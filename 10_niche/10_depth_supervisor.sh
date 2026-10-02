#!/usr/bin/env bash
# 10_depth_supervisor.sh —— M6 §5.1 深度配平守卫的无人值守守护（断点续跑 + 内存自适应）
#
# 挂法（必须**脱离会话子进程**，启动后核验 PPID=1）：
#   setsid --fork /bin/bash -c 'exec nohup /bin/bash <abs>/10_niche/10_depth_supervisor.sh' \
#     < /dev/null > <abs>/logs/_depth_supervisor.launch.log 2>&1 &
#
# 职责 —— **只管流程与资源，不碰任何口径**：
#   1) 盯驱动器 `run_depth_guard.sh`：它退出而任务没跑完就原地重挂；
#      已完成的 (切片, AGF, 度量) 由 `_status.tsv` **断点跳过**，不重算。
#   2) 内存自适应：可用 < 派活线 ⇒ 置暂停旗标（驱动器不再起新任务，**等**，不丢活）；
#                  可用 < 硬地板 ⇒ 收掉**最新起的**那个任务（投入最少，会被续跑重算）。
#   3) 全程写监视日志与内存时间线；本地电脑关掉不影响。
#
# 🔴 完成凭证 = `depth_guard/_status.tsv` 里 STATUS=DONE 的**唯一三元组**（重跑会多行）。
set -uo pipefail

WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
NICHE=/home/eto/luad_v2/results/10_niche
RLIB=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library
STATUS="$NICHE/depth_guard/_status.tsv"
PAUSE="$NICHE/depth_guard/_PAUSE"
DONE_FLAG="$NICHE/depth_guard/_DEPTH_DONE"
LOGD="$WT/logs"
SUPLOG="$LOGD/_depth_supervisor.log"
MON="$NICHE/depth_guard/memory_timeline.tsv"
SLIDES_F="$WT/10_niche/slides_all56.txt"

NPAR="${NPAR:-6}"        # 并发度只是资源旋钮，**不改变任何结果**（每张切片独立）
METRIC="${METRIC:-nUMI}" # 深-3：nUMI 为主
POLL=30
START_FLOOR_GB=60        # 派活线：低于它就不再起新任务（等）
HARD_FLOOR_GB=40         # 硬地板：低于它收掉最新起的任务
MAX_HOURS=96
MAX_RELAUNCH=20
MIN_GAP_S=120

say()       { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$SUPLOG"; }
free_gb()   { awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo; }
## 🔴 不能写 `pgrep -c ... || echo 0`：pgrep 计数为 0 时既打印 0 又返回 1 ⇒ 变量里两行。
r_running() { pgrep -cf -- 'bin/exec/R .*08_depth_guard\.R' 2>/dev/null || true; }
drv_alive() { pgrep -f 'run_depth_guard.sh' >/dev/null 2>&1 \
           || pgrep -f '08_depth_guard.R'   >/dev/null 2>&1; }
done_uniq() { awk -F'\t' 'NR>1 && $4=="DONE"{k[$1"\t"$2"\t"$3]=1} END{print length(k)+0}' "$STATUS" 2>/dev/null; }

proc_start() { sed 's/.*) //' "/proc/$1/stat" 2>/dev/null | awk '{print $20}'; }
reap_newest() {
  local p st best="" bt=""
  for p in $(pgrep -f '08_depth_guard\.R' 2>/dev/null); do
    st=$(proc_start "$p"); [[ -z "$st" ]] && continue
    if [[ -z "$bt" || "$st" -gt "$bt" ]]; then bt="$st"; best="$p"; fi
  done
  [[ -z "$best" ]] && return 1
  say "🔴 内存低于硬地板：收掉最新起的任务 pid=$best（会被断点续跑重算，不丢结果）"
  kill -TERM "$best" 2>/dev/null
  return 0
}

launch_driver() {
  say "重挂深度守卫驱动器（并行 $NPAR；已 DONE 的单元驱动会自动跳过）"
  ( cd "$WT" && setsid --fork /bin/bash -c \
      "exec nohup /bin/bash $WT/10_niche/run_depth_guard.sh $SLIDES_F BOTH $NPAR" \
      < /dev/null > "$WT/logs/_depth_driver.launch.log" 2>&1 & )
}

mkdir -p "$LOGD" "$NICHE/depth_guard"
[[ -f "$MON" ]] || printf 'ts\tfree_gb\tdone_uniq\tr_tasks\tpaused\n' > "$MON"

NS=$(grep -vc '^\s*$' "$SLIDES_F" 2>/dev/null || echo 56)
TARGET=$(( NS * 2 ))     # 2 个 AGF 版 × 1 个度量（nUMI）
say "深度守卫守护启动：目标 $TARGET（${NS} 张 × 2 AGF × 1 度量）；派活线 ${START_FLOOR_GB} GB；硬地板 ${HARD_FLOOR_GB} GB；并行 $NPAR；轮询 ${POLL}s"

t0=$(date +%s); tries=0; last=0
while :; do
  now=$(date +%s)
  f=$(free_gb); d=$(done_uniq); r=$(r_running)
  paused=0; [[ -e "$PAUSE" ]] && paused=1
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F' '%H:%M:%S)" "$f" "$d" "$r" "$paused" >> "$MON"

  # —— 1. 内存自适应 ——
  if (( f < START_FLOOR_GB )); then
    [[ -e "$PAUSE" ]] || { touch "$PAUSE"; say "⚠️ 可用 ${f} GB < 派活线 ${START_FLOOR_GB} GB ⇒ 置暂停旗标，不再起新任务（在跑的不动）"; }
  else
    [[ -e "$PAUSE" ]] && { rm -f "$PAUSE"; say "✅ 可用 ${f} GB 已恢复 ⇒ 撤暂停旗标，继续派活"; }
  fi
  if (( f < HARD_FLOOR_GB )); then reap_newest && sleep 15; fi

  # —— 2. 全跑完 ⇒ 先出汇总，再收工 ——
  if (( d >= TARGET )); then
    say "🎉 深度守卫全部完成 DONE=$d/$TARGET ⇒ 跑汇总（11_depth_summary.R）"
    env R_LIBS="$RLIB" Rscript "$WT/10_niche/11_depth_summary.R" >> "$LOGD/_depth_summary.log" 2>&1
    say "汇总 rc=$?（日志 $LOGD/_depth_summary.log）"
    rm -f "$PAUSE"; : > "$DONE_FLAG"; exit 0
  fi

  # —— 3. 驱动器守护 ——
  if ! drv_alive && (( now - last >= MIN_GAP_S )); then
    if (( tries < MAX_RELAUNCH )); then
      launch_driver; tries=$((tries+1)); last=$now
    elif (( tries == MAX_RELAUNCH )); then
      say "🔴 已重挂 $MAX_RELAUNCH 次仍没跑完 ⇒ 转只监视，不再重挂（看 _depth_driver.launch.log 与逐张日志）"
      tries=$((tries+1))
    fi
  fi

  # —— 4. 兜底超时 ——
  if (( now - t0 > MAX_HOURS * 3600 )); then
    say "⚠️ 守护运行超过 ${MAX_HOURS} 小时，退出（进度 DONE=$d/$TARGET）"
    rm -f "$PAUSE"; exit 9
  fi

  sleep "$POLL"
done
