#!/usr/bin/env bash
# 07_supervisor.sh —— M6 网格的无人值守守护进程（断点续跑 + 内存自适应）
#
# 挂法（必须**脱离会话子进程**，启动后核验 PPID=1）：
#   setsid --fork /bin/bash -c 'exec nohup /bin/bash <abs>/10_niche/07_supervisor.sh' \
#     < /dev/null > <abs>/logs/_supervisor.launch.log 2>&1 &
#
# 职责 —— **只管流程与资源，不碰任何口径**：
#   1) 盯链条：链条/网格意外退出就原地重挂；已完成的切片由 _status.tsv **断点跳过**，不重算。
#   2) 内存自适应：
#        · 可用内存 < 派活线 ⇒ 置暂停旗标，驱动器**不再起新任务**（等，不判跳过、不丢活）；
#        · 可用内存 < 硬地板 ⇒ 收掉**最新起的**那个任务（投入最少），它会被续跑重算。
#      ⚠️ 同一时刻**只有一个**收任务的进程组会动，避免和驱动器抢。
#   3) 全程写监视日志（时间 / 可用内存 / 进度 / 暂停态）；本地电脑关掉不影响。
#
# 🔴 完成凭证 = `_status.tsv` 里 (切片, AGF) 的 **DONE** 行，且**按唯一对计数**
#    （同一张切片重跑会出现多行，只有唯一对才是真进度）。
#    半截文件不会有 DONE 行 ⇒ 续跑必然重算它（02_banksy_grid.R 已改原子落盘）。
set -uo pipefail

WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
NICHE=/home/eto/luad_v2/results/10_niche
STATUS="$NICHE/banksy/_status.tsv"
PAUSE="$NICHE/banksy/_PAUSE"
DONE_FLAG="$NICHE/_DOWNSTREAM_DONE"
LOGD="$WT/logs"
SUPLOG="$LOGD/_supervisor.log"
MON="$NICHE/banksy/memory_timeline.tsv"

TARGET=112              # 56 张 × 2 个 AGF 版
POLL=30                 # 轮询秒数
START_FLOOR_GB=60       # 派活线：低于它就不再起新任务（等）
HARD_FLOOR_GB=40        # 硬地板：低于它收掉最新起的任务
MAX_HOURS=72            # 兜底：守护最多跑这么久
MAX_RELAUNCH=10         # 同一目标最多重挂几次，防止无限循环
MIN_GAP_S=300           # 两次重挂之间至少隔这么久（否则链条一挂就每 30 秒撞一次、烧掉重试名额）

say()       { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$SUPLOG"; }
free_gb()   { awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo; }
done_uniq() { awk -F'\t' 'NR>1 && $3=="DONE"{k[$1"\t"$2]=1} END{print length(k)+0}' "$STATUS" 2>/dev/null; }
## 🔴 不能写成 `pgrep -c ... || echo 0`：pgrep 在计数为 0 时**既打印 0 又返回 1**，
##    `|| echo 0` 于是再补一个 0 ⇒ 变量里是两行，printf 会多吐一行 `0\t0`（时间线里出现过）。
## 模式要求 `bin/exec/R` 与脚本名同时出现，免得把「命令行里恰好含这串字的 shell」也算进来。
r_running() { pgrep -cf -- 'bin/exec/R .*02_banksy_grid\.R' 2>/dev/null || true; }

grid_alive() { pgrep -f 'chain_hvg_then_grid.sh' >/dev/null 2>&1 \
            || pgrep -f 'run_banksy_grid.sh'   >/dev/null 2>&1 \
            || pgrep -f '02_banksy_grid.R'     >/dev/null 2>&1; }
post_alive() { pgrep -f '09_postrun_chain.sh' >/dev/null 2>&1; }

launch_chain() {
  say "重挂主链条（冻结→冒烟→全网格；已 DONE 的切片驱动会自动跳过）"
  ( cd "$WT" && setsid --fork /bin/bash -c \
      "exec nohup /bin/bash $WT/10_niche/chain_hvg_then_grid.sh" \
      < /dev/null > "$WT/logs/_chain.launch.log" 2>&1 & )
}
launch_post() {
  say "重挂下游接力（03 裁定/共识 → 04 注释 → 06 覆盖守卫 → 05 空间统计）"
  ## 🔴 2026-10-02：下游链条连续三次启动后 7–8 分钟**无声消失**（无报错、无 OOM、内存富余），
  ##    因此改走 `_run_postrun_traced.sh` 包装：它会把退出码 / 信号 / 存活秒数写进 _postrun_exit.log。
  ##    这层只记录、不改口径。
  ( cd "$WT" && setsid --fork /bin/bash -c \
      "exec nohup /bin/bash $WT/10_niche/_run_postrun_traced.sh" \
      < /dev/null > "$WT/logs/_postrun.launch.log" 2>&1 & )
}

## 取进程启动时刻（/proc/<pid>/stat 第 22 列，单位 clock tick）。
## 🔴 不能用 `ps -o etimes`：本机 ps 的这一列**恒为 0**（实测），排序会退化。
##    sed 贪心剥到**最后一个** ") " 后面，避免 comm 里带括号时取错列；starttime 是之后的第 20 列。
proc_start() { sed 's/.*) //' "/proc/$1/stat" 2>/dev/null | awk '{print $20}'; }

## 收掉**最新起**的那个网格任务（starttime 最大 = 刚起不久、投入最少）
reap_newest() {
  local p st best="" bt=""
  for p in $(pgrep -f '02_banksy_grid\.R' 2>/dev/null); do
    [[ "$(tr '\0' ' ' < "/proc/$p/cmdline" 2>/dev/null)" == *--smoke* ]] && continue
    st=$(proc_start "$p"); [[ -z "$st" ]] && continue
    if [[ -z "$bt" || "$st" -gt "$bt" ]]; then bt="$st"; best="$p"; fi
  done
  [[ -z "$best" ]] && return 1
  say "🔴 内存低于硬地板：收掉最新起的任务 pid=$best（会被断点续跑重算，不丢结果）"
  kill -TERM "$best" 2>/dev/null
  return 0
}

mkdir -p "$LOGD" "$NICHE/banksy"
[[ -f "$MON" ]] || printf 'ts\tfree_gb\tdone_uniq\tr_tasks\tpaused\n' > "$MON"
say "守护启动：目标 $TARGET；派活线 ${START_FLOOR_GB} GB；硬地板 ${HARD_FLOOR_GB} GB；轮询 ${POLL}s"

t0=$(date +%s); chain_tries=0; post_tries=0; last_chain=0; last_post=0
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
  if (( f < HARD_FLOOR_GB )); then
    reap_newest && sleep 15
  fi

  # —— 2. 下游全部完成 ⇒ 收工 ——
  if [[ -e "$DONE_FLAG" ]]; then
    say "🎉 下游全部完成，守护退出（DONE=$d/$TARGET）"
    rm -f "$PAUSE"; exit 0
  fi

  # —— 3. 链条守护 ——
  if (( d >= TARGET )); then
    if ! post_alive && (( now - last_post >= MIN_GAP_S )); then
      if (( post_tries < MAX_RELAUNCH )); then
        launch_post; post_tries=$((post_tries+1)); last_post=$now
      elif (( post_tries == MAX_RELAUNCH )); then
        say "🔴 下游已重挂 $MAX_RELAUNCH 次仍没跑完 ⇒ 转只监视，不再重挂（看 _postrun.log）"
        post_tries=$((post_tries+1))
      fi
    fi
  else
    if ! grid_alive && (( now - last_chain >= MIN_GAP_S )); then
      if (( chain_tries < MAX_RELAUNCH )); then
        launch_chain; chain_tries=$((chain_tries+1)); last_chain=$now
      elif (( chain_tries == MAX_RELAUNCH )); then
        say "🔴 主链条已重挂 $MAX_RELAUNCH 次仍没跑完 ⇒ 转只监视，不再重挂（看 _chain.log）"
        chain_tries=$((chain_tries+1))
      fi
    fi
  fi

  # —— 4. 兜底超时 ——
  if (( $(date +%s) - t0 > MAX_HOURS * 3600 )); then
    say "⚠️ 守护运行超过 ${MAX_HOURS} 小时，退出（进度 DONE=$d/$TARGET）"
    rm -f "$PAUSE"; exit 9
  fi

  sleep "$POLL"
done
