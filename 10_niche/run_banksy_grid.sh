#!/usr/bin/env bash
# run_banksy_grid.sh —— M6 §3.2/§3.4 逐切片 BANKSY 网格的驱动器
#
# 作用：把「哪些切片 × 哪个 AGF 版」跑一遍 02_banksy_grid.R，控制并行度、内存闸门与断点续跑。
# 🔴 **本脚本不含任何口径**：切片范围与并行度由命令行给；跑出来的东西是不是正式结论，
#    取决于调用它的那份「范围」是不是预注册 §3.4 允许的。**范围本身就是待签的决断**。
#
# 用法：
#   bash 10_niche/run_banksy_grid.sh <slides.txt> <TRUE|FALSE|BOTH> <并行数> [--smoke]
#   例：bash 10_niche/run_banksy_grid.sh 10_niche/slides_all56.txt BOTH 8
#   环境变量 STATUS_FILE 可换完成表（冒烟用 _status_smoke.tsv，免得污染全量的续跑集）。
#
# 产物：logs/nice_grid_<slide>_<agfT|F>.log
#       results/10_niche/banksy/<slide>/{domains,stats,profiles}_<agf>.tsv(.gz)
#       results/10_niche/banksy/_status.tsv（逐张逐版一行：DONE / FAIL / SKIP_MEM_TIMEOUT）
#
# 🔴 2026-10-01 修正：**不再用 `export -f` + `xargs bash -c`**（那套在本机根本传不下去，已实证，
#    冒烟因此 0/2 全挂、两条链路原地退出）。改为「子命令 `__run_one` 重新起本脚本」。
#    详细根因见下方「自派发入口」注释。
#
# 🔴 2026-10-01 补（用户令：断点 + 内存自适应 + 中间结果不丢）：
#    · **断点续跑**：开工前读完成表，把 (切片, AGF) 已是 DONE 的**唯一对**从任务清单里去掉；
#      重跑同一张会出现重复行，所以进度一律按**唯一对**计。
#    · **内存自适应**：可用内存 < 派活线时**等**（不判跳过、不丢活）；等超过上限才记 SKIP_MEM_TIMEOUT。
#      守护进程 07_supervisor.sh 会在内存吃紧时置暂停旗标 `banksy/_PAUSE`，这里一并尊重。
set -uo pipefail

## ——— 自派发入口 ———
## xargs 的每个任务用 `/bin/bash <本脚本> __run_one <slide> <agf>` 重新起一份本脚本，
## 靠**命令行参数**告诉子进程干什么，完全不依赖 bash 函数导出。
##
## 🔴 为什么不能 `export -f run_one`：
##    本机 PATH 里 `/home/eto/.local/bin/bash` 是个 132 字节的 `#!/bin/sh`（dash）转发壳：
##        #!/bin/sh
##        export PATH="/home/eto/.local/bin:...:$PATH"
##        exec /bin/bash "$@"
##    dash 启动时会丢掉**名字不是合法 shell 标识符**的环境变量 ⇒ `BASH_FUNC_run_one%%`
##    到不了真正的 /bin/bash。实测（2026-10-01，逐条可复现）：
##      · `env 'BASH_FUNC_f%%=() { echo OK; }' /bin/bash -c f`        → OK
##      · `env 'BASH_FUNC_f%%=() { echo OK; }' bash -c f`             → f: command not found
##        （按名字调 bash ⇒ 先经 dash 壳 ⇒ 变量已被丢）
##      · `env 'BASH_FUNC_f%%=() { echo OK; }' /bin/sh -c 'env|grep -c BASH_FUNC'` → 0
##      · 对照：普通变量 `MY_PLAIN=x` 经 dash 转发**不丢**（计数 1）
##    ⇒ 本机上**任何**按名字调 `bash` 的场合都传不了函数；只有普通环境变量能过。
##    ⇒ 凡是要并发的脚本，一律走「重新起脚本 + 子命令」，别再试 `export -f`。
DISPATCH=0
[[ "${1:-}" == "__run_one" ]] && { DISPATCH=1; shift; }

if (( DISPATCH )); then
  RUN_SLIDE="${1:?__run_one 用法: __run_one <slide> <agf>}"
  RUN_AGF="${2:?__run_one 用法: __run_one <slide> <agf>}"
  SMOKE_FLAG="${SMOKE_FLAG:-}"    # 普通环境变量，dash 不丢（只有 BASH_FUNC_* 会被丢）
else
  SLIDES_F="${1:?用法: run_banksy_grid.sh <slides.txt> <TRUE|FALSE|BOTH> <并行数> [--smoke]}"
  AGF_MODE="${2:?给 TRUE / FALSE / BOTH}"
  NPAR="${3:?给并行数（1-8）}"
  ## 🔴 硬上限 8：调用方（chain_hvg_then_grid.sh）传数进来，这里兜底。
  ##    并发度只是资源旋钮、**不改变任何结果**（每张切片独立）。
  ##    2026-10-01 首张实测后由用户裁定 3 → 8：单任务实测 **约 2 GB / 1 核**
  ##    （日志里「已用 9.7 GB」是**整机**已用量，不是单任务），机器 256 GB / 20 核。
  ##    8 路 ≈ 8 核 / 约 20 GB，对后台病理进程仍留足余量。
  [[ "$NPAR" =~ ^[0-9]+$ ]] || { echo "并行数必须是整数，收到：$NPAR" >&2; exit 2; }
  (( NPAR > 8 )) && NPAR=8
  SMOKE_FLAG=""
  [[ "${4:-}" == "--smoke" ]] && SMOKE_FLAG="--smoke"
fi

ROOT=/home/eto/luad_v2
WT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"     # 工作树根
SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
NICHE="$ROOT/results/10_niche"
LOGD="$WT/logs"
RLIB=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library

## 派活线：可用内存低于它就不起新任务（等）。01_ 的片内守卫是 25 GB，这里更严。
## 🔴 用户明令「后台还有病理进程，不要内存溢出」⇒ 宁可等，也不并发过量。
MIN_FREE_GB=60
WAIT_MAX_S=3600         # 等内存最多等 1 小时，超时才记 SKIP_MEM_TIMEOUT（不是丢活，是记账）

## 完成表：冒烟用独立表，避免它的 2 行把全量的 P25 误判成「已完成」而跳过。
STATUS="${STATUS_FILE:-$NICHE/banksy/_status.tsv}"
PAUSE="$NICHE/banksy/_PAUSE"

mkdir -p "$LOGD" "$NICHE/banksy"
[[ -f "$STATUS" ]] || printf 'slide\tagf\tstatus\tlog\n' > "$STATUS"

free_gb() { awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo; }

## 等内存。返回 1 = 等超时。
## 暂停旗标由守护进程置/撤；若守护不在了（旗标没人撤），视为失效并清掉，免得永久卡死。
wait_for_mem() {
  local waited=0 fg
  while :; do
    if [[ -e "$PAUSE" ]] && ! pgrep -f '07_supervisor.sh' >/dev/null 2>&1; then
      rm -f "$PAUSE"
    fi
    fg="$(free_gb)"
    if (( fg >= MIN_FREE_GB )) && [[ ! -e "$PAUSE" ]]; then
      (( waited > 0 )) && echo "[$(date +%H:%M:%S)] 内存已恢复（${fg} GB），开跑"
      return 0
    fi
    if (( waited % 60 == 0 )); then
      if [[ -e "$PAUSE" ]]; then
        echo "[$(date +%H:%M:%S)] 暂缓派活：守护置了暂停旗标 $PAUSE（内存吃紧，可用 ${fg} GB）；已等 ${waited}s"
      else
        echo "[$(date +%H:%M:%S)] 暂缓派活：可用 ${fg} GB < 派活线 ${MIN_FREE_GB} GB；已等 ${waited}s"
      fi
    fi
    sleep 10; waited=$((waited + 10))
    (( waited >= WAIT_MAX_S )) && return 1
  done
}

run_one() {
  local s="$1" agf="$2"
  local tag; [[ "$agf" == "TRUE" ]] && tag=agfT || tag=agfF
  local log="$LOGD/nice_grid_${s}_${tag}.log"
  ## 🔴 内存闸门是「等」不是「跳」：不再有 SKIP_MEM 把活悄悄丢掉。
  ##    旧写法 `if free_gb -lt "$MIN_FREE_GB"` 还是把 `-lt`/`60` 当**参数**喂给函数、
  ##    函数恒返回 0 ⇒ 条件恒真 ⇒ 每个任务都被误判成内存不足（2026-10-01 已实证并撤）。
  if ! wait_for_mem; then
    echo "[$(date +%H:%M:%S)] 等内存超时（${WAIT_MAX_S}s）⇒ 记 SKIP_MEM_TIMEOUT，本条需重跑" | tee -a "$log"
    printf '%s\t%s\tSKIP_MEM_TIMEOUT\t%s\n' "$s" "$agf" "$log" >> "$STATUS"
    return 3
  fi
  echo "[$(date +%H:%M:%S)] 起 $s $tag（可用 $(free_gb) GB）" | tee -a "$log"
  env R_LIBS="$RLIB" Rscript "$WT/10_niche/02_banksy_grid.R" \
      --slide "$s" --agf "$agf" $SMOKE_FLAG >> "$log" 2>&1
  local rc=$?
  if [[ $rc -eq 0 ]]; then
    printf '%s\t%s\tDONE\t%s\n' "$s" "$agf" "$log" >> "$STATUS"
  else
    printf '%s\t%s\tFAIL(rc=%d)\t%s\n' "$s" "$agf" "$rc" "$log" >> "$STATUS"
  fi
  echo "[$(date +%H:%M:%S)] 完 $s $tag rc=$rc 可用 $(free_gb) GB" | tee -a "$log"
  return $rc
}

## 派发模式：跑到这里说明本进程是被 xargs 起出来干**一件**活的，干完就走，不碰状态表以外的任何东西。
if (( DISPATCH )); then
  run_one "$RUN_SLIDE" "$RUN_AGF"; exit $?
fi

## ——— 以下是驱动器模式 ———
case "$AGF_MODE" in
  TRUE)  AGFS=(TRUE) ;;
  FALSE) AGFS=(FALSE) ;;
  BOTH)  AGFS=(TRUE FALSE) ;;
  *) echo "AGF 模式只能是 TRUE / FALSE / BOTH，收到：$AGF_MODE" >&2; exit 2 ;;
esac

mapfile -t SLIDES < <(grep -v '^\s*$' "$SLIDES_F")

## 断点：把完成表里已是 DONE 的**唯一对**读成集合，派活时跳过。
declare -A DONESET=()
while IFS=$'\t' read -r _sl _ag _st _rest; do
  [[ "$_st" == "DONE" ]] && DONESET["${_sl}|${_ag}"]=1
done < <(awk -F'\t' 'NR>1' "$STATUS")

: > /tmp/nice_grid_jobs.$$.txt
nskip=0
for s in "${SLIDES[@]}"; do for agf in "${AGFS[@]}"; do
  if [[ -n "${DONESET["${s}|${agf}"]:-}" ]]; then nskip=$((nskip + 1)); continue; fi
  echo "$s $agf" >> /tmp/nice_grid_jobs.$$.txt
done; done
NTOT=$(wc -l < /tmp/nice_grid_jobs.$$.txt)
echo "[$(date +%H:%M:%S)] 切片 ${#SLIDES[@]} 张 × AGF ${#AGFS[@]} 版；断点跳过 $nskip，本次派活 $NTOT，并行 $NPAR（smoke=${SMOKE_FLAG:-no}）"

if (( NTOT == 0 )); then
  echo "[$(date +%H:%M:%S)] 没有待跑任务（全部已完成）"
  rm -f /tmp/nice_grid_jobs.$$.txt
  exit 0
fi

## 🔴 子命令用**绝对路径 /bin/bash**（绕开 PATH 里的 dash 转发壳）跑 `__run_one`；
##    --smoke 与完成表路径走环境变量（普通变量能过 dash，函数不能）。
export SMOKE_FLAG
export STATUS_FILE="$STATUS"
xargs -a /tmp/nice_grid_jobs.$$.txt -L1 -P "$NPAR" /bin/bash "$SELF" __run_one
rc=$?
rm -f /tmp/nice_grid_jobs.$$.txt
echo "[$(date +%H:%M:%S)] 全部结束 rc=$rc；状态表 $STATUS"
exit $rc
