#!/usr/bin/env bash
# run_depth_guard.sh —— M6 §5.1 深度配平守卫的驱动器（08_depth_guard.R 的调度壳）
#
# 一个「任务」= (切片, AGF 版, 深度度量) ⇒ 56 × 2 × 1(nUMI) = **112 个单元**。
# 每个单元内部跑 **6 个条件**（五档五分位 q1–q5 + 一个共同深度稀释 down），
# 逐个原子落盘；中途被杀也保住已完成的档，重跑自动跳过。
#
# 🔴 **本脚本不含任何口径**：只跑 08_depth_guard.R；口径全在 NICHE_PREREG.md §5.1/§5.1.1。
#
# 用法：
#   bash 10_niche/run_depth_guard.sh <slides.txt> <TRUE|FALSE|BOTH> <并行数> [--smoke]
#   例：bash 10_niche/run_depth_guard.sh 10_niche/slides_all56.txt BOTH 8
#   环境变量 STATUS_FILE 可换完成表（冒烟用别的表，免得污染全量续跑集）。
#
# 产物：logs/nice_depth_<slide>_<agfT|F>.log
#       results/10_niche/depth_guard/<agf>/<slide>__<q1..q5|down>.tsv
#       results/10_niche/depth_guard/_status.tsv
#
# 🔴 自派发（不用 `export -f`）：本机 PATH 里 `/home/eto/.local/bin/bash` 是 `#!/bin/sh`(dash)
#    转发壳，会把 `BASH_FUNC_*` 全丢掉 ⇒ 一律 `/bin/bash "$SELF" __run_one` 重新起脚本。
set -uo pipefail

DISPATCH=0
[[ "${1:-}" == "__run_one" ]] && { DISPATCH=1; shift; }

if (( DISPATCH )); then
  RUN_SLIDE="${1:?__run_one 用法: __run_one <slide> <agf> <metric>}"
  RUN_AGF="${2:?__run_one 用法: __run_one <slide> <agf> <metric>}"
  RUN_METRIC="${3:?__run_one 用法: __run_one <slide> <agf> <metric>}"
else
  SLIDES_F="${1:?用法: run_depth_guard.sh <slides.txt> <TRUE|FALSE|BOTH> <并行数> [--smoke]}"
  AGF_MODE="${2:?给 TRUE / FALSE / BOTH}"
  NPAR="${3:?给并行数（1-8）}"
  [[ "$NPAR" =~ ^[0-9]+$ ]] || { echo "并行数必须是整数，收到：$NPAR" >&2; exit 2; }
  (( NPAR > 8 )) && NPAR=8
  METRIC="${METRIC:-nUMI}"          # 深-3：nUMI 为主（nFeature 列为待签的并列臂）
fi

ROOT=/home/eto/luad_v2
WT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
NICHE="$ROOT/results/10_niche"
LOGD="$WT/logs"
RLIB=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library

## 派活线：与主网格同一条线（60 GB），比片内守卫（25 GB）更严。
## 🔴 用户明令「后台还有病理进程，不要内存溢出」⇒ 宁可等，也不并发过量。
MIN_FREE_GB=60
WAIT_MAX_S=3600

STATUS="${STATUS_FILE:-$NICHE/depth_guard/_status.tsv}"
PAUSE="$NICHE/depth_guard/_PAUSE"
PAUSE_SHARED="$NICHE/banksy/_PAUSE"    # 主网格守护的旗标：只读尊重，不替它清

mkdir -p "$LOGD" "$NICHE/depth_guard"
[[ -f "$STATUS" ]] || printf 'slide\tagf\tmetric\tstatus\tlog\n' > "$STATUS"

free_gb() { awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo; }

wait_for_mem() {
  local waited=0 fg
  while :; do
    ## 自家旗标：守护不在了就清掉，免得永久卡死。共享旗标不动（那是别人的）。
    if [[ -e "$PAUSE" ]] && ! pgrep -f '10_depth_supervisor.sh' >/dev/null 2>&1; then
      rm -f "$PAUSE"
    fi
    fg="$(free_gb)"
    if (( fg >= MIN_FREE_GB )) && [[ ! -e "$PAUSE" ]] && [[ ! -e "$PAUSE_SHARED" ]]; then
      (( waited > 0 )) && echo "[$(date +%H:%M:%S)] 内存已恢复（${fg} GB），开跑"
      return 0
    fi
    if (( waited % 60 == 0 )); then
      echo "[$(date +%H:%M:%S)] 暂缓派活：可用 ${fg} GB，旗标(pause=$( [[ -e "$PAUSE" ]] && echo 1 || echo 0 ), shared=$( [[ -e "$PAUSE_SHARED" ]] && echo 1 || echo 0 ))；已等 ${waited}s"
    fi
    sleep 10; waited=$((waited + 10))
    (( waited >= WAIT_MAX_S )) && return 1
  done
}

run_one() {
  local s="$1" agf="$2" metric="$3"
  local tag; [[ "$agf" == "TRUE" ]] && tag=agfT || tag=agfF
  local log="$LOGD/nice_depth_${s}_${tag}_${metric}.log"
  if ! wait_for_mem; then
    echo "[$(date +%H:%M:%S)] 等内存超时（${WAIT_MAX_S}s）⇒ 记 SKIP_MEM_TIMEOUT，本条需重跑" | tee -a "$log"
    printf '%s\t%s\t%s\tSKIP_MEM_TIMEOUT\t%s\n' "$s" "$agf" "$metric" "$log" >> "$STATUS"
    return 3
  fi
  echo "[$(date +%H:%M:%S)] 起 $s $tag $metric（可用 $(free_gb) GB）" | tee -a "$log"
  env R_LIBS="$RLIB" Rscript "$WT/10_niche/08_depth_guard.R" \
      --slide "$s" --agf "$agf" --metric "$metric" >> "$log" 2>&1
  local rc=$?
  if [[ $rc -eq 0 ]]; then
    printf '%s\t%s\t%s\tDONE\t%s\n' "$s" "$agf" "$metric" "$log" >> "$STATUS"
  else
    printf '%s\t%s\t%s\tFAIL(rc=%d)\t%s\n' "$s" "$agf" "$metric" "$rc" "$log" >> "$STATUS"
  fi
  echo "[$(date +%H:%M:%S)] 完 $s $tag $metric rc=$rc 可用 $(free_gb) GB" | tee -a "$log"
  return $rc
}

if (( DISPATCH )); then
  run_one "$RUN_SLIDE" "$RUN_AGF" "$RUN_METRIC"; exit $?
fi

case "$AGF_MODE" in
  TRUE)  AGFS=(TRUE) ;;
  FALSE) AGFS=(FALSE) ;;
  BOTH)  AGFS=(TRUE FALSE) ;;
  *) echo "AGF 模式只能是 TRUE / FALSE / BOTH，收到：$AGF_MODE" >&2; exit 2 ;;
esac

mapfile -t SLIDES < <(grep -v '^\s*$' "$SLIDES_F")

## 断点：完成表里 STATUS=DONE 的**唯一三元组**读成集合，派活时跳过。
## （按唯一对计，不按行计——重跑会产生重复行。）
declare -A DONESET=()
while IFS=$'\t' read -r _sl _ag _me _st _rest; do
  [[ "$_st" == "DONE" ]] && DONESET["${_sl}|${_ag}|${_me}"]=1
done < <(awk -F'\t' 'NR>1' "$STATUS")

JOBS="$(mktemp /tmp/nice_depth_jobs.XXXXXX)"
nskip=0
for s in "${SLIDES[@]}"; do for agf in "${AGFS[@]}"; do
  if [[ -n "${DONESET["${s}|${agf}|${METRIC}"]:-}" ]]; then nskip=$((nskip + 1)); continue; fi
  echo "$s $agf $METRIC" >> "$JOBS"
done; done
NTOT=$(wc -l < "$JOBS")
echo "[$(date +%H:%M:%S)] 切片 ${#SLIDES[@]} 张 × AGF ${#AGFS[@]} 版 × 度量 $METRIC；断点跳过 $nskip，本次派活 $NTOT，并行 $NPAR"

if (( NTOT == 0 )); then
  echo "[$(date +%H:%M:%S)] 没有待跑任务（全部已完成）"; rm -f "$JOBS"; exit 0
fi

export STATUS_FILE="$STATUS"
xargs -a "$JOBS" -L1 -P "$NPAR" /bin/bash "$SELF" __run_one
rc=$?
rm -f "$JOBS"
ndone=$(awk -F'\t' 'NR>1 && $4=="DONE"{c++} END{print c+0}' "$STATUS")
nuniq=$(awk -F'\t' 'NR>1 && $4=="DONE"{print $1"|"$2"|"$3}' "$STATUS" | sort -u | wc -l)
echo "[$(date +%H:%M:%S)] 全部结束 rc=$rc；完成表 $STATUS（DONE 行 $ndone / 唯一三元组 $nuniq）"
exit $rc
