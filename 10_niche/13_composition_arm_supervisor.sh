#!/bin/bash
# 13_composition_arm_supervisor.sh —— M6 §14 D4 并发驱动（断点 + 内存守卫 + 逐张核进度）
#
# 🔴 必须用**绝对路径 /bin/bash** 调用：/home/eto/.local/bin/bash 是 #!/bin/sh 转发壳，
#    会丢掉 export -f 的函数（BASH_FUNC_*）⇒ 本脚本一律**自派发**，不靠导出函数。
#
# 用法（脱离式挂载）：
#   setsid --fork /bin/bash -c "exec nohup /bin/bash <abs>/13_composition_arm_supervisor.sh" \
#     < /dev/null > <abs>/logs/d4_supervisor.log 2>&1 &
#   然后核 PPID=1。

set -uo pipefail
SELF="$(readlink -f "$0")"
ROOT=/home/eto/luad_v2
NICHE=$ROOT/results/10_niche
DIAG=$NICHE/kstar_diag
PARTS=$DIAG/d4_parts
LOGD=$ROOT/logs/d4_parts
mkdir -p "$PARTS" "$LOGD"
RLIB=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library
SCRIPT=$ROOT/.claude/worktrees/vigilant-hawking-799963/10_niche/13_composition_arm.R
JOBS=8
MIN_FREE_GB=30

log(){ echo "[$(date +%H:%M:%S)] $*"; }
free_gb(){ awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo; }

__one(){
  local s="$1" t0 rc t1
  if [ -f "$PARTS/$s.tsv" ]; then log "跳过（已有）$s"; return 0; fi
  if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
    log "!! 可用内存 $(free_gb) GB < $MIN_FREE_GB GB ⇒ 本张不动，等下一轮"; return 9
  fi
  t0=$(date +%s)
  R_LIBS=$RLIB Rscript "$SCRIPT" --slide "$s" >"$LOGD/$s.log" 2>&1
  rc=$?
  t1=$(date +%s)
  if [ $rc -eq 0 ] && [ -f "$PARTS/$s.tsv" ]; then
    printf '%s\t0\t%d\t%s\n' "$s" "$((t1-t0))" "$(date +%F' '%T)" >>"$DIAG/_d4_status.tsv"
    log "OK   $s  ($((t1-t0))s)"
  else
    printf '%s\t%s\t%d\t%s\n' "$s" "$rc" "$((t1-t0))" "$(date +%F' '%T)" >>"$DIAG/_d4_status.tsv"
    log "FAIL $s  rc=$rc  (见 $LOGD/$s.log)"
  fi
  return $rc
}

__main(){
  log "==== D4 并发驱动开始；并发=$JOBS，可用内存 $(free_gb) GB ===="
  : >"$DIAG/_d4_status.tsv"
  ls "$NICHE/banksy" | grep -E '^GSM' >"$DIAG/_d4_slides.txt"
  local n; n=$(wc -l <"$DIAG/_d4_slides.txt")
  log "切片总数 $n"

  # 自派发并发（不用 export -f）
  xargs -a "$DIAG/_d4_slides.txt" -P "$JOBS" -I{} /bin/bash "$SELF" __one {}

  # 🔴 逐张核进度：应有 vs 已有（绝不只看计数）
  ls "$PARTS" 2>/dev/null | sed 's/\.tsv$//' | sort >"$DIAG/_d4_done.txt"
  sort "$DIAG/_d4_slides.txt" >"$DIAG/_d4_exp.txt"
  local miss; miss=$(comm -23 "$DIAG/_d4_exp.txt" "$DIAG/_d4_done.txt" | tr '\n' ' ')
  log "应有 $(wc -l <"$DIAG/_d4_exp.txt") / 已有 $(wc -l <"$DIAG/_d4_done.txt")"
  if [ -n "${miss// }" ]; then log "!! 缺：$miss"; else log "全部齐了"; fi

  log "---- 汇总 ----"
  R_LIBS=$RLIB Rscript "$SCRIPT" --combine
  log "==== D4 结束；可用内存 $(free_gb) GB ===="
}

case "${1:-__main}" in
  __one) shift; __one "$@" ;;
  __main|*) __main ;;
esac
