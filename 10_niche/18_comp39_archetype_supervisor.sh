#!/bin/bash
# 15_domain_archetype_supervisor.sh —— M6 §17 并发驱动（断点 + 内存守卫 + 逐张核进度）
# 🔴 必须用绝对路径 /bin/bash 调用（/home/eto/.local/bin/bash 是 dash 转发壳，丢 BASH_FUNC_*）⇒ 自派发

set -uo pipefail
SELF="$(readlink -f "$0")"
ROOT=/home/eto/luad_v2
DIAG=$ROOT/results/10_niche/kstar_diag
PARTS=$DIAG/d11_parts
LOGD=$ROOT/logs/d11_parts
mkdir -p "$PARTS" "$LOGD"
RLIB=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library
SCRIPT=$ROOT/.claude/worktrees/vigilant-hawking-799963/10_niche/18_comp39_archetype.R
JOBS=8
MIN_FREE_GB=25

log(){ echo "[$(date +%H:%M:%S)] $*"; }
free_gb(){ awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo; }

__one(){
  local s="$1" t0 rc t1
  [ -f "$PARTS/$s.tsv" ] && { log "跳过（已有）$s"; return 0; }
  if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then log "!! 内存 $(free_gb) GB < $MIN_FREE_GB ⇒ 本张不动"; return 9; fi
  t0=$(date +%s)
  R_LIBS=$RLIB Rscript "$SCRIPT" --slide "$s" >"$LOGD/$s.log" 2>&1; rc=$?
  t1=$(date +%s)
  if [ $rc -eq 0 ] && [ -f "$PARTS/$s.tsv" ]; then
    printf '%s\t0\t%d\t%s\n' "$s" "$((t1-t0))" "$(date +%F' '%T)" >>"$DIAG/_d11_status.tsv"; log "OK   $s ($((t1-t0))s)"
  else
    printf '%s\t%s\t%d\t%s\n' "$s" "$rc" "$((t1-t0))" "$(date +%F' '%T)" >>"$DIAG/_d11_status.tsv"; log "FAIL $s rc=$rc"
  fi
  return $rc
}

__main(){
  log "==== §20 comp39 归并 并发驱动；并发=$JOBS，可用内存 $(free_gb) GB ===="
  : >"$DIAG/_d11_status.tsv"
  ls "$ROOT/results/10_niche/banksy" | grep -E '^GSM' >"$DIAG/_d11_slides.txt"
  log "切片总数 $(wc -l <"$DIAG/_d11_slides.txt")"
  xargs -a "$DIAG/_d11_slides.txt" -P "$JOBS" -I{} /bin/bash "$SELF" __one {}
  ls "$PARTS" 2>/dev/null | sed 's/\.tsv$//' | sort >"$DIAG/_d11_done.txt"
  sort "$DIAG/_d11_slides.txt" >"$DIAG/_d11_exp.txt"
  local miss; miss=$(comm -23 "$DIAG/_d11_exp.txt" "$DIAG/_d11_done.txt" | tr '\n' ' ')
  log "应有 $(wc -l <"$DIAG/_d11_exp.txt") / 已有 $(wc -l <"$DIAG/_d11_done.txt")"
  [ -n "${miss// }" ] && log "!! 缺：$miss" || log "全部齐了"
  log "---- 汇总 ----"
  R_LIBS=$RLIB Rscript "$SCRIPT" --combine
  log "==== 结束；可用内存 $(free_gb) GB ===="
}

case "${1:-__main}" in
  __one) shift; __one "$@" ;;
  *) __main ;;
esac
