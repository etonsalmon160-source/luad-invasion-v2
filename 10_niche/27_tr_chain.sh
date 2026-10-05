#!/bin/bash
# 27_tr_chain.sh —— 靶点扰动逆向臂：**串行**跑 D5 / D7（等 D3 先结束；内存守卫）
# 🔴 串行而非并发：单进程峰值 ~17 GB（X + rank + sorted_vals），并发会叠乘
set -uo pipefail
SELF="$(readlink -f "$0")"
ROOT=/home/eto/luad_v2
PY=$ROOT/.claude/worktrees/vigilant-hawking-799963/10_niche/24_target_reversal.py
OUTBASE=$ROOT/results/10_niche/target_reversal
MIN_FREE_GB=40
log(){ echo "[$(date +%H:%M:%S)] $*"; }
free_gb(){ awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo; }

__wait_mem(){
  while [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; do
    log "内存 $(free_gb) GB < $MIN_FREE_GB GB ⇒ 等 60s"; sleep 60
  done
}

__one(){
  local q="$1"
  __wait_mem
  log "==== 开始 $q（可用内存 $(free_gb) GB）===="
  mkdir -p "$OUTBASE/L1a_$q"
  python3 "$PY" --lib L1a --query "$q" --outdir "$OUTBASE/L1a_$q"
  local rc=$?
  log "==== $q 结束 rc=$rc（可用内存 $(free_gb) GB）===="
  return $rc
}

__main(){
  # ① 等 D3 进程退出（若在跑）
  local d3pid; d3pid=$(pgrep -f 'query D3' | head -1 || true)
  if [ -n "$d3pid" ]; then
    log "等 D3（pid $d3pid）结束…"
    while kill -0 "$d3pid" 2>/dev/null; do sleep 60; done
    log "D3 已结束"
  fi
  __one D5
  __one D7
  log "==== 全链结束 ===="
}

case "${1:-__main}" in
  __one) shift; __one "$@" ;;
  *) __main ;;
esac
