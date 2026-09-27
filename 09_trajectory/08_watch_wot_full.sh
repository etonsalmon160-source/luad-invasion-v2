#!/usr/bin/env bash
# 盯 WOT 全量跑的真进程：记录内存轨迹 + 超阈值时保护性终止 + 结束时播报。
#
# 为什么要有这个文件（踩过的坑，别再踩）：
#   之前两次都用 `setsid nohup cmd &` 后的 `$!` 当 PID，那是**外层壳**的 PID，
#   壳一退（几秒内）监视脚本就以为"进程已结束"，于是：
#     - 内存闸从头到尾没保护过任何东西；
#     - 我据此报了一次"进程无声死掉"，纯属误判 —— 真正的 python 一直活着。
#   所以这里的 PID 是从 ps 里按命令行认出来的**真进程**，不取 $!。
#
# 阈值 120GB 的取法：后台还有 P4 inferCNV（R）常驻约 40GB；
#   120 + 40 = 160 < 256GB 物理内存，留够余量不会波及它。

set -u

TARGET_PID="${1:?用法: 08_watch_wot_full.sh <真PID> [阈值GB]}"
THRESH_GB="${2:-120}"
OUT=/home/eto/luad_v2/results/09_trajectory
RSS_FILE="$OUT/wot_full_watch_rss.tsv"
WD_LOG="$OUT/wot_full_watchdog.log"
RUN_LOG="$OUT/wot_full_run.log"
SRC=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/09_trajectory

say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$WD_LOG"; }

if ! kill -0 "$TARGET_PID" 2>/dev/null; then
  say "PID=$TARGET_PID 不存在，什么都不做"
  exit 1
fi

say "闸挂上 真PID=$TARGET_PID 阈值 ${THRESH_GB}GB（这次是真进程，不是 \$!）"
printf 'iso\trss_gb\n' > "$RSS_FILE"

peak=0
while kill -0 "$TARGET_PID" 2>/dev/null; do
  rss_kb=$(awk '/VmRSS/{print $2}' "/proc/$TARGET_PID/status" 2>/dev/null)
  if [ -n "${rss_kb:-}" ]; then
    rss_gb=$(awk -v k="$rss_kb" 'BEGIN{printf "%.2f", k/1024/1024}')
    printf '%s\t%s\n' "$(date '+%H:%M:%S')" "$rss_gb" >> "$RSS_FILE"
    awk -v a="$rss_gb" -v b="$peak" 'BEGIN{exit !(a>b)}' && peak="$rss_gb"
    if awk -v r="$rss_gb" -v t="$THRESH_GB" 'BEGIN{exit !(r>t)}'; then
      say "!!! 内存 ${rss_gb}GB 越过阈值 ${THRESH_GB}GB —— 保护性终止（先 TERM 后 KILL）"
      kill -TERM "$TARGET_PID" 2>/dev/null
      sleep 20
      kill -0 "$TARGET_PID" 2>/dev/null && kill -9 "$TARGET_PID" 2>/dev/null
      say "已终止；峰值 ${peak}GB。这是内存闸主动拦的，不是无声死。"
      break
    fi
  fi
  sleep 10
done

say "进程已结束 峰值 ${peak}GB —— 采样 ${RSS_FILE}"
echo "=== 运行日志尾部 ===" >> "$WD_LOG"
grep -v 'No GPU/TPU' "$RUN_LOG" | tail -25 >> "$WD_LOG"
echo "=== 产出目录 ===" >> "$WD_LOG"
ls -la "$OUT/wot_full/tmaps/" >> "$WD_LOG" 2>&1
echo "=== 验收结果 ===" >> "$WD_LOG"
if [ -f "$OUT/wot_full/wot_full_manifest.json" ]; then
  /home/eto/venvs/scmalig/bin/python - "$OUT/wot_full/wot_full_manifest.json" >> "$WD_LOG" 2>&1 <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
v = m.get("result", {}).get("verification", {})
print("all_checks_pass =", v.get("all_checks_pass"))
print("n_h5ad_in_dir   =", v.get("n_h5ad_in_dir"))
for r in v.get("records", []):
    print(f"  {r['pair']:12s} exists={r['exists']} shape={r.get('shape')} "
          f"expected={r['expected_shape']} X_sum={r.get('X_sum')} nan={r.get('X_has_nan')}")
print("errors =", m.get("errors"))
print("total_wall_sec =", m.get("total_wall_sec"), "peak_rss_gb =", m.get("peak_rss_gb"))
PY
else
  echo "manifest 不存在 —— 没跑完或中途出错" >> "$WD_LOG"
fi
