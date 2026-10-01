#!/usr/bin/env bash
# 20b_wait_then_report.sh —— 等全队列驱动跑完，自动出全队列上报
#
# 为什么单独挂一个守卫：驱动脚本正在被 bash 逐行读，**不能**中途编辑它去接下一步。
# 本脚本只做两件事：等驱动进程消失 → 跑只读上报 20_report_cohort.R。
#
# 🔴 挂服务器、不挂会话：setsid nohup bash 08_spatial_deconv/20b_wait_then_report.sh < /dev/null > logs/....log 2>&1 &
set -u

WL=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
cd "$WL" || exit 1

export R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library

LOGD="$WL/logs/fastcnv_cohort_20260929"
OUT=/home/eto/luad_v2/results/08_spatial_deconv/spatial_cnv/cohort_borrow
mkdir -p "$LOGD"

WAIT_MAX=43200      # 12 小时上限，超时只记一笔、不强杀任何东西
POLL=60

echo "==== 守卫开始 $(date '+%F %T')；等 19_run_cohort_all.sh 退出 ===="
t0=$(date +%s)
while pgrep -f '19_run_cohort_all\.sh' > /dev/null 2>&1; do
  now=$(date +%s)
  if [ $((now - t0)) -gt "$WAIT_MAX" ]; then
    echo "[$(date '+%T')] ⏱ 等待超过 $((WAIT_MAX/3600)) 小时，驱动仍在跑 ⇒ 本守卫退出（不干预）"
    exit 2
  fi
  sleep "$POLL"
done
echo "[$(date '+%T')] 驱动已退出（等待 $(( ($(date +%s)-t0)/60 )) 分钟）"

n=$(ls "$OUT"/P*/manifest.json 2>/dev/null | wc -l)
echo "已完成 $n / 25 例"
if [ "$n" -eq 0 ]; then
  echo "🔴 一例都没完成 ⇒ 不出上报"
  exit 1
fi

echo "[$(date '+%T')] 跑 20_report_cohort.R（只读：不筛数据、不下恶性判定）"
Rscript --vanilla 08_spatial_deconv/20_report_cohort.R > "$LOGD/report.log" 2>&1
rc=$?
echo "[$(date '+%T')] 上报结束 rc=$rc ⇒ 日志 $LOGD/report.log"
echo "==== 守卫结束 $(date '+%F %T') ===="
