#!/usr/bin/env bash
# 19_run_cohort_all.sh —— 借锚臂全队列驱动（SPATIAL_CNV_PREREG.md §19）
#
# 逐患者串行跑 19_run_fastcnv_cohort.R。**必须串行**：单对象峰值可到 ~40 GB，
# 本机还有背景病理进程 ⇒ 并发会 OOM。
#
# 🔴 挂服务器、不挂会话：用 setsid nohup 起本脚本，PPID=1。
#    例：setsid nohup bash 08_spatial_deconv/19_run_cohort_all.sh < /dev/null > logs/....log 2>&1 &
#
# 行为（逐条显式）：
#   · 有 manifest.json 的患者**跳过**（可安全重跑本脚本续跑）
#   · 单患者失败**不立即停**（患者互相独立），但**连续两个失败 ⇒ 停**（判为系统性故障）
#   · 每患者起跑前查可用内存，低于 MIN_FREE 则等 10 分钟再看，最多等 3 次
#
set -u

WL=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
cd "$WL" || exit 1

export R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library

LOGD="$WL/logs/fastcnv_cohort_20260929"
OUT=/home/eto/luad_v2/results/08_spatial_deconv/spatial_cnv/cohort_borrow
mkdir -p "$LOGD"

MIN_FREE=80          # GB；与 19_*.R 的守卫同一口径
ARMS=all

## 受体 15 例先跑（借锚臂是本轮新增的科学），再跑自带锚 10 例（含 borrow 对照）
PATIENTS=(P7 P8 P10 P12 P16 P18 P21 P3 P5 P13 P14 P15 P17 P19 P23 \
          P1 P2 P4 P6 P9 P11 P20 P22 P24 P25)

free_gb() { awk '/^MemAvailable:/ {printf "%d", $2/1048576}' /proc/meminfo; }

echo "==== 驱动开始 $(date '+%F %T')；患者 ${#PATIENTS[@]} 例；ARMS=$ARMS ===="
fails=(); consec=0
for P in "${PATIENTS[@]}"; do
  if [ -f "$OUT/$P/manifest.json" ]; then
    echo "[$(date '+%T')] SKIP  $P（已有 manifest.json）"; continue
  fi
  ## 内存守卫
  ok=0
  for try in 1 2 3; do
    fg=$(free_gb)
    if [ "$fg" -ge "$MIN_FREE" ]; then ok=1; break; fi
    echo "[$(date '+%T')] WAIT  $P：可用内存 ${fg} GB < ${MIN_FREE} GB，等 10 分钟（第 $try 次）"
    sleep 600
  done
  if [ "$ok" -ne 1 ]; then
    echo "[$(date '+%T')] 🔴 $P 内存始终不足（${fg} GB）⇒ 停止链条"; fails+=("$P:mem"); break
  fi

  echo "[$(date '+%T')] START $P（可用 $(free_gb) GB）"
  t0=$(date +%s)
  Rscript --vanilla 08_spatial_deconv/19_run_fastcnv_cohort.R --patient "$P" --arms "$ARMS" \
      > "$LOGD/$P.log" 2>&1
  rc=$?
  t1=$(date +%s)
  if [ "$rc" -eq 0 ]; then
    v=$(grep -o '"verdict": *"[A-Z]*"' "$OUT/$P/manifest.json" 2>/dev/null | head -1)
    echo "[$(date '+%T')] END   $P rc=0 用时 $(( (t1-t0)/60 )) min  $v"
    consec=0
  else
    echo "[$(date '+%T')] 🔴 END   $P rc=$rc 用时 $(( (t1-t0)/60 )) min ⇒ 日志 $LOGD/$P.log"
    fails+=("$P:rc$rc"); consec=$((consec+1))
    if [ "$consec" -ge 2 ]; then
      echo "[$(date '+%T')] 🔴🔴 连续 $consec 例失败 ⇒ 判为系统性故障，停止链条（不烧机时）"; break
    fi
  fi
done

echo "==== 驱动结束 $(date '+%F %T') ===="
printf '失败患者：%s\n' "${fails[*]:-无}"
for P in "${PATIENTS[@]}"; do
  if [ -f "$OUT/$P/manifest.json" ]; then echo "  完成 $P"; else echo "  未完成 $P"; fi
done
