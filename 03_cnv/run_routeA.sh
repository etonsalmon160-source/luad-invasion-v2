#!/usr/bin/env bash
# 03_cnv/run_routeA.sh —— 路线 A（覆盖度地板负对照）脱离会话运行
#
# 为什么要等冒烟驱动：两者都用 n.cores=1 且争同一块磁盘，并发会污染 route A 的
# 计时与 copykat 内部的并行距离计算。故先等 logs/gp1_runner.log 出现结束横幅。
# 为什么用 setsid：本地电脑关机 / app 关闭都不应影响服务器上的这个进程。
#
# 启动（脱离会话）：
#   cd /home/eto/luad_v2 && setsid nohup bash 03_cnv/run_routeA.sh \
#       > logs/routeA_runner.log 2>&1 < /dev/null &
# 校验：ps -o pid,ppid,pgid,sid,cmd -C R

set -u
cd /home/eto/luad_v2 || exit 1
mkdir -p logs

SAMPLE="${1:-P13_Normal}"
MARK='GP1 串行冒烟全部结束'
MAX_WAIT=$((12 * 3600))          # 最多等 12 h，超时则继续并由脚本自身硬停报错

echo "================================================================"
echo "路线 A 驱动启动  $(date '+%F %T')  pid=$$  样本=$SAMPLE"
echo "================================================================"

waited=0
while ! grep -q "$MARK" logs/gp1_runner.log 2>/dev/null; do
  if [ "$waited" -ge "$MAX_WAIT" ]; then
    echo "---- 等待冒烟驱动超时（${MAX_WAIT}s），继续尝试 ----"
    break
  fi
  if [ "$waited" -eq 0 ]; then
    echo "---- 等待冒烟驱动结束中 $(date '+%F %T') ----"
  fi
  sleep 60
  waited=$((waited + 60))
done
echo "---- 冒烟驱动已结束，开始路线 A $(date '+%F %T')（等待 ${waited}s）----"

Rscript 03_cnv/05_routeA_coverage_floor.R "$SAMPLE" > "logs/routeA_${SAMPLE}.log" 2>&1
rc=$?
echo "---- 路线 A 结束 rc=$rc $(date '+%F %T') ----"
[ -f "results/03_cnv/routeA/curve.csv" ] || { echo "!!! 曲线未产出，日志尾部："; tail -30 "logs/routeA_${SAMPLE}.log"; }

echo "================================================================"
echo "路线 A 全部结束  $(date '+%F %T')"
echo "================================================================"
