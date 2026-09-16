#!/usr/bin/env bash
# 03_cnv/run_gp1_smoke.sh —— GP1 三个冒烟点，**串行**跑
#
# 为什么串行：并发会互相争内存带宽与磁盘，污染 t = a + b*cells 的计时拟合。
#               单机 20 核下每样本 n.cores=1，串行的总挂钟代价可接受。
# 为什么用 setsid：本脚本须脱离发起它的会话存活（用户可能关闭 app）。
#               启动方式见下方注释；务必确认 PPID 变为 1。
#
# 启动（脱离会话，且日志落盘）：
#   cd /home/eto/luad_v2 && setsid nohup bash 03_cnv/run_gp1_smoke.sh \
#       > logs/gp1_runner.log 2>&1 < /dev/null &
#   disown
# 校验：
#   ps -o pid,ppid,sid,cmd -C R

set -u
cd /home/eto/luad_v2 || exit 1
mkdir -p logs

POINTS=(P13_Normal P14_AIS P24_LUAD)

echo "================================================================"
echo "GP1 串行冒烟启动  $(date '+%F %T')  pid=$$"
echo "================================================================"

# ---- 步骤 0: 75 样本两档基因数（GP2 排期的输入，纯读 h5ad，约 15 min）--------
echo
echo "---- [prereg_gene_tiers] 开始 $(date '+%F %T') ----"
rm -f results/03_cnv/prereg_gene_tiers.csv
Rscript 03_cnv/02_prereg_gene_tiers.R > logs/prereg_gene_tiers.log 2>&1
echo "---- [prereg_gene_tiers] 结束 rc=$? $(date '+%F %T') ----"
[ -f results/03_cnv/prereg_gene_tiers.csv ] || echo "!!! 分档 CSV 未产出，见 logs/prereg_gene_tiers.log"

for s in "${POINTS[@]}"; do
  echo
  echo "---- [$s] 开始 $(date '+%F %T') ----"
  # 清掉上一次的残留输出，避免与本次混写（曾因重复进程写同一目录出过错）
  rm -rf "results/03_cnv/smoke/$s" "results/03_cnv/smoke/$s.json"
  t0=$SECONDS
  Rscript 03_cnv/03_smoke_test_copykat.R "$s" > "logs/gp1_${s}.log" 2>&1
  rc=$?
  echo "    挂钟 $((SECONDS - t0)) s"   # 与 JSON 内的 wall_sec_* 互为佐证
  echo "---- [$s] 结束 rc=$rc $(date '+%F %T') ----"
  if [ ! -f "results/03_cnv/smoke/$s.json" ]; then
    echo "!!! [$s] 未产出 JSON，日志尾部："
    tail -20 "logs/gp1_${s}.log"
  fi
done

echo
echo "================================================================"
echo "GP1 串行冒烟全部结束  $(date '+%F %T')"
echo "================================================================"
