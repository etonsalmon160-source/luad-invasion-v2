#!/usr/bin/env bash
# 03_cnv/run_floor_ab.sh <sample_id> [...]
#
# 地板 A/B 对照臂：对每个样本跑**两次**同一份 01_copykat_gse308103.R ——
#   臂 A：LUAD_COVERAGE_FLOOR=840（主口径）
#   臂 B：LUAD_COVERAGE_FLOOR=0  （无地板；其有效 UP.DR 由源码 :57 的覆写决定）
# 两次除地板外**一切相同**（同脚本、同参数、同 n.cores）。产物分别落在
#   results/03_cnv/gp2/            （臂 A，与正式 GP2 同目录）
#   results/03_cnv/ab_floor0/      （臂 B）
#
# 用途：查证"地板是否把有效 UP.DR 从 0.05 推到 0.10，进而抬高假阳性非整倍体"。
# 并发上限 4（每样本 n.cores=1；与正在跑的 P19 共存时要留余量）。
set -u
cd /home/eto/luad_v2

N_PAR=4
QUEUE=()
for s in "$@"; do
  QUEUE+=("A:$s")
  QUEUE+=("B:$s")
done

run_one() {
  local arm="${1%%:*}" s="${1##*:}"
  if [ "$arm" = "A" ]; then
    LUAD_COVERAGE_FLOOR=840 LUAD_GP2_OUTROOT=/home/eto/luad_v2/results/03_cnv/gp2 \
      Rscript 03_cnv/01_copykat_gse308103.R "$s"
  else
    LUAD_COVERAGE_FLOOR=0 \
      LUAD_GP2_TIERS=/home/eto/luad_v2/results/03_cnv/prereg_gene_tiers.csv \
      LUAD_GP2_OUTROOT=/home/eto/luad_v2/results/03_cnv/ab_floor0 \
      Rscript 03_cnv/01_copykat_gse308103.R "$s"
  fi
  echo "[ab] 臂$arm $s 退出码=$?"
}

for item in "${QUEUE[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$N_PAR" ]; do sleep 5; done   # 简易信号量
  # ⚠️ 不可写 `nohup run_one ...`：nohup 只能执行**外部程序**，不能执行 shell 函数
  # （实测报 "nohup: failed to run command 'run_one'"）。父进程已 setsid+nohup，子进程直接继承。
  run_one "$item" >"/tmp/ab_${item%%:*}_${item##*:}.log" 2>&1 &
done
wait
echo "[ab] 全部完成"
