#!/bin/bash
# 脱离式跑"补原文那条轴（一）"（预注册 20_paper_axis_kac_prereg.md 的 21_paper_axis_kac_run.py）。
#
# 为什么脱离：读 3.1 GB h5ad + 9 次真打分 + 200 次置换零分布，估 15-20 分钟，会话断了不能死。
# 避坑：`exec` 让 bash 被 python 取代 ⇒ PID 全程不变，脚本自己写 .pid，不依赖 `$!`。
#
# 内存：133,384 × 18,069 的 CSR 子集（~2-4 GB）+ log1p + score_genes 临时，估峰值 < 20 GB。
set -u

WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
OUT=/home/eto/luad_v2/results/09_trajectory/paper_axis
LOG=$OUT/run.log
PIDF=$OUT/run.pid
mkdir -p "$OUT"

cd "$WT" || exit 1
exec >>"$LOG" 2>&1

echo "=============================================================="
echo "[$(date '+%F %T')] 脱离式启动：补原文那条轴（一）"
echo "  脚本  : 09_trajectory/21_paper_axis_kac_run.py"
echo "  预注册: 09_trajectory/20_paper_axis_kac_prereg.md（已签字）"
echo "  签字  : T0=解冻面板(限范围) / T1=D1主判+D2佐证 / T2=全8MP都打"
echo "  判据  : D1 = 8型MP的簇级argmax ∈ {Tumor cell/KAC, KAC/inflammatory}"
echo "  对照  : K1深度 / K2置换200次 / K3单一患者"
echo "  输出  : $OUT/"
echo "  本 PID: $$（exec 后 PID 不变，就是 python 的 PID）"
echo "=============================================================="

echo "$$" > "$PIDF"

exec /home/eto/venvs/scmalig/bin/python 09_trajectory/21_paper_axis_kac_run.py
