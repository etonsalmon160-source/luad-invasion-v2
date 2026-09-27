#!/bin/bash
# 脱离式跑 CytoTRACE（预注册 16_cytotrace_prereg.md 的 17_cytotrace_run.py）。
#
# 为什么脱离：两档合计约 1 小时，会话断了不能死。
# 避坑：`exec` 让 bash 被 python 取代 ⇒ PID 全程不变，脚本自己写 .pid，不依赖 `$!`。
#
# 内存：Ms 稠密化后约 12 GB（133,384 × 11,240 nnz/细胞），加 X/图/瞬时，
#       峰值估 20-25 GB。与同时段在跑的 RCTD 口径 d（~17 GB）共存无压力（机器 256 GB）。
set -u

WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
OUT=/home/eto/luad_v2/results/09_trajectory/cytotrace
LOG=$OUT/run.log
PIDF=$OUT/run.pid
mkdir -p "$OUT"

cd "$WT" || exit 1
exec >>"$LOG" 2>&1

echo "=============================================================="
echo "[$(date '+%F %T')] 脱离式启动 CytoTRACE（原文那条轴）"
echo "  脚本  : 09_trajectory/17_cytotrace_run.py"
echo "  预注册: 09_trajectory/16_cytotrace_prereg.md（已签字）"
echo "  两档  : 主=全部上皮 133,384 / 敏感性=肺泡类 118,986"
echo "  输出  : results/09_trajectory/cytotrace/"
echo "  本 PID: $$（exec 后 PID 不变，就是 python 的 PID）"
echo "=============================================================="

echo "$$" > "$PIDF"

exec /home/eto/venvs/scmalig/bin/python 09_trajectory/17_cytotrace_run.py
