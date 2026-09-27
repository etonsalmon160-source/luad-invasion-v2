#!/bin/bash
# 脱离式跑运输读数（预注册 13_transport_readout_prereg.md 的 14_transport_readout.py）。
#
# 为什么脱离：跑约 1 小时，会话断了不能死。
# 避坑：`exec` 让 bash 被 python 取代 ⇒ PID 全程不变，脚本自己写 .pid，不依赖 `$!`。
#
# 内存：峰值 <20 GB（不构造 Π），与同时段在跑的 RCTD 口径 d（~17 GB）无冲突。
set -u

WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
LOG=/home/eto/luad_v2/results/09_trajectory/readout/run.log
PIDF=/home/eto/luad_v2/results/09_trajectory/readout/run.pid
mkdir -p /home/eto/luad_v2/results/09_trajectory/readout

cd "$WT" || exit 1
exec >>"$LOG" 2>&1

echo "=============================================================="
echo "[$(date '+%F %T')] 脱离式启动运输读数"
echo "  脚本 : 09_trajectory/14_transport_readout.py"
echo "  预注册: 09_trajectory/13_transport_readout_prereg.md（已签字）"
echo "  输出 : results/09_trajectory/readout/"
echo "  本 PID: $$（exec 后 PID 不变，就是 python 的 PID）"
echo "=============================================================="

echo "$$" > "$PIDF"

exec /home/eto/venvs/scmalig/bin/python 09_trajectory/14_transport_readout.py
