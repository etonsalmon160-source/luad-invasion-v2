#!/bin/bash
# 脱离式跑 RealTimeKernel 的补充配置 self_transitions='all'（唯一让矩阵不三角的设定）。
#
# 为什么要脱离：原来的跑法是本会话的子进程，会话一断就死（用户要关机）。
#
# 🔴 避坑（上次踩过、记忆里也记过）：`setsid nohup cmd &` 的 `$!` 是**外层壳**的 PID，
#    不是 python 的。这里用 `exec` 让 bash 直接被 python 取代 ⇒ **PID 全程不变**，
#    脚本自己把这个 PID 写进 .pid 文件，不依赖任何 `$!`。
#
# 内存：本配置峰值预计与 threshold=None 那档同量级（实测 95.84 GB）。脚本自带 150 GB
# 自检闸（MEM_LIMIT_GB），保护同时段在跑的脱离队列（R，40 GB）。
set -u

WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
LOG=/home/eto/luad_v2/results/09_trajectory/smoke/realtime_kernel_all_run.log
PIDF=/home/eto/luad_v2/results/09_trajectory/smoke/realtime_kernel_all_run.pid

cd "$WT" || exit 1
exec >>"$LOG" 2>&1

echo "=============================================================="
echo "[$(date '+%F %T')] 脱离式启动 RealTimeKernel 'all' 配置"
echo "  工作树 : $WT"
echo "  脚本   : 09_trajectory/10_realtime_kernel_smoke.py"
echo "  配置   : RK_CONFIGS=all_harmony_w05_noThresh（self_transitions=all, conn_weight=0.5, threshold=None）"
echo "  输出   : results/09_trajectory/smoke/realtime_kernel_structure_all.json"
echo "  本 PID : $$（下面 exec 之后 PID 不变，就是 python 的 PID）"
echo "=============================================================="

export RK_CONFIGS=all_harmony_w05_noThresh
export RK_OUT=/home/eto/luad_v2/results/09_trajectory/smoke/realtime_kernel_structure_all.json
export RK_K_FULL=12

echo "$$" > "$PIDF"

exec /home/eto/venvs/scmalig/bin/python 09_trajectory/10_realtime_kernel_smoke.py
