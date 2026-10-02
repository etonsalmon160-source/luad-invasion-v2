#!/bin/bash
# _run_postrun_traced.sh —— 给下游接力套一层「验尸」包装
#
# 为什么要有它：2026-10-02 下游链条连续三次在启动后 7–8 分钟**无声消失**
#   （实例1 09:24:20→09:32:50；实例2 09:32:51→09:39:43），日志里既没有报错也没有
#   「失败」行，内核日志 `oom_kill 0`、无 coredump、内存全程 252 GB 富余。
#   ⇒ 说明它**不是自己退出的**，是被外部干掉的（或者 rc≠0 但没人记）。
#   所以这一层只干一件事：把**退出码 / 信号 / 存活秒数**写进 `_postrun_exit.log`。
#
# 🔴 这层不碰任何口径，只是记录。真正的链条仍是 09_postrun_chain.sh。
WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
ELOG="$WT/logs/_postrun_exit.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" >> "$ELOG"; }

START_H="$(date '+%F %H:%M:%S')"; START_S=$(date +%s)
say "===== 下游实例启动 pid=$$ 起于 $START_H ====="

on_exit() {
  local rc=$? sig=0
  (( rc > 128 )) && sig=$((rc - 128))
  if (( sig > 0 )); then
    say "🔴 链条结束 rc=$rc ⇒ 被信号 $sig（$(kill -l "$sig" 2>/dev/null)）杀掉；存活 $(( $(date +%s) - START_S )) 秒"
  else
    say "链条结束 rc=$rc（正常退出码）；存活 $(( $(date +%s) - START_S )) 秒"
  fi
}
trap on_exit EXIT
trap 'say "⚠️ 收到 SIGTERM（外部要求终止）"; exit 143' TERM
trap 'say "⚠️ 收到 SIGHUP（父会话断开？）"; exit 129' HUP
trap 'say "⚠️ 收到 SIGINT"; exit 130' INT

/bin/bash "$WT/10_niche/09_postrun_chain.sh"
rc=$?
say "内层链条返回 rc=$rc"
exit $rc
