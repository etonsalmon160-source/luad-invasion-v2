#!/usr/bin/env bash
# chain_hvg_then_grid.sh —— 无人值守链路（2026-10-01 用户签 §3.4.1 补-1：全 56 张 × 两 AGF 全跑）
#
#   等 01_hvg_freeze.R 跑完
#   → 验 hvg_3000.txt 有 3000 行
#   → 用驱动器对 1 张小切片冒烟（--smoke；顺带验证驱动器本身）
#   → 冒烟两版都 DONE ⇒ 开全量：56 张 × AGF 两版，并行 4
#
# 🔴 任一步不过就**停在那里并写清楚原因**，不静默往下走。
set -uo pipefail
WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
NICHE=/home/eto/luad_v2/results/10_niche
L="$WT/logs/_chain.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$L"; }

cd "$WT" || exit 1
say "链路启动；等待 01_hvg_freeze.R 退出"

# —— 1. 等冻结 ——
while pgrep -f "01_hvg_freeze.R" >/dev/null; do sleep 60; done
say "冻结进程已退出"
if [ ! -f "$NICHE/hvg_3000.txt" ]; then
  say "!! $NICHE/hvg_3000.txt 不存在 ⇒ 冻结没成功，链路终止。见 logs/hvg_freeze_20261001.log"; exit 1
fi
N=$(wc -l < "$NICHE/hvg_3000.txt")
if [ "$N" -ne 3000 ]; then
  say "!! hvg_3000.txt 行数 $N ≠ 3000 ⇒ 终止"; exit 1
fi
say "hvg_3000.txt 校验通过（3000 行）"

# —— 2. 冒烟（1 张最小的切片，两版 AGF）——
## 🔴 冒烟写**独立**完成表 _status_smoke.tsv：它那 2 行一旦落进正式表，
##    全量的断点续跑就会把 P25 这两版当成「已完成」而跳过。开工前先清空这张表，保证真冒烟。
SMOKE_STATUS="$NICHE/banksy/_status_smoke.tsv"
say "冒烟：驱动器 × 1 张切片 × 2 个 AGF 版（完成表独立：$SMOKE_STATUS）"
printf 'GSM9226223_P25_LUAD\n' > "$WT/10_niche/_smoke_slides.txt"
rm -f "$SMOKE_STATUS"
STATUS_FILE="$SMOKE_STATUS" bash "$WT/10_niche/run_banksy_grid.sh" "$WT/10_niche/_smoke_slides.txt" BOTH 2 --smoke >>"$L" 2>&1
ok=$(awk -F'\t' '$3=="DONE"{k[$1"\t"$2]=1} END{print length(k)+0}' "$SMOKE_STATUS" 2>/dev/null)
if [ "$ok" -lt 2 ]; then
  say "!! 冒烟只有 $ok/2 版成功 ⇒ 不开全量。逐条看 logs/nice_grid_GSM9226223_P25_LUAD_agf*.log"; exit 2
fi
say "冒烟两版都成功"

# —— 3. 开全量 ——
say "开全量：56 张 × AGF 两版，并行 8（不改结果，只压墙钟）"
nohup bash "$WT/10_niche/run_banksy_grid.sh" "$WT/10_niche/slides_all56.txt" BOTH 8 >>"$L" 2>&1 &
GW=$!
say "全量驱动器 PID=$GW（父 shell）；明细见 logs/nice_grid_*.log 与 $NICHE/banksy/_status.tsv"
wait $GW
say "全量驱动器结束 rc=$?"
## 进度一律按 (切片, AGF) 的**唯一对**计（重跑会产生重复行，按行数会虚高）
awk -F'\t' 'NR>1 && $3=="DONE"{k[$1"\t"$2]=1} END{printf "  DONE(唯一对) = %d / 112\n", length(k)}' "$NICHE/banksy/_status.tsv" | tee -a "$L"
awk -F'\t' 'NR>1{c[$3]++} END{for(k in c) printf "  %s = %d 行\n", k, c[k]}' "$NICHE/banksy/_status.tsv" | tee -a "$L"
