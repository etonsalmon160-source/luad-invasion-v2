#!/usr/bin/env bash
# 09_postrun_chain.sh —— 网格跑完后的下游接力（无人值守）
#
#   等 run_banksy_grid.sh 消失 + 02_banksy_grid.R 全部退出
#   → 数 _status.tsv 里 DONE 的行数，**不足就停在那里并写清楚**
#   → 03_consensus.R（§3.4 裁定 + §3.5 共识 + RCTD 敏感性臂）
#   → 04_annotate.R （§3.6 域注释，中性域名）
#   → 06_coverage_guard.R（§5.2/§5.3，便宜，先跑）
#   → 05_spatial_stats.py（§3.7 Squidpy，最慢，放最后）
#
# 🔴 **不含 §5.1 深度配平守卫**：其做法（逐档重跑分域 vs 共同深度下采样）未逐条签，须先裁定。
#
# 🔴 任一步失败**就地停下**，不静默往下走。
set -uo pipefail
WT=/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963
NICHE=/home/eto/luad_v2/results/10_niche
L="$WT/logs/_postrun.log"
export R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$L"; }

cd "$WT" || exit 1
say "下游接力启动；等待上游链路（chain_hvg_then_grid.sh）走完冻结→冒烟→全网格"

# —— 1. 等上游整条链路结束 ——
# 🔴 必须把 chain_hvg_then_grid.sh 也纳入等待条件：它在 `wait $GW` 上挂着，
#    网格还没起时若只等 run_banksy_grid.sh，本脚本会立刻误判"网格已结束"（竞态）。
while pgrep -f "chain_hvg_then_grid.sh" >/dev/null \
   || pgrep -f "run_banksy_grid.sh" >/dev/null \
   || pgrep -f "02_banksy_grid.R" >/dev/null; do sleep 60; done
say "上游链路与网格进程已全部退出"

# —— 2. 数结果，不足就停 ——
if [ ! -f "$NICHE/banksy/_status.tsv" ]; then
  say "!! 没有 $NICHE/banksy/_status.tsv ⇒ 网格根本没起，接力终止"; exit 1
fi
## 🔴 按 (切片, AGF) 的**唯一对**计：重跑过的切片会有多行（FAIL 行 + DONE 行），按行数会虚高/虚低
DONE=$(awk -F'\t' 'NR>1 && $3=="DONE"{k[$1"\t"$2]=1} END{print length(k)+0}' "$NICHE/banksy/_status.tsv")
## 「仍然没跑成」= 有非 DONE 行、且**没有** DONE 行的对（重跑成功的对不算失败）
FAIL=$(awk -F'\t' 'NR>1{ if($3=="DONE") ok[$1"\t"$2]=1; else bad[$1"\t"$2]=1 }
                   END{ n=0; for(k in bad) if(!(k in ok)) n++; print n+0 }' "$NICHE/banksy/_status.tsv")
say "网格状态：DONE(唯一对)=$DONE  仍未跑成的对=$FAIL（期望 112 = 56 张 × 2 个 AGF）"
if [ "$DONE" -lt 112 ]; then
  say "!! 只有 $DONE/112 成功 ⇒ **停在这里**，先看 logs/nice_grid_*.log 与 _status.tsv 里的失败行"
  awk -F'\t' '$3!="DONE"{print "   失败: "$0}' "$NICHE/banksy/_status.tsv" | head -20 | tee -a "$L"
  exit 2
fi
say "网格全绿（112/112）⇒ 开下游"

run () {   # run <标签> <命令...>
  local tag="$1"; shift
  say "▶ $tag"
  if "$@" >>"$L" 2>&1; then say "✅ $tag 完成"; else say "!! $tag 失败（rc=$?）⇒ 就地停"; exit 3; fi
}

run "03 §3.4 裁定 + §3.5 共识"  Rscript 10_niche/03_consensus.R --agf BOTH --rule frac_ge
run "04 §3.6 域注释"            Rscript 10_niche/04_annotate.R --agf BOTH
run "06 §5.2/§5.3 覆盖与患者级守卫" Rscript 10_niche/06_coverage_guard.R
run "05 §3.7 Squidpy 空间统计"   python3 10_niche/05_spatial_stats.py --agf BOTH

say "全部下游完成。"
## 完成凭据：守护进程 07_supervisor.sh 见到它就收工退出（不靠 grep 日志，避免误判）
touch "$NICHE/_DOWNSTREAM_DONE"

## —— 最后挂 §5.1 深度配平守卫（§13.4 重跑顺序的末位）——
## 🔴 08_depth_guard.R 从 consensus/<agf>/consensus_domains.tsv 读选中档 ⇒ 必须在 03 之后。
## 🔴 换口径后**断点表必须是全新的**：$NICHE/depth_guard 已随旧产物改名归档，此处天然是空表。
say "→ 挂 §5.1 深度配平守卫（10_depth_supervisor.sh；断点表全新，符合 §13.4）"
( cd "$WT" && setsid --fork /bin/bash -c \
    "exec nohup /bin/bash $WT/10_niche/10_depth_supervisor.sh" \
    < /dev/null > "$WT/logs/_depth_supervisor.launch.log" 2>&1 & )
say "深度守卫已脱离式挂起（监护日志 logs/_depth_supervisor.log）"
