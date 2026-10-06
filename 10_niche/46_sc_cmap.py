#!/usr/bin/env python3
# 46_sc_cmap.py —— 把"单细胞层面细胞类型特异签名"送进 CMap（L1），并对照已有零模型
#
# 用法: python3 46_sc_cmap.py <signature.tsv> <label>
#   signature.tsv 需含列 gene,lfc
# 零模型：直接用 results/10_niche/target_reversal/null_smatched.tsv（200 次，S 配平）
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, pandas as pd
from importlib.machinery import SourceFileLoader
HERE = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche"
TR = SourceFileLoader("tr", f"{HERE}/24_target_reversal.py").load_module()
ROOT = "/home/eto/luad_v2"; TRD = f"{ROOT}/results/10_niche/target_reversal"
SC = f"{ROOT}/results/10_niche/tr_singlecell"; DATA = f"{TRD}/data"
top = int(os.environ.get("TOP", 150))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

sig_f = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else "sc"
S = pd.read_csv(sig_f, sep="\t").dropna()
UP = S.sort_values("lfc", ascending=False).gene.head(top).tolist()
DN = S.sort_values("lfc").gene.head(top).tolist()
log(f"{label}: 签名 up {len(UP)} / down {len(DN)}")

from cmapPy.pandasGEXpress.parse_gctx import parse
G = parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df = G.data_df
genes = [str(g) for g in df.index]
X = df.values.T.astype(np.float32); sig_ids = np.array([str(c) for c in df.columns]); N = X.shape[1]
del df, G
rank, sv, _ = TR.build_rank_tables(X); del X
gpos = {g: i for i, g in enumerate(genes)}
gi = TR.load_gene_map("GSE70138")
sym2e = dict(zip(gi.pr_gene_symbol.astype(str), gi.pr_gene_id.astype(str)))
log(f"库 {len(sig_ids):,} × {N:,}")

si = TR.load_sig_info("GSE70138").set_index("sig_id").reindex(sig_ids)
sm = pd.read_csv(f"{DATA}/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",
                 sep="\t", compression="gzip").set_index("sig_id").reindex(sig_ids)
pi2 = pd.read_csv(f"{DATA}/lincs_pert_info2.tsv", sep="\t").set_index("pert_id")
lab = pd.DataFrame({"sig_id": sig_ids, "cell_id": si.cell_id.fillna("NA").values,
                    "pert_type": si.pert_type.fillna("NA").values, "pert_id": si.pert_id.values,
                    "pert_name": si.pert_iname.values})
lab["is_ts"] = pi2.is_touchstone.reindex(lab.pert_id.values).values
lab["group"] = lab.cell_id + "|" + lab.pert_type
lab["distil_ss"] = sm.distil_ss.values; lab["distil_nsample"] = sm.distil_nsample.values
qref = set(lab[(lab.is_ts == 1) & (lab.distil_nsample >= 3)]
           .sort_values(["distil_ss", "sig_id"], ascending=[False, True])
           .drop_duplicates(["cell_id", "pert_id"]).sig_id)
IDX = {s: i for i, s in enumerate(sig_ids)}
qmask = np.zeros(len(sig_ids), dtype=bool)
for s in qref: qmask[IDX[s]] = True
GRP = lab.group.values; UGRP = np.unique(GRP); GRPM = [GRP == g for g in UGRP]
log(f"Q_ref {len(qref):,}")

up_ix = np.array([gpos[sym2e[g]] for g in UP if sym2e.get(g) in gpos])
dn_ix = np.array([gpos[sym2e[g]] for g in DN if sym2e.get(g) in gpos])
log(f"映射进库：up {len(up_ix)} / down {len(dn_ix)}")
Ru, Cu = TR._rc(rank, sv, up_ix); e_u = TR.es_from_ranks(Ru, Cu, N, len(up_ix))
Rd, Cd = TR._rc(rank, sv, dn_ix); e_d = TR.es_from_ranks(Rd, Cd, N, len(dn_ix))
w = np.where(np.sign(e_u) != np.sign(e_d), (e_u - e_d) / 2.0, 0.0)
ncs = np.empty_like(w)
for m in GRPM:
    ww = w[m]
    mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
    mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
    ncs[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
tau = np.empty_like(ncs)
for m in GRPM:
    ref = np.sort(np.abs(ncs[m & qmask]))
    if ref.size == 0: continue
    tau[m] = np.sign(ncs[m]) * 100.0 / ref.size * np.searchsorted(ref, np.abs(ncs[m]), side="left")

neg = float((tau <= -90).mean()); pos = float((tau >= 90).mean())
log(f"结果：τ≤−90 = {neg:.4f}（{int((tau<=-90).sum()):,}）| τ≥+90 = {pos:.4f} | mean τ = {tau.mean():+.2f}")
NL = pd.read_csv(f"{TRD}/null_smatched.tsv", sep="\t")
v = NL.sig_le_m90.values
z = (neg - v.mean()) / v.std(ddof=1); p = (1 + (v >= neg).sum()) / (1 + len(v))
log(f"对照 S 配平零模型（n={len(v)}，均值 {v.mean():.4f}）：z = {z:+.2f}，经验 p = {p:.4f}  "
    f"{'✅ 出界' if p < 0.05 else '❌ 落随机内'}")
# 稳定性对照：阳性侧（query 自身内参照）
vp = NL.get("rate_ge_p90")
log(f"注：本查询 τ≥+90 = {pos:.4f}（对称性自检；若正负两侧都高 ⇒ 是离散度而非方向）")

# 化合物层（非 MOA）：最强的 15 个候选扰动（措辞限定）
d = pd.DataFrame({"pert_name": lab.pert_name.values, "cell_id": lab.cell_id.values,
                  "taut": tau, "ncs": ncs, "is_ts": lab.is_ts.values})
d = d.sort_values("taut")
d.to_csv(f"{TRD}/sc_{label}_ranked.tsv", sep="\t", index=False)
log("最负的 15 个签名（候选扰动）：")
print(d.head(15)[["pert_name", "cell_id", "taut", "is_ts"]].to_string(index=False))
