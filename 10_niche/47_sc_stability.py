#!/usr/bin/env python3
# 47_sc_stability.py —— 用户的核心判据：**稳定**逆转
#
# 设计：把 22 例患者随机劈成两半 → 各自**独立**推导签名 → 各自跑 CMap →
#       统计量 = 两半 τ 向量的 Spearman 相关（以及 top-500 最负签名集合的重叠）
# 零模型：**打乱 ECM_hi 标签**（同样患者数、同样细胞数，只换分组）⇒ 两半均为噪声，
#         给出"两半一致度"的零分布 ⇒ 观测若显著高于它，才是真的稳定。
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import h5py, numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from importlib.machinery import SourceFileLoader
HERE = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche"
TR = SourceFileLoader("tr", f"{HERE}/24_target_reversal.py").load_module()
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
TRD = f"{ROOT}/results/10_niche/target_reversal"; DATA = f"{TRD}/data"
TOP, K_SPLIT, SEED = 150, 6, 20261005
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ── 数据 ──
D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as f:
    genes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["genes"][:]])
    M = sparse.csr_matrix((f["counts/data"][:], f["counts/indices"][:], f["counts/indptr"][:]),
                          shape=(len(D), len(genes)))
gi = {g: i for i, g in enumerate(genes)}
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = sparse.diags(1e4 / tot) @ M; Mn = Mn.tocsr(); Mn.data = np.log1p(Mn.data)
ECM = [g for g in ["COL1A1","COL1A2","COL3A1","COL6A3","FN1","CTHRC1","POSTN","SPARC","THBS2",
                   "BGN","COL5A1","COL5A2","LUM","VCAN","MMP2","TIMP1","ASPH","COL8A1"] if g in gi]
X = np.asarray(Mn[:, [gi[g] for g in ECM]].todense())
D["ECM"] = ((X - X.mean(0)) / (X.std(0) + 1e-9)).mean(1)
fib = D[D.L2 == "Fibroblast"].copy()
Mf = Mn[fib.index.values]                      # 🔴 先取矩阵（原始行号），后 reset_index
fib = fib.reset_index(drop=True)
fib["ECM_hi"] = fib.ECM > fib.ECM.quantile(0.75)
PAT = np.array(sorted(fib.patient_id.unique())); log(f"成纤维 {len(fib):,} / 患者 {len(PAT)}")

def paired_lfc(sel_pat, lab_vec):
    """sel_pat: 患者子集；lab_vec: 布尔分组（True=高）。逐患者差 → 中位。"""
    parts = []
    for pid in sel_pat:
        pos = np.where(fib.patient_id.values == pid)[0]
        ma = pos[lab_vec[pos]]; mb = pos[~lab_vec[pos]]
        if len(ma) >= 20 and len(mb) >= 20:
            parts.append(np.asarray(Mf[ma].mean(0)).ravel() - np.asarray(Mf[mb].mean(0)).ravel())
    return np.median(np.vstack(parts), axis=0) if parts else None

# ── CMap ──
from cmapPy.pandasGEXpress.parse_gctx import parse
G = parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df = G.data_df
lgenes = [str(g) for g in df.index]; Xl = df.values.T.astype(np.float32)
sig_ids = np.array([str(c) for c in df.columns]); N = Xl.shape[1]; del df, G
rank, sv, _ = TR.build_rank_tables(Xl); del Xl
gpos = {g: i for i, g in enumerate(lgenes)}
gmap = TR.load_gene_map("GSE70138")
sym2e = dict(zip(gmap.pr_gene_symbol.astype(str), gmap.pr_gene_id.astype(str)))
si = TR.load_sig_info("GSE70138").set_index("sig_id").reindex(sig_ids)
sm = pd.read_csv(f"{DATA}/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",
                 sep="\t", compression="gzip").set_index("sig_id").reindex(sig_ids)
pi2 = pd.read_csv(f"{DATA}/lincs_pert_info2.tsv", sep="\t").set_index("pert_id")
lab = pd.DataFrame({"cell_id": si.cell_id.fillna("NA").values, "pert_type": si.pert_type.fillna("NA").values,
                    "pert_id": si.pert_id.values, "pert_name": si.pert_iname.values})
lab["is_ts"] = pi2.is_touchstone.reindex(lab.pert_id.values).values
lab["group"] = lab.cell_id + "|" + lab.pert_type
lab["ss"] = sm.distil_ss.values; lab["ns"] = sm.distil_nsample.values
qref = set(lab[(lab.is_ts == 1) & (lab.ns >= 3)].sort_values(["ss"], ascending=False)
           .drop_duplicates(["cell_id", "pert_id"]).index)
qmask = np.zeros(len(sig_ids), dtype=bool)
qmask[list(qref)] = True
GRP = lab.group.values; UGRP = np.unique(GRP); GRPM = [GRP == g for g in UGRP]
log(f"库 {len(sig_ids):,}；Q_ref {int(qmask.sum()):,}")

def tau_of(gene_up, gene_dn):
    up = np.array([gpos[sym2e[g]] for g in gene_up if sym2e.get(g) in gpos])
    dn = np.array([gpos[sym2e[g]] for g in gene_dn if sym2e.get(g) in gpos])
    Ru, Cu = TR._rc(rank, sv, up); e_u = TR.es_from_ranks(Ru, Cu, N, len(up))
    Rd, Cd = TR._rc(rank, sv, dn); e_d = TR.es_from_ranks(Rd, Cd, N, len(dn))
    w = np.where(np.sign(e_u) != np.sign(e_d), (e_u - e_d) / 2.0, 0.0)
    ncs = np.empty_like(w)
    for m in GRPM:
        ww = w[m]
        mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
        mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
        ncs[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
    t = np.empty_like(ncs)
    for m in GRPM:
        ref = np.sort(np.abs(ncs[m & qmask]))
        if ref.size == 0: continue
        t[m] = np.sign(ncs[m]) * 100.0 / ref.size * np.searchsorted(ref, np.abs(ncs[m]), side="left")
    return t, len(up), len(dn)

def sig_genes(lfc):
    r = pd.DataFrame({"gene": genes, "lfc": lfc}).replace([np.inf, -np.inf], np.nan).dropna()
    return (r.sort_values("lfc", ascending=False).gene.head(TOP).tolist(),
            r.sort_values("lfc").gene.head(TOP).tolist())

# 全队列参考（观测上限）
gu, gd = sig_genes(paired_lfc(PAT, fib.ECM_hi.values))
t_full, nu, nd = tau_of(gu, gd)
log(f"全队列：up {nu} / down {nd} ；τ≤−90 = {(t_full<=-90).mean():.4f}")
rng = np.random.default_rng(SEED)
rows = []
for k in range(K_SPLIT):
    perm = rng.permutation(PAT); h1, h2 = perm[:11], perm[11:]
    out = {"k": k}
    for tag, pats in [("h1", h1), ("h2", h2)]:
        lfc = paired_lfc(pats, fib.ECM_hi.values)
        g1, g2 = sig_genes(lfc)
        t, n1, n2 = tau_of(g1, g2)
        out[tag + "_tau"] = t; out[tag + "_n"] = (n1, n2)
        out[tag + "_frac"] = float((t <= -90).mean())
    rho = spearmanr(out["h1_tau"], out["h2_tau"]).correlation
    o1 = set(np.argsort(out["h1_tau"])[:500]); o2 = set(np.argsort(out["h2_tau"])[:500])
    out["rho_obs"] = rho; out["ovl500_obs"] = len(o1 & o2)
    # 零模型：**逐患者内**打乱分组标签（保持每患者的细胞数与高/低细胞数不变）
    lb = np.empty(len(fib), dtype=bool)
    for p in PAT:
        pos = np.where(fib.patient_id.values == p)[0]
        lb[pos] = rng.permutation(fib.ECM_hi.values[pos])
    taus = {}
    for tag, pats in [("h1", h1), ("h2", h2)]:
        g1, g2 = sig_genes(paired_lfc(pats, lb))
        t, _, _ = tau_of(g1, g2); taus[tag] = t
    out["rho_null"] = spearmanr(taus["h1"], taus["h2"]).correlation
    o1 = set(np.argsort(taus["h1"])[:500]); o2 = set(np.argsort(taus["h2"])[:500])
    out["ovl500_null"] = len(o1 & o2)
    rows.append({kk: vv for kk, vv in out.items() if not kk.endswith("_tau")})
    log(f"  split {k}: ρ_obs={out['rho_obs']:.3f} (null {out['rho_null']:.3f}) | "
        f"frac h1/h2 = {out['h1_frac']:.3f}/{out['h2_frac']:.3f} | 重叠500 = {out['ovl500_obs']} vs null {out['ovl500_null']}")

R = pd.DataFrame(rows); R.to_csv(f"{OUT}/sc_stability.tsv", sep="\t", index=False)
print("\n=== 稳定性汇总 ===")
for c in ["h1_frac", "h2_frac", "rho_obs", "rho_null", "ovl500_obs", "ovl500_null"]:
    print(f"  {c:14s} 中位 {R[c].median():.4f}   均值 {R[c].mean():.4f}")
print(f"\n观测 ρ 中位 {R.rho_obs.median():.3f} vs 零模型 {R.rho_null.median():.3f}")
print(f"观测重叠 中位 {R.ovl500_obs.median():.0f} vs 零模型 {R.ovl500_null.median():.0f}")
