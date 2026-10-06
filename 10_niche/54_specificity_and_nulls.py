#!/usr/bin/env python3
# 54_specificity_and_nulls.py —— 三件事一起做
#   ① K562 特异性闸门：基因扰动库里 K562 数据集上的得分 必须低（否则是泛化效应）
#   ② 无关程序对照：把**全部 39 个细胞型的疾病轴**各跑一遍，看候选在成纤维轴上是不是真的更靠前
#   ③ 多零模型 + BH 校正：三个零模型（打乱基因标签 / 随机高斯轴 / 按表达分位打乱），
#      要求 ≥2/3 显著且 BH 校正后 q<0.05
import os, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import h5py, numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import rankdata
ROOT = "/home/eto/luad_v2"; SC = f"{ROOT}/results/10_niche/tr_singlecell"
SEED = 20261005; N_PERM = 500; UP_DN = 200
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ── ① 基因扰动库 ──
f = h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r")
def cat(p):
    g = f[p]; c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]); return c[np.asarray(g["codes"][:])]
G = cat("var/gene_name"); X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
lab = pd.DataFrame({"gene": cat("obs/perturbed_gene_name"), "sign": np.asarray(f["obs/perturbation_sign"][:]),
                    "dataset": cat("obs/dataset")})
f.close()
Xnrm = np.linalg.norm(X.astype(np.float64), axis=1) + 1e-12
IS_K562 = lab.dataset.str.contains("K562", case=False).values
log(f"基因扰动库 {X.shape[0]:,} × {X.shape[1]:,}；K562 {int(IS_K562.sum()):,}")

# ── ② 39 个细胞型的疾病轴 ──
D = pd.read_parquet(f"{SC}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as h:
    sg = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["genes"][:]])
    M = sparse.csr_matrix((h["counts/data"][:], h["counts/indices"][:], h["counts/indptr"][:]),
                          shape=(len(D), len(sg)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
pid = D.patient_id.values; l2 = D.L2.values
PRE = np.isin(D.stage.values, ["AAH", "AIS", "MIA"]); IAC = D.stage.values == "IAC"
def axis_of(m):
    ps = []
    for p in np.unique(pid):
        ip = (pid == p) & m; a = ip & IAC; b = ip & PRE
        if a.sum() >= 20 and b.sum() >= 20:
            ps.append(np.asarray(Mn[a].mean(0)).ravel() - np.asarray(Mn[b].mean(0)).ravel())
    return np.median(np.vstack(ps), axis=0) if len(ps) >= 8 else None
AX = {t: axis_of(l2 == t) for t in np.unique(l2)}
AX = {k: v for k, v in AX.items() if v is not None}
fib_expr = np.asarray(Mn[l2 == "Fibroblast"].mean(0)).ravel()
log(f"疾病轴：{len(AX)} 个细胞型")
del M, Mn

pos = {g: i for i, g in enumerate(G)}
keep = np.array([i for i, g in enumerate(sg) if g in pos]); tgt = np.array([pos[sg[i]] for i in keep])
Xk = X[:, tgt].astype(np.float64); Xkn = Xnrm                     # 已按 tgt 重算范数
Xkn = np.linalg.norm(Xk, axis=1) + 1e-12
expr_k = fib_expr[keep]
log(f"共用基因 {len(keep):,}")

def cos_axis(dvec):
    d = dvec - dvec.mean(); dn = np.linalg.norm(d) + 1e-12
    return -((Xk @ d) / (Xkn * dn))          # 负号：越大 = 越把状态推回正常

# ── ③ 全部 39 个细胞型各跑一遍（无关程序对照）──
R = pd.DataFrame({"gene": lab.gene.values, "sign": lab.sign.values, "dataset": lab.dataset.values})
for t, ax in AX.items():
    R["ax::" + t] = cos_axis(ax[keep])
cols = [c for c in R.columns if c.startswith("ax::")]
GAGG = R.groupby(["gene", "sign"])[cols].median()
fib_col = "ax::Fibroblast"
GAGG["rank_fib"] = rankdata(-GAGG[fib_col])          # 1 = 最靠前
GAGG["rank_other_med"] = np.median([rankdata(-GAGG[c]) for c in cols if c != fib_col], axis=0)
GAGG["rank_gain"] = GAGG.rank_other_med - GAGG.rank_fib   # >0 = 在成纤维轴上更靠前
GAGG["n_axes_top1pct"] = (GAGG[cols].rank(ascending=False) <= len(GAGG) * 0.01).sum(axis=1)
CNT = R.groupby(["gene","sign"]).size()
GAGG = GAGG[GAGG.index.map(CNT) >= 3]
log(f"轴对照完成：{len(GAGG):,} 条（基因×方向）")

# ── ④ K562 闸门 ──
k5 = R[IS_K562].groupby(["gene", "sign"])[fib_col].median().rename("cos_k562")
nk = R[~IS_K562].groupby(["gene", "sign"])[fib_col].median().rename("cos_nonK562")
GAGG = GAGG.join(k5).join(nk)

# ── ⑤ 多零模型 + BH ──
d0 = AX["Fibroblast"][keep] - AX["Fibroblast"][keep].mean()
rng = np.random.default_rng(SEED)
obs = cos_axis(AX["Fibroblast"][keep])
expr_rank = rankdata(expr_k)
nbins = 20
p_mats = {}
def per_gene_p(null_cos):
    return np.array([(1 + (null_cos[i] >= obs[i]).sum()) / (1 + null_cos.shape[1])
                     for i in range(len(obs))])
# N1 打乱基因标签
NC = np.empty((len(obs), N_PERM), dtype=np.float32)
dd = d0.copy()
for k in range(N_PERM):
    rng.shuffle(dd); NC[:, k] = cos_axis(dd).astype(np.float32)
p_mats["N1打乱基因标签"] = per_gene_p(NC); log("  零模型 N1 完成")
# N2 随机高斯轴
NC = np.empty((len(obs), N_PERM), dtype=np.float32)
for k in range(N_PERM):
    NC[:, k] = cos_axis(rng.standard_normal(len(d0))).astype(np.float32)
p_mats["N2随机高斯轴"] = per_gene_p(NC); log("  零模型 N2 完成")
# N3 按表达分位打乱
binid = pd.qcut(expr_rank, nbins, labels=False)
NC = np.empty((len(obs), N_PERM), dtype=np.float32)
for k in range(N_PERM):
    dd = d0.copy()
    for b in range(nbins):
        idx = np.where(binid == b)[0]
        dd[idx] = rng.permutation(dd[idx])
    NC[:, k] = cos_axis(dd).astype(np.float32)
p_mats["N3按表达分位打乱"] = per_gene_p(NC); log("  零模型 N3 完成")

Pg = pd.DataFrame({k: pd.Series(v) for k, v in p_mats.items()})
Pg["gene"] = lab.gene.values; Pg["sign"] = lab.sign.values
Pg = Pg.groupby(["gene", "sign"])[list(p_mats)].median()
def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)
for k in p_mats: Pg["q_" + k] = bh(Pg[k].values)
Pg["q_med"] = Pg[[c for c in Pg.columns if c.startswith("q_")]].median(axis=1)
Pg["n_sig"] = (Pg[[c for c in Pg.columns if c.startswith("q_")]] < 0.05).sum(axis=1)
GAGG = GAGG.join(Pg[["q_med", "n_sig"]])
log(f"三零模型完成；BH 后 ≥2/3 显著的条目 {(GAGG.n_sig>=2).sum():,}/{len(GAGG):,}")

# ── ⑥ 合表 ──
GAGG = GAGG.reset_index()
prev = pd.read_csv(f"{SC}/CANDIDATE_TARGETS_v2.tsv", sep="\t")
prev["gene"] = prev.gene.str.upper(); GAGG["gene"] = GAGG.gene.str.upper()
for c in ["cos_k562","cos_nonK562"]:            # v2 里已有同名列，先删掉避免 merge 出 _x/_y
    if c in prev.columns: prev = prev.drop(columns=[c])
OUT = prev.merge(GAGG[["gene", "sign", "cos_k562", "cos_nonK562", "rank_fib", "rank_other_med",
                       "rank_gain", "n_axes_top1pct", "q_med", "n_sig"]],
                 on=["gene", "sign"], how="left")
def gate(r):
    if pd.isna(r.cos_k562):   k = "K562未测"
    elif r.cos_k562 > 0.05:   k = "❌泛化"
    else:                     k = "✅K562特异"
    s = "✅轴特异" if (pd.notna(r.rank_gain) and r.rank_gain > 0) else ("❌轴非特异" if pd.notna(r.rank_gain) else "轴未测")
    n = f"✅{int(r.n_sig)}/3零模型" if pd.notna(r.n_sig) and r.n_sig >= 2 else "❌零模型"
    return f"{k} | {s} | {n}"
OUT["闸门"] = OUT.apply(gate, axis=1)
OUT = OUT.sort_values(["final_tier", "cos_state"], ascending=[True, False])
OUT.round(4).to_csv(f"{SC}/CANDIDATE_TARGETS_v3.tsv", sep="\t", index=False)
log(f"落盘 CANDIDATE_TARGETS_v3.tsv（{len(OUT):,} 条）")
print("\n=== 闸门通过分布 ===")
print(OUT.闸门.value_counts().head(12).to_string())
print("\n=== T1/T2 且三闸全过 ===")
sel = OUT[(OUT.final_tier.str.startswith("T1")) | (OUT.final_tier.str.startswith("T2"))]
sel = sel[sel.闸门.str.contains("✅K562特异") & sel.闸门.str.contains("✅轴特异") & sel.闸门.str.contains("✅")]
print(sel[["final_tier", "gene", "dir_human", "cos_state", "cos_domain", "cos_k562",
           "rank_gain", "n_sig"]].round(4).to_string(index=False))
