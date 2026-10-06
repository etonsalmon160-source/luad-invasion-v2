#!/usr/bin/env python3
# 56_scmg_native_gates.py —— SCMG 作者的 causal_score：**先与官方实现逐位对拍，再跑三道闸门**
#
# 官方公式（scmg/model/causal_prediction.py，逐行）：
#   mean_shift_z = shift / gene_stds
#   pert_vec_z[i] = X[i] / gene_stds
#   pms[i] = 1 - cosine_distance(mean_shift_z, pert_vec_z[i])          # = 余弦相似度
#   gene_shift_z[i] = clip(mean_shift_z[被扰动基因], -5, 5)
#   causal_score[i] = pms[i] * perturbation_sign[i] * gene_shift_z[i]
# 本脚本把它向量化（一次 matvec），**先对拍**，再跑：
#   闸门① K562 特异   闸门② 31 条细胞型轴特异   闸门③ 多零模型 + BH
import os, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy import sparse
from scipy.stats import rankdata
sys.path.insert(0, "/home/eto/scmg_workspace/SCMG")
ROOT = "/home/eto/luad_v2"; SC = f"{ROOT}/results/10_niche/tr_singlecell"
SEED = 20261005; N_PERM = 500
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ── 数据 ──
f = h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r")
def cat(p):                                  # 兼容 categorical(categories+codes) 与普通字符串 dataset
    o = f[p]
    if isinstance(o, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]])
        return c[np.asarray(o["codes"][:])]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o[:]])
Gp = cat("var/gene_name"); X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
hid = cat("var/human_id")                       # 🔴 var.index 用 Ensembl（SCMG 自己管道的用法）
pert_ens = cat("obs/perturbed_gene"); pert_gene = cat("obs/perturbed_gene_name")
sign = np.asarray(f["obs/perturbation_sign"][:])
ds = cat("obs/dataset"); f.close()
isf = {g: i for i, g in enumerate(Gp)}
ihid = {g: i for i, g in enumerate(hid)}
keep_p = np.array([g in ihid for g in pert_ens])                # 🔴 与官方一致：按 Ensembl 过滤（否则会错位一位）
log(f"扰动库 {X.shape}；被扰动基因在 var 内 {int(keep_p.sum()):,}/{len(keep_p):,}")

D = pd.read_parquet(f"{SC}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as h:
    sg = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["genes"][:]])
    M = sparse.csr_matrix((h["counts/data"][:], h["counts/indices"][:], h["counts/indptr"][:]),
                          shape=(len(D), len(sg)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
E1 = np.asarray(Mn.mean(0)).ravel(); E2 = np.asarray(Mn.power(2).mean(0)).ravel()
std_all = np.sqrt(np.maximum(E2 - E1 ** 2, 0))
pid = D.patient_id.values; l2 = D.L2.values
PRE = np.isin(D.stage.values, ["AAH", "AIS", "MIA"]); IAC = D.stage.values == "IAC"
fib = (l2 == "Fibroblast")
def axis_of(m):
    ps = []
    for p in np.unique(pid):
        a = (pid == p) & m & IAC; b = (pid == p) & m & PRE
        if a.sum() >= 20 and b.sum() >= 20:
            ps.append(np.asarray(Mn[a].mean(0)).ravel() - np.asarray(Mn[b].mean(0)).ravel())
    return np.median(np.vstack(ps), axis=0) if len(ps) >= 8 else None
AX = {t: axis_of(l2 == t) for t in np.unique(l2)}; AX = {k: v for k, v in AX.items() if v is not None}
del M, Mn
log(f"疾病轴 {len(AX)} 条")

pos = {g: i for i, g in enumerate(sg)}
keep = np.array([i for i, g in enumerate(sg) if g in isf])
tgt = np.array([isf[sg[i]] for i in keep])   # 🔴 映射到**扰动库**的下标（原写成 pos[...] = keep，静默错位）
# 🔴 gene_stds 改用 **SCMG 官方流形**算出的（58_ 存盘）；自算版会毁掉方法（见 §9.15）
STD_F = f"{SC}/scmg_official_gene_stds.tsv"
std_p = pd.read_csv(STD_F, sep="\t")["std"].values.astype(np.float64)   # 🔴 不能写 .std（那是 DataFrame 方法）
assert len(std_p) == len(Gp), (len(std_p), len(Gp))
log(f"gene_stds：官方流形 {int((std_p>0.1).sum()):,}/{len(std_p):,} 有真值")
Xk = X[:, tgt].astype(np.float64); std_t = std_p[tgt].astype(np.float64)
Xfull = X.astype(np.float64); stdf = std_p.astype(np.float64)      # 🔴 官方在完整 var 空间算，且是 float64
Pz = Xfull / stdf; Pzf = Pz                                         # 完整空间（对拍用）
Pz = Xk / std_t                                                     # 共有基因空间（闸门用，避免结构性零）
Pn = np.linalg.norm(Pz, axis=1) + 1e-12
gidx = np.array([isf.get(g, -1) for g in pert_gene])              # 被扰动基因在 var(18,108) 里的下标
_k_of_v = -np.ones(len(Gp), dtype=np.int64); _k_of_v[tgt] = np.arange(len(keep))
gidx_k = np.where(gidx >= 0, _k_of_v[np.clip(gidx, 0, None)], -1)  # 再映射到共有基因空间

def causal_score(axis_vec, full=False):
    """full=True：在完整 var 空间(18,108)算（与官方同口径）；
       full=False：只在共有基因空间算（闸门用，避免结构性零）。"""
    if full:
        az = np.zeros(len(Gp), dtype=np.float64); az[tgt] = axis_vec[keep]
        dz = az / stdf                                  # 🔴 只除一次标准差
        Pzv = Pzf; gi = gidx
    else:
        dz = axis_vec[keep] / std_t
        Pzv = Pz;  gi = gidx_k
    dn = np.linalg.norm(dz) + 1e-12
    pms = (Pzv @ dz) / ((np.linalg.norm(Pzv, axis=1) + 1e-12) * dn)
    gs = np.zeros(len(pms), dtype=np.float32)
    ok = gi >= 0
    gs[ok] = np.clip(dz[gi[ok]], -5, 5)
    return pms, gs, pms * sign * gs

# ── 对拍：与官方实现比 ──
from scmg.model.causal_prediction import CausalGenePredictor
import anndata as ad
obs = pd.DataFrame({"perturbed_gene": pert_ens, "perturbed_gene_name": pert_gene,
                    "perturbation_sign": sign, "dataset": ds})
adata_pert = ad.AnnData(X=X.astype(np.float32), obs=obs, var=pd.DataFrame(index=pd.Index(hid)))
cgp = CausalGenePredictor(adata_pert, std_p)
axis_p = np.zeros(len(Gp), dtype=np.float64)            # 🔴 官方实现要 var 空间(18,108)的向量
axis_p[np.array([isf[g] for g in sg[keep]])] = AX["Fibroblast"][keep]
ref = cgp.calc_causal_scores(axis_p)                    # 官方实现（内部已过滤 keep_p）
pms0, gs0, cs0 = causal_score(AX["Fibroblast"], full=True)   # 自实现（完整空间，与官方同口径）
assert len(ref) == int(keep_p.sum()), (len(ref), int(keep_p.sum()))
assert len(ref) > 19000
from scipy.stats import spearmanr
d2 = np.abs(ref.pert_sim.values - pms0[keep_p])
rho = spearmanr(ref.pert_sim.values, pms0[keep_p]).correlation
sg = np.mean(np.sign(ref.pert_sim.values) == np.sign(pms0[keep_p]))
log(f"对拍（{len(ref):,} 条）：pert_sim 中位差 {np.median(d2):.2e} / 最大差 {d2.max():.3e}")
log(f"对拍：Spearman = {rho:.8f}；符号一致 {100*sg:.4f}%")
assert rho > 0.9999 and sg > 0.999, "自实现与官方不一致 ⇒ 中止"
log("✅ 自实现与 SCMG 官方在秩与符号上完全一致")

# ── 闸门② 31 轴 ──
R = pd.DataFrame({"gene": pert_gene, "sign": sign, "dataset": ds})
for t, ax in AX.items():
    pm, gs, cs = causal_score(ax)
    R["rev::" + t] = pm * sign                     # 逆转分量（去掉 gene_shift_z，那是机制权重不是逆转）
cols = [c for c in R.columns if c.startswith("rev::")]
AGG = R.groupby(["gene", "sign"])[cols].median()
fibc = "rev::Fibroblast"
AGG["rank_fib"] = rankdata(-AGG[fibc])
AGG["rank_other_med"] = np.median([rankdata(-AGG[c]) for c in cols if c != fibc], axis=0)
AGG["rank_gain"] = AGG.rank_other_med - AGG.rank_fib

# ── 闸门① K562 ──
k5 = R[R.dataset.str.contains("K562", case=False)].groupby(["gene", "sign"])[fibc].median().rename("rev_k562")
k5n = R[~R.dataset.str.contains("K562", case=False)].groupby(["gene", "sign"])[fibc].median().rename("rev_nonK562")
AGG = AGG.join(k5).join(k5n)

# ── 闸门③ 多零模型 ──
pms0, gs0, cs0 = causal_score(AX["Fibroblast"])
rng = np.random.default_rng(SEED)
d0 = AX["Fibroblast"][keep].astype(np.float64) - AX["Fibroblast"][keep].mean()
expr_rank = rankdata(E1[keep]); nb = 20; binid = pd.qcut(expr_rank, nb, labels=False)
def pvec(null):     # null: (n_pert, N_PERM) 的 rev 值
    obs = pms0 * sign
    return np.array([(1 + (null[i] >= obs[i]).sum()) / (1 + null.shape[1]) for i in range(null.shape[0])])
PS = {}
def rev_of(dvec):
    dz = dvec / std_t                              # 🔴 dvec 已是共有基因空间，不再切片
    dn = np.linalg.norm(dz) + 1e-12
    return ((Pz @ dz) / (Pn * dn)) * sign
for nm, gen in [("N1打乱基因标签", None), ("N2随机高斯轴", None), ("N3按表达分位打乱", None)]:
    NC = np.empty((len(Pz), N_PERM), dtype=np.float32)
    for k in range(N_PERM):
        if nm == "N2随机高斯轴": dd = rng.standard_normal(len(d0))
        else:
            dd = d0.copy()
            if nm == "N1打乱基因标签": rng.shuffle(dd)
            else:
                for b in range(nb):
                    ix = np.where(binid == b)[0]; dd[ix] = rng.permutation(dd[ix])
        NC[:, k] = rev_of(dd)
    PS[nm] = pvec(NC); log(f"  {nm} 完成")
Pg = pd.DataFrame({k: v for k, v in PS.items()}); Pg["gene"] = pert_gene; Pg["sign"] = sign
Pg = Pg.groupby(["gene", "sign"])[list(PS)].median()
def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]; return np.clip(q, 0, 1)
for k in PS: Pg["q_" + k] = bh(Pg[k].values)
Pg["n_sig"] = (Pg[[c for c in Pg.columns if c.startswith("q_")]] < 0.05).sum(axis=1)
AGG = AGG.join(Pg[["n_sig"]])
# 主榜：gene_shift_z>0 与 <0 分开（作者规范）
gsv = pd.Series(gs0, index=pd.MultiIndex.from_arrays([pert_gene, sign], names=["gene", "sign"])).groupby(level=[0, 1]).median().rename("gene_shift_z")
AGG = AGG.join(gsv)
AGG = AGG.reset_index()
d0df = pd.DataFrame({"gene": pert_gene, "sign": sign, "causal_score": cs0, "pms": pms0})
d0df = d0df.groupby(["gene", "sign"]).agg(causal_score=("causal_score", "median"), pms=("pms", "median")).reset_index()
d0df["gene"] = d0df.gene.str.upper(); AGG["gene"] = AGG.gene.str.upper()
OUT = AGG.merge(d0df, on=["gene", "sign"], how="left")
cnt = R.groupby(["gene", "sign"]).size().rename("n_ds"); OUT = OUT.join(cnt, on=["gene", "sign"])
OUT = OUT[OUT.n_ds >= 3]
OUT.to_csv(f"{SC}/SCMG_NATIVE_GATES.tsv", sep="\t", index=False)
log(f"落盘 SCMG_NATIVE_GATES.tsv（{len(OUT):,} 条）")
for nm, m in [("gene_shift_z>0（疾病升高，敲低逆转）", OUT.gene_shift_z > 0),
              ("gene_shift_z<0（疾病降低，过表达逆转）", OUT.gene_shift_z < 0)]:
    s = OUT[m & (OUT.rev_k562 <= 0.02) & (OUT.rank_gain > 0) & (OUT.n_sig >= 2)]
    s = s.sort_values("causal_score", ascending=False)
    print(f"\n=== {nm}：三闸全过 {len(s)} 个 · top20 ===")
    print(s.head(20)[["gene", "n_ds", "causal_score", "pms", "rev_k562", "rank_fib", "rank_gain", "n_sig"]].round(4).to_string(index=False))
