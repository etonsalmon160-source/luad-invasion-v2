#!/usr/bin/env python3
# 65_global_decomposed.py —— **全局状态逆转 · 逐细胞类型分解**
#
# 用户口径：「我要全局状态，不是一个细胞型，你当然一个个分出来更好」
# 做法：对**每一个 L2 细胞型**各算一条**患者内配对**的疾病轴（IAC − 前驱，只在该型内比），
#       各跑一次**官方 CausalGenePredictor**，再把这些结果**合并成全局榜**。
# 好处：① 覆盖全组织（全局）② 每个轴内比⇒无细胞组成污染 ③ 能看出靶点是"全局共通"还是"某型特有"
import sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy import sparse
sys.path.insert(0, "/home/eto/scmg_workspace/SCMG")
from scmg.model.causal_prediction import CausalGenePredictor
import anndata as ad
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
TOP_N = 100                      # 每个细胞型取 top N 作为"该型命中的靶点"
MIN_PAT = 8                      # 至少这么多患者两期都有
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
def cat(fh, p):
    o = fh[p]
    if isinstance(o, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]]); return c[np.asarray(o["codes"][:])]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o[:]])

# ① 官方库 + 官方 stds
f = h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r")
hid = cat(f, "var/human_id"); Gp = cat(f, "var/gene_name"); X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
pe = cat(f, "obs/perturbed_gene"); pn = cat(f, "obs/perturbed_gene_name")
sg_ = np.asarray(f["obs/perturbation_sign"][:]); ds = cat(f, "obs/dataset"); f.close()
A = ad.AnnData(X=X, obs=pd.DataFrame({"perturbed_gene": pe, "perturbed_gene_name": pn,
                                      "perturbation_sign": sg_, "dataset": ds}),
               var=pd.DataFrame(index=pd.Index(hid)))
std = pd.read_csv(f"{OUT}/scmg_official_gene_stds.tsv", sep="\t")["std"].values.astype(np.float64)
cgp = CausalGenePredictor(A, std)
gp = {g: i for i, g in enumerate(Gp)}
log(f"官方库 {A.shape}；stds 就绪")

# ② 我们的细胞 + 预后 HR（质控用）
D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
HR = pd.read_csv(f"{OUT}/tcga_pergene_HR.tsv", sep="\t"); hrd = dict(zip(HR.gene, HR.HR))
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as h:
    sgenes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["genes"][:]])
    M = sparse.csr_matrix((h["counts/data"][:], h["counts/indices"][:], h["counts/indptr"][:]),
                          shape=(len(D), len(sgenes)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
pid = D.patient_id.values; l2 = D.L2.values
PRE = np.isin(D.stage.values, ["AAH", "AIS", "MIA"]); IAC = D.stage.values == "IAC"
log(f"细胞 {Mn.shape[0]:,}；L2 类型 {len(set(l2))}")

rows, used = [], []
for t in sorted(set(l2)):
    m = l2 == t
    ps = []
    for p in np.unique(pid):
        a = (pid == p) & m & IAC; b = (pid == p) & m & PRE
        if a.sum() >= 20 and b.sum() >= 20:
            ps.append(np.asarray(Mn[a].mean(0)).ravel() - np.asarray(Mn[b].mean(0)).ravel())
    if len(ps) < MIN_PAT:
        continue
    axis = np.median(np.vstack(ps), axis=0)
    v = np.zeros(len(Gp))
    for k, g in enumerate(sgenes):
        j = gp.get(g)
        if j is not None: v[j] = axis[k]
    r = cgp.calc_causal_scores(v)
    r = r.sort_values("causal_score", ascending=False).drop_duplicates("perturbed_gene", keep="first")
    r["cell_type"] = t; r["n_pat"] = len(ps)
    rows.append(r); used.append(t)
    log(f"  {t:28s} {len(ps):2d} 例 | top: {', '.join(r.perturbed_gene_name.head(3))}")
ALL = pd.concat(rows, ignore_index=True)
ALL.to_csv(f"{OUT}/PERCELLTYPE_causal.tsv", sep="\t", index=False)

# ③ 合并成全局榜：每个 (基因, 方向) 出现在多少个细胞型的 top-N 里
res = []
for (g, s), grp in ALL.groupby(["perturbed_gene_name", "perturbation_sign"]):
    gs = grp.gene_shift_z.median()
    tops = []
    for t, gt in grp.groupby("cell_type"):
        gt2 = gt.sort_values("causal_score", ascending=False).reset_index(drop=True)
        pos = gt2.index[gt2.perturbed_gene_name == g]
        if len(pos) and pos[0] < TOP_N: tops.append(t)
    if not tops: continue
    res.append({"gene": g, "sign": s, "intent": "敲低" if s < 0 else "过表达",
                "n_celltypes_top": len(tops), "n_celltypes_total": len(used),
                "frac": len(tops) / len(used),
                "median_causal": grp.causal_score.median(),
                "median_shift_z": gs, "celltypes": ",".join(sorted(tops))})
G = pd.DataFrame(res).sort_values(["n_celltypes_top", "median_causal"], ascending=False)
# ④ 预后质控
G["HR"] = G.gene.map(hrd)
ok = ((G.intent == "敲低") & (G.HR > 1)) | ((G.intent == "过表达") & (G.HR < 1))
G["qc_keep"] = ok.fillna(False)
G.to_csv(f"{OUT}/GLOBAL_DECOMPOSED_targets.tsv", sep="\t", index=False)
log(f"跑了 {len(used)} 个细胞型；全局榜 {len(G):,} 条（质控后 {int(G.qc_keep.sum()):,}）")

print(f"\n{'='*96}\n全局榜（按'在多少个细胞型里都进 top{TOP_N}'排）· 质控通过\n{'='*96}")
print(G[G.qc_keep].head(28)[["gene", "intent", "n_celltypes_top", "n_celltypes_total", "median_causal", "HR"]].round(3).to_string(index=False))
print(f"\n{'='*96}\n被预后质控剔除的（前 12）\n{'='*96}")
print(G[~G.qc_keep].head(12)[["gene", "intent", "n_celltypes_top", "median_shift_z", "HR"]].round(3).to_string(index=False))
