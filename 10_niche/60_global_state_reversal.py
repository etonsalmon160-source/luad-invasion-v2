#!/usr/bin/env python3
# 60_global_state_reversal.py —— **全局状态逆转**（不是某一个细胞程序）
#
# 与 59_ 的唯一区别：位移不再取自"某个细胞型→某个细胞型"，
# 而是取自 **全部细胞**： mean(IAC 全部细胞) − mean(前驱期全部细胞)
# ⇒ 这是"整块组织的状态"往哪走，不是"成纤维自己怎么变"。
#
# 方法仍然 100% 是官方的：SCMG 的 CausalGenePredictor + 他们流形的 gene_stds。
import sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy import sparse
sys.path.insert(0, "/home/eto/scmg_workspace/SCMG")
from scmg.model.causal_prediction import CausalGenePredictor
import anndata as ad
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
def cat(fh, p):
    o = fh[p]
    if isinstance(o, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]]); return c[np.asarray(o["codes"][:])]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o[:]])

# ① 官方扰动库 + 官方 gene_stds
f = h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r")
hid = cat(f, "var/human_id"); X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
pe = cat(f, "obs/perturbed_gene"); pn = cat(f, "obs/perturbed_gene_name")
sg_ = np.asarray(f["obs/perturbation_sign"][:]); ds = cat(f, "obs/dataset"); f.close()
obs = pd.DataFrame({"perturbed_gene": pe, "perturbed_gene_name": pn,
                    "perturbation_sign": sg_, "dataset": ds})
A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(hid)))
std_official = pd.read_csv(f"{OUT}/scmg_official_gene_stds.tsv", sep="\t")["std"].values.astype(np.float64)
cgp = CausalGenePredictor(A, std_official)
log(f"库 {A.shape}；官方 gene_stds 就绪")

# ② 我们自己的细胞（全部，不是只取某一型）
D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as h:
    sgenes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["genes"][:]])
    M = sparse.csr_matrix((h["counts/data"][:], h["counts/indices"][:], h["counts/indptr"][:]),
                          shape=(len(D), len(sgenes)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
pid = D.patient_id.values; stg = D.stage.values
PRE = np.isin(stg, ["AAH", "AIS", "MIA"]); IAC = stg == "IAC"
log(f"我们的细胞 {Mn.shape[0]:,}（IAC {int(IAC.sum()):,} / 前驱 {int(PRE.sum()):,}）")

# ③ 全局位移：全部细胞（两种估计）
gA = np.asarray(Mn[IAC].mean(0)).ravel() - np.asarray(Mn[PRE].mean(0)).ravel()      # 全队列
ps = []
for p in np.unique(pid):
    a = (pid == p) & IAC; b = (pid == p) & PRE
    if a.sum() >= 50 and b.sum() >= 50:
        ps.append(np.asarray(Mn[a].mean(0)).ravel() - np.asarray(Mn[b].mean(0)).ravel())
gB = np.median(np.vstack(ps), axis=0)                                               # 患者内配对
log(f"全局位移：全队列版 / 患者配对版（用了 {len(ps)} 例）")
del M, Mn

# ④ 映射到官方 var（Ensembl），用 symbol→gene_name 对应
gp_pos = {g: i for i, g in enumerate(cat(h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r"), "var/gene_name"))}
for nm, sh in [("全队列", gA), ("患者配对", gB)]:
    v = np.zeros(len(hid))
    for k, g in enumerate(sgenes):
        j = gp_pos.get(g)
        if j is not None: v[j] = sh[k]
    r = cgp.calc_causal_scores(v)
    r = r.sort_values("causal_score", ascending=False).drop_duplicates("perturbed_gene", keep="first")
    r.to_csv(f"{OUT}/GLOBAL_reversal_{'all' if nm=='全队列' else 'paired'}.tsv", sep="\t", index=False)
    print(f"\n{'='*80}\n【{nm}】全局状态逆转 —— 靶点榜\n{'='*80}")
    print("① 疾病里升高的基因 → 敲低来逆转")
    print(r[r.gene_shift_z > 0].head(18)[["perturbed_gene_name", "causal_score", "pert_sim", "gene_shift_z"]].round(4).to_string(index=False))
    print("\n② 疾病里降低的基因 → 过表达来逆转")
    print(r[r.gene_shift_z < 0].head(18)[["perturbed_gene_name", "causal_score", "pert_sim", "gene_shift_z"]].round(4).to_string(index=False))
