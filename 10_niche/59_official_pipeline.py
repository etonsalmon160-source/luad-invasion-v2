#!/usr/bin/env python3
# 59_official_pipeline.py —— **照 SCMG 官方教程原样跑，不造轮子**
#
# 教程（predict_causal_genes_for_cell_state_transitions.ipynb）逐条：
#   1. 读参考流形 → CP10K+log1p
#   2. gene_stds = 每个基因在流形上的标准差
#   3. exp_shift = mean(目标细胞型) − mean(来源细胞型)
#   4. CausalGenePredictor(adata_pert, gene_stds).calc_causal_scores(exp_shift)
#   5. 按 causal_score 排序、按 perturbed_gene 去重
#   6. 按 gene_shift_z 正负**分成两张榜**
#
# 我们只是把「Epiblast → Nascent mesoderm」换成**肺的成纤维激活轴**（流形自带的细胞型）。
import sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy import sparse
sys.path.insert(0, "/home/eto/scmg_workspace/SCMG")
from scmg.model.causal_prediction import CausalGenePredictor
import anndata as ad
OUT = "/home/eto/luad_v2/results/10_niche/tr_singlecell"
MANI = "/home/eto/scmg_workspace/hf_data/ref_global_cell_state_manifold.h5ad"
PERT = "/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
def cat(fh, p):
    o = fh[p]
    if isinstance(o, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]]); return c[np.asarray(o["codes"][:])]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o[:]])

# ① 流形
h = h5py.File(MANI, "r")
gid = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["var/_index"][:]])
ct = cat(h, "obs/cell_type")
M = sparse.csr_matrix((h["X/data"][:], h["X/indices"][:], h["X/indptr"][:]),
                      shape=(len(ct), len(gid)))
h.close()
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
E1 = np.asarray(Mn.mean(0)).ravel(); E2 = np.asarray(Mn.power(2).mean(0)).ravel()
gene_stds = np.maximum(np.sqrt(np.maximum(E2 - E1 ** 2, 0)), 0.1)
log(f"流形 {M.shape[0]:,} × {M.shape[1]:,}；gene_stds 就绪")

# ② 扰动库（官方期望格式：var.index = Ensembl）
f = h5py.File(PERT, "r")
hid = cat(f, "var/human_id"); X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
pe = cat(f, "obs/perturbed_gene"); pn = cat(f, "obs/perturbed_gene_name")
sg = np.asarray(f["obs/perturbation_sign"][:]); ds = cat(f, "obs/dataset"); f.close()
obs = pd.DataFrame({"perturbed_gene": pe, "perturbed_gene_name": pn,
                    "perturbation_sign": sg, "dataset": ds})
A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(hid)))
cgp = CausalGenePredictor(A, gene_stds)

# ③ 肺的成纤维激活轴（流形自带细胞型）
PAIRS = [("pulmonary interstitial fibroblast", "myofibroblast cell", "肺间质成纤维 → 肌成纤维"),
         ("fibroblast of lung",              "myofibroblast cell", "肺成纤维 → 肌成纤维"),
         ("alveolar type 2 fibroblast cell", "myofibroblast cell", "肺泡成纤维 → 肌成纤维")]
for src, tgt, label in PAIRS:
    si = np.where(ct == src)[0]; ti = np.where(ct == tgt)[0]
    if len(si) < 20 or len(ti) < 20:
        log(f"跳过 {label}（细胞数 {len(si)}/{len(ti)}）"); continue
    shift = np.asarray(Mn[ti].mean(0)).ravel() - np.asarray(Mn[si].mean(0)).ravel()
    r = cgp.calc_causal_scores(shift)
    r = r.sort_values("causal_score", ascending=False).drop_duplicates("perturbed_gene", keep="first")
    tag = src.split()[0] + "_to_myofib"
    r.to_csv(f"{OUT}/official_{tag}.tsv", sep="\t", index=False)
    print(f"\n{'='*78}\n{label}（源 {len(si)} / 靶 {len(ti)} 细胞）\n{'='*78}")
    print(f"【疾病里升高的基因 → 敲低来逆转】top 15")
    print(r[r.gene_shift_z > 0].head(15)[["perturbed_gene_name", "causal_score", "pert_sim", "gene_shift_z"]].round(4).to_string(index=False))
    print(f"\n【疾病里降低的基因 → 过表达来逆转】top 15")
    print(r[r.gene_shift_z < 0].head(15)[["perturbed_gene_name", "causal_score", "pert_sim", "gene_shift_z"]].round(4).to_string(index=False))
