#!/usr/bin/env python3
# 55_scmg_native.py —— 用 **SCMG 作者自己的 CausalGenePredictor** 跑我们的数据（第四种口径）
#
# 与他们的官方教程逐条对齐（SCMG/docs/source/tutorials/predict_causal_genes_for_cell_state_transitions.ipynb）：
#   · exp_shift = 目标状态均值 − 来源状态均值（CP10K + log1p）
#   · gene_stds = 每个基因在参考流形上的标准差（教程里从 adata_ref 算；我们没有预计算表，用本项目数据算）
#   · causal_score = cos_sim × perturbation_sign × gene_shift_z
#   · 输出按 gene_shift_z 正负**分成两张榜**
# 两个查询版本：
#   A. 全队列平均差（教程口径：mean(IAC) − mean(前驱)）
#   B. 患者内配对中位（本项目 §9.14 的轴 B）—— 用来测"配对有没有影响"
import os, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy import sparse
sys.path.insert(0, "/home/eto/scmg_workspace/SCMG")
from scmg.model.causal_prediction import CausalGenePredictor
import anndata as ad
ROOT = "/home/eto/luad_v2"; SC = f"{ROOT}/results/10_niche/tr_singlecell"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ── ① 扰动库 → AnnData ──
H5 = "/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad"
f = h5py.File(H5, "r")
def cat(p):
    g = f[p]; c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]); return c[np.asarray(g["codes"][:])]
genes_p = cat("var/gene_name")
Xp = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
obs = pd.DataFrame({"perturbed_gene": cat("obs/perturbed_gene_name"),
                    "perturbed_gene_name": cat("obs/perturbed_gene_name"),
                    "perturbation_sign": np.asarray(f["obs/perturbation_sign"][:]),
                    "dataset": cat("obs/dataset")})
f.close()
var = pd.DataFrame(index=pd.Index(genes_p, name=None))
adata_pert = ad.AnnData(X=Xp, obs=obs, var=var)
log(f"扰动库 AnnData {adata_pert.shape}")

# ── ② 我们自己的细胞 → 基因标准差 + 两个查询向量 ──
D = pd.read_parquet(f"{SC}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as h:
    sg = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["genes"][:]])
    M = sparse.csr_matrix((h["counts/data"][:], h["counts/indices"][:], h["counts/indptr"][:]),
                          shape=(len(D), len(sg)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
log("算基因标准差（全体细胞）…")
E1 = np.asarray(Mn.mean(0)).ravel(); E2 = np.asarray(Mn.power(2).mean(0)).ravel()
std_all = np.sqrt(np.maximum(E2 - E1 ** 2, 0))
pid = D.patient_id.values; l2 = D.L2.values
PRE = np.isin(D.stage.values, ["AAH", "AIS", "MIA"]); IAC = D.stage.values == "IAC"
fib = (l2 == "Fibroblast")
# A：全队列平均差（教程口径）
shiftA = np.asarray(Mn[fib & IAC].mean(0)).ravel() - np.asarray(Mn[fib & PRE].mean(0)).ravel()
# B：患者内配对中位（我们的轴）
ps = []
for p in np.unique(pid):
    a = (pid == p) & fib & IAC; b = (pid == p) & fib & PRE
    if a.sum() >= 20 and b.sum() >= 20:
        ps.append(np.asarray(Mn[a].mean(0)).ravel() - np.asarray(Mn[b].mean(0)).ravel())
shiftB = np.median(np.vstack(ps), axis=0)
log(f"查询 A/B 就绪；配对用了 {len(ps)} 例患者")
del M, Mn

# ── ③ 对齐基因空间（按扰动库的 var 顺序排 std 与 shift）──
# 🔴🔴 本脚本**产物已作废**（见 TARGET_REVERSAL_STATUS.md §四）：
#    这里用的是**替代的 gene_stds**（缺失填 0.1），而官方口径要求用流形算出的真 std
#    ⇒ 产出 `scmg_native_A/B.tsv` 的因果分全在 1e-6 量级、全平，**「官方方法找不到」是假象**。
#    **现行版是 `56_scmg_native_gates.py`**（用官方 stds）。本脚本保留仅作审计，勿引用其产物。
#
# ⚠️ 审计补充（2026-10-06）：下面第 2 行 `np.maximum(std_p, 0.1)` 会把**真实 SD 也低于 0.1 的基因
#    抬升到 0.1**——这一步此前未单独登记，它会人为放大低方差基因。属已作废路径，不再修，
#    但任何复用这段代码的脚本都必须先删掉这一行。
pos = {g: i for i, g in enumerate(sg)}
std_p = np.array([std_all[pos[g]] if g in pos else np.nan for g in genes_p])
for nm, sh in [("A", shiftA), ("B", shiftB)]:
    v = np.array([sh[pos[g]] if g in pos else 0.0 for g in genes_p])
    globals()["shift_" + nm] = v
    log(f"  shift_{nm}: 非零 {int((v!=0).sum()):,}/{len(v):,}；std 缺 {int(np.isnan(std_p).sum()):,}")
std_p = np.nan_to_num(std_p, nan=0.1)
std_p = np.maximum(std_p, 0.1)

# ── ④ 跑作者的方法 ──
OUTS = {}
for nm in ["A", "B"]:
    cgp = CausalGenePredictor(adata_pert, std_p)
    log(f"跑 CausalGenePredictor（查询 {nm}）…")
    df = cgp.calc_causal_scores(globals()["shift_" + nm])
    df = df.sort_values("causal_score", ascending=False).drop_duplicates("perturbed_gene", keep="first")
    OUTS[nm] = df
    df.to_csv(f"{SC}/scmg_native_{nm}.tsv", sep="\t", index=False)
    up = df[df.gene_shift_z > 0].head(15); dn = df[df.gene_shift_z < 0].head(15)
    print(f"\n=== 查询 {nm}：gene_shift_z>0（在疾病里升高的基因）top15 ===")
    print(up[["perturbed_gene_name", "causal_score", "pert_sim", "gene_shift_z"]].round(4).to_string(index=False))
    print(f"=== 查询 {nm}：gene_shift_z<0（在疾病里降低的基因）top15 ===")
    print(dn[["perturbed_gene_name", "causal_score", "pert_sim", "gene_shift_z"]].round(4).to_string(index=False))

# ── ⑤ A vs B 一致性 ──
m = OUTS["A"][["perturbed_gene", "causal_score"]].merge(
    OUTS["B"][["perturbed_gene", "causal_score"]], on="perturbed_gene", suffixes=("_A", "_B"))
from scipy.stats import spearmanr
print(f"\nA(全队列) vs B(患者配对) 的 causal_score Spearman = {spearmanr(m.causal_score_A, m.causal_score_B).correlation:.3f}")
