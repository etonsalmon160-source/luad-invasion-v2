#!/usr/bin/env python3
# ============================================================================
# 🔴【已作废 · 2026-09-16】本脚本**请勿运行**，保留仅为审计留痕。
#   作废原因：本脚本实现的是旧 GP4a 口径（scran 池化 SF + log1p + HVG 2000 +
#     `sc.pp.scale(max_value=10)` + PCA 30），该口径已被用户决策「对齐论文全套」整体作废。
#     论文配方：SCTransform(vst.flavor="v2") → HVG 3000 → PCA 50 → Harmony → Louvain。
#   登记位置：docs/PARAMETERS_AND_SOURCES.md §M3-A.1「已作废登记」表。
#   现行替代：04_integration/10_seurat_traditional.R
#   ⚠️ 输入依赖 00a/00b，二者亦已作废 ⇒ 本脚本现在**也跑不通**（scran_size_factors.csv.gz 不会生成）。
# ============================================================================
# 01_preprocess_gse308103.py —— 传统分支预处理（GP4a）
#
# 输入：results/02_expression/gse308103_counts.h5ad          （GP0 产物，原始计数）
#       results/04_integration/scran_size_factors.csv.gz     （00b 产物，scran 池化 SF）
# 输出：results/04_integration/gse308103_preprocessed.h5ad
#       results/04_integration/{preprocess_manifest.json, hvg_2000.txt}
#
# ────────────────────────────────────────────────────────────────────────────
# 参数登记块（法则 3.1）—— 全部与 docs/PARAMETERS_AND_SOURCES.md M3-A.1 对齐
#
#   步              取值                              出处
#   ------------    ------------------------------    ------------------------------
#   去零基因        丢弃 count 全零的基因            ⚠️C 本项目约定（4 个基因，见 manifest）
#   归一化          **scran 池化 size factor**        Lun 2016, Genome Biol 17:75
#                   （由 00b_scran_sizefactors.R 算，此处读 CSV 应用）
#   变换            log1p                             Ahlmann-Eltze & Huber 2023,
#                                                     Nat Methods 20:665
#   HVG             2000, flavor='seurat_v3'           Stuart 2019, Cell 177:1888
#                   batch_key='sample_id'
#   ⚠️ **HVG 必须在原始计数上算**（seurat_v3 是方差稳定变换），且必须在归一化**之前**。
#   scale           max_value=10，仅 HVG 子集           🟡D scanpy 默认
#   PCA             n_comps=30, svd_solver='arpack'    Heumos 2023, Nat Rev Genet 24:550
#   kNN             n_neighbors=15, n_pcs=30           Heumos 2023
#   回归            **不做**（不回归 nCount / pct_mt）  Heumos 2023；Hafemeister & Satija 2019
#
# ⚠️ 明确不做：
#   · 不做 `var_names_make_unique`（只改名、会静默丢信号）—— 遇重名基因**硬停**
#   · 不做批次校正（Arm A）；Harmony 见 02 脚本（GP4c）
#   · 不做 `sc.pp.regress_out`
#   · **不落盘稠密 scale 矩阵**（648,945×2000 float64 = 10 GB；PCA 后即弃）
# ============================================================================
import json, time, hashlib, sys
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import pandas as pd

ROOT   = Path("/home/eto/luad_v2")
INP    = ROOT / "results/02_expression/gse308103_counts.h5ad"
SF_CSV = ROOT / "results/04_integration/scran_size_factors.csv.gz"
OUTDIR = ROOT / "results/04_integration"
H5     = OUTDIR / "gse308103_preprocessed.h5ad"
HVG_TX = OUTDIR / "hvg_2000.txt"
MAN    = OUTDIR / "preprocess_manifest.json"

N_TOP_GENES = 2000          # Stuart 2019
N_PCS       = 30            # Heumos 2023（10–50）
N_NEIGHBORS = 15            # Heumos 2023
SCALE_MAX   = 10.0          # 🟡D 常规
BATCH_KEY   = "sample_id"   # GP4c：样本是唯一真实技术批次轴

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

def sha256(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while (b := f.read(buf)):
            h.update(b)
    return h.hexdigest()

# ---- 1. 载入原始计数 ------------------------------------------------------
log(f"载入 {INP}")
adata = ad.read_h5ad(INP)
n_obs_in, n_var_in = adata.shape
log(f"  shape={adata.shape}  X={adata.X.dtype}  sparse={sp.issparse(adata.X)}")

X = adata.X
probe = X.data if sp.issparse(X) else X.ravel()
if not np.all(probe[:5_000_000] == np.floor(probe[:5_000_000])):
    sys.exit("[FAIL] 输入不是整数原始计数 —— seurat_v3 需要 raw counts，拒绝继续")
if adata.layers or adata.raw is not None or adata.obsm:
    sys.exit(f"[FAIL] 输入带归一化痕迹 layers={list(adata.layers)} "
             f"raw={adata.raw is not None} obsm={list(adata.obsm)}，拒绝继续")
log("前置断言通过：原始整数计数、无归一化痕迹")

# ---- 2. 基因名唯一性（硬停，不做 make_unique）------------------------------
dups = adata.var_names[adata.var_names.duplicated()].tolist()
if dups:
    sys.exit(f"[FAIL] 有 {len(dups)} 个重名基因，禁止 var_names_make_unique（会静默丢信号）：{dups[:10]}")
log("基因名唯一性通过")

# ---- 3. 去全零基因 --------------------------------------------------------
nz = np.asarray((X > 0).sum(axis=0)).ravel()
zero_genes = adata.var_names[nz == 0].tolist()
if zero_genes:
    log(f"丢弃 {len(zero_genes)} 个全零基因：{zero_genes}")
    adata = adata[:, nz > 0].copy()
log(f"  shape={adata.shape}")

# ---- 4. HVG：**在原始计数上**算（顺序关键，必须在归一化之前）---------------
log(f"HVG n_top_genes={N_TOP_GENES} flavor=seurat_v3 batch_key={BATCH_KEY}")
sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES,
                            flavor="seurat_v3", batch_key=BATCH_KEY)
n_hvg = int(adata.var.highly_variable.sum())
if n_hvg != N_TOP_GENES:
    sys.exit(f"[FAIL] HVG 选出 {n_hvg} 个，期望 {N_TOP_GENES}")
nbatches = adata.var.highly_variable_nbatches
log(f"  HVG={n_hvg}；逐批次命中分布 "
    f"{dict(zip(*np.unique(nbatches.values, return_counts=True)))}（共 {adata.obs[BATCH_KEY].nunique()} 批次）")
HVG_TX.write_text("\n".join(adata.var_names[adata.var.highly_variable]) + "\n")

# ---- 5. 保留原始计数层，再按 scran size factor 归一化 ---------------------
adata.layers["counts"] = adata.X.copy()

log(f"读入 scran size factor {SF_CSV}")
sfdf = pd.read_csv(SF_CSV)
if not {"cell_barcode", "size_factor", "sample_id"} <= set(sfdf.columns):
    sys.exit(f"[FAIL] scran CSV 缺列：{list(sfdf.columns)}")
if len(sfdf) != adata.n_obs:
    sys.exit(f"[FAIL] scran CSV {len(sfdf)} 行 vs h5ad {adata.n_obs} 细胞，拒绝继续")

sfdf = sfdf.set_index("cell_barcode").reindex(adata.obs_names)
n_missing = int(sfdf["size_factor"].isna().sum())
if n_missing:
    sys.exit(f"[FAIL] 有 {n_missing} 个细胞在 scran CSV 中缺失")
if not (sfdf["sample_id"].values == adata.obs[BATCH_KEY].values).all():
    sys.exit("[FAIL] scran CSV 的 sample_id 与 h5ad obs 不一致 —— 细胞顺序对不上")
sf = sfdf["size_factor"].to_numpy(np.float64)
if not (np.isfinite(sf).all() and (sf > 0).all()):
    sys.exit("[FAIL] size factor 含非法值")
log(f"  对齐通过（{adata.n_obs} 细胞）；SF 分位 min/med/max = "
    f"{np.round([sf.min(), np.median(sf), sf.max()], 4).tolist()}")

adata.obs["size_factor"] = sf.astype(np.float32)
adata.X = sp.diags(1.0 / sf) @ adata.X
sc.pp.log1p(adata)
log("counts/size_factor + log1p 完成")

# ---- 6. HVG 子集上 scale → PCA → kNN --------------------------------------
hv = adata.var.highly_variable.values
hvg = adata[:, hv].copy()
del adata.var.highly_variable
log(f"scale(max_value={SCALE_MAX}) on {hvg.shape}")
sc.pp.scale(hvg, max_value=SCALE_MAX, zero_center=True)

log(f"PCA n_comps={N_PCS} svd_solver=arpack")
sc.tl.pca(hvg, n_comps=N_PCS, svd_solver="arpack")
vr = hvg.uns["pca"]["variance_ratio"]
log(f"  前 5 PC 解释方差比 {np.round(vr[:5], 4).tolist()}；累计 {vr.sum():.4f}")

adata.obsm["X_pca"] = hvg.obsm["X_pca"]
adata.uns["pca"] = hvg.uns["pca"]
log(f"neighbors n_neighbors={N_NEIGHBORS} n_pcs={N_PCS}")
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep="X_pca")
del hvg

# ---- 7. 落盘 --------------------------------------------------------------
log(f"写出 {H5}")
adata.write_h5ad(H5, compression="lzf")
log(f"  {H5.stat().st_size/2**30:.2f} GiB；计算 SHA-256")
h5_sha = sha256(H5)

man = {
    "script": "04_integration/01_preprocess_gse308103.py",
    "checkpoint": "GP4a",
    "params": {
        "normalize": "scran pooled size factors (via 00b_scran_sizefactors.R)",
        "transform": "log1p", "n_top_genes": N_TOP_GENES, "hvg_flavor": "seurat_v3",
        "hvg_batch_key": BATCH_KEY, "hvg_on_raw_counts": True,
        "scale_max_value": SCALE_MAX, "n_comps": N_PCS, "svd_solver": "arpack",
        "n_neighbors": N_NEIGHBORS, "regress_out": None, "batch_correction": None,
    },
    "inputs": {
        "counts_h5ad": {"path": str(INP), "sha256": sha256(INP)},
        "scran_size_factors": {"path": str(SF_CSV), "sha256": sha256(SF_CSV)},
    },
    "counts": {"n_obs_in": int(n_obs_in), "n_var_in": int(n_var_in),
               "n_zero_genes_dropped": len(zero_genes), "zero_count_genes": zero_genes,
               "n_hvg": n_hvg},
    "size_factor": {"min": float(sf.min()), "median": float(np.median(sf)),
                    "max": float(sf.max()), "geometric_mean": float(np.exp(np.mean(np.log(sf))))},
    "pca_variance_ratio": [float(x) for x in vr],
    "pca_variance_ratio_cumsum": float(vr.sum()),
    "knn_graph_nnz": int(adata.obsp["distances"].nnz),
    "outputs": {
        "h5ad": {"path": str(H5), "sha256": h5_sha, "bytes": int(H5.stat().st_size),
                 "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
                 "nnz": int(adata.X.nnz), "compression": "lzf"},
        "hvg_list": {"path": str(HVG_TX), "sha256": sha256(HVG_TX), "n": n_hvg},
    },
    "versions": {"scanpy": sc.__version__, "anndata": ad.__version__,
                 "numpy": np.__version__, "scipy": scipy.__version__,
                 "pandas": pd.__version__},
    "wall_sec": round(time.time() - T0, 1),
}
MAN.write_text(json.dumps(man, indent=2, ensure_ascii=False))
log(f"manifest -> {MAN}")
log(f"完成：{adata.n_obs} 细胞 × {adata.n_vars} 基因，HVG {n_hvg}，PC {N_PCS}")
