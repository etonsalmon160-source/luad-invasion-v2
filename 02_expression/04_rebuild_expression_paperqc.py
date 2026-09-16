#!/usr/bin/env python3
# ============================================================================
# 04_rebuild_expression_paperqc.py —— 按论文 QC 分析集重建表达对象（GP0 重做）
#
# 背景：用户 2026-09-16 决策「改用论文 QC 重建」。旧的 648,945 核对象
#       （gse308103_counts.h5ad）**保留不删**，作为敏感性臂；本脚本产出新的
#       413,697 核对象，作为此后唯一分析集。
#
# 掩膜来源：results/01_qc/gse308103_analysis_mask_paperqc.csv.gz
#           = 论文计数门 ∩ scDblFinder 单细胞（见 01_qc/11_regate_paper_qc.py）
#
# 论文口径的两步（顺序照论文：先筛细胞，再筛基因）：
#   ① 细胞：nFeature>=500 & nCount>=1000 & pct_mt<=20   → 413,697 核
#   ② 基因：在**新细胞集上**检出数 <3 的基因剔除
#
# ⚠️ 基因名重名**硬停**，绝不用 var_names_make_unique（只改名、静默丢信号）。
# ⚠️ 不写归一化副本、不设 .raw —— 只存原始整数计数，下游 SCTransform 自己算。
#
# 输入：results/02_expression/gse308103_counts.h5ad            （旧 648,945 对象）
#       results/01_qc/gse308103_analysis_mask_paperqc.csv.gz   （新掩膜）
# 输出：results/02_expression/gse308103_counts_paperqc.h5ad
#       results/02_expression/rebuild_paperqc_manifest.json
#       results/04_integration/seurat_io/{data_f32,indices_i32,indptr_i32}.bin
#       results/04_integration/seurat_io/{gene_names.txt,cell_names.txt,cell_meta.csv.gz}
#       results/04_integration/seurat_io/export_manifest.json
# ============================================================================
import json, time, hashlib, sys
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse as sp
import anndata as ad
import pandas as pd

ROOT = Path("/home/eto/luad_v2")
OLD  = ROOT / "results/02_expression/gse308103_counts.h5ad"
MASK = ROOT / "results/01_qc/gse308103_analysis_mask_paperqc.csv.gz"
NEW  = ROOT / "results/02_expression/gse308103_counts_paperqc.h5ad"
MAN  = ROOT / "results/02_expression/rebuild_paperqc_manifest.json"
RIO  = ROOT / "results/04_integration/seurat_io"

EXPECT_N_CELLS = 413697      # 由 11_regate_paper_qc.py 实测得出；此处作为回归断言
GENE_MIN_CELLS = 3           # 论文：在 ≥3 个细胞中检出

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

def sha256(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while (b := f.read(buf)):
            h.update(b)
    return h.hexdigest()

# ---- 1. 读入旧对象与掩膜 ---------------------------------------------------
log(f"载入旧对象 {OLD.name}")
a = ad.read_h5ad(OLD)
n_old, n_var_old = a.shape
log(f"  {a.shape}  X={a.X.dtype}  nnz={a.X.nnz}")

if a.layers or a.raw is not None or a.obsm:
    sys.exit(f"[FAIL] 旧对象带归一化痕迹 layers={list(a.layers)} raw={a.raw is not None} "
             f"obsm={list(a.obsm)}，应只有原始计数")
if not sp.issparse(a.X):
    sys.exit("[FAIL] X 不是稀疏矩阵")

d = a.X.data
if not np.all(d[:5_000_000] == np.floor(d[:5_000_000])):
    sys.exit("[FAIL] 旧对象 X 非整数计数")

mk = pd.read_csv(MASK)
log(f"读入新掩膜 {len(mk)} 核")
if len(mk) != EXPECT_N_CELLS:
    sys.exit(f"[FAIL] 掩膜 {len(mk)} != 预期 {EXPECT_N_CELLS}")
if mk.cell_barcode.duplicated().any():
    sys.exit("[FAIL] 掩膜 cell_barcode 有重复")

# ---- 2. 细胞子集（保持旧对象原有行序，便于追溯）----------------------------
keep = pd.Index(a.obs_names).isin(set(mk.cell_barcode))
n_keep = int(keep.sum())
if n_keep != EXPECT_N_CELLS:
    sys.exit(f"[FAIL] 掩膜命中 {n_keep} != {EXPECT_N_CELLS} —— 掩膜与旧对象不同源")
log(f"掩膜命中 {n_keep} 核（占旧对象 {100*n_keep/n_old:.2f}%）")

a = a[keep].copy()
log(f"  子集后 {a.shape}  nnz={a.X.nnz}")

# obs 用掩膜口径回填 patient/stage（防旧对象里挂错）
mk_idx = mk.set_index("cell_barcode")
a.obs["patient_id"] = pd.Categorical(mk_idx.loc[a.obs_names, "patient_id"].values)
a.obs["stage"]      = pd.Categorical(mk_idx.loc[a.obs_names, "stage"].values)
a.obs["stage_token"]= pd.Categorical(mk_idx.loc[a.obs_names, "stage_token"].values)
n_bad = int((a.obs["sample_id"].astype(str).values
             != mk_idx.loc[a.obs_names, "sample_id"].values).sum())
if n_bad:
    sys.exit(f"[FAIL] {n_bad} 个核的 sample_id 在掩膜与旧对象间不一致")
log("  obs 回填通过（patient_id / stage / stage_token / sample_id 四项交叉核对）")

# ---- 3. 基因名唯一性（硬停，不做 make_unique）------------------------------
dups = a.var_names[a.var_names.duplicated()].tolist()
if dups:
    sys.exit(f"[FAIL] {len(dups)} 个重名基因，禁止 make_unique：{dups[:10]}")
log("基因名唯一性通过")

# ---- 4. 基因过滤：在新细胞集上检出 <3 的基因剔除（论文口径）----------------
det = np.asarray((a.X > 0).sum(axis=0)).ravel()
drop = det < GENE_MIN_CELLS
n_drop = int(drop.sum())
dropped_genes = a.var_names[drop].tolist()
log(f"检出 <{GENE_MIN_CELLS} 细胞的基因：{n_drop} 个 → 剔除")
if n_drop:
    log(f"  示例：{dropped_genes[:12]}")
a = a[:, ~drop].copy()
n_var_new = a.n_vars
log(f"  基因 {n_var_old} → {n_var_new}（-{n_drop}）")

# ---- 5. 写出新 h5ad --------------------------------------------------------
a.uns["build"] = {
    "note": "论文 QC 重建（nFeature>=500 & nCount>=1000 & pct_mt<=20；基因检出>=3 细胞）",
    "mask": str(MASK.name), "gate": "GP0-redo / 01_qc/11_regate_paper_qc.py",
    "doublet_caller": "scDblFinder（铁律 R2；论文用 Scrublet，见登记偏差）",
}
log(f"写出 {NEW.name}")
a.write_h5ad(NEW, compression="lzf")
log(f"  {NEW.stat().st_size/2**30:.2f} GiB；计算 SHA-256")
new_sha = sha256(NEW)

# ---- 6. 导出 R（Seurat）可读的 genes×cells CSC 二进制 ----------------------
RIO.mkdir(parents=True, exist_ok=True)
log("导出 Seurat I/O（genes×cells CSC 三元组）")
T = a.X.T.tocsc()                       # CSR 的 .T 是 CSC 视图，不复制
assert T.shape == (n_var_new, a.n_obs)
data, indices, indptr = T.data, T.indices, T.indptr
assert data.dtype == np.float32 and indices.dtype == np.int32
assert int(indptr.max()) == a.X.nnz
log(f"  CSC 就绪：nnz={a.X.nnz}")

paths = {}
for name, arr in [("data_f32", data), ("indices_i32", indices), ("indptr_i32", indptr)]:
    p = RIO / f"{name}.bin"
    arr.tofile(p)
    paths[name] = {"path": str(p.relative_to(ROOT)), "bytes": int(p.stat().st_size),
                   "sha256": sha256(p)}
    log(f"  写出 {p.name}  {p.stat().st_size/2**30:.2f} GiB")

for name, vals in [("cell_names.txt", a.obs_names), ("gene_names.txt", a.var_names)]:
    p = RIO / name
    p.write_text("\n".join(map(str, vals)) + "\n")
    paths[name] = {"path": str(p.relative_to(ROOT)), "bytes": int(p.stat().st_size),
                   "sha256": sha256(p)}
meta = a.obs[["sample_id", "patient_id", "stage", "stage_token"]].copy()
meta.insert(0, "cell_barcode", a.obs_names)
p = RIO / "cell_meta.csv.gz"
meta.to_csv(p, index=False, compression="gzip")
paths["cell_meta.csv.gz"] = {"path": str(p.relative_to(ROOT)),
                             "bytes": int(p.stat().st_size), "sha256": sha256(p)}
log("  写出 cell_names / gene_names / cell_meta")

man_io = {
    "script": "02_expression/04_rebuild_expression_paperqc.py",
    "purpose": "导出 genes×cells CSC 供 R 侧 Seurat/SCTransform 使用",
    "source_h5ad": {"path": str(NEW), "sha256": new_sha, "n_cells": int(a.n_obs),
                    "n_genes": int(n_var_new), "nnz": int(a.X.nnz)},
    "layout": "genes×cells CSC (dgCMatrix), i 为 0-based 行索引；R 侧用 readBin + new(\"dgCMatrix\")",
    "artifacts": paths,
    "versions": {"anndata": ad.__version__, "numpy": np.__version__, "scipy": scipy.__version__},
    "wall_sec": round(time.time() - T0, 1),
}
(RIO / "export_manifest.json").write_text(json.dumps(man_io, indent=2, ensure_ascii=False))
log(f"  manifest -> {RIO/'export_manifest.json'}")

# ---- 7. 重建 manifest ------------------------------------------------------
stage_counts = a.obs["stage"].value_counts().to_dict()
man = {
    "schema": "luad_v2.expression_build/2",
    "built_utc": time.strftime("%Y-%m-%dT%H:%M:%S+00:00", time.gmtime()),
    "checkpoint": "GP0-redo（论文 QC 口径）",
    "dataset": "GSE308103",
    "supersedes": {"note": "旧的 648,945 核对象保留为敏感性臂，**未删除**",
                   "old_h5ad": str(OLD.relative_to(ROOT)), "old_n_obs": int(n_old)},
    "qc_rule": {"min_nFeature": 500, "min_nCount": 1000, "max_pct_mt": 20,
                "gene_min_cells": GENE_MIN_CELLS,
                "doublet_caller": "scDblFinder (R2)"},
    "n_obs": int(a.n_obs), "n_vars": int(n_var_new), "nnz": int(a.X.nnz),
    "n_genes_dropped_lt3cells": n_drop, "genes_dropped": dropped_genes,
    "stage_counts": {k: int(v) for k, v in stage_counts.items()},
    "n_samples": int(a.obs["sample_id"].nunique()),
    "n_patients": int(a.obs["patient_id"].nunique()),
    "artifacts": {
        "h5ad": {"path": str(NEW.relative_to(ROOT)), "sha256": new_sha,
                 "bytes": int(NEW.stat().st_size)},
        "mask": {"path": str(MASK.relative_to(ROOT)), "sha256": sha256(MASK)},
        "seurat_io_dir": str(RIO.relative_to(ROOT)),
    },
    "inputs": {"old_h5ad": {"path": str(OLD.relative_to(ROOT)), "sha256": sha256(OLD)}},
    "versions": {"anndata": ad.__version__, "numpy": np.__version__,
                 "scipy": scipy.__version__, "pandas": pd.__version__},
    "wall_sec": round(time.time() - T0, 1),
}
MAN.write_text(json.dumps(man, indent=2, ensure_ascii=False))
log(f"manifest -> {MAN}")
log(f"完成：{a.n_obs} 核 × {n_var_new} 基因，nnz={a.X.nnz}")
