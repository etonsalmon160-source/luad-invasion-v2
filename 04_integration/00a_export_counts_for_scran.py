#!/usr/bin/env python3
# ============================================================================
# 🔴【已作废 · 2026-09-16】本脚本**请勿运行**，保留仅为审计留痕。
#   作废原因：用户决策「对齐论文全套」⇒ M3-A.1 的归一化主口径改为
#     SCTransform(vst.flavor="v2") → HVG 3000 → PCA 50 → Harmony(默认参数)。
#     论文配方里**没有 scran 这一步**，故本脚本登记的 scran 口径整体作废。
#   登记位置：docs/PARAMETERS_AND_SOURCES.md §M3-A.1「已作废登记」表。
#   现行替代：04_integration/10_seurat_traditional.R
#   ⚠️ 「作废」不等于「结论无效」：由本链产出的 `00c_stage_fragility_check.py`
#      的 AAH 脆弱性结论**仍然成立**（已在论文 QC 口径下复核，见
#      results/01_qc/stage_fragility_report_paperqc.json）。
# ============================================================================
# 00a_export_counts_for_scran.py —— 把原始计数导出为 R 可直接读的二进制
#
# 为什么需要这一步：
#   `docs/PARAMETERS_AND_SOURCES.md` M3-A.1 登记的归一化**主口径**是
#   scran 池化 size factor（Lun 2016, Genome Biol 17:75），scran 只有 R 实现。
#   本脚本把 GP0 的原始计数以**零拷贝语义**导出为 genes×cells 的 CSC 三元组，
#   R 侧用 readBin + new("dgCMatrix") 直接构造，无需 MatrixMarket 文本往返
#   （文本 .mtx 约 16 GB / 读写各十几分钟；二进制 3.2 GB / 分钟级）。
#
# 关键性质：
#   · h5ad 的 X 是 cells×genes 的 CSR；scipy 的 `X.T` 对 CSR 返回的是
#     **同一块内存**被解释为 CSC 的 genes×cells 矩阵（不复制），正是 dgCMatrix 的布局。
#   · dgCMatrix 的 `i` 槽是 **0-based** 行索引 —— 与 scipy 一致，不需 ±1。
#   · nnz=791,571,728 < 2^31，故 indptr 可安全存 int32。
#
# ⚠️ 本脚本**不修改任何数据**，只做只读导出 + 新写文件。
# ============================================================================
import json, time, hashlib, sys
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse as sp
import anndata as ad

ROOT = Path("/home/eto/luad_v2")
INP  = ROOT / "results/02_expression/gse308103_counts.h5ad"
OUT  = ROOT / "results/04_integration/scran_io"
OUT.mkdir(parents=True, exist_ok=True)

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

def sha256(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while (b := f.read(buf)):
            h.update(b)
    return h.hexdigest()

log(f"载入 {INP}")
a = ad.read_h5ad(INP)
n_cells, n_genes = a.shape
log(f"  cells={n_cells}  genes={n_genes}  nnz={a.X.nnz}")

# ---- 前置断言：必须是原始整数计数（scran 的池化求和前提）-------------------
d = a.X.data
if not np.all(d[: 5_000_000] == np.floor(d[: 5_000_000])):
    sys.exit("[FAIL] 输入非整数计数，拒绝导出")
if a.layers or a.raw is not None or a.obsm:
    sys.exit(f"[FAIL] 输入对象带归一化痕迹 layers={list(a.layers)} "
             f"raw={a.raw is not None} obsm={list(a.obsm)}，拒绝导出")
log("前置断言通过：原始整数计数、无归一化痕迹")

# ---- 1. 布局转换（cells×genes CSR → genes×cells CSC，共享内存）-------------
A = a.X.T.tocsc()
assert A.shape == (n_genes, n_cells)
data, indices, indptr = A.data, A.indices, A.indptr
assert data.dtype == np.float32 and indices.dtype == np.int32, (data.dtype, indices.dtype)
assert indptr.max() == a.X.nnz, (int(indptr.max()), a.X.nnz)
log(f"  CSC(gene×cell) 就绪：data={data.nbytes/2**30:.2f} GiB "
    f"indices={indices.nbytes/2**30:.2f} GiB indptr={indptr.nbytes/2**30:.2f} GiB")

# ---- 2. 落盘 --------------------------------------------------------------
paths = {}
for name, arr in [("data_f32", data), ("indices_i32", indices), ("indptr_i32", indptr)]:
    p = OUT / f"{name}.bin"
    arr.tofile(p)
    paths[name] = {"path": str(p), "bytes": int(p.stat().st_size), "sha256": sha256(p)}
    log(f"  写出 {p.name}  {p.stat().st_size/2**30:.2f} GiB")

for name, vals in [("cell_names.txt", a.obs_names), ("gene_names.txt", a.var_names)]:
    p = OUT / name
    p.write_text("\n".join(map(str, vals)) + "\n")
    paths[name] = {"path": str(p), "bytes": int(p.stat().st_size), "sha256": sha256(p)}
log("  写出 cell_names.txt / gene_names.txt")

# ---- 3. 细胞元数据（矩阵列顺序严格一致；R 侧按 sample_id 分块）------------
meta = a.obs[["sample_id", "patient_id", "stage"]].copy()
meta.insert(0, "cell_barcode", a.obs_names)
p = OUT / "cell_meta.csv.gz"
meta.to_csv(p, index=False, compression="gzip")
paths["cell_meta.csv.gz"] = {"path": str(p), "bytes": int(p.stat().st_size), "sha256": sha256(p)}
log(f"  写出 cell_meta.csv.gz（{len(meta)} 行）")

# ---- 4. manifest ---------------------------------------------------------
man = {
    "script": "04_integration/00a_export_counts_for_scran.py",
    "purpose": "导出 genes×cells CSC 供 scran::computeSumFactors 使用",
    "source": {"path": str(INP), "sha256": sha256(INP),
               "n_cells": int(n_cells), "n_genes": int(n_genes), "nnz": int(a.X.nnz)},
    "layout": "genes×cells CSC (dgCMatrix), i 为 0-based 行索引",
    "artifacts": paths,
    "versions": {"anndata": ad.__version__, "numpy": np.__version__, "scipy": scipy.__version__},
    "wall_sec": round(time.time() - T0, 1),
}
(OUT / "export_manifest.json").write_text(json.dumps(man, indent=2, ensure_ascii=False))
log(f"完成 → {OUT}")
