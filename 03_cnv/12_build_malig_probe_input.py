#!/usr/bin/env python3
"""切出上皮子集，供 scMalignantFinder 诊断用。

输入：results/02_expression/gse308103_counts_paperqc.h5ad（原始计数）
      results/05_annotation/epiA_subset_barcodes.txt（133,384 上皮细胞，未再剔除）
输出：results/03_cnv/malig_probe/malig_probe_input.h5ad（原始计数，仅上皮）

只做切片，不改数值、不做归一化（归一化交给 scMalignantFinder 的 norm_type 处理）。
"""
import hashlib
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = "/home/eto/luad_v2"
SRC = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
BARCODES = f"{ROOT}/results/05_annotation/epiA_subset_barcodes.txt"
OUTDIR = f"{ROOT}/results/03_cnv/malig_probe"
OUT = f"{OUTDIR}/malig_probe_input.h5ad"


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    barcodes_sha = sha256(BARCODES)
    with open(BARCODES) as f:
        cells = [ln.strip() for ln in f if ln.strip()]
    print(f"[in] 上皮清单 {len(cells)} 细胞  file_sha256={barcodes_sha}")
    assert barcodes_sha == "d12a115a87d0790d7fe8cf27b16170999d94daba51f3016fdb79c5aa3fdde55e", \
        "上皮清单与 run_manifest 登记不符，停下"
    assert len(cells) == len(set(cells)) == 133384, "上皮清单有重复行或行数不符，停下"

    print(f"[in] 读 {SRC} ...", flush=True)
    adata = ad.read_h5ad(SRC)
    print(f"[in] {adata.n_obs} 细胞 x {adata.n_vars} 基因")
    print(f"[in] obs 列: {list(adata.obs.columns)}")
    print(f"[in] obs_names[:2]: {list(adata.obs_names[:2])}")

    idx = pd.Index(cells)
    missing = idx.difference(adata.obs_names)
    assert len(missing) == 0, f"{len(missing)} 个清单细胞不在 h5ad，例：{list(missing[:3])}"

    sub = adata[idx].copy()
    print(f"[out] 子集 {sub.n_obs} 细胞 x {sub.n_vars} 基因")

    X = sub.X
    assert sp.issparse(X), "X 不是稀疏矩阵，停下确认"
    Xc = X.tocsc()
    mins = np.asarray(Xc.min(axis=0).todense()).ravel()
    is_int = np.allclose(X.data, np.round(X.data))
    print(f"[chk] X dtype={X.dtype} 最小={X.data.min()} 最大={X.data.max()} 全整数值={is_int}")
    assert is_int and X.data.min() >= 0, "X 不是非负整数计数 —— 不是原始计数，停下"
    assert sub.n_vars == 18069, f"基因数 {sub.n_vars} != 18069，停下"

    sub.write_h5ad(OUT, compression="gzip")
    print(f"[out] 写出 {OUT}  {os.path.getsize(OUT)/1e6:.1f} MB")

    manifest = {
        "script": "03_cnv/12_build_malig_probe_input.py",
        "source_h5ad": SRC,
        "source_h5ad_sha256": "a276cd1af69a4620ff9e3d64f4b93f33461e646537c5f2c9da1aff8c9150b9de",
        "barcodes_file": BARCODES,
        "barcodes_file_sha256": barcodes_sha,
        "n_cells_in": int(adata.n_obs),
        "n_cells_out": int(sub.n_obs),
        "n_genes": int(sub.n_vars),
        "is_raw_counts": True,
        "out": OUT,
        "out_sha256": sha256(OUT),
        "out_bytes": os.path.getsize(OUT),
    }
    with open(f"{OUTDIR}/malig_probe_input_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    sys.exit(main())
