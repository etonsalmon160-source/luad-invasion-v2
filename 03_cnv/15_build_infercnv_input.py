#!/usr/bin/env python3
"""切出 P4 的上皮核，做成 inferCNV 的输入。

输入：results/02_expression/gse308103_counts_paperqc.h5ad（原始计数）
      results/05_annotation/epiA_subset_barcodes.txt（133,384 上皮核清单）
      results/03_cnv/infercnv_smoke/gene_order_hg38.tsv（18,000 个有坐标的基因）
输出：results/03_cnv/infercnv_smoke/p4_counts_coo.tsv.gz   （gene_idx, cell_idx, count；1-based）
      p4_genes.txt / p4_cells.txt / p4_annotations.tsv / p4_infercnv_input_manifest.json

为什么要 COO 而不是稠密文本：18,000 x 19,748 的稠密文本约 500 MB，
inferCNV 用 read.table 读会极慢。COO 只有非零元，R 侧用 data.table::fread 秒读后拼 dgCMatrix。

分组口径见 INFERCNV_SMOKE_PREREG.md §四（写死，事后不许改）。
"""
import gzip
import hashlib
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

ROOT = "/home/eto/luad_v2"
SRC = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
BARCODES = f"{ROOT}/results/05_annotation/epiA_subset_barcodes.txt"
CLUSTERS = f"{ROOT}/results/04_integration/seurat_trad/epiA/clusters.csv.gz"
GENE_ORDER = f"{ROOT}/results/03_cnv/infercnv_smoke/gene_order_hg38.tsv"
OUTDIR = f"{ROOT}/results/03_cnv/infercnv_smoke"

PATIENT = "P4"
HOLDOUT_FRAC = 0.15
HOLDOUT_SEED = 0
REF_GROUP = "Normal"
HOLDOUT_GROUP = "Normal_holdout"
EXPECTED_N = 19748
BARCODES_SHA = "d12a115a87d0790d7fe8cf27b16170999d94daba51f3016fdb79c5aa3fdde55e"


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    bsha = sha256(BARCODES)
    assert bsha == BARCODES_SHA, f"上皮清单 sha 不符（{bsha}），停下"

    cl = pd.read_csv(CLUSTERS, usecols=["cell_barcode", "sample_id", "patient_id", "stage"])
    p4 = cl[cl["patient_id"] == PATIENT].copy()
    print(f"[in] {PATIENT} 上皮核 {len(p4)}")
    print(p4.groupby(["stage", "sample_id"]).size().to_string())
    assert len(p4) == EXPECTED_N, f"{PATIENT} 核数 {len(p4)} != 预注册的 {EXPECTED_N}，停下"

    # 只保留有坐标的基因（无坐标的基因 inferCNV 也放不进去；此处**显式**剔除并登记，不靠静默丢）
    go = pd.read_csv(GENE_ORDER, sep="\t", header=None, names=["gene", "chr", "start", "end"])
    keep_genes = set(go["gene"])
    print(f"[in] 有坐标的基因 {len(keep_genes)}")

    print(f"[in] 读 {SRC} ...", flush=True)
    adata = ad.read_h5ad(SRC)
    print(f"[in] {adata.n_obs} 细胞 x {adata.n_vars} 基因")

    cells = p4["cell_barcode"].tolist()
    missing = pd.Index(cells).difference(adata.obs_names)
    assert len(missing) == 0, f"{len(missing)} 个 {PATIENT} 细胞不在 h5ad，例：{list(missing[:3])}"

    sub = adata[cells].copy()
    del adata
    print(f"[out] 子集 {sub.n_obs} 细胞 x {sub.n_vars} 基因")

    var_names = list(sub.var_names)
    assert len(set(var_names)) == len(var_names), "基因名有重复，停下"
    gene_keep_mask = np.array([g in keep_genes for g in var_names])
    n_dropped = int((~gene_keep_mask).sum())
    dropped = [g for g, k in zip(var_names, gene_keep_mask) if not k]
    print(f"[gene] 保留 {int(gene_keep_mask.sum())} / {len(var_names)}；因**无坐标**剔除 {n_dropped}")

    sub = sub[:, gene_keep_mask].copy()
    genes_out = list(sub.var_names)
    X = sub.X.tocsc()
    assert np.allclose(X.data, np.round(X.data)) and X.data.min() >= 0, \
        "X 不是非负整数计数 —— 不是原始计数，停下"
    print(f"[chk] 计数稀疏：nnz={X.nnz}  密度={100.0*X.nnz/(X.shape[0]*X.shape[1]):.1f}%")

    # —— 分组（写死，见预注册 §四）——
    cell_stage = dict(zip(p4["cell_barcode"], p4["stage"]))
    cell_sample = dict(zip(p4["cell_barcode"], p4["sample_id"]))
    normal_cells = sorted([c for c in cells if cell_stage[c] == "Normal"])
    print(f"[grp] 参考候选（Normal）{len(normal_cells)}")

    rng = np.random.default_rng(HOLDOUT_SEED)
    n_hold = int(round(len(normal_cells) * HOLDOUT_FRAC))
    hold_idx = set(rng.choice(len(normal_cells), size=n_hold, replace=False).tolist())
    holdout = {normal_cells[i] for i in hold_idx}
    groups = {}
    for c in cells:
        if cell_stage[c] == "Normal":
            groups[c] = HOLDOUT_GROUP if c in holdout else REF_GROUP
        else:
            groups[c] = cell_sample[c]
    n_ref = sum(1 for v in groups.values() if v == REF_GROUP)
    print(f"[grp] {REF_GROUP}={n_ref}  {HOLDOUT_GROUP}={n_hold}  "
          f"其余观测组={sorted(set(v for v in groups.values() if v not in (REF_GROUP, HOLDOUT_GROUP)))}")

    # —— 写文件 ——
    with open(f"{OUTDIR}/p4_genes.txt", "w") as f:
        f.write("\n".join(genes_out) + "\n")
    with open(f"{OUTDIR}/p4_cells.txt", "w") as f:
        f.write("\n".join(cells) + "\n")

    with open(f"{OUTDIR}/p4_annotations.tsv", "w") as f:
        for c in cells:
            f.write(f"{c}\t{groups[c]}\n")

    # COO（1-based，R 侧方便）；要写成 genes x cells
    # ⚠️ sub.X 是**细胞 x 基因**，故 tocsc 后 coo.row = 细胞序号、coo.col = 基因序号
    coo = X.tocoo()
    coo_path = f"{OUTDIR}/p4_counts_coo.tsv.gz"
    df = pd.DataFrame({
        "gene_idx": (coo.col + 1).astype(np.int32),
        "cell_idx": (coo.row + 1).astype(np.int32),
        "count": coo.data.astype(np.int32),
    })
    df = df.sort_values(["cell_idx", "gene_idx"], kind="stable")
    df.to_csv(coo_path, sep="\t", header=False, index=False, compression="gzip")
    print(f"[out] {coo_path}  {len(df)} 非零元")
    del coo, df

    manifest = {
        "script": "03_cnv/15_build_infercnv_input.py",
        "source_h5ad": SRC,
        "source_h5ad_sha256": "a276cd1af69a4620ff9e3d64f4b93f33461e646537c5f2c9da1aff8c9150b9de",
        "barcodes_file_sha256": bsha,
        "gene_order_file": GENE_ORDER,
        "gene_order_sha256": sha256(GENE_ORDER),
        "patient": PATIENT,
        "n_cells": len(cells),
        "n_genes": len(genes_out),
        "n_genes_dropped_no_position": n_dropped,
        "genes_dropped_no_position": dropped,
        "nnz": int(X.nnz),
        "density_pct": round(100.0 * X.nnz / (X.shape[0] * X.shape[1]), 2),
        "groups": {
            "ref_group": REF_GROUP,
            "holdout_group": HOLDOUT_GROUP,
            "holdout_frac": HOLDOUT_FRAC,
            "holdout_seed": HOLDOUT_SEED,
            "n_ref": n_ref,
            "n_holdout": n_hold,
            "n_per_sample": {k: int(v) for k, v in p4["sample_id"].value_counts().items()},
        },
        "out": {
            "coo": coo_path,
            "coo_sha256": sha256(coo_path),
            "genes": f"{OUTDIR}/p4_genes.txt",
            "cells": f"{OUTDIR}/p4_cells.txt",
            "annotations": f"{OUTDIR}/p4_annotations.tsv",
        },
    }
    with open(f"{OUTDIR}/p4_infercnv_input_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in manifest.items() if k != "genes_dropped_no_position"},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    sys.exit(main())
