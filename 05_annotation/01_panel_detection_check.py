"""上皮亚型 marker 面板的检出率诊断 —— 只回答「面板有没有信号」，不产亚型标签。

为什么需要这一步
----------------
源论文的 Table S2 meta-program 是 NMF 模块，含应激/管家类共表达基因。这些基因
在本数据（10x Flex 探针法，仅 18,082 特征；核的 UMI 低）里可能几乎检不出。
若某亚型的面板基因大面积 0 检出，该亚型在 score_genes 里就只是个常数，
argmax 会系统性地偏向"基因检出好的"亚型 —— 那是技术假象，不是生物学。

本脚本**不聚类、不打标签**，只输出逐基因检出率，供决定面板是否可用。
真正产标签的是 GP8a，须另经签字。

用法:
    python3 05_annotation/01_panel_detection_check.py \
        --h5ad results/02_expression/gse308103_counts_paperqc.h5ad \
        --barcodes results/05_annotation/epiA_subset_barcodes.txt \
        --out results/05_annotation/epi_panel_detection.csv
"""

import argparse
import hashlib
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import epi_subtype_panel as P  # noqa: E402
import epi_classic_panel as C  # noqa: E402


def panels_to_check():
    """返回 [(mode, {亚型: [基因...]})]。现行口径是 classic；paper_mp 保留作历史对照。"""
    out = [("classic", C.build_panel())]
    out.append(("paper_mp", P.build_panel("paper_mp")))
    out.append(("kac_sig", P.build_panel("kac_sig")))
    return out


def _read_index(f, group):
    g = f[group]
    key = "_index" if "_index" in g else list(g.keys())[0]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[key][:]])


def load_subset(h5ad, barcodes_file):
    """返回 (cells x genes 的 CSR 子集, 基因名, 细胞名)。"""
    want = [ln for ln in open(barcodes_file).read().split("\n") if ln]
    if len(set(want)) != len(want):
        raise SystemExit("barcode 清单含重复")
    with h5py.File(h5ad, "r") as f:
        cells = _read_index(f, "obs")
        genes = _read_index(f, "var")
        pos = pd.Index(cells).get_indexer(want)
        miss = int((pos < 0).sum())
        if miss:
            raise SystemExit(f"{miss} 个 barcode 不在 h5ad 里")
        idx = np.sort(pos)

        X = f["X"]
        enc = X.attrs.get("encoding-type")
        if enc == "csr_matrix":
            indptr, indices, data = X["indptr"][:], X["indices"][:], X["data"][:]
        elif enc == "csc_matrix":
            # 按列切片便宜：直接取需要的列
            indptr, indices, data = X["indptr"][:], X["indices"][:], X["data"][:]
            keep = np.zeros(len(cells), dtype=bool)
            keep[idx] = True
            cols = np.repeat(keep, np.diff(indptr))
            new_ptr = np.concatenate([[0], np.cumsum(keep[np.repeat(np.arange(len(cells)), np.diff(indptr))])])
            # csc 的列序保持，行(基因)索引原样
            import scipy.sparse as sp
            n_keep = int(keep.sum())
            mat = sp.csc_matrix((data[cols], indices[cols], new_ptr),
                                shape=(len(genes), n_keep))
            return mat.T.tocsr(), genes, [cells[i] for i in idx]
        else:
            raise SystemExit(f"未预期的 X 编码: {enc}")

    import scipy.sparse as sp
    # csr：按行(细胞)切片
    start, end = indptr[idx], indptr[idx + 1]
    counts = end - start
    new_ptr = np.concatenate([[0], np.cumsum(counts)])
    sel = np.concatenate([np.arange(s, e) for s, e in zip(start, end)]) if len(idx) else np.array([], int)
    mat = sp.csr_matrix((data[sel], indices[sel], new_ptr), shape=(len(idx), len(genes)))
    return mat, genes, [cells[i] for i in idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="results/02_expression/gse308103_counts_paperqc.h5ad")
    ap.add_argument("--barcodes", default="results/05_annotation/epiA_subset_barcodes.txt")
    ap.add_argument("--out", default="results/05_annotation/epi_panel_detection.csv")
    a = ap.parse_args()

    print(f"读 {a.h5ad} …")
    mat, genes, cells = load_subset(a.h5ad, a.barcodes)
    print(f"子集: {mat.shape[0]:,} 细胞 × {mat.shape[1]:,} 基因, nnz={mat.nnz:,}")

    gidx = {g: i for i, g in enumerate(genes)}
    n_cells = mat.shape[0]
    det = np.asarray((mat > 0).sum(axis=0)).ravel() / n_cells  # 逐基因检出细胞比例

    rows = []
    for mode, pan in panels_to_check():
        for sub, glist in pan.items():
            present = [g for g in glist if g in gidx]
            absent = [g for g in glist if g not in gidx]
            v = np.array([det[gidx[g]] for g in present]) if present else np.array([])
            rows.append(dict(
                mode=mode, subtype=sub, n_genes=len(glist),
                n_in_matrix=len(present), n_absent=len(absent),
                absent=";".join(absent),
                det_median=float(np.median(v)) if len(v) else np.nan,
                det_p10=float(np.percentile(v, 10)) if len(v) else np.nan,
                det_max=float(v.max()) if len(v) else np.nan,
                n_det_gt5pct=int((v > 0.05).sum()),
                n_det_lt1pct=int((v < 0.01).sum()),
            ))
            # 逐基因明细
            for g in present:
                rows[-1].setdefault("_per_gene", []).append((g, float(det[gidx[g]])))

    df = pd.DataFrame([{k: v for k, v in r.items() if k != "_per_gene"} for r in rows])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    df.to_csv(a.out, index=False)
    print("\n" + df.to_string(index=False))

    # 逐基因长表
    long = []
    for r in rows:
        for g, d in r["_per_gene"]:
            long.append(dict(mode=r["mode"], subtype=r["subtype"], gene=g, det_frac=d))
    ldf = pd.DataFrame(long)
    lout = a.out.replace(".csv", "_per_gene.csv")
    ldf.to_csv(lout, index=False)

    # 最差的基因（可能在拖后腿）
    worst = ldf.nsmallest(25, "det_frac")
    print(f"\n检出率最低的 25 个基因（<1% 的基因对 score_genes 几乎无贡献）:\n"
          f"{worst.to_string(index=False)}")

    print(f"\n写出 {a.out}")
    print(f"写出 {lout}")
    print("panel sha256:", hashlib.sha256(open(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "epi_subtype_panel.py"), "rb").read()).hexdigest())
    print("barcodes sha256:", hashlib.sha256(open(a.barcodes, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
