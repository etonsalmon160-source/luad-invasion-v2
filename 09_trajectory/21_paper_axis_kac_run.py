#!/usr/bin/env python
"""
预注册 20_paper_axis_kac_prereg.md 的实现：用原文 MP 面板给现有 27 簇加"论文口径"注释。

判据（签字 T1）：
  D1 主判 = 8 型 MP 的簇级 argmax ∈ {Tumor cell/KAC, KAC/inflammatory}
  D2 佐证 = 27 簇上 MP6/KAC-inflam 与 Table S3 KAC 签名的簇级 Spearman（不当门）

对照 K1（深度）/ K2（同长度同检出率分位的置换零分布，200 次）/ K3（单一患者主导）。

边界（签字 T0/T2 与 §6 禁令）：
  · 不重聚类；打分对象 = 现有 27 簇
  · 不改任何 results/05_annotation/epiA_* 既有文件
  · "Tumor cell" 是论文给的名字，不是恶性判定（单细胞 CNV 已退场，无逐细胞恶性标签）
  · 打分参数与 GP8a 逐字相同（同一把尺子）

用法:
    python3 09_trajectory/21_paper_axis_kac_run.py [--smoke N]
"""
import argparse
import gc
import hashlib
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

import anndata as ad
import scanpy as sc

BASE = "/home/eto/luad_v2/results"
WT = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963"
H5AD = f"{BASE}/02_expression/gse308103_counts_paperqc.h5ad"
BARCODES = f"{BASE}/05_annotation/epiA_subset_barcodes.txt"
CLUSTERS = f"{BASE}/04_integration/seurat_trad/epiA/clusters.csv.gz"
OUT = f"{BASE}/09_trajectory/paper_axis"

CLUSTER_COL = "harmony_res0.7_seed0"     # 与 epiA_annotation_manifest.json 一致
CTRL_SIZE = 50                            # 与 GP8a 注册值逐字相同
N_BINS = 25
SCORE_SEED = 0
N_REP = 200                               # K2 置换次数（预注册 §5）
RANDOM_SEED = 20260927

KAC_SIDE = {"Tumor cell/KAC", "KAC/inflammatory"}   # D1 命中集
NEWLINE = "\n"

_T0 = time.time()
_HWM = 0.0


def log(msg):
    global _HWM
    try:
        for ln in open("/proc/self/status"):
            if ln.startswith("VmHWM:"):
                _HWM = max(_HWM, int(ln.split()[1]) / 1024 / 1024)
    except Exception:
        pass
    print(f"[{time.time()-_T0:7.1f}s | {_HWM:5.2f}GB] {msg}", flush=True)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def colname(subtype):
    """亚型名 -> 合法列名（避免空格与斜杠）。"""
    return "s_" + "".join(c if c.isalnum() else "_" for c in subtype)


def _read_index(f, group):
    g = f[group]
    key = "_index" if "_index" in g else list(g.keys())[0]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[key][:]])


def _read_obs_col(f, name, row_idx):
    """读 obs 的一列（自动处理 categorical 与 bytes）。"""
    node = f["obs"][name]
    if isinstance(node, h5py.Group):          # categorical
        cats = [x.decode() if isinstance(x, bytes) else str(x)
                for x in node["categories"][:]]
        codes = node["codes"][:]
        sel = codes[row_idx]
        return np.array([cats[c] if c >= 0 else "" for c in sel], dtype=object)
    v = node[:]
    if v.dtype.kind == "S":
        v = np.array([x.decode() for x in v], dtype=object)
    return v[row_idx]


def load_subset(h5ad, barcodes_file, smoke=None):
    """返回 (CSR 子集, 基因名, 细胞名, 行号, obs 列字典)。"""
    want = [ln for ln in open(barcodes_file).read().split(NEWLINE) if ln]
    if len(set(want)) != len(want):
        raise SystemExit("barcode 清单含重复")
    log(f"barcode 清单 {len(want):,} 行")

    with h5py.File(h5ad, "r") as f:
        cells = _read_index(f, "obs")
        genes = _read_index(f, "var")
        pos = pd.Index(cells).get_indexer(want)
        miss = int((pos < 0).sum())
        if miss:
            raise SystemExit(f"🔴 {miss} 个 barcode 不在 h5ad 里")
        idx = np.sort(pos)
        if smoke:
            idx = idx[:smoke]
            log(f"🔧 SMOKE 模式：只用前 {smoke} 个细胞")

        obs = {c: _read_obs_col(f, c, idx)
               for c in ["patient_id", "stage", "nFeature"] if c in f["obs"]}

        X = f["X"]
        enc = X.attrs.get("encoding-type")
        if enc != "csr_matrix":
            raise SystemExit(f"🔴 未预期的 X 编码: {enc}（本脚本只处理 csr）")
        indptr, indices, data = X["indptr"][:], X["indices"][:], X["data"][:]
        names = [cells[i] for i in idx]

    start, end = indptr[idx], indptr[idx + 1]
    counts = end - start
    new_ptr = np.concatenate([[0], np.cumsum(counts)])
    sel = np.concatenate([np.arange(s, e) for s, e in zip(start, end)]) if len(idx) else np.array([], int)
    mat = sp.csr_matrix((data[sel], indices[sel], new_ptr), shape=(len(idx), len(genes)))
    log(f"子集读入 {mat.shape[0]:,} × {mat.shape[1]:,}, nnz={mat.nnz:,}")
    return mat, genes, names, idx, obs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", type=int, default=0, help="只用前 N 个细胞（调试）")
    ap.add_argument("--nrep", type=int, default=N_REP, help="K2 置换次数（默认=预注册的 200）")
    a = ap.parse_args()
    n_rep = min(a.nrep, 5) if a.smoke else a.nrep

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(f"{OUT}/figures", exist_ok=True)

    sys.path.insert(0, f"{WT}/05_annotation")
    import epi_subtype_panel as P            # noqa: E402

    panel_mp = P.build_panel("paper_mp")     # 8 型
    panel_kac = P.build_panel("kac_sig")     # 1 型（KAC 签名的独立第二来源）
    log(f"面板 paper_mp {len(panel_mp)} 型 / kac_sig {len(panel_kac['KAC'])} 基因")

    mat, genes, cells, idx, obs = load_subset(H5AD, BARCODES, a.smoke or None)
    gidx = {g: i for i, g in enumerate(genes)}

    # ---- 检出率（K2 置换要用；在归一化之前算）----
    n_cells = mat.shape[0]
    det = np.asarray((mat > 0).sum(axis=0)).ravel() / n_cells

    # ---- 组装 AnnData ----
    adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=cells),
                       var=pd.DataFrame(index=genes))
    for k, v in obs.items():
        adata.obs[k] = v
    if "nFeature" in adata.obs:
        adata.obs["nFeature"] = pd.to_numeric(adata.obs["nFeature"], errors="coerce")

    # ---- 簇标签（与 GP8a 同一份）----
    clu = pd.read_csv(CLUSTERS, usecols=["cell_barcode", CLUSTER_COL])
    lab = clu.set_index("cell_barcode")[CLUSTER_COL].reindex(adata.obs_names)
    n_nan = int(lab.isna().sum())
    if n_nan:
        raise SystemExit(f"🔴 {n_nan} 个细胞没匹配到簇标签")
    adata.obs["cluster"] = pd.Categorical(lab.astype(str).values)
    ucl = sorted(adata.obs["cluster"].astype(str).unique(), key=lambda x: (len(x), x))
    log(f"簇 {len(ucl)} 个：{ucl}")

    # ---- 归一化 + 打分（参数与 GP8a 逐字相同）----
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    cols = {}
    for subtype, glist in panel_mp.items():
        present = [g for g in glist if g in gidx]
        if len(present) < 2:
            log(f"🔴 {subtype} 可用基因 <2，跳过")
            continue
        c = colname(subtype)
        sc.tl.score_genes(adata, present, ctrl_size=CTRL_SIZE, n_bins=N_BINS,
                          random_state=SCORE_SEED, score_name=c, use_raw=False)
        cols[c] = subtype
        log(f"打分 {subtype:<40s} ({len(present)} 基因)")

    kac_genes = [g for g in panel_kac["KAC"] if g in gidx]
    ck = colname("KACsig")
    sc.tl.score_genes(adata, kac_genes, ctrl_size=CTRL_SIZE, n_bins=N_BINS,
                      random_state=SCORE_SEED, score_name=ck, use_raw=False)
    cols[ck] = "KACsig"
    log(f"打分 {'KACsig (Table S3)':<40s} ({len(kac_genes)} 基因，"
        f"矩阵外 {len(panel_kac['KAC']) - len(kac_genes)} 个)")

    # ---- 逐簇均值 ----
    cvals = adata.obs["cluster"].astype(str).values
    allcols = list(cols.keys())
    M = adata.obs[allcols].to_numpy(dtype=np.float64)
    nfeat = adata.obs["nFeature"].to_numpy()
    pat = adata.obs["patient_id"].astype(str).to_numpy()
    rows = []
    for c in ucl:
        m = cvals == c
        if m.sum() == 0:
            continue
        vc = pd.Series(pat[m]).value_counts()
        d = dict(cluster=c, n_cells=int(m.sum()),
                 nFeature_median=float(np.nanmedian(nfeat[m])),
                 n_patients=int(len(vc)),
                 top_patient_frac=float(vc.iloc[0] / m.sum()))
        for c2, s in zip(allcols, M[m].mean(axis=0)):
            d[c2] = float(s)
        rows.append(d)
    cm = pd.DataFrame(rows)
    ucl = cm["cluster"].tolist()          # 只留有细胞的簇
    log(f"参与打分的簇 {len(ucl)} 个")
    cm.to_csv(f"{OUT}/cluster_scores_papermp.csv", index=False)
    log(f"写出 cluster_scores_papermp.csv（{len(ucl)} 簇 × {len(allcols)} 列）")

    # ---- D1：8 型 MP 的 argmax ----
    mp_cols = [c for c in allcols if cols[c] != "KACsig"]
    mp_types = [cols[c] for c in mp_cols]
    A = cm[mp_cols].to_numpy()
    am = np.array(mp_types)[A.argmax(axis=1)]
    srt = np.sort(A, axis=1)
    margin = (srt[:, -1] - srt[:, -2]) if A.shape[1] > 1 else np.zeros(len(ucl))
    cm.insert(2, "argmax_papermp", am)
    cm.insert(3, "win_margin", np.round(margin, 4))
    cm["is_kac_side"] = [t in KAC_SIDE for t in am]
    cm.to_csv(f"{OUT}/cluster_scores_papermp.csv", index=False)   # 带判定的版本覆盖

    hit = cm.loc[cm["is_kac_side"], ["cluster", "n_cells", "argmax_papermp", "win_margin",
                                     "nFeature_median", "top_patient_frac"]]
    log("D1 命中（KAC 侧）簇：" + NEWLINE + (hit.to_string(index=False) if len(hit) else "  （无）"))

    # ---- D2：簇级 Spearman（KAC 侧两 MP vs KAC 签名）----
    d2 = {}
    for t in ["Tumor cell/KAC", "KAC/inflammatory"]:
        c = colname(t)
        if c in allcols:
            rho, p = stats.spearmanr(cm[c].to_numpy(), cm[ck].to_numpy())
            d2[t] = {"spearman_rho": float(rho), "p": float(p)}
    log(f"D2 簇级 Spearman vs KAC 签名: {json.dumps(d2, ensure_ascii=False)}")

    # ---- K1 深度对照 ----
    k1 = {"spearman_nFeature_vs_KACsig": float(stats.spearmanr(
              cm["nFeature_median"].to_numpy(), cm[ck].to_numpy())[0]),
          "nFeature_median_kacside": (float(cm.loc[cm["is_kac_side"], "nFeature_median"].median())
                                      if cm["is_kac_side"].any() else None),
          "nFeature_median_rest": float(cm.loc[~cm["is_kac_side"], "nFeature_median"].median())}
    log(f"K1 {json.dumps(k1, ensure_ascii=False)}")

    # ---- K2 置换零分布：同长度、同检出率分位 ----
    # 统计量 = 27 个簇均值的 SD（真签名若切分簇结构，SD 应显著大于零分布）
    kac_arr = adata.obs[ck].to_numpy(dtype=np.float64)
    real_sd = float(np.std([kac_arr[cvals == c].mean() for c in ucl]))

    nb = 20
    edges = np.quantile(det, np.linspace(0, 1, nb + 1))
    edges[-1] += 1e-9
    binid = np.clip(np.digitize(det, edges) - 1, 0, nb - 1)
    used = set()
    for gl in list(panel_mp.values()) + [panel_kac["KAC"]]:
        used.update(g for g in gl if g in gidx)
    pool_by_bin = {b: np.array([i for i in range(len(genes))
                                if binid[i] == b and genes[i] not in used])
                   for b in range(nb)}
    rng = np.random.default_rng(RANDOM_SEED)
    kac_bins = binid[[gidx[g] for g in kac_genes]]
    null = []
    for _ in range(n_rep):
        pick = [int(rng.choice(pool_by_bin[b])) if len(pool_by_bin[b]) else int(rng.integers(len(genes)))
                for b in kac_bins]
        pick = [i for i in pick if genes[i] not in used]
        if len(pick) < 2:
            continue
        sc.tl.score_genes(adata, [genes[i] for i in pick], ctrl_size=CTRL_SIZE, n_bins=N_BINS,
                          random_state=SCORE_SEED, score_name="_null")
        v = adata.obs["_null"].to_numpy(dtype=np.float64)
        null.append(float(np.std([v[cvals == c].mean() for c in ucl])))
        del adata.obs["_null"]
    null = np.array(null)
    k2 = {"n_rep": int(len(null)), "stat": "SD of 27 cluster means",
          "real_sd": real_sd, "null_mean": float(null.mean()),
          "null_p95": float(np.percentile(null, 95)) if len(null) else None,
          "null_max": float(null.max()) if len(null) else None,
          "real_beats_p95": bool(len(null) and real_sd > np.percentile(null, 95))}
    log(f"K2 {json.dumps(k2, ensure_ascii=False)}")

    # ---- 逐细胞分数落盘 ----
    pd.DataFrame({"cell": adata.obs_names, "cluster": cvals,
                  "stage": adata.obs["stage"].astype(str).values,
                  "patient_id": adata.obs["patient_id"].astype(str).values,
                  "nFeature": adata.obs["nFeature"].to_numpy(),
                  **{c: adata.obs[c].to_numpy(dtype=np.float64) for c in allcols}}) \
        .to_csv(f"{OUT}/per_cell_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    log("写出 per_cell_scores.tsv.gz")

    # ---- K3 弱证据标记 ----
    k3 = cm.loc[cm["is_kac_side"] & (cm["top_patient_frac"] > 0.5), "cluster"].tolist()

    rep = {
        "prereg": "09_trajectory/20_paper_axis_kac_prereg.md",
        "signed": "T0=解冻(限范围) / T1=D1主判+D2佐证 / T2=全8MP都打",
        "input_h5ad": H5AD, "n_cells": int(len(cells)),
        "cluster_col": CLUSTER_COL, "n_clusters": len(ucl),
        "score_genes": {"ctrl_size": CTRL_SIZE, "n_bins": N_BINS, "random_state": SCORE_SEED,
                        "gene_pool": None, "use_raw": False},
        "panel_classes": cols,
        "kac_sig_n_in_matrix": len(kac_genes),
        "kac_sig_n_absent": int(len(panel_kac["KAC"]) - len(kac_genes)),
        "D1": {"rule": "argmax(8 MP) ∈ {Tumor cell/KAC, KAC/inflammatory}",
               "n_hit": int(cm["is_kac_side"].sum()),
               "hits": hit.to_dict(orient="records")},
        "D2": d2,
        "K1_depth": k1, "K2_permutation": k2,
        "K3_single_patient_dominated": k3,
        "sha256": {"panel": sha256(f"{WT}/05_annotation/epi_subtype_panel.py"),
                   "barcodes": sha256(BARCODES)},
    }
    with open(f"{OUT}/kac_axis_summary.json", "w") as fh:
        json.dump(rep, fh, indent=2, ensure_ascii=False)

    log(f"===== 完成；峰值 {_HWM:.2f} GB =====")
    print(json.dumps(rep, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
