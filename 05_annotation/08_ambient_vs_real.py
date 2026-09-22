"""上皮亚簇里的外来基因：**是真表达，还是环境 RNA 的污染**？

为什么要有这一步
----------------
2026-09-21 GP8a 实跑后，好几个上皮簇的 top 富集基因是**别的谱系**的经典 marker：
  簇 11 = DCN/COL1A2/COL3A1/LUM/SPARC（成纤维）
  簇 21 = IGKC/IGHG1/IGHA1/JCHAIN（浆细胞）
  簇 15 = IGKC/IGHG1 + EGFR/TIMP1
而同一个细胞按**六谱系面板**打分，最强的仍是"上皮"。

两种解释都成立，而且**后果完全相反**：
  ① 环境 RNA —— 游离 mRNA 吸附到核上（snRNA 通病）。**是假的**，该扣掉。
  ② 真生物学 —— 上皮细胞真在表达这些（EMT 就表达 COL1A1/FN1/VIM）。**是真的**，是发现。
分不出来的话，RCTD 解卷积的参考谱就带偏差（见 2026-09-21 的分期相关讨论）。

判据：**用阳性对照自校准，不设阈值**
------------------------------------
对每个簇，把待查基因和**该簇自己的真 marker**（阳性对照）放在**同一把尺子**上量：
  - `pct_in`  簇内表达该基因的细胞比例（原始计数 >0 即算表达）
  - `pct_out` 簇外表达比例
  - `mean_in` 簇内平均表达（按 nCount 归一化到 1e4 再 log1p）
  - `log2r`   簇内均值 / 簇外均值的 log2（偏在本簇的程度）
  `pct_out` 兼作**阴性对照** —— 环境 RNA 的特点是"到处都有、簇内只比簇外高一点点"。

读法（**不设阈值，看形状**）：
  - 外来基因长得像阳性对照（pct_in 高、且远高于 pct_out、log2r 大）⇒ **像真表达**
  - 外来基因 pct_in 只比 pct_out 高一点、log2r 接近 0 ⇒ **像环境 RNA**
  阳性对照给出"真"的标尺，阴性对照给出"背景"的标尺，外来基因落在哪边就是哪边。

实现
----
**不加载整个 h5ad**（413,697 × 18,069 CSR，655M 非零，约需 20+ GB）。只按列取 61 个
目标基因：读 `indices` 一次，用 `isin` 选出命中的非零元，再由 `indptr` 反推所属细胞。
峰值内存 ≈ 6 GB，且不挤占正在跑的其他任务。

⚠️ 本脚本**只诊断**：不改数据、不删细胞、不决定去留。

用法:
    python3 05_annotation/08_ambient_vs_real.py --tag epiA --lineage 上皮 --res 0.5 --seed 0
"""

import argparse
import gzip
import hashlib
import json
import os
import time

import h5py
import numpy as np
import pandas as pd

ROOT = "/home/eto/luad_v2"
H5 = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ROOT}/results/05_annotation"
SUBSETS = f"{OUT}/lineage_subsets_manifest.json"
LABELS = f"{OUT}/gp6_cell_labels.csv.gz"
TARGET_SUM = 1e4

# ---- 待查基因：只放**真正在 top 富集里出现过**的那些，不凭空加 ----
SUSPECT = {
    "成纤维": ["DCN", "COL1A1", "COL1A2", "COL3A1", "LUM", "SPARC", "FN1", "VIM", "COL6A3"],
    "B/浆":   ["IGKC", "IGHG1", "IGHA1", "JCHAIN", "IGHG3", "IGHM", "IGLC2"],
    "髓系":   ["FTL", "S100A8", "S100A9", "LYZ", "CTSB", "TYROBP", "AIF1"],
    "T/NK":   ["TRAC", "CD2", "ITGAL", "TRBC2", "CORO1A", "IKZF1", "PTPRC"],
    "内皮":   ["PECAM1", "VWF", "CLDN5"],
}

# ---- 阳性对照：该亚型**自己的**经典 marker，给"真表达"当标尺 ----
POS_CTRL = {
    "AT1": ["AGER", "CAV1", "EMP2"],
    "AT2": ["SFTPC", "SFTPA1", "ABCA3", "NAPSA"],
    "Ciliated": ["FOXJ1", "DNAH12"],
    "Basal": ["KRT15", "KRT5", "TP63"],
    "Goblet/Mucous": ["MUC5B", "TFF3"],
    "Club": ["SCGB1A1", "SCGB3A1"],
    "Serous": ["LTF", "DMBT1"],
    "Tuft": ["POU2F3"],
    "Ionocyte": ["FOXI1"],
    "Neuroendocrine": ["ASCL1"],
}

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def sha256(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def read_str_col(node, key="_index"):
    key = key if key in node else list(node.keys())[0]
    v = node[key][:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--res", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    a = ap.parse_args()

    man = json.load(open(SUBSETS, encoding="utf-8"))
    bc_file = f"{OUT}/{man['lineages'][a.lineage]['file']}"
    sub_bc = [ln for ln in open(bc_file).read().split("\n") if ln]
    log(f"{a.lineage} 子集 {len(sub_bc):,} 核（清单 {os.path.basename(bc_file)}）")

    clu_p = f"{TRAD}/{a.tag}/clusters.csv.gz"
    clu = pd.read_csv(clu_p)
    col = f"harmony_res{a.res}_seed{a.seed}"
    if col not in clu.columns:
        raise SystemExit(f"🔴 集群表没有 {col}；有 {[c for c in clu.columns if 'res' in c][:6]}")

    # 亚型标注（GP8a 的注释表，拿每簇 argmax）
    ann = pd.read_csv(f"{OUT}/{a.tag}_cluster_annotation.csv")
    amap = dict(zip(ann["cluster"].astype(str), ann["argmax"]))

    f = h5py.File(H5, "r")
    gene_names = read_str_col(f["var"])
    g2j = {g: i for i, g in enumerate(gene_names)}
    obs_bc = read_str_col(f["obs"])
    ncount = f["obs"]["nCount"][:].astype(np.float64)
    log(f"矩阵 {len(obs_bc):,} 细胞 × {len(gene_names):,} 基因")

    # ---- 目标基因 → 列号 ----
    targets, meta = [], []
    for lin, gs in SUSPECT.items():
        for g in gs:
            targets.append(g); meta.append((f"待查:{lin}", lin, g))
    pos_by_lin = {}
    for lin, gs in SUSPECT.items():
        pos_by_lin[lin] = gs
    for gs in POS_CTRL.values():
        for g in gs:
            targets.append(g); meta.append(("阳性对照", "自身", g))

    # 去掉矩阵里没有的基因，并如实登记
    keep_t, keep_m, missing = [], [], []
    for t, m in zip(targets, meta):
        if t in g2j:
            keep_t.append(t); keep_m.append(m)
        else:
            missing.append(t)
    tjs = np.array(sorted({g2j[t] for t in keep_t}), dtype=np.int32)
    log(f"目标基因 {len(keep_t)} 个（去重后 {len(tjs)} 列）"
        + (f"；矩阵里没有、已跳过：{sorted(set(missing))}" if missing else ""))

    # ---- 一次性扫 CSR，抽出命中的非零元 ----
    indptr = f["X"]["indptr"][:]
    log("读 indices / data …")
    indices = f["X"]["indices"][:]
    data = f["X"]["data"][:]
    f.close()
    log(f"非零元 {len(indices):,}；筛选目标基因 …")
    sel = np.isin(indices, tjs)
    hit_pos = np.flatnonzero(sel)
    hit_gene = indices[hit_pos].astype(np.int64)
    hit_val = data[hit_pos].astype(np.float64)
    # 非零元位置 → 所属细胞
    hit_cell = np.searchsorted(indptr, hit_pos, side="right") - 1
    del sel, hit_pos, indices, data
    log(f"命中 {len(hit_gene):,} 个非零元")

    # 反查 gene 名（列号被排序过）
    inv = {j: g for g, j in g2j.items()}

    # ---- 细胞 → 簇 ----
    order = pd.Index(obs_bc)
    ci = order.get_indexer(clu["cell_barcode"].to_numpy())
    if (ci < 0).any():
        raise SystemExit(f"🔴 有 {int((ci < 0).sum())} 个作图细胞不在矩阵里")
    cl = np.full(len(obs_bc), "", dtype=object)
    cl[ci] = clu[col].astype(str).to_numpy()
    is_epi = cl != ""
    log(f"矩阵里找到作图细胞 {is_epi.sum():,} / 清单 {len(sub_bc):,}")
    if int(is_epi.sum()) != len(sub_bc):
        raise SystemExit("🔴 作图细胞数 ≠ 子集清单，拒绝继续")

    ucl = sorted(set(cl[is_epi]), key=lambda x: (len(x), x))

    # 每个细胞的归一化尺度（log1p 前的分母）
    scale = TARGET_SUM / np.maximum(ncount, 1.0)

    rows = []
    for c in ucl:
        m = (cl == c)
        n_in, n_out = int(m.sum()), int((~m & is_epi).sum())
        targ = amap.get(str(c), "")
        want = [("阳性对照", g) for g in POS_CTRL.get(targ, []) if g in g2j] + \
               [(f"待查:{lin}", g) for lin, gs in SUSPECT.items() for g in gs if g in g2j]
        for kind, g in want:
            j = g2j[g]
            s = (hit_gene == j)
            cells = hit_cell[s]
            # ⚠️ 一律用**线性**的归一化计数（每 1e4 的拷贝数），**不取 log**。
            # 2026-09-21 修正：早先版本先 log1p 再取均值，又把两个 log 尺度的数相除当倍数，
            # 那是错的（log1p(a)/log1p(b) ≠ a/b）。倍数必须在线性尺度上算。
            vals = hit_val[s] * scale[cells]
            vin = vals[m[cells]]
            vout = vals[(~m)[cells] & is_epi[cells]]
            p_in = len(vin) / n_in
            p_out = (len(vout) / n_out) if n_out else 0.0
            # 簇内**所有**细胞的均值（未表达的记 0）—— 与 scanpy 的 mean 同口径
            mean_in = float(vin.sum()) / n_in if n_in else 0.0
            mean_out = float(vout.sum()) / n_out if n_out else 0.0
            rows.append(dict(
                cluster=str(c), argmax=targ, n_cells=n_in, kind=kind, gene=g,
                pct_in=round(p_in, 4), pct_out=round(p_out, 4),
                pct_ratio=round(p_in / p_out, 2) if p_out > 0 else np.inf,
                mean_in=round(mean_in, 4), mean_out=round(mean_out, 4),
                mean_in_pos=round(float(vin.mean()), 4) if len(vin) else 0.0,
                log2r=round(float(np.log2((mean_in + 1e-6) / (mean_out + 1e-6))), 3),
            ))
    d = pd.DataFrame(rows)
    out = f"{OUT}/{a.tag}_ambient_vs_real.csv"
    d.to_csv(out, index=False)
    log(f"写出 {os.path.relpath(out, ROOT)}（{len(d)} 行）")

    # ---- 决定性对照：同一个基因，在**真细胞**里是什么水平 ----
    # 环境 RNA 与真表达的差别不在"有没有"，在"量级"。拿 CellTypist 判为该谱系的
    # 全细胞（同一份数据、同一套流程）当标尺：上皮簇里的量若只是真细胞的一小截，
    # 那就是游离 mRNA 的底噪，不是这些核自己在转录。
    ref_rows = []
    lab = pd.read_csv(LABELS, usecols=["cell_barcode", "B_lineage", "B_celltypist_type"])
    li = order.get_indexer(lab["cell_barcode"].to_numpy())
    ok = li >= 0
    lin_all = np.full(len(obs_bc), "", dtype=object)
    lin_all[li[ok]] = lab["B_lineage"].to_numpy()[ok]
    typ_all = np.full(len(obs_bc), "", dtype=object)
    typ_all[li[ok]] = lab["B_celltypist_type"].fillna("").to_numpy()[ok]

    for lin, gs in SUSPECT.items():
        mc = (lin_all == lin)
        if mc.sum() == 0:
            continue
        for g in gs:
            if g not in g2j:
                continue
            s = (hit_gene == g2j[g])
            cells = hit_cell[s]
            vals = hit_val[s] * scale[cells]          # 线性，不取 log（同上）
            v = vals[mc[cells]]
            pct_ref = len(v) / int(mc.sum())
            mean_ref = float(v.sum()) / int(mc.sum()) if len(v) else 0.0
            epi = d[(d["gene"] == g) & (d["kind"] == f"待查:{lin}")]
            if not len(epi):
                continue
            top = epi.loc[epi["mean_in"].idxmax()]
            ref_rows.append(dict(
                gene=g, lineage=lin,
                n_cells_lineage=int(mc.sum()),
                pct_in_true=round(pct_ref, 4),
                mean_in_true=round(mean_ref, 4),
                top_epi_cluster=str(top["cluster"]), top_epi_argmax=top["argmax"],
                top_epi_pct=top["pct_in"], top_epi_mean=top["mean_in"],
                fold_epi_over_true=round(float(top["mean_in"]) / mean_ref, 3) if mean_ref > 0 else np.inf,
            ))
    ref = pd.DataFrame(ref_rows)
    out_ref = f"{OUT}/{a.tag}_foreign_in_true_lineage.csv"
    ref.to_csv(out_ref, index=False)
    log(f"写出 {os.path.relpath(out_ref, ROOT)}")
    print(f"\n\n{'#'*104}\n# 决定性对照：上皮簇里的量 vs **真细胞**里的量（CellTypist 判为该谱系的全细胞）\n"
          f"# 读法：fold ≪ 1 ⇒ 上皮簇里那点量只是游离 mRNA 的底噪（环境 RNA），不是这些核在转录\n{'#'*104}")
    print(ref.to_string(index=False))

    pd.set_option("display.width", 200)
    cols = ["gene", "pct_in", "pct_out", "pct_ratio", "mean_in", "mean_out", "log2r"]
    for c in ucl:
        s = d[d["cluster"] == c]
        su = s[s["kind"].str.startswith("待查")]
        if not len(su):
            continue
        pc = s[s["kind"] == "阳性对照"]
        print(f"\n{'='*104}\n簇 {c}（{amap.get(str(c),'')}，{int(s['n_cells'].iloc[0]):,} 核）")
        if len(pc):
            print("  【真表达的标尺 —— 本簇自己的 marker】")
            print("  " + pc[cols].to_string(index=False).replace("\n", "\n  "))
        print("  【待查的外来基因】")
        print("  " + su[cols].to_string(index=False).replace("\n", "\n  "))

    json.dump(dict(
        tag=a.tag, lineage=a.lineage, res=a.res, seed=a.seed, cluster_col=col,
        n_cells=len(sub_bc), n_clusters=len(ucl),
        suspect_genes=SUSPECT, positive_controls=POS_CTRL,
        genes_missing_in_matrix=sorted(set(missing)),
        rule="阳性对照给'真'的标尺、pct_out 给'背景'的标尺；外来基因落在哪边就是哪边。"
             "**不设阈值** —— 看形状，由人工判读（与用户'人工定'口径一致）。",
        impl_note="按列取基因，不加载整个 h5ad；pct 用原始计数 >0 判表达（与归一化无关）；"
                  "mean = 按 nCount 归一化到 1e4 的**线性**均值，未表达的细胞记 0；"
                  "倍数一律在线性尺度上算。",
        correction="2026-09-21 修正：本脚本首版把 mean 先 log1p 再取均值、又拿两个 log 尺度的"
                   "数相除当倍数，那是错的（log1p(a)/log1p(b) ≠ a/b），会系统性放大倍数。"
                   "已改为线性尺度，`fore_lineage` 表与 `epiA_foreign_in_true_lineage.csv` 均为修正后口径。",
        reference_note="`foreign_in_true_lineage.csv` 是决定性对照：同一基因在 CellTypist 判为"
                       "该谱系的**真细胞**里的水平。fold = 上皮簇最高均值 / 真细胞均值。",
        outputs={os.path.relpath(p, ROOT): sha256(p) for p in [out, out_ref]},
        inputs={os.path.relpath(p, ROOT): sha256(p) for p in
                [H5, bc_file, clu_p, f"{OUT}/{a.tag}_cluster_annotation.csv", LABELS,
                 os.path.abspath(__file__)]},
    ), open(f"{OUT}/{a.tag}_ambient_vs_real_manifest.json", "w", encoding="utf-8"),
        ensure_ascii=False, indent=2)
    log("写出 manifest")


if __name__ == "__main__":
    main()
