#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_qc/10_validate_doublets_crosslineage.py — 双体判定的【独立】验证：跨谱系共表达

原理（不依赖 scDblFinder）：
  真双体 = 两个细胞并入同一液滴 → 应出现**互斥谱系标记同时表达**。
  取一对互斥谱系：
     上皮  EPCAM   vs  免疫  PTPRC(CD45)
  对每个样本，在 QC 通过的细胞中比较：
     - 被判定为 doublet 的细胞中，**双阳性**（EPCAM>0 且 PTPRC>0）比例
     - 被判定为 singlet  的细胞中，同一比例
  期望：**doublet 的双阳性率显著高于 singlet**（若显著高 → 支持判定为真双体，
  而非"高计数的单细胞"）。

输出：results/01_qc/gse308103_doublet_crosslineage_validation.csv + 摘要
"""
import gzip, os
import numpy as np
import pandas as pd

RAW = "/home/eto/luad_invasion/data/GSE308103/extracted"
OUT = "/home/eto/luad_v2/results/01_qc"
MARKERS = ("EPCAM", "PTPRC")          # 上皮 vs 免疫（互斥）

qc = pd.read_csv(os.path.join(OUT, "gse308103_per_cell_qc.csv.gz"),
                 usecols=["sample_id", "cell_barcode", "qc_pass", "doublet_class"])
qc = qc[(qc.qc_pass) & (qc.doublet_class != "not_tested")]

files = sorted(f for f in os.listdir(RAW) if ".raw_counts" in f)
rows = []
for f in files:
    sid = f.split("_", 1)[1].split(".raw_counts")[0]
    sub = qc[qc.sample_id == sid]
    if not len(sub):
        continue
    # 读首行条码 + 仅取两个标记基因所在行（避免整矩阵）
    with gzip.open(os.path.join(RAW, f), "rt") as fh:
        bc = fh.readline().rstrip("\n").split("\t")
        vals = {}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] in MARKERS:
                vals[p[0]] = np.asarray(p[1:], dtype=np.float64)
            if len(vals) == len(MARKERS):
                break
    if len(vals) < len(MARKERS):
        print(f"  [{sid}] 缺标记，跳过"); continue
    a, b = vals[MARKERS[0]], vals[MARKERS[1]]
    dblpos = (a > 0) & (b > 0)

    m = pd.DataFrame({"cell_barcode": [f"{x}|{sid}" for x in bc], "dbl_pos": dblpos})
    m = m.merge(sub, on="cell_barcode", how="inner")
    d = m.doublet_class == "doublet"
    if d.sum() < 10 or (~d).sum() < 10:
        continue
    rows.append(dict(sample_id=sid, n=len(m), n_dbl=int(d.sum()),
                     dblpos_rate_doublet=100 * m.dbl_pos[d].mean(),
                     dblpos_rate_singlet=100 * m.dbl_pos[~d].mean()))

r = pd.DataFrame(rows)
r["ratio"] = (r.dblpos_rate_doublet / r.dblpos_rate_singlet.replace(0, np.nan)).round(2)
r.to_csv(os.path.join(OUT, "gse308103_doublet_crosslineage_validation.csv"), index=False)

print("=== 跨谱系共表达（EPCAM+ & PTPRC+）验证 ===")
print(f"  可评估样本: {len(r)}")
print(f"  doublet 双阳性率 中位: {r.dblpos_rate_doublet.median():.2f}%")
print(f"  singlet 双阳性率 中位: {r.dblpos_rate_singlet.median():.2f}%")
print(f"  doublet 更高的样本: {(r.dblpos_rate_doublet > r.dblpos_rate_singlet).sum()}/{len(r)} ({(r.dblpos_rate_doublet > r.dblpos_rate_singlet).mean()*100:.0f}%)")
print(f"  倍数比 中位: {r.ratio.median():.1f}x")
