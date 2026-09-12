#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_qc/05_doublets_scrublet.py — M1 交叉验证：用 scrublet 独立复检双体

目的：scDblFinder 是主方法，但**单一方法会被质疑**。用 **scrublet**
      （Wolock 2019, Cell Systems；独立算法：模拟双体 + 期望率阈值）逐样本复检，
      报告**两法一致率**，支撑 M1 结论。

设定（与主流程一致）：
  * 同一套 QC 掩码（读 gse308103_per_cell_qc.csv.gz 的 qc_pass）
  * scrublet `expected_doublet_rate` 用**固定默认 0.06**（不采用 scDblFinder 的估计 →
    两法**相互独立**，否则一致率被人为抬高）
  * random_state=1（确定性）

输出：gse308103_scrublet_per_cell.csv.gz
      gse308103_doublet_method_agreement.csv
"""
import gzip, os
import numpy as np
import pandas as pd
from scrublet import Scrublet

RAW = "/home/eto/luad_invasion/data/GSE308103/extracted"
OUT = "/home/eto/luad_v2/results/01_qc"
EXPECTED_RATE = 0.06
N_PRIN_COMP, N_NEIGHBORS, RANDOM_STATE = 30, 30, 1

qc = pd.read_csv(os.path.join(OUT, "gse308103_per_cell_qc.csv.gz"),
                 usecols=["sample_id", "cell_barcode", "qc_pass"])
files = sorted(f for f in os.listdir(RAW) if ".raw_counts" in f)
print(f"[scrublet] {len(files)} 样本；expected_doublet_rate={EXPECTED_RATE}", flush=True)

rows = []
for i, f in enumerate(files, 1):
    sid = f.split("_", 1)[1].split(".raw_counts")[0]
    keep = set(qc.loc[(qc.sample_id == sid) & qc.qc_pass, "cell_barcode"])

    # 首行 = 条码（无基因列名），数据行 = 基因 + 值
    # → 单独取首行做列名；再 skiprows=1 读数据（index_col=0 = 基因）；pandas 直接读 .gz
    path = os.path.join(RAW, f)
    barcodes = pd.read_csv(path, sep="\t", header=None, nrows=1).iloc[0].tolist()
    # 注意：dtype 必须**只**给数据列（1..n），否则会试图把索引列（基因名）转成 float 而报错
    dt = {i: np.float32 for i in range(1, len(barcodes) + 1)}
    df = pd.read_csv(path, sep="\t", header=None, index_col=0, skiprows=1, dtype=dt)
    if len(df.columns) != len(barcodes):
        raise RuntimeError(f"{sid}: 条码数 {len(barcodes)} != 数据列数 {len(df.columns)}")
    mask = np.fromiter((f"{b}|{sid}" in keep for b in barcodes), dtype=bool, count=len(barcodes))
    Xq = df.values[:, mask]                     # genes x cells (float32)
    del df
    bcq = [b for b, m in zip(barcodes, mask) if m]

    if Xq.shape[1] < 50:
        print(f"  [{i:2d}/{len(files)}] {sid:12s} 跳过（QC 后 <50 细胞）", flush=True)
        continue

    sc = Scrublet(Xq.T, expected_doublet_rate=EXPECTED_RATE,
                  n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)
    scores, preds = sc.scrub_doublets()   # 返回 (doublet_scores_obs_, predicted_doublets_)
    rows.append(pd.DataFrame({
        "dataset": "GSE308103", "sample_id": sid,
        "cell_barcode": [f"{b}|{sid}" for b in bcq],
        "scrublet_score": np.round(scores, 4),
        "scrublet_class": np.where(preds, "doublet", "singlet"),
    }))
    print(f"  [{i:2d}/{len(files)}] {sid:12s} n={Xq.shape[1]:6d}  scrublet 双体="
          f"{int(preds.sum()):5d} ({100*preds.mean():5.1f}%)", flush=True)

sp = pd.concat(rows, ignore_index=True)
with gzip.open(os.path.join(OUT, "gse308103_scrublet_per_cell.csv.gz"), "wt", compresslevel=6) as fh:
    sp.to_csv(fh, index=False)

# ---- 两法一致性（仅 QC 通过、且两法都判过的细胞）----
mg = qc.merge(sp, on=["sample_id", "cell_barcode"], how="inner")
mg = mg[mg.qc_pass & (mg.doublet_class != "not_tested")].copy()
mg["sd_dbl"] = mg.doublet_class == "doublet"
mg["sb_dbl"] = mg.scrublet_class == "doublet"

def _s(d):
    both = int((d.sd_dbl & d.sb_dbl).sum()); sd = int((d.sd_dbl & ~d.sb_dbl).sum())
    sb = int((~d.sd_dbl & d.sb_dbl).sum()); ne = int((~d.sd_dbl & ~d.sb_dbl).sum())
    return pd.Series({"n": len(d), "sd_rate": round(100*d.sd_dbl.mean(), 2),
                      "sb_rate": round(100*d.sb_dbl.mean(), 2), "both": both, "sd_only": sd,
                      "sb_only": sb, "neither": ne,
                      "jaccard": round(both / max(both + sd + sb, 1), 3)})
agree = mg.groupby("sample_id").apply(_s).reset_index()
agree.to_csv(os.path.join(OUT, "gse308103_doublet_method_agreement.csv"), index=False)

n = int(agree.n.sum())
print("\n=== 两法一致性（全队列）===")
print(f"  比较细胞数        : {n:,}")
print(f"  两法都判双体      : {int(agree.both.sum()):,}")
print(f"  仅 scDblFinder 判 : {int(agree.sd_only.sum()):,}")
print(f"  仅 scrublet 判    : {int(agree.sb_only.sum()):,}")
print(f"  两法都判单细胞    : {int(agree.neither.sum()):,}")
print(f"  scDblFinder 双体率: {100*mg.sd_dbl.mean():.2f}%")
print(f"  scrublet    双体率: {100*mg.sb_dbl.mean():.2f}%")
print(f"  一致率            : {100*(agree.both.sum()+agree.neither.sum())/n:.2f}%")
print(f"  Jaccard 中位      : {agree.jaccard.median():.3f}")
