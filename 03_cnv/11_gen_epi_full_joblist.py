#!/usr/bin/env python3
"""生成「全量上皮 · 锚定臂」的**作业清单**（只生成清单，不跑）。

口径来源 = `ANCHOR_PREREG.md`（2026-09-22 预注册）+ `PARAMETERS_AND_SOURCES.md` §M2。
本脚本**不引入任何新参数**，只把已签字的定义翻译成一行一个 Rscript 调用。

已签字的三个范围常数（改任何一个都必须另开预注册）：
  CALIBER_SUBTYPES : 主口径 AT2、次口径 AT1 —— 见 prereg §2.4 那两行。**不含**少见亚型
                     （Ciliated/Goblet/Basal/Serous）——它们从未被登记为 query 亚型。
  MIN_ANCHOR       : 10 —— prereg G5，取自 copykat 自身 `baseline.norm.cl(min.cells=10)`
  FLOOR            : 840 —— 用户签字，作用于 M1 的 nFeature

P4 有两个 Normal 样本（P4_Normal / P4_Normal1）。预注册禁止静默消解：
  · 主锚 = AT2 核数较多者（P4_Normal），**清单里写出来**，不是隐式取 max；
  · 对**关键路径那一条**（P4_LUAD|AT2，全批最长跑）追加 P4_Normal1 的敏感性跑。
"""
import collections
import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CNV = os.path.join(ROOT, "results", "03_cnv")
EPI = os.path.join(ROOT, "results", "05_annotation", "epiCNV_subset_barcodes.txt")
QC = os.path.join(ROOT, "results", "01_qc", "gse308103_per_cell_qc.csv.gz")
CLU = os.path.join(ROOT, "results", "04_integration", "seurat_trad", "epiA", "clusters.csv.gz")
ANN = os.path.join(ROOT, "results", "05_annotation", "epiA_cluster_annotation.csv")
EST = os.path.join(CNV, "epi_full_cost_estimate.json")
OUT = os.path.join(CNV, "epi_full_joblist.tsv")

SEED_COL = "harmony_res0.7_seed0"
FLOOR = 840
MIN_ANCHOR = 10
CALIBER_SUBTYPES = ("AT2", "AT1")     # 主口径 + 次口径；prereg §2.4
SENSITIVITY = ("P4_LUAD", "AT2")      # 唯一追加第二锚定的 combo（见模块 docstring）


def main():
    m = json.load(open(EST))["report"]["model"]["t_sec"]
    a = float(m.split(" + ")[0])
    b = float(m.split(" + ")[1].split(" * ")[0])
    c = float(m.split("^")[1])
    sec = lambda n: a + b * (n ** c)

    epi = {l.strip() for l in open(EPI) if l.strip()}
    qc = pd.read_csv(QC, usecols=["sample_id", "cell_barcode", "nFeature"])
    df = qc[qc["cell_barcode"].isin(epi) & (qc["nFeature"] >= FLOOR)].copy()
    clu = pd.read_csv(CLU, usecols=["cell_barcode", SEED_COL])
    ann = pd.read_csv(ANN, usecols=["cluster", "argmax"])
    df = df.merge(clu, on="cell_barcode", how="left")
    df["st"] = df[SEED_COL].map(dict(zip(ann["cluster"], ann["argmax"])))
    df = df[df["st"].notna()]
    idx = df.groupby(["sample_id", "st"]).size()

    pat = lambda s: s.split("_")[0]
    normals = collections.defaultdict(list)
    for s in df["sample_id"].unique():
        normals[pat(s)].append(s)
    normals = {p: sorted(x for x in v if "Normal" in x) for p, v in normals.items()}

    rows, blocked = [], []
    for (s, st), nv in idx.items():
        if "Normal" in s or st not in CALIBER_SUBTYPES:
            continue
        cand = sorted(((int(idx.get((nm, st), 0)), nm) for nm in normals.get(pat(s), [])),
                      key=lambda t: (-t[0], t[1]))
        if not cand or cand[0][0] < MIN_ANCHOR:
            blocked.append((s, st, int(nv), cand[0][0] if cand else 0))
            continue
        an, aref = cand[0]
        rows.append(dict(subtype=st, sample_id=s, anchor_ref_sample=aref,
                         anchor_n=an, lesion_n=int(nv), n_total=int(nv) + an,
                         variant="", note="primary anchor"))
        # P4 的第二个 Normal：只对唯一登记的敏感性 combo 追加
        if (s, st) == SENSITIVITY and len(cand) > 1:
            an2, aref2 = cand[1]
            if an2 >= MIN_ANCHOR:
                rows.append(dict(subtype=st, sample_id=s, anchor_ref_sample=aref2,
                                 anchor_n=an2, lesion_n=int(nv), n_total=int(nv) + an2,
                                 variant="_anchorB",
                                 note="P4 双 Normal 敏感性跑（不静默消解）"))

    rows.sort(key=lambda r: -r["n_total"])
    with open(OUT, "w") as fh:
        fh.write("subtype\tsample_id\tanchor_ref_sample\tanchor_n\tlesion_n\tn_total"
                 "\tpred_h\tvariant\tnote\n")
        for r in rows:
            fh.write(f"{r['subtype']}\t{r['sample_id']}\t{r['anchor_ref_sample']}\t"
                     f"{r['anchor_n']}\t{r['lesion_n']}\t{r['n_total']}\t"
                     f"{sec(r['n_total'])/3600:.2f}\t{r['variant']}\t{r['note']}\n")

    print(f"清单: {os.path.relpath(OUT, ROOT)}")
    print(f"  可跑 {len(rows)} 跑（口径={'+'.join(CALIBER_SUBTYPES)}），"
          f"合计 {sum(sec(r['n_total']) for r in rows)/3600:.1f} CPU·h，"
          f"最长 {rows[0]['sample_id']}|{rows[0]['subtype']} {sec(rows[0]['n_total'])/3600:.2f} h")
    print(f"  清单内被挡（锚定 <{MIN_ANCHOR}）: {len(blocked)}")
    print("\n  前 10 跑:")
    for r in rows[:10]:
        print(f"    {r['sample_id']:<12}{r['subtype']:<5}病灶{r['lesion_n']:>6} "
              f"锚定{r['anchor_n']:>6} n={r['n_total']:>6} {sec(r['n_total'])/3600:>5.2f} h"
              f"  {r['anchor_ref_sample']}{r['variant']}")
    if blocked:
        print("\n  ⚠️ 被挡（不跑，按 2026-09-22 裁定记为 cnv_not_tested）:")
        for s, st, nv, an in sorted(blocked):
            print(f"    {s:<12}{st:<5}病灶{nv:>6} 锚定{an:>3}")


if __name__ == "__main__":
    main()
