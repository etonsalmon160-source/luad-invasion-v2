#!/usr/bin/env python3
"""`TREM2+ Dendritic` 这一个标签底下到底装了什么（只读已入库产物，不重算聚类）。

背景：S4 面板下 `TREM2+ Dendritic` 赢了 11 个簇 / 20,571 核，是髓系里最大的一块。
但它并不是一群细胞 —— 各簇**自己的** top-15 富集基因把它们分成五类，其中三类不是
这个细胞类型（含两个被预注册规则放行的污染簇）。

🔴 诚实边界（必读）：
  下面 `VERDICT` 里的分类**是我逐个读各簇 top-15 基因做的判读**，属**人工观察与推断**，
  不是任何自动化检验的结果，也没有独立验证。图里每个色块都标了 `top genes` 就是让判读
  可被你自己核对/推翻。**不得**把这个分类当成计算产物引用。

输入（全部已在库、带 sha256）：
  results/05_annotation/myeloidA_s4_cluster_annotation.csv   逐簇 argmax / n_cells / excluded
  results/05_annotation/myeloidA_s4_cluster_topgenes.csv     逐簇 top-15 富集基因
  results/05_annotation/myeloidA_contamination.csv           双线污染判据逐簇留痕
  results/04_integration/seurat_trad/myeloidA/clusters.csv.gz 逐核 patient/sample/stage

输出：figures/trem2_label_triage.png（无 CJK 字体 ⇒ 全英文标签）
      results/05_annotation/trem2_label_triage_manifest.json
"""
import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "05_annotation")
FIG = os.path.join(ROOT, "figures")

LABEL = "TREM2+ Dendritic"

# 🔴 人工判读（非自动化检验）：依据 = 各簇自身 top-15 基因 + 污染判据留痕
# 值 = (类别key, 判读理由里的关键基因)
VERDICT = {
    3:  ("lam",  "APOE GPNMB PSAP PLD3 LIPA CTSD"),
    10: ("lam",  "SPP1 CTSB GPNMB TREM2 APOE CD68"),
    12: ("lam",  "CHIT1 GPNMB CYP27A1 LIPA ACP5 PLA2G7"),
    15: ("lam",  "CHI3L1 LIPA GM2A APOE CTSS GPNMB"),
    17: ("lam",  "MMP9 MMP14 ADAMDEC1 CTSB PLA2G2D CTSZ"),
    5:  ("amb",  "TNFAIP2 ITGAX NPL GPNMB (ITGAX leans DC)"),
    7:  ("ifn",  "CXCL10 GBP1 CXCL9 GBP5 STAT1 TAP1"),
    14: ("art",  "CD83 CXCL8 CCL3 FOSB ATF3 NR4A2 DUSP1"),
    6:  ("cont", "SFTPB + COL1A1 (AT2 x fibroblast doublet)"),
    24: ("cont", "IGKC IGHG1 JCHAIN + COL1A1 (IG ambient)"),
    25: ("cont", "TRAC TRBC2 CD2 ZAP70 (T cell)"),
}
CAT = {
    "lam":  ("Real lipid-associated macrophage (LAM) program", "#2e86c1"),
    "amb":  ("Ambiguous (mixed program)",                      "#95a5a6"),
    "ifn":  ("IFN-response STATE, not a cell type",            "#8e44ad"),
    "art":  ("Dissociation / immediate-early artifact",        "#e67e22"),
    "cont": ("Contaminant, HELD by the pre-registered rule",   "#c0392b"),
}
ORDER = ["lam", "amb", "ifn", "art", "cont"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    os.makedirs(FIG, exist_ok=True)

    ann = pd.read_csv(os.path.join(R, "myeloidA_s4_cluster_annotation.csv"))
    con = pd.read_csv(os.path.join(R, "myeloidA_contamination.csv")).set_index("cluster")
    cl = pd.read_csv(os.path.join(ROOT, "results", "04_integration", "seurat_trad",
                                  "myeloidA", "clusters.csv.gz"),
                     usecols=["sample_id", "patient_id", "stage", "harmony_res0.5_seed0"])
    cl = cl.rename(columns={"harmony_res0.5_seed0": "cluster"})

    base_p24 = float((cl.patient_id == "P24").mean())
    base_p24luad = float((cl.sample_id == "P24_LUAD").mean())

    mine = ann[(ann.argmax == LABEL) & (~ann.excluded.astype(bool))].set_index("cluster")
    clusters = sorted(mine.index, key=lambda c: (ORDER.index(VERDICT[c][0]), -mine.loc[c, "n_cells"]))

    rows = []
    for c in clusters:
        d = cl[cl.cluster == c]
        rows.append(dict(
            cluster=int(c), n=int(mine.loc[c, "n_cells"]),
            cat=VERDICT[c][0], genes=VERDICT[c][1],
            win_margin=float(mine.loc[c, "win_margin"]),
            n_patients=int(mine.loc[c, "n_patients"]),
            p24=float((d.patient_id == "P24").mean()),
            p24luad=float((d.sample_id == "P24_LUAD").mean()),
            iac=float((d.stage == "IAC").mean()),
            ct_main=str(con.loc[c, "ct_main"]), ct_main_frac=float(con.loc[c, "ct_main_frac"]),
            panel_top_other=str(con.loc[c, "panel_top_other"]),
            hold=str(con.loc[c, "single_line_hold"]),
        ))
    df = pd.DataFrame(rows)
    df["p24luad_fold"] = df.p24luad / base_p24luad

    # ---- 画图 ----
    fig = plt.figure(figsize=(19.5, 8.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 0.85, 0.85], wspace=0.42)

    y = np.arange(len(df))
    labels = [f"cl{int(r.cluster)}" for r in df.itertuples()]
    colors = [CAT[r.cat][1] for r in df.itertuples()]

    # (A) 每簇核数 + 自身 top 基因（判读依据）
    ax = fig.add_subplot(gs[0, 0])
    ax.barh(y, df.n, color=colors, edgecolor="white", linewidth=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=10); ax.invert_yaxis()
    ax.set_xlabel("Nuclei assigned to the label by argmax", fontsize=10)
    ax.set_xlim(0, df.n.max() * 1.62)
    for i, r in enumerate(df.itertuples()):
        ax.text(r.n + 90, i, f"{r.n:,}   {r.genes}", va="center", fontsize=7.6, color="#333")
    ax.set_title(f"(A) The label `{LABEL}` won {len(df)} clusters / {int(df.n.sum()):,} nuclei\n"
                 f"but each cluster's OWN top markers say what it is",
                 fontsize=12, loc="left")
    ax.grid(axis="x", alpha=0.22)
    ax.legend(handles=[Patch(fc=CAT[k][1], label=CAT[k][0]) for k in ORDER],
              fontsize=8.4, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2,
              framealpha=0.96)

    # (B) P24_LUAD 富集倍数
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.barh(y, df.p24luad_fold, color=colors, edgecolor="white", linewidth=0.8)
    ax2.axvline(1.0, color="#111", ls="--", lw=1.4)
    ax2.text(1.06, len(df) - 0.35, "myeloid average\n(=1x)", fontsize=7.5, va="center", color="#111")
    ax2.set_yticks(y); ax2.set_yticklabels(labels, fontsize=10); ax2.invert_yaxis()
    ax2.set_xscale("log")
    ax2.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    ax2.set_xticklabels(["0.25x", "0.5x", "1x", "2x", "4x", "8x"])
    ax2.set_xlabel("Enrichment for P24_LUAD\n(vs all 64,084 myeloid nuclei)", fontsize=9.5)
    ax2.set_xlim(0.2, 11)
    for i, r in enumerate(df.itertuples()):
        ax2.text(r.p24luad_fold * 1.09, i, f"{r.p24luad_fold:.1f}x", va="center", fontsize=8)
    ax2.set_title("(B) The held contaminants are P24-driven\n"
                  "(up to 7.7x the myeloid average)", fontsize=12, loc="left")
    ax2.grid(axis="x", alpha=0.22)

    # (C) 每类占比
    ax3 = fig.add_subplot(gs[0, 2])
    tot = int(df.n.sum())
    agg = df.groupby("cat").n.sum().reindex(ORDER).fillna(0)
    yy = np.arange(len(agg))
    ax3.barh(yy, agg.values, color=[CAT[k][1] for k in agg.index],
             edgecolor="white", linewidth=0.8)
    ax3.set_yticks(yy)
    ax3.set_yticklabels([CAT[k][0].split("(")[0].split(",")[0].strip() for k in agg.index],
                        fontsize=9)
    ax3.invert_yaxis()
    ax3.set_xlim(0, agg.max() * 1.45)
    for i, v in enumerate(agg.values):
        ax3.text(v + tot * 0.012, i, f"{int(v):,}  ({v/tot:.0%})", va="center", fontsize=9)
    ax3.set_xlabel("Nuclei", fontsize=10)
    ax3.set_title(f"(C) Split of the {tot:,} nuclei\nby what they actually are",
                  fontsize=12, loc="left")
    ax3.grid(axis="x", alpha=0.22)
    ax3.text(0.76, 0.50,
             "Only ~52% is a real\nlipid-associated macrophage.\n"
             "17% is contaminant held by the\npre-registered 2-line rule:\n"
             "single-line hits are NOT\nauto-removed (see\n06_contamination_check.py:70).",
             transform=ax3.transAxes, ha="center", va="center", fontsize=7.0,
             bbox=dict(fc="#fdf2e9", ec="#e67e22", alpha=0.95))

    out = os.path.join(FIG, "trem2_label_triage.png")
    fig.savefig(out, dpi=155, bbox_inches="tight")

    inv = {
        "figure": os.path.relpath(out, ROOT),
        "figure_sha256": sha256(out),
        "question": (f"S4 面板下 `{LABEL}` 赢了 {len(df)} 簇 / {int(df.n.sum()):,} 核，"
                     f"是髓系最大一块。这 11 簇是不是同一群细胞？"),
        "answer": "不是。按各簇自身 top-15 富集基因，至少分五类，其中三类（含两个被预注册规则放行的污染簇）不是该细胞类型。",
        "🔴_honest_boundary": (
            "VERDICT 里的分类是**人工判读**（逐簇读 top-15 基因 + 污染判据留痕），"
            "属观察与推断，**不是自动化检验的结果，也没有独立验证**。"
            "图 (A) 右侧逐簇列出 top 基因就是为了让判读可被核对或推翻。"
            "**不得当作计算产物引用。**"),
        "category_definition": {k: v[0] for k, v in CAT.items()},
        "verdict_per_cluster": {int(r.cluster): {"category": r.cat, "n": r.n,
                                                 "top_genes_as_read": r.genes}
                                for r in df.itertuples()},
        "myeloid_baseline": {"n_cells": len(cl), "p24_frac": base_p24,
                             "p24luad_frac": base_p24luad,
                             "note": "P24_LUAD 占全部 64,084 髓系核的 8.4% ⇒ 富集倍数以此为 1x"},
        "per_cluster_evidence": {int(r.cluster): {
            "n_cells": r.n, "n_patients": r.n_patients, "win_margin": r.win_margin,
            "p24_frac": r.p24, "p24luad_frac": r.p24luad, "p24luad_fold": r.p24luad_fold,
            "iac_frac": r.iac,
            "celltypist_main": r.ct_main, "celltypist_main_frac": r.ct_main_frac,
            "panel_top_other": r.panel_top_other,
            "single_line_hold": r.hold} for r in df.itertuples()},
        "held_by_rule": {
            "clusters": [int(c) for c in df[df.cat == "cont"].cluster],
            "rule": ("06_contamination_check.py:70（冻结于看结果之前）：非上皮谱系"
                     "『须两条线同时命中才剔；单线命中如实留痕，但不自动剔』。"
                     "簇 6 与 25 是单线命中（仅 CellTypist）⇒ 按规则保留。"),
            "why_not_a_bug": ("规则先定后看结果，法则 3.2 禁止事后调参 ⇒ 处置是"
                              "**登记为已登记局限**，不是回改判据。")},
        "caveats": [
            "分类是人工判读，非自动检验（见 _honest_boundary）。",
            "簇内基因是 top-15 富集，不是差异检验的显著性；未经独立验证。",
            "`TREM2+ Dendritic` 是 Travaglini 健康肺命名；源论文（彭 2026）全文 TREM2 / SPP1 "
            "0 命中，其自身髓系空间 marker 写作 APOE / GPMNB 且称为 myeloid。",
            "本图不改变任何注释结果，只做判读可视化。",
        ],
        "inputs": {f: sha256(os.path.join(R, f)) for f in (
            "myeloidA_s4_cluster_annotation.csv", "myeloidA_s4_cluster_topgenes.csv",
            "myeloidA_contamination.csv")},
        "inputs_extra": {
            "results/04_integration/seurat_trad/myeloidA/clusters.csv.gz":
                sha256(os.path.join(ROOT, "results", "04_integration", "seurat_trad",
                                    "myeloidA", "clusters.csv.gz"))},
    }
    mf = os.path.join(R, "trem2_label_triage_manifest.json")
    json.dump(inv, open(mf, "w"), ensure_ascii=False, indent=1)

    print("图:", out)
    print("清单:", mf)
    print(df[["cluster", "n", "cat", "p24luad_fold", "ct_main"]].to_string(index=False))
    print("\n分类合计:")
    for k in ORDER:
        v = int(agg.get(k, 0))
        print(f"  {k:5s} {v:6,}  {v/tot:.1%}   {CAT[k][0]}")


if __name__ == "__main__":
    main()
