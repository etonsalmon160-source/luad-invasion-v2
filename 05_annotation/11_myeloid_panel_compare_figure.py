#!/usr/bin/env python3
"""髓系面板 S1 vs S4 对比图（只读已入库产物，不重算任何东西）。

用途：给用户复核用。回答一个问题——把髓系的 marker 表从 Travaglini 2020 Table S1
的 canonical 摘要列换成同一篇论文 Table S4 的逐簇富集表，源数据本身的缺陷是否被修掉。

🔴 两张表的**亚型名不一样、粒度也不一样**（例：S1 的 `Megakaryocyte` = S4 的
`Platelet/Megakaryocyte`；S1 一个 `Basophil/Mast` 在 S4 被拆成 1 和 2；S4 多出四个
`* Dendritic` / `Proliferating Macrophage` / `OLR1+ Classical Monocyte`）。
所以配对比较必须靠一张**我们自己的**命名对照表 `MAP`。**这张表不是原文的内容**，
是本项目为了画图做的映射，登记在 manifest 的 `name_mapping_ours` 里，不得当作原文口径引用。

输入（全部已在库、带 sha256 的产物）：
  results/05_annotation/{myeloidA,myeloidA_s4}_cluster_annotation.csv   逐簇 argmax
  results/05_annotation/{myeloidA,myeloidA_s4}_annotation_manifest.json 薄面板/不可打分类型
  results/05_annotation/panel_caliber_s4.json                           S4 逐型基因数
  results/05_annotation/classic_panels.py                               S1 逐型槽位数（N_SLOTS）
  results/05_annotation/myeloidA_rstar.json                             r*（确认两面板同一对象）

输出：figures/myeloid_panel_s1_vs_s4.png
无 CJK 字体 ⇒ 全英文标签。
"""
import hashlib
import importlib.util
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "05_annotation")
FIG = os.path.join(ROOT, "figures")

# 命名对照（本项目自定，非原文口径）：S1 型名 -> S4 型名
MAP = {
    "Neutrophil": "Neutrophil",
    "Eosinophil": "Eosinophil",
    "Megakaryocyte": "Platelet/Megakaryocyte",
    "Macrophage": "Macrophage",
    "pDC": "Plasmacytoid Dendritic",
    "mDC1": "Myeloid Dendritic Type 1",
    "mDC2": "Myeloid Dendritic Type 2",
    "Classical Monocyte": "Classical Monocyte",
    "Intermediate Monocyte": "Intermediate Monocyte",
    "Nonclassical Monocyte": "Nonclassical Monocyte",
    "Basophil/Mast": "Basophil/Mast 1",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_s1_slots():
    p = os.path.join(ROOT, "05_annotation", "classic_panels.py")
    spec = importlib.util.spec_from_file_location("_cp_ro", p)
    m = importlib.util.module_from_spec(spec)
    sys.modules["_cp_ro"] = m
    spec.loader.exec_module(m)
    return dict(m.N_SLOTS["髓系"])


def main():
    os.makedirs(FIG, exist_ok=True)

    s1_ann = pd.read_csv(os.path.join(R, "myeloidA_cluster_annotation.csv"))
    s4_ann = pd.read_csv(os.path.join(R, "myeloidA_s4_cluster_annotation.csv"))
    s1_mf = json.load(open(os.path.join(R, "myeloidA_annotation_manifest.json")))
    s4_mf = json.load(open(os.path.join(R, "myeloidA_s4_annotation_manifest.json")))
    s4_cal = json.load(open(os.path.join(R, "panel_caliber_s4.json")))
    rstar = json.load(open(os.path.join(R, "myeloidA_rstar.json")))
    s1_slots = load_s1_slots()

    s4_used = {t: int(v["n_used"]) for t, v in s4_cal["inventory"]["髓系"].items()}
    s4_scorable = {t: bool(v["scorable"]) for t, v in s4_cal["inventory"]["髓系"].items()}

    # 配对行（按 S1 顺序）+ S4 独有的行
    paired = [s4 for s1, s4 in MAP.items() if s1 in s1_slots]
    s4_only = [t for t in s4_cal["inventory"]["髓系"].keys() if t not in paired]
    order = paired + s4_only
    s1_slot_of = {MAP[s1]: n for s1, n in s1_slots.items() if s1 in MAP}
    s1_name_of = {s4: s1 for s1, s4 in MAP.items()}

    s1_vals = [s1_slot_of.get(t, 0) for t in order]
    s4_vals = [s4_used.get(t, 0) for t in order]

    # 逐型赢得的簇数（排除已排除簇）
    w1 = s1_ann.loc[~s1_ann["excluded"].astype(bool), "argmax"].value_counts()
    w4 = s4_ann.loc[~s4_ann["excluded"].astype(bool), "argmax"].value_counts()
    # 赢家表是 S1 型名，映射到 S4 名以便配对
    w1_m = {}
    for nm, c in w1.items():
        w1_m[MAP.get(nm, nm)] = w1_m.get(MAP.get(nm, nm), 0) + int(c)
    v1 = [int(w1_m.get(t, 0)) for t in order]
    v4 = [int(w4.get(t, 0)) for t in order]

    # ---- 画图 ----
    fig = plt.figure(figsize=(20.5, 6.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.72], wspace=0.30)

    # (A) 逐型可用基因数
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(order))
    h = 0.38
    ax.barh(y + h / 2, s1_vals, height=h, color="#c0392b",
            label="Old panel (Table S1 canonical column)")
    ax.barh(y - h / 2, s4_vals, height=h, color="#2e86c1",
            label="New panel (same paper, Table S4 per-cluster)")
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8.5)
    ax.invert_yaxis()
    ax.axvline(2, color="gray", ls=":", lw=1)
    ax.text(2.2, 0.15, "2 slots", color="gray", fontsize=8, va="center")
    ax.axvline(20, color="#1b4f72", ls="--", lw=1)
    ax.text(19.2, 0.15, "cap 20", color="#1b4f72", fontsize=8, ha="right", va="center")
    for yi, (a, b) in enumerate(zip(s1_vals, s4_vals)):
        if a:
            ax.text(a + 0.2, yi + h / 2, str(a), va="center", fontsize=7.5)
        ax.text(b + 0.2, yi - h / 2, str(b), va="center", fontsize=7.5,
                color="#1b4f72", weight="bold" if b else "normal")
    # 标出 S1 没有对应项的行
    for yi, t in enumerate(order):
        if t in s4_only:
            ax.text(22.6, yi, "not in old panel", fontsize=7, color="#7f8c8d",
                    va="center", style="italic")
    ax.set_xlabel("Number of marker genes available to score that subtype")
    ax.set_xlim(0, 30)
    ax.set_title(f"(A) Source marker table per myeloid subtype\n"
                 f"Old: {sum(s1_vals)} genes total for {sum(1 for v in s1_vals if v)} types"
                 f"  |  New: up to {max(s4_vals)} genes, {len(order)} types", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="x", alpha=0.25)

    # (B) 逐型赢得的簇数
    ax2 = fig.add_subplot(gs[0, 1])
    x = np.arange(len(order))
    ax2.bar(x - 0.19, v1, width=0.38, color="#c0392b", label="Old panel")
    ax2.bar(x + 0.19, v4, width=0.38, color="#2e86c1", label="New panel")
    ax2.set_xticks(x)
    ax2.set_xticklabels(order, rotation=42, ha="right", fontsize=8)
    ax2.set_ylabel("Number of clusters won (argmax, non-excluded)")
    n1 = sum(1 for t in order if v1 and order.index(t) < len(paired))
    n4 = sum(1 for c in v4 if c)
    ax2.set_title(f"(B) How many myeloid subtypes ever win a cluster\n"
                  f"Old: {sum(1 for c in v1 if c)} of {len(paired)}  |  "
                  f"New: {n4} of {len(order)}", fontsize=11)
    ax2.legend(fontsize=9, loc="upper right", framealpha=0.95)
    ax2.grid(axis="y", alpha=0.25)
    ax2.text(0.73, 0.72,
             f"same r*={rstar['r_star']}, seed {s1_mf['seed']},\n"
             f"{s1_mf['n_cells']:,} nuclei, {s1_mf['n_clusters']} clusters;\n"
             f"panel is the only variable",
             transform=ax2.transAxes, ha="center", va="top", fontsize=8,
             bbox=dict(fc="#fdf2e9", ec="#e67e22", alpha=0.95))

    # (C) 质量标记
    ax3 = fig.add_subplot(gs[0, 2])
    n_thin_cl1 = int(s1_ann["flag"].notna().sum())
    n_thin_cl4 = int(s4_ann["flag"].notna().sum())
    cats = ["clusters flagged\nthin_panel", "types NOT\nscorable", "types that\nnever win"]
    c1 = [n_thin_cl1, len(s1_mf["unscorable_types"]), len(s1_mf["never_won"])]
    c4 = [n_thin_cl4, len(s4_mf["unscorable_types"]), len(s4_mf["never_won"])]
    x3 = np.arange(len(cats))
    b1 = ax3.bar(x3 - 0.19, c1, width=0.38, color="#c0392b", label="Old panel")
    b4 = ax3.bar(x3 + 0.19, c4, width=0.38, color="#2e86c1", label="New panel")
    for b in list(b1) + list(b4):
        ax3.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.12,
                 str(int(b.get_height())), ha="center", fontsize=9)
    ax3.set_xticks(x3)
    ax3.set_xticklabels(cats, fontsize=8.5)
    ax3.set_ylabel("count")
    ax3.set_ylim(0, max(c1 + c4) + 2.2)
    ax3.set_title("(C) Quality markers\n"
                  "thin_panel flag = winning type has <=2 genes", fontsize=11)
    ax3.legend(fontsize=9, loc="upper right")
    ax3.grid(axis="y", alpha=0.25)
    ax3.text(0.52, 0.80,
             "zero flags is NOT\n'evidence is now enough':\n"
             "every scorable type merely\ncrossed the <=2 boundary\n"
             "(thinnest new panels: 3-4 genes)",
             transform=ax3.transAxes, ha="center", va="top", fontsize=7.3,
             bbox=dict(fc="#fdf2e9", ec="#e67e22", alpha=0.95))

    out = os.path.join(FIG, "myeloid_panel_s1_vs_s4.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")

    inv = {
        "figure": os.path.relpath(out, ROOT),
        "figure_sha256": sha256(out),
        "comparison": (f"same object (myeloidA), same r*={rstar['r_star']}, same seed "
                       f"{s1_mf['seed']}, same {s1_mf['n_cells']} nuclei, "
                       f"{s1_mf['n_clusters']} clusters; only the marker panel differs"),
        "name_mapping_ours": MAP,
        "name_mapping_caveat": ("🔴 这张对照表是本项目为画图自定的，不是 Travaglini 原文内容。"
                                "S1 与 S4 的亚型名与粒度都不同：S1 的 Megakaryocyte = S4 的 "
                                "Platelet/Megakaryocyte；S1 的 Basophil/Mast 在 S4 拆成 1 和 2；"
                                "S1 的 mDC1/mDC2/pDC 在 S4 写成全名。"),
        "old_panel": {"source": "Travaglini 2020 Table S1 'Canonical markers' (summary column)",
                      "n_types": len(s1_slots), "slots_per_type": s1_slots,
                      "total_genes": int(sum(s1_slots.values())),
                      "thin_panel_types": s1_mf["thin_panel_types"],
                      "unscorable_types": s1_mf["unscorable_types"],
                      "never_won": s1_mf["never_won"],
                      "clusters_flagged": int(s1_ann["flag"].notna().sum())},
        "new_panel": {"source": "Travaglini 2020 Table S4 per-cluster enriched markers (same paper)",
                      "n_types": len(s4_vals), "genes_used_per_type": s4_used,
                      "scorable": s4_scorable,
                      "thin_panel_types": s4_mf["thin_panel_types"],
                      "unscorable_types": s4_mf["unscorable_types"],
                      "never_won": s4_mf["never_won"],
                      "clusters_flagged": int(s4_ann["flag"].notna().sum())},
        "clusters_won_old_by_s4_name": {t: v1[i] for i, t in enumerate(order) if v1[i]},
        "clusters_won_new": {t: v4[i] for i, t in enumerate(order) if v4[i]},
        "thin_panel_flag_rule": "flag set when the WINNING type's panel has <=2 genes",
        "honest_caveat": ("新面板零 thin_panel 标记是**阈值边界**的结果（判据是 <=2 个基因），"
                          "不等于证据变充分：最薄的新面板只有 3–4 个基因"
                          "（Capillary Intermediate 1 = 3；Capillary / Bronchial Vessel 2 / "
                          "Neutrophil = 4）。"),
        "inputs": {f: sha256(os.path.join(R, f)) for f in (
            "myeloidA_cluster_annotation.csv", "myeloidA_s4_cluster_annotation.csv",
            "myeloidA_annotation_manifest.json", "myeloidA_s4_annotation_manifest.json",
            "panel_caliber_s4.json", "myeloidA_rstar.json")},
        "inputs_extra": {"05_annotation/classic_panels.py":
                         sha256(os.path.join(ROOT, "05_annotation", "classic_panels.py"))},
        "caveat": "Panel caliber panel_caliber_s4.json is provisional until countersigned; "
                  "so is every number on the New-panel side.",
    }
    mf = os.path.join(R, "myeloid_panel_s1_vs_s4_manifest.json")
    json.dump(inv, open(mf, "w"), ensure_ascii=False, indent=1)

    print("图:", out)
    print("清单:", mf)
    print(f"旧面板赢得簇的型 ({sum(1 for c in v1 if c)}):",
          ", ".join(t for t, c in zip(order, v1) if c))
    print(f"新面板赢得簇的型 ({sum(1 for c in v4 if c)}):",
          ", ".join(t for t, c in zip(order, v4) if c))
    print(f"thin_panel 簇数 旧 {n_thin_cl1} -> 新 {n_thin_cl4}")


if __name__ == "__main__":
    main()
