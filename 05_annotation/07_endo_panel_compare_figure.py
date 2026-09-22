#!/usr/bin/env python3
"""内皮面板 S1 vs S4 对比图（只读已入库产物，不重算任何东西）。

用途：给用户签字用。回答一个问题——把内皮的 marker 表从 Travaglini 2020 Table S1
的 canonical 摘要列换成同一篇论文 Table S4 的逐簇富集表，源数据本身的缺陷是否被修掉。

输入（全部已在库、带 sha256 的产物）：
  results/05_annotation/{endoA,endoA_s4}_cluster_annotation.csv   逐簇 argmax
  results/05_annotation/{endoA,endoA_s4}_cluster_scores.csv       逐簇逐型模块分
  results/05_annotation/panel_caliber_s4.json                      S4 逐型基因数
  results/05_annotation/classic_panels.py                          S1 逐型槽位数（N_SLOTS）

输出：figures/endo_panel_s1_vs_s4.png
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
    return dict(m.N_SLOTS["内皮"])


def main():
    os.makedirs(FIG, exist_ok=True)

    s1_ann = pd.read_csv(os.path.join(R, "endoA_cluster_annotation.csv"))
    s4_ann = pd.read_csv(os.path.join(R, "endoA_s4_cluster_annotation.csv"))
    s4_cal = json.load(open(os.path.join(R, "panel_caliber_s4.json")))
    s1_slots = load_s1_slots()

    # S1 逐型槽位数；S4 逐型可用基因数（n_used）
    s4_used = {t: int(v["n_used"]) for t, v in s4_cal["inventory"]["内皮"].items()}
    s4_scorable = {t: bool(v["scorable"]) for t, v in s4_cal["inventory"]["内皮"].items()}

    # 旧面板只有这 4 型；S4 有 9 型。并集，按 S4 顺序
    order = list(s4_cal["inventory"]["内皮"].keys())
    for t in s1_slots:
        if t not in order:
            order.append(t)

    # ---- (A) 逐型基因槽位：S1 源槽位 / S4 过筛后可用基因数 ----
    s1_vals = [s1_slots.get(t, 0) for t in order]
    s4_vals = [s4_used.get(t, 0) for t in order]

    fig = plt.figure(figsize=(15, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.28, 1.0], wspace=0.28)

    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(order))
    h = 0.38
    ax.barh(y + h / 2, s1_vals, height=h, color="#c0392b", label="Old panel (Table S1 canonical column)")
    ax.barh(y - h / 2, s4_vals, height=h, color="#2e86c1", label="New panel (same paper, Table S4 per-cluster)")
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(2, color="gray", ls=":", lw=1)
    ax.text(2.15, len(order) - 0.4, "2 slots", color="gray", fontsize=8)
    for yi, (a, b) in enumerate(zip(s1_vals, s4_vals)):
        if a:
            ax.text(a + 0.15, yi + h / 2, str(a), va="center", fontsize=8)
        ax.text(b + 0.15, yi - h / 2, str(b), va="center", fontsize=8,
                color="#1b4f72", weight="bold" if b else "normal")
    ax.set_xlabel("Number of marker genes available to score that subtype")
    ax.set_title("(A) Source marker table per endothelial subtype\n"
                 "Old: 6 genes total for 4 types  |  New: up to 20 genes, 9 types", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xlim(0, 24)
    ax.grid(axis="x", alpha=0.25)

    # ---- (B) 逐簇 argmax：谁赢过 ----
    ax2 = fig.add_subplot(gs[0, 1])
    w1 = s1_ann.loc[~s1_ann["excluded"], "argmax"].value_counts()
    w4 = s4_ann.loc[~s4_ann["excluded"], "argmax"].value_counts()
    labels = list(dict.fromkeys(list(w1.index) + list(w4.index)))
    v1 = [int(w1.get(l, 0)) for l in labels]
    v4 = [int(w4.get(l, 0)) for l in labels]
    x = np.arange(len(labels))
    ax2.bar(x - 0.19, v1, width=0.38, color="#c0392b", label="Old panel")
    ax2.bar(x + 0.19, v4, width=0.38, color="#2e86c1", label="New panel")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=42, ha="right", fontsize=8.5)
    ax2.set_ylabel("Number of clusters won (argmax, non-excluded)")
    ax2.set_title("(B) How many endothelial subtypes ever win a cluster\n"
                  "Old: only 2 of 4  |  New: 8 of 9 (Cap. Interm. 2 unscorable in source)", fontsize=11)
    ax2.legend(fontsize=9, loc="upper right", framealpha=0.95)
    ax2.grid(axis="y", alpha=0.25)
    ax2.text(0.985, 0.62, "same r*=0.7, seed 0,\n43,280 nuclei;\npanel is the only\nvariable",
             transform=ax2.transAxes, ha="right", va="top", fontsize=8,
             bbox=dict(fc="#fdf2e9", ec="#e67e22", alpha=0.95))

    out = os.path.join(FIG, "endo_panel_s1_vs_s4.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")

    inv = {
        "figure": os.path.relpath(out, ROOT),
        "figure_sha256": sha256(out),
        "comparison": "same object (endoA), same r*=0.7, same seed 0, same 43,280 nuclei; only the marker panel differs",
        "old_panel": {"source": "Travaglini 2020 Table S1 'Canonical markers' (summary column)",
                      "subtypes": list(s1_slots.keys()), "slots_per_type": s1_slots,
                      "subtypes_that_ever_win": sorted(w1.index.tolist()),
                      "n_subtypes_that_ever_win": int(len(w1))},
        "new_panel": {"source": "Travaglini 2020 Table S4 per-cluster enriched markers (same paper)",
                      "subtypes": list(s4_cal["inventory"]["内皮"].keys()),
                      "genes_used_per_type": s4_used, "scorable": s4_scorable,
                      "subtypes_that_ever_win": sorted(w4.index.tolist()),
                      "n_subtypes_that_ever_win": int(len(w4))},
        "pre_existing_flags_old": json.load(open(os.path.join(R, "endoA_annotation_manifest.json")))[
            "thin_panel_types"],
        "inputs": {f: sha256(os.path.join(R, f)) for f in (
            "endoA_cluster_annotation.csv", "endoA_s4_cluster_annotation.csv",
            "panel_caliber_s4.json")},
        "inputs_extra": {"05_annotation/classic_panels.py":
                         sha256(os.path.join(ROOT, "05_annotation", "classic_panels.py"))},
        "caveat": "Panel caliber panel_caliber_s4.json is provisional until countersigned; "
                  "so is every number on the New-panel side.",
    }
    mf = os.path.join(R, "endo_panel_s1_vs_s4_manifest.json")
    json.dump(inv, open(mf, "w"), ensure_ascii=False, indent=1)

    print("图:", out)
    print("清单:", mf)
    print("旧面板赢过的型 (%d): %s" % (len(w1), ", ".join(w1.index)))
    print("新面板赢过的型 (%d): %s" % (len(w4), ", ".join(w4.index)))


if __name__ == "__main__":
    main()
