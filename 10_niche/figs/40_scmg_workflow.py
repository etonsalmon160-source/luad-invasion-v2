#!/usr/bin/env python3
# 40_scmg_workflow.py —— 靶点逆向臂总流程图（SCMG 基因层 + LINCS 化合物层 → 候选 → 三重验证）
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = "/home/eto/luad_v2/results/paper_figures"
INK, MUTE = "#1a1a1a", "#5a5a5a"
BLUE, PURP, ROSE = "#2F5597", "#7E4FA8", "#C0392B"
TEAL, SAND = "#2E8B8B", "#B08A3E"

fig = plt.figure(figsize=(9.6, 2.9))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(20, 82); ax.axis("off")


def box(x, y, w, h, title, lines, edge, fill="white", tsize=6.8, lsize=5.6, lw=1.1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5,rounding_size=1.6",
                                linewidth=lw, edgecolor=edge, facecolor=fill, zorder=2))
    ax.text(x + w / 2, y + h - 3.4, title, ha="center", va="top", fontsize=tsize,
            fontweight="bold", color=edge, zorder=3)
    for i, t in enumerate(lines):
        ax.text(x + w / 2, y + h - 8.0 - i * 4.6, t, ha="center", va="top",
                fontsize=lsize, color=MUTE, zorder=3, linespacing=1.35)


def arrow(x0, x1, y):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=9,
                                 linewidth=1.0, color="#8a8a8a", zorder=1,
                                 shrinkA=0, shrinkB=0))


YB, HB, W = 50, 27, 16.4
XS = [1.5, 20.6, 39.7, 58.8, 77.9]

box(XS[0], YB, W, HB, "Single-cell atlas",
    ["399,579 nuclei", "25 patients", "Normal→AAH→AIS→MIA→IAC"], BLUE)
box(XS[1], YB, W, HB, "Query signature",
    ["disease axis: IAC − precursor", "patient-paired, 16,104 genes"], BLUE)
box(XS[2], YB, W, HB, "Perturbation libraries",
    ["SCMG  20,345 gene perturbations", "LINCS  1,826 compounds"], PURP)
box(XS[3], YB, W, HB, "Reversal scoring",
    ["SCMG → CausalGenePredictor", "  cos × sign × gene_shift_z",
     "LINCS → signatureSearch gess", "  CMap / Cor-Spearman"], PURP, lsize=5.2)
box(XS[4], YB, W, HB, "Ranked hypotheses",
    ["causal genes (systems level)", "34 compounds", "PI3K/mTOR-centric"], ROSE)

for i in range(4):
    arrow(XS[i] + W + 0.9, XS[i + 1] - 0.9, YB + HB / 2)

# 底部验证带
ax.add_patch(FancyBboxPatch((1.5, 22.0), 96.9, 20.0,
                            boxstyle="round,pad=0.5,rounding_size=1.6",
                            linewidth=1.0, edgecolor=TEAL, facecolor="#F4F8F8", zorder=2))
ax.text(3.2, 37.4, "Validation", ha="left", va="top", fontsize=6.6, fontweight="bold",
        color=TEAL, zorder=3)
VALS = [("Positive control", "gastrulation TF screen\nTBXT #1 of 8,450"),
        ("Matched null", "1,000 × random gene sets\nsignature level"),
        ("Split-half", "patients split ×4\nρ = 0.962  (null −0.155)"),
        ("Cross-modal", "single-cell vs spatial axis\nSpearman 0.735")]
for i, (h, t) in enumerate(VALS):
    x = 4.5 + i * 23.6
    ax.text(x, 32.6, h, ha="left", va="top", fontsize=5.8, fontweight="bold", color=INK, zorder=3)
    ax.text(x, 29.2, t, ha="left", va="top", fontsize=5.0, color=MUTE, zorder=3,
            linespacing=1.4)
for i in range(3):
    ax.plot([4.5 + (i + 1) * 23.6 - 2.4] * 2, [22.5, 35.0], color="#c8d8d8",
            linewidth=0.7, zorder=3)

fig.savefig(f"{OUT}/P10_scmg_workflow.pdf", bbox_inches="tight", facecolor="white")
fig.savefig(f"{OUT}/P10_scmg_workflow.png", dpi=400, bbox_inches="tight", facecolor="white")
print("✓ P10_scmg_workflow")
