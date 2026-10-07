#!/usr/bin/env python3
# 39_domain_forest_plot.py —— 七个生态位域 marker 签名的预后森林图（OS ＋ DFS）
#   出两张独立图：
#     P9a_domain_forest_continuous —— 连续分数（HR per +1 SD）
#     P9b_domain_forest_lowhigh    —— 中位二分 Low vs High（参照 = Low）
#   数据：results/10_niche/m7a/m7a_niche_forest_lowhigh.tsv
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
D = pd.read_csv(f"{ROOT}/results/10_niche/m7a/m7a_niche_forest_lowhigh.tsv", sep="\t")

DOMC = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
DOMLAB = ["Airway (ciliated)", "iCAF", "ECM/interstitial", "Alveolar-capillary",
          "AT2", "Vascular", "Lymphoid"]
XLIM = 3.2
SIG, NS = "#C0392B", "0.35"


def pfmt(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def draw(mode, sublabel, fname, title):
    n = 7
    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 5, width_ratios=[2.15, 1.15, 0.60, 1.15, 0.60], wspace=0.10,
                          left=0.005, right=0.995, top=0.80, bottom=0.145)
    axL = fig.add_subplot(gs[0, 0]); axL.axis("off")
    fo = fig.add_subplot(gs[0, 1]); fd = fig.add_subplot(gs[0, 3])
    axPo = fig.add_subplot(gs[0, 2]); axPo.axis("off")
    axPd = fig.add_subplot(gs[0, 4]); axPd.axis("off")

    Y = np.arange(n)[::-1]
    for ax in (axL, axPo, axPd):
        ax.set_xlim(0, 1); ax.set_ylim(-0.7, n - 0.30)
    for f in (fo, fd):
        f.set_ylim(-0.7, n - 0.30); f.set_xlim(0, XLIM)
        f.axvline(1, color="0.45", linestyle=(0, (3, 2.5)), linewidth=0.8, zorder=0)
        for s in ("top", "right", "left"):
            f.spines[s].set_visible(False)
        f.spines["bottom"].set_linewidth(0.5); f.spines["bottom"].set_color("0.5")
        f.set_yticks([])
        f.tick_params(axis="x", labelsize=6, length=2, color="0.5", pad=2)
        f.set_xticks(np.arange(0, 4, 1)); f.set_xticklabels(["0", "1", "2", "3"])
        for t, v in zip(f.get_xticklabels(), np.arange(0, 4, 1)):
            if v == 1:
                t.set_fontweight("bold")
        f.set_xlabel("HRs", fontsize=6.5, labelpad=1.5)

    for ax, txt in ((fo, "OS"), (fd, "DFS"), (axPo, "P value"), (axPd, "P value")):
        ax.text(0.5, 1.03, txt, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=7, fontweight="bold")
    axL.text(0.0, 1.03, "Niche domain", transform=axL.transAxes, ha="left", va="bottom",
             fontsize=7, fontweight="bold")

    for i in range(7):
        y = Y[i]
        r = D[(D.domain == f"D{i+1}") & (D["mode"] == mode)].iloc[0]
        axL.add_patch(plt.Rectangle((0.005, y - 0.10), 0.055, 0.22, color=DOMC[i], lw=0))
        axL.text(0.085, y + 0.07, f"D{i+1}  {DOMLAB[i]}", ha="left", va="center",
                 fontsize=6.5, color="0.10")
        if sublabel:
            axL.text(0.105, y - 0.30, sublabel, ha="left", va="center", fontsize=6.0,
                     color="0.42", style="italic")
        for f, hrc, loc, hic, pc, axt in ((fo, "OS_HR", "OS_lo", "OS_hi", "OS_p", axPo),
                                          (fd, "DFS_HR", "DFS_lo", "DFS_hi", "DFS_p", axPd)):
            sig = r[pc] < 0.05
            col = SIG if sig else NS
            f.errorbar(r[hrc], y, xerr=[[r[hrc] - r[loc]], [r[hic] - r[hrc]]], fmt="o",
                       color=col, ecolor=col, elinewidth=0.9, capsize=1.8, capthick=0.9,
                       markersize=2.7 if sig else 2.2, zorder=3)
            axt.text(0.5, y, pfmt(r[pc]), ha="center", va="center", fontsize=6.3,
                     color=SIG if sig else "0.30", fontweight="bold" if sig else "normal")

    fig.text(0.5, 0.965, title, ha="center", va="top", fontsize=8.6, fontweight="bold")
    fig.savefig(f"{OUT}/{fname}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/{fname}.png", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {fname}")


draw("cont", "per +1 SD", "P9a_domain_forest_continuous",
     "Univariable Cox regression — niche-domain signatures (continuous)")
draw("lv", "Low vs High", "P9b_domain_forest_lowhigh",
     "Univariable Cox regression — niche-domain signatures (median split)")
