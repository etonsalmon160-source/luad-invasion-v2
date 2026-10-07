#!/usr/bin/env python3
# 34_domain_palette_preview.py —— 生态位域配色的候选方案对比（同一批切片）
import numpy as np, pickle, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CACHE = "/tmp/_p5_matrix_data.pkl"
REP = [("AIS", "GSM9226207_P19_AIS"), ("IAC", "GSM9226200_P15_LUAD")]
DATA = dict(zip([s for _, s in REP], pickle.load(open(CACHE, "rb"))[2:4]))
OUT = "/home/eto/luad_v2/results/paper_figures"

# 每列 = 一个方案，顺序都是 D1..D7；D4 占 45% 的 spot，决定整张图的底子
SCHEMES = {
    "E · blue-purple-pink\n(pale-pink background)":
        ["#4472A8", "#6C5CA8", "#A85CA0", "#F3CFE0", "#E8789E", "#3FA0A0", "#7E4FA8"],
    "E2 · blue-purple-pink\n(grey background)":
        ["#4472A8", "#6C5CA8", "#A85CA0", "#E8E8EA", "#E8789E", "#3FA0A0", "#7E4FA8"],
    "E3 · deeper family\n(pale-pink background)":
        ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"],
    "D · current":
        ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3"],
}
LAB = ["D1 Airway", "D2 iCAF", "D3 ECM/interstitial", "D4 Alveolar-cap.",
       "D5 AT2", "D6 Vascular", "D7 Lymphoid"]

fig, axes = plt.subplots(2, len(SCHEMES), figsize=(3.0 * len(SCHEMES), 6.4))
for r, (stage, sl) in enumerate(REP):
    D = DATA[sl]
    for c, (name, cols) in enumerate(SCHEMES.items()):
        ax = axes[r, c]
        ax.set_facecolor("black"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.4); s.set_color("0.35")
        a = D["dom"]; ok = ~pd.isna(a)
        ax.scatter(D["X"][ok], D["Y"][ok], c=[cols[int(v) - 1] for v in a[ok]],
                   s=2.2, marker="s", linewidths=0, rasterized=True)
        ax.set_xlim(0, D["W"]); ax.set_ylim(-D["H"], 0); ax.set_aspect("equal")
        if r == 0:
            ax.set_title(name, fontsize=7.5, fontweight="bold", color="white", linespacing=1.3)
        if c == 0:
            ax.set_ylabel(stage, fontsize=9, fontweight="bold", color="white")
    # 每行右侧列图例
    axes[r, -1].legend(handles=[Line2D([], [], marker="s", ls="", ms=5,
                                       mfc=SCHEMES["D · current"][i], mec="none", label=LAB[i])
                                for i in range(7)],
                       loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
                       fontsize=6.2, labelcolor="white", title="domain",
                       title_fontsize=6.5).get_title().set_color("white")

fig.patch.set_facecolor("black")
fig.tight_layout()
fig.savefig(f"{OUT}/_domain_palette_options.png", dpi=300, bbox_inches="tight",
            facecolor="black")
print("✓ _domain_palette_options")
