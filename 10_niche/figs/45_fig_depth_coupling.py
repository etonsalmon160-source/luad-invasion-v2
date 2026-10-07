#!/usr/bin/env python3
# 45_fig_depth_coupling.py —— 贯穿性发现：深度–密度耦联（草稿 §13 / Fig 4 核心面板）
#   七域各带一个深度倍数（1.00×–7.95×）；深度配平后，各域 top-150 up 基因的存活率 ∝ 1/深度
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
D = pd.read_csv(f"{ROOT}/results/10_niche/target_reversal/depth_overlap_all_domains.tsv", sep="\t")
D.columns = ["dom", "depth", "up_c", "dn_c", "old"]
D["surv"] = D.up_c / 150 * 100
D = D.sort_values("depth").reset_index(drop=True)

DOMC = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
NAME = {"D1": "Airway", "D2": "iCAF", "D3": "ECM/interstitial", "D4": "Alveolar-cap.",
        "D5": "AT2", "D6": "Vascular", "D7": "Lymphoid"}
COL = {d: DOMC[i] for i, d in enumerate(["D1", "D2", "D3", "D4", "D5", "D6", "D7"])}


def txt(c):
    """浅色（D4 的淡粉）直接当字色看不见，压暗后再用。"""
    r, g, b = [int(c[i:i + 2], 16) for i in (1, 3, 5)]
    if (0.299 * r + 0.587 * g + 0.114 * b) / 255 > 0.72:
        r, g, b = int(r * 0.45), int(g * 0.45), int(b * 0.45)
    return f"#{r:02x}{g:02x}{b:02x}"

x, y = D.depth.values, D.surv.values
rho = stats.spearmanr(x, y).correlation
b, a = np.polyfit(np.log10(x), y, 1)
xx = np.linspace(0.95, 8.4, 100)
yy = a + b * np.log10(xx)
res = y - (a + b * np.log10(x))
sd = res.std()
print(f"n={len(D)}  spearman={rho:.3f}  slope={b:.1f}  resid_sd={sd:.2f}")

fig, ax = plt.subplots(figsize=(5.6, 3.9))
ax.fill_between(xx, yy - 1.96 * sd, yy + 1.96 * sd, color="0.90", zorder=0, linewidth=0)
ax.plot(xx, yy, color="0.45", lw=1.0, zorder=1)
OFF = {"D4": (0, 11, "center"), "D2": (8, 3, "left"), "D6": (8, 3, "left"),
       "D1": (8, 2, "left"), "D7": (8, -4, "left"), "D5": (-9, -9, "right"),
       "D3": (0, -13, "center")}
for _, r in D.iterrows():
    c = COL[r.dom]
    ax.scatter(r.depth, r.surv, s=52, c=c, edgecolor="white", linewidths=0.8, zorder=4)
    dx, dy, ha = OFF[r.dom]
    ax.annotate(f"{r.dom} {NAME[r.dom]}", (r.depth, r.surv), xytext=(dx, dy),
                textcoords="offset points", fontsize=6.3, color=txt(c), ha=ha, va="center",
                fontweight="bold" if r.dom == "D3" else "normal", zorder=5,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.9))
ax.set_xscale("log")
ax.set_xticks([1, 1.5, 2, 3, 5, 8]); ax.set_xticklabels(["1", "1.5", "2", "3", "5", "8"])
ax.set_xlim(0.62, 9.4); ax.set_ylim(0, 66)
ax.set_xlabel("Domain depth relative to the shallowest domain   (×, log scale)", fontsize=7.2)
ax.set_ylabel("% of the domain's top-150 up genes\nsurviving depth matching", fontsize=7.2)
ax.tick_params(labelsize=6.5, length=2.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.6); ax.spines[s].set_color("0.45")
ax.text(0.03, 0.045, f"Spearman ρ = {rho:.2f}\nshaded = 95% band of the fit",
        transform=ax.transAxes, fontsize=6.0, color="0.35", linespacing=1.5)
ax.annotate("D3 is the sole outlier —\nand the only LUAD-only domain",
            (7.95, 10.7), xytext=(3.05, 44), fontsize=6.0, color="#9C4E93", linespacing=1.5,
            arrowprops=dict(arrowstyle="-|>", color="#9C4E93", lw=0.7,
                            connectionstyle="arc3,rad=-0.22"))
ax.set_title("A domain's signature stability is set by how deep that domain sits",
             fontsize=8.4, fontweight="bold", pad=8, loc="left")
fig.savefig(f"{OUT}/P15_depth_density_coupling.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P15_depth_density_coupling.png", dpi=400, bbox_inches="tight")
print("✓ P15_depth_density_coupling")
print(D[["dom", "depth", "surv"]].to_string(index=False))
