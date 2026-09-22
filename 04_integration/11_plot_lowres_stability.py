"""给"上皮/成纤维该不该整个作废"出一张判断用的图（决策辅助，不是产物）。

用户 2026-09-22 授权：把分辨率网格往下扩到 0.2/0.3/0.4，只重聚类、不动注释，
验证 0.5–0.8 全不达 0.90 是否只是"分辨率偏高"。本脚本只读两个新 run 的
`resolution_metrics.csv` 与已有的 r=0.7 注释表，不重算任何东西。

两个子图回答两个不同的问题：
  ① 稳定性到底是怎么掉的？—— 跨种子 ARI 对分辨率，含 0.90 门槛与各点簇数
  ② 那套不稳定的划分，注释是不是也一塌糊涂？—— r=0.7 的 27 个簇塌缩成几个型、
     赢家-亚军差有多大（用来区分"划分不稳"与"注释也不可信"）

用法: python3 04_integration/11_plot_lowres_stability.py
输出: results/05_annotation/figures/lowres_stability_decision.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = "/home/eto/luad_v2"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ROOT}/results/05_annotation/figures"
THRESH = 0.90
GRID = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
SERIES = [("epiA_lowres", "Epithelium (epiA, 133,384 nuclei)", "tab:red", "o"),
          ("fibroA_lowres", "Fibroblast (fibroA, 77,633 nuclei)", "tab:blue", "s")]

os.makedirs(OUT, exist_ok=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.6))

# ---- ① 稳定性曲线 --------------------------------------------------------
ax1.axhspan(0.70, THRESH, color="red", alpha=0.05)
ax1.axhline(THRESH, color="grey", ls="--", lw=1.6)
ax1.text(0.205, THRESH + 0.006, f"pre-registered gate  ARI >= {THRESH:.2f}",
         color="grey", fontsize=10.5, va="bottom")

for tag, label, color, mk in SERIES:
    m = pd.read_csv(f"{TRAD}/{tag}/resolution_metrics.csv").set_index("resolution")
    y = [m.loc[r, "ari_seed_mean"] for r in GRID]
    n = [int(m.loc[r, "n_clusters"]) for r in GRID]
    ok = [bool(m.loc[r, "pass_seed"]) for r in GRID]
    ax1.plot(GRID, y, "-", color=color, lw=2.2, zorder=3)
    for r, v, o in zip(GRID, y, ok):
        ax1.plot(r, v, mk, color=color, ms=11 if o else 9, zorder=4,
                 mfc=color if o else "white", mew=2.2)
        ax1.annotate(f"{v:.4f}\n{n[GRID.index(r)]} cl",
                     (r, v), textcoords="offset points",
                     xytext=(0, 13 if color == "tab:red" else -30),
                     ha="center", fontsize=9, color=color)
    ax1.plot([], [], mk, color=color, mfc=color, ms=9, label=label + "  [pass]")
    ax1.plot([], [], mk, color=color, mfc="white", ms=9, label=label + "  [fail]")

ax1.annotate("cliff 0.4 -> 0.5", xy=(0.45, 0.875), xytext=(0.52, 0.955),
             fontsize=11, color="black",
             arrowprops=dict(arrowstyle="->", color="black", lw=1.6))
ax1.set_xlabel("Leiden resolution", fontsize=12)
ax1.set_ylabel("cross-seed stability  (mean pairwise ARI over 5 seeds)", fontsize=12)
ax1.set_title("A. Stability collapses above r = 0.4\n"
              "(filled marker = passes gate, open = fails;  label = ARI / cluster count)",
              fontsize=12)
ax1.set_xticks(GRID)
ax1.set_ylim(0.70, 1.02)
ax1.grid(alpha=0.25)
ax1.legend(fontsize=8.6, loc="lower left", framealpha=0.92)

# ---- ② 注释是否也随之崩掉 -----------------------------------------------
a = pd.read_csv(f"{ROOT}/results/05_annotation/epiA_cluster_annotation.csv")
vc = a["argmax"].value_counts()
ax2.barh(range(len(vc))[::-1], vc.values, color="tab:red", alpha=0.75)
ax2.set_yticks(range(len(vc))[::-1])
ax2.set_yticklabels(vc.index, fontsize=11)
for i, (t, v) in enumerate(vc.items()):
    ax2.text(v + 0.25, len(vc) - 1 - i, f"{v} of {len(a)} clusters", va="center",
             fontsize=10, color="black")
ax2.set_xlim(0, max(vc.values) * 1.32)
ax2.set_xlabel("number of r = 0.7 clusters called as this type", fontsize=12)
ax2.set_title(f"B. Same partition, annotation side (epiA, r = 0.7)\n"
              f"{len(a)} unstable clusters collapse to {len(vc)} distinct types; "
              f"min winning margin = {a['win_margin'].min():.3f}",
              fontsize=12)
ax2.grid(alpha=0.25, axis="x")
ax2.invert_yaxis()

fig.tight_layout()
p = f"{OUT}/lowres_stability_decision.png"
fig.savefig(p, dpi=140)
print(f"写出 {os.path.relpath(p, ROOT)}")
