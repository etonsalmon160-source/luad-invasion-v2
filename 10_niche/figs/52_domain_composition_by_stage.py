#!/usr/bin/env python3
# 52_domain_composition_by_stage.py —— 七个生态位域的组分随病程变化
#   口径：**先算每张切片内的域构成，再按期别对切片取平均**
#        （若直接按 spot 汇总，张数多的期别会被组织体积主导）
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; T = f"{ROOT}/results/10_niche"
OUT = f"{ROOT}/results/paper_figures"
ST = ["Normal", "AAH", "AIS", "MIA", "IAC"]
ALIAS = {"LUAD": "IAC"}
DOMC = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
LAB = ["D1 Airway", "D2 iCAF", "D3 ECM/interstitial", "D4 Alveolar-cap.",
       "D5 AT2", "D6 Vascular", "D7 Lymphoid"]

D = pd.read_csv(f"{T}/kstar_diag/d7_domain_assign.tsv", sep="\t")
D["stage"] = D.stage.replace(ALIAS)
D = D[D.stage.isin(ST)]

# 切片内构成
D["tot"] = D.groupby("slide").n_spot.transform("sum")
D["frac"] = D.n_spot / D.tot
W = D.pivot_table(index=["slide", "stage"], columns="archetype", values="frac",
                  aggfunc="sum", fill_value=0).reindex(columns=range(1, 8), fill_value=0)
W = W.reset_index()
nsl = W.groupby("stage").size().reindex(ST)
C = W.groupby("stage")[list(range(1, 8))].mean().reindex(ST)          # 按期别对切片平均
print("逐期别的域构成（%，按切片平均）：")
print((C * 100).round(1).to_string())
print("\n各期切片数：", nsl.to_dict())

fig = plt.figure(figsize=(5.0, 3.2))
ax = fig.add_axes([0.11, 0.135, 0.60, 0.76])

xs = np.arange(len(ST))
bot = np.zeros(len(ST))
for k in range(1, 8):
    v = C[k].values * 100
    ax.bar(xs, v, bottom=bot, width=0.68, color=DOMC[k - 1], linewidth=0, label=LAB[k - 1])
    for i in range(len(ST)):
        if v[i] > 7:
            ax.text(i, bot[i] + v[i] / 2, f"{v[i]:.0f}", ha="center", va="center",
                    fontsize=5.8, color="white" if k != 4 else "0.35", fontweight="bold")
    bot += v
for i in range(len(ST)):
    ax.text(i, 101.5, f"n={nsl.iloc[i]}", ha="center", fontsize=5.9, color="0.45")
ax.set_xticks(xs); ax.set_xticklabels(ST, fontsize=7.4)
ax.set_ylim(0, 106); ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xlim(-0.6, 4.6)
ax.set_ylabel("Share of spots   (%)", fontsize=7.2)
ax.tick_params(axis="y", labelsize=6.6, length=2.4)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_linewidth(0.6); ax.spines[sp].set_color("0.45")
ax.set_title("Domain composition per stage", fontsize=8.6, fontweight="bold", pad=14, loc="left")

hs, ls = ax.get_legend_handles_labels()
# 图例右侧单列、与柱顶对齐
fig.legend(hs, ls, loc="upper left", bbox_to_anchor=(0.735, 0.885), frameon=False,
           fontsize=6.6, title="Niche domain", title_fontsize=6.9,
           labelspacing=0.62, handletextpad=0.55, handlelength=1.1)
fig.savefig(f"{OUT}/P22_domain_composition_by_stage.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P22_domain_composition_by_stage.png", dpi=400, bbox_inches="tight")
C.to_csv(f"{T}/domain_composition_by_stage.tsv", sep="\t")
print("✓ P22_domain_composition_by_stage")
