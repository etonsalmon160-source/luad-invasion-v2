#!/usr/bin/env python3
# 46_fig_cohort_design.py —— Fig 1 队列与设计：两模态 × 五期 的患者级覆盖
#   单细胞 snRNA（GSE308103）／空间 Visium（GSE307534），同一批患者的配对取材
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
ST = ["Normal", "AAH", "AIS", "MIA", "IAC"]
ALIAS = {"LUAD": "IAC"}

# ── 单细胞：患者 × 期别 的细胞数 ──
C = pd.read_csv("/tmp/_cache_cellmap.csv")
C["st"] = C.stage.replace(ALIAS)
SC = C.pivot_table(index="patient_id", columns="st", values="cell_barcode",
                   aggfunc="size", fill_value=0).reindex(columns=ST, fill_value=0)

# ── 空间：患者 × 期别 的切片数 ──
d7 = pd.read_csv(f"{ROOT}/results/10_niche/kstar_diag/d7_domain_assign.tsv", sep="\t")
d7["st"] = d7.stage.replace(ALIAS)
SP = d7.drop_duplicates("slide").pivot_table(index="patient", columns="st", values="slide",
                                             aggfunc="size", fill_value=0).reindex(columns=ST, fill_value=0)

PATS = sorted(set(SC.index) | set(SP.index), key=lambda p: (int(p[1:]) if p[1:].isdigit() else 999))
print(f"患者 {len(PATS)}（snRNA {len(SC)} / Visium {len(SP)}）")

SC = SC.reindex(PATS).fillna(0); SP = SP.reindex(PATS).fillna(0)
TOT = pd.DataFrame({"snRNA 细胞": SC.sum(1), "Visium 切片": SP.sum(1)})
print(TOT.sum().to_string())

fig = plt.figure(figsize=(6.8, 6.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.42], height_ratios=[4.0, 1.0],
                      wspace=0.06, hspace=0.045, left=0.10, right=0.985, top=0.855, bottom=0.075)
ax = fig.add_subplot(gs[0, 0]); axr = fig.add_subplot(gs[0, 1], sharey=ax)
axb = fig.add_subplot(gs[1, 0], sharex=ax)

n = len(PATS)
for j, st in enumerate(ST):
    ax.axvline(j, color="0.93", lw=0.6, zorder=0)
ax.set_xlim(-0.5, len(ST) - 0.5); ax.set_ylim(-0.5, n - 0.5); ax.invert_yaxis()
for i, p in enumerate(PATS):
    for j, st in enumerate(ST):
        cs, cp = SC.loc[p, st], SP.loc[p, st]
        if cs > 0:
            ax.scatter(j - 0.16, i, s=6 + 46 * np.sqrt(cs / SC.values.max()),
                       c="#C0392B", linewidths=0, zorder=3)
        if cp > 0:
            ax.scatter(j + 0.16, i, s=9 + 30 * np.sqrt(cp / max(SP.values.max(), 1)),
                       c="#2F5597", marker="s", linewidths=0, zorder=3)
ax.set_xticks(range(len(ST)))
ax.tick_params(axis="x", labelbottom=False, length=0)   # 期别名交给底部面板
for s in ax.spines.values():
    s.set_visible(False)

ax.set_yticks(range(n))
ax.set_yticklabels(PATS, fontsize=5.4)
ax.tick_params(axis="y", length=0, pad=1)

# 右：每患者合计
axr.barh(range(n), TOT["Visium 切片"].values, height=0.34, color="#2F5597", zorder=3)
axr.barh(range(n), -TOT["snRNA 细胞"].values / 12000, height=0.34, color="#C0392B", zorder=3)
axr.axvline(0, color="0.55", lw=0.7)
axr.set_xlim(-22, 4); axr.tick_params(axis="y", labelleft=False, length=0); axr.tick_params(axis="x", labelsize=5.6, length=2)
axr.set_xticks([-20, -10, 0, 3]); axr.set_xticklabels(["240k", "120k", "0", "3"])
axr.set_xlabel("nuclei (left)   |   sections (right)", fontsize=6.0, labelpad=2)
for s in ("top", "right", "left"):
    axr.spines[s].set_visible(False)
axr.spines["bottom"].set_linewidth(0.6); axr.spines["bottom"].set_color("0.45")

# 下：每期别合计
w = 0.36
xs = np.arange(len(ST))
axb.bar(xs - w / 2, SP[ST].sum().values, width=w, color="#2F5597", label="Visium sections")
axb.bar(xs + w / 2, SC[ST].sum().values / 12000, width=w, color="#C0392B",
        label="snRNA nuclei ÷ 12k")
for j, st in enumerate(ST):
    axb.text(j - w / 2, SP[ST].sum().values[j] + 0.7, str(int(SP[ST].sum().values[j])),
             ha="center", fontsize=5.8, color="#2F5597")
    axb.text(j + w / 2, SC[ST].sum().values[j] / 12000 + 0.7,
             f"{int(SC[ST].sum().values[j]):,}", ha="center", fontsize=5.8, color="#C0392B")
axb.set_xticks(xs); axb.set_xticklabels(ST, fontsize=7.0)
axb.set_ylim(0, 31); axb.tick_params(axis="y", labelsize=5.8, length=2)
axb.set_ylabel("counts", fontsize=6.2)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axb.spines[s].set_linewidth(0.6); axb.spines[s].set_color("0.45")
axb.legend(loc="upper left", frameon=False, fontsize=5.8, handlelength=1.2,
           labelspacing=0.25, handletextpad=0.45, bbox_to_anchor=(0.005, 0.99))

fig.text(0.10, 0.975, "Cohort and study design", fontsize=10, fontweight="bold", ha="left")
fig.text(0.10, 0.950, f"{len(PATS)} patients;  {int(SC.values.sum()):,} nuclei and "
         f"{int(SP.values.sum())} Visium sections;  paired sampling across "
         "Normal → AAH → AIS → MIA → IAC", fontsize=6.6, color="0.35", ha="left")
fig.savefig(f"{OUT}/P16_cohort_design.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P16_cohort_design.png", dpi=400, bbox_inches="tight")
print("✓ P16_cohort_design")
