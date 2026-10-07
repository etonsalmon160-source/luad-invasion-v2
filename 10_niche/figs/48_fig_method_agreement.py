#!/usr/bin/env python3
# 48_fig_method_agreement.py —— 同一库内三个打分口径的两两关系（草稿 Fig 6 核心）
#   CMap（加权 KS，只用基因成员）／Cor-Spearman／Cor-Pearson（用数值）
import numpy as np, pandas as pd, itertools
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; S = f"{ROOT}/results/10_niche/sigsearch"
OUT = f"{ROOT}/results/paper_figures"

A = {}
for f, col in [("res_CMAP.tsv", "raw_score"), ("res_COR_spearman.tsv", "cor_score"),
               ("res_COR_pearson.tsv", "cor_score")]:
    A[f] = pd.read_csv(f"{S}/{f}", sep="\t").groupby("pert")[col].median()
M = pd.DataFrame(A); M.columns = ["CMap", "Cor-Spearman", "Cor-Pearson"]
N = len(M)
LAB = list(M.columns)

RHO = pd.DataFrame(np.eye(3), index=LAB, columns=LAB)
for a, b in itertools.combinations(LAB, 2):
    RHO.loc[b, a] = RHO.loc[a, b] = M[a].corr(M[b], method="spearman")
print(RHO.round(3).to_string())

fig = plt.figure(figsize=(7.0, 3.2))
gs = fig.add_gridspec(1, 2, width_ratios=[0.82, 1.0], wspace=0.30,
                      left=0.115, right=0.985, top=0.80, bottom=0.175)

# ── A：相关矩阵 ──
ax = fig.add_subplot(gs[0, 0])
cmap = plt.get_cmap("OrRd")
im = ax.imshow(RHO.values, cmap=cmap, vmin=0, vmax=1)
for i in range(3):
    for j in range(3):
        v = RHO.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.0,
                color="white" if v > 0.62 else "0.20",
                fontweight="bold" if i != j else "normal")
ax.set_xticks(range(3)); ax.set_xticklabels(LAB, fontsize=6.4, rotation=30, ha="right",
                                            rotation_mode="anchor")
ax.set_yticks(range(3)); ax.set_yticklabels(LAB, fontsize=6.4)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("Spearman", fontsize=6.2); cb.ax.tick_params(labelsize=5.6, length=2)
cb.outline.set_linewidth(0.4)
ax.set_title("agreement between the three scores", fontsize=7.6, fontweight="bold",
             pad=6, loc="left")

# ── B：top-N 重叠 ──
axb = fig.add_subplot(gs[0, 1])
COL = {"Cor-Spearman|Cor-Pearson": "#2F5597", "CMap|Cor-Spearman": "#C0392B",
       "CMap|Cor-Pearson": "#E08A3C"}
Ns = np.unique(np.round(np.logspace(1, np.log10(600), 45)).astype(int))
axb.plot(Ns, [n * n / N for n in Ns], color="0.55", lw=1.0, ls=(0, (4, 2.5)),
         label="chance", zorder=1)
for a, b in [("Cor-Spearman", "Cor-Pearson"), ("CMap", "Cor-Spearman"),
             ("CMap", "Cor-Pearson")]:
    key = f"{a}|{b}"
    ov = [len(set(M.nsmallest(n, a).index) & set(M.nsmallest(n, b).index)) for n in Ns]
    axb.plot(Ns, ov, color=COL[key], lw=1.4, label=f"{a} × {b}", zorder=3)
axb.set_xscale("log")
axb.set_xlabel("top-N list compared   (log)", fontsize=7)
axb.set_ylabel("compounds shared", fontsize=7)
axb.set_xlim(10, 600); axb.set_ylim(0, 480)
axb.set_xticks([10, 100, 500]); axb.set_xticklabels(["10", "100", "500"])
axb.tick_params(labelsize=6.4, length=2.5)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axb.spines[s].set_linewidth(0.6); axb.spines[s].set_color("0.45")
axb.legend(loc="upper left", frameon=False, fontsize=5.9, handlelength=1.6,
           labelspacing=0.26, handletextpad=0.5)
axb.set_title("shared top-N lists", fontsize=7.6, fontweight="bold", pad=6, loc="left")

fig.text(0.115, 0.955, "One library, three scores — and they only partly agree",
         fontsize=9, fontweight="bold", ha="left")
fig.text(0.115, 0.905, "The two correlation scores are near-redundant; the CMap statistic shares "
         "only about a third of its ranking with them.",
         fontsize=6.4, color="0.35", ha="left")
fig.savefig(f"{OUT}/P18_method_agreement.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P18_method_agreement.png", dpi=400, bbox_inches="tight")
print("✓ P18_method_agreement")

pd.DataFrame([dict(pair=f"{a} × {b}", spearman=M[a].corr(M[b], method="spearman"),
                   top50=len(set(M.nsmallest(50, a).index) & set(M.nsmallest(50, b).index)),
                   top500=len(set(M.nsmallest(500, a).index) & set(M.nsmallest(500, b).index)))
              for a, b in itertools.combinations(LAB, 2)]).to_csv(
    f"{ROOT}/results/10_niche/method_agreement.tsv", sep="\t", index=False)
