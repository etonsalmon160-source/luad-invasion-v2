#!/usr/bin/env python3
# 43_fig_compound_spatial.py —— 化合物层·跨模态（单细胞轴 vs 空转轴）
#   同一个库、同一个打分口径（Cor-Spearman），只换查询轴 ⇒ 检验"结论是否依赖模态"
#   空转轴 = 逐切片按 nUMI 十分位分档、同档内比 IAC vs 前驱（74_/74b_ 产出）
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; S = f"{ROOT}/results/10_niche/sigsearch"
OUT = f"{ROOT}/results/paper_figures"

sc = pd.read_csv(f"{S}/res_COR_spearman.tsv", sep="\t").groupby("pert").cor_score.median()
sp = pd.read_csv(f"{S}/res_SPATIAL_cor.tsv", sep="\t").groupby("pert").cor_score.median()
J = pd.DataFrame({"sc": sc, "sp": sp}).dropna()
N = len(J); RK = J.rank()
CAND = list(pd.read_csv("/tmp/_comp34.tsv", sep="\t", index_col=0).index)
rho = J.sc.corr(J.sp, method="spearman")
ov100 = len(set(J.nsmallest(100, "sc").index) & set(J.nsmallest(100, "sp").index))
exp100 = 100 * 100 / N
print(f"n={N}  spearman={rho:.4f}  top100 重叠={ov100} (随机 {exp100:.1f})")

RED, GREY = "#C0392B", "#9aa0a6"
fig = plt.figure(figsize=(6.8, 3.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.26,
                      left=0.075, right=0.985, top=0.80, bottom=0.145)

# ── A：秩-秩 ──
ax = fig.add_subplot(gs[0, 0])
ax.add_patch(plt.Rectangle((1, 1), 99, 99, facecolor="#f6f1ea", edgecolor="none", zorder=0))
ax.axvline(100, color="0.72", lw=0.6, ls=(0, (3, 2.5)), zorder=1)
ax.axhline(100, color="0.72", lw=0.6, ls=(0, (3, 2.5)), zorder=1)
ax.scatter(RK.sc, RK.sp, s=2.6, c="0.78", linewidths=0, rasterized=True, zorder=2)
sub = RK.loc[[g for g in CAND if g in RK.index]]
ax.scatter(sub.sc, sub.sp, s=14, c=RED, edgecolor="white", linewidths=0.35, zorder=3)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.8, 2100); ax.set_ylim(0.8, 2100)
for a in (ax,):
    a.set_xticks([1, 10, 100, 1000]); a.set_yticks([1, 10, 100, 1000])
    a.set_xticklabels(["1", "10", "100", "1000"]); a.set_yticklabels(["1", "10", "100", "1000"])
ax.set_xlabel("rank by single-cell axis   (log)", fontsize=7)
ax.set_ylabel("rank by spatial axis   (log)", fontsize=7)
ax.tick_params(labelsize=6.2, length=2.2)
for s in ax.spines.values():
    s.set_linewidth(0.6); s.set_color("0.45")
ax.set_title(f"two independent query axes\nspearman = {rho:.2f}", fontsize=7.4,
             fontweight="bold", pad=5, loc="left")
ax.text(0.985, 0.035, f"shaded box = top 100 in both\n{ov100} compounds inside "
        f"(chance expectation ≈ {exp100:.1f})\nred = the 34 candidates from panel P12",
        transform=ax.transAxes, fontsize=5.6, color="0.30", va="bottom", ha="right",
        linespacing=1.55,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.6))

# ── B：topN 重叠富集 ──
axb = fig.add_subplot(gs[0, 1])
Ns = np.unique(np.round(np.logspace(1, np.log10(600), 40)).astype(int))
obs, exp = [], []
for n in Ns:
    obs.append(len(set(J.nsmallest(n, "sc").index) & set(J.nsmallest(n, "sp").index)))
    exp.append(n * n / N)
axb.plot(Ns, exp, color="0.55", lw=1.0, ls=(0, (4, 2.5)), label="chance expectation  (n²/N)")
axb.plot(Ns, obs, color=RED, lw=1.4, label="observed overlap")
axb.scatter([100], [ov100], s=22, facecolor="white", edgecolor=RED, linewidths=1.1, zorder=4)
axb.annotate(f"{ov100}", (100, ov100), xytext=(118, ov100 + 62), fontsize=6.4,
             color=RED, fontweight="bold",
             arrowprops=dict(arrowstyle="-", color=RED, lw=0.6, shrinkA=0, shrinkB=3))
axb.set_xscale("log")
axb.set_xlabel("top-N list compared   (log)", fontsize=7)
axb.set_ylabel("compounds shared by both top-N lists", fontsize=7)
axb.set_xlim(10, 600); axb.set_ylim(0, 330)
axb.set_xticks([10, 100, 500]); axb.set_xticklabels(["10", "100", "500"])
axb.tick_params(labelsize=6.2, length=2.2)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axb.spines[s].set_linewidth(0.6); axb.spines[s].set_color("0.45")
axb.legend(loc="upper left", frameon=False, fontsize=6.0, handlelength=1.8,
           labelspacing=0.3, handletextpad=0.5)
axb.set_title("the two axes agree far beyond chance", fontsize=7.4,
              fontweight="bold", pad=5, loc="left")
axb.grid(axis="y", color="0.92", linewidth=0.5); axb.set_axisbelow(True)

fig.text(0.065, 0.955, "The compound result is not modality-specific",
         fontsize=9, fontweight="bold", ha="left")
fig.savefig(f"{OUT}/P13_compound_reversal_crossmodal.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P13_compound_reversal_crossmodal.png", dpi=400, bbox_inches="tight")
print("✓ P13_compound_reversal_crossmodal")
