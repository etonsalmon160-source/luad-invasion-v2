#!/usr/bin/env python3
# 47_fig_split_half_stability.py —— 化合物榜的患者劈半可重复性
#   4 次劈半：患者随机分两半 → 各自推疾病轴 → 各自对 LINCS 打分 → 比较两半的化合物排名
#   零模型：打乱其中一半的排名（置换），给出 ρ 的零分布
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; S = f"{ROOT}/results/10_niche/sigsearch/splits"
OUT = f"{ROOT}/results/paper_figures"
RNG = np.random.default_rng(20261006)

rows = []
for k in range(4):
    h1 = pd.read_csv(f"{S}/cor_split{k}_h1.tsv", sep="\t").groupby("pert").cor_score.median()
    h2 = pd.read_csv(f"{S}/cor_split{k}_h2.tsv", sep="\t").groupby("pert").cor_score.median()
    j = pd.DataFrame({"h1": h1, "h2": h2}).dropna()
    rho = j.h1.corr(j.h2, method="spearman")
    # 置换零模型
    nul = np.array([j.h1.corr(pd.Series(RNG.permutation(j.h2.values), index=j.index),
                              method="spearman") for _ in range(300)])
    rows.append(dict(k=k, rho=rho, nul_med=np.median(nul), nul_q975=np.quantile(nul, .975),
                     nul_q025=np.quantile(nul, .025),
                     o500=len(set(j.nsmallest(500, "h1").index) & set(j.nsmallest(500, "h2").index)),
                     n=len(j)))
R = pd.DataFrame(rows)
print(R.round(3).to_string(index=False))
print(f"观测 ρ 中位 {R.rho.median():.3f} | 零模型中位 {R.nul_med.median():+.3f} "
      f"| top500 观测中位 {int(R.o500.median())}  随机期望 {500*500/R.n.iloc[0]:.0f}")

fig = plt.figure(figsize=(7.4, 3.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.26,
                      left=0.075, right=0.985, top=0.80, bottom=0.155)

# ── A：逐次劈半的 ρ ──
ax = fig.add_subplot(gs[0, 0])
RED, GREY = "#C0392B", "0.62"
for _, r in R.iterrows():
    ax.plot([r.k, r.k], [r.nul_q025, r.nul_q975], color=GREY, lw=5, alpha=0.55,
            solid_capstyle="butt", zorder=2)
ax.scatter(R.k, R.rho, s=58, c=RED, edgecolor="white", linewidths=0.9, zorder=4)
ax.axhline(0, color="0.75", lw=0.7, zorder=1)
ax.set_xticks(range(4)); ax.set_xticklabels([f"split {i+1}" for i in range(4)], fontsize=6.8)
ax.set_ylim(-0.25, 1.05); ax.set_xlim(-0.55, 3.55)
ax.set_ylabel("Spearman  (half 1  vs  half 2)", fontsize=7)
ax.tick_params(labelsize=6.4, length=2.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.6); ax.spines[s].set_color("0.45")
ax.text(0.03, 0.05, "grey bar = permutation null (2.5–97.5%)",
        transform=ax.transAxes, fontsize=5.8, color="0.42")
ax.set_title("rank reproducibility across patient halves", fontsize=7.6,
             fontweight="bold", pad=6, loc="left")

# ── B：top-N 重叠 ──
axb = fig.add_subplot(gs[0, 1])
h1 = pd.read_csv(f"{S}/cor_split0_h1.tsv", sep="\t").groupby("pert").cor_score.median()
h2 = pd.read_csv(f"{S}/cor_split0_h2.tsv", sep="\t").groupby("pert").cor_score.median()
j = pd.DataFrame({"h1": h1, "h2": h2}).dropna(); N = len(j)
Ns = np.unique(np.round(np.logspace(1, np.log10(600), 40)).astype(int))
obs = [len(set(j.nsmallest(n, "h1").index) & set(j.nsmallest(n, "h2").index)) for n in Ns]
axb.plot(Ns, [n * n / N for n in Ns], color="0.55", lw=1.0, ls=(0, (4, 2.5)),
         label="chance expectation  (n²/N)")
axb.plot(Ns, obs, color=RED, lw=1.4, label="observed overlap")
axb.set_xscale("log")
axb.set_xlabel("top-N list compared   (log)", fontsize=7)
axb.set_ylabel("compounds shared by both halves", fontsize=7)
axb.set_xlim(10, 600); axb.set_ylim(0, 480)
axb.set_xticks([10, 100, 500]); axb.set_xticklabels(["10", "100", "500"])
axb.tick_params(labelsize=6.4, length=2.5)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axb.spines[s].set_linewidth(0.6); axb.spines[s].set_color("0.45")
axb.legend(loc="upper left", frameon=False, fontsize=6.0, handlelength=1.8,
           labelspacing=0.28, handletextpad=0.5)
axb.set_title("top-500 overlap (split 1)", fontsize=7.6, fontweight="bold", pad=6, loc="left")

fig.text(0.075, 0.955, "The compound ranking is reproducible across independent patient halves",
         fontsize=9, fontweight="bold", ha="left")
fig.text(0.075, 0.905, f"Patients split at random, each half scored against LINCS separately; "
         f"Spearman of the two rankings.  Median across splits = {R.rho.median():.3f}.",
         fontsize=6.4, color="0.35", ha="left")
fig.savefig(f"{OUT}/P17_split_half_stability.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P17_split_half_stability.png", dpi=400, bbox_inches="tight")
print("✓ P17_split_half_stability")
