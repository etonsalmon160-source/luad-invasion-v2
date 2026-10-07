#!/usr/bin/env python3
# 44_fig_domain_cnv.py —— 生态位域的空间 CNV 分值（同切片内配对，深度自动抵消）
#   口径（NICHE_PREREG §18.1–18.2）：每个域内 spot 的 cf 减掉**同一切片**其余 spot 的 cf，
#   逐切片算一次 Δ，再跨切片汇总。cf 是连续量，**不是恶性判定**。
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
D = pd.read_csv(f"{ROOT}/results/10_niche/kstar_diag/d8_spot_archetype_cnv.tsv.gz", sep="\t")
D = D[~D.grp.isna()] if "grp" in D.columns else D   # 全部 spot（含 reference），与 §18.2 口径一致
D = D.dropna(subset=["cf", "archetype"])
print(f"spot {len(D):,} / 切片 {D.slide.nunique()}")

# 逐切片、逐域：域内均值 − 该切片其余 spot 均值
rows = []
for sl, g in D.groupby("slide"):
    if len(g) < 200:
        continue
    tot = g.cf.median()                       # 中位口径，与 §18.2 一致
    for a, h in g.groupby("archetype"):
        if len(h) < 1:
            continue
        rows.append(dict(slide=sl, arch=int(a), delta=h.cf.median() - tot, n=len(h)))
P = pd.DataFrame(rows)
P.to_csv("/home/eto/luad_v2/results/10_niche/domain_cnv_delta_perslide.tsv", sep="\t", index=False)

S = P.groupby("arch").agg(n_slide=("slide", "size"), n_pos=("delta", lambda x: int((x > 0).sum())),
                          d_med=("delta", "median"), d_mean=("delta", "mean"),
                          n_spot=("n", "sum")).reset_index()
S["frac_pos"] = S.n_pos / S.n_slide
S["p_sign"] = [stats.binomtest(int(r.n_pos), int(r.n_slide), 0.5).pvalue for r in S.itertuples()]
print(S.round(4).to_string(index=False))

DOMC = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
DOMLAB = ["D1\nAirway", "D2\niCAF", "D3\nECM/interstitial", "D4\nAlveolar-cap.",
          "D5\nAT2", "D6\nVascular", "D7\nLymphoid"]

fig, ax = plt.subplots(figsize=(6.4, 3.5))
xs = np.arange(7)
for i, a in enumerate(S.arch):
    v = P[P.arch == a].delta.values * 100
    ax.scatter(np.full(len(v), i) + np.random.default_rng(7).uniform(-0.17, 0.17, len(v)),
               v, s=7, c="0.35", alpha=0.55, linewidths=0, zorder=2)
    bar = S.loc[i, "d_med"] * 100          # 中位（§18.2 报的就是中位）
    ax.bar(i, bar, width=0.30, color=DOMC[i], zorder=3)
    p = S.loc[i, "p_sign"]
    star = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "n.s."))
    ax.text(i, max(bar, 0) + 0.35, star, ha="center", va="bottom", fontsize=7.4,
            fontweight="bold", color="#C0392B" if p < 0.05 else "0.45")
    ax.text(i, -5.6, f"{int(S.loc[i,'n_pos'])}/{int(S.loc[i,'n_slide'])}",
            ha="center", va="top", fontsize=6.0, color="0.45")

ax.axhline(0, color="0.55", lw=0.7)
ax.set_xticks(xs); ax.set_xticklabels(DOMLAB, fontsize=6.4, linespacing=1.3)
ax.set_ylabel("CNV score of the domain\n− rest of the same section   (×100)", fontsize=7)
ax.set_ylim(-6.8, 9.0); ax.set_xlim(-0.6, 6.6)
ax.tick_params(axis="y", labelsize=6.5, length=2.5); ax.tick_params(axis="x", length=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.6); ax.spines[s].set_color("0.45")
ax.set_title("Spatial CNV score per niche domain", fontsize=8.8, fontweight="bold", pad=16,
             loc="left")
fig.savefig(f"{OUT}/P14_domain_spatial_cnv.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P14_domain_spatial_cnv.png", dpi=400, bbox_inches="tight")
print("✓ P14_domain_spatial_cnv")
