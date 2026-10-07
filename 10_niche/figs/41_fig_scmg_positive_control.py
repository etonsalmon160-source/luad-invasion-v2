#!/usr/bin/env python3
# 41_fig_scmg_positive_control.py —— SCMG 阳性对照（照原文 causal_gene_prediction 的图式）
#   数据：tr_singlecell/scmg_PC_gastrulation.tsv（8,450 个扰动 × 15 个扰动数据集）
#   图式：x = gene_shift_z（被扰动基因自身的表达位移）
#         y = perturbation match score（扰动与「上胚层→中胚层」状态转移的匹配分）
#   判读：causal gene 落在「两轴同号且绝对值大」的两个象限
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
d = pd.read_csv(f"{ROOT}/results/10_niche/tr_singlecell/scmg_PC_gastrulation.tsv", sep="\t")
d = d.sort_values("causal_score", ascending=False).reset_index(drop=True)
d["rank"] = d.index + 1
N = len(d)

POS, NEG = "#2E6FB7", "#D1783C"
KNOWN = ["TBXT", "MSGN1", "MIXL1", "TBX6", "SNAI1", "EOMES", "EVX1", "NANOG", "POU5F1", "L1TD1", "PRDM14"]
lab = d[d.perturbed_gene_name.isin(KNOWN)].copy()

fig, ax = plt.subplots(figsize=(5.4, 5.0))
ax.axhline(0, color="0.75", lw=0.6, zorder=1)
ax.axvline(0, color="0.75", lw=0.6, zorder=1)
for sgn, col in ((1, POS), (-1, NEG)):
    m = d.perturbation_sign == sgn
    ax.scatter(d.gene_shift_z[m], d.pert_match_score[m], s=2.2, c=col, alpha=0.30,
               linewidths=0, rasterized=True, zorder=2)

# 已知基因：标注 + 秩位
off = {"TBXT": (0.26, 0.014), "MSGN1": (-0.98, 0.014), "MIXL1": (0.14, 0.006),
       "TBX6": (0.16, -0.030), "SNAI1": (0.12, 0.030), "EOMES": (-0.80, 0.024),
       "EVX1": (-0.74, -0.034), "NANOG": (-0.60, -0.034), "POU5F1": (0.18, -0.010),
       "L1TD1": (0.20, 0.012), "PRDM14": (0.18, -0.038)}
for _, r in lab.iterrows():
    col = POS if r.perturbation_sign > 0 else NEG
    dx, dy = off.get(r.perturbed_gene_name, (0.18, 0.016))
    ax.scatter([r.gene_shift_z], [r.pert_match_score], s=22, facecolor="white",
               edgecolor=col, linewidths=1.1, zorder=4)
    ax.annotate(f"{r.perturbed_gene_name}  #{int(r['rank'])}",
                (r.gene_shift_z, r.pert_match_score), xytext=(r.gene_shift_z + dx,
                r.pert_match_score + dy), fontsize=6.2, color=col, fontweight="bold",
                zorder=5, va="center",
                arrowprops=dict(arrowstyle="-", color=col, lw=0.6, shrinkA=0, shrinkB=2))

ax.set_xlabel("Gene expression shift  (z-score)", fontsize=7.5)
ax.set_ylabel("Perturbation match score", fontsize=7.5)
ax.set_xlim(-5.4, 6.4); ax.set_ylim(-0.315, 0.335)
ax.tick_params(labelsize=6.5, length=2.5)
for s in ax.spines.values():
    s.set_linewidth(0.6); s.set_color("0.45")
ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4.5, mfc=POS, mec="none",
                          label=f"knockdown / activation, +  (n={int((d.perturbation_sign>0).sum()):,})"),
                   Line2D([], [], marker="o", ls="", ms=4.5, mfc=NEG, mec="none",
                          label=f"knockdown, −  (n={int((d.perturbation_sign<0).sum()):,})")],
          loc="lower right", frameon=False, fontsize=6.2, handletextpad=0.3, labelspacing=0.28)
ax.set_title("Positive control: epiblast → nascent mesoderm\n"
             f"the {N:,} perturbations of the SCMG library, ranked by causal score",
             fontsize=8.6, fontweight="bold", linespacing=1.5, pad=7)
fig.savefig(f"{OUT}/P11_scmg_positive_control.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P11_scmg_positive_control.png", dpi=400, bbox_inches="tight")
print("✓ P11_scmg_positive_control")
print(lab[["rank", "perturbed_gene_name", "gene_shift_z", "pert_match_score",
           "causal_score"]].to_string(index=False))
