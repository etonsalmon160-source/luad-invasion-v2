#!/usr/bin/env python
"""预注册 20 的 §9 出图。英文标签——本机无 CJK 字体。

读 results/09_trajectory/paper_axis/ 的产物，出三张：
  1. 簇 × MP 热图（行按 argmax 分组、列按 panel 顺序），标出 KAC 侧命中
  2. KAC 签名逐簇得分条形图（降序，KAC 侧红色）
  3. MP6 / MP9 与 KAC 签名的簇级散点（D2 的可视化）
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

D = "/home/eto/luad_v2/results/09_trajectory/paper_axis"
FIG = f"{D}/figures"

cm = pd.read_csv(f"{D}/cluster_scores_papermp.csv")
score_cols = [c for c in cm.columns if c.startswith("s_")]
mp_cols = [c for c in score_cols if c != "s_KACsig"]

cm = cm.sort_values(["is_kac_side", "argmax_papermp", "win_margin"],
                    ascending=[False, True, False]).reset_index(drop=True)

# ---- 图1：热图 ----
Z = cm[score_cols].to_numpy(float)
Zz = (Z - Z.mean(0)) / (Z.std(0) + 1e-12)          # 逐列 z 便于同图比较
fig, ax = plt.subplots(figsize=(12, 10))
im = ax.imshow(Zz, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax.set_xticks(range(len(score_cols)))
ax.set_xticklabels([c[2:].replace("_", " ") for c in score_cols],
                   rotation=40, ha="right", fontsize=9)
ax.set_yticks(range(len(cm)))
ax.set_yticklabels([f"c{int(r.cluster)}  n={int(r.n_cells):,}" for r in cm.itertuples()],
                   fontsize=8)
for lab, r in zip(ax.get_yticklabels(), cm.itertuples()):
    if r.is_kac_side:
        lab.set_color("crimson")
        lab.set_fontweight("bold")
ax.set_xlim(-0.5, len(score_cols) + 2.0)
for i, r in cm.iterrows():
    ax.text(len(score_cols) - 0.35, i, r.argmax_papermp, fontsize=7.5,
            va="center", ha="left", clip_on=False,
            color=("crimson" if r.is_kac_side else "#333333"))
ax.set_title("Paper meta-program scores by cluster (column z-score)\n"
             f"n={int(cm.n_cells.sum()):,} cells, {len(cm)} clusters — "
             f'D1 hits (argmax in KAC side): {int(cm.is_kac_side.sum())} '
             "(red cluster labels = D1 hit)", fontsize=10)
cb = fig.colorbar(im, ax=ax, shrink=0.45, pad=0.04)
cb.set_label("z-score within column", fontsize=9)
plt.tight_layout()
plt.savefig(f"{FIG}/paper_mp_heatmap.png", dpi=140)
print("写好 paper_mp_heatmap.png")

# ---- 图2：KAC 签名逐簇 ----
s = cm.sort_values("s_KACsig", ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(11, 4.6))
col = ["crimson" if k else "#8fbcd4" for k in s.is_kac_side]
ax.bar(range(len(s)), s.s_KACsig, color=col)
ax.axhline(0, color="k", lw=.8)
for i, r in s.iterrows():
    ax.text(i, r.s_KACsig + (0.01 if r.s_KACsig >= 0 else -0.03),
            f"c{int(r.cluster)}", rotation=90, fontsize=7, ha="center")
ax.set_xticks([])
ax.set_ylabel("Table S3 KAC signature score")
ax.set_title("KAC signature score per cluster (red = D1 KAC-side hit)\n"
             "this is the independent 2nd source vs MP6 (only 8 genes overlap)", fontsize=10)
ax.grid(alpha=.25, axis="y")
plt.tight_layout()
plt.savefig(f"{FIG}/kac_signature_by_cluster.png", dpi=140)
print("写好 kac_signature_by_cluster.png")

# ---- 图3：D2 散点 ----
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
for ax, mp in zip(axes, ["s_Tumor_cell_KAC", "s_KAC_inflammatory"]):
    x, y = cm[mp].to_numpy(), cm["s_KACsig"].to_numpy()
    rho, p = stats.spearmanr(x, y)
    ax.scatter(x, y, c=["crimson" if k else "#8fbcd4" for k in cm.is_kac_side], s=52)
    for i, r in cm.iterrows():
        ax.annotate(f"c{int(r.cluster)}", (cm[mp][i], cm["s_KACsig"][i]),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(mp[2:].replace("_", " ") + " score")
    ax.set_ylabel("KAC signature score")
    ax.set_title(f"D2: {mp[2:].replace('_', ' ')}\nSpearman rho={rho:+.3f}  p={p:.3g}",
                 fontsize=10)
    ax.grid(alpha=.25)
plt.tight_layout()
plt.savefig(f"{FIG}/d2_mp_vs_kacsignature.png", dpi=140)
print("写好 d2_mp_vs_kacsignature.png")
print("KAC 侧命中:", cm.loc[cm.is_kac_side, "cluster"].astype(int).tolist())
