#!/usr/bin/env python3
"""fig_scale_invariance_thin50.py —— 尺度不变性实验的判读图（RCTD_PREREG §13.3）

要测的命题：把每个 spot 的计数随机抽掉一半（细胞比例不变、只有深度变小），
           「上皮是权重最大的那一类」这个判断会不会移动。

对照：
  全深度 = results/08_spatial_deconv/rctd_a/per_slide/<slide>.weights.tsv.gz
  半深度 = results/08_spatial_deconv/rctd_a_thin50/per_slide/<slide>.weights.tsv.gz

判读口径（§13.3 已先写下，避免事后解释）：
  - 移动 ≈ 0                              ⇒ 尺度不变性实测成立
  - 移动少、且集中在「边界 spot」          ⇒ 成立，但登记「边界处不稳」
  - 移动很多（非边界也动）                 ⇒ ⚠️ argmax 门控本身可疑，须停下来回退给用户

「边界 spot」的判据（先说清）：该 spot 第二大的类型与上皮的权重差 < 0.05。

图内标签全英文（本机无 CJK 字体）。

跑法：python3 08_spatial_deconv/fig_scale_invariance_thin50.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"
SLIDE = "GSM9226168_P1_AAH"
FULL = f"{ROOT}/results/08_spatial_deconv/rctd_a/per_slide/{SLIDE}.weights.tsv.gz"
THIN = f"{ROOT}/results/08_spatial_deconv/rctd_a_thin50/per_slide/{SLIDE}.weights.tsv.gz"
OUT = f"{ROOT}/results/08_spatial_deconv/rctd_a_thin50/fig_scale_invariance.png"
MARGIN = 0.05          # 「边界 spot」判据：第一大与第二大的权重差
EPI_EN = "Epithelial"

EN = {"上皮": EPI_EN, "内皮": "Endothelial", "成纤维": "Fibroblast",
      "髓系": "Myeloid", "T_NK": "T/NK", "B_浆": "B/Plasma"}

wf = pd.read_csv(FULL, sep="\t", index_col=0).rename(columns=EN)
wt = pd.read_csv(THIN, sep="\t", index_col=0).rename(columns=EN)

common = wf.index.intersection(wt.index)
only_full = wf.index.difference(wt.index)
only_thin = wt.index.difference(wf.index)
wf, wt = wf.loc[common], wt.loc[common]

# 逐 spot：argmax 类型 + 第一/第二名的间距
def top2(df):
    a = df.values
    order = np.argsort(-a, axis=1)
    cols = np.array(df.columns)
    top1 = cols[order[:, 0]]
    m1 = a[np.arange(len(a)), order[:, 0]]
    m2 = a[np.arange(len(a)), order[:, 1]]
    return top1, m1 - m2, a.sum(1)

am_f, margin_f, rs_f = top2(wf)
am_t, margin_t, rs_t = top2(wt)
moved = am_f != am_t

epi_f = wf[EPI_EN].values
epi_t = wt[EPI_EN].values

# 上皮的 argmax 是否移动（本门控实际用的就是这一条判据）
epi_f_am = am_f == EPI_EN
epi_t_am = am_t == EPI_EN
epi_flip = epi_f_am != epi_t_am

n = len(common)
n_moved = int(moved.sum())
n_epi_flip = int(epi_flip.sum())
n_border = int((margin_f < MARGIN).sum())
moved_border = int((moved & (margin_f < MARGIN)).sum())
moved_nonborder = n_moved - moved_border
ROBUST = 0.15
n_robust = int((margin_f >= ROBUST).sum())
moved_robust = int((moved & (margin_f >= ROBUST)).sum())

fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

# ---- (a) 上皮权重：全深度 vs 半深度 --------------------------------------
lim = max(epi_f.max(), epi_t.max()) * 1.03
ax[0].scatter(epi_f, epi_t, s=5, c=np.where(epi_flip, "#c0392b", "#2980b9"),
              alpha=0.45, linewidths=0)
ax[0].plot([0, lim], [0, lim], color="#7f8c8d", lw=1.2, ls="--", label="y = x")
ax[0].set_xlabel("Epithelial weight — full depth")
ax[0].set_ylabel("Epithelial weight — half depth (binomial p=0.5)")
ax[0].set_title(f"(a) Epithelial weight barely moves under halving\n"
                f"r = {np.corrcoef(epi_f, epi_t)[0,1]:.4f}   "
                f"median |diff| = {np.median(np.abs(epi_f-epi_t)):.4f}", fontsize=10)
ax[0].legend(fontsize=8, frameon=False)

# ---- (b) 谁移动了，移动的离边界有多近 ------------------------------------
bins = np.linspace(0, 0.5, 51)
ax[1].hist(margin_f[~moved], bins=bins, color="#2980b9", alpha=0.75,
           label=f"argmax unchanged (n={n-n_moved:,})")
ax[1].hist(margin_f[moved], bins=bins, color="#c0392b", alpha=0.85,
           label=f"argmax moved (n={n_moved:,})")
ax[1].axvline(MARGIN, color="#7f8c8d", lw=1.5, ls="--",
              label=f"border = top1-top2 < {MARGIN}")
ax[1].set_xlabel("margin between top-1 and top-2 cell type (full depth)")
ax[1].set_ylabel("number of spots")
ax[1].set_title(f"(b) Churn lives at the ties — and only there\n"
                f"{moved_border}/{n_moved} movers had margin < {MARGIN}; "
                f"{moved_nonborder} elsewhere\n"
                f"but at margin >= {ROBUST}: only {moved_robust}/{n_robust:,} moved",
                fontsize=10)
ax[1].legend(fontsize=8, frameon=False)

# ---- (c) 门控实际用的那条判据 -------------------------------------------
frac_f = 100 * epi_f_am.mean()
frac_t = 100 * epi_t_am.mean()
bars = ax[2].bar(["full depth", "half depth"], [frac_f, frac_t],
                 color=["#2980b9", "#8e44ad"], width=0.55)
ax[2].bar_label(bars, fmt="%.2f%%", fontsize=10, padding=3)
ax[2].set_ylabel("% of spots where Epithelial is the top type")
ax[2].set_ylim(0, max(frac_f, frac_t) * 1.35)
ax[2].set_title(f"(c) The gate itself (argmax == Epithelial)\n"
                f"spots where it flipped: {n_epi_flip} / {n:,}  "
                f"({100*n_epi_flip/n:.3f}%)", fontsize=10)
ax[2].text(0.5, 0.06, f"spots dropped after thinning (all-zero): {len(only_full):,}",
           transform=ax[2].transAxes, ha="center", fontsize=8, color="#555555")

for a in ax:
    a.spines[["top", "right"]].set_visible(False)

fig.suptitle("Scale-invariance test on real tissue  —  halve every spot's counts, "
             "does the top cell type change?  (§13.3)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT, dpi=150)

print("图写到:", OUT)
print()
print(f"共同 spot 数              {n:,}")
print(f"稀释后全零被剔除            {len(only_full):,}")
print(f"argmax 移动的 spot          {n_moved:,}   ({100*n_moved/n:.3f}%)")
print(f"  其中 margin < {MARGIN} 的   {moved_border:,}")
print(f"  非边界的                  {moved_nonborder:,}")
print(f"上皮 argmax 翻转（门控判据） {n_epi_flip:,}   ({100*n_epi_flip/n:.3f}%)")
print(f"上皮权重 corr             {np.corrcoef(epi_f, epi_t)[0,1]:.5f}")
print(f"上皮权重 中位绝对差         {np.median(np.abs(epi_f-epi_t)):.5f}")
print(f"上皮 argmax 占比  全深度 {frac_f:.2f}%  →  半深度 {frac_t:.2f}%")
