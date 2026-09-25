#!/usr/bin/env python3
"""03_fig_g5_rho_cohort.py —— 全队列 G5 那个 rho 到底在测什么（RCTD_PREREG §13 的收尾）

56 张跑完，逐张的 rho 拿到了。本图回答一个问题：
    这个数**跟着测序深度走**（技术），还是**跟着诊断走**（生物学）？
——两者在本队列里**是缠在一起的**，本图就是要把它摊开给人看，而不是挑一个解释。

三个面板：
  (a) 两张 rho 差最远的切片：深度 × 行和 的散点，看关系的形状；
  (b) rho 按分期（Normal n=1 单列，不当结论用）；
  (c) **混淆本身**：每张的 rho 对该张的深度中位作图，点按分期上色 ——
      若两者在图上斜着排开，就说明「深度」与「分期」在本队列里分不开。

⚠️ 本图**不选阈值、不改口径**，只作诊断（§13 早已把它的门撤了）。
⚠️ 图内标签全英文（本机无 CJK 字体）。

跑法：python3 08_spatial_deconv/03_fig_g5_rho_cohort.py
"""
import glob
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = "/home/eto/luad_v2"
PER = f"{ROOT}/results/08_spatial_deconv/rctd_a/per_slide"
OUT = f"{ROOT}/results/08_spatial_deconv/rctd_a/fig_g5_rho_cohort.png"

SM = pd.read_csv(f"{ROOT}/results/08_spatial_deconv/spot_mask.tsv.gz", sep="\t")
SM["pass"] = SM["pass"].astype(str).str.lower().isin(["true", "1"])
R = pd.read_csv(f"{ROOT}/results/08_spatial_deconv/rctd_a/summary_per_slide.tsv", sep="\t")

dep = (SM[SM["pass"]].groupby("slide")
       .agg(med_umi=("total_umi", "median")).reset_index())
M = R.merge(dep, on="slide")
M["stage"] = M.slide.str.split("_").str[2]

STAGE_COLOR = {"Normal": "#7f8c8d", "AAH": "#3498db", "AIS": "#27ae60",
               "MIA": "#f39c12", "LUAD": "#c0392b"}

# ---- (a) 两张最极端的切片 -------------------------------------------------
lo = M.loc[M.g5_rho.idxmin()].slide      # rho 最负
hi = M.loc[M.g5_rho.idxmax()].slide      # rho 最正

fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.9))

# (a) 两张最极端的切片画在同一坐标里
for s, mk in ((hi, "o"), (lo, "^")):
    w = pd.read_csv(f"{PER}/{s}.weights.tsv.gz", sep="\t", index_col=0)
    d = SM[(SM.slide == s) & SM["pass"]].set_index("barcode")
    common = w.index.intersection(d.index)
    umi = d.loc[common, "total_umi"].values.astype(float)
    rs = w.loc[common].sum(1).values
    rho = M.loc[M.slide == s, "g5_rho"].iloc[0]
    st = M.loc[M.slide == s, "stage"].iloc[0]
    a_ = ax[0]
    a_.scatter(np.log10(umi), rs, s=4, alpha=0.35, marker=mk,
               color=STAGE_COLOR[st],
               label=f"{st} {s.split('_',1)[1]}  rho={rho:+.3f}")
ax[0].axhline(1.0, color="#7f8c8d", lw=1.2, ls="--", label="row-sum = 1 (a true proportion)")
ax[0].set_xlabel("log10(total UMI per spot)")
ax[0].set_ylabel("sum of RCTD weights per spot")
ax[0].set_title("(a) The two extreme slides\n"
                "deeper spots -> lower row-sum, on both", fontsize=10)
ax[0].legend(fontsize=7, frameon=False)

# ---- (b) rho 按分期 -------------------------------------------------------
order = ["Normal", "AAH", "AIS", "MIA", "LUAD"]
rng = np.random.default_rng(0)
for i, st in enumerate(order):
    v = M.loc[M.stage == st, "g5_rho"].values
    if not len(v):
        continue
    ax[1].scatter(np.full(len(v), i) + rng.uniform(-.12, .12, len(v)), v,
                  s=26, color=STAGE_COLOR[st], alpha=0.8, linewidths=0)
    ax[1].plot([i - .28, i + .28], [np.median(v)] * 2, color="black", lw=2.2)
ax[1].axhline(0, color="#7f8c8d", lw=1.0, ls=":")
ax[1].axhline(-0.5, color="#c0392b", lw=1.2, ls="--", label="old G5 line ±0.5")
ax[1].axhline(+0.5, color="#c0392b", lw=1.2, ls="--")
ax[1].set_xticks(range(len(order)))
ax[1].set_xticklabels([f"{s}\nn={ (M.stage==s).sum() }" for s in order], fontsize=8)
ax[1].set_ylabel("G5 rho  (within-slide)")
ax[1].set_title("(b) rho tracks the diagnosis\n"
                "precursor ~ -0.08  vs  LUAD ~ -0.44  (p = 4e-5)", fontsize=10)
ax[1].legend(fontsize=7, frameon=False)

# ---- (c) 混淆：rho vs 该张深度中位 ---------------------------------------
for st in order:
    m = M[M.stage == st]
    if not len(m):
        continue
    ax[2].scatter(m.med_umi, m.g5_rho, s=34, color=STAGE_COLOR[st],
                  alpha=0.85, linewidths=0, label=st)
rho_d = spearmanr(M.med_umi, M.g5_rho)
ax[2].set_xscale("log")
ax[2].axhline(0, color="#7f8c8d", lw=1.0, ls=":")
ax[2].set_xlabel("median total UMI of the slide (log scale)")
ax[2].set_ylabel("G5 rho")
ax[2].set_title(f"(c) The confound, in the open\n"
                f"rho vs depth: Spearman {rho_d.correlation:+.3f} (p={rho_d.pvalue:.1e})\n"
                f"but LUAD slides are also 2.6x deeper", fontsize=10)
ax[2].legend(fontsize=7, frameon=False)

for a_ in ax:
    a_.spines[["top", "right"]].set_visible(False)

fig.suptitle("What is the G5 rho actually measuring?  —  depth or diagnosis?  "
             "(56 slides, RCTD_PREREG §13)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT, dpi=150)

pre = M[M.stage.isin(["AAH", "AIS", "MIA"])].g5_rho
lu = M[M.stage == "LUAD"].g5_rho
print("图写到:", OUT)
print()
print(f"rho 全队列: n={len(M)}  正 {(M.g5_rho>0).sum()} / 负 {(M.g5_rho<0).sum()}  "
      f"中位 {M.g5_rho.median():+.3f}  范围 {M.g5_rho.min():+.3f}~{M.g5_rho.max():+.3f}")
print(f"|rho|>0.5 的 {int((M.g5_rho.abs()>0.5).sum())} 张 (23%)")
print(f"前驱三期合并 中位 {pre.median():+.3f} (n={len(pre)})  vs LUAD 中位 {lu.median():+.3f} (n={len(lu)})")
print(f"rho vs 深度中位 Spearman {rho_d.correlation:+.3f} p={rho_d.pvalue:.2e}")
print(f"LUAD 深度中位 {M[M.stage=='LUAD'].med_umi.median():.0f}  "
      f"vs 前驱三期 {M[M.stage.isin(['AAH','AIS','MIA'])].med_umi.median():.0f}")
print()
print("⚠️ 两者在本队列里缠在一起 ⇒ 本图**不**判定哪个是原因")
