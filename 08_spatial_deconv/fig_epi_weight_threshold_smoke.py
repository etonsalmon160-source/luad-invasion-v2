#!/usr/bin/env python3
"""fig_epi_weight_threshold_smoke.py —— 门控阈值在真实数据上做了什么

输入：冒烟测试单张切片（GSM9226168_P1_AAH）的逐 spot 权重，**直接读项目目录里的正式产物**
     results/08_spatial_deconv/rctd_a/per_slide/<slide>.weights.tsv.gz
     （2026-09-25 更正：早期版本读 /tmp 的一份 9,932 spot 的中间稿，那个 spot 集合是错的
      —— 它是正确答案 9,870 的超集，多出的 62 个上皮权重中位 0.001、最大 0.043，全在分布底部，
      故当时的数字与结论不受影响；但输入必须指向正式产物，见 RCTD_PREREG.md §12.5 第 5 项。）
目的：这是**参数问题**，按项目规矩给图，不给参数对照表。
     要回答的是：阈值在这张切片上切掉多少、切下来的是哪些 spot。

图内标签全英文（本机无 CJK 字体）。

跑法：python3 08_spatial_deconv/fig_epi_weight_threshold_smoke.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SLIDE = "GSM9226168_P1_AAH"
SRC = ("/home/eto/luad_v2/results/08_spatial_deconv/rctd_a/per_slide/"
       f"{SLIDE}.weights.tsv.gz")
OUT = ("/home/eto/luad_v2/results/08_spatial_deconv/rctd_a/"
       "fig_epi_weight_threshold_smoke.png")

# 列名是中文（注释表原名）⇒ 图内必须换英文（本机无 CJK 字体）
EN = {"上皮": "Epithelial", "内皮": "Endothelial", "成纤维": "Fibroblast",
      "髓系": "Myeloid", "T_NK": "T/NK", "B_浆": "B/Plasma"}

w = pd.read_csv(SRC, sep="\t", index_col=0)
w = w.rename(columns=EN)
epi = w["Epithelial"].values
rowsum = w.sum(1).values
n = len(epi)

SIGNED = 0.5                      # 已签字口径
pcts = {p: np.percentile(epi, p) for p in (60, 70, 80)}

fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.0))

# ---- 左：上皮权重分布 + 阈值线 -------------------------------------------
ax[0].hist(epi, bins=60, range=(0, max(1.0, epi.max())), color="#8fb8de",
           edgecolor="white", linewidth=0.4)
ax[0].axvline(SIGNED, color="#c0392b", lw=2.0,
              label=f"signed gate  >{SIGNED}   keeps {100*(epi>SIGNED).mean():.1f}%")
ax[0].axvline(np.median(epi), color="#2c3e50", lw=1.6, ls="--",
              label=f"median = {np.median(epi):.3f}")
MAJORITY = 0.5 * np.median(rowsum)   # 真·「多数」线 = 行和中位的一半
ax[0].axvline(MAJORITY, color="#8e44ad", lw=1.6, ls="-.",
              label=f"true majority line = {MAJORITY:.3f}\n(0.5 x median row-sum)"
                    f" keeps {100*(epi>MAJORITY).mean():.1f}%")
for p in (60, 70, 80):
    ax[0].axvline(pcts[p], color="#7f8c8d", lw=1.1, ls=":",
                  label=f"P{p} = {pcts[p]:.3f}")
ax[0].set_xlabel("RCTD raw Epithelial weight per spot")
ax[0].set_ylabel("number of spots")
ax[0].set_title(f"Epithelial weight distribution\n{SRC.split('/')[-1]}  "
                f"({n:,} spots, caliber a)", fontsize=10)
ax[0].legend(fontsize=8, frameon=False)
ax[0].spines[["top", "right"]].set_visible(False)

# ---- 右：保留比例 vs 阈值 ------------------------------------------------
th = np.linspace(0, max(0.95, epi.max()), 400)
keep = np.array([(epi > t).mean() * 100 for t in th])
ax[1].plot(th, keep, color="#2980b9", lw=2.0, label="Epithelial weight > threshold")
argmax_frac = (w.idxmax(1) == "Epithelial").mean() * 100
ax[1].axhline(argmax_frac, color="#27ae60", lw=1.6, ls="-.",
              label=f"argmax = Epithelial  ({argmax_frac:.1f}%)")
ax[1].axvline(SIGNED, color="#c0392b", lw=2.0)
ax[1].plot([SIGNED], [(epi > SIGNED).mean() * 100], "o", color="#c0392b", ms=7)
ax[1].annotate(f" >{SIGNED}: {(epi>SIGNED).mean()*100:.1f}%",
               xy=(SIGNED, (epi > SIGNED).mean() * 100),
               xytext=(SIGNED + 0.08, (epi > SIGNED).mean() * 100 + 12),
               color="#c0392b", fontsize=9,
               arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.2))
for p in (60, 70, 80):
    ax[1].axvline(pcts[p], color="#7f8c8d", lw=1.1, ls=":")
ax[1].set_xlabel("threshold on RCTD raw Epithelial weight")
ax[1].set_ylabel("% of spots retained")
ax[1].set_title("What the threshold actually cuts\n"
                "(dotted = P60/P70/P80 of the weight itself)", fontsize=10)
ax[1].legend(fontsize=8, frameon=False)
ax[1].spines[["top", "right"]].set_visible(False)

fig.suptitle("Spatial gating threshold on raw Epithelial weight  —  "
             "smoke test, single slide (P1 AAH)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT, dpi=150)
print("图写到:", OUT)
print()
print("行和分位 (0/25/50/75/100):", np.round(np.percentile(rowsum, [0, 25, 50, 75, 100]), 3))
print("上皮权重分位 (10/25/50/75/90/95/99):",
      np.round(np.percentile(epi, [10, 25, 50, 75, 90, 95, 99]), 3))
print("上皮权重最大值:", round(epi.max(), 3))
print("“多数线” = 0.5 × 行和中位 =", round(0.5 * np.median(rowsum), 3))
for t in (0.4, 0.45, 0.5, 0.6, 0.7):
    print(f"  > {t}: {int((epi>t).sum()):>5} spot ({100*(epi>t).mean():>5.1f}%)")
print("argmax=上皮:", int((w.idxmax(1) == 'Epithelial').sum()),
      f"({argmax_frac:.1f}%)")
