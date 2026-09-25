#!/usr/bin/env python3
"""02_margin_cohort.py —— 全队列「逐 spot margin 分布」汇总（RCTD_PREREG §13.5.4 的下一步）

背景：尺度不变性实验（§13.3/§13.5，单张 GSM9226168_P1_AAH）测出——
      把计数减半重跑，`argmax == 上皮` 这条门控判据翻转 8.30%，两组上皮集合 Jaccard 0.840；
      但抖动**全部**集中在「第一名与第二名咬得很紧」的 spot 上，
      `margin >= 0.15` 的 2,516 个 spot **零移动**。
      用户裁定（§13.5.4）：**口径先不动，等 56 张跑完，先出全队列的 margin 分布再定。**

本脚本就出这个分布。**不选任何新阈值、不改口径**——只把「每张切片的 margin 长什么样、
上皮占比在「全体」与「margin 够大」两种取法下差多少」摊开。

margin 的定义：该 spot 权重最大的那一类 减去 第二大的那一类（都取原始权重，不做归一化）。

⚠️ 输入是**逐张跑出来的正式产物**，跑完几张算几张（可反复重跑，后跑的会补进来）。
⚠️ 图内标签全英文（本机无 CJK 字体）。

跑法：python3 08_spatial_deconv/02_margin_cohort.py
"""
import glob
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"
PER = f"{ROOT}/results/08_spatial_deconv/rctd_a/per_slide"
OUT_TSV = f"{ROOT}/results/08_spatial_deconv/rctd_a/margin_cohort.tsv"
OUT_PNG = f"{ROOT}/results/08_spatial_deconv/rctd_a/fig_margin_cohort.png"

EPI = "上皮"
STABLE = 0.15      # 只作**参照线**（§13.5 实测该处单张零移动），**不是**新口径
TIGHT = 0.05       # 同上，作参照

files = sorted(glob.glob(f"{PER}/*.weights.tsv.gz"))
if not files:
    raise SystemExit("还没有任何逐张权重产物")

rows, curves = [], {}
for f in files:
    slide = os.path.basename(f).replace(".weights.tsv.gz", "")
    w = pd.read_csv(f, sep="\t", index_col=0)
    a = w.values
    o = np.argsort(-a, axis=1)
    cols = np.array(w.columns)
    top1 = cols[o[:, 0]]
    margin = a[np.arange(len(a)), o[:, 0]] - a[np.arange(len(a)), o[:, 1]]
    n = len(w)
    epi = top1 == EPI
    stab = margin >= STABLE

    rows.append(dict(
        slide=slide, n_spot=n,
        epi_frac=epi.mean(),
        epi_frac_stable=epi[stab].mean() if stab.sum() else np.nan,
        n_stable=int(stab.sum()), stable_frac=stab.mean(),
        margin_p10=np.percentile(margin, 10), margin_p50=np.percentile(margin, 50),
        margin_p90=np.percentile(margin, 90),
        tight_frac=(margin < TIGHT).mean(),
        top1_endothelial_frac=(top1 == "内皮").mean(),
    ))
    # 每张的 margin 累积曲线（用于叠图看形状是否一致）
    xs = np.linspace(0, 0.4, 200)
    curves[slide] = np.searchsorted(np.sort(margin), xs) / n

df = pd.DataFrame(rows).sort_values("slide").reset_index(drop=True)
df.to_csv(OUT_TSV, sep="\t", index=False)

fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.8))

# ---- (a) 逐张 margin 累积曲线 -------------------------------------------
xs = np.linspace(0, 0.4, 200)
for s, y in curves.items():
    ax[0].plot(xs, y, lw=0.9, alpha=0.55, color="#2980b9")
med = np.median(np.vstack(list(curves.values())), axis=0)
ax[0].plot(xs, med, lw=2.6, color="#c0392b", label="cohort median")
ax[0].axvline(STABLE, color="#7f8c8d", ls="--", lw=1.4,
              label=f"reference line {STABLE} (not a caliber)")
ax[0].set_xlabel("margin between top-1 and top-2 cell type")
ax[0].set_ylabel("cumulative fraction of spots")
ax[0].set_title(f"(a) Margin CDF per slide ({len(curves)} slides)\n"
                f"how tight is the winning type?", fontsize=10)
ax[0].legend(fontsize=8, frameon=False)

# ---- (b) 逐张上皮占比：全体 vs 只在 margin 够大的 spot 上 -----------------
x = np.arange(len(df))
ax[1].bar(x - 0.2, 100 * df.epi_frac, width=0.4, color="#2980b9",
          label="all spots (argmax == Epithelial)")
ax[1].bar(x + 0.2, 100 * df.epi_frac_stable, width=0.4, color="#8e44ad",
          label=f"only spots with margin >= {STABLE}")
ax[1].set_xticks(x)
ax[1].set_xticklabels([s.split("_", 1)[1] if "_" in s else s for s in df.slide],
                      rotation=90, fontsize=5)
ax[1].set_xlabel("slide")
ax[1].set_ylabel("% of spots whose top type is Epithelial")
ax[1].set_title("(b) Two ways of picking the subsets, per slide\n"
                "if the bars diverge, the gate is depth-sensitive here", fontsize=10)
ax[1].legend(fontsize=7, frameon=False)

# ---- (c) 「并列有多普遍」 -------------------------------------------------
ax[2].hist(100 * df.stable_frac, bins=20, color="#27ae60", alpha=0.8,
           label=f"spots with margin >= {STABLE}")
ax[2].hist(100 * (1 - df.tight_frac), bins=20, color="#e67e22", alpha=0.6,
           label=f"spots with margin >= {TIGHT}")
ax[2].set_xlabel("% of spots on the slide")
ax[2].set_ylabel("number of slides")
ax[2].set_title(f"(c) How much of each slide sits away from a tie\n"
                f"{len(df)} slides", fontsize=10)
ax[2].legend(fontsize=8, frameon=False)

for a_ in ax:
    a_.spines[["top", "right"]].set_visible(False)

fig.suptitle("Cohort-wide margin distribution  —  'how close is the winning cell type?'  "
             "(RCTD_PREREG §13.5.4)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT_PNG, dpi=145)

print("表写到:", OUT_TSV)
print("图写到:", OUT_PNG)
print()
print(f"切片数 {len(df)}")
print(f"上皮占比（全体 argmax）  中位 {100*df.epi_frac.median():.2f}%  "
      f"范围 {100*df.epi_frac.min():.2f}–{100*df.epi_frac.max():.2f}%")
print(f"上皮占比（margin>={STABLE}）中位 {100*df.epi_frac_stable.median():.2f}%  "
      f"范围 {100*df.epi_frac_stable.min():.2f}–{100*df.epi_frac_stable.max():.2f}%")
print(f"margin < {TIGHT} 的 spot 占比  中位 {100*df.tight_frac.median():.1f}%  "
      f"范围 {100*df.tight_frac.min():.1f}–{100*df.tight_frac.max():.1f}%")
print(f"margin >= {STABLE} 的 spot 占比 中位 {100*df.stable_frac.median():.1f}%  "
      f"范围 {100*df.stable_frac.min():.1f}–{100*df.stable_frac.max():.1f}%")
print(f"逐张 margin 中位数          中位 {df.margin_p50.median():.3f}  "
      f"范围 {df.margin_p50.min():.3f}–{df.margin_p50.max():.3f}")
