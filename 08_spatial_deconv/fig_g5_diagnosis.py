#!/usr/bin/env python3
"""fig_g5_diagnosis.py —— G5 硬门（|Spearman(spot 总 UMI, 权重行和)| <= 0.5）为什么在第 2 张就触发

背景：56 张全跑在 GSM9226169_P1_LUAD 上硬停（rho = 0.521 > 0.5）。
     按法则 3.2 不事后调阈值；本图只为**判根因**，不给「改成多少」的建议。

要判的是：G5 的这个 rho 是不是在测它标签上写的东西（「权重随测序深度漂移」）。
三个对照：
  (a) 深度本身在空间上成不成组织样的团块 —— 技术伪影不该有平滑的空间结构；
  (b) UMI 与 n_genes 是不是同一个变量 —— 若是，rho 测的其实是「mRNA 总量 vs 分到的细胞质量」；
  (c) rho 随 x 的离散程度怎么变 —— 决定它能不能跨切片比。

图内标签全英文（本机无 CJK 字体）。

跑法：python3 08_spatial_deconv/fig_g5_diagnosis.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = "/home/eto/luad_v2"
S1 = "GSM9226168_P1_AAH"    # 过门 (rho=0.437)
S2 = "GSM9226169_P1_LUAD"   # 触发 (rho=0.521)
GATE = 0.5

SM = pd.read_csv(f"{ROOT}/results/08_spatial_deconv/spot_mask.tsv.gz", sep="\t")
SM["pass"] = SM["pass"].astype(str).str.lower().isin(["true", "1"])

W = pd.read_csv(f"{ROOT}/results/08_spatial_deconv/rctd_a/per_slide/{S1}.weights.tsv.gz",
                sep="\t", index_col=0)
D1 = SM[(SM.slide == S1) & SM["pass"]].set_index("barcode")
common = [b for b in W.index if b in D1.index]
umi = D1.loc[common, "total_umi"].values.astype(float)
ng = D1.loc[common, "n_genes"].values.astype(float)
rs = W.loc[common].sum(1).values

pos = pd.read_csv(f"{ROOT}/data/visium_spatial/{S1}/spatial/tissue_positions.csv")
pos = pos[pos.barcode.isin(common)].set_index("barcode").loc[common]
xy = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
_, idx = cKDTree(xy).query(xy, k=7)
nb = idx[:, 1:]


def moran_like(v):
    z = (v - v.mean()) / (v.std() + 1e-12)
    return float(np.mean(z * z[nb].mean(1)))


rho_all = spearmanr(umi, rs).correlation
q25, q75 = np.percentile(umi, [25, 75])
m = (umi >= q25) & (umi <= q75)
rho_mid = spearmanr(umi[m], rs[m]).correlation
I_umi, I_rs = moran_like(umi), moran_like(rs)

fig, ax = plt.subplots(2, 2, figsize=(13.0, 10.0))

# ---- (0,0) 两张切片的深度分布 --------------------------------------------
bins = np.linspace(2.0, 5.2, 70)
for s, col in ((S1, "#2980b9"), (S2, "#c0392b")):
    d = SM[(SM.slide == s) & SM["pass"]]
    ax[0, 0].hist(np.log10(d.total_umi.values), bins=bins, histtype="step",
                  lw=2.0, color=col,
                  label=f"{s.split('_',1)[1]}  n={len(d):,}  median={d.total_umi.median():.0f}")
    ax[0, 0].axvline(np.log10(d.total_umi.median()), color=col, lw=1.0, ls=":")
ax[0, 0].set_xlabel("log10(total UMI per spot)")
ax[0, 0].set_ylabel("number of spots")
ax[0, 0].set_title("(a) The two slides are not comparable in depth spread\n"
                   "the failing slide has a WIDER spread and 42% fewer spots", fontsize=10)
ax[0, 0].legend(fontsize=8, frameon=False)

# ---- (0,1) UMI vs 行和 ----------------------------------------------------
hb = ax[0, 1].hexbin(np.log10(umi), rs, gridsize=48, cmap="Blues", mincnt=1)
ax[0, 1].axhline(1.0, color="#7f8c8d", lw=1.2, ls="--")
ax[0, 1].axhline(np.median(rs), color="#e67e22", lw=1.6, ls="-.",
                 label=f"median row-sum = {np.median(rs):.3f}  (a true proportion would be 1.0)")
ax[0, 1].set_xlabel("log10(total UMI per spot)")
ax[0, 1].set_ylabel("sum of RCTD weights per spot")
ax[0, 1].set_title(f"(b) G5's statistic, slide that PASSED\n"
                   f"rho = {rho_all:+.3f} over all spots   |   "
                   f"rho = {rho_mid:+.3f} within P25-P75", fontsize=10)
ax[0, 1].legend(fontsize=8, frameon=False)
plt.colorbar(hb, ax=ax[0, 1], label="spots")
for a in ax.ravel():
    a.spines[["top", "right"]].set_visible(False)

# ---- (1,0)/(1,1) 空间图：深度与行和 -------------------------------------
for a, v, name, I in ((ax[1, 0], np.log10(umi), "log10(total UMI)", I_umi),
                      (ax[1, 1], rs, "sum of RCTD weights", I_rs)):
    sc = a.scatter(xy[:, 0], xy[:, 1], c=v, s=3.0, cmap="viridis", linewidths=0)
    a.invert_yaxis()
    a.set_aspect("equal")
    a.set_xticks([]); a.set_yticks([])
    a.set_title(f"(c) {name}, mapped on the tissue\n"
                f"neighbour correlation = {I:+.3f}  "
                f"({'tissue-structured' if I > 0.5 else 'not structured'})", fontsize=10)
    plt.colorbar(sc, ax=a, label=name.split("(")[0].strip())
    a.spines[["top", "right", "left", "bottom"]].set_visible(False)

fig.suptitle("Why G5 fired on slide 2  —  diagnosing the criterion, not re-tuning the threshold",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.955])
OUT = f"{ROOT}/results/08_spatial_deconv/rctd_a/fig_g5_diagnosis.png"
fig.savefig(OUT, dpi=145)
print("图写到:", OUT)
print()
print(f"G5 全范围        rho = {rho_all:+.3f}")
print(f"G5 P25-P75       rho = {rho_mid:+.3f}")
print(f"UMI 空间结构     I   = {I_umi:+.3f}")
print(f"行和 空间结构    I   = {I_rs:+.3f}")
print(f"rho(UMI,n_genes)     = {spearmanr(umi, ng).correlation:+.3f}")
print(f"行和 P5/P50/P95      = {np.round(np.percentile(rs, [5, 50, 95]), 3)}")
