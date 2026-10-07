#!/usr/bin/env python3
# 37_domain_dotplot.py —— 七个生态位域 × 标志基因 的 dot plot
#   颜色 = 域内平均表达（行内 z 后线性映射到 0–1）  大小 = 域内表达该基因的 spot 比例
#   比例需要回原始矩阵统计 ⇒ 56 张切片用多进程并行（否则单线程 ~20 分钟）
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, glob, os, pickle
from multiprocessing import Pool
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

ROOT = "/home/eto/luad_v2"; T = f"{ROOT}/results/10_niche"
OUT = f"{ROOT}/results/paper_figures"
KH = f"{T}/kstar_diag/d8_parts"; PARTS = f"{T}/kstar_diag/d13_parts"
FRAC_CACHE = "/tmp/_domain_frac.pkl"

DOMC = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
DOMLAB = ["Airway", "iCAF", "ECM/interstitial", "Alveolar-cap.", "AT2", "Vascular", "Lymphoid"]

d7 = pd.read_csv(f"{T}/kstar_diag/d7_domain_assign.tsv", sep="\t")
D7 = {(s, int(d)): int(a) for s, d, a in zip(d7.slide, d7.domain, d7.archetype)}

# ── ① 域 × 基因 平均表达（快：直接用 d13_parts）──
GENES = None; SUM = None; CNT = np.zeros(7)
for f in sorted(glob.glob(f"{PARTS}/*.tsv")):
    sl = os.path.basename(f)[:-4]
    D = pd.read_csv(f, sep="\t")
    if GENES is None:
        GENES = [c for c in D.columns if c not in ("domain", "slide")]
        SUM = np.zeros((7, len(GENES)))
    V = D[GENES].values.astype(float)
    for i, dm in enumerate(D["domain"].astype(int).values):
        a = D7.get((sl, int(dm)))
        if a:
            SUM[a - 1] += V[i]; CNT[a - 1] += 1
PROF = SUM / np.maximum(CNT, 1)[:, None]

eps = 1e-4
rest = (PROF.sum(0, keepdims=True) - PROF) / 6.0
lfc = np.log2((PROF + eps) / (rest + eps))
TOP, blocks = [], []
for a in range(7):
    for j in np.argsort(-lfc[a])[:6]:
        if GENES[j] not in TOP:
            TOP.append(GENES[j]); blocks.append(a)
IDX = [GENES.index(g) for g in TOP]
print(f"选中 {len(TOP)} 个基因")

# ── ② 域 × 基因 表达 spot 比例（回原始矩阵，多进程并行）──
def frac_one(sl):
    d = f"{ROOT}/data/visium_spatial/{sl}"
    feat = pd.read_csv(f"{d}/filtered_feature_bc_matrix/features.tsv.gz", sep="\t", header=None)
    bar = pd.read_csv(f"{d}/filtered_feature_bc_matrix/barcodes.tsv.gz", header=None)[0].values
    M = sio.mmread(f"{d}/filtered_feature_bc_matrix/matrix.mtx.gz").tocsr()
    if M.shape[0] == len(feat) and M.shape[1] != len(feat):
        M = M.T.tocsr()
    gi = {g: i for i, g in enumerate(feat[1].astype(str).values)}
    keep = [(k, gi[g]) for k, g in enumerate(TOP) if g in gi]
    dm = pd.read_csv(f"{KH}/{sl}.tsv", sep="\t").drop_duplicates("barcode").set_index("barcode")
    dom = pd.Series(bar).map(dm.domain).values
    arch = np.array([D7.get((sl, int(v)), 0) if not pd.isna(v) else 0 for v in dom])
    num = np.zeros((7, len(TOP))); den = np.zeros(7)
    if keep:
        B = M[:, [c for _, c in keep]]
        B = (B > 0).astype(np.int8)
        for a in range(1, 8):
            m = arch == a
            if m.sum() == 0:
                continue
            num[a - 1, [k for k, _ in keep]] += np.asarray(B[m].sum(0)).ravel()
            den[a - 1] += m.sum()
    return num, den

if os.path.exists(FRAC_CACHE):
    FRAC = pickle.load(open(FRAC_CACHE, "rb")); print("比例 cache hit")
else:
    files = sorted(glob.glob(f"{KH}/*.tsv"))
    with Pool(10) as p:
        res = p.map(frac_one, [os.path.basename(f)[:-4] for f in files])
    NUM = sum(r[0] for r in res); DEN = sum(r[1] for r in res)
    FRAC = NUM / np.maximum(DEN, 1)[:, None]
    pickle.dump(FRAC, open(FRAC_CACHE, "wb")); print("比例 cache written")
    print("每域 spot 数:", DEN.astype(int))
print("每域参与切片数:", CNT.astype(int))

# ── ③ dot plot ──
Y = PROF[:, IDX]                                    # 7 × 42
Z = (Y - Y.mean(0, keepdims=True)) / (Y.std(0, keepdims=True) + 1e-9)
Zc = (Z - (-2.2)) / 4.4                             # → 0..1，与色阶同界
F = FRAC                                          # 已是 7 × 42（按 TOP 顺序）
size = 4 + 62 * np.sqrt(np.clip(F, 0, 1))

cmap = LinearSegmentedColormap.from_list("vir", ["#440154", "#3B528B", "#21918C", "#5EC962", "#FDE725"])
fig = plt.figure(figsize=(8.8, 4.8))
gs = fig.add_gridspec(1, 2, width_ratios=[0.30, 1.0], wspace=0.02,
                      left=0.05, right=0.855, top=0.975, bottom=0.255)
axl = fig.add_subplot(gs[0, 0]); ax = fig.add_subplot(gs[0, 1])

# 左栏：域色块 + 名字
axl.set_xlim(0, 1); axl.set_ylim(6.5, -0.5); axl.axis("off")
for a in range(7):
    axl.add_patch(plt.Rectangle((0.86, a - 0.34), 0.14, 0.68, color=DOMC[a], lw=0))
    c = DOMC[a]
    if 0.299 * int(c[1:3], 16) / 255 + 0.587 * int(c[3:5], 16) / 255 + 0.114 * int(c[5:7], 16) / 255 > 0.72:
        c = "#8a7f85"                      # D4 太浅，压暗再写字
    axl.text(0.82, a, f"D{a+1}  {DOMLAB[a]}", ha="right", va="center",
             fontsize=6.8, fontweight="bold", color=c)

xx, yy = np.meshgrid(np.arange(len(TOP)), np.arange(7))
ax.scatter(xx.ravel(), yy.ravel(), s=size.ravel(), c=Zc.ravel(), cmap=cmap, vmin=0, vmax=1,
           linewidths=0)
ax.set_xlim(-0.6, len(TOP) - 0.4); ax.set_ylim(6.6, -0.6)
ax.set_yticks(range(7)); ax.set_yticklabels([])
ax.set_xticks(range(len(TOP)))
ax.set_xticklabels(TOP, rotation=90, fontsize=5.6, fontstyle="italic")
ax.tick_params(axis="x", length=0, pad=4)
for s in ax.spines.values():
    s.set_linewidth(0.5); s.set_color("0.6")
prev = blocks[0]
for i, b in enumerate(blocks):
    if b != prev:
        ax.axvline(i - 0.5, color="0.75", linewidth=0.7, zorder=0)
    prev = b
# 色标 + 大小图例都竖排在右侧：色标占上半，大小图例接在下面
cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap=cmap), ax=ax,
                  fraction=0.030, pad=0.018, shrink=0.52, anchor=(0.0, 1.0))
cb.set_label("Scaled expression", fontsize=6.5)
cb.ax.tick_params(labelsize=5.5, length=2); cb.outline.set_linewidth(0.4)
fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=np.sqrt(4 + 62 * np.sqrt(v)),
                           mfc="0.35", mec="none", label=f"{int(v*100)}%")
                    for v in (0.1, 0.4, 0.8)],
           loc="upper left", bbox_to_anchor=(0.872, 0.42), frameon=False, fontsize=6.2,
           labelspacing=1.7, handletextpad=1.0, title="Fraction of spots\nexpressing",
           title_fontsize=6.5, alignment="left")
fig.savefig(f"{OUT}/P8_domain_dotplot.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P8_domain_dotplot.png", dpi=400, bbox_inches="tight")
print("✓ P8_domain_dotplot")
