#!/usr/bin/env python3
# 36_domain_marker_heatmap.py —— 七个生态位域类型的标志基因表达热图
#   数据源：results/10_niche/kstar_diag/d13_parts/<slide>.tsv（逐片、逐域的 HVG 均值表达）
#   行序 = 该片 domain 升序（已用「行数 == d7 域数」在 56/56 上核过），
#   并会用 d13_archetype_markers.tsv 交叉验证 top marker 是否复现。
import numpy as np, pandas as pd, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT = "/home/eto/luad_v2"; T = f"{ROOT}/results/10_niche"
OUT = f"{ROOT}/results/paper_figures"
PARTS = f"{T}/kstar_diag/d13_parts"

DOMC  = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
DOMLAB = ["Airway", "iCAF", "ECM/interstitial", "Alveolar-cap.",
          "AT2", "Vascular", "Lymphoid"]

d7 = pd.read_csv(f"{T}/kstar_diag/d7_domain_assign.tsv", sep="\t")
D7 = {(s, int(d)): int(a) for s, d, a in zip(d7.slide, d7.domain, d7.archetype)}
DOMS = d7.groupby("slide").domain.apply(lambda x: sorted(set(x))).to_dict()

GENES = None
SUM = None; CNT = np.zeros(7)
for f in sorted(glob.glob(f"{PARTS}/*.tsv")):
    sl = os.path.basename(f)[:-4]
    D = pd.read_csv(f, sep="\t")                   # 列序：<基因...> domain slide
    if GENES is None:
        GENES = [c for c in D.columns if c not in ("domain", "slide")]
        SUM = np.zeros((7, len(GENES)))
    V = D[GENES].values.astype(float)
    for i, dm in enumerate(D["domain"].astype(int).values):
        a = D7.get((sl, int(dm)))
        if a is None:
            continue
        SUM[a - 1] += V[i]; CNT[a - 1] += 1
PROF = SUM / np.maximum(CNT, 1)[:, None]
print(f"域×基因画像 {PROF.shape} | 每域参与的切片数: {CNT.astype(int)}")

# ── 交叉验证：重算的 top30 marker 是否复现 d13_archetype_markers.tsv ──
eps = 1e-4
rest = (PROF.sum(0, keepdims=True) - PROF) / 6.0
lfc = np.log2((PROF + eps) / (rest + eps))
ref = pd.read_csv(f"{T}/kstar_diag/d13_archetype_markers.tsv", sep="\t")
ov = []
for a in range(1, 8):
    mine = [GENES[j] for j in np.argsort(-lfc[a - 1])[:30]]
    theirs = list(ref[ref.archetype == a].sort_values("rank").gene)
    ov.append(len(set(mine) & set(theirs)) / 30)
print("与 d13 的 top30 marker 重合率:", "  ".join(f"D{i+1}={v:.2f}" for i, v in enumerate(ov)))

# ── 选基因：每域 top6 ──
TOP, blocks = [], []
for a in range(7):
    for j in np.argsort(-lfc[a])[:6]:
        g = GENES[j]
        if g not in TOP:
            TOP.append(g); blocks.append(a)
print(f"选中 {len(TOP)} 个基因")

Y = PROF[:, [GENES.index(g) for g in TOP]]      # 7 域 × 42 基因
Z = (Y - Y.mean(0, keepdims=True)) / (Y.std(0, keepdims=True) + 1e-9)  # 每个基因在 7 域间标准化
Z = Z.T                                          # → 42 基因 × 7 域（行=基因）

cmap = LinearSegmentedColormap.from_list("bwr2", ["#4A7EBB", "#F6F6F4", "#C0392B"])
fig = plt.figure(figsize=(5.8, 0.168 * len(TOP) + 1.5))
gs = fig.add_gridspec(1, 2, width_ratios=[0.62, 1.0], wspace=0.42,
                      left=0.02, right=0.90, top=0.945, bottom=0.135)


def text_col(hexc):
    """浅色域（如 D4 的淡粉）直接当文字色看不清，压暗后再用。"""
    r, g, b = [int(hexc[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    if 0.299 * r + 0.587 * g + 0.114 * b > 0.72:
        r, g, b = r * 0.45, g * 0.45, b * 0.45
    return (r, g, b)
axb = fig.add_subplot(gs[0, 0])          # 左栏：域名字 + 色块（图例）
ax = fig.add_subplot(gs[0, 1])           # 热图主体

# ── 左栏图例：每个域的色块，左侧写域名字 ──
axb.set_xlim(0, 1); axb.set_ylim(len(TOP) - 0.5, -0.5); axb.axis("off")
runs, s0 = [], 0
for i in range(1, len(blocks) + 1):
    if i == len(blocks) or blocks[i] != blocks[s0]:
        runs.append((s0, i - 1, blocks[s0])); s0 = i
for (a, b, d) in runs:
    axb.add_patch(Rectangle((0.74, a - 0.5), 0.26, b - a + 1, color=DOMC[d], lw=0))
    axb.text(0.68, (a + b) / 2, f"D{d+1}  {DOMLAB[d]}", ha="right", va="center",
             fontsize=6.3, color=text_col(DOMC[d]), fontweight="bold")

# ── 热图主体 ──
im = ax.imshow(Z, aspect="auto", cmap=cmap, vmin=-2.2, vmax=2.2)
ax.set_xticks(range(7)); ax.set_xticklabels([f"D{i+1}" for i in range(7)], fontsize=7)
ax.set_yticks(range(len(TOP)))
ax.set_yticklabels(TOP, fontsize=5.8, fontstyle="italic", color="0.15")
ax.set_ylim(len(TOP) - 0.5, -0.5)
ax.tick_params(axis="y", length=0, pad=3); ax.tick_params(axis="x", length=0, pad=4)
for s in ax.spines.values():
    s.set_linewidth(0.5); s.set_color("0.55")
for i, (a, b, d) in enumerate(runs):
    if i:
        ax.axhline(a - 0.5, color="white", linewidth=1.2)
cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.025)
cb.set_label("Row-scaled mean expression\n(CP10K, log1p)", fontsize=6, linespacing=1.35)
cb.ax.tick_params(labelsize=5.5, length=2); cb.outline.set_linewidth(0.4)
ax.set_title("Marker genes of the seven niche domains", fontsize=8.8, fontweight="bold", pad=10)
fig.savefig(f"{OUT}/P7_domain_marker_heatmap.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P7_domain_marker_heatmap.png", dpi=400, bbox_inches="tight")
print("✓ P7_domain_marker_heatmap")

pd.DataFrame({"archetype": [blocks[i] + 1 for i in range(len(TOP))], "gene": TOP,
              "log2fc": [round(lfc[blocks[i], GENES.index(TOP[i])], 3) for i in range(len(TOP))],
              **{f"D{j+1}_mean": Y[j].round(4) for j in range(7)}}).to_csv(
    f"{T}/domain_marker_heatmap_data.tsv", sep="\t", index=False)
np.save("/tmp/_domain_prof.npy", PROF); np.save("/tmp/_domain_genes.npy", np.array(GENES))
print("saved results/10_niche/domain_marker_heatmap_data.tsv")
