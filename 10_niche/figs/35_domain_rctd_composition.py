#!/usr/bin/env python3
# 35_domain_rctd_composition.py —— 每个生态位域类型(D1..D7)的 RCTD 权重组分
#   (a) 6 谱系堆叠柱   (b) 39 亚型堆叠柱
import numpy as np, pandas as pd, colorsys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
T = f"{ROOT}/results/10_niche"

LIN_ORDER = ["上皮", "成纤维", "髓系", "内皮", "T/NK", "B/浆"]
LIN_EN = {"上皮":"Epithelial","成纤维":"Fibroblast","髓系":"Myeloid",
          "内皮":"Endothelial","T/NK":"T/NK","B/浆":"B/Plasma"}
LIN_COL = {"Epithelial":"#E64B35","Fibroblast":"#4DBBD5","Endothelial":"#00A087",
           "Myeloid":"#3C5488","T/NK":"#F39B7F","B/Plasma":"#8491B4"}
DOMLAB = ["D1 Airway (ciliated)", "D2 iCAF", "D3 ECM/interstitial",
          "D4 Alveolar-capillary", "D5 AT2", "D6 Vascular", "D7 Lymphoid"]

map6 = pd.read_csv("/tmp/_rctd2lin.csv")
r2l = dict(zip(map6.rctd, map6.L1))
REF = sorted(map6.rctd.unique())


def hx(h): return np.array(colorsys.rgb_to_hsv(*[int(h[i:i+2], 16)/255 for i in (1, 3, 5)]))
def l2_pal(hexc, n):
    h, s, v = hx(hexc)
    hs = [h] if n == 1 else [(h + o) % 1 for o in np.linspace(-0.052, 0.052, n)]
    ss = min(0.95, max(0.62, s)); vs = np.linspace(0.98, 0.68, n)
    return [colorsys.hsv_to_rgb(hs[i], ss, vs[i]) for i in range(n)]


REF_BY_LIN = [r for L in LIN_ORDER for r in REF if r2l.get(r) == L]
SUBC = {}
for L in LIN_ORDER:
    subs = [r for r in REF_BY_LIN if r2l.get(r) == L]
    for r, c in zip(subs, l2_pal(LIN_COL[LIN_EN[L]], len(subs))):
        SUBC[r] = c

nspot = pd.read_csv(f"{T}/domain_nspot.tsv", sep="\t", index_col=0)
L6 = pd.read_csv(f"{T}/domain_rctd_lineage_prop.tsv", sep="\t", index_col=0)
S39 = pd.read_csv(f"{T}/domain_rctd_subtype_prop.tsv", sep="\t", index_col=0)

x = np.arange(len(L6))
xlab = [f"D{int(a)}\n{nspot.loc[a,'pct']:.0f}%" for a in L6.index]

def stacked(df, cols, colmap, title, fname, w, h, ncol=2, fs=4.2, keyh=0.20):
    fig, ax = plt.subplots(figsize=(w, h))
    bottom = np.zeros(len(df))
    for c in cols:
        if c not in df.columns: continue
        ax.bar(x, df[c].values, bottom=bottom, width=0.72,
               color=colmap[c], edgecolor="white", linewidth=0.25, label=c)
        bottom += df[c].values
    ax.set_xticks(x); ax.set_xticklabels(xlab, fontsize=7)
    ax.set_ylim(0, 1); ax.set_ylabel("Mean RCTD weight fraction", fontsize=8)
    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.set_yticks(np.arange(0, 1.01, 0.25))
    ax.set_yticklabels([f"{int(v*100)}%" for v in np.arange(0, 1.01, 0.25)], fontsize=7)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.4); ax.set_axisbelow(True)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=fs, ncol=ncol, handleheight=keyh, labelspacing=0.25,
              handletextpad=0.4, columnspacing=0.9)
    fig.savefig(f"{OUT}/{fname}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/{fname}.png", dpi=400, bbox_inches="tight")
    print(f"  ✓ {fname}")

stacked(L6, [LIN_EN[L] for L in LIN_ORDER], LIN_COL,
        "RCTD lineage composition of each niche domain (56 sections)",
        "P6a_domain_RCTD_lineage_composition", 5.2, 3.6, ncol=1, fs=7)

stacked(S39, REF_BY_LIN, SUBC,
        "RCTD subtype composition of each niche domain (56 sections)",
        "P6b_domain_RCTD_subtype_composition", 6.4, 4.4, ncol=2, fs=4.0, keyh=0.18)
