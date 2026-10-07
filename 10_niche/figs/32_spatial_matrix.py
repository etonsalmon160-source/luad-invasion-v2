#!/usr/bin/env python3
# 32_spatial_matrix.py —— 空间矩阵：行 = 五期代表切片
#   列 = H&E | RCTD(6 谱系) | RCTD(39 亚型) | 四轴耗竭 | 生态位域类型
#   RCTD 两版都用荧光加性配色；耗竭做空间邻域平滑后按主导轴上色
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, json, colorsys, pickle, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
WHITE = "--white" in sys.argv            # 只把**整图背景**改白；面板内部保持黑
SUFFIX = "_white" if WHITE else ""
PAGE_BG = "white" if WHITE else "black"    # 整图背景（间隙/页边/图例区）
PANEL_BG = "black"                        # 面板内部**始终黑**（荧光加性配色需要）
TXT = "black" if WHITE else "white"
REP = [("Normal", "GSM9226174_P4_Normal"), ("AAH", "GSM9226222_P25_AAH"),
       ("AIS", "GSM9226207_P19_AIS"), ("MIA", "GSM9226195_P13_MIA"),
       ("IAC", "GSM9226200_P15_LUAD")]

EXH = {"T cell":  ["PDCD1","LAG3","HAVCR2","TIGIT","ENTPD1","TOX","CXCL13","LAYN","TNFRSF9"],
       "NK":      ["KLRC1","KIR2DL3","TIGIT","HAVCR2","PDCD1","CD160","KIR3DL1"],
       "B cell":  ["FCRL4","ITGAX","ZEB2","FCRL2","CD86","FCRL5"],
       "LAM":     ["TREM2","APOE","GPNMB","SPP1","LPL","LGALS3","CD9"]}   # 脂质相关巨噬（项目 Mac_LAM）
EXH_NAMES = list(EXH)
CH = np.array([[0.95,0.20,0.20], [0.20,0.85,0.35], [0.25,0.45,1.00], [1.00,0.65,0.10]])

# ── RCTD：39 亚型荧光配色 + 6 谱系配色 ──
map6 = pd.read_csv("/tmp/_rctd2lin.csv")
REF = sorted(map6.rctd.unique())
r2l = dict(zip(map6.rctd, map6.L1))
LIN_ORDER = ["上皮", "成纤维", "髓系", "内皮", "T/NK", "B/浆"]
LIN_EN = {"上皮":"Epithelial", "成纤维":"Fibroblast", "髓系":"Myeloid",
          "内皮":"Endothelial", "T/NK":"T/NK", "B/浆":"B/Plasma"}
HUE  = {"上皮":12, "成纤维":52, "髓系":125, "内皮":182, "T/NK":235, "B/浆":302}
SPAN = {"上皮":26, "成纤维":16, "髓系":62, "内皮":26, "T/NK":30, "B/浆":26}
FLUO = {}
for L in LIN_ORDER:
    sub = [r for r in REF if r2l.get(r) == L]
    n = len(sub)
    off = np.linspace(-SPAN[L]/2, SPAN[L]/2, n) if n > 1 else np.array([0.0])
    order = np.argsort(np.argsort(off))
    v = np.linspace(0.80, 1.0, n)[order] if n > 1 else np.array([1.0])
    for i, s_ in enumerate(sub):
        FLUO[s_] = np.array(colorsys.hsv_to_rgb(((HUE[L] + off[i]) % 360) / 360.0, 0.95, v[i]))
FLUO39 = np.array([FLUO[r] for r in REF])
LIN6 = np.array([colorsys.hsv_to_rgb(HUE[L]/360.0, 0.95, 0.98) for L in LIN_ORDER])
M6 = np.zeros((len(REF), 6))
for i, r in enumerate(REF):
    if r2l.get(r) in LIN_ORDER:
        M6[i, LIN_ORDER.index(r2l[r])] = 1.0

d7 = pd.read_csv(f"{ROOT}/results/10_niche/kstar_diag/d7_domain_assign.tsv", sep="\t")
D7MAP = {(s, int(d)): int(a) for s, d, a in zip(d7.slide, d7.domain, d7.archetype)}
KH = f"{ROOT}/results/10_niche/kstar_diag/d8_parts"


def load(sl):
    d = f"{ROOT}/data/visium_spatial/{sl}"
    feat = pd.read_csv(f"{d}/filtered_feature_bc_matrix/features.tsv.gz", sep="\t", header=None)
    bar  = pd.read_csv(f"{d}/filtered_feature_bc_matrix/barcodes.tsv.gz", header=None)[0].values
    M = sio.mmread(f"{d}/filtered_feature_bc_matrix/matrix.mtx.gz").tocsr()
    if M.shape[0] == len(feat) and M.shape[1] != len(feat):
        M = M.T.tocsr()
    gi = {g: i for i, g in enumerate(feat[1].astype(str).values)}
    tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
    Mn = (sp.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
    Zs, used = [], []
    for k in EXH_NAMES:
        ix = [gi[g] for g in EXH[k] if g in gi]
        used.append((k, len(ix)))
        X = np.asarray(Mn[:, ix].todense())
        Zs.append(((X - X.mean(0)) / (X.std(0) + 1e-9)).mean(1))
    Z = np.vstack(Zs).T
    W = pd.read_csv(f"{ROOT}/results/08_spatial_deconv/rctd_d/per_slide/{sl}.weights.tsv.gz",
                    sep="\t", index_col=0)
    W = W.reindex(index=bar, columns=REF).fillna(0).values
    pos = pd.read_csv(f"{d}/spatial/tissue_positions.csv")
    sf = json.load(open(f"{d}/spatial/scalefactors_json.json"))["tissue_hires_scalef"]
    he = imread(f"{d}/spatial/tissue_hires_image.png")
    pos = pos.drop_duplicates("barcode").set_index("barcode").reindex(bar)
    ok = pos.pxl_col_in_fullres.notna().values
    X = (pos.pxl_col_in_fullres.values * sf)[ok]
    Y = (-pos.pxl_row_in_fullres.values * sf)[ok]
    dm = pd.read_csv(f"{KH}/{sl}.tsv", sep="\t").drop_duplicates("barcode").set_index("barcode")
    dom = pd.Series(bar).map(dm.domain).values
    dom = np.array([D7MAP.get((sl, int(v)), np.nan) if not pd.isna(v) else np.nan for v in dom])
    return dict(he=he, W=he.shape[1], H=he.shape[0], X=X, Y=Y,
                Z=Z[ok], W39=W[ok], dom=dom.astype(float)[ok], used=used)


def smooth(Z, X, Y, k=7):
    """把每个 spot 的四轴 z 分与空间最近邻平均，去掉单 spot 噪点。"""
    t = cKDTree(np.c_[X, Y])
    _, idx = t.query(np.c_[X, Y], k=k)
    return Z[idx].mean(1)


def exh_rgb(Zs, thresh=0.15, scale=0.85):
    """每个 spot 取相对四轴均值最高的那一轴着色，亮度 = 主导程度。"""
    Zc = Zs - Zs.mean(1, keepdims=True)
    dom = Zc.argmax(1); mag = Zc.max(1)
    b = np.clip((mag - thresh) / scale, 0, 1)
    return np.clip(CH[dom] * b[:, None] + 0.05, 0, 1)


CACHE = "/tmp/_p5_matrix_data.pkl"
if os.path.exists(CACHE) and "--nocache" not in sys.argv:
    DATA = pickle.load(open(CACHE, "rb")); print("cache hit")
else:
    DATA = [load(sl) for _, sl in REP]
    pickle.dump(DATA, open(CACHE, "wb")); print("cache written")
for (stage, sl), D in zip(REP, DATA):
    print(f"{stage:7s} {sl:28s} spots={len(D['X'])}  " +
          " ".join(f"{k}={n}" for k, n in D["used"]))

DOMC = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
# 39 亚型图例要挤进单列宽度，长名字缩略（只影响图例文字，不影响数据）
ABBR = {
    "Basophil_Mast 1": "Basophil/Mast 1",
    "Vascular Smooth Muscle": "Vasc. Smooth Muscle",
    "Classical Monocyte": "Classical Mono.",
    "Nonclassical Monocyte": "Nonclass. Mono.",
    "Proliferating Macrophage": "Prolif. Macrophage",
    "Myeloid Dendritic Type 1": "Myeloid DC 1",
    "Myeloid Dendritic Type 2": "Myeloid DC 2",
    "Plasmacytoid Dendritic": "Plasmacytoid DC",
    "EREG+ Dendritic": "EREG+ DC",
    "IGSF21+ Dendritic": "IGSF21+ DC",
    "TREM2+ Dendritic": "TREM2+ DC",
    "Capillary Intermediate 1": "Capillary Interm. 1",
    "Bronchial Vessel 1": "Bronch. Vessel 1",
    "Bronchial Vessel 2": "Bronch. Vessel 2",
}
REF_BY_LIN = [r for L in LIN_ORDER for r in REF if r2l.get(r) == L]
DOMLAB = ["Airway (ciliated)", "iCAF", "ECM/interstitial", "Alveolar-capillary",
          "AT2", "Vascular", "Lymphoid"]
TITLES = ["H&E", "RCTD deconvolution\n6 lineages", "RCTD deconvolution\n39 subtypes",
          "Programs\n(T / NK / B / LAM)", "Niche domain (K*=7)"]
SPOT = 2.2


def px100(sl):
    """该切片 100 µm 对应的 hires 像素数（用相邻 spot 中位间距自标定，避开 µm/px 争议）。"""
    d = f"{ROOT}/data/visium_spatial/{sl}/spatial"
    pos = pd.read_csv(f"{d}/tissue_positions.csv")
    sf = json.load(open(f"{d}/scalefactors_json.json"))["tissue_hires_scalef"]
    pos = pos[pos.in_tissue == 1] if "in_tissue" in pos.columns else pos
    X = pos.pxl_col_in_fullres.values * sf
    Y = -pos.pxl_row_in_fullres.values * sf
    t = cKDTree(np.c_[X, Y]); dd, _ = t.query(np.c_[X, Y], k=2)
    return float(np.median(dd[:, 1]))


def draw_scalebar(ax, x0, x1, y0, y1, p100):
    """在面板左下角画比例尺；自动挑一个不超过面板宽 1/4 的整数长度。"""
    for lab, um in (("1 mm", 1000.0), ("500 µm", 500.0), ("200 µm", 200.0)):
        L = um / 100.0 * p100
        if L <= 0.24 * (x1 - x0):
            break
    xa = x0 + 0.045 * (x1 - x0)
    yb = y0 + 0.055 * (y1 - y0)
    # 白条 + 黑描边：白底 H&E 与黑底 panel 上都能看见
    ax.plot([xa, xa + L], [yb, yb], color="black", lw=3.6, solid_capstyle="butt", zorder=8)
    ax.plot([xa, xa + L], [yb, yb], color="white", lw=2.1, solid_capstyle="butt", zorder=9)
    ax.text(xa + L / 2, yb + 0.018 * (y1 - y0), lab, color="white", fontsize=5.4,
            ha="center", va="bottom", zorder=9,
            bbox=dict(facecolor="0.15", edgecolor="none", alpha=0.75, pad=0.7))

def he_box(D):
    """H&E 图里实际的组织范围（非白像素），返回 hires 坐标系下的 (x0,x1,y0,y1)。"""
    im = D["he"]
    rgb = im[..., :3] if im.ndim == 3 else np.dstack([im] * 3)
    m = (rgb.min(axis=2) < 0.88) & (rgb.max(axis=2) - rgb.min(axis=2) > 0.04)
    rows = np.where(m.sum(1) > 8)[0]; cols = np.where(m.sum(0) > 8)[0]
    if len(rows) == 0 or len(cols) == 0:
        return None
    return (float(cols.min()), float(cols.max()), float(-rows.max()), float(-rows.min()))


# 限框 = spot 外接框 ∪ H&E 组织框（只用 spot 框会把 H&E 里的组织切掉）
BOX = []
for D in DATA:
    x0, x1 = float(D["X"].min()), float(D["X"].max())
    y0, y1 = float(D["Y"].min()), float(D["Y"].max())
    hb = he_box(D)
    if hb is not None:
        x0, x1 = min(x0, hb[0]), max(x1, hb[1])
        y0, y1 = min(y0, hb[2]), max(y1, hb[3])
    m = 0.02 * max(x1 - x0, y1 - y0)
    BOX.append((x0 - m, x1 + m, y0 - m, y1 + m))
# 再把所有格子统一到同一个长宽比（只补黑边、不拉伸组织）⇒ 网格方正
_A = float(np.mean([(b[3] - b[2]) / (b[1] - b[0]) for b in BOX]))
_FIT = []
for D, (x0, x1, y0, y1) in zip(DATA, BOX):
    w, h = x1 - x0, y1 - y0
    if h / w < _A:
        d = (_A * w - h) / 2; y0 -= d; y1 += d
    else:
        d = (h / _A - w) / 2; x0 -= d; x1 += d
    # 裁到 H&E 图片本身的范围：越界会在面板边缘留出黑带（看着就像和转录组错位）
    x0 = max(x0, 0.0); x1 = min(x1, float(D["W"]))
    y0 = max(y0, -float(D["H"])); y1 = min(y1, 0.0)
    _FIT.append((x0, x1, y0, y1))
BOX = _FIT
hr = [(b[3] - b[2]) / (b[1] - b[0]) for b in BOX]
SCALE = [px100(sl) for _, sl in REP]
print(f"统一长宽比 = {_A:.3f}")
fig = plt.figure(figsize=(9.6, 9.6 * sum(hr) * 0.20 + 1.15))
fig.patch.set_facecolor(PAGE_BG)
outer = fig.add_gridspec(2, 1, height_ratios=[sum(hr), 1.00], hspace=0.002,
                         left=0.035, right=0.995, top=0.945, bottom=0.012)
gs = outer[0].subgridspec(5, 5, height_ratios=hr, wspace=0.002, hspace=0.002)
gsl = outer[1].subgridspec(1, 5, height_ratios=[1.0], hspace=0.008)

for r, ((stage, sl), D) in enumerate(zip(REP, DATA)):
    W6 = D["W39"] @ M6
    rgb6 = np.clip(W6 @ LIN6, 0, 1)
    rgb39 = np.clip(D["W39"] @ FLUO39, 0, 1)
    Zs = smooth(D["Z"], D["X"], D["Y"], k=37)
    rgb_ex = exh_rgb(Zs)
    for c in range(5):
        ax = fig.add_subplot(gs[r, c])
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.35); s.set_color("0.45" if WHITE else "0.30")
        ax.set_facecolor(PANEL_BG)
        if c == 0:
            ax.imshow(D["he"], extent=[0, D["W"], -D["H"], 0])
            _b = BOX[r]
            draw_scalebar(ax, _b[0], _b[1], _b[2], _b[3], SCALE[r])
        elif c == 1:
            ax.scatter(D["X"], D["Y"], c=rgb6, s=SPOT, marker="s", linewidths=0, rasterized=True)
        elif c == 2:
            ax.scatter(D["X"], D["Y"], c=rgb39, s=SPOT, marker="s", linewidths=0, rasterized=True)
        elif c == 3:
            ax.scatter(D["X"], D["Y"], c=rgb_ex, s=SPOT, marker="s", linewidths=0, rasterized=True)
        else:
            a = D["dom"]; ok = ~pd.isna(a)
            ax.scatter(D["X"][ok], D["Y"][ok], c=[DOMC[int(v)-1] for v in a[ok]],
                       s=SPOT, marker="s", linewidths=0, rasterized=True)
        bx0, bx1, by0, by1 = BOX[r]
        ax.set_xlim(bx0, bx1); ax.set_ylim(by0, by1); ax.set_aspect("equal")
        if r == 0:
            ax.set_title(TITLES[c], fontsize=7.2, pad=3, fontweight="bold",
                         linespacing=1.25, color=TXT)
        if c == 0:
            ax.set_ylabel(stage, fontsize=8.5, fontweight="bold", labelpad=3, color=TXT)

# ── 每列下方的图例（用独立轴排版，保证对齐）──
def leg_ax(r, c, ncol=1):
    a = fig.add_subplot(gsl[r, c]); a.axis("off"); a.set_facecolor(PAGE_BG)
    return a


def style_leg(leg):
    leg.get_frame().set_facecolor("none")
    leg.get_title().set_color(TXT)
    for t in leg.get_texts():
        t.set_color(TXT)
    return leg

LT = dict(marker="s", ls="", ms=5.2, mec="none")   # 四组图例统一色块尺寸
h_lin = [Line2D([], [], mfc=LIN6[i], label=LIN_EN[L], **LT) for i, L in enumerate(LIN_ORDER)]
ax = leg_ax(0, 1)
style_leg(ax.legend(handles=h_lin, loc="upper center", ncol=1, frameon=False,
                    fontsize=6.0, handletextpad=0.35, labelspacing=0.18,
                    title="Lineage colour", title_fontsize=6.3))

h_ex = [Line2D([], [], mfc=CH[i], label=k, **LT) for i, k in enumerate(EXH_NAMES)]
ax = leg_ax(0, 3)
style_leg(ax.legend(handles=h_ex, loc="upper center", ncol=1, frameon=False,
                    fontsize=6.0, handletextpad=0.35, labelspacing=0.18,
                    title="Programs\n(dominant one per spot)", title_fontsize=6.3))

h_dm = [Line2D([], [], mfc=DOMC[i], label=f"D{i+1} {DOMLAB[i]}", **LT) for i in range(7)]
ax = leg_ax(0, 4)
style_leg(ax.legend(handles=h_dm, loc="upper center", ncol=1, frameon=False,
                    fontsize=6.0, handletextpad=0.35, labelspacing=0.18,
                    title="Niche domain (K*=7)", title_fontsize=6.3))

# 39 亚型图例：塞进它自己那一列，2 列排布，长名缩略
LT2 = dict(marker="s", ls="", ms=5.2, mec="none")
h39 = [Line2D([], [], mfc=FLUO[r], label=ABBR.get(r, r), **LT2) for r in REF_BY_LIN]
ax = leg_ax(0, 2)
style_leg(ax.legend(handles=h39, loc="upper center", ncol=2, frameon=False,
                    fontsize=4.3, handletextpad=0.26, labelspacing=0.15,
                    columnspacing=0.55, handlelength=0.9,
                    title="39 subtypes\n(ordered by lineage, n=39)", title_fontsize=6.0))
leg_ax(0, 0)   # H&E 列无图例，占位保持列宽一致

fig.savefig(f"{OUT}/P5_spatial_matrix{SUFFIX}.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P5_spatial_matrix{SUFFIX}.png", dpi=400, bbox_inches="tight")
print("✓ P5_spatial_matrix{SUFFIX}")
