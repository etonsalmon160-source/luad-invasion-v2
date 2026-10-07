#!/usr/bin/env python3
# 51_exhaustion_by_stage.py —— 四轴耗竭随病程的变化（堆叠构成 + 绝对水平）
#   每条轴的分数 = 该基因集在 spot 上的**平均表达**（CP10K+log1p），不是 z 分
#   ⇒ 可跨切片比较（与 19_signature_overlay.R 的口径一致：colMeans）
#   56 张切片用多进程并行
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, glob, os, pickle
from multiprocessing import Pool
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; T = f"{ROOT}/results/10_niche"
OUT = f"{ROOT}/results/paper_figures"
KH = f"{T}/kstar_diag/d8_parts"; CACHE = "/tmp/_exh_stage.pkl"
ST = ["Normal", "AAH", "AIS", "MIA", "IAC"]
ALIAS = {"LUAD": "IAC"}

EXH = {"T cell":  ["PDCD1","LAG3","HAVCR2","TIGIT","ENTPD1","TOX","CXCL13","LAYN","TNFRSF9"],
       "NK":      ["KLRC1","KIR2DL3","TIGIT","HAVCR2","PDCD1","CD160","KIR3DL1"],
       "B cell":  ["FCRL4","ITGAX","ZEB2","FCRL2","CD86","FCRL5"],
       "LAM":     ["TREM2","APOE","GPNMB","SPP1","LPL","LGALS3","CD9"]}   # 脂质相关巨噬程序（项目 Mac_LAM），非 T 细胞意义的耗竭
CH = ["#C0392B", "#2E8B8B", "#2F5597", "#E08A3C"]     # 与 P5 耗竭列同色


def one(sl):
    d = f"{ROOT}/data/visium_spatial/{sl}"
    feat = pd.read_csv(f"{d}/filtered_feature_bc_matrix/features.tsv.gz", sep="\t", header=None)
    M = sio.mmread(f"{d}/filtered_feature_bc_matrix/matrix.mtx.gz").tocsr()
    if M.shape[0] == len(feat) and M.shape[1] != len(feat):
        M = M.T.tocsr()
    gi = {g: i for i, g in enumerate(feat[1].astype(str).values)}
    tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
    Mn = sp.diags(1e4 / tot) @ M; Mn.data = np.log1p(Mn.data)
    vals = []
    for k in EXH:
        ix = [gi[g] for g in EXH[k] if g in gi]
        vals.append(float(np.asarray(Mn[:, ix].mean(1)).mean()) if ix else np.nan)
    return sl, vals


if os.path.exists(CACHE):
    S = pickle.load(open(CACHE, "rb")); print("cache hit")
else:
    files = sorted(glob.glob(f"{KH}/*.tsv"))
    with Pool(10) as p:
        res = p.map(one, [os.path.basename(f)[:-4] for f in files])
    S = pd.DataFrame([(s, *v) for s, v in res], columns=["slide", "T cell", "NK", "B cell", "LAM"])
    d7 = pd.read_csv(f"{T}/kstar_diag/d7_domain_assign.tsv", sep="\t")[["slide", "stage"]].drop_duplicates()
    S = S.merge(d7, on="slide")
    S["stage"] = S.stage.replace(ALIAS)
    pickle.dump(S, open(CACHE, "wb")); print("cache written")

S = S[S.stage.isin(ST)]
KEY = list(EXH)
per_slide = S.groupby(["slide", "stage"])[KEY].mean().reset_index()
per_stage = per_slide.groupby("stage")[KEY].mean().reindex(ST)
print("\n逐期别平均（CP10K+log1p 的基因集均值）：")
print(per_stage.round(4).to_string())

# 逐 spot 的相对构成：需要逐 spot 值 ⇒ 取每张切片的构成后按切片平均
comp = per_slide[KEY].div(per_slide[KEY].sum(1), axis=0)
comp["stage"] = per_slide["stage"].values
C = comp.groupby("stage")[KEY].mean().reindex(ST)
print("\n相对构成（每张切片归一后按切片平均）：")
print((C * 100).round(1).to_string())

# ── 图 ──
fig = plt.figure(figsize=(7.2, 2.9))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.24,
                      left=0.07, right=0.86, top=0.80, bottom=0.16)

ax = fig.add_subplot(gs[0, 0])
bot = np.zeros(len(ST))
xs = np.arange(len(ST))
for k, c in zip(KEY, CH):
    ax.bar(xs, C[k].values * 100, bottom=bot * 100, width=0.68, color=c, label=k, linewidth=0)
    bot += C[k].values
ax.set_xticks(xs); ax.set_xticklabels(ST, fontsize=7)
ax.set_ylim(0, 100); ax.set_ylabel("Relative composition of the four programs  (%)", fontsize=6.8)
ax.tick_params(axis="y", labelsize=6.2, length=2.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.6); ax.spines[s].set_color("0.45")
ax.set_title("Relative make-up of the four programs", fontsize=7.6, fontweight="bold", pad=6, loc="left")

axb = fig.add_subplot(gs[0, 1])
for k, c in zip(KEY, CH):
    axb.plot(xs, per_stage[k].values, "-o", color=c, lw=1.4, ms=3.4, label=k)
axb.set_xticks(xs); axb.set_xticklabels(ST, fontsize=7)
axb.set_ylabel("Mean gene-set expression (CP10K, log1p)", fontsize=6.8)
axb.tick_params(labelsize=6.2, length=2.2)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axb.spines[s].set_linewidth(0.6); axb.spines[s].set_color("0.45")
axb.grid(axis="y", color="0.93", linewidth=0.5); axb.set_axisbelow(True)
axb.set_title("Absolute level per stage", fontsize=7.6, fontweight="bold", pad=6, loc="left")

hs, ls = ax.get_legend_handles_labels()
fig.legend(hs, ls, loc="center left", bbox_to_anchor=(0.875, 0.5), frameon=False,
           fontsize=6.6, title="Program", title_fontsize=6.8)
fig.savefig(f"{OUT}/P21_exhaustion_by_stage.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P21_exhaustion_by_stage.png", dpi=400, bbox_inches="tight")
per_stage.to_csv(f"{T}/exhaustion_by_stage.tsv", sep="\t")
print("✓ P21_exhaustion_by_stage")
