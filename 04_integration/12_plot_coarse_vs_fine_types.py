"""粗粒度(r=0.3) vs 细粒度(r=0.7) 的上皮型清单 —— 判断"细粒度多出来的是什么"。

用户 2026-09-22 授权后跑的低分辨率实验的收尾图。左图只看型清单；
右图回答关键追问：细粒度多出来的 Goblet/Mucous，是**真的一个型**，
还是被并进 AT2 的混合群？（逐细胞查 MUC 与 ATF 两侧的 marker）

用法: python3 04_integration/12_plot_coarse_vs_fine_types.py
输出: results/05_annotation/figures/epiA_coarse_vs_fine_types.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata

ROOT = "/home/eto/luad_v2"
ANN = f"{ROOT}/results/05_annotation"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ANN}/figures"
os.makedirs(OUT, exist_ok=True)

fine = pd.read_csv(f"{ANN}/epiA_cluster_annotation.csv")
coarse = pd.read_csv(f"{ANN}/epiA_lowres_r030_cluster_annotation.csv")
types = sorted(set(fine["argmax"]) | set(coarse["argmax"]))
n_fine = [int((fine["argmax"] == t).sum()) for t in types]
n_coarse = [int((coarse["argmax"] == t).sum()) for t in types]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.5, 6.6))

# ---- ① 型清单：每个型占了多少个簇 --------------------------------------
y = np.arange(len(types))
h = 0.38
ax1.barh(y + h / 2, n_fine, h, color="tab:red", alpha=0.85, label="fine  r = 0.7  (27 clusters)")
ax1.barh(y - h / 2, n_coarse, h, color="tab:blue", alpha=0.85, label="coarse  r = 0.3  (17 clusters)")
for i, (a, b) in enumerate(zip(n_fine, n_coarse)):
    if a:
        ax1.text(a + 0.3, y[i] + h / 2, str(a), va="center", fontsize=10, color="tab:red")
    if b:
        ax1.text(b + 0.3, y[i] - h / 2, str(b), va="center", fontsize=10, color="tab:blue")
ax1.set_yticks(y)
ax1.set_yticklabels(types, fontsize=11)
ax1.set_xlabel("number of clusters called as this type", fontsize=12)
ax1.set_xlim(0, max(n_fine) * 1.22)
ax1.set_title("A. Type inventory barely changes when coarsening\n"
              "coarse keeps 5 of 6 types; the only loss is Goblet/Mucous (1 cluster, 3,587 nuclei)",
              fontsize=12)
ax1.legend(fontsize=10, loc="lower right")
ax1.grid(alpha=0.25, axis="x")

# ---- ② Goblet/Mucous 那一簇到底是什么 -----------------------------------
clu = pd.read_csv(f"{TRAD}/epiA_lowres/clusters.csv.gz",
                  usecols=["cell_barcode", "harmony_res0.7_seed0"])
is_gob = (clu["harmony_res0.7_seed0"] == 12).to_numpy()
ad = anndata.read_h5ad(f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad",
                       backed="r")
sub = ad[clu["cell_barcode"].tolist()].to_memory()
tot = np.asarray(sub.X.sum(1)).ravel()
tot[tot == 0] = 1
norm = sp.diags(1e4 / tot) @ sub.X

GENS = ["MUC5AC", "MUC5B", "SPDEF", "SFTPC", "SFTPA1"]
GENS = [g for g in GENS if g in sub.var_names]
M = norm[:, [sub.var_names.get_loc(g) for g in GENS]].toarray()
det_gob, mean_gob = (M[is_gob] > 0).mean(0) * 100, M[is_gob].mean(0)
det_oth, mean_oth = (M[~is_gob] > 0).mean(0) * 100, M[~is_gob].mean(0)

x = np.arange(len(GENS))
w = 0.38
ax2.bar(x - w / 2, mean_gob, w, color="tab:red", alpha=0.85,
        label=f"the 3,587 'Goblet/Mucous' nuclei  (r = 0.7 cluster 12)")
ax2.bar(x + w / 2, mean_oth, w, color="lightgrey", edgecolor="grey",
        label="all other epithelial nuclei")
for i in range(len(GENS)):
    ax2.text(x[i] - w / 2, mean_gob[i] + 6, f"{det_gob[i]:.0f}%", ha="center",
             fontsize=9.5, color="tab:red")
    ax2.text(x[i] + w / 2, mean_oth[i] + 6, f"{det_oth[i]:.0f}%", ha="center",
             fontsize=9.5, color="dimgrey")
ax2.set_xticks(x)
ax2.set_xticklabels(GENS, fontsize=11)
ax2.set_ylabel("mean expression  (CP10K;  % above bar = detection rate)", fontsize=11.5)
ax2.axvline(2.5, color="black", ls=":", lw=1.4)
ax2.text(1.0, ax2.get_ylim()[1] * 0.93, "mucous program", ha="center", fontsize=11)
ax2.text(3.5, ax2.get_ylim()[1] * 0.93, "surfactant / AT2 program", ha="center", fontsize=11)
ax2.set_title("B. What that extra cluster actually is\n"
              "MUC genes are 100x+ enriched there, but surfactant genes are still detected "
              "in the majority\nof those same nuclei  ->  a mixed AT2/mucous population, "
              "not a clean goblet type", fontsize=12)
ax2.legend(fontsize=9.5, loc="upper left")

fig.tight_layout()
p = f"{OUT}/epiA_coarse_vs_fine_types.png"
fig.savefig(p, dpi=140)
print(f"写出 {os.path.relpath(p, ROOT)}")
print("detection%: goblet-call", dict(zip(GENS, det_gob.round(1))))
print("detection%: other      ", dict(zip(GENS, det_oth.round(1))))
print("mean      : goblet-call", dict(zip(GENS, mean_gob.round(2))))
print("mean      : other      ", dict(zip(GENS, mean_oth.round(2))))
