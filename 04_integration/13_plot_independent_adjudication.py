"""用**与我的 marker 面板不同源**的第三方注释（CellTypist / 人类肺图谱）来判两件事：

  ① 我判的每个上皮型，独立方法认不认？（纯度）
  ② 独立方法自己判出来的那些型，我的划分接不接得住？（召回）

这不是自证 —— 独立方法是**哪套分辨率都改变不了的外部锚**，故可用来定 r*。

用法: python3 04_integration/13_plot_independent_adjudication.py
输出: results/05_annotation/figures/independent_adjudication.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = "/home/eto/luad_v2"
ANN = f"{ROOT}/results/05_annotation"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ANN}/figures"
os.makedirs(OUT, exist_ok=True)

clu = pd.read_csv(f"{TRAD}/epiA_lowres/clusters.csv.gz")
lab = pd.read_csv(f"{ANN}/gp6_cell_labels_adjudicated.csv.gz",
                  usecols=["cell_barcode", "B_celltypist_type"])
ann7 = pd.read_csv(f"{ANN}/epiA_cluster_annotation.csv").set_index("cluster")["argmax"]
ann3 = pd.read_csv(f"{ANN}/epiA_lowres_r030_cluster_annotation.csv").set_index("cluster")["argmax"]
clu["ty7"] = clu["harmony_res0.7_seed0"].map(ann7)
clu["ty3"] = clu["harmony_res0.3_seed0"].map(ann3)
assert clu[["ty7", "ty3"]].notna().all().all(), "有簇映射不上"
m = clu.merge(lab, on="cell_barcode", how="left")

# 独立方法里"气道分泌/基底"这一族，以及被我合并掉的那些，统一压成少数几类便于阅读
GROUP = {
    "AT1": "AT1", "AT2": "AT2",
    "Multiciliated (non-nasal)": "Ciliated",
    "Deuterosomal": "Ciliated",
    "Basal resting": "Basal", "Suprabasal": "Basal",
    "Hillock-like": "Basal", "Basal activated": "Basal",
    "Club (non-nasal)": "Club", "Goblet (nasal)": "Goblet",
    "pre-TB secretory": "pre-TB secretory",
    "SMG mucous": "SMG mucous", "SMG serous (bronchial)": "SMG serous",
    "SMG duct": "SMG duct",
}
m["ct_grp"] = m["B_celltypist_type"].map(GROUP).fillna("other (non-epithelial / rare)")
ORDER = ["AT2", "AT1", "Ciliated", "Basal", "Club", "Goblet", "pre-TB secretory",
         "SMG mucous", "SMG serous", "SMG duct", "other (non-epithelial / rare)"]
COL = dict(zip(ORDER, ["#8c8c8c", "#4d4d4d", "#2f7fbf", "#8f4f9f", "#e08a1e",
                       "#c0392b", "#d98b8b", "#a0522d", "#7f9f3f", "#5f7f5f", "#dddddd"]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17.5, 7.4))

# ---- ① 纯度：我判的型，独立方法认不认 -----------------------------------
rows = [("ty3", t, 3) for t in ["AT1", "AT2", "Ciliated", "Basal", "Serous"]]
rows.append(("ty7", "Goblet/Mucous", 7))          # r=0.3 里已被并进 AT2，单独列出
labels, mats, ns = [], [], []
for col, t, _ in rows:
    s = m[m[col] == t]
    v = (s["ct_grp"].value_counts() / len(s) * 100)
    mats.append([v.get(k, 0.0) for k in ORDER])
    labels.append(f"{t}   ({len(s):,})")
    ns.append(len(s))
M = np.array(mats)
left = np.zeros(len(rows))
for j, k in enumerate(ORDER):
    ax1.barh(np.arange(len(rows)), M[:, j], left=left, color=COL[k],
             edgecolor="white", lw=0.6, label=k)
    for i in range(len(rows)):
        if M[i, j] >= 11:
            ax1.text(left[i] + M[i, j] / 2, i, f"{M[i, j]:.0f}", ha="center",
                     va="center", fontsize=8.6, color="white", fontweight="bold")
    left += M[:, j]
ax1.set_yticks(np.arange(len(rows)))
ax1.set_yticklabels(labels, fontsize=11)
ax1.invert_yaxis()
ax1.set_xlim(0, 100)
ax1.set_xlabel("how the INDEPENDENT cell atlas labels the nuclei in that type  (%)", fontsize=11.5)
ax1.set_title("A. Purity: does the outside method agree with each of my types?\n"
              "AT1 / AT2 / Ciliated are confirmed;  'Serous' is a mixture;\n"
              "the r=0.7 'Goblet/Mucous' is NOT confirmed (scatters into ciliated + club)",
              fontsize=11.5)
ax1.legend(fontsize=8, ncol=2, loc="lower left", framealpha=0.95)
ax1.grid(alpha=0.22, axis="x")

# ---- ② 召回：独立方法判出的型，我的划分接不接得住 ------------------------
ct_types = ["AT1", "AT2", "Basal", "Ciliated", "Club",
            "pre-TB secretory", "Goblet", "SMG mucous", "SMG serous"]
x = np.arange(len(ct_types))
w = 0.38
for off, col, lab_, color in [(-w / 2, "ty3", "my  r = 0.3", "tab:blue"),
                              (w / 2, "ty7", "my  r = 0.7", "tab:red")]:
    vals, capt = [], []
    for ct in ct_types:
        s = m[m["ct_grp"] == ct]
        if len(s) == 0:
            vals.append(0); capt.append(""); continue
        vc = s[col].value_counts()
        vals.append(vc.iloc[0] / len(s) * 100)
        capt.append(vc.index[0])
    b = ax2.bar(x + off, vals, w, color=color, alpha=0.85, label=lab_)
    for xi, v, c in zip(x + off, vals, capt):
        ax2.text(xi, v + 1.2, f"{v:.0f}", ha="center", fontsize=8.6, color=color,
                 fontweight="bold")
        ax2.text(xi, 2, c.replace("Goblet/Mucous", "Gob/Muc"), ha="center", fontsize=6.6,
                 color="white", rotation=90, va="bottom")
ax2.set_xticks(x)
ax2.set_xticklabels([t.replace(" (non-nasal)", "").replace(" (nasal)", "")
                     .replace(" (bronchial)", "") for t in ct_types],
                    fontsize=9.5, rotation=25, ha="right")
ax2.set_ylabel("share of that atlas type caught by my single largest cluster  (%)", fontsize=11)
ax2.set_ylim(0, 108)
ax2.axhline(90, color="grey", ls="--", lw=1.2)
ax2.text(len(ct_types) - 0.4, 91.5, "90%", color="grey", fontsize=9, ha="right")
ax2.set_title("B. Recall: the atlas' own cells, where do they land in my partition?\n"
              "big AT1/AT2/Ciliated/SMG-serous populations are caught (>=90%);\n"
              "the atlas' Club/Goblet cells are NOT (they fall into AT2)",
              fontsize=11.5)
ax2.legend(fontsize=10, loc="upper right")
ax2.grid(alpha=0.22, axis="y")

fig.tight_layout()
p = f"{OUT}/independent_adjudication.png"
fig.savefig(p, dpi=140)
print(f"写出 {os.path.relpath(p, ROOT)}")
print("\nA. 纯度（我判的型里，独立方法的构成 %）")
print(pd.DataFrame(M.round(1), index=labels, columns=ORDER).to_string())
print("\nB. 召回（独立方法的型，被我最大簇接住 %）")
for ct in ct_types:
    s = m[m["ct_grp"] == ct]
    a, b = s["ty7"].value_counts(), s["ty3"].value_counts()
    print(f"  {ct:<18}{len(s):>6} 核   r=0.7 {a.index[0]}({a.iloc[0]/len(s)*100:.0f}%)"
          f"   r=0.3 {b.index[0]}({b.iloc[0]/len(s)*100:.0f}%)")
