#!/usr/bin/env python3
# 49_fig_domain_reversal.py —— 状态空间之二：**哪个生态位域最容易被逆转**
#   七个域各有一条深度配平签名（target_reversal/d{1..7}_depthmatched_lfc.tsv，符号→Entrez）
#   与 LINCS 库（118,050 扰动 × 12,328 基因）逐行求相关 ⇒ 每个域得到一张化合物榜
#   口径与既有的 Cor 族一致（数值相关），只换查询轴；因此「域轴 vs 全局轴」可比
import h5py, numpy as np, pandas as pd, time, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; TR = f"{ROOT}/results/10_niche/target_reversal"
OUT = f"{ROOT}/results/paper_figures"; CACHE = "/tmp/_domain_reversal.npz"
DOMC = ["#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA", "#D95F8C", "#2E8B8B", "#6A3FA0"]
NAME = ["Airway", "iCAF", "ECM/interstitial", "Alveolar-cap.", "AT2", "Vascular", "Lymphoid"]

if os.path.exists(CACHE):
    z = np.load(CACHE, allow_pickle=True)
    COR, PERT, CAND = z["COR"], z["PERT"], list(z["CAND"])
    print("cache hit", COR.shape)
else:
    t0 = time.time()
    gi = pd.read_csv(f"{ROOT}/data/external/lincs_geneinfo.txt.gz", sep="\t") \
        if os.path.exists(f"{ROOT}/data/external/lincs_geneinfo.txt.gz") else None
    if gi is None:   # 引擎侧基因表
        gi = pd.read_csv("/home/eto/lincs_data/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
                         sep="\t")
    s2e = dict(zip(gi["pr_gene_symbol"].astype(str), gi["pr_gene_id"].astype(str)))
    with h5py.File(f"{ROOT}/results/10_niche/sigsearch/lincs_db.h5", "r") as f:
        PERT = np.array([x.decode() for x in f["colnames"][:]])
        RN = np.array([x.decode() for x in f["rownames"][:]])
        A = f["assay"][:, :].astype(np.float32)          # 118050 × 12328
    print(f"库 {A.shape}  载入 {time.time()-t0:.0f}s")
    e2i = {g: i for i, g in enumerate(RN)}
    Q = np.zeros((7, A.shape[1]), dtype=np.float32)
    for k in range(1, 8):
        d = pd.read_csv(f"{TR}/d{k}_depthmatched_lfc.tsv", sep="\t")
        d.columns = ["gene", "lfc"]
        ix = [(e2i[s2e[g]], v) for g, v in zip(d.gene, d.lfc) if g in s2e and s2e[g] in e2i]
        for i, v in ix:
            Q[k - 1, i] = v
        print(f"  D{k}: 命中 {len(ix)} 基因")
    Ac = A - A.mean(1, keepdims=True)
    An = np.linalg.norm(Ac, axis=1) + 1e-9
    Qc = Q - Q.mean(1, keepdims=True)
    Qn = np.linalg.norm(Qc, axis=1) + 1e-9
    COR = (Ac @ Qc.T) / (An[:, None] * Qn[None, :])      # 118050 × 7
    CAND = list(pd.read_csv("/tmp/_comp34.tsv", sep="\t", index_col=0).index)
    np.savez(CACHE, COR=COR, PERT=PERT, CAND=np.array(CAND))
    print(f"完成 {time.time()-t0:.0f}s")

pert_key = np.array([p.split("__")[0] for p in PERT])
MK = pd.DataFrame(COR)
MK["key"] = pert_key
MK = MK.groupby("key").mean()                             # 逐化合物跨细胞系取均值
print("化合物数:", len(MK))

CAND = [c for c in CAND if c in MK.index]
rows = []
for k in range(7):
    r = MK[k].rank()                                      # 越负越靠前 ⇒ 秩越小
    med_c = r.reindex(CAND).median()
    exp = len(MK) / 2
    rows.append(dict(dom=f"D{k+1}", dom_name=NAME[k], med_rank=med_c,
                     n_top100=int((r.reindex(CAND) <= 100).sum()), n_cand=len(CAND),
                     exp_top100=len(CAND) * 100 / len(MK)))
R = pd.DataFrame(rows)
print(R.round(1).to_string(index=False))

fig, ax = plt.subplots(figsize=(6.2, 3.4))
xs = np.arange(7)
ax.bar(xs, R.n_top100, width=0.62, color=DOMC, linewidth=0)
ax.axhline(R.exp_top100.iloc[0], color="0.45", ls=(0, (4, 2.5)), lw=1.0, zorder=3)
ax.text(6.45, R.exp_top100.iloc[0] * 1.10, f"chance ≈ {R.exp_top100.iloc[0]:.1f}",
        ha="right", va="bottom", fontsize=5.9, color="0.42")
for i, r in R.iterrows():
    ax.text(i, r.n_top100 + 0.35, str(int(r.n_top100)), ha="center", fontsize=6.6, color="0.25")
ax.set_xticks(xs); ax.set_xticklabels([f"D{i+1}\n{n}" for i, n in enumerate(NAME)],
                                      fontsize=6.2, linespacing=1.3)
ax.set_ylabel(f"of the 34 compounds, how many fall in the\n"
              f"top 100 when the query is that domain   (chance {R.exp_top100.iloc[0]:.1f})",
              fontsize=6.8)
ax.set_ylim(0, 30.5)
ax.tick_params(axis="y", labelsize=6.4, length=2.5); ax.tick_params(axis="x", length=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.6); ax.spines[s].set_color("0.45")
ax.set_title("Which niche domain is the reversible one?", fontsize=9,
             fontweight="bold", pad=16, loc="left")
fig.savefig(f"{OUT}/P19_domain_reversal.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P19_domain_reversal.png", dpi=400, bbox_inches="tight")
print("✓ P19_domain_reversal")
R.to_csv(f"{ROOT}/results/10_niche/domain_reversal_summary.tsv", sep="\t", index=False)
