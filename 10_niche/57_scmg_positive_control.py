#!/usr/bin/env python3
# 57_scmg_positive_control.py —— **用他们的一切做阳性对照**
#
# ① gene_stds **从 SCMG 自己的参考流形算**（教程口径：CP10K+log1p 后逐基因 std）
# ② 阳性对照：**逐条复现他们教程发表的例子** Epiblast → Nascent mesoderm
#    期望：EOMES/TBXT/EVX1/SNAI1 为正，POU5F1/NANOG 为负
# ③ 若阳性对照成立 ⇒ 用**他们的 stds** 重跑我们的成纤维轴，看 0/2433 会不会变
import sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy import sparse
sys.path.insert(0, "/home/eto/scmg_workspace/SCMG")
from scmg.model.causal_prediction import CausalGenePredictor
import anndata as ad
ROOT = "/home/eto/luad_v2"; SC = f"{ROOT}/results/10_niche/tr_singlecell"
MANI = "/home/eto/scmg_workspace/hf_data/ref_global_cell_state_manifold.h5ad"
PERT = "/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

def cat(fh, p):
    o = fh[p]
    if isinstance(o, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]]); return c[np.asarray(o["codes"][:])]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o[:]])

# ── ① 他们的参考流形：算 gene_stds（教程口径）──
h = h5py.File(MANI, "r")
gt = cat(h, "var/human_gene_name") if "var/human_gene_name" in h["var"] else cat(h, "var/_index")
vars_ = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["var/_index"][:]])
ct = cat(h, "obs/cell_type")
M = sparse.csr_matrix((h["X/data"][:], h["X/indices"][:], h["X/indptr"][:]),
                      shape=(len(ct), len(vars_)))
h.close()
log(f"参考流形 {M.shape[0]:,} 细胞 × {M.shape[1]:,} 基因；cell_type {len(set(ct))} 类")
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
E1 = np.asarray(Mn.mean(0)).ravel(); E2 = np.asarray(Mn.power(2).mean(0)).ravel()
std_their = np.sqrt(np.maximum(E2 - E1 ** 2, 0))
log(f"gene_stds（他们的流形）算好：范围 [{std_their[std_their>0].min():.4f}, {std_their.max():.3f}]")

# ── ② 扰动库 ──
f = h5py.File(PERT, "r")
Gp = cat(f, "var/gene_name"); hid = cat(f, "var/human_id")
X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
pert_ens = cat(f, "obs/perturbed_gene"); pn = cat(f, "obs/perturbed_gene_name")
sign = np.asarray(f["obs/perturbation_sign"][:]); ds = cat(f, "obs/dataset"); f.close()
obs = pd.DataFrame({"perturbed_gene": pert_ens, "perturbed_gene_name": pn,
                    "perturbation_sign": sign, "dataset": ds})
A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(hid)))
# gene_stds 需按 var.index(=human_id/Ensembl) 对齐
ihid = {g: i for i, g in enumerate(hid)}
std_aligned = np.array([std_their[vars_.tolist().index(g)] if g in set(vars_) else 0.1 for g in hid])
log(f"gene_stds 对齐到扰动库 var：{int((std_aligned>0.1).sum()):,}/{len(hid):,} 有真值")
cgp = CausalGenePredictor(A, std_aligned)

# ── ③ 阳性对照：复现他们的 Epiblast → Nascent mesoderm ──
src_cells = np.isin(ct, ["Epiblast"])
tgt_cells = np.isin(ct, ["Nascent mesoderm", "Nascent Mesoderm"])
log(f"Epiblast {int(src_cells.sum())} 细胞；Nascent mesoderm {int(tgt_cells.sum())} 细胞")
shift = np.asarray(Mn[tgt_cells].mean(0)).ravel() - np.asarray(Mn[src_cells].mean(0)).ravel()
# shift 在流形基因空间 → 映射到扰动库 var 顺序
gp_pos = {g: i for i, g in enumerate(Gp)}
ihid_ = {g: i for i, g in enumerate(hid)}
ov_ens = len(set(vars_) & set(hid))          # 流形 var 是否能与库的 Ensembl 对上
ov_sym = len(set(vars_) & set(Gp))
mode = "Ensembl" if ov_ens > ov_sym else "symbol"
log(f"流形 var 与库的对齐：Ensembl 重合 {ov_ens:,} / symbol 重合 {ov_sym:,} ⇒ 用 {mode}")
sh_p = np.zeros(len(hid))
if mode == "Ensembl":
    for i, e in enumerate(vars_): sh_p[ihid_[e]] = shift[i]
else:
    for i, g in enumerate(vars_): sh_p[gp_pos[g]] = shift[i]
log(f"shift 非零 {int((sh_p!=0).sum()):,}")
r = cgp.calc_causal_scores(sh_p)
r = r.sort_values("causal_score", ascending=False).drop_duplicates("perturbed_gene", keep="first")
r.to_csv(f"{SC}/scmg_PC_gastrulation.tsv", sep="\t", index=False)
EXPECT_POS = ["EOMES", "TBXT", "EVX1", "SNAI1"]; EXPECT_NEG = ["POU5F1", "NANOG"]
rr = r.set_index("perturbed_gene_name")
print("\n=== 阳性对照：他们教程的期望 ===")
print(rr.reindex(EXPECT_POS + EXPECT_NEG)[["causal_score", "pert_sim", "gene_shift_z"]].round(4).to_string())
n_up = int((r.gene_shift_z > 0).sum()); n_dn = int((r.gene_shift_z < 0).sum())
print(f"\n全库中 gene_shift_z>0 的 {n_up:,} 个；EOMES 等的排名（z>0 内）：")
zpos = r[r.gene_shift_z > 0].reset_index(drop=True)
for g in EXPECT_POS:
    if g in set(zpos.perturbed_gene_name): print(f"  {g}: 第 {zpos.index[zpos.perturbed_gene_name==g][0]+1} / {len(zpos)}")
print("POU5F1/NANOG 在 z<0 里的排名：")
zneg = r[r.gene_shift_z < 0].reset_index(drop=True)
for g in EXPECT_NEG:
    if g in set(zneg.perturbed_gene_name): print(f"  {g}: 第 {zneg.index[zneg.perturbed_gene_name==g][0]+1} / {len(zneg)}")

# ── ④ 用他们的 stds 重跑我们的成纤维轴 ──
D = pd.read_parquet(f"{SC}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as hh:
    sg = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in hh["genes"][:]])
    M2 = sparse.csr_matrix((hh["counts/data"][:], hh["counts/indices"][:], hh["counts/indptr"][:]),
                           shape=(len(D), len(sg)))
t2 = np.asarray(M2.sum(1)).ravel(); t2[t2 == 0] = 1
Mn2 = (sparse.diags(1e4 / t2) @ M2).tocsr(); Mn2.data = np.log1p(Mn2.data)
pid = D.patient_id.values; l2 = D.L2.values
PRE = np.isin(D.stage.values, ["AAH", "AIS", "MIA"]); IAC = D.stage.values == "IAC"
fib = l2 == "Fibroblast"; ps = []
for p in np.unique(pid):
    a = (pid == p) & fib & IAC; b = (pid == p) & fib & PRE
    if a.sum() >= 20 and b.sum() >= 20:
        ps.append(np.asarray(Mn2[a].mean(0)).ravel() - np.asarray(Mn2[b].mean(0)).ravel())
axis = np.median(np.vstack(ps), axis=0)
ax_p = np.zeros(len(Gp))
for k, g in enumerate(sg):
    if g in gp_pos: ax_p[gp_pos[g]] = axis[k]
log(f"成纤维轴映射：非零 {int((ax_p!=0).sum()):,}")
r2 = cgp.calc_causal_scores(ax_p)
r2 = r2.sort_values("causal_score", ascending=False).drop_duplicates("perturbed_gene", keep="first")
r2.to_csv(f"{SC}/scmg_native_theirstds_fibroblast.tsv", sep="\t", index=False)
print("\n=== 用他们的 stds 跑我们的成纤维轴：top15（z>0）===")
print(r2[r2.gene_shift_z > 0].head(15)[["perturbed_gene_name", "causal_score", "gene_shift_z"]].round(4).to_string(index=False))
print("=== top15（z<0）===")
print(r2[r2.gene_shift_z < 0].head(15)[["perturbed_gene_name", "causal_score", "gene_shift_z"]].round(4).to_string(index=False))
