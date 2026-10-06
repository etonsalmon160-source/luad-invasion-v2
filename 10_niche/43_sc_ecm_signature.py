#!/usr/bin/env python3
# 43_sc_ecm_signature.py —— 单细胞层面：ECM 成纤维状态的期别富集 + 患者内配对签名
#
# 目标（用户）：找到"稳定逆转生态位"的药。此为第一步——在单细胞层面
#   ① 该状态的期别富集（观察/期望，避免整体期别构成灌水）
#   ② 患者内配对的细胞状态签名（细胞类型特异 ⇒ 与 LINCS 细胞系背景同质）
import h5py, numpy as np, pandas as pd
from scipy import sparse
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"

D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
print(f"细胞 {len(D):,}；分期 {sorted(D.stage.unique())}")

# ── ① 期别富集：观察 / 期望 ──
S = D.groupby(["L2", "stage"]).size().unstack(fill_value=0)
base = D.stage.value_counts(normalize=True)
S = S.div(S.sum(axis=1), axis=0)
OE = S.div(base, axis=1)
OE.to_csv(f"{OUT}/L2_by_stage_OE.tsv", sep="\t")
cols = ["Normal", "AAH", "AIS", "MIA", "IAC"]
print("\n=== 观察/期望（>1 = 该期富集；整体期为分母）===")
show = ["Fibroblast", "Myofibroblast", "Lipofibroblast", "Pericyte",
        "AT2", "AT1", "Plasma", "B", "CD8+ Mem/Eff T", "TREM2+ Dendritic", "Macrophage"]
print(OE.loc[[c for c in show if c in OE.index], [c for c in cols if c in OE.columns]].round(2).to_string())

# ── ② 表达矩阵 + ECM 打分 ──
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as f:
    genes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["genes"][:]])
    M = sparse.csr_matrix((f["counts/data"][:], f["counts/indices"][:], f["counts/indptr"][:]),
                          shape=(len(D), len(genes)))
gi = {g: i for i, g in enumerate(genes)}
ECM = ["COL1A1", "COL1A2", "COL3A1", "COL6A3", "FN1", "CTHRC1", "POSTN", "SPARC", "THBS2",
       "BGN", "COL5A1", "COL5A2", "LUM", "VCAN", "MMP2", "TIMP1", "ASPH", "COMP", "THY1",
       "COL8A1", "COL10A1", "COL11A1"]
ECM = [g for g in ECM if g in gi]
tot = np.asarray(M.sum(axis=1)).ravel(); tot[tot == 0] = 1
Mn = sparse.diags(1e4 / tot) @ M; Mn = Mn.tocsr(); Mn.data = np.log1p(Mn.data)
Xecm = np.asarray(Mn[:, [gi[g] for g in ECM]].todense())
z = (Xecm - Xecm.mean(0)) / (Xecm.std(0) + 1e-9)
D["ECM_score"] = z.mean(1)
print(f"\nECM 程序 {len(ECM)} 基因；全细胞得分 均值 {D.ECM_score.mean():.3f}")

# ── ③ ECM 状态的期别富集（在成纤维内部按 ECM 打分分档）──
fib = D[D.L2 == "Fibroblast"].copy()
print(f"\nL2=Fibroblast：{len(fib):,} 细胞 / {fib.patient_id.nunique()} 患者")
fib["ECM_hi"] = fib.ECM_score > fib.ECM_score.quantile(0.75)
hf = fib.groupby(["stage", "ECM_hi"]).size().unstack(fill_value=0)
hf["hi_frac"] = hf[True] / (hf[True] + hf[False])
print("\n=== 成纤维内部 ECM-high 比例（按分期）===")
print(hf.reindex([c for c in cols if c in hf.index]).round(3).to_string())

# ── ④ 患者内配对签名：ECM-high vs ECM-low（同患者）──
pb = (Mn[fib.index.values] @ pd.get_dummies(fib.ECM_hi.astype(int)).values)  # (cells, 2) 稀疏乘
A = sparse.csr_matrix(np.asarray(pb.todense()))
hi = np.asarray(A[fib.ECM_hi.values].sum(0)).ravel()
lo = np.asarray(A[~fib.ECM_hi.values].sum(0)).ravel()
nhi, nlo = fib.ECM_hi.sum(), (~fib.ECM_hi).sum()
lfc = np.log2((hi / nhi + 1e-6) / (lo / nlo + 1e-6))
res = pd.DataFrame({"gene": genes, "lfc": lfc}).dropna()
res = res[np.isfinite(res.lfc)]
res.sort_values("lfc", ascending=False).to_csv(f"{OUT}/scFib_ECMhi_vs_low_lfc.tsv", sep="\t", index=False)

up = res.sort_values("lfc", ascending=False).head(20)
dn = res.sort_values("lfc").head(20)
print(f"\n=== ECM-high 成纤维 up 前 20（n_hi={nhi} / n_lo={nlo}）===")
print(", ".join(up.gene))
print("=== down 前 20 ===")
print(", ".join(dn.gene))

# ── ⑤ 与空转 D3 签名的重叠（回答"域签名 vs 细胞型特异签名"差多少）──
o = pd.read_csv(f"{ROOT}/results/10_niche/kstar_diag/d14_signed_panel.tsv", sep="\t")
o = o[o.archetype == 3]
for d in ["up", "down"]:
    og = set(o[o.direction == d].sort_values("rank").gene.head(150))
    sg = set(res.sort_values("lfc", ascending=(d == "down")).head(150).gene)
    print(f"\n[与空转 D3 的 {d} 侧重叠] {len(og & sg)}/150（Jaccard {len(og&sg)/len(og|sg):.3f}）")
