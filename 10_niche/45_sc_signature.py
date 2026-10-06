#!/usr/bin/env python3
# 45_sc_signature.py —— 单细胞层面：细胞类型特异签名（患者内配对）
#
# 两个签名：
#   A_状态  ：成纤维内部 ECM-high vs ECM-low（患者内配对）
#   B_疾病  ：成纤维 IAC vs 前驱(AAH/AIS/MIA)（患者内配对）—— 与"生态位"最对应的疾病程序
# 输出给 CMap（46_），并报与空转 D3 签名的重叠。
import h5py, numpy as np, pandas as pd
from scipy import sparse
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"

D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as f:
    genes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["genes"][:]])
    M = sparse.csr_matrix((f["counts/data"][:], f["counts/indices"][:], f["counts/indptr"][:]),
                          shape=(len(D), len(genes)))
gi = {g: i for i, g in enumerate(genes)}
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = sparse.diags(1e4 / tot) @ M; Mn = Mn.tocsr(); Mn.data = np.log1p(Mn.data)

ECM = [g for g in ["COL1A1","COL1A2","COL3A1","COL6A3","FN1","CTHRC1","POSTN","SPARC","THBS2",
                   "BGN","COL5A1","COL5A2","LUM","VCAN","MMP2","TIMP1","ASPH","COL8A1"] if g in gi]
X = np.asarray(Mn[:, [gi[g] for g in ECM]].todense())
D["ECM"] = ((X - X.mean(0)) / (X.std(0) + 1e-9)).mean(1)

fib = D[D.L2 == "Fibroblast"].copy()
Mf = Mn[fib.index.values]                      # 🔴 必须用**原始 D 行号**取矩阵行；
fib = fib.reset_index(drop=True)               #    先取矩阵、后 reset_index（顺序反了就静默取错行）
fib["ECM_hi"] = fib.ECM > fib.ECM.quantile(0.75)
assert Mf.shape[0] == len(fib) == 35581, (Mf.shape, len(fib))

def paired_lfc(mask_a, mask_b):
    """逐患者 (a 均值 − b 均值)，只在两组都有的患者上算，再取患者中位。"""
    A = np.asarray(Mf[mask_a].mean(axis=0)).ravel()
    Bv = np.asarray(Mf[mask_b].mean(axis=0)).ravel()
    out = np.full(len(genes), np.nan); parts = []
    for pid, idx in fib.groupby("patient_id").groups.items():
        ia = idx if isinstance(idx, np.ndarray) else np.array(idx)
        ma = mask_a[ia]; mb = mask_b[ia]
        if ma.sum() >= 20 and mb.sum() >= 20:
            va = np.asarray(Mf[ia[ma]].mean(axis=0)).ravel()
            vb = np.asarray(Mf[ia[mb]].mean(axis=0)).ravel()
            parts.append(va - vb)
    if not parts: return None, 0
    P = np.vstack(parts)
    return np.median(P, axis=0), len(parts)

SIGS = {}
d, np_ = paired_lfc(fib.ECM_hi.values, (~fib.ECM_hi).values)
SIGS["A_state_ECMhi"] = d
print(f"A 状态签名（ECM-hi vs low）：用了 {np_} 例患者")
pre = fib.stage.isin(["AAH", "AIS", "MIA"]).values
iac = (fib.stage == "IAC").values
d2, np2 = paired_lfc(iac, pre)
SIGS["B_disease_IAC"] = d2
print(f"B 疾病签名（IAC vs 前驱）：用了 {np2} 例患者")

o = pd.read_csv(f"{ROOT}/results/10_niche/kstar_diag/d14_signed_panel.tsv", sep="\t")
o = o[o.archetype == 3]
for name, v in SIGS.items():
    r = pd.DataFrame({"gene": genes, "lfc": v}).dropna()
    r = r[np.isfinite(r.lfc)]
    r.sort_values("lfc", ascending=False).to_csv(f"{OUT}/sc_{name}_lfc.tsv", sep="\t", index=False)
    print(f"\n=== {name}：up 前 15 ===")
    print(", ".join(r.sort_values("lfc", ascending=False).head(15).gene))
    print(f"=== {name}：down 前 15 ===")
    print(", ".join(r.sort_values("lfc").head(15).gene))
    for dd in ["up", "down"]:
        og = set(o[o.direction == dd].sort_values("rank").gene.head(150))
        sg = set(r.sort_values("lfc", ascending=(dd == "down")).head(150).gene)
        print(f"  [与空转 D3 {dd} 侧重叠] {len(og & sg)}/150")

# 两签名互相关（应高度一致）
a = pd.DataFrame({"gene": genes, "A": SIGS["A_state_ECMhi"], "B": SIGS["B_disease_IAC"]}).dropna()
print(f"\n两签名 Spearman（全基因）：{a.A.corr(a.B, method='spearman'):.3f}")
