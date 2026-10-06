#!/usr/bin/env python3
# 41_sc_program_home.py —— 单细胞层面：D3 的 ECM 程序落在哪些 L2 型上（只读，探索性）
import h5py, numpy as np, pandas as pd
from scipy import sparse
from scipy.io import mmwrite
ROOT="/home/eto/luad_v2"; OUT=f"{ROOT}/results/10_niche/tr_singlecell"
f=h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5","r")
genes=np.array([g.decode() if isinstance(g,bytes) else str(g) for g in f["genes"][:]])
ct   =np.array([g.decode() if isinstance(g,bytes) else str(g) for g in f["cell_types"][:]])
M=sparse.csr_matrix((f["counts/data"][:], f["counts/indices"][:], f["counts/indptr"][:]),
                    shape=(len(ct), len(genes)))
print(f"参考：{M.shape[0]} 细胞 × {M.shape[1]} 基因；{len(set(ct))} 个 L2 型")

# CP10K + log1p（与项目口径一致）
tot=np.asarray(M.sum(axis=1)).ravel(); tot[tot==0]=1
Mn=sparse.diags(1e4/tot) @ M
Mn=Mn.tocsr(); Mn.data=np.log1p(Mn.data)

# 每型均表达
types=sorted(set(ct)); rows=[]
E=np.zeros((len(types), len(genes)))
for i,t in enumerate(types):
    E[i]=np.asarray(Mn[ct==t].mean(axis=0)).ravel()
E=pd.DataFrame(E, index=types, columns=genes)

# D3 深度配平后的 ECM 程序（来自 §9.2/残差版）
ECM=["COL1A1","COL1A2","COL3A1","COL6A3","FN1","CTHRC1","POSTN","SPARC","THBS2","BGN",
     "COL5A1","COL5A2","LUM","VCAN","MMP2","TIMP1","ASPH","COMP","THY1","COL8A1","COL10A1","COL11A1"]
ECM=[g for g in ECM if g in genes]
z=(E[ECM]-E[ECM].mean())/E[ECM].std(ddof=0)
E["ECM_score"]=z.mean(axis=1)
E["n_cells"]=pd.Series(ct).value_counts().reindex(types).values
E.sort_values("ECM_score",ascending=False)[["ECM_score","n_cells"]+ECM[:8]].head(15).to_csv(
    f"{OUT}/sc_ECM_score_by_L2.tsv",sep="\t")
print("\n=== ECM 程序得分最高的 15 个 L2 型 ===")
print(E.sort_values("ECM_score",ascending=False)[["ECM_score","n_cells"]+ECM[:6]].head(15).to_string())
print("\n=== ECM 程序得分最低的 6 个 ===")
print(E.sort_values("ECM_score")[["ECM_score","n_cells"]].head(6).to_string())
np.save(f"{OUT}/_E_mean.npy", E[[c for c in E.columns if c in genes]].values)
E[[c for c in E.columns if c in genes]].to_pickle(f"{OUT}/L2_mean_expr.pkl")
print(f"\n落盘 {OUT}/sc_ECM_score_by_L2.tsv")
