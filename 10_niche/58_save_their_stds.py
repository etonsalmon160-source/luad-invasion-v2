#!/usr/bin/env python3
# 58_save_their_stds.py —— 从 SCMG 官方流形算 gene_stds 并按扰动库 var 对齐存盘
import numpy as np, h5py, pandas as pd
from scipy import sparse
MANI="/home/eto/scmg_workspace/hf_data/ref_global_cell_state_manifold.h5ad"
PERT="/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad"
OUT="/home/eto/luad_v2/results/10_niche/tr_singlecell/scmg_official_gene_stds.tsv"
def cat(fh,p):
    o=fh[p]
    if isinstance(o,h5py.Group):
        c=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in o["categories"][:]]); return c[np.asarray(o["codes"][:])]
    return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in o[:]])
h=h5py.File(MANI,"r"); vars_=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in h["var/_index"][:]])
M=sparse.csr_matrix((h["X/data"][:],h["X/indices"][:],h["X/indptr"][:]),shape=(len(h["obs/_index"]),len(vars_))); h.close()
tot=np.asarray(M.sum(1)).ravel(); tot[tot==0]=1
Mn=(sparse.diags(1e4/tot)@M).tocsr(); Mn.data=np.log1p(Mn.data)
E1=np.asarray(Mn.mean(0)).ravel(); E2=np.asarray(Mn.power(2).mean(0)).ravel()
std=np.sqrt(np.maximum(E2-E1**2,0))
f=h5py.File(PERT,"r"); hid=cat(f,"var/human_id"); Gp=cat(f,"var/gene_name"); f.close()
m=dict(zip(vars_,std))
al=np.array([m.get(g,0.1) for g in hid])
al=np.maximum(al,0.1)
pd.DataFrame({"var_id":hid,"gene":Gp,"std":al}).to_csv(OUT,sep="\t",index=False)
print(f"存盘 {OUT}；std 范围 [{al.min():.3f}, {al.max():.3f}]；用真值的 {int((al>0.1).sum()):,}/{len(al):,}")
