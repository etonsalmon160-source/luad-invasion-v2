#!/usr/bin/env python3
# 61_export_global_axis.py —— 把「全局状态轴」导出成签名文件，供化合物库 CMap 使用
import numpy as np, pandas as pd, h5py
from scipy import sparse
ROOT="/home/eto/luad_v2"; OUT=f"{ROOT}/results/10_niche/tr_singlecell"
D=pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5","r") as h:
    g=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in h["genes"][:]])
    M=sparse.csr_matrix((h["counts/data"][:],h["counts/indices"][:],h["counts/indptr"][:]),shape=(len(D),len(g)))
tot=np.asarray(M.sum(1)).ravel(); tot[tot==0]=1
Mn=(sparse.diags(1e4/tot)@M).tocsr(); Mn.data=np.log1p(Mn.data)
pid=D.patient_id.values; stg=D.stage.values
PRE=np.isin(stg,["AAH","AIS","MIA"]); IAC=stg=="IAC"
A=np.asarray(Mn[IAC].mean(0)).ravel()-np.asarray(Mn[PRE].mean(0)).ravel()
ps=[]
for p in np.unique(pid):
    a=(pid==p)&IAC; b=(pid==p)&PRE
    if a.sum()>=50 and b.sum()>=50: ps.append(np.asarray(Mn[a].mean(0)).ravel()-np.asarray(Mn[b].mean(0)).ravel())
B=np.median(np.vstack(ps),axis=0)
for nm,v in [("global_all",A),("global_paired",B)]:
    pd.DataFrame({"gene":g,"lfc":v}).to_csv(f"{OUT}/{nm}_lfc.tsv",sep="\t",index=False)
    print(nm,"up前8:",", ".join(pd.Series(g).iloc[np.argsort(-v)[:8]]))
    print(nm,"down前8:",", ".join(pd.Series(g).iloc[np.argsort(v)[:8]]))
