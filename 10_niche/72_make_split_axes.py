#!/usr/bin/env python3
# 72_make_split_axes.py —— 为稳定性检验生成"患者劈半"的查询轴
#   K 个随机劈半（11 vs 11）+ 同数目的**零模型**（逐患者内打乱分期标签）
import numpy as np, pandas as pd, h5py, os
from scipy import sparse
ROOT="/home/eto/luad_v2"; SC=f"{ROOT}/results/10_niche/tr_singlecell"; OUT=f"{ROOT}/results/10_niche/sigsearch/splits"
os.makedirs(OUT, exist_ok=True); K=4; SEED=20261006
D=pd.read_parquet(f"{SC}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5","r") as h:
    g=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in h["genes"][:]])
    M=sparse.csr_matrix((h["counts/data"][:],h["counts/indices"][:],h["counts/indptr"][:]),shape=(len(D),len(g)))
tot=np.asarray(M.sum(1)).ravel(); tot[tot==0]=1
Mn=(sparse.diags(1e4/tot)@M).tocsr(); Mn.data=np.log1p(Mn.data)
pid=D.patient_id.values; stg=D.stage.values.copy()
PAT=np.unique(pid); rng=np.random.default_rng(SEED)

def axis(sub_pat, labels=None):
    lab = stg if labels is None else labels
    PRE=np.isin(lab,["AAH","AIS","MIA"]); IAC=lab=="IAC"
    ps=[]
    for p in sub_pat:
        a=(pid==p)&IAC; b=(pid==p)&PRE
        if a.sum()>=50 and b.sum()>=50:
            ps.append(np.asarray(Mn[a].mean(0)).ravel()-np.asarray(Mn[b].mean(0)).ravel())
    return np.median(np.vstack(ps),axis=0) if len(ps)>=5 else None

rows=[]
for k in range(K):
    perm=rng.permutation(PAT); h1,h2=perm[:11],perm[11:]
    for tag,sub in [("h1",h1),("h2",h2)]:
        v=axis(sub)
        if v is not None: rows.append({"name":f"split{k}_{tag}","kind":"obs","lfc":v})
    # 零模型：逐患者内打乱分期标签
    lab=stg.copy()
    for p in PAT:
        i=np.where(pid==p)[0]; lab[i]=rng.permutation(lab[i])
    for tag,sub in [("h1",h1),("h2",h2)]:
        v=axis(sub,lab)
        if v is not None: rows.append({"name":f"split{k}_{tag}_null","kind":"null","lfc":v})
    print(f"split {k} 完成（{sum(1 for r in rows if r['name'].startswith(f'split{k}_'))} 条轴）")

E=pd.DataFrame({r["name"]:r["lfc"] for r in rows}); E.insert(0,"gene",g)
E.to_csv(f"{OUT}/split_axes.tsv",sep="\t",index=False)
pd.DataFrame([{"name":r["name"],"kind":r["kind"]} for r in rows]).to_csv(f"{OUT}/split_index.tsv",sep="\t",index=False)
print(f"\n落盘 {len(rows)} 条轴 → {OUT}/split_axes.tsv")
