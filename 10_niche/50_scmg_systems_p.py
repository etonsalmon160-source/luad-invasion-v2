#!/usr/bin/env python3
# 50_scmg_systems_p.py —— 给 49_ 的榜补**逐基因经验 p**（打乱 d 的完整零分布）
import os, time
os.environ.setdefault("OMP_NUM_THREADS","4")
import h5py, numpy as np, pandas as pd
from scipy import sparse
ROOT="/home/eto/luad_v2"; OUT=f"{ROOT}/results/10_niche/tr_singlecell"
SEED=20261005; N=300
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)
D=pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5","r") as f:
    sg=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in f["genes"][:]])
    M=sparse.csr_matrix((f["counts/data"][:],f["counts/indices"][:],f["counts/indptr"][:]),shape=(len(D),len(sg)))
tot=np.asarray(M.sum(1)).ravel(); tot[tot==0]=1
Mn=(sparse.diags(1e4/tot)@M).tocsr(); Mn.data=np.log1p(Mn.data)
pid=D.patient_id.values; l2=D.L2.values; PRE=np.isin(D.stage.values,["AAH","AIS","MIA"]); IAC=D.stage.values=="IAC"
def axis(m):
    ps=[]
    for p in np.unique(pid):
        ip=(pid==p)&m; a=ip&IAC; b=ip&PRE
        if a.sum()>=20 and b.sum()>=20: ps.append(np.asarray(Mn[a].mean(0)).ravel()-np.asarray(Mn[b].mean(0)).ravel())
    return np.median(np.vstack(ps),axis=0)
d_fib=axis(l2=="Fibroblast")
f=h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad","r")
def cat(p):
    g=f[p]; c=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in g["categories"][:]]); return c[np.asarray(g["codes"][:])]
G=cat("var/gene_name"); X=np.nan_to_num(np.asarray(f["X"][:],dtype=np.float32))
lab=pd.DataFrame({"gene":cat("obs/perturbed_gene_name"),"sign":np.asarray(f["obs/perturbation_sign"][:])}); f.close()
pos={g:i for i,g in enumerate(G)}
keep=np.array([i for i,g in enumerate(sg) if g in pos]); tgt=np.array([pos[sg[i]] for i in keep])
Xn=X[:,tgt].astype(np.float64); Xnorm=np.linalg.norm(Xn,axis=1)+1e-12
d0=d_fib[keep]-d_fib[keep].mean()
obs=-((Xn@d0)/(Xnorm*(np.linalg.norm(d0)+1e-12)))
rng=np.random.default_rng(SEED); null=np.empty((len(lab),N),dtype=np.float32)
for k in range(N):
    dd=d0.copy(); rng.shuffle(dd)
    null[:,k]=-((Xn@dd)/(Xnorm*(np.linalg.norm(dd)+1e-12))).astype(np.float32)
    if (k+1)%100==0: log(f"  null {k+1}/{N}")
lab["cos_fib"]=obs; lab["p"]=[(1+(null[i]>=obs[i]).sum())/(1+N) for i in range(len(obs))]
agg=lab.groupby(["gene","sign"]).agg(n=("cos_fib","size"),cos_fib=("cos_fib","median"),
                                     p=("p","median"),frac_pos=("cos_fib",lambda v:(v>0).mean())).reset_index()
agg=agg[agg.n>=3]
agg["方向"]=np.where(agg.sign<0,"敲低","过表达")
agg=agg.sort_values("cos_fib",ascending=False)
agg.to_csv(f"{OUT}/scmg_systems_ranked_p.tsv",sep="\t",index=False)
log(f"通过 p<0.05 的（基因×方向）：{int((agg.p<0.05).sum())} / {len(agg)}")
print("\n=== 最强 30（带 p）===")
print(agg.head(30)[["gene","方向","n","cos_fib","p","frac_pos"]].round(4).to_string(index=False))
