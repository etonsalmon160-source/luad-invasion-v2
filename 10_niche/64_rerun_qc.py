#!/usr/bin/env python3
# 64_rerun_qc.py —— 用**质控后的全局签名**重跑官方 SCMG 管道
import sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
sys.path.insert(0,"/home/eto/scmg_workspace/SCMG")
from scmg.model.causal_prediction import CausalGenePredictor
import anndata as ad
OUT="/home/eto/luad_v2/results/10_niche/tr_singlecell"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)
def cat(fh,p):
    o=fh[p]
    if isinstance(o,h5py.Group):
        c=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in o["categories"][:]]); return c[np.asarray(o["codes"][:])]
    return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in o[:]])
f=h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad","r")
hid=cat(f,"var/human_id"); Gp=cat(f,"var/gene_name"); X=np.nan_to_num(np.asarray(f["X"][:],dtype=np.float32))
pe=cat(f,"obs/perturbed_gene"); pn=cat(f,"obs/perturbed_gene_name")
sg=np.asarray(f["obs/perturbation_sign"][:]); ds=cat(f,"obs/dataset"); f.close()
std=pd.read_csv(f"{OUT}/scmg_official_gene_stds.tsv",sep="\t")["std"].values.astype(np.float64)
A=ad.AnnData(X=X,obs=pd.DataFrame({"perturbed_gene":pe,"perturbed_gene_name":pn,"perturbation_sign":sg,"dataset":ds}),
             var=pd.DataFrame(index=pd.Index(hid)))
cgp=CausalGenePredictor(A,std)
gp={g:i for i,g in enumerate(Gp)}
for nm in ["global_paired","global_all"]:
    S=pd.read_csv(f"{OUT}/{nm}_qc.tsv",sep="\t")
    v=np.zeros(len(Gp))
    for g,l in zip(S.gene,S.lfc_qc):
        j=gp.get(g)
        if j is not None: v[j]=l
    log(f"{nm}: 非零 {int((v!=0).sum()):,}")
    r=cgp.calc_causal_scores(v)
    r=r.sort_values("causal_score",ascending=False).drop_duplicates("perturbed_gene",keep="first")
    r.to_csv(f"{OUT}/QC_{nm}_causal.tsv",sep="\t",index=False)
    print(f"\n{'='*76}\n【质控后 · {nm}】\n{'='*76}")
    print("① 疾病升高 → 敲低逆转")
    print(r[r.gene_shift_z>0].head(15)[["perturbed_gene_name","causal_score","gene_shift_z"]].round(4).to_string(index=False))
    print("\n② 疾病降低 → 过表达逆转")
    print(r[r.gene_shift_z<0].head(15)[["perturbed_gene_name","causal_score","gene_shift_z"]].round(4).to_string(index=False))
