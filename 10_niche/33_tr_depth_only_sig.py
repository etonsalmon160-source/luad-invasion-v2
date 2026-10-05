#!/usr/bin/env python3
# 33_tr_depth_only_sig.py —— **决定性对照**：造一张"纯深度签名"（完全不碰 D3）跑 CMap
#   若它也出强信号 ⇒ 原 D3 签名的 CMap 信号是深度造的 ⇒ 必须改用深度配平签名
import numpy as np, pandas as pd, sys, os, time, glob
from scipy.io import mmread
import scipy.sparse as sp
sys.path.insert(0,"/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche")
from importlib.machinery import SourceFileLoader
TR=SourceFileLoader("tr","/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche/24_target_reversal.py").load_module()
ROOT="/home/eto/luad_v2"; OUT=f"{ROOT}/results/10_niche/target_reversal"; DATA=f"{OUT}/data"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)
Q=pd.read_csv(f"{ROOT}/results/08_spatial_deconv/spatial_qc_per_spot.csv.gz")
VIS=f"{ROOT}/data/visium_spatial"; slides=sorted(Q.slide.unique())
acc={}; nsm=0
for k,s in enumerate(slides):
    d=Q[Q.slide==s][["barcode","nUMI"]]
    if len(d)<500: continue
    feat=pd.read_csv(f"{VIS}/{s}/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header=None)
    bar=pd.read_csv(f"{VIS}/{s}/filtered_feature_bc_matrix/barcodes.tsv.gz",header=None)
    gi2={b:i for i,b in enumerate(bar[0])}
    ok=np.array([b in gi2 for b in d.barcode]); d=d[ok]
    idx=[gi2[b] for b in d.barcode]
    M=mmread(f"{VIS}/{s}/filtered_feature_bc_matrix/matrix.mtx.gz").tocsc()[:,idx]
    tot=np.asarray(M.sum(axis=0)).ravel(); tot[tot==0]=1
    Mn=(M@sp.diags(1e4/tot)).tocsr(); Mn.data=np.log1p(Mn.data)
    nu=d.nUMI.values
    hi=nu>=np.quantile(nu,0.8); lo=nu<=np.quantile(nu,0.2)   # 纯按深度：上/下五分位
    if hi.sum()<50 or lo.sum()<50: continue
    v=np.asarray(Mn[:,hi].mean(axis=1)).ravel()-np.asarray(Mn[:,lo].mean(axis=1)).ravel()
    for g,val in zip(feat[1].values,v): acc[g]=acc.get(g,0.0)+val
    nsm+=1
S=pd.Series(acc)/max(nsm,1); S.name="lfc_depthonly"
S.sort_values(ascending=False).to_csv(f"{OUT}/depth_only_lfc.tsv",sep="\t",header=True)
log(f"纯深度签名完成：{len(S)} 基因 / {nsm} 张切片")
hvg=set(pd.read_csv(f"{ROOT}/results/10_niche/hvg_3000.txt",header=None)[0]); Sh=S[S.index.isin(hvg)]
nu_=Sh.sort_values(ascending=False).head(150).index.tolist(); nd_=Sh.sort_values().head(150).index.tolist()

from cmapPy.pandasGEXpress.parse_gctx import parse
G=parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df=G.data_df
genes=[str(g) for g in df.index]
gi=TR.load_gene_map("GSE70138"); sym2e=dict(zip(gi.pr_gene_symbol.astype(str),gi.pr_gene_id.astype(str)))
gpos={g:i for i,g in enumerate(genes)}
up=np.array([gpos[sym2e[g]] for g in nu_ if sym2e.get(g) in gpos]); dn=np.array([gpos[sym2e[g]] for g in nd_ if sym2e.get(g) in gpos])
log(f"命中 Entrez up {len(up)} / down {len(dn)}")
X=df.values.T.astype(np.float32); sig_ids=np.array([str(c) for c in df.columns]); N=X.shape[1]; del df,G
rank,sv,_=TR.build_rank_tables(X); del X
Ru,Cu=TR._rc(rank,sv,up); e_u=TR.es_from_ranks(Ru,Cu,N,len(up))
Rd,Cd=TR._rc(rank,sv,dn); e_d=TR.es_from_ranks(Rd,Cd,N,len(dn))
w=np.where(np.sign(e_u)!=np.sign(e_d),(e_u-e_d)/2.0,0.0)
si=TR.load_sig_info("GSE70138").set_index("sig_id").reindex(sig_ids)
sm=pd.read_csv(f"{DATA}/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",sep="\t",compression="gzip").set_index("sig_id").reindex(sig_ids)
pi2=pd.read_csv(f"{DATA}/lincs_pert_info2.tsv",sep="\t").set_index("pert_id")
lab=pd.DataFrame({"sig_id":sig_ids,"pert_id":si.pert_id.values,"cell_id":si.cell_id.fillna("NA").values,
                  "pert_type":si.pert_type.fillna("NA").values,"pert_name":si.pert_iname.values})
lab["is_ts"]=pi2.is_touchstone.reindex(lab.pert_id.values).values
lab["group"]=lab.cell_id+"|"+lab.pert_type; lab["distil_ss"]=sm.distil_ss.values; lab["distil_nsample"]=sm.distil_nsample.values
qref=set(lab[(lab.is_ts==1)&(lab.distil_nsample>=3)].sort_values(["distil_ss","sig_id"],ascending=[False,True])
         .drop_duplicates(["cell_id","pert_id"]).sig_id)
QM=np.zeros(len(sig_ids),bool); IDX={s:i for i,s in enumerate(sig_ids)}
for s in qref: QM[IDX[s]]=True
GRP=lab.group.values
ncs=np.empty_like(w)
for g in np.unique(GRP):
    m=GRP==g; ww=w[m]; mp=ww[ww>0].mean() if (ww>0).any() else np.nan; mn=-ww[ww<0].mean() if (ww<0).any() else np.nan
    ncs[m]=np.where(ww>0,ww/mp,np.where(ww<0,ww/mn,0.0))
tau=np.zeros_like(ncs)
for g in np.unique(GRP):
    m=GRP==g; ref=np.sort(np.abs(ncs[m&QM]))
    if ref.size==0: continue
    tau[m]=np.sign(ncs[m])*100.0/ref.size*np.searchsorted(ref,np.abs(ncs[m]),side="left")
log("="*56)
log(f"【纯深度签名】τ≤−90 {int((tau<=-90).sum())} 条（{100*np.mean(tau<=-90):.2f}%）| ≤−95 {100*np.mean(tau<=-95):.2f}% | ≤−98 {100*np.mean(tau<=-98):.2f}%")
log(f"对照零模型：−90 4.70% | −95 2.42% | −98 1.02%")
log(f"对照原 D3 签名：−90 6.60% | −95 4.08% | −98 1.81%")
