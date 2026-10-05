#!/usr/bin/env python3
# 28_tr_null_fullqref.py —— 靶点扰动逆向臂：**全量 Q_ref 下的零模型**
#
# 决定性问题：D3/D5/D7 的 τ 分布，与"随机查询"在同一条参照下相比，**有区别吗**？
#   若没有 ⇒ 候选表是噪声，不得作结论。
#   用同一 Q_ref（§5.14 ss 档 exemplar）、同一签名全集，只换 query。
import numpy as np, pandas as pd, sys, os, time
sys.path.insert(0,"/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche")
from importlib.machinery import SourceFileLoader
TR = SourceFileLoader("tr","/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche/24_target_reversal.py").load_module()

OUT="/home/eto/luad_v2/results/10_niche/target_reversal"; DATA=f"{OUT}/data"
N_RAND=int(sys.argv[1]) if len(sys.argv)>1 else 20
SEED=20261004
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)

# ── 与生产同口径地建库 ──
from cmapPy.pandasGEXpress.parse_gctx import parse
G=parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df=G.data_df
genes=[str(g) for g in df.index]; gi=TR.load_gene_map("GSE70138")
sym2e=dict(zip(gi.pr_gene_symbol.astype(str), gi.pr_gene_id.astype(str)))
gpos={g:i for i,g in enumerate(genes)}
X=df.values.T.astype(np.float32); sig_ids=np.array([str(c) for c in df.columns]); N=X.shape[1]
del df,G
log(f"库：{X.shape[0]} 签名 × {N} 基因")
rank,sv,_=TR.build_rank_tables(X); del X
si=TR.load_sig_info("GSE70138").set_index("sig_id").reindex(sig_ids)
sm=pd.read_csv(f"{DATA}/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",sep="\t",compression="gzip").set_index("sig_id").reindex(sig_ids)
pi2=pd.read_csv(f"{DATA}/lincs_pert_info2.tsv",sep="\t").set_index("pert_id")
lab=pd.DataFrame({"sig_id":sig_ids,"cell_id":si.cell_id.fillna("NA").values,
                  "pert_type":si.pert_type.fillna("NA").values,"pert_id":si.pert_id.values,
                  "pert_name":si.pert_iname.values})
lab["is_ts"]=pi2.is_touchstone.reindex(lab.pert_id.values).values
lab["group"]=lab.cell_id+"|"+lab.pert_type
lab["distil_ss"]=sm.distil_ss.values; lab["distil_nsample"]=sm.distil_nsample.values
qref=set(lab[(lab.is_ts==1)&(lab.distil_nsample>=3)].sort_values(["distil_ss","sig_id"],ascending=[False,True])
         .drop_duplicates(["cell_id","pert_id"]).sig_id)
qmask=np.zeros(len(sig_ids),dtype=bool); IDX={s:i for i,s in enumerate(sig_ids)}
for s in qref: qmask[IDX[s]]=True
GRP=lab.group.values
log(f"Q_ref {len(qref)} 条")

def ncs_of(w):
    o=np.empty_like(w)
    for g in np.unique(GRP):
        m=GRP==g; ww=w[m]
        mp=ww[ww>0].mean() if (ww>0).any() else np.nan
        mn=-ww[ww<0].mean() if (ww<0).any() else np.nan
        o[m]=np.where(ww>0,ww/mp,np.where(ww<0,ww/mn,0.0))
    return o
def tau_of(ncs):
    o=np.zeros_like(ncs)
    for g in np.unique(GRP):
        m=GRP==g; ref=np.sort(np.abs(ncs[m&qmask]))
        if ref.size==0: continue
        o[m]=np.sign(ncs[m])*100.0/ref.size*np.searchsorted(ref,np.abs(ncs[m]),side="left")
    return o

rng=np.random.default_rng(SEED); rows=[]
for k in range(N_RAND):
    g=rng.choice(N,size=300,replace=False)
    Ru,Cu=TR._rc(rank,sv,g[:150]); e_u=TR.es_from_ranks(Ru,Cu,N,150)
    Rd,Cd=TR._rc(rank,sv,g[150:]); e_d=TR.es_from_ranks(Rd,Cd,N,150)
    w=np.where(np.sign(e_u)!=np.sign(e_d),(e_u-e_d)/2.0,0.0)
    t=tau_of(ncs_of(w))
    rows.append({"rep":k,
                 "rate_le_m90":float((t<=-90).mean()), "rate_le_m95":float((t<=-95).mean()),
                 "rate_le_m98":float((t<=-98).mean()),
                 "rate_ge_p90":float((t>=90).mean()),
                 "median":float(np.median(t)),
                 "q01":float(np.quantile(t,0.01)), "q05":float(np.quantile(t,0.05))})
    if k%5==0: log(f"  null {k}/{N_RAND}  τ≤−90 {rows[-1]['rate_le_m90']:.4f}  ≤−95 {rows[-1]['rate_le_m95']:.4f}  ≤−98 {rows[-1]['rate_le_m98']:.4f}")
R=pd.DataFrame(rows); R.to_csv(f"{OUT}/null_fullqref.tsv",sep="\t",index=False)
obs={}
for q in ["D3","D5","D7"]:
    f=f"{OUT}/L1a_{q}/ranked.tsv"
    if os.path.exists(f):
        d=pd.read_csv(f,sep="\t")
        obs[q]={-90:float((d.tau_ss<=-90).mean()),-95:float((d.tau_ss<=-95).mean()),-98:float((d.tau_ss<=-98).mean())}
log("=" * 60)
for col,thr in [("rate_le_m90",-90),("rate_le_m95",-95),("rate_le_m98",-98)]:
    log("零模型 τ≤%d：均值 %.4f  范围 [%.4f, %.4f]" % (thr, R[col].mean(), R[col].min(), R[col].max()))
log("-"*60)
for q in obs:
    d=obs[q]; parts=[]
    for col,thr in [("rate_le_m90",-90),("rate_le_m95",-95),("rate_le_m98",-98)]:
        v=d[thr]; mu=R[col].mean(); sd=R[col].std()
        parts.append("τ≤%d: %.4f (z=%+.1f)" % (thr, v, (v-mu)/sd if sd>0 else 0))
    log("  观测 %s：%s" % (q, " | ".join(parts)))
log("判据：若观测落在零模型范围内（|z|<2）⇒ **无信号**，候选表不得作结论")
