#!/usr/bin/env python3
# 31_tr_cmap_depthmatched.py —— 用**深度配平**的 D3 签名重跑 CMap，看 PCL 是否还站得住
#
# 🔴 关键对照：深度配平签名在**全部基因**上排，原签名在 **3,000 HVG** 上排 ⇒ 候选池不同会压低重叠。
#    故本脚本**在同一 3,000 HVG 内**重排，才是公平比较。
import numpy as np, pandas as pd, sys, os, time
sys.path.insert(0,"/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche")
from importlib.machinery import SourceFileLoader
TR=SourceFileLoader("tr","/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche/24_target_reversal.py").load_module()
OUT="/home/eto/luad_v2/results/10_niche/target_reversal"; KD="/home/eto/luad_v2/results/10_niche/kstar_diag"
DATA=f"{OUT}/data"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)

dm=pd.read_csv(f"{OUT}/d3_depthmatched_lfc.tsv",sep="\t",index_col=0).iloc[:,0]
hvg=set(pd.read_csv("/home/eto/luad_v2/results/10_niche/hvg_3000.txt",header=None)[0])
dm_h=dm[dm.index.isin(hvg)]
orig=pd.read_csv(f"{KD}/d14_signed_panel.tsv",sep="\t"); orig=orig[orig.archetype==3]
ou=orig[orig.direction=="up"].sort_values("rank").gene.head(150).tolist()
od=orig[orig.direction=="down"].sort_values("rank").gene.head(150).tolist()
nu=dm_h.sort_values(ascending=False).head(150).index.tolist()
nd=dm_h.sort_values().head(150).index.tolist()
log(f"【同一 3000 HVG 内重排】up ∩ = {len(set(ou)&set(nu))}/150 (J={len(set(ou)&set(nu))/len(set(ou)|set(nu)):.3f})"
    f" ; down ∩ = {len(set(od)&set(nd))}/150 (J={len(set(od)&set(nd))/len(set(od)|set(nd)):.3f})")

# ── CMap（与生产同口径）──
from cmapPy.pandasGEXpress.parse_gctx import parse
G=parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df=G.data_df
genes=[str(g) for g in df.index]
gi=TR.load_gene_map("GSE70138"); sym2e=dict(zip(gi.pr_gene_symbol.astype(str),gi.pr_gene_id.astype(str)))
gpos={g:i for i,g in enumerate(genes)}
up_ix=np.array([gpos[sym2e[g]] for g in nu if sym2e.get(g) in gpos])
dn_ix=np.array([gpos[sym2e[g]] for g in nd if sym2e.get(g) in gpos])
log(f"命中 Entrez：up {len(up_ix)} / down {len(dn_ix)}")
X=df.values.T.astype(np.float32); sig_ids=np.array([str(c) for c in df.columns]); N=X.shape[1]; del df,G
rank,sv,_=TR.build_rank_tables(X); del X
Ru,Cu=TR._rc(rank,sv,up_ix); e_u=TR.es_from_ranks(Ru,Cu,N,len(up_ix))
Rd,Cd=TR._rc(rank,sv,dn_ix); e_d=TR.es_from_ranks(Rd,Cd,N,len(dn_ix))
w=np.where(np.sign(e_u)!=np.sign(e_d),(e_u-e_d)/2.0,0.0)

si=TR.load_sig_info("GSE70138").set_index("sig_id").reindex(sig_ids)
sm=pd.read_csv(f"{DATA}/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",sep="\t",compression="gzip").set_index("sig_id").reindex(sig_ids)
pi2=pd.read_csv(f"{DATA}/lincs_pert_info2.tsv",sep="\t").set_index("pert_id")
lab=pd.DataFrame({"sig_id":sig_ids,"cell_id":si.cell_id.fillna("NA").values,"pert_type":si.pert_type.fillna("NA").values,
                  "pert_name":si.pert_iname.values,"pert_id":si.pert_id.values})
lab["is_ts"]=pi2.is_touchstone.reindex(lab.pert_id.values).values
lab["group"]=lab.cell_id+"|"+lab.pert_type
lab["distil_ss"]=sm.distil_ss.values; lab["distil_nsample"]=sm.distil_nsample.values
qref=set(lab[(lab.is_ts==1)&(lab.distil_nsample>=3)].sort_values(["distil_ss","sig_id"],ascending=[False,True])
         .drop_duplicates(["cell_id","pert_id"]).sig_id)
QM=np.zeros(len(sig_ids),bool); IDX={s:i for i,s in enumerate(sig_ids)}
for s in qref: QM[IDX[s]]=True
GRP=lab.group.values
ncs=np.empty_like(w)
for g in np.unique(GRP):
    m=GRP==g; ww=w[m]
    mp=ww[ww>0].mean() if (ww>0).any() else np.nan
    mn=-ww[ww<0].mean() if (ww<0).any() else np.nan
    ncs[m]=np.where(ww>0,ww/mp,np.where(ww<0,ww/mn,0.0))
tau=np.zeros_like(ncs)
for g in np.unique(GRP):
    m=GRP==g; ref=np.sort(np.abs(ncs[m&QM]))
    if ref.size==0: continue
    tau[m]=np.sign(ncs[m])*100.0/ref.size*np.searchsorted(ref,np.abs(ncs[m]),side="left")
lab["tau_dm"]=tau; lab["ncs_dm"]=ncs
lab.to_csv(f"{OUT}/L1a_D3_depthmatched_ranked.tsv",sep="\t",index=False)
log(f"深度配平签名 CMap：τ≤−90 {int((tau<=-90).sum())} 条（{100*np.mean(tau<=-90):.2f}%）"
    f" ; ≤−95 {int((tau<=-95).sum())} ; ≤−98 {int((tau<=-98).sum())}")
log(f"对照：原签名 τ≤−90 7790 条（6.60%）")
