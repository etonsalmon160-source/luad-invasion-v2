#!/usr/bin/env python3
# 34_tr_depth_residual_sig.py —— **修估计量**：逐切片对深度做三次回归、取残差再比域
#
# 为什么换：旧版把 spot 切十分位、每档算差、再取**中位** ⇒ 只用到十分位这一个尺度、
#   且中位丢弃信息 ⇒ **down 侧系统性失效**（六个域全 0–5/150，与域无关）。
# 新版：逐切片把**每个基因**对 log10(nUMI) 做三次回归取残差，再用**全部 spot** 比
#   D3 与非 D3 的残差均值 ⇒ 用满信息、深度已连续剔除。
#   一次线性代数算完所有基因：R = Y − X(XᵗX)⁻¹XᵗY
#
# 用法：TARGET=<1..7> python3 34_tr_depth_residual_sig.py
import numpy as np, pandas as pd, sys, os, time, glob
from scipy.io import mmread
import scipy.sparse as sp
ROOT="/home/eto/luad_v2"; OUT=f"{ROOT}/results/10_niche/target_reversal"
KD=f"{ROOT}/results/10_niche/kstar_diag"; VIS=f"{ROOT}/data/visium_spatial"
TARGET=int(os.environ.get("TARGET","3"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)

A=f"{KD}/d7_domain_assign.tsv"
D8=pd.concat([pd.read_csv(f,sep="\t") for f in glob.glob(f"{KD}/d8_parts/*.tsv")],ignore_index=True)
D8=D8.merge(pd.read_csv(A,sep="\t")[["slide","domain","archetype"]],on=["slide","domain"],how="left")
Q=pd.read_csv(f"{ROOT}/results/08_spatial_deconv/spatial_qc_per_spot.csv.gz")[["slide","barcode","nUMI"]]
D8=D8.merge(Q,on=["slide","barcode"],how="left")
log(f"D{TARGET}: spot {len(D8)}；该域 {int((D8.archetype==TARGET).sum())}")

acc={}; nsm=0
for k,s in enumerate(sorted(D8.slide.unique())):
    d=D8[D8.slide==s]
    if (d.archetype==TARGET).sum() < 100: continue
    feat=pd.read_csv(f"{VIS}/{s}/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header=None)
    bar=pd.read_csv(f"{VIS}/{s}/filtered_feature_bc_matrix/barcodes.tsv.gz",header=None)
    gi2={b:i for i,b in enumerate(bar[0])}
    ok=np.array([b in gi2 for b in d.barcode]); d=d[ok]
    idx=[gi2[b] for b in d.barcode]
    M=mmread(f"{VIS}/{s}/filtered_feature_bc_matrix/matrix.mtx.gz").tocsc()[:,idx]
    tot=np.asarray(M.sum(axis=0)).ravel(); tot[tot==0]=1
    Mn=(M@sp.diags(1e4/tot)).tocsr(); Mn.data=np.log1p(Mn.data)
    Y=np.asarray(Mn.todense(),dtype=np.float64).T          # spots × genes
    t=np.log10(d.nUMI.values+1.0); t=t-t.mean()
    X=np.column_stack([np.ones_like(t), t, t**2, t**3])    # 三次
    beta,_,_,_=np.linalg.lstsq(X,Y,rcond=None)
    R=Y-X@beta                                             # 残差 spots × genes
    is3=(d.archetype==TARGET).values
    dv=R[is3].mean(axis=0)-R[~is3].mean(axis=0)
    for g,val in zip(feat[1].values,dv): acc[g]=acc.get(g,0.0)+val
    nsm+=1
    if k%10==0: log(f"  {k+1}/{len(D8.slide.unique())} 切片；有效 {nsm}")
S=pd.Series(acc)/max(nsm,1); S.name="lfc_resid"
S.sort_values(ascending=False).to_csv(f"{OUT}/d{TARGET}_resid_lfc.tsv",sep="\t",header=True)
log(f"D{TARGET} 残差签名完成：{len(S)} 基因 / {nsm} 张切片")

hvg=set(pd.read_csv(f"{ROOT}/results/10_niche/hvg_3000.txt",header=None)[0]); Sh=S[S.index.isin(hvg)]
o=pd.read_csv(f"{KD}/d14_signed_panel.tsv",sep="\t"); o=o[o.archetype==TARGET]
for d_ in ["up","down"]:
    og=set(o[o.direction==d_].sort_values("rank").gene.head(150))
    ng=set(Sh.sort_values(ascending=(d_=="down")).head(150).index)
    log(f"  {d_}: 同一3000HVG ∩ = {len(og&ng)}/150 (J={len(og&ng)/len(og|ng):.3f})")
