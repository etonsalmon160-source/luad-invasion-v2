#!/usr/bin/env python3
# 30_tr_depth_matched_sig.py —— **决定性检验**：D3 的签名是不是深度造的？
#
# 问题：D3 域 spot 的 nUMI 中位 18,693 vs 其余域 3,946–5,642 ⇒ 原签名（域 vs 其余域）
#       可能混入"RNA 含量高"而非"缺氧侵袭"。
# 做法：**在同一切片内、按 nUMI 十分位配平**——每一档里只比同深度的 D3 spot 与 非D3 spot，
#       逐档算 log2FC，再按档取中位合并。
# 判据：新老签名重叠基因数；再用新表重跑 CMap，看逆转 PCL 是否仍是那批。
import numpy as np, pandas as pd, sys, os, time
ROOT="/home/eto/luad_v2"; import os as _os
OUT=f"{ROOT}/results/10_niche/target_reversal"
TARGET=int(_os.environ.get("TARGET","3"))
KD=f"{ROOT}/results/10_niche/kstar_diag"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)

# ① spot → D3 归属（d8）
import glob
A=pd.read_csv(f"{KD}/d8_parts/GSM9226168_P1_AAH.tsv",sep="\t")
d8=[]
for f in glob.glob(f"{KD}/d8_parts/*.tsv"):
    x=pd.read_csv(f,sep="\t"); d8.append(x)
D8=pd.concat(d8,ignore_index=True)
AR=pd.read_csv(f"{KD}/d7_domain_assign.tsv",sep="\t")[["slide","domain","archetype"]]
D8=D8.merge(AR,on=["slide","domain"],how="left")
log(f"spot 归属：{len(D8)} 条；D3 spot {int((D8.archetype==3).sum())}")

# ② nUMI
Q=pd.read_csv(f"{ROOT}/results/08_spatial_deconv/spatial_qc_per_spot.csv.gz")
D8=D8.merge(Q[["slide","barcode","nUMI"]],on=["slide","barcode"],how="left")
log(f"nUMI 缺失 {int(D8.nUMI.isna().sum())}")

# ③ 逐切片读表达，做十分位配平的 log2FC
from scipy.io import mmread
import gzip, scipy.sparse as sp
VIS=f"{ROOT}/data/visium_spatial"
slides=sorted(D8.slide.unique())
acc={}; nsm=0
for k,s in enumerate(slides):
    d=D8[D8.slide==s]
    if (d.archetype==TARGET).sum() < 100:
        continue
    feat=pd.read_csv(f"{VIS}/{s}/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header=None)
    bar=pd.read_csv(f"{VIS}/{s}/filtered_feature_bc_matrix/barcodes.tsv.gz",header=None)
    M=mmread(f"{VIS}/{s}/filtered_feature_bc_matrix/matrix.mtx.gz").tocsc()
    gi2={b:i for i,b in enumerate(bar[0])}
    ok=np.array([b in gi2 for b in d.barcode])      # 🔴 必须先过滤 d 再建 idx，否则长度不等
    d=d[ok]
    idx=[gi2[b] for b in d.barcode]
    M=M[:,idx]
    tot=np.asarray(M.sum(axis=0)).ravel(); tot[tot==0]=1
    Mn=M@sp.diags(1e4/tot)            # CP10K（🔴 M 是 基因×spot，要**右乘**缩小列）
    Mn=Mn.tocsr(); Mn.data=np.log1p(Mn.data)   # 🔴 全程稀疏，不转稠密（12k×14k 稠密要 1.4 GB/张）
    is3=(d.archetype==TARGET).values
    nu=d.nUMI.values
    dec=pd.qcut(nu,10,labels=False,duplicates="drop")
    rs=[]
    for q in np.unique(dec):
        m=dec==q
        a=is3&m; b=(~is3)&m
        if a.sum()<10 or b.sum()<10: continue
        rs.append(np.asarray(Mn[:,a].mean(axis=1)).ravel()-np.asarray(Mn[:,b].mean(axis=1)).ravel())
    if rs:
        v=np.median(np.vstack(rs),axis=0)
        for g,val in zip(feat[1].values,v):
            acc[g]=acc.get(g,0.0)+val
        nsm+=1
    if k%10==0: log(f"  {k+1}/{len(slides)} 切片；有效 {nsm}")
S=pd.Series(acc)/max(nsm,1)
S.name="lfc_dm"
S.sort_values(ascending=False).to_csv(f"{OUT}/d{TARGET}_depthmatched_lfc.tsv",sep="\t",header=True)
log(f"D{TARGET} 深度配平签名完成：{len(S)} 个基因，来自 {nsm} 张切片")

# ④ 与原始（未配平）D3 签名比
orig=pd.read_csv(f"{KD}/d14_signed_panel.tsv",sep="\t")
for d_ in ["up","down"]:
    o=set(orig[(orig.archetype==TARGET)&(orig.direction==d_)].sort_values("rank").gene.head(150))
    if d_=="up": n=set(S.sort_values(ascending=False).head(150).index)
    else:        n=set(S.sort_values().head(150).index)
    hv=set(pd.read_csv(f"{ROOT}/results/10_niche/hvg_3000.txt",header=None)[0])
    Sh=S[S.index.isin(hv)]
    nh=set(Sh.sort_values(ascending=(d_=="down")).head(150).index)
    log(f"  {d_}: 全基因池 ∩ = {len(o&n)}/150 (J={len(o&n)/len(o|n):.3f}) ; **同一3000HVG ∩ = {len(o&nh)}/150 (J={len(o&nh)/len(o|nh):.3f})**")
