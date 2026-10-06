#!/usr/bin/env python3
# 77_other_programs.py —— D-5a 补：**其他真实程序**是否同样强（ECM 的特异性终判）
import numpy as np, pandas as pd, h5py
from scipy import sparse
ROOT="/home/eto/luad_v2"; SC=f"{ROOT}/results/10_niche/tr_singlecell"
PROG={
 "ECM/胶原(我们)":["COL1A1","COL1A2","COL3A1","COL6A3","FN1","CTHRC1","POSTN","SPARC","THBS2","BGN","COL5A1","COL5A2","LUM","VCAN","MMP2","TIMP1","ASPH","COL8A1","COMP","THY1","COL10A1","COL11A1"],
 "细胞周期":["MKI67","TOP2A","CCNB1","CCNA2","CDK1","BUB1","AURKA","PLK1","CDC20","BIRC5","UBE2C","RRM2","TYMS","CCNB2","NDC80","NUSAP1","TPX2","ASPM","CENPF","KIF11"],
 "干扰素 I 型":["ISG15","IFI6","IFIT1","IFIT3","MX1","MX2","OAS1","OAS2","OASL","STAT1","STAT2","IRF7","IRF9","XAF1","IFI44L","IFI27","RSAD2","HERC5","USP18","EPSTI1"],
 "EMT":["VIM","CDH2","FN1","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","TGFB1","SPARC","COL1A1","ACTA2","MMP2","MMP9","ITGB1","CAV1","FBN1","TAGLN","SERPINE1","THBS1"],
 "缺氧":["VEGFA","LDHA","SLC2A1","PGK1","ENO1","HIF1A","BNIP3","PDK1","ADM","NDRG1","ANGPTL4","CA9","MIF","P4HA1","EGLN3","VHL","TFRC","HK2","ALDOA","GAPDH"],
 "凋亡":["BAX","BAK1","CASP3","CASP8","CASP9","BCL2L11","PMAIP1","FAS","TNFRSF10B","BBC3","APAF1","DIABLO","BID","CASP7","CYCS","TP53","PUMA" ],
 "氧化磷酸化":["NDUFA1","NDUFB1","SDHA","SDHB","UQCRC1","COX5A","COX7A2","ATP5F1A","ATP5F1B","ATP5MC1","CYCS","NDUFS1","NDUFV1","COX4I1","ATP5PF"],
 "上皮身份/AT2":["SFTPC","SFTPA1","SFTPA2","SFTPB","NAPSA","NKX2-1","ABCA3","CLDN18","EPCAM","KRT8","KRT18","AGER","HOPX","PDPN","CAV1"],
}
N_PERM=300; SEED=20261006
D=pd.read_parquet(f"{SC}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5","r") as h:
    g=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in h["genes"][:]])
    M=sparse.csr_matrix((h["counts/data"][:],h["counts/indices"][:],h["counts/indptr"][:]),shape=(len(D),len(g)))
tot=np.asarray(M.sum(1)).ravel(); tot[tot==0]=1
Mn=(sparse.diags(1e4/tot)@M).tocsr(); Mn.data=np.log1p(Mn.data)
pid=D.patient_id.values; l2=D.L2.values
PRE=np.isin(D.stage.values,["AAH","AIS","MIA"]); IAC=D.stage.values=="IAC"; fib=(l2=="Fibroblast")
mu=np.asarray(Mn[fib].mean(0)).ravel()
sd=np.sqrt(np.maximum(np.asarray(Mn[fib].power(2).mean(0)).ravel()-mu**2,1e-12))
ZT={}
for p in np.unique(pid):
    for tag,msk in [("I",(pid==p)&fib&IAC),("P",(pid==p)&fib&PRE)]:
        if msk.sum()>=20: ZT[(p,tag)]=(np.asarray(Mn[msk].mean(0)).ravel()-mu)/sd
VALID=[p for p in np.unique(pid) if (p,"I") in ZT and (p,"P") in ZT]
Z={p:(ZT[(p,"I")],ZT[(p,"P")]) for p in VALID}
gi={x:i for i,x in enumerate(g)}
del M, Mn
rng=np.random.default_rng(SEED)
# 零分布：按每程序的基因数配平（从成纤维表达池抽）
cand=np.where(mu>np.median(mu))[0]
print(f"患者 {len(VALID)} 例；成纤维表达池 {len(cand)}\n")
print(f"{'程序':16s} {'基因':>4s} {'例数':>5s} {'中位差':>8s} {'P(例数)':>9s} {'P(中位差)':>10s}")
for nm,lst in PROG.items():
    idx=np.array([gi[x] for x in lst if x in gi])
    if len(idx)<5: print(f"{nm:16s}  基因太少"); continue
    d=np.array([Z[p][0][idx].mean()-Z[p][1][idx].mean() for p in VALID])
    w=int((d>0).sum()); md=float(np.median(d))
    W=np.empty(N_PERM); Ds=np.empty(N_PERM)
    for k in range(N_PERM):
        ii=rng.choice(cand,size=len(idx),replace=False)
        dd=np.array([Z[p][0][ii].mean()-Z[p][1][ii].mean() for p in VALID])
        W[k]=(dd>0).sum(); Ds[k]=np.median(dd)
    print(f"{nm:16s} {len(idx):>4d} {w:>5d} {md:>8.4f} {(1+(W>=w).sum())/(1+N_PERM):>9.4f} {(1+(Ds>=md).sum())/(1+N_PERM):>10.4f}")
