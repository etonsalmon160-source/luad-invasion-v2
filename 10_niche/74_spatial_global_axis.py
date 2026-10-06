#!/usr/bin/env python3
# 74_spatial_global_axis.py —— 空转层面的重试：**深度配平的全局空间轴**
#
# 之前空转失败的原因：LUAD 的 spot 比前驱深 2.22 倍 ⇒ "IAC vs 前驱"的差异被深度吃掉。
# 这次的做法：**逐切片内按 nUMI 十分位分档**，只在**同档内**比 IAC vs 前驱，再按档合并。
#   ⇒ 深度在档内自动抵消，剩下的是"同等测序量下，IAC spot 与前驱 spot 的转录组差别"。
import numpy as np, pandas as pd, h5py, glob, os, time
from scipy import sparse
ROOT="/home/eto/luad_v2"; VIS=f"{ROOT}/data/visium_spatial"; OUT=f"{ROOT}/results/10_niche/sigsearch"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)
Q=pd.read_csv(f"{ROOT}/results/08_spatial_deconv/spatial_qc_per_spot.csv.gz")[["slide","barcode","nUMI"]]
slides=sorted(Q.slide.unique()); log(f"{len(slides)} 张切片")
NB=10
# ① 先建**并集**基因空间（不同切片面板不同：有的 18085、有的 10108）
feat_of={}
for s_ in slides:
    feat_of[s_]=pd.read_csv(f"{VIS}/{s_}/filtered_feature_bc_matrix/features.tsv.gz",sep="\t",header=None)[1].astype(str).values
genes=np.array(sorted(set().union(*[set(v) for v in feat_of.values()])))
gpos={g:i for i,g in enumerate(genes)}
log(f"并集基因空间 {len(genes):,} 个")

acc={b:{"a":np.zeros(len(genes)),"na":0,"b":np.zeros(len(genes)),"nb":0} for b in range(NB)}
import scipy.io as sio
for s in slides:
    st=s.split("_")[-1]
    stage="IAC" if st.startswith("LUAD") else ("PRE" if st in ("AAH","AIS","MIA") else None)
    if stage is None: continue
    fg=feat_of[s]; ix=np.array([gpos[g] for g in fg])
    Mx=sio.mmread(f"{VIS}/{s}/filtered_feature_bc_matrix/matrix.mtx.gz").tocsr()
    if Mx.shape[0]==len(fg) and Mx.shape[1]!=len(fg): Mn=Mx.T.tocsr()
    else: Mn=Mx.tocsr()
    tot=np.asarray(Mn.sum(1)).ravel()
    bnum=np.digitize(tot, np.quantile(tot, np.linspace(0,1,NB+1)[1:-1]))
    Mn=(sparse.diags(1e4/np.maximum(tot,1))@Mn).tocsr(); Mn.data=np.log1p(Mn.data)
    for b in range(NB):
        m=(bnum==b)
        if m.sum()==0: continue
        v=np.asarray(Mn[m].sum(0)).ravel()
        full=np.zeros(len(genes)); full[ix]=v
        if stage=="IAC": acc[b]["a"]+=full; acc[b]["na"]+=int(m.sum())
        else:            acc[b]["b"]+=full; acc[b]["nb"]+=int(m.sum())
    del Mx, Mn
    log(f"  {s} [{stage}] {Mn.shape if False else len(tot)} spot / 面板 {len(fg)}")
num=np.zeros(len(genes)); den=0.0
for b in range(NB):
    a,bb=acc[b]["a"],acc[b]["b"]; na,nb=acc[b]["na"],acc[b]["nb"]
    if na>=50 and nb>=50:
        num+=(a/na-bb/nb)*min(na,nb); den+=min(na,nb)
        log(f"  档{b}: IAC {na} / 前驱 {nb} ⇒ 入合并")
axis=num/den
pd.DataFrame({"gene":genes,"lfc":axis}).to_csv(f"{OUT}/spatial_global_depthmatched_lfc.tsv",sep="\t",index=False)
log(f"落盘；非零 {int((axis!=0).sum()):,}")
o=np.argsort(-axis); print("up 前10 :", ", ".join(genes[o[:10]])); print("down 前10:", ", ".join(genes[o[-10:]]))
