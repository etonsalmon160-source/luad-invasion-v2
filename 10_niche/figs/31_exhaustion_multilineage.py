#!/usr/bin/env python3
# 31_exhaustion_multilineage.py —— 多谱系耗竭综合分（逐 spot，五个代表切片）
#   T/NK/B/髓系 四套耗竭 marker，各自 z 后取平均 ⇒ 一个"耗竭"分
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, json, glob, os
ROOT="/home/eto/luad_v2"; OUT="/home/eto/luad_v2/results/paper_figures"
REP=["GSM9226174_P4_Normal","GSM9226222_P25_AAH","GSM9226207_P19_AIS","GSM9226195_P13_MIA","GSM9226200_P15_LUAD"]
EXH={
 "T cell exhaustion":  ["PDCD1","LAG3","HAVCR2","TIGIT","ENTPD1","TOX","CXCL13","LAYN","TNFRSF9"],
 "NK exhaustion":      ["KLRC1","KIR2DL3","TIGIT","HAVCR2","PDCD1","CD160","KIR3DL1"],
 "B cell exhaustion":  ["FCRL4","ITGAX","ZEB2","FCRL2","CD86","FCRL5"],
 "Myeloid exhaustion": ["TREM2","APOE","GPNMB","SPP1","LPL","LGALS3","CD9"],
}
rows=[]
for sl in REP:
    d=f"{ROOT}/data/visium_spatial/{sl}/filtered_feature_bc_matrix"
    feat=pd.read_csv(f"{d}/features.tsv.gz",sep="\t",header=None)
    bar=pd.read_csv(f"{d}/barcodes.tsv.gz",header=None)[0].values
    M=sio.mmread(f"{d}/matrix.mtx.gz").tocsr()
    if M.shape[0]==len(feat) and M.shape[1]!=len(feat): M=M.T.tocsr()
    sym=feat[1].astype(str).values
    tot=np.asarray(M.sum(1)).ravel(); tot[tot==0]=1
    Mn=(sp.diags(1e4/tot)@M).tocsr(); Mn.data=np.log1p(Mn.data)
    gi={g:i for i,g in enumerate(sym)}
    axes=[]
    for name,mk in EXH.items():
        ix=[gi[g] for g in mk if g in gi]
        if len(ix)<3: print(f"  ! {sl} {name}: 只有 {len(ix)} 个 marker"); continue
        X=np.asarray(Mn[:,ix].todense()); Z=(X-X.mean(0))/(X.std(0)+1e-9)
        axes.append(Z.mean(1))
        rows.append({"slide":sl,"axis":name,"n_markers":len(ix)})
    comp=np.mean(np.vstack(axes),axis=0)
    pd.DataFrame({"slide":sl,"barcode":bar,"exhaust_score":comp}).to_csv(f"/tmp/_exh_{sl}.tsv",sep="\t",index=False)
    print(f"  {sl}: {comp.mean():+.4f} ± {comp.std():.4f}")
print("\n各轴用到的 marker 数:"); print(pd.DataFrame(rows).to_string(index=False))
