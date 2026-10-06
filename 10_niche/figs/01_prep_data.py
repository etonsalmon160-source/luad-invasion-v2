#!/usr/bin/env python3
# 01_prep_data.py —— 为出图准备：UMAP × L2 标签 × 分期
import pandas as pd, numpy as np, gzip
ROOT="/home/eto/luad_v2"
U=pd.read_csv(f"{ROOT}/results/04_integration/seurat_trad/full/umap.csv.gz")
D=pd.read_parquet(f"{ROOT}/results/10_niche/tr_singlecell/sc_cellmap.parquet")
print("UMAP",U.shape,"| cellmap",D.shape)
M=U.merge(D,on="cell_barcode",how="inner")
print("合并后",M.shape)
# L1 谱系：由 L2 反推（用参考 manifest 的 39 型归属）
import h5py, json
man=json.load(open(f"{ROOT}/results/08_spatial_deconv/reference_d.manifest.json"))
print("L2 型:",len(man["cell_type_counts"]))
M.to_csv("/tmp/_cache_cellmap.csv",index=False)
print("缓存 /tmp/_cache_cellmap.csv")
print(M.stage.value_counts().to_dict())
print("列:",list(M.columns)[:12])
