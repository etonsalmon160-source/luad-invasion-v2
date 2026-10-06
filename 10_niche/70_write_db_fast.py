#!/usr/bin/env python3
# 70_write_db_fast.py —— 直接用 h5py 写 signatureSearch 的库（替代 R 的逐列慢写）
# 结构（照 signatureSearch:::create_empty_h5）：
#   assay    (基因 × 签名) 数值
#   rownames (基因 × 1)   固定长字符串 40
#   colnames (签名 × 1)   固定长字符串 200
#   padj 非必需（gess_cor/gess_cmap 都不用），跳过 ⇒ 省 5.8 GB
import numpy as np, pandas as pd, h5py, time
from cmapPy.pandasGEXpress.parse_gctx import parse
OUT = "/home/eto/luad_v2/results/10_niche/sigsearch/lincs_db.h5"
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

log("读 gctx …")
G = parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df = G.data_df
genes = [str(g) for g in df.index]; sigs = [str(c) for c in df.columns]
si = pd.read_csv("/home/eto/lincs_data/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
                 sep="\t", compression="gzip").set_index("sig_id").reindex(sigs)
col = (si.pert_iname.fillna("NA").astype(str) + "__" +
       si.cell_id.fillna("NA").astype(str) + "__" +
       si.pert_type.fillna("NA").astype(str)).values
assert len(col) == len(sigs) and si.pert_iname.isna().sum() == 0
M = df.values.astype(np.float32)                 # 基因 × 签名
del df, G
log(f"矩阵 {M.shape}  ({M.nbytes/1e9:.1f} GB)；列名重复数 {len(col)-len(set(col)):,}")
log("写 H5 …")
with h5py.File(OUT, "w") as f:
    # 🔴 R(rhdf5) 按列序读、h5py 按 C 序写 ⇒ 必须写转置，R 才读到 (基因 × 签名)
    Mt = np.ascontiguousarray(M.T)
    f.create_dataset("assay", data=Mt, dtype="float32", chunks=(Mt.shape[0], 128))
    f.create_dataset("rownames", data=np.array(genes, dtype=h5py.string_dtype("utf-8", 40)))
    f.create_dataset("colnames", data=np.array(col, dtype=h5py.string_dtype("utf-8", 200)))
log(f"落盘 {OUT}  {pd.io.common.file_exists(OUT) and __import__('os').path.getsize(OUT)/1e9:.2f} GB")
