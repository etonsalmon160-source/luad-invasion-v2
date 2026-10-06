#!/usr/bin/env python3
# 67_export_lincs_feather.py —— 把 LINCS 库导成 feather，供 signatureSearch 建库
# 列名格式：(pert_iname)__(cell_id)__(pert_type)，行名 = Entrez（gctx 的行就是 Entrez）
import numpy as np, pandas as pd, time
from cmapPy.pandasGEXpress.parse_gctx import parse
OUT = "/home/eto/luad_v2/results/10_niche/sigsearch"
import os; os.makedirs(OUT, exist_ok=True)
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

log("读 gctx …")
G = parse("/home/eto/lincs_data/GSE70138_Level5.gctx")
df = G.data_df                                          # 基因(Entrez) × 签名
log(f"  矩阵 {df.shape[0]:,} 基因 × {df.shape[1]:,} 签名")

si = pd.read_csv("/home/eto/lincs_data/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
                 sep="\t", compression="gzip").set_index("sig_id").reindex(df.columns)
miss = si.pert_iname.isna().sum()
log(f"  元数据对齐：缺 {miss} 条")
col = (si.pert_iname.fillna("NA").astype(str) + "__" +
       si.cell_id.fillna("NA").astype(str) + "__" +
       si.pert_type.fillna("NA").astype(str))
assert col.notna().all() and len(col) == df.shape[1]

M = df.values.astype(np.float32)                        # 基因 × 签名
log(f"  转置并写 feather（{M.nbytes/1e9:.1f} GB）…")
T = pd.DataFrame(M.T, columns=[str(g) for g in df.index], index=col.values)
T.index.name = "signature"
del M, df, G
import pyarrow as pa, pyarrow.feather as paf
paf.write_feather(pa.Table.from_pandas(T.reset_index(), preserve_index=False),
                  f"{OUT}/lincs_ref.feather", compression="uncompressed")   # 🔴 R 的 arrow 不支持 lz4
log(f"落盘 {OUT}/lincs_ref.feather")
print("列名样例:", list(T.columns[:3]), "| 行名样例:", list(T.index[:3]))
