#!/usr/bin/env python3
# 42_sc_cellmap.py —— 复原"逐细胞 → L2 型 + 患者 + 分期"的映射（只读，探索性）
#
# 依据 08_spatial_deconv/00_build_rctd_reference.py:416/463-469
#   types 按 **h5ad 条码顺序**逐行构造；keep = (type is not None)
#   ⇒ reference_d.h5 的 cell_types[i] ↔ 第 i 个"被保留"的 h5ad 条码
# 校验：候选被剔集合是否恰使保留数 = 399,579，且 L2 型与 L1 谱系是否 100% 自洽。
import h5py, numpy as np, pandas as pd
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
ANN = f"{ROOT}/results/05_annotation"

# ① L2 型（参考对象里的顺序）
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as f:
    cell_types = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["cell_types"][:]])
print(f"参考 L2 型序列：{len(cell_types):,}")

# ② h5ad 条码顺序
with h5py.File(f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad", "r") as f:
    h5bc = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/_index"][:]])
print(f"h5ad 条码：{len(h5bc):,}")

# ③ 候选被剔集合
ex = pd.read_csv(f"{ANN}/gp8c_readopt_excluded_cells.csv.gz")
exset = set(ex.cell_barcode)
print(f"候选被剔：{len(exset):,}   ⇒ 保留 {len(h5bc)-len(set(h5bc)&exset):,}（需 = {len(cell_types):,}）")

# ④ 复原
keep_bc = np.array([b for b in h5bc if b not in exset])
assert len(keep_bc) == len(cell_types), f"长度不符 {len(keep_bc)} vs {len(cell_types)}"
D = pd.DataFrame({"cell_barcode": keep_bc, "L2": cell_types})

# ⑤ 校验：L2 → L1 谱系必须一对一
gp6 = pd.read_csv(f"{ANN}/gp6_cell_labels.csv.gz",
                  usecols=["cell_barcode", "patient_id", "sample_id", "stage", "A_frozen"])
D = D.merge(gp6, on="cell_barcode", how="left")
assert D.patient_id.notna().all(), "有条码对不上 gp6"
ct = pd.crosstab(D.L2, D.A_frozen)
bad = ct.apply(lambda r: r.sum() - r.max(), axis=1)
print("\n=== 校验：每个 L2 型是否只来自一个 L1 谱系 ===")
print(f"  不纯的 L2 型：{int((bad>0).sum())} / {len(bad)}（越大越可疑）")
if (bad > 0).any(): print(ct[bad > 0].to_string())
print(f"  ✅ 分期取值：{sorted(D.stage.unique())}")
print(f"  ✅ 患者数：{D.patient_id.nunique()}")

D.to_parquet(f"{OUT}/sc_cellmap.parquet", index=False)
print(f"\n落盘 {OUT}/sc_cellmap.parquet  ({len(D):,} 行)")

# ⑥ 各 L2 型的期别构成（比例，对深度免疫）
S = D.groupby(["L2", "stage"]).size().unstack(fill_value=0)
S = S.div(S.sum(axis=1), axis=0) * 100
order = ["Normal", "AAH", "AIS", "MIA", "IAC"]
S = S[[c for c in order if c in S.columns]]
print("\n=== 各 L2 型的期别构成（%）===")
print(S.round(1).to_string())
S.to_csv(f"{OUT}/L2_by_stage_pct.tsv", sep="\t")
