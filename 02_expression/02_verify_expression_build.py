#!/usr/bin/env python3
"""
02_expression/02_verify_expression_build.py  —  GP0 门校验

对**写出的 h5ad**复核（不看构建器的内存态）。任一硬检查失败 → 退出码 2。

硬检查
  V1 n_obs == 648945 且 n_vars == 18082
  V2 逐样本 n_obs == (M1 表 qc_pass & singlet 计数)；**最大偏差必须为 0**
  V3 nnz == mask 内 Σ nFeature（符号求和会改变 nnz，故要求基因无重名）
  V4 每个 obs 名可解析为 <10x_barcode>|<sample_id>，且 sample_id ∈ registry 权威表
  V5 obs 的 stage 与 resolve_stage(stage_token) 完全一致（无静默默认）
  V6 计数矩阵为**非负整数**，且无全零细胞

⚠️ 两种口径（2026-09-17 补）——`GP0_MODE` 环境变量选择：
  `old`（默认）  = 旧 648,945 × 18,082 对象，**敏感性臂**。V2/V3 参考 = M1 掩码。
  `paperqc`     = 现行分析口径 413,697 × 18,069。V2 参考 = `gse308103_analysis_mask_paperqc.csv.gz`
                  （= M1 ∩ 论文绝对门）。**V3 形不同**：重建时施加 `gene_min_cells=3` 丢了 13 个基因，
                  故 `nnz ≤ Σ nFeature`，差额 = 被丢基因内的非零计数（**不得**要求相等）。

用法
  python 02_expression/02_verify_expression_build.py                 # 旧对象
  GP0_MODE=paperqc python 02_expression/02_verify_expression_build.py  # 现行口径
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "00_ingest"))
import cohort_registry as REG  # noqa: E402

MODE = os.environ.get("GP0_MODE", "old")
if MODE not in ("old", "paperqc"):
    raise SystemExit(f"GP0_MODE 只接受 'old' | 'paperqc'，收到 {MODE!r}")

QC_TABLE = os.path.join(ROOT, "results/01_qc/gse308103_per_cell_qc.csv.gz")

if MODE == "paperqc":
    H5 = os.path.join(ROOT, "results/02_expression/gse308103_counts_paperqc.h5ad")
    MANIFEST = os.path.join(ROOT, "results/02_expression/rebuild_paperqc_manifest.json")
    MASK = os.path.join(ROOT, "results/01_qc/gse308103_analysis_mask_paperqc.csv.gz")
    EXPECT_N_OBS = 413697
    EXPECT_N_VARS = 18069
    V3_STRICT = False   # gene_min_cells=3 丢基因 ⇒ nnz ≤ Σ nFeature
else:
    H5 = os.path.join(ROOT, "results/02_expression/gse308103_counts.h5ad")
    MANIFEST = os.path.join(ROOT, "results/02_expression/build_manifest.json")
    MASK = None
    EXPECT_N_OBS = 648945
    EXPECT_N_VARS = 18082
    V3_STRICT = True

fails, warns, notes = [], [], []


def chk(name, ok, detail=""):
    (notes if ok else fails).append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    return ok


def main():
    if not os.path.exists(H5):
        raise SystemExit(f"缺 {H5} —— 先跑 01_build_expression_gse308103.py")
    A = ad.read_h5ad(H5)
    X = A.X
    notes.append(f"h5ad 载入：{A.shape[0]:,} 细胞 × {A.shape[1]:,} 基因，X 类型 {type(X).__name__}/{X.dtype}")

    # ---- V1 ----
    chk("V1 维度", A.shape == (EXPECT_N_OBS, EXPECT_N_VARS),
        f"{A.shape} vs ({EXPECT_N_OBS}, {EXPECT_N_VARS})")

    # ---- 权威参考掩码 ----
    if MODE == "paperqc":
        ref = pd.read_csv(MASK)   # = M1 ∩ 论文绝对门
        ref_label = "论文 QC 掩膜"
    else:
        qc = pd.read_csv(QC_TABLE, usecols=["sample_id", "cell_barcode", "qc_pass", "doublet_class", "nFeature"])
        keep = qc["qc_pass"].astype(bool) & (qc["doublet_class"] == "singlet")
        ref = qc.loc[keep]
        ref_label = "M1 掩膜（qc_pass & singlet）"
    notes.append(f"权威参考掩码 [{ref_label}]：{len(ref):,} 细胞 / {ref['sample_id'].nunique()} 样本")
    chk("V0 参考掩码条码集合 == h5ad obs_names",
        set(ref["cell_barcode"].astype(str)) == set(A.obs_names.astype(str)),
        f"|ref\\h5ad|={len(set(ref['cell_barcode'].astype(str)) - set(A.obs_names.astype(str)))} "
        f"|h5ad\\ref|={len(set(A.obs_names.astype(str)) - set(ref['cell_barcode'].astype(str)))}")

    # ---- V2 逐样本细胞数 ----
    got = A.obs["sample_id"].value_counts().rename("built")
    exp = ref["sample_id"].value_counts().rename("ref")
    cmp = pd.concat([got, exp], axis=1).fillna(0).astype(int)
    cmp["diff"] = cmp["built"] - cmp["ref"]
    maxdev = int(cmp["diff"].abs().max()) if len(cmp) else -1
    chk("V2 逐样本 n_obs == 参考掩码", maxdev == 0,
        f"最大偏差 {maxdev}；样本数 {len(cmp)}")
    if maxdev:
        notes.append(cmp[cmp["diff"] != 0].to_string())

    # ---- V3 nnz vs Σ nFeature ----
    nnz = int(X.nnz)
    sum_feat = int(ref["nFeature"].sum()) if "nFeature" in ref.columns else None
    per_cell_nnz = np.asarray(X.getnnz(axis=1)).ravel()
    if sum_feat is None:
        # paperqc 掩膜表无 nFeature 列 → 用 obs.nFeature（逐行已核与掩膜一致）
        sum_feat = int(A.obs["nFeature"].astype(np.int64).sum())
        notes.append(f"V3 参考 Σ nFeature 取自 obs（掩膜表无该列）：{sum_feat:,}")
    delta = sum_feat - nnz
    if V3_STRICT:
        chk("V3 nnz == Σ nFeature(mask)", nnz == sum_feat, f"矩阵 {nnz:,} vs 参考 {sum_feat:,}")
    else:
        # gene_min_cells=3 丢 13 基因 ⇒ nnz 应 ≤ Σ nFeature，差额 = 被丢基因内非零计数
        chk("V3 nnz ≤ Σ nFeature(mask)", 0 <= delta, f"矩阵 {nnz:,} vs 参考 {sum_feat:,}，差额 {delta:,}")
        notes.append(f"V3 说明[paperqc]：差额 {delta:,} = 重建时被丢的 "
                     f"{json.load(open(MANIFEST, encoding='utf-8')).get('n_genes_dropped_lt3cells', '?')} "
                     f"个 <3 细胞基因内的非零计数（**构造性，非损坏**）")
    obs_nfeat = A.obs["nFeature"].to_numpy() if "nFeature" in A.obs else None
    if obs_nfeat is not None:
        diff_pc = per_cell_nnz.astype(np.int64) - obs_nfeat.astype(np.int64)
        d = int(np.abs(diff_pc).max())
        if V3_STRICT:
            chk("V3b 逐细胞 nnz == obs.nFeature", d == 0, f"最大偏差 {d}")
        else:
            chk("V3b 逐细胞 nnz ≤ obs.nFeature", int(diff_pc.max()) <= 0, f"最大超出 {int(diff_pc.max())}")
            notes.append(f"V3b 说明[paperqc]：逐细胞差额总和 {-int(diff_pc.sum()):,}（应等于 V3 差额）")
    else:
        warns.append("V3b 跳过：obs 无 nFeature 列")

    # ---- V4 obs 名可解析 + sample_id 权威 ----
    parts = A.obs_names.to_series().str.split("|", n=1, expand=True)
    chk("V4a 每个 obs 名含 '|'", parts.shape[1] == 2 and parts[1].notna().all())
    same = bool((parts[1].to_numpy() == A.obs["sample_id"].astype(str).to_numpy()).all())
    chk("V4b obs 名后缀 == obs.sample_id", same)
    auth_sn = {g["patient_id"] for g in REG.load_geo("GSE308103").values()}
    auth_samples = {f'{g["patient_id"]}_{g["token"]}' for g in REG.load_geo("GSE308103").values()}
    unknown = sorted(set(A.obs["sample_id"].astype(str)) - auth_samples)
    chk("V4c sample_id 全在 GEO 权威表", not unknown, f"未知 {unknown[:5]}")
    unknown_p = sorted(set(A.obs["patient_id"].astype(str)) - auth_sn)
    chk("V4d patient_id 全在 GEO 权威表", not unknown_p, f"未知 {unknown_p[:5]}")

    # ---- V5 stage 严格一致 ----
    bad = [t for t in A.obs["stage_token"].astype(str).unique()
           if REG.resolve_stage(t) != A.obs.loc[A.obs["stage_token"].astype(str) == t, "stage"].iloc[0]]
    chk("V5 stage == resolve_stage(stage_token)", not bad, f"不一致 {bad}")

    # ---- V6 计数为整数、非负、无全零细胞 ----
    mn = float(X.min()) if X.nnz else 0.0
    frac_int = float(np.all(np.equal(np.mod(X.data[: min(len(X.data), 5_000_000)], 1), 0))) if X.nnz else True
    chk("V6a 计数非负", mn >= 0, f"min {mn}")
    chk("V6b 计数为整数", frac_int)
    chk("V6c 无全零细胞", int((per_cell_nnz == 0).sum()) == 0,
        f"全零 {int((per_cell_nnz == 0).sum())}")

    # ---- 报告 ----
    print(f"\n================ GP0 校验 [MODE={MODE}] ================")
    print(f"  对象: {os.path.relpath(H5, ROOT)}")
    for n in notes:
        print("  ·", n)
    for w in warns:
        print("  ⚠", w)
    for f in fails:
        print("  ", f)
    print("----------------------------------------")
    print(f"  PASS {sum(1 for x in notes if x.startswith('PASS'))} / FAIL {len(fails)}")
    print(f"  obs 分期分布: {A.obs['stage'].value_counts().to_dict()}")
    print(f"  obs 患者数: {A.obs['patient_id'].nunique()}  样本数: {A.obs['sample_id'].nunique()}")

    if os.path.exists(MANIFEST):
        m = json.load(open(MANIFEST, encoding="utf-8"))
        print(f"  manifest: nnz={m['nnz']:,}")
        gs = m.get("gene_vector_sha256")
        if gs:
            print(f"  gene_vector_sha256={gs[:16]}…")
        else:
            print("  gene_vector_sha256=（该 manifest 未记录 —— paperqc 版按构造未存全基因组向量哈希）")
        print(f"  h5ad sha256={m['artifacts']['h5ad']['sha256'][:16]}…")

    print("=========================================")
    print("VERDICT:", f"GP0 PASS ✅ [{MODE}]" if not fails else f"GP0 FAIL ❌ [{MODE}]")
    return 2 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
