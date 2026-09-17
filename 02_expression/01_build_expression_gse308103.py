#!/usr/bin/env python3
"""
02_expression/01_build_expression_gse308103.py  —  Step 0

⚠️ **本脚本产出的是旧 648,945 × 18,082 对象 = 现在的「敏感性臂」，非现行分析口径。**
   现行口径（413,697 × 18,069，论文 QC）由 **`04_rebuild_expression_paperqc.py`** 产出。
   本脚本**未作废**（用于复现敏感性臂），但**请勿**误当作当前分析对象。
   登记：`PLAN_AND_CHECKPOINTS.md` §2.0 行 2 / §6「过门证据」；`PARAMETERS_AND_SOURCES.md` §M3-A.0。
   校验：`02_verify_expression_build.py`（**默认即校本对象**）；现行口径加 `GP0_MODE=paperqc`。

从 75 个**稠密文本**基因×细胞矩阵重建**稀疏 AnnData**。

为什么需要这步：M1 只落了汇总 CSV（per_cell_qc / per_sample），**没有持久化任何表达对象**。
M2（CopyKAT）与 M3-A（降维/聚类）都必须从原始计数开始。

铁律
  R1 分期无静默默认 —— patient/stage 一律经 cohort_registry.resolve_stage（未知 token raise）
  R3 无真实来源不计算 —— 只读真实测序矩阵；无 np.random / 无硬编码 / 无插补
  R5 产物可复现 + 哈希；patient_id 与 sample_id 分层

输入
  /home/eto/luad_invasion/data/GSE308103/extracted/<GSM>_<Pxx>_<Token>.raw_counts.mtx.txt.gz
      * 第 1 行 = barcode（tab 分隔，无行首标签）
      * 第 2.. 行 = <gene>\t<count>×n_cells   （**稠密整数**）
  results/01_qc/gse308103_per_cell_qc.csv.gz   （M1 掩码来源）
  results/00_ingest/paired_samples.csv         （权威 sample_id → patient/stage）

输出
  results/02_expression/per_sample/<sample_id>.npz   每样本 CSR（供 M2 CopyKAT 复用）
  results/02_expression/gse308103_counts.h5ad        全局稀疏（仅原始计数）
  results/02_expression/gse308103_analysis_mask.csv.gz
  results/02_expression/build_manifest.json

用法
  python 02_expression/01_build_expression_gse308103.py [--workers 8]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "00_ingest"))
import cohort_registry as REG  # noqa: E402

RAW_DIR = "/home/eto/luad_invasion/data/GSE308103/extracted"
QC_TABLE = os.path.join(ROOT, "results/01_qc/gse308103_per_cell_qc.csv.gz")
SAMPLES_CSV = os.path.join(ROOT, "results/00_ingest/paired_samples.csv")
OUT_DIR = os.path.join(ROOT, "results/02_expression")
PER_SAMPLE = os.path.join(OUT_DIR, "per_sample")

EXPECT_N_OBS = 648945           # M1 实测：qc_pass & singlet
EXPECT_N_VARS = 18082
FNAME_RE = re.compile(r"^(?P<gsm>GSM\d+)_(?P<sample_id>P\d+_[A-Za-z0-9]+)\.raw_counts\.mtx\.txt\.gz$")

log_lines = []


def log(msg):
    s = f"[{datetime.now().strftime('%H:%M:%S')}] {msg}"
    print(s, flush=True)
    log_lines.append(s)


def sha256_file(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# 1. 掩码：M1 的 (qc_pass & singlet) 逐样本 barcode 集合
# ---------------------------------------------------------------------------
def build_mask():
    qc = pd.read_csv(QC_TABLE,
                     usecols=["sample_id", "cell_barcode", "qc_pass", "doublet_class",
                              "nCount", "nFeature", "pct_mt", "doublet_score"])
    keep = qc["qc_pass"].astype(bool) & (qc["doublet_class"] == "singlet")
    sub = qc.loc[keep]
    total = len(sub)
    if total != EXPECT_N_OBS:
        raise SystemExit(f"掩码总数 {total} != 期望 {EXPECT_N_OBS} —— M1 表变了？拒绝继续。")

    # cell_barcode = <10x_barcode>|<sample_id>
    bad = ~sub["cell_barcode"].str.split("|", n=1).str[1].eq(sub["sample_id"])
    if bad.any():
        raise SystemExit(f"{int(bad.sum())} 行的 cell_barcode 后缀与 sample_id 不符 —— 拒绝继续。")
    sub = sub.assign(bc10x=sub["cell_barcode"].str.split("|", n=1).str[0])

    if sub["cell_barcode"].duplicated().any():
        raise SystemExit("cell_barcode 存在重复 —— 拒绝继续。")

    mask = {sid: set(g["bc10x"]) for sid, g in sub.groupby("sample_id")}
    log(f"掩码：{total} 细胞 / {len(mask)} 样本（qc_pass & singlet）")
    return mask, sub


# ---------------------------------------------------------------------------
# 2. 逐样本读取（worker）
# ---------------------------------------------------------------------------
def _read_one(args):
    path, sample_id, keep_bcs, gene_ref = args
    t0 = time.time()

    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        if not header:
            raise RuntimeError(f"{sample_id}: 空文件")
        bcs = header.rstrip("\n").split("\t")
        if len(bcs) != len(set(bcs)):
            raise RuntimeError(f"{sample_id}: barcode 有重复")

        dt = {i: np.int32 for i in range(1, len(bcs) + 1)}
        dt[0] = str
        df = pd.read_csv(fh, sep="\t", header=None, dtype=dt)
    n_txt_cells = df.shape[1] - 1
    if n_txt_cells != len(bcs):
        raise RuntimeError(f"{sample_id}: 表头 {len(bcs)} vs 数据列 {n_txt_cells}")
    genes = df.iloc[:, 0].to_numpy(dtype=object)
    n_genes_txt = len(genes)

    # 基因向量必须与参考逐元素一致（否则 hstack 会错位）
    if gene_ref is not None and n_genes_txt == len(gene_ref):
        if not np.array_equal(genes, gene_ref):
            raise RuntimeError(f"{sample_id}: 基因向量与参考不一致（{n_genes_txt}）")
    elif gene_ref is not None:
        raise RuntimeError(f"{sample_id}: 基因数 {n_genes_txt} != 参考 {len(gene_ref)}")

    # 按掩码取列
    col_of = {bc: i for i, bc in enumerate(bcs)}
    missing = [bc for bc in keep_bcs if bc not in col_of]
    if missing:
        raise RuntimeError(f"{sample_id}: {len(missing)} 个掩码 barcode 在该文件中缺失，例 {missing[:3]}")
    take = np.fromiter((col_of[bc] for bc in keep_bcs), dtype=np.int64, count=len(keep_bcs))

    dense = df.iloc[:, 1:].to_numpy(dtype=np.int32)[:, take]      # 18082 基因 × n_keep 细胞
    del df
    X = sp.csr_matrix(dense.T)                                    # → n_keep 细胞 × 18082 基因
    del dense
    X.sum_duplicates()
    X.eliminate_zeros()

    out = os.path.join(PER_SAMPLE, f"{sample_id}.npz")
    sp.save_npz(out, X.astype(np.float32))
    np.save(os.path.join(PER_SAMPLE, f"{sample_id}.genes.npy"), genes.astype("U"))

    return dict(
        sample_id=sample_id, gsm=FNAME_RE.match(os.path.basename(path)).group("gsm"),
        n_txt_cells=n_txt_cells, n_kept=int(X.shape[0]), n_genes=n_genes_txt,
        nnz=int(X.nnz), elapsed_s=round(time.time() - t0, 1),
        genes_sha256=hashlib.sha256("\n".join(map(str, genes)).encode()).hexdigest(),
    )


# ---------------------------------------------------------------------------
# 3. 主流程
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 4))
    ap.add_argument("--gene-ref", default=None,
                    help="参考基因向量 .npy（默认用样本表第一个文件）")
    args = ap.parse_args()

    os.makedirs(PER_SAMPLE, exist_ok=True)
    t_start = time.time()

    mask, sub = build_mask()

    # 权威 sample_id → (patient, stage)
    sm = pd.read_csv(SAMPLES_CSV)
    need_cols = {"sample_id", "patient_id", "stage", "stage_token"}
    if not need_cols.issubset(sm.columns):
        raise SystemExit(f"paired_samples.csv 缺列 {need_cols - set(sm.columns)}")
    sm = sm[sm["dataset"] == "GSE308103"] if "dataset" in sm.columns else sm
    meta = {r.sample_id: r for r in sm.itertuples()}

    # 文件清单：以掩码为准（不漏样本、不多样本）
    all_files = {}
    for fn in sorted(os.listdir(RAW_DIR)):
        m = FNAME_RE.match(fn)
        if m:
            all_files[m.group("sample_id")] = os.path.join(RAW_DIR, fn)
    log(f"原始文件 {len(all_files)} 个；掩码样本 {len(mask)} 个")

    only_mask = sorted(set(mask) - set(all_files))
    only_file = sorted(set(all_files) - set(mask))
    if only_mask:
        raise SystemExit(f"掩码里有但磁盘无文件：{only_mask}")
    if only_file:
        raise SystemExit(f"磁盘有文件但掩码里无（M1 未覆盖？）：{only_file}")
    miss_meta = sorted(set(mask) - set(meta))
    if miss_meta:
        raise SystemExit(f"掩码样本在权威样本表中缺失：{miss_meta}")

    # 参考基因向量
    if args.gene_ref:
        gene_ref = np.load(args.gene_ref, allow_pickle=False).astype("U")
    else:
        first = sorted(mask)[0]
        with gzip.open(all_files[first], "rt") as fh:
            fh.readline()
            gene_ref = pd.read_csv(fh, sep="\t", header=None, usecols=[0], dtype={0: str}).iloc[:, 0].to_numpy(dtype="U")
    log(f"参考基因向量 {len(gene_ref)}（首 {gene_ref[0]} … 末 {gene_ref[-1]}）")

    if len(gene_ref) != EXPECT_N_VARS:
        raise SystemExit(f"基因数 {len(gene_ref)} != 期望 {EXPECT_N_VARS}")

    # 并行抽取
    tasks = [(all_files[s], s, sorted(mask[s]), gene_ref) for s in sorted(mask)]
    log(f"并行抽取（workers={args.workers}）…")
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers) as pool:
        rows = []
        for i, r in enumerate(pool.imap_unordered(_read_one, tasks), 1):
            rows.append(r)
            if i % 10 == 0 or i == len(tasks):
                log(f"  {i}/{len(tasks)}  最近 {r['sample_id']} "
                    f"({r['n_kept']}/{r['n_txt_cells']} 保留, nnz={r['nnz']}, {r['elapsed_s']}s)")

    per = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)

    # 基因向量一致性
    if per["genes_sha256"].nunique() != 1:
        bad = per.loc[per["genes_sha256"] != per["genes_sha256"].mode()[0], ["sample_id", "genes_sha256"]]
        raise SystemExit(f"基因向量不一致（{len(bad)} 个样本）：\n{bad.to_string()}")
    log("✅ 全部 75 个文件的基因向量逐元素一致")

    if int(per["n_kept"].sum()) != EXPECT_N_OBS:
        raise SystemExit(f"保留细胞合计 {int(per['n_kept'].sum())} != {EXPECT_N_OBS}")

    # vstack → 全局 CSR（沿**细胞轴**拼接；各样本基因轴必须一致）
    order = per["sample_id"].tolist()
    blocks = [sp.load_npz(os.path.join(PER_SAMPLE, f"{s}.npz")).tocsr() for s in order]
    for b, s in zip(blocks, order):
        if b.shape[1] != EXPECT_N_VARS:
            raise SystemExit(f"{s}: 基因轴 {b.shape[1]} != {EXPECT_N_VARS}")
    X = sp.vstack(blocks, format="csr").astype(np.float32)
    del blocks
    if X.shape != (EXPECT_N_OBS, EXPECT_N_VARS):
        raise SystemExit(f"全局矩阵 shape {X.shape} != ({EXPECT_N_OBS}, {EXPECT_N_VARS})")
    log(f"全局稀疏矩阵 {X.shape}, nnz={X.nnz:,}, 稀疏度 {100*(1-X.nnz/(X.shape[0]*X.shape[1])):.2f}%")

    # obs：一律走 registry 的严格分期解析
    obs = sub.set_index("cell_barcode").loc[
        [f"{bc}|{s}" for s in order for bc in sorted(mask[s])]
    ].reset_index()
    if len(obs) != EXPECT_N_OBS:
        raise SystemExit("obs 行数与矩阵不符")
    obs["sample_id"] = pd.Categorical(obs["sample_id"])
    obs["patient_id"] = [meta[s].patient_id for s in obs["sample_id"]]
    obs["stage_token"] = [meta[s].stage_token for s in obs["sample_id"]]
    # R1：严格映射，未知 token 直接 raise（无静默默认）
    obs["stage"] = [REG.resolve_stage(t) for t in obs["stage_token"]]
    obs["patient_id"] = pd.Categorical(obs["patient_id"])
    obs["stage"] = pd.Categorical(obs["stage"], categories=REG.STAGES)
    obs = obs.set_index("cell_barcode")

    # mask CSV 必须在 AnnData 接管 obs **之前**固化：ad.AnnData 就地持有同一 DataFrame，
    # 而下面的 adata.obs_names.name = None 会把 obs.index.name 置为 None，
    # 之后 obs.reset_index() 给出的是 'index' 列而非 'cell_barcode'（已实测复现）。
    mask_df = obs.reset_index()[
        ["cell_barcode", "sample_id", "patient_id", "stage_token", "stage"]
    ].copy()

    import anndata as ad
    adata = ad.AnnData(X=X, obs=obs,
                       var=pd.DataFrame(index=pd.Index(gene_ref.astype(str), name=None)))
    adata.var_names.name = None
    adata.obs_names.name = None
    adata.uns["build"] = dict(
        script="02_expression/01_build_expression_gse308103.py",
        raw_dir=RAW_DIR, qc_table=QC_TABLE, source="M1 mask: qc_pass & doublet_class=='singlet'",
        built_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        note="仅原始计数；无归一化副本、无 .raw 槽（省内存）",
    )

    h5 = os.path.join(OUT_DIR, "gse308103_counts.h5ad")
    adata.write_h5ad(h5, compression="gzip")
    log(f"写出 {h5}  ({os.path.getsize(h5)/1e9:.2f} GB)")

    mask_csv = os.path.join(OUT_DIR, "gse308103_analysis_mask.csv.gz")
    if len(mask_df) != EXPECT_N_OBS:
        raise SystemExit(f"mask_df 行数 {len(mask_df)} != {EXPECT_N_OBS}")
    mask_df.to_csv(mask_csv, index=False, compression="gzip")

    per.to_csv(os.path.join(OUT_DIR, "per_sample_build_stats.csv"), index=False)

    manifest = dict(
        schema="luad_v2.expression_build/1",
        built_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        dataset="GSE308103", n_obs=int(X.shape[0]), n_vars=int(X.shape[1]), nnz=int(X.nnz),
        stage_counts={k: int(v) for k, v in obs["stage"].value_counts().items()},
        n_samples=int(len(order)),
        gene_vector_sha256=str(per["genes_sha256"].iloc[0]),
        gene_first=str(gene_ref[0]), gene_last=str(gene_ref[-1]),
        artifacts={
            "h5ad": dict(path=os.path.relpath(h5, ROOT), sha256=sha256_file(h5),
                         bytes=os.path.getsize(h5)),
            "mask_csv": dict(path=os.path.relpath(mask_csv, ROOT), sha256=sha256_file(mask_csv)),
            "per_sample_dir": os.path.relpath(PER_SAMPLE, ROOT),
        },
        script_sha256=sha256_file(os.path.abspath(__file__)),
        registry_sha256=sha256_file(os.path.join(ROOT, "00_ingest/cohort_registry.py")),
        runtime_s=round(time.time() - t_start, 1),
    )
    mf = os.path.join(OUT_DIR, "build_manifest.json")
    with open(mf, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)

    log(f"obs 分期分布: {manifest['stage_counts']}")
    log(f"清单 → {mf}")
    log(f"完成，用时 {manifest['runtime_s']}s")

    with open(os.path.join(ROOT, "logs/Step0_build_expression.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
