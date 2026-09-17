#!/usr/bin/env python3
"""
02_expression/03_audit_expression_build.py  —  GP0 对抗性审计

⚠️ **本脚本审计的是旧 648,945 × 18,082 对象（＝敏感性臂），非现行分析口径。**
   现行口径（413,697 × 18,069）的校验见 `02_verify_expression_build.py GP0_MODE=paperqc`；
   其独立重抽取类审计（对应本脚本的 A5）**尚未对 paperqc 对象做**，如实标缺。
   登记：`PLAN_AND_CHECKPOINTS.md` §6「过门证据」GP0 块。

**这不是 GP0 门校验的替代品，而是它的补强。**

`02_verify_expression_build.py` 的 V1–V6 证明的是**自洽性**：
矩阵与 M1 的表对得上。但若错误在**源头**就存在（阶段挂错、基因轴错位、
数字放错行列），自洽性检查抓不住——两方同源，会一起错。

本脚本专门查"同源检查抓不住"的部分：

  A1  基因轴：h5ad var_names 与**原始文本**的基因列**逐元素**一致；无重名、无空名
  A2  细胞轴：**逐细胞**行和 == M1 nCount（最大偏差必须为 0）
      —— 与 V3b 的 nnz==nFeature **相互独立**（一个数非零个数，一个数总和），
         两者同时精确相等几乎不可能由"行列放错"碰巧满足
  A3  逐样本行和 == M1 逐样本 Σ nCount
  A4  全局：X.data 全体有限、非负、**全为整数**（V6b 只抽查了前 500 万个元素）
  A5  **独立重抽取**：用与构建器**不同的代码路径**从原始 .txt.gz 重算 3 个样本的
      矩阵，与 h5ad 对应块**逐元素**比对。这是唯一能证伪"数字放错位置"的检查。
  A6  **R1 闭合**：每个样本的 stage_token / patient_id 与 GEO 权威表 `GSE308103_samples.tsv`
      逐行一致；且 obs.stage == resolve_stage(权威 token)
      —— V5 只查自洽（stage 与 token 内部一致），**不查 token 本身是否属于该样本**
  A7  每样本 npz ↔ h5ad 对应块：shape / nnz / **逐元素**一致
  A8  每样本 原始文件 barcode 集合 ⊇ M1 掩码集合；h5ad obs barcode 集合 == 掩码集合

用法
  python 02_expression/03_audit_expression_build.py [--n-reextract 3]
退出码 0 = 全部通过；2 = 有 FAIL
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "00_ingest"))
import cohort_registry as REG  # noqa: E402

RAW_DIR = "/home/eto/luad_invasion/data/GSE308103/extracted"
QC_TABLE = os.path.join(ROOT, "results/01_qc/gse308103_per_cell_qc.csv.gz")
OUT_DIR = os.path.join(ROOT, "results/02_expression")
PER_SAMPLE = os.path.join(OUT_DIR, "per_sample")
H5 = os.path.join(OUT_DIR, "gse308103_counts.h5ad")
GEO_TSV = os.path.join(ROOT, "00_ingest/geo_metadata/GSE308103_samples.tsv")

EXPECT_N_OBS = 648945
EXPECT_N_VARS = 18082
FNAME_RE = re.compile(r"^(?P<gsm>GSM\d+)_(?P<sample_id>P\d+_[A-Za-z0-9]+)\.raw_counts\.mtx\.txt\.gz$")

fails, passes, notes = [], [], []


def chk(name, ok, detail=""):
    (passes if ok else fails).append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-reextract", type=int, default=3,
                    help="独立重抽取的样本数（0 = 跳过，仅当时间受限）")
    args = ap.parse_args()

    # ---------------------------------------------------------------- 载入
    m1 = pd.read_csv(QC_TABLE, usecols=["sample_id", "cell_barcode", "qc_pass",
                                        "doublet_class", "nCount", "nFeature"])
    keep = m1["qc_pass"].astype(bool) & (m1["doublet_class"] == "singlet")
    m1 = m1.loc[keep].copy()
    m1["bc10x"] = m1["cell_barcode"].str.split("|", n=1).str[0]
    notes.append(f"M1 掩码 {len(m1):,} 细胞 / {m1['sample_id'].nunique()} 样本")
    if len(m1) != EXPECT_N_OBS:
        raise SystemExit(f"M1 掩码数 {len(m1)} != {EXPECT_N_OBS}")

    A = ad.read_h5ad(H5)
    X = A.X.tocsr()
    notes.append(f"h5ad {A.shape[0]:,} × {A.shape[1]:,}，nnz {X.nnz:,}")

    raw_files = {}
    for fn in sorted(os.listdir(RAW_DIR)):
        m = FNAME_RE.match(fn)
        if m:
            raw_files[m.group("sample_id")] = os.path.join(RAW_DIR, fn)

    # ================================================================ A1 基因轴
    ref_sid = sorted(raw_files)[0]
    with gzip.open(raw_files[ref_sid], "rt") as fh:
        fh.readline()
        genes_raw = pd.read_csv(fh, sep="\t", header=None, usecols=[0],
                                dtype={0: str}).iloc[:, 0].to_numpy(dtype=object)
    vn = A.var_names.to_numpy()
    chk("A1a var_names 长度 == 原始基因列", len(vn) == len(genes_raw),
        f"{len(vn)} vs {len(genes_raw)}")
    eq = bool(len(vn) == len(genes_raw) and np.array_equal(vn.astype(str), genes_raw.astype(str)))
    chk("A1b var_names 与原始文本基因列**逐元素**一致", eq,
        f"{ref_sid}；首 {vn[0]} / 末 {vn[-1]}")
    chk("A1c 基因无重名", A.var_names.is_unique,
        f"重复 {int(A.var_names.duplicated().sum())} 个")
    bad_names = [g for g in vn if (not isinstance(g, str)) or g.strip() == "" or g == "nan"]
    chk("A1d 基因无空名/NaN", not bad_names, f"异常 {bad_names[:5]}")

    # ================================================================ A2/A3 细胞轴
    m1i = m1.set_index("cell_barcode")
    aligned = m1i.reindex(A.obs_names.to_numpy())
    chk("A2a obs_names 全部能在 M1 掩码中找到", aligned["nCount"].notna().all(),
        f"缺 {int(aligned['nCount'].isna().sum())}")
    # 累加一律用 float64：float32 在**逐样本**求和的量级（可达 4e7 > 2^24）会因累加舍入
    # 产生 ±2 的**假偏差**（见 A3 的 note）。逐细胞量级（max nCount 777,447 < 2^24）本可精确，
    # 但统一用 float64 以消除歧义。
    rowsum = np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
    d2 = int(np.abs(rowsum - aligned["nCount"].to_numpy()).max())
    chk("A2b **逐细胞**行和 == M1 nCount", d2 == 0, f"最大偏差 {d2}")
    if d2:
        w = np.argmax(np.abs(rowsum - aligned["nCount"].to_numpy()))
        notes.append(f"    偏差最大处：{A.obs_names[w]} 矩阵 {rowsum[w]:.0f} vs M1 {aligned['nCount'].iloc[w]}")

    # V3b 的互补检查
    nnzc = np.asarray(X.getnnz(axis=1)).ravel()
    d2b = int(np.abs(nnzc.astype(np.int64) - aligned["nFeature"].to_numpy().astype(np.int64)).max())
    chk("A2c **逐细胞** nnz == M1 nFeature（与 V3b 同）", d2b == 0, f"最大偏差 {d2b}")

    # (nCount, nFeature) 组合唯一性 → 行列错位若发生，几乎必然被 A2b 抓到
    pair = pd.DataFrame({"s": rowsum.astype(np.int64), "f": nnzc.astype(np.int64)})
    dup = int(len(pair) - len(pair.drop_duplicates()))
    notes.append(f"    (行和, 非零数) 组合重复 {dup:,} / {len(pair):,} "
                 f"（重复越高，A2b 对'两细胞互换'的检出力越弱）")

    per_sample_sum = pd.Series(rowsum, index=A.obs["sample_id"].astype(str).to_numpy()).groupby(level=0).sum()
    m1_sum = m1.groupby("sample_id")["nCount"].sum()
    cmp = pd.concat([per_sample_sum.rename("built"), m1_sum.rename("m1")], axis=1).fillna(0).astype(np.int64)
    d3 = int((cmp["built"] - cmp["m1"]).abs().max())
    chk("A3 逐样本行和 == M1 逐样本 Σ nCount", d3 == 0, f"最大偏差 {d3}；样本 {len(cmp)}")
    # 如实记录：若用 float32 累加会出现的**假**偏差，说明为何本检查必须 float64/int64
    rs32 = np.asarray(X.sum(axis=1, dtype=np.float32)).ravel()
    s32 = (pd.Series(rs32, index=A.obs["sample_id"].astype(str).to_numpy())
           .groupby(level=0).sum().reindex(m1_sum.index).astype(np.int64))
    ref = m1_sum.astype(np.int64)
    d32 = int((s32 - ref).abs().max())
    notes.append(f"    float32 累加会给出假偏差 max|Δ|={d32}（{int((s32 != ref).sum())}/{len(ref)} 样本）"
                 f" —— 纯累加舍入（样本和量级 > 2^24），非数据错误；A3 因此用 float64")

    # ================================================================ A4 全局
    data = X.data
    n = len(data)
    bad_fin = bad_int = bad_neg = 0
    gmax, gmin = -np.inf, np.inf
    for i in range(0, n, 20_000_000):
        b = data[i:i + 20_000_000]
        bad_fin += int((~np.isfinite(b)).sum())
        bad_int += int((np.mod(b, 1) != 0).sum())
        bad_neg += int((b < 0).sum())
        gmax = max(gmax, float(b.max()))
        gmin = min(gmin, float(b.min()))
    chk("A4a X.data 全体有限（无 NaN/Inf）", bad_fin == 0, f"异常 {bad_fin} / {n:,}")
    chk("A4b X.data 全体为整数", bad_int == 0, f"非整 {bad_int} / {n:,}（V6b 只抽查前 500 万）")
    chk("A4c X.data 全体非负", bad_neg == 0, f"负值 {bad_neg}")
    chk("A4d float32 整数精度无损", gmax < 2 ** 24,
        f"max {gmax:.0f} < 2^24={2**24}，min {gmin:.0f}")

    # ================================================================ A6 R1 闭合
    geo = pd.read_csv(GEO_TSV, sep="\t")
    geo["sample_id"] = geo["patient_id"].astype(str) + "_" + geo["token"].astype(str)
    gmap = geo.set_index("sample_id")
    obs_pairs = A.obs[["sample_id", "patient_id", "stage_token", "stage"]].drop_duplicates()
    obs_pairs = obs_pairs.assign(sample_id=obs_pairs["sample_id"].astype(str))
    missing = sorted(set(obs_pairs["sample_id"]) - set(gmap.index))
    chk("A6a 每个 obs 样本都在 GEO 权威表", not missing, f"缺 {missing[:5]}")
    ok_p = ok_t = ok_s = 0
    bad_rows = []
    for r in obs_pairs.itertuples():
        if r.sample_id not in gmap.index:
            continue
        g = gmap.loc[r.sample_id]
        if str(g["patient_id"]) != str(r.patient_id):
            bad_rows.append((r.sample_id, "patient", g["patient_id"], r.patient_id))
        elif str(g["token"]) != str(r.stage_token):
            bad_rows.append((r.sample_id, "token", g["token"], r.stage_token))
        elif REG.resolve_stage(str(g["token"])) != str(r.stage):
            bad_rows.append((r.sample_id, "stage", g["token"], r.stage))
        else:
            ok_s += 1
    chk("A6b patient_id / stage_token / stage 三者与 GEO 权威表逐行一致",
        not bad_rows, f"{ok_s}/{len(obs_pairs)} 一致；不一致 {bad_rows[:5]}")

    # A6c：**本项目设计是同一患者多个分期**（同患者多病灶、多分期 = 配对设计本身，
    # 也正因此 patient_id 绝不可作 batch）。故"患者内 stage 唯一"是**错误断言**——
    # 本脚本初版误写成该断言并误报 FAIL，已更正为下面对 GUID 权威表的集合一致性检查。
    o_set = (A.obs[["sample_id", "patient_id", "stage_token"]]
             .drop_duplicates()                      # ← 必须先去重：obs 是逐细胞的，token 会重复
             .assign(pid=lambda d: d["patient_id"].astype(str))
             .groupby("pid")["stage_token"]
             .apply(lambda s: tuple(sorted(s.astype(str)))))
    g_set = geo.groupby("patient_id")["token"].apply(lambda s: tuple(sorted(s.astype(str))))
    g_set.index = g_set.index.astype(str)
    allp = sorted(set(g_set.index) | set(o_set.index))
    diff_p = {p: [list(g_set.get(p, ())), list(o_set.get(p, ()))] for p in allp
              if tuple(g_set.get(p, ())) != tuple(o_set.get(p, ()))}
    chk("A6c 每患者的 stage_token 集合 == GEO 权威表", not diff_p,
        f"患者 {len(allp)}；不一致 {diff_p}")
    chk("A6d GEO 表内无重复 (patient, token)", not geo.duplicated(["patient_id", "token"]).any(),
        f"重复 {int(geo.duplicated(['patient_id', 'token']).sum())}")
    nsamp = A.obs.assign(pid=A.obs["patient_id"].astype(str)).groupby("pid")["sample_id"].nunique()
    notes.append(f"    每患者样本数分布 {nsamp.value_counts().sort_index().to_dict()}"
                 f"（>1 = 配对设计正常）")

    # ================================================================ A7/A8 逐样本
    order = sorted(raw_files)
    bad_npz, bad_bc, bad_shape = [], [], []
    sid_all = A.obs["sample_id"].astype(str).to_numpy()
    offsets, off = {}, 0
    for s in order:
        offsets[s] = off
        off += int((sid_all == s).sum())
    chk("A7a 逐样本 obs 计数合计 == n_obs", off == EXPECT_N_OBS, f"{off} vs {EXPECT_N_OBS}")
    # A7a' 非平凡检查：每个样本的**行块必须连续**（否则 A7b 的块切片会取到别的样本）
    bad_contig = [s for s in order
                  if not np.all(sid_all[offsets[s]:offsets[s] + int((sid_all == s).sum())] == s)]
    chk("A7a' 每个样本的行块在 obs 中连续", not bad_contig, f"不连续 {bad_contig[:3]}")

    import random
    rng = random.Random(20260915)
    pool = order
    check_npz = order if len(order) <= 20 else rng.sample(pool, 12)
    for s in check_npz:
        p = os.path.join(PER_SAMPLE, f"{s}.npz")
        if not os.path.exists(p):
            bad_npz.append((s, "缺文件"))
            continue
        B = sp.load_npz(p).tocsr()
        m = (A.obs["sample_id"].astype(str) == s).to_numpy()
        blk = X[offsets[s]:offsets[s] + int(m.sum())]
        if B.shape != blk.shape:
            bad_shape.append((s, B.shape, blk.shape))
            continue
        if (B != blk).nnz != 0:
            d = (B - blk)
            bad_npz.append((s, f"元素不一致 nnz={d.nnz}"))
    chk("A7b 抽样 npz ↔ h5ad 对应块逐元素一致", not (bad_npz or bad_shape),
        f"抽 {len(check_npz)} 个；异常 {bad_npz[:3]} {bad_shape[:3]}")

    # A8 barcode 集合
    for s in check_npz:
        with gzip.open(raw_files[s], "rt") as fh:
            file_bcs = set(fh.readline().rstrip("\n").split("\t"))
        want = set(m1.loc[m1["sample_id"] == s, "bc10x"])
        if not want.issubset(file_bcs):
            bad_bc.append((s, len(want - file_bcs)))
    chk("A8a 原始文件 barcode ⊇ M1 掩码 barcode", not bad_bc, f"不足 {bad_bc[:3]}")

    obs_bc = set(A.obs_names.to_numpy())
    chk("A8b h5ad obs barcode 集合 == M1 掩码 barcode 集合",
        obs_bc == set(m1["cell_barcode"]),
        f"仅矩阵有 {len(obs_bc - set(m1['cell_barcode']))}；仅 M1 有 {len(set(m1['cell_barcode']) - obs_bc)}")

    # ================================================================ A5 独立重抽取
    if args.n_reextract > 0:
        cnt = pd.Series(A.obs["sample_id"].astype(str)).value_counts()
        cnt = cnt.reindex(order)
        # 选：最小样本、中位样本、以及一个带 '1' 后缀的样本
        picks = [cnt.idxmin(), cnt.sort_values().index[len(cnt) // 2]]
        ones = [s for s in order if s[-1] == "1"]
        if ones:
            picks.append(ones[0])
        picks = list(dict.fromkeys(picks))[:args.n_reextract]
        notes.append(f"    A5 重抽取样本：{picks}（细胞数 {[int(cnt[s]) for s in picks]}）")

        indep_fail = []
        for s in picks:
            with gzip.open(raw_files[s], "rt") as fh:
                bcs = fh.readline().rstrip("\n").split("\t")
                # ↓ 与构建器不同的路径：整表一次读入，让 pandas 自行推断 dtype
                df = pd.read_csv(fh, sep="\t", header=None, index_col=0)
            g_ind = df.index.to_numpy(dtype=object)
            if not np.array_equal(g_ind.astype(str), vn.astype(str)):
                indep_fail.append((s, "基因列与 var_names 不一致"))
                continue
            col = {bc: i for i, bc in enumerate(bcs)}
            want = sorted(m1.loc[m1["sample_id"] == s, "bc10x"])
            if any(bc not in col for bc in want):
                indep_fail.append((s, "掩码 barcode 缺失"))
                continue
            take = np.fromiter((col[bc] for bc in want), dtype=np.int64, count=len(want))
            exp = df.to_numpy(dtype=np.float32)[:, take].T     # → 细胞 × 基因
            del df
            m = (A.obs["sample_id"].astype(str) == s).to_numpy()
            got = X[offsets[s]:offsets[s] + int(m.sum())].toarray()
            if exp.shape != got.shape:
                indep_fail.append((s, f"shape {exp.shape} vs {got.shape}"))
                continue
            d = np.abs(exp - got).max()
            if d != 0:
                nz = int((exp != got).sum())
                indep_fail.append((s, f"元素不一致 max|Δ|={d} 不同元素 {nz}"))
            del exp, got
        chk(f"A5 独立重抽取 {len(picks)} 样本 × 逐元素比对 h5ad", not indep_fail,
            f"不一致 {indep_fail}")

    # ================================================================ 报告
    print("\n============== GP0 对抗性审计 ==============")
    for x in notes:
        print("  ·", x)
    print("-------------------------------------------")
    for x in passes:
        print("  ", x)
    for x in fails:
        print("  ", x)
    print("-------------------------------------------")
    print(f"  PASS {len(passes)} / FAIL {len(fails)}")
    print("VERDICT:", "GP0 审计 全绿 ✅" if not fails else "GP0 审计 有失败 ❌")
    print("===========================================")
    return 2 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
