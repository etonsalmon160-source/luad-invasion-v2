#!/usr/bin/env python3
"""指标3（谱系覆盖）补算 —— 全量对象 + 6 个谱系子集，逐个分辨率。

登记原文（PARAMETERS_AND_SOURCES.md §M3-A.3 第 237 行）
------------------------------------------------------
    「有 ≥1 个法则2 marker 集在 ≥25% 细胞检出且模块分均值 >0 的簇占比，硬约束 ≥ 0.90」

⚠️ **「marker 集在 ≥25% 细胞检出」这句有歧义**，本脚本把三种读法**全部**算出来并列上报，
   **不挑一个好看的当结论**：
      R_any  ：该 panel **任一**基因在 ≥25% 细胞检出（最松）
      R_mean ：该 panel 各基因检出率的**均值** ≥25%   ← 本次头条口径
      R_all  ：该 panel **全部**基因都在 ≥25% 细胞检出（最严）
   「模块分均值 >0」这一半无歧义：指该簇在**同一个** panel 上的模块分均值 >0。
   口径未定 ⇒ 读数并列，**由用户裁定**，本脚本不代为选择。

对象与分辨率
------------
    full r*=0.6 ；epiA 0.7 ；tnkA 0.8 ；myeloidA 0.5 ；endoA 0.7 ；bplasmaA 0.6 ；fibroA 0.8
    （2026-09-22 更新；换 r* 前为 epiA 0.5 / tnkA 0.5 / endoA 0.5）
    每对象在 {0.5,0.6,0.7,0.8} **四档全算**，每档 **5 个种子全算**（r* 由跨种子均值选出，
    这里看种子敏感性）。

模块分怎么来的
--------------
    与 GP6 标准 A **同一实现**：normalize_total(1e4) + log1p → sc.tl.score_genes(ctrl_size=50,
    random_state=0)，基因轴保留全部 18,069 个（只切细胞，不切基因），缺失基因经
    marker_panel.resolve_missing() 剔除并上报。
    ⚠️ 子集上的模块分**与全量上重算的不是同一个数**（score_genes 的对照基因分箱依赖
       传入对象的整体表达分布）。子集就该用子集自己的分，本脚本如此做。

自校验
------
    全量对象的模块分**重算一遍**，与已冻结的 `gp6_scores.csv.gz` 逐细胞比对；
    |差| 超容差即硬报错。重算对得上 ⇒ 子集部分才可信。

🔴 第一版踩过的坑（留在注释里，防复发）
    第一版对全量对象走的是 `sub_adata = A`（不拷贝），随后在其上原地 normalize+log1p，
    **污染了原始计数对象** ⇒ 后 6 个子集是从一份已 log1p 的对象上切的，被双重归一化，
    模块分全错（检出那半因零位置不变而未受影响）。
    本版起：**任何 normalize 都只在 `A[ix].copy()` 上做，原始 A 只读**，
    并在每个对象开始前用「加载时留下的 data 切片金丝雀」断言 A 未被改动。

用法:
    python3 05_annotation/04_metric3_coverage.py
"""
import hashlib
import json
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import marker_panel as MP

ROOT = "/home/eto/luad_v2"
H5 = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
SDIR = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ROOT}/results/05_annotation"
GP6_SCORES = f"{OUT}/gp6_scores.csv.gz"

CTRL_SIZE = 50
SCORE_SEED = 0
RES_GRID = [0.5, 0.6, 0.7, 0.8]
SEEDS = [0, 1, 2, 3, 4]
DET_MIN = 0.25          # 登记原文的 25%
COV_MIN = 0.90          # 登记原文的硬约束
HEADLINE = "R_mean"     # 头条口径；另两种并列上报
TOL = 1e-6              # 全量模块分重算比对容差
CANARY = 4096           # 原始计数未被改动的金丝雀长度

# ⚠️ 各对象的 r* 是**手抄**自 `results/05_annotation/<tag>_rstar.json`，不是运行时读取。
# 2026-09-22 更新：L1 标签口径由 `A_frozen` 换成 `A_adjudicated` ⇒ 谱系成员变化 ⇒
# 六谱系重聚类 ⇒ **上皮 0.5→0.7、T·NK 0.5→0.8、内皮 0.5→0.7**（髓系/B·浆/成纤维不变）。
# 换 r* 后必须重跑本脚本，否则 `coverage.csv` 里「r* 处」那几行取的是**已作废的分辨率**。
OBJECTS = [
    ("full",     None,                           0.6),
    ("epiA",     "epiA_subset_barcodes.txt",     0.7),
    ("tnkA",     "tnkA_subset_barcodes.txt",     0.8),
    ("myeloidA", "myeloidA_subset_barcodes.txt", 0.5),
    ("endoA",    "endoA_subset_barcodes.txt",    0.7),
    ("bplasmaA", "bplasmaA_subset_barcodes.txt", 0.6),
    ("fibroA",   "fibroA_subset_barcodes.txt",   0.8),
]
LIN_CN = {"full": "全量", "epiA": "上皮", "tnkA": "T/NK", "myeloidA": "髓系",
          "endoA": "内皮", "bplasmaA": "B/浆", "fibroA": "成纤维"}
# §六.6 第 2 条登记要求：除登记字面的「任意一套 panel」外，**另算「本谱系那套」**并两数列出，
# 用于检验「错认」口子有多大。全量对象无「本谱系」，该列不适用。
OWN_LIN = {"epiA": "上皮", "tnkA": "T/NK", "myeloidA": "髓系",
           "endoA": "内皮", "bplasmaA": "B/浆", "fibroA": "成纤维"}

LOG = []


def log(s):
    print(s, flush=True)
    LOG.append(s)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


# ================================================================ 载入全量计数（只读）
log(f"载入 {os.path.relpath(H5, ROOT)}")
A = ad.read_h5ad(H5)
log(f"  shape={A.shape}  nnz={A.X.nnz:,}")
assert A.X.min() >= 0 and np.allclose(A.X.data, np.round(A.X.data)), "X 不是原始整数计数"
assert not A.obs_names.duplicated().any(), "barcode 有重复"
A_X_CSR = A.X.tocsr()
_CANARY_VAL = A_X_CSR.data[:CANARY].copy()       # 金丝雀：A 若被原地改动，断言会炸
_NNZ0 = int(A_X_CSR.nnz)

USED, MISSING = MP.resolve_missing(A.var_names)
log(f"标准 A 面板：可用 {sum(len(v) for v in USED.values())} / 缺失 {MISSING}")
for lin in MP.LINEAGES:
    assert len(USED[lin]) >= 5, f"{lin} 可用 marker 过少（{len(USED[lin])}）"

# 76 个 panel 基因一次抽成小矩阵，检出全在这里算，快且不碰大对象
PANEL_GENES = [g for lin in MP.LINEAGES for g in USED[lin]]
GI = {g: i for i, g in enumerate(PANEL_GENES)}
assert len(PANEL_GENES) == len(GI), "panel 基因有重复，需去重后再定列号"
B = (A_X_CSR[:, [list(A.var_names).index(g) for g in PANEL_GENES]] > 0).tocsr()
B.data = np.ones_like(B.data, dtype=bool)        # 纯 bool，检出只看 0/1
PANEL_COL = {lin: np.array([GI[g] for g in USED[lin]]) for lin in MP.LINEAGES}
log(f"  检出矩阵 B：{B.shape}  nnz={B.nnz:,}（{len(PANEL_GENES)} 个 panel 基因）")


def content_sha256(path, chunk=1 << 22):
    """解压后内容的 sha256（`.gz` 用）。非 `.gz` 退化为同 sha256。

    ⚠️ 为什么需要它：`.gz` 的**压缩文件**哈希里含 gzip MTIME 字段 ⇒ 即使内容一字未改，
    逐次运行也会得到不同的哈希（2026-09-22 实测：08b967d8… → 794a222c…，内容相同）。
    压缩哈希因此**不可作复现锚点**（违反 R5「产物可复现且哈希」的本意），内容哈希才是。
    """
    if not path.endswith(".gz"):
        return sha256(path)
    import gzip
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def assert_raw_untouched(stage):
    assert int(A.X.nnz) == _NNZ0, f"🔴 {stage}：原始对象 nnz 变了，A 被改动"
    assert np.array_equal(A.X.tocsr().data[:CANARY], _CANARY_VAL), \
        f"🔴 {stage}：原始计数的金丝雀变了，A 被原地改动"


def cluster_stats(Bc, codes, c):
    """一个簇：6 panel 的检出三读数 + 模块分均值。Bc=该对象检出矩阵，codes=簇编号数组。"""
    m = codes == c
    Bcc = Bc[m]
    n = Bcc.shape[0]
    row = {"n_cells": int(n)}
    for lin in MP.LINEAGES:
        # 每基因检出率：该簇里 count>0 的细胞占比
        per_gene = np.asarray(Bcc[:, PANEL_COL[lin]].sum(axis=0)).ravel() / n
        row[f"any_{lin}"] = float(np.asarray(Bcc[:, PANEL_COL[lin]].sum(axis=1) > 0).mean())
        row[f"mean_{lin}"] = float(per_gene.mean())
        row[f"all_{lin}"] = float(per_gene.min())
    return row


# ================================================================ 逐对象
per_cluster_rows = []
cov_rows = []
recompute_check = None

for tag, bc_file, rstar in OBJECTS:
    t0 = time.time()
    assert_raw_untouched(f"{tag} 开始前")

    if bc_file is None:
        ix = np.arange(A.n_obs)
        cells = A.obs_names.to_numpy()
    else:
        path = os.path.join(OUT, bc_file)
        cells = np.array([l.strip() for l in open(path) if l.strip()])
        assert len(set(cells)) == len(cells), f"{bc_file} 有重复"
        ix = pd.Index(A.obs_names).get_indexer(cells)
        assert (ix >= 0).all(), f"{bc_file} 有 barcode 不在 h5ad 里"

    log(f"\n=== {tag}（{LIN_CN[tag]}）n_cells={len(ix):,}  r*={rstar} ===")

    Bsub = B[ix]                                  # 检出：只读小矩阵

    # 模块分：**在拷贝上**做归一化，原始 A 一个字节都不动
    sub = A[ix].copy()
    assert_raw_untouched(f"{tag} 取拷贝后")
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    for lin in MP.LINEAGES:
        sc.tl.score_genes(sub, list(USED[lin]), ctrl_size=CTRL_SIZE,
                          random_state=SCORE_SEED, score_name=f"scoreA_{lin}")

    if bc_file is None:
        log("  自校验：全量模块分与已冻结 gp6_scores.csv.gz 比对")
        ref = pd.read_csv(GP6_SCORES)
        assert (ref["cell_barcode"].to_numpy() == sub.obs_names.to_numpy()).all(), \
            "gp6_scores 与 h5ad 细胞顺序不一致"
        worst = 0.0
        for lin in MP.LINEAGES:
            d = np.abs(sub.obs[f"scoreA_{lin}"].to_numpy() - ref[lin].to_numpy(dtype=float))
            worst = max(worst, float(d.max()))
        log(f"  自校验最大 |差| = {worst:.3e}（容差 {TOL:g}）")
        assert worst <= TOL, f"🔴 全量模块分重算与冻结值不符（最大差 {worst:.3e}）⇒ 停"
        recompute_check = dict(max_abs_diff=worst, tol=TOL, status="PASS")

    scores = {lin: sub.obs[f"scoreA_{lin}"].to_numpy(dtype=np.float64) for lin in MP.LINEAGES}

    clu = pd.read_csv(os.path.join(SDIR, tag, "clusters.csv.gz"))
    assert (clu["cell_barcode"].to_numpy() == sub.obs_names.to_numpy()).all(), \
        f"{tag} 的分群表与细胞顺序不一致"

    for r in RES_GRID:
        for s in SEEDS:
            col = f"harmony_res{r}_seed{s}"
            assert col in clu.columns, col
            lab = pd.Categorical(clu[col].astype(str))
            codes = lab.codes
            nclu = int(codes.max()) + 1
            stats = []
            for c in range(nclu):
                row = cluster_stats(Bsub, codes, c)
                m = codes == c
                for lin in MP.LINEAGES:
                    row[f"score_{lin}"] = float(scores[lin][m].mean())
                row.update(object=tag, res=r, seed=s, cluster=str(lab.categories[c]))
                stats.append(row)
                per_cluster_rows.append(row)
            df = pd.DataFrame(stats)
            for reading, col_ in (("R_any", "any"), ("R_mean", "mean"), ("R_all", "all")):
                ok = np.zeros(len(df), dtype=bool)
                hit = []
                for lin in MP.LINEAGES:
                    pc = ((df[f"{col_}_{lin}"].to_numpy() >= DET_MIN) &
                          (df[f"score_{lin}"].to_numpy() > 0))
                    ok |= pc
                    if pc.any():
                        hit.append(f"{lin}:{int(pc.sum())}")
                cov_rows.append(dict(object=tag, res=r, seed=s, reading=reading,
                                     n_clusters=int(len(df)), n_pass=int(ok.sum()),
                                     coverage=float(ok.mean()),
                                     passes=bool(ok.mean() >= COV_MIN),
                                     panels_hit=";".join(hit) or "-"))
            # 「本谱系那套」口径（§6.6 第 2 条）：只算该对象自己那套 panel，头条读法
            if tag in OWN_LIN:
                lin = OWN_LIN[tag]
                pc_own = ((df[f"{HEADLINE.split('_')[1]}_{lin}"].to_numpy() >= DET_MIN) &
                          (df[f"score_{lin}"].to_numpy() > 0))
                cov_rows.append(dict(object=tag, res=r, seed=s,
                                     reading=f"{HEADLINE}_own",
                                     n_clusters=int(len(df)), n_pass=int(pc_own.sum()),
                                     coverage=float(pc_own.mean()),
                                     passes=bool(pc_own.mean() >= COV_MIN),
                                     panels_hit=f"{lin}:{int(pc_own.sum())}"))
    log(f"  完成，用时 {time.time()-t0:.1f}s")

assert_raw_untouched("全部对象跑完")

cov = pd.DataFrame(cov_rows)
pc = pd.DataFrame(per_cluster_rows)

# ================================================================ 落盘
os.makedirs(OUT, exist_ok=True)
cov_path = os.path.join(OUT, "metric3_coverage.csv")
cov.to_csv(cov_path, index=False)
pc_path = os.path.join(OUT, "metric3_per_cluster.csv.gz")
pc.to_csv(pc_path, index=False, compression="gzip")

# ================================================================ 汇总打印（r* 处，seed0）
log("\n" + "=" * 78)
log(f"头条口径 {HEADLINE}（panel 各基因检出率均值 ≥ {DET_MIN}）· r* 处 · seed0")
log("=" * 78)
log(f"{'对象':<9}{'r*':>5}{'簇数':>6}{'通过':>6}{'覆盖':>10}   判定")
for tag, _, rstar in OBJECTS:
    q = cov[(cov.object == tag) & (cov.res == rstar) & (cov.seed == 0) &
            (cov.reading == HEADLINE)].iloc[0]
    verdict = "✅ 过线" if q["coverage"] >= COV_MIN else "🔴 未过线"
    log(f"{LIN_CN[tag]:<9}{rstar:>5}{q['n_clusters']:>6}{q['n_pass']:>6}"
        f"{q['coverage']:>10.4f}   {verdict}")

log("\n" + "=" * 78)
log("三种读法在 r* 处（seed0）的并列读数 —— 口径未定，不代为选择")
log("=" * 78)
log(f"{'对象':<9}" + "".join(f"{rd:>17}" for rd in ("R_any", "R_mean", "R_all")))
for tag, _, rstar in OBJECTS:
    cs = []
    for rd in ("R_any", "R_mean", "R_all"):
        q = cov[(cov.object == tag) & (cov.res == rstar) & (cov.seed == 0) &
                (cov.reading == rd)].iloc[0]
        cs.append(f"{q['coverage']:.4f}({q['n_pass']}/{q['n_clusters']})")
    log(f"{LIN_CN[tag]:<9}" + "".join(f"{c:>17}" for c in cs))

log("\n" + "=" * 78)
log("§6.6 第 2 条：「本谱系那套」口径（R_mean_own，r* 处，seed0）—— 检验「错认」口子")
log("=" * 78)
for tag, _, rstar in OBJECTS:
    if tag not in OWN_LIN:
        log(f"{LIN_CN[tag]:<9}{rstar:>5}{'-':>6}{'-':>6}{'N/A':>10}   全量对象无「本谱系」")
        continue
    q = cov[(cov.object == tag) & (cov.res == rstar) & (cov.seed == 0) &
            (cov.reading == f"{HEADLINE}_own")].iloc[0]
    verdict = "✅ 仍过线" if q["coverage"] >= COV_MIN else "🔴 收紧后不过线"
    log(f"{LIN_CN[tag]:<9}{rstar:>5}{q['n_clusters']:>6}{q['n_pass']:>6}"
        f"{q['coverage']:>10.4f}   {verdict}")

log("\n" + "=" * 78)
log(f"种子敏感性（r* 处，{HEADLINE}）—— r* 由跨种子均值选出，这里看逐种子")
log("=" * 78)
for tag, _, rstar in OBJECTS:
    v = cov[(cov.object == tag) & (cov.res == rstar) &
            (cov.reading == HEADLINE)].sort_values("seed")
    log(f"{LIN_CN[tag]:<9}r*={rstar}  " +
        "  ".join(f"s{int(s)}={c:.4f}" for s, c in zip(v.seed, v.coverage)) +
        f"   最低={v.coverage.min():.4f}{'  🔴 有种子上不了线' if v.coverage.min() < COV_MIN else ''}")

log("\n" + "=" * 78)
log(f"跨分辨率（seed0，{HEADLINE}）—— 若某谱系 r* 处不过线，这里看网格内还有没有别的档")
log("=" * 78)
for tag, _, rstar in OBJECTS:
    v = cov[(cov.object == tag) & (cov.seed == 0) &
            (cov.reading == HEADLINE)].sort_values("res")
    log(f"{LIN_CN[tag]:<9}" +
        "  ".join(f"r{r}{'*' if r == rstar else ' '}={c:.4f}" for r, c in zip(v.res, v.coverage)))

# ── 两段说明改为**由 cov 现算**，不再写死 ───────────────────────────────────
# 2026-09-22 教训：旧文本写死「七个对象全部过线」。换 L1 口径后髓系子集由 54,407
# 涨到 64,084，该句立刻变成假话，而它已经写进了 manifest。⇒ 凡随数据变化的结论，
# 一律从表里算出来，不写字符串常量。
_rm = cov[(cov.reading == HEADLINE) & (cov.seed == 0)]
_at_rstar = [(t, float(_rm[(_rm.object == t) & (_rm.res == r)].coverage.iloc[0]),
              int(_rm[(_rm.object == t) & (_rm.res == r)].n_pass.iloc[0]),
              int(_rm[(_rm.object == t) & (_rm.res == r)].n_clusters.iloc[0]))
             for t, _, r in OBJECTS]
_n_fail = [x for x in _at_rstar if x[1] < COV_MIN]
_n_below = int((cov[cov.reading == HEADLINE].coverage < COV_MIN).sum())
_p = [f"{HEADLINE} 下各对象 r* 处（seed0）：",
      " · ".join(f"{LIN_CN[t]} {c:.4f}({p_}/{n})" for t, c, p_, n in _at_rstar),
      f"；最低 {min(x[1] for x in _at_rstar):.4f}。"]
if _n_fail:
    _p.append("⇒ 🔴 **%d 个对象在 r* 处不过线**（%s）⇒ 本指标**已淘汰这些谱系的 r* 候选**，"
              "须按 §M3-A.3 停在检查点升级，或由用户**显式裁定**覆盖（不得默许、不得事后放宽阈值）。"
              % (len(_n_fail), "、".join(LIN_CN[t] for t, *_ in _n_fail)))
else:
    _p.append("⇒ 全部过线 ⇒ 指标3 对 r* 无区分力，未排除任何分辨率。")
_p.append(f"全网格 {_n_below} 格低于 {COV_MIN}。不得把「指标3 过线」写成对 r* 的支持证据。")
non_binding_note = "".join(_p)

_ow = cov[(cov.reading == f"{HEADLINE}_own") & (cov.seed == 0)]
_ow_rstar = [(t, float(_ow[(_ow.object == t) & (_ow.res == r)].coverage.iloc[0]),
              int(_ow[(_ow.object == t) & (_ow.res == r)].n_pass.iloc[0]),
              int(_ow[(_ow.object == t) & (_ow.res == r)].n_clusters.iloc[0]))
             for t, _, r in OBJECTS if t in OWN_LIN]
_ow_fail = [x for x in _ow_rstar if x[1] < COV_MIN]
_ow_grid = cov[(cov.reading == f"{HEADLINE}_own") & (cov.seed == 0) & (cov.coverage < COV_MIN)]
_p = ["收紧到「本谱系那套」后，各谱系 r* 处（seed0）：",
      " · ".join(f"{LIN_CN[t]} {c:.4f}({p_}/{n})" for t, c, p_, n in _ow_rstar)]
_p.append("⇒ 🔴 **" + "、".join(LIN_CN[t] for t, *_ in _ow_fail) + " 在 r* 处不过线**。"
          if _ow_fail else "⇒ r* 处全过线。")
if len(_ow_grid):
    _p.append("全网格另有 %d 格低于 %s：" % (len(_ow_grid), COV_MIN) + "；".join(
        f"{LIN_CN[row.object]} r={row.res} s{int(row.seed)}={row.coverage:.4f}"
        for row in _ow_grid.itertuples()))
own_panel_note = "".join(_p)

man = dict(
    script="05_annotation/04_metric3_coverage.py",
    registered_definition="PARAMETERS_AND_SOURCES.md §M3-A.3 的「指标3 谱系覆盖」行（**不引行号**：该表会增行，行号会漂）",
    registered_text="有 ≥1 个法则2 marker 集在 ≥25% 细胞检出且模块分均值 >0 的簇占比，硬约束 ≥ 0.90",
    caliber_ambiguity_note=(
        "登记原文的「marker 集在 ≥25% 细胞检出」有歧义 ⇒ 三种读法并列上报："
        "R_any=任一基因检出率≥25%；R_mean=各基因检出率均值≥25%；R_all=全部基因检出率≥25%。"
        "头条取 R_mean，但**口径由用户裁定**，本脚本不代为选择。"
        "另按 §六.6 第 2 条增列 R_mean_own = 只认「本谱系那套」panel（收紧口径），"
        "用于检验「被错认也算认得出」这个口子有多大；全量对象无「本谱系」，该列不适用。"),
    headline_reading=HEADLINE, det_min=DET_MIN, coverage_min=COV_MIN,
    module_score=dict(method="sc.tl.score_genes", ctrl_size=CTRL_SIZE,
                      random_state=SCORE_SEED,
                      norm="normalize_total(1e4) + log1p",
                      gene_axis="全部 18,069 基因，仅切细胞不切基因",
                      note="子集模块分按子集自身分布重算，与全量上重算的同细胞分值不同"),
    detection=dict(rule="count>0", matrix="单独抽出的 76 基因 bool 小矩阵",
                   invariant_note="归一化不改零位置 ⇒ 检出与归一化无关"),
    panel_source="05_annotation/marker_panel.py（法则2，一次文献）",
    panel_genes_used=PANEL_GENES, panel_missing=MISSING,
    all_resolutions_computed=RES_GRID, all_seeds_computed=SEEDS,
    per_object_rstar={t: r for t, _, r in OBJECTS},
    recompute_selfcheck=recompute_check,
    known_bug_fixed=("v1 对 full 走 sub_adata=A 且原地归一化 ⇒ 污染原始对象、后 6 子集被双重归一化。"
                     "v2 起归一化只在拷贝上做，并用金丝雀断言原始对象未被改动。v1 结果作废。"),
    caliber_ruling=dict(
        date="2026-09-21", ruled_by="用户（会话内批准）",
        decision="登记口径取 R_mean（panel 各基因检出率均值 ≥25%）",
        rationale=("另两种读法可被证明是退化的，不是挑好看的那个："
                   "① R_any 恒 = 1.0 —— 面板里只要有一个基因在 ≥25% 细胞里冒头就算过"
                   "（环境 RNA 带进来一点即够），等于没测；"
                   "② R_all 结构性不可达 —— 面板故意跨亚型（上皮那套同时含 AT1 的 AGER/CAV1 与 "
                   "AT2 的 SFTPC/NAPSA），没有任何单一细胞类型能把整组基因都表达在 ≥25% 细胞里。"),
        non_binding_note=non_binding_note,
        own_panel_note=own_panel_note,
        recompute_note=("🔴 2026-09-22 重算：L1 标签口径由 `A_frozen` 换为 `A_adjudicated` ⇒ 谱系成员变化 "
                        "⇒ 按脚本头部注释的要求重跑。旧产物备份于 `results/05_annotation/.prev_metric3_20260921/`。"
                        "本次两处改动：① OBJECTS 里手抄的 r* —— epiA 0.5→0.7、tnkA 0.5→0.8、endoA 0.5→0.7"
                        "（髓系/B·浆/成纤维不变）；② `non_binding_note` 与 `own_panel_note` 改为**由 cov 现算** —— "
                        "旧文本写死「七个对象全部过线」，换子集后即成假话。")),
    inputs={os.path.relpath(p, ROOT): sha256(p) for p in
            [H5, GP6_SCORES, f"{ROOT}/05_annotation/marker_panel.py"]
            + [f"{OUT}/{b}" for _, b, _ in OBJECTS if b]},
    outputs={},
)
man["outputs"] = {os.path.relpath(p, ROOT): sha256(p) for p in [cov_path, pc_path]}
man["outputs_content_sha256"] = {os.path.relpath(p, ROOT): content_sha256(p)
                                for p in [cov_path, pc_path]}
man["hash_note"] = ("⚠️ `outputs` 里 `.gz` 项是**压缩文件**的 sha256，含 gzip MTIME ⇒ 逐次运行必变，"
                    "**不可作复现锚点**；`outputs_content_sha256` 记的是**解压后内容**的 sha256，"
                    "它才是可复现的锚点（R5）。非 `.gz` 项两者相同。")
man_path = os.path.join(OUT, "metric3_manifest.json")
with open(man_path, "w") as fh:
    json.dump(man, fh, ensure_ascii=False, indent=2)
log(f"\n写出 {os.path.relpath(cov_path, ROOT)} / {os.path.relpath(pc_path, ROOT)} / "
    f"{os.path.relpath(man_path, ROOT)}")
