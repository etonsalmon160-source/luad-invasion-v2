#!/usr/bin/env python3
# 05_spatial_stats.py —— M6 §3.7 Squidpy 空间统计（邻域富集 + 共现）
#
# 输入：
#   results/10_niche/consensus/<agf>/consensus_domains.tsv    域 → 共识类型（只含选中那档）
#   results/10_niche/banksy/<slide>/domains_<agf>.tsv.gz      逐 spot 域分配（seed=0）
#   data/visium_spatial/<slide>/spatial/tissue_positions.csv  厂商 pxl 坐标
#
# 输出（results/10_niche/spatial_stats/<agf>/）：
#   nhood_per_slide.tsv       逐切片 × 类型对：z、经验 p、观测/期望计数
#   nhood_z_median.tsv        跨切片 z 中位 + 富集切片比例（患者级加权）
#   cooccurrence_per_slide.tsv 逐切片 × 类型对 × 距离壳：共现富集比
#   manifest.json
#
# 🔴 §3.7 硬约束（逐条落实）：
#   1. `nhood_enrichment(n_perms=1000)` 只返回 z，**不返回 p** ⇒ **禁止**把结果写成 `P<0.001`。
#      经验 p 由**本脚本自己的标签置换**算：p = (b+1)/(n+1)，b = 置换中统计量 ≥ 观测的个数。
#      同时保留 squidpy 的 z 与自算 z，**两者并报**（不一致就说明零分布口径有差，须登记）。
#   2. 距离步长固定 `{100,200,400,800}` µm（S8；1000 已删）。
#   3. **逐切片算 → 汇总**（切片间没有空间关系，坐标不可拼）。
#   4. 多重检验 BH-FDR **仅作报告**；⚠️ Fisher/FDR 会被组成性 + 大 N 灌水 ⇒ **主证据用效应量**。
#   5. `co_occurrence` **不做置换** ⇒ 不得声称做过置换。
#   6. 簇数由 §3.4 决定（读 consensus_domains.tsv），**不硬编码**。
#
# 坐标：厂商 `pxl_col/pxl_row_in_fullres`（§3.2.1 实测唯一正确写法）；µm 换算
#   UM_PER_PX = 0.2304（基于标称 55 µm spot 直径，见 §3.2.2）⇒ 实测点距 ≈ 91 µm。
#   ⚠️ 这个换算只影响"距离壳"的 µm 报数，不影响 nhood 图（图只看拓扑）。
#
# 跑法：
#   python3 10_niche/05_spatial_stats.py --agf BOTH            # 全量
#   python3 10_niche/05_spatial_stats.py --agf BOTH --smoke 2  # 只跑前 2 张

import argparse, json, os, sys, time
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
import squidpy as sq

ROOT  = "/home/eto/luad_v2"
NICHE = os.path.join(ROOT, "results/10_niche")
VISIUM= os.path.join(ROOT, "data/visium_spatial")

UM_PER_PX  = 0.2304          # 标称 55 µm spot 直径换算（§3.2.2）
N_NEIGHS   = 6               # Visium 六角栅格：每个 spot 6 个近邻
N_PERMS    = 1000            # S12
SEED       = 20261001
DIST_UM    = np.array([100.0, 200.0, 400.0, 800.0])   # S8
N_HOOD_MIN = 3               # 邻域统计至少要有这么多 spot，否则跳过该切片

def log(fmt, *a):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)

def bh_fdr(p):
    """Benjamini-Hochberg；只用于报告。"""
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0: return p
    o = np.argsort(p); ranked = p[o]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1)
    return out

def read_coords(slide, barcodes):
    po = pd.read_csv(os.path.join(VISIUM, slide, "spatial", "tissue_positions.csv"))
    po = po.set_index("barcode")
    miss = [b for b in barcodes if b not in po.index]
    if miss:
        raise SystemExit("%s：%d 个 barcode 不在 tissue_positions.csv" % (slide, len(miss)))
    sub = po.loc[barcodes]
    xy = np.c_[sub["pxl_col_in_fullres"].to_numpy(float),
               sub["pxl_row_in_fullres"].to_numpy(float)]
    xy = xy - xy.min(axis=0)              # 平移非负（避负索引陷阱；距离不变）
    return xy * UM_PER_PX                 # → µm

def statistic(labels, idx, n_lab):
    """邻域计数矩阵：类型 i 邻接类型 j 的次数（无向，两个方向都记）。"""
    li, lj = labels[idx[0]], labels[idx[1]]
    M = np.zeros((n_lab, n_lab), dtype=float)
    np.add.at(M, (li, lj), 1.0)
    np.add.at(M, (lj, li), 1.0)
    return M

def nhood_one(slide, labs, xy, seed):
    """返回 (z_squidpy, z_self, p_self, obs_count, types, n_neigh_edges)。"""
    import anndata as ad
    cats = sorted(set(labs))
    code = {c: i for i, c in enumerate(cats)}
    lab_i = np.array([code[c] for c in labs])
    n_lab = len(cats)

    A = ad.AnnData(X=sparse.csr_matrix((len(labs), 1)))
    A.obs["niche"] = pd.Categorical(labs, categories=cats)
    A.obsm["spatial"] = xy
    sq.gr.spatial_neighbors(A, coord_type="generic", n_neighs=N_NEIGHS, set_diag=False)
    conn = A.obsp["spatial_connectivities"].tocoo()
    idx = np.vstack([conn.row, conn.col])

    # —— squidpy 版（预注册点名的函数；只给 z）——
    sq.gr.nhood_enrichment(A, cluster_key="niche", n_perms=N_PERMS, seed=seed,
                           show_progress_bar=False, n_jobs=1)
    z_sq = np.asarray(A.uns["niche_nhood_enrichment"]["zscore"], dtype=float)

    # —— 自算版：同一零分布定义下给 z 与经验 p ——
    obs = statistic(lab_i, idx, n_lab)
    rng = np.random.default_rng(seed)
    perms = np.empty((N_PERMS, n_lab, n_lab))
    for b in range(N_PERMS):
        perm = rng.permutation(lab_i)
        perms[b] = statistic(perm, idx, n_lab)
    mu = perms.mean(axis=0)
    sd = perms.std(axis=0, ddof=1)
    sd_safe = np.where(sd > 0, sd, np.nan)
    z_self = (obs - mu) / sd_safe
    ge = (perms >= obs[None, :, :]).sum(axis=0)
    p_self = (ge + 1.0) / (N_PERMS + 1.0)
    return z_sq, z_self, p_self, obs, np.array(cats), idx.shape[1]

def radius_pairs(xy, r_um):
    """半径图：所有距离 ≤ r_um 的点对，返回 (i, j) 两列（无向只给一次）。
    🔴 共现必须用**半径图**：6-NN 图的边长恒 ≈ 一个点距(91 µm)，拿它去切
       {100,200,400,800} µm 的壳，各壳装的是同一批边 ⇒ 距离壳完全失效。
    用 sparse_distance_matrix（C 实现），比 query_pairs 的 Python set 快得多也省内存。"""
    t = cKDTree(xy)
    S = t.sparse_distance_matrix(t, r_um, output_type="coo_matrix")
    S = sparse.triu(S, k=1).tocoo()
    return np.column_stack([S.row, S.col]).astype(np.int64)

def cooc_one(labs, xy, cats_int, cats_str, pairs):
    """在固定 µm 壳 {100,200,400,800} 上算共现富集比（不做置换）。

    🔴 2026-10-02 修：原签名 `cooc_one(labs, xy, cats, pairs)` 收的 `cats` 是**字符串**名
      （调用处传 `cat_str`），而 `labs` 是**整数**共识类型 ⇒ `code[c]` 全部 KeyError: 1
      （09:52:36 那次 rc=1 的根因）。现把**整数类目**与**报名字符串**分开传：
      索引一律用整数类目建，报出的 `type_i/type_j` 用字符串。数值口径不变。
    """
    if pairs.size == 0: return []
    rows_i, cols_i = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(xy[rows_i] - xy[cols_i], axis=1)
    code = {c: i for i, c in enumerate(cats_int)}
    li = np.array([code[c] for c in labs])
    out = []
    for D in DIST_UM:
        keep = d <= D
        if keep.sum() == 0: continue
        M = np.zeros((len(cats_int), len(cats_int)))
        np.add.at(M, (li[rows_i[keep]], li[cols_i[keep]]), 1.0)
        np.add.at(M, (li[cols_i[keep]], li[rows_i[keep]]), 1.0)
        row = M.sum(axis=1, keepdims=True); col = M.sum(axis=0, keepdims=True)
        tot = M.sum()
        exp = row @ col / tot if tot > 0 else np.full_like(M, np.nan)
        ratio = np.where(exp > 0, M / exp, np.nan)
        for i in range(len(cats_int)):
            for j in range(i, len(cats_int)):
                out.append(dict(dist_um=float(D), type_i=cats_str[i], type_j=cats_str[j],
                                n_pair=int(M[i, j]), n_expected=float(exp[i, j]),
                                enrich_ratio=float(ratio[i, j])))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--agf", default="BOTH", choices=["BOTH", "TRUE", "FALSE"])
    ap.add_argument("--smoke", type=int, default=0)
    a = ap.parse_args()
    agfs = ["agfT", "agfF"] if a.agf == "BOTH" else (["agfT"] if a.agf == "TRUE" else ["agfF"])

    for agf in agfs:
        dfile = os.path.join(NICHE, "consensus", agf, "consensus_domains.tsv")
        if not os.path.exists(dfile):
            log("！缺 %s ⇒ 跳过（先跑 03 + 04）", dfile); continue
        KD = pd.read_csv(dfile, sep="\t")
        OUT = os.path.join(NICHE, "spatial_stats", agf)
        os.makedirs(OUT, exist_ok=True)
        slides = sorted(KD["slide"].unique())
        if a.smoke: slides = slides[:a.smoke]
        log("==== AGF=%s：%d 张切片 / %d 类 ====", agf, len(slides), KD["consensus_type"].nunique())

        dom2ct = {}
        for s, d, ct in zip(KD["slide"], KD["domain"], KD["consensus_type"]):
            dom2ct[(s, int(d))] = int(ct)

        nhood_rows, cooc_rows, skipped = [], [], []
        for k, s in enumerate(slides, 1):
            bf = os.path.join(NICHE, "banksy", s, "domains_%s.tsv.gz" % agf)
            if not os.path.exists(bf):
                skipped.append((s, "no domains file")); continue
            D = pd.read_csv(bf, sep="\t")
            sub = KD[KD["slide"] == s].iloc[0]
            D = D[(D["k_geom"] == sub["k_geom"]) & (D["lambda"] == sub["lambda"]) &
                  (D["resolution"] == sub["resolution"]) & (D["seed"] == 0)]
            if len(D) < N_HOOD_MIN:
                skipped.append((s, "too few spots")); continue
            labs = [dom2ct.get((s, int(d)), -1) for d in D["domain"]]
            keep = [i for i, v in enumerate(labs) if v >= 0]
            if len(keep) < N_HOOD_MIN:
                skipped.append((s, "domain not in consensus")); continue
            D = D.iloc[keep]; labs = np.array([labs[i] for i in keep])
            ## 🔴 2026-10-02 修：某些切片域→共识类型映射后**只剩 1 个类型** ⇒ squidpy 报
            ##    `ValueError: Expected at least 2 clusters, found 1`（10:13:50 那次 rc=1 的根因）。
            ##    单类型切片**做不了**邻域富集（无类型对它）⇒ 按既有的 `skipped` 机制**逐张上报跳过**，
            ##    不是偷偷丢样本。阈值是"能不能算"，不是判据。
            if len(set(labs)) < 2:
                skipped.append((s, "single consensus type after mapping")); continue
            bc = D["barcode"].tolist()
            xy = read_coords(s, bc)

            z_sq, z_self, p_self, obs, cats, n_edge = nhood_one(s, list(labs), xy, SEED)
            # squidpy 的类别顺序 = A.obs['niche'].categories（已按 cats 排序）⇒ 与 cats 对齐
            cat_str = [str(c) for c in cats]
            for i in range(len(cats)):
                for j in range(i, len(cats)):
                    nhood_rows.append(dict(
                        slide=s, patient=s.split("_")[1], type_i=cat_str[i], type_j=cat_str[j],
                        z_squidpy=float(z_sq[i, j]), z_self=float(z_self[i, j]),
                        p_empirical=float(p_self[i, j]), obs_count=float(obs[i, j]),
                        n_edges=n_edge, n_spot=len(labs)))
            # co-occurrence：半径图（≤800 µm 的所有点对），按 µm 壳切
            pairs = radius_pairs(xy, float(DIST_UM.max()))
            for r in cooc_one(list(labs), xy, list(cats), cat_str, pairs):
                r["slide"] = s; r["patient"] = s.split("_")[1]
                cooc_rows.append(r)
            if k % 5 == 0 or k == len(slides):
                log("  %d/%d（最近 %s，%d spot，%d 类）", k, len(slides), s, len(labs), len(cats))

        if not nhood_rows:
            log("！没有任何切片产出统计 ⇒ 跳过"); continue
        NH = pd.DataFrame(nhood_rows)
        NH["q_bh_report_only"] = bh_fdr(NH["p_empirical"].to_numpy())
        NH.to_csv(os.path.join(OUT, "nhood_per_slide.tsv"), sep="\t", index=False)

        # —— 汇总：先按患者聚合（同一患者多张切片不是独立样本，§5.3），再跨患者取中位 ——
        pat = NH.groupby(["patient", "type_i", "type_j"], as_index=False).agg(
            z=("z_squidpy", "median"), z_self=("z_self", "median"),
            p_min=("p_empirical", "min"))
        agg = pat.groupby(["type_i", "type_j"], as_index=False).agg(
            z_median_across_patients=("z", "median"),
            z_self_median=("z_self", "median"),
            p_min_across_patients=("p_min", "min"),
            n_patient=("z", "size"),
            n_patient_z_gt0=("z", lambda v: int((v > 0).sum())))
        agg["frac_patient_enriched"] = agg["n_patient_z_gt0"] / agg["n_patient"]
        agg["q_bh_report_only"] = bh_fdr(agg["p_min_across_patients"].to_numpy())  # 仅报告
        agg.sort_values(["type_i", "type_j"]).to_csv(
            os.path.join(OUT, "nhood_z_median.tsv"), sep="\t", index=False)

        pd.DataFrame(cooc_rows).to_csv(os.path.join(OUT, "cooccurrence_per_slide.tsv"),
                                       sep="\t", index=False)

        # —— z 口径一致性检查（squidpy vs 自算）——
        dz = float(np.nanmax(np.abs(NH["z_squidpy"] - NH["z_self"])))
        if skipped: log("  ⚠️ 跳过 %d 张：%s", len(skipped), skipped[:5])
        if dz > 0.5:
            log("  🔴 z 口径不一致：max|z_squidpy - z_self| = %.3f ⇒ 两套零分布不同，**须登记**", dz)
        else:
            log("  ✅ z 口径一致：max|z_squidpy - z_self| = %.3f", dz)

        man = dict(step="§3.7 Squidpy 空间统计", agf=agf, n_slide=len(slides),
                   n_slide_done=NH["slide"].nunique(), n_slide_skipped=len(skipped),
                   n_perms=N_PERMS, seed=SEED, dist_um_um=DIST_UM.tolist(),
                   um_per_px=UM_PER_PX, n_neighs=N_NEIGHS,
                   squidpy=sq.__version__,
                   coord="pxl_col/pxl_row_in_fullres，平移非负 → µm（§3.2.1/§3.2.2）",
                   p_def="经验 p = (b+1)/(n+1)，b = 置换中统计量 ≥ 观测的次数（自算）",
                   caveats=["nhood_enrichment 不返回 p ⇒ 不得写 P<0.001",
                            "co_occurrence 不做置换 ⇒ 不得声称做过置换",
                            "BH-FDR 仅作报告；主证据用效应量",
                            "跨切片汇总按患者聚类，" +
                            "同一患者多张切片不是独立样本（§5.3）"],
                   max_abs_z_gap_squidpy_vs_self=dz,
                   finished=time.strftime("%Y-%m-%d %H:%M:%S"))
        with open(os.path.join(OUT, "manifest.json"), "w") as f:
            json.dump(man, f, ensure_ascii=False, indent=2)
        log("  落盘 → %s", OUT)
    log("==== 05 结束 ====")

if __name__ == "__main__":
    main()
