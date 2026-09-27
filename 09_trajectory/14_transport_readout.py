#!/usr/bin/env python
"""
运输读数：绕开谱机制，直接用 WOT 耦合做「沿期别轴」的读数。

判读**写死在** 09_trajectory/13_transport_readout_prereg.md（已签字）。
本脚本只执行；对预注册的两处偏离在文件末尾 §偏离 里登记，并在 JSON 里带出。

核心算法（预注册 §2）：不构造任何大矩阵。
    设 u 是定义在 IAC 细胞上的标量场。令 w4 = u，逐期回代
        w3 = T_34 @ w4 ;  w2 = T_23 @ w3 ;  w1 = T_12 @ w2 ;  w0 = T_01 @ w1
    则第 i 期细胞 c 的读数 s(c) = w_i(c)（i=4 时 s = u）。
    一次回代同时得到全部五期，代价 = 4 次稠密矩阵-向量乘。

🔴 实现上踩过并已修正的一个坑（写在这里防止再犯）：
    §4.2 的置换零分布必须是「**打乱 u 后重算 s**」，不能是「打乱 s」。
    SD 对置换不变 ⇒ 打乱 s 的零分布恒等于观测值，z≡0、p≡1，检验失效。
    打乱 u 改变的是**映射**（Π 的列被打乱），SD 才会变。

跑法：python 09_trajectory/14_transport_readout.py
"""
import gc
import json
import os
import time

import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp

IN = "/home/eto/luad_v2/results/09_trajectory/wot_full/epiA_wot_input.h5ad"
TMAP = "/home/eto/luad_v2/results/09_trajectory/wot_full/tmaps"
EMB = "/home/eto/luad_v2/results/04_integration/seurat_trad/epiA/embeddings"
OUT = "/home/eto/luad_v2/results/09_trajectory/readout"

STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]
PAIRS = [(0, 1), (1, 2), (2, 3), (3, 4)]
COORDS = ["U1_DC1", "U2_nFeature", "U3_pct_mt", "U4_random"]
MAIN_COORD = "U1_DC1"
NEG_COORDS = ["U2_nFeature", "U3_pct_mt", "U4_random"]
NORMS = ["rownorm", "raw"]
RANDOM_SEED = 20260927
PASS_ALPHA = 0.01
ETA2_DOWNGRADE = 0.5
KNN_K = 20

# 偏离 D1（登记）：预注册写 N_PERM=1000。用 1000 需 ~5.5 小时；而判定阈值 PASS_ALPHA=0.01
#   ⇒ 置换 p 的最小可达值是 1/(N+1)，只要 1/(N+1) < 0.01 即 N > 99 就够。
#   故取 **200**（最小可达 p = 1/201 = 0.00498 < 0.01，判定分辨率足够），省下 80% 时间。
N_PERM = 200
# 偏离 D2（登记）：IAC 期 s ≡ u（w4 = u），打乱 u 后 SD 逐位不变 ⇒ z≡0、p≡1，该期天然退化。
#   若把它留在判定里，「U1 严格压过对照」会因 0>0 恒为假 ⇒ 判定恒 FAIL。
#   故 §4.2 的置换检验与判定**只对前四期**做，IAC 期只作报告、不入判定。
PERM_STAGES = [0, 1, 2, 3]


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def rss_gb():
    with open("/proc/self/status") as fh:
        for ln in fh:
            if ln.startswith("VmRSS:"):
                return int(ln.split()[1]) / 1024 / 1024
    return 0.0


def tmap_path(i, j):
    return os.path.join(TMAP, f"tmaps_{float(i)}_{float(j)}.h5ad")


def read_tmap(i, j, want_index=False):
    with h5py.File(tmap_path(i, j), "r") as h:
        T = h["X"][:]
        if not want_index:
            return T
        oi, vi = h["obs"].attrs["_index"], h["var"].attrs["_index"]
        src = [x.decode() if isinstance(x, bytes) else x for x in h["obs"][oi][:]]
        tgt = [x.decode() if isinstance(x, bytes) else x for x in h["var"][vi][:]]
    return T, src, tgt


def pearson_rows(A, B):
    A = A - A.mean(1, keepdims=True)
    B = B - B.mean(1, keepdims=True)
    na = np.linalg.norm(A, axis=1)
    nb = np.linalg.norm(B, axis=1)
    ok = (na > 0) & (nb > 0)
    r = np.full(A.shape[0], np.nan)
    r[ok] = (A[ok] * B[ok]).sum(1) / (na[ok] * nb[ok])
    return r


def backsub(u, coups, rownorm_sums):
    """从 IAC 场 u 逐期回代，返回 {4: u, 3: w3, 2: w2, 1: w1, 0: w0}。"""
    blocks = {4: u}
    w = u
    for (i, j) in reversed(PAIRS):
        w = coups[(i, j)] @ w
        if rownorm_sums is not None:
            rs = rownorm_sums[(i, j)]
            w = w / rs
        blocks[i] = w
    return blocks


def main():
    os.makedirs(OUT, exist_ok=True)
    rep = {"prereg": "09_trajectory/13_transport_readout_prereg.md",
           "date": time.strftime("%Y-%m-%d %H:%M:%S"),
           "note": "机械执行预注册；对预注册的两处偏离见 deviations 字段"}
    rep["deviations"] = [
        {"id": "D1", "prereg_says": "N_PERM=1000", "actual": N_PERM,
         "reason": "判定阈值 α=0.01 ⇒ 置换检验只需 N>99 才有足够分辨率（最小 p=1/(N+1)）；"
                   "1000 次需 ~5.5 小时，200 次的最小 p=0.00498<0.01，够用且省 80% 时间。"},
        {"id": "D2", "prereg_says": "对五个期别都做置换检验并全部入判定",
         "actual": "只对前四期做置换检验并入判定；IAC 期只报告",
         "reason": "IAC 期 s≡u ⇒ 打乱 u 后 SD 逐位不变 ⇒ z≡0、p≡1，是构造性退化；"
                   "留在判定里会让「严格压过对照」恒为 0>0=假，判定恒 FAIL。"},
    ]

    # ─────────────── 读入 ───────────────
    log("读 obs …")
    with h5py.File(IN, "r") as h:
        o = h["obs"]
        names = np.array([x.decode() if isinstance(x, bytes) else x
                          for x in o[o.attrs["_index"]][:]])
        def cat(k):
            g = o[k]
            cats = np.array([x.decode() if isinstance(x, bytes) else x
                             for x in g["categories"][:]])
            return cats[g["codes"][:]]
        stage, pid = cat("stage"), cat("patient_id")
        nfeat = o["nFeature"][:].astype(np.float64)
        pmt = o["pct_mt"][:].astype(np.float64)
        n_var = len(h["var"][h["var"].attrs["_index"]])
        genes = np.array([x.decode() if isinstance(x, bytes) else x
                          for x in h["var"][h["var"].attrs["_index"]][:]])
    n_cells = len(names)
    si = np.array([STAGES.index(s) for s in stage])
    sizes = [int((si == i).sum()) for i in range(5)]
    assert sizes == [32251, 11891, 28525, 4583, 56134], f"期别计数与预注册不符：{sizes}"
    upats = np.unique(pid)
    pidx = pd.Index(upats).get_indexer(pid)
    log(f"  {n_cells} 细胞 × {n_var} 基因；各期 {dict(zip(STAGES, sizes))}；"
        f"{len(upats)} 患者")
    rep["input"] = {"h5ad": IN, "n_cells": n_cells, "n_vars": int(n_var),
                    "stage_counts": dict(zip(STAGES, sizes)), "n_patients": int(len(upats))}

    log("读冻结 harmony 嵌入 …")
    cells = [x.strip() for x in open(os.path.join(EMB, "cells.txt")) if x.strip()]
    H = np.fromfile(os.path.join(EMB, "harmony_f32.bin"), dtype=np.float32)
    pos = pd.Index(cells).get_indexer(names)
    assert (pos >= 0).all(), "有条码在冻结嵌入里找不到"
    harmony = H.reshape(len(cells), -1)[pos].copy()
    del H, pos, cells
    log(f"  harmony {harmony.shape}（0 缺失）")
    rep["embedding"] = {"path": EMB, "dim": int(harmony.shape[1])}

    # ─────────────── 坐标阶梯 ───────────────
    log("算 U1 = DC1（scanpy neighbors k=20 + diffmap）…")
    t = time.time()
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(X=sp.csr_matrix((n_cells, 1), dtype=np.float32))
    a.obs_names = names
    a.obsm["X_emb"] = harmony
    sc.pp.neighbors(a, n_neighbors=KNN_K, use_rep="X_emb")
    sc.tl.diffmap(a, n_comps=5)
    dc1 = np.asarray(a.obsm["X_diffmap"][:, 1], dtype=np.float64)
    log(f"  diffmap 用了 {time.time()-t:.1f}s；DC1 范围 {dc1.min():.3f}…{dc1.max():.3f}")
    del a; gc.collect()

    rng = np.random.default_rng(RANDOM_SEED)
    fields = {"U1_DC1": dc1, "U2_nFeature": nfeat, "U3_pct_mt": pmt,
              "U4_random": rng.standard_normal(n_cells)}
    for k in COORDS:                       # 四档统一定符号：与期别序号正相关
        if np.corrcoef(fields[k], si)[0, 1] < 0:
            fields[k] = -fields[k]
    rep["coord_sign_convention"] = "四档统一：使其与期别序号正相关（保证跨档可比）"
    rep["coord_stage_corr"] = {k: round(float(np.corrcoef(fields[k], si)[0, 1]), 4)
                               for k in COORDS}
    rep["coord_depth_corr"] = {k: round(float(np.corrcoef(fields[k], nfeat)[0, 1]), 4)
                               for k in COORDS}
    rep["coord_patient_eta2"] = {}
    for k in COORDS:
        gm = fields[k].mean()
        ssb = sum((pidx == q).sum() * (fields[k][pidx == q].mean() - gm) ** 2
                  for q in range(len(upats)))
        sst = ((fields[k] - gm) ** 2).sum()
        rep["coord_patient_eta2"][k] = round(float(ssb / sst), 4)
    log(f"  与期别序号相关：{rep['coord_stage_corr']}")
    log(f"  与深度相关：    {rep['coord_depth_corr']}")
    log(f"  患者 η²：       {rep['coord_patient_eta2']}")

    # ─────────────── 条码对齐（全部索引的地基）+ 耦合常驻内存 ───────────────
    log("核对耦合条码顺序并把四张耦合读进内存 …")
    coups, rownorm_sums, align = {}, {}, {}
    for (i, j) in PAIRS:
        T, src, tgt = read_tmap(i, j, want_index=True)
        ok_s = bool((np.array(src) == names[si == i]).all())
        ok_t = bool((np.array(tgt) == names[si == j]).all())
        align[f"{i}->{j}"] = {"shape": list(T.shape), "src_matches_adata_order": ok_s,
                              "tgt_matches_adata_order": ok_t}
        assert ok_s and ok_t, f"耦合 {i}->{j} 的条码顺序与 adata 不一致，停"
        rs = T.sum(1); rs[rs == 0] = 1.0
        rownorm_sums[(i, j)] = rs
        coups[(i, j)] = T
        log(f"  T_{i}{j} {T.shape} 对齐={ok_s}/{ok_t}  rss={rss_gb():.1f}GB")
    rep["tmap_alignment"] = align
    log(f"  四张耦合常驻内存，rss={rss_gb():.1f}GB")

    iac_idx = np.where(si == 4)[0]

    # ─────────────── §2 观测读数 ───────────────
    log("§2 观测读数 …")
    s_obs = {}
    for norm in NORMS:
        ssums = rownorm_sums if norm == "rownorm" else None
        for k in COORDS:
            bl = backsub(fields[k][iac_idx].copy(), coups, ssums)
            s = np.empty(n_cells)
            for i in range(5):
                s[si == i] = bl[i]
            s_obs[(k, norm)] = s
            del bl
    log("  完成")

    # ─────────────── §4.2 置换检验（打乱 u，重算 s）───────────────
    log(f"§4.2 置换检验：{N_PERM} 次 × {len(COORDS)} 档 × {len(NORMS)} 归一化 …")
    perms = [rng.permutation(len(iac_idx)) for _ in range(N_PERM)]
    ladder_rows = []
    for norm in NORMS:
        ssums = rownorm_sums if norm == "rownorm" else None
        for k in COORDS:
            t0 = time.time()
            sd_null = {i: np.empty(N_PERM) for i in PERM_STAGES}
            u0 = fields[k][iac_idx]
            for b, pm in enumerate(perms):
                bl = backsub(u0[pm].copy(), coups, ssums)
                for i in PERM_STAGES:
                    sd_null[i][b] = bl[i].std()
                del bl
                if (b + 1) % 50 == 0:
                    log(f"    {k}/{norm} 置换 {b+1}/{N_PERM}  {time.time()-t0:.0f}s")
            s = s_obs[(k, norm)]
            for i in range(5):
                g = s[si == i]
                sd_obs = float(g.std())
                if i in PERM_STAGES:
                    null = sd_null[i]
                    mu, sd = float(null.mean()), float(null.std())
                    z = (sd_obs - mu) / sd if sd > 0 else np.nan
                    p = float((null >= sd_obs).mean())
                else:
                    mu = sd = z = p = np.nan
                pp = pid[si == i]
                gm = g.mean()
                ssb = sum((pp == q).sum() * (g[pp == q].mean() - gm) ** 2
                          for q in np.unique(pp))
                sst = ((g - gm) ** 2).sum()
                ladder_rows.append(dict(
                    coord=k, norm=norm, stage=STAGES[i], n=int(g.size), sd=sd_obs,
                    perm_sd_mean=mu, perm_sd_sd=sd, z=z, p_perm=p,
                    eta2_patient=float(ssb / sst) if sst > 0 else np.nan,
                    median=float(np.median(g)),
                    in_verdict=bool(i in PERM_STAGES)))
            log(f"  {k}/{norm} 置换完成 {time.time()-t0:.0f}s")
    ladder = pd.DataFrame(ladder_rows)
    ladder.to_csv(f"{OUT}/ladder_by_coordinate.tsv", sep="\t", index=False)

    med = {}
    for k in COORDS:
        for norm in NORMS:
            m = [float(np.median(s_obs[(k, norm)][si == i])) for i in range(5)]
            med[f"{k}|{norm}"] = {"median_by_stage": [round(x, 6) for x in m],
                                  "monotone_increasing": bool(np.all(np.diff(m) > 0))}
    rep["monotonicity"] = med

    def z_of(c, norm, i):
        d = ladder[(ladder.coord == c) & (ladder.norm == norm) & (ladder.stage == STAGES[i])]
        return float(d["z"].iloc[0])

    verdict = {}
    for norm in NORMS:
        z_main = [z_of(MAIN_COORD, norm, i) for i in PERM_STAGES]
        z_neg = {c: [z_of(c, norm, i) for i in PERM_STAGES] for c in NEG_COORDS}
        sub = ladder[(ladder.coord == MAIN_COORD) & (ladder.norm == norm) & ladder.in_verdict]
        eta2_max, p_max = float(sub["eta2_patient"].max()), float(sub["p_perm"].max())
        # 判定：主坐标必须**在全部四期上**都严格压过每一个对照档（逐期比较）
        beats = {c: bool(np.all(np.array(z_main) > np.array(z_neg[c]))) for c in NEG_COORDS}
        v = {"z_main_by_stage": [round(x, 3) for x in z_main],
             "z_neg_by_stage": {c: [round(x, 3) for x in z_neg[c]] for c in NEG_COORDS},
             "max_eta2_patient": round(eta2_max, 4), "max_p_perm": round(p_max, 6),
             "median_by_stage": med[f"{MAIN_COORD}|{norm}"]["median_by_stage"],
             "monotone": med[f"{MAIN_COORD}|{norm}"]["monotone_increasing"],
             "strictly_beats_negatives_by_stage": beats}
        if not all(beats.values()):
            v["verdict"] = "FAIL：主坐标没能逐期严格压过全部阴性对照 ⇒ 运输没把连贯的东西推过去"
        elif p_max >= PASS_ALPHA:
            v["verdict"] = f"FAIL：期内离散度不显著高于置换零分布（max p={p_max:.4f}）"
        elif not v["monotone"]:
            v["verdict"] = "降级：离散度与对照可分，但 median(s|stage) 不单调"
        elif eta2_max >= ETA2_DOWNGRADE:
            v["verdict"] = "降级：患者效应过大（η²≥0.5），须按患者分层重报"
        else:
            v["verdict"] = "PASS：主坐标期内离散度显著、逐期压过三档对照、单调、患者效应可控"
        verdict[norm] = v
    rep["ladder_verdict"] = verdict
    for norm in NORMS:
        log(f"  [{norm}] {verdict[norm]['verdict']}")

    # ─────────────── §4.1 地板 ───────────────
    log("§4.1 地板检验 …")
    with h5py.File(IN, "r") as h:
        X = sp.csr_matrix((h["X"]["data"][:], h["X"]["indices"][:], h["X"]["indptr"][:]),
                          shape=(n_cells, n_var)).tocsr()
    log(f"  X nnz={X.nnz:,}  rss={rss_gb():.1f}GB")
    floor_rows = []
    rngp = np.random.default_rng(RANDOM_SEED + 2)
    for (i, j) in PAIRS:
        t = time.time()
        T = coups[(i, j)]
        cs = T.sum(0); cs[cs == 0] = 1.0
        Xi = X[si == i].toarray().astype(np.float32)
        Xj = X[si == j].toarray().astype(np.float32)
        Xh = (T.T @ Xi) / cs[:, None].astype(np.float32)
        perm = rngp.permutation(T.shape[0])
        Ts = T[perm]
        cs2 = Ts.sum(0); cs2[cs2 == 0] = 1.0
        Xn = (Ts.T @ Xi) / cs2[:, None].astype(np.float32)
        # 独立耦合零模型：x̂(d,g) = Σ_c w_c x_i(c,g) / Σ_c w_c，与 d 无关（w_c = 源行和）
        # ⇒ 每个目标细胞都是同一个「加权平均画像」，故铺 n_j 行（= T.shape[1] 是目标数）。
        a_c = T.sum(1) / T.sum()
        Xind = np.tile((a_c @ Xi) / a_c.sum(), (T.shape[1], 1))
        for nm, Xp in [("real", Xh), ("null_rowperm", Xn), ("null_indep", Xind)]:
            floor_rows.append(dict(
                pair=f"{STAGES[i]}->{STAGES[j]}", model=nm,
                median_r_across_targets=float(np.nanmedian(pearson_rows(Xp, Xj))),
                median_r_across_genes=float(np.nanmedian(pearson_rows(Xp.T, Xj.T))),
                n_targets=int(T.shape[1])))
        log(f"  {STAGES[i]}->{STAGES[j]} {time.time()-t:.0f}s  "
            f"真={floor_rows[-3]['median_r_across_targets']:.4f} "
            f"行置换={floor_rows[-2]['median_r_across_targets']:.4f} "
            f"独立={floor_rows[-1]['median_r_across_targets']:.4f}")
        del Ts, Xi, Xj, Xh, Xn, Xind; gc.collect()
    floor = pd.DataFrame(floor_rows)
    floor.to_csv(f"{OUT}/floor_test.tsv", sep="\t", index=False)
    piv = floor.pivot(index="pair", columns="model", values="median_r_across_targets")
    rep["floor_test"] = {
        "median_r_across_targets": piv.round(4).to_dict(),
        "median_r_across_genes": floor.pivot(index="pair", columns="model",
                                             values="median_r_across_genes").round(4).to_dict(),
        "real_beats_null_rowperm_all_pairs": bool((piv["real"] > piv["null_rowperm"]).all()),
        "real_beats_null_indep_all_pairs": bool((piv["real"] > piv["null_indep"]).all())}
    log("  地板表已写")

    # ─────────────── §4.3 基因趋势 ───────────────
    log("§4.3 基因趋势（自举 + 留一患者）…")
    G = int(n_var)
    dfreq = np.zeros((len(PAIRS), G))
    num = np.zeros((len(PAIRS), len(upats), G))
    den = np.zeros((len(PAIRS), len(upats)))
    base = np.zeros((len(PAIRS), len(upats), G))
    nbase = np.zeros((len(PAIRS), len(upats)))
    for pi, (i, j) in enumerate(PAIRS):
        T = coups[(i, j)]
        cs = T.sum(0); cs[cs == 0] = 1.0
        Xi = X[si == i].toarray().astype(np.float32)
        Xh = (T.T @ Xi) / cs[:, None].astype(np.float32)
        dfreq[pi] = Xh.mean(0).astype(np.float64) - Xi.mean(0).astype(np.float64)
        pt, ps = pidx[si == j], pidx[si == i]
        for q in range(len(upats)):
            m = pt == q
            num[pi, q] = Xh[m].sum(0).astype(np.float64); den[pi, q] = m.sum()
            m2 = ps == q
            base[pi, q] = Xi[m2].sum(0).astype(np.float64); nbase[pi, q] = m2.sum()
        del Xi, Xh; gc.collect()
        log(f"  pair {i}->{j} 的 Δ 已算")
    del X, coups, rownorm_sums; gc.collect()

    N_BOOT = 1000
    rngb = np.random.default_rng(RANDOM_SEED + 3)
    boot = np.empty((N_BOOT, len(PAIRS), G))
    for b in range(N_BOOT):
        pick = rngb.integers(0, len(upats), len(upats))
        for pi in range(len(PAIRS)):
            boot[b, pi] = (num[pi][pick].sum(0) / max(den[pi][pick].sum(), 1)
                           - base[pi][pick].sum(0) / max(nbase[pi][pick].sum(), 1))
        if (b + 1) % 250 == 0:
            log(f"    自举 {b+1}/{N_BOOT}")
    lo = np.percentile(boot, 2.5, axis=0)
    hi = np.percentile(boot, 97.5, axis=0)
    del boot; gc.collect()

    loo_flip = np.zeros((len(PAIRS), G), dtype=bool)
    for q in range(len(upats)):
        keep = np.ones(len(upats), bool); keep[q] = False
        for pi in range(len(PAIRS)):
            d = (num[pi][keep].sum(0) / max(den[pi][keep].sum(), 1)
                 - base[pi][keep].sum(0) / max(nbase[pi][keep].sum(), 1))
            loo_flip[pi] |= (np.sign(d) != np.sign(dfreq[pi]))
    same_sign = np.all(np.stack([np.sign(dfreq[p]) == np.sign(dfreq[0])
                                 for p in range(len(PAIRS))]), axis=0)
    beyond = np.all((dfreq <= lo) | (dfreq >= hi), axis=0)
    loo_ok = ~loo_flip.any(0)
    cand = same_sign & beyond & loo_ok

    gt = pd.DataFrame({"gene": genes[:G]})
    for pi, (i, j) in enumerate(PAIRS):
        gt[f"delta_{STAGES[i]}_{STAGES[j]}"] = dfreq[pi]
        gt[f"boot_lo_{STAGES[i]}_{STAGES[j]}"] = lo[pi]
        gt[f"boot_hi_{STAGES[i]}_{STAGES[j]}"] = hi[pi]
    gt["same_sign"] = same_sign
    gt["beyond_boot_ci"] = beyond
    gt["loo_stable"] = loo_ok
    gt["candidate"] = cand
    gt.sort_values(f"delta_{STAGES[3]}_{STAGES[4]}", key=np.abs, ascending=False) \
      .to_csv(f"{OUT}/candidate_genes.tsv", sep="\t", index=False)
    rep["gene_trends"] = {
        "n_genes": G, "n_same_sign": int(same_sign.sum()),
        "n_beyond_boot": int(beyond.sum()), "n_loo_stable": int(loo_ok.sum()),
        "n_candidate": int(cand.sum()),
        "top15_candidates_by_abs_last_delta": gt.loc[cand].reindex(
            gt.loc[cand, f"delta_{STAGES[3]}_{STAGES[4]}"].abs()
              .sort_values(ascending=False).index)["gene"].head(15).tolist()}
    log(f"  候选基因 {int(cand.sum())} 个（同号 {int(same_sign.sum())}，"
        f"越自举区间 {int(beyond.sum())}，留一稳定 {int(loo_ok.sum())}）")

    rep["peak_rss_gb"] = round(rss_gb(), 2)
    with open(f"{OUT}/readout_summary.json", "w") as fh:
        json.dump(rep, fh, ensure_ascii=False, indent=1, default=str)
    log(f"写到 {OUT}/readout_summary.json；峰值 RSS {rss_gb():.2f}GB")

    print("\n================ 摘要 ================")
    for norm in NORMS:
        print(f"[{norm}] {verdict[norm]['verdict']}")
    print(f"地板：真 > 行置换 全对？ {rep['floor_test']['real_beats_null_rowperm_all_pairs']}")
    print(f"      真 > 独立   全对？ {rep['floor_test']['real_beats_null_indep_all_pairs']}")
    print(f"候选基因：{rep['gene_trends']['n_candidate']} / {G}")
    return 0


if __name__ == "__main__":
    import traceback
    try:
        raise SystemExit(main())
    except BaseException:
        traceback.print_exc()
        raise
