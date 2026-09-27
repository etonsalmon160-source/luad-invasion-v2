#!/usr/bin/env python
"""
归一化口径的决定性检验：三种口径下，WOT 代价空间的第一主成分是不是「测序深度」？

为什么这是决定性的：
    wot/ot/ot_model.py:290  p0_x, p1_x, pca, mean = wot.ot.compute_pca(p0.X, p1.X, local_pca)
    wot/ot/ot_model.py:240  cost = pairwise_distances(p0_x, p1_x) / median(cost)
    ⇒ 代价 = 两期细胞**合起来做一次 30 维 PCA** 之后，那个 PC 空间里的欧氏距离。
    WOT 内部**没有** log、没有 scale、没有按方差挑基因（gene_filter 默认 None）。
    所以「用哪种归一化」直接决定那个 PC 空间长什么样。

    若原始计数下 PC1 就是库大小，而库大小在队列里与期别相关（已知：LUAD 中位深度是前驱的 2.22 倍），
    那么运输轴会主要沿着深度走 —— 等于把深度这个混杂编码进了「发育轴」。

只量相关性，不做任何生物学结论。不落盘结果，只出图。
"""
import os
import time

import numpy as np

H5AD = "/home/eto/luad_v2/results/02_expression/gse308103_counts_paperqc.h5ad"
EPI_BC = "/home/eto/luad_v2/results/05_annotation/epiA_nocontam_subset_barcodes.txt"
FIG = "/home/eto/luad_v2/results/09_trajectory/smoke/norm_depth_confound.png"

STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]
N_PER_STAGE = 3000
N_PC = 30
SEED = 42


def rss_gb():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return round(int(line.split()[1]) / 1024 / 1024, 2)
    return None


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] rss={rss_gb()}GB {m}", flush=True)


def make_calibers(a):
    """返回三种口径的密集表达子集（同一个细胞抽样，只换标度）。"""
    import scanpy as sc

    rng = np.random.default_rng(SEED)
    idx = []
    for st in STAGES:
        pool = np.where((a.obs["stage"].astype(str) == st).values)[0]
        idx.append(rng.choice(pool, size=min(N_PER_STAGE, len(pool)), replace=False))
    idx = np.concatenate(idx)
    sub = a[idx].copy()
    stage = sub.obs["stage"].astype(str).values
    ncount = sub.obs["nCount"].values.astype(float)

    out = {}
    out["A raw counts"] = np.asarray(sub.X.todense(), dtype=np.float64)

    b = sub.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    out["B CP10K"] = np.asarray(b.X.todense(), dtype=np.float64)

    sc.pp.log1p(b)
    out["C CP10K+log1p"] = np.asarray(b.X.todense(), dtype=np.float64)
    return out, stage, ncount


def main():
    import anndata as ad
    from sklearn.decomposition import PCA
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    log("读入 h5ad …")
    a = ad.read_h5ad(H5AD)
    want = {ln.strip() for ln in open(EPI_BC) if ln.strip()}
    a = a[a.obs_names.isin(want)].copy()
    log(f"上皮子集 {a.n_obs} × {a.n_vars}")

    cals, stage, ncount = make_calibers(a)
    log("三种口径的密集矩阵就绪")

    # 逐相邻期别对做 PCA（复刻 WOT：两期合起来拟合），看 PC_i 与 nCount 的相关
    res = {}          # caliber -> (n_pairs, n_pc) |corr|
    pc1_scatter = {}  # caliber -> (pc1, ncount, stage) 用最大的一对
    for name, X in cals.items():
        cors, biggest = [], None
        for t0 in range(4):
            m = np.isin(stage, [STAGES[t0], STAGES[t0 + 1]])
            Xp, np_p = X[m], ncount[m]
            Xc = Xp - Xp.mean(axis=0)                 # PCA 自带中心化，但显式写清
            p = PCA(n_components=N_PC, svd_solver="randomized",
                    random_state=0).fit(Xc)
            Z = p.transform(Xc)
            c = [abs(np.corrcoef(Z[:, i], np_p)[0, 1]) for i in range(5)]
            cors.append(c)
            if biggest is None or Xp.shape[0] > biggest[0]:
                biggest = (Xp.shape[0], Z[:, 0].copy(), np_p.copy(), stage[m].copy())
        res[name] = np.array(cors)
        pc1_scatter[name] = biggest
        log(f"  {name:14s} |corr(PC1,nCount)| 四对平均 = {res[name][:, 0].mean():.3f}")

    fig, ax = plt.subplots(2, 2, figsize=(15, 10))

    # 面板 1：头条数字 —— 各口径下前 5 个 PC 与深度的相关
    w, xs = 0.26, np.arange(5)
    colors = {"A raw counts": "#888888", "B CP10K": "#3b6ea5", "C CP10K+log1p": "#c0392b"}
    for i, (name, M) in enumerate(res.items()):
        ax[0, 0].bar(xs + (i - 1) * w, M.mean(axis=0), width=w,
                     label=name, color=colors[name])
    ax[0, 0].set_xticks(xs)
    ax[0, 0].set_xticklabels([f"PC{i+1}" for i in range(5)])
    ax[0, 0].set_ylabel("|corr( PC_i , total counts )|   mean over 4 stage-pairs")
    ax[0, 0].set_ylim(0, 1.02)
    ax[0, 0].axhline(0.5, color="#444444", ls=":", lw=1)
    ax[0, 0].legend(fontsize=9)
    ax[0, 0].set_title(
        "1. Is the OT cost space just a sequencing-depth axis?\n"
        "WOT's cost = Euclidean distance in a 30-dim PCA of the two stages")

    # 面板 2 & 3：PC1 vs 深度（最大的一对），两端口径对照
    for axx, name in [(ax[0, 1], "A raw counts"), (ax[1, 0], "C CP10K+log1p")]:
        n_cells, pc1, nc, st = pc1_scatter[name]
        for s, col in zip(STAGES, ["#8c8c8c", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]):
            k = st == s
            axx.scatter(nc[k], pc1[k], s=3, alpha=0.35, color=col, label=s)
        r = np.corrcoef(pc1, nc)[0, 1]
        axx.set_xlabel("total counts per cell")
        axx.set_ylabel("PC1 of the stage-pair PCA")
        axx.set_title(f"2/3. {name}:  PC1 vs depth,  r = {r:+.3f}\n"
                      f"(largest pair, {n_cells} cells)", fontsize=10)
        axx.legend(fontsize=8, markerscale=4)
        axx.set_xscale("log")

    # 面板 4：结论与读数
    ax[1, 1].axis("off")
    r1 = {k: v[:, 0].mean() for k, v in res.items()}
    txt = (
        "How to read this\n"
        "----------------\n"
        "WOT builds its cost from a joint 30-dim PCA of the two\n"
        "stages, then Euclidean distance there (ot_model.py:290,240).\n"
        "It does NOT log, scale, or variance-filter internally.\n"
        "So the caliber decides what that PCA's leading axes are.\n\n"
        "|corr(PC1, total counts)|, averaged over the 4 stage-pairs:\n"
        f"   A raw counts    : {r1['A raw counts']:.3f}\n"
        f"   B CP10K         : {r1['B CP10K']:.3f}\n"
        f"   C CP10K + log1p : {r1['C CP10K+log1p']:.3f}\n\n"
        "Why this matters HERE specifically:\n"
        "  depth is not random in this cohort. The RCTD line already\n"
        "  measured LUAD median depth = 2.22x the precursor stages.\n"
        "  So a cost axis that is mostly depth is an axis that is\n"
        "  partly stage -- i.e. the transport would partly carry\n"
        "  the confounder we are forbidden from encoding.\n\n"
        "This test does NOT pick a winner between B and C -- it\n"
        "tests A against the other two. See the other figure\n"
        "(norm_caliber_choice.png) for what B vs C changes.\n\n"
        "Caveat: this is a 3000-cells-per-stage subsample and a\n"
        "per-pair PCA fitted here, not WOT's own fit. Magnitudes\n"
        "will differ slightly from the real run; the ORDERING is\n"
        "the claim."
    )
    ax[1, 1].text(0.0, 0.98, txt, va="top", ha="left", fontsize=10,
                  family="monospace", linespacing=1.4)

    fig.suptitle("Does the normalization caliber turn WOT's cost axis into a depth axis? "
                 f"(epithelial, {N_PER_STAGE}/stage subsample)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG, dpi=115)
    log(f"图写到 {FIG}")
    for k, v in res.items():
        print(f"{k:14s} PC1..5 |corr| 四对平均: " +
              " ".join(f"{x:.3f}" for x in v.mean(axis=0)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
