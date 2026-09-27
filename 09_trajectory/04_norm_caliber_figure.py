#!/usr/bin/env python
"""
归一化口径的决策图 —— WOT 的代价矩阵直接建在表达值上，所以「用哪种归一化」会实打实
改变运输耦合。这里只做一件事：把三种候选口径下**数据的形状**画出来，让口径可看图决定。

不给任何生物学结论。不写任何结果文件，只出一张图。
"""
import os
import time

import numpy as np

H5AD = "/home/eto/luad_v2/results/02_expression/gse308103_counts_paperqc.h5ad"
EPI_BC = "/home/eto/luad_v2/results/05_annotation/epiA_nocontam_subset_barcodes.txt"
FIG = "/home/eto/luad_v2/results/09_trajectory/smoke/norm_caliber_choice.png"


def rss_gb():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return round(int(line.split()[1]) / 1024 / 1024, 2)
    return None


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] rss={rss_gb()}GB {m}", flush=True)


def main():
    import anndata as ad
    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    log("读入 h5ad …")
    a = ad.read_h5ad(H5AD)
    want = {ln.strip() for ln in open(EPI_BC) if ln.strip()}
    a = a[a.obs_names.isin(want)].copy()
    log(f"上皮子集 {a.n_obs} × {a.n_vars}")

    ncount = a.obs["nCount"].values.astype(float)

    # --- 口径 A：原始计数（交给 WOT 内部处理）---
    mean_raw = np.asarray(a.X.mean(axis=0)).ravel()

    # --- 口径 B：CP10K（每细胞缩到 1e4），不取对数 ---
    b = a.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    mean_cpm = np.asarray(b.X.mean(axis=0)).ravel()

    # --- 口径 C：CP10K + log1p（探针同款，也是 Seurat LogNormalize 的等价式）---
    sc.pp.log1p(b)
    mean_log = np.asarray(b.X.mean(axis=0)).ravel()
    log("三种口径的基因均值都算完了")

    names = np.asarray(a.var_names)
    fig, ax = plt.subplots(2, 2, figsize=(15, 10))

    # 面板 1：每细胞总计数（归一化去掉的就是这个变异）
    ax[0, 0].hist(ncount, bins=120, color="#3b6ea5")
    ax[0, 0].set_xlabel("total counts per cell (raw)")
    ax[0, 0].set_ylabel("n cells")
    ax[0, 0].set_title(
        f"1. Library size spread that normalization removes\n"
        f"median={np.median(ncount):.0f}  min={ncount.min():.0f}  "
        f"max={ncount.max():.0f}  p99/p1={np.percentile(ncount,99)/np.percentile(ncount,1):.1f}x")
    ax[0, 0].set_yscale("log")

    # 面板 2：基因均值分布（对数尺度）—— 看 log1p 把动态范围压了多少
    for v, lab, c in [(mean_raw, "A raw counts", "#888888"),
                      (mean_cpm, "B CP10K (no log)", "#3b6ea5"),
                      (mean_log, "C CP10K + log1p", "#c0392b")]:
        x = np.log10(np.clip(v, 1e-6, None))
        ax[0, 1].hist(x, bins=120, histtype="step", lw=2, label=lab, color=c)
    ax[0, 1].set_xlabel("log10( mean expression per gene )")
    ax[0, 1].set_ylabel("n genes")
    ax[0, 1].legend(fontsize=9)
    ax[0, 1].set_title(
        "2. Gene-mean distribution per caliber\n"
        "log1p squeezes the range; raw counts are dominated by a few huge genes")

    # 面板 3：榜位迁移图 —— 一个基因要在三种口径下都排前 20，才说明「主角」稳定。
    # （初版把三口径的均值画在同一根轴上比高矮，量纲差 1000 倍 ⇒ 红柱不可见，已废弃重画。）
    def rank_of(v):
        order = np.argsort(-v)
        r = np.empty(len(v), dtype=int)
        r[order] = np.arange(1, len(v) + 1)
        return r

    rA, rB, rC = rank_of(mean_raw), rank_of(mean_cpm), rank_of(mean_log)
    K, CLIP = 20, 45
    idxA = set(np.argsort(-mean_raw)[:K])
    idxB = set(np.argsort(-mean_cpm)[:K])
    idxC = set(np.argsort(-mean_log)[:K])
    toplog = list(names[np.argsort(-mean_log)[:K]])
    union_idx = sorted(idxA | idxB | idxC, key=lambda j: rC[j])

    for j in union_idx:
        ys = [min(rA[j], CLIP), min(rB[j], CLIP), min(rC[j], CLIP)]
        stable = (rA[j] <= K) and (rB[j] <= K) and (rC[j] <= K)
        ax[1, 0].plot([0, 1, 2], ys, "-o", ms=4,
                      lw=2.0 if stable else 1.0,
                      color="#1f6f3f" if stable else "#c0392b",
                      alpha=0.95 if stable else 0.8, zorder=3 if stable else 2)
        ax[1, 0].text(-0.08, ys[0], names[j], ha="right", va="center", fontsize=7)
        ax[1, 0].text(2.08, ys[2], names[j], ha="left", va="center", fontsize=7,
                      fontweight="bold")
    ax[1, 0].set_xticks([0, 1, 2])
    ax[1, 0].set_xticklabels(["A\nraw counts", "B\nCP10K", "C\nCP10K+log1p"],
                             fontsize=9)
    ax[1, 0].set_xlim(-1.5, 3.5)
    ax[1, 0].set_ylim(CLIP + 5, 0.2)
    ax[1, 0].axhspan(K + 0.5, CLIP + 5, color="#000000", alpha=0.06, zorder=0)
    ax[1, 0].axhline(K + 0.5, color="#444444", ls=":", lw=1, zorder=1)
    ax[1, 0].set_ylabel(f"gene rank by mean expression  (1 = top, clipped at {CLIP})",
                        fontsize=9)
    ax[1, 0].set_title("3. Rank migration across calibers\n"
                       "green = top-20 in all three (stable); red = enters/leaves top-20",
                       fontsize=10)

    # 面板 4：榜位重合度 —— 直接回答「换口径会不会换掉主角」
    ax[1, 1].axis("off")
    inter_raw = len(idxC & idxA)
    inter_cpm = len(idxC & idxB)
    txt = (
        "Top-20 overlap with caliber C (CP10K + log1p)\n\n"
        f"  A raw counts   : {inter_raw}/20 genes in common\n"
        f"  B CP10K        : {inter_cpm}/20 genes in common\n\n"
        "B and C share the same per-cell scaling, so their\n"
        "difference is ONLY whether values are logged.\n"
        "That is the one real choice here.\n\n"
        "What the swap looks like (panel 3):\n"
        "  stable top-20 = surfactant / secretory program\n"
        "      SFTPB SFTPC SFTPA1 SLC34A2 MUC1 NAPSA CTSH\n"
        "      RNASE1 ...  -> these drive the OT cost either way\n"
        "  enters only WITH log1p (broadly-expressed, stress/\n"
        "  housekeeping): S100A6 SAT1 FOS JUND LMNA MCL1 MT-CO3\n"
        "  enters only WITHOUT log (highly secreted, club cells):\n"
        "      SCGB1A1 SCGB3A2 SCGB3A1 PIGR SLPI SFTPD AGER\n\n"
        "Plain reading: without log1p the cost mostly asks\n"
        "\"is this cell a surfactant-secreting cell?\"; with\n"
        "log1p it compares the broader transcriptome.\n\n"
        "Panel 1 = what all three remove: the per-cell library\n"
        f"size spread ({np.percentile(ncount,99)/np.percentile(ncount,1):.0f}x between p1 and p99).\n"
        "None of these equalises PER-PATIENT depth -- separate\n"
        "question, and a real one for a 25-patient axis."
    )
    ax[1, 1].text(0.02, 0.96, txt, va="top", ha="left", fontsize=9.8,
                  family="monospace", linespacing=1.45)

    fig.suptitle("Normalization caliber for the WOT run - epithelial subset "
                 f"(n={a.n_obs} cells x {a.n_vars} genes)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG, dpi=115)
    log(f"图写到 {FIG}")

    # 控制台也打一份榜位，便于我核对图没画错
    for lab, v in [("A raw", mean_raw), ("B CP10K", mean_cpm), ("C log1p", mean_log)]:
        print(f"top20 by {lab:8s}: " + ", ".join(names[np.argsort(-v)[:20]]))
    print(f"\noverlap with C:  A {inter_raw}/20   B {inter_cpm}/20")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
