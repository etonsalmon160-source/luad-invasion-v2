#!/usr/bin/env python
"""
预注册 16_cytotrace_prereg.md 的实现：走原文那条轴（CytoTRACE 发育潜能打分）。

两档（T1 签字）：
  主    = 全部上皮 133,384
  敏感性 = 肺泡类 118,986（AT2 + AT1），各自独立建图

判据 C1-C4 事前写死（§4），阈值 η² ≥ 0.06（T2 签字，Cohen 中效应惯例）。
方向预测：CytoTRACE 分数随病程【单调下降】。方向若相反照实报，不调头重述。

零安装：cellrank 2.0.5 的 CytoTRACEKernel.compute_cytotrace() = Gulati 2020 复刻。
Ms 自己按连通度算（scvelo 在只有 X 的对象上跑不了），见 §3。
"""
import gc, json, os, sys, time
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

import anndata as ad
import scanpy as sc
from cellrank.kernels import CytoTRACEKernel

BASE = "/home/eto/luad_v2/results"
IN_H5AD = f"{BASE}/09_trajectory/wot_full/epiA_wot_input.h5ad"
EMB_DIR = f"{BASE}/04_integration/seurat_trad/epiA/embeddings"
CLUSTERS = f"{BASE}/04_integration/seurat_trad/epiA/clusters.csv.gz"
CLUSTER_ANNOT = f"{BASE}/05_annotation/epiA_cluster_annotation.csv"
OUT = f"{BASE}/09_trajectory/cytotrace"
FIGD = f"{OUT}/figures"

STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]
SIDX = {s: i for i, s in enumerate(STAGES)}
N_NEIGHBORS = 30          # 与 CellRank 教程 moments 默认一致
N_GENES = 200             # cellrank 默认
N_RAND = 50               # K2 阴性对照重复次数
CLUSTER_COL = "harmony_res0.7_seed0"   # 与 epiA_annotation_manifest.json 一致
ARMS = [("all", None), ("alveolar", {"AT2", "AT1"})]
ETA2_MIN = 0.06
P_MAX = 0.01
R_DEPTH_MAX = 0.90
RANDOM_SEED = 20260927

_T0 = time.time()
_HWM = 0.0


def log(msg):
    global _HWM
    try:
        for ln in open("/proc/self/status"):
            if ln.startswith("VmHWM:"):
                _HWM = max(_HWM, int(ln.split()[1]) / 1024 / 1024)
    except Exception:
        pass
    print(f"[{time.time()-_T0:7.1f}s | {_HWM:5.2f}GB] {msg}", flush=True)


def eta2(values, groups):
    """η² = SS_between / SS_total"""
    v = np.asarray(values, float)
    g = np.asarray(groups)
    grand = v.mean()
    ss_tot = ((v - grand) ** 2).sum()
    ss_bet = 0.0
    for lab in np.unique(g):
        m = g == lab
        ss_bet += m.sum() * (v[m].mean() - grand) ** 2
    return float(ss_bet / ss_tot) if ss_tot > 0 else 0.0


def build_ms(A, X, chunk=10000):
    """Ms = D^-1 A X（一阶 kNN 均值），分块防峰值。"""
    rows = []
    n = X.shape[0]
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        Ac = A[s:e]
        M = Ac @ X
        d = np.asarray(Ac.sum(1)).ravel().astype(np.float64)
        d[d == 0] = 1.0
        M = sp.diags(1.0 / d) @ M
        rows.append(sp.csr_matrix(M, dtype=np.float32))
        del M
    return sp.vstack(rows).tocsr()


def stage_medians(score, stage):
    return {s: float(np.median(score[stage == s])) for s in STAGES
            if (stage == s).sum() > 0}


def criteria(score, stage, patient, num_exp_genes, nfeat, label):
    """跑 C1-C4 + K3，返回 dict。"""
    med = stage_medians(score, stage)
    vals = [med[s] for s in STAGES]
    c1_monotone = bool(all(vals[i] > vals[i + 1] for i in range(len(vals) - 1)))
    c1_reverse = bool(all(vals[i] < vals[i + 1] for i in range(len(vals) - 1)))

    groups = [score[stage == s] for s in STAGES if (stage == s).sum() > 0]
    H, p_kw = stats.kruskal(*groups)
    e_stage = eta2(score, stage)
    e_pat = eta2(score, patient)

    # 与两个深度代理的相关（K1 / C4）：取绝对值更大的那个
    r_numgenes = (float(stats.pearsonr(score, num_exp_genes)[0])
                  if np.std(num_exp_genes) > 0 else float("nan"))
    r_nfeat = (float(stats.pearsonr(score, nfeat)[0])
               if np.std(nfeat) > 0 else float("nan"))
    r_worst = max(abs(r_numgenes), abs(r_nfeat))

    # K3 逐患者：ct_score 与期别序号的 Spearman
    per_pat = {}
    for pid in np.unique(patient):
        m = patient == pid
        si = np.array([SIDX[s] for s in stage[m]])
        if len(np.unique(si)) < 2:
            continue
        rho, prho = stats.spearmanr(si, score[m])
        per_pat[str(pid)] = {"n_cells": int(m.sum()), "n_stages": int(len(np.unique(si))),
                             "spearman_rho": float(rho), "p": float(prho)}
    n_neg = sum(1 for v in per_pat.values() if v["spearman_rho"] < 0)
    return {
        "label": label,
        "stage_medians": med,
        "c1_monotone_decreasing": c1_monotone,
        "c1_monotone_increasing": c1_reverse,
        "kruskal_H": float(H), "kruskal_p": float(p_kw),
        "eta2_stage": e_stage, "eta2_patient": e_pat,
        "c2_pass": bool(p_kw < P_MAX and e_stage >= ETA2_MIN),
        "c3_pass": bool(e_pat < 0.5),
        "r_score_vs_num_exp_genes": r_numgenes,
        "r_score_vs_nFeature": r_nfeat,
        "r_depth_worst": r_worst,
        "c4_pass": bool(r_worst < R_DEPTH_MAX),
        "per_patient": per_pat,
        "n_patient_negative_slope": int(n_neg),
        "n_patient_tested": int(len(per_pat)),
    }


def run_arm(arm_name, keep_labels):
    log(f"===== 档 [{arm_name}] 开始 =====")
    adata = ad.read_h5ad(IN_H5AD)
    log(f"读入 h5ad: {adata.shape}")

    # Harmony 附加（全量行序已核验与 obs._index 逐位一致）
    cells = [l.strip() for l in open(f"{EMB_DIR}/cells.txt")]
    assert list(adata.obs_names) == cells, "Harmony cells.txt 与 h5ad 行序不一致！"
    H = np.fromfile(f"{EMB_DIR}/harmony_f32.bin", dtype=np.float32).reshape(len(cells), -1)
    adata.obsm["X_harmony"] = H
    log(f"Harmony 附加完成 {H.shape}")

    # 亚型标签
    cl = pd.read_csv(CLUSTERS, usecols=["cell_barcode", CLUSTER_COL])
    an = pd.read_csv(CLUSTER_ANNOT)
    c2lab = dict(zip(an["cluster"], an["argmax"]))
    adata.obs["epi_subtype"] = (
        cl.set_index("cell_barcode")[CLUSTER_COL].reindex(adata.obs_names).map(c2lab).values)
    n_nan = int(pd.isna(adata.obs["epi_subtype"]).sum())
    log(f"亚型映射完成，未命中 {n_nan} 个细胞")
    log("亚型分布:\n" + adata.obs["epi_subtype"].value_counts().to_string())

    if keep_labels is not None:
        adata = adata[adata.obs["epi_subtype"].isin(keep_labels)].copy()
        log(f"子集到 {keep_labels}: {adata.shape}")

    # 邻接图
    sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=N_NEIGHBORS, n_pcs=None)
    log("邻接图完成")

    # Ms
    A = adata.obsp["connectivities"]
    Ms = build_ms(A, sp.csr_matrix(adata.X))
    adata.layers["Ms"] = Ms
    log(f"Ms 建好 {Ms.shape} nnz={Ms.nnz/1e6:.1f}M")
    del Ms
    gc.collect()

    # CytoTRACE
    k = CytoTRACEKernel(adata)
    k.compute_cytotrace(layer="Ms", n_genes=N_GENES)
    log("CytoTRACE 完成")

    score = np.asarray(adata.obs["ct_score"], float)
    numexp = np.asarray(adata.obs["ct_num_exp_genes"], float)
    stage = np.asarray(adata.obs["stage"].astype(str))
    patient = np.asarray(adata.obs["patient_id"].astype(str))
    nfeat = np.asarray(adata.obs["nFeature"], float)

    res = criteria(score, stage, patient, numexp, nfeat, f"{arm_name}/ct_score")
    log(f"C1 单调降={res['c1_monotone_decreasing']} 单调升={res['c1_monotone_increasing']}")
    log(f"C2 p={res['kruskal_p']:.3e} η²期别={res['eta2_stage']:.4f} pass={res['c2_pass']}")
    log(f"C3 η²患者={res['eta2_patient']:.4f} pass={res['c3_pass']}")
    log(f"C4 r(检出基因数)={res['r_score_vs_num_exp_genes']:.4f} "
        f"r(nFeature)={res['r_score_vs_nFeature']:.4f} pass={res['c4_pass']}")
    log(f"K3 逐患者: {res['n_patient_negative_slope']}/{res['n_patient_tested']} 斜率为负")

    # K2 阴性对照：随机标量场，同一套平滑 + 同一套判据
    rng = np.random.default_rng(RANDOM_SEED + 1)
    d = np.asarray(A.sum(1)).ravel().astype(np.float64); d[d == 0] = 1.0
    null = []
    for _ in range(N_RAND):
        u = rng.normal(size=adata.n_obs)
        us = np.asarray(A @ u).ravel() / d
        groups = [us[stage == s] for s in STAGES if (stage == s).sum() > 0]
        _, pn = stats.kruskal(*groups)
        null.append(eta2(us, stage))
    null = np.array(null)
    res["k2_random_field"] = {
        "n_rep": N_RAND,
        "eta2_stage_mean": float(null.mean()), "eta2_stage_p95": float(np.percentile(null, 95)),
        "observed_eta2_stage": res["eta2_stage"],
        "observed_beats_p95": bool(res["eta2_stage"] > np.percentile(null, 95)),
    }
    log(f"K2 随机场 η²均值={null.mean():.4f} p95={np.percentile(null,95):.4f} "
        f"观测={res['eta2_stage']:.4f}")

    # 逐细胞分数落盘（供后续复用）
    pd.DataFrame({"cell": adata.obs_names, "stage": stage, "patient_id": patient,
                  "epi_subtype": np.asarray(adata.obs["epi_subtype"].astype(str)),
                  "ct_score": score, "ct_pseudotime": np.asarray(adata.obs["ct_pseudotime"], float),
                  "ct_num_exp_genes": numexp, "nFeature": nfeat}) \
        .to_csv(f"{OUT}/cytotrace_per_cell_{arm_name}.tsv.gz", sep="\t",
                index=False, compression="gzip")
    log(f"逐细胞分数已落盘 ({arm_name})")

    # 基因层面：哪些基因进了 top-200
    corr = np.asarray(adata.var["ct_gene_corr"], float)
    isc = np.asarray(adata.var["ct_correlates"], bool)
    pd.DataFrame({"gene": adata.var_names, "gene_corr_with_num_exp_genes": corr,
                  "is_ct_correlate": isc}) \
        .sort_values("gene_corr_with_num_exp_genes", ascending=False) \
        .to_csv(f"{OUT}/cytotrace_gene_corr_{arm_name}.tsv", sep="\t", index=False)
    log(f"基因相关性已落盘：{int(isc.sum())} 个入 top-{N_GENES}")
    res["n_ct_correlates"] = int(isc.sum())

    del adata, A; gc.collect()
    log(f"===== 档 [{arm_name}] 结束，峰值 {_HWM:.2f} GB =====")
    return res


def main():
    os.makedirs(FIGD, exist_ok=True)
    rep = {"prereg": "09_trajectory/16_cytotrace_prereg.md",
           "signed": "T1=全部上皮主+肺泡类敏感性 / T2=η²≥0.06 / T3=monocle3先不装",
           "input": IN_H5AD, "n_neighbors": N_NEIGHBORS, "n_genes": N_GENES,
           "direction_predicted": "CytoTRACE 随病程单调下降",
           "criteria": {"P_MAX": P_MAX, "ETA2_MIN": ETA2_MIN, "R_DEPTH_MAX": R_DEPTH_MAX},
           "arms": {}}
    for name, keep in ARMS:
        rep["arms"][name] = run_arm(name, keep)
        with open(f"{OUT}/cytotrace_summary.json", "w") as f:   # 每档完就落盘
            json.dump(rep, f, indent=2, ensure_ascii=False)

    log("全部完成")
    print(json.dumps(rep, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    main()
