#!/usr/bin/env python3
"""跑 scMalignantFinder 预训练模型（诊断用），并做 MALIG_PROBE_PREREG.md 里预注册的评估。

环境：`~/venvs/scmalig`（`--system-site-packages` 复用主环境 py3.8 的
scanpy/sklearn/joblib/squidpy，仅 `pip install --no-deps` 装 17KB 的 scMalignantFinder wheel；
**主环境未被修改**）。跑法：
    ~/venvs/scmalig/bin/python 03_cnv/13_run_malig_probe.py

⚠️ 版本偏差：本机 scanpy 1.9.8 / sklearn 1.3.2，官方 pin 为 1.9.3 / 1.2.2。
本脚本只用到 `sc.pp.normalize_total` / `sc.read_h5ad` / `LogisticRegression.predict_proba`（跨版本稳定）。
pyscenic / dask 未装 —— 它们只被 `utils.py` / `spatial.py` 用到，本诊断不碰。

判定口径**写死在预注册里**，本脚本只执行、不改阈值。
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = "/home/eto/luad_v2"
PROBE = f"{ROOT}/results/03_cnv/malig_probe"
IN_H5AD = f"{PROBE}/malig_probe_input.h5ad"
PRETRAIN = f"{ROOT}/tools/scMalignantFinder_pretrain"
CLUSTERS = f"{ROOT}/results/04_integration/seurat_trad/epiA/clusters.csv.gz"
CLUSTER_ANN = f"{ROOT}/results/05_annotation/epiA_cluster_annotation.csv"
FIGDIR = f"{PROBE}/figures"

# —— 预注册常量（MALIG_PROBE_PREREG.md §四）——
POS_STAGE = "IAC"
NEG_STAGE = "Normal"
DESC_STAGES = ["AAH", "AIS", "MIA"]
AUROC_STRONG, AUROC_WEAK = 0.80, 0.60
CLUSTER_COL = "harmony_res0.7_seed0"
STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "IAC"]


def main():
    os.makedirs(FIGDIR, exist_ok=True)

    import anndata as ad
    from scMalignantFinder import classifier

    print(f"[in] {IN_H5AD}")
    adata = ad.read_h5ad(IN_H5AD)
    print(f"[in] {adata.n_obs} 细胞 x {adata.n_vars} 基因")

    model = classifier.scMalignantFinder(
        test_input=adata,
        pretrain_dir=PRETRAIN,
        norm_type=True,
        n_thread=8,
    )
    model.load()
    res = model.predict()

    print(f"[out] 结果对象 {res.n_obs} 细胞")
    print(f"[out] obs 列: {list(res.obs.columns)}")
    for c in res.obs.columns:
        v = res.obs[c]
        if v.dtype.kind in "if":
            print(f"       {c}: dtype={v.dtype} min={v.min():.4g} max={v.max():.4g}")
        else:
            print(f"       {c}: 取值 {list(pd.unique(v))[:6]}")
    print(f"[out] 打分矩阵列: {list(res.to_df().columns)[:10]} ({res.n_vars} 列)")

    # —— 找出打分列（v1.2.0 固定为 malignancy_probability）——
    score_col = None
    for cand in ["malignancy_probability", "malignant_score", "malignancy",
                 "malignant_prob", "probability", "prob", "score"]:
        if cand in res.obs.columns:
            score_col = cand
            break
    assert score_col is not None, f"找不到打分列，obs 列={list(res.obs.columns)}"
    print(f"[key] 打分列 = {score_col}")

    pred_col = None
    for cand in ["scMalignantFinder_prediction", "malignant", "prediction",
                 "predict", "pred", "label", "call"]:
        if cand in res.obs.columns:
            pred_col = cand
            break
    assert pred_col is not None, f"找不到判定列，obs 列={list(res.obs.columns)}"
    print(f"[key] 判定列 = {pred_col}")

    # —— 组装逐细胞表 ——
    per_cell = pd.DataFrame({
        "cell_barcode": res.obs_names,
        "score": res.obs[score_col].values,
    })
    if pred_col is not None:
        per_cell["pred"] = res.obs[pred_col].values

    # 标签：从 seurat 的 clusters.csv.gz 取 stage/patient/sample + 亚型
    cl = pd.read_csv(CLUSTERS, usecols=["cell_barcode", "sample_id", "patient_id",
                                        "stage", CLUSTER_COL])
    ann = pd.read_csv(CLUSTER_ANN, usecols=["cluster", "argmax"])
    cl["epi_subtype"] = cl[CLUSTER_COL].map(dict(zip(ann["cluster"], ann["argmax"])))
    assert cl["epi_subtype"].notna().all(), "有细胞的簇号在注释表里找不到，停下"

    per_cell = per_cell.merge(
        cl[["cell_barcode", "sample_id", "patient_id", "stage", "epi_subtype",
            CLUSTER_COL]].rename(columns={CLUSTER_COL: "epi_cluster"}),
        on="cell_barcode", how="left", validate="one_to_one")
    n_miss = int(per_cell["stage"].isna().sum())
    assert n_miss == 0, f"{n_miss} 个细胞匹配不到分期标签，停下"

    # 交叉核对：h5ad obs 自带 stage（来自 cohort_registry.resolve_stage）必须与 clusters 表一致
    h5_stage = adata.obs["stage"].astype(str)
    chk = per_cell.set_index("cell_barcode")["stage"].reindex(h5_stage.index)
    n_disagree = int((chk.values != h5_stage.values).sum())
    assert n_disagree == 0, f"{n_disagree} 个细胞分期在 h5ad 与 clusters 表间不一致，停下"
    print(f"[chk] 分期与 h5ad obs 逐细胞一致（{len(h5_stage)} 细胞）")

    print(f"[join] 合并后 {len(per_cell)} 行；分期分布 "
          f"{per_cell['stage'].value_counts().to_dict()}")

    per_cell.to_csv(f"{PROBE}/malig_probe_per_cell.csv.gz", index=False,
                    compression="gzip")

    # —— 预注册评估：逐患者 AUROC（Normal vs IAC）——
    rows = []
    for pid, g in per_cell.groupby("patient_id"):
        pos = g.loc[g["stage"] == POS_STAGE, "score"]
        neg = g.loc[g["stage"] == NEG_STAGE, "score"]
        if len(pos) == 0 or len(neg) == 0:
            rows.append({"patient_id": pid, "n_pos": len(pos), "n_neg": len(neg),
                         "auroc": np.nan, "note": "该患者缺一类，跳过"})
            continue
        y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
        s = np.r_[pos.values, neg.values]
        if len(np.unique(y)) < 2:
            rows.append({"patient_id": pid, "n_pos": len(pos), "n_neg": len(neg),
                         "auroc": np.nan, "note": "单类"})
            continue
        rows.append({"patient_id": pid, "n_pos": len(pos), "n_neg": len(neg),
                     "auroc": float(roc_auc_score(y, s)), "note": ""})
    per_pat = pd.DataFrame(rows).sort_values("patient_id")
    per_pat.to_csv(f"{PROBE}/malig_probe_per_patient.csv", index=False)

    valid = per_pat["auroc"].dropna()
    med = float(valid.median()) if len(valid) else float("nan")
    if med >= AUROC_STRONG:
        verdict = "表达轴有信号，可作候选第二轴（仍须另行签字才能采用）"
    elif med >= AUROC_WEAK:
        verdict = "弱信号，只能作佐证，不足以定恶性"
    else:
        verdict = "表达分不开，不采用"

    print(f"\n[VERDICT] 逐患者 AUROC 中位数 = {med:.4f} （有效患者 {len(valid)}/{len(per_pat)}）")
    print(f"[VERDICT] 预注册判读：{verdict}")
    print(per_pat.to_string(index=False))

    # —— 各分期打分分布（描述性）——
    stage_stat = per_cell.groupby("stage")["score"].describe()
    print(f"\n[描述] 各分期打分：\n{stage_stat.to_string()}")

    summary = {
        "script": "03_cnv/13_run_malig_probe.py",
        "tool": "scMalignantFinder (pretrained)",
        "pretrain_dir": PRETRAIN,
        "input_h5ad": IN_H5AD,
        "n_cells": int(len(per_cell)),
        "score_col": score_col,
        "pred_col": pred_col,
        "pos_stage": POS_STAGE,
        "neg_stage": NEG_STAGE,
        "desc_stages": DESC_STAGES,
        "auroc_strong": AUROC_STRONG,
        "auroc_weak": AUROC_WEAK,
        "per_patient_auroc_median": med,
        "n_patients_with_auroc": int(len(valid)),
        "n_patients_total": int(len(per_pat)),
        "verdict": verdict,
        "stage_score_describe": json.loads(stage_stat.to_json(orient="index")),
        "caveats": [
            "样本级批次混淆无法排除：Normal 与 IAC 样本在解离/环境RNA/批次上不同",
            "IAC 样本内的瘤旁正常上皮被误标为正类，拉低表观性能",
            "无逐细胞真值 ⇒ 只回答'能不能分开两个样本群'，不回答'每个核是否恶性'",
            "该工具及同类前辈均未在 snRNA 上验证过；无 AAH/AIS 前例",
            "上皮亚型标签来自 r*=0.7（放宽产物）⇒ 亚型相关结论为探索性",
            "版本偏差：本机 scanpy 1.9.8 / sklearn 1.3.2，官方 pin 1.9.3 / 1.2.2",
        ],
        "env": {
            "venv": "/home/eto/venvs/scmalig",
            "python": "3.8.10",
            "scMalignantFinder": "1.2.0",
            "install": "pip install --no-deps（主环境未改）",
            "scanpy": "1.9.8",
            "sklearn": "1.3.2",
            "joblib": "1.4.2",
        },
    }
    with open(f"{PROBE}/malig_probe_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"\n[out] {PROBE}/malig_probe_summary.json")

    # —— 图（图内标签全英文）——
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax = axes[0]
    data, labels = [], []
    for st in STAGE_ORDER:
        v = per_cell.loc[per_cell["stage"] == st, "score"].values
        if len(v):
            data.append(v)
            labels.append(f"{st}\n(n={len(v)})")
    bp = ax.boxplot(data, labels=labels, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=2))
    for i, box in enumerate(bp["boxes"]):
        box.set_facecolor(plt.cm.viridis(i / max(1, len(data) - 1)))
        box.set_alpha(0.75)
    ax.set_ylabel("malignancy score")
    ax.set_xlabel("stage (sample-level)")
    ax.set_title("scMalignantFinder score by stage\n(epithelial nuclei, GSE308103)")
    ax.grid(axis="y", ls=":", alpha=0.4)

    ax = axes[1]
    ok = per_pat.dropna(subset=["auroc"])
    colors = ["#c0392b" if v < AUROC_WEAK else ("#e67e22" if v < AUROC_STRONG else "#27ae60")
              for v in ok["auroc"]]
    ax.bar(ok["patient_id"], ok["auroc"], color=colors)
    ax.axhline(AUROC_STRONG, color="#27ae60", ls="--", lw=1.5,
               label=f"strong >= {AUROC_STRONG}")
    ax.axhline(AUROC_WEAK, color="#e67e22", ls="--", lw=1.5,
               label=f"weak >= {AUROC_WEAK}")
    ax.axhline(0.5, color="grey", ls=":", lw=1.2, label="chance = 0.5")
    ax.axhline(med, color="black", lw=2, label=f"median = {med:.3f}")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("AUROC (IAC vs Normal)")
    ax.set_xlabel("patient")
    ax.set_title("Per-patient AUROC on own normal-vs-tumour epithelium")
    ax.legend(fontsize=8, loc="lower left")
    ax.tick_params(axis="x", rotation=90)
    ax.grid(axis="y", ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/malig_probe_score_by_stage.png", dpi=150)
    print(f"[out] {FIGDIR}/malig_probe_score_by_stage.png")


if __name__ == "__main__":
    sys.exit(main())
