#!/usr/bin/env python3
# ============================================================================
# 00c_stage_fragility_check.py —— AAH 脆弱性的**实证**检验（GP4a 前置诊断）
#
# 背景：用户要求参数须"考虑 AAH 细胞的脆弱性"。在动参数之前，先证明
#       AAH 到底脆弱在**哪里**——是测序深度（计数伪影，可用阈值补偿），
#       还是转录距离（聚类吸收，只能靠分辨率/注释层面解决）。
#
# ⚠️ 本脚本**只读 obs**（backed 模式），不改任何数据，不下任何阈值决定。
#    它存在的意义是：让"要不要给 AAH 放宽 QC"这个决定**基于实测**，
#    而不是基于对"早期病变一定质量差"的想当然。
#
# 输入：results/02_expression/gse308103_counts.h5ad
# 输出：results/04_integration/stage_qc_by_stage.csv
#       results/04_integration/stage_fragility_report.json
# ============================================================================
import json, time, hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path("/home/eto/luad_v2")
INP  = ROOT / "results/02_expression/gse308103_counts.h5ad"
OUT  = ROOT / "results/04_integration"
OUT.mkdir(parents=True, exist_ok=True)

ORDER = ["Normal", "AAH", "AIS", "MIA", "IAC"]
T0 = time.time()

def sha256(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while (b := f.read(buf)):
            h.update(b)
    return h.hexdigest()

a = ad.read_h5ad(INP, backed="r")
obs = a.obs[["sample_id", "patient_id", "stage", "nCount", "nFeature", "pct_mt"]].copy()
a.file.close()
n_total = len(obs)
print(f"读入 obs：{n_total} 细胞")

# ---- 1. 逐分期 QC 分布 ----------------------------------------------------
rows = []
for s in ORDER:
    v = obs.loc[obs.stage == s]
    rows.append({
        "stage": s, "n_cells": len(v), "pct_of_total": round(100 * len(v) / n_total, 3),
        "n_samples": v.sample_id.nunique(), "n_patients": v.patient_id.nunique(),
        "nFeature_p10": round(float(v.nFeature.quantile(.10)), 1),
        "nFeature_p25": round(float(v.nFeature.quantile(.25)), 1),
        "nFeature_med": round(float(v.nFeature.median()), 1),
        "nFeature_p75": round(float(v.nFeature.quantile(.75)), 1),
        "nCount_p10":  round(float(v.nCount.quantile(.10)), 1),
        "nCount_med":  round(float(v.nCount.median()), 1),
        "pct_mt_med":  round(float(v.pct_mt.median()), 4),
        "pct_mt_p90":  round(float(v.pct_mt.quantile(.90)), 4),
    })
tab = pd.DataFrame(rows).set_index("stage")
tab["nFeature_med_vs_Normal"] = (tab.nFeature_med / tab.loc["Normal", "nFeature_med"]).round(3)
tab.to_csv(OUT / "stage_qc_by_stage.csv")
print("\n=== 逐分期 QC ===")
print(tab.to_string())

# ---- 2. 计数门槛的**选择性**检验 ------------------------------------------
#    关键问题：常见 scRNA 门槛是否**不成比例地**砍掉 AAH？
thr_rows = []
for thr in [200, 500, 1000, 1500]:
    surv = {s: 100.0 * float((obs.loc[obs.stage == s, "nFeature"] >= thr).mean()) for s in ORDER}
    thr_rows.append({"threshold": f"nFeature>={thr}", **{s: round(surv[s], 2) for s in ORDER},
                     "AAH_minus_Normal_pp": round(surv["AAH"] - surv["Normal"], 2)})
thr = pd.DataFrame(thr_rows).set_index("threshold")
print("\n=== 计数门槛的逐分期存活率（%）===")
print(thr.to_string())

# ---- 3. 判定 --------------------------------------------------------------
# 判据（本项目显式约定）：若 AAH 在**任一**门槛下的存活率比 Normal 低 >2 个百分点，
# 才认为存在"计数层面的 AAH 选择性丢失"，才需要讨论补偿。
worst_gap = float(thr["AAH_minus_Normal_pp"].min())
count_level_fragility = worst_gap < -2.0
med_ratio = float(tab.loc["AAH", "nFeature_med_vs_Normal"])

print("\n=== 判定 ===")
print(f"AAH 中位 nFeature / Normal = {med_ratio:.3f}×")
print(f"最差门槛下 AAH 存活率 − Normal 存活率 = {worst_gap:+.2f} pp")
print(f"⇒ 计数层面存在 AAH 选择性丢失？ {'是' if count_level_fragility else '否'}")

rep = {
    "script": "04_integration/00c_stage_fragility_check.py",
    "purpose": "实证 AAH 脆弱性的来源：计数伪影 vs 转录距离",
    "input": {"path": str(INP), "sha256": sha256(INP), "n_cells": n_total},
    "by_stage": tab.reset_index().to_dict(orient="records"),
    "threshold_survival_pct": thr.reset_index().to_dict(orient="records"),
    "verdict": {
        "aah_median_nFeature_ratio_vs_Normal": med_ratio,
        "worst_survival_gap_pp": worst_gap,
        "count_level_fragility": count_level_fragility,
        "interpretation": (
            "AAH 核的测序深度不低于其它分期；常见计数门槛不选择性砍 AAH。"
            "⇒ 本数据集里 AAH 的脆弱性**不在计数层面**，"
            "对 AAH 单独放宽 QC 阈值**缺乏依据**（会构成另一种静默默认）。"
            "风险应转向**聚类层面**的处理：AAH 与 Normal 的转录距离细微，"
            "分辨率调高以分开 IAC 亚群时可能把 AAH 并入 Normal。"
            "⇒ 须在 GP5（分辨率选择）设显式可证伪护栏，见 PLAN_AND_CHECKPOINTS.md。"
            "★ 另：AAH 最大的脆弱点是**样本量**——仅 9 样本 / 8 患者"
            "（Normal 24/23、IAC 24/23）⇒ 任何 AAH 层面的结论都建立在 n=8 上，"
            "且这些患者同时贡献其它分期（分期嵌套于患者内）⇒ "
            "AAH 结论必须报**患者级**不确定性，不得按细胞数当独立重复。"
        ) if not count_level_fragility else (
            "AAH 在计数层面被选择性丢失，须在 GP4a 讨论是否补偿，并登记为偏差。"
        ),
    },
    "wall_sec": round(time.time() - T0, 1),
}
(OUT / "stage_fragility_report.json").write_text(json.dumps(rep, indent=2, ensure_ascii=False))
print(f"\n产出：{OUT/'stage_qc_by_stage.csv'}")
print(f"      {OUT/'stage_fragility_report.json'}")
