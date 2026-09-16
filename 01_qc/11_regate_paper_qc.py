#!/usr/bin/env python3
# ============================================================================
# 11_regate_paper_qc.py —— 按**源论文 QC 口径**重建分析掩膜（GP0 重做）
#
# 背景与授权：
#   2026-09-16 用户决策「改用论文 QC 重建」+「对齐论文全套」。此前 M1 用的是
#   **逐样本自适应 MAD**（本项目自定口径），比权威口径多留 36% 的核。本脚本把
#   分析掩膜换成源论文的**固定全局阈值**，作为新的唯一分析集。
#
# 源论文（权威出处，不得改写）：
#   Peng F, Sinjab A, Dai Y, … Wang L, Kadara H. "Multimodal spatial-omics
#   reveal co-evolution of alveolar progenitors and proinflammatory niches in
#   progression of lung precursor lesions." Cancer Cell 2026;44(2):321-339.e13.
#   DOI 10.1016/j.ccell.2025.10.004 · PMC12980502 · GEO Series_pubmed_id=41202811
#   其 Methods 载明：剔除 nFeature<500 或 nCount<1000 或 pct_mt>20% 的核；
#   并保留在 ≥3 个细胞中检出的基因。论文最终保留 401,635 个核。
#
# ★ 双体口径的**显式偏差**（铁律 R2 优先于"对齐论文"）：
#   论文用 **Scrublet** 去双体；本项目铁律 R2 钉死 **scDblFinder**。
#   本脚本因此**保留 M1 的 scDblFinder 单细胞判定**，在其之上叠加论文计数门。
#   这是一条登记偏差，不是疏漏。M1 期间 Scrublet 已作为**敏感性臂**跑过
#   （01_qc/05_doublets_scrublet.py → results/01_qc/gse308103_sensitivity_doublet_rate.csv）。
#
# ⚠️ 本脚本**只做掩膜计算**，不改任何矩阵。矩阵重建见 02_expression/04_*。
# ⚠️ 同时在新细胞集上**重做** 00c 的 AAH 脆弱性检验 —— 00c 的结论是在旧口径
#    （648,945 核）上得出的，换口径必须重验，不得直接沿用。
#
# 输入：results/01_qc/gse308103_per_cell_qc.csv.gz        （M1 逐细胞 QC + 双体）
#       results/02_expression/gse308103_analysis_mask.csv.gz（M1 掩膜，供 patient/stage）
# 输出：results/01_qc/gse308103_analysis_mask_paperqc.csv.gz
#       results/01_qc/regate_paper_qc_report.json
#       results/01_qc/regate_paper_qc_cascade.csv
#       results/01_qc/regate_paper_qc_by_stage.csv
#       results/01_qc/regate_paper_qc_by_sample.csv
#       results/01_qc/stage_fragility_report_paperqc.json
# ============================================================================
import json, time, hashlib, sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT  = Path("/home/eto/luad_v2")
QC    = ROOT / "results/01_qc/gse308103_per_cell_qc.csv.gz"
MASK  = ROOT / "results/02_expression/gse308103_analysis_mask.csv.gz"
OUT   = ROOT / "results/01_qc"

# ---- 论文阈值（显式常量，禁止散落在代码里）--------------------------------
MIN_FEATURE = 500
MIN_COUNT   = 1000
MAX_PCT_MT  = 20.0          # 论文：「>20% 被剔除」⇒ 保留 pct_mt <= 20
PAPER_N_RETAINED = 401635   # 论文公布值，仅用于对照，**不作为断言目标**

ORDER = ["Normal", "AAH", "AIS", "MIA", "IAC"]
T0 = time.time()
def log(m): print(f"[{time.time()-T0:6.1f}s] {m}", flush=True)

def sha256(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while (b := f.read(buf)):
            h.update(b)
    return h.hexdigest()

# ---- 1. 读入 ---------------------------------------------------------------
log(f"读入 {QC.name}")
qc = pd.read_csv(QC)
log(f"  {len(qc)} 细胞 × {len(qc.columns)} 列：{list(qc.columns)}")
m1 = pd.read_csv(MASK)
log(f"读入 M1 掩膜 {len(m1)} 细胞（patient/stage 来源）")

if len(qc) != 798100:
    sys.exit(f"[FAIL] M1 QC 表行数 {len(qc)} != 798100（预期全集），拒绝继续")
if m1.cell_barcode.duplicated().any() or qc.cell_barcode.duplicated().any():
    sys.exit("[FAIL] cell_barcode 有重复，连接键不唯一")

# 把 patient_id / stage 从 M1 掩膜挂到 QC 表（左连接，M1 掩膜是 QC 表的子集）
df = qc.merge(m1, on="cell_barcode", how="left", suffixes=("", "_m1"))
grafted = df.stage.notna().sum()
log(f"  patient/stage 挂载：{grafted} 细胞有标注（应 == {len(m1)}）")
if grafted != len(m1):
    sys.exit(f"[FAIL] 挂载数 {grafted} != M1 掩膜 {len(m1)}")
# 交叉核对 sample_id 两表一致
bad = df.loc[df.stage.notna() & (df.sample_id != df.sample_id_m1)]
if len(bad):
    sys.exit(f"[FAIL] {len(bad)} 个细胞的 sample_id 在两表间不一致")

# ---- 2. 掩膜级联（逐步上报每一步砍了多少）--------------------------------
c_paper = (df.nFeature >= MIN_FEATURE) & (df.nCount >= MIN_COUNT) & (df.pct_mt <= MAX_PCT_MT)
c_singlet = df.doublet_class == "singlet"
c_m1mask  = df.stage.notna()

cascade = [
    ("全集（M1 输入）",                          len(df),                                   None),
    ("M1 qc_pass",                              int(df.qc_pass.sum()),                     "逐样本自适应 MAD"),
    ("M1 qc_pass & singlet = M1 掩膜",           int(c_m1mask.sum()),                       "scDblFinder 单细胞"),
    ("论文计数门（作用于全集）",                   int(c_paper.sum()),                        "nFeature>=500 & nCount>=1000 & pct_mt<=20"),
    ("M1 掩膜 ∩ 论文计数门 = **新分析集**",        int((c_m1mask & c_paper).sum()),            "双体用 scDblFinder（R2）"),
]
cas = pd.DataFrame(cascade, columns=["step", "n_cells", "口径"])
cas["pct_of_all"] = (100 * cas.n_cells / len(df)).round(3)
log("\n=== 掩膜级联 ===")
log("\n" + cas.to_string(index=False))

n_new = int((c_m1mask & c_paper).sum())

# ---- 3. 论文计数门在 M1 单细胞集内砍掉的分解 -------------------------------
m1s = c_m1mask
rem_feat  = int((m1s & (df.nFeature <  MIN_FEATURE)).sum())
rem_cnt   = int((m1s & (df.nFeature >= MIN_FEATURE) & (df.nCount < MIN_COUNT)).sum())
rem_mt    = int((m1s & (df.nFeature >= MIN_FEATURE) & (df.nCount >= MIN_COUNT)
                     & (df.pct_mt > MAX_PCT_MT)).sum())
rem_any   = int((m1s & ~c_paper).sum())
log(f"\nM1 单细胞集内，论文门剔除 {rem_any} 核：")
log(f"  nFeature<{MIN_FEATURE}            : {rem_feat}")
log(f"  仅 nCount<{MIN_COUNT}（feature 已过）: {rem_cnt}")
log(f"  仅 pct_mt>{MAX_PCT_MT}%           : {rem_mt}")
if rem_feat + rem_cnt + rem_mt != rem_any:
    sys.exit(f"[FAIL] 剔除分解不闭合：{rem_feat}+{rem_cnt}+{rem_mt} != {rem_any}")

# 边界敏感性：恰好 20% 的核有多少（若阈值是 <20 而非 <=20，会多砍这些）
n_eq20 = int((m1s & (df.pct_mt == MAX_PCT_MT)).sum())
log(f"  边界：pct_mt 恰为 20.000 的 M1 单细胞 = {n_eq20}（<=20 与 <20 的差异上界）")

# ---- 4. 新分析集的逐分期 / 逐样本构成 --------------------------------------
new = df.loc[c_m1mask & c_paper,
             ["cell_barcode", "sample_id", "patient_id", "stage_token", "stage"]].copy()
new = new.sort_values("cell_barcode").reset_index(drop=True)

by_stage = (new.groupby("stage", observed=True)
              .agg(n_cells=("cell_barcode", "size"),
                   n_samples=("sample_id", "nunique"),
                   n_patients=("patient_id", "nunique"))
              .reindex(ORDER))
# 对照 M1 掩膜的逐分期
old_stage = (df.loc[c_m1mask].groupby("stage", observed=True)
               .agg(n_cells_m1=("cell_barcode", "size")).reindex(ORDER))
by_stage = by_stage.join(old_stage)
by_stage["retention_pct"] = (100 * by_stage.n_cells / by_stage.n_cells_m1).round(2)
by_stage["pct_of_new"]    = (100 * by_stage.n_cells / n_new).round(3)
by_stage.index.name = "stage"
log("\n=== 新分析集逐分期（含相对 M1 的留存率）===")
log("\n" + by_stage.to_string())

# 关键可比性检查：AAH 的留存率是否显著低于 Normal
aah_ret = float(by_stage.loc["AAH", "retention_pct"])
nor_ret = float(by_stage.loc["Normal", "retention_pct"])
log(f"\nAAH 留存 {aah_ret:.2f}%  vs  Normal 留存 {nor_ret:.2f}%  "
    f"⇒ 差 {aah_ret-nor_ret:+.2f} pp")

by_sample = (new.groupby(["sample_id", "stage"], observed=True)
                .size().rename("n_cells").reset_index().sort_values("sample_id"))

# ---- 5. 新细胞集上的 AAH 脆弱性重验（替代 00c 的旧口径结论）----------------
qc_new = df.loc[c_m1mask & c_paper]
rows = []
for s in ORDER:
    v = qc_new.loc[qc_new.stage == s]
    rows.append({
        "stage": s, "n_cells": len(v),
        "n_samples": v.sample_id.nunique(), "n_patients": v.patient_id.nunique(),
        "nFeature_p10": round(float(v.nFeature.quantile(.10)), 1),
        "nFeature_med": round(float(v.nFeature.median()), 1),
        "nCount_p10":   round(float(v.nCount.quantile(.10)), 1),
        "nCount_med":   round(float(v.nCount.median()), 1),
        "pct_mt_med":   round(float(v.pct_mt.median()), 4),
    })
frag = pd.DataFrame(rows).set_index("stage")
frag["nFeature_med_vs_Normal"] = (frag.nFeature_med / frag.loc["Normal", "nFeature_med"]).round(3)

thr_rows = []
for thr in [200, 500, 1000, 1500]:
    surv = {s: 100.0 * float((qc_new.loc[qc_new.stage == s, "nFeature"] >= thr).mean()) for s in ORDER}
    thr_rows.append({"threshold": f"nFeature>={thr}", **{s: round(surv[s], 2) for s in ORDER},
                     "AAH_minus_Normal_pp": round(surv["AAH"] - surv["Normal"], 2)})
thr = pd.DataFrame(thr_rows).set_index("threshold")
log("\n=== [新口径] 逐分期 QC ===")
log("\n" + frag.to_string())
log("\n=== [新口径] 计数门槛逐分期存活率（%）===")
log("\n" + thr.to_string())

worst_gap = float(thr["AAH_minus_Normal_pp"].min())
count_level_fragility = worst_gap < -2.0
med_ratio = float(frag.loc["AAH", "nFeature_med_vs_Normal"])
log(f"\n[新口径] AAH 中位 nFeature/Normal = {med_ratio:.3f}×；"
    f"最差门槛存活差 = {worst_gap:+.2f} pp ⇒ 计数层面脆弱？ "
    f"{'是' if count_level_fragility else '否'}")

# ---- 6. 落盘 ---------------------------------------------------------------
mask_p = OUT / "gse308103_analysis_mask_paperqc.csv.gz"
new.to_csv(mask_p, index=False, compression="gzip")
cas.to_csv(OUT / "regate_paper_qc_cascade.csv", index=False)
by_stage.reset_index().to_csv(OUT / "regate_paper_qc_by_stage.csv", index=False)
by_sample.to_csv(OUT / "regate_paper_qc_by_sample.csv", index=False)
log(f"\n写出 {mask_p}（{len(new)} 行）")

rep = {
    "script": "01_qc/11_regate_paper_qc.py",
    "purpose": "按源论文固定全局 QC 口径重建分析掩膜（替代 M1 自适应 MAD）",
    "decided_by": "用户 2026-09-16 AskUserQuestion：「改用论文 QC 重建」",
    "source_paper": {
        "citation": "Peng F, Sinjab A, Dai Y, … Wang L, Kadara H. Cancer Cell 2026;44(2):321-339.e13",
        "doi": "10.1016/j.ccell.2025.10.004", "pmcid": "PMC12980502",
        "qc_rule": "exclude nFeature<500 | nCount<1000 | pct_mt>20%",
        "paper_n_retained": PAPER_N_RETAINED,
    },
    "thresholds": {"min_nFeature": MIN_FEATURE, "min_nCount": MIN_COUNT,
                   "max_pct_mt": MAX_PCT_MT, "gene_min_cells": 3,
                   "gene_filter_stage": "在矩阵重建时施加（02_expression），本脚本不含"},
    "deviation_from_paper": {
        "item": "doublet caller",
        "paper": "Scrublet", "ours": "scDblFinder",
        "reason": "铁律 R2 钉死 scDblFinder，优先于「对齐论文」",
        "mitigation": "M1 已跑 Scrublet 敏感性臂：results/01_qc/gse308103_sensitivity_doublet_rate.csv",
    },
    "cascade": cas.to_dict(orient="records"),
    "removed_within_m1_singlets": {"total": rem_any, "by_nFeature": rem_feat,
                                   "by_nCount_only": rem_cnt, "by_pct_mt_only": rem_mt,
                                   "n_pct_mt_exactly_20": n_eq20},
    "new_analysis_set": {
        "n_cells": n_new,
        "vs_paper_retained": {"paper": PAPER_N_RETAINED, "delta": n_new - PAPER_N_RETAINED,
                              "delta_pct": round(100 * (n_new - PAPER_N_RETAINED) / PAPER_N_RETAINED, 2)},
        "by_stage": by_stage.reset_index().to_dict(orient="records"),
    },
    "aah_fragility_recheck": {
        "note": "在新分析集上重做 00c 的检验；00c 的旧口径结论不再直接沿用",
        "by_stage": frag.reset_index().to_dict(orient="records"),
        "threshold_survival_pct": thr.reset_index().to_dict(orient="records"),
        "verdict": {
            "aah_median_nFeature_ratio_vs_Normal": med_ratio,
            "worst_survival_gap_pp": worst_gap,
            "count_level_fragility": count_level_fragility,
            "aah_retention_pct": aah_ret, "normal_retention_pct": nor_ret,
            "aah_retention_minus_normal_pp": round(aah_ret - nor_ret, 2),
        },
    },
    "inputs": {"per_cell_qc": {"path": str(QC), "sha256": sha256(QC)},
               "m1_mask":     {"path": str(MASK), "sha256": sha256(MASK)}},
    "outputs": {"mask_paperqc": {"path": str(mask_p), "sha256": sha256(mask_p),
                                 "n_cells": n_new}},
    "versions": {"pandas": pd.__version__, "numpy": np.__version__},
    "wall_sec": round(time.time() - T0, 1),
}
rep_p = OUT / "regate_paper_qc_report.json"
rep_p.write_text(json.dumps(rep, indent=2, ensure_ascii=False))
log(f"写出 {rep_p}")

frag_rep = {
    "script": "01_qc/11_regate_paper_qc.py（内嵌，重做 00c）",
    "supersedes": "results/04_integration/stage_fragility_report.json（旧 648,945 口径）",
    "cell_set": f"论文 QC ∩ scDblFinder 单细胞，n={n_new}",
    "by_stage": frag.reset_index().to_dict(orient="records"),
    "threshold_survival_pct": thr.reset_index().to_dict(orient="records"),
    "verdict": rep["aah_fragility_recheck"]["verdict"],
}
(OUT / "stage_fragility_report_paperqc.json").write_text(
    json.dumps(frag_rep, indent=2, ensure_ascii=False))
log("写出 stage_fragility_report_paperqc.json")
log(f"完成：新分析集 {n_new} 核（论文 {PAPER_N_RETAINED}，差 "
    f"{100*(n_new-PAPER_N_RETAINED)/PAPER_N_RETAINED:+.2f}%）")
