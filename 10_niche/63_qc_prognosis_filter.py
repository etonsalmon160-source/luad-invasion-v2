#!/usr/bin/env python3
# 63_qc_prognosis_filter.py —— 逆转签名的**预后质控**（回答"逆回去是好是坏"）
#
# 规则（用 TCGA-LUAD 483 例的**独立生存证据**，不是我的判断）：
#   疾病里升高（lfc>0）⇒ 我们的意图是"敲低它" ⇒ 只有当 **HR_adj > 1（高表达⇒预后差）** 才保留
#   疾病里降低（lfc<0）⇒ 我们的意图是"补回来" ⇒ 只有当 **HR_adj < 1（高表达⇒预后好）** 才保留
#   其余（HR 方向与意图相反 = 保护性程序；或无 HR）⇒ **置零剔除**
import numpy as np, pandas as pd
OUT = "/home/eto/luad_v2/results/10_niche/tr_singlecell"
HR = pd.read_csv(f"{OUT}/tcga_pergene_HR.tsv", sep="\t")
hr = dict(zip(HR.gene, HR.HR)); hq = dict(zip(HR.gene, HR.p))

for nm in ["global_paired", "global_all"]:
    S = pd.read_csv(f"{OUT}/{nm}_lfc.tsv", sep="\t")
    S["HR"] = S.gene.map(hr); S["p_HR"] = S.gene.map(hq)
    n0 = len(S)
    # 意图
    S["intent"] = np.where(S.lfc > 0, "敲低", "补回")
    # 判据
    ok = np.zeros(len(S), dtype=bool)
    m_up = S.lfc > 0; m_dn = S.lfc < 0
    ok[m_up] = (S.HR[m_up] > 1).fillna(False).values
    ok[m_dn] = (S.HR[m_dn] < 1).fillna(False).values
    S["keep"] = ok
    S["lfc_qc"] = np.where(ok, S.lfc, 0.0)
    S.to_csv(f"{OUT}/{nm}_qc.tsv", sep="\t", index=False)
    print(f"\n=== {nm} ===")
    print(f"  轴里非零基因 {int((S.lfc!=0).sum()):,}；有 HR 的 {int(S.HR.notna().sum()):,}")
    print(f"  ✅ 通过质控 {int(ok.sum()):,} ；❌ 剔除 {int(((S.lfc!=0)&~ok).sum()):,}")
    # 被剔除里最可惜的（按 |lfc| 排）
    bad = S[(S.lfc != 0) & (~ok)].copy()
    bad["abs"] = bad.lfc.abs()
    print(f"  【被剔除里 |lfc| 最大的 12 个（疑似保护性/无证据）】")
    print(bad.nlargest(12, "abs")[["gene", "intent", "lfc", "HR", "p_HR"]].round(4).to_string(index=False))
    print(f"  【保留里 |lfc| 最大的 10 个】")
    good = S[ok].copy(); good["abs"] = good.lfc.abs()
    print(good.nlargest(10, "abs")[["gene", "intent", "lfc", "HR", "p_HR"]].round(4).to_string(index=False))

# 关键抽查：免疫/干扰素那批
print("\n=== 抽查：之前担心的那批 ===")
S = pd.read_csv(f"{OUT}/global_paired_qc.tsv", sep="\t").set_index("gene")
for g in ["STAT1", "IFNGR2", "ADAR", "XBP1", "FOXP3", "IGKC", "IGHG1", "NKX2-1", "SFTPC", "TCF21", "CLDN18"]:
    if g in S.index:
        r = S.loc[g]
        print(f"  {g:8s} lfc={r.lfc:+.3f}  意图={r.intent}  HR={r.HR if pd.notna(r.HR) else float('nan'):.3f}  ⇒ {'保留' if r.keep else '剔除'}")
