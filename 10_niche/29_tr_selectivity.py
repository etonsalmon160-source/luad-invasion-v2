#!/usr/bin/env python3
# 29_tr_selectivity.py —— 靶点扰动逆向臂：**按 CMap 原义算 selectivity**
#
# 原文（Subramanian 2017，逐字）：
#   「We define the PCL selectivity s of a query q as the fraction of **PCLs** whose connectivity
#     to q is less than a given threshold τ_th.」  s_q = (1/N) Σᵢ [ |τᵢ| < τ_th ]
#   「Given a vector of normalized connectivity scores for the members of a PCL p, relative to
#     query q, in a given cell line, we apply the **maximum quantile** procedure…」
#   「**In the analyses presented here, we used Q_hi = 67, Q_lo = 33**」  ← P 档参数
#
# 🔴 档位：概念与 Q67/Q33 = **P（论文逐字）**；PCL→τ 的参照集口径**原文未给**，
#    本项目操作化为"τ_PCL = |NCS_PCL| 在全部 PCL 中的百分位" ⇒ **S 档**，须在报告标注。
#
# 用法：python3 29_tr_selectivity.py <query> [lib]
import numpy as np, pandas as pd, sys, os
OUT="/home/eto/luad_v2/results/10_niche/target_reversal"; DATA=f"{OUT}/data"
q = sys.argv[1] if len(sys.argv)>1 else "D3"
lib = sys.argv[2] if len(sys.argv)>2 else "L1a"
f=f"{OUT}/{lib}_{q}/ranked.tsv"
d=pd.read_csv(f,sep="\t"); d["pn"]=d.pert_name.astype(str).str.lower()
mo=pd.read_csv(f"{DATA}/clue_moa_list.tsv",sep="\t")
m=d.merge(mo, left_on="pn", right_on="pert_name", how="inner", suffixes=("","_moa"))
print(f"{q}: 签名 {len(d)} 条；能映射到 MOA 的 {len(m)} 条（{100*len(m)/len(d):.1f}%）；MOA 类 {m.moa.nunique()}")

# ① 每个 PCL 的 NCS 向量 → 最大分位聚合（Q67 / Q33）
rows=[]
for moa, g in m.groupby("moa"):
    v = g.ncs.values
    if v.size == 0: continue
    q67, q33 = np.quantile(v, 0.67), np.quantile(v, 0.33)
    ncs_pcl = q67 if abs(q67) >= abs(q33) else q33
    rows.append({"moa": moa, "n_member_sig": v.size, "n_compounds": g.pert_name.nunique(),
                 "q67": q67, "q33": q33, "ncs_pcl": ncs_pcl})
P=pd.DataFrame(rows)
# ② τ_PCL = |NCS_PCL| 在全部 PCL 里的百分位（**S 档：本项目操作化**）
P["tau_pcl"] = np.sign(P.ncs_pcl) * 100.0 / len(P) * np.searchsorted(np.sort(np.abs(P.ncs_pcl)),
                                                                     np.abs(P.ncs_pcl), side="left")
P=P.sort_values("tau_pcl")
P.to_csv(f"{OUT}/{lib}_{q}/selectivity_pcl.tsv", sep="\t", index=False)
s=(P.tau_pcl.abs() < 90).mean()
print(f"  选择性 s = {s:.3f}（|τ_PCL|<90 的 PCL 占比；越接近 1 越选择性）")
print(f"  τ_PCL ≤ −90 的 PCL: {(P.tau_pcl<=-90).sum()} / {len(P)}")
print("\n  逆转方向最强的 PCL（τ_PCL 最负，前 12）：")
print(P.head(12)[["moa","n_compounds","n_member_sig","ncs_pcl","tau_pcl"]].to_string(index=False))
