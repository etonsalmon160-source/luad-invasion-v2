#!/usr/bin/env python3
# 78_plip_density_corr.py —— 补落的脚本：PLIP 分与"组织占比"到底相关多少
#
# 背景：草稿 §13 引用 `corr(neoplastic, 组织占比) = 0.825` 作为"深度–密度耦联"的独立旁证。
#       该数字的生成脚本此前未落盘。2026-10-06 审计时定位到它的真实出处并核出错配。
#
# 🔴 结论：0.825 的来源是 `corr_lesion_vs_empty_within_slide / GSM9226174_P4_Normal = −0.8250`
#    —— 那是 **P4_Normal 单张切片上「lesion 分 vs 空片分」的相关**，**不是**组织占比的相关，
#    而且只来自一张切片。真正的 `corr_lesion_vs_tissuefrac` = **0.754**（P4 四张切片）。
import json, os
import pandas as pd

R = "/home/eto/luad_v2/results/07_he_pathology"
print("=" * 72)
print("PLIP 分与组织占比的相关 —— 逐来源核对")
print("=" * 72)

# ① P4 四切片对照（prompt_compare）
p = f"{R}/prompt_compare/prompt_compare_report.json"
if os.path.exists(p):
    d = json.load(open(p))
    print(f"\n[prompt_compare]  三个提示词集并列（{'═'*40}）")
    print(f"  {'set':44s} {'corr_组织占比':>12s} {'corr_空片':>10s} {'P4_Normal切片内':>15s}")
    for setname, s in d["sets"].items():
        tf = s.get("corr_lesion_vs_tissuefrac"); em = s.get("corr_lesion_vs_empty")
        w = s.get("corr_lesion_vs_empty_within_slide", {}).get("GSM9226174_P4_Normal")
        fm = lambda x: f"{x:+.4f}" if isinstance(x, (int, float)) else "—"
        star = "  <-- 草稿引的 0.825 在此" if isinstance(w, float) and abs(abs(w) - 0.825) < 1e-3 else ""
        print(f"  {setname:44s} {fm(tf):>12s} {fm(em):>10s} {fm(w):>15s}{star}")

# ② pilot 逐 spot（2 张切片）
q = f"{R}/pilot/plip_pilot_spot_scores.csv.gz"
if os.path.exists(q):
    s = pd.read_csv(q)
    print(f"\n[plip_pilot_spot_scores.gz]  {len(s)} spot / {s.gsm.nunique()} 切片")
    for c in ["lung adenocarcinoma", "normal lung tissue", "lung tissue with adenocarcinoma in situ"]:
        if c in s.columns:
            print(f"  corr(tissue_frac, {c:44s}) = {s[c].corr(s.tissue_frac):+.4f}")

print("\n" + "=" * 72)
print("⇒ 可直接引用的数：corr_lesion_vs_tissuefrac = 0.754（P4 四切片）")
print("  **0.825 不得写成「与组织占比的相关」**——它是 |corr_lesion_vs_empty| 的 P4_Normal 单切片值")
print("=" * 72)
