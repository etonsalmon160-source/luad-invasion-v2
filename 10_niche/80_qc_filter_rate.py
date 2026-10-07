#!/usr/bin/env python3
# 80_qc_filter_rate.py —— 补落的脚本：预后质控的**剔除率分母**
#
# 背景：草稿 §14.2(c) 与 STATUS 记「全局轴 46% 被剔」，但未注明分母。
#       2026-10-06 审计：三个分母给出三个数，必须写明。
import pandas as pd
D = pd.read_csv("/home/eto/luad_v2/results/10_niche/tr_singlecell/global_paired_qc.tsv", sep="\t")
n_all = len(D)
n_hr = int(D.HR.notna().sum())
n_keep = int(D.keep.sum())
n_drop = int((~D.keep).sum())
n_dir = int(((~D.keep) & D.HR.notna()).sum())   # 方向与意图相反
n_na = int(((~D.keep) & D.HR.isna()).sum())     # TCGA 无 HR
print(f"全表基因            : {n_all:,}")
print(f"  有 HR             : {n_hr:,}   无 HR: {n_na:,}")
print(f"  保留              : {n_keep:,} ({n_keep/n_all*100:.1f}%)")
print(f"  剔除              : {n_drop:,} ({n_drop/n_all*100:.1f}%)")
print(f"     ├ 方向相反      : {n_dir:,}")
print(f"     └ 无 HR         : {n_na:,}")
print(f"""
🔴 三个分母三个数：
   全表            {n_drop/n_all*100:.1f}%
   仅"有 HR"子集    {n_drop/n_hr*100:.1f}%      ← 若语境是"在能判的基因里"就用这个
   仅"方向相反"占全表 {n_dir/n_all*100:.1f}%      ← 若语境是"被生存证据否决的"
   ⇒ 引用时**必须写明分母**；草稿记的 46% 与三者都不完全对应。
""")
