#!/usr/bin/env python3
# 79_split_half_null.py —— 补落的脚本：劈半稳定性的**零模型**
#
# 背景：STATUS 记「观测 ρ 中位 0.962 vs 零模型 −0.155」，但生成脚本未落盘。
#       2026-10-06 审计：**两种零模型定义给出完全不同的数**，必须二选一并写明。
import numpy as np, pandas as pd
S = "/home/eto/luad_v2/results/10_niche/sigsearch/splits"
RNG = np.random.default_rng(20261006)
rows = []
for k in range(4):
    g = lambda t: pd.read_csv(f"{S}/cor_split{k}_{t}.tsv", sep="\t").groupby("pert").cor_score.median()
    h1, h2, n1 = g("h1"), g("h2"), g("h1_null")
    j = pd.DataFrame({"h1": h1, "h2": h2, "n": n1}).dropna()
    obs = j.h1.corr(j.h2, method="spearman")
    # 零模型 A：仓库自带的 *_null.tsv（一个"空查询"的打分）
    nulA = j.h1.corr(j.n, method="spearman")
    # 零模型 B：置换——打乱一半的排名
    nulB = np.median([j.h1.corr(pd.Series(RNG.permutation(j.h2.values), index=j.index),
                                method="spearman") for _ in range(300)])
    rows.append(dict(split=k, obs=obs, nullA_file=nulA, nullB_perm=nulB))
R = pd.DataFrame(rows)
print(R.round(4).to_string(index=False))
print(f"\n观测 ρ 中位        : {R.obs.median():.3f}")
print(f"零模型 A（自带文件）中位: {R.nullA_file.median():+.3f}   ← STATUS 记的是 −0.155，复现不出")
print(f"零模型 B（置换）    中位: {R.nullB_perm.median():+.3f}   ← 本脚本推荐")
print("""
🔴 为什么不能用零模型 A：`cor_*_null.tsv` 与真实打分在 (pert, cell) 层面
   自带 ±0.5–0.6 的相关（库结构造成），并非"无关对照"。
   ⇒ 引用时必须写清用的是哪一种；0.962 对面写「−0.155」属口径不明。
""")
