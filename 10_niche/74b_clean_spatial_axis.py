#!/usr/bin/env python3
# 74b_clean_spatial_axis.py —— 空转轴的去污染步（**事后补写的脚本**，见下）
#
# 背景：74_spatial_global_axis.py 产出 spatial_global_depthmatched_lfc.tsv（18,087 基因），
#       75_spatial_gess.R 用的却是 spatial_global_depthmatched_clean_lfc.tsv（18,060 基因）。
#       生成后者的代码当时**没有落成脚本**（临时命令），2026-10-06 从两个文件反推补写本脚本。
#
# 反推结论：clean 相对 depthmatched **只少 27 个基因、不改任何数值**
#   （两条轴在共同基因上 Spearman = 1.00000）
#
# 排除的三类：
#   ① 免疫球蛋白基因 —— 浆细胞/血清 Ig 污染，会主导差异
#   ② 线粒体基因（MT-）—— 测序质量 / 细胞状态，与疾病轴无关
#   ③ Visium 探针表里的 **抗体行**（mouse_IgG* / rat_IgG*）—— 量级 1e6 的技术行
#
# ⚠️ 已知**过度剔除**：规则 ① 的前缀 ^IGL 会一并命中 **IGLON5**，
#    它是 IgLON 家族黏附分子，**不是免疫球蛋白**。当时确实把它删了。
#    本脚本为保持与既有产物逐位一致，仍按原规则执行，但单独打印提示。
import re, sys
import pandas as pd

S = "/home/eto/luad_v2/results/10_niche/sigsearch"
SRC = f"{S}/spatial_global_depthmatched_lfc.tsv"
DST = f"{S}/spatial_global_depthmatched_clean_lfc.tsv"

PAT = [("免疫球蛋白", re.compile(r"^IG[HKL]")),
       ("线粒体",     re.compile(r"^MT-")),
       ("抗体探针行", re.compile(r"^(mouse|rat)_IgG"))]

df = pd.read_csv(SRC, sep="\t")
print(f"输入 {SRC}: {len(df)} 基因")
keep = pd.Series(True, index=df.index)
for name, pat in PAT:
    hit = df.gene.str.match(pat)
    print(f"  {name:8s} 剔除 {int(hit.sum()):3d} 个")
    keep &= ~hit
out = df[keep].copy()
print(f"输出 {len(out)} 基因（剔除 {len(df) - len(out)}）")

over = sorted(set(df.gene[~keep]) - set(df.gene[df.gene.str.match(PAT[0][1])]))
extra = sorted(set(out.gene) & {"IGLON5"})
if extra:
    print(f"⚠️  注意：{'、'.join(extra)} 是被 ^IGL 前缀误伤的黏附分子（非免疫球蛋白），"
          "此处保留原行为以便与既有产物逐位一致；若要纠正需重跑 75_spatial_gess.R")

# 自检：写盘**之前**先与既有产物比对（写后再读会被自己覆盖）
import os
import numpy as np
ref = pd.read_csv(DST, sep="\t") if os.path.exists(DST) else None
if ref is not None:
    assert list(ref.gene) == list(out.gene), "基因或顺序与既有产物不一致"
    # 既有产物落盘时小数位比源文件少 ⇒ 逐位比较会有 ~1e-21 的纯文本精度差（机器精度级），
    # 用相对容差 1e-12 比，与"数值一致"等价。
    assert np.allclose(ref.lfc.values, out.lfc.values, rtol=1e-12, atol=1e-15), \
        "数值与既有产物不一致"
    d = np.abs(ref.lfc.values - out.lfc.values).max()
    print(f"✓ 自检通过：与既有产物一致（{len(out)} 基因，最大差 {d:.2e}＝文本精度）")
else:
    print("（无既有产物可比对，跳过自检）")

out.to_csv(DST, sep="\t", index=False)
print(f"落盘 {DST}")
