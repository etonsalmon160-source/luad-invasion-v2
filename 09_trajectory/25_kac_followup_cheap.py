#!/usr/bin/env python
"""
发育谱臂 §7 任务表里的两个**便宜任务**（PLAN_AND_CHECKPOINTS.md §7 D1、D2）。

两者都**不改口径、不重跑打分、不产新判据**，只读已入库的产物做查表与统计。
⇒ 按 §7 登记，**不需签字**；结论标注为"探索（未注册）"。

  D1  **MP9 为什么反向**：D2 报 MP6 vs KAC 签名 ρ=+0.722、MP9 vs KAC 签名 ρ=−0.376。
      查 MP9 那 50 个基因是什么、MP9 高分落在哪些簇、是否与 KAC 签名共享基因。

  D2  **簇 16/21/23 是不是单一患者的亚克隆**：K3 已标 top_patient_frac 0.561/0.673/0.862。
      查这三簇的患者构成 × 该患者的分期构成，并核对三簇是否同一患者驱动。

产物：`results/09_trajectory/paper_axis/followup_cheap.md`（人读）+ 标准输出。
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

D = "/home/eto/luad_v2/results/09_trajectory/paper_axis"
WT = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963"
OUT = f"{D}/followup_cheap.md"

sys.path.insert(0, f"{WT}/05_annotation")
import epi_subtype_panel as P                                    # noqa: E402

PANEL_MP = P.build_panel("paper_mp")
PANEL_KAC = P.build_panel("kac_sig")["KAC"]
KAC_SIDE = {"Tumor cell/KAC", "KAC/inflammatory"}
STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]


def colname(subtype):
    return "s_" + "".join(c if c.isalnum() else "_" for c in subtype)


cm = pd.read_csv(f"{D}/cluster_scores_papermp.csv")
pc = pd.read_csv(f"{D}/per_cell_scores.tsv.gz", sep="\t")
print(f"簇 {len(cm)} / 细胞 {len(pc):,}")

lines = ["# 发育谱臂 §7 便宜任务（D1 / D2）—— 探索（未注册）", ""]

# ══════════════════════════════════════════════════════════════════
# D1  MP9 为什么反向
# ══════════════════════════════════════════════════════════════════
lines += ["## D1 · MP9（`KAC/inflammatory`）为什么与 KAC 签名反向", ""]

mp9 = PANEL_MP["KAC/inflammatory"]
mp6 = PANEL_MP["Tumor cell/KAC"]
lines += [
    f"MP9 共 {len(mp9)} 基因；MP9 ∩ MP6 = {sorted(set(mp9) & set(mp6))}；"
    f"MP9 ∩ Table S3 KAC 签名 = {sorted(set(mp9) & set(PANEL_KAC))}",
    "",
    "MP9 的 50 个基因（逐字）：",
    "",
    "```",
    ", ".join(mp9),
    "```",
    "",
]
# 与其它 MP 的重叠（若与 AT2/Ciliated/Club 重叠高，反向可能是共享基因造成的）
lines += ["MP9 与其余 7 个 MP 的基因重叠：", "", "| MP | 重叠数 | 基因 |", "| :--- | :---: | :--- |"]
for t, gl in PANEL_MP.items():
    if t == "KAC/inflammatory":
        continue
    ov = sorted(set(mp9) & set(gl))
    lines.append(f"| {t} | {len(ov)} | {', '.join(ov) if ov else '—'} |")
lines.append("")

c9, c6, ck = colname("KAC/inflammatory"), colname("Tumor cell/KAC"), "s_KACsig"
m9 = cm[["cluster", c9, c6, ck, "is_kac_side", "argmax_papermp", "n_cells"]].copy()
m9["mp9_rank"] = m9[c9].rank(ascending=False).astype(int)
m9 = m9.sort_values(c9, ascending=False)

lines += [
    "### 27 簇的 MP9 分数（降序）",
    "",
    "| 排名 | 簇 | MP9 | MP6 | KAC 签名 | D1 是否 KAC 侧 | 该簇 argmax |",
    "| :---: | :---: | ---: | ---: | ---: | :---: | :--- |",
]
for _, r in m9.iterrows():
    lines.append(f"| {r.mp9_rank} | c{int(r.cluster)} | {r[c9]:+.3f} | {r[c6]:+.3f} | "
                 f"{r[ck]:+.3f} | {'✅' if r.is_kac_side else ''} | {r.argmax_papermp} |")
lines.append("")

kc = m9[m9.is_kac_side]
rest = m9[~m9.is_kac_side]
lines += [
    f"- 8 个 KAC 侧簇的 MP9 中位 **{kc[c9].median():+.3f}**（排名 {kc.mp9_rank.min()}–{kc.mp9_rank.max()}）；"
    f"其余 19 簇中位 **{rest[c9].median():+.3f}**",
    f"- KAC 侧簇的 MP9 排名：{sorted(kc.mp9_rank.tolist())}（27 取 8）"
    f" ⇒ 若全挤在后段（大排名），说明 MP9 在 KAC 侧被**系统性压低**",
    f"- 簇级相关：corr(MP9, MP6) = {stats.spearmanr(cm[c9], cm[c6])[0]:+.3f}；"
    f"corr(MP9, KAC签名) = {stats.spearmanr(cm[c9], cm[ck])[0]:+.3f}",
    "",
]
# MP9 最高 / 最低的簇各自是谁
top = m9.head(5)
bot = m9.tail(5)
lines += [
    f"- MP9 最高 5 簇：{', '.join(f'c{int(r.cluster)}({r[c9]:+.2f}, {r.argmax_papermp})' for _, r in top.iterrows())}",
    f"- MP9 最低 5 簇：{', '.join(f'c{int(r.cluster)}({r[c9]:+.2f}, {r.argmax_papermp})' for _, r in bot.iterrows())}",
    "",
]

# 逐细胞层面再核一次（避免簇均值的合并假象）
r_cell = stats.spearmanr(pc[c9], pc[ck])[0]
r_cell66 = stats.spearmanr(pc[c9], pc[c6])[0]
lines += [
    f"- **逐细胞层面**（{len(pc):,} 个）：corr(MP9, KAC签名) = {r_cell:+.3f}；"
    f"corr(MP9, MP6) = {r_cell66:+.3f}",
    "",
]

# ══════════════════════════════════════════════════════════════════
# D2  簇 16/21/23 是不是单一患者的亚克隆
# ══════════════════════════════════════════════════════════════════
lines += ["---", "", "## D2 · 簇 16 / 21 / 23 是不是单一患者的亚克隆", ""]

FOCUS = [16, 21, 23]
cohort_stage = pc["stage"].value_counts(normalize=True).reindex(STAGES) * 100

lines += ["### 这三簇的患者构成", "",
          "| 簇 | 细胞数 | 患者数 | 最大患者 | 占比 | 第 2 患者 | 占比 |",
          "| :---: | ---: | :---: | :--- | ---: | :--- | ---: |"]
detail = {}
for c in FOCUS:
    g = pc[pc["cluster"] == c]
    vc = g["patient_id"].value_counts()
    detail[c] = g
    p1, p2 = vc.index[0], (vc.index[1] if len(vc) > 1 else "—")
    f2 = (vc.iloc[1] / len(g) * 100) if len(vc) > 1 else 0.0
    lines.append(f"| c{c} | {len(g):,} | {len(vc)} | {p1} | {vc.iloc[0]/len(g):.1%} | {p2} | {f2:.1f}% |")
lines.append("")

lines += ["### 每个簇的 top 患者在**全部期别**上的分布（看是不是就一个期）", "",
          "| 簇 | 主导患者 | 该患者在本簇的期别分布 | 该患者在全队列的期别分布 |",
          "| :---: | :--- | :--- | :--- |"]
for c in FOCUS:
    g = detail[c]
    tp = g["patient_id"].value_counts().index[0]
    gs = g[g.patient_id == tp]["stage"].value_counts(normalize=True).reindex(STAGES).fillna(0) * 100
    allp = pc[pc.patient_id == tp]["stage"].value_counts(normalize=True).reindex(STAGES).fillna(0) * 100
    fmt = lambda s: " / ".join(f"{st} {s[st]:.0f}%" for st in STAGES if s[st] > 0) or "—"
    lines.append(f"| c{c} | {tp} | {fmt(gs)} | {fmt(allp)} |")
lines.append("")

lines += [f"（参考）全队列期别构成：{' / '.join(f'{st} {cohort_stage[st]:.1f}%' for st in STAGES)}", ""]

# 三簇是否同一患者驱动
tops = {c: detail[c]["patient_id"].value_counts().index[0] for c in FOCUS}
lines += [
    f"三簇的 top 患者：{', '.join(f'c{c} → {p}' for c, p in tops.items())}"
    f" ⇒ **{'同一个患者' if len(set(tops.values())) == 1 else '不同患者'}**",
    "",
]

# 三簇合计：一个患者最多能解释多少
allf = pc[pc["cluster"].isin(FOCUS)]
vc_all = allf["patient_id"].value_counts()
lines += [
    f"三簇合计 {len(allf):,} 细胞；患者构成前 5：",
    "",
    "| 患者 | 细胞 | 占比 |",
    "| :--- | ---: | ---: |",
]
for p, n in vc_all.head(5).items():
    lines.append(f"| {p} | {n:,} | {n/len(allf):.1%} |")
lines.append("")

# 这三簇的患者是否也贡献了其它簇（若是，则"专属亚克隆"更弱）
tp_all = set(tops.values())
lines += ["### 这些患者是否只出现在这三簇里（越专属 ⇒ 越像私有亚克隆）", "",
          "| 患者 | 该患者总细胞 | 落在 c16/21/23 | 占其自身 | 该患者还覆盖多少其它簇 |",
          "| :--- | ---: | ---: | ---: | ---: |"]
for p in sorted(tp_all):
    gp = pc[pc.patient_id == p]
    in_f = int(gp["cluster"].isin(FOCUS).sum())
    n_other = gp[~gp["cluster"].isin(FOCUS)]["cluster"].nunique()
    lines.append(f"| {p} | {len(gp):,} | {in_f:,} | {in_f/len(gp):.1%} | {n_other} |")
lines.append("")

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\n写出 {OUT}")
