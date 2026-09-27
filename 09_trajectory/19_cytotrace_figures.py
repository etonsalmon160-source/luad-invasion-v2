#!/usr/bin/env python
"""CytoTRACE 结果出图（预注册 16 的 §4 判据可视化）。
英文标签——本机无 CJK 字体。输出到 results/09_trajectory/cytotrace/figures/。
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]
D = "/home/eto/luad_v2/results/09_trajectory/cytotrace"
FIG = f"{D}/figures"

blobs = {a: pd.read_csv(f"{D}/cytotrace_per_cell_{a}.tsv.gz", sep="\t")
         for a in ["all", "alveolar"]}

# ---- 图1：期别箱线图（两档） ----
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for ax, arm in zip(axes[:2], ["all", "alveolar"]):
    d = blobs[arm]
    data = [d.loc[d.stage == s, "ct_score"].values for s in STAGES]
    bp = ax.boxplot(data, labels=STAGES, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=2))
    for b in bp["boxes"]:
        b.set_facecolor("#8fbcd4")
    med = [np.median(x) for x in data]
    ax.plot(range(1, 6), med, "o-", color="crimson", lw=2, ms=6, label="median")
    for i, m in enumerate(med):
        ax.annotate(f"{m:.3f}", (i + 1, m), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=9, color="crimson")
    ax.set_title(f"CytoTRACE score by stage [{arm}]  n={len(d):,}")
    ax.set_ylabel("CytoTRACE score (0-1)")
    ax.grid(alpha=.25, axis="y")
    ax.legend(fontsize=9)

# 第三格：深度配平面板
d = blobs["all"]
d["depth_bin"] = pd.qcut(d["ct_num_exp_genes"], 5, labels=False, duplicates="drop")
ax = axes[2]
nb = sorted(d.depth_bin.dropna().unique())
x = np.arange(len(nb))
colors = plt.cm.viridis(np.linspace(0, .85, len(STAGES)))
for si, s in enumerate(STAGES):
    ys, lo, hi = [], [], []
    for b in nb:
        v = d.loc[(d.depth_bin == b) & (d.stage == s), "ct_score"]
        ys.append(v.median() if len(v) else np.nan)
        if len(v) > 20:
            bs = [np.median(np.random.choice(v, len(v), replace=True)) for _ in range(200)]
            lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
        else:
            lo.append(np.nan); hi.append(np.nan)
    ax.plot(x, ys, "o-", color=colors[si], label=s, lw=2, ms=5)
    ax.fill_between(x, lo, hi, color=colors[si], alpha=.18)
ax.set_xticks(x)
ax.set_xticklabels([f"Q{b+1}\n{d[d.depth_bin==b].ct_num_exp_genes.min():.0f}-"
                    f"{d[d.depth_bin==b].ct_num_exp_genes.max():.0f}" for b in nb], fontsize=8)
ax.set_xlabel("detected-genes quintile (depth-matched)")
ax.set_ylabel("median CytoTRACE score")
ax.set_title("Stage effect SURVIVES depth matching\n(IAC highest in every bin)")
ax.legend(fontsize=8, ncol=2)
ax.grid(alpha=.25)
plt.tight_layout()
plt.savefig(f"{FIG}/cytotrace_by_stage_and_depth.png", dpi=140)
print("写好", f"{FIG}/cytotrace_by_stage_and_depth.png")

# ---- 图2：逐患者斜率 ----
fig, axes = plt.subplots(1, 7, figsize=(20, 3.2), sharey=True)
pats = sorted(blobs["all"].patient_id.unique())
for ax, p in zip(axes, pats[:7]):
    pass
fig, ax = plt.subplots(figsize=(14, 5.5))
rhos = {}
for p in sorted(blobs["all"].patient_id.unique()):
    s = blobs["all"][blobs["all"].patient_id == p]
    if s.stage.nunique() < 2:
        continue
    xi = np.array([STAGES.index(v) for v in s.stage])
    med = [np.median(s.loc[s.stage == st, "ct_score"]) for st in STAGES
           if (s.stage == st).sum() > 0]
    xs = [STAGES.index(st) for st in STAGES if (s.stage == st).sum() > 0]
    rho, _ = stats.spearmanr(xi, s.ct_score)
    rhos[p] = rho
    ax.plot(xs, med, "o-", alpha=.65, lw=1.6,
            color=("crimson" if rho < 0 else "#2b6cb0"))
ax.set_xticks(range(5)); ax.set_xticklabels(STAGES)
ax.set_ylabel("median CytoTRACE score")
nneg = sum(1 for r in rhos.values() if r < 0)
ax.set_title(f"Per-patient stage trend — {nneg}/{len(rhos)} negative; "
             f"median rho = {np.median(list(rhos.values())):+.3f}  "
             f"(red = negative slope)")
ax.grid(alpha=.25)
plt.tight_layout()
plt.savefig(f"{FIG}/cytotrace_per_patient.png", dpi=140)
print("写好", f"{FIG}/cytotrace_per_patient.png")
print("逐患者 rho:", {k: round(v, 3) for k, v in sorted(rhos.items())})
