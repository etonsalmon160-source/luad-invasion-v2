#!/usr/bin/env python3
# 50_fig_signature_qc.py —— 逆转签名的预后质控：疾病方向 × 预后方向 的四象限
#   规则（63_qc_prognosis_filter.py）：疾病升高(lfc>0)⇒拟敲低⇒需 HR>1；疾病降低⇒拟补回⇒需 HR<1；
#   方向与意图相反（保护性程序）或没有 HR ⇒ 置零剔除。
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/paper_figures"
D = pd.read_csv(f"{ROOT}/results/10_niche/tr_singlecell/global_paired_qc.tsv", sep="\t")
D = D[D.HR.isna() | (D.HR > 0)]
keep = D[D.keep]
drop_dir = D[(~D.keep) & D.HR.notna()]        # 方向与意图相反（保护性程序）
drop_na  = D[(~D.keep) & D.HR.isna()]         # TCGA 里没有 HR
drop = D[~D.keep]
print(f"全表 {len(D):,}｜保留 {len(keep):,}（{len(keep)/len(D)*100:.1f}%）"
      f"｜剔除 {len(drop):,}（{len(drop)/len(D)*100:.1f}%）"
      f"  = 方向反 {len(drop_dir):,} + 无 HR {len(drop_na):,}")
print(f"  HR>1: {int((D.HR>1).sum()):,}   HR<1: {int((D.HR<1).sum()):,}")

KEEP_C, DROP_C = "#2F5597", "#C9CCD1"
fig = plt.figure(figsize=(7.0, 4.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.34], wspace=0.06,
                      left=0.075, right=0.985, top=0.84, bottom=0.135)
ax = fig.add_subplot(gs[0, 0]); axr = fig.add_subplot(gs[0, 1])

XL = 0.45
ax.axvline(0, color="0.65", lw=0.7, zorder=1)
ax.axhline(1, color="0.65", lw=0.7, zorder=1)
for m, c, lb in ((~D.keep, DROP_C, "dropped"), (D.keep, KEEP_C, "kept")):
    s = D[m]
    ax.scatter(s.lfc, s.HR, s=2.6, c=c, alpha=0.55, linewidths=0, rasterized=True,
               zorder=3 if lb == "kept" else 2, label=lb)
ax.set_yscale("log")
ax.set_xlim(-XL, XL); ax.set_ylim(0.28, 3.2)
ax.set_yticks([0.4, 0.7, 1.0, 1.5, 2.5]); ax.set_yticklabels(["0.4", "0.7", "1.0", "1.5", "2.5"])
ax.set_xlabel("Disease axis   log2 fold-change  (IAC − precursor)", fontsize=7.2)
ax.set_ylabel("Adjusted hazard ratio (TCGA-LUAD)", fontsize=7.2)
ax.tick_params(labelsize=6.5, length=2.5)
for s in ax.spines.values():
    s.set_linewidth(0.6); s.set_color("0.45")

# 象限说明
Q = [(0.435, 3.05, "disease↑ · harmful\n→ knock down   KEPT", KEEP_C, "right", "top"),
     (0.435, 0.305, "disease↑ · protective\n→ dropped", "0.55", "right", "bottom"),
     (-0.435, 0.305, "disease↓ · protective\n→ restore   KEPT", KEEP_C, "left", "bottom"),
     (-0.435, 3.05, "disease↓ · harmful\n→ dropped", "0.55", "left", "top")]
for x, y, t, c, ha, va in Q:
    ax.text(x, y, t, fontsize=6.0, color=c, ha=ha, va=va, linespacing=1.5,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.80, pad=1.2))
LAB = ["XBP1", "FOXP3", "NKX2-1", "SFTPC", "TCF21", "DKK1", "SFTPA1"]
sub = D[D.gene.isin(LAB)]
OFF = {"XBP1": (6, -10), "FOXP3": (6, 6), "NKX2-1": (-6, -10), "SFTPC": (6, -10),
       "TCF21": (6, 6), "DKK1": (6, 6), "SFTPA1": (6, -10)}
for _, r in sub.iterrows():
    c = KEEP_C if r.keep else "#8A8F98"
    dx, dy = OFF.get(r.gene, (6, 6))
    ax.annotate(r.gene, (r.lfc, r.HR), xytext=(dx, dy), textcoords="offset points",
                fontsize=6.0, color=c, ha="left" if dx > 0 else "right", va="center",
                fontweight="bold", zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.7))
hs, ls = ax.get_legend_handles_labels()
ax.set_title("the signature is filtered on independent survival evidence",
             fontsize=8.0, fontweight="bold", pad=7, loc="left")

# 右：占比
axr.bar([0], [len(keep)], width=0.55, color=KEEP_C, linewidth=0, label="kept")
axr.bar([1], [len(drop_dir)], width=0.55, color=DROP_C, linewidth=0,
        label="dropped: direction opposite")
axr.bar([1], [len(drop_na)], bottom=[len(drop_dir)], width=0.55, color="#E4E6E9",
        linewidth=0, label="dropped: no survival data")
for i, v, b in ((0, len(keep), 0), (1, len(drop_dir), 0), (1, len(drop_na), len(drop_dir))):
    axr.text(i + 0.30, b + v / 2, f"{v:,}", ha="left", va="center", fontsize=5.8, color="0.30")
axr.text(0, len(keep) + 900, f"{len(keep)/len(D)*100:.0f}%", ha="center", fontsize=7.0,
         color=KEEP_C, fontweight="bold")
axr.text(1, len(drop_dir) + len(drop_na) + 900, f"{len(drop)/len(D)*100:.0f}%", ha="center",
         fontsize=7.0, color="#8A8F98", fontweight="bold")
axr.set_xlim(-0.55, 1.75)
axr.set_xticks([0, 1]); axr.set_xticklabels(["kept", "dropped"], fontsize=6.4)
axr.legend(loc="upper left", bbox_to_anchor=(0.30, 1.0), frameon=False, fontsize=5.4,
           handlelength=0.9, labelspacing=0.22, handletextpad=0.35)
axr.set_ylim(0, len(keep) * 1.22); axr.set_yticks([])
axr.set_ylabel("genes", fontsize=6.6)
for s in ("top", "right"):
    axr.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axr.spines[s].set_linewidth(0.6); axr.spines[s].set_color("0.45")
axr.set_title("how much is removed", fontsize=8.0, fontweight="bold", pad=7, loc="left")

fig.legend(hs, ls, loc="lower center", bbox_to_anchor=(0.42, -0.035), ncol=2, frameon=False,
           fontsize=6.2, markerscale=3.4, handletextpad=0.3, columnspacing=1.6)
fig.savefig(f"{OUT}/P20_signature_qc.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P20_signature_qc.png", dpi=400, bbox_inches="tight")
print("✓ P20_signature_qc")
print(sub[["gene", "lfc", "HR", "keep"]].round(3).to_string(index=False))
