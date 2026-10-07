#!/usr/bin/env python3
# 42_fig_compound_reversal.py —— 化合物层主结果（单细胞疾病轴）
#   三法并排（CMap / Cor-Spearman / Cor-Pearson），逐化合物取**跨细胞系中位**分，
#   各取 top100 后交集 = 34 个候选（与 TARGET_REVERSAL_STATUS.md §三-5 一致）
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = "/home/eto/luad_v2"; S = f"{ROOT}/results/10_niche/sigsearch"
OUT = f"{ROOT}/results/paper_figures"

# ── 三法逐化合物中位分 ──
AGG = {}
for f, col in [("res_CMAP.tsv", "raw_score"), ("res_COR_spearman.tsv", "cor_score"),
               ("res_COR_pearson.tsv", "cor_score")]:
    AGG[f] = pd.read_csv(f"{S}/{f}", sep="\t").groupby("pert")[col].median()
M = pd.DataFrame(AGG); M.columns = ["CMAP", "SPEAR", "PEAR"]
TOPS = {c: set(M.nsmallest(100, c).index) for c in M.columns}
COMMON = sorted(set.intersection(*TOPS.values()))
CAND = M.loc[COMMON]
print(f"三法 top100 交集 = {len(COMMON)} 个化合物")

# ── 机制归类（MOAss 字段不全 ⇒ 用 STATUS/草稿里已核过的那张表）──
MOA = {
    "PI3K / mTOR": ["torin-1", "torin-2", "AZD-2014", "AZD-8055", "GSK-2126458", "GDC-0980",
                     "MLN-0128", "PKI-179", "NVP-BEZ235", "OSI-027", "WYE-125132", "PI-103",
                     "GDC-0941", "taselisib", "voxtalisib", "GSK-2110183"],
    "HSP90": ["NVP-AUY922", "tanespimycin", "AT-13387"],
    "Proteasome": ["MG-132", "bortezomib"],
    "Cell-cycle kinase": ["XL-888", "dinaciclib", "PF-03814735", "SB-939"],
    "Other / unannotated": ["BIIB-021", "NVP-TAE226", "TAK-285", "JW-7-24-1", "YM-155",
                            "lacidipine", "midostaurin", "pentobarbital", "okadaic acid"],
}
g2f = {g: fam for fam, gs in MOA.items() for g in gs}
missing = [g for g in COMMON if g not in g2f]
assert not missing, f"未归类的候选: {missing}"
assert sum(len(v) for v in MOA.values()) == 34, "机制表数量不是 34"

FAMC = {"PI3K / mTOR": "#C0392B", "HSP90": "#E08A3C", "Proteasome": "#7E4FA8",
        "Cell-cycle kinase": "#2E8B8B", "Other / unannotated": "#9aa0a6"}
ORDER = ["PI3K / mTOR", "HSP3S0" if False else "HSP90", "Proteasome", "Cell-cycle kinase",
         "Other / unannotated"]

fig = plt.figure(figsize=(8.0, 3.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.72], wspace=0.28,
                      left=0.065, right=0.985, top=0.86, bottom=0.145)

# ── A：秩-秩散点（log 轴）；方框 = top100 ──
RK = M.rank()                       # 越小越靠前（最负分）
ax = fig.add_subplot(gs[0, 0])
ax.add_patch(plt.Rectangle((1, 1), 99, 99, facecolor="#f6f1ea", edgecolor="none", zorder=0))
ax.axvline(100, color="0.72", lw=0.6, ls=(0, (3, 2.5)), zorder=1)
ax.axhline(100, color="0.72", lw=0.6, ls=(0, (3, 2.5)), zorder=1)
ax.scatter(RK.CMAP, RK.SPEAR, s=2.6, c="0.78", linewidths=0, rasterized=True, zorder=2)
for fam in ORDER:
    sub = RK.loc[[g for g in MOA[fam]]]
    ax.scatter(sub.CMAP, sub.SPEAR, s=16, c=FAMC[fam], edgecolor="white",
               linewidths=0.35, zorder=3)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.8, 2100); ax.set_ylim(0.8, 2100)
ax.set_xticks([1, 10, 100, 1000]); ax.set_yticks([1, 10, 100, 1000])
ax.set_xticklabels(["1", "10", "100", "1000"]); ax.set_yticklabels(["1", "10", "100", "1000"])
ax.set_xlabel("rank by CMap  τ   (log)", fontsize=7)
ax.set_ylabel("rank by Cor-Spearman  (log)", fontsize=7)
ax.tick_params(labelsize=6.2, length=2.2)
for s in ax.spines.values():
    s.set_linewidth(0.6); s.set_color("0.45")
exp = 100 * 100 / len(M)
ax.set_title("same library — CMap versus correlation",
             fontsize=7.4, fontweight="bold", pad=5, loc="left")
ax.text(0.985, 0.035, f"shaded box = top 100 in both\n34 compounds fall inside all three top-100 lists\n"
        f"(chance expectation ≈ {exp:.1f})",
        transform=ax.transAxes, fontsize=5.6, color="0.30", va="bottom", ha="right",
        linespacing=1.55,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.6))

# ── B：34 个的机制构成 ──
axb = fig.add_subplot(gs[0, 1])
cnt = pd.Series({f: len(MOA[f]) for f in ORDER})
axb.barh(np.arange(len(ORDER))[::-1], cnt.values,
         color=[FAMC[f] for f in ORDER], height=0.62, linewidth=0)
for i, f in enumerate(ORDER):
    y = len(ORDER) - 1 - i
    axb.text(cnt[f] + 0.35, y, str(cnt[f]), va="center", fontsize=6.6, color="0.25")
axb.set_yticks(np.arange(len(ORDER))[::-1]); axb.set_yticklabels(ORDER, fontsize=6.5)
axb.set_xlim(0, 19); axb.set_xticks([0, 5, 10, 15]); axb.tick_params(labelsize=6, length=2)
axb.set_xlabel("compounds", fontsize=7)
for s in ("top", "right", "left"):
    axb.spines[s].set_visible(False)
axb.spines["bottom"].set_linewidth(0.6); axb.spines["bottom"].set_color("0.45")
axb.set_title("mechanism of the 34 candidates", fontsize=7.4, fontweight="bold",
              pad=5, loc="left")
axb.grid(axis="x", color="0.92", linewidth=0.5); axb.set_axisbelow(True)

fig.text(0.065, 0.955, "Compound-level reversal of the single-cell disease axis",
         fontsize=9, fontweight="bold", ha="left")
fig.savefig(f"{OUT}/P12_compound_reversal_singlecell.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/P12_compound_reversal_singlecell.png", dpi=400, bbox_inches="tight")
print("✓ P12_compound_reversal_singlecell")
print(CAND.round(4).to_string())
