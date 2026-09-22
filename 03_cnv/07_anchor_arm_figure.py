#!/usr/bin/env python3
"""「锚定 + 同质」臂 —— P11 冒烟结果（只读已入库产物，不重算 CopyKAT）。

要回答的问题（ANCHOR_PREREG.md §1）：
  显式锚定 + 单亚型输入**同时**满足时，CopyKAT 的判决方向对不对、稳不稳？

图回答三件事：
  (A) 锚定的作用：**同一批 728 个核**，锚定关 vs 锚定开 —— 病灶侧 9.8% → 30.6%，
      锚定侧两次都是 0%。⚠️ "锚定关"那一次不是登记过的臂，是首跑因实现 bug
      意外没把锚定传进去而留下的对照，本图**标注为意外对照**，不得当预注册臂引用。
  (B) 两条预注册运行的分侧结果 + G1 的 0.10 / 0.50 停机线。
  (C) 细胞数级联：840 覆盖度地板 + 用户签字集 + 亚型三层收窄各砍掉多少。
  (D) 诚实边界：P11 一对是冒烟，不是"方法可用"的结论。

无 CJK 字体 ⇒ 全英文标签。
"""
import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CNV = os.path.join(ROOT, "results", "03_cnv")
FIG = os.path.join(ROOT, "figures")
OUT = os.path.join(CNV, "anchor_at2")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_runs():
    runs = {}
    for s in ("P11_LUAD", "P11_AAH"):
        with open(os.path.join(OUT, f"{s}.json")) as fh:
            runs[s] = json.load(fh)
    return runs


def accidental_control():
    """首跑（锚定被 bug 清零）留下的 P11_LUAD 判定。

    ⚠️ **不是预注册臂**，是意外对照。放在 `results/` 下（不是 `scratch/`，
    后者被 gitignore）才能让本图 (A) 的 "anchor OFF" 那根柱子可被独立复核。
    """
    p = os.path.join(ROOT, "results", "03_cnv", "anchor_firstrun_UNPLANNED_control",
                     "P11_LUAD_copykat_prediction.txt")
    if not os.path.exists(p):
        return None
    anc = les = aneu_les = aneu_anc = 0
    with open(p) as fh:
        next(fh)
        for line in fh:
            name, _, pred = line.rstrip("\n").partition("\t")
            is_anc = "P11_Normal" in name
            if is_anc:
                anc += 1; aneu_anc += (pred == "aneuploid")
            else:
                les += 1; aneu_les += (pred == "aneuploid")
    return dict(source=os.path.relpath(p, ROOT), source_sha256=sha256(p),
                n_lesion=les, n_aneuploid_lesion=aneu_les,
                frac_lesion=aneu_les / les, n_anchor=anc,
                n_aneuploid_anchor=aneu_anc, frac_anchor=aneu_anc / anc)


def main():
    os.makedirs(FIG, exist_ok=True)
    runs = load_runs()
    acc = accidental_control()

    fig = plt.figure(figsize=(19.5, 8.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.05, 1.25],
                          height_ratios=[1.0, 0.85], hspace=0.42, wspace=0.30)

    # ---- (A) 锚定的作用：同一批 728 个核 ----
    ax = fig.add_subplot(gs[0, 0])
    lu = runs["P11_LUAD"]
    conds = ["anchor OFF\n(first run,\nUNPLANNED)", "anchor ON\n(pre-registered)"]
    les_v = [acc["frac_lesion"], lu["frac_aneuploid_lesion"]] if acc else [np.nan, lu["frac_aneuploid_lesion"]]
    anc_v = [acc["frac_anchor"], lu["frac_aneuploid_anchor"]] if acc else [np.nan, lu["frac_aneuploid_anchor"]]
    les_n = [acc["n_lesion"], lu["lesion_n_judged"]] if acc else [0, lu["lesion_n_judged"]]
    anc_n = [acc["n_anchor"], lu["anchor_n_judged"]] if acc else [0, lu["anchor_n_judged"]]
    x = np.arange(2); w = 0.36
    b1 = ax.bar(x - w / 2, les_v, w, color="#c0392b", label="lesion side (P11_LUAD, IAC)")
    b2 = ax.bar(x + w / 2, anc_v, w, color="#2e86c1", label="anchor side (P11_Normal, known normal)")
    for bars, vals, ns in ((b1, les_v, les_n), (b2, anc_v, anc_n)):
        for bar, v, n in zip(bars, vals, ns):
            if np.isnan(v):
                continue
            ax.text(bar.get_x() + bar.get_width() / 2, v + 0.012,
                    f"{v*100:.1f}%\nn={n}", ha="center", fontsize=8.6)
    ax.set_xticks(x); ax.set_xticklabels(conds, fontsize=9)
    ax.set_ylabel("Fraction called 'aneuploid'", fontsize=10)
    ax.set_ylim(0, 0.58)
    ax.axhline(0.10, color="#e67e22", ls="--", lw=1.2)
    ax.text(1.46, 0.115, "G1 flag 10%", fontsize=7.6, color="#b9770e", ha="right")
    ax.axhline(0.50, color="#7d3c98", ls=":", lw=1.2)
    ax.text(1.46, 0.515, "G1 stop 50%", fontsize=7.6, color="#7d3c98", ha="right")
    ax.set_title("(A) Same 728 nuclei: adding the anchor lifts the\n"
                 "lesion 9.8% -> 30.6%, anchor stays at 0%",
                 fontsize=11, loc="left")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.22)
    ax.text(0.025, 0.60,
            "the OFF bar is an UNPLANNED control\nleft by an implementation bug -- it is\n"
            "NOT a pre-registered arm",
            transform=ax.transAxes, fontsize=7.4, va="top",
            bbox=dict(fc="#fdf2e9", ec="#e67e22", alpha=0.95))

    # ---- (B) 两条预注册运行分侧 ----
    ax2 = fig.add_subplot(gs[0, 1])
    order = ["P11_LUAD", "P11_AAH"]
    xs = np.arange(2); w2 = 0.36
    lv = [runs[s]["frac_aneuploid_lesion"] for s in order]
    av = [runs[s]["frac_aneuploid_anchor"] for s in order]
    ln = [runs[s]["lesion_n_judged"] for s in order]
    an = [runs[s]["anchor_n_judged"] for s in order]
    b1 = ax2.bar(xs - w2 / 2, lv, w2, color="#c0392b", label="lesion side")
    b2 = ax2.bar(xs + w2 / 2, av, w2, color="#2e86c1", label="anchor side")
    for bars, vals, ns in ((b1, lv, ln), (b2, av, an)):
        for bar, v, n in zip(bars, vals, ns):
            ax2.text(bar.get_x() + bar.get_width() / 2, v + 0.012,
                     f"{v*100:.1f}%\nn={n}", ha="center", fontsize=8.6)
    ax2.set_xticks(xs)
    ax2.set_xticklabels([f"{s}\n({'IAC' if 'LUAD' in s else 'AAH'} lesion)" for s in order], fontsize=9)
    ax2.set_ylabel("Fraction called 'aneuploid'", fontsize=10)
    ax2.set_ylim(0, 0.58)
    ax2.axhline(0.10, color="#e67e22", ls="--", lw=1.2)
    ax2.text(1.46, 0.115, "G1 flag 10%", fontsize=7.6, color="#b9770e", ha="right")
    ax2.axhline(0.50, color="#7d3c98", ls=":", lw=1.2)
    ax2.text(1.46, 0.515, "G1 stop 50%", fontsize=7.6, color="#7d3c98", ha="right")
    ax2.set_title("(B) The two pre-registered runs: direction is right\n"
                  "for IAC, but AAH returns exactly ZERO",
                  fontsize=11, loc="left")
    ax2.legend(fontsize=8, loc="upper left")
    ax2.grid(axis="y", alpha=0.22)
    ax2.annotate("G2's positive control\nrests on this ONE bar",
                 xy=(0 - w2 / 2, lv[0]), xytext=(0.16, 0.455),
                 fontsize=8, color="#7b241c", ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color="#7b241c", lw=1.2,
                                 connectionstyle="arc3,rad=0.18"))
    ax2.text(1.0, 0.30, "no signal:\nNOT evidence of\n'no malignancy'",
             ha="center", fontsize=8.2, color="#1a5276",
             bbox=dict(fc="#eaf2f8", ec="#2e86c1", alpha=0.95))

    # ---- (C) 细胞数级联 ----
    ax3 = fig.add_subplot(gs[0, 2])
    # 每一行必须取**它自己那个样本**的 JSON；曾把 P11_LUAD 的数套到 AAH 行上（已修）
    rows = [
        ("P11_Normal\n(anchor)", runs["P11_LUAD"],
         ["anchor_ref_n_h5ad", "anchor_ref_n_after_floor",
          "anchor_ref_n_in_signed", "anchor_ref_n_of_subtype"], "#2e86c1"),
        ("P11_LUAD\n(lesion)", runs["P11_LUAD"],
         ["n_cells_h5ad", "n_cells_ge_floor",
          "anchor_les_n_in_signed", "anchor_les_n_of_subtype"], "#c0392b"),
        ("P11_AAH\n(lesion)", runs["P11_AAH"],
         ["n_cells_h5ad", "n_cells_ge_floor",
          "anchor_les_n_in_signed", "anchor_les_n_of_subtype"], "#8e44ad"),
    ]
    stage_names = ["all nuclei in sample", "pass 840 floor",
                   "in signed 128,091 set", "are AT2 (GP8a L2)"]
    gmax = max(max(runs[s][k] for k in ("n_cells_h5ad",)) for s in runs)
    gap, bh = 4.7, 0.78
    for ri, (label, src, keys, color) in enumerate(rows):
        vals = [src[k] for k in keys]
        y0 = ri * gap
        for si, v in enumerate(vals):
            # si=3 (AT2, the smallest) 画在最上；si=0 (全样本) 在最下
            ax3.barh(y0 + si, v, color=color, alpha=0.30 + 0.22 * si, height=bh)
            ax3.text(v + gmax * 0.015, y0 + si, f"{v:,}", va="center", fontsize=8.2)
        ax3.text(-gmax * 0.03, y0 + 1.5, label, ha="right", va="center",
                 fontsize=8.8, color=color, fontweight="bold")
        ax3.text(gmax * 0.40, y0 + 3.0,
                 f"kept past 840 floor: {100*vals[1]/vals[0]:.0f}%", fontsize=7.8,
                 color=color, va="center")
    ax3.set_ylim(-0.9, (len(rows) - 1) * gap + 4.2)
    ax3.set_yticks([])
    ax3.set_xlabel("nuclei (log-free linear scale)", fontsize=10)
    ax3.set_xlim(0, gmax * 1.14)
    ax3.set_title("(C) What each narrowing step costs\n"
                  "the 840 floor is the single biggest cut, in every sample",
                  fontsize=11, loc="left")
    ax3.grid(axis="x", alpha=0.22)
    handles = [plt.Rectangle((0, 0), 1, 1, fc="#555", alpha=0.30 + 0.22 * i)
               for i in range(4)]
    ax3.legend(handles, stage_names, fontsize=7.6, loc="lower right")

    # ---- (D) 诚实边界 ----
    axd = fig.add_subplot(gs[1, :]); axd.axis("off")
    axd.set_title("(D) What this run does and does NOT establish", fontsize=11.5, loc="left")
    txt = (
        "ESTABLISHED (for P11 only):\n"
        "  - The anchor branch was actually taken (copykat logged 'baseline is from known input').\n"
        "  - The anchor cells are clean AT2: EPCAM 0.408 / PTPRC 0.000, and 0/157 of them were\n"
        "    called aneuploid in BOTH runs  ->  G1 passes, the baseline is believable.\n"
        "  - For the IAC lesion the direction is finally right: lesion 30.6% > anchor 0.0%.\n"
        "  - Counts are exact: 728 judged = 175 aneuploid + 509 diploid + 44 not.defined.\n\n"
        "NOT ESTABLISHED:\n"
        "  - That CopyKAT is usable. ANCHOR_PREREG.md 9 says P11 is a smoke test; G2's positive\n"
        "    control rests entirely on the ONE non-zero bar in panel (B).\n"
        "  - That the AAH lesion has no malignant cells. 0/428 is a NULL, not a negative result:\n"
        "    absence of a call is not evidence of absence of CNAs.\n"
        "  - Any malignancy call whatsoever. No tumour cell is named here.\n\n"
        "STILL OPEN:\n"
        "  - copykat logged '152 known normal cells found in dataset' while we passed 157: it uses\n"
        "    the SMOOTHED matrix and drops 5-6 anchor cells from the baseline median. The baseline\n"
        "    is therefore built from 152 (not 157) cells. Registered, not corrected for.\n"
        "  - SCEVAN second arm: not installed, deferred by the user -- the user's 2026-09-17 signed\n"
        "    route requires BOTH arms to agree before this caliber is released.\n"
        "  - Whether the 840 floor removes true-CNA AAH nuclei: retention differs across samples\n"
        "    (Normal 35.2% / LUAD 47.6% / AAH 44.7%) and the floor drops LOW-COVERAGE nuclei.\n\n"
        "REPRODUCIBILITY (extra run, not in the pre-registration):\n"
        "  - Identical input re-run reproduces EXACTLY: same 728 verdicts cell-by-cell, and the\n"
        "    by-cell CNA matrix is byte-identical (sha256 a28e2259... both runs). copykat() sets\n"
        "    set.seed(1234) unconditionally and has no RNG in segmentation.\n"
        "  - So the difference in panel (A) is NOT noise. It is the anchor, and only the anchor."
    )
    axd.text(0.0, 0.98, txt, transform=axd.transAxes, ha="left", va="top",
             fontsize=8.3, family="monospace",
             bbox=dict(fc="#fbfbfb", ec="#999", alpha=0.95))

    out = os.path.join(FIG, "anchor_arm_p11.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")

    inv = {
        "figure": os.path.relpath(out, ROOT),
        "figure_sha256": sha256(out),
        "question": "显式锚定 + 单亚型输入同时满足时，CopyKAT 的判决方向对不对、稳不稳？",
        "answer": "方向对了（IAC 病灶 30.6% > 锚定 0.0%）；但样本量只够冒烟，不足以判方法可用。",
        "not_a_malignancy_call": True,
        "preregistration": "results/03_cnv/ANCHOR_PREREG.md",
        "per_run": {s: {k: runs[s][k] for k in (
            "n_cells_used", "n_judged", "n_pred_aneuploid", "n_pred_diploid",
            "n_pred_not_defined", "anchor_n_judged", "lesion_n_judged",
            "anchor_n_aneuploid", "lesion_n_aneuploid", "frac_aneuploid_anchor",
            "frac_aneuploid_lesion", "G1_verdict", "norm_ref_epcam_detect",
            "norm_ref_ptprc_detect", "norm_ref_n", "norm_ref_n_usable",
            "anchor_ref_n_of_subtype", "anchor_les_n_of_subtype",
            "param_coverage_floor")} for s in runs},
        "unplanned_control": acc,
        "reproducibility_check": {
            "outroot": "results/03_cnv/anchor_at2_repeat",
            "verdict": "reproducible",
            "cell_by_cell_disagreements": 0,
            "cna_matrix_sha256_identical": True,
            "cna_matrix_sha256": "a28e2259ea4a79f56aa0274be08e9ebb1ad01e87360e1d61655df878e7bd07c5",
            "note": "不在预注册里，是额加的可复现性核实；同输入重跑逐细胞判定与 CNA 矩阵字节全同。",
        },
        "caveats": [
            "🔴 首跑那次不是登记过的臂：实现 bug 使锚定在 2c 被清零，copykat 因此走了"
            "自动猜基线。它只作**意外对照**，不得当预注册臂引用。",
            "🔴 `frac_aneuploid_of_judged` 这条 JSON 字段在锚定臂里是**混算**的"
            "（分母含锚定细胞），单独引用会误导；必须用 frac_aneuploid_anchor / _lesion 分侧值。",
            "两次运行的锚定侧结果逐档相同（151 diploid / 6 not.defined / 0 aneuploid）"
            "⇒ 自洽，但不等于证明了全局可复现（需同输入重跑，见 anchor_at2_repeat）。",
            "P11 一对是冒烟与设计验证，不是全量结果，不得据单点判『方法可用/不可用』。",
            "不产生任何恶性判定，不命名肿瘤细胞。",
        ],
        "inputs_sha256": {
            os.path.relpath(os.path.join(OUT, f"{s}.json"), ROOT): sha256(os.path.join(OUT, f"{s}.json"))
            for s in runs},
    }
    mf = os.path.join(CNV, "anchor_arm_figure_manifest.json")
    with open(mf, "w") as fh:
        json.dump(inv, fh, ensure_ascii=False, indent=1)
    print("图:", out)
    print("清单:", mf)
    for s in runs:
        r = runs[s]
        print(f"  {s}: used={r['n_cells_used']} anchor={r['anchor_n_judged']}"
              f"({r['anchor_n_aneuploid']} aneu) lesion={r['lesion_n_judged']}"
              f"({r['lesion_n_aneuploid']} aneu) G1={r['G1_verdict']}")


if __name__ == "__main__":
    main()
