#!/usr/bin/env python3
"""CopyKAT 各臂对照图：判决不是由 CNV 决定的（只读已入库 JSON，不重算）。

要回答的问题：CopyKAT 在 GSE308103 上能不能用来判"这个上皮细胞是不是肿瘤"？

图回答的是三件事：
  (A) 全细胞臂：**方向是反的** —— Normal 样本判出的"非整倍体"比 MIA 病灶还高。
  (B) 同一个样本，只换喂进去的细胞集，非整倍体率在 0%–41% 之间翻转 ⇒ 判决由输入决定。
  (C) 三条通路各自的假阳性/假阴性（供裁决走哪条路）。

数据来源：results/03_cnv/**/*.json（9 月 16 日各臂的逐样本产物，均已带 sha256）。
无 CJK 字体 ⇒ 全英文标签。
"""
import glob
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

ARM_LABEL = {
    "gp2":                      "whole-cell, floor 840, UP.DR 0.10",
    "ab_floor840_unified05":    "whole-cell, floor 840, UP.DR 0.05",
    "ab_floor0":                "whole-cell, no floor",
    "ab_epi_argmax":            "epithelium only (argmax)",
    "ab_epi_strict":            "epithelium only (strict)",
    "ab_epi_at2":               "epithelium, AT2 only",
    "trad_anchor_floor840":     "whole-cell, ANCHORED (non-epi ref)",
    "trad_default":             "whole-cell, ANCHORED, no floor",
}
WHOLE_ARMS = ["gp2", "ab_floor840_unified05", "ab_floor0"]
STAGE_COLOR = {"Normal": "#2e86c1", "AAH": "#8e44ad", "AIS": "#16a085",
               "MIA": "#e67e22", "LUAD": "#c0392b"}
# 逐样本上色（三个 Normal 样本若同色就分不开）；实线=Normal，虚线=病灶
SAMPLE_COLOR = {"P3_Normal": "#5dade2", "P11_Normal": "#1a5276",
                "P13_Normal": "#a569bd", "P14_Normal": "#154360",
                "P13_MIA": "#d35400"}
SAMPLE_MARKER = {"P3_Normal": "o", "P11_Normal": "s", "P13_Normal": "^",
                 "P14_Normal": "D", "P13_MIA": "v"}
SAMPLE_STYLE = {"P3_Normal": "-", "P11_Normal": "-", "P13_Normal": "-",
                "P14_Normal": "-", "P13_MIA": "--"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load():
    rows = []
    for f in sorted(glob.glob(os.path.join(CNV, "**", "*.json"), recursive=True)):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        if "frac_aneuploid_of_judged" not in d:
            continue
        arm = os.path.relpath(os.path.dirname(f), CNV)
        rows.append(dict(arm=arm, f=f, sample=d["sample_id"], stage=d["stage_token"],
                         n=d["n_cells_used"], aneu=d["n_pred_aneuploid"],
                         frac=d["frac_aneuploid_of_judged"],
                         epi=d.get("epi_rule"), sub=d.get("subtype"),
                         anchor=d.get("norm_ref_mode")))
    return rows


def main():
    os.makedirs(FIG, exist_ok=True)
    rows = load()

    fig = plt.figure(figsize=(19.0, 7.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 0.95], wspace=0.30)

    # ---- (A) 全细胞臂：方向是反的 ----
    ax = fig.add_subplot(gs[0, 0])
    sub = [r for r in rows if r["arm"] in WHOLE_ARMS]
    samples = ["P3_Normal", "P11_Normal", "P13_Normal", "P14_Normal", "P13_MIA"]
    x = np.arange(len(samples)); w = 0.26
    for i, arm in enumerate(WHOLE_ARMS):
        vals = []
        for s in samples:
            m = [r for r in sub if r["arm"] == arm and r["sample"] == s]
            vals.append(m[0]["frac"] if m else np.nan)
        ax.bar(x + (i - 1) * w, vals, width=w, label=ARM_LABEL[arm])
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", "\n") for s in samples], fontsize=9)
    ax.axhline(0, color="#333", lw=1)
    ax.set_ylabel("Fraction called 'aneuploid' by CopyKAT", fontsize=10)
    ax.set_ylim(0, 0.72)
    ax.set_title("(A) Whole-cell arm: the direction is INVERTED\n"
                 "3 Normal samples outrank the MIA lesion (which gets 0%)",
                 fontsize=11.5, loc="left")
    ax.legend(fontsize=7.8, loc="upper right")
    ax.grid(axis="y", alpha=0.22)
    ax.annotate("expected: ~0%", xy=(1.5, 0.015), fontsize=8, color="#1b4f72")
    ax.annotate("expected: HIGH", xy=(4.0, 0.02), fontsize=8, color="#a04000")

    # ---- (B) 同一样本跨臂翻转 ----
    ax2 = fig.add_subplot(gs[0, 1])
    show = ["P11_Normal", "P13_Normal", "P13_MIA", "P14_Normal"]
    order = ["gp2", "ab_floor840_unified05", "ab_floor0",
             "ab_epi_argmax", "ab_epi_strict", "ab_epi_at2",
             "trad_anchor_floor840", "trad_default"]
    xs = np.arange(len(order))
    for s in show:
        ys, ok = [], []
        for i, arm in enumerate(order):
            m = [r for r in rows if r["arm"] == arm and r["sample"] == s]
            ys.append(m[0]["frac"] if m else np.nan)
            ok.append(bool(m))
        st = [r["stage"] for r in rows if r["sample"] == s][0]
        ax2.plot(xs, ys, marker=SAMPLE_MARKER[s], ls=SAMPLE_STYLE[s],
                 color=SAMPLE_COLOR[s], label=f"{s} ({st})", lw=1.9, ms=8)
    ax2.set_xticks(xs)
    ax2.set_xticklabels([ARM_LABEL[a].replace(", ", ",\n") for a in order],
                        rotation=40, ha="right", fontsize=7.6)
    ax2.set_ylabel("Fraction called 'aneuploid'", fontsize=10)
    ax2.set_ylim(-0.03, 0.70)
    ax2.axhline(0, color="#333", lw=1)
    ax2.legend(fontsize=8.5, loc="upper left", ncol=2)
    ax2.grid(axis="y", alpha=0.22)
    ax2.set_title("(B) Same sample, same nuclei -- only the input cell set changes\n"
                  "P11_Normal flips 0-32%, P13_Normal 0-41%",
                  fontsize=11.5, loc="left")
    ax2.text(0.99, 0.72,
             "A verdict that only depends on\nwhich cells you feed is not a\nCNV measurement.",
             transform=ax2.transAxes, ha="right", va="top", fontsize=8,
             bbox=dict(fc="#fdf2e9", ec="#e67e22", alpha=0.95))

    # ---- (C) 三条通路各自的硬边界 ----
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.axis("off")
    ax3.set_title("(C) What is structurally blocked", fontsize=11.5, loc="left")
    txt = (
        "CopyKAT: verdict = f(input expression axis,\n"
        "guessed baseline, one binary branch).\n"
        "Baseline is GUESSED (norm.cell.names=\"\"\n"
        "-> baseline.norm.cl(); 10 of 13 arms logged\n"
        "'low confidence' then answered anyway).\n"
        "Narrowing the input RAISED the rate 0%->32%\n"
        "=> not a CNV estimate (GP1_report.md 9.7).\n\n"
        "Method availability on THIS dataset:\n"
        "  [NO] Numbat  -- needs BAM/allele. GSE308103\n"
        "              ships 75 count matrices only.\n"
        "  [NO] infercnvpy -- needs Python >=3.10; here 3.8.10\n"
        "  [??] inferCNV (R) -- feasible, needs a forced reference\n"
        "  [??] SCEVAN -- count-matrix-only; YOU deferred it\n"
        "              (waiting on the main-line architecture)\n\n"
        "The ONLY sourced way out (GP1 9.9):\n"
        "  ANCHOR (pass norm.cell.names explicitly)\n"
        "  + HOMOGENEOUS (one subtype in, malignant+normal)\n"
        "Not yet run in combination. The AT2-only arm\n"
        "satisfied 'homogeneous' but NOT 'anchor', and\n"
        "returned 0% on BOTH Normal and LUAD (n=124/160)\n"
        "=> no positive signal in a tumour sample.")
    ax3.text(0.0, 0.95, txt, transform=ax3.transAxes, ha="left", va="top",
             fontsize=8.0, family="monospace",
             bbox=dict(fc="#fbfbfb", ec="#999", alpha=0.95))

    out = os.path.join(FIG, "cnv_arm_instability.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")

    inv = {
        "figure": os.path.relpath(out, ROOT),
        "figure_sha256": sha256(out),
        "question": "CopyKAT 在 GSE308103 上能否判'这个上皮细胞是肿瘤'？",
        "answer": "当前配置下不能。判决由喂进去的细胞集决定，不由 CNV 决定。",
        "evidence": {
            "A_inverted_direction": "全细胞臂下 Normal 样本 26%/41%/59% 判非整倍体，而 MIA 病灶 0%。",
            "B_flip_within_sample": "同一个样本只换输入细胞集，非整倍体率在 0%–41% 之间翻转。",
            "C_blocked_methods": "Numbat 需 BAM（本数据集无）；infercnvpy 需 py>=3.10（本机 3.8.10）；"
                                 "CopyKAT 的唯一有出处出路是'锚定+同质'，该组合尚未跑。",
        },
        "arms_included": {a: ARM_LABEL[a] for a in ARM_LABEL},
        "per_arm_results": [{"arm": r["arm"], "sample": r["sample"], "stage": r["stage"],
                             "n_cells_used": r["n"], "n_aneuploid": r["aneu"],
                             "frac_aneuploid": r["frac"], "epi_rule": r["epi"],
                             "subtype": r["sub"], "anchor": r["anchor"]} for r in rows],
        "caveats": [
            "本图只汇总已入库的 2026-09-16 各臂 JSON，不重算任何东西。",
            "各臂的细胞集与参数不同是**故意**的（那是被测变量），所以横向比较的是'判决的稳定性'，"
            "不是'哪个臂更准'。",
            "各臂产物均为**归因实验**，按 EPI_PREREG.md §6 与 GP1 §9.9，不构成 M2 正式产物。",
            "不产生任何恶性判定。",
        ],
        "inputs_sha256": {os.path.relpath(r["f"], ROOT): sha256(r["f"]) for r in rows},
    }
    mf = os.path.join(CNV, "arm_flip_figure_manifest.json")
    json.dump(inv, open(mf, "w"), ensure_ascii=False, indent=1)
    print("图:", out)
    print("清单:", mf)
    for a in order:
        v = [r for r in rows if r["arm"] == a]
        print(f"  {a:26s} n_samples={len(v)}  fracs={[round(x['frac'],3) for x in v]}")


if __name__ == "__main__":
    main()
