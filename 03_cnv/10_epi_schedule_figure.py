#!/usr/bin/env python3
"""「全量上皮」CopyKAT 排期图：挂钟被**单个最长运行**决定 + 锚定覆盖缺口。

只读 `results/03_cnv/epi_full_cost_estimate.json`（由 `08_epi_full_cost_estimate.py` 产）。
无 CJK 字体 ⇒ 全英文标签。
"""
import collections
import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CNV = os.path.join(ROOT, "results", "03_cnv")
FIG = os.path.join(ROOT, "figures")
EST = os.path.join(CNV, "epi_full_cost_estimate.json")
QC = os.path.join(ROOT, "results", "01_qc", "gse308103_per_cell_qc.csv.gz")
EPI = os.path.join(ROOT, "results", "05_annotation", "epiCNV_subset_barcodes.txt")
CLU = os.path.join(ROOT, "results", "04_integration", "seurat_trad", "epiA",
                   "clusters.csv.gz")
ANN = os.path.join(ROOT, "results", "05_annotation", "epiA_cluster_annotation.csv")
FLOOR, MIN_ANCHOR = 840, 10
SEED_COL = "harmony_res0.7_seed0"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rebuild():
    epi = {l.strip() for l in open(EPI) if l.strip()}
    qc = pd.read_csv(QC, usecols=["sample_id", "cell_barcode", "nFeature"])
    df = qc[qc["cell_barcode"].isin(epi) & (qc["nFeature"] >= FLOOR)].copy()
    clu = pd.read_csv(CLU, usecols=["cell_barcode", SEED_COL])
    ann = pd.read_csv(ANN, usecols=["cluster", "argmax"])
    df = df.merge(clu, on="cell_barcode", how="left")
    df["st"] = df[SEED_COL].map(dict(zip(ann["cluster"], ann["argmax"])))
    return df[df["st"].notna()]


def main():
    est = json.load(open(EST))
    m = est["report"]["model"]
    a, b, c = (float(x) for x in
               (m["t_sec"].split(" + ")[0], m["t_sec"].split(" + ")[1].split(" * ")[0],
                m["t_sec"].split("^")[1]))
    sec = lambda n: a + b * (n ** c)

    df = rebuild()
    idx = df.groupby(["sample_id", "st"]).size()
    pat = lambda s: s.split("_")[0]
    norm = collections.defaultdict(list)
    for s in df["sample_id"].unique():
        norm[pat(s)].append(s)
    norm = {p: [x for x in v if "Normal" in x] for p, v in norm.items()}

    cov = collections.Counter(); bad = collections.Counter(); cells = collections.Counter()
    hours = []
    for (s, st), n in idx.items():
        if "Normal" in s:
            continue          # Normal 作锚定源，不作病灶跑；**不得**计入组合分母
        cov[st] += 1
        cells[st] += int(n)
        an = max((int(idx.get((nn, st), 0)) for nn in norm.get(pat(s), [])), default=0)
        if an < MIN_ANCHOR:
            bad[st] += 1
        else:
            hours.append((sec(int(n) + an) / 3600, f"{s}|{st}", int(n) + an))
    hours.sort(reverse=True)

    fig = plt.figure(figsize=(19.0, 6.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.30)

    # ---- (A) 201 个可跑组合，按时长排序 ----
    ax = fig.add_subplot(gs[0, 0])
    hv = [h for h, _, _ in hours]
    ax.bar(range(len(hv)), hv, color="#c0392b", width=0.85)
    ax.set_xlabel("the 201 runnable (sample x subtype) runs, longest first", fontsize=10)
    ax.set_ylabel("predicted wall time (h)", fontsize=10)
    ax.set_title("(A) Wall clock = ONE run, not the sum\n"
                 "P4_LUAD|AT2 alone is %.1f h of the %.1f h total" % (hv[0], sum(hv)),
                 fontsize=11, loc="left")
    ax.annotate(f"{hours[0][1]}\n{hours[0][2]:,} nuclei\n{hv[0]:.1f} h",
                xy=(0, hv[0]), xytext=(14, hv[0] * 0.86), fontsize=8.6,
                arrowprops=dict(arrowstyle="->", color="#7b241c", lw=1.2))
    ax.text(0.62, 0.55,
            "adding cores cannot help\nbelow the longest single run",
            transform=ax.transAxes, fontsize=8.4, color="#7b241c",
            bbox=dict(fc="#fdf2e9", ec="#e67e22", alpha=0.95))
    ax.grid(axis="y", alpha=0.22)

    # ---- (B) 锚定覆盖缺口 ----
    ax2 = fig.add_subplot(gs[0, 1])
    order = sorted(cov, key=lambda s: -cells[s])
    okv = [cov[s] - bad[s] for s in order]
    bdv = [bad[s] for s in order]
    y = np.arange(len(order))
    ax2.barh(y, okv, color="#2e86c1", label="runnable (anchor >= 10)")
    ax2.barh(y, bdv, left=okv, color="#c0392b",
             label="BLOCKED: paired Normal has <10 of this subtype")
    for i, s in enumerate(order):
        ax2.text(okv[i] + bdv[i] + 1.2, i, f"{okv[i]}/{cov[s]}", va="center", fontsize=8.4)
    ax2.set_yticks(y); ax2.set_yticklabels(order, fontsize=9)
    ax2.set_xlabel("number of (sample x subtype) combinations", fontsize=10)
    ax2.set_title("(B) Under the signed caliber the whole\n"
                  "epithelium CANNOT be run", fontsize=11, loc="left")
    ax2.legend(fontsize=7.8, loc="lower right")
    ax2.set_xlim(0, max(cov.values()) * 1.28)
    ax2.grid(axis="x", alpha=0.22)

    # ---- (C) 核数覆盖 ----
    ax3 = fig.add_subplot(gs[0, 2])
    tot = sum(cells.values())
    run = sum(cells[s] for s in order if s not in ("Serous", "Basal", "Goblet/Mucous",
                                                   "Ciliated"))
    ax3.bar([0], [run], color="#2e86c1", label="AT2 + AT1: 100%% runnable")
    ax3.bar([1], [tot], color="#7f8c8d", alpha=0.45, label="all epithelium")
    for x, v, lab in ((0, run, f"{run:,}\n({100*run/tot:.1f}%)"),
                      (1, tot, f"{tot:,}\n(100%)")):
        ax3.text(x, v * 0.5, lab, ha="center", va="center", fontsize=10,
                 color="white", fontweight="bold")
    ax3.set_xticks([0, 1])
    ax3.set_xticklabels(["AT2 + AT1\n(102/102 runs)", "all 6 subtypes\n(201/280 runs)"],
                        fontsize=9)
    ax3.set_ylabel("epithelial nuclei (post 840 floor)", fontsize=10)
    ax3.set_title("(C) The gap is rare subtypes: 10.7% of nuclei\n"
                  "sit in partially- or un-anchorable subtypes",
                  fontsize=11, loc="left")
    ax3.legend(fontsize=8, loc="upper left")
    ax3.grid(axis="y", alpha=0.22)

    fig.text(0.005, 0.012,
             "Predicted only -- no CopyKAT run was performed for this figure. Cost model "
             "t = %.1f + %.6g*n^%.3f, fitted on %d archived runs (median error %.0f%%, max %.0f%%). "
             "Floor 840 and 'anchor >= 10' (prereg G5) unchanged. P4 has TWO Normal samples "
             "(P4_Normal, P4_Normal1) -- pairing not silently resolved."
             % (a, b, c, m["n_anchors"], 100 * m["med_rel_err"], 100 * m["max_rel_err"]),
             fontsize=7.4, family="monospace")

    os.makedirs(FIG, exist_ok=True)
    out = os.path.join(FIG, "epi_full_schedule.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")

    inv = {
        "figure": os.path.relpath(out, ROOT), "figure_sha256": sha256(out),
        "question": "把 CopyKAT 扩到全量上皮要多久？",
        "answer": "忠于签字口径（同亚型锚定）挂钟约 6.2 h、CPU 43 h，但 79/280 组合无合格锚定。",
        "single_longest_run_h": hours[0][0], "wall_harmonic_bound_h": hours[0][0],
        "sum_serial_h": sum(hv),
        "coverage": {s: {"combos": cov[s], "blocked": bad[s], "lesion_nuclei": cells[s]}
                     for s in order},
        "caveats": [
            "纯预测，未为本图跑任何 CopyKAT。",
            "模型在同 n 不同样本上有 ±43% 的散布（n=393 时 130.6s vs 186.9s）⇒ 挂钟应按区间读，约 5–8 h。",
            "P4 有两个 Normal 样本，配对歧义未消解（本图锚定取两者中核数较多者）。",
            "挂钟受单个最长运行限制，加核不可缩短。",
        ],
        "inputs_sha256": {os.path.relpath(p, ROOT): sha256(p) for p in (EST, EPI, CLU, ANN)},
    }
    mf = os.path.join(CNV, "epi_full_schedule_figure_manifest.json")
    with open(mf, "w") as fh:
        json.dump(inv, fh, ensure_ascii=False, indent=1)
    print("图:", out)
    print("清单:", mf)
    print(f"  最长单跑 {hours[0][1]} {hours[0][2]} 核 {hours[0][0]:.2f} h；"
          f"串行合计 {sum(hv):.1f} h；可跑 {len(hv)} 跑")


if __name__ == "__main__":
    main()
