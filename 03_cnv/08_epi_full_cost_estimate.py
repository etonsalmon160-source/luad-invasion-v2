#!/usr/bin/env python3
"""「全量上皮」CopyKAT 排期的**成本估算**（只估时间，不跑）。

要回答的问题：把 CopyKAT 从 P11 的 AT2 单点扩到**全部 75 个样本的全上皮**，
大概要多久？

口径（全部来自已入库产物，不重算）：
  1. 上皮细胞集合 = 用户签字的 `epiCNV_subset_barcodes.txt`（128,091 核，两标准交集）。
  2. 覆盖度地板 840 照旧（用户签字，本估算**不改**）。地板作用于 M1 的 `nFeature`。
  3. 亚型标签 = GP8a 上皮 L2（`clusters.csv.gz` 的 `harmony_res0.7_seed0`
     经 `epiA_cluster_annotation.csv` 的 `argmax` 映射）；挂不上的核剔除、不填补。
  4. 成本模型：由**已入库运行**的实测 (n, 秒) 最小二乘拟合 t = a + b·n^c。

本脚本给**两种排期**，因为用户签字的口径要求锚定"同患者 Normal **同一亚型**"：
  A. **每样本一跑**（全上皮混在一起）—— 与 `ab_epi_*` 臂同形，但**违反"同一亚型"**。
  B. **每（样本 × 亚型）一跑**，锚定取同患者 Normal 的**同一亚型** —— 忠于签字口径。

⚠️ 本脚本产出**只是排期预测**，不是任何生物学结论，也不构成 GP2 产物。
"""
import collections
import glob
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CNV = os.path.join(ROOT, "results", "03_cnv")
EPI = os.path.join(ROOT, "results", "05_annotation", "epiCNV_subset_barcodes.txt")
QC = os.path.join(ROOT, "results", "01_qc", "gse308103_per_cell_qc.csv.gz")
CLU = os.path.join(ROOT, "results", "04_integration", "seurat_trad", "epiA",
                   "clusters.csv.gz")
ANN = os.path.join(ROOT, "results", "05_annotation", "epiA_cluster_annotation.csv")

FLOOR = 840          # 用户签字，不可动
CORES = 20           # 实测 cpuset 上限（GP1 §4.3）
RAM_BUDGET_GB = 200.0
MEM_A, MEM_B = 0.00928, 0.798      # GP1 §1.1 三点拟合：peak ≈ a·n^b
MIN_ANCHOR = 10                    # 预注册 G5：锚定 ≥ copykat 自身 min.cells
SEED_COL = "harmony_res0.7_seed0"  # 现行 epiA r*=0.7
UNIFIED_ARMS = ("gp2", "ab_floor840_unified05", "anchor_at2", "anchor_at2_repeat")


# ---------------------------------------------------------------- 成本模型
def anchors():
    out = []
    for p in glob.glob(os.path.join(CNV, "**", "*.json"), recursive=True):
        if os.path.relpath(p, CNV).split(os.sep)[0] not in UNIFIED_ARMS:
            continue
        try:
            d = json.load(open(p))
        except Exception:
            continue
        if not isinstance(d, dict):
            continue
        n = d.get("n_cells_used") or d.get("n_cells_in") or d.get("n_cells_ge_floor")
        t = d.get("copykat_sec") or d.get("wall_sec_copykat")
        if n and t:
            out.append((int(n), float(t), os.path.relpath(p, ROOT)))
    return sorted(out)


def fit(pts):
    n = np.array([p[0] for p in pts], float)
    t = np.array([p[1] for p in pts], float)
    best = None
    for c in np.arange(1.0, 3.0, 0.005):
        X = np.column_stack([np.ones_like(n), n ** c])
        coef, *_ = np.linalg.lstsq(X, t, rcond=None)
        pred = X @ coef
        sse = float(((t - pred) ** 2).sum())
        if best is None or sse < best[0]:
            best = (sse, c, coef[0], coef[1], pred)
    sse, c, a, b, pred = best
    rel = np.abs(pred - t) / t
    return dict(a=a, b=b, c=float(c), n=len(pts),
                max_rel_err=float(rel.max()), med_rel_err=float(np.median(rel)))


def lpt(jobs, cores):
    load = [0.0] * cores
    for v in sorted(jobs, reverse=True):
        load[int(np.argmin(load))] += v
    return max(load) if load else 0.0


# ---------------------------------------------------------------- 输入构造
def load_inputs():
    epi = {l.strip() for l in open(EPI) if l.strip()}
    qc = pd.read_csv(QC, usecols=["sample_id", "cell_barcode", "nFeature"])
    df = qc[qc["cell_barcode"].isin(epi)].copy()
    n_signed = len(df)
    df = df[df["nFeature"] >= FLOOR]
    n_floor = len(df)

    clu = pd.read_csv(CLU, usecols=["cell_barcode", SEED_COL])
    ann = pd.read_csv(ANN, usecols=["cluster", "argmax"])
    m = dict(zip(ann["cluster"], ann["argmax"]))
    df = df.merge(clu, on="cell_barcode", how="left")
    df["subtype"] = df[SEED_COL].map(m)
    n_unmapped = int(df["subtype"].isna().sum())
    df = df[df["subtype"].notna()]
    return df, dict(n_signed=n_signed, n_floor=n_floor, n_unmapped=n_unmapped,
                    n_usable=len(df))


def patient_normals(df):
    """每患者的 Normal 样本。歧义（>1）**不自行消解**，原样报出。"""
    per = collections.defaultdict(set)
    for s in df["sample_id"].unique():
        per[s.split("_")[0]].add(s)
    out = {}
    for p, ss in per.items():
        nrm = sorted(s for s in ss if "Normal" in s)
        out[p] = nrm
    return out


def main():
    pts = anchors()
    m = fit(pts)
    sec = lambda n: m["a"] + m["b"] * (n ** m["c"])
    gb = lambda n: MEM_A * (n ** MEM_B)

    df, meta = load_inputs()
    idx = df.groupby(["sample_id", "subtype"]).size()
    normals = patient_normals(df)
    ambiguous = {p: v for p, v in normals.items() if len(v) != 1}

    # ---------- 排期 A：每样本一跑（全上皮混一起）----------
    per_s = df.groupby("sample_id").size()
    # ---------- 排期 B：每（样本 × 亚型）一跑，锚定同患者 Normal 同亚型 ----------
    jobs_b, skipped_no_anchor = [], []
    for (s, st), nv in idx.items():
        p = s.split("_")[0]
        nv = int(nv)
        if "Normal" in s:
            continue                      # Normal 作锚定源，不作病灶跑
        # P4 有两个 Normal ⇒ 两种选择都算，取影响报出
        an = max((int(idx.get((nn, st), 0)) for nn in normals.get(p, [])), default=0)
        if an < MIN_ANCHOR:
            skipped_no_anchor.append((s, st, nv, an))
            continue
        jobs_b.append((f"{s}|{st}", nv + an))

    rep = {
        "model": {"t_sec": f"{m['a']:.1f} + {m['b']:.6g} * n^{m['c']:.3f}",
                  "n_anchors": m["n"], "max_rel_err": round(m["max_rel_err"], 3),
                  "med_rel_err": round(m["med_rel_err"], 3)},
        "input": meta,
        "n_samples": int(per_s.size), "n_patients": len(normals),
        "ambiguous_patients": {p: v for p, v in ambiguous.items()},
        "A_per_sample_all_epi": {
            "n_runs": int(per_s.size),
            "cpu_h": sum(sec(int(v)) for v in per_s.values) / 3600,
            "wall_h": lpt([sec(int(v)) for v in per_s.values], CORES) / 3600,
            "max_single_h": sec(int(per_s.max())) / 3600,
            "max_peak_rss_gb": max(gb(int(v)) for v in per_s.values),
        },
        "B_per_sample_subtype": {
            "n_runs": len(jobs_b),
            "cpu_h": sum(sec(v) for _, v in jobs_b) / 3600,
            "wall_h": lpt([sec(v) for _, v in jobs_b], CORES) / 3600,
            "max_single_h": (sec(max(v for _, v in jobs_b)) / 3600) if jobs_b else 0,
            "max_peak_rss_gb": (max(gb(v) for _, v in jobs_b)) if jobs_b else 0,
            "skipped_anchor_too_small": len(skipped_no_anchor),
        },
    }
    print(json.dumps(rep, ensure_ascii=False, indent=1))

    print("\n锚点（同口径臂）:")
    for n, t, p in pts:
        print(f"  n={n:6d}  实测 {t:8.1f}s  模型 {sec(n):8.1f}s  "
              f"误差 {100*(sec(n)-t)/t:+6.1f}%   {p}")

    print("\n排期 A 最长 5 个样本:")
    rows = sorted(((s, int(v), sec(int(v)) / 3600) for s, v in per_s.items()),
                  key=lambda r: -r[2])[:5]
    for s, nv, h in rows:
        print(f"    {s:<14}{nv:>7} 核   {h:>6.2f} h")
    print("\n排期 B 最长 5 个（样本|亚型）:")
    rows = sorted(((k, v, sec(v) / 3600) for k, v in jobs_b), key=lambda r: -r[2])[:5]
    for k, nv, h in rows:
        print(f"    {k:<24}{nv:>7} 核   {h:>6.2f} h")
    if skipped_no_anchor:
        print(f"\n  ⚠️ B 排期有 {len(skipped_no_anchor)} 个 (样本×亚型) 因同患者 Normal "
              f"同亚型 <{MIN_ANCHOR} 核而未排（预注册 G5）：")
        for s, st, nv, an in sorted(skipped_no_anchor)[:12]:
            print(f"      {s:<14}{st:<16} 病灶 {nv:>5} 锚定 {an:>3}")

    out = os.path.join(CNV, "epi_full_cost_estimate.json")
    with open(out, "w") as fh:
        json.dump({"report": rep,
                   "anchors": [{"n": n, "sec": t, "source": p} for n, t, p in pts],
                   "per_sample_subtype": {f"{s}|{st}": int(v)
                                          for (s, st), v in idx.items()}},
                  fh, ensure_ascii=False, indent=1)
    print("\n清单:", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
