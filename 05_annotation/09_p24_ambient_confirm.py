#!/usr/bin/env python3
"""P24 环境 RNA —— **确认版**（按已冻结口径 `p24_ambient_caliber.json` 执行）。

与探索版（08）的两处设计差别，就是为了让"看过度量"这件事不至于让结论空心化：
  ① 分析单位从**细胞**改成**病人**（细胞不独立，细胞级 p 值会虚高）；
  ② 加一组**逐基因检出率匹配的负对照基因集**（200 组 × 10 基因，seed=0），
     排除"P24 的 nFeature 本来就高 ⇒ 什么程序都高"这个伪影。

🔴 脚本开头校验口径文件存在且哈希入册；缺则 SystemExit（法则 3.1）。

产物：results/05_annotation/p24_ambient_confirm_{per_patient.csv,null_distribution.csv,manifest.json}
      figures/p24_ambient_confirm.png
"""
import hashlib
import json
import os

import numpy as np
import pandas as pd

ROOT = "/home/eto/luad_v2"
R = os.path.join(ROOT, "results", "05_annotation")
FIG = os.path.join(ROOT, "figures")
H5 = os.path.join(ROOT, "results", "02_expression", "gse308103_counts_paperqc.h5ad")
CALIBER = os.path.join(R, "p24_ambient_caliber.json")
CELLS_CSV = os.path.join(R, "p24_ambient_ig_cells.csv.gz")
CHUNK = 20000
RNG_SEED = 0


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    if not os.path.exists(CALIBER):
        raise SystemExit("口径文件不存在 —— 口径未定不许计算（法则 3.1）")
    cal = json.load(open(CALIBER))
    cal_sha = sha256(CALIBER)
    print(f"口径 {cal['caliber']}  sha256={cal_sha[:16]}…  frozen_at={cal['frozen_at']}")

    PL = cal["gene_sets"]["plasma_program"]
    IG = cal["gene_sets"]["ig_constant"]
    NONB = cal["data"]["non_b_lineages"]

    cells = pd.read_csv(CELLS_CSV)
    cells["is_nonb"] = cells["A_readopted"].isin(NONB)
    nb = cells[cells["is_nonb"]].reset_index(drop=True)
    patients = sorted(cells["patient_id"].unique())
    print(f"非 B/浆 细胞 {len(nb)}  病人 {len(patients)}")

    # ---- 一次读入：非 B/浆 细胞的**逐基因检出**布尔矩阵（对 18,069 基因） ----
    import anndata
    ad = anndata.read_h5ad(H5, backed="r")
    var = list(ad.var_names)
    vset = set(var)
    for nm, gs in (("plasma_program", PL), ("ig_constant", IG)):
        miss = [g for g in gs if g not in vset]
        if miss:
            raise SystemExit(f"{nm} 缺基因 {miss} —— 停（口径与矩阵不符）")
    # 非 B 细胞的全局行号
    row_of = {bc: i for i, bc in enumerate(ad.obs_names)}
    nb_rows = np.array([row_of[bc] for bc in nb["cell_barcode"]], dtype=np.int64)
    del row_of
    keep = np.zeros(ad.n_obs, dtype=bool)
    keep[nb_rows] = True
    n_genes = ad.n_vars
    D = np.zeros((len(nb), n_genes), dtype=bool)
    pos = 0
    for s in range(0, ad.n_obs, CHUNK):
        e = min(s + CHUNK, ad.n_obs)
        k = keep[s:e]
        if not k.any():
            continue
        ck = ad.X[s:e]
        ck = ck.toarray() if not hasattr(ck, "toarray") else np.asarray(ck.todense())
        blk = ck > 0
        m = int(k.sum())
        D[pos:pos + m] = blk[k]
        pos += m
        if (s // CHUNK) % 5 == 0:
            print(f"  读入 {e}/{ad.n_obs}", flush=True)
    del ad
    assert pos == len(nb)
    print(f"检出矩阵 {D.shape}  ≈ {D.nbytes/1e9:.2f} GB")

    rate = D.mean(0)
    gi = {g: var.index(g) for g in var}
    pl_idx = [gi[g] for g in PL]
    ig_idx = [gi[g] for g in IG]

    # ---- 主检验：逐病人 浆程序 ≥3/10 ----
    pl_hit = (D[:, pl_idx].sum(1) >= 3).astype(np.float64)
    pid = nb["patient_id"].values
    pp = (pd.DataFrame({"patient_id": pid, "hit": pl_hit, "ig_umi_frac": nb["ig_umi_frac"].values})
            .groupby("patient_id")
            .agg(n_cells=("hit", "size"), plasma_ge3=("hit", "mean"),
                 ig_umi_frac_median=("ig_umi_frac", "median"))
            .reset_index())
    others = pp.loc[pp["patient_id"] != "P24", "plasma_ge3"]
    p24v = float(pp.loc[pp["patient_id"] == "P24", "plasma_ge3"].iloc[0])
    med = float(others.median())
    mad = float(np.median(np.abs(others - med))) * 1.4826
    rank = int((pp["plasma_ge3"] > p24v).sum()) + 1
    primary = dict(
        p24_value=p24v, others_median=med, others_max=float(others.max()),
        ratio_to_median=p24v / med if med else float("nan"),
        mad_scale=mad, robust_z=(p24v - med) / mad if mad else float("nan"),
        rank_among_patients=rank, n_patients=len(pp),
        rule=cal["primary_test"]["support_H1_if"],
        verdict_supports_H1=bool(rank == 1 and (p24v / med if med else 0) >= 3.0
                                 and mad and (p24v - med) / mad >= 3.0),
    )

    # ---- 特异度负对照：200 组逐基因检出率匹配的伪基因集 ----
    panel = set(PL) | set(IG) | {"MS4A1", "CD79A", "CD79B", "CD19", "CD37", "CD74",
                                 "BANK1", "FCRL1", "CD22", "TCL1A"}
    excl = np.array([(v in panel) or v.startswith(("MT-", "RPL", "RPS", "MTRNR"))
                     for v in var])
    pool = np.where((~excl) & (rate >= 0.01))[0]
    print(f"负对照候选池 {len(pool)} 个基因")
    tol = 0.05
    rng = np.random.RandomState(RNG_SEED)
    n_null = 200
    null_vals = np.full(n_null, np.nan)
    is_p24 = pid == "P24"
    for b in range(n_null):
        pick = []
        for g in pl_idx:
            cand = pool[np.abs(rate[pool] - rate[g]) <= tol]
            cand = cand[~np.isin(cand, pick)]
            if len(cand) == 0:
                pick.append(pool[rng.randint(len(pool))])
            else:
                pick.append(int(cand[rng.randint(len(cand))]))
        pick = np.array(pick)
        null_vals[b] = float(((D[is_p24][:, pick].sum(1) >= 3)).mean())
    p_emp = float((np.sum(null_vals >= p24v) + 1) / (n_null + 1))
    q95 = float(np.percentile(null_vals, 95))
    specificity = dict(
        null=f"{n_null} 组 × {len(PL)} 基因，逐基因检出率匹配 ±{tol}，seed={RNG_SEED}",
        pool_size=int(len(pool)), tolerance=tol, n_null=n_null,
        null_mean=float(np.mean(null_vals)), null_p95=q95, null_max=float(np.max(null_vals)),
        p24_on_plasma_program=p24v, empirical_p=p_emp,
        verdict_supports_H1=bool(p24v >= q95 and p_emp <= 0.05),
    )

    # ---- 次检验：IG 恒定区 UMI 占比（不同源读数） ----
    sec = pp.sort_values("ig_umi_frac_median", ascending=False).reset_index(drop=True)
    sec_rank = int(sec.index[sec["patient_id"] == "P24"][0]) + 1
    secondary = dict(
        readout=cal["secondary_test"]["readout"],
        p24_value=float(pp.loc[pp["patient_id"] == "P24", "ig_umi_frac_median"].iloc[0]),
        others_median=float(pp.loc[pp["patient_id"] != "P24", "ig_umi_frac_median"].median()),
        rank_among_patients=sec_rank, n_patients=len(pp),
        verdict_supports_H1=bool(sec_rank == 1),
    )

    # ---- 落盘 ----
    pp.sort_values("plasma_ge3", ascending=False).to_csv(
        os.path.join(R, "p24_ambient_confirm_per_patient.csv"), index=False)
    pd.DataFrame({"null_index": np.arange(n_null), "p24_plasma_ge3": null_vals}).to_csv(
        os.path.join(R, "p24_ambient_confirm_null_distribution.csv"), index=False)

    overall = bool(primary["verdict_supports_H1"] and specificity["verdict_supports_H1"]
                   and secondary["verdict_supports_H1"])
    out = dict(
        caliber=cal["caliber"], caliber_sha256=cal_sha, caliber_frozen_at=cal["frozen_at"],
        honesty_boundary=cal["honesty_boundary"],
        primary_patient_level=primary, specificity_negative_control=specificity,
        secondary_ig_umi_fraction=secondary,
        overall_verdict=("确认支持 H1（环境 RNA）" if overall else "未通过确认"),
        interpretation=("P24 的浆程序信号在**病人级**是全局最高的离群点，且显著高过逐基因匹配的"
                        "负对照 ⇒ 不是『P24 什么基因都高』的伪影；IG 恒定区占比同样最高（不同源读数）。"
                        if overall else "三项里至少一项未过，结论不成立，须如实降级。"),
        not_independent_validation=cal["honesty_boundary"],
        inputs={"caliber": cal_sha, "cells_csv": sha256(CELLS_CSV), "h5ad": sha256(H5)},
        outputs={"per_patient": "results/05_annotation/p24_ambient_confirm_per_patient.csv",
                 "null_distribution": "results/05_annotation/p24_ambient_confirm_null_distribution.csv",
                 "figure": "figures/p24_ambient_confirm.png"},
    )

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.8))

    ax = axes[0]
    o = pp.sort_values("plasma_ge3")
    cols = ["#c0392b" if p == "P24" else "#95a5a6" for p in o["patient_id"]]
    ax.barh(o["patient_id"], o["plasma_ge3"], color=cols)
    ax.axvline(3 * med, color="#e67e22", ls="--", lw=1.4, label="3x median of others")
    ax.axvline(med, color="#7f8c8d", ls=":", lw=1.2, label="median of others")
    ax.set_xlabel("Fraction of a patient's NON-B cells with >=3 plasma-program genes")
    ax.set_title(f"(A) PRIMARY - patient level (n={len(pp)})\n"
                 f"P24 rank {rank}/{len(pp)}, {primary['ratio_to_median']:.1f}x median, "
                 f"robust z={primary['robust_z']:.1f}\n"
                 f"-> H1 {'supported' if primary['verdict_supports_H1'] else 'NOT supported'}", fontsize=10.5)
    ax.legend(fontsize=8); ax.grid(axis="x", alpha=0.25); ax.tick_params(labelsize=8)

    ax = axes[1]
    ax.hist(null_vals, bins=28, color="#95a5a6", edgecolor="white")
    ax.axvline(p24v, color="#c0392b", lw=2, label=f"P24 on plasma program ({p24v:.3f})")
    ax.axvline(q95, color="#e67e22", ls="--", lw=1.4, label="null 95th pct")
    ax.set_xlabel("P24's value using a MATCHED RANDOM 10-gene set")
    ax.set_ylabel(f"count (of {n_null} null sets)")
    ax.set_title(f"(B) SPECIFICITY - matched negative control\n"
                 f"empirical p = {p_emp:.3f} (need <=0.05)\n"
                 f"-> H1 {'supported' if specificity['verdict_supports_H1'] else 'NOT supported'}", fontsize=10.5)
    ax.legend(fontsize=8)

    ax = axes[2]
    o2 = pp.sort_values("ig_umi_frac_median")
    cols = ["#c0392b" if p == "P24" else "#95a5a6" for p in o2["patient_id"]]
    ax.barh(o2["patient_id"], o2["ig_umi_frac_median"] * 100, color=cols)
    ax.set_xlabel("Median IG constant-region UMI fraction of a patient's NON-B cells (%)")
    ax.set_title(f"(C) SECONDARY - independent readout\n"
                 f"P24 rank {sec_rank}/{len(pp)}\n"
                 f"-> H1 {'supported' if secondary['verdict_supports_H1'] else 'NOT supported'}", fontsize=10.5)
    ax.grid(axis="x", alpha=0.25); ax.tick_params(labelsize=8)

    fig.tight_layout()
    fout = os.path.join(FIG, "p24_ambient_confirm.png")
    fig.savefig(fout, dpi=160, bbox_inches="tight")
    out["figure_sha256"] = sha256(fout)
    json.dump(out, open(os.path.join(R, "p24_ambient_confirm_manifest.json"), "w"),
              ensure_ascii=False, indent=1)

    print("\n===== 主检验（病人级）=====")
    for k, v in primary.items():
        if k != "rule":
            print(f"  {k}: {v}")
    print("\n===== 特异度负对照 =====")
    for k, v in specificity.items():
        print(f"  {k}: {v}")
    print("\n===== 次检验（IG 恒定区占比）=====")
    for k, v in secondary.items():
        print(f"  {k}: {v}")
    print(f"\n总判定：{out['overall_verdict']}")
    print("图:", fout)


if __name__ == "__main__":
    main()
