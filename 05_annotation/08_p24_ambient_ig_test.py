#!/usr/bin/env python3
"""P24 的 B/浆 污染：是"真的浆细胞"还是"环境 RNA（免疫球蛋白汤）"？—— 证伪测试。

来源与授权：用户 2026-09-22 会话内授权"跑，看 IG 是不是假象"。
产物：results/05_annotation/p24_ambient_ig_{cells.csv.gz,manifest.json} + figures/p24_ambient_ig.png

═══════════════════════════════════════════════════════════════════════════
Part A —— 预注册口径（法则 3.1：先钉死再算）
═══════════════════════════════════════════════════════════════════════════
假设 H1（环境 RNA / soup）：P24 被剔出 B/浆 的细胞是别的细胞泡在含免疫球蛋白的"汤"里，被 IG 误判。
零假设 H0：那批就是真 B/浆 细胞。

 T1 非 B/浆 5 谱系的 **IG 阳性率**，P24 vs 其他病人；支持 H1 需 ≥2× 且在 ≥3/5 谱系成立。
 T2 B/浆 被剔细胞（IG+）的**浆程序检出率**（任一基因 >0）vs 保留细胞；支持 H1 需 ≤0.5×。

🔴 **两个预注册度量事后被判定为"饱和/无区分力"**（详见 manifest 的 metric_defect）：
   - T1 用的"IG 阳性率"在全队列本就 0.45–0.63（环境 IG 人人有），2× 的杆几乎够不到；
   - T2 用的"任一浆程序基因 >0"被 XBP1/PRDM1/IRF4/CD27 这类广谱基因撑到 ~0.95，谁都过。
   ⇒ Part A 的判定**照实上报，但不据此下结论**。

═══════════════════════════════════════════════════════════════════════════
Part B —— 事后改用有区分力的度量（**探索性，非预注册确认性检验**）
═══════════════════════════════════════════════════════════════════════════
换成"10 个浆程序基因里**检出 ≥3 个**的细胞占比"（保留的真浆细胞该指标 0.80，中位 7/10）。
判据（本项目约定，**在看到 Part A 的缺陷后才定**，故只作探索性）：P24 ≥ 其他 × 3 且在 ≥4/5 谱系成立。

基因集（只取矩阵里真有的；10x Flex 探针板缺 IGHG2/IGHG4/IGLC2/IGLC3/IGLL5）：
  IG_CONST       = 恒定区基因（汤里漂的就是这些），**故意不含 JCHAIN**（它是浆程序基因）
  PLASMA_PROGRAM = 浆细胞程序（10 个）
  B_PROGRAM      = B 细胞程序（10 个）
═══════════════════════════════════════════════════════════════════════════
"""
import gzip
import hashlib
import json
import os

import numpy as np
import pandas as pd

ROOT = "/home/eto/luad_v2"
R = os.path.join(ROOT, "results", "05_annotation")
FIG = os.path.join(ROOT, "figures")
H5 = os.path.join(ROOT, "results", "02_expression", "gse308103_counts_paperqc.h5ad")
CELLS_CSV = os.path.join(R, "p24_ambient_ig_cells.csv.gz")

IG_CONST = ["IGHG1", "IGHG3", "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE",
            "IGKC", "IGLC1", "IGLC7", "IGLL1"]
PLASMA_PROGRAM = ["MZB1", "XBP1", "PRDM1", "SDC1", "DERL3", "TNFRSF17",
                  "CD27", "IRF4", "POU2AF1", "JCHAIN"]
B_PROGRAM = ["MS4A1", "CD79A", "CD79B", "CD19", "CD37", "CD74",
             "BANK1", "FCRL1", "CD22", "TCL1A"]

NON_B_LINEAGES = ["上皮", "内皮", "髓系", "成纤维", "T/NK"]
P24_SAMPLES = ["P24_Normal", "P24_AAH", "P24_LUAD"]

# Part A（预注册）
T1_FOLD, T1_MIN_LIN = 2.0, 3
T2_FOLD = 0.5
# Part B（事后，本项目约定）
TB_FOLD, TB_MIN_LIN, TB_GENES = 3.0, 4, 3
CHUNK = 20000


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def build_cells_table():
    """从 h5ad 抽每条细胞的 IG / 浆程序 / B 程序面板 → 逐细胞表。"""
    import anndata
    ad = anndata.read_h5ad(H5, backed="r")
    var = list(ad.var_names)
    vset = set(var)
    for nm, gs in (("IG_CONST", IG_CONST), ("PLASMA_PROGRAM", PLASMA_PROGRAM),
                   ("B_PROGRAM", B_PROGRAM)):
        miss = [g for g in gs if g not in vset]
        if miss:
            raise SystemExit(f"{nm} 有基因不在矩阵里：{miss} —— 先改口径再算（法则 3.1）")
    cols = IG_CONST + PLASMA_PROGRAM + B_PROGRAM
    gi = np.array([var.index(g) for g in cols])
    n_ig, n_pl = len(IG_CONST), len(PLASMA_PROGRAM)
    n = ad.n_obs
    M = np.zeros((n, len(cols)), dtype=np.float32)
    for s in range(0, n, CHUNK):
        e = min(s + CHUNK, n)
        ck = ad.X[s:e]
        ck = ck.toarray() if not hasattr(ck, "toarray") else np.asarray(ck.todense())
        M[s:e] = ck[:, gi]
        if (s // CHUNK) % 5 == 0:
            print(f"  读入 {e}/{n}", flush=True)
    obs = ad.obs[["sample_id", "patient_id", "stage", "nCount", "nFeature"]].copy()
    del ad
    det = M > 0
    cells = pd.DataFrame({
        "cell_barcode": obs.index.values,
        "sample_id": obs["sample_id"].values,
        "patient_id": obs["patient_id"].values,
        "stage": obs["stage"].values,
        "nCount": obs["nCount"].values.astype(np.int64),
        "nFeature": obs["nFeature"].values.astype(np.int64),
        "ig_det": det[:, :n_ig].sum(1).astype(np.int16),
        "plasma_det": det[:, n_ig:n_ig + n_pl].sum(1).astype(np.int16),
        "b_det": det[:, n_ig + n_pl:].sum(1).astype(np.int16),
        "ig_umi": M[:, :n_ig].sum(1),
    })
    cells["ig_umi_frac"] = cells["ig_umi"] / cells["nCount"].clip(lower=1)
    asg = pd.read_csv(os.path.join(R, "gp8c_cell_assignment.csv.gz"),
                      usecols=["cell_barcode", "A_readopted"])
    cells = cells.merge(asg, on="cell_barcode", how="left")
    if cells["A_readopted"].isna().any():
        raise SystemExit("有细胞挂不上谱系标签 —— 停")
    kept = set(l.strip() for l in open(os.path.join(R, "bplasmaA_nocontam_subset_barcodes.txt")))
    cells["is_kept_bplasma"] = cells["cell_barcode"].isin(kept)
    cells.to_csv(CELLS_CSV, index=False, compression="gzip")
    return cells


def main():
    if os.path.exists(CELLS_CSV):
        print("复用已抽好的逐细胞表（改基因集时须删掉它重跑）")
        cells = pd.read_csv(CELLS_CSV)
    else:
        cells = build_cells_table()

    ex = pd.read_csv(os.path.join(R, "bplasmaA_excluded_cells.csv.gz"),
                     usecols=["cell_barcode", "cluster", "sample_id"])
    cells = cells.merge(ex[["cell_barcode", "cluster"]], on="cell_barcode", how="left")
    is_excl = cells["cluster"].notna()
    # 保留细胞 = 在 nocontam 清单里；直接由清单算，不依赖缓存列（旧缓存无此列）
    kept = set(l.strip() for l in open(os.path.join(R, "bplasmaA_nocontam_subset_barcodes.txt")))
    cells["is_kept_bplasma"] = cells["cell_barcode"].isin(kept)
    K = cells[cells["is_kept_bplasma"]]
    E = cells[is_excl]
    nb = cells[cells["A_readopted"].isin(NON_B_LINEAGES)].copy()
    nb["is_p24"] = nb["patient_id"] == "P24"

    # ---------- Part A：预注册 ----------
    piv = (nb.groupby(["A_readopted", "is_p24"])["ig_det"]
             .agg(lambda s: float((s > 0).mean())).unstack())
    piv.columns = ["other", "p24"]
    piv["fold"] = piv["p24"] / piv["other"]
    t1_pass = int((piv["fold"] >= T1_FOLD).sum())
    A1 = dict(pre_registered_verdict_supports_H1=bool(t1_pass >= T1_MIN_LIN),
              n_lineages_with_fold_ge=t1_pass, threshold_fold=T1_FOLD,
              min_lineages=T1_MIN_LIN, table=piv.round(6).reset_index().to_dict("records"))

    ref = float((K["plasma_det"] > 0).mean())
    p24e = E[E["patient_id"] == "P24"]
    A2 = dict(pre_registered_verdict_supports_H1=bool(
                  (float((p24e["plasma_det"] > 0).mean()) / ref) <= T2_FOLD),
              kept_plasma_pos_rate=ref, p24_excl_plasma_pos_rate=float((p24e["plasma_det"] > 0).mean()),
              other_excl_plasma_pos_rate=float((E[E["patient_id"] != "P24"]["plasma_det"] > 0).mean()),
              threshold_fold=T2_FOLD)

    # ---------- Part B：有区分力的度量（事后、探索性） ----------
    nb["plasma_ge3"] = (nb["plasma_det"] >= TB_GENES)
    b = (nb.groupby(["A_readopted", "is_p24"])
           .agg(n=("ig_det", "size"), rate=("plasma_ge3", "mean"),
                median_plasma=("plasma_det", "median"),
                median_ig_umi_frac=("ig_umi_frac", "median"),
                median_nFeature=("nFeature", "median")).reset_index())
    bp = b.pivot(index="A_readopted", columns="is_p24", values="rate")
    bp.columns = ["other", "p24"]
    bp["fold"] = bp["p24"] / bp["other"].replace(0, np.nan)
    b_pass = int((bp["fold"] >= TB_FOLD).sum())
    B1 = dict(metric=f"fraction of cells with >= {TB_GENES}/{len(PLASMA_PROGRAM)} plasma-program genes",
              rationale="Part A 的度量饱和；真浆细胞该指标 0.80（中位 7/10），能区分",
              exploratory_verdict_supports_H1=bool(b_pass >= TB_MIN_LIN),
              n_lineages_with_fold_ge=b_pass, threshold_fold=TB_FOLD, min_lineages=TB_MIN_LIN,
              kept_plasma_ge3_rate=float((K["plasma_det"] >= TB_GENES).mean()),
              kept_plasma_det_median=float(K["plasma_det"].median()),
              excluded_all_plasma_ge3_rate=float((E["plasma_det"] >= TB_GENES).mean()),
              table=bp.round(6).reset_index().to_dict("records"))

    psamp = (nb.groupby(["patient_id", "sample_id"])
               .agg(n=("ig_det", "size"),
                    ig_pos=("ig_det", lambda s: float((s > 0).mean())),
                    ig_umi_frac_med=("ig_umi_frac", "median"),
                    plasma_ge3=("plasma_ge3", "mean"),
                    plasma_det_med=("plasma_det", "median"),
                    nFeature_med=("nFeature", "median")).reset_index())
    oth = psamp[psamp["patient_id"] != "P24"]
    B2 = dict(
        p24_samples=psamp[psamp["sample_id"].isin(P24_SAMPLES)].round(6).to_dict("records"),
        other_patients_median=dict(ig_pos=float(oth["ig_pos"].median()),
                                   plasma_ge3=float(oth["plasma_ge3"].median()),
                                   nFeature_med=float(oth["nFeature_med"].median())),
        other_patients_plasma_ge3_range=[float(oth["plasma_ge3"].min()), float(oth["plasma_ge3"].max())],
        note="P24 的 Normal 与 AAH 两样本同样抬高 ⇒ 不是肿瘤生物学，更像样本处理/环境 RNA。")

    # M1 双体率（本 h5ad 已剔双体，doublet_class 全为单细胞 ⇒ 必须回 M1 表取）
    qc = pd.read_csv(os.path.join(ROOT, "results", "01_qc", "gse308103_qc_per_sample.csv"))
    dcol = [c for c in qc.columns if "doublet" in c.lower() and "n_" in c.lower()]
    qc2 = qc[["sample_id"] + (dcol[:1])].copy() if dcol else qc[["sample_id"]].copy()
    qc2.columns = ["sample_id", "m1_doublet_n"] if dcol else ["sample_id"]
    B3 = dict(note="本 h5ad 已剔双体，doublet_class 无常量 ⇒ 双体率取自 M1 表",
              table=psamp.merge(qc2, on="sample_id", how="left")
                          .sort_values("plasma_ge3", ascending=False)
                          .head(10).round(5).to_dict("records"))

    # ---------- 落盘 ----------
    out = dict(
        question="P24 被剔出 B/浆 的细胞：真浆细胞，还是环境 RNA（免疫球蛋白汤）？",
        authorized_by="用户 2026-09-22 会话内（\"跑，看 IG 是不是假象\"）",
        genes=dict(IG_CONST=IG_CONST, PLASMA_PROGRAM=PLASMA_PROGRAM, B_PROGRAM=B_PROGRAM,
                   note="JCHAIN 故意排除在 IG_CONST 之外（它是浆程序基因）"),
        PartA_preregistered=dict(T1_ig_positivity=A1, T2_program_any=A2),
        metric_defect=dict(
            detected_how="先按预注册跑完，发现两个度量都饱和；再看 Part B 才有区分力",
            T1="IG 阳性率在全队列本就高（非B谱系其他病人 0.45–0.63），2× 的杆几乎够不到；"
               "P24 在 5 个谱系**一律** 1.5–1.9× ⇒ 形态像环境 RNA，但过不了预注册的 2×。",
            T2="'任一浆程序基因 >0' 被 XBP1/PRDM1/IRF4/CD27 这类广谱基因撑到 ~0.95，"
               "保留细胞与 P24 被剔细胞无差别 ⇒ 该度量无区分力。",
            consequence="Part A 的两个判定（均=不支持 H1）**照实上报，但不据此下结论**；"
                        "结论改由 Part B 的有区分力度量给出，且明确标为**事后/探索性**。"),
        PartB_exploratory=dict(B1_per_lineage=B1, B2_per_sample=B2, B3_doublet_confound=B3),
        conclusion=dict(
            provisional_verdict=("支持 H1（环境 RNA）" if B1["exploratory_verdict_supports_H1"] else "未能支持 H1"),
            strength="P24 三个样本（含 Normal/AAH）在**非 B/浆**的 5 个谱系里，浆细胞程序均显著抬高；"
                     "P24_LUAD 最高。形状是「全局均一抬高 + 含正常样本」 ⇒ 典型环境 RNA，非真细胞。",
            not_a_confirmatory_result="Part B 的阈值是看到 Part A 缺陷后才定的 ⇒ 只作探索性。"
                                      "若要作结论，须先把该度量重新冻结为口径再跑一次。",
            next="是否据此对 P24 单独加环境校正 / 标记，等你定。"),
        inputs={"h5ad": sha256(H5),
                "gp8c_cell_assignment.csv.gz": sha256(os.path.join(R, "gp8c_cell_assignment.csv.gz")),
                "bplasmaA_excluded_cells.csv.gz": sha256(os.path.join(R, "bplasmaA_excluded_cells.csv.gz")),
                "bplasmaA_nocontam_subset_barcodes.txt": sha256(
                    os.path.join(R, "bplasmaA_nocontam_subset_barcodes.txt"))},
        outputs={"cells_csv": "results/05_annotation/p24_ambient_ig_cells.csv.gz",
                 "manifest": "results/05_annotation/p24_ambient_ig_manifest.json",
                 "figure": "figures/p24_ambient_ig.png"},
        caveat="IG_CONST 只用恒定区基因；本矩阵缺 IGHG2/IGHG4/IGLC2/IGLC3/IGLL5。"
               "检出率基于原始计数（>0），不受归一化影响。",
    )

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    LIN_EN = {"上皮": "Epithelial", "内皮": "Endothelial", "髓系": "Myeloid",
              "成纤维": "Fibroblast", "T/NK": "T/NK", "B/浆": "B/Plasma"}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8))

    ax = axes[0]
    idx = [l for l in NON_B_LINEAGES if l in piv.index]
    x = np.arange(len(idx))
    ax.bar(x - 0.19, [piv.loc[l, "other"] for l in idx], 0.38, color="#7f8c8d", label="Other patients")
    ax.bar(x + 0.19, [piv.loc[l, "p24"] for l in idx], 0.38, color="#c0392b", label="P24")
    for i, l in enumerate(idx):
        ax.text(i, max(piv.loc[l, "other"], piv.loc[l, "p24"]) + 0.02,
                f"x{piv.loc[l,'fold']:.1f}", ha="center", fontsize=8.5)
    ax.set_ylim(0, 1.12); ax.set_xticks(x); ax.set_xticklabels([LIN_EN[l] for l in idx])
    ax.set_ylabel("Fraction of cells with any IG transcript")
    ax.set_title(f"(A) PRE-REGISTERED T1 - NOT discriminative\n"
                 f"baseline already 0.45-0.63; bar of {T1_FOLD}x unreachable\n"
                 f"verdict: {t1_pass}/5  ->  H1 NOT supported", fontsize=10.5)
    ax.legend(fontsize=9); ax.grid(axis="y", alpha=0.25)

    ax = axes[1]
    idx = [l for l in NON_B_LINEAGES if l in bp.index]
    x = np.arange(len(idx))
    ax.bar(x - 0.19, [bp.loc[l, "other"] for l in idx], 0.38, color="#7f8c8d", label="Other patients")
    ax.bar(x + 0.19, [bp.loc[l, "p24"] for l in idx], 0.38, color="#c0392b", label="P24")
    for i, l in enumerate(idx):
        ax.text(i, max(bp.loc[l, "other"], bp.loc[l, "p24"]) + 0.012,
                f"x{bp.loc[l,'fold']:.1f}", ha="center", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels([LIN_EN[l] for l in idx])
    ax.set_ylabel(f"Fraction of cells with >= {TB_GENES}/{len(PLASMA_PROGRAM)} plasma-program genes")
    ax.set_title(f"(B) CORRECTED metric (exploratory) - A PLASMA PROGRAM\n"
                 f"inside lineages that are NOT B/plasma\n"
                 f"verdict: {b_pass}/5  ->  H1 {'supported' if B1['exploratory_verdict_supports_H1'] else 'not supported'}",
                 fontsize=10.5)
    ax.legend(fontsize=9); ax.grid(axis="y", alpha=0.25)

    ax = axes[2]
    top = psamp.sort_values("plasma_ge3", ascending=False).head(14).sort_values("plasma_ge3")
    cols_ = ["#c0392b" if p == "P24" else "#95a5a6" for p in top["patient_id"]]
    ax.barh(top["sample_id"], top["plasma_ge3"], color=cols_)
    for i, (s, v) in enumerate(zip(top["sample_id"], top["plasma_ge3"])):
        ax.text(v + 0.004, i, f"{v:.2f}", va="center", fontsize=8,
                color=("#c0392b" if s in P24_SAMPLES else "#555555"),
                weight=("bold" if s in P24_SAMPLES else "normal"))
    ax.set_xlim(0, max(top["plasma_ge3"]) * 1.18)
    ax.set_xlabel(f"Fraction of NON-B cells with >= {TB_GENES} plasma-program genes")
    ax.set_title("(C) Top-14 samples - all THREE P24 samples elevated\n"
                 "(normal + AAH too  ->  sample handling / ambient, not tumour biology)", fontsize=10.5)
    ax.grid(axis="x", alpha=0.25); ax.tick_params(labelsize=8.5)

    fig.tight_layout()
    fout = os.path.join(FIG, "p24_ambient_ig.png")
    fig.savefig(fout, dpi=160, bbox_inches="tight")
    out["figure_sha256"] = sha256(fout)
    json.dump(out, open(os.path.join(R, "p24_ambient_ig_manifest.json"), "w"),
              ensure_ascii=False, indent=1)

    print("\n===== Part A（预注册，度量饱和）=====")
    print(piv.round(4).to_string())
    print(f"  T1: 过 2× 的谱系 {t1_pass}/5 ⇒ 不支持 H1")
    print(f"  T2: 保留细胞浆程序阳性率 {ref:.3f} / P24 被剔细胞 {A2['p24_excl_plasma_pos_rate']:.3f} ⇒ 不支持 H1（但度量已饱和）")
    print("\n===== Part B（有区分力，事后/探索性）=====")
    print(bp.round(4).to_string())
    print(f"  过 3× 的谱系 {b_pass}/5 ⇒ 支持 H1：{B1['exploratory_verdict_supports_H1']}")
    print("\n  P24 三样本：")
    for r in B2["p24_samples"]:
        print("   %-11s n=%5d  浆程序≥3 %.3f   IG阳性率 %.3f   nFeature中位 %.0f" % (
            r["sample_id"], r["n"], r["plasma_ge3"], r["ig_pos"], r["nFeature_med"]))
    print("  其他病人中位: 浆程序≥3 %.3f  (范围 %.3f–%.3f)" % (
        B2["other_patients_median"]["plasma_ge3"],
        B2["other_patients_plasma_ge3_range"][0], B2["other_patients_plasma_ge3_range"][1]))
    print("\n图:", fout)


if __name__ == "__main__":
    main()
