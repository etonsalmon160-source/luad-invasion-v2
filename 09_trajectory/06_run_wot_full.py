#!/usr/bin/env python
"""
WOT 全量运输图 —— 正式跑（口径已签字，见 results/09_trajectory/WOT_PREREG.md）。

三条已签口径：
  a. 跨患者期别 Normal→AAH→AIS→MIA→IAC 当成一条运输轴        （用户「我签」2026-09-25）
  b. 归一化 = CP10K(1e4) + log1p                              （用户选定 2026-09-25）
  c. 跨患者深度差**不单独处理**，深度只作报告属性、不进筛选    （用户选定 2026-09-25）

本步只产出耦合矩阵，不给任何生物学结论。验收条款见预注册第 7 节。
"""
import json
import os
import time
import traceback

import numpy as np
import pandas as pd
import scipy.sparse as sp

H5AD = "/home/eto/luad_v2/results/02_expression/gse308103_counts_paperqc.h5ad"
EPI_BC = "/home/eto/luad_v2/results/05_annotation/epiA_nocontam_subset_barcodes.txt"
OUT = "/home/eto/luad_v2/results/09_trajectory/wot_full"
TMAP_DIR = os.path.join(OUT, "tmaps")
INPUT_H5AD = os.path.join(OUT, "epiA_wot_input.h5ad")

STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "IAC"]
DAY = {s: float(i) for i, s in enumerate(STAGE_ORDER)}
# 预注册里写定的期望形状（全项目统一上皮口径下的逐期计数）
EXPECT_COUNTS = {"Normal": 32251, "AAH": 11891, "AIS": 28525, "MIA": 4583, "IAC": 56134}


def rss_gb():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return round(int(line.split()[1]) / 1024 / 1024, 2)
    return None


def hwm_gb():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmHWM:"):
                return round(int(line.split()[1]) / 1024 / 1024, 2)
    return None


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] rss={rss_gb()}GB {m}", flush=True)


def build_input():
    import anndata as ad
    import scanpy as sc

    log(f"读入 {os.path.basename(H5AD)} …")
    a = ad.read_h5ad(H5AD)
    log(f"  全库 {a.n_obs} × {a.n_vars}")

    want = {ln.strip() for ln in open(EPI_BC) if ln.strip()}
    a = a[a.obs_names.isin(want)].copy()
    log(f"  上皮子集 {a.n_obs} × {a.n_vars}（目标 {len(want)}）")
    if a.n_obs != len(want):
        raise ValueError(f"上皮子集 {a.n_obs} != 清单 {len(want)}，口径对不上，停")

    # 期别 → float（官方 recipe：先 float 再 category）。数值只当排序标签，运输只在相邻期之间解。
    a.obs["day"] = (a.obs["stage"].astype(str).map(DAY)).astype(float).astype("category")
    if a.obs["day"].isna().any():
        bad = a.obs.loc[a.obs["day"].isna(), "stage"].value_counts().to_dict()
        raise ValueError(f"有细胞期别不在 {STAGE_ORDER} 里: {bad}")

    tp = sorted(set(a.obs["day"].astype(float)))
    log(f"  WOT 会看到的 timepoints = {tp}")
    if tp != [0.0, 1.0, 2.0, 3.0, 4.0]:
        raise ValueError(f"timepoints 排序不对: {tp} —— 运输轴顺序会错，停")

    counts = {s: int((a.obs["stage"].astype(str) == s).sum()) for s in STAGE_ORDER}
    log(f"  逐期计数 {counts}")
    for s in STAGE_ORDER:
        if counts[s] != EXPECT_COUNTS[s]:
            raise ValueError(f"期别 {s} 计数 {counts[s]} != 预注册 {EXPECT_COUNTS[s]}，停")

    # 口径 b：CP10K + log1p
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    log("  口径 b 已施加：CP10K(1e4) + log1p")

    a.write_h5ad(INPUT_H5AD)
    log(f"  预处理后的输入写到 {INPUT_H5AD}（供下游用同一对象，防口径漂移）")
    return a, counts


def depth_report(a):
    """口径 c 签的是「深度只作报告属性」⇒ 这里把跨患者/跨期深度诊断出全，供报告用。"""
    df = a.obs[["patient_id", "sample_id", "stage", "nCount"]].copy()
    df["nCount"] = df["nCount"].astype(float)
    # 全部转成原生 python 类型，免得 json 里出现 numpy 标量被 default=str 变成字符串
    by_stage = {str(s): {"count": int(g.size),
                         "median": round(float(g.median()), 1),
                         "mean": round(float(g.mean()), 1),
                         "std": round(float(g.std()), 1)}
                for s, g in df.groupby("stage")["nCount"]}
    by_pat = {str(p): {"count": int(g.size), "median": round(float(g.median()), 1)}
              for p, g in df.groupby("patient_id")["nCount"]}
    prec = df[df["stage"].isin(["Normal", "AAH", "AIS", "MIA"])]["nCount"].median()
    luad = df[df["stage"] == "IAC"]["nCount"].median()
    return {
        "note": "口径 c：深度只作报告属性，不进筛选",
        "by_stage": by_stage, "by_patient": by_pat,
        "median_precursor_vs_LUAD": {
            "precursor_median": float(prec), "LUAD_median": float(luad),
            "ratio_LUAD_over_precursor": round(float(luad) / float(prec), 3),
            "prereg_expected_ratio": 2.22,
            "reading": "RCTD 线登记的是 2.22 倍；这里核对能不能复现。对不上要记录，不改口径。",
        },
    }


def _nan_free(X):
    """稀疏查存储值；密集查全部（.ravel() 是视图，bool 临时量约 0.4 GB，可接受）。"""
    vals = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if vals.size == 0:
        return True
    return not bool(np.isnan(vals).any())


def verify_tmaps(a, counts):
    """预注册第 7 节的验收条款，逐条核。"""
    import anndata as ad

    ok_names = set(a.obs_names)
    stage_of = dict(zip(a.obs_names, a.obs["stage"].astype(str)))
    checks, records = [], []

    for i in range(len(STAGE_ORDER) - 1):
        s0, s1 = STAGE_ORDER[i], STAGE_ORDER[i + 1]
        f = os.path.join(TMAP_DIR, f"tmaps_{DAY[s0]}_{DAY[s1]}.h5ad")
        rec = {"pair": f"{s0}->{s1}", "file": f, "expected_shape":
               [EXPECT_COUNTS[s0], EXPECT_COUNTS[s1]], "exists": os.path.exists(f)}
        if rec["exists"]:
            t = ad.read_h5ad(f)
            X = t.X
            obs_src = pd.Series(t.obs_names)
            var_tgt = pd.Series(t.var_names)
            rec.update({
                "shape": list(t.shape),
                "shape_ok": list(t.shape) == rec["expected_shape"],
                "X_dtype": str(X.dtype),
                # 稀疏判据必须用 scipy.sparse.issparse：numpy 数组也有 .data（memoryview），
                # 用 hasattr(X,"data") 会把原始字节当数值查 NaN —— 静默错。
                "X_is_sparse": bool(sp.issparse(X)),
                "X_has_nan": _nan_free(X),
                "X_sum": round(float(X.sum()), 2),
                "n_obs_names_in_adata": int(obs_src.isin(ok_names).sum()),
                "n_var_names_in_adata": int(var_tgt.isin(ok_names).sum()),
                "n_obs_names_stage_ok": int((obs_src.map(stage_of) == s0).sum()),
                "n_var_names_stage_ok": int((var_tgt.map(stage_of) == s1).sum()),
                "obs_cols": [c for c in t.obs.columns],
            })
            rec["all_names_match"] = bool(
                rec["n_obs_names_in_adata"] == t.shape[0] and
                rec["n_var_names_in_adata"] == t.shape[1] and
                rec["n_obs_names_stage_ok"] == t.shape[0] and
                rec["n_var_names_stage_ok"] == t.shape[1])
            # 行和/列和：耦合的边际该与源/目标细胞数同量级
            rs = np.asarray(X.sum(axis=1)).ravel()
            rec["rowsum_median"] = round(float(np.median(rs)), 4)
            rec["total_mass"] = rec["X_sum"]
            rec["mass_vs_n_source"] = round(rec["X_sum"] / t.shape[0], 4)
        records.append(rec)
        checks.append(rec["exists"] and rec.get("shape_ok", False)
                      and rec.get("all_names_match", False)
                      and not rec.get("X_has_nan", True))

    # 目录清洁条款：from_wot 会 glob 全部 h5ad，多一个文件就会静默错数
    strays = [f for f in os.listdir(TMAP_DIR) if not f.endswith(".h5ad")]
    n_files = len([f for f in os.listdir(TMAP_DIR) if f.endswith(".h5ad")])
    return {"records": records, "n_h5ad_in_dir": n_files,
            "dir_clean_4_files": n_files == 4,
            "non_h5ad_strays": strays,
            "all_checks_pass": bool(all(checks)) and n_files == 4}


def main():
    os.makedirs(TMAP_DIR, exist_ok=True)
    t0 = time.time()
    errors, result = [], {}
    try:
        import wot
        a, counts = build_input()
        result["depth_report"] = depth_report(a)

        log(f"开跑 WOT（4 对相邻期别）… 输出目录 {TMAP_DIR}")
        om = wot.ot.OTModel(a, day_field="day")
        om.compute_all_transport_maps(tmap_out=TMAP_DIR + "/", overwrite=True)
        log("WOT 四对都算完了")

        result["verification"] = verify_tmaps(a, counts)
        log(f"验收: {result['verification']['all_checks_pass']}")
    except Exception as exc:                      # noqa: BLE001
        errors.append({"step": "main", "type": type(exc).__name__,
                       "msg": str(exc)[:600],
                       "trace_tail": traceback.format_exc().strip().splitlines()[-4:]})
        log(f"!! 失败: {exc}")

    manifest = {
        "kind": "wot_full_transport_maps",
        "scope": "只产出耦合矩阵，零生物学结论；三条口径见 WOT_PREREG.md",
        "signed_calibers": {
            "transport_axis": "跨患者 Normal->AAH->AIS->MIA->IAC 作为一条运输轴（用户「我签」2026-09-25）",
            "normalization": "CP10K(1e4) + log1p（用户选定 2026-09-25）",
            "patient_depth": "不单独处理；深度只作报告属性（用户选定 2026-09-25）",
        },
        "input": {"h5ad": H5AD, "epithelial_barcodes": EPI_BC,
                  "X": "原始计数 → 口径 b 归一化", "n_genes": 18069,
                  "gene_filter": "None（全基因，不按方差挑）"},
        "day_encoding": "Normal/AAH/AIS/MIA/IAC -> 0.0/1.0/2.0/3.0/4.0（float Categorical，"
                        "只作排序标签；运输只在相邻期之间解，数值差不进数学）",
        "wot_parameters": om.ot_config if not errors else None,
        "stage_counts": counts if not errors else None,
        "result": result,
        "errors": errors,
        "total_wall_sec": round(time.time() - t0, 1),
        "peak_rss_gb": hwm_gb(),
        "contamination": "同时段后台在跑 P4 inferCNV（R，~40GB RSS）⇒ 墙钟受污染",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(OUT, "wot_full_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False, default=str)
    log(f"总 {manifest['total_wall_sec']}s 峰值 {manifest['peak_rss_gb']}GB "
        f"{len(errors)} 个致命错误")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
