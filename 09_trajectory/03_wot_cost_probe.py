#!/usr/bin/env python
"""
WOT 成本探针 —— 只量「Waddington-OT 在我们的数据上要跑多久、吃多少内存」。

为什么这么设计：
    直接抽样 2000/期去外推不可靠（OT 的代价随 (n0×n1) 涨，且四对的规模悬殊）。
    所以**按比例抽样**：保持真实的期别大小比例，做 4 档越来越大的子样本，
    拟合出 scaling 规律，再外推到全量。

⚠️ 本次**不给任何生物学结论**：
    - 期别编码成数字 0..4 只是**排序手段**（WOT 用 sorted(set(day)) 定顺序；
      文字会按字母排成 AAH<AIS<IAC<MIA<Normal，顺序全错）。数值差不进数学。
    - 「把跨患者期别当成一条运输轴」这个口径**未签字**；本次只量代价，不解读耦合。
    - 归一化是为探针随手选的标准做法（CPM 1e4 + log1p），**正式跑的归一化口径未定**。

安全：外部 watchdog 按真 PID 监控，RSS 超阈值直接 SIGKILL。
"""
import json
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd

H5AD = "/home/eto/luad_v2/results/02_expression/gse308103_counts_paperqc.h5ad"
EPI_BC = "/home/eto/luad_v2/results/05_annotation/epiA_nocontam_subset_barcodes.txt"
OUT = "/home/eto/luad_v2/results/09_trajectory/smoke"
TMAP = os.path.join(OUT, "wot_probe_tmaps")

STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "IAC"]
# 按比例抽样的档位（分母）。
# 第一轮跑了 [64,32,16,8]：局部指数在大样本端是 0.94/0.88/0.67/0.52，明显往 1.0 涨，
# 说明小档被固定开销主导 ⇒ 用 0.466 的全局拟合外推会**低估**。
# 第二轮只跑大档 [16,8,4]，把拟合压在指数最要紧的区间。
FRACTIONS = [16, 8, 4]
SEED = 42


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


def load_epithelial():
    import anndata as ad
    import scanpy as sc

    log(f"读入 {os.path.basename(H5AD)} …")
    a = ad.read_h5ad(H5AD)
    log(f"  全库 {a.n_obs} × {a.n_vars}")

    want = {ln.strip() for ln in open(EPI_BC) if ln.strip()}
    keep = a.obs_names.isin(want)
    a = a[keep].copy()
    log(f"  上皮子集 {a.n_obs} × {a.n_vars}（目标 {len(want)}）")

    a.obs["day"] = pd.Categorical(
        a.obs["stage"].astype(str),
        categories=STAGE_ORDER, ordered=True)
    miss = int(a.obs["day"].isna().sum())
    if miss:
        raise ValueError(f"{miss} 个细胞期别不在 {STAGE_ORDER} 里")

    # 探针用的标准归一化（正式跑的口径未定）
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def probe_one(adata, denom, seed=SEED):
    """按比例抽一份子样本：每期取 max(60, n/denom) 个细胞。"""
    rng = np.random.default_rng(seed)
    idx = []
    for st in STAGE_ORDER:
        pool = np.where((adata.obs["stage"].astype(str) == st).values)[0]
        k = max(60, len(pool) // denom)
        k = min(k, len(pool))
        idx.append(rng.choice(pool, size=k, replace=False))
    idx = np.concatenate(idx)
    return idx


def run_probe(adata, denom, results):
    import wot

    idx = probe_one(adata, denom)
    sub = adata[idx].copy()
    # WOT 用 sorted(set(obs[day_field])) 定期别顺序 ⇒ 数字编码，按生物学顺序
    sub.obs["day"] = pd.Categorical(
        [STAGE_ORDER.index(s) for s in sub.obs["stage"].astype(str)],
        categories=list(range(len(STAGE_ORDER))), ordered=True)

    sizes = {s: int((sub.obs["stage"].astype(str) == s).sum()) for s in STAGE_ORDER}
    log(f"--- 档位 1/{denom} 共 {sub.n_obs} 细胞 {sizes} ---")

    for t0, t1 in zip(range(len(STAGE_ORDER) - 1), range(1, len(STAGE_ORDER))):
        m = (sub.obs["day"].astype(int).isin([t0, t1])).values
        pair = sub[m].copy()
        n0 = int((pair.obs["day"].astype(int) == t0).sum())
        n1 = int((pair.obs["day"].astype(int) == t1).sum())
        pre_rss = rss_gb()
        t_start = time.time()
        status, err = "ok", None
        try:
            om = wot.ot.OTModel(pair, day_field="day")
            om.compute_all_transport_maps(
                tmap_out=os.path.join(TMAP, f"f{denom}_{t0}_{t1}"), overwrite=True)
        except Exception as exc:            # noqa: BLE001
            status, err = "ERROR", f"{type(exc).__name__}: {str(exc)[:300]}"
            log(f"  !! {STAGE_ORDER[t0]}→{STAGE_ORDER[t1]} 失败: {err}")
        dt = round(time.time() - t_start, 1)
        rec = {"denom": denom, "pair": f"{STAGE_ORDER[t0]}->{STAGE_ORDER[t1]}",
               "n0": n0, "n1": n1, "n0xn1": n0 * n1, "wall_sec": dt,
               "rss_before_gb": pre_rss, "rss_after_gb": rss_gb(), "peak_rss_gb": hwm_gb(),
               "status": status, "error": err}
        results.append(rec)
        log(f"  {rec['pair']:16s} {n0:>5}×{n1:<5} = {n0*n1:>10}  {dt:>7.1f}s  "
            f"peak={hwm_gb()}GB  {status}")


def extrapolate(results):
    """按 log(time) ~ log(n0*n1) 拟合，外推到真实的四对规模。纯经验，误差未知。"""
    ok = [r for r in results if r["status"] == "ok" and r["wall_sec"] > 0
          and r["n0xn1"] > 0]
    real = {"Normal->AAH": 32251 * 11891, "AAH->AIS": 11891 * 28525,
            "AIS->MIA": 28525 * 4583, "MIA->IAC": 4583 * 56134}
    out = {"fit_ok": False, "caveat": "经验外推，OT 代价未必是幂律；只作量级参考"}
    if len(ok) < 4:
        out["reason"] = f"可用点只有 {len(ok)} 个，不够拟合"
        return out
    x = np.log(np.array([r["n0xn1"] for r in ok], dtype=float))
    y = np.log(np.array([r["wall_sec"] for r in ok], dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    out.update({
        "fit_ok": True,
        "n_points": len(ok),
        "slope_exponent": round(float(slope), 3),
        "r2": round(float(1 - resid.var() / y.var()), 4) if y.var() > 0 else None,
        "predicted_sec": {}, "predicted_total_sec": None,
    })
    tot = 0.0
    for name, prod in real.items():
        p = float(np.exp(slope * np.log(prod) + intercept))
        out["predicted_sec"][name] = round(p, 1)
        tot += p
    out["predicted_total_sec"] = round(tot, 1)
    out["predicted_total_human"] = f"{tot/60:.1f} 分钟" if tot < 7200 else f"{tot/3600:.1f} 小时"
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(TMAP, exist_ok=True)
    t0 = time.time()
    errors = []
    results = []
    try:
        adata = load_epithelial()
        for denom in FRACTIONS:
            run_probe(adata, denom, results)
    except Exception as exc:                # noqa: BLE001
        errors.append({"step": "main", "type": type(exc).__name__, "msg": str(exc)[:400],
                       "trace_tail": traceback.format_exc().strip().splitlines()[-3:]})
        log(f"!! 整体失败: {exc}")

    ex = extrapolate(results)
    manifest = {
        "kind": "wot_cost_probe",
        "scope": "NO BIOLOGICAL READOUT — 只量代价；方向源口径未签字；耦合不做任何解读",
        "input": {"h5ad": H5AD, "epithelial_barcodes": EPI_BC,
                  "note": "X 是原始计数；探针用 CPM 1e4 + log1p 归一（正式跑口径未定）"},
        "day_encoding": "Normal/AAH/AIS/MIA/IAC → 0..4，仅排序用；数值差不进数学",
        "wot_defaults": {"local_pca": 30, "growth_iters": 1, "epsilon": 0.05,
                         "lambda1": 1, "lambda2": 50, "tau": 10000,
                         "scaling_iter": 3000, "inner_iter_max": 50},
        "fractions": [f"1/{f}" for f in FRACTIONS],
        "results": results,
        "extrapolation": ex,
        "errors": errors,
        "total_wall_sec": round(time.time() - t0, 1),
        "peak_rss_gb": hwm_gb(),
        "contamination": "同时段后台在跑 P4 inferCNV（R，~300% CPU，38GB），墙钟受污染",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(OUT, "wot_cost_probe_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)
    log(f"总 {manifest['total_wall_sec']}s 峰值 {manifest['peak_rss_gb']}GB "
        f"{len(errors)} 个致命错误")
    log(f"外推: {ex}")


if __name__ == "__main__":
    sys.exit(main())
