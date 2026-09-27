#!/usr/bin/env python
"""
独立核验四个 WOT 耦合文件 —— 刻意**不复用** 06 里那段验收代码。

为什么必须独立：
    06 的验收报了 nan=True，我怀疑是我自己把布尔值的键名写反了（键叫 X_has_nan，
    装的其实是「无 NaN」），检查里又取了一次反 ⇒ 恒为 False。这正是本项目反复
    登记的那类「不报错但数错」。所以这里从头另写，不共享任何一行判据代码。

只核验，不产生生物学结论。
"""
import json
import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

IN = "/home/eto/luad_v2/results/09_trajectory/wot_full/epiA_wot_input.h5ad"
TDIR = "/home/eto/luad_v2/results/09_trajectory/wot_full/tmaps"
OUT = "/home/eto/luad_v2/results/09_trajectory/wot_full/independent_verify.json"

PAIRS = [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 4.0)]


def rss_gb():
    with open("/proc/self/status") as fh:
        for ln in fh:
            if ln.startswith("VmRSS:"):
                return round(int(ln.split()[1]) / 1024 / 1024, 2)
    return None


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] rss={rss_gb()}GB {m}", flush=True)


def main():
    log("读入输入 adata（拿 barcode 与期别做对照）…")
    a = ad.read_h5ad(IN)
    names = pd.Index(a.obs_names)
    stage_of = dict(zip(a.obs_names, a.obs["stage"].astype(str)))
    d = a.obs["day"].astype(float).values
    ncount = a.obs["nCount"].values.astype(float)
    log(f"  {a.n_obs} × {a.n_vars}；day 取值 {sorted(set(d))}")

    rep = {"input": {"n_obs": int(a.n_obs), "n_vars": int(a.n_vars),
                     "days": [float(x) for x in sorted(set(d))]}}

    # --- 深度事实：单核 nCount 的逐期 / 逐患者分布 ---
    st = a.obs["stage"].astype(str).values
    pat = a.obs["patient_id"].astype(str).values
    by_stage = {s: {"n": int((st == s).sum()),
                    "median": float(np.median(ncount[st == s])),
                    "p25": float(np.percentile(ncount[st == s], 25)),
                    "p75": float(np.percentile(ncount[st == s], 75))}
                for s in ["Normal", "AAH", "AIS", "MIA", "IAC"]}
    prec = ncount[np.isin(st, ["Normal", "AAH", "AIS", "MIA"])]
    luad = ncount[st == "IAC"]
    df = pd.DataFrame({"pat": pat, "st": st, "nc": ncount})
    per_pat = (df.groupby("pat")["nc"].median().round(0).astype(int).to_dict())
    rep["depth_single_nucleus"] = {
        "by_stage": by_stage,
        "precursor_median": float(np.median(prec)),
        "LUAD_median": float(np.median(luad)),
        "ratio_LUAD_over_precursor": round(float(np.median(luad) / np.median(prec)), 3),
        "patient_median_nCount": per_pat,
        "note": "这是单核 nCount；RCTD 线登记的 2.22 倍是**空转 spot** 深度，不是同一量纲。",
    }
    log(f"  单核深度：前驱中位 {np.median(prec):.0f} vs LUAD {np.median(luad):.0f} "
        f"= {np.median(luad)/np.median(prec):.3f} 倍（预注册里写的 2.22 是空转 spot 的数）")

    # --- 逐个耦合文件独立核验 ---
    per_pair = []
    for s0, s1 in PAIRS:
        f = os.path.join(TDIR, f"tmaps_{s0}_{s1}.h5ad")
        r = {"pair": f"{s0}->{s1}", "file": os.path.basename(f)}
        t = ad.read_h5ad(f)
        X = t.X
        r["shape"] = list(t.shape)
        r["X_type"] = "sparse" if sp.issparse(X) else "dense"
        r["X_dtype"] = str(X.dtype)

        # NaN：稀疏查 .data，密集查全部 —— 用 issaparse 分流，不用 hasattr(.data)
        if sp.issparse(X):
            vals = np.asarray(X.data, dtype=np.float64)
        else:
            vals = np.asarray(X, dtype=np.float64).ravel()
        r["n_stored"] = int(vals.size)
        r["n_nan"] = int(np.isnan(vals).sum())
        r["n_nonfinite"] = int((~np.isfinite(vals)).sum())
        r["min"] = float(np.nanmin(vals))
        r["max"] = float(np.nanmax(vals))

        # 名字：obs_names 必须是源期、var_names 必须是目标期、且都在输入里
        obs_in = pd.Index(t.obs_names).isin(names)
        var_in = pd.Index(t.var_names).isin(names)
        r["obs_all_in_input"] = bool(obs_in.all())
        r["var_all_in_input"] = bool(var_in.all())
        r["obs_all_src_stage"] = bool((pd.Series(t.obs_names).map(stage_of) ==
                                       {0.0: "Normal", 1.0: "AAH", 2.0: "AIS",
                                        3.0: "MIA", 4.0: "IAC"}[s0]).all())
        r["var_all_tgt_stage"] = bool((pd.Series(t.var_names).map(stage_of) ==
                                       {0.0: "Normal", 1.0: "AAH", 2.0: "AIS",
                                        3.0: "MIA", 4.0: "IAC"}[s1]).all())
        r["obs_unique"] = bool(pd.Index(t.obs_names).is_unique)
        r["var_unique"] = bool(pd.Index(t.var_names).is_unique)

        # 边际：行和应≈每个源细胞搬运出去的质量；列和≈每个目标细胞接收的质量
        rs = np.asarray(X.sum(axis=1)).ravel()
        cs = np.asarray(X.sum(axis=0)).ravel()
        r["rowsum"] = {"median": float(np.median(rs)), "min": float(rs.min()),
                       "max": float(rs.max()), "n_zero_rows": int((rs == 0).sum())}
        r["colsum"] = {"median": float(np.median(cs)), "min": float(cs.min()),
                       "max": float(cs.max()), "n_zero_cols": int((cs == 0).sum())}
        r["X_total"] = float(rs.sum())
        r["obs_cols"] = [c for c in t.obs.columns]
        r["growth_g0"] = {"median": float(np.median(t.obs["g0"])) if "g0" in t.obs else None}
        del X, rs, cs, t
        per_pair.append(r)
        log(f"  {r['pair']:10s} {r['X_type']:6s} NaN={r['n_nan']} 名字全对="
            f"{r['obs_all_in_input'] and r['var_all_in_input'] and r['obs_all_src_stage'] and r['var_all_tgt_stage']} "
            f"零行={r['rowsum']['n_zero_rows']} 零列={r['colsum']['n_zero_cols']}")

    rep["pairs"] = per_pair
    rep["_independent"] = ("刻意不复用 06 的验收代码：06 的 X_has_nan 键名与语义相反，"
                           "导致 all_checks_pass 恒为 False。这里从头另写。")
    with open(OUT, "w") as fh:
        json.dump(rep, fh, ensure_ascii=False, indent=1, default=str)
    log(f"写到 {OUT}")

    # 一句结论
    bad = [r["pair"] for r in per_pair
           if r["n_nan"] or r["n_nonfinite"] or not (r["obs_all_in_input"] and r["var_all_in_input"]
           and r["obs_all_src_stage"] and r["var_all_tgt_stage"] and r["obs_unique"] and r["var_unique"])]
    print(f"\n结论：{'全部干净，四个耦合可用' if not bad else '有问题: ' + str(bad)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
