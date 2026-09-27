#!/usr/bin/env python
"""
CellRank 成本冒烟 —— 只量「公共主干」的墙钟与内存，不产出任何生物学结论。

为什么只量公共主干：
    用户批准的第三点（用病程做"实时"方向轴）需要 cr.kernels.RealTimeKernel，
    而它要求预先算好的最优传输耦合（from_moscot / from_wot / 显式 couplings /
    重写 compute_coupling），且本 venv（Python 3.8）里 wot/moscot/ott/pot/scot
    一个都没装。方向源口径未签字 ⇒ 本脚本**故意不选任何方向源**。

    因此这里只量所有核选择都要经过的那一段：kNN 图 → 转移矩阵 → GPCCA 的
    Schur 分解与宏状态。这段的价格是任何核方案的下限，先量出来供决策。

输入（全部只读，不复算任何 PCA/HVG，避免与已签口径分叉）：
    results/04_integration/seurat_trad/epiA/embeddings/cells.txt      133,384 barcode
    results/04_integration/seurat_trad/epiA/embeddings/harmony_f32.bin 133384×50 f32
    results/04_integration/seurat_trad/epiA/embeddings/pca_f32.bin     同上（未校正臂）
    results/04_integration/seurat_io/cell_meta.csv.gz                  仅用于取 stage 做占位

不产出：不做 fate probability、不做 driver gene、不做 stage 交叉表、
        不解释宏状态。宏状态只记个数，不记它是什么。
"""
import json
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd

EMB = "/home/eto/luad_v2/results/04_integration/seurat_trad/epiA/embeddings"
META = "/home/eto/luad_v2/results/04_integration/seurat_io/cell_meta.csv.gz"
OUT = "/home/eto/luad_v2/results/09_trajectory/smoke"
N_HVG_REP = 50
K_PARAM = 20          # 与 run_manifest.json 的 FindNeighbors k.param 一致
SCHUR_N = 20
N_STATES = 10


def hwm_gb():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmHWM:"):
                return round(int(line.split()[1]) / 1024 / 1024, 2)
    return None


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class Step:
    """把每一步的墙钟、峰值内存、异常都记下来，单步失败不炸整体。"""

    def __init__(self):
        self.timings = []
        self.errors = []

    def __call__(self, name, fn):
        t0 = time.time()
        try:
            val = fn()
        except Exception as exc:            # noqa: BLE001 - 冒烟就是要抓住一切异常
            dt = round(time.time() - t0, 1)
            self.timings.append({"step": name, "wall_sec": dt, "status": "ERROR"})
            self.errors.append({"step": name, "type": type(exc).__name__,
                                "msg": str(exc)[:500],
                                "trace_tail": traceback.format_exc().strip().splitlines()[-3:]})
            log(f"  !! {name} 失败 ({type(exc).__name__}) {dt}s: {str(exc)[:200]}")
            return None
        dt = round(time.time() - t0, 1)
        self.timings.append({"step": name, "wall_sec": dt, "status": "ok",
                             "peak_rss_gb": hwm_gb()})
        log(f"  ok {name}: {dt}s  peak_rss={hwm_gb()} GB")
        return val


def load_arm(arm):
    """读 frozen 嵌入 → AnnData。arm ∈ {harmony, pca}。"""
    import anndata as ad

    cells = [ln.strip() for ln in open(os.path.join(EMB, "cells.txt")) if ln.strip()]
    fname = "harmony_f32.bin" if arm == "harmony" else "pca_f32.bin"
    X = np.fromfile(os.path.join(EMB, fname), dtype=np.float32)
    if X.size != len(cells) * N_HVG_REP:
        raise ValueError(f"{fname}: {X.size} != {len(cells)}×{N_HVG_REP}")
    X = X.reshape(len(cells), N_HVG_REP)

    meta = pd.read_csv(META)
    stage = meta.set_index("cell_barcode")["stage"]
    stage = stage.reindex(cells)
    n_missing = int(stage.isna().sum())

    obs = pd.DataFrame(index=pd.Index(cells, name="cell_barcode"))
    # 只做占位，本冒烟不消费它（交叉表属生物学结论，口径未签字）
    obs["stage"] = pd.Categorical(stage.values,
                                  categories=["Normal", "AAH", "AIS", "MIA", "IAC"])
    adata = ad.AnnData(X=X, obs=obs)
    adata.obsm[f"X_{arm}"] = X.copy()
    return adata, n_missing


def run_arm(arm, st):
    import scanpy as sc
    import cellrank as cr

    log(f"=== 臂 {arm} ===")
    adata, n_missing = load_arm(arm)
    log(f"  载入 {adata.n_obs} 细胞 × {adata.n_vars} 维；stage 缺失 {n_missing}")

    def build_graph():
        sc.pp.neighbors(adata, n_neighbors=K_PARAM, use_rep=f"X_{arm}",
                        metric="euclidean")
        return adata

    if st("knn_graph", build_graph) is None:
        return {"arm": arm, "aborted_at": "knn_graph"}

    kern = st("connectivity_kernel", lambda: cr.kernels.ConnectivityKernel(adata)
              .compute_transition_matrix())
    if kern is None:
        return {"arm": arm, "aborted_at": "connectivity_kernel"}

    # CellRank 2 里过渡矩阵写在 adata.obsp[f"T_{kernel.key}"]，估计器不会自动认领。
    # 内核须作为**位置参数**传给估计器（KernelMixin.__init__(kernel=...)），
    # 传 kernels= 会被路由到 _from_adata 上、报 TypeError。
    def init_est():
        try:
            return cr.estimators.GPCCA(kern)
        except Exception:                    # 回退：把矩阵塞进 obsp 再给 AnnData
            adata.obsp["T_fwd"] = kern.transition_matrix
            return cr.estimators.GPCCA(adata)

    g = st("init_estimator", init_est)
    if g is None:
        return {"arm": arm, "aborted_at": "init_estimator"}
    st("compute_schur", lambda: g.compute_schur(n_components=SCHUR_N))
    st("compute_macrostates", lambda: g.compute_macrostates(n_states=N_STATES))

    sizes = {}
    try:
        if "macrostates" in adata.obs:
            sizes = {str(k): int(v) for k, v in
                     adata.obs["macrostates"].value_counts().items()}
    except Exception as exc:                # noqa: BLE001
        st.errors.append({"step": "read_macrostate_sizes",
                          "type": type(exc).__name__, "msg": str(exc)[:300]})

    return {"arm": arm, "n_cells": int(adata.n_obs), "n_dims": int(adata.n_vars),
            "stage_missing": n_missing, "macrostate_sizes": sizes}


def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    st = Step()
    arms = []
    for arm in ("harmony", "pca"):
        try:
            arms.append(run_arm(arm, st))
        except Exception as exc:            # noqa: BLE001
            st.errors.append({"step": f"arm:{arm}:fatal",
                              "type": type(exc).__name__, "msg": str(exc)[:500],
                              "trace_tail": traceback.format_exc().strip().splitlines()[-3:]})
            log(f"  !! 臂 {arm} 整体失败: {exc}")

    manifest = {
        "kind": "cellrank_cost_smoke",
        "scope": "COMMON TRUNK ONLY — 不含方向源；不产出任何生物学结论",
        "why_no_direction_source": (
            "已批准的「病程做实时轴」需 RealTimeKernel，它要预计算的最优传输耦合"
            "（from_moscot/from_wot/显式 couplings/重写 compute_coupling），"
            "而 wot/moscot/ott/pot/scot 在本 venv(Python 3.8) 一个都没装。"),
        "input_embedding": f"{EMB}/{{harmony,pca}}_f32.bin",
        "graph": {"n_neighbors": K_PARAM, "metric": "euclidean",
                  "caveat": "Seurat 侧用 annoy；scanpy 侧是 umap/pynndescent NN，非同一实现"},
        "kernel": "ConnectivityKernel（无向骨架，不含任何方向/时间信息）",
        "estimator": {"name": "GPCCA", "n_schur_components": SCHUR_N,
                      "n_states_requested": N_STATES},
        "arms": arms,
        "timings": st.timings,
        "errors": st.errors,
        "total_wall_sec": round(time.time() - t_start, 1),
        "peak_rss_gb": hwm_gb(),
        "note": ("同一时刻后台在跑 P4 的 inferCNV 成本冒烟（R，~300% CPU，38 GB RSS），"
                 "本进程已 nice -n 19 且限 4 线程；墙钟**受轻度污染**，"
                 "只可作量级参考，不可当排产数字。"),
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(OUT, "cellrank_smoke_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)

    log(f"总计 {manifest['total_wall_sec']}s，峰值 {manifest['peak_rss_gb']} GB，"
        f"{len(st.errors)} 个错误")
    for e in st.errors:
        log(f"  ERROR {e['step']}: {e['type']}: {e['msg'][:200]}")


if __name__ == "__main__":
    sys.exit(main())
