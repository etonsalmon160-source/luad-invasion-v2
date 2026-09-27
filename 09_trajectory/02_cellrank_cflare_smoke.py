#!/usr/bin/env python
"""
CellRank CFLARE 成本冒烟 —— 目的只有一个：
    在**不装 PETSc/SLEPc** 的前提下，CFLARE 能不能在 133,384 细胞上把
    特征分解 + 命运概率跑完，以及花多少墙钟和内存。

为什么用 CFLARE 而不是 GPCCA：
    GPCCA 走 SchurMixin ⇒ 没 petsc4py/slepc4py 就强制降级 brandts ⇒
    稠密化 133,384² ≈ 142 GB（已实测吃到 152 GB，被 kill）。
    CFLARE **不在 SchurMixin 继承链**上，它调 scipy.sparse.linalg.eigs，
    源码原话 "Uses a sparse implementation, if possible, and only computes the
    top k eigenvectors" ⇒ 不稠密化。

⚠️ 本次**不给任何生物学解读**：
    用的是 ConnectivityKernel —— 无向、不含方向/时间信息。
    有向核（RealTimeKernel 要 OT 耦合；PseudotimeKernel 要一个签过字的根）
    都属未签字口径 ⇒ 这里**不选方向源**。
    因此本次算出的 initial/terminal states 与 fate probabilities **没有意义**，
    只是用来量代价。脚本故意不打印、不落盘它们的任何内容。

安全闸：外部 watchdog 会在 RSS 超过阈值时直接杀进程（见 run 脚本）。
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
N_DIMS = 50
K_PARAM = 20
K_EIG = 20


def hwm_gb():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmHWM:"):
                return round(int(line.split()[1]) / 1024 / 1024, 2)
    return None


def rss_gb():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return round(int(line.split()[1]) / 1024 / 1024, 2)
    return None


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] rss={rss_gb()}GB {msg}", flush=True)


class Step:
    def __init__(self):
        self.timings, self.errors = [], []

    def __call__(self, name, fn):
        t0 = time.time()
        try:
            val = fn()
        except Exception as exc:            # noqa: BLE001
            dt = round(time.time() - t0, 1)
            self.timings.append({"step": name, "wall_sec": dt, "status": "ERROR"})
            self.errors.append({"step": name, "type": type(exc).__name__,
                                "msg": str(exc)[:400],
                                "trace_tail": traceback.format_exc().strip().splitlines()[-3:]})
            log(f"  !! {name} 失败 ({type(exc).__name__}) {dt}s: {str(exc)[:200]}")
            return None
        dt = round(time.time() - t0, 1)
        self.timings.append({"step": name, "wall_sec": dt, "status": "ok",
                             "peak_rss_gb": hwm_gb(), "rss_gb": rss_gb()})
        log(f"  ok {name}: {dt}s peak={hwm_gb()}GB")
        return val


def load_arm(arm):
    import anndata as ad

    cells = [ln.strip() for ln in open(os.path.join(EMB, "cells.txt")) if ln.strip()]
    fname = "harmony_f32.bin" if arm == "harmony" else "pca_f32.bin"
    X = np.fromfile(os.path.join(EMB, fname), dtype=np.float32)
    if X.size != len(cells) * N_DIMS:
        raise ValueError(f"{fname}: {X.size} != {len(cells)}×{N_DIMS}")
    X = X.reshape(len(cells), N_DIMS)

    stage = pd.read_csv(META).set_index("cell_barcode")["stage"].reindex(cells)
    obs = pd.DataFrame(index=pd.Index(cells, name="cell_barcode"))
    obs["stage"] = pd.Categorical(stage.values,
                                  categories=["Normal", "AAH", "AIS", "MIA", "IAC"])
    adata = ad.AnnData(X=X, obs=obs)
    adata.obsm[f"X_{arm}"] = X.copy()
    return adata, int(stage.isna().sum())


def sparse_report(adata, st):
    """把转移矩阵的稀疏程度记下来——这是"会不会稠密化"的直接证据。"""
    import scipy.sparse as sp

    T = None
    for key in ("T_connectivities", "T_fwd"):
        if key in adata.obsp:
            T = adata.obsp[key]
            break
    if T is None:
        return {}
    info = {"key": key, "is_sparse": bool(sp.issparse(T))}
    if sp.issparse(T):
        info["nnz"] = int(T.nnz)
        info["dense_gb_if_materialised"] = round(
            (T.shape[0] * T.shape[1] * 8) / 1024 ** 3, 1)
    return info


def run_arm(arm, st):
    import scanpy as sc
    import cellrank as cr

    log(f"=== 臂 {arm} ===")
    adata, n_missing = load_arm(arm)
    log(f"  载入 {adata.n_obs} 细胞 × {adata.n_vars} 维；stage 缺失 {n_missing}")

    # 注意：这些函数必须**显式 return 一个非 None 的对象**。
    # 守卫写的是 `if st(...) is None: 提前返回`，若被调函数忘了 return，
    # 就会被误判成"这一步失败"而静默跳过全部重活（2026-09-25 真踩过这个坑）。
    def build_graph():
        sc.pp.neighbors(adata, n_neighbors=K_PARAM, use_rep=f"X_{arm}", metric="euclidean")
        return adata

    if st("knn_graph", build_graph) is None:
        return {"arm": arm, "aborted_at": "knn_graph"}

    def build_kernel():
        k = cr.kernels.ConnectivityKernel(adata)
        k.compute_transition_matrix()
        return k

    kern = st("connectivity_kernel", build_kernel)
    if kern is None:
        return {"arm": arm, "aborted_at": "connectivity_kernel"}

    spinfo = sparse_report(adata, st)
    log(f"  转移矩阵: {spinfo}")

    def init_est():
        try:
            return cr.estimators.CFLARE(kern)
        except Exception:
            adata.obsp["T_fwd"] = kern.transition_matrix
            return cr.estimators.CFLARE(adata)

    g = st("init_estimator", init_est)
    if g is None:
        return {"arm": arm, "aborted_at": "init_estimator", "transition_matrix": spinfo}

    def compute_eig():
        g.compute_eigendecomposition(k=K_EIG, which="LR")
        return g

    if st("compute_eigendecomposition", compute_eig) is None:
        return {"arm": arm, "aborted_at": "compute_eigendecomposition",
                "transition_matrix": spinfo}

    # 只记特征值排序后的"形状"（头尾各几个），不记细胞归属
    eig_shape = {}
    try:
        D = g.eigendecomposition["D"]
        eig_shape = {"k": int(len(D)),
                     "top5_real": [round(float(x.real), 4) for x in D[:5]],
                     "bottom3_real": [round(float(x.real), 4) for x in D[-3:]]}
    except Exception as exc:                # noqa: BLE001
        st.errors.append({"step": "read_eig", "type": type(exc).__name__, "msg": str(exc)[:300]})

    def do_fit():
        return g.fit(k=K_EIG)

    if st("fit(terminal_states)", do_fit) is None:
        return {"arm": arm, "aborted_at": "fit", "transition_matrix": spinfo,
                "eigendecomposition": eig_shape}

    def do_fate():
        g.predict()
        return g

    if st("predict+fate_probabilities", do_fate) is None:
        return {"arm": arm, "aborted_at": "fate_probabilities", "transition_matrix": spinfo,
                "eigendecomposition": eig_shape}

    # 只数个数，不看是什么
    n_term = None
    try:
        n_term = int(len(g.terminal_states.cat.categories))
    except Exception:                       # noqa: BLE001
        pass

    return {"arm": arm, "n_cells": int(adata.n_obs), "n_dims": int(adata.n_vars),
            "stage_missing": n_missing, "transition_matrix": spinfo,
            "eigendecomposition": eig_shape, "n_terminal_states": n_term}


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    st = Step()
    arms = []
    for arm in ("harmony", "pca"):
        try:
            arms.append(run_arm(arm, st))
        except Exception as exc:            # noqa: BLE001
            st.errors.append({"step": f"arm:{arm}:fatal", "type": type(exc).__name__,
                              "msg": str(exc)[:400],
                              "trace_tail": traceback.format_exc().strip().splitlines()[-3:]})
            log(f"  !! 臂 {arm} 整体失败: {exc}")

    manifest = {
        "kind": "cellrank_cflare_cost_smoke",
        "scope": "NO BIOLOGICAL READOUT — 无向核，方向源未签字；terminal states/fate probs 无意义，不落盘",
        "why_cflare": ("GPCCA 走 SchurMixin ⇒ 无 petsc4py/slepc4py 即降级 brandts ⇒ 稠密化 "
                       "133,384²≈142GB（已实测 152GB 被 kill）。CFLARE 不在 SchurMixin 链上，"
                       "调 scipy.sparse.linalg.eigs，保持稀疏。"),
        "kernel": "ConnectivityKernel（无向）—— 只为量代价，不含方向",
        "transition_matrix_caveat": "dense_gb_if_materialised 是用来说明『一旦稠密化就有多大』",
        "estimator": {"name": "CFLARE", "k_eig": K_EIG, "which": "LR"},
        "graph": {"n_neighbors": K_PARAM, "metric": "euclidean",
                  "caveat": "Seurat 用 annoy；scanpy 用 umap/pynndescent NN，非同一实现"},
        "arms": arms,
        "timings": st.timings,
        "errors": st.errors,
        "total_wall_sec": round(time.time() - t0, 1),
        "peak_rss_gb": hwm_gb(),
        "contamination": "同时段后台在跑 P4 inferCNV（R，~300% CPU，38GB），墙钟受污染，只作量级参考",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(OUT, "cellrank_cflare_smoke_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)

    log(f"总计 {manifest['total_wall_sec']}s，峰值 {manifest['peak_rss_gb']}GB，"
        f"{len(st.errors)} 个错误")
    for e in st.errors:
        log(f"  ERROR {e['step']}: {e['type']}: {e['msg'][:200]}")


if __name__ == "__main__":
    sys.exit(main())
