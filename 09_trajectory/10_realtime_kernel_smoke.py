#!/usr/bin/env python
"""
机械检验：把已签字的 WOT 耦合接成 RealTimeKernel，看它的谱结构到底是什么样。

为什么这么设计（读源码后决定的，不是猜）：
    读 cellrank/kernels/_real_time_kernel.py 的 _restich_couplings()：
      - 矩阵按期别分 5×5 块，顺序 = 期别顺序（0..4）
      - 耦合只放在**超对角**块 blocks[i][i+1]
      - 所有对角块**显式设成全零**：`blocks[src_ix][src_ix] = sp.spdiags([0]*n, 0, n, n)`
      - **只有 reference（= 最大的期别键 = IAC 4.0）那一块**被填成
        identity / uniform / connectivities
    ⇒ 这是个**块上三角**矩阵。块上三角矩阵的特征值 = 各对角块特征值的并集
      （det(A−λI) = Π det(A_ii−λI)）。
    ⇒ 推论：**谱 = {0（重数 = 前四期细胞数之和）} ∪ spec(IAC 自转块)**。
      运输耦合的**数值不进特征值**。

    ⚠️ 一处**我原先写错、后经 11_realtime_kernel_microcheck.py 实验推翻**的地方：
      我曾写「compute_transition_matrix() 结尾是 `self.transition_matrix = tmap.X`，没有行归一化」。
      事实是 `transition_matrix` 是 `Kernel` 基类的 **property，setter（`_base_kernel.py:400-436`）
      会逐行归一化**（`should_norm()` 不过就 `_normalize()`）。
      微型复现里我造的行和 0.61–0.96 出来正好变成 1.0 ⇒ 归一化确实发生。
      **但这不改变上面的结论**：行归一化是逐行乘一个正的标量（D⁻¹A），
      **保零模式** ⇒ 块上三角依旧是块上三角，对角块依旧是零块（IAC 那块本来行和就是 1）。
      代价是 WOT 的行质量（0.63–1.44，即生长率信息）被这一步丢掉了。

    这个推论必须实测，不能只靠读码。所以这里做三件事：
      1) 量出块结构：按期别重排后，严格下三角块的非零个数必须为 0
      2) 量对角块：前四期的对角块必须全零；IAC 那块非零
      3) 用两种方式求谱并对照：(a) 只对 IAC 块求，(b) 对**整个 133,384² 矩阵**求前几个
         —— 两者应当一致。这一步是推论的经验确认。
      4) （补）`self_transitions='all'` 是**唯一**让矩阵不再三角的设定：每一期的对角块
         都被填上连通性、耦合被按 (1-conn_weight) 缩放 ⇒ 谱里才含运输信息。
         那一个配置用 RK_CONFIGS 环境变量单独跑。

只做结构检验，不给任何生物学结论。
"""
import json
import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

IN = "/home/eto/luad_v2/results/09_trajectory/wot_full/epiA_wot_input.h5ad"
TMAP = "/home/eto/luad_v2/results/09_trajectory/wot_full/tmaps"
EMB = "/home/eto/luad_v2/results/04_integration/seurat_trad/epiA/embeddings"
OUT = os.environ.get("RK_OUT",
                     "/home/eto/luad_v2/results/09_trajectory/smoke/realtime_kernel_structure.json")

# 全部可选配置。为什么是这个集合：
#   diagonal / connectivities 是默认路径 —— 源码里除 reference(=IAC) 外，所有对角块被
#   显式设成全零 ⇒ 块上三角 ⇒ 谱与运输数值无关。这两个配置用来确认那个退化。
#   all 则是真正混合：**每一期**的对角块都被填上连通性、耦合按 (1-conn_weight) 缩放
#   ⇒ 矩阵不再三角，谱里才含有运输与期内的混合信息。这是唯一谱有意义的设定。
ALL_CONFIGS = {
    "diagonal_noThresh": dict(self_transitions="diagonal", threshold=None),
    "diagonal_autoThresh": dict(self_transitions="diagonal", threshold="auto"),
    "connectivities_harmony": dict(self_transitions="connectivities", threshold="auto",
                                   conn_kwargs={"use_rep": "X_pca_harmony",
                                                "n_neighbors": 20}),
    "all_harmony_w05_noThresh": dict(self_transitions="all", conn_weight=0.5, threshold=None,
                                     conn_kwargs={"use_rep": "X_pca_harmony",
                                                  "n_neighbors": 20}),
}
SELECTED = os.environ.get("RK_CONFIGS",
                          "diagonal_noThresh,diagonal_autoThresh,connectivities_harmony").split(",")
K_FULL = int(os.environ.get("RK_K_FULL", "6"))

STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]
MEM_LIMIT_GB = 150.0
K_TOP = 20


def rss_gb():
    try:
        with open("/proc/self/status") as fh:
            for ln in fh:
                if ln.startswith("VmRSS:"):
                    return int(ln.split()[1]) / 1024 / 1024
    except Exception:                                  # noqa: BLE001
        pass
    return 0.0


def hwm_gb():
    try:
        with open("/proc/self/status") as fh:
            for ln in fh:
                if ln.startswith("VmHWM:"):
                    return int(ln.split()[1]) / 1024 / 1024
    except Exception:                                  # noqa: BLE001
        pass
    return 0.0


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] rss={rss_gb():.2f}GB {m}", flush=True)


def guard(where):
    """自己看自己的内存，不用外部闸、不用 PID —— 踩过一次 PID 的坑了。"""
    r = rss_gb()
    if r > MEM_LIMIT_GB:
        raise MemoryError(f"{where}: RSS {r:.1f}GB 越过自设上限 {MEM_LIMIT_GB}GB，主动停")


def top_eigs(M, k=K_TOP, tag=""):
    """按模取前 k 个特征值。M 非对称 ⇒ 用 eigs（不是 eigsh）。

    收敛失败只登记、不中断 —— 否则一个配置失败会连带丢掉全部块结构结果。
    """
    t = time.time()
    n = M.shape[0]
    kk = min(k, n - 2)
    ncv = min(max(2 * kk + 1, 20), n - 1)
    try:
        vals = sp.linalg.eigs(M, k=kk, which="LM", ncv=ncv, tol=1e-7,
                              maxiter=3000, return_eigenvectors=False)
    except Exception as exc:                              # noqa: BLE001
        log(f"    {tag} eigs 未收敛/失败：{type(exc).__name__}: {str(exc)[:160]}"
            f"（{time.time()-t:.1f}s）")
        return {"k": int(kk), "ncv": int(ncv), "seconds": round(time.time() - t, 1),
                "failed": f"{type(exc).__name__}: {str(exc)[:200]}", "abs_sorted": []}
    mag = np.sort(np.abs(vals))[::-1]
    real_frac = float(np.mean(np.abs(vals.imag) < 1e-8))
    log(f"    {tag} eigs(k={kk}) 用了 {time.time()-t:.1f}s；"
        f"前3个 = {', '.join(f'{m:.4f}' for m in mag[:3])}")
    return {"k": int(kk), "ncv": int(ncv), "seconds": round(time.time() - t, 1),
            "abs_sorted": [round(float(m), 6) for m in mag],
            "frac_nearly_real": round(real_frac, 3)}


def block_report(M, si, tag):
    """按期别重排后量块结构：严格下三角必须为 0。"""
    perm = np.argsort(si, kind="stable")
    Mp = M[perm, :][:, perm]
    sir = si[perm]
    coo = Mp.tocoo()
    rows, cols = coo.row, coo.col
    strict_lower = int(np.sum(sir[rows] > sir[cols]))
    diag_blocks = {}
    for i, s in enumerate(STAGES):
        m = sir == i
        # 直接按行列掩码数非零（比切片省事）
        sel = (sir[rows] == i) & (sir[cols] == i)
        diag_blocks[s] = {"n_cells": int(m.sum()), "nnz_onselfblock": int(sel.sum())}
    strict_upper = int(np.sum(sir[rows] < sir[cols]))
    log(f"  {tag} 块结构：严格下三角非零 = {strict_lower}（必须 0）；"
        f"上三角非零 = {strict_upper}；对角块 = "
        f"{ {k: v['nnz_onselfblock'] for k, v in diag_blocks.items()} }")
    del perm, Mp, coo
    return {"strict_lower_nnz": strict_lower, "strict_upper_nnz": strict_upper,
            "diag_blocks": diag_blocks}


def rowsum_report(M, si, tag):
    rs = np.asarray(M.sum(axis=1)).ravel()
    out = {s: {"median": round(float(np.median(rs[si == i])), 4),
               "min": round(float(rs[si == i].min()), 4),
               "max": round(float(rs[si == i].max()), 4)}
           for i, s in enumerate(STAGES)}
    log(f"  {tag} 逐期行和：{out}")
    return out


def main():
    import anndata as ad
    import cellrank as cr
    from cellrank.kernels import RealTimeKernel

    rep = {"purpose": "机械检验 RealTimeKernel 的谱结构；零生物学结论",
           "input": IN, "tmaps": TMAP}

    log("读入已签口径的输入 adata …")
    a = ad.read_h5ad(IN)
    log(f"  {a.n_obs} × {a.n_vars}；day dtype = {a.obs['day'].dtype}；"
        f"type(iloc[0]) = {type(a.obs['day'].iloc[0]).__name__}")
    rep["adata"] = {"n_obs": int(a.n_obs), "n_vars": int(a.n_vars),
                    "day_dtype": str(a.obs["day"].dtype),
                    "day_first_value_repr": repr(a.obs["day"].iloc[0]),
                    "day_first_value_type": type(a.obs["day"].iloc[0]).__name__}

    # 冻结嵌入：给连通性自转块用（与之前无向基线同一套、同一维度）
    cells = [x.strip() for x in open(os.path.join(EMB, "cells.txt")) if x.strip()]
    H = np.fromfile(os.path.join(EMB, "harmony_f32.bin"), dtype=np.float32)
    assert H.size % len(cells) == 0, "嵌入维度对不上"
    H = H.reshape(len(cells), -1)
    emb = pd.DataFrame(H, index=cells)
    missing = a.obs_names.difference(emb.index)
    if len(missing):
        raise ValueError(f"有 {len(missing)} 个细胞在冻结嵌入里找不到，停")
    a.obsm["X_pca_harmony"] = emb.loc[a.obs_names].values.copy()
    log(f"  冻结嵌入已接上：{a.obsm['X_pca_harmony'].shape}（{len(cells)} 个条码，0 缺失）")
    rep["embedding"] = {"path": EMB, "dim": int(a.obsm["X_pca_harmony"].shape[1]),
                        "missing_barcodes": 0}
    del H, emb

    si = a.obs["stage"].astype(str).map({s: i for i, s in enumerate(STAGES)}).values
    rep["stage_counts"] = {s: int((si == i).sum()) for i, s in enumerate(STAGES)}

    # 一个 kernel 对象，多次换参数（compute_transition_matrix 按参数缓存）
    guard("构造 kernel 前")
    log("用 from_wot 接上耦合 …")
    rtk = RealTimeKernel.from_wot(a, path=TMAP, time_key="day")
    keys = sorted(rtk.couplings.keys())
    log(f"  from_wot 解析出的耦合键（{len(keys)} 个）：{keys}")
    log(f"  reference（自转块落在哪一期）= {getattr(rtk, '_reference', None)}")
    rep["parsed_couplings"] = {
        "n": len(keys),
        "keys": [[str(k[0]), str(k[1])] for k in keys],
        "key_dtype": type(keys[0][0]).__name__,
        "reference": str(getattr(rtk, "_reference", None)),
        "shapes": {f"{k[0]}->{k[1]}": [int(rtk.couplings[k].n_obs), int(rtk.couplings[k].n_vars)]
                   for k in keys},
    }

    configs = [(k, ALL_CONFIGS[k]) for k in SELECTED if k]
    log(f"将跑 {len(configs)} 个配置：{[k for k, _ in configs]}")
    rep["selected_configs"] = [k for k, _ in configs]
    results = {}
    for name, kw in configs:
        guard(f"{name} 之前")
        log(f"=== 配置 {name}：{ {k: v for k, v in kw.items() if k != 'conn_kwargs'} } ===")
        t = time.time()
        rtk.compute_transition_matrix(**kw)
        M = rtk.transition_matrix.tocsr()
        log(f"  转移矩阵 {M.shape}，非零 {M.nnz:,}，构建 {time.time()-t:.1f}s")
        d = {"params": {k: (v if k != "conn_kwargs" else dict(v)) for k, v in kw.items()},
             "shape": list(M.shape), "nnz": int(M.nnz),
             "build_seconds": round(time.time() - t, 1),
             "density": round(M.nnz / (M.shape[0] * M.shape[1]), 8)}
        d["blocks"] = block_report(M, si, name)
        d["rowsum_by_stage"] = rowsum_report(M, si, name)

        # (a) 五个对角块各自的谱。
        #     若矩阵块上三角，则整矩阵谱 = 这五个块谱的并集（含 0 度块贡献的 0）。
        guard(f"{name} 求对角块谱前")
        d["spec_of_diag_blocks"] = {}
        for i, s in enumerate(STAGES):
            idx = np.where(si == i)[0]
            sub = M[np.ix_(idx, idx)].tocsr()
            d["spec_of_diag_blocks"][s] = top_eigs(sub, k=3, tag=f"{name} 对角块 {s}")
            del sub
        iac = STAGES.index("IAC")
        idx_iac = np.where(si == iac)[0]
        D4 = M[np.ix_(idx_iac, idx_iac)].tocsr()
        d["spec_of_IAC_block"] = top_eigs(D4, k=K_TOP, tag=f"{name} IAC块")
        del D4

        # (b) 对整矩阵求前几个，作为对「特征值=对角块并集」这个推论的经验确认
        guard(f"{name} 整矩阵求谱前")
        d["spec_of_full_matrix"] = top_eigs(M, k=K_FULL, tag=f"{name} 整矩阵")
        d["_full_k"] = K_FULL

        results[name] = d
        del M
        import gc
        gc.collect()
        log(f"  {name} 完成，已释放矩阵；峰值 {hwm_gb():.2f}GB")

    rep["configs"] = results
    rep["verdict_note"] = (
        "若 blocks.strict_lower_nnz == 0 且前四期对角块 nnz == 0，则矩阵为块上三角，"
        "由 det(A-λI)=Πdet(A_ii-λI) 得 spec(A) = {0（重数=前四期细胞数）} ∪ spec(IAC自转块)。"
        "此时『有没有谱隙』完全由 IAC 期内部自转块决定，与运输耦合的数值无关。"
        "spec_of_full_matrix 与 spec_of_IAC_block 的前几个应当一致，即为该推论的经验确认。"
        "注意：装配结果会被 Kernel.transition_matrix 的 setter 逐行归一化"
        "（_base_kernel.py:400-436）；逐行正缩放保零模式，故块上三角结构不受影响。"
        "self_transitions='all' 是唯一让矩阵不再三角的设定，见 realtime_kernel_structure_all.json。")
    rep["peak_rss_gb"] = round(hwm_gb(), 2)
    rep["date"] = time.strftime("%Y-%m-%d %H:%M:%S")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(rep, fh, ensure_ascii=False, indent=1, default=str)
    log(f"写到 {OUT}；全程峰值 {hwm_gb():.2f}GB")

    print("\n================ 摘要 ================")
    for name, d in results.items():
        b = d["blocks"]
        print(f"{name}: nnz={d['nnz']:,} 严格下三角={b['strict_lower_nnz']} "
              f"对角块={ {k: v['nnz_onselfblock'] for k, v in b['diag_blocks'].items()} }")
        print(f"   IAC 块前5个|λ| = {d['spec_of_IAC_block']['abs_sorted'][:5]}")
        print(f"   整矩阵前5个|λ| = {d['spec_of_full_matrix']['abs_sorted'][:5]}")
    return 0


if __name__ == "__main__":
    import traceback
    try:
        raise SystemExit(main())
    except BaseException:
        traceback.print_exc()
        raise
