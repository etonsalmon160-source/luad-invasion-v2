#!/usr/bin/env python
"""
微型复现：用**已知行和**的假耦合，检查 RealTimeKernel 装配时到底归不归一化。

为什么要做（这是纠错，不是新发现）：
    10_realtime_kernel_smoke.py 量出「逐期行和 = 1.0（min=max=1.0）」，但我用 h5py
    直接读耦合文件，行和是 0.63–1.44，**不是 1**。两者不可能同时为真。块结构又量得
    「前四期对角块 nnz = 0」⇒ 前四期每一行只装耦合 ⇒ 行和必须等于耦合行和。
    ⇒ 必有一处是我自己的错。这个小实验把整条链路放到我能控制行和的数据上：
    耦合行和是我自己设的（0.6 + 0.4·U），如果跑完变成 1.0，就是 cellrank 归一化；
    如果保持原值，就是我 10 号脚本的测量代码有问题。

    顺带这也**独立验证**了块上三角那条结论：严格下三角必须为 0。
"""
import json
import os
import tempfile

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]
N_PER = [40, 30, 35, 20, 45]           # 故意不等，和真实数据的形态类似
OUT = "/home/eto/luad_v2/results/09_trajectory/smoke/microcheck_realtime_kernel.json"
TMAPDIR = "/home/eto/luad_v2/results/09_trajectory/smoke/microcheck_tmaps"

rng = np.random.default_rng(20260925)


def main():
    os.makedirs(TMAPDIR, exist_ok=True)
    for f in os.listdir(TMAPDIR):                     # 目录必须干净（from_wot 会 glob 全部 h5ad）
        os.remove(os.path.join(TMAPDIR, f))

    names, days, stages = [], [], []
    for i, s in enumerate(STAGES):
        for j in range(N_PER[i]):
            names.append(f"{s}_{j}")
            days.append(float(i))
            stages.append(s)
    obs = pd.DataFrame({"stage": stages,
                        "day": pd.Categorical(days, categories=[0.0, 1.0, 2.0, 3.0, 4.0])},
                       index=names)
    a = ad.AnnData(X=sp.csr_matrix((len(names), 5), dtype=np.float32), obs=obs)

    # 造 4 个耦合：行和 = 0.6 + 0.4·U(0,1)（**确定不是 1**），列向随机摊
    truth = {}
    for i in range(4):
        src = [n for n, d in zip(names, days) if d == i]
        tgt = [n for n, d in zip(names, days) if d == i + 1]
        P = rng.random((len(src), len(tgt)))
        P = P / P.sum(1, keepdims=True)
        rs = 0.6 + 0.4 * rng.random(len(src))
        P = P * rs[:, None]
        c = ad.AnnData(X=np.ascontiguousarray(P, dtype=np.float64),
                       obs=pd.DataFrame(index=src), var=pd.DataFrame(index=tgt))
        f = os.path.join(TMAPDIR, f"tmaps_{float(i)}_{float(i+1)}.h5ad")
        c.write_h5ad(f)
        truth[(float(i), float(i + 1))] = rs
        print(f"  写 {os.path.basename(f)}  shape={c.shape}  行和 {rs.min():.4f}–{rs.max():.4f}")

    import cellrank  # noqa: F401
    from cellrank.kernels import RealTimeKernel

    rtk = RealTimeKernel.from_wot(a, path=TMAPDIR, time_key="day")
    print(f"  解析出的键：{sorted(rtk.couplings.keys())}；reference={rtk._reference}")
    rtk.compute_transition_matrix(self_transitions="diagonal", threshold=None)
    M = rtk.transition_matrix.tocsr()

    si = a.obs["stage"].astype(str).map({s: i for i, s in enumerate(STAGES)}).values
    perm = np.argsort(si, kind="stable")
    coo = M[perm, :][:, perm].tocoo()
    sr = si[perm]
    strict_lower = int(np.sum(sr[coo.row] > sr[coo.col]))
    diagblk = {s: int(np.sum((sr[coo.row] == i) & (sr[coo.col] == i)))
               for i, s in enumerate(STAGES)}

    rep = {"n_per_stage": dict(zip(STAGES, N_PER)),
           "truth_rowsum_minmax": {f"{k[0]}->{k[1]}": [round(float(v.min()), 4),
                                                       round(float(v.max()), 4)]
                                   for k, v in truth.items()},
           "assembled_rowsum_by_stage": {}, "strict_lower_nnz": strict_lower,
           "diag_block_nnz": diagblk, "shape": list(M.shape), "nnz": int(M.nnz)}

    rs = np.asarray(M.sum(axis=1)).ravel()
    print(f"\n  装配后：严格下三角非零 = {strict_lower}（应 0）；对角块 = {diagblk}")
    for i, s in enumerate(STAGES):
        g = rs[si == i]
        rep["assembled_rowsum_by_stage"][s] = {"min": round(float(g.min()), 6),
                                              "max": round(float(g.max()), 6),
                                              "median": round(float(np.median(g)), 6)}
        print(f"  {s:7s} 装配后行和 {g.min():.4f}–{g.max():.4f}   "
              f"原始耦合行和 " + ("  （无耦合，自转块）" if i == 4 else
                                  f"{truth[(float(i), float(i+1))].min():.4f}–"
                                  f"{truth[(float(i), float(i+1))].max():.4f}"))

    # 前四期：装配后行和是否**逐字节等于**原始耦合行和？
    #   若 True ⇒ cellrank 不归一化；若 False ⇒ 被归一化（本次结论是 False）
    ok = True
    for i in range(4):
        t = truth[(float(i), float(i + 1))]
        g = rs[si == i]
        same = np.allclose(np.sort(t), np.sort(g), rtol=0, atol=1e-12)
        ok &= same
        print(f"  {STAGES[i]:7s} 装配行和 == 原始耦合行和（位置对应）？ {same}"
              f"   ← 预期 False（说明装配被逐行归一化）")
    rep["rowsum_preserved_exactly"] = bool(ok)

    # 第五期（IAC）自转块：diagonal ⇒ 单位阵 ⇒ 行和恰 1
    g = rs[si == 4]
    rep["IAC_selfblock_rowsum"] = {"min": round(float(g.min()), 6), "max": round(float(g.max()), 6)}
    print(f"  IAC（自转块 = 单位阵）行和 {g.min():.6f}–{g.max():.6f}")

    rep["verdict"] = ("若 rowsum_preserved_exactly=False（实测就是 False），则装配结果**被逐行归一化**；"
                      "归一化发生在 Kernel.transition_matrix 的 setter（_base_kernel.py:400-436），"
                      "不是 _restich_couplings。逐行正缩放保零模式 ⇒ 块上三角不受影响。")

    # --- 严格验证：整矩阵特征值 == 各对角块特征值之并（det(A-λI)=Πdet(A_ii-λI)）---
    # 小矩阵，直接 eigvals 全谱，不用 ARPACK
    ev_full = np.linalg.eigvals(M.toarray())
    ev_blocks = np.concatenate([np.linalg.eigvals(
        M[np.ix_(np.where(si == i)[0], np.where(si == i)[0])].toarray())
        for i in range(len(STAGES))])
    o_f = np.argsort(-np.abs(ev_full))     # 先按模排序，**趁排序前**取前几个给报告用
    o_b = np.argsort(-np.abs(ev_blocks))
    top5_full = np.abs(ev_full[o_f])[:5]
    top5_union = np.abs(ev_blocks[o_b])[:5]
    # 多重集比较：两边都按同一规则排序后再逐位比（sort_complex 会按实部重排，故只用于比较）
    a = np.sort_complex(ev_full[o_f])
    b = np.sort_complex(ev_blocks[o_b])
    maxdiff = float(np.max(np.abs(a - b)))
    n_zero = int(np.sum(np.abs(ev_full) < 1e-10))
    rep["spectrum_union_check"] = {
        "n_eigenvalues_full": int(ev_full.size),
        "max_abs_diff_vs_union_of_diag_blocks": maxdiff,
        "n_eigenvalues_near_zero": n_zero,
        "n_cells_in_first_four_stages": int(sum(N_PER[:4])),
        "top5_abs_full": [round(float(x), 6) for x in top5_full],
        "top5_abs_union": [round(float(x), 6) for x in top5_union],
    }
    print(f"\n  谱并集检验：整矩阵 {ev_full.size} 个特征值 vs 对角块并集，"
          f"最大逐位差 = {maxdiff:.3e}")
    print(f"             |λ|≈0 的个数 = {n_zero}（前四期细胞数之和 = {sum(N_PER[:4])}）")
    print(f"             前5个 |λ|（整矩阵）= {[round(float(x),4) for x in top5_full]}")
    print(f"             前5个 |λ|（块并集）= {[round(float(x),4) for x in top5_union]}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(rep, fh, ensure_ascii=False, indent=1, default=str)
    print(f"\n写到 {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
