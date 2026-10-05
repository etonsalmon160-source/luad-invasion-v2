#!/usr/bin/env python3
# 25_tr_dump_for_crosscheck.py —— 靶点扰动逆向臂：**为独立对拍导出子集**
#
# 目的：把「我们的矩阵 + 我们的 WTCS/NCS」原样导出，交给 R 侧用**官方实现**
#       （`signatureSearch:::.enrichScore` / `.lincsScores`）重算一遍 ⇒ 逐位比对。
#
# 🔴 为什么必须做这一步：本臂的 WTCS/NCS/τ 是**我们自己实现的**。
#    ES 已与 `fgsea` 对拍（1.1e-16），但 **NCS/τ 没有任何独立实现验证过**，
#    而 NCS 上确实出过一个会静默产生假阴性的符号 bug。
#
# 口径（全部显式，出处见 TARGET_REVERSAL_PARAMETERS.md）：
#   p (权重指数) = 1        ← GSEA 默认（非 CMap 明示）
#   ε            = 1e-4     ← 自定（本脚本不涉及）
#   NCS 分组      = cell_id × pert_type，**正负分开取均值**，μ⁻ 取绝对值
#   零 WTCS      = 从 μ 的计算中剔除（与官方一致）
#
# 产物（全部写 <outdir>）：
#   X_f64.bin   行主序 float64，(nsig × ngene)
#   genes.txt   每行一个基因名（列序）
#   sigs.tsv    sig_id / cell_id / pert_type（行序）
#   query.txt   q_up 一行、q_down 一行（基因名，空格分隔）
#   ours.tsv    sig_id / wtcs / ncs   ← 我们算的
#
# 用法：python3 25_tr_dump_for_crosscheck.py --ncol 2000 --outdir results/10_niche/target_reversal/xcheck

import argparse, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib.machinery import SourceFileLoader
TR = SourceFileLoader("tr", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "24_target_reversal.py")).load_module()

D_LINCS = "/home/eto/lincs_data"

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncol", type=int, default=2000)
    ap.add_argument("--query", default="D3")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    from cmapPy.pandasGEXpress.parse_gctx import parse
    G = parse(f"{D_LINCS}/GSE70138_Level5.gctx")
    df = G.data_df
    genes = [str(g) for g in df.index]                       # Entrez
    gi = TR.load_gene_map("GSE70138")
    sym2e = dict(zip(gi.pr_gene_symbol.astype(str), gi.pr_gene_id.astype(str)))
    gpos = {g: i for i, g in enumerate(genes)}

    up_sym, dn_sym = TR.load_query(int(a.query.lstrip("D")))
    up_ix = np.array([gpos[sym2e[g]] for g in up_sym if sym2e.get(g) in gpos])
    dn_ix = np.array([gpos[sym2e[g]] for g in dn_sym if sym2e.get(g) in gpos])
    log(f"查询 {a.query}：q_up {len(up_ix)} / q_down {len(dn_ix)}（Entrez 行号）")

    ncol = min(a.ncol, df.shape[1])
    sig_ids = [str(c) for c in df.columns[:ncol]]
    # 🔴 用 **float64**：源头 gctx 是 float32，但对拍必须两边同精度，否则残差 ~1e-7
    X = df.iloc[:, :ncol].values.T.astype(np.float64)        # nsig × ngene
    si = TR.load_sig_info("GSE70138").set_index("sig_id")
    meta = si.reindex(sig_ids)
    cell = meta.cell_id.fillna("NA").astype(str).values
    ptyp = meta.pert_type.fillna("NA").astype(str).values
    log(f"子集：{X.shape[0]} 签名 × {X.shape[1]} 基因")

    # —— 我们自己的 WTCS / NCS ——
    rank, sv, N = TR.build_rank_tables(X)
    R,C = TR._rc(rank, sv, up_ix); eup = TR.es_from_ranks(R, C, N, len(up_ix))
    R,C = TR._rc(rank, sv, dn_ix); edn = TR.es_from_ranks(R, C, N, len(dn_ix))
    wtcs = np.where(np.sign(eup) != np.sign(edn), (eup - edn) / 2.0, 0.0)
    ncs = TR_ncs(wtcs, cell, ptyp)
    log(f"我们：WTCS 非零 {100*(wtcs!=0).mean():.1f}% ; NCS 范围 [{ncs.min():.3f},{ncs.max():.3f}]")

    # —— 导出 ——
    d = a.outdir
    X.tofile(os.path.join(d, "X_f64.bin"))
    with open(os.path.join(d, "genes.txt"), "w") as f: f.write("\n".join(genes) + "\n")
    pd.DataFrame({"sig_id": sig_ids, "cell_id": cell, "pert_type": ptyp}).to_csv(
        os.path.join(d, "sigs.tsv"), sep="\t", index=False)
    eu = [g for g in up_sym if sym2e.get(g) in gpos]
    ed = [g for g in dn_sym if sym2e.get(g) in gpos]
    with open(os.path.join(d, "query.txt"), "w") as f:
        f.write(" ".join(TR_entrez(eu, sym2e)) + "\n" + " ".join(TR_entrez(ed, sym2e)) + "\n")
    pd.DataFrame({"sig_id": sig_ids, "wtcs": wtcs, "ncs": ncs}).to_csv(
        os.path.join(d, "ours.tsv"), sep="\t", index=False)
    log(f"导出完成 → {d}")

def TR_entrez(syms, sym2e): return [sym2e[g] for g in syms]

def TR_ncs(w, cell, ptyp):
    """与 24_target_reversal.py 的 ncs_of 完全同式（正负分取均值、μ⁻ 取绝对值、零剔除）。"""
    out = np.empty_like(w)
    key = np.array([f"{c}|{t}" for c, t in zip(cell, ptyp)])
    for g in np.unique(key):
        m = key == g
        ww = w[m]
        mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
        mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
        out[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
    return out

if __name__ == "__main__":
    main()
