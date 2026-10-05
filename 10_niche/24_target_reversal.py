#!/usr/bin/env python3
# 24_target_reversal.py —— 靶点扰动逆向臂（TARGET-REVERSAL）主实现
#
# 预注册：10_niche/TARGET_REVERSAL_PREREG.md（已签）。本脚本**严格照 §5 方法学**。
#
# 逐式对应：
#   §5.3 WTCS  = 加权 KS（p=1，权重取 |统计量值|）→ (ES_up − ES_down)/2，异号才非零
#        ✅ 对拍三重：① vs fgsea 1.24.0（1.11e-16）；② vs 官方 signatureSearch:::.enrichScore
#        （逐位相同）；③ 见 26_tr_crosscheck_official.R
#   §5.4 NCS   = w/μ⁺（w>0）否则 w/μ⁻，μ = **cell_id × pert_type** 组内带符号均值
#   §5.5 τ     = sgn·100/N·Σ[|NCS_i,r| < |NCS_q,r|]，Q_ref = 1000 随机签名（**降级**，见 §5.5）
#   §5.6 零模型 = 1000 个随机查询签名（同基因数）
#
# 数据映射（已核实 2026-10-03）：
#   gctx 列名 = sig_info 的 **sig_id**；gctx 行 = gene_info 的 **pr_gene_id（Entrez）**
#   ⇒ 查询签名（基因符号）经 gene_info 的 pr_gene_symbol → pr_gene_id 转成 Entrez
#
# 用法：
#   --lib L1a|L1b|L2        选库
#   --ncol N                只跑前 N 列（小规模验证用）
#   --query D3|D5|D7        主查询（默认 D3）
#   --outdir PATH

import argparse, os, sys, time, json
import numpy as np
import pandas as pd

D_LINCS = "/home/eto/lincs_data"
KD      = "/home/eto/luad_v2/results/10_niche/kstar_diag"
OUT_BASE= "/home/eto/luad_v2/results/10_niche/target_reversal"
N_REF   = 1000        # Q_ref 随机查询数（§5.5 降级值；出处=signatureSearch 默认，非论文）
N_NULL  = 1000        # 零模型随机签名数（§5.6）
EPS     = 1e-4        # §5.12 登记为**自定参数**
SEED    = 20261003

def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ───────────────────────── §5.3 加权 KS（p=1）─────────────────────────
def weighted_ks(values_desc, is_hit):
    """参照实现（慢，用于对拍）：values_desc 按值降序，is_hit 同序。"""
    N = values_desc.size; S = int(is_hit.sum())
    if S == 0 or S >= N: return 0.0
    Nr = values_desc[is_hit].sum()
    if Nr <= 0: return 0.0
    cs = np.cumsum(np.where(is_hit, values_desc / Nr, -1.0 / (N - S)))
    return float(cs.max()) if cs.max() > -cs.min() else float(cs.min())

def es_exact(R, C, N, S):
    """**精确向量化 ES**：在完整增量序列上 cumsum（与官方 .enrichScore 逐位等价）。
       R = (nsig, S) 命中位置的升序秩；C = (nsig, S) 那些位置的 |值|。
       🔴 早期版本只取"命中处 + 命中前一刻"的极值做快捷运算，**有 ~1e-7 的数值误差**
          （与官方实现对拍才发现）；改为全序列 cumsum 后与官方**逐位相同**。
       内存：nsig 分块调用，单块内临时数组 nsig × N。"""
    nsig = R.shape[0]
    C = np.abs(C)                                  # 🔴 权重 = |值|，不是带符号的值
                                                   #    （漏了这句 ⇒ Nr 可为负 ⇒ 全 NA/全 0；我犯过两次）
    w = np.zeros((nsig, N), dtype=np.float64)
    h = np.zeros((nsig, N), dtype=np.float64)
    np.put_along_axis(w, R, C.astype(np.float64), axis=1)
    np.put_along_axis(h, R, 1.0, axis=1)
    Nr = C.sum(axis=1, keepdims=True)
    safe = Nr > 0                                  # 官方：NR==0 ⇒ 该签名 ES 记为 0
    Nr_s = np.where(safe, Nr, 1.0)
    inc = h * w / Nr_s - (1.0 - h) / (N - S)
    cs = np.cumsum(inc, axis=1)
    idx = np.argmax(np.abs(cs), axis=1)
    out = cs[np.arange(nsig), idx]
    return np.where(safe[:, 0], out, 0.0)

def es_from_ranks(R, C, N, S, chunk=2000):
    """分块调用 es_exact，控制内存。"""
    out = np.empty(R.shape[0])
    for a in range(0, R.shape[0], chunk):
        b = min(a + chunk, R.shape[0])
        out[a:b] = es_exact(R[a:b], C[a:b], N, S)
    return out

# ───────────────────────── 元数据 ─────────────────────────
def load_gene_map(series):
    gi = pd.read_csv(f"{D_LINCS}/{series}_Broad_LINCS_gene_info"
                     + ("_2017-03-06" if series == "GSE70138" else "") + ".txt.gz",
                     sep="\t", compression="gzip")
    return gi          # 列: pr_gene_id, pr_gene_symbol, pr_gene_title, pr_is_lm, pr_is_bing

def load_sig_info(series):
    si = pd.read_csv(f"{D_LINCS}/{series}_Broad_LINCS_sig_info"
                     + ("_2017-03-06" if series == "GSE70138" else "") + ".txt.gz",
                     sep="\t", compression="gzip")
    return si          # sig_id, pert_id, pert_iname, pert_type, cell_id, ...

def load_query(archetype=3, top=150):
    """§5.2：从 d14_signed_panel.tsv 取 q_up / q_down（基因符号）。"""
    P = pd.read_csv(f"{KD}/d14_signed_panel.tsv", sep="\t")
    P = P[P.archetype == archetype]
    up = P[(P.direction == "up")].sort_values("rank").gene.tolist()[:top]
    dn = P[(P.direction == "down")].sort_values("rank").gene.tolist()[:top]
    return up, dn

# ───────────────────────── 主流程 ─────────────────────────

def _rc(rank, sv, q_idx):
    """返回 (R, C)：命中位置的升序秩 与 |值|。"""
    R = np.sort(rank[:, q_idx], axis=1)
    C = np.take_along_axis(sv, R, axis=1)
    return R, C

def build_rank_tables(X):
    """预算：order / rank / sorted_vals（一次，之后所有查询 O(1) 次向量化）。"""
    nsig, N = X.shape
    order = np.argsort(-X, axis=1, kind="stable").astype(np.int32)
    rank  = np.empty_like(order)
    np.put_along_axis(rank, order, np.broadcast_to(np.arange(N, dtype=np.int32), (nsig, N)), axis=1)
    sv = np.take_along_axis(X, order, axis=1)   # 🔴 不强制降精度：保持输入 dtype（曾硬写 float32）
    return rank, sv, N

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib", default="L1a", choices=["L1a", "L1b", "L2"])
    ap.add_argument("--ncol", type=int, default=0)
    ap.add_argument("--query", default="D3")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    arch = int(a.query.lstrip("D"))
    outdir = a.outdir or f"{OUT_BASE}/{a.lib}_{a.query}"
    os.makedirs(outdir, exist_ok=True)
    np.random.seed(SEED)

    up_sym, dn_sym = load_query(arch)
    log(f"查询 {a.query}：q_up {len(up_sym)} / q_down {len(dn_sym)}（符号）")

    if a.lib == "L2":
        import h5py
        h5 = h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r")
        genes = [x.decode() if isinstance(x, bytes) else x for x in h5["var/gene_name"][:]]
        gset = set(genes)
        up_ix = np.array([genes.index(g) for g in up_sym if g in gset])
        dn_ix = np.array([genes.index(g) for g in dn_sym if g in gset])
        ncol = h5["X"].shape[0] if a.ncol == 0 else min(a.ncol, h5["X"].shape[0])
        log(f"L2：{h5['X'].shape[0]} 扰动 × {len(genes)} 基因；命中 up {len(up_ix)} / down {len(dn_ix)}")
        X = np.asarray(h5["X"][:ncol, :], dtype=np.float32)
        lab = pd.DataFrame({
            "sig_id": [f"L2::{i}" for i in range(ncol)],
            "group": [str(d) for d in np.asarray(h5["obs/dataset"])[:ncol]],
            "pert_name": np.asarray(h5["obs/perturbed_gene_name"])[:ncol],
            "pert_sign": np.asarray(h5["obs/perturbation_sign"])[:ncol]})
    else:
        from cmapPy.pandasGEXpress.parse_gctx import parse
        gctx = f"{D_LINCS}/GSE70138_Level5.gctx" if a.lib == "L1a" else f"{D_LINCS}/GSE92742_Level5.gctx"
        if not os.path.exists(gctx): log(f"{gctx} 不存在（需先解压）"); sys.exit(2)
        G = parse(gctx); df = G.data_df
        genes = [str(g) for g in df.index]                      # Entrez
        gi = load_gene_map("GSE70138")
        sym2e = dict(zip(gi.pr_gene_symbol.astype(str), gi.pr_gene_id.astype(str)))
        gpos = {g: i for i, g in enumerate(genes)}
        up_ix = np.array([gpos[sym2e[g]] for g in up_sym if sym2e.get(g) in gpos])
        dn_ix = np.array([gpos[sym2e[g]] for g in dn_sym if sym2e.get(g) in gpos])
        ncol = df.shape[1] if a.ncol == 0 else min(a.ncol, df.shape[1])
        log(f"L1a：{len(genes)} 基因 × {df.shape[1]} 签名；取前 {ncol}；命中 up {len(up_ix)} / down {len(dn_ix)}")
        # gctx 自带元数据为空 ⇒ 用 sig_info 按 sig_id 映射
        series = "GSE70138" if a.lib == "L1a" else "GSE92742"
        si = load_sig_info(series).set_index("sig_id")
        sig_ids = np.array([str(c) for c in df.columns[:ncol]])
        meta = si.reindex(sig_ids)
        X = df.iloc[:, :ncol].values.T.astype(np.float32)        # ncol × genes
        lab = pd.DataFrame({"sig_id": sig_ids,
            "group": (meta.cell_id.fillna("NA").astype(str) + "|" + meta.pert_type.fillna("NA").astype(str)).values,
            "pert_name": meta.pert_iname.values, "pert_type": meta.pert_type.values,
            "cell_id": meta.cell_id.values})
        del df, G

    log("预算秩矩阵（order/rank/sorted_vals）…")
    t0 = time.time()
    rank, sv, N = build_rank_tables(X)
    log(f"  秩矩阵完成 {time.time()-t0:.0f}s；形状 {rank.shape}")
    del X
    nsig = rank.shape[0]

    def ncs_of_sub(w, idx):
        """在子集 idx 上按组做 NCS（与 ncs_of 同式，只作用于子集）。"""
        g_all = lab.group.values[idx]
        out = np.empty_like(w)
        for g in np.unique(g_all):
            m = g_all == g
            ww = w[m]
            mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
            mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
            out[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
        return out

    def ncs_of(w):
        out = np.empty_like(w)
        for g in np.unique(lab.group.values):
            m = lab.group.values == g
            ww = w[m]
            # 🔴 μ⁻ 必须取**负侧均值的绝对值**：CMap 原文 μ⁻ 是"negative 侧的 scale"，
            #    若直接用带符号均值（负数），w/μ⁻ 会把符号翻正 ⇒ NCS 全正、τ 全正、
            #    "τ≤−90" 恒为空（我踩过，表面像"没有东西能逆转"）
            mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
            mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
            out[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
        return out

    log("主查询 WTCS …")
    # 双臂：up 与 down 各自 ES 再合成
    def wtcs_pair(q_up, q_dn):
        Ru,Cu = _rc(rank, sv, q_up); eup = es_from_ranks(Ru, Cu, N, len(q_up))
        Rd,Cd = _rc(rank, sv, q_dn); edn = es_from_ranks(Rd, Cd, N, len(q_dn))
        return np.where(np.sign(eup) != np.sign(edn), (eup - edn) / 2.0, 0.0)
    w_main = wtcs_pair(up_ix, dn_ix)
    ncs_main = ncs_of(w_main)
    log(f"  主查询完成；WTCS 范围 [{w_main.min():.3f}, {w_main.max():.3f}]")

    # ───────── §5.14 Q_ref：Touchstone exemplar（三档） ─────────
    DATA = "/home/eto/luad_v2/results/10_niche/target_reversal/data"
    sm = pd.read_csv(f"{DATA}/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",
                     sep="\t", compression="gzip").set_index("sig_id").reindex(sig_ids)
    pi2 = pd.read_csv(f"{DATA}/lincs_pert_info2.tsv", sep="\t").set_index("pert_id")
    lab = lab.reset_index(drop=True)
    lab["pert_id"] = meta.pert_id.values
    lab["is_ts"] = pi2.is_touchstone.reindex(lab.pert_id.values).values
    lab["distil_ss"] = sm.distil_ss.values
    lab["distil_cc_q75"] = sm.distil_cc_q75.values
    lab["distil_nsample"] = sm.distil_nsample.values
    log("touchstone 命中(非NA) %d/%d ; nsample>=3 的 %d" %
        (lab.is_ts.notna().sum(), len(lab), (lab.distil_nsample >= 3).sum()))

    def build_qref(tier):
        d = lab[(lab.is_ts == 1) & (lab.distil_nsample >= 3)].copy()
        if tier == "all":
            return set(d.sig_id)
        col = "distil_ss" if tier == "ss" else "distil_cc_q75"
        d = d.sort_values([col, "sig_id"], ascending=[False, True])
        return set(d.drop_duplicates(["cell_id", "pert_id"]).sig_id)

    _qmask_sub = None   # 由 QC 段设置（子集用的整数掩码）
    def tau_of_full(ncs_sub, idx_sub, qref):
        """在子集 idx_sub 上算 τ（QC 用，避免全量开销）。"""
        sig_all = lab.sig_id.values; grp_all = lab.group.values
        out = np.zeros(len(idx_sub))
        for g in np.unique(grp_all[idx_sub]):
            loc = np.where(grp_all[idx_sub] == g)[0]
            glob = idx_sub[loc]
            # 🔴 参照必须是**该查询自己的** NCS（我曾误用 ncs_main ⇒ QC 恒 0）
            ref = np.abs(ncs_sub[loc][_qmask_sub[glob]])
            if ref.size == 0: continue
            ref = np.sort(ref)
            cnt = np.searchsorted(ref, np.abs(ncs_sub[loc]), side="left")
            out[loc] = np.sign(ncs_sub[loc]) * 100.0 / ref.size * cnt
        return out

    # 🔴 `np.isin` 对**字符串**数组走慢路径（O(n×m)），6,059 条就 16 分钟/档，
    #    "all" 档（108k）会到 ~4 小时。改用**整数索引掩码**，降为 O(n)。
    SIG = lab.sig_id.values
    GRP = lab.group.values
    _IDX = {s: i for i, s in enumerate(SIG)}
    def qmask_of(qref):
        m = np.zeros(len(SIG), dtype=bool)
        for s in qref:
            j = _IDX.get(s)
            if j is not None:
                m[j] = True
        return m

    def tau_of(ncs, qref):
        if not isinstance(qref, np.ndarray):
            qref = qmask_of(qref)
        out = np.zeros_like(ncs)
        for g in np.unique(GRP):
            m = GRP == g
            ref = np.abs(ncs[m & qref])
            if ref.size == 0:
                continue
            ref = np.sort(ref)
            cnt = np.searchsorted(ref, np.abs(ncs[m]), side="left")
            out[m] = np.sign(ncs[m]) * 100.0 / ref.size * cnt
        return out

    TAUS = {}
    for tier in ["ss", "cc", "all"]:
        qr = build_qref(tier)
        TAUS[tier] = tau_of(ncs_main, qr)
        log("  档 %-3s : Q_ref %d 个签名 ; τ≤−90 %d ; τ≤−95 %d ; τ≤−98 %d" %
            (tier, len(qr), (TAUS[tier] <= -90).sum(), (TAUS[tier] <= -95).sum(), (TAUS[tier] <= -98).sum()))

    # ── §5.6 修订：零模型改为**小规模 QC** ──
    # 🔴 2026-10-04 修订：原设计跑 1000 次全量 ES（全任务约 20 小时，不可行），
    #    且**在概念上多余**——τ 本身就是"|NCS| 在 Touchstone 参照集里的百分位"，
    #    零分布**内建在 τ 中**（官方用 τ≥90 当阈值 = 10% 尾）。故改为 QC：
    #    随机查询下 τ 应近似均匀；偏了就说明 τ 构造有偏。
    N_QC = 100
    rng = np.random.default_rng(SEED)
    SUBS = np.random.default_rng(SEED).choice(nsig, size=min(5000, nsig), replace=False)
    log(f"零模型 QC：{N_QC} 个随机查询 × {len(SUBS)} 个签名（验证 τ 近似均匀）")
    rank_s, sv_s = rank[SUBS], sv[SUBS]          # 只在子集上算，省 24 倍
    # 🔴 QC 的正确层次：随机查询下，**全部 (查询×签名) 对**的 τ 应近似均匀
    #    ⇒ 中位 ≈ 0，且 τ≤−90 的占比 ≈ 10%。曾误把"每个查询的 τ 中位数"当统计量（那是 0 才对）
    qref_ss = build_qref("ss"); all_tau = []
    _mq = qmask_of(qref_ss)
    _qmask_sub = _mq   # 全局给 tau_of_full 用（子集内的全局下标）
    for k in range(N_QC):
        g = rng.choice(N, size=300, replace=False)
        Ru, Cu = _rc(rank_s, sv_s, g[:150]); e_u = es_from_ranks(Ru, Cu, N, 150)
        Rd, Cd = _rc(rank_s, sv_s, g[150:]); e_d = es_from_ranks(Rd, Cd, N, 150)
        ww = np.where(np.sign(e_u) != np.sign(e_d), (e_u - e_d) / 2.0, 0.0)
        all_tau.append(tau_of_full(ncs_of_sub(ww, SUBS), SUBS, qref_ss))
    at = np.concatenate(all_tau)
    log("  QC：随机查询下 τ 中位 %.1f（期望≈0）；τ≤−90 占比 %.1f%%、τ≥+90 占比 %.1f%%（各期望≈10%%）；n=%d" %
        (np.median(at), 100 * np.mean(at <= -90), 100 * np.mean(at >= 90), at.size))

    tau = TAUS["ss"]
    tau = TAUS["ss"]
    # ── selectivity（CMap 2017 §"Selectivity of PCL Connections" 的精神）──
    # 原文按 PCL（perturbagen class）算"query 连不上多少比例"。此处按**化合物聚合**：
    # 对每个化合物取其所有签名的 |NCS| 最大值，再看该 query 的 τ 分布里，
    # 有多少比例的化合物"连不上"（|τ| < 90）。s 越大 = 越选择性（不是对所有化合物都强）。
    def selectivity(tau_vec):
        agg = pd.Series(np.abs(tau_vec)).groupby(lab.pert_name.values).max()
        return float((agg < 90).mean())
    SELECT = {t: selectivity(TAUS[t]) for t in TAUS}
    log("  选择性 s（|τ|<90 的化合物占比，越接近 1 越选择性）：%s" %
        " ; ".join(f"{t}={SELECT[t]:.3f}" for t in TAUS))

    res = lab.copy(); res["wtcs"] = w_main; res["ncs"] = ncs_main
    for tier in ["ss", "cc", "all"]:
        res["tau_" + tier] = TAUS[tier]
    res["n_q_up"] = len(up_ix); res["n_q_down"] = len(dn_ix)
    res = res.sort_values("tau_ss")
    res.to_csv(f"{outdir}/ranked.tsv", sep="\t", index=False)
    log(f"落盘 {outdir}/ranked.tsv（{len(res)} 行）")
    log("τ≤−90: %d ；≤−95: %d ；≤−98: %d" % ((tau<=-90).sum(), (tau<=-95).sum(), (tau<=-98).sum()))
    print(res.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
