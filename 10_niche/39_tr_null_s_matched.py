#!/usr/bin/env python3
# 39_tr_null_s_matched.py —— 诊断：**S（命中集大小）配平的零模型**
#
# 用户的质疑：「你确定我们的选药程序完备，好像选药程序本身就是降级的随机因素」
#
# 审计发现的两个口径不一致：
#   ① 真实查询的 150+150 个基因**只有一部分能映射进 LINCS 库**
#      （D3 实测 up 113 / down 76），而 28_/37_ 的零模型 `rng.choice(N,300)`
#      从**库内**抽基因 ⇒ **必然 300/300 全命中**。两边 S 不同 ⇒ ES 分布形状不同。
#   ② τ 的参照集 Q_ref 是 exemplar（高 distil_ss）⇒ τ 分布零膨胀（中位恒为 0）
#      ⇒ "低于 −90 的个数"是零膨胀分布上的尾部计数，方差大、功效低。
#
# 本脚本只做 ①：把零模型的抽样池换成**我们自己的 Visium 基因池**，
# 使随机查询的映射率与真实查询一致（S 自动配平），**其余口径逐字不变**。
#
# 用法：python3 39_tr_null_s_matched.py [N_RAND]
import os, sys, time, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np, pandas as pd
import multiprocessing as mp
from importlib.machinery import SourceFileLoader

HERE = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche"
TR = SourceFileLoader("tr", f"{HERE}/24_target_reversal.py").load_module()

ROOT = "/home/eto/luad_v2"
OUT = f"{ROOT}/results/10_niche/target_reversal"
DATA = f"{OUT}/data"
KD = f"{ROOT}/results/10_niche/kstar_diag"
VIS = f"{ROOT}/data/visium_spatial"
N_RAND = int(sys.argv[1]) if len(sys.argv) > 1 else 100
SEED = 20261005
N_WORKER = 18

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ══ 库（与 28_/37_ 同） ══
from cmapPy.pandasGEXpress.parse_gctx import parse
G = parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df = G.data_df
genes = [str(g) for g in df.index]
X = df.values.T.astype(np.float32)
sig_ids = np.array([str(c) for c in df.columns]); N = X.shape[1]
del df, G
log(f"库：{X.shape[0]} 签名 × {N} 基因")
rank, sv, _ = TR.build_rank_tables(X); del X
gpos = {g: i for i, g in enumerate(genes)}

gi = TR.load_gene_map("GSE70138")
sym2e = dict(zip(gi.pr_gene_symbol.astype(str), gi.pr_gene_id.astype(str)))

si = TR.load_sig_info("GSE70138").set_index("sig_id").reindex(sig_ids)
sm = pd.read_csv(f"{DATA}/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",
                 sep="\t", compression="gzip").set_index("sig_id").reindex(sig_ids)
pi2 = pd.read_csv(f"{DATA}/lincs_pert_info2.tsv", sep="\t").set_index("pert_id")
lab = pd.DataFrame({"sig_id": sig_ids,
                    "cell_id": si.cell_id.fillna("NA").values,
                    "pert_type": si.pert_type.fillna("NA").values,
                    "pert_id": si.pert_id.values,
                    "pert_name": si.pert_iname.values})
lab["is_ts"] = pi2.is_touchstone.reindex(lab.pert_id.values).values
lab["group"] = lab.cell_id + "|" + lab.pert_type
lab["distil_ss"] = sm.distil_ss.values
lab["distil_nsample"] = sm.distil_nsample.values
qref = set(lab[(lab.is_ts == 1) & (lab.distil_nsample >= 3)]
           .sort_values(["distil_ss", "sig_id"], ascending=[False, True])
           .drop_duplicates(["cell_id", "pert_id"]).sig_id)
IDX = {s: i for i, s in enumerate(sig_ids)}
qmask = np.zeros(len(sig_ids), dtype=bool)
for s in qref: qmask[IDX[s]] = True
GRP = lab.group.values; UGRP = np.unique(GRP)
GRPM = [GRP == g for g in UGRP]
log(f"Q_ref {len(qref)} 条；分组 {len(UGRP)} 个")

# ══ ⭐ 我们自己的 Visium 基因池（零模型的抽样池）══
# 🔴 抽样池必须是**未映射的完整符号池**；映射在抽完之后做，否则 300/300 必然全命中。
feat = np.array(sorted({p for f in glob.glob(f"{VIS}/*/filtered_feature_bc_matrix/features.tsv.gz")
                          for p in pd.read_csv(f, sep="\t", header=None)[1].astype(str)}))
POOL = np.array([gpos.get(sym2e.get(g)) if sym2e.get(g) in gpos else -1 for g in feat], dtype=np.int64)
log(f"Visium 符号池 {len(feat)} ⇒ 其中可映射进库 {int((POOL >= 0).sum())}"
    f"（{100*(POOL>=0).mean():.1f}%）")

# 真实查询的映射率（用于对照）
d3 = pd.read_csv(f"{OUT}/L1a_D3/ranked.tsv", sep="\t").iloc[0]
log(f"真实 D3 查询映射率：up {int(d3.n_q_up)}/150 ；down {int(d3.n_q_down)}/150")

def ncs_of(w, _m=GRPM):
    o = np.empty_like(w)
    for m in _m:
        ww = w[m]
        mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
        mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
        o[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
    return o

def tau_of(ncs, _m=GRPM):
    o = np.zeros_like(ncs)
    for m in _m:
        ref = np.sort(np.abs(ncs[m & qmask]))
        if ref.size == 0: continue
        o[m] = np.sign(ncs[m]) * 100.0 / ref.size * np.searchsorted(ref, np.abs(ncs[m]), side="left")
    return o

_R, _S, _POOL = rank, sv, POOL

def one_rep(k):
    rng = np.random.default_rng(SEED + 7919 * k)
    sel = rng.choice(len(_POOL), size=300, replace=False)   # 从**符号池**抽
    ix = _POOL[sel]                                          # 抽完再映射
    gu = np.unique(ix[:150][ix[:150] >= 0])                  # 去重（符号别名可能撞同一 Entrez）
    gd = np.unique(ix[150:][ix[150:] >= 0])
    Ru, Cu = TR._rc(_R, _S, gu); e_u = TR.es_from_ranks(Ru, Cu, N, len(gu))
    Rd, Cd = TR._rc(_R, _S, gd); e_d = TR.es_from_ranks(Rd, Cd, N, len(gd))
    w = np.where(np.sign(e_u) != np.sign(e_d), (e_u - e_d) / 2.0, 0.0)
    t = tau_of(ncs_of(w))
    return {"rep": k,
            "s_up": len(gu), "s_dn": len(gd),
            "sig_le_m90": float((t <= -90).mean()),
            "sig_le_m95": float((t <= -95).mean()),
            "mean_tau": float(t.mean()),
            "q05_tau": float(np.quantile(t, 0.05)),
            "min_tau": float(t.min())}

log(f"S 配平零模型：{N_RAND} 次，{N_WORKER} 进程（抽样池 = Visium 基因池）")
out_f = f"{OUT}/null_smatched.tsv"
with open(out_f, "w") as fh:
    fh.write("rep\ts_up\ts_dn\tsig_le_m90\tsig_le_m95\tmean_tau\tq05_tau\tmin_tau\n"); fh.flush()
    t0 = time.time(); done = 0
    with mp.Pool(N_WORKER) as pool:
        for r in pool.imap_unordered(one_rep, range(N_RAND)):
            fh.write("\t".join(str(r[c]) for c in
                     ["rep", "s_up", "s_dn", "sig_le_m90", "sig_le_m95", "mean_tau", "q05_tau", "min_tau"]) + "\n")
            fh.flush(); done += 1
            if done % 20 == 0 or done == N_RAND:
                el = time.time() - t0
                log(f"   {done}/{N_RAND}  已用 {el/60:.1f} min  预计总 {el/done*N_RAND/60:.1f} min")

# ══ 汇总：与旧零模型并排 ══
R = pd.read_csv(out_f, sep="\t")
OLD = pd.read_csv(f"{OUT}/null_pcl.tsv", sep="\t")
log("=" * 64)
log(f"S 配平零模型实际 S：up {R.s_up.mean():.0f}（范围 {R.s_up.min()}–{R.s_up.max()}）"
    f" ；down {R.s_dn.mean():.0f}（{R.s_dn.min()}–{R.s_dn.max()}）")
log(f"  vs 旧零模型 S = 150 / 150  ；真实 D3 = {int(d3.n_q_up)} / {int(d3.n_q_down)}")
log("-" * 64)
log("统计量                     旧零模型(库内抽样)      S配平零模型(Visium池)")
for col, name in [("sig_le_m90", "签名层 τ≤−90 占比")]:
    a, b = OLD[col].values, R[col].values
    log(f"  {name:24s} {a.mean():.4f} ± {a.std(ddof=1):.4f}   {b.mean():.4f} ± {b.std(ddof=1):.4f}")
log("-" * 64)
OBS = {"S1 原签名": 0.0660, "S2 深度配平": 0.0447, "纯深度对照": 0.0423}
for tag, ob in OBS.items():
    a, b = OLD.sig_le_m90.values, R.sig_le_m90.values
    pa = (1 + (a >= ob).sum()) / (1 + len(a)); pb = (1 + (b >= ob).sum()) / (1 + len(b))
    za = (ob - a.mean()) / a.std(ddof=1); zb = (ob - b.mean()) / b.std(ddof=1)
    log(f"  {tag:12s} 观测 {ob:.4f}   旧: z={za:+.2f} p={pa:.3f}   S配平: z={zb:+.2f} p={pb:.3f}")
log("判据：若 S 配平后 z 明显下降 ⇒ 『S 口径不一致』是原结论的一部分成因，须登记为程序缺陷。")
log("完成")
