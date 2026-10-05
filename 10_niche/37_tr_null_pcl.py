#!/usr/bin/env python3
# 37_tr_null_pcl.py —— 靶点扰动逆向臂：**PCL 层面的零模型**（补 §TR-2 的缺口）
#
# 问题：S1 有 11 个 MOA 类 τ_PCL≤−90、S2 有 18 个。但签名层面的零模型（28_）证明
#       S2 的签名层面落随机内。⇒ 那 11/18 到底是不是"类内聚合放大的假象"？
#       必须**用同一个统计量**（#PCL τ_PCL≤−90）在随机查询上标定。
#
# 做法：与 28_ 完全同一条管道，只把报告统计量从"签名占比"换成"PCL 命中数"。
#       τ_ss = τ 已有；PCL 聚合与 τ_PCL 百分位**逐字复刻 29_tr_selectivity.py**。
#
# 步骤：
#   ① 自检 —— 对已落盘的 L1a_D3 / L1a_D3_depthmatched 的 ncs 列重算 PCL，
#             必须分别复现 11 / 18。对不上就中止（说明我的复刻有误）。
#   ② 零模型 —— N_RAND 个随机 300 基因查询走全管道 ⇒ 每次的 PCL 命中数分布。
#   ③ 给出观测值的经验 p。
#
# 用法：python3 37_tr_null_pcl.py [N_RAND]   （默认 200）
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np, pandas as pd
import multiprocessing as mp
from importlib.machinery import SourceFileLoader

HERE = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche"
TR = SourceFileLoader("tr", f"{HERE}/24_target_reversal.py").load_module()

OUT = "/home/eto/luad_v2/results/10_niche/target_reversal"
DATA = f"{OUT}/data"
N_RAND = int(sys.argv[1]) if len(sys.argv) > 1 else 200
SEED = 20261005
N_WORKER = 18

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ══════════════════════════════════════════════════════════════════
# 建库（与 28_ 逐字同口径）
# ══════════════════════════════════════════════════════════════════
from cmapPy.pandasGEXpress.parse_gctx import parse
G = parse("/home/eto/lincs_data/GSE70138_Level5.gctx"); df = G.data_df
genes = [str(g) for g in df.index]
X = df.values.T.astype(np.float32)
sig_ids = np.array([str(c) for c in df.columns]); N = X.shape[1]
del df, G
log(f"库：{X.shape[0]} 签名 × {N} 基因")
rank, sv, _ = TR.build_rank_tables(X); del X

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
lab["lib_idx"] = np.arange(len(lab))

qref = set(lab[(lab.is_ts == 1) & (lab.distil_nsample >= 3)]
           .sort_values(["distil_ss", "sig_id"], ascending=[False, True])
           .drop_duplicates(["cell_id", "pert_id"]).sig_id)
IDX = {s: i for i, s in enumerate(sig_ids)}
qmask = np.zeros(len(sig_ids), dtype=bool)
for s in qref: qmask[IDX[s]] = True
GRP = lab.group.values
UGRP = np.unique(GRP)
log(f"Q_ref {len(qref)} 条；分组 {len(UGRP)} 个")

# ── MOA → 库内签名索引（与 29_ 同：先按 pert_name 小写 merge）──
mo = pd.read_csv(f"{DATA}/clue_moa_list.tsv", sep="\t")
mlab = lab[["lib_idx", "pert_name"]].copy()
mlab["pn"] = mlab.pert_name.astype(str).str.lower()
M = mlab[["lib_idx", "pn"]].merge(mo, left_on="pn", right_on="pert_name", how="inner")
MOA_NAMES = sorted(M.moa.unique())
_mid = {m: i for i, m in enumerate(MOA_NAMES)}
M["c"] = M.moa.map(_mid).astype(int)
IDX_BY_C = [M.lib_idx.values[M.c.values == c] for c in range(len(MOA_NAMES))]
N_PCL = len(MOA_NAMES)
log(f"MOA 类 {N_PCL} 个；可映射签名 {len(M)} 条（{100*len(M)/len(lab):.1f}%）")

GRPM = [GRP == g for g in UGRP]   # 预切分组掩码

def ncs_of(w, _m=GRPM):
    o = np.empty_like(w)
    for m in _m:
        ww = w[m]
        mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
        mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
        o[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
    return o

def pcl_tau(ncs, _idx=IDX_BY_C, _n=N_PCL):
    """逐字复刻 29_tr_selectivity.py 的 PCL 聚合 + τ_PCL 百分位定义。"""
    npc = np.empty(_n)
    for c, ix in enumerate(_idx):
        v = ncs[ix]
        q67 = np.quantile(v, 0.67); q33 = np.quantile(v, 0.33)
        npc[c] = q67 if abs(q67) >= abs(q33) else q33
    a = np.sort(np.abs(npc))
    return np.sign(npc) * 100.0 / len(npc) * np.searchsorted(a, np.abs(npc), side="left")

# ══════════════════════════════════════════════════════════════════
# ① 自检：重算已落盘的 S1 / S2，必须复现 11 / 18
# ══════════════════════════════════════════════════════════════════
log("① 自检：用落盘的 ncs 列重算 PCL")
# (文件, 该文件里 NCS 的列名) —— 两版落盘列名不同
CHECK = {"S1": (f"{OUT}/L1a_D3/ranked.tsv", "ncs"),
         "S2": (f"{OUT}/L1a_D3_depthmatched_ranked.tsv", "ncs_dm")}
EXPECT = {"S1": 11, "S2": 18}
ok = True
for tag, (f, ncs_col) in CHECK.items():
    d = pd.read_csv(f, sep="\t")
    # 按 sig_id 对齐到库顺序
    d = d.set_index("sig_id").reindex(sig_ids)
    if d[ncs_col].isna().any():
        log(f"  {tag}: 有 {int(d[ncs_col].isna().sum())} 条对不上库 ⇒ 中止"); ok = False; break
    t = pcl_tau(d[ncs_col].values.astype(np.float64))
    hit = int((t <= -90).sum())
    mark = "✅" if hit == EXPECT[tag] else "❌"
    log(f"  {tag}: PCL 命中 {hit} / {N_PCL}（期望 {EXPECT[tag]}）{mark}")
    if hit != EXPECT[tag]: ok = False
if not ok:
    log("自检不过 ⇒ 中止（先修脚本，不要跑零模型）"); sys.exit(2)
log("自检通过 ✅  ⇒ 进入零模型")

# ══════════════════════════════════════════════════════════════════
# ② 零模型：每个随机查询走全管道
# ══════════════════════════════════════════════════════════════════
_R = rank; _S = sv

def one_rep(k):
    g = np.random.default_rng(SEED + k).choice(N, size=300, replace=False)
    Ru, Cu = TR._rc(_R, _S, g[:150]); e_u = TR.es_from_ranks(Ru, Cu, N, 150)
    Rd, Cd = TR._rc(_R, _S, g[150:]); e_d = TR.es_from_ranks(Rd, Cd, N, 150)
    w = np.where(np.sign(e_u) != np.sign(e_d), (e_u - e_d) / 2.0, 0.0)
    ncs = ncs_of(w)
    tau_ss = np.empty_like(ncs)
    for m in GRPM:
        ref = np.sort(np.abs(ncs[m & qmask]))
        if ref.size == 0: continue
        tau_ss[m] = np.sign(ncs[m]) * 100.0 / ref.size * np.searchsorted(ref, np.abs(ncs[m]), side="left")
    tp = pcl_tau(ncs)
    return {"rep": k,
            "sig_le_m90": float((tau_ss <= -90).mean()),
            "n_pcl_le_m90": int((tp <= -90).sum()),
            "n_pcl_le_m95": int((tp <= -95).sum()),
            "n_pcl_le_m98": int((tp <= -98).sum()),
            "min_tau_pcl": float(tp.min())}

log(f"② 零模型：{N_RAND} 个随机查询，{N_WORKER} 进程")
out_f = f"{OUT}/null_pcl.tsv"
with open(out_f, "w") as fh:
    fh.write("rep\tsig_le_m90\tn_pcl_le_m90\tn_pcl_le_m95\tn_pcl_le_m98\tmin_tau_pcl\n"); fh.flush()
    t0 = time.time(); done = 0
    with mp.Pool(N_WORKER) as pool:
        for r in pool.imap_unordered(one_rep, range(N_RAND)):
            fh.write("\t".join(str(r[c]) for c in
                     ["rep", "sig_le_m90", "n_pcl_le_m90", "n_pcl_le_m95", "n_pcl_le_m98", "min_tau_pcl"]) + "\n")
            fh.flush(); done += 1
            if done % 10 == 0 or done == N_RAND:
                el = time.time() - t0
                log(f"   {done}/{N_RAND}  已用 {el/60:.1f} min  预计总 {el/done*N_RAND/60:.1f} min")

# ══════════════════════════════════════════════════════════════════
# ③ 经验 p
# ══════════════════════════════════════════════════════════════════
R = pd.read_csv(out_f, sep="\t")
log("=" * 62)
log(f"零模型 PCL 命中数（τ_PCL≤−90）：均值 {R.n_pcl_le_m90.mean():.2f}  "
    f"范围 [{R.n_pcl_le_m90.min()}, {R.n_pcl_le_m90.max()}]  中位 {R.n_pcl_le_m90.median():.0f}")
for thr, col, obs1, obs2 in [(-90, "n_pcl_le_m90", 11, 18),
                             (-95, "n_pcl_le_m95", None, None),
                             (-98, "n_pcl_le_m98", None, None)]:
    v = R[col].values
    if obs1 is not None:
        p1 = (1 + (v >= obs1).sum()) / (1 + len(v))
        p2 = (1 + (v >= obs2).sum()) / (1 + len(v))
        log(f"  τ_PCL≤{thr}: 零模型 均值 {v.mean():.2f} / q95 {np.quantile(v,.95):.0f} / max {v.max()}"
            f"  ⇒  S1(11) 经验 p = {p1:.3f} ; S2(18) 经验 p = {p2:.3f}")
    else:
        log(f"  τ_PCL≤{thr}: 零模型 均值 {v.mean():.2f} / q95 {np.quantile(v,.95):.0f} / max {v.max()}")
log("判据：经验 p ≥ 0.05 ⇒ 该 MOA 层结果与随机查询不可区分，不得作候选。")
log("完成")
