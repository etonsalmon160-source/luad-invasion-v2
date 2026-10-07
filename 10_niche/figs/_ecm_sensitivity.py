#!/usr/bin/env python3
# _ecm_sensitivity.py —— ECM 程序结论对手写基因的敏感度（审计用）
#   22 个基因里只有 14 个落在 D3 深度配平 top150 ⇒ 检验结论是否由手选的 8 个驱动
import numpy as np, pandas as pd, h5py
from scipy import sparse
ROOT = "/home/eto/luad_v2"; SC = f"{ROOT}/results/10_niche/tr_singlecell"
ALL22 = ["COL1A1","COL1A2","COL3A1","COL6A3","FN1","CTHRC1","POSTN","SPARC","THBS2","BGN",
         "COL5A1","COL5A2","LUM","VCAN","MMP2","TIMP1","ASPH","COL8A1","COMP","THY1","COL10A1","COL11A1"]
d = pd.read_csv(f"{ROOT}/results/10_niche/target_reversal/d3_depthmatched_lfc.tsv", sep="\t")
d.columns = ["gene", "lfc"]
d = d.sort_values("lfc", ascending=False).reset_index(drop=True); d["rk"] = d.index + 1
RK = dict(zip(d.gene, d.rk)); TOT = len(d)
DERIVED = [g for g in ALL22 if g in RK and RK[g] <= 150]
HAND    = [g for g in ALL22 if g not in DERIVED]
print(f"数据衍生 {len(DERIVED)}: {DERIVED}")
print(f"手工加入 {len(HAND)}: {HAND}")

D = pd.read_parquet(f"{SC}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as h:
    g = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["genes"][:]])
    M = sparse.csr_matrix((h["counts/data"][:], h["counts/indices"][:], h["counts/indptr"][:]),
                          shape=(len(D), len(g)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
pid = D.patient_id.values; l2 = D.L2.values
PRE = np.isin(D.stage.values, ["AAH", "AIS", "MIA"]); IAC = D.stage.values == "IAC"
fib = (l2 == "Fibroblast")
mu = np.asarray(Mn[fib].mean(0)).ravel()
sd = np.sqrt(np.maximum(np.asarray(Mn[fib].power(2).mean(0)).ravel() - mu ** 2, 1e-12))
mean_fib = mu.copy()
ZT = {}
for p in np.unique(pid):
    for tag, m in [("I", (pid == p) & fib & IAC), ("P", (pid == p) & fib & PRE)]:
        if m.sum() >= 20:
            ZT[(p, tag)] = (np.asarray(Mn[m].mean(0)).ravel() - mu) / sd
VALID = [p for p in np.unique(pid) if (p, "I") in ZT and (p, "P") in ZT]
Z = {p: (ZT[(p, "I")], ZT[(p, "P")]) for p in VALID}
del M, Mn
gi = {x: i for i, x in enumerate(g)}
print(f"可用配对患者 {len(VALID)} 例\n")

def stat(names):
    idx = np.array([gi[x] for x in names if x in gi])
    dd = np.array([Z[p][0][idx].mean() - Z[p][1][idx].mean() for p in VALID])
    return int((dd > 0).sum()), float(np.median(dd)), len(idx)

rng = np.random.default_rng(20261006); NP = 300
cand = np.where(mean_fib > np.median(mean_fib))[0]
print(f"{'基因集':22s} {'n':>3s} {'例数':>8s} {'中位差':>10s} {'P(例数)':>9s} {'P(中位差)':>10s}")
for nm, lst in [("全部 22（现行）", ALL22), ("仅数据衍生 14", DERIVED), ("仅手工加入 8", HAND)]:
    w0, d0, n = stat(lst)
    W = np.empty(NP, int); Ds = np.empty(NP)
    for k in range(NP):
        ix = rng.choice(cand, size=n, replace=False)
        W[k], Ds[k], _ = stat([g[i] for i in ix])
    pw = (1 + (W >= w0).sum()) / (1 + NP); pd_ = (1 + (Ds >= d0).sum()) / (1 + NP)
    print(f"{nm:22s} {n:3d} {w0:>4d}/{len(VALID):<3d} {d0:>+10.4f} {pw:>9.4f} {pd_:>10.4f}")

# 逐个剔除看稳定性
print("\n留一法（去掉某一个基因后，中位差如何变）：")
base = stat(ALL22)[1]
for gname in ALL22:
    rest = [x for x in ALL22 if x != gname]
    w1, d1, _ = stat(rest)
    flag = "  <-- 影响较大" if abs(d1 - base) > 0.15 * abs(base) else ""
    print(f"  去掉 {gname:9s}（D3 秩 {RK.get(gname,0):>5d}）: {w1}/{len(VALID)} 例, 中位差 {d1:+.4f}{flag}")
