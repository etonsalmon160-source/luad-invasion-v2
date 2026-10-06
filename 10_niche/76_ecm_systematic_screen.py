#!/usr/bin/env python3
# 76_ecm_systematic_screen.py —— D-5a：ECM 程序到底特不特殊（系统筛查）
#
# 用户批评：「这个程序不是很正常的程序吗，其实还有其他类似的系统化程序只是你没去筛」
# 做法：把 ECM 那 22 个基因当"一个程序"，用**随机基因程序**建零分布，看它排第几。
#   N1 全基因池随机 22 基因
#   N2 按**成纤维平均表达分位**配平的随机 22 基因（排除"只是表达高"）
#   N3 从**成纤维表达基因池**随机 22 基因（排除"只是成纤维表达"）
# 统计量：22 例里 IAC>前驱 的例数 + 患者内配对中位差
import numpy as np, pandas as pd, h5py, time
from scipy import sparse
ROOT = "/home/eto/luad_v2"; SC = f"{ROOT}/results/10_niche/tr_singlecell"
ECM = ["COL1A1","COL1A2","COL3A1","COL6A3","FN1","CTHRC1","POSTN","SPARC","THBS2","BGN",
       "COL5A1","COL5A2","LUM","VCAN","MMP2","TIMP1","ASPH","COL8A1","COMP","THY1","COL10A1","COL11A1"]
N_PERM = 300; SEED = 20261006
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

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
log(f"成纤维 {int(fib.sum()):,} 细胞 / 患者 {len(set(pid))}")

# ⚡ 预计算：逐患者×逐组的 z 均值 ⇒ 之后任何程序都是查表
mu = np.asarray(Mn[fib].mean(0)).ravel()
sd = np.sqrt(np.maximum(np.asarray(Mn[fib].power(2).mean(0)).ravel() - mu ** 2, 1e-12))
mean_fib = mu.copy()
ZT = {}
for p in np.unique(pid):
    for tag, msk in [("I", (pid == p) & fib & IAC), ("P", (pid == p) & fib & PRE)]:
        if msk.sum() >= 20:
            ZT[(p, tag)] = (np.asarray(Mn[msk].mean(0)).ravel() - mu) / sd
VALID = [p for p in np.unique(pid) if (p, "I") in ZT and (p, "P") in ZT]
Z = {p: (ZT[(p, "I")], ZT[(p, "P")]) for p in VALID}
log(f"可用于配对的患者 {len(VALID)} 例")
del M, Mn

gi = {x: i for i, x in enumerate(g)}
ecm = np.array([gi[x] for x in ECM if x in gi])

def stat(idx):
    d = np.array([Z[p][0][idx].mean() - Z[p][1][idx].mean() for p in VALID])
    return int((d > 0).sum()), float(np.median(d))

w0, d0 = stat(ecm)
log(f"★ ECM 程序（{len(ecm)} 基因）：{w0}/{len(VALID)} 例 IAC>前驱，中位差 {d0:+.4f}")

rng = np.random.default_rng(SEED)
nb = 20
binid = pd.qcut(pd.Series(mean_fib).rank(method="first"), nb, labels=False).values
BIN_IDX = [np.where(binid == b)[0] for b in range(nb)]
ecm_bins = [binid[i] for i in ecm]
cand_fib = np.where(mean_fib > np.median(mean_fib))[0]
log(f"分位池 {[len(x) for x in BIN_IDX]}；成纤维表达池 {len(cand_fib)}")

for nm, kind in [("N1 全基因池", 0), ("N2 表达分位配平", 1), ("N3 成纤维表达池", 2)]:
    W = np.empty(N_PERM, dtype=int); Ds = np.empty(N_PERM)
    for k in range(N_PERM):
        if kind == 0:
            idx = rng.choice(len(g), size=len(ecm), replace=False)
        elif kind == 1:
            idx = np.array([BIN_IDX[b][rng.integers(len(BIN_IDX[b]))] for b in ecm_bins])
        else:
            idx = rng.choice(cand_fib, size=len(ecm), replace=False)
        W[k], Ds[k] = stat(idx)
    p_w = (1 + (W >= w0).sum()) / (1 + N_PERM); p_d = (1 + (Ds >= d0).sum()) / (1 + N_PERM)
    print(f"\n{nm}（{N_PERM} 次；ECM = {w0}/{len(VALID)} 例，中位差 {d0:+.4f}）")
    print(f"  例数  ：均值 {W.mean():.1f}，最大 {W.max()}  ⇒ P(≥{w0}) = {p_w:.4f}")
    print(f"  中位差：均值 {Ds.mean():+.4f}，最大 {Ds.max():+.4f}  ⇒ P(≥{d0:+.4f}) = {p_d:.4f}")
    print(f"  ECM 的百分位：例数 {100*(W<w0).mean():.1f}%，中位差 {100*(Ds<d0).mean():.1f}%")
