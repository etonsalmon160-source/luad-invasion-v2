#!/usr/bin/env python3
# 44_sc_ecm_checks.py —— ECM-high 成纤维的期别跳变：**先做反证检查**，再谈签名
#
# 现象：成纤维内部 ECM-high 比例 Normal 6.0 / AAH 9.7 / AIS 6.7 / MIA 4.0 / IAC 43.9 %
# 三个必须排除的替代解释（全项目栽过的那三种）：
#   A 深度——IAC 的成纤维 UMI 更多 ⇒ ECM 高表达基因检出更多
#   B 环境 RNA——IAC 组织 ECM 转录本弥散更多，被成纤维"吸收"
#   C 患者/样本——少数患者驱动
# 外加 D 双体（成纤维-上皮）——ECM 高 + 上皮基因同时高则可疑
import h5py, numpy as np, pandas as pd
from scipy import sparse
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
ORDER = ["Normal", "AAH", "AIS", "MIA", "IAC"]

D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as f:
    genes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["genes"][:]])
    M = sparse.csr_matrix((f["counts/data"][:], f["counts/indices"][:], f["counts/indptr"][:]),
                          shape=(len(D), len(genes)))
gi = {g: i for i, g in enumerate(genes)}
nUMI = np.asarray(M.sum(1)).ravel(); nGene = np.asarray((M > 0).sum(1)).ravel()
D["nUMI"] = nUMI; D["nGene"] = nGene

ECM = [g for g in ["COL1A1","COL1A2","COL3A1","COL6A3","FN1","CTHRC1","POSTN","SPARC","THBS2",
                   "BGN","COL5A1","COL5A2","LUM","VCAN","MMP2","TIMP1","ASPH","COL8A1"] if g in gi]
EPI = [g for g in ["EPCAM","KRT8","KRT18","KRT19","SFTPC","NAPSA"] if g in gi]
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = sparse.diags(1e4 / tot) @ M; Mn = Mn.tocsr(); Mn.data = np.log1p(Mn.data)
def score(gl):
    X = np.asarray(Mn[:, [gi[g] for g in gl]].todense())
    return ((X - X.mean(0)) / (X.std(0) + 1e-9)).mean(1)
D["ECM"] = score(ECM); D["EPI"] = score(EPI)

fib = D[D.L2 == "Fibroblast"].copy().reset_index(drop=True)
print(f"L2=Fibroblast：{len(fib):,} 细胞 / {fib.patient_id.nunique()} 患者\n")

print('=== 检查 A：IAC 的成纤维是不是更深（高表达 ECM 基因检出更多）===')
g = fib.groupby("stage")[["nUMI", "nGene"]].median().reindex(ORDER)
g["ECM_median"] = fib.groupby("stage").ECM.median().reindex(ORDER)
print(g.round(2).to_string())
print(f"  ⇒ nUMI 的 IAC/前驱 中位比 = {g.loc['IAC','nUMI']/g.loc[['AAH','AIS','MIA'],'nUMI'].median():.2f}")

print("\n=== 检查 A-2：**深度配平后** ECM-high 比例还在不在 ===")
fib["ng_bin"] = pd.qcut(fib.nGene, 10, labels=False, duplicates="drop")
r = []
for b, gb in fib.groupby("ng_bin"):
    n = gb.groupby("stage").size()
    h = gb[gb.ECM > gb.ECM.quantile(0.75)].groupby("stage").size()
    r.append(pd.DataFrame({"bin": b, "stage": n.index, "n": n.values,
                           "hi": h.reindex(n.index).fillna(0).values,
                           "hi_frac": (h.reindex(n.index).fillna(0) / n).values}))
R = pd.concat(r).pivot(index="stage", columns="bin", values="hi_frac").reindex(ORDER)
print(R.round(3).to_string())
print("  （每一列 = 一个 nGene 十分位；列内比较才是深度配平的）")
print("  ⇒ 每档 IAC 是否都最高：", bool((R.loc["IAC"] == R.max()).all()))
R.to_csv(f"{OUT}/scFib_ECMhi_by_depth_bin.tsv", sep="\t")

print("\n=== 检查 B：环境 RNA——上皮基因得分在成纤维里是否也随期别涨 ===")
print(fib.groupby("stage")[["ECM", "EPI"]].mean().reindex(ORDER).round(3).to_string())

print("\n=== 检查 D：ECM-high 的成纤维是不是上皮双体 ===")
fib["hi"] = fib.ECM > fib.ECM.quantile(0.75)
print(fib.groupby(["stage", "hi"])["EPI"].mean().unstack().reindex(ORDER).round(3).to_string())

print("\n=== 检查 C：逐患者 IAC 的 ECM-high 比例（是不是靠少数患者）===")
p = fib[fib.stage == "IAC"].groupby("patient_id").agg(n=("hi", "size"), hi_frac=("hi", "mean"))
p["n_all"] = fib.groupby("patient_id").size()
print(p.round(3).to_string())

print("\n=== 逐患者：前驱(AAH/AIS/MIA) 与 IAC 各自的比例 ===")
q = fib[fib.stage.isin(["AAH","AIS","MIA"])].groupby("patient_id").hi.mean().rename("precursor")
w = fib[fib.stage == "IAC"].groupby("patient_id").hi.mean().rename("IAC")
print(pd.concat([q, w], axis=1).round(3).dropna().to_string())
