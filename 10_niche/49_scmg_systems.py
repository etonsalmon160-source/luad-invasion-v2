#!/usr/bin/env python3
# 49_scmg_systems.py —— SCMG「系统生物学」口径：把疾病状态**压回正常**（不用 CMap 的 KS/τ）
#
# 为什么换口径：CMap 的 τ 是**秩统计量 + 单基因表查询**，要先把签名截成 top150，丢信息、还依赖 Q_ref。
# 系统生物学的看法是——把"疾病"当成**状态空间里的一个方向**：
#     d = 域状态 − 正常状态        （一个 18,108 维向量）
#     每个扰动 p 自带位移向量 X_p  （SCMG 直接给，全部基因、不截断）
#     打回正常的分数 = **−cos(X_p, d)**  （位移与"疾病轴"反向的程度）
# 两个轴：d_fib = 成纤维的疾病轴；d_comp = 39 个细胞型各自疾病轴的合体
# 零模型：把 d 的基因标签打乱（保持 |d| 分布 ⇒ 只破坏"方向与基因的对应"）
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import h5py, numpy as np, pandas as pd
from scipy import sparse
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
SEED = 20261005; N_NULL = int(os.environ.get("N_NULL", 200))
ORDER = ["Normal", "AAH", "AIS", "MIA", "IAC"]
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ── ① snRNA：逐患者基因均值（成纤维 + 各细胞型）──
D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as f:
    sgenes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["genes"][:]])
    M = sparse.csr_matrix((f["counts/data"][:], f["counts/indices"][:], f["counts/indptr"][:]),
                          shape=(len(D), len(sgenes)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
pid = D.patient_id.values; stg = D.stage.values; l2 = D.L2.values
PRE = np.isin(stg, ["AAH", "AIS", "MIA"]); IAC = stg == "IAC"
log(f"snRNA {Mn.shape[0]:,} × {Mn.shape[1]:,}；患者 {len(set(pid))}")

def axis_for(mask_cells):
    """患者内配对：median_p(mean(IAC∩mask) − mean(pre∩mask))"""
    parts = []
    for p in np.unique(pid):
        ip = (pid == p) & mask_cells
        a = ip & IAC; b = ip & PRE
        if a.sum() >= 20 and b.sum() >= 20:
            parts.append(np.asarray(Mn[a].mean(0)).ravel() - np.asarray(Mn[b].mean(0)).ravel())
    return (np.median(np.vstack(parts), axis=0), len(parts)) if parts else (None, 0)

d_fib, n1 = axis_for(l2 == "Fibroblast")
log(f"A 轴 d_fib（成纤维疾病轴）：{n1} 例患者")
axes = []
for t in sorted(set(l2)):
    a, n = axis_for(l2 == t)
    if a is not None and n >= 8:
        axes.append(a / (np.linalg.norm(a) + 1e-12))
log(f"B 轴 d_comp：合了 {len(axes)} 个细胞型的疾病轴（各归一后平均）")
d_comp = np.mean(np.vstack(axes), axis=0)

# ── ② SCMG 库 ──
f = h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r")
def cat(p):
    g = f[p]; c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
    return c[np.asarray(g["codes"][:])]
G = cat("var/gene_name")
X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
lab = pd.DataFrame({"gene": cat("obs/perturbed_gene_name"), "sign": np.asarray(f["obs/perturbation_sign"][:]),
                    "dataset": cat("obs/dataset")})
f.close()
log(f"SCMG X {X.shape[0]:,} × {X.shape[1]:,}")

# ── ③ 对齐基因空间 ──
pos = {g: i for i, g in enumerate(G)}
keep = np.array([i for i, g in enumerate(sgenes) if g in pos])
tgt = np.array([pos[sgenes[i]] for i in keep])
log(f"两库共有基因 {len(keep):,} / {len(sgenes):,}；参与投影")

Xn = X[:, tgt].astype(np.float64)
Xnorm = np.linalg.norm(Xn, axis=1) + 1e-12

def scores(dvec):
    d = dvec[keep]; d = d - d.mean()
    dn = np.linalg.norm(d) + 1e-12
    cos = (Xn @ d) / (Xnorm * dn)          # 余弦
    proj = (Xn @ d) / dn                   # 有符号投影长度
    return -cos, -proj                      # 取负：越大 = 越把状态推回正常

cos_f, proj_f = scores(d_fib)
cos_c, proj_c = scores(d_comp)
log(f"d_fib: cos 范围 [{cos_f.min():.3f},{cos_f.max():.3f}] 均值 {cos_f.mean():+.4f}")
log(f"d_comp: cos 范围 [{cos_c.min():.3f},{cos_c.max():.3f}] 均值 {cos_c.mean():+.4f}")

# ── ④ 零模型：打乱 d 的基因标签 ──
rng = np.random.default_rng(SEED)
null_hi = []                            # 随机轴下"cos ≥ 观测 top 阈值"的比例
for k in range(N_NULL):
    dd = d_fib[keep].copy(); rng.shuffle(dd)
    ddn = np.linalg.norm(dd) + 1e-12
    cc = -((Xn @ dd) / (Xnorm * ddn))
    null_hi.append(np.quantile(cc, 0.99))
null_hi = np.array(null_hi)
log(f"零模型（打乱 d，n={N_NULL}）：随机轴下 top-1% cos 的中位 {np.median(null_hi):.3f}")

# ── ⑤ 汇总：按 perturbed_gene 聚合（跨 dataset 取中位 + 一致性）──
res = lab.copy(); res["cos_fib"] = cos_f; res["proj_fib"] = proj_f; res["cos_comp"] = cos_c
agg = res.groupby("gene").agg(n=("cos_fib", "size"), cos_fib=("cos_fib", "median"),
                              cos_comp=("cos_comp", "median"),
                              frac_pos=("cos_fib", lambda v: (v > 0).mean()))
agg = agg[agg.n >= 3].sort_values("cos_fib", ascending=False)
agg.to_csv(f"{OUT}/scmg_systems_ranked.tsv", sep="\t")
res.to_csv(f"{OUT}/scmg_systems_all.tsv", sep="\t", index=False)
print("\n=== 把成纤维疾病轴推回正常最强的 25 个基因扰动（cos 越大越反向）===")
print(agg.head(25).round(4).to_string())
print('\n=== 尾部（最同向 = 最像加重疾病）10 个 ===')
print(agg.tail(10).round(4).to_string())
print(f"\n观测 top-1% 的 cos 分位：{np.quantile(cos_f,0.99):.3f}（零模型中位 {np.median(null_hi):.3f}）")
