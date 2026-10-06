#!/usr/bin/env python3
# 51_domain_composite.py —— 路 1：**复合域签名**（把整个域压回正常）
#
# 思路：域不是一个细胞，是**一群细胞按组成混在一起**。
#       ⇒ 域的程序 = Σ_c w_c · axis_c
#          w_c    = 细胞型 c 在 D3 域里的占比（RCTD 39 亚型组成，56 张切片均值）
#          axis_c = 细胞型 c 自己的疾病轴（患者内配对 IAC − 前驱）
# 两个权重版本：① 域内丰度 w(D3)；② 组成变化 w(D3) − w(全部域)  ← 后者才是"域相对于常态多了什么"
# 打分：系统生物学口径 −cos(X_p, d)（与 49_ 同），先用 L2(SCMG) 跑（快、不占大内存）；
#       同时导出 top150 up/down 供 L1(CMap) 后续跑。
import os, time, glob, re
os.environ.setdefault("OMP_NUM_THREADS", "4")
import h5py, numpy as np, pandas as pd
from scipy import sparse
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
KD = f"{ROOT}/results/10_niche/kstar_diag"; TOP = 150; N_NULL = 300; SEED = 20261005
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
nz = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())

# ── ① 域组成（D3 = archetype 3）——来自 §17 的 d7_archetypes.tsv（7×39 组成表）──
A7 = pd.read_csv(f"{KD}/d7_archetypes.tsv", sep="\t")
META = {"archetype", "n_domain", "n_slide", "depth_med", "top5"}
CT = [c for c in A7.columns if c not in META]
assert len(CT) == 39, len(CT)
w3 = A7[A7.archetype == 3][CT].iloc[0].astype(float)
wall = A7[CT].astype(float).mean()
delta = (w3 - wall)
log(f"D3 域（archetype 3）；参与细胞型 {len(CT)}")
log(f"  D3 里 top5：{', '.join(f'{k}={w3[k]:.3f}' for k in w3.nlargest(5).index)}")
log(f"  组成变化 top5：{', '.join(f'{k}={delta[k]:+.3f}' for k in delta.nlargest(5).index)}")
D3MAP = {nz(k): k for k in CT}
w3v = w3.values.copy(); dlt = delta.values.copy()
names_ordered = list(CT)

# ── ② 各细胞型的疾病轴 ──
D = pd.read_parquet(f"{OUT}/sc_cellmap.parquet")
with h5py.File(f"{ROOT}/results/08_spatial_deconv/reference_d.h5", "r") as f:
    sg = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["genes"][:]])
    M = sparse.csr_matrix((f["counts/data"][:], f["counts/indices"][:], f["counts/indptr"][:]),
                          shape=(len(D), len(sg)))
tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
Mn = (sparse.diags(1e4 / tot) @ M).tocsr(); Mn.data = np.log1p(Mn.data)
pid = D.patient_id.values; l2 = D.L2.values
PRE = np.isin(D.stage.values, ["AAH", "AIS", "MIA"]); IAC = D.stage.values == "IAC"
AX = {}
for t in np.unique(l2):
    m = l2 == t; ps = []
    for p in np.unique(pid):
        ip = (pid == p) & m; a = ip & IAC; b = ip & PRE
        if a.sum() >= 20 and b.sum() >= 20:
            ps.append(np.asarray(Mn[a].mean(0)).ravel() - np.asarray(Mn[b].mean(0)).ravel())
    if len(ps) >= 8: AX[t] = np.median(np.vstack(ps), axis=0)
log(f"疾病轴：{len(AX)} 个细胞型")
del M, Mn
ref2ct = D3MAP
d_dom = np.zeros(len(sg)); d_dlc = np.zeros(len(sg)); used = 0
for i, nm in enumerate(names_ordered):
    key = next((k for k in AX if nz(k) == nz(nm)), None)
    if key is None: continue
    d_dom += w3v[i] * AX[key]; d_dlc += dlt[i] * AX[key]; used += 1
log(f"复合轴：用了 {used} 个细胞型；|d_dom|={np.linalg.norm(d_dom):.2f}  |d_dlc|={np.linalg.norm(d_dlc):.2f}")

# ── ③ L2 系统生物学投影 ──
f = h5py.File("/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad", "r")
def cat(p):
    g = f[p]; c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]); return c[np.asarray(g["codes"][:])]
G = cat("var/gene_name"); X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
lab = pd.DataFrame({"gene": cat("obs/perturbed_gene_name"), "sign": np.asarray(f["obs/perturbation_sign"][:])})
f.close()
pos = {g: i for i, g in enumerate(G)}
keep = np.array([i for i, g in enumerate(sg) if g in pos]); tgt = np.array([pos[sg[i]] for i in keep])
Xn = X[:, tgt].astype(np.float64); Xnrm = np.linalg.norm(Xn, axis=1) + 1e-12
log(f"共有基因 {len(keep):,}")

def run(axis, tag):
    d0 = axis[keep] - axis[keep].mean(); dn = np.linalg.norm(d0) + 1e-12
    obs = -((Xn @ d0) / (Xnrm * dn))
    rng = np.random.default_rng(SEED); null = np.empty((len(obs), N_NULL), dtype=np.float32)
    for k in range(N_NULL):
        dd = d0.copy(); rng.shuffle(dd)
        null[:, k] = -((Xn @ dd) / (Xnrm * (np.linalg.norm(dd) + 1e-12))).astype(np.float32)
    p = np.array([(1 + (null[i] >= obs[i]).sum()) / (1 + N_NULL) for i in range(len(obs))])
    A = lab.copy(); A["cos"] = obs; A["p"] = p
    g = A.groupby(["gene", "sign"]).agg(n=("cos", "size"), cos=("cos", "median"),
                                        p=("p", "median"),
                                        frac=("cos", lambda v: (v > 0).mean())).reset_index()
    g = g[g.n >= 3]; g["dir"] = np.where(g.sign < 0, "KD", "OE"); g = g.sort_values("cos", ascending=False)
    g.to_csv(f"{OUT}/composite_{tag}_ranked.tsv", sep="\t", index=False)
    log(f"[{tag}] top1% cos 观测 {np.quantile(obs,0.99):.4f} vs 零模型 {np.median(np.quantile(null,0.99,axis=0)):.4f}"
        f" | p<0.05 的(基因×方向) {int((g.p<0.05).sum())}/{len(g)}")
    print(g.head(15)[["gene", "dir", "n", "cos", "p", "frac"]].round(4).to_string(index=False))
    return obs

log("=== 轴 A：域内丰度加权 d_dom  ===");  o1 = run(d_dom, "ddom")
log("=== 轴 B：组成变化加权 d_dlc ===");  o2 = run(d_dlc, "ddlc")

# ── ④ 导出 top150 供 L1 CMap ──
for nm, ax in [("ddom", d_dom), ("ddlc", d_dlc)]:
    r = pd.DataFrame({"gene": sg, "lfc": ax}).replace([np.inf, -np.inf], np.nan).dropna()
    r.sort_values("lfc", ascending=False).to_csv(f"{OUT}/composite_{nm}_lfc.tsv", sep="\t", index=False)
log("已导出 composite_{ddom,ddlc}_lfc.tsv 供 L1")
