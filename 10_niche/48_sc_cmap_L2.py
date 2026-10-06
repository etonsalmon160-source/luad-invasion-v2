#!/usr/bin/env python3
# 48_sc_cmap_L2.py —— 纯净的 SCMG(L2) 逆转测试：与 L1 同款管道，逐位可并排
#
# 与 L1(46_sc_cmap.py) 的唯一差别是库；口径逐条：
#   · X = 位移谱（正=扰动后升高）⇒ build_rank_tables 按 -X 降序，与 LINCS 的"表达降序"等价
#   · 分组键 = obs/dataset（预注册 §5.2 已登记"我自定适配"）
#   · 🔴 Q_ref：预注册**未给 L2 的 Q_ref**（§5.14 的 exemplar 是 L1 专有）
#     ⇒ 本项目取 **同 dataset 内全体的 |NCS| 百分位**（自校准；登记为 S 档）
#   · 零模型：同长度随机基因集，同管道 ⇒ 经验 p
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import h5py, numpy as np, pandas as pd
from importlib.machinery import SourceFileLoader
HERE = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963/10_niche"
TR = SourceFileLoader("tr", f"{HERE}/24_target_reversal.py").load_module()
ROOT = "/home/eto/luad_v2"; OUT = f"{ROOT}/results/10_niche/tr_singlecell"
TRD = f"{ROOT}/results/10_niche/target_reversal"
TOP = int(os.environ.get("TOP", 150)); N_NULL = int(os.environ.get("N_NULL", 200))
SEED = 20261005
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

sig_f = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else "sc"
S = pd.read_csv(sig_f, sep="\t").dropna()
UP = S.sort_values("lfc", ascending=False).gene.head(TOP).tolist()
DN = S.sort_values("lfc").gene.head(TOP).tolist()

H5 = "/home/eto/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad"
f = h5py.File(H5, "r")
def cat(path):                       # h5ad 的 categorical = categories + codes 两组
    g = f[path]
    cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
    return cats[np.asarray(g["codes"][:])]
genes = cat("var/gene_name")
gpos = {g: i for i, g in enumerate(genes)}
X = np.nan_to_num(np.asarray(f["X"][:], dtype=np.float32))
lab = pd.DataFrame({"dataset": cat("obs/dataset"),
                    "gene": cat("obs/perturbed_gene_name"),
                    "sign": np.asarray(f["obs/perturbation_sign"][:])})
f.close()
N = X.shape[1]
log(f"L2 库 {X.shape[0]:,} 扰动 × {N:,} 基因；dataset {lab.dataset.nunique()} 类")
rank, sv, _ = TR.build_rank_tables(X); del X
GRP = lab.dataset.values
GRPM = [GRP == g for g in np.unique(GRP)]

def tau_of(up_ix, dn_ix, need_tau=True):
    Ru, Cu = TR._rc(rank, sv, up_ix); e_u = TR.es_from_ranks(Ru, Cu, N, len(up_ix))
    Rd, Cd = TR._rc(rank, sv, dn_ix); e_d = TR.es_from_ranks(Rd, Cd, N, len(dn_ix))
    w = np.where(np.sign(e_u) != np.sign(e_d), (e_u - e_d) / 2.0, 0.0)
    ncs = np.empty_like(w)
    for m in GRPM:
        ww = w[m]
        mp = ww[ww > 0].mean() if (ww > 0).any() else np.nan
        mn = -ww[ww < 0].mean() if (ww < 0).any() else np.nan
        ncs[m] = np.where(ww > 0, ww / mp, np.where(ww < 0, ww / mn, 0.0))
    if not need_tau: return ncs
    # 🔴 L2 的 Q_ref = 同 dataset 内全体（自校准）
    t = np.empty_like(ncs)
    for m in GRPM:
        ref = np.sort(np.abs(ncs[m]))
        if ref.size == 0: continue
        t[m] = np.sign(ncs[m]) * 100.0 / ref.size * np.searchsorted(ref, np.abs(ncs[m]), side="left")
    return t

up_ix = np.array([gpos[g] for g in UP if g in gpos])
dn_ix = np.array([gpos[g] for g in DN if g in gpos])
log(f"{label}：映射进 L2 up {len(up_ix)} / down {len(dn_ix)}")
t = tau_of(up_ix, dn_ix)
neg = float((t <= -90).mean()); pos = float((t >= 90).mean())
log(f"结果：τ≤−90 = {neg:.4f}（{int((t<=-90).sum()):,}/{len(t):,}）| τ≥+90 = {pos:.4f} | mean τ = {t.mean():+.2f}")

# 零模型：同长度随机基因集
rng = np.random.default_rng(SEED); rates = []
for k in range(N_NULL):
    g = rng.choice(N, size=300, replace=False)
    tt = tau_of(g[:150], g[150:])
    rates.append(float((tt <= -90).mean()))
    if (k + 1) % 50 == 0: log(f"  null {k+1}/{N_NULL}  均值 {np.mean(rates):.4f}")
rates = np.array(rates)
z = (neg - rates.mean()) / rates.std(ddof=1); p = (1 + (rates >= neg).sum()) / (1 + len(rates))
log(f"零模型（n={len(rates)}）：均值 {rates.mean():.4f} ± {rates.std(ddof=1):.4f}，q95 {np.quantile(rates,.95):.4f}")
log(f"⇒ z = {z:+.2f}，经验 p = {p:.4f}  {'✅ 出界' if p < 0.05 else '❌ 落随机内'}")
pd.DataFrame({"rate_le_m90": rates}).to_csv(f"{TRD}/null_L2_{label}.tsv", sep="\t", index=False)

d = pd.DataFrame({"pert_gene": lab.gene.values, "sign": lab.sign.values, "dataset": lab.dataset.values, "tau": t})
d = d.sort_values("tau")
d.to_csv(f"{TRD}/scL2_{label}_ranked.tsv", sep="\t", index=False)
log("最负的 20 个基因扰动（候选基因靶点；sign=-1 表示敲低）：")
print(d.head(20)[["pert_gene", "sign", "dataset", "tau"]].to_string(index=False))
