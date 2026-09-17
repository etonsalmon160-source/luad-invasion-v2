#!/usr/bin/env python3
"""GP6 · 双标准交叉注释（Harmony 主口径 r*=0.6）

标准 A：法则2 规范 marker（`marker_panel.py`，出处为一次文献）→ sc.tl.score_genes 模块分 → 簇均值 argmax
标准 B：CellTypist `Human_Lung_Atlas`（Sikkema 2023 HLCA）逐细胞判定
        ＋ sc.tl.rank_genes_groups('wilcoxon') top-50 × 图谱基因集 Jaccard

一致性：逐细胞 Cohen's κ ＋ 逐簇一致率。**全部分歧写入 CSV，绝不自动裁决。**
恶性身份**不**由本脚本判定（R2：须 CNV 证真）。

结构（2026-09-17 重排）：主结果（κ / 分歧清单）**先落盘**，再跑慢的 B2 wilcoxon 腿。
原因：`rank_genes_groups('wilcoxon')` 在 41 万核 × 45 簇上单线程、单簇稠密双 argsort，
实测 >24 min 未完成（run1 留档 `logs/GP6_dual_annotation.run1_killed_B2slow.log`）。
B2 是**附属腿**（只补 Jaccard 列），不得阻塞主结果。

输出（`results/05_annotation/`）：
  gp6_scores.csv.gz                 逐细胞 6 谱系模块分（标准 A 原始读数）
  gp6_cell_labels.csv.gz            逐细胞 双标准标签（供下游 GP8a 取上皮子集）
  annotation_disagreement.csv       全部分歧（用户复核用）
  gp6_cluster_labels.csv            逐簇 标准A/B1/B2 判定 + 一致率 + 驱动基因
  gp6_metrics.json                  κ / 一致率 / 各 seed 敏感性
  gp6_manifest.json                 产物哈希与溯源
"""
import os, sys, json, hashlib, time, platform

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
from sklearn.metrics import cohen_kappa_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import marker_panel as MP
from celltypist.models import Model
import celltypist

ROOT = "/home/eto/luad_v2"
H5 = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
CLU = f"{ROOT}/results/04_integration/seurat_trad/full/clusters.csv.gz"
OUT = f"{ROOT}/results/05_annotation"
os.makedirs(OUT, exist_ok=True)

RES = os.environ.get("GP6_RES", "0.6")
SEEDS = [0, 1, 2, 3, 4]
FROZEN_SEED = 0          # GP5 只给跨种子均值，未指定代表种子 ⇒ 本脚本显式指定并报全种子敏感性
SKIP_B2 = os.environ.get("GP6_SKIP_B2", "0") == "1"
CTRL_SIZE = 50
SCORE_SEED = 0
JACCARD_TOPN = 50
ATLAS_MARKER_TOPN = 20
DRIVER_TOPN = 15

T0 = time.time()
LOG = []


def log(msg):
    s = f"[{time.time()-T0:8.1f}s] {msg}"
    print(s, flush=True)
    LOG.append(s)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


# ================================================================ 1 载入 + 归一化
log(f"载入 {H5}")
adata = ad.read_h5ad(H5)
log(f"  shape={adata.shape}  nnz={adata.X.nnz:,}")
assert adata.shape == (413697, 18069), adata.shape
assert all(isinstance(g, str) and g for g in adata.var_names[:50]), "var_names 非字符串"
assert not any(str(g).startswith("b'") for g in adata.var_names), "var_names 疑似 bytes 泄漏"

log("归一化 normalize_total(1e4) + log1p")
assert adata.X.min() >= 0 and np.allclose(adata.X.data, np.round(adata.X.data))
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
log("  完成")

# ================================================================ 2 分群（r × 5 种子）
log(f"载入 r={RES} 分群")
clu = pd.read_csv(CLU)
assert len(clu) == adata.n_obs
assert (clu["cell_barcode"].to_numpy() == adata.obs_names.to_numpy()).all(), "分群表与 h5ad 细胞顺序不一致"
seed_cols = {s: f"harmony_res{RES}_seed{s}" for s in SEEDS}
for s, c in seed_cols.items():
    assert c in clu.columns, c
    adata.obs[f"clu_s{s}"] = pd.Categorical(clu[c].astype(str))
    log(f"  seed{s}: {clu[c].nunique()} 簇")

KEY = f"clu_s{FROZEN_SEED}"
codes = adata.obs[KEY].cat.codes.to_numpy()
nclu = codes.max() + 1
clu_ids = [str(x) for x in adata.obs[KEY].cat.categories]

# ================================================================ 3 标准 A · 模块分
used, missing = MP.resolve_missing(adata.var_names)
log(f"标准 A 面板：可用 {sum(len(v) for v in used.values())} / 缺失 {missing}")
score_cols = []
for lin in MP.LINEAGES:
    gl = used[lin]
    assert len(gl) >= 5, f"{lin} 可用 marker 过少（{len(gl)}）"
    name = f"scoreA_{lin}"
    sc.tl.score_genes(adata, gl, ctrl_size=CTRL_SIZE, random_state=SCORE_SEED, score_name=name)
    score_cols.append(name)
    log(f"  {lin:6s} n_genes={len(gl):2d} 均值={adata.obs[name].mean():+.4f}")

S = adata.obs[score_cols].to_numpy(dtype=np.float64)
lin_arr = np.array(MP.LINEAGES)


def cluster_argmax(seed):
    """每个簇的谱系 argmax（标准 A 的簇级判定），返回 (top1, top2, margin)。"""
    c = adata.obs[f"clu_s{seed}"].cat.codes.to_numpy()
    n = c.max() + 1
    m = np.zeros((n, len(lin_arr)))
    for j in range(n):
        mm = c == j
        if mm.any():
            m[j] = S[mm].mean(axis=0)
    o = np.argsort(-m, axis=1)
    return c, n, m, o[:, 0], (o[:, 1] if m.shape[1] > 1 else o[:, 0])


a_labels = {}
for s in SEEDS:
    c, n, m, t1, _ = cluster_argmax(s)
    a_labels[s] = lin_arr[t1][c]
    adata.obs[f"A_seed{s}"] = pd.Categorical(a_labels[s])

cF, nF, mF, t1F, t2F = cluster_argmax(FROZEN_SEED)
adata.obs["A_top1"] = pd.Categorical(lin_arr[t1F][cF])
adata.obs["A_top2_lin"] = pd.Categorical(lin_arr[t2F][cF])
adata.obs["A_margin"] = (mF[np.arange(nF), t1F] - mF[np.arange(nF), t2F])[cF]
assert codes.tolist() == cF.tolist()

# ================================================================ 4 标准 B · CellTypist
mpath = os.path.expanduser(MP.CELLTYPIST_MODEL["path"])
ct_sha = sha256(mpath)
MP.CELLTYPIST_MODEL["sha256"] = ct_sha
log(f"标准 B：CellTypist 模型 sha256={ct_sha[:16]}…")
model = Model.load(mpath)
n_feat = len(model.features)
n_have = len([f for f in model.features if f in set(adata.var_names)])
log(f"  模型特征 {n_feat}，本项目含 {n_have}（{n_have/n_feat*100:.1f}%）")

log("  CellTypist.annotate（含 majority_voting，over_clustering=seed0 簇）")
res_ct = celltypist.annotate(adata, model=model, majority_voting=True,
                             over_clustering=adata.obs[KEY].astype(str).to_numpy())
pred = res_ct.predicted_labels
ct_cell = pred["predicted_labels"].astype(str).to_numpy()
ct_clu_mv = pred["majority_voting"].astype(str).to_numpy()
conf = res_ct.probability_matrix.max(axis=1).to_numpy()

adata.obs["B_celltypist_type"] = pd.Categorical(ct_cell)
adata.obs["B_celltypist_mv"] = pd.Categorical(ct_clu_mv)
adata.obs["B_confidence"] = conf

MAPFULL = MP.lineage_map_full()
AMBIG = set(MP.CELLTYPIST_AMBIGUOUS)
b_lin = np.array([MAPFULL.get(t, "未归属") for t in ct_cell])
b_is_ambig = np.array([t in AMBIG for t in ct_cell])
adata.obs["B_lineage"] = pd.Categorical(b_lin)
adata.obs["B_is_ambiguous"] = b_is_ambig
log("  标准 B 谱系分布: " + ", ".join(f"{k}={v}" for k, v in
      pd.Series(b_lin).value_counts().items()))

# ================================================================ 5 一致性统计
log("一致性统计")
kappa = {}
keep_amb = ~b_is_ambig
for s in SEEDS:
    a = a_labels[s]
    kappa[f"seed{s}"] = dict(
        kappa_pos_cell=float(cohen_kappa_score(a, b_lin)),
        kappa_excl_ambig=float(cohen_kappa_score(a[keep_amb], b_lin[keep_amb])),
        agreement_rate_pos_cell=float((a == b_lin).mean()),
        agreement_rate_excl_ambig=float((a[keep_amb] == b_lin[keep_amb]).mean()),
    )
    log(f"  seed{s}: κ={kappa[f'seed{s}']['kappa_pos_cell']:.4f}  "
        f"一致率={kappa[f'seed{s}']['agreement_rate_pos_cell']:.4f}  "
        f"(剔判断项 κ={kappa[f'seed{s}']['kappa_excl_ambig']:.4f})")

aF = a_labels[FROZEN_SEED]
agree = aF == b_lin
kF = kappa[f"seed{FROZEN_SEED}"]["kappa_pos_cell"]
verdict = "PASS" if kF >= 0.80 else ("FLAG_人工复核" if kF >= 0.60 else "STOP")

# ---- 廉价驱动基因代理（簇均值 − 全局均值），B2 未完成前先给用户一个可读的驱动列 ----
log("计算簇富集驱动基因（mean-diff 代理）")
Gind = csr_matrix((np.ones(adata.n_obs), (np.arange(adata.n_obs), codes)),
                  shape=(adata.n_obs, nclu))
csum = np.asarray((Gind.T @ adata.X).todense())
cmean = csum / np.bincount(codes, minlength=nclu)[:, None]
gmean = np.asarray(adata.X.mean(axis=0)).ravel()
proxy_top = {}
for j in range(nclu):
    d = cmean[j] - gmean
    o = np.argsort(-d)[:DRIVER_TOPN]
    proxy_top[str(j)] = [adata.var_names[i] for i in o if d[i] > 0]
log("  完成")

# ================================================================ 6 ★ 主结果先落盘
log("★ 写出主结果（B2 wilcoxon 尚未运行）")
pd.DataFrame(index=adata.obs_names).assign(
    cell_barcode=adata.obs_names.to_numpy(),
    sample_id=adata.obs["sample_id"], patient_id=adata.obs["patient_id"],
    stage=adata.obs["stage"],
    **{f"A_seed{s}": a_labels[s] for s in SEEDS},
    A_frozen=aF, B_celltypist_type=ct_cell, B_lineage=b_lin,
    B_confidence=conf, agree=agree,
).to_csv(f"{OUT}/gp6_cell_labels.csv.gz", index=False, compression="gzip")

pd.DataFrame(S, columns=MP.LINEAGES).assign(
    cell_barcode=adata.obs_names.to_numpy(), cluster=codes
).to_csv(f"{OUT}/gp6_scores.csv.gz", index=False, compression="gzip")

idx = np.where(~agree)[0]
dis_df = pd.DataFrame({
    "cell_barcode": adata.obs_names.to_numpy()[idx],
    "sample_id": adata.obs["sample_id"].to_numpy()[idx],
    "patient_id": adata.obs["patient_id"].to_numpy()[idx],
    "stage": adata.obs["stage"].to_numpy()[idx],
    "cluster": codes[idx],
    "A_lineage(标准A)": aF[idx],
    "A_top2_lineage": adata.obs["A_top2_lin"].to_numpy()[idx],
    "A_margin": adata.obs["A_margin"].to_numpy()[idx],
    "B_celltypist_type(标准B)": ct_cell[idx],
    "B_lineage_mapped": b_lin[idx],
    "B_confidence": conf[idx],
    "B_is_ambiguous_type": b_is_ambig[idx],
    "cluster_driver_genes_meandiff": [";".join(proxy_top[str(codes[i])]) for i in idx],
})
dis_df.to_csv(f"{OUT}/annotation_disagreement.csv", index=False)
log(f"  分歧细胞 {len(dis_df):,} / {adata.n_obs:,}（{len(dis_df)/adata.n_obs*100:.2f}%）")

pair = dis_df.groupby(["A_lineage(标准A)", "B_lineage_mapped"]).size().sort_values(ascending=False)
log("  分歧对 top10:\n" + pair.head(10).to_string())

partial = dict(
    resolution=float(RES), frozen_seed=FROZEN_SEED, n_seeds=len(SEEDS),
    n_cells=int(adata.n_obs), n_clusters=int(nclu),
    kappa_by_seed=kappa, kappa_frozen=float(kF), verdict=verdict,
    n_disagreement=int(len(dis_df)), disagreement_rate=float(len(dis_df) / adata.n_obs),
    panel_missing=missing, panel_used={k: len(v) for k, v in used.items()},
    celltypist=dict(model_path=mpath, sha256=ct_sha, n_celltypes=61,
                    n_features=n_feat, n_features_present=n_have,
                    source=MP.CELLTYPIST_MODEL["source"]),
    ambiguous_celltypist_types=MP.CELLTYPIST_AMBIGUOUS,
    disagreement_pairs_top={f"{a}|{b}": int(v) for (a, b), v in pair.head(15).items()},
    b2_status="pending",
)
with open(f"{OUT}/gp6_metrics.json", "w", encoding="utf-8") as fh:
    json.dump(partial, fh, ensure_ascii=False, indent=2)
log(f"★ 主结果已落盘。冻结 seed{FROZEN_SEED} κ={kF:.4f} → 判定 {verdict}")
log(f"  分歧清单：{OUT}/annotation_disagreement.csv")

# ================================================================ 7 标准 B2（慢）
jac, clu_top = {}, {}
if SKIP_B2:
    log("GP6_SKIP_B2=1 ⇒ 跳过 B2 wilcoxon 腿")
else:
    log(f"标准 B2：rank_genes_groups(wilcoxon) top-{JACCARD_TOPN} × 图谱集 Jaccard（慢，单线程）")
    sc.tl.rank_genes_groups(adata, groupby=KEY, method="wilcoxon", n_genes=JACCARD_TOPN)
    names = adata.uns["rank_genes_groups"]["names"]
    clu_top = {str(c): [str(g) for g in names[c]] for c in names.dtype.names}
    log(f"  wilcoxon 完成，簇数 {len(clu_top)}")

    atlas_sets = {}
    for lin in MP.LINEAGES:
        types = [ct for ct, L in MAPFULL.items() if L == lin]
        genes = set()
        for ct in types:
            genes.update(str(g) for g in model.extract_top_markers(ct, top_n=ATLAS_MARKER_TOPN))
        atlas_sets[lin] = genes
        log(f"  图谱集 {lin:6s} {len(types):2d} 类 / {len(genes)} 基因")

    for j in range(nclu):
        top = set(clu_top[str(j)])
        jac[str(j)] = {lin: (len(top & atlas_sets[lin]) / len(top | atlas_sets[lin]))
                       if (top | atlas_sets[lin]) else 0.0 for lin in MP.LINEAGES}

# ================================================================ 8 逐簇表 + 最终指标
clu_tab = []
for j in range(nclu):
    m = codes == j
    jrow = jac.get(str(j))
    clu_tab.append(dict(
        cluster=j, n_cells=int(m.sum()),
        A_lineage=lin_arr[t1F[j]], A_top2=lin_arr[t2F[j]],
        A_margin=float(mF[j, t1F[j]] - mF[j, t2F[j]]),
        A_scores=";".join(f"{l}={mF[j,k]:+.3f}" for k, l in enumerate(MP.LINEAGES)),
        B_lineage_mode=pd.Series(b_lin[m]).mode().iat[0],
        B_majority_voting=ct_clu_mv[m][0],
        B_agreement_rate=float(agree[m].mean()),
        B_mean_confidence=float(conf[m].mean()),
        B_ambiguous_frac=float(b_is_ambig[m].mean()),
        B2_jaccard_argmax=(max(jrow, key=jrow.get) if jrow else None),
        B2_jaccard_value=(float(jrow[max(jrow, key=jrow.get)]) if jrow else None),
        B2_jaccard_all=(";".join(f"{l}={jrow[l]:.3f}" for l in MP.LINEAGES) if jrow else None),
        driver_genes_meandiff=";".join(proxy_top[str(j)]),
        driver_genes_wilcoxon=(";".join(clu_top[str(j)][:DRIVER_TOPN]) if clu_top else None),
    ))
clu_df = pd.DataFrame(clu_tab)
clu_df["tri_agree"] = (clu_df["A_lineage"] == clu_df["B_lineage_mode"]) & \
                      ((clu_df["A_lineage"] == clu_df["B2_jaccard_argmax"]) if jac else True)
clu_df.to_csv(f"{OUT}/gp6_cluster_labels.csv", index=False)
log(f"逐簇 A/B1/B2 三方一致: {int(clu_df['tri_agree'].sum())}/{len(clu_df)}")

if jac:
    dis_df["B2_jaccard_argmax"] = [max(jac[str(codes[i])], key=jac[str(codes[i])].get) for i in idx]
    dis_df["cluster_driver_genes_wilcoxon"] = [";".join(clu_top[str(codes[i])][:DRIVER_TOPN]) for i in idx]
    dis_df.to_csv(f"{OUT}/annotation_disagreement.csv", index=False)
    log("  分歧清单已补 B2 列并重写")

partial.update(
    b2_status=("skipped" if SKIP_B2 else "done"),
    tri_agree_clusters=int(clu_df["tri_agree"].sum()),
    cluster_tri_disagree=clu_df.loc[~clu_df["tri_agree"], "cluster"].tolist(),
    jaccard_params=dict(topn=JACCARD_TOPN, atlas_marker_topn=ATLAS_MARKER_TOPN),
)
with open(f"{OUT}/gp6_metrics.json", "w", encoding="utf-8") as fh:
    json.dump(partial, fh, ensure_ascii=False, indent=2)

man = dict(
    script=os.path.abspath(__file__), script_sha256=sha256(os.path.abspath(__file__)),
    panel_sha256=sha256(f"{ROOT}/05_annotation/marker_panel.py"),
    inputs=dict(h5ad=H5, h5ad_sha256=sha256(H5), clusters=CLU,
                clusters_sha256=sha256(CLU), celltypist_model_sha256=ct_sha),
    params=dict(res=RES, ctrl_size=CTRL_SIZE, score_seed=SCORE_SEED,
                jaccard_topn=JACCARD_TOPN, atlas_marker_topn=ATLAS_MARKER_TOPN,
                frozen_seed=FROZEN_SEED, skip_b2=SKIP_B2),
    outputs={f: sha256(f"{OUT}/{f}") for f in
             ["gp6_cell_labels.csv.gz", "gp6_scores.csv.gz", "gp6_cluster_labels.csv",
              "annotation_disagreement.csv", "gp6_metrics.json"]},
    env=dict(python=platform.python_version(), scanpy=sc.__version__,
             anndata=ad.__version__, celltypist=celltypist.__version__,
             numpy=np.__version__, pandas=pd.__version__),
    runtime_sec=round(time.time() - T0, 1),
)
with open(f"{OUT}/gp6_manifest.json", "w", encoding="utf-8") as fh:
    json.dump(man, fh, ensure_ascii=False, indent=2)

log(f"判定：{verdict}（冻结 seed{FROZEN_SEED} κ={kF:.4f}）；B2={partial['b2_status']}")
print(f"\nGP6 DONE -> {OUT}")
