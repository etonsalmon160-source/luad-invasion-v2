#!/usr/bin/env python3
# SCMG 分支 · 第一腿 · 零样本投影到全局流形
#
# 常量全部**写死在 results/06_scmg/SCMG_BRANCH_PREREG.md §四**（计算前登记），本脚本只执行。
# 流程照官方 tutorial 逐行对齐：
#   SCMG/docs/source/tutorials/project_single_cell_states_onto_the_global_cell_state_manifold.ipynb
#
# 跑法：python3 06_scmg/01_project_to_global_manifold.py
#
# ⚠️ 用**主环境**（py3.8 / numpy 1.22.4 / scanpy 1.9.8 / torch 2.4.1），不新建环境、不装包。
#    实测 import 通过（预注册 §二）。

import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, '/home/eto/scmg_workspace/SCMG')
from scmg.model.contrastive_embedding import CellEmbedder, embed_adata
from scmg.model.cell_type_search import CellTypeSearcher

import torch

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/06_scmg'
os.makedirs(D, exist_ok=True)

# ——— 冻结常量（预注册 §二 / §四）———
SCMG_DIR = '/home/eto/scmg_workspace'
STD_CSV = f'{SCMG_DIR}/SCMG/scmg/data/standard_genes.csv'
MODEL_PT = f'{SCMG_DIR}/models/embedder/model.pt'
STATE_PTH = f'{SCMG_DIR}/models/embedder/best_state_dict.pth'
REF_H5AD = f'{SCMG_DIR}/hf_data/ref_global_cell_state_manifold.h5ad'
QUERY_H5AD = f'{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad'

P_DEVICE = 'cpu'
P_BATCH_SIZE = 8192
P_PRE_NORMALIZED = False
P_N_JOBS = 20
P_EMB_KEY = 'X_scmg'
N_CELLS_EXPECTED = 413697
N_GENES_EXPECTED = 18069

t0 = time.time()
def log(msg):
    print(f'[{time.time()-t0:7.1f}s] {msg}', flush=True)

# ——— 1. 基因映射（用 SCMG 自带 standard_genes.csv；重名取首次出现）———
log('读 SCMG standard_genes.csv ...')
std = pd.read_csv(STD_CSV)
name2id = {}
n_ambiguous_seen = 0
amb_hits = []
for nm, gid in zip(std.human_name.astype(str).str.upper(), std.human_id):
    if nm in name2id:
        n_ambiguous_seen += 1
    else:
        name2id[nm] = gid          # 冻结：首次出现优先
log(f'  标准基因 {len(std)} 行，唯一符号 {len(name2id)}，重名行 {n_ambiguous_seen}')

log('读我们的表达对象 ...')
A = ad.read_h5ad(QUERY_H5AD)
assert A.shape == (N_CELLS_EXPECTED, N_GENES_EXPECTED), f'维度不符 {A.shape}'
log(f'  {A.shape[0]} 核 × {A.shape[1]} 基因；X 类型 {type(A.X).__name__} {A.X.dtype}')

ours_upper = [str(g).upper() for g in A.var.index]
eids = [name2id.get(g) for g in ours_upper]

# 撞上重名的基因（预注册 §三 定死处置）
_name_counts = std.human_name.astype(str).str.upper().value_counts()
amb_symbols = set(_name_counts[_name_counts > 1].index)
amb_hits = sorted({g for g in ours_upper if g in amb_symbols})

unmapped = [A.var.index[i] for i, e in enumerate(eids) if e is None]
with open(f'{D}/scmg_unmapped_genes.txt', 'w') as fh:
    fh.write('\n'.join(unmapped) + '\n')
coverage = 1 - len(unmapped) / len(eids)
log(f'  映射上 {len(eids)-len(unmapped)}/{len(eids)} = {coverage*100:.2f}%；未映射 {len(unmapped)} 个 → scmg_unmapped_genes.txt')
log(f'  撞上 SCMG 重名符号的基因：{amb_hits if amb_hits else "无"}（处置：取首次出现，冻结）')

keep = np.array([i for i, e in enumerate(eids) if e is not None])
eid_keep = [eids[i] for i in keep]
log(f'  建 Ensembl 索引对象：{A.shape[0]} 核 × {len(keep)} 基因')

Q = ad.AnnData(X=sp.csr_matrix(A.X[:, keep]))
Q.obs = A.obs[['patient_id', 'sample_id', 'stage']].copy()
Q.var.index = eid_keep
n_before = len(Q.var_names)
Q.var_names_make_unique()
if len(set(Q.var_names)) != n_before:
    log(f'  [warn] var_names_make_unique 改了 {n_before - len(set(Q.var_names))} 个（Ensembl 重复）')
del A
log(f'  完成：{Q.shape};stage 分布 {dict(Q.obs["stage"].value_counts())}')

# ——— 2. 载入编码器 ———
log('载入 CellEmbedder（model.pt + best_state_dict.pth）...')
model = torch.load(MODEL_PT, map_location=torch.device(P_DEVICE), weights_only=False)
model.load_state_dict(torch.load(STATE_PTH, map_location=torch.device(P_DEVICE), weights_only=False))
model.to(P_DEVICE)
model.eval()
log(f'  参数量 {sum(p.numel() for p in model.parameters())/1e6:.2f} M；device={P_DEVICE}')

# ——— 3. 零样本嵌入 ———
log(f'embed_adata（batch_size={P_BATCH_SIZE}, pre_normalized={P_PRE_NORMALIZED}）...')
embed_adata(model, Q, batch_size=P_BATCH_SIZE, pre_normalized=P_PRE_NORMALIZED)
latent = Q.obsm['X_ce_latent']
log(f'  X_ce_latent {latent.shape} {latent.dtype}；有限值 {np.isfinite(latent).all()}')
np.save(f'{D}/scmg_latent.npy', latent)
log(f'  存盘 {D}/scmg_latent.npy')
del model

# ——— 4. 投影到全局流形 ———
log('载入参照全局流形 ...')
adata_ref = ad.read_h5ad(REF_H5AD)
log(f'  {adata_ref.shape[0]} 参照细胞 × {adata_ref.shape[1]} 基因；cell_type {adata_ref.obs["cell_type"].nunique()} 类')
assert P_EMB_KEY in adata_ref.obsm, f'参照缺 obsm[{P_EMB_KEY}]'

log(f'CellTypeSearcher.search_ref_cell（n_jobs={P_N_JOBS}）...')
cts = CellTypeSearcher(adata_ref, emb_key=P_EMB_KEY)
match = cts.search_ref_cell(latent, project_umap=True, n_jobs=P_N_JOBS)
log(f'  最近参照细胞 {len(match)} 行；距离中位 {match.distance.median():.3f}')

proj = pd.DataFrame({
    'cell': Q.obs_names,
    'patient_id': Q.obs['patient_id'].values,
    'sample_id': Q.obs['sample_id'].values,
    'stage': Q.obs['stage'].values,
    'ref_cell': match['ref_cell'].values,
    'distance': match['distance'].values,
    'umap_x': match['umap_x'].values,
    'umap_y': match['umap_y'].values,
})
proj['projected_cell_type'] = adata_ref.obs['cell_type'].loc[proj['ref_cell']].values
proj['projected_tissue'] = adata_ref.obs['tissue'].loc[proj['ref_cell']].values
proj['projected_major_cell_type'] = adata_ref.obs['major_cell_type'].loc[proj['ref_cell']].values
proj.to_csv(f'{D}/scmg_projection.csv.gz', index=False, compression='gzip')
log(f'  存盘 {D}/scmg_projection.csv.gz')

# 同时把参照流形的 UMAP + 标签存一份，给分析脚本画底图（避免重复载入 3 GB）
ref_small = pd.DataFrame({
    'ref_cell': adata_ref.obs_names,
    'umap_x': adata_ref.obsm['X_umap'][:, 0],
    'umap_y': adata_ref.obsm['X_umap'][:, 1],
    'cell_type': adata_ref.obs['cell_type'].values,
    'tissue': adata_ref.obs['tissue'].values,
    'major_cell_type': adata_ref.obs['major_cell_type'].values,
})
ref_small.to_csv(f'{D}/scmg_ref_manifold.csv.gz', index=False, compression='gzip')
log(f'  存盘 {D}/scmg_ref_manifold.csv.gz')

log('[done] 第一腿投影完成')
