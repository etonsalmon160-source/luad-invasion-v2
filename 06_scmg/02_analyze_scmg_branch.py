#!/usr/bin/env python3
# SCMG 分支 · 第一腿 · 判读
# 判读**写死在 results/06_scmg/SCMG_BRANCH_PREREG.md §五**，本脚本只执行。
#
# 跑法：python3 06_scmg/02_analyze_scmg_branch.py

import os
import json
import numpy as np
import pandas as pd

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/06_scmg'
FIGDIR = f'{D}/figures'
os.makedirs(FIGDIR, exist_ok=True)

# ——— 预注册冻结常量（§三 / §五）———
COVERAGE_MIN = 0.85
N_GENES_IN = 18069
N_CELLS_EXPECTED = 413697
# §五 S2：肺相关 tissue 集（参照中占比 15,823/133,061 = 11.89%）
LUNG_TISSUES = {'Lung epithelium', 'lung', 'lung parenchyma', 'respiratory airway',
                'lingula of left lung', 'Gut and lung epithelium'}
S2_MIN_FRACTION = 0.10
# §五 S3
AUROC_STRONG = 0.80
AUROC_WEAK = 0.60
NEG_STAGE, POS_STAGE = 'Normal', 'IAC'
FIG_SEED = 0
FIG_N_REF = 20000
FIG_N_PER_STAGE = 4000

def auroc(pos, neg):
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    s = np.r_[pos, neg]
    r = pd.Series(s).rank(method='average').values
    n1, n0 = y.sum(), (1 - y).sum()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

# ——— 读入 ———
latent = np.load(f'{D}/scmg_latent.npy')
proj = pd.read_csv(f'{D}/scmg_projection.csv.gz')
ref = pd.read_csv(f'{D}/scmg_ref_manifold.csv.gz')
assert latent.shape[0] == len(proj) == N_CELLS_EXPECTED, (latent.shape, len(proj))
print(f'[in] 潜空间 {latent.shape}；投影 {len(proj)} 行；参照 {len(ref)} 行')

# ——— S1 仪器有效：覆盖率 / 有限性 / 非退化 ———
n_mapped = N_GENES_IN - sum(1 for _ in open(f'{D}/scmg_unmapped_genes.txt'))
coverage = n_mapped / N_GENES_IN
finite = bool(np.isfinite(latent).all())
dim_std = latent.std(axis=0)
nondeg = bool((dim_std > 0).all()) and bool(len(np.unique(latent[:100], axis=0)) > 1)
s1 = (coverage >= COVERAGE_MIN) and finite and nondeg
print(f'\n[S1] 覆盖率 {coverage*100:.2f}% (≥{COVERAGE_MIN*100:.0f}%) → {"过" if coverage>=COVERAGE_MIN else "不过"}')
print(f'     潜空间全有限 {finite}；逐维标准差最小 {dim_std.min():.4g} → 非退化 {nondeg}')
print(f'     S1 总判：{"过" if s1 else "不过"}')

# ——— S2 参照合理性：最近参照细胞的肺相关 tissue 占比 ———
is_lung = proj['projected_tissue'].isin(LUNG_TISSUES)
frac_lung = float(is_lung.mean())
s2 = frac_lung >= S2_MIN_FRACTION
print(f'\n[S2] 最近参照细胞落在肺相关 tissue 的占比 {frac_lung*100:.2f}% (阈值 ≥{S2_MIN_FRACTION*100:.0f}%) → {"过" if s2 else "不过"}')
print(f'     （参照流形本身的肺占比基准 = 11.89%）')
print('     各组占比：')
print(proj.assign(lung=is_lung).groupby('stage')['lung'].mean().mul(100).round(2).to_string())

if not s2:
    print('\n🔴 S2 不过 ⇒ 预注册 §五：报「投影未落在肺区」，本腿作废，不得解读任何下游分布。')

# ——— S3 分期可分性：留一患者最近质心 Normal vs IAC ———
mask = proj['stage'].isin([NEG_STAGE, POS_STAGE]).values
X = latent[mask].astype(np.float64)
pid = proj.loc[mask, 'patient_id'].values
y = (proj.loc[mask, 'stage'].values == POS_STAGE).astype(int)
print(f'\n[S3] 参与细胞 {len(X)}（{NEG_STAGE} {(y==0).sum()} / {POS_STAGE} {(y==1).sum()}）；患者 {len(np.unique(pid))}')

tot1 = X[y == 1].sum(axis=0); n1_tot = (y == 1).sum()
tot0 = X[y == 0].sum(axis=0); n0_tot = (y == 0).sum()

scores, ys, pids = [], [], []
per_patient = []
for p in np.unique(pid):
    te = pid == p
    if te.sum() == 0:
        continue
    tr = ~te
    # 训练侧质心（该患者排除在外）
    n1_tr = (y[tr] == 1).sum(); n0_tr = (y[tr] == 0).sum()
    if n1_tr == 0 or n0_tr == 0:
        continue
    pm1 = X[te & (y == 1)].sum(axis=0); pm0 = X[te & (y == 0)].sum(axis=0)
    c1 = (tot1 - pm1) / n1_tr
    c0 = (tot0 - pm0) / n0_tr
    d1 = np.linalg.norm(X[te] - c1, axis=1)
    d0 = np.linalg.norm(X[te] - c0, axis=1)
    s = d0 - d1                      # 越大越像 IAC
    scores.append(s); ys.append(y[te]); pids.append(np.repeat(p, te.sum()))
    if (y[te] == 1).sum() > 0 and (y[te] == 0).sum() > 0:
        per_patient.append((p, auroc(s[y[te] == 1], s[y[te] == 0]), int((y[te] == 1).sum()), int((y[te] == 0).sum())))

scores = np.concatenate(scores); ys = np.concatenate(ys); pids = np.concatenate(pids)
a = auroc(scores[ys == 1], scores[ys == 0])
v3 = ('分期信号被保留' if a >= AUROC_STRONG else
      '弱保留' if a >= AUROC_WEAK else '未被保留（全糊一块）')
print(f'     留一患者 pooled AUROC = {a:.4f} → {v3}')
print('\n     逐患者（有双类的）AUROC：')
for p, ap, npos, nneg in sorted(per_patient):
    print(f'       {p:12s} n_IAC={npos:6d} n_Normal={nneg:6d}  AUROC={ap:.3f}')

# ——— 描述性：各期映射到的参照细胞类型构成 ———
print('\n[描述] 各期 top-8 参照细胞类型占比（%）：')
comp = {}
for stg in ['Normal', 'AAH', 'AIS', 'MIA', 'IAC']:
    vc = proj.loc[proj['stage'] == stg, 'projected_cell_type'].value_counts(normalize=True).mul(100)
    comp[stg] = vc.head(8).round(2).to_dict()
    print(f'  --- {stg} (n={int((proj["stage"]==stg).sum())})')
    for k, v in comp[stg].items():
        print(f'      {v:6.2f}%  {k}')

# ——— 图（图内标签全英文）———
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rng = np.random.RandomState(FIG_SEED)
ri = rng.choice(len(ref), min(FIG_N_REF, len(ref)), replace=False)
fig, axes = plt.subplots(1, 2, figsize=(19, 8.2), dpi=140)
ax = axes[0]
ax.scatter(ref['umap_x'].values[ri], ref['umap_y'].values[ri], s=1.0, c='lightgrey', alpha=.5,
           label=f'SCMG reference manifold (n={len(ri):,})', rasterized=True)
ax.set_title('SCMG global manifold (reference) with our nuclei projected\n'
             '(each dot = nearest reference cell of one of our nuclei)')
ax.set_xlabel('SCMG UMAP 1'); ax.set_ylabel('SCMG UMAP 2')
colors = {'Normal': '#1f77b4', 'AAH': '#2ca02c', 'AIS': '#ff7f0e', 'MIA': '#d62728', 'IAC': '#8b0000'}
for stg in ['IAC', 'MIA', 'AIS', 'AAH', 'Normal']:
    m = (proj['stage'] == stg).values
    k = rng.choice(m.sum(), min(FIG_N_PER_STAGE, m.sum()), replace=False)
    ax.scatter(proj['umap_x'].values[m][k], proj['umap_y'].values[m][k], s=1.6,
               c=colors[stg], alpha=.75, label=f'{stg} (n={m.sum():,})', rasterized=True)
ax.legend(loc='upper right', fontsize=8, markerscale=4, framealpha=.9)

ax = axes[1]
for stg in ['IAC', 'MIA', 'AIS', 'AAH', 'Normal']:
    m = (proj['stage'] == stg).values
    k = rng.choice(m.sum(), min(FIG_N_PER_STAGE, m.sum()), replace=False)
    ax.scatter(proj['umap_x'].values[m][k], proj['umap_y'].values[m][k], s=1.6,
               c=colors[stg], alpha=.8, label=stg, rasterized=True)
ax.set_title(f'Our nuclei on the SCMG global manifold (subsampled {FIG_N_PER_STAGE:,}/stage)\n'
             f'S2 lung fraction = {frac_lung*100:.1f}%   |   S3 AUROC = {a:.3f}')
ax.set_xlabel('SCMG UMAP 1'); ax.set_ylabel('SCMG UMAP 2')
ax.legend(loc='upper right', fontsize=8, markerscale=4, framealpha=.9)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/scmg_1_global_manifold_projection.png', bbox_inches='tight')
plt.close()
print(f'\n[out] {FIGDIR}/scmg_1_global_manifold_projection.png')

# 图 2：各期参照细胞类型构成
stages = ['Normal', 'AAH', 'AIS', 'MIA', 'IAC']
top = proj['projected_cell_type'].value_counts().head(12).index.tolist()
fig, ax = plt.subplots(figsize=(13, 7.5), dpi=140)
bottom = np.zeros(len(stages))
cmap = plt.get_cmap('tab20')
for i, ct in enumerate(top):
    vals = np.array([100 * (proj.loc[proj['stage'] == s, 'projected_cell_type'] == ct).mean() for s in stages])
    ax.bar(stages, vals, bottom=bottom, label=ct[:44], color=cmap(i % 20))
    bottom += vals
ax.bar(stages, 100 - bottom, bottom=bottom, label='other', color='lightgrey')
ax.set_ylabel('share of nuclei (%)'); ax.set_ylim(0, 100)
ax.set_title('Reference cell type that our nuclei project onto, by stage\n'
             '(top-12 projected reference cell types; higher = more nuclei of that stage map there)')
ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/scmg_2_stage_composition.png', bbox_inches='tight')
plt.close()
print(f'[out] {FIGDIR}/scmg_2_stage_composition.png')

# ——— summary ———
summary = {
    'script': '06_scmg/02_analyze_scmg_branch.py',
    'prereg': 'results/06_scmg/SCMG_BRANCH_PREREG.md',
    'inputs': {
        'query_h5ad_sha256': 'a276cd1af69a4620ff9e3d64f4b93f33461e646537c5f2c9da1aff8c9150b9de',
        'embedder_model_pt_sha256': '01bd51e6d8a54486e5d3bc412c980466031c5f0e91c8ec2fe0474c4c3150f599',
        'embedder_state_sha256': 'e585faae3fb2f7805880cebdfcb18c9550a636381e8af211ee3bfcab3742119f',
        'standard_genes_sha256': '44c436d39ea978a4787132a1ddd92469aa3fbdef793b69d99f15f7c5dd09de78',
        'ref_manifold_sha256': '27669cfbf96a7274c8bf2a6ea23e16b8f744ecc22143229f8c4b34f5f3133d4b',
        'scmg_git_commit': '29f44c98c5621d575a058b549e0fdde0ba31a730'},
    'frozen_params': {'device': 'cpu', 'batch_size': 8192, 'pre_normalized': False,
                      'emb_key': 'X_scmg', 'n_jobs': 20},
    'gene_coverage': coverage, 'n_genes_mapped': int(n_mapped), 'n_genes_unmapped': int(N_GENES_IN - n_mapped),
    'latent_shape': list(latent.shape),
    'S1_coverage_pass': bool(coverage >= COVERAGE_MIN), 'S1_finite': finite,
    'S1_non_degenerate': nondeg, 'S1_pass': bool(s1),
    'S2_lung_fraction': frac_lung, 'S2_threshold': S2_MIN_FRACTION, 'S2_pass': bool(s2),
    'S2_lung_fraction_by_stage': proj.assign(l=is_lung).groupby('stage')['l'].mean().round(4).to_dict(),
    'S3_auroc': float(a), 'S3_strong': AUROC_STRONG, 'S3_weak': AUROC_WEAK, 'S3_verdict': v3,
    'S3_per_patient': [{'patient': p, 'auroc': ap, 'n_iac': npos, 'n_normal': nneg}
                       for p, ap, npos, nneg in sorted(per_patient)],
    'reference_cell_type_composition_top8_by_stage': comp,
    'figure_subsample': {'seed': FIG_SEED, 'n_ref': FIG_N_REF, 'n_per_stage': FIG_N_PER_STAGE,
                         'note': '仅影响图，不影响任何门槛'},
    'caveats': [
        'snRNA（核）投到以 scRNA（全细胞）为主的参照流形；SCMG 与 2026 基准均未在 snRNA 上验证',
        '参照流形无肿瘤态 ⇒ SCMG 不能给恶性标签，只能给"最像的正常/发育细胞状态"',
        '参照是全身泛组织图谱（331 组织），肺仅占 11.89%',
        '分期↔样本嵌套；S3 用留一患者缓解，但同患者内仍共享技术批次',
        '单数据集无 ground truth ⇒ 本腿不关闭 GP9（GP9 须与 M3-A 做 scIB 对照）',
        'SCMG 的 797 类全局状态标签 ≠ LUAD 亚型标签，报告不得混用',
        '偏离 pyproject 固定版本（用主环境现有版本），未改任何数值口径',
        '本腿不是「迁徙分析」；迁徙需亚型边界，属未签字的 R4 边界问题'],
}
with open(f'{D}/scmg_branch_summary.json', 'w') as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
print(f'[out] {D}/scmg_branch_summary.json')
print('\n[done]')
