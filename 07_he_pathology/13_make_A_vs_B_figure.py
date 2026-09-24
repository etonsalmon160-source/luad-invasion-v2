#!/usr/bin/env python3
# A 套 vs B 套：在**两套都跑完的切片**上逐张对比（英文标注，本机无中文字体）
# B 还在跑时也能出图 —— 只比较共同完成的部分，切片数写在标题里。
import glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

ROOT = '/home/eto/luad_v2'
DA = f'{ROOT}/results/07_he_pathology/spot_annotation'
DB = f'{ROOT}/results/07_he_pathology/spot_annotation_B'
ORDER = ['Normal', 'AAH', 'AIS', 'MIA', 'LUAD']
COL = {'Normal': '#1b9e77', 'AAH': '#7570b3', 'AIS': '#d95f02',
       'MIA': '#e7298a', 'LUAD': '#cc3311'}


def stage_of(g):
    s = g.split('_')[2]
    return 'AAH' if s.startswith('AAH') else s      # P4_AAH-1 归入 AAH


def load(d, cols):
    out = {}
    for f in glob.glob(f'{d}/per_slide/*.csv.gz'):
        g = os.path.basename(f)[:-7]
        x = pd.read_csv(f, usecols=cols)
        out[g] = {c: x[c].mean() for c in cols}
    return out


A = load(DA, ['neoplastic_score', 'tissue_frac'])
B = load(DB, ['lesion_score', 'p0', 'tissue_frac'])
common = sorted(set(A) & set(B), key=lambda g: (ORDER.index(stage_of(g)), g))
S = pd.DataFrame([{'slide': g, 'stage': stage_of(g),
                   'A_sum': A[g]['neoplastic_score'],
                   'B_sum': B[g]['lesion_score'],
                   'B_adeno': B[g]['p0']} for g in common])

r_sum, _ = pearsonr(S.A_sum, S.B_sum)
rho_sum, _ = spearmanr(S.A_sum, S.B_sum)
r_ad, _ = pearsonr(S.A_sum, S.B_adeno)

fig, ax = plt.subplots(1, 3, figsize=(19, 6.2))

# ——— (a) 逐张配对：A 的合成分 → B 的合成分 ———
a = ax[0]
for _, r in S.iterrows():
    a.plot([0, 1], [r.A_sum, r.B_sum], '-', color=COL[r.stage], lw=1.6, alpha=.75, zorder=2)
    a.scatter([0, 1], [r.A_sum, r.B_sum], s=26, color=COL[r.stage], alpha=.9,
              edgecolor='white', linewidth=.5, zorder=3)
for st in ORDER:
    if (S.stage == st).any():
        a.scatter([], [], color=COL[st], label=st, s=40)
a.set_xticks([0, 1]); a.set_xticklabels(['A set\n(summed lesion)', 'B set\n(summed lesion)'], fontsize=11)
a.set_ylabel('mean lesion score per slide', fontsize=12)
a.set_title(f'(a) Paired, slide by slide   (n={len(S)})\n'
            f'Pearson r = {r_sum:.3f}   Spearman = {rho_sum:.3f}',
            fontsize=13, fontweight='bold')
a.legend(fontsize=10, title='stage', title_fontsize=10)
a.grid(alpha=.3)

# ——— (b) 分期均值：A 合成 / B 合成 / B 单用腺癌条 ———
a = ax[1]
t = S.groupby('stage')[['A_sum', 'B_sum', 'B_adeno']].mean().reindex(ORDER)
x = np.arange(len(ORDER)); w = .26
for k, (c, lb, col) in enumerate([('A_sum', 'A set, summed', '#888888'),
                                  ('B_sum', 'B set, summed', '#cc3311'),
                                  ('B_adeno', 'B set, adenocarcinoma term alone', '#1b9e77')]):
    v = t[c].values
    a.bar(x + (k - 1) * w, np.nan_to_num(v), w, color=col, alpha=.9, label=lb)
    for i, y in enumerate(v):
        if not np.isnan(y):
            a.text(x[i] + (k - 1) * w, y + .012, f'{y:.2f}', ha='center', fontsize=8.5)
a.set_xticks(x); a.set_xticklabels(ORDER, fontsize=12)
a.set_ylabel('mean score per slide', fontsize=12)
a.set_title('(b) Summed vs single-term\n'
            'B\'s advantage lives in its single adenocarcinoma term, not in the sum',
            fontsize=12.5, fontweight='bold')
a.legend(fontsize=9, loc='upper left'); a.grid(axis='y', alpha=.3); a.set_ylim(0, 1.0)
for i, st in enumerate(ORDER):
    if st in ('Normal',) or i == 0:
        a.text(i, .955, f'n={int((S.stage==st).sum())}', ha='center', fontsize=8.5, color='#555555')

# ——— (c) A 合成分 vs B 单条：是否真的换了序 ———
a = ax[2]
for st in ORDER:
    v = S[S.stage == st]
    if len(v):
        a.scatter(v.A_sum, v.B_adeno, s=85, color=COL[st], alpha=.85,
                  edgecolor='white', linewidth=1, label=st, zorder=3)
lim = [0, max(S.A_sum.max(), S.B_adeno.max()) * 1.12]
a.plot(lim, lim, 'k--', lw=1.5, alpha=.7, zorder=1, label='y = x')
a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel('A set, summed lesion score', fontsize=12)
a.set_ylabel('B set, adenocarcinoma term alone', fontsize=12)
a.set_title(f'(c) Does B re-rank the slides?\nsame-axis scatter, Pearson r = {r_ad:.3f}',
            fontsize=13, fontweight='bold')
a.legend(fontsize=10, loc='upper left'); a.grid(alpha=.3)

fig.suptitle(f'PLIP prompt set A vs prompt set B — on the {len(S)} slides where both are done\n'
             f'A is complete (56/56); B is still running (14/56) — this figure only compares the overlap',
             fontsize=14, fontweight='bold', y=.99)
fig.tight_layout(rect=[0, 0, 1, .90])
out = f'{DB}/A_vs_B_figure.png'
fig.savefig(out, dpi=140, bbox_inches='tight', facecolor='white')
print(f'[out] {out}')
print(f'共同 {len(S)} 张；A合成 vs B合成 r={r_sum:.3f}；A合成 vs B单条 r={r_ad:.3f}')
