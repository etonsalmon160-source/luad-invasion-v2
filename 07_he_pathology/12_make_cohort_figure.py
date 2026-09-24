#!/usr/bin/env python3
# 04 全队列结果出图（英文标注，本机无中文字体）
# 目的：把「56 张切片 / 639,650 spot」这轮跑完的东西摊开看，包括它**没**成立的部分。
import os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/07_he_pathology/spot_annotation'
ORDER = ['Normal', 'AAH', 'AIS', 'MIA', 'LUAD']
COL = {'Normal': '#1b9e77', 'AAH': '#7570b3', 'AIS': '#d95f02',
       'MIA': '#e7298a', 'LUAD': '#cc3311'}

# ——— 读全量逐 spot 表（只取需要的列）———
parts = []
for f in sorted(glob.glob(f'{D}/per_slide/*.csv.gz')):
    parts.append(pd.read_csv(f, usecols=['gsm', 'slide_stage', 'neoplastic_score',
                                         'tissue_frac', 'argmax_label']))
A = pd.concat(parts, ignore_index=True)
S = A.groupby(['gsm', 'slide_stage'], as_index=False).agg(
        neo=('neoplastic_score', 'mean'), tf=('tissue_frac', 'mean'),
        n=('neoplastic_score', 'size'))
S.columns = ['gsm', 'stage', 'neo', 'tf', 'n']

rho_s, p_s = spearmanr(S.stage.map({v: i for i, v in enumerate(ORDER)}), S.neo)
r_tf, p_tf = pearsonr(S.tf, S.neo)

fig, ax = plt.subplots(2, 2, figsize=(15.5, 11.5))

# ——— (a) 逐切片：随分期单调 ———
a = ax[0, 0]
rng = np.random.default_rng(0)
for i, st in enumerate(ORDER):
    v = S[S.stage == st].neo.values
    a.scatter(i + rng.uniform(-.16, .16, len(v)), v, s=np.clip(S[S.stage == st].n / 220, 12, 150),
              c=COL[st], alpha=.75, edgecolor='white', linewidth=.7, zorder=3)
    if len(v) > 1:
        a.hlines(np.median(v), i - .28, i + .28, color='k', lw=2.6, zorder=4)
        a.vlines(i, np.percentile(v, 25), np.percentile(v, 75), color='k', lw=1.2, zorder=4)
a.set_xticks(range(5)); a.set_xticklabels(ORDER, fontsize=12)
a.set_ylabel('mean neoplastic score per slide', fontsize=12)
a.set_title(f'(a) The cohort-level signal is real\n'
            f'56 slides, 25 patients  ·  Spearman rho = {rho_s:.3f}  (p = {p_s:.1e})',
            fontsize=13, fontweight='bold')
a.grid(alpha=.3); a.set_ylim(0, .85)
a.annotate('only 1 Normal slide', xy=(0, .196), xytext=(.45, .40), fontsize=10, color='#1b9e77',
           arrowprops=dict(arrowstyle='->', color='#1b9e77', lw=1.4))

# ——— (b) 逐 spot：分布重叠极严重 ———
a = ax[0, 1]
data = [A[A.slide_stage == st].neoplastic_score.values for st in ORDER]
vp = a.violinplot(data, positions=range(5), widths=.78, showextrema=False, showmedians=True)
for j, b in enumerate(vp['bodies']):
    b.set_facecolor(COL[ORDER[j]]); b.set_alpha(.55); b.set_edgecolor('k'); b.set_linewidth(.6)
vp['cmedians'].set_color('k'); vp['cmedians'].set_linewidth(2)
for i, st in enumerate(ORDER):
    v = data[i]
    a.text(i, .97, f'n={len(v):,}', ha='center', fontsize=9.5, color=COL[st], fontweight='bold')
a.set_xticks(range(5)); a.set_xticklabels(ORDER, fontsize=12)
a.set_ylabel('per-spot neoplastic score', fontsize=12)
a.set_title('(b) ...but per spot the classes overlap massively\n'
            'the score shifts the bulk, it does not separate spot populations',
            fontsize=13, fontweight='bold')
a.grid(axis='y', alpha=.3); a.set_ylim(0, 1.02)

# ——— (c) 诚实面板：这仍然是"这块组织多实" ———
a = ax[1, 0]
for st in ORDER:
    v = S[S.stage == st]
    a.scatter(v.tf, v.neo, s=np.clip(v.n / 220, 22, 190), c=COL[st],
              alpha=.8, edgecolor='white', linewidth=.8, label=st, zorder=3)
z = np.polyfit(S.tf, S.neo, 1)
xs = np.linspace(S.tf.min(), S.tf.max(), 50)
a.plot(xs, np.polyval(z, xs), 'k--', lw=2, alpha=.75, zorder=4)
a.set_xlabel('mean tissue fraction per slide   (how packed the spot is)', fontsize=11.5)
a.set_ylabel('mean neoplastic score per slide', fontsize=12)
a.set_title(f'(c) The confound, at full-cohort scale\n'
            f'r = {r_tf:.3f} between "neoplastic" and tissue packing (p = {p_tf:.1e})',
            fontsize=13, fontweight='bold')
a.legend(fontsize=9.5, title='stage', title_fontsize=9.5, loc='upper left')
a.grid(alpha=.3)

# ——— (d) argmax 仍然判死 ———
a = ax[1, 1]
tab = (A.groupby(['slide_stage', 'argmax_label']).size()
       .unstack(fill_value=0).reindex(ORDER))
tab = tab.div(tab.sum(1), axis=0)
short = {'lung adenocarcinoma': 'adenocarcinoma', 'normal lung tissue': 'normal lung',
         'lung tissue with atypical adenomatous hyperplasia': 'AAH',
         'lung tissue with adenocarcinoma in situ': 'AIS',
         'fibrous stroma': 'fibrous stroma', 'lymphoid tissue': 'lymphoid',
         'an empty glass slide': 'empty slide'}
cols = list(tab.columns)
pal = ['#cc3311', '#1b9e77', '#7570b3', '#d95f02', '#999999', '#4477aa', '#dddddd']
bot = np.zeros(5)
for c, col in zip(cols, pal):
    v = tab[c].values
    a.bar(range(5), v, bottom=bot, color=col, alpha=.9, width=.66,
          edgecolor='white', linewidth=.7, label=short.get(c, c))
    bot += v
for i, st in enumerate(ORDER):
    ad = tab.loc[st].get('lung adenocarcinoma', 0)
    a.text(i, 1.015, f'adeno {ad:.1%}', ha='center', fontsize=9.5,
           color='#cc3311', fontweight='bold')
a.set_xticks(range(5)); a.set_xticklabels(ORDER, fontsize=12)
a.set_ylabel('fraction of spots by argmax label', fontsize=12)
a.set_ylim(0, 1.09)
a.set_title('(d) argmax is still unusable, cohort-wide\n'
            'true adenocarcinoma is almost never the winner; "AAH" takes the LUAD spots',
            fontsize=12.5, fontweight='bold')
a.legend(fontsize=8.5, loc='lower center', ncol=2, framealpha=.92)

fig.suptitle('PLIP on H&E — full cohort (56 slides, 639,650 spots, 442.9 min, 0 failures)\n'
             'a monotone stage signal that is largely a tissue-packing signal',
             fontsize=14.5, fontweight='bold', y=.985)
fig.tight_layout(rect=[0, 0, 1, .945])
out = f'{D}/cohort_figure.png'
fig.savefig(out, dpi=140, bbox_inches='tight', facecolor='white')
print(f'[out] {out}')
print(f'rho(slide-level) = {rho_s:.3f}   r(neo vs tissue_frac) = {r_tf:.3f}')
