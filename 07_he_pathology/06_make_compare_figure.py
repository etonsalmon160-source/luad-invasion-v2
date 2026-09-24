#!/usr/bin/env python3
# 提示词对照实验 · 出图（英文标注，本机无中文字体）
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/07_he_pathology/prompt_compare'
df = pd.read_csv(f'{D}/p4_prompt_compare_spots.csv.gz')

SETS = ['A_naive_terms+naive_template', 'C_naive_terms+PLIP_template',
        'D_naive_terms_lengthened', 'B_report_terms+PLIP_template']
SHORT = {'A_naive_terms+naive_template': 'A: naive terms\n+ naive template',
         'C_naive_terms+PLIP_template': 'C: naive terms\n+ PLIP template',
         'D_naive_terms_lengthened': 'D: naive terms\nLENGTHENED\n+ naive template',
         'B_report_terms+PLIP_template': 'B: report terms\n+ PLIP template'}
COL = {'A_naive_terms+naive_template': '#888888',
       'C_naive_terms+PLIP_template': '#d95f02',
       'D_naive_terms_lengthened': '#7570b3',
       'B_report_terms+PLIP_template': '#1b9e77'}
LC = json.load(open(f'{D}/length_control.json'))   # D 套（长度对照）
LAB = {
 'A_naive_terms+naive_template': ['lung adenocarcinoma', 'normal lung tissue', 'AAH', 'AIS',
                                  'fibrous stroma', 'lymphoid tissue', 'empty glass slide'],
 'C_naive_terms+PLIP_template': ['lung adenocarcinoma', 'normal lung tissue', 'AAH', 'AIS',
                                 'fibrous stroma', 'lymphoid tissue', 'empty glass slide'],
 'B_report_terms+PLIP_template': ['invasive adenoca +\ndesmoplastic stroma',
                                  'normal pneumocytes\nlining EMPTY alveoli',
                                  'thickened septa +\natypical type II',
                                  'neoplastic cells lining\nalveoli (lepidic)',
                                  'enlarged cuboidal\ntype II (RPII)',
                                  'fibrotic lung + collagen',
                                  'dense lymphocytic\ninfiltrate',
                                  'empty glass slide'],
}


def auc(y, s):
    y = np.asarray(y)
    r = pd.Series(np.asarray(s)).rank().values
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float('nan')
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


fig, ax = plt.subplots(2, 2, figsize=(16, 12))

# ——— (a) 密度捷径：lesion 分与「空玻璃片」提示词的相关 ———
a = ax[0, 0]
xs, ys, cs, errs, pers = [], [], [], [], []
for k in SETS:
    if k == 'D_naive_terms_lengthened':
        c_all = LC['corr_lesion_vs_empty']
        per = list(LC['corr_lesion_vs_empty_within_slide'].values())
    else:
        g = df[df.set == k]
        c_all = np.corrcoef(g.lesion_score, g.empty_score)[0, 1]
        per = [np.corrcoef(gg.lesion_score, gg.empty_score)[0, 1]
               for _, gg in g.groupby('slide')]
    xs.append(SHORT[k]); ys.append(c_all); cs.append(COL[k]); pers.append(per)
    errs.append([max(0.0, c_all - min(per)), max(0.0, max(per) - c_all)])
    a.scatter([len(xs) - 1] * len(per), per, s=55, color=COL[k],
              edgecolor='white', linewidth=1.2, zorder=3)
b = a.bar(range(len(SETS)), ys, color=cs, alpha=.75, width=.55, zorder=2)
a.errorbar(range(len(SETS)), ys, yerr=np.array(errs).T, fmt='none', ecolor='k',
           capsize=6, lw=1.3, zorder=4)
for i, v in enumerate(ys):
    a.text(i, v - .07 if v < 0 else v + .03, f'{v:+.2f}', ha='center',
           va='top' if v < 0 else 'bottom', fontsize=12, fontweight='bold')
a.axhline(0, color='k', lw=.8)
a.set_xticks(range(len(SETS))); a.set_xticklabels(xs, fontsize=9)
a.set_ylabel('corr( lesion score , "empty glass slide" probe )', fontsize=12)
a.set_title('(a) The density shortcut  —  it is the TERMS, not the length\n'
            'D matches B in sentence length but stays density-driven  '
            '(dots = per-slide)', fontsize=12, fontweight='bold')
a.set_ylim(-1.05, .45)
a.grid(axis='y', alpha=.3)
a.annotate('', xy=(3, .30), xytext=(2, .30),
           arrowprops=dict(arrowstyle='<->', color='#333', lw=1.4))
a.text(2.5, .335, 'length matched', ha='center', fontsize=9, style='italic', color='#333')

# ——— (b) 逐条提示词：判别力 vs 密度依赖 ———
a = ax[0, 1]
for k in ['A_naive_terms+naive_template', 'C_naive_terms+PLIP_template',
          'B_report_terms+PLIP_template']:
    g = df[df.set == k]
    isL = (g.stage == 'LUAD').astype(int).values
    isN = (g.stage == 'Normal').astype(int).values
    for j, lab in enumerate(LAB[k]):
        s = g[f'p{j}'].values
        x = np.corrcoef(s, g.tissue_frac)[0, 1]
        y = auc(isL, s)
        marker = 'X' if 'empty' in lab.lower() else 'o'
        a.scatter(x, y, s=150 if marker == 'o' else 190, color=COL[k],
                  marker=marker, edgecolor='white', linewidth=1.4, zorder=3,
                  alpha=.95)
        if k == 'A_naive_terms+naive_template':
            a.annotate(lab.replace('\n', ' '), (x, y), fontsize=8, color='#444',
                       xytext=(6, -11), textcoords='offset points')
        if k == 'B_report_terms+PLIP_template' and 'adenoca' in lab:
            a.annotate(lab.replace('\n', ' '), (x, y), fontsize=9, color=COL[k],
                       fontweight='bold', xytext=(-8, 9), textcoords='offset points',
                       ha='right', va='bottom')
        if k == 'B_report_terms+PLIP_template' and 'lepidic' in lab:
            a.annotate('neoplastic cells\nlining alveoli (lepidic)', (x, y), fontsize=8,
                       color=COL[k], xytext=(0, -13), textcoords='offset points',
                       ha='center', va='top')
a.axvline(0, color='k', lw=.9, ls='--', alpha=.6)
a.axhline(.5, color='k', lw=.9, ls='--', alpha=.6)
a.text(-.62, .93, 'REAL signal\n(separates, density-independent)', fontsize=9,
       ha='left', va='top', style='italic', color='#333')
a.text(.55, .10, 'density detector\n("empty" probes)', fontsize=9, ha='right',
       va='bottom', style='italic', color='#333')
a.set_xlabel('corr( probe probability , tissue fraction )   →  density-dependence', fontsize=12)
a.set_ylabel('AUC ( LUAD spots  vs  Normal spots )', fontsize=12)
a.set_title('(b) Every probe, by discriminating power and density-dependence\n'
            'X = "empty glass slide" probes', fontsize=13, fontweight='bold')
a.legend(handles=[Patch(color=COL[k], label=SHORT[k].replace('\n', ' '))
                  for k in ['A_naive_terms+naive_template',
                            'C_naive_terms+PLIP_template',
                            'B_report_terms+PLIP_template']],
         fontsize=9, loc='lower left')
a.grid(alpha=.3)
a.set_xlim(-.8, .95); a.set_ylim(0, 1.02)

# ——— (c) 分数谱：三套在 Normal/AAH/LUAD 上的 lesion 分 ———
a = ax[1, 0]
stg = ['Normal', 'AAH', 'LUAD']
w = .20
for i, k in enumerate(SETS):
    if k == 'D_naive_terms_lengthened':
        ms = [LC['mean_lesion_by_stage'][s] for s in stg]; es = [0, 0, 0]
    else:
        g = df[df.set == k]
        ms = [g[g.stage == s].lesion_score.mean() for s in stg]
        es = [g[g.stage == s].lesion_score.sem() for s in stg]
    a.bar(np.arange(3) + (i - 1.5) * w, ms, w, yerr=es, color=COL[k],
          label=SHORT[k].replace('\n', ' '), capsize=4, alpha=.9)
a.set_xticks(range(3)); a.set_xticklabels(stg, fontsize=12)
a.set_ylabel('mean lesion score  (adenoca + AAH + AIS probes)', fontsize=12)
a.set_title('(c) Does it order the progression?\n'
            'AAH should sit BETWEEN Normal and LUAD — in all three sets it does not',
            fontsize=13, fontweight='bold')
a.legend(fontsize=9)
a.grid(axis='y', alpha=.3)

# ——— (d) argmax 落在哪 ———
a = ax[1, 1]
namesA = ['adenocarcinoma', 'normal lung', 'AAH', 'AIS', 'fibrous stroma',
          'lymphoid', 'empty slide']
namesB = ['invasive adenoca', 'normal pneumocytes', 'thickened septa (AAH)',
          'lepidic (AIS)', 'RPII', 'fibrosis', 'lymphocytes', 'empty slide']


def argmax_frac(g, stage, labels):
    v = g[g.stage == stage].argmax_idx.value_counts(normalize=True)
    return np.array([v.get(i, 0.0) for i in range(len(labels))])


ypos = []
for bi, (k, names, stage) in enumerate([
        ('A_naive_terms+naive_template', namesA, 'Normal'),
        ('A_naive_terms+naive_template', namesA, 'LUAD'),
        ('B_report_terms+PLIP_template', namesB, 'Normal'),
        ('B_report_terms+PLIP_template', namesB, 'LUAD')]):
    g = df[df.set == k]
    fr = argmax_frac(g, stage, names)
    left = 0
    for i, (nm, fv) in enumerate(zip(names, fr)):
        if fv <= 0:
            continue
        a.barh(bi, fv, left=left, color=COL[k], alpha=.35 + .6 * (i / len(names)),
               edgecolor='white', height=.62)
        if fv > .07:
            a.text(left + fv / 2, bi, f'{nm}\n{fv:.0%}', ha='center', va='center',
                   fontsize=7.5, color='white' if fv > .25 else '#222', fontweight='bold')
        left += fv
    ypos.append(bi)
a.set_yticks(ypos)
a.set_yticklabels(['A  ·  Normal slide', 'A  ·  LUAD slide',
                   'B  ·  Normal slide', 'B  ·  LUAD slide'], fontsize=11)
a.set_xlim(0, 1); a.set_xlabel('fraction of spots', fontsize=12)
a.set_title('(d) Which class wins the vote (argmax)\n'
            'A: adenocarcinoma never wins · B: normal lung finally named correctly',
            fontsize=13, fontweight='bold')

fig.suptitle('PLIP prompt experiment on patient P4  ·  '
             '4,000 spots/slide, crop 448px, seed 0  ·  one patient only',
             fontsize=15, fontweight='bold', y=.985)
fig.tight_layout(rect=[0, 0, 1, .965])
out = f'{D}/prompt_compare_figure.png'
fig.savefig(out, dpi=145, bbox_inches='tight', facecolor='white')
print(f'[out] {out}')
