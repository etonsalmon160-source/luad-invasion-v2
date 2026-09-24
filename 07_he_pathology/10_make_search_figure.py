#!/usr/bin/env python3
# 提示词搜索 · 出图（英文标注，本机无中文字体）
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/07_he_pathology/prompt_compare'
R = pd.read_csv(f'{D}/prompt_search_all_pairs.csv')
C = pd.read_csv(f'{D}/density_free_ceiling.csv')

AUX = json.load(open(f'{D}/search_aux.json')) if os.path.exists(f'{D}/search_aux.json') else None
D_BASE = 0.958
BEST = R.auc_test.max()

fig, ax = plt.subplots(2, 2, figsize=(15.5, 11.5))

# ——— (a) 1344 条候选 —— 只有极少数能越过"纯密度"线 ———
a = ax[0, 0]
a.scatter(R.dens_test, R.auc_test, s=14, c='#4477aa', alpha=.5, edgecolor='none')
a.axhline(D_BASE, color='#cc3311', lw=2.2, ls='--',
          label=f'prompt-free density baseline  AUC={D_BASE:.3f}')
a.axhline(0.891, color='#ee7733', lw=1.8, ls=':',
          label='current best real prompt (B, adenocarcinoma only)  AUC=0.891')
n_beat = int((R.auc_test > D_BASE).sum())
a.scatter(R[R.auc_test > D_BASE].dens_test, R[R.auc_test > D_BASE].auc_test,
          s=26, c='#cc3311', alpha=.9, edgecolor='white', linewidth=.4,
          label=f'beats density: {n_beat} / {len(R)}  ({n_beat/len(R):.1%})')
a.axhline(R.auc_test.median(), color='#777777', lw=1.4, ls='-.',
          label=f'median candidate  AUC={R.auc_test.median():.3f}')
a.set_xlabel('density-dependence   |corr(score, tissue fraction)|', fontsize=12)
a.set_ylabel('AUC ( LUAD spots vs Normal spots )', fontsize=12)
a.set_title('(a) 1,344 searched prompt pairs\n'
            'a prompt-free density measure already beats 98.4% of them',
            fontsize=13, fontweight='bold')
a.legend(fontsize=9, loc='lower left')
a.grid(alpha=.3); a.set_ylim(.1, 1.02)

# ——— (b) 诚实的阶梯 ———
a = ax[0, 1]
names = ['A set\n(7 naive terms,\nsummed lesion)', 'B set\n(adenocarcinoma\nterm alone)',
         'prompt-free\ndensity direction d', 'best of\n1,344 prompts',
         'best prompt,\ndensity removed']
vals = [0.806, 0.891, D_BASE, BEST, 0.697]
cols = ['#888888', '#1b9e77', '#cc3311', '#4477aa', '#aa3377']
b = a.bar(range(5), vals, color=cols, alpha=.9, width=.62)
for i, v in enumerate(vals):
    a.text(i, v + .012, f'{v:.3f}', ha='center', fontsize=12, fontweight='bold')
a.axhline(.5, color='k', lw=.8, ls='--', alpha=.6)
a.set_xticks(range(5)); a.set_xticklabels(names, fontsize=9)
a.set_ylabel('AUC ( LUAD spots vs Normal spots )', fontsize=12)
a.set_title('(b) The honest ladder\n'
            'the search\'s "win" is +0.013 over a measure that knows no biology',
            fontsize=13, fontweight='bold')
a.set_ylim(0, 1.08); a.grid(axis='y', alpha=.3)

# ——— (c) 去掉密度方向后四个候选全部塌向同一平台 ———
a = ax[1, 0]
x = np.arange(len(C)); w = .36
a.bar(x - w / 2, C.auc_raw, w, color='#4477aa', alpha=.9, label='as-is')
a.bar(x + w / 2, C.auc_density_projected, w, color='#cc3311', alpha=.9,
      label='after removing the density direction')
for i, r in C.iterrows():
    a.text(i - w / 2, r.auc_raw + .012, f'{r.auc_raw:.3f}', ha='center', fontsize=9.5)
    a.text(i + w / 2, r.auc_density_projected + .012, f'{r.auc_density_projected:.3f}',
           ha='center', fontsize=9.5, fontweight='bold')
a.axhspan(.73, .80, color='#cc3311', alpha=.10)
a.set_xticks(x)
a.set_xticklabels([n.replace('-', '\n') for n in C.name], fontsize=10)
a.set_ylabel('AUC ( LUAD vs Normal )', fontsize=12)
a.set_title('(c) Same plateau regardless of wording\n'
            'once density is projected out, every candidate lands in the red band (~0.73-0.80)',
            fontsize=12.5, fontweight='bold')
a.legend(fontsize=10, loc='upper center', ncol=2)
a.grid(axis='y', alpha=.3); a.set_ylim(0, 1.22)

# ——— (d) 按组织占比分层看着没塌，但那是分层太粗 ———
a = ax[1, 1]
def auc(y, s):
    y = np.asarray(y)
    r = pd.Series(np.asarray(s)).rank().values
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float('nan')
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

f = pd.read_csv(f'{D}/p4_prompt_compare_spots.csv.gz')
S = ['GSM9226174_P4_Normal', 'GSM9226175_P4_AAH',
     'GSM9226176_P4_AAH-1', 'GSM9226177_P4_LUAD']
meta, emb = None, []
for g in S:
    e = np.load(f'{D}/embeds_{g}.npy')
    sub = f[(f.slide == g) & (f.set == 'A_naive_terms+naive_template')]
    m = sub[['tissue_frac', 'stage']].copy()
    meta = m if meta is None else pd.concat([meta, m], ignore_index=True)
    emb.append(e)
E = np.vstack(emb); meta = meta.reset_index(drop=True)
tf = meta.tissue_frac.values; isL = (meta.stage == 'LUAD').values
lo, hi = np.quantile(tf, .25), np.quantile(tf, .75)
d = E[tf >= hi].mean(0) - E[tf <= lo].mean(0); d /= np.linalg.norm(d)

import torch
from transformers import CLIPModel, CLIPProcessor
mdl = CLIPModel.from_pretrained(f'{ROOT}/models/plip').eval()
prc = CLIPProcessor.from_pretrained(f'{ROOT}/models/plip')
sc = mdl.logit_scale.exp().item()
def enc(t):
    with torch.no_grad():
        e = mdl.get_text_features(**prc(text=[t], return_tensors='pt', padding=True))
    return (e / e.norm(dim=-1, keepdim=True))[0].numpy()
ea = enc('a histopathology image of invasive adenocarcinoma with desmoplastic stroma')
eb = enc('a histopathology image of morphologically normal pneumocytes lining empty alveoli')
lg = sc * np.vstack([E @ ea, E @ eb])
ex = np.exp(lg - lg.max(0, keepdims=True)); s = ex[0] / ex.sum(0)

qs = np.quantile(tf, np.linspace(0, 1, 6))
qp, qd = [], []
for i in range(5):
    m = (tf >= qs[i]) & (tf < qs[i + 1] if i < 4 else tf <= qs[i + 1])
    qp.append(auc(isL[m], s[m])); qd.append(auc(isL[m], tf[m]))
x = np.arange(5)
a.bar(x - w / 2, qp, w, color='#4477aa', alpha=.9, label='best prompt')
a.bar(x + w / 2, qd, w, color='#cc3311', alpha=.9, label='tissue fraction alone')
a.axhline(.5, color='k', lw=.9, ls='--', alpha=.6)
a.set_xticks(x)
a.set_xticklabels([f'Q{i+1}\n(tf {qs[i]:.2f}-{qs[i+1]:.2f})' for i in range(5)], fontsize=9)
a.set_xlabel('tissue-fraction quintile  (matching density within each band)', fontsize=11)
a.set_ylabel('AUC within the band', fontsize=12)
a.set_title('(d) Why matching tissue fraction is NOT enough\n'
            '"fraction of white pixels" is one scalar; the density direction is far richer',
            fontsize=13, fontweight='bold')
a.legend(fontsize=10); a.grid(axis='y', alpha=.3); a.set_ylim(0, 1.1)

fig.suptitle('Can a better prompt be found?  —  the search mostly re-invents the density detector\n'
             'patient P4 only  ·  16,000 spots  ·  crop 448px  ·  seed 0',
             fontsize=14.5, fontweight='bold', y=.985)
fig.tight_layout(rect=[0, 0, 1, .955])
out = f'{D}/prompt_search_figure.png'
fig.savefig(out, dpi=140, bbox_inches='tight', facecolor='white')
print(f'[out] {out}')
