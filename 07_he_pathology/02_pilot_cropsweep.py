#!/usr/bin/env python3
# HE 病理模型 · 试点第二步：裁剪尺度扫描
#
# 第一步 224px 裁窗（≈1.27mm 视野）分不开 LUAD / Normal。
# 裁多大是**我的选择**，不是数据属性 ⇒ 扫一遍确认不是裁错。
# 扫完若仍无分离，结论是"分辨率不够"，不是"参数没调好"。
#
# 跑法：python3 07_he_pathology/02_pilot_cropsweep.py

import os, json, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import torch
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
from transformers import CLIPModel, CLIPProcessor

ROOT = '/home/eto/luad_v2'
MODEL_DIR = f'{ROOT}/models/plip'
EXT = '/home/eto/luad_invasion/data/GSE307534/extracted'
OUT = f'{ROOT}/results/07_he_pathology/pilot'
os.makedirs(OUT, exist_ok=True)

CROP_SIZES = [32, 64, 112, 224, 448]   # 原图像素，喂给 PLIP 前一律缩放到 224
N_SPOT = 300
SEED = 0
MIN_TISSUE_FRAC = 0.10
WHITE_THRESH = 230
PROMPT_TEMPLATE = 'a histopathology image of {}'
LABELS = [
    'lung adenocarcinoma',
    'normal lung tissue',
    'lung tissue with atypical adenomatous hyperplasia',
    'lung tissue with adenocarcinoma in situ',
    'fibrous stroma',
    'lymphoid tissue',
    'an empty glass slide',
]
TUMOR_LABEL, NORMAL_LABEL = 'lung adenocarcinoma', 'normal lung tissue'

SLIDES = [
    ('GSM9226174', 'P4_Normal', 'Normal'),
    ('GSM9226175', 'P4_AAH',    'AAH'),
    ('GSM9226177', 'P4_LUAD',   'LUAD'),
]


def slide_dir(gsm, name):
    """内层目录名不统一（P4_AAH1 / P4_AAH2 / P4_Normal），自动发现含 spatial/ 的那层。"""
    base = f'{EXT}/{gsm}_{name}'
    for inner in sorted(os.listdir(base)):
        d = f'{base}/{inner}/spatial'
        if os.path.isdir(d):
            return d
    raise FileNotFoundError(f'no spatial/ under {base}')


def load(gsm, name):
    d = slide_dir(gsm, name)
    img = Image.open(f'{d}/tissue_hires_image.png').convert('RGB')
    pf = f'{d}/tissue_positions.csv'
    if os.path.exists(pf):
        pos = pd.read_csv(pf)
    else:
        pos = pd.read_csv(f'{d}/tissue_positions_list.csv', header=None)
        pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                       'pxl_row_in_fullres', 'pxl_col_in_fullres']
    sc = json.load(open(f'{d}/scalefactors_json.json'))
    return img, pos, sc


def main():
    torch.set_num_threads(16)
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    texts = [PROMPT_TEMPLATE.format(l) for l in LABELS]

    # 预抽同一样本内固定的 spot 集合，保证跨 crop size 可比
    rng = np.random.RandomState(SEED)
    cache = {}
    for gsm, name, stage in SLIDES:
        img, pos, sf = load(gsm, name)
        sc = sf['tissue_hires_scalef']
        pos = pos[pos.in_tissue == 1].copy()
        pos['r'] = pos.pxl_row_in_fullres * sc
        pos['c'] = pos.pxl_col_in_fullres * sc
        take = pos.sample(min(N_SPOT, len(pos)), random_state=SEED)
        cache[stage] = (img, take, sc)
        print(f'[{stage}] {name} hires={img.size} spots={len(pos)} 抽样={len(take)}')

    rows = []
    for W in CROP_SIZES:
        print(f'\n--- crop {W}px ---')
        per_stage = {}
        for stage in ['Normal', 'AAH', 'LUAD']:
            img, take, sc = cache[stage]
            crops = []
            for _, x in take.iterrows():
                L, T = int(x.c) - W // 2, int(x.r) - W // 2
                L = max(0, min(L, img.size[0] - W)); T = max(0, min(T, img.size[1] - W))
                c = img.crop((L, T, L + W, T + W))
                if W != 224:
                    c = c.resize((224, 224), Image.BICUBIC)
                g = np.asarray(c.convert('L'))
                if (g <= WHITE_THRESH).mean() < MIN_TISSUE_FRAC:
                    continue
                crops.append(c)
            probs = []
            for i in range(0, len(crops), 64):
                inp = proc(text=texts, images=crops[i:i + 64], return_tensors='pt', padding=True)
                with torch.no_grad():
                    probs.append(model(**inp).logits_per_image.softmax(1).numpy())
            pr = np.vstack(probs)
            per_stage[stage] = pr
            argm = pr.argmax(1)
            print(f'  {stage:7s} n={len(crops):4d}  adenoca={pr[:,0].mean():.4f}  '
                  f'normal={pr[:,1].mean():.4f}  '
                  f'argmax: adenoca={np.mean(argm==0)*100:4.1f}% normal={np.mean(argm==1)*100:4.1f}% '
                  f'AAH={np.mean(argm==2)*100:4.1f}% empty={np.mean(argm==6)*100:4.1f}%')
            rows.append({'crop': W, 'stage': stage, 'n': len(crops),
                         'adenoca_mean': float(pr[:, 0].mean()),
                         'normal_mean': float(pr[:, 1].mean()),
                         'argmax_adenoca_frac': float(np.mean(argm == 0)),
                         'argmax_normal_frac': float(np.mean(argm == 1))})
        d_aden = per_stage['LUAD'][:, 0].mean() - per_stage['Normal'][:, 0].mean()
        d_norm = per_stage['LUAD'][:, 1].mean() - per_stage['Normal'][:, 1].mean()
        # 简单可分性：把 adenoca 分数当判别器，算 LUAD vs Normal 的 AUC
        pos = per_stage['LUAD'][:, 0]; neg = per_stage['Normal'][:, 0]
        y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
        r = pd.Series(np.r_[pos, neg]).rank().values
        n1, n0 = y.sum(), (1 - y).sum()
        auc = (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
        print(f'  >>> LUAD−Normal: adenoca {d_aden:+.4f}  normal {d_norm:+.4f}  '
              f'AUC(adenoca 分)={auc:.3f}')
        rows.append({'crop': W, 'stage': 'LUAD_MINUS_NORMAL', 'n': None,
                     'adenoca_mean': None, 'normal_mean': None,
                     'argmax_adenoca_frac': None, 'argmax_normal_frac': None,
                     'delta_adenoca': float(d_aden), 'delta_normal': float(d_norm),
                     'auc_adenoca': float(auc)})

    df = pd.DataFrame(rows)
    df.to_csv(f'{OUT}/plip_pilot_cropsweep.csv', index=False)
    print(f'\n[out] {OUT}/plip_pilot_cropsweep.csv')

    print('\n' + '=' * 60)
    print('汇总：AUC(adenoca 分数 区分 LUAD vs Normal)')
    sub = df[df.stage == 'LUAD_MINUS_NORMAL']
    for _, r in sub.iterrows():
        print(f"  crop {int(r['crop']):>4d}px  AUC={r['auc_adenoca']:.3f}  "
              f"Δadenoca={r['delta_adenoca']:+.4f}  Δnormal={r['delta_normal']:+.4f}")

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=140)
    W = CROP_SIZES
    for st, col in [('Normal', '#1f77b4'), ('AAH', '#2ca02c'), ('LUAD', '#d62728')]:
        v = [df[(df.crop == w) & (df.stage == st)]['adenoca_mean'].values[0] for w in W]
        axes[0].plot(W, v, 'o-', color=col, label=st)
    axes[0].set_xscale('log'); axes[0].set_xlabel('crop size (source px, resized to 224)')
    axes[0].set_ylabel('mean P("lung adenocarcinoma")'); axes[0].legend()
    axes[0].set_title('PLIP adenocarcinoma score vs crop size')
    axes[1].plot(W, [df[(df.crop == w) & (df.stage == 'LUAD_MINUS_NORMAL')]['auc_adenoca'].values[0] for w in W],
                 'o-', color='k')
    axes[1].axhline(.5, ls='--', c='grey'); axes[1].axhline(.8, ls=':', c='green')
    axes[1].set_xscale('log'); axes[1].set_xlabel('crop size (source px)')
    axes[1].set_ylabel('AUC (LUAD vs Normal)'); axes[1].set_ylim(0, 1)
    axes[1].set_title('Separability vs crop size (0.5 = chance)')
    plt.tight_layout()
    plt.savefig(f'{OUT}/figures/plip_pilot_cropsweep.png', bbox_inches='tight')
    plt.close()
    print(f'[out] {OUT}/figures/plip_pilot_cropsweep.png')


if __name__ == '__main__':
    main()
