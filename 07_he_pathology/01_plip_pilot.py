#!/usr/bin/env python3
# HE 病理模型 · 试点：PLIP 零样本能不能分开「肿瘤 vs 正常」
#
# 这是**试点**，不是分析：目的是判定这条路走不走得通。
# 阳性对照是内置的 —— P4 同一患者的 Normal 切片 vs LUAD 切片。
# 若 PLIP 在 LUAD 片上给的 adenocarcinoma 分数不高于 Normal 片 ⇒ 路死，如实上报。
#
# 跑法：python3 07_he_pathology/01_plip_pilot.py

import os, sys, json, tarfile, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import torch
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
from transformers import CLIPModel, CLIPProcessor

ROOT = '/home/eto/luad_v2'
MODEL_DIR = f'{ROOT}/models/plip'
EXT = '/home/eto/luad_invasion/data/GSE307534/extracted'   # 只读源
OUT = f'{ROOT}/results/07_he_pathology/pilot'
FIGDIR = f'{OUT}/figures'
os.makedirs(FIGDIR, exist_ok=True)

# ——— 冻结常量 ———
CROP = 224                 # PLIP 原生输入
N_SPOT_PER_SLIDE = 400     # 每片随机抽样 spot 数
SEED = 0
MIN_TISSUE_FRAC = 0.10     # 裁窗内非白像素<10% 则丢弃（纯空腔）
WHITE_THRESH = 230         # 灰度>230 视为背景白
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
# 主判读：tumor 分数
TUMOR_LABEL = 'lung adenocarcinoma'
NORMAL_LABEL = 'normal lung tissue'

# 三个同患者切片（GSM, 期望分期）
SLIDES = [
    ('GSM9226174', 'P4_Normal', 'Normal'),
    ('GSM9226176', 'P4_AAH-1',  'AAH'),
    ('GSM9226177', 'P4_LUAD',   'LUAD'),
]


def load_slide(gsm, name):
    """读已解压的 Visium 目录（只读；数据在旧项目 luad_invasion 下，不写入）。"""
    d = f'{EXT}/{gsm}_{name}/{name}/spatial'
    if not os.path.isdir(d):
        return None
    try:
        img = Image.open(f'{d}/tissue_hires_image.png').convert('RGB')
        pf = f'{d}/tissue_positions.csv'
        if os.path.exists(pf):
            pos = pd.read_csv(pf)
        else:
            pos = pd.read_csv(f'{d}/tissue_positions_list.csv', header=None)
            pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                           'pxl_row_in_fullres', 'pxl_col_in_fullres']
        sc = json.load(open(f'{d}/scalefactors_json.json'))
    except Exception as e:
        print(f'       [!] {name} 读取失败：{type(e).__name__} {e}')
        return None
    return img, pos, sc


def tissue_frac(crop_rgb):
    g = np.asarray(crop_rgb.convert('L'), dtype=np.uint8)
    return float((g <= WHITE_THRESH).mean())


def main():
    torch.set_num_threads(16)
    print('[load] PLIP ...')
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    print(f'       params {sum(p.numel() for p in model.parameters())/1e6:.1f}M')

    texts = [PROMPT_TEMPLATE.format(l) for l in LABELS]
    rng = np.random.RandomState(SEED)

    all_rows = []
    for gsm, name, stage in SLIDES:
        got = load_slide(gsm, name)
        if got is None:
            print(f'[skip] {name}: tar 不在')
            continue
        img, pos, sf = got
        sc = sf['tissue_hires_scalef']
        pos = pos[pos.in_tissue == 1].copy()
        pos['r'] = pos.pxl_row_in_fullres * sc
        pos['c'] = pos.pxl_col_in_fullres * sc
        print(f'\n[{name}] {stage}: hires {img.size}, in_tissue {len(pos)} spots, '
              f'spot dia {sf["spot_diameter_fullres"]*sc:.1f}px')

        # 抽样并裁窗
        take = pos.sample(min(N_SPOT_PER_SLIDE, len(pos)), random_state=SEED)
        crops, keep = [], []
        for _, x in take.iterrows():
            L, T = int(x.c) - CROP // 2, int(x.r) - CROP // 2
            L = max(0, min(L, img.size[0] - CROP)); T = max(0, min(T, img.size[1] - CROP))
            c = img.crop((L, T, L + CROP, T + CROP))
            f = tissue_frac(c)
            if f < MIN_TISSUE_FRAC:
                continue
            crops.append(c); keep.append((x.barcode, f))
        print(f'         裁窗保留 {len(crops)}/{len(take)}（组织占比 ≥{MIN_TISSUE_FRAC}）')

        # 分批打分
        probs = []
        B = 64
        for i in range(0, len(crops), B):
            inp = proc(text=texts, images=crops[i:i+B], return_tensors='pt', padding=True)
            with torch.no_grad():
                probs.append(model(**inp).logits_per_image.softmax(1).numpy())
        probs = np.vstack(probs)

        for (bc, f), p in zip(keep, probs):
            row = {'gsm': gsm, 'sample': name, 'stage': stage, 'barcode': bc,
                   'tissue_frac': f}
            row.update({l: float(v) for l, v in zip(LABELS, p)})
            all_rows.append(row)

        m = probs.mean(0)
        print('         平均概率: ' + '  '.join(f'{l[:22]}={v:.3f}' for l, v in zip(LABELS, m)))

    df = pd.DataFrame(all_rows)
    df.to_csv(f'{OUT}/plip_pilot_spot_scores.csv.gz', index=False)

    # ——— 主判读 ———
    print('\n' + '=' * 72)
    print('主判读：adenocarcinoma 分数（同患者 P4）')
    summ = {}
    for stage in ['Normal', 'AAH', 'LUAD']:
        s = df[df.stage == stage]
        if len(s) == 0:
            continue
        summ[stage] = {
            'n_spot': int(len(s)),
            'adenoca_mean': float(s[TUMOR_LABEL].mean()),
            'adenoca_median': float(s[TUMOR_LABEL].median()),
            'normal_mean': float(s[NORMAL_LABEL].mean()),
            'argmax_adenoca_frac': float((s[LABELS].values.argmax(1) == 0).mean()),
        }
        print(f'  {stage:7s} n={len(s):5d}  adenoca mean={s[TUMOR_LABEL].mean():.4f} '
              f'median={s[TUMOR_LABEL].median():.4f}  '
              f'argmax=adenoca {(s[LABELS].values.argmax(1)==0).mean()*100:5.1f}%  '
              f'normal mean={s[NORMAL_LABEL].mean():.4f}')

    if 'Normal' in summ and 'LUAD' in summ:
        d_mean = summ['LUAD']['adenoca_mean'] - summ['Normal']['adenoca_mean']
        d_argm = summ['LUAD']['argmax_adenoca_frac'] - summ['Normal']['argmax_adenoca_frac']
        print(f'\n  LUAD − Normal:  平均概率差 {d_mean:+.4f}   argmax 占比差 {d_argm*100:+.1f}pp')
        verdict = ('有分离信号（试点通过，可进预注册）' if (d_mean > 0.02 and d_argm > 0.05)
                   else '无分离信号（试点不通过）')
        print(f'  判定：{verdict}')
    else:
        d_mean = d_argm = None
        verdict = '切片不全，无法判定'
        print(f'  {verdict}')

    # ——— 图 ———
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    stages = [s for s in ['Normal', 'AAH', 'LUAD'] if s in set(df.stage)]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), dpi=140)
    ax = axes[0]
    data = [df.loc[df.stage == s, TUMOR_LABEL].values for s in stages]
    bp = ax.boxplot(data, labels=stages, patch_artist=True, showfliers=False)
    for b in bp['boxes']:
        b.set_facecolor('#d62728'); b.set_alpha(.45)
    for i, s in enumerate(stages):
        v = df.loc[df.stage == s, TUMOR_LABEL].values
        ax.scatter(np.random.RandomState(0).normal(i + 1, .05, len(v)), v, s=2,
                   c='k', alpha=.25, rasterized=True)
    ax.set_ylabel(f'PLIP zero-shot P("{TUMOR_LABEL}")')
    ax.set_title(f'PLIP adenocarcinoma score by slide (patient 4)\n'
                 f'LUAD-Normal mean diff = {d_mean:+.4f}' if d_mean is not None else 'PLIP score')

    ax = axes[1]
    mat = np.array([[df.loc[df.stage == s, l].mean() for l in LABELS] for s in stages])
    im = ax.imshow(mat, cmap='magma', aspect='auto')
    ax.set_xticks(range(len(LABELS)))
    ax.set_xticklabels([l.replace('lung tissue with ', '').replace('a histopathology image of ', '')[:26]
                        for l in LABELS], rotation=40, ha='right', fontsize=8)
    ax.set_yticks(range(len(stages))); ax.set_yticklabels(stages)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f'{mat[i,j]:.3f}', ha='center', va='center',
                    color='white' if mat[i, j] < mat.max() * .6 else 'black', fontsize=8)
    ax.set_title('Mean probability per prompt, by slide')
    plt.colorbar(im, ax=ax, fraction=.046)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/plip_pilot_p4_tumor_vs_normal.png', bbox_inches='tight')
    plt.close()
    print(f'\n[out] {FIGDIR}/plip_pilot_p4_tumor_vs_normal.png')

    summary = {
        'script': '07_he_pathology/01_plip_pilot.py',
        'status': 'PILOT — 判定路线可行性，非分析结果',
        'model': 'vinid/plip (CLIP ViT-B/32, 151.3M)',
        'image': 'spatial/tissue_hires_image.png + tissue_hires_scalef',
        'frozen': {'crop': CROP, 'n_spot_per_slide': N_SPOT_PER_SLIDE, 'seed': SEED,
                   'min_tissue_frac': MIN_TISSUE_FRAC, 'white_thresh': WHITE_THRESH,
                   'prompt_template': PROMPT_TEMPLATE, 'labels': LABELS},
        'per_stage': summ,
        'LUAD_minus_Normal_mean_prob': d_mean,
        'LUAD_minus_Normal_argmax_frac': d_argm,
        'verdict': verdict,
        'caveats': [
            'deposited H&E is hires 1951x2000 = 5.66 um/px，PLIP 原生训练尺度约 0.5 um/px ⇒ 差约 11 倍',
            'crop 224px 覆盖约 1.27mm 视野，远大于 spot 直径 55um ⇒ 这是"邻域"分类不是 spot 分类',
            '染色偏淡（全片色度均值 6.6），低对比',
            '同一患者单例对照，非队列证据',
        ],
    }
    with open(f'{OUT}/plip_pilot_summary.json', 'w') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(f'[out] {OUT}/plip_pilot_summary.json')


if __name__ == '__main__':
    main()
