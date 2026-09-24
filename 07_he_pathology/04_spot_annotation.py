#!/usr/bin/env python3
# HE 病理模型 · 第四步：全队列逐 spot 病理注释（56 张切片 / 639,816 spot）
#
# 目的：产出**逐 spot 的病理标签**，供下游用：
#   ① 分出上皮 spot（原文 SpatialInferCNV 只在上皮 spot 上跑）
#   ② 挑出"正常肺上皮" spot → 空转 CNV 的参考锚
#   ③ 与 RCTD 组成第三票，判别恶性
#
# ⚠️ 本脚本只产出**原始概率**（测量），不做阈值判定。
#    类别归属规则（哪一档算"病变"）属口径，需用户签字后才可施加 —— 见 summary 的 pending_rule。
#
# ⚠️ PLIP 给不了 WHO 细分级（AAH/AIS/MIA）。实测 AAH 会抢 LUAD 的 argmax。
#    细分级在原文由病理学家给出，我们这里**不冒充**；切片级分期来自目录名（GEO 标题）。
#
# 可断点续跑：每张切片单独落盘，已存在的跳过。
#
# 跑法：python3 07_he_pathology/04_spot_annotation.py

import os, sys, json, re, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import torch
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
from transformers import CLIPModel, CLIPProcessor

ROOT = '/home/eto/luad_v2'
MODEL_DIR = f'{ROOT}/models/plip'
DATA = f'{ROOT}/data/visium_spatial'
OUT = f'{ROOT}/results/07_he_pathology/spot_annotation'
PER_SLIDE = f'{OUT}/per_slide'
os.makedirs(PER_SLIDE, exist_ok=True)

# ——— 冻结常量（与 03_cohort_screen.py 一致，保证跨步骤可比）———
CROP = 448                     # 扫描最优（AUC 0.875）；喂 PLIP 前缩到 224
BATCH = 512                    # 基准实测：512→29.9 spot/s（全量约 5.9h）；64→14.5（12.3h）
MIN_TISSUE_FRAC = 0.10
WHITE_THRESH = 230
N_THREADS = 20
PROMPT_TEMPLATE = 'a histopathology image of {}'
# 提示词**沿用已验证的 7 个**：增删会重新归一化 softmax，使 03 的 rho=0.70 失效
LABELS = [
    'lung adenocarcinoma',
    'normal lung tissue',
    'lung tissue with atypical adenomatous hyperplasia',
    'lung tissue with adenocarcinoma in situ',
    'fibrous stroma',
    'lymphoid tissue',
    'an empty glass slide',
]
NEOPLASTIC = [0, 2, 3]         # adenoca + AAH + AIS 合成 neoplastic 分数
# 中间档：既非上皮病变、也非正常上皮（下游分 epi/non-epi 时用）
NON_EPI = [4, 5, 6]            # fibrous stroma / lymphoid / empty


def parse_dir(d):
    m = re.match(r'(GSM\d+)_(P\d+)_(Normal|AAH|AIS|MIA|LUAD)', d)
    return (m.group(1), m.group(2), m.group(3)) if m else None


def run_slide(model, proc, texts, gdir, stage):
    """逐 spot 打分；返回 DataFrame 或 None。"""
    sp = f'{DATA}/{gdir}/spatial'
    img = Image.open(f'{sp}/tissue_hires_image.png').convert('RGB')
    pos = pd.read_csv(f'{sp}/tissue_positions.csv')
    sf = json.load(open(f'{sp}/scalefactors_json.json'))
    sc = sf['tissue_hires_scalef']

    pos = pos[pos.in_tissue == 1].copy()
    pos['cx'] = pos.pxl_col_in_fullres * sc
    pos['cy'] = pos.pxl_row_in_fullres * sc

    W, H = img.size
    crops, keep = [], []
    for r in pos.itertuples():
        L = max(0, min(int(r.cx) - CROP // 2, W - CROP))
        T = max(0, min(int(r.cy) - CROP // 2, H - CROP))
        c = img.crop((L, T, L + CROP, T + CROP)).resize((224, 224), Image.BICUBIC)
        f = float((np.asarray(c.convert('L')) <= WHITE_THRESH).mean())
        if f < MIN_TISSUE_FRAC:
            continue
        crops.append(c)
        keep.append((r.barcode, r.array_row, r.array_col, r.pxl_row_in_fullres,
                     r.pxl_col_in_fullres, f))
    if not crops:
        return None

    probs = []
    for i in range(0, len(crops), BATCH):
        inp = proc(text=texts, images=crops[i:i + BATCH], return_tensors='pt', padding=True)
        with torch.no_grad():
            probs.append(model(**inp).logits_per_image.softmax(1).numpy())
    pr = np.vstack(probs)

    df = pd.DataFrame(keep, columns=['barcode', 'array_row', 'array_col',
                                     'pxl_row_in_fullres', 'pxl_col_in_fullres',
                                     'tissue_frac'])
    for j, l in enumerate(LABELS):
        df['p_' + l.replace(' ', '_')] = pr[:, j]
    df['neoplastic_score'] = pr[:, NEOPLASTIC].sum(1)
    df['nonepi_score'] = pr[:, NON_EPI].sum(1)
    df['normal_score'] = pr[:, 1]
    df['argmax_label'] = [LABELS[k] for k in pr.argmax(1)]
    df['argmax_idx'] = pr.argmax(1)
    return df


def main():
    torch.set_num_threads(N_THREADS)
    t_all = time.time()
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    texts = [PROMPT_TEMPLATE.format(l) for l in LABELS]

    slides = []
    for d in sorted(os.listdir(DATA)):
        p = parse_dir(d)
        if p and os.path.isdir(f'{DATA}/{d}'):
            slides.append((*p, d))

    done = skipped = failed = 0
    tot_spots = 0
    for i, (gsm, pat, stage, gdir) in enumerate(slides, 1):
        out_f = f'{PER_SLIDE}/{gdir}.csv.gz'
        if os.path.exists(out_f):
            skipped += 1
            continue
        try:
            df = run_slide(model, proc, texts, gdir, stage)
        except Exception as e:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 失败: {type(e).__name__} {e}', flush=True)
            failed += 1
            continue
        if df is None or len(df) == 0:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 无有效 spot', flush=True)
            failed += 1
            continue
        df.insert(0, 'gsm', gsm); df.insert(1, 'patient', pat)
        df.insert(2, 'slide_stage', stage); df.insert(3, 'gsm_dir', gdir)
        df.to_csv(out_f, index=False, compression='gzip')
        done += 1; tot_spots += len(df)
        el = time.time() - t_all
        eta = el / (done + skipped) * (len(slides) - done - skipped) / 60 if (done + skipped) else 0
        print(f'  [{i:2d}/{len(slides)}] {gdir:22s} {stage:5s} n={len(df):5d} '
              f'neoplastic={df.neoplastic_score.mean():.3f} '
              f'normal={df.normal_score.mean():.3f}  ETA {eta:.0f}min', flush=True)

    print(f'\n[done] 新跑 {done} 张 / 跳过 {skipped} 张 / 失败 {failed} 张；'
          f'本轮新评 {tot_spots:,} spot；总耗时 {(time.time()-t_all)/60:.1f} 分钟')

    summary = {
        'script': '07_he_pathology/04_spot_annotation.py',
        'status': 'RAW PROBABILITIES — 只产出测量值，未施加任何类别判定',
        'model': 'vinid/plip (CLIP ViT-B/32)',
        'frozen': {'crop': CROP, 'batch': BATCH, 'n_threads': N_THREADS,
                   'min_tissue_frac': MIN_TISSUE_FRAC, 'white_thresh': WHITE_THRESH,
                   'prompt_template': PROMPT_TEMPLATE, 'labels': LABELS,
                   'neoplastic_idx': NEOPLASTIC, 'nonepi_idx': NON_EPI},
        'why_prompt_set_unchanged':
            '沿用 03_cohort_screen.py 的 7 提示词；增删会重新归一化 softmax，'
            '使那边已验证的 Spearman rho=0.695 (p=2.8e-09) 失效',
        'n_slides_run': done, 'n_slides_skipped': skipped, 'n_slides_failed': failed,
        'n_spot_this_run': tot_spots,
        'pending_rule': {
            'note': '类别归属规则属口径，未经签字不得施加。本脚本只落盘原始概率。',
            'argmax_ruled_out': {
                'evidence': 'P4 逐 spot 实测（n=791）：腺癌在 Normal 与 LUAD 上被选为 argmax 的比例均为 0.0%；'
                            'LUAD 的 argmax 有 66.2% 落在「AAH」、22.8% 落在「空玻璃片」；'
                            '「空玻璃片」在 Normal 真组织上也占 21.5%',
                'conclusion': 'argmax 不可用于打标签 —— 腺癌永远不是 argmax，且空片提示词系统性抢票',
            },
            'candidate': {
                'epithelium_normal': 'normal_score 高（且非背景）',
                'epithelium_lesion': 'neoplastic_score 高',
                'non_epithelial': 'nonepi_score（stroma+lymphoid+empty）高',
                'background': '由 tissue_frac 与 empty 提示词共同判定',
            },
            'open_questions': [
                '三档还是五档？两两之间的分界用哪个分数、阈值多少？',
                '阈值必须用 P4_Normal 正常区与已知病变区定标后再冻结，不得回头看结果凑（法则 3.2）',
                '低组织占比 spot（tissue_frac 接近 0.10）是否单独一档？',
            ],
        },
        'caveats': [
            'PLIP 给不了 WHO 细分级（AAH/AIS/MIA）；实测 AAH 会抢 LUAD 的 argmax',
            '细分级在原文由病理学家逐 spot 给出；我们无此数据，不得冒充',
            'deposited H&E 5.66um/px，PLIP 原生约 0.5um/px（差约 11 倍）',
            'crop 448px 覆盖约 2.5mm 视野，远大于 spot 直径 55um ⇒ 是邻域分类',
            '零样本模型未在肺 FFPE 上微调；本注释是 R 侧 CNV 的锚定输入，不是金标准',
        ],
    }
    with open(f'{OUT}/spot_annotation_summary.json', 'w') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(f'[out] {OUT}/spot_annotation_summary.json')


if __name__ == '__main__':
    main()
