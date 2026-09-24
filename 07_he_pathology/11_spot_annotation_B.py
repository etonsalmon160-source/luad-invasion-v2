#!/usr/bin/env python3
# HE 病理模型 · 第十一步：B 套提示词的**全队列**逐 spot 打分（56 张切片）
#
# 与 04_spot_annotation.py 的唯一差别是提示词：
#   04 = A 套（7 条"分期名"术语 + `a histopathology image of {}`）
#   本脚本 = B 套（原文病理学家逐字判据 + PLIP 官方句式 `An H&E image patch of {}.`）
#   其余冻结常量（CROP/BATCH/MIN_TISSUE_FRAC/WHITE_THRESH）与 04 完全一致，
#   保证与 04 的 spot 集合逐字可比。
#
# ⚠️ 口径提醒（2026-09-24 判决）：B 套在 P4 上比 A 套高（0.891 vs 0.806），
#    但 1,344 条提示词搜索证明这个提升**大部分是密度假象**——一个不算任何提示词的
#    "组织密度方向" d 单独就有 AUC 0.958，去掉 d 后最好的词掉到 0.697。
#    本脚本产出的是**测量值**，不是"更准的病理标签"。判读见 prompt_search_figure.png。
#
# ⚠️ 额外产出：每张切片的**图像向量缓存**（CLIP 图像编码与提示词无关）。
#    有了它，以后换任何提示词只需重跑文本编码器（秒级），不必再过图像编码器。
#
# 可断点续跑：每张切片单独落盘，已存在的跳过。
#
# 跑法：python3 07_he_pathology/11_spot_annotation_B.py

import os, sys, json, re, time, gzip, warnings
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
OUT = f'{ROOT}/results/07_he_pathology/spot_annotation_B'
PER_SLIDE = f'{OUT}/per_slide'
EMB = f'{OUT}/embeds'
os.makedirs(PER_SLIDE, exist_ok=True)
os.makedirs(EMB, exist_ok=True)

# ——— 冻结常量（与 04_spot_annotation.py 逐字一致）———
CROP = 448
BATCH = 512
MIN_TISSUE_FRAC = 0.10
WHITE_THRESH = 230
N_THREADS = 20
PROMPT_TEMPLATE = 'An H&E image patch of {}.'

# ——— B 套：8 条，全部来自原文（Peng 2026 Cancer Cell）Methods 病理学家逐字判据 ———
LABELS = [
    'invasive adenocarcinoma with desmoplastic stroma',                                            # 原文:"invasive adenocarcinoma"
    'morphologically normal pneumocytes lining empty alveoli',                                     # 原文:正常肺判据
    'localized proliferation of thickened alveolar septa lined by atypical type II pneumocytes',   # 原文 AAH 判据
    'neoplastic epithelial cells continuously lining alveolar walls in a lepidic pattern',          # 原文 AIS 判据
    'monotonous enlarged cuboidal type II pneumocytes',                                            # 原文 RPII 判据
    'fibrotic lung tissue with collagen deposition',                                               # 原文排除标准:fibrosis
    'dense lymphocytic infiltrate',                                                                # 原文排除标准:inflammation
    'an empty glass slide',                                                                        # 对照（非病理术语）
]
LESION = [0, 2, 3]    # 腺癌 + AAH + AIS（与 04 的 NEOPLASTIC 同位）
NORMAL = 1
NON_EPI = [5, 6, 7]   # 纤维化 / 淋巴细胞 / 空片（与 04 的 NON_EPI 同位）


def parse_dir(d):
    m = re.match(r'(GSM\d+)_(P\d+)_(Normal|AAH|AIS|MIA|LUAD)', d)
    return (m.group(1), m.group(2), m.group(3)) if m else None


def csv_ok(p):
    """整条 gz 流能读完才算已完成：被杀在写盘中途会留下没有校尾的半截文件，
    只看"文件在不在"会把它当成已完成永久跳过（静默坏数据）。"""
    if not os.path.exists(p):
        return False
    try:
        with gzip.open(p, 'rb') as fh:
            while fh.read(1 << 20):
                pass
        return True
    except Exception:
        return False


def run_slide(model, proc, texts, gdir):
    """逐 spot 打分；返回 (DataFrame, 图像向量) 或 (None, None)。"""
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
        return None, None

    # ——— 图像编码（与提示词无关，缓存下来）———
    embs, probs = [], []
    for i in range(0, len(crops), BATCH):
        inp = proc(text=texts, images=crops[i:i + BATCH], return_tensors='pt', padding=True)
        with torch.no_grad():
            # 图像编码只过一次：model(**inp) 内部会再跑一遍视觉塔。
            # 这里自己复现 logits_per_image，数值等价（实测最大差 ~3e-6），实测快 2.78×。
            ie = model.get_image_features(pixel_values=inp['pixel_values'])
            te = model.get_text_features(input_ids=inp['input_ids'],
                                        attention_mask=inp['attention_mask'])
            ie = ie / ie.norm(dim=-1, keepdim=True)
            te = te / te.norm(dim=-1, keepdim=True)
            embs.append(ie.numpy())
            probs.append((model.logit_scale.exp() * (ie @ te.T)).softmax(1).numpy())
    pr = np.vstack(probs)
    emb = np.vstack(embs)

    df = pd.DataFrame(keep, columns=['barcode', 'array_row', 'array_col',
                                     'pxl_row_in_fullres', 'pxl_col_in_fullres',
                                     'tissue_frac'])
    for j in range(len(LABELS)):
        df[f'p{j}'] = pr[:, j]
    df['lesion_score'] = pr[:, LESION].sum(1)
    df['normal_score'] = pr[:, NORMAL]
    df['nonepi_score'] = pr[:, NON_EPI].sum(1)
    df['argmax_idx'] = pr.argmax(1)
    return df, emb


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
        if os.path.exists(f'{EMB}/{gdir}.npy') and csv_ok(out_f):
            skipped += 1
            continue
        try:
            df, emb = run_slide(model, proc, texts, gdir)
        except Exception as e:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 失败: {type(e).__name__} {e}', flush=True)
            failed += 1
            continue
        if df is None or len(df) == 0:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 无有效 spot', flush=True)
            failed += 1
            continue
        np.save(f'{EMB}/{gdir}.npy', emb)
        df.insert(0, 'gsm', gsm); df.insert(1, 'patient', pat)
        df.insert(2, 'slide_stage', stage); df.insert(3, 'gsm_dir', gdir)
        df.to_csv(out_f, index=False, compression='gzip')
        done += 1; tot_spots += len(df)
        el = time.time() - t_all
        eta = el / (done + skipped) * (len(slides) - done - skipped) / 60 if (done + skipped) else 0
        print(f'  [{i:2d}/{len(slides)}] {gdir:22s} {stage:5s} n={len(df):5d} '
              f'lesion={df.lesion_score.mean():.3f} '
              f'normal={df.normal_score.mean():.3f}  ETA {eta:.0f}min', flush=True)

    print(f'\n[done] 新跑 {done} 张 / 跳过 {skipped} 张 / 失败 {failed} 张；'
          f'本轮新评 {tot_spots:,} spot；总耗时 {(time.time()-t_all)/60:.1f} 分钟')

    summary = {
        'script': '07_he_pathology/11_spot_annotation_B.py',
        'status': 'RAW PROBABILITIES — 只产出测量值，未施加任何类别判定',
        'prompt_set': 'B（原文病理学家逐字判据）',
        'model': 'vinid/plip (CLIP ViT-B/32)',
        'frozen': {'crop': CROP, 'batch': BATCH, 'n_threads': N_THREADS,
                   'min_tissue_frac': MIN_TISSUE_FRAC, 'white_thresh': WHITE_THRESH,
                   'prompt_template': PROMPT_TEMPLATE, 'labels': LABELS,
                   'lesion_idx': LESION, 'normal_idx': NORMAL, 'nonepi_idx': NON_EPI},
        'columns': {f'p{j}': l for j, l in enumerate(LABELS)},
        'compare_with': 'results/07_he_pathology/spot_annotation/（A 套，04 产出）',
        'caveats': [
            'B 套在 P4 上的提升（0.891 vs 0.806）经 1,344 条搜索证明大部分是密度假象：'
            '不算提示词的密度方向 d 单独 AUC=0.958，去掉 d 后最好的词 0.697。见 prompt_search_figure.png',
            'PLIP 给不了 WHO 细分级（AAH/AIS/MIA）；实测 AAH 会抢 LUAD 的 argmax',
            'deposited H&E 5.66um/px，PLIP 原生约 0.5um/px（差约 11 倍）',
            'crop 448px 覆盖约 2.5mm 视野，远大于 spot 直径 55um ⇒ 是邻域分类',
            'Normal 只有 1 张切片（P4_Normal）⇒ 队列级"肿瘤 vs 正常"对照 n=1',
        ],
    }
    with open(f'{OUT}/spot_annotation_B_summary.json', 'w') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(f'[out] {OUT}/spot_annotation_B_summary.json')
    print(f'[out] 图像向量缓存 {EMB}/（{len(os.listdir(EMB))} 张）')


if __name__ == '__main__':
    main()
