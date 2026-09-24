#!/usr/bin/env python3
# HE 病理模型 · 第五步：提示词对照实验（三套提示词在同一批 crop 上比）
#
# 目的：回答"换更专业、来自临床病理报告术语的提示词，能不能救 PLIP"。
#   一次只变一个变量，否则说不清是谁的功劳：
#     A 套 = 原术语 + 原句式            （基线，= 04_spot_annotation.py 已跑的）
#     C 套 = 原术语 + PLIP 官方句式      →  和 A 比，看**句式**的影响
#     B 套 = 报告术语 + PLIP 官方句式    →  和 C 比，看**术语**的影响
#
# B 套术语出处：原文 Methods 里病理学家做 spot 注释时的**逐字判据**
#   （Peng et al. 2026 Cancer Cell，引用 WHO/IASLC 第 45-47 号文献）。
#   不是我编的。逐条对照见 summary 的 provenance 字段。
#
# ⚠️ 关键工程点：**图像编码与提示词无关**。同一批 crop 的三次打分，
#    只需要过一遍图像编码器，再跟三套文本向量各点积一次。省 2/3 算力。
#
# 跑法：python3 07_he_pathology/05_prompt_compare.py

import os, sys, json, time, warnings
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
EXISTING = f'{ROOT}/results/07_he_pathology/spot_annotation/per_slide'
OUT = f'{ROOT}/results/07_he_pathology/prompt_compare'
os.makedirs(OUT, exist_ok=True)

# ——— 冻结常量 ———
CROP = 448                      # 与 04 一致（扫描最优 AUC 0.875）
BATCH = 256
MIN_TISSUE_FRAC = 0.10          # 与 04 一致 → 保留的 spot 集合完全相同
WHITE_THRESH = 230
N_THREADS = 6                   # 主任务占 20 线程，这里限流
N_PER_SLIDE = 4000              # 每片抽样；够算 AUC，且不挤占主任务
SEED = 0

SLIDES = [('GSM9226174_P4_Normal', 'Normal'),
          ('GSM9226175_P4_AAH',    'AAH'),
          ('GSM9226176_P4_AAH-1',  'AAH'),
          ('GSM9226177_P4_LUAD',   'LUAD')]

# ——— 三套提示词 ———
# A：我上一版编的（已被证明是"分期名"不是"病理术语"）—— 基线
A_LABELS = [
    'lung adenocarcinoma',
    'normal lung tissue',
    'lung tissue with atypical adenomatous hyperplasia',
    'lung tissue with adenocarcinoma in situ',
    'fibrous stroma',
    'lymphoid tissue',
    'an empty glass slide',
]

# B：病理报告术语 —— 全部来自原文病理学家的逐字判据
B_LABELS = [
    'invasive adenocarcinoma with desmoplastic stroma',                                  # 原文:"invasive adenocarcinoma"
    'morphologically normal pneumocytes lining empty alveoli',                           # 原文:"morphologically normal pneumocytes lining empty alveoli"
    'localized proliferation of thickened alveolar septa lined by atypical type II pneumocytes',  # 原文 AAH 判据
    'neoplastic epithelial cells continuously lining alveolar walls in a lepidic pattern',        # 原文 AIS 判据
    'monotonous enlarged cuboidal type II pneumocytes',                                  # 原文 RPII 判据
    'fibrotic lung tissue with collagen deposition',                                     # 原文排除标准提到 fibrosis
    'dense lymphocytic infiltrate',                                                      # 原文排除标准提到 inflammation
    'an empty glass slide',                                                              # 对照（非病理术语）
]

SETS = {
    'A_naive_terms+naive_template': {
        'template': 'a histopathology image of {}',
        'labels': A_LABELS,
        'lesion': [0, 2, 3], 'normal': 1, 'empty': 6,
    },
    'C_naive_terms+PLIP_template': {
        'template': 'An H&E image patch of {}.',
        'labels': A_LABELS,
        'lesion': [0, 2, 3], 'normal': 1, 'empty': 6,
    },
    'B_report_terms+PLIP_template': {
        'template': 'An H&E image patch of {}.',
        'labels': B_LABELS,
        'lesion': [0, 2, 3], 'normal': 1, 'empty': 7,
    },
}

PROVENANCE = {
    'B_term_1': '原文:"invasive adenocarcinoma occupying at least 50% of the ST ROI" + "predominant invasive histologic pattern"',
    'B_term_2': '原文 Normal 纳入标准:"morphologically normal pneumocytes lining empty alveoli"',
    'B_term_3': '原文 AAH 判据:"localized proliferations of thickened alveolar septa lined discontinuously by mild to moderate atypical type II pneumocytes or Clara cells"',
    'B_term_4': '原文 AIS 判据:"neoplastic epithelial cells continuously lining mildly thickened alveolar walls in a lepidic pattern"',
    'B_term_5': '原文 RPII 判据:"monotonous, cuboidal-shaped, enlarged type II pneumocytes with a proportional increase of nucleus and cytoplasm"',
    'B_term_6': '原文 LUAD 排除标准提到 "extensive ... fibrosis"',
    'B_term_7': '原文正常片排除标准提到 "extensive inflammation"',
    'B_term_8': '非病理术语，保留作对照（A 套里它系统性抢票 22%）',
}


def auc(y, s):
    """秩法 AUC（免 sklearn 依赖）。"""
    y = np.asarray(y)
    r = pd.Series(np.asarray(s)).rank().values
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float('nan')
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def build_crops(gdir):
    """与 04_spot_annotation.py 完全一致的裁窗逻辑 + 抽样。"""
    sp = f'{DATA}/{gdir}/spatial'
    img = Image.open(f'{sp}/tissue_hires_image.png').convert('RGB')
    pos = pd.read_csv(f'{sp}/tissue_positions.csv')
    sf = json.load(open(f'{sp}/scalefactors_json.json'))
    sc = sf['tissue_hires_scalef']

    pos = pos[pos.in_tissue == 1].copy()
    W, H = img.size
    rows = []
    for r in pos.itertuples():
        cx = r.pxl_col_in_fullres * sc
        cy = r.pxl_row_in_fullres * sc
        L = max(0, min(int(cx) - CROP // 2, W - CROP))
        T = max(0, min(int(cy) - CROP // 2, H - CROP))
        c = img.crop((L, T, L + CROP, T + CROP)).resize((224, 224), Image.BICUBIC)
        f = float((np.asarray(c.convert('L')) <= WHITE_THRESH).mean())
        rows.append((r.barcode, r.array_row, r.array_col, r.pxl_row_in_fullres,
                     r.pxl_col_in_fullres, f, L, T))

    df = pd.DataFrame(rows, columns=['barcode', 'array_row', 'array_col',
                                     'pxl_row_in_fullres', 'pxl_col_in_fullres',
                                     'tissue_frac', 'L', 'T'])
    df = df[df.tissue_frac >= MIN_TISSUE_FRAC].reset_index(drop=True)
    n = len(df)
    if n > N_PER_SLIDE:
        idx = np.sort(np.random.default_rng(SEED).choice(n, N_PER_SLIDE, replace=False))
        df = df.iloc[idx].reset_index(drop=True)
    crops = [img.crop((int(r.L), int(r.T), int(r.L) + CROP, int(r.T) + CROP))
                .resize((224, 224), Image.BICUBIC)
             for r in df.itertuples()]
    return df.drop(columns=['L', 'T']), crops


def main():
    torch.set_num_threads(N_THREADS)
    t0 = time.time()
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    logit_scale = model.logit_scale.exp().item()

    # 文本向量：每套只算一次
    text_emb = {}
    for key, S in SETS.items():
        texts = [S['template'].replace('{}', l) for l in S['labels']]
        inp = proc(text=texts, return_tensors='pt', padding=True)
        with torch.no_grad():
            e = model.get_text_features(**inp)
        text_emb[key] = e / e.norm(dim=-1, keepdim=True)
        print(f'[text] {key}: {len(texts)} 条')
        for t in texts:
            print(f'         {t}')

    frames = []
    for gdir, stage in SLIDES:
        ts = time.time()
        df, crops = build_crops(gdir)
        n = len(crops)
        print(f'\n[slide] {gdir} ({stage}) 抽样 {n} spot', flush=True)

        # 图像向量缓存：编码器与提示词无关，缓存后重做提示词实验就免掉这步
        emb_f = f'{OUT}/embeds_{gdir}.npy'
        if os.path.exists(emb_f):
            ie_all = np.load(emb_f)
            if ie_all.shape[0] == n:
                print(f'        [cache] 复用图像向量 {ie_all.shape}')
            else:
                ie_all = None
        else:
            ie_all = None
        if ie_all is None:
            chunks = []
            for i in range(0, n, BATCH):
                inp = proc(images=crops[i:i + BATCH], return_tensors='pt')
                with torch.no_grad():
                    ie = model.get_image_features(**inp)
                chunks.append((ie / ie.norm(dim=-1, keepdim=True)).numpy())
                del inp, ie
            ie_all = np.vstack(chunks)
            np.save(emb_f, ie_all)
            print(f'        [cache] 写入图像向量 {ie_all.shape}')

        t_ie = torch.from_numpy(ie_all)
        probs_all = {}
        for key in SETS:
            lg = logit_scale * (t_ie @ text_emb[key].T)
            probs_all[key] = lg.softmax(1).numpy()

        for key, S in SETS.items():
            pr = probs_all[key]
            d = df[['barcode', 'array_row', 'array_col', 'pxl_row_in_fullres',
                    'pxl_col_in_fullres', 'tissue_frac']].copy()
            d['set'] = key
            d['stage'] = stage
            d['slide'] = gdir
            d['lesion_score'] = pr[:, S['lesion']].sum(1)
            d['normal_score'] = pr[:, S['normal']]
            d['empty_score'] = pr[:, S['empty']]
            d['argmax'] = [S['labels'][k] for k in pr.argmax(1)]
            d['argmax_idx'] = pr.argmax(1)
            for j, lab in enumerate(S['labels']):
                d[f'p{j}'] = pr[:, j]
            frames.append(d)
        print(f'        {n} spot × {len(SETS)} 套 = {n*len(SETS)} 次判定，'
              f'耗时 {time.time()-ts:.0f}s', flush=True)

    allf = pd.concat(frames, ignore_index=True)
    allf.to_csv(f'{OUT}/p4_prompt_compare_spots.csv.gz', index=False, compression='gzip')

    # ——— 一致性核验：A 套必须与 04 已落盘的结果逐 spot 相同 ———
    consistency = {}
    for gdir, _ in SLIDES:
        f_old = f'{EXISTING}/{gdir}.csv.gz'
        if not os.path.exists(f_old):
            continue
        old = pd.read_csv(f_old)[['barcode', 'neoplastic_score', 'tissue_frac']]
        new = allf[(allf.slide == gdir) & (allf.set == 'A_naive_terms+naive_template')]
        m = new.merge(old, on='barcode', suffixes=('_new', '_old'))
        if len(m):
            consistency[gdir] = {
                'n_matched': int(len(m)),
                'max_abs_diff_lesion_score': float(
                    (m.lesion_score - m.neoplastic_score).abs().max()),
            }

    # ——— 指标 ———
    report = {'consistency_vs_04': consistency, 'sets': {}, 'provenance': PROVENANCE,
              'frozen': {'crop': CROP, 'batch': BATCH, 'min_tissue_frac': MIN_TISSUE_FRAC,
                         'white_thresh': WHITE_THRESH, 'n_threads': N_THREADS,
                         'n_per_slide': N_PER_SLIDE, 'seed': SEED}}
    for key in SETS:
        g = allf[allf.set == key]
        r = {}
        r['auc_LUAD_vs_Normal'] = float(auc((g.stage == 'LUAD').astype(int), g.lesion_score))
        r['auc_AAH_vs_Normal'] = float(auc((g.stage == 'AAH').astype(int), g.lesion_score))
        r['mean_lesion_by_stage'] = {s: float(g[g.stage == s].lesion_score.mean())
                                     for s in ['Normal', 'AAH', 'LUAD']}
        r['mean_normal_by_stage'] = {s: float(g[g.stage == s].normal_score.mean())
                                     for s in ['Normal', 'AAH', 'LUAD']}
        # 密度代理检验：分数跟"空玻璃片"提示词的负相关，就是它在测组织密度的证据
        r['corr_lesion_vs_empty'] = float(np.corrcoef(g.lesion_score, g.empty_score)[0, 1])
        r['corr_lesion_vs_tissuefrac'] = float(np.corrcoef(g.lesion_score, g.tissue_frac)[0, 1])
        if 'empty' in SETS[key]:
            r['corr_lesion_vs_empty_within_slide'] = {
                sl: float(np.corrcoef(gg.lesion_score, gg.empty_score)[0, 1])
                for sl, gg in g.groupby('slide') if len(gg) > 2}
        r['argmax_by_stage'] = {}
        for s in ['Normal', 'AAH', 'LUAD']:
            vc = g[g.stage == s].argmax.value_counts(normalize=True)
            r['argmax_by_stage'][s] = {k: float(v) for k, v in vc.items()}
        report['sets'][key] = r

    with open(f'{OUT}/prompt_compare_report.json', 'w') as fh:
        json.dump(report, fh, ensure_ascii=False, indent=2)

    # 屏上摘要
    print('\n' + '=' * 78)
    print(f'{"套":<34}{"AUC(LUAD/Normal)":>17}{"AUC(AAH/Normal)":>17}')
    print('-' * 78)
    for key in SETS:
        r = report['sets'][key]
        print(f'{key:<34}{r["auc_LUAD_vs_Normal"]:>17.3f}{r["auc_AAH_vs_Normal"]:>17.3f}')
    print('=' * 78)
    for key in SETS:
        r = report['sets'][key]
        print(f'\n【{key}】')
        print('  平均 lesion 分  : ' + '  '.join(
            f'{s}={r["mean_lesion_by_stage"][s]:.3f}' for s in ['Normal', 'AAH', 'LUAD']))
        print(f'  corr(lesion, 空片) = {r["corr_lesion_vs_empty"]:+.3f}   '
              f'corr(lesion, 组织占比) = {r["corr_lesion_vs_tissuefrac"]:+.3f}')
        for s in ['Normal', 'LUAD']:
            top = sorted(r['argmax_by_stage'][s].items(), key=lambda x: -x[1])[:3]
            print(f'  argmax@{s:<7}: ' + '  '.join(f'{k} {v:.0%}' for k, v in top))
    print(f'\n[out] {OUT}/prompt_compare_report.json')
    print(f'[out] {OUT}/p4_prompt_compare_spots.csv.gz')
    print(f'[time] {(time.time()-t0)/60:.1f} 分钟')


if __name__ == '__main__':
    main()
