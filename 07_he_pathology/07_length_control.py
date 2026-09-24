#!/usr/bin/env python3
# 提示词对照 · 补一个 D 套（长度对照），回答"是不是句子变长就够了"
#
# 漏洞：B 套不仅术语更专业，句子也更长。两件事混在一起。
# D 套 = 把 A 套的诊断名**写长**，但加的话仍是"这是什么病"，
#        不是"看得见什么"（不含 empty alveoli / thickened septa / collagen 这类可见结构词）。
#
#   若 D 的密度相关也掉到 -0.2  ⇒ 是**长度**在起作用（那 B 的好消息就贬值了）
#   若 D 仍停在 -0.8 附近        ⇒ 是**术语在描述可见结构**在起作用（B 的好消息成立）
#
# 图像向量已缓存，本脚本只跑文本编码器，秒级完成。

import os, json
import numpy as np
import pandas as pd
import torch
from transformers import CLIPModel, CLIPProcessor

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/07_he_pathology/prompt_compare'
MODEL_DIR = f'{ROOT}/models/plip'
SLIDES = ['GSM9226174_P4_Normal', 'GSM9226175_P4_AAH',
          'GSM9226176_P4_AAH-1', 'GSM9226177_P4_LUAD']

D_LABELS = [   # 与 A 同序同义，只是写长；加粗部分是"名称解释"不是"形态描述"
    'lung adenocarcinoma, a malignant neoplasm of the lung',
    'normal lung tissue, without any pathological change',
    'lung tissue with atypical adenomatous hyperplasia, a precursor lesion',
    'lung tissue with adenocarcinoma in situ, a non-invasive neoplastic lesion',
    'fibrous stroma, a type of connective tissue',
    'lymphoid tissue, a type of immune tissue',
    'an empty glass slide, containing no tissue specimen',
]
TEMPLATE = 'a histopathology image of {}'   # 沿用 A 的句式，只变术语长度

if __name__ == '__main__':
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    scale = model.logit_scale.exp().item()

    texts = [TEMPLATE.replace('{}', l) for l in D_LABELS]
    with torch.no_grad():
        te = model.get_text_features(**proc(text=texts, return_tensors='pt', padding=True))
    te = te / te.norm(dim=-1, keepdim=True)

    # 对齐自检：用同一批缓存向量重算 A 套，必须复现已落盘的 lesion_score
    A_LABELS = ['lung adenocarcinoma', 'normal lung tissue',
                'lung tissue with atypical adenomatous hyperplasia',
                'lung tissue with adenocarcinoma in situ',
                'fibrous stroma', 'lymphoid tissue', 'an empty glass slide']
    with torch.no_grad():
        teA = model.get_text_features(**proc(
            text=[TEMPLATE.replace('{}', l) for l in A_LABELS],
            return_tensors='pt', padding=True))
    teA = teA / teA.norm(dim=-1, keepdim=True)
    f_all = pd.read_csv(f'{D}/p4_prompt_compare_spots.csv.gz')
    maxd = 0.0
    for gdir in SLIDES:
        emb = torch.from_numpy(np.load(f'{D}/embeds_{gdir}.npy'))
        prA = (scale * (emb @ teA.T)).softmax(1).numpy()[:, [0, 2, 3]].sum(1)
        ref = f_all[(f_all.slide == gdir) &
                    (f_all.set == 'A_naive_terms+naive_template')].lesion_score.values
        maxd = max(maxd, float(np.abs(prA - ref).max()))
    print(f'[对齐自检] 用缓存向量重算 A 套，与已落盘的最大绝对差 = {maxd:.2e}'
          f'  {"OK" if maxd < 1e-4 else "*** 不对齐，抽样顺序错 ***"}')
    assert maxd < 1e-4, '缓存向量与 CSV 行序不对齐，后续结论无效'

    frames = []
    for gdir in SLIDES:
        emb = np.load(f'{D}/embeds_{gdir}.npy')
        stage = 'Normal' if 'Normal' in gdir else ('LUAD' if 'LUAD' in gdir else 'AAH')
        pr = (scale * (torch.from_numpy(emb) @ te.T)).softmax(1).numpy()
        sp = pd.read_csv(f'{ROOT}/data/visium_spatial/{gdir}/spatial/tissue_positions.csv')
        # 抽样顺序与 05 完全一致（seed 0 + 同过滤），才能与已存的 tissue_frac 对齐
        f = pd.read_csv(f'{D}/p4_prompt_compare_spots.csv.gz')
        sub = f[(f.slide == gdir) & (f.set == 'A_naive_terms+naive_template')]
        tf = sub.tissue_frac.values
        assert len(tf) == pr.shape[0], f'{gdir}: {len(tf)} vs {pr.shape[0]}'
        frames.append(pd.DataFrame({
            'slide': gdir, 'stage': stage, 'tissue_frac': tf,
            'lesion_score': pr[:, [0, 2, 3]].sum(1),
            'normal_score': pr[:, 1], 'empty_score': pr[:, 6],
            'adenoca_score': pr[:, 0], 'argmax_idx': pr.argmax(1)}))
    g = pd.concat(frames, ignore_index=True)

    def auc(y, s):
        y = np.asarray(y)
        r = pd.Series(np.asarray(s)).rank().values
        n1 = int(y.sum()); n0 = len(y) - n1
        return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

    isL = (g.stage == 'LUAD').astype(int).values
    res = {
        'set': 'D_naive_terms_lengthened+naive_template',
        'labels': D_LABELS,
        'auc_LUAD_vs_Normal_aggregate': float(auc(isL, g.lesion_score)),
        'auc_LUAD_vs_Normal_adenoca_probe': float(auc(isL, g.adenoca_score)),
        'corr_lesion_vs_empty': float(np.corrcoef(g.lesion_score, g.empty_score)[0, 1]),
        'corr_lesion_vs_tissuefrac': float(np.corrcoef(g.lesion_score, g.tissue_frac)[0, 1]),
        'corr_lesion_vs_empty_within_slide': {
            s: float(np.corrcoef(gg.lesion_score, gg.empty_score)[0, 1])
            for s, gg in g.groupby('slide')},
        'mean_lesion_by_stage': {s: float(g[g.stage == s].lesion_score.mean())
                                 for s in ['Normal', 'AAH', 'LUAD']},
    }
    with open(f'{D}/length_control.json', 'w') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=2)

    print('=== D 套（把诊断名写长，但不描述可见结构）===')
    for t in texts:
        print('   ', t)
    print(f"\n  AUC(LUAD/Normal) 合成={res['auc_LUAD_vs_Normal_aggregate']:.3f}  "
          f"单用腺癌条={res['auc_LUAD_vs_Normal_adenoca_probe']:.3f}")
    print(f"  corr(lesion, 空片) = {res['corr_lesion_vs_empty']:+.3f}   "
          f"corr(lesion, 组织占比) = {res['corr_lesion_vs_tissuefrac']:+.3f}")
    print('  逐张：' + '  '.join(f'{k.split("_P4_")[1]}={v:+.3f}'
                                 for k, v in res['corr_lesion_vs_empty_within_slide'].items()))
    print('  平均 lesion 分：' + '  '.join(
        f'{s}={res["mean_lesion_by_stage"][s]:.3f}' for s in ['Normal', 'AAH', 'LUAD']))
