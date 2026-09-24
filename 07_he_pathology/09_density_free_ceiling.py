#!/usr/bin/env python3
# 诊断：把"密度方向"从文本向量里正交投影掉，看还剩多少判别力
#
# 动机：搜索出的高 AUC 提示词密度相关也高（0.97 AUC @ 0.57 dens），
#       怀疑它仍在读"这块白不白"。问：**去掉密度方向之后，还剩多少真信号？**
#
# 做法：在图像向量空间里定义"密度方向" d = mean(高组织占比) − mean(低组织占比)，
#       对文本向量 t 做 t' = t − (t·d)d 再归一化 ⇒ 该提示词**在数学上无法使用密度信息**。
#       看 AUC 掉多少。
#
# ⚠️ 重要解读限制：**肿瘤组织本来就比正常肺致密**，密度里有一部分是**真信号**。
#     投影掉 d 会把这部分真信号一起削掉。所以这个数不是"应该达到的目标"，
#     而是"完全不靠密度时还剩多少"的**下界**。
#
# 跑法：python3 07_he_pathology/09_density_free_ceiling.py

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

# 候选（从 08 的 Pareto 前沿上取点，覆盖整个权衡区间）
CAND = [
    ('max-AUC',       'a histopathology image of ', 'invasive adenocarcinoma with desmoplastic stroma', ''),
    ('balanced',      'a histopathology image of ', 'carcinoma cells', ''),
    ('low-density',   'a histology image showing ', 'carcinoma cells', ''),
    ('nuclei-term',   'a photomicrograph of ', 'malignant epithelium with enlarged hyperchromatic nuclei', ''),
]
NORMAL_FIXED = ('a histopathology image of ', 'morphologically normal pneumocytes lining empty alveoli', '')


def auc(y, s):
    y = np.asarray(y)
    r = pd.Series(np.asarray(s)).rank().values
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float('nan')
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


if __name__ == '__main__':
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    scale = model.logit_scale.exp().item()

    f = pd.read_csv(f'{D}/p4_prompt_compare_spots.csv.gz')
    meta, emb = None, []
    for gdir in SLIDES:
        e = np.load(f'{D}/embeds_{gdir}.npy')
        sub = f[(f.slide == gdir) & (f.set == 'A_naive_terms+naive_template')]
        m = sub[['tissue_frac', 'stage']].copy()
        meta = m if meta is None else pd.concat([meta, m], ignore_index=True)
        emb.append(e)
    E = np.vstack(emb)
    meta = meta.reset_index(drop=True)
    tf = meta.tissue_frac.values
    isL = (meta.stage == 'LUAD').values

    # ——— 密度方向：组织最多的 25% 与最少的 25% 的图像向量均值之差 ———
    lo, hi = np.quantile(tf, .25), np.quantile(tf, .75)
    d = E[tf >= hi].mean(0) - E[tf <= lo].mean(0)
    d = d / np.linalg.norm(d)
    print(f'密度方向已估计（高组织占比 {sum(tf>=hi)} 张 vs 低 {sum(tf<=lo)} 张）')
    # 自检：d 与 tissue_frac 的相关应很强
    proj_img = E @ d
    print(f'  自检 corr(<img,d>, tissue_frac) = {np.corrcoef(proj_img, tf)[0,1]:+.3f}')
    print(f'  自检 corr(<img,d>, 是LUAD)      = {np.corrcoef(proj_img, isL.astype(int))[0,1]:+.3f}')

    def encode(pre, term, suf):
        with torch.no_grad():
            e = model.get_text_features(**proc(text=[pre + term + suf],
                                               return_tensors='pt', padding=True))
        e = (e / e.norm(dim=-1, keepdim=True))[0].numpy()
        return e

    def score(e_t, e_n, project):
        a, b = e_t.copy(), e_n.copy()
        if project:
            a = a - (a @ d) * d; b = b - (b @ d) * d
            a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
        lg = scale * np.vstack([E @ a, E @ b])
        ex = np.exp(lg - lg.max(0, keepdims=True))
        return ex[0] / ex.sum(0)

    e_n = encode(*NORMAL_FIXED)
    rows = []
    for name, pre, term, suf in CAND:
        e_t = encode(pre, term, suf)
        r = {'name': name, 'tumor_term': term}
        for lbl, pj in [('raw', False), ('density_projected', True)]:
            s = score(e_t, e_n, pj)
            r[f'auc_{lbl}'] = auc(isL, s)
            r[f'dens_{lbl}'] = abs(np.corrcoef(s, tf)[0, 1])
        rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(f'{D}/density_free_ceiling.csv', index=False)

    print(f'\n正常侧固定为: {NORMAL_FIXED[1]}')
    print(f"\n{'候选':<14}{'原AUC':>8}{'原|dens|':>10}{'投影后AUC':>11}{'投影后|dens|':>13}")
    print('-' * 60)
    for _, r in R.iterrows():
        print(f"{r['name']:<14}{r.auc_raw:>8.3f}{r.dens_raw:>10.3f}"
              f"{r.auc_density_projected:>11.3f}{r.dens_density_projected:>13.3f}")
    print('\n（投影后 AUC 是"完全不靠密度时还剩多少"的下界，不是目标值）')
