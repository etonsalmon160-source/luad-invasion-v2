#!/usr/bin/env python3
# 提示词搜索：在缓存图像向量上穷举「句式 × 肿瘤术语 × 正常术语」的二元组合
#
# 为什么免费：图像编码器与提示词无关，04->05 已把 P4 的图像向量缓存到 .npy。
#             这里只跑文本编码器（几百条，秒级）+ 若干次点积。
#
# 为什么搜「二元对」而不是「多类」：PLIP 自己的肺 benchmark（WSSS4LUAD）就是 tumor/normal
#             二元。多类的 softmax 会让每条提示词的分数相互纠缠，搜不出"哪条词好"。
#
# 两个轴（不合成单一目标，交给用户选）：
#     y = AUC(LUAD spot vs Normal spot)          —— 判别力
#     x = |corr(score, tissue_frac)|             —— 密度依赖（越小越"真在看组织"）
#
# ⚠️ 防过拟合：P4 的 spot 按空间棋盘格切成 训练/留出 两半，只在训练半上挑，
#     留出半上的数字才是能看的。且 P4 只有 1 个患者，跨患者验证另做（见 09）。
#
# 跑法：python3 07_he_pathology/08_prompt_search.py

import os, json, itertools
import numpy as np
import pandas as pd
import torch
from transformers import CLIPModel, CLIPProcessor

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/07_he_pathology/prompt_compare'
MODEL_DIR = f'{ROOT}/models/plip'
SLIDES = ['GSM9226174_P4_Normal', 'GSM9226175_P4_AAH',
          'GSM9226176_P4_AAH-1', 'GSM9226177_P4_LUAD']

# ——— 句式（前 4 条来自 PLIP 官方代码里的真实写法）———
TEMPLATES = [
    ('An H&E image patch of ', '.'),
    ('a histopathology image of ', ''),
    ('An H&E image patch of ', ' tissue.'),
    ('a histopathology image of ', ' tissue'),
    ('a micrograph of ', ''),
    ('a photomicrograph of ', ''),
    ('a histology image showing ', ''),
    ('H&E stained ', ''),
]

# ——— 肿瘤侧术语 ———
# 前 3 条来自原文病理学家判据；其余为病理报告的惯用形态学描述（**描述看得见的东西**）
TUMOR = [
    'invasive adenocarcinoma with desmoplastic stroma',          # 原文判据
    'invasive adenocarcinoma',                                    # 原文判据
    'neoplastic epithelial cells continuously lining alveolar walls in a lepidic pattern',  # 原文 AIS/LUAD 判据
    'solid sheets of tumor cells replacing alveoli',
    'adenocarcinoma cells forming glands',
    'tumor cells replacing normal alveolar architecture',
    'atypical epithelial cells lining distorted alveoli',
    'malignant epithelium with enlarged hyperchromatic nuclei',
    'solid tumor growth pattern',
    'gland-forming malignant epithelium',
    'carcinoma cells',
    'malignant tumor tissue',
    'dense cellular tumor with loss of alveolar spaces',
    'invasive malignant epithelium',
]

# ——— 正常侧术语 ———
NORMAL = [
    'morphologically normal pneumocytes lining empty alveoli',    # 原文判据
    'normal lung parenchyma with patent alveolar spaces',
    'thin alveolar walls lined by normal pneumocytes',
    'aerated lung tissue with open alveolar spaces',
    'empty alveoli lined by flat pneumocytes',
    'normal lung with thin alveolar septa',
    'intact alveolar spaces',
    'unremarkable lung parenchyma',
    'normal lung tissue',
    'alveolar spaces lined by normal epithelium',
    'open alveolar spaces with thin septa',
    'healthy lung tissue',
]


def auc(y, s):
    y = np.asarray(y)
    r = pd.Series(np.asarray(s)).rank().values
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float('nan')
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    scale = model.logit_scale.exp().item()

    # ——— 载入缓存图像向量 + 元数据 ———
    f = pd.read_csv(f'{D}/p4_prompt_compare_spots.csv.gz')
    meta, emb = None, []
    for gdir in SLIDES:
        e = np.load(f'{D}/embeds_{gdir}.npy')
        sub = f[(f.slide == gdir) & (f.set == 'A_naive_terms+naive_template')]
        assert len(sub) == len(e)
        m = sub[['array_row', 'array_col', 'tissue_frac']].copy()
        m['stage'] = sub.stage.values
        meta = m if meta is None else pd.concat([meta, m], ignore_index=True)
        emb.append(e)
    E = np.vstack(emb)                       # (n_img, 512) 已 L2 归一
    meta = meta.reset_index(drop=True)
    n_img = len(E)
    print(f'图像向量 {E.shape}   Normal={sum(meta.stage=="Normal")} '
          f'AAH={sum(meta.stage=="AAH")}  LUAD={sum(meta.stage=="LUAD")}')

    # ——— 空间棋盘格切训练/留出（随机切会因空间自相关泄漏）———
    blk = ((meta.array_row // 8) + (meta.array_col // 8)) % 2
    tr = (blk == 0).values
    te = (blk == 1).values
    isL = (meta.stage == 'LUAD').values

    # ——— 把所有候选文本编码一次 ———
    texts, keys = [], []
    for ti, (pre, suf) in enumerate(TEMPLATES):
        for t in TUMOR:
            texts.append(pre + t + suf); keys.append(('T', ti, t))
        for t in NORMAL:
            texts.append(pre + t + suf); keys.append(('N', ti, t))
    with torch.no_grad():
        TE = model.get_text_features(**proc(text=texts, return_tensors='pt', padding=True))
    TE = (TE / TE.norm(dim=-1, keepdim=True)).numpy()
    print(f'文本 {TE.shape[0]} 条（{len(TEMPLATES)} 句式 × ({len(TUMOR)}+{len(NORMAL)}) 术语）')

    # 相似度矩阵 (n_text, n_img)，一次算完
    S = TE @ E.T
    IDX = {k: i for i, k in enumerate(keys)}

    # ——— 穷举二元对 ———
    rows = []
    for ti, (pre, suf) in enumerate(TEMPLATES):
        for t in TUMOR:
            i_t = IDX[('T', ti, t)]
            for nrm in NORMAL:
                i_n = IDX[('N', ti, nrm)]
                lg = scale * np.vstack([S[i_t], S[i_n]])       # (2, n_img)
                ex = np.exp(lg - lg.max(0, keepdims=True))
                p_t = ex[0] / ex.sum(0)
                rows.append({
                    'template': f'{pre}[] {suf}'.strip(),
                    'tumor_term': t, 'normal_term': nrm,
                    'auc_train': auc(isL[tr], p_t[tr]),
                    'auc_test': auc(isL[te], p_t[te]),
                    'dens_test': abs(np.corrcoef(p_t[te], meta.tissue_frac.values[te])[0, 1]),
                    'dens_all': abs(np.corrcoef(p_t, meta.tissue_frac.values)[0, 1]),
                })
    R = pd.DataFrame(rows)
    R.to_csv(f'{D}/prompt_search_all_pairs.csv', index=False)
    print(f'共评 {len(R)} 个二元对 → prompt_search_all_pairs.csv')

    # ——— Pareto 前沿（留出集上：密度依赖越小越好、AUC 越大越好）———
    R2 = R.sort_values('dens_test')
    front, best = [], -1
    for _, r in R2.iterrows():
        if r.auc_test > best:
            front.append(r); best = r.auc_test
    F = pd.DataFrame(front)

    print('\n=== 留出集上的 Pareto 前沿（按密度依赖从小到大）===')
    print(f"{'AUC(test)':>9}{'|dens|':>8}  {'AUC(tr)':>8}  肿瘤术语 / 正常术语")
    print('-' * 118)
    for _, r in F.iterrows():
        print(f"{r.auc_test:>9.3f}{r.dens_test:>8.3f}  {r.auc_train:>8.3f}  "
              f"{r.tumor_term[:46]:<46} | {r.normal_term[:40]}")
    F.to_csv(f'{D}/prompt_search_pareto.csv', index=False)

    # 与现有基线对照
    print('\n=== 对照：现有各套（同一留出集）===')
    for name, k in [('A 原版 lesion', 'A_naive_terms+naive_template'),
                    ('B 报告术语 lesion', 'B_report_terms+PLIP_template')]:
        g = f[(f.set == k)].reset_index(drop=True)
        assert len(g) == n_img and np.allclose(g.tissue_frac.values, meta.tissue_frac.values)
        s = g.lesion_score.values
        print(f'  {name:<20} AUC(test)={auc(isL[te], s[te]):.3f}  '
              f'|dens|={abs(np.corrcoef(s[te], meta.tissue_frac.values[te])[0,1]):.3f}')

    # 单用最强肿瘤条（B 套）
    g = f[(f.set == 'B_report_terms+PLIP_template')].reset_index(drop=True)
    s = g['p0'].values
    print(f'  {"B 单用腺癌条":<20} AUC(test)={auc(isL[te], s[te]):.3f}  '
          f'|dens|={abs(np.corrcoef(s[te], meta.tissue_frac.values[te])[0,1]):.3f}')

    json.dump({'n_pairs': len(R), 'n_img': int(n_img),
               'templates': [f'{p}[] {s}'.strip() for p, s in TEMPLATES],
               'tumor_terms': TUMOR, 'normal_terms': NORMAL},
              open(f'{D}/prompt_search_meta.json', 'w'), ensure_ascii=False, indent=2)


if __name__ == '__main__':
    main()
