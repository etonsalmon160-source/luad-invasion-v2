#!/usr/bin/env python3
# SCMG S3 方向审计：5 个患者 AUROC < 0.5 到底是「符号写反」还是「患者批次偏移压过类别分离」。
#
# 不产生任何新的科学结论，只做法医式核对。输入是已存在的产物，不改任何已签口径。
#
# 判据（先写死）：
#   · 符号假设：若 S3 的符号写反，则 s = d0 - d1 应改成 s = d1 - d0 后 23/23 患者都 >0.5，
#     且 pooled 也对称地变成 1 - 0.571 = 0.429。**若只是部分患者 <0.5，符号假设不成立。**
#   · 批次假设：把每个患者的潜向量减去**该患者自己的均值**（只去患者整体偏移，
#     两类共享同一个平移 ⇒ 类间差保留），再重跑同一套留一患者最近质心。
#     若翻转的患者回到 >0.5 ⇒ 是患者偏移造成的，不是方向错。

import json
import os

import numpy as np
import pandas as pd

ROOT = '/home/eto/luad_v2'
D = f'{ROOT}/results/06_scmg'
OUT = f'{D}/audit_s3'
os.makedirs(OUT, exist_ok=True)

NEG, POS = 'Normal', 'IAC'
AUROC_STRONG, AUROC_WEAK = 0.80, 0.60


def auroc(pos, neg):
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    s = np.r_[pos, neg]
    r = pd.Series(s).rank(method='average').values
    n1, n0 = y.sum(), (1 - y).sum()
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def loo_centroid_scores(X, pid, y):
    """复刻 02_analyze_scmg_branch.py 的 S3：留一患者，质心来自其他患者。"""
    tot1, tot0 = X[y == 1].sum(0), X[y == 0].sum(0)
    s_all, y_all, p_all = [], [], []
    rows = []
    for p in np.unique(pid):
        te = pid == p
        tr = ~te
        n1_tr, n0_tr = int((y[tr] == 1).sum()), int((y[tr] == 0).sum())
        if n1_tr == 0 or n0_tr == 0:
            continue
        c1 = (tot1 - X[te & (y == 1)].sum(0)) / n1_tr
        c0 = (tot0 - X[te & (y == 0)].sum(0)) / n0_tr
        d1 = np.linalg.norm(X[te] - c1, axis=1)
        d0 = np.linalg.norm(X[te] - c0, axis=1)
        s = d0 - d1
        s_all.append(s); y_all.append(y[te]); p_all.append(np.repeat(p, int(te.sum())))
        n_pos, n_neg = int((y[te] == 1).sum()), int((y[te] == 0).sum())
        a = auroc(s[y[te] == 1], s[y[te] == 0]) if n_pos and n_neg else float('nan')
        rows.append(dict(patient=p, n_iac=n_pos, n_normal=n_neg, auroc=a,
                         sep_global=float(np.linalg.norm(c1 - c0)),
                         d_ownI_to_cI=float(np.linalg.norm(X[te & (y == 1)].mean(0) - c1)),
                         d_ownI_to_cN=float(np.linalg.norm(X[te & (y == 1)].mean(0) - c0)),
                         d_ownN_to_cN=float(np.linalg.norm(X[te & (y == 0)].mean(0) - c0)),
                         d_ownN_to_cI=float(np.linalg.norm(X[te & (y == 0)].mean(0) - c1))))
    return np.concatenate(s_all), np.concatenate(y_all), np.concatenate(p_all), pd.DataFrame(rows)


def main():
    latent = np.load(f'{D}/scmg_latent.npy')
    proj = pd.read_csv(f'{D}/scmg_projection.csv.gz')
    print(f'[in] 潜空间 {latent.shape}；投影 {len(proj)} 行')

    mask = proj['stage'].isin([NEG, POS]).values
    pid_all = proj.loc[mask, 'patient_id'].values
    y_all = (proj.loc[mask, 'stage'].values == POS).astype(int)

    rep = {'note': '法医核对：S3 逐患者 AUROC<0.5 的成因。不产生新科学结论。'}

    # ---- 0) 每患者的整体偏移量（该患者全部细胞均值 vs 其余患者均值）----
    X0 = latent[mask].astype(np.float64)
    gmean = X0.mean(0)
    shift = {}
    for p in np.unique(pid_all):
        m = pid_all == p
        shift[str(p)] = float(np.linalg.norm(X0[m].mean(0) - gmean))
    rep['patient_shift_from_global_mean'] = {k: round(v, 3) for k, v in
                                             sorted(shift.items(), key=lambda kv: -kv[1])}

    # ---- 1) 原样复刻 —— 必须复现出 0.436 / 0.461 / ... 那几个值 ----
    s, y, pid, per = loo_centroid_scores(X0, pid_all, y_all)
    a_pool = auroc(s[y == 1], s[y == 0])
    print(f'\n[1] 复刻原 S3：pooled AUROC = {a_pool:.4f}（原记录 0.5709753246736938）')
    print(per[['patient', 'n_iac', 'n_normal', 'auroc', 'sep_global']].to_string(index=False))
    rep['reproduced_pooled_auroc'] = round(a_pool, 6)
    rep['per_patient_original'] = per.replace({np.nan: None}).to_dict('records')
    n_below = int((per['auroc'] < 0.5).sum())
    rep['n_patients_below_0.5'] = n_below
    print(f'    <0.5 的患者数 = {n_below} / {len(per)}')

    # ---- 2) 符号假设检验：整体取反 ----
    a_flip = auroc(-s[y == 1], -s[y == 0])
    print(f'\n[2] 符号假设：把 s 整体取反后 pooled AUROC = {a_flip:.4f}'
          f'（= 1 - {a_pool:.4f} = {1-a_pool:.4f}；恒等式，必然成立）')
    rep['sign_flip_pooled_auroc'] = round(a_flip, 6)
    rep['sign_flip_note'] = ('整体取反是把 0.571 变成 0.429，只会让更多的患者掉到 0.5 以下，'
                             '所以「符号写反」这个假设与观测不符：观测是**一部分**患者 <0.5、'
                             '另一部分高达 0.80。符号错会全体一致地反向。')

    # ---- 3) 批次假设：逐患者去均值后重跑 ----
    Xc = np.empty_like(X0)
    for p in np.unique(pid_all):
        m = pid_all == p
        Xc[m] = X0[m] - X0[m].mean(0)
    del X0
    sc, yc, pc, perc = loo_centroid_scores(Xc, pid_all, y_all)
    a_c = auroc(sc[yc == 1], sc[yc == 0])
    print(f'\n[3] 批次假设：逐患者去均值后 pooled AUROC = {a_c:.4f}')
    print(perc[['patient', 'n_iac', 'n_normal', 'auroc', 'sep_global']].to_string(index=False))
    rep['patient_centered_pooled_auroc'] = round(a_c, 6)
    rep['per_patient_centered'] = perc.replace({np.nan: None}).to_dict('records')
    rep['n_patients_below_0.5_after_centering'] = int((perc['auroc'] < 0.5).sum())

    # 翻转的患者去均值后是否回到 >0.5
    bad = set(per.loc[per['auroc'] < 0.5, 'patient'])
    fixed = {p: float(perc.loc[perc['patient'] == p, 'auroc'].iloc[0]) for p in sorted(bad)}
    rep['flipped_patients_after_centering'] = {k: round(v, 4) for k, v in fixed.items()}
    print(f'\n   原 <0.5 的患者去均值后：'
          f'{ {k: round(v,3) for k,v in fixed.items()} }')

    # ---- 4) 类别分离 vs 患者偏移：谁大 ----
    sep = per['sep_global'].values
    sh = np.array([shift[str(p)] for p in per['patient']])
    rep['median_class_separation'] = round(float(np.median(sep)), 3)
    rep['median_patient_shift'] = round(float(np.median(sh)), 3)
    rep['ratio_shift_over_separation'] = round(float(np.median(sh) / np.median(sep)), 3)
    corr = float(pd.Series(per['auroc'].values).corr(pd.Series(sh), method='spearman'))
    rep['spearman_auroc_vs_patient_shift'] = round(corr, 3)
    print(f'\n[4] 类别分离中位 {np.median(sep):.2f}  vs  患者偏移中位 {np.median(sh):.2f}'
          f'  ⇒ 比值 {np.median(sh)/np.median(sep):.3f}')
    print(f'    Spearman(逐患者 AUROC, 患者偏移) = {corr:.3f}')

    # ---- 5) 患者自己的质心落在哪边（直接看「方向」）----
    per['ownI_closer_to'] = np.where(per['d_ownI_to_cI'] < per['d_ownI_to_cN'], 'IAC质心', 'Normal质心')
    per['ownN_closer_to'] = np.where(per['d_ownN_to_cN'] < per['d_ownN_to_cI'], 'Normal质心', 'IAC质心')
    rep['own_centroid_side'] = per[['patient', 'auroc', 'ownI_closer_to',
                                    'ownN_closer_to']].to_dict('records')
    print('\n[5] 各患者「自己的 IAC 质心」离哪个全局质心更近：')
    print(per[['patient', 'auroc', 'ownI_closer_to', 'ownN_closer_to']].to_string(index=False))

    with open(f'{OUT}/s3_direction_audit.json', 'w') as fh:
        json.dump(rep, fh, ensure_ascii=False, indent=1, default=str)
    per.to_csv(f'{OUT}/s3_per_patient_original.tsv', sep='\t', index=False)
    perc.to_csv(f'{OUT}/s3_per_patient_centered.tsv', sep='\t', index=False)
    print(f'\n写到 {OUT}/')


if __name__ == '__main__':
    main()
