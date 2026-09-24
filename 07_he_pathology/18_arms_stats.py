#!/usr/bin/env python3
"""四臂图像源对照 · 合并分片 + 统计 + 出图。

输入：results/07_he_pathology/source_arms/_arms_{0..3}.csv（17_plip_source_arms.py 的分片输出）
输出：同目录 plip_source_arms.csv（合并）、arms_stats.json、arms_by_stage.png

三个要回答的问题：
  Q1  A 臂能不能复现既有基线（rho ≈ 0.695）？——不能复现则整场对照不成立。
  Q2  换更细的源（B）分数变不变？变了是**谁**造成的（C 控对比度 / D 控分辨率）？
  Q3  差异在切片间稳不稳（配对检验），还是被少数切片带偏？

图内文字全英文（本机 matplotlib 无中文字体）。
用法：python3 07_he_pathology/18_arms_stats.py
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon

OUT = '/home/eto/luad_v2/results/07_he_pathology/source_arms'
STAGE_NUM = {'Normal': 0, 'AAH': 1, 'AIS': 2, 'MIA': 3, 'LUAD': 4}
ARMS = ['A_hires', 'B_cytassist', 'C_cytassist_contrasthmatched',
        'D_cytassist_at_hires_pixel_density']
SHORT = {'A_hires': 'A hires', 'B_cytassist': 'B cytassist',
         'C_cytassist_contrasthmatched': 'C cyto+contrast',
         'D_cytassist_at_hires_pixel_density': 'D cyto@hires px'}
BASELINE_RHO = 0.695          # 03_cohort_screen.py 的既有基线，A 臂应当复现


EXPECT_PER_SHARD = 14          # 17 号脚本按 4 分片各 14 张切的；不等于这个数就是没跑完


def merge():
    parts = []
    for i in range(4):
        try:
            parts.append(pd.read_csv(f'{OUT}/_arms_{i}.csv'))
        except FileNotFoundError:
            print(f'  [警告] 缺分片 {i}')
    if not parts:
        raise SystemExit('[硬停] 一个分片都没有')
    # 🔴 防呆：分片是**边跑边写**的，跑到一半的表看起来完全正常，只是切片少。
    #    2026-09-25 我就用一张 15 切片的半成品算出了 rho=0.793 并差点当结果报出去。
    #    没跑满就绝不合并、绝不出统计 —— 半成品冒充结果的代价远高于多等一小时。
    short = [f'_arms_{i}.csv={len(p) // 4} 张' for i, p in enumerate(parts)
             if len(p) // 4 != EXPECT_PER_SHARD]
    if short:
        raise SystemExit(f'[硬停] 分片未跑满（应各 {EXPECT_PER_SHARD} 张）：{"; ".join(short)}\n'
                         f'        分片还在写就别跑统计；跑完了还缺就是真丢片，先查日志。')
    df = pd.concat(parts, ignore_index=True)
    df['stage_num'] = df.stage.map(STAGE_NUM)
    df = df.sort_values(['arm', 'stage_num', 'gsm_dir']).reset_index(drop=True)
    df.to_csv(f'{OUT}/plip_source_arms.csv', index=False)
    return df


def rho_of(g):
    return float(spearmanr(g.stage_num, g.neoplastic_mean).correlation)


def one_arm(g):
    r = {'n_slide': int(len(g)), 'n_spot': int(g.n_spot.sum()),
         'rho_neoplastic_vs_stage': rho_of(g),
         'neoplastic_mean': float(g.neoplastic_mean.mean()),
         'argmax_neoplastic_frac': float(g.argmax_neoplastic_frac.mean()),
         'empty_mean': float(g.empty_mean.mean())}
    for s in STAGE_NUM:
        sub = g[g.stage == s]
        r[f'neoplastic_mean_{s}'] = float(sub.neoplastic_mean.mean()) if len(sub) else np.nan
    return r


def main():
    df = merge()
    arms = {a: df[df.arm == a] for a in ARMS}
    res = {}
    print(f'{"臂":26s} {"切片":>4s} {"spot":>6s} {"rho":>7s} {"病变均值":>8s} {"argmax病变":>10s}')
    for a in ARMS:
        r = one_arm(arms[a])
        res[a] = r
        print(f'{SHORT[a]:26s} {r["n_slide"]:4d} {r["n_spot"]:6d} '
              f'{r["rho_neoplastic_vs_stage"]:7.3f} {r["neoplastic_mean"]:8.3f} '
              f'{r["argmax_neoplastic_frac"]:10.3f}')
        res[a]['label'] = SHORT[a]

    # Q1 基线复现
    ra = res['A_hires']['rho_neoplastic_vs_stage']
    res['baseline_check'] = {'expected_rho': BASELINE_RHO, 'got_rho': ra,
                             'reproduced': bool(abs(ra - BASELINE_RHO) < 0.05)}
    print(f'\nQ1 基线复现：期望 {BASELINE_RHO:.3f}，A 臂实得 {ra:.3f} ⇒ '
          f'{"复现" if res["baseline_check"]["reproduced"] else "未复现（差异需解释）"}')

    # Q3 配对检验：同一批切片上 A 与各臂逐切片配对
    piv = {a: arms[a].set_index('gsm_dir').neoplastic_mean for a in ARMS}
    common = sorted(set.intersection(*[set(p.index) for p in piv.values()]))
    print(f'\nQ3 配对检验（{len(common)} 张共同切片，Wilcoxon 符号秩，双侧）')
    for a in ARMS[1:]:
        x = piv['A_hires'][common].values
        y = piv[a][common].values
        d = y - x
        try:
            st, p = wilcoxon(x, y)
        except ValueError:
            st, p = np.nan, np.nan
        # ⚠️ 键名用 `arm` 不用 `B`：本循环比的是 A 与**每一个**臂（B/C/D），
        #    原先写死 `n_slide_B_gt_A` 会让 A_vs_C / A_vs_D 两块里的键名谎报自己是 B
        #    （2026-09-25 更正；数值一直是对的，只有键名误导）
        res[f'A_vs_{a}'] = {'n_paired': len(common), 'median_delta': float(np.median(d)),
                            'p': float(p),
                            'arm': a,
                            'n_slide_arm_gt_A': int((d > 0).sum()),
                            'n_slide_arm_lt_A': int((d < 0).sum())}
        print(f'  A vs {SHORT[a]:20s} 中位差 {np.median(d):+.3f}  '
              f'臂>A {int((d>0).sum()):2d} 张 / 臂<A {int((d<0).sum()):2d} 张  p={p:.2e}')

    # 敏感性：`cytassist_image.tiff` 上烧录的 spot 点阵会污染裁窗。逐切片记录了污染率，
    # 剔掉污染最重的切片看 rho 是否还稳（只有切片级污染率，故只能做切片级剔除）。
    clean = df[df.dot_frac_max <= 0.05]
    sens = {a: rho_of(clean[clean.arm == a]) for a in ARMS}
    res['sensitivity_dot_clean'] = {'rule': 'dot_frac_max <= 0.05', **sens,
                                    'n_slide_per_arm': int(len(clean[clean.arm == a]))}
    print(f'\n敏感性（剔掉点阵污染最重的切片，每臂剩 '
          f'{res["sensitivity_dot_clean"]["n_slide_per_arm"]} 张）：')
    for a in ARMS:
        print(f'  {SHORT[a]:26s} rho {res[a]["rho_neoplastic_vs_stage"]:.3f} -> {sens[a]:.3f}')

    json.dump(res, open(f'{OUT}/arms_stats.json', 'w'), indent=2, ensure_ascii=False)

    # 图：左=逐切片配对（A vs B/D），右=各期均值
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for a in ['B_cytassist', 'D_cytassist_at_hires_pixel_density']:
        ax[0].scatter(piv['A_hires'][common], piv[a][common], s=22, alpha=.7,
                      label=SHORT[a], marker='o' if a.startswith('B') else '^')
    lim = [0, max(piv['A_hires'].max(), piv['B_cytassist'].max()) * 1.05]
    ax[0].plot(lim, lim, 'k--', lw=.8)
    ax[0].set_xlabel('Arm A (hires) neoplastic prob.')
    ax[0].set_ylabel('Arm B/D neoplastic prob.')
    ax[0].set_title('Per-slide, paired (dashed = equal)')
    ax[0].legend(fontsize=8)
    w = 0.2
    stages = list(STAGE_NUM)
    for k, a in enumerate(ARMS):
        g = arms[a]
        m = [g[g.stage == s].neoplastic_mean.mean() for s in stages]
        ax[1].bar(np.arange(len(stages)) + (k - 1.5) * w, m, w, label=SHORT[a])
    ax[1].set_xticks(range(len(stages)))
    ax[1].set_xticklabels(stages)
    ax[1].set_ylabel('Mean neoplastic prob.')
    ax[1].set_title('By stage')
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f'{OUT}/arms_by_stage.png', dpi=130)
    print(f'\n[out] {OUT}/plip_source_arms.csv, arms_stats.json, arms_by_stage.png')


if __name__ == '__main__':
    main()
