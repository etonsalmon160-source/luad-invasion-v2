#!/usr/bin/env python3
# HE 病理模型 · 试点第三步：全队列筛查（56 张切片）
#
# 第二步在 P4 同一患者上发现 AUC 0.86–0.88（LUAD vs Normal）。
# 但那是**单患者**，可能只是切片批次/染色差异。
# 本步扩到全部 56 张切片，看 PLIP 的 neoplastic 分数是否随分期单调。
#
# 这是**筛查**，不是判决：目的是判断"值不值得立预注册"，不产出结论。
#
# 数据：/home/eto/luad_v2/data/visium_spatial/（已从只读源拷入，56 张）
#   源目录有两套布局（19 嵌套 / 37 扁平），拷贝时已归一为 <GSM_名>/spatial/。
#
# 跑法：python3 07_he_pathology/03_cohort_screen.py

import os, json, re, warnings
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
OUT = f'{ROOT}/results/07_he_pathology/pilot'
FIGDIR = f'{OUT}/figures'
os.makedirs(FIGDIR, exist_ok=True)

CROP = 448
N_SPOT = 250
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
NEOPLASTIC = [0, 2, 3]          # adenoca + AAH + AIS
STAGE_ORDER = ['Normal', 'AAH', 'AIS', 'MIA', 'LUAD']
COLORS = {'Normal': '#1f77b4', 'AAH': '#2ca02c', 'AIS': '#ff7f0e',
          'MIA': '#d62728', 'LUAD': '#8b0000'}


def parse_dir(d):
    """目录名 → (gsm, patient, stage)。例：GSM9226178_P5_AIS。"""
    m = re.match(r'(GSM\d+)_(P\d+)_(Normal|AAH|AIS|MIA|LUAD)', d)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3)


def spatial_of(gdir):
    """本地已归一为 <dir>/spatial/；兼容万一残留的嵌套布局。"""
    if os.path.isdir(f'{gdir}/spatial'):
        return f'{gdir}/spatial'
    for inner in sorted(os.listdir(gdir)):
        p = f'{gdir}/{inner}/spatial'
        if os.path.isdir(p):
            return p
    return None


def load_positions(sp):
    pf = f'{sp}/tissue_positions.csv'
    if os.path.exists(pf):
        return pd.read_csv(pf)
    pf = f'{sp}/tissue_positions_list.csv'
    pos = pd.read_csv(pf, header=None)
    pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                   'pxl_row_in_fullres', 'pxl_col_in_fullres']
    return pos


def main():
    torch.set_num_threads(16)
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    texts = [PROMPT_TEMPLATE.format(l) for l in LABELS]

    slides = []
    for d in sorted(os.listdir(DATA)):
        p = parse_dir(d)
        if p and os.path.isdir(f'{DATA}/{d}'):
            slides.append((*p, d))
    print(f'[清单] {len(slides)} 张切片')
    for st in STAGE_ORDER:
        n = sum(1 for s in slides if s[2] == st)
        print(f'   {st:7s} {n}')
    print(f'   患者 {len({s[1] for s in slides})}')

    rows = []
    for i, (gsm, pat, stage, gdir) in enumerate(slides, 1):
        try:
            sp = spatial_of(f'{DATA}/{gdir}')
            img = Image.open(f'{sp}/tissue_hires_image.png').convert('RGB')
            pos = load_positions(sp)
            sc = json.load(open(f'{sp}/scalefactors_json.json'))['tissue_hires_scalef']
        except Exception as e:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 跳过: {type(e).__name__} {e}')
            continue

        pos = pos[pos.in_tissue == 1].copy()
        pos['r'] = pos.pxl_row_in_fullres * sc
        pos['c'] = pos.pxl_col_in_fullres * sc
        take = pos.sample(min(N_SPOT, len(pos)), random_state=SEED)

        crops = []
        for _, x in take.iterrows():
            L = max(0, min(int(x.c) - CROP // 2, img.size[0] - CROP))
            T = max(0, min(int(x.r) - CROP // 2, img.size[1] - CROP))
            c = img.crop((L, T, L + CROP, T + CROP)).resize((224, 224), Image.BICUBIC)
            if (np.asarray(c.convert('L')) <= WHITE_THRESH).mean() < MIN_TISSUE_FRAC:
                continue
            crops.append(c)
        if len(crops) < 20:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 有效裁窗太少({len(crops)})，跳过')
            continue

        probs = []
        for k in range(0, len(crops), 64):
            inp = proc(text=texts, images=crops[k:k + 64], return_tensors='pt', padding=True)
            with torch.no_grad():
                probs.append(model(**inp).logits_per_image.softmax(1).numpy())
        pr = np.vstack(probs)
        row = {'gsm': gsm, 'gsm_dir': gdir, 'patient': pat, 'stage': stage,
               'n_spot': len(crops), 'hires_w': img.size[0], 'hires_h': img.size[1],
               'neoplastic_mean': float(pr[:, NEOPLASTIC].sum(1).mean()),
               'adenoca_mean': float(pr[:, 0].mean()),
               'normal_mean': float(pr[:, 1].mean()),
               'empty_mean': float(pr[:, 6].mean()),
               'argmax_neoplastic_frac': float(np.isin(pr.argmax(1), NEOPLASTIC).mean())}
        for j, l in enumerate(LABELS):
            row[f'p_{l.replace(" ", "_")}'] = float(pr[:, j].mean())
        rows.append(row)
        print(f'  [{i:2d}/{len(slides)}] {gdir:22s} {stage:5s} n={len(crops):3d} '
              f'neoplastic={row["neoplastic_mean"]:.3f} adenoca={row["adenoca_mean"]:.3f} '
              f'normal={row["normal_mean"]:.3f}')

    df = pd.DataFrame(rows)
    df.to_csv(f'{OUT}/plip_cohort_screen.csv', index=False)
    print(f'\n[out] {OUT}/plip_cohort_screen.csv')

    # ——— 分期趋势 ———
    print('\n' + '=' * 70)
    print('各期 neoplastic 分数（adenoca+AAH+AIS）')
    for st in STAGE_ORDER:
        s = df[df.stage == st]
        if len(s) == 0:
            continue
        print(f'  {st:7s} n={len(s):2d}  mean={s.neoplastic_mean.mean():.3f} '
              f'median={s.neoplastic_mean.median():.3f}  '
              f'range=[{s.neoplastic_mean.min():.3f},{s.neoplastic_mean.max():.3f}]')
    from scipy.stats import kruskal, spearmanr, mannwhitneyu
    groups = [df[df.stage == st].neoplastic_mean.values for st in STAGE_ORDER
              if len(df[df.stage == st]) > 0]
    if len(groups) >= 2:
        h, p = kruskal(*groups)
        print(f'\n  Kruskal-Wallis H={h:.2f} p={p:.2e}（各期是否有差异）')
    idx = {s: i for i, s in enumerate(STAGE_ORDER)}
    rho, pv = spearmanr(df.stage.map(idx), df.neoplastic_mean)
    print(f'  Spearman(分期序号, neoplastic) rho={rho:.3f} p={pv:.2e}（单调性）')

    # Normal-only-n=1 ⇒ 用 AAH 作最低非正常期做对照
    if (df.stage == 'AAH').any() and (df.stage == 'LUAD').any():
        u, pu = mannwhitneyu(df[df.stage == 'LUAD'].neoplastic_mean,
                             df[df.stage == 'AAH'].neoplastic_mean, alternative='greater')
        print(f'  LUAD vs AAH（Mann-Whitney 单侧>）U={u:.0f} p={pu:.3e}')
    if (df.stage == 'Normal').any() and (df.stage == 'LUAD').any():
        u, pu = mannwhitneyu(df[df.stage == 'LUAD'].neoplastic_mean,
                             df[df.stage == 'Normal'].neoplastic_mean, alternative='greater')
        print(f'  LUAD vs Normal（n_Normal=1，无统计意义，仅参考）p={pu:.3e}')

    # ——— 图 ———
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    stages = [s for s in STAGE_ORDER if (df.stage == s).any()]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), dpi=140)
    ax = axes[0]
    rng = np.random.RandomState(0)
    for i, st in enumerate(stages):
        v = df[df.stage == st].neoplastic_mean.values
        ax.scatter(rng.normal(i, .07, len(v)), v, s=30, alpha=.8,
                   c=COLORS[st], label=f'{st} (n={len(v)})')
    ax.set_xticks(range(len(stages))); ax.set_xticklabels(stages)
    ax.set_ylabel('mean PLIP neoplastic prob. per spot')
    ax.set_title(f'PLIP neoplastic score by lesion stage\n'
                 f'all {len(df)} Visium slides (GSE307534) | Spearman rho={rho:.2f} p={pv:.1e}\n'
                 f'crop={CROP}px  n={N_SPOT}/slide')
    ax.legend(fontsize=8)
    ax = axes[1]
    for i, st in enumerate(stages):
        v = df[df.stage == st].normal_mean.values
        ax.scatter(rng.normal(i, .07, len(v)), v, s=30, alpha=.8, c=COLORS[st],
                   label=f'{st} (n={len(v)})')
    ax.set_xticks(range(len(stages))); ax.set_xticklabels(stages)
    ax.set_ylabel('mean PLIP P("normal lung tissue")')
    ax.set_title('Same slides — normal-lung prompt (should fall)')
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/plip_cohort_screen_by_stage.png', bbox_inches='tight')
    plt.close()
    print(f'[out] {FIGDIR}/plip_cohort_screen_by_stage.png')

    summary = {
        'script': '07_he_pathology/03_cohort_screen.py',
        'status': 'SCREEN — 判断是否值得立预注册，非结论',
        'model': 'vinid/plip', 'crop': CROP, 'n_spot_per_slide': N_SPOT, 'seed': SEED,
        'data': 'data/visium_spatial（自只读源拷入；源有两套布局，已归一）',
        'n_slides': int(len(df)), 'n_patients': int(df.patient.nunique()),
        'stage_of_slide_from': '目录名 (= GEO 提交方标题分期)；Normal 仅 1 张',
        'per_stage': {st: {'n': int((df.stage == st).sum()),
                           'neoplastic_mean': float(df[df.stage == st].neoplastic_mean.mean()),
                           'neoplastic_median': float(df[df.stage == st].neoplastic_mean.median()),
                           'normal_mean': float(df[df.stage == st].normal_mean.mean())}
                      for st in stages},
        'spearman_rho': float(rho), 'spearman_p': float(pv),
        'caveats': [
            'Normal 只有 1 张切片（P4_Normal）⇒ 队列级 Normal 对照 n=1，不可作统计',
            '切片级分期来自目录名/GEO 标题，未经我们冻结的病理复核',
            '同患者多病灶 ⇒ 观测不独立，未按患者聚类',
            'PLIP 是零样本模型，未在肺 FFPE 上微调',
            'deposited H&E 5.66um/px，PLIP 原生约 0.5um/px（差约 11 倍）',
            'PLIP 无法区分 AAH/AIS/MIA —— 只能给"非正常"这个粗轴',
        ],
    }
    with open(f'{OUT}/plip_cohort_screen_summary.json', 'w') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(f'[out] {OUT}/plip_cohort_screen_summary.json')


if __name__ == '__main__':
    main()
