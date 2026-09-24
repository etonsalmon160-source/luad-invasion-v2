#!/usr/bin/env python3
# HE 病理模型 · 图像源三臂对照
#
# 问题：PLIP 一直喂的是 tissue_hires_image.png（中位 6.21 µm/px）。队列里有更细的
#       cytassist_image.tiff（4.6233 µm/px，彩色，1.343× 细）。换更细的源，分数会变吗？
#
# 四臂（同一批 spot、同一套提示词、裁框按 µm 对齐）：
#   A  hires                                    —— 复现既有基线（rho 0.695）
#   B  cytassist                                —— 换成更细的源
#   C  cytassist，对比度对齐到 hires            —— 把"分辨率"和"对比度"分开
#   D  cytassist 降到 hires 像素密度            —— 把"扫描本身"和"分辨率"分开
#
# 为什么必须有 C：两个源除了分辨率还差对比度/白平衡。不控制这一项，B-A 分不清是谁造成的。
#
# 为什么还要 D：PLIP 输入固定 224×224。B 和 A 覆盖同一物理视野 ⇒ 进网络时都是
#   ~14.6 µm/px，cytassist 的细采样优势在降采样那一步基本被丢掉。D 把 cytassist
#   先降到 hires 的像素密度再进 224：若 D≈B，则"更细的源"没起作用；差异来自扫描本身。
#
# 方向对齐：cytassist 相对 hires 是**镜像**的（53/56 张）。B/C 臂用仿射 transform 把
#   cytassist 采样回 hires 的朝向与中心，所以三臂看的是同一块组织、同一个朝向。
#
# 出帧（cytassist 3000² 装不下）的 spot 三臂一起丢，保持配对。
#
# 跑法：python3 07_he_pathology/17_plip_source_arms.py [GSM_P_Stage ...]

import os, json, re, csv, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import torch
from PIL import Image
from scipy import ndimage
Image.MAX_IMAGE_PIXELS = None
from transformers import CLIPModel, CLIPProcessor

ROOT = '/home/eto/luad_v2'
MODEL_DIR = f'{ROOT}/models/plip'
DATA = f'{ROOT}/data/visium_spatial'
RES = f'{ROOT}/results/07_he_pathology/resolution_sweep'
OUT = f'{ROOT}/results/07_he_pathology/source_arms'
os.makedirs(OUT, exist_ok=True)

# —— 与 03_cohort_screen.py 逐字一致 ——
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
NEOPLASTIC = [0, 2, 3]
STAGE_ORDER = ['Normal', 'AAH', 'AIS', 'MIA', 'LUAD']

FULLRES_UM = 0.25050
CYTO_UM = 4.6233

AFF = {r['slide']: r for r in csv.DictReader(open(f'{RES}/spot_affine.csv'))}
CSVF = ''


def parse_dir(d):
    m = re.match(r'(GSM\d+)_(P\d+)_(Normal|AAH|AIS|MIA|LUAD)', d)
    return m.groups() if m else None


def spatial_of(gdir):
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
    pos = pd.read_csv(f'{sp}/tissue_positions_list.csv', header=None)
    pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                   'pxl_row_in_fullres', 'pxl_col_in_fullres']
    return pos


def match_hist(src, ref):
    """把 src 的灰度直方图拉到 ref 上（C 臂：给 cytassist 装上 hires 的对比度）。"""
    s = np.clip(src, 0, 255).astype(np.uint8).ravel()
    r = np.clip(ref, 0, 255).astype(np.uint8).ravel()
    sh = np.bincount(s, minlength=256).cumsum().astype(np.float64); sh /= sh[-1]
    rh = np.bincount(r, minlength=256).cumsum().astype(np.float64); rh /= rh[-1]
    rh = rh + np.arange(256) * 1e-9
    lut = np.interp(sh, rh, np.arange(256))
    return lut[np.clip(src, 0, 255).astype(np.uint8)]


def match_hist_rgb(src, ref):
    """保色相的对比度对齐：算亮度 LUT，再对每个像素用同一个标量缩放三通道。
    逐通道对齐会把色相搞坏（出青绿色伪影），那种伪影自己就会 drive 模型。"""
    gl = src.mean(-1)
    rl = ref.mean(-1)
    sh = np.bincount(np.clip(gl, 0, 255).astype(np.uint8).ravel(), minlength=256).cumsum().astype(np.float64)
    sh /= sh[-1]
    rh = np.bincount(np.clip(rl, 0, 255).astype(np.uint8).ravel(), minlength=256).cumsum().astype(np.float64)
    rh /= rh[-1]
    lut = np.interp(sh, rh + np.arange(256) * 1e-9, np.arange(256))
    new_lum = lut[np.clip(gl, 0, 255).astype(np.uint8)]
    scale = new_lum / np.maximum(gl, 20.0)
    return np.clip(src * scale[..., None], 0, 255).astype(np.uint8)


def text_features(model, proc, texts):
    with torch.no_grad():
        e = model.get_text_features(**proc(text=texts, return_tensors='pt', padding=True))
    return e / e.norm(dim=-1, keepdim=True)


def score(model, proc, temb, crops):
    """与 03 的 logits_per_image.softmax(1) 数值等价，但文本只编码一次。"""
    scale = model.logit_scale.exp()
    ps = []
    with torch.no_grad():
        for k in range(0, len(crops), 64):
            e = model.get_image_features(**proc(images=crops[k:k + 64], return_tensors='pt'))
            e = e / e.norm(dim=-1, keepdim=True)
            ps.append((scale * e @ temb.T).softmax(1).numpy())
    return np.vstack(ps)


def row_of(arm, base, pr, dens):
    r = dict(base)
    r['arm'] = arm
    r['n_spot'] = len(pr)
    r['neoplastic_mean'] = float(pr[:, NEOPLASTIC].sum(1).mean())
    r['adenoca_mean'] = float(pr[:, 0].mean())
    r['normal_mean'] = float(pr[:, 1].mean())
    r['empty_mean'] = float(pr[:, 6].mean())
    r['argmax_neoplastic_frac'] = float(np.isin(pr.argmax(1), NEOPLASTIC).mean())
    for j, l in enumerate(LABELS):
        r[f'p_{l.replace(" ", "_")}'] = float(pr[:, j].mean())
    if len(dens) == len(pr) and np.std(dens) > 0:
        from scipy.stats import spearmanr
        r['dens_spearman'] = float(spearmanr(dens, pr[:, NEOPLASTIC].sum(1)).correlation)
    else:
        r['dens_spearman'] = np.nan
    return r


def one(model, proc, temb, gdir, gsm, pat, stage, i, n):
    sp = spatial_of(f'{DATA}/{gdir}')
    img = Image.open(f'{sp}/tissue_hires_image.png').convert('RGB')
    pos = load_positions(sp)
    sc = json.load(open(f'{sp}/scalefactors_json.json'))['tissue_hires_scalef']
    um_h = FULLRES_UM / sc
    ncy = int(round(CROP * um_h / CYTO_UM))
    a = AFF[gdir]
    A = np.array([[float(a['a00']), float(a['a01'])], [float(a['a10']), float(a['a11'])]])
    t = np.array([float(a['t0']), float(a['t1'])])
    M = (CYTO_UM / um_h) * A
    c0 = (ncy - 1) / 2.0

    pos = pos[pos.in_tissue == 1].copy()
    pos['r'] = pos.pxl_row_in_fullres * sc
    pos['c'] = pos.pxl_col_in_fullres * sc
    take = pos.sample(min(N_SPOT, len(pos)), random_state=SEED)

    Cyto = None
    CA, CB, CC, CD, dens, dots = [], [], [], [], [], []
    n_out = 0
    for _, x in take.iterrows():
        L = max(0, min(int(x.c) - CROP // 2, img.size[0] - CROP))
        T = max(0, min(int(x.r) - CROP // 2, img.size[1] - CROP))
        ca = img.crop((L, T, L + CROP, T + CROP))
        if (np.asarray(ca.convert('L')) <= WHITE_THRESH).mean() < MIN_TISSUE_FRAC:
            continue
        if Cyto is None:
            Cyto = np.asarray(Image.open(f'{sp}/cytassist_image.tiff').convert('RGB'), np.float32)
        ctr = A @ np.array([T + CROP / 2.0, L + CROP / 2.0]) + t
        off = ctr - M @ np.array([c0, c0])
        M3 = np.array([[M[0, 0], M[0, 1], 0.0], [M[1, 0], M[1, 1], 0.0], [0.0, 0.0, 1.0]])
        cb = ndimage.affine_transform(Cyto, M3, [off[0], off[1], 0.0],
                                      output_shape=(ncy, ncy, 3), order=1,
                                      mode='constant', cval=np.nan)
        if not np.isfinite(cb).all():
            n_out += 1
            continue                      # 出帧 -> 三臂一起丢，保持配对
        ca224 = ca.resize((224, 224), Image.BICUBIC)
        cb224 = Image.fromarray(np.clip(cb, 0, 255).astype(np.uint8)).resize((224, 224), Image.BICUBIC)
        cbc = match_hist_rgb(cb, np.asarray(ca, np.float32))
        cc224 = Image.fromarray(cbc).resize((224, 224), Image.BICUBIC)
        # D: cytassist 降到 hires 的像素密度再进 224 —— 与 B 同源同视野，只差"分辨率"
        cb448 = Image.fromarray(np.clip(cb, 0, 255).astype(np.uint8)).resize((CROP, CROP), Image.BICUBIC)
        cd224 = cb448.resize((224, 224), Image.BICUBIC)
        CA.append(ca224); CB.append(cb224); CC.append(cc224); CD.append(cd224)
        dens.append(float((np.asarray(ca224.convert('L')) <= WHITE_THRESH).mean()))
        # cytassist 图里烧录了点阵/fiducial（黑点），记录污染程度供敏感性分析
        dots.append(float((np.asarray(cb224.convert('L')) < 90).mean()))

    if len(CA) < 20:
        print(f'  [{i:2d}/{n}] {gdir} 有效裁窗太少({len(CA)})，跳过')
        return []
    del Cyto
    base = dict(gsm=gsm, gsm_dir=gdir, patient=pat, stage=stage,
                um_hires=round(um_h, 4), ncyto=ncy, n_outframe=n_out,
                dot_frac_mean=round(float(np.mean(dots)), 5),
                dot_frac_max=round(float(np.max(dots)), 5))
    out = []
    for arm, crops in (('A_hires', CA), ('B_cytassist', CB),
                       ('C_cytassist_contrasthmatched', CC),
                       ('D_cytassist_at_hires_pixel_density', CD)):
        pr = score(model, proc, temb, crops)
        out.append(row_of(arm, base, pr, dens))
    print(f'  [{i:2d}/{n}] {gdir:24s} {stage:5s} n={len(CA):3d} 出帧 {n_out:2d} | '
          + ' | '.join(f'{r["arm"][0]} {r["neoplastic_mean"]:.3f}' for r in out))
    sys.stdout.flush()
    return out


def main():
    torch.set_num_threads(int(os.environ.get('NTHREADS', '16')))
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    texts = [PROMPT_TEMPLATE.format(l) for l in LABELS]
    temb = text_features(model, proc, texts)

    slides = []
    for d in sorted(os.listdir(DATA)):
        p = parse_dir(d)
        if p and os.path.isdir(f'{DATA}/{d}') and d in AFF:
            slides.append((*p, d))
    if len(sys.argv) > 1:
        slides = [s for s in slides if s[3] in sys.argv[1:]]
    print(f'[清单] {len(slides)} 张')

    global CSVF
    tag = os.environ.get('SHARD_TAG', '')
    CSVF = f'{OUT}/_arms_{tag}.csv' if tag else f'{OUT}/plip_source_arms.csv'
    rows = []
    for i, (gsm, pat, stage, gdir) in enumerate(slides, 1):
        try:
            rows += one(model, proc, temb, gdir, gsm, pat, stage, i, len(slides))
        except Exception as e:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 出错: {type(e).__name__} {e}')
        df = pd.DataFrame(rows)
        if len(df):
            df.to_csv(CSVF, index=False)
    print(f'\n[out] {CSVF}  ({len(rows)} 行)')


if __name__ == '__main__':
    main()
