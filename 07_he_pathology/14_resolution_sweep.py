#!/usr/bin/env python3
# HE 病理模型 · 第十四步：**分辨率判决实验**
#
# 预注册（先读它）：results/07_he_pathology/resolution_sweep/PREREG.md
# 要回答的唯一问题：把 HE 图的有效分辨率提高，判别力会不会变？
#
# 关键事实：PLIP 输入恒为 224×224 ⇒ 它实际看到的分辨率 = 裁框物理边长 / 224。
#   旧基线 hires 448 px = 2550 µm 视野 ⇒ 有效 11.4 µm/px（源图的 5.69 µm/px 白扔一半）
#   本实验 F1 aligned 224 px = 425 µm 视野 ⇒ 有效 1.90 µm/px（细 6.0 倍）
#
# 五个臂（PREREG §3）：
#   F1P0 / F1P1  aligned × 224 px，P0=原图 / P1=背景压平
#   F2P0 / F2P1  aligned × 1350 px（同 2550 µm 视野，作 F1 的对照）
#   H448         hires × 448 px —— **管线锚**，必须复现 rho ∈ [0.60, 0.78]
#                （旧基线 0.696）；落在区间外即硬停查代码，不得继续判读
#
# 五臂共用**逐字相同**的 spot 集：抽样后先算 F1 裁框的 tissue_frac，按
# MIN_TISSUE_FRAC 过滤一次（PREREG §4.1）。过滤是 spot 的属性，不是臂的属性。
#
# 两套提示词从**同一批图像向量**算出（文本编码免费）。
#
# ⚠️ 本脚本产出的是**测量值**，不是病理标签、不是恶性判定。
#
# 可断点续跑：每臂每切片单独落盘，已存在的跳过。
#
# 跑法：
#   python3 07_he_pathology/14_resolution_sweep.py            # 全部五臂
#   python3 07_he_pathology/14_resolution_sweep.py --arms H448 --slides 2   # 试跑
#   python3 07_he_pathology/14_resolution_sweep.py --arms F1P0 F1P1         # 选臂

import os, sys, json, re, time, gzip, argparse, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import torch
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
from scipy.ndimage import gaussian_filter
from transformers import CLIPModel, CLIPProcessor

ROOT = '/home/eto/luad_v2'
MODEL_DIR = f'{ROOT}/models/plip'
DATA = f'{ROOT}/data/visium_spatial'
OUT = f'{ROOT}/results/07_he_pathology/resolution_sweep'
PER_SLIDE = f'{OUT}/per_slide'
EMB = f'{OUT}/embeds'

# ——— 冻结常量（PREREG §6，一经登记不得事后调整）———
RNG_SEED = 0
N_PER_SLIDE = 1000
INPUT = 224
BATCH = 512
N_THREADS = 18          # 本机 load 长期 45–50，不拉高
MIN_TISSUE_FRAC = 0.10  # 沿用 04/11
WHITE_THRESH = 230      # 沿用 04/11
FLATTEN_SIGMA = 32      # P1 高斯尺度，单位 = aligned_tissue px ≈ 61 µm

# ——— A 套（7 条，逐字抄自 04_spot_annotation.py，与冻结的 rho=0.696 同源）———
TEMPLATE_A = 'a histopathology image of {}'
LABELS_A = [
    'lung adenocarcinoma',
    'normal lung tissue',
    'lung tissue with atypical adenomatous hyperplasia',
    'lung tissue with adenocarcinoma in situ',
    'fibrous stroma',
    'lymphoid tissue',
    'an empty glass slide',
]
NEO_A, NORM_A, NONEPI_A = [0, 2, 3], 1, [4, 5, 6]

# ——— B 套（8 条，逐字抄自 11_spot_annotation_B.py → 原文病理学家判据）———
TEMPLATE_B = 'An H&E image patch of {}.'
LABELS_B = [
    'invasive adenocarcinoma with desmoplastic stroma',
    'morphologically normal pneumocytes lining empty alveoli',
    'localized proliferation of thickened alveolar septa lined by atypical type II pneumocytes',
    'neoplastic epithelial cells continuously lining alveolar walls in a lepidic pattern',
    'monotonous enlarged cuboidal type II pneumocytes',
    'fibrotic lung tissue with collagen deposition',
    'dense lymphocytic infiltrate',
    'an empty glass slide',
]
LESION_B, NORMAL_B, NONEPI_B = [0, 2, 3], 1, [5, 6, 7]

ALL_ARMS = [
    dict(name='F1P0', src='aligned', crop=224,  flat=False),
    dict(name='F1P1', src='aligned', crop=224,  flat=True),
    dict(name='F2P0', src='aligned', crop=1350, flat=False),
    dict(name='F2P1', src='aligned', crop=1350, flat=True),
    dict(name='H448', src='hires',   crop=448,  flat=False),
]
ARM_BY_NAME = {a['name']: a for a in ALL_ARMS}

IDENT = ['gsm', 'patient', 'slide_stage', 'gsm_dir']
CORE = (['barcode', 'array_row', 'array_col',
         'pxl_row_in_fullres', 'pxl_col_in_fullres', 'tissue_frac', 'clamped']
        + [f'a{j}' for j in range(len(LABELS_A))]
        + [f'b{j}' for j in range(len(LABELS_B))]
        + ['neo_A', 'norm_A', 'none_A', 'lesion_B', 'normal_B', 'nonepi_B'])
FIELDS = IDENT + CORE   # 落盘时的列序；IDENT 由 main() 插入


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


def rgb_of(path):
    return np.asarray(Image.open(path).convert('RGB'), dtype=np.uint8)


def flatten(rgb):
    """P1 预处理：减掉自身的大尺度**亮度**背景，三通道减同一个量。

    ⚠️ 必须保留颜色。实测：把 H&E 转成灰度再复制成三通道，PLIP 会把它读成
    「空玻璃片」（empty 概率 0.888 vs 真彩 0.370）—— 颜色本身就是判别信号，
    任何会毁掉色度的预处理都会把这一臂变成另一个实验。

    三通道减同一个亮度背景 ⇒ 通道间差值(R−G, R−B)**逐像素不变** ⇒ 色相/饱和度
    完全保留，只有大尺度明暗被抹平。这正是要去掉拼接块亮度差所需的。

    背景用 1/2 降采样算再升回原尺寸：背景本就平滑，省 4 倍算力，
    且比逐裁框算少了裁框边界的镜像伪影。"""
    H, W = rgb.shape[:2]
    g = np.asarray(Image.fromarray(rgb).convert('L'), dtype=np.uint8)
    small = np.asarray(Image.fromarray(g).resize((W // 2, H // 2), Image.BILINEAR),
                       dtype=np.float32)
    bg_s = gaussian_filter(small, FLATTEN_SIGMA / 2.0)
    bg = np.asarray(Image.fromarray(bg_s).resize((W, H), Image.BILINEAR), dtype=np.float32)
    return np.clip(rgb.astype(np.float32) - bg[:, :, None] + 128.0, 0, 255).astype(np.uint8)


def crop224(g, cx, cy, crop):
    """以 (cx,cy) 为中心取 crop×crop，缩到 224。返回 PIL RGB + clamped 标志。
    与 04/11 逐字一致：裁框取自 RGB 原图，tissue_frac 在**缩放后**的图上用 'L' 量。"""
    H, W = g.shape[:2]
    if W < crop or H < crop:
        return None, 0
    L, T = int(round(cx)) - crop // 2, int(round(cy)) - crop // 2
    Lc, Tc = max(0, min(L, W - crop)), max(0, min(T, H - crop))
    clamped = int(Lc != L or Tc != T)
    sub = g[Tc:Tc + crop, Lc:Lc + crop]
    return Image.fromarray(sub).resize((INPUT, INPUT), Image.BICUBIC), clamped


def encode(model, proc, imgs, teA, teB):
    """一批 PIL RGB → 图像向量 + 两套提示词的概率。"""
    inp = proc(images=list(imgs), return_tensors='pt')
    with torch.no_grad():
        ie = model.get_image_features(pixel_values=inp['pixel_values'])
        ie = ie / ie.norm(dim=-1, keepdim=True)
        s = model.logit_scale.exp()
        pA = (s * (ie @ teA.T)).softmax(1).numpy()
        pB = (s * (ie @ teB.T)).softmax(1).numpy()
    return ie.numpy(), pA, pB


def run_slide(model, proc, teA, teB, gdir, arms_todo):
    """返回 {arm_name: (DataFrame, emb)}；已存在的臂不重算。"""
    sp = f'{DATA}/{gdir}/spatial'
    sf = json.load(open(f'{sp}/scalefactors_json.json'))
    rts, ths = sf['regist_target_img_scalef'], sf['tissue_hires_scalef']

    pos = pd.read_csv(f'{sp}/tissue_positions.csv')
    pos = pos[pos.in_tissue == 1].sort_values('barcode').reset_index(drop=True)
    if len(pos) > N_PER_SLIDE:
        rng = np.random.default_rng(RNG_SEED)
        pos = pos.iloc[np.sort(rng.choice(len(pos), N_PER_SLIDE, replace=False))].reset_index(drop=True)
    pos['cx_at'] = pos.pxl_col_in_fullres * rts
    pos['cy_at'] = pos.pxl_row_in_fullres * rts
    pos['cx_hr'] = pos.pxl_col_in_fullres * ths
    pos['cy_hr'] = pos.pxl_row_in_fullres * ths

    g_at = rgb_of(f'{sp}/aligned_tissue_image.jpg')

    # ——— 统一 spot 集：一律用 F1 裁框量 tissue_frac（PREREG §4.1）———
    tf1, keep = [], []
    for r in pos.itertuples():
        c, _ = crop224(g_at, r.cx_at, r.cy_at, 224)
        if c is None:
            continue
        f = float((np.asarray(c.convert('L')) <= WHITE_THRESH).mean())
        tf1.append(f)
        keep.append(f < MIN_TISSUE_FRAC)
    pos = pos[~np.array(keep)].reset_index(drop=True)
    tf1 = np.array([v for v, k in zip(tf1, keep) if not k], dtype=np.float32)
    if len(pos) == 0:
        return {}

    g_hr, g_flat = None, None
    out = {}
    for a in arms_todo:
        name = a['name']
        if a['src'] == 'aligned':
            g = g_at
            if a['flat']:
                if g_flat is None:
                    g_flat = flatten(g_at)
                g = g_flat
            cxs, cys = pos.cx_at.values, pos.cy_at.values
        else:
            if g_hr is None:
                g_hr = rgb_of(f'{sp}/tissue_hires_image.png')
            g, cxs, cys = g_hr, pos.cx_hr.values, pos.cy_hr.values

        crops, clamped = [], np.zeros(len(pos), dtype=np.int8)
        for i in range(len(pos)):
            c, cl = crop224(g, cxs[i], cys[i], a['crop'])
            if c is None:
                continue
            crops.append(c)
            clamped[i] = cl
        if not crops:
            out[name] = (None, None)
            continue

        embs, prA, prB = [], [], []
        for i in range(0, len(crops), BATCH):
            e, pA, pB = encode(model, proc, crops[i:i + BATCH], teA, teB)
            embs.append(e); prA.append(pA); prB.append(pB)
        prA, prB = np.vstack(prA), np.vstack(prB)

        df = pos.iloc[:len(crops)].copy()
        df['tissue_frac'] = tf1[:len(crops)]
        df['clamped'] = clamped[:len(crops)]
        for j in range(len(LABELS_A)):
            df[f'a{j}'] = prA[:, j]
        for j in range(len(LABELS_B)):
            df[f'b{j}'] = prB[:, j]
        df['neo_A'] = prA[:, NEO_A].sum(1)
        df['norm_A'] = prA[:, NORM_A]
        df['none_A'] = prA[:, NONEPI_A].sum(1)
        df['lesion_B'] = prB[:, LESION_B].sum(1)
        df['normal_B'] = prB[:, NORMAL_B]
        df['nonepi_B'] = prB[:, NONEPI_B].sum(1)
        out[name] = (df[CORE], np.vstack(embs))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--arms', nargs='+', default=[a['name'] for a in ALL_ARMS])
    ap.add_argument('--slides', type=int, default=0, help='只跑前 N 张（试跑用）')
    args = ap.parse_args()
    for n in args.arms:
        assert n in ARM_BY_NAME, f'未知臂 {n}'
    for n in args.arms:
        os.makedirs(f'{PER_SLIDE}/{n}', exist_ok=True)
        os.makedirs(f'{EMB}/{n}', exist_ok=True)

    torch.set_num_threads(N_THREADS)
    t_all = time.time()
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    with torch.no_grad():   # 文本编码只做一次（两套共 15 条），与图像无关
        tA = proc(text=[TEMPLATE_A.format(l) for l in LABELS_A],
                  return_tensors='pt', padding=True)
        tB = proc(text=[TEMPLATE_B.format(l) for l in LABELS_B],
                  return_tensors='pt', padding=True)
        teA = model.get_text_features(**tA); teA = teA / teA.norm(dim=-1, keepdim=True)
        teB = model.get_text_features(**tB); teB = teB / teB.norm(dim=-1, keepdim=True)

    slides = []
    for d in sorted(os.listdir(DATA)):
        p = parse_dir(d)
        if p and os.path.isdir(f'{DATA}/{d}'):
            slides.append((*p, d))
    if args.slides:
        slides = slides[:args.slides]

    done = skipped = failed = 0
    tot = 0
    n_slides_all = len(slides)
    filtered_log = []
    for i, (gsm, pat, stage, gdir) in enumerate(slides, 1):
        todo = [ARM_BY_NAME[n] for n in args.arms
                if not (csv_ok(f'{PER_SLIDE}/{n}/{gdir}.csv.gz')
                        and os.path.exists(f'{EMB}/{n}/{gdir}.npy'))]
        if not todo:
            skipped += 1
            continue
        try:
            res = run_slide(model, proc, teA, teB, gdir, todo)
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f'  [{i:2d}/{n_slides_all}] {gdir} 失败: {type(e).__name__} {e}', flush=True)
            failed += 1
            continue
        if not res:
            print(f'  [{i:2d}/{n_slides_all}] {gdir} 无有效 spot', flush=True)
            failed += 1
            continue
        line = []
        for name, (df, emb) in res.items():
            if df is None or len(df) == 0:
                line.append(f'{name}=0')
                continue
            np.save(f'{EMB}/{name}/{gdir}.npy', emb)
            df.insert(0, 'gsm', gsm); df.insert(1, 'patient', pat)
            df.insert(2, 'slide_stage', stage); df.insert(3, 'gsm_dir', gdir)
            df.to_csv(f'{PER_SLIDE}/{name}/{gdir}.csv.gz', index=False, compression='gzip')
            tot += len(df)
            line.append(f'{name}:n={len(df)} neoA={df.neo_A.mean():.3f} '
                        f'lesB={df.lesion_B.mean():.3f} clamp={df.clamped.mean():.2f}')
            filtered_log.append(dict(gsm_dir=gdir, arm=name, n=len(df),
                                     clamped_frac=float(df.clamped.mean())))
        done += 1
        el = time.time() - t_all
        eta = el / (done + skipped) * (n_slides_all - done - skipped) / 60 if (done + skipped) else 0
        print(f'  [{i:2d}/{n_slides_all}] {gdir:22s} {stage:5s} ' + ' | '.join(line)
              + f'  ETA {eta:.0f}min', flush=True)

    print(f'\n[done] 新跑 {done} / 跳过 {skipped} / 失败 {failed}；'
          f'{tot:,} spot；{(time.time()-t_all)/60:.1f} 分钟')

    summary = {
        'script': '07_he_pathology/14_resolution_sweep.py',
        'status': 'RAW MEASUREMENTS — 只产出测量值，不施加任何类别判定、不判恶性',
        'prereg': 'results/07_he_pathology/resolution_sweep/PREREG.md',
        'arms_run': args.arms,
        'model': 'vinid/plip (CLIP ViT-B/32)',
        'frozen': dict(RNG_SEED=RNG_SEED, N_PER_SLIDE=N_PER_SLIDE, INPUT=INPUT,
                       BATCH=BATCH, N_THREADS=N_THREADS,
                       MIN_TISSUE_FRAC=MIN_TISSUE_FRAC, WHITE_THRESH=WHITE_THRESH,
                       FLATTEN_SIGMA=FLATTEN_SIGMA,
                       template_A=TEMPLATE_A, labels_A=LABELS_A,
                       template_B=TEMPLATE_B, labels_B=LABELS_B),
        'arms': [{k: v for k, v in a.items()} for a in
                 [ARM_BY_NAME[n] for n in args.arms]],
        'geometry': {
            'F1': 'aligned 224px @1.897um/px = 425um 视野 ⇒ 1.90 um/px',
            'F2': 'aligned 1350px @1.897um/px = 2560um 视野 ⇒ 11.4 um/px',
            'H448': 'hires 448px @5.69um/px = 2550um 视野 ⇒ 11.4 um/px（旧基线）',
        },
        'spot_filter_rule': '一律用 F1 裁框（aligned 224px）量 tissue_frac，五臂同一 spot 集',
        'color_caveat': 'H&E 必须保留颜色：实测把裁框转灰度再复制成三通道后，'
                        'PLIP 把它读成「空玻璃片」（empty 0.888 vs 真彩 0.370，'
                        'neoplastic 0.038 vs 0.221）。P1 的压平是三通道减同一个亮度背景，'
                        '通道间差值逐像素不变 ⇒ 色相保留。',
        'pending_reading': '判读由 15_resolution_sweep_analysis.py 出；H448 rho 必须 ∈ [0.60,0.78]',
        'filtered_log': filtered_log[:200],
    }
    with open(f'{OUT}/resolution_sweep_summary.json', 'w') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(f'[out] {OUT}/resolution_sweep_summary.json')


if __name__ == '__main__':
    main()
