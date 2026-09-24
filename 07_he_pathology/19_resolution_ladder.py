#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""HE 病理模型 · 分辨率判决实验（PREREG_v2 的执行脚本）

预注册（先读它）：results/07_he_pathology/resolution_sweep/PREREG_v2.md
要回答的唯一问题：**把 HE 图的有效分辨率提高，判别力会不会变？**

关键事实：PLIP 输入恒为 224×224 ⇒ 它实际看到的分辨率 = 裁框物理边长 / 224。
    L448 = 448 px @ 5.670 µm/px = 2540 µm 视野 ⇒ 有效 11.34 µm/px（源图一半白扔）
    L320 = 320 px                  = 1814 µm 视野 ⇒ 有效  8.10 µm/px
    L224 = 224 px                  = 1270 µm 视野 ⇒ 有效  5.670 µm/px（1:1，不重采样）
「让 PLIP 的 224 个输入像素覆盖更少的组织」—— 这就是本实验能做到的**最细**档。

六个臂 = 三档裁框 × 两种预处理（P0 原图 / P1 背景压平）：
    L448P0 ★ = **管线锚**，与旧基线（04/11）同源同裁框同提示词，只是 spot 抽得多（1000）
               判据：切片级 rho ∈ [0.60, 0.78]（旧值 0.696）—— 越界即硬停查代码
    L448P1 / L320P0 / L320P1 / L224P0 / L224P1

两个阶梯的判读分工（PREREG_v2 §4）：
    P0 阶梯（L448P0→L320P0→L224P0）：细节加进来，判别力有没有变化？
    P1 阶梯（同三档）              ：剥掉大尺度明暗后，细尺度上还剩不剩判别力？
                                     （H&E 线此前被判定更像"组织密度读数"，P1 就是去混杂那把刀）

⚠️ 架构天花板（PREREG_v2 §2.1）：ViT-B/32 ⇒ 224 px → 7×7 token ⇒ 一个 token 覆盖 32 输入像素
   ⇒ L448 一个 token = 363 µm、L224 = 181 µm。**PLIP 结构上表达不了比 ~180 µm 更细的东西。**
   ⇒ 平坦结果的正确措辞是「在 5.670 µm/px 以上，PLIP 的判别力不随细节变化」，
     **不是**「细节对病理判别无用」。

⚠️ 本脚本产出的是**测量值**，不是病理标签、不是 WHO 分型、不是恶性判定。

关于 P1 的一处实现澄清（**登记**，见 PREREG_v2 §4 与本脚本 manifest）：
    预注册写的是「裁框转灰度得亮度背景 bg」。本脚本把背景算在**整张切片**上（每次一张，
    再取裁框对应的那一块），不是逐裁框各算一份。两个理由：① 逐裁框算背景时，
    裁框自身的边界会让背景在边缘被拉偏（滤波器截断在数组边上），而背景本该是
    「这一带的大尺度明暗」，不该依赖裁框摆在哪；② 本仓库已验证过的实现
    （`14_resolution_sweep.py` 的 flatten）就是整图算一次，直接沿用同一条代码路径。
    物理尺度不变：FLATTEN_SIGMA=32 **hires px** ≈ 181 µm。

可断点续跑：每臂每切片单独落盘，已存在且 gz 能读完的跳过。

跑法：
  python3 07_he_pathology/19_resolution_ladder.py --limit 1        # G0 成本实测
  python3 07_he_pathology/19_resolution_ladder.py                  # 全量（用户已拍板）
  python3 07_he_pathology/19_resolution_ladder.py --arms L224P0 L224P1
"""

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
OUT = f'{ROOT}/results/07_he_pathology/resolution_ladder'
PER_SLIDE = f'{OUT}/per_slide'
EMB = f'{OUT}/embeds'

# ——— 冻结常量（PREREG_v2 §7，一经登记不得事后调整）———
RNG_SEED = 0
N_PER_SLIDE = 1000
CROPS = [448, 320, 224]        # hires px
INPUT = 224                    # PLIP 输入
BATCH = 512
N_THREADS = 18                 # 本机 load 长期 45–50，不拉高
MIN_TISSUE_FRAC = 0.10         # 沿用 04/11
WHITE_THRESH = 230             # 沿用 04/11
FLATTEN_SIGMA = 32             # P1 的高斯尺度（hires px ≈ 181 µm）
MIN_CHROMA_STD = 5.0           # G2 彩色守卫（PREREG_v2 §10）

# hires 的 µm/px 由 fullres 0.2505 ÷ tissue_hires_scalef 得出（逐切片实测 4.41–7.43，中位 6.21）
FULLRES_UM = 0.25050

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


def arms_all():
    """六臂：三档裁框 × 两种预处理。L448P0 是锚（名字里的顺序**不得**改）。"""
    out = []
    for c in CROPS:
        for p, flat in (('P0', False), ('P1', True)):
            out.append(dict(name=f'L{c}{p}', crop=c, flat=flat))
    return out


ALL_ARMS = arms_all()
ARM_BY_NAME = {a['name']: a for a in ALL_ARMS}
ANCHOR = 'L448P0'

IDENT = ['gsm', 'patient', 'slide_stage', 'gsm_dir']
CORE = (['barcode', 'array_row', 'array_col',
         'pxl_row_in_fullres', 'pxl_col_in_fullres',
         'tissue_frac', 'tissue_frac_224', 'tissue_frac_arm', 'clamped', 'chroma_std']
        + [f'a{j}' for j in range(len(LABELS_A))]
        + [f'b{j}' for j in range(len(LABELS_B))]
        + ['neo_A', 'norm_A', 'none_A', 'lesion_B', 'normal_B', 'nonepi_B'])


def parse_dir(d):
    m = re.match(r'(GSM\d+)_(P\d+)_(Normal|AAH|AIS|MIA|LUAD)', d)
    return (m.group(1), m.group(2), m.group(3)) if m else None


def csv_ok(p):
    """整条 gz 流能读完才算已完成：被杀在写盘中途会留下没有校尾的半截文件，
    只看「文件在不在」会把它当成已完成永久跳过（静默坏数据）。"""
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


def background_rgb(rgb):
    """P1 的背景：整图灰度的大尺度分量，三通道共用。

    ⚠️ 必须保留颜色。实测：把 H&E 转灰度再复制成三通道，PLIP 会把它读成
    「空玻璃片」（empty 0.888 vs 真彩 0.370，neoplastic 0.038 vs 0.221）——
    颜色本身就是判别信号，任何毁掉色度的预处理都会把这一臂变成另一个实验。
    三通道减同一个亮度背景 ⇒ 通道间差值（R−G、R−B）逐像素不变 ⇒ 色相/饱和度完整保留。

    背景用 1/2 降采样算再升回：背景本就平滑，省 4 倍算力（sigma 同步减半，物理尺度不变）。
    """
    H, W = rgb.shape[:2]
    g = np.asarray(Image.fromarray(rgb).convert('L'), dtype=np.uint8)
    small = np.asarray(Image.fromarray(g).resize((W // 2, H // 2), Image.BILINEAR),
                       dtype=np.float32)
    bg_s = gaussian_filter(small, FLATTEN_SIGMA / 2.0)
    return np.asarray(Image.fromarray(bg_s).resize((W, H), Image.BILINEAR), dtype=np.float32)


def crop224(g, cx, cy, crop):
    """以 (cx,cy) 为中心取 crop×crop，缩到 224。返回 PIL RGB + clamped 标志 + 原始裁框。

    ⚠️ 与 04/11 逐字一致：裁框取自 RGB 原图，tissue_frac 在**缩放后**的图上用 'L' 量。
    """
    H, W = g.shape[:2]
    if W < crop or H < crop:
        return None, 0, None
    L, T = int(round(cx)) - crop // 2, int(round(cy)) - crop // 2
    Lc, Tc = max(0, min(L, W - crop)), max(0, min(T, H - crop))
    clamped = int(Lc != L or Tc != T)
    sub = g[Tc:Tc + crop, Lc:Lc + crop]
    return Image.fromarray(sub).resize((INPUT, INPUT), Image.BICUBIC), clamped, sub


def tissue_frac_of(pil):
    return float((np.asarray(pil.convert('L')) <= WHITE_THRESH).mean())


def chroma_std_of(pil):
    """G2 守卫用的彩色度：std(R−G)。灰度图的 R=G ⇒ 恒为 0。
    在**送进网络的 224 图**上量（守卫要防的是「喂给模型的图不是彩色的」）。"""
    a = np.asarray(pil, dtype=np.float32)
    return float((a[:, :, 0] - a[:, :, 1]).std())


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


def run_slide(model, proc, teA, teB, gdir, arms_todo, log):
    """返回 {arm_name: (DataFrame, emb)}。

    🔴 六臂共用**逐字相同**的 spot 集：抽样一次（1000 个 in-tissue spot），
       过滤也只用 L448 裁框量 tissue_frac（PREREG_v2 §5）——
       **过滤是 spot 的属性，不是臂的属性**。这样六臂可直接配对比较。
    """
    sp = f'{DATA}/{gdir}/spatial'
    sf = json.load(open(f'{sp}/scalefactors_json.json'))
    um_h = FULLRES_UM / sf['tissue_hires_scalef']

    pos = pd.read_csv(f'{sp}/tissue_positions.csv')
    pos = pos[pos.in_tissue == 1].sort_values('barcode').reset_index(drop=True)
    n_pool = len(pos)
    if n_pool > N_PER_SLIDE:
        rng = np.random.default_rng(RNG_SEED)
        pos = pos.iloc[np.sort(rng.choice(n_pool, N_PER_SLIDE, replace=False))].reset_index(drop=True)
    pos['cx'] = pos.pxl_col_in_fullres * sf['tissue_hires_scalef']
    pos['cy'] = pos.pxl_row_in_fullres * sf['tissue_hires_scalef']

    g_hr = rgb_of(f'{sp}/tissue_hires_image.png')

    # ——— 统一 spot 集：一律用 L448 裁框量 tissue_frac（PREREG_v2 §5）———
    tf448, tf224, clamp448, kept = [], [], [], []
    for r in pos.itertuples():
        c4, cl4, _ = crop224(g_hr, r.cx, r.cy, 448)
        if c4 is None:
            continue
        f4 = tissue_frac_of(c4)
        c2, _, _ = crop224(g_hr, r.cx, r.cy, 224)     # 只登记，不过滤（§5 的 ⚠️）
        tf448.append(f4)
        tf224.append(tissue_frac_of(c2) if c2 is not None else np.nan)
        clamp448.append(cl4)
        kept.append(f4 >= MIN_TISSUE_FRAC)
    kept = np.array(kept, dtype=bool)
    n_sampled = len(kept)
    n_kept = int(kept.sum())
    pos = pos[kept].reset_index(drop=True)
    if len(pos) == 0:
        log.update(n_pool=n_pool, n_sampled=n_sampled, n_kept=0)
        return {}
    tf448 = np.array(tf448, dtype=np.float32)[kept]
    tf224 = np.array(tf224, dtype=np.float32)[kept]
    clamp448 = np.array(clamp448, dtype=np.int8)[kept]
    log.update(n_pool=n_pool, n_sampled=n_sampled, n_kept=n_kept,
               um_hires=round(float(um_h), 4),
               n_clamped_448=int(clamp448.sum()),
               n_filtered_out=int(n_sampled - n_kept))

    g_flat, bg = None, None
    out = {}
    for a in arms_todo:
        crop, name = a['crop'], a['name']
        if a['flat']:
            if bg is None:
                bg = background_rgb(g_hr)
                g_flat = np.clip(g_hr.astype(np.float32) - bg[:, :, None] + 128.0,
                                 0, 255).astype(np.uint8)
            g = g_flat
        else:
            g = g_hr

        crops, clamped, chroma = [], np.zeros(len(pos), dtype=np.int8), []
        tf_arm = []
        for i in range(len(pos)):
            c, cl, _ = crop224(g, pos.cx.values[i], pos.cy.values[i], crop)
            if c is None:
                continue
            crops.append(c)
            clamped[i] = cl
            chroma.append(chroma_std_of(c))
            tf_arm.append(tissue_frac_of(c))
        if not crops:
            out[name] = (None, None)
            continue

        # ——— G2 彩色守卫（PREREG_v2 §10）：任一臂不满足 ⇒ 该臂不出数 ———
        med = float(np.median(chroma))
        if med <= MIN_CHROMA_STD:
            print(f'  🔴 G2 硬停：{name} 的 std(R−G) 中位 {med:.2f} ≤ {MIN_CHROMA_STD} '
                  f'⇒ 这一臂喂进去的不是彩色图，**不出数**', flush=True)
            log.setdefault('g2_failed', []).append({'arm': name, 'chroma_median': med})
            out[name] = (None, None)
            continue

        embs, prA, prB = [], [], []
        for i in range(0, len(crops), BATCH):
            e, pA, pB = encode(model, proc, crops[i:i + BATCH], teA, teB)
            embs.append(e); prA.append(pA); prB.append(pB)
        prA, prB = np.vstack(prA), np.vstack(prB)

        df = pos.iloc[:len(crops)].copy()
        df['tissue_frac'] = tf448[:len(crops)]
        df['tissue_frac_224'] = tf224[:len(crops)]
        df['tissue_frac_arm'] = np.array(tf_arm, dtype=np.float32)
        df['clamped'] = clamped[:len(crops)]
        df['chroma_std'] = np.array(chroma, dtype=np.float32)
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
    ap.add_argument('--limit', type=int, default=0, help='只跑前 N 张（G0 成本实测用）')
    ap.add_argument('--out-tag', default='', help='给输出目录加后缀，试跑用，避免污染正式产物')
    args = ap.parse_args()
    for n in args.arms:
        assert n in ARM_BY_NAME, f'未知臂 {n}'

    global OUT, PER_SLIDE, EMB
    if args.out_tag:
        OUT = f'{OUT}_{args.out_tag}'
        PER_SLIDE, EMB = f'{OUT}/per_slide', f'{OUT}/embeds'
    for n in args.arms:
        os.makedirs(f'{PER_SLIDE}/{n}', exist_ok=True)
        os.makedirs(f'{EMB}/{n}', exist_ok=True)

    torch.set_num_threads(N_THREADS)
    t_all = time.time()
    model = CLIPModel.from_pretrained(MODEL_DIR).eval()
    proc = CLIPProcessor.from_pretrained(MODEL_DIR)
    with torch.no_grad():   # 文本编码只做一次（两套共 15 条），与图像无关
        teA = model.get_text_features(**proc(text=[TEMPLATE_A.format(l) for l in LABELS_A],
                                            return_tensors='pt', padding=True))
        teA = teA / teA.norm(dim=-1, keepdim=True)
        teB = model.get_text_features(**proc(text=[TEMPLATE_B.format(l) for l in LABELS_B],
                                            return_tensors='pt', padding=True))
        teB = teB / teB.norm(dim=-1, keepdim=True)

    slides = []
    for d in sorted(os.listdir(DATA)):
        p = parse_dir(d)
        if p and os.path.isdir(f'{DATA}/{d}'):
            slides.append((*p, d))
    if args.limit:
        slides = slides[:args.limit]

    done = skipped = failed = 0
    tot = 0
    n_slides_all = len(slides)
    per_slide_log = []
    for i, (gsm, pat, stage, gdir) in enumerate(slides, 1):
        todo = [ARM_BY_NAME[n] for n in args.arms
                if not (csv_ok(f'{PER_SLIDE}/{n}/{gdir}.csv.gz')
                        and os.path.exists(f'{EMB}/{n}/{gdir}.npy'))]
        if not todo:
            skipped += 1
            continue
        lg = {}
        t1 = time.time()
        try:
            res = run_slide(model, proc, teA, teB, gdir, todo, lg)
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
                        f'lesB={df.lesion_B.mean():.3f} chroma={df.chroma_std.median():.1f}')
        done += 1
        lg.update(gsm_dir=gdir, stage=stage, sec=round(time.time() - t1, 1))
        per_slide_log.append(lg)
        el = time.time() - t_all
        eta = el / (done + skipped) * (n_slides_all - done - skipped) / 60 if (done + skipped) else 0
        print(f'  [{i:2d}/{n_slides_all}] {gdir:22s} {stage:5s} ' + ' | '.join(line)
              + f'  {lg["sec"]}s/张  ETA {eta:.0f}min', flush=True)

    el_min = (time.time() - t_all) / 60
    print(f'\n[done] 新跑 {done} / 跳过 {skipped} / 失败 {failed}；'
          f'{tot:,} spot；{el_min:.1f} 分钟')

    summary = {
        'script': '07_he_pathology/19_resolution_ladder.py',
        'status': 'RAW MEASUREMENTS — 只产出测量值，不施加任何类别判定、不判恶性、不出 WHO',
        'prereg': 'results/07_he_pathology/resolution_sweep/PREREG_v2.md',
        'arms_run': args.arms,
        'model': 'vinid/plip (CLIP ViT-B/32)',
        'note_arch_ceiling': 'ViT-B/32 ⇒ 224px/32px = 7×7 token ⇒ 1 token = 32 输入像素；'
                             'L448 一个 token 363 µm、L224 181 µm。'
                             '平坦结果的措辞只能是「在 5.670 µm/px 以上判别力不随细节变化」。',
        'note_p1_background': 'P1 的背景算在整张切片上（登记的实现澄清，见脚本 docstring）',
        'frozen': dict(RNG_SEED=RNG_SEED, N_PER_SLIDE=N_PER_SLIDE, CROPS=CROPS, INPUT=INPUT,
                       BATCH=BATCH, N_THREADS=N_THREADS, MIN_TISSUE_FRAC=MIN_TISSUE_FRAC,
                       WHITE_THRESH=WHITE_THRESH, FLATTEN_SIGMA=FLATTEN_SIGMA,
                       MIN_CHROMA_STD=MIN_CHROMA_STD,
                       template_A=TEMPLATE_A, labels_A=LABELS_A,
                       template_B=TEMPLATE_B, labels_B=LABELS_B),
        'anchor_arm': ANCHOR,
        'spot_filter_rule': '一律用 L448 裁框（hires 448px）量 tissue_frac，六臂同一 spot 集；'
                            'tissue_frac_224 只登记不过滤',
        'per_slide_log': per_slide_log,
        'n_slide': n_slides_all, 'n_done': done, 'n_skipped': skipped, 'n_failed': failed,
        'elapsed_min': round(el_min, 1),
    }
    with open(f'{OUT}/ladder_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print(f'[out] {OUT}/ladder_summary.json')


if __name__ == '__main__':
    main()
