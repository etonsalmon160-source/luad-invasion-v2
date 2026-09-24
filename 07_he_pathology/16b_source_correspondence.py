#!/usr/bin/env python3
"""源图对应性审计（16b）—— aligned_tissue_image.jpg 与 tissue_hires_image.png 是什么关系？

要回答三件事：
  Q1  两张图是不是**同一块组织**？
  Q2  是不是**镜像/旋转**关系（用户提出：看着像轴对称）？
  Q3  如果是同一块，**视野差几倍**（⇒ aligned 的实际 µm/px 是多少）？

为什么必须重做：`16_source_audit.py` 的翻折比较是
「整幅直接对位」（把 aligned 按 1/3 缩到 hires 格子上再比），
它同时假设了 ①比例恰为 3.0 ②两者共原点。**两个假设都没验证过。**
实测：hires 是**紧贴组织的裁切**（组织占满画面），
aligned 是**宽视野**（组织只占中间一块，四周白边 + 灰色拼块），
⇒ 1/3 缩放根本没把两张图放到同一个物理尺度上，那个比较**无效**。

本脚本改用**不受比例影响**的测法：
  A. 取两边**组织轮廓**（Otsu），各自裁到外接框、归一化到同一尺寸，
     再比 8 种翻折/旋转（该测法把"视野差几倍"吸收掉，只留"形状对不对"）
  B. 在**边缘图**上做 比例+平移 搜索（边缘图不受明暗/灰阶影响，比亮度稳），
     定出**视野比** ⇒ 反推 aligned 的真实 µm/px

只做测量，不判恶性、不产标签。
"""
import os
import sys
import csv
import glob
import json
import time

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy.signal import fftconvolve

Image.MAX_IMAGE_PIXELS = None

ROOT = '/home/eto/luad_v2'
SP = os.path.join(ROOT, 'data/visium_spatial')
OUT = os.path.join(ROOT, 'results/07_he_pathology/resolution_sweep')

HIRES_UM_PER_PX = 5.670     # 由 fullres 0.25050 / tissue_hires_scalef 0.044178393
WORK = 1600                 # 轮廓提取的工作网格
NORM = 512                  # 归一化后的比较尺寸
SEARCH_R = 40               # 归一化图上的平移搜索半径
CAN = 1500                  # 边缘图比较用的画布
HIRES_GRID = 700            # 边缘图比较时 hires 载入的边长（画布够装 q 倍）


def spatial_dir(gsm):
    flat = os.path.join(SP, gsm, 'spatial')
    if os.path.isdir(flat):
        return flat
    for d in sorted(glob.glob(os.path.join(SP, gsm, '*'))):
        if os.path.isdir(os.path.join(d, 'spatial')):
            return os.path.join(d, 'spatial')
    return None


def otsu_thresh(g):
    hist, _ = np.histogram(g, bins=256, range=(0.0, 255.0))
    hist = hist.astype(np.float64)
    tot = hist.sum()
    if tot <= 0:
        return 128.0
    w = np.cumsum(hist)
    mu = np.cumsum(hist * np.arange(256, dtype=np.float64))
    denom = w * (tot - w)
    denom[denom <= 0] = np.nan
    return float(np.nanargmax((mu[-1] * w - mu * tot) ** 2 / denom))


def load_work(path, work=WORK):
    im = Image.open(path)
    w, h = im.size
    sc = work / max(w, h)
    a = np.asarray(im.convert('L').resize(
        (max(1, int(round(w * sc))), max(1, int(round(h * sc)))), Image.BOX),
        dtype=np.float32)
    return a


def silhouette(a):
    return (a < otsu_thresh(a)).astype(np.float32)


def norm_mask(m, N=NORM):
    ys, xs = np.where(m > 0.5)
    if len(ys) < 10:
        return None, None, None
    sub = m[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    asp = sub.shape[1] / sub.shape[0]
    r = np.asarray(Image.fromarray((sub * 255).astype(np.uint8)).resize((N, N), Image.BILINEAR),
                   dtype=np.float32) / 255.0
    return (r > 0.5).astype(np.float32), asp, sub.shape


TRANSFORMS = [
    ('as-is',         lambda a: a),
    ('flip-LR',       lambda a: a[:, ::-1]),
    ('flip-UD',       lambda a: a[::-1, :]),
    ('rot180',        lambda a: a[::-1, ::-1]),
    ('rot90',         lambda a: np.rot90(a, 1)),
    ('rot270',        lambda a: np.rot90(a, 3)),
    ('transpose',     lambda a: a.T),
    ('antitranspose', lambda a: np.rot90(a.T, 2)),
]


def best_shift(A, B, R):
    An = A - A.mean()
    Bn = B - B.mean()
    c = fftconvolve(An, Bn[::-1, ::-1], mode='same')
    ny, nx = A.shape
    cy, cx = ny // 2, nx // 2
    yy, xx = np.ogrid[:ny, :nx]
    c = np.where((np.abs(yy - cy) <= R) & (np.abs(xx - cx) <= R), c, -np.inf)
    i = np.unravel_index(np.argmax(c), c.shape)
    den = np.sqrt(float((An * An).sum()) * float((Bn * Bn).sum()))
    return (float(c[i] / den) if den > 0 else 0.0), int(i[0] - cy), int(i[1] - cx)


def pad_center(a, can):
    out = np.zeros((can, can), dtype=np.float32)
    h, w = a.shape
    if h > can or w > can:
        y0 = max(0, (h - can) // 2)
        x0 = max(0, (w - can) // 2)
        a = a[y0:y0 + can, x0:x0 + can]
        h, w = a.shape
    y = (can - h) // 2
    x = (can - w) // 2
    out[y:y + h, x:x + w] = a
    return out


def edge_map(a, sig=2.0):
    s = gaussian_filter(a, sig)
    gx = np.gradient(s, axis=1)
    gy = np.gradient(s, axis=0)
    return np.hypot(gx, gy)


def field_ratio_scan(pa, ph, ratios):
    """在**边缘图**上找「视野比 q」：aligned 视野 = q × hires 视野。
    同一物理尺度要求 aligned 的像素数 = q × hires 的像素数
    （q>1 ⇒ aligned 视野更宽、需要被**放大**到 q 倍才能和 hires 同尺度相比）。
    返回 [(q, ncc, dy, dx, scaled_shape)]。"""
    ah = load_work(ph, HIRES_GRID)
    ih = Image.open(ph)
    wh, hh = ih.size
    g = ah.shape[1] / wh                     # hires 载入后的缩放
    ea = edge_map(ah)
    EAn = pad_center(ea, CAN)
    EAn = (EAn - EAn.mean()) / (EAn.std() + 1e-6)
    ia = Image.open(pa)
    out = []
    for q in ratios:
        tw = int(round(wh * g * q))
        th = int(round(hh * g * q))
        if tw < 64 or th < 64 or tw > CAN or th > CAN:
            continue
        aa = np.asarray(ia.convert('L').resize((tw, th), Image.BOX), dtype=np.float32)
        EA2 = pad_center(edge_map(aa), CAN)
        EA2 = (EA2 - EA2.mean()) / (EA2.std() + 1e-6)
        ncc, dy, dx = best_shift(EAn, EA2, R=CAN // 2 - 20)
        out.append((float(q), ncc, int(dy), int(dx), (th, tw)))
    return out, dict(px_hires=wh, px_aligned=ia.size[0], px_hires_h=hh)


def main():
    limit = None
    if '--limit' in sys.argv:
        limit = int(sys.argv[sys.argv.index('--limit') + 1])
    do_scan = '--no-scan' not in sys.argv

    gsms = sorted(d for d in os.listdir(SP) if d.startswith('GSM'))
    if limit:
        gsms = gsms[:limit]

    t = np.random.default_rng(0).normal(size=(300, 300)).astype(np.float32)
    tb = np.roll(t, (7, -11), axis=(0, 1))
    ncc_t, dyt, dxt = best_shift(t, tb, R=50)
    ok = (abs(dyt + 7) <= 1) and (abs(dxt - 11) <= 1) and ncc_t > 0.9
    print(f'[自检] 已知位移 (dy=+7,dx=-11) ⇒ 找到 (dy={dyt},dx={dxt}) ncc={ncc_t:.3f} '
          f'{"OK" if ok else "🔴 符号约定错了，结果不可信"}')
    if not ok:
        return 1

    rows = []
    example_scan = None
    t0 = time.time()
    order = [n for n, _ in TRANSFORMS]
    for k, gsm in enumerate(gsms):
        sd = spatial_dir(gsm)
        if sd is None:
            continue
        ph = os.path.join(sd, 'tissue_hires_image.png')
        pa = os.path.join(sd, 'aligned_tissue_image.jpg')
        if not (os.path.exists(ph) and os.path.exists(pa)):
            print(f'[{gsm}] 缺图，跳过')
            continue
        mh = silhouette(load_work(ph))
        ma = silhouette(load_work(pa))
        H, asp_h, sh_h = norm_mask(mh)
        A, asp_a, sh_a = norm_mask(ma)
        if H is None or A is None:
            print(f'[{gsm}] 轮廓取不到，跳过')
            continue
        vals = {}
        for name, fn in TRANSFORMS:
            v, _, _ = best_shift(H, fn(A), R=SEARCH_R)
            vals[name] = round(v, 4)
        best = max(vals, key=lambda k2: vals[k2])
        row = dict(gsm=gsm, asp_hires=round(asp_h, 4), asp_aligned=round(asp_a, 4),
                   frac_hires=round(float(mh.mean()), 4), frac_aligned=round(float(ma.mean()), 4),
                   best_orient=best, ncc_best=vals[best],
                   ncc_as_is=vals['as-is'], ncc_flipLR=vals['flip-LR'],
                   ncc_flipUD=vals['flip-UD'],
                   **{f'o_{n}': vals[n] for n in order})
        if do_scan:
            sc, meta = field_ratio_scan(pa, ph, np.arange(0.60, 2.05, 0.02))
            if sc:
                bq = max(sc, key=lambda t2: t2[1])
                # aligned 视野 = q × hires 视野；两边原生像素数之比 px_a/px_h
                # ⇒ aligned 的 µm/px = hires 的 µm/px × q × px_h/px_a
                um_a = HIRES_UM_PER_PX * bq[0] * meta['px_hires'] / meta['px_aligned']
                row.update(q_best=round(bq[0], 3), q_ncc=round(bq[1], 4),
                           um_per_px_aligned=round(um_a, 4),
                           field_mm_aligned=round(HIRES_UM_PER_PX * meta['px_hires'] * bq[0] / 1000.0, 3),
                           q_scaled=f'{bq[4][0]}x{bq[4][1]}')
                if example_scan is None:
                    example_scan = (gsm, meta, sc, bq)
        rows.append(row)
        msg = (f'[{k+1}/{len(gsms)}] {gsm}  宽高比 h {asp_h:.3f}/a {asp_a:.3f}  '
               f'最佳摆法 {best} {vals[best]:.3f}  (as-is {vals["as-is"]:.3f} / '
               f'flipLR {vals["flip-LR"]:.3f} / flipUD {vals["flip-UD"]:.3f})')
        if 'q_best' in row:
            msg += f'  视野比 {row["q_best"]:.2f}'
        print(msg)

    print(f'\n{len(rows)} 张切片，{time.time()-t0:.0f} s')
    if not rows:
        return 1

    print('\n[Q2 朝向] 8 种摆法的跨切片中位数 NCC（轮廓已归一化到同尺寸 ⇒ 比例被吸收）')
    print(f'{"摆法":<16}{"NCC中位数":>10}')
    med = {}
    for n in order:
        med[n] = float(np.median([r[f'o_{n}'] for r in rows]))
    for n in sorted(order, key=lambda k2: -med[k2]):
        flag = '  ← 最好' if n == max(order, key=lambda k2: med[k2]) else ''
        print(f'{n:<16}{med[n]:>10.3f}{flag}')
    n_same = sum(1 for r in rows if r['best_orient'] == 'as-is')
    print(f'\n⇒ {n_same}/{len(rows)} 张切片的**最佳摆法是 as-is（原样）**')
    mirror = np.median([max(r['o_flip-LR'], r['o_flip-UD']) for r in rows])
    print(f'⇒ 镜像类（flip-LR / flip-UD 取大者）的跨切片中位数 = {mirror:.3f}')
    print(f'⇒ as-is {med["as-is"]:.3f} vs 镜像 {mirror:.3f}')

    if do_scan and any('q_best' in r for r in rows):
        qs = np.array([r['q_best'] for r in rows if 'q_best' in r])
        ums = np.array([r['um_per_px_aligned'] for r in rows if 'um_per_px_aligned' in r])
        flds = np.array([r['field_mm_aligned'] for r in rows if 'field_mm_aligned' in r])
        print(f'\n[Q3 视野比] aligned 视野 / hires 视野，跨切片中位数 = {np.median(qs):.3f} '
              f'(范围 {qs.min():.2f}–{qs.max():.2f}, n={len(qs)}, 不一致的切片数 '
              f'{int((qs != qs[0]).sum())})')
        print(f'⇒ aligned 视野 ≈ {np.median(flds):.2f} mm（hires {HIRES_UM_PER_PX*1995/1000:.2f} mm）')
        print(f'⇒ aligned 的实际分辨率 ≈ {np.median(ums):.3f} µm/px '
              f'⇒ 比 hires（{HIRES_UM_PER_PX} µm/px）细 {HIRES_UM_PER_PX/np.median(ums):.2f} 倍')
        print(f'   ⚠️ 旧记录写的 1.890 µm/px 来自"两边同视野、像素差 3.0 倍"的**未验证假设**；'
              f'实测视野差 {np.median(qs):.2f} 倍 ⇒ 旧值**作废**')

    # ---------------- 图 ----------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f'[warn] 无 matplotlib（{e}），跳过出图')
        return 0

    gsm = rows[0]['gsm']
    sd = spatial_dir(gsm)
    ph = os.path.join(sd, 'tissue_hires_image.png')
    pa = os.path.join(sd, 'aligned_tissue_image.jpg')
    mh = silhouette(load_work(ph))
    ma = silhouette(load_work(pa))
    H, _, _ = norm_mask(mh)
    A, _, _ = norm_mask(ma)

    def ov(a, b):
        o = np.zeros(a.shape + (3,), dtype=np.float32)
        o[..., 1] = a          # green = hires
        o[..., 0] = b          # magenta = aligned
        o[..., 2] = b
        return o

    fig, ax = plt.subplots(2, 4, figsize=(20, 10.4))
    ax[0, 0].imshow(H, cmap='gray'); ax[0, 0].set_title(f'hires silhouette (box-normalised)\naspect {rows[0]["asp_hires"]:.3f}')
    ax[0, 1].imshow(A, cmap='gray'); ax[0, 1].set_title(f'aligned silhouette (box-normalised)\naspect {rows[0]["asp_aligned"]:.3f}')
    ax[0, 2].imshow(ov(H, A)); ax[0, 2].set_title(f'as-is overlay  NCC {rows[0]["ncc_as_is"]:.3f}\ngreen=hires  magenta=aligned')
    ax[0, 3].imshow(ov(H, A[:, ::-1])); ax[0, 3].set_title(f'flip-LR overlay  NCC {rows[0]["ncc_flipLR"]:.3f}\n(the "mirror" hypothesis)')

    # 下排：按实测比例 q 把 aligned 摆到 hires 的画幅上（真·同尺度对照）
    if example_scan is not None:
        gsm2, meta, sc, bq = example_scan
        sd2 = spatial_dir(gsm2)
        ph2 = os.path.join(sd2, 'tissue_hires_image.png')
        pa2 = os.path.join(sd2, 'aligned_tissue_image.jpg')
        G = 900
        ih2 = Image.open(ph2)
        wh2, hh2 = ih2.size
        hg = np.asarray(ih2.convert('L').resize((G, max(1, int(round(G * hh2 / wh2)))), Image.BOX), dtype=np.float32)
        q = bq[0]
        tw = int(round(wh2 * (G / wh2) * q)); th = int(round(hh2 * (G / wh2) * q))
        ag = np.asarray(Image.open(pa2).convert('L').resize((tw, th), Image.BOX), dtype=np.float32)
        # 把 aligned 缩回 hires 的格子，用扫描出的平移对齐
        ag_g = np.asarray(Image.fromarray(ag.astype(np.uint8)).resize((G, hg.shape[0]), Image.BOX), dtype=np.float32)
        e1 = edge_map(hg); e2 = edge_map(ag_g)
        ncc2, dy2, dx2 = best_shift((e1 - e1.mean()) / (e1.std() + 1e-6),
                                    (e2 - e2.mean()) / (e2.std() + 1e-6), R=max(hg.shape) // 3)
        ag_g = np.roll(ag_g, (dy2, dx2), axis=(0, 1))
        s = lambda x: np.clip((x - np.percentile(x, 1)) / (np.percentile(x, 99) - np.percentile(x, 1) + 1e-6), 0, 1)
        ax[1, 0].imshow(s(hg), cmap='gray'); ax[1, 0].set_title(f'hires (as delivered)  {hg.shape}')
        ax[1, 1].imshow(s(ag_g), cmap='gray'); ax[1, 1].set_title(
            f'aligned, scaled by measured q={q:.2f} + aligned  edge NCC {ncc2:.3f}')
        ov2 = np.zeros(hg.shape + (3,), dtype=np.float32)
        ov2[..., 1] = s(hg); ov2[..., 0] = s(ag_g); ov2[..., 2] = s(ag_g)
        ax[1, 2].imshow(ov2); ax[1, 2].set_title('edge overlay at the measured scale\ngrey = agree, colour = disagree')
        ax[1, 3].plot([t[0] for t in sc], [t[1] for t in sc], marker='o', ms=3)
        ax[1, 3].axvline(q, color='r', ls='--', label=f'peak q={q:.2f}')
        ax[1, 3].set_title('edge-NCC vs candidate field ratio q'); ax[1, 3].legend(fontsize=8)
        ax[1, 3].set_xlabel('q = aligned field / hires field')
    else:
        for j in range(4):
            ax[1, j].axis('off')
    for a in (ax[0, 0], ax[0, 1], ax[0, 2], ax[0, 3]):
        a.set_xticks([]); a.set_yticks([])
    for j in range(3):
        ax[1, j].set_xticks([]); ax[1, j].set_yticks([])
    fig.suptitle(f'Top: is one a mirror of the other? (silhouettes box-normalised, scale difference absorbed)  |  '
                 f'Bottom: same-scale comparison using the measured field ratio q', fontsize=12)
    fig.tight_layout()
    fp = os.path.join(OUT, 'source_correspondence.png')
    fig.savefig(fp, dpi=100)
    plt.close(fig)
    print(f'[out] {fp}')

    fp2 = os.path.join(OUT, 'source_correspondence.csv')
    with open(fp2, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    fp3 = os.path.join(OUT, 'source_correspondence_summary.json')
    with open(fp3, 'w') as f:
        json.dump(dict(n_slides=len(rows),
                       median_ncc_by_orientation=med,
                       n_best_is_as_is=n_same,
                       median_mirror_ncc=float(mirror),
                       median_field_ratio=float(np.median(qs)) if do_scan and len(qs) else None,
                       median_field_mm_aligned=float(np.median(flds)) if do_scan and len(flds) else None,
                       median_um_per_px_aligned=float(np.median(ums)) if do_scan and len(ums) else None,
                       hires_um_per_px=HIRES_UM_PER_PX), f, indent=2, ensure_ascii=False)
    print(f'[out] {fp2}')
    print(f'[out] {fp3}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
