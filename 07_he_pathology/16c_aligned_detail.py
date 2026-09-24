#!/usr/bin/env python3
"""aligned_tissue 到底有没有**真细节**（16c）—— 回答「换这张图会不会改进」。

背景：`16_source_audit.py` 里的 `highfreq_test()` 问过同一个问题，但它的对照项
拿的是「hires 放大 **3 倍**」，而实测视野比是 **1.30**（⇒ aligned 只比 hires 细
**2.31 倍**，见 `16b_source_correspondence.py`）。**对照选错 ⇒ 那个测量作废。**
本脚本用**正确比例**重做。

三个问题，逐层收窄：
  Q1  在正确比例下把两图对齐后，**同一物理位置的像素值是否一致**？
      （rho 高 ⇒ aligned 是同一组织的真实渲染，不是插值出来的无关内容）
  Q2  aligned 的高频功率，是否**高于**「hires 按真实比例放大」这个
      **不含新信息**的对照？（高 ⇒ 可能有真细节）
  Q3  若有超出，那部分功率**长什么样**？（看图判断是真组织纹理还是压缩伪影）

只做测量。不判恶性、不产标签、不改任何既有结果。
"""
import os
import sys
import glob
import json
import csv

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, uniform_filter, shift as ndshift
from scipy.signal import fftconvolve

Image.MAX_IMAGE_PIXELS = None

ROOT = '/home/eto/luad_v2'
SP = os.path.join(ROOT, 'data/visium_spatial')
OUT = os.path.join(ROOT, 'results/07_he_pathology/resolution_sweep')

HIRES_UM_PER_PX = 5.670            # fullres 0.25050 / tissue_hires_scalef 0.044178393
Q_FIELD = 1.30                     # 实测视野比（16b，4/4 张一致）
HI_PX = 1995                       # hires 原生宽
AL_PX = 5985                       # aligned 原生宽
PATCH = 512                        # hires 格子上的取样块边长
N_PATCH = 4


def spatial_dir(gsm):
    flat = os.path.join(SP, gsm, 'spatial')
    if os.path.isdir(flat):
        return flat
    for d in sorted(glob.glob(os.path.join(SP, gsm, '*'))):
        if os.path.isdir(os.path.join(d, 'spatial')):
            return os.path.join(d, 'spatial')
    return None


def place(g, can):
    """把 g 居中放进 can×can 画布，返回 (画布, y偏移, x偏移)。"""
    out = np.zeros((can, can), dtype=np.float32)
    h, w = g.shape
    y, x = (can - h) // 2, (can - w) // 2
    out[y:y + h, x:x + w] = g
    return out, y, x


def best_shift(A, B, R):
    An, Bn = A - A.mean(), B - B.mean()
    c = fftconvolve(An, Bn[::-1, ::-1], mode='same')
    ny, nx = A.shape
    cy, cx = ny // 2, nx // 2
    yy, xx = np.ogrid[:ny, :nx]
    c = np.where((np.abs(yy - cy) <= R) & (np.abs(xx - cx) <= R), c, -np.inf)
    i = np.unravel_index(np.argmax(c), c.shape)
    den = np.sqrt(float((An * An).sum()) * float((Bn * Bn).sum()))
    return (float(c[i] / den) if den > 0 else 0.0), int(i[0] - cy), int(i[1] - cx)


def edge_map(a, sig=2.0):
    s = gaussian_filter(a, sig)
    return np.hypot(np.gradient(s, axis=1), np.gradient(s, axis=0))


def fft_upsample(g, ny, nx):
    """理想插值（FFT 补零）放大：**保留下方全部功率，Nyquist 以上功率恰为 0**。
    这是「不含任何新信息」的正确对照 —— 双线性放大会额外低通掉中频，不公平。"""
    F = np.fft.fftshift(np.fft.fft2(g))
    h, w = g.shape
    out = np.zeros((ny, nx), dtype=complex)
    oy, ox = (ny - h) // 2, (nx - w) // 2
    out[oy:oy + h, ox:ox + w] = F
    return np.real(np.fft.ifft2(np.fft.ifftshift(out))) * (ny * nx) / (h * w)


def band_rho(H, A, um_per_px, lo, hi):
    """H 与 A（同一 hires 格子）在 [lo,hi] cycles/µm 环带内的相关系数。"""
    ny, nx = H.shape
    Hw = (H - H.mean()) * np.hanning(ny)[:, None] * np.hanning(nx)[None, :]
    Aw = (A - A.mean()) * np.hanning(ny)[:, None] * np.hanning(nx)[None, :]
    FH, FA = np.fft.fftshift(np.fft.fft2(Hw)), np.fft.fftshift(np.fft.fft2(Aw))
    yy, xx = np.mgrid[0:ny, 0:nx]
    cy, cx = ny // 2, nx // 2
    k = np.sqrt(((yy - cy) / (ny * um_per_px)) ** 2 + ((xx - cx) / (nx * um_per_px)) ** 2)
    m = (k >= lo) & (k < hi)
    if m.sum() < 32:
        return float('nan')
    a, b = FH[m].real, FA[m].real
    if a.std() < 1e-9 or b.std() < 1e-9:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def jpeg_block_score(g):
    """8×8 分块指纹：横向相邻列差的均值，按列下标 mod 8 分组。
    满格比（max/median）越高 ⇒ 越像 JPEG 块效应（⇒ 高频里有压缩伪影的份）。"""
    d = np.abs(np.diff(g, axis=1)).mean(axis=0)
    n = d.size
    s = [float(d[i::8].mean()) for i in range(8)]
    med = float(np.median(s))
    return (max(s) / med if med > 0 else float('nan')), int(np.argmax(s)), [round(v, 3) for v in s]


def radial_psd(g, um_per_px, nb=48):
    ny, nx = g.shape
    gw = g * np.hanning(ny)[:, None] * np.hanning(nx)[None, :]
    F = np.fft.fftshift(np.fft.fft2(gw))
    P = (np.abs(F) ** 2) / (ny * nx)
    yy, xx = np.mgrid[0:ny, 0:nx]
    cy, cx = ny // 2, nx // 2
    k = np.sqrt(((yy - cy) / (ny * um_per_px)) ** 2 + ((xx - cx) / (nx * um_per_px)) ** 2)
    edges = np.linspace(0, k.max(), nb + 1)
    idx = np.digitize(k.ravel(), edges) - 1
    p = P.ravel()
    ok = (idx >= 0) & (idx < nb)
    num = np.bincount(idx[ok], weights=p[ok], minlength=nb)
    cnt = np.bincount(idx[ok], minlength=nb).astype(float)
    cnt[cnt == 0] = np.nan
    return 0.5 * (edges[1:] + edges[:-1]), num / cnt


def main():
    limit = 3
    if '--limit' in sys.argv:
        limit = int(sys.argv[sys.argv.index('--limit') + 1])

    um_a = HIRES_UM_PER_PX * Q_FIELD * HI_PX / AL_PX      # aligned 原生 µm/px ≈ 2.457
    f = um_a / HIRES_UM_PER_PX                            # aligned → hires 尺度 ≈ 0.4333
    up = HIRES_UM_PER_PX / um_a                           # hires 需放大几倍才同尺度 ≈ 2.308
    nyq_h = 1.0 / (2 * HIRES_UM_PER_PX)
    print(f'[比例] 视野比 q={Q_FIELD:.2f} ⇒ aligned 原生 ≈ {um_a:.3f} µm/px；'
          f'hires {HIRES_UM_PER_PX} µm/px ⇒ hires 需放大 {up:.3f} 倍才同尺度')
    print(f'[比例] hires 真实 Nyquist = {nyq_h:.4f} cycles/µm（放大后其上的功率全是插值）')

    gsms = sorted(d for d in os.listdir(SP) if d.startswith('GSM'))[:limit]
    rows = []
    example = None
    for gsm in gsms:
        sd = spatial_dir(gsm)
        if sd is None:
            continue
        ph = os.path.join(sd, 'tissue_hires_image.png')
        pa = os.path.join(sd, 'aligned_tissue_image.jpg')
        if not (os.path.exists(ph) and os.path.exists(pa)):
            continue
        H = np.asarray(Image.open(ph).convert('L'), dtype=np.float32)
        ia = Image.open(pa)
        a_nat = np.asarray(ia.convert('L'), dtype=np.float32)
        A = np.asarray(ia.convert('L').resize(
            (int(round(ia.size[0] * f)), int(round(ia.size[1] * f))), Image.BOX), dtype=np.float32)

        can = 1 << int(np.ceil(np.log2(max(H.shape[0], A.shape[0]) + 60)))
        eH, ohy, ohx = place(edge_map(H), can)
        eA, oay, oax = place(edge_map(A), can)
        ncc, dy, dx = best_shift((eH - eH.mean()) / (eH.std() + 1e-6),
                                 (eA - eA.mean()) / (eA.std() + 1e-6), R=can // 2 - 20)
        # A 平移 (dy,dx) 后与 H 对齐 ⇒ A 的画布原点变为 (oay+dy, oax+dx)
        aoy, aox = oay + dy, oax + dx
        # H 画布上的重叠区
        y0, y1 = max(ohy, aoy), min(ohy + H.shape[0], aoy + A.shape[0])
        x0, x1 = max(ohx, aox), min(ohx + H.shape[1], aox + A.shape[1])
        ov = (y1 - y0, x1 - x0)
        # 重叠区里挑组织块（H 上）
        sub = H[y0 - ohy:y1 - ohy, x0 - ohx:x1 - ohx]
        m = uniform_filter((sub < 230).astype(np.float32), 200)
        picks, tmp = [], m.copy()
        for _ in range(N_PATCH):
            iy, ix = np.unravel_index(np.argmax(tmp), tmp.shape)
            picks.append((iy, ix))
            tmp[max(0, iy - 300):iy + 300, max(0, ix - 300):ix + 300] = -1

        rhos, ex_here = [], None
        for iy, ix in picks:
            py = min(max(0, iy - PATCH // 2), sub.shape[0] - PATCH)
            px = min(max(0, ix - PATCH // 2), sub.shape[1] - PATCH)
            if py < 0 or px < 0:
                continue
            hp = sub[py:py + PATCH, px:px + PATCH]
            # sub 坐标 → H 原生坐标 → 画布 → A 索引
            cy_, cx_ = y0 - ohy + py, x0 - ohx + px
            ay, ax = ohy + cy_ - aoy, ohx + cx_ - aox
            ap = A[ay:ay + PATCH, ax:ax + PATCH]
            if ap.shape != (PATCH, PATCH):
                continue
            if hp.std() < 2 or ap.std() < 2:
                continue
            r = float(np.corrcoef(hp.ravel(), ap.ravel())[0, 1])
            rhos.append(r)
            if ex_here is None:
                ex_here = dict(py=py, px=px, r=r)
        rho = float(np.median(rhos)) if rhos else float('nan')
        rows.append(dict(gsm=gsm, edge_ncc=round(ncc, 4), dy=dy, dx=dx,
                         overlap=f'{ov[0]}x{ov[1]}', rho_pixel=round(rho, 4), n_patch=len(rhos)))
        print(f'[{gsm}] 边缘 NCC {ncc:.3f}  重叠 {ov[0]}×{ov[1]}  '
              f'同位置 rho（{len(rhos)} 块中位）= {rho:+.3f}')
        if example is None and ex_here is not None:
            # 保留 H 的这块 + aligned 原生同物理区（不缩放，供谱分析）
            cy_, cx_ = y0 - ohy + ex_here['py'], x0 - ohx + ex_here['px']
            example = dict(gsm=gsm, H=H, hc=(cy_, cx_), a_nat=a_nat,
                           map=(ohy, ohx, aoy, aox, f), r=ex_here['r'])

    if example is None:
        print('[stop] 没取到任何有效块，谱分析跳过')
        return 1

    # ——— 谱：同一物理视野，三种取法 ———
    gsm, H, a_nat, r = example['gsm'], example['H'], example['a_nat'], example['r']
    ohy, ohx, aoy, aox, fm = example['map']
    rh, rw = H.shape
    # 同一物理视野：hires 上 hs px ↔ aligned 原生 hs×up px（up = 5.670/2.457）
    hs = PATCH
    while hs >= 192:
        hy0 = min(max(0, example['hc'][0] - hs // 2), rh - hs)
        hx0 = min(max(0, example['hc'][1] - hs // 2), rw - hs)
        hcy, hcx = hy0 + hs // 2, hx0 + hs // 2
        acy = int(round((ohy + hcy - aoy) / fm))
        acx = int(round((ohx + hcx - aox) / fm))
        asz = int(round(hs * up))
        ay0, ax0 = acy - asz // 2, acx - asz // 2
        if ay0 >= 0 and ax0 >= 0 and ay0 + asz <= a_nat.shape[0] and ax0 + asz <= a_nat.shape[1]:
            break
        hs //= 2
    else:
        print('[stop] aligned 块越界（缩小到 192 px 仍不行）')
        return 1
    hy0, hx0 = hcy - hs // 2, hcx - hs // 2
    hp = H[hy0:hy0 + hs, hx0:hx0 + hs]
    ap = a_nat[ay0:ay0 + asz, ax0:ax0 + asz]
    print(f'[谱] 取块 hires ({hy0},{hx0}) {hs}px ↔ aligned 原生 ({ay0},{ax0}) {asz}px')
    hup = fft_upsample(hp, asz, asz)     # 理想插值：Nyquist 以上功率恰为 0（正确对照）
    ap_h = np.asarray(Image.fromarray(ap.astype(np.uint8)).resize((hs, hs), Image.BOX), dtype=np.float32)

    kA, pA = radial_psd(ap, um_a)
    kH, pH = radial_psd(hp, HIRES_UM_PER_PX)
    kU, pU = radial_psd(hup, um_a)

    band = (kU > nyq_h * 1.05) & (kU < kU.max() * 0.90)
    ratio = float(np.exp(np.mean(np.log(np.maximum(pU[band], 1e-12))) -
                         np.mean(np.log(np.maximum(pA[band], 1e-12)))))
    lo = (kU > 0.02) & (kU < nyq_h * 0.95)
    ratio_lo = float(np.exp(np.mean(np.log(np.maximum(pU[lo], 1e-12))) -
                            np.mean(np.log(np.maximum(pA[lo], 1e-12)))))

    bands = [(0.005, 0.02), (0.02, 0.045), (0.045, 0.070), (0.070, 0.088)]
    brho = [(a, b, band_rho(hp, ap_h, HIRES_UM_PER_PX, a, b)) for a, b in bands]
    jb, jphase, jprof = jpeg_block_score(ap)

    # 整数像素配准会在细尺度上**人为去相关**（细带周期只有 11–14 px，
    # 差 1–2 px 就够把 rho 打到 0 甚至负）。所以细带结论必须**先做亚像素配准**。
    def bandpass(g, um, lo_, hi_):
        ny, nx = g.shape
        gw = (g - g.mean()) * np.hanning(ny)[:, None] * np.hanning(nx)[None, :]
        F = np.fft.fftshift(np.fft.fft2(gw))
        yy, xx = np.mgrid[0:ny, 0:nx]
        ccy, ccx = ny // 2, nx // 2
        kk = np.sqrt(((yy - ccy) / (ny * um)) ** 2 + ((xx - ccx) / (nx * um)) ** 2)
        F = np.where((kk >= lo_) & (kk < hi_), F, 0)
        return np.real(np.fft.ifft2(np.fft.ifftshift(F)))

    def ncc_at(A, B, dy, dx):
        s = ndshift(B, (dy, dx), order=3, mode='nearest')
        m = (slice(24, -24), slice(24, -24))
        a, b = A[m], s[m]
        a, b = a - a.mean(), b - b.mean()
        d = np.sqrt(float((a * a).sum()) * float((b * b).sum()))
        return float((a * b).sum() / d) if d > 0 else 0.0

    sub = {}
    for lo_, hi_, tag in [(0.045, 0.088, 'fine'), (0.020, 0.045, 'mid'), (0.005, 0.020, 'coarse')]:
        bH = bandpass(hp, HIRES_UM_PER_PX, lo_, hi_)
        bA = bandpass(ap_h, HIRES_UM_PER_PX, lo_, hi_)
        c_int = ncc_at(bH, bA, 0, 0)
        best = (c_int, 0.0, 0.0)
        for dy in np.arange(-4, 4.01, 1.0):
            for dx in np.arange(-4, 4.01, 1.0):
                c = ncc_at(bH, bA, dy, dx)
                if c > best[0]:
                    best = (c, dy, dx)
        _, dy0, dx0 = best
        for dy in np.arange(dy0 - 1, dy0 + 1.001, 0.1):
            for dx in np.arange(dx0 - 1, dx0 + 1.001, 0.1):
                c = ncc_at(bH, bA, dy, dx)
                if c > best[0]:
                    best = (c, dy, dx)
        sub[tag] = dict(band=f'{lo_:.3f}-{hi_:.3f}', rho_integer=round(c_int, 4),
                        rho_subpixel=round(best[0], 4),
                        shift=(round(best[1], 2), round(best[2], 2)))
        print(f'[亚像素] {tag:6s} 带 {lo_:.3f}–{hi_:.3f}：整数对齐 rho = {c_int:+.3f} '
              f'⇒ 亚像素最优 rho = {best[0]:+.3f}（平移 {best[1]:+.1f},{best[2]:+.1f} px）')

    print(f'\n[谱] 例片 {gsm}（同位置 rho {r:+.3f}）；hires {hs}px ↔ aligned 原生 {asz}px，同一物理视野')
    print(f'[谱] hires 真实 Nyquist = {nyq_h:.4f} cycles/µm')
    print(f'[谱] **hires 够得到的频段**（<Nyquist）对照/实测 功率比 = {ratio_lo:.3f}'
          f'  （对照=理想插值 ⇒ ≈1 表示两图在该频段一致；>1 表示 aligned 反而更强）')
    print(f'[谱] **超出 hires 能力的频段**（>Nyquist）对照/实测 功率比 = {ratio:.3f}'
          f'  （对照按构造在此频段功率为 0）')
    print('[分频段一致性] hires 与 aligned（同一 hires 格子）逐带 rho：')
    for a, b, rr in brho:
        print(f'       {a:.3f}–{b:.3f} cyc/µm （{2*HIRES_UM_PER_PX*b:.0f}–{2*HIRES_UM_PER_PX*a:.0f} µm 尺度）  rho = {rr:+.3f}')
    print(f'[JPEG 指纹] aligned 原生 8×8 列差满格比 = {jb:.3f}（峰值在相位 {jphase}）'
          f' ⇒ ' + ('有 8 px 块效应迹象' if jb > 1.15 else '无明显 8 px 块效应'))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f'[warn] 无 matplotlib（{e}），跳过出图')
        return 0
    f2 = lambda x: np.clip((x - np.percentile(x, 1)) / (np.percentile(x, 99) - np.percentile(x, 1) + 1e-6), 0, 1)
    fig, ax = plt.subplots(1, 4, figsize=(21, 5.4))
    ax[0].imshow(f2(hp), cmap='gray', interpolation='nearest')
    ax[0].set_title(f'hires native {hs}px @5.670 µm/px\n(the only colour source we have)')
    ax[1].imshow(f2(np.asarray(Image.fromarray(ap.astype(np.uint8)).resize((hs, hs), Image.BOX),
                               dtype=np.float32)), cmap='gray', interpolation='nearest')
    ax[1].set_title(f'aligned, same field, shown at {hs}px\n(pixel rho vs hires = {r:+.3f})')
    ax[2].imshow(f2(ap), cmap='gray', interpolation='nearest')
    ax[2].set_title(f'aligned NATIVE {asz}px @{um_a:.3f} µm/px\n(greyscale — does the detail exist?)')
    ax[3].semilogy(kH, pH, lw=1.3, color='steelblue', label='hires native (5.670 µm/px)')
    ax[3].semilogy(kU, pU, lw=1.3, color='darkorange', ls='--',
                   label=f'hires, ideal {up:.2f}x upsample (zero info >Nyquist)')
    ax[3].semilogy(kA, pA, lw=1.3, color='crimson', label=f'aligned native ({um_a:.3f} µm/px)')
    ax[3].axvline(nyq_h, color='grey', ls=':', lw=1.2)
    ax[3].text(nyq_h * 1.03, max(np.nanmax(pH), np.nanmax(pA)) * 0.35,
               f'hires true\nNyquist\n{nyq_h:.3f} cyc/µm', fontsize=7.5, va='top')
    ax[3].set_xlabel('spatial frequency (cycles/µm)')
    ax[3].set_ylabel('radial power (per-pixel, same field)')
    txt = '\n'.join([f'{a:.3f}-{b:.3f}: rho {rr:+.2f}' for a, b, rr in brho])
    ax[3].text(0.02, 0.98, 'band rho (hires vs aligned)\n' + txt, transform=ax[3].transAxes,
               fontsize=7, va='top', family='monospace',
               bbox=dict(fc='white', ec='grey', alpha=0.8))
    ax[3].set_title(f'redo at the MEASURED {up:.2f}x (old test used 3x — void)\n'
                    f'>Nyquist power ratio control/aligned = {ratio:.3f}')
    ax[3].legend(fontsize=7.5, loc='lower left')
    for a in ax[:3]:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f'{gsm}: does aligned_tissue contain detail that hires cannot?  '
                 f'[greyscale, so unusable by PLIP regardless of the answer]', fontsize=12)
    fig.tight_layout()
    fp = os.path.join(OUT, 'aligned_detail.png')
    fig.savefig(fp, dpi=105)
    plt.close(fig)
    print(f'[out] {fp}')

    with open(os.path.join(OUT, 'aligned_detail.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(OUT, 'aligned_detail_summary.json'), 'w') as fh:
        json.dump(dict(q_field=Q_FIELD, um_aligned_native=round(um_a, 4),
                       up_factor=round(up, 4), hires_nyquist=round(nyq_h, 5),
                       control='fft ideal upsample (zero power > Nyquist)',
                       psd_ratio_below_nyquist=round(ratio_lo, 4),
                       psd_ratio_above_nyquist=round(ratio, 4),
                       band_rho={f'{a:.3f}-{b:.3f}': (None if np.isnan(rr) else round(rr, 4))
                                 for a, b, rr in brho},
                       jpeg_block_ratio=round(jb, 4), subpixel_band=sub,
                       median_pixel_rho=float(np.nanmedian([r2['rho_pixel'] for r2 in rows])),
                       example_slide=gsm), fh, indent=2, ensure_ascii=False)
    print(f'[out] {OUT}/aligned_detail.csv , aligned_detail_summary.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
