#!/usr/bin/env python3
# ============================================================================
# ⛔ 已作废（RETIRED，2026-09-24）—— **不要运行，不要引用它的任何输出**
# ----------------------------------------------------------------------------
# 作废理由（由 `16_source_audit.py` 给出硬证据，见
# `results/07_he_pathology/resolution_sweep/source_audit.png`）：
#
#   1. **目的消失了**：本脚本存在的唯一目的是把 `cytassist_image.tiff` 配准到 spot 坐标系，
#      好拿它当「比 hires 更细的彩色源」重跑 PLIP。但实测 cytassist 覆盖的视野**比 hires 更宽**
#      （hires 组织占画面 33%，cytassist 只占 62% 且掩码已被载物台边缘污染），
#      ⇒ cytassist **不是**更细的源，是**更粗的**。既然不更细，配准就没有意义。
#      `detected_tissue_image.jpg` 上画出的 spot 点阵自相关基周期 = 10.7 px，
#      把它换算成 µm/px 落在 4.67–9.35 区间，而 hires = 5.670 µm/px：
#      **cytassist 与 hires 同级或更粗，不存在「更细的彩色源」。**
#
#   2. **方法本身也不成立**：本脚本试过的四种配准（spot 点云 vs 组织掩码 FFT /
#      Otsu 掩码扫尺度 / 色度掩码扫尺度 / 外接框定尺度+平移）全部退化：
#      均质地毯式点阵叠在组织形状的掩码上，缩得越小重叠越"好"，NCC 归一化分母是常数
#      ⇒ 峰值永远跑到搜索边界。实测输出把 spot 推到画面外（在框率 1.5–32%），
#      亮度相关 rho ≈ 0。**没有产出过任何可用变换。**
#
#   3. deposit 里**没有**配准矩阵（`spatial/` 下只有图 + positions + scalefactors；
#      `spatial_enrichment.csv` 是基因富集表）。**变换只能自己测，而我们测不出来。**
#
# ⇒ 分辨率判决实验改为**只用 hires 的裁框阶梯**，见 `PREREG_v2.md`。
# ============================================================================
# HE 病理模型 · 第十五步：把 `cytassist_image.tiff` 配准到 Spot 坐标系
#
# 为什么需要这一步：
#   分辨率判决实验（14）想用「比 hires 更细的彩色图」重跑 PLIP。deposit 里三张候选：
#     tissue_hires_image.png   彩色  5.69 µm/px   ← 既有基线用的
#     aligned_tissue_image.jpg **灰度（色度 std = 0.00）** 1.897 µm/px  ← 🔴 废
#     cytassist_image.tiff     彩色  约 3.79 µm/px ← 唯一可用的更细彩色图
#   `aligned_tissue_image.jpg` 是灰度，而 PLIP 判别**依赖颜色**（实测：把裁框转灰度
#   再复制三通道，PLIP 判成「空玻璃片」，empty 0.888 vs 真彩 0.370）⇒ 不能用。
#
#   Space Ranger 没有在 deposit 里留下任何配准矩阵（spatial/ 下只有图 + positions
#   + scalefactors），所以 cytassist → Spot 坐标系的变换必须**我们自己测**。
#
# 方法：hires 与 cytassist 是**同一张组织的彩色 H&E**，拿组织轮廓做刚性配准 ——
#   两图各自 Otsu 出组织掩码，扫尺度 a（cyt_px / hires_px），
#   每个 a 用 FFT 互相关求最佳平移，取归一化互相关峰值最高的 a。
#   轮廓有真实结构（肺叶、裂隙），相位相关能锁住；纯点云重叠锁不住（已试过，退化）。
#
# 产出：results/07_he_pathology/resolution_sweep/cytassist_registration.csv
#       results/07_he_pathology/resolution_sweep/cytassist_registration.png（人工复核用）
#
# 跑法：python3 07_he_pathology/15_register_cytassist.py

import os, re, json, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Image.MAX_IMAGE_PIXELS = None

ROOT = '/home/eto/luad_v2'
DATA = f'{ROOT}/data/visium_spatial'
OUT = f'{ROOT}/results/07_he_pathology/resolution_sweep'
os.makedirs(OUT, exist_ok=True)

WORK = 750           # 配准搜索的工作边长（cyt 缩到 WORK）
Q = 4                # hires 缩放 = 原尺寸 × a / Q
A_RANGE = (0.75, 2.10)   # cyt_px / hires_px 的搜索范围（不看任何先验，纯扫）
N_A = 55
WIN = 60             # 局部微调的搜索半径（cyt 工作网格 px）
FIG_SLIDES = 6       # 复核图展示几张


def tissue_mask(rgb, q):
    """H&E 组织专属掩码：靠**色度**(R−G) 认出粉色组织。
    不用灰度阈值 —— cyt 图里黑色载物台背景、灰色 fiducial 点阵、白色框内空白
    在灰度上都可能和组织混在一起，但色度上分得干净（组织 R−G 明显 > 0）。"""
    a = np.asarray(rgb, dtype=np.int16)
    m = ((a[:, :, 0] - a[:, :, 1]) > 12).astype(np.float32)
    if q > 1:
        im = Image.fromarray((m * 255).astype(np.uint8))
        m = np.asarray(im.resize((m.shape[1] // q, m.shape[0] // q), Image.BILINEAR),
                       dtype=np.float32) / 255.0
    return m


def xcorr_full(A, B):
    """零填充全互相关：返回 (peak_score, dy, dx)，(dy,dx) 是 A 相对 B 需要的平移。"""
    A = A - A.mean(); B = B - B.mean()
    sA, sB = A.std(), B.std()
    if sA < 1e-9 or sB < 1e-9:
        return -1.0, 0, 0
    s = (A.shape[0] + B.shape[0] - 1, A.shape[1] + B.shape[1] - 1)
    FA = np.fft.rfft2(A, s=s); FB = np.fft.rfft2(B, s=s)
    c = np.fft.irfft2(FA * np.conj(FB), s=s)
    p = int(c.argmax())
    py, px = np.unravel_index(p, c.shape)
    dy = py - (B.shape[0] - 1)
    dx = px - (B.shape[1] - 1)
    return float(c.max() / (sA * sB * A.size)), int(dy), int(dx)


def parse_dir(d):
    m = re.match(r'(GSM\d+)_(P\d+)_(Normal|AAH|AIS|MIA|LUAD)', d)
    return (m.group(1), m.group(2), m.group(3)) if m else None


def otsu(g):
    """Otsu 阈值。组织亮、背景更亮，这里按'暗于阈值 = 组织'取。"""
    h = np.bincount(g.ravel(), minlength=256).astype(np.float64)
    tot = h.sum()
    w0 = np.cumsum(h); w1 = tot - w0
    idx = np.arange(256)
    m0 = np.cumsum(h * idx); m1 = (h * idx).sum() - m0
    with np.errstate(invalid='ignore', divide='ignore'):
        var = (m0 / w0 - m1 / w1) ** 2 * w0 * w1
    var[~np.isfinite(var)] = -1
    return int(var.argmax())


def xcorr_peak(A, B):
    """归一化互相关（FFT）。返回 (score, dy, dx)。A/B 同尺寸，零均值化。"""
    A = A - A.mean(); B = B - B.mean()
    sA, sB = A.std(), B.std()
    if sA < 1e-9 or sB < 1e-9:
        return -1.0, 0, 0
    c = np.fft.irfft2(np.fft.rfft2(A) * np.conj(np.fft.rfft2(B)), s=A.shape)
    p = int(c.argmax())
    dy, dx = np.unravel_index(p, c.shape)
    H, W = A.shape
    dy = dy if dy < H // 2 else dy - H
    dx = dx if dx < W // 2 else dx - W
    return float(c.max() / (sA * sB * A.size)), int(dy), int(dx)


def bbox_of(m, pct=0.5):
    """组织掩码的鲁棒外接框（用分位数，避免少量杂点把框撑开）。"""
    ys, xs = np.nonzero(m > 0.5)
    if len(xs) < 100:
        return None
    x0, x1 = np.percentile(xs, pct), np.percentile(xs, 100 - pct)
    y0, y1 = np.percentile(ys, pct), np.percentile(ys, 100 - pct)
    return x0, x1, y0, y1


def register(gdir):
    """hires → cyt 帧：x_c = x_h * a + dx, y_c = y_h * a + dy（各图自身像素）。
    `a` = cyt_px / hires_px = hires_µm_per_px / cyt_µm_per_px。

    两步：① 组织外接框直接给出尺度与平移（外接框是强特征，不会退化）；
         ② 在预测位置 ±40 px 内做局部互相关，报告该处的 NCC 与搜索到的最高 NCC。
    不做全局尺度扫描 —— 空白画布的零填充会让峰值跑到边界（已试过，退化）。"""
    sp = f'{DATA}/{gdir}/spatial'
    hr = Image.open(f'{sp}/tissue_hires_image.png').convert('RGB')
    cy = Image.open(f'{sp}/cytassist_image.tiff').convert('RGB')
    hw, hh = hr.size
    cw, ch = cy.size

    q_c = max(1, int(round(max(cw, ch) / WORK)))
    Mc = tissue_mask(cy, q_c)
    Mh = tissue_mask(hr, 1)
    bh, bc = bbox_of(Mh), bbox_of(Mc)
    if bh is None or bc is None:
        raise RuntimeError('组织掩码为空')

    a = float(np.median([(bc[1] - bc[0]) / (bh[1] - bh[0]) * q_c,
                         (bc[3] - bc[2]) / (bh[3] - bh[2]) * q_c]))
    dx = float((bc[0] + bc[1]) / 2 * q_c - a * (bh[0] + bh[1]) / 2)
    dy = float((bc[2] + bc[3]) / 2 * q_c - a * (bh[2] + bh[3]) / 2)

    # ② 只做**局部**验证：在 ±WIN 内微调（必须限制在预测附近 ——
    #    不加限制的全域互相关对二值掩码会退化到零填充边界，已实测）
    m = np.asarray(Image.fromarray((Mh * 255).astype(np.uint8)).resize(
        (max(1, int(round(hw / q_c))), max(1, int(round(hh / q_c)))), Image.BILINEAR),
        dtype=np.float32) / 255.0
    Hc, Wc = Mc.shape
    py, px = int(round(dy / q_c)), int(round(dx / q_c))
    y0, x0 = max(0, py - WIN), max(0, px - WIN)
    y1, x1 = min(Hc, py + m.shape[0] + WIN), min(Wc, px + m.shape[1] + WIN)
    sub_c = Mc[y0:y1, x0:x1]
    sub_h = np.zeros_like(sub_c)
    ay0, ax0 = py - y0, px - x0
    ty0, tx0 = max(0, ay0), max(0, ax0)
    ty1, tx1 = min(sub_c.shape[0], ay0 + m.shape[0]), min(sub_c.shape[1], ax0 + m.shape[1])
    if ty1 > ty0 and tx1 > tx0:
        sub_h[ty0:ty1, tx0:tx1] = m[ty0 - ay0:ty1 - ay0, tx0 - ax0:tx1 - ax0]
    sc, dy_l, dx_l = xcorr_full(sub_h, sub_c)
    return dict(gsm_dir=gdir, ncc=float(sc), a=a, cyt_px_per_hires_px=a,
                dx_cyt=float(dx + dx_l * q_c), dy_cyt=float(dy + dy_l * q_c),
                dx_bbox=dx, dy_bbox=dy,
                refine_px=float(max(abs(dx_l), abs(dy_l)) * q_c),
                hires_size=f'{hw}x{hh}', cyt_size=f'{cw}x{ch}',
                hires_um_per_px=float(HIRES_UM), cyt_um_per_px=float(HIRES_UM / a))


HIRES_UM = 5.69    # 由 448px 裁框 = 2550µm 与 scalefactors 定出的既有口径（仅供填表）


def spot_check(gdir, a, dx_c, dy_c):
    """配准自检：把 spot 中心映射到 cyt，量 425µm 视野里的组织占比，
    与 hires 同 spot 的小窗口亮度做相关 —— 配准对，两者必然强相关。"""
    sp = f'{DATA}/{gdir}/spatial'
    sf = json.load(open(f'{sp}/scalefactors_json.json'))
    ths = sf['tissue_hires_scalef']
    pos = pd.read_csv(f'{sp}/tissue_positions.csv')
    pos = pos[pos.in_tissue == 1].sort_values('barcode').reset_index(drop=True)
    if len(pos) > 1500:
        rng = np.random.default_rng(0)
        pos = pos.iloc[np.sort(rng.choice(len(pos), 1500, replace=False))].reset_index(drop=True)
    xh = pos.pxl_col_in_fullres.values * ths
    yh = pos.pxl_row_in_fullres.values * ths
    xc = xh * a + dx_c
    yc = yh * a + dy_c

    gh = np.asarray(Image.open(f'{sp}/tissue_hires_image.png').convert('L'), dtype=np.float32)
    gc = np.asarray(Image.open(f'{sp}/cytassist_image.tiff').convert('L'), dtype=np.float32)
    W_h, H_h = gh.shape[1], gh.shape[0]
    W_c, H_c = gc.shape[1], gc.shape[0]
    r = 13   # 约 100 µm 的窗口
    mh, mc = [], []
    for i in range(len(pos)):
        x0, y0 = int(xh[i]) - r, int(yh[i]) - r
        if 0 <= x0 and 0 <= y0 and x0 + 2 * r < W_h and y0 + 2 * r < H_h:
            mh.append(gh[y0:y0 + 2 * r, x0:x0 + 2 * r].mean())
        else:
            mh.append(np.nan)
        x0, y0 = int(xc[i]) - r, int(yc[i]) - r
        if 0 <= x0 and 0 <= y0 and x0 + 2 * r < W_c and y0 + 2 * r < H_c:
            mc.append(gc[y0:y0 + 2 * r, x0:x0 + 2 * r].mean())
        else:
            mc.append(np.nan)
    mh, mc = np.array(mh), np.array(mc)
    ok = np.isfinite(mh) & np.isfinite(mc)
    return float(np.corrcoef(mh[ok], mc[ok])[0, 1]), float(ok.mean())


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()

    slides = []
    for d in sorted(os.listdir(DATA)):
        p = parse_dir(d)
        if p and os.path.isfile(f'{DATA}/{d}/spatial/cytassist_image.tiff'):
            slides.append((*p, d))
    if args.limit:
        step = max(1, len(slides) // args.limit)
        slides = slides[::step][:args.limit]
    print(f'[info] 处理 {len(slides)} 张切片（共 {len(os.listdir(DATA))} 个目录）')

    t0 = time.time()
    rows = []
    for i, (gsm, pat, stage, gdir) in enumerate(slides, 1):
        try:
            r = register(gdir)
            rho, okf = spot_check(gdir, r['a'], r['dx_cyt'], r['dy_cyt'])
            r['spot_rho_100um'] = rho
            r['spot_in_bounds'] = okf
            r['patient'], r['stage'] = pat, stage
            rows.append(r)
            print(f'  [{i:2d}/{len(slides)}] {gdir:24s} '
                  f'a={r["a"]:.4f} cyt={r["cyt_um_per_px"]:.2f}µm/px '
                  f'平移=({r["dx_cyt"]:+.0f},{r["dy_cyt"]:+.0f}) NCC={r["ncc"]:.3f} '
                  f'| 亮度相关 rho={rho:.3f} 在框={okf:.1%}', flush=True)
        except Exception as e:
            print(f'  [{i:2d}/{len(slides)}] {gdir} 失败: {type(e).__name__} {e}', flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(f'{OUT}/cytassist_registration.csv', index=False)
    print(f'\n[out] {OUT}/cytassist_registration.csv  ({len(df)} 行, {(time.time()-t0)/60:.1f} min)')
    if len(df):
        print(f'[stat] a 中位 {df.a.median():.4f} 范围 {df.a.min():.3f}–{df.a.max():.3f}')
        print(f'[stat] cyt µm/px 中位 {df.cyt_um_per_px.median():.2f}  范围 '
              f'{df.cyt_um_per_px.min():.2f}–{df.cyt_um_per_px.max():.2f}')
        print(f'[stat] 亮度相关 rho 中位 {df.spot_rho_100um.median():.3f}  '
              f'最差 {df.spot_rho_100um.min():.3f}  (n={len(df)})')

    # ——— 人工复核图：spot 中心叠在两张图上 ———
    d = df.sort_values('spot_rho_100um', ascending=False)
    pick = list(d.head(FIG_SLIDES // 2).gsm_dir) + list(d.tail(FIG_SLIDES - FIG_SLIDES // 2).gsm_dir)
    pick = pick[:FIG_SLIDES]
    fig, axes = plt.subplots(2, len(pick), figsize=(3.4 * len(pick), 7.2))
    if len(pick) == 1:
        axes = axes.reshape(2, 1)
    for j, gdir in enumerate(pick):
        sp = f'{DATA}/{gdir}/spatial'
        r = df[df.gsm_dir == gdir].iloc[0]
        sf = json.load(open(f'{sp}/scalefactors_json.json'))
        pos = pd.read_csv(f'{sp}/tissue_positions.csv')
        pos = pos[pos.in_tissue == 1]
        for k, (img_f, sc, dx, dy, name) in enumerate([
                ('tissue_hires_image.png', sf['tissue_hires_scalef'], 0, 0, 'hires (baseline)'),
                ('cytassist_image.tiff', sf['tissue_hires_scalef'] * r.a, r.dx_cyt, r.dy_cyt,
                 'cytassist (registered)')]):
            im = Image.open(f'{sp}/{img_f}').convert('RGB')
            im.thumbnail((600, 600))
            s = im.size[0] / (Image.open(f'{sp}/{img_f}').size[0])
            axes[k, j].imshow(im)
            axes[k, j].scatter(pos.pxl_col_in_fullres * sc * s + dx * s,
                               pos.pxl_row_in_fullres * sc * s + dy * s,
                               s=0.4, c='lime', linewidths=0)
            axes[k, j].set_title(f'{gdir}\n{name}', fontsize=7)
            axes[k, j].axis('off')
    fig.suptitle('Spot centres overlaid — hires (top) vs registered cytassist (bottom)', fontsize=9)
    fig.tight_layout()
    fig.savefig(f'{OUT}/cytassist_registration.png', dpi=130)
    print(f'[out] {OUT}/cytassist_registration.png')


if __name__ == '__main__':
    main()
