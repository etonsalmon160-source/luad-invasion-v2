#!/usr/bin/env python3
# HE 病理模型 · 第十六步：图像源审计
#
# ⛔ 本脚本的两个函数已作废（2026-09-24 实测更正），**不得引用其结果**：
#   · `correspondence_test()` —— 它是**整幅直接对位**，内含两个**从未验证**的假设：
#     ① aligned 与 hires 像素数差 3 倍 ⇒ 同视野；② 两者共原点。
#     实测两条都错：hires 是**紧贴组织的裁切**、aligned 是**宽视野**（视野比 1.30），
#     1/3 缩放根本没把两者放到同一物理尺度 ⇒ 它给出的「内容对不上」（NCC 0.16–0.18）
#     是**方法错**，不是数据性质。**已撤回。**
#   · `highfreq_test()` —— 它拿「hires 放大 **3 倍**」当无新信息对照，实测比例是 **2.31 倍**
#     ⇒ 对照项选错，图 `source_highfreq.png` 的对照**作废**
#     （「aligned 有方块状拼接伪影」这一观察仍在，但结论不下）。
#
# ✅ 现行测法见 `07_he_pathology/16b_source_correspondence.py`（轮廓归一化 + 边缘图比例搜索）
#    → `source_correspondence.png`。结论：**同一块组织、同一朝向**（4/4 最佳摆法 as-is，
#    镜像最差），**视野比 q=1.30** ⇒ aligned ≈ **2.457 µm/px**（比 hires 细 2.31 倍）。
#    旧值「1.890 µm/px、细 3 倍」**作废**。
#
# 本脚本其余部分（灰度/色度测量、帧归属 NCC、点阵周期、径向 PSD 本身）仍然有效。
#
# 为什么有这一步：
#   分辨率判决实验（14）的原设计用的是 `aligned_tissue_image.jpg`（5985×6000，1.890 µm/px，
#   官方 `regist_target_img_scalef` 给出到 spot 的精确映射）。实测该图**是灰度的**
#   （R−G 标准差 = 0.000），而 PLIP 判别**依赖颜色**：灰度裁框会被读成「空玻璃片」
#   （empty 0.889 vs 真彩 0.370）。⇒ 原设计的 F1/F2 两个臂**跑不了**。
#
#   于是必须盘清楚：deposit 里还有没有更细的**彩色** HE 图？
#   候选只有 `cytassist_image.tiff`（3000×3000）。本脚本把它审计到能下结论为止。
#
# 本脚本**只做测量**，不判恶性、不产标签、不改任何既有结果。
#
# 跑法：python3 07_he_pathology/16_source_audit.py

import os, json, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Image.MAX_IMAGE_PIXELS = None

ROOT = '/home/eto/luad_v2'
DATA = f'{ROOT}/data/visium_spatial'
OUT = f'{ROOT}/results/07_he_pathology/resolution_sweep'
os.makedirs(OUT, exist_ok=True)

TIERS = ['cytassist_image.tiff', 'detected_tissue_image.jpg',
         'aligned_tissue_image.jpg', 'aligned_fiducials.jpg',
         'tissue_hires_image.png', 'tissue_lowres_image.png']

# 阵列几何：array_col 从 0 到 223 共 223 个「单位步」，全阵列在 fullres 里的跨度
# 实测为 44519 px。Visium 的步长是 50 µm ⇒ fullres = 50/199.6 = 0.2505 µm/px。
FULLRES_UM_PER_PX = 0.25050


def measure(f, gsm_dir):
    im = Image.open(f)
    a = np.asarray(im.convert('RGB'), dtype=np.int16)
    rg = (a[:, :, 0] - a[:, :, 1])
    return dict(gsm_dir=gsm_dir, file=os.path.basename(f), w=im.size[0], h=im.size[1],
                chroma_std=float(rg.std()), r_mean=float(a[:, :, 0].mean()),
                g_mean=float(a[:, :, 1].mean()), b_mean=float(a[:, :, 2].mean()))


def gray_norm(f, size):
    im = Image.open(f).convert('L').resize((size, size), Image.BILINEAR)
    a = np.asarray(im, dtype=np.float32)
    return (a - a.mean()) / (a.std() + 1e-9)


def lattice_pitch(f, box):
    """画出格点的自相关周期：给 detected_tissue 上的 spot 点阵量周期（px）。
    该图是 Space Ranger 把 spot 格点**画在 cytassist 图上**的产物 ⇒ 点阵周期是硬事实。"""
    g = np.asarray(Image.open(f).convert('L'), dtype=np.float32)
    sub = g[box[0]:box[1], box[2]:box[3]]
    sub = sub - sub.mean()
    F = np.fft.rfft2(sub)
    ac = np.fft.fftshift(np.fft.irfft2(F * np.conj(F), s=sub.shape))
    c = np.array(ac.shape) // 2
    prof = ac[c[0], c[1]:c[1] + 120]
    peaks = [i for i in range(3, 119) if prof[i] > prof[i - 1] and prof[i] >= prof[i + 1]]
    return prof, peaks


def tissue_bbox(f, thresh=12, pct=0.5):
    a = np.asarray(Image.open(f).convert('RGB'), dtype=np.int16)
    m = ((a[:, :, 0] - a[:, :, 1]) > thresh)
    ys, xs = np.nonzero(m)
    return (float(np.percentile(xs, pct)), float(np.percentile(xs, 100 - pct)),
            float(np.percentile(ys, pct)), float(np.percentile(ys, 100 - pct)),
            float(m.mean()))


def radial_psd(g, um_per_px, nb=48):
    """Hann 窗 + 2D FFT ⇒ 径向平均功率谱，横轴 cycles/µm。
    用途：判断一张图在它自己的像素网格上**到底含多少真实细节**。"""
    g = np.asarray(g, dtype=np.float32)
    g = g - g.mean()
    ny, nx = g.shape
    gw = g * np.hanning(ny)[:, None] * np.hanning(nx)[None, :]
    F = np.fft.fftshift(np.fft.fft2(gw))
    P = (np.abs(F) ** 2) / (ny * nx)
    cy, cx = P.shape[0] // 2, P.shape[1] // 2
    yy, xx = np.mgrid[0:P.shape[0], 0:P.shape[1]]
    k = np.sqrt(((yy - cy) / (ny * um_per_px)) ** 2 + ((xx - cx) / (nx * um_per_px)) ** 2)
    kmax = 0.5 / um_per_px
    edges = np.linspace(0, kmax, nb + 1)
    idx = np.digitize(k.ravel(), edges) - 1
    ok = (idx >= 0) & (idx < nb)
    s = np.bincount(idx[ok], weights=P.ravel()[ok], minlength=nb)
    c = np.bincount(idx[ok], minlength=nb)
    prof = np.where(c > 0, s / np.maximum(c, 1), np.nan)
    return 0.5 * (edges[:-1] + edges[1:]), prof


def correspondence_test(sp):
    """⛔ **已作废，不得引用**（见文件头横幅）。

    `aligned_tissue_image.jpg` 与 `tissue_hires_image.png` 内容对得上吗？

    它把 aligned 按 **1/3** 缩到 hires 格子上再比 —— 这**默默假设了**
    「像素数差 3 倍 ⇒ 同视野」和「两者共原点」。实测两条都错
    （hires 紧贴组织、aligned 宽视野，视野比 **1.30**），
    ⇒ 它给出的「内容对不上」是**方法错**，已撤回。
    现行测法：`16b_source_correspondence.py`。
    """
    raise NotImplementedError('已作废：整幅直接对位内含未验证的"同视野/同原点"假设；'
                              '改用 16b_source_correspondence.py')
    from scipy.signal import fftconvolve
    A0 = np.asarray(Image.open(f'{sp}/aligned_tissue_image.jpg').convert('L'), dtype=np.float32)
    H = np.asarray(Image.open(f'{sp}/tissue_hires_image.png').convert('L'), dtype=np.float32)
    h3, w3 = (A0.shape[0] // 3) * 3, (A0.shape[1] // 3) * 3
    A = A0[:h3, :w3].reshape(h3 // 3, 3, w3 // 3, 3).mean(axis=(1, 3))
    n, m = min(A.shape[0], H.shape[0]), min(A.shape[1], H.shape[1])
    A, H = A[:n, :m], H[:n, :m]

    def nx(x):
        return (x - x.mean()) / (x.std() + 1e-9)

    An, Hn = nx(A), nx(H)
    flat = float((An * Hn).mean())
    c = fftconvolve(An, Hn[::-1, ::-1], mode='same')
    shift = float(c.max() / An.size)
    Ac = nx(A[n // 4:n // 4 + 1200, m // 4:m // 4 + 1200])
    c2 = fftconvolve(Ac, Hn[::-1, ::-1], mode='valid')
    inner = float(c2.max() / Ac.size)

    dih = {}
    for name, T in [('as-is', A), ('flip-LR', A[:, ::-1]), ('flip-UD', A[::-1, :]),
                    ('rot180', A[::-1, ::-1]), ('rot90', np.rot90(A, 1)),
                    ('rot270', np.rot90(A, 3))]:
        B = np.asarray(Image.fromarray(T.astype(np.uint8)).resize((m, n), Image.BILINEAR),
                       dtype=np.float32) if T.shape != A.shape else T
        dih[name] = float((nx(B) * Hn).mean())
    return dict(ncc_asis=flat, ncc_best_shift=shift, ncc_inner_shift=inner, dihedral=dih)


def highfreq_test(sp, spot_xy_hires, um_hires=5.670):
    """⛔ **已作废，不得引用**（见文件头横幅）。

    原意：aligned_tissue（当时以为 1.890 µm/px）是不是真有比 hires 更细的细节。
    对照 B 取的是「hires 放大 **3 倍**」—— 建立在「两者像素数差 3 倍 ⇒ 同视野」这个
    **未验证假设**上。实测视野比是 **1.30** ⇒ 正确对照应是 **2.31 倍**。
    对照项选错 ⇒ 该函数的输出与 `source_highfreq.png` 的对照**作废**。
    """
    raise NotImplementedError('已作废：对照项按 3 倍选错（实测 2.31 倍）。'
                              '要重测须另跑，见 16b_source_correspondence.py')
    hx, hy = spot_xy_hires
    half_h, half_a = 128, 384     # hires 半宽 / aligned 半宽

    im_h = Image.open(f'{sp}/tissue_hires_image.png').convert('L')
    im_a = Image.open(f'{sp}/aligned_tissue_image.jpg').convert('L')
    H, W = im_h.size[1], im_h.size[0]
    hx = float(np.clip(hx, half_h + 2, W - half_h - 2))
    hy = float(np.clip(hy, half_h + 2, H - half_h - 2))

    ph = np.asarray(im_h.crop((int(hx - half_h), int(hy - half_h),
                               int(hx + half_h), int(hy + half_h))), dtype=np.float32)
    ax_, ay_ = hx * 3.0, hy * 3.0
    pa = np.asarray(im_a.crop((int(ax_ - half_a), int(ay_ - half_a),
                               int(ax_ + half_a), int(ay_ + half_a))), dtype=np.float32)

    pB = np.asarray(Image.fromarray(ph.astype(np.uint8)).resize((half_a * 2, half_a * 2),
                                                                Image.BICUBIC), dtype=np.float32)
    # A 缩回 hires 网格，与 C 逐像素比：是否只是 C 的模糊版？
    pA_small = np.asarray(Image.fromarray(pa.astype(np.uint8)).resize((half_h * 2, half_h * 2),
                                                                      Image.BICUBIC),
                          dtype=np.float32)
    m = float(pA_small.mean())
    rms_res = float(np.sqrt((((pA_small - m) - (ph - ph.mean())) ** 2).mean()))
    rho = float(np.corrcoef(pA_small.ravel(), ph.ravel())[0, 1])

    kA, pA = radial_psd(pa, 1.890)
    kB, pB_ = radial_psd(pB, 1.890)
    kC, pC = radial_psd(ph, um_hires)
    return dict(kA=kA, pA=pA, kB=kB, pB=pB_, kC=kC, pC=pC,
                rho_A3_vs_C=rho, rms_res=float(rms_res / (ph.std() + 1e-9)),
                std_A=float(pa.std()), std_B=float(pB.std()), std_C=float(ph.std()),
                nyq_hires=0.5 / um_hires, nyq_aligned=0.5 / 1.890,
                patch_hires=ph, patch_aligned=pa, patch_hires_up=pB)


def main():
    slides = sorted(d for d in os.listdir(DATA) if os.path.isdir(f'{DATA}/{d}/spatial'))
    gdir = slides[0]
    sp = f'{DATA}/{gdir}/spatial'
    sf = json.load(open(f'{sp}/scalefactors_json.json'))

    rows = [measure(f'{sp}/{f}', gdir) for f in TIERS if os.path.isfile(f'{sp}/{f}')]
    df = pd.DataFrame(rows)
    df['um_per_px'] = FULLRES_UM_PER_PX * (df.w / df.w)  # 占位，下面逐行填
    frw = df.loc[df.file == 'tissue_hires_image.png', 'w'].iloc[0] / sf['tissue_hires_scalef']
    df['field_mm'] = df.w * FULLRES_UM_PER_PX * (frw / frw)  # 占位
    df.to_csv(f'{OUT}/source_audit_tiers.csv', index=False)

    # ——— 帧归属：两两内容 NCC（同一帧的两张图 resize 到同尺寸后必然强相关）———
    keys = ['aligned_tissue_image.jpg', 'cytassist_image.tiff', 'aligned_fiducials.jpg',
            'detected_tissue_image.jpg', 'tissue_hires_image.png', 'tissue_lowres_image.png']
    keys = [k for k in keys if os.path.isfile(f'{sp}/{k}')]
    g = {k: gray_norm(f'{sp}/{k}', 800) for k in keys}
    ncc = pd.DataFrame(index=keys, columns=keys, dtype=float)
    for i in keys:
        for j in keys:
            ncc.loc[i, j] = float((g[i] * g[j]).mean())
    ncc.to_csv(f'{OUT}/source_audit_ncc.csv')
    # ⚠️ 本表是把每张图**各自归一化到 800×800** 后算的整幅相关 ⇒
    #    它对「同帧、同视野」的图对（cytassist↔detected↔fiducials、hires↔lowres）有效；
    #    对 `aligned_tissue ↔ hires` **无效**（两者视野比 1.30，见 16b），该格**作废**。

    # ——— detected_tissue 上点阵的周期 ———
    prof, peaks = lattice_pitch(f'{sp}/detected_tissue_image.jpg', (400, 1000, 500, 1100))

    # ——— 组织外接框（cyt 的色度掩码会被载物台边缘污染，必须登记污染率）———
    bh = tissue_bbox(f'{sp}/tissue_hires_image.png')
    bc = tissue_bbox(f'{sp}/cytassist_image.tiff')

    print(f'[slide] {gdir}')
    print(f'[anchor] fullres {FULLRES_UM_PER_PX} µm/px ; '
          f'hires {FULLRES_UM_PER_PX/sf["tissue_hires_scalef"]:.3f} µm/px')
    print(f'[anchor] ⚠️ aligned 的 µm/px **不能**用 '
          f'fullres/regist_target_img_scalef = {FULLRES_UM_PER_PX/sf["regist_target_img_scalef"]:.3f} 算 ——')
    print(f'         那只在"两边同视野"成立时才对，而实测视野比 1.30 ⇒ aligned ≈ **2.457 µm/px**'
          f'（16b 实测）。')
    print(df[['file', 'w', 'h', 'chroma_std', 'r_mean']].to_string(index=False))
    print(f'\n[NCC] aligned_tissue vs hires = {ncc.loc["aligned_tissue_image.jpg","tissue_hires_image.png"]:+.3f}'
          f'  ⚠️ 此值按 1/3 缩放算，比例假设已证伪 ⇒ **作废**，勿引')
    print(f'[NCC] cytassist vs detected_tissue = {ncc.loc["cytassist_image.tiff","detected_tissue_image.jpg"]:+.3f}')
    print(f'[NCC] cytassist vs aligned_fiducials = {ncc.loc["cytassist_image.tiff","aligned_fiducials.jpg"]:+.3f}')
    print(f'[NCC] hires vs lowres = {ncc.loc["tissue_hires_image.png","tissue_lowres_image.png"]:+.3f}')
    print(f'\n[点阵] detected_tissue 自相关峰位置(px) = {peaks}')
    print(f'[点阵] ⇒ 周期 {np.mean(np.diff(peaks)) if len(peaks)>1 else -1:.1f} px '
          f'⇒ 若该周期对应 100 µm 则 cyt ≈ {100/np.mean(np.diff(peaks)):.2f} µm/px')
    print(f'\n[组织框] hires   x[{bh[0]:.0f},{bh[1]:.0f}] y[{bh[2]:.0f},{bh[3]:.0f}] 掩码占比 {bh[4]:.1%}')
    print(f'[组织框] cyt     x[{bc[0]:.0f},{bc[1]:.0f}] y[{bc[2]:.0f},{bc[3]:.0f}] 掩码占比 {bc[4]:.1%}'
          f'  ← 占比过高 ⇒ 载物台边缘污染，框不可用')

    # ——— 高频检验：aligned_tissue 到底有没有比 hires 更细的真细节 ———
    pos = pd.read_csv(f'{sp}/tissue_positions.csv')
    pos = pos[pos.in_tissue == 1]
    hx = float(np.median(pos.pxl_col_in_fullres.values) * sf['tissue_hires_scalef'])
    hy = float(np.median(pos.pxl_row_in_fullres.values) * sf['tissue_hires_scalef'])
    # ⛔ [对应] 与 [高频图] 两段已作废（见文件头横幅）：
    #   correspondence_test() 用整幅直接对位，内含"同视野/同原点"两个未验证假设；
    #   highfreq_test() 的对照项按 3 倍选，实测比例是 2.31 倍。
    #   两者都**不再运行、不再出图**，留着只会误导。现行测法见 16b。
    print('\n[⛔ 作废] [对应]/[高频图] 两段已停跑（整幅直接对位 + 3 倍对照都是错的）；')
    print('         现行测法见 07_he_pathology/16b_source_correspondence.py '
          '→ source_correspondence.png')
    print('         结论：同一块组织、同一朝向；视野比 q=1.30 ⇒ aligned ≈ 2.457 µm/px（细 2.31 倍）')

    # ————————————— 审计图 —————————————
    fig = plt.figure(figsize=(15.5, 8.4))
    gs = fig.add_gridspec(2, 4)
    panels = [('cytassist_image.tiff', 'cytassist_image.tiff — whole slide, COLOUR'),
              ('detected_tissue_image.jpg', 'detected_tissue_image.jpg — same frame, spot grid drawn'),
              ('aligned_tissue_image.jpg', 'aligned_tissue_image.jpg — GRAYSCALE (chroma std = 0.000)'),
              ('tissue_hires_image.png', 'tissue_hires_image.png — the trusted source')]
    for k, (f, ttl) in enumerate(panels):
        ax = fig.add_subplot(gs[0, k])
        im = Image.open(f'{sp}/{f}').convert('RGB')
        im.thumbnail((520, 520))
        ax.imshow(im)
        ax.set_title(ttl, fontsize=7.5)
        ax.axis('off')

    ax = fig.add_subplot(gs[1, 0:2])
    ax.plot(prof, lw=0.9, color='crimson')
    for p in peaks:
        ax.axvline(p, color='grey', ls=':', lw=0.8)
    ax.set_title('detected_tissue_image.jpg — autocorrelation of the drawn spot lattice\n'
                 f'peak lags (px): {peaks}', fontsize=8)
    ax.set_xlabel('lag (px)', fontsize=8)
    ax.set_ylabel('autocorr', fontsize=8)
    ax.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[1, 2])
    a = np.asarray(Image.open(f'{sp}/aligned_fiducials.jpg').convert('RGB'), dtype=np.int16)
    m = (a[:, :, 0] > 150) & (a[:, :, 1] < 100) & (a[:, :, 2] < 100)
    ax.imshow(np.asarray(Image.fromarray((m * 255).astype(np.uint8)).resize((600, 600),
                                                                          Image.BILINEAR)),
              cmap='gray_r')
    ax.set_title(f'aligned_fiducials.jpg — red markers only (n={int(m.sum())} px)\n'
                 'whole slide, NO H&E tissue', fontsize=8)
    ax.axis('off')

    ax = fig.add_subplot(gs[1, 3])
    ax.scatter(bh[0:2], [1, 1], s=40, c='steelblue', label='hires tissue x-span')
    ax.scatter(bc[0:2], [0.4, 0.4], s=40, c='crimson', label='cyt "tissue" x-span')
    ax.annotate('', xy=(bh[0], 1), xytext=(bh[1], 1),
                arrowprops=dict(arrowstyle='<->', color='steelblue'))
    ax.annotate('', xy=(bc[0], 0.4), xytext=(bc[1], 0.4),
                arrowprops=dict(arrowstyle='<->', color='crimson'))
    ax.text(1500, 1.15, 'hires 1995 px = 11.31 mm', fontsize=7.5, ha='center')
    ax.text(1500, 0.55, 'cyt 3000 px (mask contaminated 62%)', fontsize=7.5, ha='center')
    ax.set_xlim(0, 3000)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.legend(fontsize=7, loc='lower center')
    ax.set_title('tissue extent — cyt covers a WIDER field than hires\n'
                 '⇒ cytassist is NOT the finer source', fontsize=8)

    fig.suptitle('Source audit: the resolution experiment of PREREG.md cannot be run as designed',
                 fontsize=10, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    fig.savefig(f'{OUT}/source_audit.png', dpi=125)
    print(f'\n[out] {OUT}/source_audit.png')
    print(f'[out] {OUT}/source_audit_tiers.csv , source_audit_ncc.csv')


if __name__ == '__main__':
    main()
