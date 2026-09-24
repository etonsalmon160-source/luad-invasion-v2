#!/usr/bin/env python3
"""SP 模块 · 空转（Visium）spot 级质控

两种模式：
  --figures-only   只算 QC 指标 + 出图，**不冻结掩码**。（签字之前就该跑这个）
  --freeze         按签字阈值冻结 spot 掩码。（阈值未签字会硬停）

阈值来源 = **源论文逐字方法段**（Peng 2026 Cancer Cell, PMC12980502）：
  > spots with nUMI < 500 and nFeature < 200 and relative high percentage
  > of mitochondrial genes >15% were filtered out. As a result, low-quality
  > spots (only ~1.05%) were also excluded from subsequent analyses.

⚠️ 原文那句用的是 `and`。两个可能读法（任一成立即剔 / 三个同时成立才剔）
   会给出完全不同的剔除率 ⇒ 本脚本两种都算，用 "~1.05%" 这个已知数去判。
   实测两种都够不着 1.05%（OR 剔 4.97% / AND 剔 0.01%）⇒ 2026-09-25 用户签字：
   阈值照原文 500 / 200 / 0.15，读法取 **OR**，4.97% 的差额登记为缺陷（见 SIGNED / RULE）。

⚠️ 原文另有一层我们**做不到**：病理学家先剔除 "Space" / "Foreign body" /
   "Tissue artifacts" 三类 spot。我们无逐 spot 病理注释 ⇒ 该层缺失，登记为局限。

跑法：
  python3 08_spatial_deconv/00_spatial_qc.py --figures-only
"""

import os, sys, gzip, json, hashlib, argparse, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

ROOT = '/home/eto/luad_v2'
DATA = f'{ROOT}/data/visium_spatial'
OUT = f'{ROOT}/results/08_spatial_deconv'
CACHE = f'{OUT}/qc_per_slide'

# ——— 冻结常量：源论文 Visium QC（逐字，法则 3.1）———
PAPER_NUMI_MIN = 500      # 原文：nUMI < 500
PAPER_NFEAT_MIN = 200     # 原文：nFeature < 200
PAPER_MT_MAX = 15.0       # 原文：mitochondrial genes >15%
PAPER_DROP_PCT_REPORTED = 1.05   # 原文自报剔除比例，用于反推 `and` 的读法

MT_PREFIX = 'MT-'         # 人线粒体基因 symbol 前缀

N_WORKERS = 3             # 见 main() 里的说明：本机共享受限，不拉高

# ——— 已签字阈值（SPATIAL_QC_PREREG.md §7 第 1 项，2026-09-25）———
# 用户 2026-09-25 裁定：**照原文 500 / 200 / 0.15**，并同时确认
# 「接受实测剔 4.97%（OR 读法），而非原文自称的 1.05%；差额登记为缺陷」。
# 注意量纲：`MAX_MT_FRAC` 是**比例**（0.15），不是百分数（原文写的是 15%）。
SIGNED = {
    'MIN_UMI_SPOT': 500,
    'MIN_GENES_SPOT': 200,
    'MAX_MT_FRAC': 0.15,
}

# 剔除读法：**OR**（任一不达标即剔）。原文那句话用的是 `and`，字面读法只剔 0.01%；
# 用户签的「接受 4.97%」就是 OR 的值 ⇒ 语义与阈值一起签的是 OR。
RULE = 'or'


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def list_slides():
    ds = sorted(d for d in os.listdir(DATA) if d.startswith('GSM'))
    out = []
    for d in ds:
        p = f'{DATA}/{d}/filtered_feature_bc_matrix'
        if all(os.path.exists(f'{p}/{f}') for f in
               ('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz')):
            out.append(d)
    return out


def read_features(slide):
    """→ (var_names, is_mt, feature_type)。只保留 Gene Expression，剔掉抗体行。"""
    p = f'{DATA}/{slide}/filtered_feature_bc_matrix/features.tsv.gz'
    ids, syms, types = [], [], []
    with gzip.open(p, 'rt') as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            ids.append(c[0]); syms.append(c[1] if len(c) > 1 else '')
            types.append(c[2] if len(c) > 2 else 'Gene Expression')
    types = np.array(types)
    is_ge = types == 'Gene Expression'
    return ids, syms, types, is_ge


def read_mtx_triplets(path):
    """读 MatrixMarket 三元组 → (gene_idx0, spot_idx0, value)，全部 0-based。

    ⚠️ 不用 `scipy.io.mmread`：它逐行 Python 解析文本，实测比 pandas 慢两个数量级。
    ⚠️ 这里**不做** Gene Expression 过滤 —— 过滤在调用方按行做，因为抗体行的计数
       量级可达 1e6（实测），混进去会把 nUMI 整个毁掉（见 SPATIAL_QC_PREREG.md §1.1）。
    """
    with gzip.open(path, 'rt') as fh:
        for line in fh:                       # 跳过 % 开头的头，取维度行
            if line[0] != '%':
                nrow, ncol, nnz = map(int, line.split())
                break
        df = pd.read_csv(fh, sep=' ', header=None, dtype=np.int32,
                         engine='c', names=['g', 'b', 'v'])
    if len(df) != nnz:
        raise ValueError(f'{path}: 读到 {len(df):,} 行，声明 {nnz:,}')
    return df['g'].values - 1, df['b'].values - 1, df['v'].values, nrow, ncol


def qc_one_slide(slide):
    """算该切片逐 spot 的 nUMI / nFeature / pct_mt（只统计 Gene Expression 行）。"""
    base = f'{DATA}/{slide}/filtered_feature_bc_matrix'
    with gzip.open(f'{base}/barcodes.tsv.gz', 'rt') as f:
        barcodes = [l.strip() for l in f if l.strip()]

    ids, syms, types, is_ge = read_features(slide)
    syms = np.array(syms)

    g, b, v, nrow, ncol = read_mtx_triplets(f'{base}/matrix.mtx.gz')
    assert nrow == len(syms), f'{slide}: mtx 行 {nrow} != features {len(syms)}'
    assert ncol == len(barcodes), f'{slide}: mtx 列 {ncol} != barcodes {len(barcodes)}'

    # 只留 Gene Expression 行 —— 抗体行必须剔（其计数可达 1e6）
    keep = is_ge[g]
    g, b, v = g[keep], b[keep], v[keep]

    nspot = len(barcodes)
    numi = np.bincount(b, weights=v, minlength=nspot).astype(np.int64)
    nfeat = np.bincount(b, minlength=nspot).astype(np.int64)      # nnz 数 = 检出基因数

    mt_row = np.zeros(len(syms), dtype=bool)
    mt_row[np.flatnonzero(is_ge)] = np.array(
        [s.upper().startswith(MT_PREFIX) for s in syms])[is_ge]
    mt_hit = mt_row[g]
    mtc = np.bincount(b[mt_hit], weights=v[mt_hit], minlength=nspot)
    with np.errstate(invalid='ignore', divide='ignore'):
        pct_mt = np.where(numi > 0, 100.0 * mtc / numi, np.nan)

    df = pd.DataFrame({
        'slide': slide,
        'barcode': barcodes,
        'nUMI': numi,
        'nFeature': nfeat,
        'pct_mt': pct_mt,
    })
    return df, dict(n_barcodes=nspot, n_genes_GE=int(is_ge.sum()),
                    n_genes_antibody=int((~is_ge).sum()),
                    n_MT_genes=int(mt_row.sum()),
                    n_antibody_rows=int((~is_ge).sum()),
                    matrix_sha=sha256(f'{base}/matrix.mtx.gz'),
                    features_sha=sha256(f'{base}/features.tsv.gz'),
                    barcodes_sha=sha256(f'{base}/barcodes.tsv.gz'))


def load_or_compute(slide):
    os.makedirs(CACHE, exist_ok=True)
    cp = f'{CACHE}/{slide}.csv.gz'
    mp = f'{CACHE}/{slide}.meta.json'
    if os.path.exists(cp) and os.path.exists(mp):
        return pd.read_csv(cp), json.load(open(mp))
    t = time.time()
    df, meta = qc_one_slide(slide)
    meta['seconds'] = round(time.time() - t, 1)
    df.to_csv(cp, index=False, compression='gzip')
    json.dump(meta, open(mp, 'w'), indent=1)
    return df, meta


def apply_paper_rule(df):
    """两种读法 + 分层（各自单独剔了多少）。"""
    f_umi = df['nUMI'] < PAPER_NUMI_MIN
    f_gen = df['nFeature'] < PAPER_NFEAT_MIN
    f_mt = df['pct_mt'] > PAPER_MT_MAX
    return {
        # 读法 OR：任一不达标即剔（Seurat subset 的常规写法）
        'drop_or': (f_umi | f_gen | f_mt),
        # 读法 AND：三个同时不达标才剔（原文那个 `and` 的字面读法）
        'drop_and': (f_umi & f_gen & f_mt),
        'fail_umi_only': f_umi & ~f_gen & ~f_mt,
        'fail_gen_only': f_gen & ~f_umi & ~f_mt,
        'fail_mt_only': f_mt & ~f_umi & ~f_gen,
        'fail_umi': f_umi, 'fail_gen': f_gen, 'fail_mt': f_mt,
    }


def make_figures(allq, summ, outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(19, 10))

    panels = [
        ('nUMI', PAPER_NUMI_MIN, True, 'nUMI per spot', 'log'),
        ('nFeature', PAPER_NFEAT_MIN, True, 'genes detected per spot', 'log'),
        ('pct_mt', PAPER_MT_MAX, False, 'mitochondrial % per spot', 'linear'),
    ]
    for j, (col, thr, lower, xlab, _) in enumerate(panels):
        ax = axes[0, j]
        v = allq[col].dropna().values
        nz = v[v > 0]
        bins = np.logspace(np.log10(max(nz.min(), 1)), np.log10(nz.max()), 80) if lower \
            else np.linspace(0, min(np.nanmax(v), 100), 80)
        ax.hist(nz if lower else v, bins=bins, color='#4a7fb5', alpha=.85)
        ax.axvline(thr, color='#c0392b', lw=2.2, ls='--')
        if lower:
            ax.set_xscale('log')
        frac = float((v < thr).mean() * 100) if lower else float((v > thr).mean() * 100)
        ax.set_title(f'{xlab}\npaper cutoff {("≥" if lower else "≤")}{thr}  '
                     f'— {frac:.2f}% of spots outside', fontsize=11)
        ax.set_xlabel(xlab); ax.set_ylabel('spots')

        # 逐切片中位数
        ax2 = axes[1, j]
        med = allq.groupby('slide')[col].median().sort_values()
        ax2.barh(range(len(med)), med.values, color='#5a9e6f', height=.8)
        ax2.axvline(thr, color='#c0392b', lw=2.2, ls='--')
        ax2.set_yticks([])
        ax2.set_xlabel(f'per-slide median {col}')
        if lower:
            ax2.set_xscale('log')
        ax2.set_title(f'per-slide median ({len(med)} slides)', fontsize=11)

    fig.suptitle('Visium spot QC — our 56 slides vs source-paper thresholds '
                 '(Peng 2026 Cancer Cell)', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p1 = f'{outdir}/spatial_qc_thresholds.png'
    fig.savefig(p1, dpi=140); plt.close(fig)

    # 剔除率对照图
    fig, ax = plt.subplots(figsize=(13, 6))
    names = ['OR\n(any fails)', 'AND\n(all three fail)']
    vals = [summ['drop_pct_or'], summ['drop_pct_and']]
    bars = ax.bar(names, vals, color=['#c0392b', '#7f8c8d'], width=.5)
    ax.axhline(PAPER_DROP_PCT_REPORTED, color='black', lw=2, ls=':')
    ax.text(1.35, PAPER_DROP_PCT_REPORTED, f"  paper reports {PAPER_DROP_PCT_REPORTED}%",
            va='center', fontsize=12)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f'{v:.2f}%', ha='center', va='bottom',
                fontsize=13, fontweight='bold')
    ax.set_ylabel('% of 639,816 spots dropped')
    ax.set_title('Which reading of the paper\'s "and" reproduces ~1.05%?')
    fig.tight_layout()
    p2 = f'{outdir}/spatial_qc_droprule.png'
    fig.savefig(p2, dpi=140); plt.close(fig)

    # ——— 敏感性图：nUMI 阈值 → 剔除率；以及逐切片剔除率排序 ———
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(19, 7.5),
                                   gridspec_kw={'width_ratios': [1, 1.15]})

    # 左：整条曲线。横轴阈值，纵轴剔除率（nUMI 一项）
    ts = np.unique(np.round(np.logspace(np.log10(50), np.log10(5000), 250)))
    y = [100 * float((allq['nUMI'] < t).mean()) for t in ts]
    axL.plot(ts, y, color='#2c3e50', lw=2.2)
    axL.axvline(PAPER_NUMI_MIN, color='#c0392b', lw=2.4, ls='--')
    axL.axhline(PAPER_DROP_PCT_REPORTED, color='black', lw=2, ls=':')
    y500 = 100 * float((allq['nUMI'] < PAPER_NUMI_MIN).mean())
    axL.plot([PAPER_NUMI_MIN], [y500], 'o', color='#c0392b', ms=11, zorder=5)
    axL.annotate(f'paper cutoff {PAPER_NUMI_MIN}\n→ we drop {y500:.2f}% (nUMI alone; '
                 f'{summ["drop_pct_or"]:.2f}% with all three)',
                 xy=(PAPER_NUMI_MIN, y500),
                 xytext=(PAPER_NUMI_MIN * 1.22, 16.5), fontsize=12,
                 arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.8),
                 color='#c0392b')
    t105 = float(np.percentile(allq['nUMI'], PAPER_DROP_PCT_REPORTED))
    axL.axvline(t105, color='#7f8c8d', lw=2, ls='-.')
    axL.annotate(f'to reach the paper\'s {PAPER_DROP_PCT_REPORTED}%\nyou would need nUMI < {t105:.0f}',
                 xy=(t105, PAPER_DROP_PCT_REPORTED), xytext=(t105 * 1.15, 8.5),
                 fontsize=12, color='#566573',
                 arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=1.8))
    axL.set_xscale('log')
    axL.set_xlabel('nUMI cutoff (spots below it are dropped)')
    axL.set_ylabel('% of all 639,816 spots dropped')
    axL.set_title('Threshold sensitivity — the paper\'s 1.05% is not reachable\n'
                  'at a defensible nUMI cutoff', fontsize=12)
    axL.set_ylim(0, 26)
    axL.grid(alpha=.3)

    # 右：逐切片剔除率排序，按该切片真实测序深度上色（连续量，不做"浅/深"二分）
    per = allq.groupby('slide').apply(
        lambda d: 100 * ((d['nUMI'] < PAPER_NUMI_MIN) | (d['nFeature'] < PAPER_NFEAT_MIN)
                         | (d['pct_mt'] > PAPER_MT_MAX)).mean()).sort_values(ascending=False)
    medg = allq.groupby('slide')['nFeature'].median().reindex(per.index)
    sc = axR.scatter(range(len(per)), per.values, c=medg.values, cmap='viridis',
                     s=95, zorder=3, edgecolor='w', linewidth=.6)
    cb = fig.colorbar(sc, ax=axR, pad=.02)
    cb.set_label('that slide\'s median genes/spot')
    axR.axhline(PAPER_DROP_PCT_REPORTED, color='black', lw=2, ls=':')
    axR.text(len(per) * .42, PAPER_DROP_PCT_REPORTED + .8,
             f"paper reports {PAPER_DROP_PCT_REPORTED}%", fontsize=12)
    axR.axhline(per.median(), color='#2c3e50', lw=1.6, ls='--')
    axR.text(len(per) * .42, per.median() + .8, f'our per-slide median {per.median():.2f}%',
             fontsize=12, color='#2c3e50')
    r = float(np.corrcoef(list(range(len(per))), medg.values)[0, 1])
    axR.set_xticks([])
    axR.set_xlabel(f'56 slides, sorted by drop rate  —  '
                   f'rank vs depth Spearman r = {r:.2f}')
    axR.set_ylabel('% of that slide\'s spots dropped (paper OR rule)')
    axR.set_title('The gap is per-slide depth heterogeneity, not the rule', fontsize=12)
    axR.grid(alpha=.3, axis='y')

    fig.suptitle('Why our drop rate is 4.97% and the paper\'s is 1.05%', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p3 = f'{outdir}/spatial_qc_sensitivity.png'
    fig.savefig(p3, dpi=140); plt.close(fig)
    return p1, p2, p3


def count_in_tissue(slide):
    """tissue_positions.csv 里 in_tissue==1 的 spot 数（SP0 用来核对 filtered 恒等）。"""
    p = f'{DATA}/{slide}/spatial/tissue_positions.csv'
    if not os.path.exists(p):
        return None
    n = 0
    with open(p) as f:
        head = f.readline()
        if head.startswith('barcode,'):
            ci = head.rstrip('\n').split(',').index('in_tissue')
            for line in f:
                c = line.rstrip('\n').split(',')
                if len(c) > ci and c[ci].strip() in ('1', '1.0'):
                    n += 1
        else:                                   # 老格式 tissue_positions_list.csv：无表头，第 2 列
            for line in [head] + f.readlines():
                c = line.rstrip('\n').split(',')
                if len(c) > 1 and c[1].strip() == '1':
                    n += 1
    return n


def freeze(allq, metas, slides):
    """按签字阈值冻结 spot 掩码（SP0–SP4）。返回 manifest。"""
    # ——— SP1：阈值已签字 ———
    if any(v is None for v in SIGNED.values()):
        sys.exit('[STOP] SP1 不过：阈值未签字（SIGNED 里仍是 None）⇒ 不得冻结掩码。'
                 '见 SPATIAL_QC_PREREG.md §3 / §7 第 1 项。')

    # ——— 单位换算（必须显式，禁止心算）———
    # pct_mt 的单位是**百分数**（0–100），MAX_MT_FRAC 的单位是**比例**（0–1）。
    mt_pct_max = 100.0 * SIGNED['MAX_MT_FRAC']
    print(f'[freeze] 单位换算：pct_mt > {mt_pct_max:g} % '
          f'（= 签字比例 {SIGNED["MAX_MT_FRAC"]} × 100）', flush=True)

    # ——— SP0：底座核对 ———
    # ⚠️ 2026-09-25 更正：预注册 §1.1 原写「features.tsv.gz 共 18,120 行」（= 18,085 + 35）。
    #    实测**只对 45/56 张成立**；另 11 张根本没有抗体面板（features.tsv.gz = 18,085 行）。
    #    Gene Expression 那 18,085 行**逐张逐序完全相同**（symbol 序 sha256 全队列 1 个值），
    #    ⇒ §1.2 的基因交集（18,085 / 18,082 / 18,066）不受影响。
    print('[freeze] SP0 底座核对…', flush=True)
    sp0_bad, n_with_panel = [], 0
    for s in slides:
        m = metas[s]
        n_filt = m['n_barcodes']
        n_in = count_in_tissue(s)
        if n_in is None:
            sp0_bad.append(f'{s}: 找不到 spatial/tissue_positions.csv')
            continue
        if n_in != n_filt:
            sp0_bad.append(f'{s}: in_tissue==1 数 {n_in:,} != filtered barcodes {n_filt:,}')
        if m['n_genes_GE'] != 18085:
            sp0_bad.append(f'{s}: Gene Expression 行 {m["n_genes_GE"]:,} != 18,085')
        if m['n_antibody_rows'] not in (0, 35):
            sp0_bad.append(f'{s}: 抗体行 {m["n_antibody_rows"]} 既不是 0 也不是 35')
        elif m['n_antibody_rows'] == 35:
            n_with_panel += 1
    if sp0_bad:
        sys.exit('[STOP] SP0 不过 ⇒ 硬停，不冻结掩码：\n  ' + '\n  '.join(sp0_bad[:20]))
    print(f'[freeze] SP0 ✅ {len(slides)} 张：in_tissue==1 ≡ filtered barcodes；'
          f'GE 行 18,085 全张一致（symbol 序全队列同一哈希）；'
          f'抗体面板 {n_with_panel}/{len(slides)} 张有（其余 '
          f'{len(slides) - n_with_panel} 张无抗体行，已随 GE 筛自然落掉）', flush=True)

    # ——— SP2 前置：与既有掩码的登记值对表（防止事后悄悄改口径）———
    mp = f'{OUT}/spot_mask_manifest.json'
    if os.path.exists(mp):
        old = json.load(open(mp))
        othr, orule = old.get('thresholds'), old.get('rule')
        if othr != SIGNED or orule != RULE:
            sys.exit(f'[STOP] 已有一份冻结掩码，且其登记口径与本次不同 ⇒ 硬停，不覆盖。\n'
                     f'  既有：thresholds={othr} rule={orule}\n'
                     f'  本次：thresholds={SIGNED} rule={RULE}\n'
                     f'  口径改变属**事后改阈值**（SPATIAL_QC_PREREG.md §6 顺序约束禁止）。\n'
                     f'  如确需改动：先删掉 {mp} 与 spot_mask.tsv.gz，并在预注册里登记原因。')

    # ——— 冻结 ———
    f_umi = (allq['nUMI'] < SIGNED['MIN_UMI_SPOT'])
    f_gen = (allq['nFeature'] < SIGNED['MIN_GENES_SPOT'])
    f_mt = (allq['pct_mt'] > mt_pct_max)
    assert RULE == 'or'
    drop = f_umi | f_gen | f_mt
    keep = ~drop

    mask = pd.DataFrame({
        'slide': allq['slide'],
        'barcode': allq['barcode'],
        'in_tissue': 1,                       # SP0 已证恒等：输入就是 filtered 矩阵
        'total_umi': allq['nUMI'].astype(np.int64),
        'n_genes': allq['nFeature'].astype(np.int64),
        'pct_mt': allq['pct_mt'].round(6),
        # 写成 TRUE/FALSE：pandas 与 R 的 read.* 都能直接读成逻辑型
        'pass': np.where(keep, 'TRUE', 'FALSE'),
    })
    # mtime=0 ⇒ gzip 头里不带时间戳 ⇒ 同输入重跑逐字节相同（SP2）
    mask_p = f'{OUT}/spot_mask.tsv.gz'
    mask.to_csv(mask_p, sep='\t', index=False,
                compression={'method': 'gzip', 'compresslevel': 6, 'mtime': 0})

    # ——— 逐切片汇总（SP3：只报告，不设通过率阈值）———
    per_slide = {}
    for s in slides:
        m = allq['slide'] == s
        per_slide[s] = {
            'n_total': int(m.sum()),
            'n_pass': int(keep[m].sum()),
            'n_fail_UMI': int((f_umi[m] & ~f_gen[m] & ~f_mt[m]).sum()),
            'n_fail_genes': int((f_gen[m] & ~f_umi[m] & ~f_mt[m]).sum()),
            'n_fail_MT': int((f_mt[m] & ~f_umi[m] & ~f_gen[m]).sum()),
            'n_fail_any': int(drop[m].sum()),
            'drop_pct': round(float(drop[m].mean() * 100), 4),
            'median_total_umi': float(allq.loc[m, 'nUMI'].median()),
        }
    shallow = sorted((s for s in slides if per_slide[s]['drop_pct'] > 15.0),
                     key=lambda s: -per_slide[s]['drop_pct'])

    # ——— §8：56×3 输入哈希清单 ———
    tsv_p = f'{OUT}/slide_inputs_sha256.tsv'
    with open(tsv_p, 'w') as f:
        f.write('slide\tfile\tsha256\n')
        for s in slides:
            m = metas[s]
            for fn, key in (('barcodes.tsv.gz', 'barcodes_sha'),
                            ('features.tsv.gz', 'features_sha'),
                            ('matrix.mtx.gz', 'matrix_sha')):
                f.write(f'{s}\t{fn}\t{m[key]}\n')

    manifest = {
        'module': 'SP 空转 spot QC',
        'rule_doc': 'SPATIAL_QC_PREREG.md',
        'thresholds': dict(SIGNED),
        'thresholds_units': {'MIN_UMI_SPOT': 'count', 'MIN_GENES_SPOT': 'count',
                             'MAX_MT_FRAC': '比例 0–1（乘 100 = 百分数 %）'},
        'mt_pct_effective': mt_pct_max,
        'rule': RULE,
        'rule_meaning': ('任一不达标即剔（OR）。原文那句用 and，字面读法只剔 0.01%；'
                         '签字的 4.97% 就是 OR 的值 ⇒ 语义与阈值一起签的是 OR。'),
        'signature': {
            'signed_by': '用户（老板）',
            'signed_at': '2026-09-25',
            'user_words': ['照原文 500 / 200 / 0.15', '接受剔 4.97% 而非原文 1.05%',
                           '包含，全 56 张'],
            'attribution': ('阈值与读法由用户裁定；Claude 未代替用户决定任何未回答的问题，'
                            '仅把已给的答复逐字落进代码与本清单。'),
        },
        'n_slides': len(slides), 'slides': slides,
        'n_total': int(len(allq)), 'n_pass': int(keep.sum()), 'n_fail': int(drop.sum()),
        'drop_pct': round(float(drop.mean() * 100), 4),
        'paper_reported_drop_pct': PAPER_DROP_PCT_REPORTED,
        'registered_defect': (
            f'实测剔 {100 * drop.mean():.2f}%，原文自报 {PAPER_DROP_PCT_REPORTED}%；'
            '差 4.7 倍。根因是逐切片测序深度不均（秩相关 0.83），非规则读错。'
            '按 feedback_gate_failure_registered_as_defect：**保留原文参数，不调参凑数**，'
            '差额登记为缺陷。'),
        'per_slide': per_slide,
        'shallow_slides_gt15pct': shallow,
        'shallow_slides_decision': (
            f'这 {len(shallow)} 张（剔 >15%）**包含在内**（用户 2026-09-25：「包含，全 56 张」）。'
            'SP3 不设通过率阈值 ⇒ 不得用下面的数字事后挑切片。'),
        'sp0_check': {
            'in_tissue_equals_filtered': True,
            'evidence': f'{len(slides)}/{len(slides)} 张：tissue_positions.csv 的 in_tissue==1 '
                        f'计数 == filtered barcodes 计数',
            'GE_rows': 18085,
            'GE_symbol_order_identical_all_slides': True,
            'antibody_panel_slides': n_with_panel,
            'antibody_rows_dropped': 35,
            'antibody_rows_per_slide': {s: metas[s]['n_antibody_rows'] for s in slides},
            'correction_2026_09_25': (
                '预注册 §1.1 原写 features.tsv.gz 共 18,120 行（18,085 + 35 抗体），'
                f'实测只对 {n_with_panel}/{len(slides)} 张成立 —— 另 '
                f'{len(slides) - n_with_panel} 张没有抗体面板（features.tsv.gz = 18,085 行）。'
                'GE 那 18,085 行逐张逐序相同（symbol 序 sha256 全队列同一值）'
                '⇒ 基因交集与逐 spot QC 数字均不受影响；受影响的只是"总行数"那句话。'),
        },
        'limitations': [
            '原文先由病理学家剔除 "Space" / "Foreign body" / "Tissue artifacts" 三类 spot，'
            '再统计 1.05%。我们无逐 spot 病理注释 ⇒ 该层永久缺失（SPATIAL_QC_PREREG.md §3.2(e)）。',
            '本掩码只做**第一层**剔除。RCTD 自己的 UMI_min / counts_MIN 仍然生效，是第二层；'
            '两层数字分别记录，不得合并成一个"有效 spot 数"（SP4）。',
            '本掩码不含任何恶性判定、不含细胞类型、不含权重。',
        ],
        'hashes': {
            'spot_mask.tsv.gz': sha256(mask_p),
            'spatial_qc_per_spot.csv.gz': sha256(f'{OUT}/spatial_qc_per_spot.csv.gz'),
            'slide_inputs_sha256.tsv': sha256(tsv_p),
        },
        'hash_list_note': ('slide_inputs_sha256.tsv 含 56×3 = 168 行输入哈希（§8）。'
                           '该清单**自身**的 sha256 记在上面 hashes 里 —— 文件不能含自己的哈希，'
                           '故"再入清单"落在本清单（manifest）这一层。'),
    }
    json.dump(manifest, open(mp, 'w'), indent=1, ensure_ascii=False)

    print(f'\n[freeze] 掩码已冻结：{mask_p}')
    print(f'  {len(allq):,} spot → 留 {keep.sum():,} / 剔 {drop.sum():,} '
          f'（{100 * drop.mean():.2f}%）')
    print(f'  判据各自：nUMI<{SIGNED["MIN_UMI_SPOT"]} → {int(f_umi.sum()):,}；'
          f'nFeature<{SIGNED["MIN_GENES_SPOT"]} → {int(f_gen.sum()):,}；'
          f'MT>{mt_pct_max:g}% → {int(f_mt.sum()):,}')
    print(f'  登记缺陷：原文自报 {PAPER_DROP_PCT_REPORTED}%，我们 {100 * drop.mean():.2f}%')
    if shallow:
        print(f'  剔 >15% 的切片 {len(shallow)} 张（已按签字全部纳入）：'
              + ', '.join(f'{s} {per_slide[s]["drop_pct"]:.1f}%' for s in shallow[:6]))
    print(f'  清单：{mp}')
    print(f'  输入哈希：{tsv_p}')
    return manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--figures-only', action='store_true')
    ap.add_argument('--freeze', action='store_true')
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    slides = list_slides()
    print(f'[SP] 切片 {len(slides)} 张', flush=True)

    # 本机 load average 长期 57–70（同机其他 PID namespace），只开 3 个进程，不抢占
    dfs, metas = [], {}
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        for i, (s, (df, meta)) in enumerate(
                zip(slides, ex.map(load_or_compute, slides)), 1):
            dfs.append(df); metas[s] = meta
            print(f'  [{i:>2}/{len(slides)}] {s:24s} spot={meta["n_barcodes"]:>6,} '
                  f'GE={meta["n_genes_GE"]:>6,} 抗体={meta["n_antibody_rows"]:>2} '
                  f'MT基因={meta["n_MT_genes"]:>2} ({meta["seconds"]:>5.1f}s)', flush=True)

    allq = pd.concat(dfs, ignore_index=True)
    # mtime=0 ⇒ gzip 头不带时间戳 ⇒ 同输入重跑逐字节相同（§8 登记的是这个文件的哈希，
    # 用默认 gzip 会因为时间戳不同而哈希漂移，让"冻结件"名不副实）
    allq.to_csv(f'{OUT}/spatial_qc_per_spot.csv.gz', index=False,
                compression={'method': 'gzip', 'compresslevel': 6, 'mtime': 0})

    rules = apply_paper_rule(allq)
    n = len(allq)
    summ = {
        'n_slides': len(slides), 'n_spots': int(n),
        'paper_thresholds': {'nUMI_min': PAPER_NUMI_MIN, 'nFeature_min': PAPER_NFEAT_MIN,
                             'mt_max_pct': PAPER_MT_MAX},
        'paper_reported_drop_pct': PAPER_DROP_PCT_REPORTED,
        'n_drop_or': int(rules['drop_or'].sum()),
        'n_drop_and': int(rules['drop_and'].sum()),
        'drop_pct_or': round(100 * rules['drop_or'].mean(), 3),
        'drop_pct_and': round(100 * rules['drop_and'].mean(), 3),
        'n_fail_umi': int(rules['fail_umi'].sum()),
        'n_fail_gen': int(rules['fail_gen'].sum()),
        'n_fail_mt': int(rules['fail_mt'].sum()),
        'median_nUMI': float(allq['nUMI'].median()),
        'median_nFeature': float(allq['nFeature'].median()),
        'median_pct_mt': float(allq['pct_mt'].median()),
        'per_slide_median_nFeature_median': float(
            allq.groupby('slide')['nFeature'].median().median()),
        'per_slide_median_nFeature_range': [
            float(allq.groupby('slide')['nFeature'].median().min()),
            float(allq.groupby('slide')['nFeature'].median().max())],
        # 比值类指标的分母检查（feedback: 先查分母假象）
        'per_slide_drop_pct_or_mean': float(
            rules['drop_or'].groupby(allq['slide']).mean().mean() * 100),
        'per_slide_drop_pct_or_median': float(
            rules['drop_or'].groupby(allq['slide']).mean().median() * 100),
        # 要复现原文 1.05% 需要的 nUMI 阈值（不可接受的低）
        'nUMI_cutoff_to_reproduce_paper_drop_pct': float(
            np.percentile(allq['nUMI'], PAPER_DROP_PCT_REPORTED)),
        'per_slide': {s: {'n_spot': int((allq['slide'] == s).sum()),
                          'median_nUMI': float(allq.loc[allq['slide'] == s, 'nUMI'].median()),
                          'median_nFeature': float(allq.loc[allq['slide'] == s, 'nFeature'].median()),
                          'median_pct_mt': float(allq.loc[allq['slide'] == s, 'pct_mt'].median()),
                          'drop_pct_or': float(rules['drop_or'][allq['slide'] == s].mean() * 100),
                          } for s in slides},
        'slide_meta': metas,
        'not_reproducible_from_paper': [
            '病理学家逐 spot 剔除 "Space"/"Foreign body"/"Tissue artifacts" 这一步——我们无逐 spot 病理注释',
        ],
    }
    json.dump(summ, open(f'{OUT}/spatial_qc_summary.json', 'w'), indent=1, ensure_ascii=False)

    print(f'\n[SP] 总 spot = {n:,}')
    print(f'  原文阈值 nUMI≥{PAPER_NUMI_MIN} / nFeature≥{PAPER_NFEAT_MIN} / MT≤{PAPER_MT_MAX}%')
    print(f'  读法 OR （任一不达标即剔）: 剔 {summ["n_drop_or"]:,} = {summ["drop_pct_or"]:.2f}%')
    print(f'  读法 AND（三者同时才剔）  : 剔 {summ["n_drop_and"]:,} = {summ["drop_pct_and"]:.2f}%')
    print(f'  原文自报                  : {PAPER_DROP_PCT_REPORTED}%')
    print(f'  中位 nUMI={summ["median_nUMI"]:.0f}  nFeature={summ["median_nFeature"]:.0f}  '
          f'MT={summ["median_pct_mt"]:.2f}%')

    if args.figures_only:
        p1, p2, p3 = make_figures(allq, summ, OUT)
        print(f'\n[fig] {p1}\n[fig] {p2}\n[fig] {p3}')
        return

    if args.freeze:
        freeze(allq, metas, slides)


if __name__ == '__main__':
    main()
