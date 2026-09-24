#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_build_rctd_reference.py —— 为 M5 的 RCTD 构建参考（snRNA，GSE308103）

用途
    把 GP0 的表达矩阵 + 已冻结的细胞标签，拼成 spacexr::Reference() 要的两样东西：
      ① 计数矩阵（基因 × 细胞，稀疏，整数）
      ② 每个细胞的**类型名**（粒度由 caliber 决定，见 RCTD_PREREG.md §2.2）
    再落成一个 HDF5 交给 R 侧读（R 装有 rhdf5，走 HDF5 比 MatrixMarket 文本省几十 GB）。

⚠️ 本脚本**不抽样**。抽样由 R 侧的 Reference(n_max_cells=) 完成（预注册 §2.4），
   以免抽两次。这里导出全部细胞。

⚠️ 本脚本**不判恶性**、**不改任何标签**。标签只做搬运与（在 caliber b/c 下）细化：

   基准谱系取**六谱系重聚类这一版划分**（RCTD_PREREG.md §2.2b-1），不是 gp6 的 A_frozen。
   因为 L2 子群处理时把细胞在谱系间**返还并删去**过，重聚类子集才是更晚的判定；
   六个子集恰好两两互斥、并集 = 全部 413,697 个细胞（已实测）。
   细化时上皮等按 GP8a 的 L2 亚型展开，进不了亚型的细胞**剔除并计数上报，不填补、不猜**。

🔴 未签字不跑：读同目录 calibration.json，缺失或 signed!=true 直接退出。
"""
import os, sys, json, gzip, csv, hashlib, argparse, time
from collections import Counter
import numpy as np
import h5py

ROOT = '/home/eto/luad_v2'
HERE = os.path.dirname(os.path.abspath(__file__))

H5AD     = f'{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad'
LABELS   = f'{ROOT}/results/05_annotation/gp6_cell_labels.csv.gz'
SEEDMAN  = f'{ROOT}/results/05_annotation/seed_representativeness_manifest.json'
INTEG    = f'{ROOT}/results/04_integration/seurat_trad'
ANN      = f'{ROOT}/results/05_annotation'

# 预注册 §8 冻结哈希
FROZEN = {
    H5AD:   'a276cd1af69a4620ff9e3d64f4b93f33461e646537c5f2c9da1aff8c9150b9de',
    LABELS: 'fa031a52a4ab3b7098d3fca5555b0749b2f4622711153c0b742f8ae09aaf12e2',
}

# 六谱系的「划分清单」与「污染清单」：(subset, nocontam) 两个 sha256。
# 这两组清单决定**谁进参考、谁算哪个谱系**，任一口径都读 ⇒ 冻结它们是防**事后换清单**
# （换清单 = 换口径，法则 3.1）。上皮/成纤维两行两哈希相同，与它们 excluded 列为 0 一致。
PARTITION_SHA = {
    'epiA':     ('d12a115a87d0790d7fe8cf27b16170999d94daba51f3016fdb79c5aa3fdde55e',
                 'd12a115a87d0790d7fe8cf27b16170999d94daba51f3016fdb79c5aa3fdde55e'),
    'fibroA':   ('505b82a900745cb443dd31e6ab8fdd1c544cf80d7b126b22543c03ec3606f516',
                 '505b82a900745cb443dd31e6ab8fdd1c544cf80d7b126b22543c03ec3606f516'),
    'tnkA':     ('8cb87515d5463bdbdb9c6885decef19611a73b5a8f6abbdc850507d5c620b0c5',
                 'b27ae33fa965d98e4d116c9ebc6387b621a25f8f28fa5833f9d9aeb05ff6454a'),
    'myeloidA': ('2fd34840c92f8a6ebf11ea8713903af7ecd7d0ed5a6622864836f01805b2fb1f',
                 'e93166da0922af837185507fef54a38e06a8575169dff546ca7a41e6a7b4190e'),
    'endoA':    ('9ee8d0c08a26bd14592b654d3c9578561d1143f5bb6f2a74e544e8306f2396ac',
                 'd9b5bc57254ac2b6fecea7ab0eec93dde58a98f484648514d8f5035e3dd9a786'),
    'bplasmaA': ('7e34c5c5baedb0ecfef071a8e2f3c738964c97cc29d74aa13c65c2a731a5cd0a',
                 '16e0607fb8344bf6d4ee97a5119a0a0644e8109c0535dea2c16513ba6f427c53'),
}

# 注释表哈希（决定「每个簇叫什么名字」，即粒度口径本身）。2026-09-25 登记。
# ⚠️ 两套面板并存，**不是**版本更替：`_s4_` 是 2026-09-22 裁定的新面板（髓系 8→12、
#    内皮 2→8 ⇒ 细版 29 → 39 型），旧表仍是 caliber c 的口径，两套都留着。
#    S4 与旧面板共用同一套聚类（簇 id 逐个相同、excluded 相同），差别只在 argmax ⇒
#    读错表**不会**报错、只会静默给出另一个合法的型数 ⇒ 靠 expect_n_types 自检兜住。
ANNOTATION_SHA = {
    'epiA_cluster_annotation.csv':       '589b69cc03c89afc991d18baaa19fecc935325d11b58f2dccddc2a55ab7c8e4e',
    'fibroA_cluster_annotation.csv':     'c9dff58d34b16d762f673e417431babde30d7a2f34c5ed2f97bb5941640e8f0a',
    'tnkA_cluster_annotation.csv':       'd0a9f377f0cd973a888b961b5510685cadc629402d37bfe1f953532806cff412',
    'myeloidA_cluster_annotation.csv':   'bc233579ffc36ff946c343a9bce6b26142ad63800952035ef96ffbdd2ca21678',
    'myeloidA_s4_cluster_annotation.csv': '0057f39d1fea98bf7f151ebec9fe1bdc824bd1687c9f9f42247a6401a3cbc186',
    'endoA_cluster_annotation.csv':      '673b3323b00ce5225c8f07b2467e9f5cb26e5bc497772707889815dbb1b10a34',
    'endoA_s4_cluster_annotation.csv':   '00f9fae5485e3a74d57acb894d058e3f5a3329404342b29a9263bb9fe1bc137d',
    'bplasmaA_cluster_annotation.csv':   'ad6e9a7120c328a33a55b43a4d3ee5bf0fb64854fd4781387ecbb99efe1bc094',
}


def frozen_items(cal=None):
    """全部要硬断言的输入 → {绝对路径: 期望 sha256}。

    `cal` 给定时**只**追加该口径真正会读的注释表（caliber d 读 `_s4_` 那张）。
    为什么要按口径裁：注释表决定「进参考的细胞叫什么名字」，是粒度口径本身；
    但它同时也是最容易被无关的重跑顺手改掉的产物。只冻结**这次真的会读**的表，
    才能做到「改了一张与本次口径无关的表，不该拦下这一跑」。
    """
    d = dict(FROZEN)
    for tag, (h_sub, h_noc) in PARTITION_SHA.items():
        d[f'{ANN}/{tag}_subset_barcodes.txt'] = h_sub
        d[f'{ANN}/{tag}_nocontam_subset_barcodes.txt'] = h_noc
    if cal is not None:
        panels = cal.get('panels', {}) or {}
        refine = {'a': [], 'b': ['epiA'], 'c': list(LINEAGE_TAG.values()),
                  'd': list(LINEAGE_TAG.values())}[cal['caliber']]
        for tag in refine:
            suf = '_s4' if panels.get(tag) == 's4' else ''
            p = f'{tag}{suf}_cluster_annotation.csv'
            if p not in ANNOTATION_SHA:
                die(f'注释表 {p} 没有登记哈希 ⇒ 先把它写进 ANNOTATION_SHA 再跑')
            d[f'{ANN}/{p}'] = ANNOTATION_SHA[p]
    return d

# 六谱系在 A_frozen 里的名字 → 整合产物的目录 tag
LINEAGE_TAG = {
    '上皮':   'epiA',
    '成纤维': 'fibroA',
    'T/NK':   'tnkA',
    '髓系':   'myeloidA',
    '内皮':   'endoA',
    'B/浆':   'bplasmaA',
}


def die(msg):
    print(f'\n🔴 硬停：{msg}\n', flush=True)
    sys.exit(1)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_calibration():
    p = f'{HERE}/calibration.json'
    if not os.path.exists(p):
        die(f'缺少签字文件 {p}。\n'
            f'    本臂的口径（参考粒度 / 跑哪些切片 / WHO 映射）必须先在 '
            f'RCTD_PREREG.md §9 裁定，\n'
            f'    裁定结果写进 calibration.json 后本脚本才跑。')
    with open(p, encoding='utf-8') as fh:
        c = json.load(fh)
    if not c.get('signed'):
        die(f'{p} 的 signed 不为 true ⇒ 口径未签字。')
    for k in ('caliber', 'slides', 'n_max_cells', 'who_mapping'):
        if k not in c:
            die(f'签字文件缺字段 `{k}`。')
    if c['caliber'] not in ('a', 'b', 'c', 'd'):
        die(f"caliber 只能是 a/b/c/d，收到 {c['caliber']!r}")

    # —— caliber d = S4 面板的细版（39 型）——
    # 🔴 为什么不是 caliber c：2026-09-22 裁定的 S4 面板把髓系 8→12 型、内皮 2→8 型，
    #    细版因此是 **39** 型而不是 29 型。S4 与旧面板**共用同一套簇**（簇 id 逐个相同、
    #    excluded 集合相同、只有 argmax 变了），所以只是换一张注释表 —— 但也正因为
    #    「只差一个文件名」，最容易发生的事故就是**静默沿用旧表**、跑出 29 型还以为是 39。
    #    故：panels 必须在签字文件里显式写出（值只认 's4'）。
    panels = c.get('panels', {})
    if not isinstance(panels, dict):
        die('panels 字段必须是 {谱系标签: "s4"} 形式的字典')
    bad = {k: v for k, v in panels.items() if v != 's4'}
    if bad:
        die(f"panels 的值目前只认识 's4'，收到 {bad}")
    unknown = set(panels) - set(LINEAGE_TAG.values())
    if unknown:
        die(f'panels 里的谱系标签不存在：{sorted(unknown)}（可选：{sorted(LINEAGE_TAG.values())}）')
    # 🔴 只守一个方向：**细版不许静默退回旧表**（caliber d 必须有 S4 面板）。
    #    反方向（带 panels 但跑粗版）不设闸：caliber a/b/c 不展开任何谱系，
    #    panels 根本不参与，拦它只会挡住用户明确要的粗版那一次。
    #    `expect_n_types` 是这条闸的第二道保险（39 vs 29 都是合法型数，只有它抓得住）。
    if c['caliber'] == 'd' and not panels:
        die(f"caliber='d'（S4 面板的细版）要求 panels 非空，但签字文件里没有 panels。\n"
            f"    （要出旧 L2 面板的细版请用 --caliber c；两版口径不同，不得混用。）")
    return c


def verify_frozen(cal=None):
    items = frozen_items(cal)
    print(f'校验冻结哈希（{len(items)} 个文件；h5ad 3 GB 约 15 秒）…', flush=True)
    for p, want in items.items():
        got = sha256(p)
        if got != want:
            die(f'哈希不符：{p}\n  期望 {want}\n  实得 {got}')
        print(f'  ✅ {os.path.relpath(p, ROOT)}', flush=True)


def load_partition():
    """六个重聚类子集清单 → dict(条码 -> 谱系名)。

    🔴 基准谱系用**重聚类这一版**，不用 gp6 的 A_frozen：L2 子群处理把细胞在谱系间
       返还并删去过（用户 2026-09-24 确认），因此六个子集是更晚、更全的判定。
       实测：上表六个数相加 = 并集 = 413,697 = 全部细胞，且两两互斥 ⇒ 是一个划分。
       子集不互斥说明清单串了，当场停。
    """
    lin_of = {}
    for lin, tag in LINEAGE_TAG.items():
        p = f'{ANN}/{tag}_subset_barcodes.txt'
        if not os.path.exists(p):
            die(f'缺少重聚类子集清单 {p}')
        n0 = len(lin_of)
        with open(p, encoding='utf-8') as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                if s in lin_of:
                    die(f'条码 {s} 同时落在 {lin_of[s]} 与 {lin} 的子集里 ⇒ 子集不互斥')
                lin_of[s] = lin
        print(f'  {lin:6s} {tag:10s} {len(lin_of) - n0:>7,}', flush=True)
    return lin_of


def compare_to_afrozen(lin_of):
    """与 gp6 的 A_frozen 交叉核对：**只上报差异方向，不判对错、不报错**。

    差异是预期的（返还/删去所致）；这里记下来是为了让 manifest 里留有痕迹，
    免得以后有人看到两版对不上以为是 bug。
    """
    a = {}
    with gzip.open(LABELS, 'rt') as fh:
        for row in csv.DictReader(fh):
            a[row['cell_barcode']] = row['A_frozen']
    if set(a) != set(lin_of):
        die(f'两版标签的条码集合不同：A_frozen {len(a):,} vs 重聚类 {len(lin_of):,}')
    flips = Counter((a[b], lin_of[b]) for b in lin_of if a[b] != lin_of[b])
    n = sum(flips.values())
    print(f'\n  与 A_frozen 不一致 {n:,} 个细胞（{n / len(lin_of) * 100:.1f}%，'
          f'预期：子群处理时返还/删去所致）：', flush=True)
    for (x, y), k in flips.most_common():
        print(f'     {x:8s} → {y:8s} {k:>6,}', flush=True)
    return n, {f'{x}→{y}': k for (x, y), k in flips.items()}


def signed_seed_and_r():
    """读 §12.1 判据产物：每个谱系的最终注释种子与 r*。
    规则（GP5_report.md §12.1）：gap ≤ 0.01 视为并列 ⇒ 沿用种子 0；否则用 best_seed。"""
    with open(SEEDMAN, encoding='utf-8') as fh:
        m = json.load(fh)
    per = m['verdict']['per_object']
    out = {}
    for tag in LINEAGE_TAG.values():
        if tag not in per:
            die(f'种子判据产物里没有 {tag}（{SEEDMAN}）')
        e = per[tag]
        # r* 以已签字的 <tag>_rstar.json 为准（判据产物里的值应与之一致，不一致即停）
        with open(f'{ANN}/{tag}_rstar.json', encoding='utf-8') as fh:
            rs = json.load(fh)
        if abs(rs['r_star'] - e['r_star']) > 1e-9:
            die(f'{tag} 的 r* 在两份产物里不一致：rstar.json={rs["r_star"]} vs '
                f'seed_manifest={e["r_star"]}')
        seed = 0 if e['pass_'] else int(e['best_seed'])
        out[tag] = dict(r_star=rs['r_star'], seed=seed, gap=e['gap'], pass_=e['pass_'])
    return out


def load_subset_and_nocontam(tag):
    """该谱系的重聚类子集清单与「非污染」子清单 → (subset, nocontam)。

    污染 = 被判**跨谱系污染**并移出的簇（`*_cluster_annotation.csv` 的 excluded 列，
    理由逐条写明在 excluded_reason 里）。这些簇的**谱系归属本身就不可信**，
    不只是亚型分不出 ⇒ 粗版、细版**都**要剔除它们（预注册 §2.2b-2）。
    """
    sp = f'{ANN}/{tag}_subset_barcodes.txt'
    np_ = f'{ANN}/{tag}_nocontam_subset_barcodes.txt'
    for p in (sp, np_):
        if not os.path.exists(p):
            die(f'缺少污染判据产物 {p}')

    def _rd(p):
        s = set()
        with open(p, encoding='utf-8') as fh:
            for line in fh:
                t = line.strip()
                if t:
                    s.add(t)
        return s

    sub, noc = _rd(sp), _rd(np_)
    if not noc <= sub:
        die(f'{tag} 的 nocontam 清单不完全是子集的子集（多出 {len(noc - sub)} 个）')
    return sub, noc


def load_subtype_map(tag, r_star, seed, nocontam, ann_suffix=''):
    """L2 亚型：clusters.csv.gz 的 harmony_res{r*}_seed{seed} 列 → 注释表的 argmax。

    🔴 **不能只用 argmax**：被判跨谱系污染的簇**仍有 argmax 值**，照 argmax 映射会把
       已被剔除的细胞当亚型喂进参考。故以 nocontam 清单为**唯一准入**，
       并用注释表 excluded 列反向核对（同一决定两处记录不一致即停）。

    `ann_suffix`：`'_s4'` 时读 `{tag}_s4_cluster_annotation.csv`（S4 面板，caliber d）。
    只换注释表 —— 簇本身不变（S4 与旧面板共用同一套和谐聚类与簇 id），
    所以 `clusters.csv.gz` 与 `harmony_res{r*}_seed{seed}` 列在两种面板下是同一列。

    返回 dict(条码 -> 亚型)。只含 nocontam 里的细胞；落在无 argmax 的簇里 → 值为 None。
    """
    col = f'harmony_res{r_star}_seed{seed}'
    cp = f'{INTEG}/{tag}/clusters.csv.gz'
    ap = f'{ANN}/{tag}{ann_suffix}_cluster_annotation.csv'
    for p in (cp, ap):
        if not os.path.exists(p):
            die(f'缺少 L2 产物 {p}（细版需要；若只出粗版请把 caliber 设为 a）')

    c2s, excluded_cl = {}, set()
    with open(ap, encoding='utf-8') as fh:
        r = csv.DictReader(fh)
        if 'argmax' not in r.fieldnames:
            die(f'{ap} 没有 argmax 列')
        for row in r:
            c2s[row['cluster']] = row['argmax']
            if str(row.get('excluded', '')).lower() in ('true', '1'):
                excluded_cl.add(row['cluster'])

    out = {}
    with gzip.open(cp, 'rt') as fh:
        r = csv.DictReader(fh)
        if col not in r.fieldnames:
            die(f'{cp} 没有列 {col}（可选列：{r.fieldnames[4:]}）')
        for row in r:
            b = row['cell_barcode']
            if b not in nocontam:
                continue
            cl = row[col]
            if cl in excluded_cl:
                die(f'{tag} 簇 {cl} 在注释表里标 excluded，但 nocontam 清单里还有它的'
                    f'细胞（如 {b}）⇒ 两份产物矛盾，停。')
            out[b] = c2s.get(cl)

    if len(out) != len(nocontam):
        die(f'{tag} nocontam 清单 {len(nocontam):,} 个，但只在 clusters 表里找到 '
            f'{len(out):,} 个 ⇒ 有细胞不在本次子聚类结果里')
    n_none = sum(1 for v in out.values() if v is None)
    if n_none:
        print(f'     ⚠️ {tag} 有 {n_none:,} 个细胞落在没有 argmax 的簇里 ⇒ 按未映射剔除',
              flush=True)
    return out


def build_cell_types(caliber, barcodes, lin_of, contaminated, subtypes):
    """返回 (cell_types:list[str|None], 报告 dict)。

    **粗、细两版共用同一个细胞宇宙**：六个重聚类子集的并集 **减去**跨谱系污染簇。
    这样两版的差别**只剩粒度**（颗粒度对照才干净），按预注册 §2.2 出两版。

      caliber a → 6 谱系名（粗）
      caliber b → 上皮按 L2 展开，其余五谱系用谱系名（11 型）
      caliber c → 六谱系都用**旧 L2 面板**展开（29 型）
      caliber d → 六谱系都用 **S4 面板**展开（39 型；髓系 12、内皮 8）

    剔除分两种理由，分开计数（混在一起会看不出是哪一层出的问题）：
      · `contaminated` —— 落在被判跨谱系污染的簇里（§2.2b-2）
      · `unmapped`     —— 非污染，但所在簇在注释表里没有 argmax
    """
    refine = {'a': [], 'b': ['epiA'], 'c': list(LINEAGE_TAG.values()),
              'd': list(LINEAGE_TAG.values())}[caliber]
    types, drop = [], {'contaminated': 0, 'unmapped': 0}
    n_none = 0
    for b in barcodes:
        lin = lin_of.get(b)
        if lin is None:
            die(f'细胞 {b} 不在六个重聚类子集的任何一个里 ⇒ 划分不完整，停')
        if b in contaminated:
            drop['contaminated'] += 1
            types.append(None)
            continue
        tag = LINEAGE_TAG[lin]
        if tag in refine:
            st = subtypes[tag].get(b)
            if st is None:
                drop['unmapped'] += 1
                types.append(None)
                n_none += 1
            else:
                types.append(st)
        else:
            types.append(lin)               # 不细化的谱系，保留谱系名
    if n_none:
        print(f'\n⚠️ 有 {n_none:,} 个非污染细胞所在簇没有 argmax ⇒ 剔除', flush=True)
    return types, drop


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=None,
                    help='默认 = results/08_spatial_deconv/reference_<caliber>，'
                         '保证粗/细两档各自落盘、互不覆盖')
    ap.add_argument('--caliber', choices=['a', 'b', 'c', 'd'],
                    help='覆盖 calibration.json 里的 caliber（粗/细两版要各跑一次：'
                         '粗=a 六谱系名，细=d 用 S4 面板的 39 型）')
    args = ap.parse_args()
    if args.out is None:
        # 两档 caliber 必须写不同文件，否则后跑的把先跑的覆盖掉（且看不出）
        tag = args.caliber or 'from-calibration'
        args.out = f'{ROOT}/results/08_spatial_deconv/reference_{tag}'

    t0 = time.time()
    cal = load_calibration()
    if args.caliber:
        cal['caliber'] = args.caliber
    # --caliber 覆盖后**立刻**重新核一遍：命令行覆盖是常规操作（粗/细两版各跑一次），
    # 也正是最容易把签字口径跑串的地方。放在读 h5ad（3 GB、约 30 秒）**之前**。
    panels = cal.get('panels', {}) or {}
    if cal['caliber'] == 'd' and not panels:
        die(f"命令行把 caliber 覆盖成 'd'，但签字文件里没有 panels ⇒ 无面板可换，停。")
    if panels and cal['caliber'] != 'd':
        print(f'ℹ️ 提醒：本次 caliber={cal["caliber"]!r} 不展开谱系，签字文件里的 '
              f'panels={panels} 本次**不参与**（panels 只对 caliber d 生效）。', flush=True)
    print(f'✅ 口径已签字：caliber={cal["caliber"]}  slides={cal["slides"]}  '
          f'n_max_cells={cal["n_max_cells"]}  who_mapping={cal["who_mapping"]}  '
          f'panels={panels or "空"}', flush=True)

    verify_frozen(cal)

    print('\n读 h5ad（413,697 × 18,069，约 5 GB 内存）…', flush=True)
    with h5py.File(H5AD, 'r') as f:
        data    = f['X/data'][()].astype(np.int32, copy=False)
        indices = f['X/indices'][()].astype(np.int32, copy=False)
        indptr  = f['X/indptr'][()].astype(np.int64, copy=False)
        ncell, ngene = int(f['X/indptr'].shape[0] - 1), int(f['var/_index'].shape[0])
        genes    = [x.decode() if isinstance(x, bytes) else str(x) for x in f['var/_index'][:]]
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in f['obs/_index'][:]]
    print(f'  X: {ncell} 细胞 × {ngene} 基因, nnz={len(data):,}', flush=True)

    if data.size and int(data.max()) > 2**31 - 1:
        die('计数超出 int32')

    print('\n读六个重聚类子集清单（= 基准谱系层）…', flush=True)
    lin_of = load_partition()
    print(f'  合计 {len(lin_of):,} 个细胞', flush=True)
    miss = [b for b in barcodes if b not in lin_of]
    extra = [b for b in lin_of if b not in set(barcodes)]
    if miss or extra:
        die(f'划分与 h5ad 条码不一致：h5ad 缺 {len(miss)} 个（例 {miss[:3]}），'
            f'清单多 {len(extra)} 个（例 {extra[:3]}）')
    print(f'  ✅ 与 h5ad 的 {ncell:,} 个条码完全对上（两两互斥、无遗漏）', flush=True)
    n_flip, flip_dir = compare_to_afrozen(lin_of)

    print('\n读污染清单（粗、细两版都剔除这些细胞）…', flush=True)
    contaminated, nocontam = set(), {}
    for lin, tag in LINEAGE_TAG.items():
        sub, noc = load_subset_and_nocontam(tag)
        nocontam[tag] = noc
        contaminated |= (sub - noc)
        print(f'  {tag:10s} 子集 {len(sub):>7,}  非污染 {len(noc):>7,}  '
              f'污染 {len(sub - noc):>6,}', flush=True)
    print(f'  污染合计 {len(contaminated):,}（六子集互斥，故可直接相加）', flush=True)

    seeds = signed_seed_and_r()
    print('\n各谱系签字的 r* 与注释种子（GP5_report.md §12.1）：', flush=True)
    for tag, s in seeds.items():
        print(f'  {tag:10s} r*={s["r_star"]}  seed={s["seed"]}  '
              f'gap={s["gap"]:.6f}  {"并列→种子0" if s["pass_"] else "换种子"}', flush=True)

    refine = {'a': [], 'b': ['epiA'], 'c': list(LINEAGE_TAG.values()),
              'd': list(LINEAGE_TAG.values())}[cal['caliber']]
    subtypes = {}
    if refine:
        print('\n读 L2 亚型…', flush=True)
    for tag in refine:
        s = seeds[tag]
        suf = '_s4' if panels.get(tag) == 's4' else ''
        m = load_subtype_map(tag, s['r_star'], s['seed'], nocontam[tag], suf)
        subtypes[tag] = m
        seen = sorted({v for v in m.values() if v})
        print(f'  {tag:10s}{suf or "":4s} harmony_res{s["r_star"]}_seed{s["seed"]} → '
              f'{len(seen)} 型：{", ".join(seen)}', flush=True)

    types, drop = build_cell_types(cal['caliber'], barcodes, lin_of, contaminated, subtypes)
    print(f'\n剔除（不填补、不猜）：污染 {drop["contaminated"]:,}，'
          f'无 argmax {drop["unmapped"]:,}', flush=True)

    keep = np.array([t is not None for t in types], dtype=bool)
    n_keep = int(keep.sum())
    cell_types = [t for t, k in zip(types, keep) if k]
    if n_keep == 0:
        die('细化后一个细胞都不剩')

    # —— G1：维度与类型名自检 ——
    if len(set(cell_types)) < 2:
        die(f'参考只有 {len(set(cell_types))} 种类型，RCTD 至少要 2 种')
    for t in set(cell_types):
        if t is None or str(t).startswith('UNMAPPED'):
            die(f'类型名里出现 {t!r}（说明有未映射的细胞漏过了筛选）')

    # —— G2：粒度自检（防「静默沿用旧注释表」）——
    # caliber d 若把 S4 表读成旧 L2 表，型数会**静默**变成 29 而不是 39（两者都合法、
    # 都能跑完、都不报错）。型数是这个实验唯一可自动核对的粒度指纹 ⇒ 写死在签字文件里。
    want = cal.get('expect_n_types', {}).get(cal['caliber'])
    if want is not None and len(set(cell_types)) != want:
        die(f"caliber={cal['caliber']} 应为 {want} 型，实得 {len(set(cell_types))} 型："
            f"{sorted(set(cell_types))}\n"
            f"    ⇒ 注释表读错了（caliber d 要求 <谱系>_s4_cluster_annotation.csv）。")
    if want is not None:
        print(f'  ✅ 粒度自检：{len(set(cell_types))} 型 = 签字文件登记的 {want} 型', flush=True)

    cnt = Counter(cell_types)
    print(f'\n参考细胞类型（{len(cnt)} 型，{n_keep:,} 细胞）：', flush=True)
    for t, n in cnt.most_common():
        print(f'  {t:22s} {n:>8,}', flush=True)

    # —— 按保留的细胞切稀疏矩阵（CSC：基因 × 细胞）——
    # h5ad 的 X 是 (细胞 × 基因) CSR ⇒ indptr 按**细胞**分段，段内是基因号。
    # 所以「按细胞分段、段内按基因号」正是 (基因 × 细胞) 的 CSC，重排只是删段。
    if n_keep != ncell:
        print(f'  按保留细胞重排稀疏矩阵（{ncell} → {n_keep} 列）…', flush=True)
        drop_idx = np.nonzero(~keep)[0]
        m = np.ones(len(data), dtype=bool)
        for c in drop_idx:                      # 丢掉的细胞很少，反向删段最快
            m[indptr[c]:indptr[c + 1]] = False
        data    = data[m]
        indices = indices[m]
        indptr  = np.concatenate(([0], np.cumsum(np.diff(indptr)[keep])))
        del m, drop_idx

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    hp = f'{args.out}.h5'
    print(f'\n写 {hp} …', flush=True)
    with h5py.File(hp, 'w') as f:
        # 压缩用 level 1：indices 是递增的，压得动；级别高上去只是白等时间
        g = f.create_group('counts')
        g.create_dataset('data',    data=data,    compression='gzip', compression_opts=1)
        g.create_dataset('indices', data=indices, compression='gzip', compression_opts=1)
        g.create_dataset('indptr',  data=indptr,  compression='gzip', compression_opts=1)
        g.attrs['shape'] = np.array([ngene, n_keep], dtype=np.int64)   # [基因, 细胞]
        f.create_dataset('genes', data=np.array(genes, dtype=object),
                         dtype=h5py.special_dtype(vlen=str))
        f.create_dataset('cell_types', data=np.array(cell_types, dtype=object),
                         dtype=h5py.special_dtype(vlen=str))
        f.attrs['caliber']   = cal['caliber']
        f.attrs['built_at']  = time.strftime('%Y-%m-%dT%H:%M:%S')
        f.attrs['source']    = os.path.relpath(H5AD, ROOT)
        f.attrs['labels']    = os.path.relpath(LABELS, ROOT)
        f.attrs['note']      = ('基准谱系取六谱系重聚类子集（两两互斥、并集=全部细胞），'
                                '不用 A_frozen；跨谱系污染簇已剔除。'
                                '抽样不在此处发生，由 R 侧 Reference(n_max_cells=) 完成。'
                                '本文件不含任何恶性判定。')

    man = dict(
        script=os.path.relpath(os.path.abspath(__file__), ROOT),
        built_at=time.strftime('%Y-%m-%dT%H:%M:%S'),
        caliber=cal['caliber'], slides=cal['slides'], who_mapping=cal['who_mapping'],
        signed_by=cal.get('signed_by'), signed_at=cal.get('signed_at'),
        n_genes=ngene, n_cells_h5ad=ncell, n_cells_reference=n_keep,
        nnz=int(len(data)),
        cell_type_counts=dict(cnt),
        seeds={k: v for k, v in seeds.items()},
        lineage_layer='六谱系重聚类子集（两两互斥，并集 = 全部细胞）',
        n_flip_vs_afrozen=n_flip, flip_directions=flip_dir,
        n_contaminated_excluded=drop['contaminated'],
        n_unmapped_excluded=drop['unmapped'],
        n_cell_universe=n_keep,
        frozen={os.path.relpath(p, ROOT): h for p, h in frozen_items(cal).items()},
        output=os.path.relpath(hp, ROOT),
        elapsed_sec=round(time.time() - t0, 1),
    )
    with open(f'{args.out}.manifest.json', 'w', encoding='utf-8') as fh:
        json.dump(man, fh, ensure_ascii=False, indent=1)
    print(f'[out] {hp}\n[out] {args.out}.manifest.json\n'
          f'耗时 {man["elapsed_sec"]} 秒', flush=True)


if __name__ == '__main__':
    main()
