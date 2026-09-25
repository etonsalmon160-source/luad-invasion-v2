#!/usr/bin/env python3
"""分辨率阶梯（PREREG_v2 六臂）· 判读。

输入：results/07_he_pathology/resolution_ladder/per_slide/<臂>/<切片>.csv.gz
      results/07_he_pathology/resolution_ladder/embeds/<臂>/<切片>.npy
输出：ladder_stats.json（全部数值）、ladder_by_arm.csv、ladder_reading.md、
      ladder_figure.png（图内文字全英文）

对应 PREREG_v2 §8 的四项指标 + §8.5 判读表：
  1. 切片级 mean score vs 分期 Spearman rho（n=56），两套提示词各一份
  2. 扣掉密度后的偏相关（控该臂自身的 <img, d>），置换 20,000 次给经验 p
  3. 管线锚 G1：L448P0 的 rho 必须 ∈ [0.60, 0.78]，越界**硬停**
  4. 同一 spot 上各臂 argmax 类目相对 L448P0 的改变率
  另附 G2 彩色守卫的逐臂复核（chroma_std > 5.0）

🔴 本脚本**只**把预注册里写好的量算出来，不新增任何判据。
   预注册没写死到字面的三处操作化，在此**显式登记**（见 REG里STRATION）：
   · d 的定义逐字沿用 09_density_free_ceiling.py（高/低 tissue_frac 各 25% 的均值差）
   · tissue_frac 取 L448 那一列（§5：过滤是 spot 的属性、不是臂的属性）⇒ 六臂同一分组
   · 偏相关在**切片级**算（与指标 1 同层），用秩残差法：rank 后各自对 rank(密度) 回归取残差再相关
   · 置换对象 = 切片的分期标签；种子 0；双侧经验 p
"""
import os
import sys
import json
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon

ROOT = '/home/eto/luad_v2'
LADDER = f'{ROOT}/results/07_he_pathology/resolution_ladder'
PER_SLIDE = f'{LADDER}/per_slide'
EMB = f'{LADDER}/embeds'
OUT = f'{ROOT}/results/07_he_pathology/resolution_sweep'

ARMS = ['L448P0', 'L448P1', 'L320P0', 'L320P1', 'L224P0', 'L224P1']   # 顺序 = 预注册 §4
ANCHOR = 'L448P0'
STAGE_NUM = {'Normal': 0, 'AAH': 1, 'AIS': 2, 'MIA': 3, 'LUAD': 4}
CROP_OF = {a: int(a[1:4]) for a in ARMS}
PREP_OF = {a: a[4:] for a in ARMS}

N_SLIDE_EXPECT = 56            # 全队列；不等于这个数就是没跑完
N_PERM = 20000                 # 预注册 §8 第 2 项
RNG_SEED = 0
G1_LO, G1_HI = 0.60, 0.78      # 预注册 §4 / §10
G2_MIN_CHROMA = 5.0            # 预注册 §7 / §10
FLAT_BAND = 0.05               # 「实质相同」的带宽（与 18_arms_stats 的复现判据同量级）
NEO_IDX = [0, 2, 3]            # A 套：腺癌 + AAH + AIS（与 03/17 同源）
LESION_IDX = [0, 2, 3]         # B 套：浸润 + 非典型增生 + 贴壁


def die(msg):
    print(f'\n[硬停] {msg}', flush=True)
    sys.exit(2)


def load_arm(arm):
    """一只臂的全部切片：逐 spot 表（对齐好的 spot 集）+ 图像向量。

    返回 (df, E)：df 行序与 E 行序**逐个对应**（19 号脚本同一个循环里配对写出）。
    """
    d = f'{PER_SLIDE}/{arm}'
    if not os.path.isdir(d):
        die(f'缺目录 {d} —— 阶梯还没跑过这只臂。')
    slides = sorted(f[:-7] for f in os.listdir(d) if f.endswith('.csv.gz'))
    if len(slides) != N_SLIDE_EXPECT:
        die(f'{arm} 只有 {len(slides)} 张切片（应为 {N_SLIDE_EXPECT}）⇒ 阶梯没跑完，'
            f'不判读。跑完再跑本脚本（19 号脚本会跳过已完成切片）。')
    frames, embs = [], []
    for s in slides:
        missing = f'{EMB}/{arm}/{s}.npy'
        if not os.path.exists(missing):
            die(f'{arm} 的 {s} 有逐 spot 表但缺图像向量 {missing} ⇒ 不判读。')
        with gzip.open(f'{d}/{s}.csv.gz', 'rt') as fh:
            frames.append(pd.read_csv(fh))
        embs.append(np.load(missing))
    df = pd.concat(frames, ignore_index=True)
    E = np.vstack(embs)
    if len(df) != len(E):
        die(f'{arm}：逐 spot 表 {len(df)} 行 vs 图像向量 {len(E)} 行，对不上 ⇒ 不判读。')
    return df, E


def argmax_cat(df, idx):
    """逐 spot 的 argmax 类目（在给定的标签子集或全标签上取）。"""
    return df[[f'a{j}' for j in idx]].values.argmax(1) if len(idx) == 7 \
        else df[[f'a{j}' for j in idx]].values.argmax(1)


def density_direction(E, tf):
    """密度方向 d —— 逐字沿用 09_density_free_ceiling.py 的定义：
    组织最多的 25% 与最少的 25% 的图像向量均值之差，单位化。"""
    lo, hi = np.quantile(tf, .25), np.quantile(tf, .75)
    n_hi, n_lo = int((tf >= hi).sum()), int((tf <= lo).sum())
    if n_hi < 10 or n_lo < 10:
        return None, n_hi, n_lo
    d = E[tf >= hi].mean(0) - E[tf <= lo].mean(0)
    return d / np.linalg.norm(d), n_hi, n_lo


def _rank(a):
    return pd.Series(np.asarray(a, dtype=float)).rank().values


def partial_spearman(x, y, z, rng, n_perm=N_PERM):
    """切片级偏相关（控 z）+ 置换经验双侧 p。

    x = 该臂逐切片 mean score；y = 分期序号；z = 该臂逐切片 mean <img,d>。
    秩残差法：rank 后各自对 rank(z) 作一元回归、取残差，再求 Pearson。
    置换对象是 y（分期标签）—— z 的结构固定，只打乱分期归属。
    """
    z = np.asarray(z, dtype=float)
    if z.std() == 0:
        return float('nan'), float('nan')
    rx, ry, rz = _rank(x), _rank(y), _rank(z)
    rzc = rz - rz.mean()
    vrz = (rzc ** 2).sum()

    def resid(u):
        b = (rzc * (u - u.mean())).sum() / vrz
        return u - u.mean() - b * rzc

    ex = resid(rx)
    r_obs = float(np.corrcoef(ex, resid(ry))[0, 1])

    # rank(y_perm) == rank(y)[perm]（y 的取值多重集不变 ⇒ 平均秩只由取值决定）
    perms = np.array([rng.permutation(len(y)) for _ in range(n_perm)])
    rry = ry[perms]                                    # (n_perm, n)
    rryc = rry - rry.mean(1, keepdims=True)
    b = (rryc * rzc).sum(1) / vrz                      # (n_perm,)
    ey = rryc - b[:, None] * rzc
    num = (ey * ex).sum(1)
    den = np.sqrt((ey ** 2).sum(1) * (ex ** 2).sum())
    r_perm = np.where(den > 0, num / np.maximum(den, 1e-300), 0.0)
    p = float((1 + (np.abs(r_perm) >= abs(r_obs) - 1e-12).sum()) / (n_perm + 1))
    return r_obs, p


def main():
    os.makedirs(OUT, exist_ok=True)
    print('读六臂逐 spot 表与图像向量…', flush=True)
    arms = {}
    per_slide = {}
    for a in ARMS:
        df, E = load_arm(a)
        arms[a] = (df, E)
        g = df.groupby('gsm_dir')
        ps = pd.DataFrame({
            'stage': g.slide_stage.first(),
            'n_spot': g.size(),
            'neo_A_mean': g.neo_A.mean(),
            'lesion_B_mean': g.lesion_B.mean(),
            'chroma_med': g.chroma_std.median(),
        })
        ps['stage_num'] = ps.stage.map(STAGE_NUM)
        if ps.stage_num.isna().any():
            die(f'{a} 有切片分期不在 {sorted(STAGE_NUM)}：{ps[ps.stage_num.isna()].index.tolist()}')
        per_slide[a] = ps
        print(f'  {a}: {len(df):,} spot / {len(ps)} 张', flush=True)

    # —— 六臂同一 spot 集：配对的前提，先验，别假设 ——
    ref_key = arms[ANCHOR][0][['gsm_dir', 'barcode']].values
    for a in ARMS:
        k = arms[a][0][['gsm_dir', 'barcode']].values
        if k.shape != ref_key.shape or not (k == ref_key).all():
            die(f'{a} 与 {ANCHOR} 的 spot 集/顺序不一致 ⇒ 配对前提不成立，不判读。')
    print('  ✅ 六臂逐 spot 一一对齐（同切片、同 barcode、同行序）', flush=True)

    res = {'script': '07_he_pathology/20_ladder_stats.py',
           'prereg': 'results/07_he_pathology/resolution_sweep/PREREG_v2.md',
           'arms': ARMS, 'anchor': ANCHOR,
           'registration': {
               'd_definition': 'mean(E[tissue_frac>=q75]) - mean(E[tissue_frac<=q25])，单位化；'
                               '逐字沿用 09_density_free_ceiling.py',
               'tissue_frac_used': 'L448 那一列（§5：过滤是 spot 的属性）⇒ 六臂同一分组',
               'partial_corr_level': '切片级（与指标 1 同层）；秩残差法',
               'permutation': f'打乱切片的分期标签，{N_PERM} 次，种子 {RNG_SEED}，双侧经验 p',
               'flat_band': FLAT_BAND,
               'flat_rule': f'|Δrho| <= {FLAT_BAND} 或配对 Wilcoxon p >= 0.05 ⇒ 「实质相同」',
           }}

    # ——— 指标 1：切片级 rho（两套提示词）———
    by_arm = {}
    for a in ARMS:
        ps = per_slide[a]
        rA, pA = spearmanr(ps.stage_num, ps.neo_A_mean)
        rB, pB = spearmanr(ps.stage_num, ps.lesion_B_mean)
        by_arm[a] = {
            'crop_px': CROP_OF[a], 'prep': PREP_OF[a],
            'n_slide': int(len(ps)), 'n_spot': int(ps.n_spot.sum()),
            'rho_A_neo': float(rA), 'p_A_neo': float(pA),
            'rho_B_lesion': float(rB), 'p_B_lesion': float(pB),
            'neo_A_mean': float(ps.neo_A_mean.mean()),
            'lesion_B_mean': float(ps.lesion_B_mean.mean()),
            'chroma_min_median': float(ps.chroma_med.min()),
        }
        print(f'  {a:7s} rho(A套 neo)={rA:+.3f}  rho(B套 lesion)={rB:+.3f}  '
              f'最低切片 chroma 中位={ps.chroma_med.min():.1f}', flush=True)
    res['by_arm'] = by_arm

    # ——— 指标 3：管线锚 G1（先于一切判读）———
    r_anchor = by_arm[ANCHOR]['rho_A_neo']
    g1_ok = bool(G1_LO <= r_anchor <= G1_HI)
    res['G1_anchor'] = {'arm': ANCHOR, 'rho': r_anchor, 'band': [G1_LO, G1_HI], 'pass': g1_ok}
    print(f'\nG1 管线锚：{ANCHOR} 的 rho(A 套) = {r_anchor:.3f}，要求 ∈ '
          f'[{G1_LO}, {G1_HI}] ⇒ {"过 ✅" if g1_ok else "越界 ❌"}', flush=True)
    if not g1_ok:
        res['status'] = 'G1_FAILED'
        res['note'] = ('锚越界 ⇒ 按 §10 **硬停，查代码**，不出判读表、不做任何分辨率结论。'
                       '原始测量值已逐张落盘。')
        json.dump(res, open(f'{OUT}/ladder_stats.json', 'w'), indent=2, ensure_ascii=False)
        die(f'G1 锚越界（rho={r_anchor:.3f} 不在 [{G1_LO}, {G1_HI}]）。'
            f'数值已写入 {OUT}/ladder_stats.json。\n'
            f'        可能原因：切片集与旧基线不同、裁框/缩放/提示词漂移、抽样 bug。'
            f'先查代码，**不得**据此下分辨率结论。')

    # ——— 指标 2：扣密度偏相关（逐臂，两套提示词）———
    rng = np.random.default_rng(RNG_SEED)
    print(f'\n扣密度偏相关（切片级，控该臂 <img,d>；置换 {N_PERM} 次）…', flush=True)
    for a in ARMS:
        df, E = arms[a]
        tf = df.tissue_frac.values.astype(float)
        d, n_hi, n_lo = density_direction(E, tf)
        if d is None:
            by_arm[a]['density'] = None
            print(f'  {a}: tissue_frac 分不出高低组 ⇒ 偏相关跳过')
            continue
        proj = E @ d
        ps = per_slide[a]
        pj = pd.Series(proj).groupby(df.gsm_dir.values).mean().reindex(ps.index)
        selfcheck = float(np.corrcoef(proj, tf)[0, 1])
        prA, ppA = partial_spearman(ps.neo_A_mean.values, ps.stage_num.values, pj.values, rng)
        prB, ppB = partial_spearman(ps.lesion_B_mean.values, ps.stage_num.values, pj.values, rng)
        by_arm[a]['density'] = {
            'n_hi': n_hi, 'n_lo': n_lo,
            'corr_proj_tissue_frac': selfcheck,
            'partial_rho_A_neo': prA, 'partial_p_A_neo': ppA,
            'partial_rho_B_lesion': prB, 'partial_p_B_lesion': ppB}
        print(f'  {a:7s} |<img,d> vs tissue_frac| 自检={selfcheck:+.3f}  '
              f'偏相关 A={prA:+.3f} (p={ppA:.2e})  B={prB:+.3f} (p={ppB:.2e})', flush=True)

    # ——— 指标 4：argmax 类目改变率（对各臂 vs 锚）———
    print('\nargmax 类目改变率（同一 spot，对各臂 vs 锚）：', flush=True)
    for a in ARMS:
        if a == ANCHOR:
            continue
        A0 = arms[ANCHOR][0]
        Aa = arms[a][0]
        for tag, idxs, pref in (('A', NEO_IDX + [1, 4, 5, 6], 'a'), ('B', LESION_IDX + [1, 4, 5, 6, 7], 'b')):
            c0 = A0[[f'{pref}{j}' for j in idxs]].values.argmax(1)
            ca = Aa[[f'{pref}{j}' for j in idxs]].values.argmax(1)
            ch = (c0 != ca)
            by_arm[a][f'argmax_change_{tag}'] = {
                'overall': float(ch.mean()),
                'per_slide_mean': float(pd.Series(ch).groupby(A0.gsm_dir.values).mean().mean())}
        print(f'  {a:7s} 改变率 A套={by_arm[a]["argmax_change_A"]["overall"]:.3f}  '
              f'B套={by_arm[a]["argmax_change_B"]["overall"]:.3f}', flush=True)

    # ——— G2 复核（彩色守卫）———
    g2 = {a: bool(by_arm[a]['chroma_min_median'] > G2_MIN_CHROMA) for a in ARMS}
    res['G2_chroma'] = {'min_chroma_required': G2_MIN_CHROMA, 'per_arm_pass': g2,
                        'note': '19 号脚本编码前已逐臂断言；此处用逐 spot 表复核'}
    if not all(g2.values()):
        res['G2_chroma']['failed_arms'] = [a for a, ok in g2.items() if not ok]

    # ——— §8.5 判读表 ———
    def pair(a_from, a_to):
        """a_to 相对 a_from 的 Δrho 与配对检验（逐切片 mean score）。"""
        x = per_slide[a_from].neo_A_mean.values
        y = per_slide[a_to].neo_A_mean.values
        d = y - x
        try:
            _, p = wilcoxon(x, y)
        except ValueError:
            p = float('nan')
        dr = by_arm[a_to]['rho_A_neo'] - by_arm[a_from]['rho_A_neo']
        if abs(dr) <= FLAT_BAND or (p == p and p >= 0.05):
            verdict = 'flat'
        elif dr > 0:
            verdict = 'finer_better'
        else:
            verdict = 'finer_worse'
        return {'from': a_from, 'to': a_to,
                'rho_from': by_arm[a_from]['rho_A_neo'], 'rho_to': by_arm[a_to]['rho_A_neo'],
                'delta_rho': float(dr), 'paired_p': float(p) if p == p else None,
                'median_delta_score': float(np.median(d)), 'verdict': verdict}

    table = {
        'detail_vs_anchor': [pair(ANCHOR, a) for a in ARMS if a != ANCHOR],
        'p1_vs_p0_same_crop': [pair(f'L{c}P0', f'L{c}P1') for c in (448, 320, 224)],
    }
    res['reading_table'] = table
    flat_any = any(r['verdict'] == 'flat' for r in table['detail_vs_anchor'])
    res['wording'] = {
        'if_flat': '在 5.670 µm/px 以上，PLIP 的判别力不随细节变化',
        'forbidden': '细节对病理判别无用（更细的档 ≤1.9 µm/px 本 deposit 测不了，见 PREREG_v2 §9）',
    }

    print('\n§8.5 判读表（A 套 neo，切片级 rho；带宽 %.2f + 配对 Wilcoxon）' % FLAT_BAND)
    print(f'{"对比":>16s} {"rho 起点":>9s} {"rho 终点":>9s} {"Δrho":>8s} {"配对 p":>10s}  判读')
    for r in table['detail_vs_anchor'] + table['p1_vs_p0_same_crop']:
        vd = {'flat': '实质相同', 'finer_better': '更细更好', 'finer_worse': '更细更差'}[r['verdict']]
        print(f'{r["from"]:>7s}->{r["to"]:<7s} {r["rho_from"]:9.3f} {r["rho_to"]:9.3f} '
              f'{r["delta_rho"]:+8.3f} {r["paired_p"]:10.2e}  {vd}')
    if flat_any:
        print(f'\n🔴 措辞纪律（§8.5）：平坦结果必须写成「{res["wording"]["if_flat"]}」，\n'
              f'   不得写成「{res["wording"]["forbidden"]}」。')

    # ——— 落盘 ———
    ps_all = pd.concat([per_slide[a].assign(arm=a) for a in ARMS])
    ps_all.to_csv(f'{OUT}/ladder_by_arm.csv')
    json.dump(res, open(f'{OUT}/ladder_stats.json', 'w'), indent=2, ensure_ascii=False)

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for prep, mk, col in (('P0', 'o', 'tab:blue'), ('P1', '^', 'tab:orange')):
        xs = [CROP_OF[a] for a in ARMS if PREP_OF[a] == prep]
        ys = [by_arm[a]['rho_A_neo'] for a in ARMS if PREP_OF[a] == prep]
        o = np.argsort(xs)
        ax[0].plot(np.array(xs)[o], np.array(ys)[o], mk + '-', color=col,
                   label=f'{prep} (no / background-flattened)')
    ax[0].axhspan(G1_LO, G1_HI, color='grey', alpha=.15)
    ax[0].text(450, G1_HI, 'G1 anchor band', fontsize=7, va='top', ha='right')
    ax[0].set_xlabel('Crop size (hires px)  ->  finer detail with smaller crop')
    ax[0].set_ylabel('Spearman rho (slide mean score vs stage)')
    ax[0].set_title('Detail ladder, prompt set A')
    ax[0].set_xticks([224, 320, 448])
    ax[0].legend(fontsize=8)

    others = [a for a in ARMS if a != ANCHOR]
    w = 0.38
    xx = np.arange(len(others))
    ca = [by_arm[a]['argmax_change_A']['overall'] for a in others]
    cb = [by_arm[a]['argmax_change_B']['overall'] for a in others]
    ax[1].bar(xx - w / 2, ca, w, label='prompt set A')
    ax[1].bar(xx + w / 2, cb, w, label='prompt set B')
    ax[1].set_xticks(xx)
    ax[1].set_xticklabels(others, fontsize=8)
    ax[1].set_ylabel('argmax category change rate vs L448P0')
    ax[1].set_title('Paired argmax shift (same spots)')
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f'{OUT}/ladder_figure.png', dpi=130)

    # ——— 判读文本 ———
    L = []
    L.append('# 分辨率阶梯判读（PREREG_v2 六臂）\n')
    L.append(f'> 由 `07_he_pathology/20_ladder_stats.py` 生成；口径见 `PREREG_v2.md`。')
    L.append(f'> 锚 {ANCHOR} 的 rho = **{r_anchor:.3f}**（要求 ∈ [{G1_LO}, {G1_HI}]）⇒ '
             f'**过**。\n')
    L.append('| 臂 | 裁框 | 预处理 | rho（A 套 neo） | rho（B 套 lesion） | 偏相关 A | 偏相关 p | argmax 改变率 A |')
    L.append('| :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: |')
    for a in ARMS:
        b = by_arm[a]
        dn = b.get('density') or {}
        ch = b.get('argmax_change_A', {}).get('overall')
        L.append(f'| {a} | {b["crop_px"]} | {b["prep"]} | {b["rho_A_neo"]:.3f} | '
                 f'{b["rho_B_lesion"]:.3f} | '
                 f'{dn.get("partial_rho_A_neo", float("nan")):.3f} | '
                 f'{dn.get("partial_p_A_neo", float("nan")):.2e} | '
                 f'{"" if ch is None else f"{ch:.3f}"} |')
    L.append('\n## 判读表（§8.5）\n')
    L.append('| 对比 | rho 起点 | rho 终点 | Δrho | 配对 p | 判读 |')
    L.append('| :--- | ---: | ---: | ---: | ---: | :--- |')
    for r in table['detail_vs_anchor'] + table['p1_vs_p0_same_crop']:
        vd = {'flat': '实质相同', 'finer_better': '更细更好', 'finer_worse': '更细更差'}[r['verdict']]
        L.append(f'| {r["from"]} → {r["to"]} | {r["rho_from"]:.3f} | {r["rho_to"]:.3f} | '
                 f'{r["delta_rho"]:+.3f} | {r["paired_p"]:.2e} | {vd} |')
    L.append(f'\n## 措辞纪律\n\n- 平坦时必须写：「{res["wording"]["if_flat"]}」')
    L.append(f'- **不得**写：「{res["wording"]["forbidden"]}」')
    L.append(f'- 本实验**不**产出 WHO 分型、不判恶性、不解除任何既有禁令（§11）。')
    L.append(f'\n图：`ladder_figure.png`；全部数值：`ladder_stats.json`；逐切片：`ladder_by_arm.csv`。')
    open(f'{OUT}/ladder_reading.md', 'w', encoding='utf-8').write('\n'.join(L) + '\n')

    print(f'\n[out] {OUT}/ladder_stats.json, ladder_by_arm.csv, ladder_reading.md, ladder_figure.png')


if __name__ == '__main__':
    main()
