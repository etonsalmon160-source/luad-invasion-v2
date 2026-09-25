#!/usr/bin/env python3
"""分辨率阶梯判读 · 补记（B 套判读表 / 逐期均值 / 下端脆弱点）。

为什么有这个脚本
----------------
`20_ladder_stats.py`（首跑的判读脚本）产出的 §8.5 判读表**只覆盖 A 套**提示词
（`neo_A`）。而 B 套（原文病理学家逐字形态判据，`lesion_B`）随细节的方向与
A 套**相反**；只看那张表，会把结论读成「细节对判别无用」——与实测不符。
另外那张表里 448 档出现两次（`detail_vs_anchor` 首行与
`p1_vs_p0_same_crop` 首行重叠）。

本脚本把 B 套判读表、逐期均值、以及下端脆弱点（Normal n=1）补进
`ladder_reading.md`，让冻结记录不被读成比实测更狠的结论。

纪律
----
· **只读**冻结的 `ladder_by_arm.csv` 与 `ladder_stats.json`；
  **不重算任何门值**（G1/G2）、**不改变任何既有判读**。
· 所有 rho 一律从 `ladder_stats.json` 取；同时从 CSV 重算一遍做交叉核对，
  不一致就**硬停**（不信任、要验）。
· 判读规则与 20 号脚本**逐字相同**：`|Δrho| <= 0.05` 或配对 Wilcoxon
  `p >= 0.05` ⇒「实质相同」。配对检验用逐切片 mean score（与 20 号一致）。
· 追加是**幂等**的：重复运行不会叠出第二份补记。

跑法：python3 07_he_pathology/21_ladder_supplement.py
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon

ROOT = '/home/eto/luad_v2'
OUT = f'{ROOT}/results/07_he_pathology/resolution_sweep'

ARMS = ['L448P0', 'L448P1', 'L320P0', 'L320P1', 'L224P0', 'L224P1']
ANCHOR = 'L448P0'
CROP_OF = {a: int(a[1:4]) for a in ARMS}
PREP_OF = {a: a[4:] for a in ARMS}
STAGES = ['Normal', 'AAH', 'AIS', 'MIA', 'LUAD']
N_SLIDE_EXPECT = 56
FLAT_BAND = 0.05
MARK = '\n---\n\n## 补记（'

# 两套提示词的列名与 rho 键名
SETS = {'A': ('neo_A_mean', 'rho_A_neo', 'A 套（原口径分期名术语）'),
        'B': ('lesion_B_mean', 'rho_B_lesion', 'B 套（原文病理学家逐字判据）')}


def die(msg):
    print(f'\n[硬停] {msg}', flush=True)
    sys.exit(2)


def w_p(x, y):
    try:
        r = wilcoxon(x, y)
    except ValueError:
        return float('nan')
    return float(getattr(r, 'pvalue', r[1]))


def main():
    csv_p = f'{OUT}/ladder_by_arm.csv'
    json_p = f'{OUT}/ladder_stats.json'
    md_p = f'{OUT}/ladder_reading.md'
    for p in (csv_p, json_p, md_p):
        if not os.path.exists(p):
            die(f'缺 {p} —— 先跑 20_ladder_stats.py。')

    df = pd.read_csv(csv_p)
    J = json.load(open(json_p, encoding='utf-8'))

    # ——— 前置核对：这批数就是冻结的那批 ———
    if sorted(df.arm.unique()) != sorted(ARMS):
        die(f'CSV 里的臂不是 {ARMS}：{sorted(df.arm.unique())}')
    n_slide = df.groupby('arm').gsm_dir.nunique()
    if not (n_slide == N_SLIDE_EXPECT).all():
        die(f'逐臂切片数不是 {N_SLIDE_EXPECT}：{n_slide.to_dict()}')
    if J.get('anchor') != ANCHOR:
        die(f'JSON 的 anchor 是 {J.get("anchor")}，与脚本假定的 {ANCHOR} 不符。')
    if not J.get('G1_anchor', {}).get('pass'):
        die('冻结的 G1 锚没过 ⇒ 按 §10 不该出判读表，本补记也不出。')

    piv = {k: df.pivot(index='gsm_dir', columns='arm', values=v[0])[ARMS]
           for k, v in SETS.items()}
    for k, P in piv.items():
        if P.isna().any().any():
            die(f'{k} 套：有切片在某只臂上缺值 ⇒ 六臂不是同一切片集，不配对。')
    stage_of = df.groupby('gsm_dir').stage.first()
    snum_of = df.groupby('gsm_dir').stage_num.first()
    if not (piv['A'].index == piv['B'].index).all():
        die('两套的切片顺序不一致。')

    # ——— 交叉核对：CSV 重算的 rho == JSON 里的 rho ———
    chk = []
    for a in ARMS:
        for k, (col, key, _) in SETS.items():
            r_csv = float(spearmanr(snum_of.values, piv[k][a].values)[0])
            r_json = float(J['by_arm'][a][key])
            if abs(r_csv - r_json) > 1e-9:
                die(f'{a} {k} 套：CSV 重算 rho={r_csv!r} ≠ JSON {r_json!r} ⇒ '
                    f'两份冻结产物对不上，不追加。')
            chk.append({'arm': a, 'set': k, 'rho': r_json})
    print(f'✅ 交叉核对通过：{len(chk)} 个 rho 与 ladder_stats.json 逐位一致', flush=True)

    # ——— 判读（规则逐字同 20 号脚本的 pair()）———
    def pair(k, a_from, a_to):
        _, key, _ = SETS[k]
        x = piv[k][a_from].values.astype(float)
        y = piv[k][a_to].values.astype(float)
        dr = float(J['by_arm'][a_to][key]) - float(J['by_arm'][a_from][key])
        p = w_p(x, y)
        if abs(dr) <= FLAT_BAND or (p == p and p >= 0.05):
            v = 'flat'
        elif dr > 0:
            v = 'finer_better'
        else:
            v = 'finer_worse'
        return {'set': k, 'from': a_from, 'to': a_to,
                'rho_from': float(J['by_arm'][a_from][key]),
                'rho_to': float(J['by_arm'][a_to][key]),
                'delta_rho': dr, 'paired_p': None if p != p else p,
                'median_delta_score': float(np.median(y - x)), 'verdict': v}

    VS_ANCHOR = [(ANCHOR, a) for a in ARMS if a != ANCHOR]
    STEPWISE = [(f'L{c0}P{p}', f'L{c1}P{p}') for p in ('0', '1')
                for c0, c1 in ((448, 320), (320, 224))]
    P1_VS_P0 = [(f'L{c}P0', f'L{c}P1') for c in (448, 320, 224)]

    tab = {k: {'vs_anchor': [pair(k, *pr) for pr in VS_ANCHOR],
               'stepwise': [pair(k, *pr) for pr in STEPWISE],
               'p1_vs_p0': [pair(k, *pr) for pr in P1_VS_P0]} for k in SETS}

    # ——— 逐期均值 ———
    per_stage = {}
    for k, (col, _, _) in SETS.items():
        g = df.groupby(['arm', 'stage'])[col].mean().unstack('stage')
        n = df.groupby(['arm', 'stage']).size().unstack('stage')
        per_stage[k] = {'mean': g.reindex(columns=STAGES).to_dict('index'),
                        'n': n.reindex(columns=STAGES).to_dict('index')}
    stage_n = df.groupby('stage').gsm_dir.nunique().reindex(STAGES)

    # ——— 落盘（JSON）———
    sup = {'script': '07_he_pathology/21_ladder_supplement.py',
           'date': '2026-09-25',
           'reads_only': ['ladder_by_arm.csv', 'ladder_stats.json'],
           'recomputes_gates': False,
           'flat_rule': f'|Δrho| <= {FLAT_BAND} 或配对 Wilcoxon p >= 0.05 ⇒「实质相同」',
           'cross_check': {'n_rho_checked': len(chk), 'all_match': True},
           'reading_table': tab,
           'per_stage_mean': {k: per_stage[k]['mean'] for k in SETS},
           'per_stage_n_slide': {s: int(stage_n[s]) for s in STAGES},
           'known_defect_in_20_ladder_stats': (
               '§8.5 判读表只算 A 套；且 448 档在 detail_vs_anchor 与 '
               'p1_vs_p0_same_crop 两张子表里各出现一次（表中首行与末行重复）')}
    json.dump(sup, open(f'{OUT}/ladder_supplement.json', 'w', encoding='utf-8'),
              indent=2, ensure_ascii=False)

    # ——— 图（图内文字全英文：本机无 CJK 字体）———
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for ax, k in zip(axes, ('A', 'B')):
        for a in ARMS:
            ys = [per_stage[k]['mean'][a][s] for s in STAGES]
            ax.plot(range(len(STAGES)), ys,
                    color=f'C{[448, 320, 224].index(CROP_OF[a])}',
                    ls='-' if PREP_OF[a] == 'P0' else '--',
                    marker='o' if PREP_OF[a] == 'P0' else '^',
                    ms=4, lw=1.4, alpha=.85,
                    label=f'{a}')
        ax.set_xticks(range(len(STAGES)))
        ax.set_xticklabels(STAGES)
        ax.set_xlabel('Stage (Normal -> LUAD)')
        ax.set_ylabel(f'mean score ({SETS[k][0]})')
        ax.set_title(f'Per-stage mean, prompt set {k}')
        ax.legend(fontsize=7, ncol=2)
        ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(f'{OUT}/ladder_supplement_figure.png', dpi=130)
    plt.close(fig)

    # ——— 追加到判读文本（幂等）———
    vd = {'flat': '实质相同', 'finer_better': '更细更好', 'finer_worse': '更细更差'}
    L = []
    L.append('\n---\n')
    L.append('## 补记（2026-09-25）：B 套判读表 · 逐期均值 · 下端脆弱点\n')
    L.append('> 由 `07_he_pathology/21_ladder_supplement.py` 生成。**只读**本目录冻结的')
    L.append('> `ladder_by_arm.csv` 与 `ladder_stats.json`，**不重算任何门值、不改变任何既有判读**。')
    L.append('>')
    L.append('> **为什么补**：上面的 §8.5 判读表**只有 A 套**（原口径 `neo_A`），')
    L.append('> 且 448 档出现两次（`detail_vs_anchor` 首行与 `p1_vs_p0_same_crop` 首行重叠）。')
    L.append('> **B 套**（原文病理学家逐字形态判据）随细节的方向与 A 套**相反**；')
    L.append('> 只看上面的表，会把结论读成「细节无用」—— 与实测不符。\n')

    L.append('### 1. 两套提示词：方向相反，天花板相同\n')
    L.append('| 提示词 | 448 | 320 | 224 | 变化 |')
    L.append('| :--- | ---: | ---: | ---: | :--- |')
    for k in ('A', 'B'):
        for prep in ('P0', 'P1'):
            key = SETS[k][1]
            ys = [J['by_arm'][f'L{c}{prep}'][key] for c in (448, 320, 224)]
            L.append(f'| {SETS[k][2]} · {prep} | {ys[0]:.3f} | {ys[1]:.3f} | {ys[2]:.3f} | '
                     f'{ys[2] - ys[0]:+.3f} |')
    cells = [(float(J['by_arm'][a][SETS[k][1]]), a, k) for a in ARMS for k in SETS]
    vals = [c[0] for c in cells]
    v_max, a_max, k_max = max(cells)
    L.append(f'\n**六臂 × 两套 = {len(vals)} 格全部落在 {min(vals):.3f}–{v_max:.3f}，'
             f'最高 {v_max:.3f}**（{a_max} / {k_max} 套）。')
    L.append(f'⇒ {len(vals)} 格里**没有一格达到 0.74**；可写'
             f'「换哪个旋钮都没把天花板抬起来」，**不可**写「细节对病理判别无用」。\n')

    # A 套的 vs 锚 / P1 vs P0 已在上方 §8.5 表里（含一行 448 重复），此处只补它缺的「相邻档」
    GROUPS = {'A': [('stepwise', '相邻档')],
              'B': [('vs_anchor', 'vs 锚'), ('stepwise', '相邻档'), ('p1_vs_p0', 'P1 vs P0')]}

    def vword(grp, v):
        """P1 vs P0 那组比的是**预处理**，不是细节 —— 措辞必须分开。"""
        if grp == 'p1_vs_p0':
            return {'flat': '实质相同', 'finer_better': 'P1 更好',
                    'finer_worse': 'P1 更差'}[v]
        return vd[v]

    for k in ('A', 'B'):
        L.append(f'### {"2" if k == "A" else "3"}. {SETS[k][2]} · 判读表（规则逐字同 §8.5）\n')
        if k == 'A':
            L.append('> A 套的「vs 锚」与「P1 vs P0」两组已见上方 §8.5 判读表'
                     '（该表含一行 448 重复）；下面只补它缺的**相邻档**视图。\n')
        L.append('| 对比 | rho 起点 | rho 终点 | Δrho | 配对 p | 判读 |')
        L.append('| :--- | ---: | ---: | ---: | ---: | :--- |')
        for gkey, nm in GROUPS[k]:
            for r in tab[k][gkey]:
                pp = '—' if r['paired_p'] is None else f'{r["paired_p"]:.2e}'
                L.append(f'| {r["from"]} → {r["to"]}（{nm}） | {r["rho_from"]:.3f} | '
                         f'{r["rho_to"]:.3f} | {r["delta_rho"]:+.3f} | {pp} | '
                         f'{vword(gkey, r["verdict"])} |')
        if k == 'A':
            L.append('\n**A 套在每一个相邻档上都是「平、且略降」**（−0.018 / −0.016 / '
                     '−0.007 / −0.021）⇒ 原口径下细节不带来判别力。')
        else:
            L.append('\n🔴 **B 套的升幅整个发生在 320→224 那一档**（+0.088，配对 p=5.0e-09）；'
                     '448→320 只有 +0.047（在 ±0.05 带宽内）。')
            L.append('而 **224 正是唯一不重采样的一档**（1:1），448 与 320 都要先缩图。')
            L.append('⇒ B 套的读数更像「**降采样在吃分**」，不是「越细越好」的单调规律。')
        L.append('')

    L.append('### 4. 逐期均值（切片级 mean score，再按分期取平均）\n')
    for k in ('A', 'B'):
        L.append(f'**{SETS[k][2]}** — `{SETS[k][0]}`\n')
        L.append('| 分期 | 切片 n | ' + ' | '.join(ARMS) + ' |')
        L.append('| :--- | ---: | ' + ' | '.join(['---:'] * len(ARMS)) + ' |')
        for s in STAGES:
            row = [f'{per_stage[k]["mean"][a][s]:.3f}' for a in ARMS]
            L.append(f'| {s} | {int(stage_n[s])} | ' + ' | '.join(row) + ' |')
        L.append('')

    L.append('### 5. 🔴 下端脆弱点（须随结论一起引用）\n')
    L.append('- **Normal 只有 %d 张切片**，MIA 只有 %d 张。每一条 rho 的**下端由这 %d 张 '
             'Normal 承重**；切片级 n=%d，但分期的下端 n=%d。'
             % (int(stage_n['Normal']), int(stage_n['MIA']), int(stage_n['Normal']),
                N_SLIDE_EXPECT, int(stage_n['Normal'])))
    L.append('- 故本记录能支持的是「**前驱三期（AAH/AIS/MIA）彼此分不开**」；')
    L.append('  任何依赖 Normal 具体位置的**幅度**陈述（如「良性被压到某值」）**不稳**。')
    L.append('- 切片内 spot 不独立（空间自相关），统计单位是**切片**（PREREG_v2 §9.4）。')
    L.append('- 本补记**不**产出 WHO 分型、不判恶性、不解除任何既有禁令（§11）。\n')
    L.append('数值：`ladder_supplement.json`；图：`ladder_supplement_figure.png`。')

    base = open(md_p, encoding='utf-8').read()
    if '# 分辨率阶梯判读' not in base:
        die(f'{md_p} 不像冻结的判读文件（缺标题）⇒ 不追加。')
    head = base.split(MARK)[0].rstrip('\n')
    with open(md_p, 'w', encoding='utf-8') as fh:
        fh.write(head + '\n' + '\n'.join(L) + '\n')

    print(f'\n[out] {OUT}/ladder_reading.md（已追加补记）')
    print(f'[out] {OUT}/ladder_supplement.json, ladder_supplement_figure.png')
    for k in ('A', 'B'):
        print(f'\n{k} 套 相邻档：')
        for r in tab[k]['stepwise']:
            print(f'  {r["from"]}->{r["to"]}  Δrho={r["delta_rho"]:+.3f}  '
                  f'p={r["paired_p"]:.2e}  {vd[r["verdict"]]}')


if __name__ == '__main__':
    main()
