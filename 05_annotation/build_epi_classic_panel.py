"""生成 `epi_classic_panel.py` —— 上皮亚型**经典 marker** 面板（逐基因可溯源）。

为什么有这个脚本
----------------
用户 2026-09-17 改口径：「就不按什么kac去做，就按照经典marker…按照权威经典marker做重聚类」。
⇒ 面板主干由"源论文 Table S2 meta-program"改为 **Travaglini 2020 Nature 的 Table S1
「Canonical markers」** —— 人工整编的教科书式规范表，逐细胞类型直给基因，非聚类富集输出，
不携带本数据信息（无循环性）。

基因名一律从 .xlsx 单元格**程序化取**，不手抄。

对原表的偏离（2026-09-21 用户裁定后）
------------------------------------
裁定把两类"补"明确区分开，处置不同：

**① 认人同源 —— 保留（1 类偏离）**
`ORTHOLOG_FIX`：Table S1 Club 行的 `CYP2F2` 是**小鼠**基因号（该表其余格均为人的，此格系
物种笔误）；换用人同源 `CYP2F1`。**这不叫"补"** —— `CYP2F1` 就是原表点名那个基因的人版本，
且**同篇论文自己的 Table S4**（用他们自己的人类数据算出的簇富集表）相关两行写的正是 `CYP2F1`。
⇒ 依据充足，保留。（查证过程见 `docs/PARAMETERS_AND_SOURCES.md` §上皮亚型面板 的 2026-09-17 行）

**② 找替身填缺槽 —— 删除（0 类）**
原 `TABLE_S4_FILL`（`DAPL1` 缺则取同篇 Table S4 的 `KRT15`/`KRT17` 顶替）**已按裁定删除**。
⚠️ 它此前**从未被任何代码读取**（生成器写了、产线没用，`grep` 核实）⇒ 删除它
**不改变任何已算结果**，只是清死代码并冻结该裁定。

后果（**两个型被削弱，须随结果报告**）
--------------------------------------
| 型 | Table S1 槽数 | 实际可用 | 丢的是 | 为什么 |
|---|---|---|---|---|
| Basal | 4 | **3** | `DAPL1` | 本矩阵没有；按裁定不取 `KRT15`/`KRT17` 顶替 |
| Serous | 3 | **2** | `PRR4` | 本矩阵没有；按裁定不留替身 |

⇒ **Basal 判别力下降**（用户做此裁定时已被告知）。产线（`03_subcluster_annotation.py`）
对 ≤2 基因的型自动打 `thin_panel`；本文件另设 `KNOWN_ABSENT_IN_MATRIX` 供报告逐条引用。
Club 不受影响（可用 3 个基因）。

用法:
    python3 05_annotation/build_epi_classic_panel.py \
        --tableS1 data/external/travaglini2020/TableS1_canonical_cell_types_markers.xlsx \
        --out 05_annotation/epi_classic_panel.py
"""

import argparse
import hashlib

import pandas as pd

# -------- 小鼠基因号 → 人同源号。2026-09-21 用户裁定：**此类保留**（「只认人同源」）
# 裁定把两种"偏离原表"分开：认人同源号（= 认字）**允许**；找替身填缺槽（= 找替身）**不允许**。
# 依据：Table S1 其余格均为人类基因号，唯 Club 行写 CYP2F2（小鼠），系物种笔误；
#      同篇论文 Table S4（他们自己的人类数据）该行写的正是 CYP2F1。详证见 PARAMETERS 2026-09-17 行。
ORTHOLOG_FIX = {"CYP2F2": "CYP2F1"}

# Table S1 上皮块中**本矩阵肯定没有、且无权威替身可换**的槽（供报告引用；实际缺失由矩阵过滤现算）
# 依据：DAPL1 / PRR4 经 2026-09-17 实测不在作者矩阵（18,069 基因）。
KNOWN_ABSENT_IN_MATRIX = {
    "DAPL1": "本矩阵无此基因；按 2026-09-21 裁定不从 Table S4 取 KRT15/KRT17 顶替（找替身）",
    "PRR4": "本矩阵无此基因；按 2026-09-21 裁定不留替身",
}

# Table S1 上皮块的行号（0-based，实测）
EPI_ROWS = {
    3: "Club", 4: "Ciliated", 5: "Basal", 6: "Goblet", 7: "Mucous", 8: "Serous",
    9: "Ionocyte", 10: "Neuroendocrine", 11: "Tuft",
    12: "AT1", 13: "AT2",
}

# Goblet 与 Mucous 在 Table S1 里共用 MUC5B（Goblet 多 MUC5AC/SPDEF）⇒ 用这组 marker
# **无法**把二者分开。合并为一型，如实登记，不硬拆。
MERGE = {("Goblet", "Mucous"): "Goblet/Mucous"}

# 稀有型（正常肺合计 <0.3%）—— 照算照报，但胜出簇过小时须标 rare/likely-spurious
RARE = ["Serous", "Ionocyte", "Neuroendocrine", "Tuft"]

# 非上皮污染校验（不进 argmax，仅作 QC 列）
CONTAM_CHECK = ["PTPRC", "CD3D", "COL1A1", "PECAM1", "CD68"]


def read_table_s1(path):
    """返回 {亚型: [基因...]}，Table S1 原文逐字，**只**套 `ORTHOLOG_FIX`（物种同源）。"""
    d = pd.read_excel(path, header=None)
    out = {}
    for r, name in EPI_ROWS.items():
        raw = str(d.iloc[r, 3])
        genes = [ORTHOLOG_FIX.get(g.strip(), g.strip()) for g in raw.split(",") if g.strip()]
        if not genes:
            raise ValueError(f"Table S1 第 {r} 行【{name}】未解析出基因，格式可能已变")
        out[name] = genes
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tableS1",
                    default="data/external/travaglini2020/TableS1_canonical_cell_types_markers.xlsx")
    ap.add_argument("--out", default="05_annotation/epi_classic_panel.py")
    a = ap.parse_args()

    s1 = read_table_s1(a.tableS1)

    # 合并 Goblet/Mucous
    merged = {}
    for k, v in s1.items():
        merged[k] = v
    for (x, y), newname in MERGE.items():
        merged[newname] = sorted(set(s1[x]) | set(s1[y]))
        del merged[x], merged[y]

    # 自检：KNOWN_ABSENT_IN_MATRIX 里登记的槽必须真的在原表里（防登记与现实漂移）
    all_genes = {g for gs in merged.values() for g in gs}
    stale = [g for g in KNOWN_ABSENT_IN_MATRIX if g not in all_genes]
    if stale:
        raise SystemExit(f"🔴 KNOWN_ABSENT_IN_MATRIX 里登记为缺失、但原表里根本没有的基因：{stale}")

    order = list(merged.keys())

    L = []
    A = L.append
    A('"""GP8a 上皮亚聚类 —— **经典 marker** 亚型面板（逐基因可溯源到一次文献）。')
    A('')
    A('生成方式：本文件由 `05_annotation/build_epi_classic_panel.py` 从原始补充材料 .xlsx')
    A('**程序化生成**，非手工转录 —— 保证基因名逐字一致。')
    A('')
    A('口径（用户 2026-09-17 指令）')
    A('---------------------------')
    A('「就不按什么kac去做，就按照经典marker…按照权威经典marker做重聚类和celltype那个算法比比结果」')
    A('⇒ 主干 = **Travaglini 2020 Nature 587:619 的 Table S1「Canonical markers」**；')
    A('  源论文的 meta-program / KAC 体系**不再作本面板**（存档于 docs/PARAMETERS_AND_SOURCES.md §M3-A.5b）。')
    A('')
    A('为什么 Table S1 是最合适的"权威经典"')
    A('--------------------------------------')
    A('它是**人工整编的教科书式规范表**（按 Epithelium / Endothelium / Stroma / PNS / Immune 分块），')
    A('**不是**聚类富集输出 ⇒ 不携带本数据的任何信息，作面板无循环性。')
    A('')
    A('⚠️ 三条必须随结果报告的限制')
    A('----------------------------')
    A('1. **该表面向健康人肺，不是肿瘤**。LUAD 中相当一部分上皮是恶性的，用正常 marker 给恶性细胞')
    A('   打标签本身有偏。')
    A('2. **本面板无法区分恶性/正常** —— LUAD 常保留 NKX2-1/SFTPC 等肺泡身份。⇒ 本面板只出"正常上皮')
    A('   亚型"标签，**不设 tumor 亚型**；恶性身份一律由 GP2 CNV 判（铁律 R2），两者**交叉报告**。')
    A('3. **稀有型**（Serous/Ionocyte/Neuroendocrine/Tuft）正常肺合计 <0.3%，照算照报，但胜出簇过小时')
    A('   须标 `rare/likely-spurious`。')
    A('')
    A('对原表的偏离（2026-09-21 用户裁定后）')
    A('------------------------------------')
    A('裁定把两种"偏离原表"分开，处置不同：')
    A('')
    A('**① 认人同源 —— 保留（1 类偏离）**')
    A('`ORTHOLOG_FIX`：Table S1 Club 行的 `CYP2F2` 是**小鼠**基因号（该表其余格均为人的，此格系')
    A('物种笔误）；换用人同源 `CYP2F1`。**这不叫"补"** —— `CYP2F1` 就是原表点名那个基因的人版本，')
    A('且**同篇论文自己的 Table S4**（用他们自己的人类数据算出的簇富集表）相关两行写的正是 `CYP2F1`。')
    A('⇒ 依据充足，保留。')
    A('')
    A('**② 找替身填缺槽 —— 删除（0 类）**')
    A('原 `TABLE_S4_FILL`（`DAPL1` 缺则取同篇 Table S4 的 `KRT15`/`KRT17` 顶替）**已按裁定删除**。')
    A('⚠️ 它此前**从未被任何代码读取**（生成器写了、产线没用，`grep` 核实）⇒ 删除它')
    A('**不改变任何已算结果**，只是清死代码并冻结该裁定。')
    A('')
    A('🔴 后果：两个型被削弱（**须随结果报告**）')
    A('----------------------------------------')
    A('| 型 | Table S1 槽 | 实际可用 | 丢的是 | 为什么 |')
    A('|---|---|---|---|---|')
    A('| Basal | 4 | **3** | `DAPL1` | 本矩阵没有；按裁定不取 `KRT15`/`KRT17` 顶替 |')
    A('| Serous | 3 | **2** | `PRR4` | 本矩阵没有；按裁定不留替身 |')
    A('')
    A('⇒ **Basal 判别力下降**（用户做此裁定时已被告知）。产线对 ≤2 基因的型自动打')
    A('  `thin_panel`；`KNOWN_ABSENT_IN_MATRIX` 供报告逐条引用。Club 不受影响（可用 3 个基因）。')
    A('"""')
    A('')
    A('# ---------------------------------------------------------------- 一次文献')
    A('REFERENCES = {')
    A('    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "')
    A('                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "')
    A('                      "doi:10.1038/s41586-020-2922-4；PMC7704697 "')
    A('                      "【本面板主干 Table S1 Canonical markers；Club 行另有物种核对见 Table S4】",')
    A('    "Habermann2020": "Habermann AC, Gutierrez AJ, Bui LT, et al. Single-cell RNA sequencing reveals "')
    A('                     "profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary "')
    A('                     "fibrosis. Sci Adv 2020;6(28):eaba1972. doi:10.1126/sciadv.aba1972 "')
    A('                     "【仅用于对 Table S1 存疑条目做第二意见核对，不入本面板】",')
    A('}')
    A('')
    A('# -------- 小鼠基因号 → 人同源号。2026-09-21 用户裁定：**此类保留**（「只认人同源」）')
    A('# 裁定把两种"偏离原表"分开：认人同源号（= 认字）允许；找替身填缺槽（= 找替身）不允许。')
    A('# 依据：Table S1 其余格均为人类基因号，唯 Club 行写 CYP2F2（小鼠），系物种笔误；')
    A('#      同篇论文 Table S4（他们自己的人类数据）该行写的正是 CYP2F1。')
    A('ORTHOLOG_FIX = {')
    for k, v in ORTHOLOG_FIX.items():
        A(f'    "{k}": "{v}",')
    A('}')
    A('')
    A('# Table S1 上皮块中**本矩阵肯定没有**的槽（供报告引用；实际缺失由矩阵过滤现算）')
    A('KNOWN_ABSENT_IN_MATRIX = {')
    for k, why in KNOWN_ABSENT_IN_MATRIX.items():
        A(f'    "{k}": "{why}",')
    A('}')
    A('')
    A('# ------------------------------------------------- Table S1 经典 marker（主干，逐字）')
    A('TABLE_S1_MARKERS = {')
    for k in order:
        A(f'    "{k}": [' + ", ".join(f'"{g}"' for g in merged[k]) + '],')
    A('}')
    A('')
    A('# 面板 = Table S1 主干（仅 1 类偏离：Club 行小鼠号→人同源；无 Table S4 补齐）')
    A('')
    A('# 合并说明（Goblet 与 Mucous 在 Table S1 里共用 MUC5B ⇒ 无法分开）')
    A('MERGED = {')
    for (x, y), newname in MERGE.items():
        A(f'    "{x}" + "{y}": "{newname}",')
    A('}')
    A('')
    A(f'RARE = {RARE!r}')
    A(f'CONTAM_CHECK = {CONTAM_CHECK!r}')
    A('')
    A('# 逐型基因槽数（Table S1 原文，未增补前）')
    A('N_SLOTS = {')
    for k in order:
        A(f'    "{k}": {len(merged[k])},')
    A('}')
    A('')
    A('')
    A('def build_panel(available_genes=None, table_s1=None):')
    A('    """返回 {亚型: [基因...]}。')
    A('')
    A('    available_genes : 本数据矩阵的基因集合。给了则**只保留在矩阵里的**，')
    A('                      并返回 (面板, 缺失清单)；不给则原样返回 Table S1 主干。')
    A('    """')
    A('    panel = {k: list(v) if table_s1 is None else list(table_s1[k])')
    A('             for k, v in TABLE_S1_MARKERS.items()}')
    A('    if available_genes is None:')
    A('        return panel')
    A('    have = set(available_genes)')
    A('    used, missing = {}, {}')
    A('    for k, gs in panel.items():')
    A('        used[k] = [g for g in gs if g in have]')
    A('        miss = [g for g in gs if g not in have]')
    A('        if miss:')
    A('            missing[k] = miss')
    A('    return used, missing')
    A('')
    A('')
    A('if __name__ == "__main__":')
    A('    tot = sum(N_SLOTS.values())')
    A('    print(f"经典 marker 面板：{len(N_SLOTS)} 个上皮亚型，基因槽 {tot} 个")')
    A('    for k, n in N_SLOTS.items():')
    A(f'        tag = "  [稀有]" if k in RARE else ""')
    A('        print(f"  {k:16s} {n} 个{tag}")')
    text = "\n".join(L) + "\n"
    with open(a.out, "w") as fh:
        fh.write(text)
    print(f"写出 {a.out}  ({len(text)} 字节)")
    print("out sha256:", hashlib.sha256(text.encode()).hexdigest())
    print("TableS1 sha256:", hashlib.sha256(open(a.tableS1, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
