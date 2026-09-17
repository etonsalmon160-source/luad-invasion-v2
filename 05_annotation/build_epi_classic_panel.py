"""生成 `epi_classic_panel.py` —— 上皮亚型**经典 marker** 面板（逐基因可溯源）。

为什么有这个脚本
----------------
用户 2026-09-17 改口径：「就不按什么kac去做，就按照经典marker…按照权威经典marker做重聚类」。
⇒ 面板主干由"源论文 Table S2 meta-program"改为 **Travaglini 2020 Nature 的 Table S1
「Canonical markers」** —— 人工整编的教科书式规范表，逐细胞类型直给基因，非聚类富集输出，
不携带本数据信息（无循环性）。

基因名一律从 .xlsx 单元格**程序化取**，不手抄。

两处、且仅两处允许偏离原表，均显式登记（法则0：不静默）
--------------------------------------------------------
1. `ORTHOLOG_FIX`：Table S1 写的是**小鼠**基因号时换用人同源号（CYP2F2 → CYP2F1）。
2. `TABLE_S4_FILL`：Table S1 的基因**不在本数据矩阵**时，从**同一篇 Travaglini 2020**
   的 Table S4（该亚型对应 cluster 的富集基因）取替代 —— 不另找文献、不问犄角旮旯。

用法:
    python3 05_annotation/build_epi_classic_panel.py \
        --tableS1 data/external/travaglini2020/TableS1_canonical_cell_types_markers.xlsx \
        --tableS4 data/external/travaglini2020/TableS4_cluster_enriched_markers.xlsx \
        --out 05_annotation/epi_classic_panel.py
"""

import argparse
import hashlib

import pandas as pd

# -------- 允许的偏离 #1：小鼠基因号 → 人同源号（Table S1 的 Ciliated/Basal 行是人的，唯 Club 行写了小鼠号）
ORTHOLOG_FIX = {
    "CYP2F2": "CYP2F1",   # Table S1 Club Cell 行；CYP2F2 为小鼠，人同源为 CYP2F1
}

# -------- 允许的偏离 #2：Table S1 基因不在本矩阵 → 同一篇 Table S4 同亚型 cluster 取替代
# 格式: 亚型 -> (Table S4 的 sheet 名, 该 sheet 第 0 行的细胞类型名)
TABLE_S4_SHEET = {
    "Basal": ("Cluster 4", "Basal"),      # Cluster 4/5 = Basal / Proximal Basal
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
    """返回 {亚型: [基因...]}，Table S1 原文逐字（仅套 ORTHOLOG_FIX）。"""
    d = pd.read_excel(path, header=None)
    out = {}
    for r, name in EPI_ROWS.items():
        raw = str(d.iloc[r, 3])
        genes = [ORTHOLOG_FIX.get(g.strip(), g.strip()) for g in raw.split(",") if g.strip()]
        if not genes:
            raise ValueError(f"Table S1 第 {r} 行【{name}】未解析出基因，格式可能已变")
        out[name] = genes
    return out


def read_table_s4_sheet(path, sheet):
    """返回该 sheet 的基因列（第 0 行是细胞类型名，第 1 行表头，第 2 行起是基因）。"""
    d = pd.read_excel(path, sheet_name=sheet, header=None)
    label = str(d.iloc[0, 0]).strip()
    genes = [str(v).strip() for v in d.iloc[2:, 0].dropna()]
    return label, genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tableS1",
                    default="data/external/travaglini2020/TableS1_canonical_cell_types_markers.xlsx")
    ap.add_argument("--tableS4",
                    default="data/external/travaglini2020/TableS4_cluster_enriched_markers.xlsx")
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

    # Table S4 补齐（仅对登记过的亚型）
    fills = {}
    for sub, (sheet, expect) in TABLE_S4_SHEET.items():
        label, genes = read_table_s4_sheet(a.tableS4, sheet)
        if label != expect:
            raise ValueError(f"Table S4 {sheet} 首行是【{label}】，与登记的【{expect}】不符")
        fills[sub] = genes

    # 选补齐基因：Table S1 里该亚型**不在本矩阵**的槽，用 Table S4 的**前若干个**同篇基因顶上
    # 这里不查矩阵（builder 不该依赖数据），只把候选全列出来供产线按可用性取；缺失清单由产线报。
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
    A('两处、且仅两处对原表的偏离（显式登记，不静默）')
    A('----------------------------------------------')
    A('- `ORTHOLOG_FIX`：Table S1 Club 行写的是**小鼠**基因号 `CYP2F2`，换用人同源 `CYP2F1`。')
    A('- `TABLE_S4_FILL`：Table S1 的基因若**不在本数据矩阵**，从**同一篇 Travaglini 2020** 的')
    A('  Table S4（该亚型对应 cluster 的富集基因）取替代 —— 不另找文献。')
    A('"""')
    A('')
    A('# ---------------------------------------------------------------- 一次文献')
    A('REFERENCES = {')
    A('    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "')
    A('                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "')
    A('                      "doi:10.1038/s41586-020-2922-4；PMC7704697 "')
    A('                      "【本面板主干 Table S1 Canonical markers；补齐用 Table S4】",')
    A('    "Habermann2020": "Habermann AC, Gutierrez AJ, Bui LT, et al. Single-cell RNA sequencing reveals "')
    A('                     "profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary "')
    A('                     "fibrosis. Sci Adv 2020;6(28):eaba1972. doi:10.1126/sciadv.aba1972 "')
    A('                     "【仅用于对 Table S1 存疑条目做第二意见核对，不入本面板】",')
    A('}')
    A('')
    A('# -------- 允许的偏离 #1：小鼠基因号 → 人同源号（登记，非静默）')
    A('ORTHOLOG_FIX = {')
    for k, v in ORTHOLOG_FIX.items():
        A(f'    "{k}": "{v}",   # Table S1 Club Cell 行；小鼠 → 人同源')
    A('}')
    A('')
    A('# ------------------------------------------------- Table S1 经典 marker（主干，逐字）')
    A('TABLE_S1_MARKERS = {')
    for k in order:
        A(f'    "{k}": [' + ", ".join(f'"{g}"' for g in merged[k]) + '],')
    A('}')
    A('')
    A('# 面板 = Table S1 主干')
    for sub, genes in fills.items():
        A(f'# Table S4 补齐候选（供 {sub} 在 Table S1 基因不在矩阵时顶替；同篇 Travaglini 2020）')
    A('TABLE_S4_FILL = {')
    for sub, genes in fills.items():
        A(f'    "{sub}": [')
        rows = [genes[i:i + 8] for i in range(0, 40, 8)]
        for r in rows:
            A('        ' + ", ".join(f'"{g}"' for g in r) + ',')
        A('    ],')
    A('}')
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
    for nm, p in (("TableS1", a.tableS1), ("TableS4", a.tableS4)):
        print(f"{nm} sha256:", hashlib.sha256(open(p, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
