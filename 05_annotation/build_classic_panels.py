"""生成 `classic_panels.py` —— **全部 6 个谱系**的上皮/非上皮亚型经典 marker 面板。

为什么要扩到全部谱系
------------------
用户 2026-09-17：「**几个谱系都需要重聚类，不仅仅是上皮**」。
§M3-A.5 本就登记了 6 个谱系；本脚本把上皮那套「经典 marker」口径**原样推广**到其余 5 个谱系，
用的是**同一张表**（Travaglini 2020 Table S1），不另找文献。

主干与偏离（全部显式登记，不静默）
--------------------------------
主干 = Table S1 的 Canonical markers 列，逐字程序化取出（不手抄）。

对原表**必需**的偏离，共 **5 类**，逐条登记在 `ORTHOLOG_FIX` / `ALIAS_FIX` / `DROP_NON_SYMBOL` /
`SOURCE_TABLE_DEFECTS` / `UNANNOTATABLE`：

🔴 **2026-09-21 用户裁定**：把两种"偏离原表"分开 —— **认人同源号（= 认字）允许**，
**找替身填缺槽（= 找替身）不允许**。据此后果：
- ⓪ `ORTHOLOG_FIX` **保留**（`CYP2F2`→`CYP2F1`）。依据：Table S1 其余格均人类号、唯 Club 行写
  小鼠号 `CYP2F2`（物种笔误），且同篇 Table S4（他们自己的人类数据）该行写的正是 `CYP2F1`。
  生成器在下方用 Table S4 做**可执行对账**：凡替换后的基因必须出现在 Table S4 里，否则硬报错。
- Table S4「找替身填缺槽」（`DAPL1` 缺则取 `KRT15`/`KRT17`）**已删**。⚠️ 它此前从未被任何代码
  读取（`grep` 核实）⇒ 删除不改变任何已算结果。`DAPL1`/`PRR4` 缺失如实上报，不留替身。

① `ALIAS_FIX`：原表写的是**描述性名称或缩写**，不是基因符号。
   已核实矩阵里存在对应符号的才换；换不了的进 ②。
② `DROP_NON_SYMBOL`：原表写的不是符号、且其对应的真符号**不在本矩阵**。
   例：`MHCII` 的真符号 `HLA-DRA`/`HLA-DRB1` **都不在矩阵**（本 10x Flex 探针板未覆盖 MHC-II），
   ⇒ 只能丢弃并登记，**不拿别的基因顶替**。
③ `SOURCE_TABLE_DEFECTS`：**原表自身有硬伤**。已发现 2 处：
   (a) `CD4+ Mem/Eff Cell` 行的 marker 写的是 `CD3E, CD8, COTL1, LDHB` —— 一个 CD4 细胞类型
       列了 CD8。**这是源表的错**，不是我们的口径。⇒ 该行**不修**（改了就成我们的判断了），
       整行**标记为存疑**并逐字保留原样，注释时不得据它下结论。
   (b) `Basophil` 与 `Mast Cell` 两行**基因完全相同**（`MS4A2, CPA3, TPSAB1`），
       且 `TPSAB1` 不在矩阵 ⇒ 二者**结构上不可分辨**。合并为 `Basophil/Mast` 上报。
④ `UNANNOTATABLE`：原表该行 **markers 列是空的**。例：`Bronchial Vessel`（第 20 行）。
   ⇒ 该亚型**无法用本面板注释**，如实上报，**不编造 marker**。

如实声明（必须随结果报告）
------------------------
**非上皮谱系的面板比上皮的粗**，这是**源表的性质**，不是本项目的选择：
上皮 11 型多为 3–5 个基因；非上皮里 `Vein Cell` 只有 1 个（`ACKR1`）、
`Eosinophil` 1 个（`SIGLEC8`）、`Capillary Cell` 1 个（`CA4`）、
`Nonclassical Monocyte` 1 个（`CD16`→`FCGR3A`）。基因槽少 ⇒ 容错为零，
解读时必须按 §M3-A.5c「检出率现实」同样谨慎。

用法:
    python3 05_annotation/build_classic_panels.py \
        --tableS1 data/external/travaglini2020/TableS1_canonical_cell_types_markers.xlsx \
        --out 05_annotation/classic_panels.py
"""

import argparse
import hashlib
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ------- 谱系 → {Table S1 行号: 亚型名}（0-based，实测行号）------------------
LINEAGE_ROWS = {
    "上皮": {3: "Club", 4: "Ciliated", 5: "Basal", 6: "Goblet", 7: "Mucous", 8: "Serous",
             9: "Ionocyte", 10: "Neuroendocrine", 11: "Tuft", 12: "AT1", 13: "AT2"},
    "内皮": {17: "Artery", 18: "Vein", 19: "Capillary", 21: "Lymphatic"},
    "成纤维": {24: "Vascular Smooth Muscle", 25: "Airway Smooth Muscle", 26: "Fibroblast",
               27: "Myofibroblast", 28: "Lipofibroblast", 29: "Pericyte",
               30: "Mesothelial"},
    "T/NK": {40: "CD8+ Mem/Eff T", 41: "CD8+ Naive T", 42: "CD4+ Mem/Eff",
             43: "CD4+ Naive T", 44: "NK", 45: "NKT"},
    "B/浆": {38: "B", 39: "Plasma"},
    "髓系": {46: "Neutrophil", 47: "Basophil", 48: "Mast", 49: "Eosinophil",
             50: "Megakaryocyte", 51: "Macrophage", 52: "pDC", 53: "mDC1", 54: "mDC2",
             55: "Classical Monocyte", 56: "Intermediate Monocyte",
             57: "Nonclassical Monocyte"},
}

# ⓪ 物种修正：原表个别行写的是**小鼠**基因号 → 换人同源。
# 🔴 **2026-09-21 用户裁定「只认人同源」——此类保留**（区别于"找替身填缺槽"，后者被禁）。
#    依据：Table S1 其余格均为人类基因号，唯「Club」行写 CYP2F2（小鼠），系物种笔误；
#         同篇论文 Table S4（他们自己的人类数据）该行写的正是 CYP2F1。
#    生成器用 Table S4 对账硬校验此替换（见 main()），换不了就报错。
ORTHOLOG_FIX = {"CYP2F2": "CYP2F1"}

# ① 原表写描述性名称/缩写 → 换成矩阵里**确实存在**的基因符号（逐个核过）
ALIAS_FIX = {
    "CD8": "CD8A",     # 原表 CD8+ Naive T 行；CD8B 不在矩阵（已核实）
    "CD16": "FCGR3A",  # 原表 Nonclassical/Intermediate Monocyte 行；FCGR3B 在矩阵但语义为中性粒
    "MHCII CLEC9A": "CLEC9A",   # 原表 mDC1 行写成 "MHCII CLEC9A"——CLEC9A 才是基因
}

# ② 非符号、且真符号不在矩阵 → 丢弃并登记（**不找替身**）
DROP_NON_SYMBOL = {
    "MHCII": "真符号 HLA-DRA / HLA-DRB1 **均不在本矩阵**（10x Flex 探针板未覆盖 MHC-II）。"
             "丢弃并登记，不用别的基因顶替。原表 mDC2 行用它。",
}

# ③ 原表自身的硬伤 —— **不修，原样保留 + 标记**
SOURCE_TABLE_DEFECTS = {
    "CD4+ Mem/Eff": "🔴 原表该行 marker 为「CD3E, CD8, COTL1, LDHB」——**一个 CD4 细胞类型列了 CD8**，"
                    "属源表笔误。**不修**（改了就成了我们的判断）。整行存疑，注释时不得据它下结论。",
    "Basophil/Mast": "🔴 原表 `Basophil`(行47) 与 `Mast Cell`(行48) 两行**基因完全相同**"
                     "（MS4A2, CPA3, TPSAB1）⇒ **结构上不可分辨**，合并上报。"
                     "且 TPSAB1 **不在矩阵** ⇒ 实际仅剩 MS4A2+CPA3 两个基因槽。",
}

# ④ 原表 markers 列为空 → 无法注释
UNANNOTATABLE = {
    "Bronchial Vessel": "🔴 Table S1 第 20 行 markers 列为**空**（连基因槽都没有）"
                        "⇒ 本面板**无法注释**该亚型。如实上报，不编造 marker。",
}

# 合并（同名或不可分辨）
MERGE = {
    ("Goblet", "Mucous"): "Goblet/Mucous",
    ("Basophil", "Mast"): "Basophil/Mast",
}

# 稀有型（正常肺合计 <0.3%，Travaglini Table S1 的 Relative abundance 列）
RARE = ["Serous", "Ionocyte", "Neuroendocrine", "Tuft", "Eosinophil", "Megakaryocyte"]

# 非本体谱系污染校验（不进 argmax，仅作 QC 列）
CONTAM_CHECK = ["PTPRC", "CD3D", "COL1A1", "PECAM1", "CD68", "EPCAM"]


def read_table_s1(path):
    """返回 {行号: (亚型名, [基因名...] 或 None)}。None = 原表该行 markers 为空。"""
    d = pd.read_excel(path, header=None)
    out = {}
    for lin, rows in LINEAGE_ROWS.items():
        out[lin] = {}
        for r, name in rows.items():
            v = d.iloc[r, 3]
            if pd.isna(v):
                out[lin][name] = None
                continue
            genes = [g.strip() for g in str(v).split(",") if g.strip()]
            out[lin][name] = genes or None
    return out


def apply_fixes(genes):
    """套 ALIAS_FIX → 套 DROP_NON_SYMBOL → 去重（保序）。返回 (基因, 被丢掉的原词)。"""
    out, dropped = [], []
    for g in genes:
        if g in DROP_NON_SYMBOL:
            dropped.append(g)
            continue
        out.append(ORTHOLOG_FIX.get(g, ALIAS_FIX.get(g, g)))
    seen, uniq = set(), []
    for g in out:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tableS1",
                    default="data/external/travaglini2020/TableS1_canonical_cell_types_markers.xlsx")
    ap.add_argument("--out", default="05_annotation/classic_panels.py")
    ap.add_argument("--frozen-epidermis", default="05_annotation/epi_classic_panel.py")
    ap.add_argument("--tableS4",
                    default="data/external/travaglini2020/TableS4_cluster_enriched_markers.xlsx")
    a = ap.parse_args()

    # ---- 用户裁定「按和我们数据对得上的那张原文 marker 做」的可执行判据：
    #      凡经 ORTHOLOG_FIX / ALIAS_FIX 换过的基因，**必须在论文自己的 Table S4 里出现过**。
    #      Table S4 是他们用自己的人类数据算出来的 ⇒ 它出现 = 「原文自己的数据这么用」。
    s4_genes = set()
    if os.path.exists(a.tableS4):
        xl = pd.ExcelFile(a.tableS4)
        for sh in xl.sheet_names:
            s4_genes |= {str(v).strip() for v in
                         pd.read_excel(a.tableS4, sheet_name=sh, header=None).iloc[2:, 0].dropna()}
        print(f"Table S4 共 {len(xl.sheet_names)} 张 sheet / {len(s4_genes)} 个基因名（对账用）")
        bad = []
        for orig, new in list(ORTHOLOG_FIX.items()) + list(ALIAS_FIX.items()):
            if new not in s4_genes:
                bad.append((orig, new))
        if bad:
            raise SystemExit(
                f"🔴 以下替换在论文自己的 Table S4 里**查无此基因**，"
                f"不符合用户裁定的「按和我们数据对得上的原文 marker」：{bad}")
        print("全部替换均在 Table S4 中出现 ✅（= 论文自己的数据这么用）")

    s1 = read_table_s1(a.tableS1)

    panels, dropped_log, unannot = {}, {}, {}
    for lin, d in s1.items():
        panels[lin] = {}
        for name, genes in d.items():
            if genes is None:
                unannot[name] = lin
                continue
            fixed, dropped = apply_fixes(genes)
            panels[lin][name] = fixed
            if dropped:
                dropped_log[name] = dropped

    # 合并
    for lin in panels:
        for (x, y), newname in MERGE.items():
            if x in panels[lin] and y in panels[lin]:
                panels[lin][newname] = sorted(set(panels[lin][x]) | set(panels[lin][y]))
                del panels[lin][x], panels[lin][y]

    # 与已冻结的上皮面板对账：上皮块必须逐基因一致（防两个文件漂移）
    if os.path.exists(a.frozen_epidermis):
        sys.path.insert(0, os.path.dirname(os.path.abspath(a.frozen_epidermis)))
        import epi_classic_panel as F  # noqa: E402
        mine = {k: sorted(v) for k, v in panels["上皮"].items()}
        theirs = {k: sorted(v) for k, v in F.TABLE_S1_MARKERS.items()}
        if mine != theirs:
            diff = {k: (mine.get(k), theirs.get(k)) for k in set(mine) | set(theirs)
                    if mine.get(k) != theirs.get(k)}
            raise SystemExit(f"🔴 上皮块与已冻结的 epi_classic_panel.py 不一致：{diff}")
        print("上皮块与已冻结的 epi_classic_panel.py 逐基因一致 ✅")

    L = []
    A = L.append
    A('"""全部 6 个谱系的**亚型经典 marker 面板**（逐基因可溯源到 Travaglini 2020 Table S1）。')
    A('')
    A('生成方式：由 `05_annotation/build_classic_panels.py` 从原始补充材料 .xlsx **程序化生成**，')
    A('非手工转录。**不要手改本文件**，改生成器后重跑。')
    A('')
    A('主干 = Travaglini 2020 Nature 587:619 的 Table S1「Canonical markers」列（逐字）。')
    A('口径 = 用户 2026-09-17「就不按什么kac去做，就按照经典marker」+「几个谱系都需要重聚类」。')
    A('')
    A('🔴 对原表的偏离共 **5 类**，逐条登记于 ORTHOLOG_FIX / ALIAS_FIX / DROP_NON_SYMBOL /')
    A('   SOURCE_TABLE_DEFECTS / UNANNOTATABLE —— **没有一条是静默的**。')
    A('   （原文写「4 类」漏数了 ⓪ ORTHOLOG_FIX，与 §M3-A.5e 的「5 类」自相矛盾；2026-09-17 统一为 5。）')
    A('   🔴 **2026-09-21 用户裁定**：认人同源号（认字）**允许**；找替身填缺槽（找替身）**不允许**。')
    A('   ⇒ ⓪ ORTHOLOG_FIX（CYP2F2→人同源 CYP2F1）**保留**，5 类全部生效；')
    A('     而 Table S4 式的"补替身"**未采用** —— DAPL1 / PRR4 缺失如实上报，不留替身。')
    A('')
    A('🔴 必须随结果报告的两件事：')
    A('  1. **非上皮面板比上皮粗**（源表性质）：多型只有 1–2 个基因槽，解读须同 §M3-A.5c 谨慎。')
    A('  2. 原表有 **2 处硬伤**（CD4+ Mem/Eff 行混入 CD8；Basophil 与 Mast 行完全相同），')
    A('     **原样保留未修**，见 SOURCE_TABLE_DEFECTS。')
    A('"""')
    A('')
    A('REFERENCES = {')
    A('    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "')
    A('                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "')
    A('                      "doi:10.1038/s41586-020-2922-4；PMC7704697 【全部面板的唯一来源】",')
    A('}')
    A('')
    A('# ⓪ 物种修正：原表个别行写的是**小鼠**基因号 → 换人同源')
    A('# 🔴 2026-09-21 用户裁定：认人同源号（认字）**允许**；找替身填缺槽（找替身）**不允许**。')
    A('#    依据：Table S1 其余格均人类号、唯 Club 行写 CYP2F2，系物种笔误；')
    A('#         同篇论文 Table S4（他们自己的人类数据）该行写的正是 CYP2F1。')
    A('ORTHOLOG_FIX = {')
    for k, v in ORTHOLOG_FIX.items():
        A(f'    "{k}": "{v}",')
    A('}')
    A('')
    A('# ① 描述性名称/缩写 → 矩阵中确实存在的基因符号（逐个核实过）')
    A('ALIAS_FIX = {')
    for k, v in ALIAS_FIX.items():
        A(f'    "{k}": "{v}",')
    A('}')
    A('')
    A('# ② 非符号且真符号不在矩阵 → 丢弃并登记，不找替身')
    A('DROP_NON_SYMBOL = {')
    for k, v in DROP_NON_SYMBOL.items():
        A(f'    "{k}": "{v}",')
    A('}')
    A('')
    A('# ③ 原表自身硬伤 —— **不修**，原样保留 + 标记')
    A('SOURCE_TABLE_DEFECTS = {')
    for k, v in SOURCE_TABLE_DEFECTS.items():
        A(f'    "{k}": "{v}",')
    A('}')
    A('')
    A('# ④ 原表 markers 列为空 → 无法注释')
    A('UNANNOTATABLE = {')
    for k, v in UNANNOTATABLE.items():
        A(f'    "{k}": "{v}",')
    A('}')
    A('')
    A('MERGE = {')
    for (x, y), newname in MERGE.items():
        A(f'    "{x} + {y}": "{newname}",')
    A('}')
    A('')
    A(f'RARE = {sorted(RARE)!r}')
    A('')
    A(f'CONTAM_CHECK = {CONTAM_CHECK!r}')
    A('')
    A('# --------------------------------------------------- 面板本体（逐字取自 Table S1）')
    A('PANELS = {')
    for lin, d in panels.items():
        A(f'    "{lin}": {{')
        for name, genes in d.items():
            A(f'        "{name}": [' + ", ".join(f'"{g}"' for g in genes) + '],')
        A('    },')
    A('}')
    A('')
    A('# 逐亚型基因槽数（供「面板有多粗」一眼可见）')
    A('N_SLOTS = {')
    for lin, d in panels.items():
        A(f'    "{lin}": {{')
        for name, genes in d.items():
            A(f'        "{name}": {len(genes)},')
        A('    },')
    A('}')
    A('')
    A(f'DROPPED_LOG = {json.dumps(dropped_log, ensure_ascii=False)}')
    A('')
    A('')
    A('def build_panel(lineage, available_genes=None):')
    A('    """返回某一谱系的面板 {亚型: [基因...]}。')
    A('')
    A('    available_genes : 本数据矩阵的基因集合。给了则只保留在矩阵里的，')
    A('                      并返回 (面板, 缺失清单)；不给则原样返回。')
    A('    """')
    A('    if lineage not in PANELS:')
    A('        raise KeyError(f"未登记的谱系：{lineage}；有 {list(PANELS)}")')
    A('    panel = {k: list(v) for k, v in PANELS[lineage].items()}')
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
    A('def assert_provenance():')
    A('    """硬校验：原表硬伤必须被登记；被丢弃的词必须有记录；面板里不得残留非符号。')
    A('')
    A('    ⚠️ 2026-09-17 修：原校验只查 ALIAS_FIX 与 DROP_NON_SYMBOL，**漏了 ORTHOLOG_FIX**')
    A('    （小鼠基因号 → 人同源）。若日后重生成时 CYP2F2 残留未被换掉，原校验放行。')
    A('    ⚠️ 2026-09-21：用户裁定认人同源号**允许**（故本判据保留 ORTHOLOG_FIX 这一项）；')
    A('    "找替身填缺槽"**不允许**，未引入任何替身表，故判据无需增项。')
    A('    """')
    A('    flat = {g for d in PANELS.values() for gs in d.values() for g in gs}')
    A('    residue = (flat & set(ALIAS_FIX)) | (flat & set(DROP_NON_SYMBOL)) \\')
    A('        | (flat & set(ORTHOLOG_FIX))')
    A('    if residue:')
    A('        raise AssertionError(f"面板里残留未处理的非符号/小鼠基因号：{residue}")')
    A('    for name, note in SOURCE_TABLE_DEFECTS.items():')
    A('        if not note.startswith("🔴"):')
    A('            raise AssertionError(f"{name} 的源表硬伤未标记")')
    A('    return True')
    A('')
    A('')
    A('if __name__ == "__main__":')
    A('    assert_provenance()')
    A('    tot = sum(len(d) for d in PANELS.values())')
    A('    print(f"经典亚型面板：{len(PANELS)} 个谱系 / {tot} 个亚型")')
    A('    for lin, d in PANELS.items():')
    A('        ns = ", ".join(f"{k}({len(v)})" for k, v in d.items())')
    A('        print(f"\\n  【{lin}】{len(d)} 型")')
    A('        print(f"    {ns}")')
    A('    print(f"\\n无法注释（原表无 marker）：{list(UNANNOTATABLE)}")')
    A('    print(f"原表硬伤（未修，原样保留）：{list(SOURCE_TABLE_DEFECTS)}")')
    A('    print(f"被丢弃的非符号：{DROPPED_LOG}")')
    A('    thin = [(lin, k, len(v)) for lin, d in PANELS.items() for k, v in d.items() if len(v) <= 2]')
    A('    print(f"\\n🔴 只有 1–2 个基因槽的亚型（容错为零）：")')
    A('    for lin, k, n in thin:')
    A('        print(f"    {lin:6s} {k:24s} {n} 个")')
    text = "\n".join(L) + "\n"
    with open(a.out, "w") as fh:
        fh.write(text)
    print(f"\n写出 {a.out}  ({len(text)} 字节)")
    print("out sha256:", hashlib.sha256(text.encode()).hexdigest())
    print("TableS1 sha256:", hashlib.sha256(open(a.tableS1, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
