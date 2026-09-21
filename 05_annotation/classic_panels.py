"""全部 6 个谱系的**亚型经典 marker 面板**（逐基因可溯源到 Travaglini 2020 Table S1）。

生成方式：由 `05_annotation/build_classic_panels.py` 从原始补充材料 .xlsx **程序化生成**，
非手工转录。**不要手改本文件**，改生成器后重跑。

主干 = Travaglini 2020 Nature 587:619 的 Table S1「Canonical markers」列（逐字）。
口径 = 用户 2026-09-17「就不按什么kac去做，就按照经典marker」+「几个谱系都需要重聚类」。

🔴 对原表的偏离共 **5 类**，逐条登记于 ORTHOLOG_FIX / ALIAS_FIX / DROP_NON_SYMBOL /
   SOURCE_TABLE_DEFECTS / UNANNOTATABLE —— **没有一条是静默的**。
   （原文写「4 类」漏数了 ⓪ ORTHOLOG_FIX，与 §M3-A.5e 的「5 类」自相矛盾；2026-09-17 统一为 5。）
   🔴 **2026-09-21 用户裁定**：认人同源号（认字）**允许**；找替身填缺槽（找替身）**不允许**。
   ⇒ ⓪ ORTHOLOG_FIX（CYP2F2→人同源 CYP2F1）**保留**，5 类全部生效；
     而 Table S4 式的"补替身"**未采用** —— DAPL1 / PRR4 缺失如实上报，不留替身。

🔴 必须随结果报告的两件事：
  1. **非上皮面板比上皮粗**（源表性质）：多型只有 1–2 个基因槽，解读须同 §M3-A.5c 谨慎。
  2. 原表有 **2 处硬伤**（CD4+ Mem/Eff 行混入 CD8；Basophil 与 Mast 行完全相同），
     **原样保留未修**，见 SOURCE_TABLE_DEFECTS。
"""

REFERENCES = {
    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "
                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "
                      "doi:10.1038/s41586-020-2922-4；PMC7704697 【全部面板的唯一来源】",
}

# ⓪ 物种修正：原表个别行写的是**小鼠**基因号 → 换人同源
# 🔴 2026-09-21 用户裁定：认人同源号（认字）**允许**；找替身填缺槽（找替身）**不允许**。
#    依据：Table S1 其余格均人类号、唯 Club 行写 CYP2F2，系物种笔误；
#         同篇论文 Table S4（他们自己的人类数据）该行写的正是 CYP2F1。
ORTHOLOG_FIX = {
    "CYP2F2": "CYP2F1",
}

# ① 描述性名称/缩写 → 矩阵中确实存在的基因符号（逐个核实过）
ALIAS_FIX = {
    "CD8": "CD8A",
    "CD16": "FCGR3A",
    "MHCII CLEC9A": "CLEC9A",
}

# ② 非符号且真符号不在矩阵 → 丢弃并登记，不找替身
DROP_NON_SYMBOL = {
    "MHCII": "真符号 HLA-DRA / HLA-DRB1 **均不在本矩阵**（10x Flex 探针板未覆盖 MHC-II）。丢弃并登记，不用别的基因顶替。原表 mDC2 行用它。",
}

# ③ 原表自身硬伤 —— **不修**，原样保留 + 标记
SOURCE_TABLE_DEFECTS = {
    "CD4+ Mem/Eff": "🔴 原表该行 marker 为「CD3E, CD8, COTL1, LDHB」——**一个 CD4 细胞类型列了 CD8**，属源表笔误。**不修**（改了就成了我们的判断）。整行存疑，注释时不得据它下结论。",
    "Basophil/Mast": "🔴 原表 `Basophil`(行47) 与 `Mast Cell`(行48) 两行**基因完全相同**（MS4A2, CPA3, TPSAB1）⇒ **结构上不可分辨**，合并上报。且 TPSAB1 **不在矩阵** ⇒ 实际仅剩 MS4A2+CPA3 两个基因槽。",
}

# ④ 原表 markers 列为空 → 无法注释
UNANNOTATABLE = {
    "Bronchial Vessel": "🔴 Table S1 第 20 行 markers 列为**空**（连基因槽都没有）⇒ 本面板**无法注释**该亚型。如实上报，不编造 marker。",
}

MERGE = {
    "Goblet + Mucous": "Goblet/Mucous",
    "Basophil + Mast": "Basophil/Mast",
}

RARE = ['Eosinophil', 'Ionocyte', 'Megakaryocyte', 'Neuroendocrine', 'Serous', 'Tuft']

CONTAM_CHECK = ['PTPRC', 'CD3D', 'COL1A1', 'PECAM1', 'CD68', 'EPCAM']

# --------------------------------------------------- 面板本体（逐字取自 Table S1）
PANELS = {
    "上皮": {
        "Club": ["CYP2F1", "SCGB3A2", "CCKAR"],
        "Ciliated": ["FOXJ1", "TUBB1", "TP73", "CCDC78"],
        "Basal": ["KRT5", "KRT14", "TP63", "DAPL1"],
        "Serous": ["PRR4", "LPO", "LTF"],
        "Ionocyte": ["CFTR", "FOXI1", "ASCL3"],
        "Neuroendocrine": ["CALCA", "CHGA", "ASCL1"],
        "Tuft": ["DCLK1", "ASCL2"],
        "AT1": ["AGER", "PDPN", "CLIC5"],
        "AT2": ["SFTPB", "SFTPC", "SFTPD", "MUC1", "ETV5"],
        "Goblet/Mucous": ["MUC5AC", "MUC5B", "SPDEF"],
    },
    "内皮": {
        "Artery": ["GJA5", "BMX"],
        "Vein": ["ACKR1"],
        "Capillary": ["CA4"],
        "Lymphatic": ["PROX1", "PDPN"],
    },
    "成纤维": {
        "Vascular Smooth Muscle": ["CNN1", "ACTA2", "TAGLN", "RGS5"],
        "Airway Smooth Muscle": ["CNN1", "ACTA2", "TAGLN", "DES", "LGR6"],
        "Fibroblast": ["COL1A1", "PDGFRA"],
        "Myofibroblast": ["COL1A1", "PDGFRA", "ELN", "ACTA2"],
        "Lipofibroblast": ["COL1A1", "PDGFRA", "PLIN2", "APOE"],
        "Pericyte": ["CSPG4", "TRPC6", "PDGFRB"],
        "Mesothelial": ["MSLN", "UPK3B", "WT1"],
    },
    "T/NK": {
        "CD8+ Mem/Eff T": ["CD3E", "CD8A", "GZMK", "DUSP2"],
        "CD8+ Naive T": ["CD3E", "CD8A", "GZMH", "GZMB"],
        "CD4+ Mem/Eff": ["CD3E", "CD8A", "COTL1", "LDHB"],
        "CD4+ Naive T": ["CD3E", "CD4", "CCR7", "LEF1"],
        "NK": ["KLRD1", "NKG7", "TYROBP"],
        "NKT": ["CD3E", "CD8A", "FCER1G", "TYROBP"],
    },
    "B/浆": {
        "B": ["CD79A", "CD24", "MS4A1", "CD19"],
        "Plasma": ["CD79A", "CD27", "SLAMF7"],
    },
    "髓系": {
        "Neutrophil": ["S100A8", "S100A9", "IFITM2", "FCGR3B"],
        "Eosinophil": ["SIGLEC8"],
        "Megakaryocyte": ["NRGN", "PPBP", "PF4", "OST4"],
        "Macrophage": ["MARCO", "MSR1", "MRC1"],
        "pDC": ["LILRB4", "IRF8", "LILRA4"],
        "mDC1": ["CLEC9A", "LAMP3"],
        "mDC2": ["CD1C", "PLD4"],
        "Classical Monocyte": ["CD14", "S100A8"],
        "Intermediate Monocyte": ["CD14", "S100A8", "FCGR3A"],
        "Nonclassical Monocyte": ["FCGR3A"],
        "Basophil/Mast": ["CPA3", "MS4A2", "TPSAB1"],
    },
}

# 逐亚型基因槽数（供「面板有多粗」一眼可见）
N_SLOTS = {
    "上皮": {
        "Club": 3,
        "Ciliated": 4,
        "Basal": 4,
        "Serous": 3,
        "Ionocyte": 3,
        "Neuroendocrine": 3,
        "Tuft": 2,
        "AT1": 3,
        "AT2": 5,
        "Goblet/Mucous": 3,
    },
    "内皮": {
        "Artery": 2,
        "Vein": 1,
        "Capillary": 1,
        "Lymphatic": 2,
    },
    "成纤维": {
        "Vascular Smooth Muscle": 4,
        "Airway Smooth Muscle": 5,
        "Fibroblast": 2,
        "Myofibroblast": 4,
        "Lipofibroblast": 4,
        "Pericyte": 3,
        "Mesothelial": 3,
    },
    "T/NK": {
        "CD8+ Mem/Eff T": 4,
        "CD8+ Naive T": 4,
        "CD4+ Mem/Eff": 4,
        "CD4+ Naive T": 4,
        "NK": 3,
        "NKT": 4,
    },
    "B/浆": {
        "B": 4,
        "Plasma": 3,
    },
    "髓系": {
        "Neutrophil": 4,
        "Eosinophil": 1,
        "Megakaryocyte": 4,
        "Macrophage": 3,
        "pDC": 3,
        "mDC1": 2,
        "mDC2": 2,
        "Classical Monocyte": 2,
        "Intermediate Monocyte": 3,
        "Nonclassical Monocyte": 1,
        "Basophil/Mast": 3,
    },
}

DROPPED_LOG = {"mDC2": ["MHCII"]}


def build_panel(lineage, available_genes=None):
    """返回某一谱系的面板 {亚型: [基因...]}。

    available_genes : 本数据矩阵的基因集合。给了则只保留在矩阵里的，
                      并返回 (面板, 缺失清单)；不给则原样返回。
    """
    if lineage not in PANELS:
        raise KeyError(f"未登记的谱系：{lineage}；有 {list(PANELS)}")
    panel = {k: list(v) for k, v in PANELS[lineage].items()}
    if available_genes is None:
        return panel
    have = set(available_genes)
    used, missing = {}, {}
    for k, gs in panel.items():
        used[k] = [g for g in gs if g in have]
        miss = [g for g in gs if g not in have]
        if miss:
            missing[k] = miss
    return used, missing


def assert_provenance():
    """硬校验：原表硬伤必须被登记；被丢弃的词必须有记录；面板里不得残留非符号。

    ⚠️ 2026-09-17 修：原校验只查 ALIAS_FIX 与 DROP_NON_SYMBOL，**漏了 ORTHOLOG_FIX**
    （小鼠基因号 → 人同源）。若日后重生成时 CYP2F2 残留未被换掉，原校验放行。
    ⚠️ 2026-09-21：用户裁定认人同源号**允许**（故本判据保留 ORTHOLOG_FIX 这一项）；
    "找替身填缺槽"**不允许**，未引入任何替身表，故判据无需增项。
    """
    flat = {g for d in PANELS.values() for gs in d.values() for g in gs}
    residue = (flat & set(ALIAS_FIX)) | (flat & set(DROP_NON_SYMBOL)) \
        | (flat & set(ORTHOLOG_FIX))
    if residue:
        raise AssertionError(f"面板里残留未处理的非符号/小鼠基因号：{residue}")
    for name, note in SOURCE_TABLE_DEFECTS.items():
        if not note.startswith("🔴"):
            raise AssertionError(f"{name} 的源表硬伤未标记")
    return True


if __name__ == "__main__":
    assert_provenance()
    tot = sum(len(d) for d in PANELS.values())
    print(f"经典亚型面板：{len(PANELS)} 个谱系 / {tot} 个亚型")
    for lin, d in PANELS.items():
        ns = ", ".join(f"{k}({len(v)})" for k, v in d.items())
        print(f"\n  【{lin}】{len(d)} 型")
        print(f"    {ns}")
    print(f"\n无法注释（原表无 marker）：{list(UNANNOTATABLE)}")
    print(f"原表硬伤（未修，原样保留）：{list(SOURCE_TABLE_DEFECTS)}")
    print(f"被丢弃的非符号：{DROPPED_LOG}")
    thin = [(lin, k, len(v)) for lin, d in PANELS.items() for k, v in d.items() if len(v) <= 2]
    print(f"\n🔴 只有 1–2 个基因槽的亚型（容错为零）：")
    for lin, k, n in thin:
        print(f"    {lin:6s} {k:24s} {n} 个")
