"""GP6 标准 A 的 marker 面板 —— 逐基因可溯源到**一次文献**。

设计约束（重要，不可绕过）
--------------------------------------------------
标准 A（本面板）与标准 B（CellTypist）必须是**两个独立口径**，否则一致性 κ 是自证。
因此本面板的出处**刻意排除 HLCA 整合图谱本身**（Sikkema 2023 Nat Med），
因为 CellTypist 的 `Human_Lung_Atlas.pkl` 就是该图谱训出来的。

残留的非独立性（如实声明）：HLCA 整合了本面板所引的多套一次研究数据，
故两者**并非统计独立**；本面板只能保证「marker 的定义来源不同」，不能保证「底层数据不同」。

出处审计（2026-09-17，用户指令「你给我查清楚marker到底哪里来的」）
====================================================================
审计脚本见 `docs/PARAMETERS_AND_SOURCES.md` §M3-A.5d。结论：

1. **原引用的 `VieiraBraga2019` 无效** —— 该文**没有 marker 表**（补充材料为细胞计数 /
   OMIM 基因表 / 临床元数据 / Fisher p 值 / T 细胞计数 / CellPhoneDB / 抗体表），
   上皮亚型 marker 只出现在正文图里，**不可逐基因追溯**。⇒ 已从 `src` 移除，保留在
   `UNVERIFIED_SRC` 中如实登记，**不静默删除、也不假装核对过**。

2. **80/80 基因全部可追溯到 Travaglini 2020 Table S4**（带 avg_logFC / pct_in_cluster /
   pct_out_cluster / p_val_adj）。⇒ 本面板基因**不是"犄角旮旯"来的**。

3. **谱系归属 71/80 一致**（判据：该基因被报富集的细胞类型中，至少一个属于它应属的谱系）。

4. **3 个基因偏弱，保留但标记**（用户 2026-09-17 决定：只改引用、不换基因、不重跑 GP6）：
   见 `WEAK_PROVENANCE`。

5. Travaglini Table S1 是**人工整编的参考表**，供人工簇级判读查表用（Methods 原文：
   "Clusters were assigned a canonical identity based on enriched expression of these
   marker genes."），**从不喂给自动标注器** ⇒ 作面板无循环性。
   ⚠️ 但它面向**健康人肺**，本数据含恶性细胞，用正常 marker 给恶性细胞打标签本身有偏。

替换历史：本文件取代 `docs/scientific_rigor_and_audit.md` 法则2 的旧表。
旧表**一个出处都没有**（2026-09-17 审计标缺），按用户 2026-09-17 指令「marker 要来自权威文章」
重建为下方带引用的版本。
"""

# ---------------------------------------------------------------- 一次文献清单
REFERENCES = {
    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "
                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "
                      "doi:10.1038/s41586-020-2922-4；PMC7704697 "
                      "【本面板的**可追溯主来源**：Table S1 人工整编规范 marker；Table S4 簇富集 marker】",
    "Habermann2020": "Habermann AC, Gutierrez AJ, Bui LT, et al. Single-cell RNA sequencing reveals "
                     "profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary "
                     "fibrosis. Sci Adv 2020;6(28):eaba1972. doi:10.1126/sciadv.aba1972 "
                     "【Table S3「Marker genes used for cell-type annotation」覆盖全部 6 谱系、"
                     "带正负 marker；但该表在补充 PDF 里是**图片**，无法程序化提取，故本面板**未从它取基因**，"
                     "只作人工核对的第二意见】",
    "Peng2026": "Peng F, Sinjab A, et al., Kadara H. Multimodal spatial-omics reveal co-evolution of "
                "alveolar progenitors and proinflammatory niches in progression of lung precursor "
                "lesions. Cancer Cell 2026;44(2):321-339.e13. doi:10.1016/j.ccell.2025.10.004 "
                "【本项目源论文】",
}

# 🔴 曾被本文件引用、但**本次审计未能核实其有逐基因 marker 表**的文献。
#    保留登记，不静默删除，也不假装核对过。它们的基因支持力**未知**。
UNVERIFIED_SRC = {
    "VieiraBraga2019": "Vieira Braga FA, Kar G, Berg M, et al. A cellular census of human lungs "
                       "identifies novel cell states in health and in asthma. Nat Med "
                       "2019;25(7):1153-1163. doi:10.1038/s41591-019-0468-5 "
                       "🔴 **已核实：该文没有 marker 表**（补充为细胞计数/OMIM 表/临床元数据/"
                       "Fisher p 值/T 细胞计数/CellPhoneDB/抗体表）。**不可逐基因追溯**。",
    "Reyfman2019": "Reyfman PA, Walter JM, Joshi N, et al. Am J Respir Crit Care Med "
                   "2019;199(12):1517-1536. doi:10.1164/rccm.201712-2410OC  ⚠️ 本次**未核对**其补充材料",
    "Zilionis2019": "Zilionis R, Engblom C, Pfirschke C, et al. Immunity 2019;50(5):1317-1334. "
                    "doi:10.1016/j.immuni.2019.03.009  ⚠️ 本次**未核对**其补充材料",
    "Gillich2020": "Gillich A, Zhang F, Farmer CG, et al. Nature 2020;586(7831):785-789. "
                   "doi:10.1038/s41586-020-2822-7  ⚠️ 本次**未核对**其补充材料",
    "Lambrechts2018": "Lambrechts D, Wauters E, Boeckx B, et al. Nat Med 2018;24(8):1277-1289. "
                      "doi:10.1038/s41591-018-0096-5  ⚠️ 本次**未核对**其补充材料",
    "Guo2018": "Guo X, Zhang Y, Zheng L, et al. Nat Med 2018;24(7):978-985. "
               "doi:10.1038/s41591-018-0045-3  ⚠️ 本次**未核对**其补充材料",
}

# 🔴 出处审计中检出**支持偏弱**的基因。用户 2026-09-17 决定：**保留基因、只改引用**。
#    必须随 GP6 结果一并报告，不得当作与其余基因同等强度的 marker。
WEAK_PROVENANCE = {
    "SPP1": "Travaglini Table S4 全表最高 pct_in_cluster 仅 0.21 ⇒ 该数据中不构成特异 marker。"
            "（SPP1+ 巨噬在 Habermann 2020 / Zilionis 2019 另有文献支持，但本次未逐表核对）",
    "FAP": "Travaglini Table S4 中 logFC 1.67 但 pct_in_cluster 仅 0.15 ⇒ 表达比例过低。",
    "KDR": "Travaglini Table S4 中 logFC 最高的细胞类型是 Goblet（生物学别扭）⇒ 归属不稳。",
}

# ---------------------------------------------------------------- 面板本体
# 每谱系：genes = 阳性 marker（用于 sc.tl.score_genes 模块分）
#         neg   = 应为阴性 marker（仅作校验展示，不进模块分）
#         src   = 出处键（见 REFERENCES）
PANEL = {
    "上皮": dict(
        genes=["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "NKX2-1", "SFTPC", "SFTPA1",
               "SFTPA2", "SFTPB", "NAPSA", "AGER", "CAV1", "PDPN", "SCGB3A2", "SCGB1A1"],
        neg=["PTPRC", "CD3D", "PECAM1", "COL1A1"],
        src=["Travaglini2020", "Peng2026"],
        src_crosscheck=["Habermann2020"],
        note="AT2=SFTPC/SFTPA1/SFTPA2/SFTPB/NAPSA；AT1=AGER/CAV1/PDPN；气道分泌=SCGB3A2/SCGB1A1；"
             "pan-上皮=EPCAM/KRT8/KRT18/KRT19/CDH1。NKX2-1(TTF-1) 为 LUAD 谱系决定因子（源论文）。",
    ),
    "T/NK": dict(
        genes=["CD3D", "CD3E", "CD3G", "TRAC", "CD4", "IL7R", "CD8A", "CD8B",
               "NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
        neg=["EPCAM", "COL1A1", "CD68"],
        src=["Travaglini2020"],
        src_crosscheck=["Habermann2020"],
        note="pan-T=CD3D/E/G+TRAC；CD4=CD4/IL7R；CD8=CD8A/CD8B；NK/细胞毒=NKG7/GNLY/KLRD1/PRF1/GZMB。"
             "🔴 原引 Guo2018 未核实其有 marker 表，已移入 UNVERIFIED_SRC。",
    ),
    "B/浆": dict(
        genes=["MS4A1", "CD19", "CD79A", "CD79B", "MZB1", "JCHAIN", "SDC1",
               "IGHG1", "IGKC", "XBP1", "DERL3"],
        neg=["CD3D", "EPCAM", "ACTA2"],
        src=["Travaglini2020"],
        src_crosscheck=["Habermann2020"],
        note="B=MS4A1/CD19/CD79A/CD79B；浆细胞=MZB1/JCHAIN/SDC1/IGHG1/IGKC/XBP1/DERL3。",
    ),
    "髓系": dict(
        genes=["LYZ", "AIF1", "ITGAX", "CD68", "CD163", "MSR1", "C1QA", "C1QB", "C1QC",
               "MARCO", "APOE", "FCN1", "CD14", "S100A8", "S100A9", "SPP1"],
        neg=["CD3D", "EPCAM", "PECAM1"],
        src=["Travaglini2020"],
        src_crosscheck=["Habermann2020"],
        note="pan-髓系=LYZ/AIF1/ITGAX；巨噬=CD68/CD163/MSR1/C1QA/B/C/MARCO/APOE；"
             "单核=FCN1/CD14/S100A8/S100A9；SPP1+ 巨噬=SPP1 🔴 出处偏弱，见 WEAK_PROVENANCE。"
             "🔴 原引 Zilionis2019 未核实其有 marker 表，已移入 UNVERIFIED_SRC。",
    ),
    "成纤维": dict(
        genes=["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "PDGFRA", "PDGFRB",
               "ACTA2", "TAGLN", "FAP", "CXCL12"],
        neg=["PTPRC", "EPCAM", "PECAM1"],
        src=["Travaglini2020"],
        src_crosscheck=["Habermann2020"],
        note="通用成纤维=COL1A1/A2/COL3A1/DCN/LUM/FN1；肌成纤维=ACTA2/TAGLN/FAP/PDGFRB；"
             "肺泡成纤维=PDGFRA/CXCL12。FAP 🔴 出处偏弱，见 WEAK_PROVENANCE。"
             "🔴 原引 Reyfman2019 / Lambrechts2018 未核实其有 marker 表，已移入 UNVERIFIED_SRC。",
    ),
    "内皮": dict(
        genes=["PECAM1", "CDH5", "KDR", "CD34", "VWF", "EGFL7", "RAMP2", "EMCN",
               "PLVAP", "AQP1", "CLDN5", "FLT1"],
        neg=["EPCAM", "PTPRC", "COL1A1"],
        src=["Travaglini2020"],
        src_crosscheck=["Habermann2020"],
        note="pan-内皮=PECAM1/CDH5/KDR/CD34/VWF；毛细血管 aCap/gCap 特化=EGFL7/RAMP2/EMCN/PLVAP/AQP1。"
             "KDR 🔴 出处偏弱，见 WEAK_PROVENANCE。"
             "🔴 原引 Gillich2020 / Lambrechts2018 未核实其有 marker 表，已移入 UNVERIFIED_SRC。",
    ),
}

LINEAGES = list(PANEL.keys())

# ---------------------------------------------------------------- CellTypist 61 类 → 6 谱系
# `Human_Lung_Atlas.pkl` 的真实出处是 Sikkema 2023 Nat Med 29:1563（HLCA 整合），
# **不是** PARAMETERS 原先登记的 Travaglini 2020 —— 2026-09-17 据模型自带元数据更正。
CELLTYPIST_MODEL = dict(
    path="~/.celltypist/data/models/Human_Lung_Atlas.pkl",
    cell_types=61,
    version="v2",
    source="Sikkema L, Ramírez-Suástegui C, Strobl DC, et al. An integrated cell atlas of the lung "
           "in health and disease. Nat Med 2023;29(6):1563-1577. doi:10.1038/s41591-023-02327-2",
    sha256=None,   # 运行时回填
)

# 映射为**显式登记表**，不是自动推断。标 AMBIG 的条目属判断而非共识，单列上报。
CELLTYPIST_TO_LINEAGE = {
    "AT0": "上皮", "AT1": "上皮", "AT2": "上皮", "AT2 proliferating": "上皮",
    "Basal resting": "上皮", "Club (nasal)": "上皮", "Club (non-nasal)": "上皮",
    "Deuterosomal": "上皮", "Goblet (bronchial)": "上皮", "Goblet (nasal)": "上皮",
    "Goblet (subsegmental)": "上皮", "Hillock-like": "上皮", "Ionocyte": "上皮",
    "Multiciliated (nasal)": "上皮", "Multiciliated (non-nasal)": "上皮",
    "Neuroendocrine": "上皮", "SMG duct": "上皮", "SMG mucous": "上皮",
    "SMG serous (bronchial)": "上皮", "SMG serous (nasal)": "上皮",
    "Suprabasal": "上皮", "Tuft": "上皮", "pre-TB secretory": "上皮",

    "CD4 T cells": "T/NK", "CD8 T cells": "T/NK", "NK cells": "T/NK",
    "T cells proliferating": "T/NK",

    "B cells": "B/浆", "Plasma cells": "B/浆",

    "Alveolar Mph CCL3+": "髓系", "Alveolar Mph MT-positive": "髓系",
    "Alveolar Mph proliferating": "髓系", "Alveolar macrophages": "髓系",
    "Classical monocytes": "髓系", "DC1": "髓系", "DC2": "髓系",
    "Interstitial Mph perivascular": "髓系", "Migratory DCs": "髓系",
    "Monocyte-derived Mph": "髓系", "Non-classical monocytes": "髓系",
    "Plasmacytoid DCs": "髓系",

    "Adventitial fibroblasts": "成纤维", "Alveolar fibroblasts": "成纤维",
    "Myofibroblasts": "成纤维", "Peribronchial fibroblasts": "成纤维",
    "Subpleural fibroblasts": "成纤维",

    "EC aerocyte capillary": "内皮", "EC arterial": "内皮", "EC general capillary": "内皮",
    "EC venous pulmonary": "内皮", "EC venous systemic": "内皮",
    "Lymphatic EC differentiating": "内皮", "Lymphatic EC mature": "内皮",
    "Lymphatic EC proliferating": "内皮",
}

# 判断项（非共识）—— 单列上报，不混入主 κ
CELLTYPIST_AMBIGUOUS = {
    "Hematopoietic stem cells": ("髓系", "HSC 是否归髓系无共识；HLCA 单列"),
    "Mast cells": ("髓系", "肥大细胞属髓系起源但常单列；HLCA 单列"),
    "Mesothelium": ("未归属", "间皮不属六大谱系任何一个；保持未归属"),
    "Pericytes": ("成纤维", "血管周细胞属 mural 而非成纤维；本项目归入基质"),
    "SM activated stress response": ("成纤维", "平滑肌应激态归入基质"),
    "Smooth muscle": ("成纤维", "平滑肌归入基质（非成纤维，判断项）"),
    "Smooth muscle FAM83D+": ("成纤维", "同上"),
}


def resolve_missing(available_genes):
    """返回 (可用面板, 缺失清单)。缺失按谱系列出，不静默丢弃。"""
    have = set(available_genes)
    used, missing = {}, {}
    for lin, d in PANEL.items():
        used[lin] = [g for g in d["genes"] if g in have]
        miss = [g for g in d["genes"] if g not in have]
        if miss:
            missing[lin] = miss
    return used, missing


def lineage_map_full():
    """完整 61 类映射（含判断项），用于审计打印。"""
    m = dict(CELLTYPIST_TO_LINEAGE)
    for k, (v, _) in CELLTYPIST_AMBIGUOUS.items():
        m[k] = v
    return m


def assert_provenance():
    """硬校验：`src` 里不许出现未核实出处的文献；弱出处基因必须被登记。"""
    bad = [(lin, s) for lin, d in PANEL.items() for s in d["src"] if s in UNVERIFIED_SRC]
    if bad:
        raise AssertionError(f"src 混入未核实出处：{bad}")
    unknown = [(lin, s) for lin, d in PANEL.items()
               for s in list(d["src"]) + list(d.get("src_crosscheck", [])) if s not in REFERENCES]
    if unknown:
        raise AssertionError(f"src 引用了 REFERENCES 里没有的键：{unknown}")
    allg = {g for d in PANEL.values() for g in d["genes"]}
    notmark = [g for g in WEAK_PROVENANCE if g not in allg]
    if notmark:
        raise AssertionError(f"WEAK_PROVENANCE 登记了面板里没有的基因：{notmark}")
    return True


if __name__ == "__main__":
    assert_provenance()
    n = sum(len(d["genes"]) for d in PANEL.values())
    print(f"面板：{len(PANEL)} 谱系，阳性 marker {n} 个（去重前）")
    for lin, d in PANEL.items():
        print(f"  {lin:6s} {len(d['genes']):2d} 个  "
              f"出处 {d['src']}  第二意见 {d.get('src_crosscheck', [])}")
    print(f"\n未核实出处的文献（登记备查，不在 src 内）：{list(UNVERIFIED_SRC)}")
    print(f"弱出处基因（保留，须随结果报告）：{list(WEAK_PROVENANCE)}")
    mm = lineage_map_full()
    print(f"\nCellTypist 映射覆盖 {len(mm)}/61 类；未覆盖 {61 - len(mm)} 类")
    missing_ct = [c for c in mm if c not in CELLTYPIST_TO_LINEAGE]
    print("其中属判断项:", missing_ct)
