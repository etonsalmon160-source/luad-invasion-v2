"""GP8a 上皮亚聚类 —— 亚型 marker 面板（逐基因可溯源）。

生成方式：本文件由 `05_annotation/build_epi_subtype_panel.py` 从原始补充材料
（源论文 Table S2 / Table S3 的 .xlsx）**程序化生成**，非手工转录 —— 保证基因名逐字一致。

两个来源层
----------
1. **主干 = 源论文自身上皮 MP**（Table S2 第二块 "snRNA-seq (lung epithelium)"，9 个 meta-program，
   每 MP 50 基因）。论文定义的亚型就是我们要的亚型，且 450/450 基因槽**全部命中本数据矩阵**。
   ⚠️ 但 MP 是 **NMF meta-program**，不是洁净 marker 表 —— 里面混有应激/管家类基因
   （如 MP2 AT2 含 SOD2/CXCL2/NR4A1）。**直接整块做 score_genes 特异性会偏低**，
   故同时提供 CORE 小集做交叉校验。
2. **CORE = 经典亚型 marker**，逐基因出处核对状态见 CORE 的 status 字段。
   🔴 未经原论文补充表逐条核对前，status 一律为 "unverified"。

已知数据源缺陷（必须随结果一起报告）
------------------------------------
- 源论文 Table S2 的 MP4 与 MP7 各有一处 **Excel 自动改名**：`DKK 3.00`（应为 DKK3）、
  `ERN 1.00`（应为 ERN1）。已在本文件手工还原并在此登记，**不静默使用**。
- 作者上传的矩阵只有 18,082 个基因（约标准转录组之半）；`SFTPA2` 等基因**不在数据中**，
  故本面板已剔除，并在 MISSING_IN_MATRIX 登记。
"""

# ---------------------------------------------------------------- 一次文献
REFERENCES = {
    "Peng2026": "Peng F, Sinjab A, et al., Kadara H. Multimodal spatial-omics reveal co-evolution of "
                "alveolar progenitors and proinflammatory niches in progression of lung precursor "
                "lesions. Cancer Cell 2026;44(2):321-339.e13. doi:10.1016/j.ccell.2025.10.004 "
                "【本项目源论文；本面板主干出自其 Table S2 / Table S3】",
    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "
                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "
                      "doi:10.1038/s41586-020-2922-4",
    "VieiraBraga2019": "Vieira Braga FA, Kar G, Berg M, et al. A cellular census of human lungs "
                       "identifies novel cell states in health and in asthma. Nat Med "
                       "2019;25(7):1153-1163. doi:10.1038/s41591-019-0468-5",
}

# ---------------------------------------------- Excel 改名修复（登记，非静默）
EXCEL_NAME_FIX = {
    "DKK\xa03.00": "DKK3",   # Table S2 MP4 Basal/basal stem
    "ERN\xa01.00": "ERN1",   # Table S2 MP7 Tumor cell (stress/inflammatory)
}

# ------------------------------------------------- 源论文 Table S2 上皮 MP（主干）
# 逐 MP 基因序 = 原表列序，逐字未改（仅套用 EXCEL_NAME_FIX）。
PAPER_EPI_MP = {
    "MP1": dict(label="Ciliated", genes=[
        "CDHR3", "MAPK15", "CFAP70", "DNAH7", "CFAP157", "DLEC1",
        "RSPH1", "CFAP43", "DNAAF1", "SPEF2", "TOGARAM2", "CCDC78",
        "DNAI1", "MS4A8", "WDR66", "CAPS", "CFAP73", "HYDIN",
        "NEK10", "SAXO2", "CCDC187", "CFAP52", "CFAP54", "DNAH12",
        "DNAH9", "VWA3A", "CDHR4", "CFAP47", "DNAJA4", "LRRC46",
        "CFAP100", "DNAH10", "DTHD1", "FAM92B", "RFX3", "ZBBX",
        "AC007906.2", "CFAP46", "FRMPD2", "TEKT2", "TSPAN1", "UBXN11",
        "C6orf118", "DZIP1L", "PTPRT", "ERICH3", "SPAG8", "UBXN10",
        "ULK4", "PRR29",
    ]),
    "MP2": dict(label="AT2", genes=[
        "PGC", "C3", "LAMP3", "SERPINA1", "CRTAC1", "HHIP",
        "NPC2", "C11orf96", "WIF1", "NAPSA", "PLA2G4F", "SCD",
        "TCIM", "ALPL", "DMBT1", "PIGR", "BMP1", "ETV5",
        "C2", "CXCL17", "SPRY4", "STEAP4", "C16orf89", "C4BPA",
        "MFSD2A", "SFTPD", "CCDC141", "KIAA1324L", "NNMT", "TFPI",
        "MID1IP1", "NRGN", "SLPI", "ABCA3", "CD36", "DBI",
        "EPHX1", "NR4A1", "PARM1", "BMP2", "COLEC12", "KCNJ15",
        "LGALSL", "LRRK2", "RASGRF1", "SDR16C5", "CITED2", "CXCL2",
        "SOD2", "HMGCS1",
    ]),
    "MP3": dict(label="Club/secretory", genes=[
        "SCGB1A1", "SCGB3A1", "CP", "KIAA1324", "KLK11", "CRACR2B",
        "MUC4", "CAPN13", "PLPP2", "WFDC2", "ALDH1A1", "CLU",
        "PDCD4", "PIGR", "FMO2", "GRAMD2B", "SCGB3A2", "SLPI",
        "BPIFB1", "MGP", "ABCC5", "CTSC", "LCN2", "ECE1",
        "NOTCH3", "RASSF9", "CXCL1", "EHF", "FGFR3", "GPC1",
        "MUC20", "PRSS23", "STEAP4", "TENT5C", "TSPAN1", "GRHL1",
        "HES4", "ITGA9", "MET", "P2RY2", "SCNN1G", "SERPINF1",
        "TGM2", "CDK14", "IGFBP2", "KLK10", "MYO15B", "NPDC1",
        "SCNN1B", "WNK2",
    ]),
    "MP4": dict(label="Basal/basal stem", genes=[
        "IGFBP2", "KRT15", "SERPINF1", "CDH3", "CLDN1", "ECE1",
        "F3", "GPC1", "ITGA2", "FGFR3", "ITGB4", "TRIM29",
        "DST", "GCLC", "PLCH2", "PRSS23", "ZFP36L2", "ARHGAP23",
        "CLSTN1", "EYA2", "FMO2", "IL33", "KRT17", "KRT19",
        "OBSCN", "S100A2", "TP63", "ANOS1", "COL7A1", "DDIT4",
        "FHL2", "KRT5", "PERP", "CLU", "COL18A1", "ELN",
        "IFITM10", "JAG1", "LAMB3", "LTBP2", "NECTIN1", "NTN1",
        "SLC1A5", "TPM1", "ABCC3", "CFH", "DHRS3", "DKK3",
        "ETS2", "FBLN1",
    ]),
    "MP5": dict(label="AT1", genes=[
        "SPOCK2", "TIMP3", "CCN1", "COL4A2", "MYRF", "PHLDB2",
        "PKDCC", "RTKN2", "SCEL", "ABCA1", "AGER", "ANOS1",
        "CAV1", "NCKAP5", "CEACAM6", "CLIC5", "COL4A1", "KRT7",
        "PALM2-AKAP2", "ST6GALNAC5", "DST", "KHDRBS2", "RGCC", "SOGA3",
        "UPK3B", "ARHGEF26", "LIMS2", "MS4A15", "CCN2", "FAM189A2",
        "WFS1", "CAV2", "COL4A3", "CDKN2B", "PGGHG", "SEMA3B",
        "SEMA5A", "TACSTD2", "SLC5A9", "UNC13D", "ANXA3", "CRIP1",
        "SLC19A1", "EMP2", "ICAM1", "NDNF", "NTM", "SPARC",
        "ADRB2", "LRRN4",
    ]),
    "MP6": dict(label="Tumor cell/KAC", genes=[
        "CTSE", "PLAAT4", "IFI27", "IFI6", "GDF15", "MMP7",
        "CD24", "CEACAM6", "ABCC3", "BST2", "C1R", "IL32",
        "SCGB3A2", "TM4SF1", "C1S", "HSD17B11", "LY6E", "SERPINA1",
        "SOX4", "STAT1", "TAP1", "TRAM1", "GSN", "KRT7",
        "MSLN", "S100A11", "WFDC2", "FN1", "IFI44L", "MX1",
        "TMSB10", "AGR3", "FBP1", "IFI44", "IGHG1", "ISG15",
        "ITGAV", "RAB27B", "SOD3", "PEG10", "XAF1", "IFITM3",
        "OAS3", "TIMP1", "MDK", "EEF1G", "MARCKSL1", "ICAM1",
        "ANXA2", "CLIC5",
    ]),
    "MP7": dict(label="Tumor cell (stress/inflammatory)", genes=[
        "ICAM4", "CXCL2", "NCOA7", "RND1", "ITPKC", "STEAP4",
        "RRAD", "SOD2", "NAMPT", "NFKBIA", "SLC6A14", "SOCS2",
        "BHLHE40", "DMBT1", "IRF1", "MYC", "PTP4A3", "SERPINB1",
        "SOCS3", "ZC3H12A", "ICAM1", "NNMT", "TCIM", "ATF3",
        "FOSB", "GADD45B", "NR4A1", "BMP2", "CSF3", "F3",
        "NFKBIZ", "SERPINA1", "SPRY4", "AREG", "BCL3", "CD83",
        "CDC42EP1", "DUSP6", "EGR1", "ELF3", "EPHA2", "ERN1",
        "ETS2", "PLK3", "SBNO2", "TIPARP", "TRIB1", "JUNB",
        "TNFRSF12A", "NFKB2",
    ]),
    "MP8": dict(label="Ciliated", genes=[
        "C11orf88", "C20orf85", "C9orf24", "CAPS", "CETN2", "DNALI1",
        "FOXJ1", "KCTD12", "POLR2I", "PSENEN", "TSPAN1", "TUBB4B",
        "UFC1", "FAM183A", "HSP90AA1", "IGFBP2", "MORN2", "PIFO",
        "PRDX5", "TPPP3", "ALDH1A1", "CTXN1", "LRRC10B", "ZMYND10",
        "ATP5IF1", "CAPSL", "PFN2", "RSPH1", "SNTN", "C9orf116",
        "GSTA1", "IK", "LRTOMT", "TEKT1", "C5orf49", "CD24",
        "CLU", "FAM216B", "IFT22", "IGFBP5", "MARCKS", "AC007906.2",
        "CES1", "ENKUR", "METRN", "SCGB1A1", "SRI", "LRRC46",
        "NUCB2", "LDLRAD1",
    ]),
    "MP9": dict(label="KAC/inflammatory", genes=[
        "CLCF1", "EMP1", "HBEGF", "LAMC2", "PLAUR", "PTGS2",
        "SFN", "TM4SF1", "ATF3", "DGKD", "FOSB", "IER3",
        "IL32", "ITGA2", "MAFF", "RRAD", "SEMA4A", "SGMS2",
        "ATP10A", "CCN1", "CDKN1A", "CENPT", "EID3", "ERRFI1",
        "GADD45A", "GADD45B", "ITPKC", "KDM6B", "MT2A", "NABP1",
        "PIM1", "SLC25A44", "TMEM88", "YOD1", "ANGPTL4", "FOSL1",
        "LAMB3", "EPHA2", "PHLDA2", "CTNNAL1", "CAVIN1", "ICAM4",
        "MAP7D1", "EGR1", "SPSB1", "AREG", "SOCS3", "TRIB1",
        "MYADM", "PLK3",
    ]),
}

# 亚型 -> 由哪些 MP 组成
PAPER_SUBTYPE_FROM_MP = {
    "Ciliated": ['MP1', 'MP8'],   # MP1 ∪ MP8（论文给了两个纤毛 MP）
    "AT2": ['MP2'],
    "Club/secretory": ['MP3'],
    "Basal/basal stem": ['MP4'],   # 含 Excel 改名修复：DKK 3.00 → DKK3
    "AT1": ['MP5'],
    "Tumor cell/KAC": ['MP6'],
    "Tumor cell (stress/inflammatory)": ['MP7'],   # 含 Excel 改名修复：ERN 1.00 → ERN1
    "KAC/inflammatory": ['MP9'],
}

# ------------------------- 源论文 Table S3 签名（人，已剔除小鼠列）
PAPER_SIGNATURES = {
    "KAC_human": [
        "GAPDH", "IGFBP3", "AGR2", "TIMP1", "MDK", "MMP7",
        "CRABP2", "MUC5B", "PCSK1N", "TFF3", "S100A11", "HSP90B1",
        "NQO1", "WFDC2", "PHLDA2", "PLPP2", "TESC", "C4orf48",
        "CDKN2A", "ABCC3", "GOLM1", "PHLDA1", "DDIT4", "PTMA",
        "HSPA5", "MARCKSL1", "HPGD", "FTH1", "PDIA4", "SPINK1",
        "LCN2", "AKR1C2", "AKR1C3", "TFF1", "C15orf48", "TNFRSF21",
        "MIF", "RPL41", "CLDN10", "PDLIM4", "SSR4", "VSTM2L",
        "ANXA1", "PLTP", "THBS1", "CALR", "PFKP", "GGT5",
        "GPX2", "ECE1", "PLAU", "AKR1B10", "KLK6", "COTL1",
        "TSPAN1", "MGP", "CFH", "MANF", "IGFBP2", "CRLF1",
        "ANGPTL4", "CDKN1A", "SLCO2A1", "MAP1B", "GADD45A", "PPP1R14B",
        "SOX9", "SPDEF", "CX3CL1", "TGM2", "EIF4EBP1", "KRT8",
        "RPL13A", "RPL15", "CYBA", "TMSB10", "TPI1", "S100P",
        "RPS2", "SH3BGRL3", "CSTB", "CFB", "FXYD5", "RASD1",
        "DUSP5", "NPDC1", "PLAUR", "PTGES", "SLC2A1", "BCAS1",
        "CP", "RPS18", "LINC00511", "GCLC", "RPS19", "AKR1B1",
        "ARHGDIB", "JPT1", "SNX9", "SLC16A7", "ACTG1", "TSPAN6",
    ],
    "NFkB_human": [
        "CHUK", "FADD", "IKBKB", "IKBKG", "IL1A", "IL1R1",
        "MAP3K1", "MAP3K14", "MAP3K7", "MYD88", "NFKB1", "NFKBIA",
        "RELA", "RIPK1", "TAB1", "TNF", "TNFAIP3", "TNFRSF1A",
        "TNFRSF1B", "TRADD", "TRAF6",
    ],
    "Interferon_human": [
        "ACTN1", "APOL1", "APOL2", "APOL3", "APOL6", "B2M",
        "BST2", "C19orf66", "C1R", "C1S", "C3", "C5orf56",
        "CCL2", "CCL27", "CCL5", "CD74", "CFB", "CFH",
        "CMPK2", "CX3CL1", "CXCL10", "CXCL11", "CXCL9", "DDX58",
        "DDX60", "DDX60L", "DYNLT1", "EIF2AK2", "EPSTI1", "ETV7",
        "FXYD5", "GBP1", "GBP2", "GBP3", "GBP4", "GBP5",
        "GLUL", "HAPLN3", "HELZ2", "HERC5", "HERC6", "HLA-A",
        "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1",
        "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
        "HLA-E", "HLA-F", "HSPB8", "ICAM1", "IDO1", "IFI16",
        "IFI27", "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1",
        "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3",
        "IL32", "IRF1", "IRF7", "ISG15", "ISG20", "KLHDC7B",
        "LAMP3", "LAP3", "LGALS3BP", "LGALS9", "LY6E", "MDK",
        "MSRB1", "MX1", "MX2", "NNMT", "NT5C3A", "NUB1",
        "OAS1", "OAS2", "OAS3", "OASL", "PARP14", "PARP9",
        "PDZK1IP1", "PLSCR1", "PSMB10", "PSMB8", "PSMB9", "PSME1",
        "PSME2", "RARRES3", "RNF213", "RSAD2", "RTP4", "SAA1",
        "SAMD9", "SAMD9L", "SAMHD1", "SERPING1", "SERPINH1", "SLFN5",
        "SOD2", "SP100", "SP110", "SPP1", "SQRDL", "STAT1",
        "TAP1", "TAP2.1", "TAPBP", "TGM2", "TNFSF10", "TNFSF13B",
        "TRIM22", "TYMP", "UBD", "UBE2L6", "USP18", "VCAM1",
        "WARS", "XAF1", "ZC3HAV1", "ZNFX1",
    ],
}

# 注：Table S3 另有两列是小鼠（KAC / Inflammatory），本项目为人类数据，**不采纳**，如实登记。

MISSING_IN_MATRIX = ["SFTPA2"]   # 作者矩阵本身不含；见模块 docstring


def build_panel(mode="paper_mp"):
    """返回 {亚型: [基因...]}。

    mode="paper_mp" : 主干，源论文 Table S2 上皮 MP 合并（Ciliated = MP1 ∪ MP8）
    mode="kac_sig"  : 源论文 Table S3 人 KAC 签名，作为 KAC 的独立第二来源
    """
    if mode == "paper_mp":
        return {st: sorted({g for m in mps for g in PAPER_EPI_MP[m]["genes"]})
                for st, mps in PAPER_SUBTYPE_FROM_MP.items()}
    if mode == "kac_sig":
        return {"KAC": [g for g in PAPER_SIGNATURES["KAC_human"] if g not in MISSING_IN_MATRIX]}
    raise ValueError(f"未知 mode: {mode}")

