"""🔴 本文件由 `05_annotation/build_panels_s4.py` 程序化生成，**不要手改**。

来源：Travaglini 2020 Nature 587:619，**Table S4**（per-cluster enriched markers）。
    仍是**同一篇论文**，举证只需核一张表；从 S1 的「Canonical markers」列换成 S4 的
    逐簇富集表，是因为 S1 那一列在血管内皮上不是单细胞证据（见生成器文件头）。

口径（用户 2026-09-22 签字）：
    ① 只用 10x sheet（排除 SS2/Smart-seq2）
    ② 过筛 avg_logFC ≥ 1.0 且 (pct_in − pct_out) ≥ 0.3
    ③ 按 logFC 降序，运行时取前 20 个**本矩阵里真有的**（矩阵缺则顺延补齐）
    ④ 亚型粒度照 S4 原粒度
    ⑤ S4 无 10x 表的型沿用 Table S1 行（逐条登记于 FALLBACK_FROM_S1）
    ⑥ 不做别名修正（旧符号不在矩阵者如实报 missing）

🔴 面板**只用于亚型（L2）注释**。六谱系的上位归属（L1）仍由 marker_panel.py 决定，
   不受本文件影响。
"""

SOURCE = '/home/eto/luad_v2/data/external/travaglini2020/TableS4_cluster_enriched_markers.xlsx'
SOURCE_SHA256 = 'abca9f72c5bb8d47327f543738ae9672c1bf006c92a6c0eb4b4cc2dc29031073'
SOURCE_NOTE = ('Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. doi:10.1038/s41586-020-2922-4；PMC7704697【Table S4】')

SELECTION_LOGFC = 1.0
SELECTION_SPEC = 0.3
PANEL_CAP = 20

# S4 的 10x 表缺失、沿用 Table S1 的型（逐条登记，不静默）
FALLBACK_FROM_S1 = {'髓系': {'Neutrophil': ['S100A8', 'S100A9', 'IFITM2', 'FCGR3B'], 'Eosinophil': ['SIGLEC8']}}
FALLBACK_REASON = {'Neutrophil': 'Table S4 里中性粒**只有 SS2 表**（Cluster 43 (SS2)），无 10x 表。严格走「只用 10x」会丢掉该型，而它在本数据里标注了 1,406 个核 ⇒ 沿用 S1 行。', 'Eosinophil': 'Table S4 **没有嗜酸簇**（该型稀少未被单独分出）⇒ 沿用 S1 行（其本来即 1 槽，不可打分）。'}

SOURCE_TABLE_DEFECTS = {'Cluster 21 标签不一致': '10x 表第 0 行写「Capillary Intermediate 2」，同簇 SS2 表写「Capillary Intermediate 1」。本项目只用 10x ⇒ 取「Capillary Intermediate 2」。缺陷原样登记，不修（改了就成了我们的判断）。', 'Cluster 50 标签不一致': '10x 表写「Myeloid Dendritic Type 1」，SS2 表写「Dendritic」。同上，取 10x 的标签。', '旧基因符号': 'S4 混用旧符号（CTGF/CYR61/SEPP1/KIAA0101/C10orf54/FAM26F/NAPSB 等）。本项目不做别名修正（见文件头 ⑥），不在矩阵者如实报 missing。'}

UNANNOTATABLE = {}   # S4 里每型都有 marker；原 S1 的 Bronchial Vessel 在 S4 中已有基因
UNSCORABLE_FROM_S4 = ['Capillary Intermediate 2']

RARE = ['Platelet/Megakaryocyte', 'Eosinophil', 'Bronchial Vessel 1', 'Bronchial Vessel 2',
        'Capillary Intermediate 1', 'Capillary Intermediate 2']

CONTAM_CHECK = ['PTPRC', 'CD3D', 'COL1A1', 'PECAM1', 'CD68', 'EPCAM']

# 过筛后的**完整**名单（按 logFC 降序）。运行时的封顶与矩阵存在性检查在 build_panel。
PANELS = {
    '内皮': {
        'Artery': ['CXCL12', 'IGFBP3', 'GJA5', 'DKK2', 'CXCL2', 'MGP', 'ENPP2', 'HEY1', 'SOX17', 'SERPINE2', 'IDO1', 'CTNNAL1', 'ARGLU1'],
        'Vein': ['CPE', 'ACKR1', 'PTGDS', 'C7', 'CLU', 'PRSS23', 'VWF', 'VCAM1', 'CCL23', 'CYP1B1', 'PLAT', 'PLA1A', 'IL1R1', 'PTGIS', 'RGS5', 'NCOA7', 'SRPX', 'ABI3BP', 'NNMT', 'LIFR', 'CXCL2'],
        'Capillary Aerocyte': ['SOSTDC1', 'HPGD', 'EDNRB', 'S100A3', 'IL1RL1', 'EMCN', 'S100A4', 'TSPAN12', 'FRY', 'CYP3A5', 'B3GALNT1'],
        'Capillary': ['IL7R', 'FCN3', 'SLC6A4', 'EDN1'],
        'Capillary Intermediate 1': ['IL1RL1', 'SERPINE1', 'EDNRB'],
        'Capillary Intermediate 2': [],
        'Bronchial Vessel 1': ['SPRY1', 'PLVAP', 'CTGF', 'ACKR1', 'MYC', 'CYR61', 'ADAMTS1', 'POSTN', 'SNHG5'],
        'Bronchial Vessel 2': ['SLC2A3', 'CYR61', 'EMP1', 'HBEGF', 'ABCB1'],
        'Lymphatic': ['CCL21', 'TFF3', 'FABP4', 'MMRN1', 'EFEMP1', 'RGS16', 'IGFBP5', 'PDPN', 'FABP5', 'NTS', 'HSPB1', 'APOD', 'LYVE1', 'SNCG', 'IGF1', 'CD9', 'ARL4A', 'PPP1R2', 'PDLIM4', 'STMN1', 'PROX1', 'GYPC', 'ARID5B', 'LIMS1', 'FSCN1', 'MRC1', 'LINC00152', 'TSTA3', 'RARRES2', 'OAF', 'TUBB6', 'LAPTM5', 'DKK3', 'LRRC70', 'ANGPT2', 'GGTA1P'],
    },
    '髓系': {
        'Basophil/Mast 1': ['TPSB2', 'TPSAB1', 'CPA3', 'MS4A2', 'HPGDS', 'VWA5A', 'LTC4S', 'SLC18A2', 'GATA2', 'CLU', 'RGS13', 'KIT', 'JUN', 'MAOB', 'HPGD', 'LAPTM4A', 'C1orf186', 'RGS2', 'RGS1', 'CD69', 'CD9', 'SAMSN1', 'IL1RL1', 'SDPR', 'HDC', 'TDRD3', 'CALB2', 'SIGLEC17P', 'PTGS1', 'LOC284454', 'FCER1A', 'EGR1', 'CPM', 'CTNNBL1', 'PEBP1'],
        'Basophil/Mast 2': ['LMNA', 'TPSAB1', 'TPSB2', 'CPA3', 'TSC22D1', 'MS4A2', 'BIRC3', 'HPGDS', 'LOC101927482', 'GADD45B', 'CCL4L1', 'CSF1', 'RGS13', 'SELK', 'CD69', 'IL1RL1', 'RGS2', 'HPGD', 'VWA5A', 'CPM', 'HDC', 'TNFRSF9', 'GLUL', 'GPR65', 'CREM', 'AREG', 'ID2', 'PLIN2', 'NR4A1', 'FOSB', 'C1orf186', 'CLU', 'LAPTM4A', 'GATA2', 'SLC18A2', 'SDCBP', 'PTGS2', 'RGS1', 'CD83', 'KIT', 'LTC4S', 'LEO1', 'SYAP1', 'LMO4', 'CD9', 'NR4A3', 'SIGLEC17P', 'MAOB', 'SAMSN1', 'BHLHE40', 'TNFAIP3', 'KDM6B', 'MYADM', 'CALB2', 'REL', 'SDPR', 'NFKBIZ', 'SLC2A3', 'METTL21A', 'BRE-AS1', 'DUSP6', 'SOCS3', 'TIMP3', 'LOC284454', 'SLC45A3', 'ELL2', 'CKLF', 'SGK1', 'DUSP10', 'ELF1', 'SLC26A2', 'RHOH', 'SOCS1', 'FCER1A', 'ICAM1'],
        'Platelet/Megakaryocyte': ['PPBP', 'PF4', 'GNG11', 'TUBB1', 'ACRBP', 'RGS18', 'SDPR', 'PTCRA', 'GP9', 'NRGN', 'CMTM5', 'TMEM40', 'SPARC', 'MMD', 'NCOA4', 'TREML1', 'CA2', 'TUBA4A', 'MAP3K7CL', 'HBB', 'RUFY1', 'MYL9', 'GRAP2', 'LIMS1', 'MAX', 'NT5C3A', 'TSC22D1', 'ODC1', 'PGRMC1', 'CTTN', 'RGS10', 'MPP1', 'DMTN', 'GPX1', 'RAP1B', 'CTSA', 'CLU', 'ITGA2B', 'TPM4', 'C2orf88', 'ESAM', 'FERMT3', 'NGFRAP1', 'HRAT92', 'MARCH2', 'KIF2A', 'PF4V1', 'GSTO1', 'PLA2G12A', 'CLEC1B', 'CXCR2P1', 'F13A1', 'RNF11', 'GAS2L1', 'YWHAH', 'ACTN1', 'PDZK1IP1', 'CLDN5', 'GMPR', 'TLN1', 'TPM1', 'PDLIM1', 'SMOX', 'DAPP1', 'GP1BA', 'C9orf89', 'C19orf33', 'LCN2', 'RAB11A', 'C6orf25', 'SNCA', 'ETFA', 'PTPN18', 'TLK1', 'DAB2', 'MYLK', 'VIM-AS1', 'NFE2', 'ILK', 'TPST2', 'CALM3', 'HEMGN', 'EIF2AK1', 'PTGS1', 'SNAP23', 'TUBA8', 'VCL', 'HIST1H2BK', 'SLC40A1', 'PRKAR2B', 'LYL1', 'RSU1', 'HIST2H2BE', 'HIST1H2AC', 'GFI1B', 'SMIM5', 'TNFSF4', 'LAT', 'R3HDM4', 'AMD1', 'STOM', 'ABCC3', 'STON2', 'GADD45A', 'TRIM58', 'SH3BGRL2', 'LY6G6F', 'RAB27B', 'FAM110A', 'LDLRAP1', 'RAB32', 'H3F3AP4', 'GGTA1P', 'CDKN2D', 'NAT8B', 'TBXA2R', 'ENKUR', 'NEXN', 'CAPZA2', 'LAMTOR1', 'SCGB1C1', 'SCN1B', 'FRMD3', 'ARHGAP6', 'EMC3', 'HEXIM2', 'DNAJB6', 'AKIRIN2', 'SMIM3', 'TMBIM1', 'VDAC3', 'CCDC85B', 'PARVB', 'PPM1A', 'GNAZ', 'CNST', 'CAPN1', 'ALOX12', 'HIST1H1C', 'XPNPEP1', 'HGD', 'UBL4A', 'PKM', 'MFSD1', 'PRDX6', 'TUBA1C', 'INF2', 'SPINT2', 'PDCD10', 'LYPLAL1', 'APP', 'MIR4435-2HG', 'NUTF2', 'SWI5', 'STX11', 'BIN2', 'PNMA1', 'PLEK', 'LINC00152', 'TMEM91', 'GRHL1', 'TACC3', 'WBP2', 'SNN', 'WDR1', 'RBBP6', 'ITGB1', 'ICAM2', 'TMEM55A', 'HACD4', 'TNNC2', 'PPP1R14A', 'MTURN', 'SPNS1', 'WBP5', 'TGFB1', 'MLH3', 'GTPBP2', 'RBX1', 'PYGL', 'CORO1C', 'RABGAP1L', 'FHL1', 'PIP4K2A', 'SLA2', 'GUCY1B3', 'LGALSL', 'IFRD1', 'SYMPK', 'MGLL', 'FAM63A', 'MIR6843', 'ABHD16A', 'SENCR', 'ITGB3', 'YPEL5', 'CD9', 'SELP', 'CMIP', 'H2AFJ', 'RDH11', 'PBX1', 'AHCTF1', 'FAM212A', 'ANKRD9', 'AIG1', 'DYNLL1', 'HMG20B', 'F11R', 'RPA1', 'LOC113230', 'CCL5', 'TBC1D20', 'SIAH2', 'SLC39A3', 'PDE5A', 'ARRB1', 'FRMD4B', 'UBE2E3', 'LTBP1', 'BCL2L1', 'ANAPC5', 'MPST', 'H1F0', 'TSC22D4', 'TALDO1', 'ZNF185', 'LINC00657', 'DIAPH1', 'PROS1', 'TUBB4B', 'DERA', 'MYH9', 'LEPROT', 'ANO6', 'TRAPPC5', 'TAL1', 'ZYX', 'STXBP2', 'GNAS', 'TDRP', 'ATP2A3', 'ABLIM3', 'RTN3', 'ABI1'],
        'Macrophage': ['FABP4', 'APOC1', 'CCL20', 'C1QB', 'MARCO', 'C1QA', 'MCEMP1', 'C1QC', 'ALDH2', 'FABP5', 'LGALS3', 'FBP1', 'VSIG4', 'FN1', 'ACP5', 'IFI27', 'CES1', 'SERPING1', 'GRN', 'CTSD', 'GCHFR', 'MSR1', 'SERPINA1', 'CD68', 'APOE', 'CCL18', 'RETN', 'OLR1', 'CXCL5', 'MT1X', 'LGALS3BP', 'IFI30', 'MGST3', 'ANXA2', 'CXCL3', 'TREM1', 'CYP27A1', 'IFI6', 'LPL', 'MRC1', 'CSTB', 'RGCC', 'NUPR1', 'ANXA5', 'CTSC', 'TSPO', 'GSTO1', 'SNX10', 'PDLIM1', 'HP', 'HSPB1', 'GPNMB', 'GLIPR2', 'CTSL', 'ALOX5AP', 'LTA4H', 'HLA-DRB1', 'CD9', 'STXBP2', 'RBP4', 'TXN', 'MS4A7', 'LOC731424', 'PHLDA3', 'HLA-DRB5', 'FCGRT', 'CTSZ', 'S100A13', 'PRDX1', 'MS4A4A', 'CAPG'],
        'Proliferating Macrophage': ['FABP4', 'TUBA1B', 'MARCO', 'H2AFZ', 'C1QB', 'APOC1', 'STMN1', 'C1QA', 'TUBB', 'GRN', 'FN1', 'C1QC', 'UBE2C', 'VSIG4', 'CES1', 'ALDH2', 'KIAA0101', 'LGALS3', 'CCL18', 'HMGN2', 'HP', 'CCL20', 'PTTG1', 'FABP5', 'SERPING1', 'TK1', 'SERPINA1', 'LGALS3BP', 'FBP1', 'ARL6IP1', 'CDK1', 'MCEMP1', 'MGST1', 'TYMS', 'PCNA', 'ACP5', 'CTSD', 'TUBB4B', 'CKS2', 'MS4A4A', 'MGST3', 'CKS1B', 'CTSC', 'TREM1', 'RETN', 'IFI6', 'TUBA1C', 'RBP4', 'PDLIM1', 'CD68', 'NUSAP1', 'HIST1H4C', 'GSTO1', 'APOE', 'CYP27A1', 'SNX10', 'CXCL3', 'ANXA2', 'LRPAP1', 'GCHFR', 'GPNMB', 'HEXB', 'CENPW', 'TOP2A', 'GGH', 'LYZ', 'RHEB', 'KPNA2', 'CD59', 'ANXA5', 'RRM2', 'IFI27', 'LPL', 'MT1G', 'TUBB6', 'HLA-DRB1', 'CDKN3', 'CTSL', 'HN1', 'CENPF', 'CXCL5', 'LOC731424', 'PLBD1', 'TSPO', 'BIRC5', 'MS4A7', 'MT2A'],
        'Plasmacytoid Dendritic': ['IRF7', 'GZMB', 'PLD4', 'LILRA4', 'PPP1R14B', 'IRF8', 'PTGDS', 'SERPINF1', 'TSPAN13', 'CLIC3', 'ITM2C', 'CLN8', 'SCT', 'CXCR3', 'TCF4', 'PLAC8', 'NAPSB', 'C12orf75', 'UGCG', 'PTPRS', 'TCL1A', 'LRRC26', 'CCDC50', 'GPR183', 'IL3RA', 'SEC61B', 'SPIB', 'PTCRA', 'IRF4', 'TPM2', 'MAP1A', 'APP', 'SMPD3', 'LILRB4', 'DERL3', 'BCL11A', 'TRAF4', 'HERPUD1', 'SOX4', 'CYB561A3', 'CLEC4C', 'RNASE6', 'OPN3', 'DNASE1L3', 'GAS6', 'RASD1', 'LAMP5', 'FUT7', 'RNASET2', 'MYBL2', 'MZB1', 'MPEG1', 'C12orf45', 'PRMT9', 'PLP2', 'VIMP', 'C9orf142', 'LINC00996', 'NPC2', 'SLC7A5P2', 'IGFLR1', 'P2RY14', 'PHB'],
        'Myeloid Dendritic Type 1': ['WFDC21P', 'CCL17', 'CCL22', 'C15orf48', 'G0S2', 'HLA-DPB1', 'HLA-DPA1', 'HLA-DQB1', 'HLA-DQA1', 'IDO1', 'S100B', 'TXN', 'CST3', 'LGALS2', 'CCR7', 'CSF2RA', 'HLA-DQA2', 'NAPSB', 'CPVL', 'C1orf54', 'CD83', 'GSN', 'MARCKSL1', 'HLA-DRB1', 'FSCN1', 'DUSP4', 'HLA-DRB5', 'LAMP3', 'BASP1', 'HLA-DQB2', 'IL4I1', 'IRF8', 'PPA1', 'HLA-DRB6', 'GPR183', 'SERPINB9', 'DNASE1L3', 'SYNGR2', 'CLEC9A', 'IL1R2'],
        'Myeloid Dendritic Type 2': ['FCER1A', 'CD1C', 'CLEC10A', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DRB6', 'CST3', 'HLA-DRB1', 'HLA-DRB5', 'CD1E', 'HLA-DQB2', 'CPVL', 'MS4A6A', 'HLA-DMB', 'HLA-DQA2', 'NAPSB', 'LGALS2', 'HLA-DMA', 'GSN', 'FCGR2B', 'YWHAH', 'RNASE6', 'GPR183'],
        'IGSF21+ Dendritic': ['SEPP1', 'RNASE1', 'FOLR2', 'CCL2', 'CCL13', 'SLC40A1', 'C1QC', 'F13A1', 'MRC1', 'LGMN', 'STAB1', 'C1QA', 'CD14', 'MS4A6A', 'SLCO2B1', 'HMOX1', 'PDK4', 'HLA-DQA1', 'CD163', 'HLA-DPA1', 'HLA-DRB1', 'MAFB', 'GGTA1P', 'FCGR2A', 'DAB2', 'IER3', 'HLA-DQB1', 'MS4A4A', 'HLA-DPB1', 'TMEM176B', 'HLA-DRB5', 'MS4A7', 'HLA-DMA', 'C1QB', 'EGR1', 'HLA-DRB6', 'HLA-DQA2', 'A2M', 'CFD', 'CSF1R', 'CTSZ', 'KCTD12', 'GPR34', 'CTSC', 'CEBPD', 'HLA-DMB', 'NPC2', 'CST3', 'PLTP', 'GPX1', 'FCGRT', 'TGFBI', 'FCGR2B', 'AIF1', 'TMEM176A', 'CTSB', 'PLD3', 'RNASE6', 'BLVRB', 'LILRB5', 'MFSD1'],
        'EREG+ Dendritic': ['CXCL8', 'G0S2', 'IER3', 'THBS1', 'IL1B', 'TIMP1', 'CXCL2', 'PLAUR', 'CXCL3', 'EMP1', 'EREG', 'CTSB', 'RAB31', 'C15orf48', 'CTSL', 'NAMPT', 'HLA-DQB1', 'THBD', 'CD14', 'CD163', 'FCGR2A', 'CD93', 'C5AR1', 'AREG', 'CLEC10A', 'GPR183', 'PPIF', 'HLA-DQA2', 'HLA-DPB1', 'YWHAH'],
        'TREM2+ Dendritic': ['APOE', 'CCL18', 'CHIT1', 'CHI3L1', 'APOC1', 'GPNMB', 'LIPA', 'FN1', 'NUPR1', 'ACP5', 'IFI27', 'C1QB', 'FBP1', 'TREM2', 'PLD3', 'CTSD', 'CTSZ', 'C1QC', 'MMP9', 'CTSB', 'C1QA', 'HLA-DQA1', 'LILRB4', 'CYP27A1', 'CAPG', 'GRN', 'CSTB', 'TMEM176B', 'LGMN', 'PRDX1', 'LYZ', 'CD68', 'HLA-DMB', 'TXN', 'C15orf48', 'CTSK', 'IFI30', 'ALDH2', 'GSN', 'CTSH', 'ANXA2', 'CXCL9', 'HLA-DRB1', 'GPX3', 'HLA-DMA', 'LGALS3', 'HLA-DQB1', 'BRI3', 'HLA-DRB5', 'IGSF6', 'CD9', 'TYMP', 'HLA-DPA1', 'HLA-DRB6', 'HEXB', 'UBD', 'HLA-DPB1', 'CTSS', 'LTA4H'],
        'Classical Monocyte': ['S100A8', 'S100A9', 'S100A12', 'VCAN', 'LYZ', 'MNDA', 'FCN1', 'CSTA', 'THBS1', 'CTSS', 'CD14', 'AP1S2', 'TKT', 'TNFSF13B', 'FGL2', 'GCA', 'CFP', 'LGALS2', 'MS4A6A', 'CLEC4E', 'NUP214', 'CD36', 'TYMP', 'AIF1', 'MPEG1', 'FPR1', 'CYBB'],
        'OLR1+ Classical Monocyte': ['SOD2', 'TIMP1', 'C15orf48', 'IL1RN', 'IL1B', 'G0S2', 'CCL2', 'MARCKS', 'NAMPT', 'CTSL', 'CXCL8', 'AQP9', 'S100A9', 'TNIP3', 'IER3', 'BCL2A1', 'ACSL1', 'CCL4L1', 'ATP13A3', 'THBS1', 'PPIF', 'S100A8', 'HIF1A', 'CD300E', 'VCAN', 'CD163', 'CD93', 'GK', 'SERPINB9', 'CCL20', 'IFNGR2', 'SLC43A2', 'EMP1', 'SLC39A8', 'FCN1', 'NINJ1', 'UPP1', 'SLC25A37', 'OLR1', 'WTAP', 'CXCL2', 'TNIP1', 'S100A12', 'TNFRSF1B', 'RIN2', 'CLEC4E', 'LINC00152', 'APOBEC3A', 'TLR2', 'LOC731424', 'TNFAIP3', 'CEBPB', 'PLAUR', 'MXD1', 'FPR1', 'TNFAIP6', 'SPHK1', 'PFKFB3', 'RNF19B', 'ETS2', 'VEGFA', 'DSE', 'NFKB1', 'MIR4435-2HG', 'EREG', 'MMP19', 'BACH1', 'ANPEP', 'TREM1', 'GRINA', 'RAB31', 'FPR2', 'KYNU', 'TYMP', 'EHD1', 'TNFSF8', 'IVNS1ABP', 'LCP2', 'MCEMP1', 'C5AR1'],
        'Nonclassical Monocyte': ['LST1', 'LYPD2', 'FCGR3A', 'LINC01272', 'IFITM3', 'AIF1', 'COTL1', 'LILRB2', 'FCN1', 'CFP', 'FAM26F', 'CTSS', 'CFD', 'LRRC25', 'LILRA5', 'MS4A7', 'RHOB', 'C10orf54', 'MAFB', 'FCER1G', 'SERPINA1', 'IFI30', 'RHOC', 'LYN', 'SPI1', 'PILRA', 'LILRA1', 'PECAM1', 'ZFAND5', 'STXBP2', 'CSF1R', 'FGL2', 'RNASET2', 'APOBEC3A', 'HES4'],
        'Intermediate Monocyte': ['IFITM3', 'APOBEC3A', 'LST1', 'FCN1', 'LILRB2', 'COTL1', 'LINC01272', 'CD300E', 'THBS1', 'IL1B', 'CEBPD', 'NAMPT', 'LILRA5', 'VCAN', 'CFP', 'SOCS3', 'FGR', 'PRELID1'],
        'Neutrophil': ['S100A8', 'S100A9', 'IFITM2', 'FCGR3B'],
        'Eosinophil': ['SIGLEC8'],
    },
}

# 过筛数（= 源表对该型给出的可区分基因数；不随本矩阵变化）
N_SLOTS = {
    '内皮': {
        'Artery': 13,
        'Vein': 21,
        'Capillary Aerocyte': 11,
        'Capillary': 4,
        'Capillary Intermediate 1': 3,
        'Capillary Intermediate 2': 0,
        'Bronchial Vessel 1': 9,
        'Bronchial Vessel 2': 5,
        'Lymphatic': 36,
    },
    '髓系': {
        'Basophil/Mast 1': 35,
        'Basophil/Mast 2': 75,
        'Platelet/Megakaryocyte': 245,
        'Macrophage': 71,
        'Proliferating Macrophage': 87,
        'Plasmacytoid Dendritic': 63,
        'Myeloid Dendritic Type 1': 40,
        'Myeloid Dendritic Type 2': 25,
        'IGSF21+ Dendritic': 61,
        'EREG+ Dendritic': 30,
        'TREM2+ Dendritic': 59,
        'Classical Monocyte': 27,
        'OLR1+ Classical Monocyte': 80,
        'Nonclassical Monocyte': 35,
        'Intermediate Monocyte': 18,
        'Neutrophil': 4,
        'Eosinophil': 1,
    },
}

AUDIT = {
    'Artery': {'sheet': 'Cluster 16', 'n_raw_genes': 292, 'n_passed': 13},
    'Basophil/Mast 1': {'sheet': 'Cluster 44', 'n_raw_genes': 453, 'n_passed': 35},
    'Basophil/Mast 2': {'sheet': 'Cluster 45', 'n_raw_genes': 506, 'n_passed': 75},
    'Bronchial Vessel 1': {'sheet': 'Cluster 22', 'n_raw_genes': 362, 'n_passed': 9},
    'Bronchial Vessel 2': {'sheet': 'Cluster 23', 'n_raw_genes': 333, 'n_passed': 5},
    'Capillary': {'sheet': 'Cluster 19', 'n_raw_genes': 263, 'n_passed': 4},
    'Capillary Aerocyte': {'sheet': 'Cluster 18', 'n_raw_genes': 401, 'n_passed': 11},
    'Capillary Intermediate 1': {'sheet': 'Cluster 20', 'n_raw_genes': 308, 'n_passed': 3},
    'Capillary Intermediate 2': {'sheet': 'Cluster 21', 'n_raw_genes': 410, 'n_passed': 0},
    'Classical Monocyte': {'sheet': 'Cluster 55', 'n_raw_genes': 562, 'n_passed': 27},
    'EREG+ Dendritic': {'sheet': 'Cluster 53', 'n_raw_genes': 619, 'n_passed': 30},
    'Eosinophil': {'sheet': '(S1 fallback)', 'n_raw_genes': 1, 'n_passed': 1},
    'IGSF21+ Dendritic': {'sheet': 'Cluster 52', 'n_raw_genes': 571, 'n_passed': 61},
    'Intermediate Monocyte': {'sheet': 'Cluster 58', 'n_raw_genes': 524, 'n_passed': 18},
    'Lymphatic': {'sheet': 'Cluster 24', 'n_raw_genes': 820, 'n_passed': 36},
    'Macrophage': {'sheet': 'Cluster 47', 'n_raw_genes': 704, 'n_passed': 71},
    'Myeloid Dendritic Type 1': {'sheet': 'Cluster 50', 'n_raw_genes': 411, 'n_passed': 40},
    'Myeloid Dendritic Type 2': {'sheet': 'Cluster 51', 'n_raw_genes': 279, 'n_passed': 25},
    'Neutrophil': {'sheet': '(S1 fallback)', 'n_raw_genes': 4, 'n_passed': 4},
    'Nonclassical Monocyte': {'sheet': 'Cluster 57', 'n_raw_genes': 549, 'n_passed': 35},
    'OLR1+ Classical Monocyte': {'sheet': 'Cluster 56', 'n_raw_genes': 740, 'n_passed': 80},
    'Plasmacytoid Dendritic': {'sheet': 'Cluster 49', 'n_raw_genes': 660, 'n_passed': 63},
    'Platelet/Megakaryocyte': {'sheet': 'Cluster 46', 'n_raw_genes': 1134, 'n_passed': 245},
    'Proliferating Macrophage': {'sheet': 'Cluster 48', 'n_raw_genes': 969, 'n_passed': 87},
    'TREM2+ Dendritic': {'sheet': 'Cluster 54', 'n_raw_genes': 573, 'n_passed': 59},
    'Vein': {'sheet': 'Cluster 17', 'n_raw_genes': 276, 'n_passed': 21},
}


def build_panel(lineage, available_genes=None):
    """返回某一谱系的面板 {亚型: [基因...]}。

    available_genes : 本数据矩阵的基因集合。给了则：
        ① 每型按 logFC 降序**顺延**取前 PANEL_CAP 个真在矩阵里的基因
           （矩阵里没有的跳过、继续往下取，使各型面板等长）；
        ② 返回 (面板, 缺失清单) —— 缺失清单记的是**原名单里所有**不在矩阵的基因。
    不给则原样返回（不封顶）。
    """
    if lineage not in PANELS:
        raise KeyError(f"未登记的谱系：{lineage}；有 {list(PANELS)}")
    if available_genes is None:
        return {k: list(v) for k, v in PANELS[lineage].items()}
    have = set(available_genes)
    used, missing = {}, {}
    for t, gs in PANELS[lineage].items():
        miss = [g for g in gs if g not in have]
        used[t] = [g for g in gs if g in have][:PANEL_CAP]
        if miss:
            missing[t] = miss
    return used, missing


def assert_provenance():
    """硬校验（与 classic_panels.assert_provenance 同精神）。"""
    problems = []
    for lineage, d in PANELS.items():
        if not d:
            problems.append(f"{lineage}: 面板为空")
        for t, gs in d.items():
            if len(gs) == 0 and t not in UNSCORABLE_FROM_S4 and t not in FALLBACK_FROM_S1.get(lineage, {}):
                problems.append(f"{lineage}/{t}: 名单为空且未登记")
            if len(gs) != len(set(gs)):
                problems.append(f"{lineage}/{t}: 名单内有重复基因")
            bad = [g for g in gs if not g or g.lower() == "nan" or " " in g]
            if bad:
                problems.append(f"{lineage}/{t}: 非法基因名 {bad}")
    for lineage, d in FALLBACK_FROM_S1.items():
        for t in d:
            if t not in PANELS[lineage]:
                problems.append(f"fallback {lineage}/{t} 未进入面板")
    if problems:
        raise AssertionError("面板溯源校验失败：\n  " + "\n  ".join(problems))
    return True

