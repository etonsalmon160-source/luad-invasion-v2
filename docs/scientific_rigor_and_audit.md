# 科研严谨性与审计准则 (Scientific Rigor & Audit Protocol)

> 适用于本项目（配对空间-单核图谱）。**所有脚本、报告、图表必须遵守。**
> 范围与里程碑见 [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md) 与 [`WHITEPAPER.md`](WHITEPAPER.md)。

---

## 法则 0 · 数据身份与分期真值
1. 数据集身份以 **GEO/GSA** 为准；本项目**仅两个配对数据集**：`GSE308103`(snRNA) + `GSE307534`(空间)；
2. **分期无静默默认**（严禁 `.get(x,'IAC')` 之类回退）；未知 token 一律 raise；
3. **恶性标签须 CNV 证真**（CopyKAT），**不得**以泛上皮标记（EPCAM/KRT…）argmax 代替；
4. 双体必须用**标准算法**（scDblFinder），**不得**以 nCount+nFeature 启发式替代。
5. **病灶序号无推导**：同一患者的第二病灶**不得由 token 猜** —— 两个数据集的命名约定**不同**
   （`GSE307534` 用 `AAH-1`/`AIS-1`（**带横线**）；`GSE308103` 用 `AAH1`/`AIS1`/`Normal1`/`LUAD1`（**无横线**））。
   序号**只查 GEO 权威表**（`cohort_registry.resolve_lesion_ordinal`），查不到即 raise。

---

## 法则 1 · 伪造与模拟数据零容忍
1. **严禁** `np.random` / `runif` / `rnorm` / `sample()` / `make_blobs()` / 硬编码虚拟数值来伪造细胞比例、表达或空间丰度；
2. 每个图表/统计的输入必须来自**真实测序矩阵**，并具明确来源；
3. 发现伪造脚本 → **立即隔离**，禁止留在生产主干。

---

## 法则 2 · 生物学特异性（marker 交叉校验）

> 🔴 **2026-09-17 重建（用户指令「marker 要来自权威文章，不要给我从犄角旮旯搞来」）**。
> 本表旧版（EPCAM/KRT7/19… 六行）**一个出处都没有**，属审计标缺。现按**一次文献**重建，
> 逐基因出处见 [`05_annotation/marker_panel.py`](../05_annotation/marker_panel.py)。
> **执行口径 = 该文件，不是本表**；本表为人读摘要。

| 谱系 | 阳性 marker | 应为阴性 | 出处（一次文献） |
| :--- | :--- | :--- | :--- |
| 上皮 | EPCAM, KRT8/18/19, CDH1, NKX2-1, SFTPC, SFTPA1/A2, SFTPB, NAPSA, AGER, CAV1, PDPN, SCGB3A2, SCGB1A1 | PTPRC, CD3D, PECAM1, COL1A1 | Travaglini 2020 *Nature* 587:619（AT1/AT2/Club/pan-上皮）；Vieira Braga 2019 *Nat Med* 25:1153（气道）；**彭 2026 *Cancer Cell*（源论文，LUAD 谱系 NKX2-1）** |
| T/NK | CD3D/E/G, TRAC, CD4, IL7R, CD8A/B, NKG7, GNLY, KLRD1, PRF1, GZMB | EPCAM, COL1A1, CD68 | Travaglini 2020；Vieira Braga 2019；Guo 2018 *Nat Med* 24:978（NSCLC T 细胞） |
| B/浆细胞 | MS4A1, CD19, CD79A/B, MZB1, JCHAIN, SDC1, IGHG1, IGKC, XBP1, DERL3 | CD3D, EPCAM, ACTA2 | Travaglini 2020；Vieira Braga 2019 |
| 髓系 | LYZ, AIF1, ITGAX, CD68, CD163, MSR1, C1QA/B/C, MARCO, APOE, FCN1, CD14, S100A8/9, SPP1 | CD3D, EPCAM, PECAM1 | Travaglini 2020（肺泡巨噬）；Habermann 2020 *Sci Adv* 6:eaba1972（SPP1⁺ 巨噬）；Zilionis 2019 *Immunity* 50:1317（肺肿瘤髓系） |
| 成纤维 | COL1A1/A2, COL3A1, DCN, LUM, FN1, PDGFRA/B, ACTA2, TAGLN, FAP, CXCL12 | PTPRC, EPCAM, PECAM1 | Travaglini 2020；Habermann 2020；Reyfman 2019 *AJRCCM* 199:1517；Lambrechts 2018 *Nat Med* 24:1277（肿瘤基质） |
| 内皮 | PECAM1, CDH5, KDR, CD34, VWF, EGFL7, RAMP2, EMCN, PLVAP, AQP1, CLDN5, FLT1 | EPCAM, PTPRC, COL1A1 | Travaglini 2020（EC 亚型）；Gillich 2020 *Nature* 586:785（肺泡毛细血管 aCap/gCap 特化）；Lambrechts 2018 |

**恶性身份以 CNV 为准**；marker 仅作**一致性佐证**，不作判据。

**与 GP6 标准 B 的独立性（重要）**：本表是 GP6 **标准 A** 的基因来源，**刻意不取自 HLCA 整合图谱**
（Sikkema 2023 *Nat Med* 29:1563）——因为 CellTypist 的 `Human_Lung_Atlas.pkl` 正是该图谱训出来的，
若两者同源，则 κ 变成自证。**残留非独立性如实声明**：HLCA 整合了上表所引的多套一次研究数据，
故两者**并非统计独立**，只保证「marker 定义来源不同」。

---

## 法则 3 · 参数透明与统计严谨
1. 参数**显式常量化**并登记 [`PARAMETERS_AND_SOURCES.md`](PARAMETERS_AND_SOURCES.md)（标出处与核对状态）；
2. **阈值一经确定不得事后调整**（防 p-hacking）；
3. 单细胞差异用**患者级 pseudobulk**（DESeq2），不用单细胞级 Wilcoxon（伪重复）；
4. 组成性比例数据报**效应量**，不用裸 FDR；
5. 空间统计用**类别标签**、报 **z + 经验 p**（`squidpy.nhood_enrichment` **不返回 p 值**；`co_occurrence` **不做置换**）；
6. 生存分析须输出 HR、95% CI、Log-rank P，**严禁篡改显著性**。

---

## 法则 4 · 因果与措辞边界
1. **观察性单细胞/空间不能建立因果**；合法杠杆 = **cis-MR + coloc** 或扰动实验；
2. 可说：**candidate / genetically supported candidate / prognostic association / computational hypothesis**；
3. **不可说**：`causal`、`driver`（无限定）、`validated target`、"状态逆转因子"；
4. 对接（M8）结论仅为**计算假说**；"validated" 须湿实验。

---

## 法则 5 · 跨模态（本项目核心）
1. 空间（Visium spot，多细胞混合）**必须解卷积**；参考须**模态匹配**（FFPE↔FFPE → 用 GSE308103）；
2. 跨模态一致性须过 **M4 五判据**；**不过门 → AAH 只能标为假说**；
3. 禁用 `modality 当 batch`（sc↔sn 是 system 效应且与数据集共线）。

---

## 法则 6 · 冻结与可复现
1. 产物**可复现 + 有哈希**；确定性（种子显式、无随机）；
2. `patient_id`（真患者）与 `sample_id`（切片/样本）**分层**；
3. 输入输出可溯源；冻结清单与校验报告随里程碑提交。

---

## 审计流程（产物入库前）
1. **数据源头核验**：路径有效、来源明确；
2. **计数/标签一致性核验**：与冻结表对齐；
3. **可复现性核验**：确定性脚本 + 种子显式 + 产物哈希；
4. **措辞与方法学比对**：对照本准则与顶刊标准；
5. **守卫须被证伪**：凡新增或修改的一致性校验，必须喂一份**故意做坏的输入**，确认它真的会红
   （退出码非 0、且点名具体条目）。**只"跑通"不算验证 —— 永远通过的守卫是假的。**
   实例：`00_ingest/03_verify_cohort_consistency.py` 用被删规则的坏表自证，点名 5 个样本并退出 1。
6. **同一事实不得有两份实现**：凡同一事实存在第二份推导（**尤其是"从未被调用"的那一份**），
   要么删除，要么与权威源逐字对撞。**无执行路径的缺陷不产生症状、不会自我暴露**，
   单测/流水线/产物核验全都覆盖不到它（实例：`lesion_ordinal()`，
   见 [`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §6.1）。
