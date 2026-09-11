# LUAD 多组学课题 · 最新总纲 (PROJECT_SUMMARY)

> **本文件定位**：全项目**唯一的最新事实总纲（authoritative facts ledger）**。
> 旧的报告式总结（`MASTER_SCIENTIFIC_REPORT.md` 等）已因内容与查实事实冲突而隔离至
> [`archive/legacy_docs/`](archive/legacy_docs/)。
> **叙事型白皮书**仍是 [`LUAD_SPATIOTEMPORAL_INVASION_ARCHITECTURE_AND_PROGRESS.md`](LUAD_SPATIOTEMPORAL_INVASION_ARCHITECTURE_AND_PROGRESS.md)（本总纲不替代它，仅记录已核验事实与整改状态）。
> **更新日期**：2026-09-12

---

## 0. 阅读顺序
1. 白皮书 = 研究设计与叙事；2. 本总纲 = 已核验事实 / 真值口径 / 整改状态；3. `PIPELINE_MANIFEST.md` = 脚本清单；4. `FINAL_RIGOROUS_AUDIT_VERDICT.md` = 前任终审裁决（历史）。

---

## 1. 数据底座（GEO 权威核实）

| 数据集 | 真实身份（GEO） | 用途 | 核实结论 |
| :--- | :--- | :--- | :--- |
| **GSE131907** | "Single cell RNA sequencing of lung adenocarcinoma"（Kim et al., 208,506 细胞 / 58 样本 / 44 患者；含 原发/正常肺/正常淋巴结/转移淋巴结/脑转移/胸腔积液） | 单细胞基座（Normal/IAC/LNM） | ✅ 身份正确 |
| **GSE189357** | "Spatiotemporal transcriptional atlas of LUAD from AIS to IAC [scRNA-seq]"，**9 样本 = 3 AIS + 3 MIA + 3 IAC，无 AAH** | 单细胞基座（早期） | ⚠️ 项目旧标注有误（见 §3） |
| **GSE148071** | "Comprehensive Profiling of Cancer Cells and Their Microenvironment in **Advanced NSCLC**"，42 样本 | 独立验证队列 | ⚠️ 实为**晚期 NSCLC**，非"早期 LUAD" |
| **GSE307534** | 空间转录组（Visium CytAssist FFPE，19+ 患者；每患者含 AAH/AIS/MIA/LUAD 切片） | 空间主队列 | ✅ 切片身份正确 |
| **GSE308103** | snRNA-seq（FFPE 卷片，与 GSE307534 同研究；含 normal/AAH/AIS/MIA/LUAD，75 样本） | **真 AAH 单细胞来源** | ✅ 已下载并 QC（见 §3.4） |

---

## 2. 双分支架构（关键）
- **传统分支（阶段 1–5）**：scanpy/Seurat、RCTD、PLIP、SpaGCN/Squidpy、CellRank/PAGA/DPT 等**传统算法**。
- **SCMG 分支（阶段 6，正交）**：仅用 **SCMG 官方神经网络**（`CellEmbedder` 编码 → 条件扩散模型 `generate_transition_cells` → `CausalGenePredictor` 因果基因），**不得掺入任何传统算法**，作为对传统结果的**独立正交对照**。
- 两分支最后做**互证/分歧**报告。

---

## 3. 单细胞基座：已查实的问题与权威口径

### 3.1 阶段映射（权威 = GEO 逐样本 histology）
- **GSE189357 真值**：`TD1/2/9 = IAC`、`TD3/4/6 = MIA`、`TD5/7/8 = AIS`；**无 AAH**。
  - 旧错误：`TD9` 被标 **AAH**、`TD4` 被标 IAC → 项目里 **14,312 个"AAH"细胞实为 IAC**。
- **GSE131907 真值**（按 `Sample_Origin`）：`nLung→Normal`、`tLung/tL/B→IAC`、`mLN→LNM`、`nLN→正常淋巴结`、`mBrain→脑转移`、`PE→胸腔积液`。
  - 旧错误：`nLN`(35,169) 被记为 **LNM**；`mBrain`(21,158)+`PE`(19,755) 被并入 **IAC**；`LNM` 一档约 60% 是正常淋巴结。
- **GSE148071**：应如实标注为 **Advanced NSCLC**（不做"早期 LUAD"表述）。
- **原则**：一切静默默认（如 `.get(x,'IAC')`）必须取消。

### 3.2 细胞类型映射
- 旧 `Malignant` 由 `62_super_continent.R` 用**泛上皮标记**（EPCAM/KRT8/18/19/CDH1/SFN）argmax 判定 → **不是恶性特异**，`Pure_Tumor` 可能混入正常上皮。
- **必须用 CNV 推断**（inferCNV/CopyKAT）证真恶性 vs 正常上皮。
- 参照流形无肿瘤态 → SCMG 官方映射/隐空间**均无法区分恶性 vs 正常上皮**（已实测：余弦轮廓系数 0.035、kNN 纯度 0.68）。

### 3.3 单细胞管线参数（沿用性）
- **可沿用**：QC（nFeature 500–10000 / nCount 1000–60000 / mt<10%）、LogNormalize、HVG(vst,2000)、PCA(50)、Harmony、聚类。
- **须重做**：① **双体去除**（`36` 用 nCount+nFeature z-sum 取高 2.5%，非 DoubletFinder，会误删好细胞）；② **恶性标签**（CNV）；③ **阶段映射**（GEO 真值）。
- **文档/代码不一致**：Harmony `dims` 代码为 **1:50**，白皮书写 1:30。

### 3.4 真 AAH（GSE308103，已下载 + QC）
- 9 个 AAH 样本、~62,892 细胞；**上皮占比 30–51%（最大谱系）**，纤维 19–31% → **非纤维富集，可用**。
- 口径差异：**snRNA（核）**，纳入时须注明。
- QC 结果：`results/tables/gse308103_aah_qc.csv`。

---

## 4. 空间（Visium）队列：真值
- **核心 6 阶段**：Normal `P4_Normal`、AAH `P1_AAH`、AIS `P3_AIS`、MIA **`P10_MIA`**(GSM9226189, 真 MIA)、IAC `P3_LUAD`、LNM **`PT_3_LNM`**(GSE190811)。
- **P0 级禁令**（[`.agents/rules/spatial_cohort_and_figure_prohibitions.md`](.agents/rules/spatial_cohort_and_figure_prohibitions.md)）：禁止 `P4_AAH1`/`P4_AAH` 充当 MIA；禁止把原发切片伪标 LNM。
- GEO 已证：`P1_AAH`(GSM9226168) = 真 AAH；`P10_MIA` = 真 MIA；**空转 AAH 无顶替**。
- ⚠️ `17_render_6stage_master_multimodal_matrix_complete.py` 仍用 `P4_AAH1` 当 MIA（旧版，须处理）。

---

## 5. 已查实的"伪数据/伪算法"与处置
| 对象 | 问题 | 处置 |
| :--- | :--- | :--- |
| `pipeline/11_tmb_...py` | 写死 TMB 表、HR、`np.exp` 伪 KM，无 Cox | 待重算/隔离 |
| `pipeline/107_lincs_...py` | CMap 结果反推伪造，参考数据是 28×8 stub | ✅ 已改诚实版（无真数据即拒绝）；伪表已隔离 |
| `117_scmg_master_316k_pipeline.py` | 伪时间/逆转概率写死（`stage_order` 字典 + `1-pseudotime`） | ✅ 已隔离 |
| `merge_316k_umap.py` | 传统 UMAP + 写死分期表造逆转概率 | 已识别 |
| `29/31/64_*.py` | 违禁样本映射（P4_AAH1→MIA、P3_LUAD→LNM） | ✅ 已隔离 |
| `13/44/45/53/60/40/41/85/110` 等 | 写死数值/合成曲线/伪造矩阵 | 待清理（见 `PIPELINE_MANIFEST` 与审计） |
| SCMG 全量逆转（118） | 真算，但全局相关 → 髓系污染（MCEMP1 等） | 传统口径问题，已记录 |

---

## 6. 单细胞脚本"最后一版"链（其余为旧版干扰）
`30_scRNA_qc_filtering.R` → `34_super_atlas_harmony.R` → `36_doublet_removal.R` → `62_super_continent.R` → `subcluster_tumor_paper.R` → `build_14tumor_tme_reference.R` → `01_extract_subatlas_annotations.R` → `54_rigorous_paper_harmony_subclustering.R` / `calibrate_single_cell_gating.R` / `reannotate_immune_canonical.R` / `sync_master_metadata.py` → `101_rebuild_accurate_metadata_cache.py`

产出：`results/metadata_cache/all_subatlases_barcode_to_label_mapping.csv`（316,689 细胞；**当前版本阶段/标签有误，待重建**）。

---

## 7. 整改路线（最高规格重建）
1. **阶段映射权威化**（GSE189357 TD 真值 / GSE131907 `Sample_Origin` / GSE148071 如实），取消静默默认。
2. **补真 AAH**：整合 GSE308103（snRNA，注明口径）。
3. **CNV 证真恶性标签**（inferCNV/CopyKAT）。
4. **重做双体去除**（DoubletFinder/scDblFinder）。
5. **隔离旧版脚本**（每功能只留最后一版）。
6. **下游连锁重做**：RCTD 解卷积（依赖单细胞签名）、空间 6 阶段、SCMG 分支（纯 SCMG）。
7. **重出受影响图版**。

> **原则**：不得再出现"标签/分期/数值未经权威核实即入库"；一切结果须有可追溯的真实来源。

---

## 8. 文档目录
- 白皮书（叙事，权威）：[`LUAD_SPATIOTEMPORAL_INVASION_ARCHITECTURE_AND_PROGRESS.md`](LUAD_SPATIOTEMPORAL_INVASION_ARCHITECTURE_AND_PROGRESS.md)
- 本总纲（事实）：`PROJECT_SUMMARY.md`
- 旧报告归档：[`archive/legacy_docs/`](archive/legacy_docs/)
- 脚本清单：`PIPELINE_MANIFEST.md`；历史终审：`FINAL_RIGOROUS_AUDIT_VERDICT.md`
- 规则：`.agents/rules/`（含 P0 空间样本禁令）
