# 诊断性预注册 · 表达型恶性判别器在我们数据上能不能分开

**登记时间**：2026-09-23（**计算开始之前**，见文件 mtime）
**性质**：🔬 **诊断（diagnostic）**，不是交付物。**不改变任何已签口径**，不产出可引用的分期结论。
**授权**：用户 2026-09-23 「去做」——授权执行本对照实验（此前我提出的"第 1 条零承诺对照实验"）。

---

## 一、要回答的问题

CopyKAT 的锚定判决口径已废（见 `results/03_cnv/ANCHOR_PREREG.md` §7 机械裁决 = 不可用）。
已发表的前辈工具里有若干**不用 CNV、不用正常参考**的表达型恶性判别器。
**本诊断只问一件事：它们的打分在我们这批单核数据上，能不能把「正常上皮」和「瘤上皮」分开。**

能分开 ⇒ 存在一条与 CNV 独立的第二轴（仍是候选，不是结论）。
分不开 ⇒ 佐证"表达分不开"，CNV 路线更有必要。

## 二、工具与资源（冻结）

| 项 | 值 |
|---|---|
| 工具 | `scMalignantFinder`（PyPI，conda 环境 `scmalignant`，Python 3.10.10） |
| 预训练模型 | `tools/scMalignantFinder_pretrain/model.joblib`（**22,511 B**） |
| 特征表 | `tools/scMalignantFinder_pretrain/ordered_feature.tsv`（**2707** 个基因，与论文 2,707 DEG 吻合） |
| 来源 | Zenodo record 17888140 |
| 论文 | Communications Biology 2025;8:504, DOI 10.1038/s42003-025-07942-y, PMID 40148533 |

**该模型不用 CNV**（纯 2707 基因的逻辑回归），**不需要用户提供正常参考细胞**。

## 三、输入（冻结）

- 表达矩阵：`results/02_expression/gse308103_counts_paperqc.h5ad`（sha256 `a276cd1af69a4620ff9e3d64f4b93f33461e646537c5f2c9da1aff8c9150b9de`，413,697 细胞 × 18,069 基因）
- 细胞子集：`results/05_annotation/epiA_subset_barcodes.txt`（**133,384** 个上皮细胞，sha256 `d12a115a87d0790d7fe8cf27b16170999d94daba51f3016fdb79c5aa3fdde55e`）
  —— 该清单为**未再做剔除**的原始上皮子集（**不是** nocontam 集），避免引入额外的事后筛选决定。
- 分期标签：`results/04_integration/seurat_trad/epiA/clusters.csv.gz` 的 `stage` 列（来自 `00_ingest/cohort_registry.resolve_stage`，R1 权威口径）

## 四、评估设计（**先写死，不得事后改**）

**分组只用样本级分期（弱标签），不用任何逐细胞真值**（我们没有逐细胞真值——这正是问题本身）：

- 正类（"瘤"）：`stage == IAC` 的上皮细胞
- 负类（"正常"）：`stage == Normal` 的上皮细胞
- 描述组（**只看分布，不参与打分评判**）：AAH / AIS / MIA

**主指标：逐患者 AUROC。** 对每个患者，用他本人的 Normal 上皮 vs IAC 上皮算一个 AUROC，
再取患者间中位数。**必须逐患者算**——否则一个病人的大批细胞会淹没全体。

**预注册判读（看到结果前定死）**：

| 逐患者 AUROC 中位数 | 判读 |
|---|---|
| ≥ 0.80 | 表达轴**有信号**，可作候选第二轴（仍须另行签字才能采用） |
| 0.60 – 0.80 | 弱信号，只能作佐证，不足以定恶性 |
| < 0.60 | 表达**分不开**，不采用 |

**不得**在看到结果后调阈值、换基因表、或改分组定义（法则 3.2）。

## 五、必须一并上报的已知混淆（**不许省略**）

1. **样本级批次混淆**：Normal 样本与 IAC 样本在解离、环境 RNA、测序批次上都不同。
   判别器若靠这些技术差异分开，**不是恶性信号**。本设计**无法完全排除**这一混淆——只能如实标注。
2. **瘤旁正常上皮**：IAC 样本内部也有正常上皮，会被本设计误标为正类，**拉低**表观性能。
3. **无逐细胞真值**：所以本诊断**只能**回答"能不能分开两个样本群"，**不能**回答"每个核是不是恶性"。
4. **单核未验证**：该工具及所有同类前辈**均未在 snRNA 上验证过**；本诊断是该工具在单核上的**首次**使用。
5. **无 AAH/AIS 前例**：该工具论文的前驱病变验证只到结直肠息肉与胃 NAG→CAG→IM→EGC，**从未碰过 AAH/AIS**。

## 六、产物

- `malig_probe_per_cell.csv.gz`：逐细胞打分（barcode, sample_id, patient_id, stage, epi_subtype, score, pred）
- `malig_probe_per_patient.csv`：逐患者 AUROC
- `malig_probe_summary.json`：判定结果与全部口径常量
- `figures/malig_probe_score_by_stage.png`：各分期打分分布（**图内标签全英文**）
