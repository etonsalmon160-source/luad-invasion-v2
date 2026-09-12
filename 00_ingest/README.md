# 00_ingest —— 权威输入层（第一步）

> 目标：把"用哪些数据集、每个样本属哪个阶段、什么模态"**按权威口径冻结**，供下游一致引用。

## 铁律
1. 数据身份以 **GEO/GSA 为准**；分期**无静默默认**。
2. `patient_id`（真患者）与 `sample_id`（组织/切片）**分层**。
3. 数据未到位 → **如实报缺**，不造数据。

## 单细胞队列（scRNA，**均无 AAH**）
| 数据集 | 模态 | 阶段 |
| :--- | :--- | :--- |
| GSE131907 | scRNA | 按 `Sample_Origin`：`nLung=Normal`、`tLung/tL/B=IAC`、`mLN=LNM`、**`nLN=正常淋巴结(非LNM)`**、`mBrain=脑转移`、`PE=胸腔积液` |
| GSE189357 | scRNA | `TD1/2/9=IAC`、`TD3/4/6=MIA`、`TD5/7/8=AIS`（**无 AAH**） |
| GSE148071 | scRNA | **Advanced NSCLC** |

## AAH 单细胞（决策）
- **暂用**：`GSE308103`（snRNA）→ 经 **SCMG zero-shot 跨平台并入**；
- **预留接口**：`HRA001130`（全细胞 scRNA，GSA-Human **受控，待申请**）→ 见 [`hra001130_interface.py`](hra001130_interface.py)；
- **跨模态五判据**（缺一不可，见白皮书 §3.1）：① 各阶段**可分辨** ② 重叠阶段 sn↔sc **一致** ③ 平台偏移 **δ(stage) 稳定** ④ **sn 内部配对**同向 ⑤ AAH 身份 **CNV/标记**证真。

## 文件
- [`cohort_registry.py`](cohort_registry.py) —— 权威队列登记表（single source of truth）
- [`hra001130_interface.py`](hra001130_interface.py) —— HRA001130 预留接口

## 待办（下一步）
- [x] `fetch_geo_metadata.py`：拉取并冻结 GEO 逐样本权威元数据（`geo_metadata/`）—— 提供 `patient id` / histology / origin。
- [x] `01_load_cohorts.py`：按 registry 纳入三 scRNA 队列（GEO 真值分期 + 分层），
      输出冻结的 per-cell 表 + 校验报告（细胞数/分期分布/哈希）。→ **M0 已过门**（见 `PLAN_AND_CHECKPOINTS.md` 过门证据）。
- [ ] `GSE308103` 按 snRNA 口径纳入（含 QC，组件见 tools/）—— 属 **M4**。

## 患者身份口径（重要）
- **GSE131907**：`patient_id` 取自 GEO `patient id`（权威；44 患者，形如 `P0001/P1006/P2001/P3002`）。
  **禁止**按样本名尾号推断——`LUNG_N06`(P0006) / `EBUS_06`(P1006) / `LN_06`(P2006) / `NS_06`(P3006) 是**四个不同患者**（同尾号陷阱）。
- **GSE189357 / GSE148071**：GEO 无 `patient id` 字段，按 **1 样本 = 1 患者** 显式假设（`frozen_samples.csv` 有 `patient_rule`/`patient_confidence` 标注）。

## 文件（本步产出，`results/00_ingest/`）
| 文件 | 内容 |
| :--- | :--- |
| `frozen_per_cell.csv.gz` | per-cell 表（cell_barcode / patient_id / sample_id / stage / dataset / raw_barcode） |
| `frozen_samples.csv` | 样本表（+ modality / patient_rule / patient_confidence） |
| `frozen_patients.csv` | 患者表（样本数 / 细胞数 / 涉及阶段） |
| `frozen_source_files.csv` | 源文件清单（路径 / 大小 / mtime / SHA-256） |
| `frozen_manifest.json` | 冻结清单（各产物哈希 + 计数 + 过门结果） |
| `M0_validation_report.md` | 人读校验报告（逐条过门 + OPEN ISSUE） |
