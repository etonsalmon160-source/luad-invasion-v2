# M0 输入冻结 · 校验报告

- schema: `M0-frozen-1`  · 生成时间(UTC): 2026-09-12T06:13:07Z
- 数据源(只读): `/home/eto/luad_invasion/data`
- registry: `00_ingest/cohort_registry.py` sha256=`7ad44013fd013dedeb80e77be23838b03d9a22cfab1c6d7860743ebf8d37b4cb`

## 过门判定

**✅ PASS**

| # | 检查项 | 结果 | 详情 |
| :-- | :--- | :---: | :--- |
| | C1 per-cell 字段齐全 | ✅ | 存在 ['cell_barcode', 'patient_id', 'sample_id', 'stage', 'dataset'] |
| | C2 关键字段无空值 | ✅ | 空值统计={'cell_barcode': 0, 'patient_id': 0, 'sample_id': 0, 'stage': 0, 'dataset': 0} |
| | C3 cell_barcode 全局唯一 | ✅ | 重复=0 |
| | C4 阶段 ∈ 合法词表 | ✅ | 实际阶段=['AIS', 'Adv_NSCLC', 'Brain_Met', 'IAC', 'LNM', 'MIA', 'Normal', 'Normal_LN', 'Pleural_Effusion'] |
| | C5 GSE189357 分期 == GEO 真值 (TD1/2/9=IAC, TD3/4/6=MIA, TD5/7/8=AIS) | ✅ | 实际={'IAC': ['TD1', 'TD2', 'TD9'], 'MIA': ['TD3', 'TD4', 'TD6'], 'AIS': ['TD5', 'TD7', 'TD8']} |
| | C6 GSE189357 AAH 细胞数 == 0 | ✅ | AAH=0 |
| | C7 GSE131907 nLN 样本仅映射 Normal_LN（非 LNM） | ✅ | nLN 样本=['LN_01', 'LN_02', 'LN_03', 'LN_04', 'LN_05', 'LN_06', 'LN_07', 'LN_08', 'LN_11', 'LN_12'], 其阶段=['Normal_LN'] |
| | C8 GSE131907 脑转移/胸水单列（未并入 IAC） | ✅ | 实际阶段集合=['Brain_Met', 'IAC', 'LNM', 'Normal', 'Normal_LN', 'Pleural_Effusion'] |
| | C9 GSE148071 == Adv_NSCLC | ✅ | 实际=['Adv_NSCLC'] |
| | C10 GSE189357 样本集 == registry TD 集 | ✅ | 观测=['TD1', 'TD2', 'TD3', 'TD4', 'TD5', 'TD6', 'TD7', 'TD8', 'TD9'] registry=['TD1', 'TD2', 'TD3', 'TD4', 'TD5', 'TD6', 'TD7', 'TD8', 'TD9'] |
| | C11 GSE131907 患者数 == GEO 真值 44 | ✅ | distinct=44, 前缀正确=True |
| | C12 GSE131907 GEO origin 与 registry 分期逐样本一致 | ✅ | 58/58 一致 |
| | C13 GSE189357 GEO histology == registry 分期（独立核验） | ✅ | 9/9 一致，AAH 缺席 |
| | C14 GSE131907 患者数 == 44 | ✅ | 实际=44 |
| | C14 GSE189357 患者数 == 9 | ✅ | 实际=9 |
| | C14 GSE148071 患者数 == 42 | ✅ | 实际=42 |
| | 自审：禁用模式（np.random / random.* / .get( ）| ✅ | 无 |

## 计数

- 总细胞数：**420,766**
- 样本数：109  患者数：95

### 按数据集

| 数据集 | 细胞数 | 样本数 | 患者数 | 模态 |
| :--- | ---: | ---: | ---: | :--- |
| GSE131907 | 208,506 | 58 | 44 | scRNA |
| GSE148071 | 89,887 | 42 | 42 | scRNA |
| GSE189357 | 122,373 | 9 | 9 | scRNA |

### 按阶段

| 阶段 | 细胞数 | 样本数 |
| :--- | ---: | ---: |
| AIS | 38,107 | 3 |
| Adv_NSCLC | 89,887 | 42 |
| Brain_Met | 29,060 | 10 |
| IAC | 108,075 | 18 |
| LNM | 21,479 | 7 |
| MIA | 33,413 | 3 |
| Normal | 42,995 | 11 |
| Normal_LN | 37,446 | 10 |
| Pleural_Effusion | 20,304 | 5 |

### GSE189357 逐样本（GEO 真值核对）

| sample_id | stage | 细胞数 |
| :--- | :--- | ---: |
| TD1 | IAC | 15,216 |
| TD2 | IAC | 18,064 |
| TD3 | MIA | 12,658 |
| TD4 | MIA | 11,756 |
| TD5 | AIS | 19,079 |
| TD6 | MIA | 8,999 |
| TD7 | AIS | 3,507 |
| TD8 | AIS | 15,521 |
| TD9 | IAC | 17,573 |

## 产物哈希（SHA-256）

- `frozen_per_cell.csv.gz` ： `703c5f03562d031a6aea96a68dea28f312733330db76fa0f0ed361a61657f7d2`
- `frozen_samples.csv` ： `4696a88b9c69a85d93eae742547cc3537e4328e1b7cf0c0bbe2cc690e10d980a`
- `frozen_patients.csv` ： `6945e2b71f72fcd19547e30edea118fd241bf7e8480c5f9d1b2fb9fc21a2f13a`
- `frozen_source_files.csv` ： `ed834ed6618289772573a38f22b01f11f3c59f68571ed48e2ac603955ae6a284`

## 患者身份（GEO 权威）

- **GSE131907**：GEO_VERIFIED (44 patients) — 规则：GEO: !Sample_characteristics_ch1 'patient id'（权威字段）；患者数=44
- **GSE148071**：ASSUMED (no patient field in GEO) — 规则：assumed: 1 样本 = 1 患者（GEO 无 patient id 字段）；患者数=42
- **GSE189357**：ASSUMED (histology GEO-verified) — 规则：assumed: 1 样本 = 1 患者（GEO 无 patient id；histology 由 GEO 核实）；患者数=9

## 已知显式假设（非静默默认）

- GSE131907：`patient_id` 取自 **GEO `patient id`** 字段（权威，44 患者，形如 P0001/P1006/P2001/P3002）；**不用样本名尾号推断**（同尾号跨系列并非同一患者）。
- GSE189357 / GSE148071：GEO 无 `patient id` 字段，按 **1 样本 = 1 患者** 显式假设登记（见 `frozen_samples.csv` 的 `patient_rule`/`patient_confidence`）；取得患者表后回填重跑。GSE189357 的 **histology 已由 GEO 独立核实**。
- 源文件中 `.rds`（>2 GiB）默认不做 SHA-256（表中标注跳过）；如需全量哈希，加 `--hash-large`。

## GEO 交叉核验

- GSE131907 `tissue origin abbrevation` ↔ registry `stage_by_origin`：逐样本一致（C12）。
- GSE189357 `histolgical type` ↔ registry `stage_by_TD`：逐样本一致（C13），**独立确认无 AAH**。

## 范围说明

- 本门仅冻结**三个 scRNA 队列**；AAH 单细胞（GSE308103 sn / HRA001130）属 **M4 跨模态 AAH**，不在 M0 范围。
- 未做任何 QC / 双体 / 恶性标注（属 M1/M2）。
