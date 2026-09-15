# M0 输入冻结（配对数据集）· 校验报告

- schema `M0-paired-1` · 2026-09-15T14:02:40Z · raw `/home/eto/luad_invasion/data`
- registry sha256 `228617ee52d93d1d4a5c184ccd43e87fcd52d81724070dbce974d5a47d9a60d2`
- GEO 权威表 sha256 `{'GSE308103': 'c25de6ae144e7ff3947b903499c5afb2cfaf51266ba1664edbded90b2dae2795', 'GSE307534': '65b935cd5c0332dea263c8bc842524a53cdef9ec580a253db9e53c1af2a84893'}`

## 过门判定：**❌ FAIL**

| 检查 | 结果 | 详情 |
| :-- | :--: | :-- |
| C1 两数据集均在册 | ✅ | ['GSE307534', 'GSE308103'] |
| C2 分期 ∈ 词表 | ✅ | ['AAH', 'AIS', 'IAC', 'MIA', 'Normal'] |
| C3 sample_key 全局唯一（= dataset:sample_id） | ✅ | 重复=0；唯一 sample_id=83/131（跨数据集同名 48 个=配对病灶） |
| C4 配对患者 == 两 GEO 表患者交集（且 ≥ 23 例） | ✅ | 23 例=['P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24', 'P25']；空间独有（无 snRNA，不入配对）：['P1', 'P2'] |
| C5 空转 LNM 未伪造 | ✅ | 无 LNM 行 |
| C6 样本数与 GEO 全量一致（{'GSE308103': 75, 'GSE307534': 56}） | ✅ | 实测={'GSE308103': 75, 'GSE307534': 56} |
| C7 本地名 vs GEO 标题 逐样本核对（无缺、无多、无不符） | ✅ | 131 样本全部一致 |
| C8 空转切片必需文件齐备（矩阵 + spatial 坐标/缩放） | ❌ | 1 张不完整，见报告 §切片完整性 |

## 计数

- 样本：**131**　患者：**25**　双模态配对：**23**

| 数据集 | 模态 | 样本数 | GEO 全量 | 患者数 |
| :-- | :-- | --: | --: | --: |
| GSE307534 | spatial | 56 | 56 | 25 |
| GSE308103 | snRNA | 75 | 75 | 23 |

## 配对患者核验

| 患者 | 模态 | 样本数 | 阶段 |
| :-- | :-- | --: | :-- |
| ⚠️ P1 | spatial | 2 | AAH;IAC |
| ✅ P10 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ✅ P11 | snRNA;spatial | 3 | AAH;IAC;Normal |
| ✅ P12 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P13 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ✅ P14 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P15 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ✅ P16 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P17 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P18 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ✅ P19 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ⚠️ P2 | spatial | 2 | AAH;IAC |
| ✅ P20 | snRNA;spatial | 3 | AAH;IAC;Normal |
| ✅ P21 | snRNA;spatial | 5 | AIS;IAC;Normal |
| ✅ P22 | snRNA;spatial | 4 | AAH;AIS;IAC;Normal |
| ✅ P23 | snRNA;spatial | 5 | AIS;IAC;Normal |
| ✅ P24 | snRNA;spatial | 3 | AAH;IAC;Normal |
| ✅ P25 | snRNA;spatial | 3 | AAH;IAC;Normal |
| ✅ P3 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P4 | snRNA;spatial | 6 | AAH;IAC;Normal |
| ✅ P5 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P6 | snRNA;spatial | 3 | AAH;IAC;Normal |
| ✅ P7 | snRNA;spatial | 4 | IAC;Normal |
| ✅ P8 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P9 | snRNA;spatial | 4 | AAH;AIS;IAC;Normal |

## GEO 交叉核验

两个数据集**逐 GSM** 与 GEO 样本标题核对：**无缺、无多、无不符**。

- 权威口径：`00_ingest/geo_metadata/{GSE308103,GSE307534}_samples.tsv`（由 NCBI E-utilities `esummary` 取得，含 GSM、患者号、分期、病灶序号、GEO 原始标题）。
- 交叉核验内容：本地解压目录名/文件名推导的 `(patient, stage_token)` 是否与 GEO 标题一致；本地是否存在 GEO 未登记项或缺失项。

| 数据集 | GEO 样本 | 本地命中 | 名不符 | 缺 | 多 |
| :-- | --: | --: | --: | --: | --: |
| GSE307534 | 56 | 56 | 0 | 0 | 0 |
| GSE308103 | 75 | 75 | 0 | 0 | 0 |

## 切片完整性

判据：每张空转切片必须含 `filtered_feature_bc_matrix/{matrix.mtx,barcodes.tsv,features.tsv}.gz` + `spatial/scalefactors_json.json` + `spatial/tissue_positions(.csv|_list.csv)`。

发现 **1** 张不完整：

- ⚠️ GSE307534 GSM9226176 切片文件缺失 ['spatial/scalefactors_json.json', 'spatial/tissue_positions.csv 或 spatial/tissue_positions_list.csv']（root=GSE307534/extracted/GSM9226176_P4_AAH-1/P4_AAH2）

> 处置：**重下该 GSM**（NCBI 端限速易致 tar 截断）；未补齐前 M0 不过门。

## 病灶序号（`lesion_ordinal`）

`lesion_ordinal=2` 表示**同一患者的第二个独立病灶切片**（GEO 标题作 `second … of patient N`），与首个病灶**分期相同但解剖独立**，**不得合并**：

| 数据集 | GSM | 患者 | token | GEO 标题 |
| :-- | :-- | :-- | :-- | :-- |
| GSE307534 | GSM9226176 | P4 | `AAH-1` | Lung tissue with second AAH of patient 4 [ST] |
| GSE307534 | GSM9226183 | P7 | `LUAD-1` | Lung tissue with second LUAD of patient 7 [ST] |
| GSE307534 | GSM9226212 | P21 | `AIS-1` | Lung tissue with second AIS of patient 21 [ST] |
| GSE307534 | GSM9226218 | P23 | `AIS-1` | Lung tissue with second AIS of patient 23 [ST] |
| GSE308103 | GSM9237905 | P4 | `Normal1` | Second normal lung tissue of patient 4 [snRNA-seq] |
| GSE308103 | GSM9237907 | P4 | `AAH1` | Lung tissue with second AAH of patient 4 [snRNA-seq] |
| GSE308103 | GSM9237917 | P7 | `LUAD1` | Lung tissue with second LUAD of patient 7 [snRNA-seq] |
| GSE308103 | GSM9237960 | P21 | `AIS1` | Lung tissue with second AIS of patient 21 [snRNA-seq] |
| GSE308103 | GSM9237968 | P23 | `AIS1` | Lung tissue with second AIS of patient 23 [snRNA-seq] |

## 来源与哈希

- `paired_samples.csv` `9c382df9210a2937c1e270452b6c4d74dd6ec53fccdc8833963c2824c0556f65`
- `paired_patients.csv` `c5f76baaee399677bd8f6f0747757b2612de038b2db5a21922a0ddcddbe26aeb`
- `paired_source_files.csv` `0c5e761c84a31da2fa2b880e8dc589be42898738ad7fde00ebea0945f443fdba`

## 说明
- **唯一键 = `sample_key`（`<dataset>:<sample_id>`）**。`sample_id` 跨数据集**不唯一**（同一患者同一病灶的 snRNA 与空间切片同名，如 `P3_LUAD`）→ 下游 join **必须用 `sample_key`**。
- 分期经 `cohort_registry.resolve_stage` **严格映射**；**原始 GEO 标签保留在 `stage_token`**（如 `LUAD` → 归一化 `IAC`，token 仍记 `LUAD`，便于回溯）。未知 token 直接 raise（无静默默认）。两侧命名惯例不同（snRNA 第二病灶用 `Normal1/AAH1/LUAD1/AIS1`，空间用 `Normal/AAH-1/LUAD-1/AIS-1`），**均按 GEO 原始文件名保真记录**，不强行统一。
- 患者身份取自 **GEO 样本标题**（`… of patient N`）；两数据集的 `P*` 编号即 GEO 标题中的患者号（GEO 无独立 patient 字段）—— 配对关系据此建立。
- **空转切片不靠目录名定位**：解压层级因 GEO 打包方式而异，一律由 `filtered_feature_bc_matrix/` 实际位置判定，候选 ≠ 1 即 raise（见 C8）。
- 源文件 SHA-256 未逐个计算（矩阵较大）—— 完整性以 M1 的 per-cell 产物哈希为准。
- 空转 LNM 暂缺（无合法 LUAD 数据）→ 不产出 LNM 行。
