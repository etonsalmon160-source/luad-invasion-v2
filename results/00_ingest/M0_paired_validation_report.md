# M0 输入冻结（配对数据集）· 校验报告

- schema `M0-paired-1` · 2026-09-12T12:56:59Z · raw `/home/eto/luad_invasion/data`
- registry sha256 `6843e4e86cdaef11323663c428a7c49bfb86b518edf1d4c926db7ef983fce5f6`

## 过门判定：**✅ PASS**

| 检查 | 结果 | 详情 |
| :-- | :--: | :-- |
| C1 两数据集均在册 | ✅ | ['GSE307534', 'GSE308103'] |
| C2 分期 ∈ 词表 | ✅ | ['AAH', 'AIS', 'IAC', 'MIA', 'Normal'] |
| C3 sample_key 全局唯一（= dataset:sample_id） | ✅ | 重复=0；唯一 sample_id=80/94（跨数据集同名 14 个=配对病灶） |
| C4 配对患者 == 登记表 9 例 | ✅ | 实测=['P10', 'P13', 'P15', 'P18', 'P21', 'P22', 'P25', 'P3', 'P4'] |
| C5 空转 LNM 未伪造 | ✅ | 无 LNM 行 |

## 计数

- 样本：**94**　患者：**25**　双模态配对：**9**

| 数据集 | 模态 | 样本数 | 患者数 |
| :-- | :-- | --: | --: |
| GSE307534 | spatial | 19 | 11 |
| GSE308103 | snRNA | 75 | 23 |

## 配对患者核验

| 患者 | 模态 | 样本数 | 阶段 |
| :-- | :-- | --: | :-- |
| ⚠️ P1 | spatial | 2 | AAH;IAC |
| ✅ P10 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ⚠️ P11 | snRNA | 3 | AAH;IAC;Normal |
| ⚠️ P12 | snRNA | 3 | AIS;IAC;Normal |
| ✅ P13 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ⚠️ P14 | snRNA | 3 | AIS;IAC;Normal |
| ✅ P15 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ⚠️ P16 | snRNA | 3 | AIS;IAC;Normal |
| ⚠️ P17 | snRNA | 3 | AIS;IAC;Normal |
| ✅ P18 | snRNA;spatial | 3 | IAC;MIA;Normal |
| ⚠️ P19 | snRNA | 3 | AIS;IAC;Normal |
| ⚠️ P2 | spatial | 2 | AAH;IAC |
| ⚠️ P20 | snRNA | 3 | AAH;IAC;Normal |
| ✅ P21 | snRNA;spatial | 4 | AIS;IAC;Normal |
| ✅ P22 | snRNA;spatial | 4 | AAH;AIS;IAC;Normal |
| ⚠️ P23 | snRNA | 4 | AIS;IAC;Normal |
| ⚠️ P24 | snRNA | 3 | AAH;IAC;Normal |
| ✅ P25 | snRNA;spatial | 3 | AAH;IAC;Normal |
| ✅ P3 | snRNA;spatial | 3 | AIS;IAC;Normal |
| ✅ P4 | snRNA;spatial | 6 | AAH;IAC;Normal |
| ⚠️ P5 | snRNA | 3 | AIS;IAC;Normal |
| ⚠️ P6 | snRNA | 3 | AAH;IAC;Normal |
| ⚠️ P7 | snRNA | 3 | IAC;Normal |
| ⚠️ P8 | snRNA | 3 | AIS;IAC;Normal |
| ⚠️ P9 | snRNA | 4 | AAH;AIS;IAC;Normal |

## 来源与哈希

- `paired_samples.csv` `db5dcf7f4750aa2381dba7fd9cd21063734fd098da907d23e5feaecf03aabe3a`
- `paired_patients.csv` `c9a0fc156365e49a598ac75d10b2f47be3272557e3e2bf6c14c97f85a47f0403`
- `paired_source_files.csv` `e4d7722a155857fa6199fbda24513b3647b7bfea3cf34f91ef83ba18d657ed25`

## 说明
- **唯一键 = `sample_key`（`<dataset>:<sample_id>`）**。`sample_id` 跨数据集**不唯一**（同一患者同一病灶的 snRNA 与空间切片同名，如 `P3_LUAD`）→ 下游 join **必须用 `sample_key`**。
- 分期经 `cohort_registry.resolve_stage` **严格映射**；**原始 GEO 标签保留在 `stage_token`**（如 `LUAD` → 归一化 `IAC`，token 仍记 `LUAD`，便于回溯）。未知 token 直接 raise（无静默默认）。
- 患者身份取自 **GEO 样本标题**（`… of patient N`）；两数据集的 `P*` 编号即 GEO 标题中的患者号
（GEO 无独立 patient 字段）—— 配对关系据此建立，见 §GEO 交叉核验。
- GSE307534 的 GEO 全量为 56 样本；**本地已解压**的计入本表（其余待下载）。
- 源文件 SHA-256 未逐个计算（矩阵较大）—— 完整性以 M1 的 per-cell 产物哈希为准。
- 空转 LNM 暂缺（无合法 LUAD 数据）→ 不产出 LNM 行。
