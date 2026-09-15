# 00_ingest —— 权威输入层

> **范围（2026-09-12 收窄）**：仅**两个配对数据集**——`GSE308103`(snRNA) + `GSE307534`(Visium 空间)。
> 三个 scRNA 队列（GSE131907/189357/148071）已移出范围，存档于 `/home/eto/luad_invasion/luad_v2_out_of_scope/`。

## 为什么是这两个

同一研究、**同一批患者**（23 例配对：P3–P25）、**相邻切片**，且
**模态匹配**（空间是 FFPE，参考也是 FFPE 的 snRNA）——这使空间解卷积的参考最可信，
并支撑跨模态一致性判据（AAH vs 同患者 Normal/AIS）。

| 数据集 | 模态 | 角色 | 分期 |
| :--- | :--- | :--- | :--- |
| **GSE308103** | **snRNA**（细胞核；FFPE） | 单细胞**参考**；**唯一含 AAH** 的单细胞资源 | Normal / AAH / AIS / MIA / IAC |
| **GSE307534** | **Visium spot**（FFPE CytAssist；55 µm，**非单细胞**） | **空间图谱**（原位坐标）；**解卷积对象** | Normal / AAH / AIS / MIA / IAC |

> ⚠️ **Visium spot 是多细胞混合** → 必须用单细胞参考**解卷积**才能得到每个 spot 的细胞组成。
> "配对"指**同患者/同病灶**，不是"同一细胞测了两遍"；它让参考匹配，但不取消解卷积这一步。

## 铁律
1. 分期**无静默默认**（[`cohort_registry.py`](cohort_registry.py) 的 `resolve_stage` 未知 token 一律 raise）；
2. `patient_id`（患者）与 `sample_id`（切片/样本）分层；
3. **snRNA 的 QC 阈值不得照搬整细胞 scRNA**（核的 nCount 更低、mt% 更低）——先测后定；
4. 数据未到位 → 如实报缺，不造数据。

## 文件
- [`cohort_registry.py`](cohort_registry.py) —— 权威登记表（配对患者 / 分期 token / P0 禁用法）

> HRA001130（受控库）预留接口已移出范围 → `/home/eto/luad_invasion/luad_v2_out_of_scope/00_ingest/`。

## 关键事实（已核实）
- GSE308103：**79.8 万核 / 75 样本**（实测），median nCount 1,516、median pct_mt 0.6%；
- GSE307534：GEO 56 样本 / 25 患者；**本地已解压 56 张切片，完整覆盖 23 例配对患者**；
- LNM 空转暂缺（见 registry `LNM_STATUS`）；`GSE190811` 经核实为**乳腺癌**，已废。

## 下一步
- [x] `01_qc/00_metrics_gse308103.R`：先算指标、据实定阈值（核数据）
- [ ] `01_qc/01_qc_doublets_gse308103.R`：MAD 离群 + scDblFinder 逐样本（运行中）
