# LUAD v2 · 计划指导与严格检查点 (PLAN & STRICT CHECKPOINTS)

> **性质**：本文件是**重做的计划指导（planning guide）**。每个里程碑设**硬性过门条件（gate）**，
> **不过门不得启动下一步**。
> **配套**：[`docs/WHITEPAPER.md`](docs/WHITEPAPER.md)（技术路线）· [`docs/PROJECT_SUMMARY.md`](docs/PROJECT_SUMMARY.md)（事实总纲）· [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md)（参数出处）。

---

## 0. 进度重置声明（PROGRESS RESET）
- 旧工程 `luad_invasion` 的**一切结果性产物作废**（图表/数值/表/TMB/CMap/对接），原因见 `PROJECT_SUMMARY`。
- **本工程进度归零**，从 `M0` 重新开始。旧目录仅作**只读数据源与历史参考**。
- 进度看板见本文件 §9（当前：全部 `0%`）。

---

## 1. 铁律（每步硬约束）
1. 数据身份以 **GEO/GSA** 为准；**分期无静默默认**（代码里不得出现 `.get(x,'IAC')` 之类回退）。
2. **恶性标签须 CNV 证真**；**双体用 scDblFinder**。
3. **无真实来源 = 不计算**（宁可拒绝/报缺，不伪造、不硬编码、不 np.random）。
4. 产物**可复现 + 有哈希**；`patient_id`（真患者）与 `sample_id`（组织/切片）**分层**。
5. **SCMG 分支不掺传统算法**；标准分支与 SCMG 分支**并行**，最后对照。

---

## 2. 里程碑与硬性检查点（Gates）

### M0 · 输入冻结门 (Input Freeze)
- **做**：按 [`00_ingest/cohort_registry.py`](00_ingest/cohort_registry.py) 纳入三 scRNA 队列（GEO 真值分期 + 分层）。
- **过门条件（全满足）**：
  1. per-cell 表字段齐全：`cell_barcode, patient_id, sample_id, stage, dataset`；
  2. **阶段计数与 GEO 真值一致**：GSE189357 `IAC=TD1/2/9`、`MIA=TD3/4/6`、`AIS=TD5/7/8`，**AAH 计数=0**；
  3. 代码审查**无静默默认**（grep 通过）；
  4. 产出**冻结清单 + SHA-256 + 校验报告**。
- **不过门 → 停**（不得进 M1）。

### M1 · QC / 双体门
- **做**：显式 QC；**scDblFinder 逐样本**。
- **过门**：报告 QC 前后细胞数、双体率、参数（标出处，见 PARAMETERS）；**无启发式替代**。

### M2 · 恶性证真门 (CNV)
- **做**：`inferCNV`/`CopyKAT` 判 CNV。
- **过门**：报告 CNV 阳性细胞数 + 与经典标记一致性；**恶性克隆的下游使用以此为准**。**不过门不得产出"恶性克隆"。**

### M3 · 整合门（双分支）
- **做**：**A 标准/scVI-scANVI**（modality 当 batch）与 **B 纯 SCMG**（zero-shot）并行。
- **过门**：以 **scIB 口径（LISI/kBET + 生物保守）**给数；证明**有分辨力**（各阶段可分），非"全糊一块"。

### M4 · 跨模态 AAH 门（**五判据缺一不可**）
- ① 各阶段**可分辨**；② 重叠阶段 sn↔sc **一致**；③ 平台偏移 **δ(stage) 稳定**；④ **sn 内部配对**（AAH vs 同患者 Normal/AIS）同向；⑤ AAH 身份 **CNV/标记**证真。
- **不过门 → 不下"AAH 结论"**（只能标为假说）。
- 数据：暂 `GSE308103`(sn)；`HRA001130`(sc) 留接口待申请。

### M5 · 空间解卷积门
- **做**：RCTD，reference 用 M3 冻结的签名。
- **过门**：**Σ_t P(t|s)=1**（容差内）、非负、被拒 spot 数报告；签名与 reference **哈希对齐**。

### M6 · 生态位门
- **做**：SpaGCN + GMM/BIC（无预设聚类数）+ Squidpy 置换。
- **过门**：**无硬编码聚类数/标签**；置换有 `n_perms` 记录（见 PARAMETERS）。

### M7 · 靶点门（传统轨道）
- **做**：患者级 Pseudobulk DESeq2；TCGA Cox/KM；**CMap 仅在有真实 LINCS 数据时执行**。
- **过门**：DESeq2 设计/截断有出处；**CMap 无真实数据 → 拒绝产出表**。

### M8 · 对接门
- **做**：fpocket + AutoDock Vina（靶点来自 M7 真实结果）。
- **过门**：靶点/结构文件**真实存在**；参数标出处；**未计算不虚构**。

---

## 3. 通用检查点（每步都查）
- [ ] 无 `np.random` / 写死数值 / 伪曲线；
- [ ] 无静默默认（回退到某阶段/标签）；
- [ ] 脚本可复现（确定性 + 种子显式）；
- [ ] 产物有哈希；输入输出可溯源；
- [ ] 参数在 `docs/PARAMETERS_AND_SOURCES.md` 有出处。

---

## 4. 双分支架构（贯穿 M3–M8）
- **标准/主流分支**：scVI/scANVI、Harmony、Seurat、RCTD、SpaGCN、Squidpy、DESeq2、CellRank…
- **纯 SCMG 分支**：zero-shot 整合 → 流形 → 状态 → 轨迹（条件扩散）→ 因果（`CausalGenePredictor`）。
- 两分支**对照**（scIB 口径）；工具源码见 [`tools/`](tools/)。

---

## 5. 数据来源与获取
- scRNA：GSE131907 / GSE189357 / GSE148071（**无 AAH**）；
- AAH：`GSE308103`(sn，暂用) / `HRA001130`(sc，受控，[接口](00_ingest/hra001130_interface.py))；
- 空间：GSE307534 / GSE189487 / GSE190811（核心 6 阶段切片见 P0 禁令）。

---

## 6. 里程碑进度看板
| 里程碑 | 状态 | 过门 |
| :--- | :---: | :---: |
| M0 输入冻结 | ⬜ 0% | ☐ |
| M1 QC/双体 | ⬜ 0% | ☐ |
| M2 CNV 证真 | ⬜ 0% | ☐ |
| M3 整合(双分支) | ⬜ 0% | ☐ |
| M4 跨模态 AAH | ⬜ 0% | ☐ |
| M5 空间解卷积 | ⬜ 0% | ☐ |
| M6 生态位 | ⬜ 0% | ☐ |
| M7 靶点(传统轨道) | ⬜ 0% | ☐ |
| M8 对接 | ⬜ 0% | ☐ |

> 每过一个门 → 更新本表 + 记一笔"过门证据"（产物路径 + 哈希）。
