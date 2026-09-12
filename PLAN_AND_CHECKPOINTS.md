# LUAD v2 · 计划指导与严格检查点 (PLAN & STRICT CHECKPOINTS)

> **性质**：本文件是**重做的计划指导（planning guide）**。每个里程碑设**硬性过门条件（gate）**，
> **不过门不得启动下一步**。
> **配套**：[`docs/WHITEPAPER.md`](docs/WHITEPAPER.md)（技术路线）· [`docs/PROJECT_SUMMARY.md`](docs/PROJECT_SUMMARY.md)（事实总纲）· [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md)（参数出处）。

---

## 0. 进度重置声明（PROGRESS RESET）
- 旧工程 `luad_invasion` 的**一切结果性产物作废**（图表/数值/表/TMB/CMap/对接），原因见 `PROJECT_SUMMARY`。
- **本工程进度归零**，从 `M0` 重新开始。旧目录仅作**只读数据源与历史参考**。
- 进度看板见本文件 §6（M0 已完成；M-1 起为 0%）。

---

## 1. 铁律（每步硬约束）
1. 数据身份以 **GEO/GSA** 为准；**分期无静默默认**（代码里不得出现 `.get(x,'IAC')` 之类回退）。
2. **恶性标签须 CNV 证真**；**双体用 scDblFinder**。
3. **无真实来源 = 不计算**（宁可拒绝/报缺，不伪造、不硬编码、不 np.random）。
4. 产物**可复现 + 有哈希**；`patient_id`（真患者）与 `sample_id`（组织/切片）**分层**。
5. **SCMG 分支不掺传统算法**；标准分支与 SCMG 分支**并行**，最后对照。
6. **SCMG 分支能力边界**（2026-09-12 核实）：只做 **zero-shot 跨数据集 scRNA 整合 + 流形 + 状态刻画**；
   **不得**输出「状态逆转 / 逆转因子 / 因果基因」（源码无此能力，参照流形无肿瘤态）。靶点为**候选**，不得称因果。

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

### M-1 · 计划整改门 (Plan Remediation) —— **先于 M1，必做**

> **背景**：2026-09-12 对白皮书做**逐阶段可行性/严谨性审计**（数据身份 + 8 份方法学核查），
> 发现 **3 个 🔴（方法/数据不成立）+ 7 项 🟡**。**不过此门，不得进入 M1 及其后任何阶段。**

**A. 数据层**
- **A1** 空间图谱限定 **Normal→AAH→AIS→MIA→IAC**（全部 GSE307534，同库同平台）；**删除 GSE190811**（经 GEO 核实为**乳腺癌**，非 LUAD；原 GSM5732148 在该库不存在）。规则文档已改（[`docs/spatial_cohort_and_figure_prohibitions.md`](docs/spatial_cohort_and_figure_prohibitions.md)），白皮书待改。
- **A2** LNM 暂**只保留单细胞层**（GSE131907 `mLN`）；**空转 LNM 待真实 LUAD 数据**（公开库暂无）。
- **A3** 可选：`GSE305258`（ALK+ NSCLC 淋巴结/脑转移空间，10 LNT，**GeoMx ROI 非 Visium**）仅作 **LNM 正交验证**，**不并入主 Visium 矩阵**。
- **A4** **永久禁用 WES/TMB**（TCGA 突变文件为 0 字节）。

**B. 工具/环境**
- 安装 **CopyKAT**（`MCMCpack`+`transport` from CRAN；纯 R，无 root）→ M2 **主力**。
- **infercnv 降级为可选交叉验证**（官方 README 已声明 *no longer supported*；需系统 JAGS；本机 Bioconductor 不可达）。
- **升级 scvi-tools 至 1.5.x**（现装 0.15.5，2022 版）→ 方可使用 **SysVI / scArches(scANVI surgery)**。
  ⚠️ **不能在共享环境做**（Py3.8 + root 库；见 §5b）→ 必须建 **隔离 env（micromamba, Py≥3.10）**；否则退用 0.15.5 能力并**如实标注**（无 SysVI/scArches surgery）。
- 安装 **spaGCN**、**harmonypy**、**plip(py)**（权重已缓存 1.2 GB）；**CellCharter**（若采用）。
- **GROMACS 缺失** → 装，或**砍掉 Stage-8 的 100 ns MD**（仅保留对接）。

**C. 方法学修正**（写入下方 M3/M5/M6/M7 各门）
- **M3**：弃用「modality 当 batch」（sc↔sn 是 *system* 效应，且 modality 与 dataset 共线）→ **scVI(batch=dataset)** 建 scRNA 图谱 + **scArches/scANVI 映射** sn RNA；scIB 用**完整 panel**。
- **M5**：弃用自归一化 Σ=1 作为门（full 模式无 reject，该门为空）→ 门改到 **RCTD 原生输出 + 显式 QC**；参考模态与切片匹配。
- **M6**：GMM/BIC 改**稳定性/共识选择**；Fisher 改报**效应量**；Squidpy 统计改**类别标签 + 经验 p**。
- **M7**：CMap 指标纠正为 **NCS**，或直接**砍掉**（无 LINCS 数据不产出）。
- **SCMG 分支收缩**：仅做 **zero-shot 跨数据集 scRNA 整合 + 全局流形 + 状态刻画**，与标准分支 **scIB 对照**。
  **删除所有「状态逆转 / 逆转因子 / CausalGenePredictor 因果」表述**（经核实：SCMG 无此能力；参照流形无肿瘤态）。
  靶点候选改由**标准分支**产出（SCENIC+/regulon + CellRank 命运关联 + TCGA 生存），用词限定为 **「候选调控因子/靶点」**，不得称因果。

**D. 文档修正**
- 白皮书：删作废结果（Figure 1–5「已完成」、14 克隆、107,796 spots、`Z<-5.8`、316,689 等）；LNM 改 pending；删 SCMG 逆转表述；修 CellRank 矛盾。
- PARAMETERS：修 `infercnv denoise` 默认（**FALSE**，非 TRUE）、Tau 定义矛盾（**0–1** vs ±100）；QC 阈值补文献或标「本项目约定」。

**过门条件（全满足）**：
1. 必需工具**安装且 smoke-test 通过**（CopyKAT 对 1 样本跑通；`scvi.external.SysVI` / scArches 可 import）；
2. 白皮书 / PLAN / PARAMETERS 修订**落盘**，且对作废数值（`14 克隆`、`107,796`、`Z<-5.8`、`316,689`、`PT_3_LNM`/`GSM5732148`）grep 无残留（历史勘误记录除外）；
3. GSE190811 在代码 / 文档 / 规则中**全部清除**（除勘误说明）；
4. M3 / M5 / M6 / M7 门**重写落盘**；
5. **无 🔴 遗留**。

### M1 · QC / 双体门
- **做**：显式 QC；**scDblFinder 逐样本**。
- **过门**：报告 QC 前后细胞数、双体率、参数（标出处，见 PARAMETERS）；**无启发式替代**。

### M2 · 恶性证真门 (CNV)
- **做**：**`CopyKAT` 为主力**（逐样本，**禁止合池**）；**`infercnv` 仅作 5–10k 细胞/样本的可选交叉验证**（官方已停维护）。
- **过门**：报告 CNV 阳性细胞数 + 与经典标记一致性；**恶性细胞的下游使用以此为准**。**不过门不得产出"恶性克隆"。**
- **注**：GSE189357 / GSE148071 无正常样本 → 参考用同样本免疫/基质细胞；snRNA（GSE308103）CNV 需放宽参数并谨慎解读。

### M3 · 整合门（双分支）
- **做**：**A 标准分支**：**scVI（`batch=dataset`）** 建 scRNA 图谱 → **scArches/scANVI 映射 sn RNA**（非 zero-shot）。
  **B SCMG 分支**：**zero-shot 跨数据集 scRNA 整合 + 全局流形 + 状态刻画**（**不掺传统算法**；**不输出逆转/因果**）。
- **过门**：以 **scIB 完整口径**给数——批去除（kBET + iLISI + graph-connectivity + PCR）**与** 生物保守（cLISI/ARI/NMI/ASW）**并报**；证明**有分辨力**（各阶段可分），非"全糊一块"。
- **禁用**：`modality 当 batch`（sc↔sn 为 system 效应且与 dataset 共线）；仅凭 iLISI↑ 不得判为整合成功（可被过度整合刷高）。

### M4 · 跨模态 AAH 门（**五判据缺一不可**）
- ① 各阶段**可分辨**；② 重叠阶段 sn↔sc **一致**；③ 平台偏移 **δ(stage) 稳定**；④ **sn 内部配对**（AAH vs 同患者 Normal/AIS）同向；⑤ AAH 身份 **CNV/标记**证真。
- **不过门 → 不下"AAH 结论"**（只能标为假说）。
- 数据：暂 `GSE308103`(sn)；`HRA001130`(sc) 留接口待申请。

### M5 · 空间解卷积门
- **做**：RCTD（`doublet_mode='full'`，Visium 推荐），reference 用 M3 冻结的签名。
- **过门**：以 **RCTD 原生输出**为准——权重非负、**每 spot 权重分布合理**、**参考模态与切片匹配**（FFPE↔sn，先剔除跨模态 DEG）、报告**低质量/被剔除 spot 数**（由 `UMI_min`/`counts_MIN` 显式阈值定义，**不用** self-imposed 的 Σ=1 充当门）；签名与 reference **哈希对齐**。
- **禁用**：把 `normalize_weights()` 后的 Σ=1 当作 RCTD 的性质（那是项目自己施加的）；`full` 模式**没有** reject 类别。

### M6 · 生态位门
- **做**：SpaGCN 结构域分割 → 邻域聚合 → **数据驱动聚类（稳定性/共识选择，如 CellCharter 或 bootstrap 稳定性；不用单一 BIC 极小值）** → Squidpy 空间统计。
- **过门**：**无硬编码聚类数/标签**；簇数由**稳定性/共识**决定并留档；空间统计用**类别标签**（不是解卷积比例），报 **z-score + 经验 p（(b+1)/(n+1)）**，置换 `n_perms` 有记录；富集检验报**效应量**（Fisher FDR 受组成性与大 N 灌水）。
- **禁用**：把 `nhood_enrichment` 的输出写成 `P<0.001`（该函数**不返回 p 值**）；对 `co_occurrence` 声称做过置换（它**不做**）；距离步长小于 Visium ~100 µm 点距。

### M7 · 靶点门（传统轨道）
**M7a · 生态位预后特征**：CellCharter 生态位 → **患者/切片级聚合**出签名 → TCGA-LUAD bulk 打分（ssGSEA/解卷积）→ **多变量 Cox（校正分期/年龄/性别）+ KM**。
- **过门**：签名在**患者级**聚合产生（不用 spot 级 DE）；防过拟合（惩罚/交叉验证）；措辞 = **预后关联**（`niche poor-prognosis signature`）。空间仅 11–25 例 → **不主张患者亚型分型**。

**M7b · 候选靶点池（遗传统计锚定）**：候选来源 = 恶性程序/regulon **＋** 生态位签名 → **cis-MR + coloc**（ILCCO/TRICL LUAD GWAS × 肺 eQTL）→ **Open Targets 可成药性 + DepMap 选择性依赖 + 临床期药物匹配**。
- **操作手册**：[`docs/M7B_MR_COLOC_TARGET_ANCHORING.md`](docs/M7B_MR_COLOC_TARGET_ANCHORING.md)
- **过门**：每个候选给出 `genetic_support ∈ {supported, not_supported, not_testable}`（阈值见手册，**不得事后调参**）；`not_testable` **如实标缺，不替代**；措辞 = **genetically supported candidate target**，**不得称因果**。

**M7c · CMap（可选）**：仅当获得**真实 LINCS 数据**且用对指标（**NCS**，非负 Tau）时才执行；否则**拒绝产出表**（现无 LINCS 数据 → 默认不产出）。

- **总过门**：DESeq2 设计/截断有出处；靶点措辞限定；MR 输入/输出有哈希与来源记录。

### M8 · 对接门
- **做**：fpocket + AutoDock Vina（靶点来自 M7 真实结果）；**MD 仅在 GROMACS 安装后才做**。
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
- **标准/主流分支**：scVI（`batch=dataset`）/scANVI + scArches、Harmony、Seurat、RCTD、SpaGCN、Squidpy、DESeq2、CellRank、SCENIC+…
- **纯 SCMG 分支**：**zero-shot 跨数据集 scRNA 整合 → 全局流形 → 细胞状态刻画**（**不掺传统算法**）。
  ~~条件扩散轨迹 → CausalGenePredictor 因果~~ —— **已删除**（能力不存在，见铁律 6）。
- 两分支**对照**（scIB 完整口径）；工具源码见 [`tools/`](tools/)。
- 靶点/调控因子候选**只由标准分支产出**，且措辞为「候选」。

---

## 5. 数据来源与获取
- scRNA：GSE131907 / GSE189357 / GSE148071（**无 AAH**）；
- AAH：`GSE308103`(sn，暂用) / `HRA001130`(sc，受控，[接口](00_ingest/hra001130_interface.py))；
- 空间：**GSE307534**（主，含 AAH/AIS/MIA/LUAD，25 患者）/ **GSE189487**（冻存验证）；**LNM 空转暂缺**（~~GSE190811~~ 经核实为乳腺癌，已废；待真实 LUAD 数据；可选 `GSE305258` 仅作正交验证）。
- TCGA-LUAD：仅表达 + 临床（**突变文件 0 字节，WES/TMB 永久禁用**）。

---

## 5b. 环境约束与最简决定（Environment Constraints，2026-09-12 核实）

> 本机环境**硬约束**，决定了工具选型与"哪些方法根本装不上"。**凡与之冲突的方案一律作废。**

1. **无 GPU**（无 `/dev/nvidia*`，`torch.cuda.is_available()=False`）→ GPU 依赖方法（Boltz-2、cell2location、CellCharter、SysVI…）**CPU-only 或 ~42 万细胞下不可行**。
2. **仅 Python 3.8.10、无 conda** → 现代 DL 栈装不上：
   - `cellcharter` 需 ≥3.9；`scvi-tools` 1.5.x（SysVI / scArches surgery）需 ≥3.10 → **在共享环境不可装**。
   - **解**：需要现代 DL 时用**隔离环境**（micromamba from conda-forge）；PyPI 仅经镜像（`https://pypi.tuna.tsinghua.edu.cn/simple`），conda-forge 直连。
3. **共享库 root 属主**（`/usr/local/lib/R/site-library`、`/usr/local/lib/python3.8/dist-packages`）→ **禁止把社区 DL 包 pip 进共享环境**：
   - 曾装 `cellcharter` **静默把 torch 降到 1.12.1**，搞坏 `scvi`/`pytorch-lightning`（已手工回滚 torch→2.4.1+cu118、pytorch-lightning→1.5.10.post0、torchmetrics→0.7.3）。
   - `scvi-tools 0.15.5` **保持**（pin `pytorch-lightning>=1.5,<1.6`）；装包一律进**个人库**（`~/.local/lib/python3.8/site-packages`、`~/R/.../4.2`）。
4. **网络**：`http(s)_proxy=127.0.0.1:7890` 指向**死端口**；出网须**绕代理**（`curl --noproxy '*'` / Python `ProxyHandler({})`）。`github.com` 被墙、`api.github.com`/`codeload.github.com` 通（故 `remotes::install_github` 在**清空代理变量**下可用）。
5. **R 可装**：**CopyKAT**（先装 `RcppEigen`(Eigen 4.0) 到个人库，再 `transport`）、`coloc`、`ieugwasr`、`TwoSampleMR`。

**对本计划的直接影响**：
- M-1 §B 里"升级 scvi-tools 至 1.5.x"**不能在共享环境做** → 必须建**隔离 env（micromamba）**，或改用 **scvi 0.15.5 的可用能力**（无 SysVI/scArches surgery）并如实标注。
- M2 主力 **CopyKAT** 与 M6 的 **CellCharter**：CopyKAT 可装（见上）；**CellCharter 需隔离 env**。
- 任何"GPU 加速"表述须删或降级为 CPU。

---

## 6. 里程碑进度看板
| 里程碑 | 状态 | 过门 |
| :--- | :---: | :---: |
| M0 输入冻结 | ✅ 100% | ✅ |
| M-1 计划整改 | ⬜ 0% | ☐ |
| M1 QC/双体 | ⬜ 0% | ☐ |
| M2 CNV 证真 | ⬜ 0% | ☐ |
| M3 整合(双分支) | ⬜ 0% | ☐ |
| M4 跨模态 AAH | ⬜ 0% | ☐ |
| M5 空间解卷积 | ⬜ 0% | ☐ |
| M6 生态位 | ⬜ 0% | ☐ |
| M7 靶点(传统轨道) | ⬜ 0% | ☐ |
| M8 对接 | ⬜ 0% | ☐ |

> 每过一个门 → 更新本表 + 记一笔"过门证据"（产物路径 + 哈希）。

### 过门证据 (Gate Evidence)

**M0 输入冻结 — ✅ PASS (2026-09-12)**
- 脚本：[`00_ingest/01_load_cohorts.py`](00_ingest/01_load_cohorts.py)  · registry sha256 `7ad44013fd013dedeb80e77be23838b03d9a22cfab1c6d7860743ebf8d37b4cb`
- GEO 权威元数据（冻结）：`00_ingest/geo_metadata/`（`fetch_geo_metadata.py` 拉取；含 `patient id`/histology/origin）
  - GSE131907 `8d94dcc62bd31ba6…` · GSE189357 `0f2ac9213286f0fe…` · GSE148071 `45d2fa96262eafc2…`
- 产物（`results/00_ingest/`）：
  - `frozen_per_cell.csv.gz`（420,766 细胞；109 样本；95 患者）sha256 `703c5f03562d031a6aea96a68dea28f312733330db76fa0f0ed361a61657f7d2`（确定性 gzip，跨运行稳定）
  - `frozen_samples.csv` / `frozen_patients.csv` / `frozen_source_files.csv` / `frozen_manifest.json`
  - 校验报告：[`results/00_ingest/M0_validation_report.md`](results/00_ingest/M0_validation_report.md)（17 项检查 PASS + 自审）
- 阶段计数（实测）：GSE189357 IAC=TD1/2/9、MIA=TD3/4/6、AIS=TD5/7/8（**AAH=0**）；GSE131907 nLN→Normal_LN（非 LNM）、脑转/胸水单列；GSE148071=Adv_NSCLC。
- **GEO 交叉核验**：GSE131907 `tissue origin` 逐样本 == registry（58/58）；GSE189357 `histolgical type` 逐样本 == registry（9/9，独立确认无 AAH）。
- **患者身份（GEO 权威）**：GSE131907 = **44 患者**（取自 GEO `patient id`，非样本名尾号）；GSE189357=9、GSE148071=42（1 样本=1 患者假设，GEO 无 patient id 字段）。
