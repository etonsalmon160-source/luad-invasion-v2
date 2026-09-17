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

### 2.0 执行顺序（2026-09-16 按**传统单细胞流程**重构）

> **背景**：本文件原先把 **M2（CNV）排在 M3（降维聚类）之前**，那是把 CNV 当作"前置的主力证真步骤"。
> 传统单细胞流程的实际顺序相反：**先归一化→降维→聚类→marker 注释，再对所注释出的上皮簇做 CNA 细化**。
> 本次重构后 **M 编号降级为「模块标识符」**（`M2` 只表示"CNV 模块"，不再蕴含先做后做）；
> **执行顺序以下表为准**，GP 编号为稳定的检查点编号。

> 🔴 **2026-09-16 二次重构（用户决策「改用论文 QC 重建」+「对齐论文全套」）**：第 2–5 行全部改写。
> **GP0 用论文 QC 口径重做**（413,697 核，非原 648,945）；**GP4a–4c 换成源论文的 Seurat 配方**
> （SCTransform v2 → HVG 3000 → PCA 50 → **Harmony 为主口径** → FindClusters/Louvain）。
> 旧的 scran / HVG-2000 / PCA-30 / kNN-15 / Leiden 与 `{20,30,50}` ARI 护栏**全部作废**（登记见 PARAMETERS §M3-A.1）。
> 详见 [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md) §M3-A.0–A.3。

> 🔴 **2026-09-16 三次重构（用户决策：CNV 降为"上皮亚群恶性精判"）**：第 7–10 行改写，GP8 拆为 **GP8a/GP8b**。
> **CNV 的定位从"独立证真门"收缩为"给亚群贴恶性标签的精细判别器"**，因此**必须排在上皮亚聚类之后**
> （原排 GP6 之后，只能出逐细胞读数，落不到亚群上）。**免疫/基质/内皮谱系的亚聚类（GP8b）不再被 CNV 排期阻塞**。
> **GP1（CopyKAT 冒烟）随之下沉为 GP2 的执行细节**，不再单列。

| 序 | 检查点 | 内容 | 子门 |
| :---: | :--- | :--- | :--- |
| 1 | ~~M1~~ | QC / 双体（逐样本自适应 MAD） | ✅ 已过（🔶 降级为**敏感性臂**） |
| 2 | **GP0** | 表达对象重建 —— **🔴 2026-09-16 按论文 QC 重做**：`nFeature≥500 & nCount≥1000 & pct_mt≤20` ∩ scDblFinder 单细胞 → **413,697 核 × 18,069 基因** | 🔶 已重做；**2026-09-17 审计补做门校验（13 PASS / 0 FAIL）**，待你签字 |
| 3 | **GP4a** | **预处理（论文配方）**：`SCTransform(vst.flavor="v2", variable.features.n=3000)` | ✅ 全量已跑（2026-09-17） |
| 4 | **GP4b** | `RunPCA(npcs=50)` —— **无 ARI 护栏**（原护栏随旧网格作废） | ✅ 全量已跑（2026-09-17） |
| 5 | **GP4c** | 批次：**Harmony on `sample_id`, `dims.use=1:50`（主口径，论文配方）** / 不校正（**敏感性臂，非论文配方**） | ✅ **2026-09-17 用户裁定：采用论文 Harmony 主口径**（附不可归因限定） |
| 6 | **GP5** | 分辨率选择：**`{0.5,0.6,0.7,0.8}` × 5 种子** + 指标 **1/2/3/4**（含 **AAH 吸收护栏**）+ 预注册 `r*` 规则（**注释之前**人工签字） | ✅ **2026-09-17 签字：`r*=0.6`**；指标4 承认全量无效 |
| 7 | **GP6** | **双标准注释**（法则2 marker ⊕ CellTypist）+ Cohen's κ —— **只负责定出"哪些是上皮"** | ⬜ |
| 8 | **GP8a** | **上皮谱系亚聚类**（子集内重启全流程，口径同 GP4a–c）→ **上皮亚群** | ⬜ |
| 9 | **GP2** | **CNV 精判**：逐样本 CopyKAT，输入 = GP6 定出的**上皮细胞**；结论按 GP8a 的**上皮亚群**汇总 → 逐亚群恶性判定 | ⬜ |
| 10 | **GP8b** | 其余 **5 个谱系**亚聚类（非上皮，与 CNV 无关）+ 层级自洽 | ⬜ |
| 11 | **GP3 / GP7** | 指标后端冻结（GP3）→ scIB 双 panel 评估（GP7） | ⬜ |
| 12 | **GP9** | **SCMG 对照臂**（M3-B）：zero-shot 整合 + 全局流形 + 状态刻画，**不掺传统算法**（铁律 5） | ⬜ |
| 13 | M4… | 跨模态 AAH → M5 → M6 → M7 → M8 | ⬜ |

**关键次序约束**：
- **GP2（CNV）必须晚于 GP8a（上皮亚聚类）** —— 2026-09-16 第三次调整。CNV 在本项目里的**唯一用途**
  是**给上皮亚群精判恶性**，而"亚群"这个对象在 GP8a 之前**还不存在**。停在 GP6 之后跑，
  只能给出"逐上皮**细胞**"的读数，却无法把它归到亚群上，等于把最需要的那一步留给下游猜。
  ⇒ **顺序：GP6 定上皮 → GP8a 定亚群 → GP2 定亚群恶性**。
- **GP2 不得给出"恶性克隆"级结论** —— 核 UMI 低、含内含子/环境 RNA，CNV 只能支持
  **亚群级（非整倍体 vs 否）**判定。措辞上限见 §M2。
- **GP6 之前必须 GP5 签字** —— 分辨率不可事后拟合。
- **GP8a 与 GP8b 拆开的原因** —— 只有上皮**需要**等 CNV（R2 要求恶性须 CNV 证真）；
  免疫/基质/内皮谱系的亚聚类与 CNV 无关，**不得**被 CNV 的排期阻塞。
- **GP4c 主口径 = Harmony（论文配方）**，不校正臂仅为**稳健性对照**；**GP9 的 SCMG 是另一条分支**。
  三者不是同一层对照，报告时**不得并表**。
- **GP3 须在 GP7 之前** —— 单数据集无 ground truth 的**门规格变更**必须先在 GP3 落书面批准。

### M0 · 输入冻结门 (Input Freeze)
- **做**：按 [`00_ingest/cohort_registry.py`](00_ingest/cohort_registry.py) 纳入**两个配对数据集**
  （`GSE308103` snRNA + `GSE307534` 空间），产出样本/患者/分期冻结表（含 **23 例**配对关系）。
- **过门条件（全满足）**：
  1. 样本表字段齐全：`dataset, sample_id, patient_id, stage, modality`；
  2. **分期与 GEO 真值一致**（token 严格映射，未知即 raise）；**空转 LNM 不得伪造**（暂缺如实标注）；
  3. **23 例配对患者**（P3–P25）在**两模态**均有切片；
  4. 代码审查**无静默默认**；产出**冻结清单 + SHA-256 + 校验报告**。
- **不过门 → 停**（不得进 M1）。

### M-1 · 计划整改门 (Plan Remediation) —— **先于 M1，必做**

> **背景**：2026-09-12 对白皮书做**逐阶段可行性/严谨性审计**（数据身份 + 8 份方法学核查），
> 发现 **3 个 🔴（方法/数据不成立）+ 7 项 🟡**。**不过此门，不得进入 M1 及其后任何阶段。**

**A. 数据层**
- **A1** 空间图谱限定 **Normal→AAH→AIS→MIA→IAC**（全部 GSE307534，同库同平台）；**删除 GSE190811**（经 GEO 核实为**乳腺癌**，非 LUAD；原 GSM5732148 在该库不存在）。规则文档已改（[`docs/spatial_cohort_and_figure_prohibitions.md`](docs/spatial_cohort_and_figure_prohibitions.md)），白皮书待改。
- **A2** **LNM 不在范围**（空转无合法 LUAD 数据；单细胞层随三 scRNA 队列一并移出）。日后若获得真实 LUAD LNM 空转数据再补。
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
- **M3**：弃用「modality 当 batch」（sc↔sn 是 *system* 效应，且 modality 与 dataset 共线）→ **scVI(`batch=sample_id`)**。
  **⚠️ 2026-09-15 更正**：原文 `batch=dataset` 在**单数据集**（GSE308103）下该键为**常量、无意义**；75 个样本文库才是唯一真实技术批次轴。
  **且 `patient_id` 绝不可作 batch** —— 分期嵌套于患者内，校正它会抹掉 M4 门判据 ④ 所需的患者内配对对比。scIB 用**完整 panel**。
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

> 🔴 **2026-09-16 口径反转（用户决策「改用论文 QC 重建」）**。本节原写"**严禁**照搬整细胞 scRNA 阈值
> （实测 `nCount≥1000` 会砍掉约 30% 的核）"。该结论**在数量上仍然成立**（实测确砍掉约 30%），
> 但**据此拒绝论文阈值是错的**：论文的阈值来自**这批同样的样本、同一平台（10x Fixed RNA/FFPE）**，
> 是**原生口径**；本项目的逐样本 MAD 才是无出处的自选值。故**主口径改为论文固定阈值**，
> 原 MAD 口径**降级为敏感性臂**（产物保留：`results/02_expression/gse308103_counts.h5ad`，648,945 核）。

- **做**：显式 QC；**scDblFinder 逐样本**（铁律 R2）。数据为 **GSE308103（snRNA / 细胞核）**。
- **主口径（论文，逐字）**：剔除 `nFeature<500` 或 `nCount<1000` 或 `pct_mt>20%`；基因保留「在 ≥3 个细胞中检出」。
- **过门**：报告 QC 前后细胞数、双体率、参数（标出处，见 PARAMETERS）；**无启发式替代**。
- **★ 实测（2026-09-16）**：新分析集 = **413,697 核 × 18,069 基因**（旧 M1 掩膜 648,945；论文公布 401,635，**+12,062 = +3.00%**）。
  论文门在 M1 单细胞集内剔除 **235,248** 核：`nFeature<500` 剔除 **68,694**、另 **166,554** 仅因 `nCount<1000`、
  **`pct_mt>20%` 剔除 0 个**（M1 的自适应 MAD 已先清掉）⇒ 论文的线粒体门在本数据上**形同虚设**。
  边界：`pct_mt` 恰为 20.000 的核 **0** 个 ⇒ `<=20` 与 `<20` 两读法无差异。
- **⚠️ 登记偏差**：双体口径 **scDblFinder**（R2）vs 论文 **Scrublet**。R2 优先于"对齐论文"，已有 Scrublet 敏感性臂。
- **⚠️ 不得**因"要对齐论文"而把 `pct_mt` 当分数用——本项目的 `pct_mt` 是**百分数**（0–100），
  见 `01_qc/00_metrics_gse308103.R`：`pct <- 100 * colSums(m[mt,]) / pmax(cs,1)`。
- **AAH 脆弱性**：已在**新细胞集上重做实证检验**（`results/01_qc/stage_fragility_report_paperqc.json`）：
  中位 `nFeature` 比 Normal **1.081×**，AAH 留存 65.48% vs Normal 62.05%（**+3.43 pp**）⇒ **计数层面不脆弱**。
  故**不得**给 AAH 单独放宽 QC（那是无出处的静默默认）；脆弱性改在 **GP5 的聚类吸收层面**设可证伪护栏（指标4）。

### M2 · 上皮亚群恶性精判 (CNV) —— **执行位置：GP2，在 GP8a 上皮亚聚类之后**

> **⚠️ 2026-09-16 定位两次收窄**：① 从"降维聚类之前的**主力证真门**"改为"注释之后的上皮细化"；
> ② 再收窄为**"给上皮亚群贴恶性标签的精细判别器"**，执行位置由 GP6 之后**后移到 GP8a 之后**。
> 理由见 §2.0：CNV 在本项目的**唯一用途**是判**亚群**恶性，而亚群在 GP8a 之前不存在。

- **做**：**`CopyKAT` 为主力**（逐样本，**禁止合池**），输入 = GP6 注释为**上皮谱系**、并经 **GP8a** 划出亚群的细胞；
  **`infercnv` 仅作 5–10k 细胞/样本的可选交叉验证**（官方已停维护）。
- **产出粒度 = 亚群**：逐样本得 CNV± 后，按 **GP8a 的上皮亚群**汇总（每亚群 `frac_cnv_pos`、跨样本一致性），
  据此给**每个上皮亚群**贴 `malignant` / `non-malignant` / `not_testable`。
- **过门**：上皮亚群的恶性标签**三态齐全且样本级可溯**（不出现"全样本 0 非整倍体"被当作正常证据）。
  **不过门不得产出"恶性克隆"。** 措辞上限 = **亚群级（非整倍体 vs 否）**，**禁止**克隆/亚克隆/进化树表述。
- **下游使用**：M4 判据⑤、以及一切"恶性上皮"过滤，**一律以 GP2 的亚群标签为准**，不得用泛上皮 argmax 替代。

- **★★ 2026-09-16 实测发现：`CopyKAT` 存在两个必须用守卫堵住的静默失效路径**
  （源码 `tools/copykat/R/copykat.R`，copykat 1.2.5）：
  1. **`norm.cell.names=""` ⇒ 走自猜基线**（`copykat.R:166-170` `baseline.norm.cl()`）。低置信度时
     `WNS="unclassified.prediction"` 并**静默回退** `baseline.GMM(max.normal=5, mu.cut=0.05, Nfraq.cut=0.99)`，
     该函数**重写 `WNS`** ⇒ 日志只剩 `"low confidence in classification"`。
     **守卫**：必须传 `norm.cell.names=`（免疫/基质参考，`LUAD_NORM_REF=nonepi`），
     并断言日志出现 **`baseline is from known input`**（`copykat.R:136-145` 的可证伪标记）。
  2. **`:456` 全二倍体逃生门**：`if (cor(conses.diploid, conses.aneuploid) >= 0.6) com.preN[] <- "diploid"`
     —— **没有"不确定"类别**。纯上皮输入最易撞上此门。
     **实测**：`P11_LUAD` 纯上皮输入返回 **0/160 = 0.0%** 非整倍体 ⇒ **该臂是盲的**，
     其上任何"0 非整倍体"读数**不携带信息**，不得作为正常样本的证据。
     **守卫**：每个肿瘤样本必须**同时**跑"全细胞"对照臂；若对照臂也为 0，该样本 CNV 结论标 `not_testable`。
  - **副作用提醒**：`copykat.R:18-20` 只在**不合格细胞数 > 1** 时才过滤（恰好 1 个时不过滤）——
    任何复刻 copykat 细胞过滤的链条必须逐字照搬此语义，否则 `n_judged` 守卫会误炸。

- **注**：本数据为 **snRNA（GSE308103）** → CNV 需放宽参数并谨慎解读（核、低 UMI）；
  参考用同样本免疫/基质细胞（即上条守卫 1 的 `norm.cell.names`）。
- **退化样本预注册处置**（本项目约定，无文献阈值）：`n_retained<200` → `not_testable`，**绝不**填补为 "normal"；
  `n_cnv_pos==0` → 作为**真实结果**上报，**不得**放宽 `KS.cut` 重跑（法则 3.2 禁事后调参）；
  Normal/AAH 中 `frac_cnv_pos>0.95` 或 IAC 中 `<0.01` → `implausible`。
  **弃跑判据**：`not_testable` + `implausible` 合计 > 30% → **停下上报**（R2 对该数据集不可满足），
  **禁用**泛上皮 argmax 替代。
- **亚群级传导规则（2026-09-16 新增，跑之前固定）**：样本标 `not_testable` ⇒ 该样本贡献的**所有上皮亚群细胞**
  一并标 `not_testable`（**不得**因"其他样本正常"就补一个标签）；某亚群在**全部**样本上 `not_testable`
  ⇒ 该亚群恶性状态 = **`unknown`**，**不得**进入任何"恶性 vs 非恶性"对比。报告须给出
  `亚群 × 样本` 的 `not_testable` 覆盖表。**"unknown 亚群"占比 > 30% ⇒ 停在 GP2 上报。**
- **必须上报**：`annotateGenes.hg20()` 内部位置表未覆盖而被静默丢弃的基因清单（法则 0 静默默认陷阱）。

### M3-A · 传统/主流分支（执行位置：GP4a–GP8a → GP2 → GP8b）

> **⚠️ 2026-09-16 拆分说明**：原文把两件不同层面的事塞在同一个 M3 门里，现拆开：
> **M3-A 内部的两臂**（GP4c）= **Harmony on `sample_id`（主口径，论文配方）** vs **不校正（敏感性臂，非论文配方）**；
> **M3-B 的 SCMG**（GP9）= **另一条分支**。**三者不是同一层对照，报告时不得并表。**
>
> 🔴 **2026-09-16 主/对照反转**：原为「不校正＝**主**，Harmony＝对照」。论文把 Harmony 写进**主配方**，
> 故主口径改为 Harmony；不校正臂**保留**但降级为稳健性对照（`sample_id`↔`stage` 共线风险仍在）。

- **做**（按 §2.0 的顺序）：GP4a 预处理 → GP4b 降维 → GP4c 批次两臂 → GP5 分辨率 → GP6 注释
  → **GP8a 上皮亚聚类** → **GP2 CNV 定亚群恶性**（M2 门，见上）→ **GP8b 其余 5 谱系亚聚类**。
- **传统分支的"权威参数和顺序"以 [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md) M3-A.0–A.5 为准**，
  该表为执行前登记（法则 3.1），执行时**不得偏离而不登记**。
- **过门**（GP4a/b/c）：`SCTransform(vst.flavor="v2")` / HVG `3000` / `npcs=50` / `Harmony dims=1:50` /
  `k.param=20` 全部**带哈希登记**；**参数来源逐条可溯**（论文逐字引文 or Seurat 源码默认）。
  ⚠️ **原 `{20,30,50}` ARI 护栏随旧 PCA-30 口径一并作废**，不适用于论文的固定 50 PC。
- **两臂上报规则**：逐细胞 `ARI(Harmony, 不校正)` < 0.7 ⇒ **不得择一上报**，须同时呈现两臂结果并升级裁定。
- **执行脚本**：[`04_integration/10_seurat_traditional.R`](04_integration/10_seurat_traditional.R)
  （`smoke --n_cells N` / `full` 两种模式；已过 3,000 与 60,000 细胞烟雾测试）。
- **禁用**：`modality 当 batch`（sc↔sn 为 system 效应且与 dataset 共线）；**`patient_id` 作 batch**
  （分期嵌套于患者内，校正会抹掉 M4 判据 ④ 所需的患者内配对）。

### M3-B · SCMG 对照分支（执行位置：GP9，在 M3-A 之后）

- **做**：**zero-shot 跨数据集 scRNA 整合 + 全局流形 + 状态刻画**
  （**不掺传统算法**——铁律 5；**不输出逆转/因果**——铁律 6）。
- **过门**：以 **scIB 完整口径**与 M3-A 对照——批去除（kBET + iLISI + graph-connectivity + PCR）
  **与** 生物保守（cLISI/ARI/NMI/ASW）**并报**；batch panel 与 bio panel **分开报**，附逐指标原始值。
- **判据重写（2026-09-15）**：不再是「最大化 iLISI」，而是 **batch 指标可接受 _且_ 分期可分性被保留**
  （证明**有分辨力**，非"全糊一块"）。
- **门规格变更**：单数据集**无 ground truth** ⇒ 生物保守的 ARI/NMI 只能算**未校正 vs 已校正簇标签之间**的，
  **须在 GP3 显式批准**后方可计算。
- **裁决**：Arm B 相对 Arm A 的 **ARI < 0.7** ⇒ Arm A 为主，如实报方法学局限。
- **禁用**：仅凭 iLISI↑ 判为整合成功（可被过度整合刷高）。

### M4 · 跨模态 AAH 门（**五判据缺一不可**）
- ① 各阶段**可分辨**；② 重叠阶段 sn↔sc **一致**；③ 平台偏移 **δ(stage) 稳定**；④ **sn 内部配对**（AAH vs 同患者 Normal/AIS）同向；⑤ AAH 身份 **CNV/标记**证真。
- **不过门 → 不下"AAH 结论"**（只能标为假说）。
- 数据：`GSE308103`(snRNA，主参考)；如日后获批 `HRA001130`(全细胞 scRNA，受控) 可作 AAH 跨模态验证。

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
- **标准/主流分支**：scVI（`batch=sample_id`）/scANVI + scArches、Harmony、Seurat、RCTD、SpaGCN、Squidpy、DESeq2、CellRank、SCENIC+…
- **纯 SCMG 分支**：**zero-shot 跨数据集 scRNA 整合 → 全局流形 → 细胞状态刻画**（**不掺传统算法**）。
  ~~条件扩散轨迹 → CausalGenePredictor 因果~~ —— **已删除**（能力不存在，见铁律 6）。
- 两分支**对照**（scIB 完整口径）；工具源码见 [`tools/`](tools/)。
- 靶点/调控因子候选**只由标准分支产出**，且措辞为「候选」。

> **⚠️ 2026-09-16 层级澄清**：本文档里"臂(arm)"和"分支(branch)"是**两个不同层级**，不得混用：
> - **分支** = 顶层二分（**标准/主流分支** vs **纯 SCMG 分支**），两个分支**各自独立跑完整流程**后**互相对照**。
> - **臂** = **标准分支内部**的批次处理对照（**Arm A 不校正** vs **Arm B Harmony on `sample_id`**），
>   共用同一套预处理之外的下游，只在 GP4c 分岔。
> 因此 **Arm A/B 的对照结论只属于 M3-A**，**不得**与 M3-B 的 SCMG 对照并成一张表。

---

## 5. 数据来源与获取（**范围已收窄：仅两个配对数据集**）

> **2026-09-15 修正**：配对患者数 **9 → 23**（P3–P25）。**旧值 9 是"仅下载 19/56 张空间切片"时的产物**，
> 切片下载齐后按两张 GEO 权威表求交集实为 **23 例**。见 `00_ingest/cohort_registry.py` 的 `PAIRED_PATIENTS_MIN`。
> **2026-09-12 收窄**：仅用 **`GSE308103`(snRNA) + `GSE307534`(Visium 空间)** —— 同一研究、**模态匹配（FFPE↔FFPE）**。
> 三个 scRNA 队列（GSE131907/189357/148071）**移出范围**，存档 `/home/eto/luad_invasion/luad_v2_out_of_scope/`。

- **单细胞/参考**：`GSE308103`（snRNA，75 样本 / 798,100 核实测）——**唯一含 AAH** 的单细胞资源；
- **空间**：`GSE307534`（Visium CytAssist FFPE；GEO 56 样本 / 25 患者，**本地 56 张切片齐**，覆盖全部 **23 例**配对患者；仅 P1/P2 无 snRNA 不入配对）；
- **LNM**：空转暂缺（~~GSE190811~~ 经核实为**乳腺癌**，已废）；
- （范围外，已存档）`HRA001130`（sc，受控）预留接口 → `/home/eto/luad_invasion/luad_v2_out_of_scope/`；
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

> **2026-09-16 重构**：本表原按 M 编号排序，与 §2.0 的实际执行序不一致（M2 被排在 M3 前）。
> 现**改按执行序排列**，并补入此前完全缺失的 GP1–GP9 行。

| 序 | 检查点 | 模块 | 状态 | 过门 |
| :---: | :--- | :--- | :---: | :---: |
| — | M0 输入冻结 | 输入 | ⚠️ 未过门 | ☐ |
| — | M-1 计划整改 | 计划 | 🔶 进行中 | ☐ |
| 1 | M1 QC/双体（🔶 降级为敏感性臂） | 质控 | ✅ 100% | ✅ |
| 2 | **GP0** 表达对象重建 | 表达 | 🔶 **2026-09-16 按论文 QC 重做**；09-17 审计补校验 **13/0** | ☐ |
| 3 | **GP4a** 预处理（SCTransform v2 → HVG 3000） | M3-A | ✅ **2026-09-17 全量完成** | ✅ |
| 4 | **GP4b** 降维（PCA 50，无 ARI 护栏） | M3-A | ✅ **2026-09-17 全量完成** | ✅ |
| 5 | **GP4c** 批次两臂（**Harmony 主** / 不校正对照） | M3-A | ✅ **2026-09-17 裁定：采用论文 Harmony 主口径** | ✅ |
| 6 | **GP5** 分辨率选择（区间 0.5–0.8 × 5 种子 + AAH 护栏） | M3-A | ✅ **2026-09-17 签字 `r*=0.6`**；指标4 认无效 | ✅ |
| 7 | **GP6** 双标准注释 + κ（**只定"哪些是上皮"**） | M3-A | ⬜ 0% | ☐ |
| 8 | **GP8a** **上皮**亚聚类 → 上皮亚群 | M3-A | ⬜ 0% | ☐ |
| 9 | **GP2** CNV **亚群**恶性精判 | **M2** | ⬜ 等 GP8a | ☐ |
| 10 | **GP8b** 其余 5 谱系亚聚类（不被 CNV 阻塞） | M3-A | ⬜ 0% | ☐ |
| 11 | **GP3** 指标后端冻结 → **GP7** scIB 双 panel | M3-A | ⬜ 0% | ☐ |
| 12 | **GP9** SCMG 对照臂 | **M3-B** | ⬜ 0% | ☐ |
| 13 | M4 跨模态 AAH | — | ⬜ 0% | ☐ |
| 14 | M5 空间解卷积 | — | ⬜ 0% | ☐ |
| 15 | M6 生态位 | — | ⬜ 0% | ☐ |
| 16 | M7 靶点(传统轨道) | — | ⬜ 0% | ☐ |
| 17 | M8 对接 | — | ⬜ 0% | ☐ |

> **GP1（CopyKAT 冒烟）已并入 GP2**，不再单列：它是 GP2 的**执行细节**（先 3 点外推再排期），不是独立门。
> 参考：`results/03_cnv/smoke/P19_LUAD/` —— 该冒烟**未能收敛**（13h CPU 停在 step 7），
> 是促成"CNV 收窄 + 输入改上皮亚群"的实测依据之一。

> 每过一个门 → 更新本表 + 记一笔"过门证据"（产物路径 + 哈希）。
>
> ⚠️ **GP4a/4b/4c 合并为一个执行单元**：`04_integration/10_seurat_traditional.R` 一次跑完
> （SCTransform → PCA → Harmony → 两臂 kNN → 聚类 → UMAP），但**过门证据按三步分别记**，
> 因为三者是三个可独立否决的判据，不得因"一个脚本跑通"就合并签字。

### 过门证据 (Gate Evidence)

**GP4a / GP4b / GP4c / GP5 — ✅ 全部签字 (2026-09-17)**（GP4a、GP4b 全量跑通于 00:43；GP4c、GP5 经用户裁定签字）
- 脚本：`04_integration/10_seurat_traditional.R`（`full` 模式）· 日志 `logs/seurat_full.log`
- 输入：`results/02_expression/gse308103_counts_paperqc.h5ad` sha256 `a276cd1a…`（**已现场复核与 manifest 一致**）
  = **413,697 核 × 18,069 基因**，nnz 655,222,534
- 配方：`SCTransform(vst.flavor="v2", variable.features.n=3000, rv.th=1.3, ncells=5000, seed=1448145)`
  → `RunPCA(npcs=50, seed=42)` → `Harmony(group.by.vars="sample_id", dims.use=1:50)`
  → `FindNeighbors(k.param=20, prune.SNN=1/15)` → Louvain × `{0.5,0.6,0.7,0.8}` × 种子 `{0,1,2,3,4}` + 不校正臂
- 耗时 **14,716.6 s = 4 h 05 m**；**全程峰值 VmHWM 252.5 GB / 256 GB（98.5%），无 swap，可用内存一度只剩 9 GB**
- 版本：R 4.2.2 · Seurat 4.3.0 · sctransform 0.3.5 · harmony 2.0.5（论文为 1.2.0，**不声称数值等价**）· Matrix 1.5.3
- 产物哈希：`clusters.csv.gz` `4e6e6977…` · `umap.csv.gz` `c33995dd…` · `resolution_metrics.csv` `58c35608…` ·
  `run_manifest.json` `a57c977c…` · `timings.csv` `7d8724f2…` · `embeddings/harmony_f32.bin` `93f2cf96…` · `embeddings/pca_f32.bin` `f529b2f2…`
- 报告：[`results/04_integration/seurat_trad/full/GP4_report.md`](results/04_integration/seurat_trad/full/GP4_report.md)

- ✅ **GP4c 已签字 —— 用户裁定 (2026-09-17)：采用源论文的 Harmony 主口径。**
  升级条款触发经过如实报出，用户裁定主口径＝`RunHarmony(group.by.vars="sample_id", dims.use=1:50)`（与 §M3-A.2 预注册主口径一致）；
  **不校正臂保留为敏感性/局限陈述，不得当作主结果**。
  ⚠️ **这是对"两臂分歧"的显式裁定，不是分歧消失 —— 下列限定对全部下游结果具约束力：**
  逐细胞 `ARI(Harmony, 不校正)` = r0.6 **0.6575（5/5 种子 <0.7，稳健）**、r0.7 **0.6960（1/5，边缘）**；
  kNN(k=15) 从嵌入预测分期：**校正前 0.5388（基线 0.4546，+8.42 pp）→ 校正后 0.4866（+3.19 pp）⇒ 抹掉约 60% 分期可恢复信号**；
  `sample_id`↔`stage` 一对一嵌套 ⇒ **被抹掉的部分数学上不可归因**。
  ⇒ 下游任何分期相关结论**不得声称已排除批次混淆**；M4 须把这层不可归因性写入解释边界。
  版本偏差仍在：论文 Harmony 1.2.0 vs 本机 2.0.5 ⇒ **不声称数值等价**。
- ✅ **GP5 已签字，`r* = 0.6`**（合法候选仅 r=0.6 / r=0.7）：
  指标1（跨种子 ARI ≥ 0.90 硬约束）实测 0.8963–0.9230，**剔除 r=0.5（0.8963）与 r=0.8（0.8999）** ⇒ 该护栏在全量规模**真实承重**（3,000 细胞时恒 = 1.0 空转，60,000 细胞时 0.9346–0.9834 —— **不得用小规模推断本规模**）。
  按预注册破平规则（分差 <0.01 取较低）0.9123−0.9041=0.0082 ⇒ **取 r=0.6**。
  ⚠️ **脚本原第 305 行为裸 `which.max`，行内注释却声称已破平 —— 注释与代码不符**，故全量 run 报出 0.7。
  已修脚本为显式破平实现，并同步更正 `resolution_metrics.csv` 的 `is_rstar` 列与 manifest 的 `rstar_candidate`/`metrics_table[*].is_rstar`，
  manifest 内新增 **`rstar_correction`** 块留痕；**原始指标一个字节未动**。**用户接受该更正经过程序。**
- ✅ **指标4（AAH 吸收护栏）—— 用户签字：承认其在全量规模上无效。** 4 个分辨率下**没有任何簇满足"Normal 占比最高且 >50%"**，
  分母为空 ⇒ `absorption_rate` 按定义恒为 0，**没有为 `r*` 提供任何信息**。
  ⇒ `r*` **实际由指标 1（跨种子 ARI）与指标 2（跨分辨率稳定性）承担**；指标 3 待 GP6 后补算；指标 4 不计入。
  护栏**保留登记但标记"全量无效"**，不删除、不重设计；将来若在别的数据/粒度重启，须**重新预注册**。
- ⚠️ **指标3（谱系覆盖）本轮未计算** —— 需 GP6 注释之后方可算，manifest 已标 `metric3_absent`。
  故表中 `eligible` 只反映 1/4 两条约束，**不是完整三约束判定**。
- 诚实记录：`input.n_genes` 18,069 → 变换矩阵 18,047，差 **22 个**基因，由 sctransform v2 内部过滤产生（日志可查，**非静默**）。

**Step 0 · GP0 表达对象重建（论文 QC 口径）— 🔶 已重做 (2026-09-16)，校验由 2026-09-17 审计补做**
- 脚本：`02_expression/04_rebuild_expression_paperqc.py`；来源 = 75 个稠密文本计数矩阵 → 稀疏 AnnData
- 掩码：**M1 `qc_pass & singlet` ∩ 论文 QC 门**（`nFeature≥500 & nCount≥1000 & pct_mt≤20`），逆中 `gene_min_cells=3` 在重建矩阵时施加
- 结果：**(413,697 核 × 18,069 基因)**，nnz **655,222,534**；分期 IAC 188,087 / Normal 94,506 / AIS 81,555 / AAH 35,283 / MIA 14,266（**23 患者 / 75 样本**）
- 产物：`results/02_expression/gse308103_counts_paperqc.h5ad` `a276cd1a…`（3.13 GB）· `results/01_qc/gse308103_analysis_mask_paperqc.csv.gz` `ebc74c1e…` · `rebuild_paperqc_manifest.json`
- ⚠️ **本次审计发现的缺口（已补）**：GP0-redo **当初没跑门校验、也没出报告** —— 因为 `02_verify_expression_build.py` 的 `EXPECT_N_OBS/N_VARS` **硬编码在旧对象上**（648945 / 18082），从未适配重做版。**该缺口正是"§6 曾把旧 648,945 构建记为 ✅ PASS"的根因。**
  2026-09-17 审计对 paperqc 对象**重跑全部硬检查**，独立复算 + 现场哈希核对，结果：
  V1 形状 ✅ · V2 逐样本偏差 **0** ✅ · V4 条码可解析 & `sample_id` 全在权威表 ✅ · V5 `stage == resolve_stage(token)` 无静默默认 ✅ · V6 全元素非负/整数、全零细胞 **0** ✅
  V3 修正形：`nnz == Σ nFeature − 11 = 655,222,545 − 11 = 655,222,534` ✅ —— **差 11 是构造性的**，即那 13 个 `<3 细胞` 被丢基因内残留的 11 个非零计数（旧脚本 V3 的"精确相等"形对重做版**不成立**，需带此项）
  **掩膜 ↔ h5ad ↔ manifest 三方独立对账**：条码集合完全相等（`S1 == S2`，双向差 0）、`stage/patient_id/sample_id/stage_token` 逐行一致、无 NaN、条码唯一、`stage_counts` 相等
  **下游 seurat_io 导出**（R 侧实际读入）：`cell_names`/`gene_names` 行数与**顺序**均等于 h5ad、`cell_meta` 条码集合相等、6 个产物 sha256 **逐个吻合清单**
- **未复核项（如实标缺）**：重做版**没有**旧的 A5 独立重抽取（3 样本逐元素重算）与 A4b 全元素整数性专项脚本留痕；本次 V6 已全量覆盖整数性，A5 类重抽取**仍未对 paperqc 对象做**。
- 🔎 **旧 648,945 对象保留为敏感性臂，未删除**：`gse308103_counts.h5ad` `f9dbe382…` · `gse308103_analysis_mask.csv.gz` `2f8bb0f6…` · `GP0_report.md`（该报告只描述旧对象）
- 旧对象历史（当时确已过门，**现已被论文口径取代**）：V1–V6 12/12 + 对抗审计 22/22 全绿，含 A5 独立重抽取、A2b 逐细胞行和 == nCount（偏差 0）、A6b/c R1 闭合 75/75 & 23/23；残余盲区 (行和,非零数) 重复 252,677/648,945
- 诚实记录：旧对象**首次运行崩于** `KeyError ['cell_barcode']` → 已定位修复；**崩溃未污染数据**。日志 `logs/Step0_build_expression.run1_FAILED.stdout`

**M1 QC / 双体 — ✅ PASS (2026-09-12)｜🔶 现状：按 §2.0 行 1 已**降级为敏感性臂**
- 数据集：`GSE308103`（snRNA）**75 样本 / 798,100 核**
- 脚本：`01_qc/00_metrics_gse308103.R` → `01_qc/01_qc_doublets_gse308103.R` → `02_annotate_doublet_qc.R` → `03_sensitivity_nmads.R` → `06_sensitivity_doublet_rate.R` → `08_validate_doublet_calls.R`
- 结果：pre **798,100** → pass **767,839（96.21%）**；双体 **118,894（15.48% of pass）**
- 阈值：`nCount/nFeature` 逐样本 **MAD 离群**（nmads=3, log1p）+ `pct_mt<5`（**核数据据实定**，非照搬 scRNA）
- ⚠️ **口径叠加（勿误读）**：分析掩膜 = **M1 `qc_pass & singlet` ∩ 论文 QC 门**（`nFeature≥500 & nCount≥1000 & pct_mt≤20`）。
  即**论文的绝对阈值照样施加了**，"非照搬 scRNA" 指的是 **M1 自身**的判定口径，**不是**说绝对阈值被否决。链：798,100 → 767,839(M1) → 648,945(∩单细胞) → 论文门在全集上 555,480 → **交集 413,697**
- 敏感性：nmads 3 vs 5 = +1.97 pp（**不敏感**）；双体率 vs 固定 top-10% 重合约 **62.3%**（**较敏感** → 须做下游"剔/不剔"敏感性）
- 正向验证：**A** 计数特征（75/75 样本双体 nCount 比中位 **2.39**；双体率 vs 细胞数 **r=0.921**）；**B** 跨谱系共表达（EPCAM+PTPRC+ 在 doublet 中为 singlet 的 **7.2×**，72/75 样本一致）
- 交叉验证：**本环境不可行**（scrublet 不适配稀疏核；DoubletFinder 需 Seurat 2/3 或 5）—— 已如实记录
- **已解决异常**：`P7_LUAD` 曾判 0 双体 → 根因为 **xgb 分类器塌缩**（非生物学）→ 改 `score="weighted"` 得 11.69%；主脚本已加**自动 fallback**
- 报告：[`results/01_qc/M1_validation_report.md`](results/01_qc/M1_validation_report.md)（v3）
- 产物哈希：`gse308103_per_cell_qc.csv.gz` `44c890bb…` · `gse308103_qc_per_sample.csv` `2d0bd7df…`
- ⚠️ v1 的 per-cell 表存在 `cell_barcode` 失效缺陷（fread autostart 跳过条码行）→ 已修复并全量重跑，v1/v2 报告作废

**M0 输入冻结（配对）— ⚠️ 未过门 (2026-09-15 重做)**
- 脚本：`00_ingest/01_freeze_paired.py`；registry `00_ingest/cohort_registry.py`
- 结果：**131 样本 / 25 患者 / 23 例双模态配对**（P3–P25）；唯一键 = `sample_key`（`dataset:sample_id`）
- 产物：`results/00_ingest/paired_{samples,patients,source_files}.csv` + `M0_paired_validation_report.md` + `paired_manifest.json`
- 校验：C1–C7 绿；**C8 红 1 项** → `paired_manifest.json` 的 `gate_pass=false`
- **未过门原因（唯一）**：`GSE307534/GSM9226176` 磁盘上的 tar **截断**（56,272,384 B，应为 90,677,930 B；
  `gzip -t` 报 `unexpected end of file`），缺 `spatial/scalefactors_json.json` 与 `spatial/tissue_positions.csv`。
  已实测重下载可得**完整 87 MB** tar（清单含全部必需文件，且**只有一个切片根** `P4_AAH2`）。
- **不阻塞 M2/M3-A**：该缺口在**空间**数据集（GSE307534），而 M2/M3-A 只跑 **snRNA**（GSE308103）。补下载属 M5 前置，
  且源目录 `/home/eto/luad_invasion` 为**只读**，需单独授权后另做。
- ⚠️ 旧 M0（三 scRNA 队列）**已随范围收窄作废**，其脚本/产物移至 `/home/eto/luad_invasion/luad_v2_out_of_scope/`。

> **历史记录（旧范围，已作废）**：曾冻结 3 个 scRNA 队列（GSE131907 / GSE189357 / GSE148071，420,766 细胞 / 109 样本 / 95 患者），
> 并完成 GEO 交叉核验（分期与 origin）。该产物已归档，**不再作为本项目依据**。

