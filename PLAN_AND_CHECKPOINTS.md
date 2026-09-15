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
- **做**：显式 QC；**scDblFinder 逐样本**。数据为 **GSE308103（snRNA / 细胞核）**。
- **过门**：报告 QC 前后细胞数、双体率、参数（标出处，见 PARAMETERS）；**无启发式替代**。
- **★ 必须**：QC 阈值**按核数据据实确定**（`nCount/nFeature` 逐样本 **MAD 离群** `nmads=3`, log1p, 双尾；
  `pct_mt < 5`）——**严禁照搬整细胞 scRNA 阈值**（实测 `nCount≥1000` 会砍掉约 30% 的核）。

### M2 · 恶性证真门 (CNV)
- **做**：**`CopyKAT` 为主力**（逐样本，**禁止合池**）；**`infercnv` 仅作 5–10k 细胞/样本的可选交叉验证**（官方已停维护）。
- **过门**：报告 CNV 阳性细胞数 + 与经典标记一致性；**恶性细胞的下游使用以此为准**。**不过门不得产出"恶性克隆"。**
- **注**：本数据为 **snRNA（GSE308103）** → CNV 需放宽参数并谨慎解读（核、低 UMI）；参考用同样本免疫/基质细胞。

### M3 · 整合门（双分支）
- **做**：**A 标准分支**：**scVI（`batch=sample_id`）**（原文 `batch=dataset` 见 §C 更正）。
  **双臂**：**Arm A（不校正，主）** vs **Arm B（Harmony/scVI on `sample_id`，对照）**。
  **B SCMG 分支**：**zero-shot 跨数据集 scRNA 整合 + 全局流形 + 状态刻画**（**不掺传统算法**；**不输出逆转/因果**）。
- **过门**：以 **scIB 完整口径**给数——批去除（kBET + iLISI + graph-connectivity + PCR）**与** 生物保守（cLISI/ARI/NMI/ASW）**并报**；
  **判据重写（2026-09-15）**：不再是「最大化 iLISI」，而是 **batch 指标可接受 _且_ 分期可分性被保留**（证明**有分辨力**，非"全糊一块"）。
  单数据集**无 ground truth** ⇒ 生物保守的 ARI/NMI 只能算**未校正 vs 已校正簇标签之间**的，**须在 GP3 显式批准此门规格变更**。
- **禁用**：`modality 当 batch`（sc↔sn 为 system 效应且与 dataset 共线）；**`patient_id` 作 batch**；仅凭 iLISI↑ 不得判为整合成功（可被过度整合刷高）。
- **裁决**：Arm B 相对 Arm A 的 **ARI < 0.7** ⇒ Arm A 为主，如实报方法学局限。

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
| 里程碑 | 状态 | 过门 |
| :--- | :---: | :---: |
| M0 输入冻结 | ⚠️ 未过门 | ☐ |
| M-1 计划整改 | 🔶 进行中 | ☐ |
| M1 QC/双体 | ✅ 100% | ✅ |
| Step 0 表达对象重建 (GP0) | ✅ 100% | ✅ |
| M2 CNV 证真 | ⬜ 0% | ☐ |
| M3 整合(双分支) | ⬜ 0% | ☐ |
| M4 跨模态 AAH | ⬜ 0% | ☐ |
| M5 空间解卷积 | ⬜ 0% | ☐ |
| M6 生态位 | ⬜ 0% | ☐ |
| M7 靶点(传统轨道) | ⬜ 0% | ☐ |
| M8 对接 | ⬜ 0% | ☐ |

> 每过一个门 → 更新本表 + 记一笔"过门证据"（产物路径 + 哈希）。

### 过门证据 (Gate Evidence)

**Step 0 · GP0 表达对象重建 — ✅ PASS (2026-09-15 23:27)**
- 脚本：`02_expression/01_build_expression_gse308103.py`（构建）· `02_expression/02_verify_expression_build.py`（**对写出的 h5ad 复核**，不看构建器内存态）
- 来源：75 个稠密文本计数矩阵 → 稀疏 AnnData；掩码 = M1 的 `qc_pass & doublet_class=='singlet'`
- 结果：**(648,945 细胞 × 18,082 基因)**，nnz **791,571,728**，稀疏度 93.25%；分期 IAC 294,684 / Normal 152,302 / AIS 123,988 / AAH 53,887 / MIA 24,084（23 患者 / 75 样本）
- 校验：V1–V6 **12/12 全绿**；V2 逐样本偏差 **0**；V3 nnz **精确等于** Σ nFeature；V3b **逐细胞** nnz 偏差 **0**；V5 分期无静默默认
- **对抗性审计**（`02_expression/03_audit_expression_build.py`）：**22/22 全绿**。含 **A5 独立重抽取**（不同代码路径从原始文本重算 3 样本 × 逐元素比对）、**A2b 逐细胞行和 == nCount**（偏差 0）、**A4b 全 7.9 亿元素整数性**（V6b 只抽查 500 万）、**A6b/c R1 闭合**（token/patient/stage vs GEO 权威 75/75、23/23）
- 额外独立核对：mask CSV 与 h5ad obs **逐行**四项一致；75 文件基因向量 SHA-256 唯一值数 = 1；float32 精度安全（max 52,227 ≪ 2²⁴）
- 审计残余盲区（如实标缺）：(行和,非零数) 组合重复 252,677/648,945 → 同摘要细胞互换的盲区，由 A5/A7b 兜底；A5 仅覆盖 3/75 样本
- 产物：`results/02_expression/gse308103_counts.h5ad` `f9dbe382…`（1.91 GB）· `gse308103_analysis_mask.csv.gz` `2f8bb0f6…` · `per_sample/*.npz` ×75 · `build_manifest.json`
- 报告：[`results/02_expression/GP0_report.md`](results/02_expression/GP0_report.md)
- 诚实记录：**首次运行崩于** `KeyError ['cell_barcode']`（`ad.AnnData` 就地持有 obs，`obs_names.name=None` 污染了其后的 `reset_index()`）→ 已定位、修复、重跑；**崩溃未污染数据**（h5ad 当时已成功写出）。失败日志留档 `logs/Step0_build_expression.run1_FAILED.stdout`

**M1 QC / 双体 — ✅ PASS (2026-09-12)**
- 数据集：`GSE308103`（snRNA）**75 样本 / 798,100 核**
- 脚本：`01_qc/00_metrics_gse308103.R` → `01_qc/01_qc_doublets_gse308103.R` → `02_annotate_doublet_qc.R` → `03_sensitivity_nmads.R` → `06_sensitivity_doublet_rate.R` → `08_validate_doublet_calls.R`
- 结果：pre **798,100** → pass **767,839（96.21%）**；双体 **118,894（15.48% of pass）**
- 阈值：`nCount/nFeature` 逐样本 **MAD 离群**（nmads=3, log1p）+ `pct_mt<5`（**核数据据实定**，非照搬 scRNA）
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

