# 诊断性预注册 · inferCNV 在我们单核上皮上跑不跑得通

**登记时间**：2026-09-23（**计算开始之前**，见文件 mtime）
**性质**：🔬 **冒烟诊断（smoke test）**，不是交付物。**不改变任何已签口径**，不产出可引用的分期结论。
**授权**：用户 2026-09-23 「inferCNV这个试试」。
**用户当场拍板的两项**（2026-09-23）：
1. 规模 = **单患者 P4**
2. 参考 = **同患者 Normal 全部上皮**（**注意这是对已签「同亚型」口径的偏离**，见 §六偏差登记）

---

## 一、要回答的问题

`PLAN_AND_CHECKPOINTS.md:100` 把 infercnv 登记为「**降级为可选交叉验证**」，理由是三条环境障碍。
**本冒烟只问三件事，不问"谁恶性"**：

1. 三条障碍在实测下是否真的成立？（装机层面）
2. 逐患者跑能否在可接受的时间内跑通？（工程层面）
3. 它导出的**连续残差矩阵**，形状是否合理、有没有可用于判别的信号？（方法学层面）

**是否把 inferCNV 采纳为恶性口径，不在本冒烟范围内**——那要另开预注册、另行签字。

## 二、工具与资源（冻结）

| 项 | 值 |
|---|---|
| 工具 | `infercnv` R 包 **v1.23.0**（2023-12-01），**本地源码** `/home/eto/luad_v2/tools/infercnv/` |
| 安装位置 | **隔离库** `/home/eto/Rlibs/infercnv`（照 `/home/eto/Rlibs/SCEVAN` 先例；**主 R 库未动**） |
| JAGS | **自建**，前缀 `/home/eto/local/jags`（JAGS 4.3.2 源码；**未用 sudo，未动系统**） |
| rjags | 4.17（隔离库内，链接 `~/local/jags/lib/libjags.so.4`） |
| 运行时前提 | **必须** `export LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH`（否则 rjags 装载失败） |
| 基因位置来源 | **GENCODE v44 (hg38)**，`data/external/gencode/gencode.v44.annotation.gtf.gz`（EBI 镜像，48 MB，权威注释） |

**三条注册障碍的实测处置**（`PLAN_AND_CHECKPOINTS.md:100`）：

| 注册障碍 | 实测 |
|---|---|
| 官方 README 声明 *no longer supported* | **不成立**，本地源码完整可装，`library(infercnv)` 通过 |
| 需系统 JAGS | **已解决**，源码自建，无需 sudo |
| 本机 Bioconductor 不可达 | **绕开**，从 vendored 源码 `R CMD INSTALL`，未用 Bioconductor |

## 三、输入（冻结）

- 表达矩阵：`results/02_expression/gse308103_counts_paperqc.h5ad`（sha256 `a276cd1af69a4620ff9e3d64f4b93f33461e646537c5f2c9da1aff8c9150b9de`）
- 细胞子集：`results/05_annotation/epiA_subset_barcodes.txt`（sha256 `d12a115a87d0790d7fe8cf27b16170999d94daba51f3016fdb79c5aa3fdde55e`）
- 细胞：**P4 全部上皮核 = 19,748**，分 5 个样本：

  | 样本 | stage | n |
  |---|---|---|
  | `P4_Normal` | Normal | 5,824 |
  | `P4_Normal1` | Normal | 2,711 |
  | `P4_AAH` | AAH | 1,009 |
  | `P4_AAH1` | AAH | 2,702 |
  | `P4_LUAD` | IAC | 7,502 |

  —— 选 P4 因它是本项目上皮最多的患者（次为 P23 19,130）。**不是挑结果好的患者**（选择依据在计算前写定）。
  **勘误（2026-09-23，仍在计算前）**：本文件初稿把 P4 写成「Normal 5,824 + LUAD 7,502 = 13,326」，
  漏了 `P4_AAH` / `P4_AAH1` / `P4_Normal1` 三个样本。已更正为 19,748。P4 有**两个** Normal 样本。
- 分期标签：`results/04_integration/seurat_trad/epiA/clusters.csv.gz` 的 `stage` 列（R1 权威口径）

**为什么只能逐患者跑**：inferCNV 内部是**稠密矩阵**。133,384 核 × 18,069 基因 × 8 B ≈ **19 TB**，不可能。
逐患者最大约 2.9 GB，可行。这与已签的「锚定 = 同患者」口径方向一致。

## 四、🔴 圆环性：本冒烟最要紧的设计点

**参考细胞定义基线 ⇒ 参考细胞自己的残差必然 ≈ 1.0。**
所以在同一次跑里比较「Normal 组 vs 瘤组」是**构造出来的**，**不是证据**。

**处置（计算前写死）**：从参考组（`P4_Normal` + `P4_Normal1`，共 8,535 个正常上皮）
里**抽掉 15%（种子 0）**，单独标成一个组 `Normal_holdout`，使它们**不参与基线计算**。
评估只用 `Normal_holdout` vs `P4_LUAD`——这一条**不是循环的**。

**分组写死**（5 组）：

| 组名 | 细胞 | 角色 |
|---|---|---|
| `Normal` | 8,535 × 85% | **参考组**（`ref_group_names`） |
| `Normal_holdout` | 8,535 × 15%，种子 0 | 观测，**门槛 3 的负类** |
| `P4_LUAD` | 7,502 | 观测，**门槛 3 的正类** |
| `P4_AAH` | 1,009 | 观测，**描述性**（不参与判读） |
| `P4_AAH1` | 2,702 | 观测，**描述性**（不参与判读） |

（被抽掉的细胞仍是同一患者的正常上皮，故仍共享患者的样本级批次。这层混淆**无法在本设计内排除**，如实登记。）

## 五、参数（**计算前登记，事后不许改**，法则 3.1 / 3.2）

| 参数 | 取值 | 依据 |
|---|---|---|
| `cutoff` | **0.1** | inferCNV 官方 10x 教程标准值（默认 1 是 Smart-seq2 口径） |
| `window_length` | **101** | 官方默认 |
| `denoise` | **TRUE** | 官方默认；抗单核高噪声 |
| `HMM` | **FALSE** | **本冒烟不走 JAGS 通路**；要的是连续残差，不是它的 HMM 判决 |
| `analysis_mode` | **"subclusters"** | 官方默认 |
| `cluster_by_groups` | **TRUE** | 官方默认 |
| `write_expr_matrix` | **TRUE** | 导出残差矩阵 |
| `no_plot` | **FALSE** | 顺带出一张热图 |
| `num_threads` | **8** | 本机 20 核，留余量 |
| 留出比例 / 种子 | **0.15 / 0** | 见 §四 |
| 参考组 | `Normal`（P4 自己的 Normal 上皮，扣除留出集） | 用户 2026-09-23 拍板 |

## 六、偏差登记（**必须带进报告**）

1. **参考定义偏离已签口径**：已签锚定是「同患者 Normal **同亚型**」；本次按用户 2026-09-23
   拍板改用「同患者 Normal **全部上皮**」。理由：inferCNV 一次只吃**一个**参考组，
   且部分亚型的 Normal 细胞过少。**这是冒烟口径，不构成新口径，不得被后续引用。**
2. **单核未验证**：inferCNV 及 2026 基准均**未在 snRNA 上验证**过。
3. **工具定位**：inferCNV 是**可视化**工具。2026 基准（bioRxiv 10.64898/2026.04.12.718050）
   明确警告它「把细胞分成**亚克隆**，而非**逐细胞**判决」。故本次要的是**导出的残差矩阵**，不是它的图。
4. **单患者单次跑**：不得推及全体。
5. **环境 RNA / 内含子滞留**会抬高残差；本数据无环境校正（见 P24 那条已确认的环境 RNA 结论）。

## 七、判读（**看到结果前定死**）

| # | 门槛 | 通过判据 |
|---|---|---|
| 1 | **能跑通** | 无报错跑完；**基因保留率 ≥ 90%**（相对 18,069） |
| 2 | **矩阵合理** | 参考组细胞的残差**中位数落在 [0.95, 1.05]** |
| 3 | **有没有信号** | `Normal_holdout` vs 瘤 的逐细胞打分 **AUROC** |

门槛 3 的判读：

| AUROC | 判读 |
|---|---|
| ≥ 0.80 | 有信号，值得进下一轮（**仍须另行签字**） |
| 0.60 – 0.80 | 弱信号 |
| < 0.60 | 分不开 |

**打分定义（必须先定死）**：逐细胞分数 = **常染色体**基因上 `|log2(残差)|` 的**中位数**。
用中位数而非均值，是为了抗少数极端基因的噪声。**不得事后换这个定义。**

**不得**在看到结果后调 `cutoff` / `window_length` / 阈值 / 打分定义（法则 3.2）。

## 八、产物

- `results/03_cnv/infercnv_smoke/gene_order_hg38.tsv`：本项目的基因位置表 + 覆盖报告
- `gene_order_unmatched_genes.txt`：**匹配不上的基因逐个列出**（法则 0：不得静默丢基因）
- `p4_infercnv_counts.tsv` + `p4_infercnv_annotations.tsv`：跑 inferCNV 的输入
- `expr.observations.dat` / 最终 RDS：残差矩阵
- `infercnv_smoke_per_cell.csv.gz`：逐细胞打分
- `infercnv_smoke_summary.json`：判定结果与全部口径常量
- `figures/`：热图 + 打分分布（**图内标签全英文**）

## 九、附录：首跑中止与存储类型修正（2026-09-23，**判读结果尚未产生时**）

**首跑（PID 3948218）已中止**，日志留档 `logs/infercnv_smoke_p4_ABORTED_dgeMatrix.log`，
中间产物留档 `infercnv_out_ABORTED_dgeMatrix/`。

**中止原因**：STEP 08（`subtract_ref_expr_from_obs`）跑满 22 分钟无进展。
查根因：以 `dgCMatrix` 为输入时，`CreateInfercnvObject` 把 `@expr.data` 存成 Matrix 包的
**`dgeMatrix`**；而 inferCNV 自己的 `.subtract_expr()`（`R/inferCNV_ops.R:1742`）按 **base matrix**
语义逐行抽取（`as.numeric(expr_matrix[row_idx, , drop=TRUE])`）。实测同机对照：

| `@expr.data` 类型 | 单行抽取 | 7,668 行外推 |
|---|---|---|
| `dgeMatrix`（首跑实际） | 1.229 秒 | **157 分钟** |
| base `matrix` | 0.002 秒 | **15 秒** |

**处置**：在 `03_cnv/16_run_infercnv_smoke.R` 中 `CreateInfercnvObject` 之后加一行
`infercnv_obj@expr.data <- as.matrix(infercnv_obj@expr.data)`。

**性质声明（重要）**：
1. **不是口径变更，不触碰法则 3.2**。§五 登记的 9 个参数（`cutoff=0.1` / `window_length=101` /
   `denoise=TRUE` / `HMM=FALSE` / `analysis_mode` / `cluster_by_groups` / `write_expr_matrix` /
   `no_plot` / `num_threads`）**一个未改**；§七 三条门槛与打分定义**一个未改**。
2. **不是"换参数凑门"**。改的是内存里的存储类型，不是任何阈值或算法。
   实测 `identical(as.matrix(dgeMatrix), base_matrix)` 为 **TRUE**：数值逐位相同，
   残差矩阵不会因此改变。
3. 本次修正**发生在任何判读数值被看到之前**；首跑除日志与中间对象外**未产出任何结果**，
   故不存在"看着结果改流程"。
4. 该缺陷**属于 inferCNV 自身的实现问题**（对 `dgeMatrix` 输入不自洽），不是我们数据的性质；
   但若将来采纳 inferCNV，此坑必须写进工具说明。
