# SCMG 分支 · 第一腿预注册：零样本投影到全局流形

**登记时间**：2026-09-23（**本文件 mtime 早于本腿任何分期分层结果的计算**）
**性质**：🔬 **SCMG 分支（M3-B / GP9）的第一腿**。**不掺传统算法**（铁律 R4）；**不输出逆转/因果**（铁律 R5）。
**授权**：用户 2026-09-23「先把scmg该引入了，那条路该走了，看看结果」。
**范围声明**：本腿**不是**「迁徙分析」。迁徙需用亚型边界定义起点/终点 ⇒ 那是 `docs/PARAMETERS_AND_SOURCES.md:357`
登记为 🔴 的 **R4 边界问题**（「GP8a 的亚型标签以何种形式进入 SCMG」），**尚未签字，本腿不碰**。

---

## 一、要回答的三个问题

1. **仪器/参照合理吗？**（工程层面）
2. **分期轴在 SCMG 流形上还在不在？**（= M3-B 注册判据「分期可分性被保留」）
3. **各期核落在全局流形的哪些参照细胞状态上？**（状态刻画）

**不问**「谁是恶性」——参照流形**无肿瘤态**（2026-09-12 已核实），SCMG 只能回答"最像哪个**正常/发育**细胞状态"。
恶性仍由 CNV 证真（R2），不在本腿。

## 二、工具与资产（冻结）

| 项 | 值 | sha256 |
|---|---|---|
| SCMG 源码 | `/home/eto/scmg_workspace/SCMG`（官方 `xingjiepan/SCMG`）| git `29f44c98c5621d575a058b549e0fdde0ba31a730` |
| 编码器全模型 | `models/embedder/model.pt` | `01bd51e6d8a54486e5d3bc412c980466031c5f0e91c8ec2fe0474c4c3150f599` |
| 编码器权重 | `models/embedder/best_state_dict.pth` | `e585faae3fb2f7805880cebdfcb18c9550a636381e8af211ee3bfcab3742119f` |
| 标准基因表 | `scmg/data/standard_genes.csv` | `44c436d39ea978a4787132a1ddd92469aa3fbdef793b69d99f15f7c5dd09de78` |
| 参照全局流形 | `hf_data/ref_global_cell_state_manifold.h5ad` | `27669cfbf96a7274c8bf2a6ea23e16b8f744ecc22143229f8c4b34f5f3133d4b` |

**运行环境（2026-09-23 实测，非假设）**：**用主环境**，**不新建、不装包、不碰基础库**。
主环境 py3.8.10 / numpy 1.22.4 / scanpy 1.9.8 / torch 2.4.1+cu118 / pandas 2.0.3 / anndata 0.9.2，
`sys.path.insert` 注入源码后 `from scmg.model.contrastive_embedding import CellEmbedder, embed_adata`
等**全部 import 通过**、`pkg_resources.resource_filename` 解析正常。
**无 GPU**（`torch.cuda.is_available() == False`）⇒ 全程 CPU。
⚠️ `pyproject.toml` 钉的 `numpy==1.26.4 / scanpy==1.10.4 / torch==2.3.1 / pandas==1.5.3` 要求 Python ≥3.9，
**装进主环境必然炸主库**；实测源码在**主环境现有版本**下可加载。**这是对 `pyproject` 固定版本的偏离，如实登记（§六.7）。**

## 三、输入（冻结）

| 项 | 值 |
|---|---|
| 表达对象 | `results/02_expression/gse308103_counts_paperqc.h5ad`，sha256 `a276cd1af69a4620ff9e3d64f4b93f33461e646537c5f2c9da1aff8c9150b9de` |
| 维度 | **413,697 核 × 18,069 基因**（GP0 论文 QC 口径） |
| `var.index` | **基因符号**（`SAMD11`…），**无 Ensembl 列** |
| 分期分层 | Normal 94,506 / AAH 35,283 / AIS 81,555 / MIA 14,266 / **IAC 188,087**；**23 患者 / 75 样本**（R1 权威口径） |
| 结构 | `sample_id`↔`stage` **一对一嵌套**；但**同一患者同时含 Normal 与 IAC** ⇒ 存在**患者内配对** |

**基因映射（冻结）**：**用 SCMG 自带的 `standard_genes.csv` 的 `human_name → human_id`**，
**不自造映射表**（用户 2026-09-23「亚型标签能不能用他们自己的映射器」；此处同旨）。
大写归一后实测覆盖率 **16,104 / 18,069 = 89.13%**。
未映射 1,965 个**逐个列出**（`scmg_unmapped_genes.txt`，法则 0：不得静默丢基因）。

**重名陷阱（预先定死处置）**：SCMG `human_name` 有 **11 个重名符号**（18,108 个 id 只 18,074 个唯一名）。
我们只撞上 **1 个：`MKKS`**（对应 3 个 Ensembl：`ENSG00000125863` / `ENSG00000285508` / `ENSG00000285723`）。
**处置（冻结）**：取 `standard_genes.csv` 中**首次出现**的那一行，**确定性、可复现**，并单独登记。
**禁止**改用别的映射源凑一个"更好的" MKKS。

## 四、参数（**计算前登记，事后不许改**，法则 3.1 / 3.2）

| 参数 | 取值 | 依据 |
|---|---|---|
| `device` | **cpu** | 本机无 GPU |
| `batch_size` | **8192** | 实测 512 与 8192 输出**逐位相同**；8192 使 413,697 核约 **18.3 分钟**（512 需约 62 分钟） |
| `pre_normalized` | **False** | 用 SCMG 内部的 Cp10k + `log1p`（`get_Xs_from_anndata`），**不自作归一化** |
| `CellTypeSearcher.emb_key` | **`X_scmg`** | 官方 tutorial 默认 |
| `CellTypeSearcher` 方法 | **只用 `search_ref_cell`** | 见下 |
| `n_jobs` | **20** | 本机 20 核 |

🔴 **禁用 `search_ref_cell_types`**：它的 soft 版先算 `pairwise_distances(413,697 × 133,061)`
≈ **2.2×10¹¹ 个 float ≈ 880 GB**，本机必然 OOM。**不是"慢"，是"不可能"**，故禁用。

## 五、判读（**看到结果前定死**）

| # | 门槛 | 通过判据 |
|---|---|---|
| **S1** | **仪器有效** | 基因覆盖率 **≥ 85%**（实测 89.13%）；`X_ce_latent` 全为有限值；非退化（逐维标准差 > 0，且不是所有细胞同一个向量） |
| **S2** | **参照合理性** | 我们核的**最近参照细胞**落在**肺相关 tissue** 的占比。预设肺 tissue 集 = {`Lung epithelium`, `lung`, `lung parenchyma`, `respiratory airway`, `lingula of left lung`, `Gut and lung epithelium`}（参照中占比 **15,823/133,061 = 11.89%**）。<br>**判据：占比 ≥ 10%**（宽松下界）。**若 < 10% ⇒ 报「投影未落在肺区」，本腿作废，不得解读任何下游分布** |
| **S3** | **分期可分性被保留**（M3-B 注册判据） | 在 SCMG 512 维潜空间上，**留一患者（leave-one-patient-out）最近质心分类器**区分 **Normal vs IAC**，报 pooled **AUROC** |

S3 的判读（沿本项目既有口径，保持可比）：

| AUROC | 判读 |
|---|---|
| ≥ 0.80 | 分期信号在 SCMG 流形上**被保留** |
| 0.60 – 0.80 | **弱**保留 |
| < 0.60 | **未被保留**（"全糊一块"） |

**S3 为什么用留一患者 + 最近质心（先写明理由，非看到结果后选）**：
- **留一患者**：Normal 与 IAC 来自**同一批患者**，患者级批次是最强混淆；留一患者测的是**跨患者泛化**，
  跑分器**看不到**被评患者，故不是靠患者身份蒙对。
- **最近质心**而非 kNN：留一 23 折 × 28 万细胞的 kNN 在 CPU 上不可行；质心是 O(n)，且**先于结果定死**。
- 潜空间是**零样本、无监督**产生的（**未曾用过分期标签**）⇒ 本检验**不循环**。

**描述性产出（不参与判读）**：各期核映射到的**参照细胞类型分布**（top 15）+ 逐期「肺相关参照类型」占比表。

## 六、偏差登记（**必须带进报告**）

1. **模态不匹配**：输入是 **snRNA（细胞核，FFPE）**，参照流形以 **scRNA（全细胞）** 为主。
   SCMG 与 2026 基准均**未在 snRNA 上验证**过。
2. 🔴 **参照流形无肿瘤态**（2026-09-12 已核实）⇒ SCMG **不能**给恶性标签，
   只能给「最像的**正常/发育**细胞状态」。**"最像 AT2" ≠ "是正常细胞"**。
3. **参照是全身泛组织**：331 组织 / 31 大类 / 36 数据集（含小鼠脑、人胎盘、视网膜），**肺只占 11.89%**。
4. **分期↔样本嵌套**；S3 用留一患者缓解，但同患者内仍共享技术批次。
5. **单数据集、无 ground truth** ⇒ **本腿不关闭 GP9**。GP9 的注册过门要**与 M3-A 做 scIB 完整口径对照**
   （batch panel + bio panel 并报），须 M3-A 完成后另做。**本腿只出 SCMG 一侧**。
6. 🔴 **不得把 SCMG 的 797 类全局状态标签说成 LUAD 亚型标签**。两者是**不同的标签体系**：
   SCMG 给的是「全身细胞状态里最近的那个」，**没有** KAC / AAH 亚型 / 肿瘤态。
7. **偏离 `pyproject.toml` 固定版本**（见 §二）：用的是主环境现有版本（numpy 1.22.4 / scanpy 1.9.8 / torch 2.4.1），
   非 SCMG 作者钉的 numpy 1.26.4 / scanpy 1.10.4 / torch 2.3.1。**理由**：后者要求 Python ≥3.9，装则炸主库；
   实测现有版本可加载。**未改任何数值口径**。
8. 全局流形是**预训练资产**，其训练数据（Tabula 类图谱）与本项目无关；**不得**把参照细胞当成我们队列的细胞。

## 七、产物

- `scmg_unmapped_genes.txt`：未映射到 SCMG 标准基因集的 1,965 个基因（逐行，法则 0）
- `scmg_latent.npy`：413,697 × 512 的 `X_ce_latent`
- `scmg_projection.csv.gz`：逐核 `ref_cell` / `distance` / `projected_cell_type` / `projected_tissue` / `umap_x` / `umap_y`
- `scmg_branch_summary.json`：S1/S2/S3 判定 + 全部冻结常量 + 输入哈希
- `figures/scmg_1_global_manifold_projection.png`：参照流形 UMAP 底图 + 我们各期核投影
- `figures/scmg_2_stage_composition.png`：各期映射到的参照细胞类型构成
- 脚本：`06_scmg/01_project_to_global_manifold.py`、`06_scmg/02_analyze_scmg_branch.py`
- 跑法：`python3 06_scmg/01_project_to_global_manifold.py`（本腿只有这两步；§二 已核实无需额外依赖）
