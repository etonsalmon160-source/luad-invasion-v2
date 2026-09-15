# 参数与出处对照表 (PARAMETERS & SOURCES)

> 为 [`WHITEPAPER.md`](WHITEPAPER.md) 每步提供**可溯源参数**；每行标 **核对状态** 与 **出处类型**。
>
> **出处类型**：`S`=源码实测 · `D`=工具默认 · `P`=论文推荐 · `B`=方法学基准 · `C`=本项目约定
> **核对状态**：✅=已从源码/官方文档现取确认 · 🟡=已附权威出处、值待逐条现取复核 · ⚠️=无单一文献（本项目约定）
>
> **更新**：2026-09-15（v2.0，**新增 M3-A 传统降维/重聚类/亚聚类全节**——此前该部分**零行登记**；
> 补全 M2 的 CopyKAT 缺失参数；更正 M3 的 `batch=dataset` → `sample_id` 并禁用 `patient_id` 作 batch）
> 　·　2026-09-12（v2，随范围收窄重写；并修正 v1 两处错误：`infercnv denoise` 默认、CMap `Tau` 定义）

---

## M1 · QC / 双体（**GSE308103 为 snRNA**）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| QC 阈值（**不得照搬整细胞 scRNA**） | nCount/nFeature 逐样本 **MAD 离群** `nmads=3`, `log=TRUE`, 双尾 | 🟡 D | `scuttle::isOutlier` 默认（McCarthy 2017, *Bioinformatics*） |
| pct_mt 上限 | **< 5**（核；全队列 99% 分位 6.44） | ⚠️ C | 本项目约定（核 mt 本就低） |
| 双体 | `scDblFinder(dbr=NULL, dbr.sd=NULL, aggregateFeatures=FALSE)` **逐样本** | ✅ D | Bioconductor 参考（本次现取） |
| 执行 | `SerialParam()` + `set.seed(1)` | ⚠️ C | 确定性要求 |

> ⚠️ v1 曾列 `nFeature 500–10,000 / nCount 1,000–60,000 / mt<10%`——**那是整细胞 scRNA 约定，不适用于本 snRNA 数据**（实测会砍掉约 30% 的核）。

---

## M2 · 恶性证真 (CNV)

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **CopyKAT**（主力，逐样本） | `id.type='S', cell.line='no', ngene.chr=5, min.gene.per.cell=200, LOW.DR=0.05, UP.DR=0.1, win.size=25, KS.cut=0.1, distance='euclidean', genome='hg20', n.cores=1, output.seg=FALSE` | ✅ S | `copykat` 已装源码现取（2026-09-15）；Gao 2021, *Nat Biotechnol*。**v2.0 补全**：v2 初版漏记 `min.gene.per.cell`/`LOW.DR`/`UP.DR`/`genome` |
| 退化样本处置 | `n_retained<200` → `not_testable`（**绝不填补为 Normal**）；`n_cnv_pos==0` → 作为真实结果上报，**不得放宽 `KS.cut` 重跑**；Normal/AAH `frac_cnv_pos>0.95` 或 IAC `<0.01` → `implausible` | ⚠️ C | 无文献阈值；本项目约定（防事后调参，法则 3.2） |
| **弃跑判据** | `not_testable` + `implausible` 合计 **>30%** → 停并上报 | ⚠️ C | 本项目约定；届时 R2 对该数据集不可满足，**禁用泛上皮 argmax 替代** |
| infercnv（**可选**交叉验证） | `cutoff=1, window_length=101, HMM=FALSE, HMM_type='i6', denoise=FALSE, cluster_by_groups=TRUE` | ✅ S | `infercnv/R/inferCNV_ops.R:242` — **注：v1 误写 `denoise=TRUE`；默认实为 `FALSE`** |
| 细胞数上限 | infercnv **仅用于 5–10k 细胞/样本子集** | ⚠️ C | 官方已声明 *no longer supported*；全量 200k+ 不可行 |

---

## M3 · 整合（snRNA 内部 + 跨模态）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| scVI（内部） | `n_latent=10, n_hidden=128, n_layers=1, gene_likelihood='zinb', batch_size=128`，**`batch=sample_id`** | ✅ D | scvi-tools 0.15.5 实测默认（本机版本）。**v2.0 更正**：原写 `batch=dataset`——单数据集下该键为**常量、无意义**；75 个样本文库才是唯一真实技术批次轴 |
| **禁用** | **`patient_id` 作 batch** | ⚠️ C | 分期**嵌套于**患者内；校正它会抹掉 M4 门判据 ④（患者内配对 AAH vs Normal/AIS）所需的生物学信号 |
| 跨模态 | **scArches/scANVI 标签迁移**（非 zero-shot；非 `modality 当 batch`） | 🟡 D | Lopez 2018；Lotfollahi 2022 *Nat Biotechnol*（scArches） |
| **禁用** | `modality` 作 batch（sc↔sn 为 **system** 效应且与数据集共线） | ✅ B | scvi-tools（SysVI 专设）；Hrovatin 2025 *BMC Genomics* |
| 评估 | **scIB 完整 panel**：kBET + iLISI + graph-connectivity + PCR；cLISI/ARI/NMI/ASW | 🟡 B | Luecken 2022, *Nat Methods* |
| SCMG 对照臂 | zero-shot 跨数据集整合 + 流形；**不产出逆转/因果** | ✅ S | SCMG 源码实测（无 reversal 方法） |

---

## M3-A · 传统降维与重聚类（常规分支，GSE308103）

> **本分支自原始计数重启**；不掺 SCMG（法则 4 / R4）。**执行前登记**（法则 3.1）。
> ⚠️C = 无单一权威出处，本项目约定。

### M3-A.1 预处理

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 归一化 | **scran 池化 size factor**（R 算 → CSV → 回灌 `obs['size_factor']`） | 🟡 P | Lun 2016, *Genome Biol* 17:75（为异质捕获+零膨胀设计，契合 snRNA） |
| 归一化（回退） | `sc.pp.normalize_total(target_sum=1e4)` | ⚠️ C | 预注册回退，**用则登记偏差** |
| 变换 | 归一化计数上 `log1p` | 🟡 P | Ahlmann-Eltze & Huber 2023, *Nat Methods* 20:665 |
| scVI 输入 | **原始计数**（非归一化） | ✅ S | scvi-tools NB/ZINB 似然要求 |
| HVG | `n_top_genes=2000`, **`flavor='seurat_v3'`**, `batch_key='sample_id'` | 🟡 P | Stuart 2019, *Cell* 177:1888（VST）；Heumos 2023, *Nat Rev Genet* 24:550 |
| HVG 依赖 | `scikit-misc==0.2.0`（**2026-09-15 装入 `~/.local`**，装后 `scanpy`/`scvi`/`numpy`/`torch` 版本均未变） | ✅ S | scanpy `seurat_v3` 硬依赖 |
| HVG（回退） | `flavor='seurat'`（log1p 矩阵上，无需额外包） | ⚠️ C | 仅当主依赖失效时启用，**用则登记偏差** |
| 缩放 | `sc.pp.scale(max_value=10)`，**仅 HVG 子集** | 🟡 D | scanpy 默认 |
| PCA | **`n_comps=30` 固定先验** | ⚠️ C | Heumos 2023 建议 30–50；具体值本项目定 |
| PCA 诊断 | kneedle 点（对数解释方差前 50 PC 的二阶差分 argmax）——**仅记录，不替换主值** | ⚠️ C | 防止"眼看拐点"式事后选择 |
| **PCA 护栏** | `n_pcs ∈ {20,30,50}` 三跑，两两 **ARI ≥ 0.90**；违反则三报并升级 | ⚠️ C | 本项目约定 |
| kNN | `n_neighbors=15`, `metric='euclidean'` | 🟡 D | scanpy/Seurat 默认 |
| Leiden | `method='leiden'`, `flavor='igraph'`, `n_iterations=2`, `directed=False` | 🟡 P | Traag 2019, *Sci Rep* 9:5233 |

### M3-A.2 批次校正（双臂）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **Arm A（主）** | **不校正**：未校正 PCA 空间直接聚类 | ⚠️ C | 见下方共线风险 |
| **Arm B（对照）** | `harmonypy` 0.0.10，`batch=sample_id` | 🟡 D | Korsunsky 2019, *Nat Biotechnol*（Harmony） |
| 共线风险 | 每样本＝一病灶＝一分期 ⇒ `sample_id` 与 `stage` **部分共线**；校正可能一并抹掉分期信号 | ⚠️ C | 本项目设计固有 |
| **裁决规则（预注册）** | Arm B 相对 Arm A 的 **ARI < 0.7** ⇒ Arm A 为主，如实报"scVI/Harmony 在分期嵌套设计下的方法学局限" | ⚠️ C | 本项目约定，**不可事后改** |

### M3-A.3 分辨率选择（预注册，不可事后拟合）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 网格 | `resolution ∈ {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0}` | ⚠️ C | — |
| 种子 | `random_state ∈ {1,2,3,4,5}` | ⚠️ C | 确定性要求 |
| 指标1 跨种子稳定性 | 5 种子两两 ARI 均值，**硬约束 ≥ 0.90** | ⚠️ C | Patterson-Cross 2021, *Cell Rep Methods*（ClustAssess 思路） |
| 指标2 跨分辨率稳定性 | 相邻网格 ARI（找稳定平台） | ⚠️ C | 同上 |
| 指标3 谱系覆盖 | 有 ≥1 个法则2 marker 集在 **≥25% 细胞**检出且模块分均值 >0 的簇占比，**硬约束 ≥ 0.90** | ⚠️ C | Heumos 2023（marker 评估） |
| 选择 | `r* = argmax_r [0.5·ARI_seed + 0.5·ARI_xres]`，受两硬约束 | ⚠️ C | — |
| **破平（显式）** | 两分辨率分差 **< 0.01** ⇒ 取**较低**者（少簇、省过度切分） | ⚠️ C | — |
| 无解 | **停在检查点，不放宽阈值** | ⚠️ C | — |

### M3-A.4 注释（双标准交叉）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 标准 A · 基因来源 | 法则 2 规范 marker（6 谱系） | ✅ S | `docs/scientific_rigor_and_audit.md` |
| 标准 A · 方法 | `sc.tl.score_genes` 模块分 → 簇均值 argmax | 🟡 D | scanpy |
| 标准 B · 基因来源 | **独立图谱**：`celltypist` Human_Lung_Atlas | 🟡 P | Travaglini 2020, *Nature* 587:619 |
| 标准 B · 方法 | `CellTypist.annotate` **且** `rank_genes_groups('wilcoxon')` top-50 × 图谱集 Jaccard | 🟡 D/S | 两法互不叠加，正交 |
| 一致性 | 逐细胞 **Cohen's κ** + 逐簇一致率 | 🟡 P | Landis & Koch 1977 释义；**切点本项目定** |
| **切点（预注册）** | `κ≥0.80 且 ≥90%` 通过；`0.60≤κ<0.80` 标记+**强制人工复核**；`κ<0.60` **停** | ⚠️ C | — |
| 分歧 | 全部写入 `annotation_disagreement.csv`，**绝不自动裁决** | ⚠️ C | — |
| 恶性身份 | **仍以 CNV 为准**（法则 0.3），与两标准说什么无关 | ✅ B | 法则 2 |
| 网络 | CellTypist 模型下载须走 `--noproxy` / `ProxyHandler({})` | ⚠️ C | 死代理绕过 |

### M3-A.5 亚聚类（全部主要谱系）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 谱系 | 上皮 / T·NK / B·浆 / 髓系 / 成纤维 / 内皮 | ✅ S | 法则 2 表 |
| 协议 | 每个子集**从原始计数重启全流程**：子集内重选 HVG(2000) → 重算 scran → 重缩放 → 重跑 PCA(30，同护栏) → 重建 kNN(15) → 同网格×5种子×同选择规则 → 子集内双标准注释 | ⚠️ C | **本项目约定**——未找到正式证明该协议避免 double-dipping 的论文。最接近的选择性推断文献（Gao, Bien & Witten 2020, *JASA* 115:1622）讨论的是**选择+检验复用数据**，非聚类。**不冒充有引用** |
| 层级 | `global_label`(L1) 与 `subcluster_label`(L2) **并存，绝不覆盖 L1** | ⚠️ C | — |

### M3-A.6 评估（scIB，独立 venv）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 后端 | **独立 venv** `~/venvs/scib` 装真库 `scib`（**不碰主环境** scvi 0.15.5 / scanpy 1.9.8 / numpy 1.22.4） | ✅ S | 2026-09-15 建；`scib_metrics` 需 Py≥3.9，本机 3.8.10 **不可用** |
| batch / bio | `batch=sample_id`（75 层）；`bio_label` = 冻结后的 L1 注释 | ⚠️ C | — |
| **两遍设计** | Pass 1 只算无标签指标（graph connectivity / iLISI / PCR）+ 无监督聚类；Pass 2 冻结 `bio_label` 后算依赖标签指标 | ⚠️ C | 破循环依赖（标签只在注释后存在） |
| **门规格变更** | 单数据集**无 ground truth** ⇒ Pass 2 的 ARI/NMI 是**未校正 vs 已校正簇标签之间**，非对真值。**须 GP3 批准** | ⚠️ C | 与 PLAN 原文"bio conservation ARI/NMI"含义不同，**必须显式改** |
| 报告 | batch panel 与 bio panel **分开报**，附逐指标原始值（防复合分掩盖） | ✅ B | Luecken 2022, *Nat Methods* 19:41 |

---

## M5 · 空间解卷积 (RCTD)

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| RCTD | `doublet_mode='full'`（Visium 推荐）, `UMI_min=100`, `CELL_MIN_INSTANCE=25`, `counts_MIN=10`, `max_cores` | ✅ S | `spacexr` 2.2.1（Cable 2021, *Nat Biotechnol*） |
| **门** | 用 **RCTD 原生输出**（权重非负 / 分布合理 / 剔除 spot 显式计数） | ⚠️ C | **禁用**把 `normalize_weights()` 的 Σ=1 当门（`full` 模式 `constrain=F`、**无 reject 类别**） |
| 参考 | **GSE308103（snRNA，模态匹配 FFPE↔FFPE）** | ✅ B | 基准：参考模态须与切片匹配 |

---

## M6 · 空间生态位 (BANKSY) + 统计 (Squidpy)

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| BANKSY | `lambda`、`k_geom`（空间邻域）—— **运行时据包文档现取** | 🟡 D | Singhal 2024, *Nat Genet*（运行前复核） |
| 簇数 | **稳定性/共识选择**（非单一 BIC 极小值） | ⚠️ C | BIC 在空间自相关 + 预平滑下会过选 |
| Squidpy 邻域富集 | `nhood_enrichment(n_perms=1000)` → 报 **z + 经验 p=(b+1)/(n+1)** | ✅ S | 源码：**置换标签、只返回 z、无 p 值**；`n_perms` 记录 |
| `co_occurrence` | **不做置换**（`n_perms` 代码在 `spatial_autocorr`） | ✅ S | 源码实测 |
| 距离步长 | **≥ ~100 µm**（Visium 点距；v1 的 55 µm 过小） | ⚠️ C | 平台几何 |

---

## M7 · 靶点

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| DESeq2 | `fitType='parametric'`, `test='Wald'`, `design=~patient+stage`, `alpha=0.1`；患者级 pseudobulk（细胞数≥10 聚合） | 🟡 D | Love 2014, *Genome Biol* |
| Cox/KM | `lifelines.CoxPHFitter(penalizer=0.1)`；校正年龄/性别/分期 | 🟡 D | TCGA-LUAD（表达+临床；**无突变**） |
| **cis-MR** | 工具 F>10；MR P（BH-FDR）<0.05；Egger 截距 P>0.05；Steiger 方向正确 | ⚠️ C | 见 [`M7B_MR_COLOC_TARGET_ANCHORING.md`](M7B_MR_COLOC_TARGET_ANCHORING.md) |
| **coloc** | `PP.H4 ≥ 0.8`（共享因果变异）；`PP.H3` 高 → 判 `distinct` | 🟡 B | Giambartolomei 2014；Wallace 2021（coloc v5） |
| **CMap** | 仅当有真实 LINCS 数据；指标用 **NCS** | ⚠️ | **修正 v1**：Subramanian 2017 中 **Tau 是 0–1 重现性指标，不可能为负**；v1 的 `Tau≤-90` 系误用（旧脚本 `tau=es*100` 为自造） |

---

## M8 · 结构对接

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 探袋 | `fpocket`（V≥350 Å³, Dscore≥0.55）+ **P2Rank** 共识 | 🟡 D | Le Guilloux 2009；P2Rank 2024 对比评测 |
| 对接 | `AutoDock Vina 1.2.7`，`exhaustiveness=8`, `num_modes=9`, `energy_range=3.0` | ✅ S | 源码默认（Trott & Olson 2010） |
| 重打分 | **gnina 1.3**（CNN rescoring，最好的实用精度） | 🟡 P | *J Cheminform* 2025 |
| 物理合理性 | **PoseBusters** 过滤 | 🟡 B | *Chem Sci* 2024 |
| 共折叠 | **Boltz-2**（MIT，`pip install boltz`）——**需 GPU** | 🟡 P | 本轮前沿调研 |
| MD | **缓做**（本机无 GPU/GROMACS）；如做用 **OpenMM 8.x** | ⚠️ C | 本机约束 |

---

## 规则

1. 只有 ✅ 与 🟡（已附出处）可入白皮书；⚠️ 如实标「本项目约定」。
2. 🟡 项在**实际运行前**逐条现取确认并回填本表。
3. 凡未能核实者**不写数值**，只写"按工具默认（见出处）"。
4. 阈值一经确定**不得事后调整**（防 p-hacking）；如确需调整，须记录理由与时间。
