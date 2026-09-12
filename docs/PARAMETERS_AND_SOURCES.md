# 参数与出处对照表 (PARAMETERS & SOURCES)

> 为 [`WHITEPAPER.md`](WHITEPAPER.md) 每步提供**可溯源参数**；每行标 **核对状态** 与 **出处类型**。
>
> **出处类型**：`S`=源码实测 · `D`=工具默认 · `P`=论文推荐 · `B`=方法学基准 · `C`=本项目约定
> **核对状态**：✅=已从源码/官方文档现取确认 · 🟡=已附权威出处、值待逐条现取复核 · ⚠️=无单一文献（本项目约定）
>
> **更新**：2026-09-12（v2，随范围收窄重写；并修正 v1 两处错误：`infercnv denoise` 默认、CMap `Tau` 定义）

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
| **CopyKAT**（主力，逐样本） | `ngene.chr=5, win.size=25, KS.cut=0.1, distance='euclidean'` | ✅ D | `copykat/R/copykat.R:31`（Gao 2021, *Nat Biotechnol*） |
| infercnv（**可选**交叉验证） | `cutoff=1, window_length=101, HMM=FALSE, HMM_type='i6', denoise=FALSE, cluster_by_groups=TRUE` | ✅ S | `infercnv/R/inferCNV_ops.R:242` — **注：v1 误写 `denoise=TRUE`；默认实为 `FALSE`** |
| 细胞数上限 | infercnv **仅用于 5–10k 细胞/样本子集** | ⚠️ C | 官方已声明 *no longer supported*；全量 200k+ 不可行 |

---

## M3 · 整合（snRNA 内部 + 跨模态）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| scVI（内部） | `n_latent=10, n_hidden=128, n_layers=1, gene_likelihood='zinb', batch_size=128`，`batch=dataset` | ✅ D | scvi-tools 0.15.5 实测默认（本机版本） |
| 跨模态 | **scArches/scANVI 标签迁移**（非 zero-shot；非 `modality 当 batch`） | 🟡 D | Lopez 2018；Lotfollahi 2022 *Nat Biotechnol*（scArches） |
| **禁用** | `modality` 作 batch（sc↔sn 为 **system** 效应且与数据集共线） | ✅ B | scvi-tools（SysVI 专设）；Hrovatin 2025 *BMC Genomics* |
| 评估 | **scIB 完整 panel**：kBET + iLISI + graph-connectivity + PCR；cLISI/ARI/NMI/ASW | 🟡 B | Luecken 2022, *Nat Methods* |
| SCMG 对照臂 | zero-shot 跨数据集整合 + 流形；**不产出逆转/因果** | ✅ S | SCMG 源码实测（无 reversal 方法） |

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
