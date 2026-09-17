# 参数与出处对照表 (PARAMETERS & SOURCES)

> 为 [`WHITEPAPER.md`](WHITEPAPER.md) 每步提供**可溯源参数**；每行标 **核对状态** 与 **出处类型**。
>
> **出处类型**：`S`=源码实测 · `D`=工具默认 · `P`=论文推荐 · `B`=方法学基准 · `C`=本项目约定
> **核对状态**：✅=已从源码/官方文档现取确认 · 🟡=已附权威出处、值待逐条现取复核 · ⚠️=无单一文献（本项目约定）
>
> **更新**：2026-09-16（**v3.0 — 全面对齐源论文口径**。用户当日决策「改用论文 QC 重建」+「对齐论文全套」。
> 新增 §M3-A.0 源论文与检测平台（权威出处逐字引文）；M1 的 QC 口径由"逐样本自适应 MAD"**改为论文固定全局阈值**，
> 原口径降级为敏感性臂；M3-A.1–A.5 的 scran / HVG-2000 / PC-30 / kNN-15 / Leiden 全套**已被
> SCTransform(v2) / HVG-3000 / PC-50 / kNN-20 / **Louvain** 取代**。
> 连带三处**结构性**变更（非换值）：① **M3-A.2 主口径改为 Harmony**，原本的"不校正＝主"降为**敏感性臂**；
> ② **M3-A.3 分辨率网格收窄到论文区间 0.5–0.8**，旧 `{0.2…2.0}` 仅作诊断；③ **M3-A.3 新增「指标4 · AAH 吸收护栏」**
> （可证伪、跑前固定），回应"须考虑 AAH 脆弱性"的要求——实测已排除计数层面脆弱，故护栏设在**聚类吸收**层面）
> 　·　2026-09-15（v2.0，新增 M3-A 传统降维/重聚类/亚聚类全节；补全 M2 的 CopyKAT 缺失参数；
> 更正 M3 的 `batch=dataset` → `sample_id` 并禁用 `patient_id` 作 batch）
> 　·　2026-09-12（v2，随范围收窄重写；并修正 v1 两处错误：`infercnv denoise` 默认、CMap `Tau` 定义）

---

## M3-A.0 · 源论文与检测平台（**全流程的权威锚点**，2026-09-16 建立）

> 本节存在的原因：本项目此前**没有**把源论文的原始方法学落盘，参数只能靠"通用权威文献"逐项拼凑，
> 导致口径与数据实际来源脱节。2026-09-16 重新抓取原文并**逐字核对**后补上本节。
> **凡本节与其它节冲突，以本节为准**（它是这批样本的**原生**方法学）。

| 项 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 论文 | Peng F, Sinjab A, Dai Y, Treekitkarnmongkol W, Yang S, … Wang L, Kadara H. *Multimodal spatial-omics reveal co-evolution of alveolar progenitors and proinflammatory niches in progression of lung precursor lesions.* **Cancer Cell** 2026;44(2):321–339.e13 | ✅ S | DOI `10.1016/j.ccell.2025.10.004` · PMC12980502 · PMID 41202811 · license CC BY-NC-ND |
| **检测平台（决定一切 QC 解释）** | **10x Chromium Fixed RNA Profiling（探针法）+ FFPE 卷片**，非标准 3′/5′ | ✅ S | 原文：*"using FFPE scrolls consecutive to ST-analyzed sections… Chromium Fixed RNA profiling (CG000632)"*、*"Chromium Fixed RNA Human transcriptome kits"*。**解释了 18,082 的固定基因面板** |
| 测序深度目标 | ≥ 50,000 mean read pairs / cell | ✅ S | 原文 snRNA-seq 方法节 |
| 样本 | snRNA-seq **75 样本 / 25 患者**（其中 23 例与 ST 配对） | ✅ S | 原文：*"integrated ST data with snRNA-seq profiles from a matched set of 75 samples… from 23 of the 25 patients"* |
| **逐分期样本数（论文 Figure 3H 图注）** | normal **24** / AAH **9** / AIS **14** / MIA **4** / LUAD **24** | ✅ S | 与原文字面一致。**★ 与本项目 `cohort_registry` 的 `Normal 24 / AAH 9 / AIS 14 / MIA 4 / IAC 24` 完全吻合 ⇒ 对 R1 分期权威性的独立外部验证** |
| 保留核数 | **401,635**；上皮 **139,663** | ✅ S | 原文：*"After rigorous QC and filtering, 401,635 nuclei… were retained"*、*"analyzed 139,663 nuclei from epithelial cells"* |
| 上皮亚群命名（**论文的**） | ciliated / club / basal / AT1 / AT2 / **AIC**（alveolar intermediate）/ **KAC**（KRT8-high alveolar intermediate）/ tumor | ✅ S | 原文。⚠️ **本项目 GP8a 不采用此命名体系**（用户 2026-09-17 改口径为经典 marker，见 §M3-A.5c）；此处仅记录论文原话 |
| **分期的判定方式（对 M4 门至关重要）** | **病理判读**，非无监督聚类 | ✅ S | 原文：*"H&E-stained slides… were assessed by **3 independent experienced pathologists** to determine the presence of and type of lesion (invasive, precursor)"*；ROI 用 HALO 标注，*"Paired samples were selected and classified as invasive, precursor, or normal"* |
| MIA 定义（原文逐字） | lepidic 为主、≤3 cm、浸润成分 ≤0.5 cm、无坏死与脉管侵犯 | ✅ S | 原文引用 WHO 第 45 号文献 |
| 亚群划分是否逐患者独立 | **是**（仅 **Visium ST** 节；snRNA-seq **不是**） | ✅ S | 原文 ST 节：*"Clustering analysis was performed independently for each patient"*。⚠️ 勿把此句误用到 snRNA |

> ⚠️ **对 M4 门的直接推论**：论文的 normal / AAH / AIS / MIA / LUAD 标签来自**病理学**，不是聚类产物。
> 因此本项目 **M4 判据 ⑤「AAH 身份 CNV/标记证真」不是"复现论文的聚类"**，而是在已有病理标签之上做**独立的分子证真**。
> 两者是**互补关系**，不是同一件事——不得把"聚类没分出 AAH"当作"AAH 不存在"的证据。

---

## M1 · QC / 双体（**GSE308103 为 snRNA**）

> 🔴 **2026-09-16 口径变更（用户决策「改用论文 QC 重建」）**：主口径由本项目自定的
> **逐样本自适应 MAD** 改为**源论文的固定全局阈值**。理由：论文的阈值是这批**同平台同批次**样本的
> 原生口径，且其实测保留数（401,635）与我们套用后的结果（413,697）**只差 3.0%**；
> 而旧的自适应 MAD 口径比权威口径**多留 36%的核**（648,945 vs 413,697），属系统性过宽。
> 旧口径**不删除**，降级为**敏感性臂**（`results/02_expression/gse308103_counts.h5ad` 保留）。

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **QC 主口径（新）** | `nFeature >= 500` **且** `nCount >= 1000` **且** `pct_mt < 20`（论文写 "≥20% were filtered out"）| ✅ S | **源论文 Methods 逐字**（见 §M3-A.0）；脚本 `01_qc/11_regate_paper_qc.py` |
| ↳ 低复杂度预过滤 | `nFeature >= 200`（先剔 cell debris / empty drops） | ✅ S | 同上。**被 500 门完全吸收**，单独施加不改变结果 |
| ↳ 基因过滤 | 在**新细胞集上**检出 **< 3 个细胞**的基因剔除 | ✅ S | 同上；实测剔除 **13 / 18,082** 个（`02_expression/04_rebuild_expression_paperqc.py`） |
| ↳ **实测结果** | **413,697 核**（论文 401,635，**+3.00%**）；18,069 基因；nnz 655,222,534 | ✅ S | `results/01_qc/regate_paper_qc_report.json` |
| ↳ **逐门槛实际约束力** | 线粒体门在 M1 单细胞集内砍 **0 个**核；约束力**全在** `nFeature>=500`（−68,694）与 `nCount>=1000`（−166,554） | ✅ S | 同上。⚠️ 因 snRNA 线粒体本就低（全队列 pct_mt 中位 **0.60%**、99.9 分位 30.9%），该门几近空转 |
| 双体 | `scDblFinder(dbr=NULL, dbr.sd=NULL, aggregateFeatures=FALSE)` **逐样本**；02_annotate 补未建阈值样本 | ✅ D | Bioconductor 参考（本次现取） |
| 执行 | `SerialParam()` + `set.seed(1)` | ⚠️ C | 确定性要求 |
| **🔴 与论文的登记偏差：双体调用器** | 论文 **Scrublet**（Wolock 2019）；本项目 **scDblFinder** | ⚠️ **偏差** | **铁律 R2 钉死 scDblFinder，优先于"对齐论文"。** 缓解：M1 已跑 Scrublet 敏感性臂（`01_qc/05_doublets_scrublet.py` → `results/01_qc/gse308103_sensitivity_doublet_rate.csv`），两口径并存可查 |
| **论文双体的第二道（本项目尚未做）** | *"based on cluster distribution and marker gene expression, doublets forming distinct clusters with hybrid expression features were also removed"* | 🟡 P | 原文。属**聚类后**步骤，本项目对应到 **GP6 注释之后的跨谱系共表达簇剔除**，须在 GP6 显式执行并登记 |
| ~~QC 阈值（旧，已降级）~~ | ~~nCount/nFeature 逐样本 MAD 离群 `nmads=3`, `log=TRUE`, 双尾~~ | 🔶 **降级为敏感性臂** | `scuttle::isOutlier`（McCarthy 2017, *Bioinformatics*） |
| ~~pct_mt 上限（旧）~~ | ~~`< 5`（核；全队列 99% 分位 6.44）~~ | 🔶 **降级为敏感性臂** | 本项目约定（核 mt 本就低） |

> ⚠️ v1 曾列 `nFeature 500–10,000 / nCount 1,000–60,000 / mt<10%`——**那是整细胞 scRNA 约定**。
> **注意**：论文的 snRNA 口径（500 / 1000 / 20%）与本项目 v2 的"核数据不得照搬 scRNA 阈值"结论**方向相反**，
> 但**不矛盾**：论文的 500/1000 是**下界**（剔低质量），v1 的 500–10,000 是**双边界**（上界 10,000 会砍掉高深度核）。
> 实测佐证：套用论文门后 AAH 留存率 **65.48% 反而高于** Normal 的 **62.05%**（+3.43 pp），未见对早期病变的选择性丢失。

---

## M2 · 恶性证真 (CNV)

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **CopyKAT**（主力，逐样本） | `id.type='S', cell.line='no', ngene.chr=5, min.gene.per.cell=200, LOW.DR=0.05, UP.DR=0.10, win.size=25, KS.cut=0.1, distance='euclidean', genome='hg20', n.cores=1, output.seg=FALSE, plot.genes='FALSE'` | ✅ S | `copykat` 已装源码现取（2026-09-15）；Gao 2021, *Nat Biotechnol*。**v2.0 补全**：v2 初版漏记 `min.gene.per.cell`/`LOW.DR`/`UP.DR`/`genome`。<br>🟢 **2026-09-16 18:00 定稿：`UP.DR` 取包默认 0.10**（即 `LOW.DR=0.05, UP.DR=0.10`）。曾一度改为 0.05 以求全队列统一，**该改动已随 840 地板的重算而撤回**（理由见下行末段）。撤回后本项目参数**与 `copykat` 包默认、作者 README 完全一致**，不再有任何偏离。用户 2026-09-16「按照目前最合理的情况做」授权按此定稿。<br>⚠️ **生效差异**：源码 `:57-60` 的无条件覆写仍会在判据基因集 <7000 的样本上把有效 `UP.DR` 压到 0.05。地板 840 之后这样的样本只剩 **2/75**（`P16_AIS` 6,949、`P8_Normal` 6,272；见 `prereg_gene_tiers_floor840.csv`）。GP2 对每个样本**逐条申报有效值**，不假设它是 0.10 |
| ↳ **`norm.cell.names`（锚定开关；🔴 未定，本行是 GP2 的第二道未决口）** | 当前 **`""`（不锚定）** = 让 copykat **自己猜**正常基线。已装源码 `copykat.R:136` 的**分支条件是 `length(norm.cell.names) > 1`** —— 🔴 **只给 1 个细胞名会静默落进"猜基线"分支**，不是"用这 1 个细胞当基线"。**要锚定必须一次给 ≥2 个细胞名** | 🔴 | `tools/copykat/R/copykat.R:136-184`（已装 1.2.5 源码现取）；分支后果见 [`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §9.7 |
| ↳ 两个仅影响绘图的参数（登记以免误当作结果相关） | `cell.group=""`（只画热图侧栏，源码 `:527-535`）；`test.emd="FALSE"`（**字符串比较**，源码 `:539`，只控制是否画 emd 聚类） | ✅ S | 源码现取 2026-09-17；**两者均不改变 CNV 判定** |
| ↳ **基因注释表**（`genome='hg20'` 真正用到的东西） | 包内 sysdata 两张：`full.anno` **56,051 × 7**（`abspos, chromosome_name, start_position, end_position, ensembl_gene_id, hgnc_symbol, band`）、`DNA.hg20` **12,205 × 3**（`chrom, chrompos, abspos`，分段用的 bin 网格）。🔴 **两表均无版本号、无下载日期** ⇒ 出处只能记到「copykat 1.2.5 内置」，**无法追溯到具体 biomaRt/Ensembl 版本**。`id.type='S'` 走 `hgnc_symbol` 匹配 | ⚠️ C | 🔴 **出处不完整，如实登记** |
| ↳ **静默丢弃的基因**（法则 0「静默默认陷阱」） | 本项目 18,069 基因中 **17,352 命中（96.03%）、717 被静默丢弃**。**717 个不是随机丢**：① 整条**线粒体基因组**（`MT-ND1`…`MT-CYB`，13 个）；② 大量 **HGNC 改名基因**——我们的基因轴用新符号，`full.anno` 是旧导出（`CCN1`←CYR61、`ATP5PB`←ATP5F1B、`PRUNE1`←PRUNE、`TENT5B/C`←FAM46B/C、`DIPK1A`←FAM69A、`SHISAL2A`←FAM159A、`MACO1`←TMEM57、`MINDY1`←FAM63A、`KHDC4`←KIAA0907、`TLCD4`←TMEM56、`INKA2`←FAM212B、`TAFA3`←FAM19A3、`TCTEX1D4`←DYNLT4` 等）⇒ **系统性偏向"近年改过名"的基因**。逐样本的有效基因集另受 `gene_min_cells=3` 与地板 840 影响，**GP2 须逐样本申报** | ⚠️ C | 实测 2026-09-17，脚本口径见本行下注 |
| ↳ 🔴 **Y 染色体完全缺失** | `full.anno$chromosome_name` **只取 1–22 与 23**；23 = **X**（反查：`XIST` band q13.2、`KDM6A` p11.3）。**Y 基因一个都不在表里**（`RPS4Y1` `DDX3Y` `UTY` `USP9Y` 全部缺失）⇒ **copykat 在本配置下无法检测 Y 丢失**。Y 丢失是 LUAD 常见事件 ⇒ 本项目 CNV 判定**只覆盖 chr1–22 + X**，**报告须写明，不得表述为"全基因组 CNV"**，也不得把"未见 Y 丢失"当成证据 | ⚠️ C | 实测 2026-09-17；**该限制来自 copykat 内置表，非本项目选择** |
| ↳ `!duplicated(hgnc_symbol)` 步骤 | 源码 `:19` 会去重 `hgnc_symbol`，但 `full.anno` 56,051 行**零重名** ⇒ 该步**丢 0 行**，无效但无害（登记以免日后误以为它会吞基因） | ✅ S | 实测 2026-09-17 |
| ↳ 实测脚本（可复跑） | `Rscript -e 'suppressMessages(library(copykat)); e<-new.env(); data("full.anno",package="copykat",envir=e); fa<-get("full.anno",envir=e); ours<-readLines("results/04_integration/seurat_io/gene_names.txt"); cat(sum(ours %in% fa$hgnc_symbol),"/",length(ours),"\n")'` | ✅ S | 本行全部数字由此得出 |
| 🔴 **数据源本身只有 18,082 个基因（2026-09-17 新发现）—— 成因 = 探针法（Flex）** | 作者上传的每个 `.raw_counts.mtx.txt.gz` **固定 18,082 行**（= 基因轴），按**基因组位置**排序（`SAMD11` chr1p36.33 → … → `KDM5D`/`EIF1AY`/`DAZ2` chrY → `MT-ND1`…`MT-CYB` chrM）。标准 10x `GRCh38` 参考约 **36,601** 个 feature ⇒ **约一半转录组不在文件里**。实测**完全缺失**：**GAPDH**、**全部 `RPL*`/`RPS*` 核糖体基因（0 个）**、`MALAT1`、`TPI1`、`PTMA`、`MT-CO1`、`BCAS1`、`LINC00511`，以及本项目面板基因 `KRT18`、`SFTPA2`、`CD8B`、`FCN1` | ⚠️ C | **成因已查明（探针面板，非生信过滤）**：源论文 STAR Methods 原文 —— samples processed "for Chromium Fixed RNA profiling (CG000632)"、建库用 "**Chromium Fixed RNA Human transcriptome kits** (10x Genomics, 4 rxns × 4 BC 1000475 / 4 rxns × 16 BC 100476)"、"Chromium Fixed RNA Profiling Reagent kits (CG000527)"；GEO `Series_overall_design` 亦写 "using the 10x Genomics Chromium Fixed RNA protocol"。= **10x Flex 探针捕获**，探针集本身约 1.8 万基因 ⇒ **特征集由实验方法决定**。<br>⚠️ **但论文与 GEO 均未给出该探针集的基因清单或排除规则** ⇒ 具体为何缺 GAPDH/RPL/RPS **仍无法从公开材料核实**，如实登记为**未说明的探针集内容**，不假设 |
| ↳ 佐证：基因轴是**固定面板**而非按样本过滤 | 三个样本（首个 / 第 30 个 / 末个）的基因列 `md5` **完全一致 = `ef30771526daf1d7`**，行数均 18,083（含 barcode 行）。若是"按表达过滤"，各样本基因轴必不同 | ✅ S | 实测 2026-09-17；与探针法解释一致 |
| ↳ 上游对齐（论文原文，本案未用到） | "Raw data were pre-processed (demultiplex cellular barcodes, read alignment and generation of gene count matrix) using **Cell Ranger** Single Cell Software Suite provided by 10x Genomics. For alignment of reads, human reference **GRCh38 (hg38)** was used." ⇒ **Cell Ranger 版本号、注释版本（GENCODE 版本）、探针集清单，论文与 GEO 三者均未给** | ⚠️ C | 源论文 STAR Methods「snRNA-seq data analysis」；GEO `Sample_data_processing` 逐字相同。**登记为出处缺口，不推测版本** |
| ↳ 在 / 不在（供面板取舍） | **在**：`ACTB`、`EPCAM`、`KRT8/19`、`CDH1`、`NKX2-1`、`SFTPC`、`SFTPA1/B/D`、`NAPSA`、`AGER`、`CAV1`、`PDPN`、`SCGB1A1/3A2`、`FOXJ1`、`PIFO`、`TP63`、`KRT5`、`KRT15`、`MKI67`、`TOP2A` | ✅ S | 实测 |
| ↳ 上一行「96.03%」的**正确读法** | 该比例的分母是**本项目那 18,069 个基因**，**不是全转录组**。相对标准 `GRCh38` 全转录组的真实覆盖约 **49%**（18,082 / 36,601）。⇒ **报告与图注不得表述为「CopyKAT 覆盖全转录组 96%」**，须写作「覆盖本数据可用基因的 96.03%」 | ⚠️ C | 同上；**本次更正** |
| ↳ 对 CNV 的连带影响 | ① 核糖体基因常被 CNV 工具用作**基线/参照**，本数据**一个都没有** ⇒ 基线构成与常规不同，**须在 GP2 报告中说明**；② 作者矩阵含 3 个 Y 基因（`KDM5D`/`EIF1AY`/`DAZ2`），但**copykat 注释表无任何 Y 基因**（见下行）⇒ **Y 仍不可分析**，该 3 个基因会被静默丢弃，计入 717 | ⚠️ C | 实测 2026-09-17 |
| **`UP.DR` 强制与 `LOW.DR` 相等**（**预注册：2026-09-16，GP2 开跑前**；⚠️ **书面签字待补**） | `LOW.DR = UP.DR = 0.05` | ⚠️ C | **偏离 `copykat` 包默认**（默认 `LOW.DR=0.05, UP.DR=0.1`）。理由：①`UP.DR` 名不副实，源码 **`:189`**（及 `:192`）用作**下界** `DR2 >= UP.DR`，是驱动最终分割的基因检出率门槛；②源码 **`:57-60`** `if(nrow(rawmat) < 7000){ UP.DR <- LOW.DR ; print(...) }` 是**无条件覆写**，显式传参**拦不住**，只有两值相等才能使有效阈值全队列一致；③实测 75 样本（`results/03_cnv/prereg_gene_tiers.csv` 的 `genes_ge_*` 列）：统一 **0.05** → 最终基因数 min **3,157** / p25 5,903 / 中位 **6,742** / max 9,145，**0 个**样本 <2,000；统一 **0.10** → min **1,791** / p25 3,756 / 中位 **4,494** / max 7,523，**1 个**样本 <2,000，`win.size=25` 按基因序平滑会跨极大基因组距离、分割退化；④不强制则 **23** 个样本用 0.05、**52** 个用 0.1，且覆盖率与分期相关（Normal 最高），构成沿 M4 进展轴的混杂。**此为预先决定，非事后调参**（法则 3.2）。<br>⚠️ **已刷新（2026-09-16）**：本行早前引用的 min 3,528 / 中位 7,429、min 1,231 / 中位 3,723、3 个样本 <2,000、27/48 分布，是 ③ 的**建模值**版本；`02_prereg_gene_tiers.R` 重写为调用 `copykat_chain()` 后 `genes_ge_*` 已换成**精确值**，上表为精确值。**结论方向不变**（0.05 仍更保守），但决策时依据的数字与现值不同，如实记此差异。<br>⚠️ **仍未消除的反对意见**：作者 README 明写 "I need to make sure that LOW.DR is smaller than UP.DR"，本项目此决定**与该句直接冲突**，不假装有文献支持。<br>⚠️ **签字状态**：`results/03_cnv/GP1_report.md` §7 曾把本项列为【需裁决】。**本行不主张已获签字** —— 2026-09-16 复核时用户询问"证据是否不支持"，本行据此复核确认**证据方向未变**（0.05 仍更保守），但**书面 go/no-go 仍待用户给出**<br>🔴 **2026-09-16 17:20 证据翻转（同一日内，加算 840 地板后）**：本行 ③④ 两条理由**建立在地板之前的口径上**。在地板 840 之后重算（`logs/tiers_floor840.log`、`prereg_gene_tiers_floor840.csv`）：<br>· 7000 覆写从 **23/75 塌到 2/75**（21 个样本 0.05→0.10，**无一反向**）。剩余两个是 `P16_AIS`（1,471 细胞，判据基因集 6,949）与 `P8_Normal`（126 细胞，6,272）。<br>· 因此 ④ 的"分期混杂"**基本消失**：Normal 45.8%→4.2%、MIA 50%→0%、IAC 25%→0%。<br>· 代价却**反向变大**：统一为 0.05 需要动的样本数 **52 → 73**（因为覆写后 73/75 的有效值已是 0.10），而买到的"均匀性"从 31% **缩到 2.7%**。<br>· ③ 的"0.05 更保守"仍未变（min 3,157 vs 1,791 量级不变），但**它已不是主要论据**，因为 73/75 本来就要用 0.10。<br>🔴 **本项目据此撤回原推荐**：**不做统一**，令 `LOW.DR=0.05, UP.DR=0.10` 走 `copykat` 默认 + 源码 :57 覆写；并在 GP2 报告中**逐样本申报有效 `UP.DR`** 及其与分期的关联（关联已弱化到可忽略）。**理由**：这笔交易的方向反了——代价从 52 涨到 73 个样本、收益从 31% 缩到 2.7%，且原④的混杂论据已被地板消除。**这是撤回，不是新增签字项**；原预注册从未获书面签字，故不构成事后调参。<br>⚠️ 仍存的分歧：作者 README "LOW.DR is smaller than UP.DR" **与默认参数一致**，故撤回后反而**与作者一致** |
| **`n.cores` 是结果相关参数，不得当加速旋钮调** | `n.cores=1`（全队列固定） | ✅ S | **实测证伪**（2026-09-16，R 4.2.2 / 20 核，见 [`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §4.3）：同一 `set.seed(1234)` 下 `mclapply(1:8, runif(1), mc.cores=1/2/4)` 三次结果**两两不相同**。机制：`CNA.MCMC.R:73`（MCMC 分段）与 `copykat.R:123`（dlm 平滑）用 fork 式 `mclapply`，子进程随机流按 **child 序号**分配，故 `mc.cores` 变 → 每个细胞拿到的流变 → **CNV 分段结果变**。`copykat.R:34 set.seed(1234)` 只保证**固定核数下**可复现。**处置**：GP2 全程固定该值并登记；若要改，属**换计算方法**，整批重跑，不得只重跑一部分 |
| **覆盖度地板**（**2026-09-16 由用户签字采纳**） | 负对照细胞 **nFeature ≥ 840**（`P13_Normal` 负对照中位 nFeature） | ⚠️ C | **无文献出处**，来自本项目**负对照稀释实验**（[`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §8）。做法：对纯二倍体 `P13_Normal`（801 细胞）按 f ∈ {1.0…0.2} 二项稀释 UMI 后重跑 copykat，看它**自己**在何处开始判出 `aneuploid`。**取更严的那个地板**：<br>· 预注册主口径 `f_knee = 0.70` → `C*` = **746**（判据：负对照 `rate_not_defined > 0.50`）<br>· **实测假阳性地板 `f = 0.80` → 840**（该处已有 **134** 个负对照细胞被判 `aneuploid`）<br>**采纳 840 的理由**：低覆盖度会**制造假阳性非整倍体**，而该地板比 `not.defined` 地板**更早触发**；取更严值以压假阳性。序列**非单调**（`0 → 0 → 134 → 84 → 41 → 61 → 13`），照原样记录，不做平滑。<br>⚠️ 边界：**n = 1 样本、无生物学重复**；稀释只模拟测序深度下降，真实低质量样本**只会更糟**，故本值是**下界**。<br>⚠️ **实测代价（2026-09-16 补算，签字时未知）**：按 M1 逐细胞表（`results/01_qc/gse308103_per_cell_qc.csv.gz`）实算，全队列 **648,945 → 376,906**，**弃 272,039 = 41.9%**。丢弃率**在各分期间大体均匀**（AAH 39.6% / AIS 41.8% / LUAD 42.3% / Normal 43.7% / MIA 46.8%），故**不是分期混杂**；但**样本层面极不均**（`P8_Normal` 弃 92.9%，仅留 126 细胞）。**跌破 200 细胞而转 `not_testable` 的仅 1 / 75**（`P8_Normal`）。**代价性质 = 统计功效（丢掉 42% 细胞 + 样本量不均），不是偏倚。**<br>⚠️ **仍需注意的方法学点**：840 源自**单个**样本 `P13_Normal` 的深度；对 31 个自身中位 nFeature < 840 的样本，它比 M1 的**逐样本自适应 MAD**（`nmads=3`）更严，即在本门槛下**部分样本的大部分细胞被排除**。是否改用**逐样本相对**口径，属未决<br>🔴 **实测的两个连带效应（2026-09-16 17:20 补算，签字时均未知）**：地板作用在**进 copykat 之前**，故它同时改了**分母**与**7000 覆写判据**，两个后果都必须重算（已做，`prereg_gene_tiers_floor840.csv`）：
> ⚠️ **2026-09-17 审计标缺（未修，留给 GP2）**：本行及其下两行的**全部实测数字是在旧 648,945 细胞集上算的**。
> 现行分析集已改为 **413,697**（论文 QC，2026-09-16），而 GP2 的输入是 **GP6 定出的上皮细胞**（该子集在 GP6 之前不存在）。
> 故：**地板值 840 本身**（来自 `P13_Normal` 稀释实验）**不受影响**，但**"648,945→376,906 / 弃 41.9% / `not.defined` 0.408→0.097 / 7000 覆写 23/75→2/75" 这些连带效应数字必须在新细胞集与上皮子集上重算**，
> 否则 GP2 会以过期的分母与检出率做判据。**此处如实标缺，不代为重算。**<br>① **7000 覆写 23/75 → 2/75**。机制：剔低覆盖细胞**抬高**了逐基因检出率 → 过 `LOW.DR` 的基因变多 → 判据分母 `n_after_LOWDR_fullgenes`（**全基因**、早于注释，`copykat.R:57`）中位 **7,680 → 9,136** → 7000 这条线落进更空的区域。**注意判据列是 `n_after_LOWDR_fullgenes`，不是 `n_genes_final`**（后者是注释+删周期/HLA 之后的集合，两者可差数千）。<br>② **`not.defined` 率中位 0.408 → 0.097（总量 233,582 → 30,901，÷7.6）**。原因：地板的判据与产生 `not.defined` 的判据**都是覆盖度**（`copykat.R:86-105` ToRemov2、`:194-213` ToRemov3），故地板**预先清掉了**那些细胞，而非仅缩小分母。**我曾预判此处"率会升高"，实测相反 —— 记录在案。**<br>⚠️ 因此本行与下行的"地板只改分母"隐含理解**是错的**，此处更正 |
| **传入 CopyKAT 的细胞集合**（**2026-09-16 新增；当前主口径已被实测否决**） | 🔴 **未定**。现行主口径"整个样本的全部细胞"**实测不可用**；`epi_call()`（仅上皮）**已证既非必要也非充分** | 🔴 | **无文献出处，来自本项目实测**（`results/03_cnv/GP1_report.md` §9.7–§9.9；判据预注册于 [`EPI_PREREG.md`](../results/03_cnv/EPI_PREREG.md)）。**根因**：`norm.cell.names=""`（我们的调用）令 copykat **自己猜**正常细胞基线（`copykat.R:166-170` → `baseline.norm.cl()`：`hclust(ward.D2)` 分 6 群、各群拟合三成分正态混合取 `SDM`、取 `SDM` 最小的簇当正常）。**实测 13 个臂里有 10 次**该步打印 `low confidence in classification`，之后**不拒绝作答**，静默回退 `baseline.GMM(max.normal=5,…)` 重取 `preN` 并把 `WNS` 改写成 `""`，最终标签以确定措辞输出。随后 `copykat.R:456` `if(cor(conses.diploid,conses.aneuploid)>=0.6)` 是**二值**分支：过线 → **全部判 diploid**；不过线 → 按 Wasserstein 把**每个簇**强制二分，**无"不确定"档**。<br>**实测后果（判决轴随输入漂移）**：<br>· 全细胞 → 判决 = **谱系**。`P11_Normal` 判"非整倍体"组 EPCAM 0.67/PTPRC 0.00；判"二倍体"组 LYZ 0.63/DCN 1.02。<br>· 仅上皮（`strict`）→ 判决 = **亚型**。`P11_Normal` 判"非整倍体"组 SFTPC 4.79；判"二倍体"组 AGER 2.88（免疫/基质 marker 全为 0，上皮判据确实生效）。<br>· **同一 `P11_Normal`，两条都说得过去的上皮判据给出 0.0%（argmax）与 32.3%（strict）** —— 收窄细胞集本该只降噪，却把非整倍体从 0% 抬到 32.3%。<br>**已撤回的补救方案**："改报连续 CNV burden、丢掉二值标签"**不足以**解决——基线猜错时整张 `final_results_bin_by_cell.txt` 都是相对错误基线算的。<br>**待验证的出路（"锚定 + 同质"须同时满足）**：① **锚定** = 显式给 `norm.cell.names`，不让它猜；② **同质** = 输入须为**同一亚型**的 恶性 + 正常 混合。仅锚定不够（AT1/AT2 并存时，参考之外的那一种照样被判非整倍体）。<br>🔴 **GP2 不得启动**，直至本行定稿并签字 |
| ↳ 上一行的一条已定、一条未定（2026-09-17） | **已定**：「传入的细胞集合」中的**上皮那一部分** = 两标准交集 **128,091 核**（用户 2026-09-17 签字，见 §M3-A.5 首行；它的用途就是本行的 **CNV 参考**）。清单落盘：[`results/05_annotation/epiCNV_subset_barcodes.txt`](../results/05_annotation/epiCNV_subset_barcodes.txt)，sha256 `db5368c22e1a7ec5…e36d9784`。**仍未定**：「锚定」= `norm.cell.names` 是否显式给、给谁；「同质」= 输入须为**同一亚型**的恶性+正常混合，本数据是 AT1/AT2/Club… 混装 ⇒ 该条**尚无解**。⇒ 本行**仍未定稿**，🔴 **GP2 继续禁跑** | 🔴 | 见上行引文与 §M3-A.5 |
| 退化样本处置 | `n_retained<200` → `not_testable`（**绝不填补为 Normal**）；`n_cnv_pos==0` → 作为真实结果上报，**不得放宽 `KS.cut` 重跑**；Normal/AAH `frac_cnv_pos>0.95` 或 IAC `<0.01` → `implausible` | ⚠️ C | 无文献阈值；本项目约定（防事后调参，法则 3.2） |
| **弃跑判据** | `not_testable` + `implausible` 合计 **>30%** → 停并上报 | ⚠️ C | 本项目约定；届时 R2 对该数据集不可满足，**禁用泛上皮 argmax 替代**。<br>🔴 **2026-09-16 实测评：本判据抓不住 §9.7 的失效模式。** 它的设计对象是"全部判成非整倍体"（Normal `frac > 0.95`）；而实测的失效是"在 0% 与 59% 之间随配置任意翻转"，`P13_Normal` 的 0.407 **远低于 0.95，不会触发**。本判据**必须改写**为能捕捉"判决对输入/参数的敏感性"（如：同一负对照在相邻配置下跨越某分界），否则形同虚设。改写方案待定，**不得**用当前版本充当门禁 |
| infercnv（**可选**交叉验证） | `cutoff=1, window_length=101, HMM=FALSE, HMM_type='i6', denoise=FALSE, cluster_by_groups=TRUE` | ✅ S | `infercnv/R/inferCNV_ops.R:242` — **注：v1 误写 `denoise=TRUE`；默认实为 `FALSE`** |
| 细胞数上限 | infercnv **仅用于 5–10k 细胞/样本子集** | ⚠️ C | 官方已声明 *no longer supported*；全量 200k+ 不可行 |
| 🔴 **源论文的恶性判定方法（对照；2026-09-17 核查）** | Peng 2026 *Cancer Cell* **在 snRNA 端不用 CNV 判恶性**。它的做法是：① spot 层面真值 = **病理学家在 H&E 上标注**（"含 ≥1 个肿瘤细胞即 tumor"，"normal epithelial spot 不应含肿瘤细胞"）；② CNV = **SpatialInferCNV（Erickson 2022）跑 Visium ST 的*上皮 spot***（原文："we applied SpatialInferCNV to infer genome-wide copy number alterations (CNAs) in epithelial spots from Visium ST data"）；③ 参考 spot：先 "without reference" 模式聚类分 clone，再取 CNA 信号最低者。**PMC 全文中 CopyKAT 一词从未出现**。<br>⇒ **本项目 R2 用 CopyKAT 跑 snRNA 是我们自己选的路线，源论文没这么做**；GP1 实测到的失效（判决跟踪谱系、非 CNV）在源论文里是被**换成空间模态**绕开的，不是被解决的 | ✅ S | 我本人核过 [PMC12980502](https://pmc.ncbi.nlm.nih.gov/articles/PMC12980502/) 原文（Results "Spatial landscape of copy number alterations and clonal architectures…" + STAR Methods）。**未核实项**：`cutoff=0.1 / cluster_by_groups=TRUE / HMM=TRUE / denoise=TRUE` 与参考 spot 选择的完整步骤，仅来自代理检索、我未能从原文确认，**标为未核实** |
| 🔴 **输入数据可用性（决定可选工具；2026-09-17 实测）** | `GSE308103_RAW.tar`（1.56 GB）解出 **75 个 `*.raw_counts.mtx.txt.gz`**，**别无他物**。**无 BAM、无 FASTQ、无 CRAM**。`HRA001130` 为受控库、**已移出范围**（见 [`00_ingest/README.md`](../00_ingest/README.md):29）。⇒ **凡需 BAM 做等位 pileup 的工具一律不可用** | ✅ S | 目录实测 2026-09-17 |
| ↳ 因此被排除的工具（**等位基因法，全部须 BAM**） | **Numbat**（Gao T, et al. *Nat Biotechnol* 2023;41:417-426, PMID 36163550）须 `pileup_and_phase.R --bams` ＋ 1000G VCF 定相；**HoneyBADGER**（Fan, *Genome Res* 2018;28:1217, PMID 29898899）须 BAM；**CaSpER**（Serin Harmanci, *Nat Commun* 2020;11:89, PMID 31900397）须 BAM。三者**免正常参考**的优点于是**用不上** | ⚠️ C | 出处如上；Numbat 输入要求核自其 README/vignette |
| ↳ 可用的只剩**表达矩阵法**（均为单细胞级） | **CopyKAT**（Gao R, et al. *Nat Biotechnol* 2021;39:599-608, PMID 33462507）；**SCEVAN**（De Falco, *Nat Commun* 2023;14:1074, PMID 36841879）；inferCNV（Patel 2014 *Science* 344:1396，粒度粗、须参考）。**三者都须显式给正常基线**才可靠 | ✅ S | 工具对照见下 |
| ↳ **SCEVAN 作为第二条腿**（2026-09-17 核查） | **只要表达矩阵**（README："Count matrix with genes on rows … and cells on columns"）、**不要 BAM、不要配对 DNA**、自称"completely unsupervised"；且**支持显式传入** `norm_cells`（"Vector of possible known normal cells to be used as confident normal cells (optional)"）＋ `FIXED_NORMAL_CELLS`。⇒ 与 CopyKAT **同一输入前提、不同实现**，可做互证的第二条腿 | ✅ S | 核自 SCEVAN GitHub README（AntonioDeFalco/SCEVAN）2026-09-17 |
| 🔴 **文献空白（必须如实登记）** | **没有任何论文讨论过"AAH/AIS 等癌前病变中单细胞 CNV 恶性判定是否可靠"。** CopyKAT 原文的前提是"stromal cells have diploid genomes"——**在没有真肿瘤时该前提崩塌**。⇒ 本项目对 **AAH/AIS 的恶性判定没有任何文献背书**，只能作**低置信度结果**上报，**不得**据此宣称"早期病变存在恶性克隆"（法则 4） | ⚠️ C | 检索记录 2026-09-17；**"未找到"不等于"已证伪"，也不等于"可放心用"** |
| ↳ 工具 benchmark（**仅供参考，非同行评审**） | 一项 bioRxiv 预印本（2024.09.26.615284）称 Numbat 多数指标最优、**仅表达矩阵时推荐 CopyKAT**。🔴 **该 benchmark 为 preprint，未同行评审，不得作为判据**，只作参考 | 🔴 | 预印本，**标为不可作依据** |
| ✅ **路线决策（2026-09-17 用户签字）** | **在单细胞层面解决，不走空间模态**。**双臂互证**：CopyKAT ＋ SCEVAN，喂**同一个显式锚定**（CopyKAT 走 `norm.cell.names`，SCEVAN 走 `norm_cells`），报两者一致性。**锚定来源** = **同患者 Normal 样本的同一亚型**细胞 ⇒ **GP2 依赖 GP8a 的产物**，与既有门序"GP2 必须晚于 GP8a"一致。**AAH/AIS 单独标低置信度**并写明原因 | ✅ S | 用户 2026-09-17「走这个思路很好，我们做的就是算法的边界，有问题我们就描述限制」 |
| ↳ 为什么必须"同亚型"（**本项目推断，不是文献结论**） | GP1 实测机制：输入混装 AT1/AT2 时，锚在 AT1 基线上会把 AT2 判成"非整倍体"——**两者的真实表达差异被当成了 CNV**。⇒ 输入须为**单一亚型**的恶性＋正常混装。🔴 **文献里没有写过这条**，属本项目从实测机制推出的要求，**如实标注为推断** | ⚠️ C | 由 [`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §9.7 机制推出 |
| ↳ 双臂不一致时怎么办（**预注册，跑前写死**） | 若 CopyKAT 与 SCEVAN 在同一锚定下**一致性低** ⇒ 结论即"**该数据条件下单细胞 CNV 判定不可靠**"，**照实报告**。**不得**择优取一、**不得**调参凑一致（法则 3.2 禁事后调参） | ⚠️ C | 预注册 2026-09-17 |
| ↳ 连带影响：弃跑判据须重写 | 原弃跑判据（`not_testable`+`implausible` >30%）**已实测抓不住**失效模式（见下）。新路线下它须改为能捕捉"**判决对锚定/输入配置的敏感性**"（如：同一 Normal 负对照在相邻锚定配置下跨越分界）。**改写方案待定，得签字** | 🔴 | 与下行"弃跑判据"联动 |

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

> 🔴 **2026-09-16 全节重写（用户决策「对齐论文全套」）**：本节原为 **scran 池化 / HVG-2000 / PCA-30 /
> kNN-15 / Leiden** 的**本项目自选组合**，**已全部作废**，由源论文的 Seurat 配方（见 §M3-A.0）取代。
> 旧组合移入本节末「已作废登记」以保留决策史。**凡未在下表列出的步骤，一律走工具默认，不得另设。**

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 归一化 + 缩放 | `SCTransform()`（Seurat），**`vst.flavor="v2"`** | ✅ S/P | 论文：*"Data normalization and scale transformation were performed using method 'SCTransform' in Seurat"*。Seurat v5.1.0 的 SCTransform **默认即 v2**，与本机可传值一致 |
| ↳ vst 正则化 | `do_regularize=TRUE`, `theta_regularization="od_factor"`, `residual_type="pearson"` | ✅ S | `sctransform 0.3.5` 源码现取（形参默认） |
| ↳ **v2 flavor 强制项** | `method="glmGamPoi_offset"`、`exclude_poisson=TRUE`、`min_variance="umi_median"`、`n_cells=2000` | ✅ S | `sctransform::vst` 源码：`vst.flavor=="v2"` 分支**无条件覆写**这四项。**不可单独改回**，否则不再是 v2 |
| ↳ θ 估计后端 | `glmGamPoi 1.10.0`（**已装**） | ✅ S | v2 flavor 会**硬依赖**它；缺包时 `stop()`，**不会静默回退**（源码 `if (method %in% c("glmGamPoi","glmGamPoi_offset")) …stop(…)`） |
| ↳ 参数估计子采样 | `ncells=5000`（Seurat 默认；仅用于拟合 θ～μ 曲线，不影响全部基因的残差） | ✅ S | Seurat 4.3.0 `SCTransform` 形参默认 |
| ↳ 残差裁剪 | `clip.range` = Seurat 默认 `±sqrt(ncol/30)` ⇒ 本数据 N=413,697 时 **±117.43** | ✅ S | 同上。⚠️ **不是** scanpy 惯例的 `max_value=10`；两者不可互换 |
| ↳ 回归 | `vars.to.regress=NULL`（**不回归** nCount / pct_mt） | ✅ S | Seurat 默认；也是 v2 flavor 的设计前提 |
| ↳ 随机种子 | `seed.use=1448145` | ✅ S | Seurat 4.3.0 默认 |
| **HVG** | **`variable.features.n=3000`** | ✅ S/P | 论文：*"Top 3,000 HVGs were selected"*；**恰等于 Seurat 的默认值** |
| ↳ HVG 选法 | SCTransform 内部按 `variable.features.rv.th=1.3` + 残差方差选，**不是** scanpy 的 `seurat_v3` | ✅ S | Seurat 默认。⚠️ 与旧配方的 `flavor='seurat_v3'` **不是同一个算法**，不可互相替代 |
| ↳ HVG 顺序语义 | SCTransform **直接在原始计数上**建 Pearson 残差模型 ⇒ **不存在**旧配方"必须先于归一化算 HVG"的顺序约束 | ✅ S | 机制差异，照搬旧注释即为误述 |
| **PCA** | **`npcs=50`**，`seed.use=42`，输入 = SCT `scale.data`（3000 HVG × 全部核） | ✅ S/P | 论文：*"The top 50 PCs were used to calculate the embedding"*；**50 亦为 Seurat `RunPCA` 默认** |
| ↳ SCT `scale.data` 内容 | **`return.only.var.genes=TRUE`** ⇒ `scale.data` **只含 3000 HVG**（非全 18,047 基因） | ✅ S | Seurat `SCTransform` 形参默认。⚠️ 决定 `RunPCA` 的实际输入宽度；**2026-09-17 审计补登记**（此前只在脚本与 manifest，未入表） |
| **批次校正（主口径）** | **Harmony**：`RunHarmony(group.by.vars="sample_id", reduction.use="pca", dims.use=1:50)`，其余**全部默认** | ✅ S/P | 论文：*"**Harmony (version 1.2.0)** was run with **default parameters** to remove batch effects present in the top 50 PCA space"* |
| kNN 图 | `FindNeighbors(dims=1:50, k.param=20, reduction="harmony", annoy.metric="euclidean", nn.method="annoy", prune.SNN=1/15, graph.name=c("harmony_nn","harmony_snn"))` | ✅ S/P | 论文：*"'FindNeighbors' function of Seurat was used to construct the SNN Graph"*。⚠️ **`dims` 默认 `1:10`**（不是 50），**必须显式改 1:50**；`graph.name` 默认 `NULL` ⇒ 内部推成 `paste0(reduction, c("_nn","_snn"))`，**显式写出以免改名后失联** |
| 聚类 | `FindClusters(resolution ∈ 0.5–0.8, algorithm=1, random.seed=0, modularity.fxn=1, n.start=10, n.iter=10, group.singletons=TRUE)` | ✅ S/P | 论文：*"'FindClusters' function with a resolution set to 0.5-0.8"*；`algorithm=1` 即 **Louvain**，为 Seurat 默认，与论文 Xenium 节 *"Louvain algorithm"* 字面一致 |
| UMAP | `RunUMAP(dims=1:50, reduction="harmony", n.neighbors=30L, min.dist=0.3, metric="cosine", umap.method="uwot", spread=1, n.components=2L, seed.use=42L)` | ✅ S/P | 论文：*"2-D visualization … using UMAP with Seurat function 'RunUMAP'"*。⚠️ **三处默认陷阱**：① `dims` 默认 `NULL`（→全部 50，**恰好**等于论文值，仍显式写）；② **`reduction` 默认是 `"pca"`** —— 留空则主口径的图会在**未校正空间**上画，属静默默认，**必须显式改 `"harmony"`**；③ **`metric` 默认 `"cosine"`**（非 scanpy 的 `"euclidean"`），论文未指定 ⇒ 取 Seurat 默认并**如实登记**，不假装论文说过 |
| 簇 marker | `FindAllMarkers(logfc.threshold=0.25, min.pct=0.1, only.pos=FALSE, random.seed=1)`；论文另加**在 ≥10% 组内细胞中表达**的基因才保留 | ✅ S/P | 论文：*"Genes expressed in at least 10% of each group were kept for further analysis"* |
| 环境 | Seurat **4.3.0** / SeuratObject 4.1.3 / sctransform **0.3.5** / glmGamPoi **1.10.0** / harmony **2.0.5** / future 1.29.0（R 4.2.2 Patched） | ✅ S | 2026-09-16 实测 |
| ✅ **降维参数与原文逐条复核（2026-09-17，用户质询「确认和原文一样吧…会不会画出来太稀疏或太密」）** | **复核方法**：`efetch` 取全文 XML 落盘，抽出 **`snRNA-seq data analysis`** 一节逐字比对（**注意**：该文另有 `Visium CytAssist ST data analysis` 一节，两节参数**不同**——ST 节写的是 top **30** PCs、resolution **0.8**；**勿把 ST 节点参数串到 snRNA 上**）。<br>**逐字结论**：snRNA 节原文 *"Top **3,000 HVGs** were selected for PCA and downstream unsupervised clustering. The top **50 PCs** were used to calculate the embedding. **Harmony (version 1.2.0) was run with default parameters** … **top 50 PCA space**. … **'FindNeighbors'** … **'FindClusters'** … **resolution set to 0.5-0.8**. … 2-D visualization … using **UMAP** … **'RunUMAP'"* ⇒ **HVG / PCA / Harmony / 分辨率 四项与原文一致** | ✅ S | PMC12980502 全文 XML，2026-09-17 现取 |
| ↳ 原文**未给**的两个参数（不是"不一致"，是"无从对照"） | ① `FindNeighbors` 的 `k.param`：原文只说"construct the SNN Graph"，**未给 k** ⇒ 用 Seurat 默认 20。② **`RunUMAP` 的 `n.neighbors` / `min.dist` / `spread` / `metric`：原文一个字都没写** ⇒ 全用 Seurat 默认。⚠️ **用户担心的"太稀疏/太密"正由这两个参数决定，而原文无值可对** ⇒ **任何改动都是偏离默认、不是偏离论文**，且**不得**因此声称"和原文一致" | ⚠️ S/C | 同上 |
| ↳ UMAP 疏密的正确认识（重要） | ① **UMAP 是聚完之后才画的，不参与任何聚类** ⇒ 其参数**不影响任何簇、不影响 `r*`、不影响亚型判读**；② UMAP 的点密度**不代表生物学密度**（uwot 会归一化密度）；③ 41.3 万核一张图必然呈密云状。**已实际渲染核对**：`figures/umap_full_run_by_lineage.png`（6 谱系干净分离、无过疏/过密）| ✅ S | 实测渲染，2026-09-17 |
| ↳ 矩阵往返 | 二进制 CSC 三元组 `results/04_integration/seurat_io/`（`02_expression/04_rebuild_expression_paperqc.py` 产出） | ⚠️ C | R 侧 `readBin` + `new("dgCMatrix")`；避开 hdf5r（本机未装） |

**🔴 与论文的版本偏差（三条，全部登记，不假装等价）**

| 项 | 论文 | 本机 | 影响与处置 |
| :--- | :--- | :--- | :--- |
| Seurat | **v5.1.0** | **4.3.0** | v5 的差异主要在 **layer 架构**（`JoinLayers`/`IntegrateLayers`）与 `FindClusters` 的 Leiden 选项；本流程用到的 `SCTransform`/`RunPCA`/`RunHarmony`/`FindNeighbors`/`FindClusters`/`RunUMAP` **接口与语义均相同**。**未逐项证伪差异**，如实标注 |
| sctransform | 随 Seurat 5.1 应为 **0.4.x** | **0.3.5** | **关键点**：`vst.flavor="v2"` 在 0.3.5 中**已完整可用**（已读源码确认强制项），故 v2 正则化的**语义**与论文一致；0.4.x 与 0.3.5 的 v2 实现差异**未逐行比对** |
| Harmony | **v1.2.0** | **2.0.5** | ⚠️ **数值可复现性最弱的一条**。Harmony 2.x 是重写版（新 C++ 核）。**不得声称与 1.2.0 数值一致**。缓解：GP4c 同时跑 **Arm A（不校正）** 作敏感性对照，若两臂结论分歧则升级上报 |

**🟡 已作废登记（旧本项目自选组合，保留以存档决策史，不得再执行）**

| 旧参数 | 旧值 | 作废理由 |
| :--- | :--- | :--- |
| 归一化 | scran 池化 size factor（`quickCluster(min.size=100, method="igraph", graph.fun="walktrap", use.ranks=FALSE, d=NULL, min.mean=0.1, block=NULL)` + `computeSumFactors(sizes=seq(21,101,5), max.cluster.size=3000, positive=TRUE, min.mean=0.1)`） | 被论文的 SCTransform 取代。**scran 的池化在 413,697 核上也从未跑完**（`scran_size_factors.csv.gz` 不存在），故无结果作废问题 |
| 变换 | 归一化计数 `log1p` | SCTransform 用 Pearson 残差，不做 log1p |
| HVG | `n_top_genes=2000`, `flavor='seurat_v3'`, `batch_key='sample_id'`（依赖 `scikit-misc`） | 被 3000 + SCTransform 内部选法取代 |
| 缩放 | `sc.pp.scale(max_value=10)`（仅 HVG 子集） | SCTransform 自带 scale + 自身 clip.range |
| PCA | `n_comps=30` + `{20,30,50}` ARI≥0.90 护栏 | 被固定 `npcs=50` 取代（论文值，也是 Seurat 默认）。**护栏同时作废** |
| kNN | `n_neighbors=15` | 被 Seurat `k.param=20` 取代 |
| 聚类 | Leiden（`leidenalg` 0.10.2，`n_iterations=-1`，`RBConfigurationVertexPartition`） | 被 Seurat `FindClusters(algorithm=1)` = **Louvain** 取代。论文用 Louvain |

> ⚠️ **脚本 `04_integration/00a_export_counts_for_scran.py` / `00b_scran_sizefactors.R` 随之作废**
> （面向旧的 648,945 核集与 scran 路径）；`results/04_integration/scran_io/` 为它们的中间产物，可清理。
> **`04_integration/00c_stage_fragility_check.py` 的结论不作废** —— 该检验已在**新细胞集上重做**，
> 见 `results/01_qc/stage_fragility_report_paperqc.json`。

### M3-A.2 批次校正（**Harmony 为主口径；不校正降为敏感性臂**）

> 🔴 **2026-09-16 改写（用户决策「对齐论文全套」）**：原设计为「Arm A 不校正＝**主** / Arm B Harmony＝对照」。
> 论文把 **Harmony 写进主配方**（*"Harmony (version 1.2.0) was run with default parameters to remove batch
> effects present in the top 50 PCA space"*），故**主口径改为 Harmony**。原 Arm A（不校正）**不删除**——
> 在 `sample_id`↔`stage` 共线的前提下它仍是必要的稳健性对照，但**它不属于论文配方**，
> 其结论只能作为稳健性陈述，**不得当作主结果上报**。

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **主口径** | **Harmony**：`RunHarmony(obj, group.by.vars="sample_id", reduction.use="pca", dims.use=1:50, reduction.save="harmony")`，其余**全部默认** | ✅ S/P | 论文 Methods；`harmony::RunHarmony.Seurat` 形参默认现取（2.0.5） |
| ↳ 校正轴 | **`sample_id`（75 层）** —— 唯一真实技术批次轴 | ✅ S | 本项目设计：每样本＝一病灶切片/一次建库 |
| ↳ **禁止**以 `patient_id` 作 batch | 分期**嵌套于患者内**；按患者校正会抹掉 M4 门判据④所需的**患者内配对对比** | ✅ B | 本项目设计固有约束 |
| ↳ 版本偏差 | 论文 Harmony **1.2.0** vs 本机 **2.0.5**（重写版，新 C++ 核）⇒ **不声称数值等价** | ⚠️ C | 见 §M3-A.1 版本偏差表 |
| ↳ `project.dim=TRUE` | 把 Harmony 嵌入**投影回 Seurat 对象**的 `reductions`/`embeddings` 槽（否则只返回矩阵、后续 `FindNeighbors(reduction="harmony")` 会取不到） | ✅ S | `harmony::RunHarmony.Seurat` 形参默认。**2026-09-17 审计补登记**（此前只在脚本与 manifest，未入表） |
| **敏感性臂（明确不在论文配方内）** | **不校正**：未校正 PCA 空间直接 `FindNeighbors`/`FindClusters`（同 dims/k/resolution） | ⚠️ C | 共线风险下的稳健性对照 |
| ↳ Python 路径作废 | 原登记的 `harmonypy 0.0.10` 路径**作废**——主干已在 R/Seurat 内，不走 Python 复实现 | ⚠️ C | 避免两套实现混用产生不可归因差异 |
| 共线风险 | 每样本＝一病灶＝一分期 ⇒ `sample_id` 与 `stage` **部分共线**；Harmony 可能**一并抹掉分期信号** | ⚠️ C | 本项目设计固有 |
| **两臂上报规则（预注册，不可事后改）** | 两臂各自聚类后算逐细胞 **ARI(Harmony, 不校正)**。若 **ARI < 0.7** ⇒ 两臂**结论分歧**：**不得择一上报**，须**同时呈现**两臂的簇–分期构成表，并附"分期可分性"证据（§M3-A.6），升级给你裁定 | ⚠️ C | 本项目约定 |
| ↳ **裁定结果（用户 2026-09-17 签字）** | 条款已触发（r=0.6 稳健 5/5 种子；r=0.7 边缘 1/5）。**用户裁定：采用源论文的 Harmony 主口径**；不校正臂保留为**敏感性/局限陈述**，**不得当作主结果上报**。⚠️ 这是对分歧的**显式裁定**，非分歧消失 —— 下列限定**对全部下游结果具约束力**：① 逐细胞 ARI 0.6575(r0.6)；② 抹掉约 60% 分期信号（0.5388→0.4866）；③ 嵌套致**不可归因** ⇒ **分期相关结论不得声称已排除批次混淆**，M4 须写入解释边界 | ✅ S | `results/04_integration/seurat_trad/full/GP4_report.md` §九 |
| ↳ **实测：全量 413,697 核（2026-09-17）** | 逐细胞 ARI = r0.5 **0.7072** / r0.6 **0.6575** / r0.7 **0.6960** / r0.8 **0.7251** ⇒ **r=0.6 与 r=0.7 触发升级条款**。用 5 个 Harmony 种子复核触发强度：**r=0.6 为 5/5 种子全部 <0.7（稳健）**，r=0.7 仅 1/5（0.6960，**边缘、种子依赖**，均值 0.7115 在线上） | ✅ S | `results/04_integration/seurat_trad/full/clusters.csv.gz` |
| ↳ **实测：分期可分性（升级条款要求的证据）** | kNN(k=15, 8:2 分层留出) 从嵌入预测分期：**PCA 校正前 0.5388**（多数类基线 0.4546，**+8.42 pp**）→ **Harmony 校正后 0.4866**（**+3.19 pp**）⇒ **Harmony 抹掉约 60% 的分期可恢复信号**；逐分期召回 IAC 0.733 / Normal 0.401 / AIS 0.269 / MIA 0.049 / AAH 0.082。⚠️ `sample_id`↔`stage` 嵌套 ⇒ **"批次"与"分期"不可分，被抹掉的部分无法归因**。**该风险实测非零，已如实升级** | ✅ S | 同上（复算见 GP4_report.md §4.1） |
| ↳ **实测：两臂簇–分期构成** | r=0.7 下 **98.6% 的细胞落在分期混杂簇**；分期纯度 >90% 的簇仅 5 个、合计 **1.39% 细胞**、**全部为 IAC 小簇**，其中多个 **>93% 来自单一测序样本** ⇒ 更像单样本技术产物，非分期生物学。两臂同型（48 vs 49 簇） | ✅ S | 同上（GP4_report.md §4.2） |

### M3-A.3 分辨率选择（**限定在论文区间内**，预注册，不可事后拟合）

> 🔴 **2026-09-16 改写**：原网格 `{0.2 … 2.0}`（8 点）系本项目自选，**已作废**。论文口径为
> **resolution 0.5–0.8**，新网格**限定在该区间内**。区间外的分辨率**只能作行为诊断图**，
> **不得参与 `r*` 选择**。

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 主网格 | `resolution ∈ {0.5, 0.6, 0.7, 0.8}`（＝论文区间 0.5–0.8 含端点，步长 0.1 穷举） | ✅ P | 论文：*"'FindClusters' function with a resolution set to 0.5-0.8"* |
| ↳ 步长归属 | 论文只给**区间**不给步长；0.1 为本项目在区间内的穷举，**非拟合** | ⚠️ C | 区间是论文的，步长是本项目的 |
| 算法 | `FindClusters(algorithm=1)`（Louvain），`modularity.fxn=1`，`n.start=10`，`n.iter=10` | ✅ P | 论文 Xenium 节 *"Louvain algorithm"*；Seurat 默认 |
| 种子 | `random.seed ∈ {0, 1, 2, 3, 4}`（含 Seurat 默认 0，共 5 个） | ⚠️ C | 稳定性度量要求 |
| ↳ Louvain 的确定性 | 固定 `random.seed` **单次即确定**；5 种子**用于度量稳定性，不取平均** | ✅ S | Seurat 4.3.0 源码 `set.seed(random.seed)` |
| 指标1 跨种子稳定性 | 5 种子两两 ARI 均值，**硬约束 ≥ 0.90** | ⚠️ C | Patterson-Cross 2021, *Cell Rep Methods*（ClustAssess 思路） |
| ↳ **实测：护栏随规模而变化（重要）** | **3,000 细胞**：5 种子给出**同一划分**（仅簇编号不同），两两 ARI **恒 = 1.0** ⇒ 该规模下**空约束**。**60,000 细胞**：ARI 降到 **0.9346–0.9834**（r=0.6 最低）⇒ **恢复承重**。3,000 细胞那次已用 `mclust::adjustedRandIndex` 对 40 组真实标签 + 4 组方向性对照复核，**0 处不符** ⇒ 当时恒 1.0 非本项目实现错误 | ✅ S | `results/04_integration/seurat_trad/smoke_{3000,60000}/resolution_metrics.csv` |
| ↳ **实测：全量 413,697 核（2026-09-17，判定用）** | ARI 均值 = r0.5 **0.8963** / r0.6 **0.9230** / r0.7 **0.9141** / r0.8 **0.8999** ⇒ **r=0.5 与 r=0.8 未达 0.90，被该硬约束剔除** ⇒ **护栏在全量规模上真实承重**（非空约束）。**已终结"能否用小规模推断"之争：不能。** | ✅ S | `results/04_integration/seurat_trad/full/resolution_metrics.csv` |
| ↳ 该实测的**处置** | **不得**引用 3,000 细胞那次去断言"跨种子稳定"；本护栏**只在全量 413,697 上判定**。全量若仍出现 ARI 恒 = 1.0，则在 GP5 报告里明写"指标1 空、不承重"，改由**指标2（跨分辨率）+ 指标3（谱系覆盖）+ 指标4（AAH 吸收）**承担。**全量实测非空且承重 ⇒ 上述退路不启用。** | ⚠️ C | 本项目约定 |
| 指标2 跨分辨率稳定性 | 相邻网格 ARI（找稳定平台） | ⚠️ C | 同上 |
| 指标3 谱系覆盖 | 有 ≥1 个法则2 marker 集在 **≥25% 细胞**检出且模块分均值 >0 的簇占比，**硬约束 ≥ 0.90** | ⚠️ C | Heumos 2023（marker 评估） |
| **指标4 · AAH 吸收护栏（新增，可证伪）** | 先固定定义：簇 c 中 Normal 细胞占比**最高且 >50%** ⇒ c 为 **Normal 主簇**。`absorption_rate(r) = AAH 核落在 Normal 主簇中的比例`。**硬约束 `absorption_rate ≤ 0.50`** | ⚠️ C | 本项目约定，**跑之前固定** |
| ↳ 指标4 的依据 | AAH 仅 **9 样本 / 8 患者**，且计数层面**已实证不脆弱**（`stage_fragility_report_paperqc.json`：中位 nFeature/Normal = 1.081×，最差门槛存活差 +3.43 pp）⇒ 若 AAH 确有脆弱性，**只可能在被聚类吸收** | ✅ S | `results/01_qc/stage_fragility_report_paperqc.json` |
| ↳ 违反指标4 的处置 | 该分辨率**不得选为 `r*`**；若**区间内全部**分辨率 `absorption_rate > 0.50` ⇒ **停在 GP5，不放宽阈值、不选 `r*`**，升级给你 | ⚠️ C | 本项目约定 |
| ↳ **实测：全量 413,697 核上该护栏空转（2026-09-17）** | 4 个分辨率下**没有任何簇满足"Normal 占比最高且 >50%"** ⇒ **分母为空** ⇒ `absorption_rate` 按定义恒为 **0**。**该护栏未提供任何信息，不得计入 `r*` 的通过理由。** | ✅ S | `results/04_integration/seurat_trad/full/resolution_metrics.csv` |
| ↳ **裁定（用户 2026-09-17 签字）** | **承认其在全量规模上无效**。⇒ `r*` 实际由**指标1（跨种子 ARI）+ 指标2（跨分辨率稳定性）**承担；指标3 待 GP6 后补算；指标4 不计入。护栏**保留登记但标记"全量无效"**，不删除、不重设计。⚠️ 将来若在别的数据/粒度（含 §M3-A.5 的**子集聚类**）重启，**须重新预注册**，**不得**把"全量空转"当作子集上也必然空转的理由 | ✅ S | `results/04_integration/seurat_trad/full/GP4_report.md` §九 |
| ↳ 🔴 **实测：全量（2026-09-17）本护栏空转** | 4 个分辨率下**没有任何簇满足"Normal 占比最高且 >50%"**（Normal 主导的簇存在，但占比最高仅约 31%）⇒ **分母为空，`absorption_rate` 按定义恒为 0**。**该护栏本轮未测到任何东西，不得计入 `r*` 的通过理由。** 若要使其承重，须**重新设计该护栏定义**（属改判据，须事前批准并重登记），不得事后调阈值 | ✅ S | `results/04_integration/seurat_trad/full/resolution_metrics.csv` = 0 |
| 选择 | `r* = argmax_r [0.5·ARI_seed + 0.5·ARI_xres]`，受指标 1 / 3 / 4 **三条硬约束** | ⚠️ C | — |
| **破平（显式）** | 两分辨率分差 **< 0.01** ⇒ 取**较低**者（少簇、少过度切分） | ⚠️ C | — |
| ↳ 🔴 **实测：本轮触发破平（2026-09-17）** | 合格候选 r=0.7(score 0.9123) 与 r=0.6(score 0.9041)，分差 **0.0082 < 0.01** ⇒ **`r* = 0.6`**。⚠️ 脚本 `10_seurat_traditional.R` **原第 305 行为裸 `which.max`，行内注释却声称已实现本规则 —— 注释与代码不符**，故全量 run 误报 0.7。**已修为显式破平实现**；产物同步更正（manifest 内 `rstar_correction` 块留痕），**原始指标未动** | ✅ S | 同上 |
| 无解 | **停在检查点，不放宽阈值** | ⚠️ C | — |
| 仅作诊断 | `resolution ∈ {0.2, 1.0, 1.2, 1.5, 2.0}`（论文区间外）：只出诊断图，**不参与 `r*`** | ⚠️ C | 区间外 |

### M3-A.4 注释（双标准交叉）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 标准 A · 基因来源 | 法则 2 规范 marker（6 谱系，80 个阳性）—— **2026-09-17 重建为带引用版**，见下行 | ✅ S | [`05_annotation/marker_panel.py`](../05_annotation/marker_panel.py) |
| ↳ 标准 A 面板出处（一次文献） | 上皮 Travaglini 2020 *Nature* 587:619 ＋ Vieira Braga 2019 *Nat Med* 25:1153 ＋ **彭 2026 *Cancer Cell*（源论文）**；T/NK 同前两者 ＋ Guo 2018 *Nat Med* 24:978；B/浆 同前两者；髓系 Travaglini 2020 ＋ Habermann 2020 *Sci Adv* 6:eaba1972 ＋ Zilionis 2019 *Immunity* 50:1317；成纤维 ＋ Reyfman 2019 *AJRCCM* 199:1517 ＋ Lambrechts 2018 *Nat Med* 24:1277；内皮 ＋ Gillich 2020 *Nature* 586:785 ＋ Lambrechts 2018 | ✅ S | 逐基因出处见 `marker_panel.py::REFERENCES` |
| ↳ 标准 A 与标准 B 的独立性 | 面板**刻意不取自 HLCA 整合图谱本身**（否则 κ 自证）。**残留非独立**：HLCA 整合了上述一次研究数据 ⇒ 只保证「定义来源不同」，**不保证「底层数据不同」** | ⚠️ C | 如实声明，**不声称统计独立** |
| ↳ 面板基因缺失 | 4/80 不在本对象基因轴：`KRT18` `SFTPA2` `CD8B` `FCN1`（被 `gene_min_cells=3` 滤除）⇒ 运行时按谱系剔除并上报，**不填补** | ⚠️ C | 实测 |
| 标准 A · 方法 | `sc.tl.score_genes`（`ctrl_size=50`, `random_state=0`）模块分 → 簇均值 argmax | 🟡 D | scanpy |
| ↳ 标准 A 的簇级口径 | 逐细胞标签 = **该细胞所属簇**的 argmax 谱系（GP6 规格：簇均值 argmax） | ⚠️ C | — |
| ↳ 代表种子 | `r*=0.6` 的 5 个种子**全部计算并全部上报**；代表种子 = **与其余 4 个种子平均 ARI 最高者**（最代表该分辨率共识分群）。实算 = **`seed0`，mean ARI 0.93275**。**残留问题如实登记**：GP5 原表只给跨种子均值、**未指定代表种子**，本项目是在 GP6 期间才补定此规则的 ⇒ 规则与被冻结的 `seed0` 恰好一致属**事后确认**，不是预注册。但该确认有实算支撑（5 个种子的 mean-ARI 排序），且**不改变任何已算结果**（无需重跑） | ⚠️ C | 本项目约定；补登记 2026-09-17 |
| 标准 B · 基因来源 | **独立图谱**：`celltypist` **`Human_Lung_Atlas.pkl`** | 🟡 P | 🔴 **2026-09-17 更正**：模型自带元数据为 **Sikkema 2023 *Nat Med* 29:1563（HLCA 整合图谱，doi:10.1038/s41591-023-02327-2, v2, 61 类）**——本行原登记为 Travaglini 2020，**系错误**，已改 |
| ↳ 模型特征覆盖 | 模型 5,017 特征中本项目含 **4,500（89.7%）**；缺 517 个由 CellTypist 内部处理，如实登记 | ⚠️ C | 实测 |
| 标准 B · 方法 | `CellTypist.annotate`（含 `majority_voting`，`over_clustering`=seed0 簇）**且** `rank_genes_groups('wilcoxon')` top-50 × 图谱集 Jaccard | 🟡 D/S | 两法互不叠加，正交 |
| ↳ 61 类 → 6 谱系映射 | **显式登记表**（非自动推断）；7 类判为**判断项**单列上报：`Hematopoietic stem cells` `Mast cells` `Pericytes` `SM activated stress response` `Smooth muscle` `Smooth muscle FAM83D+`；`Mesothelium` 保持**未归属** | ⚠️ C | `marker_panel.py::CELLTYPIST_TO_LINEAGE` / `CELLTYPIST_AMBIGUOUS` |
| ↳ Kappa 双口径 | 主 κ **含**判断项（阳性口径）；另报**剔除判断项**的 κ，两者并报，防单口径掩盖 | ⚠️ C | — |
| 一致性 | 逐细胞 **Cohen's κ** + 逐簇一致率 | 🟡 P | Landis & Koch 1977 释义；**切点本项目定** |
| **切点（预注册）** | `κ≥0.80 且 ≥90%` 通过；`0.60≤κ<0.80` 标记+**强制人工复核**；`κ<0.60` **停** | ⚠️ C | — |
| 分歧 | 全部写入 `annotation_disagreement.csv`，**绝不自动裁决** | ⚠️ C | — |
| 恶性身份 | **仍以 CNV 为准**（法则 0.3），与两标准说什么无关 | ✅ B | 法则 2 |
| **GP6 结果** | `r*=0.6` / 冻结 `seed0` / 45 簇：**κ = 0.8539 → PASS**；逐细胞分歧 47,636 / 413,697（11.5%）。产物 `results/05_annotation/gp6_metrics.json`、`gp6_cluster_labels.csv`、[`GP6_report.md`](../results/05_annotation/GP6_report.md) | ✅ S | 实算 2026-09-17 |
| 🔴 **B2 腿判定为不可用（重要）** | 逐簇 argmax 两两一致：**A vs B1 = 41/45**、A vs B2 = 7/45、B1 vs B2 = 8/45、三方一致 7/45。**在 38 个"三方不一致"簇里 A 与 B1 仍一致的有 34 个** ⇒ 落后的只有 B2。B1 逐簇一致率中位 **0.882**（健康）；**B2 的 Jaccard 全域仅 0.0086–0.0821，中位 0.0446**，此区间内 argmax 无区分度。⇒ **`tri_agree=7/45` 是 B2 失效的证据，不是 38 个簇标错的证据**；主判定以 A↔B1 的 κ 为准 | ✅ S | 实算 2026-09-17；**不重调 B2**（法则 3.2 禁事后调参），降级为记录性指标 |
| ↳ 面板缺基因的根因（2026-09-17 追查） | `panel_missing` = 上皮 `[KRT18, SFTPA2]`、T/NK `[CD8B]`、髓系 `[FCN1]`。**根因在数据源**：作者上传的 GEO 计数矩阵**本身只有 18,082 个基因**，上述 4 基因在原始 `.raw_counts.mtx.txt.gz` 中**不存在**（逐个 `grep -cx` = 0）。**非本项目管道丢弃** | ✅ S | 原始文件实测 |
| 网络 | CellTypist 模型下载须走 `--noproxy` / `ProxyHandler({})` | ⚠️ C | 死代理绕过 |

### M3-A.5 亚聚类（全部主要谱系）

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 谱系 | 上皮 / T·NK / B·浆 / 髓系 / 成纤维 / 内皮 | ✅ S | 法则 2 表 |
| **上皮子集口径（用途：CNV 参考）** | **两标准交集**：标准 A（法则2 marker）判上皮 **且** 标准 B（CellTypist）判上皮 → **128,091 核**（75 样本 / 23 患者；IAC 51,802 · Normal 31,903 · AIS 28,025 · AAH 11,821 · MIA 4,540）。**2026-09-17 用户签字**。**用途限定（重要）**：本口径服务于 **CNV 参考**（决定 CNV 在哪些细胞上做），**不自动继承为亚聚类（GP8a）的输入口径** | ✅ S | 本项目约定。备选未采纳：标准 A 单口径 141,105 / 标准 B 单口径 143,034 / 并集 156,048 |
| ↳ **GP8a（上皮亚聚类）的输入口径** | ✅ **标准 A 单口径**（本项目 marker 面板，法则 2）⇒ `r*=0.6` / 冻结 `seed0` 下**上皮 141,105 核**（75 样本 / 23 患者；IAC 58,750 · Normal 34,344 · AIS 30,306 · AAH 12,815 · MIA 4,890）。种子敏感性：5 个种子 141,105–141,690（极差 0.41%），**不随种子翻盘**。<br>清单落盘：[`results/05_annotation/epiA_subset_barcodes.txt`](../results/05_annotation/epiA_subset_barcodes.txt)，sha256 `a952857e3dbb119d…b0324d34`（141,105 行） | ⚠️ C | **签字范围须分清**：用户 2026-09-17 签字的是**方向**——「按照传统那一套做下去，CellTypist 只用来测试质量」（即单口径、标准 B 退为质检）。**141,105 这个数字与清单是本项目据该方向计算的产物，非用户逐字确认**。此处如实分层，不让"选了方向"冒充"确认了数字" |
| ↳ **GP8a 启动前的硬前置** | 用户 2026-09-17 指令：**亚群重聚类与 CNV 开跑前，须先确认「具体 marker 出处和参数」**。⇒ 该前置**已于 2026-09-17 满足**：① marker 出处查清并修毕（§M3-A.5c 引用行 + §M3-A.5d 审计）；② 标签口径经用户裁定为「先聚类、再按簇注释」（§M3-A.5c）。**GP8a 于 2026-09-17 21:10 启动** | ✅ S | 用户 2026-09-17 原话；两项前置均已有裁定记录 |
| 🔴 **范围更正：6 个谱系全做（用户 2026-09-17）** | 用户原话：「**几个谱系都需要重聚类，不仅仅是上皮**」。⇒ 亚聚类**不是只做上皮**，§M3-A.5 登记的**全部 6 个谱系**都要跑。本节此后按 6 谱系执行 | ✅ S | 用户 2026-09-17 指令 |
| **各谱系输入清单（构建器 `05_annotation/00_build_lineage_subsets.py`）** | 口径同 §M3-A.5（标准 A 单口径，`A_frozen` seed0）。清单经 sha256 落盘，manifest `results/05_annotation/lineage_subsets_manifest.json`。**互斥且完备**：6 谱系合计 **413,697 = 全部细胞**（断言通过）。<br>上皮 141,105（`epiA_subset_barcodes.txt`，sha256 `a952857e…`，**与已签字文件逐字节一致**，脚本会硬校验）；T·NK 62,483（`tnkA…`）；B·浆 34,385（`bplasmaA…`）；髓系 54,407（`myeloidA…`）；成纤维 77,633（`fibroA…`）；内皮 43,684（`endoA…`）。全部 **≥ 子集规模下限 1000** ⇒ 6 个都做 | ✅ S | 实测 2026-09-17 |
| **6 个亚聚类运行登记（2026-09-17 启动）** | 统一命令：`Rscript 04_integration/10_seurat_traditional.R subset --cells <清单> --tag <tag> --lineage <谱系名>`。输出 `results/04_integration/seurat_trad/<tag>/`；日志 `logs/gp8_<tag>_*.log`。**子集内的 SCTransform(v2) / HVG 3000 / PCA 50 / Harmony(`sample_id`) / FindNeighbors(50) / 网格 `{0.5,0.6,0.7,0.8}` × 5 种子 全部逐字同 full**（`10_seurat_traditional.R` 内常量，无分支差异） | ✅ S | 脚本常量现取 |
| ↳ GP8a 的 r* 须单独签字 | 子集内会算出**自己的** `r*`（同 §M3-A.3 规则与破平规则），**不得**沿用全局的 0.6。注释**必须等 r* 签字之后**才做 | ⚠️ C | §M3-A.3 |
| 🔴 **r\* 已算出，6 个谱系待签字（2026-09-17）** | 由 `05_annotation/03_subcluster_annotation.py --stage select` 按 §M3-A.3 规则算出，落盘 `results/05_annotation/<tag>_rstar.json`。**尚未签字，未做任何注释**。<br>**T·NK `0.5`**（0.5 独高 0.985）· **髓系 `0.5`**（0.5/0.6 差 0.008<0.01 并列，取低）· **内皮 `0.5`**（**四个分辨率全部落在 0.01 内**：0.978/0.980/0.980/0.975，取低）· **B·浆 `0.6`**（0.6 独高 0.943）· **上皮 `0.5`**（**仅 0.5 过种子约束**；0.6/0.7/0.8 的 `ari_seed_mean` 为 0.820/0.840/0.848 **全部 <0.90**，故三者不可选）· 🔴 **成纤维：无解** | ⚠️ C | `*_rstar.json`，2026-09-17 |
| 🔴 **成纤维亚聚类「无解」—— 预注册规则触发，不是 bug** | **事实**：四个分辨率**全部**不满足 `ari_seed_mean ≥ 0.90`（0.833 / 0.744 / 0.800 / 0.870）。**诊断**（逐种子两两 ARI，自算）：种子 0 在每个分辨率上都偏离（vs 其余 4 种子均值 0.666–0.873），但**剔掉它仍到不了 0.90**；最好的是 r=0.8（两两均值 0.870、最低 0.840）。**对照组**：内皮 r=0.5 两两均值 **0.983** ⇒ **不是流水线问题，是成纤维这一谱系本身跨种子不稳**。<br>**解释（推断，非结论）**：成纤维亚型（肌成纤维/肺泡/血管周/胸膜下）在文献中本就呈**连续谱**而非离散岛 ⇒ 换种子簇边界即漂。**未找到文献直接支持此推断，不得当证据引用** | 🔴 | 实测 + 自算 ARI，2026-09-17 |
| 🔴 **用户裁定：成纤维降级到 0.85，出 `r*=0.8`，结果只能作探索性（2026-09-17）** | **这是事后放宽预注册阈值**，已按"不得含糊"的要求显式实现与登记：① 代码层 —— `03_subcluster_annotation.py --stage select` 新增 `--relax-seed-min`，**默认不启用**；启用后产物 json 强制带 `relaxed=true` / `seed_min_applied=0.85` / `seed_min_preregistered=0.9` / `relaxed_note`。② 落盘 `results/05_annotation/fibroA_rstar.json` = `r* 0.8`（唯一入选者；其两两均值 0.870、最低 0.840）。③ **约束**：该谱系的亚型注释只能作**探索性**结论，**不得进主结论**；报告中必须标注"跨种子稳定性 <0.90"。④ **回归检查**：不带 `--relax-seed-min` 重跑其余 5 个谱系，`r*` 一字未变（0.5/0.5/0.5/0.6/0.5）⇒ 放宽开关未污染主口径 | 🔴 | 用户 2026-09-17 裁定；`fibroA_rstar.json` 现取 |
| 🔴 **首跑 6 个全部崩在最后一行（2026-09-17，已修）** | **症状**：6 个 run 全部 `Execution halted`，报 `'sha256sum' is not an exported object from 'namespace:tools'`。**根因**：manifest 里写 `tools::sha256sum(...)` —— R 的 `tools` 包**只有 `md5sum`，没有 `sha256sum`**。**影响范围**：崩在**最后一条语句**，故 `clusters.csv.gz` / `umap.csv.gz` / `resolution_metrics.csv` / `timings.csv` / `embeddings/` **全部完整**，**只有 `run_manifest.json` 缺失**。⇒ 科学结果未受影响，**受影响的只有溯源登记** | 🔴→✅ | 6 份日志 `logs/gp8_*_20260917_2112.log` 尾部 |
| ↳ 为什么会漏检（须记住的教训） | `subset` 分支是 2026-09-17 新加的，而参考全量跑 `full/`（2026-09-17 00:43）**远早于**该分支加入（脚本 mtime 19:57）⇒ **`subset` 这条代码路径从未被执行过一次**，首次执行即 6 个 run 全崩。**同类风险**：任何"新加的、未跑过的分支"都可能带着只在运行时才暴露的错。⇒ 新分支**必须先小规模冒烟再全量** | 🔴 | 脚本 mtime 与 `full/run_manifest.json` 时间对比 |
| ↳ 同批修掉的另外两处（同类漏检） | ① `lineage` 被**硬编码成 `"上皮"`** ⇒ 6 个 run 会全部登记成上皮。已改为**必填命令行参数 `--lineage`**，并在脚本内断言其属六谱系之一（写错即停）。② `caliber` 写的是**已作废的**「双标准交集」⇒ 已更正为「标准 A 单口径 = GP6 冻结 `A_frozen`（seed0, r*=0.6）」。**两处都是原样活到崩之前，没有被任何检查挡住** | ✅ S | `10_seurat_traditional.R` 参数解析块与 manifest `subset` 块 |
| ↳ 修法与重跑（2026-09-17 21:38） | `tools::sha256sum` → `digest::digest(file=, algo="sha256")`（已核对：R 侧结果与 `00_build_lineage_subsets.py` 写的 sha256 **逐字符相同**）；脚本开头加 `requireNamespace("digest")` **提前断言**，避免再"跑完 25 分钟才在最后一行报缺包"。**旧产物已留存**于 `results/04_integration/seurat_trad/.prev_nomanifest/`，用于重跑后的**确定性核对**（同输入同种子应产出相同簇） | ✅ S | 实测 2026-09-17 |
| ↳ **确定性核对结果：4/4 逐字节相同** | 重跑后与崩掉那次的产物对比：`tnkA` / `myeloidA` / `endoA` / `bplasmaA` 的 `clusters.csv.gz`、`umap.csv.gz`、`resolution_metrics.csv` **三者全部逐字节相同**；逐细胞比对 `clusters.csv.gz` 的**全部 27 个簇列**，**无一个细胞标签不同**。⇒ **那次崩溃对科学结果的影响是可验证的零**，不只是"认为没影响"。`epiA` / `fibroA` 上次未跑完（中途停下），**无旧产物可比，不声称** | ✅ S | md5sum + 逐细胞比对，2026-09-17 |
| ↳ R 侧 `sha256` 的可用实现（备查） | 本机 `digest` 与 `openssl` 两包**均已安装**；`tools` 包**无**任何 sha 系列函数（`grep("sha", ls("package:tools"))` 返回空） | ✅ S | `requireNamespace` 实测 |
| ↳ **重跑完成（2026-09-17 22:38）** | 6 个全部跑通并写出 `run_manifest.json`，`lineage` 字段逐个正确（上皮/成纤维/T·NK/髓系/内皮/B·浆），`n_actual` 与清单逐一相符，患者均 23。耗时：上皮 3677s（最长，141k 核）> 成纤维 1605s > T·NK 1199s > 髓系 1080s > 内皮 821s > B·浆 633s | ✅ S | 6 份 manifest 现取 |
| 🔴 **PR 前复核 · 又一处「共用分支写死值」（2026-09-17 深夜，已修）** | `metrics_caveat` 的 subset 分支原**写死「细胞集合已限定为上皮交集」**。subset 是六谱系**共用**的代码路径 ⇒ **五个非上皮谱系的 manifest 里都印着这句错话**；且「交集」一词本身已过时（口径已改为标准 A 单口径）。与上一行的 `lineage = "上皮"` **属同一类漏检：共用分支里写死某个谱系**。<br>**修法**：① R 源脚本改为按 `--lineage` 实参生成（并带上清单 basename 与 sha256）；② 已落盘的 6 份 manifest 按 `rstar_correction` 先例**原地更正**，每份加 `metrics_caveat_correction` 块，`original` 字段**完整保留原错文**，未静默丢弃；③ **未重跑**（六谱系约 2 小时）—— 属确定性文字替换，不触及任何随机数或计算路径，`metrics_table` / `clusters.csv.gz` / `umap.csv.gz` / `params` / `versions` / `embeddings` / `wall_sec` **一律未改** | ✅ S | 6 份 `run_manifest.json` 现取 |
| ↳ **subset 分支冒烟（2026-09-17 深夜）** | 因改了 subset 分支，按上一行「新分支必须先冒烟」的教训**真跑一次**：跨 75 样本分层抽 3,000 上皮核，`Rscript … subset --cells … --tag _smoke_subset --lineage 上皮`，**退出码 0**，manifest 文案、`subset.lineage`、清单 sha256 三项均正确。**首次冒烟其实失败了**——我误取清单前 3,000 行（按样本排布，全来自同一样本），Harmony 报 `contrasts can be applied only to factors with 2 or more levels`。⇒ **记一笔**：亚聚类跑全量时样本数必须 ≥2，本批 6 个谱系均含 75 样本，不触发。冒烟目录已删 | ✅ S | `/tmp/smoke_subset.log`，2026-09-17 |
| 🔴 **PR 前复核 · r\* 规则的实现与登记不符（2026-09-17 深夜，已修）** | §M3-A.3 第 244 行登记 r\* 受**三条**硬约束（指标1 跨种子稳定性 / 指标3 谱系覆盖 / 指标4 AAH 吸收），R 侧 `eligible` 列实现的则是**两条可算的**（跨种子 ∧ AAH 吸收）；但 `03_subcluster_annotation.py` 的 `stage_select` 连这**一条**都没用全 —— **只用跨种子稳定性**，静默丢掉了 AAH 护栏。**影响**：本批 6 谱系 `pass_aah` **全为 TRUE** ⇒ 补齐后 **r\* 一个都没变**（0.5/0.5/0.5/0.5/0.6/0.8放宽，已逐条对撞备份确认）。⇒ 属**「实现与登记不符」，不是「结果错」**，但已在 PR 前修正，并新增 `constraints_applied` 字段把实际施加的约束写进 `*_rstar.json` | 🔴→✅ | `03_subcluster_annotation.py`；6 份 `*_rstar.json` 现取 |
| ↳ **破平判据的边界不一致（同批修掉）** | R 用严格 `score > max − 0.01`，Python 用 `score >= max − 0.01`。**恰好在差 0.01 时两边会给出不同的 r\***，且 `resolution_metrics.csv` 的 `is_rstar` 列会与 `*_rstar.json` 互相矛盾。已把 R 统一为 `>=`。**本批数据未落在该边界上**（最近的是 B·浆 0.943 与 0.9323，差 0.0107），故 r\* 不变 | ✅ S | 两处代码现取；边界值逐条复核 |
| ↳ **面板溯源硬校验的漏洞（同批修掉）** | `classic_panels.py::assert_provenance()` 只查 `ALIAS_FIX` 与 `DROP_NON_SYMBOL` 的残留，**漏查 `ORTHOLOG_FIX`**（小鼠基因号 → 人同源）。若日后重生成时 `CYP2F2` 未被替换，原校验会放行。已补上。**同时验证生成器可复现**：重跑 `build_classic_panels.py` 的输出与 `classic_panels.py` **逐字节一致**（除时间戳），且全部替换词均在论文自己的 Table S4 中出现过 | ✅ S | 重跑 diff，2026-09-17 |
| 🔴 **PR 前复核 · 指标3（谱系覆盖）从未参与 r\* 选择（2026-09-17 深夜，登记不改结果）** | §M3-A.3 第 244 行登记 r\* 受**指标 1 / 3 / 4** 三条硬约束；第 242 行已裁定指标3「待 GP6 后补算」、第 241/243 行裁定指标4「全量空转、不计入」。⇒ 本批 6 个谱系的 `r*` **是在只有指标1 承重、指标4 空转、指标2 仅用于破平的条件下选出的**，指标3 **一次都没算过**。**明确后果**：GP8 注释时必须**逐谱系补算指标3**（r\* 处的覆盖率）并报告；若某谱系算得 **< 0.90**，按已登记规则该分辨率**不得作 r\***，须**停在检查点升级**，**不得**事后放宽。**本条不隐藏、不假装三条都生效** | 🔴 | §M3-A.3 第 242/244 行；待 GP8 注释补算 |
| ↳ **「先舍入再判阈」的修法在 `score` 列留下 ≤0.0001 偏离（同批，登记但**不改**产物）** | 旧 R 实现拿**未舍入**的 seed/xres 算 `score` 再舍入；修后改为**先舍入再算** ⇒ 落盘 `resolution_metrics.csv` 的 `score` 列，与「拿该行自己那两列按修后公式重算」的值在 **5/36 行**上差 **0.0001**。逐行为（前者=落盘值 / 后者=按本行列重算）`bplasmaA r=0.8` 0.8781 / 0.878 · `endoA r=0.6` 0.98 / 0.9799 · `fibroA r=0.6` 0.6776 / 0.6777 · `full r=0.5` 0.8787 / 0.8786 · `smoke_60000 r=0.5` 0.9097 / 0.9098。<br>**r\* 一个都没变**：三处落在生产谱系上，`endoA` 手算四者并列（0.9783/0.9799/0.98/0.975）仍取较低者 **0.5**，`bplasmaA` 的 r\*=0.6（0.943）与 `fibroA` 的 r\*=0.8 均不由这两行决定。<br>**决定：不重写已落盘产物** —— 它们如实对应**产出它们的那版代码**，重写只会毁掉可追溯性。**将来重跑以修后版本为准** | ✅ S | 36 行逐行重算比对，2026-09-17 |
| ↳ **确定性的诚实边界（重申）** | 崩溃前对照**只保住了 4 个谱系**（`tnkA`/`myeloidA`/`endoA`/`bplasmaA`，三者逐字节相同）；`epiA`/`fibroA` 当时**尚未跑完**，`.prev_nomanifest/` 下**只有空的 `embeddings/` 目录**（已删）。⇒ 确定性证据是 **4/6，不是 6/6**，不得简称「全部逐字节相同」 | ✅ S | `find .prev_nomanifest -type f` 实测 |
| ↳ **基线 md5 登记（可核验）** | 下列每个文件，「正式产物」与「崩溃前基线」的 md5 **完全相等** —— 这就是「逐字节相同」的证据（值相同故只列一个）。基线本身**不入库**（86 MB，`.prev_nomanifest/` 已写进 `.gitignore`，本地保留备查、可随时删）：<br>`tnkA` clusters `79c46e167eb71a7f1f35007a4fdf5a04` · umap `a3d74163f8c272182178fd4a09a329ad` · metrics `8fe9ecdc7e6c8096ebe2afb03347b0d2`<br>`myeloidA` clusters `ea0ce10abe48079666044b4ec537a8d9` · umap `f83b03ac5031c7150c4bfb4045ff2a45` · metrics `de0f2d4b72da5df177da03704ed18217`<br>`endoA` clusters `a0754562497354c8b9a9e3b94d402a62` · umap `68efccfbbfb8642412767e5fad76e8a3` · metrics `ab845de0081e912d51aaba4f1e9abc6d`<br>`bplasmaA` clusters `4fa9f308648f21e7b414bea57c4979bc` · umap `e400a51d9641229738a6728979eab003` · metrics `517f83f9e14ed5bcd154ce02699a44a4`<br>（每行的 clusters/umap/metrics **新跑值与基线值同**。正式产物已入库，基线未入库 ⇒ 本机可复核，克隆后不可） | ✅ S | md5sum 实测 2026-09-17 |
| 🔴 **簇级注释的具体做法（用户裁定后细化）** | 每个簇：① 报该簇的经典面板**逐型平均得分**（不只报 argmax 冠军）；② 报该簇 `rank_genes_groups('wilcoxon')` top 富集基因；③ 报该簇规模与跨样本/跨患者分布；④ 人工按 ①②③ 判读并**写明理由**。⑤ 与 CellTypist 在同一子集上的簇级众数**并报**。**不自动裁决**，分歧单列入 `epi_subtype_disagreement.csv` | ⚠️ C | Travaglini 2020 Methods 人工判读口径 |
| ↳ 为什么是单口径、不是交集 | **领域惯例是单套标注**：全局聚类 → 注释成谱系 → 取**某一套**标签的子集 → 重聚类。**2026-09-17 文献核查未找到任何一篇**用"两套方法都同意才纳入"来定义亚聚类输入。实例：Sinjab 2021 *Cancer Discov*（10.1158/2159-8290.CD-20-1285，EPCAM± FACS 分选 + 单套 marker，70,030 上皮 → 10 亚型）；Travaglini 2020 *Nature*（10.1038/s41586-020-2922-4，按 compartment marker `SubsetData` 再聚类）；Han/Sinjab 2024 *Nature*（10.1038/s41586-024-07113-9，EPCAM 富集 246,102 上皮核，单一通路）；Sikkema 2023 *Nat Med*（HLCA，直接全细胞分层聚到 level 2–5） | ⚠️ C | 核查记录 2026-09-17；**「未找到先例」本身如实登记，不等于「已证明错」** |
| ↳ 求交集的具体代价（**推断，非文献结论**） | HLCA 自述**过渡态细胞恰是跨研究最缺 consensus 的一类**，且把癌症数据映射到该图谱时 **22% 细胞判 unknown** ⇒ 两标准分歧的细胞很可能正是**过渡态/疾病特异态**。本项目课题是**进展**（Normal→AAH→AIS→MIA→IAC），过渡态正是要找的对象，求交集等于优先丢弃它们。**注意**：A 独有 13,014 ＋ B 独有 14,943 = 27,957 核，正是分歧所在 | ⚠️ C | 🔴 **属推断**——未找到直接研究"求交集"这一做法的文献。**不得**当作已有证据引用 |
| ↳ 标准 B（CellTypist）的角色 | **仅质量校验**：出 κ 与分歧清单，**不参与定义任何输入集合**。κ=0.8539（冻结 seed0）**PASS** ⇒ 标注质量已验收，标准 B 使命完成 | ✅ S | 用户 2026-09-17 指令 |
| ↳ 更正记录（过程留痕） | 2026-09-17 我**曾误**将 CNV 用的交集口径登记为 GP8a 输入并据此启动了一次 run；用户当场更正「交集只是用来做这个 cna 参考的啊」。该 run 已中止、空产物目录已清除（**无结果落盘**）。本行即该更正的留痕，**不静默改掉** | ⚠️ C | — |
| ↳ 为什么 CNV 取交集 | CNV 的**前提**是"这些细胞是上皮"。GP1 已实测：CopyKAT 在本数据上的非整倍体判定**跟踪的是谱系组成、不是 CNV**（见 [`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §9.7）⇒ 前提越松，基线越容易被非上皮细胞带偏。取两套办法**都**认可的最保守口径。被排除的 27,957 核（A 剔 B 13,014 ＋ B 剔 A 14,943）**正是两标准分歧最集中的那部分**，属模糊地带，**不强行拉入** | ⚠️ C | 本项目约定 |
| ↳ 未采纳并集的代价 | 并集多 27,957 核（较交集 +21.8%）；**主动放弃**，换取"上皮"前提的纯度。如实登记，**不声称是最优解** | ⚠️ C | — |
| 协议 | 每个子集**从原始计数重启全流程**（口径同 §M3-A.1，不换算法）：子集内 **SCTransform**（同 `vst.flavor="v2"`，HVG 3000 在**子集内**重选）→ **PCA(50)** → **Harmony(`sample_id`)** → FindNeighbors(50) → 同网格 `{0.5…0.8}` × 5 种子 × **同 M3-A.3 全部指标（含指标4）** → 子集内双标准注释 | ⚠️ C | **本项目约定**——未找到正式证明该协议避免 double-dipping 的论文。最接近的选择性推断文献（Gao, Bien & Witten, *JASA* **2024**;119(545):332–342, doi:10.1080/01621459.2022.2116331）讨论的是**同一份数据先聚类、再检验簇间均值差**的选择性 type I error，**不是**「取子集再聚类」的问题。**不冒充有引用**。<br>🔴 **2026-09-17 出处更正**：本行原写 `2020, JASA 115:1622`（**两个字段都错**，源自 arXiv 预印本年与错误的卷页）。已按 JASA 正式记录改正 |
| ↳ 子集内 HVG 的偏差 | 论文的 3000 HVG 是在**全部 413,697 核**上选的；子集内重选 ⇒ 亚聚类空间**不与论文等价**（论文未做亚聚类，无可对齐对象） | ⚠️ C | 如实登记，**不声称与论文一致** |
| ↳ 子集规模下限 | 子集细胞数 **< 1000** ⇒ 该谱系**不做亚聚类**，只报 L1，如实记录（`SCTransform` 的 `ncells=5000` 与 `min_cells` 类过滤在小子集上不稳） | ⚠️ C | 本项目约定 |
| 层级 | `global_label`(L1) 与 `subcluster_label`(L2) **并存，绝不覆盖 L1** | ⚠️ C | — |

### M3-A.5b 上皮亚型 marker 面板 —— ~~论文 MP 口径~~ **【已作废，勿用】**

> 🔴 **本节已于 2026-09-17 被用户改口径取代，仅存档以备追溯，不得作为 GP8a 的输入。**
> 用户原话：「**就不按什么kac去做，就按照经典marker**」。
> 现行口径见 **[§M3-A.5c](#m3-a5c-上皮亚型-marker-面板--经典-marker-口径现行)**。
> 其中第 311–312 行的「双面板决策」、第 303 行「KAC 第二来源」、第 320–321 行的
> 检出率诊断**均按论文 MP 口径算出，对现行经典口径不适用**，须在 GP8a 重算/重报。
>
> 用户 2026-09-17 指令：亚群重聚类开跑前须**确认具体 marker 出处**。本节即该确认件。
> 执行口径 = [`05_annotation/epi_subtype_panel.py`](../05_annotation/epi_subtype_panel.py)（**程序化生成，非手抄**）。

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **主干来源** | 源论文 **Table S2 第二块**「snRNA-seq (lung epithelium)」的 **9 个 meta-program（MP）**，每 MP **50 基因** | ✅ S | Peng 2026 *Cancer Cell* Table S2。**本项目已把原表冻结入库**：`data/external/peng2026_cancercell/TableS2_meta_program_gene_lists.xlsx`，sha256 `d1c2333ab6031b83…8b2158b0` |
| ↳ 为什么用论文自己的 MP | 论文**定义的亚型就是我们要的亚型**（Ciliated / AT2 / Club-secretory / Basal / AT1 / Tumor-KAC / Tumor-stress / KAC-inflammatory），且 8 个亚型的 **450 个基因槽 100% 命中**本数据矩阵（因二者是同一批数据） | ✅ S | 实测 2026-09-17 |
| ↳ ⚠️ **MP ≠ marker 表（重要限制）** | MP 是 **NMF meta-program**，含应激/管家类共表达基因（如 MP2 AT2 内含 `SOD2` `CXCL2` `NR4A1`）。**整块直接做 `score_genes` 特异性会偏低**，故同时提供 CORE 小集做交叉校验 | ⚠️ C | 如实登记限制 |
| 亚型 → MP 映射 | Ciliated = **MP1 ∪ MP8**（论文给了两个纤毛 MP）；AT2=MP2；Club/secretory=MP3；Basal/basal stem=MP4；AT1=MP5；Tumor cell/KAC=MP6；Tumor cell (stress/inflammatory)=MP7；KAC/inflammatory=MP9 | ✅ S | 原表第 3 行标签 |
| 去重后规模 | **388 个基因**（8 亚型合并去重） | ✅ S | 实测 |
| ↳ 亚型可分性 | 两两 Jaccard **最大 0.136**（MP7 × MP9），多数 <0.05 ⇒ 面板之间**天然分得开**，不是同一堆基因换名字 | ✅ S | 实测 2026-09-17 |
| **KAC 的第二独立来源** | 源论文 **Table S3**「人 KAC」签名 **102 基因**（另有 2 列小鼠，本项目为人类数据，**不采纳**）。冻结于 `data/external/peng2026_cancercell/TableS3_signature_gene_lists.xlsx`，sha256 `043ff6c684dfabd2…546ea7b0` | ✅ S | Peng 2026 *Cancer Cell* Table S3 |
| ↳ ⚠️ **论文两处 KAC 定义只重叠 8 个基因** | MP6（50 基因）∩ Table S3 签名（102 基因，剔矩阵外 11 个）= **8**（`ABCC3, MARCKSL1, MDK, MMP7, S100A11, TIMP1, TMSB10, WFDC2`）；仅 MP6 有 42、仅签名有 94。⇒ **论文自身给的两套 KAC 定义并不一致**，本项目**两套都报**，不挑一个当"真 KAC" | ⚠️ C | 实测 2026-09-17 |
| ↳ Table S3 在矩阵外的 11 个基因 | `GAPDH, PTMA, RPL41, RPL13A, RPL15, TPI1, RPS2, RPS18, RPS19, BCAS1, LINC00511` ⇒ 前 9 个是**核糖体/管家**（探针面板不含，见 M2 节），后 2 个是真实缺失 | ⚠️ C | 实测 |
| 🔴 **Excel 自动改名（2 处，已修复并登记）** | Table S2 的 **MP4** 含 `DKK 3.00`（应为 **DKK3**）、**MP7** 含 `ERN 1.00`（应为 **ERN1**）。生成脚本以 `EXCEL_NAME_FIX` 显式还原，**不静默使用** | ✅ S | 原始 .xlsx 实测，`openpyxl` 日期型单元格检查 = 0，故仅此 2 处 |
| 矩阵外基因 | `SFTPA2`（法则2 上皮 panel 成员）**不在数据中**，面板已剔除并登记于 `MISSING_IN_MATRIX` | ✅ S | 同 M2 节 |
| 🔴 **CORE 小集（经典亚型 marker）出处状态** | 候选集已列出（AT1/AT2/Club/Ciliated/Basal 各 12–14 基因），**逐基因对原论文补充表的核对仍在进行**；**核对完成前 status 一律为 `unverified`，不得进产线** | 🔴 | 调研中（Travaglini 2020 补表 / Vieira Braga 2019 补表） |
| 生成脚本（可复跑） | `python3 05_annotation/build_epi_subtype_panel.py`；产物 sha256 `6ec17e8386cdf5df…a1adf492`，**重跑逐字节一致** | ✅ S | 实测 2026-09-17 |
| ✅ **双面板决策（用户 2026-09-17 签字）** | **两套都算、都报**，不二选一：<br>· **P1 = 论文 MP 全集**（8 亚型，388 基因）—— 覆盖全但含噪声<br>· **P2 = MP ∩ 经典 marker**（更干净、覆盖窄）<br>· 两套各自出亚型标签，**并报两套之间的一致性（ARI + 混淆矩阵）**。<br>**不一致本身就是结果**——它标的正是亚型边界的模糊处 | ✅ S | 用户 2026-09-17 原话「**都做**，这对我们后面的 scmg 迁徙非常重要」 |
| ↳ **签字的理由（用户给的，须一并登记）** | 上皮亚型分辨率**直接喂给下游 SCMG 分支的迁徙分析**；尤其 **KAC**（论文定义为 AT2 与肿瘤细胞之间的中间态）正是迁徙的枢纽细胞 ⇒ 亚型边界画在哪，直接决定迁徙的起点/终点 | ✅ S | 同上 |
| 🔴 **待登记的接口问题** | GP8a 属**传统分支**，而 SCMG 分支受**铁律 R4（不得掺传统算法）**约束。⇒「GP8a 的亚型标签以何种形式进入 SCMG 迁徙」**是 R4 边界问题，须单独立项登记后才能实施**，本行只记录需求，**不预设接口形式** | 🔴 | 本项目约定；R4 |
| `score_genes` 参数（预注册） | `sc.tl.score_genes(gene_list, ctrl_size=50, n_bins=25, random_state=0, gene_pool=None, use_raw=False)` —— 前四项即 scanpy 1.9.8 默认（**已现取签名核实**）；`use_raw=False` ⇒ 在 **X = log 归一化后的表达**上打分 | ⚠️ C | scanpy 1.9.8 签名现取 2026-09-17；⚠️ `random_state=0` **必须显式固定**（该函数按表达分箱随机抽对照基因集，不固定则不可复现，违反 R5） |
| ↳ `gene_pool=None` 的后果（登记） | 对照基因从**全部 18,069 个基因**里抽，**包含表达量近 0 的基因**。这是 scanpy 默认，本项目**照用不改为自定义 pool**（改了就是偏离默认且无出处）；如实登记该点 | ⚠️ C | 如实登记 |
| ↳ 亚型标签口径 | 逐细胞得分 → **簇内均值 argmax**（与 GP6 标准 A 同口径，便于两处对齐） | ⚠️ C | 与 §M3-A.4 一致 |
| 🔴 **P2 的性质更正（重要，2026-09-17 实测）** | 实测经典核心集 **76 个基因中 55 个（72%）本来就在论文 MP 里**。⇒ **P2 = P1 ∩ 经典 不是「两套独立口径」，而是「P1 剪掉噪声基因后的高特异版本」**。P1 与 P2 的一致性**只检验特异性，不检验独立性**——两者同源，一致是应该的，不能当作交叉验证通过。**真正的独立第二来源只有 Table S3 的 KAC 签名**（与 MP6 仅重叠 8/50） | ⚠️ C | 实测 2026-09-17；**此更正须随 GP8a 结果一并报告，不得把 P1 vs P2 一致说成「两套方法互证」** |
| ↳ 实测重叠明细（逐亚型） | AT1 10/13、AT2 8/12、Club 7/9、Ciliated 12/13、Basal 10/12、Tumor-KAC 8/17；合计 **55/76 = 72%** | ✅ S | 实测 |
| ↳ P2 规模 | **55 基因**（P1 为 445 基因槽 / 去重 388） | ✅ S | 实测 |
| **面板检出率诊断（GP8a 前置，不产标签）** | 在已签字口径（标准A 上皮 **141,105 核**）上算逐基因检出率。脚本 `05_annotation/01_panel_detection_check.py`，产物 `results/05_annotation/epi_panel_detection{,_per_gene}.csv`。**实测：8 个亚型全部有信号**，逐亚型检出率中位：AT2 34.8% / Tumor-KAC 30.6% / Tumor-stress 29.9% / AT1 27.6% / KAC-inflam 20.3% / Club 18.4% / Basal 14.6% / **Ciliated 7.2%**；**无基因零检出**。复核项：脚本读到的 barcode 清单 sha256 = `a952857e…d34`，**与登记值一致**（口径未被改动） | ✅ S | 实测 2026-09-17 |
| ⚠️ 诊断出的两个风险点 | ① **Ciliated 信号最弱**（中位 7.2%、p10 4.8%），亚聚类时**可能被合并进邻近亚型**；② **Basal 的经典 marker 恰是检出最差的一批** —— `KRT5` **2.2%**、`TP63` 3.3%、`KRT17` 3.5%、`COL7A1` 2.4% ⇒ basal 身份**只能靠 MP 基因撑**，靠不了教科书 marker；basal 稀少或混入时风险高。③ KAC 签名（Table S3）中 `KLK6` 检出 0.14%，近乎无效 | ⚠️ C | 同上；**须随 GP8a 结果一并报告** |

### M3-A.5c 上皮亚型 marker 面板 —— **经典 marker 口径（现行）**

> ✅ **现行口径。**用户 2026-09-17 指令原话：「**就不按什么kac去做，就按照经典marker，kac是他提出的吗，
> 按照权威经典marker做重聚类和celltype那个算法比比结果**」。
> ⇒ ① 论文 MP / KAC 体系**不再作 GP8a 面板**（降为历史，见 §M3-A.5b 存档）；
> ② 面板改用**权威一次文献的经典 marker**；③ 亚聚类结果**须与 CellTypist 对比后并报**。

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **面板来源（权威）** | **Travaglini 2020 *Nature* 587:619 的 Table S1「Canonical markers」** —— 该表**逐细胞类型直给经典 marker**，是教科书级的权威一张表。已冻结入库：`data/external/travaglini2020/TableS1_canonical_cell_types_markers.xlsx`，sha256 `465d7e83220bf4cc…755534e7` | ✅ S | Travaglini 2020 *Nature* 587(7835):619-625；PMC7704697 |
| ↳ 为什么它是最合适的"权威经典" | 它是**人工整编的教科书式规范表**（原始 60 行 × 8 列，按 Epithelium/Endothelium/Stroma/PNS/Immune 分块），**不是聚类富集输出**，故不携带本数据的任何信息 —— 作面板无循环性 | ✅ S | 原表结构实测 |
| ↳ 交叉校验源（第二意见） | **Habermann 2020 *Sci Adv* 6:eaba1972** 补充材料（已冻结 `data/external/habermann2020/aba1972_supplementary_material.pdf`，sha256 `155e9b2128f82bcf…1b6be81`）。**仅用于核对 Travaglini 表里存疑的条目**，不替代主来源 | ⚠️ C | Habermann 2020 |
| ✅ **旧面板引用缺陷 —— 已裁定并已修** | [`05_annotation/marker_panel.py`](../05_annotation/marker_panel.py) 原 4 个谱系（上皮/T·NK/B·浆/髓系）列 `VieiraBraga2019` 作出处，但**该文没有 marker 表**（补充 XLS 为 Tables S1–S10：细胞计数 / OMIM 基因表 / 临床元数据 / Fisher p 值 / T 细胞计数 / CellPhoneDB / 抗体表），**上皮亚型 marker 只出现在正文图里，不可逐基因追溯**。**用户 2026-09-17 裁定：「只改引用、不换基因、不重跑」**。⇒ 已执行：`src` 改为 `Travaglini2020`（可逐基因追溯），另设 `src_crosscheck=Habermann2020`（第二意见）；**80 个基因一个未动**。原引文一律移入新增的 `UNVERIFIED_SRC` **登记备查**（`VieiraBraga2019`/`Reyfman2019`/`Zilionis2019`/`Gillich2020`/`Lambrechts2018`/`Guo2018`），**不静默删除**。新增 `assert_provenance()` 硬校验，`src` 混入未核实出处即报错 | ✅ S | 2026-09-17 核对 Vieira Braga 2019 补充材料后确认；用户裁定后执行 |
| 🔴 **GP6 故意不重跑（用户决定）** | 用户原话「**只改引用、不换基因、不重跑**」。理由成立：本项目标准 A 的判定完全由 `genes` 列表驱动（`src` 是出处注释，不参与计算）⇒ **基因集不变 ⇒ 现行 GP6 的 κ=0.8539 与全部簇标签逐字不变**。重跑只会得到同样的数，属空转。**本行即为此决定的登记** | ✅ S | 用户 2026-09-17 决定 |
| ↳ 但 3 个弱出处基因须随结果报告 | `SPP1`（Table S4 最高 `pct_in` 仅 0.21）/ `FAP`（logFC 1.67 但 `pct_in` 0.15）/ `KDR`（最高 logFC 落在 Goblet，归属别扭）—— 保留但**在 `marker_panel.py` 的 `WEAK_PROVENANCE` 中显式标记**，报 GP6 结果时不得与其余 77 个同等对待 | 🔴 | 2026-09-17 出处审计 |
| **亚型集合（拟）** | Travaglini Table S1 上皮全 11 型：**AT1 / AT2 / Club / Ciliated / Basal / Goblet / Mucous / Serous / Ionocyte / Neuroendocrine / Tuft** | ⚠️ C | 待签字 |
| ↳ 基因（逐型，Table S1 原文） | AT1 `AGER, PDPN, CLIC5`；AT2 `SFTPB, SFTPC, SFTPD, MUC1, ETV5`；Club `CYP2F2, SCGB3A2, CCKAR`；Ciliated `FOXJ1, TUBB1, TP73, CCDC78`；Basal `KRT5, KRT14, TP63, DAPL1`；Goblet `MUC5B, MUC5AC, SPDEF`；Mucous `MUC5B`；Serous `PRR4, LPO, LTF`；Ionocyte `CFTR, FOXI1, ASCL3`；Neuroendocrine `CALCA, CHGA, ASCL1`；Tuft `DCLK1, ASCL2` | ✅ S | 原表逐字（**未加任何一个表外基因**） |
| ↳ ⚠️ 原表面向新鲜组织，非肿瘤 | 该表是**健康人肺图谱**的规范 marker。LUAD 中 many 上皮细胞是恶性的，**用正常 marker 给恶性细胞打标签，本身就是有偏的** —— 此点须随结果报告，不得掩盖 | ⚠️ C | 如实登记 |
| 矩阵存在性（实测） | 41 个基因槽中 **38 个在本矩阵**。**缺 3 个**：`CYP2F2`（**小鼠基因**；人同源 `CYP2F1` **在**，须显式登记换用）、`DAPL1`、`PRR4` | ✅ S | 实测 2026-09-17 |
| ↳ 补齐策略（拟） | 仅对上述 3 槽，用**同一文献**的邻接信息补：Club 补 `SCGB1A1`（Travaglini 亦列为分泌细胞经典）；Basal 补 `KRT15`/`KRT17`（同表簇富集列，非另找文献）；Serous 的 `PRR4` 缺失则**该型 3 基因只剩 2**，如实记录不补。**不从犄角旮旯找替身** | ⚠️ C | 待签字；用户要求「marker 要来自权威文章」 |
| 🔴 **恶性/正常不由本面板判** | 经典 marker **无法区分恶性与正常肺泡上皮**（LUAD 常保留 NKX2-1/SFTPC 等肺泡身份）。⇒ 本面板**只出"正常上皮亚型"标签**，**不设"tumor"亚型**。恶性身份一律由 **GP2 CNV** 判（铁律 R2），两者**交叉报告**（`亚型 × CNV±` 列联表），**绝不合并** | ✅ S | R2 + 资深判断，本项目约定 |
| 稀有型处置（预注册） | Serous / Ionocyte / Neuroendocrine / Tuft 在正常肺合计 <0.3%。**照算、照报**，但若某型在任一簇的 argmax 胜出且该簇 **<20 细胞**，须标 `rare/likely-spurious` 单列上报，**不得当作发现** | ⚠️ C | 本项目约定 |
| `score_genes` 参数 | 沿用 §M3-A.5b 已核实的注册值：`ctrl_size=50, n_bins=25, random_state=0, gene_pool=None, use_raw=False`（scanpy 1.9.8 签名现取） | ✅ S | 同 §M3-A.5b 第 314–315 行 |
| 标签口径 | 逐细胞得分 → **簇内均值 argmax**（与 GP6 标准 A 同口径） | ⚠️ C | 与 §M3-A.4 一致 |
| **与 CellTypist 的对比（用户明确要求）** | GP8a 出标签后，用 **CellTypist** 在**同一上皮子集**上独立注释，两者**并报**：逐细胞 Cohen's κ + 逐簇混淆矩阵。口径同 GP6 标准 B1（`Human_Lung_Atlas`，Sikkema 2023 *Nat Med* 29:1563）。**稀有型与"未定"单列**，不混进主 κ | ⚠️ C | 用户 2026-09-17 指令「和 celltype 那个算法比比结果」 |
| 🔴 对比的定位（防误读） | CellTypist 模型 = **HLCA 整合图谱**，而 HLCA 整合了 Travaglini 2020 **自身**。⇒ 两者**不同源但同根**，一致是**应该的**、不构成独立验证；**不一致才是信息量所在**。此性质须随结果写明，**不得把"高度一致"报成"互相验证通过"** | ⚠️ C | 与 GP6 §4 残留非独立性同源 |

#### ⚠️ 经典面板的**检出率现实**（2026-09-17 实测，**先于 GP8a，重要**）

执行口径已落盘：[`05_annotation/epi_classic_panel.py`](../05_annotation/epi_classic_panel.py)
（生成器 `build_epi_classic_panel.py`，**重跑逐字节一致**，sha256 `5cedb67f6a55fb0a…40a6e730`）。
逐基因检出率诊断（同一批 141,105 核）：`results/05_annotation/epi_panel_detection{,_per_gene}.csv`。

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| 🔴 **只有 5 个经典 marker 检出率 >50%** | `SFTPB` 88.9% / `SFTPC` 80.2% / `MUC1` 80.0% / `SFTPD` 59.7% / `AGER` 55.6% —— **全部是 AT2/AT1** | ✅ S | 实测 |
| 🔴 **其余亚型的 marker 检出率极低** | Club 最好 `SCGB3A2` **36.3%**；Ciliated 最好 `FOXJ1` **10.7%**；Goblet `MUC5B` 4.7%；Basal `TP63` 3.3% / `KRT5` **2.2%** / `KRT14` 0.17%；Serous `LPO` 0.32%；Tuft `ASCL2` 0.51%；Neuroendocrine `ASCL1` 0.38% / `CHGA` 0.13% / `CALCA` 0.11%；Ionocyte `FOXI1` 0.06% / `ASCL3` **0.02%**；Club `CCKAR` **0.0099%** | ✅ S | 实测 |
| 🔴 **我此处的原推断已被推翻（自我更正，保留原话以示追溯）** | ~~「标签会被 AT1/AT2 系统性吞掉 —— 其余 8 型的 marker 表达量低到 argmax 几乎不可能胜出」~~ **此句作废**。机制上不成立：`score_genes` 的对照基因**按表达量分箱后从同一箱抽**（`ctrl_size=50, n_bins=25`）⇒ 得分对任何 marker 集**中心化在 0**，低检出集**不处于基线劣势**。已改为**实测**（见下一行） | 🔴→✅ | scanpy 1.9.8 `_score_genes.py` 源码现取；详见 §M3-A.5d |
| **实测结果（替代原推断）** | 在**已冻结的 GP6 簇**上用经典面板跑 argmax 探针（`05_annotation/02_classic_panel_argmax_probe.py`，seed0、不重聚类、不调参）：**10 型中 5 型确有簇胜出** —— AT2 9 簇/99,941 核（71%）、AT1 1 簇/30,434、Ciliated 1/6,564、Basal 1/2,551、Goblet-Mucous 1/1,615；**从未胜出**：Club / Serous / Ionocyte / Neuroendocrine / Tuft。⇒ 原判断**过强**，但 AT2 独大与 4 个稀有型全灭**确为事实**，须在设计里直面 | ✅ S | `results/05_annotation/epi_classic_argmax_probe.csv` |
| 根因（三条，均已核实） | ① **核转录组**：纤毛/分泌类基因的 mRNA 主要在胞质，核内本就低；② **10x Flex 探针法 + 核 UMI 低**（本数据中位 nFeature 945）；③ **面板只有 2–5 个基因槽，容错为零**（一个基因差就塌） | ⚠️ C | ①② 已知；③ 本项目约定 |
| 对照：论文 MP 为何"看起来好" | MP 逐型检出中位 0.15–0.35（AT2 34.8% / Ciliated 7.2%），因 **MP 是 NMF 模块，含大量高表达/管家基因** —— 这正是它的**优点也是它的噪声**。⇒ 两者是**敏感性 vs 特异性**的交换，不是谁对谁错 | ✅ S | 实测对照 |
| 矩阵缺失（实测） | 经典面板 41 槽中 **38 个在矩阵**；缺 `CYP2F2`（**小鼠基因**，已按 `ORTHOLOG_FIX` 换 `CYP2F1`，在）/ `DAPL1` / `PRR4` | ✅ S | 实测 |
| Goblet 与 Mucous **已合并** | 两型在 Table S1 里共用 `MUC5B`（Goblet 仅多 `MUC5AC`/`SPDEF`）⇒ 这组 marker **无法把二者分开**，合并为 `Goblet/Mucous` 并登记，**不硬拆** | ✅ S | Table S1 原文 |
| ✅ **标签口径已裁定（用户 2026-09-17）** | 用户选定「**先聚类，再用经典面板按簇注释**」⇒ 与 Travaglini 2020 Methods 的**原procedure一致**（原文："Clusters were assigned a canonical identity based on enriched expression of these marker genes."）。**不是**逐细胞 argmax 后投票。⇒ **GP8a 可以开跑**，但注释单位是**簇**，须逐簇报富集 marker 与判读理由，**不得**只报一个标签 | ✅ S | 用户 2026-09-17 决定；Travaglini 2020 Methods 原文 |
| ↳ 簇注释必须直面的事 | 上述探针显示 AT2 会独占多数簇、稀有型可能**一个簇都没有**。⇒ 报结果时须直说「**本数据无法解析出 Club / Serous / Ionocyte / Neuroendocrine / Tuft**」，**不得**为了凑齐 10 型把某个 AT2 簇硬掰成稀有型 | 🔴 | 实测探针 + 本项目约定 |

#### M3-A.5e 全部 6 谱系的亚型经典面板（用户 2026-09-17「几个谱系都需要重聚类」）

执行口径已落盘：[`05_annotation/classic_panels.py`](../05_annotation/classic_panels.py)
（生成器 [`build_classic_panels.py`](../05_annotation/build_classic_panels.py)；**重跑逐字节一致已实跑验证**，
产物 sha256 `46244cdc0a459dfe…17bd493`，生成器 sha256 `f5dd338d48c8721c…669cdab3`）。
**来源唯一**：Travaglini 2020 Table S1
（sha256 `465d7e83220bf4cc…755534e7`），**不另找文献**。共 **6 谱系 / 40 亚型**。

| 参数 | 值 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| ↳ 与上皮那份的一致性 | 生成器**硬校验**：本文件的上皮块与已冻结的 [`epi_classic_panel.py`](../05_annotation/epi_classic_panel.py) **逐基因一致**，不一致即 abort ⇒ 两个文件不会漂移 | ✅ S | 生成器断言 |
| ↳ 🔴 **本行 sha256 的两次更正（留痕，不静默）** | 2026-09-17 PR 前复核发现：本行先后登记过 `da9c1cbccdf6c055…bad242701` 与 `cf433ca1e7d5d479…a3c58d`，**两个都对不上实际文件**（当时是对着一个临时/中间产物算的）。现按**当前工作树实测**重新登记，并**同时给出生成器与产物两个哈希**：产物 `46244cdc0a459dfe…17bd493`、生成器 `f5dd338d48c8721c…669cdab3`。**「重跑逐字节一致」这条是当天实跑验证的**（`build_classic_panels.py --out <临时路径>` 之后 `diff` 为空），不是推断；同次运行还复核了「全部替换词均在论文 Table S4 中出现」与「上皮块与冻结版逐基因一致」两条断言 | ✅ S | `sha256sum` + `diff` 实测，2026-09-17 |
| 各谱系亚型数 | 上皮 10（Goblet/Mucous 已并）· 内皮 4 · 成纤维 7 · T·NK 6 · B·浆 2 · 髓系 11 | ✅ S | Table S1 |
| 🔴 **对原表的偏离共 5 类，逐条登记** | ⓪ `ORTHOLOG_FIX`：`CYP2F2`（小鼠）→`CYP2F1`（人同源，在矩阵）<br>① `ALIAS_FIX`：`CD8`→`CD8A`、`CD16`→`FCGR3A`、`MHCII CLEC9A`→`CLEC9A`（原表写的是描述性名称/缩写，非基因符号；**逐个核实矩阵里有该符号才换**）<br>② `DROP_NON_SYMBOL`：`MHCII` —— 真符号 `HLA-DRA`/`HLA-DRB1` **都不在矩阵**（本 10x Flex 探针板未覆盖 MHC-II）⇒ **丢弃并登记，不用别的基因顶替**<br>③ `SOURCE_TABLE_DEFECTS`：见下两行<br>④ `UNANNOTATABLE`：见下 | ✅ S | 2026-09-17 逐基因核对矩阵 |
| 🔴 **原表硬伤 (a)：CD4 行里写着 CD8** | Table S1 `CD4+ Mem/Eff Cell` 行的 marker 是 `CD3E, **CD8**, COTL1, LDHB` —— 一个 CD4 细胞类型列了 CD8。**这是源表的错**。⇒ **不修**（改了就成了我们的判断），**整行标记为存疑并原样保留**，注释时**不得据它下结论** | 🔴 | 源表第 42 行 |
| 🔴 **原表硬伤 (b)：Basophil 与 Mast 行完全相同** | 两行基因都是 `MS4A2, CPA3, TPSAB1` ⇒ **结构上不可分辨**；且 `TPSAB1` **不在矩阵** ⇒ 实际只剩 `MS4A2`+`CPA3`。合并为 `Basophil/Mast` 上报，**不硬拆** | 🔴 | 源表第 47/48 行 |
| 🔴 **`Bronchial Vessel` 无法注释** | Table S1 第 20 行 markers 列为**空** —— 连基因槽都没有 ⇒ 本面板**注释不了**该亚型。**如实上报，不编造 marker** | 🔴 | 源表第 20 行 |
| 🔴 **非上皮面板比上皮粗（源表性质）** | **11 个亚型只有 1–2 个基因槽**（容错为零）：`Vein` 1 · `Capillary` 1 · `Eosinophil` 1 · `Nonclassical Monocyte` 1 · `Tuft` 2 · `Artery` 2 · `Lymphatic` 2 · `Fibroblast` 2 · `mDC1` 2 · `mDC2` 2 · `Classical Monocyte` 2。**这不是本项目的选择，是 Table S1 本身如此**。⇒ 解读非上皮亚型时须按 §M3-A.5c「检出率现实」同样谨慎，**不得把"某型没检出"直接说成"该型不存在"**。<br>🔴 **本条按源表口径计数，不含矩阵内损耗**：`Basophil/Mast` 源表 3 槽、但 `TPSAB1` 不在矩阵（见上上行的原表硬伤 b）⇒ **运行时实为 2 槽**；`mDC2` 源表 3 槽、`MHCII` 被丢弃 ⇒ 同样降为 2 槽（`N_SLOTS` 已按丢后计为 2）。故**运行时**「1–2 槽」的亚型比 11 个**更多**；脚本每次运行都按**矩阵内实际可用槽数**重算 `thin` 清单并打印，**以运行时那个清单为准**，本行的 11 个只是源表底数 | 🔴 | Table S1 基因槽计数 + 矩阵内缺失实测 |
| 🔴 **`CYP2F2` 是鼠源 —— 源表物种笔误，已查实** | 用户 2026-09-17 追问「这个怎么还有鼠源marker」。**查证结论**：<br>① **`CYP2F2` 是小鼠基因，人的同源基因是 `CYP2F1`**。旁证：文献中成对出现「CYP2F2 knockout 与 CYP2F1 humanized mice」（把小鼠 Cyp2f2 敲掉、把人的 CYP2F1 转进去），以及「mouse-specific lung tumors from CYP2F2-mediated metabolism」——均表明 CYP2F2 无人物种对应物。<br>② 它成为小鼠 club cell 招牌 marker 的原因：小鼠**萘（naphthalene）肺损伤模型**靠 Cyp2f2 代谢萘生成毒物、专杀 club 细胞，故小鼠文献中 club cell marker 长期写作 Cyp2f2。<br>③ **Travaglini 2020 是人肺图谱**，全表其余基因均为人的 ⇒ Table S1 此格系**物种笔误**。<br>④ **决定性证据来自该论文自己**：其 **Table S4**（用他们自己的人类数据算出的簇富集 marker）在 **Cluster 8 Goblet** 与 **Cluster 9 Mucous** 两张 sheet 里写的都是 **`CYP2F1`**（人）⇒ **论文自己的分析用的就是 CYP2F1**，只有人工汇总的 Table S1 写成了鼠基因。<br>⑤ 本矩阵：`CYP2F2` **不在**，`CYP2F1` **在** | ✅ S | 2026-09-17 逐格核对 Table S1 / Table S4 + 物种关系文献核查 |
| ✅ **用户裁定：按「和我们数据对得上的那张原文 marker」做** | 用户 2026-09-17 原话：「**按照我们数据原文那张marker去做**」。⇒ 取基因以**论文自己的数据（Table S4 / 本矩阵实际存在的那个人的基因符号）**为准：`CYP2F2`→`CYP2F1` 保留，理由即上一行的 ④⑤（同篇论文自己的数据这么用，且本矩阵只有人的那个）。<br>**同一规则适用于同类情形**：字面写的是描述性名称/缩写、而矩阵里存在其正式符号的（`CD8`→`CD8A`、`CD16`→`FCGR3A`、`MHCII CLEC9A`→`CLEC9A`）照此resolve；正式符号**矩阵里没有**的（`MHCII`→`HLA-DRA`/`HLA-DRB1` 均不在）**只报缺失、不找替身**。<br>**仍然不修**：`CD4+ Mem/Eff` 行把 CD8 写进 CD4 型——那是**作者的判断错误**，不是符号问题，改了就是替作者做判断 | ✅ S | 用户 2026-09-17 原话 + 与 Table S4 逐格对账 |
| 其余（`RARE` / `CONTAM_CHECK` / `MERGE`） | `RARE` = Serous/Ionocyte/Neuroendocrine/Tuft/Eosinophil/Megakaryocyte；`CONTAM_CHECK` = PTPRC/CD3D/COL1A1/PECAM1/CD68/EPCAM（**仅作 QC 列，不进 argmax**）；`MERGE` = Goblet+Mucous、Basophil+Mast | ⚠️ C | 本项目约定 |
| `assert_provenance()` | 硬校验：面板里**不得残留**任何未处理的非符号；原表硬伤必须带 🔴 标记。任一不满足即报错 | ✅ S | — |

#### M3-A.5d marker 出处审计（用户 2026-09-17 指令「你给我查清楚marker到底哪里来的」）

审计对象：[`05_annotation/marker_panel.py`](../05_annotation/marker_panel.py) 的 6 谱系面板（80 基因）。

| 结论 | 内容 | 状态 | 证据 |
| :--- | :--- | :---: | :--- |
| 🔴 **文件名里的引用是错的** | 该文件 4 个谱系（上皮/T·NK/B·浆/髓系）列 `VieiraBraga2019` 作出处，但**该文没有 marker 表**（补充材料为 Tables S1–S10：细胞计数 / OMIM 基因表 / 临床元数据 / Fisher p 值 / T 细胞计数 / CellPhoneDB / 抗体表）。**该引用不可逐基因追溯** | ✅ S | 独立核对 Vieira Braga 2019 补充材料，2026-09-17 |
| ✅ **但基因本身不是"犄角旮旯"来的** | 80/80 基因**全部见于 Travaglini 2020 Table S4**（带 avg_logFC / pct_in_cluster / pct_out_cluster / p_val_adj），即**每一个都能追溯到一篇一次文献的表** | ✅ S | `data/external/travaglini2020/TableS4_cluster_enriched_markers.xlsx`（sha256 `abca9f72c5bb8d47…29031073`）逐基因反查 |
| ✅ **谱系归属 71/80 一致** | 判据：该基因被报富集的细胞类型中，至少有一个属于它应属的谱系（logFC≥0.25 且 pct_in≥0.5） | ✅ S | 实测 |
| ⚠️ **9 个需单独说明** | `TRAC`/`CD4`/`CD19`/`IGHG1`/`IGKC`/`PDGFRA` —— **实际合格**，只是落在 **Table S4 的 (SS2) 大表**里（该表基因数 1600–4300，被我的初筛排除）；逐条复查确认归属正确。**真正偏弱的是 3 个**：`SPP1`（全表最高 pct_in 仅 0.21）、`FAP`（logFC 1.67 但 pct_in 0.15）、`KDR`（最强富集落在 Goblet，生物学别扭） | ⚠️ C | 实测逐条复查 2026-09-17 |
| ⚠️ **Travaglini Table S4 是"簇富集 marker"表，不是"泛 marker"表** | 按构造只收 `pct_out ≤ 0.30` 的基因。**教科书级泛 marker 天生进不去**（`SFTPB` pct_out 0.78、`COL1A2` 0.75、`DCN` 0.90）⇒ 40 个基因在 Table S4 里只出现在低特异尾部。**这不是数据缺陷，是两类表用途不同** | ✅ S | 实测；两类表分工须随报告说明 |
| ✅ **Travaglini Table S1 是人工整编参考表** | Methods「Cell clustering, doublet calling, and annotation」逐字：「**Clusters were assigned a canonical identity based on enriched expression of these marker genes.**」⇒ Table S1 是作者手工整编、**供人工簇级判读查表用**，**从不喂给自动标注器**，且**不是算法输入** ⇒ 作面板**无循环性** | ✅ S | Travaglini 2020 *Nature* Methods 原文；PMC7704697 |
| ✅ **已找到"真·注释 marker 表"** | **Habermann 2020 Table S3「Marker genes used for cell-type annotation」** —— 覆盖**全部 6 谱系**（上皮/免疫/内皮/间质），带**正负 marker**。⚠️ 但它在补充 PDF 里是**图片**（`data/external/habermann2020/` p28），**无法程序化提取**，手工转录违反"非手抄"原则 ⇒ 只作**核对第二意见**用 | ✅ S | Habermann 2020 *Sci Adv* 6:eaba1972 补充材料 p28 实测 |
| 🔴 **Travaglini 是"先聚类、后人工注释"** | 工具链：Seurat v2.3 + Louvain `FindClusters` + MAST `FindMarkers` + **人工判读**；**全程无逐细胞打分、无标签迁移**。⇒ 簇级注释是**该文献的原生做法**，不是我的发明 | ✅ S | Travaglini 2020 Methods 原文 |
| ⚠️ **但"argmax 更差"没有文献基准** | 方法学综述（Heumos 2023 *Nat Rev Genet* 24:550「Mapping cell clusters to cell identities」节；Kiselev 2019 *Nat Rev Genet* 20:273；Abdelaal 2019 *Genome Biol*）**一致把注释framed 在簇级**，但**没有一篇基准研究定量证明逐细胞 argmax 更差**。⇒ 只能说"标准做法如此"，**不得说"argmax 已被证明错"** | ⚠️ C | 文献调研 2026-09-17 |
| 🔴 **我的判断被实测推翻（自我更正，须登记）** | 我曾判断「经典面板 argmax 会被 AT1/AT2 系统性吞掉」。**机制上不成立**：`score_genes` 的对照基因**按表达量分箱后从同一箱抽**（`ctrl_size=50, n_bins=25`，复刻 Seurat `AddModuleScore`，出处 Satija 2015 / Tirosh 2016）⇒ 得分对任何 marker 集**中心化在 0**，低检出集**不处于基线劣势** | ✅ S | scanpy 1.9.8 `_score_genes.py` 源码现取 |
| 实测探针结果（GP6 冻结簇，非新聚类） | 脚本 `05_annotation/02_classic_panel_argmax_probe.py`，产物 `results/05_annotation/epi_classic_argmax_probe.csv`。**10 个亚型中 5 个赢下簇**（AT2 9 簇/99,941 核；AT1 1 簇；Ciliated 1；Basal 1；Goblet/Mucous 1）；**从未胜出：Club / Serous / Ionocyte / Neuroendocrine / Tuft**。⇒ 我"会塌掉"的说法**过强，撤回**；但 **AT2 占 71% 上皮核**、**Club 一簇未赢**是真实局限 | ✅ S | 实测 2026-09-17 |
| ⚠️ 该探针的效力边界 | 只在 **GP6 的全身簇**上做，仅 13/45 簇属上皮 ⇒ **不能替代 GP8a 的检验**（GP8a 是上皮自己重聚类，AT2 那一大坨会裂开） | ⚠️ C | — |

#### KAC 出处裁定（回答用户问「kac是他提出的吗」）

| 事实 | 内容 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| **"KAC"这个名字是谁的** | **是作者自己团队的**，但**不是本篇 2026 论文首提**。原文写：「**In line with our recent report 19**, we also noted a **KRT8 high alveolar intermediate cell (KAC)** subset…」。**ref 19 = Han G, Sinjab A, Rahal Z, …, Kadara H (2024). *An atlas of epithelial cell states and plasticity in lung adenocarcinoma.* Nature 627:656–663** —— Sinjab 与 Kadara 同为 2026 本篇作者。⇒ **KAC 是同一团队在自己 2024 年 Nature 论文里先提的**，本篇是沿用 | ✅ S | 原文 ref 19 全文；本仓库 `data/external/peng2026_cancercell/` |
| **这门生物学是不是他们的** | **不是**。KRT8-high 肺泡中间态在**他组之前**已由**其他团队**独立建立：ref 32 = **Choi et al. 2020 *Cell Stem Cell* 27:366**（小鼠急性肺损伤，AT2 来源的"损伤相关瞬时祖细胞"）；ref 33 = **Strunz et al. 2020 *Nat Commun* 11:3559**（人肺纤维化中的 **Krt8+ 过渡态干细胞**）。本篇自己也写「previous studies…found that proinflammatory signaling was important for generation and/or maintenance of KRT8 high alveolar intermediates…32,33」 | ✅ S | 原文 refs 32/33 全文 + 正文引述 |
| **结论（一句话）** | **细胞状态本身是公认的、别人先发现的；"KAC"这个缩写是作者组的自家命名**，源自他们 2024 年 Nature。⇒ 把它当"已确立的独立细胞类型"用是有风险的，因为**它的定义权在作者组手里**，且本篇内两处 KAC 定义只重叠 8 个基因（见 §M3-A.5b 第 305 行） | ✅ S | 综合上两行 |

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
