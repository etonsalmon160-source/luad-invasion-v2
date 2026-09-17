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
| 上皮亚群命名 | ciliated / club / basal / AT1 / AT2 / **AIC**（alveolar intermediate）/ **KAC**（KRT8-high alveolar intermediate）/ tumor | ✅ S | 原文 |
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
| **`UP.DR` 强制与 `LOW.DR` 相等**（**预注册：2026-09-16，GP2 开跑前**；⚠️ **书面签字待补**） | `LOW.DR = UP.DR = 0.05` | ⚠️ C | **偏离 `copykat` 包默认**（默认 `LOW.DR=0.05, UP.DR=0.1`）。理由：①`UP.DR` 名不副实，源码 **`:189`**（及 `:192`）用作**下界** `DR2 >= UP.DR`，是驱动最终分割的基因检出率门槛；②源码 **`:57-60`** `if(nrow(rawmat) < 7000){ UP.DR <- LOW.DR ; print(...) }` 是**无条件覆写**，显式传参**拦不住**，只有两值相等才能使有效阈值全队列一致；③实测 75 样本（`results/03_cnv/prereg_gene_tiers.csv` 的 `genes_ge_*` 列）：统一 **0.05** → 最终基因数 min **3,157** / p25 5,903 / 中位 **6,742** / max 9,145，**0 个**样本 <2,000；统一 **0.10** → min **1,791** / p25 3,756 / 中位 **4,494** / max 7,523，**1 个**样本 <2,000，`win.size=25` 按基因序平滑会跨极大基因组距离、分割退化；④不强制则 **23** 个样本用 0.05、**52** 个用 0.1，且覆盖率与分期相关（Normal 最高），构成沿 M4 进展轴的混杂。**此为预先决定，非事后调参**（法则 3.2）。<br>⚠️ **已刷新（2026-09-16）**：本行早前引用的 min 3,528 / 中位 7,429、min 1,231 / 中位 3,723、3 个样本 <2,000、27/48 分布，是 ③ 的**建模值**版本；`02_prereg_gene_tiers.R` 重写为调用 `copykat_chain()` 后 `genes_ge_*` 已换成**精确值**，上表为精确值。**结论方向不变**（0.05 仍更保守），但决策时依据的数字与现值不同，如实记此差异。<br>⚠️ **仍未消除的反对意见**：作者 README 明写 "I need to make sure that LOW.DR is smaller than UP.DR"，本项目此决定**与该句直接冲突**，不假装有文献支持。<br>⚠️ **签字状态**：`results/03_cnv/GP1_report.md` §7 曾把本项列为【需裁决】。**本行不主张已获签字** —— 2026-09-16 复核时用户询问"证据是否不支持"，本行据此复核确认**证据方向未变**（0.05 仍更保守），但**书面 go/no-go 仍待用户给出**<br>🔴 **2026-09-16 17:20 证据翻转（同一日内，加算 840 地板后）**：本行 ③④ 两条理由**建立在地板之前的口径上**。在地板 840 之后重算（`logs/tiers_floor840.log`、`prereg_gene_tiers_floor840.csv`）：<br>· 7000 覆写从 **23/75 塌到 2/75**（21 个样本 0.05→0.10，**无一反向**）。剩余两个是 `P16_AIS`（1,471 细胞，判据基因集 6,949）与 `P8_Normal`（126 细胞，6,272）。<br>· 因此 ④ 的"分期混杂"**基本消失**：Normal 45.8%→4.2%、MIA 50%→0%、IAC 25%→0%。<br>· 代价却**反向变大**：统一为 0.05 需要动的样本数 **52 → 73**（因为覆写后 73/75 的有效值已是 0.10），而买到的"均匀性"从 31% **缩到 2.7%**。<br>· ③ 的"0.05 更保守"仍未变（min 3,157 vs 1,791 量级不变），但**它已不是主要论据**，因为 73/75 本来就要用 0.10。<br>🔴 **本项目据此撤回原推荐**：**不做统一**，令 `LOW.DR=0.05, UP.DR=0.10` 走 `copykat` 默认 + 源码 :57 覆写；并在 GP2 报告中**逐样本申报有效 `UP.DR`** 及其与分期的关联（关联已弱化到可忽略）。**理由**：这笔交易的方向反了——代价从 52 涨到 73 个样本、收益从 31% 缩到 2.7%，且原④的混杂论据已被地板消除。**这是撤回，不是新增签字项**；原预注册从未获书面签字，故不构成事后调参。<br>⚠️ 仍存的分歧：作者 README "LOW.DR is smaller than UP.DR" **与默认参数一致**，故撤回后反而**与作者一致** |
| **`n.cores` 是结果相关参数，不得当加速旋钮调** | `n.cores=1`（全队列固定） | ✅ S | **实测证伪**（2026-09-16，R 4.2.2 / 20 核，见 [`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §4.3）：同一 `set.seed(1234)` 下 `mclapply(1:8, runif(1), mc.cores=1/2/4)` 三次结果**两两不相同**。机制：`CNA.MCMC.R:73`（MCMC 分段）与 `copykat.R:123`（dlm 平滑）用 fork 式 `mclapply`，子进程随机流按 **child 序号**分配，故 `mc.cores` 变 → 每个细胞拿到的流变 → **CNV 分段结果变**。`copykat.R:34 set.seed(1234)` 只保证**固定核数下**可复现。**处置**：GP2 全程固定该值并登记；若要改，属**换计算方法**，整批重跑，不得只重跑一部分 |
| **覆盖度地板**（**2026-09-16 由用户签字采纳**） | 负对照细胞 **nFeature ≥ 840**（`P13_Normal` 负对照中位 nFeature） | ⚠️ C | **无文献出处**，来自本项目**负对照稀释实验**（[`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §8）。做法：对纯二倍体 `P13_Normal`（801 细胞）按 f ∈ {1.0…0.2} 二项稀释 UMI 后重跑 copykat，看它**自己**在何处开始判出 `aneuploid`。**取更严的那个地板**：<br>· 预注册主口径 `f_knee = 0.70` → `C*` = **746**（判据：负对照 `rate_not_defined > 0.50`）<br>· **实测假阳性地板 `f = 0.80` → 840**（该处已有 **134** 个负对照细胞被判 `aneuploid`）<br>**采纳 840 的理由**：低覆盖度会**制造假阳性非整倍体**，而该地板比 `not.defined` 地板**更早触发**；取更严值以压假阳性。序列**非单调**（`0 → 0 → 134 → 84 → 41 → 61 → 13`），照原样记录，不做平滑。<br>⚠️ 边界：**n = 1 样本、无生物学重复**；稀释只模拟测序深度下降，真实低质量样本**只会更糟**，故本值是**下界**。<br>⚠️ **实测代价（2026-09-16 补算，签字时未知）**：按 M1 逐细胞表（`results/01_qc/gse308103_per_cell_qc.csv.gz`）实算，全队列 **648,945 → 376,906**，**弃 272,039 = 41.9%**。丢弃率**在各分期间大体均匀**（AAH 39.6% / AIS 41.8% / LUAD 42.3% / Normal 43.7% / MIA 46.8%），故**不是分期混杂**；但**样本层面极不均**（`P8_Normal` 弃 92.9%，仅留 126 细胞）。**跌破 200 细胞而转 `not_testable` 的仅 1 / 75**（`P8_Normal`）。**代价性质 = 统计功效（丢掉 42% 细胞 + 样本量不均），不是偏倚。**<br>⚠️ **仍需注意的方法学点**：840 源自**单个**样本 `P13_Normal` 的深度；对 31 个自身中位 nFeature < 840 的样本，它比 M1 的**逐样本自适应 MAD**（`nmads=3`）更严，即在本门槛下**部分样本的大部分细胞被排除**。是否改用**逐样本相对**口径，属未决<br>🔴 **实测的两个连带效应（2026-09-16 17:20 补算，签字时均未知）**：地板作用在**进 copykat 之前**，故它同时改了**分母**与**7000 覆写判据**，两个后果都必须重算（已做，`prereg_gene_tiers_floor840.csv`）：
> ⚠️ **2026-09-17 审计标缺（未修，留给 GP2）**：本行及其下两行的**全部实测数字是在旧 648,945 细胞集上算的**。
> 现行分析集已改为 **413,697**（论文 QC，2026-09-16），而 GP2 的输入是 **GP6 定出的上皮细胞**（该子集在 GP6 之前不存在）。
> 故：**地板值 840 本身**（来自 `P13_Normal` 稀释实验）**不受影响**，但**"648,945→376,906 / 弃 41.9% / `not.defined` 0.408→0.097 / 7000 覆写 23/75→2/75" 这些连带效应数字必须在新细胞集与上皮子集上重算**，
> 否则 GP2 会以过期的分母与检出率做判据。**此处如实标缺，不代为重算。**<br>① **7000 覆写 23/75 → 2/75**。机制：剔低覆盖细胞**抬高**了逐基因检出率 → 过 `LOW.DR` 的基因变多 → 判据分母 `n_after_LOWDR_fullgenes`（**全基因**、早于注释，`copykat.R:57`）中位 **7,680 → 9,136** → 7000 这条线落进更空的区域。**注意判据列是 `n_after_LOWDR_fullgenes`，不是 `n_genes_final`**（后者是注释+删周期/HLA 之后的集合，两者可差数千）。<br>② **`not.defined` 率中位 0.408 → 0.097（总量 233,582 → 30,901，÷7.6）**。原因：地板的判据与产生 `not.defined` 的判据**都是覆盖度**（`copykat.R:86-105` ToRemov2、`:194-213` ToRemov3），故地板**预先清掉了**那些细胞，而非仅缩小分母。**我曾预判此处"率会升高"，实测相反 —— 记录在案。**<br>⚠️ 因此本行与下行的"地板只改分母"隐含理解**是错的**，此处更正 |
| **传入 CopyKAT 的细胞集合**（**2026-09-16 新增；当前主口径已被实测否决**） | 🔴 **未定**。现行主口径"整个样本的全部细胞"**实测不可用**；`epi_call()`（仅上皮）**已证既非必要也非充分** | 🔴 | **无文献出处，来自本项目实测**（`results/03_cnv/GP1_report.md` §9.7–§9.9；判据预注册于 [`EPI_PREREG.md`](../results/03_cnv/EPI_PREREG.md)）。**根因**：`norm.cell.names=""`（我们的调用）令 copykat **自己猜**正常细胞基线（`copykat.R:166-170` → `baseline.norm.cl()`：`hclust(ward.D2)` 分 6 群、各群拟合三成分正态混合取 `SDM`、取 `SDM` 最小的簇当正常）。**实测 13 个臂里有 10 次**该步打印 `low confidence in classification`，之后**不拒绝作答**，静默回退 `baseline.GMM(max.normal=5,…)` 重取 `preN` 并把 `WNS` 改写成 `""`，最终标签以确定措辞输出。随后 `copykat.R:456` `if(cor(conses.diploid,conses.aneuploid)>=0.6)` 是**二值**分支：过线 → **全部判 diploid**；不过线 → 按 Wasserstein 把**每个簇**强制二分，**无"不确定"档**。<br>**实测后果（判决轴随输入漂移）**：<br>· 全细胞 → 判决 = **谱系**。`P11_Normal` 判"非整倍体"组 EPCAM 0.67/PTPRC 0.00；判"二倍体"组 LYZ 0.63/DCN 1.02。<br>· 仅上皮（`strict`）→ 判决 = **亚型**。`P11_Normal` 判"非整倍体"组 SFTPC 4.79；判"二倍体"组 AGER 2.88（免疫/基质 marker 全为 0，上皮判据确实生效）。<br>· **同一 `P11_Normal`，两条都说得过去的上皮判据给出 0.0%（argmax）与 32.3%（strict）** —— 收窄细胞集本该只降噪，却把非整倍体从 0% 抬到 32.3%。<br>**已撤回的补救方案**："改报连续 CNV burden、丢掉二值标签"**不足以**解决——基线猜错时整张 `final_results_bin_by_cell.txt` 都是相对错误基线算的。<br>**待验证的出路（"锚定 + 同质"须同时满足）**：① **锚定** = 显式给 `norm.cell.names`，不让它猜；② **同质** = 输入须为**同一亚型**的 恶性 + 正常 混合。仅锚定不够（AT1/AT2 并存时，参考之外的那一种照样被判非整倍体）。<br>🔴 **GP2 不得启动**，直至本行定稿并签字 |
| 退化样本处置 | `n_retained<200` → `not_testable`（**绝不填补为 Normal**）；`n_cnv_pos==0` → 作为真实结果上报，**不得放宽 `KS.cut` 重跑**；Normal/AAH `frac_cnv_pos>0.95` 或 IAC `<0.01` → `implausible` | ⚠️ C | 无文献阈值；本项目约定（防事后调参，法则 3.2） |
| **弃跑判据** | `not_testable` + `implausible` 合计 **>30%** → 停并上报 | ⚠️ C | 本项目约定；届时 R2 对该数据集不可满足，**禁用泛上皮 argmax 替代**。<br>🔴 **2026-09-16 实测评：本判据抓不住 §9.7 的失效模式。** 它的设计对象是"全部判成非整倍体"（Normal `frac > 0.95`）；而实测的失效是"在 0% 与 59% 之间随配置任意翻转"，`P13_Normal` 的 0.407 **远低于 0.95，不会触发**。本判据**必须改写**为能捕捉"判决对输入/参数的敏感性"（如：同一负对照在相邻配置下跨越某分界），否则形同虚设。改写方案待定，**不得**用当前版本充当门禁 |
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
| 协议 | 每个子集**从原始计数重启全流程**（口径同 §M3-A.1，不换算法）：子集内 **SCTransform**（同 `vst.flavor="v2"`，HVG 3000 在**子集内**重选）→ **PCA(50)** → **Harmony(`sample_id`)** → FindNeighbors(50) → 同网格 `{0.5…0.8}` × 5 种子 × **同 M3-A.3 全部指标（含指标4）** → 子集内双标准注释 | ⚠️ C | **本项目约定**——未找到正式证明该协议避免 double-dipping 的论文。最接近的选择性推断文献（Gao, Bien & Witten 2020, *JASA* 115:1622）讨论的是**选择+检验复用数据**，非聚类。**不冒充有引用** |
| ↳ 子集内 HVG 的偏差 | 论文的 3000 HVG 是在**全部 413,697 核**上选的；子集内重选 ⇒ 亚聚类空间**不与论文等价**（论文未做亚聚类，无可对齐对象） | ⚠️ C | 如实登记，**不声称与论文一致** |
| ↳ 子集规模下限 | 子集细胞数 **< 1000** ⇒ 该谱系**不做亚聚类**，只报 L1，如实记录（`SCTransform` 的 `ncells=5000` 与 `min_cells` 类过滤在小子集上不稳） | ⚠️ C | 本项目约定 |
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
