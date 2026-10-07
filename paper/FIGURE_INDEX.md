# 正文图索引（按服务器上已排好版的三张大图）

> 来源：`/home/eto/luad_v2/results/paper_figures/图{1,2,3}.pdf`（连 `composite/` 下同名副本）
> 页面宽度均为 **180 mm** ⇒ **跨栏满宽**（full-width），LaTeX 里用 `figure*` 环境。
> 每个面板的来源脚本 / 数据在末列，便于溯源。

---

## Figure 1 — 单细胞图谱（`图1.pdf`，180×178 mm，6 面板）

| 面板 | 内容 | 来源 |
|---|---|---|
| **a** | Single-cell atlas: lineage —— 399,579 核的 UMAP，按 6 大谱系上色 | `figs/10_panels_singlecell.R` → `P1a` |
| **b** | Lineage composition —— 六谱系占比（竖条） | `figs/10_panels_singlecell.R` → `P2a` |
| **c** | Single-cell atlas: stage —— 同一 UMAP，按五期上色 | `figs/10_panels_singlecell.R` → `P1b` |
| **d** | Single-cell atlas by stage —— 按期别分面（**保留谱系配色**） | `figs/10_panels_singlecell.R` → `P1c` |
| **e** | Lineage composition by stage —— 六谱系随期别的堆叠构成 | `figs/10_panels_singlecell.R` → `P2a` |
| **f** | L2 subtype composition by stage —— **39 个亚型**随期别的堆叠构成 | `figs/10_panels_singlecell.R` → `P2b` |

---

## Figure 2 — 空间解卷积与生态位（`图2.pdf`，180×193 mm，7 面板，**内容旋转 90°**）

| 面板 | 内容 | 来源 |
|---|---|---|
| **a** | RCTD lineage composition of each niche domain —— 七个域类型的六谱系权重构成 | `figs/35_domain_rctd_composition.py` → `P6a` |
| **b** | **空间矩阵**（5 期代表切片 × 5 列）：H&E ｜ RCTD 6 谱系 ｜ RCTD 39 亚型 ｜ Programs(T/NK/B/LAM) ｜ 生态位域(K\*=7) | `figs/32_spatial_matrix.py` → `P5` |
| **c** | RCTD subtype composition of each niche domain —— 同上，**39 个亚型**粒度 | `figs/35_domain_rctd_composition.py` → `P6b` |
| **d** | Domain composition per stage —— 七域随期别的堆叠构成 | `figs/52_domain_composition_by_stage.py` → `P22` |
| **e** | 四程序随期别：相对构成（堆叠）+ 绝对水平（折线） | `figs/51_exhaustion_by_stage.py` → `P21` |
| **f** | Univariable Cox regression —— 生态位域签名的预后森林图（连续） | `figs/39_domain_forest_plot.py` → `P9a` |
| **g** | 域类型 × 标志基因 dot plot | `figs/37_domain_dotplot.py` → `P8` |

> ⚠️ **面板 e 是负面对照**（正文 §13.1）：四程序的期别差异**不随病程单调**，配平前几乎就是深度的镜像；**不得**作为病程结论引用。

---

## Figure 3 — 靶点扰动逆向臂（`图3.pdf`，180×210 mm，8 面板）

| 面板 | 内容 | 来源 |
|---|---|---|
| **a** | 逆转签名的预后质控：疾病方向 × 预后方向四象限 + 保留/剔除量 | `figs/50_fig_signature_qc.py` → `P20` |
| **b** | 阳性对照：上胚层→中胚层，8,450 个扰动的因果分排序（TBXT 第 1） | `figs/41_fig_scmg_positive_control.py` → `P11` |
| **c** | 七个域类型各自的空间 CNV 分值（同切片内配对） | `figs/44_fig_domain_cnv.py` → `P14` |
| **d** | 化合物层：CMap × Cor-Spearman 秩-秩 + 34 个候选的机制构成 | `figs/42_fig_compound_reversal.py` → `P12` |
| **e** | 跨模态：单细胞轴 vs 空转轴的化合物排名（ρ = 0.73） | `figs/43_fig_compound_spatial.py` → `P13` |
| **f** | 深度–密度耦联：域签名存活率 vs 域深度（ρ = −0.89） | `figs/45_fig_depth_coupling.py` → `P15` |
| **g** | 哪个域最可逆转：34 个候选在七个域榜里的命中数 | `figs/49_fig_domain_reversal.py` → `P19` |
| **h** | 同一库三个打分口径的两两关系（相关矩阵 + top-N 重叠） | `figs/48_fig_method_agreement.py` → `P18` |

---

## 尚未进图但已在正文的内容（需决定：补图 / 只列表 / 放补充）

| 内容 | 现状 | 建议 |
|---|---|---|
| 队列与设计（25 患者 × 五期矩阵） | `P16_cohort_design` 已单独出图 | 可作 **Fig 1 前置面板**或**补充图** |
| 阳性对照的**逐位对拍**（WTCS 差 6.6e-14） | 仅正文文字 | 表格（见 TABLES_INDEX） |
| 劈半稳定性（ρ 中位 0.962，4 次） | `P17_split_half_stability` | 正文图缺；可并入 Fig 3 或补充 |
| 域 marker 热图（P7） | 已出图 | 与 Fig 2g（dot plot）二选一 |
| 域预后森林（**中位二分**版，P9b） | 已出图 | Fig 2f 用的是连续版；二分版可放补充 |
| 各阴性臂（发育谱三条、生态位闸 FAIL、H&E 天花板） | 仅正文文字 | 补充图或表 |
