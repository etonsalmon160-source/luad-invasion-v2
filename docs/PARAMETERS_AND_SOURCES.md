# 参数与出处对照表 (PARAMETERS & SOURCES)

> 为白皮书 §3 每步提供**可溯源参数**；每行标 **核对状态** 与 **出处类型**。
>
> **出处类型**：`S`=源码实测 · `D`=工具默认 · `P`=论文推荐 · `B`=方法学基准 · `C`=本项目约定
> **核对状态**：✅=本次已从**源码/官方文档现取确认** · 🟡=已附权威出处、值待逐条现取复核 · ⚠️=无单一文献（本项目约定）

---

## 阶段 1 · 单细胞

### 1a. SCMG 分支（参数取自官方源码 → ✅ S）
| 参数 | 值 | 出处 |
| :--- | :--- | :--- |
| 标准基因数 | **18,108**（`standard_genes.csv`, human_id） | SCMG `preprocessing/data_standardization.py` |
| 标准化 | total-count 归一到 **1e4** + `log1p` | SCMG `contrastive_embedding.py:198-199` |
| 编码器 MLP | `18108 → 512`，hidden `(2048, 2048)` | SCMG `contrastive_embedding.py:58` |
| 解码器 MLP | `512+64 → 18108`，hidden `(1024, 2048)` | SCMG `contrastive_embedding.py:60` |
| 嵌入 `batch_size` | **512** | SCMG `contrastive_embedding.py:210` |
| 扩散 `min_beta/max_beta` | **1e-4 / 0.02** | SCMG `manifold_generation.py:26-27` |
| 因果 `gene_stds` 下限 / z 截断 | **0.1 / ±5**；`causal_score = 匹配分 × 基因位移z` | SCMG `causal_prediction.py:22,36,59,62` |
| 流形 kNN | `n_neighbors=15`, `metric='cosine'` | ⚠️ C（本项目；SCMG 教程示例 k=30） |
| UMAP | `min_dist=0.3`, `spread=1.0` | D（scanpy/SCMG 教程） |

### 1b. 标准 / 主流分支
| 步骤/工具 | 参数 | 状态 | 出处 |
| :--- | :--- | :---: | :--- |
| QC | nFeature 500–10,000；nCount 1,000–60,000；mt<10% | ⚠️ C | 本项目约定（需再引同区间队列文佐证） |
| **scVI** | `n_latent=10`, `n_hidden=128`, `n_layers=1`, `gene_likelihood='zinb'`, `batch_size=128` | ✅ D | scvi-tools 官方文档（本次现取；注意 `n_layers` 默认=**1**） |
| scANVI | 继承 scVI；`n_labels_per_class=5`，`unlabeled_category` 需用户指定 | 🟡 D | scvi-tools 文档 / Xu 2021, Mol Syst Biol, 10.15252/msb.20209620 |
| **scDblFinder** | `dbr=NULL`（按 10x 期望率估计）、`dbr.sd=NULL`、`aggregateFeatures=FALSE` | ✅ D | Bioconductor 参考（本次现取；**注**：`dbr.sd` 默认为 NULL，非 0.015） |
| Harmony | `theta=2`, `lambda=1`, `max_iter=10`, `nclust=NULL→min(round(N/30),100)` | 🟡 D | Korsunsky 2019, Nat Methods, 10.1038/s41592-019-0619-0 |
| Seurat v5 | `FindVariableFeatures nfeatures=2000, vst`；`RunPCA npcs=50`；`k.param=20`；`resolution=0.8` | 🟡 D | Hao 2024, Nat Biotechnol, 10.1038/s41587-023-01767-y |
| **CopyKAT** | `ngene.chr=5`, `win.size=25`, `KS.cut=0.1`, `distance='euclidean'` | 🟡 D | Gao 2021, Nat Biotechnol, 10.1038/s41587-020-00795-2 |
| inferCNV | `cutoff=1`(10x)、`HMM_type='i6'`、`denoise=TRUE`、`window_length=101` | 🟡 D | Patel 2014, Science, 10.1126/science.1254257 |

---

## 阶段 2–8（🟡 已附出处，值待逐条现取）
| 工具 | 参数 | 出处 |
| :--- | :--- | :--- |
| RCTD | `doublet_mode='full'`, `max_cores`, `UMI_min`, `CELL_MIN_INSTANCE` | Cable 2021, Nat Biotechnol, 10.1038/s41587-021-00830-w |
| PLIP | 输入 224×224；ViT-B/32；`vinid/plip` | Huang 2023, Nat Med, 10.1038/s41591-023-02504-3 |
| SpaGCN | `p`, `l`, `resolution`（无固定默认，需搜索）, `refine=True` | Hu 2021, Nat Methods, 10.1038/s41592-021-01255-8 |
| **Squidpy** | `nhood_enrichment n_perms=1000`（✅）；`spatial_neighbors n_neighs=6, n_rings=1` | ✅/🟡 D — Palla 2022, Nat Methods, 10.1038/s41592-021-01358-2 |
| DESeq2 | `fitType='parametric'`, `test='Wald'`, `design=~...`, `alpha=0.1` | Love 2014, Genome Biol, 10.1186/s13059-014-0550-8 |
| CellRank 2 | `compute_lineage_drivers n_perms=1000`；GPCCA `n_states` 无固定默认 | Weiler 2024, Nat Methods, 10.1038/s41592-024-02303-9 |
| PAGA | `n_neighbors`（scanpy 默认 15）；`thresholds` | Wolf 2019, Genome Biol, 10.1186/s13059-019-1663-x |
| **AutoDock Vina** | `exhaustiveness=8`, `num_modes=9`, `energy_range=3.0`（✅ 源码）；网格盒按位点设 ~20–30 Å | ✅/D — Trott & Olson 2010, J Comput Chem, 10.1002/jcc.21334 |
| fpocket | 无强制 Dscore/体积默认；常用 `Dscore>0.5` | Le Guilloux 2009, BMC Bioinformatics, 10.1186/1471-2105-10-168 |
| GROMACS | 2 fs 步长 + LINCS；V-rescale 300 K；PME；力场 AMBER99SB-ILDN/CHARMM36 | van der Spoel 2005, JCC, 10.1002/jcc.20291 |
| LINCS CMap | Tau 定义（τ=Σ(1−xᵢ/max x)/(n−1)，0–1）；连通性 NCS | Subramanian 2017, Cell, 10.1016/j.cell.2017.10.049 |
| scanpy | `neighbors n_neighbors=15`；`umap min_dist=0.5, spread=1.0`；HVG `flavor='seurat'` | Wolf 2018, Genome Biol, 10.1186/s13059-017-1382-0 |

---

## 规则
1. **只有 ✅ 与 🟡(已附出处) 的才可入白皮书**；⚠️ 项如实标"本项目约定"。
2. **🟡 项在写入白皮书前逐条现取确认**（本次已抓到并纠正两处错误：scVI `n_layers`=1、scDblFinder `dbr.sd`=NULL）。
3. 凡未能核实者，**不写数值**，只写"按工具默认（见出处）"。
