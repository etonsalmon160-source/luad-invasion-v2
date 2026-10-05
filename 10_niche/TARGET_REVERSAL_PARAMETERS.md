# 靶点扰动逆向臂 · 参数与出处（审稿人可复现版）

**用途**：本臂的 WTCS / NCS / τ 是**我们自己实现的**，因此每个参数都必须能追溯到出处。
本文件逐参数给出 **取值 / 出处 / 档位**，档位分三档：

| 档位 | 含义 |
|---|---|
| **P** | **论文原文**（Subramanian et al. 2017 *Cell*，CMap 原版）|
| **T** | **工具惯例/默认**——可引用，但**不是**论文规定，**不得**写成"CMap 标准" |
| **S** | **本项目自定**——是我们的选择，不是领域标准，须在报告中明示 |

---

## 1. 数据

| 资产 | 取值 | 出处 |
|---|---|---|
| 主参照库 | `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx`，**12,328 基因 × 118,050 签名**，Z 分差异谱 | NCBI GEO；完整性经 **官方 SHA512SUMS 校验通过** |
| 扩展参照库 | `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`，473,647 签名 | 同上（该集无官方 SHA512，仅核体积） |
| 元数据 | `sig_info`（sig_id→pert_id/cell_id/pert_type）、`gene_info`（Entrez↔symbol） | GEO |
| **Touchstone 标记** | `signatureSearch`（Bioc 3.16 / 1.12.0）随包的 `data/lincs_pert_info2.rda`，**逐字含列 `is_touchstone`** | **P**（CMap 2017 的 Q_ref 定义需要它） |
| 第二套库 L2 | `xingjiepan/SCMG_data` 的 `pseudo_bulk_perturbation_database.h5ad`，20,345 扰动 × 18,108 基因，**差异位移谱**（实测均值≈0、40–60% 为负）；`is_touchstone` 不适用（非化合物） | 数据集卡 **license: mit**，`gated: False` |

**关键映射（实测确认）**：gctx **列名 = `sig_info$sig_id`**；gctx **行 = `gene_info$pr_gene_id`（Entrez）**。
⇒ 查询签名（基因符号）经 `pr_gene_symbol → pr_gene_id` 转 Entrez。

---

## 2. 查询签名（§5.2）

| 参数 | 取值 | 档位 | 出处 / 说明 |
|---|---|---|---|
| 签名构造 | `log2((mean_A(g)+ε)/(mean_rest(g)+ε))` | **S** | 域 vs 其余域；ε 见下 |
| **ε** | **1e-4** | **S** | 防 log 0；**无数值出处，纯本项目选择** |
| `q_up` / `q_down` 规模 | **各 150** | **T** | `signatureSearch::rand_query_ES(set_length=150)` 的**软件默认**；**CMap 2017 正文未规定** |
| 切分方式 | 按 `log2fc` 降序 top150 / 升序 top150 | **P** | CMap 2017：「Each gene in the query carries a **sign**…」——查询是**带符号基因集** |
| 主查询 | **D3**（缺氧侵袭）；D5 / D7 为反向对照 | **S** | 由本臂 §20/§22 定出 |

---

## 3. WTCS（§5.3）

| 参数 | 取值 | 档位 | 出处 |
|---|---|---|---|
| 统计量 | 加权 Kolmogorov–Smirnov 富集分数 | **P** | CMap 2017 逐字：「based on the **weighted Kolmogorov-Smirnov enrichment statistic (ES)**」 |
| **权重指数 `p`（代码里的 `type`）** | **1** | **T** | CMap 2017 只说"weighted"，**未给 p**；`p=1` 是 **GSEA (Subramanian 2005) 默认** |
| 权重取 | **\|值\|**（绝对值） | **P** | GSEA 定义；**实测漏取 abs 会致 `Nr` 为负 ⇒ 全 NA**（本项目踩过两次） |
| 合成 | `(ES_up − ES_down)/2` **当且仅当** `sign(ES_up) ≠ sign(ES_down)`，否则 0 | **P** | CMap 2017 逐字 |
| 单臂情形 | 本项目恒为双臂 ⇒ 不触发 | — | — |

---

## 4. NCS（§5.4）

| 参数 | 取值 | 档位 | 出处 |
|---|---|---|---|
| 分组键 | **`cell_id` × `pert_type`** | **P** | CMap 2017：「normalize the values **within each cell line and perturbagen type**」 |
| 正/负分取均值 | `μ⁺` = 正 WTCS 均值；`μ⁻` = 负 WTCS 均值 | **P** | 同上 |
| **`μ⁻` 取绝对值** | `NCS = w / abs(μ⁻)` | **P** | 官方源码同处注释：「**without abs() sign of neg values would switch to pos**」 |
| 零 WTCS | **从均值计算中剔除** | **P** | 官方一致：`es_na[es_na==0] <- NA` |
| 分母为 0 | 置 `10^-12` | **T** | 官方同（避免除零） |
| L2 的分组 | 用 `dataset`（L2 无 `cell_id`） | **S** | 本项目适配；须在报告中标注 |

---

## 5. τ（§5.5）

| 参数 | 取值 | 档位 | 出处 |
|---|---|---|---|
| 公式 | `τ = sgn(NCS) × 100/N × Σᵢ[ \|NCSᵢ,r\| < \|NCS_q,r\| ]` | **P** | CMap 2017 逐字 |
| `Q_ref` 定义 | 同 **cell line × perturbation type** 的 Touchstone 签名 | **P** | CMap 2017：「**exemplar signatures of Touchstone perturbagens** that match the cell line and perturbation type」 |
| Touchstone 判定 | `lincs_pert_info2$is_touchstone` | **P/数据** | 见 §1 |
| **`exemplar` 的操作定义** | 每个 (细胞系, pert_id) 取 **`distil_ss` 最大**的签名（同分取 sig_id 字典序最小），要求 `distil_nsample ≥ 3` | **S** | 概念出自 **P**（CMap 2017「exemplar signatures」）；但**官方文档未给操作定义** ⇒ 本项目操作化。指标 `distil_ss`/`distil_nsample` 取自**官方随发的 `sig_metrics`** |
| exemplar 自洽性 | 得出的细胞系 = **A375/MCF7/HA1E/HT29/PC3/YAPC/HELA/A549/HCC515/HEPG2…**，**与 CLUE Touchstone 核心细胞系重合** ✓ | — | 规则未跑偏 |
| exemplar 规模 | A375__trt_cp **685**；共 **6,059**（29 个细胞系） | — | — |
| 敏感性① | 改用 `distil_cc_q75` 最大者 | **S** | 检验代表**指标**的影响 |
| 敏感性② | **不折叠**（该组全部 touchstone 签名，nsample≥3） | **S** | 检验"每扰动一个代表"**这一步**的影响 |
| 阈值 | **τ ≤ −90** 主档；−95 / −98 敏感性档 | **P** | **CLUE 官方文档逐字**：「we consider **tau of +90 or higher, and of −90 or lower**, as strong scores…」（Connectopedia `connectivity_scores`）；CMap 2017 结果叙述同 |

---

## 6. 零模型与报告规则

| 项 | 取值 | 档位 |
|---|---|---|
| 随机查询零模型 | **N = 1000** | **T**（`rand_query_ES` 默认） |
| 随机基因匹配 | 按表达分位 | **S** |
| 出表硬门 | ① 零模型跑完；② ~~L2 许可证~~ ✅ 已解（MIT） | — |
| 措辞 | 只准「**候选扰动**」；禁 `causal`/`driver`/`validated target`/「状态逆转因子」/「药」 | **P**（`docs/scientific_rigor_and_audit.md:64`） |

---

## 7. 尺子的独立对拍（**本臂最关键的一步**）

| 项 | 结果 |
|---|---|
| 归因 | WTCS / NCS / τ 均为**自实现** ⇒ 必须独立复算 |
| 参照实现 | `signatureSearch` 1.12.0（Girke 实验室；Bioconductor 官方 CMap/LINCS 实现）的 `.enrichScore()` / `.lincsScores()` |
| 环境 | 隔离库 `~/Rlibs/sigsearch`（`BiocManager::install(lib=…)`，R 4.2.2 / Bioc 3.16） |
| 对拍脚本 | `25_tr_dump_for_crosscheck.py`（导出子集+自算值）→ `26_tr_crosscheck_official.R`（官方复算+比对） |
| **结果（2,000 签名）** | WTCS 最大绝对差 **6.6e-14**、NCS **1.2e-13**（皆机器精度）；**符号 2000/2000 一致**；分布逐位相同 |
| 附：ES 的第三方对拍 | 另与 **`fgsea` 1.24.0** 对拍（2,000 基因 × 30 集）：最大差 **1.11e-16** |
| 🔴 **对拍逮到的自实现 bug** | ① NCS 的 `μ⁻` 未取绝对值 ⇒ NCS 全正 ⇒ **τ 全正 ⇒ 候选恒空（静默假阴性）**；② 向量化 ES 漏 `np.abs(C)` ⇒ `Nr` 为负 ⇒ 全 NA。**两个都已修并复验** |
| 复现命令 | 见 §8 |

---

## 8. 复现命令（逐条照抄即可）

```bash
R_LIBS_SS=/home/eto/Rlibs/sigsearch:/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library

# 1) 安装官方参照实现（隔离库）
R_LIBS="$R_LIBS_SS" Rscript -e '.libPaths(c("/home/eto/Rlibs/sigsearch",.libPaths())); \
  BiocManager::install(c("signatureSearch","signatureSearchData"), lib="/home/eto/Rlibs/sigsearch", ask=FALSE)'

# 2) 导出对拍子集（自算 WTCS/NCS）
python3 10_niche/25_tr_dump_for_crosscheck.py --ncol 2000 --query D3 \
        --outdir /home/eto/luad_v2/results/10_niche/target_reversal/xcheck

# 3) 官方复算 + 比对
R_LIBS="$R_LIBS_SS" Rscript 10_niche/26_tr_crosscheck_official.R \
        /home/eto/luad_v2/results/10_niche/target_reversal/xcheck
# ⇒ 期望："✅ 自实现与官方实现逐位一致，口径无误。"
```

---

## 9. 已知局限（须与结果同页）

1. **τ 的参照集口径未定**（§5 第 4 行）：官方把参照烘死在它自己的库上，我们的"全 touchstone"构造与之**不等价**。
2. `pert_id` 中约 **37%** 无 `is_touchstone` 取值 ⇒ 这些签名的 touchstone 状态不详（但 Q_ref 靠**组级**参照，影响有限）。
3. 我们**不是** CLUE 的 **TAS**。
4. L2 全是细胞系、**无肺**，与本项目组织背景有域差。
5. 无扰动实验验证 —— 本臂只出**计算排序**。
6. 自定参数（ε / 表达分位匹配 / L2 分组 / exemplar 落法）须在报告中逐条标明为**本项目选择**。
