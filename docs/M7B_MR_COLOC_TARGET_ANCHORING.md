# M7b · 靶点池的遗传统计锚定（cis-MR + coloc）操作手册

> **定位**：本文件是 [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md) **M7b** 的操作手册。
> 作用：把从单细胞/空间得到的**候选基因**，用**人类遗传学**证据升级为
> **「遗传学支持的候选靶点（genetically supported candidate target）」**——这是**无需湿实验**能达到的最高诚实台阶，
> 并且**直接喂给 M8 的分子对接**。
>
> **状态**：方法学已核实（2026-09-12）；**尚未运行**（属 M-1 之后的 M7）。
> 本文件描述方法、数据入口、判定阈值与产出 schema；**所有数值必须来自真实数据运行**（铁律 R3）。

---

## 1. 它解决什么问题（为什么需要）

单细胞/空间给出的是**关联**：某基因在 MIA/IAC 高表达。关联 ≠ 因果，三个陷阱：
**混杂**（吸烟/炎症带动）、**反向因果**（癌导致高表达）、**共调控**（一个模块齐变，分不清真凶）。

**MR（孟德尔随机化）** 用**受精时固定、随机分配**的胚系变异作工具变量，绕开这三条——
相当于一场"天然的随机对照试验"。**cis-MR + coloc** 是药物靶点验证的行业标准
（PCSK9 / IL6R / TYK2 / FXI 皆由此确认）；遗传学支持的靶点，临床成功率约翻倍。

**一句话**：MR 给**因果方向**，单细胞/空间给**细胞与位点**，两者互补。

---

## 2. 产出定义

对**每一个候选基因 G**（来自 M7a 生态位签名 与 M7 恶性程序/regulon），产出一行：

| 列 | 含义 |
| :--- | :--- |
| `gene` / `ensembl_id` | 候选基因 |
| `source` | 来源：`malignant_program` / `niche_signature` |
| `spatial_localization` | 单细胞/空间定位（如 "invasion front", "myCAF niche"） |
| `n_cis_snp` | 用作工具的 cis-SNP 数 |
| `mr_method` | `Wald`（单 SNP）/ `IVW`（多 SNP） |
| `mr_beta` / `mr_se` / `mr_or` / `mr_or_ci_low` / `mr_or_ci_high` / `mr_p` | MR 估计（OR 为**每 1 SD 表达**的风险比） |
| `mr_fdr` | BH 校正后 FDR |
| `f_stat_min` | 最小工具 F 统计量（弱工具筛查） |
| `egger_intercept_p` | 多效性检验（多 SNP 时） |
| `steiger_ok` | Steiger 方向是否正确 |
| `coloc_pph4` | 共享因果变异后验概率 |
| `coloc_verdict` | `shared`(H4) / `distinct`(H3) / `inconclusive` |
| `ot_tractability` | Open Targets 可成药性分级 |
| `known_drug` | 是否有临床期药物（ChEMBL/Open Targets） |
| `genetic_support` | **最终判定**：`supported` / `not_supported` / `not_testable` |

---

## 3. 原理速览

### 3.1 MR 三条假设（每条都要查）
1. **Relevance**：工具变量确实影响暴露（cis-eQTL 显著）→ 用 **F 统计量 > 10** 筛查。
2. **Independence**：工具变量与混杂无关 → cis 变异天然较接近此条。
3. **Exclusion restriction**：**只通过该基因**影响疾病 → 最脆的一条，靠**敏感性分析**兜底：
   - **MR-Egger 截距**（Bowden 2015）：截距显著≠0 → 存在方向性多效性；
   - **加权中位数 / 加权众数**（Bowden 2016）：对部分无效工具稳健；
   - **Cochran's Q 异质性**；
   - **Steiger 过滤**（Hemani 2017）：确保方向是"变异→表达→疾病"，而非反向。

### 3.2 cis-MR 与 coloc 为什么合起来用
- **cis-MR**（只用基因 ±1 Mb 内的变异作工具）：生物学最可信、多效性最少 → 给**方向与效应量**。
- **coloc**（Giambartolomei 2014；Wallace 2021）：贝叶斯检验 GWAS 信号与 eQTL 信号**是否同一个因果变异**
  → 排除"两个不同变异恰好在 LD 中"的假象。输出 **PP.H4**（共享）与 **PP.H3**（独立）。

### 3.3 两个合起来判
- MR 显著 **且** `PP.H4 ≥ 0.8` → `supported`（遗传学支持）
- MR 显著但 `PP.H3` 高 → `distinct`（该位点对疾病的作用不经该基因）→ `not_supported`
- 取不到 cis-SNP → `not_testable`（**如实标缺，不做替代**）

---

## 4. 数据来源（可达性已于 2026-09-12 实测）

| 角色 | 来源 | 入口 | 实测 | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| **暴露·cis-eQTL（首选）** | **eQTL Catalogue** | https://www.ebi.ac.uk/eqtl/ （直连 200） | ✅ | 跨组织**统一格式**的 cis-eQTL，含肺 |
| 暴露·cis-eQTL（备选） | GTEx v8 肺 | GTEx Portal / dbGaP phs000424 | ⚠️ 门户本机不通 | 可用 eQTL Catalogue 的 GTEx 处理版替代 |
| 暴露·cis-eQTL（血，大样本） | eQTLGen | https://www.eqtlgen.org/ （直连 200） | ✅ | n=31,684；组织不同，须注明 |
| **结局·GWAS（实操首选）** | **IEU OpenGWAS** | https://gwas.mrcieu.ac.uk/ · API `api.opengwas.io`（200） | ✅ | 需**免费 token**；TwoSampleMR 直接对接 |
| 结局·GWAS（权威） | ILCCO / TRICL 肺癌 | McKay 2017 *Nat Genet*；Byun 2022 *Nat Genet* | ⚠️ | GWAS Catalog 中不少研究 `fullPvalueSet=False`；完整汇总统计可能受控（dbGaP）→ **须先确认可得性** |
| LD 参考 | 1000 Genomes EUR | 1000G / OpenGWAS 内置 | ✅ | 祖先必须与 GWAS 匹配 |
| 可成药性 | Open Targets Platform | GraphQL API | ✅ | tractability buckets |
| 依赖/重定位 | DepMap、DGIdb、ChEMBL | 公开 | ✅ | 选择性依赖，非泛必需 |

> ⚠️ **最大不确定性**：ILCCO/TRICL 的**完整** LUAD GWAS 汇总统计可能需要申请。
> 若不可得 → 改走 **IEU OpenGWAS** 中的肺癌 GWAS，并在报告中**明确来源与样本量**；不得用别的表型顶替。

---

## 5. 环境

```bash
# R 4.2.2 已装。补齐 MR/coloc 依赖：
Rscript -e 'install.packages(c("remotes","data.table","ggplot2"))'
Rscript -e 'remotes::install_github("MRCIEU/TwoSampleMR")'
Rscript -e 'remotes::install_github("MRCIEU/ieugwasr")'
Rscript -e 'install.packages("coloc")'
# 版本记录（产物须带）：sessionInfo() 写入 run 日志
```
- `TwoSampleMR`（Hemani 2018 *eLife*）· `ieugwasr`（OpenGWAS 客户端）· `coloc`（Giambartolomei 2014；v5 Wallace 2021）
- OpenGWAS token：`ieugwasr::get_opengwas_jwt()`（首次需在 https://api.opengwas.io 注册）

---

## 6. 操作步骤

> 下为**结构骨架**，真实运行须记录输入哈希、版本、种子（项目铁律 R4/R5）。

```r
# 6.0 输入：候选基因表（来自 M7a/M7）
cand <- read.csv("results/M7/candidate_genes.csv")   # gene, source, spatial_localization

for (g in cand$gene) {
  # 6.1 cis-eQTL（肺）：取 ±1Mb 内显著 SNP（eQTL Catalogue / GTEx 处理版）
  #     → data.frame(snp, beta_exposure, se_exposure, ea, nea, p, n)
  expo <- fetch_cis_eqtl(g, tissue = "lung")

  # 6.2 结局 GWAS：取同一批 SNP 的效应（OpenGWAS / ILCCO）
  out <- fetch_outcome_gwas(expo$snp, gwas_id = LUAD_GWAS_ID)

  # 6.3 协调（等位基因对齐，避免链翻转）
  dat <- harmonise_data(expo, out)

  # 6.4 MR：单 SNP→Wald；多 SNP→IVW；敏感性分析
  res    <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  pleio  <- mr_pleiotropy_test(dat)      # Egger 截距
  het    <- mr_heterogeneity(dat)        # Cochran's Q
  steig  <- directionality_test(dat)     # Steiger
  fstat  <- min(dat$beta.exposure^2 / dat$se.exposure^2)   # 弱工具

  # 6.5 coloc：同一区域（eQTL 与 GWAS）
  col <- coloc::coloc.abf(dataset1 = eqtl_region, dataset2 = gwas_region)
  pph4 <- col$summary["PP.H4.abf"]; pph3 <- col$summary["PP.H3.abf"]

  # 6.6 判定（见 §7），写入一行
}
```

### 关键质控点
- **等位基因协调**：`harmonise_data(action=2)`，剔除模糊/回文 SNP；
- **弱工具**：`F > 10`，否则剔除并注明；
- **LD/祖先**：LD 参考与 GWAS 人群一致（EUR）；
- **样本重叠**：暴露与结局样本不应重叠（eQTL 与 GWAS 人群不同 → 通常满足）；
- **吸烟混杂**：报告是否用吸烟校正的 GWAS；不校正须在局限中声明。

---

## 7. 判定阈值

| 指标 | 阈值 | 出处/性质 |
| :--- | :--- | :--- |
| cis-SNP 数 | ≥ 1（单 SNP 用 Wald） | 方法学 |
| 工具 F 统计量 | **> 10** | 弱工具常规（Burgess & Thompson 2011） |
| MR P 值 | **BH-FDR < 0.05** | 项目约定（多重检验） |
| MR-Egger 截距 | P > 0.05（无方向性多效性） | Bowden 2015 |
| Steiger | 方向正确 | Hemani 2017 |
| coloc PP.H4 | **≥ 0.8** | 常规（Giambartolomei 2014） |
| coloc PP.H3 | 高 → 判 `distinct` | 常规 |

> ⚠️ 阈值为**常规/项目约定**，须在 `docs/PARAMETERS_AND_SOURCES.md` 登记；**不得事后调阈值**。

---

## 8. 局限与红线（措辞规范）

- ✅ 可说：**"genetically supported candidate target"**、**"prioritized target"**、**"repurposing candidate"**
- ❌ 不可说：`causal`（无条件）、`driver`、`drug target`、`validated`
- MR 是**群体水平**因果，**不等于**"该基因在某细胞里驱动了 MIA"——细胞层面结论仍来自单细胞/空间；
- 肺 eQTL 样本量有限 → 部分基因 `not_testable`，**如实标缺，不替代**；
- 目标物**结构对接**（M8）仍是**计算假说**，`validated` 须湿实验（SPR/ITC/类器官）。

---

## 9. 复现与产物要求（铁律）

1. 输入（eQTL 表、GWAS 汇总统计、候选基因表）**均记 SHA-256**；
2. 输出表：`results/M7b/mr_coloc_target_anchoring.csv` + 运行日志（`sessionInfo()` + 参数 + 种子）；
3. 每基因的中间产物（harmonised 表、coloc 区域）留存以便复核；
4. 全部脚本确定性、无 `np.random`/`runif`、无写死数值。

---

## 10. 待办（M-1 之后执行）

- [ ] 确认 LUAD GWAS 可得性：OpenGWAS 取一个可用 ID **或** 取得 ILCCO/TRICL 汇总统计（记录来源与 n）
- [ ] 下载 eQTL Catalogue 肺 cis-eQTL（记录版本）
- [ ] 注册 OpenGWAS token
- [ ] 把 §7 阈值登记进 `docs/PARAMETERS_AND_SOURCES.md`
- [ ] 用 M7a 的候选基因跑通 1 个基因（smoke test）→ 再全量
