# 🔬 科研严谨性与结果审计准则 (Scientific Rigor & Results Audit Protocol)

## 🎯 核心使命与职责定位
本准则定义了在本项目（肺腺癌侵袭与微环境演进多组学课题）中必须坚守的**最高科研诚信、生物学自洽性与数据真实性审计标准**。所有智能体（Agents）、分析脚本与报告撰写必须 100% 遵守以下铁律。

---

## 🛑 第一法则：伪造与模拟数据零容忍（Zero-Tolerance for Synthetic/Dummy Data）

1. **严禁任何非生物学随机生成**：
   - 严禁在任何分析或可视化脚本中使用 `np.random.uniform`, `np.random.normal`, `runif()`, `rnorm()`, `sample()`, `make_blobs()`, 或硬编码虚拟数值来伪造细胞比例、基因表达或空间丰度；
   - 严禁使用任何形式的虚构数据伪造 51 切片或空间微环境演进表。
2. **单一数据真实性溯源**：
   - 每一个下游图表和统计检验，其输入数据必须来自真实的原始测序矩阵（如 `*.mtx.gz`, `*.rds`, `*.h5ad`, `*.tsv.gz`），且具有明确的样本来源（GSE131907, GSE189357, GSE148071, GSE189487, GSE307534, TCGA-LUAD, Zenodo 15104582 等）。
3. **被审计出的假脚本与假数据必须立即物理销毁**：
   - 一旦发现未经真实生物学计算的废弃脚本，必须立即彻底清除，禁止在生产环境中保留。

---

## 🔍 第二法则：生物学逻辑与特征基因特异性审计（Biological Marker Specificity Audit）

在进行任何单细胞分群、空间解卷积或微环境注释时，必须严格通过经典特征基因的特异性交叉校验：

| 谱系类别 | 必须表达的经典阳性标志物（Canonical Markers） | 必须为阴性或极低的排他标志物 | 生物学自洽性底线 |
| :--- | :--- | :--- | :--- |
| **恶性/正常上皮** | *EPCAM*, *KRT7*, *KRT19*, *NKX2-1*, *SFTPC*, *SFTPB*, *AGER*, *CAV1* | *PTPRC (CD45)*, *CD3D*, *PECAM1* | 严禁与免疫细胞/内皮细胞混杂 |
| **T & NK 淋巴细胞** | *CD3D*, *CD3E*, *CD4*, *CD8A*, *GZMB*, *PRF1*, *PDCD1*, *HAVCR2*, *NKG7* | *EPCAM*, *COL1A1*, *CD68* | 严禁与上皮/髓系细胞重叠 |
| **B & 浆细胞** | *CD19*, *MS4A1 (CD20)*, *CD79A*, *SDC1 (CD138)*, *IGHG1*, *IGHA1* | *CD3D*, *EPCAM*, *ACTA2* | 严格区分 Naive B 与分化后 Plasma |
| **髓系巨噬细胞** | *CD68*, *CD163*, *C1QA*, *C1QB*, *C1QC*, *SPP1*, *MARCO*, *CLEC9A* | *CD3D*, *EPCAM*, *PECAM1* | 严格区分 C1QC+ 固有 TAM 与 SPP1+ 促侵袭 TAM |
| **成纤维基质细胞** | *COL1A1*, *COL1A2*, *ACTA2 (α-SMA)*, *PDGFRB*, *FAP*, *CXCL12* | *PTPRC*, *EPCAM*, *PECAM1* | 严格区分 myCAF（促纤维化）与 iCAF（炎症性） |
| **血管内皮细胞** | *PECAM1 (CD31)*, *VWF*, *CDH5*, *EGFL7*, *KDR* | *EPCAM*, *PTPRC*, *COL1A1* | 严格区分 Capillary 毛细与 Tip-like 肿瘤新生内皮 |

---

## 📐 第三法则：数学严谨性与参数透明度审计（Mathematical Rigor & Parameter Auditing）

1. **质控过滤阈值显式化**：
   - 必须显式声明过滤参数：500 <= nFeature_RNA <= 10,000, 1,000 <= nCount_RNA <= 60,000, percent.mt < 10.0%, 双细胞期望率 2.5%。
2. **批次矫正收敛性与高变基因**：
   - 必须记录高变基因数目（Top 3,000 HVGs, VST 法），PCA 主成分数（npcs = 30），Harmony 迭代参数（`theta=2`, `max.iter=20`）。
3. **空间流形解卷积数学公式**：
   - 必须采用 Seurat 概率流形转移算子：
     ```text
Score(s, t) = Σ_{k ∈ K} W_{sk} · [ (x_s · r_{tk}) / (||x_s||_2 * ||r_{tk}||_2) ]
```
     其中截断阈值 tau = 0.05 必须固定，严禁人为随意调整以掩盖低置信度噪点。
4. **空间富集比 Ro/e 严格计算**：
   - 必须基于真实频数列联表计算 Ro/e = (N_ij * N_total) / (N_i * N_j)，并在 chi-square 检验 P < 0.05 下判定统计显著性。
5. **生存分析统计学标准**：
   - Cox 回归必须输出风险比（HR）、95% 置信区间（95% CI）与对数秩检验（Log-rank test）精准 P 值，严禁篡改显著性。

---

## 🔬 第四法则：多组学跨平台交叉校验（Cross-Modal Verification）

1. **单细胞 ↔ 空间 Visium ↔ 亚细胞 CosMx 互证**：
   - 单细胞中发现的 14 个恶性克隆与 16 类 TME 亚型，必须在 14 张 Visium 空间切片中能观察到解剖物理分层；
   - 在 2 张 CosMx 亚细胞切片（~0.18 µm/pixel）中，必须能在肿瘤侵袭边界（Invasive Margin）清晰观察到 **myCAF 胶原物理屏障** 与 **SPP1+ TAM 促侵袭结界** 对 **CD8+ Tex 耗竭 T 细胞** 的局部排斥现象。
2. **单细胞/空间 ↔ 全外显子 WES-TMB 突变反逻辑互证（当前不可信，待重算）**：
   - ⚠️ 现有 TMB/THPP 结果**不可用**：`pipeline/11_tmb_dynamics_and_paradox_signature.py` 写死了 TMB 表、HR/P 值与 `np.exp` 伪 KM 曲线，**未做任何 Cox/生存计算**。
   - 既往文档中的"138 个 THPP 基因（TTN/MUC16/…）"与代码中的"16 个"亦互相矛盾。**在接入真实 TCGA 数据重算前，禁止引用任何 TMB/THPP 数值。**

---

## 🛑 第零法则：数据身份与分期真值（Data Identity & Stage Truth）
1. **数据集身份以 GEO 为准**：GSE131907（LUAD，含正常/原发/转移/脑转/胸水）、GSE189357（3 AIS + 3 MIA + 3 IAC，**无 AAH**）、GSE148071（**Advanced NSCLC**）。
2. **分期映射严禁静默默认**；GSE189357 以逐样本 histology 为准（TD1/2/9=IAC、TD3/4/6=MIA、TD5/7/8=AIS）。
3. **恶性标签必须经 CNV 推断证真**，禁止仅凭泛上皮标记（EPCAM/KRT…）判定。
4. **双体去除必须用标准算法**（DoubletFinder/scDblFinder），禁止以 nCount+nFeature 排序替代。
5. 详见 [`PROJECT_SUMMARY.md`](../../PROJECT_SUMMARY.md)。

---

## 🛡️ 审计执行流程（Audit Workflow）
在将任何分析结果或图表正式纳入科研成果报告（`walkthrough.md` 或论文初稿）前，审计 Agent 必须执行以下四步检查：
1. **数据源头核验**：确认输入矩阵为真实测序数据且路径有效；
2. **细胞数与标签一致性核验**：确认细胞总数（如 True_V5 的 328,028 个细胞）与各亚群计数精确无误；
3. **图表代码重现性核验**：确认图表由确定性脚本渲染生成，代码逻辑无任何随机模拟；
4. **方法学与顶刊标准比对**：对照 *Cell Reports Medicine* 2026 顶刊标准确认参数规范。
