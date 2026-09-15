# 🛑 空间转录组 6 阶段切片映射与图表严谨性最高禁令 (Spatial Cohort & Figure Prohibitions)

> **生效级别**：最高科研诚信铁律（P0 级永久约束），全流程任何脚本、报告与智能体行为必须无条件遵守，违者即视为学术不端。

---

## 一、 6 阶段空间切片映射铁律 (Golden Standard Slice Manifest)

在整个项目所有生成 Figure 1、空间矩阵、微环境解卷积及多模态图谱的脚本中，**切片与分期映射关系必须 100% 锁定如下，严禁任何形式的擅自替换或回退**：

| 阶段 (Stage) | 阶段标准名称 (Label) | 唯一合法切片 ID | 权威数据库 | GSM / 样本编号 | 组织学特征与严禁事项 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Stage 1** | `Stage 1: Normal in PT` | **`P4_Normal`** | GSE307534 | GSM9226174 | 正常肺实质对照，严禁使用任何含癌切片。 |
| **Stage 2** | `Stage 2: AAH Precursor` | **`P1_AAH`** | GSE307534 | GSM9226168 | 真实不典型腺瘤样增生（AAH）前体切片。 |
| **Stage 3** | `Stage 3: AIS In Situ` | **`P3_AIS`** | GSE307534 | GSM9226172 | 原位腺癌（AIS）纯伏壁样（Lepidic）生长切片。 |
| **Stage 4** | `Stage 4: MIA Micro-inv` | **`P10_MIA`** | GSE307534 | GSM9226189 | **【铁律禁令】必须使用真实微浸润 `P10_MIA`；绝对禁止使用 `P4_AAH1` / `P4_AAH` 充当 MIA！** |
| **Stage 5** | `Stage 5: IAC Invasive` | **`P3_LUAD`** | GSE307534 | GSM9226173 | **【铁律禁令】原发浸润性腺癌（IAC），实性/腺泡型与纤维化基质！** |
| **Stage 6** | `Stage 6: LNM Metastasis` | **⛔ 暂无合法切片（待真实 LUAD 数据）** | — | — | **【2026-09-12 勘误】原锁定的 `PT_3_LNM`(GSE190811) 经 GEO 核实为【乳腺癌】淋巴结转移（系列标题 "…breast cancer patients"），且其 GSM 号 GSM5732148 在该库并不存在（真实为 GSM5732357–2360）。该切片已【作废】，禁止用于任何 LUAD 产物；Stage 6 空间锚点暂缺，待项目方提供真实 LUAD 淋巴结转移空转数据后重锁。** |

> **⛔ 勘误记录（2026-09-12）**：本文件原将 `PT_3_LNM`(GSE190811, GSM5732148) 作为 Stage 6 永久锁定切片，
> 经 GEO 逐样本核实为**乳腺癌**数据。已作废。**LNM 阶段空间产物（解卷积 / 生态位 / PLIP / 对接）在其被替换前一律不得产出**；
> 期间空间图谱限定为 **Normal → AAH → AIS → MIA → IAC**（全部来自 GSE307534 单库单平台）。
> 单细胞层面的 LNM 仍可用 GSE131907 `mLN`（真转移淋巴结，44 患者中已 GEO 核实）。

---

## 二、 算法与图表渲染的三大红线禁令 (Rendering Prohibitions)

1. 🚫 **严禁人工取模或伪造条纹（Zero Synthetic Modulo Logic）**：
   - WHO 病理学分类（Lepidic, Acinar, Papillary, Solid, Fibrotic Stroma）与 TCGA 分子亚型（TRU, PP, PI）必须 100% 由真实的 RCTD 亚克隆及微环境解卷积权重计算得出。
   - **绝对禁止出现 `i % 2 == 0`、`i % 3 != 0` 或 `np.random` 等任何人造交替条纹代码！**

2. 🚫 **严禁在非肿瘤背景乱涂假阳性（Clean Non-tumor Background Rule）**：
   - 非肿瘤细胞（正常肺泡上皮、淋巴细胞、巨噬细胞、血管内皮等）在 WHO 与 TCGA 列中必须统一保持整洁的浅灰底色（`#E2E8F0`），与图例中的 `Non-tumor / NA` 严格一致，严禁给正常组织赋予恶性亚型颜色。

3. 🚫 **严禁 ROI 选框沦为纯黄色肿瘤块（Multilineage Frontier ROI Law）**：
   - 第 5 列的虚线选框与第 6 列的放大六边形，必须是**恶性肿瘤（黄/橙/红）、成纤维基质（深红 myCAF）、巨噬细胞（紫 TAM）、T细胞（蓝）及内皮细胞（粉）紧密接触重塑的“多系侵袭破口前沿（Multilineage Frontier）”**。
   - **绝对禁止将 ROI 中心定位于肿瘤实质内部的 100% 纯恶性细胞黄块！**
