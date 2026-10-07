# 正文表格索引（草稿）

> 原则：**能进表的就不进图**——参数、计数、对拍数值、阴性清单这些"读数字"的内容用表；
> 图只留趋势与空间格局。

## 落地状态（2026-10-07）

七张表已全部写进 LaTeX，`main.tex` 编译 16 页 0 错误、无未解析引用。

| 表 | LaTeX 文件 | 插入位置 |
|---|---|---|
| Table 1 队列与数据来源 | `tex/sections/_tab_cohort.tex` | Methods §Data |
| Table 2 参数与出处 (P/T/S) | `tex/sections/_tab_params.tex` | Methods §Software |
| Table 3 方法学对拍与阳性对照 | `tex/sections/_tab_validation.tex` | Results §"Analyses that did not work" 末 |
| Table 4 化合物候选（34） | `tex/sections/_tab_compounds.tex` | Results §"In-silico reversal" |
| Table 5 七个域类型身份 | `tex/sections/_tab_domains.tex` | Results §"Spatial deconvolution" |
| Table 6 阴性结果及可能原因（12 条） | `tex/sections/_tab_negatives.tex` | Results §"Analyses that did not work" |
| Table 7 局限性清单 | `tex/sections/04_limitations.tex` | 正文 §Limitations（**逐条排版，未做表**）|

---

## Table 1 — 队列与数据来源

| 项 | 值 |
|---|---|
| 单核 RNA（snRNA） | GSE308103 · **413,697 核 / 23 患者**（论文主体分析用 QC 后 **399,579 核**）|
| 空间转录组（Visium） | GSE307534 · **56 张切片 / 25 患者** |
| 平台 | Visium CytAssist **11 mm**，点距 **100 µm**，14,336 spot（128×112）|
| 两模态共享患者 | **23 例**（Visium 多出的 2 例无单细胞数据）|
| 期别划分 | Normal → AAH → AIS → MIA → IAC |
| 切片数/期别 | Normal **1** ｜ AAH 11 ｜ AIS 14 ｜ MIA 4 ｜ IAC 26 |

> ⚠️ Normal 仅 **1 张**、MIA 仅 **4 张** ⇒ 两端估计不可单独引用。

---

## Table 2 — 各分析臂的参数与出处（**本项目要求"逐条标 P/T/S"**）

> 图例：**P** = 论文逐字（有出处）｜**T** = 软件默认（非论文）｜**S** = 本项目自定

| 分析臂 | 参数 | 取值 | 档 | 出处 |
|---|---|---|---|---|
| 单细胞 QC | 合格核数 | 767,000 → 399,579 | S | M1 |
| 降维 | 分辨率 | 0.5–0.8 | P | 源论文 |
| 逐谱系亚聚类 | L2 亚型数 | **39** | S | GP8 |
| BANKSY 分域 | λ / k_geom / resolution | **0.2 / 18 / 0.5**（主档固定） | P+S | 官方推荐 λ=0.2（55 µm 域分割）|
| BANKSY 共识 | K\* | **7**（K=8 起稳定性 <0.814） | S | §17.3 |
| RCTD | 判据 | `argmax == 上皮` | S | 开跑前改签 |
| RCTD | 参考谱系数 | 6 谱系 / 39 亚型 | S | M5 |
| 空间 CNV | 引擎 | fastCNV | P | DOI 10.1186/s13073-026-01731-w |
| 空间 CNV | 臂 | D 臂（全切片非上皮作参考） | S | §18.1 |
| 域 CNV 比较 | 口径 | 同切片内配对（域内中位 − 同切片其余中位） | S | §18.1 |
| 查询签名 | 构造 | 患者内配对（IAC − 前驱），取患者间**中位** | S | §14.2 |
| 查询签名 | 预后质控 | 疾病↑需 HR>1；疾病↓需 HR<1 | S | §14.2c |
| 化合物库 | LINCS | 118,050 条记录 × 12,328 基因 / 1,826 化合物 | P | GSE70138 |
| 化合物打分 | 三口径 | CMap(τ) / Cor-Spearman / Cor-Pearson | P | CMap 2017 + `signatureSearch` |
| 基因库 | SCMG | 20,345 个基因扰动 | P | `xingjiepan/SCMG_data`（MIT）|
| 阳性对照 | 数据 | 官方教程 上胚层→中胚层 | P | SCMG 官方 |
| 零模型 | 置换数 | **1,000**（签名层）/ 200（MOA 层） | T+S | `signatureSearch` 默认 |
| 劈半 | 次数 | 4 | S | §N |

---

## Table 3 — 方法学对拍与阳性对照（**本项目最强的方法学证据**）

| 检验 | 对象 | 结果 |
|---|---|---|
| **WTCS 逐位对拍** | 自实现 vs `signatureSearch` 1.12.0 | 最大绝对差 **6.611e-14**；符号一致 **2000/2000** |
| NCS 对拍 | 同上 | 最大绝对差 **1.203e-13** |
| **SCMG `causal_score` 对拍** | 自实现 vs 作者 Python | Spearman **1.00000000**；符号一致 **100.0000%**（19,819 条）|
| **阳性对照** | 官方教程已知答案 | **TBXT 第 1 / 8,450**；MSGN1 2；MIXL1 3；POU5F1 7；SNAI1 9；EOMES 11；NANOG 12；EVX1 15 |
| **跨模态** | 单细胞轴 vs 空转轴 | Spearman **0.735**；top-100 重叠 **52**（随机 5.5）|
| **劈半稳定性** | 患者随机分半 ×4 | 观测 ρ 中位 **0.962**；置换零模型 **±0.045** |
| **三口径一致性** | CMap / Cor-Sp / Cor-Pr | Cor×Cor ρ=**0.950**（top500 重叠 433）；CMap×Cor ρ=**0.36–0.38** |
| **空间标定自洽** | 相邻 spot 中心距 = 100 µm | 反推 **1 px = 0.2513 µm**；整图宽 **10.9 mm**（= 11 mm 平台）|

---

## Table 4 — 化合物候选（34 个，按机制）

| 机制 | 个数 | 代表化合物 |
|---|---|---|
| **PI3K / mTOR** | **16** | torin-1 · torin-2 · AZD-2014 · AZD-8055 · GSK-2126458 · GDC-0980 · MLN-0128 · PKI-179 · NVP-BEZ235 · OSI-027 · WYE-125132 · PI-103 · GDC-0941 · taselisib · voxtalisib · GSK-2110183 |
| 细胞周期激酶 | 4 | XL-888 · dinaciclib · PF-03814735 · SB-939 |
| HSP90 | 3 | NVP-AUY922 · tanespimycin · AT-13387 |
| 蛋白酶体 | 2 | MG-132 · bortezomib |
| 其他 / 未注释 | 9 | BIIB-021 · NVP-TAE226 · TAK-285 · JW-7-24-1 · YM-155 · lacidipine · midostaurin · pentobarbital · okadaic acid |

---

## Table 5 — 七个生态位域类型的身份

| 域 | 现行名 | 依据 | 占 spot | 备注 |
|---|---|---|---|---|
| D1 | 气道上皮域 | Ciliated ×6.3、Goblet ×5.6 | 5% | |
| D2 | iCAF | IL6/CXCL8/CCL2/LIF/HAS1/HAS2 + 签名层 Δ+0.0093 | 11% | §20.2 修订 |
| D3 | ECM/间质富集型 | 配平后 up 侧转 COL1A1/COL3A1/FN1/CTHRC1/ASPH | **仅 IAC，17%** | **两次修订**：免疫浸润→缺氧侵袭→作废 |
| D4 | 肺泡壁/血气屏障域 | Capillary ×2.5、AT1 ×1.9 | **45%** | 背景域；Normal 里占 77% |
| D5 | AT2 上皮域 | AT2 = 0.480（×2.0） | 19% | 空间 CNV 显著偏高 |
| D6 | 血管/气道壁平滑肌域 | VSMC ×5.1、Artery ×2.3 | 8% | |
| D7 | 淋巴细胞聚集域 | B ×9.1、NKT ×4.7 | 4% | **五期皆有** |

---

## Table 6 — 阴性结果清单及其**可能原因**（要你确认写法）

| 阴性结果 | 观察到什么 | 可能的问题所在 |
|---|---|---|
| **染色体 CNV 在单细胞层失败** | 6 个工具族无一给出可用 ROC（唯一 0.5 是统计量退化）| ① 参考锚不可靠（同 3,890 个正常细胞三次跑出 0%/13.6%/49.8%）；② 逐细胞标签本身不可得（铁律 2）|
| **空间 CNV 队列级未建立** | 配平后 p=0.092、CI 下界=0 | 深度与组织密度不可分；**不是**"LUAD 无 CNV" |
| **域级 CNV 的 D5 细分不成立** | cor(cf, nUMI)=0.740；只按计数切复现 ~75% | 分数与 RNA 含量强耦合，**分不开"肿瘤多"与"RNA 多"** |
| **发育谱三条路线全否定** | CytoTRACE 方向与原文**相反**；WOT 谱机制废弃；运输读数仅地板弱过 | ① 测的是**分期轴**不是亚型轴（无 AIC/KAC 标签）；② CytoTRACE 高 = 检出基因多，也可读成"转录组更杂乱" |
| **H&E / PLIP 天花板 ~0.73** | 六臂×两套共 12 格全在 0.571–0.729 | 架构 181 µm/token 与最细彩色源 5.670 µm/px 的乘积；**不等于"细节无用"** |
| **CMC 提示词调优无效** | 1344 条候选只有 22 条(1.6%)超过"不算术语"的基线 | H&E 判恶性信号不在图里 |
| **基因层逆转无信号** | SCMG 官方因果分：官方对照 1.279 vs 我们轴 0.083 | 本队列疾病轴与库中扰动的表达位移**无可检出匹配**；**不是**"SCMG 找不到靶点" |
| **四程序随病程无单调趋势** | 配平前分数≈深度镜像；配平后 10 档顺序一致但与病程不符 | **深度混杂制造假趋势**（§13.1）|
| **生态位闸 1 FAIL** | 跨种子 ARI≥0.90 全档 0/40 | 生态位定义对方法敏感（BANKSY 7 / RCTD 6 / marker 5）|

---

## Table 7 — 局限性清单（正文第 5 节，逐条）

（现有 13 条，见 `PAPER_DRAFT_v2.md` §5；成稿时整理成表）
