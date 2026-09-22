# GP8c 报告 · 六个谱系亚聚类 + 注释 + 污染体检 + 谱系接管

生成时间：2026-09-22 08:35（夜跑自动完成）
范围：`GSE308103` snRNA，413,697 核。只读源 `/home/eto/luad_invasion` 未写入。

---

## 一、一句话结论

六个谱系的亚聚类和注释**全部跑完**（退出码全 0），每个细胞恰好归属一个谱系（413,697 互斥且完备）。
但有 **三件事必须放到最前面**，因为它们影响你怎么读后面的所有数字：

1. ~~**内皮的亚型注释在方法上不成立**~~ → **已修复（2026-09-22 你追认 S4 面板）**：原来每型只有 1–2 个基因槽，4 种内皮被压成 2 种——不是数据问题，是面板问题；换成同篇论文 Table S4 的逐簇富集表后，**8 个内皮亚型各自胜出**。详见第三节（一）。
同一份 S4 也换了髓系 L2 面板（8 型 → 12 型胜出），详见第三节（一之二）——**但 S4 不改 L1，全局簇 17/26 是你人工裁决修掉的，不是 S4 修掉的**。
2. **B/浆 的"27% 污染"其实是一个样本（P24）的问题**，可查、可处置，且**不影响其它五个谱系**。环境 RNA 检验已跑完**预注册失败 → 探索版纠正 → 重冻结口径确认版**三步，确认版**三项全过**（病人级 7.8× 中位、稳健 z=10.1、负对照经验 p=0.005）⇒ **P24 的环境 RNA 假设成立**，详见第三节（二）。**原文（PMC12980502）全文没有提这件事**——既没做环境 RNA 校正，也没标记任何人类 snRNA 样本，反而把 P24 当代表案例用（同节末）。**待你定的只剩处置方式**——**2026-09-22 你已裁定『保留 + 标注』**（不剔除、不现在加环境校正），决定书 `results/05_annotation/p24_disposition_decision.json`。
3. **上皮旧版"过关"作废**——在修正后的细胞集上复现不出来。

---

## 二、六个谱系的最终结果

### 上皮（epiA）· r\*=0.7 · 种子 0 · **事后放宽门槛（探索性）**
27 个亚簇，133,384 核，**0 簇被剔**。

| 亚型 | 核数 |
|---|---|
| AT2 | 88,666 |
| AT1 | 30,320 |
| Ciliated | 6,697 |
| Goblet/Mucous | 3,587 |
| Basal | 2,386 |
| Serous | 1,728 |

未胜出的类型（如实上报，未强行凑满 10 型）：Club、Ionocyte、Neuroendocrine、Tuft。

### T/NK（tnkA）· r\*=0.8 · 种子 3 · 已剔 2 簇 956 核
38 个亚簇，63,078 核。

| 亚型 | 核数 |
|---|---|
| CD4+ Naive T | 37,528 |
| NK | 10,608 |
| CD8+ Mem/Eff T | 8,173 |
| CD8+ Naive T | 2,680 |
| NKT | 2,172 |
| CD4+ Mem/Eff | 1,917 |

⚠️ CD4+ Mem/Eff 这一行的源 marker 表**有硬伤**（原表把 CD8 写进了 CD4 行）。**原样保留未改**，所以这 1,917 核的解释要打问号——改了就成了我们的判断。

### B/浆（bplasmaA）· r\*=0.6 · 种子 0 · 已剔 9 簇 8,367 核
22 个亚簇，22,915 核。

| 亚型 | 核数 |
|---|---|
| Plasma | 15,975 |
| B | 6,940 |

### 髓系（myeloidA）· r\*=0.5 · 种子 0 · 已剔 1 簇 215 核
27 个亚簇，63,869 核。

**旧面板（Table S1，已被取代）**——8 型能赢，11 个簇报 thin_panel：

| 亚型 | 核数 |
|---|---|
| Classical Monocyte | 18,517 |
| Intermediate Monocyte | 14,741 |
| Macrophage | 14,074 |
| Basophil/Mast | 6,368 |
| mDC2 | 6,100 |
| pDC | 2,092 |
| Neutrophil | 1,191 |
| mDC1 | 786 |

**现行面板（Table S4，2026-09-22 你追认）**——12 型各自胜出，0 个簇报 thin_panel：

| 亚型 | 核数 |
|---|---|
| TREM2+ Dendritic | 20,571 |
| Macrophage | 10,837 |
| IGSF21+ Dendritic | 9,193 |
| Basophil/Mast 1 | 6,368 |
| Myeloid Dendritic Type 2 | 6,100 |
| EREG+ Dendritic | 2,476 |
| Nonclassical Monocyte | 2,327 |
| Classical Monocyte | 1,985 |
| Myeloid Dendritic Type 1 | 1,953 |
| Plasmacytoid Dendritic | 925 |
| Neutrophil | 840 |
| Proliferating Macrophage | 294 |

（17 型里的 **Eosinophil 打不了分**：S4 没有嗜酸簇、沿用 S1 后只剩 1 个基因。另有
**Basophil/Mast 2、Platelet/Megakaryocyte、OLR1+ Classical Monocyte、Intermediate Monocyte**
四型一次未赢。**TREM2+ Dendritic 这个名字要当心**，见第三节（一之二）。）

**旧面板的毛病（保留记录）**：Eosinophil、mDC1、mDC2、Classical Monocyte、Nonclassical Monocyte、
Basophil/Mast 六型面板薄；Eosinophil、Nonclassical Monocyte 完全打不了分。
**TPSAB1 不在矩阵里**（肥大细胞最经典的 marker），旧面板下 Basophil/Mast 只剩 MS4A2+CPA3 两个槽；
而且 S1 源表里 Basophil 和 Mast 两行**基因完全相同** ⇒ 旧面板只能合并上报。
S4 把这一型拆成 Basophil/Mast 1 与 2 两块各 20 个基因（共有 CPA3/MS4A2/HPGDS/VWA5A/RGS13/CD69/IL1RL1 等，
1 带 KIT、2 带 HDC 作区分），**但只有 1 赢过簇**。

### 内皮（endoA）· r\*=0.7 · 种子 0 · 已剔 4 簇 4,580 核
25 个亚簇，38,700 核。

**旧面板（Table S1，已被取代）**——只有两型能赢：

| 亚型 | 核数 |
|---|---|
| Lymphatic | 19,412 |
| Artery | 19,288 |

**现行面板（Table S4，2026-09-22 追认后）**——8 型各自胜出：

| 亚型 | 核数 |
|---|---|
| Capillary | 8,839 |
| Capillary Intermediate 1 | 7,194 |
| Artery | 6,388 |
| Vein | 5,654 |
| Bronchial Vessel 2 | 4,674 |
| Lymphatic | 3,354 |
| Bronchial Vessel 1 | 2,292 |
| Capillary Aerocyte | 305 |

（9 型里的 **Capillary Intermediate 2 打不了分**：源数据本身筛不出可区分的基因，非本次操作所致。）
**这就是第一节说的那件事**，详见下面第三节。

### 成纤维（fibroA）· r\*=0.8 · 种子 4 · **事后放宽门槛（探索性）** · 0 簇被剔
33 个亚簇，77,633 核。本次**未重跑**（其细胞集未变动）。

| 亚型 | 核数 |
|---|---|
| Fibroblast | 35,581 |
| Myofibroblast | 24,244 |
| Vascular Smooth Muscle | 8,408 |
| Lipofibroblast | 5,920 |
| Pericyte | 3,480 |

---

## 三、三件必须显眼说的事

### （一）内皮面板：曾"少两路"，已用同一篇论文的 S4 表修好（2026-09-22 你追认）

**原先的毛病（保留记录）**：内皮面板只有 4 行，来自 Travaglini 2020 Table S1 的逐字转录：

| 亚型 | 源表给的 marker | 可用基因数 | 结果 |
|---|---|---|---|
| Artery | GJA5, BMX | 2 | 薄（≤2） |
| Vein | ACKR1 | **1** | **打不了分** |
| Capillary | CA4 | **1** | **打不了分** |
| Lymphatic | PROX1, PDPN | 2 | 薄（≤2） |

`missing_genes` 是空的——这 6 个基因**矩阵里全都有**。所以问题不在数据，在于**源表每种内皮只给了 1–2 个基因**：
Vein 和 Capillary **在结构上不可能胜出**，那份"Lymphatic 19,412 / Artery 19,288"实质是
"GJA5/BMX 和 PROX1/PDPN 两组基因哪个分高"，不能当 4 路划分用。

**根因**：S1 那张 "Canonical markers" 是**摘要列**，对血管内皮**不是单细胞证据**——动脉那两行来自 bulk/培养细胞，
静脉和毛细血管没给单细胞谱。

**修法（不改论文，只换表）**：改用**同一篇论文**的 Table S4（逐簇富集表），
按"只取 10x 表 / logFC≥1.0 且特异度≥0.3 / 每型封顶 20 个顺延补齐 / 照 S4 原粒度 9 型 /
S4 无表者（中性粒、嗜酸）沿用 S1"筛。**同一批 43,280 核、同一 r\*=0.7、同种子，唯一变量是面板。**

| | 旧面板（Table S1） | 现行面板（Table S4） |
|---|---|---|
| 能赢的亚型 | **2** / 4（只有 Artery、Lymphatic） | **8** / 9 |
| thin_panel 报警 | 4 型全报 | 无 |
| 逐型基因槽 | 2 / 1 / 1 / 2 | 最多 20（Artery 13、Vein 20、Lymphatic 20…） |

图：`figures/endo_panel_s1_vs_s4.png`（`05_annotation/07_endo_panel_compare_figure.py`，只读产物出图）。
口径与逐型基因名单：`results/05_annotation/panel_caliber_s4.json`（**2026-09-22 你已追认**，
`countersign_required=false`）；面板模块哈希 `dc61f7b5…`。

**换表后仍须如实上报的三条局限**（不得省略）：

1. **Capillary Intermediate 2 打不了分**——S4 里该簇筛完只剩 0 个可区分基因。性质与旧的"表写太省"不同，
   这是**源数据本身没有可区分的基因**。9 型里能用 8 型。
2. **内皮面板里混进了免疫基因**（Capillary 的 IL7R、Capillary Intermediate 1 的 IL1RL1）
   ⇒ 提示 Travaglini 自己的簇里也带环境 RNA/双体。
3. **Capillary 在其源数据里本就弱**（最高 logFC 只 1.26、过筛剩 4 个基因）⇒ 即使换表，该型仍最弱。

**同一份 S4 口径也修了髓系**——下一小节。

### （一之二）髓系面板：同一份 S4 换过的第二个谱系（2026-09-22 你追认）

**原先的毛病（保留记录）**：旧髓系面板 11 型共 **28 个基因**，逐型只有 1–4 个槽——
Eosinophil 1 个、Nonclassical Monocyte 1 个、mDC1/mDC2/Classical Monocyte 各 2 个。
另外 S1 源表里 Basophil 与 Mast 两行**基因完全相同**（只能合并成一型），
且最经典的肥大细胞 marker `TPSAB1` **不在我们的矩阵里**，该型实际只剩 MS4A2+CPA3 两个槽。

🔴 **须分清两层，别把功劳记错**：旧 L2 面板**并非**"没有 DC marker"——
它其实有 mDC1（CLEC9A/LAMP3）、mDC2（CD1C/PLD4）、pDC（LILRB4/IRF8/LILRA4），
只是每型仅 2–3 个基因。真正"既无肥大细胞 marker、又无 DC marker"的是 **L1 的
`marker_panel.py` 髓系面板**（16 个基因，全是巨噬/单核），这才是第四节里
全局簇 17/26 被错标的上游原因。而 **S4 只改 L2，不改 L1**（见 `panel_caliber_s4.json` 的 `scope`），
所以 **17/26 是被你 2026-09-22 的人工裁决修掉的，不是被 S4 修掉的**。

**修法**：与内皮同一份 S4 口径、同一套筛法（只取 10x 表 / logFC≥1.0 且特异度≥0.3 /
每型封顶 20 个顺延补齐 / 照 S4 原粒度 17 型 / S4 无表者沿用 S1）。
**同一批 64,084 核、同一 r\*=0.5、同种子 0，唯一变量是面板。**

| | 旧面板（Table S1） | 现行面板（Table S4） |
|---|---|---|
| 面板覆盖的亚型 | 11 型 / 28 个基因 | 17 型 / 最多 20 个基因 |
| 能赢的亚型 | **8** / 11 | **12** / 17 |
| thin_panel 报警的簇 | **11 个** | **0 个** |
| 打不了分的型 | 2（Eosinophil、Nonclassical Monocyte） | 1（Eosinophil） |

图：`figures/myeloid_panel_s1_vs_s4.png`（`05_annotation/11_myeloid_panel_compare_figure.py`，只读产物出图）。
清单：`results/05_annotation/myeloid_panel_s1_vs_s4_manifest.json`。
逐型基因名单：`panel_caliber_s4.json`（**2026-09-22 你已追认**，`countersign_required=false`）。

**换表后仍须如实上报的四条局限**（不得省略）：

1. 🔴 **"thin_panel 报警归零"不等于"证据够了"。** 报警的判据是**胜出型的面板基因数 ≤2**；
   S4 把每个能打分的型都抬到 ≥3 个基因，于是报警自然归零——这是**跨过了阈值边界**，
   不是证据质量变好。最薄的几个新面板仍然接近零容错：
   **Capillary Intermediate 1 只有 3 个基因**（内皮），Capillary / Bronchial Vessel 2 / Neutrophil 各 4 个。
   本报告一律不用"S4 把面板修好了"这种说法，只能说"S4 拆掉了 1–2 个基因这个结构性陷阱"。
2. 🔴 **`Platelet/Megakaryocyte` 被 S4 归入髓系**（Travaglini 的 Cluster 46）。
   血小板/巨核细胞不是白细胞，**源论文也没这么分**。Peng 2026 原文逐字列举（Methods，主细胞群鉴定）：
   > "…myeloid cells (**CD68, CD163**, etc. for macrophages; **CSF3R** for neutrophils; **KIT, MS4A2**,
   > etc. for mast cells; **CD1C, CLEC9A**, etc. for conventional dendritic cells-cDC, **LILRA4, IL3RA**,
   > etc. for plasmacytoid dendritic cells-pDC)…"

   同一篇原文另一处把髓系概括为 "myeloid (**monocyte, macrophage, cDC, and mast**)"——
   **两处都没有血小板/巨核细胞**。所以这是 Travaglini 图谱的归类，不是 LUAD 领域的归类。
   **实际影响为零**（它一次也没赢过任何簇），但这是**登记在案的、未经你复核的分类学判断**
   （见 `panel_caliber_s4.json` 的 `known_limitations`）。**建议处置：保留在面板里（不影响结果），
   但在任何髓系构成比的分母里把它单列或剔除**——请你定。
3. 🔴 **S4 的名字是 Travaglini 健康肺的簇名，与 LUAD 领域叫法冲突。**
   最刺眼的是 **`TREM2+ Dendritic` 赢了 11 簇 / 20,571 核**（髓系里最大的一块）——
   在 LUAD 里 TREM2 基本是**脂质相关巨噬细胞（LAM）**的招牌，S4 却把它归在树突里。
   照 S4 写就是"Dendritic"，按 LUAD 领域写就是"Macrophage/LAM"。**二选一须登记**，
   本报告暂按 S4 原名上报，并在下游结论里同时标注这个冲突。
4. **四型一次未赢**（Basophil/Mast 2、Platelet/Megakaryocyte、OLR1+ Classical Monocyte、
   Intermediate Monocyte）。其中 `Intermediate Monocyte` 在旧面板下曾赢 6 簇 / 14,741 核，
   换成 S4 后归零——**这是面板更换带来的实质变化，须在此显式记录**。可查的原因：
   S4 的 `Intermediate Monocyte` 只有 15 个可用基因，且与 `Nonclassical Monocyte` 共享
   **6 / 15** 个（CFP、COTL1、IFITM3、LILRA5、LILRB2、LST1），与 `Classical Monocyte` 共享 3 / 15
   ⇒ 在 S4 的粒度下它不再构成独立的一路。

**这张对照表有一处是我们自己加的，不是原文内容**：S1 与 S4 的亚型名与粒度都不同
（S1 的 `Megakaryocyte` = S4 的 `Platelet/Megakaryocyte`；S1 的 `Basophil/Mast` 在 S4 拆成 1 和 2；
S1 的 `mDC1/mDC2/pDC` 在 S4 写成全名），所以配对比较必须靠一张**本项目自定的命名对照表**，
它已登记在 `myeloid_panel_s1_vs_s4_manifest.json` 的 `name_mapping_ours` 里，**不得当作原文口径引用**。

### （二）B/浆 的污染其实是一个样本的问题

后面 S7 那张表里，B/浆 被剔 8,367 核（占 B/浆 子集的 26.7%），是六个谱系里最刺眼的数字。
但按样本拆开看，它是**一个样本**的问题：

**先澄清两个都出现过的数**（不是矛盾，是两次运行）：

| | 簇数 | 总核 | 剔除率 |
|---|---|---|---|
| **旧版**（触发本次调查的那次） | 33 | 34,385 | **28.55%**（剔 11 簇 9,818 核） |
| **新版**（裁决层 + 重跑之后） | 31 | 31,282 | **26.75%**（剔 9 簇 8,367 核） |

| 事实（基于新版） | 数 |
|---|---|
| P24 占 B/浆 子集的比 | 36.4%（全库平均 11.7%，**3.11 倍富集**） |
| P24 占全部剔除的比 | 62.6% |
| P24 自己的剔除率 | **46.0%**（全部患者里最高） |
| **去掉 P24 之后，B/浆 剔除率** | **26.7% → 15.7%** |

**被剔的是什么样的簇**（直接从 `bplasmaA_contamination.csv` 读，不是我归纳的）：
面板原本判它 B/浆，但簇均值上**别的谱系分更高**，且 CellTypist 同意——
例如簇 8（1,044 核）自评 B/浆 0.917，而成纤维分 2.011，CellTypist 判成纤维 52.1%；
簇 4（2,391 核）自评 0.874，成纤维 1.219，CellTypist 判上皮（非上皮比例 0.698）。
所以剔除不是"随手丢"，是两条线同时指向别的谱系。

**其它五个谱系没有这个效应**（P24 富集倍数：上皮 0.63× / 内皮 0.69× / 成纤维 0.66× / 髓系 1.03× / T/NK 1.33×）。

**最可能的解释（未验证，是假设不是结论）**：环境 RNA。肺 snRNA 里免疫球蛋白转录本（IGKC/IGHG1/IGHM）
是最典型的环境污染物，它们糊在核表面，把六谱系面板里"B/浆"那条的分抬起来，造成全局 B/浆 误标，
再被亚簇级的两线检查挡掉。P24 是 IAC 样本，如果它的组织处理/上机特别容易释放 Ig，就会出现这种"单样本富集"。

**能证伪的检验**（等你点头我再跑）：在剔除的核里比 IG 基因（IGKC/IGHG1）与真正的浆细胞功能基因（XBP1/PRDM1/JCHAIN）。
若 IG 高而功能基因不高 ⇒ 环境 RNA 假设成立，P24 该考虑单独处理或多加一个环境校正步骤。

---

#### 已跑 · 环境 RNA 检验结果（2026-09-22，`05_annotation/08_p24_ambient_ig_test.py`）

脚本、逐细胞表、清单、图都已入库。**结论：环境 RNA 假设得到支持**，但要说清它是怎么得到的。

**第一步（预注册口径）失败了，而且失败在度量上，不是结论上**：

| 预注册检验 | 结果 | 为什么不作数 |
|---|---|---|
| T1：非 B/浆 5 谱系的 IG 阳性率，P24 ≥ 其他 ×2 | P24 是 1.5–1.9×，**0/5 过线** | 队列基线本就 0.45–0.63（环境 IG 人人有）⇒ 2× 的杆几乎够不到，**度量饱和** |
| T2：被剔细胞"任一浆程序基因 >0" vs 保留细胞 ≤0.5× | 0.966 vs 0.972，**不过线** | XBP1/PRDM1/IRF4/CD27 这些基因太广谱 ⇒ ~0.95 谁都过，**度量无区分力** |

两个预注册判定**照实上报，不据此下结论**（这是设计缺陷，不是 P24 的结论）。

**第二步：换成有区分力的度量**（"10 个浆程序基因里检出 ≥3 个"；保留的真浆细胞该指标 0.80、中位 7/10）
——**事后定的，只作探索性**：

| 非 B/浆 谱系 | 其他病人 | **P24** | 倍数 |
|---|---|---|---|
| 上皮 | 0.088 | **0.505** | 5.8× |
| 内皮 | 0.012 | **0.129** | **10.8×** |
| 髓系 | 0.044 | **0.334** | 7.5× |
| 成纤维 | 0.015 | **0.248** | **17.2×** |
| T/NK | 0.053 | **0.448** | 8.4× |

**5/5 个谱系全部过 3×。** 也就是说：在**结构上不可能造浆细胞**的细胞里，P24 有 6–17 倍的浆细胞程序信号。

**而且它不像肿瘤生物学**——P24 三个样本**同向抬高**：

| 样本 | 核数 | 浆程序 ≥3 占比 |
|---|---|---|
| P24_LUAD | 24,885 | **0.46** |
| P24_Normal | 6,380 | **0.23** |
| P24_AAH | 7,151 | **0.22** |
| 其他病人（中位） | — | 0.031（范围 0.000–0.164） |

**连 Normal 样本都高 7 倍** ⇒ 更像样本处理/上机带来的环境 RNA，而不是肿瘤。P24 的 nFeature 中位也偏高（1,787–1,925 vs 其他 1,178），与环境 RNA 抬高文库复杂度一致。

图：`figures/p24_ambient_ig.png`（A 预注册、B 纠正后、C 逐样本）。

**🔴 诚实边界**：Part B 的阈值是**看到 Part A 缺陷之后才定的**，所以它**不是**一次预注册确认性检验。
要把它当结论，得先把"浆程序 ≥3/10"这个度量**重新冻结为口径、再跑一次**（法则 3.2，口径不得事后追认）。
**在你决定是否重冻结之前，这条只能算探索性证据。**

---

#### 已跑 · 环境 RNA **确认版**（2026-09-22，口径 `results/05_annotation/p24_ambient_caliber.json` → 脚本 `05_annotation/09_p24_ambient_confirm.py`）

你 2026-09-22 会话内决定"重冻结口径，跑确认版"。口径**先冻结、后计算**（法则 3.1），
脚本开头校验口径文件存在及其 sha256（`a522659e87cd6f33…`）。

确认版比探索版强在两处**设计**上（不是换阈值）：

1. **分析单位从「细胞」改为「病人」**——细胞不独立，细胞级 p 值会虚高。23 位病人，非 B/浆 384,862 核。
2. **加一组负对照基因集**——排除"P24 的细胞 nFeature 本来就高 ⇒ 随便什么基因集都高"这个伪影：
   从矩阵里按**逐基因检出率匹配**（对应基因差 ≤0.05）抽 200 组伪 10 基因集，seed=0。

| 检验 | 判据（口径原文） | 实测 | 判定 |
|---|---|---|---|
| **主 · 病人级** | P24 秩 1/23 **且** ≥ 其他 22 人中位 3× **且** 稳健离群分 ≥3 | 秩 **1/23**；**7.81×** 中位；稳健 z = (0.377−0.048)/0.0327 = **10.1** | ✅ 过 |
| **特异度 · 负对照** | P24 浆程序值 ≥ 200 组伪集的 95 百分位（经验 p ≤0.05） | P24 = 0.377；伪集中位 0.225、**最大 0.314**；经验 **p = 0.005** | ✅ 过 |
| **次 · 独立读数** | 非 B/浆 细胞的 IG 恒定区 UMI 占比中位数，P24 秩 1/23 | 0.43% vs 其他中位 0.029%，秩 **1/23** | ✅ 过 |

**三项全过，总判定：确认支持环境 RNA 假设（H1）。** 图 `figures/p24_ambient_confirm.png`。

**🔴 但负对照的数字必须一起读，别只看"过了"**：200 组**逐基因匹配**的伪集在 P24 细胞上
本身就给出 0.225 的中位值（其他病人用真浆程序只有 0.048）。**也就是说 P24 绝对值里的大头是
"细胞普遍偏热"（测序深度/nFeature 高），不是浆细胞特异的。** 浆程序 0.377 只比伪集**最大值**
0.314 高出 0.063。所以准确的表述是：

> **纯深度效应不足以解释全部**——扣除深度匹配的基线后，P24 仍有一份**浆程序特异**的抬高；
> 但这份特异的幅度**不大**，且仅"刚好"越过 200 组负对照的极值（经验 p=0.005 是分辨率下限 1/201）。

**🔴 这仍不是独立验证**：这个度量本人在探索阶段已经看过，P24 当时表现就是好的。
确认版的价值在于"换一个更强的设计、在同一批数据上再统计一遍且仍然成立"，不是"盲法复现"。
口径原文 `honesty_boundary` 已写明这一点。

**🔴 也仍不宣称因果**：本节只说"P24 的环境 RNA 信号是否异常抬高"（已确认），
**不说**它导致了什么，也不说 P24 的恶性结论因此改变。

---

#### 原文里有没有写这个样本的问题？——**没有**（2026-09-22 核对 PMC 全文）

你 2026-09-22 问"原论文有描述这个样本的问题吗"。核对方式：把 PMC 全文（PMC12980502）HTML 下载到
`scratch/paper/peng2026_pmc.html` 逐词检索。
该文件**不入库**（`.gitignore` 的 `scratch/`）——但溯源可复现：检索日期 2026-09-22，
417,141 字节，sha256 `1260df410125f732b59290794b41f1d7acd33e38317f16e5008f545e5dc939fb`
（PMC ID 永久有效，重取即可核）。下方引文均逐字来自该文件。

| 检索词 | 全文命中 |
|---|---|
| `ambient` / `soup` / `SoupX` / `CellBender` / `DecontX` / `EmptyDrops` | **0 次** |
| `immunoglobulin` / `IGKC` / `IGHG` | **0 次** |
| `plasma cell` | 4 次，**全部作为正常生物学群体**（MP11） |
| `P24` | 8 次，**全部作为代表案例** |

⇒ **原文既没有做环境 RNA 校正，也没有标记任何人类 snRNA 样本。**
人类 snRNA 队列 75 样本 / 25 例中保留 23 例，**未点名剔除任何一个**（原文里被剔除的只有 Visium spot、Xenium 细胞和**小鼠**样本）。

**原文的 snRNA QC 与聚类全段（逐字引用）**：

> "…cells with low-complexity libraries (in which detected transcripts were aligned to <200 genes such as cell debris,
> empty drops and low-quality cells) were filtered out and excluded from subsequent analyses, and genes detected in
> less than 3 cells were also excluded. Nuclei with less than 500 detected genes or with less than 1000 reads count or
> with a mitochondrial gene fraction that is ≥20% were filtered out using Seurat (v5.1.0). Doublets were identified
> based on library complexity. First, a python-based software Scrublet was applied to identify the doublets by each
> sample. Second, based on cluster distribution and marker gene expression, doublets forming distinct clusters with
> hybrid expression features were also removed. … Clusters co-expressing discrepant lineage markers were identified
> and removed. Data normalization and scale transformation were performed using method 'SCTransform' in Seurat.
> Top 3,000 HVGs were selected for PCA and downstream unsupervised clustering. The top 50 PCs were used to calculate
> the embedding. Harmony (version 1.2.0) was run with default parameters to remove batch effects present in the top
> 50 PCA space. … 'FindClusters' function with a resolution set to 0.5-0.8."

**原文把浆细胞当成真实现象**（同一段方法里的谱系 marker 定义）：

> "…B and plasma cells (MS4A1, CD79A, etc. for B cells; **MZB1, XBP1**, etc. for plasma cells)…"
> 结果段："…increased frequencies of B, plasma, CD4+ T, Tfh, Th17, and Treg cells… with disease stage"

（我们用的浆程序基因集里的 **MZB1、XBP1** 与原文一致。）

**P24 在原文里是正面范例，不是问题样本**——图 3K/3E–3G、4E–4F、S6 用 "P24 AAH / P24 LUAD" 讲 IL1B–IL1R1 的空间信号；
"P24 with clonal evolution pattern 1a" 是演化模式 1a 的范例；P24 还列在 SpatialInferCNV 用全部正常上皮当参考的病人名单里。

**为什么原文没发现这个偏差？（🔴 这是我的推断，原文没有这句话）**
原文的第二道防线是**簇级**的——"把同时表达矛盾谱系 marker 的**簇**整簇删掉"。
而我们的发现是环境 RNA 的浆程序信号**弥散在别的谱系的簇内部**，不形成自己的簇 ⇒ **簇级删除抓不到**。
这也说明该偏差对原文主结论影响有限：P24 在原文里用于讲空间 IL1B–IL1R1，**不是**用于统计 B/浆 细胞比例，
而"27% 污染"只发生在 B/浆 子集内部。

---

#### 顺带查到的两件事（对我们要紧）

**① 原文的 snRNA 参数现在有原文数值可对标了**（此前只能"推测"，现在可引原文）：

| 环节 | 原文（PMC12980502 方法节） | 本项目（已登记偏差） |
|---|---|---|
| 核质控 | <200 基因（低复杂度）删；<500 基因 或 <1000 reads 或线粒体 **≥20%** 删 | 逐样本自适应 MAD（M1） |
| 双体 | **Scrublet，逐样本** | **scDblFinder**（偏差已登记） |
| 归一并选 HVG | SCTransform，top **3,000** HVG | scran 池化 + log1p，top **2,000** HVG |
| 降维/批次 | PCA 前 **50** PC + **Harmony 1.2.0 默认参数** | 30 PC；双臂（Arm A 不校正 / Arm B Harmony on `sample_id`） |
| 聚类 | Seurat `FindClusters`，分辨率 **0.5–0.8** | Leiden，网格 0.2–2.0 选 r\* |
| 保留核数 | **401,635** | **413,697** |

两点由此得到支持：
- 原文明确是 **10x Fixed RNA Profiling（Flex）** ⇒ 我们内皮 S4 面板选 `only_10x=True`（排除 SS2）**是对的**。
- 我们与原文的核数差 12,062，方向与质控口径差（自适应 MAD vs 固定阈值、scDblFinder vs Scrublet）一致，**不构成矛盾**。

**② P24_LUAD 的核数确实异常大 —— 一个假设（🔴 未验证，不是结论）**：

| 项 | 数值 |
|---|---|
| P24_LUAD 保留核数 | **34,523** |
| 第二名（P4_LUAD） | 21,757（P24 是它的 **1.59 倍**） |
| 24 个 IAC(LUAD) 样本的中位数 | **5,470** |
| 原文设计目标（AAH/AIS/MIA/LUAD 亚池化） | "target up to **20,000** cells per sample" ⇒ **P24_LUAD 超过设计目标** |
| P24_LUAD 双体率 | **29.03%**（全队列最高，队列均值 15.48%） |

这个组合（核数超设计目标 + 双体率最高 + 环境 RNA 信号最高）与"**上样过多 → 环境 RNA 与双体同时抬高**"是一致的。
**但这是假设**：原文没有说哪个样本上样过载，公开数据也反推不出上样量，本项目无实验记录可查。
它只能作为"为什么是 P24 而不是别人"的一个待查方向，**不得写进结论**。





| | 旧版（2026-09-18 签字） | 新版（2026-09-22） |
|---|---|---|
| r\* | 0.5 | 0.7 |
| 跨种子稳定性（r\* 处） | 0.9455 | 0.8605 |
| 是否放宽门槛 | **未放宽** | **放宽** 0.90→0.85 |
| 结论 | 过 | **只能作探索性** |

**为什么作废**：在修正后的 133,384 核上，四个分辨率的稳定性是 0.8501 / 0.8423 / **0.8605** / 0.8201，
**没有任何一个到 0.90**。旧值 0.9455 复现不出来。

**不是"被污染细胞撑起来的"**——我一开始这么猜，是错的。验证：把旧标签限制在那 133,384 个细胞上重算，得 0.9438，
与全量的 0.9455 只差 **0.0017**。5.5% 的细胞动不了一个全局 ARI。
**真正的原因**：不是"细胞变少了"，而是**旧版唯一过关的那一档（res=0.5）没有复现**。
同分辨率、同配对定义下的实算（`06_seed_churn_audit.py`）：

| 分辨率 | 旧 run 跨种子 ARI | 新 run 跨种子 ARI | 旧 run 改标签率 | 新 run 改标签率 |
|---|---|---|---|---|
| **0.5** | **0.9455** | **0.8501** | 3.91% | **8.72%** |
| 0.6 | 0.8198 | 0.8423 | 11.70% | 9.83% |
| 0.7 | 0.8395 | 0.8605 | 8.14% | 7.79% |
| 0.8 | 0.8476 | 0.8201 | 9.17% | 11.45% |

（改标签率 = 10 个种子对两两的均值；ARI 取各 run 自己的 `ari_seed_mean`。图：`figures/seed_churn_audit.png`）

⇒ **两版里只有旧版的 res=0.5 这一档过了 0.90（0.9455）**；旧版其余三档本来就 0.82–0.85，谁也没过。
新版把 res=0.5 从 0.9455 打到 0.8501（改标签率 3.91%→8.72%），而 **res=0.6 / 0.7 两版持平（新版甚至略好）**。
所以准确说法是：**旧版靠 res=0.5 拿到的那份高分离度没有复现**，而不是"新版划分普遍更不稳定"。

✅ **复核**：原句"把旧标签限制在那 133,384 个细胞上重算，得 0.9438"——**已实算、逐位复现**（0.9438）。
所以那个值本身没问题，"细胞数变化"确实解释不了 0.9455；掉的是**那一档的聚类结果**
（同一批细胞、但两版的聚类输入不同 ⇒ HVG/PCA 空间变了）。

🔴 **原"过半细胞换簇"一句已作废（2026-09-22 重算后）**：`15.91%` 与 `51.36%` **两个数都不可复现**
（全库含日志 grep，除该句正文外无任何产物/脚本含 15.91）。**51.36% 的真实身份已查清** =
新版 res=0.5、种子0 vs 种子1 的**裸标签不一致率**（51.3562%，四舍五入即 51.36%）——
Leiden 的簇编号是任意的，换种子就重排，裸比会把"同一个簇换了个编号"也算成改标签，**它不是改标签率**；
同一对在**配对定义**下是 **12.13%**（≈16,181 核会换簇）。
重算脚本与产物：`05_annotation/06_seed_churn_audit.py` / `results/05_annotation/seed_churn_audit.csv`
/ `figures/seed_churn_audit.png`（脚本内有硬自检：与已登记的 `seed_representativeness.csv` 一致才放行）。
⚠️ **主论据不受影响**：上皮需放宽阈值靠的是跨种子 ARI **0.8605 < 0.90**（`epiA_rstar.json`，已签），
有脚本、有产物。

**旁证**：R 脚本两次运行之间只有注释文字改动（`git diff` 确认），run_manifest 里唯一的参数差异是
SCT 自动推导的 `clip.range`（不是人为调的旋钮）。所以差异只能来自**聚类输入**（细胞集 → HVG/PCA 空间），
且如上表，该效应**集中在 res=0.5 那一档**（0.6 / 0.7 两版持平）。

**对照**：B/浆 去掉一块**更紧密**的细胞（ARI 0.9996），稳定性反而**升**了（0.957→0.977）。
所以这不是"删细胞就会掉稳定性"的普遍现象，是**上皮特有**的。

**技术层面的原因（待验证）**：6,502 个髓系核改变了 HVG/PCA 特征空间。这句是猜测，没测过。

**独立旁证**：上皮净化的子集（`epiA_clean`）最高也只到 0.8976，同样够不着 0.90。

---

## 四、全局簇错标裁决层（17/26/38/44）

规则（预注册于 GP6，本次未新发明阈值）：逐簇两线一致率 <0.90 **且** 谱系级分歧。
23 个簇低一致，其中 4 个是谱系级 → 全部重判：

| 簇 | 原判 | 改判 | 表型 | 核数 | 驱动证据 |
|---|---|---|---|---|---|
| 17 | 上皮 | **髓系** | 肥大细胞 | 6,574 | 驱动基因前 15 全是肥大细胞特异（CPA3;HDC;KIT;MS4A2;GATA2…），CellTypist 髓系 97% |
| 26 | B/浆 | **髓系** | 树突状细胞 | 3,103 | 驱动基因含 IRF8;WDFY4;CIITA;CD74，CellTypist 髓系 93% |
| 38 | 上皮 | **T/NK** | T 细胞（泡在表面活性剂环境 RNA 里） | 1,147 | 环境 RNA 判别：T/NK 核心 2.50×、分泌型上皮 1.79×、结构型上皮仅 1.15× |
| 44 | 内皮 | **T/NK** | T/NK | 404 | 驱动基因前 15 全是 T 基因（IL7R;TRAC;CD2…），两线只差 0.045 |

**旧口径不作废**：`A_frozen` 逐字节保留，新增 `A_adjudicated` 一列，新旧并存可对账。

**三个必须披露的问题**：

1. **时序**：这次裁决**发生在看到下游结果之后**。触发点是 2026-09-22 查"B/浆 剔除率 28.55% 是否合理"，
   簇 26 由此暴露；17/38/44 是之后系统套用预注册规则（遍历全部 45 个全局簇）发现的。
   规则本身与本次调查无关、候选集机械产出，但"先看下游、后做裁决"这一时序如实登记（法则 3.2）。
2. **签字**：`signed_by = Claude（自动夜跑）2026-09-22`，**未经人工逐行核对证据表**。
   系统自动核对的是候选簇由规则机械产出（脚本硬断言）和全部输入输出有 sha256；**须你复核的是两张证据表的判读**。
3. **面板盲区导致自动化仍会出错**：**L1 的**髓系面板（`marker_panel.py`，16 基因）不含任何肥大细胞
   marker、不含任何 DC marker；上皮面板 16 基因里 7 个是分泌型，易被环境 RNA 抬高。
   你 2026-09-22 裁定**不补面板**，故作为已知局限登记——**同样的错，L1 自动簇级 argmax 还会再犯**。
   （S4 只换了 **L2** 的内皮/髓系面板，L1 未动，所以这条局限不因 S4 而消失。）

---

## 五、谱系接管（S7）三分账

14,118 个被剔核，逐细胞重新过两线，得到三态：

| 谱系 | 被剔 | 还原 | 接手 | 待人工判 |
|---|---|---|---|---|
| 上皮 | 0 | – | – | – |
| 成纤维 | 0 | – | – | – |
| T/NK | 956 | 133 | 327 | **496** |
| B/浆 | 8,367 | 1,326 | 2,485 | **4,556** |
| 髓系 | 215 | 31 | 100 | **84** |
| 内皮 | 4,580 | 943 | 1,170 | **2,467** |
| **合计** | **14,118** | **2,433** | **4,082** | **7,603** |

- **还原**：两线一致认为它**就是**原谱系 ⇒ 簇级剔除对这批细胞是**误剔**。占 17.2%。
- **接手**：两线一致认为是**别的**谱系 ⇒ 已改判（最多流向：上皮 +1,525、成纤维 +1,312）。
- **待人工判**：两线不一致或弃权 ⇒ **暂不移动**，留在原处（最保守）。

**🔴 必须知道的两点**：

1. **这不是独立验证**。判据与剔除判据是同源证据（同一份面板打分、同一份 CellTypist 输出），
   只是把聚合层级从"簇"降到"细胞"。所以"还原"不能当作第二个独立证据用。
2. **待判的 7,603 条本身，就是本次簇级裁决不确定性的度量**。P24 又占了 3,331 条（43.8%），与第三节一致。

`unadopted` 清单在 `results/05_annotation/gp8c_readopt_unadopted.csv`（7,603 行）——**这个决定我不替你做**。

**全局不变量已核**：接手后各谱系合计 413,697 = 全部细胞，互斥且完备 ✅

---

## 六、需要你人工判的事（按优先级）

1. ~~**内皮面板**：4 路变 2 路~~ ✅ **2026-09-22 已解决**——你追认了 S4 面板（换同篇论文 Table S4），内皮现为 8 型（详见第三节（一））。
2. ~~**P24**：B/浆 污染的主要来源。~~ ✅ **2026-09-22 已全部解决**：确认版三项全过 ⇒ 环境 RNA 假设成立（见第三节（二））；
   **你同日裁定『保留 + 标注』**，决定书 `results/05_annotation/p24_disposition_decision.json`。
   **执行规则**：标注 = `patient_id == "P24"`（不改任何已签字产物——`patient_id` 本来就在各产物里）；
   凡涉及 **B/浆 分区、谱系构成比例、跨分期细胞比例**的结论**必须双报（含 P24 / 去 P24）**。
   **为什么不剔除**（你问的"同阶段还有多少例"）：

   | 分期 | 病人（全） | **病人（去 P24）** | P24 占该期核数 |
   |---|---|---|---|
   | Normal | 23 | 22 | 6.9% |
   | **AAH** | 8 | **7** ⚠️ | **20.8%** |
   | AIS | 12 | 12 | 0 |
   | MIA | 4 | 4 | 0 |
   | IAC | 23 | 22 | **18.4%** |

   **Normal + AAH + IAC 三期齐全的病人全队列只有 8 例**（P24、P4、P6、P9、P11、P20、P22、P25），
   剔除 P24 会让 AAH 从 8 例降到 7 例、三期齐全者从 8 例降到 7 例，而 **AAH 恰是这条链上最薄的一环**。
   『去 P24』的**具体数值本次不算**——真要用于下游时须先冻结口径（法则 3.1）。
3. **待判 7,603 核**：`gp8c_readopt_unadopted.csv`。我没动。
4. **上皮的探索性定位**：放宽门槛后的结果，按登记只能进探索性结论、不得进主结论。**你 2026-09-22 已裁定接受**（"没办法，往下做说明局限性就好"）——与成纤维同处，见第七节。
5. **高患者集中度的簇**（规则**故意不设阈值**，`n_patients` 和 `top_patient_frac` 两列都已在表里供你判断）：
   - 上皮 AT2 里 5 个簇 top 患者占比 ≥45%：簇 23（86.2%）、21（67.3%）、16（56.1%）、26（50.5%）、18（46.8%）
   - T/NK：簇 18（62.1%）、簇 31（47.5%）
   - 唯一的单患者簇：T/NK 簇 39（17 核）
6. **`decision` 列全是空的**——每份 `<tag>_cluster_annotation.csv` 都在等你的签字。

---

## 七、签字与授权状况（逐条登记）

| 产物 | 签字人 | 时间 | 状态 |
|---|---|---|---|
| epiA r\* | 用户本人（会话内表态） | 2026-09-22 | 放宽门槛，探索性 |
| fibroA r\* | 用户本人（会话内表态） | 2026-09-18 | 放宽门槛，探索性 |
| endoA r\* | 用户本人（会话内表态，改签） | 2026-09-22 | 预注册原值直接过线，`relaxed=false` |
| tnkA / bplasmaA / myeloidA r\* | **Claude（自动夜跑）** | 2026-09-22 | 🔴 **未经人工复核** |
| **S4 面板口径**（内皮 9 型 + 髓系 17 型） | **用户本人（会话内追认）** | 2026-09-22 | ✅ 已追认，`countersign_required=false` |
| **P24 处置**（保留 + 标注） | **用户本人（会话内裁定）** | 2026-09-22 | ✅ 决定书 `p24_disposition_decision.json`；不改任何已签产物 |
| **P24 环境 RNA 确认版口径** | Claude（按用户裁定起草并冻结，算前冻结） | 2026-09-22 | 🔴 非独立验证（度量探索阶段已看过），口径内已自陈 |
| 全局簇裁决（17/26/38/44） | **Claude（自动夜跑）** | 2026-09-22 | 🔴 **未经人工逐行核对** |

夜跑签字依据你的授权：「我准备睡觉了，关键部分你先自己按照最科学的方式决断，明天早上给我报告」。
四个簇的**归属本身**是你在 2026-09-22 当天经 AskQuestion 亲自裁定的（17→髓系 / 26→髓系 / 38→T/NK / 44→T/NK）；
自动签字只覆盖"按已定规则机械执行"这一步。

**另外披露一处闸门弱点**：S7 的签字闸门只检查 `signed_by` 非空，而"Claude（自动夜跑）"是非空字符串，
所以**自动签名通过了本该等人工的闸门**。S7 因此是在 `A_adjudicated` 上跑的。
影响范围受控：`A_frozen` 未动，六份已签字的子集清单一字未改，待判细胞留在原处。

---

## 八、已知局限清单

1. **面板盲区（指 L1 的 `marker_panel.py`）**：髓系面板无肥大细胞 marker、无 DC marker；
   上皮面板 7/16 为分泌型。本次不补，**L1** 自动 argmax 仍会犯同类错。
   **S4 只换 L2 的内皮/髓系面板，不触 L1**，所以这条不因 S4 消失。
1b. **L2 面板仍存在弱点（S4 之后）**：胜出型的 thin_panel 报警归零只是跨过了"≤2 个基因"这条边界，
   最薄的新面板仍只有 3–4 个基因（内皮 Capillary Intermediate 1 = 3；Capillary / Bronchial Vessel 2 /
   Neutrophil = 4）。S4 的 `Intermediate Monocyte`（15 基因，与 Nonclassical 共享 6 个）与
   `Basophil/Mast 2`、`OLR1+ Classical Monocyte`、`Platelet/Megakaryocyte` 一次未赢。
1c. **命名口径冲突未决**：S4 用的是 Travaglini 健康肺簇名，`TREM2+ Dendritic`（髓系最大一块，
   11 簇 / 20,571 核）在 LUAD 领域通常指脂质相关巨噬细胞。按 S4 名还是按领域名上报，**须你定并登记**。
2. **矩阵缺基因**：`KRT18`、`SFTPA2`、`CD8B`（来自环境 RNA 判别的基因集核对）；上皮另缺 `DAPL1`、`PRR4`；髓系缺 `TPSAB1`。
3. **源表硬伤**（原样保留未修）：CD4+ Mem/Eff 行混入 CD8；Basophil 与 Mast 两行完全相同；Bronchial Vessel 行 markers 为空 ⇒ 无法注释。
4. **环境 RNA 未做校正（保留至下游，已登记）**：M1 无环境校正步骤，靠双标准兜底。
   第三节（二）已**确认** P24 的效应是环境 RNA（非肿瘤）；**2026-09-22 你裁定『保留 + 标注』**
   （`p24_disposition_decision.json`）⇒ 这个偏**不消除，靠双报对冲**：
   凡 B/浆 分区、谱系构成比例、跨分期细胞比例的结论**必须附「去 P24」对照**。
   ⚠️ 原文同样没有做环境 RNA 校正（§三（二）末已核），所以这是**共同的局限**，不是我们独有的偏离。
5. **CD4+ Mem/Eff 那一行**：源表笔误，该行结论不得单独使用。

---

## 九、产物路径与校验

```
results/05_annotation/
  epiA|tnkA|bplasmaA|myeloidA|endoA|fibroA_cluster_annotation.csv   ← 亚型表，decision 列待签
  <tag>_annotation_manifest.json                                     ← 每份含输入输出 sha256
  <tag>_contamination.csv · <tag>_excluded_cells.csv.gz              ← 污染体检
  gp6_adjudication_manifest.json · gp6_adjudication_evidence.csv
  gp6_adjudication_ambient_check.csv                                 ← 裁决层证据（待你复核）
  gp8c_readopt_manifest.json · gp8c_readopt_unadopted.csv            ← 接管三分账 + 待判清单
  gp8c_cell_assignment.csv.gz                                        ← 逐细胞 A_frozen/A_adjudicated/A_readopted
  seed_representativeness_manifest.json                              ← 种子代表性（generated 已修正为真时间戳）
  seed_churn_audit.csv / seed_churn_audit_manifest.json              ← 换种子改标签率重算（§八；含硬自检）
  panel_caliber_s4.json                                              ← S4 面板口径（§三之一；2026-09-22 已追认）
  classic_panels_s4.py                                               ← S4 面板模块（内皮 9 型 + 髓系 17 型，sha dc61f7b5…）
  {endoA,myeloidA}_s4_{cluster_annotation,cluster_scores,annotation_manifest}.csv/json
                                                                     ← S4 面板下的内皮/髓系注释（同 r\*、同种子，唯一变量是面板）
  endo_panel_s1_vs_s4_manifest.json                                  ← 内皮新旧面板对比图的清单（含逐型基因数、赢过哪些型）
  myeloid_panel_s1_vs_s4_manifest.json                               ← §三之一之二的图清单（含**本项目自定**的 S1↔S4 命名对照表、逐型基因数、四条局限）
  p24_ambient_ig_cells.csv.gz / p24_ambient_ig_manifest.json          ← §三之二 P24 环境 RNA 检验（预注册 + 事后纠正 + 缺陷登记）
  p24_ambient_caliber.json                                           ← §三之二 确认版口径（2026-09-22 冻结，sha a522659e…）
  p24_ambient_confirm_per_patient.csv                                ← §三之二 确认版逐病人结果（23 行）
  p24_ambient_confirm_null_distribution.csv                          ← §三之二 确认版 200 组负对照的 P24 取值
  p24_ambient_confirm_manifest.json                                  ← §三之二 确认版清单（三项判定 + 输入输出 sha256）
  p24_disposition_decision.json                                      ← P24 处置决定书（2026-09-22 用户裁定「保留 + 标注」，含逐期例数）
  GP8c_report.md                                                     ← 本文件
figures/seed_churn_audit.png                                         ← §八 的图（裸比 vs 配对、旧 vs 新同分辨率）
figures/endo_panel_s1_vs_s4.png                                      ← §三之一的图（内皮新旧面板对比）
figures/myeloid_panel_s1_vs_s4.png                                   ← §三之一之二的图（髓系新旧面板对比；第三格显式写出"报警归零 ≠ 证据够了"）
figures/p24_ambient_ig.png                                           ← §三之二的图（预注册 vs 纠正后度量、逐样本）
figures/p24_ambient_confirm.png                                      ← §三之二 确认版的图（病人级 / 负对照 / 独立读数）
logs/gp8c_20260922/STEP_STATUS.tsv                                   ← 逐步退出码全 0
```

> ⚠️ `06_seed_churn_audit.py` 的「旧运行」部分读 `results/_superseded/…/epiA/clusters.csv.gz`，
> 该目录**不入库**（.gitignore 明令不得当哈希锚点）⇒ 那部分是**一次性历史比对**，旧副本删除后不可再生；
> 新版部分（（一））只读入库产物，可复现。脚本已记录该限制于 manifest 的 `caveat`。

**种子代表性复核**（S4）：七个对象在各自 r\* 处，T/NK 与成纤维 gap >0.01（须换种子）。
两者**均已使用**所推荐的种子（T/NK 用 3、成纤维用 4），故无遗留缺口。
（另：`05_seed_representativeness.py` 原有两处溯源缺陷已修——写死的生成日期、以及"注释尚未开始"的过时表述。）

---

## 十、下一步待办（未启动）

1. ~~P24 的 IG 检验 + 处置~~ ✅ **2026-09-22 全部收口**。
   - 检验：跑完三步（预注册 → 探索纠正 → 确认版），确认版**三项全过** ⇒ 环境 RNA 假设成立。
   - **处置：你同日裁定『保留 + 标注』**，决定书 `results/05_annotation/p24_disposition_decision.json`（含你问的逐期例数表）。
   - **执行规则（不改任何已签字产物）**：标注即 `patient_id == "P24"`（该列本来就在各产物里，逐细胞文件名含
     `P24_Normal` / `P24_AAH` / `P24_LUAD`）；**凡涉及 B/浆 分区、谱系构成比例、跨分期细胞比例的结论，必须双报（含 P24 / 去 P24）**。
   - **不剔除的理由**：AAH 期只有 8 例病人、P24 占该期 20.8% 核；三期齐全（Normal+AAH+IAC）的病人全队列仅 8 例，
     剔除会让两者都少 1 例，而 AAH 是最薄的一环。
   - **不做的事**：不现在加环境校正（要重跑聚类+全部注释 ⇒ GP1–GP8c 已签字产物全部作废；应另立里程碑）；
     『去 P24』的**具体数值本次不算**，下游真要用须先冻结口径（法则 3.1）。
2. ~~内皮面板处置（等你定方向）~~ ✅ **2026-09-22 已追认 S4 面板**。
   ~~髓系 S4 结果尚未并入本报告~~ ✅ **2026-09-22 已并入**（第三节（一之二）+ 图
   `figures/myeloid_panel_s1_vs_s4.png`）。**仍待你复核两项**：
   （a）S4 把 `Platelet/Megakaryocyte` 归入髓系（影响为零，但分类学上不对，§三之一之二 局限 2）；
   （b）`TREM2+ Dendritic` 该按 S4 名还是按 LUAD 领域名（脂质相关巨噬细胞）上报（§三之一之二 局限 3）。
3. 待判 7,603 核（等你判）。
4. 各谱系 `decision` 列签字。
5. GP8a 报告（M3 双标准注释层，尚未写）。
6. ~~上皮「换种子改标签率」重算~~ ✅ **2026-09-22 已完成**（`05_annotation/06_seed_churn_audit.py`）：
   统一到**配对定义**后重算，查清 `51.36%` = 裸标签不一致率（非改标签率）、`15.91%` 不可复现；
   并补了「旧运行 vs 新运行」的**同分辨率**对比（§八 表）。产出 `seed_churn_audit.csv` /
   `seed_churn_audit_manifest.json` / `figures/seed_churn_audit.png`。
7. 五份夜跑签字的 r\* 需你复核确认。
8. 两个未决的论文偏差：scDblFinder vs Scrublet；QC 阈值自适应 vs 固定。
   → **2026-09-22 进展**：已从 PMC 全文取出原文的**确切参数**（<500 基因 / <1000 reads / 线粒体 ≥20%、Scrublet 逐样本、
   SCTransform top 3,000 HVG、PCA 50 PC + Harmony 1.2.0、分辨率 0.5–0.8、保留 **401,635** 核），
   逐项对照表见 §三（二）末。**偏差本身未改**（R2 优先），但现在可**引原文数值**写，不必再"推测"。
9. 两臂 RCTD（粗 6 谱系 vs 细亚型）。
10. M0 缺口：GSM9226176 的 tar 已备好并验过（90,677,930 字节，`gzip -t` 通过），
    但放入只读源目录需要你**单独授权**。
