# 空转（Visium）CNV 臂 —— 预注册（2026-09-24，**执行前**登记）

> 法则 3.1：参数必须显式常量化并登记。法则 3.2：阈值一经确定不得事后调整。
> 本文件在**首跑之前**写成。下面的输入定义、守卫与判读此后不得改动；
> 要改必须另开一版，写明改了什么、为什么、改之前跑出了什么。

---

## 0. 本臂是什么，不是什么

**是什么**：在 **Visium spot** 上调用拷贝数变异，据此判定**逐 spot 的恶性 / 非恶性**。
这是全项目**唯一**的恶性判定来源。

**不是什么**：

- **不是**单细胞 CNV 的重跑。单细胞 CNV 已由用户 2026-09-24 判定**失败**并整体退出
  （历史见 `results/03_cnv/ANCHOR_PREREG.md`、`SCEVAN_PREREG.md` 顶部横幅）。
- **不是** RCTD 的产物。RCTD 只做**构成估计与纯度 QC**，**不判恶性**。
  这是源论文自己的分工（Peng 2026 Cancer Cell 方法段：恶性由 SpatialInferCNV 判，
  RCTD 用于纯度 QC），也是用户 2026-09-24 的裁定。
- **不是**"把单细胞 CNV 结论搬到空间上"。两者不共享细胞、不共享坐标。

### 0.1 从单细胞 CNV 的失败里必须继承的两条教训

| 教训 | 具体表现（已发生的事） | 本臂的对应设计 |
| :--- | :--- | :--- |
| **锚不稳 ⇒ 全盘作废** | 同一批 3,890 正常细胞，三次跑出 **0% / 13.6% / 49.8%**；`norm.cell.names` 只算基线不参与判决；99 次里 30 次不可用 | **§4 门 SC2（锚可复算）是本臂的第一道硬门**，先于任何结论 |
| **"没跑" ≠ "没失败"** | SCEVAN 第二臂一例未跑，却差点被读成"未被证伪" | §5.3 明写：未运行的臂在本臂报告里记 **NOT RUN**，不得记作阴性 |

---

## 1. 要回答的唯一问题

> **在本队列的 Visium spot 上，能否给出一份逐 spot 的恶性判定，且该判定可复算、有可证伪的阴性对照？**

**不**回答：亚克隆结构、系统发育树、克隆与分期的关系、预后。
这些是原文的后续目标，本臂不碰（见 §6）。

---

## 2. 输入定义（逐字）

### 2.1 表达矩阵

- 来源：`data/visium_spatial/<GSM>_<P>_<期>/filtered_feature_bc_matrix/`
- **只取 `Gene Expression` 行**，剔掉 **35 行 `Antibody Capture`**
  （其计数量级可达 1e6，混入会毁掉 CNV 基线 —— 见 `SPATIAL_QC_PREREG.md` §1.1）
- 组织内 spot 集合 = `SPATIAL_QC_PREREG.md` 冻结的 `spot_mask.tsv.gz`（**待签字**）
- 合并方式：**按患者合并其全部切片**（原文按 patient 建克隆与树）

### 2.2 上皮 spot 子集 —— 🔴 **待签字**

原文只在**上皮 spot** 上跑 SpatialInferCNV。我们的上皮 spot 由 **RCTD 权重**给出
（`RCTD_PREREG.md`；病理线不得供标签，见 `SPATIAL_QC_PREREG.md` §5.3）。

**待签字项**：① 用粗版还是细版 RCTD 参考；② 上皮类目的清单；③ 权重阈值取多少；
④ 阈值未达的 spot 是剔除还是单列一档。

⚠️ **不得**在看过 CNV 结果之后再回来调这个阈值（法则 3.2）。

### 2.3 基因序（**已就绪，可复用**）

- `results/03_cnv/infercnv_smoke/gene_order_hg38.tsv`
- 来源：GENCODE **v44 (hg38)** GTF ＋ HGNC `complete_set`，经 `prev_symbol`/`alias_symbol` 解析改名
- 实测匹配率 **99.62%**（18,069 个基因里 18,000 匹配，69 未匹配）
- sha256 `1f88af4ef06710a0e60632406c7a5db8f440b1764a55e191f31c3a5c8b097194`
- ⚠️ 该文件是为 **snRNA** 建的，空转的基因名需**重新过一遍同一脚本**并**重算哈希**，
  不得直接沿用旧哈希（空转 features 里的 `HSPA14` / `TBCE` / `TMSB15B` 是**双探针**，需按 symbol 求和后取一次）

---

## 3. 方法：用 `infercnv` 实现原文记录的 SpatialInferCNV 流程

### 3.1 为什么不用 SpatialInferCNV 本体

原文用的 `SpatialInferCNV`（<https://github.com/aerickso/SpatialInferCNV>）
**仓库里没有许可证**（LICENSE 文件是模板占位符 "What license is it under?"）。
⇒ 我们**不引入**它，改为**按原文方法段逐字实现同一套流程**，跑在
`infercnv` 之上（Bioconductor，**BSD_3_clause**，许可明确）。

⚠️ **由此产生的一条登记**：本臂产出的应描述为
「**原文方法的复现**」，**不是**「原文代码的运行」。报告里不得写"用了 SpatialInferCNV"。

### 3.2 环境（**已实测可用，2026-09-24**）

```bash
export LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH
R_LIBS=/home/eto/Rlibs/infercnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library
```

| 件 | 位置 / 版本 | 实测 |
| :--- | :--- | :--- |
| JAGS | `/home/eto/local/jags/bin/jags`，**4.3.2** | ✅ 自建，无需 root |
| `rjags` | `/home/eto/Rlibs/infercnv/rjags` | ✅ 装载通过（需上面的 `LD_LIBRARY_PATH`） |
| `infercnv` | `/home/eto/Rlibs/infercnv/infercnv`，**1.23.0** | ✅ 装载通过 |
| `coda` / `RColorBrewer` | `/usr/local/lib/R/site-library` | ✅ |

⚠️ **旧记录更正**：本项目先前多处记着「infercnv 卡在 rjags ⇒ 需系统 JAGS（未装）」。
**该记录已过期。** JAGS 早已自建安装且可用 ⇒ 原文的 `HMM=TRUE` **是够得着的**
（既有的 P4 冒烟之所以没跑 HMM，是因为那份预注册当时把 `P_HMM` 设成了 `FALSE`，
不是因为环境不允许）。

### 3.3 参数 —— **逐字照抄原文，一个都不改**

原文方法段（PMC12980502，Methods「Spatial infer-CNV」）：

> ...by applying the following parameters: `cutoff = 0.1`, `cluster_by_groups = TRUE`,
> `HMM = TRUE`, and `denoise = TRUE`.

```r
# 08_spatial_deconv/<NN>_run_spatial_cnv.R
P_CUTOFF            <- 0.1     # 原文
P_CLUSTER_BY_GROUPS <- TRUE    # 原文
P_HMM               <- TRUE    # 原文 —— 本机可用（§3.2）
P_DENOISE           <- TRUE    # 原文
P_REF_GROUP         <- <待签字，见 §3.4>
P_THREADS           <- 3       # 本机 load average 长期 57–70，不拉高
```

⚠️ 与既有 P4 冒烟**不同**：那份跑的是 `P_HMM=FALSE`、`analysis_mode="subclusters"`、
参考为「同患者 Normal 全部上皮」。**本臂不继承那些取值**，按原文重设。

### 3.4 🔴 **锚（参考 spot）—— 本臂最大的缺口，必须签字**

原文选参考 spot 的**三步法**（逐字）：

> the '**without reference**' mode was first used to classify spots into different clones
> based on clustering patterns. Second, the inferred clone information for each spot was then
> mapped onto the UMAP embedding. Third, spots exhibiting the **lowest CNA signal or variation**,
> which mapped precisely to the annotated "**normal lung epithelial**" population, were selected
> as reference spots for the final round of SpatialInferCNV analysis. If no such cluster was
> identified, steps 1–3 were repeated with **finer clone clustering** until a pure normal lung
> epithelial cluster was selected.

**问题**：第三步把「最低 CNA」和「病理学家标注的正常肺上皮」**两个条件一起**要求。
前者我们能算，**后者我们没有**——病理线已判死，不得供标签（`SPATIAL_QC_PREREG.md` §5.3）。

原文另有一条**兜底**：12 例（P1,P2,P3,P4,P6,P9,P14,P15,P17,P20,P21,P24）直接用
**全部**「正常肺上皮」标注 spot 作参考 —— 这条同样依赖病理标注，我们也用不了。

#### 三个候选处置（**待用户签字，未定之前不得开跑**）

| # | 方案 | 依据 | 风险 |
| :---: | :--- | :--- | :--- |
| **A** | **纯无参**：`without reference` 模式跑完，**自动取 CNA 信号最低的克隆**作参考 | 原文第 1、3 步的前半句；**不需要任何标签** | 🔴 **循环论证风险**：无参模式下基线是**全体 spot 的均值**，若切片以肿瘤为主，最低 CNA 的克隆**可能仍是肿瘤**。单细胞臂就是这么塌的 |
| **B** | **用 RCTD 复核锚**：先按 A 取候选锚，再用 RCTD 权重检查该克隆是否由**正常上皮类**（AT2/AT1）主导；不主导则按原文"更细分克隆"重来 | 原文第三步的**意图**（要求锚是正常上皮）用**另一个模态**实现；且 RCTD 的角色仍只是**纯度 QC**，不判恶性 —— 与用户裁定一致 | 需要把 `SPATIAL_QC_PREREG.md` §5.4「CNV ↔ RCTD 无共享对象」**修订为**「共享**锚校验**这一件事，仍禁止恶性×权重的乘法」。**这是一处口径修改，必须单独签字** |
| **C** | **同患者 Normal 切片作锚**：有 Normal 切片的患者直接用其 Normal 切片的上皮 spot 作参考 | 我们队列里确有 Normal/AAH 切片；这是"同患者正常组织"的**直接**测量，不需要标签 | ⚠️ 只有部分患者有 Normal 切片（多数只有 AAH/AIS/MIA/LUAD）；**覆盖不全 ⇒ 患者间锚定义不一致**，正是单细胞锚定臂出问题的地方 |

🔴 **建议 A＋B 组合**（B 作为 A 的守卫），但**这是建议不是裁定**。
**三个方案未签字之前，本臂不跑。**

---

## 4. 门 **SC0–SC4**

| # | 守卫 | 判据 | 越界后 |
| :--- | :--- | :--- | :--- |
| **SC0** | **底座核对** | 基因序匹配率 ≥ 99% 且**重算哈希**；`Gene Expression` 行数 == 18,085 且**抗体行 35 已剔**；spot 数 == `spot_mask` 的 `n_pass` | 任一不等 → **硬停** |
| **SC1** | **参数已冻结且 HMM 真在跑** | `cutoff/cluster_by_groups/HMM/denoise` 四个值与 §3.3 逐位相同；**日志中必须出现 rjags 实际被调用的证据**（不得静默降级成 `HMM=FALSE`） | 对不上或静默降级 → **硬停** |
| **SC2** | 🔴 **锚可复算**（**本臂第一道硬门**） | 同一锚集合、同一输入，**重复 N≥3 次**，逐 spot 恶性判定的**一致率 ≥ 95%**，且**每次的锚细胞/spot 名单逐字节相同** | **不一致 → 硬停，本臂到此为止，如实上报"锚不可复算"**（这正是单细胞臂的失败模式，不许绕过） |
| **SC3** | **阴性对照** | ①锚内 spot 自己必须判为非恶性（自洽）；②锚内**留出**一部分 spot 不进参考，其恶性率须与进入参考的部分**无显著差** | 锚自己判成恶性 → 锚定义错 → **硬停** |
| **SC4** | **逐切片/逐患者上报，不设通过率阈值** | 每患者报：送进的上皮 spot 数 / 锚 spot 数 / 恶性 spot 比例 / 各染色体均值；**并单列"未定义"档** | **不硬停**；**不许**事后挑患者（法则 3.2） |

⚠️ **SC4 的"未定义"必须单列**：`infercnv` 与 CopyKAT 一样会有判不了的档。
**不得**把"未定义"并进"非恶性"——那是把测不出当成测到了。

---

## 5. 与其它臂的联动

### 5.1 依赖 RCTD（单向）

```
Visium 表达 ──[SPATIAL_QC]──→ spot 掩码
                  │
                  ├─[RCTD]──→ 逐 spot 细胞类型权重 ──→ 上皮 spot 子集 ──┐
                  │                                                    ↓
                  └───────────────────────────────────→ [infercnv] ──→ 逐 spot 恶性
                                                                        │
                  病理线（PLIP）──→ 仅切片级组织密度，**不得**供标签 ────┘
```

**要点**：CNV 吃 RCTD 的**上皮 spot 子集**（一个集合），**不吃** RCTD 的恶性判断（RCTD 没有）。

### 5.2 硬禁令（继承 `SPATIAL_QC_PREREG.md` §5.2 裁决 5）

🔴 **禁止**把「CNV 的逐 spot 恶性判定」乘「RCTD 的逐 spot 上皮权重」得到"恶性细胞数/恶性比例"。
理由不变：那是**跨定义的量纲外推**。CNV 的判定单位是 **spot**，RCTD 的权重单位是 **spot 内的比例**；
相乘得到的东西**既不是 spot 数也不是细胞数**，没有对应实体。

### 5.3 与已废弃的单细胞 CNV 臂的关系

- 单细胞 CNV 的三条线（CopyKAT 锚定臂 / SCEVAN 臂 / infercnv 冒烟）**全部退出**，
  在本臂报告里一律记 **RETIRED / NOT RUN**，**不得**记作阴性证据、**不得**与本臂结果合并。
- 唯一被继承的是**教训**（§0.1）与**工具链**（§2.3 的基因序、§3.2 的环境）。
- ⚠️ 既有 `results/03_cnv/infercnv_smoke/`（P4，snRNA，`HMM=FALSE`）的
  **AUROC 0.6524「弱信号」**是**单细胞**上的数，**不得**当作本臂的先验或预期值。

---

## 6. 本臂**不**产生的结论

- 不产生**亚克隆结构 / 系统发育树**（原文后续目标，本臂不碰）
- 不产生**克隆与分期的关系**（同上）
- 不产生**逐细胞**的恶性判定（本臂的单位是 **spot**）
- 不产生**空间邻域 / 生态位**结论（属 M6 及以后）
- 不产生**预后 / 生存**结论
- 不把**未定义**档读成"非恶性"
- 不把本臂结果与**已废弃的**单细胞 CNV 结果做一致性比较

---

## 7. 待签字（**只有这几项**）

| # | 事项 | 位置 | 现状 |
| :---: | :--- | :--- | :--- |
| 1 | **锚的处置方案 A / B / C**（建议 A＋B） | §3.4 | 🔴 **未定 ⇒ 本臂不跑** |
| 2 | **上皮 spot 子集定义**（粗/细 RCA 参考、类目清单、权重阈值、未达标 spot 的处置） | §2.2 | 🔴 **未定** |
| 3 | **`spot_mask` 的三个阈值**（`MIN_UMI_SPOT` / `MIN_GENES_SPOT` / `MAX_MT_FRAC`） | `SPATIAL_QC_PREREG.md` §3.2 | 🔴 **未定**（建议照原文 500 / 200 / 0.15，并接受 4.97% 而非原文 1.05%） |
| 4 | **若选 B：修订 §5.4「CNV ↔ RCTD 无共享对象」为「共享锚校验」** | §3.4 方案 B；`SPATIAL_QC_PREREG.md` §5.4 | 🔴 **未定** |
| 5 | **是否纳入空转 QC 的 6 张浅切片**（剔除率 15–30%，见 `SPATIAL_QC_PREREG.md` §3.2(c)） | 输入定义 | 🔴 **未定** |

---

## 8. 复现命令（**签字前跑不了**）

```bash
# 前置：SPATIAL_QC_PREREG.md 第 1 步（spot 掩码冻结）+ RCTD_PREREG.md（上皮权重）
export LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH
R_LIBS=/home/eto/Rlibs/infercnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
  Rscript 08_spatial_deconv/<NN>_run_spatial_cnv.R
# ⚠️ <NN> 脚本尚未写。本文件不声称脚本已存在（与 ANCHOR_PREREG.md §8 同一纪律）
```

---

## 9. 本文件登记时**尚未**核实的事项

1. 空转的基因名过 GENCODE/HGNC 后的实际匹配率（snRNA 是 99.62%，空转未测）
2. `infercnv` 1.23.0 在 `HMM=TRUE` 下的**内存峰值与单患者耗时**（**未实测**）——
   本机 load average 长期 57–70，且后台有 PLIP 任务在跑 ⇒ 排期前必须先测一个患者的成本
3. 上皮 spot 子集的**规模**（若 RCTD 上皮 spot 数过少，infercnv 的 HMM 可能不稳）
4. 原文那 12 例"直接用全部正常肺上皮"的兜底是否暗示**其余患者也确实能找到纯正常克隆** ——
   我们无从知晓，**不得**假设我们的患者也能

---

## 10. 本文件**不**声称

- 不声称 SpatialInferCNV 已被安装或将被安装（§3.1 明确不引入）
- 不声称本臂会成功（单细胞臂失败了；本臂的锚问题**同样存在**，只是换了个模态）
- 不声称 `HMM=TRUE` 的结果一定比 `HMM=FALSE` 好
- 不声称原文的 1.05% spot 剔除率能被我们复现（实测 4.97%，见 `SPATIAL_QC_PREREG.md` §3.2）
