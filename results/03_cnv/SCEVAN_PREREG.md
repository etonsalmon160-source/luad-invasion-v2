# SCEVAN 第二臂 —— 预注册（2026-09-22，**执行前**登记）

> 法则 3.1：参数必须显式常量化并登记。法则 3.2：阈值一经确定不得事后调整。
> 本文件在**首跑之前**写成。下面的输入定义、守卫与判读此后不得改动；
> 要改必须另开一版，写明改了什么、为什么、改之前跑出了什么。
>
> 本文件与 [`ANCHOR_PREREG.md`](ANCHOR_PREREG.md)（CopyKAT 臂）**并列**。该文件 §7 写明
> 「SCEVAN 臂的运行与判读**另开预注册**」——本文件即那一份。

---

## 0. 本臂是什么，不是什么

**是**：用户 **2026-09-17 签字**的双臂互证路线
（[`docs/PARAMETERS_AND_SOURCES.md`](../../docs/PARAMETERS_AND_SOURCES.md) §M2 末行）的
**第二条臂**——把**同一个显式锚定**喂给第二种 CNV 工具，报两臂一致性。

> 签字原文：「**双臂互证：CopyKAT ＋ SCEVAN**，喂**同一个显式锚定**
> （CopyKAT 走 `norm.cell.names`，SCEVAN 走 `norm_cell`），报两者一致性。」

**不是**：不是新路线，不是新方向，**不重新讨论**锚定来源（同患者 Normal 同亚型）。
方向由用户 2026-09-17 签字；本文件只把该方向**落到可执行的常数**上。

⚠️ **本臂仍不产生 M2 的正式恶性判定，不命名肿瘤细胞**（§6）。它检验的是
**两把独立的尺子量出的"恶性/正常"是否指向同一批细胞**——是**一致性检验**，
不是"谁是肿瘤"的裁决。裁决仍须 M2 门独立签字。

---

## 1. 要回答的唯一问题

**在同一个显式锚定下，SCEVAN 与 CopyKAT 对"哪些细胞是恶性"的判断一致到什么程度？**

三条判据（缺一不可）：

1. **两臂在**同一**输入上运行**（同患者、同亚型、同细胞集合、同锚定）——否则一致性无意义；
2. **两臂都产出了非退化的判决**（各自都有恶性细胞、都有非恶性细胞；全 0 或全 1 的臂不进入一致性统计）；
3. **逐细胞一致性 ≥ 预注册门**（§5 表），且不一致的细胞**逐条列表上报**（不自动裁决）。

### 为什么非做不可

CopyKAT 臂的**结构性局限已被实测**（`GP1_report.md` §9.7–9.9）：`copykat.R:456`
是**二值分支**——`cor >= 0.6` 则全判二倍体，否则按 Wasserstein 把**每个**簇强制分
diploid/aneuploid，**没有第三档**；同一份正常肺被亚型劈开（AT2 判非整倍体、AT1 判二倍体）。
⇒ 单一工具的判决**有可能是工具自身的聚类切分**，而非 CNV 信号。

**第二把尺子**是唯一的区分办法：若两把尺子在**同一输入**上给出同一批恶性细胞，
则判决更可能是信号；若两把尺子各指一批，则至少有一把在量自己的切分。

⚠️ **SCEVAN 有同构的局限**：`classifyTumor.R:285` 是 `cutree(hcc, 2)`——**也是强制的二聚类**。
⇒ **本臂不可能修好二值问题**，它只能**暴露**二值问题在两把尺子间是否一致。
这一点必须在报告里写明，**不得**把"两臂一致"说成"方法无偏"。

---

## 2. 输入定义（逐字 —— 与 CopyKAT 臂**逐字相同**）

### 2.1 细胞来源（同一份签字冻结件）

上皮细胞集合 = **两标准交集 128,091 核**（用户 2026-09-17 签字，
`results/05_annotation/epiCNV_subset_barcodes.txt`，
sha256 `db5368c22e1a7ec5c4557660e0f45ea1336ed6b6d3ca395d27daac7de36d9784`）。
**该文件已签字冻结，本臂不得修改、不得重新生成。**
（能挂上亚型标签的实为 127,654；437 个挂不上的**剔除并计数**，不填补 —— 与 CopyKAT 臂 §2.2 同。）

### 2.2 亚型标签

与 CopyKAT 臂 §2.2 **同源**：GP8a 上皮 L2 的 `harmony_res0.7_seed0` 列 →
`results/05_annotation/epiA_cluster_annotation.csv` 的 `argmax` 映射。
**主口径 AT2、次口径 AT1**（CopyKAT 臂 §2.4）。四份输入件哈希见 §9.1。

### 2.3 一次运行的细胞集合（**与 CopyKAT 臂逐字相同** —— 一致性的前提）

对**每个患者 × 每个病灶样本**各成一次运行：

```
count_mtx  = { 该患者 Normal 样本的 <亚型> 细胞 }  ∪  { 该病灶样本的 <亚型> 细胞 }
norm_cell  = { 该患者 Normal 样本的 <亚型> 细胞 }        ← 显式锚定，非 NULL / 非空串
```

- **同质**：矩阵里只有**一个亚型**（AT2 或 AT1）。
- **锚定**：基线不由 SCEVAN 猜。`pipelineCNA.R:43` `normalNotKnown <- length(norm_cell)==0`
  ⇒ 我们给了显式 `norm_cell` ⇒ `normalNotKnown = FALSE` ⇒ `findConfident = FALSE`
  ⇒ **不触发** SCEVAN 的自动正常细胞搜索（该搜索会引入它自己的选择），直接用我们给的锚定。
- **不合并不同患者**；同一次运行内只含同一个人的两个样本。
- **地板 840 保留**（同 CopyKAT 臂 §3，作用与代价同）：逐患者上报地板前/后两个细胞数。

### 2.4 亚型、分期、患者范围

与 CopyKAT 臂 §2.4 **逐字相同**：主口径 AT2 / 次口径 AT1；23 患者**全部跑，不挑选**；
AAH/AIS **照跑但单列低置信度**（"stromal cells have diploid genomes" 前提同样无癌前背书）。
**AT1 臂的意义同 CopyKAT 臂**：AT1 基本不是 LUAD 起源 ⇒ 天然的**近阴性对照**。

---

## 3. 参数：冻结在**实测的 `formals()` 签名**上

版本冻结：**SCEVAN 1.0.3 @ `5a49b88a`**、**yaGST 2017.8.25 @ `56227df`**，
隔离库 `/home/eto/Rlibs/SCEVAN`（见 `ANCHOR_PREREG.md` §7）。

**实测签名**（不是 README 抄来的）：

```r
pipelineCNA(count_mtx, sample="", par_cores=20, norm_cell=NULL, SUBCLONES=TRUE,
            beta_vega=0.5, ClonalCN=TRUE, plotTree=TRUE, AdditionalGeneSets=NULL,
            SCEVANsignatures=TRUE, organism="human", ngenes_chr=5, perc_genes=10,
            FIXED_NORMAL_CELLS=FALSE, output_dir="./output")
```

🔴 **参数名第二次更正**：README 写 `norm_cells`，**1.0.3 里不存在**。
`pipelineCNA` 真名 = **`norm_cell`**；`classifyTumorCells` 真名 = **`norm_cell_names`**。
照 README 字面名调用**直接报 unused argument**。

### 3.1 ✅ `FIXED_NORMAL_CELLS = FALSE`（裁定，非默认值 —— 理由在源码里）

**`FIXED_NORMAL_CELLS=TRUE` 在本臂下是退化的、必须禁用。** 源码
`classifyTumor.R:292–296`：

```r
if(FIXED_NORMAL_CELLS){
  cellType_pred <- names(hcc2)
  cellType_pred[!cellType_pred %in% norm_cell_names] <- "malignant"
  cellType_pred[cellType_pred %in% norm_cell_names] <- "non malignant"
}
```

⇒ 它把**凡不在锚定名单里的细胞一律写成 `malignant`**，**完全不看聚类、不看 CNV**。
配合我们的显式锚定，这个开关的输出是**按定义构造出来的**：
「病灶细胞 100% 恶性、锚定细胞 100% 正常」，**信息量为零**，且**恰好是我们要检验的命题本身**。

**裁定：`FIXED_NORMAL_CELLS = FALSE`**（走 `else` 分支
`classifyCluster(hcc2, norm_cell_names)`，用 `cutree` 的 2 簇 + 锚定名单来分类）。
🔴 **这是法则 0（静默默认陷阱）的一个实例：该参数默认 `FALSE`，而 README 把它描述成
"只想看克隆结构时用"——照描述理解很容易打开它，一打开结论就是废纸。**

### 3.2 其余参数

| 参数 | 取值 | 说明 |
| :--- | :--- | :--- |
| `norm_cell` | 显式锚定（§2.3） | 本臂的**核心**；非 NULL |
| `ngenes_chr` | **5** | 与 CopyKAT 臂 `ngene.chr=5` **对齐**（同一基因数门槛） |
| `perc_genes` | **10** | ⚠️ 注意：`pipelineCNA.R:45` 传的是 `perc_genes/100` ⇒ 有效 **0.1**。与 CopyKAT 的 `LOW.DR=0.05/UP.DR=0.10` **不是同一个量**，不类比、不强行对齐 |
| `SUBCLONES` | `TRUE` | 默认。**本臂只用其 normal/tumor 判决**，亚克隆树不在一致性统计内（§6） |
| `ClonalCN` / `plotTree` | `TRUE` | 默认 |
| `beta_vega` | `0.5` | 默认 |
| `SCEVANsignatures` | `TRUE` | 默认 |
| `AdditionalGeneSets` | `NULL` | 默认。**不给**——给了就是引入未登记的外部信息 |
| `organism` | `"human"` | 默认 |
| `par_cores` | **待定，见 §3.3** | 唯一未冻结的运行参数 |

### 3.3 ⚠️ `par_cores`：**尚未冻结**（因未做经验验证）

**源码级事实**（`/tmp/scevan_src/`，pinned SHA `5a49b88a`）：

- 全包**只有一处 RNG**：`classifyTumor.R:24` `rnorm(1, mean=0, sd=x)`，位于
  `removeSyntheticBaseline`，由**串行** `sapply` 调用，且 `classifyTumor.R:172 set.seed(1)`
  与 `preProcessing.R:76 set.seed(123)` **先行固定**。
- 被并行的函数（`funcCNA`、`nonLinSmooth`、`getMajorityCall`、`parDist`）**无 RNG 调用**。
- Linux 路径走 `mclapply`（fork，确定性）；`makeCluster` 仅 Windows。

⇒ **推断**：`par_cores` 是**纯速度旋钮**，不改结果。

🔴 **但这只是源码级推断，不是实测。** 残余风险明确：`parDist(threads=par_cores)`
的并行归约若存在浮点非结合性，可能改变 `hclust` 的**并列破平**，进而改 `cutree` 结果。
⇒ **冻结前置条件：必须做一次经验验证** —— 同一输入、`par_cores=20` 与 `par_cores=1`
各跑一次，逐细胞判决**逐字节比对**。一致 ⇒ 冻结（建议 20，与机器 20 核上限一致）；
不一致 ⇒ **两值都报**，并把 `par_cores` 降级为**结果相关参数**恒定为 1（放弃速度）。

⚠️ **对比提醒**：CopyKAT 的 `n.cores` 是**已实测的结果相关参数**，恒为 1。
若 SCEVAN 的 `par_cores` 验证通过而 CopyKAT 的不通过，这**不矛盾**——
两个工具的实现不同，**不得**互相类推。

### 3.4 实现约束（与既有代码的关系）

- 上游 **`preProcessing.R:163`**：`if((ncol(count_mtx_annot)-5)<15) stop("Bad sample low cells < 15")`
  ⇒ **mat 里需 ≥20 细胞**才不报错，且这是在 SCEVAN **自己的细胞过滤之前**。
  🔴 **登记：`n_total`（锚定＋病灶）< 20 的运行 → 记 `not_testable`，不填补、不放宽。**
  ⚠️ 该下限比 CopyKAT 臂的 G5（锚定 ≥10）**更紧**，两者**不必一致**，各自独立登记。
- 上游 **`preProcessing.R:156`**：`if(length(cellsFilt)==(ncol-5)) stop("all cells are filtered")`
  ⇒ 全部细胞被 SCEVAN 过滤 → 报错。处置同 §3.4 上一条：`not_testable`，**不得**改
  `ngenes_chr` 重跑（法则 3.2）。
- **不修改** `01_copykat_gse308103.R`；SCEVAN 臂是**独立脚本**（§8，待实现），
  读同一份输入件、用同一套锚定定义，**不共享 copykat 的代码路径**。
- 输出目录**同样按亚型分目录**（CopyKAT 臂 `out_dir` 只按 `sample_id` 命名，
  不分亚型会撞目录 —— 同一坑，SCEVAN 侧不得重演）。

---

## 4. 可证伪守卫（写进代码，不靠人看）

| # | 守卫 | 判据 | 越界后 |
| :--- | :--- | :--- | :--- |
| **S1** | **锚定分支确已走到** | 运行日志/返回值显示用的是给定 `norm_cell`（`normalNotKnown=FALSE`） | 缺 → **硬停**（说明锚定被静默忽略，本臂输入定义失效） |
| **S2** | **锚定纯度**（阴性对照） | 锚定细胞被判非恶性的比例 | 应 ≈ 100%。`< 0.50` → `implausible`，**停下报告**；`< 0.90` → 打旗标单列<br>⚠️ **两阈值是「本项目约定」（⚠️C）**，无文献出处；登记理由同 CopyKAT 臂 G1：必须有个**跑之前**能说清的停机点。**整条曲线一并上报** |
| **S3** | **阳性对照** | 病灶臂是否**至少有一个**非零恶性 | 全部为 0 → 结论 = **SCEVAN 在本数据集检出端亦无功效**，**停下**；**不得**降 `ngenes_chr` / 改 `perc_genes` 重跑 |
| **S4** | **细胞计数守恒** | 输入细胞数 == 输出判决行数 + SCEVAN 自己过滤掉的数 | 不等且无法归因 → **硬停**（静默丢细胞 = 判决分母不明） |
| **S5** | **非退化判决** | 每臂都须**同时**含 `malignant` 与 `non malignant` | 全 1 或全 0 的臂 → **不进一致性统计**，单列上报（退化臂的一致性数字无意义） |
| **S6** | **`FIXED_NORMAL_CELLS` 确为 FALSE** | 运行记录里该参数 == `FALSE` | 若为 TRUE → **整个结果作废**（§3.1：输出是定义构造的，信息量为零） |

🔴 **S6 是本臂最容易出的错**，因为它**不会报错、不会告警**，只会给出一份"完美"的结果
（病灶 100% 恶性）。⇒ 必须**断言**，不能靠跑的人记得。

---

## 5. 预注册判读（**执行前**定，防事后择果）

### 5.1 一致性怎么算

- **单位**：逐**细胞**（不是逐簇、不是逐样本）—— 因为两把尺子的簇结构不可比。
- **两臂判决对齐**：CopyKAT `prediction.txt` 的 aneuploid/diploid ↔ SCEVAN 的
  `malignant`/`non malignant`。⚠️ **该对齐是「本项目约定」（⚠️C）**：两个工具都不产出
  "恶性"这一生物学术语，此处是把各自的**二元非整倍体判决**视为同一命题。
  **登记理由**：这是"双臂互证"唯一可操作的形式；对齐本身是一个**未验证的假设**。
- **指标**：逐细胞 **Cohen's κ** + 原始一致率 + 2×2 表（双阳性 / CopyKAT-only / SCEVAN-only / 双阴性）。
- **判读门**（比照 `PLAN_AND_CHECKPOINTS.md` GP6 的双标准注释切点，**不新造**）：

| 观察 | 结论 |
| :--- | :--- |
| κ ≥ 0.80 且一致率 ≥ 90% | 两臂**互证**：判决更可能是 CNV 信号而非单一工具的切分。**仍不构成 M2 门通过** |
| 0.60 ≤ κ < 0.80 | **强制人工复核**：逐条看不一致的细胞（各自 CNA 图谱 + marker 佐证），**不自动裁决** |
| κ < 0.60 | **停**。结论 = 至少一把尺子在量自己的切分；**不得**挑一把当"对的" |
| 某臂退化（S5 不过） | 该运行**不计算一致性**；只有一臂可用 ⇒ 结论 = **双臂互证失败**，如实上报样本量 |
| S1/S4/S6 任一不过 | 该运行**无效**，逐条上报，不进统计 |

🔴 **不设"低于某数值就算修好"的验收阈值**（法则 3.2）。除 S2/S3 的停机点与上表的
κ 切点（**沿用既有 GP6 切点，非新造**）外，本臂只做**一致性归因**，不做验收。

### 5.2 若两臂不一致，下一步是**预先写定**的

- 不一致细胞清单 + 各自证据 → **交人工复核**（不自动裁决，同 GP6 纪律）。
- **不得**回头调 `ngenes_chr` / `perc_genes` / `KS.cut` / `LOW.DR` / 地板 / 亚型定义来凑一致。
- 第三把尺子 **inferCNV (R)**（须强制指定参考，**本机未装**）**不在本臂范围**；
  若要上，**另开预注册**。

---

## 6. 本臂**不**产生的结论

- **不**产生 M2 的正式恶性判定，**不命名肿瘤细胞**。
- **不**把"两臂一致"说成"方法无偏"——两把尺子都是**强制的二值切分**（§1），
  一致只排除"各量各的切分"，不排除"共同的结构性局限"。
- **不**使用亚克隆结构（`SUBCLONES=TRUE` 只是默认；树、克隆数**不进**一致性统计，只作附录）。
- **不**用 SCEVAN 输出推翻或替代 CopyKAT 臂的 §5 判读；两臂结论**并列**报告。
- **AAH / AIS 的恶性判定无文献背书** ⇒ 只作**低置信度**上报（同 CopyKAT 臂 §6）。
- **不**宣称因果、**不**宣称"管线修好了"。

---

## 7. 两个诚实边界（**必须随结果一同交付**）

1. 🔴 **SCEVAN 的运行成本尚未测量。** 没有任何一次 SCEVAN 运行发生过
   （`ANCHOR_PREREG.md` §7：「只装了，没跑」）。⇒ 挂钟、内存、磁盘**全是未知数**，
   **不得**用 CopyKAT 的成本模型外推（不同工具、不同算法）。
   **前置动作：先跑一次 P11 冒烟并实测三量，再谈排期。**
2. 🔴 **`par_cores` 的结果无关性目前只是源码级推断**（§3.3），**尚未经验验证**。
   冻结前必须做 20-vs-1 的重复跑比对。

---

## 8. 复现命令（**待实现 —— 本文件不声称脚本已存在**）

计划（尚未写成脚本）：

```r
# 03_cnv/<NN>_scevan_anchor.R  —— 待写
# 输入：results/03_cnv/<亚型>/<样本>/ 的同源细胞集合（与 CopyKAT 臂同定义）
# .libPaths(c("/home/eto/Rlibs/SCEVAN", .libPaths()))
# pipelineCNA(count_mtx, sample=<sample_id>, par_cores=<§3.3 待冻结>,
#             norm_cell=<显式锚定>, ngenes_chr=5, perc_genes=10,
#             FIXED_NORMAL_CELLS=FALSE, SUBCLONES=TRUE, output_dir=<按亚型分目录>)
```

⚠️ **脚本、冒烟、运行都还没做。** 本文件只登记口径。

---

## 9. 输入件哈希（与 CopyKAT 臂 §9.1 同源，本臂运行时硬断言）

| 产物 | sha256 |
| :--- | :--- |
| `results/05_annotation/epiCNV_subset_barcodes.txt` | `db5368c22e1a7ec5c4557660e0f45ea1336ed6b6d3ca395d27daac7de36d9784` |
| `results/04_integration/seurat_trad/epiA/clusters.csv.gz` | `13e58cb384d2af31cd740d46bfc4516fb0078a1efe1306a25682a7f18bde566f` |
| `results/05_annotation/epiA_cluster_annotation.csv` | `589b69cc03c89afc991d18baaa19fecc935325d11b58f2dccddc2a55ab7c8e4e` |
| `results/01_qc/gse308103_per_cell_qc.csv.gz` | 见 `results/03_cnv/epi_full_joblist.tsv` 生成器的 `inputs_sha256` |

⚠️ **只有第一份是用户签字的冻结件**；其余为 GP8a / M1 产物，此处冻结是为**本臂可复现**，
**不改变它们各自的签字状态**。

---

## 10. 待用户签字 / 待办（**只有这几项**）

| # | 项 | 我的建议 | 状态 |
| :--- | :--- | :--- | :--- |
| 1 | **SCEVAN 何时跑** | **不**与 CopyKAT 全量批跑同时跑：SCEVAN 默认 `par_cores=20` 会**吃满全部 20 核**，与 CopyKAT 的并发 8 抢核，两边的挂钟都变差且不可预测 | ⬜ 待定 |
| 2 | **`par_cores` 冻结值** | 先做 §3.3 的 20-vs-1 重复跑验证（成本 ≈ 一次小样本的 2 倍），再定 | ⬜ 待验证 |
| 3 | **`FIXED_NORMAL_CELLS=FALSE`**（§3.1） | **强制 FALSE**，理由在源码里（TRUE 的输出是定义构造的）；**须用户确认这一裁定**，因为它涉及"本臂产不产出正常/肿瘤分类" | ⬜ 待确认 |
| 4 | **单细胞判决的二分对齐**（§5.1 的对齐假设） | 登记为未验证假设，随结果交付 | ⬜ 待确认 |
| 5 | 脚本实现（§8） | 冒烟跑通后再实现全量 runner | ⬜ 待办 |

---

登记时间：**2026-09-22，SCEVAN 首跑之前。**
上游依据：[`ANCHOR_PREREG.md`](ANCHOR_PREREG.md) §7（双臂互证状态表）、
`GP1_report.md` §9.7–9.9（CopyKAT 二值结构实测）、
[`docs/PARAMETERS_AND_SOURCES.md`](../../docs/PARAMETERS_AND_SOURCES.md) §M2（2026-09-17 签字路线）。
