# LUAD v2 —— 配对空间-单核图谱（癌前→浸润轴）

> **范围（2026-09-12 收窄）**：仅**两个配对数据集**——`GSE308103`(snRNA) + `GSE307534`(Visium 空间)。
> 与旧目录 `/home/eto/luad_invasion`（结果作废）**物理隔离**；该目录仅作**只读数据源**。

---

## 项目一句话

用**同一批患者、同一病灶、相邻切片**的 **Visium 空间**与 **snRNA** 配对数据，
刻画肺腺癌 **Normal → AAH → AIS → MIA → IAC** 的细胞状态与空间生态位，
并产出**遗传学锚定的候选靶点**（→ 结构对接）。

**技术定位**：贡献是**数据独特性 + 方法严谨性（基准与验证）**，不是发明新的因果推断算法。

## 先读
1. **[`PLAN_AND_CHECKPOINTS.md`](PLAN_AND_CHECKPOINTS.md)** —— 里程碑与**硬性过门条件**（含 §5b 环境约束）
2. [`docs/WHITEPAPER.md`](docs/WHITEPAPER.md) —— 技术路线（v2）
3. [`docs/PROJECT_SUMMARY.md`](docs/PROJECT_SUMMARY.md) —— 已核实事实总纲（v2）
4. [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md) —— 参数出处（v2）

## 铁律
1. 数据身份以 **GEO/GSA** 为准；分期**无静默默认**；
2. 恶性标签**须 CNV 证真**（CopyKAT）；双体用 **scDblFinder**；
3. **模态混淆（sn/空间）显式处理**，禁 naive 合并；
4. 无真实来源 = 不计算（宁可报缺，不伪造）；
5. 措辞：**候选 / 遗传学支持的候选**，**不得称因果**。

---

## 数据底座（权威）

| 数据集 | 身份 | 模态 | 角色 |
| :--- | :--- | :--- | :--- |
| **GSE308103** | 75 样本 / **798,100 核**（实测） | **snRNA**（FFPE） | 单细胞**参考**；**唯一含 AAH** |
| **GSE307534** | GEO 56 样本 / 25 患者；本地 19 切片 | **Visium spot**（FFPE） | **空间图谱**；解卷积对象 |

**9 例配对患者**：`P3 P4 P10 P13 P15 P18 P21 P22 P25`（本地空间切片完整覆盖）。
**LNM**：空转暂缺（`GSE190811` 经核实为**乳腺癌**，已废）。

> ⚠️ Visium spot 是多细胞混合 → 必须用单细胞参考**解卷积**；"配对"指**同患者/同病灶**，不取消解卷积。

---

## 目录
```
luad_v2/
├── PLAN_AND_CHECKPOINTS.md   # ← 先读：计划 + 严格检查点
├── README.md                 # 本文件
├── docs/                     # WHITEPAPER / PROJECT_SUMMARY / PARAMETERS / M7B手册 / 规则
├── 00_ingest/                # 权威队列登记（cohort_registry.py）
├── 01_qc/                    # M1：QC + 双体（GSE308103）
├── tools/                    # 主流工具源码（不入库）
├── results/  logs/           # 产物与日志
└── data/ scmg/ traditional/ scripts/   # 预留（空）
```

## 范围外内容
三个 scRNA 队列（GSE131907/189357/148071）、HRA001130 接口、旧 M0/M1 产物 →
已移至 `/home/eto/luad_invasion/luad_v2_out_of_scope/`（**未删除**，见其 README）。

## 作废说明
旧 `luad_invasion/` 的全部结果产物（图表/数值/TMB/CMap/对接）**作废**。
一切从本目录、按 `PLAN_AND_CHECKPOINTS.md` 重新开始。
