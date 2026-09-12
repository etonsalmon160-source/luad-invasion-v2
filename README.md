# LUAD v2 —— 从头重做的干净工程根目录

> **最外层新建**，与旧目录 `/home/eto/luad_invasion`（**结果作废**）**物理隔离**。
> 旧 `luad_invasion/` 仅作**只读数据源与历史参考**，不再作为产出地。
> 创建/重置：2026-09-12

---

## 状态：🔴 **进度已重置（0%）** —— 按检查点重做
- **计划指导 + 严格检查点**：**[`PLAN_AND_CHECKPOINTS.md`](PLAN_AND_CHECKPOINTS.md)**（每个里程碑设**硬性过门条件**）。
- 白皮书：[`docs/WHITEPAPER.md`](docs/WHITEPAPER.md)（已按 GEO 真值 + 双分支更新）。
- 事实总纲：[`docs/PROJECT_SUMMARY.md`](docs/PROJECT_SUMMARY.md)；参数出处：[`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md)。

## 铁律
1. 数据身份以 GEO/GSA 为准；分期**无静默默认**。
2. 恶性标签**须 CNV 证真**；双体用 **scDblFinder**。
3. **模态混淆（sn/sc）显式处理或隔离**，禁 naive 合并。
4. 无真实来源 = 不计算（宁可报缺，不伪造）。
5. **SCMG 分支不掺传统算法**；与标准分支**并行对照**。

---

## 数据底座（权威）
| 数据集 | 身份 | 模态 | 阶段 |
| :--- | :--- | :--- | :--- |
| GSE131907 | LUAD scRNA（208,506/58 样本/44 患者） | scRNA | 按 `Sample_Origin` |
| GSE189357 | 早期 LUAD scRNA | scRNA | `TD1/2/9=IAC`、`TD3/4/6=MIA`、`TD5/7/8=AIS`（**无 AAH**） |
| GSE148071 | **Advanced NSCLC**（42 样本） | scRNA | Adv_NSCLC |
| GSE307534 / GSE189487 / GSE190811 | 空间（含真实 AAH；真 LNM） | 空间 | P0 禁令锁切片 |
| **GSE308103** | 与 GSE307534 同患者配对 | **snRNA** | 含 **AAH**（暂用，SCMG 跨平台并入） |
| （预留）HRA001130 | 全细胞 scRNA，含 AAH/AIS/MIA/IA | scRNA | **受控待申请** → [接口](00_ingest/hra001130_interface.py) |

**分期真值**：GSE131907 按 `Sample_Origin`（`nLN`=正常淋巴结，**非** LNM；脑转/胸水单列）。

---

## AAH 口径（已决策）
- **暂用** `GSE308103`（snRNA）→ **SCMG zero-shot 跨平台并入**（过**跨模态五判据**，见 PLAN M4）；
- **预留** `HRA001130`（全细胞 scRNA，GSA-Human 受控）→ 获批后**零改动替换**。

## 目录
```
luad_v2/
├── PLAN_AND_CHECKPOINTS.md  # ← 计划指导 + 严格检查点（先读这个）
├── README.md
├── docs/          # WHITEPAPER / PROJECT_SUMMARY / PARAMETERS_AND_SOURCES / 规则 / legacy
├── 00_ingest/     # 权威队列登记 + HRA001130 接口（第一步）
├── traditional/   # 标准/主流分支
├── scmg/          # 纯 SCMG 分支
├── tools/         # 18 个主流工具源码（不入库）
├── data/  results/  logs/  scripts/
```

## 作废说明
旧 `luad_invasion/` 的全部结果产物（图表/数值/TMB/CMap/对接）**作废**；旧报告已隔离至其 `archive/legacy_docs/`。一切从本目录、按 `PLAN_AND_CHECKPOINTS.md` 重新开始。
