# LUAD v2 —— 从头重做的干净工程根目录

> **最外层新建**，与旧目录 `/home/eto/luad_invasion`（已完成/作废）**物理隔离**。
> 旧的 `luad_invasion/` 只作**数据源与历史参考**，**不再作为产出地**。
> 创建：2026-09-12

---

## 状态
🟡 **计划阶段（未执行）**。白皮书见 [`docs/WHITEPAPER_plan.md`](docs/WHITEPAPER_plan.md)（**计划版，无结果**）。

## 铁律
1. 数据身份以 GEO 为准；分期**无静默默认**。
2. 恶性标签**须 CNV 证真**；双体用 **scDblFinder**。
3. **模态混淆（sn/sc）显式处理或隔离**，禁 naive 合并。
4. 无真实来源 = 不计算（宁可拒绝，不伪造）。

---

## 数据底座（权威，身份以 GEO 为准）

| 数据集 | 身份 | 模态 |
| :--- | :--- | :--- |
| GSE131907 | LUAD scRNA（208,506/58 样本/44 患者） | scRNA |
| GSE189357 | 早期 LUAD scRNA（**3 AIS + 3 MIA + 3 IAC，无 AAH**） | scRNA |
| GSE148071 | **Advanced NSCLC**（42 样本） | scRNA |
| GSE307534 | 空间（含真实 AAH/AIS/MIA/LUAD） | 空间 |
| GSE189487 / GSE190811 | 空间 | 空间 |
| GSE308103 | 与 GSE307534 **同患者配对** snRNA（含 AAH） | **snRNA** |

**分期真值**：GSE189357 `TD1/2/9=IAC、TD3/4/6=MIA、TD5/7/8=AIS`；GSE131907 按 `Sample_Origin`（`nLN`=正常淋巴结，**非** LNM；脑转/胸水单列）。

---

## 待决策（唯一阻塞项）
**AAH 口径**：(A) 用 GSE308103(snRNA) 补 / (B) 单细胞不做 AAH（仅空间+bulk 签名）/ (C) 其他。

## 目录
```
luad_v2/
├── README.md
├── docs/           # 白皮书(计划版) + 事实总纲
├── 00_ingest/      # 队列纳入 + 权威分期/层级（第一步）
├── traditional/    # 传统分支（阶段 1–5）
├── scmg/           # 纯 SCMG 分支（阶段 6，正交）
├── data/           # 只读引用旧 data/（软链，不复制大文件）
├── results/  logs/  scripts/
```

## 作废说明
旧的 `luad_invasion/` 产出（所有图表/数值/TMB/CMap/对接）**全部作废**；旧总结报告已隔离至其 `archive/legacy_docs/`。一切从本目录重新开始。
