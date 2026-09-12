# 启动提示词 (STARTUP PROMPT) — LUAD v2

> 用途：开新会话时粘贴本文件（或其内容）给 AI，即可在**正确口径**下接手"从头重做"。

---

## 提示词正文（复制以下内容）

```
你将接手一个【从头重做】的肺腺癌（LUAD）多组学课题。旧工程的一切结果已作废，本工程进度归零，按严格检查点重做。

【工作目录】 /home/eto/luad_v2/
（旧目录 /home/eto/luad_invasion 仅作【只读】数据源与历史参考，不得在其上产出。）

【先读这三个文件（必读，再动手）】
1) PLAN_AND_CHECKPOINTS.md   —— 计划指导 + 严格检查点（每个里程碑有硬性"过门条件"，不过门不得进下一步）
2) docs/PROJECT_SUMMARY.md   —— 已核实的权威事实（数据真值 / 旧错误 / 整改）
3) docs/WHITEPAPER.md        —— 技术路线（双分支：标准/主流 与 纯 SCMG）
   （参数出处见 docs/PARAMETERS_AND_SOURCES.md）

【铁律，违反即停】
R1 数据身份以 GEO/GSA 为准；【分期无静默默认】（严禁 .get(x,'IAC') 之类回退）。
R2 恶性标签【须 CNV 证真】（CopyKAT/inferCNV）；双体用 scDblFinder。
R3 无真实来源 = 不计算（不伪造、不硬编码、不用 np.random、不假生存曲线）。
R4 SCMG 分支【不掺传统算法】；与标准/主流分支【并行对照】。
R5 产物可复现 + 有哈希；patient_id（真患者）与 sample_id（组织/切片）分层。

【已知的坑（旧工程，勿再犯）】
- 单细胞"AAH"是假的：GSE189357 实为 3AIS+3MIA+3IAC，【无 AAH】；旧把 TD9 误标 AAH（实 IAC）、TD4 误标 IAC（实 MIA）。
- GSE131907 分期污染：nLN(正常淋巴结) 被记成 LNM；脑转移/胸腔积液被并入 IAC。
- GSE148071 实为【Advanced NSCLC】，非"早期 LUAD"。
- 旧脚本含写死数值/伪算法（107 CMap、11 TMB-THPP、45 PAGA、13/53/60 等）→ 不得参照其数值与方法。

【AAH 口径（已决策）】
- 暂用 GSE308103（snRNA）→ 经 SCMG zero-shot 跨平台并入，须过【跨模态五判据】（PLAN M4）。
- HRA001130（全细胞 scRNA，GSA-Human 受控）已留接口（00_ingest/hra001130_interface.py），获批后零改动替换。

【第一个任务 = M0 输入冻结】
写 00_ingest/01_load_cohorts.py：按 00_ingest/cohort_registry.py 纳入三 scRNA 队列
（GEO 真值分期 + 分层），产出【冻结清单 + SHA-256 + 校验报告】。
过门条件：阶段计数与 GEO 真值一致（GSE189357: IAC=TD1/2/9, MIA=TD3/4/6, AIS=TD5/7/8, 无AAH）；
代码无静默默认；字段齐全（cell_barcode/patient_id/sample_id/stage/dataset）。不过门就停。

【工具与数据】
- 主流工具源码：/home/eto/luad_v2/tools/（scvi-tools/harmony/CopyKAT/infercnv/RCTD/PLIP/SpaGCN/Squidpy/DESeq2/CellRank/Vina/fpocket…）
- SCMG 权重（编码器+扩散+扰动库+参照流形）：/home/eto/scmg_workspace/
- 真 AAH 数据：/home/eto/luad_invasion/data/GSE308103/extracted/（75 样本）

【汇报风格】诚实优先：报错要报得清楚；不确定就说不确定；宁可拒绝，不伪造。
```

---

## 备注
- **网络**：github.com 间歇性不通；如需推送，见 README；本地 commit 不受影响。
- 本提示词与 `PLAN_AND_CHECKPOINTS.md` 配合使用最完整。
