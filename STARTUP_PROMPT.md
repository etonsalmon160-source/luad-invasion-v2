# 启动提示词 (STARTUP PROMPT) — LUAD v2 · 配对空间-单核图谱

> 用途：开新会话时粘贴本文件（或其内容），即可在**正确范围**下接手。
> **范围（2026-09-12 收窄）**：**仅两个配对数据集**。

---

## 提示词正文（复制以下内容）

```
你接手 LUAD v2 项目：用**配对的空间+snRNA**数据刻画肺腺癌癌前→浸润轴。

【工作目录】 /home/eto/luad_v2/
（旧目录 /home/eto/luad_invasion 仅作【只读】数据源；范围外存档在其 luad_v2_out_of_scope/）

【先读（必读，再动手）】
1) PLAN_AND_CHECKPOINTS.md      —— 里程碑 + 硬性过门条件（含 §5b 环境约束）
2) docs/WHITEPAPER.md           —— 技术路线 v2
3) docs/PROJECT_SUMMARY.md      —— 已核实事实（含实测数据特征、方法学事实）
   （参数出处：docs/PARAMETERS_AND_SOURCES.md；靶点手册：docs/M7B_MR_COLOC_TARGET_ANCHORING.md）

【范围：只有两个配对数据集】
- GSE308103 —— snRNA（细胞核，FFPE；75 样本 / 798,100 核）→ 单细胞**参考**；唯一含 AAH
- GSE307534 —— Visium 空间（FFPE；本地 56 切片）→ **解卷积对象**
- 23 例配对患者：P3–P25
- ⚠️ Visium spot 是多细胞混合 → 必须解卷积；"配对"指同患者/同病灶，不取消解卷积
- 其余数据集（GSE131907/189357/148071、GSE190811、HRA001130…）**一律不在范围**

【铁律，违反即停】
R1 数据身份以 GEO/GSA 为准；分期**无静默默认**（严禁 .get(x,'IAC') 之类回退）
R2 恶性标签**须 CNV 证真**（CopyKAT）；双体用 **scDblFinder**
R3 无真实来源 = 不计算（不伪造、不硬编码、不用 np.random、不假生存曲线）
R4 产物可复现 + 有哈希；patient_id（真患者）与 sample_id（切片/样本）分层
R5 措辞：**候选 / 遗传学支持的候选**；**不得称因果**（观察性数据不能建立因果）

【已核实的关键事实（勿再犯）】
- GSE308103 是 **snRNA（细胞核）**：median nCount 1,516、median mt% 0.6 → 阈值**不可照搬整细胞 scRNA 的直觉**（盲用会砍掉约 30% 的核）
  ⚠️ **2026-09-16 口径已变，勿再引用旧说法**：分析掩膜 = **M1 自适应 MAD（`pct_mt<5`）∩ 源论文绝对门（`nFeature≥500 & nCount≥1000 & pct_mt≤20`）**
  = **413,697 核 × 18,069 基因**（旧 648,945 核 × 18,082 基因保留为**敏感性臂**）。链：798,100 → 767,839(M1) → 648,945(∩单细胞) → 论文门 555,480 → 交集 413,697
- GSE190811 经核实是**乳腺癌**（非 LUAD LNM）→ 已废
- **SCMG 无「状态逆转/因果」能力**（源码无此方法；参照流形无肿瘤态）→ 仅作整合/流形对照
- 观察性单细胞不能建立因果；合法杠杆 = **cis-MR + coloc**（M7b）
- RCTD full 模式无 reject 类别、constrain=F → M5 门不得用自归一化 Σ=1
- squidpy.nhood_enrichment 置换标签、只返回 z、**无 p 值**
- 单细胞基础模型打不过 scVI/Harmony 基线；PLIP zero-shot 判 WHO 生长模式不成立（spot 尺度不足）
- CMap 的 Tau 是 0–1 重现性指标，**不可能为负**
- 本机：**无 GPU**、**Python 3.8**、共享库 root 属主（装包进个人库）、出网须绕死代理

【环境】
- 主流工具源码：/home/eto/luad_v2/tools/
- 已装：R 侧 CopyKAT/coloc/ieugwasr/TwoSampleMR/scDblFinder/spacexr；Python 侧 harmonypy/SpaGCN/plip/vina
- 硬约束详见 PLAN §5b

【当前位置】(2026-09-17) M3-A 传统分支**全量降维已完成，GP4a/4b/4c/GP5 全部签字**：
- **GP4c 批次**：用户裁定**采用源论文的 Harmony 主口径**（`RunHarmony(group.by.vars="sample_id", dims.use=1:50)`）；
  不校正臂保留为**敏感性/局限陈述**，不作主结果。
  ⚠️ **约束性限定（须随下游结果一并陈述）**：Harmony 抹掉约 **60%** 分期可恢复信号（0.5388→0.4866，基线 0.4546），
  且 `sample_id`↔`stage` 嵌套 ⇒ **被抹掉的部分无法归因** ⇒ 分期相关结论**不得声称已排除批次混淆**。
- **GP5 分辨率**：签字 **`r* = 0.6`**；**指标4（AAH 吸收护栏）经用户签字承认在全量上无效**（分母为空），`r*` 由指标 1/2 承担。
- **GP0** 按论文 QC 重做（413,697 核）🔶 已重做、2026-09-17 审计补做门校验（13/0）、**待签字**
- 脚本 `04_integration/10_seurat_traditional.R`；报告 `results/04_integration/seurat_trad/full/GP4_report.md`（§九 裁定 / §十 签字后状态）
- **下一步 = GP6 双标准注释**（法则2 marker ⊕ CellTypist）+ Cohen's κ，口径 = **Harmony 主口径的 `r*=0.6` 聚类**

【汇报风格】诚实优先：报错要报得清楚；不确定就说不确定；宁可拒绝，不伪造。
```

---

## 备注
- 本提示词与 `PLAN_AND_CHECKPOINTS.md` 配合使用最完整。
- 发 `github.com` 被墙（`api`/`codeload` 通）；PyPI 走清华镜像。
