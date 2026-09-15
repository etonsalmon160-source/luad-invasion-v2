# 项目事实总纲 (PROJECT_SUMMARY)

> **定位**：本项目**唯一的最新事实总纲**（authoritative facts ledger）。
> **范围**：仅**两个配对数据集**——`GSE308103`(snRNA) + `GSE307534`(Visium 空间)。
> **更新日期**：2026-09-12（v2，随范围收窄重写）
> 配套：[`WHITEPAPER.md`](WHITEPAPER.md)（技术路线）· [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md)（检查点）

---

## 1. 数据底座（GEO 逐样本核实）

| 数据集 | 真实身份 | 角色 | 实测规模 |
| :--- | :--- | :--- | :--- |
| **GSE308103** | snRNA-seq（FFPE fixed RNA） | 单细胞**参考**；**唯一含 AAH** 的单细胞资源 | **75 样本 / 798,100 核**（实测） |
| **GSE307534** | Visium CytAssist FFPE 空间 | **空间图谱**（原位坐标）；解卷积对象 | GEO 56 样本/25 患者；**本地 56 张切片** |

**配对（23 例，P3–P25）**：两模态均有切片；本地空间切片**完整覆盖**这 23 例（P1/P2 仅空间、无 snRNA，不入配对）。
> **2026-09-15 修正**：旧值"9 例"是**仅下载 19/56 张切片**时的交集产物；切片补齐后按两张 GEO 权威表求交集实为 23 例。
分期（两模态一致）：**Normal / AAH / AIS / MIA / IAC**。

**LNM**：空转淋巴结转移**无合法 LUAD 数据**。`GSE190811` 经 GEO 核实系列标题为
*"…paired metastatic lymph node tumors in **breast cancer** patients"*（**乳腺癌**），
且白皮书 v1 所锁 GSM5732148 在该库**不存在**（真实为 GSM5732357–5732360）→ **该切片作废**，禁止用于任何 LUAD 产物。

**已移出范围**：GSE131907 / GSE189357 / GSE148071（三个 scRNA 队列）→ 存档 `/home/eto/luad_invasion/luad_v2_out_of_scope/`。

---

## 2. 实测数据特征（2026-09-12）

**GSE308103（snRNA，细胞核）** — 全队列分位数（1/5/25/50/75/95/99%）：

| 指标 | 1% | 5% | 25% | 50% | 75% | 95% | 99% |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| nCount | 500 | 561 | 902 | **1,516** | 2,970 | 11,296 | 31,807 |
| nFeature | 367 | 452 | 697 | **1,082** | 1,786 | 4,234 | 7,085 |
| pct_mt | 0 | 0 | 0.26 | **0.6** | 1.24 | 3.13 | 6.44 |

> ⚠️ **核数据**（低 nCount、低 mt%）→ **不得照搬整细胞 scRNA 阈值**（`nCount≥1000` 会砍掉约 30% 的核；
> `mt<10%` 几乎不筛）。QC 用**逐样本 MAD 离群** + `pct_mt < 5`。

**GSE307534（空间）** — 本地 **56/56** 张切片（2026-09-15 补齐），阶段覆盖（GEO 权威表，`IAC` = 原始 token `LUAD`）：

```
Normal ( 1): P4_Normal
AAH    (11): P1_AAH, P2_AAH, P4_AAH, P4_AAH-1, P6_AAH, P9_AAH, P11_AAH, P20_AAH, P22_AAH, P24_AAH, P25_AAH
AIS    (14): P3_AIS, P5_AIS, P8_AIS, P9_AIS, P12_AIS, P14_AIS, P16_AIS, P17_AIS, P19_AIS, P21_AIS, P21_AIS-1, P22_AIS, P23_AIS, P23_AIS-1
MIA    ( 4): P10_MIA, P13_MIA, P15_MIA, P18_MIA
IAC    (26): P1_LUAD, P2_LUAD, P3_LUAD, P4_LUAD, P5_LUAD, P6_LUAD, P7_LUAD, P7_LUAD-1, P8_LUAD, P9_LUAD, P10_LUAD, P11_LUAD, P12_LUAD, P13_LUAD, P14_LUAD, P15_LUAD, P16_LUAD, P17_LUAD, P18_LUAD, P19_LUAD, P20_LUAD, P21_LUAD, P22_LUAD, P23_LUAD, P24_LUAD, P25_LUAD
```

> ⚠️ 唯一缺陷：`GSM9226176` 的 tar **截断**（56,272,384 B / 应为 90,677,930 B；`gzip -t` 报 unexpected EOF），
> 缺 `spatial/scalefactors_json.json` 与 `spatial/tissue_positions.csv`。已实测重下载可得完整 87 MB tar（只有一个切片根 `P4_AAH2`）。

---

## 3. 方法学事实（已核实，决定选型）

| 结论 | 依据 |
| :--- | :--- |
| **Visium spot 非单细胞**（55 µm 混合多细胞）→ 必须**解卷积** | 平台定义 |
| 参考**模态需与切片匹配**（FFPE↔FFPE）→ 用 **GSE308103** 作参考 | RCTD / 基准文献 |
| **sc↔sn 是 system 效应**，不能用 `modality 当 batch` | scvi-tools（SysVI 专为此设）；Hrovatin 2025 |
| `scvi-tools 0.15.5` 的 scVI 参数默认值有效；**1.5.x 需 Py≥3.10（本机装不上）** | 实测 |
| **RCTD `doublet_mode='full'` 无 reject 类别，且 `constrain=F`**（权重非概率） | 源码实测 |
| `squidpy.nhood_enrichment` **置换标签、只返回 z、无 p 值** | 源码实测 |
| **SCMG 无「状态逆转」能力**；`generate_transition_cells` 只在参照流形**已有**类型间插值；参照流形无肿瘤态 | 源码 + 实测 |
| **观察性单细胞/空间不能建立因果**；合法杠杆 = **cis-MR + coloc** 或扰动实验 | 方法学共识 |
| 单细胞基础模型（scGPT/Geneformer…）**打不过** scVI/Harmony/PCA 基线 | Kedzierska 2025 *Genome Biol*；Ahlmann-Eltze 2025 *Nat Methods* |
| PLIP zero-shot 判 WHO 生长模式**未经验证**，且 spot 尺度**不足** | 文献 + 分辨率论证 |
| CMap 的 **Tau 是 0–1 重现性指标，不可能为负**；`Tau≤-90` 系误用 | Subramanian 2017 |

---

## 4. 环境事实（硬约束）

- **无 GPU**（无 `/dev/nvidia*`，`torch.cuda.is_available()=False`）→ 深度学习/共折叠类方法仅 CPU 或需外部节点；
- **Python 3.8.10、无 conda** → `scvi-tools 1.5`、`cellcharter` 装不上；
- **共享库 root 属主** → 装包进个人库（`~/.local/lib/python3.8/site-packages`、`~/R/.../4.2`）；
  曾因装 `cellcharter` 把 torch 降到 1.12 搞坏 scvi（已回滚 torch 2.4.1+cu118、pl 1.5.10.post0、torchmetrics 0.7.3）；
- **网络**：`http(s)_proxy` 指向死端口 → 出网须绕代理；`github.com` 被墙但 `api/codeload` 通；
  CRAN ✓、conda-forge ✓、PyPI 经清华镜像 ✓。

**已装可用工具**：
- R：`CopyKAT 1.2.5`、`coloc 5.2.3`、`ieugwasr 1.1.0`、`TwoSampleMR 0.7.9`、`scDblFinder 1.12.0`、`spacexr 2.2.1`、`DESeq2`、`scater/scuttle`
- Python：`harmonypy 0.0.10`、`SpaGCN 1.2.7`、`plip`（权重已缓存 1.2 GB）、`vina 1.2.7`、`fpocket`(CLI)
- 待装：BANKSY、PRECAST、P2Rank、gnina、PoseBusters、OpenMM（GROMACS 缺失）

---

## 5. 旧工程错误（勿再犯）

见 `/home/eto/luad_invasion/luad_v2_out_of_scope/` 与 v1 白皮书的历史记录：TD9 误标 AAH；`nLN` 记成 LNM；脑转/胸水并入 IAC；
写死数值与伪曲线；SCMG 逆转因子伪造；CMap 指标误用；GSE190811 当 LUAD LNM。

**原则**：一切标签/分期/数值须有权威来源；**未经核实不得入库**；凡无真实来源 → **拒绝产出**。

---

## 6. 相关文件

- 技术路线：[`WHITEPAPER.md`](WHITEPAPER.md)
- 检查点：`PLAN_AND_CHECKPOINTS.md`
- 参数出处：[`PARAMETERS_AND_SOURCES.md`](PARAMETERS_AND_SOURCES.md)
- 靶点 MR/coloc 手册：[`M7B_MR_COLOC_TARGET_ANCHORING.md`](M7B_MR_COLOC_TARGET_ANCHORING.md)
- 空间 P0 禁令：[`spatial_cohort_and_figure_prohibitions.md`](spatial_cohort_and_figure_prohibitions.md)
- 范围外存档：`/home/eto/luad_invasion/luad_v2_out_of_scope/`（含说明 README）
