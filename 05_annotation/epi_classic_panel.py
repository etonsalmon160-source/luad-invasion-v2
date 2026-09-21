"""GP8a 上皮亚聚类 —— **经典 marker** 亚型面板（逐基因可溯源到一次文献）。

生成方式：本文件由 `05_annotation/build_epi_classic_panel.py` 从原始补充材料 .xlsx
**程序化生成**，非手工转录 —— 保证基因名逐字一致。

口径（用户 2026-09-17 指令）
---------------------------
「就不按什么kac去做，就按照经典marker…按照权威经典marker做重聚类和celltype那个算法比比结果」
⇒ 主干 = **Travaglini 2020 Nature 587:619 的 Table S1「Canonical markers」**；
  源论文的 meta-program / KAC 体系**不再作本面板**（存档于 docs/PARAMETERS_AND_SOURCES.md §M3-A.5b）。

为什么 Table S1 是最合适的"权威经典"
--------------------------------------
它是**人工整编的教科书式规范表**（按 Epithelium / Endothelium / Stroma / PNS / Immune 分块），
**不是**聚类富集输出 ⇒ 不携带本数据的任何信息，作面板无循环性。

⚠️ 三条必须随结果报告的限制
----------------------------
1. **该表面向健康人肺，不是肿瘤**。LUAD 中相当一部分上皮是恶性的，用正常 marker 给恶性细胞
   打标签本身有偏。
2. **本面板无法区分恶性/正常** —— LUAD 常保留 NKX2-1/SFTPC 等肺泡身份。⇒ 本面板只出"正常上皮
   亚型"标签，**不设 tumor 亚型**；恶性身份一律由 GP2 CNV 判（铁律 R2），两者**交叉报告**。
3. **稀有型**（Serous/Ionocyte/Neuroendocrine/Tuft）正常肺合计 <0.3%，照算照报，但胜出簇过小时
   须标 `rare/likely-spurious`。

对原表的偏离（2026-09-21 用户裁定后）
------------------------------------
裁定把两种"偏离原表"分开，处置不同：

**① 认人同源 —— 保留（1 类偏离）**
`ORTHOLOG_FIX`：Table S1 Club 行的 `CYP2F2` 是**小鼠**基因号（该表其余格均为人的，此格系
物种笔误）；换用人同源 `CYP2F1`。**这不叫"补"** —— `CYP2F1` 就是原表点名那个基因的人版本，
且**同篇论文自己的 Table S4**（用他们自己的人类数据算出的簇富集表）相关两行写的正是 `CYP2F1`。
⇒ 依据充足，保留。

**② 找替身填缺槽 —— 删除（0 类）**
原 `TABLE_S4_FILL`（`DAPL1` 缺则取同篇 Table S4 的 `KRT15`/`KRT17` 顶替）**已按裁定删除**。
⚠️ 它此前**从未被任何代码读取**（生成器写了、产线没用，`grep` 核实）⇒ 删除它
**不改变任何已算结果**，只是清死代码并冻结该裁定。

🔴 后果：两个型被削弱（**须随结果报告**）
----------------------------------------
| 型 | Table S1 槽 | 实际可用 | 丢的是 | 为什么 |
|---|---|---|---|---|
| Basal | 4 | **3** | `DAPL1` | 本矩阵没有；按裁定不取 `KRT15`/`KRT17` 顶替 |
| Serous | 3 | **2** | `PRR4` | 本矩阵没有；按裁定不留替身 |

⇒ **Basal 判别力下降**（用户做此裁定时已被告知）。产线对 ≤2 基因的型自动打
  `thin_panel`；`KNOWN_ABSENT_IN_MATRIX` 供报告逐条引用。Club 不受影响（可用 3 个基因）。
"""

# ---------------------------------------------------------------- 一次文献
REFERENCES = {
    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "
                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "
                      "doi:10.1038/s41586-020-2922-4；PMC7704697 "
                      "【本面板主干 Table S1 Canonical markers；Club 行另有物种核对见 Table S4】",
    "Habermann2020": "Habermann AC, Gutierrez AJ, Bui LT, et al. Single-cell RNA sequencing reveals "
                     "profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary "
                     "fibrosis. Sci Adv 2020;6(28):eaba1972. doi:10.1126/sciadv.aba1972 "
                     "【仅用于对 Table S1 存疑条目做第二意见核对，不入本面板】",
}

# -------- 小鼠基因号 → 人同源号。2026-09-21 用户裁定：**此类保留**（「只认人同源」）
# 裁定把两种"偏离原表"分开：认人同源号（= 认字）允许；找替身填缺槽（= 找替身）不允许。
# 依据：Table S1 其余格均为人类基因号，唯 Club 行写 CYP2F2（小鼠），系物种笔误；
#      同篇论文 Table S4（他们自己的人类数据）该行写的正是 CYP2F1。
ORTHOLOG_FIX = {
    "CYP2F2": "CYP2F1",
}

# Table S1 上皮块中**本矩阵肯定没有**的槽（供报告引用；实际缺失由矩阵过滤现算）
KNOWN_ABSENT_IN_MATRIX = {
    "DAPL1": "本矩阵无此基因；按 2026-09-21 裁定不从 Table S4 取 KRT15/KRT17 顶替（找替身）",
    "PRR4": "本矩阵无此基因；按 2026-09-21 裁定不留替身",
}

# ------------------------------------------------- Table S1 经典 marker（主干，逐字）
TABLE_S1_MARKERS = {
    "Club": ["CYP2F1", "SCGB3A2", "CCKAR"],
    "Ciliated": ["FOXJ1", "TUBB1", "TP73", "CCDC78"],
    "Basal": ["KRT5", "KRT14", "TP63", "DAPL1"],
    "Serous": ["PRR4", "LPO", "LTF"],
    "Ionocyte": ["CFTR", "FOXI1", "ASCL3"],
    "Neuroendocrine": ["CALCA", "CHGA", "ASCL1"],
    "Tuft": ["DCLK1", "ASCL2"],
    "AT1": ["AGER", "PDPN", "CLIC5"],
    "AT2": ["SFTPB", "SFTPC", "SFTPD", "MUC1", "ETV5"],
    "Goblet/Mucous": ["MUC5AC", "MUC5B", "SPDEF"],
}

# 面板 = Table S1 主干（仅 1 类偏离：Club 行小鼠号→人同源；无 Table S4 补齐）

# 合并说明（Goblet 与 Mucous 在 Table S1 里共用 MUC5B ⇒ 无法分开）
MERGED = {
    "Goblet" + "Mucous": "Goblet/Mucous",
}

RARE = ['Serous', 'Ionocyte', 'Neuroendocrine', 'Tuft']
CONTAM_CHECK = ['PTPRC', 'CD3D', 'COL1A1', 'PECAM1', 'CD68']

# 逐型基因槽数（Table S1 原文，未增补前）
N_SLOTS = {
    "Club": 3,
    "Ciliated": 4,
    "Basal": 4,
    "Serous": 3,
    "Ionocyte": 3,
    "Neuroendocrine": 3,
    "Tuft": 2,
    "AT1": 3,
    "AT2": 5,
    "Goblet/Mucous": 3,
}


def build_panel(available_genes=None, table_s1=None):
    """返回 {亚型: [基因...]}。

    available_genes : 本数据矩阵的基因集合。给了则**只保留在矩阵里的**，
                      并返回 (面板, 缺失清单)；不给则原样返回 Table S1 主干。
    """
    panel = {k: list(v) if table_s1 is None else list(table_s1[k])
             for k, v in TABLE_S1_MARKERS.items()}
    if available_genes is None:
        return panel
    have = set(available_genes)
    used, missing = {}, {}
    for k, gs in panel.items():
        used[k] = [g for g in gs if g in have]
        miss = [g for g in gs if g not in have]
        if miss:
            missing[k] = miss
    return used, missing


if __name__ == "__main__":
    tot = sum(N_SLOTS.values())
    print(f"经典 marker 面板：{len(N_SLOTS)} 个上皮亚型，基因槽 {tot} 个")
    for k, n in N_SLOTS.items():
        tag = "  [稀有]" if k in RARE else ""
        print(f"  {k:16s} {n} 个{tag}")
