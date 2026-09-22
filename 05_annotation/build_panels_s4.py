#!/usr/bin/env python3
"""生成 `classic_panels_s4.py` —— 内皮 + 髓系两个谱系的**新版**亚型 marker 面板。

为什么要做这个（2026-09-22）
────────────────────────────────────────────────────────────────────────────
GP8c 暴露了一个**结构性缺陷**：原面板取自 Travaglini 2020 **Table S1 的
「Canonical markers」列**，而该列在血管内皮上根本不是单细胞证据 ——

    Artery Cell   : GJA5, BMX        来源 Townsley 2012 +《The Lung》第 74 章（教科书）
                                     「Extant expression profiles」栏写的是 (bulk)/(cultured)
    Vein Cell     : ACKR1            同一批教科书来源；该栏**空**
    Capillary Cell: CA4              来源 Crapo 1982（生理学）；该栏**空**
    Lymphatic Cell: PROX1, PDPN      单细胞研究；该栏写 Yes

⇒ 作者自己对动脉只有 bulk/培养证据，对静脉和毛细血管**压根没有单细胞表达谱**。
  后果：Vein（1 槽）与 Capillary（1 槽）**可用基因 <2，结构上打不了分**；
  Artery/Lymphatic（各 2 槽）触及「薄」的判定线。内皮 4 型被压成 2 型。
  髓系同理：面板不含任何肥大细胞 marker、不含任何 DC marker
  ⇒ 全局簇 17（肥大细胞）与 26（DC）被错标为上皮 / B·浆。

修法（用户 2026-09-22 选定，**路径 (a)：同一篇论文的更细表**）
────────────────────────────────────────────────────────────────────────────
换用**同一篇论文的 Table S4**（per-cluster enriched markers，101 张 sheet）。
仍是 Travaglini 2020 **一篇文献**，举证的人只需核一张表 —— 溯源不跨论文。
S4 每张 sheet 第 0 行**自带类型名**，故不需要 S2/S3 的簇→类型对照表。

🔴 但 S4 是**按 avg_logFC 排序的完整差异表达表**（每张 200–2000 个基因），
   **不是** canonical 面板。直接拿去打分，动脉/静脉/毛细血管这类相关亚型的分数会
   高度相关，簇均值 argmax 会乱跳。故**必须先按「特异性」筛**（见下）。

口径（**用户 2026-09-22 签字**，法则 3.1/3.2）
────────────────────────────────────────────────────────────────────────────
① 只取 **10x 的 sheet**（sheet 名不含 "(SS2)"/"(SS)"）。
   理由：本数据是 10x Flex。SS2 = Smart-seq2（板式全长），技术不匹配。
② 逐型过筛：  avg_logFC ≥ 1.0   且   (pct_in_cluster − pct_out_cluster) ≥ 0.3
   （即「本簇里表达高、且明显比簇外特异」——判别式 marker 的标准判据）
③ 过筛名单**按 logFC 降序**保留（不在生成器里封顶），运行时取
   **前 20 个本矩阵里真有的基因**（PANEL_CAP=20）。
   ⇒ 「矩阵里没有 → 往下顺延补齐到 20」，使各型面板等长、打分力度齐。
④ 亚型粒度**照 S4 原粒度**：内皮 9 型、髓系 15 型（不再压回原来的粗粒度）。
⑤ S4 的 10x 表缺失的型（Neutrophil 只有 SS2；Eosinophil 无此簇）
   ⇒ **沿用 Table S1 行原样**，逐条登记于 FALLBACK_FROM_S1。
⑥ **不做别名修正**。S4 里有旧符号（CTGF/CYR61/SEPP1/KIAA0101/C10orf54 等），
   本项目**不把它们换成人同源/新符号**，也不找替身 —— 它们若不在矩阵，
   就由运行时的存在性检查如实报为 missing，并被 ③ 的顺延规则跳过。
   （用户 2026-09-21：「认人同源号允许；找替身填缺槽不允许」。
     换 CTGF→CCN2 属「认字」，允许但非必需；此处选择**不做**，
     把判断量降到最低 —— 丢一个 marker 是保守误差，不是伪造结果。）

🔴 本文件**由脚本生成，不要手改**。改口径请改生成器后重跑，并重新签字。
"""

import hashlib
import json
import os
from datetime import datetime

import pandas as pd

ROOT = "/home/eto/luad_v2"
XLSX = f"{ROOT}/data/external/travaglini2020/TableS4_cluster_enriched_markers.xlsx"

# ---------------------------------------------------------------- 口径常量
SELECTION_LOGFC = 1.0
SELECTION_SPEC = 0.3
PANEL_CAP = 20

# 本谱系 → S4 的 sheet（只列 10x；标签由文件第 0 行读出，不手抄）
ENDO_SHEETS = ["Cluster 16", "Cluster 17", "Cluster 18", "Cluster 19", "Cluster 20",
               "Cluster 21", "Cluster 22", "Cluster 23", "Cluster 24"]
MYE_SHEETS = ["Cluster 44", "Cluster 45", "Cluster 46", "Cluster 47", "Cluster 48",
              "Cluster 49", "Cluster 50", "Cluster 51", "Cluster 52", "Cluster 53",
              "Cluster 54", "Cluster 55", "Cluster 56", "Cluster 57", "Cluster 58"]

# S4 的 10x 表缺失的型 → 沿用 Table S1 行（逐条登记，不静默）
FALLBACK_FROM_S1 = {
    "髓系": {
        "Neutrophil": ["S100A8", "S100A9", "IFITM2", "FCGR3B"],
        "Eosinophil": ["SIGLEC8"],
    },
}
FALLBACK_REASON = {
    "Neutrophil": "Table S4 里中性粒**只有 SS2 表**（Cluster 43 (SS2)），无 10x 表。"
                  "严格走「只用 10x」会丢掉该型，而它在本数据里标注了 1,406 个核 ⇒ 沿用 S1 行。",
    "Eosinophil": "Table S4 **没有嗜酸簇**（该型稀少未被单独分出）⇒ 沿用 S1 行（其本来即 1 槽，不可打分）。",
}

# S4 论文自身的缺陷（原样登记）
SOURCE_TABLE_DEFECTS_S4 = {
    "Cluster 21 标签不一致": "10x 表第 0 行写「Capillary Intermediate 2」，同簇 SS2 表写"
                            "「Capillary Intermediate 1」。本项目只用 10x ⇒ 取「Capillary Intermediate 2」。"
                            "缺陷原样登记，不修（改了就成了我们的判断）。",
    "Cluster 50 标签不一致": "10x 表写「Myeloid Dendritic Type 1」，SS2 表写「Dendritic」。同上，取 10x 的标签。",
    "旧基因符号": "S4 混用旧符号（CTGF/CYR61/SEPP1/KIAA0101/C10orf54/FAM26F/NAPSB 等）。"
                  "本项目不做别名修正（见文件头 ⑥），不在矩阵者如实报 missing。",
}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_sheet(xl, sheet):
    """读一张 sheet → (类型名, 过筛后的基因名单[按 logFC 降序])。"""
    lab = str(pd.read_excel(xl, sheet_name=sheet, header=None, nrows=1).iloc[0, 0]).strip()
    d = pd.read_excel(xl, sheet_name=sheet, header=1)
    d.columns = [str(c).strip() for c in d.columns]
    for c in ("avg_logFC", "pct_in_cluster", "pct_out_cluster"):
        if c not in d.columns:
            raise SystemExit(f"🔴 {sheet}: 缺列 {c}；有 {list(d.columns)}")
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["Gene"] = d["Gene"].astype(str).str.strip()
    d = d.dropna(subset=["avg_logFC", "pct_in_cluster", "pct_out_cluster"])
    n_raw = len(d)
    d["spec"] = d["pct_in_cluster"] - d["pct_out_cluster"]
    d = d[(d["avg_logFC"] >= SELECTION_LOGFC) & (d["spec"] >= SELECTION_SPEC)]
    d = d.sort_values("avg_logFC", ascending=False)
    d = d[~d["Gene"].duplicated()]
    genes = [g for g in d["Gene"].tolist() if g and g.lower() != "nan"]
    return lab, genes, n_raw


def main():
    xl = pd.ExcelFile(XLSX)
    src_sha = sha256(XLSX)

    panels, slots, audit = {}, {}, {}
    for lineage, sheets in (("内皮", ENDO_SHEETS), ("髓系", MYE_SHEETS)):
        panels[lineage] = {}
        for s in sheets:
            if s not in xl.sheet_names:
                raise SystemExit(f"🔴 {s} 不在 S4 里")
            lab, genes, n_raw = read_sheet(xl, s)
            if lab in panels[lineage]:
                raise SystemExit(f"🔴 {lineage}: 类型名重复 {lab}")
            panels[lineage][lab] = genes
            slots[lab] = len(genes)
            audit[lab] = dict(sheet=s, n_raw_genes=n_raw, n_passed=len(genes))
        for t, gs in FALLBACK_FROM_S1.get(lineage, {}).items():
            panels[lineage][t] = list(gs)
            slots[t] = len(gs)
            audit[t] = dict(sheet="(S1 fallback)", n_raw_genes=len(gs), n_passed=len(gs))

    # 硬校验：过筛为空且未被登记的型必须显式出现（不静默）
    UNSCORABLE_FROM_S4 = sorted(t for t, v in slots.items() if v == 0)
    print(f"过筛后为空（无法打分）的型：{UNSCORABLE_FROM_S4}")
    print(f"逐型过筛数：{ {k: slots[k] for k in sorted(slots)} }")

    # 生成模块文本
    L = []
    A = L.append
    A('"""🔴 本文件由 `05_annotation/build_panels_s4.py` 程序化生成，**不要手改**。')
    A("")
    A("来源：Travaglini 2020 Nature 587:619，**Table S4**（per-cluster enriched markers）。")
    A("    仍是**同一篇论文**，举证只需核一张表；从 S1 的「Canonical markers」列换成 S4 的")
    A("    逐簇富集表，是因为 S1 那一列在血管内皮上不是单细胞证据（见生成器文件头）。")
    A("")
    A("口径（用户 2026-09-22 签字）：")
    A(f"    ① 只用 10x sheet（排除 SS2/Smart-seq2）")
    A(f"    ② 过筛 avg_logFC ≥ {SELECTION_LOGFC} 且 (pct_in − pct_out) ≥ {SELECTION_SPEC}")
    A(f"    ③ 按 logFC 降序，运行时取前 {PANEL_CAP} 个**本矩阵里真有的**（矩阵缺则顺延补齐）")
    A(f"    ④ 亚型粒度照 S4 原粒度")
    A("    ⑤ S4 无 10x 表的型沿用 Table S1 行（逐条登记于 FALLBACK_FROM_S1）")
    A("    ⑥ 不做别名修正（旧符号不在矩阵者如实报 missing）")
    A("")
    A("🔴 面板**只用于亚型（L2）注释**。六谱系的上位归属（L1）仍由 marker_panel.py 决定，")
    A("   不受本文件影响。")
    A('"""')
    A("")
    A(f"SOURCE = {XLSX!r}")
    A(f"SOURCE_SHA256 = {src_sha!r}")
    A("SOURCE_NOTE = ('Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the human "
      "lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "
      "doi:10.1038/s41586-020-2922-4；PMC7704697【Table S4】')")
    A("")
    A(f"SELECTION_LOGFC = {SELECTION_LOGFC}")
    A(f"SELECTION_SPEC = {SELECTION_SPEC}")
    A(f"PANEL_CAP = {PANEL_CAP}")
    A("")
    A("# S4 的 10x 表缺失、沿用 Table S1 的型（逐条登记，不静默）")
    A(f"FALLBACK_FROM_S1 = {FALLBACK_FROM_S1!r}")
    A(f"FALLBACK_REASON = {FALLBACK_REASON!r}")
    A("")
    A(f"SOURCE_TABLE_DEFECTS = {SOURCE_TABLE_DEFECTS_S4!r}")
    A("")
    A(f"UNANNOTATABLE = {{}}   # S4 里每型都有 marker；原 S1 的 Bronchial Vessel 在 S4 中已有基因")
    A(f"UNSCORABLE_FROM_S4 = {UNSCORABLE_FROM_S4!r}")
    A("")
    A("RARE = ['Platelet/Megakaryocyte', 'Eosinophil', 'Bronchial Vessel 1', 'Bronchial Vessel 2',")
    A("        'Capillary Intermediate 1', 'Capillary Intermediate 2']")
    A("")
    A("CONTAM_CHECK = ['PTPRC', 'CD3D', 'COL1A1', 'PECAM1', 'CD68', 'EPCAM']")
    A("")
    A("# 过筛后的**完整**名单（按 logFC 降序）。运行时的封顶与矩阵存在性检查在 build_panel。")
    A("PANELS = {")
    for lineage in ("内皮", "髓系"):
        A(f"    {lineage!r}: {{")
        for t, gs in panels[lineage].items():
            A(f"        {t!r}: {gs!r},")
        A("    },")
    A("}")
    A("")
    A("# 过筛数（= 源表对该型给出的可区分基因数；不随本矩阵变化）")
    A("N_SLOTS = {")
    for lineage in ("内皮", "髓系"):
        A(f"    {lineage!r}: {{")
        for t in panels[lineage]:
            A(f"        {t!r}: {slots[t]},")
        A("    },")
    A("}")
    A("")
    A("AUDIT = {")
    for t in sorted(audit):
        A(f"    {t!r}: {audit[t]!r},")
    A("}")
    A("")
    A('''
def build_panel(lineage, available_genes=None):
    """返回某一谱系的面板 {亚型: [基因...]}。

    available_genes : 本数据矩阵的基因集合。给了则：
        ① 每型按 logFC 降序**顺延**取前 PANEL_CAP 个真在矩阵里的基因
           （矩阵里没有的跳过、继续往下取，使各型面板等长）；
        ② 返回 (面板, 缺失清单) —— 缺失清单记的是**原名单里所有**不在矩阵的基因。
    不给则原样返回（不封顶）。
    """
    if lineage not in PANELS:
        raise KeyError(f"未登记的谱系：{lineage}；有 {list(PANELS)}")
    if available_genes is None:
        return {k: list(v) for k, v in PANELS[lineage].items()}
    have = set(available_genes)
    used, missing = {}, {}
    for t, gs in PANELS[lineage].items():
        miss = [g for g in gs if g not in have]
        used[t] = [g for g in gs if g in have][:PANEL_CAP]
        if miss:
            missing[t] = miss
    return used, missing


def assert_provenance():
    """硬校验（与 classic_panels.assert_provenance 同精神）。"""
    problems = []
    for lineage, d in PANELS.items():
        if not d:
            problems.append(f"{lineage}: 面板为空")
        for t, gs in d.items():
            if len(gs) == 0 and t not in UNSCORABLE_FROM_S4 and t not in FALLBACK_FROM_S1.get(lineage, {}):
                problems.append(f"{lineage}/{t}: 名单为空且未登记")
            if len(gs) != len(set(gs)):
                problems.append(f"{lineage}/{t}: 名单内有重复基因")
            bad = [g for g in gs if not g or g.lower() == "nan" or " " in g]
            if bad:
                problems.append(f"{lineage}/{t}: 非法基因名 {bad}")
    for lineage, d in FALLBACK_FROM_S1.items():
        for t in d:
            if t not in PANELS[lineage]:
                problems.append(f"fallback {lineage}/{t} 未进入面板")
    if problems:
        raise AssertionError("面板溯源校验失败：\\n  " + "\\n  ".join(problems))
    return True
''')
    out = f"{ROOT}/05_annotation/classic_panels_s4.py"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\n写出 {out}")
    print(f"  源表 sha256 = {src_sha}")
    print(f"  内皮 {len(panels['内皮'])} 型 / 髓系 {len(panels['髓系'])} 型")
    print(f"  生成文件 sha256 = {sha256(out)}")
    register_caliber(panels, slots, src_sha, out, audit, sign=SIGN)


def register_caliber(panels, slots, src_sha, out_module, audit, sign=""):
    """写出 `results/05_annotation/panel_caliber_s4.json` —— **口径登记**（法则 3.1）。

    `sign` 为空 ⇒ `signed_by` 留空 ⇒ 03_subcluster_annotation.py 的闸门拒绝注释。
    签署后才可跑（法则 3.1/3.2：新口径先冻结签字、再计算）。
    """
    import anndata
    h5 = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
    ad = anndata.read_h5ad(h5, backed="r")
    have = set(ad.var_names)
    n_var = ad.n_vars
    del ad

    inv = {}
    for lineage in ("内皮", "髓系"):
        inv[lineage] = {}
        for t, gs in panels[lineage].items():
            used = [g for g in gs if g in have][:PANEL_CAP]
            inv[lineage][t] = dict(
                sheet=audit[t]["sheet"], n_source_genes=audit[t]["n_raw_genes"],
                n_passed=slots[t],
                n_in_matrix=len([g for g in gs if g in have]),
                n_used=len(used), scorable=len(used) >= 2, genes=used)

    j = dict(
        panel="s4",
        purpose="修 GP8c 暴露的结构性缺陷：Table S1 的 Canonical markers 列在血管内皮上"
                "不是单细胞证据（动脉 bulk/培养、静脉与毛细血管无单细胞谱），导致内皮 4 型被压成 2 型；"
                "髓系面板无肥大细胞/DC marker，导致全局簇 17/26 被错标。",
        scope="🔴 只改**内皮、髓系两个谱系**的亚型（L2）面板。其余四谱系仍用 s1（Table S1）。"
              "六谱系上位归属（L1，marker_panel.py）**不受影响**；污染体检（06）与谱系接管（10）"
              "用的是 L1 面板，故其判据不变。",
        source=dict(file=os.path.relpath(XLSX, ROOT), sha256=src_sha,
                    citation="Travaglini KJ, Nabhan AN, Penland L, et al. Nature 2020;587(7835):619-625. "
                             "doi:10.1038/s41586-020-2922-4；PMC7704697【Table S4：per-cluster enriched markers】",
                    same_paper_as_s1=True,
                    why_switch="仍是同一篇论文，举证只需核一张表；从 S1 的 Canonical 摘要列"
                               "换成 S4 的逐簇富集表，因为前者在内皮上非单细胞证据。"),
        panel_module=dict(file=os.path.relpath(out_module, ROOT), sha256=sha256(out_module)),
        generator=dict(file="05_annotation/build_panels_s4.py",
                       sha256=sha256(f"{ROOT}/05_annotation/build_panels_s4.py")),
        rule=dict(
            only_10x=True,
            only_10x_why="本数据是 10x Flex；SS2 = Smart-seq2（板式全长），技术不匹配。",
            filter=dict(avg_logFC_min=SELECTION_LOGFC, specificity_min=SELECTION_SPEC,
                        specificity_def="pct_in_cluster − pct_out_cluster",
                        why="S4 是按 logFC 排序的完整差异表（每张 200–2000 基因），不是 canonical 面板。"
                            "不筛特异性则相关亚型（动脉/静脉/毛细血管）分数高度相关、argmax 乱跳。"),
            cap=dict(n=PANEL_CAP, order="avg_logFC 降序",
                     fill="矩阵里没有的基因**跳过、继续往下取**，直到凑足 20 个真有的 ⇒ 各型面板等长"),
            granularity="照 S4 原粒度（内皮 9 型、髓系 S4 部分的 15 型）",
            fallback_from_s1=FALLBACK_FROM_S1,
            fallback_reason=FALLBACK_REASON,
            no_alias_fix=dict(decision=True,
                              why="S4 混用旧符号（CTGF/CYR61/SEPP1/KIAA0101/C10orf54 等）。"
                                  "换 CTGF→CCN2 属「认字」允许，但此处**不换也不找替身**，"
                                  "把判断量降到最低——丢一个 marker 是保守误差，不是伪造结果。"
                                  "缺失者由运行时存在性检查如实报 missing。"),
        ),
        inventory=inv,
        known_limitations=[
            "Capillary Intermediate 2：S4 过筛后**0 个基因** ⇒ 该型**无法打分**。9 型里能用 8 型。"
            "性质与 S1 不同：S1 是表写太省，此处是**源数据本身没有可区分的基因**（最高 logFC 不够）。",
            "HLA-DP/DQ/DR 整个家族不在本矩阵（10x Flex 探针板缺口）⇒ 6 个 DC 亚型最经典的 marker 拿不到，"
            "只能退用次强的。MDT2 因此只剩 13 个基因（过筛 25、矩阵缺 12）。",
            "Neutrophil 在 S4 只有 SS2 表、Eosinophil 在 S4 无此簇 ⇒ 二者沿用 S1 行（1,406 核的中性粒不丢）。"
            "⇒ 髓系内部有两个来源，本登记已逐条写明。",
            "内皮面板里混入免疫基因（Capillary 的 IL7R、Capillary Intermediate 1 的 IL1RL1）"
            "⇒ 提示 Travaglini 自己的簇里也有环境 RNA/双体。须在结果报告里如实标注。",
            "Capillary（Cluster 19）在其源数据里也弱（最高 logFC 1.26、过筛 4）⇒ 即使换表，该型仍最弱，如实上报。",
            "🔴 名与实有张力的胞型（**不改名**，但报告里不得按字面讲）：S4 名为「TREM2+ Dendritic」"
            "（Cluster 54）者，拿到的 marker 是 TREM2/APOE/C1QB/GPNMB/CHIT1/CHI3L1/FN1/MMP9 —— "
            "这是**脂质相关巨噬细胞（LAM）**的经典程序，不是 DC 程序；「IGSF21+ Dendritic」"
            "（Cluster 52）为 FOLR2/C1QC/CD163/MRC1/STAB1/F13A1，同样偏组织巨噬；"
            "「EREG+ Dendritic」（Cluster 53）为 IL1B/CXCL8/EREG/CD14/CD163，偏炎性单核。"
            "⇒ 三型必须报为「S4 标为 DC、marker 谱实为单核/巨噬程序」，"
            "不得在报告里直接当 DC 讲。",
            "🔴 我把 S4 的 Cluster 44–58 整段当作「髓系」纳入，其中 Platelet/Megakaryocyte"
            "（Cluster 46）严格说不属髓系。依据是它在 S4 簇编号里与其余髓系**连续**，"
            "推测 Travaglini 把它放在髓系超簇下。**这是我的判断，未经用户复核**；"
            "副作用：新面板首次给了血小板 marker（PPBP/PF4/GP9/TUBB1），"
            "若 myeloidA 里真有血小板群，现在才看得见（旧面板看不见）。",
        ],
        review_required=[
            "核对本登记的口径与用户 2026-09-22 的裁定是否逐条一致（5 项：10x-only / 过筛 / 封顶 / 粒度 / 回退）。",
            "签字前可先看 inventory 里的逐型基因名单是否生物学合理。",
        ],
        signed_by=sign,
        signed_at=(datetime.now().strftime("%Y-%m-%dT%H:%M:%S") if sign else ""),
        # 口径的**内容**（rule 那一整块 + inventory）是用户 2026-09-22 在会话内经
        # AskUserQuestion 亲自选的 5 个选项，不是我拟的。这个字段只说清「谁把已选口径
        # 冻结成文件、谁来追认」，避免签字被读成「Claude 自己定了口径又自己批了」。
        sign_scope=dict(
            caliber_chosen_by="用户（2026-09-22 会话内，AskUserQuestion 五项：10x-only / "
                              "过筛 logFC≥1.0 且特异度≥0.3 / 每型封顶 20 且顺延补齐 / "
                              "照 S4 原粒度 / S4 无表者沿用 S1）",
            frozen_by=sign or "（未签）",
            countersign_required=True,
            countersign_note="🔴 若 frozen_by 不是用户本人，本口径下的结果一律标"
                             "**探索性（provisional）**，直到用户追认为止。追认方式："
                             "把 countersigned_by/countersigned_at 填上。"),
        sign_note="🔴 签字方式：`python3 05_annotation/build_panels_s4.py --sign \"<署名>\"`。"
                  "未签字前 03_subcluster_annotation.py --panel s4 会 SystemExit。"
                  "🔴 签的是**这一份具体面板**（panel_module_sha256）；签字后改面板即失效。",
    )
    p = f"{ROOT}/results/05_annotation/panel_caliber_s4.json"
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(j, fh, ensure_ascii=False, indent=2)
    print(f"\n写出 {os.path.relpath(p, ROOT)}"
          + ("（已签字 ⇒ 闸门放行）" if sign else "（signed_by 为空 ⇒ 闸门会拦）"))
    print(f"  矩阵基因数 {n_var}")
    for lineage in ("内皮", "髓系"):
        n_ok = sum(1 for v in inv[lineage].values() if v["scorable"])
        print(f"  {lineage}: {n_ok}/{len(inv[lineage])} 型可打分")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--sign", default="",
                    help="把 `signed_by` 写成这个署名并盖上时刻。留空 ⇒ 只登记不签字，"
                         "闸门会拒绝注释。🔴 签的是**这一份具体面板**的哈希，改面板即失效。")
    a = ap.parse_args()
    SIGN = a.sign
    main()
    if not SIGN:
        print("\n🔴 signed_by 为空 —— 03_subcluster_annotation.py --panel s4 仍会拒绝注释。"
              "\n   确认口径无误后再跑：python3 05_annotation/build_panels_s4.py --sign \"<署名>\"")
