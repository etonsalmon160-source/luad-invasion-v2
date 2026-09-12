"""
cohort_registry.py  —  权威队列登记表（single source of truth）

用途：把"用哪些数据集、每个样本属哪个阶段、模态是什么"**集中登记**，供下游一致引用，
      杜绝散落各处的启发式/静默默认。

权威口径来源：GEO 逐样本 histology / Sample_Origin（见 README 与白皮书 §2.1）。
"""

# -----------------------------------------------------------------------------
# 1. 单细胞 scRNA 队列（三个，均【无 AAH】）
# -----------------------------------------------------------------------------
SC_RNA = {
    "GSE131907": {
        "modality": "scRNA",
        "note": "Kim et al.; 分期按 Sample_Origin（非字符串启发式）",
        "stage_by_origin": {          # Sample_Origin -> 阶段
            "nLung": "Normal", "tLung": "IAC", "tL/B": "IAC",
            "mLN": "LNM",             # 真转移淋巴结
            "nLN": "Normal_LN",       # ⚠️ 正常淋巴结（**非** LNM）
            "mBrain": "Brain_Met",    # ⚠️ 单列，勿并入 IAC
            "PE": "Pleural_Effusion", # ⚠️ 单列，勿并入 IAC
        },
    },
    "GSE189357": {
        "modality": "scRNA",
        "note": "9 样本 = 3 AIS + 3 MIA + 3 IAC；**无 AAH**",
        "stage_by_TD": {              # 逐样本 histology（GEO 核实）
            "TD1": "IAC", "TD2": "IAC", "TD3": "MIA", "TD4": "MIA",
            "TD5": "AIS", "TD6": "MIA", "TD7": "AIS", "TD8": "AIS", "TD9": "IAC",
        },
        "corrections": "旧版误标：TD9=AAH（实 IAC）、TD4=IAC（实 MIA）",
    },
    "GSE148071": {
        "modality": "scRNA",
        "note": "系列标题为 **Advanced NSCLC**（42 样本），非'早期 LUAD'",
        "stage": "Adv_NSCLC",
    },
}

# -----------------------------------------------------------------------------
# 2. AAH 单细胞来源（已决策 + 预留接口）
# -----------------------------------------------------------------------------
AAH_SOURCE = {
    "active": {
        "dataset": "GSE308103",
        "modality": "snRNA (FFPE, 10x Fixed RNA)",
        "note": "与空间 GSE307534 同患者配对；经 SCMG zero-shot 跨平台并入",
        "stages": ["Normal", "AAH", "AIS", "MIA", "LUAD"],
        "caveat": "snRNA≠scRNA；模态与队列共线，须过跨模态五判据（见白皮书 §3.1）",
    },
    "reserved": {
        "dataset": "HRA001130",
        "modality": "scRNA (whole-cell)",
        "note": "含 AAH/AIS/MIA/IA + 配对正常；**GSA-Human 受控，待申请**",
        "interface": "00_ingest/hra001130_interface.py",
        "doi": "10.1038/s41467-021-26770-2",
    },
}

# -----------------------------------------------------------------------------
# 3. 空间队列（权威切片口径；P0 禁令见 .agents/rules/...）
# -----------------------------------------------------------------------------
SPATIAL = {
    "core_6stage": {
        "Normal": "P4_Normal", "AAH": "P1_AAH", "AIS": "P3_AIS",
        "MIA": "P10_MIA", "IAC": "P3_LUAD", "LNM": "PT_3_LNM",
    },
    "forbidden": ["P4_AAH1 as MIA", "P4_AAH as MIA", "P3_LUAD as LNM"],
    "source": "GSE307534 / GSE189487 / GSE190811",
}

# -----------------------------------------------------------------------------
# 阶段词表（统一口径；AAH 由 AAH_SOURCE 提供）
# -----------------------------------------------------------------------------
STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC", "LNM"]
# 额外/单列（**不并入** IAC/LNM）：
NON_CORE = ["Normal_LN", "Brain_Met", "Pleural_Effusion", "Adv_NSCLC"]


def summary():
    print("=== 单细胞 scRNA（无 AAH）===")
    for k, v in SC_RNA.items():
        print(f"  {k}: {v['modality']} — {v['note']}")
    print("=== AAH 单细胞来源 ===")
    print(f"  active  : {AAH_SOURCE['active']['dataset']} ({AAH_SOURCE['active']['modality']})")
    print(f"  reserved: {AAH_SOURCE['reserved']['dataset']} ({AAH_SOURCE['reserved']['modality']})")
    print("=== 阶段词表 ===", STAGES, "| 单列:", NON_CORE)


if __name__ == "__main__":
    summary()
