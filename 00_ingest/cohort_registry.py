"""
cohort_registry.py  —  权威队列登记表（single source of truth）

范围（2026-09-12 收窄）：**仅两个配对数据集**
  * GSE308103 —— snRNA（细胞核），含 AAH；单细胞/解卷积参考
  * GSE307534 —— Visium CytAssist FFPE 空间；与 GSE308103 **同一研究、同一批患者**

配对价值：同一患者、同一病灶（相邻切片）+ **模态匹配（FFPE↔FFPE）**
  → 让空间解卷积有可信参考，并支撑跨模态一致性判据。

权威口径来源：GEO 逐样本（见 docs/PROJECT_SUMMARY.md §数据底座）。
铁律：**分期无静默默认**（未知 token 一律 raise）。
"""

# -----------------------------------------------------------------------------
# 1. 配对患者（经 GSE307534 与 GSE308103 样本名交叉核实）
# -----------------------------------------------------------------------------
PAIRED_PATIENTS = ["P3", "P4", "P10", "P13", "P15", "P18", "P21", "P22", "P25"]

# -----------------------------------------------------------------------------
# 2. 两个配对数据集
# -----------------------------------------------------------------------------
COHORTS = {
    "GSE308103": {
        "modality": "snRNA",                 # ⚠️ 细胞核：QC 阈值不同于整细胞
        "platform": "10x (FFPE, fixed RNA)",
        "role": "单细胞参考（含 AAH）；空间解卷积的参考签名来源",
        "n_samples_geo": 75,
        "sample_re": r"^(?P<gsm>GSM\d+)_(?P<patient>P\d+)_(?P<stage_token>[A-Za-z0-9]+)\.raw_counts",
        "note": "样本名 <GSM>_<Pxx>_<Stage>.raw_counts.mtx.txt.gz；矩阵为 **基因×细胞** 稠密文本",
    },
    "GSE307534": {
        "modality": "spatial",               # Visium spot（55 µm，非单细胞）
        "platform": "10x Visium CytAssist FFPE",
        "role": "空间图谱（原位坐标）；解卷积对象",
        "n_samples_geo": 56,
        "sample_re": r"^(?P<gsm>GSM\d+)_(?P<patient>P\d+)_(?P<stage_token>[A-Za-z0-9\-]+)$",
        "note": "样本名 <GSM>_<Pxx>_<Stage>；每个含 filtered_feature_bc_matrix + spatial/",
    },
}

# -----------------------------------------------------------------------------
# 3. 分期 token → 规范阶段（严格；无默认回退）
#    LUAD = 浸润性腺癌 = IAC；AAH1/Normal1/AIS1 等为同一患者的第二切片
# -----------------------------------------------------------------------------
STAGE_MAP = {
    "Normal": "Normal", "Normal1": "Normal",
    "AAH": "AAH", "AAH1": "AAH", "AAH-1": "AAH",   # P4 的第二个 AAH 病灶
    "AIS": "AIS", "AIS1": "AIS",
    "MIA": "MIA",
    "LUAD": "IAC", "LUAD1": "IAC",       # 数据中出现的变体
    "IAC": "IAC",
}

STAGES = ["Normal", "AAH", "AIS", "MIA", "IAC"]     # 本项目配对范围（LNM 空转暂缺）
LNM_STATUS = "空转 LNM 暂缺（公开库无合法 LUAD 淋巴结转移 Visium；待真实数据）"

# 禁用法（P0，见 docs/spatial_cohort_and_figure_prohibitions.md）
FORBIDDEN = [
    "GSE190811 作为 LNM（经 GEO 核实为**乳腺癌**）",
    "P4_AAH / P4_AAH-1 充当 MIA",
    "任何切片伪标为 LNM",
]


def resolve_stage(token: str) -> str:
    """严格分期解析：未知 token 一律 raise，**绝不静默回退**。"""
    if token in STAGE_MAP:
        return STAGE_MAP[token]
    raise KeyError(f"未知分期 token {token!r}；权威映射只含 {sorted(STAGE_MAP)} —— 拒绝静默默认。")


def summary():
    print("=== 配对数据集（本项目范围）===")
    for k, v in COHORTS.items():
        print(f"  {k}: {v['modality']:8s} | {v['platform']} | GEO 样本 {v['n_samples_geo']}")
        print(f"      role: {v['role']}")
    print(f"\n配对患者（{len(PAIRED_PATIENTS)}）: {', '.join(PAIRED_PATIENTS)}")
    print(f"阶段词表: {STAGES}")
    print(f"LNM: {LNM_STATUS}")
    print("\n禁用法（P0）:")
    for f in FORBIDDEN:
        print(f"  ⛔ {f}")


if __name__ == "__main__":
    summary()
