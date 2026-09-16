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
import os

# -----------------------------------------------------------------------------
# 1. 配对患者（= 两个 GEO 权威表患者号的交集）
#    ⚠️ 2026-09-15 修正：早期空转仅下载 19/56 张时交集为 9 例；56 张齐后实为 23 例
#    （P3–P25）。**旧值 9 是下载不全的产物**，已作废。
# -----------------------------------------------------------------------------
PAIRED_PATIENTS_MIN = 23          # 交集规模下限（不满足即视为数据不全，非静默通过）

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
    "AIS": "AIS", "AIS1": "AIS", "AIS-1": "AIS",   # P21/P23 的第二个 AIS 病灶
    "MIA": "MIA",
    "LUAD": "IAC", "LUAD1": "IAC", "LUAD-1": "IAC",  # P7 的第二个 LUAD 病灶
    "IAC": "IAC",
}

# ⚠️ "第二病灶"的命名约定**两个数据集不一样**，这一点必须是显式的。
#    （2026-09-16 修正：此处早先只写 `-1`，并据此写了一个"按后缀猜序号"的函数
#      `lesion_ordinal()`。那个函数对 GSE307534 正确、对 **GSE308103 错误**，
#      而它从未被任何代码调用 —— 详见 results/03_cnv/GP1_report.md §6.1。）
#    * GSE307534（空转）  ：`AAH-1` / `LUAD-1` / `AIS-1`            —— **带横线**
#    * GSE308103（snRNA）：`AAH1` / `Normal1` / `LUAD1` / `AIS1`    —— **无横线**
#    两者都表示"同一患者的第二个独立病灶"（GEO 标题作 "second ... of patient N"），
#    与首个病灶分期相同但**解剖独立**，**不可合并**。
#    故：**不得由 token 猜序号**。唯一权威来源是 GEO 表，读取走 resolve_lesion_ordinal()。
#    `STAGE_MAP` 已同时登记两种写法（分期可解析）；序号则一律查表。

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


def resolve_lesion_ordinal(sample_id: str, dataset: str = None) -> int:
    """病灶序号的**唯一权威来源**：查 GEO 表。

    1 = 该患者该分期的首个病灶；2 = 第二个独立病灶。
    `sample_id` 形如 `P4_AAH1`（= GEO 表的 `f"{patient_id}_{token}"`）。
    查不到即 raise —— **绝不按 token 后缀猜**（两数据集约定不同，猜必错其一）。

    ⚠️ 本函数替代了早先的 `lesion_ordinal(token)`；后者按 `-1` 后缀猜，
    对 GSE308103 的 `AAH1`/`AIS1`/… 一律返回 1（真值为 2）。
    """
    datasets = [dataset] if dataset is not None else list(COHORTS)
    for ds in datasets:
        for g in load_geo(ds).values():
            if f"{g['patient_id']}_{g['token']}" == sample_id:
                return g["lesion_ordinal"]
    raise KeyError(f"GEO 权威表中无此样本 {sample_id!r}（查过 {datasets}）"
                   f" —— 拒绝由 token 猜病灶序号。")


# -----------------------------------------------------------------------------
# 4b. GEO 权威样本表（GI 逐样本标题，provenance=00_ingest/geo_metadata/）
#     用途：与**本地目录名**推导的 (patient, stage) 交叉核验；不一致即 FAIL。
# -----------------------------------------------------------------------------
GEO_META_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geo_metadata")

# 空转切片解压后**必须存在**的文件（相对切片根目录）
SPATIAL_REQUIRED = (
    "filtered_feature_bc_matrix/matrix.mtx.gz",
    "filtered_feature_bc_matrix/barcodes.tsv.gz",
    "filtered_feature_bc_matrix/features.tsv.gz",
    "spatial/scalefactors_json.json",
)
# 必需但**文件名随 Space Ranger 版本而异**的项：每组至少命中一个
SPATIAL_REQUIRED_ANY = (
    ("spatial/tissue_positions.csv", "spatial/tissue_positions_list.csv"),
)


def load_geo(dataset: str) -> dict:
    """GEO 权威 GSM → dict(patient_id, stage, lesion_ordinal, token, title)。

    stage 为 **GEO 原始分期词**（LUAD/AIS/…），不是归一化后的 IAC；
    归一化一律走 `resolve_stage`。
    """
    if dataset not in COHORTS:
        raise KeyError(f"{dataset} 不在册；在册 {sorted(COHORTS)}")
    path = os.path.join(GEO_META_DIR, f"{dataset}_samples.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"缺 GEO 权威样本表 {path} —— 拒绝无据冻结。")
    out = {}
    with open(path, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        need = {"gsm", "patient_id", "stage", "lesion_ordinal", "token", "geo_title"}
        if set(hdr) != need:
            raise ValueError(f"GEO 表表头异常 {hdr}，期望 {sorted(need)}")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            out[c[0]] = dict(patient_id=c[1], stage=c[2], lesion_ordinal=int(c[3]),
                             token=c[4], title=c[5])
    return out


def load_spatial_geo() -> dict:
    return load_geo("GSE307534")


def summary():
    print("=== 配对数据集（本项目范围）===")
    for k, v in COHORTS.items():
        print(f"  {k}: {v['modality']:8s} | {v['platform']} | GEO 样本 {v['n_samples_geo']}")
        print(f"      role: {v['role']}")
    sn = {g["patient_id"] for g in load_geo("GSE308103").values()}
    sp = {g["patient_id"] for g in load_geo("GSE307534").values()}
    both = sorted(sn & sp, key=lambda x: int(x[1:]))
    print(f"\n配对患者（{len(both)}，两模态交集）: {', '.join(both)}")
    print(f"  仅空间（无 snRNA，不入配对）: {', '.join(sorted(sp - sn, key=lambda x: int(x[1:]))) or '无'}")
    print(f"阶段词表: {STAGES}")
    print(f"LNM: {LNM_STATUS}")
    print("\n禁用法（P0）:")
    for f in FORBIDDEN:
        print(f"  ⛔ {f}")


if __name__ == "__main__":
    summary()
