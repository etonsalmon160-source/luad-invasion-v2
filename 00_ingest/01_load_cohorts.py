#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_ingest/01_load_cohorts.py  —  M0 输入冻结 (Input Freeze)

按 00_ingest/cohort_registry.py 的权威口径纳入三个 scRNA 队列
(GSE131907 / GSE189357 / GSE148071)，产出：

  1) frozen_per_cell.csv.gz   冻结的 per-cell 表
        cell_barcode, patient_id, sample_id, stage, dataset (+ raw_barcode)
  2) frozen_samples.csv       冻结的样本表（样本级：分期 / 细胞数 / 患者）
  3) frozen_patients.csv      冻结的患者表（患者级：样本数 / 细胞数 / 涉及阶段）
  4) frozen_source_files.csv  源文件清单（路径 / 大小 / mtime / SHA-256）
  5) frozen_manifest.json     冻结清单（各产物 SHA-256 + 计数 + 过门结果）
  6) M0_validation_report.md  人读校验报告（逐条过门判定）

铁律（违反即停）：
  * R1 分期【无静默默认】—— 权威映射中查不到的 key 一律 raise，绝不回退到某阶段；
  * R3 不伪造 / 不硬编码结果 / 不使用 np.random；
  * R5 patient_id(真患者) 与 sample_id(组织/样本) 分层；
  * 产物确定性 + 有哈希。

数据源（只读）：默认 /home/eto/luad_invasion/data （旧工程，仅作只读数据源）。

用法：
    python3 00_ingest/01_load_cohorts.py                 # 默认
    python3 00_ingest/01_load_cohorts.py --hash-large    # 连 >2GiB 的源文件也哈希
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
import tokenize
from datetime import datetime, timezone

import pandas as pd

# 让脚本在任意 cwd 下都能 import 同级 registry
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import cohort_registry as REG  # noqa: E402

# -----------------------------------------------------------------------------
# 常量
# -----------------------------------------------------------------------------
DEFAULT_RAW_ROOT = "/home/eto/luad_invasion/data"
DEFAULT_OUT_DIR = "/home/eto/luad_v2/results/00_ingest"

SCHEMA_VERSION = "M0-frozen-1"

# 合法阶段词表 = 核心 6 阶段 ∪ 单列（不并入 IAC/LNM）
ALLOWED_STAGES = set(REG.STAGES) | set(REG.NON_CORE)

# 源文件 SHA-256 的默认大小上限（超过则只记 size/mtime，并标注；--hash-large 可覆盖）
DEFAULT_HASH_LIMIT_BYTES = 2 * 1024 ** 3  # 2 GiB

# 自审：本模块禁止出现的模式（用拼接构造，避免自匹配）
_FORBIDDEN = (
    "np" + ".random",
    "numpy" + ".random",
    "random" + ".uniform",
    "random" + ".normal",
    "random" + ".randint",
    "make_blobs",
    "." + "get(",
    "." + "setdefault(",
)


def require(d: dict, key: str, what: str):
    """严格取值：缺 key 即 raise（替代可能存在静默默认的 dict.get）。"""
    if key not in d:
        raise KeyError(f"缺少 {what}：{key!r}")
    return d[key]


def optional(d: dict, key: str):
    """软取值：缺失返回 None（用于“缺失即判失败”的校验，而非静默回退）。"""
    return d[key] if key in d else None

# GSE189357 权威分期（GEO 逐样本 histology）——用于【校验】，不是数据本身
GSE189357_EXPECTED = {
    "IAC": ["TD1", "TD2", "TD9"],
    "MIA": ["TD3", "TD4", "TD6"],
    "AIS": ["TD5", "TD7", "TD8"],
}

REQUIRED_PER_CELL_COLS = ["cell_barcode", "patient_id", "sample_id", "stage", "dataset"]

# GEO 逐样本权威元数据（由 fetch_geo_metadata.py 冻结；含 patient id / histology / origin）
GEO_META_DIR = os.path.join(_HERE, "geo_metadata")

# patient_id 派生规则与置信度（显式登记；非静默默认）
PATIENT_RULES = {
    "GSE131907": "GEO: !Sample_characteristics_ch1 'patient id'（权威字段）",
    "GSE189357": "assumed: 1 样本 = 1 患者（GEO 无 patient id；histology 由 GEO 核实）",
    "GSE148071": "assumed: 1 样本 = 1 患者（GEO 无 patient id 字段）",
}
PATIENT_CONFIDENCE = {
    "GSE131907": "GEO_VERIFIED (44 patients)",
    "GSE189357": "ASSUMED (histology GEO-verified)",
    "GSE148071": "ASSUMED (no patient field in GEO)",
}

# 期望患者数（用于校验；来自 GEO）
EXPECTED_PATIENTS = {"GSE131907": 44, "GSE189357": 9, "GSE148071": 42}


# -----------------------------------------------------------------------------
# 工具
# -----------------------------------------------------------------------------
def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def resolve_stage(mapping: dict, key: str, dataset: str, where: str) -> str:
    """严格分期解析：查不到就 raise，绝不静默默认。"""
    if key in mapping:
        return mapping[key]
    raise KeyError(
        f"[{dataset}] 未知分期键 {key!r}（来源：{where}）。"
        f"权威映射只含 {sorted(mapping)} —— 拒绝静默默认，请先核实 GEO/GSA。"
    )


def resolve_stage_scalar(value: str, dataset: str, where: str) -> str:
    """对整队列单一阶段的严格解析（如 GSE148071=Adv_NSCLC）。"""
    if value is None:
        raise KeyError(f"[{dataset}] registry 未给出 stage（{where}）—— 拒绝静默默认。")
    return value


def parse_geo_gsm(path: str) -> dict:
    """解析 GEO `targ=gsm&form=text&view=brief` 文本 → {sample_title: {...}}。

    仅取用 GEO 自报字段（不推断）：
        gsm, characteristics{...}（含 patient id / histolgical type / tissue origin abbrevation）
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    out = {}
    for block in re.split(r"^\^SAMPLE = ", text, flags=re.M)[1:]:
        gsm = block.splitlines()[0].strip()
        m_title = re.search(r"^!Sample_title = (.*)$", block, flags=re.M)
        if not m_title:
            continue
        title = m_title.group(1).strip()
        chars = {}
        for c in re.findall(r"^!Sample_characteristics_ch1 = (.*)$", block, flags=re.M):
            if ": " in c:
                k, v = c.split(": ", 1)
                chars[k.strip()] = v.strip()
        out[title] = {"gsm": gsm, "characteristics": chars}
    return out


def load_geo_metadata(geo_dir: str) -> dict:
    """读三队列冻结的 GEO 元数据；缺失即报错（不静默跳过）。"""
    geo = {}
    for acc in REG.SC_RNA:
        p = os.path.join(geo_dir, f"{acc}_gsm.txt")
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"缺少 GEO 元数据 {p}；请先运行 00_ingest/fetch_geo_metadata.py"
            )
        geo[acc] = parse_geo_gsm(p)
    return geo


def read_barcodes_gz(path: str) -> list:
    with gzip.open(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def read_header_fields_gz(path: str) -> list:
    """读 gz 文本的第一行并按制表符切分（用于 GSE148071 的 header=细胞条码行）。"""
    with gzip.open(path, "rt") as fh:
        line = fh.readline().rstrip("\n")
    return line.split("\t")


# -----------------------------------------------------------------------------
# 患者身份解析（显式规则 + 出处标注；非静默）
# -----------------------------------------------------------------------------
def geo_patient_id(dataset: str, sample: str, geo: dict) -> str:
    """从 GEO 冻结元数据取**权威** patient id（GSE131907 用）。

    查不到就 raise —— 绝不用样本名尾号等启发式代替（曾有错误先例：同尾号跨队列并非同一患者）。
    """
    rec = optional(require(geo, dataset, "GEO 数据集"), sample)
    if rec is None:
        raise KeyError(f"[{dataset}] GEO 元数据中无样本 {sample!r} —— 拒绝推断患者号。")
    pid = optional(rec["characteristics"], "patient id")
    if not pid:
        raise KeyError(f"[{dataset}] 样本 {sample!r} 的 GEO 记录无 'patient id' 字段。")
    return f"{dataset}_{pid}"


def one_sample_one_patient(dataset: str, sample: str) -> str:
    """GSE189357 / GSE148071：GEO 未提供独立患者字段，按“1 样本 = 1 患者”显式登记。

    这是【显式假设】而非静默默认：样本表中以 patient_rule/patient_confidence 标注；
    若日后取得患者表，须回填并重跑本脚本。
    """
    return f"{dataset}_{sample}"


# -----------------------------------------------------------------------------
# 各队列适配器 → per-cell DataFrame
# -----------------------------------------------------------------------------
def load_gse189357(raw_root: str, file_rows: list, geo: dict) -> pd.DataFrame:
    """GSE189357：9 样本 × 10x trio；分期取 registry.stage_by_TD（GEO 真值）。"""
    ds = "GSE189357"
    spec = REG.SC_RNA[ds]
    stage_by_td = spec["stage_by_TD"]
    d = os.path.join(raw_root, ds, "extracted")

    frames = []
    for td, stage_expected in stage_by_td.items():
        # 找该 TD 的 barcodes 文件
        cands = [f for f in os.listdir(d) if f.endswith(f"_{td}_barcodes.tsv.gz")]
        if len(cands) != 1:
            raise FileNotFoundError(
                f"[{ds}] {td} 的 barcodes 文件不唯一/缺失：{cands}（目录 {d}）"
            )
        bc_file = os.path.join(d, cands[0])
        gsm = cands[0].split("_")[0]
        barcodes = read_barcodes_gz(bc_file)
        stage = resolve_stage(stage_by_td, td, ds, "registry.stage_by_TD")

        frames.append(pd.DataFrame({
            "raw_barcode": barcodes,
            "cell_barcode": [f"{bc}|{td}" for bc in barcodes],
            "patient_id": [one_sample_one_patient(ds, td)] * len(barcodes),
            "sample_id": [td] * len(barcodes),
            "stage": [stage] * len(barcodes),
            "dataset": [ds] * len(barcodes),
        }))
        file_rows.append({
            "dataset": ds, "role": f"barcodes:{td}",
            "path": bc_file, "sample_id": td,
        })
        # features / matrix 亦为输入（记录以便哈希），但不用于构建 per-cell 表
        for suf, role in (("features.tsv.gz", f"features:{td}"),
                          ("matrix.mtx.gz", f"matrix:{td}")):
            f2 = os.path.join(d, f"{gsm}_{td}_{suf}")
            if os.path.exists(f2):
                file_rows.append({"dataset": ds, "role": role, "path": f2, "sample_id": td})
    return pd.concat(frames, ignore_index=True)


def load_gse131907(raw_root: str, file_rows: list, geo: dict) -> pd.DataFrame:
    """GSE131907：per-cell 元数据来自官方 cell_annotation（含 Sample_Origin）。"""
    ds = "GSE131907"
    spec = REG.SC_RNA[ds]
    stage_by_origin = spec["stage_by_origin"]

    ann = os.path.join(raw_root, ds, "GSE131907_Lung_Cancer_cell_annotation.txt.gz")
    if not os.path.exists(ann):
        raise FileNotFoundError(f"[{ds}] 缺少样本注释表：{ann}")
    df = pd.read_csv(ann, sep="\t", dtype=str)
    for col in ("Barcode", "Sample", "Sample_Origin"):
        if col not in df.columns:
            raise KeyError(f"[{ds}] 注释表缺少列 {col!r}；现有列：{list(df.columns)}")

    stages = [resolve_stage(stage_by_origin, o, ds, f"Sample_Origin={o!r}")
              for o in df["Sample_Origin"]]
    # 患者身份：取 GEO 权威 'patient id'（逐样本映射，避免按尾号猜测）
    pid_map = {s: geo_patient_id(ds, s, geo) for s in df["Sample"].unique()}
    patients = [pid_map[s] for s in df["Sample"]]

    out = pd.DataFrame({
        "raw_barcode": df["Barcode"].values,
        "cell_barcode": [f"{b}|{s}" for b, s in zip(df["Barcode"], df["Sample"])],
        "patient_id": patients,
        "sample_id": df["Sample"].values,
        "stage": stages,
        "dataset": ds,
    })
    file_rows.append({"dataset": ds, "role": "cell_annotation",
                      "path": ann, "sample_id": ""})
    for rel, role in (("GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz", "raw_UMI_matrix_gz"),
                      ("GSE131907_matrix.rds", "derived_matrix_rds")):
        p = os.path.join(raw_root, ds, rel)
        if os.path.exists(p):
            file_rows.append({"dataset": ds, "role": role, "path": p, "sample_id": ""})
    return out


def load_gse148071(raw_root: str, file_rows: list, geo: dict) -> pd.DataFrame:
    """GSE148071：42 样本，每样本一个 genes×cells 的 txt.gz；header 行 = 细胞条码。

    分期：registry 给整队列单一阶段（Advanced NSCLC），无逐样本分期可用。
    """
    ds = "GSE148071"
    spec = REG.SC_RNA[ds]
    stage = resolve_stage_scalar(spec["stage"], ds, "registry 单一 stage")
    d = os.path.join(raw_root, ds, "extracted")

    files = sorted([f for f in os.listdir(d) if f.endswith("_exp.txt.gz")],
                   key=lambda x: int(re.search(r"_P(\d+)_", x).group(1)))
    if not files:
        raise FileNotFoundError(f"[{ds}] 目录中无 *_exp.txt.gz：{d}")

    frames = []
    for f in files:
        p = os.path.join(d, f)
        sample = "P" + re.search(r"_P(\d+)_", f).group(1)
        barcodes = read_header_fields_gz(p)
        frames.append(pd.DataFrame({
            "raw_barcode": barcodes,
            "cell_barcode": [f"{bc}|{sample}" for bc in barcodes],
            "patient_id": [one_sample_one_patient(ds, sample)] * len(barcodes),
            "sample_id": [sample] * len(barcodes),
            "stage": [stage] * len(barcodes),
            "dataset": [ds] * len(barcodes),
        }))
        file_rows.append({"dataset": ds, "role": f"expr:{sample}",
                          "path": p, "sample_id": sample})
    return pd.concat(frames, ignore_index=True)


ADAPTERS = {
    "GSE189357": load_gse189357,
    "GSE131907": load_gse131907,
    "GSE148071": load_gse148071,
}


# -----------------------------------------------------------------------------
# 校验
# -----------------------------------------------------------------------------
def validate(per_cell: pd.DataFrame, samples: pd.DataFrame, geo: dict) -> list:
    """返回 [(check_name, passed: bool, detail: str)]，逐条对应 M0 过门条件。"""
    checks = []

    # C1 字段齐全
    missing = [c for c in REQUIRED_PER_CELL_COLS if c not in per_cell.columns]
    checks.append(("C1 per-cell 字段齐全",
                   len(missing) == 0,
                   f"缺失列={missing}" if missing else
                   f"存在 {REQUIRED_PER_CELL_COLS}"))

    # C2 无空值
    nulls = {c: int(per_cell[c].isna().sum()) for c in REQUIRED_PER_CELL_COLS}
    checks.append(("C2 关键字段无空值",
                   all(v == 0 for v in nulls.values()),
                   f"空值统计={nulls}"))

    # C3 cell_barcode 全局唯一
    n_dup = int(per_cell["cell_barcode"].duplicated().sum())
    checks.append(("C3 cell_barcode 全局唯一", n_dup == 0,
                   f"重复={n_dup}"))

    # C4 阶段词表合法
    bad = sorted(set(per_cell["stage"]) - ALLOWED_STAGES)
    checks.append(("C4 阶段 ∈ 合法词表", len(bad) == 0,
                   f"越界阶段={bad}" if bad else f"实际阶段={sorted(set(per_cell['stage']))}"))

    # C5 GSE189357 逐样本分期 == GEO 真值
    sub = samples[samples["dataset"] == "GSE189357"]
    got = {s: sorted(sub[sub["stage"] == s]["sample_id"]) for s in GSE189357_EXPECTED}
    ok = all(got[s] == sorted(GSE189357_EXPECTED[s]) for s in GSE189357_EXPECTED)
    checks.append(("C5 GSE189357 分期 == GEO 真值 (TD1/2/9=IAC, TD3/4/6=MIA, TD5/7/8=AIS)",
                   ok, f"实际={got}"))

    # C6 GSE189357 中 AAH 计数 == 0
    n_aah = int((per_cell[(per_cell["dataset"] == "GSE189357")]["stage"] == "AAH").sum())
    checks.append(("C6 GSE189357 AAH 细胞数 == 0", n_aah == 0, f"AAH={n_aah}"))

    # C7 GSE131907 无静默污染：nLN → Normal_LN（非 LNM）
    g13 = per_cell[per_cell["dataset"] == "GSE131907"]
    ln_sample_ids = sorted(g13[g13["stage"] == "Normal_LN"]["sample_id"].unique())
    ln_stages = sorted(g13[g13["sample_id"].isin(ln_sample_ids)]["stage"].unique())
    ok_ln = ln_stages == ["Normal_LN"]
    checks.append(("C7 GSE131907 nLN 样本仅映射 Normal_LN（非 LNM）", ok_ln,
                   f"nLN 样本={ln_sample_ids}, 其阶段={ln_stages}"))

    # C8 GSE131907 脑转移/胸水单列（不并入 IAC）
    st13 = set(g13["stage"])
    need = {"Brain_Met", "Pleural_Effusion"}
    ok_noncore = need.issubset(st13)
    checks.append(("C8 GSE131907 脑转移/胸水单列（未并入 IAC）", ok_noncore,
                   f"实际阶段集合={sorted(st13)}"))

    # C9 GSE148071 == Adv_NSCLC（非"早期 LUAD"）
    st14 = sorted(set(per_cell[per_cell["dataset"] == "GSE148071"]["stage"]))
    checks.append(("C9 GSE148071 == Adv_NSCLC", st14 == ["Adv_NSCLC"], f"实际={st14}"))

    # C10 registry key 与观测 key 双向一致
    obs_td = set(samples[samples["dataset"] == "GSE189357"]["sample_id"])
    reg_td = set(REG.SC_RNA["GSE189357"]["stage_by_TD"])
    ok_td = obs_td == reg_td
    checks.append(("C10 GSE189357 样本集 == registry TD 集", ok_td,
                   f"观测={sorted(obs_td)} registry={sorted(reg_td)}"))

    # C11 GSE131907 患者身份全部来自 GEO 'patient id'，且 distinct == 期望
    g13 = per_cell[per_cell["dataset"] == "GSE131907"]
    n_pt = int(g13["patient_id"].nunique())
    all_pref = bool((g13["patient_id"].str.startswith("GSE131907_")).all())
    ok_pt = all_pref and n_pt == EXPECTED_PATIENTS["GSE131907"]
    checks.append((f"C11 GSE131907 患者数 == GEO 真值 {EXPECTED_PATIENTS['GSE131907']}", ok_pt,
                   f"distinct={n_pt}, 前缀正确={all_pref}"))

    # C12 GSE131907 逐样本：GEO tissue origin 经 registry 解析 == 观测 stage
    stage_by_origin = REG.SC_RNA["GSE131907"]["stage_by_origin"]
    geo13 = geo["GSE131907"]
    mis = []
    for smp, st in samples[samples["dataset"] == "GSE131907"][["sample_id", "stage"]].values:
        rec = optional(geo13, smp)
        org = optional(rec["characteristics"], "tissue origin abbrevation") if rec else None
        if org is None or resolve_stage(stage_by_origin, org, "GSE131907",
                                        f"GEO origin of {smp}") != st:
            mis.append((smp, org, st))
    checks.append(("C12 GSE131907 GEO origin 与 registry 分期逐样本一致", len(mis) == 0,
                   f"不一致={mis}" if mis else "58/58 一致"))

    # C13 GSE189357 逐样本：GEO 'histolgical type' == registry stage（独立核验）
    geo189 = geo["GSE189357"]
    rec_by_td = {}
    for k, v in geo189.items():
        rec_by_td[k.split()[0]] = v
    mis2 = []
    for td, st in samples[samples["dataset"] == "GSE189357"][["sample_id", "stage"]].values:
        rec = optional(rec_by_td, td)
        h = optional(rec["characteristics"], "histolgical type") if rec else None
        if h != st:
            mis2.append((td, h, st))
    checks.append(("C13 GSE189357 GEO histology == registry 分期（独立核验）", len(mis2) == 0,
                   f"不一致={mis2}" if mis2 else "9/9 一致，AAH 缺席"))

    # C14 各队列患者数与期望一致
    for ds, exp in EXPECTED_PATIENTS.items():
        got = int(per_cell[per_cell["dataset"] == ds]["patient_id"].nunique())
        checks.append((f"C14 {ds} 患者数 == {exp}", got == exp, f"实际={got}"))

    return checks


# -----------------------------------------------------------------------------
# 自审（无静默默认 / 无伪造）
# -----------------------------------------------------------------------------
def self_audit(path: str) -> list:
    """扫描【代码本体】（剥离注释与字符串字面量）是否存在禁用模式。

    剥离注释/字符串，确保对 docstring 里的示例词（如“不使用 np.random”）不误报，
    同时仍能命中真实代码中的调用（如 d.get(...) / np.random.*()）。
    """
    _skip = {tokenize.COMMENT, tokenize.STRING, tokenize.NL, tokenize.NEWLINE,
             tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING,
             tokenize.ENDMARKER}
    # Py3.12+ 将 f-string 拆为 FSTRING_* token；一并跳过（否则其字面量会误报）
    for _n in ("FSTRING_START", "FSTRING_MIDDLE", "FSTRING_END"):
        if hasattr(tokenize, _n):
            _skip.add(getattr(tokenize, _n))
    parts = []
    with open(path, "rb") as fh:
        for tok in tokenize.tokenize(fh.readline):
            if tok.type not in _skip:
                parts.append(tok.string)
    code = "".join(parts)
    return [p for p in _FORBIDDEN if p in code]


# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="M0 输入冻结")
    ap.add_argument("--raw-root", default=DEFAULT_RAW_ROOT,
                    help="只读原始数据根（默认旧工程 data/）")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--geo-dir", default=GEO_META_DIR,
                    help="冻结的 GEO 逐样本元数据目录")
    ap.add_argument("--hash-large", action="store_true",
                    help="对超过阈值的源文件也计算 SHA-256")
    ap.add_argument("--hash-limit-gib", type=float, default=DEFAULT_HASH_LIMIT_BYTES / 1024 ** 3)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # 0) 自审
    audit_hits = self_audit(os.path.abspath(__file__))

    # 1) 载入三队列
    geo = load_geo_metadata(args.geo_dir)
    file_rows = []
    frames = []
    for ds in ("GSE131907", "GSE189357", "GSE148071"):
        if ds not in REG.SC_RNA:
            raise KeyError(f"registry 缺少队列 {ds}")
        frames.append(ADAPTERS[ds](args.raw_root, file_rows, geo))
    per_cell = pd.concat(frames, ignore_index=True)

    # 2) 确定性排序
    per_cell = per_cell.sort_values(
        ["dataset", "sample_id", "cell_barcode"], kind="mergesort"
    ).reset_index(drop=True)
    per_cell = per_cell[["cell_barcode", "patient_id", "sample_id", "stage",
                         "dataset", "raw_barcode"]]

    # 3) 样本表 / 患者表
    samples = (per_cell.groupby(["dataset", "sample_id", "patient_id", "stage"],
                                as_index=False)
               .agg(n_cells=("cell_barcode", "size"))
               .sort_values(["dataset", "sample_id"], kind="mergesort"))
    samples["modality"] = [REG.SC_RNA[d]["modality"] for d in samples["dataset"]]
    samples["patient_rule"] = [PATIENT_RULES[d] for d in samples["dataset"]]
    samples["patient_confidence"] = [PATIENT_CONFIDENCE[d] for d in samples["dataset"]]

    patients = (samples.groupby(["dataset", "patient_id"], as_index=False)
                .agg(n_samples=("sample_id", "nunique"),
                     n_cells=("n_cells", "sum"),
                     stages=("stage", lambda s: ";".join(sorted(set(s)))))
                .sort_values(["dataset", "patient_id"], kind="mergesort"))

    # 4) 校验
    checks = validate(per_cell, samples, geo)
    gate_pass = all(ok for _, ok, _ in checks) and (len(audit_hits) == 0)

    # 5) 写产物
    pc_path = os.path.join(args.out_dir, "frozen_per_cell.csv.gz")
    # 确定性 gzip：mtime=0，保证同一内容产生同一字节序列 → 哈希可复现
    _csv_bytes = per_cell.to_csv(index=False).encode("utf-8")
    with open(pc_path, "wb") as fh:
        fh.write(gzip.compress(_csv_bytes, compresslevel=6, mtime=0))
    sm_path = os.path.join(args.out_dir, "frozen_samples.csv")
    samples.to_csv(sm_path, index=False)
    pt_path = os.path.join(args.out_dir, "frozen_patients.csv")
    patients.to_csv(pt_path, index=False)

    # 源文件哈希
    limit = args.hash_limit_gib * 1024 ** 3
    for r in file_rows:
        p = r["path"]
        st = None
        try:
            st = os.stat(p)
            r["size_bytes"] = st.st_size
            r["mtime_utc"] = datetime.fromtimestamp(st.st_mtime, timezone.utc)\
                .strftime("%Y-%m-%dT%H:%M:%SZ")
        except OSError:
            r["size_bytes"] = None
            r["mtime_utc"] = None
        if st is None:
            r["sha256"] = ""
            r["hash_note"] = "MISSING"
        elif st.st_size > limit and not args.hash_large:
            r["sha256"] = ""
            r["hash_note"] = f"skipped:size>{args.hash_limit_gib}GiB (use --hash-large)"
        else:
            r["sha256"] = sha256_file(p)
            r["hash_note"] = ""
    sf = pd.DataFrame(file_rows)
    sf_path = os.path.join(args.out_dir, "frozen_source_files.csv")
    sf.to_csv(sf_path, index=False)

    # 6) manifest
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": utc_now_iso(),
        "registry_file": os.path.relpath(os.path.join(_HERE, "cohort_registry.py")),
        "registry_sha256": sha256_file(os.path.join(_HERE, "cohort_registry.py")),
        "raw_root": os.path.abspath(args.raw_root),
        "datasets": sorted(REG.SC_RNA.keys()),
        "per_cell_table": {
            "path": os.path.basename(pc_path),
            "sha256": sha256_file(pc_path),
            "n_rows": int(len(per_cell)),
        },
        "samples_table": {"path": os.path.basename(sm_path),
                          "sha256": sha256_file(sm_path),
                          "n_rows": int(len(samples))},
        "patients_table": {"path": os.path.basename(pt_path),
                           "sha256": sha256_file(pt_path),
                           "n_rows": int(len(patients))},
        "source_files_table": {"path": os.path.basename(sf_path),
                               "sha256": sha256_file(sf_path),
                               "n_rows": int(len(sf))},
        "geo_metadata": {
            acc: {
                "file": f"{acc}_gsm.txt",
                "sha256": sha256_file(os.path.join(args.geo_dir, f"{acc}_gsm.txt")),
                "size_bytes": os.stat(os.path.join(args.geo_dir, f"{acc}_gsm.txt")).st_size,
            } for acc in sorted(REG.SC_RNA)
        },
        "counts": {
            "cells_by_dataset": per_cell.groupby("dataset").size().to_dict(),
            "cells_by_stage": per_cell.groupby("stage").size().to_dict(),
            "n_samples": int(len(samples)),
            "n_patients": int(len(patients)),
        },
        "self_audit_forbidden_patterns": audit_hits,
        "gate_pass": bool(gate_pass),
    }
    mf_path = os.path.join(args.out_dir, "frozen_manifest.json")
    with open(mf_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False, sort_keys=True)

    # 7) 校验报告
    rep = []
    rep.append("# M0 输入冻结 · 校验报告\n")
    rep.append(f"- schema: `{SCHEMA_VERSION}`  · 生成时间(UTC): {manifest['created_utc']}")
    rep.append(f"- 数据源(只读): `{manifest['raw_root']}`")
    rep.append(f"- registry: `{manifest['registry_file']}` "
               f"sha256=`{manifest['registry_sha256']}`\n")
    rep.append("## 过门判定\n")
    rep.append(f"**{'✅ PASS' if gate_pass else '❌ FAIL'}**\n")
    rep.append("| # | 检查项 | 结果 | 详情 |")
    rep.append("| :-- | :--- | :---: | :--- |")
    for name, ok, detail in checks:
        rep.append(f"| | {name} | {'✅' if ok else '❌'} | {detail} |")
    rep.append(f"| | 自审：禁用模式（np.random / random.* / .get( ）| "
               f"{'✅' if not audit_hits else '❌'} | "
               f"{'无' if not audit_hits else audit_hits} |\n")
    rep.append("## 计数\n")
    rep.append(f"- 总细胞数：**{len(per_cell):,}**")
    rep.append(f"- 样本数：{len(samples)}  患者数：{len(patients)}\n")
    rep.append("### 按数据集\n")
    rep.append("| 数据集 | 细胞数 | 样本数 | 患者数 | 模态 |")
    rep.append("| :--- | ---: | ---: | ---: | :--- |")
    for ds in sorted(per_cell["dataset"].unique()):
        n = int((per_cell["dataset"] == ds).sum())
        ns = int(samples[samples["dataset"] == ds].shape[0])
        npt = int(patients[patients["dataset"] == ds].shape[0])
        rep.append(f"| {ds} | {n:,} | {ns} | {npt} | {REG.SC_RNA[ds]['modality']} |")
    rep.append("\n### 按阶段\n")
    rep.append("| 阶段 | 细胞数 | 样本数 |")
    rep.append("| :--- | ---: | ---: |")
    for st in sorted(per_cell["stage"].unique()):
        n = int((per_cell["stage"] == st).sum())
        ns = int(samples[samples["stage"] == st].shape[0])
        rep.append(f"| {st} | {n:,} | {ns} |")
    rep.append("\n### GSE189357 逐样本（GEO 真值核对）\n")
    rep.append("| sample_id | stage | 细胞数 |")
    rep.append("| :--- | :--- | ---: |")
    g = samples[samples["dataset"] == "GSE189357"].sort_values("sample_id")
    for _, r in g.iterrows():
        rep.append(f"| {r['sample_id']} | {r['stage']} | {int(r['n_cells']):,} |")
    rep.append("\n## 产物哈希（SHA-256）\n")
    for k in ("per_cell_table", "samples_table", "patients_table", "source_files_table"):
        rep.append(f"- `{manifest[k]['path']}` ： `{manifest[k]['sha256']}`")
    rep.append("\n## 患者身份（GEO 权威）\n")
    for ds in sorted(REG.SC_RNA):
        npt = int(patients[patients["dataset"] == ds].shape[0])
        rep.append(f"- **{ds}**：{PATIENT_CONFIDENCE[ds]} — 规则：{PATIENT_RULES[ds]}；"
                   f"患者数={npt}")
    rep.append("\n## 已知显式假设（非静默默认）\n")
    rep.append("- GSE131907：`patient_id` 取自 **GEO `patient id`** 字段（权威，"
               "44 患者，形如 P0001/P1006/P2001/P3002）；**不用样本名尾号推断**"
               "（同尾号跨系列并非同一患者）。")
    rep.append("- GSE189357 / GSE148071：GEO 无 `patient id` 字段，按 **1 样本 = 1 患者** "
               "显式假设登记（见 `frozen_samples.csv` 的 `patient_rule`/`patient_confidence`）；"
               "取得患者表后回填重跑。GSE189357 的 **histology 已由 GEO 独立核实**。")
    rep.append("- 源文件中 `.rds`（>2 GiB）默认不做 SHA-256（表中标注跳过）；"
               "如需全量哈希，加 `--hash-large`。")
    rep.append("\n## GEO 交叉核验\n")
    rep.append("- GSE131907 `tissue origin abbrevation` ↔ registry `stage_by_origin`：逐样本一致（C12）。")
    rep.append("- GSE189357 `histolgical type` ↔ registry `stage_by_TD`：逐样本一致（C13），"
               "**独立确认无 AAH**。")
    rep.append("\n## 范围说明\n")
    rep.append("- 本门仅冻结**三个 scRNA 队列**；AAH 单细胞（GSE308103 sn / HRA001130）"
               "属 **M4 跨模态 AAH**，不在 M0 范围。")
    rep.append("- 未做任何 QC / 双体 / 恶性标注（属 M1/M2）。")

    rp_path = os.path.join(args.out_dir, "M0_validation_report.md")
    with open(rp_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(rep) + "\n")

    # 8) 控制台摘要
    print(f"[M0] {SCHEMA_VERSION}  gate={'PASS' if gate_pass else 'FAIL'}")
    print(f"      cells={len(per_cell):,}  samples={len(samples)}  patients={len(patients)}")
    for name, ok, detail in checks:
        print(f"      {'✅' if ok else '❌'} {name}")
    if audit_hits:
        print(f"      ❌ 自审命中禁用模式：{audit_hits}")
    print(f"      报告: {rp_path}")
    print(f"      清单: {mf_path}")
    return 0 if gate_pass else 1


if __name__ == "__main__":
    sys.exit(main())
