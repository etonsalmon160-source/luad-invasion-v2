#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_ingest/01_freeze_paired.py — M0 输入冻结（两个配对数据集）

按 00_ingest/cohort_registry.py 的权威口径冻结：
  * GSE308103 —— snRNA（GEO 75 样本）
  * GSE307534 —— Visium 空间（GEO 56 样本）

枚举以 `00_ingest/geo_metadata/*_samples.tsv`（NCBI E-utilities 取得的 GEO 权威表）为准，
本地解压目录仅用于**定位文件**；两侧逐 GSM 交叉核验，不一致即计入 problems 并 FAIL。

产出：
  results/00_ingest/paired_samples.csv     样本表（dataset/sample_id/patient_id/stage/lesion_ordinal/modality）
  results/00_ingest/paired_patients.csv    患者表（双模态覆盖）
  results/00_ingest/paired_source_files.csv 源文件清单（路径/大小/mtime/SHA-256）
  results/00_ingest/M0_paired_validation_report.md  校验报告（含 GEO 交叉核验）
  results/00_ingest/paired_manifest.json   冻结清单（哈希）

铁律：分期**无静默默认**（resolve_stage 未知 token 即 raise）；切片定位**不靠目录名猜**
（候选 ≠ 1 即 raise）；不伪造缺失项。
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, os, re, sys
from datetime import datetime, timezone

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import cohort_registry as REG  # noqa: E402

DEFAULT_RAW = "/home/eto/luad_invasion/data"
DEFAULT_OUT = "/home/eto/luad_v2/results/00_ingest"
SCHEMA = "M0-paired-1"


def sha256_file(p: str) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def locate_slide(gsm_dir: str, gsm: str) -> tuple:
    """在 GSM 解压目录下定位**唯一**切片根目录（含 filtered_feature_bc_matrix/）。

    解压层级因 GEO 打包方式而异（有的直接摊平、有的多一层同名/异名目录），
    **不靠目录名猜**；候选 ≠ 1 直接 raise。返回 (root, missing_files)。
    """
    cands = []
    for dp, dns, _ in os.walk(gsm_dir):
        if "filtered_feature_bc_matrix" in dns:
            cands.append(dp)
    if len(cands) != 1:
        raise RuntimeError(
            f"{gsm}: 解压目录下 filtered_feature_bc_matrix 候选 {len(cands)} 个（应恰为 1）"
            f" → 拒绝猜测。候选={[os.path.relpath(c, gsm_dir) for c in cands]}")
    root = cands[0]
    missing = [rel for rel in REG.SPATIAL_REQUIRED if not os.path.exists(os.path.join(root, rel))]
    for grp in REG.SPATIAL_REQUIRED_ANY:
        if not any(os.path.exists(os.path.join(root, g)) for g in grp):
            missing.append(" 或 ".join(grp))
    return root, missing


def collect(raw: str, file_rows: list, problems: list, file_problems: list) -> pd.DataFrame:
    """按 **GEO 权威 GSM 表**枚举；本地解压目录仅用于定位文件。

    problems      = 命名/登记类不符（缺样本、多样本、目录名 vs GEO）
    file_problems = 切片文件完整性类问题
    """
    rows = []

    # --- GSE308103 (snRNA)：本地文件名 = <GSM>_<Pxx>_<Token>.raw_counts.mtx.txt.gz ---
    geo_sn = REG.load_geo("GSE308103")
    d = os.path.join(raw, "GSE308103", "extracted")
    local_sn = {}
    for f in sorted(os.listdir(d)):
        m = re.match(REG.COHORTS["GSE308103"]["sample_re"], f)
        if m:
            local_sn[m.group("gsm")] = (m.group("patient"), m.group("stage_token"), f)
    for gsm in sorted(geo_sn):
        g = geo_sn[gsm]
        if gsm not in local_sn:
            problems.append(f"GSE308103 {gsm}（{g['title']}）本地缺失")
            continue
        lpt, ltok, fname = local_sn[gsm]
        if (lpt, ltok) != (g["patient_id"], g["token"]):
            problems.append(f"GSE308103 {gsm} 本地名 ({lpt}_{ltok}) ≠ GEO ({g['patient_id']}_{g['token']})")
        stage = REG.resolve_stage(g["token"])                  # 未知即 raise
        sid = f"{g['patient_id']}_{g['token']}"
        rows.append(dict(dataset="GSE308103", sample_id=sid, patient_id=g["patient_id"],
                         stage=stage, stage_token=g["token"],
                         lesion_ordinal=g["lesion_ordinal"], modality="snRNA"))
        if gsm in local_sn:
            file_rows.append(dict(dataset="GSE308103", role="raw_counts", sample_id=sid,
                                  gsm=gsm, path=os.path.join(d, fname)))
    for gsm in sorted(set(local_sn) - set(geo_sn)):
        problems.append(f"GSE308103 本地多出 {gsm}（GEO 未登记）")

    # --- GSE307534 (spatial)：目录名 <GSM>_<Pxx>_<Token>，内含切片子目录 ---
    geo_sp = REG.load_geo("GSE307534")
    d = os.path.join(raw, "GSE307534", "extracted")
    local_sp = {}
    for name in sorted(os.listdir(d)):
        m = re.match(r"^(GSM\d+)_(P\d+)_(.+)$", name)
        if m and os.path.isdir(os.path.join(d, name)):
            local_sp[m.group(1)] = (name, m.group(2), m.group(3))
    for gsm in sorted(geo_sp):
        g = geo_sp[gsm]
        if gsm not in local_sp:
            problems.append(f"GSE307534 {gsm}（{g['title']}）本地缺失")
            continue
        name, lpt, ltok = local_sp[gsm]
        if (lpt, ltok) != (g["patient_id"], g["token"]):
            problems.append(f"GSE307534 {gsm} 目录名 ({lpt}_{ltok}) ≠ GEO ({g['patient_id']}_{g['token']})")
        stage = REG.resolve_stage(g["token"])
        sid = f"{g['patient_id']}_{g['token']}"
        rows.append(dict(dataset="GSE307534", sample_id=sid, patient_id=g["patient_id"],
                         stage=stage, stage_token=g["token"],
                         lesion_ordinal=g["lesion_ordinal"], modality="spatial"))
        try:
            root, missing = locate_slide(os.path.join(d, name), gsm)
        except RuntimeError as exc:
            file_problems.append(str(exc))
            continue
        if missing:
            file_problems.append(f"GSE307534 {gsm} 切片文件缺失 {missing}（root={os.path.relpath(root, raw)}）")
        file_rows.append(dict(dataset="GSE307534", role="spatial_slide", sample_id=sid,
                              gsm=gsm, path=root))
        for rel in REG.SPATIAL_REQUIRED:
            file_rows.append(dict(dataset="GSE307534", role="spatial_file", sample_id=sid,
                                  gsm=gsm, path=os.path.join(root, rel)))
        spdir = os.path.join(root, "spatial")
        if os.path.isdir(spdir):
            file_rows.append(dict(dataset="GSE307534", role="spatial_dir", sample_id=sid,
                                  gsm=gsm, path=spdir))
    for gsm in sorted(set(local_sp) - set(geo_sp)):
        problems.append(f"GSE307534 本地多出 {gsm}（GEO 未登记）")

    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-root", default=DEFAULT_RAW)
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    file_rows: list = []
    problems: list = []
    file_problems: list = []
    s = collect(args.raw_root, file_rows, problems, file_problems)
    # 唯一键：**sample_id 在本项目中跨数据集不唯一**（同一病灶的两个模态同名）→ 显式复合键
    s["sample_key"] = s["dataset"] + ":" + s["sample_id"]
    s = s[["sample_key", "dataset", "sample_id", "patient_id", "stage", "stage_token",
           "lesion_ordinal", "modality"]]
    s = s.sort_values(["dataset", "patient_id", "sample_id"])

    # 患者表：双模态覆盖
    p = (s.groupby(["patient_id"], as_index=False)
           .agg(modalities=("modality", lambda x: ";".join(sorted(set(x)))),
                n_samples=("sample_id", "nunique"),
                stages=("stage", lambda x: ";".join(sorted(set(x)))))
           .sort_values("patient_id"))
    p["paired"] = p["modalities"].apply(lambda m: set(m.split(";")) == {"snRNA", "spatial"})

    # ---- 校验 ----
    checks = []
    checks.append(("C1 两数据集均在册", set(s.dataset) == {"GSE308103", "GSE307534"},
                   str(sorted(set(s.dataset)))))
    checks.append(("C2 分期 ∈ 词表", set(s.stage) <= set(REG.STAGES), str(sorted(set(s.stage)))))
    checks.append(("C3 sample_key 全局唯一（= dataset:sample_id）",
                   not s.duplicated("sample_key").any(),
                   f"重复={int(s.duplicated('sample_key').sum())}；"
                   f"唯一 sample_id={s.sample_id.nunique()}/{len(s)}（跨数据集同名 {len(s)-s.sample_id.nunique()} 个=配对病灶）"))
    sn_pt = {g["patient_id"] for g in REG.load_geo("GSE308103").values()}
    sp_pt = {g["patient_id"] for g in REG.load_geo("GSE307534").values()}
    expect_paired = sorted(sn_pt & sp_pt, key=lambda x: int(x[1:]))
    paired = sorted(p.loc[p.paired, "patient_id"], key=lambda x: int(x[1:]))
    checks.append((f"C4 配对患者 == 两 GEO 表患者交集（且 ≥ {REG.PAIRED_PATIENTS_MIN} 例）",
                   paired == expect_paired and len(paired) >= REG.PAIRED_PATIENTS_MIN,
                   f"{len(paired)} 例={paired}"
                   + (f"；空间独有（无 snRNA，不入配对）：{sorted(sp_pt - sn_pt, key=lambda x: int(x[1:]))}")))
    checks.append(("C5 空转 LNM 未伪造", not (s.stage == "LNM").any(), "无 LNM 行"))

    # C6/C7：与 GEO 权威表逐 GSM 对齐（数量 + 逐样本名）
    cnt = {ds: int((s.dataset == ds).sum()) for ds in REG.COHORTS}
    exp = {ds: REG.COHORTS[ds]["n_samples_geo"] for ds in REG.COHORTS}
    checks.append((f"C6 样本数与 GEO 全量一致（{exp}）", cnt == exp, f"实测={cnt}"))
    checks.append(("C7 本地名 vs GEO 标题 逐样本核对（无缺、无多、无不符）",
                   len(problems) == 0,
                   "131 样本全部一致" if not problems else f"{len(problems)} 项，见报告 §GEO 交叉核验"))

    # C8：空转切片文件完整性（矩阵 + 坐标 + 缩放系数）
    checks.append(("C8 空转切片必需文件齐备（矩阵 + spatial 坐标/缩放）",
                   len(file_problems) == 0,
                   "56 切片全部齐备" if not file_problems
                   else f"{len(file_problems)} 张不完整，见报告 §切片完整性"))

    ok = all(c[1] for c in checks)
    for label, lst in (("命名/登记", problems), ("切片文件", file_problems)):
        if lst:
            print(f"[M0] {label} 问题：")
            for x in lst[:20]:
                print("   ⚠️", x)

    # ---- 输出 ----
    s.to_csv(os.path.join(args.out_dir, "paired_samples.csv"), index=False)
    p.to_csv(os.path.join(args.out_dir, "paired_patients.csv"), index=False)

    fr = pd.DataFrame(file_rows)
    sizes, hashes = [], []
    for r in fr.itertuples(index=False):
        if r.path and os.path.exists(r.path):
            st = os.stat(r.path)
            sizes.append(st.st_size if os.path.isfile(r.path) else None)
            hashes.append("")   # 大文件不逐个哈希（见 notes）
        else:
            sizes.append(None)
            hashes.append("MISSING")
    fr["size_bytes"], fr["sha256"] = sizes, hashes
    fr.to_csv(os.path.join(args.out_dir, "paired_source_files.csv"), index=False)

    manifest = {
        "schema": SCHEMA, "created_utc": utc(), "raw_root": args.raw_root,
        "registry_sha256": sha256_file(os.path.join(_HERE, "cohort_registry.py")),
        "geo_tables_sha256": {ds: sha256_file(os.path.join(REG.GEO_META_DIR, f"{ds}_samples.tsv"))
                              for ds in REG.COHORTS},
        "n_samples": int(len(s)), "n_patients": int(len(p)), "n_paired": int(p.paired.sum()),
        "gate_pass": bool(ok),
        "cross_check_problems": problems, "file_problems": file_problems,
        "tables": {k: sha256_file(os.path.join(args.out_dir, f))
                   for k, f in [("samples", "paired_samples.csv"),
                                ("patients", "paired_patients.csv"),
                                ("source_files", "paired_source_files.csv")]},
    }
    with open(os.path.join(args.out_dir, "paired_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False, sort_keys=True)

    rep = ["# M0 输入冻结（配对数据集）· 校验报告\n",
           f"- schema `{SCHEMA}` · {manifest['created_utc']} · raw `{args.raw_root}`",
           f"- registry sha256 `{manifest['registry_sha256']}`",
           f"- GEO 权威表 sha256 `{manifest['geo_tables_sha256']}`\n",
           f"## 过门判定：**{'✅ PASS' if ok else '❌ FAIL'}**\n",
           "| 检查 | 结果 | 详情 |", "| :-- | :--: | :-- |"]
    rep += [f"| {n} | {'✅' if o else '❌'} | {d} |" for n, o, d in checks]
    rep += ["\n## 计数\n",
            f"- 样本：**{len(s)}**　患者：**{len(p)}**　双模态配对：**{int(p.paired.sum())}**\n",
            "| 数据集 | 模态 | 样本数 | GEO 全量 | 患者数 |", "| :-- | :-- | --: | --: | --: |"]
    for ds, g in s.groupby("dataset"):
        rep.append(f"| {ds} | {g.modality.iloc[0]} | {len(g)} | "
                   f"{REG.COHORTS[ds]['n_samples_geo']} | {g.patient_id.nunique()} |")
    rep += ["\n## 配对患者核验\n", "| 患者 | 模态 | 样本数 | 阶段 |", "| :-- | :-- | --: | :-- |"]
    for r in p.itertuples(index=False):
        mark = "✅" if r.paired else "⚠️"
        rep.append(f"| {mark} {r.patient_id} | {r.modalities} | {r.n_samples} | {r.stages} |")

    rep += ["\n## GEO 交叉核验\n"]
    if problems:
        rep += [f"发现 **{len(problems)}** 项问题（**未静默**，逐条列出）：", ""]
        rep += [f"- ⚠️ {x}" for x in problems]
    else:
        rep += ["两个数据集**逐 GSM** 与 GEO 样本标题核对：**无缺、无多、无不符**。",
                "",
                "- 权威口径：`00_ingest/geo_metadata/{GSE308103,GSE307534}_samples.tsv`"
                "（由 NCBI E-utilities `esummary` 取得，含 GSM、患者号、分期、病灶序号、GEO 原始标题）。",
                "- 交叉核验内容：本地解压目录名/文件名推导的 `(patient, stage_token)` 是否与 GEO 标题一致；"
                "本地是否存在 GEO 未登记项或缺失项。",
                "",
                "| 数据集 | GEO 样本 | 本地命中 | 名不符 | 缺 | 多 |",
                "| :-- | --: | --: | --: | --: | --: |"]
        for ds in sorted(REG.COHORTS):
            g = REG.load_geo(ds)
            rep.append(f"| {ds} | {len(g)} | {len(g)} | 0 | 0 | 0 |")

    rep += ["\n## 切片完整性\n",
            "判据：每张空转切片必须含 "
            "`filtered_feature_bc_matrix/{matrix.mtx,barcodes.tsv,features.tsv}.gz` + "
            "`spatial/scalefactors_json.json` + `spatial/tissue_positions(.csv|_list.csv)`。\n"]
    if file_problems:
        rep += [f"发现 **{len(file_problems)}** 张不完整：", ""]
        rep += [f"- ⚠️ {x}" for x in file_problems]
        rep += ["", "> 处置：**重下该 GSM**（NCBI 端限速易致 tar 截断）；未补齐前 M0 不过门。"]
    else:
        rep += ["56 张切片全部齐备。\n"]

    rep += ["\n## 病灶序号（`lesion_ordinal`）\n",
            "`lesion_ordinal=2` 表示**同一患者的第二个独立病灶切片**（GEO 标题作 `second … of patient N`），"
            "与首个病灶**分期相同但解剖独立**，**不得合并**：\n",
            "| 数据集 | GSM | 患者 | token | GEO 标题 |", "| :-- | :-- | :-- | :-- | :-- |"]
    for ds in sorted(REG.COHORTS):
        for gsm, g in sorted(REG.load_geo(ds).items()):
            if g["lesion_ordinal"] == 2:
                rep.append(f"| {ds} | {gsm} | {g['patient_id']} | `{g['token']}` | {g['title']} |")

    rep += ["\n## 来源与哈希\n",
            f"- `paired_samples.csv` `{manifest['tables']['samples']}`",
            f"- `paired_patients.csv` `{manifest['tables']['patients']}`",
            f"- `paired_source_files.csv` `{manifest['tables']['source_files']}`\n",
            "## 说明",
            "- **唯一键 = `sample_key`（`<dataset>:<sample_id>`）**。`sample_id` 跨数据集**不唯一**"
            "（同一患者同一病灶的 snRNA 与空间切片同名，如 `P3_LUAD`）→ 下游 join **必须用 `sample_key`**。",
            "- 分期经 `cohort_registry.resolve_stage` **严格映射**；**原始 GEO 标签保留在 `stage_token`**"
            "（如 `LUAD` → 归一化 `IAC`，token 仍记 `LUAD`，便于回溯）。未知 token 直接 raise（无静默默认）。"
            "两侧命名惯例不同（snRNA 第二病灶用 `Normal1/AAH1/LUAD1/AIS1`，空间用 `Normal/AAH-1/LUAD-1/AIS-1`），"
            "**均按 GEO 原始文件名保真记录**，不强行统一。",
            "- 患者身份取自 **GEO 样本标题**（`… of patient N`）；两数据集的 `P*` 编号即 GEO 标题中的患者号"
            "（GEO 无独立 patient 字段）—— 配对关系据此建立。",
            "- **空转切片不靠目录名定位**：解压层级因 GEO 打包方式而异，"
            "一律由 `filtered_feature_bc_matrix/` 实际位置判定，候选 ≠ 1 即 raise（见 C8）。",
            "- 源文件 SHA-256 未逐个计算（矩阵较大）—— 完整性以 M1 的 per-cell 产物哈希为准。",
            "- 空转 LNM 暂缺（无合法 LUAD 数据）→ 不产出 LNM 行。"]

    with open(os.path.join(args.out_dir, "M0_paired_validation_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(rep) + "\n")

    print(f"[M0] gate={'PASS' if ok else 'FAIL'}  samples={len(s)} patients={len(p)} paired={int(p.paired.sum())}")
    for n, o, d in checks:
        print(f"     {'✅' if o else '❌'} {n}")
    print(f"     报告: {args.out_dir}/M0_paired_validation_report.md")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
