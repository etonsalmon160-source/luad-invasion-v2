#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_ingest/01_freeze_paired.py — M0 输入冻结（两个配对数据集）

按 00_ingest/cohort_registry.py 的权威口径冻结：
  * GSE308103 —— snRNA（75 样本）
  * GSE307534 —— Visium 空间（GEO 56 样本；本地已解压若干）

产出：
  results/00_ingest/paired_samples.csv     样本表（dataset/sample_id/patient_id/stage/modality）
  results/00_ingest/paired_patients.csv    患者表（双模态覆盖）
  results/00_ingest/paired_source_files.csv 源文件清单（路径/大小/mtime/SHA-256）
  results/00_ingest/M0_paired_validation_report.md  校验报告（含 9 例配对核验）
  results/00_ingest/paired_manifest.json   冻结清单（哈希）

铁律：分期**无静默默认**（resolve_stage 未知 token 即 raise）；不伪造缺失项。
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


def collect(raw: str, file_rows: list) -> pd.DataFrame:
    rows = []
    # --- GSE308103 (snRNA) ---
    d = os.path.join(raw, "GSE308103", "extracted")
    for f in sorted(os.listdir(d)):
        m = re.match(REG.COHORTS["GSE308103"]["sample_re"], f)
        if not m:
            continue
        patient, token = m.group("patient"), m.group("stage_token")
        stage = REG.resolve_stage(token)                       # 未知即 raise
        rows.append(dict(dataset="GSE308103", sample_id=f"{patient}_{token}",
                         patient_id=patient, stage=stage, stage_token=token,
                         modality="snRNA"))
        file_rows.append(dict(dataset="GSE308103", role="raw_counts", sample_id=f"{patient}_{token}",
                              gsm=m.group("gsm"), path=os.path.join(d, f)))
    # --- GSE307534 (spatial) ---
    d = os.path.join(raw, "GSE307534", "extracted")
    for name in sorted(os.listdir(d)):
        m = re.match(r"^(GSM\d+)_(P\d+)_(.+)$", name)
        if not (m and os.path.isdir(os.path.join(d, name))):
            continue
        patient, token = m.group(2), m.group(3)
        stage = REG.resolve_stage(token)
        sid = f"{patient}_{token}"
        rows.append(dict(dataset="GSE307534", sample_id=sid, patient_id=patient,
                         stage=stage, stage_token=token, modality="spatial"))
        mtx = os.path.join(d, name, token, "filtered_feature_bc_matrix", "matrix.mtx.gz")
        file_rows.append(dict(dataset="GSE307534", role="spatial_matrix", sample_id=sid,
                              gsm=m.group(1), path=mtx if os.path.exists(mtx) else ""))
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-root", default=DEFAULT_RAW)
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    file_rows: list = []
    s = collect(args.raw_root, file_rows)
    # 唯一键：**sample_id 在本项目中跨数据集不唯一**（同一病灶的两个模态同名）→ 显式复合键
    s["sample_key"] = s["dataset"] + ":" + s["sample_id"]
    s = s[["sample_key", "dataset", "sample_id", "patient_id", "stage", "stage_token", "modality"]]
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
    paired = sorted(p.loc[p.paired, "patient_id"])
    checks.append((f"C4 配对患者 == 登记表 {len(REG.PAIRED_PATIENTS)} 例",
                   paired == sorted(REG.PAIRED_PATIENTS), f"实测={paired}"))
    checks.append(("C5 空转 LNM 未伪造", not (s.stage == "LNM").any(), "无 LNM 行"))

    ok = all(c[1] for c in checks)

    # ---- 输出 ----
    s.to_csv(os.path.join(args.out_dir, "paired_samples.csv"), index=False)
    p.to_csv(os.path.join(args.out_dir, "paired_patients.csv"), index=False)

    fr = pd.DataFrame(file_rows)
    for r in fr.itertuples(index=False):
        if r.path and os.path.exists(r.path):
            st = os.stat(r.path)
            fr.loc[fr.path == r.path, "size_bytes"] = st.st_size
            fr.loc[fr.path == r.path, "sha256"] = ""   # 大文件不逐个哈希（见 notes）
        else:
            fr.loc[fr.path == r.path, "size_bytes"] = None
            fr.loc[fr.path == r.path, "sha256"] = "MISSING"
    fr.to_csv(os.path.join(args.out_dir, "paired_source_files.csv"), index=False)

    manifest = {
        "schema": SCHEMA, "created_utc": utc(), "raw_root": args.raw_root,
        "registry_sha256": sha256_file(os.path.join(_HERE, "cohort_registry.py")),
        "n_samples": int(len(s)), "n_patients": int(len(p)), "n_paired": int(p.paired.sum()),
        "gate_pass": bool(ok),
        "tables": {k: sha256_file(os.path.join(args.out_dir, f))
                   for k, f in [("samples", "paired_samples.csv"),
                                ("patients", "paired_patients.csv"),
                                ("source_files", "paired_source_files.csv")]},
    }
    with open(os.path.join(args.out_dir, "paired_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False, sort_keys=True)

    rep = ["# M0 输入冻结（配对数据集）· 校验报告\n",
           f"- schema `{SCHEMA}` · {manifest['created_utc']} · raw `{args.raw_root}`",
           f"- registry sha256 `{manifest['registry_sha256']}`\n",
           f"## 过门判定：**{'✅ PASS' if ok else '❌ FAIL'}**\n",
           "| 检查 | 结果 | 详情 |", "| :-- | :--: | :-- |"]
    rep += [f"| {n} | {'✅' if o else '❌'} | {d} |" for n, o, d in checks]
    rep += ["\n## 计数\n",
            f"- 样本：**{len(s)}**　患者：**{len(p)}**　双模态配对：**{int(p.paired.sum())}**\n",
            "| 数据集 | 模态 | 样本数 | 患者数 |", "| :-- | :-- | --: | --: |"]
    for ds, g in s.groupby("dataset"):
        rep.append(f"| {ds} | {g.modality.iloc[0]} | {len(g)} | {g.patient_id.nunique()} |")
    rep += ["\n## 配对患者核验\n", "| 患者 | 模态 | 样本数 | 阶段 |", "| :-- | :-- | --: | :-- |"]
    for r in p.itertuples(index=False):
        mark = "✅" if r.paired else "⚠️"
        rep.append(f"| {mark} {r.patient_id} | {r.modalities} | {r.n_samples} | {r.stages} |")
    rep += ["\n## 来源与哈希\n",
            f"- `paired_samples.csv` `{manifest['tables']['samples']}`",
            f"- `paired_patients.csv` `{manifest['tables']['patients']}`",
            f"- `paired_source_files.csv` `{manifest['tables']['source_files']}`\n",
            "## 说明",
            "- **唯一键 = `sample_key`（`<dataset>:<sample_id>`）**。`sample_id` 跨数据集**不唯一**"
            "（同一患者同一病灶的 snRNA 与空间切片同名，如 `P3_LUAD`）→ 下游 join **必须用 `sample_key`**。",
            "- 分期经 `cohort_registry.resolve_stage` **严格映射**；**原始 GEO 标签保留在 `stage_token`**"
            "（如 `LUAD` → 归一化 `IAC`，token 仍记 `LUAD`，便于回溯）。未知 token 直接 raise（无静默默认）。",
            "- 患者身份取自 **GEO 样本标题**（`… of patient N`）；两数据集的 `P*` 编号即 GEO 标题中的患者号",
            "（GEO 无独立 patient 字段）—— 配对关系据此建立，见 §GEO 交叉核验。",
            "- GSE307534 的 GEO 全量为 56 样本；**本地已解压**的计入本表（其余待下载）。",
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
