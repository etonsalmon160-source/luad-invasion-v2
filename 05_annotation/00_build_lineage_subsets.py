"""按 GP6 标准 A（法则2 本项 marker 面板）的单口径，切出 6 个谱系的亚聚类细胞清单。

口径出处
--------
PARAMETERS_AND_SOURCES.md §M3-A.5「GP8a 的输入口径」= **标准 A 单口径**。
即：GP6 冻结的 `A_frozen` 标签（seed0）为「某谱系」的全部细胞。
为什么用单口径不用交集，见 §M3-A.5 的两行论证（领域惯例 + 交集的代价）。

本脚本**只切清单，不产标签**。清单哈希落盘，供下游 R 脚本逐字复现输入。

用法:
    python3 05_annotation/00_build_lineage_subsets.py
"""

import hashlib
import json
import os

import pandas as pd

ROOT = "/home/eto/luad_v2"
LABELS = os.path.join(ROOT, "results/05_annotation/gp6_cell_labels.csv.gz")
OUTDIR = os.path.join(ROOT, "results/05_annotation")

# 谱系 → 清单文件名。上皮那份**已存在且在 PARAMETERS 登记了 sha256**，
# 本脚本必须复现出**逐字节相同**的文件，否则硬报错（不许悄悄换掉已签字的输入）。
LINEAGES = {
    "上皮": "epiA_subset_barcodes.txt",
    "T/NK": "tnkA_subset_barcodes.txt",
    "B/浆": "bplasmaA_subset_barcodes.txt",
    "髓系": "myeloidA_subset_barcodes.txt",
    "成纤维": "fibroA_subset_barcodes.txt",
    "内皮": "endoA_subset_barcodes.txt",
}

# 已签字的产物哈希（PARAMETERS §M3-A.5 登记）。改动即报错。
# ⚠️ 2026-09-17 修：原只冻结了 epiA 一份，其余五份**会被静默覆盖**——而它们已经是
#   六次亚聚类 run 的实际输入（各 run_manifest.json 的 cells_file_sha256 即此值）。
#   只冻一份，等于对另外五个谱系开着一扇「悄悄换输入」的后门。现补齐六份。
FROZEN_SHA256 = {
    "epiA_subset_barcodes.txt":
        "a952857e3dbb119d2a5eca6201ba8f7bbe35443d09c695f7e22e08d1b0324d34",
    "tnkA_subset_barcodes.txt":
        "b80e212a1b1e1286db9faf5056b0124674395e270fc8ab550549363bdf452bd0",
    "bplasmaA_subset_barcodes.txt":
        "b3f3ef2b07fb43bd0f0fc079d04329bc8272ab7291bcd59324ed7d073a040d04",
    "myeloidA_subset_barcodes.txt":
        "b6936588925448496afd79cd33a2fedab223cefb9e902410549d98a571e5e312",
    "fibroA_subset_barcodes.txt":
        "505b82a900745cb443dd31e6ab8fdd1c544cf80d7b126b22543c03ec3606f516",
    "endoA_subset_barcodes.txt":
        "1416797a02b41c5efa5491a017e9c49506012a2870ae80585e6d3a4d69d566c7",
}

MIN_CELLS = 1000  # §M3-A.5「子集规模下限」：低于此不做亚聚类


def main():
    d = pd.read_csv(LABELS, usecols=["cell_barcode", "A_frozen"])
    if d["cell_barcode"].duplicated().any():
        raise SystemExit("gp6_cell_labels 里 barcode 有重复")

    unknown = set(d["A_frozen"]) - set(LINEAGES)
    if unknown:
        raise SystemExit(f"A_frozen 出现未登记的标签：{unknown}")

    report = {}
    for lin, fname in LINEAGES.items():
        bcs = d.loc[d["A_frozen"] == lin, "cell_barcode"].tolist()
        path = os.path.join(OUTDIR, fname)
        text = "\n".join(bcs) + "\n"
        digest = hashlib.sha256(text.encode()).hexdigest()

        if fname in FROZEN_SHA256:
            old = open(path).read()
            if hashlib.sha256(old.encode()).hexdigest() != FROZEN_SHA256[fname]:
                raise SystemExit(f"🔴 {fname} 的现有内容与已签字哈希不符，拒绝覆盖")
            if old != text:
                raise SystemExit(f"🔴 {fname} 重算结果与已签字文件**不逐字节相同**，停")
            print(f"  {lin:6s} {len(bcs):7,d} 核  {fname}  ✅ 与已签字文件逐字节一致")
        else:
            with open(path, "w") as fh:
                fh.write(text)
            print(f"  {lin:6s} {len(bcs):7,d} 核  {fname}  sha256 {digest[:16]}…")

        if len(bcs) < MIN_CELLS:
            print(f"    ⚠️ {lin} 低于子集规模下限 {MIN_CELLS}，按 §M3-A.5 不做亚聚类")

        report[lin] = dict(file=fname, n_cells=len(bcs), sha256=digest,
                           subclusterable=len(bcs) >= MIN_CELLS)

    tot = sum(v["n_cells"] for v in report.values())
    if tot != len(d):
        raise SystemExit(f"🔴 各谱系合计 {tot} ≠ 总细胞 {len(d)}（有细胞未归属）")
    print(f"\n合计 {tot:,} = 全部细胞 {len(d):,} ✅（互斥且完备）")

    out = os.path.join(OUTDIR, "lineage_subsets_manifest.json")
    with open(out, "w") as fh:
        json.dump(dict(source=os.path.relpath(LABELS, ROOT),
                       caliber="标准 A 单口径（A_frozen, seed0, r*=0.6）",
                       min_cells=MIN_CELLS, lineages=report),
                  fh, ensure_ascii=False, indent=2)
    print(f"写出 {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
