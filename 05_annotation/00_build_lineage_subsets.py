"""按 GP6 标准 A（法则2 本项 marker 面板）的单口径，切出 6 个谱系的亚聚类细胞清单。

口径出处
--------
PARAMETERS_AND_SOURCES.md §M3-A.5「GP8a 的输入口径」= **标准 A 单口径**。
即：GP6 冻结的 `A_frozen` 标签（seed0）为「某谱系」的全部细胞。
为什么用单口径不用交集，见 §M3-A.5 的两行论证（领域惯例 + 交集的代价）。

2026-09-22 起，口径改为 **`A_adjudicated`**（`A_frozen` + 全局簇错标裁决，见
`09_adjudicate_global_labels.py`）。裁决只动了 4 个全局簇、共 11,228 核：
raw 17 上皮→髓系、raw 26 B/浆→髓系、raw 38 上皮→T/NK、raw 44 内皮→T/NK。
`A_frozen` 原样保留在裁决版文件里，新旧两版清单**都可复现**（见 `SOURCES`）。

本脚本**只切清单，不产标签**。清单哈希落盘，供下游 R 脚本逐字复现输入。

用法:
    python3 05_annotation/00_build_lineage_subsets.py
"""

import hashlib
import json
import os

import pandas as pd

ROOT = "/home/eto/luad_v2"
OUTDIR = os.path.join(ROOT, "results/05_annotation")

# ---- 标签来源：**按输入自身的哈希**选锁定表 -------------------------------------
# 2026-09-22 改。原来只有一份 `gp6_cell_labels.csv.gz` + 一张 FROZEN_SHA256 表。
# 现在多了一份裁决版（`gp6_cell_labels_adjudicated.csv.gz`，见 09_adjudicate_global_labels.py），
# 于是把两张表都钉住、由**输入文件的 sha256** 决定用哪张 ——
# 这样「新旧两版各自的六份清单」都可复现，而任何**未登记**的标签文件一律硬报错，
# 比原来只钉一版更严（原来换成任何别的文件都看不出来）。
SOURCES = {
    "gp6_cell_labels.csv.gz": dict(
        label_col="A_frozen",
        source_sha256="fa031a52a4ab3b7098d3fca5555b0749b2f4622711153c0b742f8ae09aaf12e2",
        caliber="标准 A 单口径（A_frozen, seed0, r*=0.6）",
        requires_signature=None,
        files={
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
        }),
    "gp6_cell_labels_adjudicated.csv.gz": dict(
        label_col="A_adjudicated",
        source_sha256="dfaecc05c529bbe7cfb06990fd5cd13c99de03adf906fd1e534ded7e6c743093",
        caliber="标准 A 单口径 + 2026-09-22 全局簇错标裁决（A_adjudicated, seed0, r*=0.6）",
        # 裁决版是**人工裁决**的产物，必须带签字才准用（同 r* 的闸门思路）。
        requires_signature="gp6_adjudication_manifest.json",
        files={
            "epiA_subset_barcodes.txt":
                "d12a115a87d0790d7fe8cf27b16170999d94daba51f3016fdb79c5aa3fdde55e",
            "tnkA_subset_barcodes.txt":
                "8cb87515d5463bdbdb9c6885decef19611a73b5a8f6abbdc850507d5c620b0c5",
            "bplasmaA_subset_barcodes.txt":
                "7e34c5c5baedb0ecfef071a8e2f3c738964c97cc29d74aa13c65c2a731a5cd0a",
            "myeloidA_subset_barcodes.txt":
                "2fd34840c92f8a6ebf11ea8713903af7ecd7d0ed5a6622864836f01805b2fb1f",
            "fibroA_subset_barcodes.txt":
                "505b82a900745cb443dd31e6ab8fdd1c544cf80d7b126b22543c03ec3606f516",
            "endoA_subset_barcodes.txt":
                "9ee8d0c08a26bd14592b654d3c9578561d1143f5bb6f2a74e544e8306f2396ac",
        }),
}

# 当前口径：裁决版。改这一行就切换回旧版（旧版六份清单同样可复现）。
LABELS_BASENAME = "gp6_cell_labels_adjudicated.csv.gz"

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

MIN_CELLS = 1000  # §M3-A.5「子集规模下限」：低于此不做亚聚类


def main():
    if LABELS_BASENAME not in SOURCES:
        raise SystemExit(f"🔴 {LABELS_BASENAME} 未在 SOURCES 登记 —— 拒绝用未登记的标签切子集")
    src = SOURCES[LABELS_BASENAME]
    labels = os.path.join(OUTDIR, LABELS_BASENAME)
    frozen = src["files"]
    col = src["label_col"]

    actual = hashlib.sha256(open(labels, "rb").read()).hexdigest()
    if actual != src["source_sha256"]:
        raise SystemExit(f"🔴 {LABELS_BASENAME} 的 sha256 与登记值不符：\n"
                         f"    实际 {actual}\n    登记 {src['source_sha256']}")
    print(f"标签来源 {LABELS_BASENAME}（列 {col}）  sha256 ✅")

    if src["requires_signature"]:
        sp = os.path.join(OUTDIR, src["requires_signature"])
        if not os.path.exists(sp):
            raise SystemExit(f"🔴 缺 {src['requires_signature']} —— 裁决版标签须先签字")
        j = json.load(open(sp, encoding="utf-8"))
        if not j.get("signed_by"):
            raise SystemExit(f"🔴 {src['requires_signature']} 尚无 `signed_by` —— "
                             f"拒绝用未签字的裁决切子集")
        print(f"签字闸门 ✅ {src['requires_signature']}：signed_by = {j['signed_by']}")

    d = pd.read_csv(labels, usecols=["cell_barcode", col])
    if d["cell_barcode"].duplicated().any():
        raise SystemExit("gp6_cell_labels 里 barcode 有重复")

    unknown = set(d[col]) - set(LINEAGES)
    if unknown:
        raise SystemExit(f"{col} 出现未登记的标签：{unknown}")

    report = {}
    for lin, fname in LINEAGES.items():
        bcs = d.loc[d[col] == lin, "cell_barcode"].tolist()
        path = os.path.join(OUTDIR, fname)
        text = "\n".join(bcs) + "\n"
        digest = hashlib.sha256(text.encode()).hexdigest()

        if os.path.exists(path):
            old = open(path).read()
            if hashlib.sha256(old.encode()).hexdigest() != frozen[fname]:
                raise SystemExit(f"🔴 {fname} 的现有内容与登记哈希不符，拒绝覆盖")
            if old != text:
                raise SystemExit(f"🔴 {fname} 重算结果与登记文件**不逐字节相同**，停")
            print(f"  {lin:6s} {len(bcs):7,d} 核  {fname}  ✅ 与登记文件逐字节一致")
        else:
            if digest != frozen[fname]:
                raise SystemExit(f"🔴 {fname} 重算 sha256 与登记值不符：\n"
                                 f"    实际 {digest}\n    登记 {frozen[fname]}")
            with open(path, "w") as fh:
                fh.write(text)
            print(f"  {lin:6s} {len(bcs):7,d} 核  {fname}  ✅ 新建，sha256 与登记值一致")

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
        json.dump(dict(source=os.path.relpath(labels, ROOT),
                       source_sha256=actual,
                       label_column=col,
                       # ⚠️ 2026-09-17 补：原 manifest 只记了输入**路径**不记**哈希** —— 于是
                       # 「六份清单哈希」这条溯源链在输入端是断的（换了另一版标签文件也看不出来）。
                       caliber=src["caliber"],
                       supersedes_note=(
                           "2026-09-22：本清单由 A_adjudicated 切出，取代 "
                           "results/_superseded/2026-09-22_cluster_adjudication/ 下的旧清单。"
                           "改变：上皮 141,105→133,384、髓系 54,407→64,084、B/浆 34,385→31,282、"
                           "T/NK 62,483→64,034、内皮 43,684→43,280、成纤维不变。"),
                       min_cells=MIN_CELLS, lineages=report),
                  fh, ensure_ascii=False, indent=2)
    print(f"写出 {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
