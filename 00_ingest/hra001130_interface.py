"""
hra001130_interface.py  —  HRA001130 (AAH 全细胞 scRNA) 访问接口（预留）

背景
----
本研究单细胞层面的 **AAH** 阶段缺数据：
  - GSE131907 / GSE189357 / GSE148071 三个 scRNA 队列**均无 AAH**；
  - 公开库中，含 AAH 的**全细胞 scRNA** 唯一来源是 **HRA001130**（GSA-Human，**受控，需申请**）。

HRA001130 出处（Nat Commun 2021）
  "Deciphering cell lineage specification of human lung adenocarcinoma with single-cell
   RNA sequencing".  DOI: 10.1038/s41467-021-26770-2
  数据：GSA-Human `HRA001130`，scRNA-seq（10x），含 AAH(n=3)/AIS(n=4)/MIA(n=5)/IA + 18 配对正常肺。
  申请入口：https://ngdc.cncb.ac.cn/gsa-human/browse/HRA001130

设计
----
本模块只定义**接口与期望布局**，不含任何数据、不做任何伪造：
  - 数据未到位 → `load_hra001130()` 返回 None 并打印获取指引（**不报错、不造数据**）；
  - 数据到位   → 按期望布局读取，返回统一的 per-cell 表（barcode / patient_id / sample_id / stage / histology）。

在等审批期间，AAH 暂用 **GSE308103（snRNA）** 经 **SCMG zero-shot 跨平台并入**
（见 `cohort_registry.py` 与白皮书 §3.1 的跨模态五判据）。
获批后：把数据放到 `DATA_DIR` 期望路径，本接口即自动生效，**无需改动下游**。
"""

from __future__ import annotations

import os
from typing import Optional

# 期望数据落点（获批后把解压内容放这里）
DATA_DIR = os.environ.get("HRA001130_DIR", "/home/eto/luad_v2/data/HRA001130")

# 期望的样本分期（来自论文；用于**校验**，不是硬编码结果）
EXPECTED_STAGE_COUNTS = {"AAH": 3, "AIS": 4, "MIA": 5, "IA": "≥1", "Normal": 18}

# 访问与出处
ACCESS_URL = "https://ngdc.cncb.ac.cn/gsa-human/browse/HRA001130"
PAPER = "Deciphering cell lineage specification of human LUAD with scRNA-seq, Nat Commun 2021, DOI 10.1038/s41467-021-26770-2"


def _looks_available(d: str) -> bool:
    """判断数据是否已就位（存在 barcodes/矩阵/元数据任一常见形态）。"""
    if not os.path.isdir(d):
        return False
    marks = ("barcodes.tsv", "matrix.mtx", "filtered_feature_bc_matrix",
             "metadata.csv", "cell_metadata.csv")
    return any(os.path.exists(os.path.join(d, m)) for m in marks)


def load_hra001130(data_dir: str = DATA_DIR) -> Optional[dict]:
    """返回统一的 per-cell 表（dict/DataFrame）；数据未到位时返回 None。

    返回字段（统一口径，供下游一致使用）：
        cell_barcode, patient_id, sample_id, stage, histology, dataset='HRA001130'
    """
    if not _looks_available(data_dir):
        print(
            "[HRA001130] 数据未就位（受控库，需申请）。\n"
            f"  预期路径 : {data_dir}\n"
            f"  申请入口 : {ACCESS_URL}\n"
            f"  论文      : {PAPER}\n"
            "  → 暂时使用 GSE308103(snRNA) 经 SCMG 跨平台并入 AAH；本接口留待获批后启用。"
        )
        return None

    # 数据到位后在此实现读取（不在此处编造任何内容）：
    #   1) 读矩阵 + 元数据；2) 映射到统一字段；3) 用 EXPECTED_STAGE_COUNTS 做**校验**
    #      （数量/分期对不上就 raise，不静默通过）。
    raise NotImplementedError(
        "HRA001130 已检测到数据，但读取实现待补（获批后按实际文件结构实现，"
        "并以 EXPECTED_STAGE_COUNTS 做校验）。"
    )


def status() -> dict:
    """供 registry/文档引用：当前可用性。"""
    return {
        "dataset": "HRA001130",
        "modality": "scRNA (whole-cell)",
        "role": "AAH 单细胞（首选，待申请）",
        "available": _looks_available(DATA_DIR),
        "path": DATA_DIR,
        "access_url": ACCESS_URL,
        "paper": PAPER,
    }


if __name__ == "__main__":
    print(status())
    load_hra001130()
