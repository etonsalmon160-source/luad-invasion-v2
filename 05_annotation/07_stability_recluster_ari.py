"""GP8a 稳定性校验：**剔掉跨谱系污染簇之后，聚类结构是否还站得住**。

背景（用户 2026-09-21 裁定）
----------------------------
`06_contamination_check.py` 在 141,105 个上皮核里标出 2 个非上皮簇（8,142 核 = 5.77%）。
用户选的处置是「**原地注释 + 并行跑稳定性校验**」：
  - 主分析仍用**已签字的那份 141,105 核子集**，只是注释时把被标簇剔出统计；
  - 同时**另跑一次**：只喂剩下 132,963 核重做全套聚类，看簇结构变不变。

本脚本只做**比较**，不重算聚类。判读口径**不新立阈值**：把实测 ARI 摆在
项目**自己签过的**两个参照量旁边 ——

  - `ari_seed_mean`（r=0.5 处 5 个种子的两两 ARI 均值）= **同一协议、只换随机种子**的漂移量
  - `ari_xres`（与相邻分辨率的 ARI）= 协议里**真改了参数**时的漂移量

⚠️ 诚实说明：`≥0.90` 那条阈值当初是给**跨种子**稳定性签的，**不是**给"换输入细胞集"签的。
本脚本因此**不套用**它下结论，只把三个数并排给出，由人工裁决（与用户「人工定」的口径一致）。

用法:
    python3 05_annotation/07_stability_recluster_ari.py --tag epiA --clean-tag epiA_clean \
        --lineage 上皮 --res 0.5 --seed 0
"""

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score as ari

ROOT = "/home/eto/luad_v2"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ROOT}/results/05_annotation"

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True, help="主分析 tag（已签字子集上的原地注释）")
    ap.add_argument("--clean-tag", required=True, help="干净子集重聚类 tag")
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--res", type=float, required=True, help="要比对的分辨率（用签字的 r*）")
    ap.add_argument("--seed", type=int, required=True)
    a = ap.parse_args()

    clean_bc = f"{OUT}/{a.tag}_nocontam_subset_barcodes.txt"
    if not os.path.exists(clean_bc):
        raise SystemExit(f"🔴 找不到干净 barcode 清单 {clean_bc} —— 须先用 --contamination-csv 跑过注释")
    bcs = [ln for ln in open(clean_bc).read().split("\n") if ln]
    log(f"干净子集清单 {len(bcs):,} 核")

    p_main = f"{TRAD}/{a.tag}/clusters.csv.gz"
    p_clean = f"{TRAD}/{a.clean_tag}/clusters.csv.gz"
    if not os.path.exists(p_clean):
        raise SystemExit(f"🔴 干净集重聚类还没跑：{p_clean} 不存在")
    m = pd.read_csv(p_main)
    c = pd.read_csv(p_clean)
    log(f"主分析 {len(m):,} 核 / 重聚类 {len(c):,} 核")

    if len(c) != len(bcs):
        raise SystemExit(f"🔴 重聚类细胞数 {len(c)} ≠ 清单 {len(bcs)}（输入没对齐，拒绝比较）")
    if not set(c["cell_barcode"]) == set(bcs):
        raise SystemExit("🔴 重聚类的 barcode 集合与清单不一致")

    m = m.set_index("cell_barcode")
    c = c.set_index("cell_barcode")
    bcs = [b for b in c.index if b in m.index]        # 只在两边都有的细胞上比
    if len(bcs) != len(c):
        raise SystemExit(f"🔴 主分析缺 {len(c) - len(bcs)} 个干净集细胞（主分析必须是干净集的超集）")

    ref = pd.read_csv(f"{TRAD}/{a.tag}/resolution_metrics.csv")
    ref = ref[ref["resolution"] == a.res]
    if not len(ref):
        raise SystemExit(f"🔴 参照表里没有 resolution={a.res}")
    ref = ref.iloc[0]

    rows = []
    for col in sorted([x for x in c.columns if x.startswith("harmony_res")]):
        if col not in m.columns:
            log(f"⚠️ 重聚类有而主分析没有的列，跳过：{col}")
            continue
        rows.append(dict(
            cluster_col=col,
            n_clusters_main=int(m.loc[bcs, col].nunique()),
            n_clusters_clean=int(c.loc[bcs, col].nunique()),
            ari=round(float(ari(m.loc[bcs, col].to_numpy(), c.loc[bcs, col].to_numpy())), 4),
            is_rstar=("res%.1f" % a.res in col.replace("harmony_", "")),
        ))
    d = pd.DataFrame(rows).sort_values("cluster_col")
    out = f"{OUT}/{a.tag}_vs_{a.clean_tag}_ari.csv"
    d.to_csv(out, index=False)
    log(f"写出 {os.path.relpath(out, ROOT)}")

    col = f"harmony_res{a.res}_seed{a.seed}"
    hit = d[d["cluster_col"] == col]
    if not len(hit):
        raise SystemExit(f"🔴 两个 run 都没有 {col}，无法比较")
    hit = hit.iloc[0]

    print(f"\n=== {a.tag} 原地注释 vs {a.clean_tag} 干净集重聚类（只在 {len(bcs):,} 个干净核上比）===")
    print(d.to_string(index=False))
    print(f"\n签字分辨率 r*={a.res} 处（{col}）：")
    print(f"  ARI(原地, 干净集重聚类)      = {hit['ari']:.4f}")
    print(f"  参照一 ari_seed_mean（只换种子）= {ref['ari_seed_mean']:.4f}  ← 同一协议下的正常漂移")
    print(f"  参照二 ari_xres（换了分辨率）  = {ref['ari_xres']:.4f}  ← 真改了参数的漂移")
    print(f"  簇数 原地 {hit['n_clusters_main']} → 干净集 {hit['n_clusters_clean']}（剔掉 2 个污染簇前，"
          f"主分析共 {int(ref['n_clusters'])} 个）")
    print("\n⚠️ 判读**由人工做**：本脚本不套用 ≥0.90（那条阈值当初是给跨种子签的，"
          "不是给换输入细胞集签的）。")

    man = dict(
        main_tag=a.tag, clean_tag=a.clean_tag, lineage=a.lineage,
        res=a.res, seed=a.seed, cluster_col=col,
        n_cells_compared=len(bcs),
        ari_at_rstar=float(hit["ari"]),
        reference=dict(ari_seed_mean=float(ref["ari_seed_mean"]),
                       ari_xres=float(ref["ari_xres"]),
                       n_clusters_main_run=int(ref["n_clusters"]),
                       source="主分析 resolution_metrics.csv，同一 r* 行，项目已签"),
        no_threshold_applied=True,
        threshold_note="≥0.90 是跨种子口径的签字阈值；跨协议比较**未**套用，只并排给出参照量由人工裁决",
        inputs={os.path.relpath(p, ROOT): sha256(p) for p in
                [p_main, p_clean, clean_bc, f"{TRAD}/{a.tag}/resolution_metrics.csv",
                 f"{os.path.abspath(__file__)}"]},
        outputs={os.path.relpath(out, ROOT): sha256(out)},
    )
    out_man = f"{OUT}/{a.tag}_stability_manifest.json"
    json.dump(man, open(out_man, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log(f"写出 {os.path.relpath(out_man, ROOT)}")


if __name__ == "__main__":
    main()
