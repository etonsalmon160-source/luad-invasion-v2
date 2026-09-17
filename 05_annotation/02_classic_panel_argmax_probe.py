"""经典 marker 面板在**已冻结的 GP6 簇**上的 argmax 探针 —— 只做诊断，不产新标签。

要回答的问题
------------
我先前判断「经典面板做 score_genes→argmax 会被 AT1/AT2 系统性吞掉」。
该判断**在机制上站不住**：`score_genes` 的对照基因是按表达量分箱后**从同一箱**抽的
（ctrl_size=50, n_bins=25），得分对任何 marker 集都中心化在 0，低检出集不处于基线劣势。
⇒ 不能靠推理定论，**直接测**。

本脚本在 **GP6 已冻结的簇**（seed0，不重新聚类、不调参）上，用经典面板算
逐细胞得分 → 簇均值 → argmax，看标签是否塌到 AT1/AT2。

⚠️ 本脚本**不产出任何用于下游的标签**；只是方法学探针。GP8a 的真正标签须另经签字。

用法:
    python3 05_annotation/02_classic_panel_argmax_probe.py
"""

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import epi_classic_panel as C  # noqa: E402

H5 = "results/02_expression/gse308103_counts_paperqc.h5ad"
CLU = "results/04_integration/seurat_trad/full/clusters.csv.gz"
EPI = "results/05_annotation/epiA_subset_barcodes.txt"
SEED = 0
CTRL_SIZE = 50
SCORE_SEED = 0
OUT = "results/05_annotation/epi_classic_argmax_probe.csv"


def main():
    epi = [ln for ln in open(EPI).read().split("\n") if ln]
    print(f"上皮 barcode {len(epi):,} 个")

    adata = sc.read_h5ad(H5)
    print(f"  shape={adata.shape} nnz={adata.X.nnz:,}")

    clu = pd.read_csv(CLU)
    assert len(clu) == adata.n_obs, (len(clu), adata.n_obs)
    assert (clu["cell_barcode"].to_numpy() == adata.obs_names.to_numpy()).all(), "分群表顺序不一致"
    adata.obs["clu"] = pd.Categorical(clu[f"harmony_res0.6_seed{SEED}"].astype(str))

    # 归一化与 GP6 完全一致
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 上皮子集（GP8a 的实际输入范围）
    keep = adata.obs_names.isin(set(epi))
    print(f"  命中上皮 {int(keep.sum()):,} / {len(epi):,}")
    sub = adata[keep].copy()
    del adata
    print(f"  子集 shape={sub.shape}")

    _, missing = C.build_panel(sub.var_names)
    panel = C.build_panel(sub.var_names)[0]
    print(f"面板缺失: {missing}")

    for st, gl in panel.items():
        if len(gl) < 2:
            print(f"  ⚠️ {st} 可用基因仅 {len(gl)}，跳过")
            continue
        sc.tl.score_genes(sub, gl, ctrl_size=CTRL_SIZE, random_state=SCORE_SEED,
                          score_name=f"s_{st}")

    cols = [c for c in sub.obs.columns if c.startswith("s_")]
    types = [c[2:] for c in cols]
    M = sub.obs[cols].to_numpy(dtype=np.float64)
    cl = sub.obs["clu"].astype(str).to_numpy()
    ucl = sorted(set(cl), key=lambda x: int(x) if x.isdigit() else x)

    # 簇均值
    cm = np.vstack([M[cl == c].mean(axis=0) for c in ucl])
    win = np.array(types)[cm.argmax(axis=1)]
    ncell = np.array([int((cl == c).sum()) for c in ucl])

    df = pd.DataFrame(cm, columns=types)
    df.insert(0, "cluster", ucl)
    df.insert(1, "n_cells", ncell)
    df.insert(2, "argmax", win)
    df.insert(3, "win_margin", cm.max(axis=1) - np.sort(cm, axis=1)[:, -2])
    df.to_csv(OUT, index=False)

    print("\n=== 每个 GP6 簇（上皮子集内）的经典面板 argmax ===")
    print(df[["cluster", "n_cells", "argmax", "win_margin"]].to_string(index=False))

    print("\n=== argmax 胜出的亚型分布（我原判断：会塌到 AT1/AT2）===")
    vc = pd.Series(win).value_counts()
    for k, v in vc.items():
        print(f"  {k:18s} {v:2d} 簇 / {ncell[win == k].sum():7,d} 细胞")
    never = [t for t in types if t not in set(win)]
    print(f"\n从未胜出的亚型: {never}")
    print(f"\n写成 {OUT}")


if __name__ == "__main__":
    main()
