"""GP8a 前置：上皮子集**跨谱系污染**体检（只诊断，不改数据）。

为什么需要这一步
----------------
2026-09-21 GP8a 上皮注释首次实跑，逐簇 top 富集基因暴露出子集里混着别的细胞：
eg 簇 20 的 top 基因是 `TRAC/CD2/ITGAL/TRBC2`（T 细胞），簇 6 的是 `CPA3/HDAC/KIT/MS4A2`
（肥大/嗜碱）。而上皮面板**没有"不是上皮"这个选项** ⇒ 纯 argmax 把这类簇也硬贴上 AT2。
⇒ 子集定义（GP6 标准 A 的**逐细胞** argmax，噪声大）需要一道**簇级**复核。

判据（**先定规则，后看结果**；本文件写完即冻结，不因结果调参）
--------------------------------------------------------------
两条**互相独立**的证据，命中任一即判该簇"非上皮，须剔除"：

  ① **面板线**（与 GP6 标准 A 同源、可比）：用 `marker_panel.py` 法则2 六谱系面板逐细胞打分，
     取**簇内均值**，argmax ≠ 上皮 ⇒ 非上皮。
     — 为什么用簇级而非逐细胞：逐细胞 argmax 是 GP6 已用的口径，噪声正来自它；
       簇级平均把若干细胞的随机抖动抵消掉。同一套基因、同一套 score_genes 参数，**只是聚合层级不同**。
  ② **CellTypist 线**（独立图谱，与标准 A 不同源）：GP6 已存的 `B_lineage`，簇内**非上皮占比 > 50%**
     ⇒ 非上皮。"多数细胞被独立方法判为别的谱系"是自然的多数据，不为迁就任何具体簇而调。

⚠️ 两条线**可能矛盾**（本批实测确有）。矛盾时的处理**自动完成、不人工挑**：
   - 任一命中 ⇒ 标 `contaminated`，**但 `flag` 列写明是哪条线命中的**，报告里分别陈述。
   - 面板线命中、CellTypist 未命中 ⇒ 多半是**环境 RNA**（游离 mRNA 吸附，snRNA 通病），
     须与"真的混进别的细胞"区分处理，**不得合并报**。
   - 细胞占比：`外来基因` 类（免疫球蛋白 `IGKC/IGHG*`、`COL1A1` 等）在多个簇普遍轻量出现，
     属环境 RNA 底噪；**只有当一个簇被单一外来程序整体占据**才算混入。

本脚本**只读**：不写任何清单、不改 signed 的子集文件。
**不做**：不决定"剔除后要不要重聚类"（那是检查点问题，须另行裁定）。

用法:
    python3 05_annotation/06_contamination_check.py --tag epiA --lineage 上皮 --res 0.5 --seed 0
"""

import argparse
import hashlib
import importlib
import json
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import classic_panels as CP          # noqa: E402  上皮亚型面板（判断"像哪一种上皮"）
import marker_panel as MP            # noqa: E402  法则2 六谱系面板（判断"是不是上皮"）

ROOT = "/home/eto/luad_v2"
H5 = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ROOT}/results/05_annotation"
LABELS = f"{OUT}/gp6_cell_labels.csv.gz"
SUBSETS = f"{OUT}/lineage_subsets_manifest.json"

CTRL_SIZE = 50
SCORE_SEED = 0
NONEPI_MAJORITY = 0.50   # CellTypist 线阈值：非上皮占比过半

# CellTypist 的"未归属"不计入非上皮（它是"没把握"，不是"判成别的谱系"）
CT_NA = {"未归属", None, "", "nan"}

# 合并口径按谱系分（2026-09-21 用户裁定，理由见 manifest 的 rule.combine_reason）
#   上皮   —— 沿用 GP8a 已签字的原口径："命中任一即剔"
#   非上皮 —— **须两条线同时命中**才剔；单线命中如实留痕，但不自动剔
# ⚠️ 这是**指出规则的作用域**，不是事后放宽阈值（法则 3.2）：上皮口径一字未动，
#    CellTypist 线的 0.50 与面板线的 ctrl_size/seed 也一字未动。
COMBINE_ANY_LINE = {"上皮"}

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--res", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    a = ap.parse_args()

    if a.lineage not in MP.PANEL:
        raise SystemExit(f"六谱系面板里没有 {a.lineage}；有 {MP.LINEAGES}")

    # ---- 输入（全部记哈希，可复现） ----
    man = json.load(open(SUBSETS, encoding="utf-8"))
    if a.lineage not in man["lineages"]:
        raise SystemExit(f"子集清单里没有 {a.lineage}")
    bc_file = f"{OUT}/{man['lineages'][a.lineage]['file']}"
    bcs = [ln for ln in open(bc_file).read().split("\n") if ln]
    log(f"{a.lineage} 子集 {len(bcs):,} 核（清单 {os.path.basename(bc_file)}）")

    # ⚠️ 种子必须与注释步骤（03_subcluster_annotation.py）**是同一个**。
    #    本表是按**簇编号**交给 03 去剔细胞的；用错种子 ⇒ 簇编号对不上 ⇒ 剔错细胞，
    #    而且全程静默。2026-09-21 实测踩到：髓系/B浆/成纤维 的注释须用种子 1/1/4
    #    （05_seed_representativeness.py 的预注册判据），而本体检最初一律按种子 0 跑。
    # ⚠️ 分辨率同理会踩同样的坑：体检必须在**已签字**的 r* 上做，否则簇编号同样对不上
    #    03。2026-09-21 补装。
    #    两道闸门都**直接复用 03 的**，不自己另立一套，免得两边漂移。
    AN = importlib.import_module("03_subcluster_annotation")
    rstar_ev = AN.check_rstar_signature(a.tag, a.res)
    a.seed, seed_ev = AN.enforce_annotation_seed(a.tag, a.seed)

    clu = pd.read_csv(f"{TRAD}/{a.tag}/clusters.csv.gz")
    col = f"harmony_res{a.res}_seed{a.seed}"
    if col not in clu.columns:
        raise SystemExit(f"集群表没有 {col}")
    assert len(clu) == len(bcs), (len(clu), len(bcs))

    # ---- 表达对象 ----
    big = ad.read_h5ad(H5)
    keep = big.obs_names.isin(set(bcs))
    if int(keep.sum()) != len(bcs):
        raise SystemExit(f"🔴 命中 {int(keep.sum())} ≠ 清单 {len(bcs)}")
    sub = big[keep].copy()
    del big

    if not (clu["cell_barcode"].to_numpy() == sub.obs_names.to_numpy()).all():
        log("⚠️ 集群表与子集顺序不一致，按 barcode 重排")
        clu = clu.set_index("cell_barcode").loc[sub.obs_names].reset_index()

    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    # ---- ① 面板线：法则2 六谱系，逐细胞打分 → 簇内均值 ----
    have = set(sub.var_names)
    used, miss = {}, {}
    for lin in MP.LINEAGES:
        g = [x for x in MP.PANEL[lin]["genes"] if x in have]
        lost = [x for x in MP.PANEL[lin]["genes"] if x not in have]
        if lost:
            miss[lin] = lost
        if len(g) < 2:
            raise SystemExit(f"🔴 谱系 {lin} 在矩阵里可用基因 <2：{g}")
        used[lin] = g
    if miss:
        log(f"⚠️ 六谱系面板在本矩阵缺基因（不影响本判据，但须随报告说明）：{miss}")

    for lin in MP.LINEAGES:
        sc.tl.score_genes(sub, used[lin], ctrl_size=CTRL_SIZE, random_state=SCORE_SEED,
                          score_name=f"L_{lin}")

    sub.obs["clu"] = pd.Categorical(clu.set_index("cell_barcode")
                                    .loc[sub.obs_names, col].astype(str).to_numpy())
    cl = sub.obs["clu"].astype(str).to_numpy()
    ucl = sorted(set(cl), key=lambda x: (len(x), x))
    lcols = [f"L_{lin}" for lin in MP.LINEAGES]

    M = pd.DataFrame(sub.obs[lcols].to_numpy(dtype=np.float64), columns=MP.LINEAGES)
    M["clu"] = cl
    per_clu = M.groupby("clu")[MP.LINEAGES].mean()
    per_clu["n_cells"] = M.groupby("clu").size()
    log(f"面板线：{len(ucl)} 簇 × {len(MP.LINEAGES)} 谱系打分完毕")

    # ---- ② CellTypist 线：GP6 已存的 B_lineage ----
    lab = pd.read_csv(LABELS, usecols=["cell_barcode", "A_frozen", "B_lineage"])
    m = clu[["cell_barcode"]].merge(lab, on="cell_barcode", how="left")
    if m["B_lineage"].isna().any():
        raise SystemExit(f"🔴 有 {int(m['B_lineage'].isna().sum())} 个细胞在 GP6 标签表里查不到")
    m["clu"] = cl
    m["B_非上皮"] = ~m["B_lineage"].isin(CT_NA | {a.lineage})
    ct_frac = m.groupby("clu")["B_非上皮"].mean()
    ct_main = m.groupby("clu")["B_lineage"].agg(lambda s: s.value_counts().index[0])
    ct_mainfrac = m.groupby("clu")["B_lineage"].agg(lambda s: s.value_counts().iloc[0] / len(s))

    # ---- 逐簇裁定 ----
    need_both = a.lineage not in COMBINE_ANY_LINE
    rows = []
    for c in ucl:
        sc_row = per_clu.loc[c, MP.LINEAGES]
        panel_call = str(sc_row.idxmax())
        srt = np.sort(sc_row.to_numpy())
        margin = float(srt[-1] - srt[-2])
        nonepi_panel = panel_call != a.lineage
        nonepi_ct = bool(ct_frac.loc[c] > NONEPI_MAJORITY)
        hits = []
        if nonepi_panel:
            hits.append("面板")
        if nonepi_ct:
            hits.append("CellTypist")
        # 本谱系口径：上皮任一即剔；非上皮须两线同时命中才剔
        hit = hits if (len(hits) == 2 or not need_both) else []
        # 单线命中（按本谱系口径**不剔**）单独留痕。
        # ⚠️ 不能塞进 flag —— 03_subcluster_annotation.py 认的是 `flag != "保留"` 即剔。
        # ⚠️ 无命中写 "无" 而非空串：空串落盘成空字段，pandas 读回来是 NaN，
        #    而 `NaN != ""` 恒为真 —— 会把"没留痕"误判成"留痕了"。用字面值堵死这个坑。
        hold = "+".join(hits) if (not hit and hits) else "无"
        rows.append(dict(
            cluster=c, n_cells=int(per_clu.loc[c, "n_cells"]),
            panel_call=panel_call, panel_margin=round(margin, 4),
            panel_hit=bool(nonepi_panel),
            panel_self=round(float(sc_row[a.lineage]), 4),
            panel_top_other=str(sc_row.drop(a.lineage).idxmax()),
            panel_top_other_score=round(float(sc_row.drop(a.lineage).max()), 4),
            ct_hit=bool(nonepi_ct),
            ct_main=ct_main.loc[c], ct_main_frac=round(float(ct_mainfrac.loc[c]), 3),
            ct_nonepi_frac=round(float(ct_frac.loc[c]), 4),
            flag="+".join(hit) if hit else "保留",
            single_line_hold=hold,
            agree=("两条线一致" if len(hits) == 2 else
                   "仅一条线命中" if len(hits) == 1 else "—"),
        ))
    df = pd.DataFrame(rows)
    out_csv = f"{OUT}/{a.tag}_contamination.csv"
    df.to_csv(out_csv, index=False)
    log(f"写出 {os.path.relpath(out_csv, ROOT)}")

    # 逐簇 × 六谱系 明细也留下（供报告引用，避免只给结论）
    det = per_clu[["n_cells"] + MP.LINEAGES].copy()
    det.round(4).to_csv(f"{OUT}/{a.tag}_lineage_scores_by_cluster.csv")
    log(f"写出 {a.tag}_lineage_scores_by_cluster.csv")

    # ---- 报告 ----
    bad = df[df["flag"] != "保留"]
    both = df[df["flag"] == "面板+CellTypist"]
    one = df[(df["flag"] != "保留") & (df["flag"] != "面板+CellTypist")]
    hold = df[df["single_line_hold"] != "无"]
    print(f"\n=== {a.tag} 污染体检（{len(df)} 簇 / {len(bcs):,} 核）===")
    print(df.to_string(index=False))
    print(f"\n合并口径：{'非上皮谱系 —— 须两条线同时命中才剔' if need_both else '上皮 —— 命中任一即标（GP8a 已签字原口径）'}")
    print(f"判定非{a.lineage}、须剔除：{len(bad)} 簇 / {int(bad['n_cells'].sum()):,} 核"
          f"（占子集 {bad['n_cells'].sum() / len(bcs) * 100:.2f}%）")
    print(f"  两条线一致（证据最硬）：{list(both['cluster'])} 共 {int(both['n_cells'].sum()):,} 核")
    print(f"  仅一条线命中（须单独判断）：{list(one['cluster'])} 共 {int(one['n_cells'].sum()):,} 核")
    if len(one):
        print("  ⚠️ 仅面板线命中、CellTypist 未命中 ⇒ 优先怀疑**环境 RNA**而非混入他种细胞；")
        print("     仅 CellTypist 命中 ⇒ 优先怀疑**面板基因槽太少**导致面板线失灵。")
    if len(hold):
        print(f"\n  ⚠️ 单线命中、按本谱系口径**不剔**（如实留痕，交人工判读）："
              f"{len(hold)} 簇 / {int(hold['n_cells'].sum()):,} 核")
        for _, r in hold.iterrows():
            print(f"     簇 {r['cluster']}（{int(r['n_cells']):,} 核）仅 [{r['single_line_hold']}] 命中"
                  f"｜面板自评 {r['panel_self']} / 对面最高 {r['panel_top_other']} {r['panel_top_other_score']}"
                  f"｜CellTypist 主流 {r['ct_main']} {r['ct_main_frac']}")

    out_man = f"{OUT}/{a.tag}_contamination_manifest.json"
    json.dump(dict(
        tag=a.tag, lineage=a.lineage, res=a.res, seed=a.seed, cluster_col=col,
        seed_gate=seed_ev,
        rstar_gate=dict(r_star=rstar_ev.get("r_star"), relaxed=bool(rstar_ev.get("relaxed")),
                        seed_min_applied=rstar_ev.get("seed_min_applied"),
                        signed_by=rstar_ev.get("signed_by"),
                        source=f"{a.tag}_rstar.json"),
        rstar_note=("🔴 本谱系 r* 为**事后放宽**口径 ⇒ 结果只能作探索性结论，不得进主结论"
                    if rstar_ev.get("relaxed") else "r* 为预注册口径签字"),
        n_cells=len(bcs), n_clusters=len(df),
        rule=dict(
            panel_line=f"法则2 六谱系面板逐细胞 score_genes → 簇内均值 → argmax ≠ {a.lineage} ⇒ 非本谱系",
            celltypist_line=f"GP6 已存 B_lineage 的簇内非本谱系占比 > {NONEPI_MAJORITY} ⇒ 非本谱系",
            combine=("命中任一即标" if not need_both else "须两条线同时命中才标",
                     f"（本谱系 {a.lineage} 适用；单线命中写入 single_line_hold 列留痕但不剔）"),
            combine_scope=dict(
                any_line=sorted(COMBINE_ANY_LINE),
                need_both=sorted(set(MP.LINEAGES) - COMBINE_ANY_LINE),
                applied_to=a.lineage,
                applied="命中任一即标" if not need_both else "须两条线同时命中",
            ),
            combine_reason=(
                "2026-09-21 用户裁定（GP8b 后）。面板线的前提是「面板能在本谱系打出可分辨的分」，"
                "而该前提是**为上皮冻结时成立的**。GP8b 实跑证明它在低表达谱系不成立：法则2 面板"
                "自身基因在 T/NK 里只有 7.76/1e4（上皮 54.4、B/浆 208.5、成纤维 26.4、髓系 24.9、"
                "内皮 13.7），T/NK 保留簇的面板自评中位仅 0.43，score_genes 扣掉对照集后 argmax "
                "在噪声里挑最大值 —— T/NK 被标的 11 簇有 10 簇是**单面板线**命中，而 CellTypist "
                "对同一批细胞报 98–99% T/NK。照旧口径会错删 32% 的 T 细胞。"
                "**指出规则的作用域 ≠ 事后调阈值**（法则 3.2）：上皮口径一字未改，非上皮改为两线"
                "同时命中方剔；CellTypist 线阈值 0.50、面板线 ctrl_size/seed 均一字未动，"
                "单线命中仍如实落表供人工判读。"
            ),
            frozen=True, frozen_note="规则在本脚本写定时即冻结；CellTypist 线结果先被看到，"
                                     "面板线完成前规则未改。合并口径的作用域修订为 2026-09-21 "
                                     "GP8b 之后、在**看到非上皮谱系结果之后**所做，已如实记入 "
                                     "combine_reason，不掩饰其时序。",
        ),
        score_genes=dict(ctrl_size=CTRL_SIZE, random_state=SCORE_SEED),
        lineage_genes_used=used, lineage_genes_missing=miss,
        verdict=dict(
            n_nonepi=int(len(bad)), cells_nonepi=int(bad["n_cells"].sum()),
            both_lines=list(both["cluster"]), one_line=list(one["cluster"]),
            clusters=df.set_index("cluster")["flag"].to_dict(),
            n_single_line_held=int(len(hold)),
            cells_single_line_held=int(hold["n_cells"].sum()),
            single_line_held=dict(zip(hold["cluster"].astype(str), hold["single_line_hold"])),
        ),
        inputs={os.path.relpath(p, ROOT): sha256(p) for p in
                [H5, bc_file, f"{TRAD}/{a.tag}/clusters.csv.gz", LABELS,
                 f"{HERE}/marker_panel.py", f"{HERE}/06_contamination_check.py"]},
        outputs={os.path.relpath(p, ROOT): sha256(p) for p in
                 [out_csv, f"{OUT}/{a.tag}_lineage_scores_by_cluster.csv"]},
        note="本脚本只诊断，不改数据；不决定剔除后是否重聚类（检查点问题）。"
             "下游 03_subcluster_annotation.py 只认 flag 列 —— single_line_hold 列**不会**导致剔除。",
    ), open(out_man, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log(f"写出 {os.path.relpath(out_man, ROOT)}")


if __name__ == "__main__":
    main()
