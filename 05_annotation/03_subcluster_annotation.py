#!/usr/bin/env python3
"""GP8 · 逐谱系亚聚类 —— **簇级**经典面板注释 + CellTypist 对比。

口径（用户 2026-09-17 裁定）
---------------------------
「**先聚类，再用经典面板按簇注释**」——与 Travaglini 2020 Methods 的原procedure一致
（"Clusters were assigned a canonical identity based on enriched expression of these
marker genes."）。**不是**逐细胞 argmax 后投票。

两个 stage
----------
1. `--stage select`  ：读该谱系 run 的 `resolution_metrics.csv`，按 §M3-A.3 预注册规则
                     算 `r*`（含破平规则），打印表。**不产标签**，供签字。
2. `--stage annotate --res X`：对已签字的 r* 做簇级注释。产出：
   - `<tag>_cluster_scores.csv`   逐簇 × 逐亚型的**全部**平均得分（不只 argmax）
   - `<tag>_cluster_annotation.csv` 逐簇 判读 + CellTypist 众数 + 规模 + 跨患者分布
   - `<tag>_cluster_topgenes.csv`  逐簇 wilcoxon top 富集基因
   - `<tag>_subtype_disagreement.csv` 两口径不一致的簇（**绝不自动裁决**）
   - `<tag>_annotation_manifest.json`

诚实约束（写进代码，不靠自觉）
------------------------------
- **只报不判**：脚本给出证据（得分、top 基因、患者分布），`decision` 列留空待人工填。
- **稀有型**：胜出簇 < 20 细胞 ⇒ 标 `rare/likely-spurious`。
- **面板粗的亚型**（1–2 个基因槽）在报告里单列。
- **无法注释的亚型**（`UNANNOTATABLE`）不参与 argmax，单独在报告里说明。

用法:
    python3 05_annotation/03_subcluster_annotation.py --stage select --tag epiA --lineage 上皮
    python3 05_annotation/03_subcluster_annotation.py --stage annotate --tag epiA --lineage 上皮 --res 0.6
"""

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classic_panels as CP  # noqa: E402

ROOT = "/home/eto/luad_v2"
H5 = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
OUT = f"{ROOT}/results/05_annotation"
SUBSETS = f"{OUT}/lineage_subsets_manifest.json"

# 与 GP6 / §M3-A.3 逐字一致，不得在本脚本内另立口径
CTRL_SIZE = 50
SCORE_SEED = 0
SEEDS = [0, 1, 2, 3, 4]
ARI_SEED_MIN = 0.90
TIE_BREAK_EPS = 0.01      # 得分差 < 0.01 取**较低**分辨率（§M3-A.3 破平规则）
AAH_ABSORPTION_MAX = 0.50
RARE_MIN_CELLS = 20
DRIVER_TOPN = 15

T0 = time.time()


def log(m):
    print(f"[{time.time()-T0:8.1f}s] {m}", flush=True)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


# ================================================================ stage 1 · r* 选择
def stage_select(tag, relax_seed_min=None):
    """按 §M3-A.3 预注册规则算 r*。规则本身在 R 脚本里已算好 seed/xres 两列，

    `relax_seed_min` 是**事后放宽**开关，只允许用户明确裁定后使用（如 2026-09-17
    用户对成纤维的裁定）。一旦启用，产物 json 会带 `relaxed=true` 与放宽前后阈值，
    下游注释产物必须据此标注为"低稳定性/探索性"。**绝不默认启用。**
    """
    m = pd.read_csv(f"{TRAD}/{tag}/resolution_metrics.csv")
    log(f"{tag} 分辨率表：\n{m.to_string(index=False)}")

    thr = float(relax_seed_min) if relax_seed_min is not None else float(ARI_SEED_MIN)
    if relax_seed_min is not None:
        log(f"\n🔴 **事后放宽已启用**：跨种子稳定性阈值 {ARI_SEED_MIN} → {thr}。"
            f"此放宽须有用户裁定记录，且结果只能作探索性结论。")

    # 规则：score = 0.5*seed_stab + 0.5*xres；受两硬约束；并列(≤0.01)取较低分辨率
    #
    # ⚠️ 2026-09-17 修：原实现只用了「跨种子稳定性」一条硬约束，**丢掉了 AAH 吸收护栏**
    #   （pass_aah），而 R 侧的 `eligible` 列与 PARAMETERS §M3-A.3 登记的都是**两条**。
    #   GP8 门要求「子集内套用 GP5 判据」，故此处补齐。本批六谱系 pass_aah 全为 TRUE，
    #   补齐后 r* 不变 —— 属"规则实现与登记不符"，不是"结果错"。
    ok = m[m["ari_seed_mean"] >= thr].copy()
    if "aah_absorption_rate" in m.columns:
        # 用**原始指标列**判定，不复用 R 派生的 pass_aah —— 两侧各自从同一列重算，
        # 才不会出现「R 用未舍入值、Python 用舍入值」这类边界分叉。
        na_aah = ok["aah_absorption_rate"].isna()
        if na_aah.any():
            log(f"  ⚠️ 有 {int(na_aah.sum())} 个分辨率的 AAH 吸收率为 NA（子集内无 AAH 细胞）"
                f"⇒ 该条护栏**对该分辨率不适用**，予以保留并在此声明")
        keep_aah = ok["aah_absorption_rate"] <= AAH_ABSORPTION_MAX
        keep_aah |= na_aah
        dropped = ok.loc[~keep_aah, "resolution"].tolist()
        ok = ok[keep_aah]
        if dropped:
            log(f"  AAH 吸收护栏剔除分辨率 {dropped}（阈值 ≤ {AAH_ABSORPTION_MAX}，"
                f"PARAMETERS §M3-A.3 第二条硬约束）")
    else:
        raise SystemExit("🔴 resolution_metrics.csv 无 aah_absorption_rate 列，"
                         "无法施加已登记的第二条硬约束 —— 停在检查点，不静默跳过")
    if ok.empty:
        raise SystemExit(f"🔴 {tag}：没有任何分辨率同时满足跨种子稳定性 ≥ {thr} "
                         f"与 AAH 吸收 ≤ {AAH_ABSORPTION_MAX}，停在检查点")

    best = ok["score"].max()
    tied = ok[ok["score"] >= best - TIE_BREAK_EPS].sort_values("resolution")
    r_star = float(tied.iloc[0]["resolution"])   # 取较低者
    log(f"\n最高分 {best:.4f}；并列（差 ≤{TIE_BREAK_EPS}）分辨率 {list(tied['resolution'])} "
        f"⇒ **r* = {r_star}**（取较低，§M3-A.3 破平规则）")
    if len(tied) > 1:
        log(f"  并列说明：{dict(zip(tied['resolution'], tied['score'].round(4)))}")

    out = dict(tag=tag, r_star=r_star, best_score=float(best),
               eligible=[float(x) for x in ok["resolution"]],
               tie_break="差 ≤ 0.01 取较低分辨率（§M3-A.3）",
               rule="r* = argmax_r [0.5·ARI_seed + 0.5·ARI_xres]，受 ari_seed_mean ≥ 阈值 "
                    "与 aah_absorption_rate ≤ %g 两条约束；指标3 已补算（2026-09-21，口径 R_mean）"
                    "且对 r* 无区分力" % AAH_ABSORPTION_MAX,
               constraints_applied=[
                   "ari_seed_mean ≥ %g（= 指标1 跨种子稳定性；**承重**）" % thr,
                   "aah_absorption_rate ≤ %g（= 指标4 AAH 吸收护栏；**已施加但本批实测全分辨率为 0 "
                   "⇒ 空转**，按 PARAMETERS §M3-A.3 第 241/243 行的裁定**不得计入 r* 的通过理由**）"
                   % AAH_ABSORPTION_MAX,
                   "指标3（谱系覆盖 ≥ 0.90）**已补算（2026-09-21）** —— 口径「逐基因检出率均值」"
                   "（R_mean，见 PARAMETERS §M3-A.3 的 2026-09-21 裁定）；七个对象 × 4 分辨率 × "
                   "5 种子全过线（最低 0.9310）⇒ **未淘汰任何分辨率，对 r* 无区分力**。"
                   "故本 r* 仍**只由指标1 + 指标2 承担**，不得记入指标3 的功劳。"
                   "⚠️ 原措辞「待注释后补算」是错的，见 GP5_report.md §6.1 的更正。",
                   "（前两条均自 resolution_metrics.csv 的**原始指标列**重算；"
                   "R 侧派生的 pass_seed / pass_aah / eligible 一律不复用，避免边界分叉）"],
               rstar_authoritative_over="resolution_metrics.csv 的 is_rstar 列。"
                   "R 侧按**预注册**阈值 0.90 计算，无法反映事后放宽 ⇒ 本条 relaxed=true 时，"
                   "该列的 is_rstar 必为全 FALSE，属**预期**，不是矛盾。",
               seed_min_applied=thr,
               relaxed=relax_seed_min is not None,
               seed_min_preregistered=float(ARI_SEED_MIN) if relax_seed_min is not None else None,
               relaxed_note=("🔴 事后放宽预注册阈值，经用户 2026-09-17 裁定。"
                             "该谱系结果只能作**探索性**结论，不得进主结论。"
                             if relax_seed_min is not None else None),
               signed_by=None,
               signed_note="🔴 待人工签字：确认 r* 无误后，把 signed_by 改为署名"
                           "（如 \"<姓名> 2026-09-18\"），`--stage annotate` 才会放行（§M3-A.5）。")
    p = f"{OUT}/{tag}_rstar.json"
    with open(p, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)
    log(f"写出 {os.path.relpath(p, ROOT)}  ——  🔴 须人工签字后才可进入 annotate")


# ================================================================ stage 2 · 簇级注释
def check_rstar_signature(tag, res):
    """签字闸门：r* 未签字 ⇒ 不许注释。§M3-A.5「注释必须等 r* 签字之后才做」。

    ⚠️ 2026-09-17 修：原实现只在 `select` 的日志里**口头**说「须人工签字后才可进入
       annotate」，代码里没有任何检查 —— annotate 可以拿任意 --res 直接跑，包括
       select 从没选出的分辨率。现改成硬闸门：
         ① `<tag>_rstar.json` 必须存在（select 还没跑过 ⇒ 停）；
         ② 必须带 `signed_by` 字段（人工在 json 里签名后为真 ⇒ 未签 ⇒ 停）；
         ③ `--res` 必须**逐位等于**签字文件里的 r_star（防事后换分辨率）。
       三处任一不满足即 SystemExit，不做任何静默放行。
    """
    p = f"{OUT}/{tag}_rstar.json"
    if not os.path.exists(p):
        raise SystemExit(f"🔴 没有 {os.path.relpath(p, ROOT)} —— "
                         f"r* 尚未经 --stage select 产出，拒绝注释（§M3-A.5）")
    j = json.load(open(p))
    if not j.get("signed_by"):
        raise SystemExit(f"🔴 {os.path.relpath(p, ROOT)} 尚无 `signed_by` 签字字段 —— "
                         f"按 §M3-A.5，r* 须人工签字后才可注释。停止。")
    if float(res) != float(j["r_star"]):
        raise SystemExit(f"🔴 --res {res} ≠ 签字文件里的 r* = {j['r_star']} —— "
                         f"拒绝用未签字的分辨率注释。停止。")
    if j.get("relaxed"):
        log(f"🔴 本次 r* 为**事后放宽**口径（seed_min={j.get('seed_min_applied')}，"
            f"签字人 {j['signed_by']}）⇒ 产物只能作探索性结论，已在 manifest 登记")
    return j


def stage_annotate(tag, lineage, res, seed):
    rstar = check_rstar_signature(tag, res)
    man = json.load(open(SUBSETS))
    if lineage not in man["lineages"]:
        raise SystemExit(f"谱系未登记：{lineage}")
    bc_file = f"{OUT}/{man['lineages'][lineage]['file']}"
    bcs = [ln for ln in open(bc_file).read().split("\n") if ln]
    log(f"{lineage} 子集 {len(bcs):,} 核")

    clu = pd.read_csv(f"{TRAD}/{tag}/clusters.csv.gz")
    col = f"harmony_res{res}_seed{seed}"
    if col not in clu.columns:
        raise SystemExit(f"集群表里没有 {col}；有 {[c for c in clu.columns if 'harmony' in c]}")
    log(f"簇列 {col}：{clu[col].nunique()} 个簇（{clu[col].dtype}）")
    assert len(clu) == len(bcs), (len(clu), len(bcs))

    # 面板（本谱系）
    ad = sc.read_h5ad(H5)
    assert ad.shape == (413697, 18069), ad.shape
    keep = ad.obs_names.isin(set(bcs))
    if int(keep.sum()) != len(bcs):
        raise SystemExit(f"🔴 命中 {int(keep.sum())} ≠ 清单 {len(bcs)}")
    sub = ad[keep].copy()
    del ad

    # 顺序对齐：清单顺序 = h5ad 顺序 = 集群表顺序（子集 run 保留原顺序）
    if not (clu["cell_barcode"].to_numpy() == sub.obs_names.to_numpy()).all():
        log("⚠️ 集群表与子集顺序不一致，按 barcode 重排")
        clu = clu.set_index("cell_barcode").loc[sub.obs_names].reset_index()

    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    panel, missing = CP.build_panel(lineage, sub.var_names)
    n_slot = CP.N_SLOTS[lineage]
    thin = [k for k, v in panel.items() if len(v) <= 2]
    unscorable = [k for k, v in panel.items() if len(v) < 2]
    log(f"面板 {len(panel)} 型；缺失 {missing}")
    if thin:
        log(f"🔴 只有 1–2 个基因槽的亚型（容错为零）：{thin}")
    if unscorable:
        log(f"🔴 可用基因 <2，跳过打分：{unscorable}")

    for st, gl in panel.items():
        if len(gl) < 2:
            continue
        sc.tl.score_genes(sub, gl, ctrl_size=CTRL_SIZE, random_state=SCORE_SEED,
                          score_name=f"s_{st}")

    cols = [c for c in sub.obs.columns if c.startswith("s_")]
    types = [c[2:] for c in cols]
    sub.obs["clu"] = pd.Categorical(clu.set_index("cell_barcode")
                                    .loc[sub.obs_names, col].astype(str).to_numpy())
    cl = sub.obs["clu"].astype(str).to_numpy()
    ucl = sorted(set(cl), key=lambda x: (len(x), x))

    # 逐簇 × 逐型 平均得分（**全部**，不只 argmax）
    M = sub.obs[cols].to_numpy(dtype=np.float64)
    cm = pd.DataFrame(np.vstack([M[cl == c].mean(axis=0) for c in ucl]), columns=types)
    cm.insert(0, "cluster", ucl)
    cm.insert(1, "n_cells", [int((cl == c).sum()) for c in ucl])
    cm.to_csv(f"{OUT}/{tag}_cluster_scores.csv", index=False)
    log(f"写出 {tag}_cluster_scores.csv（逐簇 × 全部 {len(types)} 型得分）")

    # argmax 与证据
    win = np.array(types)[cm[types].to_numpy().argmax(axis=1)]
    order = np.sort(cm[types].to_numpy(), axis=1)
    margin = order[:, -1] - order[:, -2] if len(types) > 1 else np.zeros(len(ucl))

    # 跨患者分布（判「一个患者撑起一个簇」这类假象）
    # ⚠️ 2026-09-17 修：原实现在 obs 缺 patient_id 时**静默**改用 barcode 字符串切分推断患者
    #   —— 这正是 R5「patient_id 与 sample_id 分层」要挡的那类静默回退，且推断错也不会报错。
    #   现改为直接硬报错；同时删掉赋值后从未被使用的 `stg`（死代码）。
    if "patient_id" not in sub.obs:
        raise SystemExit("🔴 表达对象 obs 里没有 patient_id —— 拒绝由 barcode 推断患者（R5）")
    pat = sub.obs["patient_id"].astype(str).to_numpy()
    n_pat, top_pat_frac = [], []
    for c in ucl:
        m = cl == c
        vc = pd.Series(pat[m]).value_counts()
        n_pat.append(int(vc.size))
        top_pat_frac.append(float(vc.iloc[0] / m.sum()) if m.sum() else np.nan)

    # wilcoxon top 基因（每簇）
    top_rows = []
    for c in ucl:
        sub.obs["_mk"] = (cl == c)
        try:
            sc.tl.rank_genes_groups(sub, "_mk", groups=[True], method="wilcoxon",
                                    n_genes=DRIVER_TOPN)
            g = list(sub.uns["rank_genes_groups"]["names"][True])
        except Exception as e:  # 簇极小或退化
            g = [f"<失败: {type(e).__name__}>"]
        top_rows.append(dict(cluster=c, top_genes=";".join(map(str, g))))
    pd.DataFrame(top_rows).to_csv(f"{OUT}/{tag}_cluster_topgenes.csv", index=False)
    log(f"写出 {tag}_cluster_topgenes.csv")

    ann = pd.DataFrame(dict(
        cluster=ucl, n_cells=cm["n_cells"], argmax=win, win_margin=np.round(margin, 4),
        n_patients=n_pat, top_patient_frac=np.round(top_pat_frac, 3),
    ))
    ann["flag"] = ""
    small = ann["n_cells"] < RARE_MIN_CELLS
    ann.loc[small, "flag"] = "tiny_cluster"
    ann.loc[ann["argmax"].isin(CP.RARE) & small, "flag"] = "rare/likely-spurious"
    ann.loc[ann["argmax"].isin(thin), "flag"] = (
        ann.loc[ann["argmax"].isin(thin), "flag"].radd("thin_panel;").str.strip(";"))
    ann["decision"] = ""        # 🔴 留空待人工填，脚本不自动裁决
    ann["reason"] = ""
    ann.to_csv(f"{OUT}/{tag}_cluster_annotation.csv", index=False)
    log(f"写出 {tag}_cluster_annotation.csv（decision 列留空待人工填）")

    print(f"\n=== {lineage} 簇级 argmax（{col}）===")
    print(ann.to_string(index=False))
    vc = ann.groupby("argmax")["n_cells"].agg(["count", "sum"])
    print(f"\n胜出分布:\n{vc.to_string()}")
    never = [t for t in types if t not in set(win)]
    print(f"\n从未胜出的亚型: {never}")
    if unscorable:
        print(f"无法打分（可用基因<2）: {unscorable}")
    thin_src = [k for k, v in n_slot.items() if v <= 2]
    print(f"🔴 源表面板本来就只有 1–2 槽的亚型: {thin_src}")

    mf = dict(tag=tag, lineage=lineage, res=res, seed=seed, col=col,
              n_cells=int(len(bcs)), n_clusters=len(ucl),
              rstar_signed_by=rstar["signed_by"],
              rstar_relaxed=bool(rstar.get("relaxed")),
              rstar_json_sha256=sha256(f"{OUT}/{tag}_rstar.json"),
              panel_sha256=sha256(f"{os.path.dirname(os.path.abspath(__file__))}/classic_panels.py"),
              subsets_manifest=sha256(SUBSETS),
              score_genes=dict(ctrl_size=CTRL_SIZE, random_state=SCORE_SEED),
              missing_genes=missing, thin_panel_types=thin,
              unscorable_types=unscorable,
              unannotatable=list(CP.UNANNOTATABLE),
              source_table_defects=list(CP.SOURCE_TABLE_DEFECTS))
    with open(f"{OUT}/{tag}_annotation_manifest.json", "w") as fh:
        json.dump(mf, fh, ensure_ascii=False, indent=2)
    log(f"写出 {tag}_annotation_manifest.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["select", "annotate"], required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--lineage", default=None)
    ap.add_argument("--res", type=float, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--relax-seed-min", type=float, default=None,
                    help="🔴 事后放宽跨种子稳定性阈值。仅限用户明确裁定的谱系使用；"
                         "产物会登记 relaxed=true，结果只能作探索性结论。")
    a = ap.parse_args()

    if a.stage == "select":
        stage_select(a.tag, relax_seed_min=a.relax_seed_min)
    else:
        if a.lineage is None or a.res is None:
            raise SystemExit("annotate 须给 --lineage 与 --res")
        stage_annotate(a.tag, a.lineage, a.res, a.seed)


if __name__ == "__main__":
    main()
