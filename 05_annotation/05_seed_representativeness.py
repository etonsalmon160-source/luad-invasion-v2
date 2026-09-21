#!/usr/bin/env python3
"""种子代表性复核 —— 补 GP5 报告登记的一个洞（2026-09-21）

背景（GP5_report.md §十二 登记的已知缺口）：
    「种子 0 最有代表性」这个结论**只在全量上算过**（种子 0 对其余 4 个种子的平均 ARI 最高），
    **六个谱系上并未复核**。而 GP8a 的注释恰恰是在**谱系子群**上做的
    ⇒ 若某个谱系里种子 0 恰好最差，该谱系的注释会偏，且看不出来。

本脚本只做**标签比对**，不重跑任何聚类 —— 五个种子的簇标签已随各 run 落盘
（`results/04_integration/seurat_trad/<tag>/clusters.csv.gz` 里的 `harmony_res*_seed*` 列）。

────────────────────────────────────────────────────────────────────────────
预注册判据（**跑之前写下，不得看结果再改**）
────────────────────────────────────────────────────────────────────────────
在每个对象**自己的 r\\* 处**（r\\* 取自已签字的 `<tag>_rstar.json`，不在此处重算）：
    比较「种子 0 对其余 4 个种子的平均 ARI」与「5 个种子里最好的那个的同一指标」
        gap = best - seed0
        gap ≤ TOL   ⇒ 视为**并列**，种子 0 可用（沿用本项目已有的 0.01 并列口径，
                      与 §M3-A.3 破平规则同一个数，**不新发明阈值**）
        gap >  TOL   ⇒ 🔴 标记该对象：注释**须换种子**，升级人工决定
四个分辨率**全部列出**（供查看），但**判据只在 r\\* 处** —— r\\* 之外的分辨率本来就不用于注释。

自带校验（对不上即停）：
    每个 (对象, 分辨率) 的 10 个两两 ARI 的**均值**，必须等于 R 侧
    `resolution_metrics.csv` 的 `ari_seed_mean`（容差 1e-4）。
    这同时验证了 sklearn 的 ARI 与 R 侧 mclust 的实现一致（PARAMETERS 已有先例记载）。
"""

import hashlib
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

ROOT = "/home/eto/luad_v2"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
ANN = f"{ROOT}/results/05_annotation"

RES_GRID = [0.5, 0.6, 0.7, 0.8]
SEEDS = [0, 1, 2, 3, 4]
TOL = 0.01          # 与 §M3-A.3 破平规则同一个数；此处复用，不新发明
ARI_TOL = 1e-4      # 与 R 侧 ari_seed_mean（4 位小数）比对的容差
REF_SEED = 0        # 本项目当前用于注释的种子

OBJECTS = ["full", "epiA", "tnkA", "myeloidA", "endoA", "bplasmaA", "fibroA"]
LIN_CN = {"full": "全量", "epiA": "上皮", "tnkA": "T/NK", "myeloidA": "髓系",
          "endoA": "内皮", "bplasmaA": "B/浆", "fibroA": "成纤维"}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(*a):
    print(*a, flush=True)


def pairwise_aris(cols):
    """10 个两两 ARI + 逐种子对其余 4 个的平均。"""
    codes = {s: pd.factorize(cols[s])[0] for s in SEEDS}
    pairs = {}
    per_seed = {s: [] for s in SEEDS}
    for i, a in enumerate(SEEDS):
        for b in SEEDS[i + 1:]:
            v = float(adjusted_rand_score(codes[a], codes[b]))
            pairs[(a, b)] = v
            per_seed[a].append(v)
            per_seed[b].append(v)
    mean_to_others = {s: float(np.mean(per_seed[s])) for s in SEEDS}
    return pairs, mean_to_others, codes


def matched_agreement(codes_a, codes_b):
    """两套簇标签在**最优配对**下的一致率（Hungarian 配列联表）。

    ARI 只说『像不像』；这个数说『换过去之后，到底多少细胞会改标签』——
    对『注释会不会变』是更直接的代理。
    返回 (一致细胞数, 总细胞数, 一致率)。
    """
    ct = pd.crosstab(pd.Series(codes_a, dtype="int64"), pd.Series(codes_b, dtype="int64")).to_numpy()
    ri, ci = linear_sum_assignment(-ct)
    same = int(ct[ri, ci].sum())
    return same, len(codes_a), same / len(codes_a)


def main():
    log("=" * 78)
    log("种子代表性复核（只读标签，不重跑聚类）")
    log("=" * 78 + "\n")

    rows = []
    inputs, verdicts = {}, {}

    for tag in OBJECTS:
        clu_p = f"{TRAD}/{tag}/clusters.csv.gz"
        met_p = f"{TRAD}/{tag}/resolution_metrics.csv"
        if not (os.path.exists(clu_p) and os.path.exists(met_p)):
            raise SystemExit(f"🔴 {tag}: 缺 {clu_p} 或 {met_p}")

        usecols = lambda c: c.startswith("harmony_res")
        clu = pd.read_csv(clu_p, usecols=usecols)
        met = pd.read_csv(met_p)
        # r* 来源分两处：六个谱系在已签字的 <tag>_rstar.json；
        # 全量没有签字 JSON，其 r* 记在自己 run_manifest.json 的 metrics.rstar_candidate。
        if tag == "full":
            manp = f"{TRAD}/full/run_manifest.json"
            r_star = float(json.load(open(manp, encoding="utf-8"))["metrics"]["rstar_candidate"])
            rstar_src = os.path.relpath(manp, ROOT)
        else:
            rstar_src = f"{ANN}/{tag}_rstar.json"
            r_star = float(json.load(open(rstar_src, encoding="utf-8"))["r_star"])
        log(f"=== {LIN_CN[tag]}（{tag}）r* = {r_star}  n_cells={len(clu):,} ===")

        for r in RES_GRID:
            cols = {s: clu[f"harmony_res{r}_seed{s}"].to_numpy() for s in SEEDS}
            pairs, m2o, codes = pairwise_aris(cols)
            n_clu = {s: int(pd.Series(cols[s]).nunique()) for s in SEEDS}

            # 只在 r* 处算「换种子会让多少细胞改标签」——对「注释会不会变」最直接的代理
            at_rstar = abs(r - r_star) < 1e-9
            best_seed_pre = max(SEEDS, key=lambda s: m2o[s])
            if at_rstar:
                same, tot, agree = matched_agreement(codes[REF_SEED], codes[best_seed_pre])
            else:
                same = tot = None
                agree = None

            # 自带校验：10 个两两 ARI 的均值 == R 侧 ari_seed_mean
            mine = float(np.mean(list(pairs.values())))
            row = met[met["resolution"] == r]
            if len(row) != 1:
                raise SystemExit(f"🔴 {tag} r={r}: resolution_metrics.csv 命中 {len(row)} 行")
            theirs = float(row.iloc[0]["ari_seed_mean"])
            if abs(mine - theirs) > ARI_TOL:
                raise SystemExit(f"🔴 {tag} r={r}: 自算两两均值 {mine:.6f} ≠ "
                                 f"R 侧 ari_seed_mean {theirs:.6f}（容差 {ARI_TOL}）—— 实现不一致，停")

            best_seed = max(SEEDS, key=lambda s: m2o[s])
            gap = m2o[best_seed] - m2o[REF_SEED]
            rank0 = 1 + sum(1 for s in SEEDS if m2o[s] > m2o[REF_SEED] + 1e-12)
            for s in SEEDS:
                rows.append(dict(object=tag, lineage=LIN_CN[tag], res=r, seed=s,
                                 n_cells=len(clu), n_clusters=n_clu[s],
                                 mean_ari_to_other_seeds=round(m2o[s], 6),
                                 is_rstar_res=bool(at_rstar),
                                 best_seed=best_seed, best_mean=round(m2o[best_seed], 6),
                                 ref_seed=REF_SEED, ref_mean=round(m2o[REF_SEED], 6),
                                 gap_best_minus_ref=round(gap, 6),
                                 ref_rank=rank0,
                                 tie_within_tol=bool(gap <= TOL),
                                 ref_vs_best_matched_cells=same,
                                 ref_vs_best_matched_frac=(None if agree is None
                                                           else round(agree, 6))))
            log(f"  r={r}  两两均值自算 {mine:.4f} == R 侧 {theirs:.4f} ✅   "
                f"最好种子 {best_seed}({m2o[best_seed]:.4f})  种子0({m2o[REF_SEED]:.4f}, "
                f"第 {rank0} 名)  gap={gap:+.4f}"
                + (f"   换种子后改标签 {tot - same:,}/{tot:,} 细胞"
                   f"（一致 {agree * 100:.2f}%）  ← r*" if at_rstar else ""))

            if at_rstar:
                verdicts[tag] = dict(lineage=LIN_CN[tag], r_star=r_star, n_cells=tot,
                                     best_seed=best_seed,
                                     best_mean=round(m2o[best_seed], 6),
                                     ref_mean=round(m2o[REF_SEED], 6),
                                     gap=round(gap, 6), ref_rank=rank0,
                                     ref_vs_best_matched_frac=round(agree, 6),
                                     cells_changing_label=int(tot - same),
                                     pass_=bool(gap <= TOL))

        inputs[f"{tag}/clusters.csv.gz"] = sha256(clu_p)
        inputs[f"{tag}/resolution_metrics.csv"] = sha256(met_p)
        inputs[os.path.relpath(rstar_src, ROOT)] = sha256(rstar_src)
        log("")

    df = pd.DataFrame(rows)
    at_rstar = df[df["is_rstar_res"]].drop_duplicates(subset=["object"])

    log("=" * 78)
    log(f"判据（预注册）：r* 处 gap = 最好种子 − 种子{REF_SEED}；gap ≤ {TOL} 视为并列 ⇒ 种子{REF_SEED} 可用")
    log("=" * 78)
    log(f"{'对象':<8}{'r*':>5}{'最好种子':>9}{'最好均值':>10}{'种子0均值':>11}"
        f"{'gap':>9}{'种子0排名':>9}{'换种子后改标签':>16}  判定")
    bad = []
    for _, r in at_rstar.iterrows():
        ok = r["gap_best_minus_ref"] <= TOL
        if not ok:
            bad.append(r["object"])
        frac = r["ref_vs_best_matched_frac"]
        chg = f"{frac * 100:.2f}% 一致" if frac is not None else "—"
        log(f"{r['lineage']:<8}{r['res']:>5}{int(r['best_seed']):>9}"
            f"{r['best_mean']:>10.4f}{r['ref_mean']:>11.4f}"
            f"{r['gap_best_minus_ref']:>+9.4f}{int(r['ref_rank']):>9}{chg:>16}  "
            + ("✅ 并列，种子 0 可用" if ok else "🔴 掉队，须换种子"))

    log("")
    if bad:
        log(f"🔴 结论：{len(bad)} 个对象在 r* 处，种子 {REF_SEED} 掉队超过 {TOL}："
            + "、".join(LIN_CN[t] for t in bad))
        log("   ⇒ 按预注册判据，这些对象**须换种子**。建议改用的种子：")
        for _, r in at_rstar.iterrows():
            if r["object"] in bad:
                log(f"      {r['lineage']:<8} r*={r['res']}  种子 {REF_SEED}({r['ref_mean']:.4f}) "
                    f"→ 种子 {int(r['best_seed'])}({r['best_mean']:.4f})   gap={r['gap_best_minus_ref']:+.4f}   "
                    f"两套划分最优配对下一致 {r['ref_vs_best_matched_frac'] * 100:.2f}%"
                    f"（即约 {int((1 - r['ref_vs_best_matched_frac']) * r['n_cells']):,} 个细胞会换簇）")
        log(f"   ⇒ ⚠️ 但**先别急着重跑**：GP8a 注释**尚未开始**，此刻换种子**不产生额外成本**；")
        log(f"      而 r* 是用**跨全部种子的平均** ARI 选的（与用哪个种子无关），故换种子**不动 r***。")
    else:
        log(f"✅ 结论：七个对象在各自 r* 处，种子 {REF_SEED} 与最好种子的差距都 ≤ {TOL}")
        log(f"   ⇒ 「种子 {REF_SEED} 有代表性」在**六个谱系上也成立**（此前只在全量上验过）。")
        log(f"   ⇒ GP5_report.md §十二 登记的那个洞**可以关掉**，注释继续用种子 {REF_SEED}。")

    out_csv = f"{ANN}/seed_representativeness.csv"
    df.to_csv(out_csv, index=False)

    man = dict(
        script=os.path.relpath(__file__, ROOT),
        generated="2026-09-21",
        purpose="补 GP5_report.md §十二 登记的洞：种子 0 的『代表性』原先只在全量上验过，六个谱系未验。",
        method="只读各 run 已落盘的簇标签（clusters.csv.gz 的 harmony_res*_seed* 列），"
               "算 10 个两两 ARI 与逐种子对其余 4 个的均值。**不重跑任何聚类。**",
        preregistered_rule=dict(
            where="每个对象自己的 r* 处（r* 取自已签字的 <tag>_rstar.json）",
            metric="gap = max_s(mean_ari_to_other_seeds) − mean_ari_to_other_seeds[seed 0]",
            tol=TOL,
            tol_provenance="复用 §M3-A.3 破平规则的 0.01，**不新发明阈值**",
            pass_="gap ≤ TOL 视为并列 ⇒ 种子 0 可用",
            fail="gap > TOL ⇒ 该对象须换种子并重跑该谱系注释",
        ),
        selfcheck=dict(claim="每个 (对象,分辨率) 的 10 个两两 ARI 均值 == R 侧 ari_seed_mean",
                       tol=ARI_TOL, status="PASS（否则脚本 SystemExit）"),
        verdict=dict(any_fail=bool(bad), failed_objects=bad,
                     note="种子 0 在七个对象各自的 r* 处均未掉队" if not bad else
                          "见 failed_objects，须换种子",
                     per_object=verdicts,
                     recommended=[dict(object=t, lineage=LIN_CN[t], r_star=verdicts[t]["r_star"],
                                       n_cells=verdicts[t]["n_cells"],
                                       from_seed=REF_SEED, to_seed=int(verdicts[t]["best_seed"]),
                                       from_mean=verdicts[t]["ref_mean"],
                                       to_mean=verdicts[t]["best_mean"],
                                       gap=verdicts[t]["gap"],
                                       ref_vs_best_matched_frac=verdicts[t]["ref_vs_best_matched_frac"],
                                       cells_changing_label=verdicts[t]["cells_changing_label"],
                                       ref_rank=verdicts[t]["ref_rank"]) for t in bad],
                     action_note="GP8a 注释尚未开始 ⇒ 此刻换种子不产生额外成本；"
                                 "r* 由跨全部种子的平均 ARI 选出，与种子无关 ⇒ 换种子不动 r*。"),
        scope_note="⚠️ 判据用的是**簇标签的接近程度**（ARI 与最优配对一致率），"
                   "不是**注释结果**是否相同。两套划分里若某个簇头基因没变，注释就一样；"
                   "反差出现在『换了标签但 marker 谱几乎没变』的情形。"
                   "要直接看注释差异，须用两个种子各跑一次注释（本脚本不做，成本 ~41 min/次）。",
        caveat_fibroA="成纤维的 r* 是**事后放宽**阈值后取的（relaxed=true），其种子间 ARI 本就 <0.90；"
                      "该谱系即使判为并列，产物仍只能作探索性结论。",
        inputs=inputs,
        outputs={os.path.relpath(out_csv, ROOT): sha256(out_csv)},
    )
    with open(f"{ANN}/seed_representativeness_manifest.json", "w", encoding="utf-8") as fh:
        json.dump(man, fh, ensure_ascii=False, indent=2)

    log(f"\n写出 {os.path.relpath(out_csv, ROOT)} / "
        f"{os.path.relpath(f'{ANN}/seed_representativeness_manifest.json', ROOT)}")


if __name__ == "__main__":
    main()
