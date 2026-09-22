#!/usr/bin/env python3
"""换种子改标签率 —— 重算 GP8c 报告 §八 里那个不可复现的数字（2026-09-22）

背景（GP8c_report.md §八「上皮旧版过关作废」）：
    该段写了一句「换种子后改标签的比例从 15.91%（旧）跳到 **51.36%**（新）——过半细胞换簇。」
    审计发现：① 仓库内**没有任何脚本**产出 15.91 / 51.36；② 全库（含 logs/）除那句正文外
    **不含** 15.91 这个值 ⇒ 不可复现。

本脚本做四件事
--------------
**（一）查清 51.36% 的真实身份。** 结论（可复现）：它是上皮**新版**运行在 res=0.5、
    种子 0 对种子 1 的**裸标签不一致率**（不做簇配对，直接比整数标签）＝ 51.3562%。
    Leiden 的簇编号是任意的（换种子就重排），裸比会把「同一个簇换了个编号」也算成「改标签」，
    所以**它不是改标签率**。同一对在**配对定义**下是 12.13%。

**（二）在**同一个**配对定义下重报该量。** 定义**不新发明**：直接引用
    `05_seed_representativeness.py` 的 `matched_agreement()`（匈牙利最优配对），
    并与其已登记的产物 `seed_representativeness.csv` 做**硬自检**（对不上即 SystemExit）。

**（三）补一版「旧运行 vs 新运行」的配对改标签率**（GP8c §十 待办 6 的后半句）。
    ⚠️ 输入在 `results/_superseded/`（**不入库**，.gitignore 明令不得当哈希锚点）
    ⇒ 属**一次性历史比对**，旧副本一删就不可再生。

**（四）同分辨率下的旧 vs 新对比，并复核 §八 的 0.9438。** §八 的关键句是
    「把旧标签限制在那 133,384 个细胞上重算，得 0.9438」。本脚本把它**实算一遍**，
    并把「旧 vs 新」放到**同一分辨率**下并排 —— 因为旧运行的 r*=0.5 而新运行的 r*=0.7，
    跨分辨率直接比会把「分辨率效应」误读成「细胞集效应」。

口径来源：`results/05_annotation/seed_representativeness.csv`
          （ref=种子0，best=该分辨率下对其余 4 个种子平均 ARI 最高的种子；只在各对象 r* 处落数）
本脚本**只读标签，不重跑任何聚类**。
"""

import hashlib
import importlib.util
import itertools
import json
import os
from datetime import datetime

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/eto/luad_v2"
TRAD = f"{ROOT}/results/04_integration/seurat_trad"
ANN = f"{ROOT}/results/05_annotation"
SUP = f"{ROOT}/results/_superseded/2026-09-22_cluster_adjudication/04_integration/seurat_trad"
FIGDIR = f"{ROOT}/figures"

TAG = "epiA"
SEEDS = [0, 1, 2, 3, 4]
REF_SEED = 0          # 与 05_seed_representativeness.py 同一个常量（本项目当前注释用种子）
RSTAR_JSON = f"{ANN}/{TAG}_rstar.json"           # 新版 r*（已签字）
REG_CSV = f"{ANN}/seed_representativeness.csv"   # 已登记的配对口径产物（硬自检靶）


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(*a):
    print(*a, flush=True)


def load_registered_definition():
    """把 05_seed_representativeness.py 的定义原样取来用。

    绝不在这里重写一遍 —— 重写就是换口径（本项目的诚实准则）。
    """
    p = f"{ROOT}/05_annotation/05_seed_representativeness.py"
    spec = importlib.util.spec_from_file_location("_sr", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, p


def res_grid_of(df):
    return sorted(float(x) for x in
                  df.columns.to_series().str.extract(r"harmony_res([\d.]+)_seed\d")[0]
                  .dropna().unique())


def main():
    log("=" * 78)
    log("换种子改标签率 —— 重算（只读标签，不重跑聚类）")
    log("=" * 78 + "\n")

    sr, def_src = load_registered_definition()
    log(f"配对定义取自 {os.path.relpath(def_src, ROOT)} 的 matched_agreement()（匈牙利最优配对）")
    log(f"  来源自检 sha256 = {sha256(def_src)}\n")

    def churn(a, b):
        """(配对改标签率, 裸不一致率)"""
        _, _, agree = sr.matched_agreement(pd.factorize(a)[0], pd.factorize(b)[0])
        return 1.0 - agree, float((np.asarray(a) != np.asarray(b)).mean())

    new_p = f"{TRAD}/{TAG}/clusters.csv.gz"
    old_p = f"{SUP}/{TAG}/clusters.csv.gz"
    new = pd.read_csv(new_p, compression="gzip")
    r_star = float(json.load(open(RSTAR_JSON, encoding="utf-8"))["r_star"])
    reg = pd.read_csv(REG_CSV)
    reg_epi = reg[(reg["object"] == TAG) & (reg["is_rstar_res"])].iloc[0]
    best_seed = int(reg_epi["best_seed"])
    # ⚠️ 登记列的语义是**一致率**（matched agreement），改标签率 = 1 − 该值。
    reg_agree = float(reg_epi["ref_vs_best_matched_frac"])
    log(f"新版 {TAG}: n={len(new):,}  r*={r_star}  已登记 种子{REF_SEED} vs best种子{best_seed}："
        f"一致率 {reg_agree * 100:.4f}% ⇒ 改标签率 {(1 - reg_agree) * 100:.4f}%\n")

    new_res = res_grid_of(new)
    rows = []

    # ── （一）现行运行内部：全部种子对 × 全部分辨率 ─────────────────────────
    log("（一）新版运行内 · 全部种子对 × 全部分辨率（配对定义 vs 裸定义）")
    for r in new_res:
        cols = {s: new[f"harmony_res{r}_seed{s}"].to_numpy() for s in SEEDS}
        for a, b in itertools.combinations(SEEDS, 2):
            mat, raw = churn(cols[a], cols[b])
            rows.append(dict(run="new", tag=TAG, resolution=r, seed_a=a, seed_b=b,
                             n_cells=len(new), n_clusters_a=int(pd.Series(cols[a]).nunique()),
                             n_clusters_b=int(pd.Series(cols[b]).nunique()),
                             matched_churn_frac=round(mat, 6), raw_mismatch_frac=round(raw, 6),
                             cells_changing_matched=int(round(mat * len(new))),
                             is_rstar_res=bool(abs(r - r_star) < 1e-9),
                             is_registered_pair=bool(abs(r - r_star) < 1e-9
                                                     and a == REF_SEED and b == best_seed)))
    df_new = pd.DataFrame(rows)
    hit = df_new[(df_new.is_rstar_res) & (df_new.seed_a == REF_SEED) & (df_new.seed_b == best_seed)]
    if len(hit) != 1:
        raise SystemExit(f"🔴 在 r*={r_star} 处找不到种子对 ({REF_SEED},{best_seed}) —— 自检无法进行，停")
    got = float(hit.iloc[0]["matched_churn_frac"])
    if abs((1 - got) - reg_agree) > 1e-6:
        raise SystemExit(f"🔴 自检失败：本脚本配对一致率 {1 - got:.6f} ≠ 已登记 {reg_agree:.6f}"
                         f"（ref种子{REF_SEED} vs best种子{best_seed} @ r*={r_star}）"
                         f" —— 定义不一致，停")
    log(f"  ✅ 自检通过：r*={r_star} 处 种子{REF_SEED} vs 种子{best_seed} 配对一致率 "
        f"{(1 - got) * 100:.4f}%（== 已登记 {reg_agree * 100:.4f}%）⇒ 改标签率 {got * 100:.4f}%\n")

    suspect = df_new[(df_new.resolution == 0.5) & (df_new.seed_a == 0) & (df_new.seed_b == 1)].iloc[0]
    log(f"  🔎 报告里那个 51.36% 的身份：run=new, res=0.5, 种子0 vs 种子1，**裸**不一致率 "
        f"= {suspect.raw_mismatch_frac * 100:.4f}%（四舍五入即 51.36%）")
    log(f"     同一对在**配对**定义下 = {suspect.matched_churn_frac * 100:.4f}%"
        f"（≈ {suspect.cells_changing_matched:,} 个细胞改标签）")
    log(f"     ⇒ 51.36% 不是改标签率，是 Leiden 簇编号重排造成的假象。\n")

    # ── （二）（三）（四）需要旧运行 ────────────────────────────────────────
    part_b = dict(status="skipped", reason="旧副本不存在（results/_superseded/ 未入库，可能已删）")
    same_res = None
    restricted = None
    if os.path.exists(old_p):
        old = pd.read_csv(old_p, compression="gzip")
        shared = old.merge(new, on="cell_barcode", suffixes=("_o", "_n"))
        old_res = res_grid_of(old)
        log("（二）旧运行 vs 新运行（⚠️ 输入在 _superseded、**不入库** ⇒ 一次性历史比对）")
        log(f"  旧运行 {len(old):,} 核 / 新运行 {len(new):,} 核 / 共有 {len(shared):,} 核"
            f"（旧独有 {len(old) - len(shared):,}，新独有 {len(new) - len(shared):,}）")
        for r in old_res:
            for s in SEEDS:
                co, cn = f"harmony_res{r}_seed{s}_o", f"harmony_res{r}_seed{s}_n"
                if co not in shared.columns or cn not in shared.columns:
                    continue
                mat, raw = churn(shared[co].to_numpy(), shared[cn].to_numpy())
                rows.append(dict(run="old_vs_new_same_cells", tag=TAG, resolution=r, seed_a=s, seed_b=s,
                                 n_cells=len(shared), n_clusters_a=int(shared[co].nunique()),
                                 n_clusters_b=int(shared[cn].nunique()),
                                 matched_churn_frac=round(mat, 6), raw_mismatch_frac=round(raw, 6),
                                 cells_changing_matched=int(round(mat * len(shared))),
                                 is_rstar_res=False, is_registered_pair=False))

        # （三）旧运行自身（同定义）
        log("\n（三）旧运行自身 · 种子对（配对定义，供同分辨率对比）")
        old_means = {}
        for r in old_res:
            vals = []
            for a, b in itertools.combinations(SEEDS, 2):
                col_a, col_b = old[f"harmony_res{r}_seed{a}"], old[f"harmony_res{r}_seed{b}"]
                mat, raw = churn(col_a.to_numpy(), col_b.to_numpy())
                rows.append(dict(run="old", tag=TAG, resolution=r, seed_a=a, seed_b=b,
                                 n_cells=len(old),
                                 n_clusters_a=int(col_a.nunique()), n_clusters_b=int(col_b.nunique()),
                                 matched_churn_frac=round(mat, 6), raw_mismatch_frac=round(raw, 6),
                                 cells_changing_matched=int(round(mat * len(old))),
                                 is_rstar_res=bool(abs(r - 0.5) < 1e-9), is_registered_pair=False))
                vals.append(mat)
            old_means[r] = float(np.mean(vals))
            log(f"  res={r}: 配对改标签率（10 对）均值 {np.mean(vals) * 100:.2f}%  "
                f"范围 {min(vals) * 100:.2f}–{max(vals) * 100:.2f}%")

        # （四）同分辨率并列 + 复核 §八 的 0.9438
        log("\n（四）同分辨率下 旧 vs 新（关键：两个 run 的 r* 不同，跨分辨率比会把"
            "分辨率效应误读成细胞集效应）")
        new_means = {}
        ari_check = {}
        for r in sorted(set(old_res) & set(new_res)):
            ocols = {s: old[f"harmony_res{r}_seed{s}"].to_numpy() for s in SEEDS}
            ncols = {s: new[f"harmony_res{r}_seed{s}"].to_numpy() for s in SEEDS}
            o_ari = float(np.mean(list(sr.pairwise_aris({s: pd.Series(ocols[s]) for s in SEEDS})[0].values())))
            n_ari = float(np.mean(list(sr.pairwise_aris({s: pd.Series(ncols[s]) for s in SEEDS})[0].values())))
            ov = [churn(ocols[a], ocols[b])[0] for a, b in itertools.combinations(SEEDS, 2)]
            nv = [churn(ncols[a], ncols[b])[0] for a, b in itertools.combinations(SEEDS, 2)]
            new_means[r] = float(np.mean(nv))
            log(f"  res={r}:  旧 run  跨种子 ARI {o_ari:.4f} / 配对改标签率 {np.mean(ov) * 100:.2f}%"
                f"     新 run  ARI {n_ari:.4f} / 配对改标签率 {np.mean(nv) * 100:.2f}%")
            ari_check[r] = dict(old_ari=round(o_ari, 6), new_ari=round(n_ari, 6),
                                old_churn=round(float(np.mean(ov)), 6),
                                new_churn=round(float(np.mean(nv)), 6))
        same_res = ari_check

        # 复核 §八 的「把旧标签限制在那 133,384 个细胞上重算，得 0.9438」
        log("\n  复核 GP8c §八 的 0.9438（旧标签限制到共有的 133,384 核后重算跨种子 ARI）")
        restricted = {}
        for r in sorted(set(old_res) & set(new_res)):
            oc = {s: shared[f"harmony_res{r}_seed{s}_o"] for s in SEEDS}
            o_ari = float(np.mean(list(sr.pairwise_aris({s: pd.Series(oc[s]) for s in SEEDS})[0].values())))
            restricted[r] = round(o_ari, 6)
            flag = ""
            if abs(r - 0.5) < 1e-9:
                flag = ("  ← §八 写的 0.9438 "
                        + ("✅ 复现" if abs(o_ari - 0.9438) < 1e-3 else f"🔴 对不上（差 {o_ari - 0.9438:+.4f}）"))
            log(f"    res={r}: 旧标签@共有细胞 跨种子 ARI = {o_ari:.4f}{flag}")
        part_b = dict(status="computed", shared_cells=len(shared),
                      old_only=len(old) - len(shared), new_only=len(new) - len(shared),
                      old_run=dict(path=os.path.relpath(old_p, ROOT), sha256=sha256(old_p),
                                   caliber="A_frozen（旧口径）", r_star=0.5, archived=False,
                                   note="results/_superseded/ 不入库，.gitignore 明令不得当哈希锚点"
                                        "⇒ 本段数字属一次性历史比对，旧副本删除后不可再生"))
        log("")

    df = pd.DataFrame(rows)
    out_csv = f"{ANN}/seed_churn_audit.csv"
    df.to_csv(out_csv, index=False)

    # ── 图（供肉眼判读；放项目目录 figures/，不放 /tmp）──────────────────────
    os.makedirs(FIGDIR, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.5))

    ax = axes[0]
    dn = df_new[(df_new.seed_a == 0) & (df_new.seed_b == 1)].sort_values("resolution")
    x = np.arange(len(dn))
    ax.bar(x - 0.2, dn.raw_mismatch_frac * 100, 0.4, color="#c44", edgecolor="k", linewidth=0.4,
           label="raw label mismatch (no matching)")
    ax.bar(x + 0.2, dn.matched_churn_frac * 100, 0.4, color="#48a", edgecolor="k", linewidth=0.4,
           label="matched churn (registered definition)")
    ax.annotate("51.36% = raw mismatch\n(NOT a churn rate)",
                xy=(0 - 0.2, dn.raw_mismatch_frac.iloc[0] * 100), xytext=(0.05, 60),
                arrowprops=dict(arrowstyle="->", color="#900", lw=1.1), fontsize=8.5, color="#900")
    ax.set_xticks(x); ax.set_xticklabels([f"{r}" for r in dn.resolution])
    ax.set_xlabel("Leiden resolution"); ax.set_ylabel("% of cells")
    ax.set_title(f"{TAG} new run: seed0 vs seed1", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left"); ax.grid(axis="y", alpha=0.25)

    if same_res:
        for ax, r in zip(axes[1:], sorted(same_res)[:2]):
            ob = df[(df.run == "old") & (df.resolution == r) & (df.seed_a == 0)].sort_values("seed_b")
            nb = df[(df.run == "new") & (df.resolution == r) & (df.seed_a == 0)].sort_values("seed_b")
            if ob.empty or nb.empty or list(ob.seed_b) != list(nb.seed_b):
                continue
            x = np.arange(len(ob))
            ax.bar(x - 0.2, ob.matched_churn_frac * 100, 0.4, color="#999", edgecolor="k",
                   linewidth=0.4, label="OLD run (A_frozen)")
            ax.bar(x + 0.2, nb.matched_churn_frac * 100, 0.4, color="#48a", edgecolor="k",
                   linewidth=0.4, label="NEW run (A_adjudicated)")
            ax.set_xticks(x); ax.set_xticklabels([f"{s}" for s in ob.seed_b])
            ax.set_xlabel("other seed compared to seed0")
            ax.set_ylabel("% of cells changing label")
            star = "  (OLD r*)" if abs(r - 0.5) < 1e-9 else ("  (NEW r*)" if abs(r - r_star) < 1e-9 else "")
            ax.set_title(f"matched churn at res={r}{star}\nold ARI {same_res[r]['old_ari']:.4f}"
                         f" vs new ARI {same_res[r]['new_ari']:.4f}", fontsize=10)
            ax.legend(fontsize=7.5); ax.grid(axis="y", alpha=0.25)

    fig.suptitle("epiA seed churn: the raw-vs-matched artifact, and old vs new run "
                 "at the SAME resolution", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out_fig = f"{FIGDIR}/seed_churn_audit.png"
    fig.savefig(out_fig, dpi=160)
    plt.close(fig)

    man = dict(
        script=os.path.relpath(__file__, ROOT),
        generated=datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        purpose="重算 GP8c_report.md §八 那个不可复现的『换种子改标签率』（15.91% / 51.36%），"
                "并补一版『旧运行 vs 新运行』的配对改标签率。",
        question_from="用户 2026-09-22 裁定后要求处理 GP8c §十 待办 6。",
        definition=dict(
            source=os.path.relpath(def_src, ROOT), source_sha256=sha256(def_src),
            fn="matched_agreement()（匈牙利最优配对：列联表 → linear_sum_assignment）",
            registered_caliber="results/05_annotation/seed_representativeness.csv 的 "
                               "ref_vs_best_matched_frac —— ⚠️ 该列语义是**一致率**（非改标签率），"
                               "本脚本据此硬自检：1 − matched_churn_frac == 该列。",
            note="定义**原样引用**上游脚本，未在此重写；重写即换口径。",
            raw_counterpart="raw_mismatch_frac = 直接比整数簇标签，**无配对**。"
                            "Leiden 簇编号任意 ⇒ 把『同一簇换编号』也计为改标签，**不是改标签率**。"),
        findings=dict(
            identity_of_51_36=dict(value=float(suspect.raw_mismatch_frac),
                                   is_pair=dict(run="new", resolution=0.5, seed_a=0, seed_b=1,
                                                metric="raw_mismatch_frac"),
                                   same_pair_matched=float(suspect.matched_churn_frac),
                                   statement="51.36% 是**裸**标签不一致率，非改标签率；同一对在配对"
                                             f"定义下为 {suspect.matched_churn_frac * 100:.2f}%。"),
            identity_of_15_91=dict(value=None, status="UNREPRODUCIBLE",
                                   evidence="全库（含 logs/）grep '15.91' 仅命中该句正文本身，"
                                            "无任何产物/日志/脚本含该值 ⇒ 无法复现，标记作废。"),
            registered_value_at_rstar=dict(r_star=r_star, ref_seed=REF_SEED, best_seed=best_seed,
                                           matched_churn_frac=round(got, 6))),
        selfcheck=dict(claim=f"新版 r*={r_star} 处 种子{REF_SEED} vs 种子{best_seed} 的配对**一致率** "
                             f"== seed_representativeness.csv 已登记值",
                       tol=1e-6, status="PASS（否则 SystemExit）",
                       agreement_computed=round(1 - got, 6), agreement_registered=round(reg_agree, 6),
                       churn_frac=round(got, 6)),
        same_resolution_old_vs_new=same_res,
        old_labels_restricted_to_shared_cells_gp8c_section8=restricted,
        part_b_old_vs_new=part_b,
        inputs={os.path.relpath(p, ROOT): sha256(p) for p in
                [new_p, RSTAR_JSON, REG_CSV, def_src] + ([old_p] if os.path.exists(old_p) else [])},
        outputs={os.path.relpath(out_csv, ROOT): sha256(out_csv),
                 os.path.relpath(out_fig, ROOT): sha256(out_fig)},
        caveat="⚠️ 【二】old_vs_new_same_cells 的输入在 results/_superseded/（不入库、不得作哈希锚点）"
               "⇒ 该段属一次性历史比对，不可再生；【一】与旧运行自身那两段来自同一份旧副本，同样"
               "不可再生；新版部分的输入均在库内，可复现。图只示 seed0 对比，表内含全部种子对。",
    )
    with open(f"{ANN}/seed_churn_audit_manifest.json", "w", encoding="utf-8") as fh:
        json.dump(man, fh, ensure_ascii=False, indent=2)

    log(f"\n写出 {os.path.relpath(out_csv, ROOT)}\n"
        f"     {os.path.relpath(f'{ANN}/seed_churn_audit_manifest.json', ROOT)}\n"
        f"     {os.path.relpath(out_fig, ROOT)}")


if __name__ == "__main__":
    main()
