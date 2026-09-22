"""GP8c 收尾：被剔除的细胞的**谱系接手**（逐细胞复核，三态分流）。

要解决什么问题
--------------
`06_contamination_check.py` 判"某个亚簇是不是混进了别的谱系"时，用的是**簇级**证据：
① 六谱系面板打分取**簇内均值**再 argmax；② CellTypist 的 `B_lineage` 簇内**非本谱系占比 > 50%`。
两线命中即整簇剔除（上皮只需一条命中即剔）。

簇级判断为了压噪声，代价是**一刀切**：一个 30 细胞的簇被判"混了 T 细胞"，簇里那 20 个
真上皮也跟着被剔了。这一步就是把被剔的细胞**逐个**拿回来看一遍 —— 用**同一套判据**，
只是把聚合层级从"簇"降回"细胞"。

判据（**预先写明，不因结果调**；法则 3.2）
------------------------------------------
对每一个被剔细胞，拿两条**互相独立**的线（与 06 完全同源）：

  ① **面板线（逐细胞）**：`gp6_scores.csv.gz`里六个谱系的**逐细胞**模块分 argmax。
     这份打分就是 GP6 标准 A 的原始输入（`ctrl_size=50, random_state=0`），
     **一个参数都没新加**，只是不再做簇内平均。
  ② **CellTypist 线（逐细胞）**：`B_lineage`（独立图谱，与面板不同源）。

三态（状态枚举，不是可调阈值）：

  - `restored_to_<L>`  —— 两线都说是**原来那个谱系 L** ⇒ 这条细胞是**簇级误剔**，还回去。
  - `readopted_to_<X>` —— 两线都说是 **X（≠ L）** ⇒ 跨谱系接手，改判到 X。
  - `unadopted`         —— 两线**不一致**，或 CellTypist 判 `未归属` ⇒ **不自动裁**，
                          列进 `gp8c_readopt_unadopted.csv` 交人工判（用户 2026-09-22 决定）。

⚠️ 诚实披露（写进 manifest，报告里也要提）
------------------------------------------
本步的两条线**与剔除判据同源**（同一份面板打分、同一份 CellTypist 输出），是"同一证据换个
聚合层级再看一遍"，**不是独立验证**。因此：
  - `readopted_to_X` **不能**当作"独立第二证据"去支撑 X 的结论，只能说"降一级看还成立"；
  - `unadopted` 的**条数本身就是这次裁决不确定性的度量** —— 数越小，簇级判据越干净。

`unadopted` 的临时归属
---------------------
`unadopted` 的细胞**暂不移动**（留在 `A_adjudicated` 给它的那个谱系里）。理由：剔除是主动
动作、接手也是主动动作，人工没判之前"什么都不做"最保守。故本层是**增量层**，
`A_adjudicated` 与已签字的六份子集清单**一个字都不改**（同 09 的做法）。

用法:
    python3 05_annotation/10_readopt_excluded_cells.py            # 全部六个谱系
    python3 05_annotation/10_readopt_excluded_cells.py --tags epiA tnkA
"""

import argparse
import datetime
import hashlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = "/home/eto/luad_v2"
OUT = f"{ROOT}/results/05_annotation"

# ---- 输入钉哈希（与 00_build_lineage_subsets.py 的 SOURCES 一致）------------------
SOURCES = {
    # GP6 标准 A 的**逐细胞**六谱系模块分。参数 ctrl_size=50 / seed=0，见 gp6_manifest.json。
    "gp6_scores.csv.gz":
        "c657f0e4907c625082c043653746d88911d4929efcf5ff93090a6941fdaf9a9f",
    # 裁决层产物：提供 A_adjudicated（细胞属于哪个谱系）与 B_lineage（CellTypist 线）。
    "gp6_cell_labels_adjudicated.csv.gz":
        "dfaecc05c529bbe7cfb06990fd5cd13c99de03adf906fd1e534ded7e6c743093",
    # 裁决层的签字文件 —— 与 00_build_lineage_subsets.py 同一道闸门。
    "gp6_adjudication_manifest.json": None,   # 内容会变（签字段落），只查签字、不钉哈希
}

# CellTypist 的"没把握"不算"判成别的谱系"（与 06_contamination_check.py 的 CT_NA 同一口径）
CT_NA = {"未归属", None, "", "nan"}

# 谱系 → 被剔细胞清单文件名
TAGS = {
    "epiA": "上皮",
    "tnkA": "T/NK",
    "bplasmaA": "B/浆",
    "myeloidA": "髓系",
    "fibroA": "成纤维",
    "endoA": "内皮",
}

PANEL_COLS = ["上皮", "T/NK", "B/浆", "髓系", "成纤维", "内皮"]

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_registered(basename):
    """输入文件的 sha256 必须与登记值一致，否则拒绝继续（R5：换了输入要看得见）。"""
    path = f"{OUT}/{basename}"
    if not os.path.exists(path):
        raise SystemExit(f"🔴 找不到 {basename} —— 须先跑上游")
    want = SOURCES[basename]
    got = sha256(path)
    if want is None:
        return got
    if got != want:
        raise SystemExit(f"🔴 {basename} 的 sha256 与登记值不符：\n"
                         f"    实际 {got}\n    登记 {want}")
    log(f"{basename}  sha256 ✅ {got[:16]}…")
    return got


def load_reference():
    """标签表 + 逐细胞面板分，合成一张逐细胞参照表。"""
    adj = check_registered("gp6_cell_labels_adjudicated.csv.gz")
    scr = check_registered("gp6_scores.csv.gz")

    # 签字闸门：裁决层没签字就不许用它做接手（同 00_build_lineage_subsets.py）
    sp = f"{OUT}/gp6_adjudication_manifest.json"
    if not os.path.exists(sp):
        raise SystemExit("🔴 缺 gp6_adjudication_manifest.json —— 裁决层须先签字")
    sj = json.load(open(sp, encoding="utf-8"))
    if not sj.get("signed_by"):
        raise SystemExit("🔴 裁决层尚无 `signed_by` —— 拒绝用它做谱系接手")
    log(f"签字闸门 ✅ gp6_adjudication_manifest.json：signed_by = {sj['signed_by']}")
    man_sha = sha256(sp)

    lab = pd.read_csv(f"{OUT}/gp6_cell_labels_adjudicated.csv.gz",
                      usecols=["cell_barcode", "sample_id", "patient_id", "stage",
                               "A_frozen", "A_adjudicated", "B_lineage", "B_confidence"])
    if lab["cell_barcode"].duplicated().any():
        raise SystemExit("🔴 标签表里 barcode 有重复")

    s = pd.read_csv(f"{OUT}/gp6_scores.csv.gz", usecols=PANEL_COLS + ["cell_barcode"])
    if s["cell_barcode"].duplicated().any():
        raise SystemExit("🔴 面板分表里 barcode 有重复")

    ref = lab.merge(s, on="cell_barcode", how="left", validate="one_to_one")
    if ref[PANEL_COLS].isna().all(axis=1).any():
        n = int(ref[PANEL_COLS].isna().all(axis=1).sum())
        raise SystemExit(f"🔴 有 {n} 个细胞在面板分表里查不到（两表细胞集不一致）")

    # 逐细胞面板线：六谱系模块分 argmax（这就是 GP6 标准 A 的原始口径，不做簇内平均）
    pv = ref[PANEL_COLS].to_numpy(dtype=np.float64)
    order = np.argsort(-pv, axis=1, kind="stable")
    ref["panel_top"] = [PANEL_COLS[i] for i in order[:, 0]]
    ref["panel_score_top"] = pv[np.arange(len(pv)), order[:, 0]]
    ref["panel_score_second"] = pv[np.arange(len(pv)), order[:, 1]]
    ref["panel_margin"] = ref["panel_score_top"] - ref["panel_score_second"]
    # 六个分数全是 NaN ⇒ argmax 会静默返回第 0 列（"上皮"），必须显式改成弃权，
    # 否则"没打分"会被当成"判成上皮"。
    nan_all = np.isnan(pv).all(axis=1)
    if nan_all.any():
        log(f"⚠️ 有 {int(nan_all.sum())} 个细胞六个谱系分数全为 NaN ⇒ 面板线弃权")
        ref.loc[nan_all, "panel_top"] = ""
        ref.loc[nan_all, ["panel_score_top", "panel_score_second", "panel_margin"]] = np.nan

    ref["ct_top"] = ref["B_lineage"].where(~ref["B_lineage"].isin(CT_NA), other=None)
    log(f"逐细胞参照表 {len(ref):,} 行：面板线 argmax 分布 "
        f"{ref['panel_top'].value_counts().to_dict()}")
    return ref.set_index("cell_barcode", drop=False), dict(adj=adj, scr=scr, man=man_sha)


def classify(ref_row, lineage):
    """三态分流。返回 (state, target) —— target 为 None 表示 unadopted。"""
    p, c = ref_row["panel_top"], ref_row["ct_top"]
    if p is None or c is None or pd.isna(c) or p == "":
        return "unadopted", None          # 任一线弃权 ⇒ 不自动裁
    if p != c:
        return "unadopted", None          # 两线打架 ⇒ 不自动裁
    if p == lineage:
        return "restored", lineage        # 两线都说是原谱系 ⇒ 簇级误剔
    return "readopted", p                 # 两线都说是别的谱系 ⇒ 跨谱系接手


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", nargs="*", default=list(TAGS),
                    help="只处理这些 tag（默认全部六个谱系）")
    a = ap.parse_args()
    bad = set(a.tags) - set(TAGS)
    if bad:
        raise SystemExit(f"🔴 未登记的 tag：{bad}")

    log("载入逐细胞参照表")
    ref, in_sha = load_reference()

    rows, per_tag, ex_sha = [], {}, {}
    for tag in a.tags:
        lineage = TAGS[tag]
        p = f"{OUT}/{tag}_excluded_cells.csv.gz"
        ann_mf = f"{OUT}/{tag}_annotation_manifest.json"
        if not os.path.exists(p):
            # ⚠️ 2026-09-22 修：原实现一见"清单不存在"就当成"该谱系没有污染簇"。
            # 这是把**上游没跑**误当成**上游跑了且没发现污染** —— 2026-09-22 夜跑正是这样：
            # S2 卡在检查点 → S5/S6 全跳过 → 清单根本没生成，而本脚本却报出
            # "六条谱系都没有被剔细胞"这种听上去像好消息的话。**没证据 ≠ 无罪**，
            # 必须靠上游注释 manifest 来区分。
            if not os.path.exists(ann_mf):
                raise SystemExit(
                    f"🔴 {tag}：既无剔除清单、也无注释 manifest —— 无法区分"
                    f"「上游没跑」与「上游跑了且未发现污染」。\n"
                    f"    拒绝把'没跑'报成'无污染'。\n"
                    f"    缺 {os.path.relpath(p, ROOT)}\n"
                    f"    缺 {os.path.relpath(ann_mf, ROOT)}")
            ex_sha[f"{tag}_annotation_manifest.json"] = sha256(ann_mf)
            aj = json.load(open(ann_mf, encoding="utf-8"))
            c = aj.get("contamination")
            if c and int(c.get("n_cells_excluded", 0)) > 0:
                raise SystemExit(f"🔴 {tag}：注释 manifest 记了 {c['n_cells_excluded']:,} 核被剔、"
                                 f"却找不到剔除清单 —— 上游产物自相矛盾，停")
            log(f"{tag}: 上游确认无污染簇（{tag}_annotation_manifest.json："
                f"contamination.n_cells_excluded=0）")
            per_tag[tag] = dict(lineage=lineage, n_excluded=0, states={},
                                excluded_clusters=[], note="上游确认无污染簇")
            continue
        ex_sha[f"{tag}_excluded_cells.csv.gz"] = sha256(p)
        ex = pd.read_csv(p)
        if ex.empty:
            log(f"{tag}: 剔除清单为空（0 核）")
            per_tag[tag] = dict(lineage=lineage, n_excluded=0, states={},
                                excluded_clusters=[], note="清单为空")
            continue

        miss = set(ex["cell_barcode"]) - set(ref.index)
        if miss:
            raise SystemExit(f"🔴 {tag} 有 {len(miss)} 个被剔细胞在标签表里查不到")

        sub = ref.loc[ex["cell_barcode"].to_numpy()].reset_index(drop=True)
        # 一致性硬断言：被剔细胞必然属于本谱系。不成立说明上游口径串了，拒绝继续。
        wrong = sub["A_adjudicated"] != lineage
        if wrong.any():
            raise SystemExit(f"🔴 {tag} 有 {int(wrong.sum())} 个被剔细胞的 A_adjudicated "
                             f"不是 {lineage}：{sub.loc[wrong, 'A_adjudicated'].unique()[:5]}")

        ev = ex.set_index("cell_barcode")["evidence"]
        for i in range(len(sub)):
            r = sub.iloc[i]
            state, target = classify(r, lineage)
            rows.append(dict(
                cell_barcode=r["cell_barcode"], excluded_from=lineage, tag=tag,
                cluster=ex.iloc[i]["cluster"],
                patient_id=r["patient_id"], sample_id=r["sample_id"], stage=r["stage"],
                exclusion_evidence=ev.get(r["cell_barcode"], ""),
                A_frozen=r["A_frozen"], A_adjudicated=r["A_adjudicated"],
                panel_top=r["panel_top"],
                panel_score_top=round(float(r["panel_score_top"]), 4),
                panel_margin=round(float(r["panel_margin"]), 4),
                panel_self_score=round(float(r[lineage]), 4),
                ct_top=(None if pd.isna(r["ct_top"]) else r["ct_top"]),
                ct_top_confidence=(None if pd.isna(r["B_confidence"]) else round(float(r["B_confidence"]), 4)),
                state=state, target_lineage=target,
                # 接手后的临时谱系：unadopted 不动（留在 A_adjudicated），其余按 target
                A_readopted=(lineage if target is None else target),
            ))
        n = len(sub)
        st = pd.Series([x["state"] for x in rows[-n:]])
        per_tag[tag] = dict(
            lineage=lineage, n_excluded=int(n),
            states=st.value_counts().to_dict(),
            excluded_clusters=sorted(ex["cluster"].astype(str).unique().tolist()),
            note="")
        log(f"{tag}（{lineage}）被剔 {n:,} 核 → " +
            "，".join(f"{k} {v}" for k, v in st.value_counts().items()))

    if not rows:
        log("六条谱系都没有被剔细胞 —— 本层为空（这本身是个好结果，报告里要写）")
        df = pd.DataFrame(columns=["cell_barcode", "excluded_from", "state"])
    else:
        df = pd.DataFrame(rows)
        df.sort_values(["excluded_from", "state", "cluster", "cell_barcode"], inplace=True)

    # ---- 逐谱系互斥性：同一细胞不可能被两个谱系剔（子集本身就是互斥的）----
    if len(df):
        dup = df["cell_barcode"].duplicated()
        if dup.any():
            raise SystemExit(f"🔴 有 {int(dup.sum())} 个细胞出现在多个谱系的剔除清单里"
                             f" —— 上游子集不互斥，停")

    out_detail = f"{OUT}/gp8c_readopt_excluded_cells.csv.gz"
    df.to_csv(out_detail, index=False)
    log(f"写出 {os.path.relpath(out_detail, ROOT)}（{len(df):,} 行）")

    # ---- 待人工判清单（用户 2026-09-22 决定：两线打架的不自动裁）----
    un = df[df["state"] == "unadopted"].copy() if len(df) else df
    out_un = f"{OUT}/gp8c_readopt_unadopted.csv"
    un.to_csv(out_un, index=False)
    log(f"写出 {os.path.relpath(out_un, ROOT)}（{len(un):,} 行待人工判）")

    # ---- 全局一层：每细胞恰好一个谱系，且与 A_adjudicated 可逐行对照 ----
    g = ref[["cell_barcode", "A_frozen", "A_adjudicated"]].copy()
    g["readopt_state"] = "assigned"
    g["readopt_target"] = g["A_adjudicated"]
    if len(df):
        m = df.set_index("cell_barcode")
        hit = g["cell_barcode"].isin(m.index)
        g.loc[hit, "readopt_target"] = m.loc[g.loc[hit, "cell_barcode"], "A_readopted"].to_numpy()
        g.loc[hit, "readopt_state"] = m.loc[g.loc[hit, "cell_barcode"], "state"].to_numpy()
    g.rename(columns={"readopt_target": "A_readopted"}, inplace=True)
    if g["A_readopted"].isna().any():
        raise SystemExit("🔴 有细胞没有最终谱系 —— 违反'每细胞恰好一个谱系'")
    n_move = int((df["state"] == "readopted").sum()) if len(df) else 0
    if int((g["A_readopted"] != g["A_adjudicated"]).sum()) != n_move:
        raise SystemExit("🔴 移动的细胞数 ≠ readopted 条数 —— 内部不一致，停")
    out_g = f"{OUT}/gp8c_cell_assignment.csv.gz"
    g.to_csv(out_g, index=False)

    print("\n=== 接手后各谱系细胞数（对照 A_adjudicated）===")
    cmp = pd.DataFrame({
        "A_adjudicated": g["A_adjudicated"].value_counts(),
        "A_readopted": g["A_readopted"].value_counts()}).fillna(0).astype(int)
    cmp["净变化"] = cmp["A_readopted"] - cmp["A_adjudicated"]
    print(cmp.to_string())
    assert int(cmp["A_readopted"].sum()) == len(g), "谱系合计 ≠ 总细胞"
    print(f"\n合计 {int(cmp['A_readopted'].sum()):,} = 全部细胞 {len(g):,} ✅（互斥且完备）")

    mf = dict(
        script=os.path.relpath(os.path.abspath(__file__), ROOT),
        script_sha256=sha256(os.path.abspath(__file__)),
        created_at=datetime.datetime.now().isoformat(timespec="seconds"),
        caliber=("逐细胞三态接手：面板线 = gp6_scores.csv.gz 六谱系逐细胞模块分 argmax"
                 "（ctrl_size=50, seed=0，与 GP6 标准 A 同一份打分）；"
                 "CellTypist 线 = B_lineage。两线一致且 ≠ 原谱系 ⇒ readopted；"
                 "两线一致且 = 原谱系 ⇒ restored；不一致或任一线弃权 ⇒ unadopted。"),
        no_new_parameters=True,
        no_new_parameters_note=("本步没有引入任何新阈值：两条线及其参数完全沿用已登记的"
                                "GP6 打分与已存的 CellTypist 输出，只把聚合层级由簇降为细胞。"),
        unadopted_provisional_rule=("unadopted 细胞**暂不移动**，留在 A_adjudicated 给它的谱系。"
                                    "剔除与接手都是主动动作，人工未判之前不动最保守；"
                                    "故 A_adjudicated 与已签字的六份子集清单一字未改。"),
        disclosure=(("🔴 判据与剔除判据同源（同一份面板打分、同一份 CellTypist 输出），"
                     "是同一证据换聚合层级，**不是独立验证**。readopted 不可当独立第二证据用；"
                     "unadopted 的条数本身是本次簇级裁决不确定性的度量。")),
        inputs=dict(
            sources=in_sha,
            excluded_cells_sha256=ex_sha,
        ),
        outputs=dict(
            detail=os.path.relpath(out_detail, ROOT), detail_sha256=sha256(out_detail),
            unadopted=os.path.relpath(out_un, ROOT), unadopted_sha256=sha256(out_un),
            assignment=os.path.relpath(out_g, ROOT), assignment_sha256=sha256(out_g)),
        per_tag=per_tag,
        summary=dict(
            n_excluded=int(len(df)),
            n_restored=int((df["state"] == "restored").sum()) if len(df) else 0,
            n_readopted=int((df["state"] == "readopted").sum()) if len(df) else 0,
            n_unadopted=int((df["state"] == "unadopted").sum()) if len(df) else 0,
            readopt_to=({k: int(v) for k, v in
                         df.loc[df["state"] == "readopted", "target_lineage"]
                         .value_counts().items()} if len(df) else {}),
        ),
    )
    out_mf = f"{OUT}/gp8c_readopt_manifest.json"
    with open(out_mf, "w", encoding="utf-8") as fh:
        json.dump(mf, fh, ensure_ascii=False, indent=2)
    log(f"写出 {os.path.relpath(out_mf, ROOT)}")
    print(f"\n接手小结：被剔 {mf['summary']['n_excluded']:,} 核 → "
          f"还原 {mf['summary']['n_restored']:,} / 接手 {mf['summary']['n_readopted']:,} / "
          f"待人工判 {mf['summary']['n_unadopted']:,}")
    print("⚠️ 待人工判清单在 gp8c_readopt_unadopted.csv，本脚本**不替你做这个决定**")


if __name__ == "__main__":
    main()
