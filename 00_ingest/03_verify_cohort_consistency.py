#!/usr/bin/env python3
"""00_ingest/03_verify_cohort_consistency.py

独立核验：冻结表 `results/00_ingest/paired_samples.csv` 的每一个
(patient_id, stage, lesion_ordinal) 都必须与 GEO 权威表**逐字相等**。

------------------------------------------------------------------------------
为什么要新增这个脚本（2026-09-16）

`cohort_registry.py` 曾同时存在**两份**"病灶序号"的实现：

  (a) `load_geo()`        —— 读 GEO 权威表。冻结实际走的这条路，**正确**。
  (b) `lesion_ordinal()`  —— 按 `-1` 后缀猜。**从未被任何代码调用**，且对
                             GSE308103 全错（`AAH1`/`AIS1`/… 被推成 1，真值为 2）。

后果：产物是对的，缺陷是"活的但不执行"。**单测、流水线、产物核验全都只覆盖
被执行过的代码**，所以没有任何一个环节会去读一个没人调用的函数 —— 它不产生症状，
也就无从被检出。函数 (b) 已删除，序号读取统一走 `resolve_lesion_ordinal()`。

本脚本是**防复发**的那一层：它不依赖任何推导，只把冻结表与 GEO 表对撞。
将来若有人再写一份"猜"的实现并接进冻结，这里会当场红。
------------------------------------------------------------------------------

用法：
    python3 00_ingest/03_verify_cohort_consistency.py            # 核验正式冻结表
    python3 00_ingest/03_verify_cohort_consistency.py <路径>     # 核验指定表
        （第二个用法是**证伪入口**：拿一份故意做坏的表来撞它，若不能 FAIL，
          说明这个守卫是假的。见本文件末尾的 __main__ 说明。）
"""
import csv
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
import cohort_registry as REG  # noqa: E402

FROZEN = os.path.join(_ROOT, "results", "00_ingest", "paired_samples.csv")


def main(frozen=FROZEN):
    FROZEN = frozen
    print("=" * 74)
    print("00_ingest/03_verify_cohort_consistency.py —— 冻结表 vs GEO 权威表")
    print("=" * 74)

    if not os.path.exists(FROZEN):
        print(f"[FAIL] 冻结表不存在：{FROZEN}")
        return 1

    with open(FROZEN, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    print(f"[info] 冻结表 {os.path.relpath(FROZEN, _ROOT)}：{len(rows)} 行")

    # GEO 权威：按 dataset 建 (patient_id_token -> 记录)
    geo_by_ds = {}
    for ds in REG.COHORTS:
        geo_by_ds[ds] = {f"{g['patient_id']}_{g['token']}": g
                         for g in REG.load_geo(ds).values()}
        print(f"[info] GEO {ds}：{len(geo_by_ds[ds])} 个样本")

    bad = []

    # ---- 1. 逐行对撞 ---------------------------------------------------------
    for r in rows:
        ds, sid = r["dataset"], r["sample_id"]
        if ds not in geo_by_ds:
            bad.append(f"{ds}/{sid}: 数据集不在册")
            continue
        g = geo_by_ds[ds].get(sid)
        if g is None:
            bad.append(f"{ds}/{sid}: 样本不在 GEO 权威表中（拒绝无据冻结）")
            continue

        # 序号：必须以 GEO 为准，逐字相等
        if int(r["lesion_ordinal"]) != g["lesion_ordinal"]:
            bad.append(f"{ds}/{sid}: lesion_ordinal 冻结={r['lesion_ordinal']} "
                       f"GEO={g['lesion_ordinal']}")
        # 患者
        if r["patient_id"] != g["patient_id"]:
            bad.append(f"{ds}/{sid}: patient_id 冻结={r['patient_id']} GEO={g['patient_id']}")
        # 分期：token 必须能被 resolve_stage 严格解析，且等于冻结的规范分期
        try:
            got = REG.resolve_stage(r["stage_token"])
        except KeyError as e:
            bad.append(f"{ds}/{sid}: stage_token {r['stage_token']!r} 无法解析（{e}）")
            continue
        if got != r["stage"]:
            bad.append(f"{ds}/{sid}: stage 冻结={r['stage']} "
                       f"resolve_stage({r['stage_token']!r})={got}")
        # 分期词的原始 GEO 值应与 token 一致（token 去掉尾部第二病灶标记）
        if not r["stage_token"].startswith(g["stage"]):
            bad.append(f"{ds}/{sid}: stage_token={r['stage_token']!r} "
                       f"与 GEO stage={g['stage']!r} 前缀不符")

    # ---- 2. 覆盖度：GEO 的样本是否都在冻结表里 ------------------------------
    frozen_keys = {(r["dataset"], r["sample_id"]) for r in rows}
    for ds, idx in geo_by_ds.items():
        missing = sorted(set(idx) - {s for d, s in frozen_keys if d == ds})
        if missing:
            print(f"[warn] {ds}: GEO 有 {len(missing)} 个样本不在冻结表内"
                  f"（若为有意排除，请确认）：{', '.join(missing[:8])}"
                  f"{' …' if len(missing) > 8 else ''}")

    # ---- 3. 让"两套命名约定"每次运行都可见（这正是原缺陷的藏身处）----------
    print("\n--- 第二病灶（lesion_ordinal=2）的 token 写法 ---")
    n2_frozen = sum(1 for r in rows if int(r["lesion_ordinal"]) == 2)
    n2_geo = sum(1 for idx in geo_by_ds.values()
                 for g in idx.values() if g["lesion_ordinal"] == 2)
    for ds, idx in geo_by_ds.items():
        toks = sorted(g["token"] for g in idx.values() if g["lesion_ordinal"] == 2)
        dashed = [t for t in toks if "-" in t]
        undashed = [t for t in toks if "-" not in t]
        print(f"  {ds}: 带横线 {len(dashed)} {dashed}  |  无横线 {len(undashed)} {undashed}")
    print(f"  lesion_ordinal=2 计数：冻结 {n2_frozen} vs GEO {n2_geo}")
    if n2_frozen != n2_geo:
        bad.append(f"lesion_ordinal=2 计数不符：冻结 {n2_frozen} vs GEO {n2_geo}")

    # 两套约定同时存在 → 任何"按后缀猜"的实现必然错一半。这里只报警不硬停：
    # 若将来队列真的只剩一套约定，本行会提示复核，而不是让脚本无故失败。
    kinds = {("dash" if "-" in g["token"] else "nodash")
             for idx in geo_by_ds.values() for g in idx.values()
             if g["lesion_ordinal"] == 2}
    if len(kinds) > 1:
        print("  ⚠️ 两套约定**同时存在** → 严禁由 token 猜序号；必须查 GEO 表。")

    # ---- 结论 ---------------------------------------------------------------
    print("\n" + "=" * 74)
    if bad:
        print(f"[FAIL] {len(bad)} 项不一致：")
        for b in bad:
            print(f"  ✗ {b}")
        return 1
    print("[OK] 冻结表与 GEO 权威表逐字一致（patient_id / stage / lesion_ordinal）。")
    return 0


if __name__ == "__main__":
    # 无参数 = 核验正式冻结表；带参数 = 证伪入口（拿故意做坏的表撞它）
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else FROZEN))
