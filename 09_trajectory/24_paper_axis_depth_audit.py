#!/usr/bin/env python
"""
事后（未注册）诊断：K1 报了"KAC 侧簇系统性更深"（2198.75 vs 1580.0），
按项目惯例（见 `feedback_negation_needs_assumption_audit`）——**下"这只是深度"的否定结论前，
必须先审计测法的隐含假设**，否则会自信地报错。

本脚本**不改任何判据、不重跑 D1/D2**，只回答一个问题：
  **KAC 侧那 8 个簇的 MP6-argmax，是不是深度撑起来的？**

三段，分档标准事前写死：

  A1 偏相关：27 簇上 MP6 vs KAC 签名，**扣掉 nFeature_median** 后的偏 Spearman。
             （D2 报的 +0.722 若主要是深度，偏相关应塌到 0 附近）

  A2 分档复判：**事前选定**全局 nFeature 的 [q25, q75] 为共同深度档。
             在只留档内细胞后，重算每簇 8 型 MP 的簇级均值与 argmax。
             若 8 个簇仍以 MP6 为 argmax ⇒ argmax 不是深度撑的（argmax 是簇内相对比较，
             深度对 8 个 MP 大致同向抬升）。逐簇报档内细胞数，n 太小的簇标出来。

  A3 谁最吃深度：8 个 MP 各自与 nFeature_median 的簇级 Spearman。
             若 MP6 的深度敏感性**远超** MP2(AT2) 等 ⇒ D1 的 argmax 可能被深度带偏。

边界：本脚本只产诊断，**不产任何判据**；结论只用于判断 D1/D2 该报"成立"还是"待定"。
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

D = "/home/eto/luad_v2/results/09_trajectory/paper_axis"
TSV = f"{D}/per_cell_scores.tsv.gz"
OUT = f"{D}/depth_audit.tsv"
WT = "/home/eto/luad_v2/.claude/worktrees/vigilant-hawking-799963"

KAC_SIDE = {"Tumor cell/KAC", "KAC/inflammatory"}

sys.path.insert(0, f"{WT}/05_annotation")
import epi_subtype_panel as P                                  # noqa: E402


# 列名 → 真亚型名：**必须**按名字对，不能按位置、更不能靠列名反向还原
# （`colname()` 把非字母数字换成 `_`，"Tumor cell/KAC" 的 `/` 会丢 ⇒ 反向还原造过假翻转；
#   按位置对则列序一变就静默错配，且 assert 只数个数、看不出来）
def _colname(subtype):
    return "s_" + "".join(c if c.isalnum() else "_" for c in subtype)


_PANEL_MP = list(P.build_panel("paper_mp").keys())              # 与打分时的插入顺序一致
COL2TYPE = {_colname(t): t for t in _PANEL_MP}

df = pd.read_csv(TSV, sep="\t")
mp_cols = [c for c in df.columns if c.startswith("s_") and c != "s_KACsig"]
_missing = [c for c in mp_cols if c not in COL2TYPE]
if _missing:
    raise SystemExit(f"🔴 这些列名在面板里找不到对应亚型，不敢按位置猜：{_missing}")
mp_types = [COL2TYPE[c] for c in mp_cols]
ck = "s_KACsig"

cm = df.groupby("cluster").agg(
    n_cells=("cluster", "size"),
    nFeature_median=("nFeature", "median"),
    **{c: (c, "mean") for c in mp_cols + [ck]},
).reset_index()
print(f"细胞 {len(df):,} / 簇 {len(cm)}")

prev = pd.read_csv(f"{D}/cluster_scores_papermp.csv")[["cluster", "is_kac_side", "argmax_papermp"]]
cm = cm.merge(prev, on="cluster", how="left")

# ---------------- A1 偏相关 ----------------
def partial_spearman(x, y, z):
    """扣掉 z 后的 x,y 偏 Spearman（先各自对 z 的秩回归取残差，再算 Pearson）。"""
    rx, ry, rz = (stats.rankdata(v).astype(float) for v in (x, y, z))
    Z = np.column_stack([np.ones_like(rz), rz])
    bx = np.linalg.lstsq(Z, rx, rcond=None)[0]
    by = np.linalg.lstsq(Z, ry, rcond=None)[0]
    ex, ey = rx - Z @ bx, ry - Z @ by
    return float(np.corrcoef(ex, ey)[0, 1])

rows = []
for t, c in zip(mp_types, mp_cols):
    raw = stats.spearmanr(cm[c], cm[ck])
    par = partial_spearman(cm[c].to_numpy(), cm[ck].to_numpy(), cm["nFeature_median"].to_numpy())
    rows.append(dict(mp=t, rho_raw=raw[0], p_raw=raw[1], rho_partial_ctrl_nFeature=par,
                     rho_mp_vs_nFeature=stats.spearmanr(cm[c], cm["nFeature_median"])[0]))
a3 = pd.DataFrame(rows).sort_values("rho_mp_vs_nFeature", ascending=False)

print("\n=== A3 各 MP 的簇级深度敏感性（按 |rho vs nFeature| 降序）===")
print(a3[["mp", "rho_mp_vs_nFeature", "rho_raw", "rho_partial_ctrl_nFeature"]]
      .to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
print("\n=== A1 焦点：MP6 与 KAC 签名 ===")
r6 = a3[a3.mp == "Tumor cell/KAC"].iloc[0]
print(f"  原始 rho        = {r6.rho_raw:+.3f}")
print(f"  扣 nFeature 后  = {r6.rho_partial_ctrl_nFeature:+.3f}")
r9 = a3[a3.mp == "KAC/inflammatory"].iloc[0]
print(f"  （对照 MP9：原始 {r9.rho_raw:+.3f} → 扣深度 {r9.rho_partial_ctrl_nFeature:+.3f}）")

# ---------------- A2 分档复判 ----------------
q25, q75 = np.percentile(df["nFeature"].to_numpy(), [25, 75])
band = df[(df["nFeature"] >= q25) & (df["nFeature"] <= q75)]
print(f"\n=== A2 共同深度档（事前定：全局 nFeature 的 [q25,q75] = [{q25:.0f}, {q75:.0f}]）===")
print(f"  档内细胞 {len(band):,} / {len(df):,} = {len(band)/len(df):.1%}")

bm = band.groupby("cluster").agg(n_band=("cluster", "size"),
                                 **{c: (c, "mean") for c in mp_cols}).reset_index()
B = bm[mp_cols].to_numpy()
bm["argmax_band"] = np.array(mp_types)[B.argmax(axis=1)]
bm["is_kac_side_band"] = [t in KAC_SIDE for t in bm["argmax_band"]]
j = cm[["cluster", "n_cells", "nFeature_median", "argmax_papermp", "is_kac_side"]].merge(bm, on="cluster")
j["argmax_flip"] = j.argmax_papermp != j.argmax_band
print("\n  逐簇：原 argmax vs 档内 argmax")
print(j[["cluster", "n_cells", "n_band", "nFeature_median", "argmax_papermp",
         "argmax_band", "argmax_flip", "is_kac_side", "is_kac_side_band"]]
      .sort_values(["is_kac_side", "cluster"], ascending=[False, True])
      .to_string(index=False))
n_keep = int(j.loc[j.is_kac_side, "is_kac_side_band"].sum())
print(f"\n  8 个 KAC 侧簇在档内仍判 KAC 侧：**{n_keep}/8**")
print(f"  全体 27 簇 argmax 翻转：{int(j.argmax_flip.sum())} 个")
print(f"  档内细胞数最少的 3 簇：{j.nsmallest(3, 'n_band')[['cluster','n_band']].to_dict('records')}")

j.to_csv(OUT, sep="\t", index=False)
print(f"\n写出 {OUT}")

# ---------------- B 探索（未注册）：KAC 沿期别 ----------------
print("\n=== B 探索（**未注册**，只作线索）：KAC 分数沿期别 ===")
order = ["Normal", "AAH", "AIS", "MIA", "IAC"]
b = df[df["stage"].isin(order)]
med = b.groupby("stage")[ck].median().reindex(order)
cnt = b.groupby("stage").size().reindex(order)
feat = b.groupby("stage")["nFeature"].median().reindex(order)
print("  期别        n           KACsig中位    nFeature中位")
for s in order:
    print(f"  {s:<6s} {cnt[s]:>8,}   {med[s]:+.4f}      {feat[s]:.0f}")
if len(med.dropna()) >= 3:
    rho, p = stats.spearmanr([order.index(s) for s in med.dropna().index], med.dropna().values)
    print(f"  期别序号 vs KACsig 中位：Spearman rho={rho:+.3f}  p={p:.3g}")

# 逐患者（防止合并假象）
print("\n  逐患者（KAC 分数 vs 期别序号，只收 ≥3 个期的患者）：")
prows = []
for pid, g in b.groupby("patient_id"):
    gm = g.groupby("stage")[ck].median()
    ix = [order.index(s) for s in gm.index if s in order]
    if len(ix) >= 3:
        rr, pp = stats.spearmanr(ix, gm.values)
        prows.append(dict(patient=pid, n_stages=len(ix), rho=rr, p=pp))
pr = pd.DataFrame(prows)
if len(pr):
    print(f"    {len(pr)} 例；斜率正 {(pr.rho>0).sum()} 例 / 负 {(pr.rho<0).sum()} 例；"
          f"中位 rho {pr.rho.median():+.3f}")
    print(pr.sort_values("rho", ascending=False).to_string(index=False, float_format=lambda v: f"{v:+.3f}"))

print("\n=== 簇 × 期别构成（KAC 侧 vs 其余）===")
comp = pd.crosstab(df["cluster"], df["stage"], normalize="index")[order] * 100
comp["is_kac_side"] = comp.index.map(cm.set_index("cluster")["is_kac_side"])
print(comp.round(1).sort_values(["is_kac_side"], ascending=False).to_string())
print("\n  分组均值（%）：")
print(comp.groupby("is_kac_side")[order].mean().round(1).to_string())
