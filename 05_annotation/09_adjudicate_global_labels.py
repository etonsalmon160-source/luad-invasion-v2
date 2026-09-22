"""全局簇错标的**裁决层**：新增一列 `A_adjudicated`，`A_frozen` 一字不改。

为什么要有这一步
----------------
2026-09-22 调查「B/浆 谱系为什么剔掉 28.55%」时发现：GP6 的 A 标准是**逐簇 argmax**
（簇内 module score 均值最大的谱系），一个簇整体贴一个标签、再传播给簇内每个细胞。
当某个簇**本谱系的打分接近 0** 时，六个分数会一起塌平，argmax 就由噪声决定 ——
于是**整簇**被贴错，而下游是按簇切子集的，错一个簇就是几千个细胞走错门。

已确认 4 个这样的簇（原始簇号，r*=0.6 / seed0）：

  raw 17 (6,574 核)  A=上皮 → 真身**肥大细胞**  本簇 meandiff top15 = CPA3;HDC;KIT;MS4A2;
                                                VWA5A;SLC18A2;GATA2;IL1RL1;ALOX5;SLC45A3;
                                                RGS13;TPSG1;ATP6V0A2;SIGLEC6;HPGDS（全是肥大细胞
                                                特异基因）；CellTypist 髓系 97%；A 面板六分数
                                                全塌平（上皮+0.271 / 髓系+0.061 / T/NK−0.175）
  raw 26 (3,103 核)  A=B/浆 → 真身**DC**      本簇 meandiff top15 = IRF8;WDFY4;CD74;CIITA;
                                                MPEG1;LSP1;CSF2RA;LCP1;CST3;SPI1;…（IRF8/WDFY4/
                                                CIITA/CSF2RA/SPI1 是 DC 定义基因）；CellTypist 髓系 93%
  raw 38 (1,147 核)  A=上皮 → 真身**T/NK**    本簇落在表面活性剂环境 RNA 很浓的地方。逐基因检出率
                                                对深度匹配背景的倍数：T/NK 核心 2.50×、分泌型上皮
                                                1.79×、结构型上皮仅 1.15×。A 面板给本簇上皮
                                                +1.115（高于 T/NK +0.434）**正是分泌型那 1.79× 撑的**
  raw 44   (404 核)  A=内皮 → 真身**T/NK**    本簇 meandiff top15 = IL7R;TRAC;EVL;AKNA;TMC8;
                                                CD2;TRBC2;NLRC5;CORO1A;CD96;ACAP1;ARHGAP4;ITGAL;
                                                IKZF1;CXCR4（全是 T 基因）；内皮 vs T/NK 面板分
                                                0.444 vs 0.399，只差 0.045 —— 近并列

候选簇**不是人工挑的**，由 GP6 里**预注册**的复核规则机械确定：

    逐簇两线一致率 < 0.90  **且**  A_lineage ≠ B_lineage_mode（谱系级分歧）

2026-09-22 实跑：一致率 <0.90 的簇有 23 个（174,395 核），其中**谱系级**分歧的
**恰好 4 个**，即上表。本脚本**硬断言**这一等式 —— 人工裁决表与规则输出对不上就报错。
（另外 19 个只在**同一谱系内部**的亚型层面分歧，不在本次范围。）

环境 RNA 判别（为什么簇 38 不能只看面板分）
------------------------------------------
上皮面板 16 个基因里 **7 个是分泌型**（SFTPC/SFTPA1/SFTPA2/SFTPB/NAPSA/SCGB3A2/SCGB1A1
—— 表面活性剂与分泌球蛋白，分泌到胞外、易成为游离 mRNA），另 9 个是**结构型/转录因子**
（EPCAM/KRT8/KRT18/KRT19/CDH1/NKX2-1/AGER/CAV1/PDPN）。所以「上皮分高」不等于「这个核在
转录上皮基因」—— 泡在表面活性剂里的 T 细胞也会拿到高分。

判据（**不设阈值，看形状**）：把候选簇里两类基因的**检出率**和**深度匹配背景**（全图谱中
nCount 落在本簇 IQR 内的细胞，剔除本簇自身）比。结构型检出率贴近背景 ⇒ 那个「上皮分」
是环境 RNA 撑起来的。**不做随机抽样**（法则 3：无随机数）。

⚠️ 时序披露（法则 3.2 要求，不掩饰）
------------------------------------
本次裁决**发生在看到下游结果之后**：触发点是 B/浆 剔除率异常的调查，簇 26 由该调查
直接暴露；簇 17/38/44 是在其后**系统套用**上述预注册规则时发现的。规则本身预注册于
GP6、与本次调查无关，候选集由规则机械确定、非人工挑选；但「先看到下游、后做裁决」
这一时序**必须留在案卷里**，见 manifest 的 `timing_disclosure`。

本脚本**只新增列，不改 `A_frozen`**：旧口径原样保留，新旧两版并存可对账。
下游若要改用新口径，须先在 manifest 里签字（`signed_by`）。

基因集出处
----------
结构型/分泌型上皮、T/NK 核心三组**全部取自 `marker_panel.py` 已登记面板的子集**（源：
Travaglini2020 / Peng2026），不引入新出处。逐簇的 meandiff/wilcoxon/A_scores/B 判定
一律**取自 GP6 自己的冻结产物** `gp6_cluster_labels.csv`，不重算。

输出
----
  results/05_annotation/gp6_cell_labels_adjudicated.csv.gz  逐细胞：原 14 列 + 裁决 4 列
  results/05_annotation/gp6_adjudication_evidence.csv       逐簇：23 个低一致率簇全表 + 裁决标注
  results/05_annotation/gp6_adjudication_ambient_check.csv  逐簇×基因集：环境 RNA 判别
  results/05_annotation/gp6_adjudication_manifest.json      规则原文 + 证据 + 时序披露 + 哈希

🔴 manifest 的 `signed_by` 初始为 null，须人工签字后才能被下游使用。

用法:
    python3 05_annotation/09_adjudicate_global_labels.py
"""

import hashlib
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import marker_panel as MP  # noqa: E402

ROOT = "/home/eto/luad_v2"
OUT = f"{ROOT}/results/05_annotation"
LABELS = f"{OUT}/gp6_cell_labels.csv.gz"
CLUSTER_LABELS = f"{OUT}/gp6_cluster_labels.csv"
FULL_CLU = f"{ROOT}/results/04_integration/seurat_trad/full/clusters.csv.gz"
FULL_COL = "harmony_res0.6_seed0"
H5 = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"

# GP6 预注册的「强制人工复核」阈值。**不新发明**（法则 3.1）。
AGREEMENT_TH = 0.90

# ---- 人工裁决表 --------------------------------------------------------------
# to 一律取 A 面板六谱系之一；evidence 只写**从 GP6 冻结产物或本脚本重算得出**的事实。
# 键是**原始簇号**（字符串），不是 gp6_cluster_labels.csv 里的分类编码。
DECISIONS = {
    "17": dict(to="髓系", phenotype="肥大细胞",
               evidence="driver_genes_meandiff top15 全为肥大细胞特异基因"
                        "（CPA3;HDC;KIT;MS4A2;GATA2;TPSG1;IL1RL1;ALOX5;HPGDS;SIGLEC6 等）；"
                        "CellTypist 髓系占 97%；A 面板六分数全塌平（上皮+0.271/髓系+0.061/T-NK−0.175）",
               src_panel="GATA2/KIT 属髓系范畴，但 marker_panel.py 的髓系面板未收肥大细胞 marker"),
    "26": dict(to="髓系", phenotype="树突状细胞（DC）",
               evidence="driver_genes_meandiff top15 含 IRF8;WDFY4;CIITA;CSF2RA;SPI1;CD74 —— DC 定义基因；"
                        "CellTypist 髓系占 93%；A 面板 B/浆仅 +0.282 且上皮 +0.136，六分数塌平",
               src_panel="marker_panel.py 的髓系面板未收 DC marker（无 IRF8/WDFY4/CIITA）"),
    "38": dict(to="T/NK", phenotype="T 细胞（泡在表面活性剂环境 RNA 里）",
               evidence="环境 RNA 判别（逐基因检出率对深度匹配背景的倍数）："
                        "T/NK 核心 2.50×（三组最高）、分泌型上皮 1.79×、结构型上皮仅 1.15×。"
                        "A 面板给本簇上皮 +1.115（明显高于 T/NK +0.434）**正是因为这 1.79×**——"
                        "上皮面板 7/16 是分泌型基因，会随环境 RNA 一起被抬高；"
                        "结构型上皮只比背景高 15%，撑不起「上皮」这个身份。CellTypist 族亦判 T/NK。",
               src_panel="上皮面板 7/16 为分泌型，是本簇被误判的机制"),
    "44": dict(to="T/NK", phenotype="T/NK",
               evidence="driver_genes_meandiff top15 全为 T 基因"
                        "（IL7R;TRAC;CD2;TRBC2;CORO1A;IKZF1;CD96;ITGAL;CXCR4 等）；"
                        "A 面板内皮 +0.444 vs T/NK +0.399，只差 0.045（近并列）",
               src_panel=""),
}

# ---- 环境 RNA 判别用的基因集：全部是**已登记面板的子集**，不引入新出处 ----
_SECRET_EPI = ["SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "NAPSA", "SCGB3A2", "SCGB1A1"]
_EPI_PANEL = list(MP.PANEL["上皮"]["genes"])
STRUCT_EPI = [g for g in _EPI_PANEL if g not in _SECRET_EPI]
assert set(STRUCT_EPI) | set(_SECRET_EPI) == set(_EPI_PANEL), "结构型/分泌型拆分没覆盖满上皮面板"
TNK_CORE = list(MP.PANEL["T/NK"]["genes"])

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def sha256(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def read_str_col(node, key="_index"):
    key = key if key in node else list(node.keys())[0]
    v = node[key][:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def candidate_clusters(cl):
    """预注册复核规则 ⇒ 候选簇（原始簇号）。返回 (全部低一致率簇, 谱系级分歧簇)。"""
    low = cl[cl["B_agreement_rate"] < AGREEMENT_TH].copy()
    lin = low[low["A_lineage"] != low["B_lineage_mode"]].copy()
    return low, lin


def main():
    # ---------- 1. 簇号解码：分类编码 → 原始簇号 ----------
    cl = pd.read_csv(CLUSTER_LABELS)
    full = pd.read_csv(FULL_CLU)
    if FULL_COL not in full.columns:
        raise SystemExit(f"🔴 {FULL_COL} 不在 {os.path.relpath(FULL_CLU, ROOT)} 的列里")
    cats = list(pd.Categorical(full[FULL_COL].astype(str)).categories)
    cl["raw"] = [cats[int(c)] for c in cl["cluster"]]
    vc = full[FULL_COL].astype(str).value_counts()
    bad = [(r, int(n), int(vc.get(r, -1))) for r, n in zip(cl["raw"], cl["n_cells"])
           if int(vc.get(r, -1)) != int(n)]
    if bad:
        raise SystemExit(f"🔴 簇号解码自检失败（raw,n_gp6,n_真值）：{bad[:5]}")
    log(f"簇号解码 ✅ {len(cats)} 簇，逐簇 n_cells 与 {FULL_COL} 完全一致")

    # ---------- 2. 候选簇由预注册规则机械确定 ----------
    low, lin = candidate_clusters(cl)
    cand = sorted(lin["raw"].tolist(), key=int)
    log(f"规则：一致率 <{AGREEMENT_TH} 的簇 {len(low)} 个（{low['n_cells'].sum():,} 核）；"
        f"其中谱系级分歧 {len(cand)} 个 → {cand}")
    if set(cand) != set(DECISIONS):
        raise SystemExit(f"🔴 人工裁决表 {sorted(DECISIONS, key=int)} ≠ 规则候选 "
                         f"{cand} —— 决不允许人工挑簇。停。")

    # ---------- 3. 逐细胞：新增裁决列，A_frozen 原样不动 ----------
    lab = pd.read_csv(LABELS)
    if lab["cell_barcode"].duplicated().any():
        raise SystemExit("🔴 gp6_cell_labels 里 barcode 有重复")
    clust_of = dict(zip(full["cell_barcode"], full[FULL_COL].astype(str)))
    lab["_raw"] = lab["cell_barcode"].map(clust_of)
    if lab["_raw"].isna().any():
        raise SystemExit(f"🔴 有 {int(lab['_raw'].isna().sum())} 个细胞不在 {FULL_COL} 里")
    # 自检：逐簇 A_frozen 必须唯一，且等于 gp6_cluster_labels 的 A_lineage
    x = lab.groupby("_raw")["A_frozen"].nunique()
    if (x > 1).any():
        raise SystemExit(f"🔴 这些簇里 A_frozen 不唯一：{x[x > 1].index.tolist()}")
    a_of = dict(zip(cl["raw"], cl["A_lineage"]))
    mm = lab.drop_duplicates("_raw").set_index("_raw")["A_frozen"]
    diff = [r for r in mm.index if mm[r] != a_of.get(r)]
    if diff:
        raise SystemExit(f"🔴 逐细胞与逐簇的 A 标签不符：{diff[:5]}")

    lab["A_adjudicated"] = lab["A_frozen"]
    lab["adjudicated"] = False
    lab["adj_from"] = ""
    lab["adj_reason"] = ""
    for r, dec in DECISIONS.items():
        m = lab["_raw"] == r
        if not m.any():
            raise SystemExit(f"🔴 裁决表里的簇 {r} 在逐细胞表里找不到")
        lab.loc[m, "A_adjudicated"] = dec["to"]
        lab.loc[m, "adjudicated"] = True
        lab.loc[m, "adj_from"] = a_of[r]
        lab.loc[m, "adj_reason"] = f"raw{r}:{a_of[r]}→{dec['to']}({dec['phenotype']})"
        log(f"  raw {r:>2s}  {a_of[r]:>3s} → {dec['to']:<3s}  {int(m.sum()):,} 核")
    lab = lab.drop(columns=["_raw"])

    out_cells = f"{OUT}/gp6_cell_labels_adjudicated.csv.gz"
    lab.to_csv(out_cells, index=False, compression="gzip")
    log(f"写出 {os.path.relpath(out_cells, ROOT)}")

    print("\n改变后的谱系分布（A_frozen → A_adjudicated）")
    cmp = pd.crosstab(lab["A_frozen"], lab["A_adjudicated"])
    print(cmp.to_string())
    print()
    old = lab["A_frozen"].value_counts()
    new = lab["A_adjudicated"].value_counts()
    for lin_name in MP.LINEAGES:
        o, n = int(old.get(lin_name, 0)), int(new.get(lin_name, 0))
        print(f"  {lin_name:4s} {o:7,d} → {n:7,d}  ({n - o:+,d})")

    # ---------- 4. 逐簇证据表：23 个低一致率簇全表 + 裁决标注 ----------
    ev = low[["raw", "n_cells", "A_lineage", "B_lineage_mode", "B_agreement_rate",
              "A_margin", "A_scores", "driver_genes_meandiff", "driver_genes_wilcoxon",
              "tri_agree"]].copy()
    ev["lineage_level_disagreement"] = ev["raw"].isin(cand)
    ev["adjudicated_to"] = ev["raw"].map({r: d["to"] for r, d in DECISIONS.items()}).fillna("")
    ev["adjudication_evidence"] = ev["raw"].map({r: d["evidence"] for r, d in DECISIONS.items()}).fillna("")
    ev = ev.sort_values(["lineage_level_disagreement", "B_agreement_rate"],
                        ascending=[False, True])
    out_ev = f"{OUT}/gp6_adjudication_evidence.csv"
    ev.to_csv(out_ev, index=False)
    log(f"写出 {os.path.relpath(out_ev, ROOT)}（{len(ev)} 簇）")

    # ---------- 5. 环境 RNA 判别 ----------
    cand_set = set(cand)
    amb, amb_sum, amb_missing = ambient_check(full, cand_set)

    # ---------- 6. manifest ----------
    out_man = f"{OUT}/gp6_adjudication_manifest.json"
    json.dump(dict(
        stage="GP6 全局簇错标裁决层",
        date="2026-09-22",
        caliber=("在 gp6_cell_labels.csv.gz 的逐细胞标签上**新增** A_adjudicated 一列；"
                 "A_frozen 逐字节保留，旧口径不作废，新旧两版并存可对账。"),
        rule=dict(
            text="逐簇两线一致率 < %.2f 且 A_lineage != B_lineage_mode（谱系级分歧）" % AGREEMENT_TH,
            threshold=AGREEMENT_TH,
            threshold_source="GP6 预注册的「强制人工复核」阈值，本次未新发明、未调整",
            n_low_agreement=int(len(low)),
            n_low_agreement_cells=int(low["n_cells"].sum()),
            n_lineage_level=int(len(cand)),
            candidates=cand,
            note="候选集由规则机械确定；脚本硬断言人工裁决表 == 规则输出，否则 SystemExit。",
        ),
        decisions={r: dict(from_=a_of[r], to=d["to"], phenotype=d["phenotype"],
                           n_cells=int(low.loc[low["raw"] == r, "n_cells"].iloc[0]),
                           evidence=d["evidence"], panel_blind_spot=d["src_panel"])
                   for r, d in DECISIONS.items()},
        timing_disclosure=(
            "🔴 本裁决**发生在看到下游结果之后**。触发点：2026-09-22 调查「B/浆 谱系剔除率 "
            "28.55% 是否合理」，簇 26 由该调查直接暴露；簇 17/38/44 是在其后**系统套用**上述"
            "预注册规则（遍历全部 45 个全局簇）时发现的。规则本身预注册于 GP6、与本次调查无关，"
            "候选集由规则机械确定、非人工挑选。但「先看到下游、后做裁决」这一时序如实登记，"
            "不掩饰（法则 3.2）。"),
        known_limitation=(
            "marker_panel.py 的髓系面板（16 基因）**不含任何肥大细胞 marker**"
            "（无 CPA3/TPSAB1/TPSG1/MS4A2/KIT/HDC/GATA2）、**不含任何 DC marker**"
            "（无 IRF8/WDFY4/CIITA/CSF2RA/LAMP3/CLEC9A/IL3RA）；"
            "上皮面板 16 基因中 **7 个为分泌型**（表面活性剂/分泌球蛋白），易受环境 RNA 抬高。"
            "本次经用户 2026-09-22 裁定**不补面板**，作为**已知局限**登记；"
            "意味着自动化的簇级 argmax 在这两类细胞上**仍会出错**，须靠本裁决层兜底。"),
        ambient_check=dict(
            rule="结构型上皮检出率 vs 深度匹配背景（全图谱 nCount 落在本簇 IQR 内的细胞，"
                 "剔除本簇自身）；不设阈值，看形状；不使用随机数。",
            struct_epi=STRUCT_EPI, secret_epi=_SECRET_EPI, tnk_core=TNK_CORE,
            source="三组基因集全部取自 marker_panel.py 已登记面板的子集，未引入新出处",
            gene_sets_missing_in_matrix=amb_missing,
            summary=amb_sum.to_dict(orient="records"),
            csv=os.path.relpath(f"{OUT}/gp6_adjudication_ambient_check.csv", ROOT),
        ),
        signed_by=None,
        signed_note="🔴 待人工签字：核对 gp6_adjudication_evidence.csv 与 "
                    "gp6_adjudication_ambient_check.csv 后，把 signed_by 改为署名"
                    "（如 \"<姓名> 2026-09-22\"）。未签字前下游不得改用 A_adjudicated。",
        outputs={os.path.relpath(p, ROOT): sha256(p)
                 for p in [out_cells, out_ev, f"{OUT}/gp6_adjudication_ambient_check.csv"]},
        inputs={os.path.relpath(p, ROOT): sha256(p)
                for p in [LABELS, CLUSTER_LABELS, FULL_CLU, H5, os.path.abspath(__file__),
                          os.path.join(os.path.dirname(os.path.abspath(__file__)), "marker_panel.py")]},
    ), open(out_man, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log(f"写出 {os.path.relpath(out_man, ROOT)}  ——  🔴 须人工签字")
    print(f"\n{'#' * 100}\n# 环境 RNA 判别逐基因明细（簇 38 的决定性证据；四个候选簇一律同口径报出）\n{'#' * 100}")
    print(amb.to_string(index=False))


def ambient_check(full, cand_set):
    """四个候选簇 × 基因集：检出率 vs 深度匹配背景。按列取基因，不加载整个 h5ad。"""
    genes = STRUCT_EPI + _SECRET_EPI + TNK_CORE
    f = h5py.File(H5, "r")
    gene_names = read_str_col(f["var"])
    g2j = {g: i for i, g in enumerate(gene_names)}
    obs_bc = read_str_col(f["obs"])
    ncount = f["obs"]["nCount"][:].astype(np.float64)
    missing = [g for g in genes if g not in g2j]
    if missing:
        log(f"⚠️ 矩阵里没有、已跳过（如实登记）：{missing}")
    genes = [g for g in genes if g in g2j]
    tjs = np.array(sorted({g2j[g] for g in genes}), dtype=np.int32)

    log("读 X 的 indices/data（只扫一遍列）…")
    indptr = f["X"]["indptr"][:]
    indices = f["X"]["indices"][:]
    data = f["X"]["data"][:]
    f.close()
    sel = np.isin(indices, tjs)
    hit_pos = np.flatnonzero(sel)
    hit_gene = indices[hit_pos].astype(np.int64)
    hit_val = data[hit_pos].astype(np.float64)
    hit_cell = np.searchsorted(indptr, hit_pos, side="right") - 1
    del sel, hit_pos, indices, data
    log(f"命中 {len(hit_gene):,} 个非零元")

    if not np.array_equal(obs_bc, full["cell_barcode"].to_numpy()):
        order = pd.Index(obs_bc)
        ci = order.get_indexer(full["cell_barcode"].to_numpy())
        if (ci < 0).any():
            raise SystemExit("🔴 有细胞不在矩阵里")
        raw_of = np.full(len(obs_bc), "", dtype=object)
        raw_of[ci] = full[FULL_COL].astype(str).to_numpy()
        nc = np.empty(len(obs_bc))
        nc[ci] = ncount
        ncount = nc
    else:
        raw_of = full[FULL_COL].astype(str).to_numpy()

    isn = ncount > 0
    rows = []
    for r in sorted(cand_set, key=int):
        m = raw_of == r
        if not m.any():
            raise SystemExit(f"🔴 候选簇 {r} 在矩阵里找不到细胞")
        q1, q3 = np.percentile(ncount[m], [25, 75])
        bg = (~m) & (ncount >= q1) & (ncount <= q3) & isn
        n_in, n_bg = int(m.sum()), int(bg.sum())
        if n_bg < 500:
            raise SystemExit(f"🔴 簇 {r} 的深度匹配背景只有 {n_bg} 个细胞，不足以作对照")
        for gset, gname in ((STRUCT_EPI, "结构型上皮"), (_SECRET_EPI, "分泌型上皮(易环境RNA)"),
                            (TNK_CORE, "T/NK 核心")):
            for g in gset:
                if g not in g2j:
                    continue
                s = hit_gene == g2j[g]
                c = hit_cell[s]
                p_in = float(m[c].sum()) / n_in
                p_bg = float(bg[c].sum()) / n_bg
                rows.append(dict(
                    raw=r, n_cells=n_in, gene_set=gname, gene=g,
                    pct_in=round(p_in, 4), pct_bg=round(p_bg, 4),
                    ratio=round(p_in / p_bg, 2) if p_bg > 0 else np.inf,
                    n_bg=n_bg, bg_nCount_IQR=f"{q1:.0f}~{q3:.0f}",
                ))
    d = pd.DataFrame(rows)
    gs = (d.groupby(["raw", "gene_set"], as_index=False)
            .agg(n_genes=("gene", "size"),
                 pct_in_mean=("pct_in", "mean"),
                 pct_bg_mean=("pct_bg", "mean")))
    gs["pct_in_mean"] = gs["pct_in_mean"].round(4)
    gs["pct_bg_mean"] = gs["pct_bg_mean"].round(4)
    gs["ratio_mean"] = (gs["pct_in_mean"] / gs["pct_bg_mean"]).round(2)
    gs["n_genes_over_2x"] = [
        int((d[(d["raw"] == r) & (d["gene_set"] == s)]["ratio"] > 2).sum())
        for r, s in zip(gs["raw"], gs["gene_set"])]
    d.to_csv(f"{OUT}/gp6_adjudication_ambient_check.csv", index=False)
    log(f"写出 {os.path.relpath(f'{OUT}/gp6_adjudication_ambient_check.csv', ROOT)}")
    print()
    print("# 逐簇×基因集汇总（ratio ≈1 ⇒ 与背景无异，即环境 RNA；ratio ≫1 ⇒ 本簇真在转录）")
    print(gs.to_string(index=False))
    print()
    return d, gs, missing


if __name__ == "__main__":
    main()
