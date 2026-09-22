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
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 亚型（L2）面板来源。默认 classic_panels（Table S1，保持不变）；
# --panel s4 走 2026-09-22 新登记的面板（Table S4，内皮/髓系两谱系）。
# 🔴 只影响**亚型**注释。六谱系上位归属（L1）由 marker_panel.py 决定，与此无关。
PANEL_MODULES = {"s1": "classic_panels", "s4": "classic_panels_s4"}
CP = None          # 由 main() 依 --panel 载入（保持下方各处的 CP.xxx 引用不变）
PANEL_PATH = None  # 实际使用的面板模块文件路径（写进 manifest 的 panel_sha256）
PANEL_NAME = "s1"  # 实际使用的面板名（写进 manifest 的 panel）

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
               # ⚠️ 2026-09-22 修：原来这里**写死**了"经用户 2026-09-17 裁定"——那是
               # 成纤维那次的日期。上皮若也放宽，产物会把成纤维的裁定日期冒充成自己的，
               # 属伪造溯源。改成只陈述"经用户裁定"，具体日期以 signed_by / 裁定记录为准。
               relaxed_note=("🔴 事后放宽预注册阈值，经用户裁定（裁定日期见 signed_by）。"
                             "该谱系结果只能作**探索性**结论，不得进主结论。"
                             if relax_seed_min is not None else None),
               signed_by=None,
               signed_note="🔴 待人工签字：确认 r* 无误后，把 signed_by 改为署名"
                           "（如 \"<姓名> 2026-09-18\"），`--stage annotate` 才会放行（§M3-A.5）。")
    p = f"{OUT}/{tag}_rstar.json"
    with open(p, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)
    log(f"写出 {os.path.relpath(p, ROOT)}  ——  🔴 须人工签字后才可进入 annotate")


# ==================================================== 亚型面板装载 + 面板签字闸门
def load_panel_module(panel):
    """依 --panel 载入面板模块，返回 (模块, 文件绝对路径)。"""
    if panel not in PANEL_MODULES:
        raise SystemExit(f"🔴 未登记的面板 {panel}；有 {list(PANEL_MODULES)}")
    mod = __import__(PANEL_MODULES[panel])
    path = os.path.abspath(mod.__file__)
    log(f"亚型面板 = {PANEL_MODULES[panel]}（{os.path.relpath(path, ROOT)}）")
    return mod, path


def check_panel_caliber(panel, module_path):
    """签字闸门：新面板（s4）未签字 ⇒ 不许注释。

    与 check_rstar_signature 同一精神（2026-09-22 加）：新面板的口径（来源、过筛规则、
    封顶、粒度、回退）属**新口径**，按法则 3.1/3.2 须**先冻结签字、再跑**。
    默认面板 s1 不受此闸门约束（已在 GP8a/GP8c 用过，行为不变）。
    """
    if panel == "s1":
        return None
    p = f"{OUT}/panel_caliber_{panel}.json"
    if not os.path.exists(p):
        raise SystemExit(f"🔴 没有 {os.path.relpath(p, ROOT)} —— {panel} 面板的口径尚未登记，拒绝注释")
    j = json.load(open(p, encoding="utf-8"))
    if not j.get("signed_by"):
        raise SystemExit(f"🔴 {os.path.relpath(p, ROOT)} 尚无 `signed_by` 签字字段 —— "
                         f"新口径须先签字、再跑（法则 3.1/3.2）。停止。")
    now = sha256(module_path)
    # 登记文件由 build_panels_s4.py 产出，契约是 `panel_module.sha256`（单一事实来源）。
    reg = (j.get("panel_module") or {}).get("sha256")
    if reg != now:
        raise SystemExit(f"🔴 面板模块与登记不符：登记时 {reg}，"
                         f"现在 {now} —— 签字签的是**具体那一份面板**，改动即失效。停止。")
    log(f"面板签字闸门 ✅ {panel}：signed_by = {j['signed_by']}"
        + ("（⚠️ 非手写签署，结果待追认）" if (j.get("sign_scope") or {}).get("countersign_required")
           else ""))
    return j


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


def enforce_annotation_seed(tag, seed):
    """注释用哪个种子，由 `05_seed_representativeness.py` 的**预注册判据**决定，此处只执行。

    背景（GP5_report.md §十二 登记的缺口）：原先「种子 0 有代表性」只在**全量**上验过，
    六个谱系未验。2026-09-21 补验后，三个谱系在各自 r* 处种子 0 掉队超过 0.01
    （髓系、B/浆、成纤维）⇒ 按判据须改用该谱系内最具代表性的种子。

    判据（跑之前写下，看结果前后未改，**不得在此处重新裁定**）：
        gap = max_s(mean_ari_to_other_seeds) − mean_ari_to_other_seeds[seed 0]
        gap ≤ 0.01（复用 §M3-A.3 破平口径，不新发明阈值）⇒ 并列，沿用种子 0
        gap >  0.01 ⇒ 改用 argmax 的那个种子

    种子只影响**用哪一列簇标签**（`clusters.csv.gz` 里 `harmony_res*_seed*`），
    不影响 r* —— r* 由**跨全部种子的平均** ARI 选出，与用哪个种子无关。
    """
    mf = f"{OUT}/seed_representativeness_manifest.json"
    if not os.path.exists(mf):
        raise SystemExit(f"🔴 没有 {os.path.relpath(mf, ROOT)} —— 注释的种子选择必须由"
                         f" 05_seed_representativeness.py 的预注册判据产生。"
                         f"拒绝静默默认用种子 0。")
    man = json.load(open(mf, encoding="utf-8"))
    per = man.get("verdict", {}).get("per_object", {})
    if tag not in per:
        raise SystemExit(f"🔴 {os.path.relpath(mf, ROOT)} 的 verdict.per_object 里没有 {tag}")
    po = per[tag]
    want = 0 if po["pass_"] else int(po["best_seed"])
    if int(seed) != want:
        raise SystemExit(
            f"🔴 {tag}：--seed {seed} ≠ 预注册判据决定的种子 {want} —— "
            f"（r*={po['r_star']}，种子0 均值 {po['ref_mean']:.4f}，最好种子 "
            f"{int(po['best_seed'])} 均值 {po['best_mean']:.4f}，gap {po['gap']:+.4f}，"
            f"种子0 排名 {int(po['ref_rank'])}/5）。拒绝注释。")
    log(f"种子闸门 ✅ {tag}：用种子 {want}（预注册判据；"
        + ("并列，沿用种子 0" if po["pass_"]
           else f"种子 0 掉队 gap={po['gap']:+.4f} ⇒ 改用最好种子 {want}") + "）")
    return want, dict(seed=int(want), rule=man["preregistered_rule"],
                      evidence=os.path.relpath(mf, ROOT), manifest_sha256=sha256(mf),
                      per_object=po)


def stage_annotate(tag, lineage, res, seed, contamination=None, otag=None):
    """`tag` 定**输入**（集群表、r\\* 签字、子集清单 —— 换面板不动它们）；
    `otag` 定**产物文件名**，默认 = tag。传 `--out-suffix` 时产物落到别名文件，
    旧面板的结果原样留着，两套可并排比较。"""
    ot = otag or tag
    panel_ev = check_panel_caliber(PANEL_NAME, PANEL_PATH)
    rstar = check_rstar_signature(tag, res)
    seed, seed_ev = enforce_annotation_seed(tag, seed)
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
    cm.to_csv(f"{OUT}/{ot}_cluster_scores.csv", index=False)
    log(f"写出 {ot}_cluster_scores.csv（逐簇 × 全部 {len(types)} 型得分）")

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

    # wilcoxon top 基因（每簇，one-vs-rest）
    # ⚠️ 2026-09-21 修（首次实跑暴露）：原实现写 `scanpy.tl.rank_genes_groups(sub, "_mk",
    #   groups=[True], ...)` —— `_mk` 是 **bool 列**、不是 categorical，scanpy 取组标签时
    #   对它调 `.cat` ⇒ **每一簇都抛 AttributeError**，而 except 把异常**降级成字符串写进 CSV**
    #   （`<失败: AttributeError>`），日志只报"写出成功" ⇒ **静默产出空结果**（法则 0 禁止的
    #   那类静默降级）。现改为：一次对全体簇做 one-vs-rest（统计上也是"逐簇富集"的正解）；
    #   若因个别退化解失败，再逐簇单跑，且**每次都把失败原因写进 CSV 并在日志里喊**。
    # ⚠️ 2026-09-21 补：**判簇不能只看基因名**。权威做法要求每条 marker 带
    #   **调整后 p 值 + log2 倍数变化 + 表达比例**（Heumos 2023 best-practices；本项目
    #   2026-09-21 用户裁定补上）。原实现只写基因名 —— 那会把"真高表达"和"低表达但显著"
    #   混在一起（eg 簇 15 的 top 里同时有 EGFR 与 IGKC，光看名字分不出谁是本体谁是背景）。
    #   fraction expressing 现算（不依赖 scanpy 版本的 pts 字段是否还在）。
    _X = sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
    _gi = {str(g): i for i, g in enumerate(sub.var_names)}

    def _marker_stats(group):
        """取一个簇的 rank_genes_groups 统计量 + 现算簇内/簇外表达比例。"""
        rg = sub.uns["rank_genes_groups"]
        names = rg["names"][group]
        lfc, pv, padj, scr = (rg.get(k) for k in
                              ("logfoldchanges", "pvals", "pvals_adj", "scores"))
        mask = (cl == str(group))
        rows = []
        for r, gname in enumerate(names):
            j = _gi.get(str(gname))
            if j is None:                      # 理论上不会发生；发生了也不静默
                rows.append(dict(cluster=str(group), rank=r + 1, gene=str(gname),
                                 score=np.nan, log2fc=np.nan, pval=np.nan,
                                 pval_adj=np.nan, pct_in=np.nan, pct_out=np.nan))
                continue
            col = _X[:, j]
            ine = np.asarray(col[mask].todense()).ravel()
            oute = np.asarray(col[~mask].todense()).ravel()
            rows.append(dict(
                cluster=str(group), rank=r + 1, gene=str(gname),
                score=round(float(scr[group][r]), 4),
                log2fc=round(float(lfc[group][r]), 4) if lfc is not None else np.nan,
                pval=float(pv[group][r]) if pv is not None else np.nan,
                pval_adj=float(padj[group][r]) if padj is not None else np.nan,
                pct_in=round(float((ine > 0).mean()), 4),
                pct_out=round(float((oute > 0).mean()), 4)))
        return rows

    top_rows, marker_rows = [], []
    try:
        sc.tl.rank_genes_groups(sub, "clu", method="wilcoxon", n_genes=DRIVER_TOPN)
        names = sub.uns["rank_genes_groups"]["names"]
        for c in ucl:
            top_rows.append(dict(cluster=c, top_genes=";".join(map(str, names[c]))))
            marker_rows += _marker_stats(c)
    except Exception as e:
        log(f"🔴 一次对全体簇做 rank_genes_groups 失败（{type(e).__name__}: {e}）"
            f" —— 改为逐簇单跑，失败簇会在 CSV 里标出")
        top_rows, marker_rows = [], []
        for c in ucl:
            if int((cl == c).sum()) < 2 or int((cl != c).sum()) < 2:
                g = [f"<失败: 簇 {c} 细胞数不足，无法做检验>"]
            else:
                try:
                    sc.tl.rank_genes_groups(sub, "clu", groups=[c], reference="rest",
                                            method="wilcoxon", n_genes=DRIVER_TOPN)
                    g = list(sub.uns["rank_genes_groups"]["names"][c])
                    marker_rows += _marker_stats(c)
                except Exception as e2:
                    log(f"🔴 簇 {c} 的 rank_genes_groups 失败（{type(e2).__name__}: {e2}）")
                    g = [f"<失败: {type(e2).__name__}>"]
            top_rows.append(dict(cluster=c, top_genes=";".join(map(str, g))))
    df_top = pd.DataFrame(top_rows)
    n_bad = int(df_top["top_genes"].str.startswith("<失败").sum())
    if n_bad:
        log(f"🔴 {ot}_cluster_topgenes.csv 里有 {n_bad}/{len(df_top)} 簇**取不到**富集基因"
            f" —— 该文件不可用于判读，须先修")
    pd.DataFrame(top_rows).to_csv(f"{OUT}/{ot}_cluster_topgenes.csv", index=False)
    log(f"写出 {ot}_cluster_topgenes.csv（{len(df_top)} 簇，失败 {n_bad}）")

    dm = pd.DataFrame(marker_rows)
    dm.to_csv(f"{OUT}/{ot}_cluster_markers.csv", index=False)
    n_sig = int((dm["pval_adj"] < 0.05).sum())
    log(f"写出 {tag}_cluster_markers.csv（{len(dm)} 条 = {len(ucl)} 簇 × {DRIVER_TOPN} 基因；"
        f"带 padj/log2FC/簇内外表达比例；其中 padj<0.05 的 {n_sig} 条）")

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
    # ⚠️ 2026-09-21 补：原 flag 只有 tiny/rare/thin 三条，**漏了"单一供体"** ——
    #   而单供体簇**不得当作细胞类型上报**（图谱类工作的共识）。这里只标**事实**
    #   （n_patients==1），**不设阈值**：n_patients 与 top_patient_frac 两列都已在表里，
    #   多个供体但明显偏斜的簇（eg 77% 来自一人）由人工看数判。
    ann["donor_flag"] = ""
    ann.loc[ann["n_patients"] == 1, "donor_flag"] = "single_donor"
    ann["decision"] = ""        # 🔴 留空待人工填，脚本不自动裁决
    ann["reason"] = ""

    # ---- 跨谱系污染剔除（2026-09-21）：证据来自 06_contamination_check.py，不在此处新增判据 ----
    ann["excluded"] = False
    ann["excluded_reason"] = ""
    n_excl_cells = 0
    if contamination:
        if not os.path.exists(contamination):
            raise SystemExit(f"🔴 找不到污染体检表 {contamination} —— 须先跑 06_contamination_check.py")
        con = pd.read_csv(contamination)
        need = {"cluster", "flag", "n_cells"}
        if not need <= set(con.columns):
            raise SystemExit(f"🔴 {contamination} 缺列 {need - set(con.columns)}")
        cmap = con.set_index(con["cluster"].astype(str))["flag"].to_dict()
        # 只有这三种旗子是"须剔除"；其余一律保留。未知串直接报错，不让它被静默当成剔除。
        EXCL_FLAGS = {"面板", "CellTypist", "面板+CellTypist"}
        for i, c in enumerate(ann["cluster"]):
            f = cmap.get(str(c))
            if f is None:
                raise SystemExit(f"🔴 污染体检表里没有簇 {c}（口径不一致，拒绝猜）")
            if f not in EXCL_FLAGS:
                if f != "保留":
                    raise SystemExit(f"🔴 污染体检表里出现未知旗子 {f!r}（簇 {c}）"
                                     f"—— 拒绝猜它是否该剔。已知：{sorted(EXCL_FLAGS)} 或 '保留'")
                continue
            ann.at[i, "excluded"] = True
            ann.at[i, "excluded_reason"] = (
                f"跨谱系污染（{f} 命中）；判据见 {os.path.relpath(contamination, ROOT)}")
        n_excl_cells = int(ann.loc[ann["excluded"], "n_cells"].sum())

        # 逐细胞剔除清单（用户 2026-09-21 裁定：落一份清单，逐条可查）
        bad_c = set(ann.loc[ann["excluded"], "cluster"].astype(str))
        m = np.isin(cl, list(bad_c))
        ex = pd.DataFrame({"cell_barcode": sub.obs_names[m], "cluster": cl[m]})
        for k in ("patient_id", "sample_id", "stage"):
            if k in sub.obs:
                ex[k] = sub.obs[k].to_numpy()[m]
        ex["evidence"] = ex["cluster"].map(cmap)
        ex.sort_values(["cluster", "cell_barcode"], inplace=True)
        ex.to_csv(f"{OUT}/{ot}_excluded_cells.csv.gz", index=False)
        # 剩余细胞的纯 barcode 清单（供"干净子集重聚类"稳定性校验 run 用）
        keep_bcs = sub.obs_names[~m]
        with open(f"{OUT}/{ot}_nocontam_subset_barcodes.txt", "w") as fh:
            fh.write("\n".join(keep_bcs) + "\n")
        log(f"🔴 剔除 {len(bad_c)} 个簇 / {n_excl_cells:,} 核"
            f"（占子集 {n_excl_cells / len(bcs) * 100:.2f}%）；逐细胞清单 "
            f"{os.path.relpath(f'{OUT}/{ot}_excluded_cells.csv.gz', ROOT)}")
        log(f"   剩余 {len(keep_bcs):,} 核 → {os.path.relpath(f'{OUT}/{ot}_nocontam_subset_barcodes.txt', ROOT)}"
            f"（**仅供稳定性校验**，不是签字输入）")

    ann.to_csv(f"{OUT}/{ot}_cluster_annotation.csv", index=False)
    log(f"写出 {ot}_cluster_annotation.csv（decision 列留空待人工填）")

    keep = ann[~ann["excluded"]]
    print(f"\n=== {lineage} 簇级 argmax（{col}）===")
    print(ann.to_string(index=False))
    if n_excl_cells:
        print(f"\n⚠️ 上式含 {int(ann['excluded'].sum())} 个**已剔除非上皮簇**（{n_excl_cells:,} 核）；"
              f"下列统计**只在保留簇上算**")
    vc = keep.groupby("argmax")["n_cells"].agg(["count", "sum"])
    print(f"\n胜出分布（保留簇 {len(keep)} 个 / {int(keep['n_cells'].sum()):,} 核）:\n{vc.to_string()}")
    never = [t for t in types if t not in set(keep["argmax"])]
    print(f"\n从未胜出的亚型: {never}")
    if unscorable:
        print(f"无法打分（可用基因<2）: {unscorable}")
    thin_src = [k for k, v in n_slot.items() if v <= 2]
    print(f"🔴 源表面板本来就只有 1–2 槽的亚型: {thin_src}")

    mf = dict(tag=tag, out_tag=ot, lineage=lineage, res=res, seed=seed, col=col,
              n_cells=int(len(bcs)), n_clusters=len(ucl),
              seed_selection=seed_ev,
              rstar_signed_by=rstar["signed_by"],
              rstar_relaxed=bool(rstar.get("relaxed")),
              rstar_json_sha256=sha256(f"{OUT}/{tag}_rstar.json"),
              panel=PANEL_NAME,
              panel_sha256=sha256(PANEL_PATH),
              panel_caliber=(None if panel_ev is None else dict(
                  json=f"{OUT}/panel_caliber_{PANEL_NAME}.json",
                  json_sha256=sha256(f"{OUT}/panel_caliber_{PANEL_NAME}.json"),
                  signed_by=panel_ev.get("signed_by"),
                  rule=panel_ev.get("rule"))),
              subsets_manifest=sha256(SUBSETS),
              # R5（产物可复现且哈希）：核心输入一律留哈希，否则"换了输入但 manifest 看不出来"。
              # 2026-09-21 补：原先只记了面板/子集清单/污染表，**漏了簇表、表达矩阵与脚本自身**
              # —— 而这三者恰恰决定了这份注释是什么。上皮那份正因如此留下了断链
              # （记的是改口径前的污染表哈希，与现存文件对不上）。
              inputs=dict(
                  h5ad=dict(path=os.path.relpath(H5, ROOT), sha256=sha256(H5)),
                  clusters=dict(path=os.path.relpath(f"{TRAD}/{tag}/clusters.csv.gz", ROOT),
                                sha256=sha256(f"{TRAD}/{tag}/clusters.csv.gz"),
                                cluster_col=col),
                  lineage_subsets=dict(path=os.path.relpath(bc_file, ROOT),
                                       sha256=sha256(bc_file), n_cells=len(bcs)),
                  script=dict(path=os.path.relpath(os.path.abspath(__file__), ROOT),
                              sha256=sha256(os.path.abspath(__file__))),
              ),
              score_genes=dict(ctrl_size=CTRL_SIZE, random_state=SCORE_SEED),
              missing_genes=missing, thin_panel_types=thin,
              unscorable_types=unscorable,
              unannotatable=list(CP.UNANNOTATABLE),
              source_table_defects=list(CP.SOURCE_TABLE_DEFECTS),
              summary_over=("保留簇（已剔除跨谱系污染簇）" if n_excl_cells else "全部簇"),
              n_clusters_kept=int(len(keep)),
              never_won=never, win_distribution=vc["sum"].to_dict(),
              contamination=(
                  dict(csv=os.path.relpath(contamination, ROOT),
                       csv_sha256=sha256(contamination),
                       excluded_clusters=sorted(ann.loc[ann["excluded"], "cluster"].astype(str)),
                       n_cells_excluded=n_excl_cells,
                       cells_list=f"{ot}_excluded_cells.csv.gz",
                       clean_barcodes=f"{ot}_nocontam_subset_barcodes.txt",
                       note="剔除判据完全来自 06_contamination_check.py（规则在该脚本写定时冻结）；"
                            "本脚本不新增判据。干净 barcode 清单仅供稳定性校验，不是签字子集的替代品。")
                  if n_excl_cells or contamination else None))
    with open(f"{OUT}/{ot}_annotation_manifest.json", "w") as fh:
        json.dump(mf, fh, ensure_ascii=False, indent=2)
    log(f"写出 {ot}_annotation_manifest.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["select", "annotate"], required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--lineage", default=None)
    ap.add_argument("--res", type=float, default=None)
    ap.add_argument("--seed", type=int, default=0,
                    help="注释用的重聚类种子。⚠️ 不是自由参数：由 "
                         "05_seed_representativeness.py 的预注册判据决定，"
                         "传错即 SystemExit（见 enforce_annotation_seed）。")
    ap.add_argument("--relax-seed-min", type=float, default=None,
                    help="🔴 事后放宽跨种子稳定性阈值。仅限用户明确裁定的谱系使用；"
                         "产物会登记 relaxed=true，结果只能作探索性结论。")
    ap.add_argument("--contamination-csv", default=None,
                    help="跨谱系污染体检表（06_contamination_check.py 产出）。给了就按其中的 "
                         "flag 剔除被标簇，并落逐细胞剔除清单 + 干净 barcode 清单；"
                         "不给则不做任何剔除（保留旧行为）。")
    ap.add_argument("--panel", choices=sorted(PANEL_MODULES), default="s1",
                    help="亚型（L2）面板来源。默认 s1 = classic_panels（Table S1，行为不变）；"
                         "s4 = classic_panels_s4（Table S4，仅供内皮/髓系，须先签 panel_caliber_s4.json）。")
    ap.add_argument("--out-suffix", default="",
                    help="产物文件名后缀，只改**输出**、不改输入（集群表/r*/子集清单仍按 --tag 取）。"
                         "用于换面板重跑时保留旧结果，例如 --out-suffix _s4 ⇒ endoA_s4_*。")
    a = ap.parse_args()

    global CP, PANEL_PATH, PANEL_NAME
    CP, PANEL_PATH = load_panel_module(a.panel)
    PANEL_NAME = a.panel

    if a.stage == "select":
        stage_select(a.tag, relax_seed_min=a.relax_seed_min)
    else:
        if a.lineage is None or a.res is None:
            raise SystemExit("annotate 须给 --lineage 与 --res")
        stage_annotate(a.tag, a.lineage, a.res, a.seed,
                       contamination=a.contamination_csv,
                       otag=a.tag + a.out_suffix)


if __name__ == "__main__":
    main()
