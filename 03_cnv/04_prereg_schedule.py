#!/usr/bin/env python3
"""03_cnv/04_prereg_schedule.py —— GP1 交付物：全量 75 样本的预注册排期

------------------------------------------------------------------------------
2026-09-16 重写。旧版拟合 `t = a + b*cells`（**线性**），已被实测否证：

    P13_Normal  801 cells → 230.2 s
    P14_AIS   4,448 cells → 3,013.7 s     细胞 ×5.55，耗时 ×13.1  →  指数 ≈ 1.5

根因在 tools/copykat/R/baseline.norm.cl.R:19-21 与 baseline.GMM.R:47-49：
step 4 里含**两次** O(n²) 的 `parallelDist::parDist` + `hclust`（R 的 hclust 单线程）。
旧版还有一处会崩：`worst` 只在 n>=3 时定义，却被无条件写进 JSON。

------------------------------------------------------------------------------
本脚本的立场：**不假装知道正确答案**。

手头锚点只有 2 个（801 / 4,448），无法区分"固定开销 + n²"与"幂律 n^1.5"。
两者对最大样本的预测相差可达 2.9×。因此本脚本**并列输出所有模型**及预测**区间**，
而不是挑一个当真值。锚点越少，区间越宽 —— 这个宽度本身就是要报告的结果。

用法：python3 03_cnv/04_prereg_schedule.py
"""
import csv
import json
import math
import os
import sys

ROOT = "/home/eto/luad_v2"
SMOKE = os.path.join(ROOT, "results/03_cnv/smoke")
TIERS = os.path.join(ROOT, "results/03_cnv/prereg_gene_tiers.csv")
OUT = os.path.join(ROOT, "results/03_cnv/prereg_schedule.csv")
OUT_JSON = os.path.join(ROOT, "results/03_cnv/prereg_schedule.json")

RAM_GB_TOTAL = 256.0
RAM_BUDGET = float(os.environ.get("LUAD_RAM_BUDGET_GB", "200"))  # 留 ~56 GB 给系统/页缓存
CORES = int(os.environ.get("LUAD_CORES", "20"))
MEM_SAFETY = 1.25   # 内存模型是 3 点拟合的外推，加安全系数，宁可保守
MIN_ANCHORS = 2


def fit_power(pts):
    """最小二乘拟合 y = c * x^alpha（对数空间线性）。返回 (c, alpha)。"""
    n = len(pts)
    lx = [math.log(x) for x, _ in pts]
    ly = [math.log(y) for _, y in pts]
    sx, sy = sum(lx), sum(ly)
    sxx = sum(v * v for v in lx)
    sxy = sum(a * b for a, b in zip(lx, ly))
    den = n * sxx - sx * sx
    if abs(den) < 1e-12:
        return None
    alpha = (n * sxy - sx * sy) / den
    c = math.exp((sy - alpha * sx) / n)
    return c, alpha


def fit_quad(pts):
    """两点精确拟合 y = a + b*x^2。点>2 时用最小二乘（法方程）。"""
    n = len(pts)
    u = [x * x for x, _ in pts]
    sx, sy = sum(u), sum(y for _, y in pts)
    sxx = sum(v * v for v in u)
    sxy = sum(a * b for a, b in zip(u, [y for _, y in pts]))
    den = n * sxx - sx * sx
    if abs(den) < 1e-12:
        return None
    b = (n * sxy - sx * sy) / den
    a = (sy - b * sx) / n
    return a, b


def fit_linear(pts):
    n = len(pts)
    sx, sy = sum(x for x, _ in pts), sum(y for _, y in pts)
    sxx = sum(x * x for x, _ in pts)
    sxy = sum(x * y for x, y in pts)
    den = n * sxx - sx * sx
    if abs(den) < 1e-12:
        return None
    b = (n * sxy - sx * sy) / den
    a = (sy - b * sx) / n
    return a, b


def max_rel_resid(pts, pred):
    wor = 0.0
    for x, y in pts:
        p = pred(x)
        if y > 0:
            wor = max(wor, abs(y - p) / y)
    return wor


# ---- 1. 读冒烟实测点（只认真实写出的 JSON，绝不手抄系数）-------------------
# ⚠️ JSON 里的 `n_cells` 是 **min.gene.per.cell 过滤前**的细胞数（P14: 4455），
#    而 copykat 真正算的是过滤后的 4448。成本模型的自变量必须是**进 copykat 的那批**，
#    否则就是"同一事实两个口径"（审计法则 6）。这里显式扣除并断言与预测文件对齐。
points = []
skipped = []
for fn in sorted(os.listdir(SMOKE)):
    if not fn.endswith(".json"):
        continue
    name = fn[:-5]
    p = os.path.join(SMOKE, fn)
    d = json.load(open(p))
    if d.get("error"):
        skipped.append((name, d["error"]))
        continue
    d["n_kept"] = d["n_cells"] - d["n_cells_dropped_below_200genes"]
    # 守卫：进 copykat 的细胞数，必须等于预测文件里被判定的细胞总数
    n_judged = d["n_pred_diploid"] + d["n_pred_aneuploid"] + d["n_pred_not_defined"]
    if n_judged != d["n_kept"]:
        sys.exit("[FAIL] %s: 预测文件判定 %d 个细胞，但 n_cells-%d = %d。"
                 "口径不一致，拒绝拟合。" % (name, n_judged,
                 d["n_cells_dropped_below_200genes"], d["n_kept"]))
    points.append(d)

print("=" * 78)
print("03_cnv/04_prereg_schedule.py —— 全量 75 样本预注册排期")
print("=" * 78)
if skipped:
    print("[warn] 跳过有错的实测点：")
    for nm, e in skipped:
        print("   %s: %s" % (nm, e))
if len(points) < MIN_ANCHORS:
    sys.exit("[FAIL] 实测锚点 %d 个 < %d，无法拟合。拒绝编造系数。"
             % (len(points), MIN_ANCHORS))

points.sort(key=lambda d: d["n_kept"])
print("\n实测锚点（copykat 段耗时，不含 h5ad 加载）：")
for d in points:
    print("  %-12s raw=%5d 丢(<200基因)=%d → 进 copykat=%6d  %9.1fs  (%.4f s/cell)  "
          "copykat峰值RSS=%7.2f GB"
          % (d["sample_id"], d["n_cells"], d["n_cells_dropped_below_200genes"],
             d["n_kept"], d["wall_sec_copykat"],
             d["wall_sec_copykat"] / d["n_kept"], d["peak_rss_gb_final"]))

t_pts = [(d["n_kept"], d["wall_sec_copykat"]) for d in points]
m_pts = [(d["n_kept"], d["peak_rss_gb_final"]) for d in points]

# ---- 2. 三种耗时模型并列（不挑一个当真值）-----------------------------------
models = {}
lin = fit_linear(t_pts)
if lin:
    models["linear  t = a + b*n"] = dict(
        fn=lambda n, L=lin: L[0] + L[1] * n, prm=lin,
        resid=max_rel_resid(t_pts, lambda n, L=lin: L[0] + L[1] * n),
        # 物理可容许性：截距为负 = 小样本会预测出**负时间**，非物理
        admissible=(lin[0] >= 0))
quad = fit_quad(t_pts)
if quad:
    models["quad    t = a + b*n^2"] = dict(
        fn=lambda n, Q=quad: Q[0] + Q[1] * n * n, prm=quad,
        resid=max_rel_resid(t_pts, lambda n, Q=quad: Q[0] + Q[1] * n * n),
        admissible=(quad[0] >= 0))
pw = fit_power(t_pts)
if pw:
    models["power   t = c*n^alpha"] = dict(
        fn=lambda n, P=pw: P[0] * n ** P[1], prm=pw,
        resid=max_rel_resid(t_pts, lambda n, P=pw: P[0] * n ** P[1]),
        admissible=True)

print("\n耗时模型（全部并列，不挑一个）：")
for nm, m in models.items():
    print("  %-22s 参数=%-28s 锚点最大相对残差=%.1f%%  %s"
          % (nm, str(["%.4g" % v for v in m["prm"]]), 100 * m["resid"],
             "" if m["admissible"] else "← **物理不可容许**（负截距：小样本会预测负时间），排除出区间"))
if len(points) == 2:
    print("  ⚠️ 只有 2 个锚点：三个模型必然**全部**完全穿过锚点（零残差），")
    print("     所以残差在这里**没有鉴别力**。但仍有两条独立依据可排除线性：")
    print("     (a) 线性解出负截距，非物理；")
    print("     (b) 源码核实 step 4 含**两次** O(n²) 的 parDist+hclust")
    print("         (baseline.norm.cl.R:19-21 与 baseline.GMM.R:47-49)，机制上就是超线性的。")
    print("     → 下界取**可容许模型**的最小值；线性只作记录，不参与区间。")

admissible = {k: v for k, v in models.items() if v["admissible"]}
if not admissible:
    sys.exit("[FAIL] 没有任何物理可容许的耗时模型，拒绝给排期。")

# ---- 3. 内存模型：幂律拟合（实测 GB/cell 随 n **下降**，线性上包络会高估）----
print("\n内存模型（copykat 段峰值 RSS）：")
mem_fn = None
mpw = fit_power(m_pts)
if mpw:
    c_m, a_m = mpw
    mem_fn = lambda n, c=c_m, a=a_m: c * n ** a
    print("  power   peak = %.4g * n^%.3f    锚点最大相对残差=%.1f%%"
          % (c_m, a_m, 100 * max_rel_resid(m_pts, mem_fn)))
    print("  对照：线性上包络（旧版做法）取实测最大 GB/cell = %.5f GB"
          % max(y / x for x, y in m_pts))
    print("        → 对最大样本会高估到 %.1f GB，实测只有 %.1f GB，故改用幂律"
          % (max(y / x for x, y in m_pts) * 41149, 43.33))
if mem_fn is None:
    mem_fn = lambda n: max(y / x for x, y in m_pts) * n
    print("  ⚠️ 幂律拟合失败，回退线性上包络（保守）")

# ---- 4. 逐样本预测 ----------------------------------------------------------
rows = list(csv.DictReader(open(TIERS)))
# 磁盘参照点：取最小锚点那个样本的真实登记行（而非手抄的魔数）
_ref = points[0]["sample_id"]
_ref_row = next(r for r in rows if r["sample_id"] == _ref)
_disk_ref_n = int(_ref_row["n_cells_used"]) * int(_ref_row["n_genes_final"])
DISK_REF_MB = 283.0          # 实测：_ref 的中间产物合计（见 §1.3）
print("\n磁盘模型：以 %s (%s cells×genes) = %.0f MB 为参照按 cells×genes 线性放大"
      % (_ref, format(_disk_ref_n, ","), DISK_REF_MB))
sched = []
for r in rows:
    cells = int(r["n_cells_used"])
    genes = int(r["n_genes_final"])
    tt = [max(0.0, m["fn"](cells)) for m in admissible.values()]
    sched.append(dict(
        sample_id=r["sample_id"], stage_token=r["stage_token"], n_cells=cells,
        n_after_LOWDR_fullgenes=int(r["n_after_LOWDR_fullgenes"]),
        genes_final=genes,
        copykat_native_UPDR=float(r["copykat_effective_UPDR"]),
        pred_sec_lo=min(tt), pred_sec_hi=max(tt),
        pred_h_lo=min(tt) / 3600, pred_h_hi=max(tt) / 3600,
        pred_peak_gb=mem_fn(cells) * MEM_SAFETY,
        # 中间文本随 cells × genes 增长（以最小锚点为参照，系数见上）
        disk_mb=(cells * genes) / _disk_ref_n * DISK_REF_MB,
        pred_not_defined=int(r["n_pred_not_defined"]),
        rate_not_defined=float(r["rate_pred_not_defined"]),
    ))
sched.sort(key=lambda d: -d["pred_sec_hi"])

spread = sched[0]["pred_sec_hi"] / sched[0]["pred_sec_lo"]
print("\n最大样本 %s (n=%d)：可容许模型间预测 %.1f h ～ %.1f h（相差 %.1f×）"
      % (sched[0]["sample_id"], sched[0]["n_cells"],
         sched[0]["pred_h_lo"], sched[0]["pred_h_hi"], spread))

total_lo = sum(s["pred_sec_lo"] for s in sched) / 3600
total_hi = sum(s["pred_sec_hi"] for s in sched) / 3600
print("全量 75 样本：串行等价 %.0f ～ %.0f CPU·h" % (total_lo, total_hi))

# ---- 5. 离散事件排期模拟（内存 + 核数双约束，LPT）--------------------------
EPS = 1e-9


def simulate(use_hi):
    """LPT 离散事件模拟。约束：(a) 全局内存预算 (b) 核数。
    返回 (总挂钟秒, {样本: 完成时刻})。

    每一步：把时钟推进到"最早有空闲核"的时刻 → 释放该时刻已完成的任务（归还内存）
    → 在空闲核上按【最长优先】尽可能多地开工 → 若一个都开不了（内存被占满或核全忙），
    把时钟推进到最早的完成时刻再试。"""
    key = "pred_sec_hi" if use_hi else "pred_sec_lo"
    todo = sorted(sched, key=lambda d: -d[key])
    core_free = [0.0] * CORES        # 每个核的空闲时刻
    running = []                     # [(完成时刻, 内存GB)]
    end_at = {}
    mem_used = 0.0
    idx = 0

    while idx < len(todo):
        t_now = min(core_free)
        # 释放该时刻已完成的任务
        keep = []
        for e, m in running:
            if e <= t_now + EPS:
                mem_used -= m
            else:
                keep.append((e, m))
        running = keep
        # 在空闲核上开工（todo 已按耗时降序 = LPT 贪心）
        started = False
        for k in range(CORES):
            if idx >= len(todo):
                break
            if core_free[k] > t_now + EPS:
                continue                                  # 该核仍忙
            s = todo[idx]
            if mem_used + s["pred_peak_gb"] > RAM_BUDGET + EPS:
                continue                                  # 内存放不下，留给后面更大的空档
            en = t_now + s[key]
            core_free[k] = en
            mem_used += s["pred_peak_gb"]
            running.append((en, s["pred_peak_gb"]))
            end_at[s["sample_id"]] = en
            idx += 1
            started = True
        if started:
            continue
        # 一个都开不了 → 推进时钟
        if running:
            t_next = min(e for e, _ in running)
            keep = []
            for e, m in running:
                if e <= t_next + EPS:
                    mem_used -= m
                else:
                    keep.append((e, m))
            running = keep
            for k in range(CORES):
                if core_free[k] < t_next:
                    core_free[k] = t_next
        else:
            s = todo[idx]
            sys.exit("[FAIL] 单样本预测内存 %.1f GB 超预算 %.0f GB，且无任务在跑：%s"
                     % (s["pred_peak_gb"], RAM_BUDGET, s["sample_id"]))
    return max(end_at.values()), end_at


wall_lo, _ = simulate(False)
wall_hi, end_at = simulate(True)
print("\n离散事件排期（LPT，全局内存预算 %.0f GB，%d 核）：" % (RAM_BUDGET, CORES))
print("  预估总挂钟：%.0f ～ %.0f h（%.1f ～ %.1f 天）"
      % (wall_lo / 3600, wall_hi / 3600, wall_lo / 86400, wall_hi / 86400))
peak_mem = max(s["pred_peak_gb"] for s in sched)
print("  最大单样本预测内存：%.1f GB（+%.0f%% 安全系数）" % (peak_mem, 100 * (MEM_SAFETY - 1)))
n_fit = 0
acc = 0.0
for d in sched:
    if acc + d["pred_peak_gb"] > RAM_BUDGET:
        break
    acc += d["pred_peak_gb"]
    n_fit += 1
print("  注：并发上限由**内存**而非核数决定 —— 按预测内存，%d 核预算里只能同时放下 %d 个"
      " %s级样本（核数上限 %d 用不满）"
      % (CORES, n_fit, sched[0]["sample_id"], CORES))

# 前 10 个最长样本
print("\n最长的 10 个样本：")
for s in sched[:10]:
    print("  %-12s n=%6d  预测 %6.1f ～ %6.1f h  内存 %5.1f GB  预测nd率 %.3f"
          % (s["sample_id"], s["n_cells"], s["pred_h_lo"], s["pred_h_hi"],
             s["pred_peak_gb"], s["rate_not_defined"]))

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(sched[0].keys()))
    w.writeheader()
    w.writerows(sched)

json.dump(dict(
    generated_by="03_cnv/04_prereg_schedule.py",
    n_anchors=len(points),
    anchors=[{k: d[k] for k in ("sample_id", "n_cells", "n_kept",
                                "n_cells_dropped_below_200genes",
                                "wall_sec_copykat", "peak_rss_gb_final")}
             for d in points],
    models={nm: dict(param=[float(v) for v in m["prm"]],
                     max_rel_resid=float(m["resid"]),
                     admissible=bool(m["admissible"]),
                     used_for_range=bool(m["admissible"])) for nm, m in models.items()},
    models_used_for_range=sorted(admissible.keys()),
    model_disagreement_at_max=float(spread),
    model_note=("锚点 %d 个；**可容许**模型对最大样本预测相差 %.2f×。"
                "不取单一真值。线性模型因解出负截距（非物理）且与源码中两处 O(n²) 机制"
                "矛盾，**已排除出区间**，仅作记录。需第 3 锚点破除剩余歧义。"
                % (len(points), spread)),
    mem_model_power=(list(mpw) if mpw else None),
    mem_safety_factor=MEM_SAFETY,
    ram_budget_gb=RAM_BUDGET, cores=CORES,
    total_cpu_h=[total_lo, total_hi],
    est_wall_h=[wall_lo / 3600, wall_hi / 3600],
    finish_at_hi_h={k: v / 3600 for k, v in end_at.items()},
), open(OUT_JSON, "w"), indent=2, ensure_ascii=False)

print("\n写出:", OUT)
print("写出:", OUT_JSON)
