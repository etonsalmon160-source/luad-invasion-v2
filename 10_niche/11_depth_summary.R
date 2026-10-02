#!/usr/bin/env Rscript
# 11_depth_summary.R —— M6 §5.1 深度守卫的汇总与判读
#
# 输入：results/10_niche/depth_guard/<agf>/<slide>__<cond>.tsv   （08_depth_guard.R 产出）
#       results/10_niche/depth_guard/<agf>/<slide>__manifest.json
#
# 输出（results/10_niche/depth_guard/）：
#   summary_per_slide.tsv  逐切片 × 条件 × 全局类型 的 spot 占比 + 归并相关余量
#   summary_cohort.tsv     按患者聚簇后跨患者的中位占比（§5.3：同一患者多张切片不是独立样本）
#   monotonicity.tsv       逐类型：五档单调性 + 共同深度(down) vs 全片(full) 的配对差
#   manifest.json          含 08 复现主臂的 ARI 分布（自检）
#
# 🔴 判读口径 = §5.1.1 原文：**单调 = 稳定；翻转 = 降级待定**。这里只做机械判定、只上报。
# 🔴 **全部条件（含 full 基线）都由 08 用同一套机器产生** ⇒ 档间差异不混入"换机制"的水分。
# 🔴 本臂不给任何恶性标签（§7 禁令）。跨切片汇总一律按 patient 聚类（§5.3）。

suppressMessages({ library(data.table) })
ROOT  <- "/home/eto/luad_v2"
NICHE <- file.path(ROOT, "results/10_niche")
DGD   <- file.path(NICHE, "depth_guard")
QS    <- sprintf("q%d", 1:5)
CONDS <- c(QS, "down", "full")

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
patient_of <- function(s) sub("^[^_]+_([^_]+)_.*$", "\\1", s)

all_rows <- list(); repro <- list()
for (agf in c("agfT", "agfF")) {
  dd <- file.path(DGD, agf)
  if (!dir.exists(dd)) next
  fs <- list.files(dd, pattern = "__.*\\.tsv$", full.names = TRUE)
  step("%s：%d 个条件文件", agf, length(fs))
  for (f in fs) {
    x <- fread(f)
    if (!nrow(x)) next
    ## 同一全局类型可能被多个域映到 ⇒ 占比相加（§5.1.1：按域归入后再统计占比）
    y <- x[, .(n_spot_gt = sum(n_spot_dom), frac = sum(frac_dom),
               n_spot = n_spot[1], n_domain = n_domain[1],
               depth_med = depth_med[1], nn_med_um = nn_med_um[1],
               margin_min = min(cor_margin), margin_med = median(cor_margin)),
           by = .(slide, agf, cond, global_type = best_type)]
    all_rows[[length(all_rows) + 1]] <- y
  }
  mf <- list.files(dd, pattern = "__manifest\\.json$", full.names = TRUE)
  for (f in mf) {
    j <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
    if (!is.null(j)) repro[[length(repro) + 1]] <- data.table(
      agf = agf, slide = j$slide, ari = if (is.null(j$repro_vs_main_ari)) NA_real_ else j$repro_vs_main_ari)
  }
}
if (!length(all_rows)) stop("没有任何深度守卫结果 ⇒ 先跑 08_depth_guard.R", call. = FALSE)
PS <- rbindlist(all_rows, fill = TRUE)
PS[, patient := patient_of(slide)]
setorder(PS, agf, slide, cond, global_type)
fwrite(PS, file.path(DGD, "summary_per_slide.tsv"), sep = "\t")

## ——— 跨患者汇总（§5.3：先患者内取中位，再跨患者取中位）———
pat <- PS[, .(frac_pat = median(frac)), by = .(agf, patient, cond, global_type)]
COH <- pat[, .(n_patient = .N, frac_median = median(frac_pat),
               frac_q25 = quantile(frac_pat, .25), frac_q75 = quantile(frac_pat, .75)),
           by = .(agf, cond, global_type)]
setorder(COH, agf, global_type, cond)
fwrite(COH, file.path(DGD, "summary_cohort.tsv"), sep = "\t")

## ——— 单调性判读：逐切片算 frac 对档序的 Spearman rho；再看跨患者中位曲线是否严格单调 ———
## 🔴 data.table 作用域陷阱：`DT[agf == agf]` 里两边都会被当成**列名** ⇒ 恒为 TRUE（全表）。
##    循环变量一律起成表里没有的名字（AG / GT），否则每一组算的都是同一份全表。
mono <- list()
for (AG in unique(PS$agf)) for (GT in sort(unique(PS$global_type))) {
  sub <- PS[agf == AG & global_type == GT & cond %in% QS]
  if (!nrow(sub)) next
  sub[, qi := as.integer(sub("^q", "", cond))]
  rho <- sub[, .(rho = if (.N >= 3 && length(unique(frac)) > 1)
                       cor(qi, frac, method = "spearman") else NA_real_), by = slide]
  rho <- rho[is.finite(rho)]
  curve <- COH[agf == AG & global_type == GT & cond %in% QS][order(match(cond, QS))]$frac_median
  steps <- diff(curve)
  verdict <- if (length(curve) < 5) "样本不足"
             else if (all(steps > 0)) "单调上升" else if (all(steps < 0)) "单调下降" else "非单调（翻转）"
  ## 共同深度：down vs full，**同一患者内配对**差值的中位（§5.3）
  dp <- PS[agf == AG & global_type == GT & cond %in% c("down", "full")][
           , .(f = median(frac)), by = .(patient, cond)]
  d <- dcast(dp, patient ~ cond, value.var = "f")
  delta <- if (all(c("down", "full") %in% names(d))) median(d$down - d$full, na.rm = TRUE) else NA_real_
  mono[[length(mono) + 1]] <- data.table(
    agf = AG, global_type = GT, n_slide_rho = nrow(rho),
    frac_slide_rho_gt0 = if (nrow(rho)) mean(rho$rho > 0) else NA_real_,
    rho_median = if (nrow(rho)) median(rho$rho) else NA_real_,
    curve_q1 = curve[1], curve_q5 = curve[5], quintile_verdict = verdict,
    frac_full = COH[agf == AG & global_type == GT & cond == "full"]$frac_median,
    frac_down = COH[agf == AG & global_type == GT & cond == "down"]$frac_median,
    frac_down_minus_full = delta)
}
MO <- rbindlist(mono, fill = TRUE)
setorder(MO, agf, global_type)
fwrite(MO, file.path(DGD, "monotonicity.tsv"), sep = "\t")

step("落盘 → %s", DGD)
for (AG in c("agfT", "agfF")) {
  m <- MO[agf == AG]; if (!nrow(m)) next
  step("  %s：类型 %d 个；五档判读 单调上升 %d / 单调下降 %d / 非单调 %d",
       AG, nrow(m), sum(m$quintile_verdict == "单调上升"),
       sum(m$quintile_verdict == "单调下降"), sum(m$quintile_verdict == "非单调（翻转）"))
  d <- m[is.finite(frac_down_minus_full)]
  if (nrow(d)) step("    共同深度 vs 全片：|Δ| 中位 %.4f，|Δ|>0.02 的类型 %d/%d",
                    median(abs(d$frac_down_minus_full)), sum(abs(d$frac_down_minus_full) > 0.02), nrow(d))
}
RP <- if (length(repro)) rbindlist(repro) else data.table(agf = character(), slide = character(), ari = numeric())
if (nrow(RP)) for (AG in unique(RP$agf)) {
  v <- RP[agf == AG]$ari; v <- v[is.finite(v)]
  if (length(v)) step("  %s：08 复现主臂 ARI 中位 %.4f（最小 %.4f，n=%d）",
                      AG, median(v), min(v), length(v))
}
writeLines(jsonlite::toJSON(list(
  script = "10_niche/11_depth_summary.R",
  prereg = "NICHE_PREREG.md §5.1/§5.1.1（已签）；判读：单调=稳定，翻转=降级待定",
  n_slide = length(unique(PS$slide)), n_slide_by_cond = as.list(table(PS$cond)),
  conds = CONDS, n_type_agfT = nrow(MO[agf == "agfT"]), n_type_agfF = nrow(MO[agf == "agfF"]),
  repro_vs_main_ari = if (nrow(RP)) list(agfT = median(RP[agf == "agfT"]$ari, na.rm = TRUE),
                                          agfF = median(RP[agf == "agfF"]$ari, na.rm = TRUE)) else NULL,
  harvest_rule = "域按表达画像与全局共识类型画像相关取最大归入；同类型多域占比相加；另报 cor_margin",
  cluster_rule = "跨切片汇总先按 patient 取中位，再跨患者取中位（§5.3）",
  caveats = c("本臂不给任何恶性标签（§7 禁令）",
              "本臂是探索/敏感性臂：选中档本身未过闸 1（见 §4.1）",
              "down 档目标深度取全队列最浅切片的中位 ⇒ 深切片被稀释 5–22 倍，域塌缩不等于「没有生态位」",
              "归并用 argmax 相关，共识类型 1 独大（308/1005 域）⇒ 必须同看 cor_margin"),
  finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE, pretty = TRUE),
  file.path(DGD, "manifest.json"))
step("==== 11 结束 ====")
