#!/usr/bin/env Rscript
## 20_report_cohort.R —— 借锚臂**全队列上报**（只读产物；不筛数据、不下结论）
##
## 依据（两处都是**已签字**的义务，不是本脚本新增的口径）：
##   · `SPATIAL_CNV_PREREG.md` §19.7 SC4：每患者报 上皮 spot 数 / 锚 spot 数 / cf 分位数 /
##     **四档 0.02 0.05 0.10 0.15** / 各染色体臂均值 / 🔴 **参考来源（自己 or 借自哪个池）**；
##     **不设通过率阈值、不按阈值筛 spot、不许事后挑患者**（法则 3.2）。
##   · `RCTD_PREREG.md` §13.7.3：下游**每条主要结论**须附 **margin 梯度五档复核**
##     （`margin ≥ 0.00 / 0.05 / 0.10 / 0.15 / 0.20`，**全档照报**，不挑档、不筛数据）。
##
## 跑法：
##   export R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library
##   Rscript --vanilla 08_spatial_deconv/20_report_cohort.R
##
## 读：results/08_spatial_deconv/spatial_cnv/cohort_borrow/<P>/{manifest.json,*_spot_scores.rds,SC3_holdout.rds}
##     results/08_spatial_deconv/rctd_a/per_slide/<slide>.weights.tsv.gz（只为算 margin）
## 写：results/08_spatial_deconv/spatial_cnv/cohort_borrow/report/
##
## 🔴 **本脚本只上报。** 不得把任何一格读成恶性判定；不得因某格好看就写结论。

suppressPackageStartupMessages({ library(jsonlite) })

RES  <- "/home/eto/luad_v2/results/08_spatial_deconv"
BASE <- file.path(RES, "spatial_cnv/cohort_borrow")
WTS  <- file.path(RES, "rctd_a/per_slide")
OUT  <- file.path(BASE, "report")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TIERS  <- c(0.02, 0.05, 0.10, 0.15)                     # §19.7 SC4 四档（**只报**）
MGRID  <- c(0.00, 0.05, 0.10, 0.15, 0.20)               # §13.7.3 margin 五档（**只报**）
TH_MAIN <- 0.05                                         # 与 19_*.R 同一条报告线
STAGES <- c("Normal", "AAH", "AIS", "MIA", "LUAD")
PATIENTS <- sprintf("P%d", c(1:25))

stage_of <- function(s) sub("-\\d+$", "", sub("^GSM[0-9]+_P[0-9]+_", "", s))
q3 <- function(v) unname(quantile(v, c(.25, .5, .75), na.rm = TRUE))
tf <- function(v) {
  nm <- sprintf("ge_%.2f", TIERS)
  if (!length(v)) return(setNames(rep(NA_real_, length(TIERS)), nm))
  setNames(vapply(TIERS, function(t) mean(v >= t, na.rm = TRUE), numeric(1)), nm)
}

## margin = 该 spot 六谱系权重的 top1 − top2（与 RCTD_PREREG §13 一致）；只为**已落盘的上皮 spot** 算
## 返回按 bcs 命名的向量（**同一切片内** barcode 唯一；跨切片会重名 ⇒ 调用方须逐切片调用）
margin_of <- function(slide, bcs) {
  f <- file.path(WTS, sprintf("%s.weights.tsv.gz", slide))
  if (!file.exists(f)) return(setNames(rep(NA_real_, length(bcs)), bcs))
  w <- as.matrix(read.delim(gzfile(f), row.names = 1, check.names = FALSE))
  m <- apply(w, 1, function(z) { z <- sort(z, decreasing = TRUE); z[1] - z[2] })
  m[bcs]   # apply() 已按 rownames(w) 命名 ⇒ 按 barcode 取名即对齐
}

## 逐切片取 margin，按 (slide, barcode) 对齐到向量下标 i
## 🔴 `o$bc` 落盘时**带切片前缀**（`<slide>_<barcode>`，见 19_*.R 的 build_obj）⇒ 查权重表前必须剥掉
margin_for_idx <- function(o, i) {
  mm <- rep(NA_real_, length(i))
  for (s in unique(o$slide[i])) {
    j <- which(o$slide[i] == s)
    raw <- sub(paste0("^", s, "_"), "", o$bc[i][j])
    m <- margin_of(s, raw)
    if (anyNA(m)) stop(sprintf("margin 对齐失败：%s（%d/%d 个 barcode 在权重表里找不到）",
                               s, sum(is.na(m)), length(m)), call. = FALSE)
    mm[j] <- as.numeric(m)
  }
  mm
}

## —————————————————————————————————————————————————————————————
## 一、逐患者读取
## —————————————————————————————————————————————————————————————
per_pat <- list(); obs <- list(); ctrl <- list(); darm <- list(); sc3 <- list(); miss <- character(0)

for (P in PATIENTS) {
  mf <- file.path(BASE, P, "manifest.json")
  if (!file.exists(mf)) { miss <- c(miss, P); next }
  M <- fromJSON(mf)
  per_pat[[P]] <- M
  g <- function(f) { p <- file.path(BASE, P, f); if (file.exists(p)) readRDS(p) else NULL }
  obs[[P]]  <- g(sprintf("main_%s_spot_scores.rds", M$main_kind))
  ctrl[[P]] <- g("borrow_ctrl_spot_scores.rds")
  darm[[P]] <- g("D_arm_spot_scores.rds")
  sc3[[P]]  <- g("SC3_holdout.rds")
}
cat(sprintf("读到 %d 例；缺 %s\n", length(obs), if (length(miss)) paste(miss, collapse = ",") else "无"))

## —————————————————————————————————————————————————————————————
## 二、逐患者 SC4（§19.7）
## —————————————————————————————————————————————————————————————
rows <- list()
for (P in names(obs)) {
  o <- obs[[P]]; if (is.null(o)) next
  for (s in sort(unique(o$slide))) {
    i <- which(o$slide == s)
    rows[[length(rows) + 1L]] <- data.frame(
      patient = P, sex = per_pat[[P]]$sex, main_kind = per_pat[[P]]$main_kind,
      ref_source = per_pat[[P]]$sc4$ref_source, slide = s, stage = stage_of(s),
      group = paste(unique(o$grp[i]), collapse = "|"), n_spot = length(i),
      q25 = q3(o$cf[i])[1], median = q3(o$cf[i])[2], q75 = q3(o$cf[i])[3],
      max = max(o$cf[i], na.rm = TRUE), t(tf(o$cf[i])), row.names = NULL)
  }
}
per_slide <- do.call(rbind, rows)
write.table(per_slide, file.path(OUT, "sc4_per_slide.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## —————————————————————————————————————————————————————————————
## 三、按期别汇总（**pooled 与 逐患者 并列报**；AAH/Normal 无观测切片）
## —————————————————————————————————————————————————————————————
stage_rows <- list(); pat_rows <- list()
for (P in names(obs)) {
  o <- obs[[P]]; if (is.null(o)) next
  io <- which(o$grp == "observation")
  if (!length(io)) next
  st <- stage_of(o$slide)[io]
  for (s in intersect(STAGES, unique(st))) {
    v <- o$cf[io][st == s]
    pat_rows[[length(pat_rows) + 1L]] <- data.frame(
      patient = P, stage = s, n_spot = length(v), median = median(v, na.rm = TRUE),
      t(tf(v)), row.names = NULL)
  }
}
per_pat_stage <- do.call(rbind, pat_rows)
write.table(per_pat_stage, file.path(OUT, "stage_by_patient.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

allobs <- do.call(rbind, lapply(names(obs), function(P) {
  o <- obs[[P]]; if (is.null(o)) return(NULL)
  i <- which(o$grp == "observation")
  data.frame(patient = P, slide = o$slide[i], stage = stage_of(o$slide[i]),
             bc = if (is.null(o$bc)) NA_character_ else o$bc[i], cf = o$cf[i],
             row.names = NULL)
}))
for (s in STAGES) {
  i <- allobs$stage == s
  v <- allobs$cf[i]; pts <- unique(allobs$patient[i])
  med_pt <- vapply(pts, function(p) median(allobs$cf[allobs$patient == p & allobs$stage == s],
                                           na.rm = TRUE), numeric(1))
  stage_rows[[length(stage_rows) + 1L]] <- data.frame(
    stage = s, n_slide = length(unique(allobs$slide[i])), n_patient = length(pts),
    n_spot = sum(i), median = if (sum(i)) median(v, na.rm = TRUE) else NA_real_,
    t(tf(v)),
    pt_median_min = if (length(med_pt)) min(med_pt) else NA_real_,
    pt_median_max = if (length(med_pt)) max(med_pt) else NA_real_,
    n_pt_ge_0.05 = sum(med_pt >= TH_MAIN), row.names = NULL)
}
stage_tab <- do.call(rbind, stage_rows)
write.table(stage_tab, file.path(OUT, "stage_pooled.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## —————————————————————————————————————————————————————————————
## 四、自带锚 vs 借锚对照（§19.10 第 2 项：10 例）
## —————————————————————————————————————————————————————————————
svb <- do.call(rbind, lapply(names(ctrl), function(P) {
  if (is.null(ctrl[[P]]) || is.null(obs[[P]])) return(NULL)
  f <- function(x) { i <- x$grp == "observation"
                     c(median = median(x$cf[i], na.rm = TRUE), tf(x$cf[i])) }
  data.frame(patient = P, arm = c("self", "borrow"),
             rbind(f(obs[[P]]), f(ctrl[[P]])), row.names = NULL)
}))
if (!is.null(svb)) write.table(svb, file.path(OUT, "self_vs_borrow.tsv"), sep = "\t",
                               quote = FALSE, row.names = FALSE)

## —————————————————————————————————————————————————————————————
## 五、D 对照锚（reference = 该患者全部切片的非上皮）
## —————————————————————————————————————————————————————————————
dr <- do.call(rbind, lapply(names(darm), function(P) {
  d <- darm[[P]]; if (is.null(d)) return(NULL)
  i <- d$grp == "observation"
  data.frame(patient = P, arm = "D_nonepi", n_ref = sum(d$grp == "reference"),
             median = median(d$cf[i], na.rm = TRUE), t(tf(d$cf[i])), row.names = NULL)
}))
if (!is.null(dr)) write.table(dr, file.path(OUT, "d_arm.tsv"), sep = "\t",
                              quote = FALSE, row.names = FALSE)

## —————————————————————————————————————————————————————————————
## 六、§13.7.3 margin 梯度五档（**唯一的新计算，且是已签义务**）
##    主结论的操作化：**LUAD 观测上皮 − 前驱期(AIS/MIA)观测上皮 的 cf 中位差**，按患者配对。
##    五档全报；梯度单调 ⇒「稳定」；非单调或翻向 ⇒「依赖并列判定」。
## —————————————————————————————————————————————————————————————
mg_pat <- list()
for (P in names(obs)) {
  o <- obs[[P]]; if (is.null(o) || is.null(o$bc)) next
  i <- which(o$grp == "observation")
  if (!length(i)) next
  st <- stage_of(o$slide)[i]
  if (!any(st == "LUAD") || !any(st %in% c("AIS", "MIA"))) next
  mm <- margin_for_idx(o, i)                      # 按 (slide, barcode) 对齐，**不是**按 barcode 单独对齐
  for (g in MGRID) {
    k <- which(!is.na(mm) & mm >= g)
    v <- o$cf[i][k]; s2 <- st[k]
    d <- median(v[s2 == "LUAD"], na.rm = TRUE) - median(v[s2 %in% c("AIS", "MIA")], na.rm = TRUE)
    ## 逐档**同时**记各期别自己的中位（这样换一个对比也能直接从本表读，不必重算）
    mg_pat[[length(mg_pat) + 1L]] <- data.frame(
      patient = P, margin_grid = g, n_used = length(k), delta = d,
      n_AIS = sum(s2 == "AIS"), n_MIA = sum(s2 == "MIA"), n_LUAD = sum(s2 == "LUAD"),
      med_AIS  = if (any(s2 == "AIS"))  median(v[s2 == "AIS"],  na.rm = TRUE) else NA_real_,
      med_MIA  = if (any(s2 == "MIA"))  median(v[s2 == "MIA"],  na.rm = TRUE) else NA_real_,
      med_LUAD = if (any(s2 == "LUAD")) median(v[s2 == "LUAD"], na.rm = TRUE) else NA_real_,
      row.names = NULL)
  }
}
if (length(mg_pat)) {
  mp <- do.call(rbind, mg_pat)
  write.table(mp, file.path(OUT, "margin_gradient_by_patient.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  agg <- do.call(rbind, lapply(MGRID, function(g) {
    z <- mp[mp$margin_grid == g, ]
    data.frame(margin_grid = g, n_patient = nrow(z), n_spot = sum(z$n_used),
               delta_median = median(z$delta, na.rm = TRUE),
               n_pt_positive = sum(z$delta > 0, na.rm = TRUE),
               n_pt_negative = sum(z$delta < 0, na.rm = TRUE), row.names = NULL)
  }))
  d <- agg$delta_median
  ## 单调（同向且不回头）⇒ 稳定；否则 ⇒ 待定（§13.7.3 第 3 条）
  mono <- all(diff(d) >= -1e-9) || all(diff(d) <= 1e-9)
  agg$verdict <- if (mono) "稳定（单调）" else "依赖并列判定 ⇒ 结论降级为待定"
  write.table(agg, file.path(OUT, "margin_gradient.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
} else { agg <- NULL }

## —————————————————————————————————————————————————————————————
## 七、图（**English only**：本机无 CJK 字体）
## —————————————————————————————————————————————————————————————
png(file.path(OUT, "cohort_report.png"), width = 1600, height = 700, res = 130)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3, 1), las = 1)
boxplot(cf ~ stage, data = allobs[allobs$stage %in% c("AIS", "MIA", "LUAD"), ],
        col = "steelblue", outline = FALSE, ylim = c(0, 0.4),
        main = "Observation epithelium: cnv_fraction by stage",
        xlab = "stage", ylab = "cnv_fraction")
abline(h = TH_MAIN, lty = 2, col = "grey40")
mtext("dashed = 0.05 (report line, not a filter)", side = 3, line = 0.2, cex = 0.7, adj = 1)
if (!is.null(agg)) {
  plot(agg$margin_grid, agg$delta_median, type = "b", pch = 19, col = "firebrick",
       ylim = range(c(0, mp$delta), na.rm = TRUE), xlab = "RCTD margin threshold",
       ylab = "median(LUAD) - median(AIS/MIA)", main = "Margin gradient (5 tiers, per patient)")
  for (p in unique(mp$patient)) {
    z <- mp[mp$patient == p, ]; lines(z$margin_grid, z$delta, col = "grey70", lty = 3)
  }
  lines(agg$margin_grid, agg$delta_median, type = "b", pch = 19, col = "firebrick", lwd = 2)
  abline(h = 0, lty = 2, col = "grey40")
}
par(op); dev.off()

## —————————————————————————————————————————————————————————————
## 八、本上报自己的清单（含 not_claimed）
## —————————————————————————————————————————————————————————————
rep_json <- list(
  script = "08_spatial_deconv/20_report_cohort.R",
  n_patient_reported = length(obs), missing = miss,
  tiers = TIERS, margin_grid = MGRID, th_main = TH_MAIN,
  stage_pooled = stage_tab,
  margin_gradient = agg,
  margin_verdict = if (!is.null(agg)) unique(agg$verdict) else "未算（无同时含 LUAD 与前驱期的患者）",
  not_claimed = c(
    "不产生恶性判定（cnv_fraction 是连续量；CNVClassification 只是方向分类器）",
    "观测期别只有 AIS/MIA/LUAD：AAH 全部是参考切片、Normal 全队列只有 P4 一张 ⇒ 不得读成五档病程轴",
    "参考组由构造≈0 ⇒「参考平」不是证据",
    "不声称借锚与自带锚等价（10 例对照只作并列报）",
    "margin 五档是报告网格，不是口径：不得用它们筛数据、建锚或改样本量",
    "四档 0.02/0.05/0.10/0.15 只报比例，不设通过率阈值",
    "本上报是汇总，不是结论"))
write_json(rep_json, file.path(OUT, "cohort_report.json"), pretty = TRUE, auto_unbox = TRUE)

cat("\n==== 期别汇总（观测上皮）====\n"); print(stage_tab)
if (!is.null(agg)) { cat("\n==== margin 梯度五档 ====\n"); print(agg) }
cat(sprintf("\n产物：%s\n", OUT))
