#!/usr/bin/env Rscript
# 03_cnv/05_routeA_coverage_floor.R
#
# 路线 A · 覆盖度地板的**负对照下采样实验**
#
# 目的：把 weak_cnv 的阈值 X 定在一个**不含 CNV 结果**的量上。
#
# 为什么需要它：copykat 的 not.defined 不是"算法对生物学不确定"。读源码
# tools/copykat/R/copykat.R:86-105 与 194-213 可知，它由三个**覆盖度**判据把细胞踢出分段：
#   ① 该细胞检出的（已注释）基因 < 5
#   ② 检出基因覆盖到的染色体数 < 23
#   ③ 任一染色体上的连续基因段长度 < ngene.chr(=5)
# 实测（P13_Normal）：not.defined 组 nFeature 中位 613 vs diploid 组 1007（0.61 倍），
# nCount 0.54 倍，而 pct_mt 完全一致（0.13 vs 0.12）→ 纯文库复杂度问题，非濒死细胞。
#
# 因此**样本级 not.defined 率被细胞组成混杂**（免疫/基质细胞天然低复杂度），而组成正是
# 沿 Normal→AAH→AIS→MIA→IAC 变化的那个轴。若对合并队列做离群规则，会优先命中某一分期的
# 样本 —— 最坏情况是把 Normal 对照判成失败样本，而那恰是整条梯度检验的锚。
#
# 本实验：取一个**已知二倍体**的负对照（P13_Normal 中被 copykat 判为 diploid 的细胞），
# 用二项稀释模拟不同测序深度，量出"深度降到多少时 copykat 开始把二倍体细胞判成
# not.defined"。这条曲线测的是**技术退化点**，全程没有碰癌症梯度 → 可预注册。
#
# ---------------------------------------------------------------------------
# 预注册的判定规则（**在看到曲线之前写定**，法则 3.2）
#   1. 全曲线逐点上报，不做平滑、不删点、不重跑。
#   2. f_knee := 使 not.defined 率**首次超过 0.50** 的最大稀释比例 f。
#      分母口径（在**看到任何曲线之前**写定，非事后选择）：分母 = 全部负对照细胞（494），
#      分子 = 其中被 copykat 标为 not.defined 者。被 min.gene.per.cell(=200) 整细胞丢弃的
#      负对照**不计入分子**（字面口径），其数量单独上报；若在膝点处很大，
#      须作为 OPEN ISSUE 上报，而**不得**静默改用 (not.defined + 丢弃)/N。
#      全样本 not.defined 率另列一栏仅作对照 —— 它把 307 个"全深度就已 not.defined"的
#      细胞也算进来，是一条被常数底噪抬高的曲线，不是本实验要量的东西。
#      0.50 是**本项目约定**（⚠️C），无文献出处；故整条曲线一并上报，
#      使读者能自行判断该切点移动时结论是否稳健。
#   3. 覆盖度地板 C* := 在 f_knee 处负对照细胞的 nFeature 中位数（稀释后矩阵，丢弃前）。
#   4. 样本级 weak_cnv 判据（供 GP2 队列用，本脚本不执行）：
#      observed_rate − predicted_rate（在深度对齐后用同一判据预测的率）= 生物学超额；
#      该超额须在**分期内**比较，不得对合并队列取分位数。
# ---------------------------------------------------------------------------
# 前置：03_smoke_test_copykat.R <sample_id> 已跑完 —— 需要它的 prediction 文件定义负对照。
#
# 用法：Rscript 03_cnv/05_routeA_coverage_floor.R [sample_id]     # 默认 P13_Normal

suppressMessages({library(copykat); library(rhdf5); library(Matrix); library(digest)})
source("/home/eto/luad_v2/03_cnv/00_common.R")

args <- commandArgs(trailingOnly = TRUE)
validate_only <- "--validate-only" %in% args
args <- setdiff(args, "--validate-only")
sample_id <- if (length(args) >= 1L) args[1] else "P13_Normal"

ROOT     <- "/home/eto/luad_v2"
OUT_ROOT <- file.path(ROOT, "results/03_cnv/routeA")
SMOKE    <- file.path(ROOT, "results/03_cnv/smoke")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# 与 03_smoke_test_copykat.R 完全一致的参数 —— 本实验唯一变化的量是测序深度
PARAM <- list(id.type = "S", cell.line = "no", ngene.chr = 5L, min.gene.per.cell = 200L,
              LOW.DR = 0.05, UP.DR = 0.1, win.size = 25L, KS.cut = 0.1,
              distance = "euclidean", genome = "hg20", n.cores = 1L)
# 稀释网格：1.00 是阳性对照（须复现原 run 的 diploid/not.defined 划分），其余等距下降
GRID <- c(1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20)
SEED_BASE <- 20260916L                 # 稀释随机种子基数（R5：可复现）

# ---- 1. 负对照：被 copykat 判为 diploid 的细胞 ------------------------------
pred_file <- file.path(SMOKE, sample_id, sprintf("%s_copykat_prediction.txt", sample_id))
if (!file.exists(pred_file))
  stop("找不到 ", pred_file, "\n  路线 A 需先跑完 03_cnv/03_smoke_test_copykat.R ", sample_id,
       " 以定义负对照（全深度下被判 diploid 的细胞）。")
pred_sha <- digest::digest(file = pred_file, algo = "sha256")
pr   <- read.delim(pred_file, stringsAsFactors = FALSE)
cls  <- setNames(pr[[grep("copykat.pred", colnames(pr), value = TRUE)[1]]], pr[["cell.names"]])
cat(sprintf("[routeA] 负对照来源 %s (sha256=%.12s…)\n", basename(pred_file), pred_sha))
print(table(cls))

if (any(grepl("aneuploid", names(table(cls)), ignore.case = TRUE)))
  stop("该样本原 run 存在 aneuploid 细胞 —— 它不是干净的负对照，硬停。")
neg_cells <- names(cls)[!grepl("not.defined", cls)]
if (length(neg_cells) < 200)
  stop("负对照细胞数 < 200，稀释后的子采样噪声会主导结果，硬停。")
cat(sprintf("[routeA] 负对照（全深度判 diploid）= %d 个细胞\n", length(neg_cells)))

# ---- 2. 取全深度计数 --------------------------------------------------------
M_all <- read_h5ad_sample(sample_id)
xchk  <- crosscheck_against_m1(M_all, sample_id)
if (xchk$max_abs_dev_nCount != 0 || xchk$max_abs_dev_nFeature != 0)
  stop("与 M1 交叉核对未通过，硬停")
if (!setequal(rownames(M_all), names(cls)))
  stop("h5ad 的细胞集与 prediction 文件的细胞集不等，硬停")
M0 <- M_all[match(names(cls), rownames(M_all)), , drop = FALSE]   # 行序对齐到 cls
genes_all <- colnames(M0)
rm(M_all); invisible(gc())

# ---- 3. 复现 copykat 的三个覆盖度判据（全深度）并自证还原 not.defined 集 --------
# 既是**复用**（GP2 队列要用同一判据），也是**证伪**：若重放不能精确还原
# diploid/not.defined 的划分，说明我对判据的理解有误，必须停而不是继续。
# 判定链本体在 00_common.R 的 copykat_chain()（**唯一实现**，02_prereg 共用同一份）。
n_annot <- length(intersect(genes_all, hg20_gene_fate(genes_all)$kept))
cat(sprintf("[routeA] 注释表覆盖 %d/%d 个基因（判据实际用的基因集在 LOW.DR 过滤之后才定）\n",
            n_annot, length(genes_all)))

cat("[routeA] 重放 copykat 判据（全深度）以自证…\n")
rp      <- copykat_chain(M0, low_dr = PARAM$LOW.DR, up_dr = PARAM$UP.DR,
                         ngene_chr = PARAM$ngene.chr,
                         min_gene_per_cell = PARAM$min.gene.per.cell)
obs_nd  <- names(cls)[grepl("not.defined", cls)]
pred_nd <- setdiff(names(cls), rp$survivors)
ok      <- setequal(pred_nd, obs_nd)
cat(sprintf("[routeA] 重放 not.defined=%d  实测=%d  集合相等=%s\n",
            length(pred_nd), length(obs_nd), ok))
cat(sprintf("[routeA] ToRemov2=%d  ToRemov3=%d  UP.DR有效值=%.2f\n",
            rp$n_rem2, rp$n_rem3, rp$up_dr_effective))
# 可核对量：copykat 自己写出的 *_raw_results_gene_by_cell.txt 行数 - 1 == anno.mat2 的基因数
# （= 两次基因过滤 + 注释 + 删周期/HLA 之后的**实际分段基因集**）。两者须一致，否则判据有偏。
cat(sprintf("[routeA] 判据基因集：LOW.DR 后 %d 个 → 注释并删周期/HLA 后 %d 个（%d 条染色体）→ DR2 过滤后 %d 个\n",
            rp$n_after_lowdr, nrow(rp$anno), rp$n_chromosomes, rp$n_genes_final))
# 独立证伪：与 copykat 自己落盘的 anno.mat2 行数逐字对表（存在才比，绝不拿模型数冒充实测）
gfc <- file.path(SMOKE, sample_id, sprintf("%s_copykat_raw_results_gene_by_cell.txt", sample_id))
if (file.exists(gfc)) {
  con <- file(gfc, "r"); n_meas <- 0L          # 分块计数，避免把 57 MB 整表读进内存
  while (length(ch <- readLines(con, n = 50000L))) n_meas <- n_meas + length(ch)
  close(con); n_meas <- n_meas - 1L
  cat(sprintf("[routeA] 对表 copykat 实测 anno.mat2 = %d 个基因；本实现 = %d；%s\n",
              n_meas, rp$n_genes_final,
              if (n_meas == rp$n_genes_final) "一致 ✔" else "**不一致 —— 判据有偏，硬停**"))
  if (n_meas != rp$n_genes_final)
    stop("判据基因集与 copykat 实测 anno.mat2 不一致，硬停。")
}
if (!ok) {
  writeLines(c("重放多出的细胞:", head(setdiff(pred_nd, obs_nd), 50),
               "实测多出的细胞:", head(setdiff(obs_nd, pred_nd), 50)),
             file.path(OUT_ROOT, "replay_mismatch.txt"))
  stop("判据重放未能还原 not.defined 集 —— 见 routeA/replay_mismatch.txt，硬停。")
}

# 逐细胞覆盖度指标落盘（供队列"预测率"复用，无需再跑 copykat）
# ⚠️ 基因集必须用 replay 内部定下的那一个（LOW.DR 过滤 + 注释 + 删周期/HLA），
# 否则 n_chrom 会系统性偏高，队列端的"预测 not.defined 率"随之偏低。
gi    <- match(rp$anno$hgnc_symbol, genes_all); stopifnot(!anyNA(gi))
Sub   <- M0[, gi, drop = FALSE]
chrom <- rp$anno$chromosome_name
write.csv(data.frame(
  cell_barcode  = rownames(M0),
  orig_class    = unname(cls),
  nFeature_all  = as.integer(Matrix::rowSums(M0 > 0)),
  nCount_all    = as.integer(Matrix::rowSums(M0)),
  nFeature_anno = as.integer(Matrix::rowSums(Sub > 0)),
  n_chrom       = vapply(seq_len(nrow(Sub)), function(i) length(unique(chrom[which(as.vector(Sub[i, ]) > 0)])), 0L),
  stringsAsFactors = FALSE),
  file.path(OUT_ROOT, sprintf("%s_per_cell_full_depth.csv", sample_id)), row.names = FALSE)

if (validate_only) {
  cat("[routeA] --validate-only：判据重放已自证通过，未跑稀释实验，退出。\n")
  quit(save = "no", status = 0L)
}

# ---- 4. 稀释并逐个跑 copykat ------------------------------------------------
curve <- list()
for (f in GRID) {
  tag  <- sprintf("%s_f%03d", sample_id, round(f * 100))
  fdir <- file.path(OUT_ROOT, tag)
  unlink(fdir, recursive = TRUE); dir.create(fdir, recursive = TRUE)

  set.seed(SEED_BASE + round(f * 1000))               # 各比例独立抽样、独立种子
  Ms <- M0
  Ms@x <- as.numeric(rbinom(length(M0@x), size = as.integer(M0@x), prob = f))
  n_before_drop <- nrow(Ms)
  Ms <- Ms[Matrix::rowSums(Ms) > 0, , drop = FALSE]
  n_cells_thinned <- nrow(Ms)                        # ⚠️ Md 转置后 nrow 是**基因**数，别拿它当细胞数
  Md <- as.matrix(t(Ms)); rm(Ms); invisible(gc())
  # ⚠️ Md 是 **基因×细胞**（copykat 的 rawmat 朝向）。逐细胞 nFeature 必须用 colSums；
  #    曾误用 rowSums，算出来的是"逐基因检出细胞数的中位数"，与"负对照细胞 nFeature 中位数"无关。
  nfeat_cell    <- colSums(Md > 0)
  nfeat_med     <- as.integer(median(nfeat_cell))
  nfeat_med_neg <- as.integer(median(nfeat_cell[intersect(neg_cells, names(nfeat_cell))]))

  cat(sprintf("\n[routeA] f=%.2f  cells=%d(全零丢 %d)  中位nFeature=%d(负对照 %d)  稠密=%.3f GB\n",
              f, n_cells_thinned, n_before_drop - n_cells_thinned, nfeat_med, nfeat_med_neg,
              prod(dim(Md)) * 8 / 2^30))

  old <- setwd(fdir)
  sink(file.path(fdir, "copykat.stdout.log"), split = TRUE)
  res <- tryCatch(system.time(
    copykat(rawmat = Md, id.type = PARAM$id.type, cell.line = PARAM$cell.line,
            ngene.chr = PARAM$ngene.chr, min.gene.per.cell = PARAM$min.gene.per.cell,
            LOW.DR = PARAM$LOW.DR, UP.DR = PARAM$UP.DR, win.size = PARAM$win.size,
            KS.cut = PARAM$KS.cut, distance = PARAM$distance, genome = PARAM$genome,
            n.cores = PARAM$n.cores, sam.name = tag,
            output.seg = "FALSE", plot.genes = "FALSE")),
    error = function(e) structure(conditionMessage(e), class = "copykat_error"))
  sink(); setwd(old)

  # ⚠️ 必须剥掉 class 再放进 data.frame()：data.frame(error = <class "copykat_error">)
  # 会走 as.data.frame.default 并抛 "cannot coerce class 'copykat_error' to a data.frame"，
  # 把"copykat 在这个深度失败"这件事伪装成脚本崩溃（实测 f=0.30 已踩）。
  err  <- if (inherits(res, "copykat_error")) { e <- as.character(res); class(e) <- NULL; e } else NULL
  wall <- if (is.null(err)) unname(res["elapsed"]) else NA_real_
  if (!is.null(err)) cat(sprintf("[routeA] f=%.2f  copykat 自身报错（非本脚本）：%s\n", f, err))

  # 解析成功即先归零再填充：0 与未测定语义相反，table() 不含零水平，曾把真值 0 误写成 NA
  cnt <- c(diploid = NA_integer_, aneuploid = NA_integer_, not_defined = NA_integer_)
  pf  <- file.path(fdir, sprintf("%s_copykat_prediction.txt", tag))
  cell_in_pred <- character(0)
  neg_nd <- NA_integer_
  if (is.null(err) && file.exists(pf)) {
    pr2  <- read.delim(pf, stringsAsFactors = FALSE)
    pcol <- grep("copykat.pred", colnames(pr2), value = TRUE)[1]
    t2   <- table(pr2[[pcol]])
    cnt[] <- 0L
    for (k in names(t2)) {
      kk <- gsub("[^a-z.]", ".", tolower(k))
      if (kk == "diploid")     cnt["diploid"]     <- as.integer(t2[k])
      if (kk == "aneuploid")   cnt["aneuploid"]   <- as.integer(t2[k])
      if (kk == "not.defined") cnt["not_defined"] <- as.integer(t2[k])
    }
    cell_in_pred <- as.character(pr2[["cell.names"]])
    pv           <- setNames(as.character(pr2[[pcol]]), cell_in_pred)
    neg_nd       <- sum(grepl("not.defined", pv[intersect(neg_cells, cell_in_pred)]))
  }
  n_pred <- sum(cnt, na.rm = TRUE)                    # copykat 实际尝试分段的细胞数
  n_def  <- sum(cnt[c("diploid", "aneuploid")], na.rm = TRUE)

  # 负对照受限口径 = **预注册主口径**（头文件规则 2/3：问的是"二倍体细胞何时被判成
  # not.defined"，不是"全样本 not.defined 有多少"）。全样本口径只是常数底噪的堆叠。
  n_neg          <- length(neg_cells)
  neg_absent     <- n_neg - sum(neg_cells %in% cell_in_pred)   # 被 <200 基因过滤整个丢掉
  rate_nd_neg    <- if (is.na(neg_nd)) NA_real_ else neg_nd / n_neg
  # 稳健性对照口径（不用于主判据）：把"整细胞被丢"也算作不可分段
  rate_unseg_neg <- if (is.na(neg_nd)) NA_real_ else (neg_nd + neg_absent) / n_neg

  # 中间产物（每样本约 283 MB 的 *_bin_by_cell / *_gene_by_cell）解析后即删
  keep <- c(sprintf("%s_copykat_prediction.txt", tag), "copykat.stdout.log")
  for (p in setdiff(list.files(fdir, full.names = TRUE), file.path(fdir, keep)))
    unlink(p, recursive = TRUE)

  row <- data.frame(
    dilut = f, n_cells_in = n_cells_thinned, n_cells_pred = n_pred,
    n_dropped_lowgenes = n_cells_thinned - n_pred,
    n_diploid = cnt["diploid"], n_aneuploid = cnt["aneuploid"],
    n_not_defined = cnt["not_defined"], n_defined = n_def,
    rate_not_defined = if (n_pred > 0) cnt["not_defined"] / n_pred else NA_real_,
    # ---- 预注册主口径（负对照受限，见上）----
    n_neg_control = n_neg, n_neg_in_pred = n_neg - neg_absent,
    n_neg_not_defined = neg_nd, n_neg_absent = neg_absent,
    rate_not_defined_neg = rate_nd_neg, rate_unsegmentable_neg = rate_unseg_neg,
    nFeature_median_thinned = nfeat_med, nFeature_median_neg = nfeat_med_neg,
    wall_sec = wall, error = if (is.null(err)) NA_character_ else err,
    stringsAsFactors = FALSE)
  curve[[length(curve) + 1L]] <- row
  fm <- function(v, d = 0L) if (is.na(v)) "—" else formatC(v, format = "f", digits = d)
  cat(sprintf("[routeA] f=%.2f  diploid=%s aneuploid=%s not.defined=%s | 负对照 nd=%s 丢=%d 率=%s  %s\n",
              f, fm(cnt["diploid"]), fm(cnt["aneuploid"]), fm(cnt["not_defined"]),
              fm(neg_nd), neg_absent, fm(rate_nd_neg, 3L),
              if (is.na(wall)) "—" else sprintf("%.0fs", wall)))
  rm(Md); invisible(gc())
}

df <- do.call(rbind, curve)
write.csv(df, file.path(OUT_ROOT, "curve.csv"), row.names = FALSE)

# ---- 5. 按预注册规则读出地板 ------------------------------------------------
ok_rows <- df[is.na(df$error), ]
# 主判据用**负对照受限**率（头文件规则 2）；全样本率仅作对照同时上报。
f_knee  <- if (any(ok_rows$rate_not_defined_neg > 0.5, na.rm = TRUE))
  max(ok_rows$dilut[which(ok_rows$rate_not_defined_neg > 0.5)]) else NA_real_
f_knee_allsample <- if (any(ok_rows$rate_not_defined > 0.5, na.rm = TRUE))
  max(ok_rows$dilut[which(ok_rows$rate_not_defined > 0.5)]) else NA_real_
C_star  <- if (!is.na(f_knee)) ok_rows$nFeature_median_neg[ok_rows$dilut == f_knee] else NA_real_
F1      <- ok_rows[ok_rows$dilut == 1.0, ]
# 膝点处被整细胞丢弃的负对照数：若很大，字面口径(不算 not.defined)会低估，
# 必须作为 OPEN ISSUE 上报，而不是静默换成 rate_unsegmentable_neg。
knee_absent <- if (!is.na(f_knee)) ok_rows$n_neg_absent[ok_rows$dilut == f_knee] else NA_integer_

json_write(list(
  sample_id = sample_id,
  r_version = as.character(getRversion()),
  copykat_version = as.character(packageVersion("copykat")),
  param_id_type = PARAM$id.type, param_cell_line = PARAM$cell.line,
  param_ngene_chr = PARAM$ngene.chr, param_min_gene_per_cell = PARAM$min.gene.per.cell,
  param_low_dr = PARAM$LOW.DR, param_up_dr_passed = PARAM$UP.DR,
  param_win_size = PARAM$win.size, param_ks_cut = PARAM$KS.cut,
  param_distance = PARAM$distance, param_genome = PARAM$genome,
  param_h5ad_sha256 = H5AD_SHA256_GP0,
  prediction_file_sha256 = pred_sha,
  n_negative_control_cells = length(neg_cells),
  replay_reproduced_not_defined = ok,
  replay_n_rem2 = rp$n_rem2, replay_n_rem3 = rp$n_rem3,
  replay_up_dr_effective = rp$up_dr_effective,
  replay_n_genes_after_lowdr = rp$n_after_lowdr,
  n_genes_criterion_set = nrow(rp$anno),
  n_chromosomes_criterion_set = length(unique(rp$anno$chromosome_name)),
  control_f1_rate_not_defined          = if (nrow(F1)) F1$rate_not_defined else NA_real_,
  control_f1_rate_not_defined_neg      = if (nrow(F1)) F1$rate_not_defined_neg else NA_real_,
  f_knee = f_knee, coverage_floor_nFeature = C_star,
  f_knee_allsample_secondary = f_knee_allsample,
  knee_absent_neg_at_f_knee = knee_absent,
  knee_convention = "负对照受限 rate_not_defined_neg > 0.50, 本项目约定(⚠️C)，无文献出处",
  primary_rate_denominator = "全部负对照细胞数(494)，被<200基因过滤丢弃者不计入分子(字面口径)",
  grid = GRID, seed_base = SEED_BASE,
  output_dir = OUT_ROOT
), file.path(OUT_ROOT, sprintf("%s_routeA.json", sample_id)))

cat("\n==== 路线 A 曲线 ====\n")
print(df[, c("dilut", "n_cells_in", "n_cells_pred", "n_diploid", "n_not_defined",
             "n_neg_not_defined", "n_neg_absent", "rate_not_defined_neg",
             "rate_unsegmentable_neg", "nFeature_median_neg")], row.names = FALSE)
cat(sprintf("\n阳性对照 f=1.00: 全样本 not.defined 率 = %s （原 run 实测 %d/%d = %.3f）；负对照受限率 = %s\n",
            if (nrow(F1)) sprintf("%.3f", F1$rate_not_defined) else "NA",
            sum(grepl("not.defined", cls)), length(cls), mean(grepl("not.defined", cls)),
            if (nrow(F1)) sprintf("%.3f", F1$rate_not_defined_neg) else "NA"))
cat(sprintf("f_knee(负对照受限) = %s ; 覆盖度地板 C* = %s 个基因\n",
            if (is.na(f_knee)) "未触及(全曲线率均 <= 0.50)" else sprintf("%.2f", f_knee),
            if (is.na(C_star)) "NA" else as.character(C_star)))
cat(sprintf("对照：全样本口径 f_knee = %s ; 膝点处被整细胞丢弃的负对照 = %s\n",
            if (is.na(f_knee_allsample)) "未触及" else sprintf("%.2f", f_knee_allsample),
            if (is.na(knee_absent)) "NA(无膝点)" else as.character(knee_absent)))
cat("\n写出:", file.path(OUT_ROOT, "curve.csv"), "\n")
