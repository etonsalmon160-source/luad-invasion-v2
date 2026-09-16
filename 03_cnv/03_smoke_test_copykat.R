#!/usr/bin/env Rscript
# 03_cnv/03_smoke_test_copykat.R <sample_id>
#
# GP1 · CopyKAT 冒烟测试：单样本跑一次，量出 耗时 / 峰值内存 / 保留细胞数 / CNV 判定分布，
# 供外推全量 75 样本排期。参数为本项目预注册值，不因样本大小而调整。
#
# 用法：Rscript 03_cnv/03_smoke_test_copykat.R P13_Normal

suppressMessages({library(copykat)})
source("/home/eto/luad_v2/03_cnv/00_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("用法: Rscript 03_smoke_test_copykat.R <sample_id>")
sample_id <- args[1]

OUT_ROOT <- "/home/eto/luad_v2/results/03_cnv/smoke"
out_dir  <- file.path(OUT_ROOT, sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. 取矩阵 --------------------------------------------------------------
t0 <- proc.time()
M <- read_h5ad_sample(sample_id)
n_cells <- nrow(M); n_genes_raw <- ncol(M)
rss_after_load <- peak_rss_kb()

fate <- hg20_gene_fate(colnames(M))
n_kept <- length(fate$kept); n_drop <- length(fate$dropped)

# 读取路径的独立证伪：行和/非零数须与 M1 权威表精确相等，否则宁可停也不带错数往下跑
xchk <- crosscheck_against_m1(M, sample_id)
if (xchk$max_abs_dev_nCount != 0 || xchk$max_abs_dev_nFeature != 0)
  stop(sprintf("与 M1 交叉核对未通过: nCount偏差=%s nFeature偏差=%s",
               xchk$max_abs_dev_nCount, xchk$max_abs_dev_nFeature))

# 静默丢弃清单落盘，便于逐条回溯
writeLines(c(sprintf("# sample=%s  genes_raw=%d  hg20_kept=%d  dropped=%d (%.2f%%)",
                     sample_id, n_genes_raw, n_kept, n_drop, 100 * n_drop / n_genes_raw),
             sprintf("# full.anno 唯一 hgnc_symbol = %d", fate$n_anno_symbols),
             fate$dropped),
           file.path(out_dir, "hg20_dropped_genes.txt"))

# ---- 复现 copykat 内部的静默改参判定（源码第 9-29 行）------------------------
# 必须按 copykat 自己的**顺序**复算，否则判据会错：
#   先按 min.gene.per.cell 丢细胞 → 再对**全部 18,082 基因**按 LOW.DR 丢基因
#   → 若剩下的基因数 < 7000，则 UP.DR 被**静默改成 LOW.DR**（第 27 行）。
# 注意：判据用的是 LOW.DR 过滤**之后**的基因数，不是注释前的 17,361。
cell_ngenes <- Matrix::rowSums(M > 0)
keep_cells  <- cell_ngenes >= 200
n_cells_dropped <- sum(!keep_cells)
gene_der <- if (any(keep_cells))
  Matrix::colSums(M[keep_cells, , drop = FALSE] > 0) / sum(keep_cells) else numeric(ncol(M))
n_pass_lowdr <- sum(gene_der > 0.05)
low_quality_override <- n_pass_lowdr < 7000
up_dr_effective <- if (low_quality_override) 0.05 else 0.1

# ---- 2. 转 CopyKAT 要求的 基因×细胞 稠密矩阵 --------------------------------
mat_gc <- as.matrix(t(M))          # t(dgRMatrix) → dgCMatrix（O(1)），再落稠密
rm(M); invisible(gc())
rss_after_dense <- peak_rss_kb()

# ---- 3. 调用（参数全部预注册，不随样本调整）--------------------------------
cat(sprintf("[GP1] %s: %d cells × %d genes(hg20)  dense=%.2f GB  rss=%.2f GB\n",
            sample_id, n_cells, n_kept, prod(dim(mat_gc)) * 8 / 2^30,
            rss_after_dense / 2^20))

old <- setwd(out_dir)
sink_log <- file.path(out_dir, "copykat.stdout.log")
sink(sink_log, split = TRUE)
t_copy <- system.time({
  res <- tryCatch(
    copykat(rawmat = mat_gc, id.type = "S", cell.line = "no",
            ngene.chr = 5, min.gene.per.cell = 200,
            LOW.DR = 0.05, UP.DR = 0.1, win.size = 25,
            KS.cut = 0.1, distance = "euclidean", genome = "hg20",
            n.cores = 1, sam.name = sample_id,
            output.seg = "FALSE", plot.genes = "FALSE"),
    error = function(e) structure(conditionMessage(e), class = "copykat_error"))
})
sink()   # sink(<路径字符串>) 的对应关闭就是 sink()；不可再 close(路径)
setwd(old)
wall_copykat <- unname(t_copy["elapsed"])
rss_final <- peak_rss_kb()

err <- if (inherits(res, "copykat_error")) { e <- as.character(res); class(e) <- NULL; e } else NULL

# ---- 4. 解析判定 ------------------------------------------------------------
# ⚠️ 计数 0 与"未测定"必须分清：table() 不含零水平，若某档（如 Normal 样本的 aneuploid）
# 真的一个都没有，字面量里就查不到它。故**一旦成功解析预测文件即全部归零**再按出现值填充；
# 仅解析失败时整体置 NA。曾把真值 0 误写成 NA，已改。
cnt <- list(aneuploid = NA_integer_, diploid = NA_integer_,
            not.defined = NA_integer_, NA_ = NA_integer_)
pred_file <- file.path(out_dir, sprintf("%s_copykat_prediction.txt", sample_id))
if (is.null(err) && file.exists(pred_file)) {
  pr <- read.delim(pred_file, stringsAsFactors = FALSE)
  pcol <- grep("copykat.pred", colnames(pr), value = TRUE)[1]
  tab <- table(pr[[pcol]], useNA = "always")
  cnt <- list(aneuploid = 0L, diploid = 0L, not_defined = 0L, NA_ = 0L)  # 解析成功→归零，见下方说明
  for (k in names(tab)) {
    kk <- if (is.na(k)) "NA_" else gsub("[^a-z.]", ".", tolower(k))
    if (kk == "aneuploid") cnt$aneuploid <- as.integer(tab[k])
    if (kk == "diploid") cnt$diploid <- as.integer(tab[k])
    if (kk == "not.defined") cnt$not.defined <- as.integer(tab[k])
    if (kk == "NA_") cnt$NA_ <- as.integer(tab[k])
  }
} else if (!is.null(res) && !is.null(res$prediction)) {
  tab <- table(res$prediction$copykat.pred, useNA = "always")
  cnt <- list(aneuploid = 0L, diploid = 0L, not_defined = 0L, NA_ = 0L)
  cnt$aneuploid   <- as.integer(sum(tab[grep("aneuploid", names(tab))]))
  cnt$diploid     <- as.integer(sum(tab[grep("diploid", names(tab))]))
  cnt$not_defined <- as.integer(sum(tab[grep("not\\.defined", names(tab))]))
}

n_pred <- sum(unlist(cnt)[1:3], na.rm = TRUE)
frac <- if (n_pred > 0) cnt$aneuploid / n_pred else NA_real_

json_write(list(
  sample_id = sample_id,
  r_version = as.character(getRversion()),
  copykat_version = as.character(packageVersion("copykat")),
  param_low_dr = 0.05,
  param_up_dr_passed = 0.1,
  param_genome = "hg20",
  param_ngene_chr = 5,
  param_min_gene_per_cell = 200,
  param_win_size = 25,
  param_ks_cut = 0.1,
  param_distance = "euclidean",
  param_h5ad_sha256 = H5AD_SHA256_GP0,
  n_cells = n_cells,
  n_genes_raw = n_genes_raw,
  n_genes_hg20 = n_kept,
  n_genes_dropped = n_drop,
  dropped_pct = 100 * n_drop / n_genes_raw,
  n_cells_dropped_below_200genes = n_cells_dropped,
  n_genes_pass_lowDR = n_pass_lowdr,
  low_quality_updr_override = low_quality_override,
  up_dr_effective = up_dr_effective,
  m1_max_abs_dev_nCount = xchk$max_abs_dev_nCount,
  m1_max_abs_dev_nFeature = xchk$max_abs_dev_nFeature,
  dense_gb = prod(dim(mat_gc)) * 8 / 2^30,
  wall_sec_load = unname((proc.time() - t0)["elapsed"]) - wall_copykat,
  wall_sec_copykat = wall_copykat,
  peak_rss_gb_after_load = rss_after_load / 2^20,
  peak_rss_gb_after_dense = rss_after_dense / 2^20,
  peak_rss_gb_final = rss_final / 2^20,
  n_pred_aneuploid = cnt$aneuploid,
  n_pred_diploid = cnt$diploid,
  n_pred_not_defined = cnt$not.defined,
  frac_aneuploid = frac,
  error = err,
  output_dir = out_dir
), file.path(OUT_ROOT, sprintf("%s.json", sample_id)))

cat(sprintf("[GP1] %s: 完成 copykat=%.1fs  峰值RSS=%.2f GB  aneuploid=%s/%s  err=%s\n",
            sample_id, wall_copykat, rss_final / 2^20,
            cnt$aneuploid, n_pred, if (is.null(err)) "无" else err))
