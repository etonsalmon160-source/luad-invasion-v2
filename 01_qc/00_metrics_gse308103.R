#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/00_metrics_gse308103.R — M1 第 1 步：只算 QC 指标，不设阈值
#
# 为什么先不用阈值：GSE308103 是 **snRNA（细胞核）**，其 nCount/nFeature 分布
# 与整细胞 scRNA 不同；照搬 scRNA 阈值（nFeature 500-10000 / nCount 1000-60000）
# 没有依据。故先据实统计分布，再由数据定阈值（下一步）。
#
# 输入：/home/eto/luad_invasion/data/GSE308103/extracted/<GSM>_<Pxx>_<Stage>.raw_counts.mtx.txt.gz
#       格式：基因 × 细胞（首行=条码，无基因列名）；基因列为 **符号**（可判 MT-）
# 输出：results/01_qc/gse308103_per_cell_metrics.csv.gz
#       results/01_qc/gse308103_metric_summary.csv
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
set.seed(1)

RAW <- "/home/eto/luad_invasion/data/GSE308103/extracted"
OUT <- "/home/eto/luad_v2/results/01_qc"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(RAW, pattern = "\\.raw_counts"), method = "radix")
cat(sprintf("[M1/308103] %d 个样本文件\n", length(files)))

# 读首行（条码行）。**必须单独读**：fread 遇到首行字段数(2631)与数据行(2632)不一致时
# 会 autostart 跳过首行 → 条码丢失、基因错位（曾致 cell_barcode 全为 "0"）。
read_barcode_line <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
}

res <- vector("list", length(files))
for (i in seq_along(files)) {
  f <- files[i]
  sid <- sub("^GSM[0-9]+_", "", sub("\\.raw_counts.*$", "", f))
  bc <- read_barcode_line(file.path(RAW, f))
  dt <- fread(cmd = sprintf("zcat %s/%s", RAW, f), header = FALSE, skip = 1)
  g  <- as.character(unlist(dt[[1]]))
  stopifnot(length(bc) == ncol(dt) - 1L, nrow(dt) == length(g))
  m  <- as.matrix(dt[, -1]); storage.mode(m) <- "double"   # skip=1 后首行已是基因 → 不再丢行
  cs <- colSums(m); fs <- colSums(m > 0)
  mt <- grepl("^MT-", g)
  pct <- if (any(mt)) 100 * colSums(m[mt, , drop = FALSE]) / pmax(cs, 1) else rep(0, ncol(m))
  res[[i]] <- data.table(
    dataset = "GSE308103", sample_id = sid,
    cell_barcode = paste0(bc, "|", sid),
    nCount = cs, nFeature = fs, pct_mt = round(pct, 3), n_genes_total = nrow(m)
  )
  cat(sprintf("  [%2d/%d] %-12s %5d genes x %6d cells  median(nCount)=%6.0f median(nFeature)=%5.0f median(mt%%)=%.1f\n",
              i, length(files), sid, nrow(m), ncol(m), median(cs), median(fs), median(pct))); flush.console()
  rm(dt, m); gc(verbose = FALSE)
}

per_cell <- rbindlist(res)
fwrite(per_cell, file.path(OUT, "gse308103_per_cell_metrics.csv.gz"))

per_sample <- per_cell[, .(
  n_cells = .N,
  q_nCount   = paste(round(quantile(nCount,   c(.01,.05,.25,.5,.75,.95,.99))), collapse="|"),
  q_nFeature = paste(round(quantile(nFeature, c(.01,.05,.25,.5,.75,.95,.99))), collapse="|"),
  q_pct_mt   = paste(round(quantile(pct_mt,   c(.01,.05,.25,.5,.75,.95,.99)),2), collapse="|")
), by = .(sample_id)][order(sample_id)]
fwrite(per_sample, file.path(OUT, "gse308103_metric_summary.csv"))

cat("\n=== 全队列分位数（1/5/25/50/75/95/99%）===\n")
cat("nCount  :", paste(round(quantile(per_cell$nCount, c(.01,.05,.25,.5,.75,.95,.99))), collapse=" "), "\n")
cat("nFeature:", paste(round(quantile(per_cell$nFeature,c(.01,.05,.25,.5,.75,.95,.99))), collapse=" "), "\n")
cat("pct_mt  :", paste(round(quantile(per_cell$pct_mt,  c(.01,.05,.25,.5,.75,.95,.99)),2), collapse=" "), "\n")
cat(sprintf("\n总细胞: %d   样本: %d\n", nrow(per_cell), uniqueN(per_cell$sample_id)))
cat(sprintf("产物: %s/gse308103_per_cell_metrics.csv.gz\n", OUT))
