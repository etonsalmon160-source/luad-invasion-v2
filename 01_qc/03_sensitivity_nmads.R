#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/03_sensitivity_nmads.R — M1 敏感性对照：nmads = 3 vs 5
#
# 目的：证明 QC 通过率对 MAD 阈值不过度敏感（阈值选择的稳健性检查）。
# 输入：results/01_qc/gse308103_per_cell_metrics.csv.gz（已算好的 per-cell 指标）
# 输出：results/01_qc/gse308103_sensitivity_nmads.csv + 控制台摘要
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(scuttle) })
IN  <- "/home/eto/luad_v2/results/01_qc/gse308103_per_cell_metrics.csv.gz"
OUT <- "/home/eto/luad_v2/results/01_qc"
MAX_MT <- 5

d <- fread(IN)
res <- d[, {
  p3 <- !isOutlier(nCount,   nmads = 3, type = "both", log = TRUE) &
        !isOutlier(nFeature, nmads = 3, type = "both", log = TRUE) & (pct_mt < MAX_MT)
  p5 <- !isOutlier(nCount,   nmads = 5, type = "both", log = TRUE) &
        !isOutlier(nFeature, nmads = 5, type = "both", log = TRUE) & (pct_mt < MAX_MT)
  .(n = .N, pass_nmads3 = sum(p3), pass_nmads5 = sum(p5),
    rate3 = round(100 * mean(p3), 2), rate5 = round(100 * mean(p5), 2))
}, by = sample_id][order(sample_id)]

fwrite(res, file.path(OUT, "gse308103_sensitivity_nmads.csv"))

cat("=== nmads 敏感性（全队列）===\n")
cat(sprintf("  总细胞      : %d\n", sum(res$n)))
cat(sprintf("  nmads=3 pass: %d (%.2f%%)\n", sum(res$pass_nmads3), 100*sum(res$pass_nmads3)/sum(res$n)))
cat(sprintf("  nmads=5 pass: %d (%.2f%%)\n", sum(res$pass_nmads5), 100*sum(res$pass_nmads5)/sum(res$n)))
cat(sprintf("  差异        : +%d 细胞 (%.2f 个百分点)\n",
            sum(res$pass_nmads5) - sum(res$pass_nmads3),
            100*(sum(res$pass_nmads5)-sum(res$pass_nmads3))/sum(res$n)))
cat(sprintf("\n  逐样本通过率差异：中位 %.2f pp，最大 %.2f pp\n",
            median(res$rate5 - res$rate3), max(abs(res$rate5 - res$rate3))))
