#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/06_sensitivity_doublet_rate.R — M1 敏感性：双体剔除比例的影响
#
# 动机：scDblFinder 用 dbr=NULL 按细胞数自估，大样本率高（可能被质疑）。
#       本脚本回答："若改用固定剔除比例，会改变多少细胞？"
#       —— 直接基于已保存的 doublet_score（无需重跑），比较：
#         (a) 主结果（scDblFinder 自动阈值）
#         (b) 固定 top-10% 分数
#         (c) 固定 top-5%  分数
# 输出：results/01_qc/gse308103_sensitivity_doublet_rate.csv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
OUT <- "/home/eto/luad_v2/results/01_qc"
pc <- fread(file.path(OUT, "gse308103_per_cell_qc.csv.gz"))

tested <- pc[doublet_class != "not_tested"]
res <- tested[, {
  s <- doublet_score
  o <- order(-s)
  n <- .N
  auto <- sum(doublet_class == "doublet")
  top10 <- round(0.10 * n); top5 <- round(0.05 * n)
  d10 <- length(intersect(o[seq_len(top10)], which(doublet_class == "doublet")))
  d5  <- length(intersect(o[seq_len(top5)],  which(doublet_class == "doublet")))
  .(n = n, auto_dbl = auto, auto_rate = round(100 * auto / n, 2),
    top10_dbl = top10, top5_dbl = top5,
    overlap_top10 = d10, overlap_top5 = d5)
}, by = sample_id][order(sample_id)]

fwrite(res, file.path(OUT, "gse308103_sensitivity_doublet_rate.csv"))

cat("=== 双体剔除比例敏感性（全队列）===\n")
cat(sprintf("  参与判定细胞      : %d\n", sum(res$n)))
cat(sprintf("  主结果(自动阈值)  : %d 双体 (%.2f%%)\n", sum(res$auto_dbl), 100*sum(res$auto_dbl)/sum(res$n)))
cat(sprintf("  固定 top-10%% 分数 : %d 双体 (10.00%%)\n", sum(res$top10_dbl)))
cat(sprintf("  固定 top-5%%  分数 : %d 双体 ( 5.00%%)\n", sum(res$top5_dbl)))
cat(sprintf("\n  主结果 ∩ top-10%% 重合: %d / %d (%.1f%% of 主结果)\n",
            sum(res$overlap_top10), sum(res$auto_dbl),
            100*sum(res$overlap_top10)/sum(res$auto_dbl)))
cat(sprintf("  主结果 ∩ top-5%%  重合: %d / %d (%.1f%% of 主结果)\n",
            sum(res$overlap_top5), sum(res$auto_dbl),
            100*sum(res$overlap_top5)/sum(res$auto_dbl)))
cat(sprintf("\n  即：若改用固定 10%% 剔除，约 %.1f%% 的主结果双体不在其中（差异上界）\n",
            100 - 100*sum(res$overlap_top10)/sum(res$auto_dbl)))
