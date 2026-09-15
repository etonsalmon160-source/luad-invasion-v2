#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/08_validate_doublet_calls.R — M1 双体判定的【内部正向验证】
#
# 背景：第二方法交叉验证在本环境不可行（scrublet 不适配稀疏核；DoubletFinder
#       需 Seurat 2/3 或 5，本机为 4.3.0）。改用**无需第二工具**的正向验证：
#
#   原理：真双体 = 两个细胞并成一个液滴 → **counts / 基因数 应显著高于单细胞**。
#   检查：对每个样本，比较 scDblFinder 判为 doublet vs singlet 的
#         median(nCount)、median(nFeature)，期望 **doublet > singlet**（比值 > 1）。
#   另查：双体率是否随样本细胞数上升（10x 载量关系）。
#
# 输出：results/01_qc/gse308103_doublet_call_validation.csv + 控制台摘要
# =============================================================================
suppressPackageStartupMessages(library(data.table))
OUT <- "/home/eto/luad_v2/results/01_qc"
pc <- fread(file.path(OUT, "gse308103_per_cell_qc.csv.gz"))
t <- pc[doublet_class != "not_tested" & qc_pass == TRUE]

v <- t[, {
  d <- doublet_class == "doublet"
  if (sum(d) >= 10 && sum(!d) >= 10) {
    .(n = .N, n_dbl = sum(d),
      med_cnt_dbl = as.numeric(median(nCount[d])),   med_cnt_sng = as.numeric(median(nCount[!d])),
      med_ftr_dbl = as.numeric(median(nFeature[d])), med_ftr_sng = as.numeric(median(nFeature[!d])))
  } else .(n = .N, n_dbl = sum(d), med_cnt_dbl = NA_real_, med_cnt_sng = NA_real_,
            med_ftr_dbl = NA_real_, med_ftr_sng = NA_real_)
}, by = sample_id][order(sample_id)]

v[, `:=`(cnt_ratio = round(med_cnt_dbl / med_cnt_sng, 3),
         ftr_ratio = round(med_ftr_dbl / med_ftr_sng, 3))]
fwrite(v, file.path(OUT, "gse308103_doublet_call_validation.csv"))

ok <- v[!is.na(cnt_ratio)]
cat("=== 双体判定正向验证：doublet/singlet 计数比 ===\n")
cat(sprintf("  可评估样本: %d/%d\n", nrow(ok), nrow(v)))
cat(sprintf("  nCount 比值  : 中位 %.2f  范围 %.2f–%.2f\n",
            median(ok$cnt_ratio), min(ok$cnt_ratio), max(ok$cnt_ratio)))
cat(sprintf("  nFeature 比值: 中位 %.2f  范围 %.2f–%.2f\n",
            median(ok$ftr_ratio), min(ok$ftr_ratio), max(ok$ftr_ratio)))
cat(sprintf("  比值 > 1 的样本（符合真双体预期）: %d/%d (%.0f%%)\n",
            sum(ok$cnt_ratio > 1), nrow(ok), 100 * mean(ok$cnt_ratio > 1)))
cat(sprintf("  比值 > 1.5 的样本: %d/%d\n", sum(ok$cnt_ratio > 1.5), nrow(ok)))

# 双体率 vs 细胞数（10x 载量关系）
agg <- t[, .(n_cells = .N, rate = mean(doublet_class == "doublet")), by = sample_id]
r <- cor(agg$n_cells, agg$rate)
cat(sprintf("\n  双体率 vs 样本细胞数 相关(Pearson r) = %.3f （预期为正）\n", r))
