#!/usr/bin/env Rscript
# 03_cnv/02_prereg_gene_tiers.R
#
# 为 GP2 的 UP.DR 决策提供定量依据：逐样本算出在 copykat 自身的过滤**顺序**下，
# 最终进入分割（segmentation）的基因数在两档阈值下各是多少。
#
# 为什么必须这么做：copykat 第 57 行 `if(nrow(rawmat) < 7000) UP.DR <- LOW.DR` 是
# **无条件覆写**传入的 UP.DR；UP.DR 又名不副实，实为**下界**（第 189 行 `DR2 >= UP.DR`）。
# 因此 75 个样本被切成两档：过 7000 的用 UP.DR=0.1，不足的用 UP.DR=0.05=LOW.DR。
# 想让有效值全队列统一，**唯一**办法是令 LOW.DR == UP.DR（否则覆写不可避免）。
#
# ⚠️ 2026-09-16 重写：基因数不再用**近似模型**算，改为调用 00_common.R 的
#    `copykat_chain()` —— 与 05_routeA_coverage_floor.R **同一份实现**，并且已在
#    P13_Normal 上与 copykat 自己落盘的 anno.mat2 对表通过（5903 == 5903）。
#    旧版本漏了 `:74-81` 的 cycle/HLA 删除（约 564 个基因），基因数偏高约 16%。
#
# ⚠️ 2026-09-16 追加**覆盖度地板**参数：地板（`LUAD_COVERAGE_FLOOR`，默认 840）作用在
#    进 copykat 之前，故会改变 `n_cells_used` 与（可能的）7000 覆写状态。本脚本以地板值
#    给输出文件名加后缀，**不覆盖**无地板的历史产物 `prereg_gene_tiers.csv`。
#    设 `LUAD_COVERAGE_FLOOR=0` 可复现旧行为 —— 已用作本次改动的回归检验。
#
# 用法：Rscript 03_cnv/02_prereg_gene_tiers.R
#       LUAD_COVERAGE_FLOOR=0 Rscript 03_cnv/02_prereg_gene_tiers.R   # 复现无地板版本

suppressMessages({library(copykat); library(Matrix)})
source("/home/eto/luad_v2/03_cnv/00_common.R")

H5AD_READ <- H5AD
samples <- sort(unique(as.character(h5read(H5AD_READ, "obs/sample_id/categories"))))
stopifnot(length(samples) == 75L)

genes_all <- as.character(h5read(H5AD_READ, "var/_index"))   # 基因表对所有样本相同
fate <- hg20_gene_fate(genes_all)
cat(sprintf("[tiers] 基因表 %d 个；在 full.anno 内 %d 个；被静默丢弃 %d 个（全样本相同）\n",
            length(genes_all), length(fate$kept), length(fate$dropped)))

MIN_GENE_PER_CELL <- 200L ; LOW_DR <- 0.05 ; UP_DR <- 0.10 ; NGENE_CHR <- 5L
THRESHOLDS <- c(0.05, 0.10)          # 两档候选统一值

FLOOR <- as.integer(Sys.getenv("LUAD_COVERAGE_FLOOR", "840"))
cat(sprintf("[tiers] 覆盖度地板 = %s（0 表示不施加；判据 = M1 的 nFeature >= 地板）\n",
            if (FLOOR > 0L) paste0(FLOOR, " (nFeature)") else "无"))

res <- vector("list", length(samples))
for (i in seq_along(samples)) {
  s <- samples[i]
  M <- read_h5ad_sample(s)                       # 内部已做 .assert_provenance 硬断言
  n_cells_raw <- nrow(M)                         # h5ad 里该样本的细胞数（地板**之前**）

  nf <- m1_nfeature(M, s)                        # 与 h5ad 实测非零数当场对表
  n_cells_ge_floor <- sum(nf >= FLOOR)
  M <- apply_coverage_floor(M, nf, FLOOR)
  n_floor_dropped <- n_cells_raw - n_cells_ge_floor
  rm(nf); invisible(gc())

  ch <- copykat_chain(M, low_dr = LOW_DR, up_dr = UP_DR,
                      ngene_chr = NGENE_CHR, min_gene_per_cell = MIN_GENE_PER_CELL)
  rm(M); invisible(gc())

  # 守卫（可证伪）：地板 >= 200 时它严格强于 min.gene.per.cell 判据，
  # 故"进 copykat 的细胞数"必须**恰好等于**过地板的细胞数。若子集化写错（漏施加/施加两次/对齐错位），
  # 这条会先炸，而不是静默给出一份看起来正常的表。
  if (FLOOR >= MIN_GENE_PER_CELL && ch$n_cells_used != n_cells_ge_floor)
    stop(sprintf("[FAIL] %s: 过地板 %d 个细胞，但 copykat 链路收到 %d 个 —— 地板未正确施加。",
                 s, n_cells_ge_floor, ch$n_cells_used))

  # 逐档阈值的最终基因数：dr2 已是在 ToRemov2 存活细胞上算好的检出率，
  # 直接比阈值即可，不必重跑链路（dn 与阈值无关）。
  n_final <- vapply(THRESHOLDS, function(x) sum(ch$dr2 >= x), 0L)

  # 预测的 not.defined 数（**不跑 copykat 就能得到**）：被判据踢掉的细胞 = 不能分段的细胞。
  # 这是 GP2"退化样本处置"需要的量，P13 实测 801-494 = 307 ✔
  n_pred_nd <- ch$n_cells_used - length(ch$survivors)

  res[[i]] <- data.frame(
    sample_id = s,
    stage_token = sub("^[^_]*_", "", s),
    coverage_floor = FLOOR,
    n_cells_h5ad = n_cells_raw,
    n_cells_ge_floor = n_cells_ge_floor,
    n_cells_dropped_by_floor = n_floor_dropped,
    # 地板之后、再被 200 基因判据踢掉的数（地板=0 时此列退化为旧语义）
    n_cells_dropped_below_200genes = n_cells_ge_floor - ch$n_cells_used,
    n_cells_used = ch$n_cells_used,
    # 7000 判据的分母：LOW.DR 过滤后的**全基因**数（早于注释，copykat:57）
    n_after_LOWDR_fullgenes = ch$n_after_lowdr,
    # 注释 + 删周期/HLA 之后的判据基因集
    n_after_anno_cycle_hla = nrow(ch$anno),
    # **实际进入分段的基因集**（= copykat 的 anno.mat2），这是替代旧 genes_ge_* 的那个数
    n_genes_final = ch$n_genes_final,
    genes_ge_0.05 = n_final[[1]],
    genes_ge_0.10 = n_final[[2]],
    copykat_effective_UPDR = ch$up_dr_effective,
    n_pred_not_defined = n_pred_nd,
    rate_pred_not_defined = n_pred_nd / ch$n_cells_used,
    stringsAsFactors = FALSE
  )
  cat(sprintf("[%2d/75] %-12s h5ad=%6d -floor=%5d =%6d  LOWDRfull=%5d  最终=%5d  UPDR=%.2f  not.def=%5d(%.3f)\n",
              i, s, n_cells_raw, n_floor_dropped, ch$n_cells_used, ch$n_after_lowdr,
              ch$n_genes_final, ch$up_dr_effective, n_pred_nd, n_pred_nd / ch$n_cells_used))
}

df <- do.call(rbind, res)
suffix  <- if (FLOOR > 0L) sprintf("_floor%d", FLOOR) else ""
out <- sprintf("/home/eto/luad_v2/results/03_cnv/prereg_gene_tiers%s.csv", suffix)
write.csv(df, out, row.names = FALSE)
cat("\n写出:", out, "\n")

cat("\n=== 汇总 ===\n")
if (FLOOR > 0L) {
  cat(sprintf("地板 %d：h5ad 合计 %d → 过地板 %d（丢弃 %d，%.1f%%）；过地板后又被 200 基因判据踢掉 %d\n",
              FLOOR, sum(df$n_cells_h5ad), sum(df$n_cells_ge_floor),
              sum(df$n_cells_dropped_by_floor),
              100 * sum(df$n_cells_dropped_by_floor) / sum(df$n_cells_h5ad),
              sum(df$n_cells_dropped_below_200genes)))
  cat(sprintf("  过地板后细胞数 <200 的样本（→ not_testable）: %d\n",
              sum(df$n_cells_used < 200L)))
}
cat(sprintf("7000 判据分母 (n_after_LOWDR_fullgenes)：min=%d p25=%d 中位=%d max=%d\n",
            min(df$n_after_LOWDR_fullgenes), as.integer(quantile(df$n_after_LOWDR_fullgenes, .25)),
            as.integer(median(df$n_after_LOWDR_fullgenes)), max(df$n_after_LOWDR_fullgenes)))
for (col in c("genes_ge_0.05", "genes_ge_0.10", "n_genes_final")) {
  v <- as.integer(df[[col]])
  cat(sprintf("  %-16s min=%5d p25=%5d 中位=%5d max=%5d\n",
              col, min(v), as.integer(quantile(v, .25)), as.integer(median(v)), max(v)))
}
cat(sprintf("\n被覆写（UP.DR 静默改成 LOW.DR=0.05）的样本数 = %d / %d\n",
            sum(df$copykat_effective_UPDR == LOW_DR), nrow(df)))
cat("按分期：\n"); print(table(df$stage_token, df$copykat_effective_UPDR))
cat(sprintf("\n预测 not.defined 率：中位=%.3f  max=%.3f（样本 %s）\n",
            median(df$rate_pred_not_defined), max(df$rate_pred_not_defined),
            df$sample_id[which.max(df$rate_pred_not_defined)]))
