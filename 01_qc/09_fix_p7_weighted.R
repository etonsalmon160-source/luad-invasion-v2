#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/09_fix_p7_weighted.R — P7_LUAD 双体判定的补正（xgb 塌缩 → weighted）
#
# 背景（已实测）：
#   * `P7_LUAD` 用 scDblFinder 默认 `score="xgb"` 时**分类器塌缩**（iter≥1 后排除细胞数 1、2），
#     阈值 → 0.999、全部分数 = 0、判 0 双体；
#   * 该样本**质量正常**（median nCount 1,478 / nFeature 1,084，高于队列中位）；
#   * 同患者另两个样本（P7_Normal 12.05%、P7_LUAD1 12.37%）均正常 → **非生物学原因**；
#   * 改用 `score="weighted"` → 阈值 0.325、**1,158 双体（11.7%）**，与其他样本一致。
#
# 本脚本：对该**单个样本**用 `score="weighted"` 重跑，并按 `cell_barcode` **回填**主表。
#         其余 74 个样本不变（保持 `score="xgb"` 默认）。
# 确定性：set.seed(1)；SerialParam()。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(SingleCellExperiment)
  library(scDblFinder); library(scuttle); library(BiocParallel)
})
set.seed(1); register(SerialParam())

RAW <- "/home/eto/luad_invasion/data/GSE308103/extracted"
OUT <- "/home/eto/luad_v2/results/01_qc"
SID <- "P7_LUAD"; SCORE <- "weighted"; NMADS <- 3; MAX_MT <- 5

f <- list.files(RAW, pattern = paste0(SID, "\\.raw_counts"), full.names = FALSE)
stopifnot(length(f) == 1)
con <- gzfile(file.path(RAW, f), "rt"); bc <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]; close(con)
dt <- fread(cmd = sprintf("zcat %s/%s", RAW, f), header = FALSE, skip = 1)
g <- as.character(unlist(dt[[1]])); m <- as.matrix(dt[, -1]); storage.mode(m) <- "double"
stopifnot(length(bc) == ncol(m), length(g) == nrow(m))

cs <- colSums(m); fs <- colSums(m > 0); mt <- grepl("^MT-", g)
pct <- 100 * colSums(m[mt, , drop = FALSE]) / pmax(cs, 1)
pass <- !isOutlier(cs, nmads = NMADS, type = "both", log = TRUE) &
        !isOutlier(fs, nmads = NMADS, type = "both", log = TRUE) & (pct < MAX_MT)

idx <- which(pass)
sce <- SingleCellExperiment(list(counts = m[, idx, drop = FALSE]))
sce <- scDblFinder(sce, dbr = NULL, dbr.sd = NULL, aggregateFeatures = FALSE,
                   score = SCORE, BPPARAM = SerialParam())
cls <- rep("not_tested", ncol(m)); sco <- rep(NA_real_, ncol(m))
cls[idx] <- as.character(sce$scDblFinder.class); sco[idx] <- sce$scDblFinder.score
cat(sprintf("[fix] %s score=%s: pass=%d 双体=%d (%.1f%%)\n", SID, SCORE, sum(pass),
            sum(cls == "doublet"), 100 * mean(cls == "doublet")))

fix <- data.table(cell_barcode = paste0(bc, "|", SID),
                  doublet_class = cls, doublet_score = round(sco, 4))

pc <- fread(file.path(OUT, "gse308103_per_cell_qc.csv.gz"))
stopifnot(!any(duplicated(fix$cell_barcode)))
# data.table update join：按 cell_barcode 精确替换该样本的判定
pc[fix, on = "cell_barcode", `:=`(doublet_class = i.doublet_class,
                                  doublet_score = i.doublet_score)]
fwrite(pc, file.path(OUT, "gse308103_per_cell_qc.csv.gz"))

ps <- pc[, .(n_pre = .N, n_pass = sum(qc_pass),
  n_outlier_counts = sum(outlier_counts & !outlier_features),
  n_outlier_features = sum(outlier_features & !outlier_counts),
  n_mt_filtered = sum(pct_mt >= MAX_MT & !outlier_counts & !outlier_features),
  n_doublet = sum(doublet_class == "doublet"),
  doublet_rate = round(100 * sum(doublet_class == "doublet") / max(sum(qc_pass), 1), 2)),
  by = .(dataset, sample_id)][order(sample_id)]
ps[, doublet_note := fifelse(sample_id == SID, "fixed:score=weighted(xgb 塌缩)", "")]  # 覆盖旧标记
fwrite(ps, file.path(OUT, "gse308103_qc_per_sample.csv"))

cat(sprintf("[fix] 回填完成。%s 双体率=%.2f%%  全队列 pre=%d pass=%d doublet=%d (%.2f%%)\n",
            SID, ps[sample_id == SID, doublet_rate], nrow(pc), sum(pc$qc_pass),
            sum(pc$doublet_class == "doublet"),
            100 * sum(pc$doublet_class == "doublet") / sum(pc$qc_pass)))
