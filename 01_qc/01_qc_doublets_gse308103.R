#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/01_qc_doublets_gse308103.R — M1: GSE308103 (snRNA) QC + scDblFinder 逐样本
#
# 阈值依据（先测后定，见 logs/M1_308103_metrics.log）：
#   * 本数据为 **snRNA（细胞核）**，nCount 中位数 ~1,516、mt% 中位数 0.6 →
#     **不适用** scRNA 的固定阈值（nCount>=1000 会砍掉 ~30%）。
#   * nCount / nFeature：**逐样本 MAD 离群**（scuttle::isOutlier, nmads=3, 双尾, log1p）
#     —— 标准做法（McCarthy 2017 / scater-scuttle）。
#   * pct_mt < 5：核 mt 本就低（全队列 99% 分位 = 6.44），固定低阈值；
#     出处：本项目约定（标注于报告），配合 MAD 使用。
#   * 双体：scDblFinder（dbr=NULL 按细胞数自估；dbr.sd=NULL；aggregateFeatures=FALSE）
#     —— 逐样本，禁止合池。
# 确定性：set.seed(1)；SerialParam()。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(SingleCellExperiment)
  library(scDblFinder); library(scuttle); library(BiocParallel)
})
set.seed(1); register(SerialParam())

RAW <- "/home/eto/luad_invasion/data/GSE308103/extracted"
OUT <- "/home/eto/luad_v2/results/01_qc"

NMADS   <- 3
MAX_MT  <- 5
DBL_DBR <- NULL; DBL_DBR_SD <- NULL; DBL_AGGREGATE <- FALSE

files <- sort(list.files(RAW, pattern = "\\.raw_counts"), method = "radix")
cat(sprintf("[M1/308103] %d 样本；nmads=%d, mt<%d, scDblFinder(dbr=NULL)\n",
            length(files), NMADS, MAX_MT))

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

  out_c  <- isOutlier(cs, nmads = NMADS, type = "both", log = TRUE)
  out_f  <- isOutlier(fs, nmads = NMADS, type = "both", log = TRUE)
  pass   <- !out_c & !out_f & (pct < MAX_MT)

  cls <- rep("not_tested", ncol(m)); sco <- rep(NA_real_, ncol(m))
  idx <- which(pass)
  if (length(idx) >= 50) {
    sce <- SingleCellExperiment(list(counts = m[, idx, drop = FALSE]))
    sce <- scDblFinder(sce, dbr = DBL_DBR, dbr.sd = DBL_DBR_SD,
                       aggregateFeatures = DBL_AGGREGATE, BPPARAM = SerialParam())
    # Fallback：xgb 分类器偶发塌缩（阈值→1、判 0 双体）。实测 P7_LUAD 属此情形，
    # 但用 score="weighted" 正常（11.7%）。此处自动回退并标注，避免静默产出 0。
    if (sum(sce$scDblFinder.class == "doublet") == 0) {
      cat(sprintf("    [fallback] %s: xgb 判 0 双体 → 改用 score=weighted\n", sid))
      sce <- scDblFinder(sce, dbr = DBL_DBR, dbr.sd = DBL_DBR_SD,
                         aggregateFeatures = DBL_AGGREGATE, score = "weighted",
                         BPPARAM = SerialParam())
    }
    cls[idx] <- as.character(sce$scDblFinder.class)
    sco[idx] <- sce$scDblFinder.score
  }
  res[[i]] <- data.table(
    dataset = "GSE308103", sample_id = sid,
    cell_barcode = paste0(bc, "|", sid),
    nCount = cs, nFeature = fs, pct_mt = round(pct, 3),
    outlier_counts = out_c, outlier_features = out_f,
    qc_pass = pass, doublet_class = cls, doublet_score = round(sco, 4))
  cat(sprintf("  [%2d/%d] %-12s n=%6d pass=%6d outC=%5d outF=%5d mt=%4d | doublet=%5d (%.1f%%)\n",
              i, length(files), sid, ncol(m), sum(pass), sum(out_c & !out_f), sum(out_f & !out_c),
              sum(pct >= MAX_MT), sum(cls == "doublet", na.rm = TRUE),
              100 * sum(cls == "doublet", na.rm = TRUE) / max(sum(pass), 1))); flush.console()
  rm(dt, m); gc(verbose = FALSE)
}

per_cell <- rbindlist(res)
fwrite(per_cell, file.path(OUT, "gse308103_per_cell_qc.csv.gz"))

per_sample <- per_cell[, .(
  n_pre = .N, n_pass = sum(qc_pass),
  n_outlier_counts = sum(outlier_counts & !outlier_features),
  n_outlier_features = sum(outlier_features & !outlier_counts),
  n_mt_filtered = sum(pct_mt >= MAX_MT & !outlier_counts & !outlier_features),
  n_doublet = sum(doublet_class == "doublet"),
  doublet_rate = round(100 * sum(doublet_class == "doublet") / max(sum(qc_pass), 1), 2)
), by = .(dataset, sample_id)][order(sample_id)]
fwrite(per_sample, file.path(OUT, "gse308103_qc_per_sample.csv"))

cat("\n=== M1/308103 汇总 ===\n")
cat(sprintf("pre=%d  pass=%d (%.1f%%)  doublet=%d (%.2f%% of pass)\n",
            nrow(per_cell), sum(per_cell$qc_pass), 100*mean(per_cell$qc_pass),
            sum(per_cell$doublet_class == "doublet", na.rm = TRUE),
            100*sum(per_cell$doublet_class == "doublet", na.rm=TRUE)/max(sum(per_cell$qc_pass),1)))
cat(sprintf("剔除：低质量(离群) = %d\n", sum(!per_cell$qc_pass)))
cat(sprintf("产物: %s/gse308103_per_cell_qc.csv.gz\n", OUT))
