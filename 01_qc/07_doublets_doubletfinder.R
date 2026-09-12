#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/07_doublets_doubletfinder.R — M1 第二方法交叉验证：DoubletFinder
#
# 动机：scrublet 在**稀疏 snRNA** 上失效（可检测比例 ~0.5%）→ 换 DoubletFinder
#      （基于 Seurat 的 kNN 人工双体法，对稀疏数据更稳）。
#
# 设定（与主流程一致）：
#   * 同一套 QC 掩码（读 gse308103_per_cell_qc.csv.gz 的 qc_pass）
#   * **固定 pK=0.09**（跳过 paramSweep：75 样本全跑 paramSweep 不可行；
#     pK=0.09 是 10x 数据的常用稳健值）
#   * nExp = round(rate × n)，rate=0.075（10x ~1 万细胞经验值）
#   * seed=1（确定性）
#
# 用法：
#   Rscript 01_qc/07_doublets_doubletfinder.R --n 3      # 试点 3 个样本
#   Rscript 01_qc/07_doublets_doubletfinder.R            # 全部 75
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat); library(DoubletFinder)
})
set.seed(1)

RAW <- "/home/eto/luad_invasion/data/GSE308103/extracted"
OUT <- "/home/eto/luad_v2/results/01_qc"
PK <- 0.09
RATE <- 0.075
NPCS <- 30

argv <- commandArgs(trailingOnly = TRUE)
N <- if ("--n" %in% argv) as.integer(argv[which(argv == "--n") + 1]) else Inf

qc <- fread(file.path(OUT, "gse308103_per_cell_qc.csv.gz"),
            select = c("sample_id", "cell_barcode", "qc_pass"))

read_barcode_line <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
}

files <- sort(list.files(RAW, pattern = "\\.raw_counts"), method = "radix")
if (is.finite(N)) files <- head(files, N)
cat(sprintf("[DF] %d 样本；pK=%.2f, nExp=%.3f*n, npcs=%d\n", length(files), PK, RATE, NPCS))

res <- list()
for (i in seq_along(files)) {
  f <- files[i]
  sid <- sub("^GSM[0-9]+_", "", sub("\\.raw_counts.*$", "", f))
  bc <- read_barcode_line(file.path(RAW, f))
  dt <- fread(cmd = sprintf("zcat %s/%s", RAW, f), header = FALSE, skip = 1)
  g  <- as.character(unlist(dt[[1]]))
  m  <- as.matrix(dt[, -1]); storage.mode(m) <- "double"
  rownames(m) <- make.unique(g); colnames(m) <- bc
  rm(dt); gc(verbose = FALSE)

  keep <- qc[sample_id == sid & qc_pass == TRUE, cell_barcode]
  idx <- which(paste0(bc, "|", sid) %in% keep)
  if (length(idx) < 100) { cat(sprintf("  [%s] 跳过(<100)\n", sid)); next }

  t0 <- Sys.time()
  seu <- CreateSeuratObject(counts = m[, idx, drop = FALSE])
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = NPCS, verbose = FALSE)

  nExp <- round(RATE * ncol(seu))
  seu <- doubletFinder(seu, PCs = 1:NPCS, pN = 0.25, pK = PK, nExp = nExp)
  cls_col <- grep("^DF\\.classifications", colnames(seu@meta.data), value = TRUE)[1]
  if (is.na(cls_col)) { cat(sprintf("  [%s] DF 未产出分类列\n", sid)); next }

  cls <- as.character(seu@meta.data[[cls_col]])
  res[[sid]] <- data.table(dataset = "GSE308103", sample_id = sid,
                           cell_barcode = paste0(colnames(seu), "|", sid),
                           df_class = ifelse(cls == "Doublet", "doublet", "singlet"))
  cat(sprintf("  [%2d/%d] %-12s n=%6d  DF 双体=%5d (%.1f%%)  用时 %.1f min\n",
              i, length(files), sid, ncol(seu), sum(cls == "Doublet"),
              100 * mean(cls == "Doublet"), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(seu, m); gc(verbose = FALSE)
}

if (!length(res)) { cat("[DF] 无结果\n"); quit(status = 1) }
df <- rbindlist(res)
suffix <- if (is.finite(N)) sprintf("_pilot%d", N) else ""
fwrite(df, file.path(OUT, sprintf("gse308103_doubletfinder_per_cell%s.csv.gz", suffix)))
cat(sprintf("\n[DF] 完成 %d 样本；总双体=%d (%.2f%%)\n", uniqueN(df$sample_id),
            sum(df$df_class == "doublet"), 100 * mean(df$df_class == "doublet")))
