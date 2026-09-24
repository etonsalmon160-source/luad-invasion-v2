#!/usr/bin/env Rscript
# inferCNV 冒烟跑（P4 单患者）。
#
# ⚠️ 运行前提（必须，否则 rjags 装载失败）：
#     export LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH
#   跑法：
#     LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH \
#     R_LIBS=/home/eto/Rlibs/infercnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#     Rscript 03_cnv/16_run_infercnv_smoke.R
#
# 参数**写死在 INFERCNV_SMOKE_PREREG.md §五**（计算前登记），本脚本只执行，不改。
# 用法：`Rscript 16_run_infercnv_smoke.R [--object-only]`
#   --object-only 只建对象并打印维度，不跑管线（用于验证输入管道）。

suppressMessages({
  library(infercnv)
  library(Matrix)
  library(data.table)
})

ROOT <- "/home/eto/luad_v2"
D <- file.path(ROOT, "results/03_cnv/infercnv_smoke")
OUT <- file.path(D, "infercnv_out")

# ——— 预注册常量（INFERCNV_SMOKE_PREREG.md §五）———
P_CUTOFF <- 0.1
P_WINDOW <- 101
P_DENOISE <- TRUE
P_HMM <- FALSE
P_ANALYSIS_MODE <- "subclusters"
P_CLUSTER_BY_GROUPS <- TRUE
P_WRITE_EXPR <- TRUE
P_NO_PLOT <- FALSE
P_THREADS <- 8
P_REF_GROUP <- "Normal"

args <- commandArgs(trailingOnly = TRUE)
object_only <- "--object-only" %in% args

# inferCNV 自己警告：analysis_mode="subclusters" 下不加这行会在 hclust 时报错
options(scipen = 100)

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()

cat("[in] 读 COO ...\n")
coo <- fread(file.path(D, "p4_counts_coo.tsv.gz"), header = FALSE,
             col.names = c("i", "j", "x"), data.table = FALSE)
genes <- readLines(file.path(D, "p4_genes.txt"))
cells <- readLines(file.path(D, "p4_cells.txt"))
cat(sprintf("[in] COO %d 非零元；基因 %d；细胞 %d\n", nrow(coo), length(genes), length(cells)))

mat <- sparseMatrix(i = coo$i, j = coo$j, x = coo$x,
                    dims = c(length(genes), length(cells)),
                    dimnames = list(genes, cells))
rm(coo); invisible(gc())
cat(sprintf("[in] 矩阵 %d x %d，类 %s\n", nrow(mat), ncol(mat), class(mat)[1]))

cat("[obj] CreateInfercnvObject ...\n")
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = mat,
    gene_order_file   = file.path(D, "gene_order_hg38.tsv"),
    annotations_file  = file.path(D, "p4_annotations.tsv"),
    ref_group_names   = c(P_REF_GROUP),
    delim             = "\t")

# 🔴 必须转成 base matrix（数值不变，只换存储类型）。
# 以 dgCMatrix 为输入时 CreateInfercnvObject 会把 expr.data 存成 Matrix 包的
# `dgeMatrix`；而 inferCNV 的 .subtract_expr() 是按 **base matrix** 的语义逐行抽取的：
#   as.numeric(expr_matrix[row_idx, , drop=TRUE])
# 在 dgeMatrix 上单行抽取实测 **1.23 秒/行**，base matrix 上 **0.002 秒/行**（差 ~600x）。
# 7,668 行 × 1.23 s ≈ 157 分钟 —— 实跑中 STEP 08 已耗 22 分钟仍在原地，故停跑修正。
# `identical(as.matrix(dgeMatrix), base_matrix)` 为 TRUE：数值逐位相同，非口径变更。
infercnv_obj@expr.data <- as.matrix(infercnv_obj@expr.data)
stopifnot(is.matrix(infercnv_obj@expr.data), !isS4(infercnv_obj@expr.data))

cat(sprintf("[obj] 建成：%d 基因 x %d 细胞（expr.data 类型 %s）\n",
            nrow(infercnv_obj@expr.data), ncol(infercnv_obj@expr.data),
            paste(class(infercnv_obj@expr.data), collapse = "/")))
cat(sprintf("[obj] 基因区间（chr 取值）：%s\n",
            paste(head(sort(unique(as.character(infercnv_obj@gene_order$chr))), 30), collapse = ", ")))
cat(sprintf("[obj] 参考组：%s\n", paste(infercnv_obj@reference_grouped_cell_indices |> names(), collapse = ", ")))
cat(sprintf("[obj] 观测组：%s\n",
            paste(infercnv_obj@observation_grouped_cell_indices |> names(), collapse = ", ")))
cat(sprintf("[obj] 基因保留 %d / 18000（掉了 %d）\n",
            nrow(infercnv_obj@expr.data), 18000 - nrow(infercnv_obj@expr.data)))

if (object_only) {
  cat("[obj] --object-only，停在这里。\n")
  quit(save = "no", status = 0)
}

cat("[run] infercnv::run 开始 ...\n")
infercnv_obj <- infercnv::run(
    infercnv_obj,
    out_dir            = OUT,
    cutoff             = P_CUTOFF,
    window_length      = P_WINDOW,
    denoise            = P_DENOISE,
    HMM                = P_HMM,
    analysis_mode      = P_ANALYSIS_MODE,
    cluster_by_groups  = P_CLUSTER_BY_GROUPS,
    write_expr_matrix  = P_WRITE_EXPR,
    no_plot            = P_NO_PLOT,
    num_threads        = P_THREADS,
    plot_steps         = FALSE,
    save_rds           = TRUE,
    save_final_rds     = TRUE)

el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("[run] 完成，用时 %.1f 分钟\n", el))

# 落一份最终的残差矩阵（分析脚本读它）
saveRDS(infercnv_obj@expr.data, file.path(D, "p4_final_expr_data.rds"))
cat(sprintf("[out] %s（%d x %d）\n", file.path(D, "p4_final_expr_data.rds"),
            nrow(infercnv_obj@expr.data), ncol(infercnv_obj@expr.data)))
cat(sprintf("[done] 总用时 %.1f 分钟\n", el))
