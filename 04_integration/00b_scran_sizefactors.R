#!/usr/bin/env Rscript
# ============================================================================
# 🔴【已作废 · 2026-09-16】本脚本**请勿运行**，保留仅为审计留痕。
#   作废原因：用户决策「对齐论文全套」⇒ M3-A.1 归一化主口径改为
#     SCTransform(vst.flavor="v2")（v2 内部自带偏移/正则化），**不再单独算 size factor**。
#     论文配方里没有 scran 池化这一步。
#   登记位置：docs/PARAMETERS_AND_SOURCES.md §M3-A.1「已作废登记」表。
#   现行替代：04_integration/10_seurat_traditional.R
# ============================================================================
# 00b_scran_sizefactors.R —— scran 池化 size factor（**旧口径，已作废**）
#
# 出处：Lun, Bach & Marioni 2016, Genome Biol 17:75
#       "Pooling across cells to normalize single-cell RNA sequencing data
#        with many zero counts"
# 登记于 docs/PARAMETERS_AND_SOURCES.md M3-A.1 §「已作废登记」表（原为 🟡P 主口径，2026-09-16 作废）
#
# 输入：results/04_integration/scran_io/  （由 00a_export_counts_for_scran.py 产出）
# 输出：results/04_integration/scran_size_factors.csv.gz  （cell_barcode,sample_id,size_factor）
#       results/04_integration/scran_manifest.json
#
# ────────────────────────────────────────────────────────────────────────────
# 参数（全部显式钉死，不依赖函数默认值 —— 法则 3.1 / 法则 0）
#
#   quickCluster（scran 1.26.0，实际落到 quickClusterANY）
#     min.size      = 100         （包默认 100）
#     method        = "igraph"    （包默认 c("igraph","hclust") 取首项）
#     use.ranks     = FALSE       （包默认）
#     d             = NULL        （包默认；igraph 路径不做 PCA）
#     subset.row    = NULL        （包默认，内部按 min.mean 过滤）
#     min.mean      = 0.1         （包默认 NULL ⇒ 内部取 0.1）
#     graph.fun     = "walktrap"  （包默认）
#     block         = NULL        ← **显式选整份数据聚类，见下方"设计说明"**
#     block.BPPARAM = SerialParam()
#
#   computeSumFactors → scuttle::pooledSizeFactors（scran 1.26.0）
#     sizes            = seq(21, 101, 5)  （包默认，显式写回）
#     clusters         = 上面 quickCluster 的结果
#     ref.clust        = NULL
#     max.cluster.size = 3000     （包默认）
#     positive         = TRUE     （包默认）
#     scaling          = NULL     （包默认）
#     min.mean         = 0.1      （scran 规范值；包默认 NULL）
#     subset.row       = NULL
#
# ★ 设计说明（为什么 block=NULL）
#   scran 的聚类**只用于池化**，不是生物学注释，故不做 block。理由：
#   ① 整份聚类是 Lun 2016 原文的规范用法；
#   ② quickCluster 默认走 igraph/walktrap（图方法），非 hclust，不吃 O(n²) 内存，
#      648,945 细胞可整份跑；block 在此不是为绕开规模问题而加；
#   ③ max.cluster.size=3000 已把每个池子封顶，避免单个大簇主导；
#   ④ 逐样本的规模差异由下游批次分支（Arm A 不校正 / Arm B Harmony）负责处理，
#      归一化阶段先保持单一全局口径，便于 A/B 归因。
#   ⚠️ 副作用（如实记录，不外推）：池子可能跨样本混合。若下游 Arm A 观察到
#      明显的样本驱动簇，须回到此处把 block 改为 sample_id 做敏感性对照。
# ────────────────────────────────────────────────────────────────────────────
suppressMessages({
  library(Matrix); library(SingleCellExperiment); library(S4Vectors)
  library(scran); library(scuttle); library(BiocParallel)
})

ROOT <- "/home/eto/luad_v2"
IO   <- file.path(ROOT, "results/04_integration/scran_io")
OUT  <- file.path(ROOT, "results/04_integration")

T0 <- Sys.time()
el <- function() sprintf("[%6.1fs]", as.numeric(difftime(Sys.time(), T0, units="secs")))
log <- function(m) cat(el(), m, "\n", sep=" ")

N_WORKERS <- 18L
SIZES     <- seq(21, 101, 5)

# ---- 1. 读入三元组 --------------------------------------------------------
man_in <- jsonlite::fromJSON(file.path(IO, "export_manifest.json"))
nnz    <- man_in$source$nnz
n_cells <- man_in$source$n_cells
n_genes <- man_in$source$n_genes
log(sprintf("读入 CSC 三元组：%d 基因 × %d 细胞，nnz=%d", n_genes, n_cells, nnz))

data <- readBin(file.path(IO, "data_f32.bin"),    what="numeric", n=nnz,        size=4)
idx  <- readBin(file.path(IO, "indices_i32.bin"), what="integer", n=nnz,        size=4)
ptr  <- readBin(file.path(IO, "indptr_i32.bin"),  what="integer", n=n_cells+1L, size=4)
log(sprintf("  读毕；data max=%d  ptr 末值=%d", as.integer(max(data)), ptr[length(ptr)]))

stopifnot(ptr[length(ptr)] == nnz, length(idx) == nnz, length(data) == nnz)

genes <- readLines(file.path(IO, "gene_names.txt"))
cells <- readLines(file.path(IO, "cell_names.txt"))
stopifnot(length(genes) == n_genes, length(cells) == n_cells)

mat <- new("dgCMatrix", i = idx, p = ptr, x = data, Dim = c(as.integer(n_genes), as.integer(n_cells)))
rm(data, idx, ptr); invisible(gc())
rownames(mat) <- genes; colnames(mat) <- cells
log(sprintf("  dgCMatrix 就绪；总量 %.1f GiB", as.numeric(object.size(mat))/2^30))

meta <- read.csv(gzfile(file.path(IO, "cell_meta.csv.gz")), stringsAsFactors = FALSE)
stopifnot(identical(meta$cell_barcode, cells))     # 列序必须严格一致
log(sprintf("  cell_meta 列序一致；%d 样本 / %d 分期", length(unique(meta$sample_id)), length(unique(meta$stage))))

sce <- SingleCellExperiment(list(counts = mat),
        colData = DataFrame(sample_id = meta$sample_id, stage = meta$stage,
                            row.names = meta$cell_barcode))

# ---- 2. 池化聚类（仅用于池化）--------------------------------------------
log(sprintf("quickCluster method=igraph/graph.fun=walktrap min.size=100（%d 核）", N_WORKERS))
clu <- quickCluster(sce, assay.type="counts", min.size=100L,
                    method="igraph", use.ranks=FALSE, d=NULL, subset.row=NULL,
                    min.mean=0.1, graph.fun="walktrap", block=NULL,
                    BPPARAM=MulticoreParam(workers=N_WORKERS, progressbar=FALSE),
                    block.BPPARAM=SerialParam())
log(sprintf("  聚类数=%d；簇大小 min/median/max = %s", nlevels(clu),
            paste(summary(as.integer(table(clu)))[c("Min.","Median","Max.")], collapse=" / ")))

# ---- 3. 池化 size factor --------------------------------------------------
log("computeSumFactors（= scuttle::pooledSizeFactors）")
sce <- computeSumFactors(sce, assay.type="counts", sizes=SIZES, clusters=clu,
                         ref.clust=NULL, max.cluster.size=3000, positive=TRUE,
                         scaling=NULL, min.mean=0.1, subset.row=NULL,
                         BPPARAM=MulticoreParam(workers=N_WORKERS, progressbar=FALSE))
sf <- sizeFactors(sce)
stopifnot(length(sf) == n_cells, all(is.finite(sf)), all(sf > 0))

# ---- 4. 诊断：是否需要跨样本重标定 ---------------------------------------
q <- function(v) round(as.numeric(quantile(v, c(0,.25,.5,.75,1))), 4)
by_s <- split(sf, meta$sample_id)
med_by_s <- vapply(by_s, median, numeric(1))
log(sprintf("全局 sf 分位 min/p25/med/p75/max = %s", paste(q(sf), collapse=" / ")))
log(sprintf("几何均值=%.4f 算术均值=%.4f（scran 惯例中心化到 1）",
            exp(mean(log(sf))), mean(sf)))
log(sprintf("逐样本 sf 中位数：min=%.3f p25=%.3f med=%.3f p75=%.3f max=%.3f",
            q(med_by_s)[1], q(med_by_s)[2], q(med_by_s)[3], q(med_by_s)[4], q(med_by_s)[5]))
log(sprintf("逐样本中位数跨度倍数 = %.2f×", max(med_by_s)/min(med_by_s)))

# ---- 5. 落盘 --------------------------------------------------------------
out <- data.frame(cell_barcode = cells, sample_id = meta$sample_id,
                  stage = meta$stage, size_factor = sf, stringsAsFactors = FALSE)
csv <- file.path(OUT, "scran_size_factors.csv.gz")
con <- gzfile(csv, "w"); write.csv(out, con, row.names = FALSE, quote = FALSE); close(con)
log(sprintf("写出 %s（%.1f MiB）", basename(csv), file.info(csv)$size/2^20))

man <- list(
  script     = "04_integration/00b_scran_sizefactors.R",
  method     = "scran pooled size factors (Lun 2016, Genome Biol 17:75)",
  source     = man_in$source,
  params = list(
    quickCluster = list(min.size=100L, method="igraph", use.ranks=FALSE, d=NULL,
                        min.mean=0.1, graph.fun="walktrap", block=NULL),
    pooledSizeFactors = list(sizes=SIZES, ref.clust=NULL, max.cluster.size=3000L,
                             positive=TRUE, scaling=NULL, min.mean=0.1),
    n_workers = N_WORKERS, n_clusters = nlevels(clu)),
  diagnostics = list(
    sf_quantiles_global = as.list(q(sf)),
    sf_geometric_mean   = exp(mean(log(sf))), sf_arithmetic_mean = mean(sf),
    sf_median_by_sample = as.list(med_by_s),
    per_sample_median_span = max(med_by_s)/min(med_by_s),
    cluster_size_summary = as.list(summary(as.integer(table(clu))))),
  output = list(path="results/04_integration/scran_size_factors.csv.gz",
                n_cells=n_cells, md5=tools::md5sum(csv)[[1]]),
  versions = list(R=R.version.string, scran=as.character(packageVersion("scran")),
                  scuttle=as.character(packageVersion("scuttle")),
                  SingleCellExperiment=as.character(packageVersion("SingleCellExperiment")),
                  BiocParallel=as.character(packageVersion("BiocParallel")),
                  Matrix=as.character(packageVersion("Matrix"))),
  wall_sec = as.numeric(difftime(Sys.time(), T0, units="secs")))
jsonlite::write_json(man, file.path(OUT, "scran_manifest.json"), pretty=TRUE, auto_unbox=TRUE)
log("完成")
