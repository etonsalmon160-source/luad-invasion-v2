#!/usr/bin/env Rscript
# 08_depth_guard.R —— M6 §5.1 深度配平守卫（§5.1.1 已签口径，本脚本不新开任何口径）
#
# 要回答的问题：**这些"生态位域"是不是只是测序深度差**？
#
# 已签口径（NICHE_PREREG.md §5.1.1，2026-10-01 用户逐条裁定）：
#   深-1 五档按谁切  = **逐切片切五分位**：每张切片把自己的 mask-pass spot 按深度切 5 等份
#                      ⇒ 每张都进每一档，档间差异**只能是深度**，不是切片／期别。
#   深-2 共同深度    = **降到「全队列最浅那张切片的深度中位」**，用**二项稀释**把过深 spot 的
#                      counts 降下来（**保 spot、不剔 spot**；比目标浅的 spot 不动，做不到"提升"）。
#   深-3 深度量      = 沿用 §20：**nUMI 为主**，nFeature 并列复核，两套都报。
#   执行范围         = 每次**只跑 §3.4 选中的那档 `(λ*, k_geom*, r*)`**；两版 AGF 都做。
#
# 归并法（§5.1.1 已写死，**非口径变更**）：每档重跑得到该档自己的域 → 按**域表达画像**与
#   **全局共识类型**画像的相关，取最大者归入 ⇒ 再统计各全局类型在该档的 spot 占比。
#   **不**把全局类型硬贴到 spot 上。
#
# 🔴 只上报，不改任何口径；本臂不给任何恶性标签（§7 禁令）。
#
# 跑法（一次一张切片 × 一版 AGF；条件在进程内循环、逐个原子落盘、断点可续）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 10_niche/08_depth_guard.R --slide GSM9226168_P1_AAH --agf TRUE

suppressMessages({
  library(Banksy); library(Matrix); library(Seurat); library(data.table); library(dbscan)
})

## —————————————————————————————————————————————————————————————
## 一、参数
## —————————————————————————————————————————————————————————————
.a <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, default = NULL) {
  i <- match(k, .a)
  if (is.na(i)) return(default)
  if (i == length(.a)) stop("参数 ", k, " 后面缺值", call. = FALSE)
  .a[i + 1]
}
SLIDE <- get_arg("--slide")
AGF_C <- get_arg("--agf", "TRUE")
METRIC<- get_arg("--metric", "nUMI")        # 深-3：nUMI 为主
if (is.null(SLIDE) || is.na(SLIDE)) stop("必须给 --slide <GSM..._P..._...>", call. = FALSE)
if (!toupper(AGF_C) %in% c("TRUE", "FALSE")) stop("--agf 只能是 TRUE 或 FALSE", call. = FALSE)
if (!METRIC %in% c("nUMI", "nFeature")) stop("--metric 只能是 nUMI 或 nFeature", call. = FALSE)
USE_AGF <- toupper(AGF_C) == "TRUE"
AGF_TAG <- if (USE_AGF) "agfT" else "agfF"

ROOT   <- "/home/eto/luad_v2"
VISIUM <- file.path(ROOT, "data/visium_spatial")
RES    <- file.path(ROOT, "results/08_spatial_deconv")
NICHE  <- file.path(ROOT, "results/10_niche")
CONS   <- file.path(NICHE, "consensus", AGF_TAG)
OUTD   <- file.path(NICHE, "depth_guard", AGF_TAG)
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

NRM_FACTOR  <- 1e4           # [已签 S2]
NPCS        <- 20L           # [已签]
K_NEIGHBORS <- 50L           # [已签]
SPATIAL_MODE<- "kNN_median"  # [已签]
MIN_FREE_GB <- 25
THIN_SEED   <- 20261002L     # 🔴 稀释是随机操作 ⇒ 种子必须落盘
N_QUINT     <- 5L

## —————————————————————————————————————————————————————————————
## 二、工具
## —————————————————————————————————————————————————————————————
t_start <- Sys.time()
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
free_gb <- function() {
  x <- tryCatch(readLines("/proc/meminfo"), error = function(e) NULL); if (is.null(x)) return(Inf)
  v <- x[grepl("^MemAvailable:", x)]; if (!length(v)) return(Inf)
  as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1", v[1])) / 1048576
}
guard_mem <- function(where) {
  f <- free_gb()
  if (f < MIN_FREE_GB) stop(sprintf("%s：可用内存 %.1f GB < 守卫线 %.0f GB ⇒ 硬停", where, f, MIN_FREE_GB), call. = FALSE)
  f
}
read_gz <- function(p, fn) { con <- gzfile(p, "rt"); on.exit(close(con)); fn(con) }
ir <- function(x) as.integer(round(x))
atomic_write <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("原子改名失败：", tmp, " -> ", path)
}
atomic_write_lines <- function(txt, path) {
  tmp <- paste0(path, ".tmp"); writeLines(txt, tmp)
  if (!file.rename(tmp, path)) stop("原子改名失败：", tmp, " -> ", path)
}

## —————————————————————————————————————————————————————————————
## 三、读已签的选中配置（🔴 从 consensus_domains.tsv 直接读，不从 selection.tsv 的第一行猜：
##     selection.tsv 是按 ari_med 排序的，第一行不一定是选中档）
## —————————————————————————————————————————————————————————————
cdf <- file.path(CONS, "consensus_domains.tsv")
if (!file.exists(cdf)) stop("缺 ", cdf, "（先跑 03 + 04）", call. = FALSE)
KD <- fread(cdf)
cfg <- unique(KD[, .(k_geom, lambda, resolution)])
if (nrow(cfg) != 1) stop("consensus_domains.tsv 里出现了多组配置 ⇒ 该文件不是「选中一档」的产物，硬停", call. = FALSE)
KG_SEL <- as.integer(cfg$k_geom); LAM_SEL <- as.numeric(cfg$lambda); RES_SEL <- as.numeric(cfg$resolution)

## 全局共识类型画像（§3.5 产物；一行一个共识类型）
pf_f <- file.path(CONS, "consensus_profiles.tsv")
if (!file.exists(pf_f)) stop("缺 ", pf_f, call. = FALSE)
PF <- fread(pf_f)
if (!"consensus_type" %in% names(PF)) stop("consensus_profiles.tsv 缺 consensus_type 列", call. = FALSE)
PF_types <- as.integer(PF$consensus_type)
PF_genes <- setdiff(names(PF), "consensus_type")
PFM <- as.matrix(PF[, ..PF_genes]); rownames(PFM) <- as.character(PF_types)

step("==== §5.1 深度守卫：%s（use_agf=%s，metric=%s）====", SLIDE, USE_AGF, METRIC)
step("选中档：k_geom=%d lambda=%.2f resolution=%.2f；全局共识类型 %d 个 / 画像基因 %d 个",
     KG_SEL, LAM_SEL, RES_SEL, nrow(PFM), ncol(PFM))

## —————————————————————————————————————————————————————————————
## 四、读片（与 02_banksy_grid.R **同源**，保证与主网格可比）
## —————————————————————————————————————————————————————————————
HVG <- readLines(file.path(NICHE, "hvg_3000.txt"))

sm_all <- read.delim(gzfile(file.path(RES, "spot_mask.tsv.gz")), stringsAsFactors = FALSE)
sm_all <- sm_all[as.character(sm_all$pass) %in% c("TRUE", "true", "1"), ]
bc_keep <- sm_all$barcode[sm_all$slide == SLIDE]
if (!length(bc_keep)) stop("掩码通过集里没有这张切片：", SLIDE, call. = FALSE)

d  <- file.path(VISIUM, SLIDE, "filtered_feature_bc_matrix")
ft <- read_gz(file.path(d, "features.tsv.gz"), function(con)
  read.delim(con, header = FALSE, col.names = c("id", "sym", "type"), stringsAsFactors = FALSE))
ab <- which(ft$type == "Antibody Capture")
if (length(ab)) { ft <- ft[-ab, , drop = FALSE]; step("剔抗体行 %d（剩 %d 行）", length(ab), nrow(ft)) }

po <- read.csv(file.path(VISIUM, SLIDE, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
po <- po[po$in_tissue == 1L, , drop = FALSE]
i  <- match(bc_keep, po$barcode)
if (anyNA(i)) stop(SLIDE, "：mask-pass 的 barcode 不在 tissue_positions.csv 里", call. = FALSE)
po <- po[i, , drop = FALSE]; stopifnot(identical(po$barcode, bc_keep))
xy <- cbind(x = as.numeric(po$pxl_col_in_fullres), y = as.numeric(po$pxl_row_in_fullres))
xy <- sweep(xy, 2, apply(xy, 2, min))    # 平移非负（距离不变）
rownames(xy) <- bc_keep

bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
m  <- readMM(file.path(d, "matrix.mtx.gz"))
stopifnot(nrow(m) == nrow(ft) + length(ab), ncol(m) == length(bc))
if (length(ab)) m <- m[-ab, , drop = FALSE]
colnames(m) <- bc
j <- match(bc_keep, bc)
if (anyNA(j)) stop(SLIDE, "：mask-pass barcode 不在表达矩阵里", call. = FALSE)
m <- m[, j, drop = FALSE]
## 同 symbol 合并（与 02 逐字同源）
g <- factor(ft$sym, levels = unique(ft$sym))
A <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1, dims = c(nlevels(g), length(g)))
m <- A %*% m
rownames(m) <- levels(g); colnames(m) <- bc_keep
m <- as(m, "CsparseMatrix")

DEPTH_UMI <- Matrix::colSums(m)                       # nUMI（全基因，与掩码表 total_umi 同口径）
DEPTH_GEN <- Matrix::colSums(m > 0)                   # nFeature
DEPTH <- if (METRIC == "nUMI") DEPTH_UMI else DEPTH_GEN
step("读入 %d 基因 x %d spot；%s 中位 %s、范围 %s–%s（可用 %.1f GB）",
     nrow(m), ncol(m), METRIC, format(ir(median(DEPTH)), big.mark = ","),
     format(ir(min(DEPTH)), big.mark = ","), format(ir(max(DEPTH)), big.mark = ","), free_gb())

## 全片空间近邻距离（背景值，给"档内 spot 是否被打散"这个诊断做对照）
NN_FULL_MED <- median(dbscan::kNN(xy, k = 1L)$dist)

## 从 raw counts 取 HVG 子集并做 CP10K+log1p（S2 同口径）
prep_expr <- function(M) {
  missed <- setdiff(HVG, rownames(M))
  if (length(missed)) stop(SLIDE, "：HVG 轴缺 ", length(missed), " 个基因", call. = FALSE)
  E <- as(M[HVG, , drop = FALSE], "CsparseMatrix")
  N <- E
  N@x <- log1p(N@x / rep(pmax(Matrix::colSums(E), 1e-9), diff(N@p)) * NRM_FACTOR)
  list(E = E, N = N)
}

## —————————————————————————————————————————————————————————————
## 五、condition 定义
## —————————————————————————————————————————————————————————————
## 深-2 目标：全队列最浅那张切片的深度中位（从 06 的覆盖表读，口径 = mask-pass spot 的中位 total_umi）
tgt_f <- file.path(NICHE, "guards", "coverage_per_slide.tsv")
if (!file.exists(tgt_f)) stop("缺 ", tgt_f, "（先跑 06_coverage_guard.R）", call. = FALSE)
COV <- fread(tgt_f)
TARGET_UMI <- min(COV$umi_med_pass, na.rm = TRUE)
step("深-2 共同深度目标 = %.1f UMI（全队列最浅：%s）",
     TARGET_UMI, COV$slide[which.min(COV$umi_med_pass)])

## 🔴 用**秩**切五档，不用 `quantile` + `cut`：深度是离散整数，等值点很多，
##    按数值切会出现**空档**（`cut` 的 break 落在同一个值上），空档直接让 08 硬停。
##    秩法保证每档约 n/5 个 spot 且非空；并列值（平均秩相同）整块落进同一档。
rk   <- rank(DEPTH, ties.method = "average")
qidx <- pmin(N_QUINT, ceiling(rk / (length(rk) + 1e-9) * N_QUINT))
q_lo <- vapply(seq_len(N_QUINT), function(q) min(DEPTH[qidx == q]), numeric(1))
q_hi <- vapply(seq_len(N_QUINT), function(q) max(DEPTH[qidx == q]), numeric(1))
CONDS <- c(sprintf("q%d", seq_len(N_QUINT)), "down", "full")
step("五档（按 %s 的秩逐切片等分）：%s", METRIC,
     paste(sprintf("q%d=[%d,%d]", seq_len(N_QUINT), ir(q_lo), ir(q_hi)), collapse = "  "))

## 二项稀释（保 spot）：过深的 spot 按 p = T/total 抽稀；比目标浅的不动
thin_to_target <- function(M, target, seed) {
  tot <- Matrix::colSums(M)
  p <- pmin(1, target / pmax(tot, 1e-9))
  set.seed(seed)
  ## 🔴 rbinom 返回 integer，而 dgCMatrix 的 `x` 槽必须是 double ⇒ 不 as.numeric 会报
  ##    `invalid class "dgCMatrix" object: 'x' slot is not of type "double"`。
  M@x <- as.numeric(rbinom(length(M@x), size = as.integer(round(M@x)), prob = rep(p, diff(M@p))))
  M
}

## 域 → 全局共识类型：按表达画像相关取最大
map_to_global <- function(N_sub, lab) {
  lv <- sort(unique(lab))
  G <- intersect(rownames(N_sub), PF_genes)
  if (length(G) < 100) stop("与全局画像的共有基因只有 ", length(G), " 个 ⇒ 硬停", call. = FALSE)
  Pq <- vapply(lv, function(k) Matrix::rowMeans(N_sub[G, lab == k, drop = FALSE]), numeric(length(G)))
  if (is.null(dim(Pq))) Pq <- matrix(Pq, ncol = length(lv))
  rownames(Pq) <- G; colnames(Pq) <- as.character(lv)
  R <- cor(Pq, t(PFM[, G, drop = FALSE]), method = "pearson")   # 域 x 全局类型
  if (is.null(dim(R))) R <- matrix(R, nrow = length(lv))
  best <- max.col(R, ties.method = "first")
  ord <- apply(R, 1, function(v) order(v, decreasing = TRUE))
  second <- vapply(seq_len(nrow(R)), function(i) ord[2, i], integer(1))
  ## 🔴 `cor_margin` 必须上报：共识类型 1 独大（308/1005 个域），argmax 会让大量域都倒向它，
  ##    而 best/second 常常只差 0.02–0.03 ⇒ 归并本身可能是**平局噪音**。留数不设闸，但必须能看见。
  data.table(domain = as.integer(lv), best_type = PF_types[best],
             best_cor = R[cbind(seq_len(nrow(R)), best)],
             second_type = PF_types[second],
             second_cor = R[cbind(seq_len(nrow(R)), second)],
             cor_margin = R[cbind(seq_len(nrow(R)), best)] - R[cbind(seq_len(nrow(R)), second)],
             n_spot_dom = as.integer(table(lab)[as.character(lv)]), genes_matched = length(G))
}

## `full` 档的自检：08 复现出来的域分配，必须与主网格 banksy/<slide>/domains_<agf>.tsv.gz 同档同种子一致。
## 🔴 不一致就说明 08 与主臂**不是同一套机器**，那么"同机制归并"的前提不成立 ⇒ 必须登记。
check_reproduce_main <- function(lab_cells, lab) {
  bf <- file.path(NICHE, "banksy", SLIDE, sprintf("domains_%s.tsv.gz", AGF_TAG))
  if (!file.exists(bf)) return(list(ari = NA_real_, note = "缺主网格域文件"))
  D <- fread(bf)
  D <- D[k_geom == KG_SEL & lambda == LAM_SEL & resolution == RES_SEL & seed == 0]
  j <- match(names(lab_cells), D$barcode)
  if (anyNA(j)) return(list(ari = NA_real_, note = "barcode 对不上"))
  list(ari = mclust::adjustedRandIndex(lab, D$domain[j]), note = "ok")
}

run_cond <- function(cond) {
  outf <- file.path(OUTD, sprintf("%s__%s.tsv", SLIDE, cond))
  if (file.exists(outf)) { step("  断点跳过 %s（已有 %s）", cond, basename(outf)); return(invisible(NULL)) }
  guard_mem(sprintf("%s 起跑前", cond))
  t0 <- Sys.time()

  if (cond == "down") {
    Msub <- thin_to_target(m, TARGET_UMI, THIN_SEED)
    idx  <- seq_len(ncol(Msub)); xys <- xy
    lo <- hi <- NA_real_
  } else if (cond == "full") {
    Msub <- m; idx <- seq_len(ncol(m)); xys <- xy
    lo <- hi <- NA_real_
  } else {
    q <- as.integer(sub("^q", "", cond))
    idx <- which(qidx == q); if (!length(idx)) stop(cond, " 没有 spot", call. = FALSE)
    Msub <- m[, idx, drop = FALSE]; xys <- xy[idx, , drop = FALSE]
    lo <- q_lo[q]; hi <- q_hi[q]
  }
  n <- length(idx)
  if (n < 30L) stop(cond, " 只有 ", n, " 个 spot，太小 ⇒ 硬停", call. = FALSE)
  pe <- prep_expr(Msub)
  nnq <- median(dbscan::kNN(xys, k = 1L)$dist)
  dmed <- median(if (cond == "down") Matrix::colSums(Msub) else DEPTH[idx])

  bank <- BanksyObject(own.expr = pe$E,
                       cell.locs = data.frame(sdimx = xys[, 1], sdimy = xys[, 2],
                                              row.names = colnames(pe$E)),
                       meta.data = data.frame(row.names = colnames(pe$E)))
  bank <- ComputeBanksy(bank, compute_agf = USE_AGF, spatial_mode = SPATIAL_MODE,
                        k_geom = KG_SEL, verbose = FALSE)
  bank <- NormalizeBanksy(bank, assay = "both", norm_factor = NRM_FACTOR,
                          log_norm = TRUE, pseudocount = 1, base = 10)
  bank <- ScaleBanksy(bank, assay = "both")
  ## 🔴 必须先 RunBanksyPCA 再 ClusterBanksy：直接 ClusterBanksy(pca=TRUE) 会报
  ##    `Run PCA with use_agf=NA lambda=1`（ClusterBanksy 的 checkArgs 要读已算好的 PCA）。
  ##    02_banksy_grid.R 也是这样两步走的（第 226/230 行）⇒ 与主网格同源。
  bank <- RunBanksyPCA(bank, lambda = LAM_SEL, use_agf = USE_AGF, npcs = NPCS, verbose = FALSE)
  CN_TAG <- if (USE_AGF) "M1" else "M0"
  bank <- ClusterBanksy(bank, lambda = LAM_SEL, use_agf = USE_AGF, pca = TRUE, npcs = NPCS,
                        method = "leiden", k.neighbors = K_NEIGHBORS, resolution = RES_SEL,
                        num.cores = 1, seed = 0L, verbose = FALSE)
  cn <- paste0("clust_", CN_TAG, "_lam", LAM_SEL, "_k", K_NEIGHBORS, "_res", RES_SEL)
  if (!cn %in% colnames(bank@meta.data))
    stop("找不到聚类列 ", cn, "；实际：",
         paste(grep("^clust_", colnames(bank@meta.data), value = TRUE), collapse = " "), call. = FALSE)
  lab <- as.integer(bank@meta.data[[cn]])
  names(lab) <- rownames(bank@meta.data)

  if (cond == "full") {
    chk <- check_reproduce_main(lab, lab)
    REPRO[[1]] <<- chk
    step("  full 档自检：与主网格同档同种子的域分配 ARI = %s（%s）",
         if (is.na(chk$ari)) "NA" else sprintf("%.4f", chk$ari), chk$note)
  }

  mp <- map_to_global(pe$N, lab)
  mp[, `:=`(slide = SLIDE, agf = AGF_TAG, cond = cond, metric = METRIC,
            depth_lo = lo, depth_hi = hi, depth_med = dmed,
            n_spot = n, n_domain = length(unique(lab)),
            nn_med_um = nnq * 0.2304, nn_full_med_um = NN_FULL_MED * 0.2304)]
  mp[, frac_dom := n_spot_dom / n]
  setcolorder(mp, c("slide", "agf", "cond", "metric", "n_spot", "n_domain",
                    "depth_lo", "depth_hi", "depth_med", "nn_med_um", "nn_full_med_um",
                    "domain", "n_spot_dom", "frac_dom",
                    "best_type", "best_cor", "second_type", "second_cor", "genes_matched"))
  setorder(mp, domain)
  atomic_write(mp, outf)
  step("  ✅ %s：%d spot / %d 域 / 中位深度 %s / 档内近邻 %.0fµm（全片 %.0fµm）→ %.1f min",
       cond, n, mp$n_domain[1], format(ir(dmed), big.mark = ","),
       mp$nn_med_um[1], mp$nn_full_med_um[1], elapsed(t0))
  invisible(NULL)
}

REPRO <- list(NULL)
for (cond in CONDS) run_cond(cond)

man <- list(script = "10_niche/08_depth_guard.R", prereg = "NICHE_PREREG.md §5.1/§5.1.1（已签）",
            slide = SLIDE, agf = AGF_TAG, metric = METRIC, use_agf = USE_AGF,
            k_geom = KG_SEL, lambda = LAM_SEL, resolution = RES_SEL, seed = 0L,
            n_quintile = N_QUINT, quintile_lo = as.numeric(q_lo), quintile_hi = as.numeric(q_hi),
            target_umi_downsamp = TARGET_UMI, thin_seed = THIN_SEED,
            nn_full_med_um = NN_FULL_MED * 0.2304,
            repro_vs_main_ari = if (is.null(REPRO[[1]])) NA_real_ else REPRO[[1]]$ari,
            repro_vs_main_note = if (is.null(REPRO[[1]])) "未跑" else REPRO[[1]]$note,
            note = "只上报，不改口径；不含任何生物学结论；不给恶性标签（§7 禁令）",
            finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
atomic_write_lines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE, digits = NA),
                   file.path(OUTD, sprintf("%s__manifest.json", SLIDE)))
step("==== %s %s %s 结束，用时 %.1f 分钟（可用 %.1f GB）====", SLIDE, AGF_TAG, METRIC,
     elapsed(t_start), free_gb())
