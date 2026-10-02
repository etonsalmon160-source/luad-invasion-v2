#!/usr/bin/env Rscript
# 02_banksy_grid.R —— M6 生态位臂 §3.2/§3.4：**逐切片**跑 BANKSY 全网格，落盘域分配与稳定性
#
# 口径（NICHE_PREREG.md §3.1–§3.4 已签 ＋ **§13 S4-REV/S6-REV 已签，2026-10-02**；本脚本不新开任何口径）：
#   基因轴 = §3.1 冻结的 HVG 3000（results/10_niche/hvg_3000.txt），全 56 张共用同一套。
#   归一化 = CP10K+log1p（S2）= NormalizeBanksy(norm_factor=1e4, log_norm=TRUE, base=10, pseudocount=1)
#   网格   = λ{0,0.2,0.5,0.8}(§13 S4-REV；主 0.2＝官方 Visium v1/v2 域分割推荐) × k_geom{6,18}(S5；主 18)
#            × resolution{0.5,0.6,0.8}(§13) × seed{0..4}；** λ=1.0 已退出**（旧网格曾误选它）
#   AGF    = compute_agf/use_agf 两版**分别调用**（S13：两版都留档，M6-7）
#   坐标   = 厂商 pxl_col/pxl_row_in_fullres 平移成非负（§3.2.1 实测 = 唯一正确写法）
#
# 本脚本**逐切片只跑一张**，输出该张的全部域分配 + 逐档稳定性统计；跨切片共识在 03_ 里做。
#
# 跑法（每张切片一个进程；AGF 两版各一次）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   setsid nohup Rscript 10_niche/02_banksy_grid.R --slide GSM9226190_P10_LUAD --agf TRUE \
#     < /dev/null > logs/nice_grid_P10_LUAD_agfT.log 2>&1 &

suppressMessages({
  library(Banksy)
  library(Matrix)
  library(Seurat)
  library(data.table)
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
if (is.null(SLIDE) || is.na(SLIDE)) stop("必须给 --slide <GSM..._P..._...>", call. = FALSE)
if (!toupper(AGF_C) %in% c("TRUE", "FALSE"))
  stop("--agf 只能是 TRUE 或 FALSE，收到：", AGF_C, call. = FALSE)
USE_AGF <- toupper(AGF_C) == "TRUE"

ROOT   <- "/home/eto/luad_v2"
VISIUM <- file.path(ROOT, "data/visium_spatial")
RES    <- file.path(ROOT, "results/08_spatial_deconv")
NICHE  <- file.path(ROOT, "results/10_niche")
HVG_F  <- get_arg("--hvg", file.path(NICHE, "hvg_3000.txt"))
AGF_TAG<- if (USE_AGF) "agfT" else "agfF"
SMOKE  <- "--smoke" %in% .a
OUTD   <- if (SMOKE) file.path(NICHE, "banksy_smoke_grid", SLIDE) else
                      file.path(NICHE, "banksy", SLIDE)
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

NRM_FACTOR  <- 1e4          # [已签 S2]
NPCS        <- 20L          # [已签]
K_NEIGHBORS <- 50L          # [已签]
SPATIAL_MODE<- "kNN_median" # [已签]
LAM_GRID    <- c(0, 0.2, 0.5, 0.8)          # [§13 S4-REV] 主 0.2（官方 Visium v1/v2 域分割推荐）；对照 0；敏感性 {0.5,0.8}；1.0 退出
RES_GRID    <- c(0.5, 0.6, 0.8)             # [§13 S4-REV] 官方两例 0.55/0.6 ＋ 源论文 0.5–0.8
KG_GRID     <- c(6L, 18L)                   # [已签 S5] 主 18（官方）/ 敏感性 6
SEEDS       <- 0:4                          # [已签]
N_HVG       <- 3000L
N_PERM      <- 200L         # 闸 2/C3 零分布置换数（C3 是硬门，见 §13.2/§13.6）
MIN_FREE_GB <- 25           # 内存守卫线（背景还有病理进程，见协作准则）
CN_TAG      <- if (USE_AGF) "M1" else "M0"  # ClusterBanksy 列名里的 M0/M1

## --smoke：只验证脚本能端到端跑通（缩档 + 独立输出目录），**产物不得当结论**
if (SMOKE) {
  KG_GRID <- 18L; LAM_GRID <- 0.2; RES_GRID <- c(0.5, 0.8); SEEDS <- 0:1; N_PERM <- 20L
}

## —————————————————————————————————————————————————————————————
## 二、工具
## —————————————————————————————————————————————————————————————
t_start <- Sys.time()
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
free_gb <- function() {
  x <- tryCatch(readLines("/proc/meminfo"), error = function(e) NULL)
  if (is.null(x)) return(Inf)
  v <- x[grepl("^MemAvailable:", x)]
  if (!length(v)) return(Inf)
  as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1", v[1])) / 1048576
}
peak_used_gb <- function() {
  x <- tryCatch(readLines("/proc/meminfo"), error = function(e) NULL)
  if (is.null(x)) return(NA_real_)
  g <- function(k) { v <- x[grepl(paste0("^", k, ":"), x)]
    if (!length(v)) return(NA_real_); as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1", v[1])) / 1048576 }
  g("MemTotal") - g("MemAvailable")
}
read_gz <- function(p, fn) { con <- gzfile(p, "rt"); on.exit(close(con)); fn(con) }
qq <- function(x, p) { x <- x[!is.na(x)]; if (!length(x)) NA_real_ else unname(quantile(x, p)) }
ir <- function(x) as.integer(round(x))   # 🔴 sprintf %d 只吃整数型 double；非整数会直接报错
guard_mem <- function(where) {
  f <- free_gb()
  if (f < MIN_FREE_GB)
    stop(sprintf("%s：可用内存 %.1f GB < 守卫线 %.0f GB ⇒ 硬停，不抢背景进程内存",
                 where, f, MIN_FREE_GB), call. = FALSE)
  f
}

step("==== M6 §3.2/§3.4 逐切片 BANKSY 网格：%s（compute_agf/use_agf = %s）====", SLIDE, USE_AGF)
step("Banksy %s / Seurat %s / R %s / 可用 %.1f GB",
     packageVersion("Banksy"), packageVersion("Seurat"),
     paste(R.version$major, R.version$minor, sep = "."), guard_mem("起跑"))

## 冻结 HVG 轴（硬约束：不许在本脚本里重选）
if (!file.exists(HVG_F)) stop("找不到冻结 HVG 清单：", HVG_F, "（先跑 01_hvg_freeze.R）", call. = FALSE)
HVG <- readLines(HVG_F)
if (length(HVG) != N_HVG || anyDuplicated(HVG))
  stop(sprintf("HVG 清单 %d 个（期望 %d，去重后 %d）⇒ 硬停",
               length(HVG), N_HVG, length(unique(HVG))), call. = FALSE)
step("冻结 HVG 轴：%d 个（%s）", length(HVG), HVG_F)

## —————————————————————————————————————————————————————————————
## 三、读片（口径与 00_banksy_smoke.R / 19_run_fastcnv_cohort.R 同源：N0 = 冻结掩码通过集）
## —————————————————————————————————————————————————————————————
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
po <- po[i, , drop = FALSE]
stopifnot(identical(po$barcode, bc_keep))

xy <- cbind(x = as.numeric(po$pxl_col_in_fullres), y = as.numeric(po$pxl_row_in_fullres))
xy <- sweep(xy, 2, apply(xy, 2, min))     # 平移成非负（距离不变；避开边缘负坐标陷阱）
rownames(xy) <- bc_keep

bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
m  <- readMM(file.path(d, "matrix.mtx.gz"))
stopifnot(nrow(m) == nrow(ft) + length(ab), ncol(m) == length(bc))
if (length(ab)) m <- m[-ab, , drop = FALSE]
colnames(m) <- bc
j <- match(bc_keep, bc)
if (anyNA(j)) stop(SLIDE, "：mask-pass barcode 不在表达矩阵里", call. = FALSE)
m <- m[, j, drop = FALSE]
g <- factor(ft$sym, levels = unique(ft$sym))
A <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1, dims = c(nlevels(g), length(g)))
m <- A %*% m
rownames(m) <- levels(g); colnames(m) <- bc_keep

miss <- setdiff(HVG, rownames(m))
if (length(miss)) stop(SLIDE, "：HVG 轴里有本片缺的基因 ", length(miss), " 个（轴对不上）", call. = FALSE)
EXPR <- as(m[HVG, , drop = FALSE], "CsparseMatrix")
rm(m, A, g); invisible(gc())
step("读入 HVG %d x spot %d（可用 %.1f GB）", nrow(EXPR), ncol(EXPR), free_gb())

## §3.5 画像用的表达：CP10K+log1p（与 BANKSY 同一归一化口径，S2）
NRM <- as(EXPR, "CsparseMatrix")
NRM@x <- log1p(NRM@x / rep(pmax(Matrix::colSums(EXPR), 1e-9), diff(NRM@p)) * NRM_FACTOR)
rownames(NRM) <- rownames(EXPR); colnames(NRM) <- colnames(EXPR)

## —————————————————————————————————————————————————————————————
## 四、建 bank（每 k_geom 一个）
## —————————————————————————————————————————————————————————————
build_bank <- function(k_geom) {
  bank <- BanksyObject(
    own.expr  = EXPR,
    cell.locs = data.frame(sdimx = xy[, 1], sdimy = xy[, 2], row.names = bc_keep),
    meta.data = data.frame(row.names = bc_keep))
  bank <- ComputeBanksy(bank, compute_agf = USE_AGF, spatial_mode = SPATIAL_MODE,
                        k_geom = k_geom, verbose = FALSE)
  bank <- NormalizeBanksy(bank, assay = "both", norm_factor = NRM_FACTOR,
                          log_norm = TRUE, pseudocount = 1, base = 10)
  ScaleBanksy(bank, assay = "both")
}
cluster_with <- function(bank, lambda, res, seed) {
  b <- ClusterBanksy(bank, lambda = lambda, use_agf = USE_AGF, pca = TRUE, npcs = NPCS,
                     method = "leiden", k.neighbors = K_NEIGHBORS, resolution = res,
                     num.cores = 1, seed = seed, verbose = FALSE)
  cn <- paste0("clust_", CN_TAG, "_lam", lambda, "_k", K_NEIGHBORS, "_res", res)
  if (!cn %in% colnames(b@meta.data))
    stop("找不到聚类列 ", cn, "；实际列：",
         paste(grep("^clust_", colnames(b@meta.data), value = TRUE), collapse = " "), call. = FALSE)
  as.integer(b@meta.data[[cn]])
}
## 闸 2：逐 spot 的 6 个真邻居中同域的比例
## 🔴 邻居索引只依赖坐标（Visium 六角栅格的 6 个物理近邻），**与 k_geom 无关** ⇒ 整片算一次。
##    （原先写在 coherence() 内部，会被 null_coherence 每轮重算 200 次，纯浪费。）
NB6 <- dbscan::kNN(x = xy, k = 6L)$id
coherence <- function(lab) mean(rowMeans(matrix(lab[NB6] == lab, nrow = length(lab))))
null_coherence <- function(lab, n_perm = N_PERM, seed = 1L) {
  set.seed(seed)
  v <- replicate(n_perm, coherence(sample(lab)))
  list(mean = mean(v), q05 = unname(quantile(v, 0.05)), q95 = unname(quantile(v, 0.95)),
       sd = sd(v))
}
## 域 → 平均表达画像（§3.5 第 2 步；只用表达）
domain_profile <- function(lab) {
  lv <- sort(unique(lab))
  M <- vapply(lv, function(k) Matrix::rowMeans(NRM[, lab == k, drop = FALSE]), numeric(nrow(NRM)))
  if (is.null(dim(M))) M <- matrix(M, ncol = length(lv))
  rownames(M) <- rownames(NRM); colnames(M) <- as.character(lv)
  M
}

## —————————————————————————————————————————————————————————————
## 五、主流程
## —————————————————————————————————————————————————————————————
DOM  <- list()   # 域分配（长表）
STT  <- list()   # 逐档统计（一行 = 一个 kg×λ×res）
PRF  <- list()   # 画像（一行 = 一个 kg×λ×res×域）

for (kg in KG_GRID) {
  guard_mem(sprintf("建 bank kg=%d 前", kg))
  t0 <- Sys.time()
  bank <- build_bank(kg)
  step("bank 就绪 k_geom=%d compute_agf=%s（%.1f min，已用 %.1f GB / 可用 %.1f GB）",
       kg, USE_AGF, elapsed(t0), peak_used_gb(), free_gb())

  for (lam in LAM_GRID) {
    guard_mem(sprintf("kg=%d λ=%.1f PCA 前", kg, lam))
    t0 <- Sys.time()
    bankP <- RunBanksyPCA(bank, lambda = lam, use_agf = USE_AGF, npcs = NPCS, verbose = FALSE)
    step("  PCA λ=%.1f（%.1f min，可用 %.1f GB）", lam, elapsed(t0), free_gb())

    for (r in RES_GRID) {
      labs <- lapply(SEEDS, function(sd) cluster_with(bankP, lam, r, sd))
      M <- do.call(cbind, labs)
      if (length(SEEDS) >= 2) {
        cm <- combn(length(SEEDS), 2)
        ari <- apply(cm, 2, function(j) mclust::adjustedRandIndex(M[, j[1]], M[, j[2]]))
        n_ident <- sum(apply(cm, 2, function(j) identical(M[, j[1]], M[, j[2]])))
        n_pair <- length(ari)
      } else {
        cm <- NULL; ari <- NA_real_; n_ident <- NA_integer_; n_pair <- 0L
      }
      lab0 <- labs[[1]]
      tb <- table(lab0)
      cg <- coherence(lab0)
      ng <- null_coherence(lab0)

      DOM[[length(DOM) + 1]] <- data.table(
        barcode = rep(bc_keep, length(SEEDS)), k_geom = kg, lambda = lam, resolution = r,
        seed = rep(SEEDS, each = length(bc_keep)), domain = as.integer(as.vector(M)))
      STT[[length(STT) + 1]] <- data.table(
        slide = SLIDE, k_geom = kg, lambda = lam, resolution = r,
        use_agf = USE_AGF, compute_agf = USE_AGF,
        n_domain = length(tb), n_singleton = sum(tb <= 2), n_lt10 = sum(tb < 10),
        size_min = min(tb), size_med = median(as.numeric(tb)), size_max = max(tb),
        ari_mean = mean(ari), ari_min = min(ari), ari_max = max(ari),
        n_seed_pairs = n_pair, n_seed_identical = n_ident,
        coh_mean = cg, null_coh_mean = ng$mean, null_coh_q95 = ng$q95,
        coh_z = (cg - ng$mean) / ng$sd, n_perm = N_PERM)
      PF <- domain_profile(lab0)                       # 基因 x 域
      PFd <- as.data.table(t(PF))                      # 域 x 基因
      PFd[, `:=`(k_geom = kg, lambda = lam, resolution = r,
                 domain = as.integer(colnames(PF)))]
      setcolorder(PFd, c("k_geom", "lambda", "resolution", "domain"))
      PRF[[length(PRF) + 1]] <- PFd

      step("  kg=%2d λ=%.1f res=%.1f → 域 %3d（<10 的 %2d，中位 %4s）跨种子 ARI 均 %.3f/最低 %.3f；连贯 %.3f（零分布 %.3f，z=%.1f）",
           ir(kg), lam, r, length(tb), ir(sum(tb < 10)), ir(median(as.numeric(tb))),
           mean(ari), min(ari), cg, ng$mean, (cg - ng$mean) / ng$sd)
    }
    rm(bankP); invisible(gc())
  }
  rm(bank); invisible(gc())
}

## —————————————————————————————————————————————————————————————
## 六、落盘
## —————————————————————————————————————————————————————————————
D <- rbindlist(DOM); S <- rbindlist(STT); P <- rbindlist(PRF)

## 🔴 原子落盘：先写 <file>.tmp，成功后再改名成正式名。
##    目的是「被中途杀掉（OOM/误杀/守护收任务）时，正式文件要么是完整旧版、要么是完整新版，
##    绝不出现半截文件被下游当成完整结果」。半截的只可能是 .tmp，续跑会重写它。
atomic_write <- function(x, path, gz = FALSE) {
  tmp <- paste0(path, ".tmp")
  if (gz) { con <- gzfile(tmp, "wt"); write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con) }
  else    { write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE) }
  if (!file.rename(tmp, path)) stop("原子改名失败：", tmp, " -> ", path)
}

## 按行写（给 JSON 这类**不是 data.frame** 的文本用；write.table 会在 json 对象上报
## "cannot coerce class json to a data.frame"）
atomic_write_lines <- function(txt, path) {
  tmp <- paste0(path, ".tmp")
  writeLines(txt, tmp)
  if (!file.rename(tmp, path)) stop("原子改名失败：", tmp, " -> ", path)
}

atomic_write(D, file.path(OUTD, sprintf("domains_%s.tsv.gz", AGF_TAG)), gz = TRUE)
atomic_write(S, file.path(OUTD, sprintf("stats_%s.tsv", AGF_TAG)), gz = FALSE)
atomic_write(P, file.path(OUTD, sprintf("profiles_%s.tsv.gz", AGF_TAG)), gz = TRUE)

man <- list(
  script = "10_niche/02_banksy_grid.R",
  prereg = "10_niche/NICHE_PREREG.md §3.1–§3.4（已签）",
  slide = SLIDE, n_spot = length(bc_keep),
  gene_set = sprintf("HVG %d (%s)", length(HVG), HVG_F),
  normalization = "CP10K + log1p (S2)",
  banksy_version = as.character(packageVersion("Banksy")),
  coord_system = "pxl_col/pxl_row_in_fullres, shifted non-negative (§3.2.1)",
  compute_agf = USE_AGF, use_agf = USE_AGF,
  lambda = LAM_GRID, k_geom = KG_GRID, resolution = RES_GRID, seeds = SEEDS,
  npcs = NPCS, k_neighbors = K_NEIGHBORS, spatial_mode = SPATIAL_MODE,
  n_perm_gate2 = N_PERM,
  n_domain_total_rows = nrow(S),
  files = c(sprintf("domains_%s.tsv.gz", AGF_TAG), sprintf("stats_%s.tsv", AGF_TAG),
            sprintf("profiles_%s.tsv.gz", AGF_TAG)),
  started = format(t_start, "%Y-%m-%d %H:%M:%S"),
  finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  minutes = round(elapsed(t_start), 2),
  note = "逐切片网格；跨切片共识在 03_consensus.R。本文件不含任何生物学结论。")
atomic_write_lines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE),
                   file.path(OUTD, sprintf("manifest_%s.json", AGF_TAG)))

step("落盘：%s（%d 行）/ stats（%d 行）/ profiles（%d 行）", sprintf("domains_%s.tsv.gz", AGF_TAG),
     nrow(D), nrow(S), nrow(P))
step("==== %s %s 结束，用时 %.1f 分钟 ====", SLIDE, AGF_TAG, elapsed(t_start))
