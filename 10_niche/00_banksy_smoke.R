#!/usr/bin/env Rscript
# M6 生态位臂 —— BANKSY 冒烟（NICHE_PREREG.md §3.2 前置条件 + §3.4 闸 1 + §3.2 网格）
#
# 🔴 本脚本**不产生任何生物学结论**。它只回答三个工具层面的问题：
#   A. 几何（§3.2 前置条件）
#      坐标进 kNN 之前是不是六角栅格？k_geom = 6 / 18 说的"环"是真的吗？
#      三种坐标写法对照：裸 array（错） / 各向同性化 array / 厂商 pxl（金标准）
#   B. 种子（§3.4 闸 1 是否空转）
#      leiden 跨种子 ARI 是不是恒等于 1？若是，闸 1 在本臂**不构成约束**
#      （GP5 记录里 3,000 细胞时该护栏恒 = 1.0 就是这个病，§3.4 已预警）
#   C. 粒度（用户问："组织域能做到多细"）
#      分辨率 {0.5,0.8,1.0,1.3,1.6} × k_geom {6,18} 扫一遍 → 域数 / 域尺寸谱 / 空间连贯性
#      + lambda 行（0 / 0.5 / 0.8 / 1.0）看"空间信息到底起没起作用"
#
# ⚠️ 与预注册的两处**已知偏离**（必须写进限制，不得当成正式跑）：
#   ① 基因集用**全部基因**（= §3.1 已签的"全 18,066 敏感性臂"），**不是** HVG 3000 主臂
#      （HVG 3000 须 56 张合池只选一次，见 §3.1；冒烟不做）
#   ② 只跑 2 张切片，不是 56 张
#   ⇒ 本脚本的数**不得**引用为正式结论。
#
# 跑法（只挂 fastcnv 这个 R 库，Banksy 在里面）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 10_niche/00_banksy_smoke.R

suppressMessages({
  library(Banksy)
  library(Matrix)
  library(data.table)
})

## —————————————————————————————————————————————————————————————
## 一、参数（逐项显式；来源 NICHE_PREREG.md §3.1/§3.2/§3.4）
## —————————————————————————————————————————————————————————————
ROOT   <- "/home/eto/luad_v2"
VISIUM <- file.path(ROOT, "data/visium_spatial")
RES    <- file.path(ROOT, "results/08_spatial_deconv")
OUTD   <- file.path(ROOT, "results/10_niche/smoke")

# 冒烟切片：一大一小，看 spot 数对粒度的影响
SLIDES <- c(P10_LUAD = "GSM9226190_P10_LUAD",   # 13,400 mask-pass（大队列里最大档）
            P1_LUAD  = "GSM9226169_P1_LUAD")    #  5,765 mask-pass

N_SMOKE_FULL <- "GSM9226190_P10_LUAD"  # 跑完整网格的那张

NRM_FACTOR  <- 1e4          # [已签 S2] CP10K 的 1e4；配 log_norm=TRUE = CP10K+log1p
NPCS        <- 20L          # [已签] npcs = 20
K_NEIGHBORS <- 50L          # [已签] 包默认
SPATIAL_MODE<- "kNN_median" # [已签] 包默认
LAM_GRID    <- c(0, 0.5, 0.8, 1.0)   # [已签 S4] 已删 0.2
LAM_MAIN    <- 0.8                   # [已签] 主候选来自高端
RES_GRID    <- c(0.5, 0.8, 1.0, 1.3, 1.6)   # [已签]
KG_GRID     <- c(6L, 18L)            # [已签 S5] 一环 / 两环
SEEDS       <- 0:4                   # [已签] 5 个种子
SEED_FIX    <- 0L

MIN_FREE_GB <- 40   # 内存守卫线（背景还有病理进程，见协作准则）

## —————————————————————————————————————————————————————————————
## 二、小工具
## —————————————————————————————————————————————————————————————
t_start <- Sys.time()
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
meminfo <- function() {
  x <- tryCatch(readLines("/proc/meminfo"), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  get1 <- function(k) {
    v <- x[grepl(paste0("^", k, ":"), x)]
    if (!length(v)) return(NA_real_)
    as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1", v[1])) / 1048576  # kB -> GB
  }
  c(total = get1("MemTotal"), avail = get1("MemAvailable"))
}
free_gb <- function() { m <- meminfo(); if (is.null(m) || is.na(m["avail"])) Inf else unname(m["avail"]) }
peak_gb <- function() { m <- meminfo(); if (is.null(m) || anyNA(m)) NA_real_ else unname(m["total"] - m["avail"]) }
read_gz <- function(p, fn) { con <- gzfile(p, "rt"); on.exit(close(con)); fn(con) }
`%||%` <- function(a, b) if (is.null(a)) b else a

dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

step("==== M6 生态位 BANKSY 冒烟开始 ====")
step("Banksy %s / Seurat %s / R %s", packageVersion("Banksy"),
     if (requireNamespace("Seurat", quietly = TRUE)) packageVersion("Seurat") else "n/a",
     paste(R.version$major, R.version$minor, sep = "."))
step("可用内存 %.1f GB（守卫线 %.0f GB）", free_gb(), MIN_FREE_GB)
if (free_gb() < MIN_FREE_GB)
  stop(sprintf("可用内存 %.1f GB < 守卫线 %.0f GB ⇒ 硬停，不抢背景病理进程的内存",
               free_gb(), MIN_FREE_GB), call. = FALSE)

## —————————————————————————————————————————————————————————————
## 三、读片（口径与 19_run_fastcnv_cohort.R 逐字同源：N0 = 冻结掩码通过集）
## —————————————————————————————————————————————————————————————
sm_all <- read.delim(gzfile(file.path(RES, "spot_mask.tsv.gz")), stringsAsFactors = FALSE)
sm_all <- sm_all[as.character(sm_all$pass) %in% c("TRUE", "true", "1"), ]
stopifnot(nrow(sm_all) > 0)
step("冻结掩码通过集 %d spot（56 张）", nrow(sm_all))

read_slide <- function(s, need_expr = TRUE) {
  d  <- file.path(VISIUM, s, "filtered_feature_bc_matrix")
  ft <- read_gz(file.path(d, "features.tsv.gz"), function(con)
    read.delim(con, header = FALSE, col.names = c("id", "sym", "type"),
               stringsAsFactors = FALSE))
  ab <- which(ft$type == "Antibody Capture")
  if (length(ab)) ft <- ft[-ab, , drop = FALSE]      # 抗体行计数量级 1e6，必剔
  if (length(ab)) step("  %s：剔抗体行 %d（剩 %d 行）", s, length(ab), nrow(ft))

  bc_keep <- sm_all$barcode[sm_all$slide == s]
  stopifnot(length(bc_keep) > 0)

  po <- read.csv(file.path(VISIUM, s, "spatial/tissue_positions.csv"),
                 stringsAsFactors = FALSE)
  po <- po[po$in_tissue == 1L, , drop = FALSE]
  i  <- match(bc_keep, po$barcode)
  if (anyNA(i)) stop(s, "：mask-pass 的 barcode 不在 tissue_positions.csv 里", call. = FALSE)
  po <- po[i, , drop = FALSE]
  stopifnot(identical(po$barcode, bc_keep))

  sc <- jsonlite::fromJSON(file.path(VISIUM, s, "spatial/scalefactors_json.json"))

  # 🔴 进 BANKSY 的坐标 = 厂商 pxl（全 56 张实测 top6 集合重叠 0.9995，是最优写法）。
  #    裸 array 坐标只有 0.8333（6 个邻居平均错 1 个），**不许用**；
  #    各向同性化 array 0.9827（差异是壳层边界 ties，不是几何错），留作敏感性。
  #    平移成非负：距离不变，但避免边缘负坐标（本项目已登记过的负索引陷阱）。
  xy <- cbind(x = as.numeric(po$pxl_col_in_fullres),
              y = as.numeric(po$pxl_row_in_fullres))
  xy <- sweep(xy, 2, apply(xy, 2, min))
  rownames(xy) <- bc_keep

  out <- list(slide = s, barcode = bc_keep, pos = po, xy = xy,
              spot_diam_fullres = as.numeric(sc$spot_diameter_fullres))
  if (!need_expr) return(out)

  bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
  m  <- readMM(file.path(d, "matrix.mtx.gz"))
  stopifnot(nrow(m) == nrow(ft) + length(ab), ncol(m) == length(bc))
  if (length(ab)) m <- m[-ab, , drop = FALSE]
  colnames(m) <- bc
  j <- match(bc_keep, bc)
  if (anyNA(j)) stop(s, "：mask-pass barcode 不在表达矩阵里", call. = FALSE)
  m <- m[, j, drop = FALSE]
  # 同名基因合并（与 19_ 同法）
  g <- factor(ft$sym, levels = unique(ft$sym))
  A <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1,
                    dims = c(nlevels(g), length(g)))
  m <- A %*% m
  rownames(m) <- levels(g)
  colnames(m) <- bc_keep
  out$expr <- as(m, "CsparseMatrix")
  out
}

## —————————————————————————————————————————————————————————————
## 四、Part A —— 几何：坐标进 kNN 前是不是六角栅格
## —————————————————————————————————————————————————————————————
# 三种坐标写法：
#   raw = (array_col, array_row)                   ← 天真写法（各向异性）
#   iso = (array_col/2, array_row*sqrt(3)/2)       ← 按六角几何各向同性化（无方向信息，纯几何）
#   px  = (pxl_col_in_fullres, pxl_row_in_fullres) ← 厂商像素坐标（金标准）
# 判据：在 px 空间量"物理距离"。真环成员的物理距离必须是一个尖峰（= 1 个 spot 间距）。
# 判据用**比例**，不用 all()：掩码会剔掉一部分真邻居，逐 spot 的 all() 必然 FALSE（假报警）。
# 安全的 quantile：全 NA 时返回 NA 而不是报错
qq <- function(x, p) { x <- x[!is.na(x)]; if (!length(x)) NA_real_ else unname(quantile(x, p)) }

geom_check <- function(sl) {
  po <- sl$pos
  b  <- sl$barcode
  xy_px  <- cbind(pxl_col = po$pxl_col_in_fullres, pxl_row = po$pxl_row_in_fullres)
  xy_raw <- cbind(a = po$array_col,                 b = po$array_row)
  xy_iso <- cbind(x = po$array_col / 2,             y = po$array_row * sqrt(3) / 2)
  rownames(xy_px) <- rownames(xy_raw) <- rownames(xy_iso) <- b

  KMAX <- max(KG_GRID)
  # 真值：px（厂商像素）空间的前 KMAX 名。它就是"物理上最近的那些 spot"。
  truth <- dbscan::kNN(x = xy_px, k = KMAX)
  pitch_px  <- median(truth$dist[, 1])
  um_per_px <- 55 / sl$spot_diam_fullres          # Visium spot 直径 = 55 µm
  pitch_um  <- pitch_px * um_per_px

  rows <- list(); summ <- list()
  for (sys in c("raw", "iso", "px")) {
    X  <- switch(sys, raw = xy_raw, iso = xy_iso, px = xy_px)
    kn <- dbscan::kNN(x = X, k = KMAX)
    # 该写法选出的邻居，其**物理**距离：先在真值表里查这些邻居的名次，再取该名次的距离。
    # （不能拿"邻居编号"当列号去索引 truth$dist —— 那是名次维度，只有 KMAX 列。）
    rk <- t(vapply(seq_len(nrow(kn$id)),
                   function(i) match(kn$id[i, ], truth$id[i, ]), integer(KMAX)))
    phys <- matrix(truth$dist[cbind(rep(seq_len(nrow(rk)), KMAX),
                                    as.vector(rk))], nrow = nrow(rk))
    # 与真值 top-6 / top-18 的集合重叠
    ov6  <- vapply(seq_len(nrow(kn$id)),
                   function(i) length(intersect(kn$id[i, 1:6],  truth$id[i, 1:6])),  1L)
    ov18 <- vapply(seq_len(nrow(kn$id)),
                   function(i) length(intersect(kn$id[i, ],     truth$id[i, ])),     1L)
    rows[[sys]] <- data.table(
      coord_system = sys, rank = rep(seq_len(KMAX), each = nrow(kn$id)),
      phys_dist_px = as.vector(phys))
    summ[[sys]] <- data.table(
      coord_system = sys, n_spot = nrow(kn$id),
      top6_set_overlap = mean(ov6) / 6,
      top18_set_overlap = mean(ov18) / KMAX,
      top6_p90_px = qq(phys[, 1:6], 0.90),
      top6_max_px = if (all(is.na(phys[, 1:6]))) NA_real_ else max(phys[, 1:6], na.rm = TRUE),
      rank7_18_min_px = if (all(is.na(phys[, 7:KMAX]))) NA_real_ else min(phys[, 7:KMAX], na.rm = TRUE),
      pct_neighbors_off_truth18 = mean(is.na(phys)),
      pct_spots_top6_changed = mean(ov6 < 6))
  }
  R <- rbindlist(rows)
  by_r <- R[, .(n = .N, n_off = sum(is.na(phys_dist_px)),
                med_px = qq(phys_dist_px, 0.50),
                p10_px = qq(phys_dist_px, 0.10),
                p90_px = qq(phys_dist_px, 0.90)),
            by = .(coord_system, rank)]
  S <- rbindlist(summ)

  list(by_rank = by_r, summary = S, pitch_px = pitch_px, pitch_um = pitch_um,
       um_per_px = um_per_px, spot_diam_fullres = sl$spot_diam_fullres,
       n_spot = length(b))
}

## —————————————————————————————————————————————————————————————
## 五、Part B/C —— 建 bank、PCA、聚类
## —————————————————————————————————————————————————————————————
build_bank <- function(sl, k_geom, compute_agf = TRUE) {
  bank <- BanksyObject(
    own.expr  = sl$expr,
    cell.locs = data.frame(sdimx = sl$xy[, 1],
                           sdimy = sl$xy[, 2],
                           row.names = sl$barcode),
    meta.data = data.frame(row.names = sl$barcode))
  bank <- ComputeBanksy(bank, compute_agf = compute_agf, spatial_mode = SPATIAL_MODE,
                        k_geom = k_geom, verbose = FALSE)
  bank <- NormalizeBanksy(bank, assay = "both", norm_factor = NRM_FACTOR,
                          log_norm = TRUE, pseudocount = 1, base = 10)
  bank <- ScaleBanksy(bank, assay = "both")
  bank
}

# PCA 只算一次（贵）；分辨率/种子复用同一个 reduction
prep_pca <- function(bank, lambda, use_agf = TRUE) {
  RunBanksyPCA(bank, lambda = lambda, use_agf = use_agf, npcs = NPCS, verbose = FALSE)
}
cluster_only <- function(bank, lambda, res, seed, use_agf = TRUE) {
  bank <- ClusterBanksy(bank, lambda = lambda, use_agf = use_agf, pca = TRUE, npcs = NPCS,
                        method = "leiden", k.neighbors = K_NEIGHBORS, resolution = res,
                        num.cores = 1, seed = seed, verbose = FALSE)
  cn <- paste0("clust_M", as.integer(use_agf), "_lam", lambda, "_k", K_NEIGHBORS, "_res", res)
  if (!cn %in% colnames(bank@meta.data))
    stop("找不到聚类列 ", cn, "；实际列：",
         paste(grep("^clust_", colnames(bank@meta.data), value = TRUE), collapse = " "),
         call. = FALSE)
  as.integer(bank@meta.data[[cn]])
}

# 空间连贯性：逐 spot 的 6 个真邻居中同域的比例（§3.4 闸 2）
coherence <- function(lab, sl) {
  nb <- dbscan::kNN(x = sl$xy, k = 6L)$id
  same <- rowMeans(matrix(lab[nb] == lab, nrow = length(lab)))
  list(mean = mean(same), q10 = unname(quantile(same, 0.10)),
       q50 = unname(quantile(same, 0.50)))
}
null_coherence <- function(lab, sl, n_perm = 200L, seed = 1L) {
  set.seed(seed)
  v <- replicate(n_perm, coherence(sample(lab), sl)$mean)
  list(mean = mean(v), q05 = unname(quantile(v, 0.05)), q95 = unname(quantile(v, 0.95)))
}

domain_table <- function(lab, sl) {
  tb <- table(lab)
  list(n_domain  = length(tb),
       n_singleton = sum(tb <= 2),
       n_lt10    = sum(tb < 10),
       size_min  = min(tb), size_med = median(as.numeric(tb)), size_max = max(tb))
}

## 出图：域的空间分布（英文标签；本机无 CJK 字体）
fig_domains <- function(lab_list, sl, path) {
  n <- length(lab_list)
  png(path, width = 320 * min(n, 3), height = 320 * ceiling(n / 3), res = 110)
  op <- par(mfrow = c(ceiling(n / 3), min(n, 3)), mar = c(2, 2, 3, 1))
  pal <- grDevices::hcl.colors(24, "Dark 3")
  for (i in seq_along(lab_list)) {
    lab <- lab_list[[i]]
    plot(sl$xy[, 1], -sl$xy[, 2],
         col = pal[(lab %% length(pal)) + 1L], pch = 15, cex = 0.45,
         asp = 1, xlab = "", ylab = "", axes = FALSE,
         main = names(lab_list)[i], cex.main = 1.0)
  }
  par(op); dev.off()
}

## —————————————————————————————————————————————————————————————
## 六、主流程
## —————————————————————————————————————————————————————————————
geom_out <- list(); gran_out <- list(); seed_out <- list()

for (tag in names(SLIDES)) {
  s  <- SLIDES[[tag]]
  step("---- %s（%s）----", tag, s)
  t0 <- Sys.time()
  sl <- read_slide(s, need_expr = TRUE)
  step("  读入 %d 基因 x %d spot（%.1f min，峰值 %.1f GB）",
       nrow(sl$expr), ncol(sl$expr), elapsed(t0), peak_gb())

  ## ---- Part A ----
  gm <- geom_check(sl)
  geom_out[[tag]] <- gm
  step("  A 几何：spot 间距 = %.1f px = %.1f µm（spot 直径 %.1f px ⇒ %.4f µm/px）",
       gm$pitch_px, gm$pitch_um, gm$spot_diam_fullres, gm$um_per_px)
  step("  A 前 6 名物理距离 p90 = %.1f px（应 ≈ 1 个间距 395），第 7–18 名最小 = %.1f px（应 ≈ 1.7 个间距）",
       gm$summary[coord_system == "px"]$top6_p90_px,
       gm$summary[coord_system == "px"]$rank7_18_min_px)
  for (i in seq_len(nrow(gm$summary))) {
    r <- gm$summary[i]
    step("  A %-3s 写法：top6 与真邻居集合重叠 %.4f（6 个里中 %s 个）⇒ 有 %.1f%% 的 spot 邻居表被改过",
         r$coord_system, r$top6_set_overlap, round(r$top6_set_overlap * 6, 2),
         100 * r$pct_spots_top6_changed)
    step("      其中选出的邻居里有 %.2f%% 根本不在物理最近 %d 个之内（坐标越差，这个数越大）",
         100 * r$pct_neighbors_off_truth18, max(KG_GRID))
  }
  write.table(gm$by_rank, gzfile(file.path(OUTD, sprintf("geometry_by_rank_%s.tsv.gz", tag))),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(gm$summary, gzfile(file.path(OUTD, sprintf("geometry_summary_%s.tsv.gz", tag))),
              sep = "\t", quote = FALSE, row.names = FALSE)

  ## ---- Part B/C ----
  full <- (s == N_SMOKE_FULL)
  kg_run <- if (full) KG_GRID else 6L
  lam_run <- if (full) LAM_GRID else LAM_MAIN

  for (kg in kg_run) {
    t0 <- Sys.time()
    bank <- build_bank(sl, kg, compute_agf = TRUE)
    step("  bank 就绪 k_geom=%d compute_agf=TRUE（%.1f min，峰值 %.1f GB）",
         kg, elapsed(t0), peak_gb())

    ## --- Part B：种子稳定性（只在主臂参数上做一次）---
    if (full && kg == 6L) {
      bankB <- prep_pca(bank, LAM_MAIN, use_agf = TRUE)
      labs <- list()
      for (sd in SEEDS) {
        lab <- cluster_only(bankB, LAM_MAIN, 1.0, sd, use_agf = TRUE)
        labs[[as.character(sd)]] <- lab
      }
      m <- sapply(labs, as.integer)
      cm <- combn(ncol(m), 2)
      ari <- apply(cm, 2, function(j) mclust::adjustedRandIndex(m[, j[1]], m[, j[2]]))
      ident <- apply(cm, 2, function(j) identical(m[, j[1]], m[, j[2]]))
      seed_out[[tag]] <- data.table(slide = tag, seed_a = SEEDS[cm[1, ]],
                                    seed_b = SEEDS[cm[2, ]], ari = ari, identical = ident)
      step("  B 跨种子 ARI：min=%.4f max=%.4f；逐位相同 %d/%d 对  ⇒ 闸 1 %s",
           min(ari), max(ari), sum(ident), length(ident),
           if (all(ident)) "**空转（恒等于 1）**" else "在起作用")
      # leiden 本身是否随机：同一图连跑两次
      x <- bankB@reduction[[paste0("pca_M1_lam", LAM_MAIN)]]$x[, seq_len(NPCS)]
      g1 <- Banksy:::getGraph(x, K_NEIGHBORS)
      set.seed(1); l1 <- leidenAlg::leiden.community(g1, resolution = 1)$membership
      l2 <- leidenAlg::leiden.community(g1, resolution = 1)$membership
      step("  B 同一图连跑两次 leiden：逐位相同 = %s ⇒ 随机性来源 = %s",
           identical(as.integer(l1), as.integer(l2)),
           if (identical(as.integer(l1), as.integer(l2))) "无（图是确定的）" else "leiden 自身")
    }

    ## --- Part C：粒度 ---
    for (lam in lam_run) {
      t0 <- Sys.time()
      bankC <- prep_pca(bank, lam, use_agf = TRUE)
      step("  PCA 就绪 lambda=%.1f（%.1f min，峰值 %.1f GB）", lam, elapsed(t0), peak_gb())
      lab_list <- list()
      for (r in RES_GRID) {
        lab <- cluster_only(bankC, lam, r, SEED_FIX, use_agf = TRUE)
        dt  <- domain_table(lab, sl)
        cg  <- coherence(lab, sl)
        ng  <- null_coherence(lab, sl)
        gran_out[[length(gran_out) + 1]] <- data.table(
          slide = tag, k_geom = kg, lambda = lam, resolution = r, use_agf = TRUE,
          n_domain = dt$n_domain, n_singleton = dt$n_singleton, n_lt10 = dt$n_lt10,
          size_min = dt$size_min, size_med = dt$size_med, size_max = dt$size_max,
          coh_mean = cg$mean, coh_q10 = cg$q10, coh_q50 = cg$q50,
          null_mean = ng$mean, null_q95 = ng$q95,
          coh_per_pitch = cg$mean, pitch_um = gm$pitch_um)
        lab_list[[sprintf("k_geom=%d  lambda=%.1f  resolution=%.1f  ->  %d domains",
                          kg, lam, r, dt$n_domain)]] <- lab
        step("  C %-6s kg=%2d lam=%.1f res=%.1f → 域 %3d（<10 spot 的 %2d，最小 %3d，中位 %4d）连贯 %.3f（零分布 %.3f）",
             tag, as.integer(kg), lam, r, as.integer(dt$n_domain), as.integer(dt$n_lt10),
             as.integer(round(dt$size_min)), as.integer(round(dt$size_med)),
             cg$mean, ng$mean)
      }
      if (lam == LAM_MAIN) {
        fig_domains(lab_list, sl,
                    file.path(OUTD, sprintf("fig_domains_%s_kg%d.png", tag, kg)))
      }
    }
    rm(bank); invisible(gc())
  }
  invisible(gc())
}

## —————————————————————————————————————————————————————————————
## 七、落盘
## —————————————————————————————————————————————————————————————
if (length(gran_out)) {
  G <- rbindlist(gran_out)
  write.table(G, gzfile(file.path(OUTD, "granularity.tsv.gz")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  step("落盘 granularity.tsv.gz（%d 行）", nrow(G))
}
if (length(seed_out)) {
  write.table(rbindlist(seed_out), gzfile(file.path(OUTD, "seed_ari.tsv.gz")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# 几何汇总
geo_txt <- paste0("{",
  paste(sprintf('"%s":{"n_spot":%d,"pitch_px":%.2f,"pitch_um":%.3f,"um_per_px":%.5f,"raw_top6_overlap":%.4f,"iso_top6_overlap":%.4f,"px_top6_overlap":%.4f,"raw_pct_spots_changed":%.4f}',
    names(geom_out), vapply(geom_out, function(z) as.integer(z$n_spot), 1L),
    vapply(geom_out, function(z) z$pitch_px, 1),
    vapply(geom_out, function(z) z$pitch_um, 1),
    vapply(geom_out, function(z) z$um_per_px, 1),
    vapply(geom_out, function(z) z$summary[coord_system == "raw"]$top6_set_overlap, 1),
    vapply(geom_out, function(z) z$summary[coord_system == "iso"]$top6_set_overlap, 1),
    vapply(geom_out, function(z) z$summary[coord_system == "px"]$top6_set_overlap, 1),
    vapply(geom_out, function(z) z$summary[coord_system == "raw"]$pct_spots_top6_changed, 1)),
    collapse = ","), "}")
writeLines(geo_txt, file.path(OUTD, "geometry_summary.json"))

man <- list(
  script = "10_niche/00_banksy_smoke.R",
  prerereg = "10_niche/NICHE_PREREG.md",
  banksy_version = as.character(packageVersion("Banksy")),
  run_started = format(t_start, "%Y-%m-%d %H:%M:%S"),
  run_minutes = round(elapsed(t_start), 1),
  peak_gb = round(peak_gb(), 1),
  gene_set = "ALL_GENES (signed sensitivity arm; NOT the HVG3000 main arm)",
  slides = unname(SLIDES),
  note = "smoke only - no biological conclusions")
writeLines(jsonlite::toJSON(man, pretty = TRUE, auto_unbox = TRUE),
           file.path(OUTD, "run_manifest.json"))

step("==== 完成，用时 %.1f min，峰值 %.1f GB ====", elapsed(t_start), peak_gb())
