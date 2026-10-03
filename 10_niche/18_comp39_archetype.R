#!/usr/bin/env Rscript
# 18_comp39_archetype.R —— M6 §20：RCTD comp39 路线**能不能比 K*=7 更细**
#
# 设计：**严格照搬 13_composition_arm.R 的 comp39 臂**（§14.3 D4 已签口径），
#   唯一区别 = 这次把**逐域 39 维组成画像**落盘，以便跑**与 §17 完全相同的** K 扫描归一。
#   ⇒ 两边（BANKSY 域 vs RCTD 组成域）的"跨切片能归并出几个稳定类型"才可比。
#
# 🔴 与 17_rctd_finer_k.R 的区别（那个是**非严谨探针，已停**）：
#   探针用 k-means + 跨切片 pooled + 抽样 ⇒ 聚出来主要是"切片"不是"生态位"，数不可用。
#   本脚本用 **per-slide Seurat Leiden on kNN(k=50)** + 5 种子共识，与流水线一致。
#
# --slide <name>  单张：comp39 5 种子共识划分 → 逐域 39 维组成 → 原子写 d11_parts/<slide>.tsv
# --combine       汇总：与 §17 同款 K 扫描（PCA + k-means + bootstrap 80%x100）

suppressMessages({ library(data.table); library(Seurat); library(dbscan); library(mclust) })

ROOT  <- "/home/eto/luad_v2"
BN    <- file.path(ROOT, "results/10_niche/banksy")
DEC   <- file.path(ROOT, "results/08_spatial_deconv")
VIS   <- file.path(ROOT, "data/visium_spatial")
OUT   <- file.path(ROOT, "results/10_niche/kstar_diag")
PARTS <- file.path(OUT, "d11_parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
MAIN  <- list(k_geom = 18L, lambda = 0.2, resolution = 0.5)
SEEDS <- 0:4; NPCS <- 20L; K_NEIGHBORS <- 50L; K_GEOM <- 18L; RES_MAIN <- 0.5
STAB_MIN <- 0.90

.a <- commandArgs(trailingOnly = TRUE)
ga <- function(k) { i <- match(k, .a); if (is.na(i) || i == length(.a)) NULL else .a[i + 1] }
SLIDE <- ga("--slide"); COMBINE <- "--combine" %in% .a
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())
free_gb <- function() as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1",
                grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)[1])) / 1048576

## ————————————— 单张：comp39 共识域 + 组成画像 —————————————
one_slide <- function(s) {
  if (free_gb() < 20) stop("可用内存 < 20 GB ⇒ 硬停", call. = FALSE)
  d <- fread(file.path(BN, s, "domains_agfT.tsv.gz"), sep = "\t")
  bc <- d[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution &
          seed == min(seed)]$barcode
  w <- fread(file.path(DEC, "rctd_d", "per_slide", sprintf("%s.weights.tsv.gz", s)), sep = "\t", header = TRUE)
  setnames(w, 1L, "barcode"); X <- as.matrix(w[, -1L, with = FALSE]); X[!is.finite(X)] <- 0
  rownames(X) <- w$barcode
  i <- match(bc, rownames(X)); if (anyNA(i)) stop(s, "：组成缺 spot", call. = FALSE)
  X <- X[i, , drop = FALSE]

  po <- read.csv(file.path(VIS, s, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
  po <- po[po$in_tissue == 1L, , drop = FALSE]; j <- match(bc, po$barcode)
  xy <- cbind(po$pxl_col_in_fullres[j], po$pxl_row_in_fullres[j]); xy <- sweep(xy, 2, apply(xy, 2, min))

  NN  <- dbscan::kNN(x = xy, k = K_GEOM)$id
  AGF <- t(vapply(seq_len(nrow(X)), function(r) apply(X[NN[r, ], , drop = FALSE], 2, median),
                  numeric(ncol(X))))
  F2 <- cbind(X, AGF)
  colnames(F2) <- sprintf("f%03d", seq_len(ncol(F2)))   # 纯 ASCII，防 Seurat 洗名
  F2 <- scale(F2); F2[!is.finite(F2)] <- 0; rownames(F2) <- bc

  so <- CreateSeuratObject(counts = t(as.matrix(F2)))
  so[["RNA"]]$data <- t(as.matrix(F2)); VariableFeatures(so) <- colnames(F2)
  so <- ScaleData(so, verbose = FALSE)
  np_use <- min(NPCS, ncol(F2) - 1L)
  so <- RunPCA(so, npcs = np_use, verbose = FALSE, seed.use = 1L)
  so <- FindNeighbors(so, dims = seq_len(np_use), k.param = K_NEIGHBORS, verbose = FALSE)
  L <- suppressWarnings(sapply(SEEDS, function(sd)
    as.integer(FindClusters(so, resolution = RES_MAIN, algorithm = 4, random.seed = sd,
                            verbose = FALSE)$seurat_clusters)))
  if (is.null(dim(L))) L <- matrix(L, ncol = 1)

  ## 5 种子共识（与 §16/§17 同款：共分配矩阵 + 平均连接，切到各种子域数的中位）
  C <- matrix(0, nrow(L), nrow(L))
  for (k in seq_len(ncol(L))) { O <- outer(L[, k], L[, k], "=="); C <- C + O; rm(O) }
  C <- C / ncol(L)
  K <- as.integer(median(apply(L, 2, function(x) length(unique(x)))))
  dom <- cutree(hclust(as.dist(1 - C), method = "average"), k = K)
  rm(C, L); gc(verbose = FALSE)

  M <- rowsum(X, dom) / as.numeric(table(dom))
  r <- data.table(domain = as.integer(rownames(M)), M)
  r[, `:=`(slide = s, K_consensus = K, n_spot = as.integer(table(dom)),
           stage = sub("-1$", "", sub("^GSM[0-9]+_P[0-9]+_", "", s)))]
  r
}

if (!is.null(SLIDE)) {
  r <- one_slide(SLIDE)
  tmp <- file.path(PARTS, sprintf("%s.tsv.tmp", SLIDE)); fwrite(r, tmp, sep = "\t")
  if (!file.rename(tmp, file.path(PARTS, sprintf("%s.tsv", SLIDE)))) stop("原子改名失败", call. = FALSE)
  step("落盘 %s（%d 域）", SLIDE, nrow(r)); quit(save = "no", status = 0)
}
if (!COMBINE) stop("需要 --slide 或 --combine", call. = FALSE)

## ————————————— 汇总：与 §17 同款 K 扫描 —————————————
sl <- list.files(BN, pattern = "^GSM")
pf <- file.path(PARTS, sprintf("%s.tsv", sl)); have <- file.exists(pf)
step("D11 汇总：应有 %d / 已有 %d", length(sl), sum(have))
if (any(!have)) step("  ! 缺：%s", paste(sl[!have], collapse = " "))
D <- rbindlist(lapply(pf[have], fread), fill = TRUE)
META <- c("domain", "slide", "K_consensus", "n_spot", "stage")
SUB <- setdiff(names(D), META); stopifnot(length(SUB) == 39L)
step("comp39 域总数 %d（%d 维）", nrow(D), length(SUB))

X <- as.matrix(D[, ..SUB]); X <- X / pmax(rowSums(X), 1e-9)
Z <- scale(X); Z[!is.finite(Z)] <- 0
pc <- prcomp(Z, center = FALSE, scale. = FALSE); E <- pc$x[, seq_len(min(10L, ncol(pc$x)))]
set.seed(20261003)
ks <- 2:12
stab <- sapply(ks, function(k) {
  full <- kmeans(E, centers = k, nstart = 10, iter.max = 100)$cluster
  mean(replicate(100, { i <- sample(nrow(E), floor(0.8 * nrow(E)))
    a <- kmeans(E[i, , drop = FALSE], centers = k, nstart = 10, iter.max = 100)$cluster
    mclust::adjustedRandIndex(a, full[i]) }))
})
ST <- data.table(K = ks, stab = stab)
fwrite(ST, file.path(OUT, "d11_comp39_k_stability.tsv"), sep = "\t")
ok <- ST[stab >= STAB_MIN]
step("K 稳定性：%s", paste(sprintf("K=%d %.3f", ST$K, ST$stab), collapse = "  "))
step("⇒ 阈值 %.2f：%s", STAB_MIN,
     if (nrow(ok)) sprintf("最大达标 K = %d（%.3f）", ok[which.max(K)]$K, ok[which.max(K)]$stab)
     else sprintf("无任何 K 达标；描述性最优 K=%d（%.3f）", ST[which.max(stab)]$K, max(stab)))

con <- file(file.path(OUT, "d11_summary.txt"), "wt")
w <- function(...) cat(sprintf(...), file = con, append = TRUE)
w("=== §20 comp39（RCTD 39 亚型组成域）跨切片能归并出几个类型 ===\n%s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("comp39 共识域总数 %d（%d 张切片；逐切片共识 K 中位 %g）\n\n", nrow(D), uniqueN(D$slide), median(D[, .N, by = slide]$N))
w("方法**与 §17 完全相同**（PCA10 + k-means + bootstrap 80%% x100），仅特征/来源不同\n\n")
w("K 稳定性：%s\n", paste(sprintf("K=%d %.3f", ST$K, ST$stab), collapse = "  "))
w("§17（BANKSY 域）同期曲线：K=7 0.945，K=8 0.814\n")
w("⇒ %s\n", if (nrow(ok)) sprintf("comp39 最大达标 K = %d", ok[which.max(K)]$K)
             else "comp39 **无任何 K 达标**")
close(con); cat(readLines(file.path(OUT, "d11_summary.txt")), sep = "\n")
step("产物 → %s", OUT)
