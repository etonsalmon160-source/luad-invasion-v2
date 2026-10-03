#!/usr/bin/env Rscript
# 20_domain_expression.R —— M6 §22：给每个域补上**表达画像（marker 层）**
#
# 目的（用户 2026-10-03）：
#   ① 域里"具体是什么东西" —— 现在只有组成，缺 marker；
#   ② 在 **marker 空间**重跑 K 扫描 —— 这是第三条路线（BANKSY 嵌入 / RCTD 组成之外），
#      组成空间把"初始 T"和"耗竭 T"算成同一类，表达空间能分开 ⇒ 有可能比 K*=7 更细。
#
# 用 §16 已落盘的**共识域**（d8_parts/<slide>.tsv = barcode+domain），只读 Visium 原始矩阵。
# 特征 = 项目已冻结的 3000 HVG（results/10_niche/hvg_3000.txt）。
#
# --slide <name>  单张：CP10K+log1p → 取 HVG → 逐域均值 → 原子写 d13_parts/<slide>.tsv
# --combine       ① 逐域类型 marker（域 vs 其余）② marker 空间 K 扫描（与 §17/§20 同款）

suppressMessages({ library(data.table); library(Matrix) })

ROOT  <- "/home/eto/luad_v2"
VIS   <- file.path(ROOT, "data/visium_spatial")
OUT   <- file.path(ROOT, "results/10_niche/kstar_diag")
PARTS <- file.path(OUT, "d13_parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
HVGF  <- file.path(ROOT, "results/10_niche/hvg_3000.txt")
STAB_MIN <- 0.90

.a <- commandArgs(trailingOnly = TRUE)
ga <- function(k) { i <- match(k, .a); if (is.na(i) || i == length(.a)) NULL else .a[i + 1] }
SLIDE <- ga("--slide"); COMBINE <- "--combine" %in% .a
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())

one_slide <- function(s) {
  d <- fread(file.path(OUT, "d8_parts", sprintf("%s.tsv", s)), sep = "\t")   # barcode, domain
  fea <- fread(file.path(VIS, s, "filtered_feature_bc_matrix", "features.tsv.gz"), sep = "\t", header = FALSE)
  bar <- fread(file.path(VIS, s, "filtered_feature_bc_matrix", "barcodes.tsv.gz"), header = FALSE)
  M <- as(readMM(file.path(VIS, s, "filtered_feature_bc_matrix", "matrix.mtx.gz")), "CsparseMatrix")
  rownames(M) <- fea[[2]]; colnames(M) <- bar[[1]]           # 🔴 第2列才是符号名
  hvg <- intersect(fread(HVGF, header = FALSE)[[1]], rownames(M))
  M <- M[hvg, , drop = FALSE]
  toc <- colSums(M); toc[toc == 0] <- 1
  ## 🔴 `M %*% Diagonal()` 会**丢掉列名**（我踩过：之后 match 全 NA）⇒ 必须补回
  Mn <- M %*% Diagonal(x = 1e4 / toc); colnames(Mn) <- colnames(M); Mn@x <- log1p(Mn@x)
  i <- match(d$barcode, colnames(Mn))
  if (anyNA(i)) stop(s, "：barcode 不在表达矩阵里", call. = FALSE)
  Mn <- Mn[, i, drop = FALSE]
  dom <- d$domain
  I <- sparse.model.matrix(~ 0 + factor(dom))               # spots x domains
  S <- as.matrix(Mn %*% I)                                  # genes x domains（和）
  cnt <- as.numeric(table(dom))
  S <- sweep(S, 2, cnt, "/")
  r <- as.data.table(t(S))                                  # 🔴 要**域 × 基因**，不是基因 × 域
  setnames(r, rownames(S))
  r[, domain := as.integer(sub("^factor\\(dom\\)", "", colnames(S)))]
  r[, slide := s]
  r
}

if (!is.null(SLIDE)) {
  r <- one_slide(SLIDE)
  tmp <- file.path(PARTS, sprintf("%s.tsv.tmp", SLIDE)); fwrite(r, tmp, sep = "\t")
  if (!file.rename(tmp, file.path(PARTS, sprintf("%s.tsv", SLIDE)))) stop("原子改名失败", call. = FALSE)
  step("落盘 %s（%d 域 × %d 基因）", SLIDE, nrow(r), ncol(r) - 2L)
  quit(save = "no", status = 0)
}
if (!COMBINE) stop("需要 --slide 或 --combine", call. = FALSE)

## ————————————— 汇总 —————————————
sl <- list.files(file.path(ROOT, "results/10_niche/banksy"), pattern = "^GSM")
pf <- file.path(PARTS, sprintf("%s.tsv", sl)); have <- file.exists(pf)
step("D13 汇总：应有 %d / 已有 %d", length(sl), sum(have))
if (any(!have)) step("  ! 缺：%s", paste(sl[!have], collapse = " "))
D <- rbindlist(lapply(pf[have], fread), fill = TRUE)
A <- fread(file.path(OUT, "d7_domain_assign.tsv"))[, .(slide, domain, archetype)]
D <- merge(D, A, by = c("slide", "domain"))
GENES <- setdiff(names(D), c("slide", "domain", "archetype"))
step("域 %d × 基因 %d", nrow(D), length(GENES))
fwrite(D[, c("slide", "domain", "archetype"), with = FALSE],
       file.path(OUT, "d13_domain_archetype.tsv"), sep = "\t")

## ① 逐域类型的 marker：域 vs 其余域（Welch t 近似的效应量 + 排名）
##    用 log2(均值比) 的稳健版：log2((mean+eps)/(rest+eps))
E <- as.matrix(D[, ..GENES]); eps <- 1e-4
MK <- rbindlist(lapply(sort(unique(D$archetype)), function(a) {
  inx <- D$archetype == a
  m1 <- colMeans(E[inx, , drop = FALSE]); m0 <- colMeans(E[!inx, , drop = FALSE])
  lfc <- log2((m1 + eps) / (m0 + eps))
  o <- order(-lfc)[1:30]
  data.table(archetype = a, rank = seq_along(o), gene = GENES[o], log2fc = round(lfc[o], 3))
}))
fwrite(MK, file.path(OUT, "d13_archetype_markers.tsv"), sep = "\t")

## ② marker 空间 K 扫描（与 §17/§20 完全同款）
Z <- scale(E); Z[!is.finite(Z)] <- 0
pc <- prcomp(Z, center = FALSE, scale. = FALSE); EE <- pc$x[, seq_len(min(10L, ncol(pc$x)))]
set.seed(20261003); ks <- 2:12
stab <- sapply(ks, function(k) {
  full <- kmeans(EE, centers = k, nstart = 10, iter.max = 100)$cluster
  mean(replicate(100, { i <- sample(nrow(EE), floor(0.8 * nrow(EE)))
    a <- kmeans(EE[i, , drop = FALSE], centers = k, nstart = 10, iter.max = 100)$cluster
    mclust::adjustedRandIndex(a, full[i]) }))
})
ST <- data.table(K = ks, stab = stab)
fwrite(ST, file.path(OUT, "d13_marker_k_stability.tsv"), sep = "\t")
ok <- ST[stab >= STAB_MIN]
step("marker 空间 K 稳定性：%s", paste(sprintf("K=%d %.3f", ST$K, ST$stab), collapse = "  "))
step("⇒ %s", if (nrow(ok)) sprintf("最大达标 K = %d（%.3f）", ok[which.max(K)]$K, ok[which.max(K)]$stab)
             else sprintf("无达标 K；描述性最优 K=%d（%.3f）", ST[which.max(stab)]$K, max(stab)))
step("§17（RCTD 组成空间）同期：K=7 0.945 → K=8 0.814")

con <- file(file.path(OUT, "d13_summary.txt"), "wt")
w <- function(...) cat(sprintf(...), file = con, append = TRUE)
w("=== §22 域的表达画像（marker 层）%s ===\n\n", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("域 %d × HVG %d（冻结 3000 HVG）\n\n", nrow(D), length(GENES))
w("[1] marker 空间 K 稳定性（与 §17/§20 同款 PCA10+kmeans+bootstrap 80x100）\n%s\n\n",
  paste(sprintf("K=%d %.3f", ST$K, ST$stab), collapse = "  "))
w("对照：RCTD 组成空间 K=7 0.945 / K=8 0.814\n")
w("⇒ %s\n\n", if (nrow(ok)) sprintf("marker 空间最大达标 K = %d", ok[which.max(K)]$K) else "marker 空间无达标 K")
w("[2] 每类 top15 marker\n")
for (a in sort(unique(MK$archetype))) w("  D%d : %s\n", a, paste(MK[archetype == a][rank <= 15]$gene, collapse = ", "))
close(con); cat(readLines(file.path(OUT, "d13_summary.txt")), sep = "\n")
step("产物 → %s", OUT)
