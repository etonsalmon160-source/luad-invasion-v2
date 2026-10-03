#!/usr/bin/env Rscript
# 13_composition_arm.R —— M6 §14.3 D4：RCTD 组成域平行臂（**只读已落盘 RCTD 权重，不重跑 RCTD**）
#
# 预注册：NICHE_PREREG.md §14.3 D4（2026-10-02 已签）。本脚本**不产生任何主结果**。
#   粗档 comp6  = rctd_a（6 谱系）；细档 comp39 = rctd_d（39 亚型）。
#     两者是**同一批细胞、只换标签粒度**（n_cells/nnz 逐位相同）。
#   流程与 02_banksy_grid.R **同构**（唯一差别 = 特征空间由 BANKSY 嵌入换成 RCTD 组成）：
#     own + 空间 AGF(kNN median, k=18) → z-score → PCA(npcs=20，粗档受特征数限制取 11)
#     → FindNeighbors(k=50) → FindClusters(leiden, res=0.5) × 5 seeds → 逐种子对 ARI
#     连贯性 = 同款 6 个物理近邻（与 02 一致）
#   spot 集**强制对齐** BANKSY 主档（从 domains_agfT.tsv.gz 取 k=18/λ=0.2/res=0.5/seed=1 的 barcode）。
#
# 跑法（三种模式）：
#   单张：  Rscript 10_niche/13_composition_arm.R --slide GSM9226168_P1_AAH
#   汇总：  Rscript 10_niche/13_composition_arm.R --combine
#   并发由 13_composition_arm_supervisor.sh 负责（断点 + 内存守卫 + 并发 8）

suppressMessages({
  library(data.table); library(Matrix); library(Seurat); library(dbscan); library(mclust); library(ggplot2)
})

ROOT   <- "/home/eto/luad_v2"
VISIUM <- file.path(ROOT, "data/visium_spatial")
DECONV <- file.path(ROOT, "results/08_spatial_deconv")
NICHE  <- file.path(ROOT, "results/10_niche")
OUT    <- file.path(NICHE, "kstar_diag")
PARTS  <- file.path(OUT, "d4_parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)

ARMS   <- c(comp6 = "rctd_a", comp39 = "rctd_d")
MAIN   <- list(k_geom = 18L, lambda = 0.2, resolution = 0.5)
SEEDS  <- 0:4; NPCS <- 20L; K_NEIGHBORS <- 50L; K_GEOM <- 18L; N_PERM <- 200L; RES_MAIN <- 0.5
MIN_FREE_GB <- 30

.a <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k) { i <- match(k, .a); if (is.na(i) || i == length(.a)) NULL else .a[i + 1] }
SLIDE   <- get_arg("--slide")
COMBINE <- "--combine" %in% .a

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())
free_gb <- function() {
  x <- readLines("/proc/meminfo"); v <- x[grepl("^MemAvailable:", x)]
  as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1", v[1])) / 1048576
}

## —————————————————————————————————————————————————————————————
## 单张切片：两档各一行
## —————————————————————————————————————————————————————————————
banksy_spots <- function(slide) {
  d <- fread(file.path(NICHE, "banksy", slide, "domains_agfT.tsv.gz"), sep = "\t")
  d <- d[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution & seed == 1L]
  list(bc = d$barcode, lab = d$domain)
}

one_arm <- function(slide, arm, arm_dir, bc, lab_banksy, xy) {
  wf <- file.path(DECONV, arm_dir, "per_slide", sprintf("%s.weights.tsv.gz", slide))
  if (!file.exists(wf)) { step("  ! 缺 %s", wf); return(NULL) }
  w <- fread(wf, sep = "\t", header = TRUE); setnames(w, 1L, "barcode")
  i <- match(bc, w$barcode)
  if (anyNA(i)) { step("  ! %s：%d 个 barcode 不在 %s 里", slide, sum(is.na(i)), arm_dir); return(NULL) }
  X <- as.matrix(w[i, -1L, with = FALSE]); rownames(X) <- bc; X[!is.finite(X)] <- 0

  ## 空间 AGF：kNN median，k=18（与 02 同 k、同 mode）
  NN  <- dbscan::kNN(x = xy, k = K_GEOM)$id
  AGF <- t(vapply(seq_len(nrow(X)), function(r)
    apply(X[NN[r, ], , drop = FALSE], 2, median), numeric(ncol(X))))
  rownames(AGF) <- bc

  F2 <- cbind(X, AGF)
  ## 🔴 Seurat 会把特征名里的 `_` 洗成 `-`（中文更不可靠）⇒ 一律纯 ASCII 编号，
  ##    否则 SetAssayData 报 "No feature overlap"（我踩过）
  colnames(F2) <- sprintf("f%03d", seq_len(ncol(F2)))
  F2 <- scale(F2); F2[!is.finite(F2)] <- 0; rownames(F2) <- bc

  so <- CreateSeuratObject(counts = t(as.matrix(F2)))
  so[["RNA"]]$data <- t(as.matrix(F2)); VariableFeatures(so) <- colnames(F2)
  so <- ScaleData(so, verbose = FALSE)
  ## 🔴 粗档 own+agf 只有 12 个特征 ⇒ PCA 最多 11 个 PC，不能硬要 20
  np_use <- min(NPCS, ncol(F2) - 1L)
  so <- RunPCA(so, npcs = np_use, verbose = FALSE, seed.use = 1L)
  so <- FindNeighbors(so, dims = seq_len(np_use), k.param = K_NEIGHBORS, verbose = FALSE)
  L <- suppressWarnings(sapply(SEEDS, function(sd)
    as.integer(FindClusters(so, resolution = RES_MAIN, algorithm = 4, random.seed = sd,
                            verbose = FALSE)$seurat_clusters)))
  if (is.null(dim(L))) L <- matrix(L, ncol = 1)
  cm  <- combn(length(SEEDS), 2)
  ari <- apply(cm, 2, function(k) mclust::adjustedRandIndex(L[, k[1]], L[, k[2]]))

  NB6 <- dbscan::kNN(x = xy, k = 6L)$id
  coh <- function(l) mean(rowMeans(matrix(l[NB6] == l, nrow = length(l))))
  lab0 <- L[, 1]; set.seed(1L); nul <- replicate(N_PERM, coh(sample(lab0)))

  data.table(arm = arm, slide = slide, n_domain = length(unique(lab0)),
             n_singleton = sum(table(lab0) <= 2), n_lt10 = sum(table(lab0) < 10),
             ari_mean = mean(ari), ari_min = min(ari),
             coh_mean = coh(lab0), null_coh_q95 = unname(quantile(nul, 0.95)),
             ari_vs_banksy = mclust::adjustedRandIndex(lab0, lab_banksy))
}

one_slide <- function(slide) {
  f <- free_gb()
  if (f < MIN_FREE_GB) stop(sprintf("可用内存 %.0f GB < 守卫线 %d GB ⇒ 硬停", f, MIN_FREE_GB), call. = FALSE)
  bs <- banksy_spots(slide); bc <- bs$bc
  po <- read.csv(file.path(VISIUM, slide, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
  po <- po[po$in_tissue == 1L, , drop = FALSE]
  j  <- match(bc, po$barcode)
  if (anyNA(j)) stop(slide, "：坐标缺 spot", call. = FALSE)
  xy <- cbind(x = po$pxl_col_in_fullres[j], y = po$pxl_row_in_fullres[j])
  xy <- sweep(xy, 2, apply(xy, 2, min)); rownames(xy) <- bc
  ## 🔴 arm 列存**展示键**（comp6/comp39），不是目录名（rctd_a/rctd_d）——
  ##    曾把 ARMS[[a]] 当标签传进去 ⇒ 汇总按 comp* 匹配全落空（n=0）
  r <- rbindlist(lapply(names(ARMS), function(a) one_arm(slide, a, ARMS[[a]], bc, bs$lab, xy)),
                 fill = TRUE)
  step("  %s：%s", slide,
       paste(sprintf("%s K=%s ARI=%.3f", r$arm, r$n_domain, r$ari_mean), collapse = " | "))
  r
}

## —————————————————————————————————————————————————————————————
## 模式分派
## —————————————————————————————————————————————————————————————
if (!is.null(SLIDE)) {
  r <- one_slide(SLIDE)
  tmp <- file.path(PARTS, sprintf("%s.tsv.tmp", SLIDE))
  fwrite(r, tmp, sep = "\t")
  if (!file.rename(tmp, file.path(PARTS, sprintf("%s.tsv", SLIDE))))
    stop("原子改名失败", call. = FALSE)
  step("落盘 part：%s", SLIDE)
  quit(save = "no", status = 0)
}

if (!COMBINE) stop("需要 --slide <name> 或 --combine", call. = FALSE)

## ————————————— 汇总 —————————————
slides_all <- list.files(file.path(NICHE, "banksy"), pattern = "^GSM")
pf <- file.path(PARTS, sprintf("%s.tsv", slides_all))
have <- file.exists(pf)
step("D4 汇总：应有 %d 张，已有 %d 张", length(slides_all), sum(have))
if (any(!have)) step("  ! 缺：%s", paste(slides_all[!have], collapse = " "))
RES <- rbindlist(lapply(pf[have], fread), fill = TRUE)
fwrite(RES, file.path(OUT, "d4_composition_arm.tsv"), sep = "\t")

## 🔴 括号陷阱（我写错过一次）：fill=TRUE 必须落在 **rbindlist** 上，
##    写在 lapply 的右括号外侧会被 lapply 吃掉 ⇒ "unused argument (fill = TRUE)"
ST <- rbindlist(lapply(c("agfT", "agfF"), function(agf)
  rbindlist(lapply(slides_all, function(s) {
    d <- fread(file.path(NICHE, "banksy", s, sprintf("stats_%s.tsv", agf)), sep = "\t")
    d[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution][
      , `:=`(arm = agf, slide = s)]
  }), fill = TRUE)), fill = TRUE)

con <- file(file.path(OUT, "d4_summary.txt"), "wt")
w <- function(...) cat(sprintf(...), file = con, append = TRUE)
w("=== §14 D4 RCTD 组成域平行臂 汇总（%s）===\n\n", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("%-10s %6s %8s %8s %8s %8s %11s\n", "arm", "n", "K_med", "ARI_med", "ARI_min", "coh_med", "ARI_vs_bank")
for (a in c("agfT", "agfF")) {
  x <- ST[arm == a]
  w("%-10s %6d %8.0f %8.3f %8.3f %8.3f %11s\n", a, nrow(x), median(x$n_domain),
    median(x$ari_mean), median(x$ari_min), median(x$coh_mean), "-")
}
for (a in names(ARMS)) {
  x <- RES[arm == a]
  w("%-10s %6d %8.0f %8.3f %8.3f %8.3f %11.3f\n", a, nrow(x), median(x$n_domain),
    median(x$ari_mean), median(x$ari_min), median(x$coh_mean), median(x$ari_vs_banksy))
}
w("\n判据（§14.3）：组成域 ARI 若仍 < 0.90 ⇒ 病根在数据/区域，不在方法；\n")
w("              组成域 ARI 若 >= 0.90 ⇒ BANKSY 嵌入那一步是问题所在。\n")
close(con); cat(readLines(file.path(OUT, "d4_summary.txt")), sep = "\n")

P <- rbind(ST[, .(arm = arm, ari = ari_mean)], RES[, .(arm = arm, ari = ari_mean)])
g <- ggplot(P, aes(x = arm, y = ari, fill = arm)) +
  geom_hline(yintercept = 0.90, linetype = 2, colour = "red") +
  geom_boxplot(alpha = 0.75, outlier.size = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "cross-seed ARI (mean, per slide)",
       title = "D4: reproducibility of domains under three definitions",
       subtitle = "agfT/agfF = BANKSY expression (main cell); comp6/comp39 = RCTD composition; red dashed = 0.90 gate") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(OUT, "fig_d4_ari_by_definition.png"), g, width = 8, height = 5, dpi = 150)
step("产物 → %s", OUT)
