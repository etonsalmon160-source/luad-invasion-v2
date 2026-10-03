#!/usr/bin/env Rscript
# 15_domain_archetype.R —— M6 §17：跨切片把域归并成"域类型"，如实报 K
#
# 目的（用户 2026-10-03 要求）：最终结论 = **几个域、每个域的细胞组成**。
#   这正是 §3.5 该做而静默回落掉了的那一步（§12-⑭）。本脚本重做，并遵守 §15.2：
#   **稳定性不达标就如实记「无稳定 K」**，不输出 which.max 那种退化解。
#
# 两段式（可续跑）：
#   --slide <name>  单张：5 种子共识划分 → 逐域完整组成（rctd_d 39 亚型 + rctd_a 6 谱系）
#                          → 原子写 d7_parts/<slide>.tsv      （断点单位 = 一张切片）
#   --combine       汇总：428 域 × 39 维 → PCA → K 扫描 + bootstrap 稳定性
#                          → 如实选 K（或报无稳定 K）→ 每类组成 + 图
#
# 跑法：R_LIBS=... Rscript 10_niche/15_domain_archetype.R --slide GSM... | --combine

suppressMessages({ library(data.table); library(Matrix); library(ggplot2) })

ROOT  <- "/home/eto/luad_v2"
BN    <- file.path(ROOT, "results/10_niche/banksy")
DEC   <- file.path(ROOT, "results/08_spatial_deconv")
OUT   <- file.path(ROOT, "results/10_niche/kstar_diag")
PARTS <- file.path(OUT, "d7_parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
MAIN  <- list(k_geom = 18L, lambda = 0.2, resolution = 0.5)
STAB_MIN <- 0.90          # 沿用 §3.4/§3.5 的 0.90（法则 3.2：数值不改）
set.seed(20261003)

.a <- commandArgs(trailingOnly = TRUE)
ga <- function(k) { i <- match(k, .a); if (is.na(i) || i == length(.a)) NULL else .a[i + 1] }
SLIDE <- ga("--slide"); COMBINE <- "--combine" %in% .a
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())
free_gb <- function() as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1",
                grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)[1])) / 1048576

## ————————————————— 单张 —————————————————
one_slide <- function(s) {
  if (free_gb() < 20) stop("可用内存 < 20 GB ⇒ 硬停", call. = FALSE)
  d <- fread(file.path(BN, s, "domains_agfT.tsv.gz"), sep = "\t")
  d <- d[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution]
  W <- dcast(d, barcode ~ seed, value.var = "domain"); bc <- W$barcode
  L <- as.matrix(W[, -1L, with = FALSE])
  C <- matrix(0, length(bc), length(bc))
  for (j in seq_len(ncol(L))) { O <- outer(L[, j], L[, j], "=="); C <- C + O; rm(O) }
  C <- C / ncol(L)
  K <- as.integer(median(apply(L, 2, function(x) length(unique(x)))))
  dom <- cutree(hclust(as.dist(1 - C), method = "average"), k = K)
  rm(C, L, W); gc(verbose = FALSE)

  out <- list()
  for (arm in c("rctd_d", "rctd_a")) {
    w <- fread(file.path(DEC, arm, "per_slide", sprintf("%s.weights.tsv.gz", s)), sep = "\t", header = TRUE)
    setnames(w, 1L, "barcode"); X <- as.matrix(w[, -1L, with = FALSE]); X[!is.finite(X)] <- 0
    rownames(X) <- w$barcode
    i <- match(bc, rownames(X))
    if (anyNA(i)) stop(s, "：", arm, " 缺 spot", call. = FALSE)
    X <- X[i, , drop = FALSE]
    M <- rowsum(X, dom) / as.numeric(table(dom))
    out[[arm]] <- data.table(domain = as.integer(rownames(M)), M)
    rm(w, X); gc(verbose = FALSE)
  }
  r <- merge(out$rctd_d, out$rctd_a, by = "domain", suffixes = c("", ".lin"))
  r[, `:=`(slide = s, K_consensus = K, n_spot = as.integer(table(dom)),
           stage = sub("-1$", "", sub("^GSM[0-9]+_P[0-9]+_", "", s)),
           patient = sub("^GSM[0-9]+_(P[0-9]+)_.*$", "\\1", s))]
  r
}

if (!is.null(SLIDE)) {
  r <- one_slide(SLIDE)
  tmp <- file.path(PARTS, sprintf("%s.tsv.tmp", SLIDE))
  fwrite(r, tmp, sep = "\t")
  if (!file.rename(tmp, file.path(PARTS, sprintf("%s.tsv", SLIDE)))) stop("原子改名失败", call. = FALSE)
  step("落盘 %s（%d 域）", SLIDE, nrow(r)); quit(save = "no", status = 0)
}
if (!COMBINE) stop("需要 --slide <name> 或 --combine", call. = FALSE)

## ————————————————— 汇总 —————————————————
sl <- list.files(BN, pattern = "^GSM")
pf <- file.path(PARTS, sprintf("%s.tsv", sl)); have <- file.exists(pf)
step("D7 汇总：应有 %d / 已有 %d", length(sl), sum(have))
if (any(!have)) { step("  ! 缺：%s", paste(sl[!have], collapse = " ")); }
D <- rbindlist(lapply(pf[have], fread), fill = TRUE)
## 🔴 我犯过的错：rctd_a 的 6 个粗谱系列名与 rctd_d 的 39 个细亚型**不重名**
##    ⇒ suffixes 不生效、`\.lin$` 匹配不到 ⇒ 45 个特征被当成"39 亚型"混着用（同一个东西数两遍）
LIN_NAMES <- c("B_浆", "T_NK", "上皮", "内皮", "成纤维", "髓系")
LIN <- intersect(LIN_NAMES, names(D))
SUB <- setdiff(names(D), c("domain", "slide", "K_consensus", "n_spot", "stage", "patient", LIN))
stopifnot(length(SUB) == 39L, length(LIN) == 6L)
step("域总数 %d（%d 细亚型维 + %d 粗谱系维，后者不进特征）", nrow(D), length(SUB), length(LIN))

X <- as.matrix(D[, ..SUB]); X <- X / pmax(rowSums(X), 1e-9)      # 行归一（域内组成）
Z <- scale(X); Z[!is.finite(Z)] <- 0
pc <- prcomp(Z, center = FALSE, scale. = FALSE)
nv <- min(10L, ncol(pc$x))
E <- pc$x[, seq_len(nv)]

## K 扫描 + bootstrap 稳定性（子样 80%，重抽 100 次）
ks <- 2:12
stab <- sapply(ks, function(k) {
  full <- kmeans(E, centers = k, nstart = 10, iter.max = 100)$cluster
  ari <- replicate(100, {
    i <- sample(nrow(E), floor(0.8 * nrow(E)))
    a <- kmeans(E[i, , drop = FALSE], centers = k, nstart = 10, iter.max = 100)$cluster
    mclust::adjustedRandIndex(a, full[i])
  })
  mean(ari)
})
ST <- data.table(K = ks, stab = stab)
fwrite(ST, file.path(OUT, "d7_k_stability.tsv"), sep = "\t")
ok <- ST[stab >= STAB_MIN]
step("K 稳定性：%s", paste(sprintf("K=%d %.3f", ST$K, ST$stab), collapse = "  "))
if (nrow(ok)) {
  KSTAR <- ok[which.max(K)]$K
  step("⇒ 达标 K 取最大：K* = %d（%.3f）", KSTAR, ST[K == KSTAR]$stab)
} else {
  KSTAR <- NA_integer_
  step("🔴 无任何 K 达到 %.2f ⇒ 按 §15.2 如实记「无稳定 K」，不静默挑 argmax", STAB_MIN)
}
## 无论如何都给出"描述性最优"供参考，但必须标明它不是达标解
KREF <- ST[which.max(stab)]$K
step("（参考：稳定性最高的是 K=%d，%.3f —— 仅作描述，**不作**达标口径）", KREF, max(stab))

## 取 K*；若无达标 K，仍按 KREF 出组成表但**文件名与首行都标注"未达标"**
KF <- if (is.na(KSTAR)) KREF else KSTAR
set.seed(1); cl <- kmeans(E, centers = KF, nstart = 25, iter.max = 100)$cluster
D[, archetype := cl]
## 🔴 计数一律**在同一次 by 聚合里**算出来，不要再 `as.integer(table(cl))` 按位置对齐
##    （那样一旦行序不同就整体错位——我踩过：19 个域被写成 46 张切片）
dg <- fread(file.path(ROOT, "results/10_niche/depth_guard/summary_per_slide.tsv"))
dg <- dg[agf == "agfT", .(depth_med = median(depth_med)), by = slide]
D <- merge(D, dg, by = "slide", all.x = TRUE)
comp <- D[, c(lapply(.SD, mean), list(n_domain = .N, n_slide = uniqueN(slide),
                                      depth_med = median(depth_med))),
          by = archetype, .SDcols = SUB]
setorder(comp, archetype)
comp[, top5 := apply(as.matrix(.SD[, ..SUB]), 1, function(v)
  paste(sprintf("%s:%.2f", names(sort(v, decreasing = TRUE))[1:5], sort(v, decreasing = TRUE)[1:5]),
        collapse = ";"))]
fwrite(comp, file.path(OUT, "d7_archetypes.tsv"), sep = "\t")
fwrite(D[, c("slide","stage","patient","domain","archetype","n_spot", LIN), with = FALSE],
       file.path(OUT, "d7_domain_assign.tsv"), sep = "\t")

con <- file(file.path(OUT, "d7_summary.txt"), "wt")
w <- function(...) cat(sprintf(...), file = con, append = TRUE)
w("=== §17 域类型（跨切片归并）%s ===\n\n", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("域总数 %d，来自 %d 张切片\n", nrow(D), uniqueN(D$slide))
w("K 稳定性（bootstrap 80%% x100，阈 %.2f）：%s\n", STAB_MIN,
  paste(sprintf("K=%d %.3f", ST$K, ST$stab), collapse = "  "))
w("⇒ %s\n\n", if (is.na(KSTAR)) sprintf("🔴 **无任何 K 达标** ⇒ 下面的 K=%d 只是稳定性最高者，**不是**达标解（§15.2）", KREF)
  else sprintf("K* = %d（%.3f）", KSTAR, ST[K == KSTAR]$stab))
w("[每类组成：39 亚型 top5（粗谱系不参与）]\n")
for (i in seq_len(nrow(comp))) w("  类 %2d  n=%3d 域 / %2d 张切片 / 深度中位 %5.0f : %s\n",
  comp$archetype[i], comp$n_domain[i], comp$n_slide[i], comp$depth_med[i], comp$top5[i])
w("\n[域类型 × 期别（域个数）]\n")
CT <- dcast(D[, .N, by = .(stage, archetype)], stage ~ archetype, value.var = "N", fill = 0)
CT[, stage := factor(stage, levels = c("Normal","AAH","AIS","MIA","LUAD"))]
print(CT[order(stage)], file = con)
w("\n[列合计（应等于各类域总数）] %s\n", paste(colSums(CT[, -1]), collapse = " "))
close(con); cat(readLines(file.path(OUT, "d7_summary.txt")), sep = "\n")

g <- ggplot(ST, aes(K, stab)) + geom_line(linewidth = .9) + geom_point(size = 2.4, colour = "#2C7FB8") +
  geom_hline(yintercept = STAB_MIN, linetype = 2, colour = "red") +
  labs(x = "number of domain archetypes (K)", y = "bootstrap ARI (mean, 80% x100)",
       title = "§17: stability of cross-slide domain archetypes",
       subtitle = sprintf("red dashed = %.2f gate; %s", STAB_MIN,
         if (is.na(KSTAR)) sprintf("no K passes; argmax K=%d is descriptive only", KREF) else sprintf("K*=%d", KSTAR))) +
  theme_bw(base_size = 12)
ggsave(file.path(OUT, "fig_d7_k_stability.png"), g, width = 7.5, height = 4.8, dpi = 150)
step("产物 → %s", OUT)
