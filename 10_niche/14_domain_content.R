#!/usr/bin/env Rscript
# 14_domain_content.R —— M6 §16：域内容剖析（域**到底是什么**）
#
# 背景：§4.1 闸 1 FAIL（跨种子 ARI 中位 0.765）⇒ 单一种子的划分不稳。
#   本脚本**不**挑某一个种子，而是对 5 个种子做**共识划分**：
#     共分配矩阵 C[i,j] = 5 个种子里 spot i 与 j 分在同一域的频次（0…1）
#     以 (1 - C) 为距离做平均连接层次聚类，切到 K = 该切片各种子域数的中位
#   这样每一域都要求"5 次里多数次数都待在一起"，是对不稳定性的正面处置。
#
# 主档：k_geom=18 / lambda=0.2 / res=0.5 / agfT（§13 已签固定档）。
# 域标签来源：rctd_a（6 大类，粗而稳）；精细组成另附 rctd_d（39 亚型）前 3 位。
# marker：profiles_agfT.tsv.gz 的域均值，**域内跨域 z 分数**取 top。
#
# 🔴 本脚本只读已落盘产物，不重跑 BANKSY；产物写 kstar_diag/。
# 跑法：R_LIBS=... Rscript 10_niche/14_domain_content.R

suppressMessages({ library(data.table); library(Matrix) })

ROOT  <- "/home/eto/luad_v2"
BN    <- file.path(ROOT, "results/10_niche/banksy")
DEC   <- file.path(ROOT, "results/08_spatial_deconv")
VIS   <- file.path(ROOT, "data/visium_spatial")
OUT   <- file.path(ROOT, "results/10_niche/kstar_diag")
MAIN  <- list(k_geom = 18L, lambda = 0.2, resolution = 0.5)

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())

slides <- list.files(BN, pattern = "^GSM")
step("==== §16 域内容剖析：共识划分 + 组成 + marker（%d 张）====", length(slides))

read_weights <- function(arm, slide) {
  f <- file.path(DEC, arm, "per_slide", sprintf("%s.weights.tsv.gz", slide))
  if (!file.exists(f)) return(NULL)
  w <- fread(f, sep = "\t", header = TRUE); setnames(w, 1L, "barcode")
  M <- as.matrix(w[, -1L, with = FALSE]); rownames(M) <- w$barcode
  M[!is.finite(M)] <- 0; M
}

one_slide <- function(slide) {
  ## ---- 1. 五个种子的域标签 → 共识划分 ----
  d <- fread(file.path(BN, slide, "domains_agfT.tsv.gz"), sep = "\t")
  d <- d[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution]
  if (uniqueN(d$seed) < 2L) { step("  ! %s 种子不足", slide); return(NULL) }
  ## 🔴 别用 d[match(bc, barcode)][seed == s]：match 在**含全部种子**的表里匹配，
  ##    只会命中第一个种子的块 ⇒ 其它种子全空、sapply 退化成 list（我踩过）
  W <- dcast(d, barcode ~ seed, value.var = "domain")
  bc <- W$barcode
  L  <- as.matrix(W[, -1L, with = FALSE])
  seeds <- suppressWarnings(as.integer(sub("^V?", "", colnames(L))))
  ## 共分配矩阵：同域计 1，跨域计 0，逐种子平均
  C <- matrix(0, length(bc), length(bc), dimnames = list(bc, bc))
  for (j in seq_along(seeds)) C <- C + outer(L[, j], L[, j], "==")
  C <- C / length(seeds)
  K <- as.integer(median(sapply(seeds, function(s) uniqueN(L[, s]))))
  hc <- hclust(as.dist(1 - C), method = "average")
  dom <- cutree(hc, k = K)
  ## 共识强度：每个 spot 与"同域其他 spot"的平均共分配率
  agree <- vapply(seq_along(dom), function(i) {
    same <- dom == dom[i]; same[i] <- FALSE
    if (!any(same)) NA_real_ else mean(C[i, same])
  }, numeric(1))

  ## ---- 2. 域内容：组成（rctd_a 6 类 / rctd_d 39 亚型）----
  A <- read_weights("rctd_a", slide); D <- read_weights("rctd_d", slide)
  if (is.null(A)) return(NULL)
  j <- match(bc, rownames(A)); if (anyNA(j)) { step("  ! %s 组成缺 spot", slide); return(NULL) }
  A <- A[j, , drop = FALSE]; rownames(A) <- bc
  comp <- do.call(rbind, lapply(sort(unique(dom)), function(k) colMeans(A[dom == k, , drop = FALSE])))
  comp <- as.data.table(comp); comp[, domain := sort(unique(dom))]
  comp[, n_spot := as.integer(table(dom))]
  comp[, lineage := names(.SD)[max.col(.SD, ties.method = "first")], .SDcols = setdiff(names(comp), c("domain","n_spot","lineage"))]
  comp[, purity := do.call(pmax, c(as.list(.SD), list(na.rm = TRUE))), .SDcols = setdiff(names(comp), c("domain","n_spot","lineage","purity"))]

  fine <- NULL
  if (!is.null(D)) {
    jj <- match(bc, rownames(D))
    if (!anyNA(jj)) {
      D <- D[jj, , drop = FALSE]; rownames(D) <- bc
      fm <- do.call(rbind, lapply(sort(unique(dom)), function(k) colMeans(D[dom == k, , drop = FALSE])))
      top3 <- apply(fm, 1, function(v) paste(sprintf("%s:%.2f", names(sort(v, decreasing = TRUE))[1:3], sort(v, decreasing = TRUE)[1:3]), collapse = ";"))
      fine <- data.table(domain = sort(unique(dom)), fine_top3 = top3)
    }
  }

  ## ---- 3. marker：profiles 的域均值，域内跨域 z 分数 ----
  mk <- NULL
  pf <- file.path(BN, slide, "profiles_agfT.tsv.gz")
  if (file.exists(pf)) {
    p <- fread(pf, sep = "\t")
    p <- p[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution]
    p <- p[domain %in% sort(unique(dom))]
    if (nrow(p) >= 3L) {
      g <- setdiff(names(p), c("k_geom","lambda","resolution","domain"))
      X <- as.matrix(p[, ..g]); z <- scale(X)
      mk <- data.table(domain = p$domain,
                       mark_top5 = apply(z, 1, function(v) paste(g[order(-v)][1:5], collapse = ",")))
    }
  }

  comp <- merge(comp, fine, by = "domain", all.x = TRUE)
  if (!is.null(mk)) comp <- merge(comp, mk, by = "domain", all.x = TRUE)
  comp[, `:=`(slide = slide, K_consensus = K,
              stage = sub("-1$", "", sub("^GSM[0-9]+_P[0-9]+_", "", slide)),
              agree_med = round(median(agree, na.rm = TRUE), 4))]
  comp[, .(slide, stage, K_consensus, domain, n_spot, lineage, purity = round(purity, 4),
           agree_med, fine_top3, mark_top5)]
}

RES <- rbindlist(lapply(slides, function(s) {
  r <- tryCatch(one_slide(s), error = function(e) { step("  !! %s: %s", s, conditionMessage(e)); NULL })
  if (!is.null(r)) step("  OK %s（%d 域）", s, nrow(r))
  r
}), fill = TRUE)

if (is.null(RES) || !nrow(RES)) stop("无结果", call. = FALSE)
fwrite(RES, file.path(OUT, "d6_domain_content.tsv"), sep = "\t")
step("落盘 d6_domain_content.tsv：%d 张 × 共 %d 个域", uniqueN(RES$slide), nrow(RES))

## ---- 汇报：域类型 × 期别 / 纯度 / 共识强度 ----
con <- file(file.path(OUT, "d6_summary.txt"), "wt")
w <- function(...) cat(sprintf(...), file = con, append = TRUE)
w("=== §16 域内容剖析（主档 agfT k=18 l=0.2 r=0.5；5 种子共识划分）===\n")
w("%s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("域总数 %d；每张切片域数（共识 K）中位 %g\n\n", nrow(RES), median(RES[, .N, by = slide]$N))

w("[1] 域的主导谱系 × 期别（域个数）\n")
tab <- dcast(RES[, .N, by = .(stage, lineage)], stage ~ lineage, value.var = "N", fill = 0)
print(tab, file = con)
w("\n[1b] 同一张表化成比例（行归一）\n")
pr <- RES[, .N, by = .(stage, lineage)][, frac := N / sum(N), by = stage]
print(dcast(pr, stage ~ lineage, value.var = "frac", fill = 0), file = con, digits = 3)

w("\n[2] 域纯度（主导谱系占比）按期别\n")
print(RES[, .(n = .N, purity_med = round(median(purity), 3),
              purity_q25 = round(quantile(purity, .25), 3),
              pure_ge70 = round(mean(purity >= 0.70), 2)), by = stage]
      [order(purity_med)], file = con)
w("\n[3] 共识强度（spot 与同域伙伴的平均共分配率）\n")
print(RES[, .(n = .N, agree_med = round(median(agree_med), 3)), by = stage]
      [order(agree_med)], file = con)
close(con)
cat(readLines(file.path(OUT, "d6_summary.txt")), sep = "\n")
step("产物 → %s", OUT)
