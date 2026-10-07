#!/usr/bin/env Rscript
# 23_spatial_subtype_fluorescence.R —— RCTD 全亚群（39 型）空转「荧光图」（只读）
#
# 🔴 命名禁令：RCTD 权重**绝不当恶性判读**。本图只呈现**组成比例**，
#    不产生任何「恶性度 / 肿瘤 / malignancy」读法。
#
# 口径来源
#   results/08_spatial_deconv/rctd_d/per_slide/<slide>.weights.tsv.gz  （39 亚型权重）
#   data/visium_spatial/<slide>/spatial/tissue_positions.csv           （坐标）
#   亚型→谱系对照 = 项目**自己的**亚聚类注释（不手编）：
#     内皮 / 髓系  走 S4 面板  `*_s4_cluster_annotation.csv`（9 + 15 型）
#     其余四谱系   走经典面板  `*_cluster_annotation.csv`
#
# 画法：黑底 + 加性混合（每个 spot 的颜色 = Σ_k w_k · color_k）。
#       这是**组成**的直观呈现，不是任何形式的判定。

RES <- "/home/eto/luad_v2/results/08_spatial_deconv"
DAT <- "/home/eto/luad_v2/data/visium_spatial"
ANN <- "/home/eto/luad_v2/results/05_annotation"
OUT <- file.path(RES, "spatial_cnv/cohort_borrow", "subtype_fluorescence")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## 默认跑 rctd_d 里**全部**切片；也可以把切片名当参数传进来只跑指定的
WDIR   <- file.path(RES, "rctd_d/per_slide")
SLIDES <- sort(sub("\\.weights\\.tsv\\.gz$", "",
                   list.files(WDIR, pattern = "\\.weights\\.tsv\\.gz$")))
.a <- commandArgs(trailingOnly = TRUE)
SHEET_ONLY <- "--sheet-only" %in% .a          # 调试用：只出全队列一览
.a <- setdiff(.a, "--sheet-only")
if (length(.a)) SLIDES <- intersect(SLIDES, .a)
stopifnot(length(SLIDES) > 0)
cat(sprintf("共 %d 张切片\n", length(SLIDES)))

## ——— 1. 亚型 → 谱系（取自项目自己的注释表）———
LIN_OF_FILE <- c(epiA = "上皮", fibroA = "成纤维", tnkA = "T/NK",
                 myeloidA = "髓系", endoA = "内皮", bplasmaA = "B/浆")
S4_FILES    <- c("endoA", "myeloidA")          # S4 面板只覆盖这两个谱系

sub2lin <- list()
read_argmax <- function(f) {
  if (!file.exists(f)) return(character(0))
  d <- read.csv(f, stringsAsFactors = FALSE)
  unique(d$argmax[!is.na(d$argmax)])
}
for (k in names(LIN_OF_FILE)) {
  if (k %in% S4_FILES) pat <- sprintf("%s_s4_cluster_annotation.csv", k)
  else                 pat <- sprintf("%s_cluster_annotation.csv", k)
  a <- read_argmax(file.path(ANN, pat))
  cat(sprintf("  %-9s %-34s → %d 型\n", k, pat, length(a)))
  for (x in a) sub2lin[[x]] <- LIN_OF_FILE[[k]]
}
## RCTD 输出把带 `/` 的名字改写成了 `_`，用**它自己的对照表**归一化（别自己猜）
nm <- read.delim(file.path(RES, "rctd_d/cell_type_name_map.tsv"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(nm))) names(sub2lin)[names(sub2lin) == nm$original[i]] <- nm$adapted[i]

## 参考侧 39 名以**权重文件自己的列顺序**为准（比读 manifest 更直接）
REF <- colnames(read.delim(
  gzfile(file.path(RES, "rctd_d/per_slide", paste0(SLIDES[1], ".weights.tsv.gz"))),
  row.names = 1, check.names = FALSE))
missing <- setdiff(REF, names(sub2lin))
if (length(missing)) stop("亚型→谱系 缺：", paste(missing, collapse = ", "))
cat(sprintf("对照齐：%d / %d 个亚型\n", length(REF), length(REF)))

LIN <- c("上皮", "成纤维", "髓系", "内皮", "T/NK", "B/浆")
## 🔴 本机无 CJK 字体 —— 图内一切文字必须是英文（否则渲染成方块）
LIN_EN <- c(上皮 = "Epithelial", 成纤维 = "Fibroblast", 髓系 = "Myeloid",
            内皮 = "Endothelial", "T/NK" = "T/NK", "B/浆" = "B/Plasma")
sub_lin <- unlist(sub2lin)[REF]
stopifnot(length(REF) == 39, setequal(unique(sub_lin), LIN))

## ——— 2. 颜色：谱系占色相弧，弧内按亚型铺开 ———
HUE <- c(上皮 = 12, 成纤维 = 52, 髓系 = 125, 内皮 = 182, "T/NK" = 235, "B/浆" = 302)
SPAN <- c(上皮 = 26, 成纤维 = 16, 髓系 = 62, 内皮 = 26, "T/NK" = 30, "B/浆" = 26)

sub_color <- character(length(REF)); names(sub_color) <- REF
for (L in LIN) {
  idx <- which(sub_lin == L); n <- length(idx)
  off <- if (n == 1) 0 else seq(-SPAN[[L]] / 2, SPAN[[L]] / 2, length.out = n)
  ## 明度也动一下，让同弧相邻亚型分得开
  v <- if (n == 1) 1 else seq(0.85, 1.0, length.out = n)[order(off)]
  sub_color[idx] <- hsv((HUE[[L]] + off) %% 360 / 360, s = 0.95, v = v)
}
lin_color <- setNames(hsv(HUE / 360, s = 0.95, v = 1), names(HUE))   # hsv() 不保名，必须显式 setNames

## ——— 3. 组装 ——
## 每张切片只读一次权重+坐标（39 张子图要复用，别读 39 遍）
load_slide <- function(S) {
  w <- read.delim(gzfile(file.path(RES, "rctd_d/per_slide", paste0(S, ".weights.tsv.gz"))),
                  row.names = 1, check.names = FALSE)
  pos <- read.csv(file.path(DAT, S, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
  m <- match(rownames(w), pos$barcode)
  if (anyNA(m)) stop("坐标对不上：", S)
  W <- as.matrix(w[, REF]); rownames(W) <- NULL
  list(W = W, ar = pos$array_row[m], ac = pos$array_col[m],
       R = max(pos$array_row) + 1L, C = max(pos$array_col) + 1L)
}
## RGB 向量 → array 栅格（行 = array_row，列 = array_col；array_col 双倍存储故宽 224）
raster_of <- function(sl, RGB) {
  RGB <- pmin(pmax(RGB, 0), 1) ^ 0.75                      # gamma：提亮暗部，荧光感
  M <- matrix("#000000", sl$R, sl$C)
  M[cbind(sl$ar + 1L, sl$ac + 1L)] <- rgb(RGB[, 1], RGB[, 2], RGB[, 3])
  M
}
## mode: composite = 39 型加性混合；lineage = 单谱系求和；subtype = 单个亚型
build_rgb <- function(sl, mode = c("composite", "lineage", "subtype"), key = NULL) {
  mode <- match.arg(mode)
  RGB <- if (mode == "composite") {
    sl$W %*% t(grDevices::col2rgb(sub_color[REF]) / 255)          # Σ w_k · color_k
  } else if (mode == "lineage") {
    rowSums(sl$W[, sub_lin == key, drop = FALSE]) %*%
      t(grDevices::col2rgb(lin_color[[key]]) / 255)
  } else {
    sl$W[, key] %*% t(grDevices::col2rgb(sub_color[[key]]) / 255)
  }
  raster_of(sl, RGB)
}

## ——— 补洞：Visium 六角栅格奇数格是空的，用邻格均值填上 ⇒ 连续荧光观感 ———
shiftm <- function(A, dr, dc) {
  R <- nrow(A); C <- ncol(A); out <- matrix(0, R, C)
  rs <- max(1, 1 - dr):min(R, R - dr); cs <- max(1, 1 - dc):min(C, C - dc)
  out[rs + dr, cs + dc] <- A[rs, cs]; out
}
fill_gaps <- function(M, iters = 2) {
  R <- nrow(M); C <- ncol(M)
  occ <- matrix(M != "#000000", R, C)
  if (all(occ)) return(M)
  rgbm <- col2rgb(M) / 255
  A <- lapply(1:3, function(k) matrix(rgbm[k, ], R, C))      # 三个通道各一张矩阵
  NB <- list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1),
             c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  for (it in seq_len(iters)) {
    S <- lapply(1:3, function(k) matrix(0, R, C)); N <- matrix(0, R, C)
    for (d in NB) {
      for (k in 1:3) S[[k]] <- S[[k]] + shiftm(A[[k]], d[1], d[2])
      N <- N + shiftm(occ * 1, d[1], d[2])
    }
    ok <- (!occ) & (N > 0)
    for (k in 1:3) { m <- S[[k]] / pmax(N, 1); a <- A[[k]]; a[ok] <- m[ok]; A[[k]] <- a }
    occ[ok] <- TRUE
    if (all(occ)) break
  }
  matrix(rgb(A[[1]], A[[2]], A[[3]]), R, C)
}

## ——— 几何：`array_col` 是双倍存储，且**行距 ≠ 列距** ———
## 🔴 别当方格画。实测 每 1 个 array_row = 343.8 px、每 2 个 array_col = 397.5 px
##    ⇒ 行距只有列距的 **0.865**。按方格画会把组织**纵向拉长 15.6%**（第一版就是这么错的）。
##    这里从像素坐标**现算**每张切片的比值，不写死常数。
##    朝向同理：`array_row/col → pxl_row/col` 是**轴对齐**的（全队列 56 张 |cor| ≥ 0.9996），
##    但**每个轴各自的符号不一样**：全队列 36 张 (+col, −row)、17 张 (−col, +row)、3 张 (+, +)。
##    ⇒ 必须**逐轴**判，不能整幅转 180°（那样只有 (−,−) 才对，实测一张都没有）。
geom_of <- function(S) {
  p <- read.csv(file.path(DAT, S, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
  p <- p[p$in_tissue == 1, ]
  ## 用跨度比值（取绝对值）而不是回归斜率
  spr <- diff(range(p$pxl_row_in_fullres)) / diff(range(p$array_row))
  spc <- diff(range(p$pxl_col_in_fullres)) / diff(range(p$array_col))
  list(yasp   = spr / (2 * spc),
       flip_c = cor(p$array_col, p$pxl_col_in_fullres) < 0,
       flip_r = cor(p$array_row, p$pxl_row_in_fullres) < 0)
}

## 🔴 `array_col` 双倍存储 ⇒ 横向除 2；纵向再乘 yasp 才是真实比例
draw <- function(M, main, yasp, smooth = TRUE, cex.main = 0.95, col.main = "grey80") {
  if (smooth) M <- fill_gaps(M)
  nr <- nrow(M); nc <- ncol(M) / 2
  plot.new(); plot.window(xlim = c(0, nc), ylim = c(nr * yasp, 0), asp = 1)
  rasterImage(as.raster(M), 0, nr * yasp, nc, 0, interpolate = TRUE)
  title(main = main, col.main = col.main, cex.main = cex.main)
}
wrap <- function(x, w) paste(strwrap(x, width = w), collapse = "\n")

## ——— 4. 出图 ———
geo_rows <- list(); sheet <- list()
for (S in SLIDES) {
  cat(sprintf("画 %s …\n", S))
  st <- sub("-\\d+$", "", sub("^GSM[0-9]+_P[0-9]+_", "", S))
  g  <- geom_of(S)
  cat(sprintf("  %s  y/2col = %.4f  flip_col=%s flip_row=%s\n", S, g$yasp, g$flip_c, g$flip_r))
  ## 把 array 栅格摆到**图像（H&E）朝向**：逐轴翻，坐标本身一个数都不改
  to_img <- function(M) {
    if (g$flip_r) M <- M[nrow(M):1, , drop = FALSE]
    if (g$flip_c) M <- M[, ncol(M):1, drop = FALSE]
    M
  }
  sl <- load_slide(S)
  geo_rows[[S]] <- data.frame(slide = S, n_spot = nrow(sl$W), yasp = g$yasp,
                              flip_col = g$flip_c, flip_row = g$flip_r, row.names = NULL)
  Mc <- to_img(build_rgb(sl, "composite"))
  ## 每 2 格抽 1，留着拼全队列一览
  ## 🔴 别按奇数下标抽样：`to_img` 翻过轴之后六角栅格的 (row+col) 奇偶**跟着翻**，
  ##    抽样会正好全落空（实测 53/56 张变成全黑）。直接留整幅，56 张也就 ~13 MB。
  sheet[[S]] <- Mc

  if (SHEET_ONLY) next

  ## (a) 合并图（hero）
  png(file.path(OUT, sprintf("%s_subtype_composite.png", S)),
      width = 1180, height = 1320, res = 150, bg = "black")
  par(mar = c(0.5, 0.5, 2.6, 0.5), bg = "black")
  draw(Mc, sprintf("RCTD 39-subtype composite (additive)  |  %s  (%s)", S, st), g$yasp)
  legend("topright", legend = unname(LIN_EN[LIN]), col = unname(lin_color[LIN]),
         pch = 15, pt.cex = 1.6, text.col = "grey85", bty = "n", cex = 0.8,
         title = "lineage (hue anchor)", title.col = "grey85")
  dev.off()

  ## (b) 逐谱系通道（6 格）
  png(file.path(OUT, sprintf("%s_lineage_channels.png", S)),
      width = 1600, height = 1300, res = 150, bg = "black")
  par(mfrow = c(2, 3), mar = c(0.5, 0.5, 5.0, 0.5), bg = "black", oma = c(0, 0, 3.0, 0))
  for (L in LIN) {
    sub <- REF[sub_lin == L]
    draw(to_img(build_rgb(sl, "lineage", L)),
         sprintf("%s (%d)\n%s", LIN_EN[[L]], length(sub), wrap(paste(sub, collapse = ", "), 78)),
         yasp = g$yasp, smooth = FALSE, cex.main = 0.60)
  }
  mtext(sprintf("RCTD subtype fraction by lineage channel  |  %s  (%s)", S, st),
        side = 3, line = 1.4, outer = TRUE, col = "grey85", cex = 0.95)
  dev.off()

  ## (c) 逐亚型 39 格（精细亚群；每格一个亚型，用它自己的颜色）
  png(file.path(OUT, sprintf("%s_subtype_panels.png", S)),
      width = 2000, height = 1850, res = 150, bg = "black")
  par(mfrow = c(6, 7), mar = c(0.3, 0.3, 2.0, 0.3), bg = "black", oma = c(0, 0, 2.8, 0))
  for (k in REF) {
    draw(to_img(build_rgb(sl, "subtype", k)), k, yasp = g$yasp,
         smooth = TRUE, cex.main = 0.52, col.main = sub_color[[k]])
  }
  mtext(sprintf("RCTD per-subtype fraction (39 subtypes, one panel each)  |  %s  (%s)", S, st),
        side = 3, line = 1.0, outer = TRUE, col = "grey85", cex = 0.95)
  dev.off()
}

## ——— 5. 对照表落盘（可审计）———
write.table(data.frame(subtype = REF, lineage = unname(sub_lin), color = unname(sub_color[REF])),
            file.path(OUT, "subtype_lineage_color_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, geo_rows), file.path(OUT, "slide_geometry_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ——— 6. 全队列一览：56 张合并图按期别拼成一张 ———
meta <- do.call(rbind, lapply(SLIDES, function(S) {
  m <- regmatches(S, regexec("^(GSM[0-9]+)_(P[0-9]+)_(.+)$", S))[[1]]
  if (length(m) < 4) return(NULL)
  data.frame(slide = S, patient = m[3], tag = m[4], stage = sub("-\\d+$", "", m[4]),
             row.names = NULL)
}))
meta$ord  <- match(meta$stage, c("Normal", "AAH", "AIS", "MIA", "LUAD"))
meta$pnum <- as.integer(sub("^P", "", meta$patient))
meta <- meta[order(meta$ord, meta$pnum, meta$tag), ]
NCOL <- 8; NROW <- ceiling(nrow(meta) / NCOL)
png(file.path(OUT, "cohort_overview_composite.png"),
    width = 2600, height = 2600, res = 150, bg = "black")
par(mfrow = c(NROW, NCOL), mar = c(0.4, 0.4, 1.9, 0.4), bg = "black",
    oma = c(0, 0, 3.0, 0))
for (i in seq_len(nrow(meta))) {
  S <- meta$slide[i]; M <- fill_gaps(sheet[[S]])
  nr <- nrow(M); nc <- ncol(M) / 2; ya <- geo_rows[[S]]$yasp
  plot.new(); plot.window(xlim = c(0, nc), ylim = c(nr * ya, 0), asp = 1)
  rasterImage(as.raster(M), 0, nr * ya, nc, 0, interpolate = TRUE)
  title(main = sprintf("%s %s", meta$patient[i], meta$tag[i]),
        col.main = "grey85", cex.main = 0.75)
}
mtext("P21-scale cohort: RCTD 39-subtype composite (additive)  |  ordered by stage: Normal, AAH, AIS, MIA, LUAD",
      side = 3, line = 1.2, outer = TRUE, col = "grey85", cex = 1.0)
dev.off()

cat(sprintf("\n产物目录：%s\n", OUT))
