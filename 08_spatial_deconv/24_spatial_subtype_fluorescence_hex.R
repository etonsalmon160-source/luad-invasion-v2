#!/usr/bin/env Rscript
# 24_spatial_subtype_fluorescence_hex.R
#   RCTD 39 亚型 / 6 谱系的**逐 spot 六边形**荧光图。
#
# 相对 23_ 改了三件事（23_ 保留不删，可随时回去对照）
#   ① 画法：23_ 底层只有 ~224x128 个格子，却拉宽到 1100+ 像素，再叠 fill_gaps()
#      邻居平均 + rasterImage(interpolate=TRUE) 两层平滑 ⇒ 糊；而且 fill_gaps 会
#      **给被 QC 剔除的 spot 也涂上颜色**，等于在没有数据的地方显示数据。
#      本版每个 spot 画一个正好铺满六角栅格的正六边形，不插值、不填洞 ⇒
#      未被分析的 spot 是**黑洞**（真实情况）。
#   ② 配色：39 色按**谱系内细胞数降序**排成 6 条渐变（梯度轴来自 05_annotation 的
#      cluster 表，表里 cluster 本身按细胞数排）。默认 **A = 鲜艳版**：色相沿谱系弧段
#      铺开，s=0.95，v 从 1.00 递减到 0.78。加 `--equal-luminance` 切到 **C = 等亮度版**
#      （39 色 + 6 谱系色相对亮度全配平到 Y0=0.42，粉彩）。
#   ③ 权重：合成前把每个 spot 的权重**归一到和为 1**（原权重行和中位 0.917）。
#
# 等亮度版的来龙去脉（保留备查；**当前默认不用它**）
#   P10_LUAD 上，A 版 spot 亮度与测序深度的 Spearman 达 −0.76 ~ −0.80。
#   来源不是"密度进了公式"（RCTD 权重是比例，深度不在式子里），而是：深 spot 的
#   组成更集中（有效亚型数 vs 深度 −0.58）→ 集中到上皮 → 而上皮是红色，红色天然
#   亮度 0.213（绿色 0.715）⇒ "深浅"其实在替你显示"上皮占了多少"。
#   C 版把 39 色亮度配平后，合成图亮度 = Σ(占比 × 0.42) = 常数（实测极差 0.99/255）
#   ⇒ 图上在数学上不存在明度轴；代价是红系被提亮成粉彩，观感发灰。
#   ⚠️ 用 A 版就必须知道：**明度不是纯装饰**，它与"上皮占比"同向，因而局部与测序深度
#      弱相关。要读"某类多不多"请看单通道图 / 分格图，别读合成图的深浅。
#
# 图面约定
#   合并图、全队列一览 —— **纯图，图内零文字**；标注单独出：
#     colour_key_gradient_bars.png（每个谱系一条色阶条，按顺序标亚型名 + 该谱系内细胞数）
#     colour_key_swatch_grid.png （39 个色块 6 行排开）
#     cohort_overview_order.tsv （全队列一览 56 格的排布顺序）
#   6 谱系通道图、39 亚型分格图 —— 保留每格名称（否则对不上是哪一格）。
#
# 🔴 命名禁令：RCTD 权重**绝不当恶性判读**。本图只呈现**组成比例**，
#    不产生任何「恶性度 / 肿瘤 / malignancy」读法。
# 🔴 本机无 CJK 字体 —— 图内一切文字必须是英文。

RES <- "/home/eto/luad_v2/results/08_spatial_deconv"
DAT <- "/home/eto/luad_v2/data/visium_spatial"
ANN <- "/home/eto/luad_v2/results/05_annotation"
OUT <- file.path(RES, "spatial_cnv/cohort_borrow/subtype_fluorescence")
WDIR <- file.path(RES, "rctd_d/per_slide")
stopifnot(dir.exists(OUT))

LIN <- c("上皮", "成纤维", "髓系", "内皮", "T/NK", "B/浆")
LIN_EN <- c(上皮 = "Epithelial", 成纤维 = "Fibroblast", 髓系 = "Myeloid",
            内皮 = "Endothelial", "T/NK" = "T/NK", "B/浆" = "B/Plasma")
HUE <- c(上皮 = 12, 成纤维 = 52, 髓系 = 125, 内皮 = 182, "T/NK" = 235, "B/浆" = 302)
## 每个谱系在色轮上占的弧段（度）；相邻谱系之间留空，B/浆 只有 2 型故额外收窄
ARC <- list(上皮 = c(341, 388), 成纤维 = c(36, 78), 髓系 = c(86, 149.5),
            内皮 = c(162, 204.5), "T/NK" = c(212.5, 264.5), "B/浆" = c(285, 320))
Y0 <- 0.42                      # 等亮度目标（相对亮度，Rec.709）

SLIDES <- sort(sub("\\.weights\\.tsv\\.gz$", "",
                   list.files(WDIR, pattern = "\\.weights\\.tsv\\.gz$")))
.a <- commandArgs(trailingOnly = TRUE)
SHEET_ONLY <- "--sheet-only" %in% .a
EQ_LUM     <- "--equal-luminance" %in% .a
.a <- setdiff(.a, c("--sheet-only", "--equal-luminance"))
if (length(.a)) SLIDES <- intersect(SLIDES, .a)
stopifnot(length(SLIDES) > 0)
PALETTE <- if (EQ_LUM) "C" else "A"
cat(sprintf("共 %d 张切片\n", length(SLIDES)))

## ——— 1. 梯度轴：谱系内细胞数降序（只读注释表，不改任何数）———
AFILES <- c(上皮 = "epiA_cluster_annotation.csv", 成纤维 = "fibroA_cluster_annotation.csv",
            髓系 = "myeloidA_s4_cluster_annotation.csv", 内皮 = "endoA_s4_cluster_annotation.csv",
            "T/NK" = "tnkA_cluster_annotation.csv", "B/浆" = "bplasmaA_cluster_annotation.csv")
abund_tab <- do.call(rbind, lapply(LIN, function(L) {
  a <- read.csv(file.path(ANN, AFILES[[L]]), stringsAsFactors = FALSE)
  a <- a[!(toupper(as.character(a$excluded)) %in% c("TRUE", "T")), ]   # 只数进了注释的簇
  data.frame(subtype = gsub("/", "_", a$argmax), lineage = L,
             n_cells = a$n_cells, stringsAsFactors = FALSE)
}))
abund <- tapply(abund_tab$n_cells, abund_tab$subtype, sum)
REF <- unname(unlist(lapply(LIN, function(L) {
  s <- names(abund)[abund_tab$lineage[match(names(abund), abund_tab$subtype)] == L]
  s[order(abund[s], decreasing = TRUE)]
})))
stopifnot(length(REF) == 39, !anyDuplicated(REF))
sub_lin <- setNames(unname(abund_tab$lineage[match(REF, abund_tab$subtype)]), REF)
stopifnot(all(sub_lin %in% LIN), !anyNA(sub_lin))
cat("梯度轴 = 谱系内细胞数降序（来自 05_annotation 的 cluster 表）\n")

## ——— 2. 配色：39 亚型 + 6 谱系，全部等相对亮度 ———
RELY <- c(0.2126, 0.7152, 0.0722)
relY <- function(v) sum(RELY * v)
hsv_rgb <- function(h, s, v) grDevices::col2rgb(hsv(h / 360, s, v))[, 1] / 255
## hsv() 落盘是 8 位十六进制，颜色本身带取整误差 ⇒ 最后要按**取整后**的亮度再微调一步，
## 否则 39 个色各差 ±0.5/255，合成图的亮度又会随组成轻微起伏。
eqY <- function(hue, Y0 = 0.42) {
  hue <- hue %% 360
  y1 <- relY(grDevices::col2rgb(hsv(hue / 360, 1, 1))[, 1] / 255)
  dh <- seq(-4, 4, by = 0.5)
  sc <- if (y1 >= Y0) {                                   # 纯色够亮 ⇒ 压暗
    v <- Y0 / y1; seq(max(1e-6, v * 0.97), min(1, v * 1.03), length.out = 121)
  } else {                                                # 纯色不够亮 ⇒ 加白
    s <- uniroot(function(s) relY(hsv_rgb(hue, s, 1)) - Y0, c(0, 1), tol = 1e-7)$root
    seq(max(1e-6, s - 0.02), min(1, s + 0.02), length.out = 121)
  }
  H <- rep((hue + dh) %% 360, each = length(sc))
  cand <- if (y1 >= Y0) hsv(H / 360, 1, rep(sc, times = length(dh)))
          else          hsv(H / 360, rep(sc, times = length(dh)), 1)
  cand[[which.min(abs(as.vector(RELY %*% (grDevices::col2rgb(cand) / 255)) - Y0))]]
}
sub_color <- character(0)
for (L in LIN) {
  s <- REF[sub_lin[REF] == L]; n <- length(s)
  t <- if (n == 1) 0 else seq(0, 1, length.out = n)        # 0 = 该谱系最常见
  a <- ARC[[L]]
  sub_color[s] <- if (PALETTE == "A")
    hsv(((a[1] + (a[2] - a[1]) * t) %% 360) / 360, s = 0.95, v = 1.00 - 0.22 * t)
  else
    vapply(a[1] + (a[2] - a[1]) * t, eqY, "", Y0 = Y0)
}
lin_color <- setNames(if (PALETTE == "A") hsv(HUE / 360, s = 0.95, v = 1)
                      else vapply(HUE, eqY, "", Y0 = Y0), names(HUE))
ycol <- function(x) apply(grDevices::col2rgb(x), 2, relY)
if (PALETTE == "A") {
  cat(sprintf("配色 = A 鲜艳版：39 亚型相对亮度 %.3f~%.3f，极差 %.3f ⇒ **明度轴存在**\n",
              min(ycol(sub_color[REF])), max(ycol(sub_color[REF])),
              diff(range(ycol(sub_color[REF])))))
} else {
  cat(sprintf("配色 = C 等亮度版（Y0=%.2f）：39 亚型 极差 %.4f，6 谱系 极差 %.4f\n",
              Y0, diff(range(ycol(sub_color[REF]))), diff(range(ycol(lin_color)))))
}

## ——— 3. 读取切片：权重 + 坐标 ———
load_slide <- function(S) {
  w <- read.delim(gzfile(file.path(WDIR, paste0(S, ".weights.tsv.gz"))),
                  row.names = 1, check.names = FALSE)
  keep <- setdiff(colnames(w), "spot_class")
  stopifnot(all(REF %in% keep))
  ## 只作显示：归一到和为 1（原行和中位 0.917），见文件头 ③
  W <- as.matrix(w[, REF]); W <- W / rowSums(W); rownames(W) <- NULL
  p <- read.csv(file.path(DAT, S, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
  m <- match(rownames(w), p$barcode)
  if (anyNA(m)) stop("坐标对不上：", S)
  p1 <- p[p$in_tissue == 1, ]
  ## 行距 ≠ 列距：array_col 是双倍存储 ⇒ 横向除 2、纵向乘 yasp 才是真实比例
  yasp <- (diff(range(p1$pxl_row_in_fullres)) / diff(range(p1$array_row))) /
          (2 * diff(range(p1$pxl_col_in_fullres)) / diff(range(p1$array_col)))
  flip_c <- cor(p1$array_col, p1$pxl_col_in_fullres) < 0
  flip_r <- cor(p1$array_row, p1$pxl_row_in_fullres) < 0
  x <- p$array_col[m] / 2
  y <- p$array_row[m] * yasp
  ## 摆到**图像（H&E）朝向**：逐轴翻，坐标本身一个数都不改
  if (flip_c) x <- (max(p$array_col) / 2) - x
  if (flip_r) y <- (max(p$array_row) * yasp) - y
  list(W = W, x = x, y = y, yasp = yasp, n = nrow(w),
       flip_c = flip_c, flip_r = flip_r)
}

## ——— 4. 颜色：加法混色，**不加 γ**（γ 是非线性的，会破坏等亮度）———
cols_of <- function(sl, mode = c("composite", "lineage", "subtype"), key = NULL) {
  mode <- match.arg(mode)
  RGB <- if (mode == "composite") {
    sl$W %*% t(grDevices::col2rgb(sub_color[REF]) / 255)          # Σ w̃_k · color_k
  } else if (mode == "lineage") {
    rowSums(sl$W[, sub_lin[REF] == key, drop = FALSE]) %*%
      t(grDevices::col2rgb(lin_color[[key]]) / 255)
  } else {
    sl$W[, key] %*% t(grDevices::col2rgb(sub_color[[key]]) / 255)
  }
  rgb(pmin(pmax(RGB, 0), 1))
}

## ——— 5. 画法：逐 spot 正六边形（最近邻距恒为 1 ⇒ 外接半径 1/√3 正好铺满）———
HEX_R  <- 1 / sqrt(3)
HEX_TH <- (30 + 60 * (0:5)) * pi / 180
HEX_VX <- HEX_R * cos(HEX_TH)
HEX_VY <- HEX_R * sin(HEX_TH)

draw_spots <- function(x, y, cols, yasp, main = "", cex.main = 0.55,
                       col.main = "grey80", pad = 0.9, mar = c(0.3, 0.3, 0.3, 0.3)) {
  par(mar = mar)
  plot.new()
  plot.window(xlim = range(x) + c(-pad, pad),
              ylim = range(y) + c(-pad, pad) * yasp, asp = 1)
  n <- length(x)
  X <- matrix(NA_real_, n, 7L); Y <- matrix(NA_real_, n, 7L)   # 第 7 列留 NA 分隔多边形
  X[, 1:6] <- outer(x, HEX_VX, "+")
  Y[, 1:6] <- outer(y, HEX_VY, "+")
  ## 一次 polygon 调用画全部 spot：NA 分隔 ⇒ 每行一个六边形，col 逐多边形回收
  polygon(as.vector(t(X)), as.vector(t(Y)), col = cols, border = cols, lwd = 0.25)
  if (nzchar(main)) title(main = main, col.main = col.main, cex.main = cex.main)
}
wrap <- function(x, w) paste(strwrap(x, width = w), collapse = "\n")

## ——— 6. 单独出的标注（纯图之外的唯一说明来源）———
key_bars <- function() {
  png(file.path(OUT, "colour_key_gradient_bars.png"),
      width = 2600, height = 3000, res = 150, bg = "black")
  par(mar = c(0.4, 0.4, 0.4, 0.4), bg = "black")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(10.4, 0.1))
  text(0.02, 0.30, "RCTD 39-subtype colour key (equal-luminance)", adj = 0,
       col = "grey90", cex = 1.6)
  text(0.02, 0.62, "order inside each lineage = descending cell number (05_annotation cluster tables)  ·  all 39 colours share the same relative luminance, so brightness carries no quantity",
       adj = 0, col = "grey55", cex = 0.60)
  x0 <- 0.20; x1 <- 0.90
  for (i in seq_along(LIN)) {
    L <- LIN[i]; s <- REF[sub_lin[REF] == L]; n <- length(s)
    y <- 1.1 + 1.5 * (i - 1)
    rect(x0, y - 0.115, x1, y + 0.115, col = "grey15", border = "grey30")
    xs <- seq(x0, x1, length.out = n + 1)
    for (k in seq_len(n)) rect(xs[k], y - 0.105, xs[k + 1], y + 0.105,
                               col = sub_color[s[k]], border = "grey25", lwd = 0.4)
    text(x0 - 0.012, y - 0.05, sprintf("%s  (%d)", LIN_EN[[L]], n), adj = 1,
         col = "grey90", cex = 0.95)
    text(x0 - 0.012, y + 0.10, "high to low abundance", adj = 1, col = "grey45", cex = 0.5)
    for (k in seq_len(n))
      text((xs[k] + xs[k + 1]) / 2, y + 0.16, s[k], srt = 90, adj = 1,
           col = sub_color[s[k]], cex = 0.60, xpd = NA)
    text(x1 + 0.010, y, formatC(abund[s], format = "d", big.mark = ","),
         adj = 0, col = "grey55", cex = 0.52, xpd = NA)
  }
  invisible(dev.off())
}
key_grid <- function() {
  NCOL <- 12
  png(file.path(OUT, "colour_key_swatch_grid.png"),
      width = 2600, height = 2000, res = 150, bg = "black")
  par(mar = c(0.5, 0.5, 3.0, 0.5), bg = "black")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(6.6, 0.4))
  text(0.02, 0.55, "RCTD 39-subtype swatch key (equal-luminance)", adj = 0,
       col = "grey90", cex = 1.5)
  cw <- 1 / NCOL
  for (i in seq_along(LIN)) {
    L <- LIN[i]; s <- REF[sub_lin[REF] == L]; y <- i + 0.6
    text(0.01, y - 0.28, LIN_EN[[L]], adj = 0, col = "grey90", cex = 0.9)
    for (k in seq_along(s)) {
      xa <- (k - 1) * cw
      rect(xa + 0.004, y - 0.24, xa + cw - 0.004, y + 0.10,
           col = sub_color[s[k]], border = "grey25", lwd = 0.4)
      text(xa + cw / 2, y + 0.14, s[k], srt = 90, adj = 0, col = "grey80",
           cex = 0.55, xpd = NA)
    }
  }
  invisible(dev.off())
}

## ——— 7. 出图 ———
geo_rows <- list(); sheet <- list()
for (S in SLIDES) {
  cat(sprintf("画 %s …\n", S))
  st <- sub("-\\d+$", "", sub("^GSM[0-9]+_P[0-9]+_", "", S))
  sl <- load_slide(S)
  geo_rows[[S]] <- data.frame(slide = S, n_spot = sl$n, yasp = sl$yasp,
                              flip_col = sl$flip_c, flip_row = sl$flip_r, row.names = NULL)

  sheet[[S]] <- sl            # 全队列一览要用（两种模式都要留）
  if (SHEET_ONLY) next

  ## (a) 合并图 —— **纯图，零文字**
  png(file.path(OUT, sprintf("%s_subtype_composite.png", S)),
      width = 1700, height = 1900, res = 150, bg = "black")
  draw_spots(sl$x, sl$y, cols_of(sl, "composite"), sl$yasp,
             pad = 0.6, mar = c(0, 0, 0, 0))
  invisible(dev.off())

  ## (b) 逐谱系通道（6 格，保留每格名称）
  png(file.path(OUT, sprintf("%s_lineage_channels.png", S)),
      width = 2500, height = 2000, res = 150, bg = "black")
  par(mfrow = c(2, 3), bg = "black", oma = c(0, 0, 3.0, 0))
  for (L in LIN) {
    sub <- REF[sub_lin[REF] == L]
    draw_spots(sl$x, sl$y, cols_of(sl, "lineage", L), sl$yasp,
               sprintf("%s (%d)\n%s", LIN_EN[[L]], length(sub),
                       wrap(paste(sub, collapse = ", "), 78)),
               cex.main = 0.60, mar = c(0.5, 0.5, 5.0, 0.5))
  }
  mtext(sprintf("RCTD subtype fraction by lineage channel  |  %s  (%s)", S, st),
        side = 3, line = 1.4, outer = TRUE, col = "grey85", cex = 0.95)
  invisible(dev.off())

  ## (c) 逐亚型 39 格（保留每格名称）
  png(file.path(OUT, sprintf("%s_subtype_panels.png", S)),
      width = 2800, height = 3000, res = 150, bg = "black")
  par(mfrow = c(6, 7), bg = "black", oma = c(0, 0, 2.8, 0))
  for (k in REF) {
    draw_spots(sl$x, sl$y, cols_of(sl, "subtype", k), sl$yasp,
               k, cex.main = 0.55, col.main = sub_color[[k]],
               mar = c(0.3, 0.3, 2.2, 0.3))
  }
  mtext(sprintf("RCTD per-subtype fraction (39 subtypes, one panel each)  |  %s  (%s)", S, st),
        side = 3, line = 1.0, outer = TRUE, col = "grey85", cex = 0.95)
  invisible(dev.off())
}

## ——— 8. 标注与对照表落盘 ———
key_bars(); key_grid()
write.table(data.frame(subtype = REF, lineage = unname(sub_lin[REF]),
                       abundance_cells = as.integer(abund[REF]),
                       color = unname(sub_color[REF]),
                       rel_luminance = round(ycol(sub_color[REF]), 4)),
            file.path(OUT, "subtype_lineage_color_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (length(geo_rows))
  write.table(do.call(rbind, geo_rows), file.path(OUT, "slide_geometry_audit.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

## ——— 9. 全队列一览：56 张合并图按期别拼成一张 —— **纯图，零文字** ———
if (SHEET_ONLY && !length(sheet)) {
  for (S in SLIDES) { cat(sprintf("读 %s …\n", S)); sheet[[S]] <- load_slide(S) }
}
meta <- do.call(rbind, lapply(names(sheet), function(S) {
  m <- regmatches(S, regexec("^(GSM[0-9]+)_(P[0-9]+)_(.+)$", S))[[1]]
  if (length(m) < 4) return(NULL)
  data.frame(slide = S, patient = m[3], tag = m[4], stage = sub("-\\d+$", "", m[4]),
             row.names = NULL)
}))
meta$ord  <- match(meta$stage, c("Normal", "AAH", "AIS", "MIA", "LUAD"))
meta$pnum <- as.integer(sub("^P", "", meta$patient))
meta <- meta[order(meta$ord, meta$pnum, meta$tag), ]
## 只在**队列齐全**时才写全队列一览——否则单张切片冒烟会把整张图覆盖成 1 格
ALL <- length(SLIDES) == length(list.files(WDIR, pattern = "\\.weights\\.tsv\\.gz$"))
if (SHEET_ONLY || ALL) {
  NCOL <- 8; NROW <- ceiling(nrow(meta) / NCOL)
  png(file.path(OUT, "cohort_overview_composite.png"),
      width = 4200, height = 4200, res = 150, bg = "black")
  par(mfrow = c(NROW, NCOL), bg = "black", oma = c(0, 0, 0, 0))
  for (i in seq_len(nrow(meta))) {
    sl <- sheet[[meta$slide[i]]]
    draw_spots(sl$x, sl$y, cols_of(sl, "composite"), sl$yasp,
               pad = 0.6, mar = c(0.15, 0.15, 0.15, 0.15))
  }
  invisible(dev.off())
  ## 纯图没有标题 ⇒ 位置对照单独落盘，供查格
  meta$panel <- seq_len(nrow(meta))
  write.table(meta[, c("panel", "slide", "patient", "tag", "stage")],
              file.path(OUT, "cohort_overview_order.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n全队列一览：%d 格（%d 列，按期别 Normal→AAH→AIS→MIA→LUAD 排）\n",
              nrow(meta), NCOL))
} else {
  cat("\n（本次只画了部分切片，跳过全队列一览，避免覆盖成残缺版）\n")
}
cat(sprintf("产物目录：%s\n", OUT))
