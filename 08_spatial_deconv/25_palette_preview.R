#!/usr/bin/env Rscript
# 25_palette_preview.R —— **只**做配色候选的预览，供拍板；不改任何数据、不跑全量
#
# 两套候选（39 亚型的"梯度标注"）
#   A 谱系色弧：每个谱系仍占色轮上一段弧（锚色相不变），弧内**同时**走色相与明度
#               ⇒ 同一谱系内 39 个亚型彼此可分辨，谱系身份仍一眼可辨
#   B 单一色相：色相严格等于谱系色，只用**饱和度+明度**做深浅阶梯
#               ⇒ 每个亚型都是它所属谱系色的一个浓淡，与 6 谱系通道图自洽
#
# 两套共用同一条**梯度轴**：谱系内**细胞数降序**（来自 05_annotation 的 cluster 表，
# 表里 cluster 本身按细胞数排，argmax 即亚型名）。排序是数据里的，不是我编的生物学次序。
#
# 🔴 梯度只表示"在这条色阶上的第几位"，**不表示任何每 spot 的丰度/恶性**。
#    RCTD 权重绝不当恶性判读。

RES <- "/home/eto/luad_v2/results/08_spatial_deconv"
DAT <- "/home/eto/luad_v2/data/visium_spatial"
ANN <- "/home/eto/luad_v2/results/05_annotation"
OUT <- file.path(RES, "spatial_cnv/cohort_borrow/subtype_fluorescence")
PRE <- file.path(OUT, "palette_preview")
dir.create(PRE, showWarnings = FALSE)
WDIR <- file.path(RES, "rctd_d/per_slide")

LIN <- c("上皮", "成纤维", "髓系", "内皮", "T/NK", "B/浆")
LIN_EN <- c(上皮 = "Epithelial", 成纤维 = "Fibroblast", 髓系 = "Myeloid",
            内皮 = "Endothelial", "T/NK" = "T/NK", "B/浆" = "B/Plasma")
HUE <- c(上皮 = 12, 成纤维 = 52, 髓系 = 125, 内皮 = 182, "T/NK" = 235, "B/浆" = 302)

## ——— 1. 梯度轴：谱系内细胞数降序 ———
CM <- read.delim(file.path(OUT, "subtype_lineage_color_map.tsv"), stringsAsFactors = FALSE)
REF <- CM$subtype
sub_lin <- setNames(CM$lineage, CM$subtype)
stopifnot(length(REF) == 39, !anyDuplicated(REF), all(sub_lin %in% LIN))

AFILES <- c(上皮 = "epiA_cluster_annotation.csv", 成纤维 = "fibroA_cluster_annotation.csv",
            髓系 = "myeloidA_s4_cluster_annotation.csv", 内皮 = "endoA_s4_cluster_annotation.csv",
            "T/NK" = "tnkA_cluster_annotation.csv", "B/浆" = "bplasmaA_cluster_annotation.csv")
abund <- setNames(rep(0, length(REF)), REF)
for (L in LIN) {
  a <- read.csv(file.path(ANN, AFILES[[L]]), stringsAsFactors = FALSE)
  a <- a[!(toupper(as.character(a$excluded)) %in% c("TRUE", "T")), ]   # 只数进了注释的簇
  nm <- gsub("/", "_", a$argmax)
  for (i in seq_len(nrow(a))) if (nm[i] %in% REF) abund[nm[i]] <- abund[nm[i]] + a$n_cells[i]
}
stopifnot(all(abund > 0))

ORD <- unlist(lapply(LIN, function(L) { s <- REF[sub_lin == L]; s[order(abund[s], decreasing = TRUE)] }))
stopifnot(length(ORD) == 39, setequal(ORD, REF))
cat("梯度轴（谱系内细胞数降序）：\n")
for (L in LIN) {
  s <- ORD[sub_lin[ORD] == L]
  cat(sprintf("  %-6s %s\n", L, paste(sprintf("%s(%d)", s, abund[s]), collapse = " > ")))
}

## ——— 2. 两套候选配色 ———
## A 的弧段端点（度）：相邻谱系的锚色相取中点、再各让开 8°，保证弧与弧之间有空隙；
## B/浆 只有 2 个亚型，额外收窄以保住"洋红系"的身份。
ARC <- list(上皮 = c(341, 388), 成纤维 = c(36, 78), 髓系 = c(86, 149.5),
            内皮 = c(162, 204.5), "T/NK" = c(212.5, 264.5), "B/浆" = c(285, 320))
## 相对亮度（Rec.709）。C 方案就是拿它当配平目标。
RELY <- c(0.2126, 0.7152, 0.0722)
relY  <- function(v) sum(RELY * v)
hsv_rgb <- function(h, s, v) grDevices::col2rgb(hsv(h / 360, s, v))[, 1] / 255
## 在给定色相上找"亮度 = Y0"的那个颜色：先看纯色（s=1,v=1）够不够亮
##   够亮 ⇒ 压暗（v = Y0/Y_pure）；不够亮 ⇒ 加白（v=1，解 s）
## 注意：hsv() 返回的是 8 位十六进制，颜色本身就带取整误差，所以最后要按
## **取整后**的亮度再微调一步——否则 39 个色的亮度各差 ±0.5/255，合成图的亮度
## 又会随组成轻微起伏（幅度约 0.35%，肉眼看不出来，但相关系数会把它放大）。
eqY <- function(hue, Y0 = 0.42) {
  hue <- hue %% 360                       # 上皮的弧跨 360°，色相会越界
  y1 <- relY(grDevices::col2rgb(hsv(hue / 360, 1, 1))[, 1] / 255)
  ## 两个自由度：主刻度（压暗 v / 加白 s）+ 色相 ±4° 微调；
  ## 后者是为了绕开 8 位量化的格点，把 39 个色的亮度对齐到同一个值。
  dh <- seq(-4, 4, by = 0.5)
  sc <- if (y1 >= Y0) {
    v <- Y0 / y1; seq(max(1e-6, v * 0.97), min(1, v * 1.03), length.out = 121)
  } else {
    s <- uniroot(function(s) relY(hsv_rgb(hue, s, 1)) - Y0, c(0, 1), tol = 1e-7)$root
    seq(max(1e-6, s - 0.02), min(1, s + 0.02), length.out = 121)
  }
  H <- rep((hue + dh) %% 360, each = length(sc))
  cand <- if (y1 >= Y0) hsv(H / 360, 1, rep(sc, times = length(dh)))
          else          hsv(H / 360, rep(sc, times = length(dh)), 1)
  Y <- as.vector(RELY %*% (grDevices::col2rgb(cand) / 255))
  cand[[which.min(abs(Y - Y0))]]
}

mk <- function(mode = c("A", "B", "C")) {
  mode <- match.arg(mode)
  out <- character(0)
  for (L in LIN) {
    s <- ORD[sub_lin[ORD] == L]; n <- length(s)
    t <- if (n == 1) 0 else seq(0, 1, length.out = n)     # 0 = 该谱系最常见
    a <- ARC[[L]]
    hue <- a[1] + (a[2] - a[1]) * t
    ## t = 0 是**该谱系最常见**的亚型 ⇒ 放在最浓/最饱和一端，稀有的一端才淡下去；
    ## 反过来（常见=最浅）会把主导亚型冲淡，整张合并图变灰。
    out[s] <- if (mode == "A") {
      hsv((hue %% 360) / 360, s = 0.95, v = 1.00 - 0.22 * t)
    } else if (mode == "B") {
      hsv(HUE[[L]] / 360, s = 1.00 - 0.25 * t, v = 1.00 - 0.45 * t)
    } else {
      ## C：39 色**相对亮度全部 = 0.42** ⇒ 图上不存在明度轴
      vapply(hue, eqY, "", Y0 = 0.42)
    }
  }
  out[REF]
}

## ——— 3. 画法（与 24_ 同：逐 spot 正六边形，无插值、不填洞）———
HEX_R  <- 1 / sqrt(3)
HEX_TH <- (30 + 60 * (0:5)) * pi / 180
HEX_VX <- HEX_R * cos(HEX_TH); HEX_VY <- HEX_R * sin(HEX_TH)
draw_spots <- function(x, y, cols, yasp, pad = 0.9) {
  plot.new()
  plot.window(xlim = range(x) + c(-pad, pad), ylim = range(y) + c(-pad, pad) * yasp, asp = 1)
  n <- length(x); X <- matrix(NA_real_, n, 7L); Y <- matrix(NA_real_, n, 7L)
  X[, 1:6] <- outer(x, HEX_VX, "+"); Y[, 1:6] <- outer(y, HEX_VY, "+")
  polygon(as.vector(t(X)), as.vector(t(Y)), col = cols, border = cols, lwd = 0.25)
}
load_slide <- function(S) {
  w <- read.delim(gzfile(file.path(WDIR, paste0(S, ".weights.tsv.gz"))), row.names = 1, check.names = FALSE)
  W <- as.matrix(w[, REF]); rownames(W) <- NULL
  p <- read.csv(file.path(DAT, S, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
  m <- match(rownames(w), p$barcode); if (anyNA(m)) stop("坐标对不上：", S)
  p1 <- p[p$in_tissue == 1, ]
  yasp <- (diff(range(p1$pxl_row_in_fullres)) / diff(range(p1$array_row))) /
          (2 * diff(range(p1$pxl_col_in_fullres)) / diff(range(p1$array_col)))
  x <- p$array_col[m] / 2; y <- p$array_row[m] * yasp
  if (cor(p1$array_col, p1$pxl_col_in_fullres) < 0) x <- (max(p$array_col) / 2) - x
  if (cor(p1$array_row, p1$pxl_row_in_fullres) < 0) y <- (max(p$array_row) * yasp) - y
  list(W = W, x = x, y = y, yasp = yasp, n = nrow(w))
}
cols_rgb <- function(W, pal, norm = FALSE, gamma = 0.75) {
  if (norm) W <- W / rowSums(W)          # 只作显示：把每个 spot 的构成变成占比
  RGB <- pmin(pmax(W %*% t(grDevices::col2rgb(pal[REF]) / 255), 0), 1)
  if (gamma != 1) RGB <- RGB ^ gamma
  rgb(RGB)
}

## ——— 4. 梯度标注图：每个谱系一条色阶条，条上按顺序标亚型名 ———
key_bars <- function(pal, tag) {
  png(file.path(PRE, sprintf("key_gradient_bars_%s.png", tag)),
      width = 2600, height = 3000, res = 150, bg = "black")
  par(mar = c(0.4, 0.4, 0.4, 0.4), bg = "black")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(10.4, 0.1))   # 行距 1.5，给竖排名字留够空间
  text(0.02, 0.30, sprintf("39-subtype colour ramp  —  scheme %s", tag),
       adj = 0, col = "grey90", cex = 1.6)
  text(0.02, 0.62, "order inside each lineage = descending cell number (05_annotation cluster tables)  ·  position/shade = rank only, NOT a per-spot quantity",
       adj = 0, col = "grey55", cex = 0.62)
  x0 <- 0.20; x1 <- 0.90
  for (i in seq_along(LIN)) {
    L <- LIN[i]; s <- ORD[sub_lin[ORD] == L]; n <- length(s)
    y <- 1.1 + 1.5 * (i - 1)
    rect(x0, y - 0.115, x1, y + 0.115, col = "grey15", border = "grey30")
    xs <- seq(x0, x1, length.out = n + 1)
    for (k in seq_len(n)) rect(xs[k], y - 0.105, xs[k + 1], y + 0.105,
                               col = pal[s[k]], border = "grey20", lwd = 0.4)
    text(x0 - 0.012, y - 0.05, sprintf("%s  (%d)", LIN_EN[[L]], n), adj = 1,
         col = "grey90", cex = 0.95)
    text(x0 - 0.012, y + 0.10, "high → low abundance", adj = 1,
         col = "grey45", cex = 0.5)
    ## 名字从条下方**向下**竖排（srt=90 + adj=1），不再压到上一行
    for (k in seq_len(n))
      text((xs[k] + xs[k + 1]) / 2, y + 0.16, s[k], srt = 90, adj = 1,
           col = pal[s[k]], cex = 0.60, xpd = NA)
    text(x1 + 0.010, y, formatC(abund[s], format = "d", big.mark = ","),
         adj = 0, col = "grey55", cex = 0.52, xpd = NA)
  }
  dev.off()
}

## ——— 5. 色卡：6 行（谱系）× 列（该谱系亚型），一行内按丰度降序 ———
key_grid <- function(pal, tag) {
  NCOL <- 12
  png(file.path(PRE, sprintf("key_swatch_grid_%s.png", tag)),
      width = 2600, height = 2000, res = 150, bg = "black")
  par(mar = c(0.5, 0.5, 3.0, 0.5), bg = "black")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(6.6, 0.4))
  text(0.02, 0.55, sprintf("39-subtype swatch key  —  scheme %s", tag),
       adj = 0, col = "grey90", cex = 1.5)
  cw <- 1 / NCOL
  for (i in seq_along(LIN)) {
    L <- LIN[i]; s <- ORD[sub_lin[ORD] == L]
    y <- i + 0.6
    text(0.01, y - 0.28, LIN_EN[[L]], adj = 0, col = "grey90", cex = 0.9)
    for (k in seq_along(s)) {
      xa <- (k - 1) * cw
      rect(xa + 0.004, y - 0.24, xa + cw - 0.004, y + 0.10, col = pal[s[k]],
           border = "grey25", lwd = 0.4)
      text(xa + cw / 2, y + 0.14, s[k], srt = 90, adj = 0, col = "grey80",
           cex = 0.55, xpd = NA)
    }
  }
  dev.off()
}

## ——— 6. 跑：三套各出标注图 + 一张**纯图** demo（P10_LUAD，图内零文字）———
##  A/B 沿用 23_ 的加法混色 + γ0.75（与旧图可比）；C 是等亮度，必须**关掉 γ**
## 并**把每 spot 权重归一到占比**，否则亮度轴会从"权重总量"那条路绕回来。
CFG <- list(A = list(norm = FALSE, gamma = 0.75),
            B = list(norm = FALSE, gamma = 0.75),
            C = list(norm = TRUE,  gamma = 1))
S <- "GSM9226190_P10_LUAD"
sl <- load_slide(S)
pals <- list()
for (tag in names(CFG)) {
  pal <- mk(tag); pals[[tag]] <- pal
  key_bars(pal, tag)
  key_grid(pal, tag)
  png(file.path(PRE, sprintf("scheme%s_demo_%s.png", tag, S)),
      width = 1500, height = 1700, res = 150, bg = "black")
  par(mar = c(0, 0, 0, 0), bg = "black")
  draw_spots(sl$x, sl$y, cols_rgb(sl$W, pal, CFG[[tag]]$norm, CFG[[tag]]$gamma),
             sl$yasp, pad = 0.35)
  dev.off()
  write.table(data.frame(subtype = ORD, lineage = unname(sub_lin[ORD]),
                         abundance_cells = abund[ORD], color = unname(pal[ORD]),
                         rel_luminance = round(apply(grDevices::col2rgb(pal[ORD]), 2, relY), 4)),
              file.path(PRE, sprintf("palette_%s.tsv", tag)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("scheme %s 完成（配色相对亮度 中位 %.3f，极差 %.3f）\n", tag,
              median(apply(grDevices::col2rgb(pal[ORD]), 2, relY)),
              diff(range(apply(grDevices::col2rgb(pal[ORD]), 2, relY)))))
}

## ——— 7. 体检：spot 亮度到底还跟不跟深度走 ———
##   用户提的问题——"深浅会不会和细胞密度有关"。图上说法不算数，这里直接量。
QD <- read.csv(file.path(RES, "spatial_qc_per_spot.csv.gz"), stringsAsFactors = FALSE)
QD$key <- paste(QD$slide, sub("-1$", "", QD$barcode), sep = "|")
key <- paste(S, sub("-1$", "", rownames(read.delim(
  gzfile(file.path(WDIR, paste0(S, ".weights.tsv.gz"))), row.names = 1))
), sep = "|")
numi <- QD$nUMI[match(key, QD$key)]
num_rgb <- function(W, pal, norm) {
  if (norm) W <- W / rowSums(W)
  pmin(pmax(W %*% t(grDevices::col2rgb(pal[REF]) / 255), 0), 1)
}
chk <- data.frame(scheme = character(), corr_light_vs_depth = numeric(),
                  corr_light_vs_depth_8bit = numeric(),
                  light_range_255 = numeric(), stringsAsFactors = FALSE)
rs <- rowSums(sl$W)
for (tag in names(CFG)) {
  M <- num_rgb(sl$W, pals[[tag]], CFG[[tag]]$norm)
  if (CFG[[tag]]$gamma != 1) M <- M ^ CFG[[tag]]$gamma
  lum <- as.vector(M %*% RELY)                                   # 未经 8 位量化
  lq  <- as.vector(t(grDevices::col2rgb(rgb(M)) / 255) %*% RELY) # 落盘时的 8 位量化值
  chk <- rbind(chk, data.frame(
    scheme = tag,
    corr_light_vs_depth      = round(cor(lum, numi, method = "spearman"), 3),
    corr_light_vs_depth_8bit = round(cor(lq,  numi, method = "spearman"), 3),
    light_range_255          = round(255 * diff(range(lq)), 2)))
}
write.table(chk, file.path(PRE, "brightness_vs_depth_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nspot 亮度 vs 测序深度（%s）· 权重行和 vs 深度 Spearman = %.3f\n", S,
            cor(rs, numi, method = "spearman")))
print(chk)
cat(sprintf("\n产物：%s\n", PRE))
