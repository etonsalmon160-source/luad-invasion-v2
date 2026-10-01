#!/usr/bin/env Rscript
# 22_spatial_cnv_map.R —— 逐点空间图：CNV 负荷 + RCTD 上皮权重（只读）
#
# 🔴 命名禁令（项目既定，不得违反）：
#   - `cnv_fraction` 是**连续 CNV 负荷**，**不是**恶性判定，图注/色标**不得**写 malignancy / tumor
#   - RCTD 权重**绝不当恶性判读**；本图只用 `上皮` 这一列的**占比**
#   - 两幅面板**并列**呈现、**不合成**任何单一「恶性分」
#
# 口径来源：main_borrow_spot_scores.rds（cf，只含 observation spot）
#           rctd_a/per_slide/<slide>.weights.tsv.gz（六谱系）
#           data/visium_spatial/<slide>/spatial/tissue_positions.csv（坐标）
# 坐标用 array_row / array_col：不用 pxl（边缘可为负，见 reference_visium_coordinate_conventions）

RES <- "/home/eto/luad_v2/results/08_spatial_deconv"
DAT <- "/home/eto/luad_v2/data/visium_spatial"
OUT <- file.path(RES, "spatial_cnv/cohort_borrow/spatial_map")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PATIENT <- "P21"
SLIDES  <- c("GSM9226211_P21_AIS", "GSM9226212_P21_AIS-1", "GSM9226213_P21_LUAD")

## ——— 读 cf（一次，P21 的 main 臂）———
o  <- readRDS(file.path(RES, "spatial_cnv/cohort_borrow", PATIENT, "main_borrow_spot_scores.rds"))
io <- which(o$grp == "observation")

strip_slide <- function(s, b) sub(paste0("^", s, "_"), "", b)

## ——— 逐切片组装 ———
d <- do.call(rbind, lapply(SLIDES, function(S) {
  i  <- io[o$slide[io] == S]
  bc <- strip_slide(S, o$bc[i])
  w  <- read.delim(gzfile(file.path(RES, "rctd_a/per_slide", paste0(S, ".weights.tsv.gz"))),
                   row.names = 1, check.names = FALSE)
  pos <- read.csv(file.path(DAT, S, "spatial/tissue_positions.csv"), stringsAsFactors = FALSE)
  mw <- match(bc, rownames(w)); mp <- match(bc, pos$barcode)
  if (anyNA(mw) || anyNA(mp)) stop("barcode 对不上：", S)
  epi <- w[mw, "上皮"]
  ord <- sort(epi, decreasing = TRUE)
  data.frame(slide = S, cf = o$cf[i], epi = epi, margin = ord[1] - ord[2],
             x = pos$array_col[mp] / 2, y = -pos$array_row[mp], row.names = NULL)
}))
cat(sprintf("%s：%d 张切片，%d 个观测 spot，barcode 全命中\n",
            PATIENT, length(SLIDES), nrow(d)))

## ——— 画 ———
## 🔴 `hcl.colors(..., "YlOrRd")` 在 R 4.2 里 **索引 1 = 深红 #7D0025、索引 64 = 浅黄 #FFFFC8**，
##    与 RColorBrewer 的 YlOrRd 方向**相反**（静默画反过一版，登记）。故必须 `rev = TRUE`。
pal_cf  <- hcl.colors(64, "YlOrRd", rev = TRUE)   # 1=浅黄(低) → 64=深红(高)
pal_epi <- hcl.colors(64, "Blues", rev = TRUE)    # 1=近白(低) → 64=深蓝(高)
stopifnot(pal_cf[1] == "#FFFFC8", pal_epi[1] == "#F4FAFE")

XR <- range(d$x); YR <- range(d$y)
XL <- c(XR[1] - 3, XR[2] + 26)                    # 右侧留出图例位
colbar <- function(pal, lo, hi, lab) {
  x0 <- XR[2] + 6; x1 <- x0 + 9
  y0 <- YR[1] + 12; y1 <- y0 + 74
  yy <- seq(y0, y1, length.out = length(pal) + 1)
  for (k in seq_along(pal))
    rect(x0, yy[k], x1, yy[k + 1], col = pal[k], border = NA)
  rect(x0, y0, x1, y1, border = "grey30")
  text(x0 - 2, y0, sprintf("%.2f", lo), cex = 0.62, adj = c(1, 0.5))
  text(x0 - 2, y1, sprintf("%.2f", hi), cex = 0.62, adj = c(1, 0.5))
  text(x1 + 5, (y0 + y1) / 2, lab, srt = 90, cex = 0.68)
}

png(file.path(OUT, sprintf("%s_spatial_cnv_epi.png", PATIENT)),
    width = 1560, height = 2100, res = 150)
op <- par(mfrow = c(3, 2), mar = c(0.6, 0.6, 2.6, 0.6), oma = c(1.5, 1.5, 3.0, 1.0))
cf_max <- as.numeric(quantile(d$cf, .995))

for (S in SLIDES) {
  z  <- d[d$slide == S, ]
  st <- sub("-\\d+$", "", sub("^GSM[0-9]+_P[0-9]+_", "", S))
  ## 左：CNV 负荷
  plot(z$x, z$y, pch = 15, cex = 0.42, asp = 1, axes = FALSE, xlab = "", ylab = "",
       xlim = XL, ylim = YR,
       col = pal_cf[pmax(1L, pmin(64L, as.integer(round(z$cf / cf_max * 63)) + 1L))],
       main = sprintf("%s  (%s)  |  CNV burden: cnv_fraction", S, st), cex.main = 0.8)
  colbar(pal_cf, 0, cf_max, "cnv_fraction")
  ## 右：RCTD 上皮权重
  plot(z$x, z$y, pch = 15, cex = 0.42, asp = 1, axes = FALSE, xlab = "", ylab = "",
       xlim = XL, ylim = YR,
       col = pal_epi[pmax(1L, pmin(64L, as.integer(round(z$epi * 63)) + 1L))],
       main = sprintf("%s  (%s)  |  RCTD epithelial fraction", S, st), cex.main = 0.8)
  colbar(pal_epi, 0, 1, "epi fraction")
}
mtext("P21 — per-spot spatial map (observation spots only). Left: continuous CNV burden score. Right: RCTD epithelial weight. Neither panel is a malignancy call.",
      side = 3, line = 0.8, outer = TRUE, cex = 0.72)
par(op); dev.off()
cat(sprintf("产物：%s\n", file.path(OUT, sprintf("%s_spatial_cnv_epi.png", PATIENT))))

## ——— 汇总数字（供图注引用）———
sm <- do.call(rbind, lapply(SLIDES, function(S) {
  z <- d[d$slide == S, ]
  data.frame(slide = S, stage = sub("-\\d+$", "", sub("^GSM[0-9]+_P[0-9]+_", "", S)),
             n_spot = nrow(z), cf_median = median(z$cf), cf_p90 = quantile(z$cf, .9),
             epi_median = median(z$epi), epi_ge_0.5 = mean(z$epi >= 0.5),
             cor_cf_epi = cor(z$cf, z$epi), row.names = NULL)
}))
write.table(sm, file.path(OUT, sprintf("%s_spatial_cnv_epi_summary.tsv", PATIENT)),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n==== 逐切片汇总 ====\n"); print(sm)
