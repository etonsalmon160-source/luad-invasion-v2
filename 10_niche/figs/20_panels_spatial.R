#!/usr/bin/env Rscript
# 20_panels_spatial.R —— 空间基础图：**每个签名一张，全英文**
#   P4a 生态位域类型 / P4b 耗竭 T / P4c TLS / P4d ECM 成纤维 / P4e iCAF
source("00_palette_theme.R")
ROOT <- "/home/eto/luad_v2"
S <- fread(cmd = sprintf("zcat %s/results/10_niche/kstar_diag/d12_spot_signatures.tsv.gz", ROOT))
S <- S[slide != ""]
st <- fread(sprintf("%s/results/10_niche/kstar_diag/d7_domain_assign.tsv", ROOT),
            select = c("slide", "stage"))[!duplicated(slide)]
S <- merge(S, st, by = "slide")
S[, stage := factor(norm_stage(stage), levels = STAGE_LEVELS)]

# 坐标（逐切片读 tissue_positions）
pos <- rbindlist(lapply(unique(S$slide), function(sl) {
  f <- sprintf("%s/data/visium_spatial/%s/spatial/tissue_positions.csv", ROOT, sl)
  if (!file.exists(f)) return(NULL)
  d <- fread(f, header = FALSE, skip = 1,
             col.names = c("barcode", "in_tissue", "array_row", "array_col", "px", "py"))
  d[, slide := sl][, .(slide, barcode, px, py, in_tissue)]
}))
S <- merge(S, pos, by = c("slide", "barcode"))
S <- S[is.na(in_tissue) | in_tissue == 1]
cat("spots with coords:", nrow(S), "\n")
# 坐标翻转：px 是行（向下），py 是列；画图时 y 取负以正确朝向
S[, y := -py][, x := px]

## 每个期别取 spot 数最多的一张作代表
REP <- S[, .N, by = .(stage, slide)][order(-N)][, .SD[1], by = stage]
cat("代表切片:\n"); print(REP)
W <- S[slide %in% REP$slide]
W[, stage := factor(stage, levels = STAGE_LEVELS)]

## ── 比例尺：逐切片自标定（相邻 spot 中心距 = 100 µm）──
## 每张切片的 tissue_hires_scalef 不同（0.035–0.043）⇒ 必须逐张算，不能用统一 µm/px
px100 <- function(sl) {
  f <- sprintf("%s/data/visium_spatial/%s/spatial/tissue_positions.csv", ROOT, sl)
  d <- fread(f, header = FALSE, skip = 1,
             col.names = c("barcode", "in_tissue", "array_row", "array_col", "px", "py"))
  d <- d[in_tissue == 1]
  setorder(d, array_row, array_col)                 # 同一 array_row 内相邻列 = 1 个物理步长
  dd <- d[, .(g = abs(diff(py))), by = array_row][g > 0]   # 切片朝向不同，diff 可为负 ⇒ 取绝对值
  median(dd$g)
}
RG <- W[, .(x0 = min(x), x1 = max(x), y0 = min(y), y1 = max(y)), by = stage]
SB <- merge(data.table(stage = REP$stage, slide = REP$slide), RG, by = "stage")
SB[, p100 := sapply(slide, px100)]
SB[, xa := x0 + 0.055 * (x1 - x0)]
SB[, xb := xa + 10 * p100]   # 1 mm
SB[, yb := y0 + 0.055 * (y1 - y0)]
SB[, xt := xb + 0.022 * (x1 - x0)]   # 文字放条右边
cat("比例尺（1 mm 的像素长度 = 10 × 相邻 spot 间距）:\n"); print(SB[, .(stage, slide, px_1mm = round(10 * p100, 1), panel_w = round(x1 - x0, 0))])
SEG <- list(
  geom_segment(data = SB, aes(x = xa, xend = xb, y = yb, yend = yb),
               colour = "black", linewidth = 1.5, inherit.aes = FALSE),
  geom_segment(data = SB, aes(x = xa, xend = xb, y = yb, yend = yb),
               colour = "white", linewidth = 0.8, inherit.aes = FALSE),
  geom_text(data = SB, aes(x = xt, y = yb, label = "1 mm"),
            size = 2.4, inherit.aes = FALSE, colour = "black", fontface = "bold",
            hjust = 0, vjust = 0.5))


sp_plot <- function(var, title, pal, mid = NULL, fmt = "%.2f") {
  d <- copy(W); d[, v := get(var)]
  ggplot(d, aes(x, y, colour = v)) +
    geom_point(size = 0.32, stroke = 0) + SEG +
    scale_colour_gradientn(colours = pal, name = NULL,
                           guide = guide_colourbar(barwidth = unit(2.4, "cm"),
                                                   barheight = unit(0.25, "cm"))) +
    coord_fixed() + facet_wrap(~stage, nrow = 1) +
    labs(x = NULL, y = NULL, title = title) +
    theme_paper(7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), panel.border = element_rect(linewidth = 0.3, colour = "grey60"),
          panel.spacing = unit(2, "pt"),
          strip.text = element_text(colour = "white", face = "bold", size = 7.5),
          strip.background = element_rect(fill = "grey30", linewidth = 0),
          legend.position = "bottom")
}

## P4a · 生态位域类型
d <- W[!is.na(archetype)]; d[, arch := factor(paste0("D", archetype), levels = ARCH_LEVELS)]
DOMLAB <- c(D1 = "Airway", D2 = "iCAF", D3 = "ECM/interstitial", D4 = "Alveolar-cap.",
             D5 = "AT2", D6 = "Vascular", D7 = "Lymphoid")
p <- ggplot(d, aes(x, y, colour = arch)) + geom_point(size = 0.32, stroke = 0) + SEG +
  scale_colour_manual(values = ARCH_COL, name = NULL,
                      labels = function(x) paste0(x, "  ", DOMLAB[x]),
                      guide = guide_legend(nrow = 2, byrow = TRUE,
                                           override.aes = list(size = 2.6, shape = 15))) +
  coord_fixed() +
  facet_wrap(~stage, nrow = 1) + labs(x = NULL, y = NULL, title = "Niche domain archetype (K*=7)") +
  theme_paper(7) + theme(axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), panel.border = element_rect(linewidth = 0.3, colour = "grey60"),
    panel.spacing = unit(2, "pt"), strip.text = element_text(colour = "white", face = "bold", size = 7.5),
    strip.background = element_rect(fill = "grey30", linewidth = 0),
    legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
    legend.text = element_text(size = 6.4), legend.margin = margin(t = -2))
save_fig(p, "P4a_spatial_domain_archetype", 8.6, 2.5)

## P4b–P4e · 各签名
GREY_RED  <- c("#F0F0F0", "#FDD49E", "#FC8D59", "#B30000", "#67000D")
GREY_TEAL <- c("#F0F0F0", "#C7E9E4", "#5AB4AC", "#01665E", "#003C30")
GREY_PURP <- c("#F0F0F0", "#DADAEB", "#9E9AC8", "#6A51A3", "#3F007D")
SIGS <- list(
  list("T_exhaust",  "Exhausted T-cell signature",        GREY_PURP, "P4b_spatial_T_exhaust"),
  list("TLS",        "Tertiary lymphoid structure (TLS)",  GREY_TEAL, "P4c_spatial_TLS"),
  list("Fib_ECM",    "ECM fibroblast signature",          GREY_RED,  "P4d_spatial_Fib_ECM"),
  list("iCAF",       "iCAF signature",                    GREY_TEAL, "P4e_spatial_iCAF"),
  list("Treg",       "Regulatory T-cell signature",       GREY_PURP, "P4f_spatial_Treg"),
  list("Mac_LAM",    "LAM macrophage signature",          GREY_RED,  "P4g_spatial_Mac_LAM")
)
for (s in SIGS) {
  n <- sum(W[[s[[1]]]] > 0); cat(sprintf("  %s: 非零 spot %d/%d\n", s[[1]], n, nrow(W)))
  save_fig(sp_plot(s[[1]], s[[2]], s[[3]]), s[[4]], 8.6, 2.5)
}
