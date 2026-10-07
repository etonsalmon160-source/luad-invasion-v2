# 00_palette_theme.R —— 全项目统一配色与主题（顶刊风格）
# 用法：source("00_palette_theme.R")
suppressMessages({library(ggplot2); library(data.table); library(ggsci); library(patchwork)})

OUT <- "/home/eto/luad_v2/results/paper_figures"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ── 分期：有序渐变（Normal→IAC），蓝→黄→红，色盲可辨 ──
STAGE_LEVELS <- c("Normal", "AAH", "AIS", "MIA", "IAC")
STAGE_ALIAS <- c("LUAD" = "IAC")   # 🔴 d7_domain_assign.tsv 用的是 GEO 原始 token
STAGE_COL <- c(Normal = "#2C6FA6", AAH = "#74AFD4", AIS = "#F0C24B", MIA = "#E07B39", IAC = "#B02A30")

## ── 六大谱系：NPG 配色（Nature 系列）──
LIN_LEVELS <- c("上皮", "成纤维", "内皮", "髓系", "T/NK", "B/浆")
LIN_EN <- c("上皮"="Epithelial", "成纤维"="Fibroblast", "内皮"="Endothelial",
            "髓系"="Myeloid", "T/NK"="T/NK", "B/浆"="B/Plasma")
LIN_COL <- setNames(pal_npg("nrc")(6), LIN_EN[LIN_LEVELS])   # 图例用英文名

## ── 七个生态位域类型（D1–D7）──
ARCH_LEVELS <- paste0("D", 1:7)
## 2026-10-06 定：蓝紫粉家族；D4 是占 45% spot 的肺泡毛细血管背景域 ⇒ 用最浅粉后退
ARCH_COL <- setNames(c("#2F5597", "#5B4E9E", "#9C4E93", "#F7DCEA",
                       "#D95F8C", "#2E8B8B", "#6A3FA0"), ARCH_LEVELS)

## ── 强调色（单值高亮用）──
ACCENT <- "#B02A30"; NEUTRAL <- "#BDBDBD"; HILITE <- "#E07B39"

## ── 主题：干净、无多余网格、字体统一 ──
theme_paper <- function(base_size = 8, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
      panel.border = element_rect(linewidth = 0.4, colour = "grey30"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
      axis.text = element_text(colour = "grey15"),
      axis.title = element_text(colour = "grey5"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "grey35", hjust = 0),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      legend.key.size = unit(0.32, "cm"),
      legend.margin = margin(1, 1, 1, 1),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "grey95", linewidth = 0.3, colour = "grey30"),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      plot.margin = margin(3, 3, 3, 3)
    )
}
norm_stage <- function(x) { x <- as.character(x); ifelse(x %in% names(STAGE_ALIAS), STAGE_ALIAS[x], x) }
theme_set(theme_paper())

## ── 统一的 UMAP 底图（无坐标轴数字，只留框）──
theme_umap <- function(base_size = 8) {
  theme_paper(base_size) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_text(size = base_size - 1, colour = "grey30"),
          panel.grid = element_blank())
}

save_fig <- function(p, name, w, h, dpi = 400) {
  f <- file.path(OUT, paste0(name, ".pdf"))
  ggsave(f, p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = dpi)
  cat(sprintf("  ✓ %-34s %.1f×%.1f in\n", paste0(name, ".pdf/.png"), w, h))
}

## ── 数据加载器（缓存）──
load_cellmap <- function() {
  D <- as.data.table(readRDS("/tmp/_cache_cellmap.rds"))  # 由 01_ 预转
  D[, stage := factor(stage, levels = STAGE_LEVELS)]
  D[, L1 := factor(LIN_EN[L1], levels = LIN_EN[LIN_LEVELS])]
  D
}

## ── 荧光配色（移植自 08_spatial_deconv/23_spatial_subtype_fluorescence.R，RCTD 39 亚型）──
FLUO_HUE  <- c("上皮"=12, "成纤维"=52, "髓系"=125, "内皮"=182, "T/NK"=235, "B/浆"=302)
FLUO_SPAN <- c("上皮"=26, "成纤维"=16, "髓系"=62, "内皮"=26, "T/NK"=30, "B/浆"=26)
make_fluo_palette <- function(ref, sub_lin) {          # ref=39亚型名向量; sub_lin=每条对应谱系
  col <- character(length(ref)); names(col) <- ref
  for (L in names(FLUO_HUE)) {
    idx <- which(sub_lin == L); n <- length(idx); if (!n) next
    off <- if (n == 1) 0 else seq(-FLUO_SPAN[[L]]/2, FLUO_SPAN[[L]]/2, length.out = n)
    v   <- if (n == 1) 1 else seq(0.85, 1.0, length.out = n)[order(off)]
    col[idx] <- hsv((FLUO_HUE[[L]] + off) %% 360 / 360, s = 0.95, v = v)
  }
  col
}
