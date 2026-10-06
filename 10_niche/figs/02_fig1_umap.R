#!/usr/bin/env Rscript
# 02_fig1_umap.R —— Fig 1：单细胞图谱（谱系 / 分期 / 分期拆分）
source("00_palette_theme.R")
D <- fread("/tmp/_cache_cellmap.csv")
D[, stage := factor(stage, levels = STAGE_LEVELS)]
setnames(D, "A_frozen", "L1")
D[, L1 := factor(L1, levels = LIN_LEVELS)]
cat("细胞", nrow(D), "| 谱系", paste(levels(D$L1), collapse="/"), "\n")

pts <- geom_point(size = 0.06, stroke = 0, shape = 16)
base <- ggplot(D, aes(umap_1, umap_2)) + theme_umap() + labs(x = "UMAP-1", y = "UMAP-2")

# A 谱系
pA <- base + pts + aes(colour = L1) + scale_colour_manual(values = LIN_COL, name = NULL) +
  labs(title = "Lineage") + guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
# B 分期
pB <- base + pts + aes(colour = stage) + scale_colour_manual(values = STAGE_COL, name = NULL) +
  labs(title = "Stage") + guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
# C 分期拆分
pC <- base + pts + aes(colour = stage) + scale_colour_manual(values = STAGE_COL, guide = "none") +
  facet_wrap(~stage, nrow = 1) + labs(title = "Stage (split)") +
  theme(strip.text = element_text(colour = "white", face = "bold"),
        strip.background = element_rect(fill = "grey35", linewidth = 0))
pC <- pC + theme(panel.spacing = unit(1, "pt"))

fig <- (pA | pB) / pC + plot_layout(heights = c(1, 0.95)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_fig(fig, "Fig1_singlecell_UMAP", 7.2, 5.4)
