#!/usr/bin/env Rscript
# 10_panels_singlecell.R —— 单细胞基础图：**每个子图单独出文件，标签全英文**
source("00_palette_theme.R")
D <- fread("/tmp/_cache_cellmap.csv")
setnames(D, "A_frozen", "L1")
D[, stage := factor(stage, levels = STAGE_LEVELS)]
D[, L1 := factor(LIN_EN[L1], levels = LIN_EN[LIN_LEVELS])]
cat("cells", nrow(D), "\n")

pts <- geom_point(size = 0.05, stroke = 0, shape = 16)
base <- ggplot(D, aes(umap_1, umap_2)) + theme_umap(8) + labs(x = "UMAP-1", y = "UMAP-2")
LEG <- guide_legend(override.aes = list(size = 2.4, alpha = 1), keyheight = unit(0.32, "cm"))

## ── P1a · UMAP by lineage ──
p <- base + pts + aes(colour = L1) +
  scale_colour_manual(values = LIN_COL, name = NULL) +
  guides(colour = LEG) + labs(title = "Single-cell atlas: lineage")
save_fig(p, "P1a_UMAP_lineage", 3.6, 3.3)

## ── P1b · UMAP by stage ──
p <- base + pts + aes(colour = stage) +
  scale_colour_manual(values = STAGE_COL, name = NULL) +
  guides(colour = LEG) + labs(title = "Single-cell atlas: stage")
save_fig(p, "P1b_UMAP_stage", 3.6, 3.3)

## ── P1c · UMAP split by stage ──
p <- base + pts + aes(colour = stage) + scale_colour_manual(values = STAGE_COL, guide = "none") +
  facet_wrap(~stage, nrow = 2) + labs(title = "Single-cell atlas split by stage") +
  theme(panel.spacing = unit(1.5, "pt"),
        strip.text = element_text(colour = "white", face = "bold", size = 8),
        strip.background = element_rect(fill = "grey30", linewidth = 0))
save_fig(p, "P1c_UMAP_stage_split", 4.6, 4.4)

## ── P2a · lineage composition by stage ──
s <- D[, .N, by = .(stage, L1)][, frac := N / sum(N), by = stage]
p <- ggplot(s, aes(stage, frac, fill = L1)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = LIN_COL, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0.004, 0), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Fraction of nuclei", title = "Lineage composition by stage") +
  theme_paper(8) + guides(fill = guide_legend(keyheight = unit(0.3, "cm")))
save_fig(p, "P2a_lineage_composition_by_stage", 3.2, 3.0)

## ── P2b · L2 subtype composition by stage ──
L2S <- D[, .N, by = .(L2, L1)][order(-N)][, .SD[1], by = L2][, .(L2, L1)]; setorder(L2S, L1, L2)
shade <- function(hex, f) { r <- col2rgb(hex)/255; rgb(r[1]*f+(1-f)*1, r[2]*f+(1-f)*1, r[3]*f+(1-f)*1) }
L2_COL <- unlist(lapply(levels(D$L1), function(L) {
  sub <- L2S[L1 == L]$L2; n <- max(length(sub), 1)
  setNames(sapply(seq_along(sub), function(i) shade(LIN_COL[[L]], 0.5 + 0.5*(i-1)/max(n-1,1))), sub)
}))
s <- D[, .N, by = .(stage, L2)][, frac := N / sum(N), by = stage]
s[, L2 := factor(L2, levels = names(L2_COL))]
p <- ggplot(s, aes(stage, frac, fill = L2)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = L2_COL, name = NULL, guide = guide_legend(ncol = 2, keyheight = unit(0.22, "cm"))) +
  scale_y_continuous(labels = scales::percent, expand = c(0.004, 0), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Fraction of nuclei", title = "L2 subtype composition by stage (39 types)") +
  theme_paper(8) + theme(legend.text = element_text(size = 4.4))
save_fig(p, "P2b_L2_composition_by_stage", 5.6, 4.2)
