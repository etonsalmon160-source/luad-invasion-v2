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

## 图例画在坐标框内的统一样式（半透明底，少挡数据）
LEG_IN <- theme(legend.position = c(0.015, 0.015),
                legend.justification = c(0, 0),
                legend.background = element_rect(fill = alpha("white", 0.80),
                                                 colour = "grey55", linewidth = 0.2),
                legend.margin = margin(2, 3, 2, 3),
                legend.key.height = unit(0.30, "cm"))

## ── P1a · UMAP by lineage ──
p <- base + pts + aes(colour = L1) +
  scale_colour_manual(values = LIN_COL, name = NULL) +
  guides(colour = LEG) + labs(title = "Single-cell atlas: lineage") + LEG_IN
save_fig(p, "P1a_UMAP_lineage", 4.2, 3.9)

## ── P1b · UMAP by stage ──
p <- base + pts + aes(colour = stage) +
  scale_colour_manual(values = STAGE_COL, name = NULL) +
  guides(colour = LEG) + labs(title = "Single-cell atlas: stage") + LEG_IN
save_fig(p, "P1b_UMAP_stage", 4.2, 3.9)

## ── P1c · UMAP split by stage（保留谱系配色，与参考图一致）──
p <- base + pts + aes(colour = L1) + scale_colour_manual(values = LIN_COL, guide = "none") +
  facet_wrap(~stage, nrow = 1) + labs(title = "Single-cell atlas by stage (lineage colours)") +
  theme(panel.spacing = unit(1.5, "pt"),
        strip.text = element_text(colour = "white", face = "bold", size = 7.5),
        strip.background = element_rect(fill = "grey30", linewidth = 0))
save_fig(p, "P1c_UMAP_stage_split", 8.0, 2.1)

## ── P1d · UMAP split by patient ──
D[, pat := factor(patient_id)]
p <- base + pts + aes(colour = L1) + scale_colour_manual(values = LIN_COL, guide = "none") +
  facet_wrap(~pat, nrow = 4) + labs(title = "Single-cell atlas by patient (lineage colours)") +
  theme(panel.spacing = unit(0.6, "pt"),
        strip.text = element_text(colour = "white", face = "bold", size = 5),
        strip.background = element_rect(fill = "grey35", linewidth = 0))
save_fig(p, "P1d_UMAP_patient_split", 7.2, 4.6)

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
## 39 个亚型色：在该谱系色的色相附近展开 + 保持高饱和、明度从亮到中暗
## （旧写法是往白里混 ⇒ 最浅那档只有 50% 颜色，UMAP 上发白）
l2_pal <- function(hex, n) {
  h0 <- rgb2hsv(col2rgb(substr(hex, 1, 7)))   # pal_npg 是 #RRGGBBAA，col2rgb 只吃 7 位
  hh <- h0["h", 1]; ss <- h0["s", 1]          # 矩阵要按 [行, 列] 取，单下标是线性索引
  span <- if (n == 1) 0 else 0.052
  list(h = (hh + if (n == 1) 0 else seq(-span, span, length.out = n)) %% 1,
       s = rep(min(0.95, max(0.62, ss)), n),
       v = seq(0.98, 0.68, length.out = n))
}
L2_COL <- unlist(lapply(levels(D$L1), function(L) {
  sub <- L2S[L1 == L]$L2; n <- max(length(sub), 1)
  pp <- l2_pal(LIN_COL[[L]], n)
  setNames(hsv(pp$h, pp$s, pp$v), sub)
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

## ── P2c · UMAP by 39 L2 subtypes（图例放在坐标框内）──
D[, L2f := factor(L2, levels = names(L2_COL))]
p <- base + geom_point(size = 0.28, stroke = 0, shape = 16) + aes(colour = L2f) +
  scale_colour_manual(values = L2_COL, name = NULL) +
  guides(colour = guide_legend(ncol = 3, keyheight = unit(0.17, "cm"),
                               override.aes = list(size = 1.4, alpha = 1))) +
  labs(title = "Single-cell atlas: 39 L2 subtypes") +
  theme(legend.position = c(0.010, 0.990),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 3.5),
        legend.key.size = unit(0.19, "cm"),
        legend.key.height = unit(0.19, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.80),
                                         colour = "grey55", linewidth = 0.2),
        legend.margin = margin(2, 3, 2, 3))
save_fig(p, "P2c_UMAP_L2_subtype", 7.6, 6.0)
