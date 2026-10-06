#!/usr/bin/env Rscript
# 03_fig2_composition.R —— Fig 2：细胞组成随阶段变化（谱系 + L2 亚型）
source("00_palette_theme.R")
D <- fread("/tmp/_cache_cellmap.csv"); D[, stage := factor(stage, levels = STAGE_LEVELS)]
setnames(D, "A_frozen", "L1"); D[, L1 := factor(L1, levels = LIN_LEVELS)]

mk <- function(dt, key, cols, lab) {
  s <- dt[, .N, by = c("stage", key)]
  s[, frac := N / sum(N), by = stage]
  s[[key]] <- factor(s[[key]], levels = names(cols))
  ggplot(s, aes(stage, frac, fill = .data[[key]])) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = cols, name = lab) + scale_y_continuous(labels = scales::percent, expand = c(0.005, 0)) +
    labs(x = NULL, y = "Fraction of nuclei", title = lab) + theme_paper(8)
}
# L2 按谱系着色（同谱系用明暗区分）
l2map <- fread("/tmp/_l2_map.csv")
setnames(l2map, "A_frozen", "L1")
# ⚠️ L2 名字跨谱系重名（39 个名字 / 67 个谱系×亚型组合）⇒ 每个名字只取**主谱系**着色
L2S <- D[, .N, by = .(L2, L1)][order(-N)][, .SD[1], by = L2][, .(L2, L1)]; setorder(L2S, L1, L2)
shade <- function(hex, f) { r <- col2rgb(hex)/255; rgb(r[1]*f + (1-f)*1, r[2]*f + (1-f)*1, r[3]*f + (1-f)*1) }
L2_COL <- unlist(lapply(LIN_LEVELS, function(L) {
  sub <- L2S[L1 == L]$L2; n <- max(length(sub), 1)
  setNames(sapply(seq_along(sub), function(i) shade(LIN_COL[[L]], 0.55 + 0.45 * (i-1)/max(n-1,1))), sub)
}))
pA <- mk(D, "L1", LIN_COL, "Lineage")
pB <- mk(D, "L2", L2_COL, "L2 subtype (39)")
fig <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11), legend.text = element_text(size = 5))
save_fig(fig, "Fig2_composition_by_stage", 9.6, 4.6)
