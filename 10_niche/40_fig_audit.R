#!/usr/bin/env Rscript
# 40_fig_audit.R — 选药程序完备性审计图（只读）
#
# 输入: results/10_niche/target_reversal/null_smatched.tsv   (S 配平零模型, n=200)
#       results/10_niche/target_reversal/null_fullqref.tsv    (旧零模型, n=20, 含正侧)
# 输出: results/10_niche/target_reversal/fig_audit_stats.png
#
# 三个统计量并排，看哪个把 S1/S2/纯深度分开。

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
proj <- "/home/eto/luad_v2"
RD <- file.path(proj, "results/10_niche/target_reversal")
S <- fread(file.path(RD, "null_smatched.tsv"))
O <- fread(file.path(RD, "null_fullqref.tsv"))

OBS <- data.table(lab = c("S1 original", "S2 depth-matched", "pure-depth control"),
                  col = c("#E7298A", "#1B7837", "#7570B3"))
# (统计量面板, 零模型向量, 观测值, x 轴名, 标题)
panels <- list(
  list(nm = "one-sided count\ntau<=-90",        null = S$sig_le_m90, obs = c(0.0660, 0.0447, 0.0423),
       xl = "fraction of signatures with tau <= -90", ti = "A  pre-registered (best power)"),
  list(nm = "net asymmetry\n(neg - pos)",       null = O$rate_le_m90 - O$rate_ge_p90, obs = c(0.0660-0.0483, 0.0447-0.0388, NA),
       xl = "frac(tau<=-90) - frac(tau>=+90)   [old null, n=20]", ti = "B  CMap-literal (noisier)"),
  list(nm = "mean tau",                         null = S$mean_tau, obs = c(-2.09, -3.80, NA),
       xl = "mean tau over all signatures", ti = "C  continuous (weakest)"))

mk <- function(p) {
  dt <- data.table(v = p$null)
  lim <- range(c(dt$v, p$obs), na.rm = TRUE)
  g <- ggplot(dt, aes(v)) +
    geom_histogram(bins = 30, fill = "grey78", colour = "grey50", linewidth = .2) +
    geom_vline(xintercept = mean(dt$v), linetype = 2, colour = "grey15") +
    labs(x = p$xl, y = "null reps", title = p$ti) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 10))
  for (i in seq_along(p$obs)) if (!is.na(p$obs[i]))
    g <- g + annotate("segment", x = p$obs[i], xend = p$obs[i], y = Inf, yend = Inf,
                      colour = OBS$col[i]) +
             geom_vline(xintercept = p$obs[i], colour = OBS$col[i], linewidth = .9)
  g
}
gs <- lapply(panels, mk)
g <- (gs[[1]] | gs[[2]] | gs[[3]]) +
  plot_annotation(
    title = "Which statistic separates the versions?  Only the pre-registered one.",
    subtitle = "dashed = null mean;  pink = S1 original,  green = S2 depth-matched,  purple = pure-depth control",
    caption = "A: n=200, S-matched null.  B+C: see PREREG section 9.12.  None of the three puts S2 outside the random range.")
out <- file.path(RD, "fig_audit_stats.png")
ggsave(out, g, width = 11, height = 4.2, dpi = 150)
cat("wrote", out, "\n")
