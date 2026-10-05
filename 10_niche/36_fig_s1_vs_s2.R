#!/usr/bin/env Rscript
# 36_fig_s1_vs_s2.R — S1(原签名) vs S2(深度配平签名) 的 PCL 层面并排图（只读）
#
# 输入: results/10_niche/target_reversal/pcl_S1_vs_S2.tsv   (moa, tau_S1, tau_S2, d)
# 输出: results/10_niche/target_reversal/fig_S1_vs_S2_pcl.png
#
# 论点: 两版回答两个问题（S1=总效应含密度；S2=去密度残余）。
#       本图不替读者选，只把两版并排画出来。

suppressMessages({library(data.table); library(ggplot2)})
proj <- "/home/eto/luad_v2"
D <- fread(file.path(proj, "results/10_niche/target_reversal/pcl_S1_vs_S2.tsv"))
stopifnot(all(c("moa","tau_S1","tau_S2","d") %in% names(D)))

D[, hit := fifelse(tau_S1 <= -90 & tau_S2 <= -90, "both",
             fifelse(tau_S1 <= -90, "S1 only",
              fifelse(tau_S2 <= -90, "S2 only", "neither")))]
D[, hit := factor(hit, levels = c("both","S1 only","S2 only","neither"))]

COL <- c(both = "#762A83", `S1 only` = "#E7298A", `S2 only` = "#1B7837", neither = "grey80")

n1 <- D[hit == "S1 only" | hit == "both", .N]
n2 <- D[hit == "S2 only" | hit == "both", .N]

g <- ggplot(D, aes(tau_S1, tau_S2)) +
  geom_hline(yintercept = -90, linetype = 2, colour = "grey40") +
  geom_vline(xintercept = -90, linetype = 2, colour = "grey40") +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = .4) +
  geom_point(aes(colour = hit), size = 1.8, alpha = .8) +
  scale_colour_manual(values = COL, name = "MOA class",
    labels = c(sprintf("tau<=-90 in both (n=%d)", D[hit=="both", .N]),
               sprintf("S1 (original) only (n=%d)", D[hit=="S1 only", .N]),
               sprintf("S2 (depth-matched) only (n=%d)", D[hit=="S2 only", .N]),
               "neither")) +
  annotate("label", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5, size = 3, label.size = NA,
           label = sprintf("S1 total hits %d   |   S2 total hits %d   |   overlap %d",
                           n1, n2, D[hit == "both", .N])) +
  labs(x = "tau_PCL  —  S1 original D3 signature (total effect, includes density)",
       y = "tau_PCL  —  S2 depth-matched D3 signature (density residual)",
       title = "Two versions answer two different questions — reported side by side, no winner chosen",
       caption = "dashed lines = tau -90 (strong-score threshold); the two signatures' difference is the density-carried component") +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 100)) +
  theme_bw(base_size = 11)

out <- file.path(proj, "results/10_niche/target_reversal/fig_S1_vs_S2_pcl.png")
ggsave(out, g, width = 8.4, height = 6.4, dpi = 150)
cat("wrote", out, "\n")
cat(sprintf("S1 hits %d | S2 hits %d | overlap %d | total MOA %d\n", n1, n2, D[hit=="both",.N], nrow(D)))
