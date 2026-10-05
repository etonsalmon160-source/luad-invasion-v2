#!/usr/bin/env Rscript
# 38_fig_null_pcl.R — 零模型对照图（只读）
#
# 输入: results/10_niche/target_reversal/null_pcl.tsv  (200 次随机查询)
# 输出: results/10_niche/target_reversal/fig_null_pcl.png
#
# 论点: (左) MOA 类别层 —— S1/S2 的命中数都落在随机范围内 ⇒ 无候选
#       (右) 签名层     —— 只有未调整的 S1 出界

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
proj <- "/home/eto/luad_v2"
N <- fread(file.path(proj, "results/10_niche/target_reversal/null_pcl.tsv"))
stopifnot(nrow(N) == 200L)

# ── 左：MOA/PCL 层 ──
mu <- mean(N$n_pcl_le_m90); q95 <- quantile(N$n_pcl_le_m90, .95)
lp <- (1 + sum(N$n_pcl_le_m90 >= 11)) / 201
rp <- (1 + sum(N$n_pcl_le_m90 >= 18)) / 201
gL <- ggplot(N, aes(n_pcl_le_m90)) +
  geom_histogram(bins = 26, fill = "grey75", colour = "grey45", linewidth = .2) +
  geom_vline(xintercept = mu, linetype = 2, colour = "grey20") +
  geom_vline(xintercept = 11, colour = "#E7298A", linewidth = .9) +
  geom_vline(xintercept = 18, colour = "#1B7837", linewidth = .9) +
  annotate("text", x = mu, y = Inf, vjust = 1.6, hjust = -.1, size = 3, colour = "grey20",
           label = sprintf("random mean %.1f", mu)) +
  annotate("text", x = 11, y = Inf, vjust = 3.4, hjust = -.12, size = 3, colour = "#E7298A",
           label = sprintf("S1 = 11\np = %.3f", lp)) +
  annotate("text", x = 18, y = Inf, vjust = 5.6, hjust = -.12, size = 3, colour = "#1B7837",
           label = sprintf("S2 = 18\np = %.3f", rp)) +
  labs(x = "number of MOA classes with tau_PCL <= -90", y = "random queries (of 200)",
       title = "MOA-class layer: both versions inside the random range",
       subtitle = sprintf("random genes alone yield %.1f classes on average (max %d); S1 and S2 do not stand out",
                          mu, max(N$n_pcl_le_m90))) +
  theme_bw(base_size = 10.5)

# ── 右：签名层 ──
mu2 <- mean(N$sig_le_m90); sd2 <- sd(N$sig_le_m90)
OBS <- data.table(lab = c("S1 original", "S2 depth-matched", "pure-depth control"),
                  v = c(0.0660, 0.0447, 0.0423),
                  col = c("#E7298A", "#1B7837", "#7570B3"))
OBS[, p := sapply(v, function(o) (1 + sum(N$sig_le_m90 >= o)) / 201)]
gR <- ggplot(N, aes(100 * sig_le_m90)) +
  geom_histogram(bins = 26, fill = "grey75", colour = "grey45", linewidth = .2) +
  geom_vline(xintercept = 100 * mu2, linetype = 2, colour = "grey20") +
  geom_vline(data = OBS, aes(xintercept = 100 * v, colour = col), linewidth = .9) +
  scale_colour_identity() +
  annotate("text", x = 100 * mu2, y = Inf, vjust = 1.6, hjust = -.1, size = 3, colour = "grey20",
           label = sprintf("random mean %.2f%%", 100 * mu2)) +
  annotate("text", x = 6.60, y = Inf, vjust = 3.4, hjust = 1.08, size = 3, colour = "#E7298A",
           label = sprintf("S1 6.60%%\np = %.3f", OBS[1, p])) +
  annotate("text", x = 4.47, y = Inf, vjust = 5.6, hjust = 1.08, size = 3, colour = "#1B7837",
           label = sprintf("S2 4.47%%  p = %.2f", OBS[2, p])) +
  annotate("text", x = 4.23, y = Inf, vjust = 7.4, hjust = 1.08, size = 3, colour = "#7570B3",
           label = sprintf("depth-only 4.23%%  p = %.2f", OBS[3, p])) +
  labs(x = "% of (query x signature) pairs with tau <= -90", y = "random queries (of 200)",
       title = "Signature layer: only the un-adjusted S1 clears the null",
       subtitle = "and S1 is the version that is ~89% depth displacement (see fig_depth_vs_overlap)") +
  theme_bw(base_size = 10.5)

g <- gL / gR +
  plot_annotation(caption = "null = 200 random 300-gene queries through the identical pipeline, same Q_ref, same library")
out <- file.path(proj, "results/10_niche/target_reversal/fig_null_pcl.png")
ggsave(out, g, width = 9, height = 8, dpi = 150)
cat("wrote", out, "\n")
cat(sprintf("MOA layer  S1 p=%.3f  S2 p=%.3f\nsignature  S1 p=%.3f  S2 p=%.3f  depth p=%.3f\n",
            lp, rp, OBS[1, p], OBS[2, p], OBS[3, p]))
