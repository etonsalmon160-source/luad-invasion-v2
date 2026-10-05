#!/usr/bin/env Rscript
# 35_fig_depth_overlap.R — 深度倍数 vs 域签名稳定性 诊断图（只读）
#
# 输入: results/10_niche/target_reversal/depth_overlap_all_domains.tsv
#       (由 34_tr_depth_residual_sig.py 的残差估计量产出)
# 输出: results/10_niche/target_reversal/fig_depth_vs_overlap.png
#
# 论点: 域签名的"深度配平后存活率"由该域自身深度决定（7 域秩相关 -0.893），
#       D3 是唯一极端离群；两种独立估计量给出同一答案。

suppressMessages({library(data.table); library(ggplot2)})
proj <- "/home/eto/luad_v2"
T <- fread(file.path(proj, "results/10_niche/target_reversal/depth_overlap_all_domains.tsv"))
setnames(T, c("domain", "depth_ratio", "up", "down", "old"))
stopifnot(nrow(T) == 7L, all(T$up >= 0), all(T$up <= 150))

LAB <- c(D1="D1 airway epithelium", D2="D2 peribronchovascular",
         D3="D3 ECM/stroma", D4="D4 alveolar wall",
         D5="D5 AT2", D6="D6 vessel-wall SMC", D7="D7 lymphoid aggregate")
T[, lab := LAB[domain]]
T[, pct := 100 * up / 150]

rho <- cor(T$depth_ratio, T$up, method = "spearman")

g <- ggplot(T, aes(depth_ratio, pct)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey55", fill = "grey88", linewidth = .6) +
  geom_point(aes(colour = domain == "D3"), size = 4) +
  ggrepel::geom_text_repel(aes(label = lab), size = 3.2, min.segment.length = .2,
                           seed = 1, box.padding = .6) +
  scale_colour_manual(values = c(`FALSE` = "#2C7FB8", `TRUE` = "#B2182B"), guide = "none") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 5, 8)) +
  labs(x = "domain depth relative to the shallowest domain (x, log scale)",
       y = "% of top-150 up genes surviving depth matching",
       title = "A domain's signature stability is set by how deep that domain sits",
       subtitle = sprintf("7 archetypes; two independent estimators agree (D3 up 14->16); Spearman rho = %.2f; D3 (red) is the sole extreme outlier", rho),
       caption = "depth matching: per-slide cubic regression of each gene on log10(nUMI); domain-vs-rest re-ranking on residuals") +
  theme_bw(base_size = 11)

out <- file.path(proj, "results/10_niche/target_reversal/fig_depth_vs_overlap.png")
ggsave(out, g, width = 9, height = 5.6, dpi = 150)
cat("wrote", out, "\n")
cat(sprintf("rho = %.3f\n", rho))
