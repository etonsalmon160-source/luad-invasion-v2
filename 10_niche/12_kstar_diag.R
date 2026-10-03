#!/usr/bin/env Rscript
# 12_kstar_diag.R —— M6 §14 诊断 D1 / D2a / D3（**全程只读已落盘产物，不重算 BANKSY**）
#
# 预注册：NICHE_PREREG.md §14.3（2026-10-02 已签）。本脚本**不产生任何主结果**。
#   D1  = 选中主档（k_geom=18 / λ=0.2 / res=0.5）的逐切片域大小分布
#   D2a = ARI 随"域数"的变化（用已落盘的 24 档 × 56 张 stats）
#   D3  = K 稳定性曲线（consensus/<agf>/k_stability.tsv）
#
# 跑法：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 10_niche/12_kstar_diag.R

suppressMessages({ library(data.table); library(ggplot2) })

ROOT  <- "/home/eto/luad_v2"
NICHE <- file.path(ROOT, "results/10_niche")
BANKSY<- file.path(NICHE, "banksy")
CONS  <- file.path(NICHE, "consensus")
OUT   <- file.path(NICHE, "kstar_diag")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

AGFS   <- c("agfT", "agfF")
MAIN   <- list(k_geom = 18L, lambda = 0.2, resolution = 0.5)   # §13 主档

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())

## —————————————————————————————————————————————————————————————
## 读全部 stats
## —————————————————————————————————————————————————————————————
step("==== §14 诊断 D1/D2a/D3（只读）====")
slides <- list.files(BANKSY, pattern = "^GSM")
step("切片目录 %d 个", length(slides))

ST <- rbindlist(lapply(AGFS, function(agf) {
  rbindlist(lapply(slides, function(s) {
    f <- file.path(BANKSY, s, sprintf("stats_%s.tsv", agf))
    if (!file.exists(f)) { step("  ! 缺 %s/%s", s, agf); return(NULL) }
    d <- fread(f, sep = "\t")
    d[, agf := agf]
    d
  }))
}), fill = TRUE)
ST[, slide := factor(slide)]
step("读入 stats %d 行（%d 张 × %d 档 × %d agf）", nrow(ST),
     uniqueN(ST$slide), uniqueN(ST[, .(k_geom, lambda, resolution)]), uniqueN(ST$agf))

## —————————————————————————————————————————————————————————————
## D1 —— 主档逐切片域大小
## —————————————————————————————————————————————————————————————
step("---- D1 主档（k_geom=%d λ=%.1f res=%.1f）域大小 ----",
     MAIN$k_geom, MAIN$lambda, MAIN$resolution)
D1 <- ST[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution]
D1 <- D1[order(agf, slide)]
fwrite(D1[, .(agf, slide, n_domain, n_singleton, n_lt10, size_min, size_med, size_max,
              ari_mean, ari_min, coh_mean, null_coh_q95)],
       file.path(OUT, "d1_maincell_domain_sizes.tsv"), sep = "\t")

for (a in AGFS) {
  x <- D1[agf == a]
  step("  %s：域数 中位 %g（范围 %d–%d）；单/双 spot 域 中位 %g；<10 spot 域 中位 %g；最小域 中位 %g",
       a, median(x$n_domain), min(x$n_domain), max(x$n_domain),
       median(x$n_singleton), median(x$n_lt10), median(x$size_min))
  step("      跨种子 ARI 中位 %.3f（最低侧中位 %.3f）；连贯 %.3f vs 零分布 q95 中位 %.3f",
       median(x$ari_mean), median(x$ari_min), median(x$coh_mean), median(x$null_coh_q95))
}

## —————————————————————————————————————————————————————————————
## D2a —— ARI vs 域数（只用已落盘统计）
## —————————————————————————————————————————————————————————————
step("---- D2a ARI vs 域数（24 档 × 56 张）----")
D2 <- ST[, .(agf, slide, k_geom, lambda, resolution, n_domain, ari_mean, ari_min, coh_mean)]
fwrite(D2, file.path(OUT, "d2a_ari_vs_ndomain.tsv"), sep = "\t")

## 逐切片：n_domain 与 ari_mean 的秩相关（负值 = 域越少越稳）
sp <- D2[, {
  ok <- is.finite(n_domain) & is.finite(ari_mean)
  if (sum(ok) < 4) .(rho = NA_real_, n = sum(ok))
  else .(rho = suppressWarnings(cor(n_domain[ok], ari_mean[ok], method = "spearman")), n = sum(ok))
}, by = .(agf, slide)]
step("  逐切片 Spearman(域数, ARI) 中位：%s",
     paste(sprintf("%s %.3f", AGFS,
                   sapply(AGFS, function(a) median(sp[agf == a]$rho, na.rm = TRUE))), collapse = " / "))
step("    （负 = 域越少越稳；接近 0 = 与域数无关）")

## 同一 λ 下 res 0.8→0.6→0.5 的方向
step("  同一 λ 下 res 降低（0.8→0.5）时 ARI 的配对变化（中位 ΔARI）：")
## 🔴 (slide, resolution) 不唯一（还有 k_geom 两档）⇒ key 必须带 k_geom，否则 dcast 退化成计数
for (a in AGFS) for (l in c(0, 0.2, 0.5, 0.8)) for (kg in c(6L, 18L)) {
  sub <- D2[agf == a & lambda == l & k_geom == kg]
  w  <- dcast(sub, slide ~ resolution, value.var = "ari_mean")
  nd <- dcast(sub, slide ~ resolution, value.var = "n_domain")
  if (all(c("0.5", "0.8") %in% names(w))) {
    step("    %s λ=%.1f k_geom=%2d：ΔARI(res0.5-res0.8) 中位 %+.3f；域数 %g→%g",
         a, l, kg, median(w[["0.5"]] - w[["0.8"]], na.rm = TRUE),
         median(nd[["0.8"]], na.rm = TRUE), median(nd[["0.5"]], na.rm = TRUE))
  }
}
## 总体：最低档 vs 最高档域数bin 的 ARI
step("  按域数分箱的跨种子 ARI 中位（全 24 档 pooled）：")
D2[, ndbin := cut(n_domain, breaks = c(0, 3, 5, 8, 12, 20, 50, Inf),
                  labels = c("1-3", "4-5", "6-8", "9-12", "13-20", "21-50", ">50"))]
print(D2[, .(n = .N, ari_med = round(median(ari_mean, na.rm = TRUE), 3),
             nd_med = as.numeric(median(n_domain))), by = .(agf, ndbin)][order(agf, ndbin)])

## —————————————————————————————————————————————————————————————
## D3 —— K 稳定性曲线
## —————————————————————————————————————————————————————————————
step("---- D3 K 稳定性曲线 ----")
D3 <- rbindlist(lapply(AGFS, function(a) {
  f <- file.path(CONS, a, "k_stability.tsv")
  if (!file.exists(f)) { step("  ! 缺 %s", f); return(NULL) }
  d <- fread(f); d[, agf := a]; d
}), fill = TRUE)
fwrite(D3, file.path(OUT, "d3_kstability.tsv"), sep = "\t")
K_STAB_MIN <- 0.90
step("  阈值 %.2f：%s", K_STAB_MIN,
     paste(sprintf("%s 最高 K=%d（%.3f）", AGFS,
                   sapply(AGFS, function(a) D3[agf == a][which.max(stab_mean)]$K),
                   sapply(AGFS, function(a) max(D3[agf == a]$stab_mean, na.rm = TRUE))),
           collapse = " / "))

## —————————————————————————————————————————————————————————————
## 出图（全 English 标签，无 CJK 字体）
## —————————————————————————————————————————————————————————————
step("---- 出图 ----")
TH <- theme_bw(base_size = 12) + theme(legend.position = "bottom")

## 图 1：ARI vs 域数
g1 <- ggplot(D2, aes(x = n_domain, y = ari_mean)) +
  geom_point(aes(colour = factor(lambda)), alpha = 0.45, size = 1.3) +
  geom_hline(yintercept = K_STAB_MIN, linetype = 2, colour = "red") +
  geom_smooth(aes(colour = factor(lambda)), method = "loess", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ agf, nrow = 1) +
  scale_colour_brewer(palette = "Dark2", name = "lambda") +
  labs(x = "number of domains per slide", y = "cross-seed ARI (mean)",
       title = "D2a: is instability driven by too many domains?",
       subtitle = "all 24 grid cells x 56 slides; red dashed = 0.90 gate") + TH
ggsave(file.path(OUT, "fig_d2a_ari_vs_ndomain.png"), g1, width = 10, height = 5, dpi = 150)

## 图 2：K 稳定性曲线
g2 <- ggplot(D3, aes(x = K, y = stab_mean, colour = agf)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  geom_hline(yintercept = K_STAB_MIN, linetype = 2, colour = "red") +
  geom_point(data = D3[, .SD[which.max(stab_mean)], by = agf], shape = 4, size = 4, stroke = 1.2) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "consensus K", y = "bootstrap ARI (mean)",
       title = "D3: consensus K stability vs the 0.90 threshold",
       subtitle = "cross = argmax (the silent-fallback value actually selected)") + TH
ggsave(file.path(OUT, "fig_d3_kstability.png"), g2, width = 7.5, height = 5, dpi = 150)

## 图 3：主档域数分布
g3 <- ggplot(D1, aes(x = n_domain, fill = agf)) +
  geom_histogram(binwidth = 1, position = "dodge", colour = "grey30") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "domains per slide (main cell)", y = "slides",
       title = "D1: domains per slide at the selected main cell") + TH
ggsave(file.path(OUT, "fig_d1_domains_per_slide.png"), g3, width = 7.5, height = 4.5, dpi = 150)

step("产物 → %s", OUT)
step("==== 诊断 D1/D2a/D3 结束 ====")
