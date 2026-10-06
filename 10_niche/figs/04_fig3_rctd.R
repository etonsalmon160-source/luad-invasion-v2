#!/usr/bin/env Rscript
# 04_fig3_rctd.R —— Fig 3：RCTD 生态位组分（39 亚型）随阶段变化
source("00_palette_theme.R")
DIR <- "/home/eto/luad_v2/results/08_spatial_deconv/rctd_d/per_slide"
fs <- list.files(DIR, pattern = "\\.weights\\.tsv\\.gz$", full.names = TRUE)
cat("RCTD 切片", length(fs), "\n")
L <- rbindlist(lapply(fs, function(f) {
  d <- fread(f)
  setnames(d, 1, "barcode")
  slide <- sub("\\.weights\\.tsv\\.gz$", "", basename(f))
  d[, slide := slide]
  m <- matrix(as.matrix(d[, -c(1, ncol(d)), with = FALSE]), ncol = ncol(d) - 2)
  colnames(m) <- setdiff(names(d), c("barcode", "slide"))
  data.table(slide = slide, as.data.table(m))
}))
st <- fread("/home/eto/luad_v2/results/10_niche/kstar_diag/d7_domain_assign.tsv",
            select = c("slide", "stage"))[!duplicated(slide)]
L <- merge(L, st, by = "slide")
L[, stage := factor(stage, levels = STAGE_LEVELS)]
CT <- setdiff(names(L), c("slide", "stage", "barcode"))
cat("细胞型", length(CT), "| spot", nrow(L), "\n")
S <- L[, lapply(.SD, mean), by = stage, .SDcols = CT]
S <- melt(S, id.vars = "stage", variable.name = "ct", value.name = "frac")
setorder(S, ct, stage)
# 按全体丰度排序（堆叠顺序）
ord <- L[, lapply(.SD, mean), .SDcols = CT][, sapply(.SD, mean)][order(-sapply(L[, lapply(.SD, mean), .SDcols = CT], c))]
S[, ct := factor(ct, levels = names(sort(sapply(L[, lapply(.SD, mean), .SDcols = CT], c), decreasing = TRUE)))]
# 39 型配色：按 RCTD 平均丰度用 viridis 分层
lv <- levels(S$ct)
COL39 <- setNames(viridisLite::viridis(length(lv), option = "turbo", end = 0.95), lv)
p <- ggplot(S, aes(stage, frac, fill = ct)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = COL39, name = "RCTD cell type (39)", guide = guide_legend(ncol = 1)) +
  scale_y_continuous(labels = scales::percent, expand = c(0.005, 0), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Mean RCTD weight", title = "Niche composition by stage (RCTD, 39 subtypes)") +
  theme_paper(8) + theme(legend.text = element_text(size = 4.6), legend.key.size = unit(0.26, "cm"))
save_fig(p, "Fig3_RCTD_composition_by_stage", 4.4, 6.0)
# 附：六谱系版（更易读）
map6 <- fread("/tmp/_rctd2lin.csv")
LIN6 <- intersect(LIN_LEVELS, unique(map6$L1))
S6 <- rbindlist(lapply(LIN6, function(LL) {
  cols <- intersect(map6[L1 == LL]$rctd, CT); LS <- L
  data.table(stage = LS$stage, lin = LL, frac = if (length(cols)) rowSums(LS[, ..cols]) else 0)
}), use.names = TRUE)
p6 <- ggplot(S6, aes(stage, frac, fill = lin)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = LIN_COL, name = "Lineage") +
  scale_y_continuous(labels = scales::percent, expand = c(0.005, 0), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Mean RCTD weight", title = "Niche composition by stage (6 lineages)") +
  theme_paper(8)
save_fig(p6, "Fig3b_RCTD_6lineage_by_stage", 3.4, 6.0)
