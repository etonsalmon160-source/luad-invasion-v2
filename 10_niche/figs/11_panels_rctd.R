#!/usr/bin/env Rscript
# 11_panels_rctd.R —— RCTD 生态位组分（39 亚型 / 6 谱系），分开出、全英文
source("00_palette_theme.R")
DIR <- "/home/eto/luad_v2/results/08_spatial_deconv/rctd_d/per_slide"
fs <- list.files(DIR, pattern = "\\.weights\\.tsv\\.gz$", full.names = TRUE)
L <- rbindlist(lapply(fs, function(f) {
  d <- fread(f); setnames(d, 1, "barcode")
  m <- as.matrix(d[, -1, with = FALSE]); ct <- names(d)[-1]
  data.table(slide = sub("\\.weights\\.tsv\\.gz$", "", basename(f)), as.data.table(m)[, setnames(.SD, ct)])
}))
st <- fread("/home/eto/luad_v2/results/10_niche/kstar_diag/d7_domain_assign.tsv",
            select = c("slide", "stage"))[!duplicated(slide)]
L <- merge(L, st, by = "slide"); L[, stage := factor(norm_stage(stage), levels = STAGE_LEVELS)]
CT <- setdiff(names(L), c("slide", "stage", "barcode"))
cat("slides", length(fs), "| spots", nrow(L), "| cell types", length(CT), "\n")

## ── P3a · 39 subtypes ──
S <- melt(L[, lapply(.SD, mean), by = stage, .SDcols = CT], id.vars = "stage",
          variable.name = "ct", value.name = "frac")
ab <- sapply(L[, lapply(.SD, mean), .SDcols = CT], c)
S[, ct := factor(ct, levels = names(sort(ab, decreasing = TRUE)))]
COL39 <- setNames(viridisLite::viridis(length(levels(S$ct)), option = "turbo", end = 0.95), levels(S$ct))
p <- ggplot(S, aes(stage, frac, fill = ct)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = COL39, name = "RCTD cell type", guide = guide_legend(ncol = 1, keyheight = unit(0.22, "cm"))) +
  scale_y_continuous(labels = scales::percent, expand = c(0.004, 0), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Mean RCTD weight", title = "Niche composition by stage (39 subtypes)") +
  theme_paper(8) + theme(legend.text = element_text(size = 4.4))
save_fig(p, "P3a_RCTD39_composition_by_stage", 4.6, 6.2)

## ── P3b · 6 lineages ──
map6 <- fread("/tmp/_rctd2lin.csv")
setnames(map6, "L1", "lin"); map6[, lin := LIN_EN[lin]]
S6 <- rbindlist(lapply(unique(map6$lin), function(LL) {
  cols <- intersect(map6[lin == LL]$rctd, CT)
  data.table(stage = L$stage, lin = LL, frac = if (length(cols)) rowSums(L[, ..cols]) else 0)
}))
S6 <- S6[, .(frac = mean(frac)), by = .(stage, lin)]
S6[, lin := factor(lin, levels = LIN_EN[LIN_LEVELS])]
p <- ggplot(S6, aes(stage, frac, fill = lin)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = LIN_COL, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0.004, 0), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Mean RCTD weight", title = "Niche composition by stage (6 lineages)") +
  theme_paper(8) + guides(fill = guide_legend(keyheight = unit(0.3, "cm")))
save_fig(p, "P3b_RCTD6_composition_by_stage", 3.3, 3.1)
