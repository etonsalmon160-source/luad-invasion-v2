#!/usr/bin/env Rscript
# 30_spatial_matrix.R —— 空间"矩阵图"：行=样本(五期各一)，列=H&E | RCTD | 耗竭T | TLS | ECM | 域类型
#   RCTD 列用**荧光配色**（黑底加性混合），其余列以 H&E 为底、spot 上色
source("00_palette_theme.R")
suppressMessages({library(png); library(grid)})
ROOT <- "/home/eto/luad_v2"

## ── 数据 ──
S <- fread(cmd = sprintf("zcat %s/results/10_niche/kstar_diag/d12_spot_signatures.tsv.gz", ROOT))
st <- fread(sprintf("%s/results/10_niche/kstar_diag/d7_domain_assign.tsv", ROOT),
            select = c("slide","stage"))[!duplicated(slide)]
S <- merge(S, st, by = "slide"); S[, stage := factor(norm_stage(stage), levels = STAGE_LEVELS)]
REP <- S[, .N, by = .(stage, slide)][order(-N)][, .SD[1], by = stage][order(stage)]
cat("代表切片:\n"); print(REP)

## L2 → 谱系（RCTD 名）
map6 <- fread("/tmp/_rctd2lin.csv")
REF <- setdiff(names(fread(sprintf("%s/results/08_spatial_deconv/rctd_d/per_slide/%s.weights.tsv.gz",
                                  ROOT, REP$slide[1]))), "")
REF <- REF[REF != "V1"]; REF <- REF[!grepl("^V[0-9]+$", REF)]
sub_lin <- map6$L1[match(REF, map6$rctd)]
names(sub_lin) <- REF
ok <- !is.na(sub_lin); REF <- REF[ok]; sub_lin <- sub_lin[ok]
FLUO <- make_fluo_palette(REF, sub_lin)
cat("RCTD 亚型", length(REF), "| 配色就绪\n")

load_slide <- function(SL) {
  w <- fread(sprintf("%s/results/08_spatial_deconv/rctd_d/per_slide/%s.weights.tsv.gz", ROOT, SL))
  setnames(w, 1, "barcode")
  pos <- fread(sprintf("%s/data/visium_spatial/%s/spatial/tissue_positions.csv", ROOT, SL), skip = 1,
               col.names = c("barcode","in_tissue","array_row","array_col","px","py"))
  sf  <- jsonlite::fromJSON(sprintf("%s/data/visium_spatial/%s/spatial/scalefactors_json.json", ROOT, SL))
  m <- match(w$barcode, pos$barcode)
  W <- as.matrix(w[, ..REF])
  list(bc = w$barcode, W = W, X = pos$py[m] * sf$tissue_hires_scalef,
       Y = -pos$px[m] * sf$tissue_hires_scalef,
       img = readPNG(sprintf("%s/data/visium_spatial/%s/spatial/tissue_hires_image.png", ROOT, SL)))
}

panel <- function(SL, sig, mode) {
  L <- load_slide(SL)
  sm <- S[slide == SL]; m <- match(L$bc, sm$barcode)
  H <- dim(L$img)[1]; Wd <- dim(L$img)[2]
  base <- ggplot() + annotation_raster(L$img, 0, Wd, -H, 0) +
    coord_fixed(xlim = c(0, Wd), ylim = c(-H, 0), expand = FALSE) +
    theme_void() + theme(plot.margin = margin(0,0,0,0))
  if (mode == "rctd") {
    rgbm <- L$W %*% t(col2rgb(FLUO[REF])/255)          # Σ w·color（加性）
    rgbm <- pmin(rgbm, 1)
    return(ggplot() + geom_point(aes(L$X, L$Y), colour = rgb(rgbm[,1],rgbm[,2],rgbm[,3]),
                                 size = 0.42, stroke = 0) +
      coord_fixed(xlim = c(0, Wd), ylim = c(-H, 0), expand = FALSE) +
      theme_void() + theme(plot.margin = margin(0,0,0,0),
                           panel.background = element_rect(fill = "black", colour = NA)))
  }
  if (mode == "arch") {
    d <- data.table(X = L$X, Y = L$Y, a = sm$archetype[m])
    d <- d[!is.na(a)]; d[, a := factor(paste0("D", a), levels = ARCH_LEVELS)]
    return(base + geom_point(data = d, aes(X, Y, colour = a), size = 0.42, stroke = 0) +
              scale_colour_manual(values = ARCH_COL, guide = "none"))
  }
  v <- sm[[sig]][m]
  d <- data.table(X = L$X, Y = L$Y, v = v)
  base + geom_point(data = d[!is.na(v)], aes(X, Y, colour = v), size = 0.42, stroke = 0) +
    scale_colour_gradientn(colours = c("#00000000","#FDE725","#F89441","#B02A30","#67000D"),
                           guide = "none", na.value = "#00000000")
}

## ── 组矩阵 ──
COLS <- list(
  list("he",   "H&E",                       NA),
  list("rctd", "RCTD deconvolution",        NA),
  list("sig",  "Exhausted T\n(9 markers)", "T_exhaust"),
  list("arch", "Niche domain (K*=7)",        NA)
)
## 耗竭签名（与 19_signature_overlay.R 同一套，出处 Thommen & Schumacher 2018 Nat Med）
EXH_MARKERS <- c("PDCD1","LAG3","HAVCR2","TIGIT","ENTPD1","TOX","CXCL13","LAYN","TNFRSF9")
cat("耗竭 marker:", paste(EXH_MARKERS, collapse=", "), "\n")

rows <- lapply(seq_len(nrow(REP)), function(i) {
  SL <- REP$slide[i]
  ps <- lapply(COLS, function(cc) {
    p <- if (cc[[1]] == "he") {
      img <- readPNG(sprintf("%s/data/visium_spatial/%s/spatial/tissue_hires_image.png", ROOT, SL))
      H <- dim(img)[1]; Wd <- dim(img)[2]
      ggplot() + annotation_raster(img, 0, Wd, -H, 0) + coord_fixed(xlim=c(0,Wd), ylim=c(-H,0), expand=FALSE) +
        theme_void() + theme(plot.margin = margin(0,0,0,0))
    } else if (cc[[1]] == "sig") panel(SL, cc[[3]], "sig") else panel(SL, cc[[3]], cc[[1]])
    p
  })
  wrap_plots(ps, nrow = 1)
})
fig <- wrap_plots(rows, ncol = 1) +
  plot_annotation(title = "Spatial matrix: one representative section per stage",
                  theme = theme(plot.title = element_text(face = "bold", size = 10)))
save_fig(fig, "P5_spatial_matrix", 7.2, 11.5)
