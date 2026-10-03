#!/usr/bin/env Rscript
# 19_signature_overlay.R —— M6 §21：把**功能签名**逐 spot 打上，再叠到 7 个域类型上
#
# 动机（用户 2026-10-03）：「没有什么耗竭细胞这些吗」＋「看看有没有各种耗竭细胞和这个纤维生态位」。
#   根因：RCTD 参考集（MP argmax 定的 39 类）**没有耗竭程序** ⇒ 耗竭 T 被并进 "CD8+ Naive T"
#   （实测簇 34：LAG3/TIGIT/HAVCR2/ENTPD1/GZMB/PRF1，却无 CCR7/SELL/TCF7/LEF1 ⇒ 错标）。
#
# 🔴 签名来源（逐个标注；标 [原文] 的是源论文自己用的，标 [文献] 的**需用户确认后才能正式化**）：
#   T_exhaust   [文献] Thommen & Schumacher 2018 Nat Med（NSCLC CXCL13+ 耗竭 CD8）/ Wherry & Kurachi 2015 Nat Immunol
#   T_cytotoxic [原文] 原文用 "cytotoxic CD8+ T"（小鼠治疗实验）；基因取经典细胞毒组
#   T_naive     [文献] 经典初始 T
#   Treg        [原文] 原文 "T regulatory cells ... increased with disease stage"
#   TLS         [原文] 原文 lymphoid/lymphoid aggregates 标记 = MS4A1, CXCL13
#   Mac_IL1B    [原文] 原文促炎生态位核心 = IL1B-high 巨噬 + IL1R1-high 上皮
#   Mac_LAM     [原文] 原文 myeloid 标记 = APOE, GPMNB(GPNMB)
#   Fib_ECM     [原文] 原文 fibroblast/stromal 标记 = COL1A1, ACTA2, DES（另加经典 ECM）
#   iCAF        [文献] Öhlund 2017 JEM / Elyada 2019 Cancer Discov（**胰腺来源，本癌种需确认**）
#   Myofibro    [文献+原文] ACTA2/TAGLN/DES（DES 出自原文）
#
# --slide <name>  单张：读 Visium 原始矩阵 → CP10K+log1p → 逐 spot 签名均分 → 原子写 d12_parts/
# --combine       并 d8 的 spot→archetype → 逐域签名分布 + 同切片配对

suppressMessages({ library(data.table); library(Matrix) })

ROOT  <- "/home/eto/luad_v2"
VIS   <- file.path(ROOT, "data/visium_spatial")
BN    <- file.path(ROOT, "results/10_niche/banksy")
OUT   <- file.path(ROOT, "results/10_niche/kstar_diag")
PARTS <- file.path(OUT, "d12_parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)

SIG <- list(
  T_exhaust   = c("PDCD1","LAG3","HAVCR2","TIGIT","ENTPD1","TOX","CXCL13","LAYN","TNFRSF9"),
  T_cytotoxic = c("GZMB","PRF1","NKG7","GNLY","GZMA","GZMK","KLRD1"),
  T_naive     = c("CCR7","SELL","TCF7","LEF1"),
  Treg        = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  TLS         = c("MS4A1","CXCL13","CD79A","CR2","CD19"),
  Mac_IL1B    = c("IL1B","IL1R1","NLRP3","CXCL8","S100A8","S100A9"),
  Mac_LAM     = c("APOE","GPNMB","TREM2","SPP1","LPL"),
  Fib_ECM     = c("COL1A1","COL1A2","COL3A1","FN1","POSTN","THBS2","COL11A1"),
  iCAF        = c("IL6","IL11","CXCL1","CXCL2","LIF","HAS1","PDGFRA","CFD"),
  Myofibro    = c("ACTA2","TAGLN","DES","MYH11")
)

.a <- commandArgs(trailingOnly = TRUE)
ga <- function(k) { i <- match(k, .a); if (is.na(i) || i == length(.a)) NULL else .a[i + 1] }
SLIDE <- ga("--slide"); COMBINE <- "--combine" %in% .a
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())
free_gb <- function() as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1",
                grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)[1])) / 1048576

one_slide <- function(s) {
  if (free_gb() < 20) stop("可用内存 < 20 GB ⇒ 硬停", call. = FALSE)
  mtx <- file.path(VIS, s, "filtered_feature_bc_matrix", "matrix.mtx.gz")
  fea <- fread(file.path(VIS, s, "filtered_feature_bc_matrix", "features.tsv.gz"),
               sep = "\t", header = FALSE)
  bar <- fread(file.path(VIS, s, "filtered_feature_bc_matrix", "barcodes.tsv.gz"), header = FALSE)
  M <- readMM(mtx)                       # genes x spots
  M <- as(M, "dgCMatrix")
  toc <- colSums(M); toc[toc == 0] <- 1
  Mn <- M %*% Diagonal(x = 1e4 / toc)
  Mn@x <- log1p(Mn@x)                    # CP10K + log1p（与项目其它臂一致）
  ## 🔴 第1列是 Ensembl ID，第2列才是基因符号（我踩过：用第1列 ⇒ 命中 0 个基因）
  rownames(Mn) <- fea[[2]]; colnames(Mn) <- bar[[1]]
  r <- data.table(slide = s, barcode = bar[[1]],
                  as.data.table(setNames(lapply(names(SIG), function(nm) {
                    g <- intersect(SIG[[nm]], rownames(Mn))
                    if (length(g)) Matrix::colMeans(Mn[g, , drop = FALSE]) else rep(NA_real_, ncol(Mn))
                  }), names(SIG))))
  r[, n_sig_genes_found := paste(sapply(names(SIG), function(nm) length(intersect(SIG[[nm]], rownames(Mn)))), collapse = "/")]
  r
}

if (!is.null(SLIDE)) {
  r <- one_slide(SLIDE)
  tmp <- file.path(PARTS, sprintf("%s.tsv.tmp", SLIDE)); fwrite(r, tmp, sep = "\t")
  if (!file.rename(tmp, file.path(PARTS, sprintf("%s.tsv", SLIDE)))) stop("原子改名失败", call. = FALSE)
  step("落盘 %s（%d spot；各签名命中基因 %s）", SLIDE, nrow(r), r$n_sig_genes_found[1])
  quit(save = "no", status = 0)
}
if (!COMBINE) stop("需要 --slide 或 --combine", call. = FALSE)

## ————————————— 汇总 —————————————
sl <- list.files(BN, pattern = "^GSM")
pf <- file.path(PARTS, sprintf("%s.tsv", sl)); have <- file.exists(pf)
step("D12 汇总：应有 %d / 已有 %d", length(sl), sum(have))
if (any(!have)) step("  ! 缺：%s", paste(sl[!have], collapse = " "))
S <- rbindlist(lapply(pf[have], fread), fill = TRUE)
D8 <- fread(file.path(OUT, "d8_spot_archetype_cnv.tsv.gz"))[, .(slide, barcode, domain, archetype, cf)]
S <- merge(S, D8, by = c("slide", "barcode"))
step("配对 %d / %d spot", nrow(S), sum(sapply(pf[have], function(f) nrow(fread(f, select = 1)))))
fwrite(S, file.path(OUT, "d12_spot_signatures.tsv.gz"), sep = "\t", compress = "gzip")

## 逐域均值
A <- S[, c(lapply(.SD, mean, na.rm = TRUE), list(n_spot = .N, n_slide = uniqueN(slide))),
       by = archetype, .SDcols = names(SIG)]
setorder(A, archetype)
fwrite(A, file.path(OUT, "d12_archetype_signatures.tsv"), sep = "\t")

## 同切片内配对：域内均值 − 该切片其余 spot 均值
OTH <- S[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = slide, .SDcols = names(SIG)]
setnames(OTH, names(SIG), paste0(names(SIG), "_o"))
P <- merge(S[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = .(slide, archetype), .SDcols = names(SIG)],
           OTH, by = "slide")
for (nm in names(SIG)) P[[paste0("d_", nm)]] <- P[[nm]] - P[[paste0(nm, "_o")]]
## 🔴 别写 `lapply(paste0("d_",...), function(x) median(P[[x]]))`：在 by= 里 P[[x]] 取的是**整列**
##    ⇒ 每组算出同一个全局中位、7 行完全一样（我踩过）。用 .SDcols + .SD
DCOL <- paste0("d_", names(SIG))
PS <- P[, c(lapply(.SD, function(x) median(x, na.rm = TRUE)), list(n_slide = .N)),
        by = archetype, .SDcols = DCOL]
POS <- P[, c(lapply(.SD, function(x) sum(x > 0, na.rm = TRUE))), by = archetype, .SDcols = DCOL]
setnames(POS, DCOL, paste0("pos_", names(SIG)))
PS <- merge(PS, POS, by = "archetype")
setorder(PS, archetype)
fwrite(PS, file.path(OUT, "d12_archetype_signatures_paired.tsv"), sep = "\t")

LBL <- c("D1 airway epithelium","D2 peribronchovascular","D3 immune infiltrate","D4 alveolar wall",
         "D5 AT2 epithelium","D6 vessel-wall SMC","D7 lymphoid aggregate")
con <- file(file.path(OUT, "d12_summary.txt"), "wt")
w <- function(...) cat(sprintf(...), file = con, append = TRUE)
w("=== §21 功能签名 × 7 个域类型（%s）===\n\n", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("签名来源：T_exhaust/T_naive/iCAF=[文献待确认]；其余见脚注；Myofibro 部分出自原文\n\n")
w("[1] 逐域签名均分（CP10K+log1p 均值）\n")
print(as.data.frame(A), file = con, digits = 4)
w("\n[2] 同切片内配对：域内均分 − 同切片其余 spot 均分（深度/批次自动抵消）\n")
print(as.data.frame(PS), file = con, digits = 4)
close(con); cat(readLines(file.path(OUT, "d12_summary.txt")), sep = "\n")
step("产物 → %s", OUT)
