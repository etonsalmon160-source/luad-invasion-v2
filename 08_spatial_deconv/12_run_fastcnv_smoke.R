#!/usr/bin/env Rscript
# 空转（Visium）CNV 臂 —— fastCNV **单患者体检**（SPATIAL_CNV_PREREG.md §15）
#
# 本脚本**不产生任何恶性判定**，**不替换** §3 的 infercnv 复现臂。
# 它只回答 §15.3 的三个工具层面问题：
#   V1 速度       —— 与 infercnv 同患者同输入（P4，34,333 上皮 spot），直接可比
#   V2 自平性     —— V2a 池化三张良性作参考（§11.2 C 主锚）
#                    V2b 只用 Normal 一张作参考 ⇒ AAH/AAH-1 变成**留出的**良性观测
#                        （本条是**唯一能证伪**该工具的路子）
#   V3 不崩不空   —— 无 NA、窗口非空、基因账逐项记账
#
# 🔴 本脚本产出的所有数**不得**用于 SC0–SC4 任何一门；「恶性」判据仍未签字（§13.3 第 1 项）。
#
# 跑法（**只挂 fastcnv 这个库**，别和 infercnv 的库混，两边的 Seurat 版本不同）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 08_spatial_deconv/12_run_fastcnv_smoke.R [--patient P4]

suppressMessages({
  library(fastCNV)
  library(Seurat)
  library(Matrix)
})

## —————————————————————————————————————————————————————————————
## 一、参数（逐项显式，法则 3.1；来源见 §15.4）
## —————————————————————————————————————————————————————————————
CAL_REFERENCE <- "a"        # [已签 §11 第 1 项] 粗版 6 谱系
CAL_EPI       <- "上皮"     # [已签 §11 第 2 项] 单类目
CAL_RULE      <- "argmax"   # [已签 §12.1] 行内最大类 == 上皮，零阈值

F_PREPARE_COUNTS <- FALSE   # 🔴 [必须显式覆盖工具默认] 见 §15.4：默认在 spot<60,000 时聚合，
                            #    P4 上皮 34,333 < 60,000 ⇒ 不关就把逐 spot 分辨率抹掉
F_AGGREG_FACTOR  <- 999     # <1000 ⇒ 工具自己也拒绝聚合（双保险）
F_WINDOW_SIZE    <- 150     # [工具默认，显式写出] ⚠️ 与 infercnv 的 101 不同
F_WINDOW_STEP    <- 10      # [工具默认，显式写出]
F_TOPN_GENES     <- 7000    # [工具默认，显式写出]
F_THRESH_PCT     <- 0.01    # [工具默认，显式写出] 背景噪声分位带
F_POOLED_REF     <- TRUE    # [工具默认] 多切片汇池参考（§11.2 借锚所需的机制）
F_GET_ARM        <- TRUE    # 逐臂产物（§13.3 第 1 项的候选产物）
F_GET_CLUSTERS   <- FALSE   # 体检不做聚类
F_DO_PLOT        <- FALSE   # 体检不做绘图（与 §13 关绘图同一理由）
F_ASSAY          <- "Spatial"  # 显式指定，避免二次运行时拿到刚加的 genomicScores

ROOT   <- "/home/eto/luad_v2"
RES    <- file.path(ROOT, "results/08_spatial_deconv")
VISIUM <- file.path(ROOT, "data/visium_spatial")
OUTD   <- file.path(RES, "spatial_cnv/fastcnv_smoke_P4")

## —————————————————————————————————————————————————————————————
## 二、工具函数
## —————————————————————————————————————————————————————————————
args    <- commandArgs(trailingOnly = TRUE)
patient <- { i <- which(args == "--patient"); if (length(i)) args[i + 1] else "P4" }

dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()
RUN_MIN <- list()   # 每次 CNVCalling 的墙钟（V1 分段用）
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

peak_gb <- function() {
  s <- readLines("/proc/self/status")
  hw <- grep("^VmHWM:", s, value = TRUE)
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", hw)) / 1024^2
}
step <- function(fmt, ...) cat(sprintf(
  paste0("[%s | %6.1f min | 峰值 %6.2f GB] ", fmt),
  format(Sys.time(), "%H:%M:%S"), elapsed(t_start), peak_gb(), ...), "\n")

read_gz <- function(path, fn) {
  con <- gzfile(path); on.exit(try(close(con), silent = TRUE)); fn(con)
}

## 无阈值的分离度：AUC = P(x > y)，并列各算 0.5。**只报**，不设门槛（§15.3）
auc_gt <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  r <- rank(c(x, y))
  (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y))
}

stage_of <- function(s) sub("^[^_]*_[^_]*_", "", s)
is_benign <- function(s) grepl("^Normal|^AAH", stage_of(s))

step("==== fastCNV 体检开始：患者 %s ====", patient)
step("工具 fastCNV %s / Seurat %s（**预印本方法，未经同行评审**，§15.2）",
     as.character(packageVersion("fastCNV")), as.character(packageVersion("Seurat")))

## —————————————————————————————————————————————————————————————
## 三、spot 掩码 与 RCTD 权重（口径 argmax；与 11_ 脚本逐字同源）
## —————————————————————————————————————————————————————————————
sm <- read.delim(gzfile(file.path(RES, "spot_mask.tsv.gz")), stringsAsFactors = FALSE)
sm <- sm[as.character(sm$pass) %in% c("TRUE", "true", "1"), ]
sm <- sm[grepl(paste0("_", patient, "_"), sm$slide), ]
stopifnot(nrow(sm) > 0)
slides <- sort(unique(sm$slide))
step("患者 %s 的切片 %d 张：%s", patient, length(slides), paste(slides, collapse = ", "))

epi_bc <- list()
for (s in slides) {
  wf <- file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s))
  stopifnot(file.exists(wf))
  w <- read.delim(gzfile(wf), row.names = 1, check.names = FALSE)
  stopifnot(CAL_EPI %in% colnames(w), ncol(w) == 6L)
  keep <- intersect(rownames(w), sm$barcode[sm$slide == s])
  W    <- as.matrix(w[keep, , drop = FALSE])
  mx   <- do.call(pmax, as.data.frame(W))
  n_at <- rowSums(W >= mx)
  am   <- colnames(W)[max.col(W, ties.method = "first")]
  tie  <- n_at > 1L
  epi_bc[[s]] <- keep[am == CAL_EPI & !tie]
  step("  %-26s QCpass %6d  上皮 %6d  并列剔除 %d  (%s%s)",
       s, nrow(W), length(epi_bc[[s]]), sum(tie), stage_of(s),
       if (is_benign(s)) " · 良性" else "")
}
n_epi <- sum(lengths(epi_bc))
step("上皮 spot 合计 %d（与 11_ 脚本同口径，可直接对比 infercnv 那 990 min）", n_epi)

## —————————————————————————————————————————————————————————————
## 四、读表达矩阵（只取 Gene Expression；双探针按 symbol 求和；只留上皮 spot）
##    合并成**一个** Seurat 对象 ⇒ 走单对象路径（见 15.8 第 6 项：汇池路径有 bug）
## —————————————————————————————————————————————————————————————
t_read <- Sys.time()
mats <- list(); lab_a <- character(0); lab_b <- character(0); cell_slide <- character(0)
for (s in slides) {
  d      <- file.path(VISIUM, s, "filtered_feature_bc_matrix")
  ft_all <- read_gz(file.path(d, "features.tsv.gz"), function(con)
    read.delim(con, header = FALSE, col.names = c("id", "sym", "type"), stringsAsFactors = FALSE))
  ab <- which(ft_all$type == "Antibody Capture")
  stopifnot(length(ab) %in% c(0L, 35L))
  bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
  m  <- readMM(file.path(d, "matrix.mtx.gz"))
  stopifnot(nrow(m) == nrow(ft_all), ncol(m) == length(bc))
  ## ⚠️ 陷阱：`x[-integer(0), ]` 返回**零行**（不是全部行）⇒ 必须 if 包住
  if (length(ab)) { m <- m[-ab, , drop = FALSE]; ft <- ft_all[-ab, , drop = FALSE] } else { ft <- ft_all }
  stopifnot(nrow(ft) == 18085L, all(ft$type == "Gene Expression"))

  use <- match(epi_bc[[s]], bc, nomatch = 0L)
  stopifnot(!any(use == 0L), !anyDuplicated(use))
  m  <- m[, use, drop = FALSE]
  bt <- bc[use]

  ## 双探针（HSPA14 / TBCE / TMSB15B）按 symbol 求和取一次
  g   <- factor(ft$sym, levels = unique(ft$sym))
  Agg <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1,
                      dims = c(nlevels(g), length(g)))
  rownames(Agg) <- levels(g)
  m <- Agg %*% m
  rownames(m) <- levels(g)
  ## Seurat 对 barcode 字符有限制 ⇒ 用 "_" 连接，不带 "|"
  colnames(m) <- paste(s, bt, sep = "_")

  ## 两套参考标签：V2a 池化三张良性；V2b 只用 Normal 一张（其余成为留出观测）
  ga <- if (is_benign(s)) "benign_epi" else "tumor_epi"
  gb <- if (stage_of(s) == "Normal") "ref_normal" else "heldout"
  mats[[s]]      <- m
  cell_slide     <- c(cell_slide, setNames(rep(s, ncol(m)), colnames(m)))
  lab_a          <- c(lab_a, rep(ga, ncol(m)))
  lab_b          <- c(lab_b, rep(gb, ncol(m)))
  step("  %-26s 上皮 %6d  标签 g_a=%-10s g_b=%s", s, ncol(m), ga, gb)
  rm(m, Agg); invisible(gc())
}
mat <- do.call(cbind, mats); rm(mats); invisible(gc())
stopifnot(ncol(mat) == length(lab_a), identical(colnames(mat), names(cell_slide)))
obj <- CreateSeuratObject(counts = mat, assay = F_ASSAY,
                          meta.data = data.frame(g_a = lab_a, g_b = lab_b,
                                                 slide = unname(cell_slide),
                                                 row.names = colnames(mat)))
min_load <- elapsed(t_read)   # 读矩阵＋建对象的**已用时**（此前误写为「到脚本结尾的用时」，已改）
step("读入并合并成单对象：%d 基因 x %d spot（用时 %.1f min，峰值 %.2f GB）",
     nrow(obj), ncol(obj), min_load, peak_gb())

## —————————————————————————————————————————————————————————————
## 五、基因账（§15.8 第 1 项：chrY / chrM 被静默丢弃 —— 必须逐项记账）
## —————————————————————————————————————————————————————————————
our_sym <- rownames(obj)
gm      <- getGenes()                                   # 内置 Ensembl v113，不联网
keep_ch <- c(1:22, "X")
gm_f <- gm[gm$gene_biotype %in% c("protein_coding", "lncRNA") &
           gm$chromosome_name %in% keep_ch & gm$hgnc_symbol != "", ]
tok  <- unique(gm_f$hgnc_symbol)
in_tab   <- our_sym %in% tok
gene_acct <- list(
  n_symbols_matrix = length(our_sym),
  n_matched_fastcnv = sum(in_tab),
  n_unmatched_fastcnv = sum(!in_tab),
  n_unmatched_examples = head(our_sym[!in_tab], 15),
  n_in_fastcnv_chrX = sum(in_tab & our_sym %in% unique(gm_f$hgnc_symbol[gm_f$chromosome_name == "X"])),
  rule = "只留 protein_coding+lncRNA 且 chr∈{1..22,X}；chrY/chrM 静默丢弃（§15.8 第 1 项）")
step("基因账：我方 symbol %d − 未匹配工具表 %d = 进管线候选 %d（chrY/chrM 由工具静默丢弃）",
     gene_acct$n_symbols_matrix, gene_acct$n_unmatched_fastcnv, gene_acct$n_matched_fastcnv)
if (length(gene_acct$n_unmatched_examples))
  step("  未匹配示例：%s", paste(gene_acct$n_unmatched_examples, collapse = ", "))

## —————————————————————————————————————————————————————————————
## 六、跑 fastCNV —— 两条参考口径
## —————————————————————————————————————————————————————————————
## 🔴 走**单对象**路径（`CNVCalling`），不走 `fastCNV()` 的多对象汇池路径（`CNVCallingList`）。
##    原因（§15.8 第 6 项，实测）：汇池版的窗口取均值**漏了长度 1 的保护**
##        funGenomicScore <- function(normcounts, GW) sapply(GW, function(g) colMeans(normcounts[g,]))
##    某些近端着丝粒短臂（如 13p/21p/22p）进 top-7000 的基因可能只有 1 个 ⇒ 切片退化为向量 ⇒
##    `colMeans(x)` 报 'x must be an array of at least two dimensions' 而**硬崩**。
##    单对象版 `CNVCalling` 有 `if (length(g) == 1) normCounts[g, ] else colMeans(...)` 的保护。
##    ⇒ 本体检改用单对象路径；参考／观测的切分**完全相同**（同为上皮 spot 上的同一标签），
##      且**不修改工具任何代码**。汇池路径的不可用本身作为发现登记。
run_one <- function(tag, var, lab) {
  step("---- 运行 %s：referenceVar=%s，referenceLabel=%s（单标签、单对象路径）----", tag, var, lab)
  t0 <- Sys.time()
  o <- CNVCalling(obj, assay = F_ASSAY, referenceVar = var, referenceLabel = lab,
                  scaleOnReferenceLabel = TRUE, thresholdPercentile = F_THRESH_PCT,
                  geneMetadata = gm, windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
                  saveGenomicWindows = FALSE, topNGenes = F_TOPN_GENES)
  if (F_GET_ARM) o <- CNVPerChromosomeArm(o)
  d <- elapsed(t0)
  RUN_MIN[[length(RUN_MIN) + 1L]] <<- setNames(d, tag)
  step("---- %s 完成，用时 %.2f min ----", tag, d)
  o
}
## —————————————————————————————————————————————————————————————
## 七、汇总工具（在每次运行后**立即**取数，好把大对象及时放掉，护住背景进程）
## —————————————————————————————————————————————————————————————
## ⚠️ 循环性声明（§15.3）：作参考的那一组，其 cnv_fraction 由构造必然接近 0
##    ⇒ 「参考组平」不是证据，**不得**写成「方法验证通过」
collect <- function(o, tag) {
  do.call(rbind, lapply(sort(unique(o$slide)), function(s) {
    v <- FetchData(o, vars = "cnv_fraction")[o$slide == s, 1]
    data.frame(run = tag, slide = s, stage = stage_of(s), n = length(v),
               n_na = sum(!is.finite(v)), mean = mean(v, na.rm = TRUE),
               q25 = unname(quantile(v, .25, na.rm = TRUE)), median = median(v, na.rm = TRUE),
               q75 = unname(quantile(v, .75, na.rm = TRUE)), max = max(v, na.rm = TRUE),
               row.names = NULL)
  }))
}
getv <- function(o, tag) data.frame(run = tag, slide = o$slide,
                                    stage = stage_of(o$slide),
                                    cf = FetchData(o, vars = "cnv_fraction")[[1]],
                                    row.names = NULL)

## V2a：池化三张良性作参考（§11.2 C 主锚）
o_a <- run_one("V2a 池化三张良性作参考", "g_a", "benign_epi")
tab_a <- collect(o_a, "V2a_pooled3benign")
va    <- getv(o_a, "V2a")
nwin  <- nrow(GetAssay(o_a, assay = "genomicScores"))
n_nan_win_a <- sum(is.nan(rowMeans(as.matrix(GetAssayData(o_a, assay = "genomicScores", layer = "data")))))
rm(o_a); invisible(gc())
step("V2a 取数完毕，已释放大对象（当前峰值 %.2f GB）", peak_gb())

## V2b：只用 Normal 一张作参考 ⇒ AAH/AAH-1 成为**留出的**良性观测（唯一能证伪的路子）
o_b <- run_one("V2b 只用 Normal 一张作参考", "g_b", "ref_normal")
tab_b <- collect(o_b, "V2b_Normal_only")
vb    <- getv(o_b, "V2b")
rm(o_b); invisible(gc())

tab <- rbind(tab_a, tab_b)
step("==== 逐切片 cnv_fraction（参考组在构造上必然接近 0，勿作证据）====")
for (i in seq_len(nrow(tab))) step("  %-22s %-26s n=%5d 中位 %.4f [%.4f, %.4f] 最大 %.4f NA=%d",
  tab$run[i], tab$slide[i], tab$n[i], tab$median[i], tab$q25[i], tab$q75[i], tab$max[i], tab$n_na[i])

lud  <- va$cf[va$stage == "LUAD"]
ben  <- va$cf[va$stage %in% c("Normal", "AAH", "AAH-1")]
auc_a <- auc_gt(lud, ben)
ludb <- vb$cf[vb$stage == "LUAD"]
held <- vb$cf[vb$stage %in% c("AAH", "AAH-1")]
norm <- vb$cf[vb$stage == "Normal"]
auc_b <- auc_gt(ludb, held)
step("==== 分离度（无阈值 AUC = P(tumor 得分 > 对照得分)）====")
step("  V2a  LUAD %d vs 良性参考 %d（🔴 循环，仅记录）：AUC = %.4f", length(lud), length(ben), auc_a)
step("  V2b  LUAD %d vs **留出**良性 %d（承重，不循环）  ：AUC = %.4f", length(ludb), length(held), auc_b)
step("  V2b  中位：LUAD %.4f ｜ 留出良性 %.4f ｜ 参考Normal %.4f",
     median(ludb, na.rm = TRUE), median(held, na.rm = TRUE), median(norm, na.rm = TRUE))

## V3 不崩不空（单对象 ⇒ 一套窗口，不是逐切片）
v3 <- list(n_windows = nwin, windows_all_nonempty = nwin > 0,
           n_nan_windows = n_nan_win_a, n_na_total = sum(tab$n_na))
step("V3 窗口数 %d；全 NaN 窗口 %d；cnv_fraction NA 合计 %d",
     v3$n_windows, v3$n_nan_windows, v3$n_na_total)

## —————————————————————————————————————————————————————————————
## 八、出图（英文标签：本机无 CJK 字体）
## —————————————————————————————————————————————————————————————
fig <- file.path(OUTD, "fastcnv_smoke_P4_cnv_fraction.png")
STAGES <- c("Normal", "AAH", "AAH-1", "LUAD")
COLS   <- c(Normal = "grey75", AAH = "steelblue", "AAH-1" = "skyblue", LUAD = "indianred")
png(fig, width = 1700, height = 950, res = 130)
op <- par(mfrow = c(1, 2), mar = c(6.5, 4.5, 3.5, 1))
for (k in list(list(d = va, t = "V2a: reference = Normal+AAH+AAH-1 (pooled)"),
               list(d = vb, t = "V2b: reference = Normal only (AAH/AAH-1 held out)"))) {
  lv <- STAGES[STAGES %in% k$d$stage]
  boxplot(k$d$cf ~ factor(k$d$stage, levels = lv), col = COLS[lv],
          las = 2, ylab = "cnv_fraction (fraction of non-zero windows)",
          xlab = "", main = k$t, cex.main = 0.95)
  abline(h = 0, lty = 3, col = "grey40")
  mtext(sprintf("median: %s", paste(sprintf("%s=%.3f", lv,
        tapply(k$d$cf, factor(k$d$stage, levels = lv), median, na.rm = TRUE)), collapse = "  ")),
        side = 3, line = -1.2, cex = 0.7)
}
par(op); invisible(dev.off())
step("图已出：%s", fig)

## —————————————————————————————————————————————————————————————
## 九、manifest
## —————————————————————————————————————————————————————————————
manifest <- list(
  script = "08_spatial_deconv/12_run_fastcnv_smoke.R",
  purpose = "SPATIAL_CNV_PREREG.md §15：fastCNV 单患者体检（工具层面，**不产生恶性判定**）",
  tool = list(name = "fastCNV", version = as.character(packageVersion("fastCNV")),
              license = "GPL-3",
              source = "bioRxiv 2025.10.22.683855（预印本，未经同行评审）"),
  patient = patient, slides = slides, n_spot_epi = n_epi,
  params = list(assay = F_ASSAY, windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
                topNGenes = F_TOPN_GENES, thresholdPercentile = F_THRESH_PCT,
                getCNVPerChromosomeArm = F_GET_ARM,
                path = "CNVCalling + CNVPerChromosomeArm（单对象）；**未**走 fastCNV()/CNVCallingList 汇池路径",
                avoided = list(prepareCounts = F_PREPARE_COUNTS, aggregFactor = F_AGGREG_FACTOR,
                               pooledReference = F_POOLED_REF,
                               because = paste("这几项是 fastCNV() 层参数；本体检走 CNVCalling 单对象路径，",
                                               "既不触发默认聚合（spot<60,000 会聚合，见 §15.4），",
                                               "也避开汇池路径的窗口长度 1 崩溃（§15.8 第 6 项）")),
                getCNVClusters = F_GET_CLUSTERS, doPlot = F_DO_PLOT),
  caliber = list(reference = CAL_REFERENCE, epi = CAL_EPI, rule = CAL_RULE,
                 v2a_ref_label = "benign_epi (Normal+AAH+AAH-1)",
                 v2b_ref_label = "ref_normal (Normal only)"),
  gene_account = gene_acct,
  per_slide = tab,
  separation = list(
    v2a_auc_luad_vs_benign_reference = round(auc_a, 4),
    v2a_note = "🔴 循环：良性即参考，其分数由构造接近 0，此 AUC 不是证据",
    v2b_auc_luad_vs_heldout_benign = round(auc_b, 4),
    v2b_note = "承重数：留出良性（AAH/AAH-1）非参考 ⇒ 此 AUC 不循环",
    v2b_median_luad = median(ludb, na.rm = TRUE),
    v2b_median_heldout_benign = median(held, na.rm = TRUE),
    v2b_median_ref_normal = median(norm, na.rm = TRUE)),
  v3 = v3,
  v1 = list(
    elapsed_min_load = round(min_load, 2),
    elapsed_min_compute = round(sum(unlist(RUN_MIN)), 2),
    elapsed_min_by_run = round(unlist(RUN_MIN), 2),
    peak_gb_final = round(peak_gb(), 3),
    compare = "infercnv 同患者同输入：990.24 min / 46.04 GB 峰值（§15.0）",
    why_fast = paste("快是因为删掉了聚类/BayesNet/HMM 三步（infercnv 那三步占 958 min），",
                     "不是方法更强（§15.2）")),
  figure = fig,
  elapsed_min_total = round(elapsed(t_start), 2),
  peak_gb_running = round(peak_gb(), 3),
  tool_defects_observed = c(
    paste("多对象汇池路径 CNVCallingList::funGenomicScore 缺长度 1 保护 ⇒ 某臂只有 1 个",
          "top 基因时 colMeans 硬崩（'x must be an array of at least two dimensions'）；",
          "本体检实测崩溃，改用单对象 CNVCalling（其有该保护）。见 §15.8 第 6 项"),
    "多标签参考路径（某标签只出现在一张切片时） sapply 退化致参考基线静默错数（§15.8 第 2 项）",
    "genes_by_arm 补齐循环遍历 chr_arm 而非 chr_arm_full ⇒ 凭空多出 'p'/'q' 两个各 200 基因的伪条目（§15.8 第 3 项）",
    "基因账只留 chr1-22+X，chrY/chrM 静默丢弃，与本臂 §14.3 的 chr_exclude=c('chrM') 不一致（§15.8 第 1 项）"),
  not_claimed = c("不构成恶性判定；不得用于 SC0-SC4",
                  "参考组平 = 构造必然，非方法验证",
                  "CNVClassification 默认阈值 ±0.1 是工具默认，本臂未签",
                  "预印本方法，未经同行评审",
                  "未测多患者汇池（§11.2 借锚）；且该路径实测有崩溃缺陷")
)
jsonlite::write_json(manifest, file.path(OUTD, "fastcnv_smoke_P4_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)
step("==== 体检完成：总用时 %.1f min；峰值内存 %.2f GB ====", elapsed(t_start), peak_gb())
step("产物：%s", OUTD)
