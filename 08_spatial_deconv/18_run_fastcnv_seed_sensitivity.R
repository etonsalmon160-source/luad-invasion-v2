#!/usr/bin/env Rscript
# 18_run_fastcnv_seed_sensitivity.R —— 16_（外部锚体检）的**种子敏感性复核**（SPATIAL_CNV_PREREG.md §18.11）
#
# 用户 2026-09-29 签字：
#   范围 = **只做 16_**（17_/SC 那套的锚下采样不在本检查内 ⇒ 见 §18.11 的「未做」登记）
#   判据 = **判据结论逐 seed 一致 ＋ 中位极差 ≤ 0.01**
#
# 要回答的一件事：16_ 的参考 spot 是「每张上皮抽稀至 ≤2,000」**抽**出来的，抽法带 seed。
#   把 seed 换成 1 / 2 / 3 各跑一遍运行 A，看
#     A) F1 / F2 / F3 三个判据的**结论**（过 / 不过）是否逐 seed 相同
#     B) 每张切片 cnv_fraction **中位**在 3 个 seed 间的极差是否 ≤0.01
#   seed=1 **就是 16_ 那一次**（同一条 RNG 流、同一张序）⇒ 顺带当**复现锚点 R0**：
#   本脚本 seed=1 的逐切片中位必须与 16_ 落盘的 runA_spot_scores.rds **逐位相同**。
#   R0 不过 ⇒ 本检查整体作废（说明 18_ 的数据通路与 16_ 不是同一条），不得报结论。
#
# 🔴 **结构性不做**：F4（N1↔N3 批次解耦）**不含任何抽样**（两张外部切片总是全量上皮）
#    ⇒ 对 seed 免疫，本脚本不重跑它，只在 manifest 里登记这一点（不是「测过且稳」）。
# 🔴 本脚本同样**不产生任何恶性判定**，不进 SC0–SC4，不替换 §3 的 infercnv 复现臂。
#
# 跑法（**只挂 fastcnv 这个库**，别和 infercnv 的库混）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 08_spatial_deconv/18_run_fastcnv_seed_sensitivity.R

suppressMessages({
  library(fastCNV)
  library(Seurat)
  library(Matrix)
  library(rhdf5)
})

## —————————————————————————————————————————————————————————————
## 一、参数（与 16_ **逐位相同**；法则 3.1：一个都不改）
## —————————————————————————————————————————————————————————————
CAL_REFERENCE <- "a"        # [已签 §11 第 1 项] 粗版 6 谱系
CAL_EPI       <- "上皮"     # [已签 §11 第 2 项] 单类目
CAL_RULE      <- "argmax"   # [已签 §12.1] 行内最大类 == 上皮，零阈值

F_WINDOW_SIZE    <- 150
F_WINDOW_STEP    <- 10
F_TOPN_GENES     <- 7000
F_THRESH_PCT     <- 0.01
F_GET_ARM        <- TRUE
F_ASSAY          <- "Spatial"

CAP_REF_PER_SLIDE <- 2000L  # [§17.5 执行前定] 每张参考切片上皮 spot 上限
SEEDS             <- c(1L, 2L, 3L)   # [§18.11 签字] 1 = 16_ 原跑 ⇒ 兼作复现锚点 R0

TH_F1_MEDIAN  <- 0.10
TH_F1_AUC     <- 0.75
TH_FLAT_MEDIAN<- 0.05
TH_SPREAD     <- 0.01       # [§18.11 签字] 每张切片中位在 3 seed 间的极差上限

ROOT   <- "/home/eto/luad_v2"
RES    <- file.path(ROOT, "results/08_spatial_deconv")
VISIUM <- file.path(ROOT, "data/visium_spatial")
EXTDIR <- file.path(ROOT, "data/external/GSE248082")
EXT_RES<- file.path(RES, "rctd_external_GSE248082")
OUTD   <- file.path(RES, "spatial_cnv/fastcnv_seed_sensitivity")
PREV16 <- file.path(RES, "spatial_cnv/fastcnv_anchor_flatness/runA_spot_scores.rds")

REF_SLIDES <- c("GSM9226168_P1_AAH", "GSM9226170_P2_AAH", "GSM9226180_P6_AAH",
                "GSM9226186_P9_AAH", "GSM9226191_P11_AAH", "GSM9226209_P20_AAH",
                "GSM9226214_P22_AAH", "GSM9226220_P24_AAH", "GSM9226222_P25_AAH")
P4_SLIDES  <- c("GSM9226174_P4_Normal", "GSM9226175_P4_AAH",
                "GSM9226176_P4_AAH-1", "GSM9226177_P4_LUAD")
EXT_SLIDES <- c("N1", "N3")
ALL_SLIDES <- c(REF_SLIDES, P4_SLIDES, EXT_SLIDES)

## —————————————————————————————————————————————————————————————
## 二、工具函数（与 16_ **逐字同源**）
## —————————————————————————————————————————————————————————————
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()
RUN_MIN <- list()
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

peak_gb <- function() {
  s <- readLines("/proc/self/status")
  hw <- grep("^VmHWM:", s, value = TRUE)
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", hw)) / 1024^2
}
## 🔴 `step()` 是纯日志 ⇒ **绝不能因为格式串的问题把计算打死**（16_ 第 3 跑就死在这）
step <- function(fmt, ...) {
  out <- tryCatch(
    sprintf(paste0("[%s | %6.1f min | 峰值 %6.2f GB] ", fmt),
            format(Sys.time(), "%H:%M:%S"), elapsed(t_start), peak_gb(), ...),
    error = function(e) sprintf("[LOG-ERROR] fmt nchar=%d nargs=%d：%s",
                                nchar(fmt[1]), length(list(...)),
                                conditionMessage(e)))
  cat(paste(out, collapse = "\n"), "\n")
}

read_gz <- function(path, fn) {
  con <- gzfile(path); on.exit(try(close(con), silent = TRUE)); fn(con)
}

auc_gt <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  r <- rank(c(x, y))
  (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y))
}

read_h5_counts <- function(h5) {
  shape <- as.integer(h5read(h5, "matrix/shape"))
  if (length(shape) != 2L) stop("matrix/shape 不是长度 2", call. = FALSE)
  ft_name <- as.character(h5read(h5, "matrix/features/name"))
  ft_type <- as.character(h5read(h5, "matrix/features/feature_type"))
  bc      <- as.character(h5read(h5, "matrix/barcodes"))
  m <- sparseMatrix(i = as.integer(h5read(h5, "matrix/indices")) + 1L,
                    p = as.integer(h5read(h5, "matrix/indptr")),
                    x = as.double(h5read(h5, "matrix/data")),
                    dims = shape)
  if (nrow(m) != length(ft_name) || ncol(m) != length(bc))
    stop("h5 矩阵形状与 features / barcodes 长度不符", call. = FALSE)
  rownames(m) <- ft_name; colnames(m) <- bc
  list(counts = as(m, "CsparseMatrix"), symbols = ft_name, types = ft_type,
       n_ab_rows = sum(ft_type == "Antibody Capture"))
}

agg_by_symbol <- function(m, syms) {
  us <- unique(syms)
  if (length(us) == length(syms) && !anyDuplicated(syms)) {
    rownames(m) <- us; return(m)
  }
  agg <- sparseMatrix(i = match(syms, us), j = seq_along(syms), x = 1,
                      dims = c(length(us), length(syms)))
  m2 <- as(agg %*% m, "CsparseMatrix"); rownames(m2) <- us; m2
}

## 上皮 spot（口径 argmax == 上皮，零阈值，并列剔除；与 11_/12_/16_ 逐字同源）
epi_of <- function(slide, wf) {
  w <- read.delim(gzfile(wf), row.names = 1, check.names = FALSE)
  stopifnot(CAL_EPI %in% colnames(w), ncol(w) == 6L)
  W    <- as.matrix(w)
  mx   <- do.call(pmax, as.data.frame(W))
  n_at <- rowSums(W >= mx)
  am   <- colnames(W)[max.col(W, ties.method = "first")]
  list(epi = rownames(W)[am == CAL_EPI & n_at <= 1L], n_all = nrow(W), n_tie = sum(n_at > 1L))
}

getv <- function(o, tag, grp_col) {
  md <- o@meta.data
  for (k in c("slide", grp_col)) if (!k %in% names(md)) stop("元数据缺列：", k)
  n  <- ncol(o)
  sl <- as.character(md[["slide"]]); gp <- as.character(md[[grp_col]])
  cfv <- FetchData(o, vars = "cnv_fraction")[[1]]
  if (length(sl) != n || length(gp) != n || length(cfv) != n)
    stop(sprintf("getv 长度不齐：ncol=%d slide=%d %s=%d cf=%d",
                 n, length(sl), grp_col, length(gp), length(cfv)))
  if (any(nchar(gp) > 100)) stop("getv：分组列出现超长字符串（疑似被 deparse），拒绝继续")
  data.frame(run = tag, slide = sl, grp = gp, cf = cfv, row.names = NULL)
}
collect <- function(v, tag) do.call(rbind, lapply(sort(unique(v$slide)), function(s) {
  x <- v$cf[v$slide == s]
  data.frame(run = tag, slide = as.character(s),
             group = paste(unique(as.character(v$grp[v$slide == s])), collapse = "|"),
             n = length(x),
             n_na = sum(!is.finite(x)), mean = mean(x, na.rm = TRUE),
             q25 = unname(quantile(x, .25, na.rm = TRUE)), median = median(x, na.rm = TRUE),
             q75 = unname(quantile(x, .75, na.rm = TRUE)), max = max(x, na.rm = TRUE),
             row.names = NULL)
}))

step("==== 种子敏感性复核开始（§18.11）====")
step("工具 fastCNV %s / Seurat %s", as.character(packageVersion("fastCNV")),
     as.character(packageVersion("Seurat")))
step("🔴 本脚本**不产生恶性判定**，不得用于 SC0–SC4")
step("范围=只做 16_；seed=%s；判据 = 结论逐 seed 一致 ＋ 逐切片中位极差 ≤%.2f",
     paste(SEEDS, collapse = ","), TH_SPREAD)

## —————————————————————————————————————————————————————————————
## 三、上皮名单 + 表达矩阵（**只读一次**，三个 seed 共用）
##    —— 与 16_ 的唯一差别：16_ 是「读 → 抽稀 → 建对象」三步串在一起；
##       本脚本把「读」提前做一次、「抽稀」留在 seed 循环里。
##       等价性不靠我说，靠 **R0**（seed=1 必须逐位复现 16_ 落盘的中位）。
## —————————————————————————————————————————————————————————————
t_read <- Sys.time()
ref_epi <- list(); obs_epi <- list()
step("---- 上皮名单（读已签的 RCTD 权重口径；**不重算 RCTD**）----")
for (s in REF_SLIDES) {
  e <- epi_of(s, file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s)))
  stopifnot(length(e$epi) > 0L)
  ref_epi[[s]] <- sort(e$epi)
  step("  %-24s [参考候选] 上皮 %6d", s, length(ref_epi[[s]]))
}
for (s in P4_SLIDES) {
  e <- epi_of(s, file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s)))
  stopifnot(length(e$epi) > 0L)
  obs_epi[[s]] <- sort(e$epi)
  step("  %-24s [观测·同队列] 上皮 %6d（并列剔除 %d）", s, length(obs_epi[[s]]), e$n_tie)
}
for (s in EXT_SLIDES) {
  ef <- file.path(EXT_RES, sprintf("per_slide/%s.a.epi_barcodes.txt.gz", s))
  if (!file.exists(ef))
    stop(sprintf("%s 的上皮名单尚未产出（%s 不存在）⇒ 硬停", s, ef), call. = FALSE)
  obs_epi[[s]] <- sort(read_gz(ef, readLines))
  step("  %-24s [观测·外部] 上皮 %6d", s, length(obs_epi[[s]]))
}

## 逐切片：读表达（GE 行）→ 按 symbol 求和 → **只留该切片的上皮候选列**
##   参考切片留**全部**上皮（抽稀在 seed 循环里做），观测切片留全部上皮（不抽稀）
strip <- list()
for (s in ALL_SLIDES) {
  is_ref <- s %in% REF_SLIDES
  keep   <- if (is_ref) ref_epi[[s]] else obs_epi[[s]]

  if (s %in% EXT_SLIDES) {
    gsm <- if (s == "N1") "GSM8087031" else "GSM8087033"
    h5c <- read_h5_counts(file.path(EXTDIR,
             sprintf("%s_%s_filtered_feature_bc_matrix.h5", gsm, s)))
    if (h5c$n_ab_rows != 0L)
      stop(s, ": 外部切片出现抗体行，与 §16.2 的核验不符 ⇒ 硬停", call. = FALSE)
    ge <- h5c$types == "Gene Expression"
    m  <- agg_by_symbol(as(h5c$counts[ge, , drop = FALSE], "CsparseMatrix"), h5c$symbols[ge])
    rm(h5c); invisible(gc())
  } else {
    d      <- file.path(VISIUM, s, "filtered_feature_bc_matrix")
    ft_all <- read_gz(file.path(d, "features.tsv.gz"), function(con)
      read.delim(con, header = FALSE, col.names = c("id", "sym", "type"),
                 stringsAsFactors = FALSE))
    ab <- which(ft_all$type == "Antibody Capture")
    stopifnot(length(ab) %in% c(0L, 35L))
    bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
    m  <- readMM(file.path(d, "matrix.mtx.gz"))
    stopifnot(nrow(m) == nrow(ft_all), ncol(m) == length(bc))
    ## 🔴 `readMM()` 不设 dimnames ⇒ 必须显式补列名（16_ 第 1 跑就死在这）
    colnames(m) <- bc
    if (length(ab)) { m <- m[-ab, , drop = FALSE]; ft <- ft_all[-ab, , drop = FALSE] } else { ft <- ft_all }
    stopifnot(nrow(ft) == 18085L, all(ft$type == "Gene Expression"))
    m <- agg_by_symbol(as(m, "CsparseMatrix"), ft$sym)
  }

  use <- match(keep, colnames(m), nomatch = 0L)
  stopifnot(!any(use == 0L), !anyDuplicated(use))
  strip[[s]] <- m[, use, drop = FALSE]     # 不改列名：改名放到 seed 循环里（与 16_ 同步）
  rm(m); invisible(gc())
  step("  %-24s %s 读入上皮候选 %6d", s, if (s %in% EXT_SLIDES) "[外部]" else "[队列]",
       ncol(strip[[s]]))
}

## 对齐到同一 symbol 序（**与列子集无关** ⇒ 与 16_ 的 union 相同）
sym_master <- Reduce(union, lapply(strip, rownames))
for (s in names(strip)) {
  m <- strip[[s]]; add <- setdiff(sym_master, rownames(m))
  if (length(add)) {
    z <- Matrix::Matrix(0, nrow = length(add), ncol = ncol(m), sparse = TRUE)
    rownames(z) <- add
    m <- rbind(m, z)
  }
  strip[[s]] <- m[sym_master, , drop = FALSE]
}
step("上皮候选矩阵就绪：%d symbol x %d spot（用时 %.1f min；峰值 %.2f GB）",
     length(sym_master), sum(vapply(strip, ncol, integer(1))), elapsed(t_read), peak_gb())

## —————————————————————————————————————————————————————————————
## 四、seed 循环：抽参考 spot → 建对象 → 运行 A
##    **抽法与 16_ 逐字相同**：`set.seed(seed)` 一次，然后按 `REF_SLIDES` 顺序
##    逐张 `sort(epi)[sort(sample.int(length(epi), CAP))]`（一条 RNG 流，逐张推进）
## —————————————————————————————————————————————————————————————
gm <- getGenes()
t_loop <- Sys.time()
res <- list(); picks_by_seed <- list()

for (sd in SEEDS) {
  step("============ seed = %d ============", sd)
  set.seed(sd)
  picks <- list(); sampled <- list()
  for (s in REF_SLIDES) {
    ep <- ref_epi[[s]]
    pick <- if (length(ep) > CAP_REF_PER_SLIDE)
      ep[sort(sample.int(length(ep), CAP_REF_PER_SLIDE))] else ep
    picks[[s]] <- pick
    sampled[[s]] <- list(n_epi_total = length(ep), n_used = length(pick), cap = CAP_REF_PER_SLIDE)
    step("  参考 %-22s 上皮 %6d ⇒ 取 %6d", s, length(ep), length(pick))
  }
  picks_by_seed[[as.character(sd)]] <- picks

  mats <- list(); cell_slide <- character(0); g_a <- character(0)
  for (s in ALL_SLIDES) {
    is_ref <- s %in% REF_SLIDES
    keep   <- if (is_ref) picks[[s]] else obs_epi[[s]]
    lab    <- if (is_ref) "cohort_benign_epi" else "obs_other"
    m <- strip[[s]][, match(keep, colnames(strip[[s]])), drop = FALSE]
    colnames(m) <- paste(s, keep, sep = "_")
    mats[[s]]  <- m
    cell_slide <- c(cell_slide, setNames(rep(s, ncol(m)), colnames(m)))
    g_a        <- c(g_a, rep(lab, ncol(m)))
    rm(m); invisible(gc())
  }
  n_zero <- sum(vapply(mats, function(m) sum(Matrix::colSums(m) == 0), numeric(1)))
  mat <- do.call(cbind, mats); rm(mats); invisible(gc())
  stopifnot(ncol(mat) == length(g_a), identical(colnames(mat), names(cell_slide)))
  step("  seed=%d 装配：%d symbol x %d spot（全零列 %d）", sd, nrow(mat), ncol(mat), n_zero)

  n_used <- ncol(mat)
  obj <- CreateSeuratObject(counts = mat, assay = F_ASSAY,
                            meta.data = data.frame(g_a = g_a, slide = unname(cell_slide),
                                                   row.names = colnames(mat)))
  rm(mat); invisible(gc())

  t0 <- Sys.time()
  o_a <- CNVCalling(obj, assay = F_ASSAY, referenceVar = "g_a",
                    referenceLabel = "cohort_benign_epi", scaleOnReferenceLabel = TRUE,
                    thresholdPercentile = F_THRESH_PCT, geneMetadata = gm,
                    windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
                    saveGenomicWindows = FALSE, topNGenes = F_TOPN_GENES)
  if (F_GET_ARM) o_a <- CNVPerChromosomeArm(o_a)
  d_seed <- elapsed(t0)
  ## 🔴 顶层**必须**用 `<-`：`RUN_MIN[[i]] <<- …` 在顶层（globalenv 之上没有 frame）
  ##    会报 `object 'RUN_MIN' not found` ⇒ **seed 1 的 3.7 min 白烧过一次**
  RUN_MIN[[length(RUN_MIN) + 1L]] <- setNames(d_seed, sprintf("seed_%d", sd))
  step("  seed=%d 运行 A 完成，用时 %.2f min（峰值 %.2f GB）", sd, d_seed, peak_gb())

  va    <- getv(o_a, sprintf("seed_%d", sd), "g_a")
  tab_a <- collect(va, sprintf("seed_%d", sd))
  nwin  <- nrow(GetAssay(o_a, assay = "genomicScores"))
  ## 🔴 **拿到分数先落盘**，再去做任何判据/记账（16_ 的教训：跑完才崩 ⇒ 全丢）
  saveRDS(list(va = va, tab_a = tab_a, picks = picks),
          file.path(OUTD, sprintf("seed%d_spot_scores.rds", sd)))
  step("  seed=%d 逐 spot 分数已落盘", sd)
  rm(o_a, obj); invisible(gc())

  cf <- function(s) va$cf[va$slide == s]
  ref_cf <- va$cf[va$grp == "cohort_benign_epi"]
  f1_med <- median(cf("GSM9226177_P4_LUAD"), na.rm = TRUE)
  f1_auc <- auc_gt(cf("GSM9226177_P4_LUAD"), ref_cf)
  f2 <- c(P4_Normal = median(cf("GSM9226174_P4_Normal"), na.rm = TRUE),
          P4_AAH    = median(cf("GSM9226175_P4_AAH"),    na.rm = TRUE),
          "P4_AAH-1"= median(cf("GSM9226176_P4_AAH-1"),  na.rm = TRUE))
  f3 <- c(N1 = median(cf("N1"), na.rm = TRUE), N3 = median(cf("N3"), na.rm = TRUE))
  V  <- list(F1 = (f1_med >= TH_F1_MEDIAN) && (f1_auc >= TH_F1_AUC),
             F2 = all(f2 <= TH_FLAT_MEDIAN),
             F3 = all(f3 <= TH_FLAT_MEDIAN))
  res[[as.character(sd)]] <- list(seed = sd, tab = tab_a, med = setNames(tab_a$median, tab_a$slide),
                                  f1_med = f1_med, f1_auc = f1_auc, f2 = f2, f3 = f3,
                                  verdict = V, nwin = nwin, n_used = n_used, elapsed = d_seed,
                                  ref_sampled = sampled)
  step("  seed=%d 判据：F1 %s（中位 %.4f / AUC %.4f） F2 %s（%s） F3 %s（%s）",
       sd, if (isTRUE(V$F1)) "过" else "不过", f1_med, f1_auc,
       if (isTRUE(V$F2)) "过" else "不过", paste(sprintf("%.4f", f2), collapse = ","),
       if (isTRUE(V$F3)) "过" else "不过", paste(sprintf("%.4f", f3), collapse = ","))
}
step("三个 seed 总用时 %.1f min（峰值 %.2f GB）", elapsed(t_loop), peak_gb())

## —————————————————————————————————————————————————————————————
## 五、R0 —— 复现锚点：seed=1 必须逐位复现 16_ 落盘的中位
## —————————————————————————————————————————————————————————————
step("==== R0 复现锚点（seed=1 vs 16_ 落盘 runA_spot_scores.rds）====")
R0 <- NA; r0_maxdiff <- NA_real_; r0_n <- 0L
if (!file.exists(PREV16)) {
  step("  🔴 16_ 落盘文件不存在（%s）⇒ R0 无法判，本检查整体作废", PREV16)
} else {
  prev <- readRDS(PREV16)
  m16  <- tapply(prev$va$cf, prev$va$slide, median, na.rm = TRUE)
  m1   <- res[["1"]]$med
  common <- intersect(names(m16), names(m1))
  r0_maxdiff <- max(abs(m16[common] - m1[common]))
  r0_n <- length(common)
  R0 <- (r0_n >= 15L) && isTRUE(r0_maxdiff <= 1e-9)
  step("  共同切片 %d 张；逐切片中位**最大绝对差** = %.3e（判「逐位相同」需 ≤1e-09）⇒ %s",
       r0_n, r0_maxdiff, if (isTRUE(R0)) "过" else "🔴 不过")
  if (!isTRUE(R0)) {
    bad <- common[abs(m16[common] - m1[common]) > 1e-9]
    step("  🔴 不一致的切片：%s", paste(sprintf("%s(16_=%.6f,18_=%.6f)", bad,
                                               m16[bad], m1[bad]), collapse = " "))
  }
}

## —————————————————————————————————————————————————————————————
## 六、种子敏感性判定（§18.11 已签判据）
## —————————————————————————————————————————————————————————————
lev  <- ALL_SLIDES
M    <- sapply(SEEDS, function(sd) res[[as.character(sd)]]$med[lev])
spread <- apply(M, 1, function(z) max(z) - min(z))
max_spread <- max(spread)
max_spread_slide <- names(spread)[which.max(spread)]

V <- lapply(SEEDS, function(sd) res[[as.character(sd)]]$verdict)
names(V) <- sprintf("seed_%d", SEEDS)
same_verdict <- all(vapply(c("F1", "F2", "F3"), function(k)
  length(unique(vapply(V, function(x) isTRUE(x[[k]]), logical(1)))) == 1L, logical(1)))

A_pass <- same_verdict
B_pass <- isTRUE(max_spread <= TH_SPREAD)
seed_pass <- isTRUE(A_pass && B_pass)

step("==== 判据 A：结论逐 seed 一致 ====")
for (sd in SEEDS) {
  v <- res[[as.character(sd)]]$verdict
  step("  seed=%d  F1=%s  F2=%s  F3=%s", sd,
       if (isTRUE(v$F1)) "过" else "不过", if (isTRUE(v$F2)) "过" else "不过",
       if (isTRUE(v$F3)) "过" else "不过")
}
step("  ⇒ %s", if (A_pass) "过（三条判据的结论逐 seed 相同）" else "🔴 不过（有条目换 seed 就翻）")

step("==== 判据 B：逐切片中位极差 ≤%.2f ====", TH_SPREAD)
for (s in lev)
  step("  %-24s %s  极差 %.4f", s,
       paste(sprintf("s%d=%.4f", SEEDS, M[s, ]), collapse = "  "), spread[[s]])
step("  全切片最大极差 = %.4f（在 %s）⇒ %s", max_spread, max_spread_slide,
     if (B_pass) "过" else "🔴 不过")

## 参考抽样确实换了（否则判据 B 是空转）：逐参考切片 seed1↔seed2 / seed1↔seed3 的 Jaccard
jac <- list()
for (s in REF_SLIDES) {
  a <- picks_by_seed[["1"]][[s]]
  jac[[s]] <- vapply(c("2", "3"), function(k) {
    b <- picks_by_seed[[k]][[s]]
    length(intersect(a, b)) / length(union(a, b))
  }, numeric(1))
}
step("---- 参考抽样的换手率（Jaccard，1=完全没换）----")
for (s in REF_SLIDES)
  step("  %-24s vs seed2 %.4f  vs seed3 %.4f", s, jac[[s]][1], jac[[s]][2])

step("==== 总判：%s（A %s ／ B %s；R0 %s）====",
     if (seed_pass && isTRUE(R0)) "过" else if (!isTRUE(R0)) "🔴 作废（R0 不过）" else "🔴 不过",
     if (A_pass) "过" else "不过", if (B_pass) "过" else "不过",
     if (isTRUE(R0)) "过" else "不过")

## —————————————————————————————————————————————————————————————
## 七、出图（英文标签：本机无 CJK 字体）
## —————————————————————————————————————————————————————————————
fig <- file.path(OUTD, "fastcnv_seed_sensitivity.png")
png(fig, width = 1900, height = 880, res = 130)
op <- par(mfrow = c(1, 2), mar = c(9, 4.5, 3.5, 1))
cols <- c("black", "steelblue", "darkorange")
draw <- function(ix, ylim, main) {
  plot(NA, xlim = c(0.5, length(ix) + 0.5), ylim = ylim, xaxt = "n",
       xlab = "", ylab = "median cnv_fraction", main = main, cex.main = 0.95)
  axis(1, at = seq_along(ix), labels = lev[ix], las = 2, cex.axis = 0.62)
  abline(h = TH_FLAT_MEDIAN, lty = 2, col = "grey30")
  abline(h = TH_F1_MEDIAN,   lty = 3, col = "grey30")
  for (j in seq_along(SEEDS))
    points(seq_along(ix) + (j - 2) * 0.12, M[ix, j], pch = 19, cex = 0.75, col = cols[j])
}
draw(seq_along(lev), c(0, max(M) * 1.08), "Run A medians by seed (all slides)")
legend("topleft", c(sprintf("seed %d", SEEDS), sprintf("%.2f flat threshold", TH_FLAT_MEDIAN),
                    sprintf("%.2f F1 threshold", TH_F1_MEDIAN)),
       pch = c(rep(19, 3), NA, NA), lty = c(NA, NA, NA, 2, 3),
       col = c(cols, "grey30", "grey30"), cex = 0.68, bty = "n")
ix_b <- which(!lev %in% "GSM9226177_P4_LUAD")
draw(ix_b, c(0, 0.13), "Same, excluding positive control (zoom)")
par(op); invisible(dev.off())
step("图已出：%s", fig)

## —————————————————————————————————————————————————————————————
## 八、manifest
## —————————————————————————————————————————————————————————————
manifest <- list(
  script = "08_spatial_deconv/18_run_fastcnv_seed_sensitivity.R",
  purpose = "SPATIAL_CNV_PREREG.md §18.11：16_ 外部锚体检的种子敏感性复核",
  signed_caliber = list(scope = "只做 16_（外部锚体检）",
                        criterion = "判据结论逐 seed 一致 ＋ 逐切片中位极差 ≤0.01",
                        signed_by_user_on = "2026-09-29"),
  not_producing = c("不产生恶性判定", "不进 SC0-SC4", "不替换 §3 infercnv 复现臂"),
  not_tested = list(F4 = "N1<->N3 批次解耦**不含抽样** ⇒ 对 seed 免疫，本检查未重跑；不是「测过且稳」"),
  seeds = SEEDS, cap_ref_per_slide = CAP_REF_PER_SLIDE,
  thresholds = list(f1_median = TH_F1_MEDIAN, f1_auc = TH_F1_AUC,
                    flat_median = TH_FLAT_MEDIAN, spread_max = TH_SPREAD,
                    provenance = "前三条与 16_ 逐位相同；0.01 由用户 2026-09-29 签字"),
  R0_reproduction = list(compared_to = PREV16, n_slides = r0_n,
                         max_abs_diff_median = r0_maxdiff, pass = R0,
                         rule = "seed=1 必须逐位复现 16_（≤1e-09）；不过 ⇒ 本检查作废"),
  per_seed = lapply(SEEDS, function(sd) {
    r <- res[[as.character(sd)]]
    list(seed = sd, n_used_spot = r$n_used, n_windows = r$nwin,
         elapsed_min = round(r$elapsed, 2), f1_median_luad = r$f1_med,
         f1_auc_vs_ref = r$f1_auc, f2_medians = as.list(r$f2), f3_medians = as.list(r$f3),
         verdicts = r$verdict, ref_sampled = r$ref_sampled)
  }),
  median_by_slide = as.list(as.data.frame(M)),
  spread_by_slide = as.list(spread),
  max_spread = max_spread, max_spread_slide = max_spread_slide,
  ref_sampling_jaccard_vs_seed1 = jac,
  criteria = list(A_verdicts_identical = A_pass, B_spread_le_threshold = B_pass,
                  R0_reproduction = R0, overall_pass = seed_pass && isTRUE(R0)),
  reading_limits = c(
    "本检查只证「结论对这个参考子样本不敏感」，**不**证「结论对一般情形稳」",
    "16_ 的 F2/F3 离阈值很远（0.0000-0.0093 vs 上限 0.05），故 A 判据更接近「没翻」而非「压线稳住」",
    "16_ 的 F2/F3 结论不变 ⇒ §18.10 的结论（外部两张都不得当主锚）**不因换 seed 而翻**",
    "F4 未重跑（结构上无抽样）"),
  figure = fig,
  elapsed_min_total = round(elapsed(t_start), 2)
)
jsonlite::write_json(manifest, file.path(OUTD, "fastcnv_seed_sensitivity_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)
step("==== 种子检查完成：总用时 %.1f min；峰值 %.2f GB ====", elapsed(t_start), peak_gb())
step("产物：%s", OUTD)
