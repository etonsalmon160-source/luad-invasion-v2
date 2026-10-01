#!/usr/bin/env Rscript
# 19_run_fastcnv_cohort.R —— 借锚臂**全队列**（SPATIAL_CNV_PREREG.md §19，2026-09-29）
#
# 逐患者跑，一次调用一个患者（便于崩了只重跑该患者、且内存逐患者释放）。
#   --patient P7            （必需）
#   --arms    all|main|holdout|d|borrow_ctrl   （默认 all）
#
# 本臂的门与参数**一字不改**：
#   引擎 §18.1 fastCNV v1.1.11；参数 §15.4；阈值 §18.5(a) 0.05 ＋ 四档；门 §18.9 SC0–SC4。
#   路径 = **单对象 `CNVCalling` ＋ `CNVPerChromosomeArm`**；🔴 **不用 `CNVCallingList`**（§19.4）。
#
# 🔴 本脚本**不产生恶性判定**：`cnv_fraction` 是连续量；`0.05` 只用来**报比例**，不用来筛 spot、
#    不设通过率阈值（§19.7 SC4）。`CNVClassification` 只是方向分类器。
#
# 🔴 参考组织是**良性/AAH 上皮**，不是组织学正常上皮 ⇒ 逐次写明（§19.6 第 1 条）。
#    跨性别借锚**禁止**（chrX 剂量，§19.6 第 2 条）。
#
# 跑法（**只挂 fastcnv 这个库**）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 08_spatial_deconv/19_run_fastcnv_cohort.R --patient P7

suppressMessages({ library(fastCNV); library(Seurat); library(Matrix) })

## —————————————————————————————————————————————————————————————
## 一、参数（逐项显式；法则 3.1）
## —————————————————————————————————————————————————————————————
CAL_EPI       <- "上皮"      # [已签 §12.1] 单类目
CAL_RULE      <- "argmax"    # [已签 §12.1] 行内最大类 == 上皮，零阈值，**并列剔除**
F_WINDOW_SIZE <- 150;  F_WINDOW_STEP <- 10
F_TOPN_GENES  <- 7000; F_THRESH_PCT  <- 0.01
F_GET_ARM     <- TRUE;  F_ASSAY <- "Spatial"
N_REPEAT      <- 3L          # §18.9 SC2
SC2_MIN_IDENT <- 0.95        # §18.9 SC2
HOLD_FRAC     <- 0.20        # [§19.5 🕐4] SC3 留出比例
SEED          <- 1L          # 固定种子（SC3 留出用）
TH_MAIN       <- 0.05        # [已签 §18.5(a)] 只用于**报比例**
TIERS         <- c(0.02, 0.05, 0.10, 0.15)   # [已签 §18.9 SC4] 四档敏感性

## SC0 期望值（§18.9①②；逐张实测口径）
EXP_GE_ROWS    <- 18085L
EXP_AB_SET     <- c(0L, 35L) # 🔴 抗体行**并非每张都有**：45/56 有 35 行、11/56 为 0
SC0_MIN_MATCH  <- 0.99

## 🕐 内存守卫（§19.11 建议项；本机有背景病理进程）
MIN_FREE_GB <- 80            # 可用内存低于此值 ⇒ **停并上报**，不许被 OOM 静默杀掉

ROOT   <- "/home/eto/luad_v2"
RES    <- file.path(ROOT, "results/08_spatial_deconv")
VISIUM <- file.path(ROOT, "data/visium_spatial")
GODIR  <- file.path(RES, "spatial_cnv")
OUTD   <- file.path(RES, "spatial_cnv/cohort_borrow")
EXTDIR <- file.path(ROOT, "data/external/GSE248082")   # 本臂**不引入**外部锚（§19.6 第 5 条）

## 供体/受体（§19.3，2026-09-29 从 GEO 系列元数据重取核对）
DONORS_F <- c("P1", "P4", "P9", "P11", "P22")
DONORS_M <- c("P2", "P6", "P20", "P24", "P25")
SEX <- c(P1="F", P4="F", P9="F", P11="F", P22="F",
         P2="M", P6="M", P20="M", P24="M", P25="M",
         P7="F", P8="F", P10="F", P12="F", P16="F", P18="F", P21="F",
         P3="M", P5="M", P13="M", P14="M", P15="M", P17="M", P19="M", P23="M")

## —————————————————————————————————————————————————————————————
## 二、工具函数
## —————————————————————————————————————————————————————————————
args    <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d = NULL) { i <- which(args == k); if (length(i)) args[i + 1] else d }
patient <- get_arg("--patient")
arms_wanted <- strsplit(get_arg("--arms", "all"), ",")[[1]]
if (is.null(patient) || !patient %in% names(SEX)) stop("须给 --patient，且须在 SEX 表内", call. = FALSE)

dir.create(file.path(OUTD, patient), showWarnings = FALSE, recursive = TRUE)
PDIR <- file.path(OUTD, patient)
t_start <- Sys.time()
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
peak_gb <- function() {
  s <- readLines("/proc/self/status")
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", grep("^VmHWM:", s, value = TRUE))) / 1024^2
}
free_gb <- function() {
  s <- readLines("/proc/meminfo")
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", grep("^MemAvailable:", s, value = TRUE))) / 1024^2
}
## `step()` 是纯日志 ⇒ 任何格式化异常降级成一行纯文本（888,629 字那次的教训）
step <- function(fmt, ...) {
  out <- tryCatch(
    sprintf(paste0("[%s | %6.1f min | 峰值 %6.2f GB] ", fmt),
            format(Sys.time(), "%H:%M:%S"), elapsed(t_start), peak_gb(), ...),
    error = function(e) sprintf("[LOG-ERROR] fmt nchar=%d nargs=%d：%s",
                                nchar(fmt[1]), length(list(...)), conditionMessage(e)))
  cat(paste(out, collapse = "\n"), "\n"); flush(stdout())
}
read_gz  <- function(p, fn) { con <- gzfile(p); on.exit(try(close(con), silent = TRUE)); fn(con) }
stage_of <- function(s) sub("^[^_]*_[^_]*_", "", s)
is_benign <- function(s) grepl("^Normal|^AAH", stage_of(s))
patient_of <- function(s) sub("^[^_]*_([^_]*)_?.*$", "\\1", s)
sha256 <- function(p) {
  out <- system2("sha256sum", shQuote(p), stdout = TRUE)
  if (!length(out)) stop("sha256sum 失败：", p)
  sub(" .*$", "", out[1])
}

FAILS <- character(0)
gate <- function(ok, label, detail) {
  step("  [%s] %s —— %s", if (isTRUE(ok)) "PASS" else "FAIL", label, detail)
  if (!isTRUE(ok)) FAILS <- c(FAILS, label)
  isTRUE(ok)
}

step("==== 借锚臂全队列：患者 %s（性别 %s）开始 ====", patient, SEX[[patient]])
step("引擎 fastCNV %s / Seurat %s", as.character(packageVersion("fastCNV")),
     as.character(packageVersion("Seurat")))
step("🔴 参考组织 = 良性/AAH 上皮（**不是**组织学正常上皮，§19.6 第 1 条）；本脚本不产生恶性判定")
step("可用内存 %.1f GB（守卫线 %d GB）", free_gb(), MIN_FREE_GB)
if (free_gb() < MIN_FREE_GB)
  stop(sprintf("可用内存 %.1f GB < 守卫线 %d GB ⇒ 硬停，不抢背景病理进程的内存", free_gb(), MIN_FREE_GB),
       call. = FALSE)

## —————————————————————————————————————————————————————————————
## 三、切片、掩码、上皮/非上皮名单（口径与 17_ 逐字同源）
## —————————————————————————————————————————————————————————————
sm_all <- read.delim(gzfile(file.path(RES, "spot_mask.tsv.gz")), stringsAsFactors = FALSE)
sm_all <- sm_all[as.character(sm_all$pass) %in% c("TRUE", "true", "1"), ]
stopifnot(nrow(sm_all) > 0)

epi_of <- function(s, keep) {
  wf <- file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s))
  stopifnot(file.exists(wf))
  w <- read.delim(gzfile(wf), row.names = 1, check.names = FALSE)
  stopifnot(CAL_EPI %in% colnames(w), ncol(w) == 6L)
  k <- intersect(rownames(w), keep)
  W <- as.matrix(w[k, , drop = FALSE])
  mx <- do.call(pmax, as.data.frame(W)); n_at <- rowSums(W >= mx)
  am <- colnames(W)[max.col(W, ties.method = "first")]
  list(epi = k[am == CAL_EPI & n_at <= 1L], non_epi = k[am != CAL_EPI & n_at <= 1L],
       n_all = length(k), n_tie = sum(n_at > 1L))
}

## 本患者的切片（掩码 pass）
own_slides <- sort(unique(sm_all$slide[grepl(paste0("_", patient, "_"), sm_all$slide)]))
stopifnot(length(own_slides) > 0L)
own_benign <- own_slides[vapply(own_slides, is_benign, logical(1))]
own_lesion <- setdiff(own_slides, own_benign)

epi <- list(); nonepi <- list(); n_pass <- list(); n_tie <- list()
for (s in own_slides) {
  bc <- sm_all$barcode[sm_all$slide == s]
  e  <- epi_of(s, bc)
  epi[[s]] <- e$epi; nonepi[[s]] <- e$non_epi; n_pass[[s]] <- length(bc); n_tie[[s]] <- e$n_tie
}
step("本患者 %d 张切片（良性 %d / 病灶 %d）", length(own_slides), length(own_benign), length(own_lesion))
for (s in own_slides)
  step("  %-26s pass %6d  上皮 %6d  非上皮 %6d  并列剔 %5d  %s",
       s, n_pass[[s]], length(epi[[s]]), length(nonepi[[s]]), n_tie[[s]],
       if (is_benign(s)) "良性" else "病灶")

## 主臂归属：自带锚 10 例 = 自家锚；缺锚 15 例 = 借锚（§19.5 已签）
IS_DONOR_PAT <- patient %in% c(DONORS_F, DONORS_M)
MAIN_KIND <- if (length(own_benign) > 0L && IS_DONOR_PAT) "self" else "borrow"
if (MAIN_KIND == "borrow" && length(own_benign) > 0L)
  stop("该患者有良性切片却不在供体名单 ⇒ 归属不清，硬停", call. = FALSE)

## 供体池（同性别、全取、**leave-self-out**；§19.3 规则）
pool_pats <- setdiff(if (SEX[[patient]] == "F") DONORS_F else DONORS_M, patient)
pool_slides <- character(0); pool_owner <- character(0)
for (dp in pool_pats) {
  ds <- sort(unique(sm_all$slide[grepl(paste0("_", dp, "_"), sm_all$slide)]))
  ds <- ds[vapply(ds, is_benign, logical(1))]
  pool_slides <- c(pool_slides, ds); pool_owner <- c(pool_owner, rep(dp, length(ds)))
}
names(pool_owner) <- pool_slides
for (s in pool_slides) {
  e <- epi_of(s, sm_all$barcode[sm_all$slide == s])
  epi[[s]] <- e$epi; nonepi[[s]] <- e$non_epi
  n_pass[[s]] <- sum(sm_all$slide == s); n_tie[[s]] <- e$n_tie
}
step("主臂 = %s；同性别供体池 = %s（%d 张切片）", MAIN_KIND,
     paste(pool_pats, collapse = ","), length(pool_slides))
for (s in pool_slides)
  step("  [供体 %s] %-26s 上皮 %6d", pool_owner[[s]], s, length(epi[[s]]))

## 参考构成（§19.4：逐次报每位供体贡献的 spot 数）
if (MAIN_KIND == "self") {
  ref_src_slides <- own_benign; ref_kind_lab <- "self"
} else {
  ref_src_slides <- pool_slides; ref_kind_lab <- paste0("pool_", SEX[[patient]])
}
obs_slides <- own_lesion
stopifnot(length(obs_slides) > 0L, length(ref_src_slides) > 0L)

## —————————————————————————————————————————————————————————————
## 四、需要读的切片（并集）与逐张缓存
## —————————————————————————————————————————————————————————————
## 🔴 读片集必须是**所有臂**的并集：自带锚患者的「借锚对照」要用供体池的片子，
##    只在 ref_src_slides（=自家良性）里找会直接崩（P1/P2 实测栽在这里）
ctrl_slides <- if (MAIN_KIND == "self" &&
                   any(arms_wanted %in% c("all", "borrow_ctrl"))) pool_slides else character(0)
need_slides <- unique(c(ref_src_slides, own_slides, ctrl_slides))
step("需读切片 %d 张（主臂参考 %d + 自身 %d + 借锚对照 %d，去重后）", length(need_slides),
     length(ref_src_slides), length(own_slides), length(ctrl_slides))

sc0_ge <- list(); sc0_ab <- list(); sc0_other <- list()
cache <- list()
t_read <- Sys.time()
for (s in need_slides) {
  d      <- file.path(VISIUM, s, "filtered_feature_bc_matrix")
  ft_all <- read_gz(file.path(d, "features.tsv.gz"), function(con)
    read.delim(con, header = FALSE, col.names = c("id", "sym", "type"), stringsAsFactors = FALSE))
  ab <- which(ft_all$type == "Antibody Capture")
  bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
  m  <- readMM(file.path(d, "matrix.mtx.gz"))
  stopifnot(nrow(m) == nrow(ft_all), ncol(m) == length(bc))
  colnames(m) <- bc                                   # readMM 不设 dimnames
  ## SC0②：数 **Gene Expression 行**，不是总行数
  sc0_ge[[s]] <- sum(ft_all$type == "Gene Expression")
  sc0_ab[[s]] <- length(ab)
  sc0_other[[s]] <- setdiff(unique(ft_all$type), c("Gene Expression", "Antibody Capture"))
  if (length(ab)) { m <- m[-ab, , drop = FALSE]; ft <- ft_all[-ab, , drop = FALSE] } else { ft <- ft_all }

  keep <- unique(c(epi[[s]], nonepi[[s]]))            # 该张要用的 mask-pass spot
  use  <- match(keep, bc, nomatch = 0L)
  if (any(use == 0L)) stop(s, "：mask-pass 的 barcode 不在表达矩阵里", call. = FALSE)
  m <- m[, use, drop = FALSE]
  g <- factor(ft$sym, levels = unique(ft$sym))
  Agg <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1,
                      dims = c(nlevels(g), length(g)))
  m <- Agg %*% m; rownames(m) <- levels(g)
  colnames(m) <- paste(s, keep, sep = "_")
  cache[[s]] <- as(m, "CsparseMatrix")
  rm(m); invisible(gc())
  step("  读入 %-26s %d 基因 x %6d spot（%.1f min）", s, nrow(cache[[s]]), ncol(cache[[s]]),
       elapsed(t_read))
}
step("缓存就绪，用时 %.1f min（峰值 %.2f GB）", elapsed(t_read), peak_gb())

## —————————————————————————————————————————————————————————————
## 五、建对象（参考 + 观测；标签 reference / observation）
## —————————————————————————————————————————————————————————————
build_obj <- function(ref_slides, ref_bc, obs_slides_, obs_bc, tag) {
  ms <- list(); grp <- character(0); bcs <- character(0); slv <- character(0)
  for (s in ref_slides) {
    k <- ref_bc[[s]]; if (!length(k)) next
    u <- match(paste(s, k, sep = "_"), colnames(cache[[s]]))
    if (anyNA(u)) stop(tag, "：参考 barcode 未在缓存中（", s, "）", call. = FALSE)
    ms[[length(ms) + 1L]] <- cache[[s]][, u, drop = FALSE]
    n <- length(u)
    grp <- c(grp, rep("reference", n))
    bcs <- c(bcs, paste(s, k, sep = "_")); slv <- c(slv, rep(s, n))
  }
  for (s in obs_slides_) {
    k <- obs_bc[[s]]; if (!length(k)) next
    u <- match(paste(s, k, sep = "_"), colnames(cache[[s]]))
    if (anyNA(u)) stop(tag, "：观测 barcode 未在缓存中（", s, "）", call. = FALSE)
    ms[[length(ms) + 1L]] <- cache[[s]][, u, drop = FALSE]
    n <- length(u)
    grp <- c(grp, rep("observation", n))
    bcs <- c(bcs, paste(s, k, sep = "_")); slv <- c(slv, rep(s, n))
  }
  mat <- do.call(cbind, ms)
  stopifnot(ncol(mat) == length(grp), length(bcs) == length(grp), length(slv) == length(grp))
  o <- CreateSeuratObject(counts = mat, assay = F_ASSAY,
                          meta.data = data.frame(grp = grp, slide = slv,
                                                 row.names = bcs))
  list(obj = o, n_ref = sum(grp == "reference"), n_obs = sum(grp == "observation"))
}

## 一次 CNVCalling（＋逐臂）⇒ 逐 spot cf ＋ 逐臂均值 ＋ 窗口数
run_pass <- function(o, tag) {
  t0 <- Sys.time()
  r <- CNVCalling(o, assay = F_ASSAY, referenceVar = "grp", referenceLabel = "reference",
                  scaleOnReferenceLabel = TRUE, thresholdPercentile = F_THRESH_PCT,
                  geneMetadata = getGenes(), windowSize = F_WINDOW_SIZE,
                  windowStep = F_WINDOW_STEP, saveGenomicWindows = FALSE,
                  topNGenes = F_TOPN_GENES)
  if (F_GET_ARM) r <- CNVPerChromosomeArm(r)
  md <- r@meta.data
  ## 🔴 barcode 守卫：落盘的逐 spot 分数必须能按 barcode 对回 RCTD margin（§13.7.3 的义务）
  stopifnot(!is.null(rownames(md)), length(rownames(md)) == ncol(r),
            !anyDuplicated(rownames(md)))
  cf <- md[["cnv_fraction"]]; nwin <- nrow(GetAssay(r, assay = "genomicScores"))
  ## 46 个染色体臂均值列（`<arm>_CNV`），逐 slice 逐组取均值
  arm_cols <- grep("^[0-9X]+\\.(p|q)_CNV$", names(md), value = TRUE)
  sl <- as.character(md[["slide"]]); gp <- as.character(md[["grp"]])
  per <- do.call(rbind, lapply(sort(unique(sl)), function(s) {
    i <- which(sl == s)
    data.frame(pass = tag, slide = s, group = paste(unique(gp[i]), collapse = "|"),
               n = length(i), n_na = sum(!is.finite(cf[i])),
               mean = mean(cf[i], na.rm = TRUE),
               q25 = unname(quantile(cf[i], .25, na.rm = TRUE)),
               median = median(cf[i], na.rm = TRUE),
               q75 = unname(quantile(cf[i], .75, na.rm = TRUE)),
               max = max(cf[i], na.rm = TRUE), row.names = NULL)
  }))
  arms <- if (length(arm_cols)) {
    do.call(rbind, lapply(sort(unique(sl)), function(s) {
      i <- which(sl == s)
      d <- data.frame(pass = tag, slide = s, group = paste(unique(gp[i]), collapse = "|"),
                      row.names = NULL)
      for (a in arm_cols) d[[a]] <- mean(md[[a]][i], na.rm = TRUE)
      d
    }))
  } else NULL
  d <- elapsed(t0)
  step("  [%s] 完成 %.2f min；窗口 %d；参考组中位 %.4f；观测组中位 %.4f",
       tag, d, nwin, median(cf[gp == "reference"], na.rm = TRUE),
       median(cf[gp == "observation"], na.rm = TRUE))
  ## 🔴 **不把对象带出去**（genomicScores 那个 assay 很大；带出去会把峰值钉住）
  list(cf = cf, slide = sl, grp = gp, per = per, arms = arms, nwin = nwin,
       bc = rownames(md),   # 🔴 barcode 必须落盘（§18 补过的口子；否则 §13.7.3 margin 梯度离线算不了）
       min = d, ref_names = sort(colnames(r)[gp == "reference"]))
}

## 逐位相同比例（完全相等才算；并列不算过）
frac_identical <- function(cl) {
  same <- rep(TRUE, length(cl[[1]]))
  for (k in 2:length(cl)) same <- same & (cl[[1]] == cl[[k]])
  list(frac = mean(same), same = same)
}

## 四档比例（**只报**，不筛）
tier_frac <- function(cf) setNames(vapply(TIERS, function(t) mean(cf >= t, na.rm = TRUE), numeric(1)),
                                   sprintf("ge_%.2f", TIERS))

RES_ALL <- list(); MANI <- list(patient = patient, sex = SEX[[patient]],
                                main_kind = MAIN_KIND, arms_wanted = arms_wanted)

## —————————————————————————————————————————————————————————————
## 六、主臂（SC2 重复 N≥3 次）
## —————————————————————————————————————————————————————————————
run_main <- function(kind, ref_slides, tag_prefix) {
  ref_bc <- setNames(lapply(ref_slides, function(s) sort(epi[[s]])), ref_slides)
  obs_bc <- setNames(lapply(obs_slides, function(s) sort(epi[[s]])), obs_slides)
  B <- build_obj(ref_slides, ref_bc, obs_slides, obs_bc, tag_prefix)
  step("---- %s 对象：%d 基因 x %d spot（参考 %d ＋ 观测 %d）----",
       tag_prefix, nrow(B$obj), ncol(B$obj), B$n_ref, B$n_obs)
  ## SC1：未聚合证据（逐切片 spot 数 == 我们传入的上皮数）
  exp_n <- c(vapply(ref_slides, function(s) length(ref_bc[[s]]), integer(1)),
             vapply(obs_slides, function(s) length(obs_bc[[s]]), integer(1)))
  got_n <- as.integer(table(factor(B$obj$slide, levels = names(exp_n))))
  agg_ok <- identical(unname(got_n), unname(exp_n))
  gate(agg_ok, sprintf("%s SC1 未聚合（逐 spot 分辨率保住）", tag_prefix),
       sprintf("逐切片 %s", paste(sprintf("%s=%d", names(exp_n), exp_n), collapse = " ")))
  if (!agg_ok) step("⚠️ SC1 期望 %s | 实际 %s",
                    paste(sprintf("%s=%d", names(exp_n), exp_n), collapse = " "),
                    paste(sprintf("%s=%s", names(exp_n), ifelse(is.na(got_n), "缺失", got_n)), collapse = " "))

  ## 参考构成（§19.4 逐次报）
  contrib <- if (kind == "self") {
    data.frame(donor = patient, slide = ref_slides,
               n_epi = vapply(ref_slides, function(s) length(ref_bc[[s]]), integer(1)),
               row.names = NULL)
  } else {
    d <- data.frame(donor = pool_owner[ref_slides], slide = ref_slides,
                    n_epi = vapply(ref_slides, function(s) length(ref_bc[[s]]), integer(1)),
                    row.names = NULL)
    d
  }
  for (i in seq_len(nrow(contrib)))
    step("    参考贡献 %-6s %-26s %6d spot", contrib$donor[i], contrib$slide[i], contrib$n_epi[i])

  cl <- list(); anchor_same <- TRUE; run_min <- numeric(0)
  ref_bc_sorted <- sort(colnames(B$obj)[B$obj$grp == "reference"])
  for (k in seq_len(N_REPEAT)) {
    r <- run_pass(B$obj, sprintf("%s_SC2_%d", tag_prefix, k))
    cl[[k]] <- r$cf
    if (!identical(r$ref_names, ref_bc_sorted)) anchor_same <- FALSE
    run_min <- c(run_min, r$min)
    if (k == 1L) first <- r
    rm(r); invisible(gc())
  }
  fi <- frac_identical(cl)
  gate(anchor_same, sprintf("%s SC2 锚名单逐字节相同", tag_prefix),
       sprintf("参考 %d spot，%d 次一致", length(ref_bc_sorted), N_REPEAT))
  gate(fi$frac >= SC2_MIN_IDENT, sprintf("%s SC2 逐 spot 逐位相同", tag_prefix),
       sprintf("%.6f（门槛 %.2f）", fi$frac, SC2_MIN_IDENT))

  list(first = first, fi = fi, contrib = contrib, run_min = run_min,
       ref_slides = ref_slides, ref_bc = ref_bc, obs_bc = obs_bc,
       n_ref = length(ref_bc_sorted), agg_ok = agg_ok)
}

M <- NULL; BR <- NULL
if (any(arms_wanted %in% c("all", "main", "holdout", "d"))) {
  M <- run_main(MAIN_KIND, ref_src_slides, sprintf("main_%s", MAIN_KIND))
}
if (any(arms_wanted %in% c("all", "borrow_ctrl")) && MAIN_KIND == "self") {
  step("---- 对照：自带锚 10 例**另跑一次借锚**（leave-self-out，§19.5 🕐2）----")
  BR <- run_main("borrow", pool_slides, "borrow_ctrl")
}

## —————————————————————————————————————————————————————————————
## 七、SC2 首跑产物落盘（先养活，后面再崩也不重算）
## —————————————————————————————————————————————————————————————
save_pass <- function(r, tag, extra = NULL) {
  f <- file.path(PDIR, sprintf("%s_spot_scores.rds", tag))
  saveRDS(c(list(cf = r$cf, slide = r$slide, grp = r$grp, bc = r$bc, per = r$per,
                 nwin = r$nwin, run_min = r$min), extra), f)
  f
}
if (!is.null(M)) save_pass(M$first, sprintf("main_%s", MAIN_KIND),
                           list(contrib = M$contrib, ref_slides = M$ref_slides))
if (!is.null(BR)) save_pass(BR$first, "borrow_ctrl", list(contrib = BR$contrib))
step("主臂逐 spot 分数已落盘（峰值 %.2f GB）", peak_gb())

## —————————————————————————————————————————————————————————————
## 八、SC3 留出复核（§19.5 🕐4：每张参考切片留出 20%，判据差 ≤0.05）
## —————————————————————————————————————————————————————————————
if (any(arms_wanted %in% c("all", "holdout")) && !is.null(M)) {
  step("---- SC3 留出复核：每张参考切片留出 %.0f%%（seed=%d）----", HOLD_FRAC * 100, SEED)
  set.seed(SEED)
  ref_bc_sub <- setNames(lapply(M$ref_slides, function(s) {
    v <- sort(epi[[s]])
    if (length(v) < 5L) return(v)
    sort(v[sort(sample.int(length(v), max(1L, floor((1 - HOLD_FRAC) * length(v)))))] )
  }), M$ref_slides)
  held <- sum(vapply(M$ref_slides, function(s) length(epi[[s]]), integer(1))) -
          sum(vapply(M$ref_slides, function(s) length(ref_bc_sub[[s]]), integer(1)))
  step("  全量参考 %d ⇒ 留出后 %d spot（留出 %d）",
       M$n_ref, sum(vapply(ref_bc_sub, length, integer(1))), held)
  Ho <- build_obj(M$ref_slides, ref_bc_sub, obs_slides, M$obs_bc, "holdout")
  h  <- run_pass(Ho$obj, "SC3_holdout")
  ## 判据：观测组在 0.05 上的比例，两边差 ≤0.05；且参考组自身中位 ≤0.05
  f_full <- mean(M$first$cf[M$first$grp == "observation"] >= TH_MAIN, na.rm = TRUE)
  f_hold <- mean(h$cf[h$grp == "observation"] >= TH_MAIN, na.rm = TRUE)
  d_frac <- abs(f_full - f_hold)
  ref_med_hold <- median(h$cf[h$grp == "reference"], na.rm = TRUE)
  gate(ref_med_hold <= TH_MAIN, "SC3① 锚内中位 ≤0.05",
       sprintf("%.4f（循环性警告：参考组由构造≈0，不是证据，§19.6 第 4 条）", ref_med_hold))
  gate(d_frac <= 0.05, "SC3② 留出前后观测组比例差 ≤0.05",
       sprintf("全量 %.4f vs 留出 %.4f，差 %.4f", f_full, f_hold, d_frac))
  saveRDS(list(cf = h$cf, slide = h$slide, grp = h$grp, bc = h$bc, per = h$per,
               ref_bc_sub = ref_bc_sub, f_full = f_full, f_hold = f_hold,
               d_frac = d_frac, ref_med_hold = ref_med_hold, nwin = h$nwin),
          file.path(PDIR, "SC3_holdout.rds"))
  MANI$sc3 <- list(hold_frac = HOLD_FRAC, seed = SEED, n_full = M$n_ref,
                   n_hold = sum(vapply(ref_bc_sub, length, integer(1))),
                   f_obs_full = f_full, f_obs_hold = f_hold, d_frac = d_frac,
                   ref_median_hold = ref_med_hold, nwin = h$nwin)
  rm(Ho, h); invisible(gc())
}

## —————————————————————————————————————————————————————————————
## 九、D 对照锚（§19.10 第 3 项：**跑**；spot 取自该患者**全部切片**；观测集与主臂相同）
## —————————————————————————————————————————————————————————————
if (any(arms_wanted %in% c("all", "d"))) {
  step("---- D 对照锚：参考 = 该患者**全部切片**的非上皮 spot；观测集与主臂相同 ----")
  d_ref_slides <- own_slides
  d_ref_bc <- setNames(lapply(d_ref_slides, function(s) sort(nonepi[[s]])), d_ref_slides)
  d_ref_bc <- d_ref_bc[vapply(d_ref_bc, length, integer(1)) > 0L]
  stopifnot(length(d_ref_bc) > 0L)
  obs_bc   <- setNames(lapply(obs_slides, function(s) sort(epi[[s]])), obs_slides)
  D <- build_obj(names(d_ref_bc), d_ref_bc, obs_slides, obs_bc, "d_arm")
  step("  D 对象：%d 基因 x %d spot（非上皮参考 %d ＋ 观测 %d）",
       nrow(D$obj), ncol(D$obj), D$n_ref, D$n_obs)
  dd <- run_pass(D$obj, "D_arm")
  save_pass(dd, "D_arm", list(ref_slides = names(d_ref_bc),
                              ref_kind = "patient_nonepi_all_slides"))
  MANI$d_arm <- list(n_ref = D$n_ref, n_obs = D$n_obs,
                     ref_slides = names(d_ref_bc),
                     ref_kind = "该患者全部切片的非上皮 spot（§11.2 措辞）",
                     median_obs = median(dd$cf[dd$grp == "observation"], na.rm = TRUE),
                     tiers = as.list(tier_frac(dd$cf[dd$grp == "observation"])),
                     per_slide = dd$per)
  rm(D, dd); invisible(gc())
}

## —————————————————————————————————————————————————————————————
## 十、SC0 底座核对 ＋ SC4 上报（逐患者）
## —————————————————————————————————————————————————————————————
step("---- SC0 底座核对 ----")
go_path <- file.path(GODIR, "gene_order_spatial_hg38.tsv")
go_man  <- jsonlite::fromJSON(file.path(GODIR, "gene_order_spatial_manifest.json"))
go_sym  <- read.delim(go_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
h_rec   <- sha256(go_path)
gate(identical(h_rec, go_man$out_sha256), "SC0① 基因序哈希重算",
     sprintf("%s（清单 %s）", substr(h_rec, 1, 16), substr(go_man$out_sha256, 1, 16)))
our_sym <- rownames(cache[[own_slides[1]]])
mr <- mean(our_sym %in% go_sym)
gate(mr >= SC0_MIN_MATCH, "SC0① 基因序匹配率",
     sprintf("%.4f（%d/%d）", mr, sum(our_sym %in% go_sym), length(our_sym)))
ok_ge <- all(unlist(sc0_ge[need_slides]) == EXP_GE_ROWS)
gate(ok_ge, "SC0② GE 行数", sprintf("逐张 %s（期望 %d）",
     paste(unique(unlist(sc0_ge[need_slides])), collapse = ","), EXP_GE_ROWS))
ab_v <- unlist(sc0_ab[need_slides])
gate(all(ab_v %in% EXP_AB_SET), "SC0② 抗体行 ⊆ {0,35} 且已剔",
     sprintf("逐张 %s", paste(sprintf("%s:%d", names(ab_v), ab_v), collapse = " ")))
oth <- sort(unique(unlist(sc0_other[need_slides])))
gate(!length(oth), "SC0② 无第三类 feature",
     if (length(oth)) paste(oth, collapse = ",") else "只有 GE 与 Antibody Capture")
ok_spot <- all(vapply(own_slides, function(s) length(epi[[s]]) + length(nonepi[[s]]) <= n_pass[[s]], logical(1)))
gate(ok_spot, "SC0③ 上皮/非上皮 ⊆ 掩码 pass",
     sprintf("合计 pass %d ≥ 上皮+非上皮 %d", sum(unlist(n_pass[own_slides])),
             sum(vapply(own_slides, function(s) length(epi[[s]]) + length(nonepi[[s]]), integer(1)))))
gm    <- getGenes()
gm_f  <- gm[gm$gene_biotype %in% c("protein_coding", "lncRNA") &
            gm$chromosome_name %in% c(1:22, "X") & gm$hgnc_symbol != "", ]
tok   <- unique(gm_f$hgnc_symbol); in_tab <- our_sym %in% tok
gene_acct <- list(n_our_symbols = length(our_sym), n_matched = sum(in_tab),
                  n_unmatched = sum(!in_tab), unmatched_examples = head(our_sym[!in_tab], 15),
                  dropped_by_tool_rule = "chrY 与 chrM 被工具硬编码丢弃（§15.8 第 1 项）")
step("SC0④ fastCNV 基因账：%d − %d = %d（**只报不停**）",
     gene_acct$n_our_symbols, gene_acct$n_unmatched, gene_acct$n_matched)

if (!is.null(M)) {
  step("==== SC4 上报（主臂 %s；**不设通过率阈值、不筛 spot**）====", MAIN_KIND)
  obs_cf <- M$first$cf[M$first$grp == "observation"]
  ref_cf <- M$first$cf[M$first$grp == "reference"]
  for (i in seq_len(nrow(M$first$per)))
    step("  %-26s %-12s n=%5d 中位 %.4f [%.4f, %.4f] 最大 %.4f NA=%d",
         M$first$per$slide[i], M$first$per$group[i], M$first$per$n[i],
         M$first$per$median[i], M$first$per$q25[i], M$first$per$q75[i],
         M$first$per$max[i], M$first$per$n_na[i])
  tf <- tier_frac(obs_cf)
  step("  观测组四档比例：%s", paste(sprintf("%s=%.4f", names(tf), tf), collapse = "  "))
  step("  🔴 参考来源 = %s（%s）；参考组中位 %.4f（**循环性，不是证据**）",
       ref_kind_lab, if (MAIN_KIND == "self") "该患者自己的良性上皮" else "同性别供体汇池",
       median(ref_cf, na.rm = TRUE))
  ## 各染色体臂均值（观测组）
  arms_obs <- if (!is.null(M$first$arms)) {
    a <- M$first$arms
    a[a$group == unique(a$group[grepl("observation", a$group)])[1], , drop = FALSE]
  } else NULL
  MANI$sc4 <- list(main_kind = MAIN_KIND, ref_source = ref_kind_lab,
                   n_epi_obs = length(obs_cf), n_ref = M$n_ref, n_obs = length(obs_cf),
                   main_quantiles = unname(quantile(obs_cf, c(.25, .5, .75), na.rm = TRUE)),
                   median_ref = median(ref_cf, na.rm = TRUE),
                   tiers_obs = as.list(tf), per_slide = M$first$per,
                   arms_obs = arms_obs, contrib = M$contrib)
  MANI$sc2 <- list(n_repeat = N_REPEAT, frac_bit_identical = M$fi$frac,
                   anchor_names_identical = TRUE, run_min = round(M$run_min, 3))
}

## —————————————————————————————————————————————————————————————
## 十一、图（**英文标签**：本机无 CJK 字体）
## —————————————————————————————————————————————————————————————
if (!is.null(M)) {
  png(file.path(PDIR, sprintf("%s_main_%s.png", patient, MAIN_KIND)), width = 1700, height = 900, res = 130)
  op <- par(mar = c(10, 4.5, 3.5, 1))
  d <- data.frame(cf = M$first$cf, slide = M$first$slide, grp = M$first$grp)
  d <- d[d$grp == "observation", , drop = FALSE]
  if (nrow(d)) {
    boxplot(cf ~ slide, data = d, las = 2, outline = FALSE, col = "steelblue",
            ylab = "cnv_fraction",
            main = sprintf("%s: observation (lesion) vs %s reference", patient, MAIN_KIND),
            cex.main = 1.0, cex.axis = 0.8)
    abline(h = TH_MAIN, lty = 2, col = "grey30")
    legend("topleft", sprintf("%.2f (report line, not a filter)", TH_MAIN),
           lty = 2, col = "grey30", bty = "n", cex = 0.75)
  }
  par(op); invisible(dev.off())
}

## —————————————————————————————————————————————————————————————
## 十二、manifest
## —————————————————————————————————————————————————————————————
MANI$script <- "08_spatial_deconv/19_run_fastcnv_cohort.R"
MANI$caliber <- list(reference = "a（粗版 6 谱系）", epi = CAL_EPI, rule = CAL_RULE,
                     epi_def = "argmax == 上皮，零阈值，并列剔除",
                     nonepi_def = "argmax != 上皮，并列剔除")
MANI$params <- list(assay = F_ASSAY, windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
                    topNGenes = F_TOPN_GENES, thresholdPercentile = F_THRESH_PCT,
                    getCNVPerChromosomeArm = F_GET_ARM, n_repeat = N_REPEAT,
                    path = "CNVCalling + CNVPerChromosomeArm（单对象；**未**走 CNVCallingList）")
MANI$tool <- list(name = "fastCNV", version = as.character(packageVersion("fastCNV")),
                  license = "GPL-3",
                  source = "Cabrejas et al., Genome Medicine 2026, DOI 10.1186/s13073-026-01731-w")
MANI$sc0 <- list(hash_recomputed = h_rec, hash_in_manifest = go_man$out_sha256,
                 hash_ok = identical(h_rec, go_man$out_sha256),
                 gene_order_match_rate = round(mr, 4),
                 ge_rows = as.list(unlist(sc0_ge[need_slides])),
                 ab_rows = as.list(unlist(sc0_ab[need_slides])),
                 mask_pass = as.list(unlist(n_pass[own_slides])),
                 n_tie = as.list(unlist(n_tie)), gene_account = gene_acct)
MANI$coverage <- list(own_slides = own_slides, own_benign = own_benign, own_lesion = own_lesion,
                      ref_slides = ref_src_slides, ref_kind = ref_kind_lab,
                      pool_pats = pool_pats, pool_slides = pool_slides)
MANI$verdict <- if (length(FAILS)) "FAIL" else "PASS"
MANI$fails <- FAILS
MANI$not_claimed <- c("不产生恶性判定（cnv_fraction 是连续量；CNVClassification 只是方向分类器）",
                      "参考是良性/AAH 上皮，**不是**组织学正常上皮（§19.6 第 1 条）",
                      "参考组由构造≈0 ⇒「参考平」不是证据（§19.6 第 4 条）",
                      "不声称借锚与自带锚等价（参考组织与个体都不同）",
                      "不声称 25 例覆盖了完整病程轴（Normal 仍 1 张、AAH 观测 0 新增）",
                      "0.05 只用于报比例，**不是**筛 spot 的阈值（§19.7 SC4）")
MANI$elapsed_min_total <- round(elapsed(t_start), 2)
MANI$peak_gb <- round(peak_gb(), 3)
jsonlite::write_json(MANI, file.path(PDIR, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)

step("==== 患者 %s 完成：总用时 %.1f min；峰值 %.2f GB；SC 结论 %s%s ====",
     patient, elapsed(t_start), peak_gb(), MANI$verdict,
     if (length(FAILS)) paste0("（未过：", paste(FAILS, collapse = " / "), "）") else "")
step("产物：%s", PDIR)
