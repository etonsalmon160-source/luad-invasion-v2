#!/usr/bin/env Rscript
# 16_run_fastcnv_anchor_flatness.R —— 外部锚「自己平不平」的体检（SPATIAL_CNV_PREREG.md §17.5）
#
# 🔴 本脚本**不产生任何恶性判定**，**不替换** §3 的 infercnv 复现臂，**不进** SC0–SC4。
#    它只回答 §17.4 那一件事：外部对照锚（GSE248082 的 N1/N3）自身是不是平的。
#
# 规格来源：§17.4（用户 2026-09-29 签字：用已签的 fastCNV 比）+ §17.5（判据 F1–F4）
#   运行 A  参考 = 除 P4 外 9 例自带锚患者（全 AAH）的上皮 spot（每张抽稀至 ≤2,000）
#           观测 = P4_Normal / P4_AAH / P4_AAH-1 / P4_LUAD(阳性对照) / N1 / N3
#   运行 B  批次解耦（F4）：N1 作参考 ⇒ N3 作观测，及其反向
#
# 跑法（**只挂 fastcnv 这个库**，别和 infercnv 的库混，两边的 Seurat 版本不同）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 08_spatial_deconv/16_run_fastcnv_anchor_flatness.R

suppressMessages({
  library(fastCNV)
  library(Seurat)
  library(Matrix)
  library(rhdf5)
})

## —————————————————————————————————————————————————————————————
## 一、参数（逐项显式，法则 3.1；工具项与 §15.4 逐位相同，**一个都不改**）
## —————————————————————————————————————————————————————————————
CAL_REFERENCE <- "a"        # [已签 §11 第 1 项] 粗版 6 谱系
CAL_EPI       <- "上皮"     # [已签 §11 第 2 项] 单类目
CAL_RULE      <- "argmax"   # [已签 §12.1] 行内最大类 == 上皮，零阈值

F_WINDOW_SIZE    <- 150     # [工具默认，显式写出] ⚠️ 与 infercnv 的 101 不同
F_WINDOW_STEP    <- 10      # [工具默认，显式写出]
F_TOPN_GENES     <- 7000    # [工具默认，显式写出]
F_THRESH_PCT     <- 0.01    # [工具默认，显式写出] 背景噪声分位带
F_GET_ARM        <- TRUE    # 逐臂产物
F_ASSAY          <- "Spatial"

CAP_REF_PER_SLIDE <- 2000L  # [§17.5 执行前定] 每张参考切片上皮 spot 上限，seed=1
SEED              <- 1L

## 判据阈值（§17.5，**Claude 起草、可撤回**）
TH_F1_MEDIAN  <- 0.10   # 阳性对照 median 下限
TH_F1_AUC     <- 0.75   # 阳性对照 AUC 下限
TH_FLAT_MEDIAN<- 0.05   # 良性/外部 median 上限

ROOT   <- "/home/eto/luad_v2"
RES    <- file.path(ROOT, "results/08_spatial_deconv")
VISIUM <- file.path(ROOT, "data/visium_spatial")
EXTDIR <- file.path(ROOT, "data/external/GSE248082")
EXT_RES<- file.path(RES, "rctd_external_GSE248082")
OUTD   <- file.path(RES, "spatial_cnv/fastcnv_anchor_flatness")

## 参考：除 P4 外 9 例自带锚患者的良性切片（**全是 AAH**，§17.1 事实 2）
REF_SLIDES <- c("GSM9226168_P1_AAH", "GSM9226170_P2_AAH", "GSM9226180_P6_AAH",
                "GSM9226186_P9_AAH", "GSM9226191_P11_AAH", "GSM9226209_P20_AAH",
                "GSM9226214_P22_AAH", "GSM9226220_P24_AAH", "GSM9226222_P25_AAH")
## 观测：P4 四张（同队列）+ 外部两张
P4_SLIDES  <- c("GSM9226174_P4_Normal", "GSM9226175_P4_AAH",
                "GSM9226176_P4_AAH-1", "GSM9226177_P4_LUAD")
EXT_SLIDES <- c("N1", "N3")

## —————————————————————————————————————————————————————————————
## 二、工具函数
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
## 🔴 `step()` 是纯日志 ⇒ **绝不能因为格式串的问题把 30 分钟的计算打死**
##    （第 3 次跑就死在这：`sprintf` 报 "required resulting string length 888629
##      is greater than maximal 8192"，fmt 被截断到 8191 ⇒ 打印本身出了问题）
##    所以：任何格式化异常都降级成一行纯文本，并**把 fmt 的长度打出来**定位元凶
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

## 无阈值的分离度：AUC = P(x > y)，并列各算 0.5。**只报**
auc_gt <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  r <- rank(c(x, y))
  (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y))
}

stage_of <- function(s) sub("^[^_]*_[^_]*_", "", s)

## 外部切片：10x v3 h5（CSC，**细胞为列**）—— 与 13_run_rctd_external.R 同法
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

## 按 symbol 求和（双探针 HSPA14 / TBCE / TMSB15B）
agg_by_symbol <- function(m, syms) {
  us <- unique(syms)
  if (length(us) == length(syms) && !anyDuplicated(syms)) {
    rownames(m) <- us; return(m)
  }
  agg <- sparseMatrix(i = match(syms, us), j = seq_along(syms), x = 1,
                      dims = c(length(us), length(syms)))
  m2 <- as(agg %*% m, "CsparseMatrix"); rownames(m2) <- us; m2
}

## 上皮 spot（口径 argmax == 上皮，零阈值，并列剔除；与 11_ / 12_ 脚本逐字同源）
epi_of <- function(slide, wf) {
  w <- read.delim(gzfile(wf), row.names = 1, check.names = FALSE)
  stopifnot(CAL_EPI %in% colnames(w), ncol(w) == 6L)
  W    <- as.matrix(w)
  mx   <- do.call(pmax, as.data.frame(W))
  n_at <- rowSums(W >= mx)
  am   <- colnames(W)[max.col(W, ties.method = "first")]
  list(epi = rownames(W)[am == CAL_EPI & n_at <= 1L], n_all = nrow(W), n_tie = sum(n_at > 1L))
}

step("==== 外部锚体检开始（§17.5）====")
step("工具 fastCNV %s / Seurat %s（Genome Medicine 2026-08-29 已发表，§15.2 更正）",
     as.character(packageVersion("fastCNV")), as.character(packageVersion("Seurat")))
step("🔴 本脚本**不产生恶性判定**，不得用于 SC0–SC4")

## —————————————————————————————————————————————————————————————
## 三、上皮 spot 名单（读已签的 RCTD 权重口径；**不重算 RCTD**）
## —————————————————————————————————————————————————————————————
ref_keep  <- list(); obs_keep <- list()

step("---- 参考（9 例自带锚患者，全 AAH；每张上限 %d spot，seed=%d）----", CAP_REF_PER_SLIDE, SEED)
set.seed(SEED)
ref_sampled <- list()
for (s in REF_SLIDES) {
  e <- epi_of(s, file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s)))
  stopifnot(length(e$epi) > 0L)
  pick <- if (length(e$epi) > CAP_REF_PER_SLIDE) sort(e$epi)[sort(sample.int(length(e$epi), CAP_REF_PER_SLIDE))] else sort(e$epi)
  ref_keep[[s]]     <- pick
  ref_sampled[[s]]  <- list(n_epi_total = length(e$epi), n_used = length(pick), cap = CAP_REF_PER_SLIDE)
  step("  %-24s 上皮 %6d ⇒ 取 %6d", s, length(e$epi), length(pick))
}

step("---- 观测（P4 四张，**不抽稀**）----")
for (s in P4_SLIDES) {
  e <- epi_of(s, file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s)))
  stopifnot(length(e$epi) > 0L)
  obs_keep[[s]] <- sort(e$epi)
  step("  %-24s 上皮 %6d（并列剔除 %d）", s, length(e$epi), e$n_tie)
}

step("---- 观测（外部 N1/N3，**不抽稀**；名单取自已跑完的外部 RCTD）----")
for (s in EXT_SLIDES) {
  ef <- file.path(EXT_RES, sprintf("per_slide/%s.a.epi_barcodes.txt.gz", s))
  if (!file.exists(ef))
    stop(sprintf("%s 的上皮名单尚未产出（%s 不存在）⇒ 硬停，等外部 RCTD 跑完再跑本体检",
                 s, ef), call. = FALSE)
  obs_keep[[s]] <- sort(read_gz(ef, readLines))
  step("  %-24s 上皮 %6d", s, length(obs_keep[[s]]))
}
n_ref_used <- sum(lengths(ref_keep))
step("参考 spot %d ⇒ 观测 spot %d（合计 %d）", n_ref_used,
     sum(lengths(obs_keep)), n_ref_used + sum(lengths(obs_keep)))

## —————————————————————————————————————————————————————————————
## 四、读表达矩阵（Gene Expression；双探针按 symbol 求和；只留上皮 spot）
## —————————————————————————————————————————————————————————————
t_read <- Sys.time()
mats <- list(); cell_slide <- character(0); g_a <- character(0)
all_slides <- c(REF_SLIDES, P4_SLIDES, EXT_SLIDES)

for (s in all_slides) {
  is_ref <- s %in% REF_SLIDES
  keep   <- if (is_ref) ref_keep[[s]] else obs_keep[[s]]
  lab    <- if (is_ref) "cohort_benign_epi" else "obs_other"

  if (s %in% EXT_SLIDES) {                       # —— 外部：h5 ——
    gsm <- if (s == "N1") "GSM8087031" else "GSM8087033"
    h5c <- read_h5_counts(file.path(EXTDIR, sprintf("%s_%s_filtered_feature_bc_matrix.h5", gsm, s)))
    if (h5c$n_ab_rows != 0L)
      stop(s, ": 外部切片出现抗体行，与 §16.2 的核验不符 ⇒ 硬停", call. = FALSE)
    ge  <- h5c$types == "Gene Expression"
    m   <- agg_by_symbol(as(h5c$counts[ge, , drop = FALSE], "CsparseMatrix"), h5c$symbols[ge])
    bc  <- colnames(m)
    rm(h5c); invisible(gc())
  } else {                                        # —— 本队列：mtx 目录 ——
    d      <- file.path(VISIUM, s, "filtered_feature_bc_matrix")
    ft_all <- read_gz(file.path(d, "features.tsv.gz"), function(con)
      read.delim(con, header = FALSE, col.names = c("id", "sym", "type"), stringsAsFactors = FALSE))
    ab <- which(ft_all$type == "Antibody Capture")
    stopifnot(length(ab) %in% c(0L, 35L))
    bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
    m  <- readMM(file.path(d, "matrix.mtx.gz"))
    stopifnot(nrow(m) == nrow(ft_all), ncol(m) == length(bc))
    ## 🔴 `readMM()` **不设** dimnames（.mtx 里根本没有）⇒ 必须显式补列名，
    ##    否则后面 match(keep, colnames(m)) 拿到的 colnames 是 NULL，逐张全 0 匹配
    colnames(m) <- bc
    ## ⚠️ `x[-integer(0), ]` 返回**零行** ⇒ 必须 if 包住
    if (length(ab)) { m <- m[-ab, , drop = FALSE]; ft <- ft_all[-ab, , drop = FALSE] } else { ft <- ft_all }
    stopifnot(nrow(ft) == 18085L, all(ft$type == "Gene Expression"))
    m <- agg_by_symbol(as(m, "CsparseMatrix"), ft$sym)
  }

  use <- match(keep, colnames(m), nomatch = 0L)
  stopifnot(!any(use == 0L), !anyDuplicated(use))
  m   <- m[, use, drop = FALSE]
  colnames(m) <- paste(s, keep, sep = "_")

  mats[[s]]  <- m
  cell_slide <- c(cell_slide, setNames(rep(s, ncol(m)), colnames(m)))
  g_a        <- c(g_a, rep(lab, ncol(m)))
  step("  %-24s %s 上皮 %6d  组=%s", s, if (s %in% EXT_SLIDES) "[外部]" else "[队列]", ncol(m), lab)
  rm(m); invisible(gc())
}

## 对齐到同一 symbol 序（外部与本队列探针集相同，但仍显式对齐，防顺序差异）
sym_master <- Reduce(union, lapply(mats, rownames))
for (s in names(mats)) {
  m <- mats[[s]]
  add <- setdiff(sym_master, rownames(m))
  if (length(add)) {
    z <- Matrix::Matrix(0, nrow = length(add), ncol = ncol(m), sparse = TRUE)
    rownames(z) <- add
    m <- rbind(m, z)
  }
  mats[[s]] <- m[sym_master, , drop = FALSE]
}
n_missing_tot <- sum(vapply(mats, function(m) sum(Matrix::colSums(m) == 0), numeric(1)))
mat <- do.call(cbind, mats); rm(mats); invisible(gc())
stopifnot(ncol(mat) == length(g_a), identical(colnames(mat), names(cell_slide)))
step("对齐后矩阵 %d symbol x %d spot（全零列 %d；用时 %.1f min）",
     nrow(mat), ncol(mat), n_missing_tot, elapsed(t_read))

obj <- CreateSeuratObject(counts = mat, assay = F_ASSAY,
                          meta.data = data.frame(g_a = g_a, slide = unname(cell_slide),
                                                 row.names = colnames(mat)))
rm(mat); invisible(gc())
step("Seurat 对象就绪：%d 基因 x %d spot（峰值 %.2f GB）", nrow(obj), ncol(obj), peak_gb())

## —————————————————————————————————————————————————————————————
## 五、运行 A / B（单对象 CNVCalling；**不走**汇池路径，见 §15.8 第 6 项）
## —————————————————————————————————————————————————————————————
gm <- getGenes()
run_one <- function(o, tag, var, lab) {
  step("---- %s：referenceVar=%s，referenceLabel=%s ----", tag, var, lab)
  t0 <- Sys.time()
  r <- CNVCalling(o, assay = F_ASSAY, referenceVar = var, referenceLabel = lab,
                  scaleOnReferenceLabel = TRUE, thresholdPercentile = F_THRESH_PCT,
                  geneMetadata = gm, windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
                  saveGenomicWindows = FALSE, topNGenes = F_TOPN_GENES)
  if (F_GET_ARM) r <- CNVPerChromosomeArm(r)
  d <- elapsed(t0)
  RUN_MIN[[length(RUN_MIN) + 1L]] <<- setNames(d, tag)
  step("---- %s 完成，用时 %.2f min（峰值 %.2f GB）----", tag, d, peak_gb())
  r
}
## 🔴 运行 A 的元数据列叫 `g_a`，运行 B 的叫 `g_B` ⇒ 必须把列名传进来，
##    不能写死 `o$g_a`（第 2 跑就死在这一行：'g_a' not found，Did you mean "g_B"?）
## 🔴🔴 第 3 跑的真正元凶（**已定位，不是猜**）：`o[[grp_col]]` 取回来的是
##    **一个单列 data.frame**，`as.character()` 会把**整列**deparse 成**一个**
##    888,629 字的字符串（内容形如 `c("cohort_benign_epi", "cohort_benign_epi", …)`）；
##    长度 1 又被 `data.frame()` **静默循环**贴到全部 57,167 行上。
##    这个巨串进了 `sprintf("%-18s", …)` ⇒ `required resulting string length 888629
##    is greater than maximal 8192`，把纯日志打印打崩（`o$slide` 用的是 `$` 所以是好的）。
##    修法：① 直接读 `o@meta.data`（绕开 `[[` 的分派）；② 加**长度守卫**，形状不对就硬停。
getv <- function(o, tag, grp_col = "g_a") {
  md <- o@meta.data
  for (k in c("slide", grp_col)) if (!k %in% names(md)) stop("元数据缺列：", k)
  n  <- ncol(o)
  sl <- as.character(md[["slide"]]); gp <- as.character(md[[grp_col]])
  cfv <- FetchData(o, vars = "cnv_fraction")[[1]]
  if (length(sl) != n || length(gp) != n || length(cfv) != n)
    stop(sprintf("getv 长度不齐：ncol=%d slide=%d %s=%d cf=%d", n, length(sl), grp_col, length(gp), length(cfv)))
  if (any(nchar(gp) > 100)) stop("getv：分组列出现超长字符串（疑似被 deparse），拒绝继续")
  data.frame(run = tag, slide = sl, grp = gp, cf = cfv, row.names = NULL)
}
## 分组名一律 `paste(..., collapse="|")` ⇒ **保证是长度 1 的短串**
## （原写法 `group = unique(...)` 一旦 unique 返回长度>1，`data.frame` 会静默循环/报错）
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

o_a <- run_one(obj, "A 参考=9例自带锚(全AAH)上皮", "g_a", "cohort_benign_epi")
va  <- getv(o_a, "A"); tab_a <- collect(va, "A")
nwin <- nrow(GetAssay(o_a, assay = "genomicScores"))
rm(o_a); invisible(gc())
## 🔴 运行 A 的逐 spot 分数**立刻落盘**：前两次跑都是跑完 A 才在后面崩，
##    结果全丢、白烧 10 分钟。落盘后即使后面再崩，也不用重算 A。
saveRDS(list(va = va, tab_a = tab_a, nwin = nwin),
        file.path(OUTD, "runA_spot_scores.rds"))
step("运行 A 逐 spot 分数已落盘：%s", file.path(OUTD, "runA_spot_scores.rds"))

step("==== 运行 A 逐切片 cnv_fraction ====")
for (i in seq_len(nrow(tab_a)))
  step("  %-24s %-18s n=%5d 中位 %.4f [%.4f, %.4f] 最大 %.4f NA=%d",
       tab_a$slide[i], tab_a$group[i], tab_a$n[i], tab_a$median[i],
       tab_a$q25[i], tab_a$q75[i], tab_a$max[i], tab_a$n_na[i])

## 判据 F1 / F2 / F3（§17.5）
cf <- function(s) va$cf[va$slide == s]
ref_cf <- va$cf[va$grp == "cohort_benign_epi"]
f1_med <- median(cf("GSM9226177_P4_LUAD"), na.rm = TRUE)
f1_auc <- auc_gt(cf("GSM9226177_P4_LUAD"), ref_cf)
f2 <- c(P4_Normal   = median(cf("GSM9226174_P4_Normal"), na.rm = TRUE),
        P4_AAH      = median(cf("GSM9226175_P4_AAH"),    na.rm = TRUE),
        "P4_AAH-1"  = median(cf("GSM9226176_P4_AAH-1"),  na.rm = TRUE))
f3 <- c(N1 = median(cf("N1"), na.rm = TRUE), N3 = median(cf("N3"), na.rm = TRUE))
F1 <- (f1_med >= TH_F1_MEDIAN) && (f1_auc >= TH_F1_AUC)
F2 <- all(f2 <= TH_FLAT_MEDIAN)
F3 <- all(f3 <= TH_FLAT_MEDIAN)

step("==== 判据（§17.5；阈值由 Claude 起草、可撤回）====")
step("  F1 阳性对照：median(P4_LUAD)=%.4f (需 ≥%.2f)；AUC(vs 参考)=%.4f (需 ≥%.2f) ⇒ %s",
     f1_med, TH_F1_MEDIAN, f1_auc, TH_F1_AUC, if (isTRUE(F1)) "过" else "🔴 不过")
step("     ⚠️ 该 AUC 与参考组比较，参考组分数由构造 ≈0 ⇒ **偏循环**；承重的是 median")
step("  F2 同队列良性：%s（各需 ≤%.2f）⇒ %s",
     paste(sprintf("%s=%.4f", names(f2), f2), collapse="  "), TH_FLAT_MEDIAN,
     if (isTRUE(F2)) "过" else "🔴 不过 ⇒ **全盘作废**，不得读成「外部锚不平」")
step("  F3 外部：%s（各需 ≤%.2f）⇒ %s",
     paste(sprintf("%s=%.4f", names(f3), f3), collapse="  "), TH_FLAT_MEDIAN,
     if (isTRUE(F3)) "过" else "🔴 不过 ⇒ 外部锚不得当主锚（可仍作对照）；须看 F4")

## —————————————————————————————————————————————————————————————
## 六、运行 B（F4 批次解耦）：只用外部两张，互为参考
## —————————————————————————————————————————————————————————————
rm(obj); invisible(gc())        # 运行 B 只用 obs_keep 里的名单 ⇒ 大对象可先放掉
runB <- list(); n4 <- NULL
vb_n1 <- NULL; vb_n3 <- NULL
if (all(EXT_SLIDES %in% names(obs_keep))) {
  ## 用一个只含 N1/N3 的小对象（省内存；判定对象与运行 A 里那两张**同一批 spot**）
  e_slides <- EXT_SLIDES
  matsB <- list(); gB <- character(0); sB <- character(0)
  for (s in e_slides) {
    gsm <- if (s == "N1") "GSM8087031" else "GSM8087033"
    h5c <- read_h5_counts(file.path(EXTDIR, sprintf("%s_%s_filtered_feature_bc_matrix.h5", gsm, s)))
    ge  <- h5c$types == "Gene Expression"
    m   <- agg_by_symbol(as(h5c$counts[ge, , drop = FALSE], "CsparseMatrix"), h5c$symbols[ge])
    rm(h5c); invisible(gc())
    use <- match(obs_keep[[s]], colnames(m), nomatch = 0L)
    m   <- m[, use, drop = FALSE]
    colnames(m) <- paste(s, obs_keep[[s]], sep = "_")
    matsB[[s]] <- m
    gB <- c(gB, rep(s, ncol(m)))
    sB <- c(sB, setNames(rep(s, ncol(m)), colnames(m)))
    rm(m); invisible(gc())
  }
  symB <- Reduce(union, lapply(matsB, rownames))
  for (s in names(matsB)) {
    m <- matsB[[s]]; add <- setdiff(symB, rownames(m))
    if (length(add)) {
      z <- Matrix::Matrix(0, nrow = length(add), ncol = ncol(m), sparse = TRUE)
      rownames(z) <- add; m <- rbind(m, z)
    }
    matsB[[s]] <- m[symB, , drop = FALSE]
  }
  matB <- do.call(cbind, matsB); rm(matsB); invisible(gc())
  objB <- CreateSeuratObject(counts = matB, assay = F_ASSAY,
                             meta.data = data.frame(g_B = gB, slide = unname(sB),
                                                    row.names = colnames(matB)))
  rm(matB); invisible(gc())
  step("运行 B 对象：%d 基因 x %d spot（N1 ↔ N3 互为参考）", nrow(objB), ncol(objB))

  for (lab in EXT_SLIDES) {
    oB <- run_one(objB, sprintf("B 参考=%s 上皮", lab), "g_B", lab)
    vb <- getv(oB, "B", "g_B"); vb$ref <- lab
    tabB <- collect(vb, sprintf("B_ref_%s", lab))
    step("  运行 B（参考=%s）:", lab)
    for (i in seq_len(nrow(tabB)))
      step("    %-6s n=%5d 中位 %.4f [%.4f, %.4f] 最大 %.4f",
           tabB$slide[i], tabB$n[i], tabB$median[i], tabB$q25[i], tabB$q75[i], tabB$max[i])
    runB[[lab]] <- tabB
    if (lab == "N1") vb_n1 <- vb else vb_n3 <- vb
    rm(oB, vb); invisible(gc())
    saveRDS(list(runB = runB, vb_n1 = vb_n1, vb_n3 = vb_n3),
            file.path(OUTD, "runB_spot_scores.rds"))   # 同样立刻落盘
  }
  ## F4：每一向里，被留出的那一张必须平
  f4 <- vapply(EXT_SLIDES, function(lab) {
    held <- setdiff(EXT_SLIDES, lab)
    median(runB[[lab]]$median[runB[[lab]]$slide %in% held], na.rm = TRUE)
  }, numeric(1))
  n4 <- f4
  step("==== F4 批次解耦（被留出的那一张，median 需 ≤%.2f）====", TH_FLAT_MEDIAN)
  for (lab in EXT_SLIDES)
    step("  参考=%s ⇒ 留出 %s 中位 %.4f  %s", lab, setdiff(EXT_SLIDES, lab),
         f4[[lab]], if (isTRUE(f4[[lab]] <= TH_FLAT_MEDIAN)) "过" else "🔴 不过")
} else {
  step("🔴 运行 B 跳过：外部上皮名单不全")
}
F4 <- if (is.null(n4)) NA else all(n4 <= TH_FLAT_MEDIAN)

## —————————————————————————————————————————————————————————————
## 七、出图（英文标签：本机无 CJK 字体）
## —————————————————————————————————————————————————————————————
fig <- file.path(OUTD, "fastcnv_anchor_flatness.png")
png(fig, width = 1900, height = 900, res = 130)
op <- par(mfrow = c(1, 2), mar = c(9, 4.5, 3.5, 1))
lev <- c(REF_SLIDES, P4_SLIDES, EXT_SLIDES)
d <- va[va$slide %in% lev, ]
d$slide <- factor(d$slide, levels = lev)
colv <- ifelse(lev %in% REF_SLIDES, "grey75",
        ifelse(lev %in% c("GSM9226174_P4_Normal", "GSM9226175_P4_AAH", "GSM9226176_P4_AAH-1"), "steelblue",
        ifelse(lev == "GSM9226177_P4_LUAD", "indianred", "darkorange")))
boxplot(cf ~ slide, data = d, col = colv, las = 2, outline = FALSE,
        ylab = "cnv_fraction", main = "Run A: reference = 9 own-anchor AAH epithelia",
        cex.main = 0.95, cex.axis = 0.75)
abline(h = TH_FLAT_MEDIAN, lty = 2, col = "grey30")
abline(h = TH_F1_MEDIAN,   lty = 3, col = "grey30")
legend("topleft", c("reference (by construction ~0)", "same-cohort benign",
                    "positive control (P4_LUAD)", "external (N1/N3)",
                    sprintf("%.2f flat line", TH_FLAT_MEDIAN)),
       fill = c("grey75", "steelblue", "indianred", "darkorange", NA), border = NA,
       lty = c(NA, NA, NA, NA, 2), cex = 0.65)
if (!is.null(vb_n1)) {
  d2 <- rbind(vb_n1, if (!is.null(vb_n3)) vb_n3)
  d2$grp2 <- paste(d2$slide, "| ref=", ifelse(d2$ref == "N1", "N1", "N3"))
  boxplot(cf ~ grp2, data = d2, col = rep(c("darkorange", "gold"), each = 2), las = 2,
          outline = FALSE, ylab = "cnv_fraction",
          main = "Run B (F4): external-internal, N1 <-> N3",
          cex.main = 0.95, cex.axis = 0.75)
  abline(h = TH_FLAT_MEDIAN, lty = 2, col = "grey30")
}
par(op); invisible(dev.off())
step("图已出：%s", fig)

## —————————————————————————————————————————————————————————————
## 八、manifest
## —————————————————————————————————————————————————————————————
manifest <- list(
  script = "08_spatial_deconv/16_run_fastcnv_anchor_flatness.R",
  purpose = "SPATIAL_CNV_PREREG.md §17.5：外部对照锚（GSE248082 N1/N3）自身平不平的体检",
  not_producing = c("不产生恶性判定", "不替换 §3 infercnv 复现臂", "不进 SC0-SC4"),
  tool = list(name = "fastCNV", version = as.character(packageVersion("fastCNV")),
              license = "GPL-3",
              source = "Cabrejas et al., Genome Medicine 2026, DOI 10.1186/s13073-026-01731-w (published online 2026-08-29; preprint bioRxiv 2025.10.22.683855)"),
  params = list(assay = F_ASSAY, windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
                topNGenes = F_TOPN_GENES, thresholdPercentile = F_THRESH_PCT,
                getCNVPerChromosomeArm = F_GET_ARM,
                path = "CNVCalling + CNVPerChromosomeArm（单对象）",
                cap_ref_per_slide = CAP_REF_PER_SLIDE, seed = SEED),
  caliber = list(reference = CAL_REFERENCE, epi = CAL_EPI, rule = CAL_RULE,
                 ref_slides = REF_SLIDES, p4_slides = P4_SLIDES, ext_slides = EXT_SLIDES),
  ref_sampling = ref_sampled,
  thresholds = list(f1_median = TH_F1_MEDIAN, f1_auc = TH_F1_AUC, flat_median = TH_FLAT_MEDIAN,
                    drafted_by = "Claude（可撤回）",
                    provenance = "0.05 ≈ §15.9 良性上限 0.0046-0.0138 的 4 倍；0.75 = 明显分得开；0.10 ≈ §15.9 LUAD 中位 0.5092 的 1/5"),
  per_slide_runA = tab_a,
  criteria = list(F1_pipeline_valid = F1, F2_same_cohort_benign_flat = F2,
                  F3_external_flat = F3, F4_batch_decoupled = F4),
  detail = list(f1_median_luad = f1_med, f1_auc_vs_ref = f1_auc,
                f2_medians = as.list(f2), f3_medians = as.list(f3),
                f4_heldout_medians = as.list(n4)),
  per_slide_runB = if (length(runB)) do.call(rbind, runB) else NULL,
  circularity_note = "参考组 cnv_fraction 由构造 ≈0 ⇒ 参考组平不是证据（§15.3）",
  limitation = "本设计分不开「不平」与「批次效应」；运行 B 只缩小不解除（§17.5 失效规则 4）",
  v1 = list(elapsed_min_by_run = round(unlist(RUN_MIN), 2), peak_gb_final = round(peak_gb(), 3),
            n_windows = nwin),
  figure = fig,
  elapsed_min_total = round(elapsed(t_start), 2)
)
jsonlite::write_json(manifest, file.path(OUTD, "fastcnv_anchor_flatness_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)
step("==== 体检完成：总用时 %.1f min；峰值 %.2f GB ====", elapsed(t_start), peak_gb())
step("产物：%s", OUTD)
