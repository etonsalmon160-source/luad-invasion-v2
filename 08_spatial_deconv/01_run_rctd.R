#!/usr/bin/env Rscript
## 01_run_rctd.R —— M5 空间解卷积主脚本（spacexr::RCTD，doublet_mode="full"）
## ============================================================================
## 口径来源（全部已签字，本脚本**不引入任何新参数**）
##   RCTD_PREREG.md    §2 输入定义、§3 参数、§4 守卫 G1–G5、§10/§11 签字
##   calibration.json  signed=TRUE 才准跑（未签 → 硬停）
##   SPATIAL_CNV_PREREG.md §11 门控口径（本脚本只产权重，不做门控）
##
## 本脚本**不做**的事（边界，别越）
##   - 不做恶性判定。权重不是恶性标签（RCTD_PREREG §1.1 / R2 边界）。
##   - 不做门控（取上皮 spot 子集的是下一支脚本；口径已签为 `argmax == 上皮`，
##     见 SPATIAL_CNV_PREREG §12。**不是**「上皮权重 > 0.5」，那条已作废）。
##   - 不调用 spacexr 的行归一化函数（§M5 / G4 明文禁止）：本文件**不含**该函数名。
##
## 用法
##   Rscript 08_spatial_deconv/01_run_rctd.R --caliber a                  # 全 56 张
##   Rscript 08_spatial_deconv/01_run_rctd.R --caliber a --smoke GSM9226168_P1_AAH
##   Rscript 08_spatial_deconv/01_run_rctd.R --caliber a --max-cores 4
##   # 尺度不变性实验（RCTD_PREREG §13.3）：计数减半重跑，只用于单张
##   Rscript 08_spatial_deconv/01_run_rctd.R --caliber a --smoke GSM9226168_P1_AAH \
##       --thin 0.5 --out results/08_spatial_deconv/rctd_a_thin50
## ============================================================================

suppressMessages({
  library(Matrix)
  library(spacexr)
  library(rhdf5)
  library(jsonlite)
  library(data.table)
})

## ---- 路径常量：数据与结果在原仓，脚本在工作树（口径与代码分开管） ----------
ROOT     <- "/home/eto/luad_v2"
SUP_DIR  <- file.path(ROOT, "data", "visium_spatial")
RES_DIR  <- file.path(ROOT, "results", "08_spatial_deconv")

## ---- 本脚本自身路径（G4 自检要读自己的源码） -------------------------------
self_path <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) != 1L) NA_character_ else normalizePath(f, mustWork = FALSE)
}
SELF <- self_path()
SCRIPT_DIR <- if (is.na(SELF)) file.path(ROOT, "08_spatial_deconv") else dirname(SELF)

## ============================================================================
## 0. 参数解析
## ============================================================================
parse_args <- function(a) {
  out <- list(caliber = NULL, smoke = NA_character_, max_cores = 4L,
              seed = 1L, out = NULL, slides = NULL, thin = NA_real_)
  i <- 1L
  while (i <= length(a)) {
    k <- a[i]
    if (k == "--caliber")      { out$caliber <- a[i + 1L]; i <- i + 2L }
    else if (k == "--smoke")   { out$smoke   <- a[i + 1L]; i <- i + 2L }
    else if (k == "--max-cores") { out$max_cores <- as.integer(a[i + 1L]); i <- i + 2L }
    else if (k == "--seed")    { out$seed <- as.integer(a[i + 1L]); i <- i + 2L }
    else if (k == "--out")     { out$out <- a[i + 1L]; i <- i + 2L }
    else if (k == "--slides")  { out$slides <- strsplit(a[i + 1L], ",")[[1]]; i <- i + 2L }
    else if (k == "--thin")    { out$thin <- as.numeric(a[i + 1L]); i <- i + 2L }
    else stop("未知参数: ", k)
  }
  if (!is.na(out$thin) && (out$thin <= 0 || out$thin > 1)) {
    stop("--thin 必须落在 (0, 1] 区间", call. = FALSE)
  }
  out
}
ARGS <- parse_args(commandArgs(TRUE))

## ============================================================================
## 1. G4 守卫：源码自检（在**任何**计算之前）
## ============================================================================
## §M5 禁止把「逐 spot 行归一化」当成门（会把低质量 spot 里的污染装成干净混合）。
## 守卫写成可执行的：本文件一旦出现该函数名即硬停。
## 注：断言用的字符串在此拼接，避免自检触发自身。
FORBIDDEN <- paste0("normalize", "_weights", "(", ")")
if (!is.na(SELF)) {
  src <- readLines(SELF, warn = FALSE)
  hit <- grep(FORBIDDEN, src, fixed = TRUE)
  if (length(hit)) {
    stop("G4 越界：本脚本第 ", paste(hit, collapse = ","),
         " 行出现了被 §M5 明文禁止的函数调用 ⇒ 硬停。", call. = FALSE)
  }
}
cat("[G4] 源码自检通过：未出现被禁止的行归一化调用\n")

## ============================================================================
## 2. 口径签字校验（未签不许算 —— 法则「口径未定就不许启动计算」）
## ============================================================================
CAL_PATH <- file.path(SCRIPT_DIR, "calibration.json")
if (!file.exists(CAL_PATH)) stop("找不到 calibration.json: ", CAL_PATH, call. = FALSE)
CAL <- fromJSON(CAL_PATH, simplifyVector = FALSE)
if (!isTRUE(CAL$signed)) stop("calibration.json 的 signed 不为 TRUE ⇒ 口径未签 ⇒ 硬停", call. = FALSE)

CALIBER <- if (is.null(ARGS$caliber)) CAL$coarse_caliber else ARGS$caliber
if (!CALIBER %in% names(CAL$expect_n_types))
  stop("caliber 必须是 ", paste(names(CAL$expect_n_types), collapse = "/"), call. = FALSE)
EXPECT_N_TYPES <- as.integer(CAL$expect_n_types[[CALIBER]])

REF_H5   <- file.path(RES_DIR, paste0("reference_", CALIBER, ".h5"))
REF_MANI <- file.path(RES_DIR, paste0("reference_", CALIBER, ".manifest.json"))
if (!file.exists(REF_H5))   stop("参考不存在（该 caliber 尚未构建）: ", REF_H5, call. = FALSE)
if (!file.exists(REF_MANI)) stop("参考 manifest 不存在: ", REF_MANI, call. = FALSE)

OUT_DIR <- if (is.null(ARGS$out)) file.path(RES_DIR, paste0("rctd_", CALIBER)) else ARGS$out
PER_DIR <- file.path(OUT_DIR, "per_slide")
dir.create(PER_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[口径] caliber=%s  参考=%s  期望类型数=%d  输出=%s\n",
            CALIBER, basename(REF_H5), EXPECT_N_TYPES, OUT_DIR))

## §3 参数：逐项显式常量化（法则 3.1）。一个都不改。
P <- list(
  doublet_mode     = "full",   # run.RCTD()
  UMI_min          = 100,      # create.RCTD()
  CELL_MIN_INSTANCE = 25,      # create.RCTD()
  counts_MIN       = 10,       # create.RCTD()
  max_cores        = ARGS$max_cores,
  test_mode        = FALSE,
  gene_cutoff      = 0.000125,
  fc_cutoff        = 0.5,
  gene_cutoff_reg  = 0.0002,
  fc_cutoff_reg    = 0.75,
  UMI_max          = 2e7,
  UMI_min_sigma    = 300,
  MAX_MULTI_TYPES  = 4,
  keep_reference   = FALSE,
  CONFIDENCE_THRESHOLD = 5,
  DOUBLET_THRESHOLD    = 20,
  n_max_cells      = 10000,    # Reference()；§2.4
  min_UMI          = 100,      # Reference()
  require_int      = TRUE      # Reference()
)
cat(sprintf("[参数] max_cores=%d  seed=%d  doublet_mode=%s\n",
            P$max_cores, ARGS$seed, P$doublet_mode))
## 本脚本的 n_max_cells 必须与签字文件一致，防两边悄悄分叉
if (!is.null(CAL$n_max_cells) && as.integer(CAL$n_max_cells) != P$n_max_cells)
  stop(sprintf("脚本 n_max_cells=%d ≠ calibration.json %s ⇒ 硬停",
               P$n_max_cells, CAL$n_max_cells), call. = FALSE)

set.seed(ARGS$seed)   # Reference() 抽样用；seed 一并写进 manifest

## ============================================================================
## 3. 参考读取（h5：CSC，细胞为列）
## ============================================================================
## 类型名的工具适配（🔴 2026-09-25 冒烟测试实测，见 §12）
## spacexr::Reference() 的 check_cell_types() 明文拒绝含 "/" 的类型名
## （源码：prohibited_character <- "/"）。我们的谱系/亚型名是从
## results/05_annotation/ 的注释表**逐字搬来**的，其中带 "/"（粗版 2 个、
## 细版 4 个）。h5 保持与注释表逐字一致（可审计），换名**只在送入 RCTD 这一层**做，
## 并把映射表落盘，供下游还原。
sanitize_type_names <- function(x) {
  u <- unique(x)
  bad <- grepl("/", u, fixed = TRUE)
  if (!any(bad)) return(list(x = x, map = NULL))
  u2 <- gsub("/", "_", u, fixed = TRUE)
  if (anyDuplicated(u2)) stop("类型名去斜杠后出现重名 ⇒ 映射不唯一，硬停", call. = FALSE)
  map <- data.frame(original = u[bad], adapted = u2[bad], row.names = NULL)
  list(x = u2[match(x, u)], map = map)
}

read_reference <- function(h5) {
  att       <- h5readAttributes(h5, "counts")
  shape     <- as.integer(att$shape)                 # c(基因, 细胞)
  genes_ref <- as.character(h5read(h5, "genes"))
  ctypes    <- as.character(h5read(h5, "cell_types"))
  data_i    <- h5read(h5, "counts/data")
  ind_i     <- h5read(h5, "counts/indices")
  ptr_i     <- h5read(h5, "counts/indptr")
  data_i <- as.double(data_i)
  ind_i  <- as.integer(ind_i) + 1L                   # 0-based → 1-based
  ptr_i  <- as.integer(ptr_i)
  m <- sparseMatrix(i = ind_i, p = ptr_i, x = data_i,
                    dims = shape, index1 = TRUE)
  rm(data_i, ind_i, ptr_i); gc(verbose = FALSE)
  ## h5 不存条码（RCTD 只需要 names 非空做记账）；补合成名
  if (length(ctypes) != shape[2L]) stop("cell_types 长度与矩阵列数不符", call. = FALSE)
  nm <- sanitize_type_names(ctypes)
  ctypes_orig <- ctypes                    # 与注释表逐字一致（G1 对 manifest 用这个）
  ctypes <- nm$x
  colnames(m) <- sprintf("refcell_%06d", seq_len(shape[2L]))
  rownames(m) <- genes_ref
  list(counts = m, genes = genes_ref, cell_types = ctypes,
       cell_types_orig = ctypes_orig,
       shape = shape, attrs = att, name_map = nm$map)
}

## ============================================================================
## 4. 切片读取（Visium filtered：只留 Gene Expression；双探针 symbol 求和）
## ============================================================================
read_positions <- function(slide_dir) {
  f <- file.path(slide_dir, "spatial", "tissue_positions.csv")
  if (!file.exists(f)) f <- file.path(slide_dir, "spatial", "tissue_positions_list.csv")
  if (!file.exists(f)) stop("找不到 tissue_positions: ", slide_dir, call. = FALSE)
  x <- read.csv(f, header = FALSE, stringsAsFactors = FALSE)
  ## 🔴 因为第一行是文字表头，read.csv(header=FALSE) 会把**所有**列判成字符型
  ##    （实测）⇒ 去掉表头后必须显式转数值，否则 SpatialRNA 报 "coords is not numeric"
  if (is.character(x[[1L]][1L]) && tolower(x[[1L]][1L]) == "barcode")
    x <- x[-1L, , drop = FALSE]
  if (ncol(x) != 6L) stop("tissue_positions 列数不是 6: ", ncol(x), call. = FALSE)
  names(x) <- c("barcode", "in_tissue", "array_row", "array_col",
                "pxl_row_in_fullres", "pxl_col_in_fullres")
  for (k in c("in_tissue", "array_row", "array_col",
              "pxl_row_in_fullres", "pxl_col_in_fullres"))
    x[[k]] <- as.numeric(x[[k]])
  if (anyNA(x$pxl_row_in_fullres) || anyNA(x$pxl_col_in_fullres))
    stop("tissue_positions 有坐标无法转数值", call. = FALSE)
  x
}

read_slide <- function(slide) {
  d <- file.path(SUP_DIR, slide, "filtered_feature_bc_matrix")
  feats <- fread(file.path(d, "features.tsv.gz"), header = FALSE, sep = "\t",
                 stringsAsFactors = FALSE, showProgress = FALSE)
  setnames(feats, c("gene_id", "symbol", "type"))
  is_ge <- feats$type == "Gene Expression"
  n_ab  <- sum(feats$type == "Antibody Capture")
  conb <- gzfile(file.path(d, "barcodes.tsv.gz"), "r")
  bc <- readLines(conb); close(conb)
  con <- gzfile(file.path(d, "matrix.mtx.gz"), "r"); on.exit(close(con), add = TRUE)
  m <- readMM(con)                                    # 行=feature，列=barcode
  if (nrow(m) != nrow(feats)) stop(slide, ": mtx 行数与 features 不符", call. = FALSE)
  if (ncol(m) != length(bc))  stop(slide, ": mtx 列数与 barcodes 不符", call. = FALSE)
  m <- as(m, "CsparseMatrix")[is_ge, , drop = FALSE]
  syms <- feats$symbol[is_ge]
  ## 双探针 symbol（如 TBCE/HSPA14/TMSB15B）按 symbol 求和，只留一行
  dup <- unique(syms[duplicated(syms)])
  if (length(dup)) {
    us <- unique(syms)
    agg <- sparseMatrix(i = match(syms, us), j = seq_along(syms), x = 1,
                        dims = c(length(us), length(syms)))
    m <- as(agg %*% m, "CsparseMatrix")
    syms <- us
  }
  rownames(m) <- syms; colnames(m) <- bc
  list(counts = m, symbols = syms, n_antibody = n_ab,
       n_ge_rows = sum(is_ge), dual_probe = dup)
}

## ============================================================================
## 5. 读 spot_mask（N0 的唯一来源）、参考 manifest、切片清单
## ============================================================================
SM <- fread(file.path(RES_DIR, "spot_mask.tsv.gz"), showProgress = FALSE)
SM$pass <- as.logical(SM$pass)
## 🔴 按切片取 N0 必须走 base 子集，**不能**写 SM[SM$slide == slide]：
##    data.table 的 `i` 会把右边的 `slide` 优先解析成**它自己的同名列**，
##    于是比较退化成「自己 == 自己」⇒ 全为 TRUE ⇒ N0 静默变成全队列总数
##    （2026-09-25 冒烟测试实测：报 608,043 而非本张的 ~9.9k）。这是静默错，不报错。
##
## 🔴🔴 第二个坑（同一天逮到，更隐蔽）：`split()` 的两个参数**必须一起过滤**。
##    错写成 split(SM$barcode[SM$pass], SM$slide) 时，第一个是 60.8 万、第二个是 63.98 万，
##    长度不等 ⇒ R **只给一句 warning**，然后把标签整体错位：本张会被塞进下一张的条码，
##    队列末端的切片会拿到**空组**。
##    实测：本张 N0 报 10,108（真值 **9,870**），多出的 238 恰是下一张的条码。
##    ⚠️ 这次喂给 RCTD 的 spot 恰好是对的 —— 但那是「外来条码被 ∩ 本张 filtered 筛掉」撞对的，
##    **不是设计对的**；换一张切片就会静默算错或空组崩掉。
N0_BC <- split(SM$barcode[SM$pass], SM$slide[SM$pass])
## 硬门甲：分组总数必须 == 总 pass 数（错位时会不符）
if (sum(lengths(N0_BC)) != sum(SM$pass)) {
  stop(sprintf("N0_BC 分组总数 %d ≠ 总 pass 数 %d ⇒ split 对齐错位，硬停",
               sum(lengths(N0_BC)), sum(SM$pass)), call. = FALSE)
}
## 硬门乙：另用一条**完全不同的机制**（data.table 分组，非 base split）逐张算一遍，
##        两套数必须逐张相等 —— 这是防「同一处逻辑错两遍」的独立口径
N0_TRUE <- SM[SM$pass == TRUE, .N, by = slide]
setkey(N0_TRUE, slide)
N0_DIFF <- names(N0_BC)[lengths(N0_BC) != N0_TRUE[names(N0_BC), N]]
if (length(N0_DIFF)) {
  stop(sprintf("N0 逐张与独立口径不符，涉及 %d 张：%s ⇒ 硬停",
               length(N0_DIFF), paste(head(N0_DIFF, 5), collapse = ", ")), call. = FALSE)
}
N0_SLIDE_N <- setNames(as.integer(N0_TRUE[names(N0_BC), N]), names(N0_BC))
REFM <- fromJSON(REF_MANI, simplifyVector = FALSE)

SLIDES <- if (!is.null(ARGS$slides)) ARGS$slides else unlist(CAL$slide_list)
if (!is.na(ARGS$smoke)) SLIDES <- ARGS$smoke
cat(sprintf("[排期] 本次跑 %d 张切片%s\n", length(SLIDES),
            if (!is.na(ARGS$smoke)) "（冒烟：单张）" else ""))

## ============================================================================
## 6. G1：参考已对齐（细胞数 / 类型数 / 每型计数 vs manifest）
## ============================================================================
cat("[G1] 读参考 ", basename(REF_H5), " …\n", sep = "")
REF <- read_reference(REF_H5)
n_ref_cells <- REF$shape[2L]
if (n_ref_cells != as.integer(REFM$n_cells_reference)) {
  stop(sprintf("G1 越界：参考细胞数 %d ≠ manifest %s ⇒ 硬停",
               n_ref_cells, REFM$n_cells_reference), call. = FALSE)
}
## 逐型计数用**改名前的原名**与 manifest 对（manifest 记的是注释表原名）
tab_ref <- table(REF$cell_types_orig)
n_types <- length(tab_ref)
if (n_types != EXPECT_N_TYPES) {
  stop(sprintf("G1 越界：参考类型数 %d ≠ calibration 期望 %d ⇒ 硬停",
               n_types, EXPECT_N_TYPES), call. = FALSE)
}
## 逐型计数 vs manifest（对得上的那几型）
mc <- REFM$cell_type_counts
if (!is.null(mc)) {
  bad <- character(0)
  for (k in names(mc)) {
    got <- if (k %in% names(tab_ref)) as.integer(tab_ref[[k]]) else 0L
    if (got != as.integer(mc[[k]])) bad <- c(bad, sprintf("%s: h5=%d manifest=%d", k, got, mc[[k]]))
  }
  if (length(bad)) stop("G1 越界：逐型计数与 manifest 不符 ⇒ 硬停\n  ",
                        paste(bad, collapse = "\n  "), call. = FALSE)
}
cat(sprintf("[G1] 通过：%d 细胞 × %d 基因；%d 型；逐型计数与 manifest 一致\n",
            n_ref_cells, REF$shape[1L], n_types))
cat(sprintf("[G1] 参考 gene symbol 唯一？ %s\n",
            if (anyDuplicated(REF$genes)) "否（有重复！）" else "是"))
## 类型名适配（工具约束）
if (is.null(REF$name_map)) {
  cat("[名] 类型名无需适配\n")
} else {
  cat("[名] spacexr 拒绝含 \"/\" 的类型名 ⇒ 送入 RCTD 前做适配（h5 未改动）：\n")
  print(REF$name_map, row.names = FALSE)
  write.table(REF$name_map, file.path(OUT_DIR, "cell_type_name_map.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

## ============================================================================
## 7. Reference() 构建（抽样逐型前后计数都要报）
## ============================================================================
cat("[参考] 构建 Reference 对象（n_max_cells=", P$n_max_cells, "）…\n", sep = "")
ref_types <- factor(REF$cell_types)
names(ref_types) <- colnames(REF$counts)              # Reference() 要求 names 非空
pre <- as.integer(table(ref_types))
names(pre) <- names(table(ref_types))

REF_OBJ <- Reference(REF$counts, ref_types, colSums(REF$counts),
                     n_max_cells = P$n_max_cells, min_UMI = P$min_UMI,
                     require_int = P$require_int)
post_tab <- table(REF_OBJ@cell_types)
post <- as.integer(post_tab); names(post) <- names(post_tab)
sampling <- data.frame(cell_type = names(pre),
                       n_before = pre[names(pre)],
                       n_after  = as.integer(post[names(pre)]),
                       row.names = NULL)
print(sampling, row.names = FALSE)
REF_GENES    <- REF$genes          # 逐张循环还要用来判基因交集
REF_N_GENES  <- REF$shape[1L]
REF_NAME_MAP <- REF$name_map
rm(REF); gc(verbose = FALSE)

## ============================================================================
## 8. 逐张切片：G1 基因交集 → RCTD → G2/G3/G5 → 落盘
## ============================================================================
GENE_INTERSECTION_EXPECTED <- 18066L   # RCTD_PREREG §2.1 实测值（56 张应逐张一致）
summary_rows <- list()
gene_ref_absent <- NULL
gene_extra <- NULL

t_start <- Sys.time()
for (i in seq_along(SLIDES)) {
  slide <- SLIDES[i]
  cat(sprintf("\n===== [%d/%d] %s =====\n", i, length(SLIDES), slide))
  t0 <- Sys.time()

  ## ---- G3 第一层：N0 ----
  n0_bc <- N0_BC[[slide]]
  if (is.null(n0_bc)) stop("spot_mask 里没有这张切片: ", slide, call. = FALSE)
  n0 <- length(n0_bc)
  ## 硬门丙（逐张）：N0 必须同时 > 0 且等于独立口径的数。
  ##   N0 == 0 ⇒ 该张拿到空组（split 错位的典型表现），不许继续往下算。
  if (n0 == 0L) stop(slide, ": N0 == 0（split 错位或该张全不合格）⇒ 硬停", call. = FALSE)
  n0_indep <- N0_SLIDE_N[[slide]]
  if (is.na(n0_indep) || n0 != n0_indep) {
    stop(sprintf("%s: N0 = %d ≠ 独立口径 %s ⇒ 硬停", slide, n0, format(n0_indep)),
         call. = FALSE)
  }


  ## ---- 读切片 ----
  sl <- read_slide(slide)
  ## G1 基因条款（🔴 2026-09-25 实测更正，见 §12）：
  ##   原文「参考基因 ⊇ 切片基因」**不可能成立**（两份基因宇宙互不包含）。
  ##   生效判据 = 交集数 == 18,066，且 56 张逐张一致。
  inter <- intersect(REF_GENES, sl$symbols)
  ref_absent  <- setdiff(REF_GENES, sl$symbols)       # 参考有、切片无
  slide_extra <- setdiff(sl$symbols, REF_GENES)       # 切片有、参考无
  if (length(inter) != GENE_INTERSECTION_EXPECTED) {
    stop(sprintf("G1 越界：%s 基因交集 %d ≠ 期望 %d ⇒ 硬停",
                 slide, length(inter), GENE_INTERSECTION_EXPECTED), call. = FALSE)
  }
  if (is.null(gene_ref_absent)) {
    gene_ref_absent <- ref_absent
    gene_extra <- slide_extra
  } else if (!setequal(gene_ref_absent, ref_absent)) {
    stop("G1 越界：切片间「参考有但切片无」的基因集不一致 ⇒ 硬停", call. = FALSE)
  }

  ## 先把 puck 对齐到 N0 ∩ filtered，再取交集基因（按参考基因序）
  keep_bc    <- intersect(n0_bc, colnames(sl$counts))
  n0_in_filt <- length(keep_bc)
  n0_not_filt <- n0 - n0_in_filt
  cnt <- sl$counts[inter, keep_bc, drop = FALSE]
  colnames(cnt) <- keep_bc

  ## ---- 尺度不变性实验：计数二项稀释（**仅 --thin 时启用**，见 RCTD_PREREG §13.3）----
  ## 逐 spot 独立地对**每一个 count** 以概率 p 保留 ⇒ 保持该 spot 的构成（期望不变），
  ## 只把深度按 p 缩小。用来直接检验「权重与 argmax 是否随测序深度移动」。
  ## 与 §1.3 那个合成实验**同构**（固定构成、只变深度），但发生在**真实切片**上。
  ## ⚠️ 默认关闭 ⇒ 56 张全跑的路径一字未变。
  if (!is.na(ARGS$thin)) {
    set.seed(ARGS$seed)
    n_before <- ncol(cnt); tot_before <- sum(cnt)
    cnt <- as(cnt, "CsparseMatrix")
    ## ⚠️ `rbinom()` 返回 integer，直接塞进 `@x` 会把矩阵的存储类型从 double 改成 integer，
    ##    下一步 `drop0()` 就会报 "REAL() can only be applied to a 'numeric'"。
    ##    ⇒ 必须显式转回 double（2026-09-25 实跑逮到，见 §12.5 第 6 项）。
    cnt@x <- as.numeric(rbinom(length(cnt@x), size = as.integer(cnt@x), prob = ARGS$thin))
    cnt <- drop0(cnt)
    cnt <- cnt[, colSums(cnt) > 0, drop = FALSE]   # 全零 spot 去掉（否则与坐标对不齐）
    keep_bc    <- colnames(cnt)
    n0_in_filt <- length(keep_bc)
    n0_not_filt <- n0 - n0_in_filt
    cat(sprintf("[%s] 二项稀释 p=%.2f：spot %d → %d，计数总量 %d → %d（保留 %.1f%%）\n",
                slide, ARGS$thin, n_before, ncol(cnt),
                as.integer(tot_before), as.integer(sum(cnt)),
                100 * sum(cnt) / tot_before))
  }

  ## ---- 坐标 ----
  pos <- read_positions(file.path(SUP_DIR, slide))
  rownames(pos) <- pos$barcode
  miss <- setdiff(keep_bc, rownames(pos))
  if (length(miss)) stop(slide, ": ", length(miss), " 个 spot 缺坐标", call. = FALSE)
  coords <- data.frame(x = pos[keep_bc, "pxl_col_in_fullres"],
                       y = pos[keep_bc, "pxl_row_in_fullres"],
                       row.names = keep_bc)

  nUMI <- colSums(cnt); names(nUMI) <- keep_bc
  puck <- SpatialRNA(coords, cnt, nUMI)

  ## ---- RCTD ----
  obj <- create.RCTD(puck, REF_OBJ,
                     max_cores = P$max_cores, test_mode = P$test_mode,
                     gene_cutoff = P$gene_cutoff, fc_cutoff = P$fc_cutoff,
                     gene_cutoff_reg = P$gene_cutoff_reg,
                     fc_cutoff_reg = P$fc_cutoff_reg,
                     UMI_min = P$UMI_min, UMI_max = P$UMI_max,
                     counts_MIN = P$counts_MIN, UMI_min_sigma = P$UMI_min_sigma,
                     CELL_MIN_INSTANCE = P$CELL_MIN_INSTANCE,
                     MAX_MULTI_TYPES = P$MAX_MULTI_TYPES,
                     keep_reference = P$keep_reference,
                     CONFIDENCE_THRESHOLD = P$CONFIDENCE_THRESHOLD,
                     DOUBLET_THRESHOLD = P$DOUBLET_THRESHOLD)
  obj <- run.RCTD(obj, doublet_mode = P$doublet_mode)

  w <- as.matrix(obj@results$weights)                 # spot × 类型

  ## ---- G2：权重非负 ----
  wmin <- min(w)
  if (!is.finite(wmin) || wmin < 0)
    stop(sprintf("G2 越界：%s 出现负权重 min=%.6g ⇒ 硬停", slide, wmin), call. = FALSE)

  ## ---- G3：三层显式计数 ----
  n1 <- ncol(obj@originalSpatialRNA@counts)
  n2 <- ncol(obj@spatialRNA@counts)
  if (nrow(w) != n2)
    stop(sprintf("G3 越界：%s nrow(weights)=%d ≠ N2=%d ⇒ 硬停", slide, nrow(w), n2),
         call. = FALSE)

  ## ---- G5：**已由硬门改为「只报告」**（2026-09-25 晚用户签字，诊断与裁定见 RCTD_PREREG §13）
  ## 原判据 |Spearman(输入总UMI, 权重行和)| <= 0.5 **不再是闸门**。理由：
  ##   ① rho(输入总UMI, n_genes) = +0.999 ⇒ 这个 rho 量的是「mRNA 总量 vs 分到的细胞质量」；
  ##   ② 深度在空间上高度成团（邻居相关 +0.827）、行和也成团（+0.629）⇒ 两者都是**组织密度**
  ##      结构，不是技术伪影。§1.3 的合成实验是**固定构成、只变深度**；真实切片里深度变化
  ##      恰恰**因为**密度变化 ⇒ 该判据无法区分「技术漂移」与「密处细胞多」；
  ##   ③ 0.5 落在实测区间内（P1_AAH 0.437 / P1_LUAD 0.521），且 rho 随 x 离散度上升
  ##      （同一张内 P25–P75 只 0.348，P10–P90 0.457）⇒ **跨切片不可比**，0.5 是在掷硬币；
  ##   ④ 它守的是 §1.3 支撑的**绝对阈值**门控，而该口径已废、改为 `argmax == 上皮`。
  ##      ⚠️ 这里**不能**再说「argmax 对权重整体缩放免疫 ⇒ 守卫对象不存在」——那是**恒等式**，
  ##      说的是把**产出权重**乘常数，与**输入深度**扰动是两回事。真检验（§13.5：计数减半重跑）
  ##      证明 argmax **会**动（移动 9.93%、门控判据翻转 8.30%，抖动集中在并列附近）。
  ##      ⇒ 撤 G5 = 撤掉一个**测不准的仪表**，**不是**它担心的事不存在（§13.5.1、§13.7）。
  ## ⇒ 逐张**记录** rho 与「是否越过 0.5 参考线」，供审计与后续诊断；不再据此停止。
  ## spot 集合 = N2；原始总 UMI 取 RCTD 的**输入**（originalSpatialRNA@counts），不取内部量
  bc_n2 <- rownames(w)
  umi_in <- colSums(obj@originalSpatialRNA@counts)[bc_n2]
  rho <- suppressWarnings(cor(as.numeric(umi_in), as.numeric(rowSums(w)),
                             method = "spearman"))
  ## 非有限值仍硬停 —— 那不是阈值问题，是算坏了
  if (!is.finite(rho)) {
    stop(sprintf("%s: G5 rho 非有限值 ⇒ 硬停（不是阈值问题）", slide), call. = FALSE)
  }
  g5_over_ref <- abs(rho) > 0.5
  if (g5_over_ref) {
    cat(sprintf("[%s] 注：G5 参考线 0.5 被越过（|rho| = %.3f）—— **不硬停**，见 RCTD_PREREG §13\n",
                slide, abs(rho)))
  }

  ## ---- 落盘：逐 spot 权重 ----
  out_f <- file.path(PER_DIR, paste0(slide, ".weights.tsv.gz"))
  con_o <- gzfile(out_f, "w")
  write.table(w, con_o, sep = "\t", quote = FALSE, col.names = NA)
  close(con_o)

  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] 完成 %.1f 分：N0=%d N1=%d N2=%d  基因交集=%d  G5 rho=%.3f\n",
              slide, el / 60, n0, n1, n2, length(inter), rho))
  cat(sprintf("[%s] 权重逐型 min/median/max:\n", slide))
  print(round(rbind(min = apply(w, 2, min), median = apply(w, 2, median),
                    max = apply(w, 2, max)), 3))

  summary_rows[[slide]] <- list(
    slide = slide, n0_qc_pass = n0, n0_in_filtered = n0_in_filt,
    n0_not_in_filtered = n0_not_filt, n1_entered = n1, n2_weights = n2,
    n0_minus_n1 = n0 - n1, n1_minus_n2 = n1 - n2,
    n_antibody_rows = sl$n_antibody, n_ge_rows = sl$n_ge_rows,
    dual_probe_symbols = I(sl$dual_probe),
    gene_intersection = length(inter),
    weight_min = wmin, weight_rowsum_min = min(rowSums(w)),
    weight_rowsum_max = max(rowSums(w)),
    g5_rho = rho,
    g5_over_ref_0.5 = g5_over_ref,
    weight_min_by_type = as.list(apply(w, 2, min)),
    weight_median_by_type = as.list(apply(w, 2, median)),
    weight_max_by_type = as.list(apply(w, 2, max)),
    elapsed_sec = el,
    weights_file = sub(paste0("^", ROOT, "/"), "", out_f))
  rm(obj, puck, w, cnt); gc(verbose = FALSE)
}

## ============================================================================
## 9. 运行 manifest（含所有参数、守卫、输入指纹）
## ============================================================================
sha256 <- function(p) {
  if (!file.exists(p)) return(NA_character_)
  o <- suppressWarnings(system2("sha256sum", shQuote(p), stdout = TRUE, stderr = FALSE))
  if (length(o) != 1L) NA_character_ else sub(" .*$", "", o)
}

manifest <- list(
  script = "08_spatial_deconv/01_run_rctd.R",
  ran_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  caliber = CALIBER,
  reference = list(file = sub(paste0("^", ROOT, "/"), "", REF_H5),
                   sha256 = sha256(REF_H5),
                   manifest = basename(REF_MANI),
                   n_cells = n_ref_cells, n_genes = REF_N_GENES,
                   n_types = n_types, cell_type_counts = as.list(tab_ref)),
  reference_sampling = sampling,
  cell_type_name_map = if (is.null(REF_NAME_MAP)) NULL else REF_NAME_MAP,
  caliber_signed = list(signed_by = CAL$signed_by, signed = CAL$signed),
  params = P,
  seed = ARGS$seed,
  guards = list(
    G1 = list(note = "参考细胞数/类型数/逐型计数 == manifest；基因交集 == 18066 且逐张一致",
              gene_intersection_expected = GENE_INTERSECTION_EXPECTED,
              reference_genes_absent_from_slides = I(gene_ref_absent),
              slide_genes_absent_from_reference = I(gene_extra)),
    G2 = "min(weights) >= 0，逐张判，越界硬停",
    G3 = "N0→N1→N2 三层分别计数，逐张上报；差值分别报",
    ## 注：该名也必须拼接 —— 上面 §1 的 G4 自检是 grep **本文件**找那个字面量，
    ##     这里若写全名会让自己硬停（与 §1 第 76 行的处理同理）。
    G4 = paste0("源码自检：不出现 ", "行归一化函数", " 调用（本文件不含该名）"),
    G5 = paste("**已由硬门改为只报告**（RCTD_PREREG §13，2026-09-25 用户签字）：",
               "逐张记录 rho = Spearman(输入原始总UMI, 权重行和) 与是否越过 0.5 参考线；",
               "越线不停止。原判据被证为混淆（深度=RCTD输入的同位组织密度，非技术伪影），",
               "且其支撑的绝对阈值门控已被 argmax==上皮 取代。",
               "⚠️ 撤回旧说法「后者对权重整体缩放免疫」——那是恒等式，与输入深度扰动无关；",
               "§13.5 真检验证明 argmax 会动（移动 9.93%、门控翻转 8.30%），",
               "撤 G5 是撤掉测不准的仪表，不是它担心的事不存在。"),
    thin = if (is.na(ARGS$thin)) "未启用" else sprintf("计数二项稀释 p = %.3f（尺度不变性实验专用）", ARGS$thin)
  ),
  slides = SLIDES,
  per_slide = summary_rows,
  elapsed_total_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
  machine = list(R = R.version.string, spacexr = as.character(packageVersion("spacexr")),
                 rhdf5 = as.character(packageVersion("rhdf5")),
                 Matrix = as.character(packageVersion("Matrix")),
                 cores_detected = parallel::detectCores(),
                 loadavg = paste(round(read.table("/proc/loadavg")[1, 1:3], 2), collapse = "/")),
  boundary = paste("本产物是**解卷积权重**，不是恶性判定；不做门控。",
                   "恶性×权重 的乘法仍被禁止（量纲不通）。")
)
write_json(manifest, file.path(OUT_DIR, "run_manifest.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = NA)

## 逐张汇总表（方便人看）
flat <- rbindlist(lapply(summary_rows, function(r)
  data.table(slide = r$slide, N0 = r$n0_qc_pass, N0_in_filt = r$n0_in_filtered,
             N1 = r$n1_entered, N2 = r$n2_weights,
             N0_N1 = r$n0_minus_n1, N1_N2 = r$n1_minus_n2,
             gene_int = r$gene_intersection, g5_rho = round(r$g5_rho, 3),
             g5_over_0.5 = r$g5_over_ref_0.5,
             w_min = round(r$weight_min, 4),
             rowsum_min = round(r$weight_rowsum_min, 4),
             mins = round(r$elapsed_sec / 60, 1))))
fwrite(flat, file.path(OUT_DIR, "summary_per_slide.tsv"), sep = "\t")

cat(sprintf("\n===== 全部完成：%d 张，合计 %.1f 分钟 =====\n",
            length(SLIDES), manifest$elapsed_total_sec / 60))
cat("产物：", OUT_DIR, "\n", sep = "")
print(flat)
