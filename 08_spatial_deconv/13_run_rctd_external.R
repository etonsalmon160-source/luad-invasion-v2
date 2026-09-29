#!/usr/bin/env Rscript
## 13_run_rctd_external.R —— 在**外部对照数据集**上跑 RCTD（GSE248082 的 N1 / N3）
## ============================================================================
## 口径来源（全部已签字；本脚本**不引入任何新参数**）
##   SPATIAL_CNV_PREREG.md §16        外部对照锚 GSE248082 首次登记（2026-09-29 用户签字）
##   SPATIAL_CNV_PREREG.md §12        上皮门控口径 = `argmax == 上皮`（**只存在于粗版 a**）
##   RCTD_PREREG.md §2 / §3 / §10     参考、参数、粗细各一版
##   calibration.json  signed=TRUE 才准跑（未签 → 硬停）
##
## 与 01_run_rctd.R 的**全部差别**（逐条，别再多）：
##   1. 计数从 **h5** 读（GEO 包不含 matrix.mtx.gz）；坐标从**扁平**的
##      tissue_positions_list.csv.gz 读（GEO 包不含 spatial/ 子目录），**带表头**。
##   2. spot 掩码**现算**（本队列用的是冻结掩码；外部切片不在 spot_mask.tsv.gz 里）。
##      规则照 SP §3.3 已签：nUMI<500 或 nFeature<200 或 pct_mt>15% ⇒ 剔（OR）。
##   3. G1 基因交集硬门（18,066）**放宽为只报告**（§16.3）——那是本队列 vs 参考的实测值，
##      外部切片与参考基因宇宙互不包含，硬停等于因口径错配而拒跑。
##   4. **两版参考都跑**（a 与 d）；上皮 spot 只用 a 的 `argmax == 上皮`（§16.3 登记更正）。
##
## 本脚本**不做**的事（边界，别越）
##   - 不做恶性判定；权重不是恶性标签。
##   - **不**并进主队列 56 张；**不**给细版 d 定义任何上皮映射。
##   - 不调用 spacexr 的行归一化函数（§M5 / G4 禁令）：本文件**不含**该函数名。
##
## 用法
##   Rscript 08_spatial_deconv/13_run_rctd_external.R --max-cores 4
## ============================================================================

suppressMessages({
  library(Matrix)
  library(spacexr)
  library(rhdf5)
  library(jsonlite)
  library(data.table)
})

ROOT     <- "/home/eto/luad_v2"
EXT_DIR  <- file.path(ROOT, "data", "external", "GSE248082")
RES_DIR  <- file.path(ROOT, "results", "08_spatial_deconv")

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
  out <- list(max_cores = 4L, seed = 1L, calibers = c("a", "d"), out = NULL)
  i <- 1L
  while (i <= length(a)) {
    k <- a[i]
    if (k == "--max-cores")      { out$max_cores <- as.integer(a[i + 1L]); i <- i + 2L }
    else if (k == "--seed")      { out$seed <- as.integer(a[i + 1L]); i <- i + 2L }
    else if (k == "--calibers")  { out$calibers <- strsplit(a[i + 1L], ",")[[1]]; i <- i + 2L }
    else if (k == "--out")       { out$out <- a[i + 1L]; i <- i + 2L }
    else stop("未知参数: ", k)
  }
  out
}
ARGS <- parse_args(commandArgs(TRUE))
if (!identical(ARGS$calibers, c("a", "d")))
  stop("§16.3 只授权 a 与 d 两版；--calibers 不得改成别的组合", call. = FALSE)

## ============================================================================
## 1. G4 守卫：源码自检（在**任何**计算之前）
## ============================================================================
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
## 2. 口径签字校验（未签不许算）
## ============================================================================
CAL_PATH <- file.path(SCRIPT_DIR, "calibration.json")
if (!file.exists(CAL_PATH)) stop("找不到 calibration.json: ", CAL_PATH, call. = FALSE)
CAL <- fromJSON(CAL_PATH, simplifyVector = FALSE)
if (!isTRUE(CAL$signed)) stop("calibration.json 的 signed 不为 TRUE ⇒ 硬停", call. = FALSE)

OUT_DIR <- if (is.null(ARGS$out)) file.path(RES_DIR, "rctd_external_GSE248082") else ARGS$out
PER_DIR <- file.path(OUT_DIR, "per_slide")
dir.create(PER_DIR, recursive = TRUE, showWarnings = FALSE)

## ============================================================================
## 3. 外部切片清单（逐字；N2/N4/N5 为什么不在，见 §16.2）
## ============================================================================
EXT <- data.frame(
  slide  = c("N1", "N3"),
  gsm    = c("GSM8087031", "GSM8087033"),
  mean_counts_per_spot = c(4823, 2705),   # 2026-09-25 实测，只作报告
  stringsAsFactors = FALSE)
EXT$h5   <- file.path(EXT_DIR, sprintf("%s_%s_filtered_feature_bc_matrix.h5", EXT$gsm, EXT$slide))
EXT$pos  <- file.path(EXT_DIR, sprintf("%s_%s_tissue_positions_list.csv.gz", EXT$gsm, EXT$slide))
for (k in c("h5", "pos")) if (!all(file.exists(EXT[[k]])))
  stop("外部文件缺失（", k, "）：", paste(EXT[[k]][!file.exists(EXT[[k]])], collapse = ", "),
       call. = FALSE)
cat(sprintf("[输入] 外部切片 %d 张：%s（%s）\n", nrow(EXT), paste(EXT$slide, collapse = ", "),
            paste(EXT$gsm, collapse = ", ")))

## ============================================================================
## 4. 读取：h5 计数、扁平坐标、现算掩码
## ============================================================================
## —— 4.1 h5（10x v3：CSC，**细胞为列**）——
read_h5_counts <- function(h5) {
  shape <- as.integer(h5read(h5, "matrix/shape"))          # c(基因, 细胞)
  if (length(shape) != 2L) stop("matrix/shape 不是长度 2", call. = FALSE)
  ft_id   <- as.character(h5read(h5, "matrix/features/id"))
  ft_name <- as.character(h5read(h5, "matrix/features/name"))
  ft_type <- as.character(h5read(h5, "matrix/features/feature_type"))
  bc      <- as.character(h5read(h5, "matrix/barcodes"))
  data_i  <- h5read(h5, "matrix/data")
  ind_i   <- h5read(h5, "matrix/indices")
  ptr_i   <- h5read(h5, "matrix/indptr")
  m <- sparseMatrix(i = as.integer(ind_i) + 1L, p = as.integer(ptr_i),
                    x = as.double(data_i), dims = shape, index1 = TRUE)
  rm(data_i, ind_i, ptr_i); gc(verbose = FALSE)
  if (nrow(m) != length(ft_id) || ncol(m) != length(bc))
    stop("h5 矩阵形状与 features / barcodes 长度不符", call. = FALSE)
  rownames(m) <- ft_name; colnames(m) <- bc
  list(counts = as(m, "CsparseMatrix"), symbols = ft_name, types = ft_type,
       n_ge_rows = sum(ft_type == "Gene Expression"),
       n_ab_rows = sum(ft_type == "Antibody Capture"))
}

## —— 4.2 扁平坐标（**带表头**；本队列那份在 spatial/ 下且也是带表头的 v2 格式）——
read_positions_flat <- function(p) {
  x <- read.csv(gzfile(p), header = FALSE, stringsAsFactors = FALSE)
  if (is.character(x[[1L]][1L]) && tolower(x[[1L]][1L]) == "barcode")
    x <- x[-1L, , drop = FALSE]
  if (ncol(x) != 6L) stop("tissue_positions 列数不是 6: ", ncol(x), call. = FALSE)
  names(x) <- c("barcode", "in_tissue", "array_row", "array_col",
                "pxl_row_in_fullres", "pxl_col_in_fullres")
  for (k in c("in_tissue", "array_row", "array_col",
              "pxl_row_in_fullres", "pxl_col_in_fullres")) x[[k]] <- as.numeric(x[[k]])
  if (anyNA(x$pxl_row_in_fullres) || anyNA(x$pxl_col_in_fullres))
    stop("tissue_positions 有坐标无法转数值", call. = FALSE)
  x
}

## —— 4.3 现算掩码（照 SP §3.3 已签阈值；⛔ 不是本脚本自定的新阈值）——
MIN_UMI_SPOT   <- 500L
MIN_GENES_SPOT <- 200L
MAX_MT_FRAC    <- 0.15
MT_PREFIX      <- "MT-"
qc_metrics <- function(m, syms) {
  ## 只算 Gene Expression 行（本队列 45 张带抗体行，外部 0 行，但口径必须一致）
  ge <- which(!grepl("^IgG|^CD", syms) & !is.na(syms))   # 占位：由调用方传 GE 行
  stop("internal: 请用 qc_metrics_ge()")
}
qc_metrics_ge <- function(m_ge, syms_ge) {
  numi <- Matrix::colSums(m_ge)
  ngen <- Matrix::colSums(m_ge > 0)
  is_mt <- grepl(paste0("^", MT_PREFIX), syms_ge)
  mtc <- Matrix::colSums(m_ge[is_mt, , drop = FALSE])
  pct_mt <- ifelse(numi > 0, 100 * mtc / numi, NA_real_)
  pass <- !(numi < MIN_UMI_SPOT | ngen < MIN_GENES_SPOT | pct_mt > 100 * MAX_MT_FRAC)
  data.frame(barcode = colnames(m_ge), nUMI = numi, nFeature = ngen,
             pct_mt = pct_mt, pass = pass, row.names = NULL)
}

## ============================================================================
## 5. 参考读取（与 01_run_rctd.R 逐字同源：CSC 细胞为列 + "/" 名字适配）
## ============================================================================
sanitize_type_names <- function(x) {
  u <- unique(x); bad <- grepl("/", u, fixed = TRUE)
  if (!any(bad)) return(list(x = x, map = NULL))
  u2 <- gsub("/", "_", u, fixed = TRUE)
  if (anyDuplicated(u2)) stop("类型名去斜杠后出现重名 ⇒ 硬停", call. = FALSE)
  list(x = u2[match(x, u)],
       map = data.frame(original = u[bad], adapted = u2[bad], row.names = NULL))
}
read_reference <- function(h5) {
  att       <- h5readAttributes(h5, "counts")
  shape     <- as.integer(att$shape)
  genes_ref <- as.character(h5read(h5, "genes"))
  ctypes    <- as.character(h5read(h5, "cell_types"))
  m <- sparseMatrix(i = as.integer(h5read(h5, "counts/indices")) + 1L,
                    p = as.integer(h5read(h5, "counts/indptr")),
                    x = as.double(h5read(h5, "counts/data")),
                    dims = shape, index1 = TRUE)
  gc(verbose = FALSE)
  if (length(ctypes) != shape[2L]) stop("cell_types 长度与矩阵列数不符", call. = FALSE)
  nm <- sanitize_type_names(ctypes)
  colnames(m) <- sprintf("refcell_%06d", seq_len(shape[2L])); rownames(m) <- genes_ref
  list(counts = m, genes = genes_ref, cell_types = nm$x,
       cell_types_orig = ctypes, shape = shape, name_map = nm$map)
}

## ============================================================================
## 6. 主循环：逐 caliber × 逐切片
## ============================================================================
GENE_INTERSECTION_EXPECTED <- 18066L   # 只用于「说明为什么外部切片不能套这个门」
summary_rows <- list()
t_start <- Sys.time()

for (CALIBER in ARGS$calibers) {
  EXPECT_N_TYPES <- as.integer(CAL$expect_n_types[[CALIBER]])
  REF_H5   <- file.path(RES_DIR, paste0("reference_", CALIBER, ".h5"))
  REF_MANI <- file.path(RES_DIR, paste0("reference_", CALIBER, ".manifest.json"))
  if (!file.exists(REF_H5) || !file.exists(REF_MANI))
    stop("参考缺失: ", REF_H5, call. = FALSE)
  REFM <- fromJSON(REF_MANI, simplifyVector = FALSE)

  cat(sprintf("\n############ caliber = %s（期望 %d 型） ############\n",
              CALIBER, EXPECT_N_TYPES))
  REF <- read_reference(REF_H5)
  if (REF$shape[2L] != as.integer(REFM$n_cells_reference))
    stop("G1 越界：参考细胞数 ≠ manifest ⇒ 硬停", call. = FALSE)
  tab_ref <- table(REF$cell_types_orig)
  if (length(tab_ref) != EXPECT_N_TYPES)
    stop(sprintf("G1 越界：参考类型数 %d ≠ 期望 %d ⇒ 硬停", length(tab_ref), EXPECT_N_TYPES),
         call. = FALSE)
  if (!is.null(REFM$cell_type_counts)) {
    bad <- character(0)
    for (k in names(REFM$cell_type_counts)) {
      got <- if (k %in% names(tab_ref)) as.integer(tab_ref[[k]]) else 0L
      if (got != as.integer(REFM$cell_type_counts[[k]]))
        bad <- c(bad, sprintf("%s: h5=%d manifest=%d", k, got, REFM$cell_type_counts[[k]]))
    }
    if (length(bad)) stop("G1 越界：逐型计数与 manifest 不符 ⇒ 硬停\n  ",
                          paste(bad, collapse = "\n  "), call. = FALSE)
  }
  cat(sprintf("[G1] 参考通过：%d 细胞 × %d 基因；%d 型\n",
              REF$shape[2L], REF$shape[1L], length(tab_ref)))
  if (!is.null(REF$name_map)) {
    cat("[名] 类型名含 \"/\" ⇒ 送入 RCTD 前适配（h5 未改动）\n")
    print(REF$name_map, row.names = FALSE)
  }

  set.seed(ARGS$seed)
  P <- list(doublet_mode = "full", UMI_min = 100, CELL_MIN_INSTANCE = 25, counts_MIN = 10,
            max_cores = ARGS$max_cores, test_mode = FALSE,
            gene_cutoff = 0.000125, fc_cutoff = 0.5,
            gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75,
            UMI_max = 2e7, UMI_min_sigma = 300, MAX_MULTI_TYPES = 4,
            keep_reference = FALSE, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20,
            n_max_cells = 10000, min_UMI = 100, require_int = TRUE)
  if (!is.null(CAL$n_max_cells) && as.integer(CAL$n_max_cells) != P$n_max_cells)
    stop("n_max_cells 与 calibration.json 不一致 ⇒ 硬停", call. = FALSE)

  ref_types <- factor(REF$cell_types); names(ref_types) <- colnames(REF$counts)
  pre <- as.integer(table(ref_types)); names(pre) <- names(table(ref_types))
  REF_OBJ <- Reference(REF$counts, ref_types, colSums(REF$counts),
                       n_max_cells = P$n_max_cells, min_UMI = P$min_UMI,
                       require_int = P$require_int)
  post_tab <- table(REF_OBJ@cell_types)
  sampling <- data.frame(cell_type = names(pre), n_before = pre[names(pre)],
                         n_after = as.integer(post_tab[names(pre)]), row.names = NULL)
  print(sampling, row.names = FALSE)
  REF_GENES <- REF$genes
  rm(REF); gc(verbose = FALSE)

  for (i in seq_len(nrow(EXT))) {
    slide <- EXT$slide[i]
    cat(sprintf("\n===== [%s] 外部切片 %s（%s）=====\n", CALIBER, slide, EXT$gsm[i]))
    t0 <- Sys.time()

    ## ---- 读 ----
    h5c <- read_h5_counts(EXT$h5[i])
    cat(sprintf("[%s] h5：%d 行（GE %d / 抗体 %d）× %d spot\n",
                slide, length(h5c$symbols), h5c$n_ge_rows, h5c$n_ab_rows, ncol(h5c$counts)))
    if (h5c$n_ab_rows != 0L)
      stop(slide, ": 外部切片出现抗体行，与 §16.2 的核验不符 ⇒ 硬停", call. = FALSE)
    is_ge <- h5c$types == "Gene Expression"
    m_ge  <- as(h5c$counts[is_ge, , drop = FALSE], "CsparseMatrix")
    syms  <- h5c$symbols[is_ge]
    ## 双探针 symbol 求和（与 01_run_rctd.R 同法；外部探针集相同 ⇒ 预期同样命中）
    dup <- unique(syms[duplicated(syms)])
    if (length(dup)) {
      us  <- unique(syms)
      agg <- sparseMatrix(i = match(syms, us), j = seq_along(syms), x = 1,
                          dims = c(length(us), length(syms)))
      m_ge <- as(agg %*% m_ge, "CsparseMatrix"); syms <- us
    }
    rownames(m_ge) <- syms

    pos <- read_positions_flat(EXT$pos[i]); rownames(pos) <- pos$barcode
    ## ---- R4：in_tissue==1 是否 ≡ filtered barcodes（首次在外部切片上核对）----
    n_in  <- sum(pos$in_tissue == 1)
    n_fil <- ncol(m_ge)
    r4_ok <- (n_in == n_fil)
    cat(sprintf("[%s] R4：in_tissue==1 = %d，filtered barcodes = %d ⇒ %s\n",
                slide, n_in, n_fil, if (r4_ok) "一致" else "🔴 不一致（登记）"))

    ## ---- 掩码现算 ----
    qc <- qc_metrics_ge(m_ge, syms)
    rownames(qc) <- qc$barcode
    n_pass <- sum(qc$pass)
    cat(sprintf("[%s] 掩码现算（%d/%d/%g，OR）：pass %d / %d（剔 %.2f%%）\n",
                slide, MIN_UMI_SPOT, MIN_GENES_SPOT, MAX_MT_FRAC,
                n_pass, nrow(qc), 100 * (1 - n_pass / nrow(qc))))
    cat(sprintf("[%s]   单项：nUMI<%d %d；nFeature<%d %d；pct_mt>%g%% %d\n", slide,
                MIN_UMI_SPOT, sum(qc$nUMI < MIN_UMI_SPOT),
                MIN_GENES_SPOT, sum(qc$nFeature < MIN_GENES_SPOT),
                100 * MAX_MT_FRAC, sum(qc$pct_mt > 100 * MAX_MT_FRAC)))

    n0_bc <- qc$barcode[qc$pass]
    ## ---- G1 放宽：只报告交集 ----
    inter <- intersect(REF_GENES, syms)
    cat(sprintf("[%s] G1 **只报告**（硬门已按 §16.3 放宽）：交集 %d（本队列期望 %d）；参考有切片无 %d；切片有参考无 %d\n",
                slide, length(inter), GENE_INTERSECTION_EXPECTED,
                length(setdiff(REF_GENES, syms)), length(setdiff(syms, REF_GENES))))

    keep_bc <- intersect(n0_bc, colnames(m_ge))
    if (!length(keep_bc)) stop(slide, ": 掩码后无可用 spot ⇒ 硬停", call. = FALSE)
    cnt <- m_ge[inter, keep_bc, drop = FALSE]
    miss <- setdiff(keep_bc, rownames(pos))
    if (length(miss)) stop(slide, ": ", length(miss), " 个 spot 缺坐标 ⇒ 硬停", call. = FALSE)
    coords <- data.frame(x = pos[keep_bc, "pxl_col_in_fullres"],
                         y = pos[keep_bc, "pxl_row_in_fullres"], row.names = keep_bc)
    nUMI <- colSums(cnt); names(nUMI) <- keep_bc
    puck <- SpatialRNA(coords, cnt, nUMI)

    obj <- create.RCTD(puck, REF_OBJ,
                       max_cores = P$max_cores, test_mode = P$test_mode,
                       gene_cutoff = P$gene_cutoff, fc_cutoff = P$fc_cutoff,
                       gene_cutoff_reg = P$gene_cutoff_reg, fc_cutoff_reg = P$fc_cutoff_reg,
                       UMI_min = P$UMI_min, UMI_max = P$UMI_max, counts_MIN = P$counts_MIN,
                       UMI_min_sigma = P$UMI_min_sigma,
                       CELL_MIN_INSTANCE = P$CELL_MIN_INSTANCE,
                       MAX_MULTI_TYPES = P$MAX_MULTI_TYPES,
                       keep_reference = P$keep_reference,
                       CONFIDENCE_THRESHOLD = P$CONFIDENCE_THRESHOLD,
                       DOUBLET_THRESHOLD = P$DOUBLET_THRESHOLD)
    obj <- run.RCTD(obj, doublet_mode = P$doublet_mode)
    w <- as.matrix(obj@results$weights)
    wmin <- min(w)
    if (!is.finite(wmin) || wmin < 0)
      stop(sprintf("G2 越界：%s 出现负权重 min=%.6g ⇒ 硬停", slide, wmin), call. = FALSE)
    n1 <- ncol(obj@originalSpatialRNA@counts); n2 <- ncol(obj@spatialRNA@counts)
    if (nrow(w) != n2) stop("G3 越界：nrow(weights) ≠ N2 ⇒ 硬停", call. = FALSE)

    ## ---- 上皮 spot（**只在粗版 a**；细版不给上皮映射，§16.3）----
    n_epi <- NA_integer_; epi_barcodes <- character(0)
    if (CALIBER == "a") {
      if (!("上皮" %in% colnames(w)))
        stop("粗版 a 的权重里没有「上皮」列 ⇒ 与 §12 口径不符，硬停", call. = FALSE)
      W  <- w
      mx <- do.call(pmax, as.data.frame(W))
      n_at <- rowSums(W >= mx)
      am <- colnames(W)[max.col(W, ties.method = "first")]
      tie <- n_at > 1L
      epi_barcodes <- rownames(W)[am == "上皮" & !tie]
      n_epi <- length(epi_barcodes)
      cat(sprintf("[%s] 上皮 spot = %d / %d（%.1f%%）  · 并列剔除 %d\n",
                  slide, n_epi, nrow(W), 100 * n_epi / nrow(W), sum(tie)))
      con_e <- gzfile(file.path(PER_DIR, sprintf("%s.%s.epi_barcodes.txt.gz", slide, CALIBER)), "w")
      writeLines(epi_barcodes, con_e); close(con_e)
    } else {
      cat(sprintf("[%s] 细版 d：**不**产出上皮名单（§16.3：给 d 定义上皮 = 新口径，未签）\n", slide))
    }

    out_f <- file.path(PER_DIR, sprintf("%s.%s.weights.tsv.gz", slide, CALIBER))
    con_o <- gzfile(out_f, "w"); write.table(w, con_o, sep = "\t", quote = FALSE, col.names = NA)
    close(con_o)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("[%s|%s] 完成 %.1f 分：N0=%d N1=%d N2=%d；上皮=%s\n", slide, CALIBER,
                el / 60, length(n0_bc), n1, n2,
                if (is.na(n_epi)) "不产出" else as.character(n_epi)))

    summary_rows[[paste(CALIBER, slide)]] <- list(
      caliber = CALIBER, slide = slide, gsm = EXT$gsm[i],
      n_h5_barcodes = n_fil, n_in_tissue_1 = n_in, r4_in_tissue_eq_filtered = r4_ok,
      n_mask_pass = n_pass, pct_dropped_by_mask = round(100 * (1 - n_pass / nrow(qc)), 3),
      n_umi_lt500 = sum(qc$nUMI < MIN_UMI_SPOT),
      n_gene_lt200 = sum(qc$nFeature < MIN_GENES_SPOT),
      n_mt_gt15 = sum(qc$pct_mt > 100 * MAX_MT_FRAC),
      n0_entered_rctd = n1, n2_weights = n2,
      gene_intersection = length(inter),
      gene_intersection_expected_cohort = GENE_INTERSECTION_EXPECTED,
      n_ref_genes_absent = length(setdiff(REF_GENES, syms)),
      n_slide_genes_absent = length(setdiff(syms, REF_GENES)),
      dual_probe_symbols = I(dup),
      n_epi_argmax = if (is.na(n_epi)) NA_integer_ else n_epi,
      frac_epi_argmax = if (is.na(n_epi)) NA_real_ else round(n_epi / n2, 4),
      weight_min = wmin,
      weight_rowsum_min = min(rowSums(w)), weight_rowsum_max = max(rowSums(w)),
      elapsed_sec = el,
      weights_file = sub(paste0("^", ROOT, "/"), "", out_f))
    rm(obj, puck, w, cnt, m_ge); gc(verbose = FALSE)
  }
  rm(REF_OBJ); gc(verbose = FALSE)
}

## ============================================================================
## 7. 运行 manifest
## ============================================================================
sha256 <- function(p) {
  if (!file.exists(p)) return(NA_character_)
  o <- suppressWarnings(system2("sha256sum", shQuote(p), stdout = TRUE, stderr = FALSE))
  if (length(o) != 1L) NA_character_ else sub(" .*$", "", o)
}
manifest <- list(
  script = "08_spatial_deconv/13_run_rctd_external.R",
  ran_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  purpose = "SPATIAL_CNV_PREREG.md §16：外部对照锚 GSE248082（N1/N3）首次跑 RCTD；**不产生恶性判定**",
  caliber_signed = list(signed_by = CAL$signed_by, signed = CAL$signed),
  calibers_run = ARGS$calibers,
  register_correction = paste("§16.3 登记更正：细版 d 的 39 型里**没有**「上皮」类型名，",
                              "已签的上皮门控 `argmax == 上皮` 只存在于粗版 a ⇒ 两版都跑，",
                              "上皮名单**只**由 a 产出；d 的权重仅留档，本次不给 d 定义上皮映射。"),
  external = list(dataset = "GSE248082", slides = EXT$slide, gsm = EXT$gsm,
                  source_dir = sub(paste0("^", ROOT, "/"), "", EXT_DIR),
                  platform = "10x Visium Cytassist FFPE",
                  probe_set_same_as_cohort = TRUE,
                  gene_order_expected_rows = 18085L,
                  excluded = list(N2 = paste("平均 668 counts/spot，极浅；",
                                             "且深度正是让 RCTD 门控抖 9.93% 的那条轴（§16.2）"),
                                  N4 = "10x Visium（3′ 冷冻）36,601 基因，基因空间不同",
                                  N5 = "同 N4")),
  params_identical_to_01_run_rctd = list(
    doublet_mode = "full", UMI_min = 100, counts_MIN = 10, CELL_MIN_INSTANCE = 25,
    UMI_max = 2e7, UMI_min_sigma = 300, MAX_MULTI_TYPES = 4,
    CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20,
    gene_cutoff = 0.000125, fc_cutoff = 0.5,
    gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75,
    n_max_cells = 10000, min_UMI = 100, require_int = TRUE,
    max_cores = ARGS$max_cores, seed = ARGS$seed),
  new_params_this_script_only = list(
    counts_source = "h5（matrix/{data,indices,indptr,shape,features,barcodes}）",
    positions_source = "扁平 *_tissue_positions_list.csv.gz（带表头）",
    mask_rule = sprintf("现算：nUMI<%d 或 nFeature<%d 或 pct_mt>%g%% ⇒ 剔（OR）；阈值照 SP §3.3 已签",
                        MIN_UMI_SPOT, MIN_GENES_SPOT, 100 * MAX_MT_FRAC),
    mask_rule_note = "🔴 本队列用的是冻结掩码；外部切片不在 spot_mask.tsv.gz 里 ⇒ 同规则现算，非同一条产线",
    g1_hard_gate = paste("🔴 **放宽为只报告**：硬门 18066 是本队列 vs 参考的实测值，",
                         "外部切片与参考基因宇宙互不包含 ⇒ 硬停等于因口径错配拒跑（§16.3）")),
  inputs_sha256 = setNames(lapply(c(EXT$h5, EXT$pos), sha256), basename(c(EXT$h5, EXT$pos))),
  references = lapply(ARGS$calibers, function(cc) list(
    caliber = cc, expect_n_types = as.integer(CAL$expect_n_types[[cc]]),
    file = basename(file.path(RES_DIR, paste0("reference_", cc, ".h5"))),
    sha256 = sha256(file.path(RES_DIR, paste0("reference_", cc, ".h5"))))),
  per_slide = summary_rows,
  elapsed_total_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
  machine = list(R = R.version.string, spacexr = as.character(packageVersion("spacexr")),
                 cores_detected = parallel::detectCores(),
                 loadavg = paste(round(read.table("/proc/loadavg")[1, 1:3], 2), collapse = "/")),
  not_claimed = c("不构成恶性判定；权重不是恶性标签",
                  "不并进主队列 56 张；不授权把外部锚用于全队列排期（借锚口径仍未签，§14.6）",
                  "尚未做「外部锚自身 CNV 是否平」的体检；不平则该锚作废",
                  "不声称外部切片的谱系注释可信度与本队列相同（跨实验室/跨患者/批次效应未评估）",
                  "n=2（N1/N3）很薄，只能作辅助锚",
                  "细版 d 的权重本次**没有**上皮映射，不得读成「更精细的锚」"),
  boundary = "本产物是解卷积权重，不是恶性判定；不做门控以外的任何判定。"
)
write_json(manifest, file.path(OUT_DIR, "run_manifest.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = NA)

flat <- rbindlist(lapply(summary_rows, function(r) data.table(
  caliber = r$caliber, slide = r$slide, N0_mask = r$n_mask_pass,
  N1 = r$n0_entered_rctd, N2 = r$n2_weights,
  gene_int = r$gene_intersection, epi_argmax = r$n_epi_argmax,
  frac_epi = r$frac_epi_argmax, w_min = round(r$weight_min, 4),
  mins = round(r$elapsed_sec / 60, 1))))
fwrite(flat, file.path(OUT_DIR, "summary_per_slide.tsv"), sep = "\t")

cat(sprintf("\n===== 全部完成：%d 次运行，合计 %.1f 分钟 =====\n",
            length(summary_rows), manifest$elapsed_total_sec / 60))
cat("产物：", OUT_DIR, "\n", sep = "")
print(flat, row.names = FALSE)
