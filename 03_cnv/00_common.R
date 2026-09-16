# 03_cnv/00_common.R —— CopyKAT 分支共用工具
#
# 数据源唯一权威：results/02_expression/gse308103_counts.h5ad（GP0 产物）。
# 本文件只读该 h5ad；所有写出一律落在 results/03_cnv/ 下（/home/eto/luad_invasion 为只读源）。

suppressMessages({library(rhdf5); library(Matrix); library(jsonlite); library(digest)})

H5AD <- "/home/eto/luad_v2/results/02_expression/gse308103_counts.h5ad"

# ---- 来源断言：本 h5ad 必须是 GSE308103（snRNA），绝不可与空转混淆 -----------
# 为什么需要这条：GSE308103(snRNA) 与 GSE307534(Visium) 有 **48 个同名 sample_id**
# （如 P24_LUAD 两边都叫这个名），单用 sample_id 做 key 会把空转混进本分支。
# 此处不靠"恰好"：
#   (a) h5ad 的 sha256 必须等于 GP0 审计记录值 → 确保是**已审计的同一个字节对象**；
#   (b) manifest 的 dataset 必须是 GSE308103；
#   (c) h5ad 的 sample_id 全集必须与 GSE308103 GEO 权威表**精确相等**（双向差集为空）。
H5AD_SHA256_GP0 <- "f9dbe38231b41efed5d628902155679924d7f4d13a770e3b5ecd79f4acbca168"
GP0_MANIFEST    <- "/home/eto/luad_v2/results/02_expression/build_manifest.json"
GEO_DIR         <- "/home/eto/luad_v2/00_ingest/geo_metadata"

geo_sample_ids <- function(dataset) {
  p <- file.path(GEO_DIR, sprintf("%s_samples.tsv", dataset))
  if (!file.exists(p)) stop("GEO 权威表缺失: ", p)
  t <- read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("patient_id", "token") %in% colnames(t)))
    stop("GEO 表列名不符（需 patient_id/token）: ", p)
  sort(unique(paste0(t$patient_id, "_", t$token)))
}

.assert_provenance <- local({
  checked <- FALSE
  function(h5 = H5AD) {
    if (checked) return(invisible(TRUE))
    m <- jsonlite::fromJSON(GP0_MANIFEST)
    if (!identical(m$dataset, "GSE308103"))
      stop("GP0 manifest 的 dataset 不是 GSE308103，实为: ", m$dataset)
    if (!identical(m$artifacts$h5ad$sha256, H5AD_SHA256_GP0))
      stop("GP0 manifest 内部矛盾：记录哈希与代码内钉住值不符")
    got <- digest::digest(file = h5, algo = "sha256")
    if (!identical(got, H5AD_SHA256_GP0))
      stop(sprintf("h5ad 字节与 GP0 审计对象不符，拒绝在其上跑 CNV\n  期望 %s\n  实测 %s",
                   H5AD_SHA256_GP0, got))
    cats <- sort(unique(as.character(h5read(h5, "obs/sample_id/categories"))))
    g308 <- geo_sample_ids("GSE308103")
    g534 <- geo_sample_ids("GSE307534")
    if (!identical(cats, g308))
      stop(sprintf("h5ad 的 sample_id 集与 GSE308103 权威表不等 —— 疑有跨数据集混入\n  h5ad 独有: %s\n  表独有: %s",
                   paste(setdiff(cats, g308), collapse = ","),
                   paste(setdiff(g308, cats), collapse = ",")))
    n_collide <- length(intersect(cats, g534))
    checked <<- TRUE
    cat(sprintf("[provenance] GSE308103 snRNA 已验证: %d 样本, sha256=%.12s…；与空转(GSE307534)同名 %d 个（已按集合相等排除混入）\n",
                length(cats), H5AD_SHA256_GP0, n_collide))
    invisible(TRUE)
  }
})

# ---- 读取单个样本的 细胞×基因 CSR -------------------------------------------
# GP0 审计 A7a' 已证明每个样本的行块在 obs 中连续，故可用超切片只读该块的 nnz 区间，
# 避免把 7.9 亿个非零值整体读入。
read_h5ad_sample <- function(sample_id, h5 = H5AD) {
  if (!file.exists(h5)) stop("h5ad 不存在: ", h5)
  .assert_provenance(h5)                      # 来源硬断言，失败即停
  cats <- as.character(h5read(h5, "obs/sample_id/categories"))
  si <- match(sample_id, cats)
  if (is.na(si)) stop("样本不在 h5ad 的 sample_id 类别中: ", sample_id)

  codes <- as.integer(h5read(h5, "obs/sample_id/codes"))
  rows0 <- which(codes == si - 1L) - 1L                      # 转 0-based
  n <- length(rows0)
  if (n == 0L) stop("样本无细胞: ", sample_id)
  if ((max(rows0) - min(rows0) + 1L) != n)
    stop("样本行块在 obs 中不连续 —— 与 GP0 审计 A7a' 的前提不符，硬停")
  start0 <- min(rows0)

  ip <- as.numeric(h5read(h5, "X/indptr", start = start0 + 1L, count = n + 1L))
  nnz0 <- ip[1L]; nnz <- ip[n + 1L] - nnz0
  if (nnz <= 0) stop("样本 nnz 为 0: ", sample_id)

  x  <- as.numeric(h5read(h5,   "X/data",    start = nnz0 + 1L, count = nnz))
  j  <- as.integer(h5read(h5,   "X/indices", start = nnz0 + 1L, count = nnz))

  # ⚠️ 字符串数据集**必须全量读**，不得用 start/count 超切片。
  # 实测（rhdf5 2.42.0）：对变长字符串数据集带非零 start 的 h5read 会**段错误**
  # —— start=1 通过、start=2 通过、start=100 与 start=361252 均 SIGSEGV(139)；
  # 而数值数据集（X/data、X/indices）的超切片 6/6 稳定，含 85e6 count 与大偏移。
  # 该崩溃随样本起始行号而变（P13_Normal 起始 1 侥幸通过，P24_LUAD 起始 361252 必崩），
  # 曾伪装成"随机内存损坏"。
  genes   <- as.character(h5read(h5, "var/_index"))
  bc_all  <- as.character(h5read(h5, "obs/_index"))
  bc      <- bc_all[start0 + seq_len(n)]

  # dgRMatrix 的 (j, p, x) 正是 CSR 布局本身 → O(1) 构造，无需排序
  M <- new("dgRMatrix", j = j, p = as.integer(ip - nnz0), x = x,
           Dim = c(n, length(genes)), Dimnames = list(bc, genes))
  # 读取路径的完整性自检（数值超切片若出问题，以下断言会先炸而不是静默给错数）
  if (nnz != length(M@x)) stop("nnz 自检失败")
  if (length(j) && (min(j) < 0L || max(j) >= length(genes)))
    stop("列索引越界 —— 数值超切片读取不可信")
  if (any(diff(as.integer(ip - nnz0)) < 0L)) stop("行指针非单调 —— CSR 结构损坏")
  if (!all(grepl("\\|", bc))) stop("barcode 不含 '|' 分隔符 —— 与 GP0 V4a 前提不符")
  if (anyDuplicated(bc)) stop("barcode 重复")
  M
}

# ---- 与 M1 权威逐细胞表交叉核对 ---------------------------------------------
# 这是对**整条 R 读取路径**的独立证伪：ncCount/nFeature 来自 M1 的独立产物，
# 与 h5ad 的 CSR 数值不是同一来源。行和与非零数两项同时精确相等，几乎不可能由错位凑出。
M1_QC <- "/home/eto/luad_v2/results/01_qc/gse308103_per_cell_qc.csv.gz"

crosscheck_against_m1 <- function(M, sample_id, qc_path = M1_QC) {
  qc <- read.csv(qc_path, stringsAsFactors = FALSE)
  qc <- qc[qc$sample_id == sample_id, ]
  m <- match(rownames(M), qc$cell_barcode)
  if (anyNA(m)) stop("有 barcode 在 M1 表中找不到: ", sample_id)
  qc <- qc[m, ]
  if (!all(qc$qc_pass)) stop("读入的细胞里有 M1 判为 qc_pass=FALSE 的")
  if (!all(qc$doublet_class == "singlet")) stop("读入的细胞里有 M1 判为双体的")
  list(
    n = nrow(M),
    max_abs_dev_nCount   = max(abs(Matrix::rowSums(M) - qc$nCount)),
    max_abs_dev_nFeature = max(abs(Matrix::rowSums(M > 0) - qc$nFeature))
  )
}

# ---- 覆盖度地板（2026-09-16 用户签字采纳 840；见 GP1_report §8.3）------------
# 判据定义在 **M1 的 `nFeature`** 上 —— M1 的逐细胞表是权威产物，不是 h5ad 的派生物。
# 地板作用在**进 copykat 之前**，因此它会同时改变 `n_cells_used`（分母）与可能改变
# `n_after_LOWDR_fullgenes`（进而翻转 :57 的 7000 覆写）。**两个后果都必须重算，不得沿用旧值。**
COVERAGE_FLOOR_NFEATURE <- 840L

.m1_qc_cache <- local({
  d <- NULL
  function(p = M1_QC) {
    if (is.null(d)) d <<- read.csv(p, stringsAsFactors = FALSE)
    d
  }
})

# 返回与 rownames(M) 对齐的 M1 `nFeature`，并**当场证伪**两个来源的一致性。
m1_nfeature <- function(M, sample_id, qc_path = M1_QC) {
  qc <- .m1_qc_cache(qc_path)
  qc <- qc[qc$sample_id == sample_id, ]
  m <- match(rownames(M), qc$cell_barcode)
  if (anyNA(m)) stop("有 barcode 在 M1 表中找不到: ", sample_id)
  nf  <- qc$nFeature[m]
  dev <- max(abs(Matrix::rowSums(M > 0) - nf))
  if (dev != 0)
    stop(sprintf("M1 nFeature 与 h5ad 实测非零数不符（最大偏差 %d）: %s", dev, sample_id))
  nf
}

apply_coverage_floor <- function(M, nfeature, floor = COVERAGE_FLOOR_NFEATURE) {
  if (floor <= 0L) return(M)
  M[nfeature >= floor, , drop = FALSE]
}

# ===========================================================================
# 上皮判定 —— **只用来决定"哪些细胞该进 CNV 工具"**
#
# ⚠️ 边界（务必先读）：本函数的输出**只做子集化**，不构成任何恶性判断，也绝不作
#    CNV 结果的佐证或替代。marker 判"是不是上皮"，CNV 判"这个上皮是不是恶性"——
#    两个不同的问题，不循环。法则 2 表原话："恶性身份以 CNV 为准；marker 仅作一致性佐证"。
#
# 为什么非做不可：CopyKAT 的设计前提是输入**同一谱系**，细胞间差异主要来自 CNV。
#    喂入混合谱系时，层次聚类的第一刀切的是**细胞类型**而非染色体。实测（2026-09-16，
#    见 results/03_cnv/EPI_PREREG.md）：P11_Normal 里被判"非整倍体"的细胞
#    EPCAM=0.67/PTPRC=0.00，被判"二倍体"的 LYZ=0.63/DCN=1.02 —— 判决完全由谱系决定。
#    因此混合输入的"非整倍体比例"实际是在量**上皮细胞占比**，与染色体无关。
#
# 两条规则**同跑**（执行前登记，不得事后择一），用于检验结论对定义是否稳健：
#    argmax : 6 个谱系 score 取严格最大者；与最高分差 <= min_score_margin 的归 ambiguous
#    strict : argmax==上皮 且 {EPCAM,KRT7,KRT19} 至少一个检出 且
#             {PTPRC,CD3D,PECAM1,COL1A1} 全部零检出
# 守卫（可证伪）：上皮集合的 EPCAM 检出率必须 >= 2× 集合外；否则 stop。
# marker 来源：docs/scientific_rigor_and_audit.md 法则 2 表（逐行照抄，不增删）。
# ===========================================================================
EPI_MARKERS <- list(
  epithelial  = c("EPCAM", "KRT7", "KRT19", "NKX2-1", "SFTPC", "SFTPB", "AGER"),
  tnk         = c("CD3D", "CD3E", "CD4", "CD8A", "GZMB", "PRF1", "NKG7"),
  bplasma     = c("CD19", "MS4A1", "CD79A", "SDC1", "IGHG1"),
  myeloid     = c("CD68", "CD163", "C1QA", "C1QB", "C1QC", "SPP1", "MARCO"),
  fibroblast  = c("COL1A1", "COL1A2", "ACTA2", "PDGFRB", "FAP", "CXCL12"),
  endothelial = c("PECAM1", "VWF", "CDH5", "EGFL7", "KDR")
)

epi_call <- function(M, rule = c("argmax", "strict"), min_score_margin = 1e-9) {
  rule <- match.arg(rule)
  present <- lapply(EPI_MARKERS, function(g) intersect(g, colnames(M)))
  np <- vapply(present, length, 0L)
  if (np[["epithelial"]] == 0L)
    stop("[FAIL] 法则2 上皮 marker 在矩阵里一个都不存在，无法判定")

  # 只取 marker 列，避免为打分而展开整张稠密矩阵（k ~ 35 列）
  have <- unique(unlist(present))
  lib  <- Matrix::rowSums(M); lib[lib == 0] <- 1
  Xn   <- log1p(as.matrix(M[, have, drop = FALSE]) / lib * 1e4)   # cells × have

  score <- vapply(names(present), function(L) {
    g <- present[[L]]
    if (length(g) == 0L) return(rep(-Inf, nrow(Xn)))
    rowMeans(Xn[, g, drop = FALSE])
  }, numeric(nrow(Xn)))
  sc <- score[, "epithelial"]
  mx <- apply(score, 1L, max)
  is_argmax <- (mx - sc) <= min_score_margin              # 并列 = 判据无分辨力

  if (rule == "argmax") {
    keep <- is_argmax
  } else {
    pos <- intersect(c("EPCAM", "KRT7", "KRT19"),          colnames(Xn))
    neg <- intersect(c("PTPRC", "CD3D", "PECAM1", "COL1A1"), colnames(Xn))
    has_pos  <- if (length(pos)) rowSums(Xn[, pos, drop = FALSE] > 0) > 0  else rep(TRUE,  nrow(Xn))
    neg_zero <- if (length(neg)) rowSums(Xn[, neg, drop = FALSE] > 0) == 0 else rep(FALSE, nrow(Xn))
    keep <- is_argmax & has_pos & neg_zero
  }
  names(keep) <- rownames(M)

  det <- function(g, sel) if (g %in% colnames(Xn)) mean(Xn[sel, g] > 0) else NA_real_
  e_in <- det("EPCAM", keep); e_out <- det("EPCAM", !keep)
  if (!is.na(e_in) && !is.na(e_out) && e_out > 0 && e_in < 2 * e_out)
    stop(sprintf("[FAIL] 上皮判定守卫失败：集合内 EPCAM 检出率 %.3f < 2×集合外 %.3f", e_in, e_out))
  if (sum(keep) == 0L)
    stop("[FAIL] 上皮判定后 0 个细胞 —— 不得继续")

  missing <- lapply(names(EPI_MARKERS), function(L) setdiff(EPI_MARKERS[[L]], present[[L]]))
  names(missing) <- names(EPI_MARKERS)
  list(rule = rule, keep = keep, n_input = nrow(M), n_keep = sum(keep),
       n_drop = sum(!keep), marker_present = np, marker_missing = missing,
       epcam_detect_in = e_in, epcam_detect_out = e_out)
}

# ===========================================================================
# 单一亚型收窄（AT2）—— 见 results/03_cnv/EPI_PREREG.md 附篇（执行前登记）
#
# 为什么需要第二层：上皮子集生效后（免疫/基质 marker 全零），copykat 转而把
# AT2 与 AT1 劈开（实测 P11_Normal strict：判"非整倍体"的 SFTPC=4.79、判"二倍体"的 AGER=2.88）。
# 即"表达轴换了，猜测照旧"。本函数检验把轴也拿掉之后会怎样。
#
# marker 出处：Travaglini et al. 2020 Nature 587:619（Human Lung Cell Atlas），
# 与 M3 标准 B 的 celltypist Human_Lung_Atlas 同源。
# ⚠️ 尚未逐字回原文核对（已在 EPI_PREREG.md 附篇 A2 标注），因此**只用于本归因臂**，
# 不得进入任何正式产物。
# ===========================================================================
AT2_MARKERS <- c("SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "NAPSA", "LAMP3")
AT1_MARKERS <- c("AGER", "CAV1", "PDPN", "HOPX")

subtype_at2 <- function(M) {
  a2 <- intersect(AT2_MARKERS, colnames(M))
  a1 <- intersect(AT1_MARKERS, colnames(M))
  if (length(a2) == 0L || length(a1) == 0L)
    stop(sprintf("[FAIL] AT2/AT1 marker 不足（AT2 存活 %d、AT1 存活 %d），无法判定",
                 length(a2), length(a1)))
  lib <- Matrix::rowSums(M); lib[lib == 0] <- 1
  have <- c(a2, a1)
  Xn <- log1p(as.matrix(M[, have, drop = FALSE]) / lib * 1e4)
  s2 <- rowMeans(Xn[, a2, drop = FALSE])
  s1 <- rowMeans(Xn[, a1, drop = FALSE])
  core <- intersect(c("SFTPC", "SFTPA1", "SFTPB"), colnames(Xn))
  detected <- if (length(core)) rowSums(Xn[, core, drop = FALSE] > 0) > 0 else rep(TRUE, nrow(Xn))
  keep <- (s2 > s1) & detected
  names(keep) <- rownames(M)
  if (sum(keep) == 0L) stop("[FAIL] AT2 收窄后 0 个细胞 —— 不得继续")
  list(keep = keep, n_input = nrow(M), n_keep = sum(keep), n_drop = sum(!keep),
       at2_present = a2, at1_present = a1,
       at2_missing = setdiff(AT2_MARKERS, a2), at1_missing = setdiff(AT1_MARKERS, a1),
       sftpc_detect_in = if ("SFTPC" %in% colnames(Xn)) mean(Xn[keep, "SFTPC"] > 0) else NA_real_,
       sftpc_detect_out = if ("SFTPC" %in% colnames(Xn)) mean(Xn[!keep, "SFTPC"] > 0) else NA_real_)
}

# ---- hg20 注释表的静默丢弃清单 ---------------------------------------------
# annotateGenes.hg20() 内部一句 mat <- mat[rownames(mat) %in% shar, ] 会把不在
# full.anno 里的基因**静默丢掉**。这里显式算出来上报，绝不默不作声。
hg20_gene_fate <- function(genes) {
  fa <- get("full.anno", envir = as.environment("package:copykat"))
  kept <- intersect(genes, fa$hgnc_symbol)
  list(kept = kept, dropped = setdiff(genes, fa$hgnc_symbol),
       n_anno_symbols = length(unique(fa$hgnc_symbol)))
}

# ===========================================================================
# copykat 分段链路 —— **唯一权威实现**（02/05 共用；曾各写一份，已合并）
#
# 为什么必须合并：漏掉链路里的任何一步都会让"进入分段的基因/细胞集"算错，而算错的方向
# 往往还是**系统性**的（少删细胞 → 恶性比例虚高）。已实测踩过的坑：
#   · 漏掉 LOW.DR 基因过滤（在注释之前！）→ 判据基因集从 ~6.7k 涨到 16k，少删 17 个细胞
#   · 漏掉 cycle/HLA 删除 → 多出约 1,180 个基因（约占链路基因的 17%），少删细胞
#   · 7000 判据在**全部基因**上求值，早于注释 → 用注释子集算会多报（曾报 27，真值 23）
# 行号 = tools/copykat/R/copykat.R（本项目 vendored 源码）实测行号。
# ===========================================================================

# 注释 + abspos 排序 + 删周期/HLA（copykat :65 / :70 / :74-81）
# ⚠️ :74-81 删的是 cyclegenes[[1]] 与 HLA-*，**静默**；本项目实测命中约 1,180 个基因。
hg20_anno_genes <- function(genes) {
  fa  <- get("full.anno",  envir = as.environment("package:copykat"))
  cyc <- as.vector(get("cyclegenes", envir = as.environment("package:copykat"))[[1]])
  keep <- intersect(genes, fa$hgnc_symbol)
  fa <- fa[match(keep, fa$hgnc_symbol), c("hgnc_symbol", "chromosome_name", "abspos")]
  fa <- fa[order(as.numeric(fa$abspos)), ]                     # :70 按基因组位置重排
  fa <- fa[!(fa$hgnc_symbol %in% c(cyc, grep("^HLA-", fa$hgnc_symbol, value = TRUE))), ]
  fa$chromosome_name <- as.character(fa$chromosome_name)
  fa
}

# 覆盖度判据（copykat :86-105 ToRemov2 与 :194-213 ToRemov3 的同构两段）
# 返回逐细胞的"应删除"逻辑向量。三条判据：
#   ① length(as.numeric(cbind(chromosome_name, count))) < 5
#      ⚠️ 字面语义：cbind 出来是**两列**矩阵，as.numeric 拉直后长度 = 2 × 检出基因数，
#      故字面"< 5"实际是"检出基因数 <= 2"。此处照抄**语义**而非字面（照字面会多删细胞）。
#   ② 检出基因覆盖到的染色体数 < 全部染色体数
#      ⚠️ 分母取**全集**染色体数（copykat 用 anno.mat$chromosome_name 而非 anno.mat2 的）；
#      取子集会因基因集未覆盖全部染色体而放宽判据、少删细胞（实测少删 17 个）。
#   ③ 任一染色体上的连续检出基因段长度 < ngene.chr
#   实现：基因已按 abspos 排序 → 同一染色体的基因**连续** → rle(chrom[hit]) 的段数 =
#   被检出到的染色体数、段长 = 该染色体上的检出基因数。于是判据②③可向量化为
#   一次稀疏矩阵乘（细胞 × 染色体 的检出计数矩阵），无需逐细胞 R 循环。
fails_criterion <- function(A, gene_idx, chrom, n_chr_total, ngene_chr) {
  Sub <- A[, gene_idx, drop = FALSE]
  nz  <- Matrix::rowSums(Sub > 0)
  bad <- nz <= 2L                                             # 判据①
  ok  <- !bad
  if (any(ok)) {
    S1 <- Sub[ok, , drop = FALSE]
    S1@x <- rep(1, length(S1@x))                              # 二值化（检出/未检出）
    Indic <- Matrix::sparse.model.matrix(~ 0 + factor(chrom, levels = unique(chrom)))
    Cnt   <- as.matrix(S1 %*% Indic)                          # 细胞 × 染色体 检出基因数
    nchr  <- rowSums(Cnt > 0)                                 # 判据②
    C2 <- Cnt; C2[C2 == 0L] <- NA_integer_
    mn <- apply(C2, 1L, min, na.rm = TRUE)                     # 判据③（只看被检出到的染色体）
    bad[ok] <- nchr < n_chr_total | mn < ngene_chr
  }
  bad
}

# 完整链路：细胞过滤 → 基因过滤 → 7000 覆写 → 注释/排序/删周期HLA → ToRemov2
#           → DR2 基因过滤 → ToRemov3
# 返回 genes_final（= copykat 的 anno.mat2，即**实际进入分段的基因集**）与 survivors。
copykat_chain <- function(A, low_dr = 0.05, up_dr = 0.10,
                          ngene_chr = 5L, min_gene_per_cell = 200L) {
  # :19 源码只在**失败细胞数 > 1** 时才过滤（`if(sum(genes.raw<min.gene.per.cell)>1)`）。
  # 恰好 1 个不合格时**不过滤**，该细胞照常进入。逐字复刻：否则链条会少算 1 个细胞，
  # 令 `n_judged == n_cells_used` 守卫误炸（2026-09-16 `P11_Normal` 锚定臂实测踩到）。
  ng_cell <- Matrix::rowSums(A > 0)
  keep_cell <- if (sum(ng_cell < min_gene_per_cell) > 1L) ng_cell >= min_gene_per_cell
               else rep(TRUE, nrow(A))
  A <- A[keep_cell, , drop = FALSE]
  cells_kept <- rownames(A); n_cells <- nrow(A)

  der <- Matrix::colSums(A > 0) / n_cells                         # :52
  A   <- A[, der > low_dr, drop = FALSE]                          # :53
  n_after_lowdr <- ncol(A)
  # :57 **无条件覆写**：判据用的是上一步的全基因数（不是注释后的）
  up_eff <- if (n_after_lowdr < 7000L) low_dr else up_dr

  anno <- hg20_anno_genes(colnames(A))                            # :65/:70/:74-81
  gi   <- match(anno$hgnc_symbol, colnames(A)); stopifnot(!anyNA(gi))
  n_chr_total <- length(unique(anno$chromosome_name))
  rem2 <- fails_criterion(A, gi, anno$chromosome_name, n_chr_total, ngene_chr)  # :86-105
  A2   <- A[!rem2, , drop = FALSE]

  # :191 分母是 ToRemov2 **之后**的细胞数（不是最初的）
  dr2  <- Matrix::colSums(A2[, gi, drop = FALSE] > 0) / max(nrow(A2), 1L)
  sel2 <- dr2 >= up_eff
  # ⚠️ chrom 必须与 gene_idx 同步子集化：二者按位置对齐，只子集一个会让 chrom[hit]
  # 取到别的染色体，判据③随即把所有细胞判失败（曾实测 ToRemov3=全部存活细胞）。
  rem3 <- fails_criterion(A2, gi[sel2], anno$chromosome_name[sel2], n_chr_total, ngene_chr)

  list(cells_kept = cells_kept, n_cells_used = n_cells,
       n_after_lowdr = n_after_lowdr, up_dr_effective = up_eff,
       anno = anno, dr2 = dr2, sel2 = sel2,
       genes_final = anno$hgnc_symbol[sel2], n_genes_final = sum(sel2),
       n_chromosomes = n_chr_total,
       survivors = cells_kept[!rem2][!rem3], n_rem2 = sum(rem2), n_rem3 = sum(rem3))
}

# ---- 进程峰值内存（/usr/bin/time 不可用，改读 /proc）------------------------
peak_rss_kb <- function() {
  s <- readLines("/proc/self/status", warn = FALSE)
  as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", grep("^VmHWM:", s, value = TRUE)))
}

# ---- 极简 JSON 写出（不引额外依赖）-----------------------------------------
# ⚠️ 2026-09-16 重写。旧版只能处理**标量**：`fmt` 对长度 n 的原子向量返回长度 n，
# 而 `vapply(x, fmt, "")` 要求长度 1 → 报
#   "values must be length 1, but FUN(X[[32]]) result is length 9"
# （route A 的曲线表恰有 9 行，把整份 summary 写崩）。
# 旧版还把 NA 写成 `NA`（非法 JSON 的裸 token）；现统一写 `null`。
json_write <- function(x, path) {
  esc  <- function(v) gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", v))
  fmt1 <- function(v) {                       # 标量 → JSON 片段，长度**恒为 1**
    if (length(v) != 1L) stop("fmt1 只接受标量")
    if (is.na(v))        return("null")
    if (is.logical(v))   return(if (v) "true" else "false")
    if (is.numeric(v))   return(format(v, digits = 15, scientific = FALSE))
    paste0('"', esc(as.character(v)), '"')
  }
  fmt <- function(v) {                        # 任意向量/列表 → JSON 片段，长度**恒为 1**
    if (is.null(v))      return("null")
    if (is.list(v))      return(paste0("[", paste(vapply(v, fmt, ""), collapse = ","), "]"))
    if (length(v) == 0L) return("[]")
    if (length(v) == 1L) return(fmt1(v))
    paste0("[", paste(vapply(seq_along(v), function(i) fmt1(v[[i]]), ""), collapse = ","), "]")
  }
  body <- paste(sprintf('  "%s": %s', names(x), vapply(x, fmt, "")), collapse = ",\n")
  writeLines(paste0("{\n", body, "\n}"), path)
}
