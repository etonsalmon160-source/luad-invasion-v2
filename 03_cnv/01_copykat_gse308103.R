#!/usr/bin/env Rscript
# 03_cnv/01_copykat_gse308103.R <sample_id>
#
# GP2 · 全量 75 样本的 CopyKAT 逐样本 CNV 调用。**一次只跑一个样本**
# （并发由 03_cnv/run_gp2_queue.sh 控制，R 内部不并行：n.cores 是结果相关参数，见 §4.3）。
#
# 与 GP1 冒烟脚本（03_smoke_test_copykat.R）的三处**实质差别**：
#   ① 施加 840 覆盖度地板（共享实现 apply_coverage_floor()，不另写一份）
#   ② 抽完预测即删巨大的中间产物（raw_results_*_by_cell.txt / heatmap jpeg）
#   ③ JSON 里把三档判定分开报，绝不合并
#
# ⚠️ 关于 `not.defined` 的一个**源码事实**（必须写进报告，否则会被误读）：
#    tools/copykat/R/copykat.R:495
#       ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% names(com.preN))]
#    `not.defined` 的集合 = **传进 copykat 的矩阵**里没拿到判定的细胞。因此
#    **被 840 地板剔掉的细胞根本不在 prediction.txt 里，它们不是 not.defined，是"缺席"**。
#    下游若按 prediction.txt 的行数算"恶性比例"，分母自然是地板后的细胞数 —— 这是对的；
#    但**绝不能**把缺席的细胞当成 diploid 补回去。
#
# 用法：Rscript 03_cnv/01_copykat_gse308103.R P13_Normal
#       LUAD_GP2_FORCE=1 Rscript ...          # 忽略哨兵，重跑
# 环境：LUAD_COVERAGE_FLOOR（默认 840）

suppressMessages({library(copykat); library(Matrix)})
source("/home/eto/luad_v2/03_cnv/00_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("用法: Rscript 01_copykat_gse308103.R <sample_id>")
sample_id <- args[1]

OUT_ROOT <- Sys.getenv("LUAD_GP2_OUTROOT", "/home/eto/luad_v2/results/03_cnv/gp2")
# 对撞用的 tiers 表：默认与地板配套。改 `LUAD_GP2_TIERS` 可指向另一地板口径的表
# （仅用于**敏感性分析**这一对照臂；主口径恒为 floor840）。
TIERS    <- Sys.getenv("LUAD_GP2_TIERS",
                       "/home/eto/luad_v2/results/03_cnv/prereg_gene_tiers_floor840.csv")
FORCE    <- identical(Sys.getenv("LUAD_GP2_FORCE", "0"), "1")

# 预注册参数（与 docs/PARAMETERS_AND_SOURCES.md M2 首行逐字一致；不因样本大小调整）
# ⚠️ `LUAD_GP2_LOW_DR`/`LUAD_GP2_UP_DR` 仅供**隔离实验**（把 UP.DR 的影响从细胞地板的影响里
# 拆出来）。主口径恒为 0.05/0.10；改过阈值时与 tiers 表的对撞会自动跳过（那张表只在 0.05/0.10 下有效）。
LOW_DR <- as.numeric(Sys.getenv("LUAD_GP2_LOW_DR", "0.05"))
UP_DR  <- as.numeric(Sys.getenv("LUAD_GP2_UP_DR",  "0.10"))
DEFAULT_THRESHOLDS <- (LOW_DR == 0.05 && UP_DR == 0.10)
NGENE_CHR <- 5L
MIN_GENE_PER_CELL <- 200L ; WIN_SIZE <- 25L ; KS_CUT <- 0.1
N_CORES <- 1L
FLOOR <- as.integer(Sys.getenv("LUAD_COVERAGE_FLOOR", "840"))

# ⚠️ `LUAD_EPI_RULE` 默认**关闭**（""）→ 主口径行为与之前逐字节一致。
# 设为 argmax / strict 时，只保留**上皮细胞**进 copykat（判定见 00_common.R epi_call()）。
# 依据：CopyKAT 要求输入同一谱系；混合输入时它切的是细胞类型而非染色体（2026-09-16 实测）。
# 这是**输入对象**的修正，不是参数调整；地板与阈值不变，用于隔离"谱系混杂"这一个变量。
EPI_RULE <- Sys.getenv("LUAD_EPI_RULE", "")
# 第二层收窄：`LUAD_EPI_SUBTYPE=at2` 时，在上皮子集**之上**再只留 AT2（见 EPI_PREREG.md 附篇）。
# 仅在 EPI_RULE 非空时有意义；单独设它是配置错误，直接停。
EPI_SUBTYPE <- Sys.getenv("LUAD_EPI_SUBTYPE", "")

# ⚠️ 传统锚定臂：`LUAD_NORM_REF=nonepi` 时，**不动输入细胞集**，改为把"正常细胞"
# 显式喂给 copykat 的 `norm.cell.names`（copykat.R:136 `else if(length(norm.cell.names)>1)`
# → 打印 "baseline is from known input"），走的是源码里**不猜基线**的那条分支。
# 参考细胞 = 非上皮细胞（免疫/基质/内皮），即 PLAN_AND_CHECKPOINTS.md:89
# "参考用同样本免疫/基质细胞" 的原文要求；也是 copykat 作者推荐的用法。
# 与 LUAD_EPI_RULE 互斥：一个改参考，一个改输入，不能同时上，否则归因混淆。
NORM_REF <- Sys.getenv("LUAD_NORM_REF", "")

# ⚠️ 「锚定 + 同质」臂（2026-09-22；预注册 results/03_cnv/ANCHOR_PREREG.md）。
# 与上面两个臂都不同：它**同时**改输入（只留一个亚型）与参考（同患者 Normal 的该亚型），
# 且输入里**同时含** Normal 与病灶两个样本 —— 故必须是**独立分支**，不能复用上面任何一个。
# ⚠️ 不解除 146-147 行的 NORM_REF×EPI_RULE 互斥守卫（它对归因臂仍然正确）；本臂改为
# 与新变量互斥。
ANCHOR_MODE       <- Sys.getenv("LUAD_ANCHOR_MODE", "")
ANCHOR_REF_SAMPLE <- Sys.getenv("LUAD_ANCHOR_REF_SAMPLE", "")
ANCHOR_SUBTYPE    <- Sys.getenv("LUAD_ANCHOR_SUBTYPE", "")

is_anchor <- nzchar(ANCHOR_MODE)
if (is_anchor) {
  if (ANCHOR_MODE != "sameSubtypeNormal")
    stop(sprintf("[FAIL] 未知 LUAD_ANCHOR_MODE=%s（当前仅支持 sameSubtypeNormal）", ANCHOR_MODE))
  if (!nzchar(ANCHOR_REF_SAMPLE))
    stop("[FAIL] 锚定臂必须给出 LUAD_ANCHOR_REF_SAMPLE（锚定来源样本，如 P11_Normal）")
  if (!(ANCHOR_SUBTYPE %in% names(ANCHOR_SUBTYPES)))
    stop(sprintf("[FAIL] 锚定臂必须给出 LUAD_ANCHOR_SUBTYPE ∈ {%s}",
                 paste(names(ANCHOR_SUBTYPES), collapse = ",")))
  if (nzchar(EPI_RULE) || nzchar(EPI_SUBTYPE) || nzchar(NORM_REF))
    stop("[FAIL] 锚定臂与 LUAD_EPI_RULE / LUAD_EPI_SUBTYPE / LUAD_NORM_REF 互斥：本臂同时改输入与参考，叠加会让归因无法拆开")
  if (identical(ANCHOR_REF_SAMPLE, sample_id))
    stop("[FAIL] 锚定来源样本不能等于被跑的样本（锚定必须来自同患者的**另一个**样本）")
  if (!identical(sub("_.*$", "", ANCHOR_REF_SAMPLE), sub("_.*$", "", sample_id)))
    stop(sprintf("[FAIL] 锚定来源(%s)与被跑样本(%s)不是同一患者 —— 预注册 §2.3 禁止跨患者",
                 ANCHOR_REF_SAMPLE, sample_id))
} else if (nzchar(ANCHOR_REF_SAMPLE) || nzchar(ANCHOR_SUBTYPE)) {
  stop("[FAIL] 设了 LUAD_ANCHOR_REF_SAMPLE / LUAD_ANCHOR_SUBTYPE 却没设 LUAD_ANCHOR_MODE")
}

out_dir  <- file.path(OUT_ROOT, sample_id)
json_out <- file.path(OUT_ROOT, sprintf("%s.json", sample_id))
done_out <- file.path(OUT_ROOT, sprintf("%s.done", sample_id))   # 哨兵：只在成功后落

# ---- 0. 断点续跑：哨兵存在即视为完成（哨兵是最后一步写的，不依赖 JSON 解析）----
if (file.exists(done_out) && file.exists(json_out) && !FORCE) {
  cat(sprintf("[GP2] %s: 已完成（见 %s），跳过；LUAD_GP2_FORCE=1 可重跑\n", sample_id, done_out))
  quit(save = "no", status = 0)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()

# ---- 1. 读矩阵 + 双来源交叉核对 --------------------------------------------
M <- read_h5ad_sample(sample_id)
n_cells_h5ad <- nrow(M); n_genes_raw <- ncol(M)
rss_after_load <- peak_rss_kb()

xchk <- crosscheck_against_m1(M, sample_id)      # 与 M1 权威表对 nCount/nFeature
if (xchk$max_abs_dev_nCount != 0 || xchk$max_abs_dev_nFeature != 0)
  stop(sprintf("与 M1 交叉核对未通过: nCount偏差=%s nFeature偏差=%s",
               xchk$max_abs_dev_nCount, xchk$max_abs_dev_nFeature))

fate <- hg20_gene_fate(colnames(M))
writeLines(c(sprintf("# sample=%s  genes_raw=%d  hg20_kept=%d  dropped=%d (%.2f%%)",
                     sample_id, n_genes_raw, length(fate$kept), length(fate$dropped),
                     100 * length(fate$dropped) / n_genes_raw),
             sprintf("# full.anno 唯一 hgnc_symbol = %d", fate$n_anno_symbols),
             sprintf("# coverage_floor = %d (判据 = M1 逐细胞表的 nFeature)", FLOOR),
             fate$dropped),
           file.path(out_dir, "hg20_dropped_genes.txt"))

# ---- 2. 施加覆盖度地板（共享实现，与 02_prereg_gene_tiers.R 同一份）---------
nf <- m1_nfeature(M, sample_id)
n_cells_ge_floor <- sum(nf >= FLOOR)
Mf <- apply_coverage_floor(M, nf, FLOOR)
n_floor_dropped <- n_cells_h5ad - n_cells_ge_floor
rm(nf); invisible(gc())

if (nrow(Mf) != n_cells_ge_floor)
  stop(sprintf("[FAIL] 地板施加后行数 %d != 过地板细胞数 %d", nrow(Mf), n_cells_ge_floor))
if (nrow(Mf) == 0L)
  stop(sprintf("[FAIL] %s 过地板后 0 个细胞 —— 不得继续", sample_id))

# ---- 2a. 上皮子集（仅当 LUAD_EPI_RULE 非空；默认关闭）------------------------
# 修正的是**输入对象**：CopyKAT 要求同一谱系。子集施加在地板**之后**，
# 使"全细胞 vs 仅上皮"的对照除这一个变量外完全一致。
epi <- NULL
if (nzchar(EPI_RULE)) {
  epi <- epi_call(Mf, EPI_RULE)
  Mf  <- Mf[epi$keep, , drop = FALSE]
  stopifnot(nrow(Mf) == epi$n_keep)
  cat(sprintf("[EPI] %s rule=%s: %d → %d 上皮细胞（剔除 %d）；EPCAM检出率 内=%.3f 外=%.3f；marker存活 %s\n",
              sample_id, EPI_RULE, epi$n_input, epi$n_keep, epi$n_drop,
              epi$epcam_detect_in, epi$epcam_detect_out,
              paste(sprintf("%s=%d", names(epi$marker_present), epi$marker_present), collapse = " ")))
}

sub <- NULL
if (nzchar(EPI_SUBTYPE)) {
  if (!nzchar(EPI_RULE))
    stop("[FAIL] 设了 LUAD_EPI_SUBTYPE 却没设 LUAD_EPI_RULE —— 亚型收窄必须在上皮子集之上")
  if (EPI_SUBTYPE != "at2") stop(sprintf("[FAIL] 未知 LUAD_EPI_SUBTYPE=%s", EPI_SUBTYPE))
  sub <- subtype_at2(Mf)
  Mf  <- Mf[sub$keep, , drop = FALSE]
  stopifnot(nrow(Mf) == sub$n_keep)
  cat(sprintf("[SUB] %s subtype=at2: %d → %d AT2 细胞（剔除 %d）；SFTPC检出率 内=%.3f 外=%.3f；marker存活 AT2=%d %s / AT1=%d %s\n",
              sample_id, sub$n_input, sub$n_keep, sub$n_drop,
              sub$sftpc_detect_in, sub$sftpc_detect_out,
              length(sub$at2_present), paste(sub$at2_present, collapse = ","),
              length(sub$at1_present), paste(sub$at1_present, collapse = ",")))
}

# ---- 2c. 传统锚定：显式指定正常细胞参考（仅当 LUAD_NORM_REF 非空）------------
# 不子集化 Mf —— 输入保持不变，改的是**参考**。这是与 2a/2b 的关键区别。
norm_ref <- character(0)
nref <- NULL
if (nzchar(NORM_REF)) {
  if (NORM_REF != "nonepi")
    stop(sprintf("[FAIL] 未知 LUAD_NORM_REF=%s（当前仅支持 nonepi）", NORM_REF))
  if (nzchar(EPI_RULE))
    stop("[FAIL] LUAD_NORM_REF 与 LUAD_EPI_RULE 互斥：一个改参考、一个改输入，同时上会混淆归因")
  e <- epi_call(Mf, "argmax")                      # 复用同一套法则2 marker 判定
  norm_ref <- rownames(Mf)[!e$keep]                # 非上皮 = 免疫/基质/内皮
  # copykat 会先按 min.gene.per.cell 剔列；只有活下来的参考细胞才真正参与取中位数
  nf_ref <- Matrix::rowSums(Mf[norm_ref, , drop = FALSE] > 0)
  n_usable <- sum(nf_ref >= MIN_GENE_PER_CELL)
  # 守卫取 copykat 自己的 min.cells：baseline.norm.cl(min.cells=10)（copykat.R:167）
  if (n_usable < 10L)
    stop(sprintf("[FAIL] %s: 可用参考细胞仅 %d 个（< copykat 自身的 min.cells=10），基线不可靠",
                 sample_id, n_usable))
  nref <- list(n = length(norm_ref), n_usable = n_usable,
               epcam_detect = if ("EPCAM" %in% colnames(Mf)) mean(Mf[norm_ref, "EPCAM"] > 0) else NA_real_,
               ptprc_detect = if ("PTPRC" %in% colnames(Mf)) mean(Mf[norm_ref, "PTPRC"] > 0) else NA_real_)
  cat(sprintf("[REF] %s mode=nonepi: 参考 %d 个细胞（其中 %d 个过 min.gene.per.cell=%d）；EPCAM检出率=%.3f PTPRC检出率=%.3f\n",
              sample_id, nref$n, nref$n_usable, MIN_GENE_PER_CELL,
              nref$epcam_detect, nref$ptprc_detect))
}

# ---- 2-anchor. 「锚定 + 同质」臂的输入构造（仅当 LUAD_ANCHOR_MODE 非空）-------
# 与 2a/2b/2c 都不同：本段**同时**收窄输入（一个亚型）与显式给出参考（另一个样本），
# 并把**两个样本**拼进同一张矩阵。锚定细胞在矩阵里，故也能被判 —— 这正是 G1 阴性对照。
#
# ⚠️ 本段只在锚定臂执行；不设 LUAD_ANCHOR_MODE 时整段跳过，既有各臂逐字节不变。
# 🔴 **必须放在 2c 之后**：2c 开头有一句无条件的 `norm_ref <- character(0)`，
#    放在它之前会被清零（2026-09-22 首跑踩过：锚定传递静默失效 → G1 拿到 length 0 → 崩）。
#    锚定臂下 2a/2b/2c 全不执行，故 Mf 到这里仍是地板后的病灶矩阵，位置等价。
anchor <- NULL
if (is_anchor) {
  signed <- signed_epi_cnv_barcodes()          # 硬断言签字件 sha256

  # (1) 病灶样本（= 本脚本的 sample_id，Mf 已过地板）
  A_les <- anchor_narrow(Mf, ANCHOR_SUBTYPE, signed, sprintf("病灶 %s", sample_id))

  # (2) 锚定来源样本：独立读、独立过地板、独立交叉核对（不能用病灶侧的数顶替）
  Mr <- read_h5ad_sample(ANCHOR_REF_SAMPLE)
  if (!identical(colnames(Mr), colnames(Mf)))
    stop("[FAIL] 锚定来源样本与病灶样本的基因向量不一致，无法拼成一张矩阵")
  n_ref_h5ad <- nrow(Mr)
  xr <- crosscheck_against_m1(Mr, ANCHOR_REF_SAMPLE)
  if (xr$max_abs_dev_nCount != 0 || xr$max_abs_dev_nFeature != 0)
    stop(sprintf("[FAIL] 锚定来源样本与 M1 交叉核对未通过: nCount偏差=%s nFeature偏差=%s",
                 xr$max_abs_dev_nCount, xr$max_abs_dev_nFeature))
  nfr <- m1_nfeature(Mr, ANCHOR_REF_SAMPLE)
  n_ref_ge_floor <- sum(nfr >= FLOOR)
  Mrf <- apply_coverage_floor(Mr, nfr, FLOOR)
  rm(nfr); invisible(gc())
  if (nrow(Mrf) == 0L)
    stop(sprintf("[FAIL] 锚定来源 %s 过地板后 0 个细胞 —— 基线不可能可靠", ANCHOR_REF_SAMPLE))
  A_ref <- anchor_narrow(Mrf, ANCHOR_SUBTYPE, signed, sprintf("锚定 %s", ANCHOR_REF_SAMPLE))
  rm(Mr, Mrf); invisible(gc())

  # (3) G5：锚定细胞数必须 >= copykat 自身的 min.cells=10（与 2c 同一判据、同一常数）
  if (A_ref$n_usable < 10L)
    stop(sprintf("[FAIL] %s: 可用锚定细胞仅 %d 个（< copykat 自身的 min.cells=10），基线不可靠",
                 ANCHOR_REF_SAMPLE, A_ref$n_usable))

  # (4) 拼成一张矩阵：病灶在前、锚定在后（顺序不影响 copykat，但固定下来便于复现）
  Mf <- rbind(A_les$M, A_ref$M)
  stopifnot(nrow(Mf) == A_les$n_of_subtype + A_ref$n_of_subtype)
  norm_ref <- rownames(A_ref$M)                # ← 显式锚定，不是空串

  # (5) 复用 2c 的那套 nref 结构，使 G3（日志守卫）与 JSON 的既有字段一并生效。
  #     EPCAM/PTPRC 检出率用来核对"锚定细胞确实是 AT2 而非混入免疫/基质"。
  nref <- list(mode = ANCHOR_MODE,
               n = length(norm_ref), n_usable = A_ref$n_usable,
               epcam_detect = if ("EPCAM" %in% colnames(A_ref$M)) mean(A_ref$M[, "EPCAM"] > 0) else NA_real_,
               ptprc_detect = if ("PTPRC" %in% colnames(A_ref$M)) mean(A_ref$M[, "PTPRC"] > 0) else NA_real_)

  anchor <- list(
    mode = ANCHOR_MODE, subtype_key = ANCHOR_SUBTYPE, subtype_name = A_les$subtype_name,
    ref_sample = ANCHOR_REF_SAMPLE,
    ref_n_h5ad = n_ref_h5ad, ref_n_after_floor = n_ref_ge_floor,
    ref_n_in_signed = A_ref$n_in_signed, ref_n_of_subtype = A_ref$n_of_subtype,
    ref_n_unmapped_dropped = A_ref$n_unmapped_dropped,
    ref_n_other_subtype = A_ref$n_other_subtype, ref_n_usable = A_ref$n_usable,
    les_n_after_floor = A_les$n_after_floor, les_n_in_signed = A_les$n_in_signed,
    les_n_of_subtype = A_les$n_of_subtype, les_n_unmapped_dropped = A_les$n_unmapped_dropped,
    les_n_other_subtype = A_les$n_other_subtype, les_n_usable = A_les$n_usable,
    signed_n_total = length(signed), seed_col = ANCHOR_SEED_COL,
    epicnv_sha256 = ANCHOR_EPICNV_SHA256,
    epia_clusters_sha256 = ANCHOR_EPIA_CLUSTERS_SHA256,
    epia_annot_sha256 = ANCHOR_EPIA_ANNOT_SHA256)
  cat(sprintf("[ANCHOR] %s: 拼成 %d 细胞（病灶 %d + 锚定 %d）；锚定占 %.3f；EPCAM检出率=%.3f PTPRC检出率=%.3f\n",
              sample_id, nrow(Mf), A_les$n_of_subtype, A_ref$n_of_subtype,
              A_ref$n_of_subtype / nrow(Mf), nref$epcam_detect, nref$ptprc_detect))
  rm(signed); invisible(gc())
}

# 进入后续链路的细胞数。EPI 关闭时恒等于 n_cells_ge_floor —— 三个守卫据此统一。
n_cells_in <- nrow(Mf)

# ---- 2b. 与预先算好的 tiers 表对撞（跨产物证伪：两个独立实现必须一致）-------
# ⚠️ 该表只对"全细胞 + 地板 840 + 0.05/0.10"这一条口径有效，故改了输入或参考的臂下自动跳过。
if (file.exists(TIERS) && DEFAULT_THRESHOLDS && !nzchar(EPI_RULE) && !nzchar(NORM_REF) && !is_anchor) {
  tr <- read.csv(TIERS, stringsAsFactors = FALSE)
  tr <- tr[tr$sample_id == sample_id, ]
  if (nrow(tr) != 1L) stop(sprintf("[FAIL] tiers 表里 %s 有 %d 行", sample_id, nrow(tr)))
  if (tr$coverage_floor != FLOOR)
    stop(sprintf("[FAIL] tiers 表的地板=%d，本脚本=%d", tr$coverage_floor, FLOOR))
  if (tr$n_cells_ge_floor != n_cells_ge_floor)
    stop(sprintf("[FAIL] %s: tiers 表记 %d 个细胞过地板，本脚本实测 %d —— 两份实现不一致",
                 sample_id, tr$n_cells_ge_floor, n_cells_ge_floor))
}

# ---- 3. 复算链路，取得**预测的**有效 UP.DR 与基因集 -------------------------
# copykat_chain 只做子集化，不改 Mf 本身；之后 Mf 仍用于转稠密
ch <- copykat_chain(Mf, low_dr = LOW_DR, up_dr = UP_DR,
                    ngene_chr = NGENE_CHR, min_gene_per_cell = MIN_GENE_PER_CELL)

# 地板 >= 200 时它严格强于 min.gene.per.cell 判据 → 两者必须相等（可证伪）
if (FLOOR >= MIN_GENE_PER_CELL && ch$n_cells_used != n_cells_in)
  stop(sprintf("[FAIL] %s: 送入 %d，链路收到 %d —— 地板/上皮子集未正确施加",
               sample_id, n_cells_in, ch$n_cells_used))

# 与 tiers 表的其余列对撞（同源不同产物，仍值得逐列核）
if (exists("tr")) {
  if (tr$n_cells_used != ch$n_cells_used ||
      tr$n_after_LOWDR_fullgenes != ch$n_after_lowdr ||
      tr$n_genes_final != ch$n_genes_final ||
      abs(tr$copykat_effective_UPDR - ch$up_dr_effective) > 1e-12)
    stop(sprintf("[FAIL] %s: tiers 表与 copykat_chain 复算不一致（cells_used/%s %s, LOWDRfull/%s %s, genes_final/%s %s, UPDR/%s %s）",
                 sample_id, tr$n_cells_used, ch$n_cells_used,
                 tr$n_after_LOWDR_fullgenes, ch$n_after_lowdr,
                 tr$n_genes_final, ch$n_genes_final,
                 tr$copykat_effective_UPDR, ch$up_dr_effective))
}

# ---- 4. 转 copykat 要求的 基因×细胞 稠密矩阵 -------------------------------
# 必须用**地板后**的矩阵：copykat 的 not.defined 分母就是它（源码 :495）
mat_gc <- as.matrix(t(Mf))
rm(Mf); invisible(gc())
rss_after_dense <- peak_rss_kb()

n_pred_nd <- ch$n_cells_used - length(ch$survivors)   # 链路预测的 not.defined

cat(sprintf("[GP2] %s: h5ad=%d -floor=%d%s%s%s → %d cells × %d genes  dense=%.2f GB  rss=%.2f GB  UPDR_eff=%.2f  nd预测=%d\n",
            sample_id, n_cells_h5ad, n_floor_dropped,
            if (is.null(epi)) "" else sprintf(" -非上皮=%d", epi$n_drop),
            if (is.null(sub)) "" else sprintf(" -非AT2=%d", sub$n_drop),
            if (is.null(anchor)) "" else sprintf(" [锚定臂: 病灶%s %d + 锚定%s %d]",
                                                 anchor$subtype_name, anchor$les_n_of_subtype,
                                                 anchor$ref_sample, anchor$ref_n_of_subtype),
            n_cells_in, nrow(mat_gc),
            prod(dim(mat_gc)) * 8 / 2^30, rss_after_dense / 2^20, ch$up_dr_effective, n_pred_nd))

# ---- 5. 调用（参数全部预注册，n.cores 固定）--------------------------------
old <- setwd(out_dir)
sink_log <- file.path(out_dir, "copykat.stdout.log")
sink(sink_log, split = TRUE)
t_copy <- system.time({
  res <- tryCatch(
    copykat(rawmat = mat_gc, id.type = "S", cell.line = "no",
            ngene.chr = NGENE_CHR, min.gene.per.cell = MIN_GENE_PER_CELL,
            LOW.DR = LOW_DR, UP.DR = UP_DR, win.size = WIN_SIZE,
            KS.cut = KS_CUT, distance = "euclidean", genome = "hg20",
            n.cores = N_CORES, sam.name = sample_id,
            norm.cell.names = norm_ref,
            output.seg = "FALSE", plot.genes = "FALSE"),
    error = function(e) structure(conditionMessage(e), class = "copykat_error"))
})
sink()   # sink(<路径字符串>) 的对应关闭就是 sink()，不可 close(路径)
setwd(old)
wall_copykat <- unname(t_copy["elapsed"])
rss_final <- peak_rss_kb()

# ⚠️ 必须在 rm(res) **之前**取回错误：tryCatch 是主来源（stop() 走 stderr，
# sink(type="output") 收不到；冒烟阶段曾因吞掉 copykat 的真实报错而误判）
err <- if (inherits(res, "copykat_error")) { e <- as.character(res); class(e) <- NULL; e } else NULL
rm(mat_gc, res); invisible(gc())   # 稠密矩阵立刻还给系统，否则算进后续所有操作的峰值

# copykat 自己会打印 "low data quality"（= 触发了 :57 的覆写）；从日志取回用于对撞
logged_override <- FALSE
logged_known_normal <- FALSE
if (file.exists(sink_log)) {
  lg <- readLines(sink_log, warn = FALSE)
  logged_override <- any(grepl("low data quality", lg, fixed = TRUE))
  # copykat.R:145 只在**走了锚定分支**时打印这句。用作可证伪守卫：
  # 传了 norm.cell.names 却没打这句 = 没走锚定分支，本次结果不能当锚定臂用。
  logged_known_normal <- any(grepl("baseline is from known input", lg, fixed = TRUE))
}
if (length(norm_ref) > 0L) {
  cat(sprintf("[REF] %s copykat 日志出现 \"baseline is from known input\" = %s\n",
              sample_id, logged_known_normal))
  if (!is.null(err)) {
    cat("[REF] copykat 报错，锚定分支是否走到无法判定（err 非空）\n")
  } else if (!logged_known_normal) {
    stop(sprintf("[FAIL] %s: 传了 %d 个参考细胞，但 copykat 没打印 \"baseline is from known input\" —— 未走锚定分支，结果不可用",
                 sample_id, length(norm_ref)))
  }
}

# ---- 6. 解析预测（三档分开，**绝不合档**）----------------------------------
pred_file <- file.path(out_dir, sprintf("%s_copykat_prediction.txt", sample_id))
cnt <- list(aneuploid = NA_integer_, diploid = NA_integer_,
            not_defined = NA_integer_, other = NA_integer_)
if (is.null(err) && file.exists(pred_file)) {
  pr <- read.delim(pred_file, stringsAsFactors = FALSE)
  pcol <- grep("copykat\\.pred", colnames(pr), value = TRUE)
  if (length(pcol) == 0L)
    stop(sprintf("[FAIL] %s: prediction.txt 里没有 copykat.pred 列（实际列：%s）",
                 sample_id, paste(colnames(pr), collapse = ",")))
  lab  <- tolower(as.character(pr[[pcol[1]]]))
  # 解析成功 → 全部归零再计数：table() 不含零水平，若某档真的一个都没有，
  # 字面量里就查不到它（曾把真值 0 误写成 NA）
  cnt$aneuploid   <- sum(grepl("aneuploid", lab))
  cnt$diploid     <- sum(grepl("diploid",   lab))
  cnt$not_defined <- sum(grepl("not\\.defined", lab))
  cnt$other       <- nrow(pr) - cnt$aneuploid - cnt$diploid - cnt$not_defined
  if (cnt$other > 0)
    warning(sprintf("%s: 有 %d 行的 copykat.pred 不属于三档已知值（已单列 other，未并入任何档）",
                    sample_id, cnt$other))
}

n_judged <- sum(unlist(cnt), na.rm = TRUE)

# **可证伪守卫**：prediction.txt 行数必须恰等于**进入 copykat 的细胞数**（源码 :495：
# ndef = 传入矩阵里没拿到判定的细胞，故 判定行数 ≡ ncol(rawmat)）。每个进入的细胞
# 恰好被处理一次——要么带标签，要么 not.defined。
#
# ⚠️ 2026-09-16 修：原先比的是 `n_cells_ge_floor`（过地板的细胞数）。地板 >= 200 时
# 二者恒等，所以一直没暴露；但**地板为 0 时 copykat 自己的 min.gene.per.cell=200
# 过滤器会再剔掉几个细胞**，于是合法地小于过地板数，守卫误炸
# （实测 ab_B_P14_Normal：过地板 1412、进 copykat 1410 → 误报"地板未作用于传入矩阵"）。
# 正确的比较基准是 ch$n_cells_used，它与地板、上皮子集都无关。
if (is.null(err) && !is.na(n_judged) && n_judged != ch$n_cells_used)
  stop(sprintf("[FAIL] %s: prediction.txt 判定 %d 个细胞，但进入 copykat 的是 %d",
               sample_id, n_judged, ch$n_cells_used))

frac <- if (!is.na(n_judged) && n_judged > 0) cnt$aneuploid / n_judged else NA_real_

# ---- 6b. 锚定臂分侧计分 + G1 锚定纯度（阴性对照）---------------------------
# 预注册 §4 的可证伪守卫：**锚定细胞（已知正常）被判非整倍体的比例**。
#   > 0.50 → 基线无意义，该运行 `implausible`，停下报告；> 0.10 → 打旗标单列。
# ⚠️ 两个阈值是「本项目约定」，无文献出处（预注册 §4 已登记）。整条曲线一并上报。
# 分侧计分只对锚定臂做：它把"判决"拆成 Normal 侧与病灶侧，是本臂要看的核心量。
anchor_purity <- NULL
if (is_anchor && is.null(err) && file.exists(pred_file)) {
  cn <- as.character(pr$cell.names)
  if (!all(norm_ref %in% cn))
    stop(sprintf("[FAIL] %d 个锚定细胞在 prediction.txt 里找不到（cell.names 不匹配）",
                 sum(!(norm_ref %in% cn))))
  is_anc <- cn %in% norm_ref
  lab_a  <- lab[is_anc]                       # 锚定侧
  lab_l  <- lab[!is_anc]                      # 病灶侧
  # 两侧都必须非空，否则下面的比例是 0/0 = NaN，if(NaN) 会报一句看不懂的错
  # （2026-09-22 首跑踩过：锚定传递被清零，这里 length(lab_a)=0 → 崩在 if 上）
  if (length(lab_a) == 0L || length(lab_l) == 0L)
    stop(sprintf("[FAIL] 锚定臂分侧失败：锚定侧 %d 核、病灶侧 %d 核（两侧都必须非空）",
                 length(lab_a), length(lab_l)))
  frac_anc <- sum(grepl("aneuploid", lab_a)) / length(lab_a)
  frac_les <- sum(grepl("aneuploid", lab_l)) / length(lab_l)
  verdict <- if (frac_anc > 0.50) "stop_anchor_implausible" else
             if (frac_anc > 0.10) "flag_anchor_purity"     else "ok"
  anchor_purity <- list(
    n_anchor_judged = length(lab_a), n_lesion_judged = length(lab_l),
    n_anchor_aneuploid = sum(grepl("aneuploid", lab_a)),
    n_lesion_aneuploid = sum(grepl("aneuploid", lab_l)),
    frac_aneuploid_anchor = frac_anc, frac_aneuploid_lesion = frac_les,
    lesion_minus_anchor = frac_les - frac_anc,
    G1_thresholds = "flag>0.10 stop>0.50 (本项目约定, 无文献出处)",
    G1_verdict = verdict)
  cat(sprintf("[G1] %s: 锚定侧 %d 核 判非整倍体 %d (%.3f) | 病灶侧 %d 核 判 %d (%.3f) | 差 %+.3f ⇒ %s\n",
              sample_id, length(lab_a), anchor_purity$n_anchor_aneuploid, frac_anc,
              length(lab_l), anchor_purity$n_lesion_aneuploid, frac_les,
              frac_les - frac_anc, verdict))
}

# ---- 7. 抽完即删中间产物 ----------------------------------------------------
# 保留：prediction.txt（判定）、final_results_bin_by_cell.txt（CNA 矩阵，下游 M4 要用）
# 删除：两个 raw_results_*_by_cell.txt（P19 实测 1.5 GB + 3.2 GB）+ heatmap jpeg
keep <- c(sprintf("%s_copykat_prediction.txt", sample_id),
          sprintf("%s_copykat_final_results_bin_by_cell.txt", sample_id),
          sprintf("%s_copykat_CNA_results.txt", sample_id),
          "hg20_dropped_genes.txt", "copykat.stdout.log")
freed_mb <- 0
for (fp in list.files(out_dir, full.names = TRUE)) {
  if (basename(fp) %in% keep) next
  sz <- file.size(fp)
  if (unlink(fp, force = TRUE) == 0L) freed_mb <- freed_mb + sz / 2^20
}

# ---- 8. 写 JSON（成功或**失败都写**，失败时 error 非 null，便于事后汇总）----
epi_present <- if (is.null(epi)) NULL else
  paste(sprintf("%s=%d", names(epi$marker_present), epi$marker_present), collapse = ";")
epi_missing <- if (is.null(epi)) NULL else
  paste(sprintf("%s:%s", names(epi$marker_missing),
                vapply(epi$marker_missing, function(g) if (length(g)) paste(g, collapse = ",") else "-", "")),
        collapse = ";")
json_write(list(
  stage = "GP2",
  sample_id = sample_id,
  stage_token = sub("^[^_]*_", "", sample_id),
  r_version = as.character(getRversion()),
  copykat_version = as.character(packageVersion("copykat")),
  param_low_dr = LOW_DR, param_up_dr_passed = UP_DR,
  param_ngene_chr = NGENE_CHR, param_min_gene_per_cell = MIN_GENE_PER_CELL,
  param_win_size = WIN_SIZE, param_ks_cut = KS_CUT,
  param_distance = "euclidean", param_genome = "hg20", param_n_cores = N_CORES,
  param_cell_line = "no", param_id_type = "S",
  param_h5ad_sha256 = H5AD_SHA256_GP0,
  param_coverage_floor = FLOOR,
  # ---- 细胞数三段（地板前 / 地板后 / 进 copykat）----
  n_cells_h5ad = n_cells_h5ad,
  n_cells_ge_floor = n_cells_ge_floor,
  n_cells_dropped_by_floor = n_floor_dropped,
  n_cells_used = ch$n_cells_used,
  n_cells_dropped_below_200genes = n_cells_in - ch$n_cells_used,
  # ---- 上皮子集（关闭时全为 null）----
  epi_rule = if (is.null(epi)) NULL else epi$rule,
  epi_n_input = if (is.null(epi)) NULL else epi$n_input,
  epi_n_kept = if (is.null(epi)) NULL else epi$n_keep,
  epi_n_dropped = if (is.null(epi)) NULL else epi$n_drop,
  epi_marker_present = epi_present,
  epi_marker_missing = epi_missing,
  epi_epcam_detect_in = if (is.null(epi)) NULL else epi$epcam_detect_in,
  epi_epcam_detect_out = if (is.null(epi)) NULL else epi$epcam_detect_out,
  # ---- 单一亚型收窄（关闭时为 null）----
  subtype = if (is.null(sub)) NULL else "at2",
  subtype_n_input = if (is.null(sub)) NULL else sub$n_input,
  subtype_n_kept = if (is.null(sub)) NULL else sub$n_keep,
  subtype_n_dropped = if (is.null(sub)) NULL else sub$n_drop,
  subtype_marker_present = if (is.null(sub)) NULL else
    paste(sprintf("AT2=%s;AT1=%s", paste(sub$at2_present, collapse = ","),
                  paste(sub$at1_present, collapse = ","))),
  subtype_sftpc_detect_in = if (is.null(sub)) NULL else sub$sftpc_detect_in,
  subtype_sftpc_detect_out = if (is.null(sub)) NULL else sub$sftpc_detect_out,
  # ---- 传统锚定参考（关闭时全为 null）----
  norm_ref_mode = if (is.null(nref)) NULL else nref$mode,
  norm_ref_n = if (is.null(nref)) NULL else nref$n,
  norm_ref_n_usable = if (is.null(nref)) NULL else nref$n_usable,
  norm_ref_epcam_detect = if (is.null(nref)) NULL else nref$epcam_detect,
  norm_ref_ptprc_detect = if (is.null(nref)) NULL else nref$ptprc_detect,
  norm_ref_branch_confirmed = if (is.null(nref)) NULL else logged_known_normal,
  # ---- 「锚定 + 同质」臂（未开时为 null）----
  # ⚠️ 计数口径：n_cells_h5ad / n_cells_ge_floor / n_cells_dropped_by_floor 指**病灶样本**；
  #    锚定来源样本的对应数在 anchor_ref_* 里。合并后的进 copykat 数在 n_cells_in。
  anchor_mode = if (is.null(anchor)) NULL else anchor$mode,
  anchor_subtype_key = if (is.null(anchor)) NULL else anchor$subtype_key,
  anchor_subtype_name = if (is.null(anchor)) NULL else anchor$subtype_name,
  anchor_seed_col = if (is.null(anchor)) NULL else anchor$seed_col,
  anchor_ref_sample = if (is.null(anchor)) NULL else anchor$ref_sample,
  anchor_counts_scope_note = if (is.null(anchor)) NULL else
    "n_cells_h5ad/n_cells_ge_floor/n_cells_dropped_by_floor = 病灶样本；锚定来源见 anchor_ref_*",
  anchor_ref_n_h5ad = if (is.null(anchor)) NULL else anchor$ref_n_h5ad,
  anchor_ref_n_after_floor = if (is.null(anchor)) NULL else anchor$ref_n_after_floor,
  anchor_ref_n_in_signed = if (is.null(anchor)) NULL else anchor$ref_n_in_signed,
  anchor_ref_n_of_subtype = if (is.null(anchor)) NULL else anchor$ref_n_of_subtype,
  anchor_ref_n_other_subtype = if (is.null(anchor)) NULL else anchor$ref_n_other_subtype,
  anchor_ref_n_unmapped_dropped = if (is.null(anchor)) NULL else anchor$ref_n_unmapped_dropped,
  anchor_ref_n_usable = if (is.null(anchor)) NULL else anchor$ref_n_usable,
  anchor_les_n_after_floor = if (is.null(anchor)) NULL else anchor$les_n_after_floor,
  anchor_les_n_in_signed = if (is.null(anchor)) NULL else anchor$les_n_in_signed,
  anchor_les_n_of_subtype = if (is.null(anchor)) NULL else anchor$les_n_of_subtype,
  anchor_les_n_other_subtype = if (is.null(anchor)) NULL else anchor$les_n_other_subtype,
  anchor_les_n_unmapped_dropped = if (is.null(anchor)) NULL else anchor$les_n_unmapped_dropped,
  anchor_les_n_usable = if (is.null(anchor)) NULL else anchor$les_n_usable,
  anchor_signed_n_total = if (is.null(anchor)) NULL else anchor$signed_n_total,
  anchor_epicnv_sha256 = if (is.null(anchor)) NULL else anchor$epicnv_sha256,
  anchor_epia_clusters_sha256 = if (is.null(anchor)) NULL else anchor$epia_clusters_sha256,
  anchor_epia_annot_sha256 = if (is.null(anchor)) NULL else anchor$epia_annot_sha256,
  # ---- G1 锚定纯度（阴性对照；未开锚定臂时为 null）----
  anchor_n_judged = if (is.null(anchor_purity)) NULL else anchor_purity$n_anchor_judged,
  lesion_n_judged = if (is.null(anchor_purity)) NULL else anchor_purity$n_lesion_judged,
  anchor_n_aneuploid = if (is.null(anchor_purity)) NULL else anchor_purity$n_anchor_aneuploid,
  lesion_n_aneuploid = if (is.null(anchor_purity)) NULL else anchor_purity$n_lesion_aneuploid,
  frac_aneuploid_anchor = if (is.null(anchor_purity)) NULL else anchor_purity$frac_aneuploid_anchor,
  frac_aneuploid_lesion = if (is.null(anchor_purity)) NULL else anchor_purity$frac_aneuploid_lesion,
  lesion_minus_anchor_frac = if (is.null(anchor_purity)) NULL else anchor_purity$lesion_minus_anchor,
  G1_thresholds = if (is.null(anchor_purity)) NULL else anchor_purity$G1_thresholds,
  G1_verdict = if (is.null(anchor_purity)) NULL else anchor_purity$G1_verdict,
  n_cells_in = n_cells_in,
  # ---- 基因集四段 ----
  n_genes_raw = n_genes_raw,
  n_genes_hg20 = length(fate$kept),
  n_genes_dropped_by_anno = length(fate$dropped),
  n_after_LOWDR_fullgenes = ch$n_after_lowdr,
  n_after_anno_cycle_hla = nrow(ch$anno),
  n_genes_final = ch$n_genes_final,
  # ---- 覆写状态：预测 vs copykat 自己日志里的话 ----
  up_dr_effective_predicted = ch$up_dr_effective,
  up_dr_override_predicted = (ch$up_dr_effective == LOW_DR),
  up_dr_override_logged_by_copykat = logged_override,
  # ---- 判定（三档分开，绝不合档）----
  n_judged = n_judged,
  n_pred_aneuploid = cnt$aneuploid,
  n_pred_diploid = cnt$diploid,
  n_pred_not_defined = cnt$not_defined,
  n_pred_other_label = cnt$other,
  frac_aneuploid_of_judged = frac,
  frac_not_defined_of_judged = if (!is.na(n_judged) && n_judged > 0) cnt$not_defined / n_judged else NA_real_,
  n_pred_not_defined_expected = n_pred_nd,
  rate_not_defined_expected = n_pred_nd / ch$n_cells_used,
  # ---- 成本 ----
  dense_gb = n_cells_in * n_genes_raw * 8 / 2^30,
  wall_sec_load = unname((proc.time() - t0)["elapsed"]) - wall_copykat,
  wall_sec_copykat = wall_copykat,
  peak_rss_gb_after_load = rss_after_load / 2^20,
  peak_rss_gb_after_dense = rss_after_dense / 2^20,
  peak_rss_gb_final = rss_final / 2^20,
  disk_freed_mb = freed_mb,
  # ---- 交叉核对 ----
  m1_max_abs_dev_nCount = xchk$max_abs_dev_nCount,
  m1_max_abs_dev_nFeature = xchk$max_abs_dev_nFeature,
  error = err,
  output_dir = out_dir
), json_out)

# 哨兵：**只在完全成功**（无 error 且判定数对得上）时落；断点续跑据此跳过
if (is.null(err) && !is.na(n_judged) && n_judged == ch$n_cells_used) {
  writeLines(c(sprintf("sample_id=%s", sample_id),
               sprintf("finished_at=%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
               sprintf("n_judged=%d", n_judged),
               if (is.null(anchor_purity)) NULL else
                 sprintf("G1_verdict=%s frac_aneuploid_anchor=%.4f frac_aneuploid_lesion=%.4f",
                         anchor_purity$G1_verdict, anchor_purity$frac_aneuploid_anchor,
                         anchor_purity$frac_aneuploid_lesion),
               sprintf("json=%s", json_out)), done_out)
} else {
  cat(sprintf("[GP2] %s: **未完成**（error=%s）→ 不落哨兵，队列不会跳过它\n",
              sample_id, if (is.null(err)) "判定数不符" else err))
}

cat(sprintf("[GP2] %s: copykat=%.1fs 峰值RSS=%.2f GB  判定 %d = aneuploid %s + diploid %s + not.defined %s%s  释放 %.0f MB  err=%s\n",
            sample_id, wall_copykat, rss_final / 2^20, n_judged,
            cnt$aneuploid, cnt$diploid, cnt$not_defined,
            if (is.null(anchor_purity)) "" else
              sprintf("  [锚定侧 %.3f / 病灶侧 %.3f ⇒ %s]",
                      anchor_purity$frac_aneuploid_anchor, anchor_purity$frac_aneuploid_lesion,
                      anchor_purity$G1_verdict),
            freed_mb, if (is.null(err)) "无" else err))
if (!is.null(anchor_purity) && anchor_purity$G1_verdict != "ok")
  cat(sprintf("[G1] ⚠️ %s ⇒ %s：预注册 §5 规定该情形**停下报告**，不得据此宣称任何恶性结论\n",
              sample_id, anchor_purity$G1_verdict))
