#!/usr/bin/env Rscript
# ============================================================================
# 10_seurat_traditional.R —— M3-A 传统分支主干（**逐字对齐源论文**）
#
# 权威出处（见 docs/PARAMETERS_AND_SOURCES.md §M3-A.0 / §M3-A.1–A.3）：
#   Peng F, Sinjab A, Dai Y, … Wang L, Kadara H.
#   "Multimodal spatial-omics reveal co-evolution of alveolar progenitors and
#    proinflammatory niches in progression of lung precursor lesions."
#   Cancer Cell 2026;44(2):321-339.e13 · DOI 10.1016/j.ccell.2025.10.004
#
#   论文配方：SCTransform(v2) → HVG 3000 → PCA 50 → Harmony(默认参数, 50 PC 空间)
#             → FindNeighbors(SNN) → FindClusters(resolution 0.5–0.8) → RunUMAP
#
# 用法：
#   Rscript 04_integration/10_seurat_traditional.R smoke --n_cells 20000
#   Rscript 04_integration/10_seurat_traditional.R full
#
# 输入：results/04_integration/seurat_io/
#         （02_expression/04_rebuild_expression_paperqc.py 产出，413,697 核 × 18,069 基因）
# 输出：results/04_integration/seurat_trad/<tag>/
#         clusters.csv.gz        簇标签（Harmony 臂 4 分辨率 × 5 种子；不校正臂 4 分辨率 × 种子0）
#         resolution_metrics.csv 指标1/2 表（含 r* 候选）
#         （指标4 · AAH 吸收护栏**不是单独文件**，而是 resolution_metrics.csv 的
#           aah_absorption_rate 列 —— 旧版头注释误写为 aah_absorption.csv，2026-09-17 更正）
#         umap.csv.gz            UMAP 坐标（Harmony 空间）
#         embeddings/{harmony,pca}_f32.bin + 元信息   嵌入（供下游复用，避免重跑 SCTransform）
#         run_manifest.json      参数逐条登记 + 版本 + 哈希 + 计时
#         timings.csv            逐步耗时/内存
#
# ⚠️ 本脚本**不做隐式默认**：论文未指定的参数一律显式写出并在 manifest 里登记其来源。
# ⚠️ 硬断言：细胞数/基因数与 seurat_io manifest 不符即 abort（不"尽力而为"）。
# ⚠️ 两个臂：主 = Harmony（论文配方）；敏感性 = 不校正（**不属于论文配方**）。
# ============================================================================

suppressMessages({
  library(Seurat)
  library(Matrix)
  library(harmony)
  library(jsonlite)
})
# digest 只被 `digest::digest()` 显式调用（见下方 manifest 的 cells_file_sha256）。
# 此处提前断言，避免"跑完 20 分钟才在最后一行报缺包"。
if (!requireNamespace("digest", quietly = TRUE))
  stop("缺 R 包 digest（用于 sha256）。install.packages('digest') 后重跑。")

# ---- 0. 常量（显式登记，禁止散落）------------------------------------------
ROOT     <- "/home/eto/luad_v2"
IO       <- file.path(ROOT, "results/04_integration/seurat_io")
OUTBASE  <- file.path(ROOT, "results/04_integration/seurat_trad")

VST_FLAVOR        <- "v2"          # 论文用 Seurat 5.1 的 SCTransform（默认 v2）
N_HVG             <- 3000L         # 论文："Top 3,000 HVGs"
N_PCS             <- 50L           # 论文："The top 50 PCs"
N_CELLS_MODEL     <- 5000L         # Seurat 默认（仅用于拟合 θ～μ）
RV_TH             <- 1.3           # Seurat 默认
SEED_SCT          <- 1448145L      # Seurat 默认
SEED_PCA          <- 42L           # Seurat 默认
SEED_UMAP         <- 42L           # Seurat 默认
K_PARAM           <- 20L           # Seurat 默认
PRUNE_SNN         <- 1 / 15        # Seurat 默认
RESOLUTIONS       <- c(0.5, 0.6, 0.7, 0.8)   # 论文区间 0.5–0.8，步长 0.1 穷举
SEEDS             <- c(0L, 1L, 2L, 3L, 4L)   # 含 Seurat 默认 0，共 5 个
SEED_SUBSAMPLE    <- 20260916L     # 本项目约定：smoke 抽样种子（跑前固定）
STAGES            <- c("Normal", "AAH", "AIS", "MIA", "IAC")

# 指标1 护栏阈值（跨种子稳定性，硬约束）——见 PARAMETERS §M3-A.3
ARI_SEED_MIN <- 0.90

# 指标4 护栏阈值（本项目约定，跑前固定）——见 PARAMETERS §M3-A.3
AAH_ABSORPTION_MAX <- 0.50

T0 <- Sys.time()
log <- function(...) cat(sprintf("[%7.1fs] ", as.numeric(difftime(Sys.time(), T0, units = "secs"))),
                         ..., "\n", sep = "")

rss_gb <- function() {
  l <- grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", l)) / 2^20
}
TIMINGS <- list()
step <- function(label, expr) {
  t <- Sys.time(); r0 <- rss_gb()
  v <- force(expr)
  el <- as.numeric(difftime(Sys.time(), t, units = "secs"))
  TIMINGS[[length(TIMINGS) + 1L]] <<- data.frame(
    step = label, sec = round(el, 1), rss_before_gb = round(r0, 2),
    rss_after_gb = round(rss_gb(), 2), stringsAsFactors = FALSE)
  log(sprintf("%-34s %8.1fs   RSS %.1f → %.1f GB", label, el, r0, rss_gb()))
  v
}

# ---- 1. 参数解析 -----------------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(argv) >= 1) argv[1] else ""
if (!MODE %in% c("smoke", "full", "subset"))
  stop("用法: 10_seurat_traditional.R {smoke --n_cells N | full | subset --cells FILE --tag TAG}")
N_SMOKE    <- NA_integer_
CELLS_FILE <- NA_character_
SUB_TAG    <- NA_character_
SUB_LIN    <- NA_character_   # 谱系名，写进 manifest；必须显式给，不许推断
if (MODE == "smoke") {
  i <- which(argv == "--n_cells")
  if (length(i) != 1 || length(argv) < i + 1) stop("smoke 模式须给 --n_cells N")
  N_SMOKE <- as.integer(argv[i + 1])
  if (is.na(N_SMOKE) || N_SMOKE < 2000) stop("n_cells 须 >= 2000（少于此时 HVG 与 PCA 不稳）")
}
if (MODE == "subset") {
  # GP8a/GP8b 亚聚类：细胞集合是**唯一**改动变量，其余口径逐位同 full（PARAMETERS §M3-A.5）
  i <- which(argv == "--cells"); j <- which(argv == "--tag"); k <- which(argv == "--lineage")
  if (length(i) != 1 || length(j) != 1 || length(k) != 1)
    stop("subset 模式须给 --cells FILE --tag TAG --lineage 谱系名")
  CELLS_FILE <- argv[i + 1]; SUB_TAG <- argv[j + 1]; SUB_LIN <- argv[k + 1]
  if (!file.exists(CELLS_FILE)) stop("找不到细胞清单: ", CELLS_FILE)
  if (!grepl("^[A-Za-z0-9_.-]+$", SUB_TAG)) stop("--tag 只允许 [A-Za-z0-9_.-]：", SUB_TAG)
  # 谱系名须是六谱系之一（法则2 表），挡住手滑写错
  .LIN_OK <- c("上皮", "T/NK", "B/浆", "髓系", "成纤维", "内皮")
  if (!SUB_LIN %in% .LIN_OK)
    stop(sprintf("--lineage 须是 %s 之一，收到：%s", paste(.LIN_OK, collapse = "/"), SUB_LIN))
}
TAG <- switch(MODE, full = "full", smoke = sprintf("smoke_%d", N_SMOKE), subset = SUB_TAG)
OUT <- file.path(OUTBASE, TAG)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "embeddings"), showWarnings = FALSE)
log(sprintf("模式=%s  tag=%s", MODE, TAG))

# ---- 2. 读入 CSC 三元组 ----------------------------------------------------
man_in <- fromJSON(file.path(IO, "export_manifest.json"))
src    <- man_in$source_h5ad
NNZ     <- as.numeric(src$nnz)
N_CELLS <- as.integer(src$n_cells)
N_GENES <- as.integer(src$n_genes)
log(sprintf("manifest：%d 细胞 × %d 基因，nnz=%.0f", N_CELLS, N_GENES, NNZ))

MAT <- step("读入 CSC 三元组 + 组装 dgCMatrix", {
  d  <- readBin(file.path(IO, "data_f32.bin"),    what = "numeric", n = NNZ,           size = 4)
  ix <- readBin(file.path(IO, "indices_i32.bin"), what = "integer", n = NNZ,           size = 4)
  p  <- readBin(file.path(IO, "indptr_i32.bin"),  what = "integer", n = N_CELLS + 1L,  size = 4)
  stopifnot(p[length(p)] == NNZ, length(ix) == NNZ, length(d) == NNZ)
  m <- new("dgCMatrix", i = ix, p = p, x = d, Dim = c(N_GENES, N_CELLS))
  rm(d, ix, p); invisible(gc())
  rownames(m) <- readLines(file.path(IO, "gene_names.txt"))
  colnames(m) <- readLines(file.path(IO, "cell_names.txt"))
  m
})
stopifnot(nrow(MAT) == N_GENES, ncol(MAT) == N_CELLS)

META <- read.csv(gzfile(file.path(IO, "cell_meta.csv.gz")), stringsAsFactors = FALSE)
stopifnot(identical(META$cell_barcode, colnames(MAT)))    # 列序必须严格一致
stopifnot(all(META$stage %in% STAGES))
rownames(META) <- META$cell_barcode
log(sprintf("cell_meta 列序一致；%d 样本 / %d 患者 / %d 分期",
            length(unique(META$sample_id)), length(unique(META$patient_id)),
            length(unique(META$stage))))

# ---- 3. smoke 抽样（按 sample_id 分层，跑前固定种子）-----------------------
if (MODE == "smoke") {
  set.seed(SEED_SUBSAMPLE)
  frac <- min(1, N_SMOKE / N_CELLS)
  keep <- unlist(lapply(split(seq_len(N_CELLS), META$sample_id), function(idx) {
    if (length(idx) == 0) return(integer(0))
    n <- max(1L, round(length(idx) * frac))
    sort(sample(idx, min(n, length(idx))))
  }), use.names = FALSE)
  keep <- sort(unique(keep))
  log(sprintf("smoke 分层抽样：%d → %d 细胞（目标 %d，frac=%.4f，seed=%d）",
              N_CELLS, length(keep), N_SMOKE, frac, SEED_SUBSAMPLE))
  MAT  <- MAT[, keep, drop = FALSE]
  META <- META[keep, , drop = FALSE]
  stopifnot(ncol(MAT) == nrow(META))
}

# ---- 3b. subset 过滤（GP8a/GP8b 亚聚类）------------------------------------
# 细胞集合按用户 2026-09-17 签字的口径给定（PARAMETERS §M3-A.5），本脚本只**执行**该口径，
# 不在此处重新推导谱系。清单里出现矩阵中不存在的 barcode 即硬报错 —— 不许静默丢细胞，
# 否则亚聚类的实际规模会与登记值不符而无人察觉。
if (MODE == "subset") {
  want <- readLines(CELLS_FILE)
  want <- want[nzchar(want)]
  if (anyDuplicated(want)) stop("细胞清单含重复 barcode")
  miss <- setdiff(want, colnames(MAT))
  if (length(miss)) stop(sprintf("清单 %d 个 barcode 不在矩阵中（前 3 个：%s）",
                                 length(miss), paste(head(miss, 3), collapse = ", ")))
  keep <- match(want, colnames(MAT))          # 保持清单顺序，便于与登记值对撞
  log(sprintf("subset 过滤：%d → %d 细胞（清单 %d 条，全部命中）",
              N_CELLS, length(keep), length(want)))
  MAT  <- MAT[, keep, drop = FALSE]
  META <- META[keep, , drop = FALSE]
  stopifnot(ncol(MAT) == nrow(META), ncol(MAT) == length(want))
  log(sprintf("subset 规模核对：%d 样本 / %d 患者 / %d 分期",
              length(unique(META$sample_id)), length(unique(META$patient_id)),
              length(unique(META$stage))))
}

# ---- 4. Seurat 对象 --------------------------------------------------------
# ⚠️ min.cells / min.features 显式给 0：默认值虽也为 0，但留空会被后续维护改成
#    非 0 而**静默丢基因/细胞**，破坏与论文基因集的可比性。
obj <- step("CreateSeuratObject", CreateSeuratObject(
  counts = MAT, project = "GSE308103", assay = "RNA",
  meta.data = META[, c("sample_id", "patient_id", "stage", "stage_token")],
  min.cells = 0, min.features = 0))
stopifnot(ncol(obj) == ncol(MAT), nrow(obj) == N_GENES)
rm(MAT); invisible(gc())
log(sprintf("Seurat 对象：%d 细胞 × %d 基因", ncol(obj), nrow(obj)))

# ---- 5. SCTransform（论文：归一化 + 缩放）----------------------------------
CLIP <- sqrt(ncol(obj) / 30)      # Seurat 默认 clip.range = ±sqrt(ncol/30)
# ⚠️ 该值随 ncol 变：全量 413,697 ⇒ ±117.43；小样本烟雾测试时会变小（2999 ⇒ ±10.00，
#    纯属巧合，与 scanpy 的 max_value=10 无关，勿据此认为两口径等价）。
log(sprintf("clip.range = ±%.2f（= ±sqrt(%d/30)，Seurat 默认公式；全量下应为 ±117.43）",
            CLIP, ncol(obj)))
obj <- step("SCTransform(vst.flavor=v2, HVG=3000)", SCTransform(
  obj, assay = "RNA", new.assay.name = "SCT",
  vst.flavor = VST_FLAVOR,
  do.correct.umi = TRUE,
  ncells = N_CELLS_MODEL,
  residual.features = NULL,
  variable.features.n = N_HVG,
  variable.features.rv.th = RV_TH,
  vars.to.regress = NULL,
  do.scale = FALSE, do.center = TRUE,
  clip.range = c(-CLIP, CLIP),
  # ⚠️ 实测内存代价（2026-09-16，全量 413,697 细胞，256 GB 机器）：
  #   以下两项组合使本步峰值 RSS 达 **242 GB**（可用内存一度只剩 9 GB，无 swap）。
  #   · do.correct.umi=TRUE     ⇒ sctransform::correct() 先建 18047×413697 稠密矩阵
  #                               （≈59.7 GB）再转稀疏；SCT 的 counts/data 槽用它，
  #                               下游 FindAllMarkers(slot="data") 会读到 ⇒ **不能关**
  #                               （关掉会改变 GP6 注释的 marker 检验，非等价）。
  #   · conserve.memory=FALSE   ⇒ 另留全基因残差矩阵（同尺寸）。置 TRUE 可省这一份，
  #                               但会改变 vst 路径，本 run 未采用。
  #   若机器内存 < 256 GB，本步会 OOM ⇒ 须先降采样或改用 conserve.memory=TRUE 并重登记。
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = SEED_SCT,
  verbose = TRUE))
HVG <- VariableFeatures(obj, assay = "SCT")
log(sprintf("HVG = %d（论文 3000）", length(HVG)))
stopifnot(length(HVG) == N_HVG)   # 硬断言：不足 3000 说明基因集偏小，须显式处理而非沉默

# ---- 6. PCA（论文：top 50 PCs）---------------------------------------------
obj <- step(sprintf("RunPCA(npcs=%d)", N_PCS), RunPCA(
  obj, assay = "SCT", features = HVG, npcs = N_PCS, seed.use = SEED_PCA, verbose = FALSE))
stopifnot(ncol(Embeddings(obj, "pca")) == N_PCS)

# ---- 7. Harmony（论文主口径）-----------------------------------------------
obj <- step("RunHarmony(group.by.vars=sample_id, dims=1:50)", RunHarmony(
  obj, group.by.vars = "sample_id", reduction.use = "pca",
  dims.use = seq_len(N_PCS), reduction.save = "harmony", project.dim = TRUE))
stopifnot(identical(rownames(Embeddings(obj, "harmony")), colnames(obj)))

# ---- 8. 两臂的 kNN 图 ------------------------------------------------------
run_arm <- function(obj, space, graph_names) {
  FindNeighbors(obj, reduction = space, dims = seq_len(N_PCS),
                k.param = K_PARAM, annoy.metric = "euclidean",
                nn.method = "annoy", prune.SNN = PRUNE_SNN,
                graph.name = graph_names, verbose = FALSE)
}
obj <- step("FindNeighbors(Harmony 臂)", run_arm(obj, "harmony", c("harmony_nn", "harmony_snn")))
obj <- step("FindNeighbors(不校正臂)", run_arm(obj, "pca",     c("pca_nn", "pca_snn")))

# ---- 9. 聚类：4 分辨率 × 5 种子（Harmony 臂）+ 4 分辨率 × 种子0（不校正臂）--
cluster_at <- function(obj, graph, r, seed) {
  obj <- FindClusters(obj, graph.name = graph, resolution = r, algorithm = 1,
                      modularity.fxn = 1, n.start = 10, n.iter = 10,
                      random.seed = seed, group.singletons = TRUE, verbose = FALSE)
  as.integer(as.character(obj$seurat_clusters))
}
LAB <- list()
for (r in RESOLUTIONS) {
  for (s in SEEDS) {
    LAB[[sprintf("harmony_res%.1f_seed%d", r, s)]] <-
      step(sprintf("FindClusters(harmony, r=%.1f, seed=%d)", r, s),
           cluster_at(obj, "harmony_snn", r, s))
  }
}
for (r in RESOLUTIONS) {
  LAB[[sprintf("nocorrect_res%.1f_seed0", r)]] <-
    step(sprintf("FindClusters(pca, r=%.1f, seed=0)", r), cluster_at(obj, "pca_snn", r, 0L))
}

# ---- 10. UMAP（Harmony 空间）-----------------------------------------------
obj <- step(sprintf("RunUMAP(reduction=harmony, dims=1:%d)", N_PCS), RunUMAP(
  obj, reduction = "harmony", dims = seq_len(N_PCS),
  n.neighbors = 30L, min.dist = 0.3, metric = "cosine",
  umap.method = "uwot", spread = 1, n.components = 2L,
  seed.use = SEED_UMAP, verbose = FALSE))
UMAP <- Embeddings(obj, "umap")

# ---- 11. 指标 --------------------------------------------------------------
# ARI（Hubert–Arabie），自实现以免引入依赖；与 mclust::adjustedRandIndex 同式
ARI <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  if (n < 2) return(NA_real_)
  c2 <- function(x) x * (x - 1) / 2
  sij <- sum(c2(tab)); si <- sum(c2(rowSums(tab))); sj <- sum(c2(colSums(tab)))
  exp_v <- si * sj / c2(n); mx <- (si + sj) / 2
  if (isTRUE(all.equal(mx, exp_v))) return(1)
  (sij - exp_v) / (mx - exp_v)
}

# 指标1：跨种子稳定性（每个分辨率内 5 种子两两 ARI 均值）
seed_stab <- vapply(RESOLUTIONS, function(r) {
  ls <- lapply(SEEDS, function(s) LAB[[sprintf("harmony_res%.1f_seed%d", r, s)]])
  m <- combn(length(ls), 2, function(k) ARI(ls[[k[1]]], ls[[k[2]]]))
  mean(m)
}, numeric(1))

# 指标2：跨分辨率稳定性（相邻分辨率 ARI，种子0）
adj <- vapply(seq_along(RESOLUTIONS), function(i) {
  if (i == 1) return(NA_real_)
  ARI(LAB[[sprintf("harmony_res%.1f_seed0", RESOLUTIONS[i - 1])]],
      LAB[[sprintf("harmony_res%.1f_seed0", RESOLUTIONS[i])]])
}, numeric(1))
xres <- vapply(seq_along(RESOLUTIONS), function(i) {
  v <- c(if (i > 1) adj[i] else NA_real_, if (i < length(RESOLUTIONS)) adj[i + 1] else NA_real_)
  mean(v, na.rm = TRUE)
}, numeric(1))

# 指标3：谱系覆盖（**不在本脚本计算** —— 见 05_annotation/04_metric3_coverage.py。
#        此处**只报簇数**，不得以簇数冒充覆盖分）
n_clust <- vapply(RESOLUTIONS, function(r) length(unique(LAB[[sprintf("harmony_res%.1f_seed0", r)]])),
                  integer(1))

# 指标4：AAH 吸收护栏（本项目约定）
STG <- factor(META$stage, levels = STAGES)
absorption <- function(cl) {
  tb <- table(cl, STG)
  prop_norm <- tb[, "Normal"] / rowSums(tb)
  major <- colnames(tb)[max.col(tb, ties.method = "first")]
  norm_major <- (major == "Normal") & (prop_norm > 0.5)
  aah <- STG == "AAH"
  if (!any(aah)) return(NA_real_)
  sum(cl[aah] %in% as.integer(rownames(tb))[norm_major]) / sum(aah)
}
absorb <- vapply(RESOLUTIONS, function(r) absorption(LAB[[sprintf("harmony_res%.1f_seed0", r)]]),
                 numeric(1))

# ⚠️ 2026-09-17 修：原实现拿**未舍入**的 seed_stab/absorb 判 pass_*，却把 round(...,4) 写进 CSV。
#   ⇒ CSV 里可能出现 `ari_seed_mean = 0.9000` 而 `pass_seed = FALSE` 的自相矛盾行（0.89996 舍入后
#   显示 0.9000 却仍判 FALSE）；而 03_subcluster_annotation.py 读的正是**舍入后**那一列，
#   边界样本上两侧会给出相反的 r*。
#   现改为「先舍入 → 再用舍入值算 score 与 pass_*」⇒ CSV 自洽，R 与 Python 口径逐位一致。
seed_stab_r <- round(seed_stab, 4)
xres_r      <- round(xres, 4)
adj_r       <- round(adj, 4)
absorb_r    <- round(absorb, 4)
score_r     <- round(0.5 * seed_stab_r + 0.5 * xres_r, 4)
pass_seed_r <- seed_stab_r >= ARI_SEED_MIN
pass_aah_r  <- absorb_r <= AAH_ABSORPTION_MAX
ok          <- pass_seed_r & pass_aah_r          # 指标3 不由本脚本施加（见 04_metric3_coverage.py）
METRICS <- data.frame(
  resolution = RESOLUTIONS, n_clusters = n_clust,
  ari_seed_mean = seed_stab_r,
  ari_adjacent_prev = adj_r, ari_xres = xres_r,
  aah_absorption_rate = absorb_r,
  score = score_r,
  pass_seed = pass_seed_r, pass_aah = pass_aah_r,
  eligible = ok, stringsAsFactors = FALSE)
RSTAR <- if (any(ok)) {
  # 预注册（PARAMETERS §M3-A.3）：r* = argmax score，受指标 1/3/4 三条硬约束；
  #   破平规则（显式）：两分辨率分差 < 0.01 ⇒ 取**较低**者（少簇、少过度切分）。
  # ⚠️ 2026-09-17 修复：本行原为裸 `which.max(e$score)`，注释却声称实现了破平 —— 注释与代码不符，
  #   全量 run 因此把 r*=0.7 报出，而按规则应为 r=0.6（0.9123−0.9041=0.0082 < 0.01）。
  #   已改为下方的显式实现；该 run 的产物由 run_manifest.json 的 rstar_correction 块更正。
  e <- METRICS[METRICS$eligible, ]
  e <- e[order(-e$score), ]
  # ⚠️ 2026-09-17 修：原为严格 `>`，而 Python 侧 03_subcluster_annotation.py 的同一规则用 `>=`
  #   （`score >= best - TIE_BREAK_EPS`）。恰好在差 0.01 时会分出**不同**的 r*，
  #   且 resolution_metrics.csv 的 is_rstar 列与 rstar.json 会互相矛盾。统一为 `>=`。
  tied <- e$resolution[e$score >= e$score[1] - 0.01]  # 与最高分差 ≤0.01 者视为并列
  min(tied)                                           # 取较低分辨率
} else NA_real_
METRICS$is_rstar <- !is.na(RSTAR) & METRICS$resolution == RSTAR

# ---- 12. 落盘 --------------------------------------------------------------
CL <- data.frame(cell_barcode = META$cell_barcode, sample_id = META$sample_id,
                 patient_id = META$patient_id, stage = META$stage,
                 as.data.frame(LAB, check.names = FALSE), stringsAsFactors = FALSE)
write.csv(CL, gzfile(file.path(OUT, "clusters.csv.gz")), row.names = FALSE, quote = FALSE)
write.csv(METRICS, file.path(OUT, "resolution_metrics.csv"), row.names = FALSE)
write.csv(data.frame(cell_barcode = META$cell_barcode,
                     umap_1 = round(UMAP[, 1], 5), umap_2 = round(UMAP[, 2], 5)),
          gzfile(file.path(OUT, "umap.csv.gz")), row.names = FALSE, quote = FALSE)
write.csv(do.call(rbind, TIMINGS), file.path(OUT, "timings.csv"), row.names = FALSE)

EMB <- list()
for (nm in c("harmony", "pca")) {
  e <- Embeddings(obj, nm); p <- file.path(OUT, "embeddings", sprintf("%s_f32.bin", nm))
  writeBin(as.numeric(t(e)), p, size = 4)     # 行主序展开：cell-major
  EMB[[nm]] <- list(path = file.path("embeddings", sprintf("%s_f32.bin", nm)),
                    dim = unname(dim(e)), bytes = file.info(p)$size,
                    md5 = unname(tools::md5sum(p)), seed = if (nm == "pca") SEED_PCA else NA)
}
writeLines(colnames(obj), file.path(OUT, "embeddings", "cells.txt"))

# subset 模式的清单哈希提到 man 之外算一次：manifest 的 subset 块与 metrics_caveat
# 都要用同一个值，就地算两遍容易在日后改动中走散。
MAN_SUBSET_SHA <- if (MODE == "subset")
  unname(digest::digest(file = CELLS_FILE, algo = "sha256")) else NA_character_

man <- list(
  script = "04_integration/10_seurat_traditional.R",
  tag = TAG, mode = MODE,
  source_paper = list(citation = "Peng F, Sinjab A, … Kadara H. Cancer Cell 2026;44(2):321-339.e13",
                      doi = "10.1016/j.ccell.2025.10.004", pmcid = "PMC12980502"),
  input = list(seurat_io = man_in$source_h5ad, n_cells = ncol(obj), n_genes = nrow(obj)),
  subsample = if (MODE == "smoke")
    list(n_target = N_SMOKE, n_actual = ncol(obj), stratify_by = "sample_id",
         seed = SEED_SUBSAMPLE) else NULL,
  subset = if (MODE == "subset") list(
    cells_file = CELLS_FILE,
    # R 的 tools 包**没有** sha256sum（只有 md5sum）；用 digest::digest 算 sha256。
    # 2026-09-17 首跑 6 谱系时此处崩掉 6 次（subset 分支此前从未被执行过），改此。
    cells_file_sha256 = MAN_SUBSET_SHA,
    n_listed = length(keep), n_actual = ncol(obj),
    lineage = SUB_LIN,
    caliber = "标准 A 单口径 = GP6 冻结 A_frozen 标签（seed0, r*=0.6），用户 2026-09-17 裁定；PARAMETERS §M3-A.5",
    n_samples = length(unique(META$sample_id)), n_patients = length(unique(META$patient_id)),
    stage_counts = as.list(table(META$stage))) else NULL,
  params = list(
    sctransform = list(vst.flavor = VST_FLAVOR, variable.features.n = N_HVG,
                       variable.features.rv.th = RV_TH, ncells = N_CELLS_MODEL,
                       vars.to.regress = NULL, do.scale = FALSE, do.center = TRUE,
                       clip.range = c(-CLIP, CLIP), conserve.memory = FALSE,
                       return.only.var.genes = TRUE, seed.use = SEED_SCT),
    pca      = list(npcs = N_PCS, seed.use = SEED_PCA, features = "SCT VariableFeatures (3000)"),
    harmony  = list(group.by.vars = "sample_id", reduction.use = "pca",
                    dims.use = seq_len(N_PCS), reduction.save = "harmony",
                    project.dim = TRUE, note = "其余为 harmony 2.0.5 默认"),
    neighbors = list(k.param = K_PARAM, annoy.metric = "euclidean", nn.method = "annoy",
                     prune.SNN = PRUNE_SNN, graphs = c("harmony_snn", "pca_snn")),
    clusters  = list(resolutions = RESOLUTIONS, seeds = SEEDS, algorithm = 1L,
                     modularity.fxn = 1L, n.start = 10L, n.iter = 10L,
                     group.singletons = TRUE, method = "matrix"),
    umap      = list(n.neighbors = 30L, min.dist = 0.3, metric = "cosine",
                     umap.method = "uwot", spread = 1, n.components = 2L,
                     seed.use = SEED_UMAP)),
  arms = list(
    main = "harmony（论文配方）",
    sensitivity = "不校正 pca（**不属于论文配方**，共线风险下的稳健性对照）"),
  metrics = list(
    guard = list(aah_absorption_max = AAH_ABSORPTION_MAX,
                 note = "指标4 为本项目附加，非论文配方；见 PARAMETERS §M3-A.3"),
    metric3_absent = paste0("指标3（谱系覆盖）**不在本脚本计算** —— 见 05_annotation/04_metric3_coverage.py ",
                            "与产物 results/05_annotation/metric3_coverage.csv（2026-09-21 补算：",
                            "七对象 × 4 分辨率 × 5 种子全过线，对 r* 无区分力）。",
                            "⚠️ 旧措辞「需 GP6 注释后方可计算」是错的（它不依赖注释）。",
                            "此处只报 n_clusters，不得以簇数冒充覆盖分"),
    metric1_caveat = paste0(
      # ⚠️ 本字段**由 METRICS 实算**，禁止写死结论。
      #    旧版曾把 3,000 细胞冒烟的"ARI 恒 = 1.0 ⇒ 指标1 不承重"写死在此，
      #    全量跑完后与实际 0.8963–0.9230 自相矛盾（2026-09-17 发现并改为实算）。
      sprintf("指标1（跨种子稳定性，硬约束 ARI ≥ %.2f）：本轮实测 %.4f–%.4f，",
              ARI_SEED_MIN, min(METRICS$ari_seed_mean), max(METRICS$ari_seed_mean)),
      sprintf("%d/%d 个分辨率达标", sum(METRICS$pass_seed), nrow(METRICS)),
      if (!all(METRICS$pass_seed)) sprintf(
        "；未达标：%s ⇒ 该约束在本规模上**承重**，已剔除这 %d 个候选。",
        paste(sprintf("r=%.1f(ARI=%.4f)", METRICS$resolution[!METRICS$pass_seed],
                      METRICS$ari_seed_mean[!METRICS$pass_seed]), collapse = "、"),
        sum(!METRICS$pass_seed)) else
        "；全部达标 ⇒ 该约束在本规模上未剔除任何候选。",
      " 护栏效力随规模而变（历史实测：3,000 细胞空转，ARI 恒 = 1.0；",
      "60,000 细胞 0.9346–0.9834）——不得用小规模结果推断本规模。",
      "⚠️ 本字段由 METRICS 实算，若与 METRICS 表不符以 METRICS 表为准。"),
    rstar_candidate = RSTAR,
    rstar_status = if (is.na(RSTAR)) "区间内无可行解 ⇒ 停在 GP5，不放宽阈值"
                   else "候选，**待 GP5 人工签字**"),
  metrics_interpretable = (MODE %in% c("full", "subset")),
  metrics_caveat = if (MODE == "smoke")
    paste("smoke 为**分层子样本**（每样本细胞数远小于真实值）⇒ Harmony 与 Louvain 的行为",
          "不代表全量。本轮的 resolution_metrics / aah_absorption **不得用于挑选 r\\***，",
          "仅用于验证代码路径与测时。") else if (MODE == "subset")
    # ⚠️ 2026-09-17 修：本行原写死「上皮交集」—— 但 subset 是六谱系**共用**的分支，
    #    五个非上皮谱系的 manifest 里都印着这句错话（同 `lineage = "上皮"` 那类漏检）。
    #    现改为按 --lineage 实参生成，并把登记在案的口径与清单哈希一并写出。
    sprintf(paste("本轮为**亚聚类**：谱系=%s，细胞集合 = %s（标准 A 单口径，seed0 全局 r*=0.6 的该谱系全部细胞），",
                  "清单 sha256 = %s。metric4（AAH 吸收护栏）的语义是全量口径的**投影**而非同一定义：",
                  "本谱系子集里 AAH 细胞仍存在，但簇的分期构成已随谱系而变 ⇒ 该列**不与全量值直接比较**。",
                  "其余口径（SCTransform 参数 / PCA / Harmony / kNN / 网格 / 种子）与 full 逐位相同。"),
            SUB_LIN, basename(CELLS_FILE), substr(MAN_SUBSET_SHA, 1, 16)) else NULL,
  versions = list(R = R.version.string, Seurat = as.character(packageVersion("Seurat")),
                  SeuratObject = as.character(packageVersion("SeuratObject")),
                  sctransform = as.character(packageVersion("sctransform")),
                  harmony = as.character(packageVersion("harmony")),
                  glmGamPoi = as.character(packageVersion("glmGamPoi")),
                  Matrix = as.character(packageVersion("Matrix"))),
  embeddings = EMB,
  metrics_table = METRICS,
  wall_sec = round(as.numeric(difftime(Sys.time(), T0, units = "secs")), 1))
write_json(man, file.path(OUT, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)
log(sprintf("完成 -> %s", OUT))
print(METRICS)
log(sprintf("r* 候选 = %s（待 GP5 人工签字）", ifelse(is.na(RSTAR), "无（区间内无可行解）", as.character(RSTAR))))
log(sprintf("总耗时 %.1f 分钟", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
