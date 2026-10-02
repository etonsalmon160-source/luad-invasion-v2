#!/usr/bin/env Rscript
# 01_hvg_freeze.R —— M6 生态位臂 §3.1 / S3：全队列选 HVG 3000，落盘冻结
#
# 口径（NICHE_PREREG.md §3.1，2026-10-01 已签；同日修正，见 §3.1.1）：
#   基因集 = HVG 3000，**全队列 56 张只选一次**；清单落盘 results/10_niche/hvg_3000.txt，
#   所有切片共用同一套基因。
#   ⚠️「只选一次」是硬约束：逐切片各选 HVG 会让每张切片基因轴不同，
#      §3.5 的跨切片共识（域表达画像）在数学上就不成立（向量长度都不一样）。
#
# 🔴 §3.1.1 修正（2026-10-01，用户签字）：原签做法「56 张拼成大矩阵后 FindVariableFeatures(vst)」
#   **物理上做不到**——合池矩阵非零元 ≈ 21 亿，超过 R `Matrix` 稀疏矩阵 `p` 槽的 `2^31-1` 硬上限
#   （dgCMatrix 没有 64 位索引版），`do.call(cbind, mats)` 直接 `Execution halted`。
#   ⇒ 改为 **逐切片各自 FindVariableFeatures(vst) → Seurat 官方 `SelectIntegrationFeatures`**：
#      按「一个基因在多少张切片里入选为 HVG」降序排，平局按跨切片中位排名破。
#   ⚠️ 代价（须进限制）：打分是**逐切片做的**（各自的 mean-variance loess 拟合），
#      不是全队列一次拟合 —— 与字面的"合池 vst"不同，是**另一个方法**，不是同一个的近似。
#
# 基因轴 = 56 张切片 symbol 的**交集** = 18,082（纯空转基因宇宙）。
#   🔴 别和 18,066 混：18,066 是 **RCTD 臂**的 `intersect(参考基因, 切片基因)`
#      （见 08_spatial_deconv/01_run_rctd.R G1 条款，逐张一致，落在 rctd_*/run_manifest.json）。
#      两个集合不同、都真实：主轴用纯空转交集 18,082；18,066 只作 §3.1 敏感性臂的登记参照。
#   本脚本把交集基因数打出来核对；**若 ≠ 18,082 就硬停**，不许静默带着不同的轴继续。
#
# N0 = 冻结掩码通过集 608,043 spot（§1 硬约束，与 01_run_rctd.R / 19_run_fastcnv_cohort.R 同源）。
# 🔴 本脚本**只选基因**，不产生任何生物学结论；不碰 RCTD 权重，不做任何聚类。
#
# 跑法（Banksy/fastCNV 那个 R 库栈；Seurat 在里面）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   setsid nohup Rscript 10_niche/01_hvg_freeze.R < /dev/null > logs/hvg_freeze.log 2>&1 &
#   🔴 不能用 `Rscript --vanilla`（会丢掉 /usr/local/lib/R/site-library 里的 data.table）。

suppressMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
})

ROOT   <- "/home/eto/luad_v2"
VISIUM <- file.path(ROOT, "data/visium_spatial")
RES    <- file.path(ROOT, "results/08_spatial_deconv")
OUTD   <- file.path(ROOT, "results/10_niche")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

N_HVG        <- 3000L      # [已签 S3]
N_FVF        <- 3000L      # 逐切片候选 HVG 数（喂给 SelectIntegrationFeatures 的 fvf.nfeatures）
NRM_FACTOR   <- 1e4        # [已签 S2] CP10K
N_GENE_EXPECT<- 18082L     # 56 张切片 symbol 交集（本臂主轴）；对不上就硬停
N_GENE_RCTD  <- 18066L     # RCTD 臂的 参考∩切片 交集（§3.1 敏感性臂），只登记不设门
MIN_FREE_GB  <- 60         # 内存守卫线（56 个 Seurat 对象常驻；背景还有病理进程）

t_start <- Sys.time()
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())
meminfo <- function() {
  x <- readLines("/proc/meminfo")
  g <- function(k) as.numeric(sub("^[^:]+:\\s*([0-9]+).*$", "\\1",
        x[grepl(paste0("^", k, ":"), x)][1])) / 1048576
  c(avail = g("MemAvailable"))
}
free_gb <- function() unname(meminfo()["avail"])
read_gz <- function(p, fn) { con <- gzfile(p, "rt"); on.exit(close(con)); fn(con) }

step("==== §3.1 S3：全队列 HVG %d 冻结（SelectIntegrationFeatures 版）====", N_HVG)
step("Seurat %s / R %s / 可用内存 %.1f GB",
     packageVersion("Seurat"), paste(R.version$major, R.version$minor, sep = "."), free_gb())
if (free_gb() < MIN_FREE_GB)
  stop(sprintf("可用内存 %.1f GB < 守卫线 %.0f GB ⇒ 硬停", free_gb(), MIN_FREE_GB), call. = FALSE)

## ——— 1. 冻结掩码通过集（N0）———
sm_all <- read.delim(gzfile(file.path(RES, "spot_mask.tsv.gz")), stringsAsFactors = FALSE)
sm_all <- sm_all[as.character(sm_all$pass) %in% c("TRUE", "true", "1"), ]
SLIDES <- sort(unique(sm_all$slide))
stopifnot(length(SLIDES) == 56L)
step("冻结掩码通过集 %d spot / %d 张切片", nrow(sm_all), length(SLIDES))

## 读一张切片的基因符号（剔抗体行；抗体行计数量级 1e6，必剔）
feat_of <- function(s) {
  ft <- read_gz(file.path(VISIUM, s, "filtered_feature_bc_matrix", "features.tsv.gz"),
                function(con) read.delim(con, header = FALSE,
                                         col.names = c("id", "sym", "type"),
                                         stringsAsFactors = FALSE))
  ab <- which(ft$type == "Antibody Capture")
  ## 抗体行（计数量级 1e6）必须剔；返回的 ft 已剔，ab 是**原始行号**（用于剔矩阵的同一批行）。
  ## 实测 45 张各有 35 行抗体（HLA-DRA / mouse_IgG* / rat_IgG2a），但那些符号其余切片没有，
  ## 因此剔与不剔，56 张 symbol 交集都是 18,082（已独立复核）。
  if (length(ab)) ft <- ft[-ab, , drop = FALSE]
  list(ft = ft, ab = ab)
}

## ——— 2. 基因轴 = 56 张交集 ———
step("求 56 张 symbol 交集 …")
gsets <- lapply(SLIDES, function(s) unique(feat_of(s)$ft$sym))
GENES <- Reduce(intersect, gsets)
step("交集基因数 = %d（本臂主轴登记 %d；RCTD 臂交集 %d，另计）",
     length(GENES), N_GENE_EXPECT, N_GENE_RCTD)
if (length(GENES) != N_GENE_EXPECT)
  stop(sprintf("交集基因数 %d ≠ 登记的 %d ⇒ 硬停核查（轴对不上，下游共识不成立）",
               length(GENES), N_GENE_EXPECT), call. = FALSE)
rm(gsets); invisible(gc())

## ——— 3. 逐切片读入 → 取 mask-pass spot → 对齐到冻结基因轴 ———
read_slide_on_axis <- function(s) {
  fi <- feat_of(s); ft <- fi$ft; ab <- fi$ab
  bc_keep <- sm_all$barcode[sm_all$slide == s]
  stopifnot(length(bc_keep) > 0)

  bc <- read_gz(file.path(VISIUM, s, "filtered_feature_bc_matrix", "barcodes.tsv.gz"), readLines)
  m  <- readMM(file.path(VISIUM, s, "filtered_feature_bc_matrix", "matrix.mtx.gz"))
  ## feat_of() 已把 ft 剔好，但 m 仍是**含抗体行**的原始行序 ⇒ 这里按同一批行号剔 m。
  ## （曾因把这里写成 nrow(m)==nrow(ft) 而崩过一次：两边的行数口径不一致。)
  stopifnot(nrow(m) == nrow(ft) + length(ab), ncol(m) == length(bc))
  if (length(ab)) m <- m[-ab, , drop = FALSE]
  stopifnot(nrow(m) == nrow(ft))

  j <- match(bc_keep, bc)
  if (anyNA(j)) stop(s, "：mask-pass barcode 不在表达矩阵里", call. = FALSE)
  m <- m[, j, drop = FALSE]

  ## 同名基因合并（与 00_banksy_smoke.R / 19_run_fastcnv_cohort.R 同法）
  g <- factor(ft$sym, levels = unique(ft$sym))
  A <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1,
                    dims = c(nlevels(g), length(g)))
  m <- A %*% m
  rownames(m) <- levels(g)
  colnames(m) <- bc_keep

  ## 对齐到冻结基因轴（交集，故必定全部命中）
  if (!all(GENES %in% rownames(m))) stop(s, "：冻结基因轴里有本片缺的基因", call. = FALSE)
  as(m[GENES, , drop = FALSE], "CsparseMatrix")
}

## ——— 4. 逐切片建对象 → CP10K+log1p → 逐切片 FindVariableFeatures(vst) ———
## 🔴 这里**必须**逐张独立建对象：把 56 张拼成一个大矩阵会触发 p 槽 2^31-1 上限（见文件头 §3.1.1）。
step("逐切片 NormalizeData(CP10K) + FindVariableFeatures(vst, n=%d) …", N_FVF)
obj.list  <- vector("list", length(SLIDES))
nvar_each <- integer(length(SLIDES))
for (i in seq_along(SLIDES)) {
  s <- SLIDES[i]
  m <- read_slide_on_axis(s)
  o <- CreateSeuratObject(counts = m, min.cells = 0, min.features = 0)
  rm(m); invisible(gc())
  o <- NormalizeData(o, normalization.method = "LogNormalize",
                     scale.factor = NRM_FACTOR, verbose = FALSE)
  o <- FindVariableFeatures(o, selection.method = "vst", nfeatures = N_FVF, verbose = FALSE)
  obj.list[[i]] <- o
  nvar_each[i] <- length(VariableFeatures(o))
  if (i %% 8 == 0 || i == length(SLIDES))
    step("  打分 %d/%d 张（最近一张 %d spot / %d HVG；当前可用 %.1f GB）",
         i, length(SLIDES), ncol(o), nvar_each[i], free_gb())
}
stopifnot(all(nvar_each > 0))
step("逐切片候选 HVG 数：中位 %d，范围 %d–%d（共 %d 张）",
     as.integer(median(nvar_each)), min(nvar_each), max(nvar_each), length(nvar_each))

## ——— 5. 跨切片汇总选 3000（Seurat 官方 SelectIntegrationFeatures）———
step("SelectIntegrationFeatures(nfeatures=%d, fvf.nfeatures=%d) …", N_HVG, N_FVF)
hv <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = N_HVG,
                                fvf.nfeatures = N_FVF, verbose = TRUE)
## 「该基因在多少张切片里是 HVG」——SelectIntegrationFeatures 的排序依据，留档可查
vf_counts <- table(unlist(lapply(obj.list, VariableFeatures)))
step("选得 HVG %d 个（可用 %.1f GB）", length(hv), free_gb())
if (length(hv) != N_HVG)
  stop(sprintf("选得 %d 个 ≠ 目标 %d ⇒ 硬停（平局太多？须核查，不许少给也不许多给）",
               length(hv), N_HVG), call. = FALSE)
if (anyDuplicated(hv)) stop("HVG 列表有重复 ⇒ 硬停", call. = FALSE)
if (!all(hv %in% GENES)) stop("HVG 列表里有不在冻结基因轴上的基因 ⇒ 硬停", call. = FALSE)

## ——— 6. 落盘冻结 ———
writeLines(hv, file.path(OUTD, "hvg_3000.txt"))
audit <- data.frame(gene = hv, rank = seq_along(hv),
                    n_slide_variable = as.integer(vf_counts[hv]), row.names = NULL)
fwrite(audit, file.path(OUTD, "hvg_3000_audit.tsv"), sep = "\t")
fwrite(data.table(gene = GENES), file.path(OUTD, "gene_axis_18082.txt"), sep = "\t")

man <- list(
  step = "§3.1 S3 HVG freeze",
  prereg = "10_niche/NICHE_PREREG.md §3.1 / §3.1.1",
  method = paste0('逐切片 Seurat FindVariableFeatures(selection.method="vst", nfeatures=', N_FVF,
                  ') → Seurat SelectIntegrationFeatures(nfeatures=', N_HVG,
                  ', fvf.nfeatures=', N_FVF, ')'),
  method_note = paste0("§3.1.1 修正：原签的「56 张合池后 vst」物理上做不到",
                       "（合池矩阵非零元 ≈21 亿 > R Matrix 稀疏 p 槽 2^31-1 上限）⇒ ",
                       "改为官方跨数据集选特征。打分是逐切片拟合的，非全队列一次拟合。"),
  normalization = sprintf("逐切片 LogNormalize, scale.factor=%.0f (= CP10K + log1p, S2)", NRM_FACTOR),
  pooled = TRUE, n_slides = length(SLIDES), n_spot = nrow(sm_all), n_gene_axis = length(GENES),
  gene_axis_source = "56 张切片 symbol 交集（纯空转基因宇宙）",
  sensitivity_arm_gene_axis = N_GENE_RCTD,
  sensitivity_arm_gene_axis_source =
    "RCTD 臂 intersect(参考基因, 切片基因)，见 08_spatial_deconv/01_run_rctd.R G1 条款 / rctd_*/run_manifest.json",
  n_hvg = length(hv),
  n_var_per_slide_min = min(nvar_each), n_var_per_slide_med = as.integer(median(nvar_each)),
  n_var_per_slide_max = max(nvar_each),
  hvg_min_slide_variable = as.integer(min(vf_counts[hv])),
  hvg_med_slide_variable = as.integer(median(vf_counts[hv])),
  out_hvg = "results/10_niche/hvg_3000.txt",
  seurat = as.character(packageVersion("Seurat")),
  R = paste(R.version$major, R.version$minor, sep = "."),
  bank_version = if (requireNamespace("Banksy", quietly = TRUE)) as.character(packageVersion("Banksy")) else NA,
  mask_source = "results/08_spatial_deconv/spot_mask.tsv.gz（rule = or, 500/200/0.15）",
  started = format(t_start, "%Y-%m-%d %H:%M:%S"),
  finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  minutes = round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1))
writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE),
           file.path(OUTD, "hvg_3000_manifest.json"))

step("落盘完成：hvg_3000.txt / hvg_3000_audit.tsv / gene_axis_18082.txt / hvg_3000_manifest.json")
step("（注意：HVG 清单只对「本用途」生效，受法则 3.2 保护，不得挪作他用）")
step("==== 结束，用时 %.1f 分钟 ====", man$minutes)
