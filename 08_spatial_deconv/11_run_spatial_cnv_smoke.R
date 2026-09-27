#!/usr/bin/env Rscript
# 空转（Visium）CNV 臂 —— 单患者**成本冒烟**（SPATIAL_CNV_PREREG.md §9 第 2 项）
#
# 本脚本回答的**不是**科学问题，是 §9 第 2 项那条一直没实测的事：
#     infercnv 1.23.0 在 HMM=TRUE 下，**单患者的内存峰值与墙钟耗时**是多少？
# ⇒ 排期前必须先测。所以本脚本把「建对象」「跑管线」分两段，各自计时、各自记峰值内存。
#
# ⚠️ 本脚本**不产生任何恶性判定结论**。它跑的是**成本**。
#    跑出来的 object 与终产物不落盘（除日志），除非显式 --keep。
#
# 运行前提（必须，否则 rjags 装载失败）：
#   export LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH
# 跑法：
#   LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH \
#   R_LIBS=/home/eto/Rlibs/infercnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 08_spatial_deconv/11_run_spatial_cnv_smoke.R [--object-only] [--patient P4]
#
#   --object-only  只建对象、报维度与内存，**不跑** infercnv（验证输入管道，秒级）
#   --patient ID   默认 P4（4 张切片，是队列里最大的患者 ⇒ 成本上界）
#   --keep         保留中间 .rds（默认不落盘；队列上跑会撑爆磁盘）

suppressMessages({
  library(infercnv)
  library(Matrix)
  library(rjags)
})

## —————————————————————————————————————————————————————————————
## 一、预注册常量（SPATIAL_CNV_PREREG.md §3.3「逐字照抄原文」）
## —————————————————————————————————————————————————————————————
P_CUTOFF            <- 0.1     # [已签 §3.3]  原文 methods 逐字
P_CLUSTER_BY_GROUPS <- TRUE    # [已签 §3.3]  原文
P_HMM               <- TRUE    # [已签 §3.3]  原文；本机 JAGS 4.3.2 可用 ⇒ 够得着
P_DENOISE           <- TRUE    # [已签 §3.3]  原文
P_REF_GROUP         <- "C_benign_epi"   # [已签 §11.2] 多锚并列 C 主；D 见脚本末尾说明
P_THREADS           <- 3       # [已签 §3.3]  本机 load 长期 30+，不拉高

## ⚠️ 原文**没有**记录、§3.3 也**没有**登记的参数 —— 下面显式写出取值：
##    显式写出的目的**只是**满足法则 3.1（参数可审计）。
##    🔴 那处**自相矛盾已裁定**（预注册 §13，2026-09-27 用户签字）：
##       §3.3 写「本臂不继承 analysis_mode="subclusters"」，而工具默认正是 "subclusters"。
##       原文只记了 4 个参数、未记 analysis_mode ⇒ 「未记项取默认」是一种读法，
##       但 §3.3 是本臂自己白纸黑字写了不继承，默认值那一边没有别的依据。
##       ⇒ **取 §3.3 这一边：`samples`。**
##       实测依据（§13.2）：工具默认下患者 P4（2 张切片）**24 h 未跑完**被超时杀掉，
##       峰值内存仅 7.67 GB ⇒ 不是内存问题；`subclusters` 把参考切出 181 + 观测 65 = 246 个单位
##       逐个跑 HMM/贝叶斯网，另加 9.5 h 纯绘图。取 `samples` 后单位数降到 2。
##       ⚠️ 本次**只改 analysis_mode**，**不退到 `HMM=FALSE`**（§13.4）。
P_WINDOW_LENGTH     <- 101                  # [未登记] infercnv 默认
P_ANALYSIS_MODE     <- "samples"            # [已签 §13] 不取工具默认（见上🔴）
P_NO_PLOT           <- TRUE                 # [已签 §13] 关绘图：图不是结果，不影响任何判定
P_HCLUST_METHOD     <- "ward.D2"            # [未登记] infercnv 默认
P_CHR_EXCLUDE       <- c("chrM")   # [已签 §14.3] 原文 Step1/Step3 两处**都显式写了** c("chrM")
## 🔴 2026-09-27 改签（§14.3 第 1 项）：本行**原来**是 infercnv 默认 c("chrX","chrY","chrM")，
##    被当成「未登记项取默认」处理。核原文源码后确认：原文**记了这个参数**（在代码里，不是
##    方法段里），两处 CreateInfercnvObject 都是 c("chrM")。⇒ 照原文，chrX/chrY 放回管线。
##    与 §13 那处 analysis_mode 是同一类错、方向相反：那个原文没记（取默认有理），
##    这个原文记了（取默认即偏离）。
##    ⚠️ 配套守卫（§14.3）：参考与观测**跨患者**时 chrX/Y 会把男女剂量差算成 CNV
##       ⇒ manifest 记 ref_same_patient / ref_slides；跨患者借锚者的 chrX/Y 结论降级待定。
##    改签后的基因账（逐项算，不是估的）：
##      矩阵 18,082 唯一 symbol
##        − 67   位置表里没有（未匹配上 GENCODE/HGNC，见 gene_order_spatial_unmatched.txt）
##        = 18,015 位置表行数
##        − 11   chrM（就是本行）
##        = 18,004 进 inferCNV
##    （改签前是 − 743 = chrX 716 ＋ chrY 16 ＋ chrM 11 ⇒ 进管线 17,272。差值 732。）
##    `03_cnv/16_run_infercnv_smoke.R` 当年**没写**这个参数 ⇒ 吃了默认，属另一臂的历史，不改。

## 口径（已签 §11 / §12.1）
CAL_REFERENCE  <- "a"          # [已签 §11 第 1 项] 粗版 6 谱系
CAL_EPI        <- "上皮"       # [已签 §11 第 2 项] 单类目
CAL_RULE       <- "argmax"     # [已签 §12.1] 判据 = 行内最大类 == 上皮，**零阈值**

ROOT    <- "/home/eto/luad_v2"
RES     <- file.path(ROOT, "results/08_spatial_deconv")
GENEORD <- file.path(RES, "spatial_cnv/gene_order_spatial_hg38.tsv")
SMOKE   <- file.path(RES, "spatial_cnv/smoke_P4")
VISIUM  <- file.path(ROOT, "data/visium_spatial")

## —————————————————————————————————————————————————————————————
## 二、工具
## —————————————————————————————————————————————————————————————
args         <- commandArgs(trailingOnly = TRUE)
object_only  <- "--object-only" %in% args
keep_interm  <- "--keep" %in% args
patient      <- { i <- which(args == "--patient"); if (length(i)) args[i + 1] else "P4" }

dir.create(SMOKE, showWarnings = FALSE, recursive = TRUE)
options(scipen = 100)   # inferCNV 已知问题：不设这行，subclusters 模式下 hclust 会报错

t_start <- Sys.time()
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

## 峰值常驻内存（VmHWM = high-water mark），进程级、比 gc() 可靠
peak_gb <- function() {
  s <- readLines("/proc/self/status")
  hw <- grep("^VmHWM:", s, value = TRUE)
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", hw)) / 1024^2
}
step <- function(fmt, ...) cat(sprintf(
  paste0("[%s | %6.1f min | 峰值 %6.2f GB] ", fmt),
  format(Sys.time(), "%H:%M:%S"), elapsed(t_start), peak_gb(), ...), "\n")

## gzfile 连接必须显式关，否则 R 会攒一堆未关闭连接并告警
read_gz <- function(path, fn) {
  con <- gzfile(path)
  on.exit(try(close(con), silent = TRUE))
  fn(con)
}

step("==== 冒烟开始：患者 %s，参考口径 %s，判据 %s ===", patient, CAL_REFERENCE, CAL_EPI)

## —————————————————————————————————————————————————————————————
## 三、读 spot 掩码 与 RCTD 权重（口径：argmax）
## —————————————————————————————————————————————————————————————
sm <- read.delim(gzfile(file.path(RES, "spot_mask.tsv.gz")), stringsAsFactors = FALSE)
sm <- sm[as.character(sm$pass) %in% c("TRUE", "true", "1"), ]
sm <- sm[grepl(paste0("_", patient, "_"), sm$slide), ]
stopifnot(nrow(sm) > 0)
slides <- sort(unique(sm$slide))
step("患者 %s 的切片 %d 张：%s", patient, length(slides), paste(slides, collapse = ", "))

epi_bc  <- list()   # 逐张：argmax == 上皮 的 barcode
ref_bc  <- list()   # 逐张：良性切片（Normal/AAH*）的 argmax == 上皮，作 C 主锚
stage_of <- function(s) sub("^[^_]*_[^_]*_", "", s)

for (s in slides) {
  wf <- file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s))
  stopifnot(file.exists(wf))
  w <- read.delim(gzfile(wf), row.names = 1, check.names = FALSE)
  stopifnot(CAL_EPI %in% colnames(w), ncol(w) == 6L)
  keep <- intersect(rownames(w), sm$barcode[sm$slide == s])
  w    <- w[keep, , drop = FALSE]
  ## 🔴 §12.1：判据是**原始权重列之间的大小关系**，不做逐 spot 归一化（§M5/G4 禁令）
  ## 并列最大 ⇒ 判为未达标（剔除并计数），不随机、不按名排序打破
  W    <- as.matrix(w)
  mx   <- do.call(pmax, as.data.frame(W))   # 逐行 max（pmax 是跨列的向量化）
  n_at_max <- rowSums(W >= mx)              # matrix >= 向量在 R 里是按列内逐行回收 ⇒ 逐行比较正确
  am   <- colnames(W)[max.col(W, ties.method = "first")]
  tie  <- n_at_max > 1L
  epi  <- keep[am == CAL_EPI & !tie]
  st   <- stage_of(s)
  benign <- grepl("^Normal|^AAH", st)
  step("  %-26s QCpass %6d  上皮 %6d  并列剔除 %d  (%s%s)",
       s, nrow(w), length(epi), sum(tie), st, if (benign) " · 良性⇒可作锚" else "")
  epi_bc[[s]] <- epi
  if (benign) ref_bc[[s]] <- epi
}

n_epi <- sum(lengths(epi_bc))
n_ref <- sum(lengths(ref_bc))
step("合计：上皮 spot %d；其中良性切片 %d 可作 C 主锚；观测 %d（= 上皮 − 锚）",
     n_epi, n_ref, n_epi - n_ref)
if (n_ref == 0) stop(sprintf("患者 %s 无良性切片 ⇒ C 主锚必须借队列良性端（§11.2）；本次冒烟未实现借锚，停下", patient))

## —————————————————————————————————————————————————————————————
## 四、读表达矩阵（逐切片），只取 Gene Expression + 双探针按 symbol 求和
## —————————————————————————————————————————————————————————————
t_read <- Sys.time()
mats <- list()
for (s in slides) {
  d   <- file.path(VISIUM, s, "filtered_feature_bc_matrix")
  ft_all <- read_gz(file.path(d, "features.tsv.gz"), function(con)
    read.delim(con, header = FALSE, col.names = c("id", "sym", "type"),
               stringsAsFactors = FALSE))
  ## 🔴 matrix.mtx 里**含**抗体行（45 张切片有 35 行）。必须先在**全表**上核对维度，
  ##    再连同矩阵一起剔掉抗体行 —— 若先过滤 ft 再对维度，必然对不上。
  ab <- which(ft_all$type == "Antibody Capture")
  stopifnot(length(ab) %in% c(0L, 35L))
  bc  <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
  m   <- readMM(file.path(d, "matrix.mtx.gz"))   # readMM 自己认得 .gz 路径，别包 gzfile 连接
  stopifnot(nrow(m) == nrow(ft_all), ncol(m) == length(bc))
  ## ← 35 行抗体在此剔掉（SC0）。
  ## ⚠️ 陷阱：`x[-integer(0), ]` 在 R 里返回**零行**（不是全部行）—— 11 张切片没有抗体行，
  ##    不加这个 if 就会把它们整张清空，且不报错。
  if (length(ab)) {
    m  <- m[-ab, , drop = FALSE]
    ft <- ft_all[-ab, , drop = FALSE]
  } else {
    ft <- ft_all
  }
  stopifnot(nrow(ft) == 18085L, all(ft$type == "Gene Expression"))
  ## 只留本次要用的 spot：上皮 ∪ 锚（都是上皮，锚 ⊆ 上皮）
  use  <- match(c(epi_bc[[s]]), bc, nomatch = 0L)
  stopifnot(!any(use == 0L), !anyDuplicated(use))
  m    <- m[, use, drop = FALSE]
  bc   <- bc[use]

  ## 🔴 §2.3：HSPA14 / TBCE / TMSB15B 是双探针，必须**按 symbol 求和取一次**。
  ##    不求和而直接取名，会让其中一个探针在 inferCNV 里按名 join 时静默消失。
  g    <- factor(ft$sym, levels = unique(ft$sym))
  Agg  <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1,
                       dims = c(nlevels(g), length(g)))
  rownames(Agg) <- levels(g)
  m    <- Agg %*% m
  rownames(m) <- levels(g)
  colnames(m) <- paste(s, bc, sep = "|")            # 🔴 barcode 跨切片重复，必须带切片前缀
  nd <- sum(duplicated(ft$sym))
  if (nd) step("  %-26s 双探针基因 %d 个已按 symbol 求和", s, nd)
  mats[[s]] <- m
  rm(m, Agg); invisible(gc())
}
mat <- do.call(cbind, mats)
rm(mats); invisible(gc())
step("读入并合并：%d 基因 x %d spot，nnz %s（用时 %.1f min）",
     nrow(mat), ncol(mat), format(nnz <- length(mat@x), big.mark = ","), elapsed(t_read))

## —————————————————————————————————————————————————————————————
## 五、注解：观测 / 参考（C 主锚）
## —————————————————————————————————————————————————————————————
cells_obs <- unlist(lapply(slides, function(s) paste(s, epi_bc[[s]], sep = "|")),
                    use.names = FALSE)
cells_ref <- unlist(lapply(names(ref_bc), function(s) paste(s, ref_bc[[s]], sep = "|")),
                    use.names = FALSE)
cells_ref <- setdiff(cells_ref, character(0))
stopifnot(all(cells_ref %in% cells_obs))         # C 主锚 ⊆ 上皮集合
cells_obs <- setdiff(cells_obs, cells_ref)       # inferCNV 要求观测与参考**互斥**
stopifnot(all(c(cells_obs, cells_ref) %in% colnames(mat)))

ann <- data.frame(cell = c(cells_obs, cells_ref),
                  group = c(rep(paste0("obs_", patient), length(cells_obs)),
                            rep(P_REF_GROUP, length(cells_ref))),
                  stringsAsFactors = FALSE)
annf <- file.path(SMOKE, sprintf("annotations_%s.tsv", patient))
write.table(ann, annf, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
step("注解：观测 %d（%s）＋ 参考 %d（%s 主锚）；共 %d spot",
     length(cells_obs), paste0("obs_", patient), length(cells_ref), P_REF_GROUP, nrow(ann))

## —————————————————————————————————————————————————————————————
## 六、建 infercnv 对象
## —————————————————————————————————————————————————————————————
t_obj <- Sys.time()
obj <- CreateInfercnvObject(raw_counts_matrix = mat,
                            gene_order_file   = GENEORD,
                            annotations_file  = annf,
                            ref_group_names   = c(P_REF_GROUP),
                            chr_exclude       = P_CHR_EXCLUDE,   # [已签 §14.3] 照原文 c("chrM")，见上方 gene 账
                            delim             = "\t")
## 基因账逐项在运行时**算出来**（不是照抄注释），防止哪天悄悄变了
go_tab  <- read.delim(GENEORD, header = FALSE, stringsAsFactors = FALSE)
n_sym   <- nrow(mat)
n_ord   <- nrow(go_tab)
n_sexmt <- sum(go_tab[[2]] %in% P_CHR_EXCLUDE)
stopifnot(n_ord - n_sexmt == nrow(obj@expr.data))
step("基因账：矩阵 %d − 未匹配 %d = 位置表 %d − 按 chr_exclude 剔 %d（%s）= 进管线 %d",
     n_sym, n_sym - n_ord, n_ord, n_sexmt, paste(P_CHR_EXCLUDE, collapse = "/"),
     nrow(obj@expr.data))

## 🔴 必须转 base matrix（03_cnv/16 已实测：dgeMatrix 上 .subtract_expr 单行抽取
##    1.23 s/行 vs base matrix 0.002 s/行，差 ~600x；数值逐位相同，非口径变更）
obj@expr.data <- as.matrix(obj@expr.data)
stopifnot(is.matrix(obj@expr.data), !isS4(obj@expr.data))

step("对象建成：%d 基因 x %d spot（expr.data 类 %s）用时 %.1f min",
     nrow(obj@expr.data), ncol(obj@expr.data),
     paste(class(obj@expr.data), collapse = "/"), elapsed(t_obj))
step("基因保留 %d / 18015（位置表匹配数）；参考组 %s；观测组 %s",
     nrow(obj@expr.data),
     paste(names(obj@reference_grouped_cell_indices), collapse = ","),
     paste(names(obj@observation_grouped_cell_indices), collapse = ","))

obj_gb <- as.numeric(object.size(obj@expr.data)) / 1024^3
step("expr.data 单一副本 %.2f GB；进程峰值 %.2f GB", obj_gb, peak_gb())

manifest <- list(
  script = "08_spatial_deconv/11_run_spatial_cnv_smoke.R",
  purpose = "SPATIAL_CNV_PREREG.md §9 第 2 项：单患者内存峰值与墙钟耗时（成本冒烟，非结论）",
  patient = patient, slides = slides,
  n_spot_epi = n_epi, n_spot_ref = n_ref,
  n_obs = length(cells_obs), n_ref = length(cells_ref),
  n_genes_matrix = nrow(mat), n_genes_infercnv = nrow(obj@expr.data),
  nnz = nnz, expr_data_gb = round(obj_gb, 3),
  params_signed = list(cutoff = P_CUTOFF, cluster_by_groups = P_CLUSTER_BY_GROUPS,
                       HMM = P_HMM, denoise = P_DENOISE, ref_group = P_REF_GROUP,
                       threads = P_THREADS),
  params_signed_later = list(analysis_mode = P_ANALYSIS_MODE,   # [已签 §13.1] samples，不取默认
                             no_plot = P_NO_PLOT,               # [已签 §13.1]
                             chr_exclude = P_CHR_EXCLUDE),      # [已签 §14.3] 原文代码逐字
  params_unregistered_default = list(window_length = P_WINDOW_LENGTH,
                                     hclust_method = P_HCLUST_METHOD),
  ## §14.3 的守卫：chrX/chrY 已放回管线 ⇒ 必须记下参考与观测**是否同患者**。
  ## 本脚本的锚只从**本患者**良性切片取（上面 for 循环里 ref_bc 只被本患者切片填充）
  ## ⇒ 结构上必然同患者、性别恒定、chrX/Y 无性别混淆。
  ## ⚠️ 借锚路径（§11.2 跨患者借队列良性端）**尚未实现**，见预注册 §14.6。
  chr_xy = list(kept = setdiff(c("chrX", "chrY"), P_CHR_EXCLUDE),
                dropped = P_CHR_EXCLUDE,
                n_dropped_from_gene_order = n_sexmt,
                ref_same_patient = TRUE,
                ref_slides = names(ref_bc),
                obs_slides = slides,
                rule = "跨患者借锚时 chrX/Y 结论降级待定；本脚本只做同患者锚"),
  gene_account = list(n_symbols_matrix = n_sym, n_unmatched = n_sym - n_ord,
                      n_in_gene_order = n_ord, n_dropped_sex_mt = n_sexmt,
                      n_into_infercnv = nrow(obj@expr.data)),
  caliber = list(reference = CAL_REFERENCE, epi = CAL_EPI, rule = CAL_RULE),
  peak_gb_after_object = round(peak_gb(), 3),
  elapsed_min_read = round(elapsed(t_read), 2),
  elapsed_min_object = round(elapsed(t_obj), 2)
)

if (object_only) {
  manifest$stopped_at <- "object-only"
  manifest$elapsed_min_total <- round(elapsed(t_start), 2)
  jsonlite::write_json(manifest, file.path(SMOKE, sprintf("smoke_%s_manifest.json", patient)),
                       auto_unbox = TRUE, pretty = TRUE)
  step("--object-only：停在这里。总用时 %.1f min", elapsed(t_start))
  quit(save = "no", status = 0)
}

## —————————————————————————————————————————————————————————————
## 七、跑管线（HMM=TRUE）—— SC1 要求日志里必须看得见 rjags 真被调用
## —————————————————————————————————————————————————————————————
OUT <- file.path(SMOKE, sprintf("infercnv_out_%s", patient))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
step("infercnv::run 开始（HMM=%s，analysis_mode=%s，threads=%d）", P_HMM, P_ANALYSIS_MODE, P_THREADS)
step("==== 下面会打印 inferCNV 自己的分步日志；每步耗时见其 STEP 编号 ====")

t_run <- Sys.time()
obj <- infercnv::run(obj,
                     out_dir           = OUT,
                     cutoff            = P_CUTOFF,
                     window_length     = P_WINDOW_LENGTH,
                     denoise           = P_DENOISE,
                     HMM               = P_HMM,
                     analysis_mode     = P_ANALYSIS_MODE,
                     cluster_by_groups = P_CLUSTER_BY_GROUPS,
                     hclust_method     = P_HCLUST_METHOD,
                     write_expr_matrix = FALSE,
                     no_plot           = P_NO_PLOT,
                     num_threads       = P_THREADS,
                     plot_steps        = FALSE,
                     save_rds          = keep_interm,
                     save_final_rds    = keep_interm)

manifest$elapsed_min_run <- round(elapsed(t_run), 2)
manifest$elapsed_min_total <- round(elapsed(t_start), 2)
manifest$peak_gb_final <- round(peak_gb(), 3)
manifest$out_dir <- OUT
jsonlite::write_json(manifest, file.path(SMOKE, sprintf("smoke_%s_manifest.json", patient)),
                     auto_unbox = TRUE, pretty = TRUE)

step("==== 冒烟完成：建对象 %.1f min ＋ 跑管线 %.1f min ＝ 总 %.1f min；峰值内存 %.2f GB ====",
     elapsed(t_obj), elapsed(t_run), elapsed(t_start), peak_gb())
step("产物目录 %s", OUT)
