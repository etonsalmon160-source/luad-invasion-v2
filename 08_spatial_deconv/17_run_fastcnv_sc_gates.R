#!/usr/bin/env Rscript
# 空转（Visium）CNV 臂 —— **门 SC0/SC1/SC2 首跑**（SPATIAL_CNV_PREREG.md §18.9，2026-09-29 已签）
#
# 本脚本按 §13.3 第 2 项执行：**SC2 第一炮必须是单患者**。
#   引擎 = fastCNV v1.1.11（§18.1）；参数沿用 §15.4。
#   参考 = 该患者自己的良性上皮 spot（同患者 C 主锚，§11.2）；
#   观测 = 该患者的 LUAD 上皮 spot。
#
# 三道门（逐条对照 §18.9）：
#   SC0 底座核对  —— ①基因序匹配率 ≥99% 且**重算哈希**；②GE 行 18085 / 抗体行 35 已剔；
#                    ③spot 数 == spot_mask 的 n_pass；④fastCNV 自己的基因账（**只报不停**）
#   SC1 参数冻结  —— §15.4 参数逐位相同；**必须出现「未聚合、逐 spot」的证据**
#   SC2 锚可复算  —— 同一参考同一输入，**重复 N≥3 次**，逐 spot `cnv_fraction`
#                    **≥95% 逐位相同**；每次锚名单逐字节相同
#
# 🔴 本脚本**不产生恶性判定**（阈值口径 §18.5(a) 已签，但那是 SC3/SC4 的事）；
#    它只回答「门过不过」。SC3/SC4 **本次不跑**。
#
# 跑法：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 08_spatial_deconv/17_run_fastcnv_sc_gates.R [--patient P4]

suppressMessages({
  library(fastCNV); library(Seurat); library(Matrix)
})

## —————————————————————————————————————————————————————————————
## 一、参数（逐项显式，法则 3.1；来源 §15.4 / §18.5 / §18.9）
## —————————————————————————————————————————————————————————————
CAL_EPI       <- "上皮"      # [已签 §12.1] 单类目
CAL_RULE      <- "argmax"    # [已签 §12.1] 行内最大类 == 上皮，零阈值
F_WINDOW_SIZE <- 150;  F_WINDOW_STEP <- 10
F_TOPN_GENES  <- 7000; F_THRESH_PCT  <- 0.01
F_GET_ARM     <- TRUE;  F_ASSAY <- "Spatial"
N_REPEAT      <- 3L          # §18.9 SC2
SC2_MIN_IDENT <- 0.95        # §18.9 SC2

## SC0 期望值（§4 / §18.9① ②）—— 取自 `10_build_spatial_gene_order.py` 的 manifest 口径
EXP_GE_ROWS    <- 18085L        # GE 行数**恒为** 18085（实测 56 张全同）
EXP_AB_SET     <- c(0L, 35L)    # 🔴 抗体行**并非每张都有**：实测 **45/56** 有 35 行、**11/56** 为 0
                                #    （首跑把 35 当成每张都有的期望值 ⇒ SC0② 假报警，已改）
EXP_UNIQUE_SYM <- 18082L
SC0_MIN_MATCH  <- 0.99       # §18.9 SC0①

ROOT   <- "/home/eto/luad_v2"
RES    <- file.path(ROOT, "results/08_spatial_deconv")
VISIUM <- file.path(ROOT, "data/visium_spatial")
GODIR  <- file.path(RES, "spatial_cnv")
OUTD   <- file.path(RES, "spatial_cnv/sc_gates_P4")

## —————————————————————————————————————————————————————————————
## 二、工具函数
## —————————————————————————————————————————————————————————————
args    <- commandArgs(trailingOnly = TRUE)
patient <- { i <- which(args == "--patient"); if (length(i)) args[i + 1] else "P4" }

dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
peak_gb <- function() {
  s <- readLines("/proc/self/status")
  hw <- grep("^VmHWM:", s, value = TRUE)
  as.numeric(sub(".*?([0-9]+) kB.*", "\\1", hw)) / 1024^2
}
step <- function(fmt, ...) cat(sprintf(
  paste0("[%s | %6.1f min | 峰值 %6.2f GB] ", fmt),
  format(Sys.time(), "%H:%M:%S"), elapsed(t_start), peak_gb(), ...), "\n")
read_gz  <- function(p, fn) { con <- gzfile(p); on.exit(try(close(con), silent = TRUE)); fn(con) }
stage_of <- function(s) sub("^[^_]*_[^_]*_", "", s)
is_benign <- function(s) grepl("^Normal|^AAH", stage_of(s))
## 用系统 sha256sum，保证与 `10_*.py`（Python hashlib）算出的哈希**逐位可比**
sha256 <- function(p) {
  out <- system2("sha256sum", shQuote(p), stdout = TRUE)
  if (!length(out)) stop("sha256sum 失败：", p)
  sub(" .*$", "", out[1])
}

FAILS <- character(0)
gate  <- function(ok, label, detail) {
  step("  [%s] %s —— %s", if (isTRUE(ok)) "PASS" else "FAIL", label, detail)
  if (!isTRUE(ok)) FAILS <<- c(FAILS, label)
  isTRUE(ok)
}

step("==== SC 门首跑（单患者）开始：%s ====", patient)
step("引擎 fastCNV %s / Seurat %s（Genome Medicine 2026-08-29 已发表，§15.2）",
     as.character(packageVersion("fastCNV")), as.character(packageVersion("Seurat")))

## —————————————————————————————————————————————————————————————
## 三、切片、掩码、上皮名单（口径与 11_/12_ 逐字同源）
## —————————————————————————————————————————————————————————————
sm <- read.delim(gzfile(file.path(RES, "spot_mask.tsv.gz")), stringsAsFactors = FALSE)
sm <- sm[as.character(sm$pass) %in% c("TRUE", "true", "1"), ]
sm <- sm[grepl(paste0("_", patient, "_"), sm$slide), ]
stopifnot(nrow(sm) > 0)
slides <- sort(unique(sm$slide))
step("患者 %s 的切片 %d 张：%s", patient, length(slides), paste(slides, collapse = ", "))

epi_bc <- list(); n_pass <- list()
for (s in slides) {
  wf <- file.path(RES, sprintf("rctd_a/per_slide/%s.weights.tsv.gz", s))
  stopifnot(file.exists(wf))
  w <- read.delim(gzfile(wf), row.names = 1, check.names = FALSE)
  stopifnot(CAL_EPI %in% colnames(w), ncol(w) == 6L)
  keep <- intersect(rownames(w), sm$barcode[sm$slide == s])
  W    <- as.matrix(w[keep, , drop = FALSE])
  mx   <- do.call(pmax, as.data.frame(W)); n_at <- rowSums(W >= mx)
  am   <- colnames(W)[max.col(W, ties.method = "first")]
  epi_bc[[s]] <- keep[am == CAL_EPI & n_at <= 1L]
  n_pass[[s]] <- sum(sm$slide == s)
  step("  %-26s 掩码 pass %6d  上皮 %6d  (%s%s)", s, n_pass[[s]], length(epi_bc[[s]]),
       stage_of(s), if (is_benign(s)) " · 良性" else "")
}
step("上皮 spot 合计 %d", sum(lengths(epi_bc)))

## —————————————————————————————————————————————————————————————
## 四、读表达矩阵；建**一个** Seurat 对象
## —————————————————————————————————————————————————————————————
mats <- list(); lab_a <- character(0); cell_slide <- character(0)
sc0_ge <- list(); sc0_ab <- list(); sc0_bc <- list(); sc0_other <- list()
for (s in slides) {
  d      <- file.path(VISIUM, s, "filtered_feature_bc_matrix")
  ft_all <- read_gz(file.path(d, "features.tsv.gz"), function(con)
    read.delim(con, header = FALSE, col.names = c("id", "sym", "type"), stringsAsFactors = FALSE))
  ab <- which(ft_all$type == "Antibody Capture")
  bc <- read_gz(file.path(d, "barcodes.tsv.gz"), readLines)
  m  <- readMM(file.path(d, "matrix.mtx.gz"))
  stopifnot(nrow(m) == nrow(ft_all), ncol(m) == length(bc))
  ## 🔴 `readMM()` 不设 dimnames（.mtx 里没有）⇒ 必须显式补列名（否则逐张 0 匹配）
  colnames(m) <- bc
  if (length(ab)) { m <- m[-ab, , drop = FALSE]; ft <- ft_all[-ab, , drop = FALSE] } else { ft <- ft_all }
  ## 🔴 必须数 **Gene Expression 行**，不是总行数：总行数 = GE + 抗体，
  ##    45/56 张总行数是 18,120（首跑把总行数当 GE 行数 ⇒ 那 45 张全判错，已改）
  sc0_ge[[s]] <- sum(ft_all$type == "Gene Expression")
  sc0_other[[s]] <- setdiff(unique(ft_all$type), c("Gene Expression", "Antibody Capture"))
  sc0_ab[[s]] <- length(ab); sc0_bc[[s]] <- length(bc)

  use <- match(epi_bc[[s]], bc, nomatch = 0L)
  stopifnot(!any(use == 0L), !anyDuplicated(use))

  ## SC0③ 掩码 pass 的 spot 必须**全都在**表达矩阵里（掩码是表达矩阵的子集）
  miss_mask <- setdiff(sm$barcode[sm$slide == s], bc)
  if (length(miss_mask)) step("  🔴 %s：掩码 pass 里有 %d 个 barcode 不在表达矩阵", s, length(miss_mask))

  m  <- m[, use, drop = FALSE]; bt <- bc[use]
  g   <- factor(ft$sym, levels = unique(ft$sym))
  Agg <- sparseMatrix(i = as.integer(g), j = seq_along(g), x = 1, dims = c(nlevels(g), length(g)))
  m <- Agg %*% m; rownames(m) <- levels(g)
  colnames(m) <- paste(s, bt, sep = "_")

  mats[[s]]  <- m
  cell_slide <- c(cell_slide, setNames(rep(s, ncol(m)), colnames(m)))
  lab_a      <- c(lab_a, rep(if (is_benign(s)) "benign_epi" else "tumor_epi", ncol(m)))
}
mat <- do.call(cbind, mats); rm(mats); invisible(gc())
stopifnot(ncol(mat) == length(lab_a), identical(colnames(mat), names(cell_slide)))
obj <- CreateSeuratObject(counts = mat, assay = F_ASSAY,
                          meta.data = data.frame(g_a = lab_a, slide = unname(cell_slide),
                                                 row.names = colnames(mat)))
rm(mat); invisible(gc())
step("对象就绪：%d 基因 x %d spot（峰值 %.2f GB）", nrow(obj), ncol(obj), peak_gb())

## —————————————————————————————————————————————————————————————
## 五、SC0 底座核对
## —————————————————————————————————————————————————————————————
step("---- SC0 底座核对 ----")

## ① 基因序匹配率 ≥99% ＋ **重算哈希**
go_path <- file.path(GODIR, "gene_order_spatial_hg38.tsv")
go_man  <- jsonlite::fromJSON(file.path(GODIR, "gene_order_spatial_manifest.json"))
go      <- read.delim(go_path, header = FALSE, stringsAsFactors = FALSE)
go_sym  <- go[[1]]
h_recomputed <- sha256(go_path)
gate(identical(h_recomputed, go_man$out_sha256), "SC0① 基因序哈希重算",
     sprintf("%s（清单 %s）", substr(h_recomputed, 1, 16), substr(go_man$out_sha256, 1, 16)))
our_sym <- rownames(obj)
mr <- mean(our_sym %in% go_sym)
gate(mr >= SC0_MIN_MATCH, "SC0① 基因序匹配率",
     sprintf("%.4f（%d/%d），门槛 ≥%.2f", mr, sum(our_sym %in% go_sym), length(our_sym), SC0_MIN_MATCH))
step("  ⚠️ 说明：这份基因序是 §3 infercnv 臂的底座；**新引擎 fastCNV 不消费它**（它用自带 getGenes()）")
step("     ⇒ SC0① 在新引擎下是**溯源一致性核对**，真正承重的是下面 SC0④ 的 fastCNV 基因账")

## ② GE 行数恒为 18085；抗体行 ∈ {0,35} 并**逐张报**（非通用 35，见上）
ok_ge <- all(unlist(sc0_ge) == EXP_GE_ROWS); ok_ab <- all(unlist(sc0_ab) %in% EXP_AB_SET)
gate(ok_ge, "SC0② GE 行数", sprintf("逐张 %s（期望 %d）", paste(unique(unlist(sc0_ge)), collapse = ","), EXP_GE_ROWS))
gate(ok_ab, "SC0② 抗体行 ⊆ {0,35} 且已剔",
     sprintf("逐张 %s（本患者无抗体行的切片 = %s）",
             paste(sprintf("%s:%d", names(unlist(sc0_ab)), unlist(sc0_ab)), collapse = " "),
             paste(names(unlist(sc0_ab))[unlist(sc0_ab) == 0], collapse = ",")))
## ②′ 除了 GE 与抗体，不得有第三种 type（有 ⇒ 我们的行数账不完整）
oth <- sort(unique(unlist(sc0_other)))
gate(!length(oth), "SC0② 无第三类 feature",
     if (length(oth)) sprintf("🔴 出现未预期 type：%s", paste(oth, collapse = ",")) else "只有 Gene Expression 与 Antibody Capture")

## ③ spot 数 == spot_mask 的 n_pass
ok_spot <- TRUE
for (s in slides) {
  n_epi_s <- length(epi_bc[[s]])
  step("    %-26s 掩码 pass %6d ⇒ 上皮取用 %6d", s, n_pass[[s]], n_epi_s)
  if (n_epi_s > n_pass[[s]]) ok_spot <- FALSE
}
gate(ok_spot, "SC0③ 上皮 spot ⊆ 掩码 pass", sprintf("合计 pass %d ≥ 上皮 %d",
     sum(unlist(n_pass)), sum(lengths(epi_bc))))

## ④ fastCNV 自己的基因账（**只报不停**）
gm      <- getGenes()
gm_f    <- gm[gm$gene_biotype %in% c("protein_coding", "lncRNA") &
              gm$chromosome_name %in% c(1:22, "X") & gm$hgnc_symbol != "", ]
tok     <- unique(gm_f$hgnc_symbol); in_tab <- our_sym %in% tok
gene_acct <- list(n_our_symbols = length(our_sym), n_matched = sum(in_tab),
                  n_unmatched = sum(!in_tab), unmatched_examples = head(our_sym[!in_tab], 15),
                  dropped_by_tool_rule = "chrY 与 chrM 被工具硬编码丢弃（§15.8 第 1 项 / §14.3 二次改签）")
step("SC0④ fastCNV 基因账：我方 %d − 未匹配 %d = 进管线候选 %d（**只报不停**）",
     gene_acct$n_our_symbols, gene_acct$n_unmatched, gene_acct$n_matched)
if (length(gene_acct$unmatched_examples))
  step("    未匹配示例：%s", paste(gene_acct$unmatched_examples, collapse = ", "))

## —————————————————————————————————————————————————————————————
## 六、SC1 参数已冻结且**没静默走默认聚合**
## —————————————————————————————————————————————————————————————
step("---- SC1 参数冻结 ----")
p_frozen <- c(assay = F_ASSAY, windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
              topNGenes = F_TOPN_GENES, thresholdPercentile = F_THRESH_PCT,
              getCNVPerChromosomeArm = F_GET_ARM)
step("  显式参数：%s", paste(sprintf("%s=%s", names(p_frozen), p_frozen), collapse = "  "))
step("  🔴 本脚本走 **CNVCalling 单对象路径**，`prepareCounts`/`aggregFactor` 是 `fastCNV()` 层参数、")
step("     在此路径下**不参与** ⇒ 不存在「spot<60,000 被聚合成 meta spot」的默认（§15.4 / §15.8 第 6 项）")
## 证据：逐切片 spot 数 == 掩码 pass 下的上皮数（若被聚合，会掉到约 1/6）
n_obj_per_slide <- table(obj$slide)
agg_ok <- all(as.integer(n_obj_per_slide) == unlist(lapply(slides, function(s) length(epi_bc[[s]]))))
gate(agg_ok, "SC1 未聚合证据（逐 spot 分辨率保住）",
     sprintf("对象逐切片 spot：%s", paste(sprintf("%s=%d", names(n_obj_per_slide), n_obj_per_slide), collapse = " ")))

## —————————————————————————————————————————————————————————————
## 七、SC2 锚可复算（重复 N≥3 次）
## —————————————————————————————————————————————————————————————
step("---- SC2 锚可复算：同一参考同一输入，重复 %d 次 ----", N_REPEAT)
ref_var <- "g_a"; ref_lab <- "benign_epi"
ref_bc <- sort(colnames(obj)[obj$g_a == ref_lab])
step("  参考 spot %d（%s，同患者 C 主锚 §11.2）；观测 spot %d",
     length(ref_bc), ref_lab, sum(obj$g_a != ref_lab))

cf_list <- list(); anchor_same <- TRUE; run_min <- numeric(0)
for (k in seq_len(N_REPEAT)) {
  t0 <- Sys.time()
  o  <- CNVCalling(obj, assay = F_ASSAY, referenceVar = ref_var, referenceLabel = ref_lab,
                   scaleOnReferenceLabel = TRUE, thresholdPercentile = F_THRESH_PCT,
                   geneMetadata = gm, windowSize = F_WINDOW_SIZE, windowStep = F_WINDOW_STEP,
                   saveGenomicWindows = FALSE, topNGenes = F_TOPN_GENES)
  if (F_GET_ARM) o <- CNVPerChromosomeArm(o)
  cf <- FetchData(o, vars = "cnv_fraction")[[1]]
  names(cf) <- colnames(o)
  cf_list[[k]] <- cf
  rb <- sort(colnames(o)[o$g_a == ref_lab])
  if (!identical(rb, ref_bc)) anchor_same <- FALSE
  run_min <- c(run_min, elapsed(t0))
  step("  第 %d/%d 次完成，用时 %.2f min；cnv_fraction 中位 %.4f（峰值 %.2f GB）",
       k, N_REPEAT, run_min[k], median(cf, na.rm = TRUE), peak_gb())
  rm(o); invisible(gc())
}

## 逐位相同的比例（**逐位**：完全相等；并列各算 0.5 不算过）
same <- rep(TRUE, length(cf_list[[1]]))
for (k in 2:N_REPEAT) same <- same & (cf_list[[1]] == cf_list[[k]])
frac_ident <- mean(same)
step("  逐位相同的 spot：%d / %d = %.6f", sum(same), length(same), frac_ident)
gate(anchor_same, "SC2 锚名单逐字节相同", sprintf("参考 %d spot，%d 次读数一致", length(ref_bc), N_REPEAT))
gate(frac_ident >= SC2_MIN_IDENT, "SC2 逐 spot cnv_fraction ≥95% 逐位相同",
     sprintf("%.6f（门槛 %.2f）", frac_ident, SC2_MIN_IDENT))
if (frac_ident < 1 && frac_ident >= SC2_MIN_IDENT) {
  d <- which(!same)
  step("  ⚠️ 有 %d 个 spot 不逐位相同，最大绝对差 %.3g",
       length(d), max(abs(cf_list[[1]][d] - cf_list[[2]][d])))
}

## —————————————————————————————————————————————————————————————
## 八、结论与 manifest
## —————————————————————————————————————————————————————————————
verdict <- if (length(FAILS)) "FAIL" else "PASS"
step("==== 门结论：%s%s ====", verdict,
     if (length(FAILS)) paste0("（未过：", paste(FAILS, collapse = " / "), "）") else "")
if (verdict == "FAIL") step("🔴 SC2 是本臂第一道硬门 ⇒ 门不过，本臂到此为止，不得继续 SC3/SC4")

manifest <- list(
  script = "08_spatial_deconv/17_run_fastcnv_sc_gates.R",
  purpose = "SPATIAL_CNV_PREREG.md §18.9：门 SC0/SC1/SC2 首跑（单患者，§13.3 第 2 项）",
  patient = patient, slides = slides, n_epi = sum(lengths(epi_bc)),
  tool = list(name = "fastCNV", version = as.character(packageVersion("fastCNV")), license = "GPL-3",
              source = "Cabrejas et al., Genome Medicine 2026, DOI 10.1186/s13073-026-01731-w"),
  params = as.list(p_frozen),
  reference = list(var = ref_var, label = ref_lab, n_ref = length(ref_bc),
                   caliber = "同患者良性上皮 C 主锚（§11.2）；argmax==上皮（§12.1）"),
  sc0 = list(hash_recomputed = h_recomputed, hash_in_manifest = go_man$out_sha256,
             hash_ok = identical(h_recomputed, go_man$out_sha256),
             gene_order_match_rate = round(mr, 4), gene_order_n = length(go_sym),
             ge_rows = unlist(sc0_ge), ab_rows = unlist(sc0_ab),
             mask_pass = unlist(n_pass), gene_account = gene_acct,
             note = "基因序是 §3 infercnv 臂的底座，新引擎不消费它 ⇒ SC0① 为溯源核对，非承重"),
  sc1 = list(no_aggregation_evidence = agg_ok, n_per_slide = as.list(n_obj_per_slide),
             note = "单对象 CNVCalling 路径不含 prepareCounts/aggregFactor ⇒ 不存在默认聚合"),
  sc2 = list(n_repeat = N_REPEAT, frac_bit_identical = frac_ident, min_ident = SC2_MIN_IDENT,
             anchor_names_identical = anchor_same, run_min = round(run_min, 3),
             note = "fastCNV 确定性；重复同一输入。⚠️ 本检验强度有限：重复的是**同一条**计算路径，
                    不覆盖「换一批锚 spot」的稳定性 —— 那属 §18.9 SC2 未写明的范围，不得据本条宣称锚稳"),
  verdict = verdict, fails = FAILS,
  not_claimed = c("不产生恶性判定（阈值口径 §18.5(a) 已签，但属 SC3/SC4）",
                  "未跑 SC3/SC4",
                  "SC2 只测确定性与锚名单一致性；不测「换锚 spot 后是否仍稳」",
                  "参考组 cnv_fraction 由构造接近 0，不是证据（§15.3 循环性）"),
  elapsed_min_total = round(elapsed(t_start), 2), peak_gb = round(peak_gb(), 3))

jsonlite::write_json(manifest, file.path(OUTD, "sc_gates_P4_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)
step("产物：%s", OUTD)
step("==== 完成：总用时 %.1f min；峰值 %.2f GB；结论 %s ====", elapsed(t_start), peak_gb(), verdict)
