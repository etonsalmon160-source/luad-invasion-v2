#!/usr/bin/env Rscript
# inferCNV 冒烟 v2：按染色体聚合的统计量。
# 判读**写死在 INFERCNV_SMOKE_PREREG_v2_CHRSTAT.md §五**，本脚本只执行。
#
# 跑法：
#   R_LIBS=/home/eto/Rlibs/infercnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript --vanilla 03_cnv/18_analyze_infercnv_chrstat.R

suppressMessages(library(data.table))

ROOT <- "/home/eto/luad_v2"
D <- file.path(ROOT, "results/03_cnv/infercnv_smoke")
FIGDIR <- file.path(D, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ——— 预注册 v2 冻结常量（§二 / §四 / §五）———
CHRS <- paste0("chr", 1:22)
REF_GROUP <- "Normal"
HOLDOUT_GROUP <- "Normal_holdout"
POS_GROUP <- "P4_LUAD"
N_GENES_EXPECTED <- 7668L
N_CELLS_EXPECTED <- 19748L
IQR_MIN <- 0            # 严格 > 0（§五 I1）
REF_MEDIAN_MAX <- 0.05  # §五 I2
AUROC_STRONG <- 0.80
AUROC_WEAK <- 0.60
EXPECTED_CELLS_PER_GROUP <- c(Normal = 7255L, Normal_holdout = 1280L,
                              P4_AAH = 1009L, P4_AAH1 = 2702L, P4_LUAD = 7502L)

auroc <- function(pos, neg) {
  y <- c(rep(1, length(pos)), rep(0, length(neg)))
  s <- c(pos, neg)
  r <- rank(s, ties.method = "average")
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cat("[in] 读残差矩阵 ...\n")
e <- readRDS(file.path(D, "p4_final_expr_data.rds"))
stopifnot(is.matrix(e), nrow(e) == N_GENES_EXPECTED, ncol(e) == N_CELLS_EXPECTED)
cat(sprintf("[in] %d 基因 x %d 细胞\n", nrow(e), ncol(e)))

cat("[in] 读基因位置表 ...\n")
go <- fread(file.path(D, "gene_order_hg38.tsv"), header = FALSE, sep = "\t",
            col.names = c("gene", "chr", "start", "end"), data.table = FALSE)
mi <- match(rownames(e), go$gene)
unmatched <- rownames(e)[is.na(mi)]
# 法则 0：映射不上的基因逐个列出，绝不静默丢
writeLines(unmatched, file.path(D, "v2_unmatched_genes.txt"))
cat(sprintf("[gene] 映射到位置表 %d / %d；未映射 %d（清单 → v2_unmatched_genes.txt）\n",
            sum(!is.na(mi)), nrow(e), length(unmatched)))
chr <- go$chr[mi]
off_chr <- unique(chr[!is.na(chr) & !chr %in% CHRS])
if (length(off_chr)) cat(sprintf("[gene] 落在 chr1..chr22 之外的染色体：%s\n", paste(off_chr, collapse = ", ")))

cat("[in] 读分组 ...\n")
# ⚠️ 必须显式 sep="\t"：barcode 带 `|`，fread 自动猜分隔符会猜错
ann <- fread(file.path(D, "p4_annotations.tsv"), header = FALSE, sep = "\t",
             col.names = c("cell", "group"), data.table = FALSE)
stopifnot(!anyDuplicated(ann$cell), setequal(ann$cell, colnames(e)))
grp <- setNames(ann$group, ann$cell)
g <- grp[colnames(e)]
tb <- table(g)
print(tb)
stopifnot(identical(as.integer(tb[names(EXPECTED_CELLS_PER_GROUP)]),
                    as.integer(EXPECTED_CELLS_PER_GROUP)))

# ——— §四 统计量（冻结）———
keep <- !is.na(chr) & chr %in% CHRS
cat(sprintf("[stat] 进入统计的基因 %d（chr1..chr22）\n", sum(keep)))
lg <- log2(e)
lg[!is.finite(lg)] <- NA
# 每条染色体对**细胞**求均值（列均值），得到 细胞 x 染色体
M <- sapply(CHRS, function(k) colMeans(lg[keep & chr == k, , drop = FALSE], na.rm = TRUE))
stopifnot(nrow(M) == ncol(e), ncol(M) == length(CHRS))
rownames(M) <- colnames(e); colnames(M) <- CHRS
# 主统计量：22 条染色体上取极差（对"每细胞整体减常数"不变）
score_range  <- apply(M, 1, function(v) max(v) - min(v))
score_maxabs <- apply(M, 1, function(v) max(abs(v)))

cat("\n[stat] 各组 score_range（主统计量，log2）：\n")
st <- do.call(rbind, lapply(split(score_range, g), function(v) {
  data.frame(n = length(v), median = median(v), mean = mean(v),
             q25 = quantile(v, .25, names = FALSE), q75 = quantile(v, .75, names = FALSE),
             p90 = quantile(v, .90, names = FALSE), max = max(v))
}))
st <- cbind(group = rownames(st), st); print(st, row.names = FALSE)
cat("\n[stat] 各组 score_maxabs（次统计量，仅描述）：\n")
print(round(tapply(score_maxabs, g, median), 5))

# ——— 门槛 I1：仪器有效 ———
iq <- IQR(score_range)
i1 <- iq > IQR_MIN
# ——— 门槛 I2：参考组合理性 ———
ref_med <- median(score_range[g == REF_GROUP])
i2 <- ref_med <= REF_MEDIAN_MAX
# ——— 门槛 I3：有没有信号 ———
pos <- score_range[g == POS_GROUP]
neg <- score_range[g == HOLDOUT_GROUP]
a <- auroc(pos, neg)
v3 <- if (a >= AUROC_STRONG) {
  "有信号，值得进下一轮（仍须另行签字）"
} else if (a >= AUROC_WEAK) {
  "弱信号"
} else {
  "分不开"
}

if (!i1) {
  cat("\n[I1] 统计量退化（IQR = 0）—— 本次作废，禁止记成「分不开」\n")
} else {
  cat(sprintf("\n[I1] 全体细胞 score_range IQR = %.6f  (>0)  → 过；取值数 %d / %d\n",
              iq, length(unique(round(score_range, 8))), length(score_range)))
  cat(sprintf("[I2] Normal 组 score_range 中位数 = %.4f  阈值 ≤ %.2f  → %s\n",
              ref_med, REF_MEDIAN_MAX, if (i2) "过" else "不过"))
  cat(sprintf("[I3] %s (n=%d) vs %s (n=%d)  AUROC = %.4f  → %s\n",
              HOLDOUT_GROUP, length(neg), POS_GROUP, length(pos), a, v3))
}

# ——— §六 逐染色体均值差（无论 I3 过不过都报）———
ord_groups <- c(REF_GROUP, HOLDOUT_GROUP, "P4_AAH", "P4_AAH1", POS_GROUP)
ord_groups <- ord_groups[ord_groups %in% names(tb)]
chrmean <- sapply(ord_groups, function(gg) colMeans(M[g == gg, , drop = FALSE]))
diff_vec <- chrmean[, POS_GROUP] - chrmean[, HOLDOUT_GROUP]
cat("\n[§六] LUAD − holdout 的逐染色体均值差（log2）：\n")
print(round(diff_vec, 4))
n_big <- sum(abs(diff_vec) >= 0.1)
cat(sprintf("[§六] |差| ≥ 0.1 的染色体数 = %d / 22；最大值 = %.4f（%s）\n",
            n_big, diff_vec[which.max(abs(diff_vec))], names(which.max(abs(diff_vec)))))
cat(sprintf("[§六] 22 条染色体上 |差| 的平均 = %.4f，标准差 = %.4f\n",
            mean(abs(diff_vec)), sd(diff_vec)))
if (n_big <= 5 && max(abs(diff_vec)) >= 0.1) {
  cat("[§六] 读法（预注册）：差异集中在少数染色体且幅度 ≥0.1 ⇒ 与真实染色体级 CNV 一致\n")
} else if (n_big >= 15 && max(abs(diff_vec)) < 0.1) {
  cat("[§六] 读法（预注册）：差异均匀铺在 22 条染色体且幅度小 ⇒ 更像全局技术效应，不得说成 CNV\n")
} else {
  cat("[§六] 读法（预注册）：介于两者之间（集中度/幅度都不极端）⇒ 本预注册未定其判词，如实报原始曲线\n")
}

# ——— 产物 ———
per_cell <- data.frame(cell = colnames(e), group = g, score_range = score_range,
                       score_maxabs = score_maxabs, M,
                       stringsAsFactors = FALSE, check.names = FALSE)
fwrite(per_cell, file.path(D, "infercnv_smoke_v2_chrstat_per_cell.csv.gz"))

# 图 1：分组箱线图
png(file.path(FIGDIR, "v2_1_score_range_by_group.png"), width = 1400, height = 640, res = 130)
par(mar = c(6, 5, 4, 2))
dl <- lapply(ord_groups, function(gg) score_range[g == gg])
boxplot(dl, names = sprintf("%s\n(n=%d)", ord_groups, sapply(dl, length)),
        showfliers = FALSE, col = "grey85", las = 1,
        ylab = "per-cell score  =  max_k - min_k  of  mean log2 residual per chromosome",
        xlab = "group",
        main = "inferCNV smoke v2: chromosome-level imbalance score\n(higher = more deviation from the patient's own normal baseline)")
abline(h = ref_med, col = "red", lty = 2); abline(h = REF_MEDIAN_MAX, col = "blue", lty = 3)
legend("topright", c("Normal group median", "I2 threshold 0.05"),
       col = c("red", "blue"), lty = c(2, 3), bty = "n")
dev.off()

# 图 2：逐染色体均值曲线
png(file.path(FIGDIR, "v2_2_chr_profile_by_group.png"), width = 1500, height = 640, res = 130)
par(mar = c(6, 5, 4, 2))
cols <- c("grey40", "darkorange", "purple", "darkgreen", "red")
plot(NA, xlim = c(1, 22), ylim = range(chrmean), xaxt = "n", las = 1,
     xlab = "chromosome", ylab = "mean log2 residual (group mean over cells)",
     main = "inferCNV smoke v2: per-chromosome mean residual by group")
abline(h = 0, col = "grey60", lty = 3)
axis(1, at = 1:22, labels = CHRS)
for (i in seq_along(ord_groups)) {
  lines(1:22, chrmean[, i], col = cols[i], lwd = 2, type = "o", pch = 19, cex = .6)
}
legend("topright", ord_groups, col = cols[seq_along(ord_groups)], lwd = 2, bty = "n", ncol = 2)
dev.off()
cat(sprintf("\n[out] %s\n", file.path(FIGDIR, "v2_1_score_range_by_group.png")))
cat(sprintf("[out] %s\n", file.path(FIGDIR, "v2_2_chr_profile_by_group.png")))

summary_list <- list(
  script = "03_cnv/18_analyze_infercnv_chrstat.R",
  prereg = "results/03_cnv/infercnv_smoke/INFERCNV_SMOKE_PREREG_v2_CHRSTAT.md",
  inputs = list(
    expr_rds_sha256 = "2b0671e407e27f0fd2274e23a15d1b560ed0ef64d6481d80f408613b521229e7",
    gene_order_sha256 = "1f88af4ef06710a0e60632406c7a5db8f440b1764a55e191f31c3a5c8b097194",
    annotations_sha256 = "bff80e528f193b533424cff6d117b221eb49176f53d5ef4c7730df5e716a0485"),
  dim = c(n_genes = nrow(e), n_cells = ncol(e)),
  n_genes_in_stat = sum(keep),
  n_genes_unmatched = length(unmatched),
  cells_per_group = as.list(tb),
  score_definition = "score_range = max_k - min_k over chr1..chr22 of mean(log2 residual) for that chromosome",
  secondary_definition = "score_maxabs = max_k |mean(log2 residual) on chr k| (descriptive only)",
  I1_iqr = iq, I1_pass = i1,
  I2_ref_median = ref_med, I2_threshold = REF_MEDIAN_MAX, I2_pass = i2,
  I3_auroc = a, I3_strong = AUROC_STRONG, I3_weak = AUROC_WEAK, I3_verdict = v3,
  group_stats = st,
  chrmean_by_group = chrmean,
  chrom_diff_luad_minus_holdout = diff_vec,
  chrom_diff_n_abs_ge_0.1 = n_big,
  caveats = c(
    "同一个患者、同一次 inferCNV 跑的第二次分析，不是独立验证",
    "统计量换了 ⇒ v2 的 AUROC 与 v1 的 0.5000 不可比（v1 那个 0.5 是统计量退化造成的）",
    "snRNA 未在 inferCNV 上验证过；inferCNV 是可视化工具，2026 基准警告它分亚克隆而非逐细胞判决",
    "参考定义偏离已签口径（同患者 Normal 全部上皮，而非同亚型）—— 冒烟口径，不得引用",
    "单患者单次跑，不得推及全体",
    "环境 RNA / 内含子滞留会抬高残差；本数据无环境校正",
    "仍未回答：inferCNV 是否应被采纳为恶性口径"))
jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
if (jsonlite_ok) {
  jsonlite::write_json(summary_list, file.path(D, "infercnv_smoke_v2_chrstat_summary.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("[out] %s\n", file.path(D, "infercnv_smoke_v2_chrstat_summary.json")))
} else {
  saveRDS(summary_list, file.path(D, "infercnv_smoke_v2_chrstat_summary.rds"))
  cat("[warn] jsonlite 不可用，摘要存成 .rds\n")
}
cat("[done]\n")
