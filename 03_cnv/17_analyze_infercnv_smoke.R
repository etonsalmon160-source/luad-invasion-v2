#!/usr/bin/env Rscript
# 分析 inferCNV 冒烟结果，判读**写死在 INFERCNV_SMOKE_PREREG.md §七**，本脚本只执行。
#
# 跑法：
#   LD_LIBRARY_PATH=/home/eto/local/jags/lib:$LD_LIBRARY_PATH \
#   R_LIBS=/home/eto/Rlibs/infercnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 03_cnv/17_analyze_infercnv_smoke.R

suppressMessages({
  library(data.table)
})

ROOT <- "/home/eto/luad_v2"
D <- file.path(ROOT, "results/03_cnv/infercnv_smoke")
FIGDIR <- file.path(D, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ——— 预注册常量（INFERCNV_SMOKE_PREREG.md §五 / §七）———
REF_GROUP <- "Normal"
HOLDOUT_GROUP <- "Normal_holdout"
POS_GROUP <- "P4_LUAD"
GENE_RETENTION_MIN <- 0.90
REF_MEDIAN_LO <- 0.95
REF_MEDIAN_HI <- 1.05
AUROC_STRONG <- 0.80
AUROC_WEAK <- 0.60
N_GENES_IN <- 18069

# 秩和法 AUROC（避免额外依赖）
auroc <- function(pos, neg) {
  y <- c(rep(1, length(pos)), rep(0, length(neg)))
  s <- c(pos, neg)
  r <- rank(s, ties.method = "average")
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cat("[in] 读残差矩阵 ...\n")
expr <- readRDS(file.path(D, "p4_final_expr_data.rds"))
cat(sprintf("[in] %d 基因 x %d 细胞\n", nrow(expr), ncol(expr)))

# ⚠️ 必须显式 sep="\t"：barcode 里带 `|`，fread 自动猜分隔符会猜成 `|`，
# 把 `AAACAA...1|P4_AAH<TAB>P4_AAH` 错切成 cell=`AAACAA...1`、group=`P4_AAH<TAB>P4_AAH`。
ann <- fread(file.path(D, "p4_annotations.tsv"), header = FALSE, sep = "\t",
             col.names = c("cell", "group"), data.table = FALSE)
stopifnot(!anyDuplicated(ann$cell), setequal(ann$cell, colnames(expr)))
# 按名字对齐（不依赖 inferCNV 是否重排过细胞顺序）
grp <- setNames(ann$group, ann$cell)
print(table(grp[colnames(expr)]))

cat(sprintf("[chk] 残差取值范围：min=%.4g  中位=%.4g  max=%.4g\n",
            min(expr), median(expr), max(expr)))
if (min(expr) <= 0) warning("残差里出现非正值，log2 会产生 -Inf/NaN；下面用 NA 屏蔽")

# ——— 逐细胞打分：常染色体基因上 |log2(残差)| 的中位数（预注册 §七 写死）———
# inferCNV 已按默认 chr_exclude 去掉了 chrX/chrY/chrM，故矩阵内基因即常染色体基因
lg <- log2(expr)
lg[!is.finite(lg)] <- NA
score <- apply(lg, 2, median, na.rm = TRUE)
n_na <- sum(!is.finite(score))
if (n_na) cat(sprintf("[warn] %d 个细胞打分失败（全 NA）\n", n_na))

per_cell <- data.frame(cell = colnames(expr), group = grp[colnames(expr)],
                       score = as.numeric(score), stringsAsFactors = FALSE)
fwrite(per_cell, file.path(D, "infercnv_smoke_per_cell.csv.gz"))

cat("\n[stat] 各组逐细胞打分：\n")
st <- do.call(rbind, lapply(split(per_cell$score, per_cell$group), function(v) {
  data.frame(n = length(v), median = median(v, na.rm = TRUE),
             mean = mean(v, na.rm = TRUE),
             p90 = quantile(v, 0.9, na.rm = TRUE, names = FALSE))
}))
st <- cbind(group = rownames(st), st)
print(st, row.names = FALSE)

# ——— 门槛 1：能跑通 + 基因保留率 ———
gene_ret <- nrow(expr) / N_GENES_IN
th1 <- gene_ret >= GENE_RETENTION_MIN

# ——— 门槛 2：参考组残差中位数落在 [0.95, 1.05] ———
ref_med <- median(expr[, grp[colnames(expr)] == REF_GROUP])
th2 <- ref_med >= REF_MEDIAN_LO && ref_med <= REF_MEDIAN_HI

# ——— 门槛 3：Normal_holdout vs P4_LUAD 逐细胞打分 AUROC ———
pos <- per_cell$score[per_cell$group == POS_GROUP]
neg <- per_cell$score[per_cell$group == HOLDOUT_GROUP]
a <- auroc(pos, neg)
if (a >= AUROC_STRONG) {
  v3 <- "有信号，值得进下一轮（仍须另行签字）"
} else if (a >= AUROC_WEAK) {
  v3 <- "弱信号"
} else {
  v3 <- "分不开"
}

cat(sprintf("\n[门槛1] 基因保留 %.2f%%（%d/%d） 阈值 %.0f%%  → %s\n",
            100 * gene_ret, nrow(expr), N_GENES_IN, 100 * GENE_RETENTION_MIN,
            if (th1) "过" else "不过"))
cat(sprintf("[门槛2] 参考组残差中位数 %.4f  区间 [%.2f, %.2f]  → %s\n",
            ref_med, REF_MEDIAN_LO, REF_MEDIAN_HI, if (th2) "过" else "不过"))
cat(sprintf("[门槛3] %s (n=%d) vs %s (n=%d)  AUROC = %.4f  → %s\n",
            HOLDOUT_GROUP, length(neg), POS_GROUP, length(pos), a, v3))

# ——— 图（图内标签全英文）———
ORD <- c(REF_GROUP, HOLDOUT_GROUP, "P4_AAH", "P4_AAH1", POS_GROUP)
ORD <- ORD[ORD %in% unique(per_cell$group)]
png(file.path(FIGDIR, "infercnv_smoke_score_by_group.png"),
    width = 1400, height = 620, res = 130)
par(mar = c(6, 5, 4, 2))
data <- lapply(ORD, function(g) per_cell$score[per_cell$group == g])
bp <- boxplot(data, names = sprintf("%s\n(n=%d)", ORD, sapply(data, length)),
              showfliers = FALSE, col = "grey85",
              ylab = "per-cell score  =  median |log2 residual|",
              xlab = "group", las = 1,
              main = "inferCNV smoke test on P4 epithelium\n(higher = more deviation from the patient's own normal baseline)")
abline(h = 0, col = "grey50", lty = 3)
dev.off()
cat(sprintf("\n[out] %s\n", file.path(FIGDIR, "infercnv_smoke_score_by_group.png")))

summary_list <- list(
  script = "03_cnv/17_analyze_infercnv_smoke.R",
  expr_dim = c(n_genes = nrow(expr), n_cells = ncol(expr)),
  gene_retention = gene_ret,
  threshold1_pass = th1,
  ref_median_residual = ref_med,
  threshold2_pass = th2,
  auroc_holdout_vs_luad = a,
  auroc_strong = AUROC_STRONG,
  auroc_weak = AUROC_WEAK,
  threshold3_verdict = v3,
  group_stats = st,
  score_definition = "median over autosomal genes of |log2(residual)|",
  caveats = c(
    "snRNA 从未在 inferCNV 上验证过",
    "inferCNV 是可视化工具；2026 基准警告它分亚克隆而非逐细胞判决",
    "参考定义偏离已签口径（同患者 Normal 全部上皮，而非同亚型）—— 冒烟口径，不得引用",
    "单患者单次跑，不得推及全体",
    "Normal_holdout 与 LUAD 仍共享患者级样本批次，本设计无法排除",
    "环境 RNA / 内含子滞留会抬高残差"
  )
)
jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
if (jsonlite_ok) {
  jsonlite::write_json(summary_list, file.path(D, "infercnv_smoke_summary.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("[out] %s\n", file.path(D, "infercnv_smoke_summary.json")))
} else {
  saveRDS(summary_list, file.path(D, "infercnv_smoke_summary.rds"))
  cat("[warn] jsonlite 不可用，摘要存成 .rds\n")
}
cat("[done]\n")
