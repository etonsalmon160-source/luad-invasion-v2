#!/usr/bin/env Rscript
# 62_tcga_pergene_hr.R —— TCGA-LUAD 逐基因调整 Cox，为「逆转签名质控」提供独立生存证据
#
# 目的：把「疾病里升高 → 该敲低」这一侧，按 **高表达是否真的预后差** 过滤。
#   HR_adj > 1（高表达⇒生存差）⇒ 敲低有道理 ⇒ 保留
#   HR_adj < 1（高表达⇒生存好）⇒ 敲低可能有害 ⇒ 剔除（保护性程序）
# 框架与 21_niche_prognosis.R 一致：逐基因跨样本 z → coxph(Surv ~ v + stage_num + age + sex)
suppressMessages({library(survival); library(data.table); library(parallel)})
setDTthreads(1)
ROOT <- "/home/eto/luad_v2"; ARC <- "/home/eto/luad_invasion/data/TCGA_LUAD"
OUT  <- file.path(ROOT, "results/10_niche/tr_singlecell")
N_CORE <- 16

E <- fread(cmd = sprintf("zcat %s/HiSeqV2.gz", ARC))
setnames(E, 1, "gene")
S <- fread(file.path(ROOT, "results/10_niche/m7a/tcga_clean_sample_list.tsv"))
common <- intersect(S$sampleID, names(E))
cat(sprintf("表达 %d 基因 × %d 样本；清洗后存活 %d 例\n", nrow(E), length(common), nrow(S)))
E <- E[, c("gene", common), with = FALSE]
X <- as.matrix(E[, -1]); rownames(X) <- E$gene
X <- t(scale(t(X)))                                     # 逐基因 z
D <- as.data.frame(S[match(common, sampleID)])
D$os_time_days <- as.numeric(D$os_time_days); D$os_event <- as.integer(D$os_event)
D$stage_num <- as.numeric(D$stage_num); D$age <- as.numeric(D$age); D$sex <- factor(D$sex)

one <- function(i) {
  g <- rownames(X)[i]
  d <- D; d$v <- X[i, ]
  s0 <- sd(d$v, na.rm = TRUE)
  if (is.na(s0) || s0 == 0) return(NULL)
  r <- tryCatch(summary(coxph(Surv(os_time_days, os_event) ~ v + stage_num + age + sex, data = d)),
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  data.table(gene = g, HR = r$conf.int[1, 1], lo = r$conf.int[1, 3], hi = r$conf.int[1, 4],
             p = r$coefficients[1, 5])
}
# 🔴 不用 mclapply：data.table 的 OpenMP 线程与 fork 冲突，子进程里全部报错（实测 one(1) 主进程正常、fork 后 0/20530）
raw <- vector("list", nrow(X))
for (i in seq_len(nrow(X))) {
  raw[[i]] <- one(i)
  if (i %% 4000 == 0) cat(sprintf("  %d / %d\n", i, nrow(X)))
}
raw <- Filter(function(z) !is.null(z) && is.data.frame(z), raw)
cat(sprintf("成功 %d / %d\n", length(raw), nrow(X)))
res <- rbindlist(raw)
setorder(res, p)
fwrite(res, file.path(OUT, "tcga_pergene_HR.tsv"), sep = "\t")
cat(sprintf("完成：%d 个基因有 HR；p<0.05 的 %d 个\n", nrow(res), sum(res$p < 0.05)))
cat("HR>1 且 p<0.05（高表达预后差）:", sum(res$HR > 1 & res$p < 0.05), "\n")
cat("HR<1 且 p<0.05（高表达预后好）:", sum(res$HR < 1 & res$p < 0.05), "\n")
