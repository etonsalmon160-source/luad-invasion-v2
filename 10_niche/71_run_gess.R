#!/usr/bin/env Rscript
# 71_run_gess.R —— 用 signatureSearch 官方包跑三种方法（CMap / Cor / gCMap）
suppressMessages({library(signatureSearch); library(data.table)})
OUT <- "/home/eto/luad_v2/results/10_niche/sigsearch"
DB  <- file.path(OUT, "lincs_db.h5")
SC  <- "/home/eto/luad_v2/results/10_niche/tr_singlecell"

# ① 查询：质控后的全局签名（符号 → Entrez）
gi <- fread(cmd = "zcat /home/eto/lincs_data/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
sym2e <- setNames(as.character(gi$pr_gene_id), as.character(gi$pr_gene_symbol))
S <- fread(file.path(SC, "global_paired_qc.tsv"))
S <- S[lfc_qc != 0]
S[, entrez := sym2e[gene]]
S <- S[!is.na(entrez)]
up <- S[order(-lfc_qc)][1:150]; dn <- S[order(lfc_qc)][1:150]
up <- unique(up[!is.na(entrez), entrez]); dn <- unique(dn[!is.na(entrez), entrez])
cat(sprintf("查询：UP %d / DOWN %d 个 Entrez\n", length(up), length(dn)))

# ② CMap 口径（上下调基因表）
qs1 <- qSig(query = list(upset = up, downset = dn), gess_method = "CMAP", refdb = DB)
r1  <- gess_cmap(qs1, chunk_size = 20000, workers = 8)
fwrite(as.data.frame(r1@result), file.path(OUT, "res_CMAP.tsv"), sep = "\t")
cat(sprintf("[CMAP] %d 行；最负 5 个: %s\n", nrow(r1@result),
            paste(head(r1@result$pert, 5), collapse = ", ")))

# ③ 相关口径（数值矩阵）
qm <- S[, .(val = median(lfc_qc)), by = entrez]
Q  <- matrix(qm$val, ncol = 1, dimnames = list(as.character(qm$entrez), "LUAD_global"))
for (mth in c("spearman", "pearson")) {
  qs2 <- qSig(query = Q, gess_method = "Cor", refdb = DB)
  r2  <- gess_cor(qs2, method = mth, chunk_size = 20000, workers = 8)
  fwrite(as.data.frame(r2@result), file.path(OUT, sprintf("res_COR_%s.tsv", mth)), sep = "\t")
  cat(sprintf("[Cor-%s] %d 行；最负 5: %s\n", mth, nrow(r2@result),
              paste(head(r2@result$pert, 5), collapse = ", ")))
}
cat("完成\n")
