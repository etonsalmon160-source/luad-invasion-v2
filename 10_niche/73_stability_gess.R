#!/usr/bin/env Rscript
# 73_stability_gess.R —— 对 16 条劈半轴各跑一次 gess_cor（Spearman），用于稳定性检验
suppressMessages({library(signatureSearch); library(data.table); library(ExperimentHub)})
try(setExperimentHubOption("localHub", TRUE), silent=TRUE)   # 🔴 不让它去连 bioconductor.org
OUT <- "/home/eto/luad_v2/results/10_niche/sigsearch"
SP  <- file.path(OUT, "splits")
DB  <- file.path(OUT, "lincs_db.h5")
gi <- fread(cmd = "zcat /home/eto/lincs_data/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
sym2e <- setNames(as.character(gi$pr_gene_id), as.character(gi$pr_gene_symbol))
E <- fread(file.path(SP, "split_axes.tsv"))
IDX <- fread(file.path(SP, "split_index.tsv"))
cols <- intersect(IDX$name, names(E))
cat(sprintf("待跑 %d 条轴\n", length(cols)))
for (cn in cols) {
  f_out <- file.path(SP, sprintf("cor_%s.tsv", cn))
  if (file.exists(f_out)) { cat(sprintf("  [%s] 已有，跳过\n", cn)); next }
  v <- E[[cn]]; names(v) <- E$gene
  v <- v[!is.na(v) & v != 0]
  ez <- sym2e[names(v)]; ok <- !is.na(ez)
  Q <- matrix(v[ok], ncol = 1, dimnames = list(as.character(ez[ok]), "q"))
  qs <- qSig(query = Q, gess_method = "Cor", refdb = DB)
  r  <- gess_cor(qs, method = "spearman", chunk_size = 20000, workers = 8)
  fwrite(as.data.frame(r@result), file.path(SP, sprintf("cor_%s.tsv", cn)), sep = "\t")
  cat(sprintf("  [%s] %d 基因 → %d 行\n", cn, nrow(Q), nrow(r@result)))
}
cat("全部完成\n")
