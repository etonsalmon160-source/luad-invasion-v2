#!/usr/bin/env Rscript
# 52_nichenet.R —— 路 2：NicheNet「哪个配体在维持这个状态」
#
# 设定：
#   receiver = Fibroblast（成纤维）
#   senders  = D3 域里丰度最高的非成纤维细胞型（d3_senders.txt，来自 d7_archetypes.tsv 的组成）
#   target   = **数据导出**的疾病轴 up 基因 top200（患者内配对 IAC−前驱），**不手写**
# 口径（本项目选择，登记为 S 档）：
#   「表达」判据 = 该细胞型的平均 log1p(CP10K) > 0.1（约相当于 ~10% 细胞检出）
#   potential_ligands = 配体-受体网络 ∩ sender 表达 ∩ 配体-靶矩阵行
# 参考：NicheNet 的 aupr_corrected 已内含"随机基因集"的基线校正

suppressMessages({library(nichenetr)})
OUT <- "/home/eto/luad_v2/results/10_niche/tr_singlecell"
DAT <- "/home/eto/nichenet_data"
set.seed(20261005)

ltm <- readRDS(file.path(DAT, "ligand_target_matrix_nsga2r_final.rds"))
lr  <- readRDS(file.path(DAT, "lr_network_human_21122021.rds"))
cat("配体-靶矩阵:", nrow(ltm), "×", ncol(ltm), " | 配体-受体:", nrow(lr), "\n")

E <- read.delim(file.path(OUT, "L2_mean_expr.tsv"), row.names = 1, check.names = FALSE)
cat("表达表:", nrow(E), "型 ×", ncol(E), "基因\n")
tgt <- readLines(file.path(OUT, "target_up200.txt"))
senders <- readLines(file.path(OUT, "d3_senders.txt"))
senders <- senders[senders %in% rownames(E)]
receiver <- "Fibroblast"
cat("receiver:", receiver, " | senders:", paste(senders, collapse = ", "), "\n")

THR <- 0.1
expressed_rec <- colnames(E)[E[receiver, ] > THR]
expressed_sen <- colnames(E)[apply(E[senders, , drop = FALSE], 2, max) > THR]
cat("表达基因：receiver", length(expressed_rec), " | senders", length(expressed_sen), "\n")

potential_ligands <- lr$from[lr$from %in% expressed_sen]
potential_ligands <- intersect(potential_ligands, rownames(ltm))
cat("候选配体（sender 表达 ∩ 网络 ∩ 矩阵）:", length(potential_ligands), "\n")

tgt_lm <- intersect(tgt, colnames(ltm))
cat("靶基因：给 200 个，落在配体-靶矩阵里的", length(tgt_lm), "个\n")
bg <- intersect(expressed_rec, colnames(ltm))
LTM <- t(ltm[potential_ligands, bg, drop = FALSE])  # 🔴 必须转置：nichenetr 2.2.1.1 默认 ligands_position="cols"
cat("用于打分的矩阵(转置后):", nrow(LTM), "背景 ×", ncol(LTM), "配体\n")

res <- predict_ligand_activities(
  geneset = tgt_lm,
  background_expressed_genes = bg,
  ligand_target_matrix = LTM,
  potential_ligands = potential_ligands)

res <- res[order(-res$aupr_corrected), ]
write.table(res, file.path(OUT, "nichenet_ligand_activities.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== 配体活性 top 20（按 aupr_corrected）===\n")
print(head(res[, c("test_ligand", "aupr_corrected", "aupr", "auroc", "pearson")], 20))

# 最强配体的下游靶链接
best <- head(res$test_ligand, 5)
active_ligand_target_links <- best_ligands_adapted <- NULL
links <- get_weighted_ligand_target_links(
  ligand = best, geneset = tgt_lm, ligand_target_matrix = LTM, n = 250)
if (!is.null(links) && nrow(links) > 0) {
  write.table(links, file.path(OUT, "nichenet_ligand_target_links.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n写出配体-靶链接:", nrow(links), "行（配体:", paste(unique(links$ligand), collapse = ", "), "）\n")
}
cat("\n=== 已落盘 ===\n")
cat("  nichenet_ligand_activities.tsv\n")
