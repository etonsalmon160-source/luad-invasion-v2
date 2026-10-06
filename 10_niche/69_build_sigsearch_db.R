#!/usr/bin/env Rscript
# 69_build_sigsearch_db.R —— 把 LINCS 建成 signatureSearch 的 HDF5 库
# 列名格式 (drug)__(cell)__(factor)，行名 = Entrez（包内以 character 比对）
suppressMessages({library(arrow); library(data.table); library(signatureSearch); library(rhdf5)})
setDTthreads(4)
OUT <- "/home/eto/luad_v2/results/10_niche/sigsearch"
t0 <- Sys.time()
cat(sprintf("[%s] 读 feather …\n", format(t0, "%H:%M:%S")))
T <- as.data.table(arrow::read_feather(file.path(OUT, "lincs_ref.feather")))
sig <- T$signature; T$signature <- NULL
cat(sprintf("  %d 签名 × %d 基因；列名唯一? %s\n", nrow(T), ncol(T),
            length(unique(sig)) == length(sig)))
genes <- colnames(T)
M <- t(as.matrix(T))                                    # 基因 × 签名
rm(T); gc()
rownames(M) <- genes; colnames(M) <- sig
cat(sprintf("  转置完成 %d × %d（%.1f GB）\n", nrow(M), ncol(M), as.numeric(object.size(M)) / 1e9))
cat(sprintf("[%s] build_custom_db …\n", format(Sys.time(), "%H:%M:%S")))
build_custom_db(df = M, h5file = file.path(OUT, "lincs_db.h5"))
cat(sprintf("[%s] 完成；文件 %.1f GB\n", format(Sys.time(), "%H:%M:%S"),
            file.size(file.path(OUT, "lincs_db.h5")) / 1e9))
