#!/usr/bin/env Rscript
# 🔴 已作废（superseded by 18_comp39_archetype.R）——**本脚本的数不可用，仅留档**
#    作废理由：k-means ＋ 把 56 张切片 pooled 在一起聚类 ＋ 抽样
#    ⇒ 聚出来主要是"切片"不是"生态位"，测的是跨切片混合聚类稳不稳，**不是**"组成空间能撑几个域"。
#    严谨版见 18_comp39_archetype.R（per-slide Seurat Leiden + 5 种子共识 + 与 §17 同款 K 扫描）。
#
# 17_rctd_finer_k.R —— M6 试探：39 维 RCTD 组成空间支持多少个**稳定**的 spot 簇？
# 目的：回答「RCTD 精细亚群那边能不能比 K*=7 更细」。
# 🔴 只读已落盘权重；产物写 kstar_diag/；日志**直接写文件**（不挂管道）
suppressMessages({library(data.table); library(mclust)})
ROOT <- "/home/eto/luad_v2"
DEC <- file.path(ROOT,"results/08_spatial_deconv/rctd_d/per_slide")
BN  <- file.path(ROOT,"results/10_niche/banksy")
OUT <- file.path(ROOT,"results/10_niche/kstar_diag")
PER <- 1200L; KS <- 4:16; NBOOT <- 25L
set.seed(20261003)
log <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...))); flush(stdout()) }
sl <- list.files(BN, pattern="^GSM")
X <- rbindlist(lapply(sl, function(s){
  w <- fread(file.path(DEC,sprintf("%s.weights.tsv.gz",s)),sep="\t",header=TRUE)
  setnames(w,1L,"barcode"); M <- as.matrix(w[,-1L,with=FALSE]); M[!is.finite(M)]<-0
  i <- sample(nrow(M), min(PER, nrow(M))); M <- M[i,,drop=FALSE]; M <- M/pmax(rowSums(M),1e-9)
  as.data.table(M)
}), fill=TRUE)
SUB <- names(X); Z <- scale(as.matrix(X)); Z[!is.finite(Z)] <- 0
log("抽样 spot %d，维 %d", nrow(Z), ncol(Z))
ST <- rbindlist(lapply(KS, function(k){
  t0 <- Sys.time(); full <- kmeans(Z, k, nstart=5, iter.max=100)$cluster
  s <- mean(replicate(NBOOT, { i <- sample(nrow(Z), floor(.8*nrow(Z)))
    a <- kmeans(Z[i,,drop=FALSE], k, nstart=5, iter.max=100)$cluster
    mclust::adjustedRandIndex(a, full[i]) }))
  log("  K=%2d  稳定性 %.3f  (%.0f s)", k, s, as.numeric(Sys.time()-t0, units="secs"))
  data.table(K=k, stab=s)
}))
fwrite(ST, file.path(OUT,"d10_rctd_spot_kstab.tsv"), sep="\t")
log("阈值 0.90 达标: %s", paste(ST[stab>=.90]$K, collapse=","))
log("阈值 0.85 达标: %s", paste(ST[stab>=.85]$K, collapse=","))
log("全部完成")
