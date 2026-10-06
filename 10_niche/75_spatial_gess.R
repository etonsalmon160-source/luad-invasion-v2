#!/usr/bin/env Rscript
suppressMessages({library(signatureSearch); library(data.table); library(ExperimentHub)})
try(setExperimentHubOption("localHub", TRUE), silent=TRUE)
OUT<-"/home/eto/luad_v2/results/10_niche/sigsearch"; DB<-file.path(OUT,"lincs_db.h5")
gi<-fread(cmd="zcat /home/eto/lincs_data/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
sym2e<-setNames(as.character(gi$pr_gene_id), as.character(gi$pr_gene_symbol))
S<-fread(file.path(OUT,"spatial_global_depthmatched_clean_lfc.tsv"))
S<-S[lfc!=0]; S[,entrez:=sym2e[gene]]; S<-S[!is.na(entrez)]
Q<-matrix(S$lfc, ncol=1, dimnames=list(as.character(S$entrez),"spatial")) 
cat(sprintf("空转轴：%d 基因\n", nrow(Q)))
qs<-qSig(query=Q, gess_method="Cor", refdb=DB)
r<-gess_cor(qs, method="spearman", chunk_size=20000, workers=8)
fwrite(as.data.frame(r@result), file.path(OUT,"res_SPATIAL_cor.tsv"), sep="\t")
cat(sprintf("完成 %d 行\n", nrow(r@result)))
