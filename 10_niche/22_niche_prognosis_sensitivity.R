#!/usr/bin/env Rscript
# 22_niche_prognosis_sensitivity.R —— M7a §22 的两条待办（2026-10-03）
#
#   ① 打分法敏感性：mean-z → **ssGSEA**（GSVA::gsva method="ssgsea"），看结论是否翻转
#   ② D3 归因：把 D3 的预后信号**拆成"生态位特有" vs "肿瘤含量/增殖"**
#      —— 用**正常肺结构签名**（D1 气道＋D4 肺泡壁）当**纯度反指标**，
#         再用**经典细胞周期基因**当**增殖**协变量，看 D3 的 HR 是否还站得住
#
# 🔴 探索性延续，仍不是签字主结果。所有 7 个签名一律并报，不挑选。
# 依赖：GSVA 1.46.0（已装）、survival 3.4.0

suppressMessages({ library(data.table); library(survival); library(GSVA) })

ROOT <- "/home/eto/luad_v2"; ARC <- "/home/eto/luad_invasion/data/TCGA_LUAD"
M7A  <- file.path(ROOT, "results/10_niche/m7a")
KD   <- file.path(ROOT, "results/10_niche/kstar_diag")
NTOP <- 30L

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())

## ---- 表达 + 临床 ----
SL <- fread(file.path(M7A, "tcga_clean_sample_list.tsv"))
E  <- fread(file.path(ARC, "HiSeqV2.gz"), sep = "\t", header = TRUE)
sym <- E[[1]]; M <- as.matrix(E[, -1L, with = FALSE]); rownames(M) <- sym
keep <- intersect(SL$sampleID, colnames(M)); M <- M[, keep, drop = FALSE]
CL <- SL[sampleID %in% keep][match(keep, sampleID)]
step("表达 %d 基因 × %d 样本；死亡 %d", nrow(M), ncol(M), sum(CL$os_event == 1))

## ---- 签名集 ----
MK  <- fread(file.path(KD, "d13_archetype_markers.tsv"))[rank <= NTOP]
SIG <- lapply(split(MK$gene, MK$archetype), function(g) intersect(g, rownames(M)))
names(SIG) <- paste0("D", names(SIG))
## 纯度反指标 = 正常肺结构（D1 气道 + D4 肺泡壁）
SIG[["NORMAL_LUNG"]] <- intersect(unique(c(SIG[["D1"]], SIG[["D4"]])), rownames(M))
## 增殖 = 经典细胞周期基因（通用，非本数据集拟合）
SIG[["PROLIF"]] <- intersect(c("MKI67","TOP2A","CCNB1","CCNB2","BUB1","BUB1B","AURKA","AURKB",
  "PLK1","CDC20","CDK1","CCNE1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","RRM2","TYMS",
  "CDKN3","UBE2C","TPX2","KIF11","ASPM","CENPF","ANLN"), rownames(M))
step("签名基因命中：%s", paste(sprintf("%s=%d", names(SIG), sapply(SIG, length)), collapse=" "))

## ---- ① ssGSEA ----
set.seed(1)
## 🔴 GSVA 1.46 用旧 API（`gsvaParam()` 是 1.47+ 才有的）
SS <- gsva(M, SIG, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
SS <- t(SS)
SC <- data.table(sampleID = keep, as.data.table(SS))
fwrite(SC, file.path(M7A, "m7a_niche_signature_scores_ssgsea.tsv"), sep = "\t")

## ---- ② Cox：三种模型 ----
D <- merge(CL, SC, by = "sampleID")
D[, `:=`(stage_num = as.numeric(stage_num), age = as.numeric(age))]
sigs <- grep("^D[0-9]+$", names(D), value = TRUE)
cox1 <- function(v, extra = "") {
  f <- as.formula(sprintf("Surv(os_time_days, os_event) ~ %s%s", v,
                          if (nzchar(extra)) paste0(" + ", extra) else ""))
  s <- summary(coxph(f, data = D))
  c(HR = s$conf.int[1,1], lo = s$conf.int[1,3], hi = s$conf.int[1,4], p = s$coefficients[1,5])
}
R <- rbindlist(lapply(sigs, function(v) {
  m1 <- cox1(v); m2 <- cox1(v, "stage_num + age + sex")
  m3 <- cox1(v, "stage_num + age + sex + NORMAL_LUNG + PROLIF")
  data.table(sig = v,
    HR_ssgsea = m1[1], p_ssgsea = m1[4],
    HR_adj = m2[1], p_adj = m2[4],
    HR_adj_purity_prolif = m3[1], p_adj_purity_prolif = m3[4])
}))
setorder(R, p_adj)
fwrite(R, file.path(M7A, "m7a_niche_cox_ssgsea.tsv"), sep = "\t")

cat("\n=== ① ssGSEA 打分（对 mean-z 的敏感性）===\n")
cat(sprintf("%-5s %10s %9s %11s %9s %18s %9s\n", "签名","HR_ssGSEA","p","HR_adj","p_adj","HR_adj+纯度+增殖","p"))
for (i in seq_len(nrow(R))) with(R[i,], cat(sprintf("%-5s %10.3f %9.2e %11.3f %9.2e %18.3f %9.2e\n",
  sig, HR_ssgsea, p_ssgsea, HR_adj, p_adj, HR_adj_purity_prolif, p_adj_purity_prolif)))

## ---- 对照：mean-z 的同一张表 ----
OLD <- fread(file.path(M7A, "m7a_niche_cox.tsv"))[, .(sig, HR_mz = HR_adj, p_mz = p_adj)]
CMP <- merge(R[, .(sig, HR_ssgsea, p_ssgsea, HR_adj, p_adj)], OLD, by = "sig", all.x = TRUE)
cat("\n=== 打分法一致性（mean-z vs ssGSEA，校正后 HR）===\n")
for (i in seq_len(nrow(CMP))) with(CMP[i,], cat(sprintf("  %-4s  mean-z %6.3f (p=%.1e)   ssGSEA %6.3f (p=%.1e)\n",
  sig, HR_mz, p_mz, HR_adj, p_adj)))

cat("\n=== ② D3 归因：加入「正常肺含量（纯度反指标）」＋「增殖」后 ===\n")
d3 <- R[sig == "D3"]
cat(sprintf("  D3  校正(分+龄+性)            HR %.3f  p=%.2e\n", d3$HR_adj, d3$p_adj))
cat(sprintf("  D3  + 正常肺含量 + 增殖        HR %.3f  p=%.2e\n", d3$HR_adj_purity_prolif, d3$p_adj_purity_prolif))
cat("\n  各签名之间的相关（Spearman，含协变量）：\n")
CC <- cor(D[, c(sigs, "NORMAL_LUNG", "PROLIF"), with = FALSE], method = "spearman")
print(round(CC["D3", , drop = FALSE], 3))
step("产物 → %s", M7A)
