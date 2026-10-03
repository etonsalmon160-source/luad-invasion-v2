#!/usr/bin/env Rscript
# 21_niche_prognosis.R —— M7a 试探：**生态位 marker 签名**在 TCGA-LUAD bulk 上的预后关联
#
# 用户口径（2026-10-03）：「用生态位标志物去做预后」——即
#   生态位 marker → 基因签名 → TCGA-LUAD bulk 打分 → Cox/KM。
#   （**不是**用患者级"域占比"——那个会被取材组成驱动，已否。）
#
# 🔴 本脚本是**探索性试探**，不是签字主结果。所有选择如下，逐条声明：
#   特征：每个域类型的 **top30 marker**（来自 d13_archetype_markers.tsv，域 vs 其余域）
#   打分：基因跨样本 z-score 后**取均值**（透明可审；非 ssGSEA）
#   终点：os_time_days / os_event
#   模型：单变量 Cox ＋ 多变量 Cox（校正 stage_num + age + sex）
#   分层：中位二分 + log-rank
#   **7 个域全部报出，不做挑选**（避免事后择优）
#
# 🔴 数据边界：只读 /home/eto/luad_invasion/data/TCGA_LUAD（只读源）；产物写 luad_v2。
#   样本限定在 m7a/tcga_clean_sample_list.tsv 的 483 例（177 死）。

suppressMessages({ library(data.table); library(survival) })

ROOT <- "/home/eto/luad_v2"
ARC  <- "/home/eto/luad_invasion/data/TCGA_LUAD"
M7A  <- file.path(ROOT, "results/10_niche/m7a")
OUT  <- file.path(ROOT, "results/10_niche/m7a")
NTOP <- 30L

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())

## ---- 临床 ----
SL <- fread(file.path(M7A, "tcga_clean_sample_list.tsv"))
step("清洗后样本 %d（死亡 %d）", nrow(SL), sum(SL$os_event == 1))

## ---- 表达（只读归档）----
E <- fread(file.path(ARC, "HiSeqV2.gz"), sep = "\t", header = TRUE)
step("表达矩阵 %d 基因 × %d 样本（Xena log2(RSEM+1)）", nrow(E) - 1L, ncol(E) - 1L)
sym <- E[[1]]; M <- as.matrix(E[, -1L, with = FALSE]); rownames(M) <- sym
keep <- intersect(SL$sampleID, colnames(M))
step("与临床交集 %d / %d", length(keep), nrow(SL))
M <- M[, keep, drop = FALSE]
CL <- SL[sampleID %in% keep][match(keep, sampleID)]

## ---- 签名：每个域类型的 top30 marker ----
MK <- fread(file.path(ROOT, "results/10_niche/kstar_diag/d13_archetype_markers.tsv"))
MK <- MK[rank <= NTOP]
SIG <- split(MK$gene, MK$archetype)
SC <- data.table(sampleID = keep)
for (a in names(SIG)) {
  g <- intersect(SIG[[a]], rownames(M))
  Z <- t(scale(t(M[g, , drop = FALSE])))          # 逐基因跨样本 z
  Z[!is.finite(Z)] <- 0
  SC[[paste0("D", a)]] <- colMeans(Z)
  step("  D%s：签名基因 %d/%d", a, length(g), length(SIG[[a]]))
}
fwrite(SC, file.path(OUT, "m7a_niche_signature_scores.tsv"), sep = "\t")

## ---- Cox + KM ----
D <- merge(CL, SC, by = "sampleID")
D[, `:=`(stage_num = as.numeric(stage_num), age = as.numeric(age))]
cols <- grep("^D[0-9]+$", names(D), value = TRUE)
res <- rbindlist(lapply(cols, function(v) {
  f1 <- as.formula(sprintf("Surv(os_time_days, os_event) ~ %s", v))
  f2 <- as.formula(sprintf("Surv(os_time_days, os_event) ~ %s + stage_num + age + sex", v))
  c1 <- coxph(f1, data = D); c2 <- coxph(f2, data = D)
  s1 <- summary(c1); s2 <- summary(c2)
  s <- survdiff(as.formula(sprintf("Surv(os_time_days, os_event) ~ I(%s > median(%s))", v, v)), data = D)
  p_lr <- 1 - pchisq(s$chisq, length(s$n) - 1)
  hi <- D[[v]] > median(D[[v]])
  data.table(sig = v, n_gene = nrow(MK[archetype == as.integer(sub("^D", "", v))]),
    HR_uni = s1$conf.int[1, 1], lo_uni = s1$conf.int[1, 3], hi_uni = s1$conf.int[1, 4],
    p_uni = s1$coefficients[1, 5],
    HR_adj = s2$conf.int[1, 1], lo_adj = s2$conf.int[1, 3], hi_adj = s2$conf.int[1, 4],
    p_adj = s2$coefficients[1, 5], p_logrank = p_lr,
    med_hi = median(D$os_time_days[hi]), med_lo = median(D$os_time_days[!hi]))
}))
setorder(res, p_uni)
fwrite(res, file.path(OUT, "m7a_niche_cox.tsv"), sep = "\t")

cat("\n=== M7a 生态位 marker 签名 × TCGA-LUAD 生存（483 例 / 177 死）===\n")
cat("🔴 探索性；打分 = top30 marker 的跨样本 z 均值；**7 个域全报，未挑选**\n\n")
print(as.data.frame(res[, .(sig, n_gene, HR_uni = round(HR_uni, 3), lo_uni = round(lo_uni, 3),
  hi_uni = round(hi_uni, 3), p_uni = signif(p_uni, 3), HR_adj = round(HR_adj, 3),
  lo_adj = round(lo_adj, 3), hi_adj = round(hi_adj, 3), p_adj = signif(p_adj, 3),
  p_logrank = signif(p_logrank, 3))]), row.names = FALSE)
cat("\n多重比较：7 个签名，Bonferroni 阈 = 0.05/7 =", round(0.05/7, 4), "\n")
cat("存活中位（签名高/低）：\n")
print(as.data.frame(res[, .(sig, med_hi_days = med_hi, med_lo_days = med_lo)]), row.names = FALSE)
cat("\n产物 → ", OUT, "\n", sep = "")
