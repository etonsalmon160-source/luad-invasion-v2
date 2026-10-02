#!/usr/bin/env Rscript
# 06_coverage_guard.R —— M6 §5.2 掩码覆盖守卫 + §5.3 患者级独立性登记（**只上报，不改任何口径**）
#
# 为什么必做（§5.2）：逐切片掩码剔除率 0.88%–29.6%；剔除率高的切片其**域覆盖面积本身就不完整**
#   ⇒ 报告必须**逐切片列出覆盖数**，**不得**把"域面积小"读成"这种生态位少"。
#
# §5.3：56 张来自 25 患者 ⇒ 任何跨切片汇总**必须按患者聚类**（同一患者的 2–3 张切片不是独立样本）。
#   本脚本只登记结构，不做统计（统计在各分析脚本里做）。
#
# 输入：results/08_spatial_deconv/spot_mask.tsv.gz（slide, barcode, in_tissue, total_umi, n_genes, pct_mt, pass）
#       results/10_niche/banksy/<slide>/domains_<agf>.tsv.gz（算进入分域的 spot 数）
# 输出：results/10_niche/guards/coverage_per_slide.tsv、coverage_summary.json
#
# 跑法（🔴 不能用 --vanilla：它会把 /usr/local/lib/R/site-library 剔出搜索路径，data.table 恰好在那）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 10_niche/06_coverage_guard.R

suppressMessages({ library(data.table) })

ROOT  <- "/home/eto/luad_v2"
NICHE <- file.path(ROOT, "results/10_niche")
OUT   <- file.path(NICHE, "guards"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())

sm <- fread(file.path(ROOT, "results/08_spatial_deconv/spot_mask.tsv.gz"))
sm[, pass := as.logical(pass)]
step("掩码表 %d spot / %d 张", nrow(sm), uniqueN(sm$slide))

## 逐切片覆盖
## 🔴 2026-10-02 修：median() 对整数向量**按元素个数奇偶返回不同类型**
##    （奇数个 ⇒ integer，偶数个 ⇒ double；`median(c(1L,2L))` = 1.5 是 double，
##      `median(c(1L,2L,3L))` = 2L 是 integer）⇒ data.table 分组时各切片类型不一，
##    直接报 "Column 3 of result for group 2 is type 'integer' but expecting type 'double'"
##    （09:32:50 / 09:41:29 两次 rc=1 的根因）。**只统一存储类型，数值不变。**
cov <- sm[, .(n_in_tissue = sum(in_tissue == 1, na.rm = TRUE),
              n_mask_pass = sum(pass, na.rm = TRUE),
              umi_med_in_tissue = as.numeric(median(total_umi[in_tissue == 1], na.rm = TRUE)),
              umi_med_pass     = as.numeric(median(total_umi[pass], na.rm = TRUE))),
          by = slide]
cov[, n_drop := n_in_tissue - n_mask_pass]
cov[, drop_rate := ifelse(n_in_tissue > 0, n_drop / n_in_tissue, NA_real_)]
cov[, patient := sub("^GSM[0-9]+_(P[0-9]+)_.*$", "\\1", slide)]
cov[, stage   := sub("^GSM[0-9]+_P[0-9]+_", "", slide)]

## 实际进入分域的 spot 数（取任一 AGF 版；缺就当 NA，不硬停——分域可能还没跑完）
for (agf in c("agfT", "agfF")) {
  n <- vapply(cov$slide, function(s) {
    f <- file.path(NICHE, "banksy", s, sprintf("domains_%s.tsv.gz", agf))
    if (!file.exists(f)) return(NA_integer_)
    d <- fread(f, select = "barcode"); uniqueN(d$barcode) }, integer(1))
  cov[[sprintf("n_spot_analyzed_%s", agf)]] <- n
}

setorder(cov, -drop_rate)
fwrite(cov, file.path(OUT, "coverage_per_slide.tsv"), sep = "\t")

step("==== §5.2 逐切片覆盖（按剔除率降序，全表见 coverage_per_slide.tsv）====")
for (i in seq_len(min(8L, nrow(cov))))
  step("  %-24s 组织 %5d  过掩码 %5d  剔除 %5.1f%%  深度中位 %s",
       cov$slide[i], cov$n_in_tissue[i], cov$n_mask_pass[i], 100 * cov$drop_rate[i],
       format(cov$umi_med_pass[i], big.mark = ","))
step("  剔除率范围 %.2f%% – %.2f%%（中位 %.2f%%）",
     100 * min(cov$drop_rate, na.rm = TRUE), 100 * max(cov$drop_rate, na.rm = TRUE),
     100 * median(cov$drop_rate, na.rm = TRUE))

## §5.3 结构登记：每患者几张切片、哪些期别
pat <- cov[, .(n_slide = .N, stages = paste(sort(unique(stage)), collapse = "+"),
               n_spot_pass = sum(n_mask_pass)), by = patient]
setorder(pat, patient)
fwrite(pat, file.path(OUT, "patient_slide_structure.tsv"), sep = "\t")
step("==== §5.3 患者×切片结构：%d 患者 / %d 张 ====", nrow(pat), nrow(cov))
step("  每患者切片数分布：%s", paste(sprintf("%d张=%d人", as.integer(names(table(pat$n_slide))),
                                          as.integer(table(pat$n_slide))), collapse = "  "))
step("  🔴 任何跨切片汇总必须按 patient 聚类（同一患者多张切片不是独立样本）")

## 高剔除切片登记（>15%）
hi <- cov[drop_rate > 0.15]
if (nrow(hi)) {
  step("  ⚠️ 剔除率 >15%% 的切片 %d 张（域覆盖本就不完整，不得读成「生态位少」）：", nrow(hi))
  for (i in seq_len(nrow(hi)))
    step("     %-24s %.1f%%", hi$slide[i], 100 * hi$drop_rate[i])
}

man <- list(step = "§5.2 掩码覆盖守卫 + §5.3 患者级独立性登记", n_slide = nrow(cov),
            n_patient = nrow(pat),
            drop_rate_range = c(min = min(cov$drop_rate, na.rm = TRUE),
                                max = max(cov$drop_rate, na.rm = TRUE),
                                median = median(cov$drop_rate, na.rm = TRUE)),
            n_slide_drop_gt15pct = nrow(hi),
            caveat = "域面积小 ≠ 该生态位少（§5.2）；跨切片汇总必须按患者（§5.3）",
            finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE),
           file.path(OUT, "coverage_summary.json"))
step("==== 06 结束 ====")
