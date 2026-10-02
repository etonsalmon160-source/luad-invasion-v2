#!/usr/bin/env Rscript
# 07_m7a_tcga_clean.R —— M7a：TCGA-LUAD 探针样本清洗（S11 六条规则）
#
# 🔴 **本脚本只做清洗，不做打分、不做生存**（S11 规则 ⑥ 的硬次序要求）：
#    先只跑清洗、把「保留患者数 + 各期计数 + 各步丢弃数」打出来核对，与规则一起落盘冻结，
#    **之后**才允许跑 ssGSEA 打分与 Cox。顺序反了即构成事后调参（法则 3.2）。
#
# 数据来源（只读归档，S10 已签）：/home/eto/luad_invasion/data/TCGA_LUAD/
#   HiSeqV2.gz（表达，log2(RSEM+1)，576 样本）、LUAD_clinicalMatrix（临床，706 样本）
#   🔴 读取前先落 sha256；**产物一律写 luad_v2 侧**，绝不在 luad_invasion 下产出。
#
# S11 规则（照抄，不得改）：
#   ① 只留原发瘤（sample_type = "Primary Tumor"）
#   ② 每个患者只留一个样本（重复 aliquot 按 barcode **字典序**取第一个；
#      ——— 不得按表达量或生存期挑，那是选择偏倚 ———）
#   ③ 必须有 pathologic_stage ∈ {I, II, III, IV}；缺失/[Discrepancy] 的丢
#   ④ 必须有生存时间（days_to_death 或 days_to_last_followup）且 > 0
#   ⑤ 年龄 / 性别缺失的丢
#   ⑥ 打印保留患者数 + 各期计数 + 各步丢弃数，落盘冻结；**然后停**
#
# 跑法（🔴 不能用 --vanilla：它会把 /usr/local/lib/R/site-library 剔出搜索路径，data.table 恰好在那）：
#   R_LIBS=/home/eto/Rlibs/fastcnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library \
#   Rscript 10_niche/07_m7a_tcga_clean.R

suppressMessages({ library(data.table) })

ARCH <- "/home/eto/luad_invasion/data/TCGA_LUAD"     # 只读
ROOT <- "/home/eto/luad_v2"
OUT  <- file.path(ROOT, "results/10_niche/m7a"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

F_EXPR <- file.path(ARCH, "HiSeqV2.gz")
F_CLIN <- file.path(ARCH, "LUAD_clinicalMatrix")
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())

## ——— 0. sha256（S10 要求，读之前先落）———
sha <- function(f) system2("sha256sum", shQuote(f), stdout = TRUE)
sha_expr <- sha(F_EXPR); sha_clin <- sha(F_CLIN)
step("sha256(HiSeqV2.gz)        = %s", sub(" .*$", "", sha_expr))
step("sha256(LUAD_clinicalMatrix)= %s", sub(" .*$", "", sha_clin))

## ——— 1. 表达矩阵的样本列（决定"有没有表达数据"这个上限）———
read_gz_line <- function(p, n = 1L) {
  con <- gzfile(p, "rt"); on.exit(close(con)); readLines(con, n = n)
}
hdr <- strsplit(read_gz_line(F_EXPR, 1L), "\t")[[1]]
expr_samples <- hdr[-1]
step("表达矩阵：%d 样本列（第 1 列是 sample）", length(expr_samples))

## ——— 2. 临床表 ———
cl <- fread(F_CLIN, sep = "\t", check.names = FALSE, na.strings = c("", "NA", "[Not Available]"))
step("临床表：%d 行 x %d 列", nrow(cl), ncol(cl))
need <- c("sampleID", "bcr_patient_barcode", "sample_type", "pathologic_stage",
          "days_to_death", "days_to_last_followup", "vital_status",
          "age_at_initial_pathologic_diagnosis", "gender")
miss <- setdiff(need, names(cl))
if (length(miss)) stop("临床表缺列：", paste(miss, collapse = ", "), call. = FALSE)

n0 <- nrow(cl)
step("==== S11 逐步清洗（起点 %d 行）====", n0)

## ① 只留原发瘤
cl <- cl[sample_type == "Primary Tumor"]
n1 <- nrow(cl); step("① 原发瘤         保留 %5d（丢弃 %5d）", n1, n0 - n1)

## 与表达矩阵取交集（没有表达数据的样本无法打分）
cl <- cl[sampleID %in% expr_samples]
n2 <- nrow(cl); step("①b 有表达数据     保留 %5d（丢弃 %5d）", n2, n1 - n2)

## ③ 分期 ∈ {I,II,III,IV}（罗马数字解析，先取罗马前缀再白名单）
rom <- toupper(trimws(sub("^\\s*stage\\s*", "", cl$pathologic_stage, ignore.case = TRUE)))
rom <- sub("^([IVX]+).*$", "\\1", rom)
rom[!rom %in% c("I", "II", "III", "IV")] <- NA_character_
cl[, stage_num := match(rom, c("I", "II", "III", "IV"))]
n3 <- nrow(cl[!is.na(stage_num)]); cl <- cl[!is.na(stage_num)]
step("③ 分期 ∈ I–IV    保留 %5d（丢弃 %5d）", n3, n2 - n3)

## ④ 生存时间 > 0（DECEASED 用 days_to_death，其余用 days_to_last_followup）
## 🔴 本表 vital_status 的取值是 **"DECEASED" / "LIVING"**（不是 "Dead"）⇒ 必须两个词都收，
##    否则 184 例死亡会被静默判成存活（首次跑就是这个 bug，见文件末尾登记）。
cl[, days_to_death := suppressWarnings(as.numeric(days_to_death))]
cl[, days_to_last_followup := suppressWarnings(as.numeric(days_to_last_followup))]
cl[, os_event := as.integer(tolower(trimws(vital_status)) %in% c("dead", "deceased"))]
cl[, os_time := fifelse(os_event == 1L, days_to_death, days_to_last_followup)]
n4 <- nrow(cl[!is.na(os_time) & os_time > 0]); cl <- cl[!is.na(os_time) & os_time > 0]
step("④ 生存时间 > 0    保留 %5d（丢弃 %5d）", n4, n3 - n4)

## ⑤ 年龄 / 性别缺失的丢
cl[, age := suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis))]
cl[, sex := toupper(trimws(gender))]
n5 <- nrow(cl[!is.na(age) & !is.na(sex) & sex %in% c("MALE", "FEMALE")])
cl <- cl[!is.na(age) & !is.na(sex) & sex %in% c("MALE", "FEMALE")]
step("⑤ 年龄/性别非缺   保留 %5d（丢弃 %5d）", n5, n4 - n5)

## ② 每患者只留一个样本（按 sampleID 字典序取第一个 —— **不看表达量/生存期**）
setorder(cl, sampleID)
dup <- cl[, .N, by = bcr_patient_barcode][N > 1]
cl  <- cl[, .SD[1], by = bcr_patient_barcode]      # 已按 sampleID 升序 ⇒ 取字典序第一个
setorder(cl, bcr_patient_barcode)
n6 <- nrow(cl)
step("② 每患者一例      保留 %5d（丢弃 %5d 重复；%d 位患者原本多例）", n6, n5 - n6, nrow(dup))
if (nrow(dup)) step("   多例患者：%s", paste(sprintf("%s(%d)", dup$bcr_patient_barcode, dup$N),
                                             collapse = ", "))

## ——— ⑥ 报数 + 落盘冻结，然后**停** ———
cnt <- data.table(
  stage = paste0("Stage ", c("I", "II", "III", "IV")),
  n = as.integer(table(factor(cl$stage_num, levels = 1:4))))
ev <- cl[, .(n = .N, n_death = sum(os_event == 1L),
             os_med_days = median(os_time)), by = .(stage = cl$stage_num)]
setorder(ev, stage)

step("==== ⑥ 清洗结果（**先落数，后跑生存**）====")
step("  保留患者/样本数 = %d", nrow(cl))
for (i in seq_len(nrow(cnt))) step("  %-11s n=%3d", cnt$stage[i], cnt$n[i])
step("  生存事件（死亡）合计 = %d / %d；随访中位 = %.0f 天",
     sum(cl$os_event == 1L), nrow(cl), median(cl$os_time))
step("  事件率按分期：%s", paste(sprintf("%s %.0f%%", paste0("S", ev$stage),
                                       100 * ev$n_death / ev$n), collapse = "  "))
## 守卫：TCGA-LUAD 不可能 0 死亡 ⇒ 事件解析错了必须硬停（首次跑就撞过一次）
if (sum(cl$os_event == 1L) == 0L)
  stop("死亡事件数 = 0，TCGA-LUAD 不可能 ⇒ vital_status 解析有误，硬停（别把结果放出去）",
       call. = FALSE)
step("  年龄中位 %.0f（范围 %.0f–%.0f）；女 %d / 男 %d",
     median(cl$age), min(cl$age), max(cl$age),
     sum(cl$sex == "FEMALE"), sum(cl$sex == "MALE"))

out <- cl[, .(sampleID, patient = bcr_patient_barcode, stage_num,
              stage_raw = pathologic_stage, os_time_days = os_time, os_event,
              age, sex, vital_status)]
fwrite(out, file.path(OUT, "tcga_clean_sample_list.tsv"), sep = "\t")

cnt_long <- data.table(
  step = c("起点（临床表全部行）", "① 原发瘤", "①b 有表达数据", "③ 分期 I–IV",
           "④ 生存时间>0", "⑤ 年龄/性别非缺", "② 每患者一例（终）"),
  n_kept    = c(n0, n1, n2, n3, n4, n5, n6),
  n_dropped = c(NA, n0 - n1, n1 - n2, n2 - n3, n3 - n4, n4 - n5, n5 - n6))
cnt_long[, rule := c("", "①", "①b", "③", "④", "⑤", "②")]
fwrite(cnt_long, file.path(OUT, "tcga_clean_counts.tsv"), sep = "\t")

man <- list(
  step = "M7a TCGA-LUAD 探针样本清洗（S11）",
  archive = ARCH, read_only = TRUE, products_written_to = OUT,
  sha256 = list(HiSeqV2 = sub(" .*$", "", sha_expr),
                LUAD_clinicalMatrix = sub(" .*$", "", sha_clin)),
  n_expr_sample = length(expr_samples), n_clin_row = n0, n_final = nrow(cl),
  rules = c("①=Primary Tumor", "②=每患者一例(sampleID 字典序)", "③=stage∈I..IV",
            "④=os_time>0", "⑤=age/sex 非缺", "⑥=先落数后跑生存"),
  stage_counts = as.list(setNames(cnt$n, cnt$stage)),
  n_death = sum(cl$os_event == 1L),
  os_med_days = median(cl$os_time),
  scale = "HiSeqV2 = log2(RSEM+1)（Xena）",
  frozen_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  next_step = "**本脚本到此为止**：打分(ssGSEA)与 Cox 必须另起脚本，且须核对本表后才跑")
writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE),
           file.path(OUT, "tcga_clean_manifest.json"))

step("  落盘：%s", file.path(OUT, "tcga_clean_sample_list.tsv"))
step("  🔴 S11⑥ 次序门：本脚本**只清洗**。打分与 Cox 在核对上表之前**不得**开跑（法则 3.2）")
step("==== 07 结束 ====")
