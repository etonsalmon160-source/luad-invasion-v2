#!/usr/bin/env Rscript
# 38_domain_forest_cox.R —— 七个生态位域 marker 签名的预后森林图数据
#   口径：每个域用 **top30 marker** → 基因跨样本 z → 取均值 → 按中位二分 Low/High
#         单变量 Cox，参照组 = Low（即 HR>1 表示 High 分数更差）
#   终点：OS（os_time_days/os_event）＋ DFS（新发肿瘤事件；无事件者用末次随访截尾）
#   数据边界：只读 /home/eto/luad_invasion/data/TCGA_LUAD；产物写 luad_v2
suppressMessages({ library(data.table); library(survival) })

ROOT <- "/home/eto/luad_v2"
ARC  <- "/home/eto/luad_invasion/data/TCGA_LUAD"
M7A  <- file.path(ROOT, "results/10_niche/m7a")
OUT  <- file.path(ROOT, "results/10_niche/m7a")

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())

SL  <- fread(file.path(M7A, "tcga_clean_sample_list.tsv"))
SC  <- fread(file.path(M7A, "m7a_niche_signature_scores.tsv"))
C   <- fread(file.path(ARC, "LUAD_clinicalMatrix"), sep = "\t")
C   <- C[, .(sampleID,
             dfs_event = new_tumor_event_after_initial_treatment,
             dfs_days  = days_to_new_tumor_event_after_initial_treatment,
             fu_days   = days_to_last_followup)]
D <- merge(merge(SL, SC, by = "sampleID"), C, by = "sampleID", all.x = TRUE)
step("合并后样本 %d（OS 事件 %d）", nrow(D), sum(D$os_event == 1, na.rm = TRUE))

## DFS：事件 = YES；时间 = 事件日；无事件 = 末次随访截尾
D[, ev := fifelse(toupper(trimws(dfs_event)) == "YES", 1L,
                  fifelse(toupper(trimws(dfs_event)) == "NO", 0L, NA_integer_))]
D[, tm := fifelse(ev == 1L, as.numeric(dfs_days), as.numeric(fu_days))]
DDFS <- D[!is.na(ev) & !is.na(tm) & tm > 0]
step("DFS 可用样本 %d（事件 %d，中位随访 %.0f 天）", nrow(DDFS), sum(DDFS$ev == 1),
     median(DDFS$tm))

fit_one <- function(d, tm, ev, v, mode = "lv") {
  x <- d[[v]]
  dd <- data.table(time = d[[tm]], event = d[[ev]])
  if (mode == "lv") {
    dd[, g := factor(ifelse(x > median(x, na.rm = TRUE), "High", "Low"),
                     levels = c("Low", "High"))]
    dd <- dd[!is.na(time) & !is.na(event) & time > 0 & !is.na(g)]
    if (length(unique(dd$g)) < 2) return(NULL)
    f <- coxph(Surv(time, event) ~ g, data = dd)
  } else {                                   # 连续：分数已是跨样本 z，故 HR = 每 +1 SD
    dd[, s := as.numeric(scale(x))]
    dd <- dd[!is.na(time) & !is.na(event) & time > 0 & !is.na(s)]
    f <- coxph(Surv(time, event) ~ s, data = dd)
  }
  s <- summary(f)
  data.table(HR = s$conf.int[1, "exp(coef)"], lo = s$conf.int[1, "lower .95"],
             hi = s$conf.int[1, "upper .95"], p = s$coefficients[1, "Pr(>|z|)"],
             n = nrow(dd), n_event = sum(dd$event))
}

res <- rbindlist(lapply(sprintf("D%d", 1:7), function(v) {
  rows <- list()
  for (md in c("lv", "cont")) {
    a <- fit_one(D,    "os_time_days", "os_event", v, md)
    b <- fit_one(DDFS, "tm",           "ev",       v, md)
    rows[[md]] <- data.table(domain = v, mode = md,
      OS_HR = a$HR,  OS_lo = a$lo,  OS_hi = a$hi,  OS_p = a$p,  OS_n = a$n,  OS_ev = a$n_event,
      DFS_HR = b$HR, DFS_lo = b$lo, DFS_hi = b$hi, DFS_p = b$p, DFS_n = b$n, DFS_ev = b$n_event)
  }
  rbindlist(rows)
}))
fwrite(res, file.path(OUT, "m7a_niche_forest_lowhigh.tsv"), sep = "\t")
print(res[, .(domain, mode, OS_HR = round(OS_HR, 2), OS_p = signif(OS_p, 3),
              DFS_HR = round(DFS_HR, 2), DFS_p = signif(DFS_p, 3))])
step("写出 m7a_niche_forest_lowhigh.tsv（含连续 + 二分两种口径）")
