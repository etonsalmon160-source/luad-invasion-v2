#!/usr/bin/env Rscript
## 21_depth_balance_report.R —— 深度配平审计（只读；不动任何已签口径的产物）
##
## 为什么做：借锚臂全队列结果是「LUAD 观测上皮的 cnv_fraction 高于 AIS/MIA」。
## 但 LUAD 组织本身测得更深（pooled 中位 8,765 vs 4,268/3,819 nUMI）⇒ 这个差**可能是深度造的**。
## margin 梯度管的是「并列边缘」，**管不了深度**。本脚本专门查这一条。
##
## 设计（**新口径，跑前须登记**）：
##   深度量：**nUMI**（每个 spot 的总计数）为主 —— 它是 cnv_fraction 的直接驱动量；
##           nFeature 作次要复核（同样分析再跑一遍，两套并列报）。
##   S0 深度缺口：pooled 与逐患者报 LUAD vs 前驱期(AIS/MIA) 观测上皮的深度中位与比值。
##   S1 患者内深度五分位：在**每位患者自己的**观测上皮上切五档，每档内算
##      median(cf[LUAD]) − median(cf[前驱])；跨患者取中位。⇒ 差值是否在各档内都成立。
##   S2 共同支撑：逐患者取两期别深度 [5%,95%] 的交集，只留落在里面的 spot，重算同一对比。
##   S3 随机子样对照（**关键对照**）：S2 会同时「配平深度」和「减少 spot 数」两件事。
##      为把两者分开，每位患者**随机**抽与 S2 相同数量的 LUAD / 前驱 spot（不看深度），
##      重复 N_DRAW 次取中位。⇒ 若 S3 的差 ≈ 未配平的差，则 S2 的缩水来自**深度**；
##      若 S3 也缩水，则只是样本量。
##
## 🔴 本脚本只上报。cnv_fraction 不是恶性判定；AAH/Normal 无观测切片，不得读成五档病程轴。
##
## 跑法：
##   Rscript --vanilla 08_spatial_deconv/21_depth_balance_report.R

RES  <- "/home/eto/luad_v2/results/08_spatial_deconv"
BASE <- file.path(RES, "spatial_cnv/cohort_borrow")
QC   <- file.path(RES, "spatial_qc_per_spot.csv.gz")
OUT  <- file.path(BASE, "depth_balance")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

N_BIN    <- 5      # 深度五分位
MIN_SPOT <- 10     # 档内每期别至少这么多 spot 才计入
N_DRAW   <- 500    # S3 随机子样重复次数（**同期别内**抽，只换样本量）
N_BOOT   <- 5000   # S4 患者自助法重抽次数（**重抽患者**，重复单元是患者不是 spot）
SEED     <- 1
PRE      <- c("AIS", "MIA")
set.seed(SEED)

stage_of <- function(s) sub("-\\d+$", "", sub("^GSM[0-9]+_P[0-9]+_", "", s))

## ——— 深度查表 ———
qc  <- read.csv(gzfile(QC), stringsAsFactors = FALSE)
key <- paste(qc$slide, qc$barcode, sep = "|")
depth_of <- function(slide, bc) {
  m <- match(paste(slide, bc, sep = "|"), key)
  list(nUMI = qc$nUMI[m], nFeature = qc$nFeature[m])
}

## ——— 读逐患者观测上皮 + 深度 ———
PATIENTS <- sprintf("P%d", 1:25)
obs <- list(); miss <- character(0)
for (P in PATIENTS) {
  mf <- file.path(BASE, P, "manifest.json")
  if (!file.exists(mf)) { miss <- c(miss, P); next }
  ## 主臂文件名随 main_kind 变（main_self_ / main_borrow_）⇒ 直接找，不解析 json
  fs <- list.files(file.path(BASE, P), pattern = "^main_.*_spot_scores\\.rds$", full.names = TRUE)
  if (length(fs) != 1L) { miss <- c(miss, P); next }
  o <- readRDS(fs[1])
  i <- which(o$grp == "observation")
  if (!length(i)) { miss <- c(miss, P); next }
  raw <- mapply(function(s, b) sub(paste0("^", s, "_"), "", b),
                o$slide[i], o$bc[i], USE.NAMES = FALSE)
  d <- depth_of(o$slide[i], raw)
  obs[[P]] <- data.frame(patient = P, slide = o$slide[i], stage = stage_of(o$slide[i]),
                         cf = o$cf[i], nUMI = d$nUMI, nFeature = d$nFeature,
                         row.names = NULL)
}
cat(sprintf("读到 %d 例；缺 %s\n", length(obs), if (length(miss)) paste(miss, collapse = ",") else "无"))
allobs <- do.call(rbind, obs)
cat(sprintf("观测上皮 spot 合计 %d；深度命中 %d（%.4f%%）\n",
            nrow(allobs), sum(!is.na(allobs$nUMI)), 100 * mean(!is.na(allobs$nUMI))))

## ——— S0：深度缺口 ———
s0_pool <- do.call(rbind, lapply(sort(unique(allobs$stage)), function(s) {
  v <- allobs[allobs$stage == s, ]
  data.frame(scope = "pooled", key = s, n_spot = nrow(v),
             nUMI_med = median(v$nUMI, na.rm = TRUE),
             nFeature_med = median(v$nFeature, na.rm = TRUE), row.names = NULL)
}))
s0_pat <- do.call(rbind, lapply(names(obs), function(P) {
  d <- obs[[P]]
  lu <- d$nUMI[d$stage == "LUAD"]; pr <- d$nUMI[d$stage %in% PRE]
  if (!length(lu) || !length(pr)) return(NULL)
  data.frame(scope = "patient", key = P, n_spot = length(lu) + length(pr),
             nUMI_med = median(lu) / median(pr),          # 这里放**比值**
             nFeature_med = median(d$nFeature[d$stage == "LUAD"]) /
                            median(d$nFeature[d$stage %in% PRE]),
             row.names = NULL)
}))
write.table(rbind(s0_pool, s0_pat), file.path(OUT, "S0_depth_gap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ——— 统一列结构的行构造器（🔴 三臂必须同列，否则 rbind 报 columns do not match）———
mkrow <- function(P, metric, arm, nL, nP, n_used, n_tot, kept, mlu, mpre, delta,
                  lo95 = NA_real_, hi95 = NA_real_) {
  data.frame(patient = P, metric = metric, arm = arm,
             n_LUAD = nL, n_pre = nP, n_used = n_used, n_total = n_tot,
             frac_kept = kept, med_LUAD = mlu, med_pre = mpre, delta = delta,
             lo95 = lo95, hi95 = hi95, row.names = NULL)
}

## ——— S1 / S2 / S3 主循环 ———
run_metric <- function(metric) {
  bin_rows <- list(); arm_rows <- list()
  for (P in names(obs)) {
    d <- obs[[P]]
    if (!any(d$stage == "LUAD") || !any(d$stage %in% PRE)) next
    x <- d[[metric]]; cf <- d$cf
    iL <- which(d$stage == "LUAD"); iP <- which(d$stage %in% PRE)

    ## 未配平（同一患者、全部 spot）：作为对照线
    arm_rows[[length(arm_rows) + 1L]] <- mkrow(
      P, metric, "unrestricted", length(iL), length(iP), length(iL) + length(iP),
      length(x), 1, median(cf[iL]), median(cf[iP]), median(cf[iL]) - median(cf[iP]))

    ## S1：患者内深度五分位
    if (length(unique(x[!is.na(x)])) >= N_BIN) {
      br <- unique(quantile(x, seq(0, 1, length.out = N_BIN + 1), na.rm = TRUE))
      if (length(br) >= 3) {
        br[1] <- -Inf; br[length(br)] <- Inf
        g <- cut(x, breaks = br, labels = FALSE, include.lowest = TRUE)
        for (b in seq_len(max(g, na.rm = TRUE))) {
          k <- which(g == b)
          if (sum(iL %in% k) < MIN_SPOT || sum(iP %in% k) < MIN_SPOT) next
          kL <- intersect(k, iL); kP <- intersect(k, iP)
          bin_rows[[length(bin_rows) + 1L]] <- data.frame(
            patient = P, metric = metric, bin = b, n_LUAD = length(kL), n_pre = length(kP),
            med_LUAD = median(cf[kL]), med_pre = median(cf[kP]),
            delta = median(cf[kL]) - median(cf[kP]), row.names = NULL)
        }
      }
    }

    ## S2：共同支撑（两期别深度 [5%,95%] 的交集）
    lo <- max(quantile(x[iL], .05, na.rm = TRUE), quantile(x[iP], .05, na.rm = TRUE))
    hi <- min(quantile(x[iL], .95, na.rm = TRUE), quantile(x[iP], .95, na.rm = TRUE))
    k <- which(x >= lo & x <= hi)
    kL <- intersect(k, iL); kP <- intersect(k, iP)
    if (length(kL) >= MIN_SPOT && length(kP) >= MIN_SPOT) {
      arm_rows[[length(arm_rows) + 1L]] <- mkrow(
        P, metric, "common_support", length(kL), length(kP), length(k),
        length(x), length(k) / length(x),
        median(cf[kL]), median(cf[kP]), median(cf[kL]) - median(cf[kP]))
      ## S3：随机抽**同样多**的 LUAD / 前驱 spot（**同期别内**抽、不看深度）⇒ 只留样本量效应
      dl <- replicate(N_DRAW, {
        median(cf[sample(iL, length(kL))]) - median(cf[sample(iP, length(kP))])
      })
      arm_rows[[length(arm_rows) + 1L]] <- mkrow(
        P, metric, "random_subsample", length(kL), length(kP), length(k),
        length(x), length(k) / length(x), median(dl), NA_real_, median(dl),
        quantile(dl, .05), quantile(dl, .95))

      ## S4（患者级检验）在汇总段做，见下。此处不再逐 spot 打乱标签 ——
      ## 🔴 原「逐 spot 置换期别标签」的零分布**退化**（登记为缺陷，见 §20.6）：
      ##    每组抽 ~8000 spot 求中位数，永远逼近总体中位数，两组相减恒为 0
      ##    （P21 实测 sd=0.0021、80% 恰为 0、全分布只有 4 个取值）⇒ 该检验分辨力≈0，
      ##    会给假 p=0。重复单元是**患者**不是 spot，改走符号检验 + 患者自助法。
    }
  }
  list(bin = do.call(rbind, bin_rows), arm = do.call(rbind, arm_rows))
}

res <- lapply(c("nUMI", "nFeature"), run_metric)
names(res) <- c("nUMI", "nFeature")

write.table(res$nUMI$bin, file.path(OUT, "S1_bins_by_patient.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ——— 汇总 ———
ARMS <- c("unrestricted", "random_subsample", "common_support")
summ <- function(r) {
  bins <- do.call(rbind, lapply(sort(unique(r$bin$bin)), function(x) {
    z <- r$bin[r$bin$bin == x, ]
    data.frame(bin = x, n_patient = nrow(z), n_spot = sum(z$n_LUAD + z$n_pre),
               delta_median = median(z$delta), n_pt_pos = sum(z$delta > 0),
               n_pt_neg = sum(z$delta < 0), row.names = NULL)
  }))
  p <- r$arm
  ## 只保留**三臂都齐**的患者，才可三向对比
  keep <- Reduce(intersect, lapply(ARMS, function(a) p$patient[p$arm == a]))
  pp <- p[p$patient %in% keep, ]
  cmp <- do.call(rbind, lapply(ARMS, function(a) {
    z <- pp[pp$arm == a, ]
    data.frame(arm = a, n_patient = nrow(z), delta_median = median(z$delta),
               n_pt_pos = sum(z$delta > 0), n_pt_neg = sum(z$delta < 0),
               frac_spot_kept_med = median(z$frac_kept), row.names = NULL)
  }))
  list(bins = bins, cmp = cmp, arms = pp)
}
S <- lapply(res, summ)

## ——— S4：**患者级**检验 —— 符号检验（精确二项）+ 患者自助法 CI ———
## 为什么不用逐 spot 置换：见 run_metric 内的登记。重复单元是患者。
s4_summary <- function(r) {
  do.call(rbind, lapply(c("unrestricted", "common_support"), function(a) {
    z <- r$arm[r$arm$arm == a, ]; d <- z$delta   # 🔴 字段名是 arm（单数）；写 r$arms 会 NULL 静默匹配
    nz <- d[d != 0]; n <- length(nz); npos <- sum(nz > 0)
    ## 符号检验（双侧精确二项）：零假设「正负各半」
    p_sign <- if (n == 0) NA_real_ else min(1, 2 * pbinom(min(npos, n - npos), n, .5))
    ## 患者自助法：重抽**患者**（保留患者间差异），给队列中位数的 95% CI
    bs <- replicate(N_BOOT, median(d[sample.int(length(d), length(d), replace = TRUE)]))
    data.frame(
      arm = a, level = if (a == "common_support") "depth_matched" else "all_spots",
      n_patient = length(d), n_signed = n, n_pt_pos = npos, n_pt_neg = sum(nz < 0),
      delta_median = median(d),
      boot_lo2.5 = as.numeric(quantile(bs, .025)),
      boot_hi97.5 = as.numeric(quantile(bs, .975)),
      p_sign_two_sided = p_sign, row.names = NULL)
  }))
}
s4 <- lapply(res, function(r) s4_summary(r))
write.table(do.call(rbind, lapply(names(s4), function(m) cbind(metric = m, s4[[m]]))),
            file.path(OUT, "S4_patient_level_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ——— S5：配平前后**逐期别** pooled 中位（只含同时有两期别的患者）———
pat_paired <- names(obs)[vapply(names(obs), function(P) {
  d <- obs[[P]]; any(d$stage == "LUAD") && any(d$stage %in% PRE) }, logical(1))]
s5 <- do.call(rbind, lapply(c("unrestricted", "depth_matched"), function(arm) {
  do.call(rbind, lapply(sort(unique(allobs$stage)), function(s) {
    v <- do.call(c, lapply(pat_paired, function(P) {
      d <- obs[[P]]
      if (arm == "unrestricted") return(d$cf[d$stage == s])
      x <- d$nUMI; iL <- which(d$stage == "LUAD"); iP <- which(d$stage %in% PRE)
      lo <- max(quantile(x[iL], .05, na.rm = TRUE), quantile(x[iP], .05, na.rm = TRUE))
      hi <- min(quantile(x[iL], .95, na.rm = TRUE), quantile(x[iP], .95, na.rm = TRUE))
      d$cf[d$stage == s & x >= lo & x <= hi]
    }))
    data.frame(arm = arm, stage = s, n_spot = length(v),
               median = if (length(v)) median(v) else NA_real_, row.names = NULL)
  }))
}))
write.table(s5, file.path(OUT, "S5_stage_median_matched.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(S$nUMI$bins, file.path(OUT, "S1_bins_summary_nUMI.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(S$nFeature$bins, file.path(OUT, "S1_bins_summary_nFeature.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, res), file.path(OUT, "S2S3_paired_by_patient.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(cbind(metric = "nUMI", S$nUMI$cmp),
                  cbind(metric = "nFeature", S$nFeature$cmp)),
            file.path(OUT, "S2S3_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ——— 图（English only：本机无 CJK 字体）———
png(file.path(OUT, "depth_balance.png"), width = 1700, height = 640, res = 130)
op <- par(mfrow = c(1, 3), mar = c(4.6, 5.0, 3.2, 1), las = 1)
boxplot(nUMI ~ stage, data = allobs[allobs$stage %in% c("AIS", "MIA", "LUAD"), ],
        log = "y", col = "grey80", outline = FALSE,
        main = "Observation epithelium depth", xlab = "stage", ylab = "nUMI (log10 scale)")
bb <- S$nUMI$bins
plot(bb$bin, bb$delta_median, type = "b", pch = 19, col = "firebrick", lwd = 2,
     ylim = range(c(0, res$nUMI$bin$delta), na.rm = TRUE),
     xlab = "within-patient depth quintile (1=shallowest)",
     ylab = "median(LUAD) - median(AIS/MIA)", main = "Depth-stratified")
for (p in unique(res$nUMI$bin$patient)) {
  z <- res$nUMI$bin[res$nUMI$bin$patient == p, ]
  lines(z$bin, z$delta, col = "grey70", lty = 3)
}
lines(bb$bin, bb$delta_median, type = "b", pch = 19, col = "firebrick", lwd = 2)
abline(h = 0, lty = 2, col = "grey40")
cc <- S$nUMI$cmp
bp <- barplot(cc$delta_median, names.arg = c("all\nspots", "random\nsubsample", "depth\nmatched"),
              col = c("steelblue", "grey60", "darkorange"), main = "Headline (nUMI)",
              ylab = "median paired delta")
abline(h = 0, lty = 2, col = "grey40")
text(bp, cc$delta_median, labels = sprintf("%d/%d pt\n+", cc$n_pt_pos, cc$n_patient),
     pos = 3, cex = 0.75)
par(op); dev.off()

cat("\n==== S0 深度缺口（pooled）====\n"); print(s0_pool)
cat(sprintf("\n逐患者 LUAD/前驱 深度比（nUMI）：中位 %.2f（范围 %.2f–%.2f）\n",
            median(s0_pat$nUMI_med), min(s0_pat$nUMI_med), max(s0_pat$nUMI_med)))
cat(sprintf("深度比 >1 的患者：%d/%d\n", sum(s0_pat$nUMI_med > 1), nrow(s0_pat)))
cat("\n==== S1 深度五档（nUMI）====\n"); print(bb)
cat("\n==== S1 深度五档（nFeature）====\n"); print(S$nFeature$bins)
cat("\n==== S2/S3 三向对比 ====\n"); print(S$nUMI$cmp); print(S$nFeature$cmp)
cat("\n==== S4 患者级检验（符号检验 + 患者自助法）====\n"); print(s4$nUMI); print(s4$nFeature)
cat("\n==== S5 配平前后逐期别 pooled 中位 ====\n"); print(s5)
cat(sprintf("\n随机子样 %d 次；患者自助法 %d 次；seed=%d\n", N_DRAW, N_BOOT, SEED))
cat(sprintf("\n产物：%s\n", OUT))
