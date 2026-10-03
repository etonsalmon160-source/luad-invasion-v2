#!/usr/bin/env Rscript
# 16_domain_cnv.R —— M6 §18：把空间 CNV 分值叠到 7 个域类型上
#
# 动机（用户 2026-10-03）：「这没有恶性程度给我一种怪怪的感觉」。
#   域是用 RCTD 组成命名的，而 **RCTD 的参考集里没有恶性细胞类型**
#   ⇒ 组成看不见肿瘤。CNV 是独立于 RCTD/BANKSY 的连续量，用来补这一维。
#
# 🔴 合规（铁律 2 / §11.2 / 恶性×权重禁令）：
#   * CNV 只作**连续量** `cf` 报，**不**给任何 spot / 域贴"恶性"标签；
#   * 本脚本产出的是"该域的 CNV 分值分布"，**不是**恶性判定；
#   * 不做 恶性 × 权重 之类的相乘。
#
# 两段式：
#   --slide <name>  重现该切片 5 种子共识划分 → 逐 spot 的 domain → 映射到 d7 的 archetype
#                   → 原子写 d8_parts/<slide>.tsv（barcode, domain, archetype）
#   --combine       并 CNV D 臂 spot 分数 → 逐 archetype 的 cf 分布 + 图

suppressMessages({ library(data.table); library(ggplot2) })

ROOT  <- "/home/eto/luad_v2"
BN    <- file.path(ROOT, "results/10_niche/banksy")
CNVB  <- file.path(ROOT, "results/08_spatial_deconv/spatial_cnv/cohort_borrow")
OUT   <- file.path(ROOT, "results/10_niche/kstar_diag")
PARTS <- file.path(OUT, "d8_parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
MAIN  <- list(k_geom = 18L, lambda = 0.2, resolution = 0.5)
CNV_ARM <- "D_arm_spot_scores.rds"   # 25/25 例、56/56 张

.a <- commandArgs(trailingOnly = TRUE)
ga <- function(k) { i <- match(k, .a); if (is.na(i) || i == length(.a)) NULL else .a[i + 1] }
SLIDE <- ga("--slide"); COMBINE <- "--combine" %in% .a
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())

one_slide <- function(s) {
  d <- fread(file.path(BN, s, "domains_agfT.tsv.gz"), sep = "\t")
  d <- d[k_geom == MAIN$k_geom & lambda == MAIN$lambda & resolution == MAIN$resolution]
  W <- dcast(d, barcode ~ seed, value.var = "domain"); bc <- W$barcode
  L <- as.matrix(W[, -1L, with = FALSE])
  C <- matrix(0, length(bc), length(bc))
  for (j in seq_len(ncol(L))) { O <- outer(L[, j], L[, j], "=="); C <- C + O; rm(O) }
  C <- C / ncol(L)
  K <- as.integer(median(apply(L, 2, function(x) length(unique(x)))))
  dom <- cutree(hclust(as.dist(1 - C), method = "average"), k = K)   # 确定性，可复现
  data.table(slide = s, barcode = bc, domain = as.integer(dom))
}

if (!is.null(SLIDE)) {
  r <- one_slide(SLIDE)
  tmp <- file.path(PARTS, sprintf("%s.tsv.tmp", SLIDE)); fwrite(r, tmp, sep = "\t")
  if (!file.rename(tmp, file.path(PARTS, sprintf("%s.tsv", SLIDE)))) stop("原子改名失败", call. = FALSE)
  step("落盘 %s（%d spot）", SLIDE, nrow(r)); quit(save = "no", status = 0)
}
if (!COMBINE) stop("需要 --slide 或 --combine", call. = FALSE)

## ————————————— 汇总 —————————————
sl <- list.files(BN, pattern = "^GSM")
pf <- file.path(PARTS, sprintf("%s.tsv", sl)); have <- file.exists(pf)
step("D8 汇总：应有 %d / 已有 %d", length(sl), sum(have))
if (any(!have)) step("  ! 缺：%s", paste(sl[!have], collapse = " "))
S <- rbindlist(lapply(pf[have], fread), fill = TRUE)
A <- fread(file.path(OUT, "d7_domain_assign.tsv"))[, .(slide, domain, archetype)]
S <- merge(S, A, by = c("slide", "domain"), all.x = TRUE)
if (anyNA(S$archetype)) step("  ! %d 个 spot 未映射到 archetype", sum(is.na(S$archetype)))
S[, key := paste0(slide, "_", barcode)]

## CNV：D 臂（逐患者文件）
pat <- list.files(CNVB, pattern = "^P[0-9]+$", full.names = TRUE)
CV <- rbindlist(lapply(pat, function(p) {
  f <- file.path(p, CNV_ARM); if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  ## 🔴 列名不能叫 `key`：data.table() 的 key= 是保留参数（设主键），不是列名（我踩过）
  data.table(bc_key = as.character(x$bc), cf = as.numeric(x$cf),
             slide = as.character(x$slide),
             grp = if (is.null(x$grp)) NA_character_ else as.character(x$grp))
}), fill = TRUE)
step("CNV spot 数 %d（%d 张切片，臂=%s）", nrow(CV), uniqueN(CV$slide), CNV_ARM)

M <- merge(S, CV[, .(bc_key, cf, grp)], by.x = "key", by.y = "bc_key")
step("配对成功 %d / %d（%.1f%%）", nrow(M), nrow(S), 100 * nrow(M) / nrow(S))
fwrite(M[, .(slide, barcode, domain, archetype, cf, grp)],
       file.path(OUT, "d8_spot_archetype_cnv.tsv"), sep = "\t", compress = "gzip")

dg <- fread(file.path(ROOT, "results/10_niche/depth_guard/summary_per_slide.tsv"))
dg <- dg[agf == "agfT", .(depth_med = median(depth_med)), by = slide]
K <- M[, .(n_spot = .N, cf_mean = mean(cf), cf_med = median(cf),
           cf_q75 = quantile(cf, .75), cf_ge_10 = mean(cf >= 0.10),
           n_slide = uniqueN(slide)), by = archetype]
## 逐 archetype 的深度：用该域所有 spot 所属切片的深度中位
dz <- merge(M[, .(slide, archetype)], dg, by = "slide")[, .(depth_med = median(depth_med)), by = archetype]
K <- merge(K[, .(archetype, n_spot, cf_mean, cf_med, cf_q75, cf_ge_10, n_slide)], dz, by = "archetype")
setorder(K, archetype)
fwrite(K, file.path(OUT, "d8_archetype_cnv.tsv"), sep = "\t")

## 基线对照：reference vs observation spot 的 cf
B <- M[, .(n = .N, cf_med = median(cf)), by = grp]

con <- file(file.path(OUT, "d8_summary.txt"), "wt")
w <- function(...) cat(sprintf(...), file = con, append = TRUE)
w("=== §18 域类型 × 空间 CNV（臂=%s）===\n%s\n\n", CNV_ARM, format(Sys.time(), "%Y-%m-%d %H:%M"))
w("🔴 CNV 只作**连续量**报，不是恶性判定（铁律 2 / §11.2）\n\n")
w("spot 配对 %d；切片 %d\n\n", nrow(M), uniqueN(M$slide))
w("[1] 逐域类型的 CNV 分值\n")
print(as.data.frame(K), file = con, digits = 4)
w("\n[2] 参考 spot vs 观测 spot 的基线（该臂自带区分）\n")
print(as.data.frame(B), file = con, digits = 4)
w("\n[3] 逐切片：域类型中位 cf 的排序是否与深度同向\n")
per <- M[, .(cf = median(cf), depth = NA_real_), by = .(slide, archetype)]
per <- merge(per, dg, by = "slide")
print(per[, .(n_slide = uniqueN(slide), cf_med = round(median(cf), 4),
              depth_med = round(median(depth_med))), by = archetype][order(-cf_med)],
      file = con, digits = 4)
close(con); cat(readLines(file.path(OUT, "d8_summary.txt")), sep = "\n")

LBL <- c("1 airway epithelium","2 peribronchovascular stroma","3 immune infiltrate",
         "4 alveolar wall","5 AT2 epithelium","6 vessel-wall smooth muscle","7 lymphoid aggregate")
K[, lab := factor(sprintf("D%d %s", archetype, LBL[archetype]), levels = sprintf("D%d %s", archetype, LBL[archetype]))]
g <- ggplot(K, aes(lab, cf_med)) +
  geom_col(aes(fill = depth_med), colour = "grey25") +
  geom_text(aes(label = sprintf("%.3f\nn=%d", cf_med, n_spot)), vjust = -0.2, size = 2.9) +
  scale_fill_viridis_c(option = "magma", name = "median counts\nper spot") +
  labs(x = NULL, y = "median spatial-CNV fraction (per spot)",
       title = "§18: does each domain carry a CNV signal?",
       subtitle = "bar height = median cf of the archetype's spots; fill = depth of the slides it comes from") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(OUT, "fig_d8_archetype_cnv.png"), g, width = 9.5, height = 5.2, dpi = 150)
step("产物 → %s", OUT)
