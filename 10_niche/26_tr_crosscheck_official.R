#!/usr/bin/env Rscript
# 26_tr_crosscheck_official.R —— 靶点扰动逆向臂：**用官方实现独立复算我们的尺子**
#
# 被验对象：我们自实现的 WTCS / NCS（生产管线 24_target_reversal.py）。
# 参照实现：`signatureSearch`（Bioconductor 3.16 / 1.12.0，Girke 实验室）的
#           `.enrichScore()` 与 `.lincsScores()` —— CMap/LINCS 的官方 R 实现。
#
# 🔴 为什么必须有这一步：WTCS/NCS/τ 是**我们自己写的**。ES 已与 `fgsea` 对拍，
#    但 NCS/τ 从未被独立实现验证；而 NCS 上确实出过一个**静默产生假阴性**的符号 bug
#    （μ⁻ 未取绝对值 ⇒ NCS 全正 ⇒ "τ≤−90" 恒为空，表面像"无候选"）。
#    官方源码同一处写着 `# without abs() sign of neg values would switch to pos`。
#
# 输入：由 25_tr_dump_for_crosscheck.py 导出的目录
#   X_f64.bin（行主序 float64，nsig×ngene）/ genes.txt / sigs.tsv / query.txt / ours.tsv
# 输出：crosscheck_report.txt + crosscheck_compare.tsv
#
# 口径（与生产管线逐条对齐，出处见 TARGET_REVERSAL_PARAMETERS.md）：
#   type（权重指数 p）= 1      ← GSEA 默认；**CMap 正文未明示**，属工具惯例
#   NCS 分组键        = <cell_id>__<pert_type>，**正/负分别取均值**
#   NCS 分母          = abs(组内均值)   ← 官方注释明示此 abs 必需
#   零 WTCS           = 从均值计算中剔除（同官方：es_na[es_na==0] <- NA）
#   签名命名          = <sig_id>__<cell_id>__<pert_type>（官方靠 gsub("^.*?__","") 取分组键）
#
# 用法：
#   R_LIBS=~/Rlibs/sigsearch:... Rscript 26_tr_crosscheck_official.R <xcheck_dir> [type]

suppressMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("用法: Rscript 26_tr_crosscheck_official.R <xcheck_dir> [type]")
D    <- args[1]
TYPE <- if (length(args) >= 2) as.integer(args[2]) else 1L
stopifnot(TYPE == 1L)   # 口径锁定；改这里必须同时改预注册

log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

## ── 读入 ─────────────────────────────────────────────
genes <- readLines(file.path(D, "genes.txt"))
sigs  <- fread(file.path(D, "sigs.tsv"))
qy    <- strsplit(readLines(file.path(D, "query.txt")), " ")
q_up  <- qy[[1]]; q_dn <- qy[[2]]
nsig  <- nrow(sigs); ngene <- length(genes)
log(sprintf("子集 %d 签名 × %d 基因；q_up %d / q_down %d", nsig, ngene, length(q_up), length(q_dn)))

con <- file(file.path(D, "X_f64.bin"), "rb")
X <- readBin(con, what = "numeric", n = nsig * ngene, size = 8)   # float64
close(con)
X <- matrix(X, nrow = nsig, ncol = ngene, byrow = TRUE)   # 行主序
colnames(X) <- genes
log("矩阵读入完成")

## ── 官方实现 ─────────────────────────────────────────
es_official <- signatureSearch:::.enrichScore
lincsScores <- signatureSearch:::.lincsScores

log("用官方 .enrichScore 复算 WTCS …")
t0 <- Sys.time()
res <- apply(X, 1, function(x) {
  sv <- sort(x, decreasing = TRUE)                       # 官方调用约定：降序
  eup <- es_official(sigvec = sv, Q = q_up, type = TYPE)
  edn <- es_official(sigvec = sv, Q = q_dn, type = TYPE)
  if (sign(eup) != sign(edn)) (eup - edn) / 2 else 0
})
log(sprintf("  完成 %.0f s", as.numeric(Sys.time() - t0, units = "secs")))

## 官方 NCS：名字须为 <任意>__<cell>__<pert_type>
esout <- res
names(esout) <- sprintf("%s__%s__%s", sigs$sig_id, sigs$cell_id, sigs$pert_type)
ls_out <- lincsScores(esout = esout, upset = q_up, downset = q_dn,
                      minTauRefSize = 0L, tau = FALSE)
off <- data.table(sig_id = sigs$sig_id,
                  wtcs_official = as.numeric(ls_out$WTCS),
                  ncs_official  = as.numeric(ls_out$NCS))

## ── 比对 ─────────────────────────────────────────────
ours <- fread(file.path(D, "ours.tsv"))
M <- merge(ours, off, by = "sig_id")
M[, `:=`(d_wtcs = abs(wtcs - wtcs_official), d_ncs = abs(ncs - ncs_official))]

rep_txt <- file.path(D, "crosscheck_report.txt")
w <- function(...) cat(sprintf(...), "\n", sep = "")
sink(rep_txt)
w("=== 靶点扰动逆向臂 · 尺子独立对拍报告 ===")
w("日期: %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("被验: 自实现 WTCS/NCS（24_target_reversal.py）")
w("参照: signatureSearch %s 的 .enrichScore / .lincsScores（官方 CMap/LINCS 实现）",
  as.character(packageVersion("signatureSearch")))
w("口径: type(p) = %d", TYPE)
w("样本: %d 签名（自 %s 子集）", nrow(M), basename(D))
w("")
w("[1] 逐条比对")
w("  WTCS  最大绝对差 %.3e   中位 %.3e   完全一致(<1e-12) %s",
  max(M$d_wtcs), median(M$d_wtcs), all(M$d_wtcs < 1e-12))
w("  NCS   最大绝对差 %.3e   中位 %.3e   完全一致(<1e-12) %s",
  max(M$d_ncs), median(M$d_ncs), all(M$d_ncs < 1e-12))
w("")
w("[2] 符号一致性（这是最关键的一项——曾在此处出过静默 bug）")
w("  WTCS 符号一致 %d/%d ；NCS 符号一致 %d/%d",
  sum(sign(M$wtcs) == sign(M$wtcs_official)), nrow(M),
  sum(sign(M$ncs)  == sign(M$ncs_official)),  nrow(M))
w("")
w("[3] 分布对照")
w("  自实现   WTCS [%.4f, %.4f]  NCS [%.4f, %.4f]",
  min(M$wtcs), max(M$wtcs), min(M$ncs), max(M$ncs))
w("  官方     WTCS [%.4f, %.4f]  NCS [%.4f, %.4f]",
  min(M$wtcs_official), max(M$wtcs_official), min(M$ncs_official), max(M$ncs_official))
w("")
w("结论: %s", if (max(M$d_wtcs) < 1e-12 && max(M$d_ncs) < 1e-12)
  "✅ 自实现与官方实现**逐位一致**，口径无误。" else
  "❌ 存在差异，见 crosscheck_compare.tsv，须逐条查明后才能跑正式版。")
sink()
cat(readLines(rep_txt), sep = "\n")
fwrite(M, file.path(D, "crosscheck_compare.tsv"), sep = "\t")
log(sprintf("产物 → %s", D))
