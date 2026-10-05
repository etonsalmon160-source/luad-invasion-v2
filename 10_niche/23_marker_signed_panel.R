#!/usr/bin/env Rscript
# 23_marker_signed_panel.R —— 靶点扰动逆向臂 §4.3：生成**带符号**的 150+150 查询签名
#
# CMap 的查询是**带符号基因集**（Subramanian 2017：「Each gene in the query carries a sign…」），
# 而 `d13_archetype_markers.tsv` 只存了 top30 且**只有上调侧**（`order(-lfc)[1:30]`）。
# 本脚本从 `d13_parts/`（434 域 × 3000 HVG 均值）**重算双向**：
#   q_up   = log2fc 降序 top150
#   q_down = log2fc 升序 top150
#
# 只读 d13_parts，不重跑表达提取。产物：
#   results/10_niche/kstar_diag/d14_signed_panel.tsv

suppressMessages({ library(data.table) })

ROOT  <- "/home/eto/luad_v2"
KD    <- file.path(ROOT, "results/10_niche/kstar_diag")
PARTS <- file.path(KD, "d13_parts")
NTOP  <- 150L

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...))); flush(stdout())

pf <- list.files(PARTS, pattern = "\\.tsv$", full.names = TRUE)
step("读 %d 个切片的分域表达", length(pf))
D <- rbindlist(lapply(pf, fread), fill = TRUE)
A <- fread(file.path(KD, "d13_domain_archetype.tsv"))[, .(slide, domain, archetype)]
D <- merge(D, A, by = c("slide", "domain"))
GENES <- setdiff(names(D), c("slide", "domain", "archetype"))
step("域 %d × 基因 %d", nrow(D), length(GENES))

E <- as.matrix(D[, ..GENES]); eps <- 1e-4
P <- rbindlist(lapply(sort(unique(D$archetype)), function(a) {
  inx <- D$archetype == a
  m1 <- colMeans(E[inx, , drop = FALSE]); m0 <- colMeans(E[!inx, , drop = FALSE])
  lfc <- log2((m1 + eps) / (m0 + eps))
  o   <- order(-lfc)
  up   <- o[seq_len(NTOP)]
  dn   <- order(lfc)[seq_len(NTOP)]
  rbindlist(list(
    data.table(archetype = a, direction = "up",   rank = seq_len(NTOP), gene = GENES[up], log2fc = lfc[up]),
    data.table(archetype = a, direction = "down", rank = seq_len(NTOP), gene = GENES[dn], log2fc = lfc[dn])
  ))
}))
setorder(P, archetype, direction, rank)
fwrite(P, file.path(KD, "d14_signed_panel.tsv"), sep = "\t")

step("落盘 d14_signed_panel.tsv（%d 行）", nrow(P))
cat("\n=== 各签名的 q_up / q_down top10 ===\n")
for (a in sort(unique(P$archetype))) {
  u <- P[archetype == a & direction == "up"][1:10]; d <- P[archetype == a & direction == "down"][1:10]
  cat(sprintf("  D%d  up  : %s\n", a, paste(u$gene, collapse = ", ")))
  cat(sprintf("      down: %s\n", paste(d$gene, collapse = ", ")))
}
cat("\n=== 每列 log2fc 的范围（看是否有足够动态范围）===\n")
print(P[, .(n = .N, min_lfc = round(min(log2fc), 2), max_l2fc = round(max(log2fc), 2)),
        by = .(archetype, direction)][order(archetype, direction)])
step("完成")
