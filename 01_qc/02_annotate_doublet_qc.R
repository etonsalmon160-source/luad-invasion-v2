#!/usr/bin/env Rscript
# =============================================================================
# 01_qc/02_annotate_doublet_qc.R — 标注 scDblFinder 未能建立阈值的样本
#
# 背景：个别样本 scDblFinder 报 "Threshold found: 1" → 判 0 双体。
#       这是其已知失败模式（未找到分离阈值），**不等于真的没有双体**。
#       实测 P7_LUAD 用 dbr=NULL 与固定 dbr=0.075 均无法建立阈值 →
#       **如实标注为 threshold_not_found，不伪造双体率**。
#       该样本质量本身良好（median nCount 1505 / nFeature 1100，均高于队列中位）。
#
# 依据：从主运行日志 `logs/M1_308103_qc.log` 解析（可复现），不硬编码样本名。
# 产出：给 gse308103_qc_per_sample.csv 增加 `doublet_note` 列。
# =============================================================================
OUT <- "/home/eto/luad_v2/results/01_qc"
LOG <- "/home/eto/luad_v2/logs/M1_308103_qc.log"

ps <- read.csv(file.path(OUT, "gse308103_qc_per_sample.csv"), stringsAsFactors = FALSE)

# 从日志解析：紧跟 "Threshold found:1" 之后的那条 "[i/75] <sample> ..." 行
flag <- character(0)
if (file.exists(LOG)) {
  ln <- readLines(LOG, warn = FALSE)
  for (i in seq_along(ln)) {
    if (grepl("^Threshold found:1$", ln[i])) {
      for (j in (i + 1):min(i + 3, length(ln))) {
        m <- regmatches(ln[j], regexpr("\\[\\s*[0-9]+/75\\]\\s+(\\S+)", ln[j]))
        if (length(m) == 1) { flag <- c(flag, sub(".*\\]\\s+(\\S+).*", "\\1", m)); break }
      }
    }
  }
}
ps$doublet_note <- ifelse(ps$sample_id %in% flag,
  "threshold_not_found(scDblFinder 未建立阈值; 0 doublets; 未伪造)", "")
write.csv(ps, file.path(OUT, "gse308103_qc_per_sample.csv"), row.names = FALSE)

cat(sprintf("[annotate] 标记样本数=%d: %s\n", length(flag), paste(flag, collapse = ", ")))
cat(sprintf("[annotate] 总: pre=%d pass=%d doublet=%d (%.2f%%)\n",
            sum(ps$n_pre), sum(ps$n_pass), sum(ps$n_doublet),
            100 * sum(ps$n_doublet) / sum(ps$n_pass)))
