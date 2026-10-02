#!/usr/bin/env Rscript
# 04_annotate.R —— M6 §3.6：给共识生态位类型贴「描述性」标签（RCTD 组成，**只描述不判定**）
#
# 输入：
#   results/10_niche/consensus/<agf>/consensus_domains.tsv   只含 03 选中的那一档
#   results/10_niche/banksy/<slide>/domains_<agf>.tsv.gz     逐 spot 域分配（要 barcode）
#   results/08_spatial_deconv/rctd_d/per_slide/<slide>.weights.tsv.gz  39 亚型权重，首列 barcode
#
# 输出（每个 AGF 一套，写进 results/10_niche/consensus/<agf>/）：
#   annotation_<agf>.tsv       共识类型 × 39 亚型：域内中位比例 / 切片全体中位比例 / 富集倍数
#   annotation_top_<agf>.tsv   每个共识类型富集倍数前 N 的亚型
#   domain_names_<agf>.tsv     中性域名（top 富集亚型拼）+ 通过禁令词检查
#   annotation_manifest_<agf>.json
#
# 🔴 §3.6 / §7 禁令：域名只用中性描述，**禁** tumor / malignant / cancer / KAC / invasive。
#   RCTD 权重**绝不当恶性判读**（铁律 2：本项目没有逐细胞恶性标签）。
#
# 富集倍数 = （域内该亚型中位比例）/（同批切片全体 spot 中位比例），切片内归一。
#   >1 相对富集；≈1 与切片平均无异 ⇒ 防"什么都有"的域被误读成特征域。
#
# 跑法：
#   R_LIBS=/home/eto/Rlibs/fastcnv:... Rscript 10_niche/04_annotate.R --agf BOTH

suppressMessages({ library(data.table); library(Matrix) })

.a <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d = NULL) { i <- match(k, .a); if (is.na(i)) return(d)
  if (i == length(.a)) stop("参数 ", k, " 后缺值", call. = FALSE); .a[i + 1] }

ROOT  <- "/home/eto/luad_v2"
NICHE <- file.path(ROOT, "results/10_niche")
CONS  <- file.path(NICHE, "consensus")
RCTD  <- file.path(ROOT, "results/08_spatial_deconv/rctd_d/per_slide")
AGF_MODE <- toupper(get_arg("--agf", "BOTH"))
TOP_N    <- as.integer(get_arg("--top", "5"))

BANNED <- c("tumor", "malignant", "cancer", "kac", "invasive")   # §7 禁令（小写比对）
step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())
AGFS <- if (AGF_MODE == "BOTH") c("agfT", "agfF") else
        if (AGF_MODE == "TRUE") "agfT" else
        if (AGF_MODE == "FALSE") "agfF" else
        stop("--agf 只能是 TRUE / FALSE / BOTH", call. = FALSE)

check_ban <- function(x, what) {   # 命中禁令词就硬停：宁可不产出，也不放违规名出去
  hit <- x[grepl(paste(BANNED, collapse = "|"), tolower(x))]
  if (length(hit)) stop(sprintf("%s 命中 §7 禁令词：%s ⇒ 硬停", what,
                                paste(hit, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}
read_weights <- function(s) {
  f <- file.path(RCTD, paste0(s, ".weights.tsv.gz")); if (!file.exists(f)) return(NULL)
  DT <- fread(f); bc <- DT[[1]]; M <- as.matrix(DT[, -1, with = FALSE]); rownames(M) <- bc; M
}

for (agf in AGFS) {
  dfile <- file.path(CONS, agf, "consensus_domains.tsv")
  if (!file.exists(dfile)) { step("！缺 %s ⇒ 跳过（先跑 03_consensus.R）", dfile); next }
  OUT <- file.path(CONS, agf)
  KD  <- fread(dfile)
  keys <- unique(KD[, .(k_geom, lambda, resolution)])
  if (nrow(keys) != 1L) stop("consensus_domains.tsv 里不止一个配置 ⇒ 不该发生", call. = FALSE)
  KG <- keys$k_geom; LAM <- keys$lambda; RES <- keys$resolution
  step("==== AGF=%s：%d 个域 / %d 类 / %d 张；配置 k_geom=%d λ=%.1f r=%.1f ====",
       agf, nrow(KD), uniqueN(KD$consensus_type), uniqueN(KD$slide),
       as.integer(KG), LAM, RES)

  rows <- list(); bad <- character(0); cols <- NULL
  for (s in unique(KD$slide)) {
    M <- read_weights(s); bf <- file.path(NICHE, "banksy", s, sprintf("domains_%s.tsv.gz", agf))
    if (is.null(M) || !file.exists(bf)) { bad <- c(bad, s); next }
    cols <- colnames(M)
    D <- fread(bf)[k_geom == KG & lambda == LAM & resolution == RES & seed == 0L,
                   .(barcode, domain)]
    j <- match(D$barcode, rownames(M)); ok <- !is.na(j)
    if (!nrow(D) || !any(ok)) { bad <- c(bad, s); next }
    R  <- rowsum(M[j[ok], , drop = FALSE], group = D$domain[ok], reorder = TRUE)
    ## 域内平均组成（该域所有 spot 的权重按行求和后归一）
    Rm <- R / pmax(rowSums(R), 1e-12)
    dt <- as.data.table(Rm)
    dt[, `:=`(domain = as.integer(rownames(R)),
              n_spot = as.integer(table(D$domain[ok])[rownames(R)]), slide = s)]
    rows[[length(rows) + 1L]] <- dt[, c("slide", "domain", "n_spot", cols), with = FALSE]
  }
  if (!length(rows)) { step("！没有任何切片拿到 RCTD 权重 ⇒ 跳过"); next }
  W <- rbindlist(rows, fill = TRUE)
  if (length(bad)) step("  ⚠️ %d 张缺权重/缺域文件，已跳过：%s", length(bad),
                        paste(head(bad, 8), collapse = ", "))

  ## 分母：逐切片全体 spot 平均组成，再跨切片取中位（避免被大切片主导）
  den_rows <- rbindlist(lapply(unique(W$slide), function(s) {
    M <- read_weights(s); if (is.null(M)) return(NULL)
    dt <- as.data.table(as.list(colMeans(M)))
    setnames(dt, cols); dt[, slide := s]; dt }), fill = TRUE)
  den_med <- as.numeric(den_rows[, lapply(.SD, median), .SDcols = cols])

  ## 按共识类型汇总
  W <- merge(W, KD[, .(slide, domain, consensus_type)], by = c("slide", "domain"), all.x = TRUE)
  agg <- rbindlist(lapply(sort(unique(W$consensus_type)), function(ct) {
    Z <- as.matrix(W[consensus_type == ct, ..cols])
    dm <- apply(Z, 2, median, na.rm = TRUE)
    data.table(consensus_type = ct, subtype = cols,
               prop_median_in_domain = as.numeric(dm),
               prop_median_all_spot  = den_med,
               enrichment            = as.numeric(dm) / pmax(den_med, 1e-12),
               n_domain              = nrow(Z))
  }))
  setorder(agg, consensus_type, -enrichment)
  fwrite(agg, file.path(OUT, sprintf("annotation_%s.tsv", agf)), sep = "\t")

  topN <- agg[, head(.SD, TOP_N), by = consensus_type]
  fwrite(topN, file.path(OUT, sprintf("annotation_top_%s.tsv", agf)), sep = "\t")

  check_ban(cols, "RCTD 亚型名")
  nm <- rbindlist(lapply(sort(unique(agg$consensus_type)), function(ct) {
    t1 <- topN[consensus_type == ct]
    a <- t1$subtype[1]; b <- if (nrow(t1) > 1) t1$subtype[2] else NA_character_
    lab <- if (!is.na(b) && t1$enrichment[1] - t1$enrichment[2] < 0.05)
             sprintf("%s/%s-high domain", a, b) else sprintf("%s-high domain", a)
    data.table(consensus_type = ct, domain_name = lab, top1 = a, top1_enrich = t1$enrichment[1],
               top2 = b, top2_enrich = if (nrow(t1) > 1) t1$enrichment[2] else NA_real_) }))
  check_ban(nm$domain_name, "域名")
  fwrite(nm, file.path(OUT, sprintf("domain_names_%s.tsv", agf)), sep = "\t")

  step("  中性域名（已过禁令检查）：")
  for (i in seq_len(nrow(nm)))
    cat(sprintf("    %2d -> %s  (top: %s x%.2f)\n", nm$consensus_type[i], nm$domain_name[i],
                nm$top1[i], nm$top1_enrich[i]))

  man <- list(step = "§3.6 域注释（RCTD 描述性标签）", agf = agf,
              selected = list(k_geom = KG, lambda = LAM, resolution = RES),
              n_consensus_type = nrow(nm), top_n = TOP_N,
              enrichment_def = "域内该亚型中位比例 / 同批切片全体 spot 中位比例（切片内归一）",
              ban_check = "域名与亚型名均过 §7 禁令词检查",
              caveat = "RCTD 权重只作组成描述，绝不当恶性判读（铁律 2）",
              n_slide_used = uniqueN(W$slide), n_slide_skipped = length(bad),
              finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE),
             file.path(OUT, sprintf("annotation_manifest_%s.json", agf)))
  step("  落盘：annotation_* / annotation_top_* / domain_names_* → %s", OUT)
}
step("==== 04 结束 ====")
