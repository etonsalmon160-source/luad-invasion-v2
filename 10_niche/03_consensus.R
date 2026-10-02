#!/usr/bin/env Rscript
# 03_consensus.R —— M6 §3.4 参数裁定 + §3.5 跨切片共识（域表达画像聚类）+ §3.5 RCTD 对齐敏感性臂
#
# 输入 = 02_banksy_grid.R 逐切片产物：
#   results/10_niche/banksy/<slide>/stats_<agf>.tsv      逐档稳定性统计
#   results/10_niche/banksy/<slide>/domains_<agf>.tsv.gz 逐 spot 域分配（含 5 个种子）
#   results/10_niche/banksy/<slide>/profiles_<agf>.tsv.gz 逐域平均表达画像（HVG 3000 维）
#
# 输出 = results/10_niche/consensus/<agf>/
#   selection.tsv / selection_report.txt   §3.4 逐档存活表 + 选中的 (kg*, λ*, r*)
#   consensus_domains.tsv                  逐 (slide, domain) → 共识生态位类型
#   consensus_profiles.tsv                 共识类型 × HVG 平均画像
#   k_stability.tsv                        K 候选的稳定性（§3.5 第 3 步「K 由稳定性规则定」）
#   sensitivity_rctd_ari.tsv               RCTD 对齐臂 vs 表达臂的一致性
#
# 🔴 **口径已于 2026-10-02 按 §13（S4-REV / S6-REV）重签**：
#   ① λ / k_geom **固定为官方推荐值 0.2 / 18**，不再由网格选（`--lam` / `--kg` 可覆盖，仅供对照）；
#   ② 候选门从「`跨种子 ARI >= 0.90`」换成 **C3 空间连贯性硬门**：
#      逐切片 `coh_mean > null_coh_q95`（置换零分布 95 分位）⇒ 通过；
#      跨切片通过比例须 ≥ FRAC（默认 0.90；`--rule strict_all` 则要求全过）；
#   ③ 存活档按 **C1 的跨种子 ARI 中位** (`ari_mean_med`) 降序；破平「分差 <0.01 取较低 resolution」；
#   ④ **C2 PAC 只报数值、不设硬线**（写在 pac.tsv）；C4 逐切片 ARI 绝对值照报；
#   ⑤ 🔴 **没有存活档时如实报 FAIL 并停止，不得回落到「选 ARI 最高的档」**（§13.0 旧事故根因）。
# ⚠️ 仍留 **补-3（§3.5 的 K 由 `stab_mean >= 0.90` 取最大）** 未改——它是 bootstrap 重采样稳定性、
#    与 S6 被弃的那条「跨种子 ARI」不是同一个量；**已在报告里登记为待裁定**，不私自改。
#
# 跑法：
#   R_LIBS=/home/eto/Rlibs/fastcnv:... Rscript 10_niche/03_consensus.R --agf BOTH --rule frac_ge

suppressMessages({ library(data.table); library(Matrix) })

.a <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d = NULL) { i <- match(k, .a); if (is.na(i)) return(d)
  if (i == length(.a)) stop("参数 ", k, " 后缺值", call. = FALSE); .a[i + 1] }

ROOT  <- "/home/eto/luad_v2"
NICHE <- file.path(ROOT, "results/10_niche")
BANK  <- file.path(NICHE, "banksy")
AGF_MODE <- toupper(get_arg("--agf", "BOTH"))
RULE      <- get_arg("--rule", "frac_ge")     # C3 跨切片聚合：frac_ge（≥FRAC 通过）/ strict_all（全过）
FRAC      <- as.numeric(get_arg("--frac", "0.90"))       # C3 通过比例阈（§13.6）
MIN_ARI   <- as.numeric(get_arg("--min-ari", "0.90"))    # ⚠️ 仅作 C4 报告参考，**不再当门**（§13 S6-REV）
LAM_MAIN  <- as.numeric(get_arg("--lam", "0.2"))         # [§13 S4-REV] 固定官方推荐值
KG_MAIN   <- as.integer(get_arg("--kg", "18"))           # [§13 S4-REV] 固定官方推荐值
RES_TIE   <- 0.01                                        # [已签] 分差<0.01 取较低分辨率
K_GRID    <- as.integer(strsplit(get_arg("--k-grid", "3,4,5,6,8,10,12,16"), ",")[[1]])
K_STAB_MIN<- as.numeric(get_arg("--k-stab", "0.90"))     # ⚠️ 补-3：K 的 bootstrap 稳定性阈，**未改、已登记待裁定**
N_BOOT    <- as.integer(get_arg("--n-boot", "60"))
SEED      <- 20261001L
set.seed(SEED)

step <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                      sprintf(fmt, ...))); flush(stdout())
ir <- function(x) as.integer(round(x))   # sprintf %d 只吃整数型 double
if (!RULE %in% c("strict_all", "frac_ge")) stop("未知 --rule（C3 聚合只支持 strict_all / frac_ge）: ", RULE, call. = FALSE)
AGFS <- if (AGF_MODE == "BOTH") c("agfT", "agfF") else
        if (AGF_MODE == "TRUE") "agfT" else
        if (AGF_MODE == "FALSE") "agfF" else
        stop("--agf 只能是 TRUE / FALSE / BOTH", call. = FALSE)

read_weights <- function(slide) {   # §3.6 用；RCTD 39 亚型权重（逐 spot）
  f <- file.path(ROOT, "results/08_spatial_deconv/rctd_d/per_slide",
                 paste0(slide, ".weights.tsv.gz"))
  if (!file.exists(f)) return(NULL)
  DT <- fread(f); bc <- DT[[1]]; M <- as.matrix(DT[, -1, with = FALSE])  # 首列 = barcode
  rownames(M) <- bc; M
}
sl_subdirs <- function() {
  if (!dir.exists(BANK)) stop("找不到 ", BANK, "（先跑 02_banksy_grid.R）", call. = FALSE)
  list.dirs(BANK, recursive = FALSE, full.names = FALSE)
}

## —————————————————————————————————————————————————————————————
## A. §3.4 逐档汇总 + 裁定
## —————————————————————————————————————————————————————————————
pick_config <- function(agf, slides) {
  S <- rbindlist(lapply(slides, function(s) {
    f <- file.path(BANK, s, sprintf("stats_%s.tsv", agf))
    if (!file.exists(f)) return(NULL)
    d <- fread(f); d[, slide := s]; d }))
  if (!nrow(S)) stop("没有 ", agf, " 的 stats 文件（先跑 02_banksy_grid.R）", call. = FALSE)
  setorder(S, k_geom, lambda, resolution)
  ## —— C3 逐切片：连贯性是否**显著高于该切片自己的标签置换零分布**（> 95 分位 = 单侧 α=0.05）——
  S[, pass3 := coh_mean > null_coh_q95]
  G <- S[, .(n_slide = .N,
             frac_pass3   = mean(pass3, na.rm = TRUE),      # C3 通过比例
             coh_med      = median(coh_mean, na.rm = TRUE),
             cohz_med     = median(coh_z, na.rm = TRUE),
             ari_mean_med = median(ari_mean, na.rm = TRUE), # C1 中位（跨切片）
             ari_mean_q25 = q_safe(ari_mean, .25),
             ari_mean_q75 = q_safe(ari_mean, .75),
             ari_min_med  = median(ari_min, na.rm = TRUE),  # C1 最小侧
             frac_ari_ref = mean(ari_min >= MIN_ARI, na.rm = TRUE),  # C4 报告：旧阈下通过率
             ndom_med     = median(n_domain)),
         by = .(k_geom, lambda, resolution)]
  ## —— C3 硬门（唯一门）——
  G[, survives := switch(RULE,
      strict_all = frac_pass3 >= 1,
      frac_ge    = frac_pass3 >= FRAC)]
  ## —— §13 S4-REV：主候选固定为官方推荐档（λ=0.2 / k_geom=18）——
  G[, main_lam := (lambda == LAM_MAIN)]
  G[, main_cfg := (lambda == LAM_MAIN & k_geom == KG_MAIN)]
  cand <- G[survives & main_cfg]
  ## 🔴 无存活档 ⇒ 如实 FAIL、**不回落到"选 ARI 最高的档"**（§13.0 旧事故根因）
  failed <- !nrow(cand)
  best <- NULL
  if (!failed) {
    ## 选择：C1 的跨种子 ARI 中位降序 → 破平「分差 <0.01 取较低 resolution」
    band <- cand[abs(ari_mean_med - max(ari_mean_med)) < RES_TIE]
    best <- band[which.min(resolution)]
  }
  setorder(cand, -ari_mean_med, resolution)
  list(all = G, ranked = cand, chosen = best, failed = failed, raw = S)
}

## 安全的分位数/中位（整列 NA 时返回 NA 而不是报错）
q_safe <- function(x, p) { x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_); unname(quantile(x, p)) }

## —————————————————————————————————————————————————————————————
## A2. C2：PAC（共识矩阵在 [0.1,0.9] 窗口的增量）—— **只报数值，不设硬线**
##     思路：C(i,j) = (# 种子中 i,j 同域) / H。PAC = 落在 (0.1,0.9) 的 spot 对占比。
##     用「按恰好 k 个坐标相同的对数」N_k 精确反演，避免 n×n 稠密矩阵（n≈1.4 万会爆内存）。
##     B_k = Σ_{|T|=k} (# 在 T 上全同的对)；B_k = Σ_{j>=k} C(j,k) N_j ⇒ 二项反演得 N_k。
## —————————————————————————————————————————————————————————————
pac_from_labels <- function(L) {          # L: n × H 整数矩阵（逐 spot × 种子）
  H <- ncol(L); n <- nrow(L); tot <- n * (n - 1) / 2
  if (H < 2 || tot <= 0) return(list(pac = NA_real_, Nj = rep(NA_real_, max(H, 1)), H = H))
  pairs_same <- function(cols) {          # 在这些列上全同的无序对数
    key <- do.call(paste, c(lapply(cols, function(cc) L[, cc]), sep = "\r"))
    t <- table(key); s <- as.numeric(t); sum(s * (s - 1) / 2)
  }
  Bk <- numeric(H + 1); Bk[1] <- tot
  for (k in 1:H) {
    idxs <- if (k == H) list(1:H) else combn(H, k, simplify = FALSE)
    Bk[k + 1] <- sum(vapply(idxs, pairs_same, numeric(1)))
  }
  Nj <- vapply(1:H, function(j) sum(vapply(j:H, function(k)
        (-1)^(k - j) * choose(k, j) * Bk[k + 1], numeric(1))), numeric(1))
  inwin <- which((1:H) / H > 0.1 & (1:H) / H < 0.9)
  pac <- if (length(inwin)) sum(Nj[inwin]) / tot else NA_real_
  list(pac = pac, Nj = Nj, H = H)
}

## —————————————————————————————————————————————————————————————
## B. §3.5 共识：域表达画像 → 相关 → 层次聚类 → K
## —————————————————————————————————————————————————————————————
consensus_from_profiles <- function(P) {
  # P: 行 = 域（带键），列 = HVG；只用表达（用户裁定 ③）
  key_cols  <- c("slide", "k_geom", "lambda", "resolution", "domain")
  gene_cols <- setdiff(names(P), key_cols)
  keys <- P[, ..key_cols]
  X <- as.matrix(P[, ..gene_cols])
  X <- X / sqrt(rowSums(X^2))                       # 余弦归一（等价于相关聚类）
  X[!is.finite(X)] <- 0
  D <- 1 - tcrossprod(X)                            # 1 - 余弦相似
  diag(D) <- 0
  hc <- hclust(as.dist(D), method = "average")
  ## 🔴 K 必须 ≤ 自助子样本的对象数，否则 cutree 崩（全量域数多不会触发；防的是小样本/异常切片）
  KMAX  <- floor(0.8 * nrow(X))
  KG_USE <- K_GRID[K_GRID <= KMAX]
  if (!length(KG_USE)) stop("域数太少（", nrow(X), "），K 稳定性选择做不了", call. = FALSE)
  if (length(KG_USE) < length(K_GRID))
    step("  ⚠️ 域数 %d 偏少，K 网格降为 {%s}（原 {%s}）", nrow(X),
         paste(KG_USE, collapse = ","), paste(K_GRID, collapse = ","))
  K_stab <- rbindlist(lapply(KG_USE, function(K) {
    cut0 <- cutree(hc, K)
    aris <- replicate(N_BOOT, {
      idx <- sample.int(nrow(X), size = floor(0.8 * nrow(X)), replace = FALSE)
      h <- hclust(as.dist(D[idx, idx]), method = "average")
      mclust::adjustedRandIndex(cut0[idx], cutree(h, K))
    })
    data.table(K = K, stab_mean = mean(aris, na.rm = TRUE),
               stab_q05 = if (any(is.finite(aris))) unname(quantile(aris, 0.05, na.rm = TRUE)) else NA_real_,
               n_ari_ok = sum(is.finite(aris)))
  }))
  ## 自助子样本里两划分可能都塌成一类 ⇒ ARI = NaN（0/0）⇒ 先剔除再平均（分母 n_ari_ok 已落盘备查）
  k_ok <- K_stab[is.finite(stab_mean) & stab_mean >= K_STAB_MIN]
  K_sel <- if (nrow(k_ok)) max(k_ok$K) else {
    fin <- K_stab[is.finite(stab_mean)]
    if (!nrow(fin)) stop("所有 K 的稳定性都是 NA ⇒ 共识做不了", call. = FALSE)
    fin[which.max(stab_mean)]$K }
  keys[, consensus_type := cutree(hc, K_sel)]
  list(keys = keys, K_stab = K_stab, K = K_sel, hc = hc,
       gene_cols = gene_cols,
       consensus_profile = rowsum(as.matrix(P[, ..gene_cols]), group = keys$consensus_type,
                                  reorder = TRUE))
}

## —————————————————————————————————————————————————————————————
## 主流程
## —————————————————————————————————————————————————————————————
slides <- sl_subdirs()
step("发现 %d 张切片的网格产物", length(slides))
for (agf in AGFS) {
  OUT <- file.path(NICHE, "consensus", agf); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  step("==== AGF=%s：§13 S4-REV/S6-REV 裁定（C3 硬门 rule=%s, frac=%.2f；主档 λ=%.1f k=%d）====",
       agf, RULE, FRAC, LAM_MAIN, KG_MAIN)
  pc <- pick_config(agf, slides)
  fwrite(pc$all[order(-survives, -ari_mean_med)], file.path(OUT, "selection.tsv"), sep = "\t")
  ch <- pc$chosen
  ## 🔴 C3 无存活档 ⇒ 如实报 FAIL，不产出 consensus_domains.tsv（下游会自然跳过），**不回落选 ARI 最高档**
  if (pc$failed) {
    hdr <- c(sprintf("AGF = %s", agf),
             sprintf("裁定 = **FAIL**：主档 (λ=%.1f, k_geom=%d) 下没有任何 resolution 过 C3 硬门", LAM_MAIN, KG_MAIN),
             sprintf("C3 门 = 逐切片 coh_mean > null_coh_q95（置换零分布 95 分位）；跨切片通过比例 %s%s",
                     if (RULE == "strict_all") ">= 1.00" else sprintf(">= %.2f", FRAC), ""),
             sprintf("主档各 resolution 的 C3 通过率：%s",
                     paste(sprintf("res=%.1f→%.2f", pc$all[main_cfg == TRUE]$resolution,
                                   pc$all[main_cfg == TRUE]$frac_pass3), collapse = "  ")),
             "🔴 依 §13.6：**不静默回落**。请人工裁定（放宽 C3 聚合 / 换 resolution 网格 / 换参数）。")
    writeLines(hdr, file.path(OUT, "selection_report.txt"))
    step("  🔴 %s 主档未过 C3 硬门 ⇒ FAIL，不产出共识（报告 %s）", agf, file.path(OUT, "selection_report.txt"))
    next
  }
  step("选中：k_geom=%d  lambda=%.1f  resolution=%.1f（C3 通过率 %.2f；C1 跨种子 ARI 中位 %.3f；域数中位 %d）",
       ir(ch$k_geom), ch$lambda, ch$resolution, ch$frac_pass3, ch$ari_mean_med, ir(ch$ndom_med))
  writeLines(c(
    sprintf("AGF = %s", agf),
    sprintf("口径 = §13 S4-REV/S6-REV（2026-10-02 签）；主档固定 λ=%.1f / k_geom=%d", LAM_MAIN, KG_MAIN),
    sprintf("C3 硬门 = 逐切片 coh_mean > null_coh_q95；跨切片通过比例 %s（rule=%s）",
            if (RULE == "strict_all") ">= 1.00" else sprintf(">= %.2f", FRAC), RULE),
    sprintf("选择 = 存活档里按 C1 跨种子 ARI 中位降序；破平「分差 <0.01 取较低 resolution」"),
    sprintf("选中 k_geom=%d lambda=%.1f resolution=%.1f", ir(ch$k_geom), ch$lambda, ch$resolution),
    sprintf("选中档 C3 通过率 %.2f；C1 跨种子 ARI：中位 %.3f（IQR %.3f–%.3f，最小侧中位 %.3f）",
            ch$frac_pass3, ch$ari_mean_med, ch$ari_mean_q25, ch$ari_mean_q75, ch$ari_min_med),
    sprintf("C4 报告（旧「ARI>=%.2f」阈下通过率，**仅报告不再是门**）：%.2f", MIN_ARI, ch$frac_ari_ref),
    sprintf("存活档 %d / 全部档 %d；主候选限定 (λ=%.1f, k_geom=%d)",
            sum(pc$all$survives), nrow(pc$all), LAM_MAIN, KG_MAIN),
    "⚠️ 补-3（§3.5 的 K 由 stab_mean>=0.90 取最大）**未改**：它是 bootstrap 重采样稳定性，与 S6 被弃的量不同，登记待裁定。"),
    file.path(OUT, "selection_report.txt"))

  ## —— 取选中配置的域 + 画像 ——
  keys <- rbindlist(lapply(slides, function(s) {
    f <- file.path(BANK, s, sprintf("profiles_%s.tsv.gz", agf)); if (!file.exists(f)) return(NULL)
    d <- fread(f)[k_geom == ch$k_geom & lambda == ch$lambda & resolution == ch$resolution]
    if (!nrow(d)) return(NULL); d[, slide := s] }), fill = TRUE)
  if (!nrow(keys)) { step("  ！选中配置没有画像行，跳过"); next }
  P <- keys
  step("  画像表：%d 个域 x %d 基因（来自 %d 张切片）", nrow(P), ncol(P) - 6, uniqueN(P$slide))

  ## —— C2：PAC（共识矩阵在 [0.1,0.9] 窗口的增量）—— **只报数值，不设硬线**（§13.2）——
  PAC <- rbindlist(lapply(slides, function(s) {
    f <- file.path(BANK, s, sprintf("domains_%s.tsv.gz", agf)); if (!file.exists(f)) return(NULL)
    D <- fread(f)[k_geom == ch$k_geom & lambda == ch$lambda & resolution == ch$resolution]
    if (!nrow(D) || uniqueN(D$seed) < 2) return(NULL)
    W <- dcast(D, barcode ~ seed, value.var = "domain")
    r <- pac_from_labels(as.matrix(W[, -1, with = FALSE]))
    data.table(slide = s, n_spot = nrow(W), H = r$H, pac = r$pac)
  }), fill = TRUE)
  if (nrow(PAC)) {
    fwrite(PAC, file.path(OUT, "pac.tsv"), sep = "\t")
    step("  C2 PAC（报告项，**不设线**）：切片中位 %.3f（IQR %.3f–%.3f，n=%d）",
         median(PAC$pac, na.rm = TRUE), q_safe(PAC$pac, .25), q_safe(PAC$pac, .75), nrow(PAC))
  } else step("  ⚠️ 取不到逐种子域 ⇒ PAC 跳过")

  ## —— RCTD 组成画像（敏感性臂）：同一批 (slide, domain) 上的 39 亚型平均权重 ——
  W <- rbindlist(lapply(unique(P$slide), function(s) {
    M <- read_weights(s); if (is.null(M)) return(NULL)
    f <- file.path(BANK, s, sprintf("domains_%s.tsv.gz", agf))
    D <- fread(f)[k_geom == ch$k_geom & lambda == ch$lambda & resolution == ch$resolution &
                  seed == 0L, .(barcode, domain)]
    j <- match(D$barcode, rownames(M)); ok <- !is.na(j)
    if (!any(ok)) return(NULL)
    M <- M[j[ok], , drop = FALSE]; D <- D[ok]
    R <- rowsum(M, group = D$domain, reorder = TRUE)
    R <- R / pmax(rowSums(R), 1e-12)
    dt <- as.data.table(R); dt[, domain := as.integer(rownames(R))]; dt[, slide := s]; dt }), fill = TRUE)

  ## —— 表达臂共识 ——
  C <- consensus_from_profiles(P)
  fwrite(C$K_stab, file.path(OUT, "k_stability.tsv"), sep = "\t")
  fwrite(C$keys, file.path(OUT, "consensus_domains.tsv"), sep = "\t")
  CP <- as.data.table(C$consensus_profile)
  CP[, consensus_type := as.integer(rownames(C$consensus_profile))]
  setcolorder(CP, c("consensus_type", C$gene_cols))
  fwrite(CP, file.path(OUT, "consensus_profiles.tsv"), sep = "\t")
  step("  §3.5 共识：K* = %d（稳定性阈 %.2f）；逐张域数中位 %d ⇒ 共 %d 个域归到 %d 类",
       C$K, K_STAB_MIN, ir(median(C$keys[, .N, by = slide]$N)), nrow(C$keys), C$K)
  print(C$K_stab)

  ## —— RCTD 对齐敏感性臂（只在有权重的域上比）——
  if (nrow(W)) {
    k2 <- merge(C$keys[, .(slide, domain, consensus_type)], W, by = c("slide", "domain"))
    if (nrow(k2) > 2) {
      Rm <- as.matrix(k2[, -(1:3), with = FALSE]); Rm[!is.finite(Rm)] <- 0
      Rn <- Rm / sqrt(rowSums(Rm^2)); Rn[!is.finite(Rn)] <- 0
      Dr <- 1 - tcrossprod(Rn); diag(Dr) <- 0
      hcr <- hclust(as.dist(Dr), method = "average")
      cr <- cutree(hcr, C$K)
      ari <- mclust::adjustedRandIndex(k2$consensus_type, cr)
      fwrite(data.table(n_domain = nrow(k2), K = C$K, ari_expr_vs_rctd = ari),
             file.path(OUT, "sensitivity_rctd_ari.tsv"), sep = "\t")
      step("  §3.5 敏感性臂（RCTD 组成对齐）：ARI(表达臂 vs RCTD 臂) = %.3f ⇒ %s",
           ari, if (ari >= 0.8) "两臂一致，结论稳" else "**不一致 ⇒ 须如实报「生态位定义对方法敏感」**")
    }
  } else step("  ⚠️ 没有拿到 RCTD 权重，敏感性臂跳过")
}
step("==== 03 结束 ====")
