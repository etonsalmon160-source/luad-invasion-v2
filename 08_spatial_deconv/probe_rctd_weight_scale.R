## probe_rctd_weight_scale.R —— 探针：RCTD `full` 模式权重的方向、量纲、深度不变性
##
## 为什么有这个脚本
## ----------------
## SPATIAL_CNV_PREREG.md §2.2 的门控口径写成「原始上皮权重 > 0.5」，
## 前提是「原始权重本身就在约等于比例的量纲上、且不随测序深度漂移」。
## 这个前提**不能靠读文档得出**（spacexr 的 Rd 只说 full 模式权重"不必和为 1"，
## 没说是什么量纲）。故用合成数据实测，供 RCTD_PREREG.md §11 引用、G5 硬门挂靠。
##
## 实测结论（2026-09-25，spacexr 2.2.1，doublet_mode='full'，max_cores=1）
## ---------------------------------------------------------------
## 1. `results$weights` 是 **spot × 细胞类型**（行名 = spot 条码，列名 = 类型）。
##    ⇒ `normalize_weights()` 内部 `sweep(weights, 1, rowSums(weights), "/")`
##      按行归一是**正确**的：它把每个 spot 除以自己的和。
## 2. 纯类型 spot 的自类型权重 **0.88–0.99**；逐 spot 行和 **0.85–0.99**。
##    ⇒ 权重已在「约等于比例」的量纲上，不是无标度的抽象数。
## 3. 同一批 spot 计数整体 ×2（构成一字不改），逐 spot 权重几乎不动：
##    0.972→0.963 / 0.956→0.942 / 0.957→0.948 / 0.969→0.966 / 0.930→0.920
##    ⇒ 权重**不随测序深度漂移**。
##
## 🔴 这是**合成数据**。真实切片的对应校验是硬门 G5（RCTD_PREREG.md §4），
##    跑完 56 张即自动判。若 G5 越界，本脚本的结论不成立 ⇒ §2.2 ③ 须回来重签。
##
## 跑法：Rscript 08_spatial_deconv/probe_rctd_weight_scale.R
## 开销：可忽略（3 类型 × 10 spot；即使是队列级任务的同期也不影响）
suppressMessages(library(spacexr))
set.seed(1)

mk_profile <- function(n_gene = 40, base = 0.5, spike_idx) {
  p <- rep(base, n_gene)
  p[spike_idx] <- 20
  p / sum(p)
}

n_gene <- 40
P <- list(A = mk_profile(n_gene, spike_idx = 1:8),
          B = mk_profile(n_gene, spike_idx = 9:16),
          C = mk_profile(n_gene, spike_idx = 17:24))

## 参考：每型 40 个细胞（CELL_MIN_INSTANCE=25 的下限之上）
ref_counts <- do.call(cbind, lapply(names(P), function(k)
  matrix(rpois(n_gene * 40, lambda = 500 * P[[k]]), nrow = n_gene,
         dimnames = list(paste0("G", seq_len(n_gene)), paste0(k, "_", seq_len(40))))))
ref_types <- factor(rep(names(P), each = 40))
names(ref_types) <- colnames(ref_counts)   # Reference() 要求 cell_types 带条码名

## 10 个 spot：前 5 个纯 A、后 5 个纯 B
spot_type <- rep(c("A", "B"), each = 5)
spot_counts <- sapply(spot_type, function(k) rpois(n_gene, lambda = 1000 * P[[k]]))
rownames(spot_counts) <- paste0("G", seq_len(n_gene))
colnames(spot_counts) <- paste0("spot", seq_len(ncol(spot_counts)))
coords <- data.frame(x = rep(1:5, 2), y = rep(1:2, each = 5),
                     row.names = colnames(spot_counts))

run_it <- function(counts, tag) {
  ref <- Reference(ref_counts, ref_types)
  puck <- SpatialRNA(coords, counts, colSums(counts))
  obj <- create.RCTD(puck, ref, max_cores = 1, UMI_min = 100,
                     CELL_MIN_INSTANCE = 25, test_mode = FALSE)
  obj <- run.RCTD(obj, doublet_mode = "full")
  wm <- as.matrix(obj@results$weights)
  cat("\n#### ", tag, " ####\n", sep = "")
  cat("class:", class(obj@results$weights), " dim:", paste(dim(wm), collapse = " x "), "\n")
  cat("rownames(=spot):", paste(head(rownames(wm), 5), collapse = ","), "\n")
  cat("colnames(=type):", paste(colnames(wm), collapse = ","), "\n")
  cat("rowSums:", paste(round(rowSums(wm), 3), collapse = ", "), "\n")
  cat("逐 spot 自类型权重（spot1-5 应≈A，spot6-10 应≈B）:\n")
  cat("  A:", paste(round(wm[, "A"], 3), collapse = ", "), "\n")
  cat("  B:", paste(round(wm[, "B"], 3), collapse = ", "), "\n")
  invisible(wm)
}

w1 <- run_it(spot_counts, "原始深度（每 spot ~1000 UMI）")
w2 <- run_it(spot_counts * 2, "深度翻倍（×2 UMI，构成完全相同）")

cat("\n#### 深度是否改变原始权重 ####\n")
cat("自类型权重｜原始深度:", paste(round(c(w1["spot1", "A"], w1["spot6", "B"]), 3),
                                  collapse = ", "), "\n")
cat("自类型权重｜深度×2  :", paste(round(c(w2["spot1", "A"], w2["spot6", "B"]), 3),
                                  collapse = ", "), "\n")
cat("行和｜原始深度:", paste(round(rowSums(w1), 3), collapse = ", "), "\n")
cat("行和｜深度×2  :", paste(round(rowSums(w2), 3), collapse = ", "), "\n")
