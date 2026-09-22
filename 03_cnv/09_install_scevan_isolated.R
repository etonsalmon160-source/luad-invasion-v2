#!/usr/bin/env Rscript
# SCEVAN 第二臂的**隔离安装**（只装，不跑）。
#
# 环境约束：不得污染项目依赖的 R 库。本机会话前 `SCEVAN`/`yaGST` 均未装，
# 其余依赖（scran / fgsea / ggtree / parallelDist / Rtsne ...）**已全部存在**
# ⇒ 只有两个 GitHub-only 包要装，且**都不在 CRAN**，故"冻结版本"只能是 git SHA。
#
# 装到专用库目录并用 `.libPaths()` 前置加载，主库（含 Seurat/harmony 的项目环境）不受影响。

LIB         <- "/home/eto/Rlibs/SCEVAN"
SCEVAN_REPO <- "AntonioDeFalco/SCEVAN"
SCEVAN_SHA  <- "5a49b88ac9445eeffcebb95404e3190992faac04"   # 2026-03-19
YAGST_REPO  <- "miccec/yaGST"
YAGST_SHA   <- "56227df3ae183070c9d156af11c306ee799435e6"   # 2017-11-02

dir.create(LIB, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.access(LIB, 2) == 0)
.libPaths(c(LIB, .libPaths()))
cat("[LIB] 目标库:", LIB, "\n")
cat("[LIB] libPaths:", paste(.libPaths(), collapse = " | "), "\n")

before <- .libPaths()[-1]

install_pinned <- function(repo, sha, name) {
  if (requireNamespace(name, quietly = TRUE)) {
    cat(sprintf("[SKIP] %s 已存在（版本 %s）\n", name,
                as.character(utils::packageVersion(name))))
    return(invisible(TRUE))
  }
  cat(sprintf("[INSTALL] %s @ %s -> %s\n", repo, substr(sha, 1, 8), LIB))
  remotes::install_github(paste0(repo, "@", sha), lib = LIB,
                          dependencies = FALSE,   # 依赖已齐；不碰既有包
                          upgrade = "never",
                          build_vignettes = FALSE, quiet = FALSE)
  stopifnot(requireNamespace(name, quietly = TRUE))
  invisible(TRUE)
}

# 依赖顺序：yaGST 是 SCEVAN 的 Imports
install_pinned(YAGST_REPO, YAGST_SHA, "yaGST")
install_pinned(SCEVAN_REPO, SCEVAN_SHA, "SCEVAN")

cat("\n===== 核验 =====\n")
for (p in c("yaGST", "SCEVAN")) {
  cat(sprintf("  %-8s %s  %s\n", p,
              requireNamespace(p, quietly = TRUE),
              tryCatch(as.character(utils::packageVersion(p)), error = function(e) "NA")))
}
cat("  SCEVAN 取自:", find.package("SCEVAN"), "\n")
# 🔴 2026-09-22 实测更正：README 写的 `norm_cells` **在本版本不存在**。
#    v1.0.3 的真实参数名是 pipelineCNA(norm_cell=) 与 classifyTumorCells(norm_cell_names=)。
#    若照 README 的字面名去调，会直接报 "unused argument"。此处按**代码里的真名**核验。
p_args <- names(formals(SCEVAN::pipelineCNA))
c_args <- names(formals(SCEVAN::classifyTumorCells))
cat("  pipelineCNA 有 `norm_cell`:", "norm_cell" %in% p_args, "\n")
cat("  classifyTumorCells 有 `norm_cell_names`:",
    "norm_cell_names" %in% c_args, "\n")
cat("  FIXED_NORMAL_CELLS 存在:",
    "FIXED_NORMAL_CELLS" %in% p_args, "\n")
cat("  不存在 `norm_cells`（README 的叫法）:",
    !("norm_cells" %in% c(p_args, c_args)), "\n")

# 主库完好性核对：装完后项目依赖的包仍可加载
cat("\n===== 主库未受影响的核对 =====\n")
for (p in c("Seurat", "harmony", "scran", "copykat")) {
  ok <- requireNamespace(p, quietly = TRUE)
  loc <- if (ok) dirname(find.package(p)) else "-"
  cat(sprintf("  %-9s %s  %s\n", p, ok, loc))
}
cat("\n[DONE]\n")
