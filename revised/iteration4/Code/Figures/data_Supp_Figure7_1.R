#!/usr/bin/env Rscript

Sys.setenv(
  KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
  } else {
    normalizePath(file.path(getwd(), "Code", "Figures"), mustWork = TRUE)
  }
})
source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))

data_Supp_Figure7_1 <- function(n_core = 8L, rebuild = FALSE) {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  paths <- f7r_paths(workspace_root)
  required <- c(
    file.path(paths$figure7, "response_class_smoothed_curves.tsv"),
    file.path(paths$figure7, "response_class_curve_class_by_seed.tsv"),
    file.path(
      paths$figure5,
      "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
    )
  )
  f7r_require_files(required, "Supplementary Figure 7-1 analytical input")
  invitro <- f7x_compute_separate_invitro(
    paths, n_core = as.integer(n_core), rebuild = isTRUE(rebuild)
  )
  invisible(list(required = required, invitro = invitro))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1L]])
  }
  data_Supp_Figure7_1(
    n_core = as.integer(value("n-core", "8")),
    rebuild = tolower(value("rebuild", "false")) %in% c("true", "1", "yes")
  )
}
