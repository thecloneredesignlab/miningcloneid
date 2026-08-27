#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
  } else {
    normalizePath(file.path(getwd(), "Code", "Figures"), mustWork = TRUE)
  }
})
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))

data_Supp_Figure6_1 <- function() {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  paths <- f6r_paths(workspace_root)
  required <- c(
    file.path(paths$figure6, "response_class_smoothed_curves.tsv"),
    file.path(paths$figure6, "response_class_curve_class_by_seed.tsv"),
    file.path(
      paths$figure5,
      "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
    )
  )
  f6r_require_files(required, "Supplementary Figure 6-1 analytical input")
  invisible(required)
}

if (sys.nframe() == 0L) data_Supp_Figure6_1()
