#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Supp_Figure4_1 <- function() {
  data_dir <- file.path(DATA_ROOT, "Figure4")
  panel_dir <- file.path(data_dir, "panels")
  outputs <- file.path(OUTPUT_ROOT, c(
    "supp_fig4-1_all18_cluster_prior_violins.png",
    "supp_fig4-1_all18_cluster_prior_violins.pdf"
  ))
  run_figure4_continuous_association(data_dir = data_dir)
  render_parameter_landscape(
    data_dir = data_dir,
    panel_dir = panel_dir,
    output_dir = OUTPUT_ROOT
  )
  staged <- file.path(
    panel_dir,
    c(
      "all_parameter_fitted_endpoint_distributions.png",
      "all_parameter_fitted_endpoint_distributions.pdf"
    )
  )
  require_files(staged, "Supplementary Figure 4-1 staged output")
  for (i in seq_along(staged)) {
    ok <- file.copy(staged[[i]], outputs[[i]], overwrite = TRUE)
    if (!ok) {
      stop("Failed to publish Supplementary Figure 4-1: ", staged[[i]])
    }
  }
  require_files(outputs, "Supplementary Figure 4-1 output")
  invisible(outputs[[1L]])
}

if (sys.nframe() == 0L) draw_Supp_Figure4_1()
