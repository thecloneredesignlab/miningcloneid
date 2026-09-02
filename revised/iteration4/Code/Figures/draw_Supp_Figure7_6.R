#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("This entry point must be run with Rscript.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
})
source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))
source(file.path(script_dir, "util", "analysis", "figure7_finite_time_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure7_finite_time_plots.R"))

if (sys.nframe() == 0L) {
  f7ft_draw_supplement_7_6(
    normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  )
}
