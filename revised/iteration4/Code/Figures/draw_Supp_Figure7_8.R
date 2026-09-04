#!/usr/bin/env Rscript

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
source(file.path(script_dir, "util", "analysis", "figure7_finite_time_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure7_invitro_passage_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure7_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure7_finite_time_plots.R"))

draw_Supp_Figure7_8 <- function() {
  workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  f7ft_draw_supp7_8(workspace_root = workspace_root)
}

if (sys.nframe() == 0L) draw_Supp_Figure7_8()
