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
source(file.path(script_dir, "util", "analysis", "figure6_context_extension.R"))
source(file.path(script_dir, "util", "analysis", "figure6_finite_time_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_invitro_passage_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_extended_time_o2.R"))
source(file.path(script_dir, "util", "analysis", "figure6_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_net_growth_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_net_growth_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_supplementary_b_layout.R"))
source(file.path(script_dir, "util", "analysis", "figure6_finite_time_plots.R"))

draw_Figure6 <- function() {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  f6ft_draw_main(workspace_root = workspace_root)
}

if (sys.nframe() == 0L) draw_Figure6()
