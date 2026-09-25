#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else
    normalizePath(file.path(getwd(), "Code", "Figures"), mustWork = TRUE)
})
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure6_context_extension.R"))
source(file.path(script_dir, "util", "analysis", "figure6_finite_time_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_invitro_passage_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_finite_time_plots.R"))
source(file.path(script_dir, "util", "analysis", "figure6_net_growth_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_supplementary_b_layout.R"))

draw_Supp_Figure6_14 <- function() {
  workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  f6sb_draw(
    workspace_root = workspace_root, supplement = "Supp_Figure6_14",
    filename = "supp_fig6-14_population_net_live_growth_rate",
    metric = "mean_net_growth_rate", day_limits = c(0, 1000),
    o2_limits = c(0, 2),
    title = "Population-weighted net-live growth rate"
  )
}

if (sys.nframe() == 0L) draw_Supp_Figure6_14()
