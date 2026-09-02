#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "si_figure7_eigenmodes.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))

draw_Supp_Figure7_3 <- function() {
  workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  si7_draw_weak_gap(workspace_root = workspace_root)
}

if (sys.nframe() == 0L) draw_Supp_Figure7_3()
