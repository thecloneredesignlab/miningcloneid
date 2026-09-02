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

draw_Supp_Figure7_1 <- function() {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  f7x_draw_supplement_7_1(workspace_root = workspace_root)
}

if (sys.nframe() == 0L) draw_Supp_Figure7_1()
