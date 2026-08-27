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

draw_Figure6 <- function() {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  f6r_draw_main(workspace_root = workspace_root)
}

if (sys.nframe() == 0L) draw_Figure6()
