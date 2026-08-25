#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))
source(file.path(script_dir, "util", "analysis", "si_figure6_eigenmodes.R"))

data_Supp_Figure6_3 <- function(n_core = 8L, rebuild = FALSE) {
  workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  si6_data(
    workspace_root = workspace_root,
    n_core = as.integer(n_core), rebuild = isTRUE(rebuild)
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1L]])
  }
  data_Supp_Figure6_3(
    n_core = as.integer(value("n-core", "8")),
    rebuild = tolower(value("rebuild", "false")) %in% c("true", "1", "yes")
  )
}
