#!/usr/bin/env Rscript

Sys.setenv(
  KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
  } else {
    normalizePath(file.path(getwd(), "Code", "Figures"), mustWork = TRUE)
  }
})

source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "si_figure7_eigenmodes.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))
source(file.path(script_dir, "util", "analysis", "figure7_invitro_extended_o2.R"))

data_Supp_Figure7_4 <- function(n_core = 8L, rebuild = FALSE) {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  s64_data(
    workspace_root = workspace_root,
    n_core = as.integer(n_core),
    rebuild = isTRUE(rebuild)
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1L]])
  }
  data_Supp_Figure7_4(
    n_core = as.integer(value("n-core", "8")),
    rebuild = tolower(value("rebuild", "false")) %in% c("true", "1", "yes")
  )
}
