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
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))

data_Figure7 <- function(
    n_core = 8L, rebuild = FALSE, n_resample = 100L
) {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  invivo <- f7r_data(
    workspace_root = workspace_root,
    n_core = as.integer(n_core),
    rebuild = isTRUE(rebuild),
    n_resample = as.integer(n_resample)
  )
  invitro <- f7x_data(
    workspace_root = workspace_root,
    n_core = as.integer(n_core),
    rebuild = isTRUE(rebuild)
  )
  invisible(list(invivo = invivo, invitro = invitro))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  option_value <- function(name, default = character()) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (!length(hit)) return(default)
    sub(paste0("^--", name, "="), "", hit[[1L]])
  }
  as_flag <- function(x) {
    tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
  }
  data_Figure7(
    n_core = as.integer(option_value("n-core", "8")),
    rebuild = as_flag(option_value("rebuild", "FALSE")),
    n_resample = as.integer(option_value("n-resample", "100"))
  )
}
