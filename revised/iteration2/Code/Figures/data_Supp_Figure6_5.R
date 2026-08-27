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

source(file.path(script_dir, "util", "runtime", "process_runner.R"))
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure6_context_extension.R"))
source(file.path(
  script_dir, "util", "analysis", "figure6_finite_time_consistency.R"
))

data_Supp_Figure6_5 <- function(n_core = 4L) {
  s65_data(
    workspace_root = WORKSPACE_ROOT,
    invitro_result_root = INVITRO_RESULT_ROOT,
    n_core = as.integer(n_core)
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- args[startsWith(args, "--n-core=")]
  if (length(hit) > 1L) stop("--n-core was provided more than once.")
  n_core <- if (length(hit)) {
    as.integer(sub("^--n-core=", "", hit[[1L]]))
  } else {
    4L
  }
  if (!is.finite(n_core) || n_core < 1L) {
    stop("--n-core must be a positive integer.")
  }
  data_Supp_Figure6_5(n_core = n_core)
}
