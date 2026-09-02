#!/usr/bin/env Rscript

Sys.setenv(
  KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "si_figure7_eigenmodes.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))

data_Supp_Figure7_3 <- function(n_core = 8L, rebuild = FALSE) {
  workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  # Supplementary Figure 7-3 uses only the q10 rank-1 surface and its local
  # grid behavior. Reuse the same-run q20 caches created by data_Figure7.R;
  # do not rebuild the archived ranks 2--10 eigenmode atlas, which is not part
  # of the current figure.
  paths <- f7r_paths(workspace_root)
  f7r_load_response_engine(paths)
  f7x_write_model_provenance(paths)
  objective_bundle <- f7r_objective_selection(paths)
  multiseed_cache <- f7x_compute_joint_invitro_cache(
    paths, objective_bundle,
    n_core = as.integer(n_core), rebuild = isTRUE(rebuild)
  )
  f7x_summarize_joint_invitro(
    paths, objective_bundle, multiseed_cache
  )
  f7x_supplement_7_3_context_data(workspace_root = workspace_root)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1L]])
  }
  data_Supp_Figure7_3(
    n_core = as.integer(value("n-core", "8")),
    rebuild = tolower(value("rebuild", "false")) %in% c("true", "1", "yes")
  )
}
