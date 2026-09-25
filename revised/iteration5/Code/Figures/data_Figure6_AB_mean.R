#!/usr/bin/env Rscript

# Re-summarize Figure 6A/B from validated endpoint-level caches using
# arithmetic means across the q10 optimizer-seed records. No Figure 6C-F or
# supplementary finite-time calculation is invoked here.

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
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure6_context_extension.R"))

data_Figure6_AB_mean <- function(n_core = 8L) {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  paths <- f6r_paths(workspace_root)
  objective_bundle <- f6r_objective_selection(paths)

  message("Figure 6A: re-summarizing existing in-vivo endpoint caches.")
  invivo_multiseed <- f6r_compute_multiseed(
    paths, objective_bundle, n_core = as.integer(n_core), rebuild = FALSE
  )
  message("Figure 6B: re-summarizing existing dense in-vivo endpoint caches.")
  invivo_dense <- f6r_compute_figure7d_dense(
    paths, objective_bundle, n_core = as.integer(n_core), rebuild = FALSE
  )
  invivo_inverse <- f6r_inverse_panel_data(
    paths, rebuild = FALSE, n_core = as.integer(n_core),
    dense_qc_path = file.path(paths$figure6, "figure7d_dense_endpoint_qc.tsv"),
    output_prefix = "figure7", model_context = "in vivo"
  )

  message("Figure 6A: re-summarizing existing in-vitro endpoint caches.")
  invitro_cache <- f6x_compute_joint_invitro_cache(
    paths, objective_bundle, n_core = as.integer(n_core), rebuild = FALSE
  )
  invitro_multiseed <- f6x_summarize_joint_invitro(
    paths, objective_bundle, invitro_cache
  )
  message("Figure 6B: re-summarizing existing dense in-vitro endpoint caches.")
  invitro_dense <- f6x_compute_dense_invitro(
    paths, objective_bundle, n_core = as.integer(n_core), rebuild = FALSE
  )

  invisible(list(
    paths = paths,
    invivo_multiseed = invivo_multiseed,
    invivo_dense = invivo_dense,
    invivo_inverse = invivo_inverse,
    invitro_multiseed = invitro_multiseed,
    invitro_dense = invitro_dense
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  option_value <- function(name, default) {
    hit <- args[startsWith(args, paste0("--", name, "="))]
    if (!length(hit)) return(default)
    sub(paste0("^--", name, "="), "", hit[[1L]])
  }
  data_Figure6_AB_mean(
    n_core = as.integer(option_value("n-core", "8"))
  )
}
