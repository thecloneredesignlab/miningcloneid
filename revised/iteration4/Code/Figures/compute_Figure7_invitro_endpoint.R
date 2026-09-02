#!/usr/bin/env Rscript

# Process-isolated Figure 7 in-vitro endpoint worker.
#
# macOS fork workers can inherit an already initialized OpenMP runtime and
# exhaust its shared-memory bookkeeping during the long endpoint grids.  This
# entry point is designed to be launched as several independent R processes;
# workers receive disjoint endpoint indices and use the same atomic caches as
# data_Figure7.R.  Scheduling therefore changes neither model inputs nor
# numerical calculations.

Sys.setenv(
  KMP_USE_SHM = "0", KMP_INIT_AT_FORK = "FALSE",
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
options(stringsAsFactors = FALSE, warn = 1)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("This worker must be run with Rscript.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
})
source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = character()) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
as_flag <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")

stage <- option_value("stage", "q20")
worker_id <- as.integer(option_value("worker-id", "1"))
worker_count <- as.integer(option_value("worker-count", "1"))
rebuild <- as_flag(option_value("rebuild", "FALSE"))
if (!stage %in% c("q20", "dense")) stop("--stage must be q20 or dense.")
if (!is.finite(worker_id) || !is.finite(worker_count) ||
    worker_count < 1L || worker_id < 1L || worker_id > worker_count) {
  stop("Invalid worker-id/worker-count combination.")
}

workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
paths <- f7r_paths(workspace_root)
f7r_load_response_engine(paths)
objective_bundle <- f7r_objective_selection(paths)
model_source_fingerprint <- f7r_model_source_fingerprint(paths)

if (identical(stage, "q20")) {
  manifest <- f7x_joint_context_endpoint_manifest(
    paths, objective_bundle, cutoff = "q20", displayed_only = FALSE,
    output_name = "joint_invitro_q20_endpoint_manifest.tsv"
  )
  endpoints <- manifest$endpoints
  cache_root <- file.path(paths$figure7, "multiseed_endpoint_cache_invitro")
} else {
  manifest <- f7x_joint_context_endpoint_manifest(
    paths, objective_bundle, cutoff = "q10", displayed_only = TRUE,
    output_name = "figure7_invitro_dense_endpoint_manifest.tsv"
  )
  endpoints <- manifest$endpoints
  cache_root <- file.path(paths$figure7, "figure7_invitro_dense_endpoint_cache")
}

contexts <- lapply(
  unique(endpoints$pair_id), f7r_pair_model_context,
  selected = objective_bundle$selected, paths = paths
)
names(contexts) <- unique(endpoints$pair_id)
assigned <- which(((seq_len(nrow(endpoints)) - 1L) %% worker_count) + 1L == worker_id)

message(
  "Figure 7 in-vitro ", stage, " worker ", worker_id, "/", worker_count,
  ": ", length(assigned), " assigned unique endpoints."
)

for (position in seq_along(assigned)) {
  i <- assigned[[position]]
  z <- endpoints[i, , drop = FALSE]
  pair_dir <- file.path(cache_root, z$pair_label[[1L]])
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  cache_path <- file.path(
    pair_dir, paste0("endpoint_", z$parameter_endpoint_group[[1L]], ".rds")
  )
  message(
    "worker ", worker_id, " ", stage, " endpoint ", position, "/",
    length(assigned), ": ", z$parameter_endpoint_group[[1L]]
  )
  qc <- if (identical(stage, "q20")) {
    f7r_compute_seed_cache(
      pair_id = z$pair_id[[1L]],
      seed_number = z$representative_seed_number[[1L]],
      objective = z$representative_objective[[1L]],
      master = objective_bundle$master,
      context = contexts[[z$pair_id[[1L]]]],
      cache_path = cache_path,
      parameter_source = objective_bundle$master_path,
      full_surface = z$endpoint_multiplicity_q10[[1L]] > 0L,
      force_rebuild = rebuild,
      model_context = "in vitro",
      parameter_value_column = "vitro_natural",
      simulation_mode = "invitro",
      model_source_fingerprint = model_source_fingerprint
    )
  } else {
    f7r_figure7d_compute_endpoint_cache(
      metadata = z,
      parameters = objective_bundle$parameters_invitro,
      context = contexts[[z$pair_id[[1L]]]],
      cache_path = cache_path,
      parameter_source = objective_bundle$paths[["parameters_invitro"]],
      force_rebuild = rebuild,
      model_context = "in vitro",
      simulation_mode = "invitro",
      model_source_fingerprint = model_source_fingerprint
    )
  }
  if (!isTRUE(qc$operator_qc_pass[[1L]])) {
    stop("Operator QC failed for ", z$parameter_endpoint_group[[1L]])
  }
}

message("Figure 7 in-vitro ", stage, " worker ", worker_id, " complete.")
