#!/usr/bin/env Rscript

# Independent-process endpoint worker for Supplementary Figure 6-4. Using one
# fresh R process per endpoint avoids macOS fork/OpenMP shared-memory failures.

Sys.setenv(
  KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
})

source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))
source(file.path(script_dir, "util", "analysis", "si_figure6_eigenmodes.R"))
source(file.path(script_dir, "util", "analysis", "figure6_context_extension.R"))
source(file.path(
  script_dir, "util", "analysis", "figure6_invitro_extended_o2.R"
))

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) {
    if (!is.null(default)) return(default)
    stop("Missing --", name, " option.")
  }
  sub(paste0("^--", name, "="), "", hit[[1L]])
}

profile <- match.arg(value("profile"), c("standard", "dense"))
index <- as.integer(value("index"))
rebuild <- tolower(value("rebuild", "false")) %in% c("true", "1", "yes")
workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
paths <- s64_paths(workspace_root)
f6r_load_response_engine(paths$base)
bundle <- s64_objective_bundle_from_frozen(paths)

manifest_path <- file.path(
  paths$base$figure6, "figure6_invitro_dense_endpoint_manifest.tsv"
)
endpoints <- f6r_read_tsv(manifest_path)
endpoints <- endpoints[
  endpoints$pair_label %in% s64_pair_labels() &
    endpoints$endpoint_multiplicity_q10 > 0, , drop = FALSE
]
endpoints <- endpoints[order(
  match(endpoints$display_label, s64_cluster_labels()),
  endpoints$representative_objective_rank,
  endpoints$representative_seed_number
), , drop = FALSE]
if (nrow(endpoints) != 115L || index < 1L || index > nrow(endpoints)) {
  stop("Endpoint worker index or frozen q10 endpoint manifest is invalid.")
}
z <- endpoints[index, , drop = FALSE]
selected_context <- f6r_pair_model_context(
  z$pair_id[[1L]], bundle$selected, paths$base
)
base_root <- file.path(
  paths$base$figure6,
  if (profile == "standard") {
    "multiseed_endpoint_cache_invitro"
  } else {
    "figure6_invitro_dense_endpoint_cache"
  }
)
cache_root <- if (profile == "standard") paths$joint_cache else paths$dense_cache
base_cache_path <- file.path(
  base_root, z$pair_label[[1L]],
  paste0("endpoint_", z$parameter_endpoint_group[[1L]], ".rds")
)
cache_path <- file.path(
  cache_root, z$pair_label[[1L]],
  paste0("endpoint_", z$parameter_endpoint_group[[1L]], ".rds")
)
qc <- s64_extended_endpoint_cache(
  metadata = z, parameters = bundle$parameters_invitro,
  context = selected_context, base_cache_path = base_cache_path,
  cache_path = cache_path,
  parameter_source = bundle$paths[["parameters_invitro"]],
  p_profile = profile, model_signature = s64_model_signature(paths),
  rebuild = rebuild
)
if (!isTRUE(qc$operator_qc_pass[[1L]])) stop("Endpoint QC failed.")
message(
  profile, " endpoint ", index, "/115 complete: ",
  z$display_label[[1L]], " ", z$parameter_endpoint_group[[1L]]
)
