#!/usr/bin/env Rscript

# Resume only the aggregation/publication stage after a validated targeted
# precision rescue.  This entry point never calls the stochastic propagator.

Sys.setenv(
  KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("This entry point must be run with Rscript.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
})
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure6_context_extension.R"))
source(file.path(script_dir, "util", "analysis", "figure6_finite_time_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_invitro_passage_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_precision_rescue.R"))

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = character()) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
base_run <- option_value("base-run")
run_id <- option_value("run-id")
if (!nzchar(base_run) || !nzchar(run_id)) {
  stop("Required arguments: --base-run=ABSOLUTE_PATH --run-id=ID")
}

workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
paths <- f6r_paths(workspace_root)
base_run <- f6pr_assert_base_root(paths, base_run)
run_paths <- f6g_paths(paths, run_id = run_id, create = FALSE)
run_root <- normalizePath(run_paths$run_root, mustWork = TRUE)
if (identical(run_root, base_run)) stop("Rescue run must differ from base run.")

validation_path <- file.path(
  run_root, "precision_rescue_condition_validation.tsv"
)
f6r_require_files(validation_path, "precision-rescue condition validation")
condition_validation <- f6r_read_tsv(validation_path)
if (nrow(condition_validation) != 1L ||
    !isTRUE(condition_validation$passed[[1L]]) ||
    condition_validation$maximum_mcse_N[[1L]] > 0.01 ||
    condition_validation$n_replicate_per_optimizer_endpoint[[1L]] != 200L ||
    !isTRUE(condition_validation$independent_from_base_production[[1L]])) {
  stop("Only the audited independent R=200 precision result may be aggregated.")
}

task_manifest_path <- file.path(run_root, "full_range_task_manifest.tsv")
f6r_require_files(task_manifest_path, "rescue task manifest")
tasks <- f6r_read_tsv(task_manifest_path)
if (nrow(tasks) != 420L || any(!file.exists(tasks$cache_path))) {
  stop("Rescue task manifest is incomplete.")
}
condition <- f6pr_condition()
selected <- tasks[
  tasks$model_context == condition$model_context &
    tasks$pair_label == condition$pair_label &
    abs(tasks$p_misseg - condition$p_misseg) < 1e-12 &
    tasks$task_id == "G0337", , drop = FALSE
]
if (nrow(selected) != 1L ||
    dirname(selected$cache_path[[1L]]) != file.path(run_root, "task_cache")) {
  stop("The rescue manifest does not reference its patched G0337 cache.")
}
patched <- readRDS(selected$cache_path[[1L]])
if (!isTRUE(patched$qc$passed) || !isTRUE(patched$qc$mcse_target_met) ||
    patched$qc$maximum_mcse_N > 0.01 ||
    !isTRUE(patched$precision_rescue$independent_random_streams)) {
  stop("Patched task cache did not pass its precision-rescue contract.")
}
base_cache <- file.path(
  base_run, "task_cache", basename(selected$cache_path[[1L]])
)
f6r_require_files(base_cache, "immutable base task cache")
base_fingerprint <- readRDS(base_cache)$fingerprint
rescue_fingerprint <- patched$fingerprint
if (identical(base_fingerprint, rescue_fingerprint)) {
  stop("Rescue and base fingerprints must differ.")
}

config_path <- file.path(run_root, "stochastic_config.rds")
rule_path <- file.path(run_root, "canonical_passage_rule.tsv")
f6r_require_files(c(config_path, rule_path), "rescue stochastic metadata")
stochastic_config <- readRDS(config_path)
if (stochastic_config$master_seed != 20260907L ||
    stochastic_config$replicates != 200L ||
    stochastic_config$allocation != "fixed" ||
    stochastic_config$mcse_target_N != 0.01) {
  stop("Unexpected stochastic configuration for aggregation-only rescue.")
}
passage_bundle <- list(
  rule = f6r_read_tsv(rule_path),
  stochastic = list(config = stochastic_config)
)

aggregation <- f6g_aggregate(
  tasks, run_paths, rescue_fingerprint, passage_bundle, smoke = FALSE,
  accepted_cache_fingerprints = c(base_fingerprint, rescue_fingerprint)
)
precision <- do.call(rbind, lapply(
  tolower(f6ft_family_levels()),
  function(family) f6r_read_tsv(file.path(
    run_root, paste0("stochastic_precision_", family, ".tsv")
  ))
))
if (!all(precision$passed) || max(precision$maximum_mcse_N) > 0.01) {
  stop("Global Figure 6 precision gate failed after rescue aggregation.")
}

provenance <- data.frame(
  key = c(
    "profile", "run_id", "resume_stage", "rescue_policy",
    "base_run_root", "base_fingerprint", "rescue_fingerprint",
    "model_code_root", "model_source_fingerprint", "condition",
    "replicates", "master_seed", "rng_independence", "validation"
  ),
  value = c(
    f6g_profile(), run_paths$run_id,
    "aggregation_only_after_numeric_type_normalization",
    "immutable pilot; one failed condition; fixed-R independent production",
    base_run, base_fingerprint, rescue_fingerprint,
    normalizePath(paths$oxygen_code, mustWork = TRUE),
    f6r_model_source_fingerprint(paths),
    "in vitro|C02|initial=3N|O2=0.7|p_misseg=0.01",
    "200", "20260907",
    "new L'Ecuyer-CMRG stream catalog; no pilot production draw reused",
    validation_path
  ), stringsAsFactors = FALSE
)
f6ft_atomic_write_tsv(
  provenance,
  file.path(run_root, "precision_rescue_run_provenance.tsv")
)
f6g_publish_current(paths, run_paths, rescue_fingerprint)
message(
  "Figure 6 aggregation-only rescue published: ", run_paths$run_id,
  "; global maximum MCSE=", signif(max(precision$maximum_mcse_N), 7)
)
invisible(aggregation)
