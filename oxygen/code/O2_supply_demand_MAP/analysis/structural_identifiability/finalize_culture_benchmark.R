#!/usr/bin/env Rscript
# No simulator: independently confirm external completion, hashes, constraints,
# samples, and saved score replay before accepting the measured endpoint.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: finalize_culture_benchmark.R <project-root> <study-root> <job-id>")
project <- normalizePath(args[1], mustWork = TRUE)
study <- normalizePath(args[2], mustWork = TRUE)
job <- args[3]
source(file.path(project, "oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_profile.R"))
out <- file.path(study, "benchmark/attempt_01")
if (file.exists(file.path(out, "validated_result.json"))) stop("Validated result already exists")
external <- read.delim(file.path(study, "manifests/benchmark_scheduler.tsv"), sep = "|", stringsAsFactors = FALSE, check.names = FALSE)
row <- external[as.character(external$JobID) == job, , drop = FALSE]
if (nrow(row) != 1L || row$State != "COMPLETED" || row$ExitCode != "0:0") stop("External launch did not complete successfully")
candidate <- jsonlite::read_json(file.path(out, "candidate_result.json"), simplifyVector = TRUE)
run <- readRDS(file.path(out, "run_manifest.rds"))
if (!identical(digest::digest(run$inputs, algo = "sha256"), run$fingerprint) || !identical(candidate$slurm_job_id, job)) stop("Wrong run manifest or job")
config <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
wp0 <- file.path(project, config$wp0_relative_path)
contract_path <- file.path(wp0, "manifests/production_contract.json")
if (o2prof_hash(contract_path) != run$inputs$wp0_contract_sha256) stop("Production contract changed")
contract <- jsonlite::read_json(contract_path, simplifyVector = TRUE)
input_manifest <- file.path(wp0, "manifests/input_manifest.tsv")
if (o2prof_hash(input_manifest) != contract$fingerprint_fields$input_manifest_sha256) stop("Frozen input manifest changed")
frozen <- read.delim(input_manifest, stringsAsFactors = FALSE)
if (any(vapply(file.path(wp0, frozen$path), o2prof_hash, character(1)) != frozen$sha256)) stop("Frozen production source/input changed")
modules <- file.path(project, "oxygen/code/O2_supply_demand_MAP/util", names(run$inputs$modules))
if (!identical(setNames(vapply(modules, o2prof_hash, character(1)), names(run$inputs$modules)), run$inputs$modules) ||
    o2prof_hash(file.path(study, "study.json")) != run$inputs$config_sha256 ||
    o2prof_hash(file.path(study, "manifests/benchmark_starts.rds")) != run$inputs$starts_sha256 ||
    o2prof_hash(file.path(project, "oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/run_culture_constrained_fit.R")) != run$inputs$entry_sha256) stop("Run inputs changed")
if (o2prof_hash(file.path(out, "fit_result.rds")) != candidate$fit_result_sha256 ||
    o2prof_hash(file.path(out, "native_de_initial_population.rds")) != candidate$initial_population_file_sha256) stop("Result files changed")
fit <- readRDS(file.path(out, "fit_result.rds"))
map <- readRDS(file.path(out, "parameter_map.rds"))
if (!identical(map$table, run$inputs$parameter_map)) stop("Saved parameter map changed")
actual <- readRDS(file.path(out, "native_de_initial_population.rds"))
if (!identical(actual, readRDS(file.path(out, "requested_initial_population.rds")))) stop("Actual population mismatch")
for (i in seq_len(nrow(actual))) o2prof_to_free(map, o2prof_from_free(map, setNames(actual[i, ], colnames(actual))))
manifest <- read.delim(file.path(project, config$wp0_relative_path, "manifests/observation_manifest.tsv"), stringsAsFactors = FALSE)
verified <- readRDS(file.path(out, "verified_components.rds"))
o2prof_check_samples(verified, manifest)
if (abs(candidate$independent_objective - verified$objective) > config$score_replay_tolerance ||
    abs(fit$best_components$objective - verified$objective) > config$score_replay_tolerance) stop("Saved score mismatch")
candidate$params <- fit$best_params
candidate$launch_completed <- TRUE
candidate$launch_exit_code <- 0L
candidate$launch_status <- "scheduler_completed_exit_0"
acceptance <- o2prof_accept_candidate(candidate, run$fingerprint, map, digest::digest(actual, algo = "sha256"), config$score_replay_tolerance)
candidate$accepted <- acceptance$accepted
candidate$search_status <- acceptance$search_status
candidate$optimization_flags <- character()
if (isTRUE(candidate$budget_limited)) candidate$optimization_flags <- c(candidate$optimization_flags, "declared_search_budget_reached")
if (isTRUE(candidate$optimizer_trace$local_attempted) && !isTRUE(candidate$optimizer_trace$local_convergence == 0L)) {
  candidate$optimization_flags <- c(candidate$optimization_flags, "local_optimizer_nonconverged")
  candidate$search_status <- if (isTRUE(candidate$budget_limited)) "budget_limited_and_optimizer_unresolved_provisional" else "optimizer_unresolved_provisional"
}
candidate$optimizer_convergence_certified <- FALSE
candidate$frozen_input_files_verified_unchanged <- nrow(frozen)
candidate$validated_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
jsonlite::write_json(candidate, file.path(out, "validated_result.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("Endpoint validated after external process completion.\n")
print(candidate[c("objective", "search_status", "fit_wall_seconds", "total_evaluator_invocations")])
