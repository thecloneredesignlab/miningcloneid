#!/usr/bin/env Rscript
# One explicitly configured constrained fit for engineering/cost measurement.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: run_culture_constrained_fit.R <project-root> <engineering-study-root>")
  project <- normalizePath(args[1], mustWork = TRUE)
  study <- normalizePath(args[2], mustWork = TRUE)
  config <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
  if (!identical(config$batch_profile_authorized, FALSE) || config$domain != "zero_extended" || config$benchmark$DE_NP != 256L) stop("Unsupported study or production population size")
  if (Sys.getenv("O2PROF_IMAGE_SHA256") != config$image_sha256) stop("Unverified runtime image")
  modules <- file.path(project, "oxygen/code/O2_supply_demand_MAP/util", c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R"))
  for (path in modules) source(path, local = .GlobalEnv)
  validation <- jsonlite::read_json(file.path(study, "tests/validation_status.json"), simplifyVector = TRUE)
  hashes <- setNames(vapply(modules, o2prof_hash, character(1)), basename(modules))
  if (!isTRUE(validation$pass) || !identical(unlist(validation$module_hashes), hashes)) stop("Necessary tests did not pass for the current modules")
  ctx <- o2prof_load_culture(file.path(project, config$wp0_relative_path), modules)
  map <- o2prof_map(ctx$cfg$parameter_table, ctx$saved$best_params, config$domain, fixed = as.list(config$fixed), zero_scales = unlist(config$zero_scales))
  starts_path <- file.path(study, "manifests/benchmark_starts.rds")
  starts <- readRDS(starts_path)
  fingerprint_data <- list(config = config, modules = hashes, parameter_map = map$table,
                           config_sha256 = o2prof_hash(file.path(study, "study.json")),
                           starts_sha256 = o2prof_hash(starts_path),
                           wp0_contract_sha256 = o2prof_hash(file.path(ctx$wp0_root, "manifests/production_contract.json")),
                           entry_sha256 = o2prof_hash(file.path(project, "oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/run_culture_constrained_fit.R")))
  fingerprint <- digest::digest(fingerprint_data, algo = "sha256")
  out <- file.path(study, "benchmark", "attempt_01")
  if (dir.exists(out)) stop("Existing attempt must not be overwritten; allocate and register a new attempt")
  dir.create(out, recursive = TRUE)
  saveRDS(list(fingerprint = fingerprint, inputs = fingerprint_data), file.path(out, "run_manifest.rds"))
  json_inputs <- fingerprint_data
  json_inputs$modules <- as.list(hashes)
  jsonlite::write_json(list(fingerprint = fingerprint, inputs = json_inputs,
                           stage = "bounded_cost_benchmark", inference_calibration = "absent",
                           production_source_identity_caveat = ctx$contract$runtime$historical_identity_caveat),
                       file.path(out, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  o2prof_write_tsv(map$table, file.path(out, "parameter_status.tsv"))
  saveRDS(map, file.path(out, "parameter_map.rds"))
  o2prof_install_culture_map(ctx, map, starts, out, count_evaluations = TRUE)
  cleanup <- o2prof_trace_native_population(out)
  on.exit(cleanup(), add = TRUE)
  options(o2prof.phase = "preflight")
  fit_started <- Sys.time()
  backend_time <- system.time(ctx$backend$main(list(
    parameter_table = ctx$cfg$parameter_table,
    fit_objects_dir = file.path(ctx$wp0_root, "source/inputs/oxygen/ploidyOxygen/data/fit_objects"),
    flow_density_path = file.path(ctx$wp0_root, "source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv"),
    out_dir = out, seed = config$benchmark$seed,
    itermax = config$benchmark$DE_itermax, itermax_max = config$benchmark$DE_itermax,
    NP = config$benchmark$DE_NP, n_cores = config$benchmark$n_cores,
    de_reltol = config$benchmark$de_reltol, de_steptol = config$benchmark$de_steptol,
    local_optim_maxit = config$benchmark$local_maxit, auto_viz = FALSE)))
  fit_seconds <- as.numeric(difftime(Sys.time(), fit_started, units = "secs"))
  fit <- readRDS(file.path(out, "fit_result.rds"))
  options(o2prof.phase = "verification")
  checked <- ctx$backend$ivt_objective_components(fit$best_params, ctx$objects, ctx$cfg)
  o2prof_check_samples(checked, ctx$manifest)
  independent <- o2prof_independent_score(ctx, checked)
  if (abs(checked$objective - fit$best_components$objective) > config$score_replay_tolerance ||
      abs(independent$objective - checked$objective) > config$score_replay_tolerance) stop("Independent score replay failed")
  o2prof_to_free(map, fit$best_params)
  requested <- readRDS(file.path(out, "requested_initial_population.rds"))
  actual <- readRDS(file.path(out, "native_de_initial_population.rds"))
  if (!identical(requested, actual)) stop("Native optimizer starting population mismatch")
  records <- lapply(list.files(file.path(out, "evaluation_logs"), pattern = "\\.tsv$", full.names = TRUE), read.delim, stringsAsFactors = FALSE)
  records <- do.call(rbind, records)
  o2prof_write_tsv(records, file.path(out, "evaluations.tsv"))
  phase_counts <- as.data.frame(table(records$phase), stringsAsFactors = FALSE)
  names(phase_counts) <- c("phase", "evaluations")
  o2prof_write_tsv(phase_counts, file.path(out, "evaluation_counts.tsv"))
  p <- map$table
  values <- vapply(p$parameter, function(n) as.numeric(fit$best_params[[n]]), numeric(1))
  p$estimate <- values
  p$at_lower <- abs(values - p$lower) <= 1e-8 * pmax(1, abs(p$lower))
  p$at_upper <- abs(values - p$upper) <= 1e-8 * pmax(1, abs(p$upper))
  o2prof_write_tsv(p, file.path(out, "parameter_result.tsv"))
  saveRDS(checked, file.path(out, "verified_components.rds"))
  budget_limited <- fit$optimizer_trace$iter_completed >= fit$optimizer_trace$iter_target ||
    isTRUE(fit$optimizer_trace$local_convergence == 1L)
  result <- list(
    fingerprint = fingerprint, numerical_valid = TRUE,
    launch_completed = FALSE, launch_status = "pending_external_scheduler_verification",
    initial_population_hash = digest::digest(actual, algo = "sha256"),
    fit_result_sha256 = o2prof_hash(file.path(out, "fit_result.rds")),
    initial_population_file_sha256 = o2prof_hash(file.path(out, "native_de_initial_population.rds")),
    objective = fit$best_components$objective, reevaluated_objective = checked$objective,
    independent_objective = independent$objective, samples_valid = TRUE,
    params = fit$best_params, fixed = map$fixed, active_parameters = sum(p$status == "active"),
    budget_limited = budget_limited, optimizer_trace = fit$optimizer_trace,
    fit_wall_seconds = fit_seconds, R_system_time = as.list(backend_time),
    total_evaluator_invocations = nrow(records), phase_counts = phase_counts,
    failed_evaluator_invocations = sum(!is.na(records$error) & nzchar(records$error)),
    native_de_nfeval = fit$deoptim$optim$nfeval,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
    finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(result, file.path(out, "candidate_result.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("Bounded constrained-fit benchmark completed; external launch verification remains.\n")
  print(result[c("objective", "fit_wall_seconds", "total_evaluator_invocations", "budget_limited")])
}
main()
