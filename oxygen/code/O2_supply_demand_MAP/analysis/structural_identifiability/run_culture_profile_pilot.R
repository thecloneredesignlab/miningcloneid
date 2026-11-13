#!/usr/bin/env Rscript
# One manifest-specified production-v2 fit; no scheduling or CI calculations.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: run_culture_profile_pilot.R <project> <study> <task-id>")
  project <- normalizePath(args[1], mustWork = TRUE)
  study <- normalizePath(args[2], mustWork = TRUE)
  task_id <- args[3]
  modules <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP/util",
                       c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
  for (path in modules) source(path, local = .GlobalEnv)
  o2pilot_verify_code(project, study)
  config_path <- file.path(study, "study.json")
  config <- jsonlite::read_json(config_path, simplifyVector = TRUE)
  task_path <- file.path(study, "tasks", task_id, "task.json")
  task <- jsonlite::read_json(task_path, simplifyVector = TRUE)
  if (task$task_id != task_id || !isTRUE(config$batch_profile_authorized) || config$inference_calibration != "absent") stop("Invalid approved task")
  if (Sys.getenv("O2PROF_IMAGE_SHA256") != config$image_sha256) stop("Unverified image")
  validated <- jsonlite::read_json(file.path(study, "tests/pilot_validation_status.json"), simplifyVector = TRUE)
  hashes <- setNames(vapply(modules, o2prof_hash, character(1)), basename(modules))
  if (!isTRUE(validated$pass) || !identical(unlist(validated$module_hashes), hashes) ||
      validated$code_manifest_sha256 != o2prof_hash(file.path(study, "manifests/code_manifest.tsv"))) stop("Current implementation not validated")
  ctx <- o2prof_load_culture(file.path(project, config$wp0_relative_path), modules)
  fixed <- if (length(task$fixed)) as.list(task$fixed) else list()
  map <- o2prof_map(ctx$cfg$parameter_table, ctx$saved$best_params, "zero_extended", fixed = fixed,
                    zero_scales = unlist(config$zero_scales))
  starts_path <- file.path(study, "tasks", task_id, "starts.json")
  starts <- jsonlite::read_json(starts_path, simplifyVector = FALSE)
  if (!length(starts)) stop("Empty start bank")
  out <- file.path(study, "fits", task_id)
  if (dir.exists(out)) stop("Existing attempt must not be overwritten")
  dir.create(out, recursive = TRUE)
  inputs <- list(task_id = task_id, task_sha256 = o2prof_hash(task_path),
                 starts_sha256 = o2prof_hash(starts_path), config_sha256 = o2prof_hash(config_path),
                 code_manifest_sha256 = o2prof_hash(file.path(study, "manifests/code_manifest.tsv")),
                 wp0_contract_sha256 = o2prof_hash(file.path(ctx$wp0_root, "manifests/production_contract.json")),
                 modules = as.list(hashes))
  o2pilot_json(inputs, file.path(out, "run_inputs.json"))
  fingerprint <- o2prof_hash(file.path(out, "run_inputs.json"))
  saveRDS(map, file.path(out, "parameter_map.rds"))
  o2prof_write_tsv(map$table, file.path(out, "parameter_status.tsv"))
  ctx <- o2pilot_track_context(ctx, file.path(out, "best_candidates"), fingerprint)
  o2prof_install_culture_map(ctx, map, starts, out, count_evaluations = TRUE)
  cleanup <- o2prof_trace_native_population(out)
  on.exit(cleanup(), add = TRUE)
  options(o2prof.phase = "preflight")
  started <- Sys.time()
  backend_time <- system.time(ctx$backend$main(list(
    parameter_table = ctx$cfg$parameter_table,
    fit_objects_dir = file.path(ctx$wp0_root, "source/inputs/oxygen/ploidyOxygen/data/fit_objects"),
    flow_density_path = file.path(ctx$wp0_root, "source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv"),
    out_dir = out, seed = task$seed, NP = 256L, n_cores = 8L,
    itermax = task$DE_itermax, itermax_max = task$DE_itermax,
    de_reltol = 1e-4, de_steptol = 25L, local_optim_maxit = task$local_maxit, auto_viz = FALSE)))
  fit_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  native <- readRDS(file.path(out, "fit_result.rds"))
  options(o2prof.phase = "verification")
  files <- list.files(file.path(out, "best_candidates"), pattern = "\\.rds$", full.names = TRUE)
  candidates <- lapply(files, readRDS)
  candidates <- c(list(list(params = native$best_params, objective = native$best_components$objective,
                           fingerprint = fingerprint, phase = "native_endpoint")), candidates)
  candidates <- candidates[order(vapply(candidates, function(x) x$objective, numeric(1)))]
  checks <- list()
  winner <- NULL
  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    if (!identical(candidate$fingerprint, fingerprint)) stop("Stale saved candidate")
    o2prof_to_free(map, candidate$params)
    comp <- ctx$backend$ivt_objective_components(candidate$params, ctx$objects, ctx$cfg)
    independent <- o2prof_independent_score(ctx, comp)
    if (abs(comp$objective - candidate$objective) > config$score_replay_tolerance ||
        abs(independent$objective - comp$objective) > config$score_replay_tolerance) stop("Saved candidate replay mismatch")
    checks[[i]] <- data.frame(candidate = i, phase = candidate$phase, objective = candidate$objective,
                              reevaluated = comp$objective, independent = independent$objective)
    if (is.null(winner) || comp$objective < winner$objective) {
      winner <- candidate; winner$objective <- comp$objective; winner$components <- comp
      winner$independent <- independent$objective
    }
  }
  o2prof_write_tsv(do.call(rbind, checks), file.path(out, "candidate_checks.tsv"))
  requested <- readRDS(file.path(out, "requested_initial_population.rds"))
  actual <- readRDS(file.path(out, "native_de_initial_population.rds"))
  if (!identical(requested, actual)) stop("Actual native population mismatch")
  for (i in seq_len(nrow(actual))) o2prof_from_free(map, setNames(actual[i, ], colnames(actual)))
  p <- map$table
  p$estimate <- vapply(p$parameter, function(n) as.numeric(winner$params[[n]]), numeric(1))
  p$at_lower <- abs(p$estimate - p$lower) <= 1e-8 * pmax(1, abs(p$lower))
  p$at_upper <- abs(p$estimate - p$upper) <= 1e-8 * pmax(1, abs(p$upper))
  o2prof_write_tsv(p, file.path(out, "parameter_result.tsv"))
  saveRDS(list(best_params = winner$params, best_components = winner$components,
               cfg = ctx$cfg, origin_phase = winner$phase), file.path(out, "accepted_fit.rds"))
  records <- do.call(rbind, lapply(list.files(file.path(out, "evaluation_logs"), pattern = "\\.tsv$", full.names = TRUE),
                                  read.delim, stringsAsFactors = FALSE))
  o2prof_write_tsv(records, file.path(out, "evaluations.tsv"))
  counts <- as.data.frame(table(records$phase), stringsAsFactors = FALSE)
  names(counts) <- c("phase", "evaluations")
  o2prof_write_tsv(counts, file.path(out, "evaluation_counts.tsv"))
  errors <- records$error[!is.na(records$error) & nzchar(records$error)]
  if (any(grepl("checkpoint_write_failed|Stale candidate|Immutable culture configuration|Wrong free-vector coordinates", errors))) stop("Bookkeeping/configuration errors; inspect evaluation log")
  component_means <- c(growth = mean(winner$components$growth_df$loglik),
                       karyotype = mean(winner$components$ploidy_df$mean_loglik),
                       flow = mean(winner$components$flow_df$mean_loglik))
  result_files <- c("fit_result.rds", "accepted_fit.rds", "parameter_map.rds", "parameter_status.tsv", "parameter_result.tsv",
                    "requested_initial_population.rds", "native_de_initial_population.rds", "candidate_checks.tsv", "evaluations.tsv", "evaluation_counts.tsv")
  budget_limited <- native$optimizer_trace$iter_completed >= native$optimizer_trace$iter_target ||
    isTRUE(native$optimizer_trace$local_convergence == 1L)
  result <- list(task_id = task_id, fingerprint = fingerprint, objective = winner$objective,
                 native_objective = native$best_components$objective, winner_phase = winner$phase,
                 reevaluated_objective = winner$components$objective, independent_objective = winner$independent,
                 params = winner$params, component_means = as.list(component_means), fixed = fixed,
                 active_parameters = sum(p$status == "active"), samples_valid = TRUE, population_identical = TRUE,
                 sample_counts = list(growth = nrow(winner$components$growth_df), karyotype = nrow(winner$components$ploidy_df), flow = nrow(winner$components$flow_df)),
                 files_sha256 = as.list(setNames(vapply(file.path(out, result_files), o2prof_hash, character(1)), result_files)),
                 fit_wall_seconds = fit_seconds, R_system_time = as.list(backend_time),
                 total_evaluator_invocations = nrow(records), failed_evaluator_invocations = length(errors),
                 native_de_nfeval = native$deoptim$optim$nfeval, optimizer_trace = native$optimizer_trace,
                 selected_passage_days = lapply(c("2N", "4N"), function(cohort) {
                   run <- winner$components[[paste0("run_", cohort)]]
                   setNames(lapply(run$segment_results, function(segment) segment$selection$selected_day), names(run$segment_results))
                 }),
                 budget_limited = budget_limited, optimizer_convergence_certified = FALSE,
                 launch_completed = FALSE, slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
                 finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  o2pilot_json(result, file.path(out, "candidate_result.json"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("Fit candidate ready for external scheduler validation:", task_id, "Q=", winner$objective, "\n")
}
main()
