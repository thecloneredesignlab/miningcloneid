#!/usr/bin/env Rscript
# Versioned, manifest-driven optimization review; no tumor code is loaded.
main <- function() {
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: run_culture_review.R project study mode [task]")
project <- normalizePath(args[1], mustWork = TRUE)
study <- normalizePath(args[2], mustWork = TRUE)
mode <- args[3]; tid <- if (length(args) >= 4L) args[4] else "validation"
code <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP")
modules <- file.path(code, "util", c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R",
                                    "o2_supply_demand_map_pilot.R", "o2_supply_demand_map_review.R"))
for (p in modules) source(p)
o2pilot_verify_code(project, study)
settings <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
if (!isTRUE(settings$authorized) || settings$output_language != "English" ||
    settings$qos != "xxlarge" || settings$core_hour_cap != 10000) stop("Wrong authorization/configuration")
if (Sys.getenv("O2PROF_IMAGE_SHA256") != settings$image_sha256) stop("Image fingerprint mismatch")
ctx <- o2prof_load_culture(file.path(project, settings$wp0_relative_path), modules)
if (mode == "validate") {
  source(file.path(code, "analysis/structural_identifiability/test_culture_review.R"), local = environment())
  replays <- list()
  for (id in c("fit_090", "fit_107", "fit_111", "fit_036")) {
    file <- file.path(project, settings$parent_relative_path, "fits", id, "accepted_fit.rds")
    original <- readRDS(file); fixed <- if (id == "fit_107") list(p_misseg = 0) else
      if (id == "fit_111") list(buffer_beta = 0) else if (id == "fit_036") list(p_misseg = original$best_params$p_misseg) else list()
    map <- o2review_map(ctx, "original", fixed)
    p <- o2prof_project_start(map, original$best_params)$params
    values <- replicate(3L, {
      comp <- ctx$raw_objective(p, ctx$objects, ctx$cfg); o2prof_check_samples(comp, ctx$manifest)
      independent <- o2prof_independent_score(ctx, comp)
      if (abs(comp$objective - independent$objective) > 1e-8) stop("Independent replay arithmetic failed")
      comp$objective
    })
    if (max(abs(values - original$best_components$objective)) > 1e-8 || diff(range(values)) > 1e-8)
      stop("Production replay/repeatability failed")
    replays[[id]] <- list(objective = values[1], repeats = 3L, spread = diff(range(values)), fixed = fixed)
  }
  # Domain mappings alter only two prespecified search ceilings, never cfg.
  for (domain in c("D_base", "D_sat")) {
    map <- o2review_map(ctx, domain)
    original <- readRDS(file.path(project, settings$parent_relative_path, "fits/fit_090/accepted_fit.rds"))
    o2prof_to_free(map, original$best_params)
  }
  provisional <- list(pass = TRUE, production_replays = replays, toy_tests_passed = TRUE,
      code_manifest_sha256 = o2prof_hash(file.path(study, "manifests/code_manifest.tsv")),
      image_sha256 = settings$image_sha256, integration_pass = FALSE,
      completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  o2pilot_json(provisional, file.path(study, "tests/validation_status.json"))
  fixture_ids <- c("unrestricted", "induction_zero", "survival_zero")
  for (i in seq_along(fixture_ids)) {
    id <- fixture_ids[i]; ref <- c("fit_090", "fit_107", "fit_111")[i]
    source_fit <- file.path(settings$parent_relative_path, "fits", ref, "accepted_fit.rds")
    fit <- readRDS(file.path(project, source_fit))
    fixed <- if (i == 2L) list(p_misseg = 0) else if (i == 3L) list(buffer_beta = 0) else list()
    m <- o2review_map(ctx, "original", fixed); anchor <- o2prof_project_start(m, fit$best_params)$params
    folder <- file.path(study, "tests/fixtures", id); dir.create(folder, recursive = TRUE)
    o2pilot_json(list(task_id = id, domain = "original", fixed = fixed, mode = "global", point_id = id,
        min_generations = 1L, max_generations = 2L, seed = 401L + i, pure_independent = FALSE,
        start_class = "integration_fixture", source_fit = source_fit), file.path(folder, "task.json"))
    pop_folder <- file.path(study, "populations", id); dir.create(pop_folder, recursive = TRUE)
    pop <- o2prof_initial_population(m, 256L, list(anchor = anchor))$population
    saveRDS(pop, file.path(pop_folder, "population.rds")); saveRDS(anchor, file.path(pop_folder, "anchor.rds"))
    o2pilot_json(list(pass = TRUE, population_sha256 = o2prof_hash(file.path(pop_folder, "population.rds"))), file.path(pop_folder, "screening_result.json"))
    runner <- file.path(code, "analysis/structural_identifiability/run_culture_review.R")
    for (child_mode in c("fixture", "fixture_check")) {
      status <- system2("bash", vapply(c(file.path(project, "scripts/agentRrunner.sh"), runner,
          project, study, child_mode, id), shQuote, character(1)))
      if (status != 0L) stop("Fresh-process actual-worker integration fixture failed: ", id)
    }
  }
  provisional$integration_pass <- TRUE
  provisional$actual_model_worker_fixtures <- fixture_ids
  provisional$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  o2pilot_json(provisional, file.path(study, "tests/validation_status.json"))
  cat("Review tests, repeatable replays and actual-worker fixtures passed. Fixtures are not profile evidence.\n")
  quit(status = 0L)
}
fixture <- mode %in% c("fixture", "fixture_check")
task_path <- if (fixture) file.path(study, "tests/fixtures", tid, "task.json") else file.path(study, "tasks", tid, "task.json")
task <- jsonlite::read_json(task_path, simplifyVector = TRUE)
if (task$task_id != tid) stop("Wrong task identifier")
fixed <- if (length(task$fixed)) as.list(task$fixed) else list()
map <- o2review_map(ctx, task$domain, fixed)
out <- if (fixture) file.path(study, "tests/fixtures", tid, "fit") else file.path(study, "fits", tid)
if (mode %in% c("check", "fixture_check")) {
  candidate <- jsonlite::read_json(file.path(out, "candidate_result.json"), simplifyVector = TRUE)
  inputs <- jsonlite::read_json(file.path(out, "run_inputs.json"), simplifyVector = TRUE)
  if (inputs$task_sha256 != o2prof_hash(task_path) || inputs$study_sha256 != o2prof_hash(file.path(study, "study.json")) ||
      inputs$implementation_sha256 != o2prof_hash(file.path(study, "manifests/code_manifest.tsv"))) stop("Run input fingerprints changed")
  for (name in names(candidate$files_sha256)) if (o2prof_hash(file.path(out, name)) != candidate$files_sha256[[name]]) stop("Candidate payload changed")
  fit <- readRDS(file.path(out, "accepted_fit.rds"))
  o2prof_to_free(map, fit$best_params)
  for (n in names(fixed)) if (!identical(as.numeric(fit$best_params[[n]]), as.numeric(fixed[[n]]))) stop("Target is not exact")
  comp <- ctx$raw_objective(fit$best_params, ctx$objects, ctx$cfg)
  independent <- o2prof_independent_score(ctx, comp)
  o2prof_check_samples(comp, ctx$manifest)
  if (abs(comp$objective - candidate$objective) > 1e-8 || abs(independent$objective - candidate$objective) > 1e-8) stop("Fresh-process rescoring failed")
  if (task$mode == "global") {
    actual <- readRDS(file.path(out, "native_de_initial_population.rds"))
    provided <- readRDS(file.path(study, "populations", tid, "population.rds"))
    screen <- jsonlite::read_json(file.path(study, "populations", tid, "screening_result.json"), simplifyVector = TRUE)
    if (screen$population_sha256 != o2prof_hash(file.path(study, "populations", tid, "population.rds"))) stop("Screened population fingerprint changed")
    if (!identical(actual, provided)) stop("Actual population was not injected")
    if (candidate$actual_generations < task$min_generations) stop("DE stopped before its explicit minimum")
    if (readRDS(file.path(out, "local_intercept_calls.rds")) != 1L) stop("Unexpected backend local entry count")
  }
  o2pilot_json(list(pass = TRUE, objective = comp$objective, independent_objective = independent$objective,
      candidate_result_sha256 = o2prof_hash(file.path(out, "candidate_result.json")),
      fixed_target_residual_tolerance = 0, fresh_process = TRUE,
      checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(out, "artifact_validation.json"))
  quit(status = 0L)
}
validation <- jsonlite::read_json(file.path(study, "tests/validation_status.json"), simplifyVector = TRUE)
if (!isTRUE(validation$pass) || validation$code_manifest_sha256 != o2prof_hash(file.path(study, "manifests/code_manifest.tsv"))) stop("Unvalidated code snapshot")
if (!fixture && !isTRUE(validation$integration_pass)) stop("Actual-worker integration gate not passed")
if (mode == "screen") {
  folder <- file.path(study, "populations", tid)
  if (dir.exists(folder)) stop("Existing screening attempt cannot be overwritten")
  dir.create(folder, recursive = TRUE)
  set.seed(task$seed)
  started <- as.numeric(Sys.time()); accepted <- list(); records <- list()
  for (i in seq_len(5000L)) {
    if (as.numeric(Sys.time()) - started >= 3500 || length(accepted) >= 64L) break
    proposal <- if (i %% 2L) "full_domain" else "growth_favoring"
    p <- o2review_draw(map, proposal, i); reason <- ""; q <- NA_real_
    q <- tryCatch({
      comp <- ctx$raw_objective(p, ctx$objects, ctx$cfg); o2prof_check_samples(comp, ctx$manifest); comp$objective
    }, error = function(e) { reason <<- conditionMessage(e); NA_real_ })
    records[[i]] <- list(draw = i, proposal = proposal, params = p, objective = q, error = reason)
    if (is.finite(q)) accepted[[length(accepted) + 1L]] <- p
    if (i %% 25L == 0L) o2pilot_atomic_rds(records, file.path(folder, "screening_draws.rds"))
  }
  o2pilot_atomic_rds(records, file.path(folder, "screening_draws.rds"))
  if (length(accepted) < 8L) {
    o2pilot_json(list(pass = FALSE, feasible_members = length(accepted), reason = "Independent feasible-bank threshold not met"), file.path(folder, "screening_result.json"))
    stop("Independent feasible-bank screening unresolved")
  }
  anchor <- if (isTRUE(task$pure_independent)) accepted[[1]] else {
    fit <- readRDS(file.path(project, task$source_fit)); o2prof_project_start(map, fit$best_params)$params
  }
  starts <- list(anchor = anchor)
  population <- o2prof_initial_population(map, 256L, starts)$population
  for (i in seq_along(accepted)) population[i, ] <- o2prof_to_free(map, accepted[[i]])
  if (!isTRUE(task$pure_independent)) population[1, ] <- o2prof_to_free(map, anchor)
  colnames(population) <- o2prof_spec(map, anchor)$param_name
  saveRDS(population, file.path(folder, "population.rds")); saveRDS(anchor, file.path(folder, "anchor.rds"))
  o2pilot_json(list(pass = TRUE, feasible_members = length(accepted), screening_evaluations = length(records),
      seed = task$seed, pure_independent = task$pure_independent, source_fit = task$source_fit,
      proposal_selection = "Protocol feasibility only; no objective ranking", population_sha256 = o2prof_hash(file.path(folder, "population.rds")),
      columns = colnames(population), elapsed_seconds = as.numeric(Sys.time()) - started), file.path(folder, "screening_result.json"))
  quit(status = 0L)
}
if (dir.exists(out)) stop("Existing optimization attempt cannot be overwritten")
dir.create(out, recursive = TRUE)
o2pilot_json(list(task_sha256 = o2prof_hash(task_path), study_sha256 = o2prof_hash(file.path(study, "study.json")),
    implementation_sha256 = o2prof_hash(file.path(study, "manifests/code_manifest.tsv")), source_fit = task$source_fit,
    source_fit_sha256 = if (!is.null(task$source_fit)) o2prof_hash(file.path(project, task$source_fit)) else NULL,
    parent_state_sha256 = settings$parent_state_sha256), file.path(out, "run_inputs.json"))
fingerprint <- o2prof_hash(file.path(out, "run_inputs.json")); saveRDS(map, file.path(out, "parameter_map.rds"))
if (mode %in% c("global", "fixture")) {
  population <- readRDS(file.path(study, "populations", tid, "population.rds"))
  anchor <- readRDS(file.path(study, "populations", tid, "anchor.rds"))
  screen <- jsonlite::read_json(file.path(study, "populations", tid, "screening_result.json"), simplifyVector = TRUE)
  if (!isTRUE(screen$pass) || o2prof_hash(file.path(study, "populations", tid, "population.rds")) != screen$population_sha256) stop("Screened population changed")
  ctx <- o2pilot_track_context(ctx, file.path(out, "best_candidates"), fingerprint)
  o2prof_install_culture_map(ctx, map, list(anchor = anchor), out, population, TRUE)
  intercepted <- new.env(parent = emptyenv()); intercepted$calls <- 0L
  ctx$backend$optim <- o2review_local_stub(intercepted)
  cleanup <- o2prof_trace_native_population(out); on.exit(cleanup())
  control <- o2review_de_controls(task)
  o2pilot_json(c(control, list(requested_steptol = 150L, minimum_generations = task$min_generations,
      enforcement = "Native patience at least minimum generations; actual iteration count checked")), file.path(out, "execution_controls.json"))
  ctx$backend$main(list(parameter_table = ctx$cfg$parameter_table,
      fit_objects_dir = file.path(ctx$wp0_root, "source/inputs/oxygen/ploidyOxygen/data/fit_objects"),
      flow_density_path = file.path(ctx$wp0_root, "source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv"),
      out_dir = out, seed = task$seed + 17L, NP = 256L, n_cores = 8L,
      itermax = control$itermax, itermax_max = control$itermax,
      de_reltol = control$reltol, de_steptol = control$steptol, auto_viz = FALSE))
  saveRDS(intercepted$calls, file.path(out, "local_intercept_calls.rds"))
  native <- readRDS(file.path(out, "fit_result.rds"))
  if (native$optimizer_trace$iter_completed < task$min_generations || intercepted$calls != 1L) stop("Global-only minimum/entry contract failed")
  o2review_finalize(ctx, map, out, fingerprint, task, native = native)
} else if (mode == "local") {
  source_fit <- readRDS(file.path(project, task$source_fit))
  p <- o2prof_project_start(map, source_fit$best_params)$params
  spec <- o2review_unit_spec(map, p)
  evaluation <- o2review_evaluator(ctx, map, out, fingerprint, "local", as.numeric(Sys.time()) + 21000)
  answer <- o2review_pattern(evaluation$fn, spec$unit, o2review_directions(map),
      max_evaluations = 19500L, deadline = as.numeric(Sys.time()) + 20900)
  # Endpoint boundary probes are explicit and included in the evaluation cap.
  boundary <- list(); improved <- FALSE
  for (i in seq_along(answer$par)) if (answer$par[i] < 1e-8 || answer$par[i] > 1 - 1e-8) for (h in c(.001, .01, .05)) {
    if (evaluation$counts$n >= 20000L) { answer$budget_limited <- TRUE; break }
    z <- answer$par; z[i] <- if (z[i] < .5) h else 1 - h
    q <- evaluation$fn(z)
    boundary[[length(boundary) + 1L]] <- data.frame(coordinate = spec$param_name[i], inward_fraction = h,
        endpoint_objective = answer$value, probe_objective = q)
    if (q < answer$value - 1e-5) improved <- TRUE
  }
  if (length(boundary)) o2prof_write_tsv(do.call(rbind, boundary), file.path(out, "boundary_probes.tsv"))
  if (improved) answer$tested_directional_stability <- FALSE
  answer$evaluations <- evaluation$counts$n
  o2prof_write_tsv(answer$sweeps, file.path(out, "local_sweeps.tsv"))
  saveRDS(answer, file.path(out, "local_search.rds"))
  o2review_finalize(ctx, map, out, fingerprint, task, local = answer)
} else stop("Unknown review mode")
cat("Complete candidate checkpointed for external fresh-process validation:", tid, "\n")
}
main()
