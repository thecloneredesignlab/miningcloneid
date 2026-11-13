#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[1], mustWork = TRUE)
study <- normalizePath(args[2], mustWork = TRUE)
modules <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP/util",
                     c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
for (path in modules) source(path)
config <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
ctx <- o2prof_load_culture(file.path(project, config$wp0_relative_path), modules)
checks <- list()
for (case in c("unrestricted", "no_induction", "flat_survival")) {
  fixed <- switch(case, unrestricted = list(), no_induction = list(p_misseg = 0), flat_survival = list(buffer_beta = 0))
  map <- o2prof_map(ctx$cfg$parameter_table, ctx$saved$best_params, "zero_extended", fixed = fixed)
  params <- o2prof_project_start(map, ctx$saved$best_params)$params
  directory <- file.path(study, "tests/worker_smoke", case)
  if (dir.exists(directory)) stop("Worker smoke evidence already exists")
  tracked <- o2pilot_track_context(ctx, file.path(directory, "best_candidates"), case)
  o2prof_install_culture_map(tracked, map, list(actual_start = params), directory, count_evaluations = TRUE)
  master <- tracked$backend$ivt_objective_components(params, ctx$objects, ctx$cfg)$objective
  cl <- parallel::makePSOCKcluster(2L)
  tryCatch({
    tracked$backend$initialize_invitro_deoptim_workers(cl, tracked$backend$o2simps_cpp_dll_info())
    values <- numeric()
    for (repeat_index in 1:2) values <- c(values, unlist(parallel::clusterCall(cl, function(evaluate, params, objects, cfg) {
      evaluate(params, objects, cfg)$objective
    }, tracked$backend$ivt_objective_components, params, ctx$objects, ctx$cfg)))
    stopifnot(length(values) == 4L, all(abs(values - master) <= 1e-8))
  }, finally = parallel::stopCluster(cl))
  candidates <- lapply(list.files(file.path(directory, "best_candidates"), pattern = "\\.rds$", full.names = TRUE), readRDS)
  stopifnot(length(candidates) == 3L)
  for (candidate in candidates) {
    o2prof_to_free(map, candidate$params)
    stopifnot(identical(candidate$fingerprint, case), abs(candidate$objective - master) <= 1e-8)
    for (target in names(fixed)) stopifnot(identical(candidate$params[[target]], 0))
  }
  checks[[case]] <- list(objective = master, active_parameters = sum(map$table$status == "active"),
                         saved_process_candidates = length(candidates), repeated_worker_calls = length(values), pass = TRUE)
}
o2pilot_json(list(pass = TRUE, cases = checks, optimization = FALSE,
                  module_hashes = as.list(setNames(vapply(modules, o2prof_hash, character(1)), basename(modules)))),
              file.path(study, "tests/worker_smoke_status.json"))
cat("Actual-model worker and saved-candidate smoke checks passed; no optimization was run.\n")
