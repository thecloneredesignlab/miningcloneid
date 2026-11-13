#!/usr/bin/env Rscript
# Independent saved-object verification. No model, compiler or optimizer loads.
o2pilot_validate_saved_objects <- function(out, fixed, manifest, tolerance = 1e-8, fit_name = "accepted_fit.rds") {
  fit <- readRDS(file.path(out, fit_name))
  map <- readRDS(file.path(out, "parameter_map.rds"))
  o2prof_to_free(map, fit$best_params)
  residual <- setNames(numeric(length(fixed)), names(fixed))
  for (name in names(fixed)) {
    expected <- as.numeric(fixed[[name]])
    actual <- as.numeric(fit$best_params[[name]])
    if (!identical(actual, expected)) stop("Saved full-precision target is not exact")
    residual[[name]] <- actual - expected
  }
  o2prof_check_samples(fit$best_components, manifest)
  arithmetic <- -mean(fit$best_components$growth_df$loglik) - mean(fit$best_components$ploidy_df$mean_loglik) - mean(fit$best_components$flow_df$mean_loglik)
  if (abs(arithmetic - fit$best_components$objective) > tolerance) stop("Saved score arithmetic mismatch")
  actual <- readRDS(file.path(out, "native_de_initial_population.rds"))
  requested <- readRDS(file.path(out, "requested_initial_population.rds"))
  if (!identical(actual, requested)) stop("Native starting population differs")
  for (i in seq_len(nrow(actual))) o2prof_to_free(map, o2prof_from_free(map, setNames(actual[i, ], colnames(actual))))
  list(pass = TRUE, objective = fit$best_components$objective, active_parameters = sum(map$table$status == "active"),
       numeric_params_exact = lapply(fit$best_params[vapply(fit$best_params, function(x) is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))],
                                    function(x) sprintf("%.17g", x)),
       target_residuals = as.list(residual), target_constraint_tolerance = 0,
       saved_score_arithmetic_difference = arithmetic - fit$best_components$objective,
       initial_population_identical = TRUE, initial_population_dimensions = dim(actual))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: check_culture_profile_artifact.R <project> <study> <task-id>")
  project <- normalizePath(args[1], mustWork = TRUE)
  study <- normalizePath(args[2], mustWork = TRUE)
  source(file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_profile.R"))
  config <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
  task <- jsonlite::read_json(file.path(study, "tasks", args[3], "task.json"), simplifyVector = TRUE)
  manifest <- read.delim(file.path(project, config$wp0_relative_path, "manifests/observation_manifest.tsv"), stringsAsFactors = FALSE)
  out <- file.path(study, "fits", args[3])
  result <- o2pilot_validate_saved_objects(out, as.list(task$fixed), manifest, config$score_replay_tolerance)
  jsonlite::write_json(result, file.path(out, "artifact_validation_v05.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  cat("Full-precision RDS constraints, samples, arithmetic and native population verified.\n")
}
