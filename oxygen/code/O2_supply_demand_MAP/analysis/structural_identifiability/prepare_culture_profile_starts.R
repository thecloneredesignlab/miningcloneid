#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_culture_profile_starts.R <project> <study>")
project_root <- normalizePath(args[1], mustWork = TRUE)
pilot_study <- normalizePath(args[2], mustWork = TRUE)
settings <- jsonlite::read_json(file.path(pilot_study, "study.json"), simplifyVector = TRUE)
pilot_modules <- file.path(pilot_study, "implementation/oxygen/code/O2_supply_demand_MAP/util",
                           c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
for (path in pilot_modules) source(path)
o2pilot_verify_code(project_root, pilot_study)
# One fresh process avoids stale pointers from reloaded test-only C++ contexts.
ctx <- o2prof_load_culture(file.path(project_root, settings$wp0_relative_path), pilot_modules)
map <- o2prof_map(ctx$cfg$parameter_table, ctx$saved$best_params, "zero_extended")
dir.create(file.path(pilot_study, "start_bank"), showWarnings = FALSE)
references <- list()
for (seed in c(144L, 371L, 10L)) {
  label <- paste0("seed", seed)
  file <- if (seed == 144L) file.path(ctx$wp0_root, "source/seed144/fit_result.rds") else file.path(pilot_study, "source", label, "fit_result.rds")
  fit <- readRDS(file)
  full <- o2prof_project_start(map, fit$best_params)
  comp <- ctx$raw_objective(full$params, ctx$objects, ctx$cfg)
  o2prof_check_samples(comp, ctx$manifest)
  if (abs(comp$objective - fit$best_components$objective) > 1e-8) stop("Production start endpoint does not replay: ", label)
  o2pilot_json(full$params, file.path(pilot_study, "start_bank", paste0(label, ".json")))
  o2prof_write_tsv(full$projection, file.path(pilot_study, "manifests", paste0(label, "_start_projection.tsv")))
  references[[label]] <- list(source_sha256 = o2prof_hash(file), archived_objective = fit$best_components$objective,
                              replayed_objective = comp$objective)
}
engineering_file <- file.path(project_root, "oxygen/results/structural_identifiability/20261109_culture_v2_wp1/benchmark/attempt_01/fit_result.rds")
engineering <- readRDS(engineering_file)
o2pilot_json(o2prof_project_start(map, engineering$best_params)$params, file.path(pilot_study, "start_bank/engineering_null.json"))
references$engineering_null <- list(source_sha256 = o2prof_hash(engineering_file), status = "initializer_only_not_profile_minimum")
spec <- o2prof_spec(map, o2prof_project_start(map, ctx$saved$best_params)$params)
draws <- list()
proposal_ranges <- list(lam_max = c(1, 1.5), alpha_o2 = c(0.5, 0.75), gamma_growth = c(1.5, 2),
                        mu_hp = c(0.005, 0.01), gamma_mu = c(1.2, 1.6), p_misseg = c(0, 0.003),
                        p_mis_base = c(1e-6, 1e-4), buffer_smax = c(0.8, 1))
o2pilot_json(list(method = "alternating full-domain uniform-transformed and growth-favoring initialization proposals",
                  initialization_only = TRUE, optimizer_population_domain_unchanged = TRUE,
                  proposal_natural_ranges = proposal_ranges,
                  selection = "first protocol-feasible draw; no ranking by score"),
              file.path(pilot_study, "manifests/independent_start_proposal.json"))
for (i in seq_len(8L)) {
  set.seed(settings$master_seed + i * 100003L)
  found <- FALSE
  for (j in seq_len(1000L)) {
    z <- setNames(stats::runif(length(spec$lower), spec$lower, spec$upper), spec$param_name)
    params <- o2prof_from_free(map, z)
    proposal <- "full_domain_uniform_transformed"
    if (j %% 2L == 0L) {
      proposal <- "growth_favoring_random_anchor_only"
      for (name in names(proposal_ranges)) params[[name]] <- stats::runif(1L, proposal_ranges[[name]][1], proposal_ranges[[name]][2])
      o2prof_to_free(map, params)
    }
    reason <- ""
    q <- tryCatch({
      comp <- ctx$raw_objective(params, ctx$objects, ctx$cfg)
      o2prof_check_samples(comp, ctx$manifest)
      comp$objective
    }, error = function(e) { reason <<- conditionMessage(e); NA_real_ })
    draws[[length(draws) + 1L]] <- list(start = i, draw = j, proposal = proposal, objective = q, error = reason, params = params)
    if (j %% 50L == 0L) saveRDS(draws, file.path(pilot_study, "manifests/independent_start_draws.rds"))
    if (is.finite(q)) {
      o2pilot_json(params, file.path(pilot_study, "start_bank", paste0("independent_", i, ".json")))
      references[[paste0("independent_", i)]] <- list(seed = settings$master_seed + i * 100003L, accepted_draw = j, objective = q,
                                                      proposal = proposal, selection = "first feasible mixed-proposal draw; no objective ranking")
      cat("Independent start", i, "accepted at draw", j, "using", proposal, "\n")
      found <- TRUE
      break
    }
  }
  if (!found) {
    saveRDS(draws, file.path(pilot_study, "manifests/independent_start_draws.rds"))
    stop("Independent feasible start screening exhausted")
  }
}
saveRDS(draws, file.path(pilot_study, "manifests/independent_start_draws.rds"))
o2pilot_json(references, file.path(pilot_study, "manifests/start_sources.json"))
cat("Production endpoints replayed and independent start bank prepared.\n")
