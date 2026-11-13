#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[1], mustWork = TRUE)
study <- normalizePath(args[2], mustWork = TRUE)
modules <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP/util",
                     c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
for (path in modules) source(path)
settings <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
ctx <- o2prof_load_culture(file.path(project, settings$wp0_relative_path), modules)
path <- file.path(study, "source/seed325/fit_result.rds")
out <- file.path(study, "start_bank/seed325.json")
if (file.exists(out)) stop("Additional start has already been prepared")
fit <- readRDS(path)
map <- o2prof_map(ctx$cfg$parameter_table, ctx$saved$best_params, "zero_extended")
full <- o2prof_project_start(map, fit$best_params)
comp <- ctx$raw_objective(full$params, ctx$objects, ctx$cfg)
o2prof_check_samples(comp, ctx$manifest)
independent <- o2prof_independent_score(ctx, comp)
if (abs(comp$objective - fit$best_components$objective) > 1e-8 || abs(comp$objective - independent$objective) > 1e-8)
  stop("Additional production endpoint did not replay")
jsonlite::write_json(full$params, out, auto_unbox = TRUE, pretty = TRUE, digits = 17)
o2prof_write_tsv(full$projection, file.path(study, "manifests/seed325_start_projection.tsv"))
o2pilot_json(list(pass = TRUE, seed = 325L, source_sha256 = o2prof_hash(path),
                  archived_objective = fit$best_components$objective, replayed_objective = comp$objective,
                  independent_objective = independent$objective,
                  purpose = "Independent archived low-score starting vector; not a copy of the current best fit and not a proven distinct basin"),
              file.path(study, "manifests/seed325_start_replay.json"))
cat("Additional seed325 starting vector verified against frozen production v2.\n")
