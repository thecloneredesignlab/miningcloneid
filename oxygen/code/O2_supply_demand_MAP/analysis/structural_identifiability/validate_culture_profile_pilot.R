#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: validate_culture_profile_pilot.R <project> <study>")
project_root <- normalizePath(args[1], mustWork = TRUE)
pilot_study <- normalizePath(args[2], mustWork = TRUE)
code_root <- file.path(pilot_study, "implementation")
# Existing production/kernel/mapping integration regression suite, unchanged.
source(file.path(code_root, "oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/test_culture_profile.R"))
pilot_modules <- file.path(code_root, "oxygen/code/O2_supply_demand_MAP/util",
                           c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
for (path in pilot_modules) source(path)
o2pilot_verify_code(project_root, pilot_study)
ctx <- o2prof_load_culture(file.path(project_root, settings$wp0_relative_path), pilot_modules)
extra <- as.data.frame(testthat::test_file(file.path(code_root, "oxygen/tests/testthat/test-culture-pilot.R"),
                                          reporter = "summary", env = .GlobalEnv, load_helpers = FALSE, stop_on_failure = TRUE))
o2prof_write_tsv(extra[, intersect(c("file", "test", "nb", "failed", "passed", "error"), names(extra)), drop = FALSE],
                  file.path(pilot_study, "tests/pilot_additional_results.tsv"))
# Prepare forward-model inputs in a fresh R process after tests that reload DLLs.
status <- system2("bash", c(shQuote(file.path(project_root, "scripts/agentRrunner.sh")),
                            shQuote(file.path(code_root, "oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/prepare_culture_profile_starts.R")),
                            shQuote(project_root), shQuote(pilot_study)))
if (status != 0L) stop("Fresh-process start preparation failed")
references <- jsonlite::read_json(file.path(pilot_study, "manifests/start_sources.json"), simplifyVector = FALSE)
previous <- jsonlite::read_json(file.path(pilot_study, "tests/validation_status.json"), simplifyVector = TRUE)
o2pilot_json(list(pass = TRUE, test_cases = previous$test_cases + nrow(extra),
                  passed_expectations = previous$passed_expectations + sum(extra$passed), failed = 0L,
                  module_hashes = as.list(setNames(vapply(pilot_modules, o2prof_hash, character(1)), basename(pilot_modules))),
                  code_manifest_sha256 = o2prof_hash(file.path(pilot_study, "manifests/code_manifest.tsv")),
                  production_endpoint_replays = references[c("seed144", "seed371", "seed10")],
                  independent_start_count = 8L, completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
              file.path(pilot_study, "tests/pilot_validation_status.json"))
cat("Pilot implementation and start bank validated; no biological optimization was run.\n")
