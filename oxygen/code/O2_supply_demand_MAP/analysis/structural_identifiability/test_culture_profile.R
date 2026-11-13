#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: test_culture_profile.R <project-root> <engineering-study-root>")
project_root <- normalizePath(args[1], mustWork = TRUE)
study <- normalizePath(args[2], mustWork = TRUE)
module_paths <- file.path(project_root, "oxygen/code/O2_supply_demand_MAP/util",
                          c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R"))
for (f in module_paths) source(f)
settings <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
wp0 <- file.path(project_root, settings$wp0_relative_path)
dir.create(file.path(study, "tests"), recursive = TRUE, showWarnings = FALSE)
profile_validation_root <- file.path(study, "tests")
profile_map_module <- module_paths[1]
ctx <- o2prof_load_culture(wp0, module_paths)
profile_test_table <- read.csv(ctx$cfg$parameter_table, stringsAsFactors = FALSE)
profile_test_params <- ctx$saved$best_params
# Reuse the existing independent kernel/accounting oracles without its model
# bootstrap; the production model is already loaded into ctx's environment.
oracle <- parse(file.path(project_root, "oxygen/tests/testthat/helper-oracle.R"))
for (expr in oracle) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% c("oracle_buffering_survival", "triplet_to_sparse", "shift_weight", "oracle_live_mass_per_division")) eval(expr)
}
cpp_names <- c(ls(ctx$backend, pattern = "^cpp_o2simps_"), ".lambda_eff_of_O2")
invisible(list2env(mget(cpp_names, envir = ctx$backend, inherits = TRUE), envir = .GlobalEnv))
test_files <- file.path(project_root, "oxygen/tests/testthat", c(
  "test-profile-objective.R", "test-buffer-missegregation.R", "test-event-accounting.R", "test-profile-culture-v2.R"))
all_results <- list()
for (path in test_files) {
  result <- testthat::test_file(path, reporter = "summary", env = .GlobalEnv, load_helpers = FALSE, stop_on_failure = TRUE)
  all_results[[basename(path)]] <- as.data.frame(result)
}
saveRDS(all_results, file.path(study, "tests/test_results.rds"))
summary <- do.call(rbind, all_results)
o2prof_write_tsv(summary[, intersect(c("file", "context", "test", "nb", "failed", "skipped", "error", "warning", "passed"), names(summary)), drop = FALSE], file.path(study, "tests/test_results.tsv"))
versions <- setNames(vapply(c("Rcpp", "DEoptim", "testthat", "Matrix", "jsonlite"), function(x) as.character(getNamespaceVersion(x)), character(1)), c("Rcpp", "DEoptim", "testthat", "Matrix", "jsonlite"))
jsonlite::write_json(list(pass = TRUE, files = basename(test_files), test_cases = nrow(summary),
                         passed_expectations = sum(summary$passed), failed = sum(summary$failed),
                         module_hashes = as.list(setNames(vapply(module_paths, o2prof_hash, character(1)), basename(module_paths))),
                         package_versions = as.list(versions), production_commit = ctx$contract$production_reference$commit,
                         completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(study, "tests/validation_status.json"), auto_unbox = TRUE, pretty = TRUE)
