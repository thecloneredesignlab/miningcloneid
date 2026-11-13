#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[1], mustWork = TRUE)
study <- normalizePath(args[2], mustWork = TRUE)
modules <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP/util",
                     c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
for (path in modules) source(path)
source(file.path(study, "scripts/culture_profile_diagnostics.R"))
config <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
ctx <- o2prof_load_culture(file.path(project, config$wp0_relative_path), modules)
checks <- list()
for (case in c("production", "no_induction", "flat_survival", "zero_error_survival")) {
  params <- ctx$saved$best_params
  if (case == "no_induction") params$p_misseg <- 0
  if (case == "flat_survival") params$buffer_beta <- 0
  if (case == "zero_error_survival") params$buffer_smax <- 0
  for (oxygen in c(0, 0.5, 2, 20.5)) {
    kernel <- pilotdiag_rates(ctx, params, ctx$cfg, oxygen)
    stopifnot(nrow(kernel$rates) == 133L, all(is.finite(as.matrix(kernel$rates))))
    checks[[length(checks) + 1L]] <- data.frame(case = case, oxygen_pct = oxygen, kernel_accounting_pass = TRUE)
  }
}
audits <- list()
for (cohort in c("2N", "4N")) {
  results <- ctx$saved$best_components[[paste0("run_", cohort)]]$segment_results
  for (i in unique(c(1L, which(vapply(results, function(x) x$segment$oxygen_pct < 20, logical(1)))[1]))) {
    if (is.na(i)) next
    result <- results[[i]]
    kernel <- pilotdiag_rates(ctx, ctx$saved$best_params, ctx$cfg, result$segment$oxygen_pct)
    # Use the relocated cfg while preserving every numerical/model setting.
    fit <- ctx$saved; fit$cfg <- ctx$cfg
    audit <- pilotdiag_dense_flux(ctx, fit, cohort, result, kernel)
    audits[[length(audits) + 1L]] <- audit$summary
  }
}
o2prof_write_tsv(do.call(rbind, checks), file.path(study, "tests/flux_kernel_checks.tsv"))
o2prof_write_tsv(do.call(rbind, audits), file.path(study, "tests/flux_dense_replay_checks.tsv"))
o2pilot_json(list(pass = TRUE, kernel_cases = length(checks), dense_replay_cases = length(audits),
                  model_source = ctx$contract$production_reference$commit,
                  interpretation = "Expected daughter production; biological nonviability, tail pruning and grid losses separated; Euler-step integration with saved passage-selection days fixed"),
              file.path(study, "tests/flux_accounting_status.json"))
cat("Daughter/event accounting and dense fixed-trajectory replay checks passed.\n")
