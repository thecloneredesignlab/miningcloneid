#!/usr/bin/env Rscript
# Reuse validated, read-only dynamic caches; never refit or write to fit_root.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", arg[1L])))
root <- normalizePath(file.path(here, "../.."))
source(file.path(here, "util/analysis/supp_figure4_4.R"))
output <- file.path(root, "data/Figures/Supp_Figure4_4")
fit_root <- normalizePath(Sys.getenv("SUPP4_4_FIT_ROOT"), mustWork = TRUE)
supp4_4_assert(basename(root) == "iteration5", "Only iteration5 is supported.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
baseline <- file.path(root, "data/Figures/Figure4")
reference_burden <- fread(file.path(baseline, "figure4a_seed25_burden_timecourse.tsv"))
reference_ploidy <- fread(file.path(baseline, "figure4a_seed25_weighted_mean_chromosome_timecourse.tsv"))
ranking <- fread(file.path(baseline, "invivo_best_parameters_500seeds.tsv"))
ref_provenance <- fread(file.path(baseline, "figure4a_runtime_config_provenance.tsv"))
supp4_4_assert(nrow(ranking) == 500L && identical(sort(ranking$seed_number), 1:500), "Incorrect 500-seed input table.")
key <- c("harvest", "cohort", "dose", "day")
setorderv(reference_burden, key)
setorderv(reference_ploidy, key)
supp4_4_assert(identical(reference_burden[, ..key], reference_ploidy[, ..key]), "Figure4A reference keys disagree.")
supp4_4_assert(nrow(reference_burden) == 568L && all(reference_burden$dose == 0), "Unexpected Figure4A reference scope.")
model_root <- file.path(root, "Code/vendor/O2_supply_demand_MAP_4491cd")
supp4_4_assert(unname(tools::md5sum(file.path(model_root, "model/model_O2_supply_demand_MAP.R"))) == ref_provenance$model_r_md5,
               "Model R source differs from Figure4A.")
supp4_4_assert(unname(tools::md5sum(file.path(model_root, "model/model_O2_supply_demand_MAP.cpp"))) == ref_provenance$model_cpp_md5,
               "Model C++ source differs from Figure4A.")
files <- list()
provenance <- list()
audit <- list()
seed25_errors <- numeric()
for (seed in seq_len(500L)) {
  fit <- file.path(fit_root, paste0("seed", seed))
  sim <- file.path(fit, "simulation/invivo")
  paths <- c(best_params = file.path(fit, "best_params.tsv"),
             fit_config = file.path(fit, "fit_config.rds"),
             manifest = file.path(sim, "simulation_manifest.tsv"),
             parameters_used = file.path(sim, "parameters_used.tsv"),
             burden = file.path(sim, "burden_timecourse.tsv"),
             ploidy = file.path(sim, "ploidy_weighted_mean_timecourse.tsv"))
  supp4_4_assert(all(file.exists(paths)), paste("Incomplete cache for seed", seed))
  manifest <- fread(paths[["manifest"]])
  meta <- setNames(as.character(manifest$value), manifest$key)
  cfg <- readRDS(paths[["fit_config"]])
  supp4_4_assert(meta[["status"]] == "complete" && meta[["report_dt"]] == "1" &&
                 meta[["N_UNIT"]] == "22" && meta[["start_with"]] == "chr_number" &&
                 normalizePath(meta[["fit_dir"]]) == normalizePath(fit), paste("Invalid cache manifest for seed", seed))
  supp4_4_assert(cfg$N_UNIT == 22 && cfg$DT == 0.05 && cfg$start_with == "chr_number",
                 paste("Unexpected numerical/unit configuration for seed", seed))
  if (seed == 25L) supp4_4_assert(unname(tools::md5sum(paths[["fit_config"]])) == ref_provenance$source_fit_config_md5,
                                 "Seed25 source configuration differs from Figure4A.")
  best <- fread(paths[["best_params"]])
  used <- fread(paths[["parameters_used"]])
  supp4_4_assert(!anyDuplicated(best$parameter) && !anyDuplicated(used$parameter), "Duplicate parameter definitions.")
  common <- intersect(best$parameter, used$parameter)
  b <- suppressWarnings(as.numeric(best$value[match(common, best$parameter)]))
  u <- suppressWarnings(as.numeric(used$value[match(common, used$parameter)]))
  numeric <- is.finite(b) & is.finite(u)
  supp4_4_assert(sum(numeric) >= 19L && all(abs(b[numeric] - u[numeric]) <= 1e-11 * pmax(1, abs(b[numeric]))),
                 paste("Cached parameters differ from best_params for seed", seed))
  param_names <- intersect(names(ranking), best$parameter)
  expected <- as.numeric(unlist(ranking[seed_number == seed, ..param_names]))
  actual <- as.numeric(best$value[match(param_names, best$parameter)])
  supp4_4_assert(length(param_names) >= 18L && all(abs(expected - actual) <= 1e-11 * pmax(1, abs(expected))),
                 paste("Parameters differ from Figure4 500-seed table for seed", seed))
  times <- file.info(paths)$mtime
  supp4_4_assert(max(times[1:2]) <= min(times[5:6]) + 1, paste("Cache predates fitted inputs for seed", seed))
  burden <- fread(paths[["burden"]])
  ploidy <- fread(paths[["ploidy"]])
  setorderv(burden, key)
  setorderv(ploidy, key)
  supp4_4_assert(!anyDuplicated(burden[, ..key]) && !anyDuplicated(ploidy[, ..key]) &&
                 identical(burden[, ..key], reference_burden[, ..key]) &&
                 identical(ploidy[, ..key], reference_ploidy[, ..key]), paste("Incomplete/inconsistent dynamic grid for seed", seed))
  totals <- burden$pred_burden_live_volume_mm3 + burden$pred_burden_dead_hypoxia_volume_mm3 + burden$pred_burden_dead_buffer_volume_mm3
  supp4_4_assert(all(abs(totals - burden$pred_burden_volume_mm3) <= 1e-8 * pmax(1, abs(totals))), "Burden components do not sum to total.")
  trajectory <- data.table(seed_number = seed, burden[, ..key], weighted_mean_N = ploidy$weighted_mean_N,
                           ploidy = ploidy$weighted_mean_N / 22,
                           effective_o2_pct = burden$pred_o2_pct,
                           tumor_size_mm3 = burden$pred_burden_volume_mm3)
  supp4_4_assert(all(is.finite(as.matrix(trajectory[, ..supp4_4_metrics]))) &&
                 all(trajectory$ploidy >= 1 - 1e-9 & trajectory$ploidy <= 7 + 1e-9) &&
                 all(trajectory$effective_o2_pct >= 0 & trajectory$effective_o2_pct <= 5) &&
                 all(trajectory$tumor_size_mm3 >= 0), paste("Invalid dynamic values for seed", seed))
  if (seed == 25L) {
    seed25_errors <- c(ploidy = max(abs(trajectory$ploidy - reference_ploidy$weighted_mean_N / 22)),
                       effective_o2_pct = max(abs(trajectory$effective_o2_pct - reference_burden$pred_o2_pct)),
                       tumor_size_mm3 = max(abs(trajectory$tumor_size_mm3 - reference_burden$pred_burden_volume_mm3)))
    supp4_4_assert(all(seed25_errors < 1e-8), "Seed25 cache does not reproduce packaged Figure4A.")
  }
  files[[seed]] <- trajectory
  provenance[[seed]] <- data.table(seed_number = seed, role = names(paths), path = unname(paths),
                                  md5 = unname(tools::md5sum(paths)), size_bytes = file.info(paths)$size)
  audit[[seed]] <- data.table(seed_number = seed, status = "validated_cached_dynamic_trajectory",
                             n_tumors = uniqueN(trajectory$harvest), n_rows = nrow(trajectory),
                             parameter_count_matched = sum(numeric), N_UNIT = cfg$N_UNIT, DT = cfg$DT)
  if (seed %% 25L == 0L) cat("Validated", seed, "of 500 cached endpoints\n")
}
tumor <- rbindlist(files)
result <- supp4_4_aggregate(tumor)
supp4_4_assert(nrow(tumor) == 284000L && nrow(result$cohort) == 74500L && nrow(result$summary) == 447L,
               "Unexpected output cardinality.")
supp4_4_assert(all(result$cohort[cohort == "4N" & day <= 77, n_tumors] == 4L) &&
                 all(result$cohort[cohort == "4N" & day >= 78, n_tumors] == 2L), "Unexpected 4N membership change.")
observed <- reference_burden[is.finite(obs_burden), .(observed_mean_mm3 = mean(obs_burden), n_tumors = uniqueN(harvest)), by = .(cohort, day)]
fwrite(tumor, file.path(output, "seed_tumor_trajectories.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(result$cohort, file.path(output, "seed_cohort_trajectories.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(result$summary, file.path(output, "timecourse_median_iqr.tsv"), sep = "\t")
fwrite(observed, file.path(output, "observed_tumor_size.tsv"), sep = "\t")
fwrite(rbindlist(provenance), file.path(output, "source_provenance.tsv"), sep = "\t")
fwrite(rbindlist(audit), file.path(output, "seed_validation.tsv"), sep = "\t")
validation <- data.table(metric = c("n_seeds", "n_tumor_trajectories", "n_tumor_time_rows", "n_seed_cohort_time_rows", "n_summary_rows",
                                   "seed25_ploidy_max_error", "seed25_o2_max_error", "seed25_burden_max_error", "quantile_type", "aggregation", "cache_mode"),
                          value = c(500, 4000, nrow(tumor), nrow(result$cohort), nrow(result$summary), seed25_errors, 7,
                                    "equal-tumor mean within seed/cohort/day; Q25/Q50/Q75 across 500 seeds", "read-only existing dynamic caches; no refitting"))
fwrite(validation, file.path(output, "analysis_validation.tsv"), sep = "\t")
cat("supp_figure4_4_data_validation_ok\n")
