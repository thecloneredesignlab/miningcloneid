#!/usr/bin/env Rscript
# Extend the same eight fitted untreated tumor scenarios; no new fitting.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", arg[1L])))
root <- normalizePath(file.path(here, "../.."))
source(file.path(here, "util/analysis/supp_figure4_4.R"))
supp4_4_assert(basename(root) == "iteration5", "Only iteration5 is supported.")
old_dir <- file.path(root, "data/Figures/Supp_Figure4_4")
out <- file.path(root, "data/Figures/Supp_Figure4_4_1000days")
fit_root <- normalizePath(Sys.getenv("SUPP4_4_FIT_ROOT"), mustWork = TRUE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
short_path <- file.path(old_dir, "seed_tumor_trajectories.tsv.gz")
short <- as.data.table(read.delim(gzfile(short_path), check.names = FALSE))
setkey(short, seed_number)
supp4_4_assert(nrow(short) == 284000L && uniqueN(short$seed_number) == 500L, "Incomplete original short-window data.")
old_prov <- fread(file.path(old_dir, "source_provenance.tsv"))
old_checks <- fread(file.path(old_dir, "seed_validation.tsv"))
supp4_4_assert(nrow(old_checks) == 500L && all(old_checks$status == "validated_cached_dynamic_trajectory"), "Original seed validation is incomplete.")
key <- c("harvest", "cohort", "dose", "day")
harvests <- unique(short[, .(harvest, cohort, dose)])
supp4_4_assert(nrow(harvests) == 8L && all(harvests$dose == 0), "Expected eight untreated fitted tumor scenarios.")
grid <- harvests[, .(day = 0:1000), by = .(harvest, cohort, dose)]
setorderv(grid, key)
trajectories <- provenance <- checks <- overlaps <- vector("list", 500L)
for (seed in 1:500) {
  fit <- file.path(fit_root, paste0("seed", seed))
  sim <- file.path(fit, "simulation/invivo")
  paths <- c(best_params = file.path(fit, "best_params.tsv"), fit_config = file.path(fit, "fit_config.rds"),
             manifest = file.path(sim, "simulation_manifest.tsv"), parameters_used = file.path(sim, "parameters_used.tsv"),
             burden = file.path(sim, "predict_burden_0_1000day.tsv"),
             ploidy = file.path(sim, "predict_ploidy_weighted_mean_0_1000day.tsv"))
  supp4_4_assert(all(file.exists(paths)), paste("Missing forecast cache for seed", seed))
  hashes <- unname(tools::md5sum(paths))
  for (role_name in names(paths)[1:4]) {
    old_hash <- old_prov[seed_number == seed & role == role_name, md5]
    supp4_4_assert(length(old_hash) == 1L && hashes[match(role_name, names(paths))] == old_hash,
                   paste("Fitted inputs/cache metadata changed for seed", seed, role_name))
  }
  meta_dt <- fread(paths[["manifest"]])
  meta <- setNames(as.character(meta_dt$value), meta_dt$key)
  supp4_4_assert(meta[["status"]] == "complete" && meta[["predict_report_dt"]] == "1" &&
                   "1000" %in% strsplit(meta[["predict_horizons"]], ",", fixed = TRUE)[[1]] &&
                   meta[["N_UNIT"]] == "22" && meta[["start_with"]] == "chr_number" &&
                   normalizePath(meta[["fit_dir"]]) == normalizePath(fit), "Invalid forecast manifest.")
  # Physical measures are doubles, not integer identifiers. Long-horizon
  # integer-like values can exceed fread's integer64 inference; parse their
  # original text explicitly and reject every invalid/non-finite conversion.
  burden_numeric <- c("pred_burden_live_volume_mm3", "pred_burden_dead_hypoxia_volume_mm3",
                       "pred_burden_dead_buffer_volume_mm3", "pred_burden_volume_mm3", "pred_o2_pct")
  burden <- fread(paths[["burden"]], colClasses = list(character = burden_numeric))
  ploidy <- fread(paths[["ploidy"]], colClasses = list(character = "weighted_mean_N"))
  for (column in burden_numeric) {
    value <- suppressWarnings(as.numeric(burden[[column]]))
    supp4_4_assert(all(is.finite(value)), paste("Invalid numeric forecast value for seed", seed, column))
    set(burden, j = column, value = value)
  }
  value <- suppressWarnings(as.numeric(ploidy$weighted_mean_N))
  supp4_4_assert(all(is.finite(value)), paste("Invalid numeric forecast ploidy for seed", seed))
  set(ploidy, j = "weighted_mean_N", value = value)
  setorderv(burden, key)
  setorderv(ploidy, key)
  supp4_4_assert(!anyDuplicated(burden[, ..key]) && !anyDuplicated(ploidy[, ..key]) &&
                   identical(burden[, ..key], grid) && identical(ploidy[, ..key], grid),
                 paste("Forecast grid must contain 8 tumors x 1001 daily points for seed", seed))
  totals <- burden$pred_burden_live_volume_mm3 + burden$pred_burden_dead_hypoxia_volume_mm3 + burden$pred_burden_dead_buffer_volume_mm3
  supp4_4_assert(all(is.finite(totals)) && all(abs(totals - burden$pred_burden_volume_mm3) <= 1e-8 * pmax(1, abs(totals))), "Invalid forecast burden decomposition.")
  dt <- data.table(seed_number = seed, burden[, ..key], weighted_mean_N = ploidy$weighted_mean_N,
                   ploidy = ploidy$weighted_mean_N / 22, effective_o2_pct = burden$pred_o2_pct,
                   tumor_size_mm3 = burden$pred_burden_volume_mm3)
  supp4_4_assert(all(is.finite(as.matrix(dt[, ..supp4_4_metrics]))) &&
                   all(dt$ploidy >= 1 - 1e-9 & dt$ploidy <= 7 + 1e-9) &&
                   all(dt$effective_o2_pct >= 0 & dt$effective_o2_pct <= 5) && all(dt$tumor_size_mm3 >= 0),
                 paste("Invalid forecast metric for seed", seed))
  comparison <- merge(dt, short[.(seed)], by = c("seed_number", key), suffixes = c("_long", "_short"))
  supp4_4_assert(nrow(comparison) == 568L, "Forecast does not cover all original observed-window rows.")
  errors <- vapply(supp4_4_metrics, function(metric) {
    x <- comparison[[paste0(metric, "_long")]]
    y <- comparison[[paste0(metric, "_short")]]
    max(abs(x - y) / pmax(1, abs(y)))
  }, numeric(1))
  supp4_4_assert(all(errors <= 1e-10), paste("Long-horizon simulation differs within original follow-up for seed", seed))
  trajectories[[seed]] <- dt
  provenance[[seed]] <- data.table(seed_number = seed, role = names(paths), path = unname(paths), md5 = hashes, size_bytes = file.info(paths)$size)
  checks[[seed]] <- data.table(seed_number = seed, status = "validated_1000day_cache", n_rows = nrow(dt), n_tumors = uniqueN(dt$harvest), n_days = uniqueN(dt$day))
  overlaps[[seed]] <- data.table(seed_number = seed, n_overlap_rows = nrow(comparison),
                                ploidy_max_relative_error = errors[[1]], o2_max_relative_error = errors[[2]], volume_max_relative_error = errors[[3]])
  if (seed %% 25L == 0L) cat("Validated 1000-day cache:", seed, "/500\n")
}
tumor <- rbindlist(trajectories)
result <- supp4_4_aggregate(tumor)
supp4_4_assert(nrow(tumor) == 4004000L && nrow(result$cohort) == 1001000L && nrow(result$summary) == 6006L,
               "Unexpected long-horizon output cardinality.")
supp4_4_assert(all(result$cohort$n_tumors == 4L) && all(result$cohort$segment == 1L), "Forecasts must continue all four tumors per cohort.")
diagnostic <- result$cohort[day %in% c(0, 56, 77, 91, 100, 250, 500, 750, 1000), {
  q <- quantile(ploidy, c(.05, .25, .5, .75, .95), type = 7)
  list(n_seeds = uniqueN(seed_number), minimum = min(ploidy), q05 = q[1], q25 = q[2], median = q[3], q75 = q[4], q95 = q[5], maximum = max(ploidy), iqr = q[4] - q[2])
}, by = .(cohort, day)]
fwrite(tumor, file.path(out, "seed_tumor_trajectories.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(result$cohort, file.path(out, "seed_cohort_trajectories.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(result$summary, file.path(out, "timecourse_median_iqr.tsv"), sep = "\t")
fwrite(fread(file.path(old_dir, "observed_tumor_size.tsv")), file.path(out, "observed_tumor_size.tsv"), sep = "\t")
context_paths <- c(short_path, file.path(old_dir, c("source_provenance.tsv", "seed_validation.tsv", "timecourse_median_iqr.tsv", "observed_tumor_size.tsv")))
context <- data.table(seed_number = 0L, role = paste0("short_window_", basename(context_paths)), path = context_paths,
                      md5 = unname(tools::md5sum(context_paths)), size_bytes = file.info(context_paths)$size)
fwrite(rbind(rbindlist(provenance), context), file.path(out, "source_provenance.tsv"), sep = "\t")
fwrite(rbindlist(checks), file.path(out, "seed_validation.tsv"), sep = "\t")
fwrite(rbindlist(overlaps), file.path(out, "overlap_validation.tsv"), sep = "\t")
fwrite(diagnostic, file.path(out, "ploidy_time_window_comparison.tsv"), sep = "\t")
overlap_table <- rbindlist(overlaps)
max_overlap_error <- max(as.matrix(overlap_table[, .(ploidy_max_relative_error, o2_max_relative_error, volume_max_relative_error)]))
fwrite(data.table(metric = c("n_seeds", "horizon_days", "n_tumor_time_rows", "n_cohort_time_rows", "n_summary_rows",
                             "n_tumors_per_cohort", "quantile_type", "max_original_overlap_relative_error", "extrapolation", "fit_quality_filter"),
                  value = c(500, 1000, nrow(tumor), nrow(result$cohort), nrow(result$summary), 4, 7,
                            max_overlap_error,
                            "all original tumor scenarios continued beyond observed follow-up", "none")),
       file.path(out, "analysis_validation.tsv"), sep = "\t")
cat("supp_figure4_4_1000days_data_validation_ok\n")
