# Metadata only. Never source the model, construct a kernel, or call a simulator.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
run <- normalizePath(args[[1L]], mustWork = TRUE)
stopifnot(grepl("invivo_structural_", basename(run), fixed = TRUE))
x <- readRDS(file.path(run, "inputs", "old_contract", "contract.rds"))
out <- file.path(run, "tables")
dir.create(out, showWarnings = FALSE)
write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "NA")
}
cfg <- x$cfg
scalar <- vapply(cfg, function(z) is.atomic(z) && length(z) == 1L, logical(1))
write_tsv(data.frame(setting = names(cfg)[scalar],
                     value = vapply(cfg[scalar], as.character, character(1)),
                     r_type = vapply(cfg[scalar], typeof, character(1))), "production_scalar_settings.tsv")
pt <- read.delim(file.path(run, "inputs", "old_contract", "transformed_parameters.tsv"), check.names = FALSE)
stopifnot(sum(pt$estimate) == 20L, nrow(pt) == 21L)
write_tsv(pt, "production_parameter_table.tsv")
scenarios <- do.call(rbind, lapply(x$scenarios, function(s) data.frame(
  scenario = s$harvest, cohort = s$cohort, dose = s$dose,
  harvest_day = s$harvest_day, burden_raw = length(s$obs_burden),
  burden_positive = sum(is.finite(s$obs_burden) & s$obs_burden > 0),
  burden_zero = sum(s$obs_burden == 0, na.rm = TRUE),
  burden_nonfinite = sum(!is.finite(s$obs_burden)),
  chromosome_records = length(s$endpoint_obs_z),
  necrosis_fraction = s$obs_necrosis_fraction, necrosis_record_count = s$n_necrosis_obs,
  necrosis_day = s$necrosis_day)))
write_tsv(scenarios, "scenario_inventory.tsv")
burden <- do.call(rbind, lapply(x$scenarios, function(s) data.frame(
  scenario = s$harvest, day = s$obs_days, observed_volume_mm3 = s$obs_burden,
  current_score_includes = is.finite(s$obs_burden) & s$obs_burden > 0)))
write_tsv(burden, "burden_sampling_schedule.tsv")
stopifnot(nrow(burden) == 120L, sum(scenarios$burden_positive) == x$expected_counts$burden,
          sum(scenarios$chromosome_records) == x$expected_counts$ploidy,
          sum(scenarios$necrosis_record_count) == x$expected_counts$necrosis)
# This records the exact source helper's inputs, without executing that helper.
chr <- cfg$chr_lengths_bp
write_tsv(data.frame(index = seq_along(chr), chromosome_length_bp = as.numeric(chr)), "initial_chromosome_lengths.tsv")
core <- cfg$model_core
write_tsv(data.frame(chromosome_state = core$grid_pre,
                     initial_2N = core$init_state_2N,
                     initial_4N = core$init_state_4N), "archived_initial_live_states.tsv")
stopifnot(length(core$grid_pre) == 133L, min(core$grid_pre) == 22L,
          max(core$grid_pre) == 154L,
          sum(core$init_state_2N) == cfg$init_total_size,
          sum(core$init_state_4N) == cfg$init_total_size,
          sum(core$init_state_2N > 0) == 1L,
          sum(core$init_state_4N > 0) == 1L)
write_tsv(data.frame(key = c("native_model_calls", "symbolic_calls", "estimated_coordinates", "fixed_parameter_rows", "burden_raw", "burden_included", "chromosome_records", "necrosis_records"),
                     value = c(0L, 0L, sum(pt$estimate), sum(!pt$estimate), nrow(burden), sum(burden$current_score_includes), sum(scenarios$chromosome_records), sum(scenarios$necrosis_record_count))), "contract_export_checks.tsv")
cat("Metadata export passed; zero model calls.\n")
