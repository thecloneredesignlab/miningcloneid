#!/usr/bin/env Rscript
# Read archived tumor scale settings and test the audited observation equation.
# This script performs no simulation, optimization or scheduler submission.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: audit_tumor_volume_scale.R <project> <iteration5> <output-directory>")
project <- normalizePath(args[1], mustWork = TRUE)
figures <- normalizePath(args[2], mustWork = TRUE)
output <- args[3]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
figure4 <- file.path(figures, "data/Figures/Figure4")
cfg <- readRDS(file.path(figure4, "source_fit_config/fit_config.rds"))
parameters <- read.delim(file.path(figure4, "runtime_fit/seed25/best_params.tsv"))
rho <- parameters$value[parameters$parameter == "rho_2N"]
stopifnot(length(rho) == 1L, is.finite(rho), rho > 0)
shared <- file.path(project, "oxygen/results/structural_identifiability/20261109_seed144_wp0/source/code",
                    "oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_shared.R")
# Evaluate only the two pure conversion functions from the frozen source.
functions <- c("cell_volume_mm3_by_ploidy", "cell_volume_mm3_by_chr_number")
env <- new.env(parent = baseenv())
for (expression in parse(shared)) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      is.symbol(expression[[2]]) && as.character(expression[[2]]) %in% functions)
    eval(expression, envir = env)
}
stopifnot(all(vapply(functions, exists, logical(1), envir = env, inherits = FALSE)))
N <- 22:154
v1 <- env$cell_volume_mm3_by_chr_number(N, list(rho_2N = rho), cfg)
v2 <- env$cell_volume_mm3_by_chr_number(N, list(rho_2N = 2 * rho), cfg)
stopifnot(all(v1 == 1 / rho), all(v2 == v1 / 2),
          all(env$cell_volume_mm3_by_ploidy(N / 22, list(rho_2N = rho), cfg) == v1))
live <- as.numeric(seq_along(N))
dead <- rev(live) / 3
V <- sum((live + dead) * v1)
f <- sum(dead * v1) / V
stopifnot(isTRUE(all.equal(sum(live), (1 - f) * rho * V, tolerance = 1e-12)),
          isTRUE(all.equal(sum((live + dead) * v2), V / 2, tolerance = 1e-12)))
summary <- read.delim(file.path(figure4, "runtime_fit/seed25/fit_summary.tsv"))
score <- setNames(as.list(summary$value), summary$metric)
table <- read.csv(file.path(figure4, "invivo_parameter_table_seed25.csv"))
selected <- c("init_total_size", "o2_Nref", "Crowding", "o2_burden_feedback", "harvest_init_multiplier",
              "rho_2N_init", "rho_2N_min", "rho_2N_max", "use_soft_prior", "lambda_prior",
              "prior_center_log10_rho_2N", "prior_sd_log10_rho_2N", "use_necrosis_loss")
z <- (log10(rho) - cfg$prior_center_log10_rho_2N) / cfg$prior_sd_log10_rho_2N
result <- list(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  scope = "Archived iteration5 seed25 tumor volume-scale review; no tumor baseline replay",
  model_calls = 0L, optimizations = 0L, jobs_submitted = 0L,
  archived_density_cells_per_mm3 = rho, archived_configuration = cfg[intersect(selected, names(cfg))],
  archived_density_parameter_table = table[table$param_name == "log10_rho_2N", , drop = FALSE],
  prior_center_cells_per_mm3 = 10^cfg$prior_center_log10_rho_2N,
  prior_center_below_production_lower_bound = 10^cfg$prior_center_log10_rho_2N < cfg$rho_2N_min,
  density_prior_raw_at_archived_endpoint = 0.5 * z^2,
  density_prior_weighted_at_archived_endpoint = cfg$lambda_prior * 0.5 * z^2,
  archived_score = score[intersect(names(score), c("objective", "objective_data", "objective_prior", "objective_prior_raw"))],
  observation_equation_checks = list(pass = TRUE, chromosome_states = length(N),
    constant_volume_across_states = TRUE, inverse_density_scaling = TRUE,
    live_count_identity_in_model_convention = TRUE,
    limitation = "Pure conversion and accounting checks; not a full-model structural proof or a fitted profile"),
  unchanged_future_scope = "Culture targets remain p_misseg and buffer_beta; tumor scale fitting needs a separate approved phase"
)
jsonlite::write_json(result, file.path(output, "volume_scale_audit.json"), pretty = TRUE, auto_unbox = TRUE, digits = 17)
cat("Archived tumor volume-scale audit passed; no model runs or new jobs.\n")
