#!/usr/bin/env Rscript
# Audit all optimizer endpoints with the same canonical operator used by the evaluator.
argv <- commandArgs(TRUE)
if (length(argv) != 3L) stop("Usage: audit_neighborhood_inputs.R FIT_ROOT FIGURE4_DIR OUT_JSON")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
source(file.path(dirname(normalizePath(script)), "../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R"), local = TRUE)
suppressPackageStartupMessages(library(jsonlite))
seeds <- paste0("seed", seq_len(500L))
reference <- read.delim(file.path(argv[[2L]], "fixed_o2_dominant_ploidy_201grid.tsv"))
active <- c("lam_max", "p_mis_base", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
            "buffer_n_exp", "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O")
# Includes all non-fitted settings used to construct the fixed-oxygen operator.
cfg_keys <- c("N_MIN", "N_MAX", "N_UNIT", "chr_lengths_bp", "start_with", "O2_growth",
              "ploidy_O2_death", "o2_min", "boundary", "Crowding", "K", "crowding")
audits <- vector("list", length(seeds))
baseline <- NULL
for (i in seq_along(seeds)) {
  seed <- seeds[[i]]
  directory <- file.path(argv[[1L]], seed)
  cfg_raw <- readRDS(file.path(directory, "fit_config.rds"))
  if (!identical(as.integer(cfg_raw$seed), as.integer(i))) stop("Config seed does not match directory: ", seed)
  p <- read.delim(file.path(directory, "best_params.tsv"))
  if (anyDuplicated(p$parameter) || !all(active %in% p$parameter) || any(!is.finite(p$value))) stop("Invalid endpoint: ", seed)
  params <- as.list(setNames(p$value, p$parameter))
  rp <- prepare_run_params(params, "invivo", cfg_raw, 0)
  cfg <- prepare_sim_cfg(cfg_raw, list(), 0, rp)
  selected <- cfg[cfg_keys[cfg_keys %in% names(cfg)]]
  selected$boundary <- rp$boundary
  selected$alpha <- rp$alpha
  selected$gamma <- rp$gamma
  if (cfg$N_MIN != 22L || cfg$N_MAX != 154L || cfg$N_UNIT != 22L || !isTRUE(cfg$O2_growth) ||
      !identical(cfg$ploidy_O2_death, "ploidy_related")) stop("Unexpected model configuration: ", seed)
  if (is.null(baseline)) baseline <- selected
  if (!isTRUE(all.equal(selected, baseline, tolerance = 0))) stop("Operator configuration differs: ", seed)
  max_error <- c(ploidy = 0, growth = 0)
  for (o2 in c(0, 2.5, 5)) {
    actual <- fixo2_dominant_attractor_one(seed, rp, globalenv(), cfg, o2)
    expected <- reference[reference$seed_id == seed & abs(reference$O2_pct - o2) < 1e-12, ]
    if (nrow(expected) != 1L || actual$status != "ok" || !isTRUE(actual$eigenvector_nonnegative)) stop("Invalid reference: ", seed)
    errors <- abs(c(actual$dominant_mean_ploidy - expected$dominant_mean_ploidy,
                    actual$dominant_growth_rate - expected$dominant_growth_rate))
    if (any(!is.finite(errors)) || any(errors > 1e-8)) stop("Figure 4 reference mismatch: ", seed, " O2=", o2)
    max_error <- pmax(max_error, errors)
  }
  audits[[i]] <- list(fit_seed = seed, config_seed = cfg_raw$seed,
                     reference_max_error = as.list(max_error), status = "ok")
  if (i %% 25L == 0L) message("Audited ", i, "/500 fitted seeds")
}
dir.create(dirname(argv[[3L]]), recursive = TRUE, showWarnings = FALSE)
write_json(list(n_seeds = length(seeds), operator_config = baseline, seeds = audits,
                reference_oxygen = c(0, 2.5, 5), reference_tolerance = 1e-8),
           argv[[3L]], auto_unbox = TRUE, pretty = TRUE, digits = 17)
