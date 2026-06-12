#!/usr/bin/env Rscript

# Post-fit in vivo simulations from an O2 joint/in vivo seed directory.
#
# Default design:
#   1. continuous initial ploidy distribution centered at ploidy 3.0;
#   2. O2 initial-supply sweep across five o2_S0 values and five initial burdens;
#   3. event-triggered treatment sweep: grow from initial_burden_cells until a
#      total-cell burden threshold is reached, then increase p_mis_base;
#   4. p_wgd sweep across five p_wgd values and five initial burdens.
#
# Example:
#   Rscript oxygen/code/O2_supply_demand_MAP/simulation/simulate_invivo_mixed_ploidy_perturbations.R \
#     --fit_dir=oxygen/results/softcouple_sigma1p5_warmup_seed50_seed350_500seed/seed330

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

.o2_sim_bootstrap_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})

SCRIPT_DIR <- normalizePath(.o2_sim_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())

MODEL_PATH <- file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R")
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(MODEL_PATH))
source(MODEL_PATH, local = environment())
rm(.o2_sim_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_bool <- o2sd_as_bool

resolve_path_value <- function(path_value, base_dir = getwd()) {
  txt <- path_value
  if (is.null(txt) || !length(txt)) return(NULL)
  txt <- trimws(as.character(txt[[1]]))
  if (!nzchar(txt)) return(NULL)
  if (startsWith(txt, "~")) return(normalizePath(path.expand(txt), mustWork = FALSE))
  if (grepl("^(/|[A-Za-z]:[/\\\\])", txt)) return(normalizePath(txt, mustWork = FALSE))
  normalizePath(file.path(base_dir, txt), mustWork = FALSE)
}

parse_num_list <- function(x, default) {
  if (is.null(x) || !length(x) || !nzchar(trimws(as.character(x[[1]])))) {
    return(as.numeric(default))
  }
  txt <- gsub("[; ]+", ",", trimws(as.character(x[[1]])))
  vals <- suppressWarnings(as.numeric(strsplit(txt, ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) stop("Failed to parse numeric list: ", as.character(x[[1]]))
  vals
}

safe_id <- function(...) {
  x <- paste(..., sep = "_")
  x <- gsub("\\+", "", x)
  x <- gsub("-", "m", x)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_]+", "", x)
  x
}

num_tag <- function(x) {
  safe_id(sprintf("%.4g", as.numeric(x)))
}

read_required_tsv <- function(path) {
  if (!file.exists(path)) stop("Required file was not found: ", path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

read_fit_config <- function(fit_dir) {
  path <- file.path(fit_dir, "fit_config.rds")
  if (!file.exists(path)) stop("fit_config.rds was not found: ", path)
  readRDS(path)
}

read_run_params <- function(fit_dir, cfg) {
  path <- file.path(fit_dir, "best_params.tsv")
  tab <- read_required_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain columns: parameter, value")
  }
  vals <- suppressWarnings(as.numeric(tab$value))
  if (any(!is.finite(vals))) {
    stop("Non-finite values were found in ", path)
  }
  rp <- as.list(setNames(vals, as.character(tab$parameter)))
  rp$boundary <- "drop"
  rp
}

prepare_sim_cfg <- function(cfg, argv) {
  cfg$DT <- as_num(argv$dt, as.numeric(.first_non_null_local(cfg$DT, 0.05)))
  cfg$fit_treatment <- FALSE
  cfg$dose_ref <- as.numeric(.first_non_null_local(cfg$dose_ref, 30))
  cfg$tx_mult_min <- as.numeric(.first_non_null_local(cfg$tx_mult_min, 0.05))
  cfg$Crowding <- as_bool(argv$Crowding, isTRUE(.first_non_null_local(cfg$Crowding, FALSE)))
  cfg$K <- as.numeric(.first_non_null_local(cfg$K, 1e12))
  cfg$crowding <- as.character(.first_non_null_local(cfg$crowding, "logistic"))
  cfg$min_pop <- as.numeric(.first_non_null_local(cfg$min_pop, 1e-12))
  cfg$o2_burden_feedback <- as_bool(argv$o2_burden_feedback, isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE)))
  cfg$O2_growth <- as_bool(argv$O2_growth, isTRUE(.first_non_null_local(cfg$O2_growth, TRUE)))
  cfg$o2_cache_bin_pct <- as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01))
  cfg$o2_cache_hysteresis_pct <- as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005))
  cfg$o2_cache_profile <- FALSE
  cfg$burden_log_eps <- as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12))
  cfg$N_UNIT <- as.integer(.first_non_null_local(cfg$N_UNIT, 22L))
  cfg$N_MIN <- as.integer(.first_non_null_local(cfg$N_MIN, 22L))
  cfg$N_MAX <- as.integer(.first_non_null_local(cfg$N_MAX, 154L))
  cfg$chr_lengths_bp <- .first_non_null_local(cfg$chr_lengths_bp, default_chr_lengths_bp_1to22())
  cfg$start_with <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "chr_number"))
  cfg$ploidy_O2_death <- assert_canonical_ploidy_o2_death_mode(.first_non_null_local(cfg$ploidy_O2_death, "ploidy_related"))
  cfg$o2_S0_upper_bound <- as.numeric(.first_non_null_local(cfg$o2_S0_upper_bound, cfg$o2_S0_max, 5.0))
  cfg$o2_Nref <- as.numeric(.first_non_null_local(cfg$o2_Nref, cfg$init_total_size, 1e6))
  cfg$o2_min <- as.numeric(.first_non_null_local(cfg$o2_min, 0.0))
  cfg
}

make_continuous_ploidy_state <- function(cfg,
                                         total_live_cells,
                                         mean_ploidy = 3.0,
                                         sd_ploidy = 0.4,
                                         min_ploidy = 1.5,
                                         max_ploidy = 6.0) {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  ploidy_grid <- as.numeric(grid_pre) / as.numeric(cfg$N_UNIT)
  w <- stats::dnorm(ploidy_grid, mean = as.numeric(mean_ploidy), sd = as.numeric(sd_ploidy))
  w[ploidy_grid < min_ploidy | ploidy_grid > max_ploidy] <- 0
  if (!is.finite(sum(w)) || sum(w) <= 0) {
    stop("Initial ploidy distribution has zero mass; check mean/sd/truncation range.")
  }
  w <- w / sum(w)
  x <- as.numeric(w) * as.numeric(total_live_cells)
  names(x) <- as.character(grid_pre)
  x
}

simulate_cpp_trajectory <- function(init_state,
                                    run_params,
                                    cfg,
                                    horizon_day,
                                    report_dt,
                                    init_dead_hypoxia_state = NULL,
                                    init_dead_buffer_state = NULL,
                                    start_day = 0,
                                    scenario_id,
                                    experiment,
                                    segment = "single") {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  horizon_day <- as.numeric(horizon_day)
  report_dt <- as.numeric(report_dt)
  if (!is.finite(horizon_day) || horizon_day < 0) stop("horizon_day must be finite and >= 0.")
  if (!is.finite(report_dt) || report_dt <= 0) stop("report_dt must be finite and > 0.")

  step_max <- as.integer(round(horizon_day / cfg$DT))
  keep_steps <- unique(as.integer(round(seq(0, horizon_day, by = report_dt) / cfg$DT)))
  keep_steps <- sort(unique(c(0L, keep_steps, step_max)))
  keep_steps <- keep_steps[keep_steps >= 0L & keep_steps <= step_max]
  local_days <- as.numeric(keep_steps) * as.numeric(cfg$DT)
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg))

  o2_s0_upper <- as.numeric(.first_non_null_local(cfg$o2_S0_upper_bound, 5.0))
  o2_S0_use <- as.numeric(.first_non_null_local(run_params$o2_S0, cfg$o2_S0_init, 0.5))
  o2_S0_use <- min(max(o2_S0_use, 0), o2_s0_upper)

  sim_args <- list(
    init_state = as.numeric(init_state),
    init_dead_hypoxia_state = as.numeric(init_dead_hypoxia_state %||% rep(0, length(grid_pre))),
    init_dead_buffer_state = as.numeric(init_dead_buffer_state %||% rep(0, length(grid_pre))),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(keep_steps),
    sim_end_step = as.integer(step_max),
    DT = as.numeric(cfg$DT),
    dose = 0.0,
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = as.numeric(horizon_day + 1),
    fit_treatment = FALSE,
    alpha = 0.0,
    gamma = 1.0,
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding_enabled = isTRUE(cfg$Crowding),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0)),
    o2_feedback = isTRUE(cfg$o2_burden_feedback),
    o2_S0 = as.numeric(o2_S0_use),
    kappa_O = as.numeric(.first_non_null_local(run_params$kappa_O, cfg$kappa_O_init, 1.0)),
    tau_O2 = as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0)),
    o2_Nref = as.numeric(.first_non_null_local(cfg$o2_Nref, cfg$init_total_size, 1e6)),
    o2_min = as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.0)),
    o2_S0_upper_bound = as.numeric(o2_s0_upper),
    eta_o2 = as.numeric(.first_non_null_local(run_params$eta_o2, cfg$eta_o2_init, 1.0)),
    o2_cache_bin_pct = as.numeric(cfg$o2_cache_bin_pct),
    o2_cache_hysteresis_pct = as.numeric(cfg$o2_cache_hysteresis_pct),
    o2_cache_profile = FALSE,
    O2_growth = isTRUE(cfg$O2_growth),
    lam_max = as.numeric(run_params$lam_max),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = 1e-8,
    buffer_smax = as.numeric(.first_non_null_local(run_params$buffer_smax, cfg$buffer_smax_init, 0.8)),
    buffer_beta = as.numeric(.first_non_null_local(run_params$buffer_beta, cfg$buffer_beta_init, 1.0)),
    buffer_n_exp = as.numeric(.first_non_null_local(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1.0)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = 0.0,
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu, cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0)),
    ploidy_O2_death = as.character(cfg$ploidy_O2_death),
    start_with = as.character(cfg$start_with),
    k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(cfg$burden_log_eps),
    return_full_trajectory = TRUE
  )
  sim_cpp <- cpp_o2simps_simulate_one(sim_args)

  live_state <- as.matrix(sim_cpp$live_state_obs)
  if (!identical(dim(live_state), c(length(keep_steps), length(grid_pre)))) {
    stop("Unexpected live_state_obs shape returned by cpp_o2simps_simulate_one.")
  }
  dead_hypoxia_state <- as.matrix(sim_cpp$dead_hypoxia_state_obs)
  dead_buffer_state <- as.matrix(sim_cpp$dead_buffer_state_obs)
  if (!identical(dim(dead_hypoxia_state), c(length(keep_steps), length(grid_pre)))) {
    stop("Unexpected dead_hypoxia_state_obs shape returned by cpp_o2simps_simulate_one.")
  }
  if (!identical(dim(dead_buffer_state), c(length(keep_steps), length(grid_pre)))) {
    stop("Unexpected dead_buffer_state_obs shape returned by cpp_o2simps_simulate_one.")
  }
  live_pop <- rowSums(live_state, na.rm = TRUE)
  live_frac <- live_state / pmax(live_pop, 1e-12)
  day <- as.numeric(start_day) + local_days

  burden <- data.frame(
    scenario_id = scenario_id,
    experiment = experiment,
    segment = segment,
    local_day = local_days,
    day = day,
    step = as.integer(keep_steps),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, 1e-5)),
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    o2_S0 = as.numeric(o2_S0_use),
    pred_burden_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_total_obs, sim_cpp$Ntot_obs)),
    pred_burden_live_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_live_obs, sim_cpp$Ntot_obs)),
    pred_burden_dead_hypoxia_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_hypoxia_obs, rep(0, length(keep_steps)))),
    pred_burden_dead_buffer_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_buffer_obs, rep(0, length(keep_steps)))),
    pred_burden_dead_total_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_total_obs, rep(0, length(keep_steps)))),
    pred_burden_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_total_obs, sim_cpp$Vmm3_obs)),
    pred_burden_live_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_live_obs, sim_cpp$Vmm3_obs)),
    pred_burden_dead_hypoxia_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_hypoxia_obs, rep(0, length(keep_steps)))),
    pred_burden_dead_buffer_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_buffer_obs, rep(0, length(keep_steps)))),
    pred_burden_dead_total_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_total_obs, rep(0, length(keep_steps)))),
    pred_o2_target_pct = as.numeric(.first_non_null_local(sim_cpp$O2_target_obs, rep(NA_real_, length(keep_steps)))),
    pred_o2_pct = as.numeric(.first_non_null_local(sim_cpp$O2_eff_obs, rep(NA_real_, length(keep_steps)))),
    stringsAsFactors = FALSE
  )
  burden$pred_log10_burden_cells <- log10(pmax(burden$pred_burden_cells, 1e-300))
  burden$pred_log10_burden_live_cells <- log10(pmax(burden$pred_burden_live_cells, 1e-300))
  burden$pred_log10_burden_dead_total_cells <- log10(pmax(burden$pred_burden_dead_total_cells, 1e-300))
  burden$pred_log10_burden_volume_mm3 <- log10(pmax(burden$pred_burden_volume_mm3, 1e-300))
  burden$pred_o2_lag_gap_pct <- burden$pred_o2_target_pct - burden$pred_o2_pct

  ploidy <- data.frame(
    scenario_id = scenario_id,
    experiment = experiment,
    segment = segment,
    local_day = rep(local_days, each = length(grid_pre)),
    day = rep(day, each = length(grid_pre)),
    N = rep(as.numeric(grid_pre), times = length(keep_steps)),
    ploidy = rep(as.numeric(grid_pre) / as.numeric(cfg$N_UNIT), times = length(keep_steps)),
    fraction = as.numeric(t(live_frac)),
    live_cells = rep(as.numeric(live_pop), each = length(grid_pre)),
    stringsAsFactors = FALSE
  )

  list(
    burden = burden,
    ploidy = ploidy,
    live_state = live_state,
    dead_hypoxia_state = dead_hypoxia_state,
    dead_buffer_state = dead_buffer_state,
    keep_steps = keep_steps,
    local_days = local_days
  )
}

summarise_ploidy_timecourse <- function(ploidy_df, live_min_cells = 1e-6) {
  ploidy_df %>%
    group_by(scenario_id, experiment, day) %>%
    summarise(
      live_cells = max(live_cells, na.rm = TRUE),
      frac_sum = sum(fraction, na.rm = TRUE),
      mean_N_num = sum(N * fraction, na.rm = TRUE),
      mean_ploidy_num = sum(ploidy * fraction, na.rm = TRUE),
      frac_ploidy_lt2_raw = sum(fraction[ploidy < 2], na.rm = TRUE),
      frac_ploidy_2to3_raw = sum(fraction[ploidy >= 2 & ploidy < 3], na.rm = TRUE),
      frac_ploidy_3to4_raw = sum(fraction[ploidy >= 3 & ploidy < 4], na.rm = TRUE),
      frac_ploidy_ge4_raw = sum(fraction[ploidy >= 4], na.rm = TRUE),
      frac_ploidy_ge5_raw = sum(fraction[ploidy >= 5], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      has_live_ploidy = is.finite(live_cells) & live_cells > live_min_cells & is.finite(frac_sum) & frac_sum > 1e-9,
      mean_N = ifelse(has_live_ploidy, mean_N_num / pmax(frac_sum, 1e-12), NA_real_),
      mean_ploidy = ifelse(has_live_ploidy, mean_ploidy_num / pmax(frac_sum, 1e-12), NA_real_),
      frac_ploidy_lt2 = ifelse(has_live_ploidy, frac_ploidy_lt2_raw, NA_real_),
      frac_ploidy_2to3 = ifelse(has_live_ploidy, frac_ploidy_2to3_raw, NA_real_),
      frac_ploidy_3to4 = ifelse(has_live_ploidy, frac_ploidy_3to4_raw, NA_real_),
      frac_ploidy_ge4 = ifelse(has_live_ploidy, frac_ploidy_ge4_raw, NA_real_),
      frac_ploidy_ge5 = ifelse(has_live_ploidy, frac_ploidy_ge5_raw, NA_real_)
    ) %>%
    select(
      scenario_id, experiment, day, live_cells, mean_N, mean_ploidy,
      frac_ploidy_lt2, frac_ploidy_2to3, frac_ploidy_3to4,
      frac_ploidy_ge4, frac_ploidy_ge5
    )
}

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  invisible(path)
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  invisible(path)
}

annotate_timecourses <- function(sim, design_row) {
  add_cols <- setdiff(names(design_row), names(sim$burden))
  burden <- bind_cols(as.data.frame(design_row[rep(1, nrow(sim$burden)), add_cols, drop = FALSE]), sim$burden)
  add_cols_ploidy <- setdiff(names(design_row), names(sim$ploidy))
  ploidy <- bind_cols(as.data.frame(design_row[rep(1, nrow(sim$ploidy)), add_cols_ploidy, drop = FALSE]), sim$ploidy)
  list(burden = burden, ploidy = ploidy)
}

run_static_scenario <- function(run_params,
                                cfg,
                                init_state,
                                horizon_day,
                                report_dt,
                                design_row) {
  design_row$status <- "simulated"
  rp <- run_params
  rp$o2_S0 <- as.numeric(design_row$o2_S0)
  rp$p_wgd <- as.numeric(design_row$p_wgd)
  rp$p_mis_base <- as.numeric(design_row$p_mis_base_pre)
  sim <- simulate_cpp_trajectory(
    init_state = init_state,
    run_params = rp,
    cfg = cfg,
    horizon_day = horizon_day,
    report_dt = report_dt,
    scenario_id = design_row$scenario_id,
    experiment = design_row$experiment,
    segment = "single"
  )
  annotate_timecourses(sim, design_row)
}

run_triggered_pmiss_scenario <- function(run_params,
                                         cfg,
                                         init_state,
                                         horizon_day,
                                         report_dt,
                                         trigger_check_dt,
                                         design_row) {
  rp_pre <- run_params
  rp_pre$p_mis_base <- as.numeric(design_row$p_mis_base_pre)
  rp_pre$p_wgd <- as.numeric(design_row$p_wgd)
  rp_pre$o2_S0 <- as.numeric(design_row$o2_S0)

  pre <- simulate_cpp_trajectory(
    init_state = init_state,
    run_params = rp_pre,
    cfg = cfg,
    horizon_day = horizon_day,
    report_dt = trigger_check_dt,
    scenario_id = design_row$scenario_id,
    experiment = design_row$experiment,
    segment = "pre"
  )

  trigger <- as.numeric(design_row$trigger_burden_cells)
  hit_idx <- which(pre$burden$pred_burden_cells >= trigger)
  if (!length(hit_idx)) {
    design_row$status <- "trigger_not_reached"
    design_row$trigger_day <- NA_real_
    design_row$actual_trigger_burden_cells <- NA_real_
    ann <- annotate_timecourses(pre, design_row)
    return(list(design = design_row, burden = ann$burden, ploidy = ann$ploidy))
  }

  hit <- hit_idx[[1]]
  trigger_day <- as.numeric(pre$burden$day[[hit]])
  design_row$status <- "treated"
  design_row$trigger_day <- trigger_day
  design_row$actual_trigger_burden_cells <- as.numeric(pre$burden$pred_burden_cells[[hit]])

  rp_post <- rp_pre
  rp_post$p_mis_base <- as.numeric(design_row$p_mis_base_post)
  post_horizon <- max(0, as.numeric(horizon_day) - trigger_day)
  post_init <- as.numeric(pre$live_state[hit, ])
  post_init_dead_hypoxia <- as.numeric(pre$dead_hypoxia_state[hit, ])
  post_init_dead_buffer <- as.numeric(pre$dead_buffer_state[hit, ])

  post <- simulate_cpp_trajectory(
    init_state = post_init,
    run_params = rp_post,
    cfg = cfg,
    horizon_day = post_horizon,
    report_dt = report_dt,
    init_dead_hypoxia_state = post_init_dead_hypoxia,
    init_dead_buffer_state = post_init_dead_buffer,
    start_day = trigger_day,
    scenario_id = design_row$scenario_id,
    experiment = design_row$experiment,
    segment = "post"
  )

  pre_keep <- pre
  pre_keep$burden <- pre_keep$burden %>% filter(day < trigger_day - 1e-9)
  pre_keep$ploidy <- pre_keep$ploidy %>% filter(day < trigger_day - 1e-9)

  pre_ann <- annotate_timecourses(pre_keep, design_row)
  post_ann <- annotate_timecourses(post, design_row)
  list(
    design = design_row,
    burden = bind_rows(pre_ann$burden, post_ann$burden),
    ploidy = bind_rows(pre_ann$ploidy, post_ann$ploidy)
  )
}

build_design_rows <- function(run_params,
                              initial_burden_values,
                              o2_values,
                              p_wgd_values,
                              pmis_base_values,
                              trigger_burden_values,
                              initial_burden_cells) {
  baseline_o2 <- as.numeric(run_params$o2_S0)
  baseline_pmiss <- as.numeric(run_params$p_mis_base)
  baseline_pwgd <- as.numeric(run_params$p_wgd)

  o2_design <- expand.grid(
    initial_burden_cells = initial_burden_values,
    o2_S0 = o2_values,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      experiment = "o2_initial_supply",
      scenario_id = safe_id("o2", "init", num_tag(initial_burden_cells), "s0", num_tag(o2_S0)),
      varied_parameter = "o2_S0",
      varied_value = as.numeric(o2_S0),
      trigger_burden_cells = NA_real_,
      p_mis_base_pre = baseline_pmiss,
      p_mis_base_post = baseline_pmiss,
      p_wgd = baseline_pwgd,
      status = "planned",
      trigger_day = NA_real_,
      actual_trigger_burden_cells = NA_real_
    )

  pwgd_design <- expand.grid(
    initial_burden_cells = initial_burden_values,
    p_wgd = p_wgd_values,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      experiment = "p_wgd_sweep",
      scenario_id = safe_id("pwgd", "init", num_tag(initial_burden_cells), "pwgd", num_tag(p_wgd)),
      varied_parameter = "p_wgd",
      varied_value = as.numeric(p_wgd),
      trigger_burden_cells = NA_real_,
      o2_S0 = baseline_o2,
      p_mis_base_pre = baseline_pmiss,
      p_mis_base_post = baseline_pmiss,
      status = "planned",
      trigger_day = NA_real_,
      actual_trigger_burden_cells = NA_real_
    )

  pmiss_design <- expand.grid(
    trigger_burden_cells = trigger_burden_values,
    p_mis_base_post = pmis_base_values,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      experiment = "pmiss_triggered_treatment",
      scenario_id = safe_id("pmiss", "trigger", num_tag(trigger_burden_cells), "post", num_tag(p_mis_base_post)),
      varied_parameter = "p_mis_base_post",
      varied_value = as.numeric(p_mis_base_post),
      initial_burden_cells = as.numeric(initial_burden_cells),
      o2_S0 = baseline_o2,
      p_mis_base_pre = baseline_pmiss,
      p_wgd = baseline_pwgd,
      status = "planned",
      trigger_day = NA_real_,
      actual_trigger_burden_cells = NA_real_
    )

  bind_rows(o2_design, pwgd_design, pmiss_design) %>%
    select(
      experiment, scenario_id, varied_parameter, varied_value,
      initial_burden_cells, trigger_burden_cells, trigger_day, actual_trigger_burden_cells,
      o2_S0, p_mis_base_pre, p_mis_base_post, p_wgd, status
    )
}

plot_outputs <- function(burden_all, ploidy_summary, design, out_dir) {
  burden_plot <- burden_all %>%
    mutate(
      varied_value_label = format(signif(varied_value, 4), trim = TRUE),
      initial_burden_label = paste0("init burden=", format(signif(initial_burden_cells, 3), scientific = TRUE)),
      trigger_label = ifelse(
        is.finite(trigger_burden_cells),
        paste0("trigger burden=", format(signif(trigger_burden_cells, 3), scientific = TRUE)),
        initial_burden_label
      )
    )

  ploidy_plot <- ploidy_summary %>%
    left_join(design, by = c("scenario_id", "experiment")) %>%
    mutate(
      varied_value_label = format(signif(varied_value, 4), trim = TRUE),
      initial_burden_label = paste0("init burden=", format(signif(initial_burden_cells, 3), scientific = TRUE)),
      trigger_label = ifelse(
        is.finite(trigger_burden_cells),
        paste0("trigger burden=", format(signif(trigger_burden_cells, 3), scientific = TRUE)),
        initial_burden_label
      )
    )

  static_burden <- burden_plot %>% filter(experiment %in% c("o2_initial_supply", "p_wgd_sweep"))
  if (nrow(static_burden) > 0) {
    p <- ggplot(static_burden, aes(day, pred_log10_burden_cells, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      facet_grid(experiment ~ initial_burden_label, scales = "free_y") +
      labs(x = "Day", y = "log10 Burden (total cells)", color = "Value") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "burden_static_sweeps.pdf"), p, width = 14, height = 7)
  }

  static_ploidy <- ploidy_plot %>% filter(experiment %in% c("o2_initial_supply", "p_wgd_sweep"))
  if (nrow(static_ploidy) > 0) {
    p <- ggplot(static_ploidy, aes(day, mean_ploidy, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      facet_grid(experiment ~ initial_burden_label, scales = "free_y") +
      labs(x = "Day", y = "Mean ploidy", color = "Value") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "mean_ploidy_static_sweeps.pdf"), p, width = 14, height = 7)
  }

  tx_burden <- burden_plot %>% filter(experiment == "pmiss_triggered_treatment")
  if (nrow(tx_burden) > 0) {
    p <- ggplot(tx_burden, aes(day, pred_log10_burden_cells, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      geom_vline(aes(xintercept = trigger_day), linetype = "dashed", linewidth = 0.25, alpha = 0.4) +
      facet_wrap(~ trigger_label, scales = "free_y") +
      labs(x = "Day", y = "log10 Burden (total cells)", color = "Post p_mis_base") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "burden_triggered_pmiss_treatment.pdf"), p, width = 12, height = 8)
  }

  tx_ploidy <- ploidy_plot %>% filter(experiment == "pmiss_triggered_treatment")
  if (nrow(tx_ploidy) > 0) {
    p <- ggplot(tx_ploidy, aes(day, mean_ploidy, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      geom_vline(aes(xintercept = trigger_day), linetype = "dashed", linewidth = 0.25, alpha = 0.4) +
      facet_wrap(~ trigger_label, scales = "free_y") +
      labs(x = "Day", y = "Mean ploidy", color = "Post p_mis_base") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "mean_ploidy_triggered_pmiss_treatment.pdf"), p, width = 12, height = 8)
  }
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- resolve_path_value(argv$fit_dir, getwd())
  if (is.null(fit_dir)) {
    stop("Missing required argument: --fit_dir=/path/to/seed_dir")
  }
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  cfg <- prepare_sim_cfg(read_fit_config(fit_dir), argv)
  run_params <- read_run_params(fit_dir, cfg)

  out_dir <- resolve_path_value(argv$out_dir, getwd())
  if (is.null(out_dir)) {
    out_dir <- file.path(fit_dir, "simulation", "invivo_mixed_ploidy_perturbations")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  initial_burden_cells <- as_num(argv$initial_burden_cells, as_num(argv$initial_live_cells, 1e5))
  initial_burden_values <- parse_num_list(argv$initial_burden_values, c(1e5, 2.5e5, 5e5, 1e6, 2e6))
  trigger_burden_values <- parse_num_list(argv$trigger_burden_values, c(5e5, 1e6, 2e6, 5e6, 1e7))
  o2_values <- parse_num_list(argv$o2_values, c(0.5, 1, 2, 3, 4, 5))
  pmis_base_values <- parse_num_list(argv$pmis_base_values, c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1))
  p_wgd_values <- parse_num_list(
    argv$p_wgd_values,
    as.numeric(run_params$p_wgd) * c(1e-2, 1e-1, 1, 10, 100)
  )

  if (as_bool(argv$smoke, FALSE)) {
    initial_burden_values <- head(initial_burden_values, 1)
    trigger_burden_values <- head(trigger_burden_values, 1)
    o2_values <- head(o2_values, 1)
    pmis_base_values <- head(pmis_base_values, 1)
    p_wgd_values <- head(p_wgd_values, 1)
  }

  horizon_day <- as_num(argv$horizon_day, 1000)
  report_dt <- as_num(argv$report_dt, 1.0)
  trigger_check_dt <- as_num(argv$trigger_check_dt, report_dt)
  make_plots <- as_bool(argv$make_plots, TRUE)

  init_mean <- as_num(argv$initial_ploidy_mean, 3.0)
  init_sd <- as_num(argv$initial_ploidy_sd, 0.4)
  init_min <- as_num(argv$initial_ploidy_min, 1.5)
  init_max <- as_num(argv$initial_ploidy_max, 6.0)

  design <- build_design_rows(
    run_params = run_params,
    initial_burden_values = initial_burden_values,
    o2_values = o2_values,
    p_wgd_values = p_wgd_values,
    pmis_base_values = pmis_base_values,
    trigger_burden_values = trigger_burden_values,
    initial_burden_cells = initial_burden_cells
  )
  design$initial_ploidy_mean <- init_mean
  design$initial_ploidy_sd <- init_sd
  design$initial_ploidy_min <- init_min
  design$initial_ploidy_max <- init_max
  design$fit_dir <- fit_dir

  message("Running ", nrow(design), " scenarios. Output: ", out_dir)
  burden_rows <- vector("list", nrow(design))
  ploidy_rows <- vector("list", nrow(design))
  design_rows <- vector("list", nrow(design))

  for (i in seq_len(nrow(design))) {
    dr <- design[i, , drop = FALSE]
    message("[", i, "/", nrow(design), "] ", dr$scenario_id)
    init_state <- make_continuous_ploidy_state(
      cfg = cfg,
      total_live_cells = as.numeric(dr$initial_burden_cells),
      mean_ploidy = init_mean,
      sd_ploidy = init_sd,
      min_ploidy = init_min,
      max_ploidy = init_max
    )

    if (identical(as.character(dr$experiment), "pmiss_triggered_treatment")) {
      res <- run_triggered_pmiss_scenario(
        run_params = run_params,
        cfg = cfg,
        init_state = init_state,
        horizon_day = horizon_day,
        report_dt = report_dt,
        trigger_check_dt = trigger_check_dt,
        design_row = dr
      )
      design_rows[[i]] <- res$design
      burden_rows[[i]] <- res$burden
      ploidy_rows[[i]] <- res$ploidy
    } else {
      res <- run_static_scenario(
        run_params = run_params,
        cfg = cfg,
        init_state = init_state,
        horizon_day = horizon_day,
        report_dt = report_dt,
        design_row = dr
      )
      design_rows[[i]] <- res$burden[1, names(dr), drop = FALSE]
      burden_rows[[i]] <- res$burden
      ploidy_rows[[i]] <- res$ploidy
    }
  }

  design_out <- bind_rows(design_rows)
  burden_all <- bind_rows(burden_rows)
  ploidy_all <- bind_rows(ploidy_rows)
  ploidy_summary <- summarise_ploidy_timecourse(ploidy_all)

  endpoint_burden <- burden_all %>%
    group_by(scenario_id, experiment) %>%
    arrange(day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    rename(
      endpoint_p_mis_base = p_mis_base,
      endpoint_p_wgd = p_wgd,
      endpoint_o2_S0 = o2_S0
    ) %>%
    select(
      scenario_id, experiment, day, segment, local_day, step,
      endpoint_p_mis_base, endpoint_p_wgd, endpoint_o2_S0,
      starts_with("pred_burden"),
      starts_with("pred_log10"),
      starts_with("pred_o2")
    )

  endpoint_ploidy <- ploidy_summary %>%
    group_by(scenario_id, experiment) %>%
    arrange(day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup()

  endpoint_summary <- design_out %>%
    left_join(endpoint_burden, by = c("scenario_id", "experiment")) %>%
    left_join(endpoint_ploidy, by = c("scenario_id", "experiment", "day"))

  write_tsv(design_out, file.path(out_dir, "simulation_design.tsv"))
  write_tsv(burden_all, file.path(out_dir, "burden_timecourse.tsv"))
  write_tsv(ploidy_summary, file.path(out_dir, "ploidy_summary_timecourse.tsv"))
  write_tsv(endpoint_summary, file.path(out_dir, "endpoint_summary.tsv"))
  write_tsv_gz(ploidy_all, file.path(out_dir, "ploidy_distribution_timecourse.tsv.gz"))

  if (make_plots) {
    plot_outputs(burden_all = burden_all, ploidy_summary = ploidy_summary, design = design_out, out_dir = out_dir)
  }

  message("Done. Wrote outputs to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(out_dir)
}

if (sys.nframe() == 0) {
  main()
}
