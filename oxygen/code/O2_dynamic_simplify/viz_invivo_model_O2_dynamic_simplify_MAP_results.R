#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

# -----------------------------------------------------------------------------
# Function: get_script_dir_self
# Purpose: Resolve script directory path for robust relative file loading.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
get_script_dir_self <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", farg[[1]])))
}

# -----------------------------------------------------------------------------
# Function: as_num_vec
# Purpose: Convert an input value to the target scalar/vector type with safe defaults.
# Parameters:
#   - x: Input value or vector to process.
#   - default: Fallback value used when the input is NULL or invalid.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
as_num_vec <- function(x, default = numeric(0)) {
  if (is.null(x)) return(as.numeric(default))
  s <- trimws(as.character(x))
  if (!nzchar(s)) return(as.numeric(default))
  parts <- trimws(unlist(strsplit(s, "[,;]", perl = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(as.numeric(default))
  vals <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(vals))) stop("Invalid numeric vector argument: ", x)
  vals
}

# -----------------------------------------------------------------------------
# Function: find_latest_fit_dir
# Purpose: Discover fit result directories from a root path.
# Parameters:
#   - results_root: Root directory that contains multiple fit result folders.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
find_latest_fit_dir <- function(results_root) {
  dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  if (length(dirs) == 0) {
    stop("No fit result directories found under: ", results_root)
  }
  dirs[[which.max(file.info(dirs)$mtime)]]
}

# -----------------------------------------------------------------------------
# Function: normalize_cfg_for_viz
# Purpose: Normalize stored fit configuration for visualization-time simulation calls.
# Parameters:
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
normalize_cfg_for_viz <- function(cfg) {
  cfg$N_UNIT <- as.integer(cfg$N_UNIT %||% 22L)
  cfg$N_MIN <- as.integer(cfg$N_MIN %||% 22L)
  cfg$N_MAX <- as.integer(cfg$N_MAX %||% 154L)
  cfg$DT <- as.numeric(cfg$DT %||% 0.5)
  cfg$o2_cap_pct <- as.numeric(cfg$o2_cap_pct %||% 5.0)
  if (!is.finite(cfg$o2_cap_pct) || cfg$o2_cap_pct <= 0 || cfg$o2_cap_pct > 100) {
    stop("fit_config o2_cap_pct must be in (0, 100].")
  }
  cfg$o2_curve_type <- tolower(as.character(cfg$o2_curve_type %||% "gompertz"))
  if (!cfg$o2_curve_type %in% c("gompertz", "glogistic")) {
    stop("fit_config o2_curve_type must be gompertz or glogistic.")
  }
  cfg$o2_logN_eps <- as.numeric(cfg$o2_logN_eps %||% 1.0)
  if (!is.finite(cfg$o2_logN_eps) || cfg$o2_logN_eps <= 0) cfg$o2_logN_eps <- 1.0
  cfg$o2_anchor_N <- as.numeric(cfg$o2_anchor_N %||% 1e6)
  if (!is.finite(cfg$o2_anchor_N) || cfg$o2_anchor_N < 0) cfg$o2_anchor_N <- 1e6
  cfg$tau_O2_init <- as.numeric(cfg$tau_O2_init %||% 2.0)
  cfg$tau_O2 <- as.numeric(cfg$tau_O2 %||% cfg$tau_O2_init)
  if (!is.finite(cfg$tau_O2) || cfg$tau_O2 <= 0) cfg$tau_O2 <- cfg$tau_O2_init
  cfg$K <- as.numeric(cfg$K %||% 1e12)
  cfg$crowding <- as.character(cfg$crowding %||% "logistic")
  cfg$init_total_size <- as.numeric(cfg$init_total_size %||% 1e6)
  cfg$dose_ref <- as.numeric(cfg$dose_ref %||% 30)
  cfg$tx_mult_min <- as.numeric(cfg$tx_mult_min %||% 0.05)
  cfg$min_pop <- as.numeric(cfg$min_pop %||% 1e-12)
  cfg$dose_zero_only <- isTRUE(cfg$dose_zero_only %||% TRUE)
  cfg$fit_treatment <- isTRUE(cfg$fit_treatment %||% FALSE)
  cfg$max_scenarios <- as.numeric(cfg$max_scenarios %||% Inf)
  if (is.null(cfg$truncate_at_treatment)) {
    cfg$truncate_at_treatment <- isTRUE(cfg$pretreat_only %||% FALSE)
  }
  cfg$truncate_at_treatment <- isTRUE(cfg$truncate_at_treatment)

  if (is.null(cfg$ploidy_at_harvest)) {
    cfg$ploidy_at_harvest <- TRUE
  }
  cfg$ploidy_at_harvest <- isTRUE(cfg$ploidy_at_harvest)
  cfg
}

# -----------------------------------------------------------------------------
# Function: read_run_params
# Purpose: Read fitted parameter table and reconstruct run_params list.
# Parameters:
#   - fit_dir: Directory containing fitted parameters and summary outputs.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
read_run_params <- function(fit_dir, cfg = NULL) {
  p <- file.path(fit_dir, "best_params.tsv")
  if (!file.exists(p)) stop("Missing best_params.tsv in: ", fit_dir)
  tab <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain columns: parameter, value")
  }
  vals <- setNames(as.numeric(tab$value), as.character(tab$parameter))
  needed <- c(
    "lam_min", "lam_max", "k_o", "p_misseg", "k_o_mis",
    "beta_buffer", "n_exp", "smax", "p_wgd",
    "o2_init_pct", "o2_rate", "o2_shape_v"
  )
  miss <- setdiff(needed, names(vals))
  if (length(miss) > 0) {
    stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
  }
  out <- as.list(vals[needed])
  out$o2_curve_type <- as.character(.first_non_null_local(cfg$o2_curve_type, "gompertz"))
  out$o2_cap <- as.numeric(.first_non_null_local(cfg$o2_cap_pct, 5.0))
  out$o2_anchor_N <- as.numeric(.first_non_null_local(cfg$o2_anchor_N, 1e6))
  # Observation-layer parameters for burden (new: rho_2N; legacy: c_scale).
  if ("rho_2N" %in% names(vals) && is.finite(vals[["rho_2N"]]) && vals[["rho_2N"]] > 0) {
    out$rho_2N <- vals[["rho_2N"]]
  }
  if ("c_scale" %in% names(vals) && is.finite(vals[["c_scale"]]) && vals[["c_scale"]] > 0) {
    out$c_scale <- vals[["c_scale"]]
  }
  out$beta_size <- if ("beta_size" %in% names(vals) && is.finite(vals[["beta_size"]])) vals[["beta_size"]] else default_beta_size_prior_center()
  if ("c_vol_2N_mm3" %in% names(vals) && is.finite(vals[["c_vol_2N_mm3"]]) && vals[["c_vol_2N_mm3"]] > 0) {
    out$c_vol_2N_mm3 <- vals[["c_vol_2N_mm3"]]
  }
  out$alpha <- if ("alpha" %in% names(vals) && is.finite(vals[["alpha"]])) vals[["alpha"]] else 0
  out$gamma <- if ("gamma" %in% names(vals) && is.finite(vals[["gamma"]])) vals[["gamma"]] else 1
  out$tau_O2 <- if ("tau_O2" %in% names(vals) && is.finite(vals[["tau_O2"]]) && vals[["tau_O2"]] > 0) vals[["tau_O2"]] else as.numeric(.first_non_null_local(cfg$tau_O2, cfg$tau_O2_init, 2.0))
  out
}

# -----------------------------------------------------------------------------
# Function: compute_o2_eff_from_burden
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - Ntot: Total predicted cell count (or burden proxy) at current time.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
compute_o2_eff_from_burden <- function(Ntot, run_params, cfg) {
  o2_feedback <- isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE))
  if (!isTRUE(o2_feedback)) {
    return(as.numeric(clip(as.numeric(.first_non_null_local(run_params$o2_cap, cfg$o2_cap_pct, 5.0)), 0, 100)))
  }

  o2_logN_eps <- as.numeric(.first_non_null_local(cfg$o2_logN_eps, 1.0))
  if (!is.finite(o2_logN_eps) || o2_logN_eps <= 0) o2_logN_eps <- 1.0
  o2_anchor_N <- as.numeric(.first_non_null_local(run_params$o2_anchor_N, cfg$o2_anchor_N, 1e6))
  if (!is.finite(o2_anchor_N) || o2_anchor_N < 0) o2_anchor_N <- 1e6

  if (exists(".o2_sigmoid_supply_from_burden", mode = "function", inherits = TRUE)) {
    return(as.numeric(.o2_sigmoid_supply_from_burden(
      Ntot = Ntot,
      run_params = run_params,
      O2_cap = as.numeric(.first_non_null_local(run_params$o2_cap, cfg$o2_cap_pct, 5.0)),
      o2_logN_eps = o2_logN_eps,
      o2_anchor_N = o2_anchor_N
    )))
  }

  # Fallback: explicit sigmoid formulas if helper is unavailable.
  curve_type <- as.character(.first_non_null_local(run_params$o2_curve_type, cfg$o2_curve_type, "gompertz"))
  O2_cap_use <- clip(as.numeric(.first_non_null_local(run_params$o2_cap, cfg$o2_cap_pct, 5.0)), 1e-6, 100)
  o2_init <- clip(as.numeric(.first_non_null_local(run_params$o2_init_pct, 0.5)), 1e-6, O2_cap_use - 1e-6)
  o2_rate <- as.numeric(.first_non_null_local(run_params$o2_rate, 1.0))
  if (!is.finite(o2_rate) || o2_rate <= 0) o2_rate <- 1.0
  o2_shape_v <- as.numeric(.first_non_null_local(run_params$o2_shape_v, 1.0))
  if (!is.finite(o2_shape_v) || o2_shape_v <= 0) o2_shape_v <- 1.0
  x <- log10(pmax(Ntot, 0) + o2_logN_eps) - log10(o2_anchor_N + o2_logN_eps)

  y <- if (identical(curve_type, "glogistic")) {
    b <- -log((O2_cap_use / o2_init)^o2_shape_v - 1)
    O2_cap_use / (1 + exp(-(o2_rate * x + b)))^(1 / o2_shape_v)
  } else {
    b <- -log(-log(o2_init / O2_cap_use))
    O2_cap_use * exp(-exp(-(o2_rate * x + b)))
  }
  as.numeric(clip(y, 0, 100))
}

# -----------------------------------------------------------------------------
# Function: simulate_one_full
# Purpose: Run one forward simulation trajectory for a single scenario.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenario: Single scenario data object.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - report_dt: Sampling interval for reported trajectories.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
simulate_one_full <- function(run_params, scenario, cfg, report_dt = 1.0) {
  model_core <- build_model_core(cfg = cfg)
  grid_pre <- model_core$grid_pre
  grid_post <- model_core$grid_post
  R0 <- model_core$R0
  R1 <- model_core$R1
  init_state <- if (scenario$cohort == "2N") model_core$init_state_2N else model_core$init_state_4N
  keep_days <- sort(unique(c(
    0,
    as.numeric(scenario$sim_end_day),
    as.numeric(scenario$obs_days),
    seq(0, as.numeric(scenario$sim_end_day), by = report_dt)
  )))
  keep_days <- keep_days[is.finite(keep_days) & keep_days >= 0 & keep_days <= as.numeric(scenario$sim_end_day)]

  o2_base <- as.numeric(.first_non_null_local(run_params$o2_cap, cfg$o2_cap_pct, 5.0))
  sim <- run_in_vivo_crowd(
    run_params = run_params,
    O2_schedule = list(c(t0 = 0, t1 = Inf, O2 = o2_base)),
    T_end = as.numeric(scenario$sim_end_day),
    sample_days = keep_days,
    N_UNIT = cfg$N_UNIT,
    DT = cfg$DT,
    K = cfg$K,
    crowding = cfg$crowding,
    grid_pre = grid_pre,
    grid_post = grid_post,
    init_state = init_state,
    chr_lengths_bp = cfg$chr_lengths_bp
  )

  d <- sim$all_dists
  if (is.null(d) || nrow(d) == 0) {
    return(list(burden = data.frame(), ploidy = data.frame()))
  }

  ploidy_rows <- d %>%
    group_by(day, N) %>%
    summarise(
      fraction = sum(fraction, na.rm = TRUE),
      pop = max(pop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      harvest = scenario$harvest,
      cohort = scenario$cohort,
      dose = scenario$dose
    ) %>%
    select(harvest, cohort, dose, day, N, fraction, pop)

  vol_lut <- setNames(as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg)), as.character(grid_pre))
  burden_by_day <- ploidy_rows %>%
    group_by(day) %>%
    summarise(
      pred_burden_cells = max(pop, na.rm = TRUE),
      pred_burden_volume_mm3 = max(pop, na.rm = TRUE) * sum(fraction * vol_lut[as.character(N)], na.rm = TRUE),
      .groups = "drop"
    )

  burden_by_day$step <- as.integer(round(burden_by_day$day / cfg$DT))
  obs_steps <- as.integer(round(as.numeric(scenario$obs_days) / cfg$DT))
  obs_map <- setNames(as.numeric(scenario$obs_burden), as.character(obs_steps))
  burden_by_day$obs_burden <- as.numeric(obs_map[as.character(burden_by_day$step)])
  burden_by_day$pred_o2_pct <- vapply(
    burden_by_day$pred_burden_cells,
    function(x) as.numeric(compute_o2_eff_from_burden(Ntot = x, run_params = run_params, cfg = cfg)),
    numeric(1)
  )

  burden_rows <- burden_by_day %>%
    mutate(
      harvest = scenario$harvest,
      cohort = scenario$cohort,
      dose = scenario$dose,
      treat_day = scenario$treat_day,
      pred_burden = pred_burden_volume_mm3
    ) %>%
    select(
      harvest, cohort, dose, treat_day, step, day,
      pred_burden, pred_burden_volume_mm3, pred_burden_cells,
      pred_o2_pct, obs_burden
    )

  list(
    burden = as.data.frame(burden_rows),
    ploidy = as.data.frame(select(ploidy_rows, harvest, cohort, dose, day, N, fraction))
  )
}

# -----------------------------------------------------------------------------
# Function: simulate_one_full_horizon
# Purpose: Run one forward simulation trajectory for a single scenario.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenario: Single scenario data object.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - horizon_day: Prediction horizon end time in days.
#   - report_dt: Sampling interval for reported trajectories.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
simulate_one_full_horizon <- function(run_params, scenario, cfg, horizon_day, report_dt = 1.0) {
  sc <- scenario
  sc$sim_end_day <- as.numeric(max(horizon_day, 0))
  simulate_one_full(run_params, sc, cfg, report_dt = report_dt)
}

# -----------------------------------------------------------------------------
# Function: normalize_burden_for_plot
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - burden_all: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
normalize_burden_for_plot <- function(burden_all) {
  burden_all %>%
    group_by(harvest, cohort, dose) %>%
    arrange(day, .by_group = TRUE) %>%
    group_modify(function(df, .y) {
      pred_delta <- df$pred_burden - df$pred_burden[[1]]
      pred_scale <- max(abs(pred_delta), na.rm = TRUE)
      df$pred_norm <- if (is.finite(pred_scale) && pred_scale > 0) pred_delta / pred_scale else pred_delta

      obs_vals <- df$obs_burden[is.finite(df$obs_burden)]
      if (length(obs_vals) > 0) {
        obs_delta <- df$obs_burden - obs_vals[[1]]
        obs_scale <- max(abs(obs_delta), na.rm = TRUE)
        df$obs_norm <- if (is.finite(obs_scale) && obs_scale > 0) obs_delta / obs_scale else obs_delta
      } else {
        df$obs_norm <- NA_real_
      }
      df
    }) %>%
    ungroup()
}

# -----------------------------------------------------------------------------
# Function: compute_ploidy_weighted_mean
# Purpose: Compute weighted mean with finite/positive-weight safeguards.
# Parameters:
#   - ploidy_all: Function-specific input argument.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
compute_ploidy_weighted_mean <- function(ploidy_all, cfg) {
  ploidy_all %>%
    group_by(harvest, cohort, dose, day) %>%
    summarise(
      weighted_mean_N = sum(N * fraction, na.rm = TRUE) / pmax(sum(fraction, na.rm = TRUE), 1e-12),
      .groups = "drop"
    ) %>%
    mutate(weighted_mean_ploidy = weighted_mean_N / cfg$N_UNIT)
}

# -----------------------------------------------------------------------------
# Function: plot_functional_response_curves
# Purpose: Generate and save visualization output for fitted model behavior.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - out_dir: Output directory for generated files and plots.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
plot_functional_response_curves <- function(run_params, cfg, out_dir) {
  o2_grid <- seq(0, 100, by = 0.2)
  N_ref <- as.numeric(cfg$N_UNIT * 2)
  ms_rate <- if (exists(".pmisseg_of_O2", mode = "function", inherits = TRUE)) {
    as.numeric(.pmisseg_of_O2(o2_grid, run_params))
  } else {
    k_o_mis_use <- max(as.numeric(run_params$k_o_mis), 1e-12)
    as.numeric(run_params$p_misseg) * (1 - o2_grid / (o2_grid + k_o_mis_use))
  }
  prolif_rate <- as.numeric(growth_lambda(
    O2 = o2_grid,
    N = N_ref,
    lam_min = run_params$lam_min,
    lam_max = run_params$lam_max,
    k_o = run_params$k_o
  ))

  o2_curve <- data.frame(
    oxygen_pct = o2_grid,
    ms_rate = pmax(ms_rate, 0),
    proliferation_rate = pmax(prolif_rate, 0),
    row.names = NULL
  )
  write.table(
    o2_curve,
    file = file.path(out_dir, "functional_curve_oxygen.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  p_ms <- ggplot(o2_curve, aes(x = oxygen_pct, y = ms_rate)) +
    geom_line(color = "#d62728", linewidth = 1) +
    labs(
      title = "Oxygen vs Missegregation (MS) Rate",
      x = "Oxygen (%)",
      y = "MS rate"
    ) +
    theme_bw(base_size = 11)

  p_prolif <- ggplot(o2_curve, aes(x = oxygen_pct, y = proliferation_rate)) +
    geom_line(color = "#1f77b4", linewidth = 1) +
    labs(
      title = "Oxygen vs Proliferation Rate",
      subtitle = "From fitted growth_lambda functional form",
      x = "Oxygen (%)",
      y = "Proliferation rate"
    ) +
    theme_bw(base_size = 11)

  ploidy_grid <- seq(cfg$N_MIN / cfg$N_UNIT, cfg$N_MAX / cfg$N_UNIT, by = 0.02)
  N_grid <- pmax(ploidy_grid * cfg$N_UNIT, 1e-8)
  viability <- as.numeric(run_params$smax) * exp(
    -as.numeric(run_params$beta_buffer) * ((2 * cfg$N_UNIT) / N_grid)^as.numeric(run_params$n_exp)
  )
  viability_curve <- data.frame(
    ploidy = ploidy_grid,
    viability_after_ms = pmax(viability, 0),
    row.names = NULL
  )
  write.table(
    viability_curve,
    file = file.path(out_dir, "functional_curve_ploidy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  p_viability <- ggplot(viability_curve, aes(x = ploidy, y = viability_after_ms)) +
    geom_line(color = "#2ca02c", linewidth = 1) +
    labs(
      title = "Ploidy vs Viability After MS",
      subtitle = "Viability term from fitted buffering functional form",
      x = "Ploidy (N / N_UNIT)",
      y = "Viability after MS"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "oxygen_vs_ms_rate.pdf"), p_ms, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_proliferation_rate.pdf"), p_prolif, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_viability_after_ms.pdf"), p_viability, width = 10, height = 7)
}

# -----------------------------------------------------------------------------
# Function: plot_predict_horizon
# Purpose: Generate and save visualization output for fitted model behavior.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenarios: List of scenario-specific observation data and metadata.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - out_dir: Output directory for generated files and plots.
#   - horizon_day: Prediction horizon end time in days.
#   - report_dt: Sampling interval for reported trajectories.
#   - top_n: Maximum number of scenarios selected for detailed plotting.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
plot_predict_horizon <- function(run_params, scenarios, cfg, out_dir, horizon_day, report_dt = 1.0, top_n = 6L) {
  sim_list <- lapply(scenarios, function(sc) {
    simulate_one_full_horizon(run_params, sc, cfg, horizon_day = horizon_day, report_dt = report_dt)
  })
  burden_all <- bind_rows(lapply(sim_list, `[[`, "burden"))
  ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy"))
  if (nrow(burden_all) == 0 || nrow(ploidy_all) == 0) return(invisible(NULL))

  burden_all <- burden_all %>% filter(day <= horizon_day + 1e-9)
  ploidy_all <- ploidy_all %>% filter(day <= horizon_day + 1e-9)
  burden_all <- normalize_burden_for_plot(burden_all)
  ploidy_mean <- compute_ploidy_weighted_mean(ploidy_all, cfg)

  horizon_tag <- paste0("0_", as.integer(round(horizon_day)), "day")
  # Remove deprecated multi-file prediction plot outputs for this horizon to avoid stale files.
  unlink(file.path(out_dir, c(
    paste0("predict_burden_normalized_", horizon_tag, ".pdf"),
    paste0("predict_burden_absolute_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_heatmap_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_top_states_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_weighted_mean_", horizon_tag, ".pdf"),
    paste0("forecast_burden_normalized_", horizon_tag, ".pdf"),
    paste0("forecast_burden_absolute_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_heatmap_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_top_states_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_weighted_mean_", horizon_tag, ".pdf")
  )), force = TRUE)

  write.table(burden_all, file = file.path(out_dir, paste0("predict_burden_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_all, file = file.path(out_dir, paste0("predict_ploidy_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_mean, file = file.path(out_dir, paste0("predict_ploidy_weighted_mean_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  burden_plot_df <- burden_all %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      metric = "Burden (normalized)",
      value = as.numeric(pred_norm)
    ) %>%
    bind_rows(
      burden_all %>%
        transmute(
          harvest = as.character(harvest),
          cohort = as.character(cohort),
          dose = as.numeric(dose),
          day = as.numeric(day),
          metric = "Burden (absolute)",
          value = as.numeric(pred_burden)
        )
    )

  ploidy_plot_df <- ploidy_mean %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      metric = "Weighted mean ploidy",
      value = as.numeric(weighted_mean_ploidy)
    )

  predict_plot_df <- bind_rows(burden_plot_df, ploidy_plot_df) %>%
    mutate(
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__"),
      metric = factor(metric, levels = c("Burden (normalized)", "Burden (absolute)", "Weighted mean ploidy"))
    )

  p_predict <- ggplot(
    predict_plot_df,
    aes(x = day, y = value, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    coord_cartesian(xlim = c(0, horizon_day)) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = paste0("Predict Curves: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = paste0("Single summary plot (all scenarios overlaid) | fit_dir=", basename(dirname(out_dir)), " | report_dt=", report_dt),
      x = "Day",
      y = NULL,
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, paste0("predict_curves_", horizon_tag, ".pdf")), p_predict, width = 12, height = 11)

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Function: find_fit_dirs_under
# Purpose: Discover fit result directories from a root path.
# Parameters:
#   - root_dir: Root directory used for recursive fit-folder discovery.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
find_fit_dirs_under <- function(root_dir) {
  all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  sub_dirs <- all_dirs[all_dirs != root_dir]

# -----------------------------------------------------------------------------
# Function: is_fit_dir
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - d: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  is_fit_dir <- function(d) {
    file.exists(file.path(d, "fit_config.rds")) &&
      file.exists(file.path(d, "best_params.tsv"))
  }

  # Primary mode: traverse subdirectories under fit_dir.
  fit_sub_dirs <- sub_dirs[vapply(sub_dirs, is_fit_dir, logical(1))]
  if (length(fit_sub_dirs) > 0) {
    return(sort(unique(fit_sub_dirs)))
  }

  # Fallback mode: if no fit subdirectories exist, use fit_dir itself.
  if (is_fit_dir(root_dir)) {
    return(normalizePath(root_dir, mustWork = TRUE))
  }

  character(0)
}

# -----------------------------------------------------------------------------
# Function: run_viz_for_fit_dir
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - fit_dir: Directory containing fitted parameters and summary outputs.
#   - argv: Character vector of command-line arguments in --key=value format.
#   - dt_path: Path to burden observation table (Excel file).
#   - ploidy_path: Path to terminal ploidy table (TSV file).
#   - report_dt: Sampling interval for reported trajectories.
#   - top_n: Maximum number of scenarios selected for detailed plotting.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
run_viz_for_fit_dir <- function(
  fit_dir,
  argv,
  dt_path,
  ploidy_path,
  report_dt,
  top_n
) {
  out_dir <- file.path(fit_dir, "viz")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cfg_path <- file.path(fit_dir, "fit_config.rds")
  if (!file.exists(cfg_path)) stop("Missing fit_config.rds in: ", fit_dir)
  cfg <- readRDS(cfg_path)
  cfg <- normalize_cfg_for_viz(cfg)

  if (!is.null(argv$dose_zero_only)) cfg$dose_zero_only <- as_bool(argv$dose_zero_only, cfg$dose_zero_only)
  if (!is.null(argv$truncate_at_treatment)) cfg$truncate_at_treatment <- as_bool(argv$truncate_at_treatment, cfg$truncate_at_treatment)
  if (!is.null(argv$ploidy_at_harvest)) cfg$ploidy_at_harvest <- as_bool(argv$ploidy_at_harvest, cfg$ploidy_at_harvest)
  if (!is.null(argv$max_scenarios)) cfg$max_scenarios <- as_num(argv$max_scenarios, cfg$max_scenarios)

  run_params <- read_run_params(fit_dir, cfg = cfg)
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)

  sim_list <- lapply(scenarios, function(sc) simulate_one_full(run_params, sc, cfg, report_dt = report_dt))
  burden_all <- bind_rows(lapply(sim_list, `[[`, "burden"))
  ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy"))

  if (nrow(burden_all) == 0 || nrow(ploidy_all) == 0) {
    stop("No simulation output generated; check fit/data configuration.")
  }

  burden_all <- normalize_burden_for_plot(burden_all)

  write.table(burden_all, file = file.path(out_dir, "burden_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_all, file = file.path(out_dir, "ploidy_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  ploidy_mean <- compute_ploidy_weighted_mean(ploidy_all, cfg)
  write.table(ploidy_mean, file = file.path(out_dir, "ploidy_weighted_mean_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  p_burden <- ggplot(burden_all, aes(x = day, y = pred_norm)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_line(
      data = burden_all %>% filter(!is.na(obs_norm)),
      aes(y = obs_norm),
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = burden_all %>% filter(!is.na(obs_norm)),
      aes(y = obs_norm),
      color = "black",
      size = 1
    ) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = "O2 Dynamic Simplify MAP Model: In Vivo Burden Trajectory (Normalized)",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Day",
      y = "Normalized Burden (delta / max|delta|)"
    ) +
    theme_bw(base_size = 11)

  p_burden_abs <- ggplot(burden_all, aes(x = day, y = pred_burden)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_line(
      data = burden_all %>% filter(!is.na(obs_burden)),
      aes(y = obs_burden),
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = burden_all %>% filter(!is.na(obs_burden)),
      aes(y = obs_burden),
      color = "black",
      size = 1
    ) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_y") +
    labs(
      title = "O2 Dynamic Simplify MAP Model: In Vivo Burden Trajectory (Absolute)",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Day",
      y = "Tumor burden (mm^3)"
    ) +
    theme_bw(base_size = 11)

  rho_2N_min <- as_num(argv$rho_2N_min, 3.2e4)
  rho_2N_max <- as_num(argv$rho_2N_max, 5.6e4)
  if (!is.finite(rho_2N_min) || rho_2N_min <= 0) rho_2N_min <- 3.2e4
  if (!is.finite(rho_2N_max) || rho_2N_max <= 0) rho_2N_max <- 5.6e4
  if (rho_2N_min > rho_2N_max) {
    tmp <- rho_2N_min
    rho_2N_min <- rho_2N_max
    rho_2N_max <- tmp
  }
  rho_2N_mid <- sqrt(rho_2N_min * rho_2N_max)
  pred_cell_col <- if ("pred_burden_cells" %in% names(burden_all)) "pred_burden_cells" else "pred_burden"
  burden_all_real <- burden_all %>%
    mutate(
      pred_burden_cell_number = as.numeric(.data[[pred_cell_col]]),
      obs_burden_cell_number_low = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_min, NA_real_),
      obs_burden_cell_number_mid = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_mid, NA_real_),
      obs_burden_cell_number_high = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_max, NA_real_)
    )

  p_burden_abs_real <- ggplot(burden_all_real, aes(x = day, y = pred_burden_cell_number)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_ribbon(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_low) & !is.na(obs_burden_cell_number_high)),
      aes(x = day, ymin = obs_burden_cell_number_low, ymax = obs_burden_cell_number_high),
      inherit.aes = FALSE,
      fill = "grey50",
      alpha = 0.18
    ) +
    geom_line(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_mid)),
      aes(y = obs_burden_cell_number_mid),
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_mid)),
      aes(y = obs_burden_cell_number_mid),
      color = "black",
      size = 1
    ) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_y") +
    labs(
      title = "O2 Dynamic Simplify MAP Model: In Vivo Burden Trajectory (Absolute, Real Scale)",
      subtitle = paste0(
        "fit_dir=", basename(fit_dir),
        " | report_dt=", report_dt,
        " | obs burden -> CellNumber using rho_2N range=[",
        signif(rho_2N_min, 4), ", ", signif(rho_2N_max, 4), "] cells/mm^3",
        " (mid=", signif(rho_2N_mid, 4), ")"
      ),
      x = "Day",
      y = "CellNumber (2N-equivalent range)"
    ) +
    theme_bw(base_size = 11)

  o2_burden_df <- burden_all %>%
    filter(is.finite(pred_burden), is.finite(pred_o2_pct)) %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      burden_mm3 = as.numeric(pred_burden),
      o2_pct = as.numeric(pred_o2_pct),
      sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__")
    )
  write.table(o2_burden_df, file = file.path(out_dir, "predict_burden_vs_o2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  p_burden_vs_o2 <- ggplot(o2_burden_df, aes(x = burden_mm3, y = o2_pct, color = cohort, group = sample_id)) +
    geom_path(linewidth = 0.75, alpha = 0.9) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_x") +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = "O2 Dynamic Simplify MAP Model: Predicted Oxygen vs Burden",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Tumor burden (mm^3)",
      y = "Oxygen (%)",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_heatmap <- ggplot(ploidy_all, aes(x = day, y = N, fill = fraction)) +
    geom_raster(interpolate = FALSE) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(min(ploidy_all$N, na.rm = TRUE), 100)) +
    scale_fill_gradientn(
      colours = viridisLite::viridis(3, option = "C"),
      values = c(0, 0.1, 1),
      limits = c(0, 1)
    ) +
    labs(
      title = "O2 Dynamic Simplify MAP Model: Predicted Ploidy Distribution Over Time",
      subtitle = "Heatmap of fraction by chromosome number (N)",
      x = "Day",
      y = "Chromosome Number (N)",
      fill = "Fraction"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey84", linewidth = 0.06),
      panel.grid.minor = element_blank(),
      panel.ontop = TRUE,
      panel.background = element_rect(fill = NA, color = NA)
    )

  top_states <- ploidy_all %>%
    group_by(N) %>%
    summarise(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_fraction)) %>%
    slice_head(n = top_n) %>%
    pull(N)
  ploidy_top <- ploidy_all %>%
    filter(N %in% top_states) %>%
    mutate(N = factor(N, levels = top_states))

  p_ploidy_lines <- ggplot(ploidy_top, aes(x = day, y = fraction, color = N)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ harvest, ncol = 2) +
    labs(
      title = paste0("O2 Dynamic Simplify MAP Model: Ploidy Over Time (Top ", top_n, " N States)"),
      x = "Day",
      y = "Fraction",
      color = "N"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_weighted_mean <- ggplot(ploidy_mean, aes(x = day, y = weighted_mean_ploidy)) +
    geom_line(color = "#d62728", linewidth = 0.9) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(min(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE), 5)) +
    labs(
      title = "O2 Dynamic Simplify MAP Model: Weighted Mean Ploidy Over Time",
      subtitle = "Weighted by predicted ploidy fractions",
      x = "Day",
      y = "Weighted Mean Ploidy (P = N / N_UNIT)"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "burden_trend.pdf"), p_burden, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute.pdf"), p_burden_abs, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute(real_scale).pdf"), p_burden_abs_real, width = 13, height = 9)
  ggsave(file.path(out_dir, "predict_burden_vs_o2.pdf"), p_burden_vs_o2, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_heatmap_over_time.pdf"), p_ploidy_heatmap, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_top_states_over_time.pdf"), p_ploidy_lines, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_weighted_mean_over_time.pdf"), p_ploidy_weighted_mean, width = 13, height = 9)
  plot_functional_response_curves(run_params = run_params, cfg = cfg, out_dir = out_dir)

  predict_horizons <- as_num_vec(argv$predict_horizons, c(100, 300, 1000))
  predict_horizons <- sort(unique(predict_horizons[is.finite(predict_horizons) & predict_horizons > 0]))
  predict_report_dt <- as_num(argv$predict_report_dt, report_dt)
  if (!is.finite(predict_report_dt) || predict_report_dt <= 0) predict_report_dt <- report_dt
  predict_top_n <- as_int(argv$predict_top_n, top_n)
  if (!is.finite(predict_top_n) || predict_top_n < 1) predict_top_n <- top_n
  do_predict_plots <- as_bool(argv$predict_plots, TRUE)

  if (isTRUE(do_predict_plots) && length(predict_horizons) > 0) {
    for (hz in predict_horizons) {
      message("  Predict plots: 0-", hz, " days (report_dt=", predict_report_dt, ")")
      plot_predict_horizon(
        run_params = run_params,
        scenarios = scenarios,
        cfg = cfg,
        out_dir = out_dir,
        horizon_day = hz,
        report_dt = predict_report_dt,
        top_n = predict_top_n
      )
    }
  }

  normalizePath(out_dir)
}

# -----------------------------------------------------------------------------
# Function: main
# Purpose: Entry point: parse options, run fitting/visualization workflow, and write outputs.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
main <- function() {
  script_dir <- get_script_dir_self()
  source(file.path(script_dir, "fit_invivo_model_O2_dynamic_simplify_MAP.R"))
  source(file.path(script_dir, "model_O2_dynamic_simplify.R"))

  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  results_root <- normalizePath(file.path(script_dir, "..", "..", "results"), mustWork = FALSE)
  fit_root <- if (!is.null(argv$fit_dir)) {
    normalizePath(argv$fit_dir, mustWork = TRUE)
  } else {
    normalizePath(find_latest_fit_dir(results_root), mustWork = TRUE)
  }

  if (!is.null(argv$out_dir)) {
    message("Ignoring --out_dir. Outputs are always written to each subdirectory's /viz.")
  }

  report_dt <- as_num(argv$report_dt, 1.0)
  if (!is.finite(report_dt) || report_dt <= 0) stop("report_dt must be > 0")
  top_n <- as_int(argv$top_n, 6L)
  if (!is.finite(top_n) || top_n < 1) stop("top_n must be >= 1")

  data_dir <- if (!is.null(argv$data_dir)) {
    argv$data_dir
  } else {
    normalizePath(file.path(script_dir, "..", "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
  }
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- file.path(data_dir, "all_ploidy.tsv")

  fit_dirs <- find_fit_dirs_under(fit_root)
  if (length(fit_dirs) == 0) {
    stop(
      "No valid fit results found. Need either: ",
      "(a) fit subdirectories under fit_dir, or ",
      "(b) fit_config.rds + best_params.tsv directly under fit_dir. fit_dir=", fit_root
    )
  }

  fit_dirs <- sort(unique(fit_dirs))
  message("Found ", length(fit_dirs), " fit directories under: ", fit_root)

  ok <- character(0)
  failed <- character(0)

  for (i in seq_along(fit_dirs)) {
    fit_dir <- fit_dirs[[i]]
    message("[", i, "/", length(fit_dirs), "] Processing: ", fit_dir)
    tryCatch(
      {
        out_dir <- run_viz_for_fit_dir(
          fit_dir = fit_dir,
          argv = argv,
          dt_path = dt_path,
          ploidy_path = ploidy_path,
          report_dt = report_dt,
          top_n = top_n
        )
        ok <- c(ok, fit_dir)
        message("  Done: ", out_dir)
      },
      error = function(e) {
        failed <<- c(failed, paste0(fit_dir, " :: ", conditionMessage(e)))
        message("  Failed: ", conditionMessage(e))
      }
    )
  }

  if (length(ok) == 0) {
    stop("All fit subdirectories failed.")
  }

  message("Finished. Success: ", length(ok), " | Failed: ", length(failed))
  if (length(failed) > 0) {
    message("Failed directories:")
    for (x in failed) message("  - ", x)
  }
}

if (sys.nframe() == 0) {
  main()
}
