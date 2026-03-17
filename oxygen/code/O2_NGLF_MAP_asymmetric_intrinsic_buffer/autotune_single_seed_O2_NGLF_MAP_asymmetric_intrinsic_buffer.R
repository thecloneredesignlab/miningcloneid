#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

parse_args <- function(argv) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) >= 2L) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) default else v
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v) || !is.finite(v)) default else v
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  v <- trimws(toupper(as.character(x)[1]))
  if (v %in% c("TRUE", "T", "1", "YES", "Y", "ON")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO", "N", "OFF")) return(FALSE)
  default
}

fmt_num <- function(x, digits = 8) {
  if (!is.finite(x)) return("NA")
  format(signif(x, digits), scientific = TRUE, trim = TRUE)
}

yaml_value_string <- function(x) {
  if (is.logical(x)) return(if (isTRUE(x)) "TRUE" else "FALSE")
  if (is.numeric(x)) return(fmt_num(x, digits = 10))
  x <- as.character(x)
  if (!nzchar(x)) return("\"\"")
  x
}

read_yaml_lines <- function(path) readLines(path, warn = FALSE)

yaml_get <- function(lines, key) {
  pat <- paste0("^\\s*", key, "\\s*:\\s*(.*)$")
  hit <- grep(pat, lines)
  if (!length(hit)) return(NULL)
  val <- sub(pat, "\\1", lines[hit[1]])
  val <- trimws(sub("\\s+#.*$", "", val))
  if (identical(val, "\"\"") || identical(val, "''")) return("")
  if (startsWith(val, "\"") && endsWith(val, "\"")) val <- substring(val, 2, nchar(val) - 1L)
  if (startsWith(val, "'") && endsWith(val, "'")) val <- substring(val, 2, nchar(val) - 1L)
  val
}

yaml_set <- function(lines, key, value) {
  pat <- paste0("^\\s*", key, "\\s*:")
  repl <- paste0(key, ": ", yaml_value_string(value))
  hit <- grep(pat, lines)
  if (length(hit)) {
    lines[hit[1]] <- repl
  } else {
    lines <- c(lines, repl)
  }
  lines
}

resolve_path_local <- function(path_value, base_dir) {
  if (is.null(path_value)) return(NULL)
  if (!nzchar(path_value)) return("")
  if (grepl("^~", path_value)) return(path.expand(path_value))
  if (grepl("^/", path_value)) return(path_value)
  normalizePath(file.path(base_dir, path_value), mustWork = FALSE)
}

encode_transformed <- function(param_name, value_natural) {
  if (!is.finite(value_natural)) return(NA_real_)
  if (startsWith(param_name, "log10_")) return(log10(value_natural))
  if (startsWith(param_name, "logit_")) return(qlogis(min(max(value_natural, 1e-8), 1 - 1e-8)))
  value_natural
}

decode_transformed <- function(param_name, value_transformed) {
  if (!is.finite(value_transformed)) return(NA_real_)
  if (startsWith(param_name, "log10_")) return(10^value_transformed)
  if (startsWith(param_name, "logit_")) return(plogis(value_transformed))
  value_transformed
}

read_parameter_table <- function(path) {
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("param_name", "estimate", "init_value", "lower_bound", "upper_bound")
  miss <- setdiff(required, names(tab))
  if (length(miss)) stop("parameter table missing columns: ", paste(miss, collapse = ", "))
  tab
}

set_param_natural <- function(tab, param_name, init = NULL, lower = NULL, upper = NULL, note_append = NULL) {
  idx <- match(param_name, tab$param_name)
  if (is.na(idx)) stop("parameter not found in table: ", param_name)
  cur_lower <- decode_transformed(param_name, as.numeric(tab$lower_bound[idx]))
  cur_upper <- decode_transformed(param_name, as.numeric(tab$upper_bound[idx]))
  cur_init <- decode_transformed(param_name, as.numeric(tab$init_value[idx]))
  lower_use <- lower %||% cur_lower
  upper_use <- upper %||% cur_upper
  init_use <- init %||% cur_init
  init_use <- min(max(init_use, lower_use), upper_use)
  tab$lower_bound[idx] <- encode_transformed(param_name, lower_use)
  tab$upper_bound[idx] <- encode_transformed(param_name, upper_use)
  tab$init_value[idx] <- encode_transformed(param_name, init_use)
  if (!is.null(note_append) && "note" %in% names(tab)) {
    note_cur <- trimws(tab$note[idx] %||% "")
    if (!nzchar(note_cur)) {
      tab$note[idx] <- note_append
    } else if (!grepl(note_append, note_cur, fixed = TRUE)) {
      tab$note[idx] <- paste(note_cur, note_append, sep = " | ")
    }
  }
  tab
}

get_param_natural_bounds <- function(tab, param_name) {
  idx <- match(param_name, tab$param_name)
  if (is.na(idx)) stop("parameter not found in table: ", param_name)
  list(
    init = decode_transformed(param_name, as.numeric(tab$init_value[idx])),
    lower = decode_transformed(param_name, as.numeric(tab$lower_bound[idx])),
    upper = decode_transformed(param_name, as.numeric(tab$upper_bound[idx]))
  )
}

param_hit_upper <- function(value, upper, tol_frac = 0.05) {
  if (!is.finite(value) || !is.finite(upper)) return(FALSE)
  if (upper == 0) return(abs(value - upper) < 1e-12)
  abs(value - upper) <= tol_frac * max(abs(upper), 1e-12)
}

param_hit_lower <- function(value, lower, tol_frac = 0.05) {
  if (!is.finite(value) || !is.finite(lower)) return(FALSE)
  if (lower == 0) return(abs(value - lower) < 1e-12)
  abs(value - lower) <= tol_frac * max(abs(lower), 1e-12)
}

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

dist_to_interval <- function(x, lower = 2.35, upper = 2.65) {
  pmax(0, lower - x, x - upper)
}

set_prior_center <- function(cfg_lines, key, value) yaml_set(cfg_lines, key, value)

apply_recommended_start <- function(cfg_lines, tab,
                                    growth_penalty_ploidy_on = FALSE,
                                    growth_penalty_hypoxia_on = TRUE,
                                    death_on = FALSE) {
  cfg_next <- cfg_lines
  tab_next <- tab

  cfg_next <- yaml_set(cfg_next, "growth_penalty_hypoxia", growth_penalty_hypoxia_on)
  cfg_next <- yaml_set(cfg_next, "growth_penalty_ploidy", growth_penalty_ploidy_on)
  cfg_next <- yaml_set(cfg_next, "death", death_on)
  cfg_next <- yaml_set(cfg_next, "auto_viz", TRUE)
  cfg_next <- yaml_set(cfg_next, "viz_report_dt", 1)
  cfg_next <- yaml_set(cfg_next, "viz_top_n", 6)
  # Strengthen the inner search so the autotune loop is not dominated by a weak
  # local DE run along one narrow p_misseg/beta_loss ridge.
  cfg_next <- yaml_set(cfg_next, "itermax", 80)
  cfg_next <- yaml_set(cfg_next, "NP", 160)
  cfg_next <- yaml_set(cfg_next, "n_starts", 16)
  cfg_next <- yaml_set(cfg_next, "optim_maxit", 6000)
  cfg_next <- yaml_set(cfg_next, "de_steptol", 15)

  cfg_next <- yaml_set(cfg_next, "o2_init_pct_init", 2.3)
  cfg_next <- yaml_set(cfg_next, "o2_init_pct_min", 2.0)
  cfg_next <- yaml_set(cfg_next, "o2_init_pct_max", 2.6)
  cfg_next <- yaml_set(cfg_next, "o2_rate_init", 0.70)
  cfg_next <- yaml_set(cfg_next, "o2_rate_min", 0.50)
  cfg_next <- yaml_set(cfg_next, "o2_rate_max", 1.00)
  cfg_next <- yaml_set(cfg_next, "o2_shape_v_init", 1.02)
  cfg_next <- yaml_set(cfg_next, "o2_shape_v_min", 0.95)
  cfg_next <- yaml_set(cfg_next, "o2_shape_v_max", 1.15)
  cfg_next <- yaml_set(cfg_next, "k_clear_init", 1.8e-3)
  cfg_next <- yaml_set(cfg_next, "k_clear_min", 1.0e-3)
  cfg_next <- yaml_set(cfg_next, "k_clear_max", 3.0e-3)
  cfg_next <- yaml_set(cfg_next, "sigma_burden", 0.20)
  cfg_next <- yaml_set(cfg_next, "sigma_burden_min", 0.17)
  cfg_next <- yaml_set(cfg_next, "sigma_burden_max", 0.24)
  cfg_next <- yaml_set(cfg_next, "rho_2N_min", 38000)
  cfg_next <- yaml_set(cfg_next, "rho_2N_max", 52000)

  cfg_next <- set_prior_center(cfg_next, "prior_center_log10_o2_rate", log10(0.70))
  cfg_next <- set_prior_center(cfg_next, "prior_center_log10_o2_init_pct", log10(2.3))
  cfg_next <- set_prior_center(cfg_next, "prior_center_log10_o2_shape_v", 0.0)
  cfg_next <- set_prior_center(cfg_next, "prior_center_beta_size", 0.45)
  cfg_next <- set_prior_center(cfg_next, "prior_center_log10_beta_loss", log10(0.095))
  cfg_next <- set_prior_center(cfg_next, "prior_center_logit_gain_loss_ratio", qlogis(0.15))
  cfg_next <- set_prior_center(cfg_next, "prior_center_log10_rho_2N", 4.64)
  cfg_next <- set_prior_center(cfg_next, "prior_center_log10_k_clear", log10(1.8e-3))

  tab_next <- set_param_natural(tab_next, "log10_p_misseg", init = 0.055, lower = 0.045, upper = 0.065,
                                note_append = "autotune v1: keep misseg high enough to avoid the low-pressure ridge")
  tab_next <- set_param_natural(tab_next, "log10_k_o_mis", init = 0.35, lower = 0.15, upper = 0.60,
                                note_append = "autotune v1: emphasize later hypoxia-triggered misseg")
  tab_next <- set_param_natural(tab_next, "log10_beta_loss", init = 0.095, lower = 0.08, upper = 0.12,
                                note_append = "autotune v1: keep a moderate daughter loss floor with event surcharge")
  tab_next <- set_param_natural(tab_next, "logit_gain_loss_ratio", init = 0.15, lower = 0.10, upper = 0.25,
                                note_append = "autotune v1: stronger gain penalty to discourage high-copy gain daughters")
  tab_next <- set_param_natural(tab_next, "log10_p_wgd", init = 1e-5, lower = 1e-8, upper = 1e-4,
                                note_append = "autotune v1: keep WGD weak during untreated fitting")
  tab_next <- set_param_natural(tab_next, "log10_o2_init_pct", init = 2.3, lower = 2.0, upper = 2.6,
                                note_append = "autotune v1: bring hypoxia onset slightly earlier")
  tab_next <- set_param_natural(tab_next, "log10_o2_rate", init = 0.70, lower = 0.50, upper = 1.00,
                                note_append = "autotune v1: faster oxygen decline to create later-stage convergence pressure")
  tab_next <- set_param_natural(tab_next, "log10_o2_shape_v", init = 1.02, lower = 0.95, upper = 1.15,
                                note_append = "autotune v1: moderate glogistic shape")
  tab_next <- set_param_natural(tab_next, "log10_rho_2N", init = 44000, lower = 38000, upper = 52000,
                                note_append = "autotune v1: tighter global burden scale")
  tab_next <- set_param_natural(tab_next, "beta_size", init = 0.45, lower = 0.35, upper = 0.55,
                                note_append = "autotune v1: trim chronic 4N burden overweighting")
  tab_next <- set_param_natural(tab_next, "log10_k_clear", init = 1.8e-3, lower = 1.0e-3, upper = 3.0e-3,
                                note_append = "autotune v1: keep dead-mass contribution bounded")
  tab_next <- set_param_natural(tab_next, "log10_sigma_burden", init = 0.20, lower = 0.17, upper = 0.24,
                                note_append = "autotune v1: tighten burden slack")

  if (isTRUE(growth_penalty_hypoxia_on)) {
    cfg_next <- yaml_set(cfg_next, "alpha_o2_init", 1.1)
    cfg_next <- yaml_set(cfg_next, "alpha_o2_min", 0.7)
    cfg_next <- yaml_set(cfg_next, "alpha_o2_max", 1.8)
    tab_next <- set_param_natural(tab_next, "log10_alpha_o2", init = 1.1, lower = 0.7, upper = 1.8,
                                  note_append = "autotune v1: weak background hypoxia damping under fixed eta_hypoxia=0.5")
  }

  if (isTRUE(growth_penalty_ploidy_on)) {
    cfg_next <- yaml_set(cfg_next, "gamma_growth_init", 2.0)
    cfg_next <- yaml_set(cfg_next, "gamma_growth_min", 1.0)
    cfg_next <- yaml_set(cfg_next, "gamma_growth_max", 3.0)
    tab_next <- set_param_natural(tab_next, "gamma_growth", init = 2.0, lower = 1.0, upper = 3.0,
                                  note_append = "autotune v1: only active when ploidy-size growth penalty is enabled")
  }

  if (isTRUE(death_on)) {
    cfg_next <- yaml_set(cfg_next, "mu_hp_init", 4.0e-3)
    cfg_next <- yaml_set(cfg_next, "mu_hp_min", 5.0e-4)
    cfg_next <- yaml_set(cfg_next, "mu_hp_max", 3.0e-2)
    cfg_next <- set_prior_center(cfg_next, "prior_center_log10_mu_hp", log10(4.0e-3))
    tab_next <- set_param_natural(tab_next, "log10_mu_hp", init = 4.0e-3, lower = 5.0e-4, upper = 3.0e-2,
                                  note_append = "autotune v1: moderate hypoxia-linked continuous death")
  }

  list(cfg_lines = cfg_next, param_table = tab_next)
}

analyze_run <- function(run_dir, n_unit = 22) {
  seed_dir <- file.path(run_dir, "seed1")
  fit_summary <- read.delim(file.path(seed_dir, "fit_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  best_params <- read.delim(file.path(seed_dir, "best_params.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  burden_fit <- read.delim(file.path(seed_dir, "burden_fit.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  terminal_ploidy <- read.delim(file.path(seed_dir, "terminal_ploidy_fit.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

  summary_map <- setNames(fit_summary$value, fit_summary$metric)
  best_map <- setNames(as.numeric(best_params$value), best_params$parameter)

  by_harvest <- split(burden_fit, burden_fit$harvest)
  latest_rows <- lapply(by_harvest, function(df) df[which.max(as.numeric(df$day)), , drop = FALSE])
  latest_df <- do.call(rbind, latest_rows)
  latest_df_2n <- latest_df[latest_df$cohort == "2N", , drop = FALSE]
  latest_df_4n <- latest_df[latest_df$cohort == "4N", , drop = FALSE]

  frac_safe <- function(num, den) num / pmax(den, 1e-12)

  latest_live_frac <- with(latest_df, frac_safe(as.numeric(pred_pop_live), as.numeric(pred_pop)))
  latest_dead_h_frac <- with(latest_df, frac_safe(as.numeric(pred_pop_dead_hypoxia), as.numeric(pred_pop)))
  latest_dead_b_frac <- with(latest_df, frac_safe(as.numeric(pred_pop_dead_buffer), as.numeric(pred_pop)))
  latest_dead_total_frac <- with(latest_df, frac_safe(as.numeric(pred_pop_dead_total), as.numeric(pred_pop)))

  latest_vol_live_frac <- with(latest_df, frac_safe(as.numeric(pred_burden_live_volume_mm3), as.numeric(pred_burden_volume_mm3)))
  latest_vol_dead_total_frac <- with(latest_df, frac_safe(as.numeric(pred_burden_dead_total_volume_mm3), as.numeric(pred_burden_volume_mm3)))
  latest_vol_dead_buffer_frac <- with(latest_df, frac_safe(as.numeric(pred_burden_dead_buffer_volume_mm3), as.numeric(pred_burden_volume_mm3)))

  burden_eval <- burden_fit[as.numeric(burden_fit$day) > 0, , drop = FALSE]
  burden_eval_2n <- burden_eval[burden_eval$cohort == "2N", , drop = FALSE]
  burden_eval_4n <- burden_eval[burden_eval$cohort == "4N", , drop = FALSE]
  obs <- as.numeric(burden_eval$obs_burden)
  pred <- as.numeric(burden_eval$pred_burden_volume_mm3)
  log_err <- log(pmax(pred, 1e-12)) - log(pmax(obs, 1e-12))
  log_err_2n <- log(pmax(as.numeric(burden_eval_2n$pred_burden_volume_mm3), 1e-12)) -
    log(pmax(as.numeric(burden_eval_2n$obs_burden), 1e-12))
  log_err_4n <- log(pmax(as.numeric(burden_eval_4n$pred_burden_volume_mm3), 1e-12)) -
    log(pmax(as.numeric(burden_eval_4n$obs_burden), 1e-12))
  raw_r2 <- {
    ss_tot <- sum((obs - mean(obs))^2)
    if (!is.finite(ss_tot) || ss_tot <= 0) NA_real_ else 1 - sum((pred - obs)^2) / ss_tot
  }

  ploidy_by_harvest <- split(terminal_ploidy, terminal_ploidy$harvest)
  ploidy_summary <- do.call(rbind, lapply(ploidy_by_harvest, function(df) {
    N <- as.numeric(df$N)
    pred_frac <- as.numeric(df$pred_fraction)
    obs_count <- as.numeric(df$obs_count)
    pred_mean_N <- sum(N * pred_frac, na.rm = TRUE) / max(sum(pred_frac, na.rm = TRUE), 1e-12)
    obs_mean_N <- sum(N * obs_count, na.rm = TRUE) / max(sum(obs_count, na.rm = TRUE), 1e-12)
    pred_mean_ploidy <- pred_mean_N / n_unit
    obs_mean_ploidy <- obs_mean_N / n_unit
    data.frame(
      harvest = df$harvest[1],
      cohort = df$cohort[1],
      pred_mean_N = pred_mean_N,
      obs_mean_N = obs_mean_N,
      pred_mean_ploidy = pred_mean_ploidy,
      obs_mean_ploidy = obs_mean_ploidy,
      data_fit_err_N = pred_mean_N - obs_mean_N,
      target_dist = dist_to_interval(pred_mean_ploidy, 2.35, 2.65),
      stringsAsFactors = FALSE
    )
  }))

  ploidy_2n <- ploidy_summary[ploidy_summary$cohort == "2N", , drop = FALSE]
  ploidy_4n <- ploidy_summary[ploidy_summary$cohort == "4N", , drop = FALSE]

  composite_tune_score <-
    as_num(summary_map[["objective_data"]], Inf) +
    8.0 * mean_or_na(ploidy_summary$target_dist) +
    10.0 * mean_or_na(ploidy_4n$target_dist) +
    2.0 * max(0, mean_or_na(latest_vol_dead_total_frac) - 0.35) +
    1.0 * max(0, abs(mean_or_na(log_err_2n)) - 0.30) +
    1.5 * max(0, abs(mean_or_na(log_err_4n)) - 0.30)

  list(
    fit_summary = fit_summary,
    best_params = best_params,
    burden_fit = burden_fit,
    terminal_ploidy = terminal_ploidy,
    summary_map = summary_map,
    best_map = best_map,
    ploidy_summary = ploidy_summary,
    metrics = list(
      objective = as_num(summary_map[["objective"]]),
      objective_data = as_num(summary_map[["objective_data"]]),
      objective_burden = as_num(summary_map[["objective_burden"]]),
      objective_ploidy = as_num(summary_map[["objective_ploidy"]]),
      burden_raw_r2 = raw_r2,
      burden_log_mae = mean(abs(log_err)),
      burden_log_rmse = sqrt(mean(log_err^2)),
      burden_mean_log_err = mean(log_err),
      burden_mean_log_err_2N = mean_or_na(log_err_2n),
      burden_mean_log_err_4N = mean_or_na(log_err_4n),
      latest_mean_live_frac = mean_or_na(latest_live_frac),
      latest_mean_dead_hypoxia_frac = mean_or_na(latest_dead_h_frac),
      latest_mean_dead_buffer_frac = mean_or_na(latest_dead_b_frac),
      latest_mean_dead_total_frac = mean_or_na(latest_dead_total_frac),
      latest_mean_live_volume_frac = mean_or_na(latest_vol_live_frac),
      latest_mean_dead_total_volume_frac = mean_or_na(latest_vol_dead_total_frac),
      latest_mean_dead_buffer_volume_frac = mean_or_na(latest_vol_dead_buffer_frac),
      pred_terminal_ploidy_mean_2N = mean_or_na(ploidy_2n$pred_mean_ploidy),
      pred_terminal_ploidy_mean_4N = mean_or_na(ploidy_4n$pred_mean_ploidy),
      obs_terminal_ploidy_mean_2N = mean_or_na(ploidy_2n$obs_mean_ploidy),
      obs_terminal_ploidy_mean_4N = mean_or_na(ploidy_4n$obs_mean_ploidy),
      ploidy_fit_mean_abs_err_2N = mean_or_na(abs(ploidy_2n$data_fit_err_N)),
      ploidy_fit_mean_abs_err_4N = mean_or_na(abs(ploidy_4n$data_fit_err_N)),
      target_range_mae_all = mean_or_na(ploidy_summary$target_dist),
      target_range_mae_2N = mean_or_na(ploidy_2n$target_dist),
      target_range_mae_4N = mean_or_na(ploidy_4n$target_dist),
      composite_tune_score = composite_tune_score
    )
  )
}

render_analysis <- function(version_label, analysis, previous_metrics = NULL, reasons, modifications,
                            stop_reason = NULL,
                            growth_penalty_ploidy_on = FALSE,
                            growth_penalty_hypoxia_on = TRUE,
                            death_on = FALSE) {
  m <- analysis$metrics
  bp <- analysis$best_map
  ps <- analysis$ploidy_summary
  lines <- c(
    paste("Version:", version_label),
    paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Target:",
    paste("  growth_penalty_hypoxia =", if (isTRUE(growth_penalty_hypoxia_on)) "TRUE" else "FALSE"),
    paste("  death =", if (isTRUE(death_on)) "TRUE" else "FALSE"),
    "  terminal ploidy target = approximately 2.5 for both 2N and 4N cohorts (operational band [2.35, 2.65])",
    "  burden and ploidy should both fit without burden being carried mainly by dead mass.",
    "",
    "Current metrics:",
    paste("  objective =", fmt_num(m$objective)),
    paste("  objective_data =", fmt_num(m$objective_data)),
    paste("  objective_burden =", fmt_num(m$objective_burden)),
    paste("  objective_ploidy =", fmt_num(m$objective_ploidy)),
    paste("  composite_tune_score =", fmt_num(m$composite_tune_score)),
    paste("  burden_raw_R2 =", fmt_num(m$burden_raw_r2)),
    paste("  burden_log_RMSE =", fmt_num(m$burden_log_rmse)),
    paste("  burden_mean_log_err_2N =", fmt_num(m$burden_mean_log_err_2N)),
    paste("  burden_mean_log_err_4N =", fmt_num(m$burden_mean_log_err_4N)),
    paste("  latest_mean_live_frac =", fmt_num(m$latest_mean_live_frac)),
    paste("  latest_mean_dead_hypoxia_frac =", fmt_num(m$latest_mean_dead_hypoxia_frac)),
    paste("  latest_mean_dead_total_frac =", fmt_num(m$latest_mean_dead_total_frac)),
    paste("  latest_mean_dead_total_volume_frac =", fmt_num(m$latest_mean_dead_total_volume_frac)),
    paste("  latest_mean_dead_buffer_volume_frac =", fmt_num(m$latest_mean_dead_buffer_volume_frac)),
    paste("  pred_terminal_ploidy_mean_2N =", fmt_num(m$pred_terminal_ploidy_mean_2N)),
    paste("  pred_terminal_ploidy_mean_4N =", fmt_num(m$pred_terminal_ploidy_mean_4N)),
    paste("  target_range_MAE_2N =", fmt_num(m$target_range_mae_2N)),
    paste("  target_range_MAE_4N =", fmt_num(m$target_range_mae_4N)),
    paste("  ploidy_fit_mean_abs_err_2N (N units) =", fmt_num(m$ploidy_fit_mean_abs_err_2N)),
    paste("  ploidy_fit_mean_abs_err_4N (N units) =", fmt_num(m$ploidy_fit_mean_abs_err_4N)),
    "",
    "Predicted terminal ploidy by harvest:",
    paste(apply(ps[, c("harvest", "cohort", "pred_mean_ploidy", "obs_mean_ploidy", "target_dist")], 1, function(r) {
      sprintf("  harvest=%s cohort=%s pred=%.4f obs=%.4f target_dist=%.4f", r[[1]], r[[2]], as.numeric(r[[3]]), as.numeric(r[[4]]), as.numeric(r[[5]]))
    }), collapse = "\n"),
    "",
    "Best parameters (selected):",
    paste("  p_misseg =", fmt_num(bp[["p_misseg"]])),
    paste("  k_o_mis =", fmt_num(bp[["k_o_mis"]])),
    paste("  beta_loss =", fmt_num(bp[["beta_loss"]])),
    paste("  gain_loss_ratio =", fmt_num(bp[["gain_loss_ratio"]])),
    paste("  beta_size =", fmt_num(bp[["beta_size"]])),
    paste("  alpha_o2 =", fmt_num(bp[["alpha_o2"]])),
    paste("  mu_hp =", fmt_num(bp[["mu_hp"]])),
    paste("  o2_init_pct =", fmt_num(bp[["o2_init_pct"]])),
    paste("  o2_rate =", fmt_num(bp[["o2_rate"]])),
    paste("  o2_shape_v =", fmt_num(bp[["o2_shape_v"]])),
    paste("  rho_2N =", fmt_num(bp[["rho_2N"]])),
    paste("  k_clear =", fmt_num(bp[["k_clear"]])),
    paste("  sigma_burden =", fmt_num(bp[["sigma_burden"]]))
  )
  if (isTRUE(growth_penalty_ploidy_on)) {
    lines <- c(lines, paste("  gamma_growth =", fmt_num(bp[["gamma_growth"]])))
  } else {
    lines <- c(lines, "  gamma_growth = inactive (growth_penalty_ploidy=FALSE)")
  }
  if (!is.null(previous_metrics)) {
    lines <- c(
      lines,
      "",
      "Change vs previous version:",
      paste("  delta_objective_data =", fmt_num(m$objective_data - previous_metrics$objective_data)),
      paste("  delta_composite_tune_score =", fmt_num(m$composite_tune_score - previous_metrics$composite_tune_score)),
      paste("  delta_target_range_MAE_2N =", fmt_num(m$target_range_mae_2N - previous_metrics$target_range_mae_2N)),
      paste("  delta_target_range_MAE_4N =", fmt_num(m$target_range_mae_4N - previous_metrics$target_range_mae_4N)),
      paste("  delta_burden_log_RMSE =", fmt_num(m$burden_log_rmse - previous_metrics$burden_log_rmse))
    )
  }
  lines <- c(lines, "", "Analysis reasoning:")
  if (length(reasons)) {
    lines <- c(lines, paste0("  - ", reasons))
  } else {
    lines <- c(lines, "  - No dominant failure mode was isolated from the current run.")
  }
  lines <- c(lines, "", "Next-step parameter plan:")
  if (length(modifications)) {
    lines <- c(lines, paste0("  - ", modifications))
  } else {
    lines <- c(lines, "  - No further micro-adjustments are currently justified.")
  }
  if (!is.null(stop_reason)) {
    lines <- c(lines, "", paste("Stop reason:", stop_reason))
  }
  lines
}

propose_next_version <- function(cfg_lines, tab, analysis, history_metrics,
                                 growth_penalty_hypoxia_on = TRUE,
                                 death_on = FALSE) {
  m <- analysis$metrics
  bp <- analysis$best_map
  reasons <- character()
  modifications <- character()
  tab_next <- tab
  cfg_next <- cfg_lines

  get_bp <- function(name, default = NA_real_) {
    x <- suppressWarnings(as.numeric(bp[[name]]))
    if (!is.finite(x)) default else x
  }

  p_bounds <- get_param_natural_bounds(tab_next, "log10_p_misseg")
  kmis_bounds <- get_param_natural_bounds(tab_next, "log10_k_o_mis")
  bloss_bounds <- get_param_natural_bounds(tab_next, "log10_beta_loss")
  ratio_bounds <- get_param_natural_bounds(tab_next, "logit_gain_loss_ratio")
  bsize_bounds <- get_param_natural_bounds(tab_next, "beta_size")
  alpha_bounds <- get_param_natural_bounds(tab_next, "log10_alpha_o2")
  o2init_bounds <- get_param_natural_bounds(tab_next, "log10_o2_init_pct")
  o2rate_bounds <- get_param_natural_bounds(tab_next, "log10_o2_rate")
  rho_bounds <- get_param_natural_bounds(tab_next, "log10_rho_2N")
  kclear_bounds <- get_param_natural_bounds(tab_next, "log10_k_clear")
  sigma_bounds <- get_param_natural_bounds(tab_next, "log10_sigma_burden")

  pred2 <- m$pred_terminal_ploidy_mean_2N
  pred4 <- m$pred_terminal_ploidy_mean_4N
  burden2 <- m$burden_mean_log_err_2N
  burden4 <- m$burden_mean_log_err_4N
  dead_buffer_vol <- m$latest_mean_dead_buffer_volume_frac

  target_lower <- 2.35
  target_upper <- 2.65
  split_overwide <- is.finite(pred2) && is.finite(pred4) && pred2 < target_lower && pred4 > target_upper
  both_high <- is.finite(pred2) && is.finite(pred4) && pred2 > target_upper && pred4 > target_upper
  both_low <- is.finite(pred2) && is.finite(pred4) && pred2 < target_lower && pred4 < target_lower

  if (isTRUE(growth_penalty_hypoxia_on) && !isTRUE(death_on)) {
    if (isTRUE(split_overwide)) {
      reasons <- c(reasons, "2N is over-converged while 4N remains too high; preserve global pressure floors and retime the pressure later in oxygen while penalizing high-copy gain daughters more strongly.")

      ratio_upper_new <- min(0.30, max(ratio_bounds$upper, max(get_bp("gain_loss_ratio", ratio_bounds$init) * 1.15, 0.20)))
      ratio_init_new <- min(ratio_upper_new * 0.92, max(ratio_bounds$lower, max(get_bp("gain_loss_ratio", ratio_bounds$init) * 1.08, 0.16)))
      tab_next <- set_param_natural(tab_next, "logit_gain_loss_ratio", init = ratio_init_new, upper = ratio_upper_new,
                                    note_append = "autotune: raised under 2N-low / 4N-high split")
      cfg_next <- set_prior_center(cfg_next, "prior_center_logit_gain_loss_ratio", qlogis(ratio_init_new))
      modifications <- c(modifications, sprintf("Raised gain_loss_ratio upper/init to %.4g / %.4g so high-copy gain daughters are less favored.", ratio_upper_new, ratio_init_new))

      kmis_lower_new <- max(0.12, min(kmis_bounds$lower, 0.15))
      kmis_upper_new <- min(0.50, max(kmis_bounds$lower, 0.45))
      kmis_init_new <- min(kmis_upper_new, max(kmis_lower_new, min(get_bp("k_o_mis", kmis_bounds$init) * 0.85, 0.30)))
      tab_next <- set_param_natural(tab_next, "log10_k_o_mis", init = kmis_init_new, lower = kmis_lower_new, upper = kmis_upper_new,
                                    note_append = "autotune: lower k_o_mis to shift misseg pressure later")
      modifications <- c(modifications, sprintf("Restricted k_o_mis to [%.4g, %.4g] with init %.4g to shift convergence pressure later.", kmis_lower_new, kmis_upper_new, kmis_init_new))

      o2_init_lower_new <- max(1.9, min(o2init_bounds$lower, 2.0))
      o2_init_upper_new <- min(2.5, max(o2init_bounds$upper, 2.5))
      o2_init_init_new <- min(o2_init_upper_new, max(o2_init_lower_new, min(get_bp("o2_init_pct", o2init_bounds$init) * 0.95, 2.15)))
      tab_next <- set_param_natural(tab_next, "log10_o2_init_pct", init = o2_init_init_new, lower = o2_init_lower_new, upper = o2_init_upper_new,
                                    note_append = "autotune: slightly lower initial oxygen under split failure")
      cfg_next <- yaml_set(cfg_next, "o2_init_pct_init", o2_init_init_new)
      cfg_next <- yaml_set(cfg_next, "o2_init_pct_min", o2_init_lower_new)
      cfg_next <- yaml_set(cfg_next, "o2_init_pct_max", o2_init_upper_new)
      cfg_next <- set_prior_center(cfg_next, "prior_center_log10_o2_init_pct", log10(o2_init_init_new))
      modifications <- c(modifications, sprintf("Shifted o2_init_pct to [%.4g, %.4g] with init %.4g to bring convergence pressure forward.", o2_init_lower_new, o2_init_upper_new, o2_init_init_new))

      o2_rate_lower_new <- max(0.45, min(o2rate_bounds$lower, 0.50))
      o2_rate_upper_new <- min(1.10, max(o2rate_bounds$upper, 1.00))
      o2_rate_init_new <- min(o2_rate_upper_new, max(o2_rate_lower_new, max(get_bp("o2_rate", o2rate_bounds$init) * 1.10, 0.75)))
      tab_next <- set_param_natural(tab_next, "log10_o2_rate", init = o2_rate_init_new, lower = o2_rate_lower_new, upper = o2_rate_upper_new,
                                    note_append = "autotune: faster oxygen decline under split failure")
      cfg_next <- yaml_set(cfg_next, "o2_rate_init", o2_rate_init_new)
      cfg_next <- yaml_set(cfg_next, "o2_rate_min", o2_rate_lower_new)
      cfg_next <- yaml_set(cfg_next, "o2_rate_max", o2_rate_upper_new)
      cfg_next <- set_prior_center(cfg_next, "prior_center_log10_o2_rate", log10(o2_rate_init_new))
      modifications <- c(modifications, sprintf("Shifted o2_rate to [%.4g, %.4g] with init %.4g to accelerate later hypoxia pressure.", o2_rate_lower_new, o2_rate_upper_new, o2_rate_init_new))

      alpha_lower_new <- max(alpha_bounds$lower, 0.7)
      alpha_upper_new <- min(2.1, max(alpha_bounds$upper, max(get_bp("alpha_o2", alpha_bounds$init) * 1.08, 1.3)))
      alpha_init_new <- min(alpha_upper_new, max(alpha_lower_new, max(get_bp("alpha_o2", alpha_bounds$init) * 1.05, 1.2)))
      tab_next <- set_param_natural(tab_next, "log10_alpha_o2", init = alpha_init_new, lower = alpha_lower_new, upper = alpha_upper_new,
                                    note_append = "autotune: modestly raised under split failure")
      cfg_next <- yaml_set(cfg_next, "alpha_o2_init", alpha_init_new)
      cfg_next <- yaml_set(cfg_next, "alpha_o2_min", alpha_lower_new)
      cfg_next <- yaml_set(cfg_next, "alpha_o2_max", alpha_upper_new)
      modifications <- c(modifications, sprintf("Raised alpha_o2 to [%.4g, %.4g] with init %.4g for modest high-ploidy damping.", alpha_lower_new, alpha_upper_new, alpha_init_new))

      bsize_lower_new <- max(0.30, min(bsize_bounds$lower, 0.35))
      bsize_upper_new <- min(0.52, max(bsize_bounds$lower, 0.50))
      bsize_init_new <- min(bsize_upper_new, max(bsize_lower_new, min(get_bp("beta_size", bsize_bounds$init) * 0.92, 0.42)))
      tab_next <- set_param_natural(tab_next, "beta_size", init = bsize_init_new, lower = bsize_lower_new, upper = bsize_upper_new,
                                    note_append = "autotune: trimmed under split failure to reduce 4N burden overweighting")
      cfg_next <- set_prior_center(cfg_next, "prior_center_beta_size", bsize_init_new)
      modifications <- c(modifications, sprintf("Shifted beta_size to [%.4g, %.4g] with init %.4g to reduce chronic 4N burden overweighting.", bsize_lower_new, bsize_upper_new, bsize_init_new))
    } else if (isTRUE(both_high)) {
      reasons <- c(reasons, "Both cohorts remain above target; increase late-stage convergence pressure without dropping p_misseg or beta_loss below their intended floors.")

      p_upper_new <- min(0.08, max(p_bounds$upper, max(get_bp("p_misseg", p_bounds$init) * 1.10, 0.07)))
      p_init_new <- min(p_upper_new * 0.95, max(p_bounds$lower, max(get_bp("p_misseg", p_bounds$init) * 1.05, 0.058)))
      tab_next <- set_param_natural(tab_next, "log10_p_misseg", init = p_init_new, upper = p_upper_new,
                                    note_append = "autotune: raised when both cohorts stay high")
      modifications <- c(modifications, sprintf("Raised p_misseg upper/init to %.4g / %.4g because both cohorts remain above target.", p_upper_new, p_init_new))

      bloss_upper_new <- min(0.15, max(bloss_bounds$upper, max(get_bp("beta_loss", bloss_bounds$init) * 1.10, 0.11)))
      bloss_init_new <- min(bloss_upper_new * 0.95, max(bloss_bounds$lower, max(get_bp("beta_loss", bloss_bounds$init) * 1.04, 0.095)))
      tab_next <- set_param_natural(tab_next, "log10_beta_loss", init = bloss_init_new, upper = bloss_upper_new,
                                    note_append = "autotune: raised when both cohorts stay high")
      cfg_next <- set_prior_center(cfg_next, "prior_center_log10_beta_loss", log10(bloss_init_new))
      modifications <- c(modifications, sprintf("Raised beta_loss upper/init to %.4g / %.4g to increase daughter-side convergence pressure.", bloss_upper_new, bloss_init_new))

      kmis_lower_new <- max(0.12, min(kmis_bounds$lower, 0.15))
      kmis_upper_new <- min(0.55, max(kmis_bounds$lower, 0.50))
      kmis_init_new <- min(kmis_upper_new, max(kmis_lower_new, min(get_bp("k_o_mis", kmis_bounds$init) * 0.90, 0.30)))
      tab_next <- set_param_natural(tab_next, "log10_k_o_mis", init = kmis_init_new, lower = kmis_lower_new, upper = kmis_upper_new,
                                    note_append = "autotune: lower k_o_mis when both cohorts stay high")
      modifications <- c(modifications, sprintf("Shifted k_o_mis to [%.4g, %.4g] with init %.4g to emphasize later hypoxia-driven convergence.", kmis_lower_new, kmis_upper_new, kmis_init_new))
    } else if (isTRUE(both_low)) {
      reasons <- c(reasons, "Both cohorts are below target; ease the weak hypoxia background and delay misseg pressure slightly, without collapsing the misseg floor.")

      kmis_lower_new <- max(0.20, min(kmis_bounds$lower, 0.20))
      kmis_upper_new <- min(0.70, max(kmis_bounds$upper, 0.60))
      kmis_init_new <- min(kmis_upper_new, max(kmis_lower_new, max(get_bp("k_o_mis", kmis_bounds$init) * 1.10, 0.45)))
      tab_next <- set_param_natural(tab_next, "log10_k_o_mis", init = kmis_init_new, lower = kmis_lower_new, upper = kmis_upper_new,
                                    note_append = "autotune: raise k_o_mis when both cohorts are low")
      modifications <- c(modifications, sprintf("Shifted k_o_mis to [%.4g, %.4g] with init %.4g to delay low-oxygen convergence pressure.", kmis_lower_new, kmis_upper_new, kmis_init_new))

      alpha_upper_new <- min(alpha_bounds$upper, max(alpha_bounds$lower, max(get_bp("alpha_o2", alpha_bounds$init) * 0.95, 1.0)))
      alpha_init_new <- min(alpha_upper_new * 0.95, max(alpha_bounds$lower, max(get_bp("alpha_o2", alpha_bounds$init) * 0.92, 0.85)))
      tab_next <- set_param_natural(tab_next, "log10_alpha_o2", init = alpha_init_new, upper = alpha_upper_new,
                                    note_append = "autotune: softened when both cohorts are low")
      cfg_next <- yaml_set(cfg_next, "alpha_o2_init", alpha_init_new)
      cfg_next <- yaml_set(cfg_next, "alpha_o2_max", alpha_upper_new)
      modifications <- c(modifications, sprintf("Trimmed alpha_o2 upper/init to %.4g / %.4g because both cohorts are converging too low.", alpha_upper_new, alpha_init_new))
    } else if (is.finite(pred4) && pred4 > target_upper) {
      reasons <- c(reasons, "4N terminal ploidy remains too high; push later, high-ploidy-directed pressure before relaxing core misseg or loss parameters.")

      ratio_upper_new <- min(0.30, max(ratio_bounds$upper, max(get_bp("gain_loss_ratio", ratio_bounds$init) * 1.12, 0.18)))
      ratio_init_new <- min(ratio_upper_new * 0.94, max(ratio_bounds$lower, max(get_bp("gain_loss_ratio", ratio_bounds$init) * 1.06, 0.15)))
      tab_next <- set_param_natural(tab_next, "logit_gain_loss_ratio", init = ratio_init_new, upper = ratio_upper_new,
                                    note_append = "autotune: raise gain penalty when 4N stays high")
      cfg_next <- set_prior_center(cfg_next, "prior_center_logit_gain_loss_ratio", qlogis(ratio_init_new))
      modifications <- c(modifications, sprintf("Raised gain_loss_ratio upper/init to %.4g / %.4g because 4N remains above target.", ratio_upper_new, ratio_init_new))

      if (is.finite(burden4) && burden4 > 0.20) {
        bsize_upper_new <- min(bsize_bounds$upper, max(bsize_bounds$lower, max(get_bp("beta_size", bsize_bounds$init) * 0.95, 0.50)))
        bsize_init_new <- min(bsize_upper_new, max(bsize_bounds$lower, max(get_bp("beta_size", bsize_bounds$init) * 0.90, 0.42)))
        tab_next <- set_param_natural(tab_next, "beta_size", init = bsize_init_new, upper = bsize_upper_new,
                                      note_append = "autotune: trim beta_size when 4N burden stays high")
        cfg_next <- set_prior_center(cfg_next, "prior_center_beta_size", bsize_init_new)
        modifications <- c(modifications, sprintf("Trimmed beta_size upper/init to %.4g / %.4g because 4N burden remains selectively high.", bsize_upper_new, bsize_init_new))
      }
    } else if (is.finite(pred2) && pred2 < target_lower) {
      reasons <- c(reasons, "2N terminal ploidy is below target while 4N is not persistently high; protect low-ploidy states by delaying misseg pressure rather than softening the global misseg/loss floor.")

      kmis_lower_new <- max(0.20, min(kmis_bounds$lower, 0.20))
      kmis_upper_new <- min(0.70, max(kmis_bounds$upper, 0.60))
      kmis_init_new <- min(kmis_upper_new, max(kmis_lower_new, max(get_bp("k_o_mis", kmis_bounds$init) * 1.10, 0.45)))
      tab_next <- set_param_natural(tab_next, "log10_k_o_mis", init = kmis_init_new, lower = kmis_lower_new, upper = kmis_upper_new,
                                    note_append = "autotune: delay misseg pressure when only 2N overshoots low")
      modifications <- c(modifications, sprintf("Shifted k_o_mis to [%.4g, %.4g] with init %.4g to ease early pressure on 2N states.", kmis_lower_new, kmis_upper_new, kmis_init_new))

      alpha_upper_new <- min(alpha_bounds$upper, max(alpha_bounds$lower, max(get_bp("alpha_o2", alpha_bounds$init) * 0.95, 0.9)))
      alpha_init_new <- min(alpha_upper_new * 0.95, max(alpha_bounds$lower, max(get_bp("alpha_o2", alpha_bounds$init) * 0.92, 0.9)))
      tab_next <- set_param_natural(tab_next, "log10_alpha_o2", init = alpha_init_new, upper = alpha_upper_new,
                                    note_append = "autotune: slightly soften alpha when only 2N overshoots low")
      cfg_next <- yaml_set(cfg_next, "alpha_o2_init", alpha_init_new)
      cfg_next <- yaml_set(cfg_next, "alpha_o2_max", alpha_upper_new)
      modifications <- c(modifications, sprintf("Trimmed alpha_o2 upper/init to %.4g / %.4g because low-ploidy states are overshooting downward.", alpha_upper_new, alpha_init_new))

      o2_init_lower_new <- max(2.1, min(o2init_bounds$lower, 2.1))
      o2_init_upper_new <- min(2.8, max(o2init_bounds$upper, 2.6))
      o2_init_init_new <- min(o2_init_upper_new, max(o2_init_lower_new, max(get_bp("o2_init_pct", o2init_bounds$init) * 1.03, 2.35)))
      tab_next <- set_param_natural(tab_next, "log10_o2_init_pct", init = o2_init_init_new, lower = o2_init_lower_new, upper = o2_init_upper_new,
                                    note_append = "autotune: preserve early oxygen when only 2N overshoots low")
      cfg_next <- yaml_set(cfg_next, "o2_init_pct_init", o2_init_init_new)
      cfg_next <- yaml_set(cfg_next, "o2_init_pct_min", o2_init_lower_new)
      cfg_next <- yaml_set(cfg_next, "o2_init_pct_max", o2_init_upper_new)
      cfg_next <- set_prior_center(cfg_next, "prior_center_log10_o2_init_pct", log10(o2_init_init_new))
      modifications <- c(modifications, sprintf("Shifted o2_init_pct to [%.4g, %.4g] with init %.4g to preserve low-ploidy growth earlier.", o2_init_lower_new, o2_init_upper_new, o2_init_init_new))
    }
  }

  if (is.finite(dead_buffer_vol) && dead_buffer_vol > 0.25) {
    reasons <- c(reasons, "Dead-buffer volume is too high; raise clearance first rather than softening p_misseg or beta_loss below their intended floors.")
    kclear_lower_new <- max(kclear_bounds$lower, 1.2e-03)
    kclear_upper_new <- min(4e-03, max(kclear_bounds$upper, get_bp("k_clear", kclear_bounds$init) * 1.25, 2.5e-03))
    kclear_init_new <- min(kclear_upper_new * 0.90, max(kclear_lower_new, get_bp("k_clear", kclear_bounds$init) * 1.15, 1.8e-03))
    tab_next <- set_param_natural(tab_next, "log10_k_clear", init = kclear_init_new, lower = kclear_lower_new, upper = kclear_upper_new,
                                  note_append = "autotune: increased after excessive dead-buffer contribution")
    cfg_next <- yaml_set(cfg_next, "k_clear_init", kclear_init_new)
    cfg_next <- yaml_set(cfg_next, "k_clear_min", kclear_lower_new)
    cfg_next <- yaml_set(cfg_next, "k_clear_max", kclear_upper_new)
    cfg_next <- set_prior_center(cfg_next, "prior_center_log10_k_clear", log10(kclear_init_new))
    modifications <- c(modifications, sprintf("Raised k_clear range to [%.4g, %.4g] with init %.4g to prevent dead mass from dominating burden.", kclear_lower_new, kclear_upper_new, kclear_init_new))
  }

  if (is.finite(burden2) && is.finite(burden4) && burden2 < -0.25 && burden4 < -0.25) {
    reasons <- c(reasons, "Both cohorts under-predict burden; overall live-volume support is still too weak.")
    alpha_upper_new <- min(alpha_bounds$upper, max(alpha_bounds$lower, max(get_bp("alpha_o2", alpha_bounds$init) * 0.95, 0.9)))
    alpha_init_new <- max(alpha_bounds$lower, min(alpha_upper_new, get_bp("alpha_o2", alpha_bounds$init) * 0.92))
    tab_next <- set_param_natural(tab_next, "log10_alpha_o2", init = alpha_init_new, upper = alpha_upper_new,
                                  note_append = "autotune: lowered after global burden underfit")
    cfg_next <- yaml_set(cfg_next, "alpha_o2_init", alpha_init_new)
    cfg_next <- yaml_set(cfg_next, "alpha_o2_max", alpha_upper_new)
    modifications <- c(modifications, sprintf("Lowered alpha_o2 to init %.4g with upper %.4g because burden is underfit in both cohorts.", alpha_init_new, alpha_upper_new))

    rho_lower_new <- max(32000, min(rho_bounds$lower, 36000))
    rho_init_new <- max(rho_lower_new, min(rho_bounds$upper, min(get_bp("rho_2N", rho_bounds$init) * 0.94, 42000)))
    tab_next <- set_param_natural(tab_next, "log10_rho_2N", init = rho_init_new, lower = rho_lower_new, upper = rho_bounds$upper,
                                  note_append = "autotune: lowered after global burden underfit")
    cfg_next <- set_prior_center(cfg_next, "prior_center_log10_rho_2N", log10(rho_init_new))
    modifications <- c(modifications, sprintf("Lowered rho_2N init to %.4g to allow larger volume per cell under global burden underfit.", rho_init_new))
  }

  if (is.finite(burden4) && burden4 > 0.25 && (!is.finite(burden2) || burden4 > burden2 + 0.10)) {
    reasons <- c(reasons, "4N burden remains selectively high; keep trimming beta_size because the burden map is still overweighting high-ploidy states.")
    bsize_upper_new <- min(bsize_bounds$upper, max(bsize_bounds$lower, max(get_bp("beta_size", bsize_bounds$init) * 0.95, 0.50)))
    bsize_init_new <- min(bsize_upper_new, max(bsize_bounds$lower, max(get_bp("beta_size", bsize_bounds$init) * 0.90, 0.42)))
    tab_next <- set_param_natural(tab_next, "beta_size", init = bsize_init_new, upper = bsize_upper_new,
                                  note_append = "autotune: trimmed when 4N burden stays selectively high")
    cfg_next <- set_prior_center(cfg_next, "prior_center_beta_size", bsize_init_new)
    modifications <- c(modifications, sprintf("Trimmed beta_size upper/init to %.4g / %.4g because 4N burden remains selectively high.", bsize_upper_new, bsize_init_new))
  }

  if (is.finite(bp[["sigma_burden"]]) && is.finite(sigma_bounds$upper) &&
      param_hit_upper(bp[["sigma_burden"]], sigma_bounds$upper, tol_frac = 0.05) &&
      is.finite(m$objective_burden) && m$objective_burden > 1.0) {
    reasons <- c(reasons, "sigma_burden is near its ceiling while burden mismatch remains large; tighten burden slack before giving the optimizer more room to ignore burden.")
    sigma_upper_new <- max(0.20, min(0.24, sigma_bounds$upper * 0.92))
    sigma_init_new <- min(sigma_upper_new * 0.95, max(sigma_bounds$lower, 0.20))
    tab_next <- set_param_natural(tab_next, "log10_sigma_burden", init = sigma_init_new, upper = sigma_upper_new,
                                  note_append = "autotune: tightened because burden error remained high at sigma ceiling")
    cfg_next <- yaml_set(cfg_next, "sigma_burden", sigma_init_new)
    cfg_next <- yaml_set(cfg_next, "sigma_burden_max", sigma_upper_new)
    modifications <- c(modifications, sprintf("Reduced sigma_burden upper/init to %.4g / %.4g so burden cannot be ignored as easily.", sigma_upper_new, sigma_init_new))
  }

  if (!length(modifications) && length(history_metrics) >= 2L) {
    reasons <- c(reasons, "No sharper parameter-direction signal emerged from the latest diagnostics.")
  }

  list(cfg_lines = cfg_next, param_table = tab_next, reasons = unique(reasons), modifications = unique(modifications))
}

write_history_tsv <- function(history, path) {
  if (!length(history)) return(invisible(NULL))
  keys <- unique(unlist(lapply(history, names), use.names = FALSE))
  rows <- lapply(seq_along(history), function(i) {
    vals <- history[[i]]
    out <- setNames(as.list(rep(NA_character_, length(keys))), keys)
    for (k in names(vals)) out[[k]] <- as.character(vals[[k]])
    out[["version"]] <- paste0("v", i)
    out
  })
  df <- do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
  df <- df[, c("version", setdiff(names(df), "version")), drop = FALSE]
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(df)
}

write_best_version <- function(best_version, best_score, best_obj_data, path) {
  writeLines(c(
    paste("best_version", best_version %||% "NA", sep = "\t"),
    paste("best_composite_tune_score", if (is.finite(best_score)) fmt_num(best_score) else "NA", sep = "\t"),
    paste("best_objective_data", if (is.finite(best_obj_data)) fmt_num(best_obj_data) else "NA", sep = "\t")
  ), path, useBytes = TRUE)
}

  main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  oxygen_root <- file.path(repo_root, "oxygen")

  base_config <- normalizePath(
    argv$base_config %||% file.path(oxygen_root, "config", "O2_NGLF_MAP_asymmetric_intrinsic_buffer.yaml"),
    mustWork = TRUE
  )
  seed <- as_int(argv$seed, 1L)
  max_versions <- as_int(argv$max_versions, 6L)
  min_improve <- as_num(argv$min_improve, 0.25)
  initial_params_tsv <- argv$initial_params_tsv %||% NULL
  growth_penalty_ploidy_on <- as_flag(argv$growth_penalty_ploidy, FALSE)
  growth_penalty_hypoxia_on <- as_flag(argv$growth_penalty_hypoxia, TRUE)
  death_on <- as_flag(argv$death, FALSE)
  runner <- normalizePath(
    argv$runner %||% file.path(oxygen_root, "code", "O2_NGLF_MAP_asymmetric_intrinsic_buffer", "run_fit_invivo_model_O2_NGLF_MAP_asymmetric_intrinsic_buffer.sh"),
    mustWork = TRUE
  )
  viz_script <- normalizePath(
    file.path(oxygen_root, "code", "O2_NGLF_MAP_asymmetric_intrinsic_buffer", "viz_invivo_model_O2_NGLF_MAP_asymmetric_intrinsic_buffer_results.R"),
    mustWork = TRUE
  )
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  session_root <- argv$session_root %||% file.path(oxygen_root, "results", paste0("O2_NGLF_MAP_asymmetric_intrinsic_buffer_autotune_", stamp))
  session_root <- normalizePath(session_root, mustWork = FALSE)
  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)

  base_lines <- read_yaml_lines(base_config)
  base_dir <- dirname(base_config)
  base_param_path <- resolve_path_local(yaml_get(base_lines, "parameter_table"), base_dir)
  if (is.null(base_param_path) || !file.exists(base_param_path)) stop("Could not resolve base parameter_table from ", base_config)
  data_dir_abs <- resolve_path_local(yaml_get(base_lines, "data_dir"), base_dir)
  if (is.null(data_dir_abs) || !dir.exists(data_dir_abs)) stop("Could not resolve data_dir from ", base_config)

  seeded <- apply_recommended_start(
    base_lines,
    read_parameter_table(base_param_path),
    growth_penalty_ploidy_on = growth_penalty_ploidy_on,
    growth_penalty_hypoxia_on = growth_penalty_hypoxia_on,
    death_on = death_on
  )
  working_cfg <- seeded$cfg_lines
  working_tab <- seeded$param_table

  history <- list()
  no_improve <- 0L
  best_version <- NA_character_
  best_score <- Inf
  best_objective_data <- Inf
  history_path <- file.path(session_root, "autotune_history.tsv")
  best_path <- file.path(session_root, "best_version.txt")

  message("Autotune session root: ", session_root)
  message("Base config: ", base_config)
  message("Base parameter table: ", base_param_path)

  for (v in seq_len(max_versions)) {
    version_label <- paste0("v", v)
    version_cfg <- working_cfg
    version_cfg <- yaml_set(version_cfg, "run_prefix", version_label)
    version_cfg <- yaml_set(version_cfg, "append_run_prefix_timestamp", FALSE)
    version_cfg <- yaml_set(version_cfg, "out_root", ".")
    version_cfg <- yaml_set(version_cfg, "data_dir", data_dir_abs)
    version_cfg <- yaml_set(version_cfg, "seeds_file", "")
    version_cfg <- yaml_set(version_cfg, "seeds_csv", as.character(seed))
    version_cfg <- yaml_set(version_cfg, "auto_viz", FALSE)
    version_cfg <- yaml_set(version_cfg, "parameter_table", paste0(version_label, "_parameter_table.csv"))

    cfg_path <- file.path(session_root, paste0(version_label, "_config.yaml"))
    tab_path <- file.path(session_root, paste0(version_label, "_parameter_table.csv"))
    analysis_path <- file.path(session_root, paste0(version_label, "_analysis.txt"))
    run_dir <- file.path(session_root, version_label)

    writeLines(version_cfg, cfg_path, useBytes = TRUE)
    write.table(working_tab, file = tab_path, sep = ",", quote = TRUE, row.names = FALSE, col.names = TRUE)

    cmd <- c(
      runner,
      paste0("--config=", cfg_path),
      "--seeds_file=",
      paste0("--seeds_csv=", seed),
      "--auto_viz=FALSE",
      "--append_run_prefix_timestamp=FALSE"
    )
    if (v == 1L && !is.null(initial_params_tsv) && nzchar(initial_params_tsv)) {
      init_tsv_abs <- normalizePath(initial_params_tsv, mustWork = TRUE)
      cmd <- c(cmd, paste0("--init_params_tsv=", init_tsv_abs))
    } else if (v > 1L) {
      prev_ckpt <- file.path(session_root, paste0("v", v - 1L), "seed1", "checkpoints", "best_params_transformed_latest.tsv")
      if (file.exists(prev_ckpt)) cmd <- c(cmd, paste0("--init_params_tsv=", prev_ckpt))
    }

    message("[autotune] Running ", version_label, " ...")
    status <- system2("bash", cmd)
    if (!identical(status, 0L)) stop("Autotune run failed for ", version_label, " with exit status ", status)
    if (!dir.exists(run_dir)) stop("Expected run directory missing: ", run_dir)

    viz_cmd <- c(
      viz_script,
      paste0("--fit_dir=", file.path(run_dir, "seed1")),
      paste0("--data_dir=", data_dir_abs),
      "--report_dt=1",
      "--top_n=6",
      "--predict_horizons=100",
      "--predict_report_dt=1",
      "--predict_top_n=6",
      "--n_cores=1"
    )
    message("[autotune] Running light viz for ", version_label, " ...")
    viz_status <- system2("Rscript", viz_cmd)
    if (!identical(viz_status, 0L)) stop("Autotune viz failed for ", version_label, " with exit status ", viz_status)

    analysis <- analyze_run(run_dir)
    metrics <- c(list(version = version_label), analysis$metrics)
    history[[length(history) + 1L]] <- metrics

    previous_metrics <- if (length(history) >= 2L) history[[length(history) - 1L]] else NULL
    proposal <- propose_next_version(
      working_cfg,
      working_tab,
      analysis,
      history,
      growth_penalty_hypoxia_on = growth_penalty_hypoxia_on,
      death_on = death_on
    )

    score <- as_num(analysis$metrics$composite_tune_score, Inf)
    data_obj <- as_num(analysis$metrics$objective_data, Inf)
    if (is.finite(score) && score + min_improve < best_score) {
      best_score <- score
      best_objective_data <- data_obj
      best_version <- version_label
      no_improve <- 0L
    } else if (is.finite(score)) {
      no_improve <- no_improve + 1L
    }

    write_history_tsv(history, history_path)
    write_best_version(best_version, best_score, best_objective_data, best_path)

    stop_reason <- NULL
    if (!length(proposal$modifications)) {
      stop_reason <- "No further parameter micro-adjustments were indicated by the current diagnostics."
    } else if (no_improve >= 2L) {
      stop_reason <- paste0("Composite tune score failed to improve by at least ", min_improve, " for two consecutive versions.")
    } else if (v >= max_versions) {
      stop_reason <- paste0("Reached max_versions=", max_versions, ".")
    }

    analysis_lines <- render_analysis(
      version_label = version_label,
      analysis = analysis,
      previous_metrics = previous_metrics,
      reasons = proposal$reasons,
      modifications = proposal$modifications,
      stop_reason = stop_reason,
      growth_penalty_ploidy_on = growth_penalty_ploidy_on,
      growth_penalty_hypoxia_on = growth_penalty_hypoxia_on,
      death_on = death_on
    )
    writeLines(analysis_lines, analysis_path, useBytes = TRUE)
    message("[autotune] Wrote analysis: ", analysis_path)

    if (!is.null(stop_reason)) break

    working_cfg <- proposal$cfg_lines
    working_tab <- proposal$param_table
  }

  write_history_tsv(history, history_path)
  write_best_version(best_version, best_score, best_objective_data, best_path)
  message("Autotune complete. Session root: ", session_root)
  message("Best version: ", best_version %||% "NA", " (composite_tune_score=", fmt_num(best_score), ")")
}

if (sys.nframe() == 0L) main()
