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

read_yaml_lines <- function(path) {
  readLines(path, warn = FALSE)
}

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

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

analyze_run <- function(run_dir) {
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
  latest_live_frac <- with(latest_df, as.numeric(pred_pop_live) / pmax(as.numeric(pred_pop), 1e-12))
  latest_dead_buffer_frac <- with(latest_df, as.numeric(pred_pop_dead_buffer) / pmax(as.numeric(pred_pop), 1e-12))
  latest_dead_hyp_frac <- with(latest_df, as.numeric(pred_pop_dead_hypoxia) / pmax(as.numeric(pred_pop), 1e-12))
  latest_df_2n <- latest_df[latest_df$cohort == "2N", , drop = FALSE]
  latest_df_4n <- latest_df[latest_df$cohort == "4N", , drop = FALSE]
  latest_live_frac_2n <- with(latest_df_2n, as.numeric(pred_pop_live) / pmax(as.numeric(pred_pop), 1e-12))
  latest_live_frac_4n <- with(latest_df_4n, as.numeric(pred_pop_live) / pmax(as.numeric(pred_pop), 1e-12))
  latest_dead_buffer_frac_2n <- with(latest_df_2n, as.numeric(pred_pop_dead_buffer) / pmax(as.numeric(pred_pop), 1e-12))
  latest_dead_buffer_frac_4n <- with(latest_df_4n, as.numeric(pred_pop_dead_buffer) / pmax(as.numeric(pred_pop), 1e-12))
  day25_df <- burden_fit[abs(as.numeric(burden_fit$day) - 25.0) < 1e-9, , drop = FALSE]
  day25_live_frac <- with(day25_df, as.numeric(pred_pop_live) / pmax(as.numeric(pred_pop), 1e-12))

  burden_eval <- burden_fit[as.numeric(burden_fit$day) > 0, , drop = FALSE]
  obs <- as.numeric(burden_eval$obs_burden)
  pred <- as.numeric(burden_eval$pred_burden_volume_mm3)
  obs_log <- log(pmax(obs, 1e-12))
  pred_log <- log(pmax(pred, 1e-12))
  log_err <- pred_log - obs_log
  burden_eval_2n <- burden_eval[burden_eval$cohort == "2N", , drop = FALSE]
  burden_eval_4n <- burden_eval[burden_eval$cohort == "4N", , drop = FALSE]
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
    pred_mean <- sum(N * pred_frac, na.rm = TRUE) / max(sum(pred_frac, na.rm = TRUE), 1e-12)
    obs_mean <- sum(N * obs_count, na.rm = TRUE) / max(sum(obs_count, na.rm = TRUE), 1e-12)
    data.frame(
      harvest = df$harvest[1],
      cohort = df$cohort[1],
      pred_mean_N = pred_mean,
      obs_mean_N = obs_mean,
      mean_err_N = pred_mean - obs_mean,
      stringsAsFactors = FALSE
    )
  }))

  list(
    fit_summary = fit_summary,
    best_params = best_params,
    burden_fit = burden_fit,
    terminal_ploidy = terminal_ploidy,
    summary_map = summary_map,
    best_map = best_map,
    metrics = list(
      objective = as_num(summary_map[["objective"]]),
      objective_data = as_num(summary_map[["objective_data"]]),
      objective_burden = as_num(summary_map[["objective_burden"]]),
      objective_ploidy = as_num(summary_map[["objective_ploidy"]]),
      latest_mean_live_frac = mean_or_na(latest_live_frac),
      latest_mean_live_frac_2N = mean_or_na(latest_live_frac_2n),
      latest_mean_live_frac_4N = mean_or_na(latest_live_frac_4n),
      latest_mean_dead_buffer_frac = mean_or_na(latest_dead_buffer_frac),
      latest_mean_dead_buffer_frac_2N = mean_or_na(latest_dead_buffer_frac_2n),
      latest_mean_dead_buffer_frac_4N = mean_or_na(latest_dead_buffer_frac_4n),
      latest_mean_dead_hypoxia_frac = mean_or_na(latest_dead_hyp_frac),
      day25_mean_live_frac = mean_or_na(day25_live_frac),
      burden_raw_r2 = raw_r2,
      burden_log_mae = mean(abs(log_err)),
      burden_log_rmse = sqrt(mean(log_err^2)),
      burden_mean_log_err = mean(log_err),
      burden_mean_log_err_2N = mean_or_na(log_err_2n),
      burden_mean_log_err_4N = mean_or_na(log_err_4n),
      ploidy_mean_abs_err_2N = mean(abs(ploidy_summary$mean_err_N[ploidy_summary$cohort == "2N"])),
      ploidy_mean_abs_err_4N = mean(abs(ploidy_summary$mean_err_N[ploidy_summary$cohort == "4N"])),
      ploidy_mean_signed_err_2N = mean(ploidy_summary$mean_err_N[ploidy_summary$cohort == "2N"]),
      ploidy_mean_signed_err_4N = mean(ploidy_summary$mean_err_N[ploidy_summary$cohort == "4N"])
    ),
    ploidy_summary = ploidy_summary
  )
}

render_analysis <- function(version_label, analysis, previous_metrics = NULL, reasons, modifications, stop_reason = NULL) {
  m <- analysis$metrics
  bp <- analysis$best_map
  lines <- c(
    paste("Version:", version_label),
    paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Current metrics:",
    paste("  objective =", fmt_num(m$objective)),
    paste("  objective_data =", fmt_num(m$objective_data)),
    paste("  objective_burden =", fmt_num(m$objective_burden)),
    paste("  objective_ploidy =", fmt_num(m$objective_ploidy)),
    paste("  latest_mean_live_frac =", fmt_num(m$latest_mean_live_frac)),
    paste("  latest_mean_live_frac_2N =", fmt_num(m$latest_mean_live_frac_2N)),
    paste("  latest_mean_live_frac_4N =", fmt_num(m$latest_mean_live_frac_4N)),
    paste("  latest_mean_dead_buffer_frac =", fmt_num(m$latest_mean_dead_buffer_frac)),
    paste("  latest_mean_dead_buffer_frac_2N =", fmt_num(m$latest_mean_dead_buffer_frac_2N)),
    paste("  latest_mean_dead_buffer_frac_4N =", fmt_num(m$latest_mean_dead_buffer_frac_4N)),
    paste("  latest_mean_dead_hypoxia_frac =", fmt_num(m$latest_mean_dead_hypoxia_frac)),
    paste("  burden_raw_R2 =", fmt_num(m$burden_raw_r2)),
    paste("  burden_log_MAE =", fmt_num(m$burden_log_mae)),
    paste("  burden_log_RMSE =", fmt_num(m$burden_log_rmse)),
    paste("  burden_mean_log_err =", fmt_num(m$burden_mean_log_err)),
    paste("  burden_mean_log_err_2N =", fmt_num(m$burden_mean_log_err_2N)),
    paste("  burden_mean_log_err_4N =", fmt_num(m$burden_mean_log_err_4N)),
    paste("  ploidy_mean_abs_err_2N =", fmt_num(m$ploidy_mean_abs_err_2N)),
    paste("  ploidy_mean_abs_err_4N =", fmt_num(m$ploidy_mean_abs_err_4N)),
    paste("  ploidy_mean_signed_err_2N =", fmt_num(m$ploidy_mean_signed_err_2N)),
    paste("  ploidy_mean_signed_err_4N =", fmt_num(m$ploidy_mean_signed_err_4N)),
    "",
    "Best parameters (selected):",
    paste("  p_misseg =", fmt_num(bp[["p_misseg"]])),
    paste("  beta_buffer =", fmt_num(bp[["beta_buffer"]])),
    paste("  smax =", fmt_num(bp[["smax"]])),
    paste("  beta_loss =", fmt_num(bp[["beta_loss"]])),
    paste("  gain_loss_ratio =", fmt_num(bp[["gain_loss_ratio"]])),
    paste("  k_o_mis =", fmt_num(bp[["k_o_mis"]])),
    paste("  beta_size =", fmt_num(bp[["beta_size"]])),
    paste("  alpha_o2 =", fmt_num(bp[["alpha_o2"]])),
    paste("  o2_init_pct =", fmt_num(bp[["o2_init_pct"]])),
    paste("  o2_rate =", fmt_num(bp[["o2_rate"]])),
    paste("  k_clear =", fmt_num(bp[["k_clear"]])),
    paste("  sigma_burden =", fmt_num(bp[["sigma_burden"]]))
  )
  if (!is.null(previous_metrics)) {
    lines <- c(
      lines,
      "",
      "Change vs previous version:",
      paste("  delta_objective_data =", fmt_num(m$objective_data - previous_metrics$objective_data)),
      paste("  delta_objective_burden =", fmt_num(m$objective_burden - previous_metrics$objective_burden)),
      paste("  delta_objective_ploidy =", fmt_num(m$objective_ploidy - previous_metrics$objective_ploidy)),
      paste("  delta_latest_mean_live_frac =", fmt_num(m$latest_mean_live_frac - previous_metrics$latest_mean_live_frac))
    )
  }
  lines <- c(lines, "", "Analysis reasons:")
  if (length(reasons)) {
    lines <- c(lines, paste0("  - ", reasons))
  } else {
    lines <- c(lines, "  - No strong failure mode detected.")
  }
  lines <- c(lines, "", "Proposed modifications for next version:")
  if (length(modifications)) {
    lines <- c(lines, paste0("  - ", modifications))
  } else {
    lines <- c(lines, "  - No further parameter changes proposed.")
  }
  if (!is.null(stop_reason)) {
    lines <- c(lines, "", paste("Stop reason:", stop_reason))
  }
  lines
}

propose_next_version <- function(cfg_lines, tab, analysis, history_metrics) {
  m <- analysis$metrics
  bp <- analysis$best_map
  reasons <- character()
  modifications <- character()
  tab_next <- tab
  cfg_next <- cfg_lines

  p_bounds <- get_param_natural_bounds(tab_next, "log10_p_misseg")
  bbuf_bounds <- get_param_natural_bounds(tab_next, "beta_buffer")
  smax_bounds <- get_param_natural_bounds(tab_next, "log10_smax")
  bloss_bounds <- get_param_natural_bounds(tab_next, "log10_beta_loss")
  kmis_bounds <- get_param_natural_bounds(tab_next, "log10_k_o_mis")
  bsize_bounds <- get_param_natural_bounds(tab_next, "beta_size")
  o2init_bounds <- get_param_natural_bounds(tab_next, "log10_o2_init_pct")
  o2rate_bounds <- get_param_natural_bounds(tab_next, "log10_o2_rate")
  alpha_bounds <- get_param_natural_bounds(tab_next, "log10_alpha_o2")
  sigma_bounds <- get_param_natural_bounds(tab_next, "log10_sigma_burden")

  if (is.finite(m$ploidy_mean_signed_err_4N) && m$ploidy_mean_signed_err_4N > 20 &&
      is.finite(m$latest_mean_live_frac_4N) && m$latest_mean_live_frac_4N > 0.70) {
    reasons <- c(reasons, "4N terminal ploidy remains too high while 4N cells stay viable; the next version should preserve early 4N burden and shift stronger ploidy movement later into lower-oxygen phases.")
    upper_new <- min(0.07, max(p_bounds$upper, bp[["p_misseg"]] * 1.10, 0.06))
    init_new <- min(upper_new * 0.90, max(bp[["p_misseg"]] * 1.05, p_bounds$init, 0.045))
    tab_next <- set_param_natural(
      tab_next, "log10_p_misseg",
      init = init_new,
      upper = upper_new,
      note_append = "autotune: kept high so late hypoxia-coupled missegregation can still drive 4N convergence"
    )
    modifications <- c(modifications, sprintf("Raised p_misseg upper to %.4g and init to %.4g so late ploidy mobility remains available.", upper_new, init_new))

    kmis_lower_new <- max(kmis_bounds$lower, 0.2)
    kmis_upper_new <- min(1.0, max(kmis_lower_new, 0.8))
    kmis_init_new <- min(kmis_upper_new * 0.80, max(kmis_lower_new, 0.5))
    tab_next <- set_param_natural(
      tab_next, "log10_k_o_mis",
      init = kmis_init_new,
      lower = kmis_lower_new,
      upper = kmis_upper_new,
      note_append = "autotune: lowered so missegregation stays weaker earlier and turns on later as oxygen falls"
    )
    modifications <- c(modifications, sprintf("Lowered k_o_mis range to [%.4g, %.4g] with init %.4g to delay strong missegregation until lower oxygen.", kmis_lower_new, kmis_upper_new, kmis_init_new))

    alpha_lower_new <- max(alpha_bounds$lower, 8e-05)
    alpha_upper_new <- min(alpha_bounds$upper, 2.0e-04)
    alpha_init_new <- min(alpha_upper_new * 0.90, max(alpha_lower_new, 1.4e-04))
    tab_next <- set_param_natural(
      tab_next, "log10_alpha_o2",
      init = alpha_init_new,
      lower = alpha_lower_new,
      upper = alpha_upper_new,
      note_append = "autotune: capped lower so 4N burden is not damped too early"
    )
    cfg_next <- yaml_set(cfg_next, "alpha_o2_init", alpha_init_new)
    cfg_next <- yaml_set(cfg_next, "alpha_o2_min", alpha_lower_new)
    cfg_next <- yaml_set(cfg_next, "alpha_o2_max", alpha_upper_new)
    modifications <- c(modifications, sprintf("Restricted alpha_o2 to [%.4g, %.4g] with init %.4g so 4N growth damping stays weaker.", alpha_lower_new, alpha_upper_new, alpha_init_new))

    bsize_lower_new <- max(bsize_bounds$lower, 0.45)
    bsize_upper_new <- max(bsize_bounds$upper, 0.70)
    bsize_init_new <- min(bsize_upper_new, max(bp[["beta_size"]] * 1.25, 0.55))
    tab_next <- set_param_natural(
      tab_next, "beta_size",
      init = bsize_init_new,
      lower = bsize_lower_new,
      upper = bsize_upper_new,
      note_append = "autotune: increased so higher-ploidy cells contribute more burden volume"
    )
    modifications <- c(modifications, sprintf("Raised beta_size range to [%.4g, %.4g] with init %.4g to support 4N burden.", bsize_lower_new, bsize_upper_new, bsize_init_new))

    o2_init_lower_new <- max(o2init_bounds$lower, 2.2)
    o2_init_upper_new <- min(o2init_bounds$upper, 3.0)
    o2_init_init_new <- min(o2_init_upper_new, max(bp[["o2_init_pct"]] * 1.10, 2.5))
    tab_next <- set_param_natural(
      tab_next, "log10_o2_init_pct",
      init = o2_init_init_new,
      lower = o2_init_lower_new,
      upper = o2_init_upper_new,
      note_append = "autotune: shifted upward to preserve early oxygen and let 4N burden accumulate"
    )
    cfg_next <- yaml_set(cfg_next, "o2_init_pct_init", o2_init_init_new)
    cfg_next <- yaml_set(cfg_next, "o2_init_pct_min", o2_init_lower_new)
    cfg_next <- yaml_set(cfg_next, "o2_init_pct_max", o2_init_upper_new)
    modifications <- c(modifications, sprintf("Raised o2_init_pct to [%.4g, %.4g] with init %.4g to preserve early 4N growth.", o2_init_lower_new, o2_init_upper_new, o2_init_init_new))

    o2_rate_lower_new <- max(o2rate_bounds$lower, 0.25)
    o2_rate_upper_new <- min(o2rate_bounds$upper, 0.70)
    o2_rate_init_new <- min(o2_rate_upper_new, max(o2_rate_lower_new, min(bp[["o2_rate"]] * 0.75, 0.45)))
    tab_next <- set_param_natural(
      tab_next, "log10_o2_rate",
      init = o2_rate_init_new,
      lower = o2_rate_lower_new,
      upper = o2_rate_upper_new,
      note_append = "autotune: slowed oxygen decline so hypoxia-driven convergence is pushed later"
    )
    cfg_next <- yaml_set(cfg_next, "o2_rate_init", o2_rate_init_new)
    cfg_next <- yaml_set(cfg_next, "o2_rate_min", o2_rate_lower_new)
    cfg_next <- yaml_set(cfg_next, "o2_rate_max", o2_rate_upper_new)
    modifications <- c(modifications, sprintf("Reduced o2_rate to [%.4g, %.4g] with init %.4g so hypoxia-driven convergence happens later.", o2_rate_lower_new, o2_rate_upper_new, o2_rate_init_new))
  }

  if (is.finite(m$latest_mean_dead_buffer_frac) && m$latest_mean_dead_buffer_frac > 0.18) {
    reasons <- c(reasons, "Dead-buffer routing is still too large; nonviable daughter mass is dominating instead of redistributing viable daughters.")
    upper_new <- max(0.005, min(p_bounds$upper, bp[["p_misseg"]] * 1.20))
    init_new <- min(upper_new * 0.80, max(1e-4, bp[["p_misseg"]] * 0.85))
    tab_next <- set_param_natural(
      tab_next, "log10_p_misseg",
      init = init_new,
      upper = upper_new,
      note_append = "autotune: reduced after excessive dead-buffer routing"
    )
    modifications <- c(modifications, sprintf("Reduced p_misseg upper to %.4g and init to %.4g.", upper_new, init_new))

    bbuf_upper_new <- max(0.01, min(bbuf_bounds$upper, max(bp[["beta_buffer"]] * 0.80, 0.01)))
    bbuf_init_new <- min(bbuf_upper_new * 0.50, max(0.0, bp[["beta_buffer"]] * 0.70))
    tab_next <- set_param_natural(
      tab_next, "beta_buffer",
      init = bbuf_init_new,
      upper = bbuf_upper_new,
      note_append = "autotune: weakened background buffering after excessive dead-buffer routing"
    )
    modifications <- c(modifications, sprintf("Reduced beta_buffer upper to %.4g and init to %.4g.", bbuf_upper_new, bbuf_init_new))

    smax_lower_new <- min(0.997, max(smax_bounds$lower, bp[["smax"]] * 0.998, 0.995))
    smax_init_new <- min(0.999, max(bp[["smax"]], 0.997))
    tab_next <- set_param_natural(
      tab_next, "log10_smax",
      init = smax_init_new,
      lower = smax_lower_new,
      note_append = "autotune: kept extremely close to 1 to avoid extra background daughter death"
    )
    modifications <- c(modifications, sprintf("Raised smax lower to %.4g and init to %.4g.", smax_lower_new, smax_init_new))

    bloss_upper_new <- max(0.008, min(bloss_bounds$upper, max(bp[["beta_loss"]] * 0.90, 0.01)))
    bloss_init_new <- min(bloss_upper_new * 0.80, max(0.003, bp[["beta_loss"]] * 0.85))
    tab_next <- set_param_natural(
      tab_next, "log10_beta_loss",
      init = bloss_init_new,
      upper = bloss_upper_new,
      note_append = "autotune: lightened loss penalty after excessive dead-buffer routing"
    )
    modifications <- c(modifications, sprintf("Reduced beta_loss upper to %.4g and init to %.4g.", bloss_upper_new, bloss_init_new))
  }

  if (is.finite(m$burden_mean_log_err_2N) && is.finite(m$burden_mean_log_err_4N) &&
      m$burden_mean_log_err_2N > 0.25 && m$burden_mean_log_err_4N < 0.00) {
    reasons <- c(reasons, "2N burden is over-predicted while 4N burden remains under-supported; the burden mapping still underweights high-ploidy cells.")
    bsize_lower_new <- max(bsize_bounds$lower, 0.50)
    bsize_upper_new <- max(bsize_bounds$upper, 0.75)
    bsize_init_new <- min(bsize_upper_new, max(bp[["beta_size"]] * 1.15, 0.60))
    tab_next <- set_param_natural(
      tab_next, "beta_size",
      init = bsize_init_new,
      lower = bsize_lower_new,
      upper = bsize_upper_new,
      note_append = "autotune: increased after 2N-over / 4N-under burden split"
    )
    modifications <- c(modifications, sprintf("Raised beta_size range to [%.4g, %.4g] with init %.4g to shift more burden volume toward 4N states.", bsize_lower_new, bsize_upper_new, bsize_init_new))
  }

  if (is.finite(bp[["sigma_burden"]]) &&
      is.finite(sigma_bounds$upper) &&
      param_hit_upper(bp[["sigma_burden"]], sigma_bounds$upper, tol_frac = 0.05) &&
      is.finite(m$objective_burden) && m$objective_burden > 1.2) {
    reasons <- c(reasons, "sigma_burden is near its upper bound while burden error remains large; the optimizer still has too much room to ignore burden mismatch.")
    sigma_upper_new <- max(0.24, min(0.30, sigma_bounds$upper * 0.90))
    sigma_init_new <- min(sigma_upper_new * 0.90, max(sigma_bounds$lower, 0.22))
    tab_next <- set_param_natural(
      tab_next, "log10_sigma_burden",
      init = sigma_init_new,
      upper = sigma_upper_new,
      note_append = "autotune: tightened because burden mismatch stayed large while sigma_burden sat near its ceiling"
    )
    cfg_next <- yaml_set(cfg_next, "sigma_burden", sigma_init_new)
    cfg_next <- yaml_set(cfg_next, "sigma_burden_max", sigma_upper_new)
    modifications <- c(modifications, sprintf("Reduced sigma_burden upper to %.4g and init to %.4g so burden fit cannot be relaxed as easily.", sigma_upper_new, sigma_init_new))
  }

  if (is.finite(m$burden_mean_log_err_2N) && is.finite(m$burden_mean_log_err_4N) &&
      m$burden_mean_log_err_2N > 0.40 && m$burden_mean_log_err_4N > 0.40 &&
      is.finite(m$ploidy_mean_signed_err_4N) && m$ploidy_mean_signed_err_4N < 15 &&
      is.finite(m$latest_mean_live_frac) && m$latest_mean_live_frac > 0.70) {
    reasons <- c(reasons, "Both 2N and 4N burdens are over-predicted after 4N ploidy error has been reduced; hypoxia growth damping can be increased modestly.")
    alpha_upper_new <- min(3.0e-04, max(alpha_bounds$upper, bp[["alpha_o2"]] * 1.15))
    alpha_init_new <- min(alpha_upper_new * 0.90, max(bp[["alpha_o2"]] * 1.05, alpha_bounds$init))
    tab_next <- set_param_natural(
      tab_next, "log10_alpha_o2",
      init = alpha_init_new,
      upper = alpha_upper_new,
      note_append = "autotune: modestly increased only after both cohorts remained burden-high"
    )
    cfg_next <- yaml_set(cfg_next, "alpha_o2_init", alpha_init_new)
    cfg_next <- yaml_set(cfg_next, "alpha_o2_max", alpha_upper_new)
    modifications <- c(modifications, sprintf("Raised alpha_o2 upper to %.4g and init to %.4g only after both cohorts remained burden-high.", alpha_upper_new, alpha_init_new))
  }

  if (!length(modifications) && length(history_metrics) >= 2L) {
    reasons <- c(reasons, "No clear parameter-bound failure mode remains after the latest run.")
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

write_best_version <- function(best_version, best_data_obj, path) {
  writeLines(c(
    paste("best_version", best_version %||% "NA", sep = "\t"),
    paste("objective_data", if (is.finite(best_data_obj)) fmt_num(best_data_obj) else "NA", sep = "\t")
  ), path, useBytes = TRUE)
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  oxygen_root <- file.path(repo_root, "oxygen")

  base_config <- normalizePath(
    argv$base_config %||% file.path(oxygen_root, "config", "O2_NGLF_MAP_asymmetric_constrained_v2.yaml"),
    mustWork = TRUE
  )
  seed <- as_int(argv$seed, 1L)
  max_versions <- as_int(argv$max_versions, 5L)
  min_improve <- as_num(argv$min_improve, 0.10)
  runner <- normalizePath(
    argv$runner %||% file.path(oxygen_root, "code", "O2_NGLF_MAP_asymmetric", "run_fit_invivo_model_O2_NGLF_MAP_asymmetric.sh"),
    mustWork = TRUE
  )
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  session_root <- argv$session_root %||% file.path(oxygen_root, "results", paste0("O2_NGLF_MAP_asymmetric_autotune_", stamp))
  session_root <- normalizePath(session_root, mustWork = FALSE)
  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)

  base_lines <- read_yaml_lines(base_config)
  base_dir <- dirname(base_config)
  base_param_path <- resolve_path_local(yaml_get(base_lines, "parameter_table"), base_dir)
  if (is.null(base_param_path) || !file.exists(base_param_path)) stop("Could not resolve base parameter_table from ", base_config)
  data_dir_abs <- resolve_path_local(yaml_get(base_lines, "data_dir"), base_dir)
  if (is.null(data_dir_abs) || !dir.exists(data_dir_abs)) stop("Could not resolve data_dir from ", base_config)

  working_cfg <- base_lines
  working_tab <- read_parameter_table(base_param_path)
  history <- list()
  no_improve <- 0L
  best_version <- NA_character_
  best_data_obj <- Inf
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
    if (v > 1L) {
      prev_ckpt <- file.path(session_root, paste0("v", v - 1L), "seed1", "checkpoints", "best_params_transformed_latest.tsv")
      if (file.exists(prev_ckpt)) cmd <- c(cmd, paste0("--init_params_tsv=", prev_ckpt))
    }
    message("[autotune] Running ", version_label, " ...")
    status <- system2("bash", cmd)
    if (!identical(status, 0L)) stop("Autotune run failed for ", version_label, " with exit status ", status)
    if (!dir.exists(run_dir)) stop("Expected run directory missing: ", run_dir)

    analysis <- analyze_run(run_dir)
    metrics <- c(list(version = version_label), analysis$metrics)
    history[[length(history) + 1L]] <- metrics

    previous_metrics <- if (length(history) >= 2L) history[[length(history) - 1L]] else NULL
    proposal <- propose_next_version(working_cfg, working_tab, analysis, history)

    data_obj <- as_num(analysis$metrics$objective_data, Inf)
    if (is.finite(data_obj) && data_obj + min_improve < best_data_obj) {
      best_data_obj <- data_obj
      best_version <- version_label
      no_improve <- 0L
    } else if (is.finite(data_obj) && data_obj < Inf) {
      no_improve <- no_improve + 1L
    }

    write_history_tsv(history, history_path)
    write_best_version(best_version, best_data_obj, best_path)

    stop_reason <- NULL
    if (!length(proposal$modifications)) {
      stop_reason <- "No further parameter micro-adjustments were indicated by the current diagnostics."
    } else if (no_improve >= 2L) {
      stop_reason <- paste0("objective_data failed to improve by at least ", min_improve, " for two consecutive versions.")
    } else if (v >= max_versions) {
      stop_reason <- paste0("Reached max_versions=", max_versions, ".")
    }

    analysis_lines <- render_analysis(
      version_label = version_label,
      analysis = analysis,
      previous_metrics = previous_metrics,
      reasons = proposal$reasons,
      modifications = proposal$modifications,
      stop_reason = stop_reason
    )
    writeLines(analysis_lines, analysis_path, useBytes = TRUE)
    message("[autotune] Wrote analysis: ", analysis_path)

    if (!is.null(stop_reason)) break

    working_cfg <- proposal$cfg_lines
    working_tab <- proposal$param_table
  }

  write_history_tsv(history, history_path)
  write_best_version(best_version, best_data_obj, best_path)
  message("Autotune complete. Session root: ", session_root)
  message("Best version: ", best_version %||% "NA", " (objective_data=", fmt_num(best_data_obj), ")")
}

if (sys.nframe() == 0L) main()
