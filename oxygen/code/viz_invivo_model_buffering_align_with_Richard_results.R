#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

get_script_dir_self <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", farg[[1]])))
}

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

find_latest_fit_dir <- function(results_root) {
  dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  if (length(dirs) == 0) {
    stop("No fit result directories found under: ", results_root)
  }
  dirs[[which.max(file.info(dirs)$mtime)]]
}

normalize_cfg_for_viz <- function(cfg) {
  cfg$N_UNIT <- as.integer(cfg$N_UNIT %||% 22L)
  cfg$N_MIN <- as.integer(cfg$N_MIN %||% 22L)
  cfg$N_MAX <- as.integer(cfg$N_MAX %||% 154L)
  cfg$DT <- as.numeric(cfg$DT %||% 0.5)
  cfg$O2_fixed <- as.numeric(cfg$O2_fixed %||% 1.0)
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

read_run_params <- function(fit_dir) {
  p <- file.path(fit_dir, "best_params.tsv")
  if (!file.exists(p)) stop("Missing best_params.tsv in: ", fit_dir)
  tab <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain columns: parameter, value")
  }
  vals <- setNames(as.numeric(tab$value), as.character(tab$parameter))
  needed <- c("lam_min", "lam_max", "k_o", "p_misseg", "k_o_mis", "beta_buffer", "n_exp", "smax", "p_wgd")
  miss <- setdiff(needed, names(vals))
  if (length(miss) > 0) {
    stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
  }
  out <- as.list(vals[needed])
  out$alpha <- if ("alpha" %in% names(vals) && is.finite(vals[["alpha"]])) vals[["alpha"]] else 0
  out$gamma <- if ("gamma" %in% names(vals) && is.finite(vals[["gamma"]])) vals[["gamma"]] else 1
  out
}

simulate_one_full <- function(run_params, scenario, cfg, report_dt = 1.0) {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  grid_post <- cfg$N_MIN:cfg$N_MAX
  R0 <- length(grid_pre)
  R1 <- length(grid_post)

  init_state <- if (scenario$cohort == "2N") {
    make_init_state(
      grid_pre = grid_pre,
      grid_post = grid_post,
      ploidy = 2,
      layer = "pre",
      N_UNIT = cfg$N_UNIT,
      total_size = cfg$init_total_size
    )
  } else {
    make_init_state(
      grid_pre = grid_pre,
      grid_post = grid_post,
      ploidy = 4,
      layer = "post",
      N_UNIT = cfg$N_UNIT,
      total_size = cfg$init_total_size
    )
  }

  lambda0 <- growth_lambda(
    cfg$O2_fixed, grid_pre,
    lam_min = run_params$lam_min,
    lam_max = run_params$lam_max,
    k_o = run_params$k_o,
    R = run_params$lam_min,
    beta = run_params$beta_buffer,
    eta = 0,
    N_unit = cfg$N_UNIT
  )
  lambda1 <- growth_lambda(
    cfg$O2_fixed, grid_post,
    lam_min = run_params$lam_min,
    lam_max = run_params$lam_max,
    k_o = run_params$k_o,
    R = run_params$lam_min,
    beta = run_params$beta_buffer,
    eta = 0,
    N_unit = cfg$N_UNIT
  )
  p_mis <- if (exists(".pmisseg_of_O2", mode = "function", inherits = TRUE)) {
    as.numeric(.pmisseg_of_O2(cfg$O2_fixed, run_params))
  } else {
    k_o_mis_use <- max(as.numeric(run_params$k_o_mis), 1e-12)
    as.numeric(run_params$p_misseg) * (1 - (cfg$O2_fixed / (cfg$O2_fixed + k_o_mis_use)))
  }
  p_mis <- clip(p_mis, 0, 1)

  G <- .build_G_with_WGD(
    N0min = cfg$N_MIN, N0max = cfg$N_MAX, lambda0_vec = lambda0,
    p0_vec = p_mis, wgd_prob_vec = run_params$p_wgd,
    N1min = cfg$N_MIN, N1max = cfg$N_MAX, lambda1_vec = lambda1,
    p1_vec = p_mis,
    boundary = "drop",
    N_unit = cfg$N_UNIT,
    beta_buffer = run_params$beta_buffer,
    n_exp = run_params$n_exp,
    smax = run_params$smax
  )

  sim_end_step <- as.integer(round(scenario$sim_end_day / cfg$DT))
  obs_steps <- as.integer(round(scenario$obs_days / cfg$DT))
  report_steps <- as.integer(round(seq(0, scenario$sim_end_day, by = report_dt) / cfg$DT))
  keep_steps <- sort(unique(c(0L, sim_end_step, obs_steps, report_steps)))
  keep_steps <- keep_steps[keep_steps >= 0L & keep_steps <= sim_end_step]

  obs_map <- setNames(as.numeric(scenario$obs_burden), as.character(obs_steps))

  v <- as.numeric(init_state)
  I <- Diagonal(n = length(v))
  dose_scaled <- scenario$dose / cfg$dose_ref
  if (!is.finite(dose_scaled) || dose_scaled < 0) dose_scaled <- 0

  burden_rows <- vector("list", length(keep_steps))
  ploidy_rows <- vector("list", length(keep_steps))
  k <- 0L

  for (step in 0:sim_end_step) {
    if (step %in% keep_steps) {
      k <- k + 1L
      frac_pre <- v[seq_len(R0)]
      frac_post <- v[R0 + seq_len(R1)]
      frac_N <- frac_pre + frac_post
      if (sum(frac_N) > 0) {
        frac_N <- frac_N / sum(frac_N)
      } else {
        frac_N <- rep(1 / length(frac_N), length(frac_N))
      }

      burden_rows[[k]] <- data.frame(
        harvest = scenario$harvest,
        cohort = scenario$cohort,
        dose = scenario$dose,
        treat_day = scenario$treat_day,
        step = step,
        day = step * cfg$DT,
        pred_burden = sum(v),
        obs_burden = as.numeric(obs_map[as.character(step)]),
        row.names = NULL
      )
      ploidy_rows[[k]] <- data.frame(
        harvest = scenario$harvest,
        cohort = scenario$cohort,
        dose = scenario$dose,
        day = step * cfg$DT,
        N = as.integer(grid_pre),
        fraction = as.numeric(frac_N),
        row.names = NULL
      )
    }

    if (step >= sim_end_step) break

    t <- step * cfg$DT
    tx_mult <- 1.0
    if (isTRUE(cfg$fit_treatment)) {
      tx_mult <- if (t < scenario$treat_day || scenario$dose <= 0) {
        1.0
      } else {
        exp(-run_params$alpha * (dose_scaled^run_params$gamma))
      }
      tx_mult <- clip(tx_mult, cfg$tx_mult_min, 1.0)
    }

    Ntot <- sum(v)
    crowd <- if (cfg$crowding == "logistic") max(0, 1 - Ntot / cfg$K) else exp(-Ntot / cfg$K)
    A <- I + cfg$DT * (crowd * tx_mult * G)
    v <- as.numeric(A %*% v)
    v[!is.finite(v)] <- 0
    v[v < 0] <- 0
    if (sum(v) <= cfg$min_pop) break
  }

  list(
    burden = bind_rows(burden_rows),
    ploidy = bind_rows(ploidy_rows)
  )
}

simulate_one_full_horizon <- function(run_params, scenario, cfg, horizon_day, report_dt = 1.0) {
  sc <- scenario
  sc$sim_end_day <- as.numeric(max(horizon_day, 0))
  simulate_one_full(run_params, sc, cfg, report_dt = report_dt)
}

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

compute_ploidy_weighted_mean <- function(ploidy_all, cfg) {
  ploidy_all %>%
    group_by(harvest, cohort, dose, day) %>%
    summarise(
      weighted_mean_N = sum(N * fraction, na.rm = TRUE) / pmax(sum(fraction, na.rm = TRUE), 1e-12),
      .groups = "drop"
    ) %>%
    mutate(weighted_mean_ploidy = weighted_mean_N / cfg$N_UNIT)
}

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

find_fit_dirs_under <- function(root_dir) {
  all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  sub_dirs <- all_dirs[all_dirs != root_dir]

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

  run_params <- read_run_params(fit_dir)
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
      title = "Richard Model: In Vivo Burden Trajectory (Normalized)",
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
      title = "Richard Model: In Vivo Burden Trajectory (Absolute)",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Day",
      y = "Predicted burden / observed burden (raw units)"
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
  burden_all_real <- burden_all %>%
    mutate(
      pred_burden_cell_number = as.numeric(pred_burden),
      obs_burden_cell_number_low = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_min, NA_real_),
      obs_burden_cell_number_mid = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_mid, NA_real_),
      obs_burden_cell_number_high = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_max, NA_real_)
    )

  p_burden_abs_real <- ggplot(burden_all_real, aes(x = day, y = pred_burden_cell_number)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_ribbon(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_low) & !is.na(obs_burden_cell_number_high)),
      aes(ymin = obs_burden_cell_number_low, ymax = obs_burden_cell_number_high),
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
      title = "Richard Model: In Vivo Burden Trajectory (Absolute, Real Scale)",
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
      title = "Richard Model: Predicted Ploidy Distribution Over Time",
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
      title = paste0("Richard Model: Ploidy Over Time (Top ", top_n, " N States)"),
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
      title = "Richard Model: Weighted Mean Ploidy Over Time",
      subtitle = "Weighted by predicted ploidy fractions",
      x = "Day",
      y = "Weighted Mean Ploidy (P = N / N_UNIT)"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "burden_trend.pdf"), p_burden, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute.pdf"), p_burden_abs, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute(real_scale).pdf"), p_burden_abs_real, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_heatmap_over_time.pdf"), p_ploidy_heatmap, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_top_states_over_time.pdf"), p_ploidy_lines, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_weighted_mean_over_time.pdf"), p_ploidy_weighted_mean, width = 13, height = 9)

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

main <- function() {
  script_dir <- get_script_dir_self()
  source(file.path(script_dir, "fit_invivo_model_buffering_align_with_Richard.R"))
  source(file.path(script_dir, "model_buffering_align_with_Richard.R"))

  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  results_root <- normalizePath(file.path(script_dir, "..", "results"), mustWork = FALSE)
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
    normalizePath(file.path(script_dir, "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
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
