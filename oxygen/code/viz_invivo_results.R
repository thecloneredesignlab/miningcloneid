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
  needed <- c("R", "beta", "eta", "pwgd", "mr_lethality0", "mr_lethality1", "pmis_O2_1", "pmis_O2_0", "alpha", "gamma")
  miss <- setdiff(needed, names(vals))
  if (length(miss) > 0) {
    stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
  }
  as.list(vals[needed])
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

  lambda0 <- growth_lambda(cfg$O2_fixed, grid_pre, R = run_params$R, beta = run_params$beta, eta = run_params$eta, N_unit = cfg$N_UNIT)
  lambda1 <- growth_lambda(cfg$O2_fixed, grid_post, R = run_params$R, beta = run_params$beta, eta = run_params$eta, N_unit = cfg$N_UNIT)
  logp <- (1 - cfg$O2_fixed) * log10(run_params$pmis_O2_0) + cfg$O2_fixed * log10(run_params$pmis_O2_1)
  p_mis <- 10^logp

  G <- .build_G_with_WGD(
    N0min = cfg$N_MIN, N0max = cfg$N_MAX, lambda0_vec = lambda0,
    p0_vec = p_mis, wgd_prob_vec = run_params$pwgd,
    N1min = cfg$N_MIN, N1max = cfg$N_MAX, lambda1_vec = lambda1,
    p1_vec = p_mis, mr_lethality0 = run_params$mr_lethality0,
    mr_lethality1 = run_params$mr_lethality1
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

main <- function() {
  script_dir <- get_script_dir_self()
  source(file.path(script_dir, "fit_invivo_ploidy_buffer.R"))
  source(file.path(script_dir, "model_functions_ploidy_buffer.R"))

  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  results_root <- normalizePath(file.path(script_dir, "..", "results"), mustWork = FALSE)
  fit_dir <- if (!is.null(argv$fit_dir)) {
    normalizePath(argv$fit_dir, mustWork = TRUE)
  } else {
    normalizePath(find_latest_fit_dir(results_root), mustWork = TRUE)
  }

  out_dir <- if (!is.null(argv$out_dir)) {
    argv$out_dir
  } else {
    file.path(fit_dir, "viz")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cfg_path <- file.path(fit_dir, "fit_config.rds")
  if (!file.exists(cfg_path)) stop("Missing fit_config.rds in: ", fit_dir)
  cfg <- readRDS(cfg_path)
  cfg <- normalize_cfg_for_viz(cfg)

  if (!is.null(argv$dose_zero_only)) cfg$dose_zero_only <- as_bool(argv$dose_zero_only, cfg$dose_zero_only)
  if (!is.null(argv$truncate_at_treatment)) cfg$truncate_at_treatment <- as_bool(argv$truncate_at_treatment, cfg$truncate_at_treatment)
  if (!is.null(argv$ploidy_at_harvest)) cfg$ploidy_at_harvest <- as_bool(argv$ploidy_at_harvest, cfg$ploidy_at_harvest)
  if (!is.null(argv$max_scenarios)) cfg$max_scenarios <- as_num(argv$max_scenarios, cfg$max_scenarios)

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

  run_params <- read_run_params(fit_dir)
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)

  sim_list <- lapply(scenarios, function(sc) simulate_one_full(run_params, sc, cfg, report_dt = report_dt))
  burden_all <- bind_rows(lapply(sim_list, `[[`, "burden"))
  ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy"))

  if (nrow(burden_all) == 0 || nrow(ploidy_all) == 0) {
    stop("No simulation output generated; check fit/data configuration.")
  }

  write.table(burden_all, file = file.path(out_dir, "burden_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_all, file = file.path(out_dir, "ploidy_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  ploidy_mean <- ploidy_all %>%
    group_by(harvest, cohort, dose, day) %>%
    summarise(
      weighted_mean_N = sum(N * fraction, na.rm = TRUE) / pmax(sum(fraction, na.rm = TRUE), 1e-12),
      .groups = "drop"
    ) %>%
    mutate(weighted_mean_ploidy = weighted_mean_N / cfg$N_UNIT)
  write.table(ploidy_mean, file = file.path(out_dir, "ploidy_weighted_mean_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  p_burden <- ggplot(burden_all, aes(x = day, y = pred_burden)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_point(
      data = burden_all %>% filter(!is.na(obs_burden)),
      aes(y = obs_burden),
      color = "black",
      size = 1
    ) +
    facet_wrap(~ harvest, scales = "free_y", ncol = 2) +
    labs(
      title = "In Vivo Burden Trajectory: Predicted vs Observed",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Day",
      y = "Tumor Volume"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_heatmap <- ggplot(ploidy_all, aes(x = day, y = N, fill = fraction)) +
    geom_tile() +
    facet_wrap(~ harvest, ncol = 2) +
    scale_fill_viridis_c(option = "C", trans = "sqrt") +
    labs(
      title = "Predicted Ploidy Distribution Over Time",
      subtitle = "Heatmap of fraction by chromosome number (N)",
      x = "Day",
      y = "Chromosome Number (N)",
      fill = "Fraction"
    ) +
    theme_bw(base_size = 11)

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
      title = paste0("Predicted Ploidy Over Time (Top ", top_n, " N States)"),
      x = "Day",
      y = "Fraction",
      color = "N"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_weighted_mean <- ggplot(ploidy_mean, aes(x = day, y = weighted_mean_ploidy)) +
    geom_line(color = "#d62728", linewidth = 0.9) +
    facet_wrap(~ harvest, ncol = 2) +
    labs(
      title = "Weighted Mean Ploidy Over Time",
      subtitle = "Weighted by predicted ploidy fractions",
      x = "Day",
      y = "Weighted Mean Ploidy (P = N / N_UNIT)"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "burden_trend.png"), p_burden, width = 13, height = 9, dpi = 150)
  ggsave(file.path(out_dir, "ploidy_heatmap_over_time.png"), p_ploidy_heatmap, width = 13, height = 9, dpi = 150)
  ggsave(file.path(out_dir, "ploidy_top_states_over_time.png"), p_ploidy_lines, width = 13, height = 9, dpi = 150)
  ggsave(file.path(out_dir, "ploidy_weighted_mean_over_time.png"), p_ploidy_weighted_mean, width = 13, height = 9, dpi = 150)

  message("Done. Visualization outputs written to: ", normalizePath(out_dir))
}

if (sys.nframe() == 0) {
  main()
}
