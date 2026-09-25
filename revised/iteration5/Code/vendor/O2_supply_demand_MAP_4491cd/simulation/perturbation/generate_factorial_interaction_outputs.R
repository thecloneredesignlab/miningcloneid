#!/usr/bin/env Rscript

# Generate factorial post-fit interaction trajectories only.  Endpoint analysis
# and figures are produced by separate consumers.

suppressPackageStartupMessages(library(dplyr))

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_perturbation_utils.R"), local = environment())
rm(.script_dir)
load_o2sd_perturbation_model(environment())

make_interaction_design <- function(run_params,
                                    o2_values,
                                    pmis_base_values,
                                    p_wgd_values,
                                    trigger_burden_values,
                                    initial_burden_cells,
                                    experiment_name = "triggered_drug",
                                    scenario_prefix = "drug") {
  baseline_p_mis_base <- as.numeric(run_params$p_mis_base)
  baseline_p_wgd <- as.numeric(run_params$p_wgd)
  expand.grid(
    o2_S0 = as.numeric(o2_values),
    p_wgd_post = as.numeric(p_wgd_values),
    trigger_burden_cells = as.numeric(trigger_burden_values),
    p_mis_base_post = as.numeric(pmis_base_values),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      experiment = experiment_name,
      scenario_id = safe_id(
        scenario_prefix, "o2", num_tag(o2_S0), "pwgd", num_tag(p_wgd_post),
        "trigger", num_tag(trigger_burden_cells), "pmiss", num_tag(p_mis_base_post)
      ),
      initial_burden_cells = as.numeric(initial_burden_cells),
      p_mis_base_pre = baseline_p_mis_base,
      p_wgd_pre = baseline_p_wgd,
      p_wgd = p_wgd_post,
      varied_parameter = "factorial",
      varied_value = NA_real_, status = "planned", trigger_day = NA_real_,
      actual_trigger_burden_cells = NA_real_, post_treatment_duration_day = NA_real_,
      sim_end_day = NA_real_
    ) %>%
    select(
      experiment, scenario_id, varied_parameter, varied_value,
      initial_burden_cells, trigger_burden_cells, trigger_day,
      actual_trigger_burden_cells, o2_S0, p_mis_base_pre,
      p_mis_base_post, p_wgd_pre, p_wgd_post, p_wgd,
      post_treatment_duration_day, sim_end_day, status
    )
}

make_untreated_design <- function(run_params,
                                  o2_values,
                                  pmis_base_values,
                                  p_wgd_values,
                                  initial_burden_cells,
                                  horizon_day) {
  expand.grid(
    o2_S0 = as.numeric(o2_values),
    p_wgd_post = as.numeric(p_wgd_values),
    p_mis_base_post = as.numeric(pmis_base_values),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      experiment = "untreated_factorial",
      scenario_id = safe_id(
        "untreated", "o2", num_tag(o2_S0), "pwgd", num_tag(p_wgd_post),
        "pmiss", num_tag(p_mis_base_post)
      ),
      initial_burden_cells = as.numeric(initial_burden_cells),
      trigger_burden_cells = NA_real_, trigger_day = NA_real_,
      actual_trigger_burden_cells = NA_real_,
      p_mis_base_pre = as.numeric(run_params$p_mis_base),
      p_wgd_pre = as.numeric(run_params$p_wgd), p_wgd = p_wgd_post,
      varied_parameter = "factorial", varied_value = NA_real_, status = "untreated",
      post_treatment_duration_day = NA_real_, sim_end_day = as.numeric(horizon_day)
    ) %>%
    select(
      experiment, scenario_id, varied_parameter, varied_value,
      initial_burden_cells, trigger_burden_cells, trigger_day,
      actual_trigger_burden_cells, o2_S0, p_mis_base_pre,
      p_mis_base_post, p_wgd_pre, p_wgd_post, p_wgd,
      post_treatment_duration_day, sim_end_day, status
    )
}

make_exact_ploidy_mixture_state <- function(cfg, total_live_cells, ploidy_values, fractions) {
  ploidy_values <- as.numeric(ploidy_values)
  fractions <- as.numeric(fractions)
  if (length(ploidy_values) != length(fractions)) stop("ploidy_values and fractions must have the same length.")
  if (!length(ploidy_values) || any(!is.finite(ploidy_values)) || any(!is.finite(fractions))) {
    stop("Initial ploidy mixture values must be finite.")
  }
  if (any(fractions < 0) || sum(fractions) <= 0) {
    stop("Initial ploidy mixture fractions must be non-negative and have positive total mass.")
  }
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  target_N <- as.integer(round(ploidy_values * as.numeric(cfg$N_UNIT)))
  if (any(target_N < cfg$N_MIN | target_N > cfg$N_MAX)) stop("Initial ploidy mixture maps outside cfg N range.")
  fractions <- fractions / sum(fractions)
  x <- rep(0, length(grid_pre)); names(x) <- as.character(grid_pre)
  for (i in seq_along(target_N)) {
    x[as.character(target_N[[i]])] <- x[as.character(target_N[[i]])] + as.numeric(total_live_cells) * fractions[[i]]
  }
  x
}

make_initial_state_for_mode <- function(cfg, total_live_cells, initial_state_mode,
                                        init_mean, init_sd, init_min, init_max) {
  initial_state_mode <- tolower(as.character(initial_state_mode))
  if (initial_state_mode %in% c("2n4n", "2n_4n", "two_four_equal", "two_four_50_50")) {
    return(make_exact_ploidy_mixture_state(cfg, total_live_cells, c(2, 4), c(0.5, 0.5)))
  }
  if (initial_state_mode %in% c("continuous", "normal", "truncated_normal")) {
    return(make_continuous_ploidy_state(cfg, total_live_cells, init_mean, init_sd, init_min, init_max))
  }
  stop("Unknown initial_state_mode=", initial_state_mode, ". Use continuous or 2n4n.")
}

summarise_factorial_ploidy <- function(sim, cfg, live_min_cells = 1e-6) {
  live_state <- sim$live_state
  live_pop <- rowSums(live_state, na.rm = TRUE)
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  ploidy_grid <- as.numeric(grid_pre) / as.numeric(cfg$N_UNIT)
  metrics <- matrix(NA_real_, nrow = nrow(live_state), ncol = 7L)
  colnames(metrics) <- c("mean_N", "mean_ploidy", "frac_ploidy_lt2", "frac_ploidy_2to3",
                         "frac_ploidy_3to4", "frac_ploidy_ge4", "frac_ploidy_ge5")
  for (i in which(is.finite(live_pop) & live_pop > live_min_cells)) {
    frac <- as.numeric(live_state[i, ]) / live_pop[[i]]
    s <- sum(frac, na.rm = TRUE)
    if (is.finite(s) && s > 1e-9) {
      metrics[i, ] <- c(
        sum(as.numeric(grid_pre) * frac, na.rm = TRUE) / s,
        sum(ploidy_grid * frac, na.rm = TRUE) / s,
        sum(frac[ploidy_grid < 2], na.rm = TRUE),
        sum(frac[ploidy_grid >= 2 & ploidy_grid < 3], na.rm = TRUE),
        sum(frac[ploidy_grid >= 3 & ploidy_grid < 4], na.rm = TRUE),
        sum(frac[ploidy_grid >= 4], na.rm = TRUE),
        sum(frac[ploidy_grid >= 5], na.rm = TRUE)
      )
    }
  }
  data.frame(
    scenario_id = sim$burden$scenario_id,
    experiment = sim$burden$experiment,
    segment = sim$burden$segment,
    local_day = sim$burden$local_day,
    day = sim$burden$day,
    live_cells = live_pop,
    as.data.frame(metrics, check.names = FALSE),
    stringsAsFactors = FALSE
  )
}

simulate_factorial_trajectory <- function(init_state, run_params, cfg, horizon_day, report_dt,
                                          init_dead_hypoxia_state = NULL,
                                          init_dead_buffer_state = NULL,
                                          start_day = 0, scenario_id = "scenario",
                                          experiment = "factorial_interaction",
                                          segment = "single", live_min_cells = 1e-6) {
  sim <- simulate_cpp_trajectory(
    init_state, run_params, cfg, horizon_day, report_dt,
    init_dead_hypoxia_state, init_dead_buffer_state, start_day,
    scenario_id, experiment, segment
  )
  keep_burden <- c(
    "scenario_id", "experiment", "segment", "local_day", "day", "step",
    "p_mis_base", "p_wgd", "o2_S0", "pred_burden_cells",
    "pred_burden_live_cells", "pred_burden_dead_hypoxia_cells",
    "pred_burden_dead_buffer_cells", "pred_burden_dead_total_cells",
    "pred_burden_volume_mm3", "pred_burden_live_volume_mm3",
    "pred_burden_dead_total_volume_mm3", "pred_o2_target_pct", "pred_o2_pct",
    "pred_log10_burden_cells", "pred_log10_burden_live_cells",
    "pred_log10_burden_dead_total_cells", "pred_log10_burden_volume_mm3",
    "pred_o2_lag_gap_pct"
  )
  sim$burden <- sim$burden[, keep_burden, drop = FALSE]
  sim$ploidy_summary <- summarise_factorial_ploidy(sim, cfg, live_min_cells)
  sim$ploidy <- NULL
  sim
}

annotate_factorial_timecourses <- function(sim, design_row) {
  sim$burden$scenario_id <- as.character(design_row$scenario_id)
  sim$burden$experiment <- as.character(design_row$experiment)
  sim$ploidy_summary$scenario_id <- as.character(design_row$scenario_id)
  sim$ploidy_summary$experiment <- as.character(design_row$experiment)
  add_cols <- setdiff(names(design_row), names(sim$burden))
  burden <- bind_cols(as.data.frame(design_row[rep(1, nrow(sim$burden)), add_cols, drop = FALSE]), sim$burden)
  add_cols_ploidy <- setdiff(names(design_row), names(sim$ploidy_summary))
  ploidy <- bind_cols(as.data.frame(design_row[rep(1, nrow(sim$ploidy_summary)), add_cols_ploidy, drop = FALSE]), sim$ploidy_summary)
  list(burden = burden, ploidy_summary = ploidy)
}

simulate_intermittent_post <- function(init_state,
                                       init_dead_hypoxia_state,
                                       init_dead_buffer_state,
                                       run_params_on,
                                       run_params_off,
                                       cfg,
                                       start_day,
                                       followup_day,
                                       on_day = 7,
                                       off_day = 21,
                                       report_dt = 1,
                                       scenario_id,
                                       experiment,
                                       live_min_cells = 1e-6) {
  current_live <- as.numeric(init_state)
  current_dead_hypoxia <- as.numeric(init_dead_hypoxia_state)
  current_dead_buffer <- as.numeric(init_dead_buffer_state)
  current_day <- as.numeric(start_day)
  remaining <- as.numeric(followup_day)
  cycle <- 1L
  segment_index <- 1L
  burden_rows <- list()
  ploidy_rows <- list()

  run_segment <- function(rp, duration, segment_name, treatment_on, cycle_index, keep_start) {
    seg <- simulate_factorial_trajectory(
      init_state = current_live,
      init_dead_hypoxia_state = current_dead_hypoxia,
      init_dead_buffer_state = current_dead_buffer,
      run_params = rp, cfg = cfg, horizon_day = duration, report_dt = report_dt,
      start_day = current_day, scenario_id = scenario_id, experiment = experiment,
      segment = segment_name, live_min_cells = live_min_cells
    )
    current_live <<- as.numeric(seg$live_state[nrow(seg$live_state), ])
    current_dead_hypoxia <<- as.numeric(seg$dead_hypoxia_state[nrow(seg$dead_hypoxia_state), ])
    current_dead_buffer <<- as.numeric(seg$dead_buffer_state[nrow(seg$dead_buffer_state), ])
    current_day <<- current_day + as.numeric(duration)
    keep <- rep(TRUE, nrow(seg$burden))
    if (!keep_start) keep <- seg$burden$local_day > 1e-9
    seg$burden <- seg$burden[keep, , drop = FALSE]
    seg$burden$cycle <- cycle_index
    seg$burden$treatment_on <- treatment_on
    keep_p <- rep(TRUE, nrow(seg$ploidy_summary))
    if (!keep_start) keep_p <- seg$ploidy_summary$local_day > 1e-9
    seg$ploidy_summary <- seg$ploidy_summary[keep_p, , drop = FALSE]
    seg$ploidy_summary$cycle <- cycle_index
    seg$ploidy_summary$treatment_on <- treatment_on
    seg
  }

  while (remaining > 1e-9) {
    duration_on <- min(as.numeric(on_day), remaining)
    if (duration_on > 1e-9) {
      seg <- run_segment(run_params_on, duration_on, "drug_on", TRUE, cycle, segment_index == 1L)
      burden_rows[[segment_index]] <- seg$burden
      ploidy_rows[[segment_index]] <- seg$ploidy_summary
      segment_index <- segment_index + 1L
      remaining <- remaining - duration_on
    }
    duration_off <- min(as.numeric(off_day), remaining)
    if (duration_off > 1e-9) {
      seg <- run_segment(run_params_off, duration_off, "drug_off", FALSE, cycle, FALSE)
      burden_rows[[segment_index]] <- seg$burden
      ploidy_rows[[segment_index]] <- seg$ploidy_summary
      segment_index <- segment_index + 1L
      remaining <- remaining - duration_off
    }
    cycle <- cycle + 1L
  }
  list(
    burden = bind_rows(burden_rows), ploidy_summary = bind_rows(ploidy_rows),
    live_state = current_live, dead_hypoxia_state = current_dead_hypoxia,
    dead_buffer_state = current_dead_buffer
  )
}

run_factorial_interaction_simulation <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- resolve_path_value(argv$fit_dir %||% argv$run_dir, getwd())
  if (is.null(fit_dir)) stop("Missing required argument: --fit_dir=/path/to/seed_dir")
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  cfg <- prepare_sim_cfg(read_fit_config(fit_dir), argv)
  run_params <- read_run_params(fit_dir, cfg)
  mode <- normalize_interaction_mode(argv$mode %||% "triggered_drug")

  out_dir <- resolve_path_value(argv$simulation_dir %||% argv$out_dir, getwd())
  if (is.null(out_dir)) {
    out_dir <- file.path(fit_dir, "simulation", "interactions", interaction_output_stem(mode))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  default_initial_burden_cells <- if (identical(mode, "untreated_factorial")) 1e6 else 1e5
  initial_burden_cells <- as_num(argv$initial_burden_cells, as_num(argv$initial_live_cells, default_initial_burden_cells))
  o2_values <- parse_num_list(argv$o2_values, c(0.5, 1, 2, 3, 4, 5))
  sweep_ratios <- c(1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 10, 50, 100)
  pmis_base_values <- parse_num_list(argv$pmis_base_values, 5e-3 * sweep_ratios)
  p_wgd_values <- parse_num_list(argv$p_wgd_values, as.numeric(run_params$p_wgd) * sweep_ratios)
  trigger_burden_values <- parse_num_list(argv$trigger_burden_values, c(5e5, 1e6, 2e6, 5e6, 1e7))
  if (as_bool(argv$smoke, FALSE)) {
    o2_values <- head(o2_values, 2); pmis_base_values <- head(pmis_base_values, 2)
    p_wgd_values <- head(p_wgd_values, 2); trigger_burden_values <- head(trigger_burden_values, 2)
  }

  pre_horizon_day <- as_num(argv$pre_horizon_day, 1000)
  post_treatment_day <- as_num(argv$post_treatment_day, 1000)
  untreated_horizon_day <- as_num(argv$horizon_day, if (identical(mode, "untreated_factorial")) 100 else 1000)
  intermittent_followup_day <- as_num(argv$intermittent_followup_day, as_num(argv$intermittent_followup_weeks, 156) * 7)
  intermittent_on_day <- as_num(argv$intermittent_on_day, as_num(argv$drug_on_day, 7))
  intermittent_off_day <- as_num(argv$intermittent_off_day, as_num(argv$drug_off_day, 21))
  report_dt <- as_num(argv$report_dt, 1.0)
  trigger_check_dt <- as_num(argv$trigger_check_dt, report_dt)
  live_min_cells <- as_num(argv$live_min_cells, 1e-6)
  init_mean <- as_num(argv$initial_ploidy_mean, 3.0)
  init_sd <- as_num(argv$initial_ploidy_sd, 0.4)
  init_min <- as_num(argv$initial_ploidy_min, 1.5)
  init_max <- as_num(argv$initial_ploidy_max, 6.0)
  initial_state_mode <- as.character(argv$initial_state_mode %||%
    if (identical(mode, "untreated_factorial")) "2n4n" else "continuous")

  design <- if (mode %in% c("triggered_drug", "intermittent_drug")) {
    make_interaction_design(
      run_params, o2_values, pmis_base_values, p_wgd_values, trigger_burden_values,
      initial_burden_cells, mode, if (identical(mode, "intermittent_drug")) "intermittent" else "drug"
    )
  } else {
    make_untreated_design(run_params, o2_values, pmis_base_values, p_wgd_values,
                          initial_burden_cells, untreated_horizon_day)
  }
  design$initial_ploidy_mean <- init_mean
  design$initial_ploidy_sd <- init_sd
  design$initial_ploidy_min <- init_min
  design$initial_ploidy_max <- init_max
  design$initial_state_mode <- initial_state_mode
  design$fit_dir <- fit_dir

  message("Running ", mode, " design: ", nrow(design), " scenarios. Output: ", out_dir)
  init_state <- make_initial_state_for_mode(
    cfg, initial_burden_cells, initial_state_mode, init_mean, init_sd, init_min, init_max
  )
  burden_rows <- vector("list", nrow(design))
  ploidy_rows <- vector("list", nrow(design))
  design_rows <- vector("list", nrow(design))

  if (identical(mode, "untreated_factorial")) {
    for (i in seq_len(nrow(design))) {
      dr <- design[i, , drop = FALSE]
      message("[", i, "/", nrow(design), "] ", dr$scenario_id)
      rp <- run_params; rp$o2_S0 <- as.numeric(dr$o2_S0)
      rp$p_wgd <- as.numeric(dr$p_wgd_post); rp$p_mis_base <- as.numeric(dr$p_mis_base_post)
      sim <- simulate_factorial_trajectory(
        init_state, rp, cfg, untreated_horizon_day, report_dt,
        scenario_id = dr$scenario_id, experiment = dr$experiment,
        segment = "untreated", live_min_cells = live_min_cells
      )
      ann <- annotate_factorial_timecourses(sim, dr)
      burden_rows[[i]] <- ann$burden; ploidy_rows[[i]] <- ann$ploidy_summary; design_rows[[i]] <- dr
    }
  } else {
    pre_cache <- list()
    pre_grid <- unique(design[, "o2_S0", drop = FALSE])
    for (i in seq_len(nrow(pre_grid))) {
      key <- format_value_label(pre_grid$o2_S0[[i]])
      rp <- run_params; rp$o2_S0 <- as.numeric(pre_grid$o2_S0[[i]])
      rp$p_wgd <- as.numeric(run_params$p_wgd); rp$p_mis_base <- as.numeric(run_params$p_mis_base)
      message("[pre ", i, "/", nrow(pre_grid), "] O2=", signif(rp$o2_S0, 6))
      pre_cache[[key]] <- simulate_factorial_trajectory(
        init_state, rp, cfg, pre_horizon_day, trigger_check_dt,
        scenario_id = paste0("trajectory_", i), experiment = "factorial_interaction_cache",
        segment = "pre", live_min_cells = live_min_cells
      )
    }
    for (i in seq_len(nrow(design))) {
      dr <- design[i, , drop = FALSE]
      pre <- pre_cache[[format_value_label(dr$o2_S0)]]
      message("[", i, "/", nrow(design), "] ", dr$scenario_id)
      hit_idx <- which(pre$burden$pred_burden_cells >= as.numeric(dr$trigger_burden_cells))
      if (!length(hit_idx)) {
        dr$status <- "trigger_not_reached"; dr$trigger_day <- NA_real_
        dr$actual_trigger_burden_cells <- NA_real_; dr$post_treatment_duration_day <- NA_real_
        dr$sim_end_day <- pre_horizon_day
        ann <- annotate_factorial_timecourses(pre, dr)
        burden_rows[[i]] <- ann$burden; ploidy_rows[[i]] <- ann$ploidy_summary; design_rows[[i]] <- dr
        next
      }
      hit <- hit_idx[[1]]
      trigger_day <- as.numeric(pre$burden$day[[hit]])
      dr$status <- if (identical(mode, "intermittent_drug")) "intermittent_treated" else "treated"
      dr$trigger_day <- trigger_day
      dr$actual_trigger_burden_cells <- as.numeric(pre$burden$pred_burden_cells[[hit]])
      followup_day <- if (identical(mode, "intermittent_drug")) intermittent_followup_day else post_treatment_day
      dr$post_treatment_duration_day <- followup_day; dr$sim_end_day <- trigger_day + followup_day
      if (identical(mode, "intermittent_drug")) {
        dr$intermittent_on_day <- intermittent_on_day; dr$intermittent_off_day <- intermittent_off_day
      }
      rp_post <- run_params; rp_post$o2_S0 <- as.numeric(dr$o2_S0)
      rp_post$p_wgd <- as.numeric(dr$p_wgd_post); rp_post$p_mis_base <- as.numeric(dr$p_mis_base_post)
      if (identical(mode, "intermittent_drug")) {
        rp_off <- run_params; rp_off$o2_S0 <- as.numeric(dr$o2_S0)
        rp_off$p_wgd <- as.numeric(run_params$p_wgd); rp_off$p_mis_base <- as.numeric(run_params$p_mis_base)
        post <- simulate_intermittent_post(
          as.numeric(pre$live_state[hit, ]), as.numeric(pre$dead_hypoxia_state[hit, ]),
          as.numeric(pre$dead_buffer_state[hit, ]), rp_post, rp_off, cfg, trigger_day,
          intermittent_followup_day, intermittent_on_day, intermittent_off_day, report_dt,
          dr$scenario_id, dr$experiment, live_min_cells
        )
      } else {
        post <- simulate_factorial_trajectory(
          as.numeric(pre$live_state[hit, ]), rp_post, cfg, post_treatment_day, report_dt,
          as.numeric(pre$dead_hypoxia_state[hit, ]), as.numeric(pre$dead_buffer_state[hit, ]),
          trigger_day, dr$scenario_id, dr$experiment, "post", live_min_cells
        )
      }
      pre_keep <- pre
      pre_keep$burden <- pre_keep$burden %>% filter(day < trigger_day - 1e-9)
      pre_keep$ploidy_summary <- pre_keep$ploidy_summary %>% filter(day < trigger_day - 1e-9)
      pre_ann <- annotate_factorial_timecourses(pre_keep, dr)
      post_ann <- annotate_factorial_timecourses(post, dr)
      burden_rows[[i]] <- bind_rows(pre_ann$burden, post_ann$burden)
      ploidy_rows[[i]] <- bind_rows(pre_ann$ploidy_summary, post_ann$ploidy_summary)
      design_rows[[i]] <- dr
    }
  }

  design_out <- bind_rows(design_rows)
  burden_all <- bind_rows(burden_rows)
  ploidy_summary <- bind_rows(ploidy_rows)
  write_tsv(design_out, file.path(out_dir, "interaction_design.tsv"))
  write_tsv(burden_all, file.path(out_dir, "burden_timecourse.tsv"))
  write_tsv(ploidy_summary, file.path(out_dir, "ploidy_summary_timecourse.tsv"))
  write_artifact_manifest(
    out_dir,
    list(
      artifact_manifest_row("design", "interaction_design.tsv", "simulation_design", design_out),
      artifact_manifest_row("burden", "burden_timecourse.tsv", "simulated_trajectory", burden_all),
      artifact_manifest_row("ploidy", "ploidy_summary_timecourse.tsv", "simulated_trajectory", ploidy_summary)
    ),
    "simulation_manifest.tsv"
  )
  message("Done. Wrote simulation outputs to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_factorial_interaction_simulation()
}
