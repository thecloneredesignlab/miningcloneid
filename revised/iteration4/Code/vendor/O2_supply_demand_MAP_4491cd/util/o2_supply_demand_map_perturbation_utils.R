#!/usr/bin/env Rscript

# Shared post-fit perturbation simulation utilities.  This file defines no
# executable workflow: simulation, analysis, and visualization entrypoints live
# in their respective functional folders.

suppressPackageStartupMessages(library(dplyr))

.o2_sim_bootstrap_workflow_root <- local({
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
  file_args <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  )
  starts <- unique(c(dirname(frame_files), dirname(file_args), getwd()))
  for (start in starts[nzchar(starts)]) {
    current <- normalizePath(start, mustWork = FALSE)
    for (depth in 0:10) {
      candidates <- c(
        current,
        file.path(current, "oxygen", "code", "O2_supply_demand_MAP"),
        file.path(current, "code", "O2_supply_demand_MAP")
      )
      hits <- candidates[file.exists(file.path(
        candidates,
        "util",
        "o2_supply_demand_map_shared.R"
      ))]
      if (length(hits)) return(normalizePath(hits[[1L]], mustWork = TRUE))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate the O2_supply_demand_MAP workflow root.")
})

WORKFLOW_ROOT <- .o2_sim_bootstrap_workflow_root
SCRIPT_DIR <- file.path(WORKFLOW_ROOT, "util")
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())

rm(.o2_sim_bootstrap_workflow_root)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_bool <- o2sd_as_bool

load_o2sd_perturbation_model <- function(envir = parent.frame()) {
  model_path <- file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R")
  if (!file.exists(model_path)) stop("O2 supply-demand model was not found: ", model_path)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
  sys.source(model_path, envir = envir)
  invisible(model_path)
}

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

normalize_interaction_mode <- function(mode) {
  mode <- tolower(as.character(mode %||% "triggered_drug"))
  if (mode %in% c("drug", "triggered", "triggered-drug")) mode <- "triggered_drug"
  if (mode %in% c("untreated", "factorial", "untreated-factorial")) mode <- "untreated_factorial"
  if (mode %in% c("intermittent", "intermittent-drug", "cycle", "cyclic")) mode <- "intermittent_drug"
  allowed <- c("triggered_drug", "untreated_factorial", "intermittent_drug")
  if (!mode %in% allowed) {
    stop("Unknown interaction mode: ", mode, ". Use ", paste(allowed, collapse = ", "), ".")
  }
  mode
}

interaction_output_stem <- function(mode) {
  switch(
    normalize_interaction_mode(mode),
    triggered_drug = "invivo_triggered_drug_fixed_pre_1000post",
    untreated_factorial = "invivo_untreated_factorial_2N4N_1e6_100day",
    intermittent_drug = "invivo_intermittent_drug_7on21off_156week"
  )
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

format_value_label <- function(x, digits = 4) {
  format(signif(as.numeric(x), digits), trim = TRUE, scientific = TRUE)
}

read_required_tsv <- o2sd_read_required_tsv

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

write_artifact_manifest <- function(out_dir, artifact_rows, manifest_name) {
  manifest <- bind_rows(artifact_rows)
  required <- c("artifact", "relative_path", "role", "rows", "columns")
  missing_cols <- setdiff(required, names(manifest))
  if (length(missing_cols)) {
    stop("Artifact manifest rows are missing columns: ", paste(missing_cols, collapse = ", "))
  }
  manifest$relative_path <- as.character(manifest$relative_path)
  manifest$path <- normalizePath(file.path(out_dir, manifest$relative_path), mustWork = FALSE)
  manifest$exists <- file.exists(manifest$path)
  if (any(!manifest$exists)) {
    stop("Manifest references missing artifacts: ", paste(manifest$relative_path[!manifest$exists], collapse = ", "))
  }
  write_tsv(manifest, file.path(out_dir, manifest_name))
}

artifact_manifest_row <- function(artifact, relative_path, role, value) {
  data.frame(
    artifact = as.character(artifact),
    relative_path = as.character(relative_path),
    role = as.character(role),
    rows = if (is.data.frame(value)) nrow(value) else NA_integer_,
    columns = if (is.data.frame(value)) ncol(value) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

require_artifact_files <- function(base_dir, relative_paths, stage_label) {
  paths <- file.path(base_dir, relative_paths)
  missing <- relative_paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      stage_label, " requires pre-existing upstream artifacts in ", base_dir,
      ": ", paste(missing, collapse = ", "),
      ". Run the simulation producer first."
    )
  }
  invisible(paths)
}

validate_artifact_manifest <- function(base_dir,
                                       manifest_name,
                                       required_relative_paths,
                                       stage_label) {
  require_artifact_files(
    base_dir,
    c(manifest_name, required_relative_paths),
    stage_label
  )
  manifest <- read_required_tsv(file.path(base_dir, manifest_name))
  required_columns <- c("artifact", "relative_path", "role", "rows", "columns")
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns)) {
    stop(stage_label, " found an invalid ", manifest_name, "; missing columns: ",
         paste(missing_columns, collapse = ", "))
  }
  missing_entries <- setdiff(required_relative_paths, as.character(manifest$relative_path))
  if (length(missing_entries)) {
    stop(stage_label, " found an incomplete ", manifest_name, "; missing entries: ",
         paste(missing_entries, collapse = ", "))
  }
  invisible(manifest)
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
