#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))

.script_dir_local <- local({
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

SCRIPT_DIR <- normalizePath(.script_dir_local, mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(OXYGEN_ROOT, "code", "O2_supply_demand_MAP"), mustWork = FALSE)
MODEL_DIR <- normalizePath(file.path(WORKFLOW_ROOT, "model"), mustWork = FALSE)
UTIL_SCRIPT <- normalizePath(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), mustWork = FALSE)

source(UTIL_SCRIPT, local = environment())
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = MODEL_DIR)
source(file.path(MODEL_DIR, "model_O2_supply_demand_MAP.R"), local = environment())

parse_args <- o2sd_parse_args
as_num <- o2sd_as_num

normalize_cfg_for_viz <- function(cfg) {
  normalize_sim_cfg_common(cfg, context = "viz")
}

read_run_params <- function(fit_dir, cfg = NULL) {
  p <- file.path(fit_dir, "best_params.tsv")
  if (!file.exists(p)) stop("Missing best_params.tsv in: ", fit_dir)
  tab <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain columns: parameter, value")
  }
  vals <- setNames(as.numeric(tab$value), as.character(tab$parameter))
  needed <- c(
    "lam_max", "p_misseg", "k_o_mis",
    "gamma_loss", "p_wgd",
    "o2_S0", "kappa_O", "eta_o2",
    "alpha_o2", "gamma_growth",
    "mu_hp", "gamma_mu", "O2_crit", "n_O", "k_clear"
  )
  miss <- setdiff(needed, names(vals))
  if (length(miss) > 0) {
    stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
  }
  out <- as.list(vals[needed])
  p_mis_base_val <- if ("p_mis_base" %in% names(vals)) vals[["p_mis_base"]] else NULL
  o2_min_val <- if ("o2_min" %in% names(vals)) vals[["o2_min"]] else NULL
  if (!is.null(p_mis_base_val) && is.finite(p_mis_base_val)) out$p_mis_base <- as.numeric(p_mis_base_val)
  if (!is.null(o2_min_val) && is.finite(o2_min_val)) out$o2_min <- as.numeric(o2_min_val)
  if ("rho_2N" %in% names(vals) && is.finite(vals[["rho_2N"]]) && vals[["rho_2N"]] > 0) {
    out$rho_2N <- vals[["rho_2N"]]
  }
  if ("c_vol_2N_mm3" %in% names(vals) && is.finite(vals[["c_vol_2N_mm3"]]) && vals[["c_vol_2N_mm3"]] > 0) {
    out$c_vol_2N_mm3 <- vals[["c_vol_2N_mm3"]]
  }
  fit_treatment_use <- isTRUE(.first_non_null_local(cfg$fit_treatment, FALSE))
  if (fit_treatment_use) {
    miss_tx <- setdiff(c("alpha", "gamma"), names(vals))
    if (length(miss_tx) > 0) {
      stop("best_params.tsv missing treatment parameters while fit_treatment=TRUE: ", paste(miss_tx, collapse = ", "))
    }
    out$alpha <- vals[["alpha"]]
    out$gamma <- vals[["gamma"]]
  } else {
    out$alpha <- if ("alpha" %in% names(vals) && is.finite(vals[["alpha"]])) vals[["alpha"]] else 0
    out$gamma <- if ("gamma" %in% names(vals) && is.finite(vals[["gamma"]])) vals[["gamma"]] else 1
  }
  out$tau_O2 <- if ("tau_O2" %in% names(vals) && is.finite(vals[["tau_O2"]]) && vals[["tau_O2"]] > 0) vals[["tau_O2"]] else NULL
  normalize_run_params_common(out, cfg = cfg)
}

canonical_cohort_from_harvest <- function(harvest) {
  h <- as.character(harvest)
  if (grepl("2N", h, fixed = TRUE)) return("2N")
  if (grepl("4N", h, fixed = TRUE)) return("4N")
  NA_character_
}

read_input_table <- function(dt_path, sheet = 1) {
  dt <- readxl::read_excel(dt_path, sheet = sheet)
  required_cols <- c("Mouse ID", "Dose", "harvest", "Day of 1st treatment")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in treatment sheet: ", paste(missing_cols, collapse = ", "))
  }
  day_cols <- grep("^Day_[0-9]+$", names(dt), value = TRUE)
  if (length(day_cols) == 0L) {
    stop("No Day_* columns found in treatment sheet: ", dt_path)
  }
  dt
}

build_prediction_rows <- function(dt) {
  out <- dt %>%
    transmute(
      mouse_id = as.character(`Mouse ID`),
      dose = suppressWarnings(as.numeric(Dose)),
      harvest = as.character(harvest),
      treatment_day = suppressWarnings(as.numeric(`Day of 1st treatment`)),
      cohort = vapply(as.character(harvest), canonical_cohort_from_harvest, character(1))
    ) %>%
    filter(
      nzchar(harvest),
      is.finite(treatment_day),
      treatment_day >= 0,
      !is.na(cohort)
    )
  if (!nrow(out)) {
    stop("No valid harvest rows with finite treatment day and recognizable 2N/4N cohort.")
  }
  out
}

integerize_cell_count <- function(x, method = "round") {
  x_num <- suppressWarnings(as.numeric(x))
  x_num[!is.finite(x_num)] <- 0
  x_num[x_num < 0] <- 0
  method <- tolower(trimws(as.character(method[[1]])))
  vals <- switch(
    method,
    floor = floor(x_num),
    ceiling = ceiling(x_num),
    round = round(x_num),
    stop("Unsupported count_rounding method: ", method)
  )
  as.integer(pmax(vals, 0))
}

open_output_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "wt")
  } else {
    file(path, open = "wt")
  }
}

write_expanded_cell_rows <- function(count_tab, out_path, chunk_size = 100000L) {
  needed <- c("harvest", "treatment_day", "ploidy", "label", "cell_count_int")
  missing_cols <- setdiff(needed, names(count_tab))
  if (length(missing_cols) > 0L) {
    stop("Expanded writer missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  con <- open_output_connection(out_path)
  on.exit(close(con), add = TRUE)
  writeLines("harvest\ttreatment_day\tploidy\tlabel", con = con)
  total_rows <- 0L
  chunk_size <- as.integer(max(1L, chunk_size))

  for (i in seq_len(nrow(count_tab))) {
    n_rep <- as.integer(count_tab$cell_count_int[[i]])
    if (!is.finite(n_rep) || n_rep <= 0L) next
    line_i <- paste(
      as.character(count_tab$harvest[[i]]),
      format(as.numeric(count_tab$treatment_day[[i]]), trim = TRUE, scientific = FALSE),
      format(as.numeric(count_tab$ploidy[[i]]), trim = TRUE, scientific = FALSE),
      as.character(count_tab$label[[i]]),
      sep = "\t"
    )
    remaining <- n_rep
    while (remaining > 0L) {
      n_chunk <- min(chunk_size, remaining)
      writeLines(rep.int(line_i, n_chunk), con = con)
      remaining <- remaining - n_chunk
    }
    total_rows <- total_rows + n_rep
  }

  invisible(total_rows)
}

simulate_state_at_treatment_day <- function(run_params, cfg, harvest_row, min_count = 1e-6) {
  model_core <- build_model_core(cfg = cfg)
  grid_pre <- as.numeric(model_core$grid_pre)
  cohort_value <- canonical_cohort_from_harvest(.first_non_null_local(harvest_row$harvest[[1]], harvest_row$harvest))
  if (is.na(cohort_value)) {
    cohort_value <- as.character(.first_non_null_local(harvest_row$cohort[[1]], harvest_row$cohort))
  }
  init_state <- if (identical(cohort_value, "2N")) model_core$init_state_2N else model_core$init_state_4N

  treat_day <- as.numeric(harvest_row$treatment_day)
  treat_step <- as.integer(round(treat_day / cfg$DT))
  obs_steps_cpp <- as.integer(treat_step)

  o2_s0_upper_use <- as.numeric(.first_non_null_local(cfg$o2_S0_upper_bound, 5.0))
  if (!is.finite(o2_s0_upper_use) || o2_s0_upper_use <= 0) o2_s0_upper_use <- 5.0
  o2_S0_use <- as.numeric(.first_non_null_local(run_params$o2_S0, cfg$o2_S0_init, 0.5))
  if (!is.finite(o2_S0_use) || o2_S0_use < 0) o2_S0_use <- 0.5
  o2_S0_use <- min(max(o2_S0_use, 0), o2_s0_upper_use)
  kappa_O_use <- as.numeric(.first_non_null_local(run_params$kappa_O, cfg$kappa_O_init, 1.0))
  if (!is.finite(kappa_O_use) || kappa_O_use <= 0) kappa_O_use <- 1.0
  eta_o2_use <- as.numeric(.first_non_null_local(run_params$eta_o2, cfg$eta_o2_init, 1.0))
  if (!is.finite(eta_o2_use) || eta_o2_use <= 0) eta_o2_use <- 1.0
  O2_crit_use <- as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  o2_Nref_use <- as.numeric(.first_non_null_local(cfg$o2_Nref, cfg$init_total_size, 1e6))
  if (!is.finite(o2_Nref_use) || o2_Nref_use <= 0) o2_Nref_use <- 1e6
  o2_min_use <- as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.5))
  if (!is.finite(o2_min_use) || o2_min_use < 0) o2_min_use <- 0.5
  o2_min_use <- min(max(o2_min_use, 0), o2_s0_upper_use)
  tau_O2_use <- as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)
  o2_growth_use <- isTRUE(.first_non_null_local(cfg$O2_growth, TRUE))

  sim_cpp <- cpp_o2simps_simulate_one(
    init_state = as.numeric(init_state),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = obs_steps_cpp,
    sim_end_step = as.integer(treat_step),
    DT = as.numeric(cfg$DT),
    dose = as.numeric(harvest_row$dose),
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = as.numeric(treat_day),
    fit_treatment = isTRUE(cfg$fit_treatment),
    alpha = as.numeric(.first_non_null_local(run_params$alpha, 0.0)),
    gamma = as.numeric(.first_non_null_local(run_params$gamma, 1.0)),
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding_enabled = isTRUE(cfg$Crowding),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(O2_crit_use),
    o2_feedback = isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE)),
    o2_S0 = as.numeric(o2_S0_use),
    kappa_O = as.numeric(kappa_O_use),
    tau_O2 = as.numeric(tau_O2_use),
    o2_Nref = as.numeric(o2_Nref_use),
    o2_min = as.numeric(o2_min_use),
    eta_o2 = as.numeric(eta_o2_use),
    o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01)),
    o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005)),
    o2_cache_profile = isTRUE(.first_non_null_local(cfg$o2_cache_profile, FALSE)),
    O2_growth = isTRUE(o2_growth_use),
    lam_max = as.numeric(run_params$lam_max),
    has_p_misseg = !is.null(run_params$p_misseg),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = as.numeric(1e-8),
    gamma_loss = as.numeric(.first_non_null_local(run_params$gamma_loss, 0.1)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = 0.0,
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu, cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0)),
    ploidy_O2_death = assert_canonical_ploidy_o2_death_mode(
      .first_non_null_local(cfg$ploidy_O2_death, run_params$ploidy_O2_death, "diploid_NULL")
    ),
    start_with = assert_canonical_start_with_mode(
      .first_non_null_local(cfg$start_with, "ploidy")
    ),
    k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor),
    return_full_trajectory = TRUE
  )

  live_state <- as.numeric(as.matrix(sim_cpp$live_state_obs)[1, ])
  dead_h_state <- as.numeric(as.matrix(sim_cpp$dead_hypoxia_state_obs)[1, ])
  dead_b_state <- as.numeric(as.matrix(sim_cpp$dead_buffer_state_obs)[1, ])

  out <- bind_rows(
    data.frame(ploidy = grid_pre / as.numeric(cfg$N_UNIT), label = "live", cell_count = live_state, stringsAsFactors = FALSE),
    data.frame(ploidy = grid_pre / as.numeric(cfg$N_UNIT), label = "O2_dead", cell_count = dead_h_state, stringsAsFactors = FALSE),
    data.frame(ploidy = grid_pre / as.numeric(cfg$N_UNIT), label = "MS_dead", cell_count = dead_b_state, stringsAsFactors = FALSE)
  ) %>%
    mutate(
      harvest = as.character(harvest_row$harvest),
      treatment_day = as.numeric(treat_day)
    ) %>%
    select(harvest, treatment_day, ploidy, label, cell_count) %>%
    filter(is.finite(cell_count), cell_count > min_count) %>%
    arrange(match(label, c("live", "O2_dead", "MS_dead")), ploidy)

  out
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  fit_dir_arg <- .first_non_null_local(argv$fit_dir, NULL)
  if (is.null(fit_dir_arg) || !nzchar(trimws(as.character(fit_dir_arg)))) {
    stop(
      "Missing required argument --fit_dir. Example: ",
      "Rscript oxygen/scripts/predict_treatment_day_ploidy_from_fit.R ",
      "--fit_dir=oxygen/results/fit_invivo_o2_supply_demand_MAP/seed28"
    )
  }

  fit_dir <- normalizePath(
    fit_dir_arg,
    mustWork = TRUE
  )
  dt_path <- normalizePath(
    .first_non_null_local(
      argv$dt_path,
      file.path(dirname(OXYGEN_ROOT), "data", "InVivoData_Gemcitabine", "dt_Gem_VT_20260209_v5.xlsx")
    ),
    mustWork = TRUE
  )
  out_path <- normalizePath(
    .first_non_null_local(
      argv$out_path,
      file.path(OXYGEN_ROOT, "results", paste0("treatment_day_ploidy_counts_", basename(dirname(fit_dir)), "_", basename(fit_dir), ".tsv"))
    ),
    mustWork = FALSE
  )
  expanded_out_path <- normalizePath(
    .first_non_null_local(
      argv$expanded_out_path,
      file.path(OXYGEN_ROOT, "results", paste0("treatment_day_ploidy_cells_", basename(dirname(fit_dir)), "_", basename(fit_dir), ".tsv.gz"))
    ),
    mustWork = FALSE
  )
  min_count <- as_num(argv$min_count, 1e-6)
  if (!is.finite(min_count) || min_count < 0) min_count <- 1e-6
  sheet <- .first_non_null_local(argv$sheet, 1)
  max_harvests <- suppressWarnings(as.integer(.first_non_null_local(argv$max_harvests, NA_integer_)))
  count_rounding <- .first_non_null_local(argv$count_rounding, "round")
  write_expanded <- o2sd_as_bool(.first_non_null_local(argv$write_expanded, TRUE), TRUE)
  expand_chunk_size <- suppressWarnings(as.integer(.first_non_null_local(argv$expand_chunk_size, 100000L)))
  if (!is.finite(expand_chunk_size) || expand_chunk_size < 1L) expand_chunk_size <- 100000L

  cfg <- readRDS(file.path(fit_dir, "fit_config.rds"))
  cfg <- normalize_cfg_for_viz(cfg)
  run_params <- read_run_params(fit_dir, cfg = cfg)

  dt <- read_input_table(dt_path, sheet = sheet)
  pred_rows <- build_prediction_rows(dt)
  if (is.finite(max_harvests) && max_harvests > 0L) {
    pred_rows <- pred_rows %>% slice_head(n = max_harvests)
  }

  day_cols <- grep("^Day_[0-9]+$", names(dt), value = TRUE)
  valid_day_lookup <- unique(as.numeric(sub("^Day_", "", day_cols)))
  if (any(!pred_rows$treatment_day %in% valid_day_lookup)) {
    bad_days <- sort(unique(pred_rows$treatment_day[!pred_rows$treatment_day %in% valid_day_lookup]))
    warning("Some treatment days are not present as Day_* columns in the workbook: ", paste(bad_days, collapse = ", "))
  }

  result <- bind_rows(lapply(seq_len(nrow(pred_rows)), function(i) {
    simulate_state_at_treatment_day(run_params, cfg, pred_rows[i, , drop = FALSE], min_count = min_count)
  }))
  result$cell_count_int <- integerize_cell_count(result$cell_count, method = count_rounding)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write.table(result, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

  expanded_rows_written <- NA_integer_
  if (isTRUE(write_expanded)) {
    expanded_rows_written <- write_expanded_cell_rows(
      count_tab = result,
      out_path = expanded_out_path,
      chunk_size = expand_chunk_size
    )
  }

  message("Wrote treatment-day ploidy counts to: ", out_path)
  if (isTRUE(write_expanded)) {
    message("Wrote treatment-day expanded cell rows to: ", expanded_out_path)
    message("Expanded rows written: ", format(expanded_rows_written, scientific = TRUE))
  }
  message("Harvests: ", length(unique(result$harvest)))
  message("Rows: ", nrow(result))
}

if (sys.nframe() == 0) {
  main()
}
