#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEoptim))
suppressPackageStartupMessages(library(dplyr))

.o2g_joint_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2g_joint_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = FALSE)

source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2g_joint_bootstrap_script_dir)

parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
.first_non_null_local <- o2sd_first_non_null
clip <- o2sd_clip

load_backend_env_local <- function(path) {
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env, chdir = TRUE)
  env
}

INVIVO_ENV <- load_backend_env_local(
  file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_fit_invivo_backend.R")
)
INVITRO_ENV <- load_backend_env_local(
  file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_fit_invitro_backend.R")
)

trim_cli_scalar_local <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  txt <- trimws(as.character(x[[1]]))
  if (!nzchar(txt)) return(NULL)
  txt
}

path_or_default <- function(x, default) {
  out <- trim_cli_scalar_local(x)
  if (is.null(out)) default else out
}

resolve_joint_restriction_flags <- function(cfg_raw) {
  master_present <- !is.null(cfg_raw$joint_restriction)
  legacy_default <- as_bool(cfg_raw$joint_require_biological_constraints, FALSE)

  if (isTRUE(master_present)) {
    restriction <- as_bool(cfg_raw$joint_restriction, FALSE)
    return(list(
      joint_restriction = restriction,
      joint_require_biological_constraints = restriction,
      joint_require_invivo_pred1000_ploidy_gt2 = restriction,
      joint_require_invitro_growth_nonnegative = restriction,
      joint_require_invitro_ploidy_phenotype = restriction
    ))
  }

  require_invivo <- as_bool(
    cfg_raw$joint_require_invivo_pred1000_ploidy_gt2,
    legacy_default
  )
  require_growth <- as_bool(
    cfg_raw$joint_require_invitro_growth_nonnegative,
    legacy_default
  )
  require_ploidy <- as_bool(
    cfg_raw$joint_require_invitro_ploidy_phenotype,
    legacy_default
  )
  restriction <- isTRUE(require_invivo) || isTRUE(require_growth) || isTRUE(require_ploidy)

  list(
    joint_restriction = restriction,
    joint_require_biological_constraints = restriction,
    joint_require_invivo_pred1000_ploidy_gt2 = require_invivo,
    joint_require_invitro_growth_nonnegative = require_growth,
    joint_require_invitro_ploidy_phenotype = require_ploidy
  )
}

default_joint_out_dir <- function(cfg_raw) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_root <- path_or_default(cfg_raw$out_root, file.path(OXYGEN_ROOT, "results"))
  run_prefix <- path_or_default(cfg_raw$run_prefix, "fit_joint_O2G_supply_demand_MAP")
  append_ts <- as_bool(.first_non_null_local(cfg_raw$append_run_prefix_timestamp, TRUE), TRUE)
  ts_format <- path_or_default(cfg_raw$run_prefix_timestamp_format, "%Y%m%d_%H%M%S")
  if (isTRUE(append_ts)) {
    stamp <- format(Sys.time(), ts_format)
    run_prefix <- paste0(run_prefix, "_joint_", stamp)
  } else {
    run_prefix <- paste0(run_prefix, "_joint")
  }
  normalizePath(file.path(out_root, run_prefix), mustWork = FALSE)
}

safe_log10 <- function(x, fallback = 1e-12) {
  x <- as.numeric(x)
  if (!is.finite(x) || x <= 0) x <- fallback
  log10(x)
}

safe_qlogis <- function(x, name) {
  x <- as.numeric(x)
  if (!is.finite(x)) stop("Non-finite probability for ", name)
  qlogis(clip(x, 1e-12, 1 - 1e-12))
}

normalize_joint_n_cores <- function(x) {
  n <- suppressWarnings(as.integer(x))
  if (!is.finite(n) || is.na(n) || n < 1L) 1L else n
}

start_joint_deoptim_cluster <- function(n_cores) {
  n_use <- normalize_joint_n_cores(n_cores)
  if (n_use <= 1L) return(NULL)
  if (.Platform$OS.type == "unix" &&
      exists("makeForkCluster", envir = asNamespace("parallel"), inherits = FALSE)) {
    return(parallel::makeForkCluster(n_use))
  }
  parallel::makePSOCKcluster(n_use)
}

resolve_joint_raw_config <- function(argv) {
  parsed <- INVIVO_ENV$.runner_resolve_config(
    argv = argv,
    script_dir = SCRIPT_DIR,
    caller_wd = getwd()
  )
  parsed$cfg
}

build_joint_invivo_context <- function(cfg_raw) {
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg_raw$glucose, TRUE),
    default = TRUE
  ))

  parameter_table <- trim_cli_scalar_local(cfg_raw$parameter_table)
  if (is.null(parameter_table)) {
    parameter_table <- default_o2g_parameter_table_path_common(
      script_dir = SCRIPT_DIR,
      glucose = glucose_use,
      must_exist = TRUE
    )
  }
  if (!file.exists(parameter_table)) {
    stop("Joint fit in vivo parameter table not found: ", parameter_table)
  }

  model_path <- file.path(WORKFLOW_ROOT, "model", "model_O2G_supply_demand_MAP.R")
  if (!file.exists(model_path)) stop("Cannot find model_O2G_supply_demand_MAP.R at ", model_path)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
  sys.source(model_path, envir = INVIVO_ENV, chdir = TRUE)
  cpp_dll <- INVIVO_ENV$o2simps_cpp_dll_info()

  o2_S0_upper_arg <- read_o2_S0_natural_upper_bound_common(parameter_table, fallback = NA_real_)
  if (!is.finite(o2_S0_upper_arg) || o2_S0_upper_arg <= 0) {
    stop("Failed to read natural-scale o2_S0 upper bound from parameter_table.")
  }

  cfg <- list(
    model_path = model_path,
    cpp_dll_name = as.character(cpp_dll$name),
    cpp_dll_path = as.character(cpp_dll$path),
    cpp_wrapper_path = as.character(cpp_dll$wrapper_path),
    N_UNIT = as_int(cfg_raw$N_UNIT, 22L),
    chr_lengths_bp = INVIVO_ENV$default_chr_lengths_bp_1to22(),
    N_MIN = as_int(cfg_raw$N_MIN, 22L),
    N_MAX = as_int(cfg_raw$N_MAX, 154L),
    DT = as_num(.first_non_null_local(cfg_raw$dt, cfg_raw$DT), 0.5),
    o2_S0_upper_bound = o2_S0_upper_arg,
    ploidy_O2_death = canonical_ploidy_o2_death_mode(cfg_raw$ploidy_O2_death, "diploid_NULL"),
    glucose = glucose_use,
    start_with = canonical_start_with_mode(cfg_raw$start_with, "ploidy"),
    o2_burden_feedback = as_bool(cfg_raw$o2_burden_feedback, TRUE),
    O2_growth = as_bool(cfg_raw$O2_growth, TRUE),
    o2_cache_bin_pct = as_num(cfg_raw$o2_cache_bin_pct, 0.01),
    o2_cache_hysteresis_pct = as_num(cfg_raw$o2_cache_hysteresis_pct, 0.005),
    o2_cache_profile = as_bool(cfg_raw$o2_cache_profile, FALSE),
    o2_Nref = as_num(cfg_raw$o2_Nref, as_num(cfg_raw$init_total_size, 1e6)),
    o2_min = as_num(cfg_raw$o2_min, 0.0),
    Crowding = o2sd_as_bool_scalar(cfg_raw$Crowding, TRUE),
    fit_tau_O2 = as_bool(cfg_raw$fit_tau_O2, FALSE),
    parameter_table = normalizePath(parameter_table, mustWork = FALSE),
    K = as_num(cfg_raw$K, 1e12),
    crowding = if (!is.null(cfg_raw$crowding)) cfg_raw$crowding else "logistic",
    init_total_size = as_num(cfg_raw$init_total_size, 1e6),
    dose_ref = as_num(cfg_raw$dose_ref, 30),
    tx_mult_min = as_num(cfg_raw$tx_mult_min, 0.05),
    min_pop = as_num(cfg_raw$min_pop, 1e-12),
    sigma_ploidy = as_num(cfg_raw$sigma_ploidy, 0.08),
    burden_log_eps = as_num(cfg_raw$burden_log_eps, 1e-12),
    burden_exclude_day0 = as_bool(cfg_raw$burden_exclude_day0, TRUE),
    use_soft_prior = as_bool(cfg_raw$use_soft_prior, TRUE),
    lambda_prior = as_num(cfg_raw$lambda_prior, 0.1),
    harvest_init_multiplier = as_bool(cfg_raw$harvest_init_multiplier, FALSE),
    prior_center_log_init_mult = as_num(cfg_raw$prior_center_log_init_mult, 0.0),
    prior_sd_log_init_mult = as_num(cfg_raw$prior_sd_log_init_mult, 0.35),
    log_init_mult_lower = as_num(cfg_raw$log_init_mult_lower, -1.0),
    log_init_mult_upper = as_num(cfg_raw$log_init_mult_upper, 1.0),
    prior_center_log10_k_o = as_num(cfg_raw$prior_center_log10_k_o, log10(50)),
    prior_sd_log10_k_o = as_num(cfg_raw$prior_sd_log10_k_o, 1.0),
    prior_center_log10_kappa_O = as_num(cfg_raw$prior_center_log10_kappa_O, NA_real_),
    prior_sd_log10_kappa_O = as_num(cfg_raw$prior_sd_log10_kappa_O, 1.0),
    prior_center_log10_o2_S0 = as_num(cfg_raw$prior_center_log10_o2_S0, NA_real_),
    prior_sd_log10_o2_S0 = as_num(cfg_raw$prior_sd_log10_o2_S0, 0.5),
    prior_center_log10_eta_o2 = as_num(cfg_raw$prior_center_log10_eta_o2, NA_real_),
    prior_sd_log10_eta_o2 = as_num(cfg_raw$prior_sd_log10_eta_o2, 0.5),
    prior_center_buffer_smax = as_num(cfg_raw$prior_center_buffer_smax, 0.8),
    prior_sd_buffer_smax = as_num(cfg_raw$prior_sd_buffer_smax, 0.25),
    prior_center_log10_buffer_beta = as_num(cfg_raw$prior_center_log10_buffer_beta, 0.0),
    prior_sd_log10_buffer_beta = as_num(cfg_raw$prior_sd_log10_buffer_beta, 0.75),
    prior_center_log10_buffer_n_exp = as_num(cfg_raw$prior_center_log10_buffer_n_exp, 0.0),
    prior_sd_log10_buffer_n_exp = as_num(cfg_raw$prior_sd_log10_buffer_n_exp, 0.75),
    prior_center_log10_rho_2N = as_num(cfg_raw$prior_center_log10_rho_2N, NA_real_),
    prior_sd_log10_rho_2N = as_num(cfg_raw$prior_sd_log10_rho_2N, 0.35),
    prior_center_log10_mu_hp = as_num(cfg_raw$prior_center_log10_mu_hp, NA_real_),
    prior_sd_log10_mu_hp = as_num(cfg_raw$prior_sd_log10_mu_hp, 1.0),
    prior_center_gamma_mu = as_num(cfg_raw$prior_center_gamma_mu, NA_real_),
    prior_sd_gamma_mu = as_num(cfg_raw$prior_sd_gamma_mu, 0.5),
    prior_center_log10_k_clear = as_num(cfg_raw$prior_center_log10_k_clear, NA_real_),
    prior_sd_log10_k_clear = as_num(cfg_raw$prior_sd_log10_k_clear, 1.0),
    fit_treatment = as_bool(cfg_raw$fit_treatment, FALSE),
    dose_zero_only = as_bool(cfg_raw$dose_zero_only, TRUE),
    paired_only = as_bool(cfg_raw$paired_only, TRUE),
    truncate_at_treatment = as_bool(cfg_raw$truncate_at_treatment, FALSE),
    ploidy_at_harvest = as_bool(cfg_raw$ploidy_at_harvest, TRUE),
    use_deoptim = TRUE,
    deoptim_parallel = FALSE,
    de_init_mode = tolower(trimws(as.character(.first_non_null_local(cfg_raw$de_init_mode, "hybrid")))),
    de_init_uniform_frac = as_num(cfg_raw$de_init_uniform_frac, 0.3),
    de_init_sigma_frac = as_num(cfg_raw$de_init_sigma_frac, 0.1),
    de_reltol = as_num(cfg_raw$de_reltol, 1e-4),
    de_steptol = as_int(cfg_raw$de_steptol, 25L),
    itermax = as_int(cfg_raw$itermax, 40L),
    NP = as_int(cfg_raw$NP, 80L),
    n_cores = 1L,
    predict_n_cores = 1L,
    seed = as_int(cfg_raw$seed, 1L),
    max_scenarios = as_num(cfg_raw$max_scenarios, Inf),
    trace_obj = as_bool(cfg_raw$trace_obj, FALSE),
    optim_trace = as_bool(cfg_raw$optim_trace, TRUE),
    optim_trace_every = as_int(cfg_raw$optim_trace_every, 1L)
  )
  cfg <- normalize_sim_cfg_common(cfg, context = "fit")

  data_dir <- path_or_default(cfg_raw$data_dir, file.path(OXYGEN_ROOT, "data", "InVivoData_Gemcitabine"))
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- INVIVO_ENV$resolve_terminal_ploidy_path(data_dir)
  scenarios <- INVIVO_ENV$prepare_data(dt_path, ploidy_path, cfg)
  cfg$harvest_param_ids <- unique(vapply(scenarios, function(sc) as.character(sc$harvest), character(1)))

  param_bundle <- INVIVO_ENV$build_transformed_parameter_table(
    path = cfg$parameter_table,
    fit_treatment = cfg$fit_treatment,
    fit_tau_O2 = cfg$fit_tau_O2,
    O2_growth = cfg$O2_growth,
    glucose = cfg$glucose,
    harvest_init_multiplier = cfg$harvest_init_multiplier,
    harvest_ids = cfg$harvest_param_ids,
    prior_center_log_init_mult = cfg$prior_center_log_init_mult,
    log_init_mult_lower = cfg$log_init_mult_lower,
    log_init_mult_upper = cfg$log_init_mult_upper
  )
  cfg <- INVIVO_ENV$sync_cfg_from_natural_parameter_table(cfg, param_bundle$natural)
  cfg <- INVIVO_ENV$finalize_prior_defaults(cfg)
  cfg$model_core <- INVIVO_ENV$build_model_core(cfg = cfg)
  cfg$scenario_cpp <- INVIVO_ENV$prepare_cpp_scenarios(scenarios, cfg)

  list(cfg = cfg, scenarios = scenarios, param_bundle = param_bundle)
}

build_joint_invitro_context <- function(cfg_raw) {
  parameter_table <- trim_cli_scalar_local(.first_non_null_local(
    cfg_raw$invitro_parameter_table,
    cfg_raw$parameter_table_invitro
  ))
  if (is.null(parameter_table)) {
    parameter_table <- INVITRO_ENV$default_parameter_table_path(
      must_exist = TRUE
    )
  }
  fit_objects_dir <- path_or_default(
    cfg_raw$fit_objects_dir,
    INVITRO_ENV$default_fit_objects_dir(must_exist = TRUE)
  )
  flow_density_path <- INVITRO_ENV$resolve_optional_flow_density_path(cfg_raw$flow_density_path)
  cfg <- INVITRO_ENV$build_invitro_cfg(
    parameter_table = parameter_table,
    dt = as_num(.first_non_null_local(cfg_raw$invitro_dt, cfg_raw$dt, cfg_raw$DT), 0.1),
    init_total_size = as_num(.first_non_null_local(cfg_raw$invitro_init_total_size, cfg_raw$init_total_size), 1e6),
    o2_upper_bound = as_num(cfg_raw$invitro_o2_upper_bound, 21),
    fixed_oxygen = TRUE
  )
  fit_objects <- INVITRO_ENV$ivt_load_fit_objects_compat(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path
  )
  spec <- INVITRO_ENV$ivt_optimizer_spec(cfg)
  list(
    cfg = cfg,
    fit_objects = fit_objects,
    spec = spec,
    natural = read.csv(parameter_table, stringsAsFactors = FALSE),
    parameter_table = normalizePath(parameter_table, mustWork = FALSE),
    fit_objects_dir = normalizePath(fit_objects_dir, mustWork = FALSE),
    flow_density_path = normalizePath(flow_density_path, mustWork = FALSE)
  )
}

shared_invitro_param_names <- function(invivo_glucose) {
  loss_shared <- c("buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp")
  growth_shared <- c("log10_lam_max")
  out <- c(
    growth_shared, "log10_p_misseg",
    "log10_p_mis_base", "log10_k_o_mis", loss_shared, "log10_alpha_o2", "gamma_growth",
    "log10_mu_hp", "gamma_mu", "log10_O2_crit", "n_O", "log10_p_wgd"
  )
  out
}

joint_shared_natural_param_names <- function(invivo_glucose) {
  loss_shared <- c("buffer_smax", "buffer_beta", "buffer_n_exp")
  growth_shared <- c("lam_max")
  out <- c(
    growth_shared, "p_mis_base", "p_misseg", "k_o_mis",
    loss_shared, "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
    "O2_crit", "n_O", "p_wgd"
  )
  out
}

invitro_shared_param_name_for_natural <- function(symbol) {
  map <- c(
    lam_max = "log10_lam_max",
    p_mis_base = "log10_p_mis_base",
    p_misseg = "log10_p_misseg",
    k_o_mis = "log10_k_o_mis",
    buffer_smax = "buffer_smax",
    buffer_beta = "log10_buffer_beta",
    buffer_n_exp = "log10_buffer_n_exp",
    alpha_o2 = "log10_alpha_o2",
    gamma_growth = "gamma_growth",
    mu_hp = "log10_mu_hp",
    gamma_mu = "gamma_mu",
    O2_crit = "log10_O2_crit",
    n_O = "n_O",
    p_wgd = "log10_p_wgd"
  )
  unname(map[[as.character(symbol)]])
}

transform_invitro_shared_slot <- function(value, symbol, slot_label) {
  value <- as.numeric(value)
  if (!is.finite(value)) {
    stop("Non-finite ", slot_label, " for shared in vitro parameter ", symbol, call. = FALSE)
  }
  if (symbol %in% c("buffer_smax", "gamma_growth", "gamma_mu", "n_O")) {
    return(value)
  }
  if (value <= 0) {
    stop("Shared in vitro parameter ", symbol, " ", slot_label, " must be > 0 for log10 transform.", call. = FALSE)
  }
  log10(value)
}

merge_joint_shared_optimizer_bounds <- function(invivo,
                                                invitro,
                                                invivo_glucose) {
  shared_names <- joint_shared_natural_param_names(
    invivo_glucose = invivo_glucose
  )
  invivo_nat <- invivo$param_bundle$natural
  invitro_nat <- invitro$natural
  merged_nat <- invivo_nat
  summary_rows <- list()

  row_for <- function(tab, symbol, label) {
    out <- tab[as.character(tab$param_symbol) == as.character(symbol), , drop = FALSE]
    if (nrow(out) != 1L) {
      stop("Joint shared parameter '", symbol, "' missing or duplicated in ", label, " parameter table.", call. = FALSE)
    }
    out
  }

  for (symbol in shared_names) {
    inv_row <- row_for(invivo_nat, symbol, "in vivo")
    ivt_row <- row_for(invitro_nat, symbol, "in vitro")
    lower_joint <- min(as.numeric(inv_row$lower_bound[[1]]), as.numeric(ivt_row$lower_bound[[1]]), na.rm = TRUE)
    upper_joint <- max(as.numeric(inv_row$upper_bound[[1]]), as.numeric(ivt_row$upper_bound[[1]]), na.rm = TRUE)
    if (!is.finite(lower_joint) || !is.finite(upper_joint) || lower_joint > upper_joint) {
      stop("Invalid merged joint bounds for shared parameter '", symbol, "'.", call. = FALSE)
    }
    idx <- which(as.character(merged_nat$param_symbol) == as.character(symbol))
    merged_nat$lower_bound[idx] <- lower_joint
    merged_nat$upper_bound[idx] <- upper_joint
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      parameter = symbol,
      invivo_lower = as.numeric(inv_row$lower_bound[[1]]),
      invivo_upper = as.numeric(inv_row$upper_bound[[1]]),
      invitro_lower = as.numeric(ivt_row$lower_bound[[1]]),
      invitro_upper = as.numeric(ivt_row$upper_bound[[1]]),
      joint_lower = lower_joint,
      joint_upper = upper_joint,
      stringsAsFactors = FALSE
    )
  }

  invivo_init <- invivo$param_bundle$optimizer$init
  invivo_lower <- invivo$param_bundle$optimizer$lower
  invivo_upper <- invivo$param_bundle$optimizer$upper
  invitro_clip_lower <- setNames(as.numeric(invitro$spec$lower), as.character(invitro$spec$param_name))
  invitro_clip_upper <- setNames(as.numeric(invitro$spec$upper), as.character(invitro$spec$param_name))
  invivo_specs <- INVIVO_ENV$parameter_table_specs()

  merged_row <- function(symbol) {
    row_for(merged_nat, symbol, "merged")
  }

  for (symbol in shared_names) {
    spec_row <- invivo_specs[as.character(invivo_specs$param_symbol) == as.character(symbol), , drop = FALSE]
    if (nrow(spec_row) == 1L) {
      pname <- as.character(spec_row$param_name[[1]])
      if (pname %in% names(invivo_lower) && !identical(as.character(spec_row$transform[[1]]), "delta_lam")) {
        row <- merged_row(symbol)
        invivo_lower[[pname]] <- INVIVO_ENV$transform_param_slot(
          as.numeric(row$lower_bound[[1]]),
          as.character(spec_row$transform[[1]]),
          symbol,
          "lower_bound"
        )
        invivo_upper[[pname]] <- INVIVO_ENV$transform_param_slot(
          as.numeric(row$upper_bound[[1]]),
          as.character(spec_row$transform[[1]]),
          symbol,
          "upper_bound"
        )
      }
    }

    ivt_pname <- invitro_shared_param_name_for_natural(symbol)
    if (!is.null(ivt_pname) && ivt_pname %in% names(invitro_clip_lower)) {
      row <- merged_row(symbol)
      invitro_clip_lower[[ivt_pname]] <- transform_invitro_shared_slot(
        row$lower_bound[[1]],
        symbol,
        "lower_bound"
      )
      invitro_clip_upper[[ivt_pname]] <- transform_invitro_shared_slot(
        row$upper_bound[[1]],
        symbol,
        "upper_bound"
      )
    }
  }

  list(
    invivo_optimizer = list(init = invivo_init, lower = invivo_lower, upper = invivo_upper),
    invitro_clip = list(lower = invitro_clip_lower, upper = invitro_clip_upper),
    natural = merged_nat,
    summary = dplyr::bind_rows(summary_rows)
  )
}

split_joint_natural_parameter_tables <- function(invivo_param_df,
                                                 invitro_param_df,
                                                 invivo_glucose) {
  shared_names <- joint_shared_natural_param_names(
    invivo_glucose = invivo_glucose
  )
  invivo_shared <- invivo_param_df[invivo_param_df$parameter %in% shared_names, , drop = FALSE]
  invitro_shared <- invitro_param_df[invitro_param_df$parameter %in% shared_names, , drop = FALSE]
  shared_df <- dplyr::full_join(
    dplyr::rename(invivo_shared, in_vivo_value = value),
    dplyr::rename(invitro_shared, in_vitro_value = value),
    by = "parameter"
  )
  if (nrow(shared_df) > 0L) {
    shared_df$parameter <- as.character(shared_df$parameter)
    shared_df$shared_value <- ifelse(
      is.finite(shared_df$in_vivo_value),
      shared_df$in_vivo_value,
      shared_df$in_vitro_value
    )
    shared_df$value_difference <- shared_df$in_vivo_value - shared_df$in_vitro_value
    shared_df <- shared_df[order(match(shared_df$parameter, shared_names), shared_df$parameter), , drop = FALSE]
    shared_df <- shared_df[, c("parameter", "shared_value", "in_vivo_value", "in_vitro_value", "value_difference"), drop = FALSE]
  }
  invivo_only <- invivo_param_df[!(invivo_param_df$parameter %in% shared_names), , drop = FALSE]
  invitro_only <- invitro_param_df[!(invitro_param_df$parameter %in% shared_names), , drop = FALSE]
  list(
    shared = shared_df,
    invivo_only = invivo_only,
    invitro_only = invitro_only
  )
}

join_invitro_path_map <- function(df, ctx) {
  if (!is.data.frame(df) || !nrow(df) || !"segment_id" %in% names(df) || !"cohort" %in% names(df)) {
    return(df)
  }
  path_map <- dplyr::bind_rows(
    INVITRO_ENV$ivt_terminal_path_map(ctx$invitro$fit_objects$jobs_2N, cohort = "2N"),
    INVITRO_ENV$ivt_terminal_path_map(ctx$invitro$fit_objects$jobs_4N, cohort = "4N")
  )
  suppressWarnings(dplyr::left_join(df, path_map, by = c("cohort", "segment_id")))
}

build_invitro_transformed_from_joint <- function(invivo_run_params,
                                                 invitro_extra_t,
                                                 invitro_spec,
                                                 invivo_cfg,
                                                 clip_lower = NULL,
                                                 clip_upper = NULL) {
  par_t <- setNames(as.numeric(invitro_spec$init), invitro_spec$param_name)
  set_if_present <- function(name, value) {
    if (name %in% names(par_t)) par_t[[name]] <<- as.numeric(value)
  }
  set_if_present("log10_lam_max", safe_log10(invivo_run_params$lam_max))
  set_if_present("log10_p_mis_base", safe_log10(.first_non_null_local(invivo_run_params$p_mis_base, invivo_cfg$p_mis_base, 1e-5)))
  set_if_present("log10_p_misseg", safe_log10(invivo_run_params$p_misseg))
  set_if_present("log10_k_o_mis", safe_log10(invivo_run_params$k_o_mis))
  set_if_present("buffer_smax", invivo_run_params$buffer_smax)
  set_if_present("log10_buffer_beta", safe_log10(invivo_run_params$buffer_beta))
  set_if_present("log10_buffer_n_exp", safe_log10(invivo_run_params$buffer_n_exp))
  set_if_present("log10_p_wgd", safe_log10(invivo_run_params$p_wgd))
  set_if_present("log10_alpha_o2", safe_log10(invivo_run_params$alpha_o2))
  set_if_present("gamma_growth", invivo_run_params$gamma_growth)
  set_if_present("log10_mu_hp", safe_log10(invivo_run_params$mu_hp))
  set_if_present("gamma_mu", invivo_run_params$gamma_mu)
  set_if_present("log10_O2_crit", safe_log10(invivo_run_params$O2_crit))
  set_if_present("n_O", invivo_run_params$n_O)
  if (length(invitro_extra_t) > 0L) {
    par_t[names(invitro_extra_t)] <- as.numeric(invitro_extra_t)
  }
  lower_use <- setNames(as.numeric(invitro_spec$lower), as.character(invitro_spec$param_name))
  upper_use <- setNames(as.numeric(invitro_spec$upper), as.character(invitro_spec$param_name))
  if (!is.null(clip_lower)) lower_use[names(clip_lower)] <- as.numeric(clip_lower)
  if (!is.null(clip_upper)) upper_use[names(clip_upper)] <- as.numeric(clip_upper)
  pmin(pmax(par_t, lower_use[names(par_t)]), upper_use[names(par_t)])
}

build_joint_context <- function(argv) {
  cfg_raw <- resolve_joint_raw_config(argv)
  restriction_flags <- resolve_joint_restriction_flags(cfg_raw)
  invivo <- build_joint_invivo_context(cfg_raw)
  invitro <- build_joint_invitro_context(cfg_raw)

  invivo_names <- names(invivo$param_bundle$optimizer$init)
  shared_ivt <- shared_invitro_param_names(invivo_glucose = invivo$cfg$glucose)
  ivt_extra_names <- setdiff(invitro$spec$param_name, shared_ivt)
  ivt_extra_prefixed <- paste0("ivt__", ivt_extra_names)
  joint_bounds <- merge_joint_shared_optimizer_bounds(
    invivo = invivo,
    invitro = invitro,
    invivo_glucose = invivo$cfg$glucose
  )

  ivt_init <- setNames(as.numeric(invitro$spec$init[match(ivt_extra_names, invitro$spec$param_name)]), ivt_extra_prefixed)
  ivt_lower <- setNames(as.numeric(invitro$spec$lower[match(ivt_extra_names, invitro$spec$param_name)]), ivt_extra_prefixed)
  ivt_upper <- setNames(as.numeric(invitro$spec$upper[match(ivt_extra_names, invitro$spec$param_name)]), ivt_extra_prefixed)

  init <- c(joint_bounds$invivo_optimizer$init, ivt_init)
  lower <- c(joint_bounds$invivo_optimizer$lower, ivt_lower)
  upper <- c(joint_bounds$invivo_optimizer$upper, ivt_upper)

  list(
    raw = cfg_raw,
    invivo = invivo,
    invitro = invitro,
    joint_shared_bounds = joint_bounds$summary,
    invitro_clip_lower = joint_bounds$invitro_clip$lower,
    invitro_clip_upper = joint_bounds$invitro_clip$upper,
    invivo_names = invivo_names,
    ivt_extra_names = ivt_extra_names,
    ivt_extra_prefixed = ivt_extra_prefixed,
    init = init,
    lower = lower,
    upper = upper,
    joint_weight_invivo = as_num(cfg_raw$joint_weight_invivo, 1.0),
    joint_weight_invitro = as_num(cfg_raw$joint_weight_invitro, 1.0),
    seed = as_int(cfg_raw$seed, 1L),
    itermax = as_int(cfg_raw$itermax, 40L),
    NP = as_int(cfg_raw$NP, 80L),
    joint_np_min_factor = as_int(cfg_raw$joint_np_min_factor, 10L),
    joint_invitro_growth_weight = as_num(cfg_raw$joint_invitro_growth_weight, 1.0),
    joint_invitro_ploidy_weight = as_num(cfg_raw$joint_invitro_ploidy_weight, 1.0),
    joint_invitro_flow_weight = as_num(cfg_raw$joint_invitro_flow_weight, 1.0),
    joint_restriction = isTRUE(restriction_flags$joint_restriction),
    joint_require_biological_constraints = isTRUE(restriction_flags$joint_require_biological_constraints),
    joint_constraint_penalty = as_num(cfg_raw$joint_constraint_penalty, 1e9),
    joint_require_invivo_pred1000_ploidy_gt2 = isTRUE(restriction_flags$joint_require_invivo_pred1000_ploidy_gt2),
    joint_invivo_ploidy_horizon_day = as_num(cfg_raw$joint_invivo_ploidy_horizon_day, 1000),
    joint_invivo_min_ploidy_fold = as_num(cfg_raw$joint_invivo_min_ploidy_fold, 2),
    joint_require_invitro_growth_nonnegative = isTRUE(restriction_flags$joint_require_invitro_growth_nonnegative),
    joint_require_invitro_ploidy_phenotype = isTRUE(restriction_flags$joint_require_invitro_ploidy_phenotype),
    joint_invitro_2N_wgd_min_N = as_num(cfg_raw$joint_invitro_2N_wgd_min_N, 80),
    joint_invitro_2N_wgd_min_fraction = as_num(cfg_raw$joint_invitro_2N_wgd_min_fraction, 0.01),
    joint_invitro_4N_min_chr_drop = as_num(cfg_raw$joint_invitro_4N_min_chr_drop, 2),
    n_cores_requested = normalize_joint_n_cores(.first_non_null_local(cfg_raw$joint_n_cores, cfg_raw$n_cores, 1L)),
    n_cores_used = NA_integer_,
    out_dir = path_or_default(cfg_raw$out_dir, default_joint_out_dir(cfg_raw))
  )
}

read_joint_init_candidates <- function(path, ctx) {
  path <- trim_cli_scalar_local(path)
  if (is.null(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("joint_init_candidates_tsv not found: ", path, call. = FALSE)
  }
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  req <- c("candidate", "parameter", "transformed_value")
  if (!all(req %in% names(tab))) {
    stop(
      "joint_init_candidates_tsv must contain columns: ",
      paste(req, collapse = ", "),
      call. = FALSE
    )
  }

  tab$candidate <- as.character(tab$candidate)
  tab$parameter <- as.character(tab$parameter)
  tab$transformed_value <- suppressWarnings(as.numeric(tab$transformed_value))
  tab <- tab[nzchar(tab$candidate) & nzchar(tab$parameter), , drop = FALSE]
  retired_parameters <- c("ivt__log10_p_mis_base")
  tab <- tab[!(tab$parameter %in% retired_parameters), , drop = FALSE]
  if (nrow(tab) == 0L) {
    stop("joint_init_candidates_tsv contains no usable candidate rows: ", path, call. = FALSE)
  }
  if (any(!is.finite(tab$transformed_value))) {
    bad <- tab$parameter[!is.finite(tab$transformed_value)]
    stop(
      "joint_init_candidates_tsv contains non-finite transformed values for: ",
      paste(unique(bad), collapse = ", "),
      call. = FALSE
    )
  }

  full_names <- names(ctx$init)
  candidates <- unique(tab$candidate)
  mat <- matrix(NA_real_, nrow = length(candidates), ncol = length(full_names))
  rownames(mat) <- candidates
  colnames(mat) <- full_names
  used_rows <- vector("list", length(candidates))

  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    rows <- tab[tab$candidate == candidate, , drop = FALSE]
    dup <- rows$parameter[duplicated(rows$parameter)]
    if (length(dup) > 0L) {
      stop(
        "joint_init_candidates_tsv has duplicate rows for candidate '",
        candidate, "': ", paste(unique(dup), collapse = ", "),
        call. = FALSE
      )
    }
    if ("log10_p_mis_base" %in% full_names && !("log10_p_mis_base" %in% rows$parameter)) {
      compat_row <- rows[1L, , drop = FALSE]
      compat_row[,] <- NA
      compat_row$candidate <- candidate
      compat_row$parameter <- "log10_p_mis_base"
      compat_row$transformed_value <- as.numeric(ctx$init[["log10_p_mis_base"]])
      rows <- rbind(rows, compat_row)
    }
    missing_names <- setdiff(full_names, rows$parameter)
    extra_names <- setdiff(rows$parameter, full_names)
    if (length(missing_names) > 0L || length(extra_names) > 0L) {
      msg <- character(0)
      if (length(missing_names) > 0L) {
        msg <- c(msg, paste0("missing=", paste(missing_names, collapse = ",")))
      }
      if (length(extra_names) > 0L) {
        msg <- c(msg, paste0("extra=", paste(extra_names, collapse = ",")))
      }
      stop(
        "joint_init_candidates_tsv candidate '", candidate,
        "' does not match current joint parameter space (",
        paste(msg, collapse = "; "), ").",
        call. = FALSE
      )
    }
    vals <- setNames(rows$transformed_value, rows$parameter)[full_names]
    clipped <- clip(vals, ctx$lower, ctx$upper)
    mat[i, ] <- clipped
    used_rows[[i]] <- data.frame(
      candidate = candidate,
      parameter = full_names,
      transformed_value_input = as.numeric(vals),
      transformed_value_used = as.numeric(clipped),
      lower = as.numeric(ctx$lower),
      upper = as.numeric(ctx$upper),
      clipped = as.logical(clipped != vals),
      stringsAsFactors = FALSE
    )
  }

  list(
    path = normalizePath(path, mustWork = FALSE),
    matrix = mat,
    used = bind_rows(used_rows)
  )
}

build_joint_de_initialpop <- function(np, lower, upper, candidates = NULL) {
  n_par <- length(lower)
  pop <- matrix(runif(np * n_par), nrow = np, ncol = n_par)
  for (j in seq_len(n_par)) {
    pop[, j] <- lower[[j]] + pop[, j] * (upper[[j]] - lower[[j]])
  }
  colnames(pop) <- names(lower)
  if (!is.null(candidates) && nrow(candidates) > 0L) {
    n_candidates <- nrow(candidates)
    if (n_candidates > np) {
      stop("Number of joint init candidates exceeds DEoptim NP.", call. = FALSE)
    }
    pop[seq_len(n_candidates), ] <- candidates
  }
  pop
}

score_joint_init_candidates <- function(candidate_mat, ctx) {
  if (is.null(candidate_mat) || nrow(candidate_mat) == 0L) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(candidate_mat)), function(i) {
    candidate <- rownames(candidate_mat)[[i]]
    par_t <- as.numeric(candidate_mat[i, ])
    names(par_t) <- colnames(candidate_mat)
    comp <- tryCatch(
      joint_objective_components(par_t, ctx),
      error = function(e) e
    )
    if (inherits(comp, "error")) {
      return(data.frame(
        candidate = candidate,
        objective_joint = 1e9,
        objective_invivo = NA_real_,
        objective_invitro = NA_real_,
        objective_invivo_data = NA_real_,
        objective_invivo_prior = NA_real_,
        objective_invivo_burden = NA_real_,
        objective_invivo_ploidy = NA_real_,
        error = conditionMessage(comp),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      candidate = candidate,
      objective_joint = as.numeric(comp$objective),
      objective_invivo = as.numeric(comp$invivo$L),
      objective_invitro = as.numeric(comp$invitro$objective),
      objective_invivo_data = as.numeric(comp$invivo$L_data),
      objective_invivo_prior = as.numeric(comp$invivo$L_prior),
      objective_invivo_burden = as.numeric(comp$invivo$L_b),
      objective_invivo_ploidy = as.numeric(comp$invivo$L_p),
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

deduplicate_joint_rows <- function(df, cols) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  key_cols <- intersect(cols, names(df))
  if (!length(key_cols)) return(df)
  key <- do.call(paste, c(df[key_cols], sep = "\r"))
  df[!duplicated(key), , drop = FALSE]
}

joint_invivo_pred1000_ploidy_gate <- function(invivo_run_params, ctx) {
  horizon <- as_num(ctx$joint_invivo_ploidy_horizon_day, 1000)
  threshold <- as_num(ctx$joint_invivo_min_ploidy_fold, 2)
  out <- list(
    pass = FALSE,
    invivo_pred1000_min_ploidy_fold = NA_real_,
    invivo_pred1000_threshold_ploidy_fold = threshold,
    invivo_pred1000_n_rows = 0L,
    invivo_pred1000_status = "not_evaluated"
  )
  if (!is.finite(horizon) || horizon <= 0 || !is.finite(threshold)) {
    out$invivo_pred1000_status <- "invalid_threshold"
    return(out)
  }

  vals <- tryCatch({
    cfg <- ctx$invivo$cfg
    model_core <- if (!is.null(cfg$model_core)) cfg$model_core else INVIVO_ENV$build_model_core(cfg = cfg)
    dplyr::bind_rows(lapply(ctx$invivo$scenarios, function(sc) {
      sc$obs_days <- horizon
      sc$sim_end_day <- horizon
      sim <- INVIVO_ENV$simulate_one(
        run_params = invivo_run_params,
        scenario = sc,
        cfg = cfg,
        model_core = model_core
      )
      frac <- as.numeric(sim$frac_N)
      n_grid <- suppressWarnings(as.numeric(names(sim$frac_N)))
      keep <- is.finite(frac) & is.finite(n_grid) & frac >= 0
      if (!any(keep)) {
        return(data.frame(
          cohort = as.character(sc$cohort),
          weighted_mean_N = NA_real_,
          weighted_mean_ploidy_fold = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      frac <- frac[keep]
      n_grid <- n_grid[keep]
      weighted_mean_N <- sum(n_grid * frac, na.rm = TRUE) / pmax(sum(frac, na.rm = TRUE), 1e-12)
      data.frame(
        cohort = as.character(sc$cohort),
        weighted_mean_N = weighted_mean_N,
        weighted_mean_ploidy_fold = weighted_mean_N / pmax(as.numeric(cfg$N_UNIT), 1e-12),
        stringsAsFactors = FALSE
      )
    }))
  }, error = function(e) {
    structure(list(message = conditionMessage(e)), class = "joint_gate_error")
  })

  if (inherits(vals, "joint_gate_error")) {
    out$invivo_pred1000_status <- paste0("error: ", vals$message)
    return(out)
  }
  vals <- vals[is.finite(vals$weighted_mean_ploidy_fold), , drop = FALSE]
  out$invivo_pred1000_n_rows <- nrow(vals)
  if (!nrow(vals)) {
    out$invivo_pred1000_status <- "missing_prediction"
    return(out)
  }
  out$invivo_pred1000_min_ploidy_fold <- min(vals$weighted_mean_ploidy_fold, na.rm = TRUE)
  out$pass <- isTRUE(all(vals$weighted_mean_ploidy_fold > threshold))
  out$invivo_pred1000_status <- if (isTRUE(out$pass)) "pass" else "fail"
  out
}

joint_invitro_ploidy_phenotype_gate <- function(invitro_comp, ctx) {
  wgd_min_N <- as_num(ctx$joint_invitro_2N_wgd_min_N, 80)
  wgd_min_fraction <- as_num(ctx$joint_invitro_2N_wgd_min_fraction, 0.01)
  drop_min <- as_num(ctx$joint_invitro_4N_min_chr_drop, 2)
  out <- list(
    pass = FALSE,
    invitro_2N_deprived_wgd_pass = FALSE,
    invitro_2N_deprived_wgd_max_fraction = NA_real_,
    invitro_2N_deprived_wgd_min_N = wgd_min_N,
    invitro_2N_deprived_wgd_min_fraction = wgd_min_fraction,
    invitro_2N_deprived_max_mean_N = NA_real_,
    invitro_4N_deprived_chr_drop_pass = FALSE,
    invitro_4N_deprived_chr_drop = NA_real_,
    invitro_4N_deprived_min_chr_drop_required = drop_min,
    invitro_4N_deprived_initial_mean_N = NA_real_,
    invitro_4N_deprived_min_mean_N = NA_real_
  )

  summary_df <- join_invitro_path_map(invitro_comp$summary, ctx)
  if (!is.data.frame(summary_df) || !nrow(summary_df)) return(out)
  if (!"lineage_label" %in% names(summary_df)) summary_df$lineage_label <- NA_character_
  if (!"lineage_terminal_key" %in% names(summary_df)) summary_df$lineage_terminal_key <- NA_character_
  if (!"lineage_passage_index" %in% names(summary_df)) summary_df$lineage_passage_index <- summary_df$passage_index
  summary_df <- deduplicate_joint_rows(
    summary_df,
    c("cohort", "segment_id", "passage_index", "oxygen_pct", "lineage_label", "lineage_terminal_key")
  )
  summary_df$predicted_mean_kary_N <- suppressWarnings(as.numeric(summary_df$predicted_mean_kary_N))
  summary_df$lineage_passage_index <- suppressWarnings(as.numeric(summary_df$lineage_passage_index))
  summary_df$passage_index <- suppressWarnings(as.numeric(summary_df$passage_index))

  d2_mean <- summary_df[
    as.character(summary_df$cohort) == "2N" &
      as.character(summary_df$lineage_label) == "deprived" &
      is.finite(summary_df$predicted_mean_kary_N),
    , drop = FALSE
  ]
  if (nrow(d2_mean)) {
    out$invitro_2N_deprived_max_mean_N <- max(d2_mean$predicted_mean_kary_N, na.rm = TRUE)
  }

  dist <- tryCatch(
    INVITRO_ENV$ivt_collect_distribution_summary(invitro_comp$run_2N),
    error = function(e) NULL
  )
  if (is.data.frame(dist) && nrow(dist)) {
    dist <- join_invitro_path_map(dist, ctx)
    if (!"lineage_label" %in% names(dist)) dist$lineage_label <- NA_character_
    dist <- deduplicate_joint_rows(
      dist,
      c("cohort", "segment_id", "passage_index", "oxygen_pct", "lineage_label", "N")
    )
    dist$N <- suppressWarnings(as.numeric(dist$N))
    dist$fraction <- suppressWarnings(as.numeric(dist$fraction))
    d2 <- dist[
      as.character(dist$cohort) == "2N" &
        as.character(dist$lineage_label) == "deprived" &
        is.finite(dist$N) &
        is.finite(dist$fraction) &
        dist$N >= wgd_min_N,
      , drop = FALSE
    ]
    if (nrow(d2)) {
      group_cols <- intersect(
        c("cohort", "segment_id", "passage_index", "oxygen_pct", "lineage_label"),
        names(d2)
      )
      if (length(group_cols)) {
        formula_use <- stats::as.formula(paste("fraction ~", paste(group_cols, collapse = " + ")))
        d2_frac <- stats::aggregate(formula_use, data = d2, FUN = sum, na.rm = TRUE)
        out$invitro_2N_deprived_wgd_max_fraction <- max(d2_frac$fraction, na.rm = TRUE)
      } else {
        out$invitro_2N_deprived_wgd_max_fraction <- sum(d2$fraction, na.rm = TRUE)
      }
    }
  }
  if (is.finite(out$invitro_2N_deprived_wgd_max_fraction)) {
    out$invitro_2N_deprived_wgd_pass <- isTRUE(
      out$invitro_2N_deprived_wgd_max_fraction >= wgd_min_fraction
    )
  } else {
    out$invitro_2N_deprived_wgd_pass <- isTRUE(
      is.finite(out$invitro_2N_deprived_max_mean_N) &&
        out$invitro_2N_deprived_max_mean_N >= wgd_min_N
    )
  }

  d4 <- summary_df[
    as.character(summary_df$cohort) == "4N" &
      as.character(summary_df$lineage_label) == "deprived" &
      is.finite(summary_df$predicted_mean_kary_N),
    , drop = FALSE
  ]
  if (nrow(d4)) {
    lineage_key <- if ("lineage_terminal_key" %in% names(d4)) {
      as.character(d4$lineage_terminal_key)
    } else {
      rep("deprived", nrow(d4))
    }
    drops <- lapply(split(d4, lineage_key), function(x) {
      x <- x[order(x$lineage_passage_index, x$passage_index, na.last = TRUE), , drop = FALSE]
      vals <- x$predicted_mean_kary_N[is.finite(x$predicted_mean_kary_N)]
      if (!length(vals)) return(NULL)
      data.frame(
        initial_mean_N = vals[[1]],
        min_mean_N = min(vals, na.rm = TRUE),
        chr_drop = min(vals, na.rm = TRUE) - vals[[1]],
        stringsAsFactors = FALSE
      )
    })
    drop_df <- dplyr::bind_rows(Filter(Negate(is.null), drops))
    if (nrow(drop_df)) {
      idx <- which.min(drop_df$chr_drop)
      out$invitro_4N_deprived_initial_mean_N <- drop_df$initial_mean_N[[idx]]
      out$invitro_4N_deprived_min_mean_N <- drop_df$min_mean_N[[idx]]
      out$invitro_4N_deprived_chr_drop <- drop_df$chr_drop[[idx]]
      out$invitro_4N_deprived_chr_drop_pass <- isTRUE(
        out$invitro_4N_deprived_chr_drop <= -drop_min
      )
    }
  }

  out$pass <- isTRUE(out$invitro_2N_deprived_wgd_pass) &&
    isTRUE(out$invitro_4N_deprived_chr_drop_pass)
  out
}

joint_constraint_metrics <- function(invivo_run_params, invitro_comp, ctx) {
  n_growth_negative <- as_num(invitro_comp$n_growth_negative_pred, NA_real_)
  growth_pass <- is.finite(n_growth_negative) && n_growth_negative == 0
  metrics <- list(
    joint_restriction = isTRUE(ctx$joint_restriction),
    joint_constraint_penalty = as_num(ctx$joint_constraint_penalty, 1e9),
    joint_require_invivo_pred1000_ploidy_gt2 = isTRUE(ctx$joint_require_invivo_pred1000_ploidy_gt2),
    joint_require_invitro_growth_nonnegative = isTRUE(ctx$joint_require_invitro_growth_nonnegative),
    joint_require_invitro_ploidy_phenotype = isTRUE(ctx$joint_require_invitro_ploidy_phenotype),
    invitro_growth_nonnegative_pass = growth_pass,
    invitro_n_growth_negative_pred = n_growth_negative
  )

  invivo_gate <- if (isTRUE(ctx$joint_require_invivo_pred1000_ploidy_gt2)) {
    joint_invivo_pred1000_ploidy_gate(invivo_run_params, ctx)
  } else {
    list(
      pass = NA,
      invivo_pred1000_min_ploidy_fold = NA_real_,
      invivo_pred1000_threshold_ploidy_fold = as_num(ctx$joint_invivo_min_ploidy_fold, 2),
      invivo_pred1000_n_rows = NA_integer_,
      invivo_pred1000_status = "disabled"
    )
  }
  names(invivo_gate)[names(invivo_gate) == "pass"] <- "invivo_pred1000_ploidy_pass"
  metrics <- c(metrics, invivo_gate)

  ploidy_gate <- if (isTRUE(ctx$joint_require_invitro_ploidy_phenotype)) {
    joint_invitro_ploidy_phenotype_gate(invitro_comp, ctx)
  } else {
    list(
      pass = NA,
      invitro_2N_deprived_wgd_pass = NA,
      invitro_2N_deprived_wgd_max_fraction = NA_real_,
      invitro_2N_deprived_wgd_min_N = as_num(ctx$joint_invitro_2N_wgd_min_N, 80),
      invitro_2N_deprived_wgd_min_fraction = as_num(ctx$joint_invitro_2N_wgd_min_fraction, 0.01),
      invitro_2N_deprived_max_mean_N = NA_real_,
      invitro_4N_deprived_chr_drop_pass = NA,
      invitro_4N_deprived_chr_drop = NA_real_,
      invitro_4N_deprived_min_chr_drop_required = as_num(ctx$joint_invitro_4N_min_chr_drop, 2),
      invitro_4N_deprived_initial_mean_N = NA_real_,
      invitro_4N_deprived_min_mean_N = NA_real_
    )
  }
  names(ploidy_gate)[names(ploidy_gate) == "pass"] <- "invitro_ploidy_phenotype_pass"
  metrics <- c(metrics, ploidy_gate)

  failures <- 0L
  if (isTRUE(metrics$joint_require_invivo_pred1000_ploidy_gt2) &&
      !isTRUE(metrics$invivo_pred1000_ploidy_pass)) failures <- failures + 1L
  if (isTRUE(metrics$joint_require_invitro_growth_nonnegative) &&
      !isTRUE(metrics$invitro_growth_nonnegative_pass)) failures <- failures + 1L
  if (isTRUE(metrics$joint_require_invitro_ploidy_phenotype) &&
      !isTRUE(metrics$invitro_ploidy_phenotype_pass)) failures <- failures + 1L

  penalty <- as_num(metrics$joint_constraint_penalty, 1e9)
  if (!is.finite(penalty) || penalty < 0) penalty <- 1e9
  metrics$joint_constraint_failures <- failures
  metrics$joint_constraint_penalty_total <- penalty * failures
  metrics$joint_constraints_pass <- failures == 0L
  metrics
}

joint_constraint_component_df <- function(metrics) {
  if (is.null(metrics) || !length(metrics)) {
    return(data.frame(component = character(0), value = character(0), stringsAsFactors = FALSE))
  }
  data.frame(
    component = names(metrics),
    value = vapply(metrics, function(x) {
      if (is.null(x) || !length(x)) return(NA_character_)
      as.character(x[[1]])
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

joint_objective_components <- function(par_t, ctx) {
  par_t <- as.numeric(par_t)
  names(par_t) <- names(ctx$init)
  invivo_par <- par_t[ctx$invivo_names]
  invivo_comp <- INVIVO_ENV$evaluate_objective_components(
    invivo_par,
    scenarios = ctx$invivo$scenarios,
    cfg = ctx$invivo$cfg
  )
  invivo_run_params <- INVIVO_ENV$decode_params(
    invivo_par,
    fit_treatment = ctx$invivo$cfg$fit_treatment,
    fit_tau_O2 = ctx$invivo$cfg$fit_tau_O2,
    cfg = ctx$invivo$cfg
  )
  ivt_extra <- numeric(0)
  if (length(ctx$ivt_extra_prefixed) > 0L) {
    ivt_extra <- par_t[ctx$ivt_extra_prefixed]
    names(ivt_extra) <- ctx$ivt_extra_names
  }
  ivt_par <- build_invitro_transformed_from_joint(
    invivo_run_params = invivo_run_params,
    invitro_extra_t = ivt_extra,
    invitro_spec = ctx$invitro$spec,
    invivo_cfg = ctx$invivo$cfg,
    clip_lower = ctx$invitro_clip_lower,
    clip_upper = ctx$invitro_clip_upper
  )
  invitro_run_params <- INVITRO_ENV$ivt_optim_par_to_run_params(ivt_par, cfg = ctx$invitro$cfg)
  invitro_run_params$p_mis_base <- as.numeric(.first_non_null_local(
    invivo_run_params$p_mis_base,
    ctx$invivo$cfg$p_mis_base,
    1e-5
  ))
  invitro_run_params$glucose <- FALSE
  invitro_comp <- tryCatch(
    INVITRO_ENV$ivt_objective_components(
      run_params = invitro_run_params,
      fit_objects = ctx$invitro$fit_objects,
      cfg = ctx$invitro$cfg,
      fallback_max_passage_days = 14,
      growth_weight = ctx$joint_invitro_growth_weight,
      ploidy_weight = ctx$joint_invitro_ploidy_weight,
      flow_weight = ctx$joint_invitro_flow_weight
    ),
    error = function(e) {
      INVITRO_ENV$make_penalty_components(reason = paste0("invitro_error: ", conditionMessage(e)))
    }
  )
  joint <- ctx$joint_weight_invivo * invivo_comp$L +
    ctx$joint_weight_invitro * invitro_comp$objective
  constraint_metrics <- joint_constraint_metrics(
    invivo_run_params = invivo_run_params,
    invitro_comp = invitro_comp,
    ctx = ctx
  )
  joint <- joint + as_num(constraint_metrics$joint_constraint_penalty_total, 0)
  if (!is.finite(joint)) joint <- 1e9
  list(
    objective = joint,
    objective_unpenalized = ctx$joint_weight_invivo * invivo_comp$L +
      ctx$joint_weight_invitro * invitro_comp$objective,
    constraint_metrics = constraint_metrics,
    invivo = invivo_comp,
    invitro = invitro_comp,
    invivo_run_params = invivo_run_params,
    invitro_run_params = invitro_run_params,
    invitro_transformed = ivt_par
  )
}

write_tsv_if_nonempty <- function(df, path) {
  if (is.data.frame(df) && nrow(df) > 0L) {
    write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

write_joint_outputs <- function(best_par_t, best_comp, ctx, out_dir, de_fit, local_fit, optimizer_trace = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    data.frame(
      transformed_parameter = names(best_par_t),
      transformed_value = as.numeric(best_par_t),
      stringsAsFactors = FALSE
    ),
    file = file.path(out_dir, "best_params_transformed.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invivo_params <- best_comp$invivo_run_params[vapply(best_comp$invivo_run_params, is.numeric, logical(1))]
  invivo_params <- filter_family_specific_run_params_for_output_common(
    invivo_params,
    glucose = ctx$invivo$cfg$glucose
  )
  invitro_params <- best_comp$invitro_run_params[vapply(best_comp$invitro_run_params, is.numeric, logical(1))]
  invitro_params <- filter_family_specific_run_params_for_output_common(
    invitro_params,
    glucose = FALSE
  )
  invivo_param_df <- data.frame(
    parameter = names(invivo_params),
    value = as.numeric(unlist(invivo_params)),
    stringsAsFactors = FALSE
  )
  invitro_param_df <- data.frame(
    parameter = names(invitro_params),
    value = as.numeric(unlist(invitro_params)),
    stringsAsFactors = FALSE
  )
  param_long_df <- bind_rows(
    data.frame(scope = "shared_invivo", invivo_param_df, stringsAsFactors = FALSE),
    data.frame(scope = "invitro_effective", invitro_param_df, stringsAsFactors = FALSE)
  )
  param_tables <- split_joint_natural_parameter_tables(
    invivo_param_df = invivo_param_df,
    invitro_param_df = invitro_param_df,
    invivo_glucose = ctx$invivo$cfg$glucose
  )
  write.table(invivo_param_df, file = file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(invitro_param_df, file = file.path(out_dir, "invitro_effective_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_long_df, file = file.path(out_dir, "joint_best_params_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_tables$shared, file = file.path(out_dir, "joint_params_shared.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_tables$invivo_only, file = file.path(out_dir, "joint_params_invivo_only.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_tables$invitro_only, file = file.path(out_dir, "joint_params_invitro_only.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ctx$joint_shared_bounds, file = file.path(out_dir, "joint_shared_bounds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  preds <- INVIVO_ENV$collect_predictions(best_comp$invivo_run_params, ctx$invivo$scenarios, ctx$invivo$cfg)
  write.table(preds$burden, file = file.path(out_dir, "invivo_burden_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$ploidy, file = file.path(out_dir, "invivo_terminal_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$burden, file = file.path(out_dir, "burden_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$ploidy, file = file.path(out_dir, "terminal_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  write_tsv_if_nonempty(join_invitro_path_map(best_comp$invitro$summary, ctx), file.path(out_dir, "invitro_lineage_summary.tsv"))
  write_tsv_if_nonempty(join_invitro_path_map(best_comp$invitro$growth_df, ctx), file.path(out_dir, "invitro_growth_loglik.tsv"))
  write_tsv_if_nonempty(join_invitro_path_map(best_comp$invitro$ploidy_df, ctx), file.path(out_dir, "invitro_ploidy_loglik.tsv"))
  write_tsv_if_nonempty(join_invitro_path_map(best_comp$invitro$flow_df, ctx), file.path(out_dir, "invitro_flow_loglik.tsv"))
  write_tsv_if_nonempty(join_invitro_path_map(best_comp$invitro$flow_overlay_df, ctx), file.path(out_dir, "invitro_flow_overlay.tsv"))

  dist_summary <- dplyr::bind_rows(
    INVITRO_ENV$ivt_collect_distribution_summary(best_comp$invitro$run_2N),
    INVITRO_ENV$ivt_collect_distribution_summary(best_comp$invitro$run_4N)
  )
  write_tsv_if_nonempty(join_invitro_path_map(dist_summary, ctx), file.path(out_dir, "invitro_distribution_summary.tsv"))

  ploidy_quantile_probs <- seq(0.01, 0.99, length.out = 50L)
  dist_quantiles <- dplyr::bind_rows(
    INVITRO_ENV$ivt_collect_distribution_quantiles(best_comp$invitro$run_2N, probs = ploidy_quantile_probs),
    INVITRO_ENV$ivt_collect_distribution_quantiles(best_comp$invitro$run_4N, probs = ploidy_quantile_probs)
  )
  write_tsv_if_nonempty(join_invitro_path_map(dist_quantiles, ctx), file.path(out_dir, "invitro_distribution_quantiles.tsv"))

  daily_counts <- dplyr::bind_rows(
    INVITRO_ENV$ivt_collect_daily_counts(best_comp$invitro$run_2N),
    INVITRO_ENV$ivt_collect_daily_counts(best_comp$invitro$run_4N)
  )
  write_tsv_if_nonempty(join_invitro_path_map(daily_counts, ctx), file.path(out_dir, "invitro_daily_counts.tsv"))

  observed_kary <- dplyr::bind_rows(
    INVITRO_ENV$ivt_collect_observed_kary_summary(best_comp$invitro$run_2N, ctx$invitro$fit_objects$fit_data),
    INVITRO_ENV$ivt_collect_observed_kary_summary(best_comp$invitro$run_4N, ctx$invitro$fit_objects$fit_data)
  )
  write_tsv_if_nonempty(join_invitro_path_map(observed_kary, ctx), file.path(out_dir, "invitro_observed_kary.tsv"))

  observed_flow <- dplyr::bind_rows(
    INVITRO_ENV$ivt_collect_observed_flow_summary(best_comp$invitro$run_2N, ctx$invitro$fit_objects$fit_data),
    INVITRO_ENV$ivt_collect_observed_flow_summary(best_comp$invitro$run_4N, ctx$invitro$fit_objects$fit_data)
  )
  write_tsv_if_nonempty(join_invitro_path_map(observed_flow, ctx), file.path(out_dir, "invitro_observed_flow.tsv"))

  joint_components <- data.frame(
    component = c(
      "objective_joint",
      "objective_joint_unpenalized",
      "weight_invivo",
      "weight_invitro",
      "objective_invivo",
      "objective_invivo_data",
      "objective_invivo_prior",
      "objective_invivo_burden",
      "objective_invivo_ploidy",
      "objective_invitro",
      "invitro_growth_loglik",
      "invitro_ploidy_loglik",
      "invitro_flow_loglik"
    ),
    value = as.character(c(
      best_comp$objective,
      best_comp$objective_unpenalized,
      ctx$joint_weight_invivo,
      ctx$joint_weight_invitro,
      best_comp$invivo$L,
      best_comp$invivo$L_data,
      best_comp$invivo$L_prior,
      best_comp$invivo$L_b,
      best_comp$invivo$L_p,
      best_comp$invitro$objective,
      best_comp$invitro$growth_loglik,
      best_comp$invitro$ploidy_loglik,
      best_comp$invitro$flow_loglik
    )),
    stringsAsFactors = FALSE
  )
  joint_components <- dplyr::bind_rows(
    joint_components,
    joint_constraint_component_df(best_comp$constraint_metrics)
  )
  write.table(joint_components, file = file.path(out_dir, "joint_components.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (is.null(optimizer_trace)) optimizer_trace <- list()
  optimizer_method <- as.character(.first_non_null_local(
    optimizer_trace$method,
    if (isTRUE(ctx$n_cores_used > 1L)) "DEoptim_parallel_plus_LBFGSB_serial" else "DEoptim_serial_plus_LBFGSB_serial"
  ))
  optimizer_deoptim_objective <- as.numeric(.first_non_null_local(
    optimizer_trace$deoptim_objective,
    if (is.list(de_fit)) de_fit$optim$bestval else NULL,
    NA_real_
  ))
  optimizer_local_objective <- as.numeric(.first_non_null_local(
    optimizer_trace$local_objective,
    if (is.list(local_fit)) local_fit$value else NULL,
    NA_real_
  ))
  optimizer_local_attempted <- isTRUE(.first_non_null_local(optimizer_trace$local_attempted, is.list(local_fit)))
  optimizer_local_accepted <- isTRUE(.first_non_null_local(optimizer_trace$local_accepted, FALSE))
  optimizer_local_convergence <- as.integer(.first_non_null_local(
    optimizer_trace$local_convergence,
    if (is.list(local_fit)) local_fit$convergence else NULL,
    NA_integer_
  ))
  optimizer_local_maxit <- as.integer(.first_non_null_local(optimizer_trace$local_maxit, NA_integer_))

  summary_df <- data.frame(
    metric = c(
      "fit_mode",
      "optimizer_method",
      "optimizer_deoptim_objective",
      "optimizer_local_objective",
      "optimizer_local_attempted",
      "optimizer_local_accepted",
      "optimizer_local_convergence",
      "optimizer_local_maxit",
      "objective",
      "objective_invivo",
      "objective_invitro",
      "joint_weight_invivo",
      "joint_weight_invitro",
      "joint_invitro_growth_weight",
      "joint_invitro_ploidy_weight",
      "joint_invitro_flow_weight",
      "joint_restriction",
      "glucose",
      "seed",
      "itermax",
      "NP_requested",
      "NP_used",
      "joint_np_min_factor",
      "n_cores_requested",
      "n_cores_used",
      "n_parameters",
      "n_invivo_scenarios",
      "joint_init_candidates_tsv",
      "n_joint_init_candidates",
      "invitro_forced_glucose"
    ),
    value = c(
      "fit_joint",
      optimizer_method,
      as.character(optimizer_deoptim_objective),
      as.character(optimizer_local_objective),
      as.character(optimizer_local_attempted),
      as.character(optimizer_local_accepted),
      as.character(optimizer_local_convergence),
      as.character(optimizer_local_maxit),
      as.character(best_comp$objective),
      as.character(best_comp$invivo$L),
      as.character(best_comp$invitro$objective),
      as.character(ctx$joint_weight_invivo),
      as.character(ctx$joint_weight_invitro),
      as.character(ctx$joint_invitro_growth_weight),
      as.character(ctx$joint_invitro_ploidy_weight),
      as.character(ctx$joint_invitro_flow_weight),
      as.character(ctx$joint_restriction),
      as.character(ctx$invivo$cfg$glucose),
      as.character(ctx$seed),
      as.character(ctx$itermax),
      as.character(ctx$NP),
      as.character(max(ctx$NP, ctx$joint_np_min_factor * length(ctx$init))),
      as.character(ctx$joint_np_min_factor),
      as.character(ctx$n_cores_requested),
      as.character(ctx$n_cores_used),
      as.character(length(ctx$init)),
      as.character(length(ctx$invivo$scenarios)),
      as.character(.first_non_null_local(ctx$joint_init_candidates_path, NA_character_)),
      as.character(.first_non_null_local(ctx$n_joint_init_candidates, 0L)),
      "TRUE"
    ),
    stringsAsFactors = FALSE
  )
  constraint_summary_df <- joint_constraint_component_df(best_comp$constraint_metrics)
  if (nrow(constraint_summary_df)) {
    summary_df <- dplyr::bind_rows(
      summary_df,
      data.frame(
        metric = constraint_summary_df$component,
        value = constraint_summary_df$value,
        stringsAsFactors = FALSE
      )
    )
  }
  write.table(summary_df, file = file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  file.copy(ctx$invivo$cfg$parameter_table, file.path(out_dir, "parameter_table_input_invivo.csv"), overwrite = TRUE)
  file.copy(ctx$invitro$parameter_table, file.path(out_dir, "parameter_table_input_invitro.csv"), overwrite = TRUE)
  write.table(
    ctx$invivo$param_bundle$transformed_output,
    file = file.path(out_dir, "parameter_table_invivo_transformed.csv"),
    sep = ",",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    ctx$invitro$spec,
    file = file.path(out_dir, "parameter_table_invitro_transformed.csv"),
    sep = ",",
    quote = FALSE,
    row.names = FALSE
  )
  saveRDS(
    list(
      ctx = list(
        invivo_cfg = INVIVO_ENV$sanitize_cfg_for_persistence(ctx$invivo$cfg),
        invitro_cfg = ctx$invitro$cfg,
        ivt_extra_names = ctx$ivt_extra_names,
        joint_weight_invivo = ctx$joint_weight_invivo,
        joint_weight_invitro = ctx$joint_weight_invitro,
        joint_invitro_growth_weight = ctx$joint_invitro_growth_weight,
        joint_invitro_ploidy_weight = ctx$joint_invitro_ploidy_weight,
        joint_invitro_flow_weight = ctx$joint_invitro_flow_weight,
        joint_restriction = ctx$joint_restriction,
        joint_require_biological_constraints = ctx$joint_require_biological_constraints,
        joint_constraint_penalty = ctx$joint_constraint_penalty,
        joint_require_invivo_pred1000_ploidy_gt2 = ctx$joint_require_invivo_pred1000_ploidy_gt2,
        joint_invivo_ploidy_horizon_day = ctx$joint_invivo_ploidy_horizon_day,
        joint_invivo_min_ploidy_fold = ctx$joint_invivo_min_ploidy_fold,
        joint_require_invitro_growth_nonnegative = ctx$joint_require_invitro_growth_nonnegative,
        joint_require_invitro_ploidy_phenotype = ctx$joint_require_invitro_ploidy_phenotype,
        joint_invitro_2N_wgd_min_N = ctx$joint_invitro_2N_wgd_min_N,
        joint_invitro_2N_wgd_min_fraction = ctx$joint_invitro_2N_wgd_min_fraction,
        joint_invitro_4N_min_chr_drop = ctx$joint_invitro_4N_min_chr_drop,
        n_cores_requested = ctx$n_cores_requested,
        n_cores_used = ctx$n_cores_used,
        optimizer_trace = optimizer_trace
      ),
      deoptim = de_fit,
      local_optim = local_fit,
      best_components = best_comp
    ),
    file = file.path(out_dir, "fit_result.rds")
  )
  saveRDS(
    INVIVO_ENV$sanitize_cfg_for_persistence(ctx$invivo$cfg),
    file = file.path(out_dir, "fit_config.rds")
  )
}

validate_fit_joint_inputs <- function(argv) {
  cfg_raw <- resolve_joint_raw_config(argv)
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg_raw$glucose, TRUE),
    default = TRUE
  ))

  invivo_parameter_table <- trim_cli_scalar_local(cfg_raw$parameter_table)
  if (is.null(invivo_parameter_table)) {
    invivo_parameter_table <- default_o2g_parameter_table_path_common(
      script_dir = SCRIPT_DIR,
      glucose = glucose_use,
      must_exist = TRUE
    )
  }
  if (!file.exists(invivo_parameter_table)) {
    stop("Joint fit in vivo parameter table not found: ", invivo_parameter_table, call. = FALSE)
  }
  INVIVO_ENV$build_transformed_parameter_table(
    path = invivo_parameter_table,
    fit_treatment = isTRUE(as_bool(cfg_raw$fit_treatment, FALSE)),
    fit_tau_O2 = isTRUE(as_bool(cfg_raw$fit_tau_O2, FALSE)),
    O2_growth = isTRUE(as_bool(cfg_raw$O2_growth, TRUE)),
    glucose = glucose_use,
    harvest_init_multiplier = isTRUE(as_bool(cfg_raw$harvest_init_multiplier, FALSE))
  )

  data_dir <- path_or_default(cfg_raw$data_dir, file.path(OXYGEN_ROOT, "data", "InVivoData_Gemcitabine"))
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  if (!file.exists(dt_path)) {
    stop("Joint fit missing in vivo tumor-burden workbook: ", dt_path, call. = FALSE)
  }
  INVIVO_ENV$resolve_terminal_ploidy_path(data_dir)

  invitro_parameter_table <- trim_cli_scalar_local(.first_non_null_local(
    cfg_raw$invitro_parameter_table,
    cfg_raw$parameter_table_invitro
  ))
  if (is.null(invitro_parameter_table)) {
    invitro_parameter_table <- INVITRO_ENV$default_parameter_table_path(
      must_exist = TRUE
    )
  }
  INVITRO_ENV$validate_invitro_parameter_table(
    parameter_table = invitro_parameter_table,
    dt = as_num(.first_non_null_local(cfg_raw$invitro_dt, cfg_raw$dt, cfg_raw$DT), 0.1),
    init_total_size = as_num(.first_non_null_local(cfg_raw$invitro_init_total_size, cfg_raw$init_total_size), 1e6),
    o2_upper_bound = as_num(cfg_raw$invitro_o2_upper_bound, 21),
    fixed_oxygen = TRUE
  )

  fit_objects_dir <- path_or_default(
    cfg_raw$fit_objects_dir,
    INVITRO_ENV$default_fit_objects_dir(must_exist = TRUE)
  )
  flow_density_path <- INVITRO_ENV$resolve_optional_flow_density_path(cfg_raw$flow_density_path)
  INVITRO_ENV$validate_invitro_fit_objects(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path
  )
  invisible(TRUE)
}

main_fit_seed_joint <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  ctx <- build_joint_context(argv)
  set.seed(ctx$seed)
  out_dir <- ctx$out_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  objective_value <- function(par_t) {
    tryCatch(
      joint_objective_components(par_t, ctx)$objective,
      error = function(e) {
        if (isTRUE(as_bool(ctx$raw$trace_obj, FALSE))) {
          message("[fit_joint] objective error: ", conditionMessage(e))
        }
        1e9
      }
    )
  }

  NP_use <- max(ctx$NP, ctx$joint_np_min_factor * length(ctx$init))
  joint_init <- read_joint_init_candidates(ctx$raw$joint_init_candidates_tsv, ctx)
  if (!is.null(joint_init)) {
    NP_use <- max(NP_use, nrow(joint_init$matrix))
    ctx$joint_init_candidates_path <- joint_init$path
    ctx$n_joint_init_candidates <- nrow(joint_init$matrix)
    write.table(
      joint_init$used,
      file = file.path(out_dir, "joint_init_candidates_used.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    candidate_scores <- score_joint_init_candidates(joint_init$matrix, ctx)
    write.table(
      candidate_scores,
      file = file.path(out_dir, "joint_init_candidate_scores.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message(
      "[fit_joint] Using joint init candidates: n=", nrow(joint_init$matrix),
      ", path=", joint_init$path
    )
  } else {
    ctx$joint_init_candidates_path <- NA_character_
    ctx$n_joint_init_candidates <- 0L
  }
  de_ctrl <- list(
    trace = TRUE,
    NP = NP_use,
    itermax = max(1L, ctx$itermax)
  )
  de_cluster <- NULL
  de_active_cores <- 1L
  if (ctx$n_cores_requested > 1L) {
    message("[fit_joint] DEoptim parallel requested with n_cores=", ctx$n_cores_requested, ".")
    de_cluster <- start_joint_deoptim_cluster(ctx$n_cores_requested)
    on.exit(try(parallel::stopCluster(de_cluster), silent = TRUE), add = TRUE)
    de_ctrl$cluster <- de_cluster
    de_active_cores <- length(de_cluster)
    message("[fit_joint] DEoptim parallel enabled: workers=", de_active_cores, ".")
  } else {
    de_ctrl$parallelType <- "none"
    message("[fit_joint] DEoptim running in serial mode (n_cores=1).")
  }
  ctx$n_cores_used <- de_active_cores
  message(
    "[fit_joint] Starting DEoptim: parameters=", length(ctx$init),
    ", NP=", NP_use,
    ", itermax=", ctx$itermax,
    ", NP_min_factor=", ctx$joint_np_min_factor,
    ", n_cores=", ctx$n_cores_used,
    ", weights(invivo,invitro)=(",
    signif(ctx$joint_weight_invivo, 6), ", ",
    signif(ctx$joint_weight_invitro, 6), ")"
  )
  if (!is.null(joint_init)) {
    de_ctrl$initialpop <- build_joint_de_initialpop(
      np = NP_use,
      lower = ctx$lower,
      upper = ctx$upper,
      candidates = joint_init$matrix
    )
  }
  de_fit <- DEoptim::DEoptim(
    fn = objective_value,
    lower = ctx$lower,
    upper = ctx$upper,
    control = de_ctrl
  )
  best_t <- as.numeric(de_fit$optim$bestmem)
  names(best_t) <- names(ctx$init)
  de_best_objective <- suppressWarnings(as.numeric(de_fit$optim$bestval))
  local_maxit <- as_int(.first_non_null_local(ctx$raw$local_optim_maxit, ctx$raw$optim_maxit, 200L), 200L)
  if (!is.finite(local_maxit) || is.na(local_maxit) || local_maxit < 1L) local_maxit <- 200L
  local_attempted <- FALSE
  local_accepted <- FALSE
  local_fit <- NULL
  local_best_objective <- NA_real_
  local_convergence <- NA_integer_
  local_message <- NA_character_
  if (is.finite(de_best_objective)) {
    local_attempted <- TRUE
    message("[fit_joint] Starting L-BFGS-B local refinement from DEoptim best; maxit=", local_maxit, ".")
    local_fit <- tryCatch(
      suppressWarnings(
        optim(
          par = best_t,
          fn = objective_value,
          method = "L-BFGS-B",
          lower = ctx$lower,
          upper = ctx$upper,
          control = list(maxit = local_maxit)
        )
      ),
      error = function(e) {
        warning("[fit_joint] L-BFGS-B local refinement failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.list(local_fit)) {
      local_best_objective <- suppressWarnings(as.numeric(local_fit$value))
      local_convergence <- suppressWarnings(as.integer(local_fit$convergence))
      local_message <- as.character(.first_non_null_local(local_fit$message, NA_character_))
      if (is.finite(local_best_objective) && local_best_objective < de_best_objective) {
        best_t <- as.numeric(local_fit$par)
        names(best_t) <- names(ctx$init)
        best_t <- clip(best_t, ctx$lower, ctx$upper)
        local_accepted <- TRUE
        message(
          "[fit_joint] L-BFGS-B local refinement improved objective: ",
          signif(de_best_objective, 8), " -> ", signif(local_best_objective, 8), "."
        )
      } else {
        message("[fit_joint] L-BFGS-B local refinement did not improve objective; keeping DEoptim best.")
      }
    }
  }
  de_method <- if (isTRUE(ctx$n_cores_used > 1L)) "DEoptim_parallel" else "DEoptim_serial"
  optimizer_method <- if (isTRUE(local_attempted)) paste0(de_method, "_plus_LBFGSB_serial") else de_method
  optimizer_trace <- list(
    method = optimizer_method,
    deoptim_objective = as.numeric(de_best_objective),
    local_objective = as.numeric(local_best_objective),
    local_attempted = isTRUE(local_attempted),
    local_accepted = isTRUE(local_accepted),
    local_convergence = as.integer(local_convergence),
    local_maxit = as.integer(local_maxit),
    local_message = local_message
  )
  best_comp <- joint_objective_components(best_t, ctx)
  write_joint_outputs(
    best_t,
    best_comp,
    ctx,
    out_dir = out_dir,
    de_fit = de_fit,
    local_fit = local_fit,
    optimizer_trace = optimizer_trace
  )
  message("Done. Joint results written to: ", normalizePath(out_dir, mustWork = FALSE))
  message("Best joint objective: ", signif(best_comp$objective, 6))
}

main_run_from_config_joint <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  parsed <- INVIVO_ENV$.runner_resolve_config(
    argv = argv,
    script_dir = SCRIPT_DIR,
    caller_wd = getwd()
  )
  cfg <- parsed$cfg

  require_nonempty <- function(key) {
    val <- INVIVO_ENV$.runner_cli_string(cfg[[key]])
    if (is.null(val) || !nzchar(trimws(val))) {
      stop("Config key '", key, "' is required for --fit_joint --mode=run.", call. = FALSE)
    }
    val
  }

  require_nonempty("out_root")
  require_nonempty("data_dir")

  seed_plan <- INVIVO_ENV$.runner_resolve_seeds(parsed)
  append_ts <- as_bool(cfg$append_run_prefix_timestamp, FALSE)
  ts_format <- INVIVO_ENV$.runner_cli_string(.first_non_null_local(cfg$run_prefix_timestamp_format, "%Y%m%d_%H%M%S"))
  run_prefix <- INVIVO_ENV$.runner_cli_string(.first_non_null_local(cfg$run_prefix, "fit_joint_O2G_supply_demand_MAP"))
  if (is.null(run_prefix) || !nzchar(trimws(run_prefix))) {
    run_prefix <- "fit_joint_O2G_supply_demand_MAP"
  }
  if (append_ts) {
    run_prefix <- paste0(run_prefix, "_", format(Sys.time(), ts_format))
  }

  out_root <- normalizePath(INVIVO_ENV$.runner_cli_string(cfg$out_root), mustWork = FALSE)
  data_dir <- normalizePath(INVIVO_ENV$.runner_cli_string(cfg$data_dir), mustWork = FALSE)
  auto_viz <- as_bool(cfg$auto_viz, TRUE)
  viz_report_dt <- as_num(cfg$viz_report_dt, 1)
  viz_top_n <- as_int(cfg$viz_top_n, 6L)

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(out_root, run_prefix)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  run_log <- file.path(run_dir, "run_status.log")
  if (!file.exists(run_log)) file.create(run_log)
  log_line <- function(...) {
    line <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    cat(line, "\n", sep = "")
    cat(line, "\n", sep = "", file = run_log, append = TRUE)
  }

  fit_script <- normalizePath(file.path(WORKFLOW_ROOT, "optimizer", "fit_model_O2G_supply_demand_MAP.R"), mustWork = FALSE)
  viz_script <- normalizePath(file.path(WORKFLOW_ROOT, "vis", "viz_invivo_model_O2G_supply_demand_MAP_results.R"), mustWork = FALSE)
  invitro_viz_script <- normalizePath(file.path(WORKFLOW_ROOT, "vis", "viz_invitro_model_O2G_supply_demand_MAP_results.R"), mustWork = FALSE)
  report_script <- normalizePath(file.path(WORKFLOW_ROOT, "report", "render_fit_report.R"), mustWork = FALSE)
  fit_base <- INVIVO_ENV$.runner_build_fit_base_args(cfg)
  snapshots <- INVIVO_ENV$.runner_write_config_snapshots(
    run_dir = run_dir,
    config_path = parsed$config_path,
    resolved_cfg = cfg,
    fit_arg_values = fit_base$values,
    seed_plan = seed_plan,
    fit_script = fit_script,
    viz_script = viz_script,
    report_script = report_script
  )

  has_parameter_snapshot <- isTRUE(INVIVO_ENV$.runner_copy_parameter_table_snapshot(cfg$parameter_table, run_dir))

  log_line("Running O2G_supply_demand_MAP joint fit")
  log_line("Config: ", parsed$config_path)
  log_line("Config input snapshot: ", snapshots$input)
  log_line("Config resolved snapshot: ", snapshots$resolved)
  log_line("Fit script: ", fit_script)
  log_line("Viz script: ", viz_script)
  log_line("In vitro viz script: ", invitro_viz_script)
  log_line("Report script: ", report_script)
  log_line("Data dir: ", data_dir)
  log_line("Seeds: ", seed_plan$seeds_csv, " (", seed_plan$seed_source, ")")
  log_line("Run dir: ", run_dir)
  log_line("Run log: ", run_log)
  if (has_parameter_snapshot) {
    log_line("Parameter table input snapshot: ", file.path(run_dir, "parameter_table_input.csv"))
  } else {
    log_line("Parameter table snapshot: (missing --parameter_table or file not found)")
  }
  log_line("Run prefix timestamp suffix: ", append_ts, " (format=", ts_format, ")")
  log_line("Auto viz/report: ", auto_viz, " (report_dt=", viz_report_dt, ", top_n=", viz_top_n, ")")

  for (seed in seed_plan$seeds) {
    if (!is.finite(seed)) next
    seed <- as.integer(seed)
    seed_dir <- file.path(run_dir, paste0("seed", seed))
    dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)
    INVIVO_ENV$.runner_copy_parameter_table_snapshot(cfg$parameter_table, seed_dir)

    fit_log <- file.path(seed_dir, "fit_status.log")
    viz_log <- file.path(seed_dir, "viz_status.log")
    invitro_viz_log <- file.path(seed_dir, "invitro_viz_status.log")
    report_log <- file.path(seed_dir, "report_status.log")
    fit_args <- c(
      fit_script,
      "--fit_joint",
      paste0("--seed=", seed),
      paste0("--out_dir=", seed_dir),
      paste0("--data_dir=", data_dir),
      fit_base$args
    )

    log_line("seed=", seed, ": start")
    log_line("seed=", seed, ": fit_log=", fit_log)
    log_line("Fit command: Rscript ", paste(fit_args, collapse = " "))
    fit_status <- INVIVO_ENV$.runner_exec_to_log("Rscript", fit_args, fit_log, run_log_path = run_log)
    if (!identical(fit_status, 0L)) {
      INVIVO_ENV$.runner_stop_with_log_tail(paste0("seed=", seed, " fit"), fit_log, fit_status)
    }
    log_line("seed=", seed, ": fit done")

    if (auto_viz) {
      viz_args <- c(
        viz_script,
        paste0("--fit_dir=", seed_dir),
        paste0("--out_dir=", file.path(seed_dir, "viz", "invivo")),
        paste0("--data_dir=", data_dir),
        paste0("--report_dt=", viz_report_dt),
        paste0("--top_n=", viz_top_n),
        "--n_cores=1"
      )
      log_line("seed=", seed, ": viz start")
      log_line("seed=", seed, ": viz_log=", viz_log)
      log_line("Viz command: Rscript ", paste(viz_args, collapse = " "))
      viz_status <- INVIVO_ENV$.runner_exec_to_log("Rscript", viz_args, viz_log, run_log_path = run_log)
      if (!identical(viz_status, 0L)) {
        INVIVO_ENV$.runner_stop_with_log_tail(paste0("seed=", seed, " viz"), viz_log, viz_status)
      }
      log_line("seed=", seed, ": viz done")

      invitro_viz_args <- c(
        invitro_viz_script,
        paste0("--fit_dir=", seed_dir),
        paste0("--out_dir=", file.path(seed_dir, "viz", "invitro"))
      )
      log_line("seed=", seed, ": in vitro viz start")
      log_line("seed=", seed, ": invitro_viz_log=", invitro_viz_log)
      log_line("In vitro viz command: Rscript ", paste(invitro_viz_args, collapse = " "))
      invitro_viz_status <- INVIVO_ENV$.runner_exec_to_log("Rscript", invitro_viz_args, invitro_viz_log, run_log_path = run_log)
      if (!identical(invitro_viz_status, 0L)) {
        INVIVO_ENV$.runner_stop_with_log_tail(paste0("seed=", seed, " in vitro viz"), invitro_viz_log, invitro_viz_status)
      }
      log_line("seed=", seed, ": in vitro viz done")

      report_args <- c(
        report_script,
        paste0("--fit_dir=", seed_dir)
      )
      log_line("seed=", seed, ": report start")
      log_line("seed=", seed, ": report_log=", report_log)
      log_line("Report command: Rscript ", paste(report_args, collapse = " "))
      report_status <- INVIVO_ENV$.runner_exec_to_log("Rscript", report_args, report_log, run_log_path = run_log)
      if (!identical(report_status, 0L)) {
        INVIVO_ENV$.runner_stop_with_log_tail(paste0("seed=", seed, " report"), report_log, report_status)
      }
      log_line("seed=", seed, ": report done")
    }
  }

  message("All done. Joint run directory: ", normalizePath(run_dir, mustWork = FALSE))
  invisible(normalizePath(run_dir, mustWork = FALSE))
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  mode_raw <- trim_cli_scalar_local(argv$mode)
  if (is.null(mode_raw) || !nzchar(mode_raw)) {
    return(main_fit_seed_joint(argv))
  }
  mode_use <- tolower(mode_raw)
  if (mode_use %in% c("run", "runner")) {
    return(main_run_from_config_joint(argv))
  }
  if (mode_use %in% c("fit_seed", "fit", "single_seed")) {
    argv$mode <- NULL
    return(main_fit_seed_joint(argv))
  }
  stop("Unsupported --mode for --fit_joint: ", mode_use, ". Expected one of: run, fit_seed.", call. = FALSE)
}

if (sys.nframe() == 0) {
  main()
}
