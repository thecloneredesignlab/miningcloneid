#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEoptim))
suppressPackageStartupMessages(library(dplyr))

.o2_joint_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2_joint_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = FALSE)

source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2_joint_bootstrap_script_dir)

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
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_invivo_backend.R")
)
INVITRO_ENV <- load_backend_env_local(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_invitro_backend.R")
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
  run_prefix <- path_or_default(cfg_raw$run_prefix, "fit_joint_O2_supply_demand_MAP")
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

joint_deoptim_control <- function(ctx, NP_use) {
  itermax <- as.integer(.first_non_null_local(ctx$itermax, 1L))
  if (!is.finite(itermax) || is.na(itermax) || itermax < 1L) itermax <- 1L
  de_reltol <- as.numeric(.first_non_null_local(ctx$de_reltol, 1e-4))
  if (!is.finite(de_reltol) || is.na(de_reltol) || de_reltol <= 0) de_reltol <- 1e-4
  de_steptol <- as.integer(.first_non_null_local(ctx$de_steptol, 25L))
  if (!is.finite(de_steptol) || is.na(de_steptol) || de_steptol < 1L) de_steptol <- 25L
  list(
    trace = TRUE,
    NP = NP_use,
    itermax = itermax,
    initialpop = joint_deoptim_initial_population(ctx, NP_use),
    reltol = de_reltol,
    steptol = de_steptol
  )
}

joint_deoptim_iter_completed <- function(de_fit, iter_target = NA_integer_) {
  iter_completed <- as.integer(.first_non_null_local(
    if (is.list(de_fit) && !is.null(de_fit$optim$iter)) as.integer(de_fit$optim$iter) else NULL,
    if (is.list(de_fit) &&
        !is.null(de_fit$member$bestmemit) &&
        (is.matrix(de_fit$member$bestmemit) || is.data.frame(de_fit$member$bestmemit))) {
      as.integer(nrow(de_fit$member$bestmemit))
    } else {
      NULL
    },
    if (is.list(de_fit) && !is.null(de_fit$member$bestvalit)) {
      as.integer(length(de_fit$member$bestvalit))
    } else {
      NULL
    },
    iter_target,
    NA_integer_
  ))
  if (!is.finite(iter_completed) || is.na(iter_completed) || iter_completed < 0L) {
    return(NA_integer_)
  }
  iter_target <- as.integer(iter_target)
  if (is.finite(iter_target) && !is.na(iter_target) && iter_completed > iter_target) {
    iter_completed <- iter_target
  }
  iter_completed
}

joint_deoptim_stop_reason <- function(iter_completed, iter_target, interrupted = FALSE) {
  if (isTRUE(interrupted)) {
    return("user_interrupt")
  }
  iter_completed <- as.integer(iter_completed)
  iter_target <- as.integer(iter_target)
  if (is.finite(iter_completed) && is.finite(iter_target) && iter_completed < iter_target) {
    return("early_stop_reltol_or_steptol")
  }
  if (is.finite(iter_completed) && is.finite(iter_target) && iter_completed >= iter_target) {
    return("itermax_reached")
  }
  NA_character_
}

resolve_joint_raw_config <- function(argv) {
  parsed <- INVIVO_ENV$.runner_resolve_config(
    argv = argv,
    script_dir = SCRIPT_DIR,
    caller_wd = getwd()
  )
  cfg <- parsed$cfg
  cfg$.config_dir <- parsed$config_dir
  cfg$.caller_wd <- parsed$caller_wd
  cfg
}

build_joint_invivo_context <- function(cfg_raw) {
  parameter_table <- trim_cli_scalar_local(cfg_raw$parameter_table)
  if (is.null(parameter_table)) {
    parameter_table <- default_o2_parameter_table_path_common(
      script_dir = SCRIPT_DIR,
      must_exist = TRUE
    )
  }
  if (!file.exists(parameter_table)) {
    stop("Joint fit in vivo parameter table not found: ", parameter_table)
  }

  model_path <- file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R")
  if (!file.exists(model_path)) stop("Cannot find model_O2_supply_demand_MAP.R at ", model_path)
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
    use_necrosis_loss = as_bool(cfg_raw$use_necrosis_loss, FALSE),
    necrosis_mapping_csv = trim_cli_scalar_local(cfg_raw$necrosis_mapping_csv),
    sigma_necrosis_logit = as_num(cfg_raw$sigma_necrosis_logit, 0.75),
    lambda_necrosis = as_num(cfg_raw$lambda_necrosis, 1.0),
    necrosis_fraction_eps = as_num(cfg_raw$necrosis_fraction_eps, 1e-4),
    use_soft_prior = as_bool(cfg_raw$use_soft_prior, TRUE),
    lambda_prior = as_num(cfg_raw$lambda_prior, 0.1),
    harvest_init_multiplier = as_bool(cfg_raw$harvest_init_multiplier, FALSE),
    prior_center_log_init_mult = as_num(cfg_raw$prior_center_log_init_mult, 0.0),
    prior_sd_log_init_mult = as_num(cfg_raw$prior_sd_log_init_mult, 0.35),
    log_init_mult_lower = as_num(cfg_raw$log_init_mult_lower, -1.0),
    log_init_mult_upper = as_num(cfg_raw$log_init_mult_upper, 1.0),
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

shared_invitro_param_names <- function() {
  loss_shared <- c("buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp")
  growth_shared <- c("log10_lam_max")
  out <- c(
    growth_shared, "log10_p_misseg",
    "log10_p_mis_base", "log10_k_o_mis", loss_shared, "log10_alpha_o2", "gamma_growth",
    "log10_mu_hp", "gamma_mu", "log10_O2_crit", "n_O", "log10_p_wgd"
  )
  out
}

joint_shared_natural_param_names <- function() {
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

joint_parse_soft_param_list <- function(x, default = character()) {
  if (is.null(x) || !length(x)) return(unique(as.character(default)))
  vals <- unlist(x, recursive = TRUE, use.names = FALSE)
  vals <- as.character(vals)
  vals <- unlist(strsplit(vals, "[,;]", perl = TRUE), use.names = FALSE)
  vals <- trimws(vals)
  unique(vals[nzchar(vals)])
}

joint_soft_split_delta_name <- function(center_name) {
  paste0("delta__", as.character(center_name))
}

joint_default_soft_coupling_params <- function() {
  c(
    "O2_crit", "mu_hp", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O",
    "alpha_o2", "gamma_growth", "lam_max", "p_mis_base",
    "p_wgd", "gamma_mu"
  )
}

joint_soft_split_natural_param_names <- function(cfg_raw) {
  if (!isTRUE(as_bool(cfg_raw$joint_soft_coupling_enable, TRUE))) {
    return(character(0))
  }
  default_params <- joint_default_soft_coupling_params()
  params <- joint_parse_soft_param_list(cfg_raw$joint_soft_coupling_params, default = default_params)
  shared <- joint_shared_natural_param_names()
  bad <- setdiff(params, shared)
  if (length(bad) > 0L) {
    stop(
      "joint_soft_coupling_params contains parameters that are not joint-shared: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  params
}

joint_soft_split_transformed_name_map <- function(split_params) {
  out <- vapply(split_params, invitro_shared_param_name_for_natural, character(1))
  names(out) <- split_params
  out
}

joint_transformed_to_natural <- function(value, symbol) {
  value <- as.numeric(value)
  if (!is.finite(value)) return(NA_real_)
  if (symbol %in% c("buffer_smax", "gamma_growth", "gamma_mu", "n_O")) {
    return(value)
  }
  10^value
}

joint_soft_probability_parameter <- function(symbol) {
  as.character(symbol) %in% c("p_mis_base", "p_misseg", "p_wgd")
}

joint_safe_logit <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p) & p > 0 & p < 1
  out[ok] <- log(p[ok] / (1 - p[ok]))
  out
}

joint_soft_coupling_metadata <- function(split_params,
                                         joint_bounds,
                                         invivo,
                                         invitro,
                                         cfg_raw) {
  if (!length(split_params)) return(data.frame())
  invivo_specs <- INVIVO_ENV$parameter_table_specs()
  bounds_tab <- joint_bounds$summary
  sigma_default <- as_num(cfg_raw$joint_soft_coupling_sigma_default, 1.5)
  if (!is.finite(sigma_default) || sigma_default <= 0) {
    stop("joint_soft_coupling_sigma_default must be > 0.", call. = FALSE)
  }

  rows <- lapply(split_params, function(symbol) {
    bound_row <- bounds_tab[as.character(bounds_tab$parameter) == as.character(symbol), , drop = FALSE]
    if (nrow(bound_row) != 1L) {
      stop("Missing merged joint-bound row for soft-coupled parameter: ", symbol, call. = FALSE)
    }
    spec_row <- invivo_specs[as.character(invivo_specs$param_symbol) == as.character(symbol), , drop = FALSE]
    if (nrow(spec_row) != 1L) {
      stop("Missing in vivo parameter-table spec for soft-coupled parameter: ", symbol, call. = FALSE)
    }
    center_name <- as.character(spec_row$param_name[[1]])
    invitro_name <- invitro_shared_param_name_for_natural(symbol)
    if (is.null(invitro_name) || !nzchar(invitro_name)) {
      stop("No in vitro transformed name is known for soft-coupled parameter: ", symbol, call. = FALSE)
    }
    if (!identical(center_name, invitro_name)) {
      stop(
        "Soft coupling currently requires matching transformed names between backends. ",
        symbol, " maps to in vivo '", center_name, "' but in vitro '", invitro_name, "'.",
        call. = FALSE
      )
    }
    if (!(center_name %in% names(joint_bounds$invivo_optimizer$init))) {
      stop("Soft-coupled center parameter missing from joint optimizer: ", center_name, call. = FALSE)
    }
    if (!(invitro_name %in% as.character(invitro$spec$param_name))) {
      stop("Soft-coupled in vitro parameter missing from in vitro optimizer spec: ", invitro_name, call. = FALSE)
    }

    transform <- as.character(spec_row$transform[[1]])
    joint_union_lower_t <- as.numeric(bound_row$joint_union_lower_t[[1]])
    joint_union_upper_t <- as.numeric(bound_row$joint_union_upper_t[[1]])
    if (!is.finite(joint_union_lower_t) ||
        !is.finite(joint_union_upper_t) ||
        joint_union_lower_t > joint_union_upper_t) {
      stop("Invalid joint union transformed bounds for soft-coupled parameter: ", symbol, call. = FALSE)
    }

    sigma_key <- paste0("joint_soft_coupling_sigma_", symbol)
    sigma <- as_num(cfg_raw[[sigma_key]], sigma_default)
    if (!is.finite(sigma) || sigma <= 0) {
      stop(sigma_key, " must be > 0.", call. = FALSE)
    }
    center_init_t <- as.numeric(joint_bounds$invivo_optimizer$init[[center_name]])
    if (!is.finite(center_init_t) ||
        center_init_t < joint_union_lower_t ||
        center_init_t > joint_union_upper_t) {
      stop("Soft-coupled center init is outside the joint union transformed bounds for: ", symbol, call. = FALSE)
    }
    delta_abs <- joint_union_upper_t - joint_union_lower_t

    data.frame(
      parameter = symbol,
      center_name = center_name,
      invitro_name = invitro_name,
      delta_name = joint_soft_split_delta_name(center_name),
      transform = transform,
      center_init_t = center_init_t,
      joint_union_lower_t = joint_union_lower_t,
      joint_union_upper_t = joint_union_upper_t,
      delta_lower_t = -delta_abs,
      delta_upper_t = delta_abs,
      sigma_delta = sigma,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

joint_is_absolute_path <- function(path) {
  path <- as.character(path[[1]])
  grepl("^/", path) || grepl("^[A-Za-z]:[\\/]", path)
}

joint_resolve_path_local <- function(path, cfg_raw, must_exist = FALSE) {
  path <- trim_cli_scalar_local(path)
  if (is.null(path)) return(NULL)
  candidates <- if (joint_is_absolute_path(path)) {
    path
  } else {
    unique(c(
      file.path(trim_cli_scalar_local(cfg_raw$.config_dir), path),
      file.path(trim_cli_scalar_local(cfg_raw$.caller_wd), path),
      file.path(getwd(), path),
      path
    ))
  }
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  out <- if (length(hit)) hit[[1]] else candidates[[1]]
  if (isTRUE(must_exist) && !file.exists(out)) {
    stop("Path not found: ", out, call. = FALSE)
  }
  normalizePath(out, mustWork = FALSE)
}

joint_soft_coupling_start_table_path <- function(cfg_raw, invivo_parameter_table) {
  explicit <- trim_cli_scalar_local(.first_non_null_local(
    cfg_raw$joint_soft_coupling_parameters_table,
    cfg_raw$joint_soft_coupling_parameters_table_path
  ))
  base_dir <- dirname(normalizePath(invivo_parameter_table, mustWork = FALSE))
  if (is.null(explicit)) {
    return(normalizePath(file.path(base_dir, "joint_soft_coupling_parameters_table.csv"), mustWork = FALSE))
  }
  path <- explicit
  if (!joint_is_absolute_path(path)) {
    candidate <- file.path(base_dir, path)
    if (file.exists(candidate)) path <- candidate
  }
  normalizePath(path, mustWork = FALSE)
}

joint_read_soft_coupling_start_table <- function(path, required = FALSE) {
  if (is.null(path) || !nzchar(as.character(path))) return(data.frame())
  if (!file.exists(path)) {
    if (isTRUE(required)) {
      stop("joint_soft_coupling_parameters_table not found: ", path, call. = FALSE)
    }
    return(data.frame())
  }
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  param_col <- intersect(c("param_name", "parameter", "name"), names(tab))
  value_col <- intersect(c("value", "init_value", "init_transformed", "transformed_value"), names(tab))
  scale_col <- intersect(c("scale", "value_scale"), names(tab))
  if (!length(param_col) || !length(value_col) || !length(scale_col)) {
    stop(
      "joint_soft_coupling_parameters_table must contain columns: param_name, value, scale.",
      call. = FALSE
    )
  }
  out <- data.frame(
    param_name = trimws(as.character(tab[[param_col[[1]]]])),
    value = suppressWarnings(as.numeric(tab[[value_col[[1]]]])),
    scale = tolower(trimws(as.character(tab[[scale_col[[1]]]]))),
    stringsAsFactors = FALSE
  )
  bad_name <- !nzchar(out$param_name)
  if (any(bad_name)) stop("Empty param_name in joint_soft_coupling_parameters_table.", call. = FALSE)
  bad_value <- !is.finite(out$value)
  if (any(bad_value)) {
    stop(
      "Non-finite value in joint_soft_coupling_parameters_table for: ",
      paste(out$param_name[bad_value], collapse = ", "),
      call. = FALSE
    )
  }
  out$scale[out$scale %in% c("optimizer", "optim", "optimizer_scale", "transformed_scale")] <- "transformed"
  allowed <- c("transformed", "natural", "log10", "identity")
  bad_scale <- !(out$scale %in% allowed)
  if (any(bad_scale)) {
    stop(
      "Unsupported scale in joint_soft_coupling_parameters_table for: ",
      paste(out$param_name[bad_scale], collapse = ", "),
      ". Allowed scales: ", paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  dup <- out$param_name[duplicated(out$param_name)]
  if (length(dup)) {
    stop(
      "Duplicate param_name in joint_soft_coupling_parameters_table: ",
      paste(unique(dup), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

joint_optimizer_param_metadata <- function(optimizer_name) {
  optimizer_name <- as.character(optimizer_name)
  base_name <- sub("^ivt__", "", optimizer_name)
  invivo_specs <- INVIVO_ENV$parameter_table_specs()
  spec_row <- invivo_specs[as.character(invivo_specs$param_name) == base_name, , drop = FALSE]
  if (nrow(spec_row) == 1L) {
    return(list(
      parameter = as.character(spec_row$param_symbol[[1]]),
      transform = as.character(spec_row$transform[[1]])
    ))
  }
  if (startsWith(base_name, "log10_")) {
    return(list(parameter = sub("^log10_", "", base_name), transform = "log10"))
  }
  list(parameter = o2sd_parameter_natural_name(base_name), transform = "identity")
}

joint_soft_start_lookup <- function(soft_meta, optimizer_names = character()) {
  rows <- list()
  if (is.data.frame(soft_meta) && nrow(soft_meta)) {
    center_rows <- data.frame(
      input_name = as.character(soft_meta$center_name),
      optimizer_name = as.character(soft_meta$center_name),
      center_name = as.character(soft_meta$center_name),
      parameter = as.character(soft_meta$parameter),
      transform = as.character(soft_meta$transform),
      is_delta = FALSE,
      stringsAsFactors = FALSE
    )
    natural_alias_rows <- center_rows
    natural_alias_rows$input_name <- as.character(soft_meta$parameter)
    natural_alias_rows <- natural_alias_rows[natural_alias_rows$input_name != natural_alias_rows$optimizer_name, , drop = FALSE]
    delta_rows <- data.frame(
      input_name = as.character(soft_meta$delta_name),
      optimizer_name = as.character(soft_meta$delta_name),
      center_name = as.character(soft_meta$center_name),
      parameter = as.character(soft_meta$parameter),
      transform = as.character(soft_meta$transform),
      is_delta = TRUE,
      stringsAsFactors = FALSE
    )
    rows <- c(rows, list(center_rows, natural_alias_rows, delta_rows))
  }
  optimizer_names <- unique(as.character(optimizer_names))
  optimizer_names <- optimizer_names[nzchar(optimizer_names) & !startsWith(optimizer_names, "delta__")]
  if (length(optimizer_names)) {
    ordinary_rows <- lapply(optimizer_names, function(nm) {
      meta <- joint_optimizer_param_metadata(nm)
      data.frame(
        input_name = nm,
        optimizer_name = nm,
        center_name = NA_character_,
        parameter = as.character(meta$parameter),
        transform = as.character(meta$transform),
        is_delta = FALSE,
        stringsAsFactors = FALSE
      )
    })
    ordinary_rows <- dplyr::bind_rows(ordinary_rows)
    natural_alias_rows <- ordinary_rows
    natural_alias_rows$input_name <- as.character(ordinary_rows$parameter)
    natural_alias_rows <- natural_alias_rows[natural_alias_rows$input_name != natural_alias_rows$optimizer_name, , drop = FALSE]
    rows <- c(rows, list(ordinary_rows, natural_alias_rows))
  }
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(out)
  out <- out[!duplicated(out$input_name), , drop = FALSE]
  row.names(out) <- NULL
  out
}

joint_soft_delta_natural_to_optimizer <- function(delta_natural, center_t, transform, symbol) {
  delta_natural <- as.numeric(delta_natural)
  center_t <- as.numeric(center_t)
  if (!is.finite(delta_natural) || !is.finite(center_t)) return(NA_real_)
  if (identical(transform, "identity")) {
    return(delta_natural)
  }
  if (transform %in% c("log10", "log10_nonnegative")) {
    center_natural <- 10^center_t
    if (!is.finite(center_natural) || center_natural <= 0) return(NA_real_)
    ratio <- delta_natural / center_natural
    return(2 * asinh(ratio / 2) / log(10))
  }
  stop("scale='natural' is not supported for delta parameter ", symbol, " with transform ", transform, ".", call. = FALSE)
}

joint_soft_start_value_to_optimizer <- function(value, scale, lookup_row, center_t = NA_real_) {
  value <- as.numeric(value)
  scale <- as.character(scale)
  transform <- as.character(lookup_row$transform[[1]])
  symbol <- as.character(lookup_row$parameter[[1]])
  if (isTRUE(lookup_row$is_delta[[1]])) {
    if (identical(scale, "transformed")) {
      return(value)
    }
    if (identical(scale, "log10")) {
      if (!(transform %in% c("log10", "log10_nonnegative"))) {
        stop("scale='log10' is incompatible with delta parameter ", lookup_row$optimizer_name[[1]], " transform ", transform, ".", call. = FALSE)
      }
      return(value)
    }
    if (identical(scale, "identity")) {
      if (!identical(transform, "identity")) {
        stop("scale='identity' is incompatible with delta parameter ", lookup_row$optimizer_name[[1]], " transform ", transform, ".", call. = FALSE)
      }
      return(value)
    }
    if (identical(scale, "natural")) {
      return(joint_soft_delta_natural_to_optimizer(value, center_t, transform, symbol))
    }
    stop("Unsupported scale for delta parameter ", lookup_row$optimizer_name[[1]], ": ", scale, call. = FALSE)
  }
  if (identical(scale, "transformed")) {
    return(value)
  }
  if (identical(scale, "log10")) {
    if (!(transform %in% c("log10", "log10_nonnegative"))) {
      stop("scale='log10' is incompatible with parameter ", symbol, " transform ", transform, ".", call. = FALSE)
    }
    return(value)
  }
  if (identical(scale, "identity")) {
    if (!identical(transform, "identity")) {
      stop("scale='identity' is incompatible with parameter ", symbol, " transform ", transform, ".", call. = FALSE)
    }
    return(value)
  }
  if (identical(scale, "natural")) {
    return(INVIVO_ENV$transform_param_slot(value, transform, symbol, "joint_soft_coupling_parameters_table value"))
  }
  stop("Unsupported scale: ", scale, call. = FALSE)
}

joint_check_soft_reconstruction_feasibility <- function(par_t, soft_meta, tol = 1e-10) {
  if (!is.data.frame(soft_meta) || !nrow(soft_meta)) {
    return(list(feasible = TRUE, parameters = character(0)))
  }
  par_names <- names(par_t)
  par_t <- as.numeric(par_t)
  names(par_t) <- par_names
  bad <- character(0)
  for (i in seq_len(nrow(soft_meta))) {
    center_name <- as.character(soft_meta$center_name[[i]])
    delta_name <- as.character(soft_meta$delta_name[[i]])
    if (!(center_name %in% names(par_t)) || !(delta_name %in% names(par_t))) {
      bad <- c(bad, as.character(soft_meta$parameter[[i]]))
      next
    }
    center <- as.numeric(par_t[[center_name]])
    delta <- as.numeric(par_t[[delta_name]])
    lower <- as.numeric(soft_meta$joint_union_lower_t[[i]])
    upper <- as.numeric(soft_meta$joint_union_upper_t[[i]])
    vivo <- center + delta / 2
    vitro <- center - delta / 2
    ok <- is.finite(center) && is.finite(delta) &&
      is.finite(lower) && is.finite(upper) &&
      vivo >= lower - tol && vivo <= upper + tol &&
      vitro >= lower - tol && vitro <= upper + tol
    if (!isTRUE(ok)) bad <- c(bad, as.character(soft_meta$parameter[[i]]))
  }
  list(feasible = length(bad) == 0L, parameters = unique(bad))
}

joint_apply_soft_coupling_start_table <- function(init,
                                                  lower,
                                                  upper,
                                                  soft_meta,
                                                  cfg_raw,
                                                  invivo_parameter_table) {
  path <- joint_soft_coupling_start_table_path(cfg_raw, invivo_parameter_table)
  required <- !is.null(cfg_raw$joint_soft_coupling_parameters_table) ||
    !is.null(cfg_raw$joint_soft_coupling_parameters_table_path)
  tab <- joint_read_soft_coupling_start_table(path, required = required)
  if (!is.data.frame(tab) || !nrow(tab) || !is.data.frame(soft_meta) || !nrow(soft_meta)) {
    return(list(init = init, lower = lower, upper = upper, metadata = soft_meta, path = path, applied = data.frame()))
  }

  lookup <- joint_soft_start_lookup(soft_meta, names(init))
  init_out <- init
  lower_out <- lower
  upper_out <- upper
  meta_out <- soft_meta
  tab_is_delta <- vapply(seq_len(nrow(tab)), function(i) {
    hit <- lookup[lookup$input_name == as.character(tab$param_name[[i]]), , drop = FALSE]
    if (nrow(hit) != 1L) FALSE else isTRUE(hit$is_delta[[1]])
  }, logical(1))
  tab_order <- c(which(!tab_is_delta), which(tab_is_delta))
  applied <- lapply(tab_order, function(i) {
    pname <- as.character(tab$param_name[[i]])
    hit <- lookup[lookup$input_name == pname, , drop = FALSE]
    if (nrow(hit) != 1L) {
      stop(
        "joint_soft_coupling_parameters_table param_name is not an active soft-coupling parameter: ",
        pname,
        call. = FALSE
      )
    }
    opt_name <- as.character(hit$optimizer_name[[1]])
    if (!(opt_name %in% names(init_out))) {
      stop("Start-table optimizer parameter missing from joint init: ", opt_name, call. = FALSE)
    }
    center_t <- if (isTRUE(hit$is_delta[[1]])) {
      center_name <- as.character(hit$center_name[[1]])
      if (!(center_name %in% names(init_out))) {
        stop("Start-table delta center parameter missing from joint init: ", center_name, call. = FALSE)
      }
      as.numeric(init_out[[center_name]])
    } else {
      NA_real_
    }
    opt_value <- joint_soft_start_value_to_optimizer(tab$value[[i]], tab$scale[[i]], hit, center_t = center_t)
    if (!is.finite(opt_value)) {
      stop("Non-finite optimizer-scale value after converting start-table row: ", pname, call. = FALSE)
    }
    old_lower <- as.numeric(lower_out[[opt_name]])
    old_upper <- as.numeric(upper_out[[opt_name]])
    tol <- 1e-12
    if (opt_value < old_lower - tol) {
      stop(
        "joint_soft_coupling_parameters_table value for ", pname,
        " is below the optimizer lower bound (", signif(old_lower, 8), ").",
        call. = FALSE
      )
    }
    if (opt_value > old_upper + tol) {
      stop(
        "joint_soft_coupling_parameters_table value for ", pname,
        " is above the optimizer upper bound (", signif(old_upper, 8), ").",
        call. = FALSE
      )
    }
    init_out[[opt_name]] <<- opt_value
    data.frame(
      param_name = pname,
      optimizer_name = opt_name,
      parameter = as.character(hit$parameter[[1]]),
      is_delta = isTRUE(hit$is_delta[[1]]),
      scale = as.character(tab$scale[[i]]),
      input_value = as.numeric(tab$value[[i]]),
      optimizer_value = opt_value,
      natural_value = if (isTRUE(hit$is_delta[[1]])) NA_real_ else joint_transformed_to_natural(opt_value, hit$parameter[[1]]),
      optimizer_lower_before = old_lower,
      optimizer_upper_before = old_upper,
      optimizer_lower_after = as.numeric(lower_out[[opt_name]]),
      optimizer_upper_after = as.numeric(upper_out[[opt_name]]),
      bound_action = "inside",
      stringsAsFactors = FALSE
    )
  })
  check <- joint_check_soft_reconstruction_feasibility(init_out, meta_out)
  if (!isTRUE(check$feasible)) {
    stop(
      "joint_soft_coupling_parameters_table values reconstruct outside joint union bounds for: ",
      paste(check$parameters, collapse = ", "),
      call. = FALSE
    )
  }
  list(
    init = init_out,
    lower = lower_out,
    upper = upper_out,
    metadata = meta_out,
    path = path,
    applied = dplyr::bind_rows(applied)
  )
}

joint_warmup_default_sigmaN <- function() {
  0.05 / stats::qnorm(0.95)
}

joint_read_table_auto <- function(path) {
  tab <- tryCatch(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(tab) || ncol(tab) <= 1L) {
    tab <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  tab
}

joint_table_map <- function(tab, name_cols, value_cols) {
  name_col <- intersect(name_cols, names(tab))
  value_col <- intersect(value_cols, names(tab))
  if (!length(name_col) || !length(value_col)) return(NULL)
  nm <- trimws(as.character(tab[[name_col[[1]]]]))
  val <- suppressWarnings(as.numeric(tab[[value_col[[1]]]]))
  ok <- nzchar(nm) & is.finite(val)
  out <- val[ok]
  names(out) <- nm[ok]
  out[!duplicated(names(out))]
}

joint_invitro_parameter_specs <- function() {
  data.frame(
    param_symbol = c(
      "lam_max", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
      "buffer_n_exp", "p_wgd", "alpha_o2", "gamma_growth", "mu_hp",
      "gamma_mu", "O2_crit", "n_O", "p_mis_base", "sigma_growth",
      "sigma_kary", "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"
    ),
    param_name = c(
      "log10_lam_max", "log10_p_misseg", "log10_k_o_mis", "buffer_smax",
      "log10_buffer_beta", "log10_buffer_n_exp", "log10_p_wgd",
      "log10_alpha_o2", "gamma_growth", "log10_mu_hp", "gamma_mu",
      "log10_O2_crit", "n_O", "log10_p_mis_base", "log10_sigma_growth",
      "log10_sigma_kary", "init_mean_2N", "log10_init_sd_2N",
      "init_mean_4N", "log10_init_sd_4N"
    ),
    transform = c(
      "log10", "log10", "log10", "identity", "log10", "log10", "log10",
      "log10", "identity", "log10", "identity", "log10", "identity",
      "log10", "log10", "log10", "identity", "log10", "identity", "log10"
    ),
    stringsAsFactors = FALSE
  )
}

joint_transform_natural_best_map <- function(natural_map, kind) {
  if (is.null(natural_map) || !length(natural_map)) return(numeric(0))
  specs <- if (identical(kind, "invivo")) {
    INVIVO_ENV$parameter_table_specs()
  } else {
    joint_invitro_parameter_specs()
  }
  out <- numeric(0)
  for (symbol in names(natural_map)) {
    row <- specs[as.character(specs$param_symbol) == as.character(symbol), , drop = FALSE]
    if (nrow(row) != 1L) next
    val <- tryCatch(
      INVIVO_ENV$transform_param_slot(
        as.numeric(natural_map[[symbol]]),
        as.character(row$transform[[1]]),
        as.character(symbol),
        "warm-up best_params.tsv value"
      ),
      error = function(e) NA_real_
    )
    if (is.finite(val)) out[[as.character(row$param_name[[1]])]] <- val
  }
  out
}

joint_read_best_seed_transformed_map <- function(seed_path, kind) {
  seed_path <- normalizePath(seed_path, mustWork = FALSE)
  if (!file.exists(seed_path)) {
    stop("Warm-up ", kind, " seed path not found: ", seed_path, call. = FALSE)
  }
  candidates <- if (dir.exists(seed_path)) {
    c(
      file.path(seed_path, "best_params_transformed.tsv"),
      file.path(seed_path, "fit_parameter_stages.tsv"),
      file.path(seed_path, "checkpoint", "best_params_transformed_latest.tsv"),
      file.path(seed_path, "best_params.tsv")
    )
  } else {
    seed_path
  }
  for (path in candidates) {
    if (!file.exists(path)) next
    tab <- joint_read_table_auto(path)
    transformed <- joint_table_map(
      tab,
      name_cols = c("transformed_parameter", "param_name", "parameter"),
      value_cols = c("transformed_value", "value", "init_value")
    )
    if (!is.null(transformed) && length(transformed) &&
        (!identical(basename(path), "best_params.tsv") ||
         any(grepl("^log10_|^ivt__|^delta__", names(transformed))))) {
      attr(transformed, "source_path") <- normalizePath(path, mustWork = FALSE)
      attr(transformed, "source_scale") <- "transformed"
      return(transformed)
    }
    if (identical(basename(path), "best_params.tsv")) {
      natural <- joint_table_map(tab, name_cols = c("parameter", "param_symbol"), value_cols = c("value"))
      transformed <- joint_transform_natural_best_map(natural, kind = kind)
      if (length(transformed)) {
        attr(transformed, "source_path") <- normalizePath(path, mustWork = FALSE)
        attr(transformed, "source_scale") <- "natural"
        return(transformed)
      }
    }
  }
  stop("Could not read warm-up transformed best parameters from: ", seed_path, call. = FALSE)
}

joint_warmup_seed_path <- function(cfg_raw, keys) {
  raw <- NULL
  for (key in keys) {
    if (!is.null(cfg_raw[[key]])) {
      raw <- cfg_raw[[key]]
      break
    }
  }
  joint_resolve_path_local(raw, cfg_raw, must_exist = TRUE)
}

joint_set_warmup_init_value <- function(init,
                                        lower,
                                        upper,
                                        name,
                                        value,
                                        source,
                                        source_detail = "") {
  name <- as.character(name)
  if (!(name %in% names(init)) || !is.finite(value)) {
    return(list(init = init, lower = lower, upper = upper, row = NULL))
  }
  old <- as.numeric(init[[name]])
  old_lower <- as.numeric(lower[[name]])
  old_upper <- as.numeric(upper[[name]])
  if (!is.finite(old_lower) || !is.finite(old_upper) || old_lower > old_upper) {
    return(list(init = init, lower = lower, upper = upper, row = NULL))
  }
  if (old_lower == old_upper) {
    return(list(init = init, lower = lower, upper = upper, row = NULL))
  }
  value_use <- as.numeric(value)
  bound_action <- "inside"
  tol <- 1e-12
  if (value_use < old_lower - tol || value_use > old_upper + tol) {
    stop(
      "Warm-up value for ", name, " from ", source,
      " is outside optimizer bounds [", signif(old_lower, 8), ", ",
      signif(old_upper, 8), "].",
      call. = FALSE
    )
  }
  init[[name]] <- value_use
  row <- data.frame(
    optimizer_name = name,
    source = source,
    source_detail = as.character(source_detail),
    old_init = old,
    warmup_init = value_use,
    optimizer_lower_before = old_lower,
    optimizer_upper_before = old_upper,
    optimizer_lower_after = as.numeric(lower[[name]]),
    optimizer_upper_after = as.numeric(upper[[name]]),
    bound_action = bound_action,
    stringsAsFactors = FALSE
  )
  list(init = init, lower = lower, upper = upper, row = row)
}

joint_named_num <- function(x, name, default = NA_real_) {
  if (is.null(x) || is.null(names(x)) || !(name %in% names(x))) return(default)
  val <- suppressWarnings(as.numeric(x[[name]]))
  if (!length(val) || !is.finite(val[[1]])) default else val[[1]]
}

joint_apply_warmup_initial_values <- function(init,
                                              lower,
                                              upper,
                                              soft_meta,
                                              cfg_raw,
                                              invivo,
                                              invitro,
                                              ivt_extra_names) {
  enabled <- isTRUE(as_bool(cfg_raw$joint_warmup_enable, FALSE))
  sigmaN <- as_num(cfg_raw$joint_warmup_sigmaN, joint_warmup_default_sigmaN())
  if (!is.finite(sigmaN) || sigmaN <= 0) {
    stop("joint_warmup_sigmaN must be > 0.", call. = FALSE)
  }
  if (!enabled) {
    return(list(
      init = init, lower = lower, upper = upper,
      enabled = FALSE, sigmaN = sigmaN,
      invivo_seed_dir = NA_character_, invitro_seed_dir = NA_character_,
      invivo_source_path = NA_character_, invitro_source_path = NA_character_,
      applied = data.frame()
    ))
  }

  invivo_seed <- joint_warmup_seed_path(cfg_raw, c("joint_warmup_invivo_seed_dir", "joint_warmup_invivo_best_seed_dir"))
  invitro_seed <- joint_warmup_seed_path(cfg_raw, c("joint_warmup_invitro_seed_dir", "joint_warmup_invitro_best_seed_dir", "joint_warmup_vitro_seed_dir"))
  invivo_map <- joint_read_best_seed_transformed_map(invivo_seed, kind = "invivo")
  invitro_map <- joint_read_best_seed_transformed_map(invitro_seed, kind = "invitro")

  init_out <- init
  lower_out <- lower
  upper_out <- upper
  records <- list()
  apply_value <- function(name, value, source, detail = "") {
    res <- joint_set_warmup_init_value(
      init = init_out,
      lower = lower_out,
      upper = upper_out,
      name = name,
      value = value,
      source = source,
      source_detail = detail
    )
    init_out <<- res$init
    lower_out <<- res$lower
    upper_out <<- res$upper
    if (!is.null(res$row)) records[[length(records) + 1L]] <<- res$row
  }

  if (is.data.frame(soft_meta) && nrow(soft_meta)) {
    for (i in seq_len(nrow(soft_meta))) {
      center_name <- as.character(soft_meta$center_name[[i]])
      delta_name <- as.character(soft_meta$delta_name[[i]])
      invitro_name <- as.character(soft_meta$invitro_name[[i]])
      vivo <- joint_named_num(invivo_map, center_name)
      vitro <- joint_named_num(invitro_map, invitro_name)
      if (!is.finite(vivo) || !is.finite(vitro)) next
      center <- (vivo + vitro) / 2
      delta <- vivo - vitro
      apply_value(center_name, center, "soft_center_from_best_seed_mean", paste(center_name, invitro_name, sep = "|"))
      apply_value(delta_name, delta, "soft_delta_from_best_seed_difference", paste(center_name, invitro_name, sep = "|"))
    }
  }

  soft_params <- if (is.data.frame(soft_meta) && nrow(soft_meta)) as.character(soft_meta$parameter) else character(0)
  shared_symbols <- setdiff(joint_shared_natural_param_names(), soft_params)
  for (symbol in shared_symbols) {
    pname <- invitro_shared_param_name_for_natural(symbol)
    if (is.null(pname) || !(pname %in% names(init_out))) next
    vivo <- joint_named_num(invivo_map, pname)
    vitro <- joint_named_num(invitro_map, pname)
    if (!is.finite(vivo) || !is.finite(vitro)) next
    apply_value(pname, (vivo + vitro) / 2, "shared_mean_from_best_seeds", paste(symbol, pname, sep = "|"))
  }

  shared_ivt <- shared_invitro_param_names()
  invivo_only <- setdiff(names(invivo$param_bundle$optimizer$init), shared_ivt)
  for (pname in invivo_only) {
    val <- joint_named_num(invivo_map, pname)
    if (is.finite(val)) apply_value(pname, val, "invivo_best_seed", pname)
  }

  for (ivt_name in as.character(ivt_extra_names)) {
    opt_name <- paste0("ivt__", ivt_name)
    val <- joint_named_num(invitro_map, ivt_name)
    if (is.finite(val)) apply_value(opt_name, val, "invitro_best_seed", ivt_name)
  }

  check <- joint_check_soft_reconstruction_feasibility(init_out, soft_meta)
  if (!isTRUE(check$feasible)) {
    stop(
      "Warm-up values reconstruct outside joint union bounds for: ",
      paste(check$parameters, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    init = init_out,
    lower = lower_out,
    upper = upper_out,
    enabled = TRUE,
    sigmaN = sigmaN,
    invivo_seed_dir = invivo_seed,
    invitro_seed_dir = invitro_seed,
    invivo_source_path = as.character(attr(invivo_map, "source_path")),
    invitro_source_path = as.character(attr(invitro_map, "source_path")),
    applied = dplyr::bind_rows(records)
  )
}

merge_joint_shared_optimizer_bounds <- function(invivo,
                                                invitro) {
  shared_names <- joint_shared_natural_param_names()
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
    invivo_spec_row <- INVIVO_ENV$parameter_table_specs()
    invivo_spec_row <- invivo_spec_row[as.character(invivo_spec_row$param_symbol) == as.character(symbol), , drop = FALSE]
    if (nrow(invivo_spec_row) != 1L) {
      stop("Missing in vivo transformed spec for shared parameter '", symbol, "'.", call. = FALSE)
    }
    vivo_bound_lower_t <- INVIVO_ENV$transform_param_slot(
      as.numeric(inv_row$lower_bound[[1]]),
      as.character(invivo_spec_row$transform[[1]]),
      symbol,
      "invivo_lower"
    )
    vivo_bound_upper_t <- INVIVO_ENV$transform_param_slot(
      as.numeric(inv_row$upper_bound[[1]]),
      as.character(invivo_spec_row$transform[[1]]),
      symbol,
      "invivo_upper"
    )
    vitro_bound_lower_t <- transform_invitro_shared_slot(ivt_row$lower_bound[[1]], symbol, "invitro_lower")
    vitro_bound_upper_t <- transform_invitro_shared_slot(ivt_row$upper_bound[[1]], symbol, "invitro_upper")
    joint_union_lower_t <- min(vivo_bound_lower_t, vitro_bound_lower_t)
    joint_union_upper_t <- max(vivo_bound_upper_t, vitro_bound_upper_t)
    if (!is.finite(joint_union_lower_t) ||
        !is.finite(joint_union_upper_t) ||
        joint_union_lower_t > joint_union_upper_t) {
      stop("Invalid joint union transformed bounds for shared parameter '", symbol, "'.", call. = FALSE)
    }
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
      joint_union_lower_t = joint_union_lower_t,
      joint_union_upper_t = joint_union_upper_t,
      stringsAsFactors = FALSE
    )
  }

  invivo_init <- invivo$param_bundle$optimizer$init
  invivo_lower <- invivo$param_bundle$optimizer$lower
  invivo_upper <- invivo$param_bundle$optimizer$upper
  invivo_specs <- INVIVO_ENV$parameter_table_specs()
  bounds_summary <- dplyr::bind_rows(summary_rows)

  for (symbol in shared_names) {
    spec_row <- invivo_specs[as.character(invivo_specs$param_symbol) == as.character(symbol), , drop = FALSE]
    if (nrow(spec_row) == 1L) {
      pname <- as.character(spec_row$param_name[[1]])
      if (pname %in% names(invivo_lower)) {
        row <- bounds_summary[as.character(bounds_summary$parameter) == as.character(symbol), , drop = FALSE]
        invivo_lower[[pname]] <- as.numeric(row$joint_union_lower_t[[1]])
        invivo_upper[[pname]] <- as.numeric(row$joint_union_upper_t[[1]])
      }
    }
  }

  list(
    invivo_optimizer = list(init = invivo_init, lower = invivo_lower, upper = invivo_upper),
    natural = merged_nat,
    summary = bounds_summary
  )
}

split_joint_natural_parameter_tables <- function(invivo_param_df,
                                                 invitro_param_df,
                                                 soft_split_params = character()) {
  all_shared_names <- joint_shared_natural_param_names()
  soft_split_params <- intersect(as.character(soft_split_params), all_shared_names)
  shared_names <- setdiff(all_shared_names, soft_split_params)
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
  invivo_only <- invivo_param_df[!(invivo_param_df$parameter %in% all_shared_names), , drop = FALSE]
  invitro_only <- invitro_param_df[!(invitro_param_df$parameter %in% all_shared_names), , drop = FALSE]
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
                                                 invivo_cfg) {
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
  par_t
}

build_joint_context <- function(argv) {
  cfg_raw <- resolve_joint_raw_config(argv)
  restriction_flags <- resolve_joint_restriction_flags(cfg_raw)
  invivo <- build_joint_invivo_context(cfg_raw)
  invitro <- build_joint_invitro_context(cfg_raw)

  invivo_names <- names(invivo$param_bundle$optimizer$init)
  shared_ivt <- shared_invitro_param_names()
  ivt_extra_names <- setdiff(invitro$spec$param_name, shared_ivt)
  ivt_extra_prefixed <- paste0("ivt__", ivt_extra_names)
  joint_bounds <- merge_joint_shared_optimizer_bounds(
    invivo = invivo,
    invitro = invitro
  )
  soft_split_params <- joint_soft_split_natural_param_names(
    cfg_raw = cfg_raw
  )
  soft_meta <- joint_soft_coupling_metadata(
    split_params = soft_split_params,
    joint_bounds = joint_bounds,
    invivo = invivo,
    invitro = invitro,
    cfg_raw = cfg_raw
  )
  if (nrow(soft_meta) > 0L) {
    for (i in seq_len(nrow(soft_meta))) {
      center_name <- soft_meta$center_name[[i]]
      joint_bounds$invivo_optimizer$init[[center_name]] <- soft_meta$center_init_t[[i]]
      joint_bounds$invivo_optimizer$lower[[center_name]] <- soft_meta$joint_union_lower_t[[i]]
      joint_bounds$invivo_optimizer$upper[[center_name]] <- soft_meta$joint_union_upper_t[[i]]
    }
  }

  ivt_init <- setNames(as.numeric(invitro$spec$init[match(ivt_extra_names, invitro$spec$param_name)]), ivt_extra_prefixed)
  ivt_lower <- setNames(as.numeric(invitro$spec$lower[match(ivt_extra_names, invitro$spec$param_name)]), ivt_extra_prefixed)
  ivt_upper <- setNames(as.numeric(invitro$spec$upper[match(ivt_extra_names, invitro$spec$param_name)]), ivt_extra_prefixed)
  delta_init <- setNames(rep(0.0, nrow(soft_meta)), soft_meta$delta_name)
  delta_lower <- setNames(as.numeric(soft_meta$delta_lower_t), soft_meta$delta_name)
  delta_upper <- setNames(as.numeric(soft_meta$delta_upper_t), soft_meta$delta_name)

  init <- c(joint_bounds$invivo_optimizer$init, ivt_init, delta_init)
  lower <- c(joint_bounds$invivo_optimizer$lower, ivt_lower, delta_lower)
  upper <- c(joint_bounds$invivo_optimizer$upper, ivt_upper, delta_upper)
  warmup <- joint_apply_warmup_initial_values(
    init = init,
    lower = lower,
    upper = upper,
    soft_meta = soft_meta,
    cfg_raw = cfg_raw,
    invivo = invivo,
    invitro = invitro,
    ivt_extra_names = ivt_extra_names
  )
  init <- warmup$init
  lower <- warmup$lower
  upper <- warmup$upper
  soft_start <- joint_apply_soft_coupling_start_table(
    init = init,
    lower = lower,
    upper = upper,
    soft_meta = soft_meta,
    cfg_raw = cfg_raw,
    invivo_parameter_table = invivo$cfg$parameter_table
  )
  init <- soft_start$init
  lower <- soft_start$lower
  upper <- soft_start$upper
  soft_meta <- soft_start$metadata

  list(
    raw = cfg_raw,
    invivo = invivo,
    invitro = invitro,
    joint_shared_bounds = joint_bounds$summary,
    invivo_names = invivo_names,
    ivt_extra_names = ivt_extra_names,
    ivt_extra_prefixed = ivt_extra_prefixed,
    joint_soft_coupling = list(
      enabled = nrow(soft_meta) > 0L,
      params = soft_split_params,
      metadata = soft_meta,
      sigma_default = as_num(cfg_raw$joint_soft_coupling_sigma_default, 1.5)
    ),
    joint_soft_coupling_start_table = list(
      path = soft_start$path,
      applied = soft_start$applied
    ),
    joint_warmup = warmup,
    init = init,
    lower = lower,
    upper = upper,
    joint_weight_invivo = as_num(cfg_raw$joint_weight_invivo, 1.0),
    joint_weight_invitro = as_num(cfg_raw$joint_weight_invitro, 1.0),
    seed = as_int(cfg_raw$seed, 1L),
    itermax = as_int(cfg_raw$itermax, 40L),
    de_reltol = as_num(cfg_raw$de_reltol, 1e-4),
    de_steptol = as_int(cfg_raw$de_steptol, 25L),
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

joint_deoptim_initial_population <- function(ctx, NP_use) {
  init <- as.numeric(ctx$init)
  names(init) <- names(ctx$init)
  lower <- as.numeric(ctx$lower)
  upper <- as.numeric(ctx$upper)
  if (any(!is.finite(init)) || any(!is.finite(lower)) || any(!is.finite(upper))) {
    stop("Joint optimizer init/lower/upper must be finite before DEoptim starts.", call. = FALSE)
  }
  if (any(lower > upper)) {
    bad <- names(ctx$init)[lower > upper]
    stop("Joint optimizer lower > upper for: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (any(init < lower - 1e-12 | init > upper + 1e-12)) {
    bad <- names(ctx$init)[init < lower - 1e-12 | init > upper + 1e-12]
    stop("Joint optimizer init outside bounds for: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  n_par <- length(init)
  NP_use <- as.integer(NP_use)
  pop <- matrix(NA_real_, nrow = NP_use, ncol = n_par)
  if (isTRUE(ctx$joint_warmup$enabled)) {
    sigmaN <- as_num(ctx$joint_warmup$sigmaN, joint_warmup_default_sigmaN())
    if (!is.finite(sigmaN) || sigmaN <= 0) stop("joint_warmup_sigmaN must be > 0.", call. = FALSE)
    for (j in seq_len(n_par)) {
      if (lower[[j]] == upper[[j]]) {
        pop[, j] <- lower[[j]]
        next
      }
      scale_ref <- abs(init[[j]])
      if (!is.finite(scale_ref) || scale_ref <= 1e-12) {
        scale_ref <- upper[[j]] - lower[[j]]
      }
      sd_use <- as.numeric(scale_ref) * sigmaN
      if (!is.finite(sd_use) || sd_use <= 0) {
        pop[, j] <- init[[j]]
        next
      }
      for (i in seq_len(NP_use)) {
        repeat {
          draw <- stats::rnorm(1L, mean = init[[j]], sd = sd_use)
          if (is.finite(draw) && draw >= lower[[j]] && draw <= upper[[j]]) {
            pop[i, j] <- draw
            break
          }
        }
      }
    }
  } else {
    pop <- matrix(stats::runif(NP_use * n_par), nrow = NP_use, ncol = n_par)
    pop <- sweep(pop, 2L, upper - lower, "*")
    pop <- sweep(pop, 2L, lower, "+")
  }
  pop[1L, ] <- init
  colnames(pop) <- names(ctx$init)
  pop
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

joint_soft_coupling_active <- function(ctx) {
  !is.null(ctx$joint_soft_coupling) &&
    isTRUE(ctx$joint_soft_coupling$enabled) &&
    is.data.frame(ctx$joint_soft_coupling$metadata) &&
    nrow(ctx$joint_soft_coupling$metadata) > 0L
}

joint_soft_coupling_penalty_components <- function(par_t, ctx) {
  if (!joint_soft_coupling_active(ctx)) {
    return(list(total = 0.0, terms = data.frame()))
  }
  par_t <- as.numeric(par_t)
  names(par_t) <- names(ctx$init)
  meta <- ctx$joint_soft_coupling$metadata
  terms <- lapply(seq_len(nrow(meta)), function(i) {
    delta_name <- meta$delta_name[[i]]
    delta <- as.numeric(par_t[[delta_name]])
    sigma <- as.numeric(meta$sigma_delta[[i]])
    penalty <- delta^2 / (2 * sigma^2)
    data.frame(
      parameter = meta$parameter[[i]],
      delta_name = delta_name,
      delta_transformed = delta,
      regularization_sigma = sigma,
      penalty_paid = penalty,
      stringsAsFactors = FALSE
    )
  })
  terms <- dplyr::bind_rows(terms)
  total <- sum(as.numeric(terms$penalty_paid), na.rm = TRUE)
  list(total = total, terms = terms)
}

joint_build_context_specific_transformed_vectors <- function(par_t, ctx) {
  par_t <- as.numeric(par_t)
  names(par_t) <- names(ctx$init)
  invivo_par <- par_t[ctx$invivo_names]
  if (!joint_soft_coupling_active(ctx)) {
    return(list(
      invivo_par = invivo_par,
      soft_derived = data.frame(),
      feasible = TRUE,
      infeasible_parameters = character(0)
    ))
  }

  meta <- ctx$joint_soft_coupling$metadata
  rows <- lapply(seq_len(nrow(meta)), function(i) {
    center_name <- meta$center_name[[i]]
    delta_name <- meta$delta_name[[i]]
    center <- as.numeric(par_t[[center_name]])
    delta <- as.numeric(par_t[[delta_name]])
    vivo_raw <- center + delta / 2
    vitro_raw <- center - delta / 2
    joint_lower <- as.numeric(meta$joint_union_lower_t[[i]])
    joint_upper <- as.numeric(meta$joint_union_upper_t[[i]])
    feasible <- is.finite(vivo_raw) && is.finite(vitro_raw) &&
      is.finite(joint_lower) && is.finite(joint_upper) &&
      vivo_raw >= joint_lower - 1e-10 && vivo_raw <= joint_upper + 1e-10 &&
      vitro_raw >= joint_lower - 1e-10 && vitro_raw <= joint_upper + 1e-10
    invivo_par[[center_name]] <<- vivo_raw
    data.frame(
      parameter = meta$parameter[[i]],
      center_name = center_name,
      invitro_name = meta$invitro_name[[i]],
      delta_name = delta_name,
      transform = meta$transform[[i]],
      center_transformed = center,
      delta_transformed = delta,
      vivo_transformed = vivo_raw,
      vitro_transformed = vitro_raw,
      joint_union_lower_transformed = joint_lower,
      joint_union_upper_transformed = joint_upper,
      feasible_at_point = isTRUE(feasible),
      stringsAsFactors = FALSE
    )
  })
  soft_derived <- dplyr::bind_rows(rows)
  feasible <- !nrow(soft_derived) || all(soft_derived$feasible_at_point)
  list(
    invivo_par = invivo_par,
    soft_derived = soft_derived,
    feasible = feasible,
    infeasible_parameters = if (feasible) character(0) else as.character(soft_derived$parameter[!soft_derived$feasible_at_point])
  )
}

joint_apply_soft_coupling_to_invitro <- function(ivt_par, soft_derived) {
  if (!is.data.frame(soft_derived) || nrow(soft_derived) == 0L) {
    return(ivt_par)
  }
  for (i in seq_len(nrow(soft_derived))) {
    nm <- soft_derived$invitro_name[[i]]
    if (nm %in% names(ivt_par)) {
      ivt_par[[nm]] <- as.numeric(soft_derived$vitro_transformed[[i]])
    }
  }
  ivt_par
}

joint_soft_coupling_summary_table <- function(par_t, ctx) {
  if (!joint_soft_coupling_active(ctx)) return(data.frame())
  derived <- joint_build_context_specific_transformed_vectors(par_t, ctx)$soft_derived
  penalty <- joint_soft_coupling_penalty_components(par_t, ctx)$terms
  if (!nrow(derived)) return(data.frame())
  out <- dplyr::left_join(
    derived,
    penalty[, c("parameter", "regularization_sigma", "penalty_paid"), drop = FALSE],
    by = "parameter"
  )
  out$split_enabled <- TRUE
  out$center_natural <- mapply(joint_transformed_to_natural, out$center_transformed, out$parameter)
  out$vivo_natural <- mapply(joint_transformed_to_natural, out$vivo_transformed, out$parameter)
  out$vitro_natural <- mapply(joint_transformed_to_natural, out$vitro_transformed, out$parameter)
  out$joint_union_lower_bound <- mapply(joint_transformed_to_natural, out$joint_union_lower_transformed, out$parameter)
  out$joint_union_upper_bound <- mapply(joint_transformed_to_natural, out$joint_union_upper_transformed, out$parameter)
  out$feasible_at_solution <- out$feasible_at_point
  out$ratio_vivo_to_vitro <- ifelse(
    is.finite(out$vivo_natural) & is.finite(out$vitro_natural) & out$vitro_natural != 0,
    out$vivo_natural / out$vitro_natural,
    NA_real_
  )
  log_scale <- grepl("^log10", out$center_name)
  probability_scale <- joint_soft_probability_parameter(out$parameter)
  out$natural_difference_vivo_to_vitro <- out$vivo_natural - out$vitro_natural
  out$transformed_difference_vivo_to_vitro <- out$vivo_transformed - out$vitro_transformed
  out$log10_ratio_vivo_to_vitro <- ifelse(
    log_scale,
    out$transformed_difference_vivo_to_vitro,
    NA_real_
  )
  out$fold_change_vivo_to_vitro <- out$ratio_vivo_to_vitro
  out$logit_vivo <- ifelse(probability_scale, joint_safe_logit(out$vivo_natural), NA_real_)
  out$logit_vitro <- ifelse(probability_scale, joint_safe_logit(out$vitro_natural), NA_real_)
  out$logit_difference_vivo_to_vitro <- out$logit_vivo - out$logit_vitro
  out$odds_ratio_vivo_to_vitro <- ifelse(
    probability_scale & is.finite(out$logit_difference_vivo_to_vitro),
    exp(out$logit_difference_vivo_to_vitro),
    NA_real_
  )
  out$delta_interpretable <- ifelse(
    log_scale,
    out$transformed_difference_vivo_to_vitro,
    out$natural_difference_vivo_to_vitro
  )
  out[, c(
    "parameter",
    "split_enabled",
    "center_name",
    "delta_name",
    "center_transformed",
    "delta_transformed",
    "vivo_transformed",
    "vitro_transformed",
    "center_natural",
    "vivo_natural",
    "vitro_natural",
    "delta_interpretable",
    "natural_difference_vivo_to_vitro",
    "transformed_difference_vivo_to_vitro",
    "log10_ratio_vivo_to_vitro",
    "fold_change_vivo_to_vitro",
    "ratio_vivo_to_vitro",
    "logit_vivo",
    "logit_vitro",
    "logit_difference_vivo_to_vitro",
    "odds_ratio_vivo_to_vitro",
    "regularization_sigma",
    "penalty_paid",
    "joint_union_lower_transformed",
    "joint_union_upper_transformed",
    "joint_union_lower_bound",
    "joint_union_upper_bound",
    "feasible_at_solution"
  ), drop = FALSE]
}

joint_objective_components <- function(par_t, ctx) {
  par_t <- as.numeric(par_t)
  names(par_t) <- names(ctx$init)
  context_vectors <- joint_build_context_specific_transformed_vectors(par_t, ctx)
  if (!isTRUE(context_vectors$feasible)) {
    penalty <- as_num(ctx$joint_constraint_penalty, 1e9)
    if (!is.finite(penalty) || penalty <= 0) penalty <- 1e9
    soft_penalty <- joint_soft_coupling_penalty_components(par_t, ctx)
    return(list(
      objective = penalty,
      objective_unpenalized = penalty,
      objective_soft_coupling = as_num(soft_penalty$total, 0),
      constraint_metrics = list(
        joint_constraint_penalty_total = penalty,
        joint_constraints_pass = FALSE,
        joint_soft_coupling_feasible = FALSE,
        joint_soft_coupling_infeasible_parameters = paste(context_vectors$infeasible_parameters, collapse = ",")
      ),
      soft_coupling = soft_penalty,
      soft_coupling_derived = context_vectors$soft_derived,
      invivo = list(L = penalty, L_data = NA_real_, L_prior = NA_real_, L_b = NA_real_, L_p = NA_real_),
      invitro = INVITRO_ENV$make_penalty_components(
        objective = penalty,
        reason = paste0("joint_soft_coupling_infeasible: ", paste(context_vectors$infeasible_parameters, collapse = ","))
      ),
      invivo_run_params = list(),
      invitro_run_params = list(),
      invitro_transformed = numeric(0)
    ))
  }
  invivo_par <- context_vectors$invivo_par
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
    invivo_cfg = ctx$invivo$cfg
  )
  ivt_par <- joint_apply_soft_coupling_to_invitro(
    ivt_par = ivt_par,
    soft_derived = context_vectors$soft_derived
  )
  invitro_run_params <- INVITRO_ENV$ivt_optim_par_to_run_params(ivt_par, cfg = ctx$invitro$cfg)
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
  soft_penalty <- joint_soft_coupling_penalty_components(par_t, ctx)
  joint <- joint + as_num(soft_penalty$total, 0)
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
    objective_soft_coupling = as_num(soft_penalty$total, 0),
    constraint_metrics = constraint_metrics,
    soft_coupling = soft_penalty,
    soft_coupling_derived = context_vectors$soft_derived,
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
  invivo_params <- filter_family_specific_run_params_for_output_common(invivo_params)
  invitro_params <- best_comp$invitro_run_params[vapply(best_comp$invitro_run_params, is.numeric, logical(1))]
  invitro_params <- filter_family_specific_run_params_for_output_common(invitro_params)
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
    soft_split_params = ctx$joint_soft_coupling$params
  )
  soft_coupling_df <- joint_soft_coupling_summary_table(best_par_t, ctx)
  write.table(invivo_param_df, file = file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(invitro_param_df, file = file.path(out_dir, "invitro_effective_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_long_df, file = file.path(out_dir, "joint_best_params_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_tables$shared, file = file.path(out_dir, "joint_params_shared.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_tables$invivo_only, file = file.path(out_dir, "joint_params_invivo_only.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(param_tables$invitro_only, file = file.path(out_dir, "joint_params_invitro_only.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (joint_soft_coupling_active(ctx)) {
    write.table(soft_coupling_df, file = file.path(out_dir, "joint_soft_coupling.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  write.table(ctx$joint_shared_bounds, file = file.path(out_dir, "joint_shared_bounds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  preds <- INVIVO_ENV$collect_predictions(best_comp$invivo_run_params, ctx$invivo$scenarios, ctx$invivo$cfg)
  write.table(preds$burden, file = file.path(out_dir, "invivo_burden_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$ploidy, file = file.path(out_dir, "invivo_terminal_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$necrosis, file = file.path(out_dir, "invivo_necrosis_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$burden, file = file.path(out_dir, "burden_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$ploidy, file = file.path(out_dir, "terminal_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$necrosis, file = file.path(out_dir, "necrosis_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

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
      "objective_soft_coupling",
      "objective_constraints",
      "joint_soft_coupling_enabled",
      "joint_soft_coupling_n_params",
      "weight_invivo",
      "weight_invitro",
      "objective_invivo",
      "objective_invivo_data",
      "objective_invivo_prior",
      "objective_invivo_burden",
      "objective_invivo_ploidy",
      "objective_invivo_necrosis",
      "objective_invivo_necrosis_raw",
      "objective_invitro",
      "invitro_growth_loglik",
      "invitro_ploidy_loglik",
      "invitro_flow_loglik"
    ),
    value = as.character(c(
      best_comp$objective,
      best_comp$objective_unpenalized,
      as_num(best_comp$objective_soft_coupling, 0),
      as_num(best_comp$constraint_metrics$joint_constraint_penalty_total, 0),
      joint_soft_coupling_active(ctx),
      length(ctx$joint_soft_coupling$params),
      ctx$joint_weight_invivo,
      ctx$joint_weight_invitro,
      best_comp$invivo$L,
      best_comp$invivo$L_data,
      best_comp$invivo$L_prior,
      best_comp$invivo$L_b,
      best_comp$invivo$L_p,
      best_comp$invivo$L_n,
      best_comp$invivo$L_n_raw,
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
  if (joint_soft_coupling_active(ctx) && nrow(best_comp$soft_coupling$terms) > 0L) {
    soft_component_df <- data.frame(
      component = paste0("joint_soft_coupling_penalty__", best_comp$soft_coupling$terms$parameter),
      value = as.character(best_comp$soft_coupling$terms$penalty_paid),
      stringsAsFactors = FALSE
    )
    joint_components <- dplyr::bind_rows(joint_components, soft_component_df)
  }
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
  optimizer_interrupted <- isTRUE(.first_non_null_local(optimizer_trace$interrupted, FALSE))
  optimizer_iter_target <- as.integer(.first_non_null_local(optimizer_trace$iter_target, ctx$itermax, NA_integer_))
  if (!is.finite(optimizer_iter_target) || is.na(optimizer_iter_target)) optimizer_iter_target <- NA_integer_
  optimizer_iter_completed <- as.integer(.first_non_null_local(optimizer_trace$iter_completed, NA_integer_))
  if (!is.finite(optimizer_iter_completed) || is.na(optimizer_iter_completed)) {
    optimizer_iter_completed <- joint_deoptim_iter_completed(de_fit, optimizer_iter_target)
  }
  deoptim_stop_iteration <- optimizer_iter_completed
  deoptim_iter_target <- optimizer_iter_target
  deoptim_stop_reason <- as.character(.first_non_null_local(optimizer_trace$deoptim_stop_reason, NA_character_))
  if (!length(deoptim_stop_reason) ||
      is.na(deoptim_stop_reason) ||
      !nzchar(trimws(deoptim_stop_reason)) ||
      tolower(trimws(deoptim_stop_reason)) %in% c("na", "nan", "null")) {
    deoptim_stop_reason <- joint_deoptim_stop_reason(
      iter_completed = optimizer_iter_completed,
      iter_target = optimizer_iter_target,
      interrupted = optimizer_interrupted
    )
  }
  optimizer_de_reltol <- as.numeric(.first_non_null_local(optimizer_trace$de_reltol, ctx$de_reltol, NA_real_))
  optimizer_de_steptol <- as.integer(.first_non_null_local(optimizer_trace$de_steptol, ctx$de_steptol, NA_integer_))

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
      "optimizer_interrupted",
      "optimizer_iter_completed",
      "optimizer_iter_target",
      "deoptim_stop_iteration",
      "deoptim_iter_target",
      "deoptim_stop_reason",
      "objective",
      "objective_invivo",
      "objective_invivo_data",
      "objective_invivo_prior",
      "objective_invivo_burden",
      "objective_invivo_ploidy",
      "objective_invivo_necrosis",
      "objective_invivo_necrosis_raw",
      "objective_invivo_necrosis_neg2loglik_raw",
      "use_necrosis_loss",
      "lambda_necrosis",
      "sigma_necrosis_logit",
      "necrosis_fraction_eps",
      "necrosis_mapping_csv",
      "n_necrosis",
      "n_necrosis_obs_total",
      "objective_invitro",
      "objective_soft_coupling",
      "objective_constraints",
      "joint_weight_invivo",
      "joint_weight_invitro",
      "joint_invitro_growth_weight",
      "joint_invitro_ploidy_weight",
      "joint_invitro_flow_weight",
      "joint_soft_coupling_enabled",
	      "joint_soft_coupling_params",
	      "joint_soft_coupling_n_params",
	      "joint_soft_coupling_sigma_default",
	      "joint_warmup_enabled",
	      "joint_warmup_sigmaN",
	      "joint_warmup_invivo_seed_dir",
	      "joint_warmup_invitro_seed_dir",
	      "joint_warmup_invivo_source_path",
	      "joint_warmup_invitro_source_path",
      "joint_restriction",
      "seed",
      "itermax",
      "de_reltol",
      "de_steptol",
      "NP_requested",
      "NP_used",
      "joint_np_min_factor",
      "n_cores_requested",
      "n_cores_used",
      "n_parameters",
      "n_invivo_scenarios"
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
      as.character(optimizer_interrupted),
      as.character(optimizer_iter_completed),
      as.character(optimizer_iter_target),
      as.character(deoptim_stop_iteration),
      as.character(deoptim_iter_target),
      as.character(deoptim_stop_reason),
      as.character(best_comp$objective),
      as.character(best_comp$invivo$L),
      as.character(best_comp$invivo$L_data),
      as.character(best_comp$invivo$L_prior),
      as.character(best_comp$invivo$L_b),
      as.character(best_comp$invivo$L_p),
      as.character(best_comp$invivo$L_n),
      as.character(best_comp$invivo$L_n_raw),
      as.character(best_comp$invivo$objective_necrosis_neg2loglik_raw),
      as.character(ctx$invivo$cfg$use_necrosis_loss),
      as.character(ctx$invivo$cfg$lambda_necrosis),
      as.character(ctx$invivo$cfg$sigma_necrosis_logit),
      as.character(ctx$invivo$cfg$necrosis_fraction_eps),
      as.character(if (is.null(ctx$invivo$cfg$necrosis_mapping_csv)) NA_character_ else normalizePath(ctx$invivo$cfg$necrosis_mapping_csv, mustWork = FALSE)),
      as.character(.first_non_null_local(best_comp$invivo$n_necrosis, NA_integer_)),
      as.character(.first_non_null_local(best_comp$invivo$n_necrosis_obs_total, NA_integer_)),
      as.character(best_comp$invitro$objective),
      as.character(as_num(best_comp$objective_soft_coupling, 0)),
      as.character(as_num(best_comp$constraint_metrics$joint_constraint_penalty_total, 0)),
      as.character(ctx$joint_weight_invivo),
      as.character(ctx$joint_weight_invitro),
      as.character(ctx$joint_invitro_growth_weight),
      as.character(ctx$joint_invitro_ploidy_weight),
      as.character(ctx$joint_invitro_flow_weight),
      as.character(joint_soft_coupling_active(ctx)),
	      paste(ctx$joint_soft_coupling$params, collapse = ","),
	      as.character(length(ctx$joint_soft_coupling$params)),
	      as.character(ctx$joint_soft_coupling$sigma_default),
	      as.character(isTRUE(ctx$joint_warmup$enabled)),
	      as.character(ctx$joint_warmup$sigmaN),
	      as.character(ctx$joint_warmup$invivo_seed_dir),
	      as.character(ctx$joint_warmup$invitro_seed_dir),
	      as.character(ctx$joint_warmup$invivo_source_path),
	      as.character(ctx$joint_warmup$invitro_source_path),
      as.character(ctx$joint_restriction),
      as.character(ctx$seed),
      as.character(ctx$itermax),
      as.character(optimizer_de_reltol),
      as.character(optimizer_de_steptol),
      as.character(ctx$NP),
      as.character(max(ctx$NP, ctx$joint_np_min_factor * length(ctx$init))),
      as.character(ctx$joint_np_min_factor),
      as.character(ctx$n_cores_requested),
      as.character(ctx$n_cores_used),
      as.character(length(ctx$init)),
      as.character(length(ctx$invivo$scenarios))
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
  if (!is.null(ctx$joint_soft_coupling_start_table$path) &&
      file.exists(ctx$joint_soft_coupling_start_table$path)) {
    file.copy(
      ctx$joint_soft_coupling_start_table$path,
      file.path(out_dir, "joint_soft_coupling_parameters_table_input.csv"),
      overwrite = TRUE
    )
  }
  if (is.data.frame(ctx$joint_soft_coupling_start_table$applied) &&
      nrow(ctx$joint_soft_coupling_start_table$applied) > 0L) {
    write.table(
      ctx$joint_soft_coupling_start_table$applied,
      file = file.path(out_dir, "joint_soft_coupling_initial_values.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  if (is.data.frame(ctx$joint_warmup$applied) &&
      nrow(ctx$joint_warmup$applied) > 0L) {
    write.table(
      ctx$joint_warmup$applied,
      file = file.path(out_dir, "joint_warmup_initial_values.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
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
        joint_soft_coupling = ctx$joint_soft_coupling,
        joint_soft_coupling_start_table = ctx$joint_soft_coupling_start_table,
        joint_warmup = ctx$joint_warmup,
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
  invivo_parameter_table <- trim_cli_scalar_local(cfg_raw$parameter_table)
  if (is.null(invivo_parameter_table)) {
    invivo_parameter_table <- default_o2_parameter_table_path_common(
      script_dir = SCRIPT_DIR,
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
  de_ctrl <- joint_deoptim_control(ctx, NP_use)
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
    ", reltol=", signif(de_ctrl$reltol, 6),
    ", steptol=", de_ctrl$steptol,
    ", NP_min_factor=", ctx$joint_np_min_factor,
    ", n_cores=", ctx$n_cores_used,
    ", weights(invivo,invitro)=(",
    signif(ctx$joint_weight_invivo, 6), ", ",
    signif(ctx$joint_weight_invitro, 6), ")"
  )
  de_fit <- DEoptim::DEoptim(
    fn = objective_value,
    lower = ctx$lower,
    upper = ctx$upper,
    control = de_ctrl
  )
  de_iter_target <- as.integer(de_ctrl$itermax)
  de_iter_completed <- joint_deoptim_iter_completed(de_fit, de_iter_target)
  deoptim_stop_reason <- joint_deoptim_stop_reason(
    iter_completed = de_iter_completed,
    iter_target = de_iter_target,
    interrupted = FALSE
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
        candidate_t <- as.numeric(local_fit$par)
        names(candidate_t) <- names(ctx$init)
        candidate_outside_bounds <- any(!is.finite(candidate_t) | candidate_t < ctx$lower - 1e-8 | candidate_t > ctx$upper + 1e-8)
        if (isTRUE(candidate_outside_bounds)) {
          message("[fit_joint] L-BFGS-B local refinement improved objective but returned an out-of-bounds point; keeping DEoptim best.")
        } else {
          best_t <- candidate_t
          local_accepted <- TRUE
          message(
            "[fit_joint] L-BFGS-B local refinement improved objective: ",
            signif(de_best_objective, 8), " -> ", signif(local_best_objective, 8), "."
          )
        }
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
    interrupted = FALSE,
    iter_completed = as.integer(de_iter_completed),
    iter_target = as.integer(de_iter_target),
    deoptim_stop_reason = deoptim_stop_reason,
    de_reltol = as.numeric(de_ctrl$reltol),
    de_steptol = as.integer(de_ctrl$steptol),
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
  run_prefix <- INVIVO_ENV$.runner_cli_string(.first_non_null_local(cfg$run_prefix, "fit_joint_O2_supply_demand_MAP"))
  if (is.null(run_prefix) || !nzchar(trimws(run_prefix))) {
    run_prefix <- "fit_joint_O2_supply_demand_MAP"
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

  fit_script <- normalizePath(file.path(WORKFLOW_ROOT, "optimizer", "fit_model_O2_supply_demand_MAP.R"), mustWork = FALSE)
  viz_script <- normalizePath(file.path(WORKFLOW_ROOT, "vis", "viz_invivo_model_O2_supply_demand_MAP_results.R"), mustWork = FALSE)
  invitro_viz_script <- normalizePath(file.path(WORKFLOW_ROOT, "vis", "viz_invitro_model_O2_supply_demand_MAP_results.R"), mustWork = FALSE)
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

  log_line("Running O2_supply_demand_MAP joint fit")
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
