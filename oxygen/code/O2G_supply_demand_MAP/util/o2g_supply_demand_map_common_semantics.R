if (!exists("o2sd_first_non_null", mode = "function", inherits = TRUE)) {
  .o2sd_common_script_dir <- local({
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
  source(file.path(.o2sd_common_script_dir, "o2g_supply_demand_map_shared.R"), local = FALSE)
  rm(.o2sd_common_script_dir)
}

# -----------------------------------------------------------------------------
# Function: canonical_ploidy_o2_death_mode
# Purpose: Canonicalize ploidy_O2_death from config/CLI aliases to one mode.
# -----------------------------------------------------------------------------
canonical_ploidy_o2_death_mode <- function(x, default = "diploid_NULL") {
  val <- o2sd_first_non_null(x, default)
  if (is.logical(val) && length(val) > 0L && !is.na(val[[1]])) {
    return(if (isTRUE(val[[1]])) "diploid_NULL" else "uniform")
  }
  s <- tolower(trimws(as.character(val[[1]])))
  if (!nzchar(s)) s <- tolower(trimws(as.character(default[[1]])))
  if (s %in% c("uniform", "false", "f", "0", "no", "n")) return("uniform")
  if (s %in% c("diploid_null", "diploid-null", "diploidnull", "true", "t", "1", "yes", "y")) {
    return("diploid_NULL")
  }
  if (s %in% c("ploidy_related", "ploidy-related", "ploidyrelated")) return("ploidy_related")
  stop(
    "Invalid ploidy_O2_death mode: '", as.character(val[[1]]),
    "'. Allowed values: uniform, diploid_NULL, ploidy_related."
  )
}

# -----------------------------------------------------------------------------
# Function: canonical_start_with_mode
# Purpose: Canonicalize endpoint observation mode for fit/viz dispatch.
# -----------------------------------------------------------------------------
canonical_start_with_mode <- function(x, default = "ploidy") {
  val <- o2sd_first_non_null(x, default)
  s <- tolower(trimws(as.character(val[[1]])))
  if (!nzchar(s)) s <- tolower(trimws(as.character(default[[1]])))
  if (s %in% c("ploidy", "p")) return("ploidy")
  if (s %in% c("chr_number", "chr-number", "chrnumber", "chromosome_number", "chromosome-number", "chromosomenumber", "chromosome", "chromosomes", "n")) {
    return("chr_number")
  }
  stop(
    "Invalid start_with mode: '", as.character(val[[1]]),
    "'. Allowed values: ploidy, chr_number."
  )
}

# -----------------------------------------------------------------------------
# Function: assert_canonical_start_with_mode
# Purpose: Enforce that runtime start_with is already canonical.
# -----------------------------------------------------------------------------
assert_canonical_start_with_mode <- function(x) {
  val <- o2sd_first_non_null(x, NA_character_)
  s <- trimws(as.character(val[[1]]))
  if (!nzchar(s)) {
    stop("start_with must be provided as one of: ploidy, chr_number.")
  }
  if (identical(s, "ploidy") || identical(s, "chr_number")) {
    return(s)
  }
  stop(
    "start_with must already be canonical before runtime dispatch. ",
    "Allowed values: ploidy, chr_number; got '", s, "'."
  )
}

# -----------------------------------------------------------------------------
# Function: assert_canonical_ploidy_o2_death_mode
# Purpose: Enforce that runtime ploidy_O2_death is already canonical.
# -----------------------------------------------------------------------------
assert_canonical_ploidy_o2_death_mode <- function(x) {
  val <- o2sd_first_non_null(x, NA_character_)
  s <- trimws(as.character(val[[1]]))
  if (!nzchar(s)) {
    stop("ploidy_O2_death must be provided as one of: uniform, diploid_NULL, ploidy_related.")
  }
  if (identical(s, "uniform") || identical(s, "diploid_NULL") || identical(s, "ploidy_related")) {
    return(s)
  }
  stop(
    "ploidy_O2_death must already be canonical before runtime dispatch. ",
    "Allowed values: uniform, diploid_NULL, ploidy_related; got '", s, "'."
  )
}

# -----------------------------------------------------------------------------
# Function: canonical_glucose_enabled
# Purpose: Canonicalize the top-level glucose family switch to a scalar boolean.
# -----------------------------------------------------------------------------
canonical_glucose_enabled <- function(x, default = TRUE) {
  val <- o2sd_first_non_null(x, default)
  if (is.logical(val) && length(val) > 0L && !is.na(val[[1]])) {
    return(isTRUE(val[[1]]))
  }
  s <- tolower(trimws(as.character(val[[1]])))
  if (!nzchar(s)) s <- tolower(trimws(as.character(default[[1]])))
  if (s %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
  stop(
    "Invalid glucose value: '", as.character(val[[1]]),
    "'. Allowed values: TRUE/FALSE."
  )
}

# -----------------------------------------------------------------------------
# Function: canonical_glucose_dynamic
# Purpose: Canonicalize the glucose_dynamic config switch to a scalar boolean.
# -----------------------------------------------------------------------------
canonical_glucose_dynamic <- function(x, default = FALSE) {
  val <- o2sd_first_non_null(x, default)
  if (is.logical(val) && length(val) > 0L && !is.na(val[[1]])) {
    return(isTRUE(val[[1]]))
  }
  s <- tolower(trimws(as.character(val[[1]])))
  if (!nzchar(s)) s <- tolower(trimws(as.character(default[[1]])))
  if (s %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
  stop(
    "Invalid glucose_dynamic value: '", as.character(val[[1]]),
    "'. Allowed values: TRUE/FALSE."
  )
}

# -----------------------------------------------------------------------------
# Function: canonical_glucose_stress_mode
# Purpose: Canonicalize glucose stress coupling mode for runtime dispatch.
# -----------------------------------------------------------------------------
canonical_glucose_stress_mode <- function(x, default = "coupled_to_O2") {
  val <- o2sd_first_non_null(x, default)
  if (is.logical(val) && length(val) > 0L && !is.na(val[[1]])) {
    return(if (isTRUE(val[[1]])) "coupled_to_O2" else "off")
  }
  s <- tolower(trimws(as.character(val[[1]])))
  if (!nzchar(s)) s <- tolower(trimws(as.character(default[[1]])))
  if (s %in% c("off", "none", "false", "f", "0", "no", "n")) return("off")
  if (s %in% c("coupled_to_o2", "coupled-to-o2", "coupledtoo2", "coupled", "true", "t", "1", "yes", "y")) {
    return("coupled_to_O2")
  }
  if (s %in% c("dynamic", "dyn", "explicit", "stateful")) return("dynamic")
  stop(
    "Invalid glucose_stress_mode: '", as.character(val[[1]]),
    "'. Allowed values: off, coupled_to_O2, dynamic."
  )
}

# -----------------------------------------------------------------------------
# Function: assert_canonical_glucose_stress_mode
# Purpose: Enforce that runtime glucose_stress_mode is already canonical.
# -----------------------------------------------------------------------------
assert_canonical_glucose_stress_mode <- function(x) {
  val <- o2sd_first_non_null(x, NA_character_)
  s <- trimws(as.character(val[[1]]))
  if (!nzchar(s)) {
    stop("glucose_stress_mode must be provided as one of: off, coupled_to_O2, dynamic.")
  }
  if (identical(s, "off") || identical(s, "coupled_to_O2") || identical(s, "dynamic")) {
    return(s)
  }
  stop(
    "glucose_stress_mode must already be canonical before runtime dispatch. ",
    "Allowed values: off, coupled_to_O2, dynamic; got '", s, "'."
  )
}

# -----------------------------------------------------------------------------
# Function: resolve_glucose_runtime_mode
# Purpose: Derive the runtime glucose-stress mode from the dynamic switch plus
#   any legacy fallback mode.
# -----------------------------------------------------------------------------
resolve_glucose_runtime_mode <- function(glucose_dynamic = NULL,
                                         glucose_stress_mode = NULL,
                                         default_dynamic = FALSE,
                                         default_static_mode = "coupled_to_O2") {
  dynamic_use <- canonical_glucose_dynamic(
    o2sd_first_non_null(glucose_dynamic, default_dynamic),
    default = default_dynamic
  )
  if (isTRUE(dynamic_use)) return("dynamic")

  mode <- canonical_glucose_stress_mode(
    o2sd_first_non_null(glucose_stress_mode, default_static_mode),
    default = default_static_mode
  )
  if (identical(mode, "dynamic")) {
    mode <- canonical_glucose_stress_mode(default_static_mode, default = "coupled_to_O2")
  }
  mode
}

# -----------------------------------------------------------------------------
# Function: default_glucose_ref_mM
# Purpose: Return the normal-blood-glucose reference concentration used to
#   interpret glucose percentages on a 0-100 scale.
# -----------------------------------------------------------------------------
default_glucose_ref_mM <- function() {
  5.0
}

# -----------------------------------------------------------------------------
# Function: default_glucose_pct_scale
# Purpose: Return the default percent-of-normal glucose scale parameters.
# -----------------------------------------------------------------------------
default_glucose_pct_scale <- function() {
  list(
    G_S0 = 100.0,
    kappa_G = 20.0,
    eta_G = 1.2,
    G_c = 30.0
  )
}

# -----------------------------------------------------------------------------
# Function: normalize_glucose_ref_mM
# Purpose: Canonicalize glucose_ref_mM to a positive scalar.
# -----------------------------------------------------------------------------
normalize_glucose_ref_mM <- function(x, default = default_glucose_ref_mM()) {
  val <- suppressWarnings(as.numeric(o2sd_first_non_null(x, default)))
  if (!is.finite(val) || val <= 0) val <- as.numeric(default_glucose_ref_mM())
  val
}

# -----------------------------------------------------------------------------
# Function: glucose_pct_to_mM
# Purpose: Convert glucose as percent-of-normal (0-100) to absolute mM.
# -----------------------------------------------------------------------------
glucose_pct_to_mM <- function(G_pct, glucose_ref_mM = default_glucose_ref_mM()) {
  ref_use <- normalize_glucose_ref_mM(glucose_ref_mM)
  pct <- suppressWarnings(as.numeric(G_pct))
  pct[!is.finite(pct)] <- NA_real_
  ref_use * pct / 100
}

# -----------------------------------------------------------------------------
# Function: read_param_table_prototype_slot_common
# Purpose: Read one natural-scale prototype slot from parameter_table.csv.
# -----------------------------------------------------------------------------
read_param_table_prototype_slot_common <- function(path, param_prototype, slot = c("init", "lower", "upper")) {
  slot <- match.arg(slot)
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NA_real_)
  tab <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab)) return(NA_real_)

  if (all(c("param_symbol", "init_value", "lower_bound", "upper_bound") %in% names(tab))) {
    col_map <- c(
      init = "init_value",
      lower = "lower_bound",
      upper = "upper_bound"
    )
    proto <- trimws(as.character(tab$param_symbol))
  } else {
    col_map <- c(
      init = "prototype_init_value",
      lower = "prototype_lower_bound",
      upper = "prototype_upper_bound"
    )
    if (!all(c("param_prototype", unname(col_map)) %in% names(tab))) return(NA_real_)
    proto <- trimws(as.character(tab$param_prototype))
  }

  idx <- match(param_prototype, proto)
  if (is.na(idx)) return(NA_real_)
  suppressWarnings(as.numeric(tab[[col_map[[slot]]]][idx]))
}

# -----------------------------------------------------------------------------
# Function: read_o2_S0_natural_upper_bound_common
# Purpose: Read natural-scale o2_S0 upper bound from parameter_table.csv.
# -----------------------------------------------------------------------------
read_o2_S0_natural_upper_bound_common <- function(path, fallback = 5.0) {
  ub <- read_param_table_prototype_slot_common(path, "o2_S0", slot = "upper")
  if (!is.finite(ub) || ub <= 0) ub <- as.numeric(fallback)
  if (!is.finite(ub) || ub <= 0) ub <- 5.0
  ub
}

# -----------------------------------------------------------------------------
# Function: normalize_sim_cfg_common
# Purpose: Canonicalize shared simulation config semantics for fit/viz.
# -----------------------------------------------------------------------------
normalize_sim_cfg_common <- function(cfg, context = c("fit", "viz")) {
  context <- match.arg(context)
  if (is.null(cfg)) cfg <- list()

  cfg$N_UNIT <- as.integer(o2sd_first_non_null(cfg$N_UNIT, 22L))
  cfg$N_MIN <- as.integer(o2sd_first_non_null(cfg$N_MIN, 22L))
  cfg$N_MAX <- as.integer(o2sd_first_non_null(cfg$N_MAX, 154L))
  cfg$DT <- as.numeric(o2sd_first_non_null(cfg$DT, 0.5))

  cfg$o2_S0_upper_bound <- as.numeric(o2sd_first_non_null(
    cfg$o2_S0_upper_bound,
    read_o2_S0_natural_upper_bound_common(cfg$parameter_table, fallback = 5.0)
  ))
  if (!is.finite(cfg$o2_S0_upper_bound) || cfg$o2_S0_upper_bound <= 0) {
    stop(context, "_config o2_S0_upper_bound must be > 0.")
  }

  cfg$o2_min <- as.numeric(o2sd_first_non_null(cfg$o2_min, 0.5))
  if (!is.finite(cfg$o2_min) || cfg$o2_min < 0) cfg$o2_min <- 0.5
  cfg$o2_min <- min(max(cfg$o2_min, 0), cfg$o2_S0_upper_bound)

  cfg$init_total_size <- as.numeric(o2sd_first_non_null(cfg$init_total_size, 1e6))
  cfg$o2_Nref <- as.numeric(o2sd_first_non_null(cfg$o2_Nref, cfg$init_total_size, 1e6))
  if (!is.finite(cfg$o2_Nref) || cfg$o2_Nref <= 0) cfg$o2_Nref <- 1e6

  cfg$tau_O2_init <- as.numeric(o2sd_first_non_null(cfg$tau_O2_init, 2.0))
  cfg$tau_O2 <- as.numeric(o2sd_first_non_null(cfg$tau_O2, cfg$tau_O2_init))
  if (!is.finite(cfg$tau_O2) || cfg$tau_O2 <= 0) cfg$tau_O2 <- cfg$tau_O2_init

  cfg$K <- as.numeric(o2sd_first_non_null(cfg$K, 1e12))
  if (!is.finite(cfg$K) || cfg$K <= 0) cfg$K <- 1e12
  cfg$crowding <- as.character(o2sd_first_non_null(cfg$crowding, "logistic"))
  cfg$Crowding <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$Crowding, TRUE), TRUE)
  if (!cfg$crowding %in% c("logistic", "gompertz")) {
    stop(context, "_config crowding must be logistic or gompertz.")
  }

  cfg$dose_ref <- as.numeric(o2sd_first_non_null(cfg$dose_ref, 30))
  cfg$tx_mult_min <- as.numeric(o2sd_first_non_null(cfg$tx_mult_min, 0.05))
  cfg$min_pop <- as.numeric(o2sd_first_non_null(cfg$min_pop, 1e-12))
  cfg$alpha_o2_init <- as.numeric(o2sd_first_non_null(cfg$alpha_o2_init, 0.5))
  cfg$gamma_growth_init <- as.numeric(o2sd_first_non_null(cfg$gamma_growth_init, 2.0))
  cfg$mu_hp_init <- as.numeric(o2sd_first_non_null(cfg$mu_hp_init, 1e-3))
  cfg$gamma_mu_init <- as.numeric(o2sd_first_non_null(cfg$gamma_mu_init, 1.0))
  cfg$o2_crit_init <- as.numeric(o2sd_first_non_null(cfg$o2_crit_init, 1.0))
  cfg$n_O_init <- as.numeric(o2sd_first_non_null(cfg$n_O_init, 1.0))
  cfg$k_clear_init <- as.numeric(o2sd_first_non_null(cfg$k_clear_init, 1e-3))
  cfg$p_wgd_init <- as.numeric(o2sd_first_non_null(cfg$p_wgd_init, 1e-4))
  cfg$p_wgd_max_init <- as.numeric(o2sd_first_non_null(cfg$p_wgd_max_init, 1e-3))
  cfg$O2_wgd_init <- as.numeric(o2sd_first_non_null(cfg$O2_wgd_init, 0.1))
  cfg$harvest_init_multiplier <- o2sd_as_bool_scalar(
    o2sd_first_non_null(cfg$harvest_init_multiplier, TRUE),
    TRUE
  )
  cfg$prior_center_log_init_mult <- as.numeric(o2sd_first_non_null(cfg$prior_center_log_init_mult, 0.0))
  cfg$prior_sd_log_init_mult <- as.numeric(o2sd_first_non_null(cfg$prior_sd_log_init_mult, 0.35))
  cfg$log_init_mult_lower <- as.numeric(o2sd_first_non_null(cfg$log_init_mult_lower, -1.0))
  cfg$log_init_mult_upper <- as.numeric(o2sd_first_non_null(cfg$log_init_mult_upper, 1.0))
  cfg$glucose <- canonical_glucose_enabled(
    o2sd_first_non_null(cfg$glucose, TRUE),
    default = TRUE
  )
  cfg$glucose_dynamic <- canonical_glucose_dynamic(
    o2sd_first_non_null(cfg$glucose_dynamic, FALSE),
    default = FALSE
  )
  if (!isTRUE(cfg$glucose)) {
    cfg$glucose_dynamic <- FALSE
  }
  cfg$ploidy_O2_death <- canonical_ploidy_o2_death_mode(
    o2sd_first_non_null(cfg$ploidy_O2_death, "diploid_NULL"),
    default = "diploid_NULL"
  )
  cfg$glucose_stress_mode <- resolve_glucose_runtime_mode(
    glucose_dynamic = cfg$glucose_dynamic,
    glucose_stress_mode = if (isTRUE(cfg$glucose)) {
      o2sd_first_non_null(cfg$glucose_stress_mode, "coupled_to_O2")
    } else {
      "off"
    },
    default_dynamic = FALSE,
    default_static_mode = if (isTRUE(cfg$glucose)) "coupled_to_O2" else "off"
  )
  cfg$start_with <- canonical_start_with_mode(
    o2sd_first_non_null(cfg$start_with, "ploidy"),
    default = "ploidy"
  )

  if (!is.finite(cfg$mu_hp_init) || cfg$mu_hp_init <= 0) cfg$mu_hp_init <- 1e-3
  if (!is.finite(cfg$gamma_mu_init) || cfg$gamma_mu_init <= 0) cfg$gamma_mu_init <- 1.0
  if (!is.finite(cfg$o2_crit_init) || cfg$o2_crit_init < 0) cfg$o2_crit_init <- 1.0
  if (!is.finite(cfg$n_O_init) || cfg$n_O_init < 0) cfg$n_O_init <- 1.0
  if (!is.finite(cfg$k_clear_init) || cfg$k_clear_init <= 0) cfg$k_clear_init <- 1e-3
  if (!is.finite(cfg$p_wgd_init) || cfg$p_wgd_init <= 0) cfg$p_wgd_init <- 1e-4
  if (!is.finite(cfg$p_wgd_max_init) || cfg$p_wgd_max_init <= 0) cfg$p_wgd_max_init <- 1e-3
  if (!is.finite(cfg$O2_wgd_init) || cfg$O2_wgd_init <= 0) cfg$O2_wgd_init <- 0.1
  if (!is.finite(cfg$prior_center_log_init_mult)) cfg$prior_center_log_init_mult <- 0.0
  if (!is.finite(cfg$prior_sd_log_init_mult) || cfg$prior_sd_log_init_mult <= 0) cfg$prior_sd_log_init_mult <- 0.35
  if (!is.finite(cfg$log_init_mult_lower)) cfg$log_init_mult_lower <- -1.0
  if (!is.finite(cfg$log_init_mult_upper)) cfg$log_init_mult_upper <- 1.0
  if (cfg$log_init_mult_upper < cfg$log_init_mult_lower) {
    tmp <- cfg$log_init_mult_lower
    cfg$log_init_mult_lower <- cfg$log_init_mult_upper
    cfg$log_init_mult_upper <- tmp
  }

  glucose_defaults <- default_glucose_pct_scale()
  cfg$glucose_ref_mM <- normalize_glucose_ref_mM(o2sd_first_non_null(cfg$glucose_ref_mM, default_glucose_ref_mM()))
  cfg$G_S0_init <- as.numeric(o2sd_first_non_null(cfg$G_S0_init, glucose_defaults$G_S0))
  cfg$kappa_G_init <- as.numeric(o2sd_first_non_null(cfg$kappa_G_init, glucose_defaults$kappa_G))
  cfg$eta_G_init <- as.numeric(o2sd_first_non_null(cfg$eta_G_init, glucose_defaults$eta_G))
  cfg$G_c_init <- as.numeric(o2sd_first_non_null(cfg$G_c_init, glucose_defaults$G_c))
  cfg$tau_G_init <- as.numeric(o2sd_first_non_null(cfg$tau_G_init, cfg$tau_O2_init, cfg$tau_O2, 0.1))
  cfg$tau_G <- as.numeric(o2sd_first_non_null(cfg$tau_G, cfg$tau_G_init))
  if (!is.finite(cfg$G_S0_init) || cfg$G_S0_init <= 0) cfg$G_S0_init <- glucose_defaults$G_S0
  if (!is.finite(cfg$kappa_G_init) || cfg$kappa_G_init <= 0) cfg$kappa_G_init <- glucose_defaults$kappa_G
  if (!is.finite(cfg$eta_G_init) || cfg$eta_G_init <= 0) cfg$eta_G_init <- glucose_defaults$eta_G
  if (!is.finite(cfg$G_c_init) || cfg$G_c_init <= 0) cfg$G_c_init <- glucose_defaults$G_c
  if (!is.finite(cfg$tau_G_init) || cfg$tau_G_init <= 0) cfg$tau_G_init <- 0.1
  if (!is.finite(cfg$tau_G) || cfg$tau_G <= 0) cfg$tau_G <- cfg$tau_G_init

  cfg$dose_zero_only <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$dose_zero_only, TRUE), TRUE)
  cfg$fit_treatment <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$fit_treatment, FALSE), FALSE)
  cfg$max_scenarios <- as.numeric(o2sd_first_non_null(cfg$max_scenarios, Inf))
  cfg$o2_burden_feedback <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$o2_burden_feedback, TRUE), TRUE)
  cfg$O2_growth <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$O2_growth, TRUE), TRUE)
  cfg$o2_cache_profile <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$o2_cache_profile, FALSE), FALSE)

  if (is.null(cfg$truncate_at_treatment)) {
    cfg$truncate_at_treatment <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$pretreat_only, FALSE), FALSE)
  }
  cfg$truncate_at_treatment <- o2sd_as_bool_scalar(cfg$truncate_at_treatment, FALSE)

  if (is.null(cfg$ploidy_at_harvest)) cfg$ploidy_at_harvest <- TRUE
  cfg$ploidy_at_harvest <- o2sd_as_bool_scalar(cfg$ploidy_at_harvest, TRUE)

  cfg
}

# -----------------------------------------------------------------------------
# Function: normalize_run_params_common
# Purpose: Fill shared run_params defaults/canonical values after reconstruction.
# -----------------------------------------------------------------------------
normalize_run_params_common <- function(run_params, cfg = NULL) {
  if (is.null(run_params)) run_params <- list()
  if (is.null(cfg)) cfg <- list()

  run_params$p_mis_base <- as.numeric(o2sd_first_non_null(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5))
  run_params$harvest_init_multiplier <- o2sd_as_bool_scalar(
    o2sd_first_non_null(run_params$harvest_init_multiplier, cfg$harvest_init_multiplier, TRUE),
    TRUE
  )
  run_params$glucose <- canonical_glucose_enabled(
    o2sd_first_non_null(run_params$glucose, cfg$glucose, TRUE),
    default = TRUE
  )
  run_params$glucose_dynamic <- canonical_glucose_dynamic(
    o2sd_first_non_null(run_params$glucose_dynamic, cfg$glucose_dynamic, FALSE),
    default = FALSE
  )
  if (!isTRUE(run_params$glucose)) {
    run_params$glucose_dynamic <- FALSE
  }
  run_params$ploidy_O2_death <- canonical_ploidy_o2_death_mode(
    o2sd_first_non_null(run_params$ploidy_O2_death, cfg$ploidy_O2_death, "diploid_NULL"),
    default = "diploid_NULL"
  )
  run_params$glucose_stress_mode <- resolve_glucose_runtime_mode(
    glucose_dynamic = run_params$glucose_dynamic,
    glucose_stress_mode = if (isTRUE(run_params$glucose)) {
      o2sd_first_non_null(run_params$glucose_stress_mode, cfg$glucose_stress_mode, "coupled_to_O2")
    } else {
      "off"
    },
    default_dynamic = cfg$glucose_dynamic,
    default_static_mode = if (isTRUE(run_params$glucose)) "coupled_to_O2" else "off"
  )
  run_params$start_with <- canonical_start_with_mode(
    o2sd_first_non_null(run_params$start_with, cfg$start_with, "ploidy"),
    default = "ploidy"
  )

  run_params$o2_S0_upper_bound <- as.numeric(o2sd_first_non_null(
    run_params$o2_S0_upper_bound,
    cfg$o2_S0_upper_bound,
    read_o2_S0_natural_upper_bound_common(cfg$parameter_table, fallback = 5.0)
  ))
  if (!is.finite(run_params$o2_S0_upper_bound) || run_params$o2_S0_upper_bound <= 0) {
    run_params$o2_S0_upper_bound <- 5.0
  }

  run_params$o2_min <- as.numeric(o2sd_first_non_null(run_params$o2_min, cfg$o2_min, 0.5))
  if (!is.finite(run_params$o2_min) || run_params$o2_min < 0) run_params$o2_min <- 0.5
  run_params$o2_min <- min(max(run_params$o2_min, 0), run_params$o2_S0_upper_bound)

  run_params$o2_Nref <- as.numeric(o2sd_first_non_null(run_params$o2_Nref, cfg$o2_Nref, cfg$init_total_size, 1e6))
  if (!is.finite(run_params$o2_Nref) || run_params$o2_Nref <= 0) run_params$o2_Nref <- 1e6

  fit_treatment_use <- isTRUE(o2sd_first_non_null(cfg$fit_treatment, FALSE))
  if (fit_treatment_use) {
    if (is.null(run_params$alpha) || !is.finite(as.numeric(run_params$alpha))) {
      stop("run_params$alpha must be present when fit_treatment=TRUE.")
    }
    if (is.null(run_params$gamma) || !is.finite(as.numeric(run_params$gamma))) {
      stop("run_params$gamma must be present when fit_treatment=TRUE.")
    }
    run_params$alpha <- as.numeric(run_params$alpha)
    run_params$gamma <- as.numeric(run_params$gamma)
  } else {
    run_params$alpha <- as.numeric(o2sd_first_non_null(run_params$alpha, 0))
    run_params$gamma <- as.numeric(o2sd_first_non_null(run_params$gamma, 1))
  }

  tau_use <- as.numeric(run_params$tau_O2)
  if (!is.finite(tau_use) || tau_use <= 0) {
    tau_use <- as.numeric(o2sd_first_non_null(cfg$tau_O2, cfg$tau_O2_init, 2.0))
  }
  if (!is.finite(tau_use) || tau_use <= 0) tau_use <- 2.0
  run_params$tau_O2 <- tau_use

  glucose_defaults <- default_glucose_pct_scale()
  run_params$p_wgd <- as.numeric(o2sd_first_non_null(run_params$p_wgd, cfg$p_wgd_init, 1e-4))
  if (!is.finite(run_params$p_wgd) || run_params$p_wgd < 0) run_params$p_wgd <- 0.0
  run_params$glucose_ref_mM <- normalize_glucose_ref_mM(
    o2sd_first_non_null(run_params$glucose_ref_mM, cfg$glucose_ref_mM, default_glucose_ref_mM())
  )
  run_params$G_S0 <- as.numeric(o2sd_first_non_null(run_params$G_S0, cfg$G_S0_init, glucose_defaults$G_S0))
  run_params$kappa_G <- as.numeric(o2sd_first_non_null(run_params$kappa_G, cfg$kappa_G_init, glucose_defaults$kappa_G))
  run_params$eta_G <- as.numeric(o2sd_first_non_null(run_params$eta_G, cfg$eta_G_init, glucose_defaults$eta_G))
  run_params$G_c <- as.numeric(o2sd_first_non_null(run_params$G_c, cfg$G_c_init, glucose_defaults$G_c))
  run_params$tau_G <- as.numeric(o2sd_first_non_null(run_params$tau_G, cfg$tau_G, cfg$tau_G_init, cfg$tau_O2, cfg$tau_O2_init, 0.1))
  if (!is.finite(run_params$G_S0) || run_params$G_S0 <= 0) run_params$G_S0 <- glucose_defaults$G_S0
  if (!is.finite(run_params$kappa_G) || run_params$kappa_G <= 0) run_params$kappa_G <- glucose_defaults$kappa_G
  if (!is.finite(run_params$eta_G) || run_params$eta_G <= 0) run_params$eta_G <- glucose_defaults$eta_G
  if (!is.finite(run_params$G_c) || run_params$G_c <= 0) run_params$G_c <- glucose_defaults$G_c
  if (!is.finite(run_params$tau_G) || run_params$tau_G <= 0) run_params$tau_G <- 0.1

  run_params
}
