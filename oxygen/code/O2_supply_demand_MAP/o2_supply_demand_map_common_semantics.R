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
  source(file.path(.o2sd_common_script_dir, "o2_supply_demand_map_shared.R"), local = FALSE)
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
  cfg$ploidy_O2_death <- canonical_ploidy_o2_death_mode(
    o2sd_first_non_null(cfg$ploidy_O2_death, "diploid_NULL"),
    default = "diploid_NULL"
  )

  if (!is.finite(cfg$mu_hp_init) || cfg$mu_hp_init <= 0) cfg$mu_hp_init <- 1e-3
  if (!is.finite(cfg$gamma_mu_init) || cfg$gamma_mu_init <= 0) cfg$gamma_mu_init <- 1.0
  if (!is.finite(cfg$o2_crit_init) || cfg$o2_crit_init < 0) cfg$o2_crit_init <- 1.0
  if (!is.finite(cfg$n_O_init) || cfg$n_O_init < 0) cfg$n_O_init <- 1.0
  if (!is.finite(cfg$k_clear_init) || cfg$k_clear_init <= 0) cfg$k_clear_init <- 1e-3

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
  run_params$ploidy_O2_death <- canonical_ploidy_o2_death_mode(
    o2sd_first_non_null(run_params$ploidy_O2_death, cfg$ploidy_O2_death, "diploid_NULL"),
    default = "diploid_NULL"
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

  run_params
}
