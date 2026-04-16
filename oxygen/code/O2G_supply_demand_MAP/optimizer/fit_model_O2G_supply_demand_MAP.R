#!/usr/bin/env Rscript

.o2g_fit_model_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2g_fit_model_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2g_fit_model_bootstrap_script_dir)

parse_args <- o2sd_parse_args
as_bool <- o2sd_as_bool
as_num <- o2sd_as_num
as_int <- o2sd_as_int
.first_non_null_local <- o2sd_first_non_null

fail_fit_input <- function(mode_label, ...) {
  stop(paste0(mode_label, " input error: ", paste0(..., collapse = "")), call. = FALSE)
}

trim_cli_scalar <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  txt <- trimws(as.character(x[[1]]))
  if (!nzchar(txt)) return(NULL)
  txt
}

require_exactly_one_fit_mode <- function(argv) {
  has_invivo <- isTRUE(as_bool(argv$fit_invivo, FALSE))
  has_invitro <- isTRUE(as_bool(argv$fit_invitro, FALSE))
  if (has_invivo && has_invitro) {
    stop("Ambiguous mode flags: pass exactly one of --fit_invivo or --fit_invitro, not both.", call. = FALSE)
  }
  if (!has_invivo && !has_invitro) {
    stop("Missing required mode flag: pass exactly one of --fit_invivo or --fit_invitro.", call. = FALSE)
  }
  if (has_invivo) "fit_invivo" else "fit_invitro"
}

strip_fit_mode_flags <- function(argv) {
  argv$fit_invivo <- NULL
  argv$fit_invitro <- NULL
  argv
}

load_backend_env <- local({
  cache <- new.env(parent = emptyenv())

  function(mode_name) {
    key <- if (identical(mode_name, "fit_invivo")) "invivo" else "invitro"
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    script_path <- if (identical(mode_name, "fit_invivo")) {
      file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_fit_invivo_backend.R")
    } else {
      file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_fit_invitro_backend.R")
    }
    if (!file.exists(script_path)) {
      stop("Backend script not found: ", script_path, call. = FALSE)
    }

    env <- new.env(parent = globalenv())
    sys.source(script_path, envir = env, chdir = TRUE)
    assign(key, env, envir = cache)
    env
  }
})

default_invivo_parameter_table_path <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "data", "O2G_supply_demand", "parameter_table.csv"), mustWork = FALSE)
}

default_invivo_data_dir <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
}

resolve_invivo_backend_mode <- function(argv) {
  mode_raw <- trim_cli_scalar(argv$mode)
  if (is.null(mode_raw)) {
    if (!is.null(trim_cli_scalar(argv$config))) {
      return("run")
    }
    return("fit_seed")
  }
  mode_use <- tolower(mode_raw)
  if (!mode_use %in% c("run", "runner", "fit_seed", "fit", "single_seed")) {
    stop(
      "Unsupported --mode for --fit_invivo: ", mode_use,
      ". Expected one of: run, fit_seed.",
      call. = FALSE
    )
  }
  if (mode_use %in% c("run", "runner")) "run" else "fit_seed"
}

validate_parameter_table_for_mode <- function(backend_env,
                                              parameter_table,
                                              mode_label,
                                              fit_treatment = FALSE,
                                              fit_tau_O2 = FALSE,
                                              O2_growth = TRUE,
                                              glucose = TRUE,
                                              glucose_dynamic = FALSE) {
  if (!file.exists(parameter_table)) {
    fail_fit_input(mode_label, "parameter_table not found: ", parameter_table)
  }
  tryCatch(
    backend_env$build_transformed_parameter_table(
      path = parameter_table,
      fit_treatment = fit_treatment,
      fit_tau_O2 = fit_tau_O2,
      O2_growth = O2_growth,
      glucose = glucose,
      glucose_dynamic = glucose_dynamic
    ),
    error = function(e) {
      fail_fit_input(mode_label, "parameter_table is invalid: ", conditionMessage(e))
    }
  )
  invisible(normalizePath(parameter_table, mustWork = FALSE))
}

validate_invivo_dataset <- function(backend_env,
                                    data_dir,
                                    cfg_probe,
                                    mode_label = "fit_invivo") {
  if (!dir.exists(data_dir)) {
    fail_fit_input(mode_label, "data_dir not found: ", data_dir)
  }
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  if (!file.exists(dt_path)) {
    fail_fit_input(
      mode_label,
      "missing tumor-burden workbook: ", dt_path,
      ". Expected file name: dt_Gem_VT_20260209_v5.xlsx"
    )
  }
  ploidy_path <- tryCatch(
    backend_env$resolve_terminal_ploidy_path(data_dir),
    error = function(e) {
      fail_fit_input(mode_label, conditionMessage(e))
    }
  )
  tryCatch(
    backend_env$prepare_data(dt_path, ploidy_path, cfg_probe),
    error = function(e) {
      fail_fit_input(mode_label, conditionMessage(e))
    }
  )
  invisible(list(
    data_dir = normalizePath(data_dir, mustWork = FALSE),
    dt_path = normalizePath(dt_path, mustWork = FALSE),
    ploidy_path = normalizePath(ploidy_path, mustWork = FALSE)
  ))
}

validate_fit_invivo_inputs <- function(argv, backend_env) {
  invivo_mode <- resolve_invivo_backend_mode(argv)

  if (identical(invivo_mode, "run")) {
    parsed <- tryCatch(
      backend_env$.runner_resolve_config(argv = argv, script_dir = SCRIPT_DIR, caller_wd = getwd()),
      error = function(e) {
        fail_fit_input("fit_invivo", conditionMessage(e))
      }
    )
    cfg_raw <- parsed$cfg
    parameter_table <- trim_cli_scalar(cfg_raw$parameter_table)
    if (is.null(parameter_table)) {
      parameter_table <- default_invivo_parameter_table_path()
    }
    glucose_dynamic_use <- isTRUE(canonical_glucose_dynamic(
      .first_non_null_local(cfg_raw$glucose_dynamic, FALSE),
      default = FALSE
    ))
    glucose_use <- isTRUE(canonical_glucose_enabled(
      .first_non_null_local(cfg_raw$glucose, TRUE),
      default = TRUE
    ))
    validate_parameter_table_for_mode(
      backend_env = backend_env,
      parameter_table = parameter_table,
      mode_label = "fit_invivo",
      fit_treatment = isTRUE(as_bool(cfg_raw$fit_treatment, FALSE)),
      fit_tau_O2 = isTRUE(as_bool(cfg_raw$fit_tau_O2, FALSE)),
      O2_growth = isTRUE(as_bool(cfg_raw$O2_growth, TRUE)),
      glucose = glucose_use,
      glucose_dynamic = glucose_dynamic_use
    )

    data_dir <- trim_cli_scalar(cfg_raw$data_dir)
    if (is.null(data_dir)) {
      fail_fit_input("fit_invivo", "data_dir is required for --fit_invivo --mode=run.")
    }
    cfg_probe <- list(
      dose_zero_only = isTRUE(as_bool(cfg_raw$dose_zero_only, TRUE)),
      truncate_at_treatment = isTRUE(as_bool(cfg_raw$truncate_at_treatment, FALSE)),
      ploidy_at_harvest = isTRUE(as_bool(cfg_raw$ploidy_at_harvest, TRUE)),
      paired_only = isTRUE(as_bool(cfg_raw$paired_only, TRUE)),
      start_with = canonical_start_with_mode(cfg_raw$start_with, "ploidy"),
      N_UNIT = as_int(cfg_raw$N_UNIT, 22L),
      max_scenarios = as_num(cfg_raw$max_scenarios, Inf)
    )
    validate_invivo_dataset(backend_env, data_dir = data_dir, cfg_probe = cfg_probe, mode_label = "fit_invivo")
    return(invisible(list(mode = "run", parsed = parsed)))
  }

  parameter_table <- trim_cli_scalar(.first_non_null_local(argv$parameter_table, argv$parameters))
  if (is.null(parameter_table)) {
    parameter_table <- default_invivo_parameter_table_path()
  }
  glucose_dynamic_use <- isTRUE(canonical_glucose_dynamic(
    .first_non_null_local(argv$glucose_dynamic, FALSE),
    default = FALSE
  ))
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(argv$glucose, TRUE),
    default = TRUE
  ))
  validate_parameter_table_for_mode(
    backend_env = backend_env,
    parameter_table = parameter_table,
    mode_label = "fit_invivo",
    fit_treatment = isTRUE(as_bool(argv$fit_treatment, FALSE)),
    fit_tau_O2 = isTRUE(as_bool(argv$fit_tau_O2, FALSE)),
    O2_growth = isTRUE(as_bool(argv$O2_growth, TRUE)),
    glucose = glucose_use,
    glucose_dynamic = glucose_dynamic_use
  )

  data_dir <- trim_cli_scalar(argv$data_dir)
  if (is.null(data_dir)) {
    data_dir <- default_invivo_data_dir()
  }
  cfg_probe <- list(
    dose_zero_only = isTRUE(as_bool(argv$dose_zero_only, TRUE)),
    truncate_at_treatment = isTRUE(as_bool(argv$truncate_at_treatment, FALSE)),
    ploidy_at_harvest = isTRUE(as_bool(argv$ploidy_at_harvest, TRUE)),
    paired_only = isTRUE(as_bool(argv$paired_only, TRUE)),
    start_with = canonical_start_with_mode(argv$start_with, "ploidy"),
    N_UNIT = as_int(argv$N_UNIT, 22L),
    max_scenarios = as_num(argv$max_scenarios, Inf)
  )
  validate_invivo_dataset(backend_env, data_dir = data_dir, cfg_probe = cfg_probe, mode_label = "fit_invivo")
  invisible(list(mode = "fit_seed"))
}

validate_invitro_observation_tables <- function(x_data, growth_data) {
  required_ploidy_passages <- expand.grid(
    ploidy = c("2N", "4N"),
    passage = c(0L, 7L, 17L),
    stringsAsFactors = FALSE
  )
  observed_pp <- unique(x_data[, c("ploidy", "passage")])
  missing_pp <- required_ploidy_passages[
    !paste(required_ploidy_passages$ploidy, required_ploidy_passages$passage) %in%
      paste(observed_pp$ploidy, observed_pp$passage),
    ,
    drop = FALSE
  ]
  if (nrow(missing_pp) > 0L) {
    missing_txt <- paste(paste0(missing_pp$ploidy, "@passage", missing_pp$passage), collapse = ", ")
    fail_fit_input("fit_invitro", "ploidy distribution data missing required ploidy/passage observations: ", missing_txt)
  }

  required_growth_pairs <- expand.grid(
    ploidy = c("2N", "4N"),
    o2 = c(0, 1),
    stringsAsFactors = FALSE
  )
  observed_growth_pairs <- unique(growth_data[, c("ploidy", "o2")])
  missing_growth_pairs <- required_growth_pairs[
    !paste(required_growth_pairs$ploidy, required_growth_pairs$o2) %in%
      paste(observed_growth_pairs$ploidy, observed_growth_pairs$o2),
    ,
    drop = FALSE
  ]
  if (nrow(missing_growth_pairs) > 0L) {
    missing_txt <- paste(paste0(missing_growth_pairs$ploidy, "@o2=", missing_growth_pairs$o2), collapse = ", ")
    fail_fit_input("fit_invitro", "growth data missing required ploidy/o2 combinations: ", missing_txt)
  }

  if (!any(is.finite(growth_data$passage) & growth_data$passage >= 1L)) {
    fail_fit_input("fit_invitro", "growth data must contain at least one finite passage >= 1.")
  }
  invisible(TRUE)
}

validate_fit_invitro_inputs <- function(argv, backend_env) {
  parameter_table <- trim_cli_scalar(argv$parameter_table)
  if (is.null(parameter_table)) {
    parameter_table <- backend_env$default_parameter_table_path()
  }
  glucose_dynamic_use <- isTRUE(canonical_glucose_dynamic(
    .first_non_null_local(argv$glucose_dynamic, FALSE),
    default = FALSE
  ))
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(argv$glucose, TRUE),
    default = TRUE
  ))
  validate_parameter_table_for_mode(
    backend_env = backend_env,
    parameter_table = parameter_table,
    mode_label = "fit_invitro",
    fit_treatment = FALSE,
    fit_tau_O2 = FALSE,
    O2_growth = TRUE,
    glucose = glucose_use,
    glucose_dynamic = glucose_dynamic_use
  )

  x_data_path <- trim_cli_scalar(argv$x_data)
  if (is.null(x_data_path)) {
    x_data_path <- backend_env$default_ploidy_data_path()
  }
  growth_data_path <- trim_cli_scalar(argv$growth_data)
  if (is.null(growth_data_path)) {
    growth_data_path <- backend_env$default_growth_data_path()
  }
  if (!file.exists(x_data_path)) {
    fail_fit_input("fit_invitro", "x_data not found: ", x_data_path)
  }
  if (!file.exists(growth_data_path)) {
    fail_fit_input("fit_invitro", "growth_data not found: ", growth_data_path)
  }

  x_data <- tryCatch(
    backend_env$load_ploidy_distribution_data(x_data_path),
    error = function(e) {
      fail_fit_input("fit_invitro", conditionMessage(e))
    }
  )
  growth_data <- tryCatch(
    backend_env$load_growth_data(growth_data_path),
    error = function(e) {
      fail_fit_input("fit_invitro", conditionMessage(e))
    }
  )
  validate_invitro_observation_tables(x_data, growth_data)
  invisible(list(
    parameter_table = normalizePath(parameter_table, mustWork = FALSE),
    x_data = normalizePath(x_data_path, mustWork = FALSE),
    growth_data = normalizePath(growth_data_path, mustWork = FALSE)
  ))
}

validate_dispatch_argv <- function(argv, fit_mode, backend_env) {
  if (identical(fit_mode, "fit_invivo")) {
    return(validate_fit_invivo_inputs(argv, backend_env))
  }
  validate_fit_invitro_inputs(argv, backend_env)
}

dispatch_backend_main <- function(argv, fit_mode, backend_env) {
  dispatch_argv <- strip_fit_mode_flags(argv)
  if (!identical(fit_mode, "fit_invivo")) {
    dispatch_argv$mode <- NULL
  }
  if (!exists("main", envir = backend_env, inherits = FALSE) || !is.function(backend_env$main)) {
    stop("Backend entrypoint 'main' not found for ", fit_mode, ".", call. = FALSE)
  }
  backend_env$main(dispatch_argv)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_mode <- require_exactly_one_fit_mode(argv)
  backend_env <- load_backend_env(fit_mode)
  validate_dispatch_argv(strip_fit_mode_flags(argv), fit_mode, backend_env)
  dispatch_backend_main(argv, fit_mode, backend_env)
}

if (sys.nframe() == 0) {
  main()
}
