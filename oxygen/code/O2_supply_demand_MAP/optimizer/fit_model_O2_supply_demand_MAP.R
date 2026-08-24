#!/usr/bin/env Rscript

.o2_fit_model_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2_fit_model_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2_fit_model_bootstrap_script_dir)

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
  has_joint <- isTRUE(as_bool(argv$fit_joint, FALSE))
  n_modes <- sum(c(has_invivo, has_invitro, has_joint))
  if (n_modes > 1L) {
    stop("Ambiguous mode flags: pass exactly one of --fit_invivo, --fit_invitro, or --fit_joint.", call. = FALSE)
  }
  if (n_modes < 1L) {
    stop("Missing required mode flag: pass exactly one of --fit_invivo, --fit_invitro, or --fit_joint.", call. = FALSE)
  }
  if (has_invivo) {
    "fit_invivo"
  } else if (has_invitro) {
    "fit_invitro"
  } else {
    "fit_joint"
  }
}

strip_fit_mode_flags <- function(argv) {
  argv$fit_invivo <- NULL
  argv$fit_invitro <- NULL
  argv$fit_joint <- NULL
  argv
}

load_backend_env <- local({
  cache <- new.env(parent = emptyenv())

  function(mode_name) {
    key <- switch(
      mode_name,
      fit_invivo = "invivo",
      fit_invitro = "invitro",
      fit_joint = "joint",
      stop("Unsupported fit mode: ", mode_name, call. = FALSE)
    )
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    script_path <- switch(
      mode_name,
      fit_invivo = file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_invivo_backend.R"),
      fit_invitro = file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_invitro_backend.R"),
      fit_joint = file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_joint_backend.R")
    )
    if (!file.exists(script_path)) {
      stop("Backend script not found: ", script_path, call. = FALSE)
    }

    env <- new.env(parent = globalenv())
    sys.source(script_path, envir = env, chdir = TRUE)
    assign(key, env, envir = cache)
    env
  }
})

default_invivo_parameter_table_path <- function(script_dir = SCRIPT_DIR, must_exist = FALSE) {
  default_o2_parameter_table_path_common(
    script_dir = script_dir,
    must_exist = must_exist
  )
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
                                              O2_growth = TRUE) {
  if (!file.exists(parameter_table)) {
    fail_fit_input(mode_label, "parameter_table not found: ", parameter_table)
  }
  tryCatch(
    backend_env$build_transformed_parameter_table(
      path = parameter_table,
      fit_treatment = fit_treatment,
      fit_tau_O2 = fit_tau_O2,
      O2_growth = O2_growth
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
      parameter_table <- default_invivo_parameter_table_path(must_exist = TRUE)
    }
    validate_parameter_table_for_mode(
      backend_env = backend_env,
      parameter_table = parameter_table,
      mode_label = "fit_invivo",
      fit_treatment = isTRUE(as_bool(cfg_raw$fit_treatment, FALSE)),
      fit_tau_O2 = isTRUE(as_bool(cfg_raw$fit_tau_O2, FALSE)),
      O2_growth = isTRUE(as_bool(cfg_raw$O2_growth, TRUE))
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
    parameter_table <- default_invivo_parameter_table_path(must_exist = TRUE)
  }
  validate_parameter_table_for_mode(
    backend_env = backend_env,
    parameter_table = parameter_table,
    mode_label = "fit_invivo",
    fit_treatment = isTRUE(as_bool(argv$fit_treatment, FALSE)),
    fit_tau_O2 = isTRUE(as_bool(argv$fit_tau_O2, FALSE)),
    O2_growth = isTRUE(as_bool(argv$O2_growth, TRUE))
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
  backend_env$ivt_reject_removed_passage_mode(
    argv$passage_mode,
    source = "the in vitro CLI"
  )
  parameter_table <- trim_cli_scalar(argv$parameter_table)
  if (is.null(parameter_table)) {
    parameter_table <- backend_env$default_parameter_table_path(
      must_exist = TRUE
    )
  }
  tryCatch(
    backend_env$validate_invitro_parameter_table(
      parameter_table = parameter_table,
      dt = as_num(argv$dt, 0.05),
      init_total_size = as_num(argv$init_total_size, 1e6),
      o2_upper_bound = as_num(argv$o2_upper_bound, 21),
      fixed_oxygen = TRUE
    ),
    error = function(e) {
      fail_fit_input("fit_invitro", conditionMessage(e))
    }
  )

  fit_objects_dir <- trim_cli_scalar(argv$fit_objects_dir)
  if (is.null(fit_objects_dir)) {
    fit_objects_dir <- backend_env$default_fit_objects_dir(must_exist = TRUE)
  }
  flow_density_path <- trim_cli_scalar(argv$flow_density_path)
  flow_density_use <- tryCatch(
    backend_env$resolve_optional_flow_density_path(flow_density_path),
    error = function(e) {
      fail_fit_input("fit_invitro", conditionMessage(e))
    }
  )
  tryCatch(
    backend_env$validate_invitro_fit_objects(
      fit_objects_dir = fit_objects_dir,
      flow_density_path = flow_density_use
    ),
    error = function(e) {
      fail_fit_input("fit_invitro", conditionMessage(e))
    }
  )
  invisible(list(
    parameter_table = normalizePath(parameter_table, mustWork = FALSE),
    fit_objects_dir = normalizePath(fit_objects_dir, mustWork = FALSE),
    flow_density_path = normalizePath(flow_density_use, mustWork = FALSE)
  ))
}

validate_dispatch_argv <- function(argv, fit_mode, backend_env) {
  if (identical(fit_mode, "fit_invivo")) {
    return(validate_fit_invivo_inputs(argv, backend_env))
  }
  if (identical(fit_mode, "fit_invitro")) {
    return(validate_fit_invitro_inputs(argv, backend_env))
  }
  if (!exists("validate_fit_joint_inputs", envir = backend_env, inherits = FALSE) ||
      !is.function(backend_env$validate_fit_joint_inputs)) {
    stop("Backend validator 'validate_fit_joint_inputs' not found for fit_joint.", call. = FALSE)
  }
  backend_env$validate_fit_joint_inputs(argv)
}

dispatch_backend_main <- function(argv, fit_mode, backend_env) {
  dispatch_argv <- strip_fit_mode_flags(argv)
  if (identical(fit_mode, "fit_invitro")) {
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
