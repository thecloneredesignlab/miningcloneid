#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEoptim))
suppressPackageStartupMessages(library(dplyr))

.o2g_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2g_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = FALSE)
HELPER_DIR <- normalizePath(file.path(OXYGEN_ROOT, "code", "in-vitro-utils"), mustWork = FALSE)

source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "model", "model_O2G_supply_demand_MAP.R"), local = environment())
helper_files <- c("io.R", "lineage_adapter.R", "runner.R", "summaries.R", "objective.R")
for (helper_name in helper_files) {
  helper_path <- file.path(HELPER_DIR, helper_name)
  if (!file.exists(helper_path)) {
    stop("Missing in vitro helper file: ", helper_path)
  }
  sys.source(helper_path, envir = environment(), chdir = TRUE)
}
rm(.o2g_bootstrap_script_dir)

parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
.first_non_null_local <- o2sd_first_non_null

default_parameter_table_path <- function(script_dir = SCRIPT_DIR,
                                         must_exist = FALSE) {
  path <- ivt_parameter_table_path(
    repo_root = OXYGEN_ROOT
  )
  if (isTRUE(must_exist) && !file.exists(path)) {
    stop("Default in vitro parameter table not found: ", path)
  }
  normalizePath(path, mustWork = FALSE)
}

default_fit_objects_dir <- function(script_dir = SCRIPT_DIR, must_exist = FALSE) {
  path <- normalizePath(
    file.path(OXYGEN_ROOT, "ploidyOxygen", "data", "fit_objects"),
    mustWork = FALSE
  )
  if (isTRUE(must_exist) && !dir.exists(path)) {
    stop("Default in vitro fit_objects directory not found: ", path)
  }
  path
}

default_flow_density_path <- function(script_dir = SCRIPT_DIR) {
  normalizePath(
    file.path(OXYGEN_ROOT, "data", "g0g1_ploidy_density_grid.csv"),
    mustWork = FALSE
  )
}

default_out_dir <- function(script_dir = SCRIPT_DIR) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  normalizePath(
    file.path(OXYGEN_ROOT, "results", paste0("fit_model_O2G_supply_demand_MAP_invitro_", stamp)),
    mustWork = FALSE
  )
}

normalize_invitro_n_cores <- function(x) {
  n <- suppressWarnings(as.integer(x))
  if (!is.finite(n) || is.na(n) || n < 1L) 1L else n
}

start_invitro_deoptim_cluster <- function(n_cores) {
  n_use <- normalize_invitro_n_cores(n_cores)
  if (n_use <= 1L) return(NULL)
  if (.Platform$OS.type == "unix" && exists("makeForkCluster", envir = asNamespace("parallel"), inherits = FALSE)) {
    return(parallel::makeForkCluster(n_use))
  }
  parallel::makePSOCKcluster(n_use)
}

resolve_optional_flow_density_path <- function(raw_path = NULL) {
  if (!is.null(raw_path)) {
    path <- normalizePath(raw_path, mustWork = FALSE)
    if (!file.exists(path)) {
      stop("flow_density_path not found: ", path)
    }
    return(path)
  }
  default_path <- default_flow_density_path()
  if (file.exists(default_path)) default_path else default_path
}

ivt_load_fit_objects_compat <- function(fit_objects_dir,
                                        flow_density_path = NULL) {
  load_formals <- names(formals(ivt_load_fit_objects))
  call_args <- list(repo_root = OXYGEN_ROOT)
  if ("fit_objects_dir" %in% load_formals) {
    call_args$fit_objects_dir <- fit_objects_dir
  }
  if ("flow_csv_path" %in% load_formals) {
    call_args$flow_csv_path <- flow_density_path
  }
  do.call(ivt_load_fit_objects, call_args)
}

build_invitro_cfg <- function(parameter_table,
                              dt = 0.05,
                              init_total_size = 1e6,
                              o2_upper_bound = 21,
                              fixed_oxygen = TRUE) {
  cfg <- ivt_build_default_cfg(
    repo_root = OXYGEN_ROOT,
    dt = dt,
    init_total_size = init_total_size,
    o2_upper_bound = o2_upper_bound,
    fixed_oxygen = fixed_oxygen
  )
  cfg$parameter_table <- normalizePath(parameter_table, mustWork = FALSE)
  cfg <- normalize_sim_cfg_common(cfg, context = "fit")
  cfg
}

validate_invitro_parameter_table <- function(parameter_table,
                                             dt = 0.05,
                                             init_total_size = 1e6,
                                             o2_upper_bound = 21,
                                             fixed_oxygen = TRUE) {
  if (!file.exists(parameter_table)) {
    stop("In vitro parameter table not found: ", parameter_table)
  }
  cfg <- build_invitro_cfg(
    parameter_table = parameter_table,
    dt = dt,
    init_total_size = init_total_size,
    o2_upper_bound = o2_upper_bound,
    fixed_oxygen = fixed_oxygen
  )
  ivt_optimizer_spec(cfg)
  invisible(cfg)
}

validate_invitro_fit_objects <- function(fit_objects_dir,
                                         flow_density_path = NULL) {
  ivt_load_fit_objects_compat(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path
  )
  invisible(TRUE)
}

make_penalty_components <- function(objective = 1e9, reason = "penalty") {
  empty_summary <- data.frame()
  empty_result <- list(
    adapter = NULL,
    model_core = NULL,
    grid_pre = integer(0),
    segment_results = list()
  )
  list(
    objective = as.numeric(objective),
    total_loglik = -as.numeric(objective),
    growth_loglik = -as.numeric(objective),
    ploidy_loglik = 0.0,
    flow_loglik = 0.0,
    growth_loglik_sum = -as.numeric(objective),
    ploidy_loglik_sum = 0.0,
    flow_loglik_sum = 0.0,
    sigma_growth = NA_real_,
    sigma_kary = NA_real_,
    sigma_flow_ploidy = NA_real_,
    n_growth = 0L,
    n_growth_observed = 0L,
    n_growth_missing_pred = 0L,
    n_growth_negative_pred = 0L,
    n_ploidy_passages = 0L,
    n_kary_cells = 0L,
    n_flow_passages = 0L,
    n_flow_samples = 0L,
    summary = empty_summary,
    growth_df = data.frame(),
    ploidy_df = data.frame(),
    flow_df = data.frame(),
    flow_overlay_df = data.frame(),
    run_2N = empty_result,
    run_4N = empty_result,
    penalty_reason = as.character(reason)
  )
}

write_tsv_if_nonempty <- function(df, path) {
  if (is.data.frame(df) && nrow(df) > 0L) {
    write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

run_invitro_rscript_to_log <- function(label, script_path, args, log_path) {
  if (!file.exists(script_path)) {
    stop("Missing ", label, " script: ", script_path, call. = FALSE)
  }
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript",
    args = c(normalizePath(script_path, mustWork = TRUE), args),
    stdout = log_path,
    stderr = log_path
  )
  if (!identical(status, 0L)) {
    tail_txt <- if (file.exists(log_path)) {
      paste(utils::tail(readLines(log_path, warn = FALSE), 25L), collapse = "\n")
    } else {
      ""
    }
    stop(
      label,
      " failed with exit status ",
      status,
      ". See ",
      log_path,
      if (nzchar(tail_txt)) paste0("\nLast log lines:\n", tail_txt) else "",
      call. = FALSE
    )
  }
  invisible(status)
}

run_invitro_auto_viz_report <- function(out_dir) {
  viz_script <- file.path(WORKFLOW_ROOT, "vis", "viz_invitro_model_O2G_supply_demand_MAP_results.R")
  report_script <- file.path(WORKFLOW_ROOT, "report", "render_invitro_fit_report.R")
  run_invitro_rscript_to_log(
    label = "invitro viz",
    script_path = viz_script,
    args = paste0("--fit_dir=", normalizePath(out_dir, mustWork = FALSE)),
    log_path = file.path(out_dir, "viz_status.log")
  )
  run_invitro_rscript_to_log(
    label = "invitro report",
    script_path = report_script,
    args = paste0("--fit_dir=", normalizePath(out_dir, mustWork = FALSE)),
    log_path = file.path(out_dir, "report_status.log")
  )
  invisible(TRUE)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  parameter_table <- if (!is.null(argv$parameter_table)) {
    argv$parameter_table
  } else {
    default_parameter_table_path(must_exist = TRUE)
  }
  fit_objects_dir <- if (!is.null(argv$fit_objects_dir)) {
    argv$fit_objects_dir
  } else {
    default_fit_objects_dir(must_exist = TRUE)
  }
  flow_density_path <- resolve_optional_flow_density_path(argv$flow_density_path)
  out_dir <- if (!is.null(argv$out_dir)) argv$out_dir else default_out_dir()

  seed <- as.integer(.first_non_null_local(argv$seed, 1L))
  itermax_requested <- as.integer(.first_non_null_local(argv$itermax, 500L))
  itermax_max <- as.integer(.first_non_null_local(argv$itermax_max, 500L))
  if (!is.finite(itermax_requested) || is.na(itermax_requested)) itermax_requested <- 500L
  if (!is.finite(itermax_max) || is.na(itermax_max) || itermax_max < 1L) itermax_max <- 500L
  itermax <- min(max(itermax_requested, 1L), itermax_max)
  NP_requested <- as.integer(.first_non_null_local(argv$NP, 80L))
  n_cores_requested <- normalize_invitro_n_cores(.first_non_null_local(argv$n_cores, 1L))
  de_reltol <- as.numeric(.first_non_null_local(argv$de_reltol, 1e-4))
  de_steptol <- as.integer(.first_non_null_local(argv$de_steptol, 25L))
  if (!is.finite(de_reltol) || de_reltol <= 0) de_reltol <- 1e-4
  if (!is.finite(de_steptol) || is.na(de_steptol) || de_steptol < 1L) de_steptol <- 25L
  dt_use <- as.numeric(.first_non_null_local(argv$dt, 0.05))
  init_total_size_use <- as.numeric(.first_non_null_local(argv$init_total_size, 1e6))
  o2_upper_bound_use <- as.numeric(.first_non_null_local(argv$o2_upper_bound, 21))
  fixed_oxygen_use <- TRUE
  auto_viz <- as_bool(.first_non_null_local(argv$auto_viz, TRUE), TRUE)

  validate_invitro_parameter_table(
    parameter_table = parameter_table,
    dt = dt_use,
    init_total_size = init_total_size_use,
    o2_upper_bound = o2_upper_bound_use,
    fixed_oxygen = fixed_oxygen_use
  )
  validate_invitro_fit_objects(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(parameter_table, file.path(out_dir, "parameter_table_input.csv"), overwrite = TRUE)
  set.seed(seed)

  cfg_local <- build_invitro_cfg(
    parameter_table = parameter_table,
    dt = dt_use,
    init_total_size = init_total_size_use,
    o2_upper_bound = o2_upper_bound_use,
    fixed_oxygen = fixed_oxygen_use
  )
  fit_objects <- ivt_load_fit_objects_compat(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path
  )

  optim_spec <- ivt_optimizer_spec(cfg_local)
  free_names <- optim_spec$param_name
  init_free <- setNames(as.numeric(optim_spec$init), free_names)
  lower_free <- setNames(as.numeric(optim_spec$lower), free_names)
  upper_free <- setNames(as.numeric(optim_spec$upper), free_names)

  objective_from_free <- function(par_free_t) {
    par_t <- init_free
    par_t[free_names] <- as.numeric(par_free_t)
    run_params <- ivt_optim_par_to_run_params(par_t = par_t, cfg = cfg_local)
    comp <- tryCatch(
      ivt_objective_components(
        run_params = run_params,
        fit_objects = fit_objects,
        cfg = cfg_local,
        fallback_max_passage_days = 14,
        growth_weight = 1,
        ploidy_weight = 1,
        flow_weight = 1
      ),
      error = function(e) {
        make_penalty_components(reason = paste0("simulation_error: ", conditionMessage(e)))
      }
    )
    comp$run_params <- run_params
    comp$full_t <- par_t
    comp
  }

  objective_value <- function(par_free_t) {
    objective_from_free(par_free_t)$objective
  }

  NP_use <- max(NP_requested, 10L * length(free_names))
  de_ctrl <- list(
    trace = TRUE,
    itermax = max(itermax, 1L),
    NP = NP_use,
    reltol = de_reltol,
    steptol = de_steptol
  )
  de_cluster <- NULL
  de_active_cores <- 1L
  if (n_cores_requested > 1L) {
    message("[fit_invitro] DEoptim parallel requested with n_cores=", n_cores_requested, ".")
    de_cluster <- tryCatch(
      start_invitro_deoptim_cluster(n_cores_requested),
      error = function(e) {
        stop("[fit_invitro] Could not start DEoptim workers: ", conditionMessage(e), call. = FALSE)
      }
    )
    on.exit(try(parallel::stopCluster(de_cluster), silent = TRUE), add = TRUE)
    parallel::clusterExport(
      de_cluster,
      varlist = c("objective_value"),
      envir = environment()
    )
    de_ctrl$cluster <- de_cluster
    de_active_cores <- length(de_cluster)
    message("[fit_invitro] DEoptim parallel enabled: workers=", de_active_cores, ".")
  } else {
    de_ctrl$parallelType <- "none"
    message("[fit_invitro] DEoptim running in serial mode (n_cores=1).")
  }
  de_fit <- DEoptim::DEoptim(
    fn = objective_value,
    lower = lower_free,
    upper = upper_free,
    control = de_ctrl
  )
  de_best_free_t <- as.numeric(de_fit$optim$bestmem)
  names(de_best_free_t) <- free_names
  de_best_objective <- suppressWarnings(as.numeric(de_fit$optim$bestval))
  best_free_t <- de_best_free_t

  local_maxit <- as_int(.first_non_null_local(argv$local_optim_maxit, argv$optim_maxit, 200L), 200L)
  if (!is.finite(local_maxit) || is.na(local_maxit) || local_maxit < 1L) local_maxit <- 200L
  local_attempted <- FALSE
  local_accepted <- FALSE
  local_fit <- NULL
  local_best_objective <- NA_real_
  local_convergence <- NA_integer_
  local_message <- NA_character_
  if (is.finite(de_best_objective)) {
    local_attempted <- TRUE
    message("[fit_invitro] Starting L-BFGS-B local refinement from DEoptim best; maxit=", local_maxit, ".")
    local_fit <- tryCatch(
      suppressWarnings(
        optim(
          par = best_free_t,
          fn = objective_value,
          method = "L-BFGS-B",
          lower = lower_free,
          upper = upper_free,
          control = list(maxit = local_maxit)
        )
      ),
      error = function(e) {
        warning("[fit_invitro] L-BFGS-B local refinement failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.list(local_fit)) {
      local_best_objective <- suppressWarnings(as.numeric(local_fit$value))
      local_convergence <- suppressWarnings(as.integer(local_fit$convergence))
      local_message <- as.character(.first_non_null_local(local_fit$message, NA_character_))
      if (is.finite(local_best_objective) && local_best_objective < de_best_objective) {
        best_free_t <- as.numeric(local_fit$par)
        names(best_free_t) <- free_names
        best_free_t <- o2sd_clip(best_free_t, lower_free, upper_free)
        local_accepted <- TRUE
        message(
          "[fit_invitro] L-BFGS-B local refinement improved objective: ",
          signif(de_best_objective, 8), " -> ", signif(local_best_objective, 8), "."
        )
      } else {
        message("[fit_invitro] L-BFGS-B local refinement did not improve objective; keeping DEoptim best.")
      }
    }
  }
  de_method <- if (de_active_cores > 1L) "DEoptim_parallel" else "DEoptim_serial"
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

  best_comp <- objective_from_free(best_free_t)
  best_run_params <- best_comp$run_params
  best_full_t <- best_comp$full_t

  best_numeric_params <- best_run_params[vapply(best_run_params, is.numeric, logical(1))]
  best_numeric_params <- filter_family_specific_run_params_for_output_common(best_numeric_params)
  best_numeric_params <- best_numeric_params[!vapply(best_numeric_params, is.null, logical(1))]

  best_params_df <- data.frame(
    parameter = names(best_numeric_params),
    value = as.numeric(unlist(best_numeric_params)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(best_params_df, file = file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  best_transformed_df <- data.frame(
    transformed_parameter = names(best_full_t),
    transformed_value = as.numeric(best_full_t),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(best_transformed_df, file = file.path(out_dir, "best_params_transformed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  write_tsv_if_nonempty(best_comp$summary, file.path(out_dir, "invitro_lineage_summary.tsv"))
  write_tsv_if_nonempty(best_comp$growth_df, file.path(out_dir, "invitro_growth_loglik.tsv"))
  write_tsv_if_nonempty(best_comp$ploidy_df, file.path(out_dir, "invitro_ploidy_loglik.tsv"))
  write_tsv_if_nonempty(best_comp$flow_df, file.path(out_dir, "invitro_flow_loglik.tsv"))
  write_tsv_if_nonempty(best_comp$flow_overlay_df, file.path(out_dir, "invitro_flow_overlay.tsv"))

  dist_summary <- dplyr::bind_rows(
    ivt_collect_distribution_summary(best_comp$run_2N),
    ivt_collect_distribution_summary(best_comp$run_4N)
  )
  write_tsv_if_nonempty(dist_summary, file.path(out_dir, "invitro_distribution_summary.tsv"))

  ploidy_quantile_probs <- seq(0.01, 0.99, length.out = 50L)
  dist_quantiles <- dplyr::bind_rows(
    ivt_collect_distribution_quantiles(best_comp$run_2N, probs = ploidy_quantile_probs),
    ivt_collect_distribution_quantiles(best_comp$run_4N, probs = ploidy_quantile_probs)
  )
  write_tsv_if_nonempty(dist_quantiles, file.path(out_dir, "invitro_distribution_quantiles.tsv"))

  daily_counts <- dplyr::bind_rows(
    ivt_collect_daily_counts(best_comp$run_2N),
    ivt_collect_daily_counts(best_comp$run_4N)
  )
  write_tsv_if_nonempty(daily_counts, file.path(out_dir, "invitro_daily_counts.tsv"))

  observed_kary <- dplyr::bind_rows(
    ivt_collect_observed_kary_summary(best_comp$run_2N, fit_objects$fit_data),
    ivt_collect_observed_kary_summary(best_comp$run_4N, fit_objects$fit_data)
  )
  write_tsv_if_nonempty(observed_kary, file.path(out_dir, "invitro_observed_kary.tsv"))

  observed_flow <- dplyr::bind_rows(
    ivt_collect_observed_flow_summary(best_comp$run_2N, fit_objects$fit_data),
    ivt_collect_observed_flow_summary(best_comp$run_4N, fit_objects$fit_data)
  )
  write_tsv_if_nonempty(observed_flow, file.path(out_dir, "invitro_observed_flow.tsv"))

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
      "objective_total",
      "total_loglik",
      "growth_loglik",
      "ploidy_loglik",
      "flow_loglik",
      "growth_loglik_sum",
      "ploidy_loglik_sum",
      "flow_loglik_sum",
      "sigma_growth",
      "sigma_kary",
      "sigma_flow_ploidy",
      "n_growth",
      "n_growth_observed",
      "n_growth_missing_pred",
      "n_growth_negative_pred",
      "n_ploidy_passages",
      "n_kary_cells",
      "n_flow_passages",
      "n_flow_samples",
      "seed",
      "itermax",
      "itermax_requested",
      "itermax_max",
      "de_reltol",
      "de_steptol",
      "NP_requested",
      "NP_used",
      "n_cores_requested",
      "n_cores_used",
      "dt",
      "init_total_size",
      "parameter_table",
      "fit_objects_dir",
      "flow_density_path"
    ),
    value = c(
      "fit_invitro",
      optimizer_method,
      as.character(de_best_objective),
      as.character(local_best_objective),
      as.character(local_attempted),
      as.character(local_accepted),
      as.character(local_convergence),
      as.character(local_maxit),
      as.character(best_comp$objective),
      as.character(best_comp$total_loglik),
      as.character(best_comp$growth_loglik),
      as.character(best_comp$ploidy_loglik),
      as.character(best_comp$flow_loglik),
      as.character(best_comp$growth_loglik_sum),
      as.character(best_comp$ploidy_loglik_sum),
      as.character(best_comp$flow_loglik_sum),
      as.character(best_comp$sigma_growth),
      as.character(best_comp$sigma_kary),
      as.character(best_comp$sigma_flow_ploidy),
      as.character(best_comp$n_growth),
      as.character(best_comp$n_growth_observed),
      as.character(best_comp$n_growth_missing_pred),
      as.character(best_comp$n_growth_negative_pred),
      as.character(best_comp$n_ploidy_passages),
      as.character(best_comp$n_kary_cells),
      as.character(best_comp$n_flow_passages),
      as.character(best_comp$n_flow_samples),
      as.character(seed),
      as.character(itermax),
      as.character(itermax_requested),
      as.character(itermax_max),
      as.character(de_reltol),
      as.character(de_steptol),
      as.character(NP_requested),
      as.character(NP_use),
      as.character(n_cores_requested),
      as.character(de_active_cores),
      as.character(dt_use),
      as.character(init_total_size_use),
      normalizePath(parameter_table, mustWork = FALSE),
      normalizePath(fit_objects_dir, mustWork = FALSE),
      normalizePath(flow_density_path, mustWork = FALSE)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  summary_df <- filter_fit_summary_metrics_for_output_common(summary_df)
  write.table(summary_df, file = file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  saveRDS(
    list(
      cfg = cfg_local,
      deoptim = de_fit,
      local_optim = local_fit,
      optimizer_trace = optimizer_trace,
      best_components = best_comp,
      best_params = best_run_params,
      fit_objects_dir = normalizePath(fit_objects_dir, mustWork = FALSE),
      flow_density_path = normalizePath(flow_density_path, mustWork = FALSE)
    ),
    file = file.path(out_dir, "fit_result.rds")
  )

  if (isTRUE(auto_viz)) {
    run_invitro_auto_viz_report(out_dir)
  }

  message("Done. Results written to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(normalizePath(out_dir, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  main()
}
