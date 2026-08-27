optimizer_extra_results_environment <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    extra_results <- file.path(
      MODEL_CODE_ROOT, "analysis", "fit_results", "extra_results.R"
    )
    required <- c(
      extra_results,
      file.path(MODEL_CODE_ROOT, "util", "o2_supply_demand_map_shared.R"),
      file.path(
        MODEL_CODE_ROOT, "util", "o2_supply_demand_map_common_semantics.R"
      )
    )
    require_files(required, "optimizer diagnostic model helper")

    env <- new.env(parent = globalenv())
    env$SCRIPT_DIR <- dirname(extra_results)
    env$WORKFLOW_ROOT <- MODEL_CODE_ROOT
    sys.source(required[[2L]], envir = env)
    sys.source(required[[3L]], envir = env)
    env$`%||%` <- env$o2sd_null_coalesce
    env$parse_args <- env$o2sd_parse_args
    env$as_num <- env$o2sd_as_num
    env$as_bool <- env$o2sd_as_bool

    lines <- readLines(extra_results, warn = FALSE)
    start <- grep("^rm\\(\\.o2sd_bootstrap_script_dir\\)", lines)[[1L]] + 1L
    stop_line <- grep("^if \\(sys\\.nframe\\(\\) == 0\\)", lines)[[1L]] - 1L
    eval(
      parse(text = lines[start:stop_line], keep.source = FALSE),
      envir = env
    )
    cache <<- env
    cache
  }
})

optimizer_metric_value <- function(metrics, name, default = NA_character_) {
  if (name %in% names(metrics)) metrics[[name]] else default
}

optimizer_collect_seed_metrics <- function(
    run_dir,
    objective_metric,
    expected_seed_count = 500L
) {
  run_dir <- assert_allowed_input(run_dir)
  seed_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
  rows <- lapply(seed_dirs, function(seed_dir) {
    summary_path <- file.path(seed_dir, "fit_summary.tsv")
    if (!file.exists(summary_path)) return(NULL)
    metrics <- read_metric_map(summary_path)
    objective <- suppressWarnings(as.numeric(metrics[[objective_metric]]))
    if (!is.finite(objective)) return(NULL)
    data.frame(
      seed = basename(seed_dir),
      objective = objective,
      fit_status = as.character(
        optimizer_metric_value(metrics, "fit_status", "ok")
      ),
      optimizer_local_attempted = as_boolean(
        optimizer_metric_value(metrics, "optimizer_local_attempted"), FALSE
      ),
      optimizer_local_accepted = as_boolean(
        optimizer_metric_value(metrics, "optimizer_local_accepted"), FALSE
      ),
      optimizer_local_convergence = suppressWarnings(as.numeric(
        optimizer_metric_value(metrics, "optimizer_local_convergence")
      )),
      optimizer_deoptim_objective = suppressWarnings(as.numeric(
        optimizer_metric_value(metrics, "optimizer_deoptim_objective")
      )),
      optimizer_local_objective = suppressWarnings(as.numeric(
        optimizer_metric_value(metrics, "optimizer_local_objective")
      )),
      optimizer_iter_completed = suppressWarnings(as.numeric(
        optimizer_metric_value(
          metrics,
          "optimizer_iter_completed",
          optimizer_metric_value(metrics, "deoptim_stop_iteration")
        )
      )),
      optimizer_iter_target = suppressWarnings(as.numeric(
        optimizer_metric_value(
          metrics,
          "optimizer_iter_target",
          optimizer_metric_value(metrics, "deoptim_iter_target")
        )
      )),
      deoptim_stop_reason = as.character(
        optimizer_metric_value(metrics, "deoptim_stop_reason")
      ),
      source_file = normalizePath(summary_path, mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) stop("No objective-bearing seeds under: ", run_dir)
  out <- do.call(rbind, rows)
  out <- out[order(
    out$objective,
    as.integer(sub("^seed", "", out$seed))
  ), , drop = FALSE]
  out$objective_rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  if (!is.null(expected_seed_count) && nrow(out) != expected_seed_count) {
    stop(
      "Expected ", expected_seed_count,
      " objective-bearing seeds under: ", run_dir,
      "; observed ", nrow(out)
    )
  }
  out
}

optimizer_boundary_summary_for_seed <- function(run_dir, seed) {
  extra <- optimizer_extra_results_environment()
  seed_dir <- file.path(run_dir, seed)
  required <- c(
    file.path(seed_dir, "fit_summary.tsv"),
    file.path(seed_dir, "best_params.tsv")
  )
  require_files(required, "optimizer boundary input")
  fit_summary_vals <- extra$read_metric_map(
    required[[1L]], "metric", "value"
  )
  fit_summary_vals <- extra$supplement_joint_invitro_metrics(
    fit_summary_vals, seed_dir
  )
  best_vals <- extra$read_metric_map(
    required[[2L]], "parameter", "value"
  )
  parameter_path <- extra$find_parameter_table_for_seed(
    run_dir, seed_dir, fit_summary_vals
  )
  parameter_table <- extra$read_parameter_table_checked(parameter_path)
  objective <- extra$as_num(extra$summary_metric_value(
    fit_summary_vals,
    "objective",
    extra$summary_metric_value(fit_summary_vals, "objective_total", NA_real_)
  ))
  parameter_long <- extra$build_parameter_long_table(
    seed = seed,
    objective = objective,
    fit_summary_vals = fit_summary_vals,
    best_vals = best_vals,
    param_table = parameter_table,
    near_thresh = 0.05
  )
  record <- extra$build_seed_summary_record(
    seed = seed,
    fit_summary_vals = fit_summary_vals,
    best_vals = best_vals,
    parameter_long = parameter_long,
    pred_gate_metrics = NULL
  )
  out <- data.frame(
    seed = seed,
    optimizer_local_accepted = as.logical(record$optimizer_local_accepted),
    optimizer_local_convergence = as.numeric(
      record$optimizer_local_convergence
    ),
    n_at_bound_active = as.integer(record$n_at_bound_active),
    n_near_bound_only_active = as.integer(
      record$n_near_bound_only_active
    ),
    n_active_params = as.integer(record$n_active_params),
    parameter_table_source = normalizePath(parameter_path, mustWork = TRUE),
    best_params_source = normalizePath(required[[2L]], mustWork = TRUE),
    stringsAsFactors = FALSE
  )
  attr(out, "source_files") <- c(
    required[[1L]], required[[2L]], parameter_path
  )
  out
}

optimizer_collect_boundary_summaries <- function(run_dir, metrics) {
  rows <- lapply(
    metrics$seed,
    function(seed) optimizer_boundary_summary_for_seed(run_dir, seed)
  )
  sources <- unique(unlist(
    lapply(rows, attr, which = "source_files"),
    use.names = FALSE
  ))
  out <- do.call(rbind, rows)
  out$parameter_table_source <- NULL
  out$best_params_source <- NULL
  out <- merge(
    metrics[, c("seed", "objective_rank", "objective"), drop = FALSE],
    out,
    by = "seed",
    all.x = TRUE,
    sort = FALSE
  )
  out <- out[order(out$objective_rank), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "source_files") <- sources
  out
}

optimizer_write_bundle <- function(
    metrics,
    output_dir,
    fit_label,
    selected_seed,
    run_dir,
    include_all_boundaries = FALSE
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  objective <- metrics[, c(
    "seed", "objective_rank", "objective"
  ), drop = FALSE]
  write_intermediate_tsv(
    objective, file.path(output_dir, "seed_objective_simple.tsv")
  )
  diagnostic_columns <- setdiff(names(metrics), "source_file")
  write_intermediate_tsv(
    metrics[, diagnostic_columns, drop = FALSE],
    file.path(output_dir, "seed_optimizer_diagnostics.tsv")
  )

  selected <- optimizer_boundary_summary_for_seed(
    run_dir = run_dir, seed = selected_seed
  )
  selected_sources <- attr(selected, "source_files")
  selected$parameter_table_source <- NULL
  selected$best_params_source <- NULL
  write_intermediate_tsv(
    selected, file.path(output_dir, "seed_summary.tsv")
  )

  deoptim_flagged <- sum(
    !is.na(metrics$deoptim_stop_reason) &
      nzchar(metrics$deoptim_stop_reason) &
      metrics$deoptim_stop_reason != "NA"
  )
  convergence <- data.frame(
    Fit = fit_label,
    `Total seeds` = nrow(metrics),
    `DEoptim convergence flag` = deoptim_flagged,
    `Local attempted` = sum(metrics$optimizer_local_attempted, na.rm = TRUE),
    `Local accepted` = sum(metrics$optimizer_local_accepted, na.rm = TRUE),
    `Local code 0` = sum(
      is.finite(metrics$optimizer_local_convergence) &
        metrics$optimizer_local_convergence == 0,
      na.rm = TRUE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write_intermediate_tsv(
    convergence, file.path(output_dir, "convergence_summary.tsv")
  )

  all_boundary_sources <- character()
  if (isTRUE(include_all_boundaries)) {
    boundaries <- optimizer_collect_boundary_summaries(run_dir, metrics)
    all_boundary_sources <- attr(boundaries, "source_files")
    write_intermediate_tsv(
      boundaries, file.path(output_dir, "seed_boundary_summary.tsv")
    )
  }

  source_files <- unique(c(
    metrics$source_file,
    selected_sources,
    all_boundary_sources
  ))
  attr(selected, "source_files") <- source_files
  invisible(selected)
}
