#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

si3_extra_results_environment <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    extra_results <- file.path(
      MODEL_CODE_ROOT, "analysis", "fit_results", "extra_results.R"
    )
    env <- new.env(parent = globalenv())
    env$SCRIPT_DIR <- dirname(extra_results)
    env$WORKFLOW_ROOT <- MODEL_CODE_ROOT
    sys.source(
      file.path(
        MODEL_CODE_ROOT, "util", "o2_supply_demand_map_shared.R"
      ),
      envir = env
    )
    sys.source(
      file.path(
        MODEL_CODE_ROOT, "util",
        "o2_supply_demand_map_common_semantics.R"
      ),
      envir = env
    )
    env$`%||%` <- env$o2sd_null_coalesce
    env$parse_args <- env$o2sd_parse_args
    env$as_num <- env$o2sd_as_num
    env$as_bool <- env$o2sd_as_bool

    lines <- readLines(extra_results, warn = FALSE)
    start <- grep(
      "^rm\\(\\.o2sd_bootstrap_script_dir\\)", lines
    )[[1L]] + 1L
    stop_line <- grep(
      "^if \\(sys\\.nframe\\(\\) == 0\\)", lines
    )[[1L]] - 1L
    eval(
      parse(text = lines[start:stop_line], keep.source = FALSE),
      envir = env
    )
    cache <<- env
    cache
  }
})

si3_collect_seed_metrics <- function(run_dir, objective_metric) {
  metric_value <- function(metrics, name, default = NA_character_) {
    if (name %in% names(metrics)) metrics[[name]] else default
  }
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
        metric_value(metrics, "fit_status", "ok")
      ),
      optimizer_local_accepted =
        as_boolean(
          metric_value(metrics, "optimizer_local_accepted"), FALSE
        ),
      optimizer_local_convergence =
        suppressWarnings(as.numeric(metric_value(
          metrics, "optimizer_local_convergence"
        ))),
      deoptim_stop_reason = as.character(
        metric_value(metrics, "deoptim_stop_reason")
      ),
      source_file = normalizePath(summary_path, mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  out <- do.call(rbind, rows)
  out <- out[order(
    out$objective,
    as.integer(sub("^seed", "", out$seed))
  ), , drop = FALSE]
  out$objective_rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  if (nrow(out) != 500L) {
    stop("Expected 500 objective-bearing seeds under: ", run_dir)
  }
  out
}

si3_selected_boundary_summary <- function(run_dir, seed) {
  extra <- si3_extra_results_environment()
  seed_dir <- file.path(run_dir, seed)
  fit_summary_vals <- extra$read_metric_map(
    file.path(seed_dir, "fit_summary.tsv"), "metric", "value"
  )
  fit_summary_vals <- extra$supplement_joint_invitro_metrics(
    fit_summary_vals, seed_dir
  )
  best_vals <- extra$read_metric_map(
    file.path(seed_dir, "best_params.tsv"), "parameter", "value"
  )
  parameter_path <- extra$find_parameter_table_for_seed(
    run_dir, seed_dir, fit_summary_vals
  )
  parameter_table <- extra$read_parameter_table_checked(parameter_path)
  objective <- extra$as_num(extra$summary_metric_value(
    fit_summary_vals,
    "objective",
    extra$summary_metric_value(
      fit_summary_vals, "objective_total", NA_real_
    )
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
  data.frame(
    seed = seed,
    optimizer_local_accepted =
      as.logical(record$optimizer_local_accepted),
    optimizer_local_convergence =
      as.numeric(record$optimizer_local_convergence),
    n_at_bound_active = as.integer(record$n_at_bound_active),
    n_active_params = as.integer(record$n_active_params),
    parameter_table_source = normalizePath(parameter_path, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

si3_write_bundle <- function(
    metrics,
    output_dir,
    fit_label,
    selected_seed,
    run_dir
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  objective <- metrics[, c("seed", "objective_rank", "objective"), drop = FALSE]
  write_intermediate_tsv(
    objective, file.path(output_dir, "seed_objective_simple.tsv")
  )
  selected <- si3_selected_boundary_summary(
    run_dir = run_dir, seed = selected_seed
  )
  boundary_sources <- c(
    file.path(run_dir, selected_seed, "best_params.tsv"),
    selected$parameter_table_source[[1L]]
  )
  selected$parameter_table_source <- NULL
  write_intermediate_tsv(
    selected, file.path(output_dir, "seed_summary.tsv")
  )
  deoptim_converged <- sum(
    !is.na(metrics$deoptim_stop_reason) &
      nzchar(metrics$deoptim_stop_reason) &
      metrics$deoptim_stop_reason != "NA"
  )
  convergence <- data.frame(
    Fit = fit_label,
    `Total seeds` = nrow(metrics),
    `DEoptim converged` = deoptim_converged,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (identical(fit_label, "joint fitting")) {
    convergence$`L-BFGS-B accepted` <-
      sum(metrics$optimizer_local_accepted, na.rm = TRUE)
    convergence$`Converged and accepted` <- sum(
      metrics$optimizer_local_accepted &
        is.finite(metrics$optimizer_local_convergence) &
        metrics$optimizer_local_convergence == 0,
      na.rm = TRUE
    )
    convergence$`Converged only` <-
      convergence$`DEoptim converged` -
      convergence$`Converged and accepted`
    convergence$`Accepted only` <-
      convergence$`L-BFGS-B accepted` -
      convergence$`Converged and accepted`
  }
  write_intermediate_tsv(
    convergence, file.path(output_dir, "convergence_summary.tsv")
  )
  attr(selected, "boundary_sources") <- boundary_sources
  invisible(selected)
}

data_SI_Figure3 <- function() {
  output_root <- file.path(DATA_ROOT, "SI_Figure3")
  selection_path <- file.path(
    DATA_ROOT, "Figure5", "figure5_frozen_inputs",
    "selected_results.tsv"
  )
  require_files(selection_path, "Figure 5 selected-results table")
  selection <- utils::read.delim(
    selection_path, check.names = FALSE, stringsAsFactors = FALSE
  )
  joint_selection <- selection[
    selection$record_type == "joint_pair_best", , drop = FALSE
  ]
  if (nrow(joint_selection) != 6L) {
    stop("SI Figure 3 requires six selected joint winners.")
  }

  invitro_metrics <- si3_collect_seed_metrics(
    INVITRO_RESULT_ROOT, "objective_total"
  )
  invivo_metrics <- si3_collect_seed_metrics(
    INVIVO_RESULT_ROOT, "objective"
  )
  invitro_selected <- si3_write_bundle(
    invitro_metrics,
    file.path(output_root, "separate", "invitro"),
    "in vitro", "seed10", INVITRO_RESULT_ROOT
  )
  invivo_selected <- si3_write_bundle(
    invivo_metrics,
    file.path(output_root, "separate", "invivo"),
    "in vivo", "seed25", INVIVO_RESULT_ROOT
  )

  joint_source_files <- character()
  boundary_source_files <- c(
    attr(invitro_selected, "boundary_sources"),
    attr(invivo_selected, "boundary_sources")
  )
  for (i in seq_len(nrow(joint_selection))) {
    row <- joint_selection[i, , drop = FALSE]
    pair_id <- paste0("fit_joint_", row$warmup_label[[1L]])
    pair_dir <- file.path(JOINT_RESULT_ROOT, pair_id)
    metrics <- si3_collect_seed_metrics(pair_dir, "objective")
    selected_summary <- si3_write_bundle(
      metrics,
      file.path(output_root, "joint", row$warmup_label[[1L]]),
      "joint fitting",
      row$selected_seed[[1L]],
      pair_dir
    )
    joint_source_files <- c(joint_source_files, metrics$source_file)
    boundary_source_files <- c(
      boundary_source_files,
      attr(selected_summary, "boundary_sources")
    )
  }

  sources <- c(
    invitro_metrics$source_file,
    invivo_metrics$source_file,
    joint_source_files,
    boundary_source_files
  )
  contract <- data.frame(
    role = c(
      rep(
        "optimizer diagnostic fit summary",
        length(invitro_metrics$source_file) +
          length(invivo_metrics$source_file) +
          length(joint_source_files)
      ),
      rep(
        "selected-fit boundary diagnostic input",
        length(boundary_source_files)
      )
    ),
    source = sources,
    local_file = NA_character_,
    source_md5 = unname(tools::md5sum(sources)),
    local_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  write_data_contract("SI_Figure3", contract)
  invisible(TRUE)
}

if (sys.nframe() == 0L) data_SI_Figure3()
