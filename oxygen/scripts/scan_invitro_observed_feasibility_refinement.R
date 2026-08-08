#!/usr/bin/env Rscript

# Second-round, parameter-only observed-data feasibility scan for the in-vitro
# model.  This driver reuses the first-round simulator/evaluator without
# changing the model, objective, likelihoods, or adding functional penalties.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

.ref_script_dir <- local({
  frame_files <- unlist(lapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value) || !length(value)) character() else as.character(value[[1L]])
  }), use.names = FALSE)
  frame_files <- frame_files[
    basename(frame_files) == "scan_invitro_observed_feasibility_refinement.R"
  ]
  if (length(frame_files)) {
    return(dirname(normalizePath(frame_files[[length(frame_files)]], mustWork = TRUE)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  file_path <- sub("^--file=", "", file_arg)
  file_path <- file_path[
    basename(file_path) == "scan_invitro_observed_feasibility_refinement.R"
  ]
  if (length(file_path)) {
    dirname(normalizePath(file_path[[1L]], mustWork = TRUE))
  } else {
    cursor <- normalizePath(getwd(), mustWork = TRUE)
    repeat {
      candidate <- file.path(
        cursor, "oxygen", "scripts", "scan_invitro_observed_feasibility_refinement.R"
      )
      if (file.exists(candidate)) {
        return(normalizePath(dirname(candidate), mustWork = TRUE))
      }
      parent <- dirname(cursor)
      if (identical(parent, cursor)) break
      cursor <- parent
    }
    getwd()
  }
})

if (!nzchar(Sys.getenv("O2SD_PROJECT_ROOT", unset = ""))) {
  Sys.setenv(O2SD_PROJECT_ROOT = normalizePath(
    file.path(.ref_script_dir, "..", ".."), mustWork = TRUE
  ))
}

# Sourcing the first-round driver supplies the tested fit-object loading,
# objective replay, flow/karyotype metric collection, fork-cluster preflight,
# and common CLI/path helpers.  Its main() is not executed when sourced.
source(
  file.path(.ref_script_dir, "scan_invitro_observed_feasibility.R"),
  local = environment(),
  chdir = TRUE
)

.ivt_ref_base_evaluate_task <- .ivt_scan_evaluate_task
.ivt_ref_parameter_names <- c(
  "p_wgd", "p_mis_base", "gamma_growth", "mu_hp", "gamma_mu",
  "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
  "buffer_n_exp", "O2_crit", "n_O"
)
.ivt_ref_schema_version <- 2L

.ivt_ref_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Failed to atomically install: ", path, call. = FALSE)
  invisible(path)
}

.ivt_ref_digest_object <- function(object) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3)
  digest <- .ivt_scan_sha256(tmp)
  if (!is.character(digest) || length(digest) != 1L || is.na(digest) || !nzchar(digest)) {
    stop("Unable to compute a SHA-256 design digest.", call. = FALSE)
  }
  digest
}

.ivt_ref_add_metrics <- function(result, task) {
  summary <- result$summary
  if (!nrow(summary)) return(result)
  flow <- result$flow
  control_wgd <- NA_real_
  n_low_groups <- 0L
  if (is.data.frame(flow) && nrow(flow)) {
    control_keep <- as.character(flow$segment_id) == "2N-C-A12" &
      as.character(flow$series) == "Predicted"
    values <- suppressWarnings(as.numeric(flow$wgd_mass[control_keep]))
    values <- values[is.finite(values)]
    if (length(values)) control_wgd <- mean(values)
    low_keep <- as.character(flow$series) == "Predicted" &
      grepl("^2N-O[12]-A(12|18|23)$", as.character(flow$segment_id))
    n_low_groups <- length(unique(as.character(flow$segment_id[low_keep])))
  }
  summary$control_flow_wgd_mass <- control_wgd
  summary$n_low_o2_flow_groups <- as.integer(n_low_groups)
  initial_4n <- suppressWarnings(as.numeric(.ivt_scan_context$base_params$init_mean_4N))
  terminal_4n <- suppressWarnings(as.numeric(summary$pooled_4N_A17_mean_N))
  summary$decline_4N_to_A17_N <- if (
    length(initial_4n) == 1L && is.finite(initial_4n)
  ) initial_4n - terminal_4n else NA_real_
  for (passage in c(12L, 18L, 23L)) {
    predicted <- suppressWarnings(as.numeric(
      summary[[paste0("flow_2N_low_A", passage, "_wgd_mass")]]
    ))
    observed <- suppressWarnings(as.numeric(
      summary[[paste0("observed_flow_2N_low_A", passage, "_wgd_mass")]]
    ))
    summary[[paste0("flow_2N_low_A", passage, "_wgd_mass_abs_error")]] <-
      abs(predicted - observed)
  }
  result$summary <- summary
  result
}

# Override only the scan evaluator symbol used by the inherited cluster/task
# machinery.  The underlying objective replay remains the first-round one.
.ivt_scan_evaluate_task <- function(task) {
  .ivt_ref_add_metrics(.ivt_ref_base_evaluate_task(task), task)
}

.ivt_ref_read_anchors <- function(path) {
  if (is.null(path) || !length(path) || is.na(path[[1L]]) || !nzchar(path[[1L]])) {
    return(NULL)
  }
  path <- normalizePath(as.character(path[[1L]]), mustWork = TRUE)
  anchors <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("anchor_id", .ivt_ref_parameter_names)
  missing <- setdiff(required, names(anchors))
  if (length(missing)) {
    stop("Anchor table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invalid <- .ivt_ref_parameter_names[vapply(.ivt_ref_parameter_names, function(name) {
    value <- suppressWarnings(as.numeric(anchors[[name]]))
    any(!is.finite(value))
  }, logical(1))]
  if (length(invalid)) {
    stop("Anchor table has non-finite parameter(s): ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  anchors
}

.ivt_ref_resolved_anchor_parameters <- function(anchors) {
  if (is.null(anchors) || !nrow(anchors)) return(data.frame())
  rows <- lapply(seq_len(nrow(anchors)), function(i) {
    inherited <- character()
    if ("inherited_fields" %in% names(anchors)) {
      text <- as.character(anchors$inherited_fields[[i]])
      if (!is.na(text) && nzchar(text)) inherited <- strsplit(text, ";", fixed = TRUE)[[1L]]
    }
    scanned <- data.frame(
      anchor_id = as.character(anchors$anchor_id[[i]]),
      parameter = .ivt_ref_parameter_names,
      value = vapply(.ivt_ref_parameter_names, function(name) {
        as.numeric(anchors[[name]][[i]])
      }, numeric(1)),
      value_origin = ifelse(
        .ivt_ref_parameter_names %in% inherited,
        "current_baseline",
        "source_row"
      ),
      unit = "model_parameter",
      is_scanned = TRUE,
      is_locked = FALSE,
      stringsAsFactors = FALSE
    )
    locked <- data.frame(
      anchor_id = as.character(anchors$anchor_id[[i]]),
      parameter = c("lam_max", "alpha_o2", "sigma_growth"),
      value = c(1.48, 0.5, 0.032),
      value_origin = "locked_setting",
      unit = "model_parameter",
      is_scanned = FALSE,
      is_locked = TRUE,
      stringsAsFactors = FALSE
    )
    rbind(scanned, locked)
  })
  do.call(rbind, rows)
}

.ivt_ref_design_tasks <- function(design) {
  lapply(seq_len(nrow(design)), function(i) {
    .ivt_scan_task(
      as.character(design$candidate_id[[i]]),
      as.character(design$stage[[i]]),
      unlist(design[i, .ivt_ref_parameter_names, drop = FALSE], use.names = TRUE)
    )
  })
}

.ivt_ref_attach_design_metadata <- function(summary, design_rows) {
  if (!identical(as.character(summary$candidate_id), as.character(design_rows$candidate_id))) {
    stop("Evaluator output order does not match the batch design.", call. = FALSE)
  }
  metadata <- c(
    "parent_candidate_id", "region", "is_anchor",
    "conditional_p_mis_base_upper", "parent_pool", "local_draw_index",
    "local_radius", "perturbed_parameter", "perturbation_direction",
    "requested_relative_change", "effective_relative_change"
  )
  metadata <- intersect(metadata, names(design_rows))
  for (name in metadata) summary[[name]] <- design_rows[[name]]
  summary
}

.ivt_ref_validate_checkpoint <- function(object, design, idx, batch_id, stage, digest) {
  if (!is.list(object) || !identical(as.integer(object$schema_version), .ivt_ref_schema_version)) {
    stop("Checkpoint batch ", batch_id, " has an invalid schema version.", call. = FALSE)
  }
  if (!identical(as.integer(object$batch_id), as.integer(batch_id)) ||
      !identical(as.character(object$stage), as.character(stage)) ||
      !identical(as.character(object$design_digest), as.character(digest))) {
    stop("Checkpoint batch ", batch_id, " has incompatible provenance.", call. = FALSE)
  }
  expected_ids <- as.character(design$candidate_id[idx])
  expected_parent <- as.character(design$parent_candidate_id[idx])
  if (!identical(as.character(object$candidate_ids), expected_ids) ||
      !identical(as.character(object$parent_candidate_ids), expected_parent) ||
      !is.data.frame(object$summary) ||
      !identical(as.character(object$summary$candidate_id), expected_ids) ||
      !identical(as.character(object$summary$parent_candidate_id), expected_parent)) {
    stop("Checkpoint batch ", batch_id, " membership/order differs from the design.", call. = FALSE)
  }
  TRUE
}

.ivt_ref_run_stage <- function(design,
                               stage,
                               out_dir,
                               cluster,
                               batch_size,
                               resume) {
  if (!is.data.frame(design) || !nrow(design)) {
    stop(stage, " design must contain at least one candidate.", call. = FALSE)
  }
  if (anyDuplicated(design$candidate_id) ||
      !all(c("candidate_id", "stage", "parent_candidate_id", .ivt_ref_parameter_names) %in% names(design))) {
    stop(stage, " design has an invalid schema or duplicate IDs.", call. = FALSE)
  }
  stage_dir <- file.path(out_dir, "checkpoints", stage)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  design_path <- file.path(out_dir, paste0(stage, "_design.rds"))
  design_digest <- .ivt_ref_digest_object(design)
  if (file.exists(design_path)) {
    prior <- readRDS(design_path)
    if (!isTRUE(resume) || !identical(prior, design)) {
      stop("Existing ", stage, " design differs from this invocation.", call. = FALSE)
    }
  } else {
    .ivt_ref_atomic_save_rds(design, design_path)
  }
  .ivt_scan_write_tsv(design, file.path(out_dir, paste0(stage, "_design.tsv")))

  batches <- split(seq_len(nrow(design)), ceiling(seq_len(nrow(design)) / batch_size))
  paths <- file.path(stage_dir, sprintf("batch_%04d.rds", seq_along(batches)))
  for (batch_id in seq_along(batches)) {
    idx <- batches[[batch_id]]
    checkpoint <- paths[[batch_id]]
    if (file.exists(checkpoint)) {
      if (!isTRUE(resume)) stop("Checkpoint already exists: ", checkpoint, call. = FALSE)
      .ivt_ref_validate_checkpoint(
        readRDS(checkpoint), design, idx, batch_id, stage, design_digest
      )
      message("[invitro_refinement] Resume verified ", stage, " batch ", batch_id,
              "/", length(batches))
      next
    }
    started <- Sys.time()
    rows <- .ivt_scan_evaluate_tasks(.ivt_ref_design_tasks(design[idx, , drop = FALSE]), cluster)
    summary <- .ivt_ref_attach_design_metadata(
      .ivt_scan_rbind(rows, "summary"), design[idx, , drop = FALSE]
    )
    object <- list(
      schema_version = .ivt_ref_schema_version,
      stage = stage,
      design_digest = design_digest,
      batch_id = batch_id,
      candidate_ids = as.character(design$candidate_id[idx]),
      parent_candidate_ids = as.character(design$parent_candidate_id[idx]),
      summary = summary,
      flow = .ivt_scan_rbind(rows, "flow"),
      kary = .ivt_scan_rbind(rows, "kary"),
      started_at = format(started, "%Y-%m-%d %H:%M:%S %Z"),
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
    .ivt_ref_atomic_save_rds(object, checkpoint)
    message("[invitro_refinement] Completed ", stage, " batch ", batch_id,
            "/", length(batches), " (", length(idx), " candidates)")
  }

  objects <- lapply(paths, readRDS)
  for (i in seq_along(objects)) {
    .ivt_ref_validate_checkpoint(objects[[i]], design, batches[[i]], i, stage, design_digest)
  }
  # Also exercise the shared membership validator against the complete ordered
  # checkpoint set.
  ivt_scan_validate_checkpoint_membership(
    objects, design, batches, stage = stage, design_digest = design_digest
  )
  summary <- do.call(rbind, lapply(objects, `[[`, "summary"))
  if (nrow(summary) != nrow(design) || anyDuplicated(summary$candidate_id) ||
      !identical(as.character(summary$candidate_id), as.character(design$candidate_id))) {
    stop(stage, " checkpoint aggregation is incomplete or out of order.", call. = FALSE)
  }
  list(
    summary = summary,
    flow = .ivt_scan_rbind(objects, "flow"),
    kary = .ivt_scan_rbind(objects, "kary"),
    design_digest = design_digest
  )
}

.ivt_ref_stage0_pass <- function(df, control_limit) {
  gate <- ivt_scan_hard_gate(df, control_limit)
  exact <- c(
    n_summary_rows = 114, n_growth = 219, n_death = 90,
    n_flow_samples = 20, n_ploidy_samples = 12,
    n_invalid_distributions = 0, n_low_o2_flow_groups = 6
  )
  for (name in names(exact)) {
    value <- suppressWarnings(as.numeric(df[[name]]))
    gate <- gate & is.finite(value) & value == exact[[name]]
  }
  wgd <- suppressWarnings(as.numeric(df$control_flow_wgd_mass))
  gate & is.finite(wgd) & wgd <= control_limit
}

.ivt_ref_stage0_design <- function(p_wgd_values, p_mis_values, stage, prefix) {
  grid <- expand.grid(
    p_wgd = as.numeric(p_wgd_values),
    p_mis_base = as.numeric(p_mis_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design <- data.frame(
    candidate_id = sprintf("%s_%05d", prefix, seq_len(nrow(grid))),
    stage = stage,
    region = prefix,
    parent_candidate_id = NA_character_,
    is_anchor = FALSE,
    conditional_p_mis_base_upper = NA_real_,
    grid,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (name in setdiff(.ivt_ref_parameter_names, c("p_wgd", "p_mis_base"))) {
    design[[name]] <- as.numeric(.ivt_scan_context$base_params[[name]])
  }
  design[, c(
    "candidate_id", "stage", "region", "parent_candidate_id", "is_anchor",
    "conditional_p_mis_base_upper", .ivt_ref_parameter_names
  ), drop = FALSE]
}

.ivt_ref_refinement_grid <- function(coarse_summary, p_wgd_values, n_refine) {
  rows <- list()
  for (i in seq_along(p_wgd_values)) {
    use <- coarse_summary[
      abs(as.numeric(coarse_summary$p_wgd) - p_wgd_values[[i]]) <
        .Machine$double.eps^0.5 * max(1, p_wgd_values[[i]]),
      , drop = FALSE
    ]
    use <- use[order(as.numeric(use$p_mis_base)), , drop = FALSE]
    pass <- as.logical(use$stage0_calibration_pass)
    first_fail <- which(!pass)[1L]
    if (!length(first_fail) || is.na(first_fail) || first_fail <= 1L || n_refine <= 0L) next
    lower <- as.numeric(use$p_mis_base[[first_fail - 1L]])
    upper <- as.numeric(use$p_mis_base[[first_fail]])
    values <- exp(seq(log(lower), log(upper), length.out = n_refine + 2L))
    values <- values[-c(1L, length(values))]
    rows[[length(rows) + 1L]] <- data.frame(
      p_wgd = p_wgd_values[[i]], p_mis_base = values,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

.ivt_ref_robustness_design <- function(candidates,
                                       specs,
                                       envelope,
                                       n_parents = 10L,
                                       relative_change = 0.05) {
  required <- c(
    "candidate_id", "full_target_pass", "pareto_front", "objective",
    .ivt_ref_parameter_names
  )
  missing <- setdiff(required, names(candidates))
  if (!is.data.frame(candidates) || length(missing)) {
    stop(
      "candidates is missing robustness column(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  n_parents <- as.integer(n_parents)
  if (length(n_parents) != 1L || is.na(n_parents) || n_parents <= 0L) {
    stop("n_parents must be a positive integer.", call. = FALSE)
  }
  relative_change <- .ivt_scan_as_number(
    relative_change, "relative_change", 0, 1
  )
  if (relative_change <= 0) stop("relative_change must be positive.", call. = FALSE)
  eligible_rows <- !is.na(candidates$full_target_pass) &
    as.logical(candidates$full_target_pass)
  if ("is_anchor" %in% names(candidates)) {
    eligible_rows <- eligible_rows &
      (is.na(candidates$is_anchor) | !as.logical(candidates$is_anchor))
  }
  eligible <- candidates[eligible_rows, , drop = FALSE]
  eligible_envelope <- envelope[
    !is.na(envelope$eligible) & as.logical(envelope$eligible),
    ,
    drop = FALSE
  ]
  if (!nrow(eligible_envelope)) stop("No eligible robustness envelope.", call. = FALSE)
  p_wgd_range <- range(as.numeric(eligible_envelope$p_wgd))
  spec_index <- setNames(seq_len(nrow(specs)), as.character(specs$param_symbol))
  robustness_eligible <- vapply(seq_len(nrow(eligible)), function(i) {
    base_values <- vapply(.ivt_ref_parameter_names, function(parameter) {
      as.numeric(eligible[[parameter]][[i]])
    }, numeric(1))
    all(vapply(.ivt_ref_parameter_names, function(name) {
      j <- spec_index[[name]]
      lower <- as.numeric(specs$lower[[j]])
      upper <- as.numeric(specs$upper[[j]])
      all(vapply(c(1 - relative_change, 1 + relative_change), function(factor) {
        values <- base_values
        values[[name]] <- values[[name]] * factor
        target_in_bounds <- is.finite(values[[name]]) &&
          values[[name]] >= lower && values[[name]] <= upper
        p_wgd_in_envelope <- is.finite(values[["p_wgd"]]) &&
          values[["p_wgd"]] >= p_wgd_range[[1L]] &&
          values[["p_wgd"]] <= p_wgd_range[[2L]]
        if (!target_in_bounds || !p_wgd_in_envelope) return(FALSE)
        cap <- .ivt_scan_conditional_cap(values[["p_wgd"]], envelope)
        p_mis_idx <- spec_index[["p_mis_base"]]
        is.finite(values[["p_mis_base"]]) &&
          values[["p_mis_base"]] >= as.numeric(specs$lower[[p_mis_idx]]) &&
          values[["p_mis_base"]] <= as.numeric(specs$upper[[p_mis_idx]]) &&
          values[["p_mis_base"]] <= cap
      }, logical(1)))
    }, logical(1)))
  }, logical(1))
  eligible <- eligible[robustness_eligible, , drop = FALSE]
  if (nrow(eligible) < n_parents) return(data.frame())
  eligible <- eligible[order(
    !as.logical(eligible$pareto_front),
    suppressWarnings(as.numeric(eligible$objective)),
    as.character(eligible$candidate_id),
    na.last = TRUE
  ), , drop = FALSE]
  parents <- eligible[seq_len(n_parents), , drop = FALSE]

  rows <- vector("list", n_parents * length(.ivt_ref_parameter_names) * 2L)
  out_index <- 0L
  for (i in seq_len(nrow(parents))) {
    for (name in .ivt_ref_parameter_names) {
      for (direction in c("minus", "plus")) {
        out_index <- out_index + 1L
        factor <- if (direction == "minus") 1 - relative_change else 1 + relative_change
        values <- vapply(.ivt_ref_parameter_names, function(parameter) {
          as.numeric(parents[[parameter]][[i]])
        }, numeric(1))
        original <- values[[name]]
        values[[name]] <- original * factor
        cap <- .ivt_scan_conditional_cap(values[["p_wgd"]], envelope)
        safe_parent <- gsub("[^A-Za-z0-9]+", "_", as.character(parents$candidate_id[[i]]))
        rows[[out_index]] <- data.frame(
          candidate_id = paste("robust", safe_parent, name, direction, sep = "_"),
          stage = "robustness",
          region = "one_at_a_time_robustness",
          parent_candidate_id = as.character(parents$candidate_id[[i]]),
          parent_pool = "target_pass",
          local_draw_index = NA_integer_,
          local_radius = relative_change,
          is_anchor = FALSE,
          conditional_p_mis_base_upper = cap,
          perturbed_parameter = name,
          perturbation_direction = direction,
          requested_relative_change = if (direction == "minus") {
            -relative_change
          } else {
            relative_change
          },
          effective_relative_change = if (original != 0) {
            values[[name]] / original - 1
          } else {
            NA_real_
          },
          as.list(values),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  design <- do.call(rbind, rows)
  rownames(design) <- NULL
  design
}

.ivt_ref_mark_source_bounds <- function(df, parameter_table) {
  pass <- rep(TRUE, nrow(df))
  for (name in .ivt_ref_parameter_names) {
    row <- parameter_table[as.character(parameter_table$param_symbol) == name, , drop = FALSE]
    value <- suppressWarnings(as.numeric(df[[name]]))
    inside <- rep(FALSE, nrow(df))
    if (nrow(row) == 1L) {
      lower <- suppressWarnings(as.numeric(row$lower_bound[[1L]]))
      upper <- suppressWarnings(as.numeric(row$upper_bound[[1L]]))
      inside <- is.finite(value) & is.finite(lower) & is.finite(upper) &
        value >= lower & value <= upper
    }
    df[[paste0("within_source_bounds_", name)]] <- inside
    pass <- pass & inside
  }
  df$within_source_parameter_bounds <- pass
  df
}

.ivt_ref_annotate <- function(summary, parameter_table, control_limit) {
  summary <- .ivt_ref_mark_source_bounds(summary, parameter_table)
  summary <- ivt_scan_target_gate(summary, control_limit = control_limit)
  summary$hard_pass <- summary$technical_pass
  metrics <- c(
    "abs_error_2N_A19_mean_N", "abs_error_2N_A19_mass_ge70",
    "abs_error_4N_A17_mean_N", "flow_low_o2_mass_ge3_rmse",
    "flow_low_o2_wgd_mass_rmse", "flow_low_o2_bimodal_count_abs_error",
    "flow_nll", "growth_log_residual_rmse", "death_logit_residual_rmse"
  )
  summary$pareto_front <- FALSE
  eligible <- which(summary$technical_pass)
  if (length(eligible)) {
    summary$pareto_front[eligible] <-
      ivt_scan_pareto_front(summary[eligible, , drop = FALSE], metrics)
  }
  order_idx <- order(
    !summary$full_target_pass, !summary$technical_pass, !summary$pareto_front,
    summary$flow_low_o2_wgd_mass_rmse, summary$abs_error_2N_A19_mean_N,
    summary$objective, na.last = TRUE
  )
  summary[order_idx, , drop = FALSE]
}

.ivt_ref_write_stage_outputs <- function(stage, annotated, flow, kary, out_dir) {
  prefix <- file.path(out_dir, stage)
  .ivt_scan_write_tsv(annotated, paste0(prefix, "_summary.tsv"))
  .ivt_scan_write_tsv(
    annotated[annotated$technical_pass, , drop = FALSE], paste0(prefix, "_hard_pass.tsv")
  )
  .ivt_scan_write_tsv(
    annotated[annotated$technical_pass & annotated$pareto_front, , drop = FALSE],
    paste0(prefix, "_pareto.tsv")
  )
  .ivt_scan_write_tsv(
    annotated[annotated$full_target_pass, , drop = FALSE], paste0(prefix, "_target.tsv")
  )
  failure_cols <- intersect(c(
    "candidate_id", "stage", "status", "error", "failed_cohort",
    "failed_scenario_id", "failed_segment_id", "failed_segment_ordinal",
    "failed_threshold_target_cells", "failed_max_live_cells",
    "failed_relative_threshold_shortfall", "passage_failure_reason"
  ), names(annotated))
  .ivt_scan_write_tsv(
    annotated[as.character(annotated$status) != "OK", failure_cols, drop = FALSE],
    paste0(prefix, "_failures.tsv")
  )
  .ivt_scan_write_tsv(flow, paste0(prefix, "_flow_metrics.tsv"))
  .ivt_scan_write_tsv(kary, paste0(prefix, "_kary_metrics.tsv"))
  invisible(annotated)
}

.ivt_ref_contract <- function(paths, settings, observed_contract) {
  code_paths <- c(
    file.path(.ref_script_dir, "scan_invitro_observed_feasibility_refinement.R"),
    file.path(.ref_script_dir, "scan_invitro_observed_feasibility.R"),
    file.path(.ref_script_dir, "invitro_observed_feasibility_utils.R")
  )
  workflow_files <- c(
    list.files(file.path(WORKFLOW_ROOT, "model"), pattern = "\\.(R|cpp)$",
               full.names = TRUE, recursive = TRUE),
    list.files(file.path(WORKFLOW_ROOT, "util"), pattern = "\\.R$",
               full.names = TRUE, recursive = TRUE)
  )
  workflow_files <- workflow_files[!grepl(
    paste0("[\\/]", "\\.rcpp_cache_"),
    workflow_files
  )]
  fit_files <- list.files(paths$fit_objects_dir, full.names = TRUE, recursive = TRUE,
                          all.files = TRUE, no.. = TRUE)
  fit_files <- fit_files[file.info(fit_files)$isdir %in% FALSE]
  inputs <- rbind(
    .ivt_scan_file_inventory("scan_code", code_paths),
    .ivt_scan_file_inventory("workflow_code", workflow_files),
    .ivt_scan_file_inventory("fit_result", paths$fit_result),
    .ivt_scan_file_inventory("parameter_table", paths$parameter_table),
    .ivt_scan_file_inventory("fit_object", fit_files),
    .ivt_scan_file_inventory("flow_density", paths$flow_density),
    .ivt_scan_file_inventory("death_data", paths$death_data)
  )
  if (!is.null(paths$anchors)) {
    inputs <- rbind(inputs, .ivt_scan_file_inventory("anchors", paths$anchors))
  }
  list(
    schema_version = .ivt_ref_schema_version,
    inputs = inputs,
    settings = data.frame(
      key = names(settings),
      value = vapply(settings, function(x) paste(as.character(x), collapse = ","), character(1)),
      stringsAsFactors = FALSE
    ),
    observed_contract = observed_contract
  )
}

.ivt_ref_prepare_contract <- function(contract, out_dir, resume) {
  path <- file.path(out_dir, "run_contract.rds")
  reusable <- c(
    list.files(file.path(out_dir, "checkpoints"), pattern = "\\.rds$",
               recursive = TRUE, full.names = TRUE),
    list.files(out_dir, pattern = "_(design|envelope)\\.rds$", full.names = TRUE)
  )
  if (file.exists(path)) {
    if (!isTRUE(resume)) stop("resume=FALSE requires a new output directory.", call. = FALSE)
    if (!identical(readRDS(path), contract)) {
      stop("Run contract differs; refusing to reuse outputs.", call. = FALSE)
    }
  } else {
    if (isTRUE(resume) && length(reusable)) {
      stop("Reusable outputs exist without run_contract.rds.", call. = FALSE)
    }
    .ivt_ref_atomic_save_rds(contract, path)
  }
  .ivt_scan_write_tsv(contract$inputs, file.path(out_dir, "run_contract_inputs.tsv"))
  .ivt_scan_write_tsv(contract$settings, file.path(out_dir, "run_contract_settings.tsv"))
  .ivt_scan_write_tsv(contract$observed_contract, file.path(out_dir, "acceptance_contract.tsv"))
}

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_result_path <- normalizePath(
    .ivt_scan_first(argv$fit_result, stop("--fit_result is required.", call. = FALSE)),
    mustWork = TRUE
  )
  fit_dir <- dirname(fit_result_path)
  out_dir <- normalizePath(
    .ivt_scan_first(argv$out_dir, stop("--out_dir is required.", call. = FALSE)),
    mustWork = FALSE
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  n_focused <- .ivt_scan_as_integer(.ivt_scan_first(argv$n_focused, 15000L), "n_focused", 0L)
  n_outer <- .ivt_scan_as_integer(.ivt_scan_first(argv$n_outer, 5000L), "n_outer", 0L)
  n_stage2 <- .ivt_scan_as_integer(.ivt_scan_first(argv$n_stage2, 5000L), "n_stage2", 0L)
  if (n_focused + n_outer < 1L) stop("Stage 1 requires at least one sampled candidate.", call. = FALSE)
  n_cores <- .ivt_scan_as_integer(.ivt_scan_first(argv$n_cores, 1L), "n_cores")
  batch_size <- .ivt_scan_as_integer(
    .ivt_scan_first(argv$batch_size, max(5L * n_cores, n_cores)), "batch_size"
  )
  seed <- .ivt_scan_as_integer(.ivt_scan_first(argv$design_seed, 10L), "design_seed", 0L)
  stage0_nodes <- .ivt_scan_as_integer(
    .ivt_scan_first(argv$stage0_p_wgd_points, 17L), "stage0_p_wgd_points", 2L
  )
  stage0_coarse <- .ivt_scan_as_integer(
    .ivt_scan_first(argv$stage0_p_mis_points, 31L), "stage0_p_mis_points", 3L
  )
  stage0_refine <- .ivt_scan_as_integer(
    .ivt_scan_first(argv$stage0_refine_points, 15L), "stage0_refine_points", 0L
  )
  max_parents <- .ivt_scan_as_integer(
    .ivt_scan_first(argv$max_stage2_parents, 100L), "max_stage2_parents"
  )
  local_fraction <- .ivt_scan_as_number(
    .ivt_scan_first(argv$stage2_neighborhood_fraction, 0.05),
    "stage2_neighborhood_fraction", 0, 1
  )
  if (local_fraction <= 0) {
    stop("stage2_neighborhood_fraction must be strictly positive.", call. = FALSE)
  }
  control_limit <- .ivt_scan_as_number(
    .ivt_scan_first(argv$control_high_mass_limit, 0.05),
    "control_high_mass_limit", 0, 1
  )
  retries <- .ivt_scan_as_integer(
    .ivt_scan_first(argv$max_preflight_retries, 5L), "max_preflight_retries", 0L
  )
  resume <- .ivt_scan_as_bool(argv$resume, TRUE)

  fit_result <- readRDS(fit_result_path)
  if (is.null(fit_result$cfg) || is.null(fit_result$best_params)) {
    stop("fit_result must contain cfg and best_params.", call. = FALSE)
  }
  parameter_table_path <- .ivt_scan_resolve_input(
    argv$parameter_table, file.path(fit_dir, "parameter_table_input.csv"),
    file.path(OXYGEN_ROOT, "results", "fit_invitro_O2_buffering_500seed", "seed10",
              "parameter_table_input.csv"), "parameter table"
  )
  fit_objects_dir <- .ivt_scan_resolve_input(
    argv$fit_objects_dir, fit_result$fit_objects_dir,
    file.path(OXYGEN_ROOT, "ploidyOxygen", "data", "fit_objects"),
    "fit objects directory", directory = TRUE
  )
  flow_density_path <- .ivt_scan_resolve_input(
    argv$flow_density_path, fit_result$flow_density_path,
    file.path(OXYGEN_ROOT, "data", "g0g1_ploidy_density_grid.csv"), "flow-density data"
  )
  death_data_path <- .ivt_scan_resolve_input(
    argv$death_data_path, fit_result$death_data_path,
    file.path(PROJECT_ROOT, "data", "InVitroData",
              "sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"),
    "death-likelihood data"
  )
  anchors_path <- NULL
  if (!is.null(argv$anchors) && length(argv$anchors) && nzchar(argv$anchors[[1L]])) {
    anchors_path <- normalizePath(argv$anchors[[1L]], mustWork = TRUE)
  }
  anchors <- .ivt_ref_read_anchors(anchors_path)
  parameter_table <- utils::read.csv(
    parameter_table_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  fit_objects <- ivt_load_fit_objects(
    repo_root = OXYGEN_ROOT,
    fit_objects_dir = fit_objects_dir,
    flow_csv_path = flow_density_path,
    death_data_path = death_data_path
  )

  base_params <- fit_result$best_params
  base_params$lam_max <- 1.48
  base_params$alpha_o2 <- 0.5
  base_params$sigma_growth <- 0.032
  required_base_params <- unique(c(
    .ivt_ref_parameter_names, "lam_max", "alpha_o2", "sigma_growth",
    "init_mean_4N"
  ))
  invalid_base_params <- required_base_params[!vapply(required_base_params, function(name) {
    value <- suppressWarnings(as.numeric(base_params[[name]]))
    length(value) == 1L && is.finite(value)
  }, logical(1))]
  if (length(invalid_base_params)) {
    stop(
      "fit_result/base overrides are missing finite parameter(s): ",
      paste(invalid_base_params, collapse = ", "), ".",
      call. = FALSE
    )
  }
  .ivt_scan_context <<- list(
    fit_objects = fit_objects,
    cfg = fit_result$cfg,
    base_params = base_params,
    targets = list(
      pooled_2N_A19_mean_N = 77.45,
      pooled_2N_A19_mass_ge70 = 0.625,
      pooled_4N_A17_mean_N = 81.475
    ),
    control_limit = control_limit,
    ploidy_weight = as.numeric(fit_result$ploidy_weight %||% 1),
    death_weight = as.numeric(fit_result$death_weight %||% 1),
    passage_time_weight = as.numeric(fit_result$passage_time_weight %||% 0.25),
    passage_time_tolerance_days = as.numeric(fit_result$passage_time_tolerance_days %||% 1),
    passage_time_sigma_days = as.numeric(fit_result$passage_time_sigma_days %||% 1),
    passage_time_df = as.numeric(fit_result$passage_time_df %||% 4),
    sigma_death_logit = as.numeric(fit_result$sigma_death_logit %||% 0.75),
    death_fraction_eps = as.numeric(fit_result$death_fraction_eps %||% 1e-4)
  )

  # Master-only observed replay and acceptance contract happen before the fork
  # cluster exists, so every fork inherits the exact observed targets.
  observed_comp <- do.call(ivt_objective_components, .ivt_scan_objective_args(base_params))
  observed_targets <- .ivt_scan_observed_targets(observed_comp, fit_objects)
  .ivt_scan_context$targets <- observed_targets$values
  observed_contract <- .ivt_scan_contract_table(observed_comp, observed_targets$table)
  if (!all(observed_contract$pass)) {
    .ivt_scan_write_tsv(observed_contract, file.path(out_dir, "acceptance_contract.tsv"))
    stop("Observed-data contract failed: ",
         paste(observed_contract$check[!observed_contract$pass], collapse = ", "),
         call. = FALSE)
  }

  focused_specs <- ivt_scan_extended_parameter_specs("focused", parameter_table)
  outer_specs <- ivt_scan_extended_parameter_specs("outer", parameter_table)
  render_specs <- function(specs) {
    paste(
      paste(
        specs$param_symbol, specs$lower, specs$upper, specs$transform,
        sep = ":"
      ),
      collapse = ";"
    )
  }
  settings <- list(
    n_focused = n_focused, n_outer = n_outer, n_stage2 = n_stage2,
    n_cores = n_cores, batch_size = batch_size, design_seed = seed,
    stage1_seed = seed, stage2_seed = seed + 2L,
    stage0_p_wgd_points = stage0_nodes, stage0_p_mis_points = stage0_coarse,
    stage0_p_wgd_lower = 1e-4, stage0_p_wgd_upper = 2.5e-3,
    stage0_p_mis_lower = 1e-7, stage0_p_mis_upper = 5e-3,
    stage0_refine_points = stage0_refine, max_stage2_parents = max_parents,
    stage2_neighborhood_fraction = local_fraction,
    robustness_min_target_parents = 10,
    robustness_relative_change = 0.05,
    robustness_parameters = paste(.ivt_ref_parameter_names, collapse = ";"),
    control_high_mass_limit = control_limit, max_preflight_retries = retries,
    expected_n_scenarios = 6, expected_n_summary_rows = 114,
    expected_n_growth = 219, expected_n_death = 90,
    expected_n_flow_samples = 20, expected_n_ploidy_samples = 12,
    expected_n_low_o2_flow_groups = 6,
    target_wgd_mass_abs_tolerance = 0.15,
    target_2N_A19_mean = 77.45, target_2N_A19_mean_tolerance = 8,
    target_2N_A19_mass_ge70 = 0.625,
    target_2N_A19_mass_ge70_tolerance = 0.20,
    target_4N_A17_mean = 81.475, target_4N_A17_mean_tolerance = 4,
    target_4N_decline_lower = 0, target_4N_decline_upper = 6,
    focused_parameter_specs = render_specs(focused_specs),
    outer_parameter_specs = render_specs(outer_specs),
    conditional_envelope_rule = "minimum_adjacent_eligible_caps",
    fixed_lam_max = 1.48, fixed_alpha_o2 = 0.5, fixed_sigma_growth = 0.032,
    no_de = TRUE, no_new_likelihood = TRUE, no_soft_constraint = TRUE,
    container_image = Sys.getenv("O2SD_CONTAINER_IMAGE", unset = "")
  )
  contract <- .ivt_ref_contract(
    list(
      fit_result = fit_result_path, parameter_table = parameter_table_path,
      fit_objects_dir = fit_objects_dir, flow_density = flow_density_path,
      death_data = death_data_path, anchors = anchors_path
    ),
    settings,
    observed_contract
  )
  .ivt_ref_prepare_contract(contract, out_dir, resume)
  .ivt_scan_write_tsv(observed_targets$table, file.path(out_dir, "observed_targets.tsv"))
  if (!is.null(anchors)) {
    .ivt_scan_write_tsv(anchors, file.path(out_dir, "anchors_input.tsv"))
    .ivt_scan_write_tsv(
      .ivt_ref_resolved_anchor_parameters(anchors),
      file.path(out_dir, "anchor_resolved_parameters.tsv")
    )
  }

  probe_task <- .ivt_scan_task(
    "preflight_anchor", "preflight",
    unlist(as.list(base_params[.ivt_ref_parameter_names]), use.names = TRUE)
  )
  cluster_state <- .ivt_scan_start_cluster(
    n_cores, probe_task, file.path(out_dir, "worker_preflight.tsv"), retries
  )
  cluster <- cluster_state$cluster
  if (!is.null(cluster)) on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)

  # Stage 0: joint normoxic admissibility envelope in p_wgd x p_mis_base.
  p_wgd_values <- exp(seq(log(1e-4), log(2.5e-3), length.out = stage0_nodes))
  p_mis_values <- exp(seq(log(1e-7), log(5e-3), length.out = stage0_coarse))
  coarse_design <- .ivt_ref_stage0_design(
    p_wgd_values, p_mis_values, "stage0_coarse", "stage0_coarse"
  )
  coarse_result <- .ivt_ref_run_stage(
    coarse_design, "stage0_coarse", out_dir, cluster, batch_size, resume
  )
  coarse_summary <- coarse_result$summary
  coarse_summary$stage0_calibration_pass <- .ivt_ref_stage0_pass(coarse_summary, control_limit)

  refine_grid <- .ivt_ref_refinement_grid(coarse_summary, p_wgd_values, stage0_refine)
  refined_summary <- coarse_summary[FALSE, , drop = FALSE]
  if (nrow(refine_grid)) {
    refine_design <- coarse_design[FALSE, , drop = FALSE]
    refine_design <- refine_design[rep(NA_integer_, nrow(refine_grid)), , drop = FALSE]
    refine_design$candidate_id <- sprintf("stage0_refine_%05d", seq_len(nrow(refine_grid)))
    refine_design$stage <- "stage0_refine"
    refine_design$region <- "stage0_refine"
    refine_design$parent_candidate_id <- NA_character_
    refine_design$is_anchor <- FALSE
    refine_design$conditional_p_mis_base_upper <- NA_real_
    refine_design$p_wgd <- refine_grid$p_wgd
    refine_design$p_mis_base <- refine_grid$p_mis_base
    for (name in setdiff(.ivt_ref_parameter_names, c("p_wgd", "p_mis_base"))) {
      refine_design[[name]] <- as.numeric(base_params[[name]])
    }
    refine_result <- .ivt_ref_run_stage(
      refine_design, "stage0_refine", out_dir, cluster, batch_size, resume
    )
    refined_summary <- refine_result$summary
    refined_summary$stage0_calibration_pass <-
      .ivt_ref_stage0_pass(refined_summary, control_limit)
  }
  calibration <- rbind(coarse_summary, refined_summary)
  calibration <- calibration[order(as.numeric(calibration$p_wgd),
                                     as.numeric(calibration$p_mis_base)), , drop = FALSE]
  # The shared selector derives contiguous-prefix caps from hard_pass.  For the
  # joint calibration, hard_pass is the stricter Stage-0 technical/control gate.
  calibration$hard_pass <- calibration$stage0_calibration_pass
  calibration_for_envelope <- calibration
  stage0_fail <- is.na(calibration_for_envelope$stage0_calibration_pass) |
    !calibration_for_envelope$stage0_calibration_pass
  calibration_for_envelope$status[stage0_fail] <-
    "STAGE0_GATE_FAIL"
  envelope <- ivt_scan_select_conditional_envelope(
    calibration_for_envelope, control_limit
  )
  .ivt_scan_write_tsv(calibration, file.path(out_dir, "stage0_calibration.tsv"))
  .ivt_scan_write_tsv(envelope, file.path(out_dir, "stage0_conditional_envelope.tsv"))
  .ivt_ref_atomic_save_rds(envelope, file.path(out_dir, "stage0_conditional_envelope.rds"))
  if (!any(as.logical(envelope$eligible))) {
    writeLines(c("status=NO_ADMISSIBLE_ENVELOPE", paste0("completed_at=", Sys.time())),
               file.path(out_dir, "status.log"))
    stop("Stage 0 found no contiguous admissible envelope.", call. = FALSE)
  }

  # Stage 1: 15k focused + 5k outer LHS by default, plus forced audit anchors.
  stage1_design <- ivt_scan_extended_design(
    n_focused, n_outer, seed, envelope, anchors = anchors
  )
  stage1_design$parent_pool <- NA_character_
  stage1_design$local_draw_index <- NA_integer_
  stage1_design$local_radius <- NA_real_
  stage1_result <- .ivt_ref_run_stage(
    stage1_design, "stage1", out_dir, cluster, batch_size, resume
  )
  stage1 <- .ivt_ref_annotate(stage1_result$summary, parameter_table, control_limit)
  .ivt_ref_write_stage_outputs(
    "stage1", stage1, stage1_result$flow, stage1_result$kary, out_dir
  )

  # Stage 2: local transformed-space refinement.  No-parent is a valid,
  # explicitly recorded diagnostic result rather than a crash.
  stage2_status <- "COMPLETE"
  stage2_result <- list(
    summary = stage1_result$summary[FALSE, , drop = FALSE],
    flow = stage1_result$flow[FALSE, , drop = FALSE],
    kary = stage1_result$kary[FALSE, , drop = FALSE]
  )
  stage2 <- stage1[FALSE, , drop = FALSE]
  if (n_stage2 > 0L && any(stage1$technical_pass)) {
    specs <- ivt_scan_extended_parameter_specs("outer")
    stage2_design <- tryCatch(
      ivt_scan_stage2_design(
        stage1, specs, n_stage2, seed + 2L, envelope,
        max_parents = max_parents, neighborhood_fraction = local_fraction
      ),
      error = function(e) {
        if (grepl("No refinement-eligible parent", conditionMessage(e), fixed = TRUE)) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (is.null(stage2_design)) {
      stage2_status <- "NO_REFINEMENT_ELIGIBLE"
    } else {
      parent_ids <- unique(as.character(stage2_design$parent_candidate_id))
      refinement_parents <- stage1[
        match(parent_ids, as.character(stage1$candidate_id)),
        ,
        drop = FALSE
      ]
      refinement_parents$refinement_parent_pool <- vapply(parent_ids, function(id) {
        as.character(stage2_design$parent_pool[match(id, stage2_design$parent_candidate_id)])
      }, character(1))
      .ivt_scan_write_tsv(
        refinement_parents,
        file.path(out_dir, "refinement_parents.tsv")
      )
      stage2_result <- .ivt_ref_run_stage(
        stage2_design, "stage2", out_dir, cluster, batch_size, resume
      )
      stage2 <- .ivt_ref_annotate(stage2_result$summary, parameter_table, control_limit)
    }
  } else {
    stage2_status <- if (n_stage2 <= 0L) {
      "SKIPPED_N_STAGE2_ZERO"
    } else {
      "NO_REFINEMENT_ELIGIBLE"
    }
  }
  if (!identical(stage2_status, "COMPLETE")) {
    empty_refinement_parents <- stage1[FALSE, , drop = FALSE]
    empty_refinement_parents$refinement_parent_pool <- character()
    .ivt_scan_write_tsv(
      empty_refinement_parents,
      file.path(out_dir, "refinement_parents.tsv")
    )
    empty_stage2_design <- stage1_design[FALSE, , drop = FALSE]
    .ivt_scan_write_tsv(empty_stage2_design, file.path(out_dir, "stage2_design.tsv"))
    empty_stage2_rds <- file.path(out_dir, "stage2_design.rds")
    if (file.exists(empty_stage2_rds)) {
      if (!isTRUE(resume) || !identical(readRDS(empty_stage2_rds), empty_stage2_design)) {
        stop("Existing empty Stage-2 design differs from this invocation.", call. = FALSE)
      }
    } else {
      .ivt_ref_atomic_save_rds(empty_stage2_design, empty_stage2_rds)
    }
    .ivt_scan_write_tsv(
      data.frame(
        status = stage2_status,
        reason = if (n_stage2 <= 0L) {
          "n_stage2_is_zero"
        } else {
          "no_refinement_parent_pool"
        },
        stringsAsFactors = FALSE
      ),
      file.path(out_dir, "stage2_status.tsv")
    )
  }
  .ivt_ref_write_stage_outputs(
    "stage2", stage2, stage2_result$flow, stage2_result$kary, out_dir
  )

  combined_raw <- rbind(stage1_result$summary, stage2_result$summary)
  combined <- .ivt_ref_annotate(combined_raw, parameter_table, control_limit)
  combined_flow <- rbind(stage1_result$flow, stage2_result$flow)
  combined_kary <- rbind(stage1_result$kary, stage2_result$kary)
  .ivt_ref_write_stage_outputs("combined", combined, combined_flow, combined_kary, out_dir)

  # Conditional one-at-a-time robustness audit.  It is diagnostic only and
  # never creates a DE table or launches an optimizer.
  robustness_status <- "NOT_TRIGGERED_INSUFFICIENT_TARGET_PASS"
  robustness_gate_pass <- FALSE
  robustness_result <- list(
    summary = stage1_result$summary[FALSE, , drop = FALSE],
    flow = stage1_result$flow[FALSE, , drop = FALSE],
    kary = stage1_result$kary[FALSE, , drop = FALSE]
  )
  robustness <- combined[FALSE, , drop = FALSE]
  target_nonanchor <- !is.na(combined$full_target_pass) &
    as.logical(combined$full_target_pass) &
    (is.na(combined$is_anchor) | !as.logical(combined$is_anchor))
  robustness_design <- .ivt_ref_robustness_design(
    combined, outer_specs, envelope, n_parents = 10L, relative_change = 0.05
  )
  if (nrow(robustness_design)) {
    robustness_result <- .ivt_ref_run_stage(
      robustness_design, "robustness", out_dir, cluster, batch_size, resume
    )
    robustness <- .ivt_ref_annotate(
      robustness_result$summary, parameter_table, control_limit
    )
    by_parent <- lapply(split(
      seq_len(nrow(robustness)),
      as.character(robustness$parent_candidate_id)
    ), function(idx) {
      data.frame(
        parent_candidate_id = as.character(robustness$parent_candidate_id[idx[[1L]]]),
        n_perturbations = length(idx),
        n_exact_requested_perturbations = sum(
          abs(
            as.numeric(robustness$effective_relative_change[idx]) -
              as.numeric(robustness$requested_relative_change[idx])
          ) <= 1e-12
        ),
        n_technical_pass = sum(robustness$technical_pass[idx]),
        n_full_target_pass = sum(robustness$full_target_pass[idx]),
        all_perturbations_pass = all(robustness$full_target_pass[idx]),
        stringsAsFactors = FALSE
      )
    })
    robustness_by_parent <- do.call(rbind, by_parent)
    rownames(robustness_by_parent) <- NULL
    robustness_gate_pass <- nrow(robustness_by_parent) >= 10L &&
      all(robustness_by_parent$n_perturbations == 2L * length(.ivt_ref_parameter_names)) &&
      all(
        robustness_by_parent$n_exact_requested_perturbations ==
          robustness_by_parent$n_perturbations
      ) &&
      all(robustness_by_parent$all_perturbations_pass)
    robustness_status <- if (robustness_gate_pass) {
      "COMPLETE_PASS"
    } else {
      "COMPLETE_FAIL"
    }
  } else {
    robustness_by_parent <- data.frame(
      parent_candidate_id = character(),
      n_perturbations = integer(),
      n_exact_requested_perturbations = integer(),
      n_technical_pass = integer(),
      n_full_target_pass = integer(),
      all_perturbations_pass = logical(),
      stringsAsFactors = FALSE
    )
    empty_robustness_design <- stage1_design[FALSE, , drop = FALSE]
    empty_robustness_design$perturbed_parameter <- character()
    empty_robustness_design$perturbation_direction <- character()
    empty_robustness_design$requested_relative_change <- numeric()
    empty_robustness_design$effective_relative_change <- numeric()
    .ivt_scan_write_tsv(
      empty_robustness_design,
      file.path(out_dir, "robustness_design.tsv")
    )
    empty_robustness_rds <- file.path(out_dir, "robustness_design.rds")
    if (file.exists(empty_robustness_rds)) {
      if (!isTRUE(resume) ||
          !identical(readRDS(empty_robustness_rds), empty_robustness_design)) {
        stop("Existing empty robustness design differs from this invocation.", call. = FALSE)
      }
    } else {
      .ivt_ref_atomic_save_rds(empty_robustness_design, empty_robustness_rds)
    }
  }
  .ivt_ref_write_stage_outputs(
    "robustness", robustness, robustness_result$flow,
    robustness_result$kary, out_dir
  )
  .ivt_scan_write_tsv(
    robustness_by_parent,
    file.path(out_dir, "robustness_by_parent.tsv")
  )
  .ivt_scan_write_tsv(
    data.frame(
      status = robustness_status,
      n_available_nonanchor_target_pass = sum(target_nonanchor),
      n_parent_candidates = nrow(robustness_by_parent),
      n_perturbation_evaluations = nrow(robustness),
      robustness_gate_pass = robustness_gate_pass,
      de_started = FALSE,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "robustness_status.tsv")
  )

  specs <- rbind(focused_specs, outer_specs)
  specs$exploratory_bounds <- TRUE
  .ivt_scan_write_tsv(specs, file.path(out_dir, "extended_parameter_specs.tsv"))

  checksum_paths <- c(
    file.path(.ref_script_dir, "scan_invitro_observed_feasibility_refinement.R"),
    file.path(.ref_script_dir, "scan_invitro_observed_feasibility.R"),
    file.path(.ref_script_dir, "invitro_observed_feasibility_utils.R"),
    fit_result_path, parameter_table_path, flow_density_path, death_data_path
  )
  checksum_names <- c(
    "refinement_driver", "base_scan_driver", "scan_helpers", "fit_result",
    "parameter_table", "flow_density", "death_data"
  )
  if (!is.null(anchors_path)) {
    checksum_paths <- c(checksum_paths, anchors_path)
    checksum_names <- c(checksum_names, "anchors")
  }
  checksums <- data.frame(
    input = checksum_names,
    path = normalizePath(checksum_paths, mustWork = TRUE),
    sha256 = vapply(checksum_paths, .ivt_scan_sha256, character(1)),
    stringsAsFactors = FALSE
  )
  .ivt_scan_write_tsv(checksums, file.path(out_dir, "input_checksums.tsv"))
  manifest <- data.frame(
    key = c(
      "schema_version", "generated_at", "project_root", "git_branch", "git_head",
      "git_status", "fit_result", "parameter_table", "fit_objects_dir",
      "flow_density_path", "death_data_path", "anchors_path", "n_stage0",
      "n_stage1", "n_stage2", "stage2_status", "n_technical_pass",
      "n_full_target_pass", "n_pareto", "robustness_status",
      "n_robustness_evaluations", "n_robustness_parents",
      "robustness_gate_pass", "fixed_lam_max", "fixed_alpha_o2",
      "fixed_sigma_growth", "de_started"
    ),
    value = c(
      .ivt_ref_schema_version, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      PROJECT_ROOT, .ivt_scan_git_value(c("branch", "--show-current")),
      .ivt_scan_git_value(c("rev-parse", "HEAD")),
      .ivt_scan_git_value(c("status", "--short")), fit_result_path,
      parameter_table_path, fit_objects_dir, flow_density_path, death_data_path,
      anchors_path %||% "", nrow(calibration), nrow(stage1), nrow(stage2),
      stage2_status, sum(combined$technical_pass), sum(combined$full_target_pass),
      sum(combined$technical_pass & combined$pareto_front), robustness_status,
      nrow(robustness), nrow(robustness_by_parent), robustness_gate_pass,
      1.48, 0.5, 0.032, FALSE
    ),
    stringsAsFactors = FALSE
  )
  .ivt_scan_write_tsv(manifest, file.path(out_dir, "manifest.tsv"))
  writeLines(
    c(
      "status=COMPLETE", paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste0("stage2_status=", stage2_status), paste0("n_stage1=", nrow(stage1)),
      paste0("n_stage2=", nrow(stage2)),
      paste0("n_technical_pass=", sum(combined$technical_pass)),
      paste0("n_full_target_pass=", sum(combined$full_target_pass)),
      paste0("robustness_status=", robustness_status),
      paste0("robustness_gate_pass=", robustness_gate_pass),
      "de_started=FALSE"
    ),
    file.path(out_dir, "status.log")
  )
  message("[invitro_refinement] COMPLETE: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
