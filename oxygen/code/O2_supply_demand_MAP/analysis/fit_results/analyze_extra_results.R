#!/usr/bin/env Rscript

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE), character(1)))
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
  }
})
source(file.path(.script_dir, "extra_results_analysis_utils.R"), local = TRUE)
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
write_tsv <- o2fr_write_tsv

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  run_dir <- argv$run_dir
  if (is.null(run_dir) || !nzchar(trimws(run_dir))) stop("Missing --run_dir.", call. = FALSE)
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  simulation_dir <- normalizePath(argv$simulation_dir, mustWork = TRUE)
  if (!file.exists(file.path(simulation_dir, "simulation_manifest.tsv"))) stop("Missing fit-result simulation manifest: ", simulation_dir, call. = FALSE)
  gate_path <- file.path(simulation_dir, "prediction_gate_metrics.tsv")
  if (!file.exists(gate_path)) stop("Missing prediction gate metrics: ", gate_path, call. = FALSE)
  out_dir <- normalizePath(argv$analysis_dir %||% argv$out_dir %||% file.path(run_dir, "extra_results"), mustWork = FALSE)
  allow_partial_seed_dirs <- as_bool(argv$allow_partial_seed_dirs, FALSE)
  near_thresh <- as_num(argv$near_thresh, 0.05)
  if (!is.finite(near_thresh) || near_thresh <= 0 || near_thresh >= 0.5) stop("near_thresh must be in (0, 0.5).")

  seed_dirs <- find_seed_dirs(run_dir)
  if (!length(seed_dirs)) stop("No valid seed directories found under: ", run_dir)
  existing_seed_summary_path <- file.path(out_dir, "seed_summary.tsv")
  if (!allow_partial_seed_dirs && file.exists(existing_seed_summary_path)) {
    old <- tryCatch(utils::read.delim(existing_seed_summary_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.data.frame(old) && "seed" %in% names(old) && length(unique(old$seed)) > length(seed_dirs)) {
      stop("Refusing to overwrite an existing all-seed summary with fewer seed directories; set --allow_partial_seed_dirs=TRUE if intentional.")
    }
  }
  gates <- utils::read.delim(gate_path, check.names = FALSE, stringsAsFactors = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  long_rows <- vector("list", length(seed_dirs))
  summary_records <- vector("list", length(seed_dirs))
  joint_invitro_long_rows <- vector("list", length(seed_dirs))
  joint_soft_rows <- vector("list", length(seed_dirs))

  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    seed <- basename(seed_dir)
    fit_summary_vals <- read_metric_map(file.path(seed_dir, "fit_summary.tsv"), "metric", "value")
    fit_summary_vals <- supplement_joint_invitro_metrics(fit_summary_vals, seed_dir)
    best_vals <- read_metric_map(file.path(seed_dir, "best_params.tsv"), "parameter", "value")
    gate <- gates[as.character(gates$seed) == seed, , drop = FALSE]
    pred_gate_metrics <- if (nrow(gate)) as.list(gate[1L, c("pred1000_2N", "pred1000_4N", "pred1000_both_gt44"), drop = FALSE]) else list(pred1000_2N = NA_real_, pred1000_4N = NA_real_, pred1000_both_gt44 = FALSE)
    param_table <- read_parameter_table_checked(find_parameter_table_for_seed(run_dir, seed_dir, fit_summary_vals))
    objective <- as_num(summary_metric_value(fit_summary_vals, "objective", summary_metric_value(fit_summary_vals, "objective_total", NA_real_)))
    long_df <- build_parameter_long_table(seed, objective, fit_summary_vals, best_vals, param_table, near_thresh)
    long_rows[[i]] <- long_df
    summary_records[[i]] <- build_seed_summary_record(seed, fit_summary_vals, best_vals, long_df, pred_gate_metrics)
    if (is_joint_fit_summary(fit_summary_vals)) {
      soft_path <- file.path(seed_dir, "joint_soft_coupling.tsv")
      if (file.exists(soft_path)) {
        soft <- tryCatch(utils::read.delim(soft_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.data.frame(soft) && nrow(soft)) { soft$seed <- seed; joint_soft_rows[[i]] <- soft }
      }
      invitro_path <- find_joint_invitro_parameter_table_for_seed(run_dir, seed_dir)
      if (!is.na(invitro_path) && nzchar(invitro_path)) {
        invitro_summary <- prepare_joint_invitro_fit_summary_values(fit_summary_vals)
        joint_invitro_long_rows[[i]] <- build_parameter_long_table(
          seed, as_num(summary_metric_value(invitro_summary, "objective", NA_real_)), invitro_summary,
          best_vals, read_parameter_table_checked(invitro_path), near_thresh
        )
      }
    }
  }

  parameter_long <- do.call(rbind, long_rows)
  joint_invitro_long_rows <- joint_invitro_long_rows[vapply(joint_invitro_long_rows, function(x) !is.null(x) && nrow(x), logical(1))]
  joint_invitro_parameter_long <- if (length(joint_invitro_long_rows)) do.call(rbind, joint_invitro_long_rows) else NULL
  joint_soft_rows <- joint_soft_rows[vapply(joint_soft_rows, function(x) !is.null(x) && nrow(x), logical(1))]
  joint_soft_coupling_all <- if (length(joint_soft_rows)) do.call(rbind, joint_soft_rows) else NULL
  seed_summary <- bind_records(summary_records)
  required_optional <- c(
    "fit_mode", "objective_total", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik",
    "growth_loglik_sum", "ploidy_loglik_sum", "flow_loglik_sum", "sigma_growth", "sigma_kary", "sigma_flow_ploidy",
    "n_growth", "n_ploidy_passages", "n_kary_cells", "n_flow_passages", "n_flow_samples",
    "objective_invivo", "objective_invitro", "objective_soft_coupling", "objective_constraints",
    "joint_weight_invivo", "joint_weight_invitro", "joint_invitro_growth_weight", "joint_invitro_ploidy_weight", "joint_invitro_flow_weight",
    "joint_soft_coupling_enabled", "joint_soft_coupling_params", "joint_soft_coupling_n_params", "joint_soft_coupling_sigma_default", "joint_soft_coupling_welsch_c",
    "n_cores_requested", "n_cores_used", "n_parameters", "n_invivo_scenarios", "itermax",
    "optimizer_deoptim_objective", "optimizer_local_objective", "optimizer_local_attempted", "optimizer_local_accepted",
    "optimizer_local_convergence", "optimizer_local_maxit", "optimizer_interrupted", "optimizer_iter_completed", "optimizer_iter_target", "deoptim_stop_reason"
  )
  for (col in required_optional) if (!(col %in% names(seed_summary))) seed_summary[[col]] <- NA
  seed_summary <- supplement_optimizer_fields_from_refinement_csv(seed_summary, run_dir)
  numeric_cols <- c(
    "objective", "objective_total", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik",
    "growth_loglik_sum", "ploidy_loglik_sum", "flow_loglik_sum", "sigma_growth", "sigma_kary", "sigma_flow_ploidy",
    "n_growth", "n_ploidy_passages", "n_kary_cells", "n_flow_passages", "n_flow_samples",
    "objective_invivo", "objective_invitro", "objective_soft_coupling", "objective_constraints",
    "joint_weight_invivo", "joint_weight_invitro", "joint_invitro_growth_weight", "joint_invitro_ploidy_weight", "joint_invitro_flow_weight",
    "joint_soft_coupling_n_params", "joint_soft_coupling_sigma_default", "joint_soft_coupling_welsch_c",
    "n_cores_requested", "n_cores_used", "n_parameters", "n_invivo_scenarios", "itermax",
    "optimizer_deoptim_objective", "optimizer_local_objective", "optimizer_local_convergence", "optimizer_local_maxit",
    "optimizer_iter_completed", "optimizer_iter_target", "objective_ploidy", "objective_burden",
    "objective_ploidy_neg2loglik_raw", "objective_burden_neg2loglik_raw", "boundary_penalty_active", "min_rel_dist_active",
    "boundary_penalty_active_excl_sigma_burden", "min_rel_dist_active_excl_sigma_burden", "pred1000_2N", "pred1000_4N"
  )
  for (col in numeric_cols) seed_summary[[col]] <- suppressWarnings(as.numeric(seed_summary[[col]]))
  for (col in c("joint_soft_coupling_enabled", "optimizer_local_attempted", "optimizer_local_accepted", "pred1000_both_gt44")) seed_summary[[col]] <- as.logical(seed_summary[[col]])
  for (col in c("fit_mode", "optimizer_interrupted", "deoptim_stop_reason")) seed_summary[[col]] <- as.character(seed_summary[[col]])

  seed_summary$objective_rank <- rank(seed_summary$objective, ties.method = "first", na.last = "keep")
  seed_summary$objective_ploidy_rank <- rank(seed_summary$objective_ploidy, ties.method = "first", na.last = "keep")
  seed_summary$objective_burden_rank <- rank(seed_summary$objective_burden, ties.method = "first", na.last = "keep")
  is_joint_run <- any(seed_summary$fit_mode == "fit_joint", na.rm = TRUE) || any(is.finite(seed_summary$objective_invivo) | is.finite(seed_summary$objective_invitro))
  is_invitro_run <- any(seed_summary$fit_mode == "fit_invitro", na.rm = TRUE) || any(is.finite(seed_summary$objective_total) | is.finite(seed_summary$growth_loglik) | is.finite(seed_summary$ploidy_loglik))
  is_invitro_only_run <- is_invitro_run && !is_joint_run
  boundary_order <- if (is_joint_run || is_invitro_run) {
    order(seed_summary$boundary_penalty_active, -seed_summary$min_rel_dist_active, seed_summary$objective, seed_summary$seed, na.last = TRUE)
  } else {
    order(seed_summary$boundary_penalty_active, -seed_summary$min_rel_dist_active, seed_summary$objective_burden, seed_summary$objective, seed_summary$seed, na.last = TRUE)
  }
  seed_summary$boundary_rank_active_support <- NA_integer_
  seed_summary$boundary_rank_active_support[boundary_order] <- seq_len(nrow(seed_summary))
  if (is_joint_run) {
    joint_order <- function(idx) idx[order(seed_summary$objective[idx], seed_order_key(seed_summary$seed[idx]), seed_summary$seed[idx], na.last = TRUE)]
    eligible <- which(!is.na(seed_summary$pred1000_both_gt44) & seed_summary$pred1000_both_gt44 & is.finite(seed_summary$objective))
    recommend_order <- if (length(eligible)) c(joint_order(eligible), joint_order(setdiff(seq_len(nrow(seed_summary)), eligible))) else joint_order(seq_len(nrow(seed_summary)))
    seed_summary$recommend_score_burden_ploidy_boundary <- NA_real_
    seed_summary$recommend_score_burden_ploidy_boundary[recommend_order] <- seq_along(recommend_order)
  } else if (is_invitro_run) {
    seed_summary$recommend_score_burden_ploidy_boundary <- with(seed_summary, objective_rank + 0.1 * boundary_rank_active_support)
    recommend_order <- order(seed_summary$objective, seed_summary$boundary_penalty_active, -seed_summary$min_rel_dist_active, seed_summary$seed, na.last = TRUE)
  } else {
    seed_summary$recommend_score_burden_ploidy_boundary <- with(seed_summary, objective_burden_rank + 0.2 * objective_ploidy_rank + 0.1 * boundary_rank_active_support)
    recommend_order <- order(seed_summary$objective_burden, seed_summary$objective_ploidy, seed_summary$boundary_penalty_active, -seed_summary$min_rel_dist_active, seed_summary$objective, seed_summary$seed, na.last = TRUE)
  }
  seed_summary$recommend_score_ploidy_burden_boundary <- seed_summary$recommend_score_burden_ploidy_boundary
  seed_summary$recommend_score_ploidy_boundary <- seed_summary$recommend_score_burden_ploidy_boundary
  for (col in c("recommend_rank_burden_ploidy_boundary_first", "recommend_rank_ploidy_burden_boundary_first", "recommend_rank_ploidy_boundary_first")) {
    seed_summary[[col]] <- NA_integer_; seed_summary[[col]][recommend_order] <- seq_len(nrow(seed_summary))
  }
  seed_summary$recommend_rank_ploidy_first <- seed_summary$recommend_rank_burden_ploidy_boundary_first
  seed_summary$forest_plot_rank_simple <- suppressWarnings(as.integer(seed_summary$objective_rank))
  seed_summary$forest_plot_rank_plus_ploidy_simple <- NA_integer_
  eligible_plot <- which(!is.na(seed_summary$pred1000_both_gt44) & seed_summary$pred1000_both_gt44)
  if (length(eligible_plot)) {
    ord <- eligible_plot[order(seed_summary$objective[eligible_plot], seed_summary$seed[eligible_plot], na.last = TRUE)]
    seed_summary$forest_plot_rank_plus_ploidy_simple[ord] <- seq_along(ord)
  }
  seed_summary <- seed_summary[order(seed_summary$objective, seed_summary$seed), , drop = FALSE]
  row.names(seed_summary) <- NULL
  convergence_summary <- build_convergence_summary(seed_summary)

  objective_cols <- c("seed", "objective")
  if (is_joint_run) objective_cols <- c(objective_cols, "objective_invivo", "objective_invitro")
  if (is_invitro_run) objective_cols <- c(objective_cols, "objective_total", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik")
  objective_cols <- c(objective_cols, intersect(c("objective_burden", "objective_ploidy"), names(seed_summary)))
  objective_simple <- seed_summary[, objective_cols, drop = FALSE]
  objective_simple$objective_rank <- suppressWarnings(as.integer(seed_summary$forest_plot_rank_simple))
  objective_simple$objective_rank_plus_ploidy <- suppressWarnings(as.integer(seed_summary$forest_plot_rank_plus_ploidy_simple))
  objective_simple <- objective_simple[, c("seed", "objective_rank", "objective_rank_plus_ploidy", setdiff(names(objective_simple), c("seed", "objective_rank", "objective_rank_plus_ploidy"))), drop = FALSE]
  objective_simple <- objective_simple[order(objective_simple$objective_rank, objective_simple$objective, objective_simple$seed, na.last = TRUE), , drop = FALSE]
  row.names(objective_simple) <- NULL

  outputs <- c("seed_summary.tsv", "convergence_summary.tsv", "parameter_boundary_long.tsv", "seed_objective_simple.tsv")
  write_tsv(seed_summary, file.path(out_dir, outputs[[1]]))
  write_tsv(convergence_summary, file.path(out_dir, outputs[[2]]))
  write_tsv(parameter_long, file.path(out_dir, outputs[[3]]))
  write_tsv(objective_simple, file.path(out_dir, outputs[[4]]))
  if (!is.null(joint_invitro_parameter_long) && nrow(joint_invitro_parameter_long)) { write_tsv(joint_invitro_parameter_long, file.path(out_dir, "invitro_parameter_boundary_long.tsv")); outputs <- c(outputs, "invitro_parameter_boundary_long.tsv") }
  if (!is.null(joint_soft_coupling_all) && nrow(joint_soft_coupling_all)) { write_tsv(joint_soft_coupling_all, file.path(out_dir, "joint_soft_coupling_all.tsv")); outputs <- c(outputs, "joint_soft_coupling_all.tsv") }
  if (is_joint_run) {
    cols <- intersect(c("seed", "objective_rank", "objective", "objective_invivo", "objective_invitro", "objective_soft_coupling", "objective_constraints", "joint_soft_coupling_enabled", "joint_soft_coupling_params", "joint_soft_coupling_n_params", "joint_soft_coupling_sigma_default", "joint_soft_coupling_welsch_c", "joint_weight_invivo", "joint_weight_invitro", "joint_invitro_growth_weight", "joint_invitro_ploidy_weight", "joint_invitro_flow_weight", "boundary_penalty_active", "min_rel_dist_active"), names(seed_summary))
    joint_simple <- seed_summary[, cols, drop = FALSE]
    joint_simple <- joint_simple[order(joint_simple$objective, joint_simple$seed, na.last = TRUE), , drop = FALSE]
    write_tsv(joint_simple, file.path(out_dir, "joint_objective_simple.tsv")); outputs <- c(outputs, "joint_objective_simple.tsv")
  }
  if (is_invitro_run) {
    invitro_summary <- if (is_joint_run) prepare_invitro_summary_for_plot(seed_summary, TRUE) else seed_summary
    cols <- intersect(c("seed", "objective_rank", "objective", "objective_total", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik", "sigma_growth", "sigma_kary", "sigma_flow_ploidy", "n_growth", "n_ploidy_passages", "n_kary_cells", "n_flow_passages", "n_flow_samples", "boundary_penalty_active", "min_rel_dist_active"), names(invitro_summary))
    invitro_simple <- invitro_summary[, cols, drop = FALSE]
    invitro_simple <- invitro_simple[order(invitro_simple$objective, invitro_simple$seed, na.last = TRUE), , drop = FALSE]
    write_tsv(invitro_simple, file.path(out_dir, "invitro_objective_simple.tsv")); outputs <- c(outputs, "invitro_objective_simple.tsv")
  }
  objective_long <- build_long_objective_table(seed_summary)
  write_tsv(objective_long, file.path(out_dir, "objective_components_long.tsv")); outputs <- c(outputs, "objective_components_long.tsv")
  run_mode <- data.frame(run_dir = run_dir, run_label = basename(run_dir), is_joint_run = is_joint_run, is_invitro_run = is_invitro_run, is_invitro_only_run = is_invitro_only_run, near_thresh = near_thresh, stringsAsFactors = FALSE)
  write_tsv(run_mode, file.path(out_dir, "fit_results_analysis_context.tsv")); outputs <- c(outputs, "fit_results_analysis_context.tsv")
  write_tsv(data.frame(stage = "analysis", file = outputs, stringsAsFactors = FALSE), file.path(out_dir, "analysis_manifest.tsv"))
  message("Wrote fit-results analysis contract: ", out_dir)
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
