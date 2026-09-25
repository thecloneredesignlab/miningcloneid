#!/usr/bin/env Rscript

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
parse_args <- o2fr_parse_args
`%||%` <- o2fr_null_coalesce
as_num <- o2fr_as_num

main <- function(argv = parse_args()) {
  simulation_dir <- normalizePath(argv$simulation_dir, mustWork = TRUE)
  manifest <- file.path(simulation_dir, "simulation_manifest.tsv")
  metrics_path <- file.path(simulation_dir, "invivo_long_ploidy_metrics.tsv")
  if (!file.exists(manifest) || !file.exists(metrics_path)) stop("Missing long-ploidy simulation contract: ", simulation_dir, call. = FALSE)
  threshold_n <- as_num(argv$threshold_N, 44)
  out_tsv <- normalizePath(argv$out_tsv %||% file.path(simulation_dir, "best_long_ploidy_gt2_seed.tsv"), mustWork = FALSE)
  best_dir_file <- normalizePath(argv$best_dir_file %||% file.path(simulation_dir, "best_long_ploidy_gt2_seed.dir"), mustWork = FALSE)
  metrics <- utils::read.delim(metrics_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("seed", "seed_dir", "objective", "long_term_weighted_mean_N", "long_term_status", "prediction_path")
  missing <- setdiff(required, names(metrics))
  if (length(missing)) stop("Long-ploidy metrics missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  res <- metrics
  res$eligible <- is.finite(res$objective) & res$long_term_status == "ok" & is.finite(res$long_term_weighted_mean_N) & res$long_term_weighted_mean_N > threshold_n
  res$selected <- FALSE
  res$selection_threshold_N <- threshold_n
  legacy_order <- c(
    "seed", "seed_dir", "objective", "long_term_day", "long_term_weighted_mean_N",
    "long_term_rows", "long_term_status", "eligible", "prediction_path", "selected",
    "selection_horizon", "selection_threshold_N", "selection_cohort", "selection_dose"
  )
  eligible <- res[res$eligible %in% TRUE, , drop = FALSE]
  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  if (!nrow(eligible)) {
    utils::write.table(res[, legacy_order, drop = FALSE], out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    stop("No eligible seed found with long-term weighted_mean_N > ", threshold_n, ". Summary written to: ", out_tsv, call. = FALSE)
  }
  eligible <- eligible[order(eligible$objective, -eligible$long_term_weighted_mean_N, eligible$seed), , drop = FALSE]
  best_seed <- eligible$seed[[1]]
  res$selected[res$seed == best_seed] <- TRUE
  res <- res[, legacy_order, drop = FALSE]
  utils::write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(res$seed_dir[res$selected][[1]], best_dir_file)
  manifest_out <- data.frame(stage = "analysis", file = c(basename(out_tsv), basename(best_dir_file)), stringsAsFactors = FALSE)
  utils::write.table(manifest_out, file.path(dirname(out_tsv), "analysis_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Selected in vivo seed: seed", best_seed)
  invisible(out_tsv)
}

if (sys.nframe() == 0L) main()
