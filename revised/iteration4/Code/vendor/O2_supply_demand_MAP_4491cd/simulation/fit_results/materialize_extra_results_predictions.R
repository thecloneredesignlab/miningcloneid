#!/usr/bin/env Rscript

# Materialize concrete O2-model prediction values from completed fit outputs.
# This stage performs no statistical comparison, plotting, or report rendering.

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
parse_args <- o2fr_parse_args
`%||%` <- o2fr_null_coalesce
seed_order_key <- function(x) suppressWarnings(as.numeric(sub("^seed", "", basename(as.character(x)))))

find_seed_dirs <- function(run_dir) {
  dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^seed[0-9]+$", basename(dirs))]
  dirs[order(seed_order_key(dirs), basename(dirs))]
}

seed_prediction_path <- function(seed_dir, filename) {
  candidates <- c(
    file.path(seed_dir, "simulation", "invivo", filename),
    file.path(seed_dir, "viz", filename),
    file.path(seed_dir, "viz", "invivo", filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else NA_character_
}

read_prediction_tsv <- function(path) {
  if (!length(path) || is.na(path) || !file.exists(path)) return(NULL)
  tryCatch(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}
first_existing_col <- o2fr_first_existing_col

prediction_seed_day_mean <- function(tab, seed, value_col, value_name) {
  if (is.null(tab) || !nrow(tab) || is.na(value_col) || !all(c("cohort", "day", value_col) %in% names(tab))) return(data.frame())
  value <- suppressWarnings(as.numeric(tab[[value_col]]))
  day <- suppressWarnings(as.numeric(tab$day))
  cohort <- as.character(tab$cohort)
  keep <- cohort %in% c("2N", "4N") & is.finite(day) & is.finite(value)
  if (!any(keep)) return(data.frame())
  raw <- data.frame(seed = as.character(seed), cohort = cohort[keep], day = day[keep], value = value[keep], stringsAsFactors = FALSE)
  out <- stats::aggregate(value ~ seed + cohort + day, raw, function(x) mean(x[is.finite(x)], na.rm = TRUE))
  names(out)[names(out) == "value"] <- value_name
  out[order(seed_order_key(out$seed), out$seed, out$cohort, out$day), , drop = FALSE]
}

summarise_seed_prediction_ci <- function(seed_day_df, value_col) {
  if (!nrow(seed_day_df)) return(data.frame())
  seed_day_df[[value_col]] <- suppressWarnings(as.numeric(seed_day_df[[value_col]]))
  seed_day_df <- seed_day_df[is.finite(seed_day_df[[value_col]]) & !is.na(seed_day_df$cohort) & is.finite(seed_day_df$day), , drop = FALSE]
  parts <- split(seq_len(nrow(seed_day_df)), interaction(seed_day_df$cohort, seed_day_df$day, drop = TRUE, sep = "\r"))
  rows <- lapply(parts, function(idx) {
    sub <- seed_day_df[idx, , drop = FALSE]
    vals <- sub[[value_col]]
    n_seed <- length(vals)
    mean_value <- mean(vals, na.rm = TRUE)
    sd_value <- if (n_seed > 1L) stats::sd(vals, na.rm = TRUE) else 0
    se_value <- if (n_seed > 1L) sd_value / sqrt(n_seed) else 0
    data.frame(
      cohort = sub$cohort[[1]], day = sub$day[[1]], n_seed = n_seed,
      mean_value = mean_value, median_value = stats::median(vals, na.rm = TRUE),
      sd_value = sd_value, se_value = se_value,
      ci_low = mean_value - 1.96 * se_value, ci_high = mean_value + 1.96 * se_value,
      min_value = min(vals, na.rm = TRUE), max_value = max(vals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out[order(out$cohort, out$day), , drop = FALSE]
}

read_ploidy_predictions <- function(seed_dirs, horizon_tag) {
  rows <- lapply(seed_dirs, function(seed_dir) {
    tab <- read_prediction_tsv(seed_prediction_path(seed_dir, paste0("predict_ploidy_weighted_mean_", horizon_tag, ".tsv")))
    value_col <- if (is.null(tab)) NA_character_ else first_existing_col(tab, c("weighted_mean_endpoint", "weighted_mean_N", "weighted_mean_ploidy"))
    prediction_seed_day_mean(tab, basename(seed_dir), value_col, "ploidy_value")
  })
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

read_burden_predictions <- function(seed_dirs, horizon_tag) {
  rows <- lapply(seed_dirs, function(seed_dir) {
    tab <- read_prediction_tsv(seed_prediction_path(seed_dir, paste0("predict_burden_", horizon_tag, ".tsv")))
    value_col <- if (is.null(tab)) NA_character_ else first_existing_col(tab, c("pred_burden_volume_mm3", "pred_burden", "burden_total"))
    prediction_seed_day_mean(tab, basename(seed_dir), value_col, "burden_value")
  })
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

gate_metrics <- function(seed, ploidy, target_day = 1000, threshold = 44) {
  tab <- ploidy[as.character(ploidy$seed) == seed, , drop = FALSE]
  out <- data.frame(seed = seed, pred1000_2N = NA_real_, pred1000_4N = NA_real_, pred1000_both_gt44 = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab)) return(out)
  target <- tab[abs(tab$day - target_day) <= 1e-8, , drop = FALSE]
  if (!nrow(target)) {
    day_use <- max(tab$day[tab$day <= target_day], na.rm = TRUE)
    if (!is.finite(day_use)) day_use <- max(tab$day, na.rm = TRUE)
    target <- tab[abs(tab$day - day_use) <= 1e-8, , drop = FALSE]
  }
  values <- lapply(c("2N", "4N"), function(cohort) target$ploidy_value[target$cohort == cohort & is.finite(target$ploidy_value)])
  names(values) <- c("2N", "4N")
  if (length(values[["2N"]])) out$pred1000_2N <- min(values[["2N"]])
  if (length(values[["4N"]])) out$pred1000_4N <- min(values[["4N"]])
  out$pred1000_both_gt44 <- length(values[["2N"]]) > 0L && length(values[["4N"]]) > 0L && all(values[["2N"]] > threshold) && all(values[["4N"]] > threshold)
  out
}

write_tsv <- o2fr_write_tsv

main <- function(argv = parse_args()) {
  run_dir <- normalizePath(argv$run_dir, mustWork = TRUE)
  out_dir <- normalizePath(argv$simulation_dir %||% argv$out_dir %||% file.path(run_dir, "simulation", "fit_results", "extra_results"), mustWork = FALSE)
  horizon_tag <- argv$horizon_tag %||% "0_1000day"
  seed_dirs <- find_seed_dirs(run_dir)
  if (!length(seed_dirs)) stop("No valid seed directories found under: ", run_dir, call. = FALSE)
  ploidy <- read_ploidy_predictions(seed_dirs, horizon_tag)
  burden <- read_burden_predictions(seed_dirs, horizon_tag)
  ploidy_ci <- summarise_seed_prediction_ci(ploidy, "ploidy_value")
  burden_ci <- summarise_seed_prediction_ci(burden, "burden_value")
  gates <- do.call(rbind, lapply(basename(seed_dirs), gate_metrics, ploidy = ploidy))
  status <- data.frame(
    seed = basename(seed_dirs), seed_dir = normalizePath(seed_dirs, mustWork = TRUE),
    ploidy_available = basename(seed_dirs) %in% unique(as.character(ploidy$seed)),
    burden_available = basename(seed_dirs) %in% unique(as.character(burden$seed)), stringsAsFactors = FALSE
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- list(
    predict1000_ploidy_seed_day_mean.tsv = ploidy,
    predict1000_ploidy_mean_ci.tsv = ploidy_ci,
    predict1000_burden_total_seed_day_mean.tsv = burden,
    predict1000_burden_total_mean_ci.tsv = burden_ci,
    prediction_gate_metrics.tsv = gates,
    prediction_input_status.tsv = status
  )
  for (name in names(outputs)) write_tsv(outputs[[name]], file.path(out_dir, name))
  manifest <- data.frame(stage = "simulation", file = names(outputs), stringsAsFactors = FALSE)
  write_tsv(manifest, file.path(out_dir, "simulation_manifest.tsv"))
  message("Wrote fit-result prediction simulation contract: ", out_dir)
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
