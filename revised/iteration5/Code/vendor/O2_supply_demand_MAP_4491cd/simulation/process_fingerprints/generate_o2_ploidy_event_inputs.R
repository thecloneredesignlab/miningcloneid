#!/usr/bin/env Rscript

# Materialize fitted-seed metadata plus existing per-seed O2/ploidy trajectories
# for downstream event-coupling and medium-O2 analyses.  No classification,
# hypothesis test, plotting, or report generation occurs here.

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
source(file.path(.script_dir, "process_fingerprint_simulation_legacy_utils.R"), local = TRUE)
source(file.path(.script_dir, "process_fingerprint_simulation_utils.R"), local = TRUE)

oes_safe_read <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 1L) return(data.frame())
  tryCatch(o2ipa_read_tsv(path), error = function(e) data.frame())
}

oes_read_ploidy_curve <- function(seed_id, seed_dir, n_unit = 22) {
  sim_dir <- file.path(seed_dir, "simulation", "invivo")
  candidates <- file.path(sim_dir, c(
    "predict_ploidy_weighted_mean_0_1000day.tsv",
    "ploidy_weighted_mean_timecourse.tsv",
    "predict_ploidy_weighted_mean_0_300day.tsv",
    "predict_ploidy_0_1000day.tsv",
    "ploidy_timecourse.tsv"
  ))
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) return(data.frame())
  path <- hits[[1]]
  tab <- oes_safe_read(path)
  if (!nrow(tab) || !all(c("cohort", "day") %in% names(tab))) return(data.frame())
  vcol <- intersect(c("weighted_mean_N", "weighted_mean_endpoint", "ploidy_value", "weighted_mean_ploidy"), names(tab))
  if (!length(vcol)) {
    if (!all(c("N", "fraction") %in% names(tab))) return(data.frame())
    tab$N <- suppressWarnings(as.numeric(tab$N))
    tab$fraction <- suppressWarnings(as.numeric(tab$fraction))
    key <- paste(tab$cohort, tab$day, sep = "\r")
    tab <- do.call(rbind, lapply(split(seq_len(nrow(tab)), key), function(idx) {
      d <- tab[idx, , drop = FALSE]
      denom <- sum(d$fraction, na.rm = TRUE)
      data.frame(cohort = d$cohort[[1]], day = d$day[[1]], weighted_mean_N = if (denom > 0) sum(d$N * d$fraction, na.rm = TRUE) / denom else NA_real_)
    }))
    vcol <- "weighted_mean_N"
  }
  value <- suppressWarnings(as.numeric(tab[[vcol[[1]]]]))
  ploidy <- if (length(value) && stats::median(value, na.rm = TRUE) > 10) value / n_unit else value
  out <- data.frame(
    seed_id = seed_id, cohort = as.character(tab$cohort),
    day = suppressWarnings(as.numeric(tab$day)), ploidy = ploidy,
    source_file = normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE
  )
  out <- out[out$cohort %in% c("2N", "4N") & is.finite(out$day) & is.finite(out$ploidy), , drop = FALSE]
  if (!nrow(out)) return(out)
  result <- aggregate(out$ploidy, by = list(seed_id = out$seed_id, cohort = out$cohort, day = out$day, source_file = out$source_file), FUN = mean)
  names(result)[5] <- "ploidy"
  result
}

oes_read_o2_curve <- function(seed_id, seed_dir) {
  sim_dir <- file.path(seed_dir, "simulation", "invivo")
  candidates <- file.path(sim_dir, c(
    "predict_burden_0_1000day.tsv", "o2_lag_timecourse.tsv",
    "predict_burden_vs_o2.tsv", "predict_burden_0_300day.tsv"
  ))
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) return(data.frame())
  path <- hits[[1]]
  tab <- oes_safe_read(path)
  if (!nrow(tab) || !all(c("cohort", "day") %in% names(tab))) return(data.frame())
  o2_col <- intersect(c("pred_o2_pct", "o2_eff_pct", "o2_pct"), names(tab))
  target_col <- intersect(c("pred_o2_target_pct", "o2_target_pct"), names(tab))
  burden_col <- intersect(c(
    "pred_burden_cells", "pred_burden_total_cells", "pred_burden",
    "pred_burden_volume_mm3", "burden_mm3", "burden_total", "burden_value"
  ), names(tab))
  if (!length(o2_col)) return(data.frame())
  out <- data.frame(
    seed_id = seed_id, cohort = as.character(tab$cohort),
    day = suppressWarnings(as.numeric(tab$day)),
    o2_pct = suppressWarnings(as.numeric(tab[[o2_col[[1]]]])),
    o2_target_pct = if (length(target_col)) suppressWarnings(as.numeric(tab[[target_col[[1]]]])) else NA_real_,
    burden = if (length(burden_col)) suppressWarnings(as.numeric(tab[[burden_col[[1]]]])) else NA_real_,
    source_file = normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE
  )
  out <- out[out$cohort %in% c("2N", "4N") & is.finite(out$day) & is.finite(out$o2_pct), , drop = FALSE]
  if (!nrow(out)) return(out)
  rows <- lapply(split(seq_len(nrow(out)), paste(out$cohort, out$day, sep = "\r")), function(idx) {
    d <- out[idx, , drop = FALSE]
    mean_or_na <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
    data.frame(
      seed_id = seed_id, cohort = d$cohort[[1]], day = d$day[[1]],
      o2_pct = mean_or_na(d$o2_pct), o2_target_pct = mean_or_na(d$o2_target_pct),
      burden = mean_or_na(d$burden), source_file = d$source_file[[1]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

oes_read_labels <- function(path) {
  schema <- data.frame(
    seed_id = character(), trajectory_regime = character(), mode_label = character(),
    objective = numeric(), delta_objective = numeric(), stringsAsFactors = FALSE
  )
  if (!nzchar(path) || !file.exists(path)) return(schema)
  tab <- oes_safe_read(path)
  if (!nrow(tab) || !"seed_id" %in% names(tab)) return(schema)
  keep <- intersect(c(
    "seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective",
    "mean_o2_day0_pct", "mean_time_o2_near_floor", "mean_terminal_ploidy",
    "mean_late_drop_amplitude"
  ), names(tab))
  tab[, keep, drop = FALSE]
}

run_o2_ploidy_event_input_simulation <- function(argv = o2ipa_parse_args()) {
  run_dir <- o2ipa_as_chr(argv$run_dir)
  if (!nzchar(run_dir) || !dir.exists(run_dir)) stop("Missing or invalid --run_dir.")
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_dir <- o2ipa_as_chr(o2ipa_null_coalesce(argv$simulation_dir, argv$out_dir), file.path(run_dir, "simulation", "process_fingerprints", "o2_ploidy_events"))
  if (!nzchar(out_dir)) stop("Missing --simulation_dir (or --out_dir).")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  n_unit <- o2ipa_as_num(argv$n_unit, 22)
  max_seeds <- o2ipa_as_int(argv$max_seeds, NA_integer_)
  label_file <- o2ipa_as_chr(argv$label_file)

  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  manifest <- inputs$manifest
  if (!is.na(max_seeds) && max_seeds > 0L && nrow(manifest) > max_seeds) {
    manifest <- manifest[order(o2ipa_seed_number(manifest$seed_id)), , drop = FALSE]
    manifest <- manifest[seq_len(max_seeds), , drop = FALSE]
    inputs$params_long <- inputs$params_long[inputs$params_long$seed_id %in% manifest$seed_id, , drop = FALSE]
  }
  if (!nrow(manifest)) stop("No fitted seeds discovered in: ", run_dir)

  ploidy_rows <- list()
  o2_rows <- list()
  status_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    seed <- manifest$seed_id[[i]]
    seed_dir <- as.character(manifest$seed_dir[[i]])
    ploidy <- oes_read_ploidy_curve(seed, seed_dir, n_unit)
    oxygen <- oes_read_o2_curve(seed, seed_dir)
    if (nrow(ploidy)) ploidy_rows[[seed]] <- ploidy
    if (nrow(oxygen)) o2_rows[[seed]] <- oxygen
    status_rows[[i]] <- data.frame(
      seed_id = seed, ploidy_available = nrow(ploidy) > 0L,
      o2_available = nrow(oxygen) > 0L,
      ploidy_source_file = if (nrow(ploidy)) ploidy$source_file[[1]] else NA_character_,
      o2_source_file = if (nrow(oxygen)) oxygen$source_file[[1]] else NA_character_,
      stringsAsFactors = FALSE
    )
  }
  ploidy <- if (length(ploidy_rows)) do.call(rbind, ploidy_rows) else data.frame()
  oxygen <- if (length(o2_rows)) do.call(rbind, o2_rows) else data.frame()
  status <- do.call(rbind, status_rows)
  if (!nrow(ploidy) || !nrow(oxygen)) {
    stop("O2/ploidy event simulation requires materialized per-seed simulation/invivo O2 and ploidy trajectories.")
  }

  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  o2ipa_write_tsv(manifest, file.path(out_dir, "tables", "event_seed_manifest.tsv"))
  o2ipa_write_tsv(inputs$params_long, file.path(out_dir, "tables", "event_parameter_values_long.tsv"))
  o2ipa_write_tsv(ploidy, file.path(out_dir, "tables", "event_ploidy_timecourses.tsv"))
  o2ipa_write_tsv(oxygen, file.path(out_dir, "tables", "event_o2_timecourses.tsv"))
  o2ipa_write_tsv(status, file.path(out_dir, "tables", "event_input_status.tsv"))
  o2ipa_write_tsv(oes_read_labels(label_file), file.path(out_dir, "tables", "event_input_labels.tsv"))
  o2ipa_write_tsv(data.frame(
    argument = c("run_dir", "n_unit", "max_seeds", "label_file"),
    value = c(run_dir, n_unit, max_seeds, label_file), stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "event_simulation_arguments.tsv"))
  files <- file.path("tables", c(
    "event_seed_manifest.tsv", "event_parameter_values_long.tsv",
    "event_ploidy_timecourses.tsv", "event_o2_timecourses.tsv",
    "event_input_status.tsv", "event_input_labels.tsv", "event_simulation_arguments.tsv"
  ))
  pfs_write_manifest(out_dir, files)
  message("Completed O2/ploidy event input simulation: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_o2_ploidy_event_input_simulation()
