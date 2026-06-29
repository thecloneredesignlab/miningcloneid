#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
WORKFLOW_DIR <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)
REPO_ROOT <- normalizePath(file.path(WORKFLOW_DIR, "..", "..", ".."), mustWork = FALSE)
FIXO2_SCRIPT_PATH <- file.path(WORKFLOW_DIR, "simulation", "fix_o2_simulation.R")
CURVE_UTILS_PATH <- normalizePath(file.path(SCRIPT_DIR, "curve_classification_utils.R"), mustWork = TRUE)

source(CURVE_UTILS_PATH)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[[1]], fixed = TRUE)
    val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x)) return(default)
  val <- tolower(as.character(x[[1]]))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(x[[1]] %||% default))
  if (length(val) && is.finite(val)) val else default
}

as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(x[[1]] %||% default))
  if (length(val) && is.finite(val)) val else default
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(as.character(x[[1]]))) return(default)
  vals <- suppressWarnings(as.numeric(trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

format_num <- function(x) {
  format(as.numeric(x), scientific = FALSE, trim = TRUE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

rbind_nonempty <- function(rows) {
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

seed_number <- function(seed_id) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

normalize_seed_ids <- function(seed_ids) {
  seed_ids <- as.character(seed_ids)
  ifelse(grepl("^seed", seed_ids), seed_ids, paste0("seed", seed_ids))
}

seed_file_stem <- function(seed_id) {
  sn <- seed_number(seed_id)
  if (is.finite(sn)) sprintf("seed%03d", sn) else as.character(seed_id)
}

num_key <- function(x) {
  vapply(x, function(xx) format(signif(as.numeric(xx), 12), scientific = FALSE, trim = TRUE), character(1))
}

num_slug <- function(x) {
  key <- num_key(x)
  key <- gsub("-", "minus", key, fixed = TRUE)
  key <- gsub("[^0-9A-Za-z]+", "p", key)
  key <- gsub("^p+|p+$", "", key)
  ifelse(nzchar(key), key, "NA")
}

default_fit_root <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
}

default_output_root <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "analysis", "dense-grid_initial-ploidy_trajectory")
}

analysis_paths <- function(output_root) {
  table_dir <- file.path(output_root, "tables")
  daily_dir <- file.path(table_dir, "daily_trajectories")
  figure_dir <- file.path(output_root, "figures")
  list(
    selected = file.path(table_dir, "fixed_o2_initial_ploidy_selected_time_curves.tsv"),
    delta = file.path(table_dir, "fixed_o2_initial_ploidy_delta_by_seed_o2_time.tsv"),
    class_by_seed_time = file.path(table_dir, "fixed_o2_initial_ploidy_curve_class_by_seed_time.tsv"),
    convergence = file.path(table_dir, "fixed_o2_initial_ploidy_convergence_summary.tsv"),
    run_arguments = file.path(table_dir, "analysis_run_arguments.tsv"),
    validation = file.path(table_dir, "validation.tsv"),
    daily_manifest = file.path(daily_dir, "manifest.tsv"),
    figures = figure_dir,
    figure_manifest = file.path(figure_dir, "figure_manifest.tsv")
  )
}

load_fixo2_env <- function(script_path = FIXO2_SCRIPT_PATH) {
  script_path <- normalizePath(path.expand(script_path), mustWork = TRUE)
  env <- new.env(parent = globalenv())
  env$commandArgs <- function(trailingOnly = FALSE) {
    if (isTRUE(trailingOnly)) character(0) else paste0("--file=", script_path)
  }
  old_error <- getOption("error")
  on.exit(options(error = old_error), add = TRUE)
  sys.source(script_path, envir = env, chdir = TRUE)
  options(error = old_error)
  required <- c(
    "cf2_fixed_matrix", "cf2_init_vector", "cf2_eigen_trajectory",
    "fixo2_normalize_state_matrix", "fixo2_eigen_trajectory_cached",
    "fixo2_expm_trajectory", "fixo2_trajectory_with_fallback", "fixo2_dominant_from_eig",
    "o2ipa_collect_seed_inputs", "o2ipa_params_wide", "o2pr_first_seed_cfg",
    "o2pr_run_params_from_vec", "o2ipa_source_model"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = TRUE)]
  if (length(missing)) stop("FixO2 helper environment is missing: ", paste(missing, collapse = ", "))
  env
}

fixo2_helper <- local({
  cache <- NULL

  function(name) {
    if (is.null(cache)) cache <<- load_fixo2_env()
    get(name, envir = cache, inherits = TRUE)
  }
})

normalize_state_matrix <- function(state_mat) {
  fixo2_helper("fixo2_normalize_state_matrix")(state_mat)
}

eigen_trajectory_cached <- function(eig, ngrid, init, time_grid, n_unit) {
  fixo2_helper("fixo2_eigen_trajectory_cached")(eig, ngrid, init, time_grid, n_unit)
}

expm_trajectory_daily <- function(M, ngrid, init, time_grid, n_unit) {
  fixo2_helper("fixo2_expm_trajectory")(M, ngrid, init, time_grid, n_unit)
}

trajectory_with_fallback <- function(M, eig, ngrid, init, time_grid, n_unit) {
  fixo2_helper("fixo2_trajectory_with_fallback")(M, eig, ngrid, init, time_grid, n_unit)
}

dominant_from_eig <- function(eig, ngrid, n_unit) {
  fixo2_helper("fixo2_dominant_from_eig")(eig, ngrid, n_unit)
}

initial_specs <- function(n_unit) {
  data.frame(
    initial_condition = c("init_2N", "init_4N"),
    initial_ploidy = c(2, 4),
    requested_initial_N = c(2, 4) * n_unit,
    stringsAsFactors = FALSE
  )
}

process_one_seed <- function(seed_id, param_mat, cfg, model_env, fixed_matrix_fn, init_vector_fn,
                             run_params_fn, o2_grid, time_grid, selected_times, daily_dir, force) {
  message("Processing ", seed_id)
  seed_num <- seed_number(seed_id)
  daily_path <- file.path(daily_dir, paste0(seed_file_stem(seed_id), ".tsv.gz"))
  expected_rows <- length(o2_grid) * 2L * length(time_grid)
  if (file.exists(daily_path) && !force) {
    stop("Daily trajectory file already exists; rerun with --force=true to overwrite: ", daily_path)
  }
  pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
  names(pvec) <- colnames(param_mat)
  run_params <- run_params_fn(pvec, cfg)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  init_specs <- initial_specs(n_unit)
  rows <- list()
  init_audit_rows <- list()
  row_i <- 0L
  audit_i <- 0L
  status_values <- character()

  for (O2 in o2_grid) {
    fm <- tryCatch(fixed_matrix_fn(model_env, cfg, run_params, O2), error = function(e) e)
    if (inherits(fm, "error")) {
      status_values <- c(status_values, paste0("matrix_error:", conditionMessage(fm)))
      next
    }
    eig <- tryCatch(eigen(fm$M, only.values = FALSE), error = function(e) e)
    if (inherits(eig, "error")) {
      status_values <- c(status_values, paste0("eigen_error:", conditionMessage(eig)))
      next
    }
    dom <- dominant_from_eig(eig, fm$ngrid, n_unit)
    for (j in seq_len(nrow(init_specs))) {
      init <- init_vector_fn(fm$ngrid, init_specs$requested_initial_N[[j]])
      audit_i <- audit_i + 1L
      init_audit_rows[[audit_i]] <- data.frame(
        seed_id = seed_id,
        seed_number = seed_num,
        initial_condition = init_specs$initial_condition[[j]],
        requested_initial_ploidy = init_specs$initial_ploidy[[j]],
        requested_initial_N = init_specs$requested_initial_N[[j]],
        used_initial_N = init$used_N,
        used_initial_ploidy = init$used_ploidy,
        O2_pct = O2,
        stringsAsFactors = FALSE
      )
      sim <- trajectory_with_fallback(fm$M, eig, fm$ngrid, init$vector, time_grid, n_unit)
      tr <- sim$trajectory
      status_values <- c(status_values, sim$status)
      if (!nrow(tr)) next
      tr$seed_id <- seed_id
      tr$seed_number <- seed_num
      tr$O2_pct <- O2
      tr$O2_key <- num_key(O2)
      tr$initial_condition <- init_specs$initial_condition[[j]]
      tr$initial_ploidy <- init_specs$initial_ploidy[[j]]
      tr$requested_initial_N <- init_specs$requested_initial_N[[j]]
      tr$used_initial_N <- init$used_N
      tr$used_initial_ploidy <- init$used_ploidy
      tr$status <- sim$status
      tr$trajectory_method <- sim$trajectory_method
      tr$dominant_mean_N <- dom$dominant_mean_N[[1L]]
      tr$dominant_mean_ploidy <- dom$dominant_mean_ploidy[[1L]]
      tr$dominant_fraction_N_le_25 <- dom$dominant_fraction_N_le_25[[1L]]
      tr$dominant_fraction_N_below_44 <- dom$dominant_fraction_N_below_44[[1L]]
      tr$dominant_fraction_N_ge_44 <- dom$dominant_fraction_N_ge_44[[1L]]
      tr$dominant_growth_rate <- dom$dominant_growth_rate[[1L]]
      tr$second_growth_rate <- dom$second_growth_rate[[1L]]
      tr$spectral_gap <- dom$spectral_gap[[1L]]
      tr$relative_spectral_gap <- dom$relative_spectral_gap[[1L]]
      tr$relax_time_days <- dom$relax_time_days[[1L]]
      tr$time_to_10x_days <- dom$time_to_10x_days[[1L]]
      tr$time_to_100x_days <- dom$time_to_100x_days[[1L]]
      tr$log10_advantage_1000d <- dom$log10_advantage_1000d[[1L]]
      tr$dominance_class <- dom$dominance_class[[1L]]
      row_i <- row_i + 1L
      rows[[row_i]] <- tr[, c(
        "seed_id", "seed_number", "O2_pct", "O2_key", "initial_condition", "initial_ploidy",
        "requested_initial_N", "used_initial_N", "used_initial_ploidy", "status", "day",
        "trajectory_method", "mean_N", "mean_ploidy", "fraction_N_le_25", "fraction_N_below_44",
        "fraction_N_ge_44", "fraction_N_ge_66", "fraction_N_ge_88", "dominant_mean_N",
        "dominant_mean_ploidy", "dominant_fraction_N_le_25", "dominant_fraction_N_below_44",
        "dominant_fraction_N_ge_44", "dominant_growth_rate", "second_growth_rate", "spectral_gap",
        "relative_spectral_gap", "relax_time_days", "time_to_10x_days", "time_to_100x_days",
        "log10_advantage_1000d", "dominance_class"
      )]
    }
  }

  daily <- rbind_nonempty(rows)
  if (nrow(daily)) {
    daily <- daily[order(daily$seed_number, daily$O2_pct, daily$initial_ploidy, daily$day), , drop = FALSE]
    write_tsv_gz(daily, daily_path)
  }
  selected <- if (nrow(daily)) daily[daily$day %in% selected_times, , drop = FALSE] else data.frame()
  selected <- selected[order(selected$seed_number, selected$O2_pct, selected$initial_ploidy, selected$day), , drop = FALSE]
  status <- if (!nrow(daily)) {
    "failed_no_rows"
  } else if (nrow(daily) != expected_rows) {
    "partial"
  } else if (all(unique(status_values) == "ok")) {
    "ok"
  } else {
    "ok_with_warnings"
  }
  init_audit <- rbind_nonempty(init_audit_rows)
  init_audit_unique <- if (nrow(init_audit)) {
    unique(init_audit[, c(
      "seed_id", "seed_number", "initial_condition", "requested_initial_ploidy",
      "requested_initial_N", "used_initial_N", "used_initial_ploidy"
    ), drop = FALSE])
  } else {
    data.frame()
  }
  manifest <- data.frame(
    seed_id = seed_id,
    seed_number = seed_num,
    daily_file = daily_path,
    expected_rows = expected_rows,
    observed_rows = nrow(daily),
    selected_rows = nrow(selected),
    file_exists = file.exists(daily_path),
    status = status,
    status_values = paste(sort(unique(status_values)), collapse = ";"),
    n_expm_fallback_trajectories = if (nrow(daily) && "trajectory_method" %in% names(daily)) {
      length(unique(paste(daily$O2_pct, daily$initial_condition, sep = "\r")[daily$trajectory_method == "expm_fallback"]))
    } else {
      0L
    },
    init_2N_used_initial_N = if (nrow(init_audit_unique)) init_audit_unique$used_initial_N[init_audit_unique$initial_condition == "init_2N"][[1L]] else NA_real_,
    init_4N_used_initial_N = if (nrow(init_audit_unique)) init_audit_unique$used_initial_N[init_audit_unique$initial_condition == "init_4N"][[1L]] else NA_real_,
    stringsAsFactors = FALSE
  )
  list(selected = selected, manifest = manifest, init_audit = init_audit_unique)
}

rename_value_cols <- function(d, suffix, keys) {
  value_cols <- setdiff(names(d), keys)
  names(d)[match(value_cols, names(d))] <- paste0(value_cols, suffix)
  d
}

build_delta_table <- function(selected) {
  keys <- c("seed_id", "seed_number", "O2_pct", "O2_key", "day")
  meta_cols <- c(
    keys, "dominant_mean_N", "dominant_mean_ploidy", "dominant_growth_rate",
    "second_growth_rate", "spectral_gap", "relative_spectral_gap", "relax_time_days",
    "time_to_10x_days", "time_to_100x_days", "log10_advantage_1000d", "dominance_class"
  )
  keep_value <- c(
    keys, "mean_N", "mean_ploidy", "fraction_N_le_25", "fraction_N_below_44",
    "fraction_N_ge_44", "fraction_N_ge_66", "fraction_N_ge_88", "used_initial_N",
    "used_initial_ploidy", "status"
  )
  d2 <- selected[selected$initial_condition == "init_2N", intersect(keep_value, names(selected)), drop = FALSE]
  d4 <- selected[selected$initial_condition == "init_4N", intersect(keep_value, names(selected)), drop = FALSE]
  d2 <- rename_value_cols(d2, "_2N", keys)
  d4 <- rename_value_cols(d4, "_4N", keys)
  out <- merge(d2, d4, by = keys, all = TRUE, sort = FALSE)
  meta <- unique(selected[, intersect(meta_cols, names(selected)), drop = FALSE])
  out <- merge(out, meta, by = keys, all.x = TRUE, sort = FALSE)
  out$delta_initial <- out$mean_ploidy_4N - out$mean_ploidy_2N
  out$abs_delta_initial <- abs(out$delta_initial)
  out$distance_to_dominant_ploidy_2N <- out$mean_ploidy_2N - out$dominant_mean_ploidy
  out$distance_to_dominant_ploidy_4N <- out$mean_ploidy_4N - out$dominant_mean_ploidy
  out$abs_distance_to_dominant_ploidy_2N <- abs(out$distance_to_dominant_ploidy_2N)
  out$abs_distance_to_dominant_ploidy_4N <- abs(out$distance_to_dominant_ploidy_4N)
  out[order(out$seed_number, out$O2_pct, out$day), , drop = FALSE]
}

classify_selected_curves <- function(selected, flat_range_threshold, step_epsilon_abs,
                                     step_epsilon_fraction, reverse_fraction_tolerance,
                                     plateau_min_points) {
  group_key <- paste(selected$seed_id, selected$day, selected$initial_condition, sep = "\r")
  groups <- split(selected, group_key)
  class_rows <- list()
  diff_rows <- list()
  for (i in seq_along(groups)) {
    curve <- groups[[i]]
    curve <- curve[order(curve$O2_pct), , drop = FALSE]
    res <- classify_o2_ploidy_curve(
      curve,
      value_col = "mean_ploidy",
      x_col = "O2_pct",
      id_col = "seed_id",
      flat_range_threshold = flat_range_threshold,
      step_epsilon_abs = step_epsilon_abs,
      step_epsilon_fraction = step_epsilon_fraction,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      plateau_min_points = plateau_min_points
    )
    sm <- res$summary
    class_rows[[i]] <- data.frame(
      seed_id = curve$seed_id[[1L]],
      seed_number = curve$seed_number[[1L]],
      day = curve$day[[1L]],
      initial_condition = curve$initial_condition[[1L]],
      initial_ploidy = curve$initial_ploidy[[1L]],
      n_o2 = nrow(curve),
      o2_min = min(curve$O2_pct, na.rm = TRUE),
      o2_max = max(curve$O2_pct, na.rm = TRUE),
      sm,
      min_mean_ploidy = min(curve$mean_ploidy, na.rm = TRUE),
      max_mean_ploidy = max(curve$mean_ploidy, na.rm = TRUE),
      min_spectral_gap = min(curve$spectral_gap, na.rm = TRUE),
      median_spectral_gap = stats::median(curve$spectral_gap, na.rm = TRUE),
      fraction_o2_gap_below_0p005 = mean(curve$spectral_gap < 0.005, na.rm = TRUE),
      fraction_o2_gap_below_0p01 = mean(curve$spectral_gap < 0.01, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    diffs <- res$differences
    diffs$day <- curve$day[[1L]]
    diffs$initial_condition <- curve$initial_condition[[1L]]
    diffs$initial_ploidy <- curve$initial_ploidy[[1L]]
    diff_rows[[i]] <- diffs[, c(
      "seed_id", "day", "initial_condition", "initial_ploidy", "O2_pct",
      "finite_difference_next", "local_slope_sign", "step_epsilon"
    )]
  }
  class_tab <- rbind_nonempty(class_rows)
  diff_tab <- rbind_nonempty(diff_rows)
  pair_rows <- lapply(split(class_tab, paste(class_tab$seed_id, class_tab$day, sep = "\r")), function(d) {
    c2 <- d$curve_class[d$initial_condition == "init_2N"]
    c4 <- d$curve_class[d$initial_condition == "init_4N"]
    s2 <- d$sign_sequence[d$initial_condition == "init_2N"]
    s4 <- d$sign_sequence[d$initial_condition == "init_4N"]
    data.frame(
      seed_id = d$seed_id[[1L]],
      day = d$day[[1L]],
      curve_class_2N = if (length(c2)) c2[[1L]] else NA_character_,
      curve_class_4N = if (length(c4)) c4[[1L]] else NA_character_,
      curve_class_consistent = if (length(c2) && length(c4)) identical(c2[[1L]], c4[[1L]]) else NA,
      sign_sequence_2N = if (length(s2)) s2[[1L]] else NA_character_,
      sign_sequence_4N = if (length(s4)) s4[[1L]] else NA_character_,
      sign_sequence_consistent = if (length(s2) && length(s4)) identical(s2[[1L]], s4[[1L]]) else NA,
      stringsAsFactors = FALSE
    )
  })
  pair_tab <- rbind_nonempty(pair_rows)
  class_tab <- merge(class_tab, pair_tab, by = c("seed_id", "day"), all.x = TRUE, sort = FALSE)
  class_tab <- class_tab[order(class_tab$seed_number, class_tab$day, class_tab$initial_ploidy), , drop = FALSE]
  diff_tab <- diff_tab[order(seed_number(diff_tab$seed_id), diff_tab$day, diff_tab$initial_ploidy, diff_tab$O2_pct), , drop = FALSE]
  list(class_table = class_tab, differences = diff_tab)
}

attach_curve_differences <- function(selected, differences) {
  if (!nrow(selected) || !nrow(differences)) return(selected)
  out <- merge(
    selected,
    differences,
    by = c("seed_id", "day", "initial_condition", "initial_ploidy", "O2_pct"),
    all.x = TRUE,
    sort = FALSE
  )
  out[order(out$seed_number, out$O2_pct, out$initial_ploidy, out$day), , drop = FALSE]
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

summarize_convergence <- function(delta, time_end) {
  terminal <- delta[delta$day == time_end, , drop = FALSE]
  if (!nrow(terminal)) return(data.frame())
  finite_gap <- is.finite(terminal$spectral_gap)
  gap_rank <- rep(NA_real_, nrow(terminal))
  gap_rank[finite_gap] <- rank(terminal$spectral_gap[finite_gap], ties.method = "average") / sum(finite_gap)
  delta_rank <- rep(NA_real_, nrow(terminal))
  finite_delta <- is.finite(terminal$abs_delta_initial)
  delta_rank[finite_delta] <- rank(-terminal$abs_delta_initial[finite_delta], ties.method = "average") / sum(finite_delta)
  terminal$spectral_gap_rank_fraction <- gap_rank
  terminal$abs_delta_top_rank_fraction <- delta_rank
  terminal$small_spectral_gap_bottom_10pct <- is.finite(gap_rank) & gap_rank <= 0.10
  terminal$large_delta_top_10pct <- is.finite(delta_rank) & delta_rank <= 0.10
  terminal$spectral_gap_decile <- cut(gap_rank, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)

  summarize_group <- function(d, scope, group_label) {
    data.frame(
      scope = scope,
      group = group_label,
      terminal_day = time_end,
      n_seed_o2 = nrow(d),
      n_seed = length(unique(d$seed_id)),
      n_o2 = length(unique(d$O2_pct)),
      spectral_gap_median = stats::median(d$spectral_gap, na.rm = TRUE),
      spectral_gap_min = suppressWarnings(min(d$spectral_gap, na.rm = TRUE)),
      spectral_gap_q25 = as.numeric(stats::quantile(d$spectral_gap, 0.25, na.rm = TRUE, names = FALSE)),
      spectral_gap_q75 = as.numeric(stats::quantile(d$spectral_gap, 0.75, na.rm = TRUE, names = FALSE)),
      abs_delta_initial_median = stats::median(d$abs_delta_initial, na.rm = TRUE),
      abs_delta_initial_q95 = as.numeric(stats::quantile(d$abs_delta_initial, 0.95, na.rm = TRUE, names = FALSE)),
      abs_delta_initial_max = suppressWarnings(max(d$abs_delta_initial, na.rm = TRUE)),
      abs_distance_to_dominant_ploidy_2N_median = stats::median(d$abs_distance_to_dominant_ploidy_2N, na.rm = TRUE),
      abs_distance_to_dominant_ploidy_4N_median = stats::median(d$abs_distance_to_dominant_ploidy_4N, na.rm = TRUE),
      fraction_2N_abs_distance_le_0p01 = mean(d$abs_distance_to_dominant_ploidy_2N <= 0.01, na.rm = TRUE),
      fraction_4N_abs_distance_le_0p01 = mean(d$abs_distance_to_dominant_ploidy_4N <= 0.01, na.rm = TRUE),
      fraction_2N_abs_distance_le_0p05 = mean(d$abs_distance_to_dominant_ploidy_2N <= 0.05, na.rm = TRUE),
      fraction_4N_abs_distance_le_0p05 = mean(d$abs_distance_to_dominant_ploidy_4N <= 0.05, na.rm = TRUE),
      fraction_abs_delta_le_0p01 = mean(d$abs_delta_initial <= 0.01, na.rm = TRUE),
      fraction_abs_delta_le_0p05 = mean(d$abs_delta_initial <= 0.05, na.rm = TRUE),
      fraction_small_gap_bottom_10pct = mean(d$small_spectral_gap_bottom_10pct, na.rm = TRUE),
      fraction_large_delta_top_10pct = mean(d$large_delta_top_10pct, na.rm = TRUE),
      fraction_large_delta_in_small_gap = if (any(d$large_delta_top_10pct %in% TRUE, na.rm = TRUE)) {
        mean(d$small_spectral_gap_bottom_10pct[d$large_delta_top_10pct %in% TRUE], na.rm = TRUE)
      } else {
        NA_real_
      },
      spearman_abs_delta_vs_spectral_gap = safe_cor(d$abs_delta_initial, d$spectral_gap, "spearman"),
      stringsAsFactors = FALSE
    )
  }

  rows <- list(summarize_group(terminal, "global", "all_seed_o2"))
  deciles <- split(terminal, terminal$spectral_gap_decile)
  if (length(deciles)) {
    rows <- c(rows, lapply(names(deciles), function(nm) summarize_group(deciles[[nm]], "spectral_gap_decile", nm)))
  }
  out <- rbind_nonempty(rows)
  out[order(match(out$scope, c("global", "spectral_gap_decile")), out$group), , drop = FALSE]
}

save_plot_pair <- function(name, figure_dir, expr, width = 8, height = 5) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(figure_dir, paste0(name, ".pdf"))
  png_path <- file.path(figure_dir, paste0(name, ".png"))
  expr <- substitute(expr)
  grDevices::pdf(pdf_path, width = width, height = height, onefile = FALSE)
  tryCatch(eval(expr, envir = parent.frame()), finally = grDevices::dev.off())
  grDevices::png(png_path, width = width, height = height, units = "in", res = 170)
  tryCatch(eval(expr, envir = parent.frame()), finally = grDevices::dev.off())
  c(pdf = pdf_path, png = png_path)
}

init_color_map <- function() {
  c(init_2N = "#0072B2", init_4N = "#D55E00")
}

plot_median_curves <- function(selected, figure_dir) {
  d <- selected[is.finite(selected$mean_ploidy), c("day", "initial_condition", "O2_pct", "mean_ploidy"), drop = FALSE]
  med <- aggregate(
    d$mean_ploidy,
    by = list(day = d$day, initial_condition = d$initial_condition, O2_pct = d$O2_pct),
    FUN = stats::median,
    na.rm = TRUE
  )
  names(med)[4] <- "median_mean_ploidy"
  plot_days <- sort(unique(med$day))
  cols <- init_color_map()
  save_plot_pair("median_mean_ploidy_curves_by_selected_time", figure_dir, {
    op <- par(mfrow = c(4, 4), mar = c(3.2, 3.5, 2.0, 0.8), oma = c(1, 1, 0, 0))
    on.exit(par(op), add = TRUE)
    yr <- range(med$median_mean_ploidy, na.rm = TRUE)
    for (day_i in plot_days) {
      z <- med[med$day == day_i, , drop = FALSE]
      plot(range(z$O2_pct), yr, type = "n", xlab = "O2 (%)", ylab = "Mean ploidy",
           main = paste0("Day ", day_i))
      for (ic in names(cols)) {
        ss <- z[z$initial_condition == ic, , drop = FALSE]
        ss <- ss[order(ss$O2_pct), , drop = FALSE]
        lines(ss$O2_pct, ss$median_mean_ploidy, col = cols[[ic]], lwd = 2)
      }
      abline(h = c(1, 2, 4), col = "grey85", lty = 3)
    }
    if (length(plot_days) < 16L) {
      for (i in seq_len(16L - length(plot_days))) plot.new()
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("bottom", legend = c("2N start", "4N start"), col = cols, lwd = 2, horiz = TRUE, bty = "n")
  }, width = 12, height = 9)
}

plot_delta_heatmap <- function(delta, figure_dir, terminal_day) {
  d <- delta[delta$day == terminal_day & is.finite(delta$delta_initial), , drop = FALSE]
  if (!nrow(d)) return(c(pdf = NA_character_, png = NA_character_))
  o2 <- sort(unique(d$O2_pct))
  score <- aggregate(d$abs_delta_initial, by = list(seed_id = d$seed_id), FUN = max, na.rm = TRUE)
  names(score)[2] <- "max_abs_delta"
  score$seed_number <- seed_number(score$seed_id)
  order_seed <- score$seed_id[order(score$max_abs_delta, score$seed_number)]
  mat <- matrix(NA_real_, nrow = length(order_seed), ncol = length(o2), dimnames = list(order_seed, num_key(o2)))
  for (i in seq_len(nrow(d))) {
    mat[d$seed_id[[i]], num_key(d$O2_pct[[i]])] <- d$delta_initial[[i]]
  }
  zlim <- as.numeric(stats::quantile(abs(mat), 0.99, na.rm = TRUE, names = FALSE))
  if (!is.finite(zlim) || zlim <= 0) zlim <- max(abs(mat), na.rm = TRUE)
  save_plot_pair("day1000_delta_initial_heatmap", figure_dir, {
    op <- par(mar = c(4.5, 4.5, 2.2, 1.2))
    on.exit(par(op), add = TRUE)
    image(
      x = o2,
      y = seq_along(order_seed),
      z = t(mat[rev(seq_along(order_seed)), , drop = FALSE]),
      zlim = c(-zlim, zlim),
      col = grDevices::hcl.colors(101, "Blue-Red 3"),
      xlab = "Fixed O2 (%)",
      ylab = "Seeds ordered by max |delta|",
      main = paste0("Day ", terminal_day, ": 4N-start minus 2N-start mean ploidy")
    )
    box()
    mtext(paste0("Color clipped at +/-", signif(zlim, 3)), side = 3, line = 0.2, cex = 0.8)
  }, width = 8.5, height = 7)
}

plot_delta_vs_gap <- function(delta, figure_dir, terminal_day) {
  d <- delta[delta$day == terminal_day & is.finite(delta$abs_delta_initial) & is.finite(delta$spectral_gap), , drop = FALSE]
  if (!nrow(d)) return(c(pdf = NA_character_, png = NA_character_))
  x <- pmax(d$spectral_gap, .Machine$double.eps)
  y <- pmax(d$abs_delta_initial, 1e-12)
  gap_q10 <- as.numeric(stats::quantile(d$spectral_gap, 0.10, na.rm = TRUE, names = FALSE))
  delta_q90 <- as.numeric(stats::quantile(d$abs_delta_initial, 0.90, na.rm = TRUE, names = FALSE))
  hi <- d$spectral_gap <= gap_q10 & d$abs_delta_initial >= delta_q90
  save_plot_pair("day1000_abs_delta_vs_spectral_gap", figure_dir, {
    op <- par(mar = c(4.8, 4.8, 2.0, 1))
    on.exit(par(op), add = TRUE)
    plot(x, y, log = "xy", pch = 16, cex = 0.45, col = grDevices::adjustcolor("black", alpha.f = 0.18),
         xlab = "Spectral gap", ylab = "|4N-start - 2N-start mean ploidy|",
         main = paste0("Day ", terminal_day, " initial-ploidy dependence vs spectral gap"))
    points(x[hi], y[hi], pch = 16, cex = 0.55, col = grDevices::adjustcolor("#D55E00", alpha.f = 0.8))
    abline(v = gap_q10, lty = 2, col = "#0072B2", lwd = 1.5)
    abline(h = max(delta_q90, 1e-12), lty = 2, col = "#D55E00", lwd = 1.5)
  }, width = 7.5, height = 5.8)
}

plot_curve_class_consistency <- function(class_table, figure_dir) {
  d <- unique(class_table[, c("seed_id", "day", "curve_class_consistent"), drop = FALSE])
  d <- d[!is.na(d$curve_class_consistent), , drop = FALSE]
  frac <- aggregate(as.numeric(d$curve_class_consistent), by = list(day = d$day), FUN = mean)
  names(frac)[2] <- "fraction_consistent"
  frac <- frac[order(frac$day), , drop = FALSE]
  save_plot_pair("curve_class_consistency_by_time", figure_dir, {
    op <- par(mar = c(4.5, 4.5, 1.5, 1))
    on.exit(par(op), add = TRUE)
    plot(frac$day, frac$fraction_consistent, type = "b", pch = 16, lwd = 2,
         ylim = c(0, 1), xlab = "Day", ylab = "Fraction of seeds with same curve class",
         main = "2N-start vs 4N-start O2-curve class consistency")
    grid(col = "grey85")
  }, width = 7.5, height = 5.2)
}

plot_convergence_boxplot <- function(delta, figure_dir, terminal_day) {
  d <- delta[delta$day == terminal_day, , drop = FALSE]
  vals <- data.frame(
    metric = rep(c("|2N - attractor|", "|4N - attractor|", "|4N - 2N|"), each = nrow(d)),
    value = pmax(c(
      d$abs_distance_to_dominant_ploidy_2N,
      d$abs_distance_to_dominant_ploidy_4N,
      d$abs_delta_initial
    ), 1e-12),
    stringsAsFactors = FALSE
  )
  vals <- vals[is.finite(vals$value), , drop = FALSE]
  save_plot_pair("day1000_convergence_distance_boxplot", figure_dir, {
    op <- par(mar = c(5.5, 4.8, 1.5, 1))
    on.exit(par(op), add = TRUE)
    boxplot(value ~ metric, data = vals, log = "y", las = 2, outline = FALSE,
            ylab = "Absolute mean-ploidy difference", xlab = "",
            col = c("#56B4E9", "#E69F00", "#999999"),
            main = paste0("Day ", terminal_day, " convergence and initial-ploidy difference"))
    stripchart(value ~ metric, data = vals[sample.int(nrow(vals), min(4000L, nrow(vals))), ],
               vertical = TRUE, method = "jitter", pch = 16, cex = 0.18,
               col = grDevices::adjustcolor("black", alpha.f = 0.08), add = TRUE)
  }, width = 7.5, height = 5.8)
}

plot_fallback_counts <- function(daily_manifest, figure_dir) {
  if (!"n_expm_fallback_trajectories" %in% names(daily_manifest)) {
    return(c(pdf = NA_character_, png = NA_character_))
  }
  d <- daily_manifest
  d <- d[order(d$seed_number), , drop = FALSE]
  save_plot_pair("expm_fallback_trajectory_counts_by_seed", figure_dir, {
    op <- par(mar = c(4.5, 4.5, 1.5, 1))
    on.exit(par(op), add = TRUE)
    plot(d$seed_number, d$n_expm_fallback_trajectories, type = "h", lwd = 1.2,
         xlab = "Seed", ylab = "Number of O2 x initial-ploidy trajectories",
         main = "Matrix::expm fallback use for numerically singular eigen trajectories")
    points(d$seed_number[d$n_expm_fallback_trajectories > 0],
           d$n_expm_fallback_trajectories[d$n_expm_fallback_trajectories > 0],
           pch = 16, cex = 0.55, col = "#D55E00")
    grid(col = "grey88")
  }, width = 8, height = 4.8)
}

curve_palette <- function(classes) {
  base <- c(
    approximately_flat = "#999999",
    monotone_increasing = "#0072B2",
    monotone_decreasing = "#D55E00",
    single_transition_increase_then_plateau = "#56B4E9",
    single_transition_decrease_then_plateau = "#E69F00",
    u_shaped = "#009E73",
    inverted_u_shaped = "#CC79A7",
    complex_nonmonotone = "#000000",
    insufficient_data = "#666666"
  )
  missing <- setdiff(classes, names(base))
  if (length(missing)) {
    extra <- grDevices::rainbow(length(missing), s = 0.65, v = 0.75)
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[classes]
}

finite_time_mode_label <- function(mean_ploidy) {
  out <- ifelse(is.finite(mean_ploidy) & mean_ploidy >= 2, "mode1", "mode2")
  out[!is.finite(mean_ploidy)] <- NA_character_
  out
}

monotonicity_class_curve_summary <- function(curves, by_seed) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class")], by = "seed_id", all.x = TRUE, sort = FALSE)
  rows <- lapply(split(d, list(d$curve_class, d$O2_pct), drop = TRUE), function(x) {
    data.frame(
      curve_class = x$curve_class[[1L]],
      O2_pct = x$O2_pct[[1L]],
      n_seed = length(unique(x$seed_id)),
      median_mean_ploidy = stats::median(x$mean_ploidy, na.rm = TRUE),
      q25_mean_ploidy = as.numeric(stats::quantile(x$mean_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
      q75_mean_ploidy = as.numeric(stats::quantile(x$mean_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
      median_spectral_gap = stats::median(x$spectral_gap, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- rbind_nonempty(rows)
  out[order(out$curve_class, out$O2_pct), , drop = FALSE]
}

monotonicity_representative_seeds <- function(curves, by_seed) {
  rows <- lapply(split(by_seed, by_seed$curve_class), function(seed_tab) {
    cls <- seed_tab$curve_class[[1L]]
    class_curves <- curves[curves$seed_id %in% seed_tab$seed_id, , drop = FALSE]
    med <- aggregate(class_curves$mean_ploidy, by = list(O2_pct = class_curves$O2_pct), FUN = stats::median, na.rm = TRUE)
    names(med)[2] <- "class_median"
    scores <- lapply(split(class_curves, class_curves$seed_id), function(curve) {
      d <- merge(curve[, c("seed_id", "O2_pct", "mean_ploidy")], med, by = "O2_pct", all.x = TRUE)
      data.frame(
        seed_id = curve$seed_id[[1L]],
        median_curve_rmse = sqrt(mean((d$mean_ploidy - d$class_median)^2, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    })
    score <- rbind_nonempty(scores)
    hit <- score[which.min(score$median_curve_rmse), , drop = FALSE]
    data.frame(
      curve_class = cls,
      representative_seed_id = hit$seed_id,
      representative_seed_number = seed_number(hit$seed_id),
      median_curve_rmse = hit$median_curve_rmse,
      class_n_seed = nrow(seed_tab),
      stringsAsFactors = FALSE
    )
  })
  out <- rbind_nonempty(rows)
  out[order(out$curve_class), , drop = FALSE]
}

plot_initial_all_curves_by_class <- function(curves, by_seed, figure_dir, title_suffix) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class")], by = "seed_id", all.x = TRUE, sort = FALSE)
  classes <- sort(unique(d$curve_class))
  cols <- curve_palette(classes)
  save_plot_pair("fixed_o2_initial_ploidy_mean_ploidy_all_seed_curves_by_class", figure_dir, {
    op <- par(mar = c(4.5, 4.5, 2.0, 1))
    on.exit(par(op), add = TRUE)
    plot(range(d$O2_pct), range(d$mean_ploidy, na.rm = TRUE), type = "n",
         xlab = "Fixed O2 (%)", ylab = "Mean ploidy",
         main = paste("All seed curves by class", title_suffix))
    for (seed in unique(d$seed_id[order(seed_number(d$seed_id))])) {
      z <- d[d$seed_id == seed, , drop = FALSE]
      z <- z[order(z$O2_pct), , drop = FALSE]
      cls <- z$curve_class[[1L]]
      lines(z$O2_pct, z$mean_ploidy, col = grDevices::adjustcolor(cols[[cls]], alpha.f = 0.18), lwd = 0.7)
    }
    legend("topright", legend = classes, col = cols, lwd = 2, bty = "n", cex = 0.7)
  }, width = 9, height = 5.5)
}

plot_initial_class_iqr <- function(summary, figure_dir, title_suffix) {
  classes <- sort(unique(summary$curve_class))
  cols <- curve_palette(classes)
  save_plot_pair("fixed_o2_initial_ploidy_mean_ploidy_median_iqr_by_class", figure_dir, {
    op <- par(mar = c(4.5, 4.5, 2.0, 1))
    on.exit(par(op), add = TRUE)
    plot(range(summary$O2_pct), range(c(summary$q25_mean_ploidy, summary$q75_mean_ploidy), na.rm = TRUE),
         type = "n", xlab = "Fixed O2 (%)", ylab = "Mean ploidy",
         main = paste("Median/IQR by class", title_suffix))
    for (cls in classes) {
      z <- summary[summary$curve_class == cls, , drop = FALSE]
      z <- z[order(z$O2_pct), , drop = FALSE]
      polygon(
        c(z$O2_pct, rev(z$O2_pct)),
        c(z$q25_mean_ploidy, rev(z$q75_mean_ploidy)),
        col = grDevices::adjustcolor(cols[[cls]], alpha.f = 0.20),
        border = NA
      )
      lines(z$O2_pct, z$median_mean_ploidy, col = cols[[cls]], lwd = 2)
    }
    legend("topright", legend = classes, col = cols, lwd = 2, bty = "n", cex = 0.7)
  }, width = 9, height = 5.5)
}

plot_initial_heatmap <- function(curves, by_seed, figure_dir, metric, file_name, zlab, transform = identity) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class", "seed_number")], by = "seed_id", all.x = TRUE, sort = FALSE)
  order_seed <- by_seed$seed_id[order(by_seed$curve_class, by_seed$seed_number)]
  o2 <- sort(unique(d$O2_pct))
  mat <- matrix(NA_real_, nrow = length(order_seed), ncol = length(o2), dimnames = list(order_seed, num_key(o2)))
  for (i in seq_len(nrow(d))) {
    mat[d$seed_id[[i]], num_key(d$O2_pct[[i]])] <- suppressWarnings(as.numeric(d[[metric]][[i]]))
  }
  mat <- transform(mat)
  save_plot_pair(file_name, figure_dir, {
    op <- par(mar = c(4.5, 4.5, 1.5, 5))
    on.exit(par(op), add = TRUE)
    z <- t(mat[nrow(mat):1L, , drop = FALSE])
    image(
      x = o2,
      y = seq_len(nrow(mat)),
      z = z,
      xlab = "Fixed O2 (%)",
      ylab = "Seeds ordered by class",
      col = grDevices::hcl.colors(80, "Viridis")
    )
    box()
    mtext(zlab, side = 4, line = 3)
  }, width = 8, height = 7)
}

plot_initial_mode_heatmap <- function(curves, by_seed, figure_dir) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class", "seed_number")], by = "seed_id", all.x = TRUE, sort = FALSE)
  d$finite_time_mode_label <- finite_time_mode_label(d$mean_ploidy)
  order_seed <- by_seed$seed_id[order(by_seed$curve_class, by_seed$seed_number)]
  o2 <- sort(unique(d$O2_pct))
  vals <- c(mode2 = 1, mode1 = 2)
  mat <- matrix(NA_real_, nrow = length(order_seed), ncol = length(o2), dimnames = list(order_seed, num_key(o2)))
  for (i in seq_len(nrow(d))) {
    mat[d$seed_id[[i]], num_key(d$O2_pct[[i]])] <- vals[[as.character(d$finite_time_mode_label[[i]])]]
  }
  save_plot_pair("fixed_o2_initial_ploidy_mode_label_heatmap_ordered_by_class", figure_dir, {
    op <- par(mar = c(4.5, 4.5, 1.5, 4))
    on.exit(par(op), add = TRUE)
    z <- t(mat[nrow(mat):1L, , drop = FALSE])
    image(
      x = o2,
      y = seq_len(nrow(mat)),
      z = z,
      xlab = "Fixed O2 (%)",
      ylab = "Seeds ordered by class",
      col = c("#D55E00", "#0072B2"),
      breaks = c(0.5, 1.5, 2.5)
    )
    box()
    legend("right", inset = -0.18, xpd = TRUE, legend = c("mean ploidy < 2", "mean ploidy >= 2"),
           fill = c("#D55E00", "#0072B2"), bty = "n", cex = 0.8)
  }, width = 8, height = 7)
}

plot_initial_gap_scatter <- function(by_seed, figure_dir, title_suffix) {
  classes <- sort(unique(by_seed$curve_class))
  cols <- curve_palette(classes)
  save_plot_pair("fixed_o2_initial_ploidy_min_spectral_gap_vs_ploidy_range", figure_dir, {
    op <- par(mar = c(4.5, 4.5, 2.0, 1))
    on.exit(par(op), add = TRUE)
    plot(by_seed$ploidy_range, by_seed$min_spectral_gap, log = "y",
         pch = 19, col = grDevices::adjustcolor(cols[by_seed$curve_class], alpha.f = 0.75),
         xlab = "Mean-ploidy range across O2", ylab = "Minimum spectral gap",
         main = paste("Range vs min spectral gap", title_suffix))
    abline(h = c(0.005, 0.01), lty = 2, col = "grey40")
    legend("topright", legend = classes, col = cols, pch = 19, bty = "n", cex = 0.7)
  }, width = 8, height = 5.5)
}

plot_initial_representatives <- function(curves, reps, figure_dir, title_suffix) {
  reps <- reps[order(reps$curve_class), , drop = FALSE]
  if (!nrow(reps)) return(c(pdf = NA_character_, png = NA_character_))
  save_plot_pair("fixed_o2_initial_ploidy_representative_curves_by_class", figure_dir, {
    op <- par(mfrow = c(nrow(reps), 2), mar = c(3, 4, 2, 1), oma = c(1, 0, 1.2, 0))
    on.exit(par(op), add = TRUE)
    for (i in seq_len(nrow(reps))) {
      seed <- reps$representative_seed_id[[i]]
      cls <- reps$curve_class[[i]]
      z <- curves[curves$seed_id == seed, , drop = FALSE]
      z <- z[order(z$O2_pct), , drop = FALSE]
      plot(z$O2_pct, z$mean_ploidy, type = "l", lwd = 2,
           xlab = "O2 (%)", ylab = "Mean ploidy", main = paste(cls, seed))
      abline(h = c(1.5, 2.0), lty = 2, col = "grey60")
      plot(z$O2_pct, z$spectral_gap, type = "l", lwd = 2, log = "y",
           xlab = "O2 (%)", ylab = "Spectral gap", main = paste(cls, seed))
      abline(h = c(0.005, 0.01), lty = 2, col = "grey60")
    }
    mtext(paste("Representative finite-time curves", title_suffix), outer = TRUE, cex = 1)
  }, width = 10, height = max(6, 2.1 * nrow(reps)))
}

write_monotonicity_like_figures <- function(selected, class_table, output_root) {
  paths <- analysis_paths(output_root)
  root_dir <- file.path(paths$figures, "monotonicity_like")
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
  keep_cols <- c("seed_id", "seed_number", "day", "initial_condition", "initial_ploidy", "O2_pct",
                 "mean_ploidy", "spectral_gap", "dominant_mean_ploidy")
  selected <- selected[, intersect(keep_cols, names(selected)), drop = FALSE]
  groups <- unique(class_table[, c("day", "initial_condition", "initial_ploidy"), drop = FALSE])
  groups <- groups[order(groups$day, groups$initial_ploidy), , drop = FALSE]
  manifests <- list()
  k <- 0L
  for (i in seq_len(nrow(groups))) {
    day_i <- groups$day[[i]]
    init_i <- groups$initial_condition[[i]]
    init_ploidy_i <- groups$initial_ploidy[[i]]
    curves <- selected[selected$day == day_i & selected$initial_condition == init_i, , drop = FALSE]
    by_seed <- class_table[class_table$day == day_i & class_table$initial_condition == init_i, , drop = FALSE]
    if (!nrow(curves) || !nrow(by_seed)) next
    curves <- curves[order(curves$seed_number, curves$O2_pct), , drop = FALSE]
    by_seed <- by_seed[order(by_seed$curve_class, by_seed$seed_number), , drop = FALSE]
    summary <- monotonicity_class_curve_summary(curves, by_seed)
    reps <- monotonicity_representative_seeds(curves, by_seed)
    group_dir <- file.path(root_dir, paste0("day_", sprintf("%04d", as.integer(day_i))), init_i)
    title_suffix <- paste0("(day ", day_i, ", ", init_i, ")")
    fig_map <- list(
      all_seed_curves_by_class = plot_initial_all_curves_by_class(curves, by_seed, group_dir, title_suffix),
      median_iqr_by_class = plot_initial_class_iqr(summary, group_dir, title_suffix),
      mean_ploidy_heatmap_ordered_by_class = plot_initial_heatmap(
        curves, by_seed, group_dir, "mean_ploidy",
        "fixed_o2_initial_ploidy_mean_ploidy_heatmap_ordered_by_class",
        "Mean ploidy"
      ),
      mode_label_heatmap_ordered_by_class = plot_initial_mode_heatmap(curves, by_seed, group_dir),
      spectral_gap_heatmap_ordered_by_class = plot_initial_heatmap(
        curves, by_seed, group_dir, "spectral_gap",
        "fixed_o2_initial_ploidy_spectral_gap_heatmap_ordered_by_class",
        "log10 spectral gap",
        transform = function(x) log10(pmax(x, .Machine$double.eps))
      ),
      min_spectral_gap_vs_ploidy_range = plot_initial_gap_scatter(by_seed, group_dir, title_suffix),
      representative_curves_by_class = plot_initial_representatives(curves, reps, group_dir, title_suffix)
    )
    for (nm in names(fig_map)) {
      k <- k + 1L
      manifests[[k]] <- data.frame(
        figure_id = paste("monotonicity_like", paste0("day_", day_i), init_i, nm, sep = "_"),
        day = day_i,
        initial_condition = init_i,
        initial_ploidy = init_ploidy_i,
        pdf = unname(fig_map[[nm]][["pdf"]]),
        png = unname(fig_map[[nm]][["png"]]),
        description = paste("Monotonicity-classification-style finite-time figure:", nm),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- rbind_nonempty(manifests)
  if (nrow(out)) {
    write_tsv(out, file.path(root_dir, "monotonicity_like_figure_manifest.tsv"))
  }
  out
}

time_grid_overlay_days <- function() {
  c(100, 200, 300, 500, 700, 1000)
}

overlay_base_line_style <- function(n_seed) {
  if (n_seed >= 200) return(list(alpha = 0.10, lwd = 0.30))
  if (n_seed >= 100) return(list(alpha = 0.12, lwd = 0.32))
  if (n_seed >= 50) return(list(alpha = 0.18, lwd = 0.38))
  if (n_seed >= 25) return(list(alpha = 0.25, lwd = 0.50))
  if (n_seed >= 10) return(list(alpha = 0.35, lwd = 0.65))
  if (n_seed >= 5) return(list(alpha = 0.48, lwd = 0.85))
  list(alpha = 0.65, lwd = 1.05)
}

overlay_curve_overlap_density <- function(z, y_col, y_range, transform = identity, n_y_bins = 24L) {
  y <- transform(z[[y_col]])
  y_range <- transform(y_range)
  ok <- is.finite(z$O2_pct) & is.finite(y)
  if (!any(ok) || !all(is.finite(y_range))) return(NA_real_)
  y_min <- min(y_range, na.rm = TRUE)
  y_max <- max(y_range, na.rm = TRUE)
  if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) return(NA_real_)
  breaks <- seq(y_min, y_max, length.out = n_y_bins + 1L)
  x_key <- sprintf("%.8f", z$O2_pct[ok])
  y_split <- split(y[ok], x_key)
  density <- vapply(y_split, function(vals) {
    bins <- cut(vals, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    bins <- bins[is.finite(bins)]
    n_bins <- length(unique(bins))
    if (!n_bins) return(NA_real_)
    length(vals) / n_bins
  }, numeric(1L))

  stats::median(density[is.finite(density)], na.rm = TRUE)
}

overlay_line_style <- function(z, y_col, y_range, transform = identity) {
  n_seed <- length(unique(z$seed_id))
  style <- overlay_base_line_style(n_seed)
  density <- overlay_curve_overlap_density(z, y_col, y_range, transform = transform)
  if (!is.finite(density) || density <= 0) return(style)

  density_target <- 12
  visibility_boost <- max(1, sqrt(density_target / density))
  style$alpha <- min(0.75, style$alpha * visibility_boost)
  style$lwd <- min(1.25, style$lwd * sqrt(visibility_boost))
  style
}

plot_empty_time_grid_panel <- function(panel_title, x_range, y_range, y_log = FALSE, panel_bg = NA_character_) {
  if (isTRUE(y_log)) {
    plot(x_range, y_range, type = "n", log = "y", xlab = "", ylab = "", main = panel_title,
         axes = FALSE, cex.main = 0.78)
  } else {
    plot(x_range, y_range, type = "n", xlab = "", ylab = "", main = panel_title,
         axes = FALSE, cex.main = 0.78)
  }
  if (!is.na(panel_bg)) {
    rect(x_range[[1L]], y_range[[1L]], x_range[[2L]], y_range[[2L]], col = panel_bg, border = NA)
  }
  box(col = "grey80")
  y_mid <- if (isTRUE(y_log)) exp(mean(log(y_range))) else mean(y_range)
  text(mean(x_range), y_mid, labels = "n=0", col = "grey45", cex = 0.8)
}

plot_time_grid_seed_overlays_by_class <- function(selected,
                                                  class_table,
                                                  output_root,
                                                  initial_condition,
                                                  metric = c("mean_ploidy", "spectral_gap")) {
  metric <- match.arg(metric)
  paths <- analysis_paths(output_root)
  figure_dir <- file.path(paths$figures, "all_seed_curves_by_class_time_grid")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  days <- time_grid_overlay_days()
  selected_cols <- c("seed_id", "day", "initial_condition", "O2_pct", metric)
  curves <- selected[
    selected$day %in% days & selected$initial_condition == initial_condition,
    intersect(selected_cols, names(selected)),
    drop = FALSE
  ]
  by_seed <- class_table[
    class_table$day %in% days & class_table$initial_condition == initial_condition,
    c("seed_id", "seed_number", "day", "initial_condition", "initial_ploidy", "curve_class"),
    drop = FALSE
  ]
  if (!nrow(curves) || !nrow(by_seed)) {
    return(c(pdf = NA_character_, png = NA_character_))
  }

  d <- merge(
    curves,
    by_seed,
    by = c("seed_id", "day", "initial_condition"),
    all = FALSE,
    sort = FALSE
  )
  d <- d[is.finite(d[[metric]]) & is.finite(d$O2_pct), , drop = FALSE]
  if (!nrow(d)) {
    return(c(pdf = NA_character_, png = NA_character_))
  }

  classes <- sort(unique(by_seed$curve_class))
  classes <- classes[!is.na(classes)]
  cols <- curve_palette(classes)
  init_ploidy <- unique(by_seed$initial_ploidy)
  init_ploidy <- init_ploidy[is.finite(init_ploidy)]
  init_ploidy <- if (length(init_ploidy)) init_ploidy[[1L]] else NA_real_
  x_range <- range(d$O2_pct, na.rm = TRUE)

  if (identical(metric, "spectral_gap")) {
    y_range <- range(d[[metric]][is.finite(d[[metric]]) & d[[metric]] > 0], c(0.005, 0.01), na.rm = TRUE)
    y_range <- pmax(y_range, .Machine$double.eps)
    y_label <- "Spectral gap"
    file_stub <- paste0("fixed_o2_initial_ploidy_spectral_gap_by_class_time_grid_", initial_condition)
    description_metric <- "spectral gap"
    y_log <- TRUE
    y_transform <- function(x) log10(pmax(x, .Machine$double.eps))
  } else {
    y_range <- range(d[[metric]][is.finite(d[[metric]])], c(1.5, 2.0), na.rm = TRUE)
    y_label <- "Finite-time mean ploidy"
    file_stub <- paste0("fixed_o2_initial_ploidy_mean_ploidy_by_class_time_grid_", initial_condition)
    description_metric <- "finite-time mean ploidy"
    y_log <- FALSE
    y_transform <- identity
  }

  if (!all(is.finite(y_range)) || y_range[[1L]] == y_range[[2L]]) {
    y_range <- if (isTRUE(y_log)) c(1e-6, 1) else c(1, 4)
  }

  save_plot_pair(file_stub, figure_dir, {
    op <- par(
      mfrow = c(length(classes), length(days)),
      mar = c(1.7, 2.1, 1.7, 0.45),
      oma = c(3.2, 3.4, 2.6, 0.3),
      mgp = c(1.8, 0.45, 0),
      tcl = -0.22
    )
    on.exit(par(op), add = TRUE)
    for (cls in classes) {
      for (day_i in days) {
        z <- d[d$curve_class == cls & d$day == day_i, , drop = FALSE]
        n_seed <- length(unique(z$seed_id))
        panel_title <- paste0(cls, "\nDay ", day_i, "  n=", n_seed)
        if (!n_seed) {
          plot_empty_time_grid_panel(panel_title, x_range, y_range, y_log = y_log)
          next
        }

        z <- z[order(z$seed_number, z$O2_pct), , drop = FALSE]
        seed_ids <- unique(z$seed_id[order(z$seed_number)])
        style <- overlay_line_style(z, metric, y_range, transform = y_transform)
        line_col <- grDevices::adjustcolor(cols[[cls]], alpha.f = style$alpha)
        if (isTRUE(y_log)) {
          plot(x_range, y_range, type = "n", log = "y", xlab = "", ylab = "",
               main = panel_title, cex.main = 0.78, cex.axis = 0.65)
          abline(h = c(0.005, 0.01), lty = 2, col = "grey72")
        } else {
          plot(x_range, y_range, type = "n", xlab = "", ylab = "",
               main = panel_title, cex.main = 0.78, cex.axis = 0.65)
          abline(h = c(1.5, 2.0), lty = 2, col = "grey72")
        }
        for (seed in seed_ids) {
          zz <- z[z$seed_id == seed, , drop = FALSE]
          zz <- zz[order(zz$O2_pct), , drop = FALSE]
          yy <- zz[[metric]]
          if (isTRUE(y_log)) yy <- pmax(yy, .Machine$double.eps)
          lines(zz$O2_pct, yy, col = line_col, lwd = style$lwd)
        }
        box()
      }
    }
    mtext("Fixed O2 (%)", side = 1, outer = TRUE, line = 1.8)
    mtext(y_label, side = 2, outer = TRUE, line = 2.1)
    mtext(paste0(initial_condition, " start: ", description_metric, " curves by class and selected time"),
          side = 3, outer = TRUE, line = 0.7, cex = 1.05)
  }, width = 2.15 * length(days), height = max(6.5, 1.85 * length(classes)))
}

write_time_grid_overlay_figures <- function(selected, class_table, output_root) {
  manifests <- list()
  k <- 0L
  for (initial_condition in c("init_2N", "init_4N")) {
    init_rows <- class_table[class_table$initial_condition == initial_condition, , drop = FALSE]
    init_ploidy <- unique(init_rows$initial_ploidy)
    init_ploidy <- init_ploidy[is.finite(init_ploidy)]
    init_ploidy <- if (length(init_ploidy)) init_ploidy[[1L]] else NA_real_
    for (metric in c("mean_ploidy", "spectral_gap")) {
      fig <- plot_time_grid_seed_overlays_by_class(
        selected = selected,
        class_table = class_table,
        output_root = output_root,
        initial_condition = initial_condition,
        metric = metric
      )
      k <- k + 1L
      manifests[[k]] <- data.frame(
        figure_id = paste("class_time_grid", initial_condition, metric, sep = "_"),
        day = NA_real_,
        initial_condition = initial_condition,
        initial_ploidy = init_ploidy,
        pdf = unname(fig[["pdf"]]),
        png = unname(fig[["png"]]),
        description = paste(
          "Rows are curve classifications and columns are days 100, 200, 300, 500, 700, and 1000;",
          metric,
          "curves use density-adaptive alpha/line width."
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  rbind_nonempty(manifests)
}

initial_ploidy_panel_label <- function(initial_condition) {
  ifelse(initial_condition == "init_2N", "2N", "4N")
}

plot_time_grid_initial_ploidy_comparison_by_class <- function(selected,
                                                              class_table,
                                                              output_root,
                                                              metric = c("mean_ploidy", "spectral_gap")) {
  metric <- match.arg(metric)
  paths <- analysis_paths(output_root)
  figure_dir <- file.path(paths$figures, "all_seed_curves_by_class_time_grid_initial_ploidy_comparison")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  days <- time_grid_overlay_days()
  initials <- c("init_2N", "init_4N")
  selected_cols <- c("seed_id", "day", "initial_condition", "O2_pct", metric)
  curves <- selected[
    selected$day %in% days & selected$initial_condition %in% initials,
    intersect(selected_cols, names(selected)),
    drop = FALSE
  ]
  by_seed <- class_table[
    class_table$day %in% days & class_table$initial_condition %in% initials,
    c("seed_id", "seed_number", "day", "initial_condition", "initial_ploidy", "curve_class"),
    drop = FALSE
  ]
  if (!nrow(curves) || !nrow(by_seed)) {
    return(c(pdf = NA_character_, png = NA_character_))
  }

  d <- merge(
    curves,
    by_seed,
    by = c("seed_id", "day", "initial_condition"),
    all = FALSE,
    sort = FALSE
  )
  d <- d[is.finite(d[[metric]]) & is.finite(d$O2_pct), , drop = FALSE]
  if (!nrow(d)) {
    return(c(pdf = NA_character_, png = NA_character_))
  }

  classes <- sort(unique(by_seed$curve_class))
  classes <- classes[!is.na(classes)]
  cols <- curve_palette(classes)
  x_range <- range(d$O2_pct, na.rm = TRUE)
  panel_bg_4n <- grDevices::adjustcolor("#D55E00", alpha.f = 0.08)

  if (identical(metric, "spectral_gap")) {
    global_y_range <- range(d[[metric]][is.finite(d[[metric]]) & d[[metric]] > 0], c(0.005, 0.01), na.rm = TRUE)
    global_y_range <- pmax(global_y_range, .Machine$double.eps)
    y_label <- "Spectral gap"
    file_stub <- "fixed_o2_initial_ploidy_spectral_gap_by_class_time_grid_2N_vs_4N"
    description_metric <- "spectral gap"
    y_log <- TRUE
    y_transform <- function(x) log10(pmax(x, .Machine$double.eps))
    ref_lines <- c(0.005, 0.01)
  } else {
    global_y_range <- range(d[[metric]][is.finite(d[[metric]])], c(1.5, 2.0), na.rm = TRUE)
    y_label <- "Finite-time mean ploidy"
    file_stub <- "fixed_o2_initial_ploidy_mean_ploidy_by_class_time_grid_2N_vs_4N"
    description_metric <- "finite-time mean ploidy"
    y_log <- FALSE
    y_transform <- identity
    ref_lines <- c(1.5, 2.0)
  }
  if (!all(is.finite(global_y_range)) || global_y_range[[1L]] == global_y_range[[2L]]) {
    global_y_range <- if (isTRUE(y_log)) c(1e-6, 1) else c(1, 4)
  }

  row_y_range <- setNames(vector("list", length(classes)), classes)
  for (cls in classes) {
    zz <- d[d$curve_class == cls, , drop = FALSE]
    if (identical(metric, "spectral_gap")) {
      yr <- range(zz[[metric]][is.finite(zz[[metric]]) & zz[[metric]] > 0], ref_lines, na.rm = TRUE)
      yr <- pmax(yr, .Machine$double.eps)
    } else {
      yr <- range(zz[[metric]][is.finite(zz[[metric]])], ref_lines, na.rm = TRUE)
    }
    if (!all(is.finite(yr)) || yr[[1L]] == yr[[2L]]) yr <- global_y_range
    row_y_range[[cls]] <- yr
  }

  save_plot_pair(file_stub, figure_dir, {
    op <- par(
      mfrow = c(length(classes), length(days) * length(initials)),
      mar = c(1.65, 1.75, 1.8, 0.35),
      oma = c(3.2, 3.4, 2.8, 0.3),
      mgp = c(1.6, 0.42, 0),
      tcl = -0.20
    )
    on.exit(par(op), add = TRUE)
    for (cls in classes) {
      y_range <- row_y_range[[cls]]
      for (day_i in days) {
        for (init_i in initials) {
          init_label <- initial_ploidy_panel_label(init_i)
          z <- d[d$curve_class == cls & d$day == day_i & d$initial_condition == init_i, , drop = FALSE]
          n_seed <- length(unique(z$seed_id))
          panel_title <- paste0(
            if (day_i == days[[1L]] && init_i == initials[[1L]]) paste0(cls, "\n") else "",
            "Day ", day_i, " / ", init_label, "\n",
            "n=", n_seed
          )
          panel_bg <- if (init_i == "init_4N") panel_bg_4n else NA_character_
          if (!n_seed) {
            plot_empty_time_grid_panel(panel_title, x_range, y_range, y_log = y_log, panel_bg = panel_bg)
            next
          }

          z <- z[order(z$seed_number, z$O2_pct), , drop = FALSE]
          seed_ids <- unique(z$seed_id[order(z$seed_number)])
          style <- overlay_line_style(z, metric, y_range, transform = y_transform)
          line_col <- grDevices::adjustcolor(cols[[cls]], alpha.f = style$alpha)
          if (isTRUE(y_log)) {
            plot(x_range, y_range, type = "n", log = "y", xlab = "", ylab = "",
                 main = panel_title, cex.main = 0.58, cex.axis = 0.52)
          } else {
            plot(x_range, y_range, type = "n", xlab = "", ylab = "",
                 main = panel_title, cex.main = 0.58, cex.axis = 0.52)
          }
          if (!is.na(panel_bg)) {
            rect(x_range[[1L]], y_range[[1L]], x_range[[2L]], y_range[[2L]], col = panel_bg, border = NA)
          }
          abline(h = ref_lines, lty = 2, col = "grey72")
          for (seed in seed_ids) {
            zz <- z[z$seed_id == seed, , drop = FALSE]
            zz <- zz[order(zz$O2_pct), , drop = FALSE]
            yy <- zz[[metric]]
            if (isTRUE(y_log)) yy <- pmax(yy, .Machine$double.eps)
            lines(zz$O2_pct, yy, col = line_col, lwd = style$lwd)
          }
          box()
        }
      }
    }
    mtext("Fixed O2 (%)", side = 1, outer = TRUE, line = 1.8)
    mtext(y_label, side = 2, outer = TRUE, line = 2.1)
    mtext(paste0("2N vs 4N starts: ", description_metric, " curves by class and selected time"),
          side = 3, outer = TRUE, line = 0.7, cex = 1.05)
  }, width = 2.15 * length(days) * length(initials), height = max(6.5, 1.85 * length(classes)))
}

write_time_grid_initial_ploidy_comparison_figures <- function(selected, class_table, output_root) {
  manifests <- list()
  for (metric in c("mean_ploidy", "spectral_gap")) {
    fig <- plot_time_grid_initial_ploidy_comparison_by_class(
      selected = selected,
      class_table = class_table,
      output_root = output_root,
      metric = metric
    )
    manifests[[metric]] <- data.frame(
      figure_id = paste("class_time_grid_initial_ploidy_comparison", metric, sep = "_"),
      day = NA_real_,
      initial_condition = "init_2N_vs_init_4N",
      initial_ploidy = NA_real_,
      pdf = unname(fig[["pdf"]]),
      png = unname(fig[["png"]]),
      description = paste(
        "Rows are curve classifications; columns are days 100, 200, 300, 500, 700, and 1000,",
        "with adjacent 2N and 4N panels per day; each row shares the y scale and 4N panels use a transparent orange background;",
        metric,
        "curves use density-adaptive alpha/line width."
      ),
      stringsAsFactors = FALSE
    )
  }
  rbind_nonempty(manifests)
}

write_figures <- function(selected, delta, class_table, convergence_summary, daily_manifest, output_root) {
  paths <- analysis_paths(output_root)
  figure_dir <- paths$figures
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  terminal_day <- if (nrow(convergence_summary) && "terminal_day" %in% names(convergence_summary)) {
    suppressWarnings(as.numeric(convergence_summary$terminal_day[convergence_summary$scope == "global"][[1L]]))
  } else {
    suppressWarnings(max(delta$day, na.rm = TRUE))
  }
  figures <- list(
    median_curves = plot_median_curves(selected, figure_dir),
    delta_heatmap = plot_delta_heatmap(delta, figure_dir, terminal_day),
    delta_vs_gap = plot_delta_vs_gap(delta, figure_dir, terminal_day),
    class_consistency = plot_curve_class_consistency(class_table, figure_dir),
    convergence_boxplot = plot_convergence_boxplot(delta, figure_dir, terminal_day),
    fallback_counts = plot_fallback_counts(daily_manifest, figure_dir)
  )
  descriptions <- c(
    median_curves = "Median O2 -> mean ploidy curves by initial ploidy across selected days.",
    delta_heatmap = "Day-1000 heatmap of mean_ploidy_4N - mean_ploidy_2N by seed and fixed O2.",
    delta_vs_gap = "Day-1000 initial-ploidy dependence versus spectral gap.",
    class_consistency = "Fraction of seeds where 2N-start and 4N-start O2-curve classifications agree.",
    convergence_boxplot = "Day-1000 distance to dominant attractor and residual initial-ploidy difference.",
    fallback_counts = "Per-seed count of Matrix::expm fallback trajectories used for numerically singular eigen trajectories."
  )
  manifest <- do.call(rbind, lapply(names(figures), function(nm) {
    data.frame(
      figure_id = nm,
      day = NA_real_,
      initial_condition = NA_character_,
      initial_ploidy = NA_real_,
      pdf = unname(figures[[nm]][["pdf"]]),
      png = unname(figures[[nm]][["png"]]),
      description = descriptions[[nm]],
      stringsAsFactors = FALSE
    )
  }))
  monotonicity_manifest <- write_monotonicity_like_figures(selected, class_table, output_root)
  if (nrow(monotonicity_manifest)) {
    manifest <- rbind(manifest, monotonicity_manifest[, names(manifest), drop = FALSE])
  }
  time_grid_manifest <- write_time_grid_overlay_figures(selected, class_table, output_root)
  if (nrow(time_grid_manifest)) {
    manifest <- rbind(manifest, time_grid_manifest[, names(manifest), drop = FALSE])
  }
  comparison_manifest <- write_time_grid_initial_ploidy_comparison_figures(selected, class_table, output_root)
  if (nrow(comparison_manifest)) {
    manifest <- rbind(manifest, comparison_manifest[, names(manifest), drop = FALSE])
  }
  write_tsv(manifest, paths$figure_manifest)
  manifest
}

write_figures_from_tables <- function(output_root) {
  paths <- analysis_paths(output_root)
  needed <- c(paths$selected, paths$delta, paths$class_by_seed_time, paths$convergence, paths$daily_manifest)
  missing <- needed[!file.exists(needed)]
  if (length(missing)) stop("Cannot plot; missing table(s): ", paste(missing, collapse = ", "))
  selected <- read_tsv(paths$selected)
  delta <- read_tsv(paths$delta)
  class_table <- read_tsv(paths$class_by_seed_time)
  convergence_summary <- read_tsv(paths$convergence)
  daily_manifest <- read_tsv(paths$daily_manifest)
  write_figures(selected, delta, class_table, convergence_summary, daily_manifest, output_root)
}

validate_cached_eigen <- function(fixo2_env, model_env, cfg, param_mat, seed_id, O2, time_grid) {
  fixed_matrix_fn <- get("cf2_fixed_matrix", envir = fixo2_env, inherits = TRUE)
  init_vector_fn <- get("cf2_init_vector", envir = fixo2_env, inherits = TRUE)
  cf2_eigen_fn <- get("cf2_eigen_trajectory", envir = fixo2_env, inherits = TRUE)
  run_params_fn <- get("o2pr_run_params_from_vec", envir = fixo2_env, inherits = TRUE)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
  names(pvec) <- colnames(param_mat)
  run_params <- run_params_fn(pvec, cfg)
  fm <- fixed_matrix_fn(model_env, cfg, run_params, O2)
  init <- init_vector_fn(fm$ngrid, 2 * n_unit)
  eig <- eigen(fm$M, only.values = FALSE)
  cached <- eigen_trajectory_cached(eig, fm$ngrid, init$vector, time_grid, n_unit)$trajectory
  original <- cf2_eigen_fn(fm$M, fm$ngrid, init$vector, time_grid, n_unit)$trajectory
  max(abs(cached$mean_ploidy - original$mean_ploidy), na.rm = TRUE)
}

build_validation <- function(daily_manifest, selected, delta, class_table, convergence_summary,
                             seeds, o2_grid, time_grid, selected_times, init_audit,
                             fixo2_env, model_env, cfg, param_mat, curve_utils_path) {
  expected_daily_rows <- length(o2_grid) * 2L * length(time_grid)
  expected_selected_rows <- length(seeds) * length(o2_grid) * 2L * length(selected_times)
  expected_delta_rows <- length(seeds) * length(o2_grid) * length(selected_times)
  expected_class_rows <- length(seeds) * length(selected_times) * 2L
  helper_rows <- run_curve_classification_validation()
  helper_validation <- data.frame(
    test_case = paste0("curve_helper_", helper_rows$test_case),
    expected = helper_rows$expected_class,
    observed = helper_rows$observed_class,
    passed = helper_rows$passed,
    stringsAsFactors = FALSE
  )
  compare_error <- tryCatch(
    validate_cached_eigen(fixo2_env, model_env, cfg, param_mat, seeds[[1L]], o2_grid[[1L]], selected_times),
    error = function(e) NA_real_
  )
  audit_2N <- init_audit[init_audit$initial_condition == "init_2N", , drop = FALSE]
  audit_4N <- init_audit[init_audit$initial_condition == "init_4N", , drop = FALSE]
  rows <- list(
    data.frame(test_case = "curve_classification_utils_sourced", expected = curve_utils_path, observed = curve_utils_path,
               passed = file.exists(curve_utils_path) && exists("classify_o2_ploidy_curve") && exists("finite_diff_curve"), stringsAsFactors = FALSE),
    data.frame(test_case = "classification_rule_version", expected = "global_range_v2", observed = curve_classification_rule_version(),
               passed = identical(curve_classification_rule_version(), "global_range_v2"), stringsAsFactors = FALSE),
    data.frame(test_case = "initial_ploidy_2N_maps_to_N44", expected = "44", observed = paste(sort(unique(audit_2N$used_initial_N)), collapse = ","),
               passed = nrow(audit_2N) > 0L && all(audit_2N$used_initial_N == 44), stringsAsFactors = FALSE),
    data.frame(test_case = "initial_ploidy_4N_maps_to_N88", expected = "88", observed = paste(sort(unique(audit_4N$used_initial_N)), collapse = ","),
               passed = nrow(audit_4N) > 0L && all(audit_4N$used_initial_N == 88), stringsAsFactors = FALSE),
    data.frame(test_case = "daily_gz_files_exist", expected = as.character(length(seeds)), observed = as.character(sum(file.exists(daily_manifest$daily_file))),
               passed = all(file.exists(daily_manifest$daily_file)), stringsAsFactors = FALSE),
    data.frame(test_case = "daily_gz_row_count_per_seed", expected = as.character(expected_daily_rows),
               observed = paste(sort(unique(daily_manifest$observed_rows)), collapse = ","), passed = all(daily_manifest$observed_rows == expected_daily_rows), stringsAsFactors = FALSE),
    data.frame(test_case = "daily_manifest_status_all_ok", expected = "ok",
               observed = paste(names(table(daily_manifest$status)), as.integer(table(daily_manifest$status)), collapse = ","),
               passed = all(daily_manifest$status == "ok"), stringsAsFactors = FALSE),
    data.frame(test_case = "selected_time_curve_rows", expected = as.character(expected_selected_rows), observed = as.character(nrow(selected)),
               passed = nrow(selected) == expected_selected_rows, stringsAsFactors = FALSE),
    data.frame(test_case = "selected_mean_ploidy_no_NA", expected = "0", observed = as.character(sum(is.na(selected$mean_ploidy))),
               passed = sum(is.na(selected$mean_ploidy)) == 0L, stringsAsFactors = FALSE),
    data.frame(test_case = "delta_by_seed_o2_time_rows", expected = as.character(expected_delta_rows), observed = as.character(nrow(delta)),
               passed = nrow(delta) == expected_delta_rows, stringsAsFactors = FALSE),
    data.frame(test_case = "delta_initial_no_NA", expected = "0", observed = as.character(sum(is.na(delta$delta_initial))),
               passed = sum(is.na(delta$delta_initial)) == 0L, stringsAsFactors = FALSE),
    data.frame(test_case = "curve_class_by_seed_time_rows", expected = as.character(expected_class_rows), observed = as.character(nrow(class_table)),
               passed = nrow(class_table) == expected_class_rows, stringsAsFactors = FALSE),
    data.frame(test_case = "selected_times_complete", expected = paste(selected_times, collapse = ","),
               observed = paste(sort(unique(selected$day)), collapse = ","), passed = setequal(sort(unique(selected$day)), selected_times), stringsAsFactors = FALSE),
    data.frame(test_case = "finite_difference_columns_present", expected = "finite_difference_next,local_slope_sign",
               observed = paste(intersect(c("finite_difference_next", "local_slope_sign"), names(selected)), collapse = ","),
               passed = all(c("finite_difference_next", "local_slope_sign") %in% names(selected)), stringsAsFactors = FALSE),
    data.frame(test_case = "cached_eigen_matches_cf2_eigen_trajectory", expected = "<=1e-10",
               observed = format(compare_error, scientific = TRUE), passed = is.finite(compare_error) && compare_error <= 1e-10, stringsAsFactors = FALSE),
    data.frame(test_case = "convergence_summary_nonempty", expected = ">0", observed = as.character(nrow(convergence_summary)),
               passed = nrow(convergence_summary) > 0L, stringsAsFactors = FALSE)
  )
  out <- rbind_nonempty(c(rows, list(helper_validation)))
  out$status <- ifelse(out$passed, "PASS", "FAIL")
  out[, c("test_case", "status", "expected", "observed", "passed")]
}

run_analysis <- function(args = parse_args()) {
  fit_root <- normalizePath(path.expand(args$fit_root %||% args$run_dir %||% default_fit_root()), mustWork = FALSE)
  output_root <- normalizePath(path.expand(args$output_root %||% args$out_dir %||% default_output_root()), mustWork = FALSE)
  o2_min <- as_num(args$o2_min, 0)
  o2_max <- as_num(args$o2_max, 5)
  o2_by <- as_num(args$o2_by, 0.025)
  o2_grid <- sort(unique(as_num_vec(args$o2_grid, seq(o2_min, o2_max, by = o2_by))))
  time_end <- as_int(args$time_end, 1000L)
  time_grid <- seq.int(0L, time_end)
  selected_times <- sort(unique(as_num_vec(args$selected_times, c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300, 500, 700, 1000))))
  selected_times <- selected_times[selected_times %in% time_grid]
  max_seeds <- as_int(args$max_seeds, NA_integer_)
  n_workers <- as_int(args$n_workers, 1L)
  force <- as_bool(args$force %||% args$overwrite, TRUE)
  flat_range_threshold <- as_num(args$flat_range_threshold, 0.05)
  step_epsilon_abs <- as_num(args$step_epsilon_abs, 1e-6)
  step_epsilon_fraction <- as_num(args$step_epsilon_fraction, 1e-4)
  reverse_fraction_tolerance <- as_num(args$reverse_fraction_tolerance, 0.05)
  plateau_min_points <- as_int(args$plateau_min_points, 3L)
  run_validation <- as_bool(args$run_validation, TRUE)
  generate_figures <- as_bool(args$generate_figures, TRUE)
  plot_only <- as_bool(args$plot_only, FALSE)

  table_dir <- file.path(output_root, "tables")
  daily_dir <- file.path(table_dir, "daily_trajectories")
  paths <- analysis_paths(output_root)

  if (plot_only) {
    message("Generating figures from existing tables under: ", output_root)
    write_figures_from_tables(output_root)
    message("Completed figure generation: ", paths$figures)
    return(invisible(paths))
  }

  if (!dir.exists(fit_root)) stop("fit_root does not exist: ", fit_root)
  if (!length(o2_grid)) stop("O2 grid is empty.")
  if (!length(selected_times)) stop("selected_times is empty after intersecting with time_grid.")

  dir.create(daily_dir, recursive = TRUE, showWarnings = FALSE)

  fixo2_env <- load_fixo2_env()
  collect_inputs_fn <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)
  params_wide_fn <- get("o2ipa_params_wide", envir = fixo2_env, inherits = TRUE)
  first_seed_cfg_fn <- get("o2pr_first_seed_cfg", envir = fixo2_env, inherits = TRUE)
  source_model_fn <- get("o2ipa_source_model", envir = fixo2_env, inherits = TRUE)
  fixed_matrix_fn <- get("cf2_fixed_matrix", envir = fixo2_env, inherits = TRUE)
  init_vector_fn <- get("cf2_init_vector", envir = fixo2_env, inherits = TRUE)
  run_params_fn <- get("o2pr_run_params_from_vec", envir = fixo2_env, inherits = TRUE)

  inputs <- collect_inputs_fn(fit_root, objective_source = "auto")
  param_mat <- params_wide_fn(inputs$params_long, "value")
  seeds <- rownames(param_mat)
  seeds <- seeds[order(seed_number(seeds))]
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seeds <- seeds[seq_len(min(max_seeds, length(seeds)))]
  }
  if (!length(seeds)) stop("No fitted seed parameters found under: ", fit_root)
  cfg <- first_seed_cfg_fn(inputs$manifest)
  model_env <- source_model_fn(dirname(FIXO2_SCRIPT_PATH))

  message(
    "Computing fixed-O2 initial-ploidy analytical trajectories: ",
    length(seeds), " seeds x ", length(o2_grid), " O2 x 2 initial ploidies x ",
    length(time_grid), " days"
  )

  worker <- function(seed_id) {
    process_one_seed(
      seed_id = seed_id,
      param_mat = param_mat,
      cfg = cfg,
      model_env = model_env,
      fixed_matrix_fn = fixed_matrix_fn,
      init_vector_fn = init_vector_fn,
      run_params_fn = run_params_fn,
      o2_grid = o2_grid,
      time_grid = time_grid,
      selected_times = selected_times,
      daily_dir = daily_dir,
      force = force
    )
  }
  n_workers <- max(1L, min(n_workers, length(seeds)))
  results <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    parallel::mclapply(seeds, worker, mc.cores = n_workers)
  } else {
    lapply(seeds, worker)
  }

  daily_manifest <- rbind_nonempty(lapply(results, `[[`, "manifest"))
  init_audit <- rbind_nonempty(lapply(results, `[[`, "init_audit"))
  selected_raw <- rbind_nonempty(lapply(results, `[[`, "selected"))
  if (!nrow(selected_raw)) stop("No selected trajectory rows were generated.")
  selected_raw <- selected_raw[order(selected_raw$seed_number, selected_raw$O2_pct, selected_raw$initial_ploidy, selected_raw$day), , drop = FALSE]

  class_info <- classify_selected_curves(
    selected_raw,
    flat_range_threshold = flat_range_threshold,
    step_epsilon_abs = step_epsilon_abs,
    step_epsilon_fraction = step_epsilon_fraction,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    plateau_min_points = plateau_min_points
  )
  selected <- attach_curve_differences(selected_raw, class_info$differences)
  delta <- build_delta_table(selected_raw)
  convergence_summary <- summarize_convergence(delta, time_end = time_end)

  run_args <- data.frame(
    argument = c(
      "script", "fit_root", "output_root", "fixo2_script", "curve_classification_utils",
      "n_seed", "max_seeds", "o2_grid", "n_o2", "time_grid", "n_time",
      "selected_times", "n_selected_time", "initial_ploidy_values",
      "expected_daily_rows_per_seed", "expected_selected_rows", "expected_delta_rows",
      "expected_class_rows", "force", "n_workers", "classification_rule_version",
      "flat_range_threshold", "step_epsilon_abs", "step_epsilon_fraction",
      "reverse_fraction_tolerance", "plateau_min_points", "generate_figures", "analytical_method",
      "daily_file_pattern"
    ),
    value = c(
      normalizePath(file.path(SCRIPT_DIR, "fixed_o2_initial_ploidy_trajectory.R"), mustWork = FALSE),
      fit_root,
      output_root,
      normalizePath(FIXO2_SCRIPT_PATH, mustWork = FALSE),
      CURVE_UTILS_PATH,
      as.character(length(seeds)),
      as.character(max_seeds),
      paste(format_num(o2_grid), collapse = ","),
      as.character(length(o2_grid)),
      paste(range(time_grid), collapse = ":"),
      as.character(length(time_grid)),
      paste(format_num(selected_times), collapse = ","),
      as.character(length(selected_times)),
      "2N,4N",
      as.character(length(o2_grid) * 2L * length(time_grid)),
      as.character(length(seeds) * length(o2_grid) * 2L * length(selected_times)),
      as.character(length(seeds) * length(o2_grid) * length(selected_times)),
      as.character(length(seeds) * length(selected_times) * 2L),
      as.character(force),
      as.character(n_workers),
      curve_classification_rule_version(),
      as.character(flat_range_threshold),
      as.character(step_epsilon_abs),
      as.character(step_epsilon_fraction),
      as.character(reverse_fraction_tolerance),
      as.character(plateau_min_points),
      as.character(generate_figures),
      "cf2_fixed_matrix + cf2_init_vector + cached eigen propagation equivalent to cf2_eigen_trajectory with Matrix::expm fallback for numerically singular eigen trajectories; no stochastic simulation; no refitting",
      file.path(daily_dir, "seed%03d.tsv.gz")
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(daily_manifest, paths$daily_manifest)
  write_tsv(selected, paths$selected)
  write_tsv(delta, paths$delta)
  write_tsv(class_info$class_table, paths$class_by_seed_time)
  write_tsv(convergence_summary, paths$convergence)
  write_tsv(run_args, paths$run_arguments)
  if (run_validation) {
    validation <- build_validation(
      daily_manifest = daily_manifest,
      selected = selected,
      delta = delta,
      class_table = class_info$class_table,
      convergence_summary = convergence_summary,
      seeds = seeds,
      o2_grid = o2_grid,
      time_grid = time_grid,
      selected_times = selected_times,
      init_audit = init_audit,
      fixo2_env = fixo2_env,
      model_env = model_env,
      cfg = cfg,
      param_mat = param_mat,
      curve_utils_path = CURVE_UTILS_PATH
    )
    write_tsv(validation, paths$validation)
    if (!all(validation$passed)) {
      stop("Validation failed: ", paste(validation$test_case[!validation$passed], collapse = ", "))
    }
  }
  if (generate_figures) {
    write_figures(selected, delta, class_info$class_table, convergence_summary, daily_manifest, output_root)
  }

  message("Completed fixed-O2 initial-ploidy trajectory analysis: ", output_root)
  invisible(paths)
}

if (identical(environment(), globalenv())) {
  run_analysis(parse_args())
}
