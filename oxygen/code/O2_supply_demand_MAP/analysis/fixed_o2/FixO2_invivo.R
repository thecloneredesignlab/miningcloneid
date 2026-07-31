#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(
      nzchar,
      vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      }, character(1))
    )
    if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
  }
})

source(file.path(SCRIPT_DIR, "..", "process_fingerprints", "process_fingerprint_utils.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "..", "process_fingerprints", "ploidy_regime_utils.R"), local = TRUE)
FIXO2_SIMULATION_SCRIPT <- normalizePath(
  file.path(SCRIPT_DIR, "..", "..", "simulation", "fix_o2_simulation.R"),
  mustWork = TRUE
)
FIXO2_SIMULATION_ENV <- new.env(parent = globalenv())
FIXO2_SIMULATION_ENV$commandArgs <- function(trailingOnly = FALSE) {
  if (isTRUE(trailingOnly)) character(0) else paste0("--file=", FIXO2_SIMULATION_SCRIPT)
}
sys.source(FIXO2_SIMULATION_SCRIPT, envir = FIXO2_SIMULATION_ENV, chdir = TRUE)
for (nm in c(
  "fixo2_fixed_matrix",
  "fixo2_init_vector",
  "fixo2_normalize_state",
  "fixo2_eigen_states",
  "fixo2_eigen_trajectory",
  "fixo2_expm_states",
  "fixo2_trajectory_with_fallback",
  "fixo2_dominant_from_eig",
  "fixo2_dominant_one",
  "fixo2_simulation_output_complete",
  "fixo2_filter_missing_simulation_tasks"
)) {
  assign(nm, get(nm, envir = FIXO2_SIMULATION_ENV, inherits = TRUE), envir = environment())
}

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required")
})

options(error = function() {
  traceback(2)
  quit(status = 1)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

# Attractor analysis helpers.
fo2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

fixo2_attractor_o2_grid <- function(args) {
  sort(unique(fixo2_as_num_vec(args$attractor_o2_grid, seq(0, 5, by = 0.05))))
}

fixo2_format_o2_list <- function(x, max_n = 18L) {
  x <- sort(unique(as.numeric(x)))
  labs <- format(x, scientific = FALSE, trim = TRUE)
  if (length(labs) > max_n) {
    labs <- c(labs[seq_len(max_n)], paste0("... (", length(x), " total)"))
  }
  paste(labs, collapse = ",")
}

fixo2_validate_mode_reference_o2 <- function(mode_reference_o2, attractor_o2_grid) {
  mode_reference_o2 <- suppressWarnings(as.numeric(mode_reference_o2))
  attractor_o2_grid <- sort(unique(suppressWarnings(as.numeric(attractor_o2_grid))))
  attractor_o2_grid <- attractor_o2_grid[is.finite(attractor_o2_grid)]
  if (!is.finite(mode_reference_o2)) {
    stop("--mode_reference_o2 must be a finite numeric O2 value.", call. = FALSE)
  }
  if (!length(attractor_o2_grid)) {
    stop("--attractor_o2_grid must contain at least one finite numeric O2 value.", call. = FALSE)
  }
  if (!any(abs(attractor_o2_grid - mode_reference_o2) < 1e-9)) {
    stop(
      "--mode_reference_o2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " is invalid. It must exactly match one value in --attractor_o2_grid. Available attractor O2 values: ",
      fixo2_format_o2_list(attractor_o2_grid),
      call. = FALSE
    )
  }
  invisible(mode_reference_o2)
}

fo2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

fo2_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

fo2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs", "report")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

fo2_metric_map <- function(path) {
  if (!file.exists(path)) return(list())
  tab <- tryCatch(fo2_read_tsv(path), error = function(e) data.frame())
  if (!all(c("metric", "value") %in% names(tab))) return(list())
  vals <- as.list(tab$value)
  names(vals) <- tab$metric
  vals
}

fo2_map_num <- function(map, keys) {
  for (key in keys) {
    val <- suppressWarnings(as.numeric(map[[key]]))
    if (length(val) && is.finite(val[[1]])) return(val[[1]])
  }
  NA_real_
}

fo2_map_chr <- function(map, keys) {
  for (key in keys) {
    val <- map[[key]]
    if (!is.null(val) && length(val) && !is.na(val[[1]]) && nzchar(as.character(val[[1]]))) return(as.character(val[[1]]))
  }
  NA_character_
}

fo2_seed_manifest_extras <- function(manifest) {
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    seed <- manifest$seed_id[[i]]
    seed_dir <- manifest$seed_dir[[i]]
    m <- if (!is.na(seed_dir) && nzchar(seed_dir)) fo2_metric_map(file.path(seed_dir, "fit_summary.tsv")) else list()
    data.frame(
      seed_id = seed,
      optimizer_iter_completed = fo2_map_num(m, c("optimizer_iter_completed", "deoptim_iter_completed", "itermax")),
      optimizer_stop_reason = fo2_map_chr(m, c("deoptim_stop_reason", "optimizer_stop_reason", "optimizer_local_convergence")),
      optimizer_local_objective = fo2_map_num(m, c("optimizer_local_objective")),
      optimizer_deoptim_objective = fo2_map_num(m, c("optimizer_deoptim_objective")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fixo2_mode_threshold <- function() 2

fixo2_mode_regimes <- function() {
  c(
    mode1 = "mode1_attractor_dominant_ploidy_ge_2",
    mode2 = "mode2_attractor_dominant_ploidy_lt_2"
  )
}

fixo2_mode_panel_names <- function() {
  c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
}

fixo2_mode_colors <- function(alpha = 1) {
  cols <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "#0072B2",
    mode2_attractor_dominant_ploidy_lt_2 = "#D55E00"
  )
  if (is.finite(alpha) && alpha < 1) cols <- grDevices::adjustcolor(cols, alpha.f = alpha)
  cols
}

fixo2_mode_labels_from_regime <- function(regime) {
  regimes <- fixo2_mode_regimes()
  out <- names(regimes)[match(as.character(regime), unname(regimes))]
  out[is.na(out)] <- NA_character_
  out
}

fixo2_o2_key <- function(x) {
  vapply(x, function(xx) format(signif(as.numeric(xx), 12), scientific = FALSE, trim = TRUE), character(1))
}

fixo2_mode_reference_o2 <- function(args = NULL) {
  val <- if (is.null(args)) NA_real_ else o2ipa_as_num(args$mode_reference_o2, NA_real_)
  if (!is.finite(val)) val <- 2
  val
}

fixo2_mode_fields <- function(dominant_ploidy) {
  dominant_ploidy <- suppressWarnings(as.numeric(dominant_ploidy))
  threshold <- fixo2_mode_threshold()
  regime <- ifelse(
    !is.finite(dominant_ploidy),
    NA_character_,
    ifelse(dominant_ploidy >= threshold, fixo2_mode_regimes()[["mode1"]], fixo2_mode_regimes()[["mode2"]])
  )
  data.frame(
    trajectory_regime = regime,
    mode_label = fixo2_mode_labels_from_regime(regime),
    mode_source = "fixed_o2_attractor_dominant_ploidy",
    mode_rule = "dominant_mean_ploidy >= 2 => mode1; dominant_mean_ploidy < 2 => mode2",
    mode_threshold_dominant_ploidy = threshold,
    stringsAsFactors = FALSE
  )
}

fixo2_assign_attractor_modes <- function(tab, ploidy_col = "dominant_mean_ploidy") {
  if (!nrow(tab)) return(tab)
  if (!ploidy_col %in% names(tab)) stop("Cannot assign FixO2 modes; missing column: ", ploidy_col)
  if ("trajectory_regime" %in% names(tab) && !"source_trajectory_regime" %in% names(tab)) tab$source_trajectory_regime <- tab$trajectory_regime
  if ("mode_label" %in% names(tab) && !"source_mode_label" %in% names(tab)) tab$source_mode_label <- tab$mode_label
  fields <- fixo2_mode_fields(tab[[ploidy_col]])
  replace_cols <- intersect(names(fields), names(tab))
  tab[, replace_cols] <- NULL
  cbind(tab, fields, stringsAsFactors = FALSE)
}

fixo2_attractor_mode_table <- function(attractors) {
  if (!nrow(attractors)) return(data.frame())
  d <- fixo2_assign_attractor_modes(attractors, "dominant_mean_ploidy")
  d$O2_key <- fixo2_o2_key(d$O2_pct)
  keep <- intersect(c(
    "seed_id", "O2_pct", "O2_key", "dominant_mean_ploidy", "trajectory_regime",
    "mode_label", "mode_source", "mode_rule", "mode_threshold_dominant_ploidy",
    "status", "dominant_growth_rate", "spectral_gap", "objective", "delta_objective",
    "in_attractor_o2_grid", "is_mode_reference_o2"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d[order(o2ipa_seed_number(d$seed_id), d$O2_pct), , drop = FALSE]
}

fixo2_attractor_mode_summary_by_seed <- function(mode_by_seed_o2, standard_o2 = c(0, 0.1, 0.5, 1, 2, 5)) {
  if (!nrow(mode_by_seed_o2)) return(data.frame())
  rows <- lapply(split(mode_by_seed_o2, mode_by_seed_o2$seed_id), function(d) {
    d <- d[order(d$O2_pct), , drop = FALSE]
    out <- data.frame(
      seed_id = d$seed_id[[1]],
      n_o2 = nrow(d),
      n_o2_mode1 = sum(d$mode_label == "mode1", na.rm = TRUE),
      n_o2_mode2 = sum(d$mode_label == "mode2", na.rm = TRUE),
      fraction_o2_mode1 = mean(d$mode_label == "mode1", na.rm = TRUE),
      fraction_o2_mode2 = mean(d$mode_label == "mode2", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (O2 in standard_o2) {
      key <- paste0("mode_at_o2_", gsub("[^0-9A-Za-z]+", "p", format(O2, scientific = FALSE, trim = TRUE)))
      hit <- d$mode_label[abs(as.numeric(d$O2_pct) - O2) < 1e-9]
      out[[key]] <- if (length(hit)) hit[[1]] else NA_character_
    }
    out
  })
  out <- do.call(rbind, rows)
  out[order(o2ipa_seed_number(out$seed_id)), , drop = FALSE]
}

fixo2_reference_mode_table <- function(mode_by_seed_o2, mode_reference_o2) {
  if (!nrow(mode_by_seed_o2)) return(data.frame())
  d <- mode_by_seed_o2[abs(as.numeric(mode_by_seed_o2$O2_pct) - mode_reference_o2) < 1e-9, , drop = FALSE]
  if (!nrow(d)) {
    stop(
      "No FixO2 attractor mode rows matched --mode_reference_o2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      ". Include this O2 value in the mode table or allow the workflow to compute it."
    )
  }
  d <- d[order(o2ipa_seed_number(d$seed_id)), , drop = FALSE]
  d <- d[!duplicated(d$seed_id), , drop = FALSE]
  threshold <- if ("mode_threshold_dominant_ploidy" %in% names(d)) d$mode_threshold_dominant_ploidy else rep(fixo2_mode_threshold(), nrow(d))
  out <- data.frame(
    seed_id = d$seed_id,
    mode_reference_o2_pct = as.numeric(d$O2_pct),
    mode_reference_o2_key = fixo2_o2_key(d$O2_pct),
    mode_reference_dominant_mean_ploidy = suppressWarnings(as.numeric(d$dominant_mean_ploidy)),
    trajectory_regime = d$trajectory_regime,
    mode_label = d$mode_label,
    mode_source = "fixed_o2_attractor_dominant_ploidy_at_reference_o2",
    mode_rule = paste0(
      "dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " >= 2 => mode1; dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " < 2 => mode2"
    ),
    mode_threshold_dominant_ploidy = threshold,
    stringsAsFactors = FALSE
  )
  optional_cols <- c(
    status = "mode_reference_status",
    dominant_growth_rate = "mode_reference_dominant_growth_rate",
    spectral_gap = "mode_reference_spectral_gap",
    objective = "objective",
    delta_objective = "delta_objective"
  )
  for (src in names(optional_cols)) {
    if (src %in% names(d)) out[[optional_cols[[src]]]] <- d[[src]]
  }
  out
}

fixo2_apply_reference_modes <- function(tab, reference_modes) {
  if (!nrow(tab) || !is.data.frame(reference_modes) || !nrow(reference_modes) || !"seed_id" %in% names(tab)) return(tab)
  mode_cols <- intersect(c(
    "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
    "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
    "mode_reference_dominant_mean_ploidy", "mode_reference_status",
    "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
  ), names(reference_modes))
  if (!all(c("seed_id", "trajectory_regime") %in% mode_cols)) return(tab)
  meta <- reference_modes[, mode_cols, drop = FALSE]
  meta <- meta[!duplicated(meta$seed_id), , drop = FALSE]
  replace_cols <- setdiff(mode_cols, "seed_id")
  tab[, intersect(replace_cols, names(tab))] <- NULL
  merge(tab, meta, by = "seed_id", all.x = TRUE, sort = FALSE)
}

fo2_load_optional_seed_metadata <- function(path, manifest) {
  if (nzchar(path) && file.exists(path)) {
    tab <- fo2_read_tsv(path)
    if (!"seed_id" %in% names(tab)) stop("label_file must contain seed_id: ", path)
    keep <- intersect(c(
      "seed_id", "trajectory_regime", "mode_label", "mode_reason", "mean_terminal_ploidy",
      "mean_late_drop_amplitude", "4N__terminal_median_ploidy",
      "4N__late_drop_amplitude", "4N__time_o2_near_floor"
    ), names(tab))
    out <- tab[, keep, drop = FALSE]
    names(out) <- sub("^trajectory_regime$", "source_trajectory_regime", names(out))
    names(out) <- sub("^mode_label$", "source_mode_label", names(out))
    names(out) <- sub("^mode_reason$", "source_mode_reason", names(out))
    return(out)
  }
  data.frame(seed_id = manifest$seed_id, stringsAsFactors = FALSE)
}

fo2_dominant_attractor_one <- function(seed_id, run_params, model_env, cfg, O2) {
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22L), as.integer(cfg$N_MAX %||% 154L))
  G <- o2pr_build_G(model_env, cfg, run_params, O2)
  mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
  M <- G - Matrix::Diagonal(x = mu_all)
  eig <- tryCatch(eigen(as.matrix(M), only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) {
    return(data.frame(
      seed_id = seed_id,
      O2_pct = O2,
      status = "eigen_failed",
      dominant_mean_N = NA_real_,
      dominant_mean_ploidy = NA_real_,
      dominant_fraction_N_le_25 = NA_real_,
      dominant_fraction_N_below_44 = NA_real_,
      dominant_fraction_N_ge_44 = NA_real_,
      dominant_growth_rate = NA_real_,
      spectral_gap = NA_real_,
      eigenvector_nonnegative = NA,
      selection_22_vs_44 = NA_real_,
      selection_44_vs_88 = NA_real_,
      eff_growth_N22 = NA_real_,
      eff_growth_N44 = NA_real_,
      eff_growth_N88 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  idx <- which.max(Re(eig$values))
  v <- Re(eig$vectors[, idx])
  if (sum(v, na.rm = TRUE) < 0) v <- -v
  nonneg <- all(v >= -1e-8, na.rm = TRUE)
  v <- pmax(v, 0)
  status <- "ok"
  if (!sum(v, na.rm = TRUE) > 0) {
    v <- rep(NA_real_, length(ngrid))
    status <- "empty_dominant_vector_after_truncation"
  } else {
    v <- v / sum(v, na.rm = TRUE)
  }
  lambda1 <- Re(eig$values[[idx]])
  lambda2 <- sort(Re(eig$values), decreasing = TRUE)[min(2L, length(eig$values))]
  eff <- vapply(c(22L, 44L, 88L), function(N) {
    col <- as.integer(N - (cfg$N_MIN %||% 22L) + 1L)
    if (col < 1L || col > ncol(G)) return(NA_real_)
    sum(G[, col]) - mu_all[[col]]
  }, numeric(1))
  names(eff) <- c("22", "44", "88")
  data.frame(
    seed_id = seed_id,
    O2_pct = O2,
    status = status,
    dominant_mean_N = sum(ngrid * v, na.rm = TRUE),
    dominant_mean_ploidy = sum(ngrid * v, na.rm = TRUE) / as.numeric(cfg$N_UNIT %||% 22),
    dominant_fraction_N_le_25 = sum(v[ngrid <= 25], na.rm = TRUE),
    dominant_fraction_N_below_44 = sum(v[ngrid < 44], na.rm = TRUE),
    dominant_fraction_N_ge_44 = sum(v[ngrid >= 44], na.rm = TRUE),
    dominant_growth_rate = lambda1,
    spectral_gap = lambda1 - lambda2,
    eigenvector_nonnegative = nonneg,
    selection_22_vs_44 = eff[["22"]] - eff[["44"]],
    selection_44_vs_88 = eff[["44"]] - eff[["88"]],
    eff_growth_N22 = eff[["22"]],
    eff_growth_N44 = eff[["44"]],
    eff_growth_N88 = eff[["88"]],
    stringsAsFactors = FALSE
  )
}

fo2_wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
                           a = "mode1_attractor_dominant_ploidy_ge_2", b = "mode2_attractor_dominant_ploidy_lt_2",
                           feature_name = value_col) {
  d <- tab[tab[[group_col]] %in% c(a, b) & is.finite(tab[[value_col]]), , drop = FALSE]
  if (!nrow(d) || length(unique(d[[group_col]])) < 2L) {
    return(data.frame(feature = feature_name, n_mode1 = 0L, n_mode2 = 0L, median_mode1 = NA_real_, median_mode2 = NA_real_, median_diff_mode2_minus_mode1 = NA_real_, wilcox_p_value = NA_real_))
  }
  x <- d[[value_col]][d[[group_col]] == a]
  y <- d[[value_col]][d[[group_col]] == b]
  p <- if (length(x) >= 1L && length(y) >= 1L) suppressWarnings(stats::wilcox.test(x, y)$p.value) else NA_real_
  data.frame(
    feature = feature_name,
    n_mode1 = length(x),
    n_mode2 = length(y),
    median_mode1 = stats::median(x, na.rm = TRUE),
    median_mode2 = stats::median(y, na.rm = TRUE),
    median_diff_mode2_minus_mode1 = stats::median(y, na.rm = TRUE) - stats::median(x, na.rm = TRUE),
    wilcox_p_value = p,
    stringsAsFactors = FALSE
  )
}

fo2_attractor_regime_tests <- function(attractors) {
  metrics <- c(
    "dominant_mean_ploidy", "dominant_fraction_N_le_25",
    "dominant_fraction_N_below_44", "dominant_growth_rate",
    "spectral_gap", "selection_22_vs_44", "selection_44_vs_88"
  )
  rows <- list()
  for (o2 in sort(unique(attractors$O2_pct))) {
    d <- attractors[attractors$O2_pct == o2, , drop = FALSE]
    for (m in metrics) {
      r <- fo2_wilcox_row(d, m, feature_name = paste0(m, "_O2_", gsub("[^0-9]+", "p", o2)))
      r$O2_pct <- o2
      r$metric <- m
      rows[[paste(o2, m, sep = "||")]] <- r
    }
  }
  out <- do.call(rbind, rows)
  out$BH_adjusted_p_value <- stats::p.adjust(out$wilcox_p_value, method = "BH")
  out[, c("O2_pct", "metric", "feature", "n_mode1", "n_mode2", "median_mode1", "median_mode2", "median_diff_mode2_minus_mode1", "wilcox_p_value", "BH_adjusted_p_value")]
}

fo2_parameter_correlations <- function(attractors, params_long) {
  if (!nrow(params_long)) return(data.frame())
  params_wide <- o2ipa_params_wide(params_long, "value")
  params_wide$seed_id <- rownames(params_wide)
  rows <- list()
  for (o2 in sort(unique(attractors$O2_pct))) {
    d0 <- attractors[attractors$O2_pct == o2, , drop = FALSE]
    d <- merge(d0, params_wide, by = "seed_id", all.x = TRUE, sort = FALSE)
    for (p in setdiff(names(params_wide), "seed_id")) {
      x <- suppressWarnings(as.numeric(d[[p]]))
      xtr <- if (all(x > 0, na.rm = TRUE)) log10(x) else x
      for (metric in c("dominant_mean_ploidy", "dominant_fraction_N_le_25", "selection_22_vs_44", "selection_44_vs_88")) {
        y <- suppressWarnings(as.numeric(d[[metric]]))
        ok <- is.finite(xtr) & is.finite(y)
        rows[[paste(o2, p, metric, sep = "||")]] <- data.frame(
          O2_pct = o2,
          parameter = p,
          metric = metric,
          n = sum(ok),
          spearman_rho = if (sum(ok) >= 3L) suppressWarnings(stats::cor(xtr[ok], y[ok], method = "spearman")) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$abs_spearman_rho <- abs(out$spearman_rho)
  out
}

fo2_mode_seed_stack_table <- function(attractors) {
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  ok_status <- if ("status" %in% names(attractors)) attractors$status == "ok" else rep(TRUE, nrow(attractors))
  d <- attractors[
    ok_status &
      attractors$trajectory_regime %in% regimes &
      is.finite(attractors$dominant_mean_ploidy),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(data.frame())
  d$seed_number <- o2ipa_seed_number(d$seed_id)
  d$stack_panel <- ifelse(
    d$trajectory_regime == "mode1_attractor_dominant_ploidy_ge_2",
    "mode1",
    "mode2"
  )
  keep <- intersect(c(
    "stack_panel", "seed_id", "seed_number", "trajectory_regime", "mode_label",
    "O2_pct", "dominant_mean_ploidy", "dominant_mean_N",
    "dominant_fraction_N_le_25", "dominant_fraction_N_below_44",
    "dominant_fraction_N_ge_44", "dominant_growth_rate", "spectral_gap",
    "selection_22_vs_44", "selection_44_vs_88",
    "objective", "delta_objective", "mean_terminal_ploidy", "mean_late_drop_amplitude"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d$trajectory_regime <- factor(d$trajectory_regime, levels = regimes)
  d <- d[order(d$trajectory_regime, d$seed_number, d$O2_pct), , drop = FALSE]
  d$trajectory_regime <- as.character(d$trajectory_regime)
  d
}

fo2_write_mode_comparison_tables <- function(attractors, out_dir) {
  regimes <- c(mode1 = "mode1_attractor_dominant_ploidy_ge_2", mode2 = "mode2_attractor_dominant_ploidy_lt_2")
  ok_status <- if ("status" %in% names(attractors)) attractors$status == "ok" else rep(TRUE, nrow(attractors))
  d <- attractors[
    ok_status &
      attractors$trajectory_regime %in% unname(regimes) &
      is.finite(attractors$dominant_mean_ploidy),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(invisible(FALSE))
  o2_vals <- sort(unique(d$O2_pct))
  rows <- lapply(o2_vals, function(O2) {
    row <- data.frame(O2_pct = O2, stringsAsFactors = FALSE)
    for (nm in names(regimes)) {
      vals <- d$dominant_mean_ploidy[d$O2_pct == O2 & d$trajectory_regime == regimes[[nm]]]
      row[[paste0(nm, "_n")]] <- length(vals)
      row[[paste0(nm, "_median")]] <- stats::median(vals, na.rm = TRUE)
      row[[paste0(nm, "_mean")]] <- mean(vals, na.rm = TRUE)
      row[[paste0(nm, "_q25")]] <- as.numeric(stats::quantile(vals, 0.25, na.rm = TRUE, names = FALSE))
      row[[paste0(nm, "_q75")]] <- as.numeric(stats::quantile(vals, 0.75, na.rm = TRUE, names = FALSE))
    }
    row
  })
  summary <- do.call(rbind, rows)
  fo2_write_tsv(summary, file.path(out_dir, "tables", "dominant_mean_ploidy_summary_mode1_mode2.tsv"))

  median_table <- data.frame(
    O2 = formatC(summary$O2_pct, format = "f", digits = 2),
    mode1 = round(summary$mode1_median, 3),
    mode2 = round(summary$mode2_median, 3),
    stringsAsFactors = FALSE
  )
  fo2_write_tsv(median_table, file.path(out_dir, "tables", "dominant_mean_ploidy_median_mode1_mode2.tsv"))

  display <- summary
  names(display)[names(display) == "O2_pct"] <- "O2"
  display$O2 <- formatC(display$O2, format = "f", digits = 2)
  numeric_cols <- setdiff(names(display), c("O2", "mode1_n", "mode2_n"))
  for (col in numeric_cols) display[[col]] <- round(display[[col]], 3)
  fo2_write_tsv(display, file.path(out_dir, "tables", "dominant_mean_ploidy_summary_mode1_mode2_display.tsv"))
  invisible(TRUE)
}

fo2_finite_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
}

fo2_finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

fo2_finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

fo2_spectral_gap_by_seed <- function(attractors) {
  ok_status <- if ("status" %in% names(attractors)) attractors$status == "ok" else rep(TRUE, nrow(attractors))
  d <- attractors[
    ok_status &
      is.finite(attractors$dominant_growth_rate) &
      is.finite(attractors$spectral_gap),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(data.frame())
  d$seed_number <- o2ipa_seed_number(d$seed_id)
  d$lambda1 <- d$dominant_growth_rate
  d$lambda2 <- d$dominant_growth_rate - d$spectral_gap
  d$relative_spectral_gap <- d$spectral_gap / pmax(abs(d$lambda1), .Machine$double.eps)
  d$relax_time_days <- ifelse(d$spectral_gap > 0, 1 / d$spectral_gap, NA_real_)
  d$time_to_10x_days <- ifelse(d$spectral_gap > 0, log(10) / d$spectral_gap, NA_real_)
  d$time_to_100x_days <- ifelse(d$spectral_gap > 0, log(100) / d$spectral_gap, NA_real_)
  d$log10_advantage_1000d <- d$spectral_gap * 1000 / log(10)
  d$dominance_class <- ifelse(
    !is.finite(d$spectral_gap) | d$spectral_gap <= 0, "non_positive",
    ifelse(d$spectral_gap < 0.001, "very_weak",
           ifelse(d$spectral_gap < 0.005, "weak",
                  ifelse(d$spectral_gap < 0.01, "moderate", "strong")))
  )
  d$dominance_class <- factor(d$dominance_class, levels = c("non_positive", "very_weak", "weak", "moderate", "strong"))
  keep <- intersect(c(
    "seed_id", "seed_number", "O2_pct", "trajectory_regime", "mode_label",
    "dominant_mean_ploidy", "lambda1", "lambda2", "spectral_gap",
    "relative_spectral_gap", "relax_time_days", "time_to_10x_days",
    "time_to_100x_days", "log10_advantage_1000d", "dominance_class",
    "dominant_fraction_N_le_25", "dominant_fraction_N_below_44",
    "objective", "delta_objective", "mean_terminal_ploidy", "mean_late_drop_amplitude"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d <- d[order(d$trajectory_regime, d$seed_number, d$O2_pct), , drop = FALSE]
  d$dominance_class <- as.character(d$dominance_class)
  d
}

fo2_spectral_gap_summary <- function(gap_by_seed) {
  if (!nrow(gap_by_seed)) return(data.frame())
  rows <- lapply(split(gap_by_seed, list(gap_by_seed$O2_pct, gap_by_seed$trajectory_regime), drop = TRUE), function(d) {
    data.frame(
      O2_pct = d$O2_pct[[1]],
      trajectory_regime = d$trajectory_regime[[1]],
      n_seed = length(unique(d$seed_id)),
      dominant_mean_ploidy_median = fo2_finite_median(d$dominant_mean_ploidy),
      spectral_gap_median = fo2_finite_median(d$spectral_gap),
      spectral_gap_q25 = fo2_finite_quantile(d$spectral_gap, 0.25),
      spectral_gap_q75 = fo2_finite_quantile(d$spectral_gap, 0.75),
      spectral_gap_min = suppressWarnings(min(d$spectral_gap, na.rm = TRUE)),
      relative_gap_median = fo2_finite_median(d$relative_spectral_gap),
      relax_time_days_median = fo2_finite_median(d$relax_time_days),
      time_to_10x_days_median = fo2_finite_median(d$time_to_10x_days),
      time_to_10x_days_q25 = fo2_finite_quantile(d$time_to_10x_days, 0.25),
      time_to_10x_days_q75 = fo2_finite_quantile(d$time_to_10x_days, 0.75),
      time_to_100x_days_median = fo2_finite_median(d$time_to_100x_days),
      log10_advantage_1000d_median = fo2_finite_median(d$log10_advantage_1000d),
      log10_advantage_1000d_q25 = fo2_finite_quantile(d$log10_advantage_1000d, 0.25),
      log10_advantage_1000d_q75 = fo2_finite_quantile(d$log10_advantage_1000d, 0.75),
      fraction_gap_lt_0p001 = mean(d$spectral_gap < 0.001, na.rm = TRUE),
      fraction_gap_lt_0p005 = mean(d$spectral_gap < 0.005, na.rm = TRUE),
      fraction_gap_lt_0p01 = mean(d$spectral_gap < 0.01, na.rm = TRUE),
      fraction_gap_ge_0p005 = mean(d$spectral_gap >= 0.005, na.rm = TRUE),
      fraction_gap_ge_0p01 = mean(d$spectral_gap >= 0.01, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$O2_pct, out$trajectory_regime), , drop = FALSE]
}

fo2_ploidy_gap_reliability_composite_table <- function(gap_by_seed) {
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  d <- gap_by_seed[
    gap_by_seed$trajectory_regime %in% regimes &
      is.finite(gap_by_seed$dominant_mean_ploidy) &
      is.finite(gap_by_seed$spectral_gap),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(data.frame())
  keep <- intersect(c(
    "seed_id", "seed_number", "trajectory_regime", "mode_label", "O2_pct",
    "dominant_mean_ploidy", "lambda1", "lambda2", "spectral_gap",
    "relative_spectral_gap", "relax_time_days", "time_to_10x_days",
    "time_to_100x_days", "log10_advantage_1000d", "dominance_class"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d$trajectory_regime <- factor(d$trajectory_regime, levels = regimes)
  d <- d[order(d$trajectory_regime, d$seed_number, d$O2_pct), , drop = FALSE]
  d$trajectory_regime <- as.character(d$trajectory_regime)
  d
}

fo2_shade_uncertain_o2 <- function(summary, reg, y_rng, o2_vals) {
  if (!nrow(summary)) return(invisible(FALSE))
  sub <- summary[
    summary$trajectory_regime == reg &
      (
        (!is.na(summary$spectral_gap_median) & summary$spectral_gap_median < 0.005) |
          (!is.na(summary$time_to_10x_days_median) & summary$time_to_10x_days_median > 500)
      ),
    ,
    drop = FALSE
  ]
  if (!nrow(sub)) return(invisible(FALSE))
  step <- suppressWarnings(min(diff(sort(unique(o2_vals))), na.rm = TRUE))
  if (!is.finite(step) || step <= 0) step <- 0.05
  half_step <- step / 2
  rect(
    xleft = pmax(min(o2_vals, na.rm = TRUE), sub$O2_pct - half_step),
    ybottom = y_rng[[1]],
    xright = pmin(max(o2_vals, na.rm = TRUE), sub$O2_pct + half_step),
    ytop = y_rng[[2]],
    col = grDevices::adjustcolor("grey80", alpha.f = 0.35),
    border = NA
  )
  invisible(TRUE)
}

fo2_plot_ploidy_gap_reliability_composite <- function(gap_by_seed, summary, fig_dir) {
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  d <- fo2_ploidy_gap_reliability_composite_table(gap_by_seed)
  d <- d[d$trajectory_regime %in% regimes, , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  d_gap <- d[d$spectral_gap > 0, , drop = FALSE]
  d_time <- d[is.finite(d$time_to_10x_days) & d$time_to_10x_days > 0, , drop = FALSE]
  if (!nrow(d_gap) || !nrow(d_time)) return(invisible(FALSE))

  panel_names <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
  o2_vals <- sort(unique(d$O2_pct))
  ploidy_rng <- range(d$dominant_mean_ploidy, 1, 1.5, 2, 2.5, na.rm = TRUE)
  ploidy_pad <- diff(ploidy_rng) * 0.04
  if (!is.finite(ploidy_pad) || ploidy_pad == 0) ploidy_pad <- 0.1
  ploidy_rng <- ploidy_rng + c(-ploidy_pad, ploidy_pad)
  gap_rng <- range(d_gap$spectral_gap, 0.001, 0.005, 0.01, na.rm = TRUE)
  time_rng <- range(d_time$time_to_10x_days, 100, 500, 1000, na.rm = TRUE)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.pdf"), width = 11, height = 10)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(4.1, 4.4, 2.1, 0.8), oma = c(0.6, 0.4, 2, 0.4))
  seed_lwd <- 0.95
  for (metric_row in c("ploidy", "gap", "time10x")) {
    for (j in seq_along(regimes)) {
      reg <- regimes[[j]]
      sub <- d[d$trajectory_regime == reg, , drop = FALSE]
      seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
      seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
      main_label <- if (metric_row == "ploidy") panel_names[[reg]] else ""
      x_axis <- if (metric_row == "time10x") "s" else "n"
      y_axis <- if (j == 1L) "s" else "n"
      if (metric_row == "ploidy") {
        plot(NA, xlim = range(o2_vals), ylim = ploidy_rng,
             xlab = "", ylab = if (j == 1L) "Dominant mean ploidy" else "",
             main = main_label, xaxt = x_axis, yaxt = y_axis)
        fo2_shade_uncertain_o2(summary, reg, ploidy_rng, o2_vals)
        for (seed in seed_ids) {
          ss <- sub[sub$seed_id == seed, , drop = FALSE]
          ss <- ss[order(ss$O2_pct), , drop = FALSE]
          lines(ss$O2_pct, ss$dominant_mean_ploidy, col = seed_col, lwd = seed_lwd)
        }
        abline(h = c(1, 1.5, 2, 2.5), col = "grey62", lty = 3, lwd = 1.4)
        if (j == 1L) {
          legend("topright", legend = c("uncertain O2 band", "single seed", "reference line"),
                 fill = c(grDevices::adjustcolor("grey80", alpha.f = 0.35), NA, NA),
                 border = c(NA, NA, NA), lty = c(NA, 1, 3),
                 col = c(NA, seed_col, "grey62"), lwd = c(NA, seed_lwd, 1.4),
                 bty = "n", cex = 0.75)
        }
      } else if (metric_row == "gap") {
        plot(NA, xlim = range(o2_vals), ylim = gap_rng, log = "y",
             xlab = "", ylab = if (j == 1L) "Spectral gap" else "",
             main = main_label, xaxt = x_axis, yaxt = y_axis)
        fo2_shade_uncertain_o2(summary, reg, gap_rng, o2_vals)
        for (seed in seed_ids) {
          ss <- sub[sub$seed_id == seed & sub$spectral_gap > 0, , drop = FALSE]
          ss <- ss[order(ss$O2_pct), , drop = FALSE]
          lines(ss$O2_pct, ss$spectral_gap, col = seed_col, lwd = seed_lwd)
        }
        abline(h = c(0.001, 0.005, 0.01), col = "grey55", lty = 3, lwd = 1.4)
      } else {
        plot(NA, xlim = range(o2_vals), ylim = time_rng, log = "y",
             xlab = "", ylab = if (j == 1L) "Days to 10x dominance" else "",
             main = main_label, xaxt = x_axis, yaxt = y_axis)
        fo2_shade_uncertain_o2(summary, reg, time_rng, o2_vals)
        for (seed in seed_ids) {
          ss <- sub[sub$seed_id == seed & is.finite(sub$time_to_10x_days) & sub$time_to_10x_days > 0, , drop = FALSE]
          ss <- ss[order(ss$O2_pct), , drop = FALSE]
          lines(ss$O2_pct, ss$time_to_10x_days, col = seed_col, lwd = seed_lwd)
        }
        abline(h = c(100, 500, 1000), col = "grey55", lty = 3, lwd = 1.4)
        mtext("Fixed O2 pct", side = 1, line = 2.7, cex = 0.9)
      }
    }
  }
  mtext("Fixed-O2 ploidy attractor and spectral-gap reliability", outer = TRUE, cex = 1.1, font = 2)
  par(oldpar)
  grDevices::dev.off()
  invisible(TRUE)
}

fo2_ploidy_gap_reliability_violin_table <- function(gap_by_seed,
                                                    o2_values = c(0, 0.1, 0.2, 0.5, 0.75, 1, 2, 3, 4, 5)) {
  regimes <- c(mode1 = "mode1_attractor_dominant_ploidy_ge_2", mode2 = "mode2_attractor_dominant_ploidy_lt_2")
  d <- fo2_ploidy_gap_reliability_composite_table(gap_by_seed)
  if (!nrow(d)) return(data.frame())
  target <- round(o2_values, 8)
  key <- round(d$O2_pct, 8)
  keep <- key %in% target & d$trajectory_regime %in% unname(regimes)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  d$O2_index <- match(round(d$O2_pct, 8), target)
  d$O2_label <- c("0", "0.1", "0.2", "0.5", "0.75", "1", "2", "3", "4", "5")[d$O2_index]
  d$mode_group <- ifelse(d$trajectory_regime == regimes[["mode1"]], "mode1", "mode2")
  d <- d[order(d$O2_index, d$mode_group, d$seed_number), , drop = FALSE]
  keep_cols <- intersect(c(
    "O2_index", "O2_label", "O2_pct", "mode_group", "seed_id", "seed_number",
    "trajectory_regime", "mode_label", "dominant_mean_ploidy", "spectral_gap",
    "time_to_10x_days", "time_to_100x_days", "log10_advantage_1000d",
    "dominance_class", "lambda1", "lambda2", "relative_spectral_gap",
    "relax_time_days"
  ), names(d))
  d[, keep_cols, drop = FALSE]
}

fo2_draw_violin_box <- function(y, x, width, fill_col, border_col = "grey25") {
  y <- y[is.finite(y)]
  if (!length(y)) return(invisible(FALSE))
  if (length(unique(y)) >= 2L) {
    dens <- stats::density(y, n = 128)
    scale <- if (max(dens$y, na.rm = TRUE) > 0) dens$y / max(dens$y, na.rm = TRUE) else dens$y
    graphics::polygon(
      c(x - scale * width, rev(x + scale * width)),
      c(dens$x, rev(dens$x)),
      col = fill_col,
      border = border_col,
      lwd = 0.7
    )
  } else {
    graphics::segments(x - width * 0.55, y[[1]], x + width * 0.55, y[[1]], col = border_col, lwd = 1.2)
  }
  graphics::boxplot(
    y, at = x, add = TRUE, axes = FALSE, outline = FALSE,
    boxwex = width * 0.55, col = grDevices::adjustcolor("white", alpha.f = 0.75),
    border = border_col, staplewex = 0.55
  )
  invisible(TRUE)
}

fo2_axis_ticks_log10 <- function(vals, rng) {
  vals <- vals[is.finite(vals) & vals > 0]
  vals <- vals[log10(vals) >= rng[[1]] & log10(vals) <= rng[[2]]]
  vals
}

fo2_plot_ploidy_gap_reliability_violin <- function(gap_by_seed, fig_dir,
                                                   o2_values = c(0, 0.1, 0.2, 0.5, 0.75, 1, 2, 3, 4, 5)) {
  d <- fo2_ploidy_gap_reliability_violin_table(gap_by_seed, o2_values = o2_values)
  if (!nrow(d)) return(invisible(FALSE))
  o2_labels <- c("0", "0.1", "0.2", "0.5", "0.75", "1", "2", "3", "4", "5")
  x_centers <- seq_along(o2_labels)
  offsets <- c(mode1 = -0.18, mode2 = 0.18)
  cols <- c(mode1 = "#0072B2", mode2 = "#D55E00")
  fill_cols <- stats::setNames(grDevices::adjustcolor(cols, alpha.f = 0.45), names(cols))
  violin_width <- 0.16

  ploidy_rng <- range(d$dominant_mean_ploidy, c(1, 1.5, 2, 2.5), na.rm = TRUE)
  ploidy_pad <- diff(ploidy_rng) * 0.05
  if (!is.finite(ploidy_pad) || ploidy_pad == 0) ploidy_pad <- 0.1
  ploidy_rng <- ploidy_rng + c(-ploidy_pad, ploidy_pad)

  gap_vals <- d$spectral_gap[d$spectral_gap > 0]
  time_vals <- d$time_to_10x_days[is.finite(d$time_to_10x_days) & d$time_to_10x_days > 0]
  gap_rng <- range(log10(gap_vals), log10(c(0.001, 0.005, 0.01)), na.rm = TRUE)
  time_rng <- range(log10(time_vals), log10(c(100, 500, 1000)), na.rm = TRUE)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.pdf"), width = 11, height = 9.5)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(3, 1), mar = c(2.5, 5, 2.1, 1), oma = c(3.3, 0.5, 1.5, 0.5))

  panel_specs <- list(
    list(metric = "dominant_mean_ploidy", ylab = "Dominant mean ploidy", ylim = ploidy_rng, log10 = FALSE,
         ref = c(1, 1.5, 2, 2.5), ref_label = "ploidy reference"),
    list(metric = "spectral_gap", ylab = "Spectral gap", ylim = gap_rng, log10 = TRUE,
         ref = c(0.001, 0.005, 0.01), ref_label = "gap threshold"),
    list(metric = "time_to_10x_days", ylab = "Days to 10x dominance", ylim = time_rng, log10 = TRUE,
         ref = c(100, 500, 1000), ref_label = "time reference")
  )

  for (i in seq_along(panel_specs)) {
    spec <- panel_specs[[i]]
    xaxt <- if (i == length(panel_specs)) "s" else "n"
    plot(NA, xlim = c(0.5, length(o2_labels) + 0.5), ylim = spec$ylim,
         axes = FALSE, xlab = "", ylab = spec$ylab, main = "")
    if (spec$log10) {
      ticks <- if (spec$metric == "spectral_gap") {
        fo2_axis_ticks_log10(c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1), spec$ylim)
      } else {
        fo2_axis_ticks_log10(c(10, 30, 100, 300, 500, 1000, 3000, 10000), spec$ylim)
      }
      axis(2, at = log10(ticks), labels = format(ticks, trim = TRUE, scientific = FALSE), las = 1)
    } else {
      axis(2, las = 1)
    }
    if (xaxt == "s") axis(1, at = x_centers, labels = o2_labels)
    box()
    for (idx in x_centers) {
      for (mode in names(offsets)) {
        vals <- d[d$O2_index == idx & d$mode_group == mode, spec$metric]
        vals <- suppressWarnings(as.numeric(vals))
        if (spec$log10) vals <- log10(vals[is.finite(vals) & vals > 0])
        fo2_draw_violin_box(vals, idx + offsets[[mode]], violin_width, fill_cols[[mode]], cols[[mode]])
      }
    }
    ref_y <- if (spec$log10) log10(spec$ref) else spec$ref
    graphics::abline(h = ref_y, col = "grey45", lty = 3, lwd = 1.5)
    if (i == 1L) {
      legend("topright",
             legend = c("mode1", "mode2", spec$ref_label),
             fill = c(fill_cols[["mode1"]], fill_cols[["mode2"]], NA),
             border = c(cols[["mode1"]], cols[["mode2"]], NA),
             lty = c(NA, NA, 3), col = c(NA, NA, "grey45"),
             lwd = c(NA, NA, 1.5), bty = "n", cex = 0.85)
    }
  }
  mtext("Fixed O2 pct", side = 1, outer = TRUE, line = 1.8, cex = 0.95)
  mtext("Fixed-O2 ploidy attractor and reliability distributions", side = 3, outer = TRUE, line = 0.3, cex = 1.1, font = 2)
  par(oldpar)
  grDevices::dev.off()
  invisible(TRUE)
}

fo2_plot_spectral_gap_outputs <- function(gap_by_seed, summary, fig_dir) {
  if (!nrow(gap_by_seed)) return(invisible(FALSE))
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  panel_names <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
  d <- gap_by_seed[gap_by_seed$trajectory_regime %in% regimes & gap_by_seed$spectral_gap > 0, , drop = FALSE]
  if (nrow(d)) {
    y_rng <- range(d$spectral_gap, 0.001, 0.005, 0.01, na.rm = TRUE)
    grDevices::pdf(file.path(fig_dir, "fixed_o2_spectral_gap_all_seed_stack_mode1_mode2.pdf"), width = 8.5, height = 8)
    oldpar <- par(no.readonly = TRUE)
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
    for (reg in regimes) {
      sub <- d[d$trajectory_regime == reg, , drop = FALSE]
      plot(NA, xlim = range(d$O2_pct), ylim = y_rng, log = "y",
           xlab = "Fixed O2 pct", ylab = "Spectral gap",
           main = panel_names[[reg]])
      abline(h = c(0.001, 0.005, 0.01), col = "grey80", lty = 3)
      seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
      seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
      for (seed in seed_ids) {
        ss <- sub[sub$seed_id == seed, , drop = FALSE]
        ss <- ss[order(ss$O2_pct), , drop = FALSE]
        lines(ss$O2_pct, ss$spectral_gap, col = seed_col, lwd = 0.7)
      }
      legend("topright", legend = c("0.001 / 0.005 / 0.01 thresholds", "single seed"),
             col = c("grey80", seed_col), lty = c(3, 1), lwd = c(1, 1), bty = "n")
    }
    par(oldpar)
    grDevices::dev.off()
  }

  d10 <- gap_by_seed[gap_by_seed$trajectory_regime %in% regimes & is.finite(gap_by_seed$time_to_10x_days) & gap_by_seed$time_to_10x_days > 0, , drop = FALSE]
  if (nrow(d10)) {
    y_rng <- range(d10$time_to_10x_days, 10, 100, 1000, na.rm = TRUE)
    grDevices::pdf(file.path(fig_dir, "fixed_o2_time_to_10x_all_seed_stack_mode1_mode2.pdf"), width = 8.5, height = 8)
    oldpar <- par(no.readonly = TRUE)
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
    for (reg in regimes) {
      sub <- d10[d10$trajectory_regime == reg, , drop = FALSE]
      plot(NA, xlim = range(d10$O2_pct), ylim = y_rng, log = "y",
           xlab = "Fixed O2 pct", ylab = "Days for dominant mode to reach 10x",
           main = panel_names[[reg]])
      abline(h = c(10, 100, 1000), col = "grey80", lty = 3)
      seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
      seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
      for (seed in seed_ids) {
        ss <- sub[sub$seed_id == seed, , drop = FALSE]
        ss <- ss[order(ss$O2_pct), , drop = FALSE]
        lines(ss$O2_pct, ss$time_to_10x_days, col = seed_col, lwd = 0.7)
      }
      legend("topright", legend = c("10 / 100 / 1000 days", "single seed"),
             col = c("grey80", seed_col), lty = c(3, 1), lwd = c(1, 1), bty = "n")
    }
    par(oldpar)
    grDevices::dev.off()
  }

  s <- summary[summary$trajectory_regime %in% regimes, , drop = FALSE]
  if (nrow(s)) {
    grDevices::pdf(file.path(fig_dir, "fixed_o2_gap_reliability_fraction_mode1_mode2.pdf"), width = 8.5, height = 6)
    oldpar <- par(no.readonly = TRUE)
    plot(NA, xlim = range(s$O2_pct), ylim = c(0, 1),
         xlab = "Fixed O2 pct", ylab = "Fraction of seeds",
         main = "Spectral gap reliability fractions")
    ltys <- c(fraction_gap_ge_0p005 = 2, fraction_gap_ge_0p01 = 1)
    for (reg in regimes) {
      sub <- s[s$trajectory_regime == reg, , drop = FALSE]
      sub <- sub[order(sub$O2_pct), , drop = FALSE]
      lines(sub$O2_pct, sub$fraction_gap_ge_0p005, col = cols[[reg]], lty = ltys[["fraction_gap_ge_0p005"]], lwd = 2)
      lines(sub$O2_pct, sub$fraction_gap_ge_0p01, col = cols[[reg]], lty = ltys[["fraction_gap_ge_0p01"]], lwd = 2)
    }
    legend("bottomright",
           legend = c("mode1 gap >= 0.005", "mode1 gap >= 0.01", "mode2 gap >= 0.005", "mode2 gap >= 0.01"),
           col = c(cols[["mode1_attractor_dominant_ploidy_ge_2"]], cols[["mode1_attractor_dominant_ploidy_ge_2"]], cols[["mode2_attractor_dominant_ploidy_lt_2"]], cols[["mode2_attractor_dominant_ploidy_lt_2"]]),
           lty = c(2, 1, 2, 1), lwd = 2, bty = "n")
    par(oldpar)
    grDevices::dev.off()
  }
  invisible(TRUE)
}

fo2_write_spectral_gap_outputs <- function(attractors, out_dir, generate_figures = TRUE) {
  gap_by_seed <- fo2_spectral_gap_by_seed(attractors)
  gap_summary <- fo2_spectral_gap_summary(gap_by_seed)
  fo2_write_tsv(gap_by_seed, file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_by_seed.tsv"))
  fo2_write_tsv(gap_summary, file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_regime_summary.tsv"))
  fo2_write_tsv(fo2_ploidy_gap_reliability_composite_table(gap_by_seed), file.path(out_dir, "tables", "fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv"))
  fo2_write_tsv(fo2_ploidy_gap_reliability_violin_table(gap_by_seed), file.path(out_dir, "tables", "fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv"))
  if (isTRUE(generate_figures)) {
    dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
    fo2_plot_spectral_gap_outputs(gap_by_seed, gap_summary, file.path(out_dir, "figures"))
    fo2_plot_ploidy_gap_reliability_composite(gap_by_seed, gap_summary, file.path(out_dir, "figures"))
    fo2_plot_ploidy_gap_reliability_violin(gap_by_seed, file.path(out_dir, "figures"))
  }
  invisible(list(by_seed = gap_by_seed, summary = gap_summary))
}

fo2_plot_mode_seed_stack <- function(attractors, fig_dir) {
  d <- fo2_mode_seed_stack_table(attractors)
  if (!nrow(d)) return(invisible(FALSE))
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  panel_names <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
  x_vals <- sort(unique(d$O2_pct))
  y_rng <- range(d$dominant_mean_ploidy, na.rm = TRUE)
  y_pad <- diff(y_rng) * 0.04
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.1
  y_rng <- y_rng + c(-y_pad, y_pad)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.pdf"), width = 8.5, height = 8)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
  for (reg in regimes) {
    sub <- d[d$trajectory_regime == reg, , drop = FALSE]
    plot(NA,
         xlim = range(x_vals, na.rm = TRUE), ylim = y_rng, xaxt = "n",
         xlab = "Fixed O2 pct", ylab = "Dominant mean ploidy",
         main = panel_names[[reg]])
    axis(1, at = x_vals, labels = format(x_vals, trim = TRUE, scientific = FALSE))
    abline(h = c(1, 1.5, 2, 2.5), col = "grey85", lty = 3)
    seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
    seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
    for (seed in seed_ids) {
      ss <- sub[sub$seed_id == seed, , drop = FALSE]
      ss <- ss[order(ss$O2_pct), , drop = FALSE]
      lines(ss$O2_pct, ss$dominant_mean_ploidy, col = seed_col, lwd = 0.7)
    }
    legend("topright", legend = "single seed", col = seed_col, lwd = 1, bty = "n")
  }
  invisible(TRUE)
}

fo2_plot_outputs <- function(attractors, out_dir) {
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00", unknown = "grey60")
  d <- attractors[attractors$status == "ok" & is.finite(attractors$dominant_mean_ploidy), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  grDevices::pdf(file.path(fig_dir, "fixed_o2_dominant_ploidy_by_regime.pdf"), width = 8, height = 6)
  plot(NA, xlim = range(d$O2_pct, na.rm = TRUE), ylim = range(d$dominant_mean_ploidy, na.rm = TRUE),
       xlab = "Fixed O2 pct", ylab = "Dominant mean ploidy", main = "Fixed-O2 attractor by regime")
  for (reg in unique(d$trajectory_regime)) {
    sub <- d[d$trajectory_regime == reg, , drop = FALSE]
    if (!nrow(sub)) next
    by_o2 <- aggregate(sub$dominant_mean_ploidy, by = list(O2_pct = sub$O2_pct), FUN = median, na.rm = TRUE)
    lines(by_o2$O2_pct, by_o2$x, col = cols[reg] %||% "grey60", lwd = 2, type = "b", pch = 16)
  }
  legend("topright", legend = names(cols), col = cols, lwd = 2, pch = 16, bty = "n")
  grDevices::dev.off()

  low <- d[d$O2_pct %in% sort(unique(d$O2_pct))[seq_len(min(4L, length(unique(d$O2_pct))))], , drop = FALSE]
  if (nrow(low)) {
    grDevices::pdf(file.path(fig_dir, "low_o2_dominant_ploidy_boxplot.pdf"), width = 9, height = 6)
    boxplot(dominant_mean_ploidy ~ interaction(O2_pct, trajectory_regime, drop = TRUE), data = low, las = 2,
            ylab = "Dominant mean ploidy", main = "Low fixed-O2 attractor distribution")
    grDevices::dev.off()
  }
  fo2_plot_mode_seed_stack(attractors, fig_dir)
  invisible(TRUE)
}

fo2_write_report <- function(out_dir, run_dir, label_file, attractors, regime_summary, tests, correlations) {
  counts <- if (nrow(attractors)) {
    rows <- lapply(split(attractors, attractors$trajectory_regime, drop = TRUE), function(d) {
      data.frame(
        trajectory_regime = d$trajectory_regime[[1]],
        n_seed_o2 = nrow(d),
        n_seed = length(unique(d$seed_id)),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  } else {
    data.frame()
  }
  low_o2 <- regime_summary[regime_summary$O2_pct %in% c(0, 0.01, 0.05, 0.1), , drop = FALSE]
  top_tests <- tests[is.finite(tests$BH_adjusted_p_value), , drop = FALSE]
  top_tests <- top_tests[order(top_tests$BH_adjusted_p_value, top_tests$wilcox_p_value), , drop = FALSE]
  top_corr <- correlations[is.finite(correlations$abs_spearman_rho), , drop = FALSE]
  top_corr <- top_corr[order(-top_corr$abs_spearman_rho), , drop = FALSE]
  gap_summary_path <- file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_regime_summary.tsv")
  gap_summary <- if (file.exists(gap_summary_path)) fo2_read_tsv(gap_summary_path) else data.frame()
  key_gap <- gap_summary[
    nrow(gap_summary) > 0 &
      gap_summary$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2") &
      gap_summary$O2_pct %in% c(0, 0.5, 1, 2, 5),
    ,
    drop = FALSE
  ]
  key_gap <- key_gap[, intersect(c(
    "O2_pct", "trajectory_regime", "n_seed", "spectral_gap_median",
    "time_to_10x_days_median", "time_to_100x_days_median",
    "log10_advantage_1000d_median", "fraction_gap_ge_0p005",
    "fraction_gap_ge_0p01"
  ), names(key_gap)), drop = FALSE]
  lines <- c(
    "# In vivo fixed-O2 ploidy attractor analysis",
    "",
    paste0("- run_dir: ", run_dir),
    paste0("- optional_label_file: ", if (nzchar(label_file)) label_file else "<none>"),
    "- FixO2 mode rule: dominant_mean_ploidy >= 2 is mode1; dominant_mean_ploidy < 2 is mode2.",
    paste0("- analyzed seeds: ", length(unique(attractors$seed_id))),
    "",
    "## Regime counts",
    "",
    paste(utils::capture.output(print(counts, row.names = FALSE)), collapse = "\n"),
    "",
    "## Low-O2 Attractor Summary",
    "",
    paste(utils::capture.output(print(low_o2, row.names = FALSE)), collapse = "\n"),
    "",
    "## Strongest Mode1-vs-Mode2 Attractor Differences",
    "",
    paste(utils::capture.output(print(utils::head(top_tests, 20), row.names = FALSE)), collapse = "\n"),
    "",
    "## Strongest Parameter-Attractor Correlations",
    "",
    paste(utils::capture.output(print(utils::head(top_corr, 30), row.names = FALSE)), collapse = "\n"),
    "",
    "## Seed-Level Mode1/Mode2 Stack Outputs",
    "",
    "- Table: `tables/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv`",
    "- Figure: `figures/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.pdf`",
    "",
    "## Spectral Gap Reliability Outputs",
    "",
    "- Seed table: `tables/fixed_o2_attractor_spectral_gap_by_seed.tsv`",
    "- Regime summary: `tables/fixed_o2_attractor_spectral_gap_regime_summary.tsv`",
    "- Figures: `figures/fixed_o2_spectral_gap_all_seed_stack_mode1_mode2.pdf`, `figures/fixed_o2_time_to_10x_all_seed_stack_mode1_mode2.pdf`, `figures/fixed_o2_gap_reliability_fraction_mode1_mode2.pdf`",
    "- Composite table: `tables/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv`",
    "- Composite figure: `figures/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.pdf`",
    "- Violin table: `tables/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv`",
    "- Violin figure: `figures/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.pdf`",
    "",
    paste(utils::capture.output(print(key_gap, row.names = FALSE)), collapse = "\n"),
    "",
    "## Interpretation Boundary",
    "",
    "These are fixed-O2 attractors from the fitted model parameters. They test model-internal mechanism consistency under standardized O2; they do not by themselves prove biological causality."
  )
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}


# Counterfactual trajectory helpers.
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

cf2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

cf2_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

cf2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs", "report")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

cf2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

cf2_default_time_grid <- function() {
  sort(unique(c(seq(0, 100, by = 1), 125, 150, 175, 200, 250, 300, 400, 500, 700, 1000)))
}

cf2_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  fixo2_fixed_matrix(model_env, cfg, run_params, O2)
}

cf2_init_vector <- function(ngrid, init_N) {
  fixo2_init_vector(ngrid, init_N, n_unit = 22)
}

cf2_normalize_state <- function(x) {
  fixo2_normalize_state(x)
}

cf2_eigen_trajectory <- function(M, ngrid, init, time_grid, n_unit) {
  fixo2_eigen_trajectory(M, ngrid, init, time_grid, n_unit)
}

cf2_crossing_time <- function(day, value, threshold, direction = "down") {
  ok <- is.finite(day) & is.finite(value)
  day <- day[ok]
  value <- value[ok]
  if (length(day) < 2L) return(NA_real_)
  ord <- order(day)
  day <- day[ord]
  value <- value[ord]
  hit <- if (identical(direction, "down")) which(value <= threshold) else which(value >= threshold)
  if (!length(hit)) return(NA_real_)
  i <- hit[[1]]
  if (i == 1L) return(day[[1]])
  x0 <- day[[i - 1L]]
  x1 <- day[[i]]
  y0 <- value[[i - 1L]]
  y1 <- value[[i]]
  if (!is.finite(y1 - y0) || abs(y1 - y0) < 1e-12) return(x1)
  x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
}

cf2_summarize_trajectory <- function(traj, max_time) {
  if (!nrow(traj)) {
    return(data.frame(
      terminal_mean_ploidy = NA_real_,
      terminal_fraction_N_le_25 = NA_real_,
      terminal_fraction_N_below_44 = NA_real_,
      min_mean_ploidy = NA_real_,
      day_min_mean_ploidy = NA_real_,
      time_crossing_ploidy_2p0_down = NA_real_,
      time_crossing_ploidy_1p8_down = NA_real_,
      time_crossing_ploidy_1p5_down = NA_real_,
      time_crossing_ploidy_1p3_down = NA_real_,
      time_crossing_ploidy_1p1_down = NA_real_,
      time_crossing_ploidy_1p5_down_censored = max_time + 1,
      reaches_ploidy_1p5 = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  last <- traj[which.max(traj$day), , drop = FALSE]
  i_min <- which.min(traj$mean_ploidy)
  t15 <- cf2_crossing_time(traj$day, traj$mean_ploidy, 1.5, "down")
  data.frame(
    terminal_mean_ploidy = last$mean_ploidy[[1]],
    terminal_fraction_N_le_25 = last$fraction_N_le_25[[1]],
    terminal_fraction_N_below_44 = last$fraction_N_below_44[[1]],
    min_mean_ploidy = traj$mean_ploidy[[i_min]],
    day_min_mean_ploidy = traj$day[[i_min]],
    time_crossing_ploidy_2p0_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 2.0, "down"),
    time_crossing_ploidy_1p8_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 1.8, "down"),
    time_crossing_ploidy_1p5_down = t15,
    time_crossing_ploidy_1p3_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 1.3, "down"),
    time_crossing_ploidy_1p1_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 1.1, "down"),
    time_crossing_ploidy_1p5_down_censored = if (is.finite(t15)) t15 else max_time + 1,
    reaches_ploidy_1p5 = is.finite(t15),
    stringsAsFactors = FALSE
  )
}

cf2_dominant_one <- function(M, ngrid, n_unit) {
  fixo2_dominant_one(M, ngrid, n_unit)
}

cf2_wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
                           a = "mode1_attractor_dominant_ploidy_ge_2", b = "mode2_attractor_dominant_ploidy_lt_2",
                           feature_name = value_col) {
  d <- tab[tab[[group_col]] %in% c(a, b) & is.finite(tab[[value_col]]), , drop = FALSE]
  if (!nrow(d) || length(unique(d[[group_col]])) < 2L) {
    return(data.frame(feature = feature_name, n_mode1 = 0L, n_mode2 = 0L, median_mode1 = NA_real_, mean_mode1 = NA_real_, median_mode2 = NA_real_, mean_mode2 = NA_real_, median_diff_mode2_minus_mode1 = NA_real_, mean_diff_mode2_minus_mode1 = NA_real_, wilcox_p_value = NA_real_))
  }
  x <- d[[value_col]][d[[group_col]] == a]
  y <- d[[value_col]][d[[group_col]] == b]
  p <- if (length(x) >= 1L && length(y) >= 1L) suppressWarnings(stats::wilcox.test(x, y)$p.value) else NA_real_
  data.frame(
    feature = feature_name,
    n_mode1 = length(x),
    n_mode2 = length(y),
    median_mode1 = stats::median(x, na.rm = TRUE),
    mean_mode1 = mean(x, na.rm = TRUE),
    median_mode2 = stats::median(y, na.rm = TRUE),
    mean_mode2 = mean(y, na.rm = TRUE),
    median_diff_mode2_minus_mode1 = stats::median(y, na.rm = TRUE) - stats::median(x, na.rm = TRUE),
    mean_diff_mode2_minus_mode1 = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
    wilcox_p_value = p,
    stringsAsFactors = FALSE
  )
}

cf2_regime_summary <- function(summary_by_seed) {
  rows <- list()
  counter <- 0L
  regs <- unname(fixo2_mode_regimes())
  for (O2 in sort(unique(summary_by_seed$O2_pct))) {
    for (init in unique(summary_by_seed$initial_condition)) {
      for (reg in regs) {
        d <- summary_by_seed[summary_by_seed$O2_pct == O2 & summary_by_seed$initial_condition == init & summary_by_seed$trajectory_regime == reg, , drop = FALSE]
        if (!nrow(d)) next
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          O2_pct = O2,
          initial_condition = init,
          trajectory_regime = reg,
          n_seed = nrow(d),
          terminal_mean_ploidy_median = stats::median(d$terminal_mean_ploidy, na.rm = TRUE),
          terminal_mean_ploidy_mean = mean(d$terminal_mean_ploidy, na.rm = TRUE),
          terminal_mean_ploidy_q25 = as.numeric(stats::quantile(d$terminal_mean_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
          terminal_mean_ploidy_q75 = as.numeric(stats::quantile(d$terminal_mean_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
          reaches_ploidy_1p5_fraction = mean(d$reaches_ploidy_1p5 %in% TRUE, na.rm = TRUE),
          time_crossing_1p5_censored_median = stats::median(d$time_crossing_ploidy_1p5_down_censored, na.rm = TRUE),
          time_crossing_1p5_censored_mean = mean(d$time_crossing_ploidy_1p5_down_censored, na.rm = TRUE),
          terminal_minus_dominant_median = stats::median(d$terminal_minus_dominant_ploidy, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

cf2_regime_tests <- function(summary_by_seed) {
  metrics <- c(
    "terminal_mean_ploidy", "terminal_fraction_N_le_25", "terminal_fraction_N_below_44",
    "min_mean_ploidy", "time_crossing_ploidy_1p5_down_censored",
    "time_crossing_ploidy_1p3_down", "terminal_minus_dominant_ploidy"
  )
  rows <- list()
  counter <- 0L
  for (O2 in sort(unique(summary_by_seed$O2_pct))) {
    for (init in unique(summary_by_seed$initial_condition)) {
      d <- summary_by_seed[summary_by_seed$O2_pct == O2 & summary_by_seed$initial_condition == init, , drop = FALSE]
      for (metric in metrics) {
        counter <- counter + 1L
        r <- cf2_wilcox_row(d, metric, feature_name = paste(metric, "O2", O2, init, sep = "__"))
        r$O2_pct <- O2
        r$initial_condition <- init
        r$metric <- metric
        rows[[counter]] <- r
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) out$BH_adjusted_p_value <- stats::p.adjust(out$wilcox_p_value, method = "BH")
  out[order(out$BH_adjusted_p_value, out$O2_pct, out$initial_condition, out$metric), , drop = FALSE]
}

cf2_parameter_correlations <- function(summary_by_seed, params_long) {
  if (!nrow(summary_by_seed) || !nrow(params_long)) return(data.frame())
  params <- params_long
  params$transformed_value <- mapply(o2ipa_transform_parameter_value, params$parameter, params$value)
  params_wide <- o2ipa_params_wide(params, value_col = "transformed_value")
  params_wide$seed_id <- rownames(params_wide)
  feature_cols <- c(
    "terminal_mean_ploidy", "terminal_fraction_N_le_25", "terminal_fraction_N_below_44",
    "time_crossing_ploidy_1p5_down_censored", "terminal_minus_dominant_ploidy"
  )
  merged <- merge(summary_by_seed[, c("seed_id", "trajectory_regime", "O2_pct", "initial_condition", feature_cols), drop = FALSE], params_wide, by = "seed_id", all.x = TRUE)
  scopes <- list(
    all_seeds = merged,
    mode1_mode2_only = merged[merged$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  )
  rows <- list()
  counter <- 0L
  for (scope_name in names(scopes)) {
    ds <- scopes[[scope_name]]
    for (O2 in sort(unique(ds$O2_pct))) {
      for (init in unique(ds$initial_condition)) {
        dd <- ds[ds$O2_pct == O2 & ds$initial_condition == init, , drop = FALSE]
        for (feature in feature_cols) {
          y <- dd[[feature]]
          for (param in o2ipa_target_params()) {
            x <- dd[[param]]
            ok <- is.finite(x) & is.finite(y)
            if (sum(ok) < 5L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) next
            ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
            counter <- counter + 1L
            rows[[counter]] <- data.frame(
              correlation_scope = scope_name,
              O2_pct = O2,
              initial_condition = init,
              feature = feature,
              parameter = param,
              parameter_transform = o2ipa_optimizer_transform(param),
              n = sum(ok),
              spearman_rho = unname(ct$estimate),
              abs_spearman_rho = abs(unname(ct$estimate)),
              p_value = ct$p.value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$BH_adjusted_p_value <- ave(out$p_value, out$correlation_scope, FUN = function(p) stats::p.adjust(p, method = "BH"))
  out[order(out$correlation_scope, -out$abs_spearman_rho, out$BH_adjusted_p_value), , drop = FALSE]
}

cf2_plot <- function(trajectory, summary_by_seed, fig_dir) {
  d <- trajectory[trajectory$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_counterfactual_median_trajectories.pdf"), width = 14, height = 12, onefile = FALSE, bg = "white")
  median_pdf_open <- TRUE
  on.exit(if (isTRUE(median_pdf_open)) grDevices::dev.off(), add = TRUE)

  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  n_o2 <- length(o2_use)
  o2_cols <- if (n_o2 >= 3L) 3L else max(1L, n_o2)
  o2_rows <- ceiling(n_o2 / o2_cols)
  if (n_o2 == 6L) {
    o2_rows <- 2L
    o2_cols <- 3L
  }
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
  line_cols <- setNames(grDevices::adjustcolor(unname(cols), alpha.f = 0.20), names(cols))
  reg_labels <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1 dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2 dominant ploidy < 2"
  )
  x_rng <- range(d$day, na.rm = TRUE)
  y_rng <- range(d$mean_ploidy, na.rm = TRUE)

  par(mfrow = c(o2_rows * length(init_use), o2_cols), mar = c(2.0, 4.2, 2.7, 1.0), oma = c(3.0, 3.2, 3.6, 0.5))
  for (o2_row in seq_len(o2_rows)) {
    row_o2 <- o2_use[((o2_row - 1L) * o2_cols + 1L):min(o2_row * o2_cols, n_o2)]
    for (init_idx in seq_along(init_use)) {
      init <- init_use[[init_idx]]
      for (col_idx in seq_len(o2_cols)) {
        if (col_idx > length(row_o2)) {
          plot.new()
          next
        }
        O2 <- row_o2[[col_idx]]
        d_panel <- d[d$O2_pct == O2 & d$initial_condition == init, , drop = FALSE]
        show_x <- init_idx == length(init_use)
        show_y <- col_idx == 1L
        plot(
          NA,
          xlim = x_rng,
          ylim = y_rng,
          xlab = "",
          ylab = "",
          xaxt = if (show_x) "s" else "n",
          yaxt = if (show_y) "s" else "n",
          main = if (init_idx == 1L) paste0("O2 = ", O2, "%; initial condition: ", init) else paste("Initial condition:", init),
          cex.main = 0.9
        )
        abline(h = c(1.5, 2), col = "grey80", lty = 2)
        for (reg in names(cols)) {
          sub <- d_panel[d_panel$trajectory_regime == reg, , drop = FALSE]
          if (!nrow(sub)) next
          seed_traces <- split(sub, sub$seed_id)
          for (trace in seed_traces) {
            trace <- trace[order(trace$day), , drop = FALSE]
            lines(trace$day, trace$mean_ploidy, col = line_cols[[reg]], lwd = 0.25)
          }
        }
        if (o2_row == 1L && init_idx == 1L && col_idx == 1L) {
          legend(
            "topright",
            legend = unname(reg_labels[names(cols)]),
            col = unname(cols),
            lty = 1,
            lwd = 2,
            bty = "n",
            cex = 0.72
          )
        }
      }
    }
  }
  mtext("Day", side = 1, outer = TRUE, line = 1.5)
  mtext("Mean ploidy", side = 2, outer = TRUE, line = 1.8)
  mtext("Fixed-O2 seed trajectories", side = 3, outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  median_pdf_open <- FALSE

  grDevices::pdf(file.path(fig_dir, "fixed_o2_counterfactual_terminal_boxplots.pdf"), width = 10, height = 7, bg = "white")
  d2 <- summary_by_seed[summary_by_seed$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  boxplot(terminal_mean_ploidy ~ interaction(O2_pct, initial_condition, trajectory_regime, drop = TRUE),
          data = d2, las = 2, ylab = "Terminal mean ploidy", main = "Fixed-O2 counterfactual terminal ploidy")
  grDevices::dev.off()
  invisible(TRUE)
}

cf2_report <- function(out_dir, args, regime_summary, tests, correlations) {
  key_summary <- regime_summary[regime_summary$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2") &
                                  regime_summary$O2_pct %in% c(1, 2, 5), , drop = FALSE]
  key_tests <- tests[tests$metric %in% c("terminal_mean_ploidy", "time_crossing_ploidy_1p5_down_censored") &
                       tests$O2_pct %in% c(1, 2, 5), , drop = FALSE]
  key_tests <- key_tests[order(key_tests$BH_adjusted_p_value), , drop = FALSE]
  top_corr <- data.frame()
  if (is.data.frame(correlations) && nrow(correlations) &&
      all(c("correlation_scope", "O2_pct", "abs_spearman_rho", "BH_adjusted_p_value") %in% names(correlations))) {
    top_corr <- correlations[correlations$correlation_scope == "mode1_mode2_only" & correlations$O2_pct %in% c(1, 2, 5), , drop = FALSE]
    if (nrow(top_corr)) {
      top_corr <- head(top_corr[order(-top_corr$abs_spearman_rho, top_corr$BH_adjusted_p_value), , drop = FALSE], 30L)
    }
  }
  lines <- c(
    "# In vivo fixed-O2 trajectory counterfactual analysis",
    "",
    paste0("- run_dir: `", o2ipa_as_chr(args$run_dir, ""), "`"),
    paste0("- optional_label_file: `", o2ipa_as_chr(args$label_file, ""), "`"),
    "- FixO2 mode rule: dominant_mean_ploidy >= 2 is mode1; dominant_mean_ploidy < 2 is mode2.",
    paste0("- O2 grid: `", o2ipa_as_chr(args$o2_grid, "0,0.5,1,2,5"), "`"),
    paste0("- time grid: default dense early grid unless `--time_grid` was supplied"),
    "",
    "## Key high-O2 terminal summaries",
    "",
    paste(capture.output(print(key_summary, row.names = FALSE)), collapse = "\n"),
    "",
    "## Key mode1-vs-mode2 tests",
    "",
    paste(capture.output(print(head(key_tests, 30L), row.names = FALSE)), collapse = "\n"),
    "",
    "## Top mode1/mode2 parameter correlations at O2=1/2/5",
    "",
    paste(capture.output(print(top_corr, row.names = FALSE)), collapse = "\n"),
    "",
    "## Notes",
    "",
    "- This counterfactual uses the same fixed-O2 matrix as the attractor analysis.",
    "- Initial conditions are point masses at N=44 (`init_2N`) and N=88 (`init_4N`) unless clipped by the configured N grid.",
    "- State vectors are propagated through the exact eigen-decomposition of the fixed linear system and normalized to composition at each time point.",
    "- `time_crossing_ploidy_1p5_down_censored` is set to max_time + 1 when the trajectory never crosses 1.5."
  )
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}


# Simulation validation helpers.
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

vf2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

vf2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

vf2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

vf2_as_seed_vec <- function(x, default) {
  vals <- o2ipa_split_csv(x, default)
  vals <- unique(o2ipa_norm_seed(vals))
  if (length(vals)) vals else default
}

fixo2_arg_nonempty <- function(x) {
  !is.null(x) && length(x) && !is.na(x[[1]]) && nzchar(trimws(as.character(x[[1]])))
}

vf2_seed_modes <- function(seed_ids, mode_arg = NULL) {
  modes <- o2ipa_split_csv(mode_arg, character())
  if (length(modes) == length(seed_ids)) {
    mode_map <- setNames(modes, seed_ids)
  } else {
    mode_map <- character()
  }
  out <- rep("", length(seed_ids))
  names(out) <- seed_ids
  has_mode <- seed_ids %in% names(mode_map) & nzchar(mode_map[seed_ids])
  out[has_mode] <- unname(mode_map[seed_ids[has_mode]])
  out
}

vf2_seed_labels <- function(seed_ids, mode_arg = NULL) {
  seed_modes <- vf2_seed_modes(seed_ids, mode_arg)
  labels <- seed_ids
  has_mode <- nzchar(seed_modes)
  labels[has_mode] <- paste0(seed_ids[has_mode], " (", seed_modes[has_mode], ")")
  labels
}

fixo2_simulation_mode_objective_metadata <- function(analysis_dir, run_dir, o2_grid,
                                                     mode_reference_o2, seed_ids = NULL,
                                                     n_workers = 1L) {
  dat <- read_seed_objectives(
    analysis_dir = analysis_dir,
    fit_dir = run_dir,
    o2_values = o2_grid,
    seed_ids = seed_ids,
    n_workers = n_workers,
    mode_reference_o2 = mode_reference_o2
  )
  if (!nrow(dat)) return(data.frame())
  dat$seed_id <- o2ipa_norm_seed(dat$seed_id)
  if (!"objective" %in% names(dat)) dat$objective <- NA_real_
  dat$objective <- suppressWarnings(as.numeric(dat$objective))
  if (!"objective_source" %in% names(dat)) dat$objective_source <- NA_character_
  if (!"mode_label" %in% names(dat)) dat$mode_label <- NA_character_
  if (!"trajectory_regime" %in% names(dat)) dat$trajectory_regime <- NA_character_
  dat <- dat[!duplicated(dat$seed_id), , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids)) dat <- dat[dat$seed_id %in% o2ipa_norm_seed(seed_ids), , drop = FALSE]
  dat
}

fixo2_select_best_objective_seed_by_mode <- function(metadata, mode_reference_o2) {
  required_modes <- c("mode1", "mode2")
  if (!nrow(metadata)) stop("Cannot select simulation representative seeds; no seed mode/objective metadata was available.")
  rows <- lapply(required_modes, function(mode) {
    d <- metadata[metadata$mode_label == mode & is.finite(metadata$objective), , drop = FALSE]
    if (!nrow(d)) {
      stop(
        "Cannot select a simulation representative seed for ",
        mode,
        "; no finite final objective was available among seeds assigned to this mode at O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        "."
      )
    }
    d <- d[order(d$objective, o2ipa_seed_number(d$seed_id)), , drop = FALSE]
    d <- d[1L, , drop = FALSE]
    data.frame(
      selection_mode = "best_final_objective_by_reference_mode",
      selection_rule = paste0(
        "Within ",
        mode,
        ", choose the seed with the smallest final objective among seeds assigned by the fixed-O2 attractor at O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE)
      ),
      mode_label = mode,
      trajectory_regime = d$trajectory_regime[[1]],
      seed_id = d$seed_id[[1]],
      final_objective = d$objective[[1]],
      objective_source = d$objective_source[[1]],
      mode_reference_o2_pct = if ("mode_reference_o2_pct" %in% names(d)) d$mode_reference_o2_pct[[1]] else mode_reference_o2,
      mode_reference_dominant_mean_ploidy = if ("mode_reference_dominant_mean_ploidy" %in% names(d)) d$mode_reference_dominant_mean_ploidy[[1]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$selection_rank <- seq_len(nrow(out))
  out[, c("selection_rank", setdiff(names(out), "selection_rank")), drop = FALSE]
}

fixo2_manual_simulation_seed_table <- function(seed_ids, seed_modes, metadata, mode_reference_o2) {
  seed_ids <- o2ipa_norm_seed(seed_ids)
  rows <- lapply(seq_along(seed_ids), function(i) {
    seed <- seed_ids[[i]]
    meta <- metadata[metadata$seed_id == seed, , drop = FALSE]
    mode <- seed_modes[[i]]
    if (!nzchar(mode) && nrow(meta) && "mode_label" %in% names(meta)) mode <- as.character(meta$mode_label[[1]])
    data.frame(
      selection_rank = i,
      selection_mode = "manual_seed_ids",
      selection_rule = "Use the explicit seed_ids supplied to the simulation validation step; mode and objective are attached from the FixO2 reference mode/objective table when available.",
      mode_label = mode,
      trajectory_regime = if (nrow(meta) && "trajectory_regime" %in% names(meta)) meta$trajectory_regime[[1]] else NA_character_,
      seed_id = seed,
      final_objective = if (nrow(meta) && "objective" %in% names(meta)) meta$objective[[1]] else NA_real_,
      objective_source = if (nrow(meta) && "objective_source" %in% names(meta)) meta$objective_source[[1]] else NA_character_,
      mode_reference_o2_pct = if (nrow(meta) && "mode_reference_o2_pct" %in% names(meta)) meta$mode_reference_o2_pct[[1]] else mode_reference_o2,
      mode_reference_dominant_mean_ploidy = if (nrow(meta) && "mode_reference_dominant_mean_ploidy" %in% names(meta)) meta$mode_reference_dominant_mean_ploidy[[1]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

vf2_as_int_vec <- function(x, default) {
  vals <- suppressWarnings(as.integer(o2ipa_split_csv(x, as.character(default))))
  vals <- vals[!is.na(vals) & vals > 0L]
  if (length(vals)) unique(vals) else default
}

vf2_num_path_tag <- function(x) {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) != 1L || !is.finite(val)) return("NA")
  txt <- format(val, scientific = FALSE, trim = TRUE)
  txt <- sub("^-", "m", txt)
  txt <- gsub("\\.", "p", txt)
  txt <- gsub("[^A-Za-z0-9]+", "", txt)
  if (!nzchar(txt)) "NA" else txt
}

vf2_default_time_grid <- function() {
  sort(unique(c(seq(0, 100, by = 1), 125, 150, 175, 200, 250, 300, 400, 500, 700, 1000)))
}

vf2_normalize_state <- function(x) {
  x <- as.numeric(Re(x))
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  x <- pmax(x, 0)
  s <- sum(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  x / s
}

vf2_init_vector <- function(ngrid, init_N) {
  idx <- which.min(abs(ngrid - init_N))
  v <- numeric(length(ngrid))
  v[[idx]] <- 1
  list(vector = v, used_N = ngrid[[idx]], used_ploidy = ngrid[[idx]] / 22)
}

vf2_state_metrics <- function(state, ngrid, n_unit) {
  w <- as.numeric(Re(state))
  w[!is.finite(w)] <- 0
  w <- pmax(w, 0)
  total <- sum(w, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(data.frame(
      mean_N = NA_real_,
      mean_ploidy = NA_real_,
      sd_ploidy = NA_real_,
      fraction_N_le_25 = NA_real_,
      fraction_N_below_44 = NA_real_,
      fraction_N_ge_44 = NA_real_,
      fraction_N_ge_66 = NA_real_,
      fraction_N_ge_88 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  w <- w / total
  ploidy_grid <- ngrid / n_unit
  mean_N <- sum(ngrid * w, na.rm = TRUE)
  mean_ploidy <- sum(ploidy_grid * w, na.rm = TRUE)
  second_ploidy <- sum(ploidy_grid^2 * w, na.rm = TRUE)
  sd_ploidy <- sqrt(max(0, second_ploidy - mean_ploidy^2))
  data.frame(
    mean_N = mean_N,
    mean_ploidy = mean_ploidy,
    sd_ploidy = sd_ploidy,
    fraction_N_le_25 = sum(w[ngrid <= 25], na.rm = TRUE),
    fraction_N_below_44 = sum(w[ngrid < 44], na.rm = TRUE),
    fraction_N_ge_44 = sum(w[ngrid >= 44], na.rm = TRUE),
    fraction_N_ge_66 = sum(w[ngrid >= 66], na.rm = TRUE),
    fraction_N_ge_88 = sum(w[ngrid >= 88], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

vf2_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22L), as.integer(cfg$N_MAX %||% 154L))
  G <- o2pr_build_G(model_env, cfg, run_params, O2)
  mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
  M <- G - Matrix::Diagonal(x = mu_all)
  list(M = M, G = G, mu_all = mu_all, ngrid = ngrid)
}

vf2_eigen_states <- function(M, init, time_grid) {
  Mdense <- as.matrix(M)
  eig <- tryCatch(eigen(Mdense, only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) stop("eigen decomposition failed")
  coef <- tryCatch(solve(eig$vectors, init), error = function(e) NULL)
  if (is.null(coef)) stop("eigen coefficient solve failed")
  lambda_ref <- max(Re(eig$values), na.rm = TRUE)
  states <- vector("list", length(time_grid))
  for (i in seq_along(time_grid)) {
    tt <- time_grid[[i]]
    weights <- exp((eig$values - lambda_ref) * tt) * coef
    states[[i]] <- vf2_normalize_state(eig$vectors %*% weights)
  }
  list(states = states, lambda_ref = lambda_ref)
}

vf2_expm_states <- function(M, init, time_grid) {
  time_grid <- sort(unique(as.numeric(time_grid)))
  if (!length(time_grid) || time_grid[[1]] < -1e-12) stop("time_grid must start at a non-negative time")
  x <- as.numeric(init)
  t_now <- 0
  states <- vector("list", length(time_grid))
  expm_cache <- new.env(parent = emptyenv())
  get_step_expm <- function(delta) {
    key <- format(signif(delta, 15), scientific = FALSE, trim = TRUE)
    if (!exists(key, envir = expm_cache, inherits = FALSE)) {
      assign(key, Matrix::expm(M * delta), envir = expm_cache)
    }
    get(key, envir = expm_cache, inherits = FALSE)
  }
  for (i in seq_along(time_grid)) {
    target <- time_grid[[i]]
    if (target < t_now - 1e-12) stop("time_grid must be sorted")
    delta <- target - t_now
    if (delta > 1e-12) {
      x <- as.numeric(get_step_expm(delta) %*% x)
      scale <- max(abs(x), na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) {
        x[] <- NA_real_
      } else {
        x <- x / scale
      }
      t_now <- target
    }
    states[[i]] <- vf2_normalize_state(x)
  }
  states
}

vf2_euler_states <- function(M, init, time_grid, dt) {
  if (!is.finite(dt) || dt <= 0) stop("dt must be positive")
  time_grid <- sort(unique(as.numeric(time_grid)))
  x <- as.numeric(init)
  t_now <- 0
  states <- vector("list", length(time_grid))
  for (i in seq_along(time_grid)) {
    target <- time_grid[[i]]
    if (target < t_now - 1e-12) stop("time_grid must be sorted")
    while (t_now < target - 1e-12) {
      h <- min(dt, target - t_now)
      x <- x + h * as.numeric(M %*% x)
      scale <- max(abs(x), na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) {
        x[] <- NA_real_
        break
      }
      x <- x / scale
      t_now <- t_now + h
    }
    states[[i]] <- vf2_normalize_state(x)
  }
  states
}

vf2_compare_states <- function(eigen_states, expm_states, euler_states, ngrid, n_unit, time_grid) {
  rows <- vector("list", length(time_grid))
  for (i in seq_along(time_grid)) {
    e_state <- eigen_states[[i]]
    x_state <- expm_states[[i]]
    u_state <- euler_states[[i]]
    em <- vf2_state_metrics(e_state, ngrid, n_unit)
    xm <- vf2_state_metrics(x_state, ngrid, n_unit)
    um <- vf2_state_metrics(u_state, ngrid, n_unit)
    rows[[i]] <- data.frame(
      day = time_grid[[i]],
      eigen_mean_N = em$mean_N,
      expm_mean_N = xm$mean_N,
      euler_mean_N = um$mean_N,
      diff_euler_expm_mean_N = um$mean_N - xm$mean_N,
      abs_diff_euler_expm_mean_N = abs(um$mean_N - xm$mean_N),
      diff_eigen_expm_mean_N = em$mean_N - xm$mean_N,
      abs_diff_eigen_expm_mean_N = abs(em$mean_N - xm$mean_N),
      eigen_mean_ploidy = em$mean_ploidy,
      expm_mean_ploidy = xm$mean_ploidy,
      euler_mean_ploidy = um$mean_ploidy,
      diff_euler_expm_mean_ploidy = um$mean_ploidy - xm$mean_ploidy,
      abs_diff_euler_expm_mean_ploidy = abs(um$mean_ploidy - xm$mean_ploidy),
      diff_eigen_expm_mean_ploidy = em$mean_ploidy - xm$mean_ploidy,
      abs_diff_eigen_expm_mean_ploidy = abs(em$mean_ploidy - xm$mean_ploidy),
      eigen_sd_ploidy = em$sd_ploidy,
      expm_sd_ploidy = xm$sd_ploidy,
      euler_sd_ploidy = um$sd_ploidy,
      diff_euler_expm_sd_ploidy = um$sd_ploidy - xm$sd_ploidy,
      abs_diff_euler_expm_sd_ploidy = abs(um$sd_ploidy - xm$sd_ploidy),
      diff_eigen_expm_sd_ploidy = em$sd_ploidy - xm$sd_ploidy,
      abs_diff_eigen_expm_sd_ploidy = abs(em$sd_ploidy - xm$sd_ploidy),
      eigen_fraction_N_le_25 = em$fraction_N_le_25,
      expm_fraction_N_le_25 = xm$fraction_N_le_25,
      euler_fraction_N_le_25 = um$fraction_N_le_25,
      eigen_fraction_N_below_44 = em$fraction_N_below_44,
      expm_fraction_N_below_44 = xm$fraction_N_below_44,
      euler_fraction_N_below_44 = um$fraction_N_below_44,
      eigen_fraction_N_ge_44 = em$fraction_N_ge_44,
      expm_fraction_N_ge_44 = xm$fraction_N_ge_44,
      euler_fraction_N_ge_44 = um$fraction_N_ge_44,
      eigen_fraction_N_ge_66 = em$fraction_N_ge_66,
      expm_fraction_N_ge_66 = xm$fraction_N_ge_66,
      euler_fraction_N_ge_66 = um$fraction_N_ge_66,
      eigen_fraction_N_ge_88 = em$fraction_N_ge_88,
      expm_fraction_N_ge_88 = xm$fraction_N_ge_88,
      euler_fraction_N_ge_88 = um$fraction_N_ge_88,
      max_abs_state_diff_euler_expm = max(abs(u_state - x_state), na.rm = TRUE),
      l1_state_diff_euler_expm = sum(abs(u_state - x_state), na.rm = TRUE),
      max_abs_state_diff_eigen_expm = max(abs(e_state - x_state), na.rm = TRUE),
      l1_state_diff_eigen_expm = sum(abs(e_state - x_state), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

vf2_error_summary_one <- function(comp) {
  terminal <- comp[which.max(comp$day), , drop = FALSE]
  data.frame(
    n_timepoint = nrow(comp),
    max_abs_diff_euler_expm_mean_ploidy = max(comp$abs_diff_euler_expm_mean_ploidy, na.rm = TRUE),
    terminal_abs_diff_euler_expm_mean_ploidy = terminal$abs_diff_euler_expm_mean_ploidy[[1]],
    max_abs_diff_eigen_expm_mean_ploidy = max(comp$abs_diff_eigen_expm_mean_ploidy, na.rm = TRUE),
    terminal_abs_diff_eigen_expm_mean_ploidy = terminal$abs_diff_eigen_expm_mean_ploidy[[1]],
    max_abs_diff_euler_expm_sd_ploidy = max(comp$abs_diff_euler_expm_sd_ploidy, na.rm = TRUE),
    terminal_abs_diff_euler_expm_sd_ploidy = terminal$abs_diff_euler_expm_sd_ploidy[[1]],
    max_abs_diff_eigen_expm_sd_ploidy = max(comp$abs_diff_eigen_expm_sd_ploidy, na.rm = TRUE),
    terminal_abs_diff_eigen_expm_sd_ploidy = terminal$abs_diff_eigen_expm_sd_ploidy[[1]],
    max_abs_diff_euler_expm_mean_N = max(comp$abs_diff_euler_expm_mean_N, na.rm = TRUE),
    terminal_abs_diff_euler_expm_mean_N = terminal$abs_diff_euler_expm_mean_N[[1]],
    max_abs_state_diff_euler_expm = max(comp$max_abs_state_diff_euler_expm, na.rm = TRUE),
    terminal_max_abs_state_diff_euler_expm = terminal$max_abs_state_diff_euler_expm[[1]],
    max_l1_state_diff_euler_expm = max(comp$l1_state_diff_euler_expm, na.rm = TRUE),
    terminal_l1_state_diff_euler_expm = terminal$l1_state_diff_euler_expm[[1]],
    stringsAsFactors = FALSE
  )
}

vf2_read_simulation_state_mean <- function(path) {
  if (!file.exists(path)) stop("Simulation state file not found: ", path)
  header_con <- gzfile(path, open = "rt")
  hdr <- strsplit(readLines(header_con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  close(header_con)
  needed <- c("simulation_id", "day", "N", "ploidy", "status", "cell_count")
  missing <- setdiff(needed, hdr)
  if (length(missing)) stop("Simulation state file missing column(s): ", paste(missing, collapse = ", "), " in ", path)
  classes <- rep("NULL", length(hdr))
  names(classes) <- hdr
  classes["simulation_id"] <- "integer"
  classes["day"] <- "numeric"
  classes["N"] <- "integer"
  classes["ploidy"] <- "numeric"
  classes["status"] <- "character"
  classes["cell_count"] <- "numeric"
  data_con <- gzfile(path, open = "rt")
  on.exit(close(data_con), add = TRUE)
  tab <- utils::read.delim(data_con, check.names = FALSE, stringsAsFactors = FALSE, colClasses = unname(classes))
  live <- tab[tab$status == "live", , drop = FALSE]
  if (!nrow(live)) stop("No live rows found in simulation state file: ", path)
  split_live <- split(live, live$day)
  rows <- lapply(split_live, function(sub) {
    total <- sum(sub$cell_count, na.rm = TRUE)
    mean_ploidy <- sum(sub$ploidy * sub$cell_count, na.rm = TRUE) / total
    second_ploidy <- sum(sub$ploidy^2 * sub$cell_count, na.rm = TRUE) / total
    data.frame(
      day = as.numeric(sub$day[[1]]),
      simulation_mean_N = sum(sub$N * sub$cell_count, na.rm = TRUE) / total,
      simulation_mean_ploidy = mean_ploidy,
      simulation_sd_ploidy = sqrt(max(0, second_ploidy - mean_ploidy^2)),
      simulation_live_cells = total,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$day), , drop = FALSE]
}

vf2_read_simulation_trajectories <- function(sim_root, simulation_mode, seed_ids, seed_label_map, seed_mode_map, o2_grid, init_specs, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    for (O2 in o2_grid) {
      o2_tag <- vf2_num_path_tag(O2)
      for (j in seq_len(nrow(init_specs))) {
        ploidy_tag <- vf2_num_path_tag(init_specs$initial_ploidy[[j]])
        for (sim_id in simulation_ids) {
          state_path <- file.path(
            sim_root,
            simulation_mode,
            paste0("O2_", o2_tag),
            paste0("ploidy", ploidy_tag),
            seed_id,
            paste0("sim", sim_id),
            "state_trajectory.tsv.gz"
          )
          sim <- vf2_read_simulation_state_mean(state_path)
          sim$seed_id <- seed_id
          sim$seed_label <- seed_label_map[[seed_id]]
          sim$seed_mode <- seed_mode_map[[seed_id]]
          sim$O2_pct <- O2
          sim$initial_condition <- init_specs$initial_condition[[j]]
          sim$initial_ploidy <- init_specs$initial_ploidy[[j]]
          sim$simulation_id <- sim_id
          sim$state_file <- state_path
          k <- k + 1L
          rows[[k]] <- sim[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy",
            "simulation_id", "day", "simulation_mean_N", "simulation_mean_ploidy", "simulation_sd_ploidy",
            "simulation_live_cells", "state_file"
          )]
        }
      }
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

vf2_simulation_summary <- function(sim_traj) {
  if (is.null(sim_traj) || !nrow(sim_traj)) return(data.frame())
  keys <- c("seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "day")
  split_key <- interaction(sim_traj[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sim_traj, split_key), function(sub) {
    data.frame(
      seed_id = sub$seed_id[[1]],
      seed_mode = sub$seed_mode[[1]],
      seed_label = sub$seed_label[[1]],
      O2_pct = sub$O2_pct[[1]],
      initial_condition = sub$initial_condition[[1]],
      initial_ploidy = sub$initial_ploidy[[1]],
      day = sub$day[[1]],
      simulation_n = length(unique(sub$simulation_id)),
      simulation_mean_mean_N = mean(sub$simulation_mean_N, na.rm = TRUE),
      simulation_sd_mean_N = stats::sd(sub$simulation_mean_N, na.rm = TRUE),
      simulation_min_mean_N = min(sub$simulation_mean_N, na.rm = TRUE),
      simulation_max_mean_N = max(sub$simulation_mean_N, na.rm = TRUE),
      simulation_mean_mean_ploidy = mean(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_mean_ploidy = stats::sd(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_min_mean_ploidy = min(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_max_mean_ploidy = max(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(sub$simulation_sd_ploidy, na.rm = TRUE),
      simulation_sd_sd_ploidy = stats::sd(sub$simulation_sd_ploidy, na.rm = TRUE),
      simulation_min_sd_ploidy = min(sub$simulation_sd_ploidy, na.rm = TRUE),
      simulation_max_sd_ploidy = max(sub$simulation_sd_ploidy, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$seed_id, out$O2_pct, out$initial_condition, out$day), , drop = FALSE]
}

vf2_compare_solution_to_simulation <- function(solution_traj, sim_summary, solution_name) {
  if (is.null(sim_summary) || !nrow(sim_summary)) return(data.frame())
  mean_N_col <- paste0(solution_name, "_mean_N")
  mean_ploidy_col <- paste0(solution_name, "_mean_ploidy")
  sd_ploidy_col <- paste0(solution_name, "_sd_ploidy")
  if (!all(c(mean_N_col, mean_ploidy_col, sd_ploidy_col) %in% names(solution_traj))) {
    stop("solution_traj missing columns for method: ", solution_name)
  }
  sol <- solution_traj[, c(
    "seed_id", "seed_label", "O2_pct", "initial_condition", "day",
    mean_N_col, mean_ploidy_col, sd_ploidy_col
  ), drop = FALSE]
  names(sol)[names(sol) == mean_N_col] <- "analytical_mean_N"
  names(sol)[names(sol) == mean_ploidy_col] <- "analytical_mean_ploidy"
  names(sol)[names(sol) == sd_ploidy_col] <- "analytical_sd_ploidy"
  comp <- merge(
    sim_summary,
    sol,
    by = c("seed_id", "seed_label", "O2_pct", "initial_condition", "day"),
    all = FALSE,
    sort = FALSE
  )
  comp$solution_method <- solution_name
  comp$diff_simulation_analytical_mean_ploidy <- comp$simulation_mean_mean_ploidy - comp$analytical_mean_ploidy
  comp$abs_diff_simulation_analytical_mean_ploidy <- abs(comp$diff_simulation_analytical_mean_ploidy)
  comp$diff_simulation_analytical_sd_ploidy <- comp$simulation_mean_sd_ploidy - comp$analytical_sd_ploidy
  comp$abs_diff_simulation_analytical_sd_ploidy <- abs(comp$diff_simulation_analytical_sd_ploidy)
  comp$diff_simulation_analytical_mean_N <- comp$simulation_mean_mean_N - comp$analytical_mean_N
  comp$abs_diff_simulation_analytical_mean_N <- abs(comp$diff_simulation_analytical_mean_N)
  comp
}

vf2_simulation_comparison_summary <- function(comp) {
  if (is.null(comp) || !nrow(comp)) return(data.frame())
  split_key <- interaction(comp$seed_id, comp$O2_pct, comp$initial_condition, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(comp, split_key), function(sub) {
    terminal <- sub[which.max(sub$day), , drop = FALSE]
    sd_vals <- sub$simulation_sd_mean_ploidy
    max_sd <- if (all(is.na(sd_vals))) 0 else max(sd_vals, na.rm = TRUE)
    data.frame(
      seed_id = sub$seed_id[[1]],
      seed_label = sub$seed_label[[1]],
      O2_pct = sub$O2_pct[[1]],
      initial_condition = sub$initial_condition[[1]],
      solution_method = sub$solution_method[[1]],
      n_timepoint = nrow(sub),
      simulation_n = max(sub$simulation_n, na.rm = TRUE),
      max_abs_diff_simulation_analytical_mean_ploidy = max(sub$abs_diff_simulation_analytical_mean_ploidy, na.rm = TRUE),
      terminal_abs_diff_simulation_analytical_mean_ploidy = terminal$abs_diff_simulation_analytical_mean_ploidy[[1]],
      max_abs_diff_simulation_analytical_sd_ploidy = max(sub$abs_diff_simulation_analytical_sd_ploidy, na.rm = TRUE),
      terminal_abs_diff_simulation_analytical_sd_ploidy = terminal$abs_diff_simulation_analytical_sd_ploidy[[1]],
      max_abs_diff_simulation_analytical_mean_N = max(sub$abs_diff_simulation_analytical_mean_N, na.rm = TRUE),
      terminal_abs_diff_simulation_analytical_mean_N = terminal$abs_diff_simulation_analytical_mean_N[[1]],
      max_simulation_sd_mean_ploidy = max_sd,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

vf2_plot <- function(traj, out_path, plot_dt) {
  d <- traj[abs(traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  seed_use <- unique(d$seed_id)
  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  pal <- c("#0072B2", "#D55E00", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  eigen_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(d$seed_label[d$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init <- d[d$initial_condition == init, , drop = FALSE]
    yr <- range(c(sub_init$expm_mean_ploidy, sub_init$euler_mean_ploidy), na.rm = TRUE)
    xr <- range(sub_init$day, na.rm = TRUE)
    for (seed in seed_use) {
      sub_seed <- sub_init[sub_init$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Day",
        ylab = "Mean ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed[sub_seed$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub$expm_mean_ploidy, col = eigen_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub <- sub_seed[sub_seed$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub$euler_mean_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.6, lty = 2)
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), "Expm analytical", paste0("Euler dt=", plot_dt)),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    "Fixed-O2 expm analytical vs Euler integration",
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}

vf2_plot_solution_vs_simulation <- function(solution_traj, simulation_traj, out_path, plot_dt, solution_name, solution_label, title) {
  d <- solution_traj[abs(solution_traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d) || is.null(simulation_traj) || !nrow(simulation_traj)) return(invisible(FALSE))
  solution_col <- paste0(solution_name, "_mean_ploidy")
  if (!solution_col %in% names(d)) stop("solution_traj missing column: ", solution_col)
  seed_use <- unique(d$seed_id)
  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  pal <- c("#0072B2", "#D55E00", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  analytical_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(d$seed_label[d$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init_sol <- d[d$initial_condition == init, , drop = FALSE]
    sub_init_sim <- simulation_traj[simulation_traj$initial_condition == init, , drop = FALSE]
    yr <- range(c(sub_init_sol[[solution_col]], sub_init_sim$simulation_mean_ploidy), na.rm = TRUE)
    xr <- range(c(sub_init_sol$day, sub_init_sim$day), na.rm = TRUE)
    for (seed in seed_use) {
      sub_seed_sol <- sub_init_sol[sub_init_sol$seed_id == seed, , drop = FALSE]
      sub_seed_sim <- sub_init_sim[sub_init_sim$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Day",
        ylab = "Mean ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed_sol[sub_seed_sol$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub[[solution_col]], col = analytical_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub_o2_sim <- sub_seed_sim[sub_seed_sim$O2_pct == O2, , drop = FALSE]
        reps <- split(sub_o2_sim, sub_o2_sim$simulation_id)
        for (rep_dat in reps) {
          rep_dat <- rep_dat[order(rep_dat$day), , drop = FALSE]
          lines(rep_dat$day, rep_dat$simulation_mean_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.2, lty = 2)
        }
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), solution_label, "Simulation reps"),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    title,
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}

vf2_range_pad <- function(x, frac = 0.04) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  xr <- range(x)
  span <- diff(xr)
  if (!is.finite(span) || span <= 0) {
    pad <- max(abs(xr[[1]]) * frac, 0.05)
  } else {
    pad <- span * frac
  }
  xr + c(-pad, pad)
}

vf2_solution_phase_plane <- function(solution_traj, plot_dt, solution_name, seed_mode_map) {
  d <- solution_traj[abs(solution_traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  mean_col <- paste0(solution_name, "_mean_ploidy")
  sd_col <- paste0(solution_name, "_sd_ploidy")
  if (!all(c(mean_col, sd_col) %in% names(d))) stop("solution_traj missing phase-plane columns for method: ", solution_name)
  out <- data.frame(
    seed_id = d$seed_id,
    seed_mode = unname(seed_mode_map[d$seed_id]),
    seed_label = d$seed_label,
    initial_ploidy = d$initial_ploidy,
    initial_condition = d$initial_condition,
    o2_pct = d$O2_pct,
    simulation_id = NA_integer_,
    solution_method = solution_name,
    time_days = d$day,
    mean_ploidy = d[[mean_col]],
    sd_ploidy = d[[sd_col]],
    dt = d$dt,
    stringsAsFactors = FALSE
  )
  seed_rank <- match(out$seed_id, names(seed_mode_map))
  out[order(seed_rank, out$initial_ploidy, out$o2_pct, out$time_days), , drop = FALSE]
}

vf2_simulation_phase_plane <- function(simulation_traj, seed_mode_map) {
  if (is.null(simulation_traj) || !nrow(simulation_traj)) return(data.frame())
  out <- data.frame(
    seed_id = simulation_traj$seed_id,
    seed_mode = if ("seed_mode" %in% names(simulation_traj)) simulation_traj$seed_mode else unname(seed_mode_map[simulation_traj$seed_id]),
    seed_label = simulation_traj$seed_label,
    initial_ploidy = simulation_traj$initial_ploidy,
    initial_condition = simulation_traj$initial_condition,
    o2_pct = simulation_traj$O2_pct,
    simulation_id = simulation_traj$simulation_id,
    solution_method = "simulation",
    time_days = simulation_traj$day,
    mean_ploidy = simulation_traj$simulation_mean_ploidy,
    sd_ploidy = simulation_traj$simulation_sd_ploidy,
    dt = NA_real_,
    stringsAsFactors = FALSE
  )
  seed_rank <- match(out$seed_id, names(seed_mode_map))
  out[order(seed_rank, out$initial_ploidy, out$o2_pct, out$simulation_id, out$time_days), , drop = FALSE]
}

vf2_plot_phase_plane_solution_vs_simulation <- function(solution_phase, simulation_phase, out_path, solution_label, title) {
  if (is.null(solution_phase) || !nrow(solution_phase) || is.null(simulation_phase) || !nrow(simulation_phase)) {
    return(invisible(FALSE))
  }
  seed_use <- unique(solution_phase$seed_id)
  o2_use <- sort(unique(solution_phase$o2_pct))
  init_use <- unique(solution_phase$initial_condition)
  pal <- c("#0072B2", "#D55E00", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  analytical_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(solution_phase$seed_label[solution_phase$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4.4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init_sol <- solution_phase[solution_phase$initial_condition == init, , drop = FALSE]
    sub_init_sim <- simulation_phase[simulation_phase$initial_condition == init, , drop = FALSE]
    xr <- vf2_range_pad(c(sub_init_sol$mean_ploidy, sub_init_sim$mean_ploidy))
    yr <- vf2_range_pad(c(sub_init_sol$sd_ploidy, sub_init_sim$sd_ploidy))
    for (seed in seed_use) {
      sub_seed_sol <- sub_init_sol[sub_init_sol$seed_id == seed, , drop = FALSE]
      sub_seed_sim <- sub_init_sim[sub_init_sim$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Mean ploidy",
        ylab = "SD ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed_sol[sub_seed_sol$o2_pct == O2, , drop = FALSE]
        sub <- sub[order(sub$time_days), , drop = FALSE]
        lines(sub$mean_ploidy, sub$sd_ploidy, col = analytical_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub_o2_sim <- sub_seed_sim[sub_seed_sim$o2_pct == O2, , drop = FALSE]
        reps <- split(sub_o2_sim, sub_o2_sim$simulation_id)
        for (rep_dat in reps) {
          rep_dat <- rep_dat[order(rep_dat$time_days), , drop = FALSE]
          lines(rep_dat$mean_ploidy, rep_dat$sd_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.1, lty = 2)
        }
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), solution_label, "Simulation reps"),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    title,
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}


# Unified FixO2 in vivo workflow.
fixo2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

fixo2_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

fixo2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

fixo2_default_run_dir <- function(repo_root) {
  default_hpc_run_dir <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed"
  local_run_dir <- file.path(repo_root, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
  if (dir.exists(default_hpc_run_dir)) default_hpc_run_dir else local_run_dir
}

fixo2_default_label_file <- function(repo_root) {
  ""
}

fixo2_resolve_repo_path <- function(path, repo_root, mustWork = FALSE) {
  if (is.null(path) || !length(path) || is.na(path[[1]]) || !nzchar(as.character(path[[1]]))) return(path)
  path <- as.character(path[[1]])
  if (identical(path, "~")) {
    path <- repo_root
  } else if (startsWith(path, "~/")) {
    path <- file.path(repo_root, substring(path, 3L))
  } else if (!grepl("^/", path)) {
    path <- file.path(repo_root, path)
  }
  normalizePath(path, mustWork = mustWork)
}

fixo2_default_simulation_dir <- function(repo_root) {
  file.path(repo_root, "oxygen", "results", "O2_fixed_simulation")
}

fixo2_default_html_report_script <- function(repo_root) {
  file.path(repo_root, "oxygen", "code", "O2_supply_demand_MAP", "report", "render_fixo2_invivo_report.R")
}

fixo2_prepare_dirs <- function(out_dir) {
  dirs <- c(
    out_dir,
    file.path(out_dir, "tables"),
    file.path(out_dir, "logs"),
    file.path(out_dir, "report"),
    file.path(out_dir, "attractors"),
    file.path(out_dir, "counterfactual_trajectories"),
    file.path(out_dir, "simulation"),
    file.path(out_dir, "simulation", "analytical_simulation_agreement")
  )
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

fixo2_write_output_manifest <- function(out_dir) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) {
    manifest <- data.frame()
  } else {
    rel <- substring(files, nchar(normalizePath(out_dir, mustWork = FALSE)) + 2L)
    manifest <- data.frame(
      section = ifelse(grepl("^attractors/", rel), "attractors",
        ifelse(grepl("^counterfactual_trajectories/", rel), "counterfactual_trajectories",
          ifelse(grepl("^simulation/analytical_simulation_agreement/", rel), "analytical_simulation_agreement",
            ifelse(grepl("^simulation/", rel), "simulation", "top_level")))),
      relative_path = rel,
      absolute_path = normalizePath(files, mustWork = FALSE),
      size_bytes = as.numeric(file.info(files)$size),
      stringsAsFactors = FALSE
    )
  }
  fixo2_write_tsv(manifest, file.path(out_dir, "tables", "FixO2_invivo_output_manifest.tsv"))
}

fixo2_render_html_report <- function(args, repo_root, out_dir) {
  render_html <- o2ipa_as_bool(fixo2_arg_value(args, "render_html_report", "generate_html_report", TRUE), TRUE)
  if (!isTRUE(render_html)) {
    message("Skipping FixO2 HTML report rendering because --render_html_report=FALSE.")
    return(invisible(NULL))
  }

  render_script <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "html_report_script", "report_script", fixo2_default_html_report_script(repo_root)),
      fixo2_default_html_report_script(repo_root)
    ),
    repo_root,
    mustWork = FALSE
  )
  if (!file.exists(render_script)) {
    stop("FixO2 HTML report script does not exist: ", render_script, call. = FALSE)
  }

  html_out_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "html_report_dir", "report_out_dir", file.path(out_dir, "html_report")), file.path(out_dir, "html_report")),
    repo_root,
    mustWork = FALSE
  )
  report_basename <- o2ipa_as_chr(fixo2_arg_value(args, "html_report_basename", "report_basename", "index"), "index")
  rscript <- file.path(R.home("bin"), "Rscript")
  if (!file.exists(rscript)) rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) stop("Rscript executable was not found for FixO2 HTML report rendering.", call. = FALSE)

  cmd_args <- c(
    render_script,
    paste0("--analysis_dir=", out_dir),
    paste0("--out_dir=", html_out_dir),
    paste0("--report_basename=", report_basename)
  )
  message("Rendering FixO2 HTML report: ", file.path(html_out_dir, paste0(report_basename, ".html")))
  status <- system2(rscript, cmd_args)
  if (!identical(as.integer(status), 0L)) {
    stop("FixO2 HTML report rendering failed with exit status ", status, call. = FALSE)
  }
  out_path <- file.path(html_out_dir, paste0(report_basename, ".html"))
  if (!file.exists(out_path) || is.na(file.info(out_path)$size) || file.info(out_path)$size <= 0) {
    stop("FixO2 HTML report renderer completed but did not create a non-empty report: ", out_path, call. = FALSE)
  }
  message("Completed FixO2 HTML report: ", normalizePath(out_path, mustWork = FALSE))
  invisible(out_path)
}

fixo2_normalize_parts <- function(x) {
  parts <- tolower(o2ipa_split_csv(x, "all"))
  if ("all" %in% parts) return(c("attractors", "counterfactual_trajectories", "simulation", "analytical_simulation_agreement"))
  parts <- sub("^counterfactual$", "counterfactual_trajectories", parts)
  parts <- sub("^trajectories$", "counterfactual_trajectories", parts)
  parts <- sub("^validation$", "simulation", parts)
  parts <- sub("^agreement$", "analytical_simulation_agreement", parts)
  parts <- sub("^analytical_agreement$", "analytical_simulation_agreement", parts)
  parts <- sub("^analytical_simulation$", "analytical_simulation_agreement", parts)
  parts <- sub("^scatter$", "analytical_simulation_agreement", parts)
  parts <- sub("^scatters$", "analytical_simulation_agreement", parts)
  valid <- c("attractors", "counterfactual_trajectories", "simulation", "analytical_simulation_agreement")
  bad <- setdiff(parts, valid)
  if (length(bad)) stop("Unknown run_parts value(s): ", paste(bad, collapse = ","))
  unique(parts)
}

fixo2_run_attractors <- function(args, repo_root, run_dir, label_file, out_dir) {
  if (!dir.exists(run_dir)) stop("run_dir does not exist: ", run_dir)
  fo2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "FixO2_invivo_attractors.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  objective_source <- o2ipa_as_chr(args$objective_source, "auto")
  o2_grid <- fixo2_attractor_o2_grid(args)
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, o2_grid)
  mode_o2_grid <- o2_grid
  max_seeds <- o2ipa_as_int(args$max_seeds, NA_integer_)
  generate_figures <- o2ipa_as_bool(args$generate_figures, TRUE)
  random_seed <- o2ipa_as_int(args$random_seed, 20260623L)
  set.seed(random_seed)

  run_args <- data.frame(
    argument = c("run_dir", "out_dir", "optional_label_file", "mode_source", "mode_rule", "mode_reference_o2", "objective_source", "attractor_o2_grid", "mode_o2_grid", "max_seeds", "generate_figures", "random_seed"),
    value = c(
      run_dir, out_dir, label_file,
      "fixed_o2_attractor_dominant_ploidy",
      "dominant_mean_ploidy >= 2 => mode1; dominant_mean_ploidy < 2 => mode2",
      mode_reference_o2, objective_source, paste(o2_grid, collapse = ","), paste(mode_o2_grid, collapse = ","),
      max_seeds, generate_figures, random_seed
    ),
    stringsAsFactors = FALSE
  )
  fo2_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  message("optional_label_file: ", if (nzchar(label_file)) label_file else "<none>")
  message("mode_reference_o2: ", mode_reference_o2)

  message("Collecting seed inputs.")
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = objective_source)
  if (!is.na(max_seeds) && max_seeds > 0L) {
    valid <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
    valid <- valid[order(o2ipa_seed_number(valid$seed_id)), , drop = FALSE]
    keep <- valid$seed_id[seq_len(min(max_seeds, nrow(valid)))]
    inputs$manifest <- inputs$manifest[inputs$manifest$seed_id %in% keep, , drop = FALSE]
    inputs$params_long <- inputs$params_long[inputs$params_long$seed_id %in% keep, , drop = FALSE]
    message("Using max_seeds subset: ", paste(keep, collapse = ", "))
  }
  extras <- fo2_seed_manifest_extras(inputs$manifest)
  manifest <- merge(inputs$manifest, extras, by = "seed_id", all.x = TRUE, sort = FALSE)
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  optional_metadata <- fo2_load_optional_seed_metadata(label_file, manifest)
  if (nrow(optional_metadata)) {
    manifest <- merge(manifest, optional_metadata, by = "seed_id", all.x = TRUE, sort = FALSE)
  }
  fo2_write_tsv(manifest, file.path(out_dir, "tables", "fixo2_seed_metadata.tsv"))
  fo2_write_tsv(manifest, file.path(out_dir, "tables", "seed_manifest_with_labels.tsv"))
  fo2_write_tsv(inputs$params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))

  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  cfg_meta <- o2pr_cfg_metadata(cfg)
  fo2_write_tsv(cfg_meta, file.path(out_dir, "tables", "fit_config_schema_summary.tsv"))
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  model_env <- o2ipa_source_model(SCRIPT_DIR)

  message("Computing fixed-O2 attractors.")
  rows <- list()
  valid_manifest <- manifest[manifest$fit_success %in% TRUE & manifest$seed_id %in% rownames(param_mat), , drop = FALSE]
  valid_manifest <- valid_manifest[order(o2ipa_seed_number(valid_manifest$seed_id)), , drop = FALSE]
  counter <- 0L
  for (i in seq_len(nrow(valid_manifest))) {
    seed <- valid_manifest$seed_id[[i]]
    if (i %% 25L == 0L || i == 1L) message("Processing ", seed, " (", i, "/", nrow(valid_manifest), ")")
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    for (O2 in mode_o2_grid) {
      counter <- counter + 1L
      rows[[counter]] <- fo2_dominant_attractor_one(seed, run_params, model_env, cfg, O2)
      rows[[counter]]$in_attractor_o2_grid <- any(abs(O2 - o2_grid) < 1e-9)
      rows[[counter]]$is_mode_reference_o2 <- abs(O2 - mode_reference_o2) < 1e-9
    }
  }
  all_mode_attractors <- do.call(rbind, rows)
  all_mode_attractors <- merge(all_mode_attractors, manifest[, intersect(c(
    "seed_id", "objective", "delta_objective", "mean_terminal_ploidy",
    "mean_late_drop_amplitude", "source_trajectory_regime",
    "source_mode_label", "source_mode_reason"
  ), names(manifest)), drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  all_mode_attractors <- fixo2_assign_attractor_modes(all_mode_attractors, "dominant_mean_ploidy")
  mode_by_seed_o2 <- fixo2_attractor_mode_table(all_mode_attractors)
  mode_reference_by_seed <- fixo2_reference_mode_table(mode_by_seed_o2, mode_reference_o2)
  fo2_write_tsv(mode_by_seed_o2, file.path(out_dir, "tables", "fixed_o2_attractor_mode_by_seed_o2.tsv"))
  fo2_write_tsv(mode_reference_by_seed, file.path(out_dir, "tables", "fixed_o2_attractor_mode_reference_by_seed.tsv"))
  fo2_write_tsv(mode_reference_by_seed, file.path(out_dir, "tables", "fixed_o2_attractor_mode_by_seed.tsv"))
  fo2_write_tsv(fixo2_attractor_mode_summary_by_seed(mode_by_seed_o2), file.path(out_dir, "tables", "fixed_o2_attractor_mode_summary_by_seed.tsv"))

  attractors <- all_mode_attractors[all_mode_attractors$in_attractor_o2_grid %in% TRUE, , drop = FALSE]
  attractors <- fixo2_assign_attractor_modes(attractors, "dominant_mean_ploidy")
  fo2_write_tsv(attractors, file.path(out_dir, "tables", "fixed_o2_attractors_by_seed.tsv"))
  fo2_write_tsv(fo2_mode_seed_stack_table(attractors), file.path(out_dir, "tables", "fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv"))
  fo2_write_mode_comparison_tables(attractors, out_dir)
  fo2_write_spectral_gap_outputs(attractors, out_dir, generate_figures = generate_figures)

  regime_summary <- do.call(rbind, lapply(split(attractors, list(attractors$O2_pct, attractors$trajectory_regime), drop = TRUE), function(d) {
    data.frame(
      O2_pct = d$O2_pct[[1]],
      trajectory_regime = d$trajectory_regime[[1]],
      n_seed = length(unique(d$seed_id)),
      dominant_mean_ploidy_median = stats::median(d$dominant_mean_ploidy, na.rm = TRUE),
      dominant_mean_ploidy_q25 = as.numeric(stats::quantile(d$dominant_mean_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
      dominant_mean_ploidy_q75 = as.numeric(stats::quantile(d$dominant_mean_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
      fraction_N_le_25_median = stats::median(d$dominant_fraction_N_le_25, na.rm = TRUE),
      fraction_N_below_44_median = stats::median(d$dominant_fraction_N_below_44, na.rm = TRUE),
      selection_22_vs_44_median = stats::median(d$selection_22_vs_44, na.rm = TRUE),
      selection_44_vs_88_median = stats::median(d$selection_44_vs_88, na.rm = TRUE),
      spectral_gap_median = stats::median(d$spectral_gap, na.rm = TRUE),
      eigenvector_nonnegative_fraction = mean(d$eigenvector_nonnegative %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  regime_summary <- regime_summary[order(regime_summary$O2_pct, regime_summary$trajectory_regime), , drop = FALSE]
  fo2_write_tsv(regime_summary, file.path(out_dir, "tables", "fixed_o2_attractor_regime_summary.tsv"))

  tests <- fo2_attractor_regime_tests(attractors)
  fo2_write_tsv(tests, file.path(out_dir, "tables", "fixed_o2_attractor_regime_tests.tsv"))

  corr <- fo2_parameter_correlations(attractors, inputs$params_long)
  fo2_write_tsv(corr[order(corr$O2_pct, corr$metric, -corr$abs_spearman_rho), , drop = FALSE], file.path(out_dir, "tables", "parameter_attractor_correlations.tsv"))

  if (generate_figures) fo2_plot_outputs(attractors, out_dir)
  fo2_write_report(out_dir, run_dir, label_file, attractors, regime_summary, tests, corr)
  message("Completed fixed-O2 attractor analysis: ", out_dir)
  invisible(list(out_dir = out_dir, o2_grid = o2_grid, n_seed = nrow(valid_manifest)))
}

fixo2_run_counterfactual <- function(args, repo_root, run_dir, label_file, out_dir, o2_grid) {
  cf2_mkdirs(out_dir)
  time_grid <- sort(unique(cf2_as_num_vec(args$time_grid, cf2_default_time_grid())))
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  max_seeds <- o2ipa_as_int(args$max_seeds, 0L)
  generate_figures <- o2ipa_as_bool(args$generate_figures, TRUE)
  max_time <- max(time_grid, na.rm = TRUE)

  log_file <- file.path(out_dir, "logs", "FixO2_invivo_counterfactual_trajectories.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
  message("run_dir: ", run_dir)
  message("optional_label_file ignored for FixO2 mode assignment: ", if (nzchar(label_file)) label_file else "<none>")
  message("out_dir: ", out_dir)
  message("O2 grid: ", paste(o2_grid, collapse = ","))
  message("mode_reference_o2: ", mode_reference_o2)
  message("time grid length: ", length(time_grid), "; max_time=", max_time)

  seed_inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  manifest <- seed_inputs$manifest
  params_long <- seed_inputs$params_long
  if (max_seeds > 0L && nrow(manifest) > max_seeds) {
    manifest <- manifest[seq_len(max_seeds), , drop = FALSE]
    params_long <- params_long[params_long$seed_id %in% manifest$seed_id, , drop = FALSE]
  }
  cfg <- o2pr_first_seed_cfg(manifest)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  param_mat <- o2ipa_params_wide(params_long, "value")

  traj_rows <- list()
  summary_rows <- list()
  matrix_rows <- list()
  mode_reference_rows <- list()
  counter_traj <- 0L
  counter_summary <- 0L
  counter_matrix <- 0L
  counter_mode_reference <- 0L
  init_specs <- data.frame(
    initial_condition = c("init_2N", "init_4N"),
    requested_N = c(44, 88),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(manifest))) {
    seed <- manifest$seed_id[[i]]
    if (i == 1L || i %% 25L == 0L) message("Processing ", seed, " (", i, "/", nrow(manifest), ")")
    if (!seed %in% rownames(param_mat)) next
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    ref_fm <- tryCatch(cf2_fixed_matrix(model_env, cfg, run_params, mode_reference_o2), error = function(e) NULL)
    if (is.null(ref_fm)) {
      ref_dom <- data.frame(
        status = "matrix_error",
        dominant_growth_rate = NA_real_,
        second_growth_rate = NA_real_,
        spectral_gap = NA_real_,
        dominant_mean_N = NA_real_,
        dominant_mean_ploidy = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      ref_dom <- cf2_dominant_one(ref_fm$M, ref_fm$ngrid, n_unit)
    }
    mode_fields <- fixo2_mode_fields(ref_dom$dominant_mean_ploidy[[1]])
    mode_fields$mode_source <- "fixed_o2_attractor_dominant_ploidy_at_reference_o2"
    mode_fields$mode_rule <- paste0(
      "dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " >= 2 => mode1; dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " < 2 => mode2"
    )
    mode_fields$mode_reference_o2_pct <- mode_reference_o2
    mode_fields$mode_reference_o2_key <- fixo2_o2_key(mode_reference_o2)
    mode_fields$mode_reference_dominant_mean_ploidy <- ref_dom$dominant_mean_ploidy[[1]]
    mode_fields$mode_reference_status <- ref_dom$status[[1]]
    mode_fields$mode_reference_dominant_growth_rate <- ref_dom$dominant_growth_rate[[1]]
    mode_fields$mode_reference_spectral_gap <- ref_dom$spectral_gap[[1]]
    counter_mode_reference <- counter_mode_reference + 1L
    mode_reference_rows[[counter_mode_reference]] <- data.frame(seed_id = seed, mode_fields, stringsAsFactors = FALSE)
    for (O2 in o2_grid) {
      fm <- tryCatch(cf2_fixed_matrix(model_env, cfg, run_params, O2), error = function(e) NULL)
      if (is.null(fm)) next
      dom <- cf2_dominant_one(fm$M, fm$ngrid, n_unit)
      counter_matrix <- counter_matrix + 1L
      matrix_rows[[counter_matrix]] <- data.frame(seed_id = seed, O2_pct = O2, dom, mode_fields, stringsAsFactors = FALSE)
      for (j in seq_len(nrow(init_specs))) {
        init <- cf2_init_vector(fm$ngrid, init_specs$requested_N[[j]])
        sim <- cf2_eigen_trajectory(fm$M, fm$ngrid, init$vector, time_grid, n_unit)
        tr <- sim$trajectory
        if (nrow(tr)) {
          tr$seed_id <- seed
          tr$O2_pct <- O2
          tr$initial_condition <- init_specs$initial_condition[[j]]
          tr$requested_initial_N <- init_specs$requested_N[[j]]
          tr$used_initial_N <- init$used_N
          tr$status <- sim$status
          tr$trajectory_regime <- mode_fields$trajectory_regime[[1]]
          tr$mode_label <- mode_fields$mode_label[[1]]
          tr$mode_source <- mode_fields$mode_source[[1]]
          tr$mode_rule <- mode_fields$mode_rule[[1]]
          tr$mode_threshold_dominant_ploidy <- mode_fields$mode_threshold_dominant_ploidy[[1]]
          tr$mode_reference_o2_pct <- mode_fields$mode_reference_o2_pct[[1]]
          tr$mode_reference_o2_key <- mode_fields$mode_reference_o2_key[[1]]
          tr$mode_reference_dominant_mean_ploidy <- mode_fields$mode_reference_dominant_mean_ploidy[[1]]
          counter_traj <- counter_traj + 1L
          traj_rows[[counter_traj]] <- tr[, c(
            "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
            "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
            "mode_reference_dominant_mean_ploidy", "O2_pct", "initial_condition",
            "requested_initial_N", "used_initial_N", "status", "day",
            "mean_N", "mean_ploidy", "fraction_N_le_25", "fraction_N_below_44",
            "fraction_N_ge_44", "fraction_N_ge_66", "fraction_N_ge_88"
          )]
        }
        sm <- cf2_summarize_trajectory(tr, max_time)
        counter_summary <- counter_summary + 1L
        summary_rows[[counter_summary]] <- data.frame(
          seed_id = seed,
          mode_fields,
          O2_pct = O2,
          initial_condition = init_specs$initial_condition[[j]],
          requested_initial_N = init_specs$requested_N[[j]],
          used_initial_N = init$used_N,
          status = sim$status,
          sm,
          dom,
          terminal_minus_dominant_ploidy = sm$terminal_mean_ploidy - dom$dominant_mean_ploidy,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  trajectory <- if (length(traj_rows)) do.call(rbind, traj_rows) else data.frame()
  summary_by_seed <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
  dominant_by_seed <- if (length(matrix_rows)) do.call(rbind, matrix_rows) else data.frame()
  mode_reference_by_seed <- if (length(mode_reference_rows)) do.call(rbind, mode_reference_rows) else data.frame()
  counterfactual_mode_by_seed_o2 <- if (nrow(dominant_by_seed)) {
    dominant_by_seed$O2_key <- fixo2_o2_key(dominant_by_seed$O2_pct)
    dominant_by_seed[, intersect(c(
      "seed_id", "O2_pct", "O2_key", "dominant_mean_ploidy", "trajectory_regime",
      "mode_label", "mode_source", "mode_rule", "mode_threshold_dominant_ploidy",
      "mode_reference_o2_pct", "mode_reference_o2_key", "mode_reference_dominant_mean_ploidy",
      "dominant_growth_rate", "spectral_gap"
    ), names(dominant_by_seed)), drop = FALSE]
  } else {
    data.frame()
  }
  regime_summary <- cf2_regime_summary(summary_by_seed)
  tests <- cf2_regime_tests(summary_by_seed)
  correlations <- cf2_parameter_correlations(summary_by_seed, params_long)

  cf2_write_tsv(data.frame(
    argument = c("run_dir", "optional_label_file", "mode_source", "mode_rule", "mode_reference_o2", "out_dir", "o2_grid", "time_grid", "max_seeds", "generate_figures"),
    value = c(
      run_dir, label_file,
      "fixed_o2_attractor_dominant_ploidy_at_reference_o2",
      paste0(
        "dominant_mean_ploidy at fixed O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        " >= 2 => mode1; dominant_mean_ploidy at fixed O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        " < 2 => mode2"
      ),
      mode_reference_o2,
      out_dir, paste(o2_grid, collapse = ","), paste(time_grid, collapse = ","),
      as.character(max_seeds), as.character(generate_figures)
    ),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  cf2_write_tsv(trajectory, file.path(out_dir, "tables", "fixed_o2_counterfactual_trajectories.tsv"))
  cf2_write_tsv(summary_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_summary_by_seed.tsv"))
  cf2_write_tsv(dominant_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_dominant_by_seed.tsv"))
  cf2_write_tsv(counterfactual_mode_by_seed_o2, file.path(out_dir, "tables", "fixed_o2_counterfactual_mode_by_seed_o2.tsv"))
  cf2_write_tsv(mode_reference_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_mode_reference_by_seed.tsv"))
  cf2_write_tsv(regime_summary, file.path(out_dir, "tables", "fixed_o2_counterfactual_regime_summary.tsv"))
  cf2_write_tsv(tests, file.path(out_dir, "tables", "fixed_o2_counterfactual_regime_tests.tsv"))
  cf2_write_tsv(correlations, file.path(out_dir, "tables", "fixed_o2_counterfactual_parameter_correlations.tsv"))
  cf2_write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  if (isTRUE(generate_figures)) cf2_plot(trajectory, summary_by_seed, file.path(out_dir, "figures"))
  cf2_report(out_dir, args, regime_summary, tests, correlations)
  message("Completed fixed-O2 counterfactual trajectory analysis: ", out_dir)
  invisible(list(out_dir = out_dir, o2_grid = o2_grid, n_seed = nrow(manifest)))
}

fixo2_simulation_state_path <- function(simulation_dir, simulation_mode, O2, initial_ploidy, seed_id, simulation_id) {
  file.path(
    simulation_dir,
    simulation_mode,
    paste0("O2_", vf2_num_path_tag(O2)),
    paste0("ploidy", vf2_num_path_tag(initial_ploidy)),
    seed_id,
    paste0("sim", simulation_id),
    "state_trajectory.tsv.gz"
  )
}

fixo2_expected_simulation_tasks <- function(run_dir, simulation_dir, simulation_mode, seed_ids, o2_grid, init_specs, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    seed_dir <- file.path(run_dir, seed_id)
    for (O2 in o2_grid) {
      for (j in seq_len(nrow(init_specs))) {
        initial_ploidy <- init_specs$initial_ploidy[[j]]
        leaf_dir <- dirname(dirname(fixo2_simulation_state_path(simulation_dir, simulation_mode, O2, initial_ploidy, seed_id, 1L)))
        for (sim_id in simulation_ids) {
          state_file <- fixo2_simulation_state_path(simulation_dir, simulation_mode, O2, initial_ploidy, seed_id, sim_id)
          k <- k + 1L
          rows[[k]] <- data.frame(
            task_id = k,
            seed_id = seed_id,
            seed_dir = seed_dir,
            O2_pct = O2,
            initial_ploidy = initial_ploidy,
            initial_condition = init_specs$initial_condition[[j]],
            simulation_id = as.integer(sim_id),
            output_dir = dirname(state_file),
            state_file = state_file,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    info <- file.info(out$state_file)
    out$complete_before <- file.exists(out$state_file) & !is.na(info$size) & info$size > 0
  }
  out
}

fixo2_simulation_forward_args <- function(args) {
  keys <- c("dt", "report_dt", "initial_cells", "joint_scope", "Crowding", "O2_growth", "start_with", "ploidy_O2_death")
  out <- character(0)
  for (key in keys) {
    if (!is.null(args[[key]]) && length(args[[key]]) > 0L) {
      out <- c(out, paste0("--", key, "=", as.character(args[[key]][[1]])))
    }
  }
  if (!is.null(args$save_every_days) && length(args$save_every_days) > 0L) {
    out <- c(out, paste0("--save_every_days=", as.character(args$save_every_days[[1]])))
  } else if (!is.null(args$time_step_days) && length(args$time_step_days) > 0L) {
    out <- c(out, paste0("--save_every_days=", as.character(args$time_step_days[[1]])))
  }
  out
}

fixo2_simulation_n_core <- function(args, task_count = Inf) {
  raw <- args$simulation_n_core %||% args$simulation_workers %||% args$n_core
  if (is.null(raw) || !length(raw)) raw <- 1L
  n_core <- o2ipa_as_int(raw, 1L)
  if (!is.finite(n_core) || is.na(n_core) || n_core < 1L) n_core <- 1L
  if (is.finite(task_count)) n_core <- min(as.integer(n_core), as.integer(task_count))
  max(1L, as.integer(n_core))
}

fixo2_simulation_worker_threads <- function(args) {
  raw <- args$simulation_worker_threads %||% args$simulation_threads_per_worker
  if (is.null(raw) || !length(raw)) raw <- 1L
  n_threads <- o2ipa_as_int(raw, 1L)
  if (!is.finite(n_threads) || is.na(n_threads) || n_threads < 1L) n_threads <- 1L
  as.integer(n_threads)
}

fixo2_simulation_worker_env <- function(args) {
  n_threads <- fixo2_simulation_worker_threads(args)
  c(
    paste0(
      c(
        "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS"
      ),
      "=",
      n_threads
    ),
    "OMP_WAIT_POLICY=PASSIVE"
  )
}

fixo2_task_log_completed <- function(log_file) {
  if (!file.exists(log_file)) return(NA)
  lines <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) character(0))
  if (!length(lines)) return(FALSE)
  any(grepl("^Done\\.$", lines))
}

fixo2_run_missing_simulation_task <- function(idx, tasks, rscript, sim_script, simulation_mode, time_days, extra_args, out_dir, worker_env) {
  task <- tasks[idx, , drop = FALSE]
  log_file <- file.path(out_dir, "logs", sprintf("fix_o2_simulation_task_%04d.log", task$task_id[[1]]))
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  cmd_args <- c(
    sim_script,
    paste0("--fit_dir=", task$seed_dir[[1]]),
    paste0("--simulation=", simulation_mode),
    paste0("--initial_ploidy=", task$initial_ploidy[[1]]),
    paste0("--time_days=", time_days),
    "--n_sim=1",
    paste0("--simulation_id=", task$simulation_id[[1]]),
    paste0("--o2=", task$O2_pct[[1]]),
    paste0("--seed_id=", task$seed_id[[1]]),
    paste0("--out_dir=", task$output_dir[[1]]),
    extra_args
  )
  error_message <- ""
  status <- tryCatch(
    system2(rscript, cmd_args, stdout = log_file, stderr = log_file, env = worker_env),
    error = function(e) {
      error_message <<- conditionMessage(e)
      cat("ERROR launching fixed-O2 simulation task: ", error_message, "\n", file = log_file, append = TRUE, sep = "")
      1L
    }
  )
  status <- if (is.null(status)) 0L else as.integer(status)
  data.frame(
    task_index = as.integer(idx),
    task_id = as.integer(task$task_id[[1]]),
    run_status = status,
    generated = TRUE,
    log_file = log_file,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

fixo2_ensure_simulation <- function(args, run_dir, simulation_dir, simulation_mode, seed_ids, o2_grid, init_specs, simulation_ids, out_dir) {
  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
  tasks <- fixo2_expected_simulation_tasks(run_dir, simulation_dir, simulation_mode, seed_ids, o2_grid, init_specs, simulation_ids)
  if (!nrow(tasks)) return(tasks)

  generate_missing <- o2ipa_as_bool(args$generate_missing_simulation, TRUE)
  time_days <- o2ipa_as_num(args$time_days, 1000)
  sim_script <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "simulation", "fix_o2_simulation.R"), mustWork = FALSE)
  if (!file.exists(sim_script)) stop("Simulation script was not found: ", sim_script)
  rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
  tasks$run_status <- NA_integer_
  tasks$generated <- FALSE
  tasks$log_file <- file.path(out_dir, "logs", sprintf("fix_o2_simulation_task_%04d.log", tasks$task_id))
  tasks$log_completed_before <- vapply(tasks$log_file, fixo2_task_log_completed, logical(1))
  tasks$error_message <- ""
  tasks$complete_before <- tasks$complete_before & (is.na(tasks$log_completed_before) | tasks$log_completed_before)

  missing <- which(!tasks$complete_before)
  if (length(missing) && !generate_missing) {
    fixo2_write_tsv(tasks, file.path(out_dir, "tables", "fixed_o2_simulation_file_status.tsv"))
    stop("Missing ", length(missing), " fixed-O2 simulation file(s); set --generate_missing_simulation=TRUE to create them.")
  }

  if (length(missing)) {
    n_core <- fixo2_simulation_n_core(args, length(missing))
    worker_threads <- fixo2_simulation_worker_threads(args)
    worker_env <- fixo2_simulation_worker_env(args)
    message("Generating missing fixed-O2 simulation task(s): ", length(missing))
    message("simulation_n_core: ", n_core)
    message("simulation_worker_threads: ", worker_threads)
    extra_args <- fixo2_simulation_forward_args(args)
    bad_seed_dir <- missing[!dir.exists(tasks$seed_dir[missing])]
    if (length(bad_seed_dir)) {
      stop("Seed directory does not exist for simulation task(s): ", paste(tasks$task_id[bad_seed_dir], collapse = ","))
    }
    runner <- function(idx) {
      tryCatch(
        fixo2_run_missing_simulation_task(
          idx = idx,
          tasks = tasks,
          rscript = rscript,
          sim_script = sim_script,
          simulation_mode = simulation_mode,
          time_days = time_days,
          extra_args = extra_args,
          out_dir = out_dir,
          worker_env = worker_env
        ),
        error = function(e) {
          log_file <- file.path(out_dir, "logs", sprintf("fix_o2_simulation_task_%04d.log", tasks$task_id[[idx]]))
          cat("ERROR running fixed-O2 simulation task: ", conditionMessage(e), "\n", file = log_file, append = TRUE, sep = "")
          data.frame(
            task_index = as.integer(idx),
            task_id = as.integer(tasks$task_id[[idx]]),
            run_status = 1L,
            generated = FALSE,
            log_file = log_file,
            error_message = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        }
      )
    }
    if (n_core > 1L && .Platform$OS.type != "windows") {
      run_results <- parallel::mclapply(missing, runner, mc.cores = n_core, mc.preschedule = FALSE)
    } else {
      if (n_core > 1L) warning("Parallel fixed-O2 simulation generation is not supported on Windows; falling back to sequential.")
      run_results <- lapply(missing, runner)
    }
    run_results <- do.call(rbind, run_results)
    tasks$run_status[run_results$task_index] <- run_results$run_status
    tasks$generated[run_results$task_index] <- run_results$generated
    tasks$log_file[run_results$task_index] <- run_results$log_file
    tasks$error_message[run_results$task_index] <- run_results$error_message
    failed <- run_results[run_results$run_status != 0L | nzchar(run_results$error_message), , drop = FALSE]
    if (nrow(failed)) {
      fixo2_write_tsv(tasks, file.path(out_dir, "tables", "fixed_o2_simulation_file_status.tsv"))
      stop(
        "Fixed-O2 simulation task(s) failed: task_id=",
        paste(failed$task_id, collapse = ","),
        "; log=",
        paste(failed$log_file, collapse = ",")
      )
    }
  } else {
    message("All requested fixed-O2 simulation files already exist; skipping simulation generation.")
  }

  info_after <- file.info(tasks$state_file)
  tasks$complete_after <- file.exists(tasks$state_file) & !is.na(info_after$size) & info_after$size > 0
  fixo2_write_tsv(tasks, file.path(out_dir, "tables", "fixed_o2_simulation_file_status.tsv"))
  if (any(!tasks$complete_after)) {
    failed <- tasks[!tasks$complete_after, , drop = FALSE]
    stop("Missing fixed-O2 simulation output(s) after generation: ", paste(failed$task_id, collapse = ","))
  }
  invisible(tasks)
}

fixo2_run_simulation_validation <- function(args, repo_root, run_dir, out_dir, simulation_dir, o2_grid) {
  simulation_mode <- o2ipa_as_chr(args$simulation_mode, "invivo")
  simulation_ids <- vf2_as_int_vec(args$simulation_ids, seq_len(10L))
  include_simulation <- o2ipa_as_bool(args$include_simulation, TRUE)
  generate_missing <- o2ipa_as_bool(args$generate_missing_simulation, TRUE)
  simulation_n_core <- fixo2_simulation_n_core(args)
  simulation_worker_threads <- fixo2_simulation_worker_threads(args)
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  seed_selection_n_workers <- o2ipa_as_int(args$simulation_seed_selection_n_workers %||% args$seed_selection_n_workers, simulation_n_core)
  if (!is.finite(seed_selection_n_workers) || seed_selection_n_workers < 1L) seed_selection_n_workers <- 1L
  seed_selection_n_workers <- as.integer(seed_selection_n_workers)
  manual_seed_selection <- fixo2_arg_nonempty(args$seed_ids) || fixo2_arg_nonempty(args$seed_id)
  simulation_analysis_dir <- normalizePath(file.path(out_dir, ".."), mustWork = FALSE)
  seed_mode_objective_metadata <- fixo2_simulation_mode_objective_metadata(
    analysis_dir = simulation_analysis_dir,
    run_dir = run_dir,
    o2_grid = o2_grid,
    mode_reference_o2 = mode_reference_o2,
    seed_ids = NULL,
    n_workers = seed_selection_n_workers
  )
  if (manual_seed_selection) {
    seed_ids <- vf2_as_seed_vec(args$seed_ids %||% args$seed_id, character())
    if (!length(seed_ids)) stop("Manual simulation seed selection was requested, but no valid seed_ids were supplied.")
    seed_mode_arg <- args$seed_modes %||% args$seed_labels
    seed_modes <- vf2_seed_modes(seed_ids, seed_mode_arg)
    seed_selection <- fixo2_manual_simulation_seed_table(
      seed_ids = seed_ids,
      seed_modes = seed_modes,
      metadata = seed_mode_objective_metadata,
      mode_reference_o2 = mode_reference_o2
    )
  } else {
    seed_selection <- fixo2_select_best_objective_seed_by_mode(
      metadata = seed_mode_objective_metadata,
      mode_reference_o2 = mode_reference_o2
    )
    seed_ids <- seed_selection$seed_id
  }
  seed_modes <- seed_selection$mode_label
  seed_modes[is.na(seed_modes)] <- ""
  seed_labels <- seed_ids
  has_mode <- nzchar(seed_modes)
  seed_labels[has_mode] <- paste0(seed_ids[has_mode], " (", seed_modes[has_mode], ")")
  seed_mode_map <- stats::setNames(seed_modes, seed_ids)
  seed_label_map <- stats::setNames(seed_labels, seed_ids)
  seed_selection_mode <- if (nrow(seed_selection)) seed_selection$selection_mode[[1]] else NA_character_
  default_time_grid <- if (include_simulation) {
    seq(0, o2ipa_as_num(args$time_days, 1000), by = o2ipa_as_num(args$time_step_days, 1))
  } else {
    vf2_default_time_grid()
  }
  time_grid <- sort(unique(vf2_as_num_vec(args$time_grid, default_time_grid)))
  dt_grid <- sort(unique(vf2_as_num_vec(args$dt_grid %||% args$dt, c(0.05))))
  plot_dt <- o2ipa_as_num(args$plot_dt, min(dt_grid))

  vf2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "FixO2_invivo_simulation_validation.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  message("simulation_dir: ", simulation_dir)
  message("include_simulation: ", include_simulation)
  message("generate_missing_simulation: ", generate_missing)
  message("simulation_n_core: ", simulation_n_core)
  message("simulation_worker_threads: ", simulation_worker_threads)
  message("mode_reference_o2: ", mode_reference_o2)
  message("seed_selection_mode: ", seed_selection_mode)
  message("seed_selection_n_workers: ", seed_selection_n_workers)
  if (include_simulation) message("simulation_ids: ", paste(simulation_ids, collapse = ","))
  message("seed_ids: ", paste(seed_ids, collapse = ","))
  message("seed_modes: ", paste(seed_modes, collapse = ","))
  message("seed_labels: ", paste(seed_labels, collapse = ","))
  message("O2 grid: ", paste(o2_grid, collapse = ","))
  message("dt grid: ", paste(dt_grid, collapse = ","))

  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  init_specs <- data.frame(
    initial_condition = c("init_2N", "init_4N"),
    requested_N = c(44, 88),
    initial_ploidy = c(2, 4),
    stringsAsFactors = FALSE
  )
  vf2_write_tsv(seed_selection, file.path(out_dir, "tables", "fixed_o2_simulation_representative_seeds.tsv"))

  if (include_simulation) {
    fixo2_ensure_simulation(
      args = args,
      run_dir = run_dir,
      simulation_dir = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      o2_grid = o2_grid,
      init_specs = init_specs,
      simulation_ids = simulation_ids,
      out_dir = out_dir
    )
  }

  trajectory_rows <- list()
  state_error_rows <- list()
  summary_rows <- list()
  tc <- 0L
  ec <- 0L
  sc <- 0L
  for (seed_id in seed_ids) {
    message("Processing seed=", seed_id, " label=", seed_label_map[[seed_id]])
    manifest <- inputs$manifest[inputs$manifest$seed_id == seed_id, , drop = FALSE]
    if (!nrow(manifest)) stop("seed_id not found in run_dir: ", seed_id)
    params_long <- inputs$params_long[inputs$params_long$seed_id == seed_id, , drop = FALSE]
    if (!nrow(params_long)) stop("No parameters found for ", seed_id)
    cfg <- o2pr_first_seed_cfg(manifest)
    n_unit <- as.numeric(cfg$N_UNIT %||% 22)
    if (!seed_id %in% rownames(param_mat)) stop("seed_id missing from parameter matrix: ", seed_id)
    pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    for (O2 in o2_grid) {
      message("  O2=", O2)
      fm <- vf2_fixed_matrix(model_env, cfg, run_params, O2)
      for (j in seq_len(nrow(init_specs))) {
        init <- vf2_init_vector(fm$ngrid, init_specs$requested_N[[j]])
        eig <- vf2_eigen_states(fm$M, init$vector, time_grid)
        expm <- vf2_expm_states(fm$M, init$vector, time_grid)
        for (dt in dt_grid) {
          message("    ", init_specs$initial_condition[[j]], ", dt=", dt)
          eu <- vf2_euler_states(fm$M, init$vector, time_grid, dt)
          comp <- vf2_compare_states(eig$states, expm, eu, fm$ngrid, n_unit, time_grid)
          comp$seed_id <- seed_id
          comp$seed_mode <- seed_mode_map[[seed_id]]
          comp$seed_label <- seed_label_map[[seed_id]]
          comp$O2_pct <- O2
          comp$initial_condition <- init_specs$initial_condition[[j]]
          comp$requested_initial_N <- init_specs$requested_N[[j]]
          comp$initial_ploidy <- init_specs$initial_ploidy[[j]]
          comp$used_initial_N <- init$used_N
          comp$dt <- dt
          comp$dominant_growth_rate <- eig$lambda_ref
          comp <- comp[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "requested_initial_N", "used_initial_N",
            "dt", "day", "dominant_growth_rate",
            "eigen_mean_N", "expm_mean_N", "euler_mean_N",
            "diff_euler_expm_mean_N", "abs_diff_euler_expm_mean_N",
            "diff_eigen_expm_mean_N", "abs_diff_eigen_expm_mean_N",
            "eigen_mean_ploidy", "expm_mean_ploidy", "euler_mean_ploidy",
            "diff_euler_expm_mean_ploidy", "abs_diff_euler_expm_mean_ploidy",
            "diff_eigen_expm_mean_ploidy", "abs_diff_eigen_expm_mean_ploidy",
            "eigen_sd_ploidy", "expm_sd_ploidy", "euler_sd_ploidy",
            "diff_euler_expm_sd_ploidy", "abs_diff_euler_expm_sd_ploidy",
            "diff_eigen_expm_sd_ploidy", "abs_diff_eigen_expm_sd_ploidy",
            "eigen_fraction_N_le_25", "expm_fraction_N_le_25", "euler_fraction_N_le_25",
            "eigen_fraction_N_below_44", "expm_fraction_N_below_44", "euler_fraction_N_below_44",
            "eigen_fraction_N_ge_44", "expm_fraction_N_ge_44", "euler_fraction_N_ge_44",
            "eigen_fraction_N_ge_66", "expm_fraction_N_ge_66", "euler_fraction_N_ge_66",
            "eigen_fraction_N_ge_88", "expm_fraction_N_ge_88", "euler_fraction_N_ge_88",
            "max_abs_state_diff_euler_expm", "l1_state_diff_euler_expm",
            "max_abs_state_diff_eigen_expm", "l1_state_diff_eigen_expm"
          )]
          tc <- tc + 1L
          trajectory_rows[[tc]] <- comp
          serr <- comp[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "dt", "day",
            "max_abs_state_diff_euler_expm", "l1_state_diff_euler_expm",
            "max_abs_state_diff_eigen_expm", "l1_state_diff_eigen_expm"
          ), drop = FALSE]
          ec <- ec + 1L
          state_error_rows[[ec]] <- serr
          sm <- vf2_error_summary_one(comp)
          sm$seed_id <- seed_id
          sm$seed_mode <- seed_mode_map[[seed_id]]
          sm$seed_label <- seed_label_map[[seed_id]]
          sm$O2_pct <- O2
          sm$initial_condition <- init_specs$initial_condition[[j]]
          sm$initial_ploidy <- init_specs$initial_ploidy[[j]]
          sm$requested_initial_N <- init_specs$requested_N[[j]]
          sm$used_initial_N <- init$used_N
          sm$dt <- dt
          sm$dominant_growth_rate <- eig$lambda_ref
          sm <- sm[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "requested_initial_N", "used_initial_N",
            "dt", "dominant_growth_rate", "n_timepoint",
            "max_abs_diff_euler_expm_mean_ploidy", "terminal_abs_diff_euler_expm_mean_ploidy",
            "max_abs_diff_eigen_expm_mean_ploidy", "terminal_abs_diff_eigen_expm_mean_ploidy",
            "max_abs_diff_euler_expm_sd_ploidy", "terminal_abs_diff_euler_expm_sd_ploidy",
            "max_abs_diff_eigen_expm_sd_ploidy", "terminal_abs_diff_eigen_expm_sd_ploidy",
            "max_abs_diff_euler_expm_mean_N", "terminal_abs_diff_euler_expm_mean_N",
            "max_abs_state_diff_euler_expm", "terminal_max_abs_state_diff_euler_expm",
            "max_l1_state_diff_euler_expm", "terminal_l1_state_diff_euler_expm"
          )]
          sc <- sc + 1L
          summary_rows[[sc]] <- sm
        }
      }
    }
  }

  trajectories <- do.call(rbind, trajectory_rows)
  state_errors <- do.call(rbind, state_error_rows)
  summaries <- do.call(rbind, summary_rows)
  summaries <- summaries[order(match(summaries$seed_id, seed_ids), summaries$O2_pct, summaries$initial_condition, summaries$dt), , drop = FALSE]
  seed_slug <- paste(seed_ids, collapse = "_")

  simulation_traj <- data.frame()
  simulation_summary <- data.frame()
  expm_simulation_comparison <- data.frame()
  expm_simulation_comparison_summary <- data.frame()
  eigen_simulation_comparison <- data.frame()
  eigen_simulation_comparison_summary <- data.frame()
  expm_phase_plane <- data.frame()
  eigen_phase_plane <- data.frame()
  simulation_phase_plane <- data.frame()
  if (include_simulation) {
    message("Reading fixed-O2 simulation trajectories")
    simulation_traj <- vf2_read_simulation_trajectories(
      sim_root = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      seed_label_map = seed_label_map,
      seed_mode_map = seed_mode_map,
      o2_grid = o2_grid,
      init_specs = init_specs,
      simulation_ids = simulation_ids
    )
    simulation_summary <- vf2_simulation_summary(simulation_traj)
    expm_simulation_comparison <- vf2_compare_solution_to_simulation(trajectories, simulation_summary, "expm")
    expm_simulation_comparison_summary <- vf2_simulation_comparison_summary(expm_simulation_comparison)
    eigen_simulation_comparison <- vf2_compare_solution_to_simulation(trajectories, simulation_summary, "eigen")
    eigen_simulation_comparison_summary <- vf2_simulation_comparison_summary(eigen_simulation_comparison)
  }
  expm_phase_plane <- vf2_solution_phase_plane(trajectories, plot_dt, "expm", seed_mode_map)
  eigen_phase_plane <- vf2_solution_phase_plane(trajectories, plot_dt, "eigen", seed_mode_map)
  if (include_simulation) {
    simulation_phase_plane <- vf2_simulation_phase_plane(simulation_traj, seed_mode_map)
  }

  run_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "simulation_dir", "simulation_mode", "include_simulation",
      "generate_missing_simulation", "simulation_n_core", "simulation_worker_threads",
      "mode_reference_o2", "seed_selection_mode", "seed_selection_n_workers",
      "simulation_ids", "seed_ids", "seed_modes", "seed_labels", "o2_grid", "time_grid", "dt_grid", "plot_dt"
    ),
    value = c(
      run_dir, out_dir, simulation_dir, simulation_mode, as.character(include_simulation),
      as.character(generate_missing), as.character(simulation_n_core), as.character(simulation_worker_threads),
      as.character(mode_reference_o2), seed_selection_mode, as.character(seed_selection_n_workers),
      paste(simulation_ids, collapse = ","), paste(seed_ids, collapse = ","), paste(seed_modes, collapse = ","), paste(seed_labels, collapse = ","),
      paste(o2_grid, collapse = ","), paste(time_grid, collapse = ","), paste(dt_grid, collapse = ","), as.character(plot_dt)
    ),
    stringsAsFactors = FALSE
  )
  vf2_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  vf2_write_tsv(trajectories, file.path(out_dir, "tables", "eigen_vs_euler_trajectories.tsv"))
  vf2_write_tsv(state_errors, file.path(out_dir, "tables", "eigen_vs_euler_state_errors.tsv"))
  vf2_write_tsv(summaries, file.path(out_dir, "tables", "eigen_vs_euler_error_summary.tsv"))
  vf2_write_tsv(expm_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_expm_analytical.tsv"))
  vf2_write_tsv(eigen_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_eigen_analytical.tsv"))
  if (include_simulation) {
    vf2_write_tsv(simulation_traj, file.path(out_dir, "tables", "simulation_replicate_mean_ploidy_trajectories.tsv"))
    vf2_write_tsv(simulation_summary, file.path(out_dir, "tables", "simulation_summary_mean_ploidy_trajectories.tsv"))
    vf2_write_tsv(expm_simulation_comparison, file.path(out_dir, "tables", "expm_vs_simulation_mean_ploidy_comparison.tsv"))
    vf2_write_tsv(expm_simulation_comparison_summary, file.path(out_dir, "tables", "expm_vs_simulation_error_summary.tsv"))
    vf2_write_tsv(eigen_simulation_comparison, file.path(out_dir, "tables", "eigen_vs_simulation_mean_ploidy_comparison.tsv"))
    vf2_write_tsv(eigen_simulation_comparison_summary, file.path(out_dir, "tables", "eigen_vs_simulation_error_summary.tsv"))
    vf2_write_tsv(simulation_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_simulation_replicates.tsv"))
    vf2_write_tsv(rbind(expm_phase_plane, simulation_phase_plane), file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_expm_vs_simulation.tsv"))
    vf2_write_tsv(rbind(eigen_phase_plane, simulation_phase_plane), file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_eigen_vs_simulation.tsv"))
  }
  vf2_plot(trajectories, file.path(out_dir, "figures", sprintf("expm_vs_euler_mean_ploidy_%s.pdf", seed_slug)), plot_dt = plot_dt)
  if (include_simulation) {
    vf2_plot_solution_vs_simulation(
      trajectories,
      simulation_traj,
      file.path(out_dir, "figures", sprintf("expm_vs_simulation_replicates_mean_ploidy_%s.pdf", seed_slug)),
      plot_dt = plot_dt,
      solution_name = "expm",
      solution_label = "Expm analytical",
      title = "Fixed-O2 expm analytical vs simulation replicates"
    )
    vf2_plot_solution_vs_simulation(
      trajectories,
      simulation_traj,
      file.path(out_dir, "figures", sprintf("eigen_vs_simulation_replicates_mean_ploidy_%s.pdf", seed_slug)),
      plot_dt = plot_dt,
      solution_name = "eigen",
      solution_label = "Eigen analytical",
      title = "Fixed-O2 eigen analytical vs simulation replicates"
    )
    vf2_plot_phase_plane_solution_vs_simulation(
      expm_phase_plane,
      simulation_phase_plane,
      file.path(out_dir, "figures", sprintf("expm_vs_simulation_phase_plane_mean_sd_ploidy_%s.pdf", seed_slug)),
      solution_label = "Expm analytical",
      title = "Fixed-O2 phase plane: expm analytical vs simulation replicates"
    )
    vf2_plot_phase_plane_solution_vs_simulation(
      eigen_phase_plane,
      simulation_phase_plane,
      file.path(out_dir, "figures", sprintf("eigen_vs_simulation_phase_plane_mean_sd_ploidy_%s.pdf", seed_slug)),
      solution_label = "Eigen analytical",
      title = "Fixed-O2 phase plane: eigen analytical vs simulation replicates"
    )
  }
  message("Completed simulation validation: ", out_dir)
  invisible(list(out_dir = out_dir, o2_grid = o2_grid, seed_ids = seed_ids, simulation_dir = simulation_dir))
}

# Analytical-simulation agreement helpers.
path_results_suffix <- function(path) {
  path <- normalizePath(as.character(path), mustWork = FALSE)
  suffix <- sub("^.*(/oxygen/results/.*)$", "\\1", path)
  ifelse(identical(suffix, path), path, suffix)
}

split_csv <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) {
    return(default)
  }
  vals <- trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

normalize_seed_ids <- function(x) {
  if (is.null(x) || !length(x)) return(character())
  o2ipa_norm_seed(x)
}

as_num_vec <- function(x, default) {
  vals <- suppressWarnings(as.numeric(split_csv(x, as.character(default))))
  vals <- vals[is.finite(vals)]
  if (length(vals)) unique(vals) else default
}

as_int_vec <- function(x, default) {
  vals <- suppressWarnings(as.integer(split_csv(x, as.character(default))))
  vals <- vals[!is.na(vals)]
  if (length(vals)) unique(vals) else default
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  tolower(as.character(x[[1]])) %in% c("true", "1", "yes", "y", "on")
}

num_key <- function(x) {
  vapply(as.numeric(x), function(val) {
    if (!is.finite(val)) return(NA_character_)
    format(signif(val, 12), scientific = FALSE, trim = TRUE)
  }, character(1))
}

num_path_tag <- function(x) {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) != 1L || !is.finite(val)) return("NA")
  txt <- format(val, scientific = FALSE, trim = TRUE)
  txt <- sub("^-", "m", txt)
  txt <- gsub("\\.", "p", txt)
  txt <- gsub("[^A-Za-z0-9]+", "", txt)
  if (!nzchar(txt)) "NA" else txt
}

read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing file: ", path)
    return(data.frame())
  }
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (is_gz) {
      gzip <- Sys.which("gzip")
      if (!nzchar(gzip)) stop("Reading gzip-compressed tables requires gzip on PATH: ", path)
      return(as.data.frame(data.table::fread(cmd = paste(shQuote(gzip), "-cd", shQuote(path)), sep = "\t", data.table = FALSE, showProgress = FALSE)))
    }
    return(as.data.frame(data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE)))
  }
  if (is_gz) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(utils::read.table(
      con,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ))
  }
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    wrote <- FALSE
    if (requireNamespace("data.table", quietly = TRUE)) {
      wrote <- tryCatch({
        data.table::fwrite(x, file = path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
        TRUE
      }, error = function(e) FALSE)
    }
    if (!wrote) {
      con <- gzfile(path, open = "wt")
      on.exit(close(con), add = TRUE)
      utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    }
  } else {
    utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  invisible(path)
}

format_o2_label <- function(x) {
  paste0("O2 = ", format(as.numeric(x), scientific = FALSE, trim = TRUE))
}

initial_condition_from_ploidy <- function(x) {
  paste0("init_", format(as.numeric(x), scientific = FALSE, trim = TRUE), "N")
}

filter_by_values <- function(df, col, values) {
  if (!nrow(df) || !(col %in% names(df))) return(df)
  key <- num_key(df[[col]])
  keep <- key %in% num_key(values)
  df[keep, , drop = FALSE]
}

method_slug <- function(x) {
  x <- tolower(trimws(as.character(x[[1]])))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "analytical" else x
}

method_label <- function(x) {
  x <- method_slug(x)
  if (identical(x, "eigen")) return("Eigen analytical")
  if (identical(x, "expm")) return("Expm analytical")
  paste0(x, " analytical")
}

analytical_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  fm <- fixo2_fixed_matrix(model_env, cfg, run_params, O2)
  list(M = fm$M, ngrid = fm$ngrid)
}

analytical_init_vector <- function(ngrid, init_N) {
  fixo2_init_vector(ngrid, init_N, n_unit = 22)
}

analytical_normalize_state <- function(x) {
  fixo2_normalize_state(x)
}

analytical_state_metrics <- function(state, ngrid, n_unit) {
  w <- analytical_normalize_state(state)
  if (all(is.na(w))) {
    return(data.frame(
      analytical_mean_N = NA_real_,
      analytical_mean_ploidy = NA_real_,
      analytical_sd_ploidy = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  ploidy_grid <- ngrid / n_unit
  mean_N <- sum(ngrid * w, na.rm = TRUE)
  mean_ploidy <- sum(ploidy_grid * w, na.rm = TRUE)
  second_ploidy <- sum(ploidy_grid^2 * w, na.rm = TRUE)
  data.frame(
    analytical_mean_N = mean_N,
    analytical_mean_ploidy = mean_ploidy,
    analytical_sd_ploidy = sqrt(max(0, second_ploidy - mean_ploidy^2)),
    stringsAsFactors = FALSE
  )
}

analytical_eigen_states <- function(M, init, time_grid) {
  fixo2_eigen_states(M, init, time_grid)
}

analytical_expm_states <- function(M, init, time_grid) {
  fixo2_expm_states(M, init, time_grid)
}

generate_seed_analytical_trajectories <- function(seed_id, inputs, param_mat, model_env,
                                                  time_points, o2_values, initial_ploidy_values,
                                                  analytical_methods) {
  manifest <- inputs$manifest[inputs$manifest$seed_id == seed_id, , drop = FALSE]
  if (!nrow(manifest) || !seed_id %in% rownames(param_mat)) return(data.frame())
  cfg <- o2pr_first_seed_cfg(manifest)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
  names(pvec) <- colnames(param_mat)
  run_params <- o2pr_run_params_from_vec(pvec, cfg)
  init_specs <- data.frame(
    initial_condition = paste0("init_", format(initial_ploidy_values, scientific = FALSE, trim = TRUE), "N"),
    initial_ploidy = initial_ploidy_values,
    requested_N = initial_ploidy_values * n_unit,
    stringsAsFactors = FALSE
  )

  rows <- list()
  k <- 0L
  for (O2 in o2_values) {
    fm <- analytical_fixed_matrix(model_env, cfg, run_params, O2)
    for (j in seq_len(nrow(init_specs))) {
      init <- analytical_init_vector(fm$ngrid, init_specs$requested_N[[j]])
      states_by_method <- list()
      if ("eigen" %in% analytical_methods) {
        states_by_method$eigen <- tryCatch(analytical_eigen_states(fm$M, init$vector, time_points), error = function(e) NULL)
      }
      if ("expm" %in% analytical_methods) {
        states_by_method$expm <- tryCatch(analytical_expm_states(fm$M, init$vector, time_points), error = function(e) NULL)
      }
      for (method in names(states_by_method)) {
        states <- states_by_method[[method]]
        if (is.null(states)) next
        for (i in seq_along(time_points)) {
          met <- analytical_state_metrics(states[[i]], fm$ngrid, n_unit)
          k <- k + 1L
          rows[[k]] <- data.frame(
            seed_id = seed_id,
            O2_pct = O2,
            O2_key = num_key(O2),
            initial_condition = init_specs$initial_condition[[j]],
            initial_ploidy = init_specs$initial_ploidy[[j]],
            requested_initial_N = init_specs$requested_N[[j]],
            used_initial_N = init$used_N,
            day = time_points[[i]],
            day_key = num_key(time_points[[i]]),
            analytical_method = method,
            analytical_method_label = method_label(method),
            met,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame()
  out
}

generate_analytical_trajectories <- function(run_dir, time_points, o2_values, initial_ploidy_values,
                                             analytical_methods, n_workers = 1L, seed_ids = NULL) {
  if (!dir.exists(run_dir)) stop("run_dir does not exist and is required to generate analytical trajectories: ", run_dir)
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  seeds <- if (is.null(seed_ids) || !length(seed_ids)) rownames(param_mat) else intersect(seed_ids, rownames(param_mat))
  if (!length(seeds)) stop("No seed parameters were found for analytical trajectory generation.")
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, length(seeds))
  message("Generating analytical trajectories from fitted parameters: ", length(seeds), " seeds, methods=", paste(analytical_methods, collapse = ","), ", workers=", n_workers)
  worker <- function(seed_id) {
    generate_seed_analytical_trajectories(
      seed_id = seed_id,
      inputs = inputs,
      param_mat = param_mat,
      model_env = model_env,
      time_points = time_points,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods
    )
  }
  rows <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    parallel::mclapply(seeds, worker, mc.cores = n_workers)
  } else {
    lapply(seeds, worker)
  }
  out <- do.call(rbind, rows[vapply(rows, nrow, integer(1)) > 0L])
  if (is.null(out)) out <- data.frame()
  if (nrow(out)) out$analytical_source_run_dir <- normalizePath(run_dir, mustWork = FALSE)
  out
}

filter_analytical_trajectories <- function(analytical, time_points, o2_values, initial_ploidy_values, analytical_methods, seed_ids = NULL) {
  if (!nrow(analytical)) return(analytical)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "analytical_mean_N", "analytical_mean_ploidy", "analytical_sd_ploidy"),
    names(analytical)
  )
  for (col in numeric_cols) analytical[[col]] <- suppressWarnings(as.numeric(analytical[[col]]))
  analytical$analytical_method <- vapply(analytical$analytical_method, method_slug, character(1))
  out <- filter_by_values(analytical, "day", time_points)
  out <- filter_by_values(out, "O2_pct", o2_values)
  out <- filter_by_values(out, "initial_ploidy", initial_ploidy_values)
  out <- out[out$analytical_method %in% analytical_methods, , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids)) out <- out[out$seed_id %in% seed_ids, , drop = FALSE]
  out$O2_key <- num_key(out$O2_pct)
  out$day_key <- num_key(out$day)
  out
}

expected_analytical_seed_ids <- function(run_dir, seed_ids = NULL) {
  seed_ids <- normalize_seed_ids(seed_ids)
  if (length(seed_ids)) return(seed_ids)
  if (is.null(run_dir) || !nzchar(run_dir) || !dir.exists(run_dir)) return(character())
  seeds <- tryCatch(o2ipa_discover_seeds(run_dir)$seed_id, error = function(e) character())
  normalize_seed_ids(seeds)
}

analytical_cache_missing_keys <- function(analytical, time_points, o2_values, initial_ploidy_values,
                                          analytical_methods, expected_seed_ids = character(), expected_run_dir = NULL) {
  if (!nrow(analytical)) return("all")
  required_cols <- c("seed_id", "O2_pct", "initial_ploidy", "day", "analytical_method", "analytical_mean_ploidy")
  missing_cols <- setdiff(required_cols, names(analytical))
  if (length(missing_cols)) return(paste0("missing column: ", paste(missing_cols, collapse = ",")))

  if (!is.null(expected_run_dir) && nzchar(expected_run_dir)) {
    expected_run_dir <- normalizePath(expected_run_dir, mustWork = FALSE)
    if (!"analytical_source_run_dir" %in% names(analytical)) return("missing analytical source run_dir")
    observed_run_dirs <- unique(normalizePath(as.character(analytical$analytical_source_run_dir), mustWork = FALSE))
    observed_ok <- observed_run_dirs %in% expected_run_dir |
      path_results_suffix(observed_run_dirs) %in% path_results_suffix(expected_run_dir)
    if (!length(observed_run_dirs) || any(is.na(observed_run_dirs)) || !all(observed_ok)) {
      return(paste0("analytical source run_dir mismatch: ", paste(head(observed_run_dirs, 3L), collapse = ",")))
    }
  }

  seeds <- normalize_seed_ids(expected_seed_ids)
  if (!length(seeds)) seeds <- normalize_seed_ids(unique(analytical$seed_id))
  if (!length(seeds)) return("no seeds")

  cache <- data.frame(
    seed_id = normalize_seed_ids(analytical$seed_id),
    O2_key = num_key(analytical$O2_pct),
    initial_ploidy_key = num_key(analytical$initial_ploidy),
    day_key = num_key(analytical$day),
    analytical_method = vapply(analytical$analytical_method, method_slug, character(1)),
    analytical_mean_ploidy = suppressWarnings(as.numeric(analytical$analytical_mean_ploidy)),
    stringsAsFactors = FALSE
  )
  available_keys <- unique(paste(cache$seed_id, cache$O2_key, cache$initial_ploidy_key, cache$day_key, cache$analytical_method, sep = "\r"))
  required <- expand.grid(
    seed_id = seeds,
    O2_key = num_key(o2_values),
    initial_ploidy_key = num_key(initial_ploidy_values),
    day_key = num_key(time_points),
    analytical_method = analytical_methods,
    stringsAsFactors = FALSE
  )
  required_keys <- paste(required$seed_id, required$O2_key, required$initial_ploidy_key, required$day_key, required$analytical_method, sep = "\r")
  missing <- setdiff(required_keys, available_keys)
  if (!length(missing)) return(character())
  head(missing, 5L)
}

generate_fixo2_attractor_mode_table <- function(run_dir, o2_values, seed_ids = NULL, n_workers = 1L) {
  if (is.null(run_dir) || !nzchar(run_dir) || !dir.exists(run_dir)) {
    stop("run_dir is required to generate FixO2 attractor mode table: ", run_dir)
  }
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  seeds <- if (is.null(seed_ids) || !length(seed_ids)) rownames(param_mat) else intersect(normalize_seed_ids(seed_ids), rownames(param_mat))
  if (!length(seeds)) stop("No seed parameters were found for FixO2 attractor mode generation.")
  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  o2_values <- sort(unique(as.numeric(o2_values)))
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, length(seeds))
  message("Generating FixO2 attractor mode table: ", length(seeds), " seeds, ", length(o2_values), " O2 values, workers=", n_workers)
  worker <- function(seed) {
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    rows <- lapply(o2_values, function(O2) fo2_dominant_attractor_one(seed, run_params, model_env, cfg, O2))
    do.call(rbind, rows)
  }
  rows <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    parallel::mclapply(seeds, worker, mc.cores = n_workers)
  } else {
    lapply(seeds, worker)
  }
  attractors <- do.call(rbind, rows[vapply(rows, nrow, integer(1)) > 0L])
  if (is.null(attractors)) attractors <- data.frame()
  if (!nrow(attractors)) return(data.frame())
  manifest <- inputs$manifest
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  attractors <- merge(
    attractors,
    manifest[, intersect(c("seed_id", "objective", "delta_objective"), names(manifest)), drop = FALSE],
    by = "seed_id",
    all.x = TRUE,
    sort = FALSE
  )
  attractors <- fixo2_assign_attractor_modes(attractors, "dominant_mean_ploidy")
  fixo2_attractor_mode_table(attractors)
}

read_seed_objectives <- function(analysis_dir, fit_dir = NULL, o2_values = NULL, seed_ids = NULL, n_workers = 1L, mode_reference_o2 = 2) {
  mode_reference_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_mode_reference_by_seed.tsv"
  )
  mode_seed_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_mode_by_seed.tsv"
  )
  mode_o2_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_mode_by_seed_o2.tsv"
  )
  mode_tab <- data.frame(seed_id = character(), trajectory_regime = character(), mode_label = character(), stringsAsFactors = FALSE)
  if (file.exists(mode_reference_path) || file.exists(mode_seed_path)) {
    mode_path <- if (file.exists(mode_reference_path)) mode_reference_path else mode_seed_path
    tab <- read_tsv(mode_path)
    cols <- intersect(c(
      "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
      "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
      "mode_reference_dominant_mean_ploidy", "mode_reference_status",
      "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
    ), names(tab))
    if ("seed_id" %in% cols) {
      mode_tab <- tab[, cols, drop = FALSE]
      mode_tab <- mode_tab[!duplicated(mode_tab$seed_id), , drop = FALSE]
      if ("mode_reference_o2_pct" %in% names(mode_tab)) {
        ref_vals <- suppressWarnings(as.numeric(mode_tab$mode_reference_o2_pct))
        if (any(is.finite(ref_vals)) && !any(abs(ref_vals - mode_reference_o2) < 1e-9, na.rm = TRUE)) {
          warning("Existing reference mode table was generated for a different mode_reference_o2; regenerating if fit_dir is available.")
          mode_tab <- data.frame(seed_id = character(), trajectory_regime = character(), mode_label = character(), stringsAsFactors = FALSE)
        }
      } else {
        warning("Existing seed-level mode table does not record mode_reference_o2; regenerating if fit_dir is available.")
        mode_tab <- data.frame(seed_id = character(), trajectory_regime = character(), mode_label = character(), stringsAsFactors = FALSE)
      }
    }
  } else if (file.exists(mode_o2_path)) {
    tab <- read_tsv(mode_o2_path)
    if ("O2_pct" %in% names(tab)) {
      reference_tab <- tryCatch(
        fixo2_reference_mode_table(tab, mode_reference_o2),
        error = function(e) {
          warning(conditionMessage(e))
          data.frame()
        }
      )
      if (nrow(reference_tab)) {
        cols <- intersect(c(
          "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
          "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
          "mode_reference_dominant_mean_ploidy", "mode_reference_status",
          "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
        ), names(reference_tab))
        mode_tab <- reference_tab[, cols, drop = FALSE]
        dir.create(dirname(mode_reference_path), recursive = TRUE, showWarnings = FALSE)
        write_tsv(reference_tab, mode_reference_path)
        write_tsv(reference_tab, mode_seed_path)
      }
    }
  }
  if (!nrow(mode_tab) && !is.null(fit_dir) && nzchar(fit_dir) && dir.exists(fit_dir) && !is.null(o2_values) && length(o2_values)) {
    mode_by_seed_o2 <- generate_fixo2_attractor_mode_table(
      run_dir = fit_dir,
      o2_values = sort(unique(c(o2_values, mode_reference_o2))),
      seed_ids = seed_ids,
      n_workers = n_workers
    )
    mode_reference_tab <- fixo2_reference_mode_table(mode_by_seed_o2, mode_reference_o2)
    mode_cols <- intersect(c(
      "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
      "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
      "mode_reference_dominant_mean_ploidy", "mode_reference_status",
      "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
    ), names(mode_reference_tab))
    mode_tab <- mode_reference_tab[, mode_cols, drop = FALSE]
    dir.create(dirname(mode_o2_path), recursive = TRUE, showWarnings = FALSE)
    write_tsv(mode_by_seed_o2, mode_o2_path)
    write_tsv(mode_reference_tab, mode_reference_path)
    write_tsv(mode_reference_tab, mode_seed_path)
  }

  if (is.null(fit_dir) || !nzchar(fit_dir) || !dir.exists(fit_dir)) {
    warning("Final seed objective values require fit_dir, but fit_dir is unavailable. Objective-colored plots will have missing objective values.")
    mode_tab$objective <- NA_real_
    mode_tab$objective_source <- NA_character_
    return(mode_tab)
  }

  read_summary_objectives <- function(fit_dir) {
    summary_path <- file.path(fit_dir, "extra_results", "seed_summary.tsv")
    if (!file.exists(summary_path)) return(data.frame())
    tab <- tryCatch(read_tsv(summary_path), error = function(e) data.frame())
    if (!nrow(tab)) return(data.frame())
    if (!"seed_id" %in% names(tab)) {
      if ("seed" %in% names(tab)) {
        tab$seed_id <- normalize_seed_ids(tab$seed)
      } else {
        return(data.frame())
      }
    }
    final_cols <- intersect(
      c("objective_total", "objective", "optimizer_local_objective", "optimizer_deoptim_objective", "objective_data"),
      names(tab)
    )
    if (!length(final_cols)) return(data.frame())
    objective <- rep(NA_real_, nrow(tab))
    source <- rep(NA_character_, nrow(tab))
    for (col in final_cols) {
      vals <- suppressWarnings(as.numeric(tab[[col]]))
      fill <- !is.finite(objective) & is.finite(vals)
      objective[fill] <- vals[fill]
      source[fill] <- col
    }
    out <- data.frame(
      seed_id = normalize_seed_ids(tab$seed_id),
      objective = objective,
      objective_source = source,
      stringsAsFactors = FALSE
    )
    out <- out[is.finite(out$objective), , drop = FALSE]
    out[!duplicated(out$seed_id), , drop = FALSE]
  }

  objectives <- read_summary_objectives(fit_dir)
  if (nrow(objectives)) {
    out <- if (nrow(mode_tab)) merge(mode_tab, objectives, by = "seed_id", all = TRUE) else objectives
    return(out)
  }

  seed_dirs <- list.dirs(fit_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("seed[0-9]+$", basename(seed_dirs))]
  rows <- lapply(seed_dirs, function(seed_dir) {
    candidates <- c(
      file.path(seed_dir, "fit_summary.tsv"),
      file.path(seed_dir, "best_params.tsv"),
      file.path(seed_dir, "best_parameters.tsv")
    )
    hits <- candidates[file.exists(candidates)]
    path <- if (length(hits)) hits[[1]] else NA_character_
    if (is.na(path)) return(NULL)
    tab <- tryCatch(read_tsv(path), error = function(e) data.frame())
    if (!nrow(tab)) return(NULL)
    if (all(c("metric", "value") %in% names(tab))) {
      vals <- as.list(tab$value)
      names(vals) <- tab$metric
      objective_cols <- c("objective_total", "objective", "optimizer_local_objective", "optimizer_deoptim_objective", "objective_data")
      hit <- objective_cols[vapply(objective_cols, function(col) {
        is.finite(suppressWarnings(as.numeric(vals[[col]])))
      }, logical(1))]
      if (!length(hit)) return(NULL)
      objective <- suppressWarnings(as.numeric(vals[[hit[[1]]]]))
      source <- hit[[1]]
    } else {
      objective_cols <- intersect(
        c("objective_total", "objective", "optimizer_local_objective", "optimizer_deoptim_objective", "objective_data"),
        names(tab)
      )
      if (!length(objective_cols)) return(NULL)
      objective <- NA_real_
      source <- NA_character_
      for (col in objective_cols) {
        val <- suppressWarnings(as.numeric(tab[[col]][[1]]))
        if (is.finite(val)) {
          objective <- val
          source <- col
          break
        }
      }
      if (!is.finite(objective)) return(NULL)
    }
    data.frame(
      seed_id = basename(seed_dir),
      objective = objective,
      objective_source = source,
      stringsAsFactors = FALSE
    )
  })
  objectives <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(objectives)) objectives <- data.frame(seed_id = character(), objective = numeric(), objective_source = character(), stringsAsFactors = FALSE)
  if (nrow(mode_tab)) merge(mode_tab, objectives, by = "seed_id", all = TRUE) else objectives
}

task_table_path <- function(simulation_dir, simulation_mode) {
  file.path(simulation_dir, simulation_mode, "task_list.tsv")
}

build_task_table_from_paths <- function(simulation_dir, simulation_mode, seed_ids, o2_values, initial_ploidy_values, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    for (O2 in o2_values) {
      for (initial_ploidy in initial_ploidy_values) {
        for (simulation_id in simulation_ids) {
          k <- k + 1L
          output_dir <- file.path(
            simulation_dir,
            simulation_mode,
            paste0("O2_", num_path_tag(O2)),
            paste0("ploidy", num_path_tag(initial_ploidy)),
            seed_id,
            paste0("sim", simulation_id)
          )
          rows[[k]] <- data.frame(
            task_id = k,
            seed_id = seed_id,
            fixed_o2_pct = O2,
            initial_ploidy = initial_ploidy,
            initial_condition = initial_condition_from_ploidy(initial_ploidy),
            simulation_id = simulation_id,
            output_dir = output_dir,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

read_simulation_tasks <- function(simulation_dir, simulation_mode, analytical, o2_values, initial_ploidy_values, simulation_ids) {
  path <- task_table_path(simulation_dir, simulation_mode)
  if (file.exists(path)) {
    tasks <- read_tsv(path)
    if (!"fixed_o2_pct" %in% names(tasks) && "O2_pct" %in% names(tasks)) {
      names(tasks)[names(tasks) == "O2_pct"] <- "fixed_o2_pct"
    }
    required <- c("seed_id", "fixed_o2_pct", "initial_ploidy", "simulation_id", "output_dir")
    missing <- setdiff(required, names(tasks))
    if (length(missing)) stop("Simulation task table is missing column(s): ", paste(missing, collapse = ", "))
    if (!"initial_condition" %in% names(tasks)) {
      tasks$initial_condition <- initial_condition_from_ploidy(tasks$initial_ploidy)
    }
  } else {
    seed_ids <- sort(unique(analytical$seed_id))
    tasks <- build_task_table_from_paths(
      simulation_dir = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      simulation_ids = simulation_ids
    )
  }
  tasks$fixed_o2_pct <- as.numeric(tasks$fixed_o2_pct)
  tasks$initial_ploidy <- as.numeric(tasks$initial_ploidy)
  tasks$simulation_id <- as.integer(tasks$simulation_id)
  expected_output_dir <- file.path(
    simulation_dir,
    simulation_mode,
    paste0("O2_", vapply(tasks$fixed_o2_pct, num_path_tag, character(1))),
    paste0("ploidy", vapply(tasks$initial_ploidy, num_path_tag, character(1))),
    tasks$seed_id,
    paste0("sim", tasks$simulation_id)
  )
  use_expected <- !dir.exists(tasks$output_dir) & dir.exists(expected_output_dir)
  tasks$output_dir[use_expected] <- expected_output_dir[use_expected]
  tasks <- filter_by_values(tasks, "fixed_o2_pct", o2_values)
  tasks <- filter_by_values(tasks, "initial_ploidy", initial_ploidy_values)
  tasks <- tasks[tasks$simulation_id %in% simulation_ids, , drop = FALSE]
  tasks <- tasks[tasks$seed_id %in% unique(analytical$seed_id), , drop = FALSE]
  tasks$state_file <- file.path(tasks$output_dir, "state_trajectory.tsv.gz")
  tasks
}

read_state_metrics_awk <- function(path, time_points = NULL) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    return(data.frame())
  }
  gzip <- Sys.which("gzip")
  awk <- Sys.which("awk")
  if (!nzchar(gzip) || !nzchar(awk)) {
    stop("Reading state trajectories requires gzip and awk on PATH.")
  }
  all_days <- is.null(time_points) || !length(time_points)
  days_arg <- if (all_days) "" else paste(num_key(time_points), collapse = ",")
  awk_script <- paste(
    'BEGIN {',
    '  keep_all = all_days + 0;',
    '  if (!keep_all) {',
    '    split(days, d, ",");',
    '    for (i in d) keep[sprintf("%.10g", d[i] + 0)] = 1;',
    '  }',
    '  OFS = "\t";',
    '}',
    'NR == 1 {',
    '  for (i = 1; i <= NF; i++) idx[$i] = i;',
    '  next;',
    '}',
    '{',
    '  day = sprintf("%.10g", $(idx["day"]) + 0);',
    '  if ((keep_all || (day in keep)) && $(idx["status"]) == "live") {',
    '    c = $(idx["cell_count"]) + 0;',
    '    n = $(idx["N"]) + 0;',
    '    p = $(idx["ploidy"]) + 0;',
    '    sumc[day] += c;',
    '    sumn[day] += n * c;',
    '    sump[day] += p * c;',
    '    sump2[day] += p * p * c;',
    '  }',
    '}',
    'END {',
    '  for (day in sumc) {',
    '    if (sumc[day] > 0) {',
    '      meanp = sump[day] / sumc[day];',
    '      varp = sump2[day] / sumc[day] - meanp * meanp;',
    '      if (varp < 0 && varp > -1e-9) varp = 0;',
    '      print day, sumn[day] / sumc[day], meanp, sqrt(varp), sumc[day];',
    '    }',
    '  }',
    '}',
    sep = "\n"
  )
  cmd <- paste(
    shQuote(gzip), "-cd", shQuote(path), "|",
    shQuote(awk),
    "-v", shQuote(paste0("days=", days_arg)),
    "-v", shQuote(paste0("all_days=", if (all_days) "1" else "0")),
    shQuote(awk_script)
  )
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
  if (!length(out)) return(data.frame())
  con <- textConnection(out)
  on.exit(close(con), add = TRUE)
  tab <- utils::read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(tab) <- c("day", "simulation_mean_N", "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells")
  tab$day <- as.numeric(tab$day)
  tab <- tab[order(tab$day), , drop = FALSE]
  tab
}

read_simulation_metrics_chunk <- function(tasks, idx, time_points = NULL, progress_every = 100L, worker_label = NULL, total_tasks = nrow(tasks)) {
  rows <- vector("list", length(idx))
  missing <- 0L
  for (j in seq_along(idx)) {
    i <- idx[[j]]
    if (progress_every > 0L && (j == 1L || j %% progress_every == 0L || j == length(idx))) {
      label <- if (is.null(worker_label) || !nzchar(worker_label)) "" else paste0(worker_label, ": ")
      message(label, "Reading simulation state metrics: ", j, "/", length(idx), " (task ", i, "/", total_tasks, ")")
    }
    task <- tasks[i, , drop = FALSE]
    metric <- read_state_metrics_awk(task$state_file[[1]], time_points)
    if (!nrow(metric)) {
      missing <- missing + 1L
      next
    }
    metric$seed_id <- task$seed_id[[1]]
    metric$O2_pct <- task$fixed_o2_pct[[1]]
    metric$O2_key <- num_key(metric$O2_pct)
    metric$initial_ploidy <- task$initial_ploidy[[1]]
    metric$initial_condition <- task$initial_condition[[1]]
    metric$simulation_id <- task$simulation_id[[1]]
    metric$day_key <- num_key(metric$day)
    rows[[j]] <- metric
  }
  if (missing) {
    label <- if (is.null(worker_label) || !nzchar(worker_label)) "" else paste0(worker_label, ": ")
    warning(label, "Missing or unreadable simulation state files: ", missing)
  }
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) out <- data.frame()
  out
}

read_simulation_metrics <- function(tasks, time_points = NULL, progress_every = 100L, n_workers = 1L) {
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, nrow(tasks))
  if (n_workers <= 1L || !identical(.Platform$OS.type, "unix")) {
    return(read_simulation_metrics_chunk(
      tasks = tasks,
      idx = seq_len(nrow(tasks)),
      time_points = time_points,
      progress_every = progress_every,
      total_tasks = nrow(tasks)
    ))
  }

  message("Reading simulation state metrics with ", n_workers, " workers: ", nrow(tasks), " tasks")
  idx <- seq_len(nrow(tasks))
  chunks <- split(idx, ((idx - 1L) %% n_workers) + 1L)
  res <- parallel::mclapply(
    seq_along(chunks),
    function(worker_id) {
      read_simulation_metrics_chunk(
        tasks = tasks,
        idx = chunks[[worker_id]],
        time_points = time_points,
        progress_every = progress_every,
        worker_label = sprintf("worker %02d", worker_id),
        total_tasks = nrow(tasks)
      )
    },
    mc.cores = n_workers
  )
  out <- do.call(rbind, res[vapply(res, nrow, integer(1)) > 0L])
  if (is.null(out)) out <- data.frame()
  out
}

filter_simulation_metrics <- function(sim_metrics, time_points, o2_values, initial_ploidy_values, simulation_ids, seed_ids = NULL) {
  if (!nrow(sim_metrics)) return(sim_metrics)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "simulation_id", "simulation_mean_N",
      "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells"),
    names(sim_metrics)
  )
  for (col in numeric_cols) sim_metrics[[col]] <- suppressWarnings(as.numeric(sim_metrics[[col]]))
  out <- filter_by_values(sim_metrics, "day", time_points)
  out <- filter_by_values(out, "O2_pct", o2_values)
  out <- filter_by_values(out, "initial_ploidy", initial_ploidy_values)
  if ("simulation_id" %in% names(out)) out <- out[out$simulation_id %in% simulation_ids, , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids) && "seed_id" %in% names(out)) {
    out <- out[out$seed_id %in% seed_ids, , drop = FALSE]
  }
  out$O2_key <- num_key(out$O2_pct)
  out$day_key <- num_key(out$day)
  out
}

aggregate_replicates <- function(sim_metrics) {
  if (!nrow(sim_metrics)) return(sim_metrics)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "simulation_id", "simulation_mean_N",
      "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells"),
    names(sim_metrics)
  )
  for (col in numeric_cols) sim_metrics[[col]] <- suppressWarnings(as.numeric(sim_metrics[[col]]))
  keys <- c("seed_id", "O2_pct", "O2_key", "initial_condition", "initial_ploidy", "day", "day_key")
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(sim_metrics)
    out <- dt[, .(
      simulation_n = data.table::uniqueN(simulation_id),
      simulation_mean_N = mean(simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(simulation_live_cells, na.rm = TRUE)
    ), by = keys]
    return(as.data.frame(out))
  }
  split_key <- interaction(sim_metrics[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sim_metrics, split_key), function(x) {
    data.frame(
      x[1, keys, drop = FALSE],
      simulation_n = length(unique(x$simulation_id)),
      simulation_mean_N = mean(x$simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(x$simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(x$simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(x$simulation_live_cells, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

merge_scatter_data <- function(analytical, sim_summary, objectives) {
  by_cols <- c("seed_id", "O2_key", "initial_condition", "day_key")
  dat <- merge(
    analytical,
    sim_summary,
    by = by_cols,
    all = FALSE,
    suffixes = c("_analytical", "_simulation")
  )
  if ("O2_pct_analytical" %in% names(dat)) dat$O2_pct <- dat$O2_pct_analytical
  if ("day_analytical" %in% names(dat)) dat$day <- dat$day_analytical
  objective_by <- "seed_id"
  if ("O2_key" %in% names(dat) && "O2_key" %in% names(objectives)) objective_by <- c("seed_id", "O2_key")
  dat <- merge(dat, objectives, by = objective_by, all.x = TRUE, suffixes = c("", "_objective"))
  if ("mode_label_objective" %in% names(dat) && "mode_label" %in% names(dat)) {
    fill <- !nzchar(as.character(dat$mode_label)) | is.na(dat$mode_label)
    dat$mode_label[fill] <- dat$mode_label_objective[fill]
  }
  dat$O2_factor <- factor(format_o2_label(dat$O2_pct), levels = format_o2_label(sort(unique(dat$O2_pct))))
  dat$day_factor <- factor(paste0("Day ", format(as.numeric(dat$day), scientific = FALSE, trim = TRUE)),
                           levels = paste0("Day ", format(sort(unique(as.numeric(dat$day))), scientific = FALSE, trim = TRUE)))
  dat$initial_condition <- factor(dat$initial_condition, levels = sort(unique(dat$initial_condition)))
  mode <- mode_values(dat)
  dat$mode_factor <- factor(mode, levels = mode_levels(mode))
  dat$objective <- as.numeric(dat$objective)
  dat[is.finite(dat$analytical_mean_ploidy) & is.finite(dat$simulation_mean_ploidy), , drop = FALSE]
}

plot_limits <- function(dat) {
  vals <- c(dat$analytical_mean_ploidy, dat$simulation_mean_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(0, 1))
  rng <- range(vals)
  pad <- diff(rng) * 0.04
  if (!is.finite(pad) || pad <= 0) pad <- 0.05
  rng + c(-pad, pad)
}

objective_aesthetic <- function(dat, transform = "identity") {
  objective <- as.numeric(dat$objective)
  label <- "Final objective"
  if (identical(transform, "log10")) {
    objective <- ifelse(is.finite(objective) & objective > 0, log10(objective), NA_real_)
    label <- "log10(final objective)"
  }
  list(value = objective, label = label)
}

mode_values <- function(dat) {
  mode <- rep("unknown", nrow(dat))
  if ("mode_label" %in% names(dat)) {
    mode <- as.character(dat$mode_label)
  } else if ("trajectory_regime" %in% names(dat)) {
    mode <- as.character(dat$trajectory_regime)
  }
  mode[is.na(mode) | !nzchar(mode)] <- "unknown"
  mode
}

mode_levels <- function(mode) {
  preferred <- c("mode1", "mode2", "unknown")
  c(intersect(preferred, unique(mode)), sort(setdiff(unique(mode), preferred)))
}

mode_palette <- function(levels) {
  base <- c(
    mode1 = "#0072B2",
    mode2 = "#D55E00",
    unknown = "#C9C9C9"
  )
  missing <- setdiff(levels, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), "Dark 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[levels]
}

base_scatter <- function(dat, limits, point_size = 0.9, alpha = 0.55) {
  ggplot2::ggplot(dat, ggplot2::aes(x = analytical_mean_ploidy, y = simulation_mean_ploidy)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey55", alpha = 0.28, linewidth = 0.25) +
    ggplot2::coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    ggplot2::labs(
      x = "Analytical solution mean ploidy",
      y = "Simulation-inferred mean ploidy"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
}

discrete_guides <- function() {
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(shape = 21, color = "grey35", stroke = 0.25, alpha = 1, size = 2.6)),
    shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
  )
}

save_plot <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", device = "pdf")
  invisible(path)
}

plot_time_facets_color_o2 <- function(dat, path, limits) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = O2_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::labs(fill = "Fixed O2", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

plot_time_facets_color_mode <- function(dat, path, limits) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

plot_time_facets_color_objective <- function(dat, path, limits, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

plot_o2_facets_color_mode <- function(dat, path, limits, title = NULL) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition", title = title) +
    discrete_guides()
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}

plot_o2_facets_color_objective <- function(dat, path, limits, title = NULL, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition", title = title) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}

augment_comparison_data <- function(dat) {
  dat$residual_ploidy <- as.numeric(dat$simulation_mean_ploidy) - as.numeric(dat$analytical_mean_ploidy)
  dat$agreement_mean_ploidy <- (as.numeric(dat$simulation_mean_ploidy) + as.numeric(dat$analytical_mean_ploidy)) / 2
  dat[is.finite(dat$residual_ploidy) & is.finite(dat$agreement_mean_ploidy), , drop = FALSE]
}

residual_limits <- function(dat) {
  vals <- as.numeric(dat$residual_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(-1, 1))
  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 0.05
  pad <- max_abs * 0.06
  c(-(max_abs + pad), max_abs + pad)
}

bland_altman_stats <- function(dat) {
  vals <- as.numeric(dat$residual_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(bias = 0, lower = NA_real_, upper = NA_real_))
  bias <- mean(vals)
  sd_val <- stats::sd(vals)
  if (!is.finite(sd_val)) sd_val <- 0
  c(bias = bias, lower = bias - 1.96 * sd_val, upper = bias + 1.96 * sd_val)
}

base_comparison_plot <- function(dat, comparison = c("residual", "bland_altman"), title = NULL) {
  comparison <- match.arg(comparison)
  dat <- augment_comparison_data(dat)
  if (identical(comparison, "residual")) {
    dat$comparison_x <- as.numeric(dat$analytical_mean_ploidy)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = comparison_x, y = residual_ploidy)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey50", alpha = 0.35, linewidth = 0.3) +
      ggplot2::labs(
        x = "Analytical solution mean ploidy",
        y = "Simulation - analytical mean ploidy",
        title = title
      )
  } else {
    dat$comparison_x <- as.numeric(dat$agreement_mean_ploidy)
    ba <- bland_altman_stats(dat)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = comparison_x, y = residual_ploidy)) +
      ggplot2::geom_hline(yintercept = ba[["bias"]], color = "grey35", alpha = 0.45, linewidth = 0.35) +
      ggplot2::geom_hline(yintercept = ba[c("lower", "upper")], color = "grey45", alpha = 0.35, linetype = 2, linewidth = 0.3) +
      ggplot2::labs(
        x = "Mean of analytical and simulation mean ploidy",
        y = "Simulation - analytical mean ploidy",
        title = title
      )
  }
  p +
    ggplot2::coord_cartesian(ylim = residual_limits(dat), expand = FALSE) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
}

plot_time_facets_color_o2_comparison <- function(dat, path, comparison) {
  p <- base_comparison_plot(dat, comparison = comparison) +
    ggplot2::geom_point(
      ggplot2::aes(fill = O2_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::labs(fill = "Fixed O2", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

plot_time_facets_color_objective_comparison <- function(dat, path, comparison, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_comparison_plot(dat, comparison = comparison) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

plot_time_facets_color_mode_comparison <- function(dat, path, comparison) {
  p <- base_comparison_plot(dat, comparison = comparison) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

plot_o2_facets_color_objective_comparison <- function(dat, path, comparison, title = NULL, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_comparison_plot(dat, comparison = comparison, title = title) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}

plot_o2_facets_color_mode_comparison <- function(dat, path, comparison, title = NULL) {
  p <- base_comparison_plot(dat, comparison = comparison, title = title) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}

make_comparison_outputs <- function(dat, fig_dir, comparison, prefix, objective_transform = "identity") {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  plot_time_facets_color_o2_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_time_color_o2.pdf")),
    comparison = comparison
  )
  plot_time_facets_color_objective_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_time_color_objective.pdf")),
    comparison = comparison,
    objective_transform = objective_transform
  )
  plot_time_facets_color_mode_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_time_color_mode.pdf")),
    comparison = comparison
  )
  plot_o2_facets_color_objective_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_o2_color_objective_all_times.pdf")),
    comparison = comparison,
    title = "All selected time points",
    objective_transform = objective_transform
  )
  plot_o2_facets_color_mode_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_o2_color_mode_all_times.pdf")),
    comparison = comparison,
    title = "All selected time points"
  )
  for (day in sort(unique(as.numeric(dat$day)))) {
    day_dat <- dat[abs(as.numeric(dat$day) - day) < 1e-9, , drop = FALSE]
    day_label <- format(day, scientific = FALSE, trim = TRUE)
    plot_o2_facets_color_objective_comparison(
      day_dat,
      file.path(fig_dir, paste0(prefix, "_by_o2_color_objective_day", day_label, ".pdf")),
      comparison = comparison,
      title = paste0("Day ", day_label),
      objective_transform = objective_transform
    )
    plot_o2_facets_color_mode_comparison(
      day_dat,
      file.path(fig_dir, paste0(prefix, "_by_o2_color_mode_day", day_label, ".pdf")),
      comparison = comparison,
      title = paste0("Day ", day_label)
    )
  }
  invisible(fig_dir)
}

plot_analytical_solution_vs_fixed_o2 <- function(analytical, path, analytical_method = "analytical") {
  dat <- analytical
  if (!nrow(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))
  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  dat$O2_factor <- factor(format_o2_label(dat$O2_pct), levels = format_o2_label(o2_levels))
  dat$day_factor <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  )
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$seed_initial_group <- interaction(dat$seed_id, dat$initial_condition, drop = TRUE)
  dat <- order_initial_condition_for_drawing(dat)

  summary <- stats::aggregate(
    analytical_mean_ploidy ~ O2_factor + day_factor + initial_condition,
    data = dat,
    FUN = stats::median,
    na.rm = TRUE
  )
  names(summary)[names(summary) == "analytical_mean_ploidy"] <- "median_analytical_mean_ploidy"
  summary <- order_initial_condition_for_drawing(summary)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = O2_factor, y = analytical_mean_ploidy)) +
    ggplot2::geom_line(
      ggplot2::aes(group = seed_initial_group, color = initial_condition),
      alpha = 0.08,
      linewidth = 0.22
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = initial_condition),
      alpha = 0.22,
      size = 0.35,
      stroke = 0
    ) +
    ggplot2::geom_line(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, group = initial_condition, color = initial_condition),
      linewidth = 0.75,
      alpha = 0.95
    ) +
    ggplot2::geom_point(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, color = initial_condition),
      size = 1.15,
      alpha = 0.95
    ) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::labs(
      x = "Fixed O2",
      y = "Analytical solution mean ploidy",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution across seeds")
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

analytical_initial_condition_levels <- function(x) {
  vals <- unique(as.character(x))
  preferred <- c("init_2N", "init_4N")
  c(preferred[preferred %in% vals], sort(setdiff(vals, preferred)))
}

analytical_initial_condition_palette <- function(levels) {
  base <- c(
    init_2N = "#1B9E77",
    init_4N = "#7B3294"
  )
  missing <- setdiff(levels, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[levels]
}

order_initial_condition_for_drawing <- function(dat) {
  if (!"initial_condition" %in% names(dat)) return(dat)
  draw_order <- match(as.character(dat$initial_condition), c("init_4N", "init_2N"))
  draw_order[is.na(draw_order)] <- 0L
  dat[order(draw_order), , drop = FALSE]
}

add_analytical_mode_labels <- function(analytical, seed_metadata = NULL) {
  dat <- analytical
  if (!nrow(dat) || !"seed_id" %in% names(dat)) return(dat)
  if (!"trajectory_regime" %in% names(dat) && is.data.frame(seed_metadata) && nrow(seed_metadata) && "seed_id" %in% names(seed_metadata)) {
    cols <- intersect(c(
      "seed_id", "O2_key", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
      "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
      "mode_reference_dominant_mean_ploidy", "mode_reference_status",
      "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
    ), names(seed_metadata))
    if (all(c("seed_id", "trajectory_regime") %in% cols)) {
      meta <- seed_metadata[, cols, drop = FALSE]
      if ("O2_key" %in% names(meta) && "O2_key" %in% names(dat)) {
        meta <- meta[!duplicated(meta[, c("seed_id", "O2_key"), drop = FALSE]), , drop = FALSE]
        dat <- merge(dat, meta, by = c("seed_id", "O2_key"), all.x = TRUE, sort = FALSE)
      } else {
        meta <- meta[!duplicated(meta$seed_id), , drop = FALSE]
        dat <- merge(dat, meta, by = "seed_id", all.x = TRUE, sort = FALSE)
      }
    }
  }
  dat
}

plot_analytical_solution_vs_fixed_o2_by_mode <- function(analytical, path, analytical_method = "analytical", seed_metadata = NULL) {
  dat <- add_analytical_mode_labels(analytical, seed_metadata = seed_metadata)
  if (!nrow(dat) || !"trajectory_regime" %in% names(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  dat <- dat[dat$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))

  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  day_levels_label <- paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  mode_labels <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "mode1",
    mode2_attractor_dominant_ploidy_lt_2 = "mode2"
  )
  dat$O2_factor <- factor(format_o2_label(dat$O2_pct), levels = format_o2_label(o2_levels))
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$mode_group <- factor(mode_labels[as.character(dat$trajectory_regime)], levels = c("mode1", "mode2"))
  dat$day_label <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = day_levels_label
  )
  dat$facet_label <- factor(
    paste(as.character(dat$day_label), as.character(dat$mode_group), sep = "\n"),
    levels = as.vector(rbind(
      paste(day_levels_label, "mode1", sep = "\n"),
      paste(day_levels_label, "mode2", sep = "\n")
    ))
  )
  dat$seed_initial_group <- interaction(dat$seed_id, dat$initial_condition, dat$mode_group, drop = TRUE)
  dat <- order_initial_condition_for_drawing(dat)

  summary <- stats::aggregate(
    analytical_mean_ploidy ~ O2_factor + facet_label + initial_condition,
    data = dat,
    FUN = stats::median,
    na.rm = TRUE
  )
  names(summary)[names(summary) == "analytical_mean_ploidy"] <- "median_analytical_mean_ploidy"
  summary <- order_initial_condition_for_drawing(summary)
  panel_background <- unique(dat[, c("facet_label", "mode_group"), drop = FALSE])

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = O2_factor, y = analytical_mean_ploidy)) +
    ggplot2::geom_rect(
      data = panel_background,
      ggplot2::aes(fill = mode_group),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE,
      alpha = 0.35,
      color = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = seed_initial_group, color = initial_condition),
      alpha = 0.08,
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = initial_condition),
      alpha = 0.22,
      size = 0.32,
      stroke = 0
    ) +
    ggplot2::geom_line(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, group = initial_condition, color = initial_condition),
      linewidth = 0.7,
      alpha = 0.95
    ) +
    ggplot2::geom_point(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, color = initial_condition),
      size = 1.05,
      alpha = 0.95
    ) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::scale_fill_manual(
      values = c(mode1 = "#DCEEFF", mode2 = "#FFE8CC"),
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::facet_wrap(~facet_label, nrow = 2, ncol = 8) +
    ggplot2::labs(
      x = "Fixed O2",
      y = "Analytical solution mean ploidy",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution by mode")
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 16, height = 7)
}

half_violin_polygon_data <- function(dat, x_col, y_col, side_col, group_cols, width = 0.36, n = 128L) {
  if (!nrow(dat)) return(data.frame())
  split_cols <- unique(c(group_cols, side_col))
  split_key <- interaction(dat[, split_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(dat, split_key), function(d) {
    y <- suppressWarnings(as.numeric(d[[y_col]]))
    y <- y[is.finite(y)]
    if (length(y) < 2L || length(unique(y)) < 2L) return(NULL)
    dens <- stats::density(y, n = n, from = min(y), to = max(y), na.rm = TRUE)
    if (!length(dens$y) || !is.finite(max(dens$y)) || max(dens$y) <= 0) return(NULL)
    side <- as.character(d[[side_col]][[1]])
    direction <- if (identical(side, "init_2N")) 1 else -1
    center <- suppressWarnings(as.numeric(d[[x_col]][[1]]))
    if (!is.finite(center)) return(NULL)
    outer_x <- center + direction * width * dens$y / max(dens$y)
    out <- data.frame(
      x = c(rep(center, length(dens$x)), rev(outer_x)),
      y = c(dens$x, rev(dens$x)),
      stringsAsFactors = FALSE
    )
    for (col in split_cols) out[[col]] <- d[[col]][[1]]
    out$violin_group <- paste(vapply(split_cols, function(col) as.character(d[[col]][[1]]), character(1)), collapse = "\r")
    out
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}

half_violin_median_data <- function(dat, x_col, y_col, side_col, group_cols, offset = 0.16) {
  if (!nrow(dat)) return(data.frame())
  split_cols <- unique(c(group_cols, side_col))
  rows <- lapply(split(dat, interaction(dat[, split_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)), function(d) {
    y <- suppressWarnings(as.numeric(d[[y_col]]))
    y <- y[is.finite(y)]
    if (!length(y)) return(NULL)
    side <- as.character(d[[side_col]][[1]])
    direction <- if (identical(side, "init_2N")) 1 else -1
    center <- suppressWarnings(as.numeric(d[[x_col]][[1]]))
    if (!is.finite(center)) return(NULL)
    out <- data.frame(
      x = center + direction * offset,
      x0 = center + direction * offset - 0.06,
      x1 = center + direction * offset + 0.06,
      median_y = stats::median(y, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (col in split_cols) out[[col]] <- d[[col]][[1]]
    out
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}

plot_analytical_solution_vs_fixed_o2_init_half_violin <- function(analytical, path, analytical_method = "analytical") {
  dat <- analytical
  if (!nrow(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))

  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)
  dat$O2_factor <- factor(format_o2_label(dat$O2_pct), levels = format_o2_label(o2_levels))
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$day_factor <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  )
  x_levels <- data.frame(
    O2_factor = levels(dat$O2_factor),
    x_pos = seq_along(levels(dat$O2_factor)),
    stringsAsFactors = FALSE
  )
  dat <- merge(dat, x_levels, by = "O2_factor", all.x = TRUE, sort = FALSE)
  dat <- order_initial_condition_for_drawing(dat)

  violin <- half_violin_polygon_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor")
  )
  if (!nrow(violin)) return(invisible(NULL))
  violin$initial_condition <- factor(violin$initial_condition, levels = init_levels)
  med <- half_violin_median_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor")
  )
  med$initial_condition <- factor(med$initial_condition, levels = init_levels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = violin,
      ggplot2::aes(x = x, y = y, group = violin_group, fill = initial_condition),
      color = NA,
      alpha = 0.62
    ) +
    ggplot2::geom_segment(
      data = med,
      ggplot2::aes(x = x0, xend = x1, y = median_y, yend = median_y, color = initial_condition),
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      data = med,
      ggplot2::aes(x = x, y = median_y, color = initial_condition),
      fill = "white",
      shape = 21,
      size = 1.15,
      stroke = 0.45
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_levels$x_pos,
      labels = x_levels$O2_factor,
      expand = ggplot2::expansion(mult = c(0.03, 0.03))
    ) +
    ggplot2::scale_fill_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::labs(
      x = "Fixed O2",
      y = "Analytical solution mean ploidy",
      fill = "Initial condition",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution distribution across fixed O2")
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}

plot_analytical_solution_vs_fixed_o2_mode_init_half_violin <- function(analytical, path, analytical_method = "analytical", seed_metadata = NULL) {
  dat <- add_analytical_mode_labels(analytical, seed_metadata = seed_metadata)
  if (!nrow(dat) || !"trajectory_regime" %in% names(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  dat <- dat[dat$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))

  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  day_levels_label <- paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  mode_labels <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "mode1",
    mode2_attractor_dominant_ploidy_lt_2 = "mode2"
  )
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)

  dat$O2_factor <- factor(format_o2_label(dat$O2_pct), levels = format_o2_label(o2_levels))
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$mode_group <- factor(mode_labels[as.character(dat$trajectory_regime)], levels = c("mode1", "mode2"))
  dat$day_factor <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = day_levels_label
  )
  x_levels <- expand.grid(
    mode_group = levels(dat$mode_group),
    O2_factor = levels(dat$O2_factor),
    stringsAsFactors = FALSE
  )
  x_levels <- x_levels[, c("O2_factor", "mode_group"), drop = FALSE]
  x_levels$x_pos <- seq_len(nrow(x_levels))
  x_levels$x_label <- paste(x_levels$O2_factor, x_levels$mode_group, sep = "\n")
  dat <- merge(dat, x_levels[, c("O2_factor", "mode_group", "x_pos"), drop = FALSE],
               by = c("O2_factor", "mode_group"), all.x = TRUE, sort = FALSE)
  dat <- order_initial_condition_for_drawing(dat)

  violin <- half_violin_polygon_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor", "mode_group")
  )
  if (!nrow(violin)) return(invisible(NULL))
  violin$initial_condition <- factor(violin$initial_condition, levels = init_levels)
  med <- half_violin_median_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor", "mode_group")
  )
  med$initial_condition <- factor(med$initial_condition, levels = init_levels)

  bg <- unique(dat[, c("day_factor", "x_pos", "mode_group"), drop = FALSE])
  bg$xmin <- bg$x_pos - 0.5
  bg$xmax <- bg$x_pos + 0.5
  bg_mode1 <- bg[bg$mode_group == "mode1", , drop = FALSE]
  bg_mode2 <- bg[bg$mode_group == "mode2", , drop = FALSE]

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = bg_mode1,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#DCEEFF",
      alpha = 0.42,
      color = NA
    ) +
    ggplot2::geom_rect(
      data = bg_mode2,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#FFE8CC",
      alpha = 0.42,
      color = NA
    ) +
    ggplot2::geom_polygon(
      data = violin,
      ggplot2::aes(x = x, y = y, group = violin_group, fill = initial_condition),
      color = NA,
      alpha = 0.58
    ) +
    ggplot2::geom_segment(
      data = med,
      ggplot2::aes(x = x0, xend = x1, y = median_y, yend = median_y, color = initial_condition),
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      data = med,
      ggplot2::aes(x = x, y = median_y, color = initial_condition),
      fill = "white",
      shape = 21,
      size = 1.05,
      stroke = 0.45
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_levels$x_pos,
      labels = x_levels$x_label,
      expand = ggplot2::expansion(mult = c(0.015, 0.015))
    ) +
    ggplot2::scale_fill_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::labs(
      x = "Fixed O2 and mode",
      y = "Analytical solution mean ploidy",
      fill = "Initial condition",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution distribution by fixed O2, mode, and initial condition")
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6.5),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 16, height = 8.5)
}

make_scatter_outputs <- function(dat, out_dir, objective_transform = "identity", analytical_method = "analytical") {
  method_dir <- file.path(out_dir, method_slug(analytical_method))
  scatter_dir <- file.path(method_dir, "scatter")
  residual_dir <- file.path(method_dir, "residual")
  bland_altman_dir <- file.path(method_dir, "bland_altman")
  dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
  limits <- plot_limits(dat)

  plot_time_facets_color_o2(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_time_color_o2.pdf"),
    limits = limits
  )
  plot_time_facets_color_objective(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_time_color_objective.pdf"),
    limits = limits,
    objective_transform = objective_transform
  )
  plot_time_facets_color_mode(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_time_color_mode.pdf"),
    limits = limits
  )
  plot_o2_facets_color_objective(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_o2_color_objective_all_times.pdf"),
    limits = limits,
    title = "All selected time points",
    objective_transform = objective_transform
  )
  plot_o2_facets_color_mode(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_o2_color_mode_all_times.pdf"),
    limits = limits,
    title = "All selected time points"
  )

  for (day in sort(unique(as.numeric(dat$day)))) {
    day_dat <- dat[abs(as.numeric(dat$day) - day) < 1e-9, , drop = FALSE]
    day_label <- format(day, scientific = FALSE, trim = TRUE)
    plot_o2_facets_color_objective(
      day_dat,
      file.path(scatter_dir, paste0("scatter_analytical_vs_simulation_by_o2_color_objective_day", day_label, ".pdf")),
      limits = limits,
      title = paste0("Day ", day_label),
      objective_transform = objective_transform
    )
    plot_o2_facets_color_mode(
      day_dat,
      file.path(scatter_dir, paste0("scatter_analytical_vs_simulation_by_o2_color_mode_day", day_label, ".pdf")),
      limits = limits,
      title = paste0("Day ", day_label)
    )
  }
  make_comparison_outputs(
    dat = dat,
    fig_dir = residual_dir,
    comparison = "residual",
    prefix = "residual_simulation_minus_analytical",
    objective_transform = objective_transform
  )
  make_comparison_outputs(
    dat = dat,
    fig_dir = bland_altman_dir,
    comparison = "bland_altman",
    prefix = "bland_altman_simulation_vs_analytical",
    objective_transform = objective_transform
  )
  invisible(method_dir)
}

make_analytical_solution_outputs <- function(analytical, out_dir, analytical_method = "analytical", seed_metadata = NULL) {
  method_dir <- file.path(out_dir, method_slug(analytical_method))
  plot_dir <- file.path(method_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_analytical_solution_vs_fixed_o2(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time.pdf"),
    analytical_method = analytical_method
  )
  plot_analytical_solution_vs_fixed_o2_init_half_violin(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time_init_half_violin_median.pdf"),
    analytical_method = analytical_method
  )
  plot_analytical_solution_vs_fixed_o2_by_mode(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time_mode1_mode2.pdf"),
    analytical_method = analytical_method,
    seed_metadata = seed_metadata
  )
  plot_analytical_solution_vs_fixed_o2_mode_init_half_violin(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time_mode_init_half_violin_median.pdf"),
    analytical_method = analytical_method,
    seed_metadata = seed_metadata
  )
  invisible(plot_dir)
}

fixo2_arg_value <- function(args, primary, fallback = NULL, default = NULL) {
  for (key in c(primary, fallback)) {
    if (!is.null(key) && !is.null(args[[key]]) && length(args[[key]]) && !is.na(args[[key]][[1]])) {
      return(args[[key]])
    }
  }
  default
}

fixo2_agreement_cache_read_path <- function(path, legacy_path, recompute = FALSE) {
  if (isTRUE(recompute) || file.exists(path)) return(path)
  legacy_path <- legacy_path[nzchar(legacy_path)]
  hit <- legacy_path[file.exists(legacy_path)]
  if (length(hit)) hit[[1]] else path
}

fixo2_run_analytical_simulation_agreement <- function(args, repo_root, run_dir, analysis_dir, out_dir, simulation_dir, o2_grid) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package is required for analytical-simulation agreement plots.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(out_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  log_dir <- file.path(out_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(log_dir, "FixO2_invivo_analytical_simulation_agreement.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  fit_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_fit_dir", "fit_dir", run_dir), run_dir),
    repo_root,
    mustWork = FALSE
  )
  agreement_run_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_run_dir", "run_dir", fit_dir), fit_dir),
    repo_root,
    mustWork = FALSE
  )
  agreement_analysis_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_analysis_dir", "analysis_dir", analysis_dir), analysis_dir),
    repo_root,
    mustWork = FALSE
  )
  agreement_simulation_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_simulation_dir", "simulation_dir", simulation_dir), simulation_dir),
    repo_root,
    mustWork = FALSE
  )
  simulation_mode <- o2ipa_as_chr(args$simulation_mode, "invivo")
  time_points <- sort(as_num_vec(fixo2_arg_value(args, "agreement_time_points", "time_points", "25,50,100,200,300,500,700,1000"), c(25, 50, 100, 200, 300, 500, 700, 1000)))
  agreement_o2 <- sort(as_num_vec(fixo2_arg_value(args, "agreement_o2_values", "o2_values", paste(o2_grid, collapse = ",")), o2_grid))
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  initial_ploidy_values <- sort(as_num_vec(fixo2_arg_value(args, "agreement_initial_ploidy_values", "initial_ploidy_values", "2,4"), c(2, 4)))
  simulation_ids <- sort(as_int_vec(fixo2_arg_value(args, "agreement_simulation_ids", "simulation_ids", "1,2,3"), c(1L, 2L, 3L)))
  objective_transform <- o2ipa_as_chr(fixo2_arg_value(args, "agreement_objective_transform", "objective_transform", "identity"), "identity")
  if (!objective_transform %in% c("identity", "log10")) stop("--agreement_objective_transform must be identity or log10.")
  recompute <- o2ipa_as_bool(fixo2_arg_value(args, "agreement_recompute", "recompute", FALSE), FALSE)
  cache_all_times <- o2ipa_as_bool(fixo2_arg_value(args, "agreement_cache_all_times", "cache_all_times", TRUE), TRUE)
  recompute_analytical <- o2ipa_as_bool(fixo2_arg_value(args, "agreement_recompute_analytical", "recompute_analytical", recompute), recompute)
  progress_every <- o2ipa_as_int(fixo2_arg_value(args, "agreement_progress_every", "progress_every", 100L), 100L)
  n_workers <- o2ipa_as_int(
    fixo2_arg_value(args, "agreement_n_workers", "n_workers", Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
    1L
  )
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  analytical_methods <- split_csv(fixo2_arg_value(args, "agreement_analytical_methods", "analytical_methods", "eigen,expm"), default = c("eigen", "expm"))
  analytical_methods <- unique(vapply(analytical_methods, method_slug, character(1)))
  seed_ids <- normalize_seed_ids(split_csv(fixo2_arg_value(args, "agreement_seed_ids", "seed_ids", ""), default = character()))

  legacy_table_dir <- file.path(agreement_analysis_dir, "simulation", "scatters", "tables")
  analytical_cache_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_analytical_cache_table", "analytical_cache_table", file.path(table_dir, "agreement_generated_analytical_trajectories.tsv.gz")),
      file.path(table_dir, "agreement_generated_analytical_trajectories.tsv.gz")
    ),
    repo_root,
    mustWork = FALSE
  )
  all_time_sim_metric_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_all_time_simulation_metric_table", "all_time_simulation_metric_table", file.path(table_dir, "agreement_simulation_all_time_metrics_by_replicate.tsv.gz")),
      file.path(table_dir, "agreement_simulation_all_time_metrics_by_replicate.tsv.gz")
    ),
    repo_root,
    mustWork = FALSE
  )
  sim_metric_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_simulation_metric_table", "simulation_metric_table", file.path(table_dir, "agreement_simulation_selected_time_metrics_by_replicate.tsv")),
      file.path(table_dir, "agreement_simulation_selected_time_metrics_by_replicate.tsv")
    ),
    repo_root,
    mustWork = FALSE
  )
  sim_summary_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_simulation_summary_table", "simulation_summary_table", file.path(table_dir, "agreement_simulation_selected_time_metrics.tsv")),
      file.path(table_dir, "agreement_simulation_selected_time_metrics.tsv")
    ),
    repo_root,
    mustWork = FALSE
  )
  agreement_data_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_data_table", "scatter_data_table", file.path(table_dir, "agreement_analytical_vs_simulation_data.tsv")),
      file.path(table_dir, "agreement_analytical_vs_simulation_data.tsv")
    ),
    repo_root,
    mustWork = FALSE
  )

  analytical_cache_read_path <- fixo2_agreement_cache_read_path(
    analytical_cache_path,
    c(
      file.path(table_dir, "scatter_generated_analytical_trajectories.tsv.gz"),
      file.path(legacy_table_dir, "scatter_generated_analytical_trajectories.tsv.gz")
    ),
    recompute_analytical
  )
  all_time_sim_metric_read_path <- fixo2_agreement_cache_read_path(
    all_time_sim_metric_path,
    c(
      file.path(table_dir, "scatter_simulation_all_time_metrics_by_replicate.tsv.gz"),
      file.path(legacy_table_dir, "scatter_simulation_all_time_metrics_by_replicate.tsv.gz")
    ),
    recompute
  )
  sim_metric_read_path <- fixo2_agreement_cache_read_path(
    sim_metric_path,
    c(
      file.path(table_dir, "scatter_simulation_selected_time_metrics_by_replicate.tsv"),
      file.path(legacy_table_dir, "scatter_simulation_selected_time_metrics_by_replicate.tsv")
    ),
    recompute
  )
  sim_summary_read_path <- fixo2_agreement_cache_read_path(
    sim_summary_path,
    c(
      file.path(table_dir, "scatter_simulation_selected_time_metrics.tsv"),
      file.path(legacy_table_dir, "scatter_simulation_selected_time_metrics.tsv")
    ),
    recompute
  )

  message("run_dir: ", agreement_run_dir)
  message("analysis_dir: ", agreement_analysis_dir)
  message("out_dir: ", out_dir)
  message("simulation_dir: ", agreement_simulation_dir)
  message("simulation_mode: ", simulation_mode)
  message("time_points: ", paste(time_points, collapse = ","))
  message("O2 values: ", paste(agreement_o2, collapse = ","))
  message("mode_reference_o2: ", mode_reference_o2)
  message("initial_ploidy_values: ", paste(initial_ploidy_values, collapse = ","))
  message("simulation_ids: ", paste(simulation_ids, collapse = ","))
  message("analytical_methods: ", paste(analytical_methods, collapse = ","))
  message("n_workers: ", n_workers)

  expected_seed_ids <- expected_analytical_seed_ids(agreement_run_dir, seed_ids)
  message("Generating or reading analytical trajectories.")
  if (!recompute_analytical && file.exists(analytical_cache_read_path)) {
    message("Reading cached analytical trajectories: ", analytical_cache_read_path)
    analytical <- read_tsv(analytical_cache_read_path)
    missing_keys <- analytical_cache_missing_keys(
      analytical = analytical,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      expected_seed_ids = expected_seed_ids,
      expected_run_dir = agreement_run_dir
    )
    analytical <- filter_analytical_trajectories(
      analytical = analytical,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      seed_ids = seed_ids
    )
    missing_methods <- setdiff(analytical_methods, unique(analytical$analytical_method))
    if (length(missing_methods) || length(missing_keys)) analytical <- data.frame()
  } else {
    analytical <- data.frame()
  }
  if (!nrow(analytical)) {
    analytical <- generate_analytical_trajectories(
      run_dir = agreement_run_dir,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      n_workers = n_workers,
      seed_ids = seed_ids
    )
    write_tsv(analytical, analytical_cache_path)
    analytical <- filter_analytical_trajectories(
      analytical = analytical,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      seed_ids = seed_ids
    )
  } else if (!identical(normalizePath(analytical_cache_read_path, mustWork = FALSE), normalizePath(analytical_cache_path, mustWork = FALSE))) {
    write_tsv(analytical, analytical_cache_path)
  }
  if (!nrow(analytical)) stop("No analytical trajectory rows were found for the requested O2/time grid.")

  message("Reading seed objective values.")
  objectives <- read_seed_objectives(
    analysis_dir = agreement_analysis_dir,
    fit_dir = fit_dir,
    o2_values = agreement_o2,
    seed_ids = seed_ids,
    n_workers = n_workers,
    mode_reference_o2 = mode_reference_o2
  )

  if (!cache_all_times && !recompute && file.exists(sim_summary_read_path)) {
    message("Reading cached simulation summary: ", sim_summary_read_path)
    sim_summary <- read_tsv(sim_summary_read_path)
    if (!identical(normalizePath(sim_summary_read_path, mustWork = FALSE), normalizePath(sim_summary_path, mustWork = FALSE))) {
      write_tsv(sim_summary, sim_summary_path)
    }
  } else {
    if (cache_all_times) {
      if (!recompute && file.exists(all_time_sim_metric_read_path)) {
        message("Reading cached all-time simulation metrics: ", all_time_sim_metric_read_path)
        all_time_sim_metrics <- read_tsv(all_time_sim_metric_read_path)
      } else {
        if (!dir.exists(agreement_simulation_dir)) stop("Simulation directory is required when all-time cache is missing: ", agreement_simulation_dir)
        tasks <- read_simulation_tasks(
          simulation_dir = agreement_simulation_dir,
          simulation_mode = simulation_mode,
          analytical = analytical,
          o2_values = agreement_o2,
          initial_ploidy_values = initial_ploidy_values,
          simulation_ids = simulation_ids
        )
        if (!nrow(tasks)) stop("No simulation tasks matched the requested seed/O2/initial ploidy/simulation id filters.")
        message("Matched simulation tasks: ", nrow(tasks))
        all_time_sim_metrics <- read_simulation_metrics(
          tasks = tasks,
          time_points = NULL,
          progress_every = progress_every,
          n_workers = n_workers
        )
        if (!nrow(all_time_sim_metrics)) stop("No all-time simulation metrics were read from state trajectories.")
        write_tsv(all_time_sim_metrics, all_time_sim_metric_path)
      }
      sim_metrics <- filter_simulation_metrics(
        sim_metrics = all_time_sim_metrics,
        time_points = time_points,
        o2_values = agreement_o2,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids,
        seed_ids = unique(analytical$seed_id)
      )
    } else if (!recompute && file.exists(sim_metric_read_path)) {
      message("Reading cached selected-time simulation metrics: ", sim_metric_read_path)
      sim_metrics <- read_tsv(sim_metric_read_path)
      sim_metrics <- filter_simulation_metrics(
        sim_metrics = sim_metrics,
        time_points = time_points,
        o2_values = agreement_o2,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids,
        seed_ids = unique(analytical$seed_id)
      )
    } else {
      if (!dir.exists(agreement_simulation_dir)) stop("Simulation directory is required when selected-time cache is missing: ", agreement_simulation_dir)
      tasks <- read_simulation_tasks(
        simulation_dir = agreement_simulation_dir,
        simulation_mode = simulation_mode,
        analytical = analytical,
        o2_values = agreement_o2,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids
      )
      if (!nrow(tasks)) stop("No simulation tasks matched the requested seed/O2/initial ploidy/simulation id filters.")
      message("Matched simulation tasks: ", nrow(tasks))
      sim_metrics <- read_simulation_metrics(
        tasks = tasks,
        time_points = time_points,
        progress_every = progress_every,
        n_workers = n_workers
      )
    }
    if (!nrow(sim_metrics)) stop("No simulation metrics were read from state trajectories.")
    write_tsv(sim_metrics, sim_metric_path)
    sim_summary <- aggregate_replicates(sim_metrics)
    write_tsv(sim_summary, sim_summary_path)
  }

  message("Merging analytical and simulation summaries.")
  agreement_rows <- list()
  fig_dirs <- character()
  plot_dirs <- character()
  for (method in sort(unique(analytical$analytical_method))) {
    method_analytical <- analytical[analytical$analytical_method %in% method, , drop = FALSE]
    plot_dirs[[method]] <- make_analytical_solution_outputs(
      analytical = method_analytical,
      out_dir = out_dir,
      analytical_method = method,
      seed_metadata = objectives
    )
    agreement_data <- merge_scatter_data(method_analytical, sim_summary, objectives)
    if (!nrow(agreement_data)) {
      warning("No merged analytical-vs-simulation rows were produced for method: ", method)
      next
    }
    method_table_dir <- file.path(table_dir, method_slug(method))
    dir.create(method_table_dir, recursive = TRUE, showWarnings = FALSE)
    method_agreement_data_path <- file.path(method_table_dir, "agreement_analytical_vs_simulation_data.tsv")
    write_tsv(agreement_data, method_agreement_data_path)
    agreement_rows[[method]] <- agreement_data
    message("Drawing analytical-simulation agreement plots for method: ", method)
    fig_dirs[[method]] <- make_scatter_outputs(
      dat = agreement_data,
      out_dir = out_dir,
      objective_transform = objective_transform,
      analytical_method = method
    )
  }
  if (!length(fig_dirs)) stop("No analytical-simulation agreement figures were produced for any analytical method.")
  combined_agreement_data <- do.call(rbind, agreement_rows)
  if (!is.null(combined_agreement_data) && nrow(combined_agreement_data)) write_tsv(combined_agreement_data, agreement_data_path)

  manifest <- data.frame(
    field = c(
      "simulation_dir", "analysis_dir", "fit_dir", "run_dir", "out_dir", "simulation_mode",
      "time_points", "o2_values", "mode_reference_o2", "initial_ploidy_values", "simulation_ids",
      "objective_transform", "cache_all_times", "n_workers", "seed_ids",
      "analytical_methods", "recompute", "recompute_analytical", "analytical_cache_table",
      "all_time_simulation_metric_table", "simulation_metric_table", "simulation_summary_table",
      "agreement_data_table", "figure_dirs", "analytical_solution_plot_dirs"
    ),
    value = c(
      agreement_simulation_dir, agreement_analysis_dir, fit_dir, agreement_run_dir, out_dir, simulation_mode,
      paste(time_points, collapse = ","), paste(agreement_o2, collapse = ","), mode_reference_o2,
      paste(initial_ploidy_values, collapse = ","), paste(simulation_ids, collapse = ","),
      objective_transform, as.character(cache_all_times), as.character(n_workers), paste(seed_ids, collapse = ","),
      paste(analytical_methods, collapse = ","), as.character(recompute), as.character(recompute_analytical), analytical_cache_path,
      all_time_sim_metric_path, sim_metric_path, sim_summary_path, agreement_data_path, paste(fig_dirs, collapse = ","), paste(plot_dirs, collapse = ",")
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(table_dir, "analytical_simulation_agreement_run_manifest.tsv"))
  message("Completed analytical-simulation agreement analysis: ", out_dir)
  invisible(list(out_dir = out_dir, figure_dirs = fig_dirs, analytical_methods = analytical_methods))
}

fixo2_main <- function(args = o2ipa_parse_args()) {
  repo_root <- o2ipa_repo_root(SCRIPT_DIR)
  run_dir <- fixo2_resolve_repo_path(o2ipa_as_chr(args$run_dir, fixo2_default_run_dir(repo_root)), repo_root, mustWork = FALSE)
  out_dir <- fixo2_resolve_repo_path(o2ipa_as_chr(args$out_dir, file.path(repo_root, "oxygen", "results", "analysis", "FixO2_invivo_500seed")), repo_root, mustWork = FALSE)
  label_file <- fixo2_resolve_repo_path(o2ipa_as_chr(args$label_file, fixo2_default_label_file(repo_root)), repo_root, mustWork = FALSE)
  simulation_dir <- fixo2_resolve_repo_path(o2ipa_as_chr(args$simulation_dir, fixo2_default_simulation_dir(repo_root)), repo_root, mustWork = FALSE)
  o2_grid <- sort(unique(fixo2_as_num_vec(args$o2_grid, c(0, 0.1, 0.5, 1, 2, 5))))
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  parts <- fixo2_normalize_parts(args$run_parts)
  fixo2_prepare_dirs(out_dir)

  render_html_report <- o2ipa_as_bool(fixo2_arg_value(args, "render_html_report", "generate_html_report", TRUE), TRUE)
  html_report_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "html_report_dir", "report_out_dir", file.path(out_dir, "html_report")), file.path(out_dir, "html_report")),
    repo_root,
    mustWork = FALSE
  )
  html_report_basename <- o2ipa_as_chr(fixo2_arg_value(args, "html_report_basename", "report_basename", "index"), "index")

  top_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "optional_label_file", "mode_source", "mode_rule", "mode_reference_o2",
      "simulation_dir", "o2_grid", "run_parts", "simulation_n_core", "simulation_worker_threads",
      "render_html_report", "html_report_dir", "html_report_basename"
    ),
    value = c(
      run_dir, out_dir, label_file,
      "fixed_o2_attractor_dominant_ploidy_at_reference_o2",
      paste0(
        "downstream analyses use the seed-level mode assigned at fixed O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        "; the complete attractor output still reports O2-dependent modes for all attractor O2 values"
      ),
      mode_reference_o2,
      simulation_dir, paste(o2_grid, collapse = ","), paste(parts, collapse = ","),
      as.character(fixo2_simulation_n_core(args)), as.character(fixo2_simulation_worker_threads(args)),
      as.character(render_html_report), html_report_dir, html_report_basename
    ),
    stringsAsFactors = FALSE
  )
  fixo2_write_tsv(top_args, file.path(out_dir, "tables", "FixO2_invivo_run_arguments.tsv"))

  results <- list()
  if ("attractors" %in% parts) {
    results$attractors <- fixo2_run_attractors(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      label_file = label_file,
      out_dir = file.path(out_dir, "attractors")
    )
  }
  if ("counterfactual_trajectories" %in% parts) {
    results$counterfactual_trajectories <- fixo2_run_counterfactual(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      label_file = label_file,
      out_dir = file.path(out_dir, "counterfactual_trajectories"),
      o2_grid = o2_grid
    )
  }
  if ("simulation" %in% parts) {
    results$simulation <- fixo2_run_simulation_validation(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      out_dir = file.path(out_dir, "simulation"),
      simulation_dir = simulation_dir,
      o2_grid = o2_grid
    )
  }
  if ("analytical_simulation_agreement" %in% parts) {
    results$analytical_simulation_agreement <- fixo2_run_analytical_simulation_agreement(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      analysis_dir = out_dir,
      out_dir = file.path(out_dir, "simulation", "analytical_simulation_agreement"),
      simulation_dir = simulation_dir,
      o2_grid = o2_grid
    )
  }

  html_report <- fixo2_render_html_report(args, repo_root = repo_root, out_dir = out_dir)
  if (!is.null(html_report)) results$html_report <- html_report
  fixo2_write_output_manifest(out_dir)
  message("Completed FixO2 in vivo workflow: ", out_dir)
  invisible(results)
}

if (identical(environment(), globalenv())) {
  fixo2_main()
}
