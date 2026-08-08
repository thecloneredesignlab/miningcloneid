# Pure helpers for the in-vitro observed-data feasibility scan.
#
# This file intentionally has no package dependencies and performs no work when
# sourced.  The scan entrypoint owns model evaluation, file IO, and parallelism.

.ivt_scan_assert_scalar_number <- function(x, name, positive = FALSE) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || !is.finite(x) || (positive && x <= 0)) {
    qualifier <- if (positive) " one finite positive number" else " one finite number"
    stop(name, " must be", qualifier, ".", call. = FALSE)
  }
  x
}

.ivt_scan_restore_rng <- function(had_seed, old_seed, old_kind) {
  do.call(RNGkind, as.list(old_kind))
  if (had_seed) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

ivt_scan_lhs_design <- function(n, specs, seed) {
  n_num <- .ivt_scan_assert_scalar_number(n, "n", positive = TRUE)
  if (n_num > .Machine$integer.max || n_num != as.integer(n_num)) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n_num)
  seed_num <- .ivt_scan_assert_scalar_number(seed, "seed")
  if (seed_num != as.integer(seed_num) || seed_num < 0 ||
      seed_num > .Machine$integer.max) {
    stop("seed must be an integer between 0 and .Machine$integer.max.", call. = FALSE)
  }
  if (!is.data.frame(specs)) {
    stop("specs must be a data frame.", call. = FALSE)
  }
  required <- c("param_symbol", "lower", "upper", "transform")
  missing <- setdiff(required, names(specs))
  if (length(missing)) {
    stop("specs is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!nrow(specs)) {
    stop("specs must contain at least one parameter.", call. = FALSE)
  }

  symbols <- as.character(specs$param_symbol)
  if (anyNA(symbols) || any(!nzchar(symbols)) || anyDuplicated(symbols)) {
    stop("specs$param_symbol must contain unique non-empty names.", call. = FALSE)
  }
  lower <- suppressWarnings(as.numeric(specs$lower))
  upper <- suppressWarnings(as.numeric(specs$upper))
  transforms <- tolower(as.character(specs$transform))
  transforms[transforms == "log10"] <- "log"
  if (any(!is.finite(lower)) || any(!is.finite(upper)) || any(lower > upper)) {
    stop("Every parameter must have finite bounds with lower <= upper.", call. = FALSE)
  }
  if (any(!transforms %in% c("linear", "log"))) {
    stop("specs$transform must contain only 'linear' or 'log'.", call. = FALSE)
  }
  if (any(transforms == "log" & lower <= 0)) {
    stop("Log-transformed parameters require strictly positive bounds.", call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  on.exit(.ivt_scan_restore_rng(had_seed, old_seed, old_kind), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed_num))

  values <- matrix(NA_real_, nrow = n, ncol = nrow(specs))
  for (j in seq_len(nrow(specs))) {
    if (lower[[j]] == upper[[j]]) {
      values[, j] <- lower[[j]]
      next
    }
    strata <- ((seq_len(n) - 1L) + stats::runif(n)) / n
    u <- strata[sample.int(n, size = n, replace = FALSE)]
    if (transforms[[j]] == "log") {
      values[, j] <- exp(log(lower[[j]]) + u * (log(upper[[j]]) - log(lower[[j]])))
    } else {
      values[, j] <- lower[[j]] + u * (upper[[j]] - lower[[j]])
    }
  }
  colnames(values) <- symbols
  id_width <- max(4L, nchar(as.character(n)))
  data.frame(
    candidate_id = sprintf(paste0("candidate_%0", id_width, "d"), seq_len(n)),
    values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

ivt_scan_parameter_specs <- function(parameter_table, p_mis_base_upper) {
  if (is.character(parameter_table) && length(parameter_table) == 1L) {
    parameter_table <- utils::read.csv(
      parameter_table,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (!is.data.frame(parameter_table)) {
    stop("parameter_table must be a data frame or a CSV path.", call. = FALSE)
  }
  required <- c("param_symbol", "lower_bound", "upper_bound")
  missing <- setdiff(required, names(parameter_table))
  if (length(missing)) {
    stop(
      "parameter_table is missing column(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  symbols <- c(
    "p_mis_base", "p_misseg", "k_o_mis", "gamma_mu",
    "buffer_smax", "buffer_beta", "buffer_n_exp", "O2_crit", "n_O"
  )
  transforms <- c("log", "log", "log", "linear", "linear", "log", "log", "log", "linear")
  table_symbols <- as.character(parameter_table$param_symbol)
  duplicated_targets <- symbols[vapply(
    symbols,
    function(x) sum(table_symbols == x, na.rm = TRUE) != 1L,
    logical(1)
  )]
  if (length(duplicated_targets)) {
    stop(
      "parameter_table must contain exactly one row for: ",
      paste(duplicated_targets, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  idx <- match(symbols, table_symbols)
  lower <- suppressWarnings(as.numeric(parameter_table$lower_bound[idx]))
  source_upper <- suppressWarnings(as.numeric(parameter_table$upper_bound[idx]))
  if (any(!is.finite(lower)) || any(!is.finite(source_upper)) || any(lower > source_upper)) {
    stop("Scan parameter bounds must be finite with lower_bound <= upper_bound.", call. = FALSE)
  }

  calibrated_upper <- .ivt_scan_assert_scalar_number(
    p_mis_base_upper,
    "p_mis_base_upper",
    positive = TRUE
  )
  p_idx <- match("p_mis_base", symbols)
  tolerance <- 1e-12 * max(1, abs(source_upper[[p_idx]]))
  if (calibrated_upper < lower[[p_idx]] - tolerance ||
      calibrated_upper > source_upper[[p_idx]] + tolerance) {
    stop(
      "p_mis_base_upper must lie within the parameter-table p_mis_base bounds.",
      call. = FALSE
    )
  }
  if (abs(calibrated_upper - lower[[p_idx]]) <= tolerance) {
    calibrated_upper <- lower[[p_idx]]
  }
  if (abs(calibrated_upper - source_upper[[p_idx]]) <= tolerance) {
    calibrated_upper <- source_upper[[p_idx]]
  }
  upper <- source_upper
  upper[[p_idx]] <- min(max(calibrated_upper, lower[[p_idx]]), source_upper[[p_idx]])
  if (any(transforms == "log" & lower <= 0)) {
    stop("Log-transformed scan parameters require strictly positive bounds.", call. = FALSE)
  }

  data.frame(
    param_symbol = symbols,
    lower = lower,
    upper = upper,
    transform = transforms,
    stringsAsFactors = FALSE
  )
}

.ivt_scan_density_grid <- function(df) {
  if (!is.data.frame(df) || !all(c("ploidy", "density") %in% names(df))) {
    stop("df must contain ploidy and density columns.", call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(df$ploidy))
  y <- suppressWarnings(as.numeric(df$density))
  if (length(x) < 2L || any(!is.finite(x)) || any(!is.finite(y)) || any(y < 0)) {
    stop("A density grid needs at least two finite points and non-negative density.", call. = FALSE)
  }
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  dx_values <- diff(x)
  if (any(dx_values <= 0)) {
    stop("Density-grid ploidy values must be unique and increasing.", call. = FALSE)
  }
  dx <- stats::median(dx_values)
  tolerance <- max(1, abs(dx)) * 1e-8
  if (any(abs(dx_values - dx) > tolerance)) {
    stop("Density-grid ploidy values must be equally spaced.", call. = FALSE)
  }
  total <- sum(y) * dx
  if (!is.finite(total) || total <= 0) {
    stop("Density-grid total mass must be strictly positive.", call. = FALSE)
  }
  list(x = x, y = y, dx = dx, total = total)
}

.ivt_scan_same_value <- function(x, value) {
  if (is.na(value)) is.na(x) else !is.na(x) & x == value
}

ivt_scan_flow_group_metrics <- function(df, lower = 3) {
  lower <- .ivt_scan_assert_scalar_number(lower, "lower")
  required <- c("segment_id", "passage_id", "sample_name", "series", "ploidy", "density")
  missing <- setdiff(required, names(df))
  if (!is.data.frame(df) || length(missing)) {
    stop("df is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  keep <- as.character(df$series) == "Predicted" &
    as.character(df$segment_id) == "2N-C-A12"
  pred <- df[!is.na(keep) & keep, , drop = FALSE]
  empty <- data.frame(
    segment_id = character(), passage_id = character(), sample_name = character(),
    n_grid = integer(), grid_step = numeric(), total_mass = numeric(),
    high_mass_ge_lower = numeric(), stringsAsFactors = FALSE
  )
  if (!nrow(pred)) return(empty)

  keys <- unique(pred[c("segment_id", "passage_id", "sample_name")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    in_group <- .ivt_scan_same_value(pred$segment_id, key$segment_id[[1L]]) &
      .ivt_scan_same_value(pred$passage_id, key$passage_id[[1L]]) &
      .ivt_scan_same_value(pred$sample_name, key$sample_name[[1L]])
    grid <- .ivt_scan_density_grid(pred[in_group, , drop = FALSE])
    data.frame(
      segment_id = as.character(key$segment_id[[1L]]),
      passage_id = as.character(key$passage_id[[1L]]),
      sample_name = as.character(key$sample_name[[1L]]),
      n_grid = length(grid$x),
      grid_step = grid$dx,
      total_mass = grid$total,
      high_mass_ge_lower = sum(grid$y[grid$x >= lower]) * grid$dx / grid$total,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.ivt_scan_interval_mass <- function(x, probability, lower, upper = Inf) {
  keep <- x >= lower & x < upper
  sum(probability[keep])
}

ivt_scan_flow_shape_metrics <- function(df) {
  grid <- .ivt_scan_density_grid(df)
  x <- grid$x
  probability <- grid$y * grid$dx / grid$total
  mean_ploidy <- sum(x * probability)
  sd_ploidy <- sqrt(max(0, sum((x - mean_ploidy)^2 * probability)))

  local_idx <- integer()
  if (length(x) >= 3L) {
    local_idx <- which(
      probability[2:(length(x) - 1L)] > probability[1:(length(x) - 2L)] &
        probability[2:(length(x) - 1L)] >= probability[3:length(x)]
    ) + 1L
  }
  low_candidates <- local_idx[x[local_idx] >= 1.4 & x[local_idx] <= 2.8]
  high_candidates <- local_idx[x[local_idx] >= 3.3 & x[local_idx] <= 5.3]
  low_idx <- if (length(low_candidates)) {
    low_candidates[[which.max(probability[low_candidates])]]
  } else {
    NA_integer_
  }
  high_idx <- if (length(high_candidates)) {
    high_candidates[[which.max(probability[high_candidates])]]
  } else {
    NA_integer_
  }
  low_peak <- if (!is.na(low_idx)) x[[low_idx]] else NA_real_
  high_peak <- if (!is.na(high_idx)) x[[high_idx]] else NA_real_
  peak_ratio <- if (!is.na(low_idx) && !is.na(high_idx) && low_peak > 0) {
    high_peak / low_peak
  } else {
    NA_real_
  }
  valley_ratio <- NA_real_
  high_relative_height <- NA_real_
  if (!is.na(low_idx) && !is.na(high_idx) && high_idx > low_idx) {
    between <- seq.int(low_idx, high_idx)
    valley_ratio <- min(probability[between]) /
      min(probability[[low_idx]], probability[[high_idx]])
    high_relative_height <- probability[[high_idx]] / max(probability)
  }
  has_target_bimodality <- is.finite(peak_ratio) && peak_ratio >= 1.75 &&
    peak_ratio <= 2.25 && is.finite(valley_ratio) && valley_ratio <= 0.5 &&
    is.finite(high_relative_height) && high_relative_height >= 0.08

  data.frame(
    mean_ploidy = mean_ploidy,
    sd_ploidy = sd_ploidy,
    bridge_mass = .ivt_scan_interval_mass(x, probability, 2.5, 3.5),
    wgd_mass = .ivt_scan_interval_mass(x, probability, 3.5, Inf),
    mass_3_5 = .ivt_scan_interval_mass(x, probability, 3, 5),
    tail_gt4_5 = .ivt_scan_interval_mass(x, probability, 4.5, Inf),
    tail_gt5 = .ivt_scan_interval_mass(x, probability, 5, Inf),
    low_peak = low_peak,
    high_peak = high_peak,
    peak_ratio = peak_ratio,
    valley_ratio = valley_ratio,
    high_relative_height = high_relative_height,
    has_target_bimodality = has_target_bimodality,
    stringsAsFactors = FALSE
  )
}

.ivt_scan_segment_id <- function(seg_res) {
  if (!is.null(seg_res$segment$segment_id)) as.character(seg_res$segment$segment_id) else NA_character_
}

.ivt_scan_weighted_quantile <- function(x, probability, target = 0.5) {
  x[[which(cumsum(probability) >= target)[[1L]]]]
}

ivt_scan_segment_distribution_metrics <- function(run, segment_id, sigma_kary = NULL) {
  if (!is.list(run) || is.null(run$grid_pre) || is.null(run$segment_results)) {
    stop("run must contain grid_pre and segment_results.", call. = FALSE)
  }
  if (length(segment_id) != 1L || is.na(segment_id) || !nzchar(segment_id)) {
    stop("segment_id must be one non-empty string.", call. = FALSE)
  }
  ids <- vapply(run$segment_results, .ivt_scan_segment_id, character(1))
  idx <- which(ids == as.character(segment_id))
  if (length(idx) != 1L) {
    stop("Expected exactly one segment result for segment_id: ", segment_id, ".", call. = FALSE)
  }
  seg_res <- run$segment_results[[idx]]
  grid <- suppressWarnings(as.numeric(run$grid_pre))
  probability <- suppressWarnings(as.numeric(seg_res$selection$selected_frac))
  if (!length(grid) || length(grid) != length(probability) || any(!is.finite(grid)) ||
      any(!is.finite(probability)) || any(probability < 0)) {
    stop("The selected segment distribution is invalid.", call. = FALSE)
  }
  total <- sum(probability)
  if (!is.finite(total) || total <= 0) {
    stop("The selected segment distribution has no positive mass.", call. = FALSE)
  }
  probability <- probability / total
  basis <- "latent"
  sigma_value <- NA_real_
  if (!is.null(sigma_kary)) {
    sigma_value <- .ivt_scan_assert_scalar_number(sigma_kary, "sigma_kary", positive = TRUE)
    kernel <- outer(
      X = grid,
      Y = grid,
      FUN = function(observed_n, true_n) {
        stats::dnorm(observed_n, mean = true_n, sd = sigma_value)
      }
    )
    kernel_sums <- colSums(kernel)
    kernel_sums[!is.finite(kernel_sums) | kernel_sums <= 0] <- 1
    kernel <- sweep(kernel, 2L, kernel_sums, "/", check.margin = FALSE)
    probability <- as.numeric(kernel %*% probability)
    probability <- probability / sum(probability)
    basis <- "observation_smoothed"
  }

  ord <- order(grid)
  grid <- grid[ord]
  probability <- probability[ord]
  mean_n <- sum(grid * probability)
  segment <- seg_res$segment
  field <- function(name, default) {
    value <- segment[[name]]
    if (is.null(value) || !length(value)) default else value[[1L]]
  }
  data.frame(
    segment_id = as.character(segment_id),
    cohort = as.character(field("cohort", NA_character_)),
    lineage_id = as.character(field("lineage_id", NA_character_)),
    scenario_id = as.character(field("scenario_id", NA_character_)),
    oxygen_pct = suppressWarnings(as.numeric(field("oxygen_pct", NA_real_))),
    distribution_basis = basis,
    sigma_kary = sigma_value,
    mean_N = mean_n,
    sd_N = sqrt(max(0, sum((grid - mean_n)^2 * probability))),
    median_N = .ivt_scan_weighted_quantile(grid, probability, 0.5),
    mass_N_ge_70 = sum(probability[grid >= 70]),
    mass_N_ge_80 = sum(probability[grid >= 80]),
    mass_N_ge_88 = sum(probability[grid >= 88]),
    n_grid_states = length(grid),
    stringsAsFactors = FALSE
  )
}

.ivt_scan_first_existing_column <- function(df, candidates, label) {
  found <- intersect(candidates, names(df))
  if (!length(found)) {
    stop(label, " column is required; accepted names: ", paste(candidates, collapse = ", "), ".", call. = FALSE)
  }
  found[[1L]]
}

ivt_scan_select_calibration_cap <- function(df, limit) {
  if (!is.data.frame(df) || !"p_mis_base" %in% names(df)) {
    stop("df must contain p_mis_base.", call. = FALSE)
  }
  limit <- .ivt_scan_assert_scalar_number(limit, "limit")
  if (limit < 0 || limit > 1) {
    stop("limit must lie in [0, 1].", call. = FALSE)
  }
  mass_col <- .ivt_scan_first_existing_column(
    df,
    c(
      "high_mass_ge_lower", "control_high_mass", "control_high_mass_ge_3",
      "control_flow_mass_ge_3", "control_flow_mass_ge3"
    ),
    "Control high-mass"
  )
  p <- suppressWarnings(as.numeric(df$p_mis_base))
  mass <- suppressWarnings(as.numeric(df[[mass_col]]))
  if (!length(p) || any(!is.finite(p)) || any(p <= 0)) {
    stop("p_mis_base values must be finite and strictly positive.", call. = FALSE)
  }

  unique_p <- sort(unique(p))
  summarized_mass <- numeric(length(unique_p))
  valid <- logical(length(unique_p))
  for (i in seq_along(unique_p)) {
    rows <- which(p == unique_p[[i]])
    row_valid <- all(is.finite(mass[rows]))
    if ("status" %in% names(df)) {
      row_valid <- row_valid && all(as.character(df$status[rows]) == "OK")
    }
    if ("protocol_feasibility_status" %in% names(df)) {
      row_valid <- row_valid && all(as.character(df$protocol_feasibility_status[rows]) == "PASS")
    }
    valid[[i]] <- row_valid
    summarized_mass[[i]] <- if (row_valid) max(mass[rows]) else NA_real_
  }
  passes <- valid & summarized_mass <= limit
  prefix <- cumprod(as.integer(passes)) == 1L
  n_prefix <- sum(prefix)
  first_fail_idx <- which(!passes)[1L]
  finite_mass <- all(is.finite(summarized_mass))
  monotone <- finite_mass && all(diff(summarized_mass) >= -sqrt(.Machine$double.eps))
  reentry <- if (is.na(first_fail_idx)) {
    FALSE
  } else {
    any(passes[seq.int(first_fail_idx, length(passes))])
  }

  data.frame(
    p_mis_base_upper = if (n_prefix) unique_p[[n_prefix]] else NA_real_,
    limit = limit,
    n_evaluated = length(unique_p),
    n_prefix_pass = n_prefix,
    first_failing_p_mis_base = if (is.na(first_fail_idx)) NA_real_ else unique_p[[first_fail_idx]],
    monotone_nondecreasing = monotone,
    pass_reentry_after_failure = reentry,
    all_pass = all(passes),
    stringsAsFactors = FALSE
  )
}

.ivt_scan_zero_count_gate <- function(df, candidates, label) {
  found <- intersect(candidates, names(df))
  if (!length(found)) {
    stop(label, " count column is required; accepted names: ", paste(candidates, collapse = ", "), ".", call. = FALSE)
  }
  Reduce(`&`, lapply(found, function(name) {
    value <- suppressWarnings(as.numeric(df[[name]]))
    is.finite(value) & value == 0
  }))
}

ivt_scan_hard_gate <- function(row, control_limit = 0.05) {
  if (is.list(row) && !is.data.frame(row)) row <- as.data.frame(row, stringsAsFactors = FALSE)
  if (!is.data.frame(row) || !nrow(row)) {
    stop("row must be a non-empty data frame or named list.", call. = FALSE)
  }
  required <- c("status", "protocol_feasibility_status", "n_scenarios", "n_insufficient_boundaries")
  missing <- setdiff(required, names(row))
  if (length(missing)) {
    stop("row is missing hard-gate column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  control_limit <- .ivt_scan_assert_scalar_number(control_limit, "control_limit")
  if (control_limit < 0 || control_limit > 1) {
    stop("control_limit must lie in [0, 1].", call. = FALSE)
  }
  control_col <- .ivt_scan_first_existing_column(
    row,
    c(
      "control_high_mass", "control_high_mass_ge_3", "control_flow_mass_ge_3",
      "control_flow_mass_ge3", "high_mass_ge_lower"
    ),
    "Control high-mass"
  )
  control_mass <- suppressWarnings(as.numeric(row[[control_col]]))
  n_scenarios <- suppressWarnings(as.numeric(row$n_scenarios))
  n_insufficient <- suppressWarnings(as.numeric(row$n_insufficient_boundaries))
  no_missing <- .ivt_scan_zero_count_gate(
    row,
    c("n_missing_predictions", "n_growth_missing_pred"),
    "Missing-prediction"
  )
  no_negative <- .ivt_scan_zero_count_gate(
    row,
    c("n_negative_predictions", "n_growth_negative_pred"),
    "Negative-prediction"
  )
  gate <- !is.na(row$status) & as.character(row$status) == "OK" &
    !is.na(row$protocol_feasibility_status) &
    as.character(row$protocol_feasibility_status) == "PASS" &
    is.finite(n_scenarios) & n_scenarios == 6 &
    is.finite(n_insufficient) & n_insufficient == 0 &
    no_missing & no_negative &
    is.finite(control_mass) & control_mass <= control_limit
  if ("all_passage_boundaries_feasible" %in% names(row)) {
    gate <- gate & !is.na(row$all_passage_boundaries_feasible) &
      as.logical(row$all_passage_boundaries_feasible)
  }
  for (flag in c(
    "all_passage_executed",
    "all_selected_at_or_after_last_observation",
    "all_selected_within_tolerance"
  )) {
    if (flag %in% names(row)) {
      gate <- gate & !is.na(row[[flag]]) & as.logical(row[[flag]])
    }
  }
  as.logical(gate)
}

ivt_scan_pareto_front <- function(df, metric_cols) {
  if (!is.data.frame(df)) {
    stop("df must be a data frame.", call. = FALSE)
  }
  metric_cols <- as.character(metric_cols)
  if (!length(metric_cols) || anyNA(metric_cols) || any(!nzchar(metric_cols)) ||
      anyDuplicated(metric_cols)) {
    stop("metric_cols must contain unique non-empty column names.", call. = FALSE)
  }
  missing <- setdiff(metric_cols, names(df))
  if (length(missing)) {
    stop("df is missing metric column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  metrics <- do.call(cbind, lapply(metric_cols, function(name) {
    suppressWarnings(as.numeric(df[[name]]))
  }))
  if (!is.matrix(metrics)) metrics <- matrix(metrics, ncol = length(metric_cols))
  valid <- apply(metrics, 1L, function(x) all(is.finite(x)))
  result <- rep(FALSE, nrow(df))
  valid_idx <- which(valid)
  if (!length(valid_idx)) return(result)

  ord <- valid_idx[order(metrics[valid_idx, 1L], valid_idx)]
  front <- integer()
  for (idx in ord) {
    dominated <- FALSE
    if (length(front)) {
      dominated <- any(vapply(front, function(other) {
        all(metrics[other, ] <= metrics[idx, ]) && any(metrics[other, ] < metrics[idx, ])
      }, logical(1)))
    }
    if (dominated) next
    if (length(front)) {
      keep <- !vapply(front, function(other) {
        all(metrics[idx, ] <= metrics[other, ]) && any(metrics[idx, ] < metrics[other, ])
      }, logical(1))
      front <- front[keep]
    }
    front <- c(front, idx)
  }
  result[front] <- TRUE
  result
}

# ----------------------------------------------------------------------------
# Extended observed-data refinement helpers
# ----------------------------------------------------------------------------

.ivt_scan_extended_parameter_names <- c(
  "p_wgd", "p_mis_base", "gamma_growth", "mu_hp", "gamma_mu",
  "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
  "buffer_n_exp", "O2_crit", "n_O"
)

ivt_scan_extended_parameter_specs <- function(region = c("focused", "outer"),
                                              parameter_table = NULL) {
  region <- match.arg(region)
  if (!is.null(parameter_table) && !is.data.frame(parameter_table) &&
      !(is.character(parameter_table) && length(parameter_table) == 1L)) {
    stop("parameter_table must be NULL, a data frame, or one CSV path.", call. = FALSE)
  }

  focused_lower <- c(
    2e-4, 1e-7, 0.01, 2e-4, 0.3, 0.001, 1e-5, 0.90, 0.2, 4, 1e-4, 0.5
  )
  focused_upper <- c(
    2.5e-3, 5e-3, 0.35, 0.015, 3.5, 0.03, 0.005, 0.9999, 2, 10, 2.5, 5
  )
  outer_lower <- c(
    1e-4, 1e-7, 0.01, 2e-4, 0.1, 0.001, 1e-5, 0.80, 0.2, 4, 1e-4, 0.5
  )
  outer_upper <- c(
    2.5e-3, 5e-3, 0.50, 0.020, 3.5, 0.03, 0.020, 0.9999, 2, 10, 2.5, 5
  )
  transforms <- c(
    "log", "conditional_log", "log", "log", "linear", "log", "log",
    "complement_log", "log", "log", "log", "linear"
  )

  data.frame(
    param_symbol = .ivt_scan_extended_parameter_names,
    lower = if (region == "focused") focused_lower else outer_lower,
    upper = if (region == "focused") focused_upper else outer_upper,
    transform = transforms,
    region = region,
    stringsAsFactors = FALSE
  )
}

ivt_scan_select_conditional_envelope <- function(calibration, control_limit = 0.05) {
  if (!is.data.frame(calibration) || !"p_wgd" %in% names(calibration)) {
    stop("calibration must be a data frame containing p_wgd.", call. = FALSE)
  }
  p_wgd <- suppressWarnings(as.numeric(calibration$p_wgd))
  if (!length(p_wgd) || any(!is.finite(p_wgd)) || any(p_wgd <= 0)) {
    stop("calibration$p_wgd must be finite and strictly positive.", call. = FALSE)
  }
  groups <- split(seq_len(nrow(calibration)), p_wgd)
  rows <- lapply(groups, function(idx) {
    cap <- ivt_scan_select_calibration_cap(
      calibration[idx, , drop = FALSE],
      control_limit
    )
    data.frame(
      p_wgd = unique(p_wgd[idx])[[1L]],
      cap,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  envelope <- do.call(rbind, rows)
  envelope <- envelope[order(envelope$p_wgd), , drop = FALSE]
  rownames(envelope) <- NULL
  local_eligible <- is.finite(envelope$p_mis_base_upper) &
    envelope$n_prefix_pass > 0L
  envelope$eligible_local <- local_eligible
  # A higher-p_wgd re-entry is not part of the admissible contiguous envelope.
  envelope$eligible <- cumprod(as.integer(local_eligible)) == 1L
  envelope$cross_p_wgd_reentry <- !envelope$eligible & local_eligible
  envelope
}

.ivt_scan_conditional_cap <- function(p_wgd, envelope) {
  required <- c("p_wgd", "p_mis_base_upper", "eligible")
  if (!is.data.frame(envelope) || !all(required %in% names(envelope))) {
    stop(
      "envelope must contain p_wgd, p_mis_base_upper, and eligible.",
      call. = FALSE
    )
  }
  nodes <- envelope[
    !is.na(envelope$eligible) & as.logical(envelope$eligible) &
      is.finite(as.numeric(envelope$p_wgd)) &
      is.finite(as.numeric(envelope$p_mis_base_upper)),
    ,
    drop = FALSE
  ]
  if (!nrow(nodes)) {
    stop("The conditional envelope has no eligible node.", call. = FALSE)
  }
  nodes$p_wgd <- as.numeric(nodes$p_wgd)
  nodes$p_mis_base_upper <- as.numeric(nodes$p_mis_base_upper)
  nodes <- nodes[order(nodes$p_wgd), , drop = FALSE]
  if (anyDuplicated(nodes$p_wgd)) {
    stop("Eligible envelope p_wgd nodes must be unique.", call. = FALSE)
  }

  values <- suppressWarnings(as.numeric(p_wgd))
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop("p_wgd values must be finite and strictly positive.", call. = FALSE)
  }
  vapply(values, function(value) {
    if (value <= nodes$p_wgd[[1L]]) return(nodes$p_mis_base_upper[[1L]])
    if (value >= nodes$p_wgd[[nrow(nodes)]]) {
      return(nodes$p_mis_base_upper[[nrow(nodes)]])
    }
    right <- which(nodes$p_wgd >= value)[[1L]]
    left <- right - 1L
    # Use the smaller adjacent cap rather than interpolating upward.  The full
    # candidate control replay remains the authoritative hard gate.
    min(nodes$p_mis_base_upper[c(left, right)])
  }, numeric(1))
}

.ivt_scan_lhs_unit <- function(n, d, seed) {
  n_num <- .ivt_scan_assert_scalar_number(n, "n")
  d_num <- .ivt_scan_assert_scalar_number(d, "d")
  seed_num <- .ivt_scan_assert_scalar_number(seed, "seed")
  if (n_num < 0 || n_num != as.integer(n_num) ||
      d_num <= 0 || d_num != as.integer(d_num) ||
      seed_num < 0 || seed_num != as.integer(seed_num)) {
    stop("n, d, and seed must be non-negative/positive integer values.", call. = FALSE)
  }
  n <- as.integer(n_num)
  d <- as.integer(d_num)
  if (!n) return(matrix(numeric(), nrow = 0L, ncol = d))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  on.exit(.ivt_scan_restore_rng(had_seed, old_seed, old_kind), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed_num))

  out <- matrix(NA_real_, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    strata <- ((seq_len(n) - 1L) + stats::runif(n)) / n
    out[, j] <- strata[sample.int(n, n, replace = FALSE)]
  }
  out
}

.ivt_scan_from_unit <- function(u, lower, upper, transform) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  u <- pmin(pmax(as.numeric(u), 0), 1)
  if (!is.finite(lower) || !is.finite(upper) || lower > upper) {
    stop("Invalid transformed bounds.", call. = FALSE)
  }
  if (lower == upper) return(rep(lower, length(u)))
  if (transform %in% c("log", "conditional_log")) {
    if (lower <= 0 || upper <= 0) stop("Log bounds must be positive.", call. = FALSE)
    return(exp(log(lower) + u * (log(upper) - log(lower))))
  }
  if (transform == "complement_log") {
    if (lower < 0 || upper >= 1) {
      stop("complement_log bounds must satisfy 0 <= lower <= upper < 1.", call. = FALSE)
    }
    gap_low <- 1 - upper
    gap_high <- 1 - lower
    return(1 - exp(log(gap_low) + u * (log(gap_high) - log(gap_low))))
  }
  if (transform == "linear") return(lower + u * (upper - lower))
  stop("Unsupported transform: ", transform, call. = FALSE)
}

.ivt_scan_to_coordinate <- function(x, transform) {
  x <- as.numeric(x)
  if (transform %in% c("log", "conditional_log")) return(log(x))
  if (transform == "complement_log") return(log(1 - x))
  if (transform == "linear") return(x)
  stop("Unsupported transform: ", transform, call. = FALSE)
}

.ivt_scan_region_design <- function(n, specs, seed, envelope, region) {
  n <- as.integer(n)
  if (n < 0L) stop("n must be non-negative.", call. = FALSE)
  if (!n) return(data.frame())
  eligible <- envelope[as.logical(envelope$eligible), , drop = FALSE]
  if (!nrow(eligible)) stop("The conditional envelope has no eligible p_wgd range.", call. = FALSE)
  p_idx <- match("p_wgd", specs$param_symbol)
  specs$lower[[p_idx]] <- max(specs$lower[[p_idx]], min(eligible$p_wgd))
  specs$upper[[p_idx]] <- min(specs$upper[[p_idx]], max(eligible$p_wgd))
  if (specs$lower[[p_idx]] > specs$upper[[p_idx]]) {
    stop("The conditional envelope does not overlap the requested p_wgd range.", call. = FALSE)
  }

  u <- .ivt_scan_lhs_unit(n, nrow(specs), seed)
  out <- data.frame(matrix(NA_real_, nrow = n, ncol = nrow(specs)))
  names(out) <- specs$param_symbol
  for (j in seq_len(nrow(specs))) {
    name <- specs$param_symbol[[j]]
    if (name == "p_mis_base") next
    out[[name]] <- .ivt_scan_from_unit(
      u[, j], specs$lower[[j]], specs$upper[[j]], specs$transform[[j]]
    )
  }
  cap <- .ivt_scan_conditional_cap(out$p_wgd, envelope)
  mis_idx <- match("p_mis_base", specs$param_symbol)
  cap <- pmin(cap, specs$upper[[mis_idx]])
  if (any(!is.finite(cap)) || any(cap < specs$lower[[mis_idx]])) {
    stop("The conditional p_mis_base cap falls below its scan lower bound.", call. = FALSE)
  }
  out$p_mis_base <- vapply(seq_len(n), function(i) {
    .ivt_scan_from_unit(
      u[i, mis_idx], specs$lower[[mis_idx]], cap[[i]], "conditional_log"
    )
  }, numeric(1))
  out <- out[, specs$param_symbol, drop = FALSE]
  id_width <- max(5L, nchar(as.character(n)))
  data.frame(
    candidate_id = sprintf(paste0("stage1_", region, "_%0", id_width, "d"), seq_len(n)),
    stage = "stage1",
    region = region,
    parent_candidate_id = NA_character_,
    is_anchor = FALSE,
    conditional_p_mis_base_upper = cap,
    out,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

ivt_scan_extended_design <- function(n_focused,
                                     n_outer,
                                     seed,
                                     envelope,
                                     anchors = NULL) {
  focused <- .ivt_scan_region_design(
    n_focused,
    ivt_scan_extended_parameter_specs("focused"),
    seed,
    envelope,
    "focused"
  )
  outer <- .ivt_scan_region_design(
    n_outer,
    ivt_scan_extended_parameter_specs("outer"),
    as.integer(seed) + 1L,
    envelope,
    "outer"
  )
  rows <- Filter(function(x) is.data.frame(x) && nrow(x), list(focused, outer))
  design <- if (length(rows)) do.call(rbind, rows) else data.frame()

  if (!is.null(anchors) && nrow(anchors)) {
    if (!is.data.frame(anchors) || !"anchor_id" %in% names(anchors)) {
      stop("anchors must be a data frame containing anchor_id.", call. = FALSE)
    }
    anchor_id <- as.character(anchors$anchor_id)
    if (anyNA(anchor_id) || any(!nzchar(anchor_id)) || anyDuplicated(anchor_id)) {
      stop("anchor_id values must be unique non-empty strings.", call. = FALSE)
    }
    missing <- setdiff(.ivt_scan_extended_parameter_names, names(anchors))
    if (length(missing)) {
      stop("anchors is missing parameter(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    clean_id <- gsub("[^A-Za-z0-9]+", "_", anchor_id)
    clean_id <- gsub("^_+|_+$", "", clean_id)
    candidate_id <- paste0("anchor_", clean_id)
    if (anyDuplicated(candidate_id)) {
      stop("Sanitized anchor candidate IDs must be unique.", call. = FALSE)
    }
    anchor_params <- anchors[, .ivt_scan_extended_parameter_names, drop = FALSE]
    for (name in names(anchor_params)) anchor_params[[name]] <- as.numeric(anchor_params[[name]])
    anchor_rows <- data.frame(
      candidate_id = candidate_id,
      stage = "stage1",
      region = "anchor",
      parent_candidate_id = NA_character_,
      is_anchor = TRUE,
      conditional_p_mis_base_upper = .ivt_scan_conditional_cap(anchor_params$p_wgd, envelope),
      anchor_params,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    design <- if (nrow(design)) rbind(design, anchor_rows) else anchor_rows
  }
  if (anyDuplicated(design$candidate_id)) {
    stop("Stage-1 candidate IDs must be unique.", call. = FALSE)
  }
  rownames(design) <- NULL
  design
}

ivt_scan_target_gate <- function(df,
                                 control_limit = 0.05,
                                 wgd_abs_tolerance = 0.15,
                                 two_n_mean_target = 77.45,
                                 two_n_mean_tolerance = 8,
                                 two_n_mass_target = 0.625,
                                 two_n_mass_tolerance = 0.20,
                                 four_n_mean_target = 81.475,
                                 four_n_mean_tolerance = 4,
                                 four_n_decline_lower = 0,
                                 four_n_decline_upper = 6) {
  if (!is.data.frame(df) || !nrow(df)) {
    stop("df must be a non-empty data frame.", call. = FALSE)
  }
  required <- c(
    "status", "protocol_feasibility_status", "n_scenarios",
    "n_summary_rows", "n_growth", "n_death", "n_flow_samples",
    "n_ploidy_samples", "n_insufficient_boundaries", "n_invalid_distributions",
    "n_missing_predictions", "n_negative_predictions",
    "all_passage_boundaries_feasible", "all_passage_executed",
    "all_selected_at_or_after_last_observation", "all_selected_within_tolerance",
    "control_flow_mass_ge3", "control_flow_wgd_mass", "n_low_o2_flow_groups",
    "pooled_2N_A19_mean_N", "pooled_2N_A19_mass_ge70",
    "pooled_4N_A17_mean_N", "decline_4N_to_A17_N",
    paste0("flow_2N_low_A", c(12, 18, 23), "_wgd_mass"),
    paste0("observed_flow_2N_low_A", c(12, 18, 23), "_wgd_mass"),
    paste0("flow_2N_low_A", c(12, 18, 23), "_bimodal_count")
  )
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("df is missing target-gate column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }

  technical <- ivt_scan_hard_gate(df, control_limit = control_limit)
  expected_counts <- c(
    n_summary_rows = 114, n_growth = 219, n_death = 90,
    n_flow_samples = 20, n_ploidy_samples = 12,
    n_invalid_distributions = 0, n_low_o2_flow_groups = 6
  )
  for (name in names(expected_counts)) {
    value <- suppressWarnings(as.numeric(df[[name]]))
    technical <- technical & is.finite(value) & value == expected_counts[[name]]
  }
  control_wgd <- suppressWarnings(as.numeric(df$control_flow_wgd_mass))
  technical <- technical & is.finite(control_wgd) & control_wgd <= control_limit

  bimodal_cols <- paste0("flow_2N_low_A", c(12, 18, 23), "_bimodal_count")
  bimodal_counts <- do.call(cbind, lapply(bimodal_cols, function(name) {
    suppressWarnings(as.numeric(df[[name]]))
  }))
  bimodality_pass <- Reduce(`&`, lapply(bimodal_cols, function(name) {
    value <- suppressWarnings(as.numeric(df[[name]]))
    is.finite(value) & value == 2
  }))
  wgd_errors <- do.call(cbind, lapply(c(12, 18, 23), function(passage) {
    predicted <- suppressWarnings(as.numeric(df[[paste0("flow_2N_low_A", passage, "_wgd_mass")]]))
    observed <- suppressWarnings(as.numeric(df[[paste0("observed_flow_2N_low_A", passage, "_wgd_mass")]]))
    abs(predicted - observed)
  }))
  wgd_mass_pass <- apply(wgd_errors, 1L, function(x) {
    all(is.finite(x)) && all(x <= wgd_abs_tolerance)
  })
  two_n_mean <- suppressWarnings(as.numeric(df$pooled_2N_A19_mean_N))
  two_n_mass <- suppressWarnings(as.numeric(df$pooled_2N_A19_mass_ge70))
  kary_2N_pass <- is.finite(two_n_mean) & is.finite(two_n_mass) &
    abs(two_n_mean - two_n_mean_target) <= two_n_mean_tolerance &
    abs(two_n_mass - two_n_mass_target) <= two_n_mass_tolerance
  four_n_mean <- suppressWarnings(as.numeric(df$pooled_4N_A17_mean_N))
  four_n_decline <- suppressWarnings(as.numeric(df$decline_4N_to_A17_N))
  kary_4N_pass <- is.finite(four_n_mean) & is.finite(four_n_decline) &
    abs(four_n_mean - four_n_mean_target) <= four_n_mean_tolerance &
    four_n_decline >= four_n_decline_lower & four_n_decline <= four_n_decline_upper
  bimodality_gap <- apply(bimodal_counts, 1L, function(x) {
    if (!all(is.finite(x))) return(NA_real_)
    mean(abs(x - 2) / 2)
  })
  wgd_mass_gap <- apply(wgd_errors, 1L, function(x) {
    if (!all(is.finite(x))) return(NA_real_)
    mean(pmax(x - wgd_abs_tolerance, 0) / wgd_abs_tolerance)
  })
  kary_2N_gap <- rowMeans(cbind(
    pmax(abs(two_n_mean - two_n_mean_target) - two_n_mean_tolerance, 0) /
      two_n_mean_tolerance,
    pmax(abs(two_n_mass - two_n_mass_target) - two_n_mass_tolerance, 0) /
      two_n_mass_tolerance
  ))
  four_n_decline_span <- four_n_decline_upper - four_n_decline_lower
  if (four_n_decline_span <= 0) four_n_decline_span <- 1
  kary_4N_gap <- rowMeans(cbind(
    pmax(abs(four_n_mean - four_n_mean_target) - four_n_mean_tolerance, 0) /
      four_n_mean_tolerance,
    pmax(
      four_n_decline_lower - four_n_decline,
      four_n_decline - four_n_decline_upper,
      0
    ) / four_n_decline_span
  ))
  target_gap_score <- rowSums(cbind(
    bimodality_gap, wgd_mass_gap, kary_2N_gap, kary_4N_gap
  ))
  target_gap_score[
    !is.finite(bimodality_gap) |
      !is.finite(wgd_mass_gap) |
      !is.finite(kary_2N_gap) |
      !is.finite(kary_4N_gap)
  ] <- NA_real_
  target_gate_failure_count <- rowSums(!cbind(
    bimodality_pass, wgd_mass_pass, kary_2N_pass, kary_4N_pass
  ))
  phenotype <- bimodality_pass & wgd_mass_pass & kary_2N_pass & kary_4N_pass
  failure_class <- ifelse(
    !technical,
    "technical_fail",
    ifelse(
      phenotype,
      "target_pass",
      ifelse(
        bimodality_pass & wgd_mass_pass & kary_2N_pass & !kary_4N_pass,
        "two_n_close_four_n_fail",
        ifelse(
          kary_4N_pass & !(bimodality_pass & wgd_mass_pass & kary_2N_pass),
          "four_n_close_two_n_fail",
          "multiple_target_fail"
        )
      )
    )
  )
  cbind(
    df,
    data.frame(
      technical_pass = as.logical(technical),
      bimodality_pass = as.logical(bimodality_pass),
      wgd_mass_pass = as.logical(wgd_mass_pass),
      kary_2N_pass = as.logical(kary_2N_pass),
      kary_4N_pass = as.logical(kary_4N_pass),
      bimodality_gap = as.numeric(bimodality_gap),
      wgd_mass_gap = as.numeric(wgd_mass_gap),
      kary_2N_gap = as.numeric(kary_2N_gap),
      kary_4N_gap = as.numeric(kary_4N_gap),
      target_gap_score = as.numeric(target_gap_score),
      target_gate_failure_count = as.integer(target_gate_failure_count),
      phenotype_pass = as.logical(phenotype),
      full_target_pass = as.logical(technical & phenotype),
      target_failure_class = failure_class,
      stringsAsFactors = FALSE
    )
  )
}

ivt_scan_select_refinement_parents <- function(stage1_summary,
                                               specs,
                                               max_parents = 100L) {
  required <- c(
    "candidate_id", "hard_pass", "full_target_pass", "target_failure_class",
    "pareto_front", "target_gap_score", "objective",
    "bimodality_pass", "wgd_mass_pass", "kary_2N_pass", "kary_4N_pass",
    .ivt_scan_extended_parameter_names
  )
  missing <- setdiff(required, names(stage1_summary))
  if (!is.data.frame(stage1_summary) || length(missing)) {
    stop(
      "stage1_summary is missing refinement-parent column(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(as.character(stage1_summary$candidate_id))) {
    stop("stage1_summary candidate_id values must be unique.", call. = FALSE)
  }
  max_parents <- as.integer(max_parents)
  if (length(max_parents) != 1L || is.na(max_parents) || max_parents <= 0L) {
    stop("max_parents must be a positive integer.", call. = FALSE)
  }

  hard <- !is.na(stage1_summary$hard_pass) & as.logical(stage1_summary$hard_pass)
  full <- !is.na(stage1_summary$full_target_pass) & as.logical(stage1_summary$full_target_pass)
  failure_class <- as.character(stage1_summary$target_failure_class)
  pareto <- !is.na(stage1_summary$pareto_front) & as.logical(stage1_summary$pareto_front)
  priority <- rep(NA_integer_, nrow(stage1_summary))
  pool_name <- rep(NA_character_, nrow(stage1_summary))
  priority[hard & pareto] <- 4L
  pool_name[hard & pareto] <- "global_pareto"
  four_n_close <- hard & !is.na(failure_class) & failure_class == "four_n_close_two_n_fail"
  two_n_close <- hard & !is.na(failure_class) & failure_class == "two_n_close_four_n_fail"
  priority[four_n_close] <- 3L
  pool_name[four_n_close] <- "four_n_close_two_n_fail"
  priority[two_n_close] <- 2L
  pool_name[two_n_close] <- "two_n_close_four_n_fail"
  priority[hard & full] <- 1L
  pool_name[hard & full] <- "target_pass"
  priority[hard & !is.finite(priority)] <- 5L
  pool_name[hard & is.na(pool_name)] <- "target_gap_only"
  keep <- is.finite(priority)
  pool <- stage1_summary[keep, , drop = FALSE]
  priority <- priority[keep]
  pool_name <- pool_name[keep]
  if (!nrow(pool)) {
    stop("No refinement-eligible parent candidate was found.", call. = FALSE)
  }

  target_gap <- suppressWarnings(as.numeric(pool$target_gap_score))
  objective <- suppressWarnings(as.numeric(pool$objective))
  if (any(!is.finite(target_gap))) {
    stop("Refinement-eligible candidates must have finite target_gap_score values.", call. = FALSE)
  }
  objective[!is.finite(objective)] <- Inf
  target_gate_failure_count <- rowSums(!cbind(
    as.logical(pool$bimodality_pass),
    as.logical(pool$wgd_mass_pass),
    as.logical(pool$kary_2N_pass),
    as.logical(pool$kary_4N_pass)
  ))
  if (any(!is.finite(target_gate_failure_count))) {
    stop("Refinement-eligible candidates must have finite target-gate flags.", call. = FALSE)
  }

  coordinates <- matrix(NA_real_, nrow = nrow(pool), ncol = nrow(specs))
  for (j in seq_len(nrow(specs))) {
    name <- as.character(specs$param_symbol[[j]])
    transform <- as.character(specs$transform[[j]])
    lower <- as.numeric(specs$lower[[j]])
    upper <- as.numeric(specs$upper[[j]])
    value <- pmin(pmax(as.numeric(pool[[name]]), lower), upper)
    coordinate <- .ivt_scan_to_coordinate(value, transform)
    limits <- .ivt_scan_to_coordinate(c(lower, upper), transform)
    span <- abs(diff(range(limits)))
    coordinates[, j] <- if (span > 0) {
      (coordinate - min(limits)) / span
    } else {
      0
    }
  }

  selected <- integer()
  target_focus_fraction <- 0.10
  n_target_focused <- min(
    nrow(pool),
    max_parents,
    max(1L, as.integer(ceiling(max_parents * target_focus_fraction)))
  )
  target_order <- order(
    target_gate_failure_count,
    target_gap,
    priority,
    objective,
    as.character(pool$candidate_id),
    na.last = TRUE
  )
  selected <- target_order[seq_len(n_target_focused)]
  for (tier in sort(unique(priority))) {
    available <- setdiff(which(priority == tier), selected)
    while (length(available) && length(selected) < max_parents) {
      min_distance <- vapply(available, function(idx) {
        delta <- sweep(
          coordinates[selected, , drop = FALSE],
          2L,
          coordinates[idx, ],
          "-"
        )
        min(sqrt(rowSums(delta^2)))
      }, numeric(1))
      best <- available[min_distance == max(min_distance)]
      chosen <- best[[order(
        target_gap[best],
        objective[best],
        as.character(pool$candidate_id[best])
      )[[1L]]]]
      selected <- c(selected, chosen)
      available <- setdiff(available, chosen)
    }
    if (length(selected) >= max_parents) break
  }
  out <- pool[selected, , drop = FALSE]
  out$refinement_parent_pool <- pool_name[selected]
  rownames(out) <- NULL
  out
}

ivt_scan_stage2_design <- function(stage1_summary,
                                   specs,
                                   n,
                                   seed,
                                   envelope,
                                   max_parents = 100L,
                                   neighborhood_fraction = 0.05) {
  if (!is.data.frame(stage1_summary) || !nrow(stage1_summary)) {
    stop("stage1_summary must be a non-empty data frame.", call. = FALSE)
  }
  required <- c("candidate_id", .ivt_scan_extended_parameter_names)
  missing <- setdiff(required, names(stage1_summary))
  if (length(missing)) {
    stop("stage1_summary is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  max_parents <- as.integer(max_parents)
  n <- as.integer(n)
  if (n <= 0L || max_parents <= 0L) stop("n and max_parents must be positive.", call. = FALSE)
  neighborhood_fraction <- .ivt_scan_assert_scalar_number(
    neighborhood_fraction, "neighborhood_fraction", positive = TRUE
  )
  if (neighborhood_fraction > 1) stop("neighborhood_fraction must be <= 1.", call. = FALSE)

  eligible_envelope <- envelope[
    !is.na(envelope$eligible) & as.logical(envelope$eligible),
    ,
    drop = FALSE
  ]
  if (!nrow(eligible_envelope)) {
    stop("The conditional envelope has no eligible p_wgd range.", call. = FALSE)
  }
  p_wgd_idx <- match("p_wgd", as.character(specs$param_symbol))
  if (is.na(p_wgd_idx)) stop("specs must contain p_wgd.", call. = FALSE)
  specs$lower[[p_wgd_idx]] <- max(
    as.numeric(specs$lower[[p_wgd_idx]]),
    min(as.numeric(eligible_envelope$p_wgd))
  )
  specs$upper[[p_wgd_idx]] <- min(
    as.numeric(specs$upper[[p_wgd_idx]]),
    max(as.numeric(eligible_envelope$p_wgd))
  )
  if (specs$lower[[p_wgd_idx]] > specs$upper[[p_wgd_idx]]) {
    stop("The conditional envelope does not overlap the Stage-2 p_wgd bounds.", call. = FALSE)
  }

  parents <- ivt_scan_select_refinement_parents(
    stage1_summary = stage1_summary,
    specs = specs,
    max_parents = max_parents
  )

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  on.exit(.ivt_scan_restore_rng(had_seed, old_seed, old_kind), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))

  parent_index <- rep(seq_len(nrow(parents)), length.out = n)
  parent_index <- parent_index[sample.int(n, n, replace = FALSE)]
  out <- data.frame(matrix(NA_real_, nrow = n, ncol = nrow(specs)))
  names(out) <- as.character(specs$param_symbol)
  for (j in seq_len(nrow(specs))) {
    name <- as.character(specs$param_symbol[[j]])
    transform <- as.character(specs$transform[[j]])
    lower <- as.numeric(specs$lower[[j]])
    upper <- as.numeric(specs$upper[[j]])
    parent_value <- pmin(pmax(as.numeric(parents[[name]][parent_index]), lower), upper)
    if (name == "p_mis_base") next
    center <- .ivt_scan_to_coordinate(parent_value, transform)
    limits <- .ivt_scan_to_coordinate(c(lower, upper), transform)
    radius <- neighborhood_fraction * abs(diff(range(limits)))
    coordinate <- center + stats::runif(n, -radius, radius)
    coordinate <- pmin(pmax(coordinate, min(limits)), max(limits))
    if (transform %in% c("log", "conditional_log")) {
      out[[name]] <- exp(coordinate)
    } else if (transform == "complement_log") {
      out[[name]] <- 1 - exp(coordinate)
    } else {
      out[[name]] <- coordinate
    }
    out[[name]] <- pmin(pmax(out[[name]], lower), upper)
  }
  cap <- .ivt_scan_conditional_cap(out$p_wgd, envelope)
  mis_idx <- match("p_mis_base", specs$param_symbol)
  cap <- pmin(cap, specs$upper[[mis_idx]])
  parent_mis <- pmin(
    pmax(as.numeric(parents$p_mis_base[parent_index]), specs$lower[[mis_idx]]),
    cap
  )
  lower_log <- log(specs$lower[[mis_idx]])
  cap_log <- log(cap)
  radius <- neighborhood_fraction * (cap_log - lower_log)
  out$p_mis_base <- exp(
    pmin(
      pmax(log(parent_mis) + stats::runif(n, -radius, radius), lower_log),
      cap_log
    )
  )
  out$p_mis_base <- pmin(pmax(out$p_mis_base, specs$lower[[mis_idx]]), cap)
  out <- out[, specs$param_symbol, drop = FALSE]
  local_draw_index <- stats::ave(
    seq_len(n),
    as.character(parents$candidate_id[parent_index]),
    FUN = seq_along
  )
  data.frame(
    candidate_id = sprintf("stage2_%05d", seq_len(n)),
    stage = "stage2",
    region = "local_refinement",
    parent_candidate_id = as.character(parents$candidate_id[parent_index]),
    parent_pool = as.character(parents$refinement_parent_pool[parent_index]),
    local_draw_index = as.integer(local_draw_index),
    local_radius = neighborhood_fraction,
    is_anchor = FALSE,
    conditional_p_mis_base_upper = cap,
    out,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

ivt_scan_validate_checkpoint_membership <- function(checkpoint_objects,
                                                     design,
                                                     batches,
                                                     stage = NULL,
                                                     design_digest = NULL) {
  if (!is.list(checkpoint_objects) || length(checkpoint_objects) != length(batches)) {
    stop("Checkpoint count does not match the requested batches.", call. = FALSE)
  }
  if (!is.data.frame(design) || !all(c("candidate_id", "parent_candidate_id") %in% names(design))) {
    stop("design must contain candidate_id and parent_candidate_id.", call. = FALSE)
  }
  for (i in seq_along(batches)) {
    object <- checkpoint_objects[[i]]
    idx <- batches[[i]]
    expected_ids <- as.character(design$candidate_id[idx])
    expected_parent <- as.character(design$parent_candidate_id[idx])
    fail <- function(reason) {
      stop("Checkpoint batch ", i, " ", reason, ".", call. = FALSE)
    }
    if (!is.list(object) || !identical(as.integer(object$batch_id), as.integer(i))) {
      fail("has an invalid batch id")
    }
    if (!identical(as.integer(object$schema_version), 2L)) {
      fail("has an invalid schema version")
    }
    if (!is.null(stage) && !identical(as.character(object$stage), as.character(stage))) {
      fail("has an invalid stage")
    }
    if (!is.null(design_digest) &&
        !identical(as.character(object$design_digest), as.character(design_digest))) {
      fail("has an invalid design digest")
    }
    if (!identical(as.character(object$candidate_ids), expected_ids) ||
        !is.data.frame(object$summary) ||
        !identical(as.character(object$summary$candidate_id), expected_ids)) {
      fail("candidate membership/order does not match")
    }
    if (!"parent_candidate_id" %in% names(object$summary) ||
        !identical(as.character(object$summary$parent_candidate_id), expected_parent)) {
      fail("parent metadata does not match")
    }
  }
  TRUE
}
