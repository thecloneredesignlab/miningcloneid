#!/usr/bin/env Rscript

# Derive the continuous fitted-parameter/ploidy association layer used by
# Figure 4B and Supplementary Figure 4-1. The 500 rows are optimizer-derived fitted
# endpoints. They are not posterior samples or biological replicates. O2-window
# comparisons therefore quantify stability across fitted endpoints and are not
# interpreted as population-level biological inference.

suppressPackageStartupMessages(library(data.table))

parameter_display_dictionary <- function() {
  data.table(
    parameter = c(
      "lam_max", "p_mis_base", "p_wgd", "p_misseg", "k_o_mis",
      "O2_crit", "n_O", "o2_S0", "kappa_O", "eta_o2",
      "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "k_clear",
      "buffer_smax", "buffer_beta", "buffer_n_exp"
    ),
    parameter_plot_label = c(
      "lam_max | Maximum division rate",
      "p_mis_base | Baseline missegregation probability",
      "p_wgd | Whole-genome-doubling probability",
      "p_misseg | Maximum stress-linked missegregation increment",
      "k_o_mis | Stress-linked missegregation half-saturation",
      "O2_crit | Critical O2 level",
      "n_O | O2-response Hill exponent",
      "o2_S0 | Low-burden effective O2 supply",
      "kappa_O | O2-drop amplitude",
      "eta_o2 | Chromosome-weighted O2-demand exponent",
      "alpha_o2 | High-ploidy growth damping",
      "gamma_growth | Ploidy exponent of growth penalty",
      "mu_hp | Stress-associated death scale",
      "gamma_mu | Ploidy exponent of stress-associated death",
      "k_clear | Dead-biomass clearance rate",
      "buffer_smax | Maximum post-missegregation survival",
      "buffer_beta | Post-missegregation viability-loss strength",
      "buffer_n_exp | Ploidy exponent of post-missegregation survival"
    ),
    parameter_display_label = c(
      "Maximum division rate (lam_max)",
      "Baseline per-chromosome missegregation probability (p_mis_base)",
      "Per-division whole-genome-doubling probability (p_wgd)",
      "Maximum stress-linked per-chromosome missegregation increment (p_misseg)",
      "Stress-linked missegregation half-saturation scale (k_o_mis)",
      "Critical oxygen level (O2_crit)",
      "Oxygen-response Hill exponent (n_O)",
      "Low-burden effective oxygen supply (o2_S0)",
      "Oxygen-drop amplitude (kappa_O)",
      "Chromosome-number-weighted oxygen-demand exponent (eta_o2)",
      "Resource-stress-dependent high-ploidy growth damping (alpha_o2)",
      "Ploidy exponent of the resource-stress growth penalty (gamma_growth)",
      "Stress-associated death scale (mu_hp)",
      "Ploidy exponent of stress-associated death (gamma_mu)",
      "Dead-biomass clearance rate (k_clear)",
      "Maximum post-missegregation survival (buffer_smax)",
      "Post-missegregation viability-loss strength (buffer_beta)",
      "Ploidy exponent of post-missegregation survival (buffer_n_exp)"
    )
  )
}

figure4_o2_windows <- function() {
  data.table(
    o2_window = c("Low O2", "High O2"),
    window_order = 1:2,
    lower_o2 = c(0, 3),
    upper_o2 = c(1, 5)
  )
}

figure4_min_global_peak_abs_rho <- function() 0.3

figure4_o2_association_decision <- function(delta, bh_q, global_peak_abs_rho) {
  if (length(delta) != 1L || length(bh_q) != 1L ||
      length(global_peak_abs_rho) != 1L ||
      any(!is.finite(c(delta, bh_q, global_peak_abs_rho))) ||
      bh_q < 0 || bh_q > 1 || global_peak_abs_rho < 0 ||
      global_peak_abs_rho > 1) {
    stop("O2 association decisions require finite scalar delta, q, and peak |rho|.")
  }
  minimum_peak <- figure4_min_global_peak_abs_rho()
  peak_passes <- global_peak_abs_rho > minimum_peak
  if (!peak_passes) {
    group <- "O2-independent"
    rule <- sprintf("Global peak |rho| <= %g; fails the magnitude requirement", minimum_peak)
  } else if (bh_q >= 0.05) {
    group <- "O2-independent"
    rule <- "No detected Low-High difference at BH q < 0.05"
  } else if (delta < 0) {
    group <- "High O2"
    rule <- sprintf("High [3,5] exceeds Low [0,1] at BH q < 0.05 and global peak |rho| > %g", minimum_peak)
  } else if (delta > 0) {
    group <- "Low O2"
    rule <- sprintf("Low [0,1] exceeds High [3,5] at BH q < 0.05 and global peak |rho| > %g", minimum_peak)
  } else {
    stop("A significant eligible Low-High contrast cannot have zero observed delta.")
  }
  list(group = group, rule = rule, minimum_peak = minimum_peak, peak_passes = peak_passes)
}

normalized_trapezoid_weights <- function(o2_grid, windows) {
  weights <- matrix(
    0,
    nrow = length(o2_grid),
    ncol = nrow(windows),
    dimnames = list(NULL, windows$o2_window)
  )
  for (window_index in seq_len(nrow(windows))) {
    lower <- windows$lower_o2[[window_index]]
    upper <- windows$upper_o2[[window_index]]
    grid_index <- which(
      o2_grid >= lower - 1e-12 & o2_grid <= upper + 1e-12
    )
    grid <- o2_grid[grid_index]
    if (length(grid) < 2L ||
        abs(grid[[1L]] - lower) > 1e-12 ||
        abs(grid[[length(grid)]] - upper) > 1e-12) {
      stop("The fixed-O2 grid does not span an exact O2 window boundary.")
    }
    delta <- diff(grid)
    local_weights <- numeric(length(grid))
    local_weights[[1L]] <- delta[[1L]] / 2
    local_weights[[length(grid)]] <- delta[[length(delta)]] / 2
    if (length(grid) > 2L) {
      local_weights[2:(length(grid) - 1L)] <-
        (delta[-length(delta)] + delta[-1L]) / 2
    }
    weights[grid_index, window_index] <- local_weights / (upper - lower)
  }
  if (any(abs(colSums(weights) - 1) > 1e-12)) {
    stop("Normalized O2-window trapezoid weights do not sum to one.")
  }
  weights
}

derive_o2_window_statistics <- function(
    parameter_matrix,
    outcome_matrix,
    parameter_names,
    o2_grid,
    bootstrap_reps = 5000L,
    bootstrap_seed = 5826L,
    bootstrap_cores = 1L
) {
  if (nrow(parameter_matrix) != 500L ||
      nrow(outcome_matrix) != 500L ||
      ncol(parameter_matrix) != 18L ||
      ncol(outcome_matrix) != 201L) {
    stop("O2-window bootstrap requires the 500 x 18 / 500 x 201 matrices.")
  }
  if (!is.finite(bootstrap_reps) || bootstrap_reps < 1000L) {
    stop("At least 1000 endpoint bootstrap replicates are required.")
  }
  if (!is.finite(bootstrap_cores) || bootstrap_cores < 1L) {
    stop("The O2-window bootstrap core count must be positive.")
  }
  windows <- figure4_o2_windows()
  window_weights <- normalized_trapezoid_weights(o2_grid, windows)
  observed_rho <- suppressWarnings(stats::cor(
    parameter_matrix,
    outcome_matrix,
    method = "spearman"
  ))
  if (!all(dim(observed_rho) == c(18L, 201L)) ||
      any(!is.finite(observed_rho))) {
    stop("Observed parameter/O2 Spearman matrix is incomplete.")
  }
  # Global peak over all 201 O2 values in [0,5], not either window's mean/peak.
  global_peak_abs_rho <- apply(abs(observed_rho), 1L, max)
  observed_scores <- abs(observed_rho) %*% window_weights
  dimnames(observed_scores) <- list(parameter_names, windows$o2_window)

  set.seed(bootstrap_seed)
  bootstrap_index <- matrix(
    sample.int(
      nrow(parameter_matrix),
      nrow(parameter_matrix) * bootstrap_reps,
      replace = TRUE
    ),
    nrow = nrow(parameter_matrix),
    ncol = bootstrap_reps
  )
  score_one_bootstrap <- function(bootstrap_index_number) {
    row_index <- bootstrap_index[, bootstrap_index_number]
    rho <- suppressWarnings(stats::cor(
      parameter_matrix[row_index, , drop = FALSE],
      outcome_matrix[row_index, , drop = FALSE],
      method = "spearman"
    ))
    score <- abs(rho) %*% window_weights
    if (any(!is.finite(score))) {
      stop("A bootstrap replicate produced a non-finite O2-window score.")
    }
    score
  }
  worker <- if (bootstrap_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      seq_len(bootstrap_reps),
      score_one_bootstrap,
      mc.cores = bootstrap_cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(bootstrap_reps), score_one_bootstrap)
  }
  bootstrap_scores <- simplify2array(worker)
  expected_dimensions <- c(18L, nrow(windows), bootstrap_reps)
  if (!identical(dim(bootstrap_scores), expected_dimensions)) {
    stop("Unexpected O2-window bootstrap score dimensions.")
  }

  score_rows <- vector("list", length(parameter_names) * nrow(windows))
  output_index <- 0L
  for (parameter_index in seq_along(parameter_names)) {
    for (window_index in seq_len(nrow(windows))) {
      output_index <- output_index + 1L
      values <- bootstrap_scores[parameter_index, window_index, ]
      interval <- stats::quantile(
        values,
        probs = c(0.025, 0.5, 0.975),
        names = FALSE,
        type = 8
      )
      score_rows[[output_index]] <- data.table(
        parameter = parameter_names[[parameter_index]],
        o2_window = windows$o2_window[[window_index]],
        window_order = windows$window_order[[window_index]],
        lower_o2 = windows$lower_o2[[window_index]],
        upper_o2 = windows$upper_o2[[window_index]],
        observed_mean_abs_rho = observed_scores[
          parameter_index, window_index
        ],
        bootstrap_median = interval[[2L]],
        bootstrap_ci_lower = interval[[1L]],
        bootstrap_ci_upper = interval[[3L]],
        bootstrap_reps = bootstrap_reps,
        bootstrap_seed = bootstrap_seed
      )
    }
  }
  window_scores <- rbindlist(score_rows)

  contrast_definition <- data.table(
    contrast = "Low - High",
    window_a = "Low O2",
    window_b = "High O2",
    window_a_index = 1L,
    window_b_index = 2L
  )
  contrast_rows <- vector(
    "list",
    length(parameter_names) * nrow(contrast_definition)
  )
  output_index <- 0L
  for (parameter_index in seq_along(parameter_names)) {
    for (contrast_index in seq_len(nrow(contrast_definition))) {
      output_index <- output_index + 1L
      definition <- contrast_definition[contrast_index]
      delta_bootstrap <-
        bootstrap_scores[parameter_index, definition$window_a_index, ] -
        bootstrap_scores[parameter_index, definition$window_b_index, ]
      delta_observed <-
        observed_scores[parameter_index, definition$window_a_index] -
        observed_scores[parameter_index, definition$window_b_index]
      interval <- stats::quantile(
        delta_bootstrap,
        probs = c(0.025, 0.5, 0.975),
        names = FALSE,
        type = 8
      )
      lower_tail <-
        (sum(delta_bootstrap <= 0) + 1) / (bootstrap_reps + 1)
      upper_tail <-
        (sum(delta_bootstrap >= 0) + 1) / (bootstrap_reps + 1)
      contrast_rows[[output_index]] <- data.table(
        parameter = parameter_names[[parameter_index]],
        contrast = definition$contrast,
        contrast_order = contrast_index,
        window_a = definition$window_a,
        window_b = definition$window_b,
        observed_delta_mean_abs_rho = delta_observed,
        bootstrap_median_delta = interval[[2L]],
        bootstrap_ci_lower = interval[[1L]],
        bootstrap_ci_upper = interval[[3L]],
        bootstrap_sign_p_value = min(1, 2 * min(lower_tail, upper_tail)),
        bootstrap_reps = bootstrap_reps,
        bootstrap_seed = bootstrap_seed
      )
    }
  }
  pairwise_tests <- rbindlist(contrast_rows)
  pairwise_tests[, bh_adjusted_p_value := stats::p.adjust(
    bootstrap_sign_p_value,
    method = "BH"
  )]
  pairwise_tests[, significant_bh_0p05 := bh_adjusted_p_value < 0.05]

  group_levels <- c("High O2", "Low O2", "O2-independent")
  classification_rows <- vector("list", length(parameter_names))
  for (parameter_index in seq_along(parameter_names)) {
    current_parameter <- parameter_names[[parameter_index]]
    current_tests <- pairwise_tests[parameter == current_parameter]
    if (nrow(current_tests) != 1L) {
      stop("Each parameter must have exactly one Low-High contrast.")
    }
    significant_count <- sum(current_tests$significant_bh_0p05)
    delta <- current_tests$observed_delta_mean_abs_rho[[1L]]
    decision <- figure4_o2_association_decision(
      delta, current_tests$bh_adjusted_p_value[[1L]],
      global_peak_abs_rho[[parameter_index]]
    )
    assigned_group <- decision$group
    decision_rule <- decision$rule
    classification_rows[[parameter_index]] <- data.table(
      parameter = current_parameter,
      o2_association_group = assigned_group,
      o2_association_group_order = match(assigned_group, group_levels),
      o2_window_significant_contrast_count = significant_count,
      o2_association_global_peak_abs_rho = global_peak_abs_rho[[parameter_index]],
      o2_association_min_peak_abs_rho = decision$minimum_peak,
      o2_association_peak_abs_rho_passes = decision$peak_passes,
      o2_window_decision_rule = decision_rule
    )
  }
  classification <- rbindlist(classification_rows)
  if (anyNA(classification$o2_association_group_order)) {
    stop("An O2-window association group is outside the configured order.")
  }
  if (classification[
        o2_window_significant_contrast_count == 0L | !o2_association_peak_abs_rho_passes,
        any(o2_association_group != "O2-independent")
      ] ||
      classification[
        o2_window_significant_contrast_count > 0L & o2_association_peak_abs_rho_passes,
        any(o2_association_group == "O2-independent")
      ]) {
    stop("O2-window significance, global-peak magnitude gate, and groups disagree.")
  }

  list(
    observed_rho = observed_rho,
    window_scores = window_scores,
    pairwise_tests = pairwise_tests,
    classification = classification,
    windows = windows,
    bootstrap_reps = bootstrap_reps,
    bootstrap_seed = bootstrap_seed,
    bootstrap_cores = bootstrap_cores
  )
}

derive_figure4_continuous_ploidy_association <- function(data_dir) {
  data_dir <- normalizePath(data_dir, mustWork = TRUE)

  input_paths <- file.path(data_dir, c(
    "fixed_o2_dominant_ploidy_201grid.tsv",
    "invivo_best_parameters_500seeds.tsv",
    "parameter_function_groups.tsv",
    "parameter_function_group_palette.tsv",
    "invivo_parameter_table_seed25.csv",
    "invivo_best_tsne_cluster_coordinates.tsv",
    "figure4a_seed25_burden_timecourse.tsv",
    "invivo_fit_objective_ranking_500seeds.tsv"
  ))
  missing <- input_paths[!file.exists(input_paths)]
  if (length(missing)) {
    stop(
      "Missing Figure 4 continuous-association input(s): ",
      paste(missing, collapse = ", ")
    )
  }

  fixed <- fread(input_paths[[1L]])
  parameters <- fread(input_paths[[2L]])
  parameter_meta <- fread(input_paths[[3L]])
  parameter_palette <- fread(input_paths[[4L]])
  bounds <- fread(input_paths[[5L]])
  clusters <- fread(input_paths[[6L]])[dataset == "invivo"]
  burden <- fread(input_paths[[7L]])
  objective_ranking <- fread(input_paths[[8L]])

  required_fixed <- c("seed_number", "O2_pct", "dominant_mean_ploidy")
  if (!all(required_fixed %in% names(fixed))) {
    stop("Fixed-O2 table lacks continuous ploidy analysis columns.")
  }
  all_parameters <- parameter_meta[order(parameter_order), parameter]
  if (length(all_parameters) != 18L || !all(all_parameters %in% names(parameters))) {
    stop("Continuous association requires the configured 18 fitted parameters.")
  }
  if (nrow(parameters) != 500L || uniqueN(parameters$seed_number) != 500L ||
      nrow(fixed) != 500L * 201L || uniqueN(fixed$O2_pct) != 201L ||
      uniqueN(fixed$seed_number) != 500L) {
    stop("Expected 500 fitted endpoints evaluated on exactly 201 O2 values.")
  }
  required_objective <- c("seed_number", "objective_rank", "objective")
  if (!all(required_objective %in% names(objective_ranking)) ||
      nrow(objective_ranking) != 500L ||
      uniqueN(objective_ranking$seed_number) != 500L ||
      !identical(sort(objective_ranking$objective_rank), seq_len(500L))) {
    stop("The canonical in-vivo objective ranking is incomplete.")
  }
  best_seed <- objective_ranking[objective_rank == 1L, seed_number]
  if (length(best_seed) != 1L || best_seed != 25L) {
    stop(
      "The canonical lowest-objective in-vivo fit must be seed25; observed: ",
      paste(best_seed, collapse = ",")
    )
  }

  display <- parameter_display_dictionary()
  if (!setequal(display$parameter, all_parameters)) {
    stop("The standardized parameter display dictionary is incomplete.")
  }
  display <- merge(
    parameter_meta,
    display,
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(display, parameter_order)
  fwrite(display, file.path(data_dir, "parameter_display_labels.tsv"), sep = "\t")

  parameter_wide <- parameters[, c("seed_number", all_parameters), with = FALSE]
  parameter_long_natural <- melt(
    parameter_wide,
    id.vars = "seed_number",
    variable.name = "parameter",
    value.name = "natural_value"
  )
  parameter_long_natural[, parameter := as.character(parameter)]

  fixed_min <- fixed[, .(
    seed_number = as.integer(seed_number),
    O2_pct = as.numeric(O2_pct),
    dominant_mean_ploidy = as.numeric(dominant_mean_ploidy)
  )]
  association_input <- merge(
    fixed_min,
    parameter_long_natural,
    by = "seed_number",
    allow.cartesian = TRUE,
    sort = FALSE
  )
  association <- association_input[, {
    keep <- is.finite(natural_value) & is.finite(dominant_mean_ploidy)
    n_complete <- sum(keep)
    rho <- if (n_complete >= 3L &&
               uniqueN(natural_value[keep]) > 1L &&
               uniqueN(dominant_mean_ploidy[keep]) > 1L) {
      suppressWarnings(cor(
        natural_value[keep],
        dominant_mean_ploidy[keep],
        method = "spearman"
      ))
    } else {
      NA_real_
    }
    .(n_complete = n_complete, spearman_rho = rho)
  }, by = .(parameter, O2_pct)]
  association <- merge(
    association,
    display[, .(
      parameter, parameter_group, parameter_order,
      parameter_plot_label, parameter_display_label
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )

  seed_order <- sort(as.integer(parameter_wide$seed_number))
  o2_grid <- sort(unique(fixed_min$O2_pct))
  parameter_matrix <- as.matrix(
    parameter_wide[match(seed_order, seed_number), ..all_parameters]
  )
  fixed_matrix_table <- copy(fixed_min)
  setorder(fixed_matrix_table, seed_number, O2_pct)
  if (!identical(unique(fixed_matrix_table$seed_number), seed_order) ||
      !identical(unique(fixed_matrix_table$O2_pct), o2_grid)) {
    stop("The fixed-O2 table cannot be aligned to the fitted endpoints.")
  }
  outcome_matrix <- matrix(
    fixed_matrix_table$dominant_mean_ploidy,
    nrow = length(seed_order),
    ncol = length(o2_grid),
    byrow = TRUE
  )
  bootstrap_reps <- suppressWarnings(as.integer(Sys.getenv(
    "FIGURE4_O2_WINDOW_BOOTSTRAP_REPS",
    unset = "5000"
  )))
  bootstrap_seed <- suppressWarnings(as.integer(Sys.getenv(
    "FIGURE4_O2_WINDOW_BOOTSTRAP_SEED",
    unset = "5826"
  )))
  bootstrap_cores <- suppressWarnings(as.integer(Sys.getenv(
    "FIGURE4_O2_WINDOW_BOOTSTRAP_CORES",
    unset = "1"
  )))
  o2_window_statistics <- derive_o2_window_statistics(
    parameter_matrix = parameter_matrix,
    outcome_matrix = outcome_matrix,
    parameter_names = all_parameters,
    o2_grid = o2_grid,
    bootstrap_reps = bootstrap_reps,
    bootstrap_seed = bootstrap_seed,
    bootstrap_cores = bootstrap_cores
  )
  observed_association_check <- association[
    order(match(parameter, all_parameters), O2_pct),
    spearman_rho
  ]
  if (max(abs(
    observed_association_check -
      as.vector(t(o2_window_statistics$observed_rho))
  )) > 1e-12) {
    stop("Bootstrap input does not reproduce the canonical Spearman matrix.")
  }

  ranking <- association[, {
    finite_rows <- which(is.finite(spearman_rho))
    if (!length(finite_rows)) {
      stop("A parameter has no finite fixed-O2 Spearman correlations.")
    }
    finite_rho <- spearman_rho[finite_rows]
    finite_o2 <- O2_pct[finite_rows]
    max_abs_rho <- max(abs(finite_rho))
    peak_candidates <- which(abs(abs(finite_rho) - max_abs_rho) <= 1e-12)
    peak_index <- peak_candidates[which.min(finite_o2[peak_candidates])]
    .(
      n_o2 = length(finite_rho),
      mean_rho = mean(finite_rho),
      median_rho = median(finite_rho),
      rms_rho = sqrt(mean(finite_rho^2)),
      fraction_positive = mean(finite_rho > 0),
      fraction_negative = mean(finite_rho < 0),
      max_abs_rho = max_abs_rho,
      rho_at_max_abs = finite_rho[[peak_index]],
      O2_at_max_abs = finite_o2[[peak_index]]
    )
  }, by = .(
    parameter, parameter_group, parameter_order,
    parameter_plot_label, parameter_display_label
  )]
  if (any(ranking$n_o2 != 201L) ||
      any(!is.finite(ranking$max_abs_rho)) ||
      any(abs(ranking$max_abs_rho - abs(ranking$rho_at_max_abs)) > 1e-12)) {
    stop("The max-|rho| parameter ranking is incomplete or internally inconsistent.")
  }
  if (any(abs(ranking$rho_at_max_abs) <= 1e-12)) {
    stop("A zero peak rho cannot be assigned to a positive/negative point color.")
  }
  setorder(ranking, -max_abs_rho, parameter_order)
  ranking[, importance_rank := seq_len(.N)]
  ranking[, peak_direction := fifelse(
    rho_at_max_abs > 0,
    "Positive peak rho",
    "Negative peak rho"
  )]
  ranking <- merge(
    ranking,
    o2_window_statistics$classification,
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(ranking$o2_association_group) ||
      anyNA(ranking$o2_association_group_order)) {
    stop("O2-window association classifications are incomplete.")
  }
  setorder(
    ranking,
    o2_association_group_order, -max_abs_rho, parameter_order
  )
  ranking[, within_o2_association_group_rank := seq_len(.N),
          by = o2_association_group_order]
  ranking[, display_order := seq_len(.N)]
  association <- merge(
    association,
    ranking[, .(
      parameter, display_order, importance_rank,
      peak_direction,
      o2_association_group, o2_association_group_order,
      within_o2_association_group_rank,
      o2_window_significant_contrast_count, o2_window_decision_rule,
      o2_association_global_peak_abs_rho, o2_association_min_peak_abs_rho,
      o2_association_peak_abs_rho_passes,
      max_abs_rho, rho_at_max_abs, O2_at_max_abs
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(association, display_order, O2_pct)

  bound_meta <- bounds[
    estimate == TRUE & param_prototype %in% all_parameters,
    .(
      parameter = param_prototype,
      initial_value = as.numeric(prototype_init_value),
      lower_bound = as.numeric(prototype_lower_bound),
      upper_bound = as.numeric(prototype_upper_bound),
      transformation = fifelse(grepl("^log10_", param_name), "log10", "identity")
    )
  ]
  if (nrow(bound_meta) != 18L || !setequal(bound_meta$parameter, all_parameters)) {
    stop("Could not recover optimizer bounds for all 18 fitted parameters.")
  }
  if (any(!is.finite(unlist(bound_meta[, .(
        initial_value, lower_bound, upper_bound
      )]))) ||
      any(bound_meta$lower_bound >= bound_meta$upper_bound) ||
      any(bound_meta$initial_value < bound_meta$lower_bound) ||
      any(bound_meta$initial_value > bound_meta$upper_bound)) {
    stop("Configured in-vivo parameter ranges or initial values are invalid.")
  }
  prior_sd <- sqrt(1 / 12)
  parameter_prior <- merge(
    parameter_long_natural,
    bound_meta,
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  parameter_prior[, transformed_value := fifelse(
    transformation == "log10",
    log10(natural_value),
    natural_value
  )]
  parameter_prior[, lower_transformed := fifelse(
    transformation == "log10",
    log10(lower_bound),
    lower_bound
  )]
  parameter_prior[, upper_transformed := fifelse(
    transformation == "log10",
    log10(upper_bound),
    upper_bound
  )]
  parameter_prior[, prior_unit :=
    (transformed_value - lower_transformed) /
      (upper_transformed - lower_transformed)]
  parameter_prior[, prior_referenced_z := (prior_unit - 0.5) / prior_sd]
  parameter_prior <- merge(
    parameter_prior,
    clusters[, .(
      seed_number = as.integer(seed), cluster_id, tSNE1, tSNE2
    )],
    by = "seed_number",
    all.x = TRUE,
    sort = FALSE
  )
  parameter_prior <- merge(
    parameter_prior,
    display[, .(
      parameter, parameter_group, parameter_order,
      parameter_plot_label, parameter_display_label
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  parameter_prior <- merge(
    parameter_prior,
    ranking[, .(
      parameter, display_order, importance_rank,
      peak_direction,
      o2_association_group, o2_association_group_order,
      within_o2_association_group_rank,
      max_abs_rho, rho_at_max_abs, O2_at_max_abs
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  parameter_prior[, is_lowest_objective_fit := seed_number == best_seed]
  setorder(parameter_prior, display_order, seed_number)
  if (nrow(parameter_prior) != 18L * 500L ||
      any(!is.finite(parameter_prior$prior_referenced_z)) ||
      anyNA(parameter_prior$cluster_id) ||
      parameter_prior[, sum(is_lowest_objective_fit)] != 18L) {
    stop("The all-parameter prior-referenced endpoint table is incomplete.")
  }

  positive_log_values <- c(
    parameter_prior$natural_value[parameter_prior$natural_value > 0],
    bound_meta$initial_value[bound_meta$initial_value > 0],
    bound_meta$lower_bound[bound_meta$lower_bound > 0],
    bound_meta$upper_bound[bound_meta$upper_bound > 0]
  )
  if (!length(positive_log_values) || any(!is.finite(positive_log_values))) {
    stop("No finite positive values are available for the Figure 4B log10 axis.")
  }
  log_floor_raw <- 10^(floor(log10(min(positive_log_values))) - 1)
  log_floor_plot <- log10(log_floor_raw)
  to_log10_plot <- function(x) {
    fifelse(x > 0, log10(x), log_floor_plot)
  }
  parameter_prior[, log10_plot_value := to_log10_plot(natural_value)]

  endpoint_ranges <- merge(
    bound_meta,
    ranking[, .(
      parameter, display_order, importance_rank,
      peak_direction,
      o2_association_group, o2_association_group_order,
      within_o2_association_group_rank,
      max_abs_rho, rho_at_max_abs, O2_at_max_abs
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  endpoint_ranges[, `:=`(
    parameter_y = 19 - display_order,
    lower_plot = to_log10_plot(lower_bound),
    upper_plot = to_log10_plot(upper_bound),
    initial_plot = to_log10_plot(initial_value),
    log_floor_raw = log_floor_raw
  )]
  best_log_rows <- parameter_prior[
    is_lowest_objective_fit == TRUE,
    .(
      parameter,
      best_seed_number = seed_number,
      best_value = natural_value,
      best_plot = log10_plot_value
    )
  ]
  endpoint_ranges <- merge(
    endpoint_ranges,
    best_log_rows,
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  endpoint_ranges[, `:=`(
    range_ymin = parameter_y - 0.34,
    range_ymax = parameter_y - 0.08,
    range_ycenter = parameter_y - 0.21,
    density_ybaseline = parameter_y + 0.06
  )]
  setorder(endpoint_ranges, display_order)

  endpoint_density_rows <- vector("list", nrow(endpoint_ranges))
  best_point_y <- numeric(nrow(endpoint_ranges))
  for (index in seq_len(nrow(endpoint_ranges))) {
    current_parameter <- endpoint_ranges$parameter[[index]]
    values <- parameter_prior[
      parameter == current_parameter,
      log10_plot_value
    ]
    span <- diff(range(values))
    if (!is.finite(span) || span <= 1e-10) {
      half_width <- max(
        0.015,
        0.0125 * (
          endpoint_ranges$upper_plot[[index]] -
            endpoint_ranges$lower_plot[[index]]
        )
      )
      grid_x <- values[[1L]] + c(-half_width, 0, half_width)
      density_scaled <- c(0, 0.38, 0)
      density_method <- "point_mass"
    } else {
      estimate <- density(
        values,
        n = 256L,
        from = min(values),
        to = max(values),
        adjust = 0.75,
        cut = 0,
        na.rm = TRUE
      )
      grid_x <- estimate$x
      density_scaled <- 0.38 * estimate$y / max(estimate$y)
      density_method <- "kernel_density_log10"
    }
    baseline <- endpoint_ranges$density_ybaseline[[index]]
    endpoint_density_rows[[index]] <- data.table(
      parameter = current_parameter,
      display_order = endpoint_ranges$display_order[[index]],
      parameter_y = endpoint_ranges$parameter_y[[index]],
      x_plot = grid_x,
      density_scaled = density_scaled,
      y_baseline = baseline,
      y_density = baseline + density_scaled,
      density_method = density_method
    )
    height_at_best <- approx(
      grid_x,
      density_scaled,
      xout = endpoint_ranges$best_plot[[index]],
      rule = 2
    )$y
    best_point_y[[index]] <- baseline + min(
      0.22,
      max(0.05, 0.55 * height_at_best)
    )
  }
  endpoint_density <- rbindlist(endpoint_density_rows)
  endpoint_ranges[, best_point_y := best_point_y]
  if (nrow(endpoint_density) < 18L * 3L ||
      any(!is.finite(endpoint_density$x_plot)) ||
      any(!is.finite(endpoint_density$y_density)) ||
      nrow(best_log_rows) != 18L ||
      any(endpoint_ranges$best_seed_number != best_seed)) {
    stop("The Figure 4B log10 range/density summaries are incomplete.")
  }

  pooled_summary <- parameter_prior[, {
    q <- quantile(prior_referenced_z, c(0.25, 0.50, 0.75), names = FALSE)
    iqr <- q[[3L]] - q[[1L]]
    lower_fence <- q[[1L]] - 1.5 * iqr
    upper_fence <- q[[3L]] + 1.5 * iqr
    lower_whisker <- min(prior_referenced_z[prior_referenced_z >= lower_fence])
    upper_whisker <- max(prior_referenced_z[prior_referenced_z <= upper_fence])
    best_rows <- which(is_lowest_objective_fit)
    if (length(best_rows) != 1L) {
      stop("Each parameter must contain exactly one lowest-objective marker.")
    }
    .(
      n = .N,
      q25_prior_z = q[[1L]],
      median_prior_z = q[[2L]],
      q75_prior_z = q[[3L]],
      lower_whisker_prior_z = lower_whisker,
      upper_whisker_prior_z = upper_whisker,
      best_seed_number = seed_number[[best_rows]],
      best_natural_value = natural_value[[best_rows]],
      best_prior_referenced_z = prior_referenced_z[[best_rows]]
    )
  }, by = .(
    parameter, parameter_group, parameter_plot_label,
    parameter_display_label, display_order, importance_rank,
    o2_association_group, o2_association_group_order,
    within_o2_association_group_rank,
    peak_direction,
    max_abs_rho, rho_at_max_abs, O2_at_max_abs
  )]
  setorder(pooled_summary, display_order)

  o2_window_scores <- merge(
    o2_window_statistics$window_scores,
    ranking[, .(
      parameter, parameter_order, parameter_plot_label,
      parameter_display_label, display_order, importance_rank,
      o2_association_group, o2_association_group_order,
      within_o2_association_group_rank,
      max_abs_rho, rho_at_max_abs, O2_at_max_abs
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(o2_window_scores, display_order, window_order)
  o2_window_pairwise_tests <- merge(
    o2_window_statistics$pairwise_tests,
    ranking[, .(
      parameter, parameter_order, parameter_plot_label,
      parameter_display_label, display_order, importance_rank,
      o2_association_group, o2_association_group_order,
      within_o2_association_group_rank,
      max_abs_rho, rho_at_max_abs, O2_at_max_abs
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(o2_window_pairwise_tests, display_order, contrast_order)
  o2_window_classification <- merge(
    o2_window_statistics$classification,
    ranking[, .(
      parameter, parameter_order, parameter_plot_label,
      parameter_display_label, display_order, importance_rank,
      within_o2_association_group_rank,
      max_abs_rho, rho_at_max_abs, O2_at_max_abs
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(o2_window_classification, display_order)
  if (nrow(o2_window_scores) != 18L * 2L ||
      nrow(o2_window_pairwise_tests) != 18L ||
      nrow(o2_window_classification) != 18L ||
      anyNA(o2_window_scores$display_order) ||
      anyNA(o2_window_pairwise_tests$display_order) ||
      anyNA(o2_window_classification$display_order)) {
    stop("The O2-window statistical audit tables are incomplete.")
  }

  fwrite(
    association,
    file.path(data_dir, "continuous_ploidy_spearman_by_o2.tsv"),
    sep = "\t"
  )
  fwrite(
    ranking,
    file.path(data_dir, "continuous_ploidy_parameter_ranking.tsv"),
    sep = "\t"
  )
  fwrite(
    parameter_prior,
    file.path(data_dir, "all_parameter_fitted_endpoint_values.tsv"),
    sep = "\t"
  )
  fwrite(
    pooled_summary,
    file.path(data_dir, "all_parameter_pooled_distribution_summary.tsv"),
    sep = "\t"
  )
  fwrite(
    endpoint_ranges,
    file.path(data_dir, "all_parameter_log10_range_summary.tsv"),
    sep = "\t"
  )
  fwrite(
    endpoint_density,
    file.path(data_dir, "all_parameter_log10_density.tsv"),
    sep = "\t"
  )
  fwrite(
    o2_window_scores,
    file.path(data_dir, "continuous_ploidy_o2_window_absrho_scores.tsv"),
    sep = "\t"
  )
  fwrite(
    o2_window_pairwise_tests,
    file.path(data_dir, "continuous_ploidy_o2_window_pairwise_tests.tsv"),
    sep = "\t"
  )
  fwrite(
    o2_window_classification,
    file.path(data_dir, "continuous_ploidy_o2_window_classification.tsv"),
    sep = "\t"
  )

  required_burden_columns <- c("harvest", "cohort", "day", "obs_burden")
  if (!all(required_burden_columns %in% names(burden))) {
    stop("Figure 4A burden table lacks measurement-audit columns.")
  }
  burden_obs <- burden[
    is.finite(obs_burden),
    .(
      harvest, cohort, day = as.numeric(day),
      figure4_obs_burden_mm3 = as.numeric(obs_burden)
    )
  ]
  figure1_path <- file.path(dirname(data_dir), "Figure1", "invivo_burden_long.tsv")
  if (!file.exists(figure1_path)) {
    stop("Missing package-internal Figure 1 burden table for Figure 4A audit: ", figure1_path)
  }
  figure1 <- fread(figure1_path)[
    burden_present == TRUE,
    .(
      harvest, cohort, day = as.numeric(day),
      figure1_obs_burden_mm3 = as.numeric(burden)
    )
  ]
  burden_audit <- merge(
    burden_obs,
    figure1,
    by = c("harvest", "cohort", "day"),
    all.x = TRUE,
    sort = FALSE
  )
  burden_audit[, exact_package_internal_match :=
    is.finite(figure1_obs_burden_mm3) &
      abs(figure4_obs_burden_mm3 - figure1_obs_burden_mm3) < 1e-10]
  burden_audit[, `:=`(
    cohort_day_n = .N,
    cohort_day_mean_mm3 = mean(figure4_obs_burden_mm3),
    cohort_day_sd_mm3 = if (.N > 1L) sd(figure4_obs_burden_mm3) else NA_real_,
    cohort_day_min_mm3 = min(figure4_obs_burden_mm3),
    cohort_day_max_mm3 = max(figure4_obs_burden_mm3)
  ), by = .(cohort, day)]
  burden_audit[, observed_line_segment := fifelse(
    cohort == "4N" & day >= 81,
    "4N late segment (n=2)",
    fifelse(
      cohort == "4N",
      "4N early segment (n=4)",
      "2N segment (n=4)"
    )
  )]
  burden_audit[, audit_interpretation := fifelse(
    cohort == "4N" & day == 77,
    paste0(
      "High variance is driven by four verified measurements, including ",
      "5899.278 mm3; retain the point and sample-SD ribbon."
    ),
    fifelse(
      cohort == "4N" & day == 81,
      paste0(
        "Only two tumors remain; break the observed-mean line after day 77 ",
        "because the apparent drop is not a within-tumor trajectory."
      ),
      "Verified package-internal measurement used in the cohort summary."
    )
  )]
  if (any(!burden_audit$exact_package_internal_match)) {
    stop("Figure 4A and Figure 1 packaged burden measurements disagree.")
  }
  setorder(burden_audit, cohort, day, harvest)
  fwrite(
    burden_audit,
    file.path(data_dir, "figure4a_burden_measurement_audit.tsv"),
    sep = "\t"
  )

  o2_group_counts <- ranking[
    order(o2_association_group_order),
    .N,
    by = .(o2_association_group, o2_association_group_order)
  ][, paste0(o2_association_group, "=", N, collapse = ";")]
  validation <- data.table(
    metric = c(
      "n_fitted_endpoints", "n_fixed_o2_values", "n_parameters",
      "n_parameter_o2_correlations", "association_metric",
      "outcome_is_continuous", "binary_ploidy_class_used_in_figure4b",
      "o2_window_score", "o2_window_boundaries",
      "o2_window_bootstrap_unit", "o2_window_bootstrap_reps",
      "o2_window_bootstrap_seed", "o2_window_bootstrap_cores",
      "o2_window_pairwise_test_count", "o2_window_multiple_testing",
      "o2_window_significance_threshold", "o2_association_global_peak_minimum",
      "o2_association_magnitude_gate", "o2_association_group_counts",
      "parameter_sort_primary", "parameter_sort_secondary",
      "parameter_sort_tertiary", "parameter_sort_tie_break",
      "ranking_magnitude_field", "ranking_signed_color_field",
      "canonical_lowest_objective_seed", "pooled_endpoint_n_per_parameter",
      "lowest_objective_markers", "endpoint_display_scale",
      "endpoint_zero_floor_raw", "endpoint_samples_are_posterior",
      "figure4a_observations_audited", "figure4a_package_internal_exact_matches",
      "figure4a_4N_day77_n", "figure4a_4N_day81_n",
      "figure4a_observed_mean_line_break_after_day77"
    ),
    value = c(
      500, 201, 18, nrow(association), "Spearman rho",
      "TRUE", "FALSE",
      "normalized trapezoid AUC of absolute Spearman rho",
      "Low [0,1]; High [3,5]; (1,3) excluded from window statistics",
      "complete fitted-endpoint row with its full 201-point O2 curve",
      o2_window_statistics$bootstrap_reps,
      o2_window_statistics$bootstrap_seed,
      o2_window_statistics$bootstrap_cores,
      nrow(o2_window_pairwise_tests),
      "Benjamini-Hochberg across all 18 Low-High window contrasts",
      "BH q < 0.05", figure4_min_global_peak_abs_rho(),
      "global max |rho| over all 201 O2 values in [0,5] must be strictly > 0.3",
      o2_group_counts,
      paste(
        "O2 association group: High; Low; O2-independent"
      ),
      "descending maximum absolute Spearman rho within O2 association group",
      "configured parameter order for exact max-|rho| ties",
      "none", "max_abs_rho", "rho_at_max_abs", best_seed,
      paste(sort(unique(pooled_summary$n)), collapse = ","),
      sum(parameter_prior$is_lowest_objective_fit),
      "original parameter value on shared log10 axis",
      format(log_floor_raw, scientific = TRUE),
      "FALSE", nrow(burden_audit), sum(burden_audit$exact_package_internal_match),
      burden_audit[cohort == "4N" & day == 77, unique(cohort_day_n)],
      burden_audit[cohort == "4N" & day == 81, unique(cohort_day_n)],
      "TRUE"
    )
  )
  fwrite(
    validation,
    file.path(data_dir, "continuous_ploidy_analysis_validation.tsv"),
    sep = "\t"
  )

  provenance <- data.table(
    source = c(
      basename(input_paths),
      basename(figure1_path)
    ),
    path = c(input_paths, figure1_path),
    md5 = unname(tools::md5sum(c(input_paths, figure1_path))),
    role = c(
      "continuous fixed-O2 outcome", "500 fitted endpoints",
      "parameter function metadata", "parameter function palette",
      "optimizer parameter bounds", "exploratory t-SNE assignments",
      "Figure 4A selected burden measurements",
      "canonical 500-seed objective ranking",
      "package-internal burden cross-check"
    )
  )
  fwrite(
    provenance,
    file.path(data_dir, "continuous_ploidy_analysis_source_provenance.tsv"),
    sep = "\t"
  )

  invisible(list(
    association = association,
    ranking = ranking,
    pooled_summary = pooled_summary,
    o2_window_scores = o2_window_scores,
    o2_window_pairwise_tests = o2_window_pairwise_tests,
    o2_window_classification = o2_window_classification,
    burden_audit = burden_audit
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  data_arg <- sub("^--data-dir=", "", args[grepl("^--data-dir=", args)])
  if (!length(data_arg)) stop("Usage: --data-dir=PATH")
  derive_figure4_continuous_ploidy_association(data_dir = data_arg[[1L]])
}
