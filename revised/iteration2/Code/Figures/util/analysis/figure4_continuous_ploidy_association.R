#!/usr/bin/env Rscript

# Derive the continuous fitted-parameter/ploidy association layer used by
# Figure 4B and Supplementary Figure 4-1. The 500 rows are optimizer-derived fitted
# endpoints.  They are not posterior samples or biological replicates.

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
    stop("A zero peak rho cannot be assigned to the requested positive/negative display groups.")
  }
  setorder(ranking, -max_abs_rho, parameter_order)
  ranking[, importance_rank := seq_len(.N)]
  ranking[, `:=`(
    peak_direction = fifelse(
      rho_at_max_abs > 0,
      "Positive peak rho",
      "Negative peak rho"
    ),
    peak_direction_order = fifelse(rho_at_max_abs > 0, 1L, 2L)
  )]
  setorder(ranking, peak_direction_order, -max_abs_rho, parameter_order)
  ranking[, within_direction_rank := seq_len(.N), by = peak_direction_order]
  ranking[, display_order := seq_len(.N)]
  association <- merge(
    association,
    ranking[, .(
      parameter, display_order, importance_rank,
      peak_direction, peak_direction_order, within_direction_rank,
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
      peak_direction, peak_direction_order, within_direction_rank,
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
      peak_direction, peak_direction_order, within_direction_rank,
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
    max_abs_rho, rho_at_max_abs, O2_at_max_abs
  )]
  setorder(pooled_summary, display_order)

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

  validation <- data.table(
    metric = c(
      "n_fitted_endpoints", "n_fixed_o2_values", "n_parameters",
      "n_parameter_o2_correlations", "association_metric",
      "outcome_is_continuous", "binary_ploidy_class_used_in_figure4b",
      "parameter_sort_primary", "parameter_sort_secondary",
      "parameter_sort_tie_break", "positive_peak_parameter_count",
      "negative_peak_parameter_count",
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
      "signed peak rho group: positive first, negative second",
      "descending maximum absolute Spearman rho within sign group",
      "configured parameter order", sum(ranking$rho_at_max_abs > 0),
      sum(ranking$rho_at_max_abs < 0),
      "max_abs_rho", "rho_at_max_abs", best_seed,
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
    burden_audit = burden_audit
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  data_arg <- sub("^--data-dir=", "", args[grepl("^--data-dir=", args)])
  if (!length(data_arg)) stop("Usage: --data-dir=PATH")
  derive_figure4_continuous_ploidy_association(data_dir = data_arg[[1L]])
}
