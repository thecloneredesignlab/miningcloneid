#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))
source(resolve_utility_file("analysis/optimizer_diagnostics.R"))

weighted_quantile <- function(values, weights, probability) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[keep]
  weights <- weights[keep]
  if (!length(values)) return(NA_real_)
  order_index <- order(values)
  values <- values[order_index]
  weights <- weights[order_index]
  cumulative <- cumsum(weights) / sum(weights)
  values[[which(cumulative >= probability)[[1L]]]]
}

figure3_parameter_order <- c(
  "alpha_o2",
  "buffer_smax",
  "gamma_growth",
  "gamma_mu",
  "init_mean_2N",
  "init_mean_4N",
  "init_sd_2N",
  "init_sd_4N",
  "k_o_mis",
  "mu_hp",
  "O2_crit",
  "p_wgd",
  "n_O",
  "p_mis_base",
  "sigma_kary",
  "lam_max",
  "buffer_n_exp",
  "buffer_beta",
  "sigma_growth",
  "p_misseg"
)

figure3_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

figure3_flag_true <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

figure3_read_best_parameters <- function(path, parameters, seed) {
  table <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  figure3_require_columns(
    table,
    c("parameter", "value"),
    paste0("Figure 3G ", seed, " best-parameter table")
  )
  table$parameter <- trimws(as.character(table$parameter))
  if (anyDuplicated(table$parameter)) {
    stop("Figure 3G found duplicate parameters in: ", path)
  }
  missing <- setdiff(parameters, table$parameter)
  if (length(missing)) {
    stop(
      "Figure 3G ", seed, " lacks fitted parameter(s): ",
      paste(missing, collapse = ", ")
    )
  }
  values <- suppressWarnings(as.numeric(
    table$value[match(parameters, table$parameter)]
  ))
  if (any(!is.finite(values))) {
    stop("Figure 3G found non-finite fitted values in: ", path)
  }
  values
}

build_figure3g_parameter_ensemble <- function(run_root, metrics, output_dir) {
  expected_seeds <- paste0("seed", seq_len(500L))
  seed_numbers <- suppressWarnings(as.integer(sub("^seed", "", metrics$seed)))
  metrics_by_seed <- metrics[order(seed_numbers), , drop = FALSE]
  if (!identical(metrics_by_seed$seed, expected_seeds)) {
    stop("Figure 3G requires exactly seed1 through seed500 once each.")
  }
  best_row <- metrics[which.min(metrics$objective), , drop = FALSE]
  if (nrow(best_row) != 1L ||
      best_row$seed[[1L]] != INVITRO_VISUALIZATION_SEED ||
      best_row$objective_rank[[1L]] != 1L) {
    stop(
      "Figure 3G requires ", INVITRO_VISUALIZATION_SEED,
      " to be the unique objective-rank-1 fit."
    )
  }

  parameter_path <- assert_allowed_input(file.path(
    run_root, INVITRO_VISUALIZATION_SEED, "parameter_table_input.csv"
  ))
  parameter_table <- utils::read.csv(
    parameter_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  figure3_require_columns(
    parameter_table,
    c("param_symbol", "init_value", "lower_bound", "upper_bound"),
    "Figure 3G parameter-range table"
  )
  fitted_column <- if ("use_invitro_fit" %in% names(parameter_table)) {
    "use_invitro_fit"
  } else if ("estimate" %in% names(parameter_table)) {
    "estimate"
  } else {
    stop("Figure 3G parameter table has no fitted-parameter flag.")
  }
  parameter_table$param_symbol <- trimws(as.character(
    parameter_table$param_symbol
  ))
  fitted <- figure3_flag_true(parameter_table[[fitted_column]])
  fitted_table <- parameter_table[fitted, , drop = FALSE]
  if (!setequal(fitted_table$param_symbol, figure3_parameter_order) ||
      nrow(fitted_table) != length(figure3_parameter_order)) {
    stop(
      "Figure 3G fitted-parameter contract changed; expected exactly: ",
      paste(figure3_parameter_order, collapse = ", ")
    )
  }
  fitted_table <- fitted_table[
    match(figure3_parameter_order, fitted_table$param_symbol),
    ,
    drop = FALSE
  ]
  fitted_table$init_value <- suppressWarnings(as.numeric(
    fitted_table$init_value
  ))
  fitted_table$lower_bound <- suppressWarnings(as.numeric(
    fitted_table$lower_bound
  ))
  fitted_table$upper_bound <- suppressWarnings(as.numeric(
    fitted_table$upper_bound
  ))
  if (any(!is.finite(unlist(fitted_table[
    c("init_value", "lower_bound", "upper_bound")
  ]))) ||
      any(fitted_table$lower_bound >= fitted_table$upper_bound) ||
      any(fitted_table$init_value < fitted_table$lower_bound) ||
      any(fitted_table$init_value > fitted_table$upper_bound)) {
    stop("Figure 3G found an invalid configured parameter range or initial value.")
  }

  best_paths <- file.path(run_root, expected_seeds, "best_params.tsv")
  require_files(best_paths, "Figure 3G 500-seed best-parameter input")
  best_paths <- vapply(best_paths, assert_allowed_input, character(1L))
  endpoint_rows <- lapply(seq_along(expected_seeds), function(index) {
    seed <- expected_seeds[[index]]
    metric_row <- metrics_by_seed[index, , drop = FALSE]
    values <- figure3_read_best_parameters(
      best_paths[[index]], figure3_parameter_order, seed
    )
    data.frame(
      seed = seed,
      seed_number = index,
      objective = metric_row$objective[[1L]],
      objective_rank = metric_row$objective_rank[[1L]],
      parameter = figure3_parameter_order,
      parameter_order = seq_along(figure3_parameter_order),
      value = values,
      is_best_seed = seed == best_row$seed[[1L]],
      source_file = best_paths[[index]],
      stringsAsFactors = FALSE
    )
  })
  endpoints <- do.call(rbind, endpoint_rows)
  row.names(endpoints) <- NULL
  if (nrow(endpoints) != 10000L ||
      length(unique(endpoints$seed)) != 500L ||
      any(table(endpoints$parameter) != 500L)) {
    stop("Figure 3G endpoint table does not satisfy the 500 x 20 contract.")
  }

  ranges <- data.frame(
    parameter = fitted_table$param_symbol,
    parameter_order = seq_along(figure3_parameter_order),
    row_y = rev(seq_along(figure3_parameter_order)),
    init_value = fitted_table$init_value,
    lower_bound = fitted_table$lower_bound,
    upper_bound = fitted_table$upper_bound,
    fitted_flag_column = fitted_column,
    stringsAsFactors = FALSE
  )
  bound_index <- match(endpoints$parameter, ranges$parameter)
  bound_width <- ranges$upper_bound[bound_index] - ranges$lower_bound[bound_index]
  bound_tolerance <- pmax(1e-12, 1e-8 * bound_width)
  outside_bound <-
    endpoints$value < ranges$lower_bound[bound_index] - bound_tolerance |
    endpoints$value > ranges$upper_bound[bound_index] + bound_tolerance
  if (any(outside_bound)) {
    stop(
      "Figure 3G found ", sum(outside_bound),
      " fitted endpoints outside configured bounds."
    )
  }

  positive_values <- c(
    ranges$lower_bound[ranges$lower_bound > 0],
    ranges$upper_bound[ranges$upper_bound > 0],
    ranges$init_value[ranges$init_value > 0],
    endpoints$value[endpoints$value > 0]
  )
  if (!length(positive_values)) {
    stop("Figure 3G has no positive values for the requested log10 scale.")
  }
  log_floor_raw <- 10^(floor(log10(min(positive_values))) - 1)
  log_floor_plot <- log10(log_floor_raw)
  to_plot <- function(x) {
    ifelse(x > 0, log10(x), log_floor_plot)
  }
  ranges$lower_plot <- to_plot(ranges$lower_bound)
  ranges$upper_plot <- to_plot(ranges$upper_bound)
  ranges$init_plot <- to_plot(ranges$init_value)
  endpoints$plot_value <- to_plot(endpoints$value)
  ranges$log_floor_raw <- log_floor_raw
  ranges$range_ymin <- ranges$row_y - 0.34
  ranges$range_ymax <- ranges$row_y - 0.08
  ranges$range_ycenter <- ranges$row_y - 0.21
  ranges$density_ybaseline <- ranges$row_y + 0.06

  best_endpoints <- endpoints[endpoints$is_best_seed, , drop = FALSE]
  best_endpoints <- best_endpoints[
    match(ranges$parameter, best_endpoints$parameter),
    ,
    drop = FALSE
  ]
  if (nrow(best_endpoints) != nrow(ranges) || anyNA(best_endpoints$parameter)) {
    stop(
      "Figure 3G could not resolve one ", INVITRO_VISUALIZATION_SEED,
      " value per fitted parameter."
    )
  }
  ranges$best_seed <- best_row$seed[[1L]]
  ranges$best_seed_value <- best_endpoints$value
  ranges$best_seed_plot <- best_endpoints$plot_value

  density_rows <- vector("list", nrow(ranges))
  best_point_y <- numeric(nrow(ranges))
  for (index in seq_len(nrow(ranges))) {
    parameter <- ranges$parameter[[index]]
    values <- endpoints$plot_value[endpoints$parameter == parameter]
    span <- diff(range(values))
    if (!is.finite(span) || span <= 1e-10) {
      half_width <- max(
        0.015,
        0.0125 * (ranges$upper_plot[[index]] - ranges$lower_plot[[index]])
      )
      grid_x <- values[[1L]] + c(-half_width, 0, half_width)
      density_scaled <- c(0, 0.38, 0)
      density_method <- "point_mass"
    } else {
      estimate <- stats::density(
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
    baseline <- ranges$density_ybaseline[[index]]
    density_rows[[index]] <- data.frame(
      parameter = parameter,
      parameter_order = index,
      row_y = ranges$row_y[[index]],
      x_plot = grid_x,
      density_scaled = density_scaled,
      y_baseline = baseline,
      y_density = baseline + density_scaled,
      density_method = density_method,
      stringsAsFactors = FALSE
    )
    height_at_best <- stats::approx(
      grid_x,
      density_scaled,
      xout = ranges$best_seed_plot[[index]],
      rule = 2
    )$y
    best_point_y[[index]] <- baseline + min(
      0.22,
      max(0.05, 0.55 * height_at_best)
    )
  }
  density_table <- do.call(rbind, density_rows)
  row.names(density_table) <- NULL
  ranges$best_seed_y <- best_point_y

  validation <- data.frame(
    metric = c(
      "seed_universe",
      "seed_count",
      "fitted_parameter_count",
      "endpoint_rows",
      "finite_endpoint_rows",
      "endpoints_outside_bounds",
      "best_seed",
      "best_seed_objective_rank",
      "range_source_seed",
      "range_source_md5",
      "display_scale",
      "zero_floor_raw",
      "upper_glyph_semantics",
      "lower_glyph_semantics",
      "uncertainty_interpretation"
    ),
    value = c(
      "seed1_through_seed500_exact",
      500L,
      nrow(ranges),
      nrow(endpoints),
      sum(is.finite(endpoints$value)),
      sum(outside_bound),
      best_row$seed[[1L]],
      best_row$objective_rank[[1L]],
      INVITRO_VISUALIZATION_SEED,
      unname(tools::md5sum(parameter_path)),
      "original_raw_value_on_log10_axis",
      format(log_floor_raw, scientific = TRUE),
      "kernel_density_of_500_optimizer_endpoints",
      "configured_lower_to_upper_range_not_empirical_quartiles",
      "numerical_search_endpoints_not_posterior_or_biological_replicates"
    ),
    stringsAsFactors = FALSE
  )

  output_paths <- c(
    ranges = write_intermediate_tsv(
      ranges,
      file.path(output_dir, "figure3g_parameter_ranges.tsv")
    ),
    endpoints = write_intermediate_tsv(
      endpoints,
      file.path(output_dir, "figure3g_parameter_endpoints_500seeds.tsv")
    ),
    density = write_intermediate_tsv(
      density_table,
      file.path(output_dir, "figure3g_violin_density.tsv")
    ),
    validation = write_intermediate_tsv(
      validation,
      file.path(output_dir, "figure3g_validation.tsv")
    )
  )
  list(
    source_files = c(best_paths, parameter_path),
    output_files = output_paths,
    best_seed = best_row$seed[[1L]],
    fitted_parameters = ranges$parameter
  )
}

data_Figure3 <- function() {
  source_root <- file.path(INVITRO_RESULT_ROOT, INVITRO_VISUALIZATION_SEED)
  local_root <- file.path(
    DATA_ROOT, "Figure3", paste0("source_", INVITRO_VISUALIZATION_SEED)
  )
  files <- c(
    "invitro_flow_overlay.tsv",
    "invitro_distribution_summary.tsv",
    "invitro_daily_counts.tsv",
    "invitro_lineage_summary.tsv",
    "invitro_growth_loglik.tsv",
    "invitro_observed_kary.tsv",
    "best_params.tsv"
  )
  sources <- file.path(source_root, files)
  require_files(
    sources, paste("Figure 3 approved", INVITRO_VISUALIZATION_SEED, "input")
  )
  destinations <- file.path(local_root, files)
  copied <- Map(copy_input, sources, destinations)
  copied_paths <- unlist(copied, use.names = FALSE)

  distribution <- utils::read.delim(
    file.path(local_root, "invitro_distribution_summary.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  observed <- utils::read.delim(
    file.path(local_root, "invitro_observed_kary.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  figure3_require_columns(
    distribution, c("N", "fraction"), "Figure 3 distribution"
  )
  figure3_require_columns(
    observed, "observed_kary_N", "Figure 3 observed karyotype"
  )
  predicted_q01 <- weighted_quantile(
    distribution$N, distribution$fraction, 0.01
  )
  predicted_q99 <- weighted_quantile(
    distribution$N, distribution$fraction, 0.99
  )
  observed_values <- observed$observed_kary_N[
    is.finite(observed$observed_kary_N)
  ]
  observed_q01 <- unname(stats::quantile(
    observed_values, 0.01, names = FALSE, type = 1
  ))
  observed_q99 <- unname(stats::quantile(
    observed_values, 0.99, names = FALSE, type = 1
  ))
  display_lower <- floor(min(predicted_q01, observed_q01)) - 2
  display_upper <- ceiling(max(predicted_q99, observed_q99)) + 2
  total_mass <- sum(distribution$fraction, na.rm = TRUE)
  display_ranges <- data.frame(
    variable = "chromosome_state_N",
    lower_quantile = 0.01,
    upper_quantile = 0.99,
    predicted_q01 = predicted_q01,
    predicted_q99 = predicted_q99,
    observed_q01 = observed_q01,
    observed_q99 = observed_q99,
    display_lower = display_lower,
    display_upper = display_upper,
    predicted_mass_below_display = sum(
      distribution$fraction[distribution$N < display_lower], na.rm = TRUE
    ) / total_mass,
    predicted_mass_above_display = sum(
      distribution$fraction[distribution$N > display_upper], na.rm = TRUE
    ) / total_mass,
    observed_count_below_display = sum(observed_values < display_lower),
    observed_count_above_display = sum(observed_values > display_upper),
    display_rule = "pooled weighted predicted and observed 1st-99th percentiles plus 2 chromosomes",
    clipping_rule = "coordinate zoom only; source rows retained",
    stringsAsFactors = FALSE
  )
  write_intermediate_tsv(
    display_ranges,
    file.path(DATA_ROOT, "Figure3", "figure3_display_ranges.tsv")
  )

  invitro_metrics <- optimizer_collect_seed_metrics(
    INVITRO_RESULT_ROOT, "objective_total"
  )
  selected_optimizer <- optimizer_write_bundle(
    metrics = invitro_metrics,
    output_dir = file.path(DATA_ROOT, "Figure3", "optimizer"),
    fit_label = "in vitro",
    selected_seed = INVITRO_VISUALIZATION_SEED,
    run_dir = INVITRO_RESULT_ROOT,
    include_all_boundaries = FALSE
  )
  optimizer_sources <- attr(selected_optimizer, "source_files")
  parameter_ensemble <- build_figure3g_parameter_ensemble(
    run_root = INVITRO_RESULT_ROOT,
    metrics = invitro_metrics,
    output_dir = file.path(DATA_ROOT, "Figure3")
  )

  input_contract <- data.frame(
    role = paste(
      "approved separate in-vitro", INVITRO_VISUALIZATION_SEED,
      "post-fit input"
    ),
    source = normalizePath(sources, mustWork = TRUE),
    local_file = copied_paths,
    source_md5 = unname(tools::md5sum(sources)),
    local_md5 = unname(tools::md5sum(copied_paths)),
    stringsAsFactors = FALSE
  )
  optimizer_contract <- data.frame(
    role = "in-vitro optimizer diagnostic input",
    source = optimizer_sources,
    local_file = NA_character_,
    source_md5 = unname(tools::md5sum(optimizer_sources)),
    local_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  parameter_contract <- data.frame(
    role = "in-vitro 500-seed fitted-parameter endpoint input",
    source = parameter_ensemble$source_files,
    local_file = NA_character_,
    source_md5 = unname(tools::md5sum(parameter_ensemble$source_files)),
    local_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  contract <- rbind(
    input_contract,
    optimizer_contract,
    parameter_contract
  )
  contract <- contract[!duplicated(contract$source), , drop = FALSE]
  write_data_contract("Figure3", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) data_Figure3()
