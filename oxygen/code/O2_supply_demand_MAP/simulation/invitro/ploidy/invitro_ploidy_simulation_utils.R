# Ploidy-distribution outputs derived from fitted in-vitro simulations.

ivt_sim_collect_ploidy_tables <- function(components) {
  run_2n <- components$run_2N
  run_4n <- components$run_4N
  if (!is.list(run_2n) || !is.list(run_4n)) {
    stop("In-vitro best_components must contain run_2N and run_4N.")
  }
  list(
    invitro_distribution_summary = dplyr::bind_rows(
      ivt_collect_distribution_summary(run_2n),
      ivt_collect_distribution_summary(run_4n)
    ),
    invitro_distribution_quantiles = dplyr::bind_rows(
      ivt_collect_distribution_quantiles(
        run_2n,
        probs = seq(0.01, 0.99, length.out = 50L)
      ),
      ivt_collect_distribution_quantiles(
        run_4n,
        probs = seq(0.01, 0.99, length.out = 50L)
      )
    )
  )
}

ivt_sim_o2_reference_levels <- function(o2_crit, grid_min = 0, grid_max = 5) {
  crit <- suppressWarnings(as.numeric(o2_crit[[1L]]))
  levels <- c(seq(grid_min, grid_max, by = 0.5), grid_min, grid_max)
  if (is.finite(crit) && crit > 0) {
    levels <- c(
      levels,
      crit * c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 25)
    )
  }
  sort(unique(signif(
    levels[is.finite(levels) & levels >= grid_min & levels <= grid_max],
    12
  )))
}

ivt_sim_compute_ploidy_response_tables <- function(run_params, cfg, context = NULL) {
  if (is.null(context)) {
    context <- ivt_sim_build_functional_response_context(run_params, cfg)
  }
  n_states <- seq.int(as.integer(cfg$N_MIN), as.integer(cfg$N_MAX))
  ratio <- (2 * context$n_unit) / pmax(n_states, 1e-12)
  viability <- pmax(pmin(
    context$buffer_smax * exp(
      -context$buffer_beta * pmax(ratio, 0)^context$buffer_n_exp
    ),
    1
  ), 0)
  viability_table <- data.frame(
    N = n_states,
    ploidy = n_states / context$n_unit,
    endpoint_value = if (identical(context$start_mode, "chr_number")) {
      n_states
    } else {
      n_states / context$n_unit
    },
    viability_after_ms = viability,
    stringsAsFactors = FALSE
  )

  state_dense <- if (identical(context$start_mode, "chr_number")) {
    seq(as.numeric(cfg$N_MIN), as.numeric(cfg$N_MAX), by = 1)
  } else {
    seq(0, 10, by = 0.05)
  }
  o2_levels <- ivt_sim_o2_reference_levels(context$o2_crit)
  ploidy_o2 <- dplyr::bind_rows(lapply(o2_levels, function(o2) {
    n_grid <- if (identical(context$start_mode, "chr_number")) {
      state_dense
    } else {
      state_dense * context$n_unit
    }
    rates <- context$compute_rates(rep(o2, length(n_grid)), n_grid)
    data.frame(
      oxygen_pct = o2,
      ploidy = if (identical(context$start_mode, "chr_number")) {
        n_grid / context$n_unit
      } else {
        state_dense
      },
      N = n_grid,
      endpoint_value = if (identical(context$start_mode, "chr_number")) {
        n_grid
      } else {
        state_dense
      },
      proliferation_rate = rates$proliferation_rate,
      death_rate = rates$death_rate,
      net_growth_rate = rates$net_growth_rate,
      stringsAsFactors = FALSE
    )
  }))

  list(
    functional_curve_ploidy = viability_table,
    functional_curve_ploidy_by_o2 = ploidy_o2
  )
}
