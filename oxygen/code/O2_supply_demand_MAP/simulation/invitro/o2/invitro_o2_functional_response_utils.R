# Oxygen-domain response tables evaluated from fitted parameters.

ivt_sim_o2_grid <- function(o2_crit, grid_min = 0, grid_max = 5) {
  base <- seq(grid_min, grid_max, by = 0.02)
  crit <- suppressWarnings(as.numeric(o2_crit[[1L]]))
  if (!is.finite(crit) || crit <= 0) return(base)
  dense <- c(
    seq(grid_min, min(grid_max, max(0.02, 25 * crit)), length.out = 220L),
    crit * 10^seq(-4, 2, length.out = 220L),
    crit * c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 25, 50, 100)
  )
  grid <- c(base, dense, grid_min, grid_max)
  sort(unique(signif(
    grid[is.finite(grid) & grid >= grid_min & grid <= grid_max],
    14
  )))
}

ivt_sim_reference_states <- function(cfg) {
  n_unit <- suppressWarnings(as.numeric(cfg$N_UNIT))
  if (!is.finite(n_unit) || n_unit <= 0) n_unit <- 22
  multiple <- seq(1.5, 5, by = 0.5)
  label <- ifelse(
    abs(multiple - round(multiple)) < 1e-8,
    paste0(as.integer(round(multiple)), "N"),
    paste0(format(multiple, trim = TRUE, nsmall = 1), "N")
  )
  data.frame(
    cohort = label,
    ploidy_multiple = multiple,
    N_ref = multiple * n_unit,
    stringsAsFactors = FALSE
  )
}

ivt_sim_compute_o2_response_tables <- function(run_params, cfg, context = NULL) {
  if (is.null(context)) {
    context <- ivt_sim_build_functional_response_context(run_params, cfg)
  }
  o2_grid <- ivt_sim_o2_grid(context$o2_crit)
  ref_states <- ivt_sim_reference_states(cfg)
  ref_primary <- data.frame(
    cohort = c("2N", "4N"),
    N_ref = c(2, 4) * as.numeric(cfg$N_UNIT),
    stringsAsFactors = FALSE
  )

  build_o2_table <- function(ref_df) {
    dplyr::bind_rows(lapply(seq_len(nrow(ref_df)), function(i) {
      rates <- context$compute_rates(
        o2_grid,
        rep(ref_df$N_ref[[i]], length(o2_grid))
      )
      out <- data.frame(
        oxygen_pct = o2_grid,
        cohort = ref_df$cohort[[i]],
        N_ref = ref_df$N_ref[[i]],
        ms_rate = pmax(as.numeric(.pmisseg_of_O2(
          O2 = o2_grid,
          run_params = run_params,
          N = ref_df$N_ref[[i]],
          O2_crit = context$o2_crit
        )), 0),
        rates[, setdiff(names(rates), c("O2", "N")), drop = FALSE],
        row.names = NULL,
        check.names = FALSE
      )
      if ("ploidy_multiple" %in% names(ref_df)) {
        out$ploidy_multiple <- ref_df$ploidy_multiple[[i]]
      }
      out
    }))
  }

  list(
    functional_curve_oxygen = build_o2_table(ref_primary),
    functional_curve_oxygen_multi_ploidy = build_o2_table(ref_states)
  )
}
