.load_invitro_observed_feasibility_helpers <- function() {
  env <- new.env(parent = baseenv())
  sys.source(
    file.path(
      repo_info$root,
      "oxygen",
      "scripts",
      "invitro_observed_feasibility_utils.R"
    ),
    envir = env
  )
  env
}

testthat::test_that("scan parameter specs use canonical order and calibrated baseline cap", {
  env <- .load_invitro_observed_feasibility_helpers()
  symbols <- c(
    "n_O", "p_mis_base", "buffer_smax", "O2_crit", "gamma_mu",
    "p_misseg", "buffer_beta", "k_o_mis", "buffer_n_exp", "unused"
  )
  table <- data.frame(
    param_symbol = symbols,
    lower_bound = c(1.5, 1e-6, 0, 0.5, 1.2, 0.005, 0.01, 0.001, 0.1, -1),
    upper_bound = c(3, 0.005, 1, 2, 3.5, 0.4, 10, 0.02, 10, 1),
    stringsAsFactors = FALSE
  )
  specs <- env$ivt_scan_parameter_specs(table, p_mis_base_upper = 2e-4)

  testthat::expect_identical(
    specs$param_symbol,
    c(
      "p_mis_base", "p_misseg", "k_o_mis", "gamma_mu", "buffer_smax",
      "buffer_beta", "buffer_n_exp", "O2_crit", "n_O"
    )
  )
  testthat::expect_equal(specs$upper[specs$param_symbol == "p_mis_base"], 2e-4)
  testthat::expect_identical(
    specs$transform,
    c("log", "log", "log", "linear", "linear", "log", "log", "log", "linear")
  )
  testthat::expect_error(
    env$ivt_scan_parameter_specs(table, p_mis_base_upper = 0.01),
    "within the parameter-table"
  )
  testthat::expect_error(
    env$ivt_scan_parameter_specs(table[table$param_symbol != "n_O", ], 2e-4),
    "exactly one row"
  )
  snapped <- env$ivt_scan_parameter_specs(table, p_mis_base_upper = 1e-6 * (1 + 1e-12))
  testthat::expect_identical(
    snapped$upper[snapped$param_symbol == "p_mis_base"],
    snapped$lower[snapped$param_symbol == "p_mis_base"]
  )
})

testthat::test_that("Latin hypercube design is deterministic, stratified, and RNG-neutral", {
  env <- .load_invitro_observed_feasibility_helpers()
  specs <- data.frame(
    param_symbol = c("linear_x", "log_x", "fixed_x"),
    lower = c(0, 1e-4, 7),
    upper = c(1, 1e-1, 7),
    transform = c("linear", "log", "linear"),
    stringsAsFactors = FALSE
  )

  set.seed(91)
  before <- .Random.seed
  design_a <- env$ivt_scan_lhs_design(20, specs, seed = 17)
  testthat::expect_identical(.Random.seed, before)
  design_b <- env$ivt_scan_lhs_design(20, specs, seed = 17)
  testthat::expect_identical(design_a, design_b)
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  design_other_kind <- env$ivt_scan_lhs_design(20, specs, seed = 17)
  testthat::expect_identical(design_a, design_other_kind)
  testthat::expect_identical(design_a$candidate_id, sprintf("candidate_%04d", 1:20))
  testthat::expect_true(all(design_a$linear_x > 0 & design_a$linear_x < 1))
  testthat::expect_true(all(design_a$log_x > 1e-4 & design_a$log_x < 1e-1))
  testthat::expect_equal(design_a$fixed_x, rep(7, 20))
  testthat::expect_equal(sort(floor(design_a$linear_x * 20)), 0:19)
  log_u <- (log(design_a$log_x) - log(1e-4)) / (log(1e-1) - log(1e-4))
  testthat::expect_equal(sort(floor(log_u * 20)), 0:19)
})

testthat::test_that("control flow high mass is grouped and uses equal-grid probability mass", {
  env <- .load_invitro_observed_feasibility_helpers()
  one <- data.frame(
    segment_id = "2N-C-A12",
    passage_id = "passage_a",
    sample_name = "sample_a",
    series = "Predicted",
    ploidy = 1:4,
    density = c(1, 1, 2, 2),
    stringsAsFactors = FALSE
  )
  two <- transform(one, passage_id = "passage_b", sample_name = "sample_b", density = c(3, 1, 0, 0))
  distractors <- rbind(
    transform(one, series = "Observed"),
    transform(one, segment_id = "2N-O1-A12")
  )
  out <- env$ivt_scan_flow_group_metrics(rbind(one, two, distractors), lower = 3)

  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_equal(out$high_mass_ge_lower[out$passage_id == "passage_a"], 4 / 6)
  testthat::expect_equal(out$high_mass_ge_lower[out$passage_id == "passage_b"], 0)
  testthat::expect_equal(out$grid_step, c(1, 1))
  testthat::expect_error(
    env$ivt_scan_flow_group_metrics(transform(one, ploidy = c(1, 2, 4, 5))),
    "equally spaced"
  )
})

testthat::test_that("flow shape metrics detect a separated approximately-WGD peak", {
  env <- .load_invitro_observed_feasibility_helpers()
  x <- seq(1, 6, by = 0.5)
  y <- rep(0.1, length(x))
  y[x == 2] <- 10
  y[x == 4] <- 5
  metrics <- env$ivt_scan_flow_shape_metrics(data.frame(ploidy = x, density = y))

  testthat::expect_equal(metrics$low_peak, 2)
  testthat::expect_equal(metrics$high_peak, 4)
  testthat::expect_equal(metrics$peak_ratio, 2)
  testthat::expect_true(metrics$valley_ratio < 0.5)
  testthat::expect_true(metrics$has_target_bimodality)
  testthat::expect_equal(
    metrics$bridge_mass + metrics$wgd_mass +
      sum(y[x < 2.5]) / sum(y),
    1
  )
})

testthat::test_that("segment distribution metrics summarize latent and smoothed states", {
  env <- .load_invitro_observed_feasibility_helpers()
  run <- list(
    grid_pre = c(40, 70, 90),
    segment_results = list(list(
      segment = list(
        segment_id = "2N-O1-A19", cohort = "2N", lineage_id = "O1",
        scenario_id = "2N-O1", oxygen_pct = 0.2
      ),
      selection = list(selected_frac = c(0.2, 0.3, 0.5))
    ))
  )
  latent <- env$ivt_scan_segment_distribution_metrics(run, "2N-O1-A19")
  smoothed <- env$ivt_scan_segment_distribution_metrics(run, "2N-O1-A19", sigma_kary = 1)

  testthat::expect_equal(latent$mean_N, 74)
  testthat::expect_equal(latent$mass_N_ge_70, 0.8)
  testthat::expect_equal(latent$median_N, 70)
  testthat::expect_identical(latent$distribution_basis, "latent")
  testthat::expect_identical(smoothed$distribution_basis, "observation_smoothed")
  testthat::expect_true(is.finite(smoothed$mean_N))
  testthat::expect_error(
    env$ivt_scan_segment_distribution_metrics(run, "missing"),
    "exactly one segment"
  )
})

testthat::test_that("calibration cap stops at the first failure and reports non-monotonic reentry", {
  env <- .load_invitro_observed_feasibility_helpers()
  calibration <- data.frame(
    p_mis_base = c(1e-6, 2e-6, 3e-6, 4e-6),
    high_mass_ge_lower = c(0.02, 0.04, 0.07, 0.03),
    status = "OK",
    protocol_feasibility_status = "PASS",
    stringsAsFactors = FALSE
  )
  cap <- env$ivt_scan_select_calibration_cap(calibration, limit = 0.05)

  testthat::expect_equal(cap$p_mis_base_upper, 2e-6)
  testthat::expect_equal(cap$n_prefix_pass, 2L)
  testthat::expect_equal(cap$first_failing_p_mis_base, 3e-6)
  testthat::expect_false(cap$monotone_nondecreasing)
  testthat::expect_true(cap$pass_reentry_after_failure)

  calibration$high_mass_ge_lower <- 0.08
  no_cap <- env$ivt_scan_select_calibration_cap(calibration, limit = 0.05)
  testthat::expect_true(is.na(no_cap$p_mis_base_upper))
  testthat::expect_equal(no_cap$n_prefix_pass, 0L)
})

testthat::test_that("hard gate enforces protocol, prediction, scenario, and control checks", {
  env <- .load_invitro_observed_feasibility_helpers()
  good <- data.frame(
    status = "OK",
    protocol_feasibility_status = "PASS",
    n_scenarios = 6,
    n_insufficient_boundaries = 0,
    n_growth_missing_pred = 0,
    n_growth_negative_pred = 0,
    control_high_mass = 0.05,
    stringsAsFactors = FALSE
  )
  testthat::expect_true(env$ivt_scan_hard_gate(good))
  testthat::expect_false(env$ivt_scan_hard_gate(transform(good, control_high_mass = 0.051)))
  testthat::expect_false(env$ivt_scan_hard_gate(transform(good, n_scenarios = 5)))
  testthat::expect_false(env$ivt_scan_hard_gate(transform(good, n_growth_missing_pred = 1)))
  alias <- good
  names(alias)[names(alias) == "control_high_mass"] <- "control_flow_mass_ge3"
  testthat::expect_true(env$ivt_scan_hard_gate(alias))
  testthat::expect_false(env$ivt_scan_hard_gate(transform(good, all_passage_executed = FALSE)))
  testthat::expect_error(
    env$ivt_scan_hard_gate(good[, names(good) != "n_growth_negative_pred"]),
    "Negative-prediction"
  )
})

testthat::test_that("Pareto front keeps ties and excludes dominated or invalid rows", {
  env <- .load_invitro_observed_feasibility_helpers()
  candidates <- data.frame(
    metric_a = c(1, 2, 3, 3, 1, NA),
    metric_b = c(3, 2, 1, 3, 3, 0)
  )
  front <- env$ivt_scan_pareto_front(candidates, c("metric_a", "metric_b"))
  testthat::expect_identical(front, c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE))
})
