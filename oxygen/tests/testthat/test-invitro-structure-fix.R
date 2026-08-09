.load_invitro_structure_api <- function() {
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  env <- new.env(parent = globalenv())
  sys.source(
    file.path(workflow_root, "util", "o2_supply_demand_map_shared.R"),
    envir = env,
    chdir = TRUE
  )
  sys.source(
    file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R"),
    envir = env,
    chdir = TRUE
  )
  source(
    file.path(workflow_root, "util", "o2_supply_demand_map_invitro_utils.R"),
    local = env,
    chdir = TRUE
  )
  env
}

.make_invitro_diagnostic_sim_args <- function(init_cells = 100,
                                               n_steps = 20L,
                                               dt = 0.1,
                                               lam_max = 0,
                                               p_mis_base = 0,
                                               mu_hp = 0,
                                               k_clear = 0,
                                               assigned_o2 = 0) {
  list(
    init_state = as.numeric(init_cells),
    N0min = 44L,
    N0max = 44L,
    N1min = 44L,
    N1max = 44L,
    obs_steps = 0:as.integer(n_steps),
    sim_end_step = as.integer(n_steps),
    DT = as.numeric(dt),
    dose = 0,
    dose_ref = 30,
    treat_day = Inf,
    fit_treatment = FALSE,
    alpha = 0,
    gamma = 1,
    tx_mult_min = 0,
    crowding_enabled = FALSE,
    crowding = "logistic",
    K = 1e12,
    min_pop = 1e-12,
    O2_crit = 1,
    o2_feedback = FALSE,
    o2_S0 = as.numeric(assigned_o2),
    kappa_O = 1,
    tau_O2 = 2,
    o2_Nref = 1e6,
    o2_min = as.numeric(assigned_o2),
    o2_S0_upper_bound = 21,
    eta_o2 = 1,
    o2_cache_bin_pct = 0.001,
    o2_cache_hysteresis_pct = 0,
    o2_cache_profile = FALSE,
    lam_max = as.numeric(lam_max),
    p_mis_base = as.numeric(p_mis_base),
    p_misseg = 0,
    k_o_mis = 1,
    p_wgd = 0,
    boundary = "drop",
    eps_tail = 0,
    buffer_smax = 1,
    buffer_beta = 0,
    buffer_n_exp = 1,
    N_unit = 22L,
    beta_size = 0,
    O2_growth = FALSE,
    alpha_o2 = 0,
    gamma_growth = 1,
    mu_hp = as.numeric(mu_hp),
    gamma_mu = 1,
    n_O = 1,
    ploidy_O2_death = "ploidy_related",
    start_with = "chr_number",
    k_clear = as.numeric(k_clear),
    vol_by_N = 1,
    burden_floor = 1e-12,
    return_full_trajectory = TRUE
  )
}

testthat::test_that("formal passage adapters form six independent true-duration scenarios", {
  env <- .load_invitro_structure_api()
  oxygen_root <- file.path(repo_info$root, "oxygen")
  fit_objects <- env$ivt_load_fit_objects(oxygen_root)
  adapters <- list(
    env$ivt_build_job_table_adapter(fit_objects$jobs_2N, fit_objects$fit_data, "2N"),
    env$ivt_build_job_table_adapter(fit_objects$jobs_4N, fit_objects$fit_data, "4N")
  )
  segments <- c(adapters[[1]]$segments, adapters[[2]]$segments)
  passage_map <- dplyr::bind_rows(lapply(adapters, env$ivt_make_passage_map))

  testthat::expect_length(segments, 114L)
  testthat::expect_identical(
    sort(unique(passage_map$scenario_id)),
    sort(c("2N-C", "2N-O1", "2N-O2", "4N-C", "4N-O1", "4N-O2"))
  )
  cumulative <- tapply(passage_map$passage_duration, passage_map$scenario_id, sum)
  testthat::expect_equal(
    as.numeric(cumulative[c("2N-C", "4N-C", "2N-O1", "2N-O2", "4N-O1", "4N-O2")]),
    c(29, 29, 122, 122, 125, 121)
  )
  testthat::expect_true(all(
    passage_map$duration_days == passage_map$observed_passage_duration + 1
  ))
  testthat::expect_true(all(
    passage_map$endpoint_day == passage_map$search_horizon_day
  ))
  testthat::expect_true(all(passage_map$passage_time_tolerance_days == 1))
  testthat::expect_true(all(vapply(segments, function(seg) length(seg$data_ids) == 1L, logical(1))))

  initial_records <- c(
    adapters[[1]]$initial_observations,
    adapters[[2]]$initial_observations
  )
  testthat::expect_setequal(
    vapply(initial_records, `[[`, character(1), "passage_id"),
    c("SUM-159_NLS_2N_A5_seed", "SUM-159_NLS_4N_A2_seed")
  )
  testthat::expect_true(all(vapply(
    initial_records,
    function(record) record$source_depth == 0L,
    logical(1)
  )))

  landmark_records <- c(
    adapters[[1]]$landmark_observations,
    adapters[[2]]$landmark_observations
  )
  landmark_targets <- stats::setNames(
    vapply(landmark_records, `[[`, character(1), "target_segment_id"),
    vapply(landmark_records, `[[`, character(1), "passage_id")
  )
  expected_landmark_targets <- c(
    "SUM-159_NLS_2N_A6M_seed" = "2N-C-A1",
    "SUM-159_NLS_4N_A4M_seed" = "4N-C-A2"
  )
  testthat::expect_identical(
    unname(landmark_targets[names(expected_landmark_targets)]),
    unname(expected_landmark_targets)
  )
  testthat::expect_true(all(vapply(
    landmark_records,
    function(record) record$source_depth > 0L && record$lineage_id == "C",
    logical(1)
  )))

  child_segments <- Filter(
    function(seg) !is.na(seg$parent_segment_id) && nzchar(seg$parent_segment_id),
    segments
  )
  testthat::expect_true(all(vapply(child_segments, function(seg) {
    startsWith(seg$parent_segment_id, paste0(seg$scenario_id, "-A"))
  }, logical(1))))
  testthat::expect_false(any(vapply(child_segments, function(seg) {
    seg$lineage_id == "O1" && grepl("-O2-", seg$parent_segment_id, fixed = TRUE)
  }, logical(1))))
  testthat::expect_false(any(vapply(child_segments, function(seg) {
    seg$lineage_id == "O2" && grepl("-O1-", seg$parent_segment_id, fixed = TRUE)
  }, logical(1))))
})

testthat::test_that("custom observation grids are clamped to each passage search horizon", {
  env <- .load_invitro_structure_api()
  oxygen_root <- file.path(repo_info$root, "oxygen")
  fit_objects <- env$ivt_load_fit_objects(oxygen_root)
  build_adapter <- get(
    ".ivt_build_independent_scenario_adapter",
    envir = env,
    inherits = FALSE
  )
  adapter <- build_adapter(
    jobs = fit_objects$jobs_2N,
    fit_data = fit_objects$fit_data,
    cohort = "2N",
    obs_days_local = c(-2, 0:14, Inf, NA_real_)
  )

  testthat::expect_true(all(vapply(adapter$segments, function(segment) {
    days <- segment$obs_days_local
    all(is.finite(days)) &&
      min(days) == 0 &&
      max(days) == segment$endpoint_day &&
      all(days <= segment$endpoint_day) &&
      sum(abs(days - segment$endpoint_day) <= 1e-10) == 1L
  }, logical(1))))
  testthat::expect_true(any(vapply(
    adapter$segments,
    function(segment) segment$search_horizon_day < 14 && !14 %in% segment$obs_days_local,
    logical(1)
  )))
})

testthat::test_that("non-root landmarks use their source segment endpoint distribution", {
  env <- .load_invitro_structure_api()
  segment <- list(
    segment_id = "2N-C-A1",
    cohort = "2N",
    lineage_id = "C",
    lineage_group = "control",
    lineage_label = "C",
    scenario_id = "2N-C",
    lineage_terminal_key = "2N-C-A1",
    passage_index = 1L,
    lineage_passage_index = 1L,
    oxygen_pct = 20.5,
    data_ids = "formal"
  )
  run <- list(
    grid_pre = c(22, 44),
    segment_results = stats::setNames(
      list(list(
        segment = segment,
        selection = list(selected_frac = c(0, 1))
      )),
      segment$segment_id
    ),
    initial_fraction = c(1, 0),
    initial_observations = list(),
    landmark_observations = list(list(
      passage_id = "landmark",
      target_segment_id = segment$segment_id
    ))
  )
  fit_data <- list(
    formal = list(kary = numeric()),
    landmark = list(kary = 44)
  )

  loglik <- env$ivt_ploidy_loglik_df(
    run = run,
    fit_data = fit_data,
    sigma_kary = 0.1
  )
  landmark_row <- loglik[loglik$passage_id == "landmark", , drop = FALSE]
  testthat::expect_identical(nrow(landmark_row), 1L)
  testthat::expect_identical(landmark_row$segment_id, segment$segment_id)
  testthat::expect_identical(landmark_row$lineage_id, "C")
  testthat::expect_gt(landmark_row$mean_loglik, -1e-8)
})

testthat::test_that("passage selection never precedes the last observation and never upsamples", {
  env <- .load_invitro_structure_api()
  endpoint_state <- c(3, 5)
  sim <- list(
    Ntot_live_obs = c(10, 12, 9, sum(endpoint_state)),
    live_state_obs = rbind(
      c(6, 4),
      c(7, 5),
      c(4, 5),
      endpoint_state
    )
  )
  closest_day_one <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = 9,
    grid_pre = c(22, 44),
    target_live_cells = 12,
    obs_days_local = 0:3,
    observed_passage_day = 1
  )
  closest_endpoint <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = 9,
    grid_pre = c(22, 44),
    target_live_cells = 8,
    obs_days_local = 0:3,
    observed_passage_day = 2
  )

  testthat::expect_identical(closest_day_one$closest_day_diagnostic, 1)
  testthat::expect_identical(closest_endpoint$closest_day_diagnostic, 3)
  testthat::expect_identical(closest_day_one$selected_day, 1)
  testthat::expect_identical(closest_endpoint$selected_day, 2)
  testthat::expect_true(closest_day_one$selected_day_after_last_observation)
  testthat::expect_true(closest_endpoint$selected_day_after_last_observation)
  testthat::expect_identical(closest_day_one$endpoint_state, c(7, 5))
  testthat::expect_equal(closest_day_one$reseeded_state, c(7, 5) * 0.75)
  testthat::expect_identical(closest_endpoint$reseeded_state, c(4, 5))
  testthat::expect_identical(closest_day_one$reseed_mode, "downsample_to_observed_inoculum")
  testthat::expect_equal(closest_day_one$boundary_scale, 0.75)
  testthat::expect_identical(closest_day_one$cell_number_before, 12)
  testthat::expect_identical(closest_day_one$cell_number_after, 9)

  downsampled <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = 4,
    grid_pre = c(22, 44),
    target_live_cells = 12,
    obs_days_local = 0:3,
    observed_passage_day = 1
  )
  testthat::expect_identical(downsampled$reseed_mode, "downsample_to_observed_inoculum")
  testthat::expect_equal(downsampled$boundary_scale, 1 / 3)
  testthat::expect_equal(downsampled$reseeded_state, c(7, 5) / 3)
  testthat::expect_true(downsampled$boundary_scale <= 1)

  first_eligible <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 15, 12, 10),
      live_state_obs = matrix(c(5, 15, 12, 10), ncol = 1)
    ),
    reseed_live_cells = 9,
    grid_pre = 44,
    target_live_cells = 10,
    obs_days_local = 0:3,
    observed_passage_day = 1,
    passage_time_tolerance_days = 2
  )
  testthat::expect_identical(first_eligible$closest_day_diagnostic, 3)
  testthat::expect_identical(first_eligible$selected_day, 1)
  testthat::expect_identical(first_eligible$selected_live_cells, 15)
  testthat::expect_identical(first_eligible$effective_threshold_cells, 10)

  testthat::expect_identical(closest_endpoint$effective_threshold_cells, 9)
  testthat::expect_match(
    closest_endpoint$threshold_target_source,
    "next_initial"
  )

  reached_early <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 10, 12, 15, 18, 20, 22),
      live_state_obs = matrix(c(5, 10, 12, 15, 18, 20, 22), ncol = 1)
    ),
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 10,
    obs_days_local = 0:6,
    observed_passage_day = 5,
    passage_time_tolerance_days = 1
  )
  testthat::expect_identical(reached_early$closest_day_diagnostic, 1)
  testthat::expect_identical(reached_early$last_observation_day, 5)
  testthat::expect_identical(reached_early$selected_day, 5)
  testthat::expect_identical(
    reached_early$predicted_live_cells_at_observation,
    20
  )

  reached_after_observation <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 8, 10, 12, 15, 18, 20),
      live_state_obs = matrix(c(5, 8, 10, 12, 15, 18, 20), ncol = 1)
    ),
    reseed_live_cells = 20,
    grid_pre = 44,
    target_live_cells = 20,
    obs_days_local = 0:6,
    observed_passage_day = 5,
    passage_time_tolerance_days = 1
  )
  testthat::expect_identical(reached_after_observation$selected_day, 6)
  testthat::expect_identical(
    reached_after_observation$predicted_live_cells_at_observation,
    18
  )

  early_only <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 20, 18, 15),
      live_state_obs = matrix(c(5, 20, 18, 15), ncol = 1)
    ),
    reseed_live_cells = 20,
    grid_pre = 44,
    target_live_cells = 20,
    obs_days_local = 0:3,
    observed_passage_day = 3,
    passage_time_tolerance_days = 0
  )
  testthat::expect_identical(early_only$closest_day_diagnostic, 1)
  testthat::expect_false(early_only$passage_executed)
  testthat::expect_match(
    early_only$passage_failure_reason,
    "after_last_observation"
  )

  infeasible <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = 9,
    grid_pre = c(22, 44),
    target_live_cells = 20,
    obs_days_local = 0:3,
    observed_passage_day = 2,
    passage_time_tolerance_days = 1
  )
  testthat::expect_false(infeasible$passage_executed)
  testthat::expect_identical(
    infeasible$reseed_mode,
    "no_passage_threshold_not_reached"
  )
  testthat::expect_null(infeasible$reseeded_state)
  testthat::expect_true(is.na(infeasible$boundary_scale))
  testthat::expect_match(infeasible$passage_failure_reason, "threshold=20")
})

testthat::test_that("growth likelihood uses every measured post-seeding live-cell count", {
  env <- .load_invitro_structure_api()
  summary_df <- data.frame(
    observed_live_cells_at_observation = c(100, 400),
    predicted_live_cells_at_observation = c(125, 500),
    last_observation_day = c(5, 5),
    observation_day = c(2, 5),
    stringsAsFactors = FALSE
  )
  out <- env$ivt_growth_loglik_df(summary_df, sigma_growth = 0.2)

  testthat::expect_identical(
    out$growth_likelihood_scale,
    rep("constant_sd_log_absolute_live_cells_at_all_measured_timepoints", 2)
  )
  testthat::expect_equal(out$sigma_log_live_cells, c(0.2, 0.2))
  testthat::expect_equal(
    out$loglik,
    stats::dnorm(
      log(c(100, 400)),
      mean = log(c(125, 500)),
      sd = c(0.2, 0.2),
      log = TRUE
    )
  )
  testthat::expect_equal(out$loglik[[1]], out$loglik[[2]])
})

testthat::test_that("growth likelihood weights equal log-count errors equally across time", {
  env <- .load_invitro_structure_api()
  passage_fn <- get(
    ".ivt_growth_passage_loglik_df",
    envir = env,
    inherits = FALSE
  )
  observation_day <- c(2, 5, 1, 2, 3, 4, 5)
  passage_id <- c(rep("two_points", 2), rep("five_points", 5))
  observed <- rep(100, length(observation_day))
  summary_df <- data.frame(
    passage_id = passage_id,
    cohort = "2N",
    lineage_id = "O1",
    scenario_id = "2N-O1",
    observed_live_cells_at_observation = observed,
    predicted_live_cells_at_observation = observed * exp(0.1),
    last_observation_day = c(rep(5, 2), rep(5, 5)),
    observation_day = observation_day,
    stringsAsFactors = FALSE
  )

  growth <- env$ivt_growth_loglik_df(summary_df, sigma_growth = 0.2)
  passage <- passage_fn(growth)

  testthat::expect_equal(length(unique(growth$loglik)), 1L)
  testthat::expect_equal(
    passage$mean_loglik[passage$passage_id == "two_points"],
    passage$mean_loglik[passage$passage_id == "five_points"]
  )
})

testthat::test_that("growth measurement matching excludes day zero and preserves passage weights", {
  env <- .load_invitro_structure_api()
  measurement_fn <- get(
    ".ivt_growth_measurement_summary",
    envir = env,
    inherits = FALSE
  )
  passage_fn <- get(
    ".ivt_growth_passage_loglik_df",
    envir = env,
    inherits = FALSE
  )
  aggregate_fn <- get(".ivt_hierarchical_loglik", envir = env, inherits = FALSE)
  summary_df <- data.frame(
    passage_id = "p1",
    segment_id = "s1",
    cohort = "2N",
    lineage_id = "O1",
    scenario_id = "2N-O1",
    last_observation_day = 5,
    selected_day = 6,
    observed_live_cells_at_observation = 400,
    predicted_live_cells_at_observation = 500,
    stringsAsFactors = FALSE
  )
  fit_data <- list(p1 = list(growth_timecourse = data.frame(
    growth_observation_id = c("p1-day0", "p1-day2", "p1-day5"),
    observation_day = c(0, 2, 5),
    observed_live_cells = c(10, 100, 400),
    growth_data_source = "unit_test",
    stringsAsFactors = FALSE
  )))
  daily_counts <- data.frame(
    segment_id = rep("s1", 4),
    passage_id = rep("p1", 4),
    day = c(0, 2, 5, 6),
    live_cells = c(10, 125, 500, 600),
    stringsAsFactors = FALSE
  )
  matched <- measurement_fn(summary_df, fit_data, daily_counts)
  testthat::expect_identical(matched$growth_observation_id, c("p1-day2", "p1-day5"))
  testthat::expect_equal(matched$observation_day, c(2, 5))
  testthat::expect_equal(matched$observed_live_cells_at_observation, c(100, 400))
  testthat::expect_equal(matched$predicted_live_cells_at_observation, c(125, 500))
  testthat::expect_identical(matched$is_last_observation, c(FALSE, TRUE))
  selected_too_early <- summary_df
  selected_too_early$selected_day <- 4
  testthat::expect_error(
    measurement_fn(selected_too_early, fit_data, daily_counts),
    "selected before its last growth observation day"
  )

  synthetic <- data.frame(
    passage_id = c("p1", "p1", "p2"),
    cohort = rep("2N", 3),
    lineage_id = rep("O1", 3),
    scenario_id = rep("2N-O1", 3),
    loglik = c(-2, -4, -9),
    stringsAsFactors = FALSE
  )
  passage <- passage_fn(synthetic)
  testthat::expect_equal(
    passage$mean_loglik[match(c("p1", "p2"), passage$passage_id)],
    c(-3, -9)
  )
  hierarchy <- aggregate_fn(passage, value_col = "mean_loglik", modality = "growth")
  testthat::expect_equal(hierarchy$value, -6)
})

testthat::test_that("fit-object loader restores complete measured growth timecourses", {
  env <- .load_invitro_structure_api()
  oxygen_root <- file.path(repo_info$root, "oxygen")
  fit_objects <- env$ivt_load_fit_objects(oxygen_root)
  a2 <- fit_objects$fit_data[["SUM-159_NLS_2N_O1_A2_seed"]]$growth_timecourse
  testthat::expect_equal(a2$observation_day, c(0, 2, 5))
  testthat::expect_equal(a2$observed_live_cells, c(404570, 832980, 6311159))

  formal_ids <- grep(
    "_(2N|4N)_(C|O1|O2)_A[0-9]+_seed$",
    names(fit_objects$fit_data),
    value = TRUE
  )
  testthat::expect_length(formal_ids, 114L)
  n_positive <- sum(vapply(formal_ids, function(passage_id) {
    sum(fit_objects$fit_data[[passage_id]]$growth_timecourse$observation_day > 0)
  }, integer(1)))
  n_terminal <- sum(vapply(formal_ids, function(passage_id) {
    entry <- fit_objects$fit_data[[passage_id]]
    any(abs(
      entry$growth_timecourse$observation_day - entry$passage_duration
    ) <= 1e-10)
  }, logical(1)))
  testthat::expect_identical(n_positive, 219L)
  testthat::expect_identical(n_terminal, 114L)
  testthat::expect_identical(
    unique(fit_objects$fit_data[["SUM-159_NLS_2N_O1_A23_seed"]]$growth_timecourse$growth_data_source),
    "fit_data_endpoint_fallback"
  )
})

testthat::test_that("threshold crossing handles all diagnostic edge cases", {
  env <- .load_invitro_structure_api()

  reached <- env$ivt_first_threshold_crossing(
    days = 0:3,
    live_cells = c(5, 8, 12, 11),
    threshold_target_cells = 10
  )
  testthat::expect_true(reached$threshold_reached_by_endpoint)
  testthat::expect_equal(reached$predicted_threshold_crossing_day, 1.5)
  testthat::expect_equal(reached$threshold_time_grid_resolution_days, 1)
  testthat::expect_equal(
    reached$threshold_crossing_interval_width_days,
    1
  )

  reached_at_start <- env$ivt_first_threshold_crossing(
    days = 0:2,
    live_cells = c(10, 8, 9),
    threshold_target_cells = 10
  )
  testthat::expect_true(
    reached_at_start$threshold_reached_by_endpoint
  )
  testthat::expect_identical(
    reached_at_start$predicted_threshold_crossing_day,
    0
  )

  not_reached <- env$ivt_first_threshold_crossing(
    days = 0:2,
    live_cells = c(2, 4, 6),
    threshold_target_cells = 10
  )
  testthat::expect_false(not_reached$threshold_reached_by_endpoint)
  testthat::expect_true(
    is.na(not_reached$predicted_threshold_crossing_day)
  )

  nonmonotonic <- env$ivt_first_threshold_crossing(
    days = 0:3,
    live_cells = c(5, 12, 4, 13),
    threshold_target_cells = 10
  )
  testthat::expect_true(nonmonotonic$threshold_reached_by_endpoint)
  testthat::expect_equal(
    nonmonotonic$predicted_threshold_crossing_day,
    5 / 7
  )

  missing_target <- env$ivt_first_threshold_crossing(
    days = 0:2,
    live_cells = c(2, 20, 30),
    threshold_target_cells = NA_real_
  )
  testthat::expect_false(
    missing_target$threshold_reached_by_endpoint
  )
  testthat::expect_true(
    is.na(missing_target$predicted_threshold_crossing_day)
  )
})

testthat::test_that("insufficient boundaries do not create a child state", {
  env <- .load_invitro_structure_api()
  endpoint_state <- as.numeric(seq_len(133L))
  initial_state <- rev(endpoint_state)
  sim <- list(
    Ntot_live_obs = c(sum(initial_state), sum(endpoint_state)),
    live_state_obs = rbind(initial_state, endpoint_state)
  )
  selection <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = sum(endpoint_state) + 1,
    grid_pre = seq_len(133L),
    target_live_cells = sum(endpoint_state) + 100,
    obs_days_local = c(0, 1)
  )
  testthat::expect_identical(
    selection$reseed_mode,
    "no_passage_threshold_not_reached"
  )
  testthat::expect_false(selection$passage_executed)
  testthat::expect_true(is.na(selection$boundary_scale))
  testthat::expect_null(selection$endpoint_state)
  testthat::expect_null(selection$reseeded_state)

  protocol <- env$.ivt_collect_protocol_feasibility(data.frame(
    cohort = "2N",
    lineage_id = "O1",
    scenario_id = "2N-O1",
    passage_id = "2N-O1-A1",
    passage_index = 1L,
    lineage_passage_index = 1L,
    passage_duration = 1,
    endpoint_day = 1,
    reseed_mode = "no_passage_threshold_not_reached",
    insufficient_boundary = TRUE,
    boundary_feasible = FALSE,
    available_cells = selection$available_cells,
    required_cells = selection$required_cells,
    supply_ratio = selection$supply_ratio,
    boundary_scale = selection$boundary_scale,
    cell_number_before = selection$cell_number_before,
    cell_number_after = selection$cell_number_after,
    stringsAsFactors = FALSE
  ))
  testthat::expect_identical(
    protocol$protocol_feasibility_status,
    "FAIL"
  )
  testthat::expect_false(
    protocol$all_passage_boundaries_feasible
  )
  testthat::expect_identical(
    protocol$n_insufficient_boundaries,
    1L
  )
})

testthat::test_that("C++ event diagnostics are cumulative and mechanism-specific", {
  null_sim <- cpp_o2simps_simulate_one(
    .make_invitro_diagnostic_sim_args(
      lam_max = 0,
      p_mis_base = 0,
      mu_hp = 0
    )
  )
  cumulative_names <- c(
    "cumulative_gross_divisions_obs",
    "cumulative_hypoxia_deaths_obs",
    "cumulative_dead_buffer_inflow_obs",
    "cumulative_nonlive_inflow_obs"
  )
  for (field in cumulative_names) {
    values <- as.numeric(null_sim[[field]])
    testthat::expect_true(all(is.finite(values)))
    testthat::expect_true(all(values == 0))
    testthat::expect_true(all(diff(values) >= 0))
  }

  division_sim <- cpp_o2simps_simulate_one(
    .make_invitro_diagnostic_sim_args(
      lam_max = 0.2,
      p_mis_base = 0,
      mu_hp = 0
    )
  )
  gross_divisions <- as.numeric(
    division_sim$cumulative_gross_divisions_obs
  )
  testthat::expect_equal(gross_divisions[[1L]], 0)
  testthat::expect_equal(gross_divisions[[2L]], 2, tolerance = 1e-12)
  testthat::expect_true(all(diff(gross_divisions) >= 0))
  testthat::expect_gt(tail(gross_divisions, 1L), 0)

  death_sim <- cpp_o2simps_simulate_one(
    .make_invitro_diagnostic_sim_args(
      n_steps = 200L,
      dt = 0.05,
      lam_max = 0,
      p_mis_base = 0,
      mu_hp = 1,
      k_clear = 2
    )
  )
  cumulative_deaths <- as.numeric(
    death_sim$cumulative_hypoxia_deaths_obs
  )
  dead_stock <- as.numeric(death_sim$Ntot_dead_hypoxia_obs)
  testthat::expect_true(all(diff(cumulative_deaths) >= -1e-12))
  testthat::expect_gt(tail(cumulative_deaths, 1L), 0)
  testthat::expect_true(any(diff(dead_stock) < 0))
  testthat::expect_equal(
    as.numeric(death_sim$cumulative_dead_buffer_inflow_obs),
    rep(0, length(cumulative_deaths)),
    tolerance = 1e-12
  )

  buffer_sim <- cpp_o2simps_simulate_one(
    .make_invitro_diagnostic_sim_args(
      lam_max = 0.2,
      p_mis_base = 0.15,
      mu_hp = 0
    )
  )
  cumulative_buffer <- as.numeric(
    buffer_sim$cumulative_dead_buffer_inflow_obs
  )
  testthat::expect_true(all(diff(cumulative_buffer) >= -1e-12))
  testthat::expect_gt(tail(cumulative_buffer, 1L), 0)
  testthat::expect_equal(
    as.numeric(buffer_sim$cumulative_hypoxia_deaths_obs),
    rep(0, length(cumulative_buffer)),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    as.numeric(buffer_sim$cumulative_nonlive_inflow_obs),
    cumulative_buffer,
    tolerance = 1e-12
  )
  testthat::expect_match(
    buffer_sim$cumulative_dead_buffer_inflow_definition,
    "missegregation-linked nonviable"
  )
})

testthat::test_that("fixed external O2 is invariant to initial cell count", {
  assigned_o2 <- 0.3
  low_count <- cpp_o2simps_simulate_one(
    .make_invitro_diagnostic_sim_args(
      init_cells = 100,
      lam_max = 0.2,
      assigned_o2 = assigned_o2
    )
  )
  high_count <- cpp_o2simps_simulate_one(
    .make_invitro_diagnostic_sim_args(
      init_cells = 1e8,
      lam_max = 0.2,
      assigned_o2 = assigned_o2
    )
  )
  for (field in c("O2_target_obs", "O2_eff_obs")) {
    testthat::expect_equal(
      as.numeric(low_count[[field]]),
      rep(assigned_o2, length(low_count[[field]])),
      tolerance = 1e-12
    )
    testthat::expect_equal(
      as.numeric(high_count[[field]]),
      rep(assigned_o2, length(high_count[[field]])),
      tolerance = 1e-12
    )
    testthat::expect_identical(
      as.numeric(low_count[[field]]),
      as.numeric(high_count[[field]])
    )
  }
})

testthat::test_that("segment simulation cannot silently rescale an explicit parent state", {
  env <- .load_invitro_structure_api()
  testthat::expect_error(
    env$ivt_run_segment_fixed_o2(
      segment = list(duration_days = 1, obs_days_local = c(0, 2)),
      cfg = list(),
      run_params = list(),
      model_core = list(),
      vol_by_N = numeric()
    ),
    "cannot extend beyond the segment search horizon"
  )
  testthat::expect_error(
    env$ivt_run_segment_fixed_o2(
      segment = list(
        oxygen_pct = 1,
        duration_days = 1,
        obs_days_local = c(0, 1)
      ),
      cfg = list(o2_burden_feedback = FALSE, o2_min = 1),
      run_params = list(o2_S0 = 1),
      model_core = list(grid_pre = c(22, 44)),
      vol_by_N = c(1, 1),
      init_state_override = c(2, 3),
      init_cells_override = 10
    ),
    "cannot rescale an explicit parent state"
  )
})

.mock_independent_adapter <- function(o1_duration = 2, identical_inputs = FALSE) {
  make_segment <- function(lineage_id, passage_index, duration, parent = NA_character_) {
    scenario_id <- paste("2N", lineage_id, sep = "-")
    initial <- if (identical_inputs) 10 else if (lineage_id == "O1") 12 else 10
    list(
      segment_id = paste0(scenario_id, "-A", passage_index),
      parent_segment_id = parent,
      cohort = "2N",
      lineage_id = lineage_id,
      lineage_group = "deprived",
      lineage_label = lineage_id,
      scenario_id = scenario_id,
      lineage_terminal_key = paste0(scenario_id, "-A2"),
      passage_id = paste0(scenario_id, "-A", passage_index),
      passage_index = passage_index,
      lineage_passage_index = passage_index,
      oxygen_pct = if (identical_inputs) 1 else if (lineage_id == "O1") 1 else 2,
      duration_days = duration,
      passage_duration = duration,
      endpoint_day = duration,
      initial_cells = initial,
      final_cells = initial,
      obs_days_local = c(0, duration),
      data_ids = paste0(scenario_id, "-A", passage_index),
      observed = list()
    )
  }
  o1_a1 <- make_segment("O1", 1L, o1_duration)
  o1_a2 <- make_segment("O1", 2L, 2, parent = o1_a1$segment_id)
  o2_a1 <- make_segment("O2", 1L, if (identical_inputs) o1_duration else 2)
  o2_a2 <- make_segment("O2", 2L, 2, parent = o2_a1$segment_id)
  list(
    cohort = "2N",
    n_segments = 4L,
    n_scenarios = 2L,
    scenario_ids = c("2N-O1", "2N-O2"),
    segments = list(o1_a1, o1_a2, o2_a1, o2_a2),
    initial_observations = list()
  )
}

testthat::test_that("an unreachable threshold stops before any child simulation", {
  env <- .load_invitro_structure_api()
  adapter <- .mock_independent_adapter(o1_duration = 2, identical_inputs = TRUE)
  adapter$segments <- adapter$segments[1:2]
  adapter$n_segments <- 2L
  adapter$n_scenarios <- 1L
  adapter$scenario_ids <- "2N-O1"
  adapter$segments[[1L]]$final_cells <- 100
  calls <- character()
  env$cell_volume_mm3_by_N <- function(grid_pre, run_params, cfg) rep(1, length(grid_pre))
  env$ivt_run_segment_fixed_o2 <- function(segment,
                                           cfg,
                                           run_params,
                                           model_core,
                                           vol_by_N,
                                           init_state_override = NULL,
                                           init_cells_override = NULL) {
    calls <<- c(calls, segment$segment_id)
    init_state <- c(as.numeric(init_cells_override), 0)
    list(
      segment = segment,
      sim = list(
        Ntot_live_obs = c(sum(init_state), 20),
        live_state_obs = rbind(init_state, c(20, 0))
      ),
      initial_state = init_state,
      initial_cells = sum(init_state)
    )
  }
  testthat::expect_error(
    env$ivt_run_lineage(
      adapter,
      cfg = list(init_total_size = 10),
      run_params = list(init_mean_2N = NA_real_, init_sd_2N = NA_real_),
      model_core = list(
        grid_pre = c(22, 44),
        init_state_2N = c(1, 0),
        init_state_4N = c(0, 1)
      )
    ),
    "protocol_infeasible:.*2N-O1-A1",
    class = "invitro_protocol_infeasible"
  )
  testthat::expect_identical(calls, "2N-O1-A1")
})

testthat::test_that("scenario propagation is independent and uses one shared parameter vector", {
  env <- .load_invitro_structure_api()
  env$cell_volume_mm3_by_N <- function(grid_pre, run_params, cfg) rep(1, length(grid_pre))
  env$ivt_run_segment_fixed_o2 <- function(segment,
                                           cfg,
                                           run_params,
                                           model_core,
                                           vol_by_N,
                                           init_state_override = NULL,
                                           init_cells_override = NULL) {
    init_state <- if (is.null(init_state_override)) {
      c(as.numeric(init_cells_override), 0)
    } else {
      as.numeric(init_state_override)
    }
    multiplier <- 1 + 0.1 * as.numeric(segment$duration_days) + 0.01 * as.numeric(segment$oxygen_pct)
    endpoint_state <- init_state * multiplier
    list(
      segment = segment,
      sim = list(
        Ntot_live_obs = c(sum(init_state), sum(endpoint_state)),
        live_state_obs = rbind(init_state, endpoint_state)
      ),
      initial_state = init_state,
      initial_cells = sum(init_state)
    )
  }
  cfg <- list(init_total_size = 10)
  run_params <- list(
    lam_max = 0.4,
    p_misseg = 0.01,
    init_mean_2N = NA_real_,
    init_sd_2N = NA_real_
  )
  model_core <- list(
    grid_pre = c(22, 44),
    init_state_2N = c(1, 0),
    init_state_4N = c(0, 1)
  )

  base_run <- env$ivt_run_lineage(
    .mock_independent_adapter(o1_duration = 2),
    cfg = cfg,
    run_params = run_params,
    model_core = model_core
  )
  changed_o1_run <- env$ivt_run_lineage(
    .mock_independent_adapter(o1_duration = 5),
    cfg = cfg,
    run_params = run_params,
    model_core = model_core
  )
  o2_keys <- c("2N-O2-A1", "2N-O2-A2")
  for (key in o2_keys) {
    testthat::expect_identical(
      base_run$segment_results[[key]]$sim$live_state_obs,
      changed_o1_run$segment_results[[key]]$sim$live_state_obs
    )
  }
  testthat::expect_false(identical(
    base_run$segment_results[["2N-O1-A1"]]$sim$live_state_obs,
    changed_o1_run$segment_results[["2N-O1-A1"]]$sim$live_state_obs
  ))

  identical_run <- env$ivt_run_lineage(
    .mock_independent_adapter(o1_duration = 2, identical_inputs = TRUE),
    cfg = cfg,
    run_params = run_params,
    model_core = model_core
  )
  for (passage_index in 1:2) {
    testthat::expect_identical(
      identical_run$segment_results[[paste0("2N-O1-A", passage_index)]]$sim$live_state_obs,
      identical_run$segment_results[[paste0("2N-O2-A", passage_index)]]$sim$live_state_obs
    )
  }
  testthat::expect_identical(identical_run$shared_run_params, run_params)
  testthat::expect_true(all(vapply(
    identical_run$adapter$segments,
    function(seg) !"run_params" %in% names(seg),
    logical(1)
  )))
})

testthat::test_that("likelihood hierarchy gives lineages equal weight within cohort", {
  env <- .load_invitro_structure_api()
  synthetic <- data.frame(
    cohort = c("2N", "2N", "2N", "2N"),
    lineage_id = c("C", "O1", "O1", "O2"),
    scenario_id = c("2N-C", "2N-O1", "2N-O1", "2N-O2"),
    loglik = c(-3, 0, 0, -6),
    stringsAsFactors = FALSE
  )
  aggregate_fn <- get(".ivt_hierarchical_loglik", envir = env, inherits = FALSE)
  result <- aggregate_fn(synthetic, value_col = "loglik", modality = "growth")
  testthat::expect_equal(result$value, -3)
  testthat::expect_equal(
    result$lineage$mean_loglik[match(c("C", "O1", "O2"), result$lineage$lineage_id)],
    c(-3, 0, -6)
  )
})

testthat::test_that("passage-time likelihood is zero inside tolerance and robust outside", {
  env <- .load_invitro_structure_api()
  summary_df <- data.frame(
    passage_id = paste0("p", 1:4),
    cohort = rep("2N", 4),
    lineage_id = rep("O1", 4),
    scenario_id = rep("2N-O1", 4),
    selected_day = c(4, 5, 2, 7),
    observed_passage_day = rep(4, 4),
    stringsAsFactors = FALSE
  )
  out <- env$ivt_passage_time_loglik_df(
    summary_df,
    passage_time_tolerance_days = 1,
    passage_time_sigma_days = 1,
    passage_time_df = 4
  )
  testthat::expect_equal(out$passage_time_residual_days, c(0, 1, -2, 3))
  testthat::expect_equal(out$passage_time_excess_days, c(0, 0, 1, 2))
  testthat::expect_equal(out$passage_time_loglik[1:2], c(0, 0))
  testthat::expect_lt(out$passage_time_loglik[[3L]], 0)
  testthat::expect_lt(out$passage_time_loglik[[4L]], out$passage_time_loglik[[3L]])
})

testthat::test_that("death likelihood cannot use a passage state before its observation", {
  env <- .load_invitro_structure_api()
  passage_id <- "synthetic_2N_O1_A1_seed"
  segment <- list(
    segment_id = "2N-O1-A1",
    cohort = "2N",
    lineage_id = "O1",
    lineage_group = "deprived",
    lineage_label = "O1",
    scenario_id = "2N-O1",
    lineage_terminal_key = "2N-O1-A1",
    passage_index = 1L,
    lineage_passage_index = 1L,
    oxygen_pct = 1,
    passage_id = passage_id,
    data_ids = passage_id,
    duration_days = 3,
    passage_duration = 2,
    observed_passage_day = 2,
    search_horizon_day = 3,
    passage_time_tolerance_days = 1,
    endpoint_day = 3,
    obs_days_local = 0:3
  )
  sim <- list(
    Ntot_live_obs = c(100, 150, 180, 200),
    Ntot_dead_total_obs = c(0, 15, 9, 40),
    Ntot_total_obs = c(100, 165, 189, 240),
    live_state_obs = matrix(c(100, 150, 180, 200), ncol = 1)
  )
  selection <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = NA_real_,
    grid_pre = 44,
    target_live_cells = 150,
    obs_days_local = 0:3,
    observed_passage_day = 2,
    passage_time_tolerance_days = 1
  )
  run <- list(
    grid_pre = 44,
    segment_results = list(list(segment = segment, sim = sim, selection = selection))
  )
  death_data <- data.frame(
    observation_id = "death-p1",
    cohort = "2N",
    lineage_id = "O1",
    scenario_id = "2N-O1",
    model_passage_id = passage_id,
    model_segment_id = "2N-O1-A1",
    likelihood_observation_day = 2,
    dead_count = 10,
    eligible_denominator = 200,
    observed_dead_fraction = 0.05,
    stringsAsFactors = FALSE
  )
  out <- env$ivt_death_loglik_df(
    run_2N = run,
    run_4N = list(grid_pre = 44, segment_results = list()),
    death_data = death_data
  )
  testthat::expect_identical(selection$selected_day, 2)
  testthat::expect_identical(out$likelihood_observation_day, 2)
  testthat::expect_identical(out$model_day, 2)
  testthat::expect_equal(out$predicted_dead_fraction, 9 / 189)
})

testthat::test_that("selected-time count plot keeps independent lineages renderable", {
  env <- .load_invitro_structure_api()
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  sys.source(
    file.path(
      workflow_root,
      "vis",
      "invitro",
      "o2_supply_demand_map_invitro_plot_utils.R"
    ),
    envir = env,
    chdir = TRUE
  )
  summary_df <- data.frame(
    cohort = "2N",
    lineage_id = rep(c("O1", "O2"), each = 3),
    lineage_label = rep(c("O1", "O2"), each = 3),
    scenario_id = rep(c("2N-O1", "2N-O2"), each = 3),
    passage_index = rep(1:3, 2),
    oxygen_pct = c(2, 1, 0.5, 2, 1, 0.5),
    predicted_live_cells = c(10, 12, 15, 10, 11, 13),
    selected_day = rep(4, 6),
    stringsAsFactors = FALSE
  )
  plot <- env$ivt_plot_lineage_counts(summary_df)
  testthat::expect_s3_class(plot, "ggplot")
  testthat::expect_silent(
    withr::with_pdf(
      tempfile(fileext = ".pdf"),
      ggplot2::ggplotGrob(plot)
    )
  )
})

testthat::test_that("post-fit tables preserve selected passage time", {
  env <- .load_invitro_structure_api()
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  sys.source(
    file.path(
      workflow_root,
      "util",
      "o2_supply_demand_map_invitro_postfit_io_utils.R"
    ),
    envir = env,
    chdir = TRUE
  )
  legacy <- data.frame(
    segment_id = "legacy-segment",
    cohort = "2N",
    lineage_label = "deprived",
    passage_duration = 4,
    selected_day = 2,
    stringsAsFactors = FALSE
  )
  normalized <- env$ivt_sim_normalize_lineage_columns(legacy)
  testthat::expect_identical(normalized$endpoint_day, 2)
  testthat::expect_identical(normalized$selected_day, 2)
  testthat::expect_true(is.na(normalized$closest_day_diagnostic))

  fixed_endpoint <- legacy
  fixed_endpoint$selected_day <- NULL
  normalized_fixed <- env$ivt_sim_normalize_lineage_columns(fixed_endpoint)
  testthat::expect_identical(normalized_fixed$endpoint_day, 4)
  testthat::expect_identical(normalized_fixed$selected_day, 4)
  testthat::expect_true(is.na(normalized_fixed$closest_day_diagnostic))
})

testthat::test_that("legacy run collectors normalize missing segment and endpoint fields", {
  env <- .load_invitro_structure_api()
  legacy_run <- list(
    grid_pre = c(22, 44),
    segment_results = list(list(
      segment = list(
        segment_id = "20.5_20.5_2",
        parent_segment_id = "20.5_20.5",
        cohort = "2N",
        oxygen_pct = 2,
        duration_days = 4,
        initial_cells = 10,
        final_cells = 18,
        obs_days_local = 0:4,
        passage_index = 4L,
        depth = 2L,
        data_ids = c("legacy_O1", "legacy_O2")
      ),
      sim = list(
        Ntot_live_obs = c(10, 12, 14, 16, 18),
        Ntot_dead_hypoxia_obs = c(0, 0, 1, 1, 2),
        Ntot_dead_buffer_obs = c(0, 1, 1, 2, 2),
        live_state_obs = rbind(
          c(5, 5),
          c(4, 8),
          c(3.5, 10.5),
          c(6, 10),
          c(9, 9)
        )
      ),
      selection = list(
        selected_index = 3L,
        selected_day = 2,
        selected_live_cells = 14,
        target_live_cells = 15,
        selected_frac = c(0.25, 0.75),
        predicted_mean_kary_N = 38.5
      )
    ))
  )

  daily <- env$ivt_collect_daily_counts(legacy_run)
  distribution <- env$ivt_collect_distribution_summary(legacy_run)
  quantiles <- env$ivt_collect_distribution_quantiles(
    legacy_run,
    probs = c(0.25, 0.5, 0.75)
  )

  testthat::expect_identical(nrow(daily), 5L)
  testthat::expect_identical(nrow(distribution), 2L)
  testthat::expect_identical(nrow(quantiles), 3L)
  testthat::expect_true(all(daily$lineage_id == "deprived"))
  testthat::expect_true(all(daily$scenario_id == "2N-deprived"))
  testthat::expect_true(all(daily$passage_id == "legacy_O1;legacy_O2"))
  testthat::expect_true(all(daily$lineage_passage_index == 3L))
  testthat::expect_true(all(daily$passage_duration == 4))
  testthat::expect_true(all(daily$endpoint_day == 2))
  testthat::expect_true(all(daily$selected_day == 2))
  testthat::expect_true(all(is.na(daily$closest_day_diagnostic)))
  testthat::expect_true(all(distribution$endpoint_day == 2))
  testthat::expect_equal(distribution$fraction, c(0.25, 0.75))
  testthat::expect_true(all(quantiles$selected_day == 2))

  normalize_segment <- get(
    ".ivt_normalize_segment_result",
    envir = env,
    inherits = FALSE
  )
  normalized <- normalize_segment(
    legacy_run$segment_results[[1L]],
    grid_pre = legacy_run$grid_pre
  )
  testthat::expect_identical(normalized$selection$selected_index, 3L)
  testthat::expect_identical(normalized$selection$selected_live_cells, 14)
  testthat::expect_equal(normalized$selection$selected_frac, c(0.25, 0.75))
  testthat::expect_identical(normalized$selection$predicted_mean_kary_N, 38.5)
  testthat::expect_true(is.na(normalized$selection$closest_index_diagnostic))
  testthat::expect_true(is.na(normalized$selection$closest_live_cells_diagnostic))

  no_endpoint_trajectory <- legacy_run$segment_results[[1L]]
  no_endpoint_trajectory$sim$live_state_obs <- NULL
  fallback <- normalize_segment(
    no_endpoint_trajectory,
    grid_pre = legacy_run$grid_pre
  )
  testthat::expect_identical(fallback$selection$selected_day, 2)
  testthat::expect_true(is.na(fallback$selection$closest_day_diagnostic))
  testthat::expect_equal(fallback$selection$selected_frac, c(0.25, 0.75))

  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  sys.source(
    file.path(
      workflow_root,
      "util",
      "o2_supply_demand_map_invitro_postfit_io_utils.R"
    ),
    envir = env,
    chdir = TRUE
  )
  sys.source(
    file.path(
      workflow_root,
      "simulation",
      "invitro",
      "population",
      "invitro_population_simulation_utils.R"
    ),
    envir = env,
    chdir = TRUE
  )
  population_tables <- env$ivt_sim_collect_population_tables(list(
    summary = data.frame(
      segment_id = "20.5_20.5_2",
      cohort = "2N",
      lineage_label = "deprived",
      duration_days = 4,
      selected_day = 2,
      predicted_live_cells = 14,
      predicted_mean_kary_N = 38.5,
      stringsAsFactors = FALSE
    ),
    run_2N = legacy_run,
    run_4N = list(grid_pre = c(22, 44), segment_results = list())
  ))
  refreshed_summary <- population_tables$invitro_lineage_summary
  testthat::expect_identical(refreshed_summary$endpoint_day, 2)
  testthat::expect_identical(refreshed_summary$selected_day, 2)
  testthat::expect_true(is.na(refreshed_summary$closest_day_diagnostic))
  testthat::expect_identical(refreshed_summary$predicted_live_cells, 18)
  testthat::expect_identical(
    refreshed_summary$predicted_passage_live_cells,
    14
  )
  testthat::expect_identical(refreshed_summary$predicted_mean_kary_N, 38.5)
})

testthat::test_that("post-fit diagnostics use exact definitions and one shared schema", {
  env <- .load_invitro_structure_api()
  passage_id <- "synthetic_2N_C_A1_seed"
  segment <- list(
    segment_id = "2N-C-A1",
    parent_segment_id = NA_character_,
    cohort = "2N",
    lineage_id = "C",
    lineage_group = "control",
    lineage_label = "C",
    scenario_id = "2N-C",
    lineage_terminal_key = "2N-C-A1",
    passage_index = 1L,
    lineage_passage_index = 1L,
    oxygen_pct = 1,
    duration_days = 2,
    passage_duration = 2,
    endpoint_day = 2,
    initial_cells = 100,
    final_cells = 180,
    obs_days_local = 0:2,
    data_ids = passage_id,
    passage_id = passage_id
  )
  sim <- list(
    Ntot_live_obs = c(100, 150, 180),
    Ntot_dead_hypoxia_obs = c(0, 10, 5),
    Ntot_dead_buffer_obs = c(0, 5, 4),
    Ntot_dead_total_obs = c(0, 15, 9),
    Ntot_total_obs = c(100, 165, 189),
    Vmm3_live_obs = c(100, 150, 180),
    Vmm3_dead_hypoxia_obs = c(0, 10, 5),
    Vmm3_dead_buffer_obs = c(0, 5, 4),
    Vmm3_dead_total_obs = c(0, 15, 9),
    Vmm3_total_obs = c(100, 165, 189),
    O2_target_obs = rep(1, 3),
    O2_eff_obs = rep(1, 3),
    live_state_obs = matrix(c(100, 150, 180), ncol = 1),
    cumulative_gross_divisions_obs = c(0, 60, 150),
    cumulative_hypoxia_deaths_obs = c(0, 12, 20),
    cumulative_dead_buffer_inflow_obs = c(0, 5, 10),
    cumulative_nonlive_inflow_obs = c(0, 17, 30),
    cumulative_gross_divisions_terminal = 150,
    cumulative_hypoxia_deaths_terminal = 20,
    cumulative_dead_buffer_inflow_terminal = 10,
    cumulative_nonlive_inflow_terminal = 30,
    cumulative_dead_buffer_inflow_definition = paste(
      "missegregation-linked nonviable daughters plus",
      "grid boundary-routed loss"
    )
  )
  selection <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = NA_real_,
    grid_pre = 44,
    target_live_cells = 180,
    obs_days_local = 0:2
  )
  run <- list(
    grid_pre = 44,
    segment_results = list(list(
      segment = segment,
      sim = sim,
      selection = selection
    )),
    initial_observations = list(),
    landmark_observations = list()
  )
  fit_data <- stats::setNames(list(list(
    g = log(180 / 100) / 2,
    passage_duration = 2,
    initial_cells = 100,
    final_cells = 180,
    kary = numeric(),
    flow = NULL
  )), passage_id)
  summary_df <- env$ivt_collect_lineage_summary(run, fit_data)

  death_data <- data.frame(
    observation_id = "synthetic_death_A1",
    cohort = "2N",
    lineage_id = "C",
    scenario_id = "2N-C",
    model_passage_id = passage_id,
    model_segment_id = "2N-C-A1",
    likelihood_observation_day = 2,
    dead_count = 10,
    eligible_denominator = 200,
    observed_dead_fraction = 0.05,
    stringsAsFactors = FALSE
  )
  death_df <- env$ivt_death_loglik_df(
    run_2N = run,
    run_4N = list(grid_pre = 44, segment_results = list()),
    death_data = death_data,
    sigma_death_logit = 0.75,
    death_fraction_eps = 1e-4
  )
  testthat::expect_equal(death_df$predicted_dead_fraction, 9 / 189)
  testthat::expect_equal(
    death_df$loglik,
    stats::dnorm(stats::qlogis(0.05), stats::qlogis(9 / 189), 0.75, log = TRUE)
  )
  testthat::expect_equal(
    death_df$logit_residual,
    stats::qlogis(0.05) - stats::qlogis(9 / 189)
  )

  testthat::expect_equal(
    summary_df$observed_net_population_doublings,
    log2(1.8)
  )
  testthat::expect_equal(
    summary_df$predicted_net_population_doublings,
    log2(1.8)
  )
  testthat::expect_equal(
    summary_df$observed_minimum_division_events,
    80
  )
  testthat::expect_equal(
    summary_df$predicted_gross_division_events,
    150
  )
  testthat::expect_equal(
    summary_df$predicted_cumulative_hypoxia_deaths,
    20
  )
  testthat::expect_equal(
    summary_df$predicted_cumulative_dead_buffer_inflow,
    10
  )
  testthat::expect_equal(
    summary_df$predicted_cumulative_nonlive_inflow,
    30
  )
  testthat::expect_equal(
    summary_df$predicted_divisions_per_initial_cell,
    1.5
  )
  testthat::expect_equal(
    summary_df$predicted_nonlive_inflow_to_division_ratio,
    0.2
  )
  testthat::expect_equal(
    summary_df$predicted_gross_division_to_net_gain_ratio,
    150 / 80
  )
  testthat::expect_identical(
    summary_df$protocol_feasibility_status,
    "PASS"
  )
  testthat::expect_equal(
    summary_df$cumulative_experimental_time,
    2
  )
  testthat::expect_equal(summary_df$cumulative_gross_divisions, 150)

  components <- list(
    summary = summary_df,
    growth_df = summary_df,
    death_df = data.frame(),
    run_2N = run,
    run_4N = list(grid_pre = 44, segment_results = list())
  )
  tables <- env$ivt_collect_postfit_tables(components)
  expected_tables <- c(
    "invitro_lineage_summary",
    "invitro_passage_audit",
    "invitro_growth_loglik",
    "invitro_passage_time_loglik",
    "invitro_death_loglik",
    "invitro_daily_counts",
    "invitro_division_death_diagnostics",
    "invitro_protocol_feasibility",
    "invitro_threshold_crossing_diagnostics"
  )
  testthat::expect_setequal(names(tables), expected_tables)
  testthat::expect_true(all(c(
    "predicted_gross_division_events",
    "predicted_cumulative_hypoxia_deaths",
    "predicted_cumulative_dead_buffer_inflow"
  ) %in% names(tables$invitro_division_death_diagnostics)))
  testthat::expect_true(all(c(
    "threshold_target_cells",
    "threshold_target_source",
    "predicted_threshold_crossing_day"
  ) %in% names(tables$invitro_threshold_crossing_diagnostics)))

  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  standalone_text <- paste(readLines(file.path(
    workflow_root,
    "util",
    "o2_supply_demand_map_fit_invitro_backend.R"
  ), warn = FALSE), collapse = "\n")
  joint_text <- paste(readLines(file.path(
    workflow_root,
    "util",
    "o2_supply_demand_map_fit_joint_backend.R"
  ), warn = FALSE), collapse = "\n")
  testthat::expect_match(
    standalone_text,
    "ivt_collect_postfit_tables(best_comp)",
    fixed = TRUE
  )
  testthat::expect_match(
    joint_text,
    "ivt_collect_postfit_tables(",
    fixed = TRUE
  )
})

.load_invitro_viz_api <- function() {
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  env <- new.env(parent = globalenv())
  assign("commandArgs", function(trailingOnly = FALSE) character(), envir = env)
  sys.source(
    file.path(
      workflow_root,
      "vis",
      "viz_invitro_model_O2_supply_demand_MAP_results.R"
    ),
    envir = env,
    chdir = TRUE
  )
  env
}

testthat::test_that("initial karyotypes are excluded only from scenario composite rows", {
  env <- .load_invitro_viz_api()
  observed <- data.frame(
    segment_id = c("2N-INITIAL-1", "2N-C-A1"),
    cohort = "2N",
    lineage_id = c("INITIAL", "C"),
    lineage_group = c("initial", "control"),
    lineage_label = c("INITIAL", "C"),
    scenario_id = c("2N-INITIAL", "2N-C"),
    passage_index = c(0L, 1L),
    observed_kary_N = c(44, 45),
    stringsAsFactors = FALSE
  )
  scenario_rows <- env$.invitro_scenario_kary_rows(observed)

  testthat::expect_identical(nrow(observed), 2L)
  testthat::expect_identical(nrow(scenario_rows), 1L)
  testthat::expect_identical(scenario_rows$lineage_id, "C")
  testthat::expect_identical(observed$lineage_id[[1L]], "INITIAL")
})

testthat::test_that("objective diagnostics use weighted hierarchical components", {
  env <- .load_invitro_viz_api()
  summary_df <- data.frame(
    metric = c(
      "objective_total",
      "growth_loglik", "ploidy_loglik", "flow_loglik", "death_loglik",
      "growth_weight", "ploidy_weight", "flow_weight", "death_weight",
      "growth_loglik_sum", "ploidy_loglik_sum", "flow_loglik_sum", "death_loglik_sum"
    ),
    value = c(8.5, -2, -3, -4, -2, 2, 0.5, 0.25, 1, -20, -30, -40, -20),
    stringsAsFactors = FALSE
  )
  metrics <- env$.invitro_objective_component_metrics(summary_df)

  testthat::expect_identical(attr(metrics, "scale_mode"), "hierarchical_objective")
  testthat::expect_equal(metrics$value[-1L], c(4, 1.5, 1, 2))
  testthat::expect_equal(sum(metrics$value[-1L]), metrics$value[[1L]])
  testthat::expect_false(any(abs(metrics$value) %in% c(20, 30, 40)))

  legacy_summary <- summary_df[grepl("_sum$", summary_df$metric), , drop = FALSE]
  legacy_metrics <- env$.invitro_objective_component_metrics(legacy_summary)
  testthat::expect_identical(attr(legacy_metrics, "scale_mode"), "legacy_raw_sum")
  testthat::expect_false("Total objective" %in% legacy_metrics$component)
  testthat::expect_match(attr(legacy_metrics, "y_label"), "not objective scale", fixed = TRUE)
})

testthat::test_that("objective diagnostics include passage time and buffer-prior penalties", {
  env <- .load_invitro_viz_api()
  summary_df <- data.frame(
    metric = c(
      "objective_total",
      "growth_loglik", "ploidy_loglik", "flow_loglik", "death_loglik",
      "passage_time_loglik",
      "growth_weight", "ploidy_weight", "flow_weight", "death_weight",
      "passage_time_weight", "buffer_prior_penalty"
    ),
    value = c(8.7, -2, -3, -4, -2, -1, 1, 0.5, 0.25, 1, 0.5, 1.7),
    stringsAsFactors = FALSE
  )
  metrics <- env$.invitro_objective_component_metrics(summary_df)

  testthat::expect_identical(attr(metrics, "scale_mode"), "hierarchical_objective")
  testthat::expect_true("Buffer soft-prior penalty" %in% metrics$component)
  testthat::expect_true(any(grepl("Passage time", metrics$component, fixed = TRUE)))
  testthat::expect_equal(sum(metrics$value[-1L]), metrics$value[[1L]])
})
