.load_invitro_passage_mode_api <- function() {
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  support <- new.env(parent = globalenv())
  sys.source(
    file.path(workflow_root, "util", "o2_supply_demand_map_shared.R"),
    envir = support,
    chdir = TRUE
  )
  sys.source(
    file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R"),
    envir = support,
    chdir = TRUE
  )
  env <- new.env(parent = support)
  source(
    file.path(workflow_root, "util", "o2_supply_demand_map_invitro_utils.R"),
    local = env,
    chdir = TRUE
  )
  env
}

testthat::test_that("org freezes upward reseeding while v1 selects from above", {
  env <- .load_invitro_passage_mode_api()
  sim <- list(
    Ntot_live_obs = c(5, 8, 12),
    live_state_obs = matrix(c(5, 8, 12), ncol = 1)
  )

  org <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 8,
    obs_days_local = 0:2,
    passage_mode = "org"
  )
  v1 <- env$ivt_extract_passage_end_state(
    sim = sim,
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 8,
    obs_days_local = 0:2,
    passage_mode = "v1"
  )

  testthat::expect_identical(org$selected_day, 1)
  testthat::expect_identical(org$selected_live_cells, 8)
  testthat::expect_equal(sum(org$reseeded_state), 10)
  testthat::expect_gt(org$boundary_scale, 1)
  testthat::expect_identical(org$reseed_mode, "org_rescale_to_requested_inoculum")

  testthat::expect_identical(v1$selected_day, 2)
  testthat::expect_identical(v1$selected_live_cells, 12)
  testthat::expect_equal(sum(v1$reseeded_state), 10)
  testthat::expect_lte(v1$boundary_scale, 1)
  testthat::expect_identical(v1$reseed_mode, "downsample_to_observed_inoculum")
})

testthat::test_that("v1 leaves sufficient org selections unchanged", {
  env <- .load_invitro_passage_mode_api()
  sim <- list(
    Ntot_live_obs = c(5, 12, 15),
    live_state_obs = cbind(c(2, 7, 9), c(3, 5, 6))
  )
  args <- list(
    sim = sim,
    reseed_live_cells = 10,
    grid_pre = c(22, 44),
    target_live_cells = 12,
    obs_days_local = 0:2
  )
  org <- do.call(env$ivt_extract_passage_end_state, c(args, passage_mode = "org"))
  v1 <- do.call(env$ivt_extract_passage_end_state, c(args, passage_mode = "v1"))

  testthat::expect_identical(v1$selected_index, org$selected_index)
  testthat::expect_identical(v1$selected_day, org$selected_day)
  testthat::expect_equal(v1$selected_frac, org$selected_frac, tolerance = 0)
  testthat::expect_equal(v1$reseeded_state, org$reseeded_state, tolerance = 0)
})

testthat::test_that("v1 eligibility is based on the state that will be reseeded", {
  env <- .load_invitro_passage_mode_api()
  selection <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 12),
      live_state_obs = matrix(c(5, 9), ncol = 1)
    ),
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 5,
    obs_days_local = 0:1,
    passage_mode = "v1"
  )

  testthat::expect_false(selection$passage_executed)
  testthat::expect_identical(selection$available_cells, 9)
  testthat::expect_null(selection$reseeded_state)
})

testthat::test_that("v1 rejects an unreachable inoculum without creating cells", {
  env <- .load_invitro_passage_mode_api()
  selection <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 8, 9),
      live_state_obs = matrix(c(5, 8, 9), ncol = 1)
    ),
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 8,
    obs_days_local = 0:2,
    passage_mode = "v1"
  )

  testthat::expect_false(selection$passage_executed)
  testthat::expect_identical(
    selection$reseed_mode,
    "no_passage_threshold_not_reached"
  )
  testthat::expect_null(selection$reseeded_state)
  testthat::expect_true(is.na(selection$boundary_scale))
  testthat::expect_identical(selection$available_cells, 9)
  testthat::expect_identical(selection$required_cells, 10)
  testthat::expect_match(selection$passage_failure_reason, "required_cells=10")
})

.two_segment_passage_adapter <- function() {
  list(
    cohort = "2N",
    n_segments = 2L,
    segments = list(
      list(
        segment_id = "p1",
        parent_segment_id = NA_character_,
        cohort = "2N",
        oxygen_pct = 2,
        duration_days = 2,
        initial_cells = 5,
        final_cells = 8,
        obs_days_local = 0:2,
        passage_index = 1L,
        data_ids = "p1"
      ),
      list(
        segment_id = "p2",
        parent_segment_id = "p1",
        cohort = "2N",
        oxygen_pct = 2,
        duration_days = 1,
        initial_cells = 10,
        final_cells = 11,
        obs_days_local = 0:1,
        passage_index = 2L,
        data_ids = "p2"
      )
    )
  )
}

.install_passage_runner_stubs <- function(env, first_endpoint = 9) {
  calls <- character()
  initial_totals <- numeric()
  env$cell_volume_mm3_by_N <- function(grid_pre, run_params, cfg) {
    rep(1, length(grid_pre))
  }
  env$ivt_run_segment_fixed_o2 <- function(segment,
                                           cfg,
                                           run_params,
                                           model_core,
                                           vol_by_N,
                                           init_state_override = NULL,
                                           init_cells_override = NULL) {
    calls <<- c(calls, segment$segment_id)
    initial_totals <<- c(initial_totals, as.numeric(init_cells_override))
    cells <- if (identical(segment$segment_id, "p1")) {
      c(5, 8, first_endpoint)
    } else {
      c(as.numeric(init_cells_override), as.numeric(init_cells_override) + 1)
    }
    list(
      segment = segment,
      sim = list(
        Ntot_live_obs = cells,
        live_state_obs = matrix(cells, ncol = 1)
      )
    )
  }
  list(
    calls = function() calls,
    initial_totals = function() initial_totals
  )
}

testthat::test_that("org keeps the historical insufficient boundary behavior", {
  env <- .load_invitro_passage_mode_api()
  observed <- .install_passage_runner_stubs(env, first_endpoint = 9)
  run <- env$ivt_run_lineage(
    adapter = .two_segment_passage_adapter(),
    cfg = list(init_total_size = 5, passage_mode = "org"),
    run_params = list(),
    model_core = list(grid_pre = 44)
  )

  testthat::expect_identical(observed$calls(), c("p1", "p2"))
  testthat::expect_identical(observed$initial_totals(), c(5, 10))
  testthat::expect_length(run$segment_results, 2L)
})

testthat::test_that("v1 stops before a child when its inoculum is unreachable", {
  env <- .load_invitro_passage_mode_api()
  observed <- .install_passage_runner_stubs(env, first_endpoint = 9)

  testthat::expect_error(
    env$ivt_run_lineage(
      adapter = .two_segment_passage_adapter(),
      cfg = list(init_total_size = 5, passage_mode = "v1"),
      run_params = list(),
      model_core = list(grid_pre = 44)
    ),
    "protocol_infeasible:.*parent_segment=p1; child_segment=p2",
    class = "invitro_protocol_infeasible"
  )
  testthat::expect_identical(observed$calls(), "p1")
})

testthat::test_that("passage mode validation rejects unknown modes", {
  env <- .load_invitro_passage_mode_api()
  testthat::expect_error(
    env$ivt_extract_passage_end_state(
      sim = list(
        Ntot_live_obs = c(5, 10),
        live_state_obs = matrix(c(5, 10), ncol = 1)
      ),
      reseed_live_cells = 5,
      grid_pre = 44,
      passage_mode = "future"
    ),
    "passage_mode must be one of: org, v1"
  )
})

testthat::test_that("the unified runner validates and forwards passage_mode", {
  runner <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "runner",
    "run_o2_fit.sh"
  )
  text <- paste(readLines(runner, warn = FALSE), collapse = "\n")
  testthat::expect_identical(
    system2("bash", c("-n", runner), stdout = FALSE, stderr = FALSE),
    0L
  )
  testthat::expect_match(text, "DEFAULT_PASSAGE_MODE=\"org\"", fixed = TRUE)
  testthat::expect_match(text, "--passage_mode=*) PASSAGE_MODE=", fixed = TRUE)
  testthat::expect_gte(
    lengths(regmatches(text, gregexpr("--passage_mode=${PASSAGE_MODE}", text, fixed = TRUE))),
    2L
  )
})
