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

.v2_branch_adapter_fixture <- function() {
  jobs <- data.frame(
    parent_key = c(NA_character_, "root", "root", "low_a1"),
    oxygen = c(20.5, 20.5, 2, 1),
    sim_key = c("root", "control_a1", "low_a1", "low_a2"),
    depth = 0:3,
    stringsAsFactors = FALSE
  )
  jobs$data_ids <- I(list(
    "fixture_2N_root_seed",
    "fixture_2N_C_A1_seed",
    c("fixture_2N_O1_A1_seed", "fixture_2N_O2_A1_seed"),
    c("fixture_2N_O1_A2_seed", "fixture_2N_O2_A2_seed")
  ))
  values <- list(
    fixture_2N_root_seed = c(1, 100, 200),
    fixture_2N_C_A1_seed = c(1, 100, 120),
    fixture_2N_O1_A1_seed = c(2, 10, 20),
    fixture_2N_O2_A1_seed = c(3, 11, 21),
    fixture_2N_O1_A2_seed = c(4, 12, 22),
    fixture_2N_O2_A2_seed = c(5, 13, 23)
  )
  fit_data <- lapply(values, function(x) {
    list(
      g = 0,
      passage_duration = x[[1]],
      initial_cells = x[[2]],
      final_cells = x[[3]],
      kary = numeric(),
      flow = NULL
    )
  })
  list(jobs = jobs, fit_data = fit_data)
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

testthat::test_that("v2 inherits the v1 no-upscaling boundary rule", {
  env <- .load_invitro_passage_mode_api()
  args <- list(
    sim = list(
      Ntot_live_obs = c(5, 8, 12),
      live_state_obs = matrix(c(5, 8, 12), ncol = 1)
    ),
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 8,
    obs_days_local = 0:2
  )
  v1 <- do.call(env$ivt_extract_passage_end_state, c(args, passage_mode = "v1"))
  v2 <- do.call(env$ivt_extract_passage_end_state, c(args, passage_mode = "v2"))

  testthat::expect_identical(v2$passage_mode, "v2")
  testthat::expect_identical(
    v2[setdiff(names(v2), "passage_mode")],
    v1[setdiff(names(v1), "passage_mode")]
  )
})

testthat::test_that("v2 splits only paired O1 and O2 adapter nodes", {
  env <- .load_invitro_passage_mode_api()
  fixture <- .v2_branch_adapter_fixture()

  org <- env$ivt_build_job_table_adapter(
    jobs = fixture$jobs,
    fit_data = fixture$fit_data,
    cohort = "2N"
  )
  v1_jobs <- fixture$jobs
  attr(v1_jobs, "passage_mode") <- "v1"
  v1 <- env$ivt_build_job_table_adapter(
    jobs = v1_jobs,
    fit_data = fixture$fit_data,
    cohort = "2N"
  )
  v2_jobs <- fixture$jobs
  attr(v2_jobs, "passage_mode") <- "v2"
  v2 <- env$ivt_build_job_table_adapter(
    jobs = v2_jobs,
    fit_data = fixture$fit_data,
    cohort = "2N"
  )

  testthat::expect_identical(v1, org)
  testthat::expect_identical(v2$n_segments, 6L)
  testthat::expect_identical(
    vapply(v2$segments, `[[`, character(1), "segment_id"),
    c(
      "root", "control_a1", "low_a1::O1", "low_a1::O2",
      "low_a2::O1", "low_a2::O2"
    )
  )
  testthat::expect_identical(v2$segments[[1]], org$segments[[1]])
  testthat::expect_identical(v2$segments[[2]], org$segments[[2]])
  testthat::expect_identical(
    vapply(v2$segments, function(x) as.character(x$parent_segment_id), character(1)),
    c(NA_character_, "root", "root", "root", "low_a1::O1", "low_a1::O2")
  )
  testthat::expect_identical(v2$segments[[3]]$data_ids, "fixture_2N_O1_A1_seed")
  testthat::expect_identical(v2$segments[[4]]$data_ids, "fixture_2N_O2_A1_seed")
  testthat::expect_identical(v2$segments[[3]]$duration_days, 2)
  testthat::expect_identical(v2$segments[[4]]$duration_days, 3)
  testthat::expect_identical(v2$segments[[5]]$initial_cells, 12)
  testthat::expect_identical(v2$segments[[6]]$initial_cells, 13)
  testthat::expect_setequal(
    unlist(lapply(v2$segments, `[[`, "data_ids"), use.names = FALSE),
    unlist(fixture$jobs$data_ids, use.names = FALSE)
  )
})

testthat::test_that("v2 propagates O1 and O2 states only within their own branches", {
  env <- .load_invitro_passage_mode_api()
  fixture <- .v2_branch_adapter_fixture()
  v2_jobs <- fixture$jobs
  attr(v2_jobs, "passage_mode") <- "v2"
  adapter <- env$ivt_build_job_table_adapter(
    jobs = v2_jobs,
    fit_data = fixture$fit_data,
    cohort = "2N"
  )

  initial_states <- list()
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
    initial_total <- as.numeric(init_cells_override)
    initial_state <- if (is.null(init_state_override)) {
      c(initial_total / 2, initial_total / 2)
    } else {
      as.numeric(init_state_override)
    }
    initial_states[[segment$segment_id]] <<- initial_state
    endpoint_total <- as.numeric(segment$final_cells)
    if (!is.finite(endpoint_total) || endpoint_total <= 0) {
      endpoint_total <- initial_total * 2
    }
    endpoint_state <- if (grepl("::O1$", segment$segment_id)) {
      c(endpoint_total, 0)
    } else if (grepl("::O2$", segment$segment_id)) {
      c(0, endpoint_total)
    } else {
      c(endpoint_total / 2, endpoint_total / 2)
    }
    n_obs <- length(segment$obs_days_local)
    state_mat <- vapply(
      seq_len(n_obs),
      function(i) if (i == n_obs) endpoint_state else initial_state,
      numeric(2)
    )
    state_mat <- t(state_mat)
    list(
      segment = segment,
      sim = list(
        Ntot_live_obs = rowSums(state_mat),
        live_state_obs = state_mat
      )
    )
  }

  env$ivt_run_lineage(
    adapter = adapter,
    cfg = list(init_total_size = 100, passage_mode = "v2"),
    run_params = list(),
    model_core = list(grid_pre = c(44, 88))
  )

  testthat::expect_equal(initial_states[["low_a2::O1"]], c(12, 0), tolerance = 0)
  testthat::expect_equal(initial_states[["low_a2::O2"]], c(0, 13), tolerance = 0)
})

testthat::test_that("v2 splits the committed O1 and O2 trees without touching other nodes", {
  env <- .load_invitro_passage_mode_api()
  fit_dir <- file.path(
    repo_info$root,
    "oxygen",
    "ploidyOxygen",
    "data",
    "fit_objects"
  )
  fit_data <- readRDS(file.path(fit_dir, "fit_data.Rds"))
  expected_segments <- c("2N" = 62L, "4N" = 60L)

  for (cohort in c("2N", "4N")) {
    jobs <- readRDS(file.path(fit_dir, paste0("jobs_", cohort, ".Rds")))
    org <- env$ivt_build_job_table_adapter(jobs, fit_data, cohort)
    v2_jobs <- jobs
    attr(v2_jobs, "passage_mode") <- "v2"
    v2 <- env$ivt_build_job_table_adapter(v2_jobs, fit_data, cohort)

    testthat::expect_identical(v2$n_segments, expected_segments[[cohort]])
    testthat::expect_setequal(
      unlist(lapply(v2$segments, `[[`, "data_ids"), use.names = FALSE),
      unlist(jobs$data_ids, use.names = FALSE)
    )
    testthat::expect_false(anyDuplicated(
      unlist(lapply(v2$segments, `[[`, "data_ids"), use.names = FALSE)
    ) > 0L)

    org_nonbranch <- Filter(
      function(x) !any(grepl("_O[12]_", x$data_ids)),
      org$segments
    )
    v2_nonbranch <- Filter(
      function(x) !any(grepl("_O[12]_", x$data_ids)),
      v2$segments
    )
    testthat::expect_identical(v2_nonbranch, org_nonbranch)

    branch_segments <- Filter(
      function(x) any(grepl("_O[12]_", x$data_ids)),
      v2$segments
    )
    testthat::expect_true(all(vapply(branch_segments, function(x) {
      length(unique(sub(".*_(O[12])_.*", "\\1", x$data_ids))) == 1L
    }, logical(1))))
    testthat::expect_true(all(vapply(branch_segments, function(x) {
      parent <- as.character(x$parent_segment_id)
      if (!grepl("::O[12]$", parent)) return(TRUE)
      sub(".*::", "", parent) == sub(".*::", "", x$segment_id)
    }, logical(1))))
  }
})

testthat::test_that("v2 preserves the 45-by-2 Death passage mapping contract", {
  env <- .load_invitro_passage_mode_api()
  fit_dir <- file.path(
    repo_info$root,
    "oxygen",
    "ploidyOxygen",
    "data",
    "fit_objects"
  )
  fit_data <- readRDS(file.path(fit_dir, "fit_data.Rds"))
  segment_maps <- lapply(c("2N", "4N"), function(cohort) {
    jobs <- readRDS(file.path(fit_dir, paste0("jobs_", cohort, ".Rds")))
    attr(jobs, "passage_mode") <- "v2"
    adapter <- env$ivt_build_job_table_adapter(jobs, fit_data, cohort)
    do.call(rbind, lapply(adapter$segments, function(segment) {
      data.frame(
        model_passage_id = as.character(segment$data_ids),
        cohort = cohort,
        matched_old_segment_id = sub(
          "::O[12]$",
          "",
          as.character(segment$segment_id)
        ),
        stringsAsFactors = FALSE
      )
    }))
  })
  segment_map <- do.call(rbind, segment_maps)
  death <- read.delim(
    file.path(
      repo_info$root,
      "data",
      "InVitroData",
      "sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  death <- death[death$include_in_current_endpoint_likelihood, , drop = FALSE]
  hit <- match(death$model_passage_id, segment_map$model_passage_id)

  testthat::expect_identical(nrow(death), 90L)
  testthat::expect_false(anyNA(hit))
  testthat::expect_identical(segment_map$cohort[hit], death$cohort)
  group_key <- paste(
    death$cohort,
    segment_map$matched_old_segment_id[hit],
    sep = "::"
  )
  group_counts <- table(group_key)
  testthat::expect_identical(length(group_counts), 45L)
  testthat::expect_true(all(group_counts == 2L))
  lineage_by_group <- split(death$lineage_id, group_key)
  testthat::expect_true(all(vapply(
    lineage_by_group,
    function(x) identical(sort(as.character(x)), c("O1", "O2")),
    logical(1)
  )))
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
    "passage_mode must be one of: org, v1, v2"
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
  testthat::expect_match(text, "--passage_mode=org|v1|v2", fixed = TRUE)
  testthat::expect_match(text, "org|v1|v2) ;;", fixed = TRUE)
  testthat::expect_match(text, "--passage_mode=*) PASSAGE_MODE=", fixed = TRUE)
  testthat::expect_gte(
    lengths(regmatches(text, gregexpr("--passage_mode=${PASSAGE_MODE}", text, fixed = TRUE))),
    2L
  )
})

testthat::test_that("HPC in-vitro submitters and array workers forward passage_mode", {
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  scripts <- c(
    file.path(workflow_root, "hpc", "submit", "submit_o2_fit.sh"),
    file.path(workflow_root, "Docker", "hpc", "submit", "submit_o2_fit.sh"),
    file.path(workflow_root, "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub"),
    file.path(workflow_root, "Docker", "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub")
  )
  testthat::expect_true(all(file.exists(scripts)))
  for (script in scripts) {
    text <- paste(readLines(script, warn = FALSE), collapse = "\n")
    testthat::expect_identical(
      system2("bash", c("-n", script), stdout = FALSE, stderr = FALSE),
      0L,
      info = script
    )
    testthat::expect_match(text, "PASSAGE_MODE", fixed = TRUE, info = script)
    testthat::expect_match(text, "org|v1|v2) ;;", fixed = TRUE, info = script)
  }
  submitter_text <- paste(readLines(scripts[[2]], warn = FALSE), collapse = "\n")
  worker_text <- paste(readLines(scripts[[4]], warn = FALSE), collapse = "\n")
  testthat::expect_match(submitter_text, 'export_arg+=",PASSAGE_MODE=', fixed = TRUE)
  testthat::expect_match(worker_text, '"--passage_mode=${PASSAGE_MODE}"', fixed = TRUE)
  testthat::expect_match(worker_text, '--passage_mode="${PASSAGE_MODE}"', fixed = TRUE)
  testthat::expect_match(
    worker_text,
    'O2SD_CONTAINER_RCPP_CACHE="${O2SD_CONTAINER_RCPP_CACHE:-/tmp/o2sd-rcpp-cache-',
    fixed = TRUE
  )
  testthat::expect_match(worker_text, '${SLURM_ARRAY_TASK_ID:-0}', fixed = TRUE)
})
