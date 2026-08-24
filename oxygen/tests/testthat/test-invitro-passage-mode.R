.load_invitro_passage_api <- function() {
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

.fixed_branch_adapter_fixture <- function() {
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

testthat::test_that("the fixed passage boundary never creates cells", {
  env <- .load_invitro_passage_api()
  selection <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 8, 12),
      live_state_obs = matrix(c(5, 8, 12), ncol = 1)
    ),
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 8,
    obs_days_local = 0:2
  )

  testthat::expect_identical(env$INVITRO_PASSAGE_IMPLEMENTATION, "v2")
  testthat::expect_identical(selection$selected_index, 3L)
  testthat::expect_identical(selection$selected_day, 2)
  testthat::expect_identical(selection$selected_live_cells, 12)
  testthat::expect_equal(sum(selection$reseeded_state), 10, tolerance = 0)
  testthat::expect_identical(selection$reseed_mode, "downsample_to_observed_inoculum")
  testthat::expect_equal(selection$boundary_scale, 10 / 12, tolerance = 0)
  testthat::expect_lte(selection$boundary_scale, 1)
  testthat::expect_false("passage_mode" %in% names(selection))
})

testthat::test_that("passage eligibility uses the state that will be reseeded", {
  env <- .load_invitro_passage_api()
  selection <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 12),
      live_state_obs = matrix(c(5, 9), ncol = 1)
    ),
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 5,
    obs_days_local = 0:1
  )

  testthat::expect_false(selection$passage_executed)
  testthat::expect_identical(selection$available_cells, 9)
  testthat::expect_null(selection$reseeded_state)
})

testthat::test_that("an unreachable inoculum is reported without creating cells", {
  env <- .load_invitro_passage_api()
  selection <- env$ivt_extract_passage_end_state(
    sim = list(
      Ntot_live_obs = c(5, 8, 9),
      live_state_obs = matrix(c(5, 8, 9), ncol = 1)
    ),
    reseed_live_cells = 10,
    grid_pre = 44,
    target_live_cells = 8,
    obs_days_local = 0:2
  )

  testthat::expect_false(selection$passage_executed)
  testthat::expect_identical(selection$reseed_mode, "no_passage_threshold_not_reached")
  testthat::expect_null(selection$reseeded_state)
  testthat::expect_true(is.na(selection$boundary_scale))
  testthat::expect_identical(selection$available_cells, 9)
  testthat::expect_identical(selection$required_cells, 10)
  testthat::expect_match(selection$passage_failure_reason, "required_cells=10")
})

testthat::test_that("the fixed adapter splits only paired O1 and O2 nodes", {
  env <- .load_invitro_passage_api()
  fixture <- .fixed_branch_adapter_fixture()
  adapter <- env$ivt_build_job_table_adapter(
    jobs = fixture$jobs,
    fit_data = fixture$fit_data,
    cohort = "2N"
  )

  testthat::expect_identical(adapter$n_segments, 6L)
  testthat::expect_identical(
    vapply(adapter$segments, `[[`, character(1), "segment_id"),
    c(
      "root", "control_a1", "low_a1::O1", "low_a1::O2",
      "low_a2::O1", "low_a2::O2"
    )
  )
  testthat::expect_identical(
    vapply(adapter$segments, function(x) as.character(x$parent_segment_id), character(1)),
    c(NA_character_, "root", "root", "root", "low_a1::O1", "low_a1::O2")
  )
  testthat::expect_identical(adapter$segments[[3]]$data_ids, "fixture_2N_O1_A1_seed")
  testthat::expect_identical(adapter$segments[[4]]$data_ids, "fixture_2N_O2_A1_seed")
  testthat::expect_identical(adapter$segments[[3]]$duration_days, 2)
  testthat::expect_identical(adapter$segments[[4]]$duration_days, 3)
  testthat::expect_identical(adapter$segments[[5]]$initial_cells, 12)
  testthat::expect_identical(adapter$segments[[6]]$initial_cells, 13)
  testthat::expect_setequal(
    unlist(lapply(adapter$segments, `[[`, "data_ids"), use.names = FALSE),
    unlist(fixture$jobs$data_ids, use.names = FALSE)
  )
})

testthat::test_that("the committed adapters are byte-equivalent to pre-removal v2", {
  env <- .load_invitro_passage_api()
  fit_dir <- file.path(
    repo_info$root,
    "oxygen",
    "ploidyOxygen",
    "data",
    "fit_objects"
  )
  fit_data <- readRDS(file.path(fit_dir, "fit_data.Rds"))
  expected <- list(
    `2N` = list(
      n_segments = 62L,
      sha256 = "2549af2337809c38208797fb7406432a6c4d9dbb10f424d783d5217e144a7b20"
    ),
    `4N` = list(
      n_segments = 60L,
      sha256 = "49b6b4f71b3db2e0744d0aa8d3dd3dd9fc6250863ef61aeccc05b61564c81338"
    )
  )

  for (cohort in names(expected)) {
    jobs <- readRDS(file.path(fit_dir, paste0("jobs_", cohort, ".Rds")))
    adapter <- env$ivt_build_job_table_adapter(jobs, fit_data, cohort)
    testthat::expect_identical(adapter$n_segments, expected[[cohort]]$n_segments)
    testthat::expect_identical(
      digest::digest(adapter, algo = "sha256", serialize = TRUE),
      expected[[cohort]]$sha256,
      info = cohort
    )
    testthat::expect_setequal(
      unlist(lapply(adapter$segments, `[[`, "data_ids"), use.names = FALSE),
      unlist(jobs$data_ids, use.names = FALSE)
    )
    testthat::expect_false(anyDuplicated(
      unlist(lapply(adapter$segments, `[[`, "data_ids"), use.names = FALSE)
    ) > 0L, info = cohort)
  }
})

testthat::test_that("O1 and O2 states propagate only within their own branches", {
  env <- .load_invitro_passage_api()
  fixture <- .fixed_branch_adapter_fixture()
  adapter <- env$ivt_build_job_table_adapter(
    jobs = fixture$jobs,
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
    endpoint_state <- if (grepl("::O1$", segment$segment_id)) {
      c(endpoint_total, 0)
    } else if (grepl("::O2$", segment$segment_id)) {
      c(0, endpoint_total)
    } else {
      c(endpoint_total / 2, endpoint_total / 2)
    }
    n_obs <- length(segment$obs_days_local)
    state_mat <- t(vapply(
      seq_len(n_obs),
      function(i) if (i == n_obs) endpoint_state else initial_state,
      numeric(2)
    ))
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
    cfg = list(init_total_size = 100),
    run_params = list(),
    model_core = list(grid_pre = c(44, 88))
  )

  testthat::expect_equal(initial_states[["low_a2::O1"]], c(12, 0), tolerance = 0)
  testthat::expect_equal(initial_states[["low_a2::O2"]], c(0, 13), tolerance = 0)
})

testthat::test_that("the fixed adapter preserves the Death passage mapping", {
  env <- .load_invitro_passage_api()
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
    adapter <- env$ivt_build_job_table_adapter(jobs, fit_data, cohort)
    do.call(rbind, lapply(adapter$segments, function(segment) {
      data.frame(
        model_passage_id = as.character(segment$data_ids),
        cohort = cohort,
        matched_old_segment_id = sub("::O[12]$", "", as.character(segment$segment_id)),
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
  group_key <- paste(death$cohort, segment_map$matched_old_segment_id[hit], sep = "::")
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

testthat::test_that("in-vitro workers prefer their process-local C++ wrapper", {
  env <- .load_invitro_passage_api()
  symbol <- "cpp_o2simps_simulate_one"
  had_binding <- exists(symbol, envir = .GlobalEnv, inherits = FALSE)
  if (had_binding) old_binding <- get(symbol, envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_binding) {
      assign(symbol, old_binding, envir = .GlobalEnv)
    } else if (exists(symbol, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = symbol, envir = .GlobalEnv)
    }
  }, add = TRUE)

  worker_binding <- function(sim_args) "worker-local"
  assign(symbol, worker_binding, envir = .GlobalEnv)
  runner_env <- environment(env$ivt_run_segment_fixed_o2)
  testthat::expect_identical(
    get(".ivt_cpp_backend_function", envir = runner_env, inherits = FALSE)(symbol),
    worker_binding
  )
  testthat::expect_identical(
    get(".ivt_cpp_simulate_one", envir = runner_env, inherits = FALSE)(list()),
    "worker-local"
  )
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
  function() calls
}

testthat::test_that("the fixed lineage stops before an unreachable child", {
  env <- .load_invitro_passage_api()
  calls <- .install_passage_runner_stubs(env, first_endpoint = 9)

  condition <- tryCatch(
    {
      env$ivt_run_lineage(
        adapter = .two_segment_passage_adapter(),
        cfg = list(init_total_size = 5),
        run_params = list(),
        model_core = list(grid_pre = 44)
      )
      NULL
    },
    invitro_protocol_infeasible = identity
  )
  testthat::expect_s3_class(condition, "invitro_protocol_infeasible")
  testthat::expect_match(
    conditionMessage(condition),
    "protocol_infeasible:.*parent_segment=p1; child_segment=p2"
  )
  testthat::expect_identical(condition$segment$segment_id, "p1")
  testthat::expect_identical(condition$selection$required_cells, 10)
  testthat::expect_identical(condition$selection$available_cells, 9)
  testthat::expect_identical(calls(), "p1")
})

testthat::test_that("protocol-infeasible fits receive a graded DE penalty", {
  backend_path <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "util",
    "o2_supply_demand_map_fit_invitro_backend.R"
  )
  backend <- new.env(parent = globalenv())
  assign(
    "commandArgs",
    function(trailingOnly = FALSE) {
      if (isTRUE(trailingOnly)) character(0) else c("R", paste0("--file=", backend_path))
    },
    envir = backend
  )
  source(backend_path, local = backend, chdir = TRUE)

  make_condition <- function(available_cells, ordinal = 1L, cohort = "2N") {
    structure(
      list(
        message = "protocol_infeasible",
        call = NULL,
        segment = list(cohort = cohort),
        selection = list(required_cells = 100, available_cells = available_cells),
        segment_ordinal = ordinal,
        segment_count = 4L
      ),
      class = c("invitro_protocol_infeasible", "error", "condition")
    )
  }

  large_shortfall <- backend$invitro_protocol_penalty_objective(make_condition(10))
  small_shortfall <- backend$invitro_protocol_penalty_objective(make_condition(90))
  later_failure <- backend$invitro_protocol_penalty_objective(make_condition(90, ordinal = 3L))
  testthat::expect_lt(small_shortfall, large_shortfall)
  testthat::expect_lt(later_failure, small_shortfall)
  testthat::expect_lt(large_shortfall, backend$INVITRO_DE_PENALTY_OBJECTIVE)
  testthat::expect_identical(
    backend$invitro_protocol_penalty_objective(simpleError("other")),
    backend$INVITRO_DE_PENALTY_OBJECTIVE
  )
})

testthat::test_that("fixed-v2 DE initialization anchors parameter-table values", {
  backend_path <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "util",
    "o2_supply_demand_map_fit_invitro_backend.R"
  )
  backend <- new.env(parent = globalenv())
  assign(
    "commandArgs",
    function(trailingOnly = FALSE) {
      if (isTRUE(trailingOnly)) character(0) else c("R", paste0("--file=", backend_path))
    },
    envir = backend
  )
  source(backend_path, local = backend, chdir = TRUE)

  set.seed(10)
  population <- backend$build_invitro_de_initial_population(
    NP = 5,
    lower = c(a = -2, b = 0),
    upper = c(a = 2, b = 10),
    init = c(a = 0.5, b = 12)
  )
  testthat::expect_equal(population[1, ], c(0.5, 10), tolerance = 0)
  testthat::expect_true(all(population[, 1] >= -2 & population[, 1] <= 2))
  testthat::expect_true(all(population[, 2] >= 0 & population[, 2] <= 10))

  backend_text <- paste(readLines(backend_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(backend_text, "de_include_parameter_init <- TRUE", fixed = TRUE)
  testthat::expect_false(grepl("identical(passage_mode_use", backend_text, fixed = TRUE))
})

testthat::test_that("passage_mode input has been removed from config and launchers", {
  env <- .load_invitro_passage_api()
  testthat::expect_error(
    env$ivt_reject_removed_passage_mode("org", source = "test input"),
    "passage_mode has been removed from test input"
  )
  testthat::expect_invisible(env$ivt_reject_removed_passage_mode(NULL))

  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  runner <- file.path(workflow_root, "runner", "run_o2_fit.sh")
  submitters <- c(
    file.path(workflow_root, "hpc", "submit", "submit_o2_fit.sh"),
    file.path(workflow_root, "Docker", "hpc", "submit", "submit_o2_fit.sh")
  )
  workers <- c(
    file.path(workflow_root, "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub"),
    file.path(workflow_root, "Docker", "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub")
  )
  scripts <- c(runner, submitters, workers)
  testthat::expect_true(all(file.exists(scripts)))

  for (script in scripts) {
    text <- paste(readLines(script, warn = FALSE), collapse = "\n")
    testthat::expect_identical(
      system2("bash", c("-n", script), stdout = FALSE, stderr = FALSE),
      0L,
      info = script
    )
    testthat::expect_match(text, "passage_mode has been removed", ignore.case = TRUE, info = script)
    testthat::expect_false(grepl("DEFAULT_PASSAGE_MODE", text, fixed = TRUE), info = script)
    testthat::expect_false(grepl("org|v1|v2) ;;", text, fixed = TRUE), info = script)
    testthat::expect_false(grepl("--passage_mode=${PASSAGE_MODE}", text, fixed = TRUE), info = script)
    testthat::expect_false(grepl('--passage_mode="${PASSAGE_MODE}"', text, fixed = TRUE), info = script)
  }

  for (submitter in submitters) {
    text <- paste(readLines(submitter, warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl('export_arg+=",PASSAGE_MODE=', text, fixed = TRUE))
  }
  for (worker in workers) {
    text <- paste(readLines(worker, warn = FALSE), collapse = "\n")
    testthat::expect_match(text, "passage_implementation", fixed = TRUE)
  }

  config_text <- paste(readLines(
    file.path(repo_info$root, "oxygen", "config", "O2_supply_demand.yaml"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_false(grepl("passage_mode", config_text, fixed = TRUE))
})

testthat::test_that("container workers retain isolated provenance and Rcpp caches", {
  worker <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "Docker",
    "hpc",
    "array_workers",
    "submit_fit_seed_array_invitro_buffering.sub"
  )
  text <- paste(readLines(worker, warn = FALSE), collapse = "\n")
  testthat::expect_match(
    text,
    'TASK_PROVENANCE_DIR="${RUN_DIR}/.array_task_provenance/task${TASK_ID}"',
    fixed = TRUE
  )
  testthat::expect_false(grepl('o2sd_prov_write_standard "${RUN_DIR}"', text, fixed = TRUE))
  testthat::expect_match(text, 'O2SD_CONTAINER_RCPP_CACHE="/tmp/o2sd-rcpp-cache-', fixed = TRUE)
  testthat::expect_match(text, '${SLURM_ARRAY_JOB_ID:-manual}', fixed = TRUE)
  testthat::expect_match(text, '${SLURM_ARRAY_TASK_ID:-0}', fixed = TRUE)
})

testthat::test_that("all fitting entrypoints default DEoptim to 1000 iterations", {
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  runner <- file.path(workflow_root, "runner", "run_o2_fit.sh")
  joint_runner <- file.path(workflow_root, "runner", "run_multi_warmup_joint.sh")
  submitters <- c(
    file.path(workflow_root, "hpc", "submit", "submit_o2_fit.sh"),
    file.path(workflow_root, "Docker", "hpc", "submit", "submit_o2_fit.sh")
  )
  workers <- c(
    file.path(workflow_root, "hpc", "array_workers", "submit_fit_seed_array_buffering.sub"),
    file.path(workflow_root, "Docker", "hpc", "array_workers", "submit_fit_seed_array_buffering.sub"),
    file.path(workflow_root, "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub"),
    file.path(workflow_root, "Docker", "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub"),
    file.path(workflow_root, "hpc", "array_workers", "submit_fit_seed_array_joint_buffering.sub"),
    file.path(workflow_root, "Docker", "hpc", "array_workers", "submit_fit_seed_array_joint_buffering.sub")
  )
  scripts <- c(runner, joint_runner, submitters, workers)

  testthat::expect_true(all(file.exists(scripts)))
  for (script in scripts) {
    testthat::expect_identical(
      system2("bash", c("-n", script), stdout = FALSE, stderr = FALSE),
      0L,
      info = script
    )
    text <- paste(readLines(script, warn = FALSE), collapse = "\n")
    testthat::expect_match(text, 'DEFAULT_ITERMAX="1000"', fixed = TRUE, info = script)
  }

  runner_text <- paste(readLines(runner, warn = FALSE), collapse = "\n")
  joint_runner_text <- paste(readLines(joint_runner, warn = FALSE), collapse = "\n")
  testthat::expect_match(runner_text, 'DEFAULT_ITERMAX_MAX="1000"', fixed = TRUE)
  testthat::expect_match(runner_text, '--itermax_max=*) ITERMAX_MAX=', fixed = TRUE)
  testthat::expect_match(runner_text, '"--itermax_max=${ITERMAX_MAX}"', fixed = TRUE)
  combined_runner_text <- paste(runner_text, joint_runner_text, sep = "\n")
  testthat::expect_gte(
    lengths(regmatches(
      combined_runner_text,
      gregexpr('"--itermax=${ITERMAX}"', combined_runner_text, fixed = TRUE)
    )),
    3L
  )

  for (submitter in submitters) {
    text <- paste(readLines(submitter, warn = FALSE), collapse = "\n")
    testthat::expect_match(text, 'DEFAULT_ITERMAX_MAX="1000"', fixed = TRUE)
    testthat::expect_match(text, '--itermax_max=*) ITERMAX_MAX=', fixed = TRUE)
    testthat::expect_match(text, 'export_arg+=",ITERMAX_MAX=${ITERMAX_MAX}"', fixed = TRUE)
    testthat::expect_gte(
      lengths(regmatches(
        text,
        gregexpr('export_arg+=",ITERMAX=${ITERMAX}"', text, fixed = TRUE)
      )),
      2L
    )
    testthat::expect_match(text, '"--itermax=${ITERMAX}"', fixed = TRUE)
    testthat::expect_match(text, 'optimizer itermax_max "${ITERMAX_MAX:-NA}"', fixed = TRUE)
  }

  for (worker in workers) {
    text <- paste(readLines(worker, warn = FALSE), collapse = "\n")
    testthat::expect_match(text, '"--itermax=${ITERMAX}"', fixed = TRUE)
    testthat::expect_match(text, '--itermax="${ITERMAX}"', fixed = TRUE)
    testthat::expect_match(text, 'optimizer itermax "${ITERMAX}"', fixed = TRUE)
  }

  config_text <- paste(readLines(
    file.path(repo_info$root, "oxygen", "config", "O2_supply_demand.yaml"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_match(config_text, "itermax: 1000", fixed = TRUE)

  backend_defaults <- c(
    invivo = file.path(workflow_root, "util", "o2_supply_demand_map_fit_invivo_backend.R"),
    invitro = file.path(workflow_root, "util", "o2_supply_demand_map_fit_invitro_backend.R"),
    joint = file.path(workflow_root, "util", "o2_supply_demand_map_fit_joint_backend.R")
  )
  backend_text <- lapply(backend_defaults, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  })
  testthat::expect_match(
    backend_text$invivo,
    "itermax = as_int(argv$itermax, 1000L)",
    fixed = TRUE
  )
  testthat::expect_match(
    backend_text$invitro,
    ".first_non_null_local(argv$itermax, 1000L)",
    fixed = TRUE
  )
  testthat::expect_match(
    backend_text$invitro,
    ".first_non_null_local(argv$itermax_max, 1000L)",
    fixed = TRUE
  )
  testthat::expect_match(
    backend_text$joint,
    "itermax = as_int(cfg_raw$itermax, 1000L)",
    fixed = TRUE
  )

  multi_warmup_patterns <- c(
    "runner/run_multi_warmup_joint.sh" = 'ITERMAX="${ITERMAX:-${DEFAULT_ITERMAX}}"',
    "hpc/submit/submit_multi_warmup_joint.sh" = 'ITERMAX="${ITERMAX:-1000}"',
    "Docker/hpc/submit/submit_multi_warmup_joint.sh" = 'ITERMAX="${ITERMAX:-1000}"',
    "hpc/submit/build_multi_warmup_task_table.R" = 'itermax = setting(settings, "itermax", "1000")',
    "Docker/hpc/submit/build_multi_warmup_task_table.R" = 'itermax = setting(settings, "itermax", "1000")',
    "hpc/array_workers/run_multi_warmup_task_table_array.sub" = 'emit("ITERMAX", value("itermax", "1000"))',
    "Docker/hpc/array_workers/run_multi_warmup_task_table_array.sub" = 'emit("ITERMAX", value("itermax", "1000"))'
  )
  for (relative_path in names(multi_warmup_patterns)) {
    path <- file.path(workflow_root, relative_path)
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(
      text,
      multi_warmup_patterns[[relative_path]],
      fixed = TRUE,
      info = path
    )
  }
})
