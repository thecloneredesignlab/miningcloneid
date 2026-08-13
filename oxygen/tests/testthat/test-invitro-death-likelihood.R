.invitro_death_paths <- function() {
  workflow <- file.path(repo_info$root, "oxygen", "code", "O2_supply_demand_MAP")
  list(
    workflow = workflow,
    standalone_backend = file.path(
      workflow,
      "util",
      "o2_supply_demand_map_fit_invitro_backend.R"
    ),
    joint_backend = file.path(
      workflow,
      "util",
      "o2_supply_demand_map_fit_joint_backend.R"
    ),
    death_data = file.path(
      repo_info$root,
      "data",
      "InVitroData",
      "sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"
    )
  )
}

.load_invitro_death_backend <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  had_global_command_args <- exists("commandArgs", envir = .GlobalEnv, inherits = FALSE)
  if (had_global_command_args) {
    old_global_command_args <- get("commandArgs", envir = .GlobalEnv, inherits = FALSE)
  }
  assign(
    "commandArgs",
    function(trailingOnly = FALSE) {
      if (isTRUE(trailingOnly)) character(0) else c("R", paste0("--file=", path))
    },
    envir = .GlobalEnv
  )
  on.exit({
    if (had_global_command_args) {
      assign("commandArgs", old_global_command_args, envir = .GlobalEnv)
    } else if (exists("commandArgs", envir = .GlobalEnv, inherits = FALSE)) {
      rm("commandArgs", envir = .GlobalEnv)
    }
  }, add = TRUE)
  source(path, local = env, chdir = TRUE)
  env
}

.invitro_death_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    paths <- .invitro_death_paths()
    standalone <- .load_invitro_death_backend(paths$standalone_backend)
    parameter_table <- standalone$default_parameter_table_path(must_exist = TRUE)
    fit_objects_dir <- standalone$default_fit_objects_dir(must_exist = TRUE)
    flow_density_path <- standalone$default_flow_density_path()
    cfg <- standalone$build_invitro_cfg(
      parameter_table = parameter_table,
      dt = 0.05,
      init_total_size = 1e6,
      o2_upper_bound = 21,
      fixed_oxygen = TRUE
    )
    fit_objects <- standalone$ivt_load_fit_objects_compat(
      fit_objects_dir = fit_objects_dir,
      flow_density_path = flow_density_path
    )
    run_params <- standalone$ivt_load_default_run_params(cfg)
    comp <- standalone$ivt_objective_components(
      run_params = run_params,
      fit_objects = fit_objects,
      cfg = cfg,
      fallback_max_passage_days = 14,
      growth_weight = 1,
      ploidy_weight = 1,
      flow_weight = 1
    )
    cached <<- list(
      paths = paths,
      standalone = standalone,
      parameter_table = parameter_table,
      fit_objects_dir = fit_objects_dir,
      flow_density_path = flow_density_path,
      cfg = cfg,
      fit_objects = fit_objects,
      run_params = run_params,
      comp = comp
    )
    cached
  }
})

testthat::test_that("Death loader enforces the fixed 90-row input contract", {
  fixture <- .invitro_death_fixture()
  fit_objects <- fixture$fit_objects
  death <- fit_objects$death_data

  testthat::expect_true(isTRUE(fit_objects$death_enabled))
  testthat::expect_identical(
    fit_objects$death_data_path,
    normalizePath(fixture$paths$death_data, mustWork = TRUE)
  )
  testthat::expect_identical(
    fit_objects$death_data_md5,
    "f346d7803596e15968b90b8213f15d95"
  )
  testthat::expect_identical(fit_objects$death_data_n_file_rows, 90L)
  testthat::expect_identical(nrow(death), 90L)
  testthat::expect_true(all(death$include_in_current_endpoint_likelihood))
  testthat::expect_false(anyDuplicated(death$observation_id) > 0L)
  testthat::expect_false(anyDuplicated(death$model_passage_id) > 0L)
  testthat::expect_true(all(death$dead_count >= 0))
  testthat::expect_true(all(death$dead_count <= death$eligible_denominator))
  testthat::expect_identical(
    death$observed_live_count,
    death$eligible_denominator - death$dead_count
  )
  testthat::expect_true(all(death$observed_live_count > 0))
  testthat::expect_lt(
    max(abs(death$observed_dead_fraction - death$dead_count / death$eligible_denominator)),
    5e-9
  )
  observed_groups <- table(paste(death$cohort, death$lineage_id, sep = "-"))
  expected_groups <- c("2N-O1" = 23L, "2N-O2" = 23L, "4N-O1" = 22L, "4N-O2" = 22L)
  testthat::expect_identical(
    as.integer(observed_groups[names(expected_groups)]),
    unname(expected_groups)
  )
})

testthat::test_that("Death rows map to 45 old segments and selected remaining stocks", {
  comp <- .invitro_death_fixture()$comp
  death <- comp$death_df
  testthat::expect_identical(nrow(death), 90L)
  testthat::expect_false(anyDuplicated(death$observation_id) > 0L)
  testthat::expect_true(all(is.finite(death$loglik)))

  segment_key <- paste(death$cohort, death$segment_id, sep = "::")
  segment_counts <- table(segment_key)
  testthat::expect_identical(length(segment_counts), 45L)
  testthat::expect_true(all(segment_counts == 2L))
  lineage_by_segment <- split(death$lineage_id, segment_key)
  testthat::expect_true(all(vapply(
    lineage_by_segment,
    function(x) identical(sort(x), c("O1", "O2")),
    logical(1)
  )))
  testthat::expect_true(all(vapply(
    split(death, segment_key),
    function(x) {
      length(unique(x$selected_index)) == 1L &&
        length(unique(x$selected_day)) == 1L &&
        length(unique(x$predicted_live_count)) == 1L &&
        length(unique(x$predicted_dead_count)) == 1L
    },
    logical(1)
  )))

  collect_selected_stock <- function(run) {
    dplyr::bind_rows(lapply(run$segment_results, function(seg_res) {
      idx <- as.integer(seg_res$selection$selected_index)
      data.frame(
        cohort = seg_res$segment$cohort,
        segment_id = seg_res$segment$segment_id,
        selected_index = idx,
        selected_day = as.numeric(seg_res$selection$selected_day),
        predicted_live_count = as.numeric(seg_res$sim$Ntot_live_obs[[idx]]),
        predicted_dead_hypoxia_count = as.numeric(seg_res$sim$Ntot_dead_hypoxia_obs[[idx]]),
        predicted_dead_buffer_count = as.numeric(seg_res$sim$Ntot_dead_buffer_obs[[idx]]),
        stringsAsFactors = FALSE
      )
    }))
  }
  selected_stock <- dplyr::bind_rows(
    collect_selected_stock(comp$run_2N),
    collect_selected_stock(comp$run_4N)
  )
  selected_key <- paste(selected_stock$cohort, selected_stock$segment_id, sep = "::")
  selected <- selected_stock[match(segment_key, selected_key), , drop = FALSE]
  testthat::expect_false(anyNA(match(segment_key, selected_key)))
  testthat::expect_equal(death$selected_index, selected$selected_index, tolerance = 0)
  testthat::expect_equal(death$selected_day, selected$selected_day, tolerance = 1e-15)
  testthat::expect_equal(death$predicted_live_count, selected$predicted_live_count, tolerance = 1e-12)
  testthat::expect_equal(
    death$predicted_dead_hypoxia_count,
    selected$predicted_dead_hypoxia_count,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    death$predicted_dead_buffer_count,
    selected$predicted_dead_buffer_count,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    death$predicted_dead_count,
    selected$predicted_dead_hypoxia_count + selected$predicted_dead_buffer_count,
    tolerance = 1e-12
  )
})

testthat::test_that("Death uses the reference logit Gaussian likelihood", {
  comp <- .invitro_death_fixture()$comp
  death <- comp$death_df
  eps <- 1e-4
  clamp_probability <- function(x) pmin(pmax(x, eps), 1 - eps)
  observed_fraction_from_counts <- death$observed_dead_count / death$eligible_denominator
  observed_fraction <- death$observed_dead_fraction
  predicted_fraction <- death$predicted_dead_count /
    (death$predicted_live_count + death$predicted_dead_count)
  observed_logit <- stats::qlogis(clamp_probability(observed_fraction))
  predicted_logit <- stats::qlogis(clamp_probability(predicted_fraction))
  row_loglik <- stats::dnorm(observed_logit, mean = predicted_logit, sd = 0.75, log = TRUE)

  testthat::expect_equal(observed_fraction, observed_fraction_from_counts, tolerance = 5e-9)
  testthat::expect_equal(death$predicted_dead_fraction, predicted_fraction, tolerance = 1e-15)
  testthat::expect_equal(death$observed_dead_logit, observed_logit, tolerance = 1e-12)
  testthat::expect_equal(death$predicted_dead_logit, predicted_logit, tolerance = 1e-12)
  testthat::expect_equal(death$logit_residual, observed_logit - predicted_logit, tolerance = 1e-12)
  testthat::expect_equal(death$loglik, row_loglik, tolerance = 1e-12)
  testthat::expect_identical(unique(death$sigma_death_logit), 0.75)
  testthat::expect_identical(unique(death$death_fraction_eps), 1e-4)
  testthat::expect_equal(comp$death_loglik_sum, sum(row_loglik), tolerance = 1e-12)
  testthat::expect_equal(comp$death_loglik, mean(row_loglik), tolerance = 1e-12)
  testthat::expect_equal(
    comp$total_loglik,
    comp$growth_loglik + comp$ploidy_loglik + comp$flow_loglik + comp$death_loglik,
    tolerance = 1e-15
  )
  testthat::expect_identical(comp$death_weight, 1.0)
  testthat::expect_identical(comp$n_death_observations, 90L)
})

testthat::test_that("production opt-in leaves the legacy direct objective contract intact", {
  fixture <- .invitro_death_fixture()
  direct_fit_objects <- fixture$standalone$ivt_load_fit_objects(
    repo_root = file.path(repo_info$root, "oxygen"),
    fit_objects_dir = fixture$fit_objects_dir,
    flow_csv_path = fixture$flow_density_path
  )
  testthat::expect_false(isTRUE(direct_fit_objects$death_enabled))
  testthat::expect_identical(nrow(direct_fit_objects$death_data), 90L)

  death_env <- environment(fixture$standalone$ivt_objective_components)
  empty_death <- get("ivt_empty_death_loglik_df", envir = death_env, inherits = FALSE)()
  testthat::expect_identical(nrow(empty_death), 0L)
  testthat::expect_identical(
    fixture$comp$total_loglik - fixture$comp$death_loglik,
    fixture$comp$growth_loglik + fixture$comp$ploidy_loglik + fixture$comp$flow_loglik
  )
})

testthat::test_that("standalone and joint backends produce identical Death likelihoods", {
  fixture <- .invitro_death_fixture()
  joint <- .load_invitro_death_backend(fixture$paths$joint_backend)
  joint_api <- joint$INVITRO_ENV
  joint_cfg <- joint_api$build_invitro_cfg(
    parameter_table = fixture$parameter_table,
    dt = 0.05,
    init_total_size = 1e6,
    o2_upper_bound = 21,
    fixed_oxygen = TRUE
  )
  joint_fit_objects <- joint_api$ivt_load_fit_objects_compat(
    fit_objects_dir = fixture$fit_objects_dir,
    flow_density_path = fixture$flow_density_path
  )
  joint_run_params <- joint_api$ivt_load_default_run_params(joint_cfg)
  joint_comp <- joint_api$ivt_objective_components(
    run_params = joint_run_params,
    fit_objects = joint_fit_objects,
    cfg = joint_cfg,
    fallback_max_passage_days = 14,
    growth_weight = 1,
    ploidy_weight = 1,
    flow_weight = 1
  )

  testthat::expect_equal(joint_comp$death_loglik, fixture$comp$death_loglik, tolerance = 1e-12)
  testthat::expect_equal(joint_comp$death_loglik_sum, fixture$comp$death_loglik_sum, tolerance = 1e-12)
  testthat::expect_equal(joint_comp$death_df$loglik, fixture$comp$death_df$loglik, tolerance = 1e-12)
  testthat::expect_identical(
    joint_comp$death_df$model_passage_id,
    fixture$comp$death_df$model_passage_id
  )
})

testthat::test_that("standalone and joint outputs enumerate only Death additions", {
  paths <- .invitro_death_paths()
  standalone_text <- paste(readLines(paths$standalone_backend, warn = FALSE), collapse = "\n")
  joint_text <- paste(readLines(paths$joint_backend, warn = FALSE), collapse = "\n")
  required_summary_fields <- c(
    "death_loglik", "death_loglik_sum", "sigma_death_logit",
    "death_fraction_eps", "death_weight", "n_death_observations",
    "death_data_path", "death_data_md5", "death_data_n_file_rows"
  )
  testthat::expect_match(standalone_text, "invitro_death_loglik.tsv", fixed = TRUE)
  testthat::expect_match(joint_text, "invitro_death_loglik.tsv", fixed = TRUE)
  for (field in required_summary_fields) {
    testthat::expect_match(standalone_text, paste0('"', field, '"'), fixed = TRUE)
    testthat::expect_match(joint_text, field, fixed = TRUE)
  }
})
