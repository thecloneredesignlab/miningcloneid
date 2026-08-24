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
    default_fit_objects <- standalone$ivt_load_fit_objects_compat(
      fit_objects_dir = fit_objects_dir,
      flow_density_path = flow_density_path
    )
    run_params <- standalone$ivt_load_default_run_params(cfg)
    # The frozen v2 boundary forbids upscaling. Use a deliberately growth-rich
    # fixture so every observed inoculum is reachable while these tests isolate
    # the Death-likelihood contract.
    run_params$lam_max <- 5
    default_comp <- standalone$ivt_objective_components(
      run_params = run_params,
      fit_objects = default_fit_objects,
      cfg = cfg,
      fallback_max_passage_days = 14,
      growth_weight = 1,
      ploidy_weight = 1,
      flow_weight = 1
    )
    load_death_data <- get(
      "load_death_data",
      envir = environment(standalone$ivt_load_fit_objects),
      inherits = FALSE
    )
    death_input <- load_death_data(file.path(repo_info$root, "oxygen"))
    fit_objects <- default_fit_objects
    fit_objects$death_enabled <- TRUE
    fit_objects$death_data <- death_input$data
    fit_objects$death_data_path <- death_input$path
    fit_objects$death_data_md5 <- death_input$md5
    fit_objects$death_data_n_file_rows <- death_input$n_file_rows
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
      default_fit_objects = default_fit_objects,
      default_comp = default_comp,
      death_input = death_input,
      fit_objects = fit_objects,
      run_params = run_params,
      comp = comp
    )
    cached
  }
})

testthat::test_that("production loading leaves Death disabled and does not read its file", {
  fixture <- .invitro_death_fixture()
  fit_objects <- fixture$default_fit_objects
  comp <- fixture$default_comp

  testthat::expect_false(isTRUE(fit_objects$death_enabled))
  testthat::expect_true(is.data.frame(fit_objects$death_data))
  testthat::expect_identical(nrow(fit_objects$death_data), 0L)
  testthat::expect_true(is.na(fit_objects$death_data_path))
  testthat::expect_true(is.na(fit_objects$death_data_md5))
  testthat::expect_identical(fit_objects$death_data_n_file_rows, 0L)
  testthat::expect_identical(nrow(comp$death_df), 0L)
  testthat::expect_true(is.na(comp$death_loglik))
  testthat::expect_true(is.na(comp$death_loglik_sum))
  testthat::expect_identical(comp$death_weight, 0.0)
  testthat::expect_identical(comp$n_death_observations, 0L)
  testthat::expect_true(is.finite(comp$total_loglik))
  testthat::expect_true(is.finite(comp$objective))
  testthat::expect_equal(
    comp$total_loglik,
    comp$growth_loglik + comp$ploidy_loglik + comp$flow_loglik,
    tolerance = 1e-15
  )

  fake_root <- file.path(tempdir(), "invitro_death_default_disabled")
  fake_death_dir <- file.path(fake_root, "data", "InVitroData")
  dir.create(fake_death_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("this file must not be read", file.path(fake_death_dir, basename(fixture$paths$death_data)))
  loaded_with_invalid_death_file <- fixture$standalone$ivt_load_fit_objects(
    repo_root = fake_root,
    fit_objects_dir = fixture$fit_objects_dir,
    flow_csv_path = fixture$flow_density_path
  )
  testthat::expect_false(isTRUE(loaded_with_invalid_death_file$death_enabled))
  testthat::expect_identical(nrow(loaded_with_invalid_death_file$death_data), 0L)
})

testthat::test_that("explicit Death loading enforces the fixed 90-row input contract", {
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

testthat::test_that("Death rows map to 45 old segments and branch-specific selected stocks", {
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
  collect_selected_stock <- function(run) {
    dplyr::bind_rows(lapply(run$segment_results, function(seg_res) {
      idx <- as.integer(seg_res$selection$selected_index)
      data.frame(
        model_passage_id = as.character(seg_res$segment$data_ids),
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
  selected <- selected_stock[
    match(death$model_passage_id, selected_stock$model_passage_id),
    ,
    drop = FALSE
  ]
  testthat::expect_false(anyNA(match(
    death$model_passage_id,
    selected_stock$model_passage_id
  )))
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

testthat::test_that("empty Death data contributes nothing even when enabled is requested", {
  fixture <- .invitro_death_fixture()
  empty_fit_objects <- fixture$default_fit_objects
  empty_fit_objects$death_enabled <- TRUE
  empty_comp <- fixture$standalone$ivt_objective_components(
    run_params = fixture$run_params,
    fit_objects = empty_fit_objects,
    cfg = fixture$cfg,
    fallback_max_passage_days = 14,
    growth_weight = 1,
    ploidy_weight = 1,
    flow_weight = 1
  )

  death_env <- environment(fixture$standalone$ivt_objective_components)
  empty_death <- get("ivt_empty_death_loglik_df", envir = death_env, inherits = FALSE)()
  testthat::expect_identical(nrow(empty_death), 0L)
  testthat::expect_identical(nrow(empty_comp$death_df), 0L)
  testthat::expect_true(is.na(empty_comp$death_loglik))
  testthat::expect_true(is.na(empty_comp$death_loglik_sum))
  testthat::expect_identical(empty_comp$death_weight, 0.0)
  testthat::expect_true(is.finite(empty_comp$total_loglik))
  testthat::expect_true(is.finite(empty_comp$objective))
  testthat::expect_equal(
    empty_comp$total_loglik,
    empty_comp$growth_loglik + empty_comp$ploidy_loglik + empty_comp$flow_loglik,
    tolerance = 1e-15
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
  joint_run_params$lam_max <- 5
  joint_default_comp <- joint_api$ivt_objective_components(
    run_params = joint_run_params,
    fit_objects = joint_fit_objects,
    cfg = joint_cfg,
    fallback_max_passage_days = 14,
    growth_weight = 1,
    ploidy_weight = 1,
    flow_weight = 1
  )
  testthat::expect_true(is.na(joint_default_comp$death_loglik))
  testthat::expect_true(is.na(joint_default_comp$death_loglik_sum))
  testthat::expect_identical(joint_default_comp$death_weight, 0.0)
  testthat::expect_true(is.finite(joint_default_comp$objective))
  testthat::expect_equal(
    joint_default_comp$total_loglik,
    joint_default_comp$growth_loglik +
      joint_default_comp$ploidy_loglik +
      joint_default_comp$flow_loglik,
    tolerance = 1e-15
  )

  joint_fit_objects$death_enabled <- TRUE
  joint_fit_objects$death_data <- fixture$death_input$data
  joint_fit_objects$death_data_path <- fixture$death_input$path
  joint_fit_objects$death_data_md5 <- fixture$death_input$md5
  joint_fit_objects$death_data_n_file_rows <- fixture$death_input$n_file_rows
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

testthat::test_that("standalone and joint in-vitro contexts use fixed v2 passage", {
  fixture <- .invitro_death_fixture()
  standalone <- .load_invitro_death_backend(fixture$paths$standalone_backend)
  standalone_cfg <- standalone$build_invitro_cfg(
    parameter_table = fixture$parameter_table,
    dt = 0.05,
    init_total_size = 1e6,
    o2_upper_bound = 21,
    fixed_oxygen = TRUE
  )
  testthat::expect_false("passage_mode" %in% names(standalone_cfg))
  testthat::expect_identical(standalone$INVITRO_PASSAGE_IMPLEMENTATION, "v2")

  joint <- .load_invitro_death_backend(fixture$paths$joint_backend)
  joint_ctx <- joint$build_joint_invitro_context(list(
    invitro_parameter_table = fixture$parameter_table,
    fit_objects_dir = fixture$fit_objects_dir,
    flow_density_path = fixture$flow_density_path
  ))
  testthat::expect_false("passage_mode" %in% names(joint_ctx$cfg))
  testthat::expect_identical(joint$INVITRO_ENV$INVITRO_PASSAGE_IMPLEMENTATION, "v2")
  testthat::expect_error(
    joint$build_joint_invitro_context(list(passage_mode = "org")),
    "passage_mode has been removed from the joint fit configuration"
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
