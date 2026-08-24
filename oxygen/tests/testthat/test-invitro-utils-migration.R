.invitro_migration_paths <- function() {
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  list(
    oxygen_root = file.path(repo_info$root, "oxygen"),
    workflow_root = workflow_root,
    shared = file.path(workflow_root, "util", "o2_supply_demand_map_shared.R"),
    common = file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R"),
    loader = file.path(workflow_root, "util", "o2_supply_demand_map_invitro_utils.R"),
    plot = file.path(
      workflow_root,
      "vis",
      "invitro",
      "o2_supply_demand_map_invitro_plot_utils.R"
    ),
    legacy_dir = file.path(repo_info$root, "oxygen", "code", "in-vitro-utils")
  )
}

.load_canonical_invitro_api <- function(include_plots = TRUE) {
  paths <- .invitro_migration_paths()
  required <- c(paths$shared, paths$common, paths$loader)
  if (isTRUE(include_plots)) required <- c(required, paths$plot)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Missing canonical in-vitro file(s): ", paste(missing, collapse = ", "))
  }

  support_env <- new.env(parent = globalenv())
  sys.source(paths$shared, envir = support_env, chdir = TRUE)
  sys.source(paths$common, envir = support_env, chdir = TRUE)

  api_env <- new.env(parent = support_env)
  # Use the public loader exactly as production callers do. An absolute path is
  # intentional: source(chdir=TRUE) otherwise makes a relative ofile ambiguous.
  source(paths$loader, local = api_env, chdir = TRUE)
  if (isTRUE(include_plots)) {
    sys.source(paths$plot, envir = api_env, chdir = TRUE)
  }
  api_env
}

.expected_invitro_plot_api <- c(
  "ivt_ploidy_fraction_fill_scale",
  "ivt_plot_daily_counts",
  "ivt_plot_distribution_heatmap",
  "ivt_plot_lineage_counts",
  "ivt_plot_lineage_flow_density",
  "ivt_plot_lineage_growth",
  "ivt_plot_lineage_ploidy"
)

# Frozen from the stage-0 manifest at
# /tmp/oxygen_stage0_invitro/ivt_api_manifest.tsv. The test deliberately embeds
# the contract and never reads the /tmp baseline at runtime.
.expected_invitro_formals <- list(
  ivt_attach_g0g1_flow_data = as.pairlist(alist(
    fit_data = , repo_root = ,
    csv_path = file.path(repo_root, "data", "g0g1_ploidy_density_grid.csv")
  )),
  ivt_build_default_cfg = as.pairlist(alist(
    repo_root = , dt = 0.05, init_total_size = 1e+06,
    o2_upper_bound = 21, fixed_oxygen = TRUE
  )),
  ivt_build_job_table_adapter = as.pairlist(alist(
    jobs = , fit_data = , cohort = , fallback_max_passage_days = 14
  )),
  ivt_build_lineage_adapter = as.pairlist(alist(
    jobs = , fit_data = , terminal_sim_key = , cohort = ,
    max_segment_days = 14, obs_days_local = NULL
  )),
  ivt_build_segment_o2_profile = as.pairlist(alist(
    target_o2_pct = , cfg = , run_params = , burden_grid =
  )),
  ivt_canonical_flow_sample_name = as.pairlist(alist(x = )),
  ivt_choose_demo_terminal_key = as.pairlist(alist(
    jobs = , fit_data = , min_kary = 1L
  )),
  ivt_collect_daily_counts = as.pairlist(alist(run = )),
  ivt_collect_distribution_quantiles = as.pairlist(alist(run = , probs = )),
  ivt_collect_distribution_summary = as.pairlist(alist(run = )),
  ivt_collect_lineage_summary = as.pairlist(alist(run = , fit_data = )),
  ivt_collect_observed_flow_summary = as.pairlist(alist(run = , fit_data = )),
  ivt_collect_observed_kary_summary = as.pairlist(alist(run = , fit_data = )),
  ivt_default_best_param_path = as.pairlist(alist(repo_root = )),
  ivt_default_flow_kernel_sd_ploidy = as.pairlist(alist(run = , fit_data = )),
  ivt_extract_passage_end_state = as.pairlist(alist(
    sim = , reseed_live_cells = , grid_pre = , target_live_cells = NA_real_,
    obs_days_local = NULL
  )),
  ivt_flow_loglik_df = as.pairlist(alist(
    run = , fit_data = , n_unit = , sigma_flow_ploidy = ,
    density_floor = 1e-12
  )),
  ivt_flow_overlay_df = as.pairlist(alist(
    run = , fit_data = , n_unit = , sigma_flow_ploidy =
  )),
  ivt_gaussian_initial_state = as.pairlist(alist(
    grid_pre = , mean_N = , sd_N =
  )),
  ivt_growth_loglik_df = as.pairlist(alist(summary_df = , sigma_growth = )),
  ivt_list_terminal_keys = as.pairlist(alist(jobs = )),
  ivt_load_best_sampled_run_params = as.pairlist(alist(
    cfg = , repo_root = ,
    tsv_path = file.path(
      repo_root, "workflow", "sampled_objective_draws",
      "invitro_objective_array.tsv"
    ),
    best_param_path = ivt_default_best_param_path(repo_root)
  )),
  ivt_load_committed_best_run_params = as.pairlist(alist(
    cfg = , repo_root = ,
    best_param_path = ivt_default_best_param_path(repo_root)
  )),
  ivt_load_default_run_params = as.pairlist(alist(cfg = )),
  ivt_load_fit_objects = as.pairlist(alist(
    repo_root = ,
    fit_objects_dir = file.path(repo_root, "ploidyOxygen", "data", "fit_objects"),
    flow_csv_path = file.path(repo_root, "data", "g0g1_ploidy_density_grid.csv")
  )),
  ivt_load_run_params_from_row = as.pairlist(alist(best_row = , cfg = )),
  ivt_locate_repo_root = as.pairlist(alist(start = getwd())),
  ivt_log_growth_rate = as.pairlist(alist(
    initial_cells = , final_cells = , duration_days = , eps = 1e-12
  )),
  ivt_make_passage_map = as.pairlist(alist(adapter = )),
  ivt_nested_observed_median = as.pairlist(alist(
    observed_nested_list = , field = , default = NA_real_
  )),
  ivt_normalize_flow_slot = as.pairlist(alist(flow_entry = )),
  ivt_objective = as.pairlist(alist(
    run_params = , fit_objects = , cfg = , fallback_max_passage_days = 14,
    growth_weight = 1, ploidy_weight = 1, flow_weight = 1
  )),
  ivt_objective_components = as.pairlist(alist(
    run_params = , fit_objects = , cfg = , fallback_max_passage_days = 14,
    growth_weight = 1, ploidy_weight = 1, flow_weight = 1,
    ploidy_prob_floor = 1e-12, flow_density_floor = 1e-12,
    flow_kernel_sd_ploidy = NULL
  )),
  ivt_observed_passage_summary = as.pairlist(alist(fit_entry = )),
  ivt_optim_par_to_run_params = as.pairlist(alist(par_t = , cfg = )),
  ivt_optimizer_spec = as.pairlist(alist(cfg = )),
  ivt_parameter_table_path = as.pairlist(alist(repo_root = )),
  ivt_ploidy_fraction_fill_scale = as.pairlist(alist(
    fill_max = , name = "Fraction"
  )),
  ivt_ploidy_loglik_df = as.pairlist(alist(
    run = , fit_data = , sigma_kary = , prob_floor = 1e-12
  )),
  ivt_plot_daily_counts = as.pairlist(alist(count_df = )),
  ivt_plot_distribution_heatmap = as.pairlist(alist(dist_df = , max_N = 110)),
  ivt_plot_lineage_counts = as.pairlist(alist(summary_df = )),
  ivt_plot_lineage_flow_density = as.pairlist(alist(
    flow_overlay_df = , max_facets = 20L
  )),
  ivt_plot_lineage_growth = as.pairlist(alist(
    summary_df = , comparison_summary_df = NULL,
    primary_label = "Loaded best", comparison_label = "Comparison"
  )),
  ivt_plot_lineage_ploidy = as.pairlist(alist(
    summary_df = , quantile_df = NULL, observed_kary_df = NULL,
    primary_label = "Loaded best", quantile_alpha = 0.5
  )),
  ivt_predicted_flow_density = as.pairlist(alist(
    observed_grid = , run_grid_pre = , selected_frac = ,
    sigma_flow_ploidy = , n_unit =
  )),
  ivt_resolve_sampled_objective_tsv = as.pairlist(alist(
    repo_root = ,
    tsv_path = file.path(
      repo_root, "workflow", "sampled_objective_draws",
      "invitro_objective_array.tsv"
    )
  )),
  ivt_reject_removed_passage_mode = as.pairlist(alist(
    value = , source = "input"
  )),
  ivt_run_lineage = as.pairlist(alist(
    adapter = , cfg = , run_params = , model_core = NULL
  )),
  ivt_run_params_to_optim_par = as.pairlist(alist(run_params = , cfg = )),
  ivt_run_segment_fixed_o2 = as.pairlist(alist(
    segment = , cfg = , run_params = , model_core = , vol_by_N = ,
    init_state_override = NULL, init_cells_override = NULL
  )),
  ivt_segment_median_value = as.pairlist(alist(
    observed_list = , field = , default = NA_real_
  )),
  ivt_select_segment_observation = as.pairlist(alist(
    sim = , reseed_live_cells = , grid_pre = , target_live_cells = NA_real_,
    obs_days_local = NULL
  )),
  ivt_set_segment_o2 = as.pairlist(alist(
    target_o2_pct = , cfg = , run_params =
  )),
  ivt_smooth_kary_distribution = as.pairlist(alist(
    grid = , probs = , sigma_kary =
  )),
  ivt_source_map_model = as.pairlist(alist(repo_root = )),
  ivt_terminal_path_map = as.pairlist(alist(jobs = , cohort = )),
  ivt_trace_lineage = as.pairlist(alist(jobs = , terminal_sim_key = )),
  ivt_weighted_mean_kary_N = as.pairlist(alist(weights = , grid_pre = )),
  ivt_weighted_quantile_kary_N = as.pairlist(alist(
    weights = , grid_pre = , probs =
  ))
)

testthat::test_that("canonical in-vitro loader exports only its public API", {
  env <- .load_canonical_invitro_api(include_plots = FALSE)
  expected <- sort(setdiff(
    names(.expected_invitro_formals),
    .expected_invitro_plot_api
  ))
  exported <- ls(env, pattern = "^ivt_", all.names = TRUE)

  testthat::expect_length(exported, 53L)
  testthat::expect_identical(exported, expected)
  testthat::expect_identical(
    ls(env, all.names = TRUE),
    sort(c("INVITRO_PASSAGE_IMPLEMENTATION", expected))
  )
  testthat::expect_false(any(grepl("^\\.o2sd_invitro_", ls(env, all.names = TRUE))))
})

testthat::test_that("canonical loader and plot utils expose the fixed-v2 APIs and formals", {
  env <- .load_canonical_invitro_api(include_plots = TRUE)
  expected_names <- sort(names(.expected_invitro_formals))
  actual_names <- ls(env, pattern = "^ivt_", all.names = TRUE)

  testthat::expect_length(actual_names, 60L)
  testthat::expect_identical(actual_names, expected_names)
  for (fn_name in expected_names) {
    testthat::expect_true(
      is.function(get(fn_name, envir = env, inherits = FALSE)),
      info = fn_name
    )
    testthat::expect_identical(
      formals(get(fn_name, envir = env, inherits = FALSE)),
      .expected_invitro_formals[[fn_name]],
      info = fn_name
    )
  }
})

testthat::test_that("canonical loader resolves cwd and dispatcher call paths", {
  paths <- .invitro_migration_paths()
  support_env <- new.env(parent = globalenv())
  sys.source(paths$shared, envir = support_env, chdir = TRUE)
  sys.source(paths$common, envir = support_env, chdir = TRUE)

  relative_env <- new.env(parent = support_env)
  local({
    original_wd <- setwd(repo_info$root)
    on.exit(setwd(original_wd), add = TRUE)
    source(
      file.path(
        "oxygen",
        "code",
        "O2_supply_demand_MAP",
        "util",
        "o2_supply_demand_map_invitro_utils.R"
      ),
      local = relative_env,
      chdir = TRUE
    )
  })
  testthat::expect_length(ls(relative_env, pattern = "^ivt_"), 53L)

  tmp_env <- new.env(parent = support_env)
  local({
    original_wd <- setwd(tempdir())
    on.exit(setwd(original_wd), add = TRUE)
    source(paths$loader, local = tmp_env, chdir = TRUE)
  })
  testthat::expect_length(ls(tmp_env, pattern = "^ivt_"), 53L)

  dispatcher_path <- normalizePath(
    file.path(paths$workflow_root, "optimizer", "fit_model_O2_supply_demand_MAP.R"),
    mustWork = TRUE
  )
  dispatcher_expr <- paste0(
    "e <- new.env(parent = globalenv()); ",
    "source(", deparse(dispatcher_path), ", local = e, chdir = TRUE); ",
    "iv <- e$load_backend_env('fit_invitro'); ",
    "stopifnot(exists('ivt_objective_components', envir = iv, inherits = FALSE)); ",
    "joint <- e$load_backend_env('fit_joint'); ",
    "stopifnot(exists('INVITRO_ENV', envir = joint, inherits = FALSE)); ",
    "stopifnot(exists(",
    "'ivt_objective_components', envir = joint$INVITRO_ENV, inherits = FALSE",
    "));"
  )
  dispatcher_output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(dispatcher_expr)),
    stdout = TRUE,
    stderr = TRUE
  )
  dispatcher_status <- attr(dispatcher_output, "status")
  if (is.null(dispatcher_status)) dispatcher_status <- 0L
  testthat::expect_identical(
    as.integer(dispatcher_status),
    0L,
    info = paste(dispatcher_output, collapse = "\n")
  )
})

testthat::test_that("legacy in-vitro-utils directory has been removed", {
  paths <- .invitro_migration_paths()
  testthat::expect_false(dir.exists(paths$legacy_dir))
  testthat::expect_false(file.exists(file.path(paths$legacy_dir, "io.R")))
  testthat::expect_false(file.exists(file.path(paths$legacy_dir, "plotting.R")))
})

testthat::test_that("in-vitro defaults and canonical paths preserve stage-0 behavior", {
  env <- .load_canonical_invitro_api(include_plots = FALSE)
  oxygen_root <- normalizePath(.invitro_migration_paths()$oxygen_root, mustWork = TRUE)
  cfg_default <- env$ivt_build_default_cfg(repo_root = oxygen_root)

  testthat::expect_identical(
    cfg_default$parameter_table,
    file.path(oxygen_root, "data", "O2_supply_demand", "parameter_table_invitro.csv")
  )
  testthat::expect_identical(
    env$ivt_parameter_table_path(oxygen_root),
    file.path(
      oxygen_root,
      "data",
      "O2_supply_demand",
      "parameter_table_invitro_buffering.csv"
    )
  )
  testthat::expect_identical(
    env$ivt_default_best_param_path(oxygen_root),
    file.path(oxygen_root, "data", "O2_supply_demand", "invitro_best_parameters.tsv")
  )

  cfg_default$parameter_table <- env$ivt_parameter_table_path(oxygen_root)
  expected_cfg <- list(
    parameter_table = file.path(
      oxygen_root,
      "data",
      "O2_supply_demand",
      "parameter_table_invitro_buffering.csv"
    ),
    N_UNIT = 22L,
    N_MIN = 22L,
    N_MAX = 154L,
    start_with = "chr_number",
    DT = 0.05,
    init_total_size = 1e+06,
    o2_Nref = 1e+06,
    o2_min = 0,
    o2_burden_feedback = FALSE,
    O2_growth = TRUE,
    Crowding = TRUE,
    crowding = "logistic",
    K = 1e+12,
    dose_ref = 30,
    tx_mult_min = 0.05,
    min_pop = 1e-12,
    fit_treatment = FALSE,
    o2_cache_profile = FALSE,
    burden_log_eps = 1e-12,
    ploidy_O2_death = "ploidy_related",
    o2_S0_upper_bound = 21,
    tau_O2_init = 2,
    tau_O2 = 2,
    alpha_o2_init = 0.5,
    gamma_growth_init = 2,
    mu_hp_init = 0.001,
    gamma_mu_init = 1,
    o2_crit_init = 1,
    n_O_init = 1,
    k_clear_init = 0.001,
    p_wgd_init = 1e-04,
    harvest_init_multiplier = FALSE,
    use_necrosis_loss = FALSE,
    sigma_necrosis_logit = 0.75,
    lambda_necrosis = 1,
    necrosis_fraction_eps = 1e-04,
    prior_center_log_init_mult = 0,
    prior_sd_log_init_mult = 0.35,
    log_init_mult_lower = -1,
    log_init_mult_upper = 1,
    buffer_smax_init = 0.8,
    buffer_beta_init = 1,
    buffer_n_exp_init = 1,
    prior_center_buffer_smax = 0.8,
    prior_sd_buffer_smax = 0.25,
    prior_center_log10_buffer_beta = 0,
    prior_sd_log10_buffer_beta = 0.75,
    prior_center_log10_buffer_n_exp = 0,
    prior_sd_log10_buffer_n_exp = 0.75,
    dose_zero_only = TRUE,
    max_scenarios = Inf,
    truncate_at_treatment = FALSE,
    ploidy_at_harvest = TRUE
  )
  testthat::expect_identical(cfg_default, expected_cfg)

  feedback_cfg <- env$ivt_build_default_cfg(
    repo_root = oxygen_root,
    dt = 0.1,
    init_total_size = 2e+06,
    o2_upper_bound = 18,
    fixed_oxygen = FALSE
  )
  testthat::expect_identical(feedback_cfg$DT, 0.1)
  testthat::expect_identical(feedback_cfg$init_total_size, 2e+06)
  testthat::expect_identical(feedback_cfg$o2_Nref, 2e+06)
  testthat::expect_identical(feedback_cfg$o2_S0_upper_bound, 18)
  testthat::expect_true(feedback_cfg$o2_burden_feedback)
})

testthat::test_that("canonical fit-object IO preserves flow normalization", {
  env <- .load_canonical_invitro_api(include_plots = FALSE)
  temp_root <- tempfile("invitro-migration-io-")
  fit_dir <- file.path(temp_root, "fit_objects")
  dir.create(fit_dir, recursive = TRUE)
  on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)

  fit_id <- "experiment_SUM-159_NLS_A_rep1"
  fit_data <- setNames(
    list(list(flow = data.frame(bin = 1:2, count = c(4, 6)))),
    fit_id
  )
  jobs_2n <- data.frame(sim_key = "two_n", stringsAsFactors = FALSE)
  jobs_4n <- data.frame(sim_key = "four_n", stringsAsFactors = FALSE)
  saveRDS(fit_data, file.path(fit_dir, "fit_data.Rds"))
  saveRDS(jobs_2n, file.path(fit_dir, "jobs_2N.Rds"))
  saveRDS(jobs_4n, file.path(fit_dir, "jobs_4N.Rds"))

  flow_path <- file.path(temp_root, "flow.csv")
  utils::write.csv(
    data.frame(
      sample_name = c("Sample_SUM159_A", "Sample_SUM159_A"),
      grid_index = c(2L, 1L),
      ploidy = c(2.1, 1.9),
      log_density = log(c(0.2, 0.8)),
      stringsAsFactors = FALSE
    ),
    flow_path,
    row.names = FALSE
  )

  loaded <- env$ivt_load_fit_objects(
    repo_root = temp_root,
    fit_objects_dir = fit_dir,
    flow_csv_path = flow_path
  )
  testthat::expect_identical(loaded$jobs_2N, jobs_2n)
  testthat::expect_identical(loaded$jobs_4N, jobs_4n)
  flow <- loaded$fit_data[[fit_id]]$flow
  testthat::expect_identical(flow$type, "g0g1_ploidy_density")
  testthat::expect_identical(flow$sample_name_canonical, "SUM-159_NLS_A")
  testthat::expect_identical(flow$g0g1_ploidy_density$grid_index, 1:2)
  testthat::expect_equal(flow$g0g1_ploidy_density$density, c(0.8, 0.2))
  testthat::expect_identical(
    env$ivt_canonical_flow_sample_name(c("Sample_SUM159_A", "SUM159_B", "other")),
    c("SUM-159_NLS_A", "SUM-159_NLS_B", "other")
  )

  unchanged <- env$ivt_attach_g0g1_flow_data(
    fit_data = fit_data,
    repo_root = temp_root,
    csv_path = file.path(temp_root, "missing.csv")
  )
  testthat::expect_identical(unchanged, fit_data)
})

testthat::test_that("weighted summary helpers preserve edge-case contracts", {
  env <- .load_canonical_invitro_api(include_plots = FALSE)
  weights <- c(0, 1, 3)
  grid <- c(22, 44, 66)

  testthat::expect_equal(env$ivt_weighted_mean_kary_N(weights, grid), 60.5)
  testthat::expect_equal(
    env$ivt_weighted_quantile_kary_N(
      weights,
      grid,
      probs = c(-0.1, 0, 0.25, 0.5, 0.75, 1, 1.1)
    ),
    c(44, 44, 44, 66, 66, 66, 66)
  )
  testthat::expect_true(is.na(env$ivt_weighted_mean_kary_N(c(0, 0), c(22, 44))))
  testthat::expect_true(all(is.na(
    env$ivt_weighted_quantile_kary_N(c(0, NA), c(22, 44), c(0.25, 0.5))
  )))
  testthat::expect_length(
    env$ivt_weighted_quantile_kary_N(weights, grid, numeric(0)),
    0L
  )
  testthat::expect_error(
    env$ivt_weighted_mean_kary_N(1:2, grid),
    "length mismatch",
    fixed = TRUE
  )
  testthat::expect_equal(
    env$ivt_log_growth_rate(100, 400, 2),
    log(4) / 2
  )
  testthat::expect_true(is.na(env$ivt_log_growth_rate(0, 400, 2)))
})

testthat::test_that("optimizer parameters preserve stage-0 values and round-trip", {
  env <- .load_canonical_invitro_api(include_plots = FALSE)
  oxygen_root <- normalizePath(.invitro_migration_paths()$oxygen_root, mustWork = TRUE)
  cfg <- env$ivt_build_default_cfg(repo_root = oxygen_root)
  cfg$parameter_table <- env$ivt_parameter_table_path(oxygen_root)
  testthat::expect_true(file.exists(cfg$parameter_table))

  expected_optim <- c(
    log10_lam_max = -0.454452912438915,
    log10_p_misseg = -1.82390874094432,
    log10_k_o_mis = -3,
    buffer_smax = 0.8,
    log10_buffer_beta = 0,
    log10_buffer_n_exp = 0,
    log10_p_wgd = -5,
    log10_alpha_o2 = -0.301029995663981,
    gamma_growth = 1.5,
    log10_mu_hp = -2.30102999566398,
    gamma_mu = 2,
    log10_O2_crit = 0,
    n_O = 0.1,
    log10_p_mis_base = -5,
    log10_sigma_growth = -1,
    log10_sigma_kary = 0,
    init_mean_2N = 47.725,
    log10_init_sd_2N = 0.0538441249319371,
    init_mean_4N = 87.666667,
    log10_init_sd_4N = 0.562158985935193
  )
  spec <- env$ivt_optimizer_spec(cfg)
  testthat::expect_identical(spec$param_name, names(expected_optim))
  testthat::expect_equal(
    stats::setNames(spec$init, spec$param_name),
    expected_optim,
    tolerance = 1e-13
  )

  run_params <- env$ivt_load_default_run_params(cfg)
  optim_par <- env$ivt_run_params_to_optim_par(run_params, cfg)
  testthat::expect_equal(optim_par, expected_optim, tolerance = 1e-13)

  restored <- env$ivt_optim_par_to_run_params(optim_par, cfg)
  natural_names <- c(
    "lam_max", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
    "buffer_n_exp", "p_wgd", "alpha_o2", "gamma_growth", "mu_hp",
    "gamma_mu", "O2_crit", "n_O", "p_mis_base", "sigma_growth",
    "sigma_kary", "init_mean_2N", "init_sd_2N", "init_mean_4N",
    "init_sd_4N"
  )
  testthat::expect_equal(
    unlist(restored[natural_names], use.names = TRUE),
    unlist(run_params[natural_names], use.names = TRUE),
    tolerance = 1e-12
  )
})

testthat::test_that("seed10 objective and all exported tables match the fixed-v2 golden", {
  env <- .load_canonical_invitro_api(include_plots = FALSE)
  oxygen_root <- normalizePath(.invitro_migration_paths()$oxygen_root, mustWork = TRUE)
  parameter_table <- file.path(
    oxygen_root,
    "data",
    "O2_supply_demand",
    "parameter_table_invitro_buffering.csv"
  )
  fit_objects_dir <- file.path(oxygen_root, "ploidyOxygen", "data", "fit_objects")
  flow_density_path <- file.path(oxygen_root, "data", "g0g1_ploidy_density_grid.csv")

  cfg <- env$ivt_build_default_cfg(
    repo_root = oxygen_root,
    dt = 0.05,
    init_total_size = 1e6,
    o2_upper_bound = 21,
    fixed_oxygen = TRUE
  )
  cfg$parameter_table <- normalizePath(parameter_table, mustWork = TRUE)
  cfg <- get(
    "normalize_sim_cfg_common",
    envir = env,
    inherits = TRUE
  )(cfg, context = "fit")

  fit_objects <- env$ivt_load_fit_objects(
    repo_root = oxygen_root,
    fit_objects_dir = fit_objects_dir,
    flow_csv_path = flow_density_path
  )
  run_params <- env$ivt_load_default_run_params(cfg)
  seed10_best <- c(
    lam_max = 1.3485151936293,
    p_misseg = 0.0105535150770604,
    k_o_mis = 0.000112898682522252,
    buffer_smax = 1,
    buffer_beta = 1.59113027351929,
    buffer_n_exp = 3.15756611375965,
    p_wgd = 0.00033120539227592,
    o2_S0 = 3.5,
    kappa_O = 0.5,
    eta_o2 = 1.2,
    rho_2N = 56000,
    beta_size = 0.35,
    alpha_o2 = 0.5,
    gamma_growth = 1.5,
    mu_hp = 0.005,
    gamma_mu = 1.2,
    O2_crit = 1.58620678578554,
    n_O = 1.86857483366499,
    k_clear = 1e-04,
    sigma_burden = 0.2,
    sigma_growth = 0.203840002499712,
    sigma_kary = 0.0233545733387599,
    init_mean_2N = 46.3459118478185,
    init_sd_2N = 2.43406581031514,
    init_mean_4N = 83.8021620039257,
    init_sd_4N = 3.61719096780439,
    tau_O2 = 0.1,
    alpha = 1e-04,
    gamma = 1,
    p_mis_base = 8.54181703118544e-06,
    o2_S0_upper_bound = 21,
    o2_min = 0,
    o2_Nref = 1e+06
  )
  for (nm in names(seed10_best)) {
    run_params[[nm]] <- unname(seed10_best[[nm]])
  }
  run_params <- get(
    "normalize_run_params_common",
    envir = env,
    inherits = TRUE
  )(run_params, cfg = cfg)

  comp <- env$ivt_objective_components(
    run_params = run_params,
    fit_objects = fit_objects,
    cfg = cfg,
    fallback_max_passage_days = 14,
    growth_weight = 1,
    ploidy_weight = 1,
    flow_weight = 1
  )
  testthat::expect_equal(comp$objective, 3.8559429273636243, tolerance = 1e-12)
  testthat::expect_equal(comp$total_loglik, -3.8559429273636243, tolerance = 1e-12)
  testthat::expect_equal(comp$growth_loglik, 0.17068361334059332, tolerance = 1e-12)
  testthat::expect_equal(comp$ploidy_loglik, -3.0975000547409732, tolerance = 1e-12)
  testthat::expect_equal(comp$flow_loglik, -0.92912648596324454, tolerance = 1e-12)

  tables <- list(
    invitro_lineage_summary = comp$summary,
    invitro_growth_loglik = comp$growth_df,
    invitro_ploidy_loglik = comp$ploidy_df,
    invitro_flow_loglik = comp$flow_df,
    invitro_flow_overlay = comp$flow_overlay_df,
    invitro_distribution_summary = dplyr::bind_rows(
      env$ivt_collect_distribution_summary(comp$run_2N),
      env$ivt_collect_distribution_summary(comp$run_4N)
    ),
    invitro_distribution_quantiles = dplyr::bind_rows(
      env$ivt_collect_distribution_quantiles(
        comp$run_2N,
        probs = seq(0.01, 0.99, length.out = 50L)
      ),
      env$ivt_collect_distribution_quantiles(
        comp$run_4N,
        probs = seq(0.01, 0.99, length.out = 50L)
      )
    ),
    invitro_daily_counts = dplyr::bind_rows(
      env$ivt_collect_daily_counts(comp$run_2N),
      env$ivt_collect_daily_counts(comp$run_4N)
    ),
    invitro_observed_kary = dplyr::bind_rows(
      env$ivt_collect_observed_kary_summary(comp$run_2N, fit_objects$fit_data),
      env$ivt_collect_observed_kary_summary(comp$run_4N, fit_objects$fit_data)
    ),
    invitro_observed_flow = dplyr::bind_rows(
      env$ivt_collect_observed_flow_summary(comp$run_2N, fit_objects$fit_data),
      env$ivt_collect_observed_flow_summary(comp$run_4N, fit_objects$fit_data)
    )
  )
  expected_rows <- c(
    invitro_lineage_summary = 131L,
    invitro_growth_loglik = 114L,
    invitro_ploidy_loglik = 12L,
    invitro_flow_loglik = 20L,
    invitro_flow_overlay = 8000L,
    invitro_distribution_summary = 16226L,
    invitro_distribution_quantiles = 6100L,
    invitro_daily_counts = 782L,
    invitro_observed_kary = 220L,
    invitro_observed_flow = 4000L
  )
  testthat::expect_identical(
    vapply(tables, nrow, integer(1)),
    expected_rows
  )

  expected_md5 <- c(
    invitro_lineage_summary = "c2b2be0ac09f3d79b1d276357d7ea9f5",
    invitro_growth_loglik = "351f78dc1c9be9d9a96500d75ebb1807",
    invitro_ploidy_loglik = "1e3d0723874958648319d24f0087c6a2",
    invitro_flow_loglik = "3a5a6bcf51508f283deaa25e6a532b2a",
    invitro_flow_overlay = "5640b597ac3534d544f71bb521896287",
    invitro_distribution_summary = "e3c5c37d76a1e8583e0c5c13c2c27a8d",
    invitro_distribution_quantiles = "f06e1d21e3bb2b26e8140268494bf48b",
    invitro_daily_counts = "2078e9a5dbca2a09581292a9f8d3c7ed",
    invitro_observed_kary = "ebb47e0e1f4da49504d3112849af9478",
    invitro_observed_flow = "68d0915d8b1e7997d2793afe7748dce6"
  )
  output_dir <- tempfile("invitro-seed10-golden-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  output_paths <- vapply(
    names(tables),
    function(nm) {
      path <- file.path(output_dir, paste0(nm, ".tsv"))
      utils::write.table(
        tables[[nm]],
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = "NA"
      )
      path
    },
    character(1)
  )
  actual_md5 <- unname(tools::md5sum(output_paths))
  names(actual_md5) <- names(output_paths)
  testthat::expect_identical(actual_md5, expected_md5)
})
