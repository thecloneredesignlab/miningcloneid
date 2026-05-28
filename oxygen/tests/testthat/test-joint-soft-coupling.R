testthat::skip_if_not_installed("DEoptim")
testthat::skip_if_not_installed("dplyr")

joint_backend_env <- local({
  env <- new.env(parent = globalenv())
  joint_path <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2G_supply_demand_MAP",
    "util",
    "o2g_supply_demand_map_fit_joint_backend.R"
  )
  had_command_args <- exists("commandArgs", envir = globalenv(), inherits = FALSE)
  old_command_args <- if (had_command_args) get("commandArgs", envir = globalenv(), inherits = FALSE) else NULL
  assign("commandArgs", function(trailingOnly = FALSE) character(0), envir = globalenv())
  on.exit({
    if (had_command_args) {
      assign("commandArgs", old_command_args, envir = globalenv())
    } else {
      rm("commandArgs", envir = globalenv())
    }
  }, add = TRUE)
  sys.source(joint_path, envir = env, chdir = TRUE)
  env
})

make_soft_ctx <- function(delta = 0, center = 0, sigma = 0.35) {
  meta <- data.frame(
    parameter = "O2_crit",
    center_name = "log10_O2_crit",
    invitro_name = "log10_O2_crit",
    delta_name = "delta__log10_O2_crit",
    transform = "log10",
    center_init_t = 0,
    center_lower_t = -1,
    center_upper_t = log10(2.5),
    invivo_lower_t = -1,
    invivo_upper_t = log10(2.5),
    invitro_lower_t = -4,
    invitro_upper_t = log10(2.5),
    delta_lower_t = -0.5,
    delta_upper_t = 0.5,
    sigma_delta = sigma,
    stringsAsFactors = FALSE
  )
  init <- c(
    log10_O2_crit = center,
    delta__log10_O2_crit = delta
  )
  list(
    init = init,
    invivo_names = "log10_O2_crit",
    joint_soft_coupling = list(
      enabled = TRUE,
      params = "O2_crit",
      metadata = meta
    )
  )
}

testthat::test_that("default soft-coupling list matches the joint split policy", {
  expected <- c(
    "O2_crit", "mu_hp", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O",
    "lam_max", "p_mis_base", "p_wgd", "gamma_mu"
  )
  actual <- joint_backend_env$joint_soft_split_natural_param_names(
    cfg_raw = list(joint_soft_coupling_enable = TRUE),
    invivo_glucose = FALSE
  )

  testthat::expect_identical(actual, expected)
  testthat::expect_false("alpha_o2" %in% actual)
  testthat::expect_false("gamma_growth" %in% actual)
})

testthat::test_that("soft coupling maps zero delta to the center", {
  ctx <- make_soft_ctx(delta = 0, center = 0)
  derived <- joint_backend_env$joint_build_context_specific_transformed_vectors(ctx$init, ctx)

  testthat::expect_equal(unname(derived$invivo_par[["log10_O2_crit"]]), 0, tolerance = 1e-12)
  testthat::expect_equal(derived$soft_derived$vivo_transformed, 0, tolerance = 1e-12)
  testthat::expect_equal(derived$soft_derived$vitro_transformed, 0, tolerance = 1e-12)
  testthat::expect_false(derived$soft_derived$vivo_clipped)
  testthat::expect_false(derived$soft_derived$vitro_clipped)
})

testthat::test_that("soft coupling is symmetric on the transformed scale", {
  ctx <- make_soft_ctx(delta = 0.2, center = -0.2)
  derived <- joint_backend_env$joint_build_context_specific_transformed_vectors(ctx$init, ctx)

  testthat::expect_equal(
    derived$soft_derived$vivo_transformed - derived$soft_derived$vitro_transformed,
    0.2,
    tolerance = 1e-12
  )
})

testthat::test_that("soft coupling penalty uses delta squared over two sigma squared", {
  ctx <- make_soft_ctx(delta = 0.21, center = 0, sigma = 0.35)
  penalty <- joint_backend_env$joint_soft_coupling_penalty_components(ctx$init, ctx)

  testthat::expect_equal(
    penalty$total,
    0.21^2 / (2 * 0.35^2),
    tolerance = 1e-12
  )
})

testthat::test_that("soft coupling reports backend-bound clipping", {
  ctx <- make_soft_ctx(delta = 0.5, center = log10(2.5))
  derived <- joint_backend_env$joint_build_context_specific_transformed_vectors(ctx$init, ctx)

  testthat::expect_true(derived$soft_derived$vivo_clipped)
  testthat::expect_equal(derived$soft_derived$boundary_status_vivo, "clipped_upper")
  testthat::expect_equal(derived$soft_derived$vivo_transformed, log10(2.5), tolerance = 1e-12)
})

testthat::test_that("soft-coupling start table converts scale-aware values into optimizer init", {
  tmp <- tempfile("soft-start-")
  dir.create(tmp)
  table_path <- file.path(tmp, "joint_soft_coupling_parameters_table.csv")
  utils::write.csv(
    data.frame(
      param_name = c(
        "log10_O2_crit",
        "delta__log10_O2_crit",
        "mu_hp",
        "delta__log10_mu_hp",
        "buffer_smax",
        "delta__buffer_smax"
      ),
      value = c(-0.30103, -1.39794, 1e-4, 1e-4, 0.8, -0.1),
      scale = c("log10", "log10", "natural", "natural", "identity", "identity")
    ),
    file = table_path,
    row.names = FALSE,
    quote = FALSE
  )
  meta <- data.frame(
    parameter = c("O2_crit", "mu_hp", "buffer_smax"),
    center_name = c("log10_O2_crit", "log10_mu_hp", "buffer_smax"),
    invitro_name = c("log10_O2_crit", "log10_mu_hp", "buffer_smax"),
    delta_name = c("delta__log10_O2_crit", "delta__log10_mu_hp", "delta__buffer_smax"),
    transform = c("log10_nonnegative", "log10", "identity"),
    delta_lower_t = c(-0.5, -0.5, -0.5),
    delta_upper_t = c(0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  init <- c(
    log10_O2_crit = 0,
    delta__log10_O2_crit = 0,
    log10_mu_hp = -3,
    delta__log10_mu_hp = 0,
    buffer_smax = 0.5,
    delta__buffer_smax = 0
  )
  lower <- c(
    log10_O2_crit = -1,
    delta__log10_O2_crit = -0.5,
    log10_mu_hp = -5,
    delta__log10_mu_hp = -0.5,
    buffer_smax = 0,
    delta__buffer_smax = -0.5
  )
  upper <- c(
    log10_O2_crit = log10(2.5),
    delta__log10_O2_crit = 0.5,
    log10_mu_hp = -2,
    delta__log10_mu_hp = 0.5,
    buffer_smax = 1,
    delta__buffer_smax = 0.5
  )

  out <- joint_backend_env$joint_apply_soft_coupling_start_table(
    init = init,
    lower = lower,
    upper = upper,
    soft_meta = meta,
    cfg_raw = list(),
    invivo_parameter_table = file.path(tmp, "parameter_table_O2.csv")
  )

  testthat::expect_equal(out$init[["log10_O2_crit"]], -0.30103, tolerance = 1e-12)
  testthat::expect_equal(out$init[["delta__log10_O2_crit"]], -1.39794, tolerance = 1e-12)
  testthat::expect_equal(out$lower[["delta__log10_O2_crit"]], -1.39794, tolerance = 1e-12)
  testthat::expect_equal(out$metadata$delta_lower_t[out$metadata$parameter == "O2_crit"], -1.39794, tolerance = 1e-12)
  testthat::expect_equal(out$init[["log10_mu_hp"]], -4, tolerance = 1e-12)
  testthat::expect_equal(out$init[["delta__log10_mu_hp"]], 2 * asinh(0.5) / log(10), tolerance = 1e-12)
  testthat::expect_equal(out$init[["buffer_smax"]], 0.8, tolerance = 1e-12)
  testthat::expect_equal(out$init[["delta__buffer_smax"]], -0.1, tolerance = 1e-12)
})

testthat::test_that("warm-up init combines best seed transformed parameters by ownership", {
  tmp <- tempfile("joint-warmup-")
  invivo_dir <- file.path(tmp, "invivo_seed")
  invitro_dir <- file.path(tmp, "invitro_seed")
  dir.create(invivo_dir, recursive = TRUE)
  dir.create(invitro_dir, recursive = TRUE)
  utils::write.table(
    data.frame(
      transformed_parameter = c("log10_O2_crit", "log10_alpha_o2", "log10_k_clear"),
      transformed_value = c(-0.2, 0.1, -1.5)
    ),
    file = file.path(invivo_dir, "fit_parameter_stages.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      transformed_parameter = c("log10_O2_crit", "log10_alpha_o2", "log10_sigma_growth"),
      transformed_value = c(-0.6, 0.5, -1.2)
    ),
    file = file.path(invitro_dir, "best_params_transformed.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  meta <- data.frame(
    parameter = "O2_crit",
    center_name = "log10_O2_crit",
    invitro_name = "log10_O2_crit",
    delta_name = "delta__log10_O2_crit",
    transform = "log10",
    stringsAsFactors = FALSE
  )
  init <- c(
    log10_O2_crit = 0,
    delta__log10_O2_crit = 0,
    log10_alpha_o2 = 0,
    log10_k_clear = -2,
    ivt__log10_sigma_growth = -2
  )
  lower <- c(
    log10_O2_crit = -1,
    delta__log10_O2_crit = -0.1,
    log10_alpha_o2 = -1,
    log10_k_clear = -3,
    ivt__log10_sigma_growth = -3
  )
  upper <- c(
    log10_O2_crit = 1,
    delta__log10_O2_crit = 0.1,
    log10_alpha_o2 = 1,
    log10_k_clear = 0,
    ivt__log10_sigma_growth = 0
  )
  invivo <- list(param_bundle = list(optimizer = list(init = init[!startsWith(names(init), "ivt__")])))

  out <- joint_backend_env$joint_apply_warmup_initial_values(
    init = init,
    lower = lower,
    upper = upper,
    soft_meta = meta,
    cfg_raw = list(
      joint_warmup_enable = TRUE,
      joint_warmup_invivo_seed_dir = invivo_dir,
      joint_warmup_invitro_seed_dir = invitro_dir
    ),
    invivo = invivo,
    invitro = list(),
    invivo_glucose = FALSE,
    ivt_extra_names = "log10_sigma_growth"
  )

  testthat::expect_true(out$enabled)
  testthat::expect_equal(out$init[["log10_O2_crit"]], -0.4, tolerance = 1e-12)
  testthat::expect_equal(out$init[["delta__log10_O2_crit"]], 0.4, tolerance = 1e-12)
  testthat::expect_equal(out$upper[["delta__log10_O2_crit"]], 0.4, tolerance = 1e-12)
  testthat::expect_equal(out$init[["log10_alpha_o2"]], 0.3, tolerance = 1e-12)
  testthat::expect_equal(out$init[["log10_k_clear"]], -1.5, tolerance = 1e-12)
  testthat::expect_equal(out$init[["ivt__log10_sigma_growth"]], -1.2, tolerance = 1e-12)
  testthat::expect_true(all(c(
    "soft_center_from_best_seed_mean",
    "soft_delta_from_best_seed_difference",
    "shared_mean_from_best_seeds",
    "invivo_best_seed",
    "invitro_best_seed"
  ) %in% out$applied$source))
})

testthat::test_that("warm-up DEoptim population uses bounded normal samples around init", {
  ctx <- list(
    init = c(a = 2, b = 0),
    lower = c(a = 1, b = -1),
    upper = c(a = 3, b = 1),
    joint_warmup = list(enabled = TRUE, sigmaN = joint_backend_env$joint_warmup_default_sigmaN())
  )
  set.seed(10)
  pop <- joint_backend_env$joint_deoptim_initial_population(ctx, NP_use = 20)

  testthat::expect_equal(unname(pop[1, ]), unname(ctx$init), tolerance = 1e-12)
  testthat::expect_true(all(pop[, "a"] >= 1 & pop[, "a"] <= 3))
  testthat::expect_true(all(pop[, "b"] >= -1 & pop[, "b"] <= 1))
  testthat::expect_true(stats::sd(pop[-1, "a"]) > 0)
  testthat::expect_true(stats::sd(pop[-1, "b"]) > 0)
})

testthat::test_that("soft-coupled p_mis_base reaches the in vitro objective", {
  old_invivo <- joint_backend_env$INVIVO_ENV
  old_invitro <- joint_backend_env$INVITRO_ENV
  on.exit({
    assign("INVIVO_ENV", old_invivo, envir = joint_backend_env)
    assign("INVITRO_ENV", old_invitro, envir = joint_backend_env)
  }, add = TRUE)

  mock_invivo <- new.env(parent = globalenv())
  mock_invivo$evaluate_objective_components <- function(par_t, scenarios, cfg) {
    list(L = 0)
  }
  mock_invivo$decode_params <- function(par_t, fit_treatment, fit_tau_O2, cfg) {
    list(p_mis_base = 1e-4)
  }

  seen_p_mis_base <- NA_real_
  mock_invitro <- new.env(parent = globalenv())
  mock_invitro$ivt_optim_par_to_run_params <- function(par_t, cfg) {
    list(p_mis_base = 10^par_t[["log10_p_mis_base"]])
  }
  mock_invitro$ivt_objective_components <- function(run_params,
                                                    fit_objects,
                                                    cfg,
                                                    fallback_max_passage_days,
                                                    growth_weight,
                                                    ploidy_weight,
                                                    flow_weight) {
    seen_p_mis_base <<- run_params$p_mis_base
    list(objective = 0, n_growth_negative_pred = 0)
  }
  assign("INVIVO_ENV", mock_invivo, envir = joint_backend_env)
  assign("INVITRO_ENV", mock_invitro, envir = joint_backend_env)

  meta <- data.frame(
    parameter = "p_mis_base",
    center_name = "log10_p_mis_base",
    invitro_name = "log10_p_mis_base",
    delta_name = "delta__log10_p_mis_base",
    transform = "log10",
    center_lower_t = -8,
    center_upper_t = -2,
    invivo_lower_t = -8,
    invivo_upper_t = -2,
    invitro_lower_t = -8,
    invitro_upper_t = -2,
    sigma_delta = 0.35,
    stringsAsFactors = FALSE
  )
  ctx <- list(
    init = c(log10_p_mis_base = -5, delta__log10_p_mis_base = 2),
    invivo_names = "log10_p_mis_base",
    joint_soft_coupling = list(enabled = TRUE, metadata = meta),
    invivo = list(scenarios = list(), cfg = list(fit_treatment = FALSE, fit_tau_O2 = FALSE)),
    invitro = list(
      spec = data.frame(
        param_name = "log10_p_mis_base",
        init = -5,
        lower = -8,
        upper = -2,
        stringsAsFactors = FALSE
      ),
      cfg = list(),
      fit_objects = list()
    ),
    ivt_extra_prefixed = character(0),
    ivt_extra_names = character(0),
    invitro_clip_lower = c(log10_p_mis_base = -8),
    invitro_clip_upper = c(log10_p_mis_base = -2),
    joint_weight_invivo = 1,
    joint_weight_invitro = 1,
    joint_invitro_growth_weight = 1,
    joint_invitro_ploidy_weight = 1,
    joint_invitro_flow_weight = 1,
    joint_restriction = FALSE,
    joint_constraint_penalty = 1e9,
    joint_require_invivo_pred1000_ploidy_gt2 = FALSE,
    joint_require_invitro_growth_nonnegative = FALSE,
    joint_require_invitro_ploidy_phenotype = FALSE
  )

  comp <- joint_backend_env$joint_objective_components(ctx$init, ctx)

  testthat::expect_equal(comp$invivo_run_params$p_mis_base, 1e-4, tolerance = 1e-12)
  testthat::expect_equal(comp$soft_coupling_derived$vitro_transformed, -6, tolerance = 1e-12)
  testthat::expect_equal(comp$invitro_run_params$p_mis_base, 1e-6, tolerance = 1e-12)
  testthat::expect_equal(seen_p_mis_base, 1e-6, tolerance = 1e-12)
})
