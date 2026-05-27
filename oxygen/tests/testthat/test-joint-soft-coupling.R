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

testthat::test_that("expanded soft split parameters map to optimizer names", {
  split_params <- joint_backend_env$joint_parse_soft_param_list(
    "O2_crit,alpha_o2,mu_hp,p_misseg,k_o_mis,buffer_smax,buffer_beta,buffer_n_exp,n_O,gamma_growth"
  )
  name_map <- joint_backend_env$joint_soft_split_transformed_name_map(split_params)

  testthat::expect_equal(
    name_map,
    c(
      O2_crit = "log10_O2_crit",
      alpha_o2 = "log10_alpha_o2",
      mu_hp = "log10_mu_hp",
      p_misseg = "log10_p_misseg",
      k_o_mis = "log10_k_o_mis",
      buffer_smax = "buffer_smax",
      buffer_beta = "log10_buffer_beta",
      buffer_n_exp = "log10_buffer_n_exp",
      n_O = "n_O",
      gamma_growth = "gamma_growth"
    )
  )
})

testthat::test_that("soft coupling reports backend-bound clipping", {
  ctx <- make_soft_ctx(delta = 0.5, center = log10(2.5))
  derived <- joint_backend_env$joint_build_context_specific_transformed_vectors(ctx$init, ctx)

  testthat::expect_true(derived$soft_derived$vivo_clipped)
  testthat::expect_equal(derived$soft_derived$boundary_status_vivo, "clipped_upper")
  testthat::expect_equal(derived$soft_derived$vivo_transformed, log10(2.5), tolerance = 1e-12)
})
