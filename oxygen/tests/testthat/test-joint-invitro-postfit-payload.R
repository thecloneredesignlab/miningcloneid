joint_postfit_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

joint_postfit_env <- new.env(parent = globalenv())
for (source_path in c(
  file.path(joint_postfit_root, "util", "o2_supply_demand_map_common_semantics.R"),
  file.path(joint_postfit_root, "util", "o2_supply_demand_map_invitro_postfit_io_utils.R")
)) {
  sys.source(source_path, envir = joint_postfit_env, chdir = TRUE)
}

testthat::test_that("joint post-fit payload reads the parent best-components parameter list", {
  fit_dir <- tempfile("joint-postfit-")
  dir.create(fit_dir)
  fit_result <- list(
    ctx = list(invitro_cfg = list(N_UNIT = 22L)),
    best_components = list(
      invitro_run_params = list(tau_O2 = 0.25, p_wgd = 0.002),
      invitro = list(objective = 123)
    )
  )

  payload <- joint_postfit_env$ivt_sim_extract_payload(fit_result, fit_dir)

  testthat::expect_equal(payload$run_params$tau_O2, 0.25)
  testthat::expect_equal(payload$run_params$p_wgd, 0.002)
  testthat::expect_equal(payload$components$objective, 123)
})

testthat::test_that("best_params.tsv remains a fallback when a legacy RDS has no parameter list", {
  fit_dir <- tempfile("joint-postfit-tsv-")
  dir.create(fit_dir)
  utils::write.table(
    data.frame(parameter = c("tau_O2", "p_wgd"), value = c(0.4, 0.003)),
    file.path(fit_dir, "best_params.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  fit_result <- list(
    ctx = list(invitro_cfg = list(N_UNIT = 22L)),
    best_components = list(invitro = list(objective = 456))
  )

  payload <- joint_postfit_env$ivt_sim_extract_payload(fit_result, fit_dir)

  testthat::expect_equal(payload$run_params$tau_O2, 0.4)
  testthat::expect_equal(payload$run_params$p_wgd, 0.003)
})
