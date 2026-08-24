testthat::test_that("in vitro DEoptim metadata distinguishes early stop from itermax", {
  backend_path <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "util",
    "o2_supply_demand_map_fit_invitro_backend.R"
  )
  expression_list <- parse(backend_path)
  env <- new.env(parent = globalenv())
  env$`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
  wanted <- c("invitro_deoptim_iter_completed", "invitro_deoptim_stop_reason")
  for (expr in expression_list) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% wanted) {
      eval(expr, envir = env)
    }
  }

  early_fit <- list(optim = list(iter = 321L))
  capped_fit <- list(optim = list(iter = 500L))
  testthat::expect_identical(env$invitro_deoptim_iter_completed(early_fit, 500L), 321L)
  testthat::expect_identical(env$invitro_deoptim_iter_completed(capped_fit, 500L), 500L)
  testthat::expect_identical(
    env$invitro_deoptim_stop_reason(321L, 500L),
    "early_stop_reltol_or_steptol"
  )
  testthat::expect_identical(
    env$invitro_deoptim_stop_reason(500L, 500L),
    "itermax_reached"
  )
})

testthat::test_that("extra-results backfills legacy in vitro DEoptim metadata from fit_result.rds", {
  utils_path <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "analysis",
    "fit_results",
    "extra_results_analysis_utils.R"
  )
  expression_list <- parse(utils_path)
  env <- new.env(parent = globalenv())
  env$`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
  env$as_num <- function(x, default = NA_real_) {
    value <- suppressWarnings(as.numeric(x))
    if (length(value) != 1L || !is.finite(value)) default else value
  }
  env$o2sd_missing_summary_value <- function(x) {
    text <- trimws(as.character(x))
    is.na(x) | !nzchar(text) | tolower(text) %in% c("na", "nan", "null")
  }
  wanted <- c("summary_metric_value", "supplement_deoptim_metrics_from_fit_result")
  for (expr in expression_list) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% wanted) {
      eval(expr, envir = env)
    }
  }

  seed_dir <- tempfile("legacy-invitro-seed-")
  dir.create(seed_dir)
  on.exit(unlink(seed_dir, recursive = TRUE, force = TRUE), add = TRUE)
  saveRDS(
    list(deoptim = list(optim = list(iter = 321L))),
    file.path(seed_dir, "fit_result.rds")
  )
  summary_vals <- c(
    fit_mode = "fit_invitro",
    itermax = "500",
    optimizer_deoptim_objective = "12.5"
  )
  supplemented <- env$supplement_deoptim_metrics_from_fit_result(summary_vals, seed_dir)

  testthat::expect_identical(supplemented[["optimizer_interrupted"]], "FALSE")
  testthat::expect_identical(supplemented[["optimizer_iter_completed"]], "321")
  testthat::expect_identical(supplemented[["optimizer_iter_target"]], "500")
  testthat::expect_identical(
    supplemented[["deoptim_stop_reason"]],
    "early_stop_reltol_or_steptol"
  )
})
