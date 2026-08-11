testthat::test_that("in vitro default death mode is ploidy-related", {
  env <- new.env(parent = globalenv())
  source(
    file.path(
      repo_info$root,
      "oxygen",
      "code",
      "O2_supply_demand_MAP",
      "util",
      "o2_supply_demand_map_invitro_utils.R"
    ),
    local = env,
    chdir = TRUE
  )

  cfg <- env$ivt_build_default_cfg(
    repo_root = file.path(repo_info$root, "oxygen")
  )

  testthat::expect_identical(cfg$ploidy_O2_death, "ploidy_related")
})

testthat::test_that("standalone in-vitro YAML controls the effective simulation config", {
  testthat::skip_if_not_installed("yaml")
  env <- new.env(parent = globalenv())
  source(
    file.path(
      repo_info$root,
      "oxygen",
      "code",
      "O2_supply_demand_MAP",
      "util",
      "o2_supply_demand_map_invitro_utils.R"
    ),
    local = env,
    chdir = TRUE
  )

  oxygen_root <- file.path(repo_info$root, "oxygen")
  config_path <- file.path(
    oxygen_root,
    "config",
    "O2_supply_demand_invitro.yaml"
  )
  cfg <- env$.ivt_build_cfg_from_config(
    repo_root = oxygen_root,
    config_path = config_path,
    parameter_table = env$ivt_parameter_table_path(oxygen_root)
  )

  testthat::expect_false(cfg$Crowding)
  testthat::expect_identical(cfg$crowding, "logistic")
  testthat::expect_equal(cfg$K, 1e12)
  testthat::expect_equal(cfg$DT, 0.05)
  testthat::expect_true(cfg$fixed_oxygen)
  testthat::expect_false(cfg$o2_burden_feedback)
  testthat::expect_identical(cfg$config_fit_mode, "fit_invitro")
  testthat::expect_match(cfg$config_sha256, "^[0-9a-f]{64}$")
  testthat::expect_identical(
    cfg$config_path,
    normalizePath(config_path, mustWork = TRUE)
  )

  override_cfg <- env$.ivt_build_cfg_from_config(
    repo_root = oxygen_root,
    config_path = config_path,
    parameter_table = env$ivt_parameter_table_path(oxygen_root),
    cli_overrides = list(Crowding = TRUE, K = 2e7)
  )
  testthat::expect_true(override_cfg$Crowding)
  testthat::expect_equal(override_cfg$K, 2e7)
})

testthat::test_that("standalone in-vitro YAML rejects wrong modes and unknown keys", {
  testthat::skip_if_not_installed("yaml")
  env <- new.env(parent = globalenv())
  source(
    file.path(
      repo_info$root,
      "oxygen",
      "code",
      "O2_supply_demand_MAP",
      "util",
      "o2_supply_demand_map_invitro_utils.R"
    ),
    local = env,
    chdir = TRUE
  )

  wrong_mode <- tempfile("invitro_wrong_mode_", fileext = ".yaml")
  unknown_key <- tempfile("invitro_unknown_key_", fileext = ".yaml")
  on.exit(unlink(c(wrong_mode, unknown_key), force = TRUE), add = TRUE)
  writeLines(c("config_version: 1", "fit_mode: fit_invivo"), wrong_mode)
  writeLines(c(
    "config_version: 1",
    "fit_mode: fit_invitro",
    "unsupported_setting: 1"
  ), unknown_key)

  testthat::expect_error(
    env$.ivt_read_invitro_config(wrong_mode),
    "fit_mode must be 'fit_invitro'",
    fixed = TRUE
  )
  testthat::expect_error(
    env$.ivt_read_invitro_config(unknown_key),
    "Unsupported in-vitro config key"
  )
})

testthat::test_that("in-vitro launchers default to the standalone YAML", {
  workflow_root <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  array_workers <- c(
    file.path(workflow_root, "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub"),
    file.path(workflow_root, "Docker", "hpc", "array_workers", "submit_fit_seed_array_invitro_buffering.sub")
  )
  unified_launchers <- c(
    file.path(workflow_root, "runner", "run_o2_fit.sh"),
    file.path(workflow_root, "hpc", "submit", "submit_o2_fit.sh"),
    file.path(workflow_root, "Docker", "hpc", "submit", "submit_o2_fit.sh")
  )

  for (path in array_workers) {
    content <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(
      content,
      "DEFAULT_CONFIG_PATH=.*O2_supply_demand_invitro\\.yaml",
      info = path
    )
  }
  for (path in unified_launchers) {
    content <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(content, "FITTING_MODE.*invitro", info = path)
    testthat::expect_match(content, "O2_supply_demand_invitro\\.yaml", info = path)
    testthat::expect_match(content, "O2_supply_demand\\.yaml", info = path)
    testthat::expect_match(content, "INVITRO_CONFIG_PATH", info = path)
    testthat::expect_match(content, "invitro_config_path", info = path)
  }

  local_content <- readLines(unified_launchers[[1L]], warn = FALSE)
  testthat::expect_true(any(grepl(
    '"--config=${INVITRO_CONFIG_PATH}"',
    local_content,
    fixed = TRUE
  )))
  for (path in unified_launchers[-1L]) {
    lines <- readLines(path, warn = FALSE)
    invivo_start <- grep("^submit_invivo_array\\(\\)", lines)
    invitro_start <- grep("^submit_invitro_array\\(\\)", lines)
    joint_start <- grep("^submit_joint_array\\(\\)", lines)
    invivo_block <- lines[invivo_start:(invitro_start - 1L)]
    invitro_block <- lines[invitro_start:(joint_start - 1L)]
    testthat::expect_true(any(grepl(
      "CONFIG_PATH=${CONFIG_PATH}",
      invivo_block,
      fixed = TRUE
    )), info = path)
    testthat::expect_true(any(grepl(
      "CONFIG_PATH=${INVITRO_CONFIG_PATH}",
      invitro_block,
      fixed = TRUE
    )), info = path)
  }
})

testthat::test_that("in vitro parameter tables share passage-rate uncertainty", {
  table_names <- c(
    "parameter_table_invitro.csv",
    "parameter_table_invitro_buffering.csv",
    "parameter_table_invitro_wgd_bimodal.csv"
  )
  rows <- lapply(table_names, function(table_name) {
    path <- file.path(
      repo_info$root,
      "oxygen",
      "data",
      "O2_supply_demand",
      table_name
    )
    table <- utils::read.csv(path, stringsAsFactors = FALSE)
    row <- table[table$param_symbol == "sigma_growth", , drop = FALSE]
    testthat::expect_equal(nrow(row), 1L, info = table_name)
    row$parameter_table <- table_name
    row
  })
  rows <- do.call(rbind, rows)

  testthat::expect_equal(rows$init_value, rep(0.2, length(table_names)))
  testthat::expect_equal(rows$lower_bound, rep(0.05, length(table_names)))
  testthat::expect_equal(rows$upper_bound, rep(0.5, length(table_names)))
  testthat::expect_true(all(grepl("passage-average", rows$description, fixed = TRUE)))
  testthat::expect_true(all(grepl("day^-1", rows$description, fixed = TRUE)))
})
