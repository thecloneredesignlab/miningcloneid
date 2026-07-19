profile_stage_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

profile_stage_files <- list(
  simulation = file.path(
    profile_stage_root,
    "simulation",
    "invivo",
    "cin",
    "generate_live_effective_pms_outputs.R"
  ),
  analysis = file.path(
    profile_stage_root,
    "analysis",
    "profile_likelihood",
    "compare_sigma_burden_live_effective_pms.R"
  ),
  legacy_wrapper = file.path(
    profile_stage_root,
    "analysis",
    "profile_likelihood",
    "estimate_live_effective_pms.R"
  ),
  visualization = file.path(
    profile_stage_root,
    "vis",
    "profile_likelihood",
    "plot_sigma_burden_live_effective_pms.R"
  ),
  runner = file.path(
    profile_stage_root,
    "runner",
    "profile_likelihood",
    "run_live_effective_pms_comparison.R"
  )
)

read_profile_code <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  paste(lines, collapse = "\n")
}

write_profile_tsv <- function(tab, path) {
  write.table(
    tab,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

make_profile_simulation_fixture <- function(root) {
  seed_dir <- file.path(root, "seed1")
  simulation_dir <- file.path(seed_dir, "simulation", "invivo")
  dir.create(simulation_dir, recursive = TRUE)

  write_profile_tsv(
    data.frame(
      parameter = c(
        "p_mis_base",
        "p_misseg",
        "k_o_mis",
        "mu_hp",
        "gamma_mu",
        "O2_crit",
        "n_O"
      ),
      value = c(0.001, 0.1, 0.02, 0.04, 1.5, 1.0, 2.0)
    ),
    file.path(seed_dir, "best_params.tsv")
  )
  saveRDS(
    list(N_UNIT = 22, ploidy_O2_death = "ploidy_related"),
    file.path(seed_dir, "fit_config.rds")
  )
  write_profile_tsv(
    data.frame(
      key = c("schema_version", "status"),
      value = c("o2sd-invivo-simulation-v1", "complete")
    ),
    file.path(simulation_dir, "simulation_manifest.tsv")
  )
  write_profile_tsv(
    data.frame(
      harvest = c("h1", "h1", "h1", "h1"),
      cohort = c("2N", "2N", "4N", "4N"),
      dose = c(0, 0, 1, 1),
      day = c(0, 1, 0, 1),
      N = c(44, 44, 88, 88),
      fraction = 1
    ),
    file.path(simulation_dir, "ploidy_timecourse.tsv")
  )
  write_profile_tsv(
    data.frame(
      harvest = c("h1", "h1", "h1", "h1"),
      cohort = c("2N", "2N", "4N", "4N"),
      dose = c(0, 0, 1, 1),
      day = c(0, 1, 0, 1),
      o2_pct = c(5, 2, 4, 1)
    ),
    file.path(simulation_dir, "predict_burden_vs_o2.tsv")
  )
  write_profile_tsv(
    data.frame(
      N = c(44, 88),
      viability_after_ms = c(0.9, 0.8)
    ),
    file.path(simulation_dir, "functional_curve_ploidy.tsv")
  )
  list(seed_dir = seed_dir, simulation_dir = simulation_dir)
}

testthat::test_that("profile-likelihood stage entrypoints exist and parse", {
  testthat::expect_true(all(file.exists(unlist(profile_stage_files))))
  for (path in unlist(profile_stage_files)) {
    testthat::expect_silent(parse(path))
  }
})

testthat::test_that("live-effective-p_ms responsibilities are physically separated", {
  simulation_text <- read_profile_code(profile_stage_files$simulation)
  analysis_text <- read_profile_code(profile_stage_files$analysis)
  vis_text <- read_profile_code(profile_stage_files$visualization)
  wrapper_text <- paste(
    readLines(profile_stage_files$legacy_wrapper, warn = FALSE),
    collapse = "\n"
  )

  testthat::expect_false(grepl(
    "ggplot|ggsave|grDevices::|pdf\\s*\\(|png\\s*\\(",
    simulation_text,
    ignore.case = TRUE,
    perl = TRUE
  ))
  testthat::expect_match(simulation_text, "live_effective_pms_manifest.tsv", fixed = TRUE)

  testthat::expect_false(grepl(
    "ggplot|ggsave|pdf\\s*\\(|png\\s*\\(|system2\\s*\\(",
    analysis_text,
    ignore.case = TRUE,
    perl = TRUE
  ))
  testthat::expect_match(
    analysis_text,
    "sigma_burden_p_misseg_vs_live_cell_plot.tsv",
    fixed = TRUE
  )

  testthat::expect_false(grepl(
    "readRDS\\s*\\(|best_params|fit_result|model_O2_supply_demand_MAP|source\\([^\\n]*simulation",
    vis_text,
    ignore.case = TRUE,
    perl = TRUE
  ))
  testthat::expect_match(vis_text, "live_effective_pms_visualization_manifest.tsv", fixed = TRUE)
  testthat::expect_match(wrapper_text, "Deprecated compatibility wrapper", fixed = TRUE)
  testthat::expect_match(
    wrapper_text,
    "generate_live_effective_pms_outputs.R",
    fixed = TRUE
  )
  runner_text <- read_profile_code(profile_stage_files$runner)
  analysis_pos <- regexpr(
    "compare_sigma_burden_live_effective_pms.R",
    runner_text,
    fixed = TRUE
  )[[1L]]
  vis_pos <- regexpr(
    "plot_sigma_burden_live_effective_pms.R",
    runner_text,
    fixed = TRUE
  )[[1L]]
  testthat::expect_gt(analysis_pos, 0L)
  testthat::expect_gt(vis_pos, analysis_pos)
})

testthat::test_that("live-effective-p_ms simulation-analysis-vis chain honors artifacts", {
  fixture_root <- tempfile("profile_stage_fixture_")
  dir.create(fixture_root)
  on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)
  fixture <- make_profile_simulation_fixture(fixture_root)

  sim_env <- new.env(parent = globalenv())
  source(profile_stage_files$simulation, local = sim_env)
  live_dir <- file.path(
    fixture$simulation_dir,
    "cin",
    "live_effective_pms"
  )
  testthat::expect_output(
    sim_env$main(list(
      seed_dir = fixture$seed_dir,
      simulation_dir = fixture$simulation_dir,
      out_dir = live_dir
    )),
    "Overall live-weighted effective p_ms",
    fixed = TRUE
  )
  expected_simulation_files <- c(
    "live_effective_pms_context.tsv",
    "live_effective_pms_sample_day.tsv",
    "live_effective_pms_overall.tsv",
    "live_effective_pms_harvest_only.tsv",
    "live_effective_pms_cohort_all_days.tsv",
    "live_effective_pms_cohort_harvest_only.tsv",
    "live_effective_pms_schema.tsv",
    "live_effective_pms_manifest.tsv"
  )
  testthat::expect_true(all(file.exists(file.path(live_dir, expected_simulation_files))))
  sample_day <- read.delim(
    file.path(live_dir, "live_effective_pms_sample_day.tsv"),
    check.names = FALSE
  )
  testthat::expect_equal(nrow(sample_day), 4L)
  testthat::expect_true(all(sample_day$live_weighted_effective_p_ms > 0.001))
  testthat::expect_true(all(sample_day$live_weighted_effective_p_ms < 0.101))

  for (sigma_cap in c("0.05", "0.15")) {
    for (seed in c("seed1", "seed2")) {
      target <- file.path(
        fixture_root,
        paste0("run_", sigma_cap),
        seed,
        "simulation",
        "invivo",
        "cin",
        "live_effective_pms"
      )
      dir.create(target, recursive = TRUE)
      testthat::expect_true(file.copy(
        file.path(live_dir, list.files(live_dir)),
        target,
        overwrite = TRUE
      ) |> all())
    }
  }

  analysis_env <- new.env(parent = globalenv())
  source(profile_stage_files$analysis, local = analysis_env)
  analysis_dir <- file.path(fixture_root, "analysis")
  testthat::expect_message(
    analysis_env$main(list(
      run_dir_template = file.path(fixture_root, "run_{sigma}"),
      sigma_caps = "0.05,0.15",
      out_dir = analysis_dir,
      n_cores = 1L
    )),
    "Wrote analysis manifest"
  )
  by_seed <- read.delim(
    file.path(analysis_dir, "sigma_burden_live_effective_pms_by_seed.tsv"),
    check.names = FALSE
  )
  plot_tab <- read.delim(
    file.path(analysis_dir, "sigma_burden_p_misseg_vs_live_cell_plot.tsv"),
    check.names = FALSE
  )
  testthat::expect_equal(nrow(by_seed), 4L)
  testthat::expect_equal(nrow(plot_tab), 8L)
  testthat::expect_false(any(grepl("[.]pdf$|[.]png$", list.files(analysis_dir))))

  vis_env <- new.env(parent = globalenv())
  source(profile_stage_files$visualization, local = vis_env)
  figure_dir <- file.path(fixture_root, "figures")
  testthat::expect_message(
    vis_env$main(list(analysis_dir = analysis_dir, out_dir = figure_dir)),
    "Wrote visualization manifest"
  )
  testthat::expect_true(file.exists(file.path(
    figure_dir,
    "sigma_burden_p_misseg_vs_live_cell_violin.pdf"
  )))
  testthat::expect_true(file.exists(file.path(
    figure_dir,
    "live_effective_pms_visualization_manifest.tsv"
  )))
})

testthat::test_that("profile analysis fails when simulation artifacts are absent", {
  analysis_env <- new.env(parent = globalenv())
  source(profile_stage_files$analysis, local = analysis_env)
  missing_dir <- tempfile("missing_live_effective_")
  dir.create(missing_dir)
  testthat::expect_error(
    analysis_env$validate_live_effective_contract(missing_dir),
    "Run simulation/invivo/cin/generate_live_effective_pms_outputs.R",
    fixed = TRUE
  )
  testthat::expect_false(file.exists(file.path(
    missing_dir,
    "live_effective_pms_manifest.tsv"
  )))
})
