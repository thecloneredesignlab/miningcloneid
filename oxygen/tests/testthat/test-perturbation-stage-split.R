perturbation_root <- file.path(
  repo_info$root, "oxygen", "code", "O2_supply_demand_MAP"
)

perturbation_files <- list(
  shared = file.path(perturbation_root, "util", "o2_supply_demand_map_perturbation_utils.R"),
  mixed_wrapper = file.path(perturbation_root, "simulation", "simulate_invivo_mixed_ploidy_perturbations.R"),
  factorial_wrapper = file.path(perturbation_root, "simulation", "simulate_invivo_factorial_interaction.R"),
  mixed_sim = file.path(perturbation_root, "simulation", "perturbation", "generate_mixed_ploidy_perturbation_outputs.R"),
  factorial_sim = file.path(perturbation_root, "simulation", "perturbation", "generate_factorial_interaction_outputs.R"),
  mixed_analysis = file.path(perturbation_root, "analysis", "perturbation", "analyze_mixed_ploidy_perturbations.R"),
  factorial_analysis = file.path(perturbation_root, "analysis", "interactions", "analyze_factorial_interactions.R"),
  mixed_vis = file.path(perturbation_root, "vis", "perturbation", "plot_mixed_ploidy_perturbations.R"),
  factorial_vis = file.path(perturbation_root, "vis", "interactions", "plot_factorial_interactions.R")
)

read_perturbation_file <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

testthat::test_that("perturbation split entrypoints exist and parse", {
  testthat::expect_true(all(file.exists(unlist(perturbation_files))))
  for (path in unlist(perturbation_files)) {
    testthat::expect_silent(parse(path))
  }
})

testthat::test_that("canonical perturbation producers write numeric artifacts only", {
  forbidden <- "ggplot|ggsave|grDevices::pdf|pdf\\s*\\(|png\\s*\\(|[.]pdf|[.]png"
  for (path in c(perturbation_files$mixed_sim, perturbation_files$factorial_sim)) {
    text <- read_perturbation_file(path)
    testthat::expect_false(grepl(forbidden, text, ignore.case = TRUE, perl = TRUE), info = path)
    testthat::expect_match(text, "simulation_manifest.tsv", fixed = TRUE)
    testthat::expect_match(text, "burden_timecourse.tsv", fixed = TRUE)
  }
  factorial_text <- read_perturbation_file(perturbation_files$factorial_sim)
  testthat::expect_false(
    grepl("simulate_invivo_mixed_ploidy_perturbations.R", factorial_text, fixed = TRUE)
  )
})

testthat::test_that("perturbation analysis and visualization are strict consumers", {
  consumer_paths <- c(
    perturbation_files$mixed_analysis, perturbation_files$factorial_analysis,
    perturbation_files$mixed_vis, perturbation_files$factorial_vis
  )
  forbidden <- "best_params|fit_config[.]rds|fit_result[.]rds|model[/]|optimizer[/]|cpp_o2|run_.*simulation\\s*\\("
  for (path in consumer_paths) {
    text <- read_perturbation_file(path)
    testthat::expect_false(grepl(forbidden, text, ignore.case = TRUE, perl = TRUE), info = path)
    testthat::expect_match(text, "validate_artifact_manifest", fixed = TRUE)
  }
  for (path in c(perturbation_files$mixed_analysis, perturbation_files$factorial_analysis)) {
    text <- read_perturbation_file(path)
    testthat::expect_false(grepl("ggplot|ggsave|[.]pdf|[.]png", text, ignore.case = TRUE, perl = TRUE), info = path)
    testthat::expect_match(text, "analysis_manifest.tsv", fixed = TRUE)
  }
  for (path in c(perturbation_files$mixed_vis, perturbation_files$factorial_vis)) {
    testthat::expect_match(read_perturbation_file(path), "visualization_manifest.tsv", fixed = TRUE)
  }
})

testthat::test_that("legacy perturbation CLI paths are compatibility wrappers", {
  mixed <- read_perturbation_file(perturbation_files$mixed_wrapper)
  factorial <- read_perturbation_file(perturbation_files$factorial_wrapper)
  testthat::expect_match(mixed, "run_legacy_mixed_ploidy_pipeline", fixed = TRUE)
  testthat::expect_match(mixed, "generate_mixed_ploidy_perturbation_outputs.R", fixed = TRUE)
  testthat::expect_match(factorial, "run_legacy_factorial_interaction_pipeline", fixed = TRUE)
  testthat::expect_match(factorial, "generate_factorial_interaction_outputs.R", fixed = TRUE)
})

testthat::test_that("analysis and visualization fail rather than materialize missing simulations", {
  empty <- tempfile("missing_perturbation_sim_")
  dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE, force = TRUE), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  consumers <- c(
    perturbation_files$mixed_analysis, perturbation_files$factorial_analysis,
    perturbation_files$mixed_vis, perturbation_files$factorial_vis
  )
  for (path in consumers) {
    destination <- tempfile("missing_perturbation_output_")
    out <- suppressWarnings(system2(
      rscript,
      c(
        shQuote(path), paste0("--simulation_dir=", shQuote(empty)),
        paste0("--analysis_dir=", shQuote(empty)), paste0("--out_dir=", shQuote(destination))
      ),
      stdout = TRUE,
      stderr = TRUE
    ))
    testthat::expect_false(is.null(attr(out, "status")), info = path)
    testthat::expect_true(any(grepl("requires pre-existing upstream artifacts", out, fixed = TRUE)), info = path)
    testthat::expect_false(file.exists(file.path(empty, "simulation_manifest.tsv")), info = path)
    testthat::expect_false(dir.exists(destination), info = path)
  }
})
