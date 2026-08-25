invivo_stage2_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

read_invivo_stage2 <- function(...) {
  paste(readLines(file.path(invivo_stage2_root, ...), warn = FALSE), collapse = "\n")
}

testthat::test_that("in-vivo visualization is a pure simulation-table consumer", {
  entry_text <- read_invivo_stage2(
    "vis",
    "viz_invivo_model_O2_supply_demand_MAP_results.R"
  )
  helper_text <- read_invivo_stage2(
    "vis",
    "invivo",
    "o2_supply_demand_map_invivo_plot_utils.R"
  )
  text <- paste(entry_text, helper_text, sep = "\n")
  forbidden <- c(
    "readRDS\\(",
    "best_params",
    "fit_config",
    "dt_Gem",
    "all_ploidy",
    "model[/]",
    "cpp_o2",
    "prepare_data",
    "simulate_one",
    "generate_invivo_simulation_outputs",
    "source\\([^)]*simulation"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("pure in-vivo vis contains forbidden dependency:", pattern)
    )
  }
  testthat::expect_false(grepl("write[.]table", helper_text, perl = TRUE))
  testthat::expect_match(entry_text, 'file.path\\(fit_dir, "simulation", "invivo"\\)')
  testthat::expect_match(entry_text, 'file.path\\(fit_dir, "viz", "invivo"\\)')
  testthat::expect_match(entry_text, "viz_manifest.tsv", fixed = TRUE)
  testthat::expect_match(entry_text, 'options\\(bitmapType = "cairo"\\)')
})

testthat::test_that("in-vivo simulation producer is data-only and domain-split", {
  simulation_root <- file.path(invivo_stage2_root, "simulation", "invivo")
  producer_paths <- list.files(
    simulation_root,
    pattern = "[.]R$",
    recursive = TRUE,
    full.names = TRUE
  )
  testthat::expect_true(length(producer_paths) >= 7L)
  text <- paste(unlist(lapply(producer_paths, readLines, warn = FALSE)), collapse = "\n")
  for (pattern in c("ggplot\\s*\\(", "ggsave\\s*\\(", "grDevices::pdf\\s*\\(", "grDevices::png\\s*\\(")) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("in-vivo simulation contains plotting behavior:", pattern)
    )
  }
  for (domain in c("o2", "ploidy", "cin", "population", "functional_response")) {
    testthat::expect_true(dir.exists(file.path(simulation_root, domain)))
    testthat::expect_true(length(list.files(file.path(simulation_root, domain), pattern = "[.]R$")) >= 1L)
  }
})

testthat::test_that("in-vivo table and CLI contracts are explicit", {
  producer_text <- read_invivo_stage2(
    "simulation",
    "invivo",
    "generate_invivo_simulation_outputs.R"
  )
  functional_text <- read_invivo_stage2(
    "simulation",
    "invivo",
    "functional_response",
    "invivo_functional_response_simulation.R"
  )
  vis_text <- read_invivo_stage2(
    "vis",
    "viz_invivo_model_O2_supply_demand_MAP_results.R"
  )
  for (name in c(
    "functional_curve_oxygen",
    "functional_curve_oxygen_multi_ploidy",
    "functional_curve_ploidy",
    "functional_curve_ploidy_by_o2"
  )) {
    testthat::expect_match(functional_text, name, fixed = TRUE)
    testthat::expect_match(vis_text, paste0(name, ".tsv"), fixed = TRUE)
  }
  for (name in c("simulation_manifest.tsv", "output_schema.tsv", "observations.tsv")) {
    testthat::expect_match(producer_text, name, fixed = TRUE)
  }
  testthat::expect_match(producer_text, "argv$run_dir", fixed = TRUE)
  testthat::expect_match(producer_text, "argv$simulation_dir", fixed = TRUE)
  testthat::expect_match(vis_text, "argv$run_dir", fixed = TRUE)
  testthat::expect_match(vis_text, "argv$simulation_manifest", fixed = TRUE)
})

testthat::test_that("in-vivo consumer reports incomplete materialization clearly", {
  vis_path <- file.path(
    invivo_stage2_root,
    "vis",
    "viz_invivo_model_O2_supply_demand_MAP_results.R"
  )
  env <- new.env(parent = globalenv())
  source(vis_path, local = env)

  missing_path <- tempfile("missing_invivo_manifest_", fileext = ".tsv")
  testthat::expect_error(
    env$.read_key_value_manifest(missing_path),
    "Run the in-vivo simulation producer before visualization",
    fixed = TRUE
  )

  simulation_dir <- tempfile("incomplete_invivo_simulation_")
  dir.create(simulation_dir)
  on.exit(unlink(simulation_dir, recursive = TRUE, force = TRUE), add = TRUE)
  manifest <- c(
    schema_version = "o2sd-invivo-simulation-v1",
    status = "complete",
    fit_dir = simulation_dir,
    report_dt = "1",
    predict_report_dt = "1",
    predict_horizons = "100",
    start_with = "chr_number",
    N_UNIT = "22",
    N_MIN = "1",
    N_MAX = "133"
  )
  write.table(
    data.frame(file = character(), role = character(), n_rows = integer(), columns = character()),
    file.path(simulation_dir, "output_schema.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  testthat::expect_error(
    env$.validate_simulation_contract(simulation_dir, manifest, horizons = 100),
    "Simulation output is incomplete; missing table(s)",
    fixed = TRUE
  )
})
