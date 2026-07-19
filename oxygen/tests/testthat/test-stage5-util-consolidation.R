stage5_map_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

stage5_source_functions <- function(path) {
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env, chdir = TRUE)
  names <- ls(env, all.names = TRUE)
  names <- names[vapply(
    names,
    function(name) is.function(get(name, envir = env, inherits = FALSE)),
    logical(1)
  )]
  stats::setNames(lapply(names, function(name) get(name, envir = env, inherits = FALSE)), names)
}

stage5_expect_same_functions <- function(actual, expected) {
  testthat::expect_identical(names(actual), names(expected))
  for (name in names(expected)) {
    testthat::expect_identical(formals(actual[[name]]), formals(expected[[name]]))
    testthat::expect_identical(body(actual[[name]]), body(expected[[name]]))
  }
}

testthat::test_that("best-fit shared utilities have canonical top-level implementations", {
  old_root <- file.path(
    stage5_map_root,
    "analysis",
    "best_fit_parameter_feature",
    "util"
  )
  mappings <- c(
    cli_utils.R = "o2_supply_demand_map_bpf_cli_utils.R",
    path_utils.R = "o2_supply_demand_map_bpf_path_utils.R",
    table_io_utils.R = "o2_supply_demand_map_bpf_table_io_utils.R",
    report_utils.R = "o2_supply_demand_map_bpf_report_utils.R",
    curve_classification_utils.R = "o2_supply_demand_map_curve_classification_utils.R"
  )
  for (old_name in names(mappings)) {
    old_path <- file.path(old_root, old_name)
    canonical_path <- file.path(stage5_map_root, "util", mappings[[old_name]])
    testthat::expect_true(file.exists(old_path))
    testthat::expect_true(file.exists(canonical_path))
    testthat::expect_lt(length(readLines(old_path, warn = FALSE)), 40L)
    testthat::expect_match(
      paste(readLines(old_path, warn = FALSE), collapse = "\n"),
      mappings[[old_name]],
      fixed = TRUE
    )
    stage5_expect_same_functions(
      stage5_source_functions(old_path),
      stage5_source_functions(canonical_path)
    )
  }
})

testthat::test_that("all fixed-O2 compatibility paths share one canonical implementation", {
  canonical_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_fixed_o2_utils.R"
  )
  analysis_wrapper <- file.path(
    stage5_map_root,
    "analysis",
    "best_fit_parameter_feature",
    "util",
    "fixed_o2_shared_utils.R"
  )
  simulation_wrappers <- c(
    file.path(stage5_map_root, "simulation", "fix_o2_simulation_shared_utils.R"),
    file.path(
      stage5_map_root,
      "simulation",
      "o2",
      "fixed_o2",
      "fixed_o2_simulation_utils.R"
    )
  )
  canonical <- stage5_source_functions(canonical_path)
  analysis_api <- stage5_source_functions(analysis_wrapper)
  stage5_expect_same_functions(analysis_api, canonical)

  simulation_names <- setdiff(
    names(canonical),
    c(
      "bpf_default_fixed_o2_grid",
      "bpf_default_dense_attractor_o2_grid",
      "bpf_o2_key",
      "bpf_o2_slug",
      "fixed_o2_shared_utils_dir"
    )
  )
  expected_simulation_api <- canonical[simulation_names]
  for (wrapper in simulation_wrappers) {
    testthat::expect_lt(length(readLines(wrapper, warn = FALSE)), 50L)
    stage5_expect_same_functions(
      stage5_source_functions(wrapper),
      expected_simulation_api
    )
  }

  run_script <- paste(
    readLines(
      file.path(
        stage5_map_root,
        "simulation",
        "o2",
        "fixed_o2",
        "run_fixed_o2_simulation.R"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  testthat::expect_match(run_script, "o2_supply_demand_map_fixed_o2_utils.R", fixed = TRUE)
  testthat::expect_false(grepl('source\\([^)]*"fixed_o2_simulation_utils.R"', run_script, perl = TRUE))
})

testthat::test_that("fixed-O2 mode semantics are defined once in canonical util", {
  canonical_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_fixed_o2_mode_utils.R"
  )
  analysis_path <- file.path(
    stage5_map_root,
    "analysis",
    "fixed_o2",
    "fixed_o2_analysis_utils.R"
  )
  simulation_path <- file.path(
    stage5_map_root,
    "simulation",
    "o2",
    "fixed_o2",
    "run_fixed_o2_simulation.R"
  )
  canonical <- stage5_source_functions(canonical_path)
  analysis_api <- stage5_source_functions(analysis_path)
  stage5_expect_same_functions(
    analysis_api[names(canonical)],
    canonical
  )

  for (path in c(analysis_path, simulation_path)) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(
      text,
      "o2_supply_demand_map_fixed_o2_mode_utils.R",
      fixed = TRUE
    )
    for (name in names(canonical)) {
      definition_pattern <- paste0(
        "(?m)^\\Q",
        name,
        "\\E\\s*<-\\s*function"
      )
      testthat::expect_false(
        grepl(definition_pattern, text, perl = TRUE),
        info = paste(path, name)
      )
    }
  }
})

testthat::test_that("parameter-landscape canonical analysis and deprecated bridge keep directional boundaries", {
  canonical <- file.path(
    stage5_map_root,
    "analysis",
    "parameter_landscape_clustering",
    "parameter_landscape_utils.R"
  )
  canonical_text <- paste(readLines(canonical, warn = FALSE), collapse = "\n")
  testthat::expect_match(canonical_text, "parameter_landscape_analysis_utils.R", fixed = TRUE)
  testthat::expect_false(grepl("simulation/", canonical_text, fixed = TRUE))
  testthat::expect_false(grepl("viz_invivo_model_O2_supply_demand_MAP_results", canonical_text, fixed = TRUE))

  canonical_env <- new.env(parent = globalenv())
  sys.source(canonical, envir = canonical_env, chdir = TRUE)
  testthat::expect_true(exists("o2pl_materialize_embedding", envir = canonical_env, inherits = TRUE))
  testthat::expect_false(exists("invivo_simulation_env", envir = canonical_env, inherits = TRUE))

  bridge <- file.path(
    stage5_map_root,
    "analysis",
    "best_fit_parameter_feature",
    "util",
    "parameter_landscape_shared_utils.R"
  )
  bridge_text <- paste(readLines(bridge, warn = FALSE), collapse = "\n")
  testthat::expect_match(bridge_text, "Deprecated compatibility", fixed = TRUE)
  testthat::expect_match(bridge_text, "parameter_landscape_simulation_utils.R", fixed = TRUE)
  testthat::expect_match(bridge_text, "parameter_landscape_invivo_feature_simulation.R", fixed = TRUE)
  testthat::expect_match(bridge_text, "parameter_landscape_analysis_utils.R", fixed = TRUE)
  testthat::expect_false(grepl("viz_invivo_model_O2_supply_demand_MAP_results", bridge_text, fixed = TRUE))

  bridge_env <- new.env(parent = globalenv())
  sys.source(bridge, envir = bridge_env, chdir = TRUE)
  bridge_env$default_oxygen_dir <- function() file.path(repo_info$root, "oxygen")
  simulation_env <- bridge_env$invivo_simulation_env()
  for (name in c("normalize_cfg_for_simulation", "read_run_params", "simulate_one_full_horizon")) {
    testthat::expect_true(exists(name, envir = simulation_env, inherits = TRUE))
  }
  testthat::expect_false(exists("render_invivo_visualizations", envir = simulation_env, inherits = TRUE))
})

testthat::test_that("legacy best-fit runner CLI still resolves migrated wrappers", {
  runner <- file.path(
    stage5_map_root,
    "analysis",
    "best_fit_parameter_feature",
    "runner.R"
  )
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(runner, "--workflow=fixed_o2", "--dry_run=true"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  testthat::expect_identical(as.integer(status), 0L)
  testthat::expect_true(any(grepl("FixO2_invivo.R", output, fixed = TRUE)))
})

testthat::test_that("process-fingerprint helpers have one canonical util implementation", {
  canonical_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_process_fingerprint_utils.R"
  )
  analysis_loader <- file.path(
    stage5_map_root,
    "analysis",
    "process_fingerprints",
    "process_fingerprint_utils.R"
  )
  simulation_helpers <- file.path(
    stage5_map_root,
    "simulation",
    "process_fingerprints",
    "process_fingerprint_simulation_legacy_utils.R"
  )
  testthat::expect_true(all(file.exists(c(
    canonical_path,
    analysis_loader,
    simulation_helpers
  ))))
  testthat::expect_lt(length(readLines(analysis_loader, warn = FALSE)), 45L)
  testthat::expect_match(
    paste(readLines(analysis_loader, warn = FALSE), collapse = "\n"),
    "o2_supply_demand_map_process_fingerprint_utils.R",
    fixed = TRUE
  )

  canonical <- stage5_source_functions(canonical_path)
  analysis_api <- stage5_source_functions(analysis_loader)
  stage5_expect_same_functions(analysis_api, canonical)

  simulation_api <- stage5_source_functions(simulation_helpers)
  shared_names <- setdiff(
    names(canonical),
    c("o2ipa_find_workflow_root", "o2ipa_extract_param")
  )
  stage5_expect_same_functions(
    simulation_api[shared_names],
    canonical[shared_names]
  )
  simulation_text <- paste(
    readLines(simulation_helpers, warn = FALSE),
    collapse = "\n"
  )
  for (name in shared_names) {
    definition_pattern <- paste0(
      "(?m)^",
      gsub(".", "\\\\.", name, fixed = TRUE),
      "\\s*<-\\s*function"
    )
    testthat::expect_false(
      grepl(definition_pattern, simulation_text, perl = TRUE),
      info = name
    )
  }
})

testthat::test_that("ploidy-regime helpers have one canonical util implementation", {
  canonical_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_ploidy_regime_utils.R"
  )
  analysis_loader <- file.path(
    stage5_map_root,
    "analysis",
    "process_fingerprints",
    "ploidy_regime_utils.R"
  )
  simulation_helpers <- file.path(
    stage5_map_root,
    "simulation",
    "process_fingerprints",
    "ploidy_regime_simulation_utils.R"
  )
  testthat::expect_true(all(file.exists(c(
    canonical_path,
    analysis_loader,
    simulation_helpers
  ))))
  testthat::expect_lt(length(readLines(analysis_loader, warn = FALSE)), 45L)
  testthat::expect_match(
    paste(readLines(analysis_loader, warn = FALSE), collapse = "\n"),
    "o2_supply_demand_map_ploidy_regime_utils.R",
    fixed = TRUE
  )

  canonical <- stage5_source_functions(canonical_path)
  analysis_api <- stage5_source_functions(analysis_loader)
  stage5_expect_same_functions(analysis_api, canonical)

  simulation_api <- stage5_source_functions(simulation_helpers)
  stage5_expect_same_functions(
    simulation_api[names(canonical)],
    canonical
  )
  simulation_text <- paste(
    readLines(simulation_helpers, warn = FALSE),
    collapse = "\n"
  )
  for (name in names(canonical)) {
    definition_pattern <- paste0(
      "(?m)^\\Q",
      name,
      "\\E",
      "\\s*<-\\s*function"
    )
    testthat::expect_false(
      grepl(definition_pattern, simulation_text, perl = TRUE),
      info = name
    )
  }
})

testthat::test_that("post-fit input and model-probe helpers are not duplicated across layers", {
  input_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_postfit_input_utils.R"
  )
  probe_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_postfit_probe_utils.R"
  )
  fixed_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_fixed_o2_utils.R"
  )
  process_path <- file.path(
    stage5_map_root,
    "simulation",
    "process_fingerprints",
    "process_fingerprint_simulation_legacy_utils.R"
  )
  ploidy_path <- file.path(
    stage5_map_root,
    "simulation",
    "process_fingerprints",
    "ploidy_regime_simulation_utils.R"
  )

  input_api <- stage5_source_functions(input_path)
  probe_api <- stage5_source_functions(probe_path)
  fixed_api <- stage5_source_functions(fixed_path)
  process_api <- stage5_source_functions(process_path)
  ploidy_api <- stage5_source_functions(ploidy_path)
  stage5_expect_same_functions(fixed_api[names(input_api)], input_api)
  stage5_expect_same_functions(process_api[names(input_api)], input_api)
  stage5_expect_same_functions(fixed_api[names(probe_api)], probe_api)
  stage5_expect_same_functions(ploidy_api[names(probe_api)], probe_api)

  consumers <- list(
    input = c(fixed_path, process_path),
    probe = c(fixed_path, ploidy_path)
  )
  canonical_apis <- list(input = input_api, probe = probe_api)
  canonical_files <- c(
    input = "o2_supply_demand_map_postfit_input_utils.R",
    probe = "o2_supply_demand_map_postfit_probe_utils.R"
  )
  for (kind in names(consumers)) {
    for (path in consumers[[kind]]) {
      text <- paste(readLines(path, warn = FALSE), collapse = "\n")
      testthat::expect_match(text, canonical_files[[kind]], fixed = TRUE)
      for (name in names(canonical_apis[[kind]])) {
        definition_pattern <- paste0(
          "(?m)^\\Q",
          name,
          "\\E\\s*<-\\s*function"
        )
        testthat::expect_false(
          grepl(definition_pattern, text, perl = TRUE),
          info = paste(path, name)
        )
      }
    }
  }
})

testthat::test_that("remaining fixed-O2 and perturbation shared helpers have canonical homes", {
  fixed_format <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_fixed_o2_format_utils.R"
  )
  fixed_table <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_fixed_o2_table_utils.R"
  )
  fixed_validation <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_fixed_o2_validation_utils.R"
  )
  perturbation <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_perturbation_utils.R"
  )
  testthat::expect_true("num_path_tag" %in% names(stage5_source_functions(fixed_format)))
  testthat::expect_setequal(
    names(stage5_source_functions(fixed_table)),
    c("fo2_read_tsv", "fo2_write_tsv")
  )
  testthat::expect_identical(
    names(stage5_source_functions(fixed_validation)),
    "fixo2_validate_mode_reference_o2"
  )
  testthat::expect_true(
    "format_value_label" %in% names(stage5_source_functions(perturbation))
  )

  fixed_simulation <- file.path(
    stage5_map_root,
    "simulation",
    "o2",
    "fixed_o2",
    "run_fixed_o2_simulation.R"
  )
  fixed_analysis <- file.path(
    stage5_map_root,
    "analysis",
    "fixed_o2",
    "run_fixed_o2_analysis.R"
  )
  fixed_runner_loader <- file.path(
    stage5_map_root,
    "runner",
    "fixed_o2",
    "fixed_o2_pipeline_loader.R"
  )
  fixed_legacy <- file.path(
    stage5_map_root,
    "runner",
    "fixed_o2",
    "fixed_o2_legacy_pipeline_functions.R"
  )
  fixed_simulation_text <- paste(readLines(fixed_simulation, warn = FALSE), collapse = "\n")
  fixed_analysis_text <- paste(readLines(fixed_analysis, warn = FALSE), collapse = "\n")
  fixed_runner_loader_text <- paste(readLines(fixed_runner_loader, warn = FALSE), collapse = "\n")
  fixed_legacy_text <- paste(readLines(fixed_legacy, warn = FALSE), collapse = "\n")
  testthat::expect_match(fixed_simulation_text, basename(fixed_format), fixed = TRUE)
  testthat::expect_match(fixed_simulation_text, basename(fixed_validation), fixed = TRUE)
  testthat::expect_match(fixed_analysis_text, basename(fixed_table), fixed = TRUE)
  testthat::expect_match(fixed_runner_loader_text, basename(fixed_table), fixed = TRUE)
  testthat::expect_match(fixed_runner_loader_text, basename(fixed_validation), fixed = TRUE)
  for (name in c("num_path_tag", "fixo2_validate_mode_reference_o2")) {
    pattern <- paste0("(?m)^\\Q", name, "\\E\\s*<-\\s*function")
    testthat::expect_false(grepl(pattern, fixed_simulation_text, perl = TRUE))
    testthat::expect_false(grepl(pattern, fixed_legacy_text, perl = TRUE))
  }
  for (name in c("fo2_read_tsv", "fo2_write_tsv")) {
    pattern <- paste0("(?m)^\\Q", name, "\\E\\s*<-\\s*function")
    testthat::expect_false(grepl(pattern, fixed_analysis_text, perl = TRUE))
    testthat::expect_false(grepl(pattern, fixed_legacy_text, perl = TRUE))
  }

  perturbation_consumers <- c(
    file.path(
      stage5_map_root,
      "simulation",
      "perturbation",
      "generate_factorial_interaction_outputs.R"
    ),
    file.path(
      stage5_map_root,
      "vis",
      "interactions",
      "plot_factorial_interactions.R"
    )
  )
  for (path in perturbation_consumers) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(text, basename(perturbation), fixed = TRUE)
    testthat::expect_false(grepl(
      "(?m)^format_value_label\\s*<-\\s*function",
      text,
      perl = TRUE
    ))
  }
})

testthat::test_that("report-only shared helpers have one canonical util implementation", {
  canonical_path <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_report_utils.R"
  )
  testthat::expect_true(file.exists(canonical_path))
  canonical <- stage5_source_functions(canonical_path)
  expected <- c(
    "o2sd_report_read_table_optional",
    "o2sd_report_safe_read_delim",
    "o2sd_report_parse_equals_args",
    "o2sd_report_null_coalesce_simple",
    "o2sd_report_null_coalesce",
    "o2sd_report_truthy",
    "o2sd_report_format_numeric_like",
    "o2sd_report_transformed_param_to_natural",
    "o2sd_report_annotate_parameter_table_for_report",
    "o2sd_report_pandoc_available",
    "o2sd_report_pdflatex_available",
    "o2sd_report_magick_available",
    "o2sd_report_ghostscript_bin",
    "o2sd_report_gs_available",
    "o2sd_report_base64enc_available",
    "o2sd_report_render_pdf_preview_png_gs",
    "o2sd_report_file_to_data_uri"
  )
  testthat::expect_setequal(names(canonical), expected)

  testthat::expect_identical(
    canonical$o2sd_report_truthy(c(TRUE, FALSE, NA)),
    c(TRUE, FALSE, FALSE)
  )
  testthat::expect_identical(
    canonical$o2sd_report_truthy(c("yes", "0", "ON", NA_character_)),
    c(TRUE, FALSE, TRUE, FALSE)
  )
  testthat::expect_identical(
    canonical$o2sd_report_transformed_param_to_natural(
      c("ivt__log10_alpha", "logit_beta", "gamma")
    ),
    c("alpha", "beta", "gamma")
  )
  testthat::expect_identical(
    canonical$o2sd_report_format_numeric_like(c("1", "1.25", "0.0001", "label")),
    c("1", "1.250", "1.000e-04", "label")
  )
  testthat::expect_identical(
    canonical$o2sd_report_parse_equals_args(c("--alpha=1", "ignored", "--label=a=b")),
    list(alpha = "1", label = "a=b")
  )
  testthat::expect_identical(canonical$o2sd_report_null_coalesce(NA_character_, "fallback"), "fallback")
  testthat::expect_null(
    canonical$o2sd_report_read_table_optional(file.path(tempdir(), "o2-missing-report-table.tsv"))
  )
  testthat::expect_equal(
    canonical$o2sd_report_safe_read_delim(file.path(tempdir(), "o2-missing-safe-table.tsv")),
    data.frame()
  )

  consumers <- list(
    render_fit_report.R = c(
      annotate_parameter_table_for_report = "o2sd_report_annotate_parameter_table_for_report",
      report_truthy = "o2sd_report_truthy",
      transformed_param_to_natural = "o2sd_report_transformed_param_to_natural",
      report_pandoc_available = "o2sd_report_pandoc_available",
      report_pdflatex_available = "o2sd_report_pdflatex_available",
      report_magick_available = "o2sd_report_magick_available",
      report_gs_available = "o2sd_report_gs_available",
      report_base64enc_available = "o2sd_report_base64enc_available",
      render_pdf_preview_png_gs = "o2sd_report_render_pdf_preview_png_gs",
      file_to_data_uri = "o2sd_report_file_to_data_uri"
    ),
    render_invitro_fit_report.R = c(
      read_table_optional = "o2sd_report_read_table_optional",
      format_numeric_like = "o2sd_report_format_numeric_like",
      annotate_parameter_table_for_report = "o2sd_report_annotate_parameter_table_for_report",
      report_truthy = "o2sd_report_truthy",
      transformed_param_to_natural = "o2sd_report_transformed_param_to_natural"
    ),
    render_fixo2_invivo_report.R = c(
      read_table_optional = "o2sd_report_read_table_optional",
      format_numeric_like = "o2sd_report_format_numeric_like"
    ),
    `fit_results/render_extra_results_report.R` = c(
      report_magick_available = "o2sd_report_magick_available",
      report_gs_available = "o2sd_report_gs_available",
      report_base64enc_available = "o2sd_report_base64enc_available",
      file_to_data_uri = "o2sd_report_file_to_data_uri",
      render_pdf_preview_png_gs = "o2sd_report_render_pdf_preview_png_gs",
      read_table_optional = "o2sd_report_read_table_optional",
      report_truthy = "o2sd_report_truthy"
    ),
    `multi_warmup/render_multi_warmup_results_report.R` = c(
      magick_available = "o2sd_report_magick_available",
      ghostscript_bin = "o2sd_report_ghostscript_bin",
      ghostscript_available = "o2sd_report_gs_available",
      render_pdf_preview_png_gs = "o2sd_report_render_pdf_preview_png_gs"
    ),
    `process_fingerprints/render_medium_o2_window_report.R` = c(
      parse_args = "o2sd_report_parse_equals_args",
      safe_read = "o2sd_report_safe_read_delim"
    ),
    `process_fingerprints/render_o2_ploidy_event_coupling_report.R` = c(
      parse_args = "o2sd_report_parse_equals_args",
      safe_read = "o2sd_report_safe_read_delim"
    )
  )
  for (relative_path in names(consumers)) {
    path <- file.path(stage5_map_root, "report", relative_path)
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(text, basename(canonical_path), fixed = TRUE)
    for (legacy_name in names(consumers[[relative_path]])) {
      canonical_name <- consumers[[relative_path]][[legacy_name]]
      testthat::expect_match(
        text,
        paste0(legacy_name, " <- ", canonical_name),
        fixed = TRUE,
        info = paste(relative_path, legacy_name)
      )
      definition_pattern <- paste0(
        "(?m)^\\Q",
        legacy_name,
        "\\E\\s*<-\\s*function"
      )
      testthat::expect_false(
        grepl(definition_pattern, text, perl = TRUE),
        info = paste(relative_path, legacy_name)
      )
    }
  }

  fixed_o2_report_paths <- c(
    file.path(stage5_map_root, "report", "render_fixo2_invivo_report.R"),
    file.path(stage5_map_root, "report", "fixed_o2", "render_fixed_o2_report.R")
  )
  for (path in fixed_o2_report_paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(text, "o2sd_report_null_coalesce", fixed = TRUE)
    testthat::expect_false(grepl(
      "(?m)^`%\\|\\|%`\\s*<-\\s*function",
      text,
      perl = TRUE
    ))
  }
})

testthat::test_that("runner and HPC shell helpers have one canonical util implementation", {
  canonical <- file.path(
    stage5_map_root,
    "util",
    "o2_supply_demand_map_shell_utils.sh"
  )
  testthat::expect_true(file.exists(canonical))
  testthat::expect_identical(
    unname(system2("bash", c("-n", canonical), stdout = TRUE, stderr = TRUE)),
    character()
  )

  expected_functions <- c(
    "truthy",
    "is_null_value",
    "normalize_fitting_mode",
    "require_positive_int",
    "require_nonnegative_int",
    "log_msg",
    "print_command",
    "run_or_print",
    "shell_join",
    "load_r_module",
    "derive_joint_warmup_seed_label",
    "label_joint_run_prefix",
    "resolve_existing_dir",
    "first_line",
    "o2sd_prov_shell_join",
    "o2sd_prov_cell",
    "o2sd_prov_init",
    "o2sd_prov_append",
    "o2sd_prov_write_many",
    "o2sd_prov_write_text",
    "o2sd_prov_write_standard",
    "o2sd_prov_record_sbatch"
  )
  declare_command <- paste(
    "source",
    shQuote(canonical),
    "; declare -F",
    paste(expected_functions, collapse = " ")
  )
  declared <- system2(
    "bash",
    c("-c", shQuote(declare_command)),
    stdout = TRUE,
    stderr = TRUE
  )
  testthat::expect_identical(attr(declared, "status"), NULL)
  declared_names <- sub("^.* ", "", declared)
  testthat::expect_setequal(declared_names, expected_functions)

  consumers <- c(
    file.path(stage5_map_root, "runner", "run_multi_warmup_joint.sh"),
    file.path(stage5_map_root, "runner", "run_o2_fit.sh"),
    file.path(stage5_map_root, "runner", "run_fit_joint_model_O2_supply_demand_MAP.sh"),
    file.path(stage5_map_root, "runner", "run_fit_model_O2_supply_demand_MAP.sh"),
    file.path(stage5_map_root, "hpc", "warm_up_joint_fitting_results_extra", "run_warm_up_joint_fitting_results_extra_hpc.sh"),
    file.path(stage5_map_root, "hpc", "warm_up_joint_fitting_results_extra", "submit_warm_up_joint_curve_array_hpc.sh"),
    file.path(stage5_map_root, "hpc", "submit", "submit_multi_warmup_task_table.sh"),
    file.path(stage5_map_root, "hpc", "submit", "submit_o2_fit.sh"),
    file.path(stage5_map_root, "hpc", "submit", "submit_multi_warmup_joint.sh"),
    file.path(stage5_map_root, "hpc", "fix_o2_simulation", "submit_fix_o2_simulation_array.sh"),
    file.path(stage5_map_root, "hpc", "best_fit_parameter_feature", "submit_best_fit_parameter_feature.sh"),
    file.path(stage5_map_root, "hpc", "fixo2_eigen_attractor", "submit_fixo2_eigen_attractor_array.sh"),
    file.path(stage5_map_root, "hpc", "postprocess", "submit_rerender_joint_seed_reports.sh"),
    file.path(stage5_map_root, "hpc", "postprocess", "postprocess_extra_results.sh"),
    file.path(stage5_map_root, "hpc", "array_workers", "run_fix_o2_simulation_array.sub"),
    file.path(stage5_map_root, "hpc", "array_workers", "run_multi_warmup_task_table_array.sub"),
    file.path(stage5_map_root, "hpc", "array_workers", "run_landscape_seed_space_task.sub"),
    file.path(stage5_map_root, "hpc", "parameter_landscape", "submit_parameter_landscape_full.sh"),
    file.path(stage5_map_root, "hpc", "dense_grid_monotonicity_classification", "submit_dense_grid_monotonicity_classification.sh")
  )
  testthat::expect_true(all(file.exists(consumers)))
  for (path in consumers) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_match(
      text,
      "o2_supply_demand_map_shell_utils.sh",
      fixed = TRUE,
      info = path
    )
    syntax <- system2("bash", c("-n", path), stdout = TRUE, stderr = TRUE)
    testthat::expect_identical(attr(syntax, "status"), NULL, info = path)
  }

  canonicalized <- list(
    truthy = consumers[c(1:16, 18:19)],
    is_null_value = consumers[c(1:2, 7:9)],
    require_nonnegative_int = consumers[c(1:2, 8:9)],
    log_msg = consumers[c(1, 7, 9)],
    print_command = consumers[c(1:2, 8:9)],
    run_or_print = consumers[c(1:2, 9)],
    shell_join = consumers[c(3:11, 14, 19)],
    load_r_module = consumers[c(2, 5, 7:9, 13:14, 16:19)],
    normalize_fitting_mode = consumers[c(2, 8)],
    require_positive_int = consumers[c(2, 8)],
    derive_joint_warmup_seed_label = consumers[c(2, 8)],
    label_joint_run_prefix = consumers[c(2, 8)],
    resolve_existing_dir = consumers[c(2, 8)],
    first_line = consumers[c(2, 8)]
  )
  for (name in names(canonicalized)) {
    pattern <- paste0("(?m)^", name, "\\s*\\(\\)\\s*\\{")
    for (path in canonicalized[[name]]) {
      text <- paste(readLines(path, warn = FALSE), collapse = "\n")
      testthat::expect_false(
        grepl(pattern, text, perl = TRUE),
        info = paste(name, path)
      )
    }
  }

  provenance_wrapper <- file.path(
    stage5_map_root,
    "hpc",
    "util",
    "write_run_provenance.sh"
  )
  provenance_text <- paste(readLines(provenance_wrapper, warn = FALSE), collapse = "\n")
  testthat::expect_lt(length(readLines(provenance_wrapper, warn = FALSE)), 15L)
  testthat::expect_match(provenance_text, basename(canonical), fixed = TRUE)
  testthat::expect_false(grepl("(?m)^o2sd_prov_[a-z_]+\\s*\\(\\)\\s*\\{", provenance_text, perl = TRUE))

  provenance_consumers <- list.files(
    file.path(stage5_map_root, "hpc"),
    pattern = "[.](sh|sub)$",
    recursive = TRUE,
    full.names = TRUE
  )
  provenance_consumer_text <- paste(
    unlist(lapply(provenance_consumers, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  testthat::expect_false(grepl("hpc/util/write_run_provenance.sh", provenance_consumer_text, fixed = TRUE))
  testthat::expect_match(provenance_consumer_text, basename(canonical), fixed = TRUE)
})
