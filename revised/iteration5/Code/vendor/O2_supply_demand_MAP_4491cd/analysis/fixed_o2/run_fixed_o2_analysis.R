#!/usr/bin/env Rscript

# Table-only fixed-O2 analyses. All numerical trajectories, attractors, and
# fitted-parameter values must already have been materialized by an upstream
# producer before this entrypoint is run.

.fixo2_analysis_script_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      path <- env$ofile
      if (is.null(path)) "" else normalizePath(path, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[basename(frame_files) == "run_fixed_o2_analysis.R"]
  if (length(own)) {
    dirname(own[[length(own)]])
  } else {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    own_file_arg <- file_arg[
      basename(sub("^--file=", "", file_arg)) == "run_fixed_o2_analysis.R"
    ]
    if (length(own_file_arg)) {
      dirname(normalizePath(sub("^--file=", "", own_file_arg[[1L]]), mustWork = FALSE))
    } else {
      normalizePath(getwd(), mustWork = FALSE)
    }
  }
})

.fixo2_analysis_workflow_root <- normalizePath(
  file.path(.fixo2_analysis_script_dir, "..", ".."),
  mustWork = TRUE
)
sys.source(
  file.path(
    .fixo2_analysis_workflow_root,
    "util",
    "o2_supply_demand_map_fixed_o2_utils.R"
  ),
  envir = environment(),
  chdir = TRUE
)
sys.source(
  file.path(
    .fixo2_analysis_workflow_root,
    "util",
    "o2_supply_demand_map_fixed_o2_table_utils.R"
  ),
  envir = environment(),
  chdir = TRUE
)

fixo2_analysis_parse_args <- o2ipa_parse_args

fixo2_analysis_usage <- function() {
  paste(
    "Usage:",
    "  Rscript run_fixed_o2_analysis.R --simulation_dir DIR [options]",
    "",
    "Required:",
    "  --simulation_dir DIR  Root containing materialized fixed-O2 tables.",
    "",
    "Options:",
    "  --analysis_dir DIR    Output root (default: DIR/analysis).",
    "  --parts CSV           attractors,counterfactual,agreement (default: all).",
    "  --mode_reference_o2 N Fixed-O2 value used for seed mode assignment (default: 2).",
    "  --help                Print this message.",
    sep = "\n"
  )
}



sys.source(
  file.path(.fixo2_analysis_script_dir, "fixed_o2_analysis_utils.R"),
  envir = environment(),
  chdir = TRUE
)

fixo2_analysis_parts <- function(x) {
  values <- tolower(trimws(strsplit(as.character(x %||% "all"), ",", fixed = TRUE)[[1L]]))
  if ("all" %in% values) {
    return(c("attractors", "counterfactual", "agreement"))
  }
  values <- sub("^counterfactual_trajectories$", "counterfactual", values)
  allowed <- c("attractors", "counterfactual", "agreement")
  unknown <- setdiff(values, allowed)
  if (length(unknown)) {
    stop("Unknown --parts value(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  unique(values)
}

fixo2_analysis_find_input <- function(root, candidates, required = TRUE) {
  paths <- normalizePath(file.path(root, candidates), mustWork = FALSE)
  hit <- paths[file.exists(paths)]
  if (length(hit)) return(hit[[1L]])
  if (required) {
    stop(
      "Required materialized table was not found. Checked: ",
      paste(paths, collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

fixo2_analysis_validate_columns <- function(tab, required, label) {
  missing <- setdiff(required, names(tab))
  if (length(missing)) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

fixo2_analysis_attractors <- function(simulation_dir, analysis_dir, mode_reference_o2) {
  source_path <- fixo2_analysis_find_input(
    simulation_dir,
    c(
      "attractors/tables/fixed_o2_attractors_by_seed.tsv",
      "tables/fixed_o2_attractors_by_seed.tsv"
    )
  )
  params_path <- fixo2_analysis_find_input(
    simulation_dir,
    c(
      "attractors/tables/parameter_values_long.tsv",
      "tables/parameter_values_long.tsv"
    ),
    required = FALSE
  )
  out_dir <- file.path(analysis_dir, "attractors")
  table_dir <- file.path(out_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  attractors <- fo2_read_tsv(source_path)
  fixo2_analysis_validate_columns(
    attractors,
    c("seed_id", "O2_pct", "dominant_mean_ploidy"),
    "Fixed-O2 attractor table"
  )
  attractors <- fixo2_assign_attractor_modes(
    attractors,
    ploidy_col = "dominant_mean_ploidy"
  )
  mode_by_seed_o2 <- fixo2_attractor_mode_table(attractors)
  mode_reference_by_seed <- fixo2_reference_mode_table(
    mode_by_seed_o2,
    mode_reference_o2
  )
  mode_summary <- fixo2_attractor_mode_summary_by_seed(mode_by_seed_o2)
  regime_summary <- fo2_attractor_regime_summary(attractors)
  regime_tests <- fo2_attractor_regime_tests(attractors)
  gap_by_seed <- fo2_spectral_gap_by_seed(attractors)
  gap_summary <- fo2_spectral_gap_summary(gap_by_seed)

  fo2_write_tsv(
    mode_by_seed_o2,
    file.path(table_dir, "fixed_o2_attractor_mode_by_seed_o2.tsv")
  )
  fo2_write_tsv(
    mode_reference_by_seed,
    file.path(table_dir, "fixed_o2_attractor_mode_reference_by_seed.tsv")
  )
  fo2_write_tsv(
    mode_reference_by_seed,
    file.path(table_dir, "fixed_o2_attractor_mode_by_seed.tsv")
  )
  fo2_write_tsv(
    mode_summary,
    file.path(table_dir, "fixed_o2_attractor_mode_summary_by_seed.tsv")
  )
  fo2_write_tsv(
    fo2_mode_seed_stack_table(attractors),
    file.path(table_dir, "fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv")
  )
  fo2_write_mode_comparison_tables(attractors, out_dir)
  fo2_write_tsv(
    regime_summary,
    file.path(table_dir, "fixed_o2_attractor_regime_summary.tsv")
  )
  fo2_write_tsv(
    regime_tests,
    file.path(table_dir, "fixed_o2_attractor_regime_tests.tsv")
  )
  fo2_write_tsv(
    gap_by_seed,
    file.path(table_dir, "fixed_o2_attractor_spectral_gap_by_seed.tsv")
  )
  fo2_write_tsv(
    gap_summary,
    file.path(table_dir, "fixed_o2_attractor_spectral_gap_regime_summary.tsv")
  )
  fo2_write_tsv(
    fo2_ploidy_gap_reliability_composite_table(gap_by_seed),
    file.path(table_dir, "fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv")
  )
  fo2_write_tsv(
    fo2_ploidy_gap_reliability_violin_table(gap_by_seed),
    file.path(table_dir, "fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv")
  )

  if (!is.na(params_path)) {
    params <- fo2_read_tsv(params_path)
    correlations <- fo2_parameter_correlations(attractors, params)
    if (nrow(correlations)) {
      correlations <- correlations[
        order(
          correlations$O2_pct,
          correlations$metric,
          -correlations$abs_spearman_rho
        ),
        ,
        drop = FALSE
      ]
    }
    fo2_write_tsv(
      correlations,
      file.path(table_dir, "parameter_attractor_correlations.tsv")
    )
  }

  data.frame(
    part = "attractors",
    input_table = source_path,
    output_dir = normalizePath(out_dir, mustWork = FALSE),
    input_rows = nrow(attractors),
    input_seeds = length(unique(attractors$seed_id)),
    stringsAsFactors = FALSE
  )
}

fixo2_analysis_counterfactual <- function(simulation_dir, analysis_dir) {
  source_path <- fixo2_analysis_find_input(
    simulation_dir,
    c(
      "counterfactual_trajectories/tables/fixed_o2_counterfactual_summary_by_seed.tsv",
      "counterfactual/tables/fixed_o2_counterfactual_summary_by_seed.tsv",
      "tables/fixed_o2_counterfactual_summary_by_seed.tsv"
    )
  )
  params_path <- fixo2_analysis_find_input(
    simulation_dir,
    c(
      "counterfactual_trajectories/tables/parameter_values_long.tsv",
      "counterfactual/tables/parameter_values_long.tsv",
      "attractors/tables/parameter_values_long.tsv",
      "tables/parameter_values_long.tsv"
    ),
    required = FALSE
  )
  out_dir <- file.path(analysis_dir, "counterfactual_trajectories")
  table_dir <- file.path(out_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  summary_by_seed <- fo2_read_tsv(source_path)
  fixo2_analysis_validate_columns(
    summary_by_seed,
    c("seed_id", "O2_pct", "trajectory_regime", "initial_condition"),
    "Fixed-O2 counterfactual summary table"
  )
  regime_summary <- cf2_regime_summary(summary_by_seed)
  regime_tests <- cf2_regime_tests(summary_by_seed)
  fo2_write_tsv(
    regime_summary,
    file.path(table_dir, "fixed_o2_counterfactual_regime_summary.tsv")
  )
  fo2_write_tsv(
    regime_tests,
    file.path(table_dir, "fixed_o2_counterfactual_regime_tests.tsv")
  )
  if (!is.na(params_path)) {
    params <- fo2_read_tsv(params_path)
    correlations <- cf2_parameter_correlations(summary_by_seed, params)
    fo2_write_tsv(
      correlations,
      file.path(table_dir, "fixed_o2_counterfactual_parameter_correlations.tsv")
    )
  }

  data.frame(
    part = "counterfactual",
    input_table = source_path,
    output_dir = normalizePath(out_dir, mustWork = FALSE),
    input_rows = nrow(summary_by_seed),
    input_seeds = length(unique(summary_by_seed$seed_id)),
    stringsAsFactors = FALSE
  )
}

fixo2_analysis_agreement <- function(simulation_dir, analysis_dir) {
  source_path <- fixo2_analysis_find_input(
    simulation_dir,
    c(
      "analytical_simulation_agreement/tables/agreement_analytical_vs_simulation_data.tsv",
      "simulation/analytical_simulation_agreement/tables/agreement_analytical_vs_simulation_data.tsv",
      "tables/agreement_analytical_vs_simulation_data.tsv"
    )
  )
  out_dir <- file.path(analysis_dir, "analytical_simulation_agreement")
  table_dir <- file.path(out_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  agreement <- fo2_read_tsv(source_path)
  agreement <- fixo2_add_agreement_limits(agreement)
  fo2_write_tsv(
    agreement,
    file.path(table_dir, "agreement_augmented_data.tsv")
  )
  fo2_write_tsv(
    fixo2_agreement_bland_altman_summary(agreement),
    file.path(table_dir, "agreement_bland_altman_summary.tsv")
  )
  data.frame(
    part = "agreement",
    input_table = source_path,
    output_dir = normalizePath(out_dir, mustWork = FALSE),
    input_rows = nrow(agreement),
    input_seeds = if ("seed_id" %in% names(agreement)) {
      length(unique(agreement$seed_id))
    } else {
      NA_integer_
    },
    stringsAsFactors = FALSE
  )
}

fixo2_analysis_main <- function(args = fixo2_analysis_parse_args()) {
  if (isTRUE(args$help)) {
    cat(fixo2_analysis_usage(), "\n")
    return(invisible(NULL))
  }
  simulation_dir <- as.character(args$simulation_dir %||% "")
  if (!nzchar(simulation_dir)) {
    stop("--simulation_dir is required.\n", fixo2_analysis_usage(), call. = FALSE)
  }
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  analysis_dir <- normalizePath(
    as.character(args$analysis_dir %||% file.path(simulation_dir, "analysis")),
    mustWork = FALSE
  )
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  parts <- fixo2_analysis_parts(args$parts)
  mode_reference_o2 <- o2ipa_as_num(args$mode_reference_o2, 2)
  rows <- list()
  if ("attractors" %in% parts) {
    rows[[length(rows) + 1L]] <- fixo2_analysis_attractors(
      simulation_dir,
      analysis_dir,
      mode_reference_o2
    )
  }
  if ("counterfactual" %in% parts) {
    rows[[length(rows) + 1L]] <- fixo2_analysis_counterfactual(
      simulation_dir,
      analysis_dir
    )
  }
  if ("agreement" %in% parts) {
    rows[[length(rows) + 1L]] <- fixo2_analysis_agreement(
      simulation_dir,
      analysis_dir
    )
  }
  manifest <- if (length(rows)) do.call(rbind, rows) else data.frame()
  fo2_write_tsv(manifest, file.path(analysis_dir, "fixed_o2_analysis_manifest.tsv"))
  message("Completed table-only fixed-O2 analyses: ", analysis_dir)
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  fixo2_analysis_main()
}
