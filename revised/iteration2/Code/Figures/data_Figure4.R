#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

data_Figure4 <- function(
    n_core = 8L,
    recompute_fixed_o2 = FALSE,
    recompute_tsne = FALSE
) {
  data_dir <- file.path(DATA_ROOT, "Figure4")
  viz_dir <- file.path(data_dir, "seed25_viz")
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

  seed25 <- file.path(INVIVO_RESULT_ROOT, "seed25")
  raw_files <- file.path(
    GEMCITABINE_DATA_ROOT,
    c(
      "dt_Gem_VT_20260209_v5.xlsx",
      "all_ploidy.csv",
      "histology_to_dt_Gem_VT_20260209_v5_mapping.csv"
    )
  )
  necrosis_source <- raw_files[[3L]]
  fit_files <- file.path(
    seed25,
    c("best_params.tsv", "fit_config.rds", "fit_summary.tsv")
  )
  require_files(c(raw_files, fit_files), "Figure 4A approved input")

  source_config_dir <- file.path(data_dir, "source_fit_config")
  runtime_fit_dir <- file.path(data_dir, "runtime_fit", "seed25")
  dir.create(source_config_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(runtime_fit_dir, recursive = TRUE, showWarnings = FALSE)

  source_fit_config <- fit_files[[2L]]
  source_fit_config_md5 <- unname(
    tools::md5sum(source_fit_config)[[1L]]
  )
  source_fit_config_copy <- copy_input(
    source_fit_config,
    file.path(source_config_dir, "fit_config.rds")
  )
  runtime_best_params <- copy_input(
    fit_files[[1L]],
    file.path(runtime_fit_dir, "best_params.tsv")
  )
  runtime_fit_summary <- copy_input(
    fit_files[[3L]],
    file.path(runtime_fit_dir, "fit_summary.tsv")
  )
  runtime_fit_config <- file.path(runtime_fit_dir, "fit_config.rds")
  if (!file.copy(
    source_fit_config,
    runtime_fit_config,
    overwrite = TRUE,
    copy.mode = TRUE
  )) {
    stop("Failed to create Figure 4A runtime fit configuration.")
  }

  runtime_cfg <- readRDS(runtime_fit_config)
  original_necrosis_mapping_csv <- as.character(
    runtime_cfg$necrosis_mapping_csv
  )
  runtime_cfg$necrosis_mapping_csv <- normalizePath(
    necrosis_source, mustWork = TRUE
  )
  saveRDS(runtime_cfg, runtime_fit_config)
  runtime_cfg_check <- readRDS(runtime_fit_config)
  if (!isTRUE(runtime_cfg_check$use_necrosis_loss)) {
    stop("Figure 4A runtime configuration does not enable necrosis loss.")
  }
  if (!identical(
    normalizePath(runtime_cfg_check$necrosis_mapping_csv, mustWork = TRUE),
    normalizePath(necrosis_source, mustWork = TRUE)
  )) {
    stop("Figure 4A runtime necrosis mapping was not applied.")
  }

  necrosis_table <- utils::read.csv(
    necrosis_source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_necrosis_columns <- c("dt_harvest", "percent_necrosis")
  missing_necrosis_columns <- setdiff(
    required_necrosis_columns, names(necrosis_table)
  )
  if (length(missing_necrosis_columns)) {
    stop(
      "Figure 4A necrosis mapping is missing column(s): ",
      paste(missing_necrosis_columns, collapse = ", ")
    )
  }
  necrosis_harvest <- trimws(as.character(necrosis_table$dt_harvest))
  necrosis_percent <- suppressWarnings(
    as.numeric(necrosis_table$percent_necrosis)
  )
  necrosis_keep <- nzchar(necrosis_harvest) & is.finite(necrosis_percent)
  if ("mapping_status" %in% names(necrosis_table)) {
    necrosis_keep <- necrosis_keep &
      tolower(trimws(as.character(necrosis_table$mapping_status))) == "mapped"
  }
  if (any(
    necrosis_keep & (necrosis_percent < 0 | necrosis_percent > 100),
    na.rm = TRUE
  )) {
    stop("Figure 4A percent_necrosis contains a value outside [0, 100].")
  }
  n_necrosis_rows <- sum(necrosis_keep)
  n_necrosis_harvests <- length(unique(necrosis_harvest[necrosis_keep]))
  if (n_necrosis_rows != 7L || n_necrosis_harvests != 7L) {
    stop(
      "Figure 4A expected 7 usable necrosis rows across 7 harvests; found ",
      n_necrosis_rows, " rows across ", n_necrosis_harvests, " harvests."
    )
  }

  run_process(
    "Rscript",
    args = c(
      file.path(
        MODEL_CODE_ROOT, "vis",
        "viz_invivo_model_O2_supply_demand_MAP_results.R"
      ),
      paste0("--fit_dir=", runtime_fit_dir),
      paste0("--data_dir=", GEMCITABINE_DATA_ROOT),
      paste0("--out_dir=", viz_dir),
      "--report_dt=1",
      "--n_cores=1"
    )
  )

  runtime_provenance <- data.frame(
    source_fit_config = normalizePath(source_fit_config, mustWork = TRUE),
    source_fit_config_md5 = source_fit_config_md5,
    source_fit_config_copy = normalizePath(
      source_fit_config_copy, mustWork = TRUE
    ),
    original_necrosis_mapping_csv = original_necrosis_mapping_csv,
    runtime_fit_config = normalizePath(runtime_fit_config, mustWork = TRUE),
    runtime_fit_config_md5 = unname(
      tools::md5sum(runtime_fit_config)[[1L]]
    ),
    runtime_necrosis_mapping_csv = normalizePath(
      runtime_cfg_check$necrosis_mapping_csv, mustWork = TRUE
    ),
    necrosis_source = normalizePath(necrosis_source, mustWork = TRUE),
    necrosis_source_md5 = unname(
      tools::md5sum(necrosis_source)[[1L]]
    ),
    use_necrosis_loss = isTRUE(runtime_cfg_check$use_necrosis_loss),
    usable_necrosis_rows = n_necrosis_rows,
    usable_necrosis_harvests = n_necrosis_harvests,
    stringsAsFactors = FALSE
  )
  write_intermediate_tsv(
    runtime_provenance,
    file.path(data_dir, "figure4a_runtime_config_provenance.tsv")
  )

  viz_map <- c(
    burden_timecourse.tsv = "figure4a_seed25_burden_timecourse.tsv",
    ploidy_weighted_mean_timecourse.tsv =
      "figure4a_seed25_weighted_mean_chromosome_timecourse.tsv",
    terminal_ploidy_observed_vs_predicted.tsv =
      "figure4a_seed25_terminal_chromosome_observed_vs_predicted.tsv"
  )
  generated <- file.path(viz_dir, names(viz_map))
  require_files(generated, "Figure 4A regenerated table")
  destinations <- file.path(data_dir, unname(viz_map))
  for (i in seq_along(generated)) {
    ok <- file.copy(generated[[i]], destinations[[i]], overwrite = TRUE)
    if (!ok) stop("Failed to stage Figure 4A intermediate: ", generated[[i]])
  }

  fixed_stage <- file.path(data_dir, "fixed_o2_analysis")
  fixed_data <- file.path(fixed_stage, "data")
  run_fixed_o2_auc(
    args = c(
      paste0("--fit_root=", INVIVO_RESULT_ROOT),
      paste0("--output_dir=", fixed_stage),
      paste0(
        "--fixed_o2_script=",
        file.path(MODEL_CODE_ROOT, "simulation", "fix_o2_simulation.R")
      ),
      paste0("--n_core=", as.integer(n_core)),
      paste0("--overwrite=", if (isTRUE(recompute_fixed_o2)) "TRUE" else "FALSE"),
      "--analysis_only=TRUE"
    ),
    env = paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT)
  )
  fixed_outputs <- c(
    "invivo_best_parameters_500seeds.tsv",
    "parameter_auc_by_o2.tsv",
    "fixed_o2_dominant_ploidy_201grid.tsv",
    "ploidy_class_counts_by_o2.tsv",
    "validation_summary.tsv",
    "run_arguments.tsv",
    "source_file_provenance.tsv"
  )
  require_files(
    file.path(fixed_data, fixed_outputs),
    "Figure 4 fixed-O2 intermediate"
  )
  for (name in fixed_outputs) {
    ok <- file.copy(
      file.path(fixed_data, name),
      file.path(data_dir, name),
      overwrite = TRUE
    )
    if (!ok) stop("Failed to stage fixed-O2 intermediate: ", name)
  }

  parameter_table_source <- file.path(
    INVIVO_RESULT_ROOT, "seed25", "parameter_table.csv"
  )
  parameter_table_local <- copy_input(
    parameter_table_source,
    file.path(data_dir, "invivo_parameter_table_seed25.csv")
  )
  config_sources <- file.path(
    AUDIT_ROOT, "parameters",
    c("parameter_function_groups.tsv", "parameter_function_group_palette.tsv")
  )
  for (path in config_sources) {
    ok <- file.copy(path, file.path(data_dir, basename(path)), overwrite = TRUE)
    if (!ok) stop("Failed to stage Figure 4 plotting configuration: ", path)
  }

  run_invivo_tsne_landscape(
    args = c(
      paste0("--fit_root=", INVIVO_RESULT_ROOT),
      paste0("--data_dir=", data_dir),
      paste0("--n_core=", as.integer(n_core)),
      paste0("--overwrite=", if (isTRUE(recompute_tsne)) "TRUE" else "FALSE")
    ),
    env = paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT)
  )
  tsne_outputs <- c(
    "invivo_initial_population_18params.tsv.gz",
    "invivo_initial_best_tsne_coordinates.tsv.gz",
    "invivo_best_tsne_cluster_coordinates.tsv",
    "invivo_tsne_cluster_summary.tsv",
    "invivo_tsne_clustering_diagnostics.tsv",
    "invivo_tsne_preprocessing_metadata.tsv",
    "invivo_tsne_run_metadata.tsv",
    "invivo_tsne_source_provenance.tsv"
  )
  require_files(
    file.path(data_dir, tsne_outputs),
    "Figure 4B in vivo-only 18-parameter t-SNE output"
  )

  source(resolve_utility_file("analysis/optimizer_diagnostics.R"))
  objective_metrics <- optimizer_collect_seed_metrics(
    INVIVO_RESULT_ROOT,
    objective_metric = "objective"
  )
  objective_ranking <- objective_metrics[, c(
    "seed", "objective_rank", "objective"
  ), drop = FALSE]
  objective_ranking$seed_number <- as.integer(sub(
    "^seed", "", objective_ranking$seed
  ))
  objective_ranking <- objective_ranking[, c(
    "seed_number", "seed", "objective_rank", "objective"
  )]
  if (nrow(objective_ranking) != 500L ||
      objective_ranking$seed_number[[1L]] != 25L ||
      objective_ranking$objective_rank[[1L]] != 1L) {
    stop("Figure 4B requires seed25 to be rank 1 by canonical fit objective.")
  }
  write_intermediate_tsv(
    objective_ranking,
    file.path(data_dir, "invivo_fit_objective_ranking_500seeds.tsv")
  )

  source(resolve_utility_file(
    "analysis/figure4_continuous_ploidy_association.R"
  ))
  derive_figure4_continuous_ploidy_association(data_dir = data_dir)
  continuous_outputs <- c(
    "figure4a_burden_measurement_audit.tsv",
    "invivo_fit_objective_ranking_500seeds.tsv",
    "continuous_ploidy_spearman_by_o2.tsv",
    "continuous_ploidy_parameter_ranking.tsv",
    "all_parameter_fitted_endpoint_values.tsv",
    "all_parameter_pooled_distribution_summary.tsv",
    "all_parameter_log10_range_summary.tsv",
    "all_parameter_log10_density.tsv",
    "parameter_display_labels.tsv",
    "continuous_ploidy_analysis_validation.tsv",
    "continuous_ploidy_analysis_source_provenance.tsv"
  )
  require_files(
    file.path(data_dir, continuous_outputs),
    "Figure 4 continuous-ploidy output"
  )

  fixed_seed_dirs <- sort(list.dirs(
    INVIVO_RESULT_ROOT, recursive = FALSE, full.names = TRUE
  ))
  fixed_seed_dirs <- fixed_seed_dirs[grepl("^seed[0-9]+$", basename(fixed_seed_dirs))]
  fixed_seed_inputs <- c(
    file.path(fixed_seed_dirs, "fit_summary.tsv"),
    file.path(fixed_seed_dirs, "best_params.tsv"),
    file.path(fixed_seed_dirs[[1L]], "fit_config.rds")
  )
  tsne_seed_inputs <- c(
    file.path(fixed_seed_dirs, "fit_config.rds"),
    file.path(fixed_seed_dirs, "parameter_table.csv"),
    file.path(fixed_seed_dirs, "fit_status.log")
  )
  require_files(
    c(fixed_seed_inputs, tsne_seed_inputs),
    "Figure 4 separate in vivo seed input"
  )
  source_rows <- c(
    raw_files,
    fit_files,
    parameter_table_source,
    fixed_seed_inputs,
    tsne_seed_inputs
  )
  local_rows <- c(
    rep(NA_character_, length(raw_files)),
    runtime_best_params,
    source_fit_config_copy,
    runtime_fit_summary,
    parameter_table_local,
    rep(
      NA_character_,
      length(fixed_seed_inputs) + length(tsne_seed_inputs)
    )
  )
  contract <- data.frame(
    role = c(
      rep("Figure 4A raw/fit input", length(raw_files) + length(fit_files)),
      "Figure 4B parameter-bounds input",
      rep("Figure 4B fixed-O2 seed input", length(fixed_seed_inputs)),
      rep(
        "Figure 4B in vivo-only t-SNE initialization input",
        length(tsne_seed_inputs)
      )
    ),
    source = normalizePath(source_rows, mustWork = TRUE),
    local_file = local_rows,
    source_md5 = unname(tools::md5sum(source_rows)),
    local_md5 = c(
      rep(NA_character_, length(raw_files)),
      unname(tools::md5sum(c(
        runtime_best_params,
        source_fit_config_copy,
        runtime_fit_summary
      ))),
      unname(tools::md5sum(parameter_table_local)),
      rep(
        NA_character_,
        length(fixed_seed_inputs) + length(tsne_seed_inputs)
      )
    ),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure4", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_core_arg <- sub("^--n-core=", "", args[grepl("^--n-core=", args)])
  recompute_arg <- sub(
    "^--recompute-fixed-o2=", "",
    args[grepl("^--recompute-fixed-o2=", args)]
  )
  recompute_tsne_arg <- sub(
    "^--recompute-tsne=", "",
    args[grepl("^--recompute-tsne=", args)]
  )
  data_Figure4(
    n_core = if (length(n_core_arg)) as.integer(n_core_arg[[1L]]) else 8L,
    recompute_fixed_o2 = if (length(recompute_arg)) {
      as_boolean(recompute_arg[[1L]])
    } else {
      FALSE
    },
    recompute_tsne = if (length(recompute_tsne_arg)) {
      as_boolean(recompute_tsne_arg[[1L]])
    } else {
      FALSE
    }
  )
}
