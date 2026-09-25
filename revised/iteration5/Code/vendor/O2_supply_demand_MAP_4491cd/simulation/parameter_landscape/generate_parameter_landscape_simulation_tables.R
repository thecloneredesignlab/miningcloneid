#!/usr/bin/env Rscript

.o2pl_generator_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "generate_parameter_landscape_simulation_tables.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(.o2pl_generator_dir, "parameter_landscape_simulation_utils.R"), local = environment(), chdir = TRUE)

o2pl_materialize_growth_turnover <- function(input_dir,
                                             root_dir,
                                             data_dir = NULL,
                                             max_seeds = NA_integer_,
                                             horizon_day = 100,
                                             report_dt = 1) {
  feature_file <- file.path(.o2pl_generator_dir, "parameter_landscape_invivo_feature_simulation.R")
  sys.source(feature_file, envir = environment(), chdir = TRUE)
  canonical_dir <- o2pl_simulation_tables_dir(root_dir, "invivo")
  canonical_csv <- file.path(canonical_dir, "growth_turnover_100d.csv")
  best_input <- file.path(canonical_dir, "fitted_parameter_features.tsv")
  if (!file.exists(best_input)) stop("Missing materialized fitted-parameter table: ", best_input, call. = FALSE)
  best_csv <- file.path(canonical_dir, "fitted_parameter_features.csv")
  write_csv(read_tsv(best_input), best_csv)
  paper_generate_invivo_growth_turnover_table(
    input_dir = input_dir,
    best_csv = best_csv,
    output_csv = canonical_csv,
    data_dir = data_dir %||% default_invivo_data_dir(),
    max_seeds = max_seeds,
    horizon_day = horizon_day,
    report_dt = report_dt
  )
  data <- read_csv_plain(canonical_csv)
  canonical_tsv <- file.path(canonical_dir, "growth_turnover_100d.tsv")
  schema <- file.path(canonical_dir, "growth_turnover_100d.schema.tsv")
  legacy <- file.path(paper_tables_dir("invivo", root_dir), "invivo_best_params_growth_turnover_100d.csv")
  write_tsv_plain(data, canonical_tsv)
  write_csv(data, legacy)
  o2pl_write_schema(data, "invivo_growth_turnover_100d", schema)
  o2pl_record_artifact(root_dir, "simulation", "invivo_growth_turnover_100d", canonical_tsv, data, schema, basename(sys.frame(1)$ofile %||% "generate_parameter_landscape_simulation_tables.R"), "invivo", input_dir)
  invisible(canonical_tsv)
}

o2pl_simulation_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  dataset <- tolower(trimws(as.character(argv$dataset %||% argv$analysis_part %||% "both")))
  dataset <- switch(dataset,
    invivo_tables = "invivo",
    invitro_tables = "invitro",
    all = "both",
    dataset
  )
  if (!dataset %in% c("invivo", "invitro", "both")) stop("--dataset must be invivo, invitro, or both.", call. = FALSE)
  root_dir <- normalizePath(path.expand(argv$result_root %||% o2pl_default_result_root()), mustWork = FALSE)
  max_seeds <- as_int(argv$max_seeds, NA_integer_)
  datasets <- if (dataset == "both") c("invivo", "invitro") else dataset
  outputs <- list()
  for (current in datasets) {
    input_dir <- if (current == "invivo") argv$invivo_input %||% argv$input_dir %||% DEFAULT_INVIVO_INPUT_DIR else argv$invitro_input %||% argv$input_dir %||% DEFAULT_INVITRO_INPUT_DIR
    outputs[[current]] <- o2pl_materialize_seed_parameter_tables(
      dataset = current,
      input_dir = input_dir,
      root_dir = root_dir,
      max_seeds = max_seeds,
      write_best = as_bool(argv$write_best, TRUE),
      write_initial = as_bool(argv$write_initial, TRUE)
    )
    if (current == "invivo" && as_bool(argv$write_growth_turnover, FALSE)) {
      outputs$growth_turnover <- o2pl_materialize_growth_turnover(
        input_dir = input_dir,
        root_dir = root_dir,
        data_dir = argv$data_dir,
        max_seeds = max_seeds,
        horizon_day = as_num(argv$growth_turnover_horizon_day, 100),
        report_dt = as_num(argv$growth_turnover_report_dt, 1)
      )
    }
    if (current == "invivo" && as_bool(argv$write_modes, TRUE)) {
      outputs$fixed_o2_modes <- o2pl_materialize_fixed_o2_modes(
        input_dir = input_dir,
        root_dir = root_dir,
        attractor_o2_grid = as_num_vec(argv$attractor_o2_grid, paper_default_attractor_o2_grid()),
        summary_o2 = as_num_vec(argv$summary_o2, paper_default_mode_summary_o2()),
        reference_o2 = as_num_vec(argv$mode_reference_o2_values, as_num(argv$mode_reference_o2, 2)),
        max_seeds = max_seeds,
        n_workers = as_int(argv$n_workers, 1L),
        overwrite = as_bool(argv$overwrite_modes %||% argv$force_modes, FALSE)
      )
    }
  }
  args_table <- data.frame(
    argument = names(argv),
    value = vapply(argv, function(x) paste(x, collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  args_path <- file.path(o2pl_simulation_tables_dir(root_dir), "run_arguments.tsv")
  write_tsv_plain(args_table, args_path)
  o2pl_record_artifact(root_dir, "simulation", "parameter_landscape_simulation_run_arguments", args_path, args_table, NA_character_, "generate_parameter_landscape_simulation_tables.R", dataset, paste(datasets, collapse = ","))
  message("Parameter-landscape simulation tables complete: ", root_dir)
  invisible(outputs)
}

if (sys.nframe() == 0L) o2pl_simulation_main()
