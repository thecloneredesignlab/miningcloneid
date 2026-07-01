#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
source(file.path(SCRIPT_DIR, "parameter_landscape_utils.R"))

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript clustering_runner.R --run_parts=all [options]",
      "",
      "run_parts values:",
      "  all                         Run tables, reductions, and clustering report.",
      "  tables                      Run invivo_tables and invitro_tables.",
      "  invivo_tables               Write in vivo tables and fixed-O2 mode tables.",
      "  invitro_tables              Write in vitro tables.",
      "  reductions,umap             Run invivo_reductions, invitro_reductions, and pooled_reductions.",
      "  invivo_reductions           Generate in vivo reductions and clustered plots.",
      "  invitro_reductions          Generate in vitro reductions and clustered plots.",
      "  pooled_reductions           Generate pooled in vivo/in vitro reductions and distance tables.",
      "  reports,clustering_report   Generate only the clustering HTML reports.",
      "",
      "Common options:",
      "  --result_root=/path/to/oxygen/results/analysis/best_fit_parameter_feature/02_parameter_landscape_clustering/parameter_landscape",
      "  --invivo_input=/path/to/fit_invivo_O2_buffering_500seed",
      "  --invitro_input=/path/to/fit_invitro_O2_buffering_500seed",
      "  --mode_reference_o2=2",
      "  --mode_reference_o2_values=0,0.1,0.5,1,2,5",
      "  --attractor_o2_grid=0,0.1,0.5,1,2,5",
      "  --summary_o2=0,0.1,0.5,1,2,5",
      "  --attractor_feature_o2_values=0,0.1,0.5,1,2,5",
      "  --max_seeds=N",
      "  --n_threads=N",
      "  --run_report=TRUE|FALSE",
      "  --run_umap_cluster_report=TRUE|FALSE",
      "  --run_pca_cluster_report=TRUE|FALSE",
      "  --dry_run=TRUE|FALSE",
      sep = "\n"
    ),
    "\n"
  )
}

append_arg <- function(args, key, value) {
  if (is.null(value) || length(value) == 0L || all(is.na(value))) return(args)
  value <- as.character(value[[1L]])
  if (!nzchar(value)) return(args)
  c(args, paste0("--", key, "=", value))
}

append_args <- function(args, values) {
  for (key in names(values)) {
    args <- append_arg(args, key, values[[key]])
  }
  args
}

run_child_script <- function(label, script, args, dry_run = FALSE) {
  script_path <- file.path(SCRIPT_DIR, script)
  if (!file.exists(script_path)) stop("Missing child script: ", script_path)
  rscript <- file.path(R.home("bin"), "Rscript")
  cmd <- c("--vanilla", script_path, args)
  message("")
  message("[", format(Sys.time(), "%F %T"), "] ", label)
  message("Command: ", paste(shQuote(c(rscript, cmd)), collapse = " "))
  if (isTRUE(dry_run)) return(invisible(0L))
  status <- system2(rscript, cmd)
  if (!identical(status, 0L)) {
    stop("Child script failed with status ", status, ": ", script)
  }
  invisible(status)
}

normalize_run_parts <- function(parts) {
  raw <- tolower(as_char_vec(parts, "all"))
  aliases <- list(
    all = c(
      "invivo_tables", "invitro_tables",
      "invivo_reductions", "invitro_reductions", "pooled_reductions", "clustering_report"
    ),
    table = c("invivo_tables", "invitro_tables"),
    tables = c("invivo_tables", "invitro_tables"),
    umap = c("invivo_reductions", "invitro_reductions", "pooled_reductions"),
    reduction = c("invivo_reductions", "invitro_reductions", "pooled_reductions"),
    reductions = c("invivo_reductions", "invitro_reductions", "pooled_reductions"),
    report = "clustering_report",
    reports = "clustering_report",
    umap_report = "clustering_report"
  )
  expanded <- character()
  for (part in raw) {
    expanded <- c(expanded, aliases[[part]] %||% part)
  }
  valid <- c(
    "invivo_tables", "invitro_tables",
    "invivo_reductions", "invitro_reductions", "pooled_reductions",
    "clustering_report"
  )
  unknown <- setdiff(expanded, valid)
  if (length(unknown)) stop("Unknown run_parts value(s): ", paste(unknown, collapse = ", "))
  ordered <- c(
    "invivo_tables", "invitro_tables",
    "invivo_reductions", "invitro_reductions", "pooled_reductions",
    "clustering_report"
  )
  ordered[ordered %in% unique(expanded)]
}

common_max_seed_args <- function(argv) {
  args <- character()
  append_arg(args, "max_seeds", argv$max_seeds)
}

run_invivo_tables <- function(argv, result_root, invivo_input, dry_run) {
  args <- append_args(character(), list(
    analysis_part = "invivo_tables",
    input_dir = invivo_input,
    result_root = result_root,
    write_modes = "TRUE",
    overwrite_modes = argv$overwrite_modes %||% argv$force_modes %||% "FALSE",
    n_workers = argv$n_workers %||% argv$n_threads,
    mode_reference_o2 = argv$mode_reference_o2 %||% "2",
    mode_reference_o2_values = argv$mode_reference_o2_values %||% "0,0.1,0.5,1,2,5",
    attractor_o2_grid = argv$attractor_o2_grid,
    summary_o2 = argv$summary_o2
  ))
  args <- c(args, common_max_seed_args(argv))
  run_child_script("Write in vivo tables and fixed-O2 mode tables", "clustering_analysis.R", args, dry_run = dry_run)
}

run_invitro_tables <- function(argv, result_root, invitro_input, dry_run) {
  args <- append_args(character(), list(
    analysis_part = "invitro_tables",
    input_dir = invitro_input,
    result_root = result_root
  ))
  args <- c(args, common_max_seed_args(argv))
  run_child_script("Write in vitro tables", "clustering_analysis.R", args, dry_run = dry_run)
}

run_invivo_reductions <- function(argv, result_root, invivo_input, dry_run) {
  args <- append_args(character(), list(
    analysis_part = "invivo_reductions",
    result_root = result_root,
    objective_seed_dir = invivo_input,
    reductions = argv$reductions %||% "umap,pca,tsne",
    preprocess_modes = argv$preprocess_modes %||% "zscore,prior_unit",
    run_pca_umap = "TRUE",
    run_full_tsne = argv$run_full_tsne %||% "FALSE",
    mode_tables_dir = argv$mode_tables_dir %||% file.path(result_root, "FixO2Modes"),
    mode_reference_o2 = argv$mode_reference_o2 %||% "2",
    shape_by_mode = "TRUE",
    drop_parameter_table_initial = "TRUE",
    umap_seed = argv$umap_seed %||% "123",
    cluster_seed = argv$cluster_seed %||% "123",
    n_threads = argv$n_threads,
    attractor_feature_o2_values = argv$attractor_feature_o2_values
  ))
  run_child_script("Generate in vivo reductions", "clustering_analysis.R", args, dry_run = dry_run)
}

run_invitro_reductions <- function(argv, result_root, invitro_input, dry_run) {
  args <- append_args(character(), list(
    analysis_part = "invitro_reductions",
    result_root = result_root,
    objective_seed_dir = invitro_input,
    reductions = argv$reductions %||% "umap,pca,tsne",
    preprocess_modes = argv$preprocess_modes %||% "zscore,prior_unit",
    run_pca_umap = "TRUE",
    run_full_tsne = argv$run_full_tsne %||% "FALSE",
    drop_parameter_table_initial = "TRUE",
    umap_seed = argv$umap_seed %||% "123",
    cluster_seed = argv$cluster_seed %||% "123",
    n_threads = argv$n_threads
  ))
  run_child_script("Generate in vitro reductions", "clustering_analysis.R", args, dry_run = dry_run)
}

run_pooled_reductions <- function(argv, result_root, invivo_input, invitro_input, dry_run) {
  args <- append_args(character(), list(
    analysis_part = "pooled_reductions",
    result_root = result_root,
    invivo_objective_seed_dir = invivo_input,
    invitro_objective_seed_dir = invitro_input,
    reductions = argv$reductions %||% "umap,pca,tsne",
    preprocess_modes = argv$pooled_preprocess_modes %||% "zscore,context_prior_unit,common_prior_unit",
    run_pca_umap = "TRUE",
    run_full = argv$run_pooled_full %||% argv$run_full %||% "TRUE",
    run_sampled = argv$run_pooled_sampled %||% argv$run_sampled %||% "TRUE",
    drop_parameter_table_initial = "TRUE",
    drop_invitro_parameter_table_initial = "TRUE",
    umap_seed = argv$umap_seed %||% "123",
    cluster_seed = argv$cluster_seed %||% "123",
    n_threads = argv$n_threads
  ))
  run_child_script("Generate pooled in vivo/in vitro reductions", "clustering_analysis.R", args, dry_run = dry_run)
}

run_cluster_report <- function(label, argv, result_root, dry_run, reductions, output_html) {
  args <- append_args(character(), list(
    result_root = result_root,
    reductions = reductions,
    output_html = output_html
  ))
  run_child_script(label, "clustering_report.R", args, dry_run = dry_run)
}

run_umap_report <- function(argv, result_root, dry_run) {
  if (!as_bool(argv$run_report %||% "TRUE", TRUE)) {
    message("Skipping clustering reports because --run_report=FALSE")
    return(invisible(NULL))
  }
  if (as_bool(argv$run_umap_cluster_report, TRUE)) {
    run_cluster_report(
      "Render full clustering HTML report",
      argv,
      result_root,
      dry_run,
      reductions = argv$cluster_report_reductions %||% "pca,umap,tsne",
      output_html = argv$umap_report_html %||% file.path(result_root, "parameter_landscape_clustering_umap_cluster_report.html")
    )
  } else {
    message("Skipping full clustering report because --run_umap_cluster_report=FALSE")
  }
  if (as_bool(argv$run_pca_cluster_report, TRUE)) {
    run_cluster_report(
      "Render PCA-only clustering HTML report",
      argv,
      result_root,
      dry_run,
      reductions = "pca",
      output_html = argv$pca_report_html %||% file.path(result_root, "parameter_landscape_clustering_pca_cluster_report.html")
    )
  } else {
    message("Skipping PCA-only clustering report because --run_pca_cluster_report=FALSE")
  }
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  if (as_bool(argv$help, FALSE) || as_bool(argv$h, FALSE)) {
    usage()
    return(invisible(NULL))
  }
  parts <- normalize_run_parts(argv$run_parts %||% argv$parts %||% "all")
  result_root <- normalizePath(
    path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()),
    mustWork = FALSE
  )
  invivo_input <- normalizePath(
    path.expand(argv$invivo_input %||% argv$input_dir %||% default_dataset_input_dir("invivo")),
    mustWork = FALSE
  )
  invitro_input <- normalizePath(
    path.expand(argv$invitro_input %||% default_dataset_input_dir("invitro")),
    mustWork = FALSE
  )
  dry_run <- as_bool(argv$dry_run, FALSE)

  message("Parameter landscape clustering runner")
  message("Result root: ", result_root)
  message("Run parts: ", paste(parts, collapse = ", "))
  message("Dry run: ", dry_run)

  for (part in parts) {
    switch(
      part,
      invivo_tables = run_invivo_tables(argv, result_root, invivo_input, dry_run),
      invitro_tables = run_invitro_tables(argv, result_root, invitro_input, dry_run),
      invivo_reductions = run_invivo_reductions(argv, result_root, invivo_input, dry_run),
      invitro_reductions = run_invitro_reductions(argv, result_root, invitro_input, dry_run),
      pooled_reductions = run_pooled_reductions(argv, result_root, invivo_input, invitro_input, dry_run),
      clustering_report = {
        run_umap_report(argv, result_root, dry_run)
      },
      stop("Unhandled run part: ", part)
    )
  }
  message("Parameter landscape clustering runner complete.")
}

main()
