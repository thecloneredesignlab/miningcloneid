.process_runner_path <- local({
  frames <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(frame) {
        path <- frame$ofile
        if (is.null(path)) "" else normalizePath(path, mustWork = FALSE)
      },
      character(1L)
    )
  )
  if (length(frames)) {
    frames[[length(frames)]]
  } else {
    file.path(getwd(), "process_runner.R")
  }
})
source(file.path(dirname(.process_runner_path), "workspace_paths.R"))
rm(.process_runner_path)

resolve_utility_file <- function(name) {
  path <- file.path(CODE_ROOT, "util", name)
  if (!file.exists(path)) stop("Missing utility implementation: ", path)
  normalizePath(path, mustWork = TRUE)
}

compose_three_row_raster <- function(
    data_dir,
    panel_dir,
    output_dir,
    panels,
    output_basename,
    validation_basename
) {
  expected_labels <- LETTERS[1:7]
  if (!identical(names(panels), expected_labels)) {
    stop("Three-row compositor requires named panels A through G.")
  }
  run_process(
    "python3",
    resolve_utility_file("graphics/raster_compositor.py"),
    env = c(
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("COMPOSITOR_DATA_DIR=", data_dir),
      paste0("COMPOSITOR_PANEL_DIR=", panel_dir),
      paste0("COMPOSITOR_OUTPUT_DIR=", output_dir),
      paste0("COMPOSITOR_PANEL_", names(panels), "=", panels),
      paste0("COMPOSITOR_OUTPUT_BASENAME=", output_basename),
      paste0("COMPOSITOR_VALIDATION_BASENAME=", validation_basename)
    )
  )
}

run_fixed_o2_auc <- function(args, env = character()) {
  run_process(
    "Rscript",
    args = c(resolve_utility_file("analysis/fixed_o2_auc.R"), args),
    env = env
  )
}

run_figure4_continuous_association <- function(data_dir) {
  run_process(
    "Rscript",
    args = c(
      resolve_utility_file(
        "analysis/figure4_continuous_ploidy_association.R"
      ),
      paste0("--data-dir=", data_dir)
    )
  )
}

run_invivo_tsne_landscape <- function(args, env = character()) {
  run_process(
    "Rscript",
    args = c(
      resolve_utility_file("analysis/invivo_tsne_landscape.R"),
      args
    ),
    env = env
  )
}

render_parameter_landscape <- function(
    data_dir,
    panel_dir,
    output_dir
) {
  run_process(
    "Rscript",
    resolve_utility_file("analysis/parameter_landscape.R"),
    env = c(
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("ANALYSIS_DATA_DIR=", data_dir),
      paste0("PLOT_OUTPUT_DIR=", panel_dir),
      paste0("DELIVERABLE_OUTPUT_DIR=", output_dir)
    )
  )
}

compose_stacked_panels <- function(
    data_dir,
    panel_dir,
    output_dir,
    top_panel,
    bottom_panel,
    supplementary_source,
    supplementary_destination,
    output_basename,
    validation_basename
) {
  run_process(
    "python3",
    resolve_utility_file("graphics/stacked_panel_compositor.py"),
    env = c(
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("COMPOSITOR_DATA_DIR=", data_dir),
      paste0("COMPOSITOR_PANEL_DIR=", panel_dir),
      paste0("COMPOSITOR_OUTPUT_DIR=", output_dir),
      paste0("COMPOSITOR_TOP_PANEL=", top_panel),
      paste0("COMPOSITOR_BOTTOM_PANEL=", bottom_panel),
      paste0("COMPOSITOR_SUPPLEMENTARY_SOURCE=", supplementary_source),
      paste0(
        "COMPOSITOR_SUPPLEMENTARY_DESTINATION=",
        supplementary_destination
      ),
      paste0("COMPOSITOR_OUTPUT_BASENAME=", output_basename),
      paste0("COMPOSITOR_VALIDATION_BASENAME=", validation_basename)
    )
  )
}

compose_top_with_left_and_stacked_right_panels <- function(
    data_dir,
    panel_dir,
    output_dir,
    top_panel,
    bottom_left_panel,
    bottom_right_top_panel,
    bottom_right_bottom_panel,
    supplementary_source,
    supplementary_destination,
    output_basename,
    validation_basename
) {
  run_process(
    "python3",
    resolve_utility_file(
      "graphics/top_with_two_column_bottom_compositor.py"
    ),
    env = c(
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("COMPOSITOR_DATA_DIR=", data_dir),
      paste0("COMPOSITOR_PANEL_DIR=", panel_dir),
      paste0("COMPOSITOR_OUTPUT_DIR=", output_dir),
      paste0("COMPOSITOR_TOP_PANEL=", top_panel),
      paste0("COMPOSITOR_BOTTOM_LEFT_PANEL=", bottom_left_panel),
      paste0(
        "COMPOSITOR_BOTTOM_RIGHT_TOP_PANEL=",
        bottom_right_top_panel
      ),
      paste0(
        "COMPOSITOR_BOTTOM_RIGHT_BOTTOM_PANEL=",
        bottom_right_bottom_panel
      ),
      paste0("COMPOSITOR_SUPPLEMENTARY_SOURCE=", supplementary_source),
      paste0(
        "COMPOSITOR_SUPPLEMENTARY_DESTINATION=",
        supplementary_destination
      ),
      paste0("COMPOSITOR_OUTPUT_BASENAME=", output_basename),
      paste0("COMPOSITOR_VALIDATION_BASENAME=", validation_basename)
    )
  )
}

compose_three_stacked_panels <- function(
    data_dir,
    panel_dir,
    output_dir,
    panel_a,
    panel_b,
    panel_c,
    supplementary_source,
    supplementary_destination,
    output_basename,
    validation_basename
) {
  run_process(
    "python3",
    resolve_utility_file("graphics/three_stacked_panel_compositor.py"),
    env = c(
      paste0("COMPOSITOR_DATA_DIR=", data_dir),
      paste0("COMPOSITOR_PANEL_DIR=", panel_dir),
      paste0("COMPOSITOR_OUTPUT_DIR=", output_dir),
      paste0("COMPOSITOR_PANEL_A=", panel_a),
      paste0("COMPOSITOR_PANEL_B=", panel_b),
      paste0("COMPOSITOR_PANEL_C=", panel_c),
      paste0("COMPOSITOR_SUPPLEMENTARY_SOURCE=", supplementary_source),
      paste0(
        "COMPOSITOR_SUPPLEMENTARY_DESTINATION=",
        supplementary_destination
      ),
      paste0("COMPOSITOR_OUTPUT_BASENAME=", output_basename),
      paste0("COMPOSITOR_VALIDATION_BASENAME=", validation_basename)
    )
  )
}

run_oxygen_response_pipeline <- function(args, env = character()) {
  run_process(
    "Rscript",
    args = c(resolve_utility_file("oxygen/response_pipeline.R"), args),
    env = env
  )
}

stage_output_file <- function(
    source,
    destination_name = basename(source)
) {
  if (!file.exists(source)) stop("Missing generated figure: ", source)
  destination <- file.path(OUTPUT_ROOT, destination_name)
  ok <- file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE)
  if (!ok) stop("Failed to stage figure: ", source, " -> ", destination)
  normalizePath(destination, mustWork = TRUE)
}

require_files <- function(paths, label = "required file") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing ", label, "(s):\n", paste(missing, collapse = "\n"))
  }
  invisible(normalizePath(paths, mustWork = TRUE))
}

write_data_contract <- function(artifact_id, rows) {
  path <- file.path(
    DATA_ROOT, artifact_id, "data_contract.tsv"
  )
  rows$figure_id <- artifact_id
  rows <- rows[, c("figure_id", setdiff(names(rows), "figure_id")), drop = FALSE]
  write_intermediate_tsv(rows, path)
}
