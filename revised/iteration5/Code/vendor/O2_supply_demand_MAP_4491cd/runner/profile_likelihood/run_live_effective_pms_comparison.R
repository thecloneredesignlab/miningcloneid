#!/usr/bin/env Rscript

# Orchestrates the separated live-effective-p_ms analysis and visualization.
# Upstream per-seed simulations must already be materialized.
#
# Usage:
#   Rscript oxygen/code/O2_supply_demand_MAP/runner/profile_likelihood/run_live_effective_pms_comparison.R \
#     --run_dir_template=/path/to/run_sigma_{sigma} \
#     --sigma_caps=0.05,0.15,0.3,0.6 \
#     --out_dir=/path/to/comparison

.o2sd_profile_runner_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_profile_runner_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_profile_runner_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_bool <- o2sd_as_bool
resolve_path_value <- o2sd_resolve_path
write_tsv <- o2sd_write_tsv

add_cli_arg <- function(args, name, value) {
  if (is.null(value) || !length(value)) return(args)
  txt <- trimws(as.character(value[[1]]))
  if (!nzchar(txt)) return(args)
  c(args, paste0("--", name, "=", txt))
}

run_rscript <- function(script, args, label) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = shQuote(c(script, args))
  )
  if (!identical(status, 0L)) {
    stop(label, " failed with exit status ", status)
  }
  invisible(status)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  out_dir <- resolve_path_value(argv$out_dir, getwd())
  if (is.null(out_dir)) stop("--out_dir is required.")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_script <- file.path(
    WORKFLOW_ROOT,
    "analysis",
    "profile_likelihood",
    "compare_sigma_burden_live_effective_pms.R"
  )
  vis_script <- file.path(
    WORKFLOW_ROOT,
    "vis",
    "profile_likelihood",
    "plot_sigma_burden_live_effective_pms.R"
  )
  if (!file.exists(analysis_script)) stop("Missing analysis entrypoint: ", analysis_script)
  if (!file.exists(vis_script)) stop("Missing visualization entrypoint: ", vis_script)

  analysis_args <- paste0("--out_dir=", out_dir)
  for (name in c(
    "run_dir_template",
    "sigma_caps",
    "live_effective_subdir",
    "max_seeds",
    "n_cores"
  )) {
    analysis_args <- add_cli_arg(analysis_args, name, argv[[name]])
  }
  run_rscript(
    analysis_script,
    analysis_args,
    "Live-effective-p_ms comparison analysis"
  )

  skip_visualization <- as_bool(argv$skip_visualization, FALSE)
  figure_dir <- resolve_path_value(argv$figure_dir, getwd()) %||%
    file.path(out_dir, "figures")
  figure_dir <- normalizePath(figure_dir, mustWork = FALSE)
  if (!skip_visualization) {
    vis_args <- c(
      paste0("--analysis_dir=", out_dir),
      paste0("--out_dir=", figure_dir)
    )
    vis_args <- add_cli_arg(vis_args, "sigma_caps", argv$sigma_caps)
    run_rscript(
      vis_script,
      vis_args,
      "Live-effective-p_ms comparison visualization"
    )
  }

  manifest <- data.frame(
    key = c(
      "schema_version",
      "status",
      "analysis_dir",
      "visualization_status",
      "figure_dir",
      "generated_at"
    ),
    value = c(
      "o2sd-live-effective-pms-pipeline-v1",
      "complete",
      out_dir,
      if (skip_visualization) "skipped" else "complete",
      if (skip_visualization) "" else figure_dir,
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(out_dir, "live_effective_pms_pipeline_manifest.tsv")
  write_tsv(manifest, manifest_path)
  message("Wrote pipeline manifest: ", manifest_path)
}

if (sys.nframe() == 0) {
  main()
}
