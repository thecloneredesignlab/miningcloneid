#!/usr/bin/env Rscript

# Sequential orchestration for the joint coupling workflow:
# analysis tables -> visualization -> report.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

run_stage <- function(label, script, args, dry_run = FALSE) {
  command <- c(normalizePath(script, mustWork = TRUE), args)
  message("[", label, "] Rscript ", paste(shQuote(command), collapse = " "))
  if (dry_run) return(invisible(0L))
  status <- system2("Rscript", command)
  if (!identical(as.integer(status), 0L)) stop(label, " failed with exit status ", status, call. = FALSE)
  invisible(status)
}

run_joint_coupling_pipeline <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  args <- o2jca_parse_args(raw_args)
  result_root <- args$result_root %||% stop("--result_root is required", call. = FALSE)
  output_root <- args$output_root %||% o2jca_default_output_root(result_root, WORKFLOW_ROOT)
  o2jca_assert_separate_output_root(result_root, output_root)
  ratio_dir <- file.path(output_root, "tables", "soft_coupling")
  ploidy_dir <- file.path(output_root, "tables", "ploidy_coupling")
  overview_fig_dir <- file.path(output_root, "figures", "overview")
  ratio_fig_dir <- file.path(output_root, "figures", "soft_coupling")
  ploidy_fig_dir <- file.path(output_root, "figures", "ploidy_categories")
  association_fig_dir <- file.path(output_root, "figures", "category_association")
  appendix_fig_dir <- file.path(output_root, "figures", "appendix")
  report_dir <- file.path(output_root, "report")
  manifest_dir <- file.path(output_root, "manifests")
  log_dir <- file.path(output_root, "logs")
  dry_run <- o2jca_as_bool(args$dry_run, FALSE)
  run_vis <- o2jca_as_bool(args$run_vis, TRUE)
  run_report <- o2jca_as_bool(args$run_report, TRUE)
  common <- c(
    paste0("--result_root=", result_root), paste0("--class_threshold=", o2jca_as_num(args$class_threshold, 1.1)),
    paste0("--pair_pattern=", args$pair_pattern %||% "^fit_joint_.*_vt_seed[0-9]+$")
  )
  if (is.finite(o2jca_as_int(args$max_pairs, NA_integer_))) common <- c(common, paste0("--max_pairs=", o2jca_as_int(args$max_pairs)))
  if (is.finite(o2jca_as_int(args$max_seeds, NA_integer_))) common <- c(common, paste0("--max_seeds=", o2jca_as_int(args$max_seeds)))
  ratio_analysis <- file.path(WORKFLOW_ROOT, "analysis", "joint_soft_coupling_stability")
  ploidy_analysis <- file.path(WORKFLOW_ROOT, "analysis", "joint_ploidy_coupling_association")
  ratio_vis <- file.path(WORKFLOW_ROOT, "vis", "joint_soft_coupling_stability")
  ploidy_vis <- file.path(WORKFLOW_ROOT, "vis", "joint_ploidy_coupling_association")

  if (!dry_run) {
    dirs <- c(ratio_dir, ploidy_dir, overview_fig_dir, ratio_fig_dir, ploidy_fig_dir, association_fig_dir, appendix_fig_dir, report_dir, manifest_dir, log_dir)
    invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
    config <- data.frame(
      key = c("result_root", "output_root", "run_name", "class_threshold", "created_at", "input_is_read_only"),
      value = c(
        normalizePath(result_root, mustWork = TRUE), normalizePath(output_root, mustWork = TRUE),
        basename(sub("/+$", "", result_root)), o2jca_as_num(args$class_threshold, 1.1),
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "TRUE"
      ),
      stringsAsFactors = FALSE
    )
    o2jca_write_tsv(config, file.path(manifest_dir, "pipeline_config.tsv"))
  }

  run_stage("ratio master", file.path(ratio_analysis, "build_joint_soft_coupling_master_table.R"), c(common, paste0("--out_dir=", ratio_dir)), dry_run)
  run_stage("within-pair ratio", file.path(ratio_analysis, "analyze_within_pair_soft_coupling.R"), c(paste0("--out_dir=", ratio_dir), paste0("--thresholds=", args$thresholds %||% "1.05,1.075,1.1,1.15,1.2")), dry_run)
  run_stage("between-pair ratio", file.path(ratio_analysis, "analyze_between_pair_soft_coupling.R"), c(paste0("--out_dir=", ratio_dir), paste0("--n_boot=", o2jca_as_int(args$n_boot, 2000L))), dry_run)
  run_stage("process summary", file.path(ratio_analysis, "summarize_soft_coupling_processes.R"), c(paste0("--out_dir=", ratio_dir)), dry_run)

  run_stage("ploidy categories", file.path(ploidy_analysis, "classify_joint_ploidy_trajectories.R"), c(common, paste0("--out_dir=", ploidy_dir), paste0("--trajectory_step=", o2jca_as_int(args$trajectory_step, 10L))), dry_run)
  run_stage("within-pair ploidy coupling", file.path(ploidy_analysis, "analyze_within_pair_ploidy_coupling.R"), c(paste0("--out_dir=", ploidy_dir), paste0("--ratio_analysis_dir=", ratio_dir)), dry_run)
  run_stage("between-pair ploidy coupling", file.path(ploidy_analysis, "analyze_between_pair_ploidy_coupling.R"), c(paste0("--out_dir=", ploidy_dir), paste0("--n_boot=", o2jca_as_int(args$n_boot, 2000L))), dry_run)
  run_stage("Cat x Class association", file.path(ploidy_analysis, "analyze_ploidy_ratio_class_association.R"), c(paste0("--out_dir=", ploidy_dir), paste0("--permutations=", o2jca_as_int(args$permutations, 999L))), dry_run)

  if (run_vis) {
    run_stage("joint overview figures", file.path(ratio_vis, "plot_joint_coupling_overview.R"), c(paste0("--ratio_analysis_dir=", ratio_dir), paste0("--ploidy_analysis_dir=", ploidy_dir), paste0("--out_dir=", overview_fig_dir)), dry_run)
    run_stage("within-pair ratio figures", file.path(ratio_vis, "plot_within_pair_soft_coupling.R"), c(paste0("--analysis_dir=", ratio_dir), paste0("--out_dir=", ratio_fig_dir)), dry_run)
    run_stage("between-pair ratio figures", file.path(ratio_vis, "plot_between_pair_soft_coupling.R"), c(paste0("--analysis_dir=", ratio_dir), paste0("--out_dir=", ratio_fig_dir), paste0("--appendix_dir=", appendix_fig_dir)), dry_run)
    run_stage("process figures", file.path(ratio_vis, "plot_soft_coupling_processes.R"), c(paste0("--analysis_dir=", ratio_dir), paste0("--out_dir=", ratio_fig_dir)), dry_run)
    run_stage("ploidy category figures", file.path(ploidy_vis, "plot_joint_ploidy_categories.R"), c(paste0("--analysis_dir=", ploidy_dir), paste0("--out_dir=", ploidy_fig_dir), paste0("--appendix_dir=", appendix_fig_dir)), dry_run)
    run_stage("ploidy coupling figures", file.path(ploidy_vis, "plot_ploidy_coupling_associations.R"), c(paste0("--analysis_dir=", ploidy_dir), paste0("--out_dir=", association_fig_dir)), dry_run)
  }
  if (run_report) {
    run_stage(
      "joint coupling report",
      file.path(WORKFLOW_ROOT, "report", "joint_coupling", "render_joint_coupling_report.R"),
      c(
        paste0("--ratio_analysis_dir=", ratio_dir), paste0("--ploidy_analysis_dir=", ploidy_dir),
        paste0("--figure_root=", file.path(output_root, "figures")),
        paste0("--out_dir=", report_dir)
      ),
      dry_run
    )
  }
  if (!dry_run) {
    manifests <- list.files(output_root, pattern = "manifest[.]tsv$", recursive = TRUE, full.names = TRUE)
    manifests <- manifests[!o2jca_path_is_within(manifests, manifest_dir)]
    if (length(manifests)) {
      manifest_index <- data.frame(
        manifest = basename(manifests),
        relative_path = substring(normalizePath(manifests, mustWork = TRUE), nchar(normalizePath(output_root, mustWork = TRUE)) + 2L),
        stringsAsFactors = FALSE
      )
      o2jca_write_tsv(manifest_index, file.path(manifest_dir, "manifest_index.tsv"))
    }
  }
  invisible(list(output_root = output_root, ratio_dir = ratio_dir, ploidy_dir = ploidy_dir, report_dir = report_dir))
}

if (identical(environment(), globalenv())) run_joint_coupling_pipeline()
