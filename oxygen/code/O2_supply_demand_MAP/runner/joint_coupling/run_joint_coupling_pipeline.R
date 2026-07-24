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

reuse_ploidy_classification <- function(source_dir, out_dir, dry_run = FALSE) {
  manifest_path <- file.path(source_dir, "ploidy_classification_manifest.tsv")
  if (!file.exists(manifest_path)) stop("Missing source ploidy classification manifest: ", manifest_path, call. = FALSE)
  manifest <- o2jca_read_tsv(manifest_path)
  o2jca_assert_columns(manifest, "file", "source ploidy classification manifest")
  source_files <- file.path(source_dir, manifest$file)
  missing <- source_files[!file.exists(source_files)]
  if (length(missing)) stop("Missing source ploidy classification files: ", paste(missing, collapse = ", "), call. = FALSE)
  message("[ploidy categories] reuse ", length(source_files), " validated materialized tables from ", source_dir)
  if (dry_run) return(invisible(character()))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  target_files <- file.path(out_dir, basename(source_files))
  copied <- file.copy(source_files, target_files, overwrite = TRUE, copy.mode = TRUE)
  if (!all(copied)) stop("Failed to reuse ploidy classification files", call. = FALSE)
  o2jca_write_manifest(target_files, file.path(out_dir, "ploidy_classification_manifest.tsv"))
  invisible(target_files)
}

resolve_fixed_o2_source_root <- function(args, result_root) {
  explicit <- args$fixed_o2_source_root %||% ""
  if (nzchar(explicit)) return(normalizePath(path.expand(explicit), mustWork = TRUE))
  search_root <- file.path(
    o2jca_repo_root(WORKFLOW_ROOT), "oxygen", "results", "analysis",
    "warm_up_joint_fitting_results_extra"
  )
  if (!dir.exists(search_root)) {
    stop("--fixed_o2_source_root is required because the default search root is missing: ", search_root, call. = FALSE)
  }
  fit_name <- basename(sub("/+$", "", result_root))
  candidates <- list.dirs(search_root, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[startsWith(basename(candidates), fit_name)]
  curve_rel <- file.path(
    "curve_classification", "dense-grid_monotonicity_regression_classification",
    "tables", "fixed_o2_ploidy_monotonicity_regression_curves.tsv"
  )
  class_rel <- file.path(
    "curve_classification", "dense-grid_monotonicity_regression_classification",
    "tables", "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"
  )
  manifest_rel <- file.path("tables", "joint_best_curve_synthetic_run_manifest.tsv")
  candidates <- candidates[
    file.exists(file.path(candidates, curve_rel)) &
      file.exists(file.path(candidates, class_rel)) &
      file.exists(file.path(candidates, manifest_rel))
  ]
  if (length(candidates) != 1L) {
    stop(
      "Could not resolve exactly one complete fixed-O2 source for ", fit_name,
      "; supply --fixed_o2_source_root. Candidates: ", paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

fixed_o2_source_paths <- function(source_root) {
  table_dir <- file.path(
    source_root, "curve_classification", "dense-grid_monotonicity_regression_classification", "tables"
  )
  paths <- list(
    curve_table = file.path(table_dir, "fixed_o2_ploidy_monotonicity_regression_curves.tsv"),
    class_table = file.path(table_dir, "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"),
    seed_manifest = file.path(source_root, "tables", "joint_best_curve_synthetic_run_manifest.tsv")
  )
  missing <- unlist(paths, use.names = FALSE)[!file.exists(unlist(paths, use.names = FALSE))]
  if (length(missing)) stop("Missing fixed-O2 source files: ", paste(missing, collapse = ", "), call. = FALSE)
  lapply(paths, normalizePath, mustWork = TRUE)
}

run_joint_coupling_pipeline <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  args <- o2jca_parse_args(raw_args)
  source_analysis_root <- args$source_analysis_root %||% NULL
  source_mode <- !is.null(source_analysis_root) && nzchar(source_analysis_root)
  if (source_mode) {
    source_analysis_root <- normalizePath(source_analysis_root, mustWork = TRUE)
    source_ratio_dir <- file.path(source_analysis_root, "tables", "soft_coupling")
    source_ploidy_dir <- file.path(source_analysis_root, "tables", "ploidy_coupling")
    source_config <- o2jca_read_tsv(file.path(source_ratio_dir, "analysis_config.tsv"))
    source_result_root <- source_config$value[source_config$key == "result_root"]
    result_root <- args$result_root %||% if (length(source_result_root)) source_result_root[[1L]] else source_analysis_root
  } else {
    result_root <- args$result_root %||% stop("--result_root is required unless --source_analysis_root is supplied", call. = FALSE)
  }
  output_root <- args$output_root %||% o2jca_default_output_root(result_root, WORKFLOW_ROOT)
  o2jca_assert_separate_output_root(result_root, output_root)
  if (source_mode && o2jca_path_is_within(output_root, source_analysis_root)) {
    stop("output_root must be outside source_analysis_root", call. = FALSE)
  }
  ratio_dir <- file.path(output_root, "tables", "soft_coupling")
  ploidy_dir <- file.path(output_root, "tables", "ploidy_coupling")
  fixed_o2_dir <- file.path(output_root, "tables", "fixed_o2_ploidy_classification")
  overview_fig_dir <- file.path(output_root, "figures", "overview")
  ratio_fig_dir <- file.path(output_root, "figures", "soft_coupling")
  ploidy_fig_dir <- file.path(output_root, "figures", "ploidy_categories")
  fixed_o2_fig_dir <- file.path(output_root, "figures", "fixed_o2_ploidy_classification")
  association_fig_dir <- file.path(output_root, "figures", "category_association")
  appendix_fig_dir <- file.path(output_root, "figures", "appendix")
  report_dir <- file.path(output_root, "report")
  manifest_dir <- file.path(output_root, "manifests")
  log_dir <- file.path(output_root, "logs")
  comparison_dir <- file.path(output_root, "tables", "classification_comparison")
  dry_run <- o2jca_as_bool(args$dry_run, FALSE)
  run_vis <- o2jca_as_bool(args$run_vis, TRUE)
  run_report <- o2jca_as_bool(args$run_report, TRUE)
  class_spec <- o2jca_classification_spec(
    threshold = o2jca_as_num(args$class_threshold, 1.1),
    lower_bound = o2jca_as_num(args$class_lower_bound, NA_real_),
    upper_bound = o2jca_as_num(args$class_upper_bound, NA_real_),
    boundary_rule = args$class_boundary_rule %||% "classb_inclusive"
  )
  fixed_o2_source_root <- resolve_fixed_o2_source_root(args, result_root)
  fixed_o2_sources <- fixed_o2_source_paths(fixed_o2_source_root)
  common <- c(
    paste0("--result_root=", result_root), paste0("--class_threshold=", class_spec$threshold),
    paste0("--class_lower_bound=", class_spec$lower_bound),
    paste0("--class_upper_bound=", class_spec$upper_bound),
    paste0("--class_boundary_rule=", class_spec$boundary_rule),
    paste0("--pair_pattern=", args$pair_pattern %||% "^fit_joint_.*_vt_seed[0-9]+$")
  )
  if (is.finite(o2jca_as_int(args$max_pairs, NA_integer_))) common <- c(common, paste0("--max_pairs=", o2jca_as_int(args$max_pairs)))
  if (is.finite(o2jca_as_int(args$max_seeds, NA_integer_))) common <- c(common, paste0("--max_seeds=", o2jca_as_int(args$max_seeds)))
  ratio_analysis <- file.path(WORKFLOW_ROOT, "analysis", "joint_soft_coupling_stability")
  ploidy_analysis <- file.path(WORKFLOW_ROOT, "analysis", "joint_ploidy_coupling_association")
  fixed_o2_analysis <- file.path(WORKFLOW_ROOT, "analysis", "joint_fixed_o2_ploidy_classification")
  ratio_vis <- file.path(WORKFLOW_ROOT, "vis", "joint_soft_coupling_stability")
  ploidy_vis <- file.path(WORKFLOW_ROOT, "vis", "joint_ploidy_coupling_association")
  fixed_o2_vis <- file.path(WORKFLOW_ROOT, "vis", "joint_fixed_o2_ploidy_classification")

  if (!dry_run) {
    dirs <- c(ratio_dir, ploidy_dir, fixed_o2_dir, comparison_dir, overview_fig_dir, ratio_fig_dir, ploidy_fig_dir, fixed_o2_fig_dir, association_fig_dir, appendix_fig_dir, report_dir, manifest_dir, log_dir)
    invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
    config <- data.frame(
      key = c(
        "result_root", "source_analysis_root", "output_root", "run_name",
        "class_threshold", "class_lower_bound", "class_upper_bound",
        "class_boundary_rule", "class_scheme", "class_rule", "created_at",
        "input_is_read_only", "ploidy_classification_reused", "fixed_o2_source_root"
      ),
      value = c(
        if (source_mode) result_root else normalizePath(result_root, mustWork = TRUE),
        if (source_mode) source_analysis_root else NA_character_,
        normalizePath(output_root, mustWork = TRUE), basename(sub("/+$", "", output_root)),
        class_spec$threshold, class_spec$lower_bound, class_spec$upper_bound,
        class_spec$boundary_rule, class_spec$scheme, class_spec$class_rule,
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "TRUE", source_mode,
        fixed_o2_source_root
      ),
      stringsAsFactors = FALSE
    )
    o2jca_write_tsv(config, file.path(manifest_dir, "pipeline_config.tsv"))
  }

  if (source_mode) {
    run_stage(
      "ratio master reclassification",
      file.path(ratio_analysis, "reclassify_joint_soft_coupling_master_table.R"),
      c(
        paste0("--source_analysis_dir=", source_ratio_dir),
        paste0("--out_dir=", ratio_dir),
        paste0("--analysis_label=", basename(sub("/+$", "", output_root))),
        paste0("--class_threshold=", class_spec$threshold),
        paste0("--class_lower_bound=", class_spec$lower_bound),
        paste0("--class_upper_bound=", class_spec$upper_bound),
        paste0("--class_boundary_rule=", class_spec$boundary_rule)
      ),
      dry_run
    )
  } else {
    run_stage("ratio master", file.path(ratio_analysis, "build_joint_soft_coupling_master_table.R"), c(common, paste0("--out_dir=", ratio_dir)), dry_run)
  }
  run_stage("within-pair ratio", file.path(ratio_analysis, "analyze_within_pair_soft_coupling.R"), c(paste0("--out_dir=", ratio_dir), paste0("--thresholds=", args$thresholds %||% "1.05,1.075,1.1,1.15,1.2")), dry_run)
  run_stage("between-pair ratio", file.path(ratio_analysis, "analyze_between_pair_soft_coupling.R"), c(paste0("--out_dir=", ratio_dir), paste0("--n_boot=", o2jca_as_int(args$n_boot, 2000L))), dry_run)
  run_stage("process summary", file.path(ratio_analysis, "summarize_soft_coupling_processes.R"), c(paste0("--out_dir=", ratio_dir)), dry_run)
  if (source_mode) {
    run_stage(
      "classification version comparison",
      file.path(ratio_analysis, "compare_joint_coupling_classification_versions.R"),
      c(
        paste0("--source_analysis_dir=", source_ratio_dir),
        paste0("--current_analysis_dir=", ratio_dir),
        paste0("--out_dir=", comparison_dir)
      ),
      dry_run
    )
  }

  if (source_mode) {
    reuse_ploidy_classification(source_ploidy_dir, ploidy_dir, dry_run)
  } else {
    run_stage("ploidy categories", file.path(ploidy_analysis, "classify_joint_ploidy_trajectories.R"), c(common, paste0("--out_dir=", ploidy_dir), paste0("--trajectory_step=", o2jca_as_int(args$trajectory_step, 10L))), dry_run)
  }
  run_stage("within-pair ploidy coupling", file.path(ploidy_analysis, "analyze_within_pair_ploidy_coupling.R"), c(paste0("--out_dir=", ploidy_dir), paste0("--ratio_analysis_dir=", ratio_dir)), dry_run)
  run_stage("between-pair ploidy coupling", file.path(ploidy_analysis, "analyze_between_pair_ploidy_coupling.R"), c(paste0("--out_dir=", ploidy_dir), paste0("--n_boot=", o2jca_as_int(args$n_boot, 2000L))), dry_run)
  run_stage("Cat x Class association", file.path(ploidy_analysis, "analyze_ploidy_ratio_class_association.R"), c(paste0("--out_dir=", ploidy_dir), paste0("--permutations=", o2jca_as_int(args$permutations, 999L))), dry_run)
  run_stage(
    "fixed-O2 steady-state ploidy classification",
    file.path(fixed_o2_analysis, "analyze_joint_fixed_o2_ploidy_classes.R"),
    c(
      paste0("--curve_table=", fixed_o2_sources$curve_table),
      paste0("--class_table=", fixed_o2_sources$class_table),
      paste0("--seed_manifest=", fixed_o2_sources$seed_manifest),
      paste0("--ploidy_category_table=", file.path(ploidy_dir, "ploidy_pair_category_assignment.tsv")),
      paste0("--out_dir=", fixed_o2_dir)
    ),
    dry_run
  )

  if (run_vis) {
    run_stage("joint overview figures", file.path(ratio_vis, "plot_joint_coupling_overview.R"), c(paste0("--ratio_analysis_dir=", ratio_dir), paste0("--ploidy_analysis_dir=", ploidy_dir), paste0("--out_dir=", overview_fig_dir)), dry_run)
    run_stage("within-pair ratio figures", file.path(ratio_vis, "plot_within_pair_soft_coupling.R"), c(paste0("--analysis_dir=", ratio_dir), paste0("--out_dir=", ratio_fig_dir)), dry_run)
    run_stage("between-pair ratio figures", file.path(ratio_vis, "plot_between_pair_soft_coupling.R"), c(paste0("--analysis_dir=", ratio_dir), paste0("--out_dir=", ratio_fig_dir), paste0("--appendix_dir=", appendix_fig_dir)), dry_run)
    run_stage("process figures", file.path(ratio_vis, "plot_soft_coupling_processes.R"), c(paste0("--analysis_dir=", ratio_dir), paste0("--out_dir=", ratio_fig_dir)), dry_run)
    run_stage("ploidy category figures", file.path(ploidy_vis, "plot_joint_ploidy_categories.R"), c(paste0("--analysis_dir=", ploidy_dir), paste0("--out_dir=", ploidy_fig_dir), paste0("--appendix_dir=", appendix_fig_dir)), dry_run)
    run_stage("ploidy coupling figures", file.path(ploidy_vis, "plot_ploidy_coupling_associations.R"), c(paste0("--analysis_dir=", ploidy_dir), paste0("--out_dir=", association_fig_dir)), dry_run)
    run_stage("fixed-O2 steady-state ploidy figures", file.path(fixed_o2_vis, "plot_joint_fixed_o2_ploidy_curves.R"), c(paste0("--analysis_dir=", fixed_o2_dir), paste0("--out_dir=", fixed_o2_fig_dir)), dry_run)
  }
  if (run_report) {
    run_stage(
      "joint coupling report",
      file.path(WORKFLOW_ROOT, "report", "joint_coupling", "render_joint_coupling_report.R"),
      c(
        paste0("--ratio_analysis_dir=", ratio_dir), paste0("--ploidy_analysis_dir=", ploidy_dir),
        paste0("--fixed_o2_analysis_dir=", fixed_o2_dir),
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
  invisible(list(output_root = output_root, ratio_dir = ratio_dir, ploidy_dir = ploidy_dir, fixed_o2_dir = fixed_o2_dir, report_dir = report_dir))
}

if (identical(environment(), globalenv())) run_joint_coupling_pipeline()
