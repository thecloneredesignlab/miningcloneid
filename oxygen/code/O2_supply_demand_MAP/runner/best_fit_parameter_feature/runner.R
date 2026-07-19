#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
UTIL_DIR <- file.path(WORKFLOW_ROOT, "analysis", "best_fit_parameter_feature", "util")
source(file.path(UTIL_DIR, "path_utils.R"))
source(file.path(UTIL_DIR, "cli_utils.R"))

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript runner.R --workflow=all [options]",
      "  Rscript runner.R --workflow=fixed_o2 [options]",
      "  Rscript runner.R --workflow=parameter_landscape [options]",
      "  Rscript runner.R --workflow=dense_grid_monotonicity [options]",
      "  Rscript runner.R --workflow=combine [options]",
      "  Rscript runner.R --workflow=combine_report [options]",
      "",
      "Workflow selection:",
      "  --workflow=NAME                 all, fixed_o2, parameter_landscape, dense_grid_monotonicity, combine, combine_report.",
      "  --parameter_parts=PARTS         clustering, mode_contribution, dominant_ploidy_contribution, or all.",
      "  --dense_parts=PARTS             monotonicity, regression, initial_ploidy, report, compare, or all.",
      "  --combine_parts=PARTS           pooled_curve_class, average_slope, report, fixo2_eigen_attractor, fixo2_eigen_report, or all.",
      "",
      "Shared path defaults:",
      "  fixed_o2 out_dir                oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo_500seed",
      "  parameter result_root           oxygen/results/analysis/best_fit_parameter_feature/02_parameter_landscape_clustering",
      "  dense-grid result_root          oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification",
      "  combine out_dir                 oxygen/results/analysis/best_fit_parameter_feature/04_combine_parameter_landscape/pooled_embedding_curve_class",
      "  combine_report output_html      oxygen/results/analysis/best_fit_parameter_feature/04_combine_parameter_landscape/pooled_embedding_curve_class/pooled_embedding_curve_class_report.html",
      "  fixo2_eigen_attractor out_dir   oxygen/results/analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/pooled_embedding_curve_class",
      "",
      "Other --key=value options are forwarded to the selected child workflow scripts.",
      sep = "\n"
    ),
    "\n"
  )
}

normalize_workflows <- function(x) {
  vals <- tolower(gsub("-", "_", bpf_split_csv(x, "all")))
  out <- character()
  for (val in vals) {
    if (val %in% c("all", "full")) {
      out <- c(out, "fixed_o2", "parameter_landscape", "dense_grid_monotonicity", "combine")
    } else if (val %in% c("fixed_o2", "fixo2", "fixed")) {
      out <- c(out, "fixed_o2")
    } else if (val %in% c("parameter_landscape", "parameter_landscape_clustering", "landscape")) {
      out <- c(out, "parameter_landscape")
    } else if (val %in% c("dense_grid_monotonicity", "dense_grid", "dense_grid_monotonicity_classification", "monotonicity")) {
      out <- c(out, "dense_grid_monotonicity")
    } else if (val %in% c("combine", "combined", "integration", "integrate")) {
      out <- c(out, "combine")
    } else if (val %in% c("combine_report", "report_combine", "integration_report", "combined_report")) {
      out <- c(out, "combine_report")
    } else if (nzchar(val)) {
      stop("Unknown workflow: ", val, call. = FALSE)
    }
  }
  unique(out)
}

normalize_parameter_parts <- function(x) {
  vals <- tolower(gsub("-", "_", bpf_split_csv(x, "clustering")))
  out <- character()
  for (val in vals) {
    if (val %in% c("all", "full")) {
      out <- c(out, "clustering", "mode_contribution", "dominant_ploidy_contribution")
    } else if (val %in% c("clustering", "reduction", "reductions", "tables")) {
      out <- c(out, "clustering")
    } else if (val %in% c("mode", "mode_contribution", "discrete")) {
      out <- c(out, "mode_contribution")
    } else if (val %in% c("dominant_ploidy", "dominant_ploidy_contribution", "continuous")) {
      out <- c(out, "dominant_ploidy_contribution")
    } else if (nzchar(val)) {
      stop("Unknown parameter part: ", val, call. = FALSE)
    }
  }
  unique(out)
}

normalize_dense_parts <- function(x) {
  vals <- tolower(gsub("-", "_", bpf_split_csv(x, "monotonicity")))
  out <- character()
  for (val in vals) {
    if (val %in% c("all", "full")) {
      out <- c(out, "monotonicity", "regression", "initial_ploidy", "report", "compare")
    } else if (val %in% c("monotonicity", "classification", "pointwise")) {
      out <- c(out, "monotonicity")
    } else if (val %in% c("regression", "smooth", "smoothed")) {
      out <- c(out, "regression")
    } else if (val %in% c("initial_ploidy", "initial", "trajectory")) {
      out <- c(out, "initial_ploidy")
    } else if (val %in% c("report", "html_report")) {
      out <- c(out, "report")
    } else if (val %in% c("compare", "compare_1000d_10000d")) {
      out <- c(out, "compare")
    } else if (nzchar(val)) {
      stop("Unknown dense-grid part: ", val, call. = FALSE)
    }
  }
  unique(out)
}

normalize_combine_parts <- function(x) {
  vals <- tolower(gsub("-", "_", bpf_split_csv(x, "pooled_curve_class")))
  out <- character()
  for (val in vals) {
    if (val %in% c("all", "full")) {
      out <- c(out, "average_slope", "pooled_curve_class", "report")
    } else if (val %in% c("pooled_curve_class", "curve_class", "pooled_embeddings", "embedding_curve_class")) {
      out <- c(out, "pooled_curve_class")
    } else if (val %in% c("average_slope", "slope", "regression_slope", "curve_average_slope")) {
      out <- c(out, "average_slope")
    } else if (val %in% c("report", "html_report", "combine_report", "pooled_curve_class_report")) {
      out <- c(out, "report")
    } else if (val %in% c("fixo2_eigen_attractor", "fixo2_eigen", "fixo2_eigen_curve_class", "eigen_attractor")) {
      out <- c(out, "fixo2_eigen_attractor")
    } else if (val %in% c("fixo2_eigen_report", "fixo2_eigen_attractor_report", "eigen_attractor_report")) {
      out <- c(out, "fixo2_eigen_report")
    } else if (nzchar(val)) {
      stop("Unknown combine part: ", val, call. = FALSE)
    }
  }
  unique(out)
}

drop_forward_keys <- function(args, keys) {
  if (!length(args)) return(args)
  keep <- vapply(args, function(arg) {
    if (!startsWith(arg, "--")) return(TRUE)
    key <- sub("^--", "", arg)
    key <- sub("=.*$", "", key)
    key <- gsub("-", "_", key, fixed = TRUE)
    !(key %in% keys)
  }, logical(1))
  args[keep]
}

has_arg <- function(args, key) {
  key <- gsub("-", "_", key, fixed = TRUE)
  any(vapply(args, function(arg) {
    startsWith(arg, "--") && identical(gsub("-", "_", sub("=.*$", "", sub("^--", "", arg)), fixed = TRUE), key)
  }, logical(1)))
}

append_default_arg <- function(args, key, value) {
  if (has_arg(args, key)) return(args)
  c(args, paste0("--", key, "=", value))
}

run_rscript <- function(label, script, args, dry_run = FALSE) {
  script <- normalizePath(script, mustWork = TRUE)
  cmd <- c(script, args)
  message("[", label, "] Rscript ", paste(shQuote(cmd), collapse = " "))
  if (isTRUE(dry_run)) return(invisible(0L))
  status <- system2(file.path(R.home("bin"), "Rscript"), cmd)
  if (!identical(as.integer(status), 0L)) {
    stop(label, " failed with exit status ", status, call. = FALSE)
  }
  invisible(status)
}

run_fixed_o2 <- function(args, repo_root, dry_run) {
  script <- file.path(bpf_fixed_o2_code_dir(repo_root), "FixO2_invivo.R")
  args <- append_default_arg(args, "out_dir", bpf_fixed_o2_result_dir(repo_root))
  run_rscript("fixed_o2", script, args, dry_run)
}

run_parameter_landscape <- function(args, argv, repo_root, dry_run) {
  parts <- normalize_parameter_parts(argv$parameter_parts)
  args <- append_default_arg(args, "result_root", bpf_parameter_landscape_result_dir(repo_root))
  canonical_runner <- file.path(WORKFLOW_ROOT, "runner", "parameter_landscape", "run_parameter_landscape.R")
  for (part in parts) {
    if (identical(part, "clustering")) {
      run_rscript("parameter_landscape_clustering", canonical_runner, append_default_arg(args, "run_parts", "all"), dry_run)
    } else if (identical(part, "mode_contribution")) {
      run_rscript("mode_parameter_contribution", canonical_runner, append_default_arg(args, "run_parts", "mode_contribution"), dry_run)
    } else if (identical(part, "dominant_ploidy_contribution")) {
      run_rscript("dominant_ploidy_parameter_contribution", canonical_runner, append_default_arg(args, "run_parts", "dominant_ploidy_contribution"), dry_run)
    }
  }
}

run_dense_grid <- function(args, argv, repo_root, dry_run) {
  parts <- normalize_dense_parts(argv$dense_parts)
  code_dir <- bpf_dense_grid_code_dir(repo_root)
  canonical_runner <- file.path(WORKFLOW_ROOT, "runner", "dense_grid_monotonicity", "run_dense_grid_monotonicity.R")
  result_root <- bpf_dense_grid_result_root(repo_root)
  monotonicity_dir <- file.path(result_root, "dense-grid_monotonicity_classification")
  initial_dir <- file.path(result_root, "dense-grid_initial-ploidy_trajectory")
  regression_dir <- file.path(result_root, "dense-grid_monotonicity_regression_classification")
  compare_dir <- file.path(result_root, "compare")

  for (part in parts) {
    if (identical(part, "monotonicity")) {
      part_args <- append_default_arg(args, "out_dir", monotonicity_dir)
      part_args <- append_default_arg(part_args, "part", "monotonicity")
      part_args <- append_default_arg(part_args, "mode", "run")
      run_rscript("dense_grid_monotonicity", canonical_runner, part_args, dry_run)
    } else if (identical(part, "regression")) {
      part_args <- append_default_arg(args, "out_dir", monotonicity_dir)
      part_args <- append_default_arg(part_args, "part", "monotonicity")
      part_args <- append_default_arg(part_args, "mode", "analyze")
      run_rscript("dense_grid_monotonicity_regression", canonical_runner, part_args, dry_run)
    } else if (identical(part, "initial_ploidy")) {
      part_args <- append_default_arg(args, "out_dir", initial_dir)
      part_args <- append_default_arg(part_args, "part", "initial_ploidy")
      part_args <- append_default_arg(part_args, "mode", "run")
      run_rscript("dense_grid_initial_ploidy", canonical_runner, part_args, dry_run)
    } else if (identical(part, "report")) {
      part_args <- append_default_arg(args, "out_dir", monotonicity_dir)
      part_args <- append_default_arg(part_args, "part", "monotonicity")
      part_args <- append_default_arg(part_args, "mode", "report")
      run_rscript("dense_grid_monotonicity_report", canonical_runner, part_args, dry_run)
    } else if (identical(part, "compare")) {
      part_args <- append_default_arg(args, "compare_root", compare_dir)
      run_rscript("dense_grid_compare_report", file.path(code_dir, "compare_1000D_10000D_report.R"), part_args, dry_run)
    }
  }
}

run_combine <- function(args, argv, repo_root, dry_run) {
  parts <- normalize_combine_parts(argv$combine_parts)
  code_dir <- bpf_combine_code_dir(repo_root)
  fixo2_code_dir <- bpf_combine_fixo2_eigen_attractor_code_dir(repo_root)
  canonical_combined_runner <- file.path(WORKFLOW_ROOT, "runner", "combined_parameter_landscape", "run_combined_parameter_landscape.R")
  canonical_fixo2_runner <- file.path(WORKFLOW_ROOT, "runner", "combined_fixo2_eigen", "run_fixo2_eigen_curve_class_pipeline.R")
  for (part in parts) {
    if (identical(part, "pooled_curve_class")) {
      part_args <- append_default_arg(
        args,
        "out_dir",
        file.path(bpf_combine_result_dir(repo_root), "pooled_embedding_curve_class")
      )
      part_args <- append_default_arg(part_args, "mode", "run")
      part_args <- append_default_arg(part_args, "run_report", "false")
      run_rscript(
        "combine_pooled_embedding_curve_class",
        canonical_combined_runner,
        part_args,
        dry_run
      )
    } else if (identical(part, "average_slope")) {
      part_args <- append_default_arg(
        args,
        "out_dir",
        file.path(bpf_combine_result_dir(repo_root), "pooled_embedding_curve_class")
      )
      part_args <- append_default_arg(part_args, "mode", "analyze")
      part_args <- append_default_arg(part_args, "slope_only", "true")
      run_rscript(
        "combine_regression_curve_average_slope",
        canonical_combined_runner,
        part_args,
        dry_run
      )
    } else if (identical(part, "report")) {
      part_args <- append_default_arg(args, "out_dir", file.path(bpf_combine_result_dir(repo_root), "pooled_embedding_curve_class"))
      part_args <- append_default_arg(part_args, "mode", "report")
      run_rscript(
        "combine_pooled_embedding_curve_class_report",
        canonical_combined_runner,
        part_args,
        dry_run
      )
    } else if (identical(part, "fixo2_eigen_attractor")) {
      out_dir <- file.path(bpf_combine_fixo2_eigen_attractor_result_dir(repo_root), "pooled_embedding_curve_class")
      part_args <- append_default_arg(args, "out_dir", out_dir)
      part_args <- append_default_arg(part_args, "run_report", "false")
      run_rscript(
        "combine_fixo2_eigen_attractor_curve_class",
        canonical_fixo2_runner,
        part_args,
        dry_run
      )
    } else if (identical(part, "fixo2_eigen_report")) {
      part_args <- append_default_arg(
        args,
        "embedding_dir",
        file.path(bpf_combine_fixo2_eigen_attractor_result_dir(repo_root), "pooled_embedding_curve_class")
      )
      part_args <- append_default_arg(
        part_args,
        "classification_dir",
        file.path(bpf_dense_grid_result_root(repo_root), "dense-grid_monotonicity_classification")
      )
      run_rscript(
        "combine_fixo2_eigen_attractor_report",
        file.path(WORKFLOW_ROOT, "report", "combined_fixo2_eigen", "render_fixo2_eigen_attractor_embedding_curve_class_report.R"),
        part_args,
        dry_run
      )
    }
  }
}

run_combine_report <- function(args, repo_root, dry_run) {
  canonical_combined_runner <- file.path(WORKFLOW_ROOT, "runner", "combined_parameter_landscape", "run_combined_parameter_landscape.R")
  args <- append_default_arg(args, "out_dir", file.path(bpf_combine_result_dir(repo_root), "pooled_embedding_curve_class"))
  args <- append_default_arg(args, "mode", "report")
  run_rscript(
    "combine_pooled_embedding_curve_class_report",
    canonical_combined_runner,
    args,
    dry_run
  )
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  argv <- bpf_parse_args(raw_args)
  if (bpf_as_bool(argv$help %||% argv$h, FALSE)) {
    usage()
    return(invisible(NULL))
  }

  repo_root <- bpf_repo_root(SCRIPT_DIR)
  dry_run <- bpf_as_bool(argv$dry_run, FALSE)
  workflows <- normalize_workflows(argv$workflow)
  forward_args <- drop_forward_keys(
    raw_args,
    c("workflow", "parameter_parts", "dense_parts", "combine_parts", "help", "h")
  )

  for (workflow in workflows) {
    if (identical(workflow, "fixed_o2")) {
      run_fixed_o2(forward_args, repo_root, dry_run)
    } else if (identical(workflow, "parameter_landscape")) {
      run_parameter_landscape(forward_args, argv, repo_root, dry_run)
    } else if (identical(workflow, "dense_grid_monotonicity")) {
      run_dense_grid(forward_args, argv, repo_root, dry_run)
    } else if (identical(workflow, "combine")) {
      run_combine(forward_args, argv, repo_root, dry_run)
    } else if (identical(workflow, "combine_report")) {
      run_combine_report(forward_args, repo_root, dry_run)
    }
  }
  invisible(NULL)
}

if (identical(environment(), globalenv())) {
  main()
}
