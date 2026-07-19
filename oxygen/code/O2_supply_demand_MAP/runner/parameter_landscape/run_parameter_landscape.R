#!/usr/bin/env Rscript

.o2pl_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "run_parameter_landscape.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
})
.o2pl_runner_root <- normalizePath(file.path(.o2pl_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_runner_root, "util", "o2_supply_demand_map_parameter_landscape_io_utils.R"), local = environment(), chdir = TRUE)

o2pl_cli_arg <- function(name, value) {
  if (is.null(value) || !length(value) || all(is.na(value))) return(character())
  paste0("--", name, "=", paste(value, collapse = ","))
}

o2pl_forward_args <- function(argv, exclude = character()) {
  names_to_use <- setdiff(names(argv), exclude)
  unlist(lapply(names_to_use, function(name) o2pl_cli_arg(name, argv[[name]])), use.names = FALSE)
}

o2pl_run_child <- function(label, script, args = character(), dry_run = FALSE) {
  script <- normalizePath(script, mustWork = FALSE)
  if (!file.exists(script)) stop("Missing workflow script: ", script, call. = FALSE)
  command <- c("--vanilla", script, args)
  message("[", label, "] ", paste(shQuote(c(file.path(R.home("bin"), "Rscript"), command)), collapse = " "))
  if (dry_run) return(invisible(0L))
  status <- system2(file.path(R.home("bin"), "Rscript"), command)
  if (!identical(status, 0L)) stop(label, " failed with status ", status, call. = FALSE)
  invisible(status)
}

o2pl_normalize_run_parts <- function(value) {
  parts <- tolower(as_char_vec(value, "all"))
  aliases <- c(
    invivo_table = "invivo_tables",
    invitro_table = "invitro_tables",
    invivo_reduction = "invivo_reductions",
    invitro_reduction = "invitro_reductions",
    pooled_reduction = "pooled_reductions",
    clustering_report = "reports",
    report = "reports",
    contribution = "mode_contribution"
  )
  parts <- vapply(parts, function(part) if (part %in% names(aliases)) aliases[[part]] else part, character(1))
  if ("all" %in% parts) parts <- c("invivo_tables", "invitro_tables", "invivo_reductions", "invitro_reductions", "pooled_reductions", "reports")
  unique(parts)
}

o2pl_parameter_landscape_runner_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root_dir <- normalizePath(path.expand(argv$result_root %||% o2pl_default_result_root()), mustWork = FALSE)
  dry_run <- as_bool(argv$dry_run, FALSE)
  requested <- argv$run_parts %||% argv$analysis_part %||% argv$part %||% "all"
  parts <- o2pl_normalize_run_parts(requested)
  allowed <- c("invivo_tables", "invitro_tables", "invivo_reductions", "invitro_reductions", "pooled_reductions", "mode_contribution", "dominant_ploidy_contribution", "reports")
  unknown <- setdiff(parts, allowed)
  if (length(unknown)) stop("Unknown run part(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  forwarded <- o2pl_forward_args(argv, c("run_parts", "analysis_part", "part", "dry_run", "result_root", "invivo_input", "invitro_input", "input_dir"))
  simulation_script <- file.path(.o2pl_runner_root, "simulation", "parameter_landscape", "generate_parameter_landscape_simulation_tables.R")
  analysis_script <- file.path(.o2pl_runner_root, "analysis", "parameter_landscape_clustering", "analyze_parameter_landscape.R")
  contribution_script <- file.path(.o2pl_runner_root, "analysis", "parameter_landscape_clustering", "parameter_contribution_analysis.R")
  visualization_script <- file.path(.o2pl_runner_root, "vis", "parameter_landscape", "visualize_parameter_landscape.R")
  report_dir <- file.path(.o2pl_runner_root, "report", "parameter_landscape")
  common <- c(o2pl_cli_arg("result_root", root_dir), forwarded)

  if ("invivo_tables" %in% parts) {
    o2pl_run_child("materialize invivo parameter landscape", simulation_script, c(
      o2pl_cli_arg("dataset", "invivo"),
      o2pl_cli_arg("input_dir", argv$invivo_input %||% argv$input_dir %||% DEFAULT_INVIVO_INPUT_DIR),
      o2pl_cli_arg("result_root", root_dir),
      o2pl_cli_arg("write_growth_turnover", argv$write_growth_turnover %||% "TRUE"),
      o2pl_cli_arg("write_modes", argv$write_modes %||% "TRUE"),
      forwarded
    ), dry_run)
  }
  if ("invitro_tables" %in% parts) {
    o2pl_run_child("materialize invitro parameter landscape", simulation_script, c(
      o2pl_cli_arg("dataset", "invitro"),
      o2pl_cli_arg("input_dir", argv$invitro_input %||% argv$input_dir %||% DEFAULT_INVITRO_INPUT_DIR),
      o2pl_cli_arg("result_root", root_dir),
      o2pl_cli_arg("write_modes", "FALSE"),
      forwarded
    ), dry_run)
  }
  for (scope in intersect(c("invivo", "invitro", "pooled"), sub("_reductions$", "", parts[grepl("_reductions$", parts)]))) {
    o2pl_run_child(paste("analyze", scope, "parameter landscape"), analysis_script, c(common, o2pl_cli_arg("part", scope)), dry_run)
    o2pl_run_child(paste("visualize", scope, "parameter landscape"), visualization_script, c(common, o2pl_cli_arg("part", scope)), dry_run)
  }
  for (target in c("mode", "dominant_ploidy")) {
    token <- paste0(target, "_contribution")
    if (!token %in% parts) next
    contribution_args <- c(common, o2pl_cli_arg("mode_contribution_target", target))
    o2pl_run_child(paste("analyze", target, "parameter contribution"), contribution_script, contribution_args, dry_run)
    o2pl_run_child(paste("visualize", target, "parameter contribution"), visualization_script, c(contribution_args, o2pl_cli_arg("part", "contribution")), dry_run)
    if (as_bool(argv$render_html_report %||% argv$run_report, TRUE)) {
      report_script <- file.path(report_dir, if (target == "mode") "mode_parameter_contribution_report.R" else "dominant_ploidy_parameter_contribution_report.R")
      o2pl_run_child(paste("assemble", target, "parameter report"), report_script, contribution_args, dry_run)
    }
  }
  if ("reports" %in% parts) {
    o2pl_run_child("assemble parameter-landscape clustering report", file.path(report_dir, "clustering_report.R"), common, dry_run)
  }
  message("Parameter-landscape runner complete: ", root_dir)
  invisible(TRUE)
}

if (sys.nframe() == 0L) o2pl_parameter_landscape_runner_main()
