args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1]])
} else {
  "workflow/render_joint_fit_report.R"
}

parse_args <- function(args) {
  out <- list()
  if (!length(args)) return(out)
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
argv <- parse_args(commandArgs(trailingOnly = TRUE))

as_logical_arg <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(default)
  tolower(trimws(as.character(x[[1]]))) %in% c("true", "t", "1", "yes", "y")
}

as_num_arg <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) return(default)
  out <- suppressWarnings(as.numeric(x[[1]]))
  if (is.finite(out)) out else default
}

render_params <- list(
  parameter_object_path = argv$parameter_object_path %||% "workflow/joint_fit_parameter_object.R",
  selected_parameter_set = argv$selected_parameter_set %||% NULL,
  invivo_fit_dir = argv$invivo_fit_dir %||% NULL,
  render_invivo_viz = as_logical_arg(argv$render_invivo_viz, TRUE),
  invivo_report_dt = as_num_arg(argv$invivo_report_dt, 1),
  invivo_top_n = as.integer(as_num_arg(argv$invivo_top_n, 6)),
  invivo_predict_horizons = argv$invivo_predict_horizons %||% "1000",
  invivo_max_scenarios = as_num_arg(argv$invivo_max_scenarios, Inf)
)

out_dir <- normalizePath(argv$output_dir %||% file.path(repo_root, "workflow"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- argv$output_file %||% "joint_fit_report.html"

rmarkdown::render(
  input = file.path(repo_root, "workflow", "joint_fit_report.Rmd"),
  output_dir = out_dir,
  output_file = out_file,
  params = render_params,
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
