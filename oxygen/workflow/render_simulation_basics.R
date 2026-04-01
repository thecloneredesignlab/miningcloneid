args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1]])
} else {
  "workflow/render_simulation_basics.R"
}

repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
rmarkdown::render(
  input = file.path(repo_root, "workflow", "simulation_basics.Rmd"),
  output_dir = file.path(repo_root, "workflow"),
  quiet = TRUE
)
