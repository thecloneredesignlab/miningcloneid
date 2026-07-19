#!/usr/bin/env Rscript

.bpf_util_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "o2_supply_demand_map_bpf_report_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(.bpf_util_dir, "o2_supply_demand_map_html_utils.R"), local = environment())
rm(.bpf_util_dir)

bpf_rel_path <- function(path, from_dir) {
  path <- normalizePath(path, mustWork = FALSE)
  from_dir <- normalizePath(from_dir, mustWork = FALSE)
  if (requireNamespace("tools", quietly = TRUE)) {
    return(tools::file_path_as_absolute(path))
  }
  path
}

bpf_html_escape <- o2sd_html_escape_standard
