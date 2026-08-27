#!/usr/bin/env Rscript

# Backward-compatible loader. Canonical implementation lives in workflow util/.
.bpf_cli_wrapper_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frame_files[basename(frame_files) == "cli_utils.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
    if (basename(path) == "cli_utils.R") return(dirname(path))
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.bpf_cli_canonical <- normalizePath(
  file.path(.bpf_cli_wrapper_dir, "..", "..", "..", "util", "o2_supply_demand_map_bpf_cli_utils.R"),
  mustWork = TRUE
)
source(.bpf_cli_canonical, local = environment(), chdir = TRUE)
rm(.bpf_cli_wrapper_dir, .bpf_cli_canonical)
