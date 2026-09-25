#!/usr/bin/env Rscript
# Deprecated compatibility loader; slope statistics live in analysis/combined_parameter_landscape.
.legacy_file <- local({ frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1))); own <- frames[basename(frames) == "calculate_regression_curve_average_slope.R"]; if (length(own)) own[[length(own)]] else normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE) })
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
.combined_analysis_source_override <- deprecated_combined_analysis()
sys.source(deprecated_combined_analysis(), envir = environment(), chdir = FALSE)
if (sys.nframe() == 0L) {
  raw <- commandArgs(trailingOnly = TRUE)
  out_arg <- grep("^--out[_-]dir=", raw, value = TRUE)
  if (length(out_arg)) {
    old_out <- sub("^--out[_-]dir=", "", out_arg[[1L]])
    if (identical(basename(normalizePath(old_out, mustWork = FALSE)), "tables")) {
      raw <- raw[!grepl("^--out[_-]dir=", raw)]
      raw <- c(raw, paste0("--out_dir=", dirname(old_out)))
      if (!any(grepl("^--output[_-]file=", raw))) raw <- c(raw, paste0("--output_file=", file.path(old_out, "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv")))
    }
  }
  deprecated_combined_exec(raw, "analyze", "--slope_only=true")
}
