#!/usr/bin/env Rscript
# Deprecated compatibility loader; classification functions live in util.
.legacy_file <- local({ frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1))); own <- frames[basename(frames) == "curve_classification_utils.R"]; if (length(own)) own[[length(own)]] else normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE) })
.root <- normalizePath(file.path(dirname(.legacy_file), "..", ".."), mustWork = TRUE)
source(file.path(.root, "util", "o2_supply_demand_map_curve_classification_utils.R"), local = environment(), chdir = TRUE)
