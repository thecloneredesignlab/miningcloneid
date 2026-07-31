#!/usr/bin/env Rscript

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!grepl("^--[^=]+=", arg)) next
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    out[[key]] <- sub("^[^=]+=", "", arg)
  }
  out
}

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) y else x

main <- function(argv = parse_args()) {
  repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  source_path <- argv$source %||% paste0(
    "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/analysis/",
    "best_fit_parameter_feature/04_combine_parameter_landscape/pooled_embedding_curve_class/",
    "tables/TSNEs/Full/",
    "pooled_invivo_invitro_initial_vs_best_context_prior_unit_tsne_best_points_curve_class.csv"
  )
  output_path <- argv$output %||% file.path(
    repo_root, "oxygen", "results", "figure6_fixed_o2_cin_map",
    "invivo_parameter_landscape_tsne_curve_class.csv"
  )
  if (!file.exists(source_path)) stop("Missing parameter-landscape source: ", source_path, call. = FALSE)

  tab <- utils::read.csv(source_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c(
    "tSNE1", "tSNE2", "dataset", "point_type", "source_group", "seed",
    "objective", "reduction", "seed_number", "curve_class", "average_slope"
  )
  if (!all(required %in% names(tab))) {
    stop("Parameter-landscape source is missing: ", paste(setdiff(required, names(tab)), collapse = ", "), call. = FALSE)
  }
  out <- tab[tab$dataset == "invivo" & tab$point_type == "best", required, drop = FALSE]
  out <- out[order(as.integer(out$seed)), , drop = FALSE]
  if (nrow(out) != 500L || anyDuplicated(out$seed) || !all(is.finite(out$tSNE1)) ||
      !all(is.finite(out$tSNE2)) || !all(is.finite(out$objective)) || anyNA(out$curve_class)) {
    stop("Expected 500 unique, fully annotated in vivo best-fit landscape points.", call. = FALSE)
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, output_path, row.names = FALSE, quote = TRUE)
  message("Wrote ", output_path, " [", nrow(out), " in vivo best-fit points]")
}

if (sys.nframe() == 0L) main()
