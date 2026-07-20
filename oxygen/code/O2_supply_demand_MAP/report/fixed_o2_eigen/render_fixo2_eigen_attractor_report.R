#!/usr/bin/env Rscript

# Report-only assembly from materialized feature, analysis, and figure files.

.fixo2ea_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_fixo2_eigen_attractor_report.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  workflow_frames <- frames[
    grepl("/O2_supply_demand_MAP/", frames, fixed = TRUE)
  ]
  if (length(workflow_frames)) {
    root <- dirname(workflow_frames[[length(workflow_frames)]])
    while (!identical(basename(root), "O2_supply_demand_MAP")) {
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (identical(basename(root), "O2_supply_demand_MAP")) {
      return(file.path(root, "report", "fixed_o2_eigen"))
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.fixo2ea_report_workflow_root <- normalizePath(file.path(.fixo2ea_report_dir, "..", ".."), mustWork = TRUE)
source(file.path(.fixo2ea_report_workflow_root, "util", "o2_supply_demand_map_eigen_attractor_utils.R"), local = environment())
source(file.path(.fixo2ea_report_workflow_root, "util", "o2_supply_demand_map_html_utils.R"), local = environment())

fixo2ea_html_escape <- o2sd_html_escape_standard

fixo2ea_file_shape <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(path = path, exists = FALSE, rows = NA_integer_, cols = NA_integer_, stringsAsFactors = FALSE))
  }
  df <- tryCatch(utils::read.csv(path, nrows = 5L, check.names = FALSE), error = function(e) NULL)
  cols <- if (is.null(df)) NA_integer_ else ncol(df)
  rows <- tryCatch({
    max(0L, length(readLines(path, warn = FALSE)) - 1L)
  }, error = function(e) NA_integer_)
  data.frame(path = path, exists = TRUE, rows = rows, cols = cols, stringsAsFactors = FALSE)
}

fixo2ea_html_table <- function(df) {
  if (is.null(df) || !nrow(df)) return("<p>No rows.</p>")
  header <- paste0("<tr>", paste0("<th>", fixo2ea_html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  body <- apply(df, 1L, function(row) {
    paste0("<tr>", paste0("<td>", fixo2ea_html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table>", header, paste(body, collapse = "\n"), "</table>")
}

fixo2ea_write_report <- function(result_root = fixo2ea_default_result_root(),
                                 output = file.path(result_root, "fixo2_eigen_attractor_based_clustering_report.html")) {
  result_root <- normalizePath(path.expand(result_root), mustWork = FALSE)
  feature_paths <- file.path(
    result_root,
    "FixO2EigenAttractors", "Tables",
    c(
      "invivo_best_fixo2_eigen_attractor_wide.csv",
      "invivo_initial_fixo2_eigen_attractor_wide.csv",
      "invitro_best_fixo2_eigen_attractor_wide.csv",
      "invitro_initial_fixo2_eigen_attractor_wide.csv",
      "fixo2_eigen_feature_manifest.csv"
    )
  )
  feature_shapes <- do.call(rbind, lapply(feature_paths, fixo2ea_file_shape))
  feature_shapes$path <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", result_root), "/?"), "", feature_shapes$path)

  coord_files <- list.files(result_root, pattern = "_coordinates\\.csv$", recursive = TRUE, full.names = TRUE)
  coord_shapes <- if (length(coord_files)) {
    out <- do.call(rbind, lapply(coord_files, fixo2ea_file_shape))
    out$path <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", result_root), "/?"), "", out$path)
    out
  } else {
    data.frame(path = character(), exists = logical(), rows = integer(), cols = integer(), stringsAsFactors = FALSE)
  }

  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\">",
    "<title>FixO2 Eigen Attractor Based Clustering</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:32px;line-height:1.45;color:#1f2933}",
    "h1,h2{font-weight:650} table{border-collapse:collapse;margin:12px 0 28px 0;font-size:13px}",
    "th,td{border:1px solid #d7dde5;padding:6px 8px;text-align:left} th{background:#eef2f6}",
    "code{background:#f4f6f8;padding:1px 4px;border-radius:3px}",
    "</style></head><body>",
    "<h1>FixO2 Eigen Attractor Based Clustering</h1>",
    "<p>Feature input is the wide matrix of dominant eigen-attractor mean ploidy across fixed-O2 values. Initial rows are full DEoptim initial populations, not sampled subsets.</p>",
    "<h2>Feature Tables</h2>",
    fixo2ea_html_table(feature_shapes),
    "<h2>Reduction Coordinate Tables</h2>",
    fixo2ea_html_table(coord_shapes),
    "</body></html>"
  )
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output, useBytes = TRUE)
  o2sd_inject_report_image_lightbox(output)
  message("Wrote ", output)
  invisible(output)
}
