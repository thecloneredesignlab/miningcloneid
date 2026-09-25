#!/usr/bin/env Rscript

# Shared file-contract helpers for multi-warmup stages.  No workflow execution,
# numerical simulation, statistical analysis, plotting, or report assembly.

o2mw_util_dir <- function() {
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "o2_supply_demand_map_multi_warmup_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
}

o2mw_workflow_root <- function() normalizePath(file.path(o2mw_util_dir(), ".."), mustWork = FALSE)

source(
  file.path(o2mw_util_dir(), "o2_supply_demand_map_shared.R"),
  local = TRUE,
  chdir = TRUE
)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) return(y)
  if (length(x) == 1L && (is.na(x) || !nzchar(trimws(as.character(x))))) return(y)
  x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pair[[1L]], fixed = TRUE)
    out[[key]] <- if (length(pair) > 1L) paste(pair[-1L], collapse = "=") else "TRUE"
  }
  out
}

as_chr <- function(x, default = "") as.character((x %||% default)[[1L]])

as_bool <- function(x, default = FALSE) {
  value <- tolower(trimws(as_chr(x, if (default) "true" else "false")))
  if (value %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  value <- suppressWarnings(as.integer(as_chr(x, as.character(default))))
  if (is.na(value)) default else value
}

as_num <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(as_chr(x, as.character(default))))
  if (!is.finite(value)) default else value
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(default)
  values <- suppressWarnings(as.numeric(trimws(strsplit(
    paste(as.character(x), collapse = ","),
    ",",
    fixed = TRUE
  )[[1L]])))
  values <- values[is.finite(values)]
  if (length(values)) values else default
}

as_char_vec <- function(x, default = character()) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(default)
  values <- trimws(strsplit(
    paste(as.character(x), collapse = ","),
    ",",
    fixed = TRUE
  )[[1L]])
  values <- values[nzchar(values)]
  if (length(values)) values else default
}

num_col <- o2sd_numeric

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

read_table_optional <- function(path, sep = "\t") {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    utils::read.table(path, sep = sep, header = TRUE, quote = "", comment.char = "",
                      check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

rbind_fill <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x), rows)
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

read_csv_plain <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

rbind_fill_plain <- rbind_fill

integrated_extra_results_dir <- function(root) {
  file.path(root, "multi_warmup_integrated_joint_run", "extra_results")
}

warmup_levels_from_manifest <- function(root, fallback = character(0)) {
  manifest <- read_table_optional(file.path(root, "multi_warmup_manifest.tsv"))
  if (is.data.frame(manifest) && nrow(manifest) && "warmup_label" %in% names(manifest)) {
    out <- unique(toupper(sub("_.*$", "", as.character(manifest$warmup_label))))
    out <- out[nzchar(out)]
    if (length(out)) return(out)
  }
  unique(as.character(fallback[nzchar(fallback)]))
}

add_invivo_warmup_from_seed <- function(tab, root, seed_col = "seed") {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab) || !(seed_col %in% names(tab))) return(tab)
  tab$warmup_label <- sub("__.*$", "", as.character(tab[[seed_col]]))
  tab$invivo_warmup <- toupper(sub("_.*$", "", tab$warmup_label))
  levels <- warmup_levels_from_manifest(root, fallback = tab$invivo_warmup)
  tab$invivo_warmup <- factor(tab$invivo_warmup, levels = levels)
  tab
}

reduction_file_suffix <- function(reduction) normalize_reduction(reduction)

format_heatmap_value <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.finite(x),
    ifelse(
      abs(x) >= 1000 | abs(x) < 0.01,
      formatC(x, format = "e", digits = 2),
      format(signif(x, 3), trim = TRUE)
    ),
    ""
  )
}

sanitize_label <- function(x, fallback = "warmup") {
  value <- gsub("[^A-Za-z0-9_.-]+", "_", trimws(as.character(x)))
  value <- gsub("^_+|_+$", "", value)
  if (!length(value) || !nzchar(value[[1L]])) fallback else value
}

seed_order_key <- function(x) {
  x <- as.character(x)
  value <- suppressWarnings(as.numeric(sub("^seed", "", x)))
  value[!is.finite(value)] <- Inf
  value
}

family_palette <- function(levels) {
  levels <- as.character(levels)
  levels <- levels[nzchar(levels)]
  palette <- c(
    "#009E73",
    "#6A3D9A",
    "#E69F00",
    "#000000",
    "#8C6D31",
    "#F0E442",
    "#66A61E",
    "#7F7F7F",
    "#B07AA1",
    "#A6761D"
  )
  if (length(levels) > length(palette)) {
    palette <- grDevices::colorRampPalette(palette)(length(levels))
  } else {
    palette <- palette[seq_along(levels)]
  }
  stats::setNames(palette, levels)
}

safe_component <- function(x, fallback = "item") {
  value <- gsub("[^A-Za-z0-9_.-]+", "_", trimws(as.character(x)))
  value <- gsub("^_+|_+$", "", value)
  ifelse(nzchar(value), value, fallback)
}
