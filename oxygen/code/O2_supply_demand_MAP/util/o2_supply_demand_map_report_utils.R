#!/usr/bin/env Rscript

# Side-effect-free report support shared by report assemblers.  These helpers
# read already-materialized tables, normalize display values, and prepare
# existing files for embedding.  They do not calculate scientific results,
# draw figures, or assemble a report.

o2sd_report_read_table_optional <- function(path, sep = "\t") {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.table(
      path,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
}

o2sd_report_safe_read_delim <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 1L) return(data.frame())
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

o2sd_report_parse_equals_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    m <- regmatches(arg, regexec("^--([^=]+)=(.*)$", arg))[[1]]
    if (length(m)) out[[m[[2]]]] <- m[[3]]
  }
  out
}

o2sd_report_null_coalesce_simple <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

o2sd_report_null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

o2sd_report_truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y", "on")
}

o2sd_report_format_numeric_like <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  out <- as.character(x)
  x_trim <- trimws(out)
  numeric_pattern <- "^[-+]?((\\d+\\.?\\d*)|(\\.\\d+))([eE][-+]?\\d+)?$"
  numeric_like <- !is.na(out) & nzchar(x_trim) & grepl(numeric_pattern, x_trim)
  if (!any(numeric_like)) return(out)

  num <- suppressWarnings(as.numeric(x_trim[numeric_like]))
  keep <- is.finite(num)
  if (!any(keep)) return(out)

  out_num <- as.character(num[keep])
  int_like <- abs(num[keep] - round(num[keep])) < 1e-9
  decimal_text <- formatC(num[keep], format = "f", digits = 3)
  sci_nonzero <- !int_like & num[keep] != 0 & grepl("\\.000$", decimal_text)
  decimal <- !int_like & !sci_nonzero

  out_num[int_like] <- format(round(num[keep][int_like]), scientific = FALSE, trim = TRUE, digits = 15)
  out_num[sci_nonzero] <- formatC(num[keep][sci_nonzero], format = "e", digits = 3)
  out_num[decimal] <- formatC(num[keep][decimal], format = "f", digits = 3)

  idx <- which(numeric_like)[keep]
  out[idx] <- out_num
  out
}

o2sd_report_transformed_param_to_natural <- function(x) {
  x <- as.character(x)
  x <- sub("^ivt__", "", x)
  x <- sub("^log10_", "", x)
  x <- sub("^logit_", "", x)
  x
}

o2sd_report_annotate_parameter_table_for_report <- function(tab, fit_dir) {
  o2sd_add_parameter_descriptions(
    tab,
    parameter_tables = parameter_description_table_paths(fit_dir)
  )
}

o2sd_report_pandoc_available <- function() {
  isTRUE(rmarkdown::pandoc_available("1.12.3"))
}

o2sd_report_pdflatex_available <- function() {
  if (identical(Sys.getenv("O2_REPORT_FORCE_NO_PDFLATEX", unset = ""), "TRUE")) {
    return(FALSE)
  }
  nzchar(Sys.which("pdflatex"))
}

o2sd_report_magick_available <- function() {
  if (identical(Sys.getenv("O2_REPORT_FORCE_NO_MAGICK", unset = ""), "TRUE")) {
    return(FALSE)
  }
  requireNamespace("magick", quietly = TRUE)
}

o2sd_report_ghostscript_bin <- function() {
  candidates <- c(Sys.which("gs"), Sys.which("gswin64c"), Sys.which("gswin32c"))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

o2sd_report_gs_available <- function() {
  nzchar(o2sd_report_ghostscript_bin())
}

o2sd_report_base64enc_available <- function() {
  requireNamespace("base64enc", quietly = TRUE)
}

o2sd_report_render_pdf_preview_png_gs <- function(src_pdf, dest_png, density = 180) {
  gs_bin <- o2sd_report_ghostscript_bin()
  if (!nzchar(gs_bin)) {
    stop("Ghostscript ('gs') was requested for PDF preview rendering but is not available in PATH.")
  }
  src_pdf_use <- normalizePath(src_pdf, mustWork = TRUE)
  dest_png_use <- normalizePath(dest_png, mustWork = FALSE)
  density_use <- suppressWarnings(as.integer(density))
  if (!is.finite(density_use) || density_use <= 0L) density_use <- 180L
  args <- c(
    "-dSAFER",
    "-dBATCH",
    "-dNOPAUSE",
    "-sDEVICE=pngalpha",
    sprintf("-r%d", density_use),
    sprintf("-sOutputFile=%s", shQuote(dest_png_use)),
    shQuote(src_pdf_use)
  )
  out <- suppressWarnings(system2(gs_bin, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "Ghostscript failed while rendering PDF preview for ", src_pdf, ": ",
      paste(out, collapse = "\n")
    )
  }
  if (!file.exists(dest_png_use)) {
    stop("Ghostscript did not create expected PNG preview: ", dest_png_use)
  }
  normalizePath(dest_png_use, mustWork = TRUE)
}

o2sd_report_file_to_data_uri <- function(path, mime) {
  if (o2sd_report_base64enc_available()) {
    return(base64enc::dataURI(file = path, mime = mime))
  }
  base64_bin <- Sys.which("base64")
  if (nzchar(base64_bin)) {
    enc <- tryCatch(
      suppressWarnings(system2(base64_bin, c("-w", "0", path), stdout = TRUE, stderr = TRUE)),
      error = function(e) character()
    )
    if (!length(enc)) {
      enc <- tryCatch(
        suppressWarnings(system2(base64_bin, path, stdout = TRUE, stderr = TRUE)),
        error = function(e) character()
      )
    }
    if (length(enc) > 0L) {
      return(sprintf("data:%s;base64,%s", mime, paste(enc, collapse = "")))
    }
  }
  stop(
    "HTML report fallback requires either the R package 'base64enc' or a system 'base64' command ",
    "when 'magick' is unavailable."
  )
}
