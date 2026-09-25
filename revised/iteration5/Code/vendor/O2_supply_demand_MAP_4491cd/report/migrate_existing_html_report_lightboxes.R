#!/usr/bin/env Rscript

# Add the shared image-viewing layer to existing HTML reports without changing
# their scientific content, figures, captions, navigation, or embedded images.

.o2sd_migration_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    file_path <- sub("^--file=", "", file_arg)
    own_arg <- file_path[basename(file_path) == "migrate_existing_html_report_lightboxes.R"]
    if (length(own_arg)) {
      return(dirname(normalizePath(own_arg[[1L]], mustWork = FALSE)))
    }
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
    }, character(1L))
  )
  own <- frame_files[basename(frame_files) == "migrate_existing_html_report_lightboxes.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2sd_migration_workflow_root <- normalizePath(file.path(.o2sd_migration_script_dir, ".."), mustWork = TRUE)
source(
  file.path(.o2sd_migration_workflow_root, "util", "o2_supply_demand_map_html_utils.R"),
  local = environment()
)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

o2sd_lightbox_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
  }
  out
}

o2sd_lightbox_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || !nzchar(trimws(x[[1L]]))) return(default)
  value <- tolower(trimws(x[[1L]]))
  if (value %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Expected a Boolean value, received: ", x[[1L]], call. = FALSE)
}

o2sd_repo_root_from_report_dir <- function() {
  normalizePath(file.path(.o2sd_migration_workflow_root, "..", "..", ".."), mustWork = TRUE)
}

o2sd_migrate_existing_html_report_lightboxes <- function(report_root, manifest_path, dry_run = FALSE) {
  report_root <- normalizePath(path.expand(report_root), mustWork = TRUE)
  html_paths <- sort(list.files(report_root, pattern = "[.]html$", full.names = TRUE, recursive = TRUE))
  rows <- vector("list", length(html_paths))

  for (i in seq_along(html_paths)) {
    path <- html_paths[[i]]
    before <- as.numeric(file.info(path)$size)
    result <- tryCatch({
      if (isTRUE(dry_run)) {
        text <- local({
          con <- file(path, open = "rb")
          on.exit(close(con), add = TRUE)
          readChar(con, nchars = before, useBytes = TRUE)
        })
        status <- if (grepl("id=\"o2sd-report-image-lightbox-script\"", text, fixed = TRUE)) {
          "already_current"
        } else if (!grepl("<img", text, fixed = TRUE)) {
          "skipped_no_images"
        } else {
          "would_inject"
        }
        list(status = status, changed = FALSE)
      } else {
        o2sd_inject_report_image_lightbox(path)
      }
    }, error = function(e) list(status = "failed", changed = FALSE, error = conditionMessage(e)))
    after <- as.numeric(file.info(path)$size)
    rows[[i]] <- data.frame(
      report_path = substring(path, nchar(report_root) + 2L),
      status = result$status,
      changed = isTRUE(result$changed),
      bytes_before = before,
      bytes_after = after,
      error = result$error %||% "",
      stringsAsFactors = FALSE
    )
  }

  manifest <- if (length(rows)) do.call(rbind, rows) else data.frame(
    report_path = character(), status = character(), changed = logical(),
    bytes_before = numeric(), bytes_after = numeric(), error = character()
  )
  manifest_path <- normalizePath(path.expand(manifest_path), mustWork = FALSE)
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  manifest
}

main <- function() {
  args <- o2sd_lightbox_parse_args()
  repo_root <- o2sd_repo_root_from_report_dir()
  report_root <- args$report_root %||% file.path(repo_root, "oxygen", "results")
  manifest_path <- args$manifest %||% file.path(
    repo_root, "oxygen", "results", "analysis", "report_image_lightbox_migration.tsv"
  )
  result <- o2sd_migrate_existing_html_report_lightboxes(
    report_root = report_root,
    manifest_path = manifest_path,
    dry_run = o2sd_lightbox_as_bool(args$dry_run, FALSE)
  )
  counts <- sort(table(result$status), decreasing = TRUE)
  message("HTML reports scanned: ", nrow(result))
  if (length(counts)) message(paste(names(counts), as.integer(counts), sep = "=", collapse = ", "))
  message("Migration manifest: ", normalizePath(manifest_path, mustWork = TRUE))
  if (any(result$status == "failed")) stop("One or more HTML reports could not be migrated.", call. = FALSE)
}

if (sys.nframe() == 0L) main()
