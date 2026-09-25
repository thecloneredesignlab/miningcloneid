#!/usr/bin/env Rscript

# Freeze/verify Figure 6 assets outside the requested A/B arithmetic-mean
# change. The main C-F panel objects and every Supplementary Figure 6 data and
# rendered asset must remain byte-identical during the A/B-only run.

options(stringsAsFactors = FALSE, warn = 1)

f6abg_workspace_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "Code", "Figures")) &&
        dir.exists(file.path(current, "data", "Figures"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Cannot locate iteration4 workspace.")
    current <- parent
  }
}

f6abg_source_freeze_helpers <- function(workspace_root) {
  helper <- file.path(
    workspace_root, "Code", "Figures", "util", "analysis",
    "figure6_checksum_helpers.R"
  )
  if (!file.exists(helper)) stop("Missing Figure 6 freeze helper: ", helper)
  sys.source(helper, envir = globalenv())
}

f6abg_current_run_root <- function(figure6_root, manifest_name) {
  manifest_path <- file.path(figure6_root, manifest_name)
  if (!file.exists(manifest_path)) stop("Missing current-run manifest: ", manifest_path)
  manifest <- utils::read.delim(manifest_path, check.names = FALSE)
  if (nrow(manifest) != 1L || !"relative_run_path" %in% names(manifest)) {
    stop("Invalid current-run manifest: ", manifest_path)
  }
  run_root <- file.path(figure6_root, manifest$relative_run_path[[1L]])
  if (!dir.exists(run_root)) stop("Missing current run directory: ", run_root)
  normalizePath(run_root, mustWork = TRUE)
}

f6abg_non_ab_files <- function(workspace_root = f6abg_workspace_root()) {
  figure6_root <- file.path(
    workspace_root, "data", "Figures", "Figure6", "fixed_pmisseg_v1"
  )
  finite_root <- f6abg_current_run_root(
    figure6_root, "finite_time_q10_current.tsv"
  )
  passage_root <- f6abg_current_run_root(
    figure6_root, "invitro_passage_q10_current.tsv"
  )
  required <- c(
    file.path(figure6_root, "finite_time_q10_current.tsv"),
    file.path(finite_root, "finite_time_panel_c.rds"),
    file.path(finite_root, "finite_time_panel_d.rds"),
    file.path(figure6_root, "invitro_passage_q10_current.tsv"),
    file.path(passage_root, "passage_time_panel_e.rds"),
    file.path(passage_root, "passage_time_panel_f.rds")
  )
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Missing Figure 6C-F freeze input: ", paste(missing, collapse = ", "))
  }
  supplementary_data_roots <- list.dirs(
    file.path(workspace_root, "data", "Figures"),
    recursive = FALSE, full.names = TRUE
  )
  supplementary_data_roots <- supplementary_data_roots[
    grepl("/Supp_Figure6_[0-9]+$", supplementary_data_roots)
  ]
  supplementary_data <- unlist(lapply(
    supplementary_data_roots,
    list.files, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE, include.dirs = FALSE
  ), use.names = FALSE)
  rendered_roots <- file.path(
    workspace_root, c("Figures", "manuscript/Figures")
  )
  rendered <- unlist(lapply(rendered_roots[dir.exists(rendered_roots)], function(root) {
    files <- list.files(root, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    files[grepl("^supp_fig6-", basename(files)) & !dir.exists(files)]
  }), use.names = FALSE)
  sort(unique(c(required, supplementary_data, rendered)))
}

f6abg_inventory <- function(workspace_root = f6abg_workspace_root()) {
  f6abg_source_freeze_helpers(workspace_root)
  repository_root <- f6g_repository_root(workspace_root)
  files <- f6abg_non_ab_files(workspace_root)
  if (!length(files)) stop("No non-A/B Figure 6 files found for freeze inventory.")
  relative_paths <- vapply(
    files, f6g_relative, character(1L), root = repository_root
  )
  info <- file.info(files)
  data.frame(
    relative_path = relative_paths,
    size_bytes = as.numeric(info$size),
    sha256 = f6g_sha256(files),
    stringsAsFactors = FALSE
  )
}

f6abg_freeze <- function(workspace_root = f6abg_workspace_root(), output) {
  f6abg_source_freeze_helpers(workspace_root)
  inventory <- f6abg_inventory(workspace_root)
  f6g_write_tsv_atomic(inventory, output)
  message("Frozen ", nrow(inventory), " non-A/B Figure 6 files: ", output)
  invisible(output)
}

f6abg_verify <- function(workspace_root = f6abg_workspace_root(), baseline,
                         output) {
  f6abg_source_freeze_helpers(workspace_root)
  if (!file.exists(baseline)) stop("Missing Figure 6 non-A/B baseline: ", baseline)
  before <- utils::read.delim(baseline, check.names = FALSE)
  after <- f6abg_inventory(workspace_root)
  comparison <- merge(
    before, after, by = "relative_path", all = TRUE,
    suffixes = c("_before", "_after"), sort = TRUE
  )
  comparison$status <- ifelse(
    is.na(comparison$sha256_before), "added_to_non_ab_scope",
    ifelse(
      is.na(comparison$sha256_after), "missing",
      ifelse(
        comparison$sha256_before == comparison$sha256_after &
          comparison$size_bytes_before == comparison$size_bytes_after,
        "unchanged", "changed"
      )
    )
  )
  f6g_write_tsv_atomic(comparison, output)
  failures <- comparison$status != "unchanged"
  if (any(failures)) {
    stop(
      "Figure 6 non-A/B freeze guard failed for: ",
      paste(comparison$relative_path[failures], collapse = ", ")
    )
  }
  message(
    "Figure 6 non-A/B freeze guard passed for ", nrow(comparison), " files."
  )
  invisible(output)
}
