#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  stop("MININGCLONEID_REPO_ROOT is not set. Run with scripts/agentRrunner.sh.")
}
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
draft_rel <- "agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting"
draft_root <- file.path(repo_root, draft_rel)

fail <- function(...) stop(..., call. = FALSE)

artifacts <- data.frame(
  id = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Supplement"),
  stem = c(
    "final_figures/recommended/figure1",
    "final_figures/recommended/figure2",
    "final_figures/recommended/figure3",
    "final_figures/recommended/figure4",
    "final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions",
    "supplementary/fit_quality_optimizer_diagnostics"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(artifacts))) {
  png_path <- file.path(draft_root, paste0(artifacts$stem[[i]], ".png"))
  pdf_path <- file.path(draft_root, paste0(artifacts$stem[[i]], ".pdf"))
  if (!file.exists(png_path) || !file.exists(pdf_path)) {
    fail("Missing PNG/PDF pair for ", artifacts$id[[i]])
  }

  png_info <- attributes(png::readPNG(png_path, info = TRUE))$info
  if (is.null(png_info$dpi) || any(abs(png_info$dpi - 300) > 0.2)) {
    fail(artifacts$id[[i]], " is not a 300-dpi PNG.")
  }
  width_in <- png_info$dim[[1]] / png_info$dpi[[1]]
  height_in <- png_info$dim[[2]] / png_info$dpi[[2]]
  if (width_in > 7.1001 || height_in > 9.0001) {
    fail(
      artifacts$id[[i]], " exceeds the 7.1 x 9.0 inch style envelope: ",
      sprintf("%.3f x %.3f", width_in, height_in)
    )
  }

  pdf_metadata <- system2("pdfinfo", pdf_path, stdout = TRUE, stderr = TRUE)
  pdf_status <- attr(pdf_metadata, "status")
  if (!is.null(pdf_status) && pdf_status != 0L) {
    fail("pdfinfo failed for ", pdf_path)
  }
  pages <- sub("^Pages:[[:space:]]+", "", grep("^Pages:", pdf_metadata, value = TRUE))
  if (length(pages) != 1L || as.integer(pages) != 1L) {
    fail(artifacts$id[[i]], " PDF is not one page.")
  }
}

feedback <- readLines(file.path(draft_root, "feedback_intake.md"), warn = FALSE)
if (any(grepl("\\| pending \\|", feedback, fixed = FALSE))) {
  fail("feedback_intake.md still contains pending directive statuses.")
}

fidelity <- read.csv(
  file.path(draft_root, "prior_code_fidelity.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (any(!nzchar(fidelity$fidelity_status)) || any(fidelity$fidelity_status == "pending")) {
  fail("prior_code_fidelity.csv has blank or pending statuses.")
}
for (column in c("copied_baseline_code_path", "active_local_script_path", "diff_path")) {
  inherited <- fidelity$inheritance_class %in% c(
    "inherited_preserve", "inherited_targeted_fix", "inherited_move",
    "inherited_replace"
  )
  paths <- fidelity[[column]][inherited]
  if (any(!file.exists(file.path(repo_root, paths)))) {
    fail("Missing inherited-panel ", column, " artifact.")
  }
}

required_docs <- c(
  "feedback_manager_context.md",
  "feedback_intake.md",
  "prior_panel_disposition.csv",
  "prior_code_fidelity.csv",
  "source_manifest.csv",
  "drafting_panels.md",
  "not_drafted.md",
  "validation_report.md",
  "final_figures/recommended/legend.md",
  "report_manifest.csv",
  "review_report.html"
)
missing_docs <- required_docs[
  !file.exists(file.path(draft_root, required_docs))
]
if (length(missing_docs)) {
  fail("Missing required package document(s): ", paste(missing_docs, collapse = ", "))
}

report_manifest <- read.csv(
  file.path(draft_root, "report_manifest.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (any(report_manifest$include_status == "missing")) {
  fail("review_report.html omits one or more drafted PNGs.")
}
if (!all(report_manifest$include_status == "included")) {
  fail("Report manifest contains an unexpected non-included PNG status.")
}

source_manifest <- read.csv(
  file.path(draft_root, "source_manifest.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (anyDuplicated(source_manifest[c("figure", "role", "path")])) {
  fail("source_manifest.csv contains duplicate figure/role/path rows.")
}
real_sources <- !grepl("^git:", source_manifest$path)
missing_sources <- source_manifest$path[
  real_sources & !file.exists(file.path(repo_root, source_manifest$path))
]
if (length(missing_sources)) {
  fail("source_manifest.csv references missing files: ", paste(missing_sources, collapse = ", "))
}
sha256_file <- function(path) {
  result <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) fail("sha256sum failed for ", path)
  sub("[[:space:]].*$", "", result[[1]])
}
real_rows <- which(real_sources)
actual_sha256 <- vapply(
  file.path(repo_root, source_manifest$path[real_rows]),
  sha256_file,
  character(1)
)
hash_mismatch <- real_rows[actual_sha256 != source_manifest$sha256[real_rows]]
if (length(hash_mismatch)) {
  fail(
    "source_manifest.csv has stale SHA-256 values for: ",
    paste(source_manifest$path[hash_mismatch], collapse = ", ")
  )
}

if (any(grepl("figure3e|figure6", list.files(
  file.path(draft_root, "final_figures"),
  recursive = TRUE,
  ignore.case = TRUE
)))) {
  fail("An excluded Figure 3E or Figure 6 artifact exists in final_figures/.")
}

cat(
  "Drafting package validation: PASS\n",
  "- ", nrow(artifacts), " PNG/PDF pairs within 7.1 x 9.0 inches\n",
  "- ", nrow(report_manifest), " drafted PNGs represented in review_report.html\n",
  "- ", nrow(source_manifest), " provenance rows with existing real-file sources\n",
  sep = ""
)
