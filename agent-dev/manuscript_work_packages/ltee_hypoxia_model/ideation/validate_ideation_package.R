#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
package_root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation"
}

required <- file.path(
  package_root,
  c(
    "ideation_report.html",
    "ideas.md",
    "feedback_manager_context.md",
    "existing_panel_disposition.csv",
    "tex_figure_requirements.csv",
    "source_inventory.csv",
    "feedback_decisions.csv",
    "fit_quality_evidence.csv",
    "existing_context_manifest.tsv",
    "baseline_build.md",
    "build_ideation_report.R"
  )
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing required artifact(s): ", paste(missing, collapse = ", "))
}

disposition <- read.csv(
  file.path(package_root, "existing_panel_disposition.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
requirements <- read.csv(
  file.path(package_root, "tex_figure_requirements.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
inventory <- read.csv(
  file.path(package_root, "source_inventory.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
decisions <- read.csv(
  file.path(package_root, "feedback_decisions.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fit_quality <- read.csv(
  file.path(package_root, "fit_quality_evidence.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manifest <- read.delim(
  file.path(package_root, "existing_context_manifest.tsv"),
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expected_disposition_columns <- c(
  "artifact_id", "path", "current_role", "disposition", "target_role",
  "rationale", "required_action"
)
stopifnot(identical(names(disposition), expected_disposition_columns))
stopifnot(nrow(disposition) == 16L)
stopifnot(!anyDuplicated(disposition$artifact_id))
allowed_dispositions <- c(
  "preserve", "targeted_fix", "move_or_duplicate", "replace", "drop",
  "uncertain", "not_applicable", "defer"
)
unexpected <- setdiff(unique(disposition$disposition), allowed_dispositions)
if (length(unexpected) > 0L) {
  stop("Unexpected disposition value(s): ", paste(unexpected, collapse = ", "))
}

stopifnot(nrow(requirements) == 28L)
stopifnot(!anyDuplicated(requirements$panel_id))
stopifnot(nrow(inventory) == 24L)
stopifnot(!anyDuplicated(inventory$source_id))
stopifnot(nrow(decisions) == 18L)
stopifnot(!anyDuplicated(decisions$decision_id))
stopifnot(nrow(fit_quality) == 3L)
stopifnot(!anyDuplicated(fit_quality$figure))
stopifnot(nrow(manifest) == 6L)
stopifnot(sum(grepl("\\.png$", manifest$local_review_copy)) == 5L)

html_path <- file.path(package_root, "ideation_report.html")
html <- paste(readLines(html_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
stopifnot(grepl("<!doctype html>", html, fixed = TRUE))
stopifnot(grepl("<title>Approved LTEE hypoxia manuscript figure plan</title>", html, fixed = TRUE))
stopifnot(lengths(regmatches(html, gregexpr("data:image/png;base64,", html, fixed = TRUE))) == 5L)
stopifnot(!grepl("<script", html, ignore.case = TRUE))
stopifnot(!grepl("@import", html, ignore.case = TRUE))
stopifnot(!grepl("url\\(['\"]?https?://", html, ignore.case = TRUE))

src_matches <- regmatches(
  html,
  gregexpr("(?:src|href)=['\"][^'\"]+['\"]", html, perl = TRUE)
)[[1L]]
src_values <- sub("^(?:src|href)=['\"]", "", src_matches, perl = TRUE)
src_values <- sub("['\"]$", "", src_values)
allowed_reference <- grepl("^(data:|#)", src_values)
if (!all(allowed_reference)) {
  stop(
    "Non-self-contained src/href value(s): ",
    paste(src_values[!allowed_reference], collapse = ", ")
  )
}

required_sections <- c(
  "1. Existing Figures And Panels",
  "2. Relevant Feedback And Recorded Decisions",
  "3. Existing-Panel Disposition",
  "4. Revised Working Figure Architecture",
  "5. Gate A Resolution",
  "6. Appendix: Sources, Caveats, And Deferred Work"
)
section_positions <- vapply(
  required_sections,
  function(section) regexpr(section, html, fixed = TRUE)[[1L]],
  integer(1L)
)
stopifnot(all(section_positions > 0L))
stopifnot(all(diff(section_positions) > 0L))

feedback_context <- paste(
  readLines(
    file.path(package_root, "feedback_manager_context.md"),
    warn = FALSE,
    encoding = "UTF-8"
  ),
  collapse = "\n"
)
stopifnot(grepl("SPAN-20260724-001", feedback_context, fixed = TRUE))
stopifnot(grepl("SPAN-20260703-001", feedback_context, fixed = TRUE))
stopifnot(grepl("ITEM-20260723T122502-001-manuscript-figure-workflow", feedback_context, fixed = TRUE))
stopifnot(grepl("Gate A is closed", feedback_context, fixed = TRUE))

required_report_claims <- c(
  "5</strong><span>active main figures",
  "0</strong><span>refits authorized",
  "Figure 3E is rejected",
  "unchanged pooled in-vivo/in-vitro embedding",
  "must not be filtered or recomputed as an in-vivo-only embedding",
  "fit adequacy",
  "not confidence intervals",
  "retain all three already displayed in-vivo regions",
  "Drafting authorized within bounds",
  "no drafting decision remains open"
)
for (claim in required_report_claims) {
  if (!grepl(claim, html, fixed = TRUE)) {
    stop("Required revised-report statement missing: ", claim)
  }
}

cat("PASS required ideation artifacts:", length(required), "\n")
cat("PASS disposition rows:", nrow(disposition), "\n")
cat("PASS TeX panel roles:", nrow(requirements), "\n")
cat("PASS source records:", nrow(inventory), "\n")
cat("PASS recorded feedback decisions:", nrow(decisions), "\n")
cat("PASS fit-quality evidence rows:", nrow(fit_quality), "\n")
cat("PASS embedded historical PNGs: 5\n")
cat("PASS report section order and self-contained references\n")
cat("PASS current and earlier feedback-manager context\n")
