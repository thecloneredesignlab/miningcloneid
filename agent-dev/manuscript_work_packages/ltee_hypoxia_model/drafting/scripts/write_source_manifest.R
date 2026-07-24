#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  stop("MININGCLONEID_REPO_ROOT is not set. Run with scripts/agentRrunner.sh.")
}
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
draft_rel <- "agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting"
draft_root <- file.path(repo_root, draft_rel)

sha256_file <- function(path) {
  result <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, ": ", paste(result, collapse = "\n"))
  }
  sub("[[:space:]].*$", "", result[[1]])
}

records <- list()
record_index <- 0L

add_record <- function(figure, role, relative_path, provenance_note = "") {
  absolute_path <- file.path(repo_root, relative_path)
  if (!file.exists(absolute_path)) {
    stop("Manifest source does not exist: ", relative_path)
  }
  info <- file.info(absolute_path)
  record_index <<- record_index + 1L
  records[[record_index]] <<- data.frame(
    figure = figure,
    role = role,
    path = gsub("\\\\", "/", relative_path),
    sha256 = sha256_file(absolute_path),
    bytes = as.numeric(info$size),
    provenance_note = provenance_note,
    stringsAsFactors = FALSE
  )
}

add_many <- function(figure, role, paths, provenance_note = "") {
  for (path in paths) add_record(figure, role, path, provenance_note)
}

# Gate, manuscript, and style references are frozen under source_context/ and
# are recorded by the local-file walk below. This avoids making manifest
# regeneration depend on ignored or environment-owned files.

# Figure 2 is a semantic diagram of the implemented model.
add_many(
  "Figure 2", "semantic_code_source",
  c(
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R",
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp"
  )
)

# Figures 1 and 3--5, plus the optimizer supplement, read only the immutable
# copies under source_tables/frozen_inputs/.  The local-file walk below records
# every frozen file and the frozen_input_manifest.csv provenance map, so this
# manifest can be regenerated without access to the upstream ignored bundle.

# Add every local drafting artifact except self-referential report products.
local_files <- list.files(
  draft_root,
  recursive = TRUE,
  full.names = FALSE,
  all.files = FALSE
)
local_files <- local_files[
  !local_files %in% c(
    "source_manifest.csv",
    "report_manifest.csv",
    "review_report.html"
  )
]

infer_figure <- function(path) {
  if (grepl("(^|/)F1|figure1", path, ignore.case = TRUE)) return("Figure 1")
  if (grepl("(^|/)F2|figure2", path, ignore.case = TRUE)) return("Figure 2")
  if (grepl("(^|/)F3|figure3", path, ignore.case = TRUE)) return("Figure 3")
  if (grepl("(^|/)F4|figure4", path, ignore.case = TRUE)) return("Figure 4")
  if (grepl("(^|/)F5|figure5", path, ignore.case = TRUE)) return("Figure 5")
  if (grepl("diagnostic|supplementary", path, ignore.case = TRUE)) {
    return("Supplementary diagnostics")
  }
  "package"
}

infer_role <- function(path) {
  if (grepl("^scripts/.*\\.R$", path)) return("local_generator")
  if (grepl("^prior_code/", path)) return("copied_prior_code")
  if (grepl("^code_diffs/", path)) return("fidelity_diff")
  if (grepl("^source_context/", path)) return("frozen_scope_reference")
  if (identical(path, "source_tables/frozen_input_manifest.csv")) {
    return("frozen_input_provenance")
  }
  if (grepl("^source_tables/frozen_inputs/", path)) {
    return("frozen_generator_input")
  }
  if (grepl("\\.(png|pdf)$", path, ignore.case = TRUE)) return("rendered_draft")
  if (grepl("^source_tables/|\\.tsv$", path, ignore.case = TRUE)) return("derived_table")
  if (grepl("\\.csv$", path, ignore.case = TRUE)) return("ledger")
  if (grepl("\\.md$", path, ignore.case = TRUE)) return("documentation")
  "package_artifact"
}

for (path in local_files) {
  add_record(
    infer_figure(path),
    infer_role(path),
    file.path(draft_rel, path)
  )
}

# The pooled embedding raster is copied byte-for-byte from this pinned git blob.
pooled_copy <- file.path(
  draft_rel,
  "initial_subpanels/F4D/figure4D_pooled_embedding_preserved_source.png"
)
pooled_row <- records[[
  which(vapply(
    records,
    function(x) identical(x$path[[1]], pooled_copy),
    logical(1)
  ))
]]
record_index <- record_index + 1L
records[[record_index]] <- data.frame(
  figure = "Figure 4",
  role = "pinned_git_blob",
  path = paste0(
    "git:7e72dca88caf784fc61271d87a1c0dfb564b8303:",
    "oxygen/figures/iteration1/fig4f_landscape_tsne_clusters.png"
  ),
  sha256 = pooled_row$sha256,
  bytes = pooled_row$bytes,
  provenance_note = "Byte-identical to the preserved local copy; generator also asserts MD5 14cecff29ab4884823b84d83f0379119.",
  stringsAsFactors = FALSE
)

manifest <- do.call(rbind, records)
manifest <- manifest[!duplicated(manifest[c("figure", "role", "path")]), , drop = FALSE]
manifest <- manifest[order(manifest$figure, manifest$role, manifest$path), , drop = FALSE]
write.csv(
  manifest,
  file.path(draft_root, "source_manifest.csv"),
  row.names = FALSE,
  quote = TRUE
)
message("Wrote ", nrow(manifest), " source-manifest rows.")
