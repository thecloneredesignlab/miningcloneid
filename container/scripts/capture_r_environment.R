#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: capture_r_environment.R OUT_DIR DIRECT_PACKAGES_TSV USAGE_TSV", call. = FALSE)
}

out_dir <- args[[1L]]
direct_path <- args[[2L]]
usage_path <- args[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, name) {
  x[] <- lapply(x, function(column) {
    if (!is.character(column)) return(column)
    gsub("[\t\r\n]+", " ", column, perl = TRUE)
  })
  utils::write.table(
    x,
    file = file.path(out_dir, name),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

scalar_or_na <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  as.character(x[[1L]])
}

capture_output_file <- function(expr, name) {
  writeLines(utils::capture.output(force(expr)), file.path(out_dir, name))
}

capture_output_file(sessionInfo(), "sessionInfo.txt")
capture_output_file(capabilities(), "capabilities.txt")
capture_output_file(extSoftVersion(), "extSoftVersion.txt")
writeLines(R.version.string, file.path(out_dir, "R-version.txt"))

write_tsv(
  data.frame(
    order = seq_along(.libPaths()),
    library_path = .libPaths(),
    exists = file.exists(.libPaths()),
    stringsAsFactors = FALSE
  ),
  "r-libpaths.tsv"
)

r_detail <- data.frame(
  field = c(
    "version.string", "major", "minor", "platform", "arch", "os",
    "system", "status", "svn.rev", "language", "nickname", "R.home"
  ),
  value = c(
    scalar_or_na(R.version.string), scalar_or_na(R.version$major),
    scalar_or_na(R.version$minor), scalar_or_na(R.version$platform),
    scalar_or_na(R.version$arch), scalar_or_na(R.version$os),
    scalar_or_na(R.version$system), scalar_or_na(R.version$status),
    scalar_or_na(R.version$svn.rev), scalar_or_na(R.version$language),
    scalar_or_na(R.version$nickname), scalar_or_na(R.home())
  ),
  stringsAsFactors = FALSE
)
write_tsv(r_detail, "r-version-details.tsv")

extra_fields <- c(
  "Depends", "Imports", "LinkingTo", "Suggests", "Enhances", "Repository",
  "RemoteType", "RemoteHost", "RemoteRepo", "RemoteUsername", "RemoteRef",
  "RemoteSha", "GithubRepo", "GithubUsername", "GithubRef", "GithubSHA1",
  "Built", "NeedsCompilation", "SystemRequirements", "Packaged"
)

manifest_rows <- list()
for (lib_order in seq_along(.libPaths())) {
  lib <- .libPaths()[[lib_order]]
  if (!dir.exists(lib)) next
  ip <- tryCatch(
    installed.packages(lib.loc = lib, fields = extra_fields, noCache = TRUE),
    error = function(e) NULL
  )
  if (is.null(ip) || !nrow(ip)) next
  frame <- as.data.frame(ip, stringsAsFactors = FALSE, check.names = FALSE)
  frame$library_order <- lib_order
  frame$library_path <- lib
  frame$package_path <- file.path(lib, frame$Package)
  manifest_rows[[length(manifest_rows) + 1L]] <- frame
}

manifest <- if (length(manifest_rows)) do.call(rbind, manifest_rows) else data.frame()
if (!nrow(manifest)) stop("No installed R packages were discovered", call. = FALSE)
manifest <- manifest[order(manifest$library_order, manifest$Package), , drop = FALSE]
write_tsv(manifest, "r-installed-packages-all.tsv")

selected <- manifest[!duplicated(manifest$Package), , drop = FALSE]
selected$selected_by_libpaths <- TRUE
write_tsv(selected, "r-installed-packages-selected.tsv")

direct <- utils::read.delim(direct_path, check.names = FALSE, stringsAsFactors = FALSE)
usage <- utils::read.delim(usage_path, check.names = FALSE, stringsAsFactors = FALSE)
runtime_scopes <- c("analysis_direct", "o2_runtime", "o2_auxiliary")
direct_packages <- sort(unique(direct$package[direct$dependency_scope %in% runtime_scopes]))
direct_packages <- direct_packages[nzchar(direct_packages)]

resolution <- lapply(direct_packages, function(package) {
  hits <- manifest[manifest$Package == package, , drop = FALSE]
  selected_hit <- selected[selected$Package == package, , drop = FALSE]
  data.frame(
    package = package,
    installed = nrow(selected_hit) == 1L,
    selected_version = if (nrow(selected_hit)) selected_hit$Version[[1L]] else NA_character_,
    selected_library = if (nrow(selected_hit)) selected_hit$library_path[[1L]] else NA_character_,
    all_installed_versions = if (nrow(hits)) paste(unique(hits$Version), collapse = ";") else NA_character_,
    all_installed_libraries = if (nrow(hits)) paste(unique(hits$library_path), collapse = ";") else NA_character_,
    requireNamespace_ok = isTRUE(tryCatch(requireNamespace(package, quietly = TRUE), error = function(e) FALSE)),
    stringsAsFactors = FALSE
  )
})
resolution <- if (length(resolution)) do.call(rbind, resolution) else data.frame()
write_tsv(resolution, "r-package-resolution.tsv")

split_dependencies <- function(value) {
  if (is.na(value) || !nzchar(value)) return(character())
  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  pieces <- trimws(sub("[[:space:]]*\\(.*\\)[[:space:]]*$", "", pieces))
  setdiff(pieces[nzchar(pieces)], "R")
}

dep_fields <- c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")
edge_rows <- list()
for (i in seq_len(nrow(selected))) {
  parent <- selected$Package[[i]]
  for (field in dep_fields) {
    deps <- split_dependencies(selected[[field]][[i]])
    if (!length(deps)) next
    edge_rows[[length(edge_rows) + 1L]] <- data.frame(
      parent = parent,
      dependency = deps,
      dependency_type = field,
      stringsAsFactors = FALSE
    )
  }
}
edges <- if (length(edge_rows)) unique(do.call(rbind, edge_rows)) else data.frame(
  parent = character(), dependency = character(), dependency_type = character()
)
write_tsv(edges, "r-installed-package-dependency-edges.tsv")

closure <- function(roots, allowed_types, label) {
  queue <- data.frame(
    root = roots,
    package = roots,
    depth = 0L,
    via_parent = NA_character_,
    via_type = "direct",
    stringsAsFactors = FALSE
  )
  out <- queue
  cursor <- 1L
  seen <- paste(queue$root, queue$package, sep = "\r")
  while (cursor <= nrow(queue)) {
    row <- queue[cursor, , drop = FALSE]
    cursor <- cursor + 1L
    children <- edges[
      edges$parent == row$package[[1L]] & edges$dependency_type %in% allowed_types,
      ,
      drop = FALSE
    ]
    if (!nrow(children)) next
    additions <- data.frame(
      root = row$root[[1L]],
      package = children$dependency,
      depth = row$depth[[1L]] + 1L,
      via_parent = row$package[[1L]],
      via_type = children$dependency_type,
      stringsAsFactors = FALSE
    )
    keys <- paste(additions$root, additions$package, sep = "\r")
    additions <- additions[!keys %in% seen, , drop = FALSE]
    if (nrow(additions)) {
      seen <- c(seen, paste(additions$root, additions$package, sep = "\r"))
      queue <- rbind(queue, additions)
      out <- rbind(out, additions)
    }
  }
  out$closure_scope <- label
  hit <- match(out$package, selected$Package)
  out$version <- selected$Version[hit]
  out$library_path <- selected$library_path[hit]
  out$installed <- !is.na(hit)
  out
}

runtime_closure <- closure(direct_packages, c("Depends", "Imports", "LinkingTo"), "runtime")
full_closure <- closure(direct_packages, dep_fields, "full_with_optional")
write_tsv(runtime_closure, "r-runtime-dependency-closure.tsv")
write_tsv(full_closure, "r-full-dependency-closure.tsv")

analysis_roots <- sort(unique(usage$package[usage$dependency_scope == "analysis_direct"]))
analysis_runtime <- closure(analysis_roots, c("Depends", "Imports", "LinkingTo"), "analysis_runtime")
analysis_full <- closure(analysis_roots, dep_fields, "analysis_full_with_optional")
write_tsv(analysis_runtime, "r-analysis-runtime-dependency-closure.tsv")
write_tsv(analysis_full, "r-analysis-full-dependency-closure.tsv")

description_hash <- lapply(seq_len(nrow(selected)), function(i) {
  path <- file.path(selected$package_path[[i]], "DESCRIPTION")
  data.frame(
    Package = selected$Package[[i]],
    Version = selected$Version[[i]],
    description_path = path,
    description_md5 = if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_,
    stringsAsFactors = FALSE
  )
})
write_tsv(do.call(rbind, description_hash), "r-package-description-hashes.tsv")
