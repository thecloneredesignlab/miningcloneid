#!/usr/bin/env Rscript
# Recoverable publication-only archive. Immutable calculation runs are excluded.
root <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT", getwd()), mustWork = TRUE)
stopifnot(basename(root) == "iteration4", dir.exists(file.path(root, "Code", "Figures")))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, grepl("^[A-Za-z0-9_.-]+$", args))
archive <- file.path(root, "audit", "figure7_publication_archive", args)
directories <- c(file.path(root, "Figures"), file.path(root, "manuscript", "Figures"),
  file.path(root, "data", "Figures", paste0("Supp_Figure7_", 8:12)))
names <- c("assembled_fig7", "supp_fig7-8_inverse_response",
  "supp_fig7-9_invivo_continuous_full_range", "supp_fig7-10_invitro_continuous_full_range",
  "supp_fig7-11_invitro_passage_full_range")
records <- list()
for (directory in directories) {
  files <- list.files(directory, full.names = TRUE)
  for (path in files) {
    if (!any(vapply(names, function(name) startsWith(basename(path), paste0(name, ".")) ||
      startsWith(basename(path), paste0(name, "_render_validation.")), logical(1)))) next
    if (dir.exists(path)) next
    relative <- substring(path, nchar(root) + 2L)
    target <- file.path(archive, relative)
    if (file.exists(target)) stop("Archive already contains ", target)
    hash <- unname(tools::md5sum(path))
    size <- file.info(path)$size
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(path, target)) stop("Cannot archive ", path)
    stopifnot(identical(unname(tools::md5sum(target)), hash))
    records[[length(records) + 1L]] <- data.frame(path = relative, size = size,
      md5 = hash, archived_path = substring(target, nchar(root) + 2L))
  }
}
dir.create(archive, recursive = TRUE, showWarnings = FALSE)
if (length(records)) write.table(do.call(rbind, records), file.path(archive, "manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
cat("Archived", length(records), "previous publication files to", archive, "\n")
