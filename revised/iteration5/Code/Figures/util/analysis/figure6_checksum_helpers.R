#!/usr/bin/env Rscript

# Shared path, hashing, inventory, and atomic-TSV helpers used by the current
# Figure 6 workflow.  The obsolete pre-promotion whole-figure freeze command
# has been removed.

options(stringsAsFactors = FALSE, warn = 1)

f6g_workspace_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "Code", "Figures")) &&
        dir.exists(file.path(current, "data", "Figures"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Cannot locate iteration4 workspace.")
    current <- parent
  }
}

f6g_repository_root <- function(workspace_root) {
  normalizePath(file.path(workspace_root, "..", ".."), mustWork = TRUE)
}

f6g_relative <- function(path, root) {
  path <- normalizePath(path, mustWork = TRUE)
  root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(path, prefix)) stop("Path outside repository: ", path)
  substring(path, nchar(prefix) + 1L)
}

f6g_figure6_files <- function(workspace_root) {
  roots <- file.path(
    workspace_root,
    c("Code/Figures", "Code/hpc", "data/Figures", "Figures", "manuscript/Figures")
  )
  roots <- roots[dir.exists(roots)]
  files <- unique(unlist(lapply(roots, function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
               no.. = TRUE, include.dirs = FALSE)
  }), use.names = FALSE))
  relative_to_workspace <- substring(
    normalizePath(files, mustWork = TRUE),
    nchar(normalizePath(workspace_root, mustWork = TRUE)) + 2L
  )
  keep <- grepl(
    "(^|/)([^/]*(Figure6|figure6|Fig6|fig6)[^/]*)(/|$)|(^|/)Supp_Figure6_[^/]+(/|$)|(^|/)supp_fig6-[^/]+$",
    relative_to_workspace,
    perl = TRUE
  )
  sort(files[keep])
}

f6g_sha256 <- function(files) {
  executable <- Sys.which("shasum")
  flavor <- "shasum"
  if (!nzchar(executable)) {
    executable <- Sys.which("sha256sum")
    flavor <- "sha256sum"
  }
  if (!nzchar(executable)) {
    stop("shasum or sha256sum is required for Figure 6 freeze verification.")
  }
  # Batch files to avoid thousands of process launches while staying below the
  # operating system command-line limit.
  groups <- split(files, ceiling(seq_along(files) / 100L))
  hashes <- unlist(lapply(groups, function(group) {
    arguments <- if (identical(flavor, "shasum")) {
      c("-a", "256", group)
    } else {
      group
    }
    output <- system2(executable, arguments, stdout = TRUE)
    if (length(output) != length(group)) {
      stop("Failed to hash one or more Figure 6 files in a batch.")
    }
    sub("[[:space:]].*$", "", output)
  }), use.names = FALSE)
  if (length(hashes) != length(files)) stop("Incomplete SHA256 inventory.")
  hashes
}

f6g_git_status <- function(repository_root, relative_paths) {
  scope <- sub("^(([^/]+/){2}).*$", "\\1", relative_paths[[1L]])
  command <- c(
    "-C", repository_root, "status", "--porcelain=v1", "--untracked-files=all",
    "--", scope
  )
  output <- suppressWarnings(system2(
    "git", command, stdout = TRUE, stderr = TRUE
  ))
  status <- rep("clean", length(relative_paths))
  command_status <- attr(output, "status")
  if (!is.null(command_status) && command_status != 0L) {
    return(rep("unavailable", length(relative_paths)))
  }
  if (!length(output)) return(status)
  for (line in output) {
    if (nchar(line) < 4L) next
    raw_path <- substring(line, 4L)
    if (grepl(" -> ", raw_path, fixed = TRUE)) {
      raw_path <- sub("^.* -> ", "", raw_path)
    }
    index <- match(raw_path, relative_paths)
    if (!is.na(index)) status[[index]] <- substring(line, 1L, 2L)
  }
  status
}

f6g_inventory <- function(workspace_root = f6g_workspace_root()) {
  repository_root <- f6g_repository_root(workspace_root)
  files <- f6g_figure6_files(workspace_root)
  if (!length(files)) stop("No Figure 6 files found for freeze inventory.")
  relative_paths <- vapply(files, f6g_relative, character(1L), root = repository_root)
  info <- file.info(files)
  data.frame(
    relative_path = relative_paths,
    size_bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    sha256 = f6g_sha256(files),
    git_status_at_freeze = f6g_git_status(repository_root, relative_paths),
    stringsAsFactors = FALSE
  )
}

f6g_write_tsv_atomic <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(data, temporary, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
  if (!file.rename(temporary, path)) stop("Cannot publish freeze manifest: ", path)
  normalizePath(path, mustWork = TRUE)
}
