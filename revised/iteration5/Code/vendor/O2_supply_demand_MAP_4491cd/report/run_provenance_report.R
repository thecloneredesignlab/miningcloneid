`%||%` <- if (exists("%||%", inherits = TRUE)) get("%||%", inherits = TRUE) else function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

o2sd_prov_escape_html <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

o2sd_prov_table_html <- function(tab, max_rows = 200L) {
  if (!is.data.frame(tab) || nrow(tab) == 0L) {
    return('<p class="report-empty">No run provenance records were found.</p>')
  }
  tab <- tab[seq_len(min(nrow(tab), as.integer(max_rows))), , drop = FALSE]
  tab[] <- lapply(tab, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  header <- paste(sprintf("<th>%s</th>", o2sd_prov_escape_html(names(tab))), collapse = "")
  rows <- vapply(seq_len(nrow(tab)), function(i) {
    vals <- vapply(tab[i, , drop = FALSE], function(x) o2sd_prov_escape_html(x[[1]]), character(1))
    paste0("<tr>", paste(sprintf("<td>%s</td>", vals), collapse = ""), "</tr>")
  }, character(1))
  paste0(
    '<table class="report-table"><thead><tr>', header, '</tr></thead><tbody>',
    paste(rows, collapse = ""),
    "</tbody></table>"
  )
}

o2sd_prov_read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.delim(
      path,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
}

o2sd_prov_read_text <- function(path, max_lines = 200L) {
  if (!file.exists(path)) return(NULL)
  txt <- tryCatch(readLines(path, warn = FALSE, n = max_lines), error = function(e) NULL)
  if (is.null(txt)) return(NULL)
  paste(txt, collapse = "\n")
}

o2sd_prov_is_seed_dir <- function(path) {
  grepl("^seed[0-9]+$", basename(normalizePath(path, mustWork = FALSE)), ignore.case = TRUE)
}

o2sd_prov_parent_run_dir <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (basename(path) == "extra_results") return(dirname(path))
  if (o2sd_prov_is_seed_dir(path)) return(dirname(path))
  path
}

o2sd_prov_find_first_seed_file <- function(run_dir, filename) {
  seed_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs), ignore.case = TRUE)]
  if (!length(seed_dirs)) return(NA_character_)
  paths <- file.path(seed_dirs, filename)
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1]] else NA_character_
}

o2sd_prov_flatten_yaml <- function(x, prefix = character(0)) {
  if (is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    out <- Map(function(item, nm) o2sd_prov_flatten_yaml(item, c(prefix, nm)), x, nms)
    return(do.call(rbind, out))
  }
  key <- paste(prefix, collapse = ".")
  value <- if (length(x) == 0L) {
    ""
  } else if (length(x) == 1L) {
    as.character(x)
  } else {
    paste(as.character(x), collapse = ",")
  }
  data.frame(key = key, value = value, stringsAsFactors = FALSE)
}

o2sd_prov_read_yaml_table <- function(path, source_label) {
  if (!file.exists(path) || !requireNamespace("yaml", quietly = TRUE)) return(NULL)
  y <- tryCatch(yaml::read_yaml(path), error = function(e) NULL)
  if (is.null(y)) return(NULL)
  tab <- o2sd_prov_flatten_yaml(y)
  if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
  tab <- tab[nzchar(tab$key), , drop = FALSE]
  data.frame(source = source_label, tab, stringsAsFactors = FALSE)
}

o2sd_prov_command_table <- function(run_dir) {
  parent_dir <- o2sd_prov_parent_run_dir(run_dir)
  command_files <- c(
    run_command = file.path(parent_dir, "run_command.txt"),
    sbatch_command = file.path(parent_dir, "sbatch_command.txt"),
    fit_command = file.path(run_dir, "fit_command.txt")
  )
  if (!file.exists(command_files[["fit_command"]])) {
    first_fit_command <- o2sd_prov_find_first_seed_file(parent_dir, "fit_command.txt")
    if (!is.na(first_fit_command)) command_files[["fit_command_sample_seed"]] <- first_fit_command
  }
  rows <- lapply(names(command_files), function(nm) {
    txt <- o2sd_prov_read_text(command_files[[nm]])
    if (is.null(txt)) return(NULL)
    data.frame(command_type = nm, path = command_files[[nm]], command = txt, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

o2sd_prov_metadata_table <- function(run_dir) {
  parent_dir <- o2sd_prov_parent_run_dir(run_dir)
  candidate_dirs <- unique(c(run_dir, parent_dir))
  prov_rows <- lapply(candidate_dirs, function(path) {
    tab <- o2sd_prov_read_tsv(file.path(path, "run_provenance.tsv"))
    if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
    tab$provenance_dir <- path
    tab
  })
  prov_rows <- Filter(Negate(is.null), prov_rows)
  if (length(prov_rows)) {
    return(do.call(rbind, prov_rows))
  }
  manifest <- o2sd_prov_read_tsv(file.path(parent_dir, "fit_array_manifest.tsv"))
  if (is.data.frame(manifest) && all(c("key", "value") %in% names(manifest))) {
    return(data.frame(
      section = "legacy_fit_array_manifest",
      key = manifest$key,
      value = manifest$value,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    section = "provenance",
    key = "status",
    value = "No run_provenance.tsv or fit_array_manifest.tsv found.",
    stringsAsFactors = FALSE
  )
}

o2sd_prov_effective_args_table <- function(run_dir) {
  parent_dir <- o2sd_prov_parent_run_dir(run_dir)
  paths <- c(
    file.path(run_dir, "run_effective_args.tsv"),
    file.path(parent_dir, "run_effective_args.tsv")
  )
  first_seed <- o2sd_prov_find_first_seed_file(parent_dir, "run_effective_args.tsv")
  if (!is.na(first_seed)) paths <- c(paths, first_seed)
  paths <- unique(paths[file.exists(paths)])
  if (!length(paths)) return(NULL)
  rows <- lapply(paths, function(path) {
    tab <- o2sd_prov_read_tsv(path)
    if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
    tab$path <- path
    tab
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

o2sd_prov_yaml_tables <- function(run_dir) {
  parent_dir <- o2sd_prov_parent_run_dir(run_dir)
  paths <- c(
    config_input = file.path(parent_dir, "config.input.yaml"),
    config_resolved = file.path(parent_dir, "config.resolved.yaml")
  )
  rows <- lapply(names(paths), function(nm) o2sd_prov_read_yaml_table(paths[[nm]], nm))
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

o2sd_run_provenance_html <- function(run_dir, title = "Run Provenance", section_id = "run-provenance") {
  run_dir <- normalizePath(run_dir, mustWork = FALSE)
  metadata <- o2sd_prov_metadata_table(run_dir)
  commands <- o2sd_prov_command_table(run_dir)
  args <- o2sd_prov_effective_args_table(run_dir)
  yaml <- o2sd_prov_yaml_tables(run_dir)

  paste0(
    '<section class="report-section" id="', o2sd_prov_escape_html(section_id), '">',
    '<h2 class="report-figure-title">', o2sd_prov_escape_html(title), '</h2>',
    '<p class="report-figure-legend">Execution metadata, input configuration snapshots, and the effective runtime arguments recorded for this result.</p>',
    '<h3>Execution Metadata</h3>',
    o2sd_prov_table_html(metadata, max_rows = 300L),
    '<h3>Commands</h3>',
    o2sd_prov_table_html(commands, max_rows = 20L),
    '<h3>Effective Runtime Arguments</h3>',
    o2sd_prov_table_html(args, max_rows = 400L),
    '<h3>YAML Configuration Snapshot</h3>',
    o2sd_prov_table_html(yaml, max_rows = 400L),
    '</section>'
  )
}
