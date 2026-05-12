#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0) {
      out[[kv]] <- TRUE
    } else {
      key <- substr(kv, 1, pos - 1)
      val <- substr(kv, pos + 1, nchar(kv))
      out[[key]] <- val
    }
  }
  out
}

first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

as_chr <- function(x, default = "") {
  val <- as.character(first_non_null(x, default))
  if (!length(val) || !nzchar(val[[1]])) default else val[[1]]
}

split_csv <- function(x, default = character()) {
  txt <- trimws(as_chr(x, paste(default, collapse = ",")))
  if (!nzchar(txt)) return(default)
  vals <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

seed_number <- function(seed_label) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_label))))
}

seed_dir_for <- function(run_dir, seed_label) {
  label <- as.character(seed_label)
  if (grepl("^seed[0-9]+$", label)) {
    return(file.path(run_dir, label))
  }
  if (grepl("^[0-9]+$", label)) {
    return(file.path(run_dir, paste0("seed", label)))
  }
  file.path(run_dir, label)
}

choose_objective <- function(tab, objective_columns) {
  objective <- rep(NA_real_, nrow(tab))
  objective_source <- rep(NA_character_, nrow(tab))
  for (col in objective_columns) {
    if (!col %in% names(tab)) next
    vals <- suppressWarnings(as.numeric(tab[[col]]))
    fill <- !is.finite(objective) & is.finite(vals)
    objective[fill] <- vals[fill]
    objective_source[fill] <- col
  }
  list(value = objective, source = objective_source)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  run_dir <- normalizePath(as_chr(argv$run_dir), mustWork = FALSE)
  if (!dir.exists(run_dir)) {
    stop("run_dir does not exist: ", run_dir, call. = FALSE)
  }

  summary_tsv <- normalizePath(
    as_chr(argv$summary_tsv, file.path(run_dir, "extra_results", "seed_summary.tsv")),
    mustWork = FALSE
  )
  if (!file.exists(summary_tsv)) {
    stop("summary_tsv does not exist: ", summary_tsv, call. = FALSE)
  }

  out_tsv <- normalizePath(
    as_chr(argv$out_tsv, file.path(run_dir, "best_seed_from_summary.tsv")),
    mustWork = FALSE
  )
  best_dir_file <- normalizePath(
    as_chr(argv$best_dir_file, file.path(run_dir, "best_seed_from_summary.dir")),
    mustWork = FALSE
  )
  objective_columns <- split_csv(argv$objective_columns, c("objective", "objective_total"))
  required_files <- split_csv(argv$required_files, c("best_params_transformed.tsv"))

  tab <- utils::read.delim(summary_tsv, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab)) {
    stop("summary_tsv has no rows: ", summary_tsv, call. = FALSE)
  }
  if (!"seed" %in% names(tab)) {
    stop("summary_tsv is missing required column 'seed': ", summary_tsv, call. = FALSE)
  }

  obj <- choose_objective(tab, objective_columns)
  seed_dirs <- vapply(tab$seed, seed_dir_for, character(1), run_dir = run_dir)
  seed_nums <- vapply(tab$seed, seed_number, integer(1))
  missing_files <- vapply(
    seed_dirs,
    function(seed_dir) {
      missing <- required_files[!file.exists(file.path(seed_dir, required_files))]
      paste(missing, collapse = ",")
    },
    character(1)
  )

  res <- tab
  res$seed_dir <- normalizePath(seed_dirs, mustWork = FALSE)
  res$selection_objective <- obj$value
  res$selection_objective_source <- obj$source
  res$selection_missing_required_files <- missing_files
  res$eligible <- is.finite(obj$value) & dir.exists(seed_dirs) & !nzchar(missing_files)
  res$selected <- FALSE

  eligible <- res[res$eligible %in% TRUE, , drop = FALSE]
  if (!nrow(eligible)) {
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    stop(
      "No eligible seed found in ", summary_tsv,
      ". Required objective columns: ", paste(objective_columns, collapse = ","),
      ". Required files: ", paste(required_files, collapse = ","),
      call. = FALSE
    )
  }

  eligible_seed_nums <- vapply(eligible$seed, seed_number, integer(1))
  eligible <- eligible[order(eligible$selection_objective, eligible_seed_nums, eligible$seed), , drop = FALSE]
  selected_seed <- eligible$seed[[1]]
  res$selected[as.character(res$seed) == as.character(selected_seed)] <- TRUE

  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(res$seed_dir[res$selected][[1]], con = best_dir_file)

  message("Selected seed: ", selected_seed)
  message("  seed_dir: ", res$seed_dir[res$selected][[1]])
  message("  objective: ", signif(res$selection_objective[res$selected][[1]], 8))
  message("  objective_source: ", res$selection_objective_source[res$selected][[1]])
  message("  summary: ", out_tsv)
  message("  best_dir_file: ", best_dir_file)
}

main()
