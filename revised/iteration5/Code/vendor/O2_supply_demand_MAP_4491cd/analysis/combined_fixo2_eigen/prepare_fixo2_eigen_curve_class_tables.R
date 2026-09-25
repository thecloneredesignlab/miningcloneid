#!/usr/bin/env Rscript

# Pure analysis consumer: annotate already-materialized FixO2 eigen-attractor
# coordinates with already-materialized dense-grid curve classes and slopes.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

WORKFLOW_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
UTIL_DIR <- file.path(WORKFLOW_DIR, "util")
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_path_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_cli_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_combined_fixo2_eigen_utils.R"))

read_csv_plain <- o2cfe_read_csv_plain
read_tsv_plain <- o2cfe_read_tsv_plain
write_tsv <- o2cfe_write_tsv

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = path, row.names = FALSE, na = "")
}

normalize_reductions <- o2cfe_normalize_reductions
normalize_variants <- o2cfe_normalize_variants

find_latest_class_table <- function(dense_grid_dir, explicit_path = NULL) {
  if (!is.null(explicit_path) && length(explicit_path) &&
      !is.na(explicit_path[[1L]]) && nzchar(explicit_path[[1L]])) {
    return(normalizePath(explicit_path[[1L]], mustWork = TRUE))
  }
  candidates <- list.files(
    dense_grid_dir,
    pattern = "^fixed_o2_ploidy_monotonicity_regression_by_seed\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(candidates)) {
    stop("No regression by-seed class table found under: ", dense_grid_dir, call. = FALSE)
  }
  candidates[[which.max(file.info(candidates)$mtime)]]
}

read_curve_class_map <- function(path, class_col = "curve_class") {
  tab <- read_tsv_plain(path)
  if (!"seed_number" %in% names(tab)) {
    if (!"seed_id" %in% names(tab)) {
      stop("Class table must contain seed_number or seed_id: ", path, call. = FALSE)
    }
    tab$seed_number <- suppressWarnings(as.integer(sub("^seed", "", as.character(tab$seed_id))))
  }
  if (!class_col %in% names(tab)) {
    stop("Class table is missing --class_col=", class_col, ": ", path, call. = FALSE)
  }
  tab$seed_number <- suppressWarnings(as.integer(tab$seed_number))
  tab[[class_col]] <- trimws(as.character(tab[[class_col]]))
  tab <- tab[is.finite(tab$seed_number) & nzchar(tab[[class_col]]), , drop = FALSE]
  if (!nrow(tab)) stop("Class table has no usable seed/class rows: ", path, call. = FALSE)
  dup <- duplicated(tab$seed_number)
  if (any(dup)) {
    dup_seeds <- unique(tab$seed_number[dup])
    inconsistent <- vapply(dup_seeds, function(seed) {
      length(unique(tab[[class_col]][tab$seed_number == seed])) > 1L
    }, logical(1))
    if (any(inconsistent)) {
      stop(
        "Class table has duplicate seed_number rows with conflicting classes: ",
        paste(dup_seeds[inconsistent], collapse = ", "),
        call. = FALSE
      )
    }
    tab <- tab[!duplicated(tab$seed_number), , drop = FALSE]
  }
  stats::setNames(tab[[class_col]], as.character(tab$seed_number))
}

read_average_slope_map <- function(path) {
  if (is.null(path) || !length(path) || is.na(path[[1L]]) || !nzchar(path[[1L]])) {
    return(numeric())
  }
  path <- normalizePath(path[[1L]], mustWork = TRUE)
  tab <- read_tsv_plain(path)
  if (!"seed_number" %in% names(tab)) {
    if (!"seed_id" %in% names(tab)) {
      stop("Average-slope table must contain seed_number or seed_id: ", path, call. = FALSE)
    }
    tab$seed_number <- suppressWarnings(as.integer(sub("^seed", "", as.character(tab$seed_id))))
  }
  if (!"average_slope" %in% names(tab)) {
    stop("Average-slope table is missing average_slope: ", path, call. = FALSE)
  }
  tab$seed_number <- suppressWarnings(as.integer(tab$seed_number))
  tab$average_slope <- suppressWarnings(as.numeric(tab$average_slope))
  tab <- tab[is.finite(tab$seed_number) & is.finite(tab$average_slope), , drop = FALSE]
  if (!nrow(tab)) stop("Average-slope table has no usable seed/slope rows: ", path, call. = FALSE)
  tab <- tab[!duplicated(tab$seed_number), , drop = FALSE]
  stats::setNames(tab$average_slope, as.character(tab$seed_number))
}

default_average_slope_table <- function(out_dir) {
  file.path(out_dir, "tables", "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv")
}

original_png_path <- function(coordinate_table) {
  path <- gsub("/Tables/", "/Figures/", coordinate_table, fixed = TRUE)
  sub("_coordinates\\.csv$", ".png", path)
}

discover_coordinate_tables <- function(pooled_root, reductions, variants) {
  rows <- list()
  for (reduction in reductions) {
    for (variant in variants) {
      table_dir <- file.path(pooled_root, reduction, "Tables", variant)
      if (!dir.exists(table_dir)) next
      files <- sort(list.files(table_dir, pattern = "_coordinates\\.csv$", full.names = TRUE))
      for (path in files) {
        rows[[length(rows) + 1L]] <- data.frame(
          reduction = reduction,
          variant = variant,
          coordinate_table = normalizePath(path, mustWork = FALSE),
          original_png = original_png_path(path),
          stub = sub("_coordinates\\.csv$", "", basename(path)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      reduction = character(), variant = character(), coordinate_table = character(),
      original_png = character(), stub = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

add_curve_class <- function(plot_data, class_map) {
  required <- c("dataset", "point_type", "seed")
  missing <- setdiff(required, names(plot_data))
  if (length(missing)) {
    stop("Coordinate table is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  plot_data$seed_number <- suppressWarnings(as.integer(plot_data$seed))
  plot_data$curve_class <- NA_character_
  idx <- plot_data$dataset == "invivo" & plot_data$point_type == "best" & is.finite(plot_data$seed_number)
  plot_data$curve_class[idx] <- unname(class_map[as.character(plot_data$seed_number[idx])])
  plot_data
}

add_average_slope <- function(plot_data, slope_map) {
  if (!"seed_number" %in% names(plot_data)) {
    plot_data$seed_number <- suppressWarnings(as.integer(plot_data$seed))
  }
  plot_data$average_slope <- NA_real_
  if (!length(slope_map)) return(plot_data)
  idx <- plot_data$dataset == "invivo" & plot_data$point_type == "best" & is.finite(plot_data$seed_number)
  plot_data$average_slope[idx] <- suppressWarnings(as.numeric(
    unname(slope_map[as.character(plot_data$seed_number[idx])])
  ))
  plot_data
}

analyze_embedding <- function(row, class_map, slope_map, average_slope_table, class_table, class_col, out_dir) {
  plot_data <- add_curve_class(read_csv_plain(row$coordinate_table), class_map)
  plot_data <- add_average_slope(plot_data, slope_map)
  table_dir <- file.path(out_dir, "tables", row$reduction, row$variant)
  annotated_path <- file.path(table_dir, paste0(row$stub, "_coordinates_curve_class.csv"))
  best_points_path <- file.path(table_dir, paste0(row$stub, "_best_points_curve_class.csv"))
  counts_path <- file.path(table_dir, paste0(row$stub, "_curve_class_counts.tsv"))
  write_csv(plot_data, annotated_path)
  best_points <- plot_data[plot_data$point_type == "best", , drop = FALSE]
  write_csv(best_points, best_points_path)
  counts <- as.data.frame(table(
    dataset = best_points$dataset,
    curve_class = ifelse(
      is.na(best_points$curve_class) | !nzchar(best_points$curve_class),
      "not_applicable_or_missing",
      best_points$curve_class
    ),
    useNA = "no"
  ), stringsAsFactors = FALSE)
  counts <- counts[counts$Freq > 0, , drop = FALSE]
  write_tsv(counts, counts_path)
  data.frame(
    reduction = row$reduction,
    variant = row$variant,
    stub = row$stub,
    coordinate_table = row$coordinate_table,
    annotated_coordinate_table = annotated_path,
    original_png = row$original_png,
    best_points_table = best_points_path,
    curve_class_counts_table = counts_path,
    average_slope_table = average_slope_table,
    class_table = class_table,
    class_col = class_col,
    stringsAsFactors = FALSE
  )
}

prepare_fixo2_eigen_curve_class_tables <- function(pooled_root,
                                                    dense_grid_dir,
                                                    out_dir,
                                                    class_table = NULL,
                                                    class_col = "curve_class",
                                                    average_slope_table = default_average_slope_table(out_dir),
                                                    reductions = c("PCAs", "UMAPs", "TSNEs"),
                                                    variants = c("Full", "BestOnly"),
                                                    dry_run = FALSE) {
  pooled_root <- normalizePath(path.expand(pooled_root), mustWork = TRUE)
  dense_grid_dir <- normalizePath(path.expand(dense_grid_dir), mustWork = TRUE)
  out_dir <- normalizePath(path.expand(out_dir), mustWork = FALSE)
  average_slope_table <- normalizePath(path.expand(average_slope_table), mustWork = !isTRUE(dry_run))
  class_table <- find_latest_class_table(dense_grid_dir, class_table)
  embeddings <- discover_coordinate_tables(pooled_root, reductions, variants)
  if (!nrow(embeddings)) stop("No pooled coordinate tables found under: ", pooled_root, call. = FALSE)
  if (isTRUE(dry_run)) return(invisible(embeddings))
  class_map <- read_curve_class_map(class_table, class_col)
  slope_map <- read_average_slope_map(average_slope_table)
  manifest <- do.call(rbind, lapply(seq_len(nrow(embeddings)), function(i) {
    analyze_embedding(
      embeddings[i, , drop = FALSE], class_map, slope_map, average_slope_table,
      class_table, class_col, out_dir
    )
  }))
  manifest_path <- file.path(out_dir, "tables", "pooled_embedding_curve_class_analysis_manifest.tsv")
  write_tsv(manifest, manifest_path)
  message("Wrote analysis manifest: ", manifest_path)
  invisible(manifest)
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  argv <- bpf_parse_args(raw_args)
  repo_root <- bpf_repo_root(SCRIPT_DIR)
  out_dir <- bpf_resolve_repo_path(
    argv$out_dir %||% file.path(bpf_combine_fixo2_eigen_attractor_result_dir(repo_root), "pooled_embedding_curve_class"),
    repo_root
  )
  prepare_fixo2_eigen_curve_class_tables(
    pooled_root = bpf_resolve_repo_path(
      argv$pooled_root %||% file.path(bpf_fixo2_eigen_attractor_result_dir(repo_root), "pooled_invivo_invitro"),
      repo_root,
      mustWork = TRUE
    ),
    dense_grid_dir = bpf_resolve_repo_path(argv$dense_grid_dir %||% bpf_dense_grid_result_root(repo_root), repo_root, mustWork = TRUE),
    out_dir = out_dir,
    class_table = argv$class_table,
    class_col = argv$class_col %||% "curve_class",
    average_slope_table = bpf_resolve_repo_path(
      argv$average_slope_table %||% default_average_slope_table(out_dir),
      repo_root,
      mustWork = !bpf_as_bool(argv$dry_run, FALSE)
    ),
    reductions = normalize_reductions(argv$reductions),
    variants = normalize_variants(argv$variants),
    dry_run = bpf_as_bool(argv$dry_run, FALSE)
  )
}

if (identical(environment(), globalenv())) {
  main()
}
