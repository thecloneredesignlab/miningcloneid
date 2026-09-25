#!/usr/bin/env Rscript

# Cross-layer path, CLI, and materialized-table contracts for the parameter
# landscape workflow. This file intentionally contains no model or plotting
# code so simulation, analysis, visualization, and reporting can share it.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) return(y)
  if (length(x) > 1L) return(x)
  if (is.na(x) || !nzchar(trimws(as.character(x)))) y else x
}

o2pl_source_file <- function() {
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "o2_supply_demand_map_parameter_landscape_io_utils.R"]
  if (length(own)) return(own[[length(own)]])
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE) else ""
}

o2pl_workflow_root <- function() {
  override <- Sys.getenv("O2PL_WORKFLOW_ROOT", unset = "")
  if (nzchar(override)) return(normalizePath(path.expand(override), mustWork = FALSE))
  starts <- unique(c(dirname(o2pl_source_file()), normalizePath(getwd(), mustWork = FALSE)))
  for (start in starts[nzchar(starts)]) {
    current <- start
    for (i in seq_len(10L)) {
      if (dir.exists(file.path(current, "model")) && dir.exists(file.path(current, "analysis"))) {
        return(normalizePath(current, mustWork = FALSE))
      }
      candidate <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP")
      if (dir.exists(candidate)) return(normalizePath(candidate, mustWork = FALSE))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate O2_supply_demand_MAP workflow root.", call. = FALSE)
}

o2pl_oxygen_root <- function() {
  normalizePath(file.path(o2pl_workflow_root(), "..", ".."), mustWork = FALSE)
}

o2pl_default_result_root <- function(layout = Sys.getenv("O2PL_DEFAULT_RESULT_LAYOUT", unset = "canonical")) {
  explicit <- Sys.getenv("O2PL_DEFAULT_RESULT_ROOT", unset = "")
  if (nzchar(explicit)) return(normalizePath(path.expand(explicit), mustWork = FALSE))
  layout <- tolower(trimws(layout %||% "canonical"))
  if (layout %in% c("legacy", "legacy_bpf", "best_fit_parameter_feature", "bpf")) {
    return(file.path(
      o2pl_oxygen_root(), "results", "analysis", "best_fit_parameter_feature",
      "02_parameter_landscape_clustering"
    ))
  }
  file.path(o2pl_oxygen_root(), "results", "analysis", "parameter_landscape")
}

default_oxygen_dir <- o2pl_oxygen_root
default_parameter_landscape_clustering_dir <- o2pl_default_result_root
default_paperfigures_dir <- o2pl_default_result_root

DEFAULT_INVIVO_INPUT_DIR <- file.path(o2pl_oxygen_root(), "results", "fit_invivo_O2_buffering_500seed")
DEFAULT_INVITRO_INPUT_DIR <- file.path(o2pl_oxygen_root(), "results", "fit_invitro_O2_buffering_500seed")
DEFAULT_OUTPUT_DIR <- o2pl_default_result_root()

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", kv[[1L]], fixed = TRUE)
    out[[key]] <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- tolower(trimws(as.character(x[[1L]])))
  if (value %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- suppressWarnings(as.integer(x[[1L]]))
  if (length(value) && is.finite(value)) value else default
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- suppressWarnings(as.numeric(x[[1L]]))
  if (length(value) && is.finite(value)) value else default
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(default)
  value <- suppressWarnings(as.numeric(trimws(unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE)))))
  value <- value[is.finite(value)]
  if (length(value)) value else default
}

as_char_vec <- function(x, default = character()) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(default)
  value <- trimws(unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE)))
  value <- value[nzchar(value)]
  if (length(value)) value else default
}

o2pl_normalize_dataset <- function(dataset) {
  dataset <- tolower(trimws(as.character(dataset %||% "invivo")))
  if (!dataset %in% c("invivo", "invitro")) stop("dataset must be invivo or invitro.", call. = FALSE)
  dataset
}

normalize_dataset <- o2pl_normalize_dataset
dataset_label <- function(dataset) if (identical(o2pl_normalize_dataset(dataset), "invivo")) "in vivo" else "in vitro"

o2pl_normalize_reduction <- function(reduction) {
  reduction <- gsub("-", "_", tolower(trimws(as.character(reduction %||% "umap"))), fixed = TRUE)
  if (identical(reduction, "t_sne")) reduction <- "tsne"
  if (!reduction %in% c("umap", "pca", "tsne")) stop("reduction must be umap, pca, or tsne.", call. = FALSE)
  reduction
}

normalize_reduction <- o2pl_normalize_reduction
reduction_dir_name <- function(reduction) switch(o2pl_normalize_reduction(reduction), umap = "UMAPs", pca = "PCAs", tsne = "TSNEs")
reduction_file_suffix <- function(reduction) o2pl_normalize_reduction(reduction)
reduction_coordinate_names <- function(reduction) switch(o2pl_normalize_reduction(reduction), umap = c("UMAP1", "UMAP2"), pca = c("PCA1", "PCA2"), tsne = c("tSNE1", "tSNE2"))
reduction_axis_labels <- function(reduction) switch(o2pl_normalize_reduction(reduction), umap = c("UMAP 1", "UMAP 2"), pca = c("PCA 1", "PCA 2"), tsne = c("t-SNE 1", "t-SNE 2"))
reduction_coordinate_cluster_dir <- function(reduction) switch(o2pl_normalize_reduction(reduction), umap = "ByUMAPCoordinates", pca = "ByPCACoordinates", tsne = "ByTSNECoordinates")
reduction_coordinate_cluster_label <- function(reduction) paste(reduction_coordinate_names(reduction), collapse = "_")

paper_dataset_dir <- function(dataset, root_dir = o2pl_default_result_root()) file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), o2pl_normalize_dataset(dataset))
paper_reduction_dir <- function(dataset, reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_dataset_dir(dataset, root_dir), reduction_dir_name(reduction))
paper_umap_dir <- function(dataset, root_dir = o2pl_default_result_root()) paper_reduction_dir(dataset, "umap", root_dir)
paper_reduction_tables_dir <- function(dataset, reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_reduction_dir(dataset, reduction, root_dir), "Tables")
paper_reduction_figures_dir <- function(dataset, reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_reduction_dir(dataset, reduction, root_dir), "Figures")
paper_reduction_tables_wclusters_dir <- function(dataset, reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_reduction_dir(dataset, reduction, root_dir), "TablesWclusters")
paper_reduction_figures_wclusters_dir <- function(dataset, reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_reduction_dir(dataset, reduction, root_dir), "FiguresWclusters")
paper_tables_dir <- function(dataset, root_dir = o2pl_default_result_root()) paper_reduction_tables_dir(dataset, "umap", root_dir)
paper_figures_dir <- function(dataset, root_dir = o2pl_default_result_root()) paper_reduction_figures_dir(dataset, "umap", root_dir)
paper_tables_wclusters_dir <- function(dataset, root_dir = o2pl_default_result_root()) paper_reduction_tables_wclusters_dir(dataset, "umap", root_dir)
paper_figures_wclusters_dir <- function(dataset, root_dir = o2pl_default_result_root()) paper_reduction_figures_wclusters_dir(dataset, "umap", root_dir)
paper_seed_parameter_tables_dir <- function(root_dir = o2pl_default_result_root()) normalizePath(path.expand(root_dir), mustWork = FALSE)
paper_seed_parameter_table_path <- function(dataset, table = c("best", "initial"), root_dir = o2pl_default_result_root()) {
  dataset <- o2pl_normalize_dataset(dataset)
  table <- match.arg(table)
  file.path(paper_seed_parameter_tables_dir(root_dir), if (table == "best") paste0(dataset, "_best_params_by_seed.csv") else paste0(dataset, "_deoptim_initial_population.csv"))
}
paper_best_params_csv <- function(dataset, root_dir = o2pl_default_result_root()) paper_seed_parameter_table_path(dataset, "best", root_dir)
paper_initial_population_csv <- function(dataset, root_dir = o2pl_default_result_root()) paper_seed_parameter_table_path(dataset, "initial", root_dir)
paper_fixo2_mode_tables_dir <- function(root_dir = o2pl_default_result_root()) file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), "FixO2Modes")
paper_pooled_dataset_dir <- function(root_dir = o2pl_default_result_root()) file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), "pooled_invivo_invitro")
paper_pooled_reduction_dir <- function(reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_pooled_dataset_dir(root_dir), reduction_dir_name(reduction))
paper_pooled_umap_dir <- function(root_dir = o2pl_default_result_root()) paper_pooled_reduction_dir("umap", root_dir)
paper_pooled_reduction_tables_dir <- function(reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_pooled_reduction_dir(reduction, root_dir), "Tables")
paper_pooled_reduction_figures_dir <- function(reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_pooled_reduction_dir(reduction, root_dir), "Figures")
paper_pooled_reduction_tables_wclusters_dir <- function(reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_pooled_reduction_dir(reduction, root_dir), "TablesWclusters")
paper_pooled_reduction_figures_wclusters_dir <- function(reduction = "umap", root_dir = o2pl_default_result_root()) file.path(paper_pooled_reduction_dir(reduction, root_dir), "FiguresWclusters")
paper_pooled_tables_dir <- function(root_dir = o2pl_default_result_root()) paper_pooled_reduction_tables_dir("umap", root_dir)
paper_pooled_figures_dir <- function(root_dir = o2pl_default_result_root()) paper_pooled_reduction_figures_dir("umap", root_dir)
paper_pooled_tables_wclusters_dir <- function(root_dir = o2pl_default_result_root()) paper_pooled_reduction_tables_wclusters_dir("umap", root_dir)
paper_pooled_figures_wclusters_dir <- function(root_dir = o2pl_default_result_root()) paper_pooled_reduction_figures_wclusters_dir("umap", root_dir)
paper_default_attractor_o2_grid <- function() seq(0, 5, by = 0.05)
paper_default_mode_summary_o2 <- function() c(0, 0.1, 0.5, 1, 2, 5)
default_dataset_input_dir <- function(dataset) if (o2pl_normalize_dataset(dataset) == "invivo") DEFAULT_INVIVO_INPUT_DIR else DEFAULT_INVITRO_INPUT_DIR

# Shared parameter metadata used to validate materialized feature tables and to
# label report tables. Numerical transformations remain in the analysis layer.
umap_parameter_set <- function(dataset) {
  dataset <- o2pl_normalize_dataset(dataset)
  if (dataset == "invitro") {
    return(c(
      "O2_crit", "mu_hp", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
      "buffer_n_exp", "n_O", "alpha_o2", "gamma_growth", "lam_max", "p_mis_base",
      "p_wgd", "gamma_mu"
    ))
  }
  c(
    "lam_max", "p_mis_base", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
    "buffer_n_exp", "p_wgd", "o2_S0", "kappa_O", "eta_o2", "alpha_o2",
    "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O", "k_clear"
  )
}
umap_log10_parameter_set <- function(dataset) intersect(c(
  "lam_max", "p_mis_base", "p_misseg", "k_o_mis", "buffer_beta", "buffer_n_exp",
  "p_wgd", "o2_S0", "kappa_O", "eta_o2", "alpha_o2", "mu_hp", "O2_crit", "k_clear"
), umap_parameter_set(dataset))

o2pl_simulation_tables_dir <- function(root_dir = o2pl_default_result_root(), dataset = NULL) {
  path <- file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), "simulation_tables")
  if (!is.null(dataset)) path <- file.path(path, o2pl_normalize_dataset(dataset))
  path
}
o2pl_analysis_tables_dir <- function(root_dir = o2pl_default_result_root(), scope = NULL) {
  path <- file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), "analysis_tables")
  if (!is.null(scope)) path <- file.path(path, as.character(scope))
  path
}
o2pl_vis_figures_dir <- function(root_dir = o2pl_default_result_root(), scope = NULL) {
  path <- file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), "figures")
  if (!is.null(scope)) path <- file.path(path, as.character(scope))
  path
}

read_tsv <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
read_csv_plain <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df %||% data.frame(), path, quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}
write_tsv_plain <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df %||% data.frame(), path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}
rbind_fill_plain <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x), rows)
  if (!length(rows)) return(data.frame())
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned <- lapply(rows, function(df) {
    for (name in setdiff(all_names, names(df))) df[[name]] <- rep(NA, nrow(df))
    df[, all_names, drop = FALSE]
  })
  do.call(rbind, aligned)
}

o2pl_require_columns <- function(df, required, artifact = "table") {
  missing <- setdiff(required, names(df))
  if (length(missing)) stop(artifact, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(df)
}

o2pl_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (!length(output)) return(NA_character_)
  sub("[[:space:]].*$", "", output[[1L]])
}

o2pl_schema_for <- function(df, artifact_id, descriptions = NULL) {
  descriptions <- descriptions %||% stats::setNames(rep("", ncol(df)), names(df))
  data.frame(
    artifact_id = artifact_id,
    column_index = seq_along(df),
    column_name = names(df),
    r_type = vapply(df, function(x) paste(class(x), collapse = "/"), character(1)),
    required = TRUE,
    description = unname(descriptions[names(df)] %||% rep("", ncol(df))),
    stringsAsFactors = FALSE
  )
}

o2pl_write_schema <- function(df, artifact_id, path, descriptions = NULL) {
  write_tsv_plain(o2pl_schema_for(df, artifact_id, descriptions), path)
}

o2pl_manifest_path <- function(root_dir = o2pl_default_result_root(), layer = c("simulation", "analysis", "visualization", "report")) {
  layer <- match.arg(layer)
  file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), paste0(layer, "_artifact_manifest.tsv"))
}

o2pl_record_artifact <- function(root_dir,
                                 layer,
                                 artifact_id,
                                 path,
                                 data = NULL,
                                 schema_path = NA_character_,
                                 producer,
                                 dataset = NA_character_,
                                 source = NA_character_,
                                 contract_version = "parameter_landscape_v1") {
  path <- normalizePath(path, mustWork = FALSE)
  row <- data.frame(
    artifact_id = artifact_id,
    layer = layer,
    dataset = dataset,
    path = path,
    rows = if (is.data.frame(data)) nrow(data) else NA_integer_,
    columns = if (is.data.frame(data)) ncol(data) else NA_integer_,
    sha256 = o2pl_sha256(path),
    schema_path = if (!is.na(schema_path)) normalizePath(schema_path, mustWork = FALSE) else NA_character_,
    producer = producer,
    source = source,
    contract_version = contract_version,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  manifest_path <- o2pl_manifest_path(root_dir, layer)
  existing <- if (file.exists(manifest_path)) read_tsv(manifest_path) else data.frame()
  if (nrow(existing) && "artifact_id" %in% names(existing)) existing <- existing[existing$artifact_id != artifact_id, , drop = FALSE]
  write_tsv_plain(rbind_fill_plain(list(existing, row)), manifest_path)
  invisible(row)
}

fixed_o2_o2_slug <- function(x) {
  text <- format(as.numeric(x), scientific = FALSE, trim = TRUE, digits = 12)
  text <- sub("0+$", "", text)
  text <- sub("\\.$", "", text)
  if (!nzchar(text)) text <- "0"
  gsub("-", "m", gsub("\\.", "p", text))
}
fixed_o2_reference_dir <- function(mode_tables_dir, mode_reference_o2) file.path(mode_tables_dir, paste0("reference_o2_", fixed_o2_o2_slug(mode_reference_o2)))
