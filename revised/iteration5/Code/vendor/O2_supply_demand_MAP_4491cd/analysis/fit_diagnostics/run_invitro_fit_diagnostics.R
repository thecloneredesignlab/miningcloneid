#!/usr/bin/env Rscript

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_postfit_io_utils.R"),
  local = environment()
)
rm(.script_dir)

write_tsv <- ivt_sim_write_tsv

read_optimizer_population <- function(path) {
  if (!file.exists(path)) return(NULL)
  population <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(population) || !is.data.frame(population)) return(NULL)
  population$optimizer_sample_id <- NULL
  if (!ncol(population)) return(NULL)
  population[] <- lapply(population, function(x) suppressWarnings(as.numeric(x)))
  keep <- vapply(
    population,
    function(x) sum(is.finite(x)) >= 3L && stats::var(x, na.rm = TRUE) > 0,
    logical(1)
  )
  population <- population[, keep, drop = FALSE]
  population <- population[stats::complete.cases(population), , drop = FALSE]
  if (nrow(population) >= 3L && ncol(population) >= 2L) population else NULL
}

empty_diagnostic_tables <- function() {
  list(
    invitro_identifiability_variance = data.frame(
      direction = character(),
      variance = numeric(),
      variance_fraction = numeric(),
      weakest_rank = integer(),
      stringsAsFactors = FALSE
    ),
    invitro_identifiability_loadings = data.frame(
      direction = character(),
      parameter = character(),
      loading = numeric(),
      absolute_loading = numeric(),
      loading_rank = integer(),
      stringsAsFactors = FALSE
    ),
    invitro_parameter_correlation = data.frame(
      parameter_x = character(),
      parameter_y = character(),
      correlation = numeric(),
      stringsAsFactors = FALSE
    ),
    invitro_fit_diagnostics_summary = data.frame(
      metric = character(),
      value = character(),
      stringsAsFactors = FALSE
    )
  )
}

compute_diagnostic_tables <- function(population) {
  tables <- empty_diagnostic_tables()
  if (is.null(population)) {
    tables$invitro_fit_diagnostics_summary <- data.frame(
      metric = "status",
      value = "optimizer_population_unavailable",
      stringsAsFactors = FALSE
    )
    return(tables)
  }
  scaled <- scale(population)
  pca <- stats::prcomp(scaled, center = FALSE, scale. = FALSE)
  variance <- pca$sdev^2
  variance_fraction <- variance / sum(variance)
  direction <- paste0("PC", seq_along(variance))
  weakest_order <- order(variance, decreasing = FALSE)
  weakest_rank <- match(seq_along(variance), weakest_order)
  tables$invitro_identifiability_variance <- data.frame(
    direction = direction,
    variance = variance,
    variance_fraction = variance_fraction,
    weakest_rank = weakest_rank,
    stringsAsFactors = FALSE
  )

  weak_index <- utils::head(weakest_order, min(5L, length(weakest_order)))
  loading_rows <- lapply(weak_index, function(index) {
    loading <- pca$rotation[, index]
    ordered <- order(abs(loading), decreasing = TRUE)
    ordered <- utils::head(ordered, min(8L, length(ordered)))
    data.frame(
      direction = paste0("PC", index),
      parameter = names(loading)[ordered],
      loading = loading[ordered],
      absolute_loading = abs(loading[ordered]),
      loading_rank = seq_along(ordered),
      stringsAsFactors = FALSE
    )
  })
  tables$invitro_identifiability_loadings <- do.call(rbind, loading_rows)

  correlation <- stats::cor(population, use = "pairwise.complete.obs")
  correlation_grid <- expand.grid(
    parameter_x = colnames(correlation),
    parameter_y = colnames(correlation),
    stringsAsFactors = FALSE
  )
  correlation_grid$correlation <- as.numeric(correlation)
  tables$invitro_parameter_correlation <- correlation_grid

  positive_variance <- variance[is.finite(variance) & variance > 0]
  condition_proxy <- if (length(positive_variance)) {
    max(positive_variance) / min(positive_variance)
  } else {
    NA_real_
  }
  tables$invitro_fit_diagnostics_summary <- data.frame(
    metric = c(
      "status", "optimizer_population_rows", "optimizer_parameter_count",
      "pca_direction_count", "variance_condition_proxy"
    ),
    value = c(
      "ok",
      as.character(nrow(population)),
      as.character(ncol(population)),
      as.character(length(variance)),
      format(condition_proxy, digits = 17)
    ),
    stringsAsFactors = FALSE
  )
  tables
}

schema_for_tables <- function(tables) {
  do.call(rbind, lapply(names(tables), function(name) {
    table <- tables[[name]]
    data.frame(
      table = paste0(name, ".tsv"),
      column = names(table),
      class = vapply(table, function(x) paste(class(x), collapse = ","), character(1)),
      rows = nrow(table),
      stringsAsFactors = FALSE
    )
  }))
}

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% stop(
    "Usage: run_invitro_fit_diagnostics.R --fit_dir=/abs/seed [--simulation_dir=/abs/simulation/invitro] [--analysis_dir=/abs/output]",
    call. = FALSE
  )
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  simulation_dir <- argv$simulation_dir %||% file.path(fit_dir, "simulation", "invitro")
  simulation_dir <- normalizePath(simulation_dir, mustWork = FALSE)
  analysis_dir <- argv$analysis_dir %||% argv$out_dir %||%
    file.path(fit_dir, "analysis", "invitro")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)

  population_path <- file.path(simulation_dir, "invitro_optimizer_population.tsv")
  population <- read_optimizer_population(population_path)
  tables <- compute_diagnostic_tables(population)

  for (name in names(tables)) {
    write_tsv(tables[[name]], file.path(analysis_dir, paste0(name, ".tsv")))
  }
  write_tsv(
    schema_for_tables(tables),
    file.path(analysis_dir, "invitro_fit_diagnostics_schema.tsv")
  )
  manifest <- data.frame(
    key = c(
      "fit_dir", "simulation_dir", "analysis_dir", "optimizer_population", "generated_at", "status",
      "table_count", paste0("rows.", names(tables))
    ),
    value = c(
      fit_dir,
      simulation_dir,
      analysis_dir,
      normalizePath(population_path, mustWork = FALSE),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      tables$invitro_fit_diagnostics_summary$value[
        match("status", tables$invitro_fit_diagnostics_summary$metric)
      ],
      as.character(length(tables)),
      as.character(vapply(tables, nrow, integer(1)))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(analysis_dir, "invitro_fit_diagnostics_manifest.tsv"))
  message("In-vitro fit diagnostics written to: ", analysis_dir)
  invisible(analysis_dir)
}

if (sys.nframe() == 0L) main()
