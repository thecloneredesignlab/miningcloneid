#!/usr/bin/env Rscript

# Reconstruct the actual differential-evolution initial populations used by the
# 500 separate in vivo fits, remove the parameter-table row from each seed, and
# embed the remaining population together with the 500 best endpoints.
#
# The fit configuration records de_init_mode=hybrid, but the single-stage fit
# backend explicitly overrides that setting to uniform before DEoptim. The fit
# logs are therefore treated as the authority for the effective initialization
# mode. All 20 optimized variables are reconstructed in their original order so
# that the RNG stream matches the optimizer; rho_2N and sigma_burden are removed
# only after reconstruction, leaving the 18 mechanistic landscape parameters.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[pair[[1L]]]] <- paste(pair[-1L], collapse = "=")
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || !nzchar(x)) return(default)
  tolower(trimws(as.character(x[[1L]]))) %in% c("1", "true", "t", "yes", "y")
}

required_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  normalizePath(path, mustWork = TRUE)
}

sample_uniform_box <- function(n, lower, upper) {
  n <- as.integer(n)
  d <- length(lower)
  if (is.na(n) || n <= 0L || d <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  u <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  out <- sweep(u, 2L, as.numeric(upper - lower), `*`)
  out <- sweep(out, 2L, as.numeric(lower), `+`)
  colnames(out) <- names(lower)
  out
}

build_uniform_initial_population <- function(np, lower, upper, init_use) {
  # This intentionally preserves the otherwise redundant first uniform draw in
  # the fit backend before candidate 1 and candidates 2:NP are overwritten.
  pop <- sample_uniform_box(np, lower, upper)
  init_vec <- as.numeric(init_use[names(lower)])
  names(init_vec) <- names(lower)
  if (any(!is.finite(init_vec))) {
    stop("Initial vector contains missing/non-finite values.", call. = FALSE)
  }
  pop[1L, ] <- pmin(pmax(init_vec, lower), upper)
  if (np > 1L) {
    pop[2L:np, ] <- sample_uniform_box(np - 1L, lower, upper)
  }
  pop
}

inverse_transform <- function(param_name, values) {
  if (startsWith(param_name, "log10_")) 10^as.numeric(values) else as.numeric(values)
}

relabel_clusters_by_embedding <- function(cluster_num, plot_data) {
  centers <- stats::aggregate(
    cbind(tSNE1, tSNE2) ~ cluster,
    data = data.frame(
      cluster = as.integer(cluster_num),
      tSNE1 = as.numeric(plot_data$tSNE1),
      tSNE2 = as.numeric(plot_data$tSNE2)
    ),
    FUN = stats::median
  )
  centers <- centers[order(centers$tSNE1, centers$tSNE2), , drop = FALSE]
  relabel_map <- stats::setNames(seq_len(nrow(centers)), as.character(centers$cluster))
  as.integer(relabel_map[as.character(cluster_num)])
}

select_tsne_clusters <- function(best_coordinates, seed = 123L) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Required R package is not installed: cluster", call. = FALSE)
  }
  mat <- as.matrix(best_coordinates[, c("tSNE1", "tSNE2")])
  storage.mode(mat) <- "double"
  if (nrow(mat) != 500L || any(!is.finite(mat))) {
    stop("Clustering requires 500 finite best-endpoint coordinates.", call. = FALSE)
  }

  k_values <- 2:8
  distance <- stats::dist(mat)
  scores <- rep(NA_real_, length(k_values))
  centers_by_k <- vector("list", length(k_values))
  for (i in seq_along(k_values)) {
    k <- k_values[[i]]
    set.seed(as.integer(seed) + k)
    km <- try(
      stats::kmeans(mat, centers = k, nstart = 10L, iter.max = 100L),
      silent = TRUE
    )
    if (inherits(km, "try-error") || length(unique(km$cluster)) < 2L) next
    silhouette <- try(cluster::silhouette(km$cluster, distance), silent = TRUE)
    if (inherits(silhouette, "try-error")) next
    scores[[i]] <- mean(silhouette[, "sil_width"])
    centers_by_k[[i]] <- km$centers
  }
  diagnostics <- data.table(
    cluster_source = "invivo_best_tSNE1_tSNE2",
    k = k_values,
    average_silhouette = scores,
    sample_n = nrow(mat),
    selected = FALSE
  )
  valid <- which(is.finite(diagnostics$average_silhouette))
  if (!length(valid)) stop("No valid t-SNE k-means solution.", call. = FALSE)
  selected_i <- valid[[which.max(diagnostics$average_silhouette[valid])]]
  diagnostics[selected_i, selected := TRUE]
  selected_k <- diagnostics$k[[selected_i]]

  final_km <- try(
    stats::kmeans(
      mat,
      centers = centers_by_k[[selected_i]],
      iter.max = 100L,
      algorithm = "Lloyd"
    ),
    silent = TRUE
  )
  if (inherits(final_km, "try-error")) {
    set.seed(as.integer(seed) + selected_k + 1000L)
    final_km <- stats::kmeans(
      mat, centers = selected_k, nstart = 10L, iter.max = 100L
    )
  }
  cluster_num <- relabel_clusters_by_embedding(final_km$cluster, best_coordinates)
  list(
    cluster_num = cluster_num,
    cluster_id = sprintf("vi_C%02d", cluster_num),
    selected_k = length(unique(cluster_num)),
    selected_silhouette = diagnostics$average_silhouette[[selected_i]],
    diagnostics = diagnostics
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
fit_root <- normalizePath(
  path.expand(args$fit_root),
  mustWork = TRUE
)
data_dir <- normalizePath(
  path.expand(args$data_dir),
  mustWork = TRUE
)
overwrite <- as_bool(args$overwrite, FALSE)
n_core <- suppressWarnings(as.integer(if (is.null(args$n_core)) 1L else args$n_core))
if (!is.finite(n_core) || is.na(n_core) || n_core < 1L) n_core <- 1L

parameter_meta_path <- required_file(
  file.path(data_dir, "parameter_function_groups.tsv"),
  "18-parameter metadata"
)
best_path <- required_file(
  file.path(data_dir, "invivo_best_parameters_500seeds.tsv"),
  "500-seed best-parameter table"
)

output_paths <- c(
  initial_population = file.path(data_dir, "invivo_initial_population_18params.tsv.gz"),
  full_coordinates = file.path(data_dir, "invivo_initial_best_tsne_coordinates.tsv.gz"),
  best_coordinates = file.path(data_dir, "invivo_best_tsne_cluster_coordinates.tsv"),
  cluster_summary = file.path(data_dir, "invivo_tsne_cluster_summary.tsv"),
  cluster_diagnostics = file.path(data_dir, "invivo_tsne_clustering_diagnostics.tsv"),
  preprocessing = file.path(data_dir, "invivo_tsne_preprocessing_metadata.tsv"),
  run_metadata = file.path(data_dir, "invivo_tsne_run_metadata.tsv"),
  source_provenance = file.path(data_dir, "invivo_tsne_source_provenance.tsv")
)
if (!overwrite && all(file.exists(output_paths))) {
  message("Reusing complete in vivo t-SNE landscape outputs in ", data_dir)
  quit(save = "no", status = 0L)
}

parameter_meta <- fread(parameter_meta_path)
required_meta <- c("parameter", "transformation", "parameter_order")
if (!all(required_meta %in% names(parameter_meta))) {
  stop(
    "Parameter metadata is missing: ",
    paste(setdiff(required_meta, names(parameter_meta)), collapse = ", "),
    call. = FALSE
  )
}
setorder(parameter_meta, parameter_order)
target_parameters <- as.character(parameter_meta$parameter)
if (nrow(parameter_meta) != 18L ||
    !identical(parameter_meta$parameter_order, seq_len(18L)) ||
    any(c("rho_2N", "sigma_burden") %in% target_parameters)) {
  stop("The t-SNE target must be the ordered 18-parameter set excluding rho_2N and sigma_burden.", call. = FALSE)
}

seed_dirs <- file.path(fit_root, paste0("seed", seq_len(500L)))
seed_inputs <- unlist(lapply(
  seed_dirs,
  function(seed_dir) {
    file.path(
      seed_dir,
      c("fit_config.rds", "parameter_table.csv", "fit_status.log")
    )
  }
), use.names = FALSE)
missing_inputs <- seed_inputs[!file.exists(seed_inputs)]
if (length(missing_inputs)) {
  stop(
    "Missing separate in vivo fit input(s): ",
    paste(head(missing_inputs, 20L), collapse = ", "),
    call. = FALSE
  )
}

initial_rows <- vector("list", 500L)
initialization_rows <- vector("list", 500L)
reference_optimizer_names <- NULL
reference_parameter_signature <- NULL

for (seed in seq_len(500L)) {
  seed_dir <- seed_dirs[[seed]]
  cfg <- readRDS(file.path(seed_dir, "fit_config.rds"))
  parameter_table <- fread(file.path(seed_dir, "parameter_table.csv"))
  status_text <- paste(
    readLines(file.path(seed_dir, "fit_status.log"), warn = FALSE),
    collapse = "\n"
  )
  if (!grepl(
    "init_mode=uniform, warm_start=TRUE, init_local=0, init_uniform=255",
    status_text,
    fixed = TRUE
  )) {
    stop("Seed ", seed, " does not document the expected effective uniform initialization.", call. = FALSE)
  }
  if (!identical(as.integer(cfg$seed), as.integer(seed)) ||
      !identical(as.integer(cfg$NP), 256L)) {
    stop("Unexpected seed or NP in fit_config.rds for seed ", seed, ".", call. = FALSE)
  }

  opt <- parameter_table[estimate == TRUE]
  optimizer_names <- as.character(opt$param_name)
  prototypes <- as.character(opt$param_prototype)
  if (nrow(opt) != 20L ||
      !setequal(prototypes, c(target_parameters, "rho_2N", "sigma_burden"))) {
    stop("Seed ", seed, " does not contain the expected 20 optimized variables.", call. = FALSE)
  }
  parameter_signature <- paste(
    optimizer_names,
    signif(as.numeric(opt$init_value), 17L),
    signif(as.numeric(opt$lower_bound), 17L),
    signif(as.numeric(opt$upper_bound), 17L),
    prototypes,
    sep = ":",
    collapse = "|"
  )
  if (is.null(reference_optimizer_names)) {
    reference_optimizer_names <- optimizer_names
    reference_parameter_signature <- parameter_signature
  } else if (!identical(optimizer_names, reference_optimizer_names) ||
             !identical(parameter_signature, reference_parameter_signature)) {
    stop("Optimizer parameter order or bounds differ at seed ", seed, ".", call. = FALSE)
  }

  init <- setNames(as.numeric(opt$init_value), optimizer_names)
  lower <- setNames(as.numeric(opt$lower_bound), optimizer_names)
  upper <- setNames(as.numeric(opt$upper_bound), optimizer_names)
  set.seed(as.integer(seed))
  population_t <- build_uniform_initial_population(
    np = 256L,
    lower = lower,
    upper = upper,
    init_use = init
  )

  # The first row is the parameter-table initial vector and is excluded before
  # feature construction, as specified for the fitted-landscape input.
  population_t <- population_t[-1L, , drop = FALSE]
  natural <- data.table(
    seed_number = rep.int(as.integer(seed), nrow(population_t)),
    population_row = seq.int(2L, 256L)
  )
  for (parameter in target_parameters) {
    j <- match(parameter, prototypes)
    natural[[parameter]] <- inverse_transform(
      optimizer_names[[j]],
      population_t[, j]
    )
  }
  initial_rows[[seed]] <- natural
  initialization_rows[[seed]] <- data.table(
    seed_number = seed,
    configured_mode = as.character(cfg$de_init_mode),
    effective_mode = "uniform",
    NP = as.integer(cfg$NP),
    parameter_table_rows = 1L,
    retained_initial_rows = nrow(population_t),
    optimized_variables_reconstructed = nrow(opt),
    embedding_variables_retained = length(target_parameters),
    excluded_variables = "rho_2N,sigma_burden"
  )
}

initial <- rbindlist(initial_rows, use.names = TRUE)
initialization_audit <- rbindlist(initialization_rows, use.names = TRUE)
if (nrow(initial) != 127500L ||
    uniqueN(initial$seed_number) != 500L ||
    any(initial[, .N, by = seed_number]$N != 255L)) {
  stop("Expected 127,500 retained initial-population rows (500 x 255).", call. = FALSE)
}

best <- fread(best_path)
if (nrow(best) != 500L ||
    !identical(sort(as.integer(best$seed_number)), seq_len(500L)) ||
    anyDuplicated(best$seed_number)) {
  stop("Best-parameter table must contain seeds 1 through 500 exactly once.", call. = FALSE)
}
missing_best <- setdiff(target_parameters, names(best))
if (length(missing_best)) {
  stop("Best-parameter table is missing: ", paste(missing_best, collapse = ", "), call. = FALSE)
}
setorder(best, seed_number)

initial_features <- as.data.frame(initial[, ..target_parameters])
best_features <- as.data.frame(best[, ..target_parameters])
for (parameter_name in target_parameters) {
  transformation <- parameter_meta[
    parameter == parameter_name,
    transformation
  ][[1L]]
  if (identical(transformation, "log10")) {
    if (any(initial_features[[parameter_name]] <= 0) ||
        any(best_features[[parameter_name]] <= 0)) {
      stop("Cannot log10-transform non-positive values for ", parameter_name, ".", call. = FALSE)
    }
    initial_features[[parameter_name]] <- log10(initial_features[[parameter_name]])
    best_features[[parameter_name]] <- log10(best_features[[parameter_name]])
  } else if (!identical(transformation, "identity")) {
    stop(
      "Unsupported transformation for ",
      parameter_name,
      ": ",
      transformation,
      call. = FALSE
    )
  }
}
feature_matrix <- rbind(
  as.matrix(initial_features),
  as.matrix(best_features)
)
storage.mode(feature_matrix) <- "double"
scaled_matrix <- scale(feature_matrix)
centers <- attr(scaled_matrix, "scaled:center")
scales <- attr(scaled_matrix, "scaled:scale")
if (nrow(scaled_matrix) != 128000L ||
    ncol(scaled_matrix) != 18L ||
    any(!is.finite(scaled_matrix)) ||
    any(!is.finite(scales) | scales == 0)) {
  stop("The standardized t-SNE matrix must be finite with shape 128,000 x 18.", call. = FALSE)
}

if (!requireNamespace("Rtsne", quietly = TRUE)) {
  stop("Required R package is not installed: Rtsne", call. = FALSE)
}
Sys.setenv(OMP_NUM_THREADS = as.character(n_core))
set.seed(123L)
message("Running in vivo-only t-SNE: 128,000 rows x 18 parameters.")
embedding <- Rtsne::Rtsne(
  scaled_matrix,
  dims = 2L,
  perplexity = 30,
  theta = 0.5,
  max_iter = 1000L,
  pca = FALSE,
  check_duplicates = FALSE,
  num_threads = n_core,
  verbose = TRUE
)$Y
colnames(embedding) <- c("tSNE1", "tSNE2")
if (!identical(dim(embedding), c(128000L, 2L)) || any(!is.finite(embedding))) {
  stop("t-SNE returned an invalid embedding.", call. = FALSE)
}

full_coordinates <- rbind(
  data.table(
    tSNE1 = embedding[seq_len(nrow(initial)), 1L],
    tSNE2 = embedding[seq_len(nrow(initial)), 2L],
    dataset = "invivo",
    point_type = "initial",
    source_group = "invivo_initial",
    seed = initial$seed_number,
    population_row = initial$population_row,
    objective = NA_real_,
    reduction = "invivo_full_18param_tsne"
  ),
  data.table(
    tSNE1 = embedding[nrow(initial) + seq_len(nrow(best)), 1L],
    tSNE2 = embedding[nrow(initial) + seq_len(nrow(best)), 2L],
    dataset = "invivo",
    point_type = "best",
    source_group = "invivo_best",
    seed = best$seed_number,
    population_row = NA_integer_,
    objective = if ("objective_data" %in% names(best)) {
      as.numeric(best$objective_data)
    } else {
      as.numeric(best$objective)
    },
    reduction = "invivo_full_18param_tsne"
  )
)
best_coordinates <- full_coordinates[point_type == "best"]
cluster_result <- select_tsne_clusters(best_coordinates, seed = 123L)
best_coordinates[, `:=`(
  dataset_label = "in vivo",
  cluster_scope = "invivo_best",
  cluster_source = "invivo_best_tSNE1_tSNE2",
  cluster_prefix = "vi",
  cluster_base_id = sprintf("C%02d", cluster_result$cluster_num),
  cluster_id = cluster_result$cluster_id,
  cluster_num = cluster_result$cluster_num,
  cluster_k = cluster_result$selected_k,
  cluster_silhouette_avg = cluster_result$selected_silhouette,
  cluster_silhouette_sample_n = 500L
)]
setorder(best_coordinates, seed)

cluster_summary <- best_coordinates[, .(
  dataset = unique(dataset),
  dataset_label = unique(dataset_label),
  cluster_prefix = unique(cluster_prefix),
  cluster_base_id = unique(cluster_base_id),
  cluster_num = unique(cluster_num),
  n_seeds = .N,
  seed_min = min(seed),
  seed_max = max(seed),
  objective_mean = mean(objective),
  objective_median = median(objective),
  objective_min = min(objective),
  objective_max = max(objective),
  tSNE1_median = median(tSNE1),
  tSNE2_median = median(tSNE2),
  tSNE1_min = min(tSNE1),
  tSNE1_max = max(tSNE1),
  tSNE2_min = min(tSNE2),
  tSNE2_max = max(tSNE2),
  cluster_source = unique(cluster_source),
  selected_k = unique(cluster_k),
  selected_average_silhouette = unique(cluster_silhouette_avg)
), by = cluster_id]
setorder(cluster_summary, cluster_num)

preprocessing <- data.table(
  parameter = target_parameters,
  parameter_order = parameter_meta$parameter_order,
  transformation = parameter_meta$transformation,
  zscore_center = as.numeric(centers[target_parameters]),
  zscore_scale = as.numeric(scales[target_parameters]),
  used_in_tsne = TRUE
)
run_metadata <- rbindlist(list(
  data.table(
    run_key = c(
      "fit_scope", "seed_count", "optimizer_variables_reconstructed",
      "excluded_variables", "configured_initialization_mode",
      "effective_initialization_mode", "initial_population_rows_per_seed",
      "dropped_parameter_table_rows_per_seed", "retained_initial_rows",
      "best_endpoint_rows", "total_tsne_rows", "tsne_parameters",
      "tsne_seed", "tsne_perplexity", "tsne_theta", "tsne_max_iter",
      "cluster_seed", "cluster_k_candidates", "selected_cluster_k",
      "selected_average_silhouette"
    ),
    value = as.character(c(
      "separate_in_vivo_only", 500L, 20L,
      "rho_2N,sigma_burden", "hybrid", "uniform",
      256L, 1L, 127500L, 500L, 128000L, 18L,
      123L, 30, 0.5, 1000L, 123L, "2,3,4,5,6,7,8",
      cluster_result$selected_k, signif(cluster_result$selected_silhouette, 16L)
    ))
  ),
  initialization_audit[, .(
    run_key = paste0("seed", seed_number, "_initialization"),
    value = paste0(
      "configured=", configured_mode,
      ";effective=", effective_mode,
      ";NP=", NP,
      ";dropped=", parameter_table_rows,
      ";retained=", retained_initial_rows
    )
  )]
))
setnames(run_metadata, "run_key", "key")

source_paths <- c(parameter_meta_path, best_path, seed_inputs)
source_provenance <- data.table(
  role = c(
    "18-parameter embedding definition",
    "500 separate in vivo best endpoints",
    rep(
      c(
        "separate in vivo fit configuration",
        "separate in vivo optimizer parameter table",
        "separate in vivo effective initialization log"
      ),
      times = 500L
    )
  ),
  source = normalizePath(source_paths, mustWork = TRUE),
  source_md5 = unname(tools::md5sum(source_paths))
)

fwrite(initial, output_paths[["initial_population"]], sep = "\t", compress = "gzip")
fwrite(
  full_coordinates,
  output_paths[["full_coordinates"]],
  sep = "\t",
  compress = "gzip"
)
fwrite(best_coordinates, output_paths[["best_coordinates"]], sep = "\t")
fwrite(cluster_summary, output_paths[["cluster_summary"]], sep = "\t")
fwrite(
  cluster_result$diagnostics,
  output_paths[["cluster_diagnostics"]],
  sep = "\t"
)
fwrite(preprocessing, output_paths[["preprocessing"]], sep = "\t")
fwrite(run_metadata, output_paths[["run_metadata"]], sep = "\t")
fwrite(source_provenance, output_paths[["source_provenance"]], sep = "\t")

message(
  "Completed in vivo-only t-SNE: 127,500 initial + 500 best = 128,000 rows; ",
  cluster_result$selected_k, " best-endpoint clusters."
)
