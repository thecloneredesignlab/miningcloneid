#!/usr/bin/env Rscript

# Build a frozen proposal preconditioner for the Figure 5F generalized-posterior
# sampler. The 500 optimizer endpoints from each selected warm-start pair are
# used only to estimate proposal geometry. They do not enter the target density,
# prior, posterior draw table, or posterior weights.

resolve_own_file <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Unable to resolve this script path.")
  normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE)
}

script_path <- resolve_own_file()
iteration_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
selection_path <- file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv")
output_dir <- file.path(figure5_dir, "generalized_posterior")
output_rds <- file.path(output_dir, "figure5f_pilot2_proposal_covariance.rds")
output_tsv <- file.path(output_dir, "figure5f_endpoint_proposal_diagnostics.tsv")
output_manifest <- file.path(output_dir, "figure5f_endpoint_proposal_sources.tsv")
result_root <- Sys.getenv(
  "FIGURE5F_RESULT_ROOT",
  unset = paste0(
    "/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/",
    "fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633"
  )
)
result_root <- normalizePath(result_root, mustWork = TRUE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

hash_lines <- function(x) {
  path <- tempfile("figure5f_hash_lines_")
  on.exit(unlink(path), add = TRUE)
  writeLines(as.character(x), path, useBytes = TRUE)
  unname(tools::md5sum(path))
}

selection <- utils::read.delim(selection_path, check.names = FALSE, stringsAsFactors = FALSE)
selection <- selection[selection$selected_for_figure5f %in% TRUE, , drop = FALSE]
families <- sprintf("C%02d", 1:6)
selection <- selection[order(match(selection$family, families)), , drop = FALSE]
if (nrow(selection) != 6L || !identical(selection$family, families)) {
  stop("Expected one selected warm-start pair for each of C01-C06.")
}

selected_pair_dir <- function(row) {
  file.path(result_root, paste0("fit_joint_", row$warmup_label[[1L]]))
}

template_dir <- file.path(selected_pair_dir(selection[1L, , drop = FALSE]), "seed1")
required_template <- file.path(
  template_dir,
  c(
    "best_params_transformed.tsv",
    "joint_soft_coupling_parameters_table_input.csv",
    "joint_shared_bounds.tsv",
    "parameter_table_invivo_transformed.csv",
    "parameter_table_invitro_transformed.csv",
    "run_effective_args.tsv"
  )
)
missing_template <- required_template[!file.exists(required_template)]
if (length(missing_template)) stop("Missing proposal template input(s):\n", paste(missing_template, collapse = "\n"))

template_par <- utils::read.delim(required_template[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
soft_input <- utils::read.csv(required_template[[2L]], check.names = FALSE, stringsAsFactors = FALSE)
shared_bounds <- utils::read.delim(required_template[[3L]], check.names = FALSE, stringsAsFactors = FALSE)
vivo_bounds <- utils::read.csv(required_template[[4L]], check.names = FALSE, stringsAsFactors = FALSE)
vitro_bounds <- utils::read.csv(required_template[[5L]], check.names = FALSE, stringsAsFactors = FALSE)
effective_args <- utils::read.delim(required_template[[6L]], check.names = FALSE, stringsAsFactors = FALSE)

parameter_from_center <- function(x) sub("^log10_", "", x)
available_centers <- soft_input$param_name[!startsWith(soft_input$param_name, "delta__")]
soft_order_value <- effective_args$value[
  effective_args$key == "joint_soft_coupling_params"
][[1L]]
soft_parameter_order <- strsplit(soft_order_value, ",", fixed = TRUE)[[1L]]
center_lookup <- setNames(available_centers, parameter_from_center(available_centers))
center_names <- unname(center_lookup[soft_parameter_order])
delta_names <- paste0("delta__", center_names)
if (length(center_names) != 14L || anyNA(center_names) ||
    any(!delta_names %in% soft_input$param_name)) {
  stop("Expected 14 center/delta soft-coupled parameter pairs.")
}
ordinary_names <- setdiff(template_par$transformed_parameter, c(center_names, delta_names))
if (length(ordinary_names) != 12L) stop("Expected 12 ordinary joint coordinates.")

coord_rows <- list()
row_i <- 1L
for (name in ordinary_names) {
  if (startsWith(name, "ivt__")) {
    lookup <- sub("^ivt__", "", name)
    hit <- vitro_bounds[vitro_bounds$param_name == lookup, , drop = FALSE]
    role <- "ordinary_invitro"
  } else {
    hit <- vivo_bounds[vivo_bounds$param_name == name, , drop = FALSE]
    role <- "ordinary_invivo"
  }
  if (nrow(hit) != 1L) stop("Missing ordinary-coordinate bounds for ", name)
  coord_rows[[row_i]] <- data.frame(
    coordinate = name,
    role = role,
    parameter = name,
    original_name = name,
    lower = as.numeric(hit$lower_bound %||% hit$lower),
    upper = as.numeric(hit$upper_bound %||% hit$upper),
    stringsAsFactors = FALSE
  )
  row_i <- row_i + 1L
}
for (center_name in center_names) {
  parameter <- parameter_from_center(center_name)
  hit <- shared_bounds[shared_bounds$parameter == parameter, , drop = FALSE]
  if (nrow(hit) != 1L) stop("Missing joint-union bounds for ", parameter)
  for (role in c("vivo", "vitro")) {
    coord_rows[[row_i]] <- data.frame(
      coordinate = paste0(role, "__", center_name),
      role = role,
      parameter = parameter,
      original_name = center_name,
      lower = as.numeric(hit$joint_union_lower_t),
      upper = as.numeric(hit$joint_union_upper_t),
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
}
coord <- do.call(rbind, coord_rows)
coord$width <- coord$upper - coord$lower
if (nrow(coord) != 40L || any(!is.finite(coord$width)) || any(coord$width <= 0)) {
  stop("Invalid 40-dimensional proposal-coordinate contract.")
}

read_endpoint_unit <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  par <- setNames(as.numeric(tab$transformed_value), tab$transformed_parameter)
  if (!setequal(names(par), template_par$transformed_parameter)) {
    stop("Endpoint parameter schema differs: ", path)
  }
  actual <- setNames(numeric(nrow(coord)), coord$coordinate)
  for (name in ordinary_names) actual[[name]] <- par[[name]]
  for (center_name in center_names) {
    center <- par[[center_name]]
    delta <- par[[paste0("delta__", center_name)]]
    actual[[paste0("vivo__", center_name)]] <- center + delta / 2
    actual[[paste0("vitro__", center_name)]] <- center - delta / 2
  }
  unit <- (actual[coord$coordinate] - coord$lower) / coord$width
  if (any(!is.finite(unit)) || any(unit < -1e-8) || any(unit > 1 + 1e-8)) {
    stop("Endpoint lies outside the reconstructed joint box: ", path)
  }
  pmin(pmax(unit, 0), 1)
}

endpoint_blocks <- list()
source_rows <- list()
for (i in seq_len(nrow(selection))) {
  family <- selection$family[[i]]
  pair_dir <- selected_pair_dir(selection[i, , drop = FALSE])
  paths <- file.path(pair_dir, paste0("seed", seq_len(500L)), "best_params_transformed.tsv")
  if (any(!file.exists(paths))) {
    stop("Selected pair does not contain all 500 endpoint tables: ", pair_dir)
  }
  block <- t(vapply(paths, read_endpoint_unit, numeric(nrow(coord))))
  colnames(block) <- coord$coordinate
  center <- apply(block, 2L, stats::median)
  endpoint_blocks[[family]] <- sweep(block, 2L, center, "-")
  source_rows[[i]] <- data.frame(
    family = family,
    warmup_label = selection$warmup_label[[i]],
    n_endpoint_files = length(paths),
    source_directory = normalizePath(pair_dir, mustWork = TRUE),
    endpoint_md5_manifest = hash_lines(unname(tools::md5sum(paths))),
    stringsAsFactors = FALSE
  )
}

residual <- do.call(rbind, endpoint_blocks)
raw_cov <- stats::cov(residual)
minimum_variance <- 0.01^2
diagonal_target <- diag(pmax(diag(raw_cov), minimum_variance), nrow(raw_cov))
shrinkage <- 0.50
base_cov <- (1 - shrinkage) * raw_cov + shrinkage * diagonal_target
ridge <- max(diag(base_cov)) * 1e-6
base_cov <- base_cov + diag(ridge, nrow(base_cov))
dimnames(base_cov) <- list(coord$coordinate, coord$coordinate)

temperatures <- c(0.5, 1, 2, 4, 8)
covariance <- setNames(lapply(temperatures, function(temp) base_cov * temp), as.character(temperatures))
minimum_eigenvalue <- min(eigen(base_cov, symmetric = TRUE, only.values = TRUE)$values)
if (!is.finite(minimum_eigenvalue) || minimum_eigenvalue <= 0) {
  stop("Endpoint-derived proposal covariance is not positive definite.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source_table <- do.call(rbind, source_rows)
object <- list(
  covariance = covariance,
  coordinate_names = coord$coordinate,
  temperatures = temperatures,
  shrinkage_to_diagonal = shrinkage,
  minimum_unit_sd = sqrt(minimum_variance),
  ridge = ridge,
  endpoint_source_manifest = source_table,
  selection_md5 = unname(tools::md5sum(selection_path)),
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  scientific_draws_included = FALSE,
  role = "frozen endpoint-derived proposal geometry only"
)
temporary <- paste0(output_rds, ".tmp_", Sys.getpid())
saveRDS(object, temporary)
if (!file.rename(temporary, output_rds)) stop("Failed to save proposal covariance.")

diagnostics <- data.frame(
  n_families = length(endpoint_blocks),
  n_endpoints = nrow(residual),
  n_coordinates = ncol(residual),
  shrinkage_to_diagonal = shrinkage,
  minimum_unit_sd = sqrt(minimum_variance),
  ridge = ridge,
  minimum_eigenvalue = minimum_eigenvalue,
  maximum_eigenvalue = max(eigen(base_cov, symmetric = TRUE, only.values = TRUE)$values),
  proposal_covariance_md5 = unname(tools::md5sum(output_rds)),
  scientific_draws_included = FALSE,
  stringsAsFactors = FALSE
)
utils::write.table(diagnostics, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
utils::write.table(source_table, output_manifest, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Saved frozen endpoint-derived proposal preconditioner:\n", output_rds, "\n", sep = "")
print(diagnostics)
