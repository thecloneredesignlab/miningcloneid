#!/usr/bin/env Rscript

# Figure 6 cluster-selection, objective-eligibility, multi-seed response-surface,
# robustness, and drawing utilities.
#
# Scientific-input boundary: fitting roots are read-only inputs; all executable
# code and generated outputs are below revised/iteration4. Model evaluation
# loads the required read-only implementation from soft_couping_org.

options(stringsAsFactors = FALSE, warn = 1)

f6r_family_levels <- function() {
  if (exists("JOINT_FAMILY_LEVELS", inherits = TRUE)) {
    levels <- get("JOINT_FAMILY_LEVELS", inherits = TRUE)
  } else {
    joint_root <- Sys.getenv("FIGURE_JOINT_RESULT_ROOT", unset = "")
    manifest_path <- file.path(joint_root, "multi_warmup_manifest.tsv")
    if (!nzchar(joint_root) || !file.exists(manifest_path)) {
      stop("Cannot resolve the joint primary-family manifest.")
    }
    manifest <- utils::read.delim(
      manifest_path, check.names = FALSE, stringsAsFactors = FALSE
    )
    levels <- regmatches(
      manifest$warmup_label, regexpr("C[0-9]{2}", manifest$warmup_label)
    )
    levels <- levels[order(as.integer(sub("^C", "", levels)))]
  }
  expected <- sprintf("C%02d", seq_along(levels))
  if (!length(levels) || !identical(as.character(levels), expected)) {
    stop("Joint primary-family labels are missing or nonconsecutive.")
  }
  as.character(levels)
}

f6r_family_count <- function() length(f6r_family_levels())

f6r_find_workspace_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "manager.sh")) &&
        dir.exists(file.path(current, "Code", "Figures")) &&
        dir.exists(file.path(current, "data", "Figures"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the iteration4 figure workspace from: ", start)
    }
    current <- parent
  }
}

f6r_paths <- function(workspace_root = f6r_find_workspace_root()) {
  root <- normalizePath(workspace_root, mustWork = TRUE)
  repository_root <- root
  default_model_code_root <-
    "/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP"
  model_code_root <- trimws(Sys.getenv(
    "FIGURE_MODEL_CODE_ROOT", unset = default_model_code_root
  ))
  if (!nzchar(model_code_root)) model_code_root <- default_model_code_root
  list(
    root = root,
    repository_root = repository_root,
    oxygen_code = normalizePath(model_code_root, mustWork = TRUE),
    code = file.path(root, "Code", "Figures"),
    figure6 = file.path(root, "data", "Figures", "Figure6"),
    figure4 = file.path(root, "data", "Figures", "Figure4"),
    figure5 = file.path(root, "data", "Figures", "Figure5"),
    supp5_1 = file.path(root, "data", "Figures", "Supp_Figure5_1"),
    supp5_2 = file.path(root, "data", "Figures", "Supp_Figure5_2"),
    supp6_2 = file.path(root, "data", "Figures", "Supp_Figure6_2"),
    figures = file.path(root, "Figures"),
    manuscript_figures = file.path(root, "manuscript", "Figures"),
    manuscript = file.path(root, "manuscript", "ltee_hypoxia_model.tex")
  )
}

f6r_require_packages <- function(packages) {
  missing <- packages[!vapply(
    packages, requireNamespace, logical(1L), quietly = TRUE
  )]
  if (length(missing)) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

f6r_resilient_lapply <- function(X, FUN, n_core = 1L) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
  )
  result <- if (n_core > 1L) {
    f6r_require_packages(c("future", "future.apply"))
    # Multisession workers are independent R processes rather than forked
    # children. This permits the multi-seed, dense-grid, and inverse stages to
    # create successive pools without exhausting macOS/OpenMP shared memory.
    previous_plan <- future::plan()
    on.exit(future::plan(previous_plan), add = TRUE)
    previous_worker_limit <- getOption("parallelly.maxWorkers.localhost")
    on.exit(
      options(parallelly.maxWorkers.localhost = previous_worker_limit),
      add = TRUE
    )
    options(
      parallelly.maxWorkers.localhost = max(3, as.integer(n_core))
    )
    future::plan(
      future::multisession,
      workers = min(as.integer(n_core), length(X))
    )
    future.apply::future_lapply(
      X, FUN,
      future.seed = TRUE,
      future.scheduling = Inf
    )
  } else {
    lapply(X, FUN)
  }
  missing <- which(vapply(result, is.null, logical(1L)))
  if (length(missing)) {
    message(
      "Retrying ", length(missing),
      " fork worker result(s) sequentially after an empty return."
    )
    result[missing] <- lapply(X[missing], FUN)
  }
  result
}

f6r_require_files <- function(paths, label = "required input") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing ", label, " file(s):\n", paste(missing, collapse = "\n"))
  }
  invisible(normalizePath(paths, mustWork = TRUE))
}

f6r_read_tsv <- function(path) {
  f6r_require_files(path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

f6r_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, file = path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  normalizePath(path, mustWork = TRUE)
}

f6r_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)[[1L]])
}

f6r_model_source_files <- function(paths) {
  files <- file.path(
    paths$oxygen_code,
    c(
      "model/model_O2_supply_demand_MAP.R",
      "model/model_O2_supply_demand_MAP.cpp",
      "util/o2_supply_demand_map_shared.R",
      "util/o2_supply_demand_map_common_semantics.R",
      "simulation/fix_o2_simulation.R",
      "simulation/fix_o2_simulation_shared_utils.R",
      "simulation/o2/fixed_o2/run_fixed_o2_simulation.R"
    )
  )
  f6r_require_files(files, "external model implementation")
  normalizePath(files, mustWork = TRUE)
}

f6r_model_source_fingerprint <- function(paths) {
  files <- f6r_model_source_files(paths)
  paste(
    paste(basename(files), unname(tools::md5sum(files)), sep = "="),
    collapse = ";"
  )
}

f6r_adjusted_rand <- function(a, b) {
  if (length(a) != length(b) || length(a) < 2L) return(NA_real_)
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  nij <- sum(choose2(tab))
  ai <- sum(choose2(rowSums(tab)))
  bj <- sum(choose2(colSums(tab)))
  total <- choose2(sum(tab))
  if (!is.finite(total) || total <= 0) return(NA_real_)
  expected <- ai * bj / total
  maximum <- 0.5 * (ai + bj)
  denom <- maximum - expected
  if (abs(denom) < .Machine$double.eps) return(1)
  (nij - expected) / denom
}

f6r_silhouette <- function(labels, distance_object) {
  mean(cluster::silhouette(as.integer(labels), distance_object)[, "sil_width"])
}

f6r_shared_parameters <- function() {
  c(
    "lam_max", "p_mis_base", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp", "p_wgd",
    "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O"
  )
}

f6r_endpoint_signatures <- function(
    master, value_column = "vivo_natural",
    signature_column = "parameter_signature"
) {
  f6r_require_packages("data.table")
  shared <- f6r_shared_parameters()
  x <- data.table::as.data.table(master)[parameter %in% shared]
  if (!value_column %in% names(x)) {
    stop("Missing endpoint-signature value column: ", value_column)
  }
  signatures <- x[, {
    ordered_values <- get(value_column)[match(shared, parameter)]
    if (length(ordered_values) != length(shared) || any(!is.finite(ordered_values))) {
      stop("Cannot construct a finite 14-parameter endpoint signature.")
    }
    list(parameter_signature = paste(sprintf("%.17g", ordered_values), collapse = "|"))
  }, by = .(pair_id, seed_number)]
  expected_joint_seed_count <- 500L * f6r_family_count()
  if (nrow(signatures) != expected_joint_seed_count ||
      anyDuplicated(signatures[, .(pair_id, seed_number)])) {
    stop(
      "Expected one exact endpoint signature for each of ",
      expected_joint_seed_count, " joint seeds."
    )
  }
  data.table::setnames(signatures, "parameter_signature", signature_column)
  as.data.frame(signatures)
}

f6r_endpoint_parameter_matrix <- function(paths) {
  long <- f6r_read_tsv(file.path(
    paths$figure4, "all_parameter_fitted_endpoint_values.tsv"
  ))
  shared <- f6r_shared_parameters()
  long <- long[
    long$parameter %in% shared,
    c("seed_number", "parameter", "transformed_value"),
    drop = FALSE
  ]
  if (nrow(long) != 500L * length(shared) ||
      any(!is.finite(long$transformed_value))) {
    stop("The local 14-parameter endpoint matrix is incomplete.")
  }
  wide <- reshape(
    long, idvar = "seed_number", timevar = "parameter", direction = "wide"
  )
  names(wide) <- sub("transformed_value.", "", names(wide), fixed = TRUE)
  wide <- wide[order(wide$seed_number), c("seed_number", shared), drop = FALSE]
  rownames(wide) <- NULL
  wide
}

# The superseded secondary-cluster diagnostic implementation was removed here.
# The active implementation below reconstructs only the primary families used
# by the current joint fitting.
f6r_cluster_diagnostics <- function(paths, n_resample = 100L) {
  f6r_require_packages(c("cluster", "data.table"))
  dir.create(paths$figure6, recursive = TRUE, showWarnings = FALSE)

  endpoint_path <- file.path(
    paths$figure5,
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
  )
  endpoints <- utils::read.csv(
    endpoint_path, check.names = FALSE, stringsAsFactors = FALSE
  )
  endpoints <- endpoints[
    endpoints$dataset == "invivo",
    c("seed", "objective", "cluster_base_id", "cluster_k",
      "cluster_silhouette_avg", "tSNE1", "tSNE2"),
    drop = FALSE
  ]
  endpoints <- endpoints[order(endpoints$seed), , drop = FALSE]
  expected_regions <- f6r_family_levels()
  selected_k <- length(expected_regions)
  observed_counts <- table(factor(
    endpoints$cluster_base_id, levels = expected_regions
  ))
  if (nrow(endpoints) != 500L || anyDuplicated(endpoints$seed) ||
      any(observed_counts < 1L) ||
      !all(endpoints$cluster_k == selected_k)) {
    stop("Expected the saved 500-endpoint primary partition at k=", selected_k, ".")
  }

  coords <- as.matrix(endpoints[, c("tSNE1", "tSNE2")])
  coord_dist <- stats::dist(coords)
  saved_labels <- match(endpoints$cluster_base_id, expected_regions)
  saved_silhouette <- f6r_silhouette(saved_labels, coord_dist)
  primary_rows <- lapply(2:8, function(k) {
    set.seed(6600L + k)
    km <- stats::kmeans(coords, centers = k, nstart = 100L, iter.max = 200L)
    sizes <- as.integer(table(km$cluster))
    data.frame(
      analysis_level = "primary pooled t-SNE regions",
      region = "all in-vivo endpoints",
      feature_space = "two-dimensional pooled 14-parameter t-SNE coordinates",
      k = k,
      average_silhouette = f6r_silhouette(km$cluster, coord_dist),
      n_solution = nrow(coords),
      minimum_cluster_size = min(sizes),
      maximum_cluster_size = max(sizes),
      selected_for_warm_starts = k == selected_k,
      stringsAsFactors = FALSE
    )
  })
  k_selection <- do.call(rbind, primary_rows)
  k_selection$maximum_silhouette <- max(k_selection$average_silhouette)
  k_selection$difference_from_maximum <-
    k_selection$maximum_silhouette - k_selection$average_silhouette
  k_selection$within_0p001_of_maximum <-
    k_selection$difference_from_maximum <= 0.001 + 1e-12
  k_selection$selection_basis <- ifelse(
    k_selection$selected_for_warm_starts,
    paste0(
      "Saved primary k=", selected_k,
      " solution used to define the joint-fit warm starts"
    ),
    "Candidate primary partition"
  )

  selected <- f6r_read_tsv(file.path(
    paths$figure5, "figure5_frozen_inputs", "selected_results.tsv"
  ))
  selected <- selected[selected$record_type == "joint_pair_best", , drop = FALSE]
  if (nrow(selected) != selected_k ||
      !identical(sort(selected$family), expected_regions)) {
    stop("Expected one selected joint fit for each declared primary region.")
  }
  representatives <- stats::setNames(
    as.integer(sub("^seed", "", selected$invivo_seed)), selected$family
  )
  assignments <- data.frame(
    seed = endpoints$seed,
    primary_region = endpoints$cluster_base_id,
    subcluster = "",
    warm_start_stratum = endpoints$cluster_base_id,
    objective = endpoints$objective,
    representative = endpoints$seed == unname(
      representatives[endpoints$cluster_base_id]
    ),
    stringsAsFactors = FALSE
  )
  representative_summary <- assignments[assignments$representative, , drop = FALSE]
  representative_summary <- representative_summary[
    order(representative_summary$primary_region), , drop = FALSE
  ]
  if (nrow(representative_summary) != selected_k ||
      !identical(representative_summary$primary_region, expected_regions)) {
    stop("Could not recover all primary-region warm-start representatives.")
  }
  representative_summary$selected_k <- selected_k
  representative_summary$selected_average_silhouette <- saved_silhouette

  repeated_primary <- vapply(seq_len(n_resample), function(run) {
    set.seed(71000L + run)
    labels <- stats::kmeans(
      coords, centers = selected_k, nstart = 1L, iter.max = 100L
    )$cluster
    f6r_adjusted_rand(saved_labels, labels)
  }, numeric(1L))

  bootstrap_primary <- matrix(
    NA_real_, nrow = n_resample, ncol = 2L,
    dimnames = list(NULL, c("ari", "selected_k"))
  )
  for (run in seq_len(n_resample)) {
    set.seed(72000L + run)
    idx <- sort(sample(seq_len(nrow(coords)), size = floor(0.8 * nrow(coords))))
    x <- coords[idx, , drop = FALSE]
    d <- stats::dist(x)
    set.seed(73000L + run)
    km6 <- stats::kmeans(x, centers = selected_k, nstart = 25L, iter.max = 200L)
    bootstrap_primary[run, "ari"] <- f6r_adjusted_rand(
      saved_labels[idx], km6$cluster
    )
    sil <- vapply(2:8, function(k) {
      set.seed(74000L + 100L * run + k)
      km <- stats::kmeans(x, centers = k, nstart = 25L, iter.max = 200L)
      f6r_silhouette(km$cluster, d)
    }, numeric(1L))
    bootstrap_primary[run, "selected_k"] <- (2:8)[which.max(sil)]
  }

  parameter_wide <- f6r_endpoint_parameter_matrix(paths)
  joined <- merge(
    endpoints, parameter_wide,
    by.x = "seed", by.y = "seed_number", sort = TRUE
  )
  global_parameter_x <- scale(as.matrix(
    joined[, f6r_shared_parameters(), drop = FALSE]
  ))
  set.seed(75001L)
  parameter_k6 <- stats::kmeans(
    global_parameter_x, centers = selected_k, nstart = 100L, iter.max = 200L
  )$cluster
  parameter_space_ari <- f6r_adjusted_rand(
    match(joined$cluster_base_id, expected_regions), parameter_k6
  )

  stability <- do.call(rbind, list(
    data.frame(
      analysis_level = "primary pooled t-SNE regions",
      region = "all in-vivo endpoints",
      perturbation = paste0(
        n_resample, " one-start k-means initializations at k=", selected_k
      ),
      n_runs = n_resample,
      ari_median = stats::median(repeated_primary),
      ari_q05 = unname(stats::quantile(repeated_primary, 0.05)),
      ari_q95 = unname(stats::quantile(repeated_primary, 0.95)),
      proportion_ari_ge_0p9 = mean(repeated_primary >= 0.9),
      stringsAsFactors = FALSE
    ),
    data.frame(
      analysis_level = "primary pooled t-SNE regions",
      region = "all in-vivo endpoints",
      perturbation = paste0(
        n_resample, " 80%-endpoint subsamples without replacement at k=",
        selected_k
      ),
      n_runs = n_resample,
      ari_median = stats::median(bootstrap_primary[, "ari"]),
      ari_q05 = unname(stats::quantile(bootstrap_primary[, "ari"], 0.05)),
      ari_q95 = unname(stats::quantile(bootstrap_primary[, "ari"], 0.95)),
      proportion_ari_ge_0p9 = mean(bootstrap_primary[, "ari"] >= 0.9),
      stringsAsFactors = FALSE
    ),
    data.frame(
      analysis_level = "primary pooled t-SNE regions",
      region = "all in-vivo endpoints",
      perturbation = paste0(
        "k=", selected_k,
        " clustering in standardized 14-parameter endpoint space"
      ),
      n_runs = 1L,
      ari_median = parameter_space_ari,
      ari_q05 = parameter_space_ari,
      ari_q95 = parameter_space_ari,
      proportion_ari_ge_0p9 = as.numeric(parameter_space_ari >= 0.9),
      stringsAsFactors = FALSE
    )
  ))
  bootstrap_frequency <- as.data.frame(table(
    selected_k = as.integer(bootstrap_primary[, "selected_k"])
  ))
  names(bootstrap_frequency)[[2L]] <- "n_subsample"
  bootstrap_frequency$proportion <-
    bootstrap_frequency$n_subsample / n_resample

  robust_selected_k <- k_selection$average_silhouette[
    k_selection$k == selected_k
  ]
  paths_out <- c(
    k_selection = f6r_write_tsv(
      k_selection, file.path(paths$figure6, "cluster_k_selection.tsv")
    ),
    stability = f6r_write_tsv(
      stability, file.path(paths$figure6, "cluster_stability.tsv")
    ),
    resample_frequency = f6r_write_tsv(
      bootstrap_frequency,
      file.path(paths$figure6, "cluster_k_resample_selection_frequency.tsv")
    ),
    bootstrap_frequency_legacy_alias = f6r_write_tsv(
      bootstrap_frequency,
      file.path(paths$figure6, "cluster_k_bootstrap_selection_frequency.tsv")
    ),
    assignments = f6r_write_tsv(
      assignments, file.path(paths$figure6, "cluster_subcluster_assignments.tsv")
    ),
    representatives = f6r_write_tsv(
      representative_summary,
      file.path(paths$figure6, "cluster_warm_start_representatives.tsv")
    ),
    audit = f6r_write_tsv(
      data.frame(
        metric = c(
          paste0("saved_primary_k", selected_k, "_average_silhouette"),
          paste0("robust_k", selected_k, "_average_silhouette"),
          paste0("saved_minus_robust_k", selected_k, "_silhouette"),
          paste0("primary_parameter_space_k", selected_k, "_ARI"),
          "tSNE_seed_perplexity_sensitivity"
        ),
        value = c(
          saved_silhouette, robust_selected_k,
          saved_silhouette - robust_selected_k, parameter_space_ari, NA
        ),
        interpretation = c(
          paste0("Silhouette of the saved k=", selected_k, " assignment"),
          paste0("High-restart audit of k=", selected_k),
          paste0(
            "Difference between saved and independently reconstructed k=",
            selected_k, " partitions"
          ),
          paste0(
            "Agreement between saved t-SNE regions and endpoint-space k=",
            selected_k
          ),
          paste0(
            "Not estimated from the frozen embedding because the source ",
            "parameter vectors required to re-estimate t-SNE are not vendored"
          )
        ),
        stringsAsFactors = FALSE
      ),
      file.path(paths$figure6, "cluster_selection_audit.tsv")
    )
  )
  invisible(list(
    k_selection = k_selection,
    stability = stability,
    assignments = assignments,
    representatives = representative_summary,
    bootstrap_frequency = bootstrap_frequency,
    paths = paths_out
  ))
}

f6r_selected_results <- function(paths) {
  selected <- f6r_read_tsv(file.path(
    paths$figure5, "figure5_frozen_inputs", "selected_results.tsv"
  ))
  selected <- selected[selected$record_type == "joint_pair_best", , drop = FALSE]
  family_levels <- f6r_family_levels()
  if (nrow(selected) != length(family_levels)) {
    stop("Expected one selected joint pair per primary family.")
  }
  selected$pair_id <- paste0("fit_joint_", selected$warmup_label)
  selected$selected_seed_number <- as.integer(sub("^seed", "", selected$selected_seed))
  selected$pair_label <- sub(
    ".*_(C[0-9]{2})_vt.*", "\\1", selected$warmup_label
  )
  if (!identical(selected$pair_label, family_levels)) {
    stop("Selected Figure 6 fits must resolve to the declared primary regions.")
  }
  selected$bundle_path <- file.path(
    paths$figure5, "figure5_frozen_inputs", "winners", selected$warmup_label
  )
  selected <- selected[order(selected$pair_label), , drop = FALSE]
  rownames(selected) <- NULL
  selected
}

f6r_objective_selection <- function(paths) {
  f6r_require_packages("data.table")
  selected <- f6r_selected_results(paths)
  family_levels <- f6r_family_levels()
  objective_rows <- lapply(seq_len(nrow(selected)), function(i) {
    path <- file.path(
      paths$supp5_2, "joint", selected$warmup_label[[i]],
      "seed_objective_simple.tsv"
    )
    x <- f6r_read_tsv(path)
    x$pair_id <- selected$pair_id[[i]]
    x$pair_label <- selected$pair_label[[i]]
    x$seed_number <- as.integer(sub("^seed", "", x$seed))
    x$source_objective_file <- normalizePath(path, mustWork = TRUE)
    x
  })
  objectives <- do.call(rbind, objective_rows)
  rownames(objectives) <- NULL
  expected_joint_seed_count <- 500L * length(family_levels)
  if (nrow(objectives) != expected_joint_seed_count ||
      any(table(objectives$pair_id) != 500L) ||
      anyDuplicated(objectives[, c("pair_id", "seed_number")])) {
    stop(
      "Joint objective cache must contain ", length(family_levels),
      " groups of 500 unique seeds."
    )
  }

  master_path <- file.path(paths$supp5_1, "soft_coupling_master_long.tsv")
  master <- data.table::as.data.table(f6r_read_tsv(master_path))
  parameter_qc <- master[, .(
    n_parameter_record = .N,
    n_parameter = data.table::uniqueN(parameter),
    all_parameter_finite = all(is.finite(vivo_natural)),
    all_feasible_at_solution = all(as.logical(feasible_at_solution)),
    all_feasible_before_projection = all(as.logical(feasible_before_projection)),
    any_projection_applied = any(as.logical(projection_applied))
  ), by = .(pair_id, seed_number)]
  if (nrow(parameter_qc) != expected_joint_seed_count ||
      any(parameter_qc$n_parameter_record != length(f6r_shared_parameters())) ||
      any(parameter_qc$n_parameter != length(f6r_shared_parameters()))) {
    stop(
      "The local joint parameter cache is not complete for all ",
      expected_joint_seed_count, " endpoints."
    )
  }
  endpoint_signatures <- f6r_endpoint_signatures(master)
  invitro_endpoint_signatures <- f6r_endpoint_signatures(
    master,
    value_column = "vitro_natural",
    signature_column = "parameter_signature_invitro"
  )
  parameter_qc <- merge(
    as.data.frame(parameter_qc), endpoint_signatures,
    by = c("pair_id", "seed_number"), all = TRUE, sort = FALSE
  )
  parameter_qc <- merge(
    parameter_qc, invitro_endpoint_signatures,
    by = c("pair_id", "seed_number"), all = TRUE, sort = FALSE
  )
  objectives <- merge(
    objectives, parameter_qc,
    by = c("pair_id", "seed_number"), all.x = TRUE, sort = FALSE
  )
  objectives <- objectives[order(
    objectives$pair_label, objectives$objective_rank
  ), , drop = FALSE]

  config_map <- stats::setNames(
    file.path(selected$bundle_path, "fit_config.rds"), selected$pair_id
  )
  objectives$config_available <- file.exists(config_map[objectives$pair_id])
  objectives$objective_finite <- is.finite(objectives$objective)
  objectives$seed_id_valid <- objectives$seed_number %in% 1:500
  objectives$hard_qc_pass <- with(
    objectives,
    objective_finite & seed_id_valid & n_parameter_record == 14L &
      n_parameter == 14L &
      all_parameter_finite & all_feasible_at_solution &
      all_feasible_before_projection & !any_projection_applied &
      config_available
  )

  objectives$delta_objective <- NA_real_
  objectives$within_pair_empirical_quantile <- NA_real_
  objectives$eligible_q05 <- FALSE
  objectives$eligible_q10 <- FALSE
  objectives$eligible_q20 <- FALSE
  for (pair in unique(objectives$pair_id)) {
    idx <- which(objectives$pair_id == pair & objectives$hard_qc_pass)
    idx <- idx[order(objectives$objective[idx], objectives$seed_number[idx])]
    n <- length(idx)
    objectives$delta_objective[idx] <-
      objectives$objective[idx] - min(objectives$objective[idx])
    objectives$within_pair_empirical_quantile[idx] <- seq_len(n) / n
    objectives$eligible_q05[idx[seq_len(ceiling(0.05 * n))]] <- TRUE
    objectives$eligible_q10[idx[seq_len(ceiling(0.10 * n))]] <- TRUE
    objectives$eligible_q20[idx[seq_len(ceiling(0.20 * n))]] <- TRUE
  }
  objectives$primary_acceptable <- objectives$eligible_q10

  objectives$parameter_endpoint_group <- NA_character_
  objectives$parameter_endpoint_group_invitro <- NA_character_
  for (pair in unique(objectives$pair_id)) {
    idx <- which(objectives$pair_id == pair)
    signatures_in_rank_order <- unique(objectives$parameter_signature[idx])
    objectives$parameter_endpoint_group[idx] <- paste0(
      objectives$pair_label[idx], "_E",
      sprintf("%03d", match(objectives$parameter_signature[idx], signatures_in_rank_order))
    )
    invitro_signatures_in_rank_order <- unique(
      objectives$parameter_signature_invitro[idx]
    )
    objectives$parameter_endpoint_group_invitro[idx] <- paste0(
      objectives$pair_label[idx], "_VT_E",
      sprintf(
        "%03d",
        match(
          objectives$parameter_signature_invitro[idx],
          invitro_signatures_in_rank_order
        )
      )
    )
  }
  objectives$endpoint_multiplicity_all500 <- ave(
    objectives$seed_number,
    objectives$pair_id, objectives$parameter_endpoint_group,
    FUN = length
  )
  for (cutoff_name in c("q05", "q10", "q20")) {
    eligible_field <- paste0("eligible_", cutoff_name)
    multiplicity_field <- paste0("endpoint_multiplicity_", cutoff_name)
    objectives[[multiplicity_field]] <- 0L
    idx <- which(objectives[[eligible_field]])
    objectives[[multiplicity_field]][idx] <- ave(
      objectives$seed_number[idx],
      objectives$pair_id[idx], objectives$parameter_endpoint_group[idx],
      FUN = length
    )
  }
  objectives$operator_qc_pass <- NA
  objectives$final_primary_acceptable <- NA
  objectives$exclusion_reason <- ifelse(
    !objectives$hard_qc_pass, "failed local objective/parameter/config QC",
    ifelse(
      objectives$eligible_q10, "included: pair-specific lowest 10% full MAP objective",
      "excluded from primary ensemble: objective rank above pair-specific lowest 10%"
    )
  )

  cutoff_definitions <- c(q05 = 0.05, q10 = 0.10, q20 = 0.20)
  cutoff_rows <- do.call(rbind, lapply(unique(objectives$pair_id), function(pair) {
    z <- objectives[objectives$pair_id == pair & objectives$hard_qc_pass, , drop = FALSE]
    z <- z[order(z$objective, z$seed_number), , drop = FALSE]
    do.call(rbind, lapply(names(cutoff_definitions), function(cutoff_name) {
      q <- unname(cutoff_definitions[[cutoff_name]])
      n <- ceiling(q * nrow(z))
      data.frame(
        pair_id = pair,
        pair_label = z$pair_label[[1L]],
        cutoff = cutoff_name,
        retained_fraction = q,
        n_hard_qc = nrow(z),
        n_accepted = n,
        objective_cutoff = z$objective[[n]],
        delta_objective_cutoff = z$delta_objective[[n]],
        stringsAsFactors = FALSE
      )
    }))
  }))

  endpoint_multiplicity <- objectives[, c(
    "pair_id", "pair_label", "seed", "seed_number", "objective_rank",
    "objective", "parameter_signature", "parameter_endpoint_group",
    "endpoint_multiplicity_all500", "endpoint_multiplicity_q05",
    "endpoint_multiplicity_q10", "endpoint_multiplicity_q20",
    "eligible_q05", "eligible_q10", "eligible_q20"
  )]
  endpoint_counts <- do.call(rbind, lapply(
    unique(objectives$pair_id),
    function(pair) {
      z <- objectives[objectives$pair_id == pair, , drop = FALSE]
      do.call(rbind, lapply(c(q05 = "eligible_q05", q10 = "eligible_q10", q20 = "eligible_q20"), function(field) {
        keep <- z[[field]]
        data.frame(
          pair_id = pair,
          pair_label = z$pair_label[[1L]],
          cutoff = sub("eligible_", "", field),
          n_seed = sum(keep),
          n_unique_parameter_endpoint = length(unique(z$parameter_endpoint_group[keep])),
          maximum_endpoint_multiplicity = max(table(z$parameter_endpoint_group[keep])),
          stringsAsFactors = FALSE
        )
      }))
    }
  ))
  rownames(endpoint_counts) <- NULL

  selected_master <- as.data.frame(master[
    objectives[objectives$eligible_q20, c("pair_id", "seed_number")],
    on = .(pair_id, seed_number), nomatch = 0L
  ])
  parameter_long <- selected_master[, c(
    "pair_id", "seed_number", "parameter", "vivo_natural"
  )]
  names(parameter_long)[names(parameter_long) == "vivo_natural"] <- "value"
  parameter_long$value_role <- "joint in-vivo effective value"
  rho_rows <- lapply(seq_len(nrow(selected)), function(i) {
    bp <- f6r_read_tsv(file.path(selected$bundle_path[[i]], "best_params.tsv"))
    rho <- bp$value[bp$parameter == "rho_2N"]
    seeds <- objectives$seed_number[
      objectives$pair_id == selected$pair_id[[i]] & objectives$eligible_q20
    ]
    data.frame(
      pair_id = selected$pair_id[[i]],
      seed_number = seeds,
      parameter = "rho_2N",
      value = as.numeric(rho[[1L]]),
      value_role = paste0(
        "pair-reference nuisance value; required by the generic run-parameter ",
        "schema but absent from the fixed-O2 transition operator"
      ),
      stringsAsFactors = FALSE
    )
  })
  parameter_long <- rbind(parameter_long, do.call(rbind, rho_rows))
  parameter_long <- parameter_long[order(
    parameter_long$pair_id, parameter_long$seed_number, parameter_long$parameter
  ), , drop = FALSE]
  parameter_long_invitro <- selected_master[, c(
    "pair_id", "seed_number", "parameter", "vitro_natural"
  )]
  names(parameter_long_invitro)[
    names(parameter_long_invitro) == "vitro_natural"
  ] <- "value"
  parameter_long_invitro$value_role <- "joint in-vitro effective value"
  parameter_long_invitro <- parameter_long_invitro[order(
    parameter_long_invitro$pair_id,
    parameter_long_invitro$seed_number,
    parameter_long_invitro$parameter
  ), , drop = FALSE]

  output_paths <- c(
    fit_qc = f6r_write_tsv(
      objectives[, c(
        "pair_id", "pair_label", "seed", "seed_number", "objective_rank",
        "objective", "objective_finite", "seed_id_valid",
        "n_parameter_record", "n_parameter",
        "all_parameter_finite", "all_feasible_at_solution",
        "all_feasible_before_projection", "any_projection_applied",
        "config_available", "hard_qc_pass", "operator_qc_pass"
      )],
      file.path(paths$figure6, "joint_seed_fit_qc.tsv")
    ),
    acceptance = f6r_write_tsv(
      objectives,
      file.path(paths$figure6, "joint_seed_acceptance.tsv")
    ),
    cutoff_summary = f6r_write_tsv(
      cutoff_rows,
      file.path(paths$figure6, "joint_seed_acceptance_summary.tsv")
    ),
    parameters = f6r_write_tsv(
      parameter_long,
      file.path(paths$figure6, "joint_acceptable_seed_invivo_parameters.tsv")
    ),
    parameters_invitro = f6r_write_tsv(
      parameter_long_invitro,
      file.path(paths$figure6, "joint_acceptable_seed_invitro_parameters.tsv")
    ),
    endpoint_multiplicity = f6r_write_tsv(
      endpoint_multiplicity,
      file.path(paths$figure6, "joint_seed_parameter_multiplicity.tsv")
    ),
    endpoint_counts = f6r_write_tsv(
      endpoint_counts,
      file.path(paths$figure6, "joint_unique_parameter_endpoint_counts.tsv")
    )
  )
  invisible(list(
    objectives = objectives,
    cutoffs = cutoff_rows,
    master = as.data.frame(master),
    selected = selected,
    parameters = parameter_long,
    parameters_invitro = parameter_long_invitro,
    endpoint_multiplicity = endpoint_multiplicity,
    endpoint_counts = endpoint_counts,
    paths = output_paths,
    master_path = master_path
  ))
}

f6r_load_response_engine <- local({
  loaded <- FALSE
  function(paths) {
    if (loaded && exists(
      "cpp_o2simps_build_G_for_o2_triplet",
      envir = globalenv(), inherits = FALSE
    )) {
      return(invisible(TRUE))
    }
    engine <- file.path(paths$code, "util", "oxygen", "response_pipeline.R")
    current_model_engine <- file.path(
      paths$oxygen_code, "simulation", "o2", "fixed_o2",
      "run_fixed_o2_simulation.R"
    )
    f6r_require_files(
      c(engine, current_model_engine),
      "analytical response engine and current-branch model"
    )
    sys.source(engine, envir = globalenv())
    dominant_attractor_with_population <- get(
      "fixo2_dominant_attractor_one", envir = globalenv(), inherits = FALSE
    )
    # Override the packaged analytical helpers with the required external
    # model implementation. Figure-specific analysis and plotting remain in
    # revised/iteration4, while the model source is read-only.
    source(current_model_engine, local = globalenv(), chdir = TRUE)
    # The latest canonical simulator omits the population-average
    # missegregation column from this otherwise identical attractor producer.
    # Retain the figure-workspace producer so the latest model implementation
    # is evaluated while preserving the established Figure 6 table contract.
    assign(
      "fixo2_dominant_attractor_one",
      dominant_attractor_with_population,
      envir = globalenv()
    )
    loaded <<- TRUE
    invisible(TRUE)
  }
})

f6r_pair_model_context <- function(pair_id, selected, paths) {
  row <- selected[selected$pair_id == pair_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Cannot resolve selected pair metadata for ", pair_id)
  config_path <- file.path(row$bundle_path, "fit_config.rds")
  best_path <- file.path(row$bundle_path, "best_params.tsv")
  f6r_require_files(c(config_path, best_path), "pair model context")
  best <- f6r_read_tsv(best_path)
  rho <- as.numeric(best$value[best$parameter == "rho_2N"])
  if (length(rho) != 1L || !is.finite(rho) || rho <= 0) {
    stop("Cannot recover pair-reference rho_2N for ", pair_id)
  }
  list(
    pair_id = pair_id,
    pair_label = row$pair_label[[1L]],
    config = readRDS(config_path),
    config_path = normalizePath(config_path, mustWork = TRUE),
    best_path = normalizePath(best_path, mustWork = TRUE),
    rho_2N = rho,
    selected_seed = row$selected_seed_number[[1L]],
    warmup_label = row$warmup_label[[1L]]
  )
}

f6r_compute_seed_cache <- function(
    pair_id, seed_number, objective, master, context, cache_path,
    parameter_source, full_surface = TRUE, force_rebuild = FALSE,
    model_context = "in vivo", parameter_value_column = "vivo_natural",
    simulation_mode = "joint", model_source_fingerprint
) {
  requested_surface_profile <- if (isTRUE(full_surface)) {
    "full_201x60"
  } else {
    "claim_diagnostic_union"
  }
  existing <- NULL
  if (file.exists(cache_path) && !isTRUE(force_rebuild)) {
    existing <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    expected_surface_rows <- if (isTRUE(full_surface)) 201L * 60L else 774L
    existing_context <- if (!is.null(existing$metadata$model_context)) {
      as.character(existing$metadata$model_context)
    } else {
      "in vivo"
    }
    if (!is.null(existing) &&
        identical(as.character(existing$metadata$pair_id), as.character(pair_id)) &&
        identical(existing_context, as.character(model_context)) &&
        identical(
          as.character(existing$metadata$model_source_fingerprint),
          as.character(model_source_fingerprint)
        ) &&
        identical(as.integer(existing$metadata$seed_number), as.integer(seed_number)) &&
        identical(
          as.character(existing$metadata$surface_profile),
          requested_surface_profile
        ) &&
        isTRUE(all.equal(
          as.numeric(existing$metadata$objective), as.numeric(objective), tolerance = 1e-12
        )) &&
        nrow(existing$trajectory) == 201L &&
        nrow(existing$surface) == expected_surface_rows &&
        isTRUE(existing$qc$operator_qc_pass[[1L]])) {
      return(existing$qc)
    }
  }

  rows <- master[
    master$pair_id == pair_id & master$seed_number == seed_number,
    , drop = FALSE
  ]
  shared <- f6r_shared_parameters()
  if (nrow(rows) != length(shared) || !setequal(rows$parameter, shared)) {
    stop("Incomplete parameters for ", pair_id, " seed", seed_number)
  }
  if (!parameter_value_column %in% names(rows)) {
    stop("Missing context parameter column: ", parameter_value_column)
  }
  params <- stats::setNames(
    as.numeric(rows[[parameter_value_column]]), rows$parameter
  )
  if (identical(simulation_mode, "joint")) {
    params[["rho_2N"]] <- context$rho_2N
  }
  run_params <- prepare_run_params(
    param_values = params,
    simulation = simulation_mode,
    cfg = context$config,
    fixed_o2 = 5
  )
  config <- prepare_sim_cfg(
    context$config, argv = list(), fixed_o2 = 5, run_params = run_params
  )
  run_params$O2_growth <- isTRUE(config$O2_growth)
  run_params$ploidy_O2_death <- config$ploidy_O2_death

  o2_values <- seq(0, 5, length.out = 201L)
  cin_values <- 10^seq(log10(0.005), log10(0.5), length.out = 60L)
  seed_id <- paste0("seed", seed_number)

  trajectory_rows <- lapply(o2_values, function(o2) {
    fixo2_dominant_attractor_one(
      seed_id = seed_id,
      run_params = run_params,
      model_env = globalenv(),
      cfg = config,
      O2 = o2
    )
  })
  trajectory <- do.call(rbind, trajectory_rows)
  trajectory <- trajectory[, c(
    "O2_pct", "status", "population_average_p_misseg",
    "dominant_mean_ploidy", "spectral_gap", "dominant_growth_rate"
  ), drop = FALSE]

  if (isTRUE(full_surface)) {
    surface_grid <- expand.grid(
      O2_pct = o2_values,
      effective_p_misseg = cin_values,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    cin_targets <- unique(vapply(c(0.005, 0.05, 0.5), function(target) {
      cin_values[[which.min(abs(cin_values - target))]]
    }, numeric(1L)))
    surface_grid <- unique(rbind(
      expand.grid(
        O2_pct = o2_values,
        effective_p_misseg = cin_targets,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      ),
      expand.grid(
        O2_pct = c(0, 1, 5),
        effective_p_misseg = cin_values,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    ))
  }
  surface_rows <- lapply(
    sort(unique(surface_grid$effective_p_misseg)),
    function(cin) {
    o2_subset <- surface_grid$O2_pct[
      abs(surface_grid$effective_p_misseg - cin) < 1e-12
    ]
    forced <- response_force_effective_p_misseg(run_params, cin)
    rows <- lapply(o2_subset, function(o2) {
      z <- fixo2_dominant_attractor_one(
        seed_id = seed_id,
        run_params = forced,
        model_env = globalenv(),
        cfg = config,
        O2 = o2
      )
      data.frame(
        O2_pct = as.numeric(z$O2_pct[[1L]]),
        effective_p_misseg = cin,
        status = as.character(z$status[[1L]]),
        actual_effective_p_misseg =
          as.numeric(z$population_average_p_misseg[[1L]]),
        dominant_mean_ploidy = as.numeric(z$dominant_mean_ploidy[[1L]]),
        spectral_gap = as.numeric(z$spectral_gap[[1L]]),
        dominant_growth_rate = as.numeric(z$dominant_growth_rate[[1L]]),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  surface <- do.call(rbind, surface_rows)
  surface <- surface[order(
    surface$effective_p_misseg, surface$O2_pct
  ), , drop = FALSE]
  trajectory$pair_id <- pair_id
  trajectory$pair_label <- context$pair_label
  trajectory$seed_number <- seed_number
  trajectory$objective <- objective
  trajectory$model_context <- model_context
  surface$pair_id <- pair_id
  surface$pair_label <- context$pair_label
  surface$seed_number <- seed_number
  surface$objective <- objective
  surface$model_context <- model_context
  surface$surface_profile <- requested_surface_profile

  valid_trajectory <- all(trajectory$status == "ok") && all(is.finite(
    trajectory$population_average_p_misseg
  )) && all(is.finite(trajectory$dominant_mean_ploidy)) &&
    all(is.finite(trajectory$spectral_gap))
  valid_surface <- all(surface$status == "ok") && all(is.finite(
    surface$dominant_mean_ploidy
  )) && all(is.finite(surface$spectral_gap))
  max_error <- max(abs(
    surface$actual_effective_p_misseg - surface$effective_p_misseg
  ), na.rm = TRUE)
  qc <- data.frame(
    pair_id = pair_id,
    pair_label = context$pair_label,
    model_context = model_context,
    seed_number = seed_number,
    objective = objective,
    surface_profile = requested_surface_profile,
    n_trajectory = nrow(trajectory),
    n_surface = nrow(surface),
    trajectory_all_status_ok = valid_trajectory,
    surface_all_status_ok = valid_surface,
    max_abs_actual_minus_requested_p_misseg = max_error,
    operator_qc_pass = valid_trajectory && valid_surface &&
      is.finite(max_error) && max_error <= 1e-8,
    cache_path = normalizePath(cache_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(cache_path, ".tmp-", Sys.getpid())
  saveRDS(
    list(
      metadata = list(
        pair_id = pair_id, pair_label = context$pair_label,
        model_context = model_context,
        seed_number = seed_number, objective = objective,
        config_path = context$config_path,
        surface_profile = requested_surface_profile,
        model_source_fingerprint = model_source_fingerprint,
        parameter_source = normalizePath(parameter_source, mustWork = TRUE)
      ),
      trajectory = trajectory,
      surface = surface,
      qc = qc
    ),
    temporary,
    compress = "gzip"
  )
  if (!file.copy(temporary, cache_path, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Failed to publish seed cache: ", cache_path)
  }
  unlink(temporary)
  qc
}

f6r_figure6d_p_values <- function() {
  seq.int(5L, 500L, by = 1L) / 1000
}

f6r_figure6d_o2_values <- function() {
  seq(0, 5, length.out = 201L)
}

f6r_figure6d_endpoint_manifest <- function(paths, objective_bundle) {
  acceptance <- objective_bundle$objectives
  display_manifest <- f6r_display_pair_manifest(acceptance$pair_id, "D")
  selected <- acceptance[
    acceptance$pair_id %in% display_manifest$pair_id &
      acceptance$eligible_q10 & acceptance$hard_qc_pass,
    , drop = FALSE
  ]
  selected$display_label <- display_manifest$display_label[
    match(selected$pair_id, display_manifest$pair_id)
  ]
  observed_seed_counts <- table(factor(
    selected$pair_label, levels = display_manifest$pair_label
  ))
  if (nrow(selected) != 50L * f6r_family_count() ||
      any(observed_seed_counts != 50L)) {
    stop("Figure 6D requires exactly 50 q10 endpoints per displayed pair.")
  }

  split_groups <- split(
    selected,
    interaction(
      selected$pair_id, selected$parameter_endpoint_group,
      drop = TRUE, lex.order = TRUE
    )
  )
  endpoint_rows <- lapply(split_groups, function(z) {
    z <- z[order(z$objective_rank, z$seed_number), , drop = FALSE]
    data.frame(
      display_label = z$display_label[[1L]],
      pair_label = z$pair_label[[1L]],
      pair_id = z$pair_id[[1L]],
      parameter_endpoint_group = z$parameter_endpoint_group[[1L]],
      representative_seed_number = z$seed_number[[1L]],
      representative_objective_rank = z$objective_rank[[1L]],
      representative_objective = z$objective[[1L]],
      endpoint_multiplicity_q10 = nrow(z),
      represented_seed_numbers = paste(z$seed_number, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  endpoints <- do.call(rbind, endpoint_rows)
  endpoints <- endpoints[order(
    match(endpoints$display_label, display_manifest$display_label),
    endpoints$representative_objective_rank,
    endpoints$representative_seed_number
  ), , drop = FALSE]
  rownames(endpoints) <- NULL
  unique_counts <- table(factor(
    endpoints$pair_label, levels = display_manifest$pair_label
  ))
  multiplicity_counts <- tapply(
    endpoints$endpoint_multiplicity_q10,
    factor(endpoints$pair_label, levels = display_manifest$pair_label),
    sum
  )
  if (length(unique_counts) != f6r_family_count() || any(unique_counts < 1L) ||
      any(multiplicity_counts != 50L)) {
    stop(
      "Unexpected Figure 6D endpoint multiplicity: unique counts ",
      paste(as.integer(unique_counts), collapse = ","),
      "; represented counts ",
      paste(as.integer(multiplicity_counts), collapse = ",")
    )
  }
  manifest_path <- f6r_write_tsv(
    endpoints,
    file.path(paths$figure6, "figure6d_dense_endpoint_manifest.tsv")
  )
  invisible(list(
    endpoints = endpoints,
    display_manifest = display_manifest,
    path = manifest_path
  ))
}

f6r_figure6d_compute_endpoint_cache <- function(
    metadata, parameters, context, cache_path, parameter_source,
    force_rebuild = FALSE, model_context = "in vivo",
    simulation_mode = "joint", model_source_fingerprint
) {
  p_values <- f6r_figure6d_p_values()
  o2_values <- f6r_figure6d_o2_values()
  expected_rows <- length(p_values) * length(o2_values)
  profile <- "dense_201x496_step0p001_v1"
  if (file.exists(cache_path) && !isTRUE(force_rebuild)) {
    existing <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    existing_context <- if (!is.null(existing$metadata$model_context)) {
      as.character(existing$metadata$model_context)
    } else {
      "in vivo"
    }
    if (!is.null(existing) &&
        identical(as.character(existing$metadata$profile), profile) &&
        identical(existing_context, as.character(model_context)) &&
        identical(
          as.character(existing$metadata$model_source_fingerprint),
          as.character(model_source_fingerprint)
        ) &&
        identical(
          as.character(existing$metadata$pair_id),
          as.character(metadata$pair_id)
        ) &&
        identical(
          as.character(existing$metadata$parameter_endpoint_group),
          as.character(metadata$parameter_endpoint_group)
        ) &&
        identical(
          as.integer(existing$metadata$representative_seed_number),
          as.integer(metadata$representative_seed_number)
        ) &&
        nrow(existing$surface) == expected_rows &&
        isTRUE(existing$qc$operator_qc_pass[[1L]])) {
      return(existing$qc)
    }
  }

  seed_number <- metadata$representative_seed_number[[1L]]
  z <- parameters[
    parameters$pair_id == metadata$pair_id[[1L]] &
      parameters$seed_number == seed_number,
    , drop = FALSE
  ]
  required_parameters <- if (identical(simulation_mode, "joint")) {
    c(f6r_shared_parameters(), "rho_2N")
  } else {
    f6r_shared_parameters()
  }
  if (nrow(z) != length(required_parameters) ||
      !setequal(z$parameter, required_parameters) ||
      any(!is.finite(z$value))) {
    stop(
      "Incomplete Figure 6D parameters for ", metadata$pair_label,
      " representative seed ", seed_number
    )
  }
  params <- stats::setNames(as.numeric(z$value), z$parameter)
  run_params <- prepare_run_params(
    param_values = params,
    simulation = simulation_mode,
    cfg = context$config,
    fixed_o2 = 5
  )
  config <- prepare_sim_cfg(
    context$config, argv = list(), fixed_o2 = 5, run_params = run_params
  )
  run_params$O2_growth <- isTRUE(config$O2_growth)
  run_params$ploidy_O2_death <- config$ploidy_O2_death
  seed_id <- paste0("seed", seed_number)

  # Preallocate the dense table.  The numerical operator and evaluation order
  # are unchanged; avoiding one single-row data.frame allocation per grid cell
  # materially reduces R bookkeeping for the 99,696-condition endpoint grid.
  n_o2 <- length(o2_values)
  n_surface <- length(p_values) * n_o2
  surface_o2 <- rep(o2_values, times = length(p_values))
  surface_p <- rep(p_values, each = n_o2)
  surface_status <- character(n_surface)
  surface_actual_p <- numeric(n_surface)
  surface_ploidy <- numeric(n_surface)
  surface_gap <- numeric(n_surface)
  surface_growth <- numeric(n_surface)
  for (p_index in seq_along(p_values)) {
    p_fixed <- p_values[[p_index]]
    forced <- response_force_effective_p_misseg(run_params, p_fixed)
    offset <- (p_index - 1L) * n_o2
    for (o2_index in seq_along(o2_values)) {
      row_index <- offset + o2_index
      result <- fixo2_dominant_attractor_one(
        seed_id = seed_id,
        run_params = forced,
        model_env = globalenv(),
        cfg = config,
        O2 = o2_values[[o2_index]]
      )
      surface_o2[[row_index]] <- as.numeric(result$O2_pct[[1L]])
      surface_status[[row_index]] <- as.character(result$status[[1L]])
      surface_actual_p[[row_index]] <-
        as.numeric(result$population_average_p_misseg[[1L]])
      surface_ploidy[[row_index]] <-
        as.numeric(result$dominant_mean_ploidy[[1L]])
      surface_gap[[row_index]] <- as.numeric(result$spectral_gap[[1L]])
      surface_growth[[row_index]] <-
        as.numeric(result$dominant_growth_rate[[1L]])
    }
  }
  surface <- data.frame(
    O2_pct = surface_o2,
    effective_p_misseg = surface_p,
    status = surface_status,
    actual_effective_p_misseg = surface_actual_p,
    dominant_mean_ploidy = surface_ploidy,
    spectral_gap = surface_gap,
    dominant_growth_rate = surface_growth,
    stringsAsFactors = FALSE
  )
  surface <- surface[order(
    surface$effective_p_misseg, surface$O2_pct
  ), , drop = FALSE]
  surface$model_context <- model_context
  forcing_error <- abs(
    surface$actual_effective_p_misseg - surface$effective_p_misseg
  )
  finite_forcing_error <- forcing_error[is.finite(forcing_error)]
  max_error <- if (length(finite_forcing_error)) {
    max(finite_forcing_error)
  } else {
    NA_real_
  }
  valid <- nrow(surface) == expected_rows &&
    all(surface$status == "ok") &&
    all(is.finite(surface$dominant_mean_ploidy)) &&
    all(is.finite(surface$spectral_gap)) &&
    all(is.finite(surface$dominant_growth_rate)) &&
    is.finite(max_error) && max_error <= 1e-8
  qc <- data.frame(
    display_label = metadata$display_label,
    pair_label = metadata$pair_label,
    pair_id = metadata$pair_id,
    model_context = model_context,
    parameter_endpoint_group = metadata$parameter_endpoint_group,
    representative_seed_number = seed_number,
    endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10,
    n_fixed_p = length(unique(surface$effective_p_misseg)),
    n_o2_per_fixed_p = if (nrow(surface)) {
      min(table(surface$effective_p_misseg))
    } else {
      0L
    },
    n_surface = nrow(surface),
    all_status_ok = all(surface$status == "ok"),
    all_numeric_outputs_finite =
      all(is.finite(surface$dominant_mean_ploidy)) &&
      all(is.finite(surface$spectral_gap)) &&
      all(is.finite(surface$dominant_growth_rate)),
    max_abs_actual_minus_requested_p_misseg = max_error,
    operator_qc_pass = valid,
    cache_path = normalizePath(cache_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(cache_path, ".tmp-", Sys.getpid())
  saveRDS(
    list(
      metadata = list(
        profile = profile,
        model_source_fingerprint = model_source_fingerprint,
        model_context = model_context,
        display_label = metadata$display_label[[1L]],
        pair_label = metadata$pair_label[[1L]],
        pair_id = metadata$pair_id[[1L]],
        parameter_endpoint_group =
          metadata$parameter_endpoint_group[[1L]],
        representative_seed_number = seed_number,
        endpoint_multiplicity_q10 =
          metadata$endpoint_multiplicity_q10[[1L]],
        config_path = context$config_path,
        parameter_source =
          normalizePath(parameter_source, mustWork = TRUE)
      ),
      surface = surface,
      qc = qc
    ),
    temporary,
    compress = "gzip"
  )
  if (!file.copy(temporary, cache_path, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Failed to publish Figure 6D dense cache: ", cache_path)
  }
  unlink(temporary)
  qc
}

f6r_figure6d_summarize_dense_caches <- function(
    paths, endpoint_manifest, cache_paths
) {
  f6r_require_packages("matrixStats")
  summaries <- list()
  for (pair in endpoint_manifest$display_manifest$pair_id) {
    metadata <- endpoint_manifest$endpoints[
      endpoint_manifest$endpoints$pair_id == pair,
      , drop = FALSE
    ]
    objects <- lapply(cache_paths[metadata$parameter_endpoint_group], readRDS)
    if (any(vapply(
      objects, function(x) !isTRUE(x$qc$operator_qc_pass[[1L]]), logical(1L)
    ))) {
      stop("Figure 6D encountered a non-passing dense endpoint cache for ", pair)
    }
    reference <- objects[[1L]]$surface
    key <- paste(reference$effective_p_misseg, reference$O2_pct, sep = "|")
    same_grid <- vapply(objects, function(x) {
      identical(
        paste(x$surface$effective_p_misseg, x$surface$O2_pct, sep = "|"),
        key
      )
    }, logical(1L))
    if (!all(same_grid)) {
      stop("Figure 6D dense endpoint caches do not share one grid for ", pair)
    }
    multiplicity <- metadata$endpoint_multiplicity_q10[
      match(
        vapply(
          objects,
          function(x) x$metadata$parameter_endpoint_group,
          character(1L)
        ),
        metadata$parameter_endpoint_group
      )
    ]
    if (sum(multiplicity) != 50L) {
      stop("Figure 6D endpoint multiplicity does not sum to 50 for ", pair)
    }
    expand_index <- rep(seq_along(objects), times = multiplicity)
    ploidy <- do.call(cbind, lapply(
      objects, function(x) x$surface$dominant_mean_ploidy
    ))
    gap <- do.call(cbind, lapply(
      objects, function(x) x$surface$spectral_gap
    ))
    growth <- do.call(cbind, lapply(
      objects, function(x) x$surface$dominant_growth_rate
    ))
    ploidy_seed_weighted <- ploidy[, expand_index, drop = FALSE]
    gap_seed_weighted <- gap[, expand_index, drop = FALSE]
    growth_seed_weighted <- growth[, expand_index, drop = FALSE]
    ploidy_quantiles <- matrixStats::rowQuantiles(
      ploidy_seed_weighted, probs = c(0.10, 0.25, 0.75, 0.90)
    )
    proportion_high <- rowMeans(ploidy_seed_weighted >= 2)
    summaries[[pair]] <- data.frame(
      pair_id = pair,
      pair_label = metadata$pair_label[[1L]],
      O2_pct = reference$O2_pct,
      effective_p_misseg = reference$effective_p_misseg,
      n_seed = 50L,
      n_unique_parameter_endpoint = nrow(metadata),
      dominant_mean_ploidy_median =
        matrixStats::rowMedians(ploidy_seed_weighted),
      dominant_mean_ploidy_q10 = ploidy_quantiles[, 1L],
      dominant_mean_ploidy_q25 = ploidy_quantiles[, 2L],
      dominant_mean_ploidy_q75 = ploidy_quantiles[, 3L],
      dominant_mean_ploidy_q90 = ploidy_quantiles[, 4L],
      proportion_dominant_mean_ploidy_ge_2 = proportion_high,
      ploidy_regime_consensus = pmax(proportion_high, 1 - proportion_high),
      spectral_gap_median = matrixStats::rowMedians(gap_seed_weighted),
      proportion_spectral_gap_below_0p005 =
        rowMeans(gap_seed_weighted < 0.005),
      dominant_growth_rate_median =
        matrixStats::rowMedians(growth_seed_weighted),
      max_abs_actual_minus_requested_p_misseg = max(vapply(
        objects,
        function(x) x$qc$max_abs_actual_minus_requested_p_misseg[[1L]],
        numeric(1L)
      )),
      cutoff = "q10",
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, summaries)
  summary$display_label <- endpoint_manifest$display_manifest$display_label[
    match(summary$pair_id, endpoint_manifest$display_manifest$pair_id)
  ]
  summary <- summary[order(
    match(
      summary$display_label,
      endpoint_manifest$display_manifest$display_label
    ),
    summary$effective_p_misseg, summary$O2_pct
  ), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

f6r_compute_figure6d_dense <- function(
    paths, objective_bundle, n_core = 8L, rebuild = FALSE
) {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  f6r_require_packages(c("Matrix", "Rcpp", "matrixStats"))
  f6r_load_response_engine(paths)
  endpoint_manifest <- f6r_figure6d_endpoint_manifest(
    paths, objective_bundle
  )
  endpoints <- endpoint_manifest$endpoints
  parameters <- objective_bundle$parameters
  cache_root <- file.path(paths$figure6, "figure6d_dense_endpoint_cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_paths <- stats::setNames(
    file.path(
      cache_root, endpoints$pair_label,
      paste0(
        "endpoint_",
        gsub("[^A-Za-z0-9_.-]", "_", endpoints$parameter_endpoint_group),
        "_seed", endpoints$representative_seed_number, ".rds"
      )
    ),
    endpoints$parameter_endpoint_group
  )
  contexts <- lapply(
    endpoint_manifest$display_manifest$pair_id,
    f6r_pair_model_context,
    selected = objective_bundle$selected,
    paths = paths
  )
  names(contexts) <- endpoint_manifest$display_manifest$pair_id
  compute_one <- function(i) {
    f6r_load_response_engine(paths)
    metadata <- endpoints[i, , drop = FALSE]
    message(
      "Figure 6D dense endpoint start: ", metadata$pair_label,
      " seed", metadata$representative_seed_number,
      " (", i, "/", nrow(endpoints), ")"
    )
    result <- f6r_figure6d_compute_endpoint_cache(
      metadata = metadata,
      parameters = parameters,
      context = contexts[[metadata$pair_id[[1L]]]],
      cache_path = cache_paths[[metadata$parameter_endpoint_group[[1L]]]],
      parameter_source = objective_bundle$paths[["parameters"]],
      force_rebuild = rebuild,
      model_source_fingerprint = f6r_model_source_fingerprint(paths)
    )
    message(
      "Figure 6D dense endpoint complete: ", metadata$pair_label,
      " seed", metadata$representative_seed_number,
      " (", i, "/", nrow(endpoints), ")"
    )
    result
  }
  index <- seq_len(nrow(endpoints))
  qc_rows <- f6r_resilient_lapply(index, compute_one, n_core = n_core)
  if (any(vapply(qc_rows, function(x) inherits(x, "try-error"), logical(1L)))) {
    stop("One or more Figure 6D dense endpoint workers failed.")
  }
  qc <- do.call(rbind, qc_rows)
  if (!all(qc$operator_qc_pass)) {
    stop("Figure 6D dense endpoint operator QC failed.")
  }
  summary <- f6r_figure6d_summarize_dense_caches(
    paths, endpoint_manifest, cache_paths
  )
  summary_path <- f6r_write_tsv(
    summary,
    file.path(paths$figure6, "figure6d_fixed_p_curve_family.tsv")
  )
  qc_path <- f6r_write_tsv(
    qc,
    file.path(paths$figure6, "figure6d_dense_endpoint_qc.tsv")
  )

  original <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  original <- original[
    original$cutoff == "q10" &
      original$pair_id %in% endpoint_manifest$display_manifest$pair_id &
      original$effective_p_misseg %in% c(0.005, 0.5),
    c(
      "pair_id", "O2_pct", "effective_p_misseg",
      "dominant_mean_ploidy_median"
    ), drop = FALSE
  ]
  overlap <- merge(
    summary[, c(
      "pair_id", "O2_pct", "effective_p_misseg",
      "dominant_mean_ploidy_median"
    )],
    original,
    by = c("pair_id", "O2_pct", "effective_p_misseg"),
    suffixes = c("_dense", "_original"),
    all = FALSE, sort = FALSE
  )
  overlap_difference <- max(abs(
    overlap$dominant_mean_ploidy_median_dense -
      overlap$dominant_mean_ploidy_median_original
  ))
  validation <- data.frame(
    check = c(
      "unique_endpoint_cache_count",
      "represented_seed_count_per_pair",
      "dense_summary_row_count",
      "fixed_p_count_per_pair",
      "oxygen_count_per_fixed_p",
      "fixed_p_range_and_interval",
      "all_endpoint_operator_qc_pass",
      "maximum_forcing_error",
      "overlap_row_count_with_original_surface",
      "overlap_matches_original_surface"
    ),
    observed = c(
      nrow(qc),
      paste(tapply(
        qc$endpoint_multiplicity_q10,
        factor(
          qc$pair_label,
          levels = endpoint_manifest$display_manifest$pair_label
        ),
        sum
      ), collapse = ","),
      nrow(summary),
      paste(as.integer(table(summary$pair_label) / 201L), collapse = ","),
      paste(sort(unique(as.integer(table(interaction(
        summary$pair_id, summary$effective_p_misseg, drop = TRUE
      ))))), collapse = ","),
      paste0(
        min(summary$effective_p_misseg), "--",
        max(summary$effective_p_misseg), " by ",
        min(diff(sort(unique(summary$effective_p_misseg))))
      ),
      all(qc$operator_qc_pass),
      max(qc$max_abs_actual_minus_requested_p_misseg),
      nrow(overlap),
      overlap_difference
    ),
    expected = c(
      as.character(nrow(qc)), paste(rep(50L, f6r_family_count()), collapse = ","),
      as.character(f6r_family_count() * 496L * 201L),
      paste(rep(496L, f6r_family_count()), collapse = ","),
      "201", "0.005--0.5 by 0.001", "TRUE", "<=1e-8",
      as.character(f6r_family_count() * 2L * 201L), "<=1e-12"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- c(
    nrow(qc) > 0L,
    identical(
      as.integer(tapply(
        qc$endpoint_multiplicity_q10,
        factor(
          qc$pair_label,
          levels = endpoint_manifest$display_manifest$pair_label
        ),
        sum
      )),
      rep(50L, f6r_family_count())
    ),
    nrow(summary) == f6r_family_count() * 496L * 201L,
    all(table(summary$pair_label) / 201L == 496L),
    all(table(interaction(
      summary$pair_id, summary$effective_p_misseg, drop = TRUE
    )) == 201L),
    isTRUE(all.equal(
      range(summary$effective_p_misseg), c(0.005, 0.5),
      tolerance = 1e-12
    )) &&
      isTRUE(all.equal(
        unique(round(diff(sort(unique(summary$effective_p_misseg))), 12)),
        0.001,
        tolerance = 1e-12
      )),
    all(qc$operator_qc_pass),
    max(qc$max_abs_actual_minus_requested_p_misseg) <= 1e-8,
    nrow(overlap) == f6r_family_count() * 2L * 201L,
    is.finite(overlap_difference) && overlap_difference <= 1e-12
  )
  validation_path <- f6r_write_tsv(
    validation,
    file.path(paths$figure6, "figure6d_dense_model_validation.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Figure 6D dense model validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(list(
    summary = summary,
    qc = qc,
    endpoint_manifest = endpoint_manifest,
    paths = c(
      summary = summary_path,
      qc = qc_path,
      endpoint_manifest = endpoint_manifest$path,
      validation = validation_path
    )
  ))
}

f6r_direction <- function(delta, tolerance = 0.05) {
  ifelse(
    !is.finite(delta), NA_character_,
    ifelse(delta > tolerance, "increase",
           ifelse(delta < -tolerance, "decrease", "approximately flat"))
  )
}

f6r_seed_claim_row <- function(cache_object) {
  tr <- cache_object$trajectory
  sf <- cache_object$surface
  get_o2 <- function(x, o2) x[which.min(abs(x$O2_pct - o2)), , drop = FALSE]
  tr0 <- get_o2(tr, 0)
  tr1 <- get_o2(tr, 1)
  tr5 <- get_o2(tr, 5)
  cin_targets <- c(low = 0.005, mid = 0.05, high = 0.5)
  diagnostic_cin <- unique(vapply(cin_targets, function(cin) {
    sf$effective_p_misseg[[which.min(abs(sf$effective_p_misseg - cin))]]
  }, numeric(1L)))
  cin_delta <- vapply(cin_targets, function(cin) {
    sub <- sf[which.min(abs(sf$effective_p_misseg - cin)), , drop = FALSE]
    cin_used <- sub$effective_p_misseg[[1L]]
    z <- sf[abs(sf$effective_p_misseg - cin_used) < 1e-12, , drop = FALSE]
    z5 <- get_o2(z, 5)$dominant_mean_ploidy[[1L]]
    z0 <- get_o2(z, 0)$dominant_mean_ploidy[[1L]]
    z5 - z0
  }, numeric(1L))
  domain_values <- lapply(c(0, 1, 5), function(o2) {
    z <- sf[abs(sf$O2_pct - o2) < 1e-12, , drop = FALSE]
    c(
      both = any(z$dominant_mean_ploidy >= 2) &&
        any(z$dominant_mean_ploidy < 2),
      high_fraction = mean(z$dominant_mean_ploidy >= 2)
    )
  })
  names(domain_values) <- paste0("o2_", c(0, 1, 5))
  diagnostic_surface <- sf[
    vapply(sf$O2_pct, function(x) any(abs(x - c(0, 1, 5)) < 1e-12), logical(1L)) |
      vapply(
        sf$effective_p_misseg,
        function(x) any(abs(x - diagnostic_cin) < 1e-12),
        logical(1L)
      ),
    , drop = FALSE
  ]
  data.frame(
    pair_id = cache_object$metadata$pair_id,
    pair_label = cache_object$metadata$pair_label,
    seed_number = cache_object$metadata$seed_number,
    objective = cache_object$metadata$objective,
    trajectory_ploidy_o2_0 = tr0$dominant_mean_ploidy[[1L]],
    trajectory_ploidy_o2_1 = tr1$dominant_mean_ploidy[[1L]],
    trajectory_ploidy_o2_5 = tr5$dominant_mean_ploidy[[1L]],
    trajectory_delta_ploidy_o2_5_minus_0 =
      tr5$dominant_mean_ploidy[[1L]] - tr0$dominant_mean_ploidy[[1L]],
    trajectory_direction_o2_0_to_5 = f6r_direction(
      tr5$dominant_mean_ploidy[[1L]] - tr0$dominant_mean_ploidy[[1L]]
    ),
    trajectory_p_misseg_o2_0 = tr0$population_average_p_misseg[[1L]],
    trajectory_p_misseg_o2_1 = tr1$population_average_p_misseg[[1L]],
    trajectory_p_misseg_o2_5 = tr5$population_average_p_misseg[[1L]],
    surface_both_regimes_o2_0 = as.logical(domain_values$o2_0[["both"]]),
    surface_both_regimes_o2_1 = as.logical(domain_values$o2_1[["both"]]),
    surface_both_regimes_o2_5 = as.logical(domain_values$o2_5[["both"]]),
    surface_high_fraction_o2_0 = as.numeric(domain_values$o2_0[["high_fraction"]]),
    surface_high_fraction_o2_1 = as.numeric(domain_values$o2_1[["high_fraction"]]),
    surface_high_fraction_o2_5 = as.numeric(domain_values$o2_5[["high_fraction"]]),
    surface_delta_ploidy_o2_5_minus_0_cin_low = cin_delta[["low"]],
    surface_delta_ploidy_o2_5_minus_0_cin_mid = cin_delta[["mid"]],
    surface_delta_ploidy_o2_5_minus_0_cin_high = cin_delta[["high"]],
    surface_direction_o2_0_to_5_cin_low = f6r_direction(cin_delta[["low"]]),
    surface_direction_o2_0_to_5_cin_mid = f6r_direction(cin_delta[["mid"]]),
    surface_direction_o2_0_to_5_cin_high = f6r_direction(cin_delta[["high"]]),
    surface_fraction_spectral_gap_ge_0p005 =
      mean(diagnostic_surface$spectral_gap >= 0.005),
    stringsAsFactors = FALSE
  )
}

f6r_summarize_multiseed <- function(paths, objective_bundle, cache_paths) {
  f6r_require_packages("data.table")
  acceptance <- objective_bundle$objectives
  cutoffs <- list(q05 = "eligible_q05", q10 = "eligible_q10", q20 = "eligible_q20")
  surface_summary <- list()
  surface_unique_summary <- list()
  trajectory_summary <- list()
  trajectory_unique_summary <- list()
  seed_claims <- list()
  qc_rows <- list()

  for (pair in unique(acceptance$pair_id)) {
    pair_cache <- cache_paths[[pair]]
    objects <- lapply(pair_cache, readRDS)
    qc <- do.call(rbind, lapply(objects, `[[`, "qc"))
    qc_rows[[pair]] <- qc
    good <- vapply(objects, function(z) isTRUE(z$qc$operator_qc_pass), logical(1L))
    objects <- objects[good]
    if (!length(objects)) stop("No operator-QC-passing seed caches for ", pair)
    seed_claims[[pair]] <- do.call(rbind, lapply(objects, f6r_seed_claim_row))
    surface_all <- data.table::rbindlist(lapply(objects, `[[`, "surface"))
    trajectory_all <- data.table::rbindlist(lapply(objects, `[[`, "trajectory"))
    pair_acceptance <- acceptance[acceptance$pair_id == pair, , drop = FALSE]
    for (cutoff in names(cutoffs)) {
      eligible_rows <- pair_acceptance[
        pair_acceptance[[cutoffs[[cutoff]]]] & pair_acceptance$hard_qc_pass,
        , drop = FALSE
      ]
      eligible <- eligible_rows$seed_number
      eligible_unique <- eligible_rows$seed_number[
        !duplicated(eligible_rows$parameter_endpoint_group)
      ]
      trajectory <- trajectory_all[seed_number %in% eligible]
      trajectory_unique <- trajectory_all[seed_number %in% eligible_unique]
      if (cutoff %in% c("q05", "q10")) {
        surface <- surface_all[
          seed_number %in% eligible & surface_profile == "full_201x60"
        ]
        surface_unique <- surface_all[
          seed_number %in% eligible_unique & surface_profile == "full_201x60"
        ]
        surface_summary[[paste(pair, cutoff)]] <- surface[, .(
          n_seed = data.table::uniqueN(seed_number),
          dominant_mean_ploidy_median = stats::median(dominant_mean_ploidy),
          dominant_mean_ploidy_q10 = stats::quantile(dominant_mean_ploidy, 0.10),
          dominant_mean_ploidy_q25 = stats::quantile(dominant_mean_ploidy, 0.25),
          dominant_mean_ploidy_q75 = stats::quantile(dominant_mean_ploidy, 0.75),
          dominant_mean_ploidy_q90 = stats::quantile(dominant_mean_ploidy, 0.90),
          proportion_dominant_mean_ploidy_ge_2 = mean(dominant_mean_ploidy >= 2),
          ploidy_regime_consensus = max(
            mean(dominant_mean_ploidy >= 2),
            mean(dominant_mean_ploidy < 2)
          ),
          spectral_gap_median = stats::median(spectral_gap),
          proportion_spectral_gap_below_0p005 = mean(spectral_gap < 0.005),
          dominant_growth_rate_median = stats::median(dominant_growth_rate),
          max_abs_actual_minus_requested_p_misseg = max(abs(
            actual_effective_p_misseg - effective_p_misseg
          ))
        ), by = .(pair_id, pair_label, O2_pct, effective_p_misseg)][
          , cutoff := cutoff
        ]
        surface_unique_summary[[paste(pair, cutoff)]] <- surface_unique[, .(
          n_seed = data.table::uniqueN(seed_number),
          dominant_mean_ploidy_median = stats::median(dominant_mean_ploidy),
          dominant_mean_ploidy_q10 = stats::quantile(dominant_mean_ploidy, 0.10),
          dominant_mean_ploidy_q90 = stats::quantile(dominant_mean_ploidy, 0.90),
          proportion_dominant_mean_ploidy_ge_2 = mean(dominant_mean_ploidy >= 2),
          ploidy_regime_consensus = max(
            mean(dominant_mean_ploidy >= 2),
            mean(dominant_mean_ploidy < 2)
          ),
          spectral_gap_median = stats::median(spectral_gap),
          proportion_spectral_gap_below_0p005 = mean(spectral_gap < 0.005),
          dominant_growth_rate_median = stats::median(dominant_growth_rate)
        ), by = .(pair_id, pair_label, O2_pct, effective_p_misseg)][
          , cutoff := cutoff
        ]
      }
      trajectory_summary[[paste(pair, cutoff)]] <- trajectory[, .(
        n_seed = data.table::uniqueN(seed_number),
        population_average_p_misseg_median = stats::median(
          population_average_p_misseg
        ),
        population_average_p_misseg_q10 = stats::quantile(
          population_average_p_misseg, 0.10
        ),
        population_average_p_misseg_q90 = stats::quantile(
          population_average_p_misseg, 0.90
        ),
        dominant_mean_ploidy_median = stats::median(dominant_mean_ploidy),
        dominant_mean_ploidy_q10 = stats::quantile(dominant_mean_ploidy, 0.10),
        dominant_mean_ploidy_q90 = stats::quantile(dominant_mean_ploidy, 0.90),
        spectral_gap_median = stats::median(spectral_gap),
        proportion_spectral_gap_below_0p005 = mean(spectral_gap < 0.005)
      ), by = .(pair_id, pair_label, O2_pct)][, cutoff := cutoff]
      trajectory_unique_summary[[paste(pair, cutoff)]] <- trajectory_unique[, .(
        n_seed = data.table::uniqueN(seed_number),
        population_average_p_misseg_median = stats::median(
          population_average_p_misseg
        ),
        dominant_mean_ploidy_median = stats::median(dominant_mean_ploidy),
        spectral_gap_median = stats::median(spectral_gap),
        proportion_spectral_gap_below_0p005 = mean(spectral_gap < 0.005)
      ), by = .(pair_id, pair_label, O2_pct)][, cutoff := cutoff]
    }
  }
  surface_summary <- as.data.frame(data.table::rbindlist(surface_summary))
  surface_unique_summary <- as.data.frame(
    data.table::rbindlist(surface_unique_summary)
  )
  trajectory_summary <- as.data.frame(data.table::rbindlist(trajectory_summary))
  trajectory_unique_summary <- as.data.frame(
    data.table::rbindlist(trajectory_unique_summary)
  )

  surface_compare_grid <- merge(
    surface_summary,
    surface_unique_summary,
    by = c("pair_id", "pair_label", "O2_pct", "effective_p_misseg", "cutoff"),
    suffixes = c("_seed_weighted", "_unique_endpoint"),
    all = TRUE, sort = FALSE
  )
  surface_comparison <- data.table::as.data.table(surface_compare_grid)[, .(
    n_optimizer_seed = unique(n_seed_seed_weighted),
    n_unique_parameter_endpoint = unique(n_seed_unique_endpoint),
    median_absolute_ploidy_difference = stats::median(abs(
      dominant_mean_ploidy_median_seed_weighted -
        dominant_mean_ploidy_median_unique_endpoint
    )),
    maximum_absolute_ploidy_difference = max(abs(
      dominant_mean_ploidy_median_seed_weighted -
        dominant_mean_ploidy_median_unique_endpoint
    )),
    proportion_grid_same_median_ploidy_regime = mean(
      (dominant_mean_ploidy_median_seed_weighted >= 2) ==
        (dominant_mean_ploidy_median_unique_endpoint >= 2)
    ),
    proportion_grid_same_low_consensus_flag = mean(
      (ploidy_regime_consensus_seed_weighted < 0.80) ==
        (ploidy_regime_consensus_unique_endpoint < 0.80)
    )
  ), by = .(pair_id, pair_label, cutoff)]
  surface_comparison <- as.data.frame(surface_comparison)

  trajectory_compare_grid <- merge(
    trajectory_summary,
    trajectory_unique_summary,
    by = c("pair_id", "pair_label", "O2_pct", "cutoff"),
    suffixes = c("_seed_weighted", "_unique_endpoint"),
    all = TRUE, sort = FALSE
  )
  trajectory_comparison <- data.table::as.data.table(trajectory_compare_grid)[, .(
    n_optimizer_seed = unique(n_seed_seed_weighted),
    n_unique_parameter_endpoint = unique(n_seed_unique_endpoint),
    median_absolute_ploidy_difference = stats::median(abs(
      dominant_mean_ploidy_median_seed_weighted -
        dominant_mean_ploidy_median_unique_endpoint
    )),
    maximum_absolute_ploidy_difference = max(abs(
      dominant_mean_ploidy_median_seed_weighted -
        dominant_mean_ploidy_median_unique_endpoint
    )),
    median_absolute_p_misseg_difference = stats::median(abs(
      population_average_p_misseg_median_seed_weighted -
        population_average_p_misseg_median_unique_endpoint
    )),
    maximum_absolute_p_misseg_difference = max(abs(
      population_average_p_misseg_median_seed_weighted -
        population_average_p_misseg_median_unique_endpoint
    ))
  ), by = .(pair_id, pair_label, cutoff)]
  trajectory_comparison <- as.data.frame(trajectory_comparison)
  seed_claims <- do.call(rbind, seed_claims)
  qc <- do.call(rbind, qc_rows)
  rownames(seed_claims) <- NULL
  rownames(qc) <- NULL

  endpoint_map <- acceptance[, c(
    "pair_id", "seed_number", "parameter_endpoint_group",
    "endpoint_multiplicity_q05", "endpoint_multiplicity_q10",
    "endpoint_multiplicity_q20"
  )]
  seed_claims <- merge(
    seed_claims, endpoint_map,
    by = c("pair_id", "seed_number"), all.x = TRUE, sort = FALSE
  )
  seed_claims <- seed_claims[order(
    seed_claims$pair_label, seed_claims$objective, seed_claims$seed_number
  ), , drop = FALSE]

  acceptance$operator_qc_pass <- qc$operator_qc_pass[
    match(
      paste(acceptance$pair_id, acceptance$seed_number),
      paste(qc$pair_id, qc$seed_number)
    )
  ]
  acceptance$final_primary_acceptable <-
    acceptance$primary_acceptable & acceptance$operator_qc_pass %in% TRUE
  failed_primary_operator <- acceptance$primary_acceptable &
    !(acceptance$operator_qc_pass %in% TRUE)
  acceptance$exclusion_reason[failed_primary_operator] <-
    "failed fixed-O2 operator status/finite-value/forced-input QC"
  f6r_write_tsv(
    acceptance, file.path(paths$figure6, "joint_seed_acceptance.tsv")
  )
  f6r_write_tsv(
    acceptance[, c(
      "pair_id", "pair_label", "seed", "seed_number", "objective_rank",
      "objective", "objective_finite", "seed_id_valid",
      "n_parameter_record", "n_parameter",
      "all_parameter_finite", "all_feasible_at_solution",
      "all_feasible_before_projection", "any_projection_applied",
      "config_available", "hard_qc_pass", "operator_qc_pass"
    )],
    file.path(paths$figure6, "joint_seed_fit_qc.tsv")
  )

  claim_specs <- list(
    both_o2_0 = list(
      label = "Both ploidy regimes present at 0% O2",
      type = "boolean", field = "surface_both_regimes_o2_0"
    ),
    both_o2_1 = list(
      label = "Both ploidy regimes present at 1% O2",
      type = "boolean", field = "surface_both_regimes_o2_1"
    ),
    both_o2_5 = list(
      label = "Both ploidy regimes present at 5% O2",
      type = "boolean", field = "surface_both_regimes_o2_5"
    ),
    trajectory_direction = list(
      label = "Unmodified trajectory: O2 response direction",
      type = "direction", field = "trajectory_direction_o2_0_to_5"
    ),
    surface_direction_low_cin = list(
      label = "Surface O2 response at low effective missegregation",
      type = "direction", field = "surface_direction_o2_0_to_5_cin_low"
    ),
    surface_direction_mid_cin = list(
      label = "Surface O2 response at intermediate effective missegregation",
      type = "direction", field = "surface_direction_o2_0_to_5_cin_mid"
    ),
    surface_direction_high_cin = list(
      label = "Surface O2 response at high effective missegregation",
      type = "direction", field = "surface_direction_o2_0_to_5_cin_high"
    ),
    spectral_gap_majority_reliable = list(
      label = "At least half of the prespecified diagnostic union has spectral gap >= 0.005",
      type = "boolean_numeric",
      field = "surface_fraction_spectral_gap_ge_0p005"
    )
  )

  summarize_pair_claims <- function(
      z, cutoff, weighting, n_optimizer_seed, n_unique_endpoint
  ) {
    rows <- list()
    for (claim_id in names(claim_specs)) {
      spec <- claim_specs[[claim_id]]
      value <- z[[spec$field]]
      if (spec$type == "direction") {
        counts <- table(value, useNA = "no")
        modal <- names(counts)[which.max(counts)]
        support <- mean(value == modal)
        median_value <- NA_real_
        support_definition <- paste0("modal direction: ", modal)
      } else if (spec$type == "boolean_numeric") {
        logical_value <- value >= 0.5
        fraction_true <- mean(logical_value)
        modal <- if (fraction_true >= 0.5) "TRUE" else "FALSE"
        support <- max(fraction_true, 1 - fraction_true)
        median_value <- stats::median(value)
        support_definition <- paste0("modal state: ", modal)
      } else {
        logical_value <- as.logical(value)
        fraction_true <- mean(logical_value)
        modal <- if (fraction_true >= 0.5) "TRUE" else "FALSE"
        support <- max(fraction_true, 1 - fraction_true)
        median_value <- fraction_true
        support_definition <- paste0("modal state: ", modal)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        cutoff = cutoff,
        pair_id = z$pair_id[[1L]],
        pair_label = z$pair_label[[1L]],
        claim_id = claim_id,
        claim_label = spec$label,
        analysis_weighting = weighting,
        n_seed = n_optimizer_seed,
        n_unique_parameter_endpoint = n_unique_endpoint,
        n_analysis_unit = nrow(z),
        modal_result = modal,
        modal_support_fraction = support,
        median_estimand = median_value,
        support_definition = support_definition,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  robustness_rows <- list()
  unique_endpoint_rows <- list()
  for (cutoff in names(cutoffs)) {
    eligible_map <- acceptance[
      acceptance[[cutoffs[[cutoff]]]] & acceptance$operator_qc_pass %in% TRUE,
      c("pair_id", "seed_number", "parameter_endpoint_group"), drop = FALSE
    ]
    claims <- merge(
      seed_claims, eligible_map,
      by = c("pair_id", "seed_number", "parameter_endpoint_group"), all = FALSE
    )
    for (pair in unique(claims$pair_id)) {
      z <- claims[claims$pair_id == pair, , drop = FALSE]
      z <- z[order(z$objective, z$seed_number), , drop = FALSE]
      z_unique <- z[!duplicated(z$parameter_endpoint_group), , drop = FALSE]
      n_unique <- nrow(z_unique)
      robustness_rows[[length(robustness_rows) + 1L]] <- summarize_pair_claims(
        z, cutoff, "optimizer-seed endpoints", nrow(z), n_unique
      )
      unique_endpoint_rows[[length(unique_endpoint_rows) + 1L]] <-
        summarize_pair_claims(
          z_unique, cutoff, "unique 14-parameter endpoints",
          nrow(z), n_unique
        )
    }
  }
  robustness <- do.call(rbind, robustness_rows)
  robustness_unique_endpoint <- do.call(rbind, unique_endpoint_rows)
  dedup_comparison <- merge(
    robustness,
    robustness_unique_endpoint,
    by = c("cutoff", "pair_id", "pair_label", "claim_id", "claim_label"),
    suffixes = c("_seed_weighted", "_unique_endpoint"),
    all = TRUE, sort = FALSE
  )
  dedup_comparison$same_modal_result <- with(
    dedup_comparison,
    modal_result_seed_weighted == modal_result_unique_endpoint
  )
  dedup_comparison$absolute_support_difference <- with(
    dedup_comparison,
    abs(modal_support_fraction_seed_weighted -
          modal_support_fraction_unique_endpoint)
  )

  claim_fields <- unique(vapply(claim_specs, `[[`, character(1L), "field"))
  claim_signature <- do.call(
    paste,
    c(lapply(seed_claims[, claim_fields, drop = FALSE], as.character), sep = "|")
  )
  endpoint_claim_consistency <- data.table::as.data.table(transform(
    seed_claims, claim_signature = claim_signature
  ))[, .(
    n_optimizer_seed = .N,
    n_distinct_claim_profile = data.table::uniqueN(claim_signature),
    all_duplicate_seed_claims_identical =
      data.table::uniqueN(claim_signature) == 1L
  ), by = .(pair_id, pair_label, parameter_endpoint_group)]
  endpoint_claim_consistency <- as.data.frame(endpoint_claim_consistency)
  cutoff_consistency <- do.call(rbind, lapply(
    split(robustness, interaction(robustness$pair_id, robustness$claim_id)),
    function(z) {
      z <- z[match(c("q05", "q10", "q20"), z$cutoff), , drop = FALSE]
      data.frame(
        pair_id = z$pair_id[[1L]],
        pair_label = z$pair_label[[1L]],
        claim_id = z$claim_id[[1L]],
        claim_label = z$claim_label[[1L]],
        n_cutoff = nrow(z),
        same_modal_result_all_cutoffs =
          length(unique(z$modal_result)) == 1L,
        modal_result_q05 = z$modal_result[z$cutoff == "q05"],
        modal_result_q10 = z$modal_result[z$cutoff == "q10"],
        modal_result_q20 = z$modal_result[z$cutoff == "q20"],
        minimum_modal_support_across_cutoffs =
          min(z$modal_support_fraction),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(cutoff_consistency) <- NULL
  overall <- do.call(rbind, lapply(
    split(robustness, interaction(robustness$cutoff, robustness$claim_id)),
    function(z) {
      mode_counts <- table(z$modal_result)
      mode <- names(mode_counts)[which.max(mode_counts)]
      data.frame(
        cutoff = z$cutoff[[1L]],
        claim_id = z$claim_id[[1L]],
        claim_label = z$claim_label[[1L]],
        n_pair = nrow(z),
        cross_pair_modal_result = mode,
        n_pair_with_cross_pair_mode = sum(z$modal_result == mode),
        pair_support_fraction_mean = mean(z$modal_support_fraction),
        pair_support_fraction_min = min(z$modal_support_fraction),
        pair_support_fraction_max = max(z$modal_support_fraction),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(overall) <- NULL

  response_delta_fields <- c(
    trajectory = "trajectory_delta_ploidy_o2_5_minus_0",
    fixed_cin_low = "surface_delta_ploidy_o2_5_minus_0_cin_low",
    fixed_cin_intermediate = "surface_delta_ploidy_o2_5_minus_0_cin_mid",
    fixed_cin_high = "surface_delta_ploidy_o2_5_minus_0_cin_high"
  )
  response_delta_rows <- list()
  for (cutoff in names(cutoffs)) {
    eligible_map <- acceptance[
      acceptance[[cutoffs[[cutoff]]]] & acceptance$operator_qc_pass %in% TRUE,
      c("pair_id", "seed_number"), drop = FALSE
    ]
    claims <- merge(
      seed_claims, eligible_map,
      by = c("pair_id", "seed_number"), all = FALSE
    )
    for (pair in unique(claims$pair_id)) {
      z <- claims[claims$pair_id == pair, , drop = FALSE]
      for (estimand in names(response_delta_fields)) {
        value <- z[[response_delta_fields[[estimand]]]]
        response_delta_rows[[length(response_delta_rows) + 1L]] <- data.frame(
          cutoff = cutoff,
          pair_id = pair,
          pair_label = z$pair_label[[1L]],
          estimand = estimand,
          definition = "dominant mean ploidy at 5% O2 minus value at 0% O2",
          n_seed = length(value),
          median = stats::median(value),
          q10 = unname(stats::quantile(value, 0.10)),
          q90 = unname(stats::quantile(value, 0.90)),
          minimum = min(value),
          maximum = max(value),
          modal_direction = names(which.max(table(f6r_direction(value)))),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  response_delta_summary <- do.call(rbind, response_delta_rows)
  rownames(response_delta_summary) <- NULL

  q10_unique_counts <- vapply(
    split(acceptance[acceptance$eligible_q10, ],
          acceptance$pair_id[acceptance$eligible_q10]),
    function(z) length(unique(z$parameter_endpoint_group)),
    integer(1L)
  )
  robustness_audit_summary <- data.frame(
    metric = c(
      "operator_qc_passing_endpoints",
      "pair_claim_combinations",
      "pair_claim_modal_results_invariant_q05_q10_q20",
      "minimum_modal_support_across_all_pairs_claims_cutoffs",
      "seed_vs_unique_endpoint_claim_comparisons",
      "seed_vs_unique_endpoint_same_modal_result",
      "minimum_q10_grid_same_median_ploidy_regime_after_deduplication",
      "minimum_q10_unique_parameter_endpoints_per_pair",
      "maximum_q10_unique_parameter_endpoints_per_pair"
    ),
    value = c(
      paste0(sum(qc$operator_qc_pass), "/", nrow(qc)),
      nrow(cutoff_consistency),
      sum(cutoff_consistency$same_modal_result_all_cutoffs),
      min(cutoff_consistency$minimum_modal_support_across_cutoffs),
      nrow(dedup_comparison),
      sum(dedup_comparison$same_modal_result),
      min(surface_comparison$proportion_grid_same_median_ploidy_regime[
        surface_comparison$cutoff == "q10"
      ]),
      min(q10_unique_counts),
      max(q10_unique_counts)
    ),
    interpretation = c(
      "All q20 endpoints passed operator status, finite-value, and forced-input checks",
      paste0(
        f6r_family_count(),
        " warm-start families times eight prespecified qualitative diagnostics"
      ),
      "Modal result unchanged when the objective set widened from 5% to 10% to 20%",
      "Smallest exact within-set modal support over all qualitative audits",
      "All pair-by-claim-by-cutoff comparisons",
      "Modal result unchanged after exact in-vivo response-parameter endpoint deduplication",
      "Smallest pointwise regime agreement for the primary complete surfaces",
      "Smallest number of exact unique in-vivo response-parameter endpoints among 50 q10 seeds",
      "Largest number of exact unique in-vivo response-parameter endpoints among 50 q10 seeds"
    ),
    stringsAsFactors = FALSE
  )

  output_paths <- c(
    surface = f6r_write_tsv(
      surface_summary,
      file.path(paths$figure6, "joint_multiseed_surface_summary.tsv")
    ),
    trajectory = f6r_write_tsv(
      trajectory_summary,
      file.path(paths$figure6, "joint_multiseed_trajectory_summary.tsv")
    ),
    surface_unique_endpoint = f6r_write_tsv(
      surface_unique_summary,
      file.path(paths$figure6, "joint_unique_parameter_surface_summary.tsv")
    ),
    trajectory_unique_endpoint = f6r_write_tsv(
      trajectory_unique_summary,
      file.path(paths$figure6, "joint_unique_parameter_trajectory_summary.tsv")
    ),
    surface_dedup_comparison = f6r_write_tsv(
      surface_comparison,
      file.path(paths$figure6, "joint_seed_vs_unique_parameter_surface_comparison.tsv")
    ),
    trajectory_dedup_comparison = f6r_write_tsv(
      trajectory_comparison,
      file.path(paths$figure6, "joint_seed_vs_unique_parameter_trajectory_comparison.tsv")
    ),
    seed_claims = f6r_write_tsv(
      seed_claims,
      file.path(paths$figure6, "joint_seed_biological_robustness.tsv")
    ),
    robustness = f6r_write_tsv(
      robustness,
      file.path(paths$figure6, "joint_seed_claim_robustness.tsv")
    ),
    robustness_overall = f6r_write_tsv(
      overall,
      file.path(paths$figure6, "joint_seed_claim_robustness_overall.tsv")
    ),
    robustness_unique_endpoint = f6r_write_tsv(
      robustness_unique_endpoint,
      file.path(paths$figure6, "joint_unique_parameter_claim_robustness.tsv")
    ),
    dedup_comparison = f6r_write_tsv(
      dedup_comparison,
      file.path(paths$figure6, "joint_seed_vs_unique_parameter_robustness.tsv")
    ),
    endpoint_claim_consistency = f6r_write_tsv(
      endpoint_claim_consistency,
      file.path(paths$figure6, "joint_duplicate_endpoint_claim_consistency.tsv")
    ),
    cutoff_consistency = f6r_write_tsv(
      cutoff_consistency,
      file.path(paths$figure6, "joint_seed_cutoff_consistency.tsv")
    ),
    response_delta_summary = f6r_write_tsv(
      response_delta_summary,
      file.path(paths$figure6, "joint_seed_response_delta_summary.tsv")
    ),
    robustness_audit_summary = f6r_write_tsv(
      robustness_audit_summary,
      file.path(paths$figure6, "joint_seed_robustness_audit_summary.tsv")
    ),
    qc = f6r_write_tsv(
      qc, file.path(paths$figure6, "joint_multiseed_operator_qc.tsv")
    )
  )
  invisible(list(
    surface = surface_summary,
    trajectory = trajectory_summary,
    surface_unique_endpoint = surface_unique_summary,
    trajectory_unique_endpoint = trajectory_unique_summary,
    surface_dedup_comparison = surface_comparison,
    trajectory_dedup_comparison = trajectory_comparison,
    seed_claims = seed_claims,
    robustness = robustness,
    robustness_unique_endpoint = robustness_unique_endpoint,
    dedup_comparison = dedup_comparison,
    endpoint_claim_consistency = endpoint_claim_consistency,
    cutoff_consistency = cutoff_consistency,
    response_delta_summary = response_delta_summary,
    robustness_audit_summary = robustness_audit_summary,
    overall = overall,
    qc = qc,
    paths = output_paths
  ))
}

f6r_compute_multiseed <- function(
    paths, objective_bundle, n_core = 8L, rebuild = FALSE
) {
  # Each numerical seed is already parallelized at the process level. Limit
  # nested BLAS/OpenMP threading to avoid oversubscription and nondeterministic
  # throughput when several endpoint operators are evaluated concurrently.
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  f6r_require_packages(c("Matrix", "Rcpp", "data.table"))
  f6r_load_response_engine(paths)
  acceptance <- objective_bundle$objectives
  selected <- objective_bundle$selected
  master <- objective_bundle$master
  cache_root <- file.path(paths$figure6, "multiseed_seed_cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_paths <- list()
  qc_rows <- list()
  jobs <- list()
  job_pairs <- character()

  for (pair in unique(acceptance$pair_id)) {
    context <- f6r_pair_model_context(pair, selected, paths)
    eligible <- acceptance[
      acceptance$pair_id == pair & acceptance$eligible_q20 &
        acceptance$hard_qc_pass,
      , drop = FALSE
    ]
    eligible <- eligible[order(eligible$objective_rank), , drop = FALSE]
    if (nrow(eligible) != 100L) {
      stop("Expected 100 q20 seeds for ", pair, "; found ", nrow(eligible))
    }
    pair_cache_dir <- file.path(cache_root, context$pair_label)
    dir.create(pair_cache_dir, recursive = TRUE, showWarnings = FALSE)
    pair_cache <- stats::setNames(
      file.path(pair_cache_dir, paste0("seed", eligible$seed_number, ".rds")),
      eligible$seed_number
    )
    if (isTRUE(rebuild)) {
      # Rebuild is implemented by ignoring existing caches; each seed cache is
      # atomically replaced after a successful calculation.
      existing_ok <- rep(FALSE, length(pair_cache))
    } else {
      existing_ok <- vapply(pair_cache, function(path) {
        if (!file.exists(path)) return(FALSE)
        cached <- tryCatch(readRDS(path), error = function(e) NULL)
        !is.null(cached) && isTRUE(cached$qc$operator_qc_pass[[1L]])
      }, logical(1L))
    }
    message(
      "Figure 6 multi-seed ", context$pair_label, ": ",
      sum(existing_ok), "/", length(pair_cache), " seed caches already available."
    )
    pair_jobs <- lapply(seq_len(nrow(eligible)), function(j) {
      list(
        pair_id = pair,
        pair_label = context$pair_label,
        seed_number = eligible$seed_number[[j]],
        objective = eligible$objective[[j]],
        context = context,
        cache_path = pair_cache[[as.character(eligible$seed_number[[j]])]],
        full_surface = eligible$eligible_q10[[j]]
      )
    })
    jobs <- c(jobs, pair_jobs)
    job_pairs <- c(job_pairs, rep(pair, length(pair_jobs)))
    cache_paths[[pair]] <- unname(pair_cache)
  }

  # Use one fixed fork pool for the complete primary-family task set. On macOS,
  # repeatedly constructing a new fork pool for each family eventually causes
  # Intel/LLVM OpenMP shared-memory initialization to fail, even with nested
  # OpenMP threads disabled.
  compute_job <- function(index) {
    f6r_load_response_engine(paths)
    job <- jobs[[index]]
    f6r_compute_seed_cache(
      pair_id = job$pair_id,
      seed_number = job$seed_number,
      objective = job$objective,
      master = master,
      context = job$context,
      cache_path = job$cache_path,
      parameter_source = objective_bundle$master_path,
      full_surface = job$full_surface,
      force_rebuild = rebuild,
      model_source_fingerprint = f6r_model_source_fingerprint(paths)
    )
  }
  job_qc <- f6r_resilient_lapply(
    seq_along(jobs), compute_job, n_core = n_core
  )
  if (any(vapply(
    job_qc, function(x) inherits(x, "try-error"), logical(1L)
  ))) {
    stop("One or more primary-family multi-seed workers failed.")
  }
  for (pair in unique(job_pairs)) {
    pair_indices <- which(job_pairs == pair)
    pair_qc <- job_qc[pair_indices]
    pair_label <- jobs[[pair_indices[[1L]]]]$pair_label
    if (length(pair_qc) != 100L) {
      stop("Expected 100 completed multi-seed jobs for ", pair_label, ".")
    }
    message("  ", pair_label, " complete (100/100).")
    qc_rows[[pair]] <- do.call(rbind, pair_qc)
  }

  summary <- f6r_summarize_multiseed(paths, objective_bundle, cache_paths)
  invisible(list(
    cache_paths = cache_paths,
    computation_qc = do.call(rbind, qc_rows),
    summary = summary
  ))
}

f6r_chart_contract <- function(paths) {
  contract <- data.frame(
    artifact = c(
      "Supplementary Figure 6-1A", "Supplementary Figure 6-1B",
      "Supplementary Figure 6-1C", "Figure 6A", "Figure 6B",
      "Supplementary Figure 6-2A", "Supplementary Figure 6-2B",
      "Supplementary Figure 6-2C", "Supplementary Figure 6-2D"
    ),
    analytical_question = c(
      "Which fixed-O2 response shapes occur across the 500 separate in-vivo fitted endpoints?",
      "Do fixed-O2 response classes occupy distinct regions of pooled parameter space?",
      "Do fixed-O2 response classes differ in full MAP fit-quality distributions?",
      "Does the oxygen-CIN-ploidy topology persist across objective-eligible joint endpoints?",
      "At each oxygen concentration and target dominant ploidy, is there a stable unique effective missegregation probability that reproduces that target?",
      "How well supported is the saved primary warm-start partition?",
      "How stable is the saved primary partition under numerical and subsampling perturbations?",
      "Which numerical endpoints satisfy the objective-based eligibility rule?",
      "Do qualitative oxygen-response diagnostics agree across objective-eligible endpoints and cutoffs?"
    ),
    chart_family = c(
      "four-by-two small-multiple response curves",
      "parameter-space scatter",
      "half-violin, boxplot, and seed-level points",
      "uncertainty heatmap plus trajectory interval",
      "small-multiple inverse-response heatmap with fixed-input reference trajectories",
      "silhouette profile", "faceted silhouette profile",
      "faceted ranked objective curve", "robustness heatmap"
    ),
    data_grain = c(
      "eight response classes summarizing 500 separate in-vivo fitted endpoints across 201 oxygen values",
      "500 separate in-vivo fitted endpoints",
      "500 separate in-vivo fitted endpoints",
      paste0(
        f6r_family_count(),
        " primary-region pairs x oxygen x effective missegregation, summarized over 50 seeds"
      ),
      paste0(
        f6r_family_count(),
        " primary-region pairs x 201 oxygen values x target-ploidy grid (1.000-7.000 by 0.025), inverted endpoint-wise over 50 seed-weighted endpoints"
      ),
      "candidate k", "primary-partition perturbation",
      "pair x numerical seed", "claim x pair x objective cutoff"
    ),
    primary_encoding = c(
      "individual fitted-endpoint curves by response-class color and pointwise class median in black",
      "response-class color; class-best black ring; warm-start-region outlines",
      "full-MAP objective density, quartiles, all seed-level endpoints, and global lowest-decile cutoff",
      "original-unit ploidy labels with log-scaled seed-weighted median fill; trajectory median and 10-90% band; low-consensus marks; unique-parameter endpoint sensitivity in source tables",
      "log-scaled color for the seed-weighted median unique inverse; gray for no stable unique inverse; hatching for endpoint-level multiple solutions; four black fixed-input trajectories and the red arithmetic mean across all 496 fixed-input trajectories",
      "average silhouette by k", "adjusted Rand index under perturbation",
      "delta full-MAP objective by rank",
      "minimum modal support across cutoffs; primary modal result and cutoff-consistency flag"
    ),
    caveat = c(
      "Response classes pool reliable, caution, and unreliable spectral-gap strata and are descriptive",
      "Pooled t-SNE axes are unitless and embedding-dependent; response-class separation is descriptive",
      "Objective distributions are descriptive; numerical endpoints are not biological replicates",
      "Post-fit asymptotic diagnostic; optimizer seeds are not biological replicates",
      "Numerical inverse of the standardized post-fit response surface; target ploidy is not clamped in the model, gray or hatched regions do not define a unique required probability, and the red curve is an unweighted visualization summary rather than a fitted trajectory",
      "t-SNE axes are unitless and embedding-dependent",
      "Weak silhouettes imply coverage strata, not biological subtypes",
      "Pair-specific empirical quantiles are numerical eligibility rules, not confidence sets",
      "Exact support fractions describe numerical endpoints, not inferential uncertainty; repeated exact parameter endpoints are also analyzed once each"
    ),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(contract, file.path(paths$figure6, "figure6_chart_contract.tsv"))
}

f6r_local_data_contract <- function(paths) {
  sources <- c(
    file.path(paths$figure6, "response_class_smoothed_curves.tsv"),
    file.path(paths$figure6, "response_class_curve_class_by_seed.tsv"),
    file.path(paths$figure6, "response_class_class_counts.tsv"),
    file.path(
      paths$figure5,
      "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
    ),
    file.path(
      paths$figure5,
      "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_full_coordinates.csv"
    ),
    file.path(paths$figure4, "all_parameter_fitted_endpoint_values.tsv"),
    file.path(paths$supp5_1, "soft_coupling_master_long.tsv"),
    unlist(lapply(
      f6r_selected_results(paths)$warmup_label,
      function(label) file.path(paths$supp5_2, "joint", label, "seed_objective_simple.tsv")
    )),
    unlist(lapply(
      f6r_selected_results(paths)$bundle_path,
      function(bundle) c(
        file.path(bundle, "best_params.tsv"), file.path(bundle, "fit_config.rds")
      )
    )),
    file.path(paths$code, "util", "oxygen", "response_pipeline.R"),
    file.path(paths$code, "util", "analysis", "figure6_robustness.R")
  )
  sources <- unique(sources)
  f6r_require_files(sources, "Figure 6 local-contract source")
  root_prefix <- paste0(paths$root, .Platform$file.sep)
  relative <- ifelse(
    startsWith(normalizePath(sources), root_prefix),
    substring(normalizePath(sources), nchar(root_prefix) + 1L),
    normalizePath(sources)
  )
  contract <- data.frame(
    figure_id = "Figure6",
    role = "local scientific or analytical input",
    source = relative,
    local_file = normalizePath(sources),
    source_md5 = vapply(sources, f6r_md5, character(1L)),
    local_md5 = vapply(sources, f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(contract, file.path(paths$figure6, "data_contract.tsv"))
}

f6r_prepare_invivo_response_classes <- function(paths) {
  input_path <- file.path(
    paths$figure4, "fixed_o2_dominant_ploidy_201grid.tsv"
  )
  f6r_require_files(input_path, "fresh in-vivo fixed-O2 grid")
  curves <- f6r_read_tsv(input_path)
  required <- c(
    "seed_id", "seed_number", "O2_pct", "dominant_mean_ploidy",
    "spectral_gap", "objective"
  )
  if (!all(required %in% names(curves))) {
    stop(
      "Fresh in-vivo fixed-O2 grid lacks: ",
      paste(setdiff(required, names(curves)), collapse = ", ")
    )
  }
  if (nrow(curves) != 500L * 201L ||
      length(unique(curves$seed_id)) != 500L ||
      length(unique(curves$O2_pct)) != 201L ||
      !isTRUE(all.equal(range(curves$O2_pct), c(0, 5))) ||
      anyDuplicated(curves[, c("seed_id", "O2_pct")]) ||
      any(!is.finite(curves$dominant_mean_ploidy)) ||
      any(!is.finite(curves$spectral_gap))) {
    stop("Fresh in-vivo fixed-O2 response grid failed its 500 x 201 contract.")
  }

  response_environment <- new.env(parent = globalenv())
  sys.source(
    file.path(paths$code, "util", "oxygen", "response_pipeline.R"),
    envir = response_environment,
    chdir = TRUE
  )
  classified <- response_environment$response_classify_all_curves(curves)
  by_seed <- classified$summary
  by_seed$seed_number <- as.integer(sub("^seed", "", by_seed$seed_id))
  ranking_path <- file.path(
    paths$figure4, "invivo_fit_objective_ranking_500seeds.tsv"
  )
  f6r_require_files(ranking_path, "canonical in-vivo objective ranking")
  objective_map <- f6r_read_tsv(ranking_path)[, c("seed", "objective")]
  names(objective_map)[[1L]] <- "seed_id"
  by_seed <- merge(by_seed, objective_map, by = "seed_id", all.x = TRUE)

  gap_rows <- lapply(split(curves, curves$seed_id), function(z) {
    data.frame(
      seed_id = z$seed_id[[1L]],
      fraction_gap_below_0p005 = mean(z$spectral_gap < 0.005),
      fraction_gap_below_0p01 = mean(z$spectral_gap < 0.01),
      any_gap_below_0p005 = any(z$spectral_gap < 0.005),
      stringsAsFactors = FALSE
    )
  })
  gap_summary <- do.call(rbind, gap_rows)
  gap_summary$spectral_gap_class <- ifelse(
    gap_summary$fraction_gap_below_0p005 >= 0.25, "unreliable",
    ifelse(
      gap_summary$any_gap_below_0p005 |
        gap_summary$fraction_gap_below_0p01 >= 0.10,
      "caution", "reliable"
    )
  )
  by_seed <- merge(by_seed, gap_summary, by = "seed_id", all.x = TRUE)
  by_seed$context <- "in vivo"
  classified$smooth$context <- "in vivo"
  if (nrow(classified$segments)) {
    classified$segments$context <- "in vivo"
  }

  class_order <- response_environment$response_curve_class_order
  class_counts <- data.frame(
    smooth_curve_class = class_order,
    n_seed = as.integer(table(factor(
      by_seed$smooth_curve_class, levels = class_order
    ))),
    stringsAsFactors = FALSE
  )
  class_counts$fraction_seed <- class_counts$n_seed / 500L
  reliability_order <- c("reliable", "caution", "unreliable")
  reliability_counts <- data.frame(
    spectral_gap_class = reliability_order,
    n_seed = as.integer(table(factor(
      by_seed$spectral_gap_class, levels = reliability_order
    ))),
    stringsAsFactors = FALSE
  )
  reliability_counts$fraction_seed <- reliability_counts$n_seed / 500L

  outputs <- c(
    smoothed = f6r_write_tsv(
      classified$smooth,
      file.path(paths$figure6, "response_class_smoothed_curves.tsv")
    ),
    by_seed = f6r_write_tsv(
      by_seed,
      file.path(paths$figure6, "response_class_curve_class_by_seed.tsv")
    ),
    segments = f6r_write_tsv(
      classified$segments,
      file.path(paths$figure6, "response_class_persistent_segments.tsv")
    ),
    class_counts = f6r_write_tsv(
      class_counts,
      file.path(paths$figure6, "response_class_class_counts.tsv")
    ),
    reliability_counts = f6r_write_tsv(
      reliability_counts,
      file.path(paths$figure6, "response_class_reliability_counts.tsv")
    )
  )
  validation <- data.frame(
    check = c(
      "seed_count", "oxygen_count", "curve_row_count",
      "class_count_total", "reliability_count_total", "best_seed"
    ),
    observed = c(
      length(unique(curves$seed_id)), length(unique(curves$O2_pct)),
      nrow(curves), sum(class_counts$n_seed), sum(reliability_counts$n_seed),
      by_seed$seed_number[[which.min(by_seed$objective)]]
    ),
    expected = c(500L, 201L, 500L * 201L, 500L, 500L, 25L),
    stringsAsFactors = FALSE
  )
  validation$passed <- validation$observed == validation$expected
  validation_path <- f6r_write_tsv(
    validation,
    file.path(paths$figure6, "response_class_validation.tsv")
  )
  provenance <- data.frame(
    role = c(
      "fresh_fixed_o2_grid", "canonical_objective_ranking",
      "classification_code"
    ),
    source_file = c(
      normalizePath(input_path, mustWork = TRUE),
      normalizePath(ranking_path, mustWork = TRUE),
      normalizePath(
        file.path(paths$code, "util", "oxygen", "response_pipeline.R"),
        mustWork = TRUE
      )
    ),
    md5 = unname(tools::md5sum(c(
      input_path,
      ranking_path,
      file.path(paths$code, "util", "oxygen", "response_pipeline.R")
    ))),
    stringsAsFactors = FALSE
  )
  provenance_path <- f6r_write_tsv(
    provenance,
    file.path(paths$figure6, "response_class_source_provenance.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Fresh in-vivo response classification failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(c(
    outputs, validation = validation_path, provenance = provenance_path
  ))
}

f6r_data <- function(
    workspace_root = f6r_find_workspace_root(), n_core = 8L,
    rebuild = FALSE, n_resample = 100L
) {
  paths <- f6r_paths(workspace_root)
  dir.create(paths$figure6, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$supp6_2, recursive = TRUE, showWarnings = FALSE)
  message("Figure 6 response classes: deriving fresh in-vivo classes from Figure 4.")
  response_class_bundle <- f6r_prepare_invivo_response_classes(paths)
  f6r_chart_contract(paths)
  message("Figure 6 robustness: reconstructing cluster-selection evidence.")
  cluster_bundle <- f6r_cluster_diagnostics(paths, n_resample = n_resample)
  message("Figure 6 robustness: applying pair-specific objective eligibility.")
  objective_bundle <- f6r_objective_selection(paths)
  message("Figure 6 robustness: computing/reusing q20 multi-seed response caches.")
  multiseed_bundle <- f6r_compute_multiseed(
    paths, objective_bundle, n_core = n_core, rebuild = rebuild
  )
  message("Figure 6D: computing/reusing exact 0.001-step fixed-missegregation caches.")
  dense_figure6d_bundle <- f6r_compute_figure6d_dense(
    paths, objective_bundle, n_core = n_core, rebuild = rebuild
  )
  message("Figure 6B: computing/reusing in-vivo inverse-response caches.")
  inverse_bundle <- f6r_inverse_panel_data(
    paths, rebuild = rebuild, n_core = n_core,
    dense_qc_path = file.path(paths$figure6, "figure6d_dense_endpoint_qc.tsv"),
    output_prefix = "figure6", model_context = "in vivo"
  )
  f6r_local_data_contract(paths)
  invisible(list(
    paths = paths,
    clusters = cluster_bundle,
    objective = objective_bundle,
    multiseed = multiseed_bundle,
    dense_figure6d = dense_figure6d_bundle,
    inverse = inverse_bundle,
    response_classes = response_class_bundle
  ))
}

f6r_theme <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#222222"),
      axis.text = ggplot2::element_text(size = 7, colour = "#333333"),
      axis.title = ggplot2::element_text(size = 8),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8, face = "bold"),
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.2, colour = "#555555"),
      plot.caption = ggplot2::element_text(size = 6.6, colour = "#555555"),
      plot.tag = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 7.5),
      legend.text = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}

f6r_save_plot <- function(plot, base_path, width, height, dpi = 300) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  png_path <- paste0(base_path, ".png")
  pdf_path <- paste0(base_path, ".pdf")
  ggplot2::ggsave(
    png_path, plot = plot, width = width, height = height,
    units = "in", dpi = dpi, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf_path, plot = plot, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf, bg = "white",
    limitsize = FALSE
  )
  c(
    png = normalizePath(png_path, mustWork = TRUE),
    pdf = normalizePath(pdf_path, mustWork = TRUE)
  )
}

f6r_compose_three_panel_row <- function(
    plots, title, legend_rel_width = 0.78, title_size = 10.5,
    panel_width_mm = 43, panel_height_mm = 43
) {
  family_count <- f6r_family_count()
  if (length(plots) != family_count) {
    stop("An iteration4 Figure 6 row requires one plot per primary region.")
  }
  legend <- cowplot::get_legend(
    plots[[1L]] + ggplot2::theme(legend.position = "right")
  )
  panel_plots <- lapply(plots, function(p) {
    egg::set_panel_size(
      p + ggplot2::theme(legend.position = "none"),
      width = grid::unit(panel_width_mm, "mm"),
      height = grid::unit(panel_height_mm, "mm"),
      margin = grid::unit(0.5, "mm")
    )
  })
  panel_row <- cowplot::plot_grid(
    plotlist = panel_plots, nrow = 1L,
    align = "hv", axis = "tblr", rel_widths = rep(1, family_count)
  )
  body <- cowplot::plot_grid(
    panel_row, legend, nrow = 1L,
    rel_widths = c(family_count, legend_rel_width)
  )
  cowplot::ggdraw() +
    cowplot::draw_plot(body, x = 0, y = 0, width = 1, height = 0.91) +
    cowplot::draw_label(
      title, x = 0.006, y = 0.995,
      hjust = 0, vjust = 1, fontface = "bold",
      fontfamily = "Helvetica", size = title_size
    )
}

f6r_primary_row_width <- function(
    reference_width = 17.2, legend_rel_width = 0.72,
    reference_family_count = 6L
) {
  reference_width *
    (f6r_family_count() + legend_rel_width) /
    (reference_family_count + legend_rel_width)
}

f6r_display_pair_manifest <- function(pair_ids, panel_label) {
  family_levels <- f6r_family_levels()
  display_pair_labels <- stats::setNames(family_levels, family_levels)
  pair_ids <- unique(as.character(pair_ids))
  pair_id_to_label <- stats::setNames(
    sub(".*_(C[0-9]{2})_vt.*", "\\1", pair_ids),
    pair_ids
  )
  ordered_pairs <- names(pair_id_to_label)[
    match(unname(display_pair_labels), pair_id_to_label)
  ]
  if (anyNA(ordered_pairs) ||
      length(unique(ordered_pairs)) != length(family_levels)) {
    stop("Cannot resolve the requested Figure 6", panel_label,
         " pair labels.")
  }
  data.frame(
    display_label = names(display_pair_labels),
    pair_label = unname(display_pair_labels),
    pair_id = ordered_pairs,
    stringsAsFactors = FALSE
  )
}

f6r_weak_gap_hatch_data <- function(
    surface, threshold = 0.5, spacing = 0.14
) {
  x <- sort(unique(surface$O2_pct))
  y <- sort(unique(surface$log10_effective_p_misseg))
  z <- matrix(NA_real_, nrow = length(y), ncol = length(x))
  z[cbind(
    match(surface$log10_effective_p_misseg, y),
    match(surface$O2_pct, x)
  )] <- surface$proportion_spectral_gap_below_0p005
  if (anyNA(z)) stop("Incomplete Figure 6D weak-gap surface grid.")

  bands <- isoband::isobands(
    x, y, z,
    levels_low = threshold,
    levels_high = max(z, na.rm = TRUE) + 1e-8
  )
  geometries <- isoband::iso_to_sfg(bands)
  if (!length(geometries)) {
    return(data.frame(
      O2_pct = numeric(), log10_effective_p_misseg = numeric(),
      hatch_group = integer()
    ))
  }
  region <- sf::st_sfc(geometries[[1L]])
  bounds <- sf::st_bbox(region)
  slope <- diff(range(y)) / diff(range(x))
  intercepts <- seq(
    bounds[["ymin"]] - slope * bounds[["xmax"]],
    bounds[["ymax"]] - slope * bounds[["xmin"]],
    by = spacing
  )
  hatch_lines <- sf::st_sfc(lapply(intercepts, function(intercept) {
    sf::st_linestring(matrix(
      c(
        bounds[["xmin"]], intercept + slope * bounds[["xmin"]],
        bounds[["xmax"]], intercept + slope * bounds[["xmax"]]
      ),
      ncol = 2L, byrow = TRUE
    ))
  }))
  clipped <- suppressWarnings(sf::st_intersection(hatch_lines, region))
  clipped <- suppressWarnings(sf::st_collection_extract(
    clipped, "LINESTRING"
  ))
  clipped <- suppressWarnings(sf::st_cast(clipped, "LINESTRING"))
  if (!length(clipped)) {
    return(data.frame(
      O2_pct = numeric(), log10_effective_p_misseg = numeric(),
      hatch_group = integer()
    ))
  }
  coordinates <- sf::st_coordinates(clipped)
  data.frame(
    O2_pct = coordinates[, "X"],
    log10_effective_p_misseg = coordinates[, "Y"],
    hatch_group = coordinates[, "L1"],
    stringsAsFactors = FALSE
  )
}

f6r_draw_surface_panel <- function(paths) {
  f6r_require_packages(c(
    "ggplot2", "cowplot", "egg", "scales", "isoband", "sf"
  ))
  surface <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  trajectory <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_trajectory_summary.tsv"
  ))
  surface <- surface[surface$cutoff == "q10", , drop = FALSE]
  trajectory <- trajectory[trajectory$cutoff == "q10", , drop = FALSE]
  if (length(unique(surface$pair_id)) != f6r_family_count() ||
      any(table(surface$pair_id) != 201L * 60L) ||
      any(surface$n_seed != 50L) || any(trajectory$n_seed != 50L)) {
    stop("Primary Figure 6A summary must contain 50 seeds per pair.")
  }

  display_manifest <- f6r_display_pair_manifest(surface$pair_id, "A")
  display_pair_labels <- stats::setNames(
    display_manifest$pair_label, display_manifest$display_label
  )
  ordered_pairs <- display_manifest$pair_id
  f6r_write_tsv(
    display_manifest,
    file.path(paths$figure6, "figure6a_displayed_pairs.tsv")
  )
  surface <- surface[surface$pair_id %in% ordered_pairs, , drop = FALSE]
  trajectory <- trajectory[trajectory$pair_id %in% ordered_pairs, , drop = FALSE]
  surface$log10_effective_p_misseg <- log10(surface$effective_p_misseg)
  trajectory$log10_p_median <- log10(
    trajectory$population_average_p_misseg_median
  )
  trajectory$log10_p_q10 <- log10(trajectory$population_average_p_misseg_q10)
  trajectory$log10_p_q90 <- log10(trajectory$population_average_p_misseg_q90)

  fill_limits <- range(
    surface$dominant_mean_ploidy_median, c(1, 4), na.rm = TRUE
  )
  if (fill_limits[[1L]] <= 0) {
    stop("Figure 6A log-scaled ploidy fill requires strictly positive values.")
  }
  fill_breaks <- c(1, 1.5, 2, 3, 4, 6)
  fill_breaks <- fill_breaks[
    fill_breaks >= fill_limits[[1L]] & fill_breaks <= fill_limits[[2L]]
  ]
  y_breaks <- pretty(range(surface$log10_effective_p_misseg), n = 5)
  y_labels <- function(x) formatC(10^x, format = "e", digits = 0)

  plots <- lapply(seq_along(ordered_pairs), function(i) {
    pair <- ordered_pairs[[i]]
    s <- surface[surface$pair_id == pair, , drop = FALSE]
    t <- trajectory[trajectory$pair_id == pair, , drop = FALSE]
    hatch <- f6r_weak_gap_hatch_data(s)
    low_consensus <- s[
      s$ploidy_regime_consensus < 0.80 &
        (match(s$O2_pct, sort(unique(s$O2_pct))) %% 8L == 1L) &
        (match(s$effective_p_misseg,
               sort(unique(s$effective_p_misseg))) %% 3L == 1L),
      , drop = FALSE
    ]
    p <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = s,
        ggplot2::aes(
          x = O2_pct, y = log10_effective_p_misseg,
          fill = dominant_mean_ploidy_median
        )
      ) +
      ggplot2::geom_path(
        data = hatch,
        ggplot2::aes(
          x = O2_pct, y = log10_effective_p_misseg,
          group = hatch_group
        ),
        inherit.aes = FALSE, colour = "#9B59B6",
        linewidth = 0.20, alpha = 0.72,
        lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_contour(
        data = s,
        ggplot2::aes(
          x = O2_pct, y = log10_effective_p_misseg,
          z = proportion_spectral_gap_below_0p005,
          colour = "Weak spectral-gap boundary",
          linetype = "Weak spectral-gap boundary"
        ),
        breaks = 0.5, linewidth = 0.38, show.legend = TRUE
      ) +
      ggplot2::geom_point(
        data = low_consensus,
        ggplot2::aes(x = O2_pct, y = log10_effective_p_misseg),
        shape = 4, size = 0.42, stroke = 0.28,
        colour = "white", show.legend = FALSE
      ) +
      ggplot2::geom_ribbon(
        data = t,
        ggplot2::aes(x = O2_pct, ymin = log10_p_q10, ymax = log10_p_q90),
        fill = "#E6E6E6", alpha = 0.65, inherit.aes = FALSE
      ) +
      ggplot2::geom_path(
        data = t,
        ggplot2::aes(
          x = O2_pct, y = log10_p_median,
          colour = "Model-implied median trajectory",
          linetype = "Model-implied median trajectory"
        ),
        linewidth = 0.58
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          "Model-implied median trajectory" = "#111111",
          "Weak spectral-gap boundary" = "#9B59B6"
        ),
        breaks = c(
          "Model-implied median trajectory",
          "Weak spectral-gap boundary"
        ),
        labels = c(
          "Median fitted\ntrajectory",
          "Weak-gap region\n(hatched)"
        ),
        name = "Overlays"
      ) +
      ggplot2::scale_linetype_manual(
        values = c(
          "Model-implied median trajectory" = "solid",
          "Weak spectral-gap boundary" = "dotted"
        ),
        breaks = c(
          "Model-implied median trajectory",
          "Weak spectral-gap boundary"
        ),
        labels = c(
          "Median fitted\ntrajectory",
          "Weak-gap region\n(hatched)"
        ),
        name = "Overlays"
      ) +
      ggplot2::scale_fill_gradientn(
        colours = c("#2C7BB6", "#74ADD1", "#FEE08B", "#F46D43", "#A50026"),
        limits = fill_limits, breaks = fill_breaks,
        labels = scales::label_number(accuracy = 0.1),
        trans = "log10", oob = scales::squish,
        name = "Median dominant\nploidy\n(log colors)"
      ) +
      ggplot2::scale_x_continuous(breaks = 0:5, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(
        breaks = y_breaks, labels = y_labels, expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(
        ylim = range(s$log10_effective_p_misseg), expand = FALSE
      ) +
      ggplot2::labs(
        title = names(display_pair_labels)[[i]],
        x = "Fixed oxygen (%)",
        y = "Effective per-chromosome\nmissegregation probability"
      ) +
      f6r_theme(8) +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "right",
        legend.box = "vertical",
        legend.title = ggplot2::element_text(size = 7.2, lineheight = 0.9),
        legend.text = ggplot2::element_text(size = 6.6, lineheight = 0.9),
        legend.key.height = grid::unit(2.8, "mm"),
        legend.spacing.y = grid::unit(0.7, "mm"),
        legend.box.spacing = grid::unit(1.0, "mm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title.position = "panel",
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_colourbar(
          order = 1,
          barheight = grid::unit(22, "mm"),
          barwidth = grid::unit(3.2, "mm")
        ),
        colour = ggplot2::guide_legend(
          order = 2, override.aes = list(linewidth = 0.7),
          keywidth = grid::unit(8, "mm"),
          keyheight = grid::unit(3.2, "mm")
        ),
        linetype = ggplot2::guide_legend(order = 2)
      )
    col <- i
    if (col != 1L) p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
    p
  })
  composite <- f6r_compose_three_panel_row(
    plots,
    title = "A. Joint-fit oxygen-CIN-ploidy surfaces",
    legend_rel_width = 0.72
  )
  output <- f6r_save_plot(
    composite,
    file.path(
      paths$figure6, "panels",
      "pair_surface_o2_effective_p_misseg_primary_family_grid"
    ),
    width = f6r_primary_row_width(), height = 3.55
  )
  invisible(list(plot = composite, paths = output))
}

f6r_bilinear_grid_value <- function(x_grid, y_grid, z, x, y) {
  if (nrow(z) != length(y_grid) || ncol(z) != length(x_grid)) {
    stop("Bilinear interpolation grid dimensions do not match coordinates.")
  }
  ix <- findInterval(x, x_grid, all.inside = TRUE)
  iy <- findInterval(y, y_grid, all.inside = TRUE)
  ix <- pmin(pmax(ix, 1L), length(x_grid) - 1L)
  iy <- pmin(pmax(iy, 1L), length(y_grid) - 1L)
  tx <- (x - x_grid[ix]) / (x_grid[ix + 1L] - x_grid[ix])
  ty <- (y - y_grid[iy]) / (y_grid[iy + 1L] - y_grid[iy])
  z00 <- z[cbind(iy, ix)]
  z10 <- z[cbind(iy, ix + 1L)]
  z01 <- z[cbind(iy + 1L, ix)]
  z11 <- z[cbind(iy + 1L, ix + 1L)]
  (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 +
    (1 - tx) * ty * z01 + tx * ty * z11
}

f6r_map_panel_c_gap_contour_to_panel_d <- function(
    paths, curve_data, display_manifest
) {
  surface <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  surface <- surface[
    surface$cutoff == "q10" &
      surface$pair_id %in% display_manifest$pair_id,
    , drop = FALSE
  ]
  boundaries <- list()
  boundary_index <- 0L
  for (pair in display_manifest$pair_id) {
    s <- surface[surface$pair_id == pair, , drop = FALSE]
    d <- curve_data[curve_data$pair_id == pair, , drop = FALSE]
    surface_x <- sort(unique(s$O2_pct))
    surface_y <- sort(unique(log10(s$effective_p_misseg)))
    s <- s[order(log10(s$effective_p_misseg), s$O2_pct), , drop = FALSE]
    gap_matrix <- matrix(
      s$proportion_spectral_gap_below_0p005,
      nrow = length(surface_y), ncol = length(surface_x), byrow = TRUE
    )
    contours <- grDevices::contourLines(
      x = surface_x, y = surface_y, z = gap_matrix, levels = 0.5
    )

    dense_x <- sort(unique(d$O2_pct))
    dense_y <- sort(unique(log10(d$effective_p_misseg)))
    d <- d[order(log10(d$effective_p_misseg), d$O2_pct), , drop = FALSE]
    ploidy_matrix <- matrix(
      d$dominant_mean_ploidy_median,
      nrow = length(dense_y), ncol = length(dense_x), byrow = TRUE
    )
    for (j in seq_along(contours)) {
      boundary_index <- boundary_index + 1L
      contour <- contours[[j]]
      boundaries[[boundary_index]] <- data.frame(
        pair_id = pair,
        pair_label = d$pair_label[[1L]],
        display_label = d$display_label[[1L]],
        boundary_id = paste0(pair, "_gap_boundary_", sprintf("%03d", j)),
        O2_pct = contour$x,
        effective_p_misseg = 10^contour$y,
        dominant_mean_ploidy_median = f6r_bilinear_grid_value(
          dense_x, dense_y, ploidy_matrix, contour$x, contour$y
        ),
        threshold = 0.5,
        criterion = "proportion_spectral_gap_below_0p005",
        source_grid = "Figure6C_q10_201x60",
        stringsAsFactors = FALSE
      )
    }
  }
  boundary <- do.call(rbind, boundaries)
  rownames(boundary) <- NULL
  boundary
}

f6r_panel_d_data <- function(paths) {
  curve_path <- file.path(
    paths$figure6, "figure6d_fixed_p_curve_family.tsv"
  )
  dense_validation_path <- file.path(
    paths$figure6, "figure6d_dense_model_validation.tsv"
  )
  f6r_require_files(
    c(curve_path, dense_validation_path),
    "exact Figure 6D dense-grid data"
  )
  curve_data <- f6r_read_tsv(curve_path)
  if (length(unique(curve_data$pair_id)) != f6r_family_count() ||
      any(table(curve_data$pair_id) != 201L * 496L) ||
      any(curve_data$n_seed != 50L)) {
    stop("Primary Figure 6D source must contain a complete 201 x 496 grid and 50 endpoints per displayed pair.")
  }

  display_manifest <- f6r_display_pair_manifest(curve_data$pair_id, "D")
  display_manifest_path <- f6r_write_tsv(
    display_manifest,
    file.path(paths$figure6, "figure6d_displayed_pairs.tsv")
  )
  curve_data$display_label <- display_manifest$display_label[
    match(curve_data$pair_id, display_manifest$pair_id)
  ]
  curve_data <- curve_data[order(
    match(curve_data$display_label, display_manifest$display_label),
    curve_data$effective_p_misseg, curve_data$O2_pct
  ), , drop = FALSE]
  rownames(curve_data) <- NULL

  highlighted_p <- c(0.01, 0.10, 0.20, 0.30)
  highlighted <- curve_data[
    curve_data$effective_p_misseg %in% highlighted_p,
    , drop = FALSE
  ]
  highlighted$highlight_label <- factor(
    sprintf("%.2f", highlighted$effective_p_misseg),
    levels = sprintf("%.2f", highlighted_p)
  )
  highlighted_path <- f6r_write_tsv(
    highlighted,
    file.path(paths$figure6, "figure6d_highlighted_trajectories.tsv")
  )
  mean_curve <- stats::aggregate(
    curve_data$dominant_mean_ploidy_median,
    by = list(
      pair_id = curve_data$pair_id,
      pair_label = curve_data$pair_label,
      display_label = curve_data$display_label,
      O2_pct = curve_data$O2_pct
    ),
    FUN = mean
  )
  names(mean_curve)[names(mean_curve) == "x"] <-
    "dominant_mean_ploidy_mean_across_fixed_p"
  mean_curve$n_fixed_p <- 496L
  mean_curve <- mean_curve[order(
    match(mean_curve$display_label, display_manifest$display_label),
    mean_curve$O2_pct
  ), , drop = FALSE]
  rownames(mean_curve) <- NULL
  mean_curve_path <- f6r_write_tsv(
    mean_curve,
    file.path(paths$figure6, "figure6b_mean_ploidy_across_fixed_p.tsv")
  )
  gap_boundary <- f6r_map_panel_c_gap_contour_to_panel_d(
    paths, curve_data, display_manifest
  )
  gap_boundary_path <- f6r_write_tsv(
    gap_boundary,
    file.path(paths$figure6, "figure6d_spectral_gap_boundary.tsv")
  )

  p_counts <- table(factor(
    curve_data$pair_label,
    levels = f6r_family_levels()
  )) / 201L
  o2_counts <- table(interaction(
    curve_data$pair_id, curve_data$effective_p_misseg, drop = TRUE
  ))
  grid_key <- paste(
    curve_data$pair_id, curve_data$effective_p_misseg, curve_data$O2_pct,
    sep = "|"
  )
  p_range <- range(curve_data$effective_p_misseg)
  p_interval <- unique(round(diff(sort(unique(
    curve_data$effective_p_misseg
  ))), 12))
  max_forcing_error <- max(
    curve_data$max_abs_actual_minus_requested_p_misseg, na.rm = TRUE
  )
  validation <- data.frame(
    check = c(
      "displayed_pair_order", "curve_family_row_count",
      "fixed_p_count_per_pair", "oxygen_count_per_fixed_p_curve",
      "endpoint_count_per_grid_point", "grid_rows_unique",
      "all_numeric_outputs_finite", "fixed_p_positive",
      "fixed_p_range", "fixed_p_interval",
      "forcing_error_within_tolerance",
      "dense_model_validation_passed",
      "highlighted_fixed_p_values",
      "highlighted_oxygen_count_per_curve",
      "mean_curve_row_count", "mean_curve_pair_order",
      "mean_curve_oxygen_count_per_pair", "mean_curve_fixed_p_count",
      "mean_curve_values_finite",
      "spectral_gap_boundary_all_pairs",
      "spectral_gap_boundary_values_finite",
      "spectral_gap_boundary_criterion"
    ),
    observed = c(
      paste(display_manifest$pair_label, collapse = ","), nrow(curve_data),
      paste(as.integer(p_counts), collapse = ","),
      paste(sort(unique(as.integer(o2_counts))), collapse = ","),
      paste(sort(unique(curve_data$n_seed)), collapse = ","),
      !anyDuplicated(grid_key),
      all(is.finite(curve_data$dominant_mean_ploidy_median)) &&
        all(is.finite(curve_data$effective_p_misseg)) &&
        all(is.finite(curve_data$spectral_gap_median)),
      all(curve_data$effective_p_misseg > 0),
      paste(format(p_range, scientific = FALSE, trim = TRUE), collapse = ","),
      paste(format(p_interval, scientific = FALSE, trim = TRUE), collapse = ","),
      max_forcing_error,
      all(f6r_read_tsv(dense_validation_path)$passed),
      paste(sort(unique(highlighted$effective_p_misseg)), collapse = ","),
      paste(sort(unique(as.integer(table(interaction(
        highlighted$pair_id, highlighted$effective_p_misseg, drop = TRUE
      ))))), collapse = ","),
      nrow(mean_curve), paste(unique(mean_curve$pair_label), collapse = ","),
      paste(sort(unique(as.integer(table(mean_curve$pair_id)))), collapse = ","),
      paste(sort(unique(mean_curve$n_fixed_p)), collapse = ","),
      all(is.finite(mean_curve$dominant_mean_ploidy_mean_across_fixed_p)),
      paste(unique(gap_boundary$pair_label), collapse = ","),
      all(is.finite(gap_boundary$O2_pct)) &&
        all(is.finite(gap_boundary$effective_p_misseg)) &&
        all(is.finite(gap_boundary$dominant_mean_ploidy_median)),
      paste(unique(gap_boundary$criterion), collapse = ",")
    ),
    expected = c(
      paste(f6r_family_levels(), collapse = ","),
      as.character(f6r_family_count() * 201L * 496L),
      paste(rep(496L, f6r_family_count()), collapse = ","), "201",
      "50", "TRUE", "TRUE", "TRUE", "0.005,0.5", "0.001",
      "<=1e-8", "TRUE", "0.01,0.1,0.2,0.3", "201",
      as.character(f6r_family_count() * 201L),
      paste(f6r_family_levels(), collapse = ","),
      "201", "496", "TRUE", paste(f6r_family_levels(), collapse = ","), "TRUE",
      "proportion_spectral_gap_below_0p005"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- c(
    identical(display_manifest$pair_label, f6r_family_levels()),
    nrow(curve_data) == f6r_family_count() * 201L * 496L,
    all(p_counts == 496L),
    all(o2_counts == 201L), all(curve_data$n_seed == 50L),
    !anyDuplicated(grid_key),
    all(is.finite(curve_data$dominant_mean_ploidy_median)) &&
      all(is.finite(curve_data$effective_p_misseg)) &&
      all(is.finite(curve_data$spectral_gap_median)),
    all(curve_data$effective_p_misseg > 0),
    isTRUE(all.equal(unname(p_range), c(0.005, 0.5), tolerance = 1e-12)),
    isTRUE(all.equal(p_interval, 0.001, tolerance = 1e-12)),
    is.finite(max_forcing_error) && max_forcing_error <= 1e-8,
    all(f6r_read_tsv(dense_validation_path)$passed),
    identical(sort(unique(highlighted$effective_p_misseg)), highlighted_p),
    nrow(highlighted) == 4L * f6r_family_count() * 201L &&
      all(table(interaction(
        highlighted$pair_id, highlighted$effective_p_misseg, drop = TRUE
      )) == 201L),
    nrow(mean_curve) == f6r_family_count() * 201L,
    identical(unique(mean_curve$pair_label),
              f6r_family_levels()),
    all(table(mean_curve$pair_id) == 201L),
    all(mean_curve$n_fixed_p == 496L),
    all(is.finite(mean_curve$dominant_mean_ploidy_mean_across_fixed_p)),
    identical(unique(gap_boundary$pair_label),
              f6r_family_levels()),
    nrow(gap_boundary) > 0L &&
      all(is.finite(gap_boundary$O2_pct)) &&
      all(is.finite(gap_boundary$effective_p_misseg)) &&
      all(is.finite(gap_boundary$dominant_mean_ploidy_median)),
    identical(
      unique(gap_boundary$criterion),
      "proportion_spectral_gap_below_0p005"
    )
  )
  validation_path <- f6r_write_tsv(
    validation,
    file.path(paths$figure6, "figure6d_validation.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Figure 6D data validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(list(
    curve_data = curve_data, highlighted = highlighted,
    mean_curve = mean_curve,
    gap_boundary = gap_boundary,
    manifest = display_manifest,
    paths = c(
      displayed_pairs = display_manifest_path,
      curve_family = curve_path,
      highlighted_trajectories = highlighted_path,
      mean_ploidy_across_fixed_p = mean_curve_path,
      spectral_gap_boundary = gap_boundary_path,
      dense_model_validation = dense_validation_path,
      validation = validation_path
    )
  ))
}

f6r_draw_fixed_p_curve_panel <- function(paths) {
  f6r_require_packages(c("ggplot2", "cowplot", "egg", "scales"))
  bundle <- f6r_panel_d_data(paths)
  curve_data <- bundle$curve_data
  highlighted <- bundle$highlighted
  mean_curve <- bundle$mean_curve
  display_manifest <- bundle$manifest
  reference_levels <- c(
    "0.01", "0.10", "0.20", "0.30",
    "Mean across 496 fixed p_miss,eff values"
  )
  mean_curve$reference_label <- factor(
    "Mean across 496 fixed p_miss,eff values",
    levels = reference_levels
  )
  color_limits <- range(curve_data$effective_p_misseg, finite = TRUE)
  if (any(!is.finite(color_limits)) || color_limits[[1L]] <= 0) {
    stop("Figure 6D log-scaled colors require positive finite values.")
  }
  color_breaks <- scales::breaks_log(n = 5)(color_limits)
  color_breaks <- color_breaks[
    color_breaks >= color_limits[[1L]] & color_breaks <= color_limits[[2L]]
  ]
  y_limits <- c(1, 7)

  plots <- lapply(seq_len(nrow(display_manifest)), function(i) {
    pair <- display_manifest$pair_id[[i]]
    d <- curve_data[curve_data$pair_id == pair, , drop = FALSE]
    d <- d[order(d$effective_p_misseg, d$O2_pct), , drop = FALSE]
    h <- highlighted[highlighted$pair_id == pair, , drop = FALSE]
    h <- h[order(h$effective_p_misseg, h$O2_pct), , drop = FALSE]
    m <- mean_curve[mean_curve$pair_id == pair, , drop = FALSE]
    m <- m[order(m$O2_pct), , drop = FALSE]
    p <- ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = O2_pct, y = dominant_mean_ploidy_median,
        group = effective_p_misseg, colour = effective_p_misseg
      )
    ) +
      ggplot2::geom_path(
        linewidth = 0.18, alpha = 0.72, lineend = "round"
      ) +
      ggplot2::geom_path(
        data = h,
        ggplot2::aes(linetype = highlight_label),
        colour = "#111111", linewidth = 0.62, alpha = 1,
        lineend = "round"
      ) +
      ggplot2::geom_path(
        data = m,
        ggplot2::aes(
          x = O2_pct,
          y = dominant_mean_ploidy_mean_across_fixed_p,
          linetype = reference_label
        ),
        inherit.aes = FALSE,
        colour = "#D62728", linewidth = 0.82, alpha = 1,
        lineend = "round"
      ) +
      ggplot2::scale_colour_viridis_c(
        option = "D", begin = 0.05, end = 0.95,
        limits = color_limits, breaks = color_breaks, trans = "log10",
        labels = function(x) formatC(x, format = "e", digits = 1),
        name = paste0(
          "Fixed p_miss,eff\n",
          "(log-scaled colors)"
        )
      ) +
      ggplot2::scale_linetype_manual(
        values = c(
          "0.01" = "solid", "0.10" = "F28282",
          "0.20" = "dotdash", "0.30" = "dotted",
          "Mean across 496 fixed p_miss,eff values" = "solid"
        ),
        breaks = reference_levels,
        labels = c(
          "p_miss,eff = 0.01", "p_miss,eff = 0.10",
          "p_miss,eff = 0.20", "p_miss,eff = 0.30",
          "Mean across 496 fixed\np_miss,eff values"
        ),
        name = "Reference curves"
      ) +
      ggplot2::scale_x_continuous(
        breaks = 0:5, limits = c(0, 5), expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(ylim = y_limits, expand = FALSE) +
      ggplot2::labs(
        title = display_manifest$display_label[[i]],
        x = "Fixed oxygen (%)", y = "Dominant mean ploidy"
      ) +
      f6r_theme(8) +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "right",
        legend.box = "vertical",
        legend.title = ggplot2::element_text(size = 7.2, lineheight = 0.9),
        legend.text = ggplot2::element_text(size = 6.6, lineheight = 0.9),
        legend.key.height = grid::unit(2.8, "mm"),
        legend.spacing.y = grid::unit(0.7, "mm"),
        legend.box.spacing = grid::unit(1.0, "mm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title.position = "panel",
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_colourbar(
          barheight = grid::unit(22, "mm"),
          barwidth = grid::unit(3.2, "mm"), order = 1
        ),
        linetype = ggplot2::guide_legend(
          keywidth = grid::unit(8, "mm"),
          keyheight = grid::unit(3.2, "mm"), order = 2,
          override.aes = list(
            colour = c(rep("#111111", 4L), "#D62728"),
            linewidth = c(rep(0.62, 4L), 0.82),
            alpha = 1
          )
        )
      )
    if (i != 1L) {
      p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
    p
  })
  composite <- f6r_compose_three_panel_row(
    plots,
    title = "B. Ploidy responses at fixed missegregation probability",
    legend_rel_width = 0.72
  )
  output <- f6r_save_plot(
    composite,
    file.path(
      paths$figure6, "panels",
      "pair_fixed_p_miss_eff_o2_ploidy_curve_primary_family_grid"
    ),
    width = f6r_primary_row_width(), height = 3.55
  )
  invisible(list(plot = composite, paths = output, data = bundle))
}

f6r_inverse_curve_solutions <- function(
    effective_p_misseg, dominant_mean_ploidy, target_ploidy,
    numerical_tolerance = 1e-10
) {
  p <- as.numeric(effective_p_misseg)
  y <- as.numeric(dominant_mean_ploidy)
  target <- as.numeric(target_ploidy)
  keep <- is.finite(p) & is.finite(y)
  p <- p[keep]
  y <- y[keep]
  if (length(p) < 2L || any(!is.finite(target))) {
    stop("Inverse response requires at least two finite forward points and finite targets.")
  }
  ord <- order(p)
  p <- p[ord]
  y <- y[ord]
  if (anyDuplicated(p) || any(diff(p) <= 0)) {
    stop("Effective-missegregation inputs must be strictly increasing for inversion.")
  }

  y_scale <- max(1, diff(range(y, finite = TRUE)))
  tol <- max(as.numeric(numerical_tolerance), 1e-10 * y_scale)
  dy <- diff(y)
  raw_sign <- ifelse(abs(dy) <= tol, 0L, ifelse(dy > 0, 1L, -1L))
  flat_segment <- raw_sign == 0L
  work_sign <- raw_sign
  if (all(work_sign == 0L)) {
    hit <- abs(target - mean(y)) <= tol
    return(data.frame(
      target_ploidy = target,
      n_solution = ifelse(hit, 2L, 0L),
      p_unique = NA_real_,
      p_solution_min = ifelse(hit, min(p), NA_real_),
      p_solution_max = ifelse(hit, max(p), NA_real_),
      forward_reconstruction_error = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  for (i in seq_along(work_sign)) {
    if (work_sign[[i]] == 0L && i > 1L) work_sign[[i]] <- work_sign[[i - 1L]]
  }
  for (i in rev(seq_along(work_sign))) {
    if (work_sign[[i]] == 0L && i < length(work_sign)) {
      work_sign[[i]] <- work_sign[[i + 1L]]
    }
  }

  change <- which(diff(work_sign) != 0L)
  run_start <- c(1L, change + 1L)
  run_end <- c(change + 1L, length(y))
  solutions <- matrix(
    NA_real_, nrow = length(target), ncol = length(run_start)
  )
  for (run_index in seq_along(run_start)) {
    idx <- run_start[[run_index]]:run_end[[run_index]]
    yr <- y[idx]
    pr <- p[idx]
    run_order <- order(yr, pr)
    yr <- yr[run_order]
    pr <- pr[run_order]
    unique_y <- !duplicated(yr)
    yr <- yr[unique_y]
    pr <- pr[unique_y]
    if (length(yr) >= 2L && diff(range(yr)) > tol) {
      solutions[, run_index] <- stats::approx(
        x = yr, y = pr, xout = target,
        rule = 1, ties = "ordered"
      )$y
    }
  }

  if (ncol(solutions) > 1L) {
    for (later in 2:ncol(solutions)) {
      for (earlier in seq_len(later - 1L)) {
        duplicate_solution <- is.finite(solutions[, later]) &
          is.finite(solutions[, earlier]) &
          abs(solutions[, later] - solutions[, earlier]) <= 1e-10
        solutions[duplicate_solution, later] <- NA_real_
      }
    }
  }
  n_solution <- rowSums(is.finite(solutions))
  p_min <- rep(Inf, length(target))
  p_max <- rep(-Inf, length(target))
  for (column in seq_len(ncol(solutions))) {
    p_min <- pmin(p_min, solutions[, column], na.rm = TRUE)
    p_max <- pmax(p_max, solutions[, column], na.rm = TRUE)
  }

  if (any(flat_segment)) {
    flat_indices <- which(flat_segment)
    for (segment_index in flat_indices) {
      hit <- abs(target - y[[segment_index]]) <= tol
      if (any(hit)) {
        n_solution[hit] <- pmax(n_solution[hit], 2L)
        p_min[hit] <- pmin(p_min[hit], p[[segment_index]], na.rm = TRUE)
        p_max[hit] <- pmax(p_max[hit], p[[segment_index + 1L]], na.rm = TRUE)
      }
    }
  }
  p_min[!is.finite(p_min)] <- NA_real_
  p_max[!is.finite(p_max)] <- NA_real_
  p_unique <- ifelse(n_solution == 1L, p_min, NA_real_)
  reconstruction_error <- rep(NA_real_, length(target))
  unique_index <- which(is.finite(p_unique))
  if (length(unique_index)) {
    reconstructed <- stats::approx(
      x = p, y = y, xout = p_unique[unique_index],
      rule = 1, ties = "ordered"
    )$y
    reconstruction_error[unique_index] <- abs(
      reconstructed - target[unique_index]
    )
  }
  data.frame(
    target_ploidy = target,
    n_solution = as.integer(n_solution),
    p_unique = p_unique,
    p_solution_min = p_min,
    p_solution_max = p_max,
    forward_reconstruction_error = reconstruction_error,
    stringsAsFactors = FALSE
  )
}

f6r_inverse_endpoint_cache <- function(
    dense_cache_path, inverse_cache_path, target_ploidy,
    force_rebuild = FALSE
) {
  profile <- "inverse_target_ploidy_1to7_step0p025_v1"
  dense_md5 <- f6r_md5(dense_cache_path)
  if (!force_rebuild && file.exists(inverse_cache_path)) {
    cached <- readRDS(inverse_cache_path)
    if (identical(cached$metadata$profile, profile) &&
        identical(cached$metadata$dense_source_md5, dense_md5) &&
        identical(as.numeric(cached$metadata$target_ploidy),
                  as.numeric(target_ploidy))) {
      return(cached$qc)
    }
  }

  dense <- readRDS(dense_cache_path)
  surface <- dense$surface
  required <- c("O2_pct", "effective_p_misseg", "dominant_mean_ploidy")
  if (!all(required %in% names(surface))) {
    stop("Dense Figure 6 endpoint cache is missing inverse-response fields: ",
         dense_cache_path)
  }
  oxygen_values <- sort(unique(surface$O2_pct))
  inverse_rows <- lapply(oxygen_values, function(o2) {
    current <- surface[surface$O2_pct == o2, , drop = FALSE]
    current <- current[order(current$effective_p_misseg), , drop = FALSE]
    result <- f6r_inverse_curve_solutions(
      effective_p_misseg = current$effective_p_misseg,
      dominant_mean_ploidy = current$dominant_mean_ploidy,
      target_ploidy = target_ploidy
    )
    result$O2_pct <- o2
    result[, c(
      "O2_pct", "target_ploidy", "n_solution", "p_unique",
      "p_solution_min", "p_solution_max",
      "forward_reconstruction_error"
    ), drop = FALSE]
  })
  inverse <- do.call(rbind, inverse_rows)
  rownames(inverse) <- NULL
  metadata <- c(
    dense$metadata,
    list(
      profile = profile,
      dense_source_path = normalizePath(dense_cache_path, mustWork = TRUE),
      dense_source_md5 = dense_md5,
      target_ploidy = as.numeric(target_ploidy)
    )
  )
  finite_error <- inverse$forward_reconstruction_error[
    is.finite(inverse$forward_reconstruction_error)
  ]
  max_error <- if (length(finite_error)) max(finite_error) else NA_real_
  qc <- data.frame(
    display_label = metadata$display_label,
    pair_label = metadata$pair_label,
    pair_id = metadata$pair_id,
    parameter_endpoint_group = metadata$parameter_endpoint_group,
    representative_seed_number = metadata$representative_seed_number,
    endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10,
    n_o2 = length(oxygen_values),
    n_target_ploidy = length(target_ploidy),
    n_inverse_grid = nrow(inverse),
    n_grid_with_any_solution = sum(inverse$n_solution >= 1L),
    n_grid_with_unique_solution = sum(inverse$n_solution == 1L),
    n_grid_with_multiple_solutions = sum(inverse$n_solution >= 2L),
    max_forward_reconstruction_error = max_error,
    inverse_qc_pass = nrow(inverse) == length(oxygen_values) * length(target_ploidy) &&
      all(inverse$n_solution >= 0L) &&
      all(!is.finite(inverse$p_unique) | inverse$n_solution == 1L) &&
      (is.na(max_error) || max_error <= 1e-8),
    cache_path = normalizePath(inverse_cache_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(inverse_cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(metadata = metadata, inverse = inverse, qc = qc),
    inverse_cache_path, compress = "xz"
  )
  qc$cache_path <- normalizePath(inverse_cache_path, mustWork = TRUE)
  saved <- readRDS(inverse_cache_path)
  saved$qc <- qc
  saveRDS(saved, inverse_cache_path, compress = "xz")
  qc
}

f6r_weighted_empirical_quantile <- function(values, weights, probs) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[keep]
  weights <- weights[keep]
  if (!length(values)) return(rep(NA_real_, length(probs)))
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord]
  cumulative <- cumsum(weights)
  total <- sum(weights)
  vapply(probs, function(probability) {
    values[[which(cumulative >= probability * total)[[1L]]]]
  }, numeric(1L))
}

f6r_inverse_panel_data <- function(
    paths, rebuild = FALSE, n_core = 1L,
    dense_qc_path = file.path(
      paths$figure6, "figure6d_dense_endpoint_qc.tsv"
    ),
    output_prefix = "figure6", model_context = "in vivo"
) {
  f6r_require_packages("data.table")
  f6r_require_files(dense_qc_path, "Figure 6 dense endpoint QC")
  dense_qc <- f6r_read_tsv(dense_qc_path)
  dense_qc <- dense_qc[
    dense_qc$pair_label %in% f6r_family_levels(),
    , drop = FALSE
  ]
  if (nrow(dense_qc) < f6r_family_count() ||
      !all(dense_qc$operator_qc_pass) ||
      any(!file.exists(dense_qc$cache_path)) ||
      !identical(
        as.integer(tapply(
          dense_qc$endpoint_multiplicity_q10,
          factor(dense_qc$pair_label, levels = f6r_family_levels()), sum
        )),
        rep(50L, f6r_family_count())
      )) {
    stop("Figure 6B inversion requires validated dense caches representing 50 endpoints in each primary family.")
  }
  target_ploidy <- seq(1, 7, by = 0.025)
  inverse_root <- file.path(
    paths$figure6, paste0(output_prefix, "_inverse_endpoint_cache")
  )
  inverse_paths <- file.path(
    inverse_root, dense_qc$pair_label,
    paste0("inverse_", dense_qc$parameter_endpoint_group, ".rds")
  )
  compute_one <- function(index) {
    message(
      "Figure 6B inverse endpoint: ", dense_qc$pair_label[[index]], " ",
      dense_qc$parameter_endpoint_group[[index]], " (", index, "/",
      nrow(dense_qc), ")"
    )
    f6r_inverse_endpoint_cache(
      dense_cache_path = dense_qc$cache_path[[index]],
      inverse_cache_path = inverse_paths[[index]],
      target_ploidy = target_ploidy,
      force_rebuild = rebuild
    )
  }
  index <- seq_len(nrow(dense_qc))
  qc_rows <- f6r_resilient_lapply(index, compute_one, n_core = n_core)
  if (any(vapply(qc_rows, function(x) inherits(x, "try-error"), logical(1L)))) {
    stop("One or more Figure 6B inverse endpoint workers failed.")
  }
  inverse_qc <- do.call(rbind, qc_rows)
  inverse_qc$model_context <- model_context
  if (!all(inverse_qc$inverse_qc_pass)) {
    stop("Figure 6B inverse endpoint QC failed.")
  }
  inverse_qc_path <- f6r_write_tsv(
    inverse_qc,
    file.path(paths$figure6, paste0(output_prefix, "_inverse_endpoint_qc.tsv"))
  )

  pair_order <- f6r_family_levels()
  pair_summaries <- lapply(pair_order, function(pair_label) {
    pair_qc <- inverse_qc[inverse_qc$pair_label == pair_label, , drop = FALSE]
    endpoint_tables <- lapply(seq_len(nrow(pair_qc)), function(index) {
      cached <- readRDS(pair_qc$cache_path[[index]])
      d <- data.table::as.data.table(cached$inverse)
      d[, `:=`(
        pair_label = pair_label,
        pair_id = pair_qc$pair_id[[index]],
        display_label = pair_qc$display_label[[index]],
        parameter_endpoint_group = pair_qc$parameter_endpoint_group[[index]],
        endpoint_multiplicity_q10 = pair_qc$endpoint_multiplicity_q10[[index]]
      )]
      d
    })
    endpoint_data <- data.table::rbindlist(endpoint_tables, use.names = TRUE)
    summary <- endpoint_data[, {
      multiplicity <- as.integer(endpoint_multiplicity_q10)
      any_solution <- n_solution >= 1L
      unique_solution <- n_solution == 1L
      multiple_solution <- n_solution >= 2L
      quantiles <- f6r_weighted_empirical_quantile(
        p_unique[unique_solution], multiplicity[unique_solution],
        probs = c(0.10, 0.50, 0.90)
      )
      finite_min <- p_solution_min[any_solution & is.finite(p_solution_min)]
      finite_max <- p_solution_max[any_solution & is.finite(p_solution_max)]
      finite_error <- forward_reconstruction_error[
        unique_solution & is.finite(forward_reconstruction_error)
      ]
      n_seed <- sum(multiplicity)
      list(
        n_unique_parameter_endpoint = .N,
        n_seed = n_seed,
        n_seed_any_solution = sum(multiplicity[any_solution]),
        n_seed_unique_solution = sum(multiplicity[unique_solution]),
        n_seed_multiple_solutions = sum(multiplicity[multiple_solution]),
        fraction_any_solution = sum(multiplicity[any_solution]) / n_seed,
        fraction_unique_solution = sum(multiplicity[unique_solution]) / n_seed,
        fraction_multiple_solutions = sum(multiplicity[multiple_solution]) / n_seed,
        p_unique_q10 = quantiles[[1L]],
        p_unique_median = quantiles[[2L]],
        p_unique_q90 = quantiles[[3L]],
        p_solution_min = if (length(finite_min)) min(finite_min) else NA_real_,
        p_solution_max = if (length(finite_max)) max(finite_max) else NA_real_,
        max_forward_reconstruction_error = if (length(finite_error)) {
          max(finite_error)
        } else {
          NA_real_
        }
      )
    }, by = .(pair_id, pair_label, display_label, O2_pct, target_ploidy)]
    as.data.frame(summary)
  })
  inverse_summary <- do.call(rbind, pair_summaries)
  rownames(inverse_summary) <- NULL
  inverse_summary$model_context <- model_context
  inverse_summary$inverse_class <- ifelse(
    inverse_summary$fraction_multiple_solutions >= 0.20,
    "multiple solutions",
    ifelse(
      inverse_summary$fraction_any_solution < 0.80,
      "no stable unique inverse",
      ifelse(
        inverse_summary$fraction_unique_solution >= 0.80,
        "stable unique inverse",
        "no stable unique inverse"
      )
    )
  )
  inverse_summary$p_display <- ifelse(
    inverse_summary$inverse_class == "stable unique inverse",
    inverse_summary$p_unique_median,
    NA_real_
  )
  inverse_summary <- inverse_summary[order(
    match(inverse_summary$pair_label, pair_order),
    inverse_summary$target_ploidy, inverse_summary$O2_pct
  ), , drop = FALSE]
  summary_path <- f6r_write_tsv(
    inverse_summary,
    file.path(
      paths$figure6, paste0(output_prefix, "_inverse_response_summary.tsv")
    )
  )

  class_summary <- data.table::as.data.table(inverse_summary)[, .(
    n_inverse_grid = .N,
    n_stable_unique_inverse = sum(inverse_class == "stable unique inverse"),
    n_multiple_solutions = sum(inverse_class == "multiple solutions"),
    n_no_stable_unique_inverse = sum(
      inverse_class == "no stable unique inverse"
    ),
    fraction_stable_unique_inverse = mean(
      inverse_class == "stable unique inverse"
    ),
    fraction_multiple_solutions = mean(
      inverse_class == "multiple solutions"
    ),
    fraction_no_stable_unique_inverse = mean(
      inverse_class == "no stable unique inverse"
    ),
    minimum_displayed_p = if (any(is.finite(p_display))) {
      min(p_display, na.rm = TRUE)
    } else {
      NA_real_
    },
    maximum_displayed_p = if (any(is.finite(p_display))) {
      max(p_display, na.rm = TRUE)
    } else {
      NA_real_
    }
  ), by = .(pair_id, pair_label, display_label)]
  class_summary <- as.data.frame(class_summary)
  class_summary$model_context <- model_context
  class_summary_path <- f6r_write_tsv(
    class_summary,
    file.path(
      paths$figure6, paste0(output_prefix, "_inverse_class_summary.tsv")
    )
  )
  anchor_summary <- inverse_summary[
    inverse_summary$target_ploidy == 4 &
      inverse_summary$O2_pct %in% c(0, 1, 5),
    c(
      "pair_id", "pair_label", "display_label", "O2_pct",
      "target_ploidy", "inverse_class", "n_seed",
      "fraction_any_solution", "fraction_unique_solution",
      "fraction_multiple_solutions", "p_unique_q10",
      "p_unique_median", "p_unique_q90"
    ), drop = FALSE
  ]
  anchor_summary <- anchor_summary[order(
    match(anchor_summary$pair_label, pair_order), anchor_summary$O2_pct
  ), , drop = FALSE]
  anchor_summary$model_context <- model_context
  anchor_summary_path <- f6r_write_tsv(
    anchor_summary,
    file.path(
      paths$figure6, paste0(output_prefix, "_inverse_ploidy4_anchor_summary.tsv")
    )
  )
  manuscript_class_summary_path <- f6r_write_tsv(
    class_summary,
    file.path(
      paths$root, "manuscript", "tables", "data", "fixed_o2",
      if (identical(model_context, "in vivo")) {
        "joint_inverse_grid_class_summary.tsv"
      } else {
        paste0("joint_invitro_inverse_grid_class_summary.tsv")
      }
    )
  )
  manuscript_anchor_summary_path <- f6r_write_tsv(
    anchor_summary,
    file.path(
      paths$root, "manuscript", "tables", "data", "fixed_o2",
      if (identical(model_context, "in vivo")) {
        "joint_inverse_ploidy4_anchor_summary.tsv"
      } else {
        paste0("joint_invitro_inverse_ploidy4_anchor_summary.tsv")
      }
    )
  )

  endpoint_manifest <- inverse_qc[, c(
    "display_label", "pair_label", "pair_id", "parameter_endpoint_group",
    "representative_seed_number", "endpoint_multiplicity_q10", "cache_path"
  ), drop = FALSE]
  endpoint_manifest$cache_md5 <- vapply(
    endpoint_manifest$cache_path, f6r_md5, character(1L)
  )
  endpoint_manifest_path <- f6r_write_tsv(
    endpoint_manifest,
    file.path(
      paths$figure6, paste0(output_prefix, "_inverse_endpoint_manifest.tsv")
    )
  )

  target_interval <- unique(round(diff(sort(unique(
    inverse_summary$target_ploidy
  ))), 12))
  max_reconstruction_error <- max(
    inverse_summary$max_forward_reconstruction_error, na.rm = TRUE
  )
  validation <- data.frame(
    check = c(
      "displayed_pair_order", "inverse_summary_row_count",
      "oxygen_count_per_pair_target", "target_ploidy_count_per_pair_oxygen",
      "target_ploidy_range", "target_ploidy_interval",
      "represented_seed_count_per_pair", "inverse_fraction_bounds",
      "inverse_classes_valid", "display_values_only_for_stable_unique_cells",
      "display_probability_range", "forward_reconstruction_error",
      "inverse_endpoint_qc_pass", "ploidy4_anchor_row_count"
    ),
    observed = c(
      paste(class_summary$pair_label, collapse = ","),
      nrow(inverse_summary),
      paste(sort(unique(as.integer(table(interaction(
        inverse_summary$pair_id, inverse_summary$target_ploidy, drop = TRUE
      ))))), collapse = ","),
      paste(sort(unique(as.integer(table(interaction(
        inverse_summary$pair_id, inverse_summary$O2_pct, drop = TRUE
      ))))), collapse = ","),
      paste(range(inverse_summary$target_ploidy), collapse = ","),
      paste(target_interval, collapse = ","),
      paste(tapply(
        inverse_qc$endpoint_multiplicity_q10,
        factor(inverse_qc$pair_label, levels = pair_order), sum
      ), collapse = ","),
      all(inverse_summary$fraction_any_solution >= 0 &
            inverse_summary$fraction_any_solution <= 1 &
            inverse_summary$fraction_unique_solution >= 0 &
            inverse_summary$fraction_unique_solution <= 1 &
            inverse_summary$fraction_multiple_solutions >= 0 &
            inverse_summary$fraction_multiple_solutions <= 1),
      all(inverse_summary$inverse_class %in% c(
        "stable unique inverse", "multiple solutions",
        "no stable unique inverse"
      )),
      all(is.finite(inverse_summary$p_display) ==
            (inverse_summary$inverse_class == "stable unique inverse")),
      all(!is.finite(inverse_summary$p_display) |
            (inverse_summary$p_display >= 0.005 &
               inverse_summary$p_display <= 0.500)),
      max_reconstruction_error,
      all(inverse_qc$inverse_qc_pass), nrow(anchor_summary)
    ),
    expected = c(
      paste(pair_order, collapse = ","),
      as.character(f6r_family_count() * 201L * 241L),
      "201", "241", "1,7", "0.025",
      paste(rep(50L, f6r_family_count()), collapse = ","),
      "TRUE", "TRUE", "TRUE", "TRUE", "<=1e-8", "TRUE",
      as.character(3L * f6r_family_count())
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == validation$expected
  validation$passed[validation$check == "forward_reconstruction_error"] <-
    is.finite(max_reconstruction_error) && max_reconstruction_error <= 1e-8
  validation_path <- f6r_write_tsv(
    validation,
    file.path(
      paths$figure6, paste0(output_prefix, "_inverse_validation.tsv")
    )
  )
  if (!all(validation$passed)) {
    stop(
      "Figure 6B inverse-response validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(list(
    summary = inverse_summary,
    class_summary = class_summary,
    endpoint_qc = inverse_qc,
    paths = c(
      response_summary = summary_path,
      class_summary = class_summary_path,
      ploidy4_anchor_summary = anchor_summary_path,
      manuscript_class_summary = manuscript_class_summary_path,
      manuscript_ploidy4_anchor_summary = manuscript_anchor_summary_path,
      endpoint_manifest = endpoint_manifest_path,
      endpoint_qc = inverse_qc_path,
      validation = validation_path
    )
  ))
}

f6r_inverse_multivalue_hatch_data <- function(
    inverse_summary, threshold = 0.20, spacing = 0.20
) {
  x <- sort(unique(inverse_summary$O2_pct))
  y <- sort(unique(inverse_summary$target_ploidy))
  z <- matrix(NA_real_, nrow = length(y), ncol = length(x))
  z[cbind(
    match(inverse_summary$target_ploidy, y),
    match(inverse_summary$O2_pct, x)
  )] <- inverse_summary$fraction_multiple_solutions
  if (anyNA(z)) stop("Incomplete Figure 6B inverse-response grid.")
  if (max(z, na.rm = TRUE) < threshold) {
    return(data.frame(
      O2_pct = numeric(), target_ploidy = numeric(), hatch_group = integer()
    ))
  }
  bands <- isoband::isobands(
    x, y, z,
    levels_low = threshold,
    levels_high = max(z, na.rm = TRUE) + 1e-8
  )
  geometries <- isoband::iso_to_sfg(bands)
  if (!length(geometries)) {
    return(data.frame(
      O2_pct = numeric(), target_ploidy = numeric(), hatch_group = integer()
    ))
  }
  region <- sf::st_sfc(geometries[[1L]])
  bounds <- sf::st_bbox(region)
  slope <- diff(range(y)) / diff(range(x))
  intercepts <- seq(
    bounds[["ymin"]] - slope * bounds[["xmax"]],
    bounds[["ymax"]] - slope * bounds[["xmin"]],
    by = spacing
  )
  hatch_lines <- sf::st_sfc(lapply(intercepts, function(intercept) {
    sf::st_linestring(matrix(
      c(
        bounds[["xmin"]], intercept + slope * bounds[["xmin"]],
        bounds[["xmax"]], intercept + slope * bounds[["xmax"]]
      ),
      ncol = 2L, byrow = TRUE
    ))
  }))
  clipped <- suppressWarnings(sf::st_intersection(hatch_lines, region))
  clipped <- suppressWarnings(sf::st_collection_extract(
    clipped, "LINESTRING"
  ))
  clipped <- suppressWarnings(sf::st_cast(clipped, "LINESTRING"))
  if (!length(clipped)) {
    return(data.frame(
      O2_pct = numeric(), target_ploidy = numeric(), hatch_group = integer()
    ))
  }
  coordinates <- sf::st_coordinates(clipped)
  data.frame(
    O2_pct = coordinates[, "X"],
    target_ploidy = coordinates[, "Y"],
    hatch_group = coordinates[, "L1"],
    stringsAsFactors = FALSE
  )
}

f6r_draw_inverse_response_panel <- function(paths) {
  f6r_require_packages(c(
    "ggplot2", "cowplot", "egg", "scales", "isoband", "sf", "data.table"
  ))
  inverse_bundle <- f6r_inverse_panel_data(paths)
  overlay_bundle <- f6r_panel_d_data(paths)
  inverse_summary <- inverse_bundle$summary
  highlighted <- overlay_bundle$highlighted
  mean_curve <- overlay_bundle$mean_curve
  display_manifest <- unique(inverse_summary[, c(
    "display_label", "pair_label", "pair_id"
  )])
  display_manifest <- display_manifest[match(
    f6r_family_levels(), display_manifest$display_label
  ), , drop = FALSE]
  color_limits <- c(0.005, 0.500)
  color_breaks <- c(0.005, 0.01, 0.05, 0.10, 0.50)
  status_key <- data.frame(
    O2_pct = -1, target_ploidy = -1,
    status_label = "No stable unique inverse",
    stringsAsFactors = FALSE
  )
  reference_levels <- c(
    "0.01", "0.10", "0.20", "0.30",
    "Mean across 496 fixed p_miss,eff values"
  )
  mean_curve$reference_label <- factor(
    "Mean across 496 fixed p_miss,eff values",
    levels = reference_levels
  )

  plots <- lapply(seq_len(nrow(display_manifest)), function(index) {
    pair <- display_manifest$pair_id[[index]]
    current <- inverse_summary[inverse_summary$pair_id == pair, , drop = FALSE]
    hatch <- f6r_inverse_multivalue_hatch_data(current)
    h <- highlighted[highlighted$pair_id == pair, , drop = FALSE]
    h <- h[order(h$effective_p_misseg, h$O2_pct), , drop = FALSE]
    m <- mean_curve[mean_curve$pair_id == pair, , drop = FALSE]
    m <- m[order(m$O2_pct), , drop = FALSE]
    p <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = current,
        ggplot2::aes(x = O2_pct, y = target_ploidy, fill = p_display)
      ) +
      ggplot2::geom_hline(
        yintercept = c(2, 4), colour = "#555555",
        linewidth = 0.25, linetype = "longdash"
      ) +
      ggplot2::geom_path(
        data = hatch,
        ggplot2::aes(
          x = O2_pct, y = target_ploidy, group = hatch_group
        ),
        inherit.aes = FALSE, colour = "#7B3294",
        linewidth = 0.20, alpha = 0.75,
        lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_contour(
        data = current,
        ggplot2::aes(
          x = O2_pct, y = target_ploidy,
          z = fraction_multiple_solutions,
          colour = "Multiple inverse solutions"
        ),
        breaks = 0.20, linewidth = 0.35, linetype = "dotted",
        show.legend = TRUE
      ) +
      ggplot2::geom_point(
        data = status_key,
        ggplot2::aes(
          x = O2_pct, y = target_ploidy, shape = status_label
        ),
        inherit.aes = FALSE, fill = "#EFEFEF", colour = "#888888",
        size = 2.1, stroke = 0.4, show.legend = TRUE
      ) +
      ggplot2::geom_path(
        data = h,
        ggplot2::aes(
          x = O2_pct, y = dominant_mean_ploidy_median,
          group = highlight_label, linetype = highlight_label
        ),
        inherit.aes = FALSE, colour = "#111111",
        linewidth = 0.62, alpha = 1, lineend = "round"
      ) +
      ggplot2::geom_path(
        data = m,
        ggplot2::aes(
          x = O2_pct,
          y = dominant_mean_ploidy_mean_across_fixed_p,
          linetype = reference_label
        ),
        inherit.aes = FALSE, colour = "#D62728",
        linewidth = 0.82, alpha = 1, lineend = "round"
      ) +
      ggplot2::scale_fill_viridis_c(
        option = "D", begin = 0.05, end = 0.95,
        limits = color_limits, breaks = color_breaks,
        trans = "log10", oob = scales::squish,
        na.value = "#EFEFEF",
        labels = function(x) formatC(x, format = "f", digits = 3),
        name = paste0(
          "Median required\np_miss,eff\n(log colors)"
        )
      ) +
      ggplot2::scale_colour_manual(
        values = c("Multiple inverse solutions" = "#7B3294"),
        labels = c("Multiple inverse\nsolutions (hatched)"),
        name = "Inverse status"
      ) +
      ggplot2::scale_linetype_manual(
        values = c(
          "0.01" = "solid", "0.10" = "F28282",
          "0.20" = "dotdash", "0.30" = "dotted",
          "Mean across 496 fixed p_miss,eff values" = "solid"
        ),
        breaks = reference_levels,
        labels = c(
          "p_miss,eff = 0.01", "p_miss,eff = 0.10",
          "p_miss,eff = 0.20", "p_miss,eff = 0.30",
          "Mean across 496 fixed\np_miss,eff values"
        ),
        name = "Reference curves"
      ) +
      ggplot2::scale_shape_manual(
        values = c("No stable unique inverse" = 22),
        labels = c("No stable unique\ninverse (gray)"),
        name = "Availability"
      ) +
      ggplot2::scale_x_continuous(
        breaks = 0:5, expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        breaks = 1:7, expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 5), ylim = c(1, 7), clip = "on", expand = FALSE
      ) +
      ggplot2::labs(
        title = display_manifest$display_label[[index]],
        x = "Fixed oxygen (%)", y = "Target dominant mean ploidy"
      ) +
      f6r_theme(8) +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "right",
        legend.box = "vertical",
        legend.title = ggplot2::element_text(size = 7.2, lineheight = 0.9),
        legend.text = ggplot2::element_text(size = 6.6, lineheight = 0.9),
        legend.key.height = grid::unit(2.8, "mm"),
        legend.spacing.y = grid::unit(0.7, "mm"),
        legend.box.spacing = grid::unit(1.0, "mm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title.position = "panel",
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_colourbar(
          barheight = grid::unit(22, "mm"),
          barwidth = grid::unit(3.2, "mm"), order = 1
        ),
        colour = ggplot2::guide_legend(
          keywidth = grid::unit(8, "mm"),
          keyheight = grid::unit(3.2, "mm"), order = 2,
          override.aes = list(linetype = "dotted", linewidth = 0.35)
        ),
        linetype = ggplot2::guide_legend(
          keywidth = grid::unit(8, "mm"),
          keyheight = grid::unit(3.2, "mm"), order = 3,
          override.aes = list(
            colour = c(rep("#111111", 4L), "#D62728"),
            linewidth = c(rep(0.62, 4L), 0.82),
            alpha = 1
          )
        ),
        shape = ggplot2::guide_legend(
          order = 4,
          override.aes = list(
            alpha = 1, fill = "#EFEFEF", colour = "#888888",
            linetype = 0
          )
        )
      )
    if (index != 1L) {
      p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
    p
  })
  composite <- f6r_compose_three_panel_row(
    plots,
    title = "B. Effective missegregation required for target ploidy",
    legend_rel_width = 0.72
  )
  output <- f6r_save_plot(
    composite,
    file.path(
      paths$figure6, "panels",
      "pair_inverse_o2_target_ploidy_required_p_miss_eff_primary_family_grid"
    ),
    width = f6r_primary_row_width(), height = 3.55
  )
  invisible(list(
    plot = composite, paths = output,
    data = list(inverse = inverse_bundle, overlay = overlay_bundle)
  ))
}

f6r_draw_supp6_2 <- function(paths) {
  f6r_require_packages(c("ggplot2", "patchwork", "scales", "magick"))
  ksel <- f6r_read_tsv(file.path(paths$figure6, "cluster_k_selection.tsv"))
  bootstrap_frequency <- f6r_read_tsv(file.path(
    paths$figure6, "cluster_k_resample_selection_frequency.tsv"
  ))
  stability <- f6r_read_tsv(file.path(
    paths$figure6, "cluster_stability.tsv"
  ))
  acceptance <- f6r_read_tsv(file.path(
    paths$figure6, "joint_seed_acceptance.tsv"
  ))
  robustness <- f6r_read_tsv(file.path(
    paths$figure6, "joint_seed_claim_robustness.tsv"
  ))
  cutoff_consistency <- f6r_read_tsv(file.path(
    paths$figure6, "joint_seed_cutoff_consistency.tsv"
  ))
  endpoint_counts <- f6r_read_tsv(file.path(
    paths$figure6, "joint_unique_parameter_endpoint_counts.tsv"
  ))
  dedup_comparison <- f6r_read_tsv(file.path(
    paths$figure6, "joint_seed_vs_unique_parameter_robustness.tsv"
  ))
  primary <- ksel[ksel$analysis_level == "primary pooled t-SNE regions", ]
  selected_k_frequency <- paste0(
    "k=", bootstrap_frequency$selected_k, ": ",
    bootstrap_frequency$n_subsample, "/", sum(bootstrap_frequency$n_subsample),
    collapse = "; "
  )

  p_a <- ggplot2::ggplot(primary, ggplot2::aes(k, average_silhouette)) +
    ggplot2::geom_line(colour = "#555555", linewidth = 0.55) +
    ggplot2::geom_point(
      ggplot2::aes(fill = selected_for_warm_starts),
      shape = 21, size = 2.4, colour = "#222222", stroke = 0.45
    ) +
    ggplot2::annotate(
      "text", x = f6r_family_count(),
      y = primary$average_silhouette[primary$k == f6r_family_count()] - 0.015,
      label = paste0("saved k=", f6r_family_count()),
      size = 2.45, hjust = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "white", `TRUE` = "#0072B2"), guide = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = 2:8) +
    ggplot2::labs(
      tag = "A", title = "Primary warm-start-region selection",
      subtitle = paste0(
        "The joint workflow uses the saved k=", f6r_family_count(),
        " primary partition.\n",
        "80% subsample silhouette maxima: ", selected_k_frequency, "."
      ),
      x = "Number of regions (k)", y = "Average silhouette"
    ) + f6r_theme()

  stability$perturbation_label <- c(
    "One-start\ninitialization",
    "80% endpoint\nsubsample",
    "14-parameter\nspace"
  )
  p_b <- ggplot2::ggplot(
    stability, ggplot2::aes(perturbation_label, ari_median)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ari_q05, ymax = ari_q95),
      width = 0.14, linewidth = 0.5, colour = "#555555"
    ) +
    ggplot2::geom_point(
      shape = 21, size = 2.7, fill = "#0072B2", colour = "#222222",
      stroke = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = 0.9, linetype = "dashed", linewidth = 0.35,
      colour = "#777777"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      tag = "B", title = "Saved primary-partition stability",
      subtitle = "Median adjusted Rand index; bars show 5th-95th percentiles when replicated.",
      x = NULL, y = "Adjusted Rand index"
    ) + f6r_theme()

  objective_floor <- 1e-4
  acceptance$delta_display <- pmax(acceptance$delta_objective, objective_floor)
  q10_endpoint_counts <- endpoint_counts[
    endpoint_counts$cutoff == "q10",
    c("pair_label", "n_seed", "n_unique_parameter_endpoint")
  ]
  q10_endpoint_counts$facet_label <- with(
    q10_endpoint_counts,
    paste0(
      pair_label,
      "\n", n_seed, " seeds\n", n_unique_parameter_endpoint,
      " unique endpoints"
    )
  )
  acceptance$facet_label <- q10_endpoint_counts$facet_label[
    match(acceptance$pair_label, q10_endpoint_counts$pair_label)
  ]
  p_c <- ggplot2::ggplot(
    acceptance, ggplot2::aes(objective_rank, delta_display)
  ) +
    ggplot2::annotate(
      "rect", xmin = 0.5, xmax = 50.5,
      ymin = objective_floor, ymax = Inf,
      fill = "#0072B2", alpha = 0.10
    ) +
    ggplot2::geom_line(colour = "#666666", linewidth = 0.45) +
    ggplot2::geom_vline(
      data = data.frame(
        xintercept = c(25.5, 50.5, 100.5),
        cutoff = factor(c("5%", "10%", "20%"), levels = c("5%", "10%", "20%"))
      ),
      ggplot2::aes(xintercept = xintercept, linetype = cutoff),
      colour = "#333333", linewidth = 0.35, show.legend = FALSE
    ) +
    ggplot2::scale_linetype_manual(values = c(`5%` = "dotted", `10%` = "solid", `20%` = "dashed")) +
    ggplot2::facet_wrap(~facet_label, nrow = 2L, scales = "free_y") +
    ggplot2::scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
    ggplot2::scale_x_continuous(breaks = c(1, 100, 250, 500)) +
    ggplot2::labs(
      tag = "C", title = "Pair-specific full-MAP objective eligibility",
      subtitle = paste0(
        "Blue region: primary lowest 10% (50/500); lines mark 5%, 10%, ",
        "and 20%.\nExact zero is displayed at 1e-4."
      ),
      x = "Objective rank within warm-start pair",
      y = expression(Delta*" full-MAP objective (log10)")
    ) + f6r_theme() +
    ggplot2::theme(
      panel.spacing = grid::unit(4, "mm"),
      strip.text = ggplot2::element_text(size = 6.5, lineheight = 0.9)
    )

  robust <- robustness[robustness$cutoff == "q10", , drop = FALSE]
  robust <- merge(
    robust,
    cutoff_consistency[, c(
      "pair_id", "claim_id", "same_modal_result_all_cutoffs",
      "minimum_modal_support_across_cutoffs"
    )],
    by = c("pair_id", "claim_id"), all.x = TRUE, sort = FALSE
  )
  dedup_q10 <- dedup_comparison[
    dedup_comparison$cutoff == "q10",
    c("pair_id", "claim_id", "same_modal_result")
  ]
  robust <- merge(
    robust, dedup_q10,
    by = c("pair_id", "claim_id"), all.x = TRUE, sort = FALSE
  )
  claim_order <- unique(robustness$claim_label)
  robust$claim_label <- factor(
    robust$claim_label, levels = rev(claim_order)
  )
  robust$pair_label <- factor(
    robust$pair_label,
    levels = f6r_family_levels()
  )
  robust$display_result <- ifelse(
    robust$modal_result == "TRUE", "yes",
    ifelse(robust$modal_result == "FALSE", "no", robust$modal_result)
  )
  robust$display <- paste0(
    robust$display_result, "\n",
    sprintf("%.0f%%", 100 * robust$modal_support_fraction), " ",
    ifelse(robust$same_modal_result_all_cutoffs, "all q", "varies")
  )
  p_d <- ggplot2::ggplot(
    robust,
    ggplot2::aes(
      pair_label, claim_label,
      fill = minimum_modal_support_across_cutoffs
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = display), size = 2.0, lineheight = 0.85) +
    ggplot2::geom_point(
      data = robust[!robust$same_modal_result, , drop = FALSE],
      shape = 4, size = 2.0, stroke = 0.55, colour = "#B2182B",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F2F2F2", high = "#0072B2", limits = c(0.5, 1),
      name = "Minimum modal support\nacross 5%, 10%, 20% sets"
    ) +
    ggplot2::scale_y_discrete(labels = scales::label_wrap(38)) +
    ggplot2::labs(
      tag = "D", title = "Seed-, cutoff-, and endpoint-multiplicity robustness",
      subtitle = paste0(
        "Primary 10% modal result/support; 'all q' = unchanged at 5%, 10%, and 20%.\n",
        "Red cross = changed after exact-parameter endpoint deduplication."
      ),
      x = "Warm-start pair", y = NULL
    ) + f6r_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6.4)
    )

  composite <- ((p_a | p_b) / (p_c | p_d) +
    patchwork::plot_layout(heights = c(0.78, 1.35), widths = c(1, 1.35)) +
    patchwork::plot_annotation(
      caption = paste0(
        "Primary regions are optimizer-derived search strata. Numerical seeds are ",
        "search endpoints, not biological replicates or confidence intervals."
      )
    )) & ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 7, colour = "#555555", hjust = 0)
    )
  output <- f6r_save_plot(
    composite,
    file.path(paths$supp6_2, "supp_fig6-2_joint_ensemble_robustness"),
    width = 11.0, height = 8.4
  )
  invisible(list(plot = composite, paths = output))
}

f6r_publish <- function(source, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Failed to publish ", source, " -> ", destination)
  }
  normalizePath(destination, mustWork = TRUE)
}

f6r_draw_main <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "patchwork", "magick", "scales"))
  paths <- f6r_paths(workspace_root)
  f6r_chart_contract(paths)
  panel_dir <- file.path(paths$figure6, "panels")
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

  # Main Figure 6 contains the joint-fit response surface and the endpoint-wise
  # inverse-response map with selected fixed-input and mean-ploidy overlays.
  # The former response-class panels are published independently as
  # Supplementary Figure 6-1.
  panel_a <- f6r_draw_surface_panel(paths)
  panel_b <- f6r_draw_inverse_response_panel(paths)
  f6r_require_files(
    c(panel_a$paths[["png"]], panel_b$paths[["png"]]),
    "Figure 6 panel"
  )
  dir.create(file.path(paths$figure6, "rendered"), recursive = TRUE, showWarnings = FALSE)
  output_png <- file.path(paths$figure6, "rendered", "assembled_oxygen_response.png")
  output_pdf <- file.path(paths$figure6, "rendered", "assembled_oxygen_response.pdf")
  image_a <- magick::image_read(panel_a$paths[["png"]])
  image_b <- magick::image_read(panel_b$paths[["png"]])
  main_width <- max(
    magick::image_info(image_a)$width, magick::image_info(image_b)$width
  )
  image_a <- magick::image_resize(image_a, paste0(main_width, "x"))
  image_b <- magick::image_resize(image_b, paste0(main_width, "x"))
  assembled <- magick::image_append(c(image_a, image_b), stack = TRUE)
  magick::image_write(assembled, output_png, format = "png", density = 300)
  magick::image_write(assembled, output_pdf, format = "pdf", density = 300)
  published <- c(
    figures_png = f6r_publish(
      output_png, file.path(paths$figures, "assembled_fig6.png")
    ),
    figures_pdf = f6r_publish(
      output_pdf, file.path(paths$figures, "assembled_fig6.pdf")
    ),
    manuscript_png = f6r_publish(
      output_png, file.path(paths$manuscript_figures, "assembled_fig6.png")
    ),
    manuscript_pdf = f6r_publish(
      output_pdf, file.path(paths$manuscript_figures, "assembled_fig6.pdf")
    )
  )
  info <- magick::image_info(magick::image_read(output_png))
  expected_main_width <- main_width
  expected_main_height <- magick::image_info(image_a)$height[[1L]] +
    magick::image_info(image_b)$height[[1L]]
  inverse_bundle <- panel_b$data$inverse
  overlay_bundle <- panel_b$data$overlay
  curve_data <- overlay_bundle$curve_data
  highlighted <- overlay_bundle$highlighted
  mean_curve <- overlay_bundle$mean_curve
  curve_o2_counts <- table(interaction(
    curve_data$pair_id, curve_data$effective_p_misseg, drop = TRUE
  ))
  curve_p_counts <- table(factor(
    curve_data$pair_label,
    levels = f6r_family_levels()
  )) / 201L
  validation <- data.frame(
    check = c(
      "main_figure_exists", "main_figure_width", "main_figure_height",
      "primary_family_seed_count", "figure6a_displayed_pair_labels",
      "figure6b_displayed_pair_labels", "figure6b_inverse_grid_rows",
      "figure6b_inverse_validation", "figure6b_fixed_p_reference_values",
      "figure6b_fixed_p_reference_oxygen_count",
      "figure6b_mean_curve_rows", "figure6b_mean_curve_values_finite",
      "figure6b_curve_family_rows", "figure6b_fixed_p_count_per_pair",
      "figure6b_oxygen_count_per_fixed_p",
      "figure6b_endpoint_count_per_grid_point",
      "figure6b_effective_missegregation_positive",
      "figure6b_dense_model_validation",
      "surface_grid_per_pair",
      "trajectory_grid_per_pair", "operator_qc_all_primary",
      "primary_cutoff_pair_specific"
    ),
    observed = c(
      as.character(file.exists(output_png)), info$width[[1L]], info$height[[1L]],
      paste(unique(f6r_read_tsv(file.path(
        paths$figure6, "joint_multiseed_surface_summary.tsv"
      ))$n_seed), collapse = ","),
      paste(f6r_read_tsv(file.path(
        paths$figure6, "figure6a_displayed_pairs.tsv"
      ))$pair_label, collapse = ","),
      paste(inverse_bundle$class_summary$pair_label, collapse = ","),
      nrow(inverse_bundle$summary),
      all(f6r_read_tsv(file.path(
        paths$figure6, "figure6_inverse_validation.tsv"
      ))$passed),
      paste(sort(unique(highlighted$effective_p_misseg)), collapse = ","),
      paste(sort(unique(as.integer(table(interaction(
        highlighted$pair_id, highlighted$effective_p_misseg, drop = TRUE
      ))))), collapse = ","),
      nrow(mean_curve),
      all(is.finite(mean_curve$dominant_mean_ploidy_mean_across_fixed_p)),
      nrow(curve_data),
      paste(as.integer(curve_p_counts), collapse = ","),
      paste(sort(unique(as.integer(curve_o2_counts))), collapse = ","),
      paste(sort(unique(curve_data$n_seed)), collapse = ","),
      all(curve_data$effective_p_misseg > 0),
      all(f6r_read_tsv(file.path(
        paths$figure6, "figure6d_dense_model_validation.tsv"
      ))$passed),
      201L * 60L, 201L,
      all(f6r_read_tsv(file.path(
        paths$figure6, "joint_multiseed_operator_qc.tsv"
      ))$operator_qc_pass),
      all(f6r_read_tsv(file.path(
        paths$figure6, "joint_seed_acceptance_summary.tsv"
      ))$n_accepted == rep(c(25L, 50L, 100L), f6r_family_count()))
    ),
    expected = c(
      "TRUE", as.character(expected_main_width),
      as.character(expected_main_height), "25,50",
      paste(f6r_family_levels(), collapse = ","),
      paste(f6r_family_levels(), collapse = ","),
      as.character(f6r_family_count() * 201L * 241L),
      "TRUE", "0.01,0.1,0.2,0.3",
      "201", as.character(f6r_family_count() * 201L), "TRUE",
      as.character(f6r_family_count() * 201L * 496L),
      paste(rep(496L, f6r_family_count()), collapse = ","), "201",
      "50", "TRUE", "TRUE",
      "12060", "201", "TRUE", "TRUE"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == validation$expected
  # The multi-cutoff table contains 25/50/100; the main panel itself is q10=50.
  validation$passed[validation$check == "primary_family_seed_count"] <-
    grepl("50", validation$observed[validation$check == "primary_family_seed_count"])
  f6r_write_tsv(
    validation, file.path(paths$figure6, "figure6_multiseed_validation.tsv")
  )
  output_files <- c(
    panel_a_png = panel_a$paths[["png"]],
    panel_a_pdf = panel_a$paths[["pdf"]],
    panel_b_png = panel_b$paths[["png"]],
    panel_b_pdf = panel_b$paths[["pdf"]],
    panel_b_inverse_response_summary_tsv =
      inverse_bundle$paths[["response_summary"]],
    panel_b_inverse_class_summary_tsv =
      inverse_bundle$paths[["class_summary"]],
    panel_b_inverse_ploidy4_anchor_summary_tsv =
      inverse_bundle$paths[["ploidy4_anchor_summary"]],
    panel_b_manuscript_class_summary_tsv =
      inverse_bundle$paths[["manuscript_class_summary"]],
    panel_b_manuscript_ploidy4_anchor_summary_tsv =
      inverse_bundle$paths[["manuscript_ploidy4_anchor_summary"]],
    panel_b_inverse_endpoint_manifest_tsv =
      inverse_bundle$paths[["endpoint_manifest"]],
    panel_b_inverse_endpoint_qc_tsv =
      inverse_bundle$paths[["endpoint_qc"]],
    panel_b_inverse_validation_tsv =
      inverse_bundle$paths[["validation"]],
    panel_b_displayed_pairs_tsv = overlay_bundle$paths[["displayed_pairs"]],
    panel_b_curve_family_tsv = overlay_bundle$paths[["curve_family"]],
    panel_b_highlighted_trajectories_tsv =
      overlay_bundle$paths[["highlighted_trajectories"]],
    panel_b_mean_ploidy_across_fixed_p_tsv =
      overlay_bundle$paths[["mean_ploidy_across_fixed_p"]],
    panel_b_dense_endpoint_manifest_tsv = file.path(
      paths$figure6, "figure6d_dense_endpoint_manifest.tsv"
    ),
    panel_b_dense_endpoint_qc_tsv = file.path(
      paths$figure6, "figure6d_dense_endpoint_qc.tsv"
    ),
    panel_b_dense_model_validation_tsv =
      overlay_bundle$paths[["dense_model_validation"]],
    panel_b_overlay_validation_tsv = overlay_bundle$paths[["validation"]],
    assembled_png = output_png,
    assembled_pdf = output_pdf,
    published
  )
  output_manifest <- data.frame(
    artifact = names(output_files),
    path = unname(output_files),
    exists = file.exists(output_files),
    size_bytes = file.info(output_files)$size,
    md5 = vapply(output_files, f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(
    output_manifest, file.path(paths$figure6, "figure6_output_manifest.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Figure 6 validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(list(output = c(png = output_png, pdf = output_pdf), published = published))
}

f6r_draw_supplement_6_1 <- function(
    workspace_root = f6r_find_workspace_root()
) {
  f6r_require_packages(c(
    "ggplot2", "cowplot", "magick", "shadowtext", "scales"
  ))
  paths <- f6r_paths(workspace_root)
  f6r_chart_contract(paths)
  f6r_load_response_engine(paths)
  panel_dir <- file.path(paths$figure6, "panels")
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

  smooth <- f6r_read_tsv(file.path(
    paths$figure6, "response_class_smoothed_curves.tsv"
  ))
  by_seed <- f6r_read_tsv(file.path(
    paths$figure6, "response_class_curve_class_by_seed.tsv"
  ))
  panel_a_paths <- response_draw_panel_a(
    smooth, by_seed,
    file.path(panel_dir, "response_class_fixed_o2_response_classes_4x2.png"),
    file.path(panel_dir, "response_class_fixed_o2_response_classes_4x2.pdf")
  )
  panel_bc_bundle <- response_build_panel_c(
    file.path(
      paths$figure5,
      "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
    ),
    by_seed, paths$figure6, panel_dir, panel_tag = "B"
  )
  f6r_require_files(
    c(panel_a_paths, panel_bc_bundle$figure_paths),
    "Supplementary Figure 6-1 panel"
  )

  rendered_dir <- file.path(paths$figure6, "rendered")
  dir.create(rendered_dir, recursive = TRUE, showWarnings = FALSE)
  output_png <- file.path(
    rendered_dir, "supp_fig6-1_response_class_diagnostics.png"
  )
  output_pdf <- file.path(
    rendered_dir, "supp_fig6-1_response_class_diagnostics.pdf"
  )
  image_a <- magick::image_read(panel_a_paths[["png"]])
  image_bc <- magick::image_read(panel_bc_bundle$figure_paths[["png"]])
  output_width <- max(
    magick::image_info(image_a)$width, magick::image_info(image_bc)$width
  )
  image_a <- magick::image_resize(image_a, paste0(output_width, "x"))
  image_bc <- magick::image_resize(image_bc, paste0(output_width, "x"))
  assembled <- magick::image_append(c(image_a, image_bc), stack = TRUE)
  magick::image_write(assembled, output_png, format = "png", density = 300)
  magick::image_write(assembled, output_pdf, format = "pdf", density = 300)

  published <- c(
    figures_png = f6r_publish(
      output_png,
      file.path(paths$figures, "supp_fig6-1_response_class_diagnostics.png")
    ),
    figures_pdf = f6r_publish(
      output_pdf,
      file.path(paths$figures, "supp_fig6-1_response_class_diagnostics.pdf")
    ),
    manuscript_png = f6r_publish(
      output_png,
      file.path(
        paths$manuscript_figures,
        "supp_fig6-1_response_class_diagnostics.png"
      )
    ),
    manuscript_pdf = f6r_publish(
      output_pdf,
      file.path(
        paths$manuscript_figures,
        "supp_fig6-1_response_class_diagnostics.pdf"
      )
    )
  )
  info <- magick::image_info(magick::image_read(output_png))
  class_counts <- table(factor(
    by_seed$smooth_curve_class,
    levels = response_curve_class_order
  ))
  declared_class_counts <- f6r_read_tsv(file.path(
    paths$figure6, "response_class_class_counts.tsv"
  ))
  declared_class_counts <- declared_class_counts[
    match(response_curve_class_order, declared_class_counts$smooth_curve_class),
    , drop = FALSE
  ]
  if (anyNA(declared_class_counts$smooth_curve_class) ||
      sum(declared_class_counts$n_seed) != 500L) {
    stop("The fresh in-vivo response-class count table is incomplete.")
  }
  validation <- data.frame(
    check = c(
      "figure_exists", "figure_width", "figure_height",
      "panel_a_grid_rows", "panel_a_grid_columns", "response_class_count",
      "response_class_seed_counts", "response_class_total_seed_count",
      "parameter_space_endpoint_count", "class_best_ring_count",
      "publication_png_md5_match", "publication_pdf_md5_match"
    ),
    observed = c(
      file.exists(output_png), info$width[[1L]], info$height[[1L]],
      4L, 2L, length(class_counts),
      paste(as.integer(class_counts), collapse = ","), sum(class_counts),
      nrow(panel_bc_bundle$data),
      sum(panel_bc_bundle$data$class_best_objective),
      f6r_md5(output_png) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output_pdf) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(
      "TRUE", as.character(output_width),
      as.character(
        magick::image_info(image_a)$height[[1L]] +
          magick::image_info(image_bc)$height[[1L]]
      ),
      "4", "2", "8",
      paste(declared_class_counts$n_seed, collapse = ","),
      "500", "500", "8",
      "TRUE", "TRUE"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == validation$expected
  validation_path <- f6r_write_tsv(
    validation, file.path(paths$figure6, "supp_fig6-1_validation.tsv")
  )
  output_files <- c(
    panel_a_png = panel_a_paths[["png"]],
    panel_a_pdf = panel_a_paths[["pdf"]],
    panel_bc_png = panel_bc_bundle$figure_paths[["png"]],
    panel_bc_pdf = panel_bc_bundle$figure_paths[["pdf"]],
    selection_data = panel_bc_bundle$data_paths[["selection"]],
    selection_validation = panel_bc_bundle$data_paths[["validation"]],
    assembled_png = output_png, assembled_pdf = output_pdf,
    validation = validation_path, published
  )
  manifest <- data.frame(
    artifact = names(output_files), path = unname(output_files),
    exists = file.exists(output_files), size_bytes = file.info(output_files)$size,
    md5 = vapply(output_files, f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(
    manifest, file.path(paths$figure6, "supp_fig6-1_output_manifest.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Supplementary Figure 6-1 validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(list(
    output = c(png = output_png, pdf = output_pdf),
    published = published, validation = validation
  ))
}

f6r_draw_supplement_6_2 <- function(workspace_root = f6r_find_workspace_root()) {
  paths <- f6r_paths(workspace_root)
  family_count <- f6r_family_count()
  drawn <- f6r_draw_supp6_2(paths)
  published <- c(
    figures_png = f6r_publish(
      drawn$paths[["png"]],
      file.path(paths$figures, "supp_fig6-2_joint_ensemble_robustness.png")
    ),
    figures_pdf = f6r_publish(
      drawn$paths[["pdf"]],
      file.path(paths$figures, "supp_fig6-2_joint_ensemble_robustness.pdf")
    ),
    manuscript_png = f6r_publish(
      drawn$paths[["png"]],
      file.path(paths$manuscript_figures, "supp_fig6-2_joint_ensemble_robustness.png")
    ),
    manuscript_pdf = f6r_publish(
      drawn$paths[["pdf"]],
      file.path(paths$manuscript_figures, "supp_fig6-2_joint_ensemble_robustness.pdf")
    )
  )
  info <- magick::image_info(magick::image_read(drawn$paths[["png"]]))
  validation <- data.frame(
    check = c(
      "supplement_exists", "supplement_width", "supplement_height",
      "cluster_k_rows", "cluster_stability_rows", "objective_seed_rows",
      "robustness_pair_claim_cutoff_rows", "cutoff_consistency_rows",
      "unique_endpoint_count_rows", "seed_dedup_claim_comparison_rows",
      "seed_dedup_surface_comparison_rows",
      "seed_dedup_trajectory_comparison_rows",
      "response_delta_summary_rows", "robustness_audit_summary_rows"
    ),
    observed = c(
      as.character(file.exists(drawn$paths[["png"]])), info$width[[1L]], info$height[[1L]],
      nrow(f6r_read_tsv(file.path(paths$figure6, "cluster_k_selection.tsv"))),
      nrow(f6r_read_tsv(file.path(paths$figure6, "cluster_stability.tsv"))),
      nrow(f6r_read_tsv(file.path(paths$figure6, "joint_seed_acceptance.tsv"))),
      nrow(f6r_read_tsv(file.path(paths$figure6, "joint_seed_claim_robustness.tsv"))),
      nrow(f6r_read_tsv(file.path(paths$figure6, "joint_seed_cutoff_consistency.tsv"))),
      nrow(f6r_read_tsv(file.path(
        paths$figure6, "joint_unique_parameter_endpoint_counts.tsv"
      ))),
      nrow(f6r_read_tsv(file.path(
        paths$figure6, "joint_seed_vs_unique_parameter_robustness.tsv"
      ))),
      nrow(f6r_read_tsv(file.path(
        paths$figure6, "joint_seed_vs_unique_parameter_surface_comparison.tsv"
      ))),
      nrow(f6r_read_tsv(file.path(
        paths$figure6, "joint_seed_vs_unique_parameter_trajectory_comparison.tsv"
      ))),
      nrow(f6r_read_tsv(file.path(
        paths$figure6, "joint_seed_response_delta_summary.tsv"
      ))),
      nrow(f6r_read_tsv(file.path(
        paths$figure6, "joint_seed_robustness_audit_summary.tsv"
      )))
    ),
    expected = c(
      "TRUE", "3300", "2520", "22", "9",
      as.character(500L * family_count),
      as.character(24L * family_count),
      as.character(8L * family_count),
      as.character(3L * family_count),
      as.character(24L * family_count),
      as.character(2L * family_count),
      as.character(3L * family_count),
      as.character(12L * family_count), "9"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == validation$expected
  f6r_write_tsv(
    validation, file.path(paths$supp6_2, "supp_figure6-2_validation.tsv")
  )
  output_files <- c(drawn$paths, published)
  output_manifest <- data.frame(
    artifact = names(output_files),
    path = unname(output_files),
    exists = file.exists(output_files),
    size_bytes = file.info(output_files)$size,
    md5 = vapply(output_files, f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(
    output_manifest, file.path(paths$supp6_2, "supp_figure6-2_output_manifest.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Supplementary Figure 6-2 validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(list(drawn = drawn, published = published))
}
