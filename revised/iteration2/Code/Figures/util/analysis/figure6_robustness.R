#!/usr/bin/env Rscript

# Figure 6 cluster-selection, objective-eligibility, multi-seed response-surface,
# robustness, and drawing utilities.
#
# Scientific-input boundary: every input resolved by this file is below the
# portable revised/iteration2 workspace.  No fit-result directory outside this
# workspace is read.  The embedded analytical response engine is the packaged
# Code/Figures/util/oxygen/response_pipeline.R implementation.

options(stringsAsFactors = FALSE, warn = 1)

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
      stop("Could not locate the iteration2 figure workspace from: ", start)
    }
    current <- parent
  }
}

f6r_paths <- function(workspace_root = f6r_find_workspace_root()) {
  root <- normalizePath(workspace_root, mustWork = TRUE)
  list(
    root = root,
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

f6r_endpoint_signatures <- function(master) {
  f6r_require_packages("data.table")
  shared <- f6r_shared_parameters()
  x <- data.table::as.data.table(master)[parameter %in% shared]
  signatures <- x[, {
    ordered_values <- vivo_natural[match(shared, parameter)]
    if (length(ordered_values) != length(shared) || any(!is.finite(ordered_values))) {
      stop("Cannot construct a finite 14-parameter endpoint signature.")
    }
    list(parameter_signature = paste(sprintf("%.17g", ordered_values), collapse = "|"))
  }, by = .(pair_id, seed_number)]
  if (nrow(signatures) != 3000L || anyDuplicated(signatures[, .(pair_id, seed_number)])) {
    stop("Expected one exact endpoint signature for each of 3,000 joint seeds.")
  }
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
    c("seed", "objective", "cluster_base_id", "tSNE1", "tSNE2"),
    drop = FALSE
  ]
  endpoints <- endpoints[order(endpoints$seed), , drop = FALSE]
  if (nrow(endpoints) != 500L || anyDuplicated(endpoints$seed)) {
    stop("Expected 500 unique in-vivo endpoints in the pooled embedding.")
  }
  expected_counts <- c(C01 = 99L, C02 = 385L, C03 = 16L)
  if (!identical(as.integer(table(endpoints$cluster_base_id)[names(expected_counts)]),
                 as.integer(expected_counts))) {
    stop("Primary warm-start-region counts do not match 99/385/16.")
  }

  coords <- as.matrix(endpoints[, c("tSNE1", "tSNE2")])
  coord_dist <- stats::dist(coords)
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
      selected_for_warm_starts = k == 3L,
      stringsAsFactors = FALSE
    )
  })
  primary <- do.call(rbind, primary_rows)
  primary$maximum_silhouette <- max(primary$average_silhouette)
  primary$difference_from_maximum <-
    primary$maximum_silhouette - primary$average_silhouette
  primary$within_0p001_of_maximum <-
    primary$difference_from_maximum <= 0.001 + 1e-12
  primary$selection_basis <- ifelse(
    primary$selected_for_warm_starts,
    paste0(
      "Retained as the smaller near-tied solution; k=3 is within 0.001 ",
      "of the numerical maximum, whereas k=4 creates a singleton region"
    ),
    "Candidate partition"
  )

  saved_silhouette <- f6r_silhouette(
    match(endpoints$cluster_base_id, c("C01", "C02", "C03")), coord_dist
  )

  parameter_wide <- f6r_endpoint_parameter_matrix(paths)
  joined <- merge(
    endpoints, parameter_wide,
    by.x = "seed", by.y = "seed_number", sort = TRUE
  )
  shared <- f6r_shared_parameters()
  expected_representatives <- list(
    C01 = c(Sc01 = 366L, Sc02 = 290L),
    C02 = c(Sc01 = 25L, Sc02 = 322L),
    C03 = c(Sc01 = 138L, Sc02 = 311L)
  )
  subcluster_rows <- list()
  assignment_rows <- list()
  reference_subcluster <- list()

  primary_order <- c("C01", "C02", "C03")
  for (i in seq_along(primary_order)) {
    region <- primary_order[[i]]
    z <- joined[joined$cluster_base_id == region, , drop = FALSE]
    x <- scale(as.matrix(z[, shared, drop = FALSE]))
    d <- stats::dist(x)
    selected_km <- NULL
    for (k in 2:min(6L, nrow(z) - 1L)) {
      set.seed(1123L + i - 1L)
      km <- stats::kmeans(x, centers = k, nstart = 10L, iter.max = 100L)
      sizes <- as.integer(table(km$cluster))
      subcluster_rows[[length(subcluster_rows) + 1L]] <- data.frame(
        analysis_level = "within-primary warm-start strata",
        region = region,
        feature_space = "14 transformed parameters z-scored within primary region",
        k = k,
        average_silhouette = f6r_silhouette(km$cluster, d),
        n_solution = nrow(z),
        minimum_cluster_size = min(sizes),
        maximum_cluster_size = max(sizes),
        selected_for_warm_starts = k == 2L,
        maximum_silhouette = NA_real_,
        difference_from_maximum = NA_real_,
        within_0p001_of_maximum = NA,
        selection_basis = if (k == 2L) {
          "Maximum average silhouette; used only to distribute two warm starts"
        } else {
          "Candidate partition"
        },
        stringsAsFactors = FALSE
      )
      if (k == 2L) selected_km <- km
    }
    if (is.null(selected_km)) stop("Failed to construct selected within-region strata.")
    minima <- vapply(
      sort(unique(selected_km$cluster)),
      function(g) z$seed[selected_km$cluster == g][
        which.min(z$objective[selected_km$cluster == g])
      ],
      integer(1L)
    )
    expected <- expected_representatives[[region]]
    if (!setequal(as.integer(minima), as.integer(expected))) {
      stop(
        "Reconstructed ", region, " representatives differ from the saved design: ",
        paste(minima, collapse = ",")
      )
    }
    label_map <- stats::setNames(
      names(expected)[match(minima, expected)], as.character(sort(unique(selected_km$cluster)))
    )
    assigned <- unname(label_map[as.character(selected_km$cluster)])
    assignment_rows[[region]] <- data.frame(
      seed = z$seed,
      primary_region = region,
      subcluster = assigned,
      warm_start_stratum = paste0(region, assigned),
      objective = z$objective,
      representative = z$seed %in% expected,
      stringsAsFactors = FALSE
    )
    reference_subcluster[[region]] <- list(
      x = x, labels = assigned, seeds = z$seed
    )
  }
  subclusters <- do.call(rbind, subcluster_rows)
  for (region in primary_order) {
    idx <- subclusters$region == region
    maximum <- max(subclusters$average_silhouette[idx])
    subclusters$maximum_silhouette[idx] <- maximum
    subclusters$difference_from_maximum[idx] <-
      maximum - subclusters$average_silhouette[idx]
    subclusters$within_0p001_of_maximum[idx] <-
      subclusters$difference_from_maximum[idx] <= 0.001 + 1e-12
  }
  k_selection <- rbind(primary, subclusters)
  assignments <- do.call(rbind, assignment_rows)
  rownames(assignments) <- NULL

  repeated_primary <- vapply(seq_len(n_resample), function(run) {
    set.seed(71000L + run)
    labels <- stats::kmeans(
      coords, centers = 3L, nstart = 1L, iter.max = 100L
    )$cluster
    f6r_adjusted_rand(endpoints$cluster_base_id, labels)
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
    km3 <- stats::kmeans(x, centers = 3L, nstart = 25L, iter.max = 200L)
    bootstrap_primary[run, "ari"] <- f6r_adjusted_rand(
      endpoints$cluster_base_id[idx], km3$cluster
    )
    sil <- vapply(2:8, function(k) {
      set.seed(74000L + 100L * run + k)
      km <- stats::kmeans(x, centers = k, nstart = 25L, iter.max = 200L)
      f6r_silhouette(km$cluster, d)
    }, numeric(1L))
    bootstrap_primary[run, "selected_k"] <- (2:8)[which.max(sil)]
  }

  global_parameter_x <- scale(as.matrix(joined[, shared, drop = FALSE]))
  set.seed(75001L)
  parameter_k3 <- stats::kmeans(
    global_parameter_x, centers = 3L, nstart = 100L, iter.max = 200L
  )$cluster
  parameter_space_ari <- f6r_adjusted_rand(
    joined$cluster_base_id, parameter_k3
  )

  stability_values <- list(
    data.frame(
      analysis_level = "primary pooled t-SNE regions",
      region = "all in-vivo endpoints",
      perturbation = paste0(n_resample, " one-start k-means initializations"),
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
        n_resample, " 80%-endpoint subsamples without replacement"
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
      perturbation = "k=3 clustering in standardized 14-parameter endpoint space",
      n_runs = 1L,
      ari_median = parameter_space_ari,
      ari_q05 = parameter_space_ari,
      ari_q95 = parameter_space_ari,
      proportion_ari_ge_0p9 = as.numeric(parameter_space_ari >= 0.9),
      stringsAsFactors = FALSE
    )
  )

  for (i in seq_along(primary_order)) {
    region <- primary_order[[i]]
    ref <- reference_subcluster[[region]]
    repeated <- vapply(seq_len(n_resample), function(run) {
      set.seed(76000L + 1000L * i + run)
      labels <- stats::kmeans(
        ref$x, centers = 2L, nstart = 1L, iter.max = 100L
      )$cluster
      f6r_adjusted_rand(ref$labels, labels)
    }, numeric(1L))
    bootstrap <- vapply(seq_len(n_resample), function(run) {
      set.seed(77000L + 1000L * i + run)
      idx <- sort(sample(
        seq_len(nrow(ref$x)), size = max(4L, floor(0.8 * nrow(ref$x)))
      ))
      set.seed(78000L + 1000L * i + run)
      labels <- stats::kmeans(
        ref$x[idx, , drop = FALSE], centers = 2L,
        nstart = 25L, iter.max = 200L
      )$cluster
      f6r_adjusted_rand(ref$labels[idx], labels)
    }, numeric(1L))
    stability_values[[length(stability_values) + 1L]] <- data.frame(
      analysis_level = "within-primary warm-start strata",
      region = region,
      perturbation = paste0(n_resample, " one-start k-means initializations"),
      n_runs = n_resample,
      ari_median = stats::median(repeated),
      ari_q05 = unname(stats::quantile(repeated, 0.05)),
      ari_q95 = unname(stats::quantile(repeated, 0.95)),
      proportion_ari_ge_0p9 = mean(repeated >= 0.9),
      stringsAsFactors = FALSE
    )
    stability_values[[length(stability_values) + 1L]] <- data.frame(
      analysis_level = "within-primary warm-start strata",
      region = region,
      perturbation = paste0(
        n_resample, " 80%-endpoint subsamples without replacement"
      ),
      n_runs = n_resample,
      ari_median = stats::median(bootstrap),
      ari_q05 = unname(stats::quantile(bootstrap, 0.05)),
      ari_q95 = unname(stats::quantile(bootstrap, 0.95)),
      proportion_ari_ge_0p9 = mean(bootstrap >= 0.9),
      stringsAsFactors = FALSE
    )
  }
  stability <- do.call(rbind, stability_values)
  bootstrap_frequency <- as.data.frame(table(
    selected_k = as.integer(bootstrap_primary[, "selected_k"])
  ))
  names(bootstrap_frequency)[[2L]] <- "n_subsample"
  bootstrap_frequency$proportion <-
    bootstrap_frequency$n_subsample / n_resample

  representative_summary <- assignments[assignments$representative, , drop = FALSE]
  representative_summary <- representative_summary[
    order(representative_summary$primary_region,
          representative_summary$subcluster), , drop = FALSE
  ]
  representative_summary$selected_k <- 2L
  representative_summary$selected_average_silhouette <- vapply(
    representative_summary$primary_region,
    function(region) subclusters$average_silhouette[
      subclusters$region == region & subclusters$k == 2L
    ],
    numeric(1L)
  )

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
          "saved_primary_k3_average_silhouette",
          "robust_k3_average_silhouette",
          "robust_k4_average_silhouette",
          "robust_k4_minus_k3_silhouette",
          "primary_parameter_space_k3_ARI",
          "tSNE_seed_perplexity_sensitivity"
        ),
        value = c(
          saved_silhouette,
          primary$average_silhouette[primary$k == 3L],
          primary$average_silhouette[primary$k == 4L],
          primary$average_silhouette[primary$k == 4L] -
            primary$average_silhouette[primary$k == 3L],
          parameter_space_ari,
          NA
        ),
        interpretation = c(
          "Silhouette of the saved C01/C02/C03 assignment",
          "High-restart audit of k=3",
          "High-restart numerical maximum; near-tied with k=3",
          "Absolute improvement is below 0.001",
          "Agreement between frozen-embedding regions and endpoint-space k=3",
          paste0(
            "Not estimable from the frozen iteration2 embedding: the 228,000 ",
            "source parameter vectors required to re-estimate t-SNE are not vendored"
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
  if (nrow(selected) != 6L) stop("Expected six selected joint pairs.")
  selected$pair_id <- paste0("fit_joint_", selected$warmup_label)
  selected$selected_seed_number <- as.integer(sub("^seed", "", selected$selected_seed))
  selected$pair_label <- sub(
    ".*_(C[0-9]+Sc[0-9]+)_vt.*", "\\1", selected$warmup_label
  )
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
  if (nrow(objectives) != 3000L ||
      any(table(objectives$pair_id) != 500L) ||
      anyDuplicated(objectives[, c("pair_id", "seed_number")])) {
    stop("Joint objective cache must contain six groups of 500 unique seeds.")
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
  if (nrow(parameter_qc) != 3000L ||
      any(parameter_qc$n_parameter_record != length(f6r_shared_parameters())) ||
      any(parameter_qc$n_parameter != length(f6r_shared_parameters()))) {
    stop("The local joint parameter cache is not complete for all 3,000 endpoints.")
  }
  endpoint_signatures <- f6r_endpoint_signatures(master)
  parameter_qc <- merge(
    as.data.frame(parameter_qc), endpoint_signatures,
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
  for (pair in unique(objectives$pair_id)) {
    idx <- which(objectives$pair_id == pair)
    signatures_in_rank_order <- unique(objectives$parameter_signature[idx])
    objectives$parameter_endpoint_group[idx] <- paste0(
      objectives$pair_label[idx], "_E",
      sprintf("%03d", match(objectives$parameter_signature[idx], signatures_in_rank_order))
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
    endpoint_multiplicity = endpoint_multiplicity,
    endpoint_counts = endpoint_counts,
    paths = output_paths,
    master_path = master_path
  ))
}

f6r_load_response_engine <- local({
  loaded <- FALSE
  function(paths) {
    if (loaded) return(invisible(TRUE))
    engine <- file.path(paths$code, "util", "oxygen", "response_pipeline.R")
    f6r_require_files(engine, "packaged analytical response engine")
    sys.source(engine, envir = globalenv())
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
    parameter_source, full_surface = TRUE, force_rebuild = FALSE
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
    if (!is.null(existing) &&
        identical(as.character(existing$metadata$pair_id), as.character(pair_id)) &&
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
  params <- stats::setNames(as.numeric(rows$vivo_natural), rows$parameter)
  params[["rho_2N"]] <- context$rho_2N
  run_params <- prepare_run_params(
    param_values = params,
    simulation = "joint",
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
  surface$pair_id <- pair_id
  surface$pair_label <- context$pair_label
  surface$seed_number <- seed_number
  surface$objective <- objective
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
        seed_number = seed_number, objective = objective,
        config_path = context$config_path,
        surface_profile = requested_surface_profile,
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
  if (nrow(selected) != 150L || any(observed_seed_counts != 50L)) {
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
  if (!identical(as.integer(unique_counts), c(50L, 15L, 50L)) ||
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
    force_rebuild = FALSE
) {
  p_values <- f6r_figure6d_p_values()
  o2_values <- f6r_figure6d_o2_values()
  expected_rows <- length(p_values) * length(o2_values)
  profile <- "dense_201x496_step0p001_v1"
  if (file.exists(cache_path) && !isTRUE(force_rebuild)) {
    existing <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(existing) &&
        identical(as.character(existing$metadata$profile), profile) &&
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
  required_parameters <- c(f6r_shared_parameters(), "rho_2N")
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
    simulation = "joint",
    cfg = context$config,
    fixed_o2 = 5
  )
  config <- prepare_sim_cfg(
    context$config, argv = list(), fixed_o2 = 5, run_params = run_params
  )
  run_params$O2_growth <- isTRUE(config$O2_growth)
  run_params$ploidy_O2_death <- config$ploidy_O2_death
  seed_id <- paste0("seed", seed_number)

  surface_rows <- lapply(p_values, function(p_fixed) {
    forced <- response_force_effective_p_misseg(run_params, p_fixed)
    rows <- lapply(o2_values, function(o2) {
      result <- fixo2_dominant_attractor_one(
        seed_id = seed_id,
        run_params = forced,
        model_env = globalenv(),
        cfg = config,
        O2 = o2
      )
      data.frame(
        O2_pct = as.numeric(result$O2_pct[[1L]]),
        effective_p_misseg = p_fixed,
        status = as.character(result$status[[1L]]),
        actual_effective_p_misseg =
          as.numeric(result$population_average_p_misseg[[1L]]),
        dominant_mean_ploidy =
          as.numeric(result$dominant_mean_ploidy[[1L]]),
        spectral_gap = as.numeric(result$spectral_gap[[1L]]),
        dominant_growth_rate =
          as.numeric(result$dominant_growth_rate[[1L]]),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  surface <- do.call(rbind, surface_rows)
  surface <- surface[order(
    surface$effective_p_misseg, surface$O2_pct
  ), , drop = FALSE]
  max_error <- max(abs(
    surface$actual_effective_p_misseg - surface$effective_p_misseg
  ), na.rm = TRUE)
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
      force_rebuild = rebuild
    )
    message(
      "Figure 6D dense endpoint complete: ", metadata$pair_label,
      " seed", metadata$representative_seed_number,
      " (", i, "/", nrow(endpoints), ")"
    )
    result
  }
  index <- seq_len(nrow(endpoints))
  if (.Platform$OS.type != "windows" && n_core > 1L) {
    qc_rows <- parallel::mclapply(
      index, compute_one,
      mc.cores = min(as.integer(n_core), length(index)),
      mc.preschedule = FALSE
    )
  } else {
    qc_rows <- lapply(index, compute_one)
  }
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
      "115", "50,50,50", "299088", "496,496,496", "201",
      "0.005--0.5 by 0.001", "TRUE", "<=1e-8", "1206", "<=1e-12"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- c(
    nrow(qc) == 115L,
    identical(
      as.integer(tapply(
        qc$endpoint_multiplicity_q10,
        factor(
          qc$pair_label,
          levels = endpoint_manifest$display_manifest$pair_label
        ),
        sum
      )),
      c(50L, 50L, 50L)
    ),
    nrow(summary) == 3L * 496L * 201L,
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
    nrow(overlap) == 3L * 2L * 201L,
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
      "Six warm-start pairs times eight prespecified qualitative diagnostics",
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
      existing_ok <- file.exists(pair_cache)
    }
    message(
      "Figure 6 multi-seed ", context$pair_label, ": ",
      sum(existing_ok), "/", length(pair_cache), " seed caches already available."
    )
    batches <- split(
      seq_len(nrow(eligible)),
      ceiling(seq_len(nrow(eligible)) / max(1L, as.integer(n_core)))
    )
    pair_qc <- list()
    for (batch_index in seq_along(batches)) {
      idx <- batches[[batch_index]]
      compute_one <- function(j) {
        cache_path <- pair_cache[[as.character(eligible$seed_number[[j]])]]
        f6r_compute_seed_cache(
          pair_id = pair,
          seed_number = eligible$seed_number[[j]],
          objective = eligible$objective[[j]],
          master = master,
          context = context,
          cache_path = cache_path,
          parameter_source = objective_bundle$master_path,
          full_surface = eligible$eligible_q10[[j]],
          force_rebuild = rebuild
        )
      }
      rows <- if (length(idx) > 1L && n_core > 1L &&
                  .Platform$OS.type != "windows") {
        parallel::mclapply(
          idx, compute_one,
          mc.cores = min(as.integer(n_core), length(idx)),
          mc.preschedule = FALSE
        )
      } else {
        lapply(idx, compute_one)
      }
      if (any(vapply(rows, function(x) inherits(x, "try-error"), logical(1L)))) {
        stop("One or more multi-seed workers failed for ", pair)
      }
      pair_qc <- c(pair_qc, rows)
      message(
        "  ", context$pair_label, " batch ", batch_index, "/",
        length(batches), " complete (", max(idx), "/", nrow(eligible), ")."
      )
    }
    cache_paths[[pair]] <- unname(pair_cache)
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
      "At each fixed effective missegregation probability, how does dominant ploidy respond to oxygen?",
      "Why were three primary warm-start regions retained?",
      "Why were two within-region warm-start strata retained?",
      "Which numerical endpoints satisfy the objective-based eligibility rule?",
      "Do qualitative oxygen-response diagnostics agree across objective-eligible endpoints and cutoffs?"
    ),
    chart_family = c(
      "four-by-two small-multiple response curves",
      "parameter-space scatter",
      "half-violin, boxplot, and seed-level points",
      "uncertainty heatmap plus trajectory interval",
      "small-multiple fixed-input curve family",
      "silhouette profile", "faceted silhouette profile",
      "faceted ranked objective curve", "robustness heatmap"
    ),
    data_grain = c(
      "eight response classes summarizing 500 separate in-vivo fitted endpoints across 201 oxygen values",
      "500 separate in-vivo fitted endpoints",
      "500 separate in-vivo fitted endpoints",
      "three displayed pairs x oxygen x effective missegregation, summarized over 50 seeds; six pairs retained analytically",
      "three displayed pairs x 496 fixed effective missegregation probabilities (0.005-0.500 by 0.001) x 201 oxygen values, summarized over 50 endpoints",
      "candidate k", "primary region x candidate k",
      "pair x numerical seed", "claim x pair x objective cutoff"
    ),
    primary_encoding = c(
      "individual fitted-endpoint curves by response-class color and pointwise class median in black",
      "response-class color; class-best black ring; warm-start-region outlines",
      "full-MAP objective density, quartiles, all seed-level endpoints, and global lowest-decile cutoff",
      "original-unit ploidy labels with log-scaled seed-weighted median fill; trajectory median and 10-90% band; low-consensus marks; unique-parameter endpoint sensitivity in source tables",
      "one oxygen-ploidy curve per fixed effective missegregation probability; shared axes and log-scaled curve colors",
      "average silhouette by k", "average silhouette by k",
      "delta full-MAP objective by rank",
      "minimum modal support across cutoffs; primary modal result and cutoff-consistency flag"
    ),
    caveat = c(
      "Response classes pool reliable, caution, and unreliable spectral-gap strata and are descriptive",
      "Pooled t-SNE axes are unitless and embedding-dependent; response-class separation is descriptive",
      "Objective distributions are descriptive; numerical endpoints are not biological replicates",
      "Post-fit asymptotic diagnostic; optimizer seeds are not biological replicates",
      "Alternative line representation of the standardized post-fit surface in Figure 6A; not independently measured CIN or biological uncertainty",
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
    file.path(paths$figure6, "selection_diagnostic_selection_data.tsv"),
    file.path(
      paths$figure5,
      "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
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

f6r_data <- function(
    workspace_root = f6r_find_workspace_root(), n_core = 8L,
    rebuild = FALSE, n_resample = 100L
) {
  paths <- f6r_paths(workspace_root)
  dir.create(paths$figure6, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$supp6_2, recursive = TRUE, showWarnings = FALSE)
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
  f6r_local_data_contract(paths)
  invisible(list(
    paths = paths,
    clusters = cluster_bundle,
    objective = objective_bundle,
    multiseed = multiseed_bundle,
    dense_figure6d = dense_figure6d_bundle
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
  if (length(plots) != 3L) {
    stop("A three-panel Figure 6 row requires exactly three plots.")
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
    align = "hv", axis = "tblr", rel_widths = rep(1, 3)
  )
  body <- cowplot::plot_grid(
    panel_row, legend, nrow = 1L,
    rel_widths = c(3, legend_rel_width)
  )
  cowplot::ggdraw() +
    cowplot::draw_plot(body, x = 0, y = 0, width = 1, height = 0.91) +
    cowplot::draw_label(
      title, x = 0.006, y = 0.995,
      hjust = 0, vjust = 1, fontface = "bold",
      fontfamily = "Helvetica", size = title_size
    )
}

f6r_display_pair_manifest <- function(pair_ids, panel_label) {
  display_pair_labels <- c(
    C01 = "C01Sc01",
    C02 = "C02Sc01",
    C03 = "C03Sc02"
  )
  pair_ids <- unique(as.character(pair_ids))
  pair_id_to_label <- stats::setNames(
    sub(".*_(C[0-9]+Sc[0-9]+)_vt.*", "\\1", pair_ids),
    pair_ids
  )
  ordered_pairs <- names(pair_id_to_label)[
    match(unname(display_pair_labels), pair_id_to_label)
  ]
  if (anyNA(ordered_pairs) || length(unique(ordered_pairs)) != 3L) {
    stop("Cannot resolve the three requested Figure 6", panel_label,
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
  if (length(unique(surface$pair_id)) != 6L ||
      any(table(surface$pair_id) != 201L * 60L) ||
      any(surface$n_seed != 50L) || any(trajectory$n_seed != 50L)) {
    stop("Primary Figure 6C summary must contain 50 seeds per pair.")
  }

  display_manifest <- f6r_display_pair_manifest(surface$pair_id, "C")
  display_pair_labels <- stats::setNames(
    display_manifest$pair_label, display_manifest$display_label
  )
  ordered_pairs <- display_manifest$pair_id
  f6r_write_tsv(
    display_manifest,
    file.path(paths$figure6, "figure6c_displayed_pairs.tsv")
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
    stop("Figure 6C log-scaled ploidy fill requires strictly positive values.")
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
      "pair_surface_o2_effective_p_misseg_three_pair_grid"
    ),
    width = 8.6, height = 3.55
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
  if (length(unique(curve_data$pair_id)) != 3L ||
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
  gap_boundary <- f6r_map_panel_c_gap_contour_to_panel_d(
    paths, curve_data, display_manifest
  )
  gap_boundary_path <- f6r_write_tsv(
    gap_boundary,
    file.path(paths$figure6, "figure6d_spectral_gap_boundary.tsv")
  )

  p_counts <- table(factor(
    curve_data$pair_label,
    levels = c("C01Sc01", "C02Sc01", "C03Sc02")
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
      paste(unique(gap_boundary$pair_label), collapse = ","),
      all(is.finite(gap_boundary$O2_pct)) &&
        all(is.finite(gap_boundary$effective_p_misseg)) &&
        all(is.finite(gap_boundary$dominant_mean_ploidy_median)),
      paste(unique(gap_boundary$criterion), collapse = ",")
    ),
    expected = c(
      "C01Sc01,C02Sc01,C03Sc02", "299088", "496,496,496", "201",
      "50", "TRUE", "TRUE", "TRUE", "0.005,0.5", "0.001",
      "<=1e-8", "TRUE", "0.01,0.1,0.2,0.3", "201",
      "C01Sc01,C02Sc01,C03Sc02", "TRUE",
      "proportion_spectral_gap_below_0p005"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- c(
    identical(display_manifest$pair_label, c("C01Sc01", "C02Sc01", "C03Sc02")),
    nrow(curve_data) == 299088L, all(p_counts == 496L),
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
    nrow(highlighted) == 4L * 3L * 201L &&
      all(table(interaction(
        highlighted$pair_id, highlighted$effective_p_misseg, drop = TRUE
      )) == 201L),
    identical(unique(gap_boundary$pair_label),
              c("C01Sc01", "C02Sc01", "C03Sc02")),
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
    gap_boundary = gap_boundary,
    manifest = display_manifest,
    paths = c(
      displayed_pairs = display_manifest_path,
      curve_family = curve_path,
      highlighted_trajectories = highlighted_path,
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
  display_manifest <- bundle$manifest
  color_limits <- range(curve_data$effective_p_misseg, finite = TRUE)
  if (any(!is.finite(color_limits)) || color_limits[[1L]] <= 0) {
    stop("Figure 6D log-scaled colors require positive finite values.")
  }
  color_breaks <- scales::breaks_log(n = 5)(color_limits)
  color_breaks <- color_breaks[
    color_breaks >= color_limits[[1L]] & color_breaks <= color_limits[[2L]]
  ]
  y_limits <- range(curve_data$dominant_mean_ploidy_median, finite = TRUE)
  y_padding <- max(0.04 * diff(y_limits), 0.03)
  y_limits <- y_limits + c(-y_padding, y_padding)

  plots <- lapply(seq_len(nrow(display_manifest)), function(i) {
    pair <- display_manifest$pair_id[[i]]
    d <- curve_data[curve_data$pair_id == pair, , drop = FALSE]
    d <- d[order(d$effective_p_misseg, d$O2_pct), , drop = FALSE]
    h <- highlighted[highlighted$pair_id == pair, , drop = FALSE]
    h <- h[order(h$effective_p_misseg, h$O2_pct), , drop = FALSE]
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
          "0.20" = "dotdash", "0.30" = "dotted"
        ),
        name = "Highlighted p_miss,eff"
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
          keyheight = grid::unit(3.2, "mm"), order = 2
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
      "pair_fixed_p_miss_eff_o2_ploidy_curve_family_three_pair_grid"
    ),
    width = 8.6, height = 3.55
  )
  invisible(list(plot = composite, paths = output, data = bundle))
}

f6r_draw_supp6_2 <- function(paths) {
  f6r_require_packages(c("ggplot2", "patchwork", "scales", "magick"))
  ksel <- f6r_read_tsv(file.path(paths$figure6, "cluster_k_selection.tsv"))
  bootstrap_frequency <- f6r_read_tsv(file.path(
    paths$figure6, "cluster_k_resample_selection_frequency.tsv"
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
  sub <- ksel[ksel$analysis_level == "within-primary warm-start strata", ]

  p_a <- ggplot2::ggplot(primary, ggplot2::aes(k, average_silhouette)) +
    ggplot2::geom_line(colour = "#555555", linewidth = 0.55) +
    ggplot2::geom_point(
      ggplot2::aes(fill = selected_for_warm_starts),
      shape = 21, size = 2.4, colour = "#222222", stroke = 0.45
    ) +
    ggplot2::annotate(
      "text", x = 3, y = primary$average_silhouette[primary$k == 3] - 0.015,
      label = "retained k=3", size = 2.45, hjust = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "white", `TRUE` = "#0072B2"), guide = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = 2:8) +
    ggplot2::labs(
      tag = "A", title = "Primary warm-start-region selection",
      subtitle = paste0(
        "k=3 and k=4 differ by <0.001; k=4 has a singleton.\n",
        "80% subsample re-selection: ",
        bootstrap_frequency$n_subsample[bootstrap_frequency$selected_k == 3],
        "/", sum(bootstrap_frequency$n_subsample), " versus ",
        bootstrap_frequency$n_subsample[bootstrap_frequency$selected_k == 4],
        "/", sum(bootstrap_frequency$n_subsample), "."
      ),
      x = "Number of regions (k)", y = "Average silhouette"
    ) + f6r_theme()

  p_b <- ggplot2::ggplot(
    sub, ggplot2::aes(k, average_silhouette, colour = region)
  ) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(
      ggplot2::aes(fill = selected_for_warm_starts),
      shape = 21, size = 2.0, stroke = 0.4
    ) +
    ggplot2::facet_wrap(~region, nrow = 1L) +
    ggplot2::scale_colour_manual(
      values = c(C01 = "#FFD700", C02 = "#9B59B6", C03 = "#39FF14"),
      guide = "none"
    ) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "white", `TRUE` = "#333333"), guide = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = 2:6) +
    ggplot2::labs(
      tag = "B", title = "Within-region warm-start strata",
      subtitle = paste0(
        "k=2 maximizes silhouette in each region.\n",
        "Low values limit interpretation to search coverage."
      ),
      x = "Number of strata (k)", y = "Average silhouette"
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
      sub("(C[0-9]+)(Sc[0-9]+)", "\\1 \\2", pair_label),
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
    levels = c("C01Sc01", "C01Sc02", "C02Sc01", "C02Sc02", "C03Sc01", "C03Sc02")
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
          "Primary regions and within-region partitions are optimizer-derived search strata. Numerical seeds are ",
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

  # Main Figure 6 now contains only the joint-fit response surface and its
  # fixed-missegregation curve-family representation. The former panels A-C
  # are published independently as Supplementary Figure 6-1.
  panel_a <- f6r_draw_surface_panel(paths)
  panel_b <- f6r_draw_fixed_p_curve_panel(paths)
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
  figure6d_curve_data <- f6r_read_tsv(file.path(
    paths$figure6, "figure6d_fixed_p_curve_family.tsv"
  ))
  figure6d_o2_counts <- table(interaction(
    figure6d_curve_data$pair_id, figure6d_curve_data$effective_p_misseg,
    drop = TRUE
  ))
  figure6d_p_counts <- table(factor(
    figure6d_curve_data$pair_label,
    levels = c("C01Sc01", "C02Sc01", "C03Sc02")
  )) / 201L
  validation <- data.frame(
    check = c(
      "main_figure_exists", "main_figure_width", "main_figure_height",
      "six_pair_primary_seed_count", "figure6a_displayed_pair_labels",
      "figure6b_displayed_pair_labels", "figure6b_curve_family_rows",
      "figure6b_fixed_p_count_per_pair", "figure6b_oxygen_count_per_fixed_p",
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
        paths$figure6, "figure6c_displayed_pairs.tsv"
      ))$pair_label, collapse = ","),
      paste(f6r_read_tsv(file.path(
        paths$figure6, "figure6d_displayed_pairs.tsv"
      ))$pair_label, collapse = ","),
      nrow(figure6d_curve_data),
      paste(as.integer(figure6d_p_counts), collapse = ","),
      paste(sort(unique(as.integer(figure6d_o2_counts))), collapse = ","),
      paste(sort(unique(figure6d_curve_data$n_seed)), collapse = ","),
      all(figure6d_curve_data$effective_p_misseg > 0),
      all(f6r_read_tsv(file.path(
        paths$figure6, "figure6d_dense_model_validation.tsv"
      ))$passed),
      201L * 60L, 201L,
      all(f6r_read_tsv(file.path(
        paths$figure6, "joint_multiseed_operator_qc.tsv"
      ))$operator_qc_pass),
      all(f6r_read_tsv(file.path(
        paths$figure6, "joint_seed_acceptance_summary.tsv"
      ))$n_accepted == rep(c(25L, 50L, 100L), 6L))
    ),
    expected = c(
      "TRUE", as.character(expected_main_width),
      as.character(expected_main_height), "25,50",
      "C01Sc01,C02Sc01,C03Sc02", "C01Sc01,C02Sc01,C03Sc02",
      "299088", "496,496,496", "201", "50", "TRUE", "TRUE",
      "12060", "201", "TRUE", "TRUE"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == validation$expected
  # The multi-cutoff table contains 25/50/100; the main panel itself is q10=50.
  validation$passed[validation$check == "six_pair_primary_seed_count"] <-
    grepl("50", validation$observed[validation$check == "six_pair_primary_seed_count"])
  f6r_write_tsv(
    validation, file.path(paths$figure6, "figure6_multiseed_validation.tsv")
  )
  output_files <- c(
    panel_a_png = panel_a$paths[["png"]],
    panel_a_pdf = panel_a$paths[["pdf"]],
    panel_b_png = panel_b$paths[["png"]],
    panel_b_pdf = panel_b$paths[["pdf"]],
    panel_b_displayed_pairs_tsv = panel_b$data$paths[["displayed_pairs"]],
    panel_b_curve_family_tsv = panel_b$data$paths[["curve_family"]],
    panel_b_highlighted_trajectories_tsv =
      panel_b$data$paths[["highlighted_trajectories"]],
    panel_b_spectral_gap_boundary_tsv =
      panel_b$data$paths[["spectral_gap_boundary"]],
    panel_b_dense_endpoint_manifest_tsv = file.path(
      paths$figure6, "figure6d_dense_endpoint_manifest.tsv"
    ),
    panel_b_dense_endpoint_qc_tsv = file.path(
      paths$figure6, "figure6d_dense_endpoint_qc.tsv"
    ),
    panel_b_dense_model_validation_tsv =
      panel_b$data$paths[["dense_model_validation"]],
    panel_b_validation_tsv = panel_b$data$paths[["validation"]],
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
      "4", "2", "8", "316,50,50,47,26,5,4,2", "500", "500", "8",
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
      "TRUE", "3300", "2520", "22", "9", "3000", "144", "48",
      "18", "144", "12", "18", "72", "9"
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
