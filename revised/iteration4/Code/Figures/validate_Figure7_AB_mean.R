#!/usr/bin/env Rscript

# Independent anchor checks for the Figure 7A/B arithmetic-mean estimands.
# The checks read the saved endpoint-level caches, retain q10 optimizer-seed
# multiplicity, and compare direct means with the published summary tables.

options(stringsAsFactors = FALSE, warn = 1)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
  } else {
    normalizePath(file.path(getwd(), "Code", "Figures"), mustWork = TRUE)
  }
})
source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))

f7ab_resolve_cache_path <- function(recorded_path, paths) {
  recorded_path <- as.character(recorded_path)
  if (file.exists(recorded_path)) {
    return(normalizePath(recorded_path, mustWork = TRUE))
  }
  marker <- "/data/Figures/"
  marker_position <- regexpr(marker, recorded_path, fixed = TRUE)[[1L]]
  if (marker_position < 1L) {
    stop("Cannot localize recorded Figure 7 cache path: ", recorded_path)
  }
  relative <- substring(recorded_path, marker_position + 1L)
  candidate <- file.path(paths$root, relative)
  if (!file.exists(candidate)) {
    stop("Localized Figure 7 cache does not exist: ", candidate)
  }
  normalizePath(candidate, mustWork = TRUE)
}

f7ab_exact_rows <- function(data, oxygen, p_misseg = NULL,
                            target_ploidy = NULL) {
  keep <- abs(data$O2_pct - oxygen) <= 1e-12
  if (!is.null(p_misseg)) {
    keep <- keep & abs(data$p_misseg - p_misseg) <= 1e-12
  }
  if (!is.null(target_ploidy)) {
    keep <- keep & abs(data$target_ploidy - target_ploidy) <= 1e-12
  }
  data[keep, , drop = FALSE]
}

f7ab_in_numeric_set <- function(values, reference, tolerance = 1e-12) {
  vapply(
    values,
    function(value) any(abs(value - reference) <= tolerance),
    logical(1L)
  )
}

f7ab_surface_anchor_rows <- function(paths, model_context) {
  suffix <- if (identical(model_context, "in vivo")) "" else "_invitro"
  surface_summary <- f7r_read_tsv(file.path(
    paths$figure7, paste0("joint_multiseed_surface_summary", suffix, ".tsv")
  ))
  trajectory_summary <- f7r_read_tsv(file.path(
    paths$figure7, paste0("joint_multiseed_trajectory_summary", suffix, ".tsv")
  ))
  endpoint_qc <- f7r_read_tsv(file.path(
    paths$figure7, paste0("joint_multiseed_operator_qc", suffix, ".tsv")
  ))
  if (identical(model_context, "in vivo")) {
    acceptance <- f7r_read_tsv(file.path(
      paths$figure7, "joint_seed_acceptance.tsv"
    ))
    eligible <- acceptance[
      acceptance$eligible_q10 & acceptance$hard_qc_pass,
      c("pair_id", "seed_number"), drop = FALSE
    ]
    eligible$endpoint_multiplicity_q10 <- 1L
    endpoint_qc <- merge(
      endpoint_qc, eligible, by = c("pair_id", "seed_number"), sort = FALSE
    )
  } else {
    eligible <- f7r_read_tsv(file.path(
      paths$figure7, "joint_invitro_q20_endpoint_manifest.tsv"
    ))
    eligible <- eligible[
      eligible$endpoint_multiplicity_q10 > 0L,
      c(
        "pair_id", "pair_label", "model_context",
        "representative_seed_number", "endpoint_multiplicity_q10"
      ), drop = FALSE
    ]
    names(eligible)[names(eligible) == "representative_seed_number"] <-
      "seed_number"
    endpoint_qc <- merge(
      endpoint_qc, eligible,
      by = c("pair_id", "pair_label", "model_context", "seed_number"),
      sort = FALSE
    )
  }
  endpoint_qc <- endpoint_qc[
    endpoint_qc$model_context == model_context & endpoint_qc$operator_qc_pass,
    , drop = FALSE
  ]
  endpoint_qc$cache_path <- vapply(
    endpoint_qc$cache_path, f7ab_resolve_cache_path, character(1L), paths = paths
  )
  represented <- tapply(
    endpoint_qc$endpoint_multiplicity_q10, endpoint_qc$pair_id, sum
  )
  if (length(represented) != f7r_family_count() || any(represented != 50L)) {
    stop("Figure 7A validation requires exactly 50 q10 seed records per pair in ",
         model_context)
  }

  oxygen_anchors <- c(0, 0.5, 5)
  p_grid <- sort(unique(surface_summary$p_misseg))
  if (length(p_grid) != 60L) {
    stop("Figure 7A validation expected a 60-point p_misseg grid in ",
         model_context)
  }
  p_anchors <- p_grid[c(1L, ceiling(length(p_grid) / 2), length(p_grid))]
  surface_rows <- list()
  trajectory_rows <- list()
  surface_index <- 0L
  trajectory_index <- 0L
  for (index in seq_len(nrow(endpoint_qc))) {
    metadata <- endpoint_qc[index, , drop = FALSE]
    object <- readRDS(metadata$cache_path[[1L]])
    for (oxygen in oxygen_anchors) {
      trajectory <- f7ab_exact_rows(object$trajectory, oxygen = oxygen)
      if (nrow(trajectory) != 1L) stop("Missing Figure 7A trajectory anchor.")
      trajectory_index <- trajectory_index + 1L
      trajectory_rows[[trajectory_index]] <- data.frame(
        model_context = model_context,
        pair_id = metadata$pair_id,
        pair_label = metadata$pair_label,
        seed_number = metadata$seed_number,
        endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10,
        O2_pct = oxygen,
        fitted_p_misseg = trajectory$fitted_p_misseg,
        stringsAsFactors = FALSE
      )
      for (p_value in p_anchors) {
        surface <- f7ab_exact_rows(
          object$surface, oxygen = oxygen, p_misseg = p_value
        )
        if (nrow(surface) != 1L) stop("Missing Figure 7A surface anchor.")
        surface_index <- surface_index + 1L
        surface_rows[[surface_index]] <- data.frame(
          model_context = model_context,
          pair_id = metadata$pair_id,
          pair_label = metadata$pair_label,
          seed_number = metadata$seed_number,
          endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10,
          O2_pct = oxygen,
          p_misseg = p_value,
          dominant_mean_ploidy = surface$dominant_mean_ploidy,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  surface_rows <- do.call(rbind, surface_rows)
  trajectory_rows <- do.call(rbind, trajectory_rows)
  direct_surface <- do.call(rbind, lapply(split(
    surface_rows,
    interaction(
      surface_rows$model_context, surface_rows$pair_id,
      surface_rows$O2_pct, surface_rows$p_misseg, drop = TRUE
    )
  ), function(group) data.frame(
    model_context = group$model_context[[1L]],
    pair_id = group$pair_id[[1L]], pair_label = group$pair_label[[1L]],
    O2_pct = group$O2_pct[[1L]], p_misseg = group$p_misseg[[1L]],
    n_seed_direct = sum(group$endpoint_multiplicity_q10),
    direct_arithmetic_mean = stats::weighted.mean(
      group$dominant_mean_ploidy, group$endpoint_multiplicity_q10
    ),
    stringsAsFactors = FALSE
  )))
  published_surface <- surface_summary[
    surface_summary$cutoff == "q10" &
      f7ab_in_numeric_set(surface_summary$O2_pct, oxygen_anchors) &
      f7ab_in_numeric_set(surface_summary$p_misseg, p_anchors),
    c(
      "pair_id", "pair_label", "O2_pct", "p_misseg", "n_seed",
      "dominant_mean_ploidy_mean", "dominant_mean_ploidy_median"
    ), drop = FALSE
  ]
  checked_surface <- merge(
    direct_surface, published_surface,
    by = c("pair_id", "pair_label", "O2_pct", "p_misseg"),
    all = TRUE, sort = FALSE
  )
  checked_surface$model_context <- model_context
  checked_surface$statistic <- "Figure 7A tile"
  checked_surface$published_arithmetic_mean <-
    checked_surface$dominant_mean_ploidy_mean
  checked_surface$absolute_error <- abs(
    checked_surface$direct_arithmetic_mean -
      checked_surface$published_arithmetic_mean
  )
  checked_surface$passed <- checked_surface$n_seed == 50L &
    checked_surface$n_seed_direct == 50L &
    is.finite(checked_surface$absolute_error) &
    checked_surface$absolute_error <= 1e-12

  direct_trajectory <- do.call(rbind, lapply(split(
    trajectory_rows,
    interaction(
      trajectory_rows$model_context, trajectory_rows$pair_id,
      trajectory_rows$O2_pct, drop = TRUE
    )
  ), function(group) data.frame(
    model_context = group$model_context[[1L]],
    pair_id = group$pair_id[[1L]], pair_label = group$pair_label[[1L]],
    O2_pct = group$O2_pct[[1L]],
    n_seed_direct = sum(group$endpoint_multiplicity_q10),
    direct_arithmetic_mean = stats::weighted.mean(
      group$fitted_p_misseg, group$endpoint_multiplicity_q10
    ),
    stringsAsFactors = FALSE
  )))
  published_trajectory <- trajectory_summary[
    trajectory_summary$cutoff == "q10" &
      f7ab_in_numeric_set(trajectory_summary$O2_pct, oxygen_anchors),
    c(
      "pair_id", "pair_label", "O2_pct", "n_seed",
      "fitted_p_misseg_mean", "fitted_p_misseg_median"
    ), drop = FALSE
  ]
  checked_trajectory <- merge(
    direct_trajectory, published_trajectory,
    by = c("pair_id", "pair_label", "O2_pct"),
    all = TRUE, sort = FALSE
  )
  checked_trajectory$model_context <- model_context
  checked_trajectory$statistic <- "Figure 7A fitted-p_misseg line"
  checked_trajectory$p_misseg <- NA_real_
  checked_trajectory$dominant_mean_ploidy_median <- NA_real_
  checked_trajectory$published_arithmetic_mean <-
    checked_trajectory$fitted_p_misseg_mean
  checked_trajectory$absolute_error <- abs(
    checked_trajectory$direct_arithmetic_mean -
      checked_trajectory$published_arithmetic_mean
  )
  checked_trajectory$passed <- checked_trajectory$n_seed == 50L &
    checked_trajectory$n_seed_direct == 50L &
    is.finite(checked_trajectory$absolute_error) &
    checked_trajectory$absolute_error <= 1e-12

  common <- c(
    "statistic", "model_context", "pair_id", "pair_label", "O2_pct",
    "p_misseg", "n_seed", "direct_arithmetic_mean",
    "published_arithmetic_mean", "absolute_error", "passed"
  )
  rbind(
    checked_surface[, common, drop = FALSE],
    checked_trajectory[, common, drop = FALSE]
  )
}

f7ab_inverse_anchor_rows <- function(paths, model_context) {
  prefix <- if (identical(model_context, "in vivo")) {
    "figure7"
  } else {
    "figure7_invitro"
  }
  inverse_summary <- f7r_read_tsv(file.path(
    paths$figure7, paste0(prefix, "_inverse_response_summary.tsv")
  ))
  inverse_qc <- f7r_read_tsv(file.path(
    paths$figure7, paste0(prefix, "_inverse_endpoint_qc.tsv")
  ))
  inverse_qc <- inverse_qc[inverse_qc$inverse_qc_pass, , drop = FALSE]
  inverse_qc$cache_path <- vapply(
    inverse_qc$cache_path, f7ab_resolve_cache_path, character(1L), paths = paths
  )
  represented <- tapply(
    inverse_qc$endpoint_multiplicity_q10, inverse_qc$pair_id, sum
  )
  if (length(represented) != f7r_family_count() || any(represented != 50L)) {
    stop("Figure 7B validation requires multiplicity 50 per pair in ",
         model_context)
  }
  oxygen_anchors <- c(0, 1, 5)
  target_anchors <- c(2, 4, 6)
  endpoint_rows <- list()
  row_index <- 0L
  for (index in seq_len(nrow(inverse_qc))) {
    metadata <- inverse_qc[index, , drop = FALSE]
    object <- readRDS(metadata$cache_path[[1L]])
    for (oxygen in oxygen_anchors) for (target in target_anchors) {
      inverse <- f7ab_exact_rows(
        object$inverse, oxygen = oxygen, target_ploidy = target
      )
      if (nrow(inverse) != 1L) stop("Missing Figure 7B inverse anchor.")
      row_index <- row_index + 1L
      endpoint_rows[[row_index]] <- data.frame(
        model_context = model_context,
        pair_id = metadata$pair_id,
        pair_label = metadata$pair_label,
        O2_pct = oxygen,
        target_ploidy = target,
        endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10,
        n_solution = inverse$n_solution,
        p_unique = inverse$p_unique,
        stringsAsFactors = FALSE
      )
    }
  }
  endpoint_rows <- do.call(rbind, endpoint_rows)
  groups <- split(
    endpoint_rows,
    interaction(
      endpoint_rows$pair_id, endpoint_rows$O2_pct,
      endpoint_rows$target_ploidy, drop = TRUE
    )
  )
  direct <- do.call(rbind, lapply(groups, function(group) {
    unique_solution <- group$n_solution == 1L & is.finite(group$p_unique)
    data.frame(
      model_context = model_context,
      pair_id = group$pair_id[[1L]],
      pair_label = group$pair_label[[1L]],
      O2_pct = group$O2_pct[[1L]],
      target_ploidy = group$target_ploidy[[1L]],
      n_seed = sum(group$endpoint_multiplicity_q10),
      n_seed_unique_solution = sum(
        group$endpoint_multiplicity_q10[unique_solution]
      ),
      direct_arithmetic_mean = if (any(unique_solution)) {
        stats::weighted.mean(
          group$p_unique[unique_solution],
          group$endpoint_multiplicity_q10[unique_solution]
        )
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  published <- inverse_summary[
    f7ab_in_numeric_set(inverse_summary$O2_pct, oxygen_anchors) &
      f7ab_in_numeric_set(inverse_summary$target_ploidy, target_anchors),
    c(
      "pair_id", "pair_label", "O2_pct", "target_ploidy",
      "inverse_class", "n_seed", "n_seed_unique_solution",
      "p_unique_mean", "p_unique_median", "p_display"
    ), drop = FALSE
  ]
  checked <- merge(
    direct, published,
    by = c("pair_id", "pair_label", "O2_pct", "target_ploidy"),
    suffixes = c("_direct", "_published"), all = TRUE, sort = FALSE
  )
  checked$statistic <- "Figure 7B required p_misseg"
  checked$model_context <- model_context
  checked$p_misseg <- NA_real_
  checked$n_seed <- checked$n_seed_published
  checked$published_arithmetic_mean <- checked$p_unique_mean
  checked$absolute_error <- abs(
    checked$direct_arithmetic_mean - checked$published_arithmetic_mean
  )
  both_missing <- !is.finite(checked$direct_arithmetic_mean) &
    !is.finite(checked$published_arithmetic_mean)
  mean_matches <- both_missing |
    (is.finite(checked$absolute_error) & checked$absolute_error <= 1e-12)
  display_matches <- ifelse(
    checked$inverse_class == "stable unique inverse",
    is.finite(checked$p_display) &
      abs(checked$p_display - checked$p_unique_mean) <= 1e-12,
    !is.finite(checked$p_display)
  )
  checked$passed <- checked$n_seed_direct == 50L &
    checked$n_seed_published == 50L &
    checked$n_seed_unique_solution_direct ==
      checked$n_seed_unique_solution_published &
    mean_matches & display_matches
  checked[, c(
    "statistic", "model_context", "pair_id", "pair_label", "O2_pct",
    "p_misseg", "target_ploidy", "inverse_class", "n_seed",
    "n_seed_unique_solution_published", "direct_arithmetic_mean",
    "published_arithmetic_mean", "p_display", "absolute_error", "passed"
  ), drop = FALSE]
}

validate_Figure7_AB_mean <- function() {
  workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  paths <- f7r_paths(workspace_root)
  anchor_rows <- f7x_rbind_fill(
    f7ab_surface_anchor_rows(paths, "in vivo"),
    f7ab_surface_anchor_rows(paths, "in vitro"),
    f7ab_inverse_anchor_rows(paths, "in vivo"),
    f7ab_inverse_anchor_rows(paths, "in vitro")
  )
  anchor_path <- f7r_write_tsv(
    anchor_rows,
    file.path(paths$figure7, "figure7_ab_mean_anchor_validation.tsv")
  )

  surface_vivo <- f7r_read_tsv(file.path(
    paths$figure7, "joint_multiseed_surface_summary.tsv"
  ))
  surface_vitro <- f7r_read_tsv(file.path(
    paths$figure7, "joint_multiseed_surface_summary_invitro.tsv"
  ))
  inverse_vivo <- f7r_read_tsv(file.path(
    paths$figure7, "figure7_inverse_response_summary.tsv"
  ))
  inverse_vitro <- f7r_read_tsv(file.path(
    paths$figure7, "figure7_invitro_inverse_response_summary.tsv"
  ))
  surface_fields <- c(
    "cutoff", "n_seed", "dominant_mean_ploidy_mean",
    "dominant_mean_ploidy_median"
  )
  primary_surface <- rbind(
    surface_vivo[
      surface_vivo$cutoff == "q10", surface_fields, drop = FALSE
    ],
    surface_vitro[
      surface_vitro$cutoff == "q10", surface_fields, drop = FALSE
    ]
  )
  primary_inverse <- rbind(inverse_vivo, inverse_vitro)
  stable <- primary_inverse$inverse_class == "stable unique inverse"
  checks <- data.frame(
    check = c(
      "all_anchor_means_match_endpoint_level_recomputation",
      "primary_surface_rows_use_50_q10_seed_records",
      "primary_surface_mean_is_finite",
      "primary_surface_median_audit_column_is_retained",
      "inverse_rows_represent_50_q10_seed_records",
      "inverse_display_uses_mean_only_for_stable_unique_cells",
      "inverse_mean_available_exactly_when_unique_solution_exists"
    ),
    observed = c(
      all(anchor_rows$passed),
      all(primary_surface$n_seed == 50L),
      all(is.finite(primary_surface$dominant_mean_ploidy_mean)),
      "dominant_mean_ploidy_median" %in% names(primary_surface),
      all(primary_inverse$n_seed == 50L),
      all(
        is.finite(primary_inverse$p_display) == stable &
          (!stable | abs(
            primary_inverse$p_display - primary_inverse$p_unique_mean
          ) <= 1e-12)
      ),
      all(
        is.finite(primary_inverse$p_unique_mean) ==
          (primary_inverse$n_seed_unique_solution > 0L)
      )
    ),
    expected = TRUE,
    stringsAsFactors = FALSE
  )
  checks$passed <- checks$observed == checks$expected
  validation_path <- f7r_write_tsv(
    checks, file.path(paths$figure7, "figure7_ab_mean_validation.tsv")
  )
  if (!all(checks$passed)) {
    stop(
      "Figure 7A/B arithmetic-mean validation failed: ",
      paste(checks$check[!checks$passed], collapse = ", ")
    )
  }
  message("Figure 7A/B arithmetic-mean validation passed.")
  invisible(list(validation = validation_path, anchors = anchor_path))
}

if (sys.nframe() == 0L) validate_Figure7_AB_mean()
