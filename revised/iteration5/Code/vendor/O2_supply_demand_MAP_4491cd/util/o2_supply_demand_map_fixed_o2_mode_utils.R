#!/usr/bin/env Rscript

# Canonical fixed-O2 attractor-mode semantics shared by simulation and analysis.

fixo2_mode_threshold <- function() 2


fixo2_mode_regimes <- function() {
  c(
    mode1 = "mode1_attractor_dominant_ploidy_ge_2",
    mode2 = "mode2_attractor_dominant_ploidy_lt_2"
  )
}


fixo2_o2_key <- function(x) {
  vapply(x, function(xx) format(signif(as.numeric(xx), 12), scientific = FALSE, trim = TRUE), character(1))
}


fixo2_assign_attractor_modes <- function(tab, ploidy_col = "dominant_mean_ploidy") {
  if (!nrow(tab)) return(tab)
  if (!ploidy_col %in% names(tab)) stop("Cannot assign FixO2 modes; missing column: ", ploidy_col)
  if ("trajectory_regime" %in% names(tab) && !"source_trajectory_regime" %in% names(tab)) tab$source_trajectory_regime <- tab$trajectory_regime
  if ("mode_label" %in% names(tab) && !"source_mode_label" %in% names(tab)) tab$source_mode_label <- tab$mode_label
  fields <- fixo2_mode_fields(tab[[ploidy_col]])
  replace_cols <- intersect(names(fields), names(tab))
  tab[, replace_cols] <- NULL
  cbind(tab, fields, stringsAsFactors = FALSE)
}


fixo2_attractor_mode_table <- function(attractors) {
  if (!nrow(attractors)) return(data.frame())
  d <- fixo2_assign_attractor_modes(attractors, "dominant_mean_ploidy")
  d$O2_key <- fixo2_o2_key(d$O2_pct)
  keep <- intersect(c(
    "seed_id", "O2_pct", "O2_key", "dominant_mean_ploidy", "trajectory_regime",
    "mode_label", "mode_source", "mode_rule", "mode_threshold_dominant_ploidy",
    "status", "dominant_growth_rate", "spectral_gap", "objective", "delta_objective",
    "in_attractor_o2_grid", "is_mode_reference_o2"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d[order(o2ipa_seed_number(d$seed_id), d$O2_pct), , drop = FALSE]
}


fixo2_attractor_mode_summary_by_seed <- function(mode_by_seed_o2, standard_o2 = c(0, 0.1, 0.5, 1, 2, 5)) {
  if (!nrow(mode_by_seed_o2)) return(data.frame())
  rows <- lapply(split(mode_by_seed_o2, mode_by_seed_o2$seed_id), function(d) {
    d <- d[order(d$O2_pct), , drop = FALSE]
    out <- data.frame(
      seed_id = d$seed_id[[1]],
      n_o2 = nrow(d),
      n_o2_mode1 = sum(d$mode_label == "mode1", na.rm = TRUE),
      n_o2_mode2 = sum(d$mode_label == "mode2", na.rm = TRUE),
      fraction_o2_mode1 = mean(d$mode_label == "mode1", na.rm = TRUE),
      fraction_o2_mode2 = mean(d$mode_label == "mode2", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (O2 in standard_o2) {
      key <- paste0("mode_at_o2_", gsub("[^0-9A-Za-z]+", "p", format(O2, scientific = FALSE, trim = TRUE)))
      hit <- d$mode_label[abs(as.numeric(d$O2_pct) - O2) < 1e-9]
      out[[key]] <- if (length(hit)) hit[[1]] else NA_character_
    }
    out
  })
  out <- do.call(rbind, rows)
  out[order(o2ipa_seed_number(out$seed_id)), , drop = FALSE]
}


fixo2_reference_mode_table <- function(mode_by_seed_o2, mode_reference_o2) {
  if (!nrow(mode_by_seed_o2)) return(data.frame())
  d <- mode_by_seed_o2[abs(as.numeric(mode_by_seed_o2$O2_pct) - mode_reference_o2) < 1e-9, , drop = FALSE]
  if (!nrow(d)) {
    stop(
      "No FixO2 attractor mode rows matched --mode_reference_o2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      ". Include this O2 value in the mode table or allow the workflow to compute it."
    )
  }
  d <- d[order(o2ipa_seed_number(d$seed_id)), , drop = FALSE]
  d <- d[!duplicated(d$seed_id), , drop = FALSE]
  threshold <- if ("mode_threshold_dominant_ploidy" %in% names(d)) d$mode_threshold_dominant_ploidy else rep(fixo2_mode_threshold(), nrow(d))
  out <- data.frame(
    seed_id = d$seed_id,
    mode_reference_o2_pct = as.numeric(d$O2_pct),
    mode_reference_o2_key = fixo2_o2_key(d$O2_pct),
    mode_reference_dominant_mean_ploidy = suppressWarnings(as.numeric(d$dominant_mean_ploidy)),
    trajectory_regime = d$trajectory_regime,
    mode_label = d$mode_label,
    mode_source = "fixed_o2_attractor_dominant_ploidy_at_reference_o2",
    mode_rule = paste0(
      "dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " >= 2 => mode1; dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " < 2 => mode2"
    ),
    mode_threshold_dominant_ploidy = threshold,
    stringsAsFactors = FALSE
  )
  optional_cols <- c(
    status = "mode_reference_status",
    dominant_growth_rate = "mode_reference_dominant_growth_rate",
    spectral_gap = "mode_reference_spectral_gap",
    objective = "objective",
    delta_objective = "delta_objective"
  )
  for (src in names(optional_cols)) {
    if (src %in% names(d)) out[[optional_cols[[src]]]] <- d[[src]]
  }
  out
}
