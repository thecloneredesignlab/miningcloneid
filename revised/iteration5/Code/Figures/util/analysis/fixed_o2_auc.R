#!/usr/bin/env Rscript

# Build a dense fixed-O2, single-parameter AUC heatmap from all 500 separate
# in vivo fitting endpoints. This script performs no refitting.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(scales)
})

script_path <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  } else {
    normalizePath("build_invivo_fixed_o2_auc_heatmap.R", mustWork = FALSE)
  }
})
script_dir <- dirname(script_path)
repo_root <- normalizePath(
  Sys.getenv(
    "HYPOXIA_REPO_ROOT",
    unset = file.path(script_dir, "..", "..", "..", "..", "..", "..")
  ),
  mustWork = TRUE
)

parse_cli <- function(x) {
  out <- list()
  for (arg in x) {
    arg <- sub("^--", "", arg)
    bits <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- bits[[1L]]
    value <- if (length(bits) > 1L) paste(bits[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || !nzchar(as.character(x[[1L]]))) return(default)
  tolower(trimws(as.character(x[[1L]]))) %in% c("true", "t", "1", "yes", "y")
}

as_int <- function(x, default) {
  value <- suppressWarnings(as.integer(x %||% default))
  if (!is.finite(value) || is.na(value)) default else value
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (is.null(args$fit_root) || !nzchar(trimws(args$fit_root))) {
  stop("--fit_root=PATH is required.")
}
fit_root <- normalizePath(args$fit_root, mustWork = TRUE)
output_dir <- normalizePath(args$output_dir %||% script_dir, mustWork = FALSE)
data_dir <- file.path(output_dir, "data")
figure_dir <- file.path(output_dir, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

fixed_o2_script <- normalizePath(
  args$fixed_o2_script %||%
    file.path(
      repo_root,
      "oxygen", "code", "O2_supply_demand_MAP", "simulation",
      "fix_o2_simulation.R"
    ),
  mustWork = TRUE
)
legacy_analysis_root <- normalizePath(
  args$legacy_analysis_root %||%
    file.path(output_dir, "disabled_legacy_validation"),
  mustWork = FALSE
)
n_core_default <- max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
n_core <- max(1L, as_int(args$n_core, n_core_default))
overwrite <- as_bool(args$overwrite, FALSE)
analysis_only <- as_bool(args$analysis_only, TRUE)

o2_grid <- seq(0, 5, length.out = 201L)
stopifnot(length(o2_grid) == 201L, isTRUE(all.equal(diff(o2_grid), rep(0.025, 200L))))

parameter_meta <- data.frame(
  parameter = c(
    "O2_crit", "n_O", "o2_S0", "kappa_O", "eta_o2",
    "lam_max", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "k_clear",
    "p_mis_base", "p_misseg", "k_o_mis", "p_wgd",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  ),
  parameter_group = c(
    rep("Oxygen response and supply-demand", 5L),
    rep("Growth, stress-associated death, and clearance", 6L),
    rep("Chromosome instability and WGD", 4L),
    rep("Post-missegregation survival", 3L)
  ),
  parameter_description = c(
    "Critical oxygen level for the oxygen-linked resource-stress function",
    "Hill exponent controlling oxygen-response steepness",
    "Low-burden effective oxygen supply level",
    "Oxygen-drop amplitude in the supply-demand target",
    "Chromosome-number-weighted oxygen-demand exponent",
    "Maximum division rate",
    "Resource-stress-dependent high-chromosome-count growth damping",
    "Chromosome-number dependence of the resource-stress growth penalty",
    "Stress-associated death scale",
    "Chromosome-number penalty on stress-associated death",
    "Dead-biomass clearance rate",
    "Baseline per-chromosome missegregation probability",
    "Maximum stress-induced per-chromosome missegregation increment",
    "Death-hazard half-saturation scale for stress-induced missegregation",
    "Per-division whole-genome-doubling probability",
    "Maximum per-copy post-missegregation survival factor",
    "Ploidy-dependent post-missegregation viability-loss strength",
    "Exponent controlling ploidy dependence of post-missegregation survival"
  ),
  transformation = c(
    "log10", "identity", "log10", "log10", "log10",
    "log10", "log10", "identity", "log10", "identity", "log10",
    "log10", "log10", "log10", "log10",
    "identity", "log10", "log10"
  ),
  stringsAsFactors = FALSE
)
parameter_meta$parameter_order <- seq_len(nrow(parameter_meta))
group_levels <- unique(parameter_meta$parameter_group)
parameter_meta$parameter_group <- factor(parameter_meta$parameter_group, levels = group_levels)
parameter_names <- parameter_meta$parameter
fwrite(parameter_meta, file.path(data_dir, "parameter_function_groups.tsv"), sep = "\t")

message("Loading fixed-O2 analytical engine: ", fixed_o2_script)
fix_env <- new.env(parent = globalenv())
source(fixed_o2_script, local = fix_env)

message("Reading 500 fitting endpoints directly from: ", fit_root)
inputs <- fix_env$o2ipa_collect_seed_inputs(fit_root, objective_source = "auto")
manifest <- as.data.table(inputs$manifest)
manifest[, seed_number := suppressWarnings(as.integer(sub("^seed", "", seed_id)))]
setorder(manifest, seed_number)
expected_seeds <- paste0("seed", seq_len(500L))
if (nrow(manifest) != 500L || !identical(manifest$seed_id, expected_seeds)) {
  missing <- setdiff(expected_seeds, manifest$seed_id)
  extra <- setdiff(manifest$seed_id, expected_seeds)
  stop(
    "Expected exactly seed1 through seed500. Missing: ",
    paste(missing, collapse = ","),
    "; extra: ",
    paste(extra, collapse = ",")
  )
}
if (any(!manifest$fit_success)) {
  stop("Non-successful fits found: ", paste(manifest$seed_id[!manifest$fit_success], collapse = ","))
}

read_best_parameter_row <- function(seed_id, parameter_file) {
  tab <- fread(parameter_file)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("Malformed best-parameter file: ", parameter_file)
  }
  values <- suppressWarnings(as.numeric(tab$value))
  names(values) <- as.character(tab$parameter)
  missing <- setdiff(parameter_meta$parameter, names(values))
  if (length(missing)) {
    stop("Missing parameters for ", seed_id, ": ", paste(missing, collapse = ", "))
  }
  row <- as.list(values[parameter_meta$parameter])
  row$seed_id <- seed_id
  as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
}

parameter_rows <- lapply(
  seq_len(nrow(manifest)),
  function(i) read_best_parameter_row(manifest$seed_id[[i]], manifest$parameter_file[[i]])
)
best_params <- as.data.table(do.call(rbind, parameter_rows))
for (parameter in parameter_meta$parameter) {
  set(best_params, j = parameter, value = suppressWarnings(as.numeric(best_params[[parameter]])))
}
best_params[, seed_number := suppressWarnings(as.integer(sub("^seed", "", seed_id)))]
best_params <- merge(
  best_params,
  manifest[, .(
    seed_id, objective, objective_source, objective_data, objective_burden,
    objective_ploidy, convergence_status, parameter_file
  )],
  by = "seed_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(best_params, seed_number)
if (any(!is.finite(as.matrix(best_params[, ..parameter_names])))) {
  stop("At least one of the 18 fitted parameters is missing or non-finite.")
}
fwrite(best_params, file.path(data_dir, "invivo_best_parameters_500seeds.tsv"), sep = "\t")

attractor_path <- file.path(data_dir, "fixed_o2_dominant_ploidy_201grid.tsv")
if (file.exists(attractor_path) && !overwrite) {
  message("Reusing existing dense attractor table: ", attractor_path)
  attractors <- fread(attractor_path)
} else {
  message(
    "Calculating analytical dominant attractors for 500 seeds x 201 O2 values ",
    "with ", n_core, " workers."
  )
  attractors <- as.data.table(
    fix_env$generate_fixo2_attractor_mode_table(
      run_dir = fit_root,
      o2_values = o2_grid,
      seed_ids = expected_seeds,
      n_workers = n_core
    )
  )
  drop_abstract_columns <- intersect(
    c(
      "trajectory_regime", "mode_label", "mode_source", "mode_rule",
      "mode_threshold_dominant_ploidy", "is_mode_reference_o2"
    ),
    names(attractors)
  )
  if (length(drop_abstract_columns)) attractors[, (drop_abstract_columns) := NULL]
  attractors[, ploidy_class := fifelse(
    is.finite(dominant_mean_ploidy) & dominant_mean_ploidy >= 2,
    "dominant_mean_ploidy \u2265 2",
    fifelse(
      is.finite(dominant_mean_ploidy) & dominant_mean_ploidy < 2,
      "dominant_mean_ploidy < 2",
      NA_character_
    )
  )]
  attractors[, seed_number := suppressWarnings(as.integer(sub("^seed", "", seed_id)))]
  setorder(attractors, seed_number, O2_pct)
  fwrite(attractors, attractor_path, sep = "\t")
}

required_attractor_cols <- c("seed_id", "O2_pct", "status", "dominant_mean_ploidy", "ploidy_class")
missing_attractor_cols <- setdiff(required_attractor_cols, names(attractors))
if (length(missing_attractor_cols)) {
  stop("Dense attractor table is missing: ", paste(missing_attractor_cols, collapse = ", "))
}
if (nrow(attractors) != 500L * 201L) {
  stop("Expected 100,500 attractor rows; observed ", nrow(attractors), ".")
}
if (uniqueN(attractors$seed_id) != 500L || uniqueN(attractors$O2_pct) != 201L) {
  stop("Dense attractor table does not contain 500 seeds and 201 O2 values.")
}

transform_parameter <- function(x, transformation) {
  x <- suppressWarnings(as.numeric(x))
  if (identical(transformation, "log10")) {
    if (any(!is.finite(x) | x <= 0)) return(rep(NA_real_, length(x)))
    return(log10(x))
  }
  x
}

rank_auc <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- as.integer(y[keep])
  score <- as.numeric(score[keep])
  n_positive <- sum(y == 1L)
  n_negative <- sum(y == 0L)
  if (n_positive == 0L || n_negative == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

best_params_index <- match(expected_seeds, best_params$seed_id)
parameter_natural <- as.matrix(best_params[best_params_index, ..parameter_names])
storage.mode(parameter_natural) <- "double"
parameter_transformed <- parameter_natural
for (j in seq_len(nrow(parameter_meta))) {
  parameter_transformed[, j] <- transform_parameter(
    parameter_natural[, j],
    parameter_meta$transformation[[j]]
  )
}

auc_rows <- vector("list", length(o2_grid) * nrow(parameter_meta))
class_rows <- vector("list", length(o2_grid))
row_index <- 0L
for (o2_index in seq_along(o2_grid)) {
  o2 <- o2_grid[[o2_index]]
  current <- attractors[abs(O2_pct - o2) < 1e-10]
  current <- current[match(expected_seeds, seed_id)]
  y <- fifelse(
    is.finite(current$dominant_mean_ploidy),
    as.integer(current$dominant_mean_ploidy >= 2),
    NA_integer_
  )
  n_positive <- sum(y == 1L, na.rm = TRUE)
  n_negative <- sum(y == 0L, na.rm = TRUE)
  n_missing <- sum(is.na(y))
  class_rows[[o2_index]] <- data.frame(
    O2_pct = o2,
    n_dominant_mean_ploidy_ge_2 = n_positive,
    n_dominant_mean_ploidy_lt_2 = n_negative,
    n_missing_or_failed = n_missing,
    auc_defined = n_positive > 0L && n_negative > 0L,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(nrow(parameter_meta))) {
    row_index <- row_index + 1L
    raw_auc <- rank_auc(y, parameter_transformed[, j])
    auc_strength <- if (is.finite(raw_auc)) max(raw_auc, 1 - raw_auc) else NA_real_
    signed_auc <- if (!is.finite(raw_auc)) {
      NA_real_
    } else if (abs(raw_auc - 0.5) < 1e-12) {
      0
    } else if (raw_auc > 0.5) {
      auc_strength
    } else {
      -auc_strength
    }
    direction <- if (!is.finite(raw_auc)) {
      "AUC undefined: only one dominant-mean-ploidy class present"
    } else if (abs(raw_auc - 0.5) < 1e-12) {
      "No directional separation"
    } else if (raw_auc > 0.5) {
      "Higher parameter values associated with dominant_mean_ploidy \u2265 2"
    } else {
      "Higher parameter values associated with dominant_mean_ploidy < 2"
    }
    natural_values <- parameter_natural[, j]
    auc_rows[[row_index]] <- data.frame(
      O2_pct = o2,
      parameter = parameter_meta$parameter[[j]],
      parameter_group = as.character(parameter_meta$parameter_group[[j]]),
      transformation = parameter_meta$transformation[[j]],
      n_total = length(y),
      n_dominant_mean_ploidy_ge_2 = n_positive,
      n_dominant_mean_ploidy_lt_2 = n_negative,
      n_missing_or_failed = n_missing,
      raw_auc_ge_2_positive = raw_auc,
      discriminative_auc = auc_strength,
      signed_auc_for_color = signed_auc,
      direction = direction,
      mean_parameter_ge_2 = if (n_positive) mean(natural_values[y == 1L], na.rm = TRUE) else NA_real_,
      mean_parameter_lt_2 = if (n_negative) mean(natural_values[y == 0L], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

auc_table <- as.data.table(do.call(rbind, auc_rows))
class_counts <- as.data.table(do.call(rbind, class_rows))
fwrite(auc_table, file.path(data_dir, "parameter_auc_by_o2.tsv"), sep = "\t")
fwrite(class_counts, file.path(data_dir, "ploidy_class_counts_by_o2.tsv"), sep = "\t")

validation_rows <- list()
legacy_attractor_path <- file.path(
  legacy_analysis_root,
  "01_fixed_o2", "FixO2_invivo_500seed", "attractors", "tables",
  "fixed_o2_attractor_mode_by_seed_o2.tsv"
)
if (file.exists(legacy_attractor_path)) {
  legacy <- fread(legacy_attractor_path)
  dense_validation <- attractors[, .(
    seed_id,
    O2_key_validation = sprintf("%.6f", O2_pct),
    dense_value = dominant_mean_ploidy
  )]
  legacy_validation <- legacy[, .(
    seed_id,
    O2_key_validation = sprintf("%.6f", O2_pct),
    legacy_value = dominant_mean_ploidy
  )]
  overlap <- merge(
    dense_validation,
    legacy_validation,
    by = c("seed_id", "O2_key_validation"),
    all = FALSE
  )
  diff <- abs(overlap$dense_value - overlap$legacy_value)
  validation_rows[[length(validation_rows) + 1L]] <- data.frame(
    check = "Dense attractor values reproduce legacy 0.05-grid values",
    expected = "50,500 overlapping rows; max absolute difference <= 1e-10",
    observed = paste0(
      nrow(overlap), " rows; max absolute difference=",
      format(max(diff, na.rm = TRUE), scientific = TRUE)
    ),
    passed = nrow(overlap) == 50500L && max(diff, na.rm = TRUE) <= 1e-10,
    source = legacy_attractor_path,
    stringsAsFactors = FALSE
  )
}

for (reference_o2 in c(0, 1, 5)) {
  reference_slug <- format(reference_o2, scientific = FALSE, trim = TRUE)
  legacy_auc_path <- file.path(
    legacy_analysis_root,
    "02_parameter_landscape_clustering", "mode_parameter_contribution",
    paste0("reference_o2_", reference_slug), "auc_plots",
    "mode_parameter_auc_by_feature.csv"
  )
  if (!file.exists(legacy_auc_path)) next
  legacy_auc <- fread(legacy_auc_path)
  legacy_auc <- legacy_auc[feature_type == "main", .(parameter = feature_name, legacy_auc = auc_mode1_higher)]
  current_auc <- auc_table[
    abs(O2_pct - reference_o2) < 1e-10,
    .(parameter, current_auc = raw_auc_ge_2_positive)
  ]
  comparison <- merge(current_auc, legacy_auc, by = "parameter", all = FALSE)
  auc_diff <- abs(comparison$current_auc - comparison$legacy_auc)
  validation_rows[[length(validation_rows) + 1L]] <- data.frame(
    check = paste0("AUC values reproduce legacy reference O2=", reference_o2),
    expected = "18 parameters; max absolute difference <= 1e-12",
    observed = paste0(
      nrow(comparison), " parameters; max absolute difference=",
      format(max(auc_diff, na.rm = TRUE), scientific = TRUE)
    ),
    passed = nrow(comparison) == 18L && max(auc_diff, na.rm = TRUE) <= 1e-12,
    source = legacy_auc_path,
    stringsAsFactors = FALSE
  )
}

validation_rows[[length(validation_rows) + 1L]] <- data.frame(
  check = "Output dimensions",
  expected = "500 seeds; 201 O2 values; 100,500 attractor rows; 3,618 AUC rows",
  observed = paste0(
    uniqueN(attractors$seed_id), " seeds; ",
    uniqueN(attractors$O2_pct), " O2 values; ",
    nrow(attractors), " attractor rows; ",
    nrow(auc_table), " AUC rows"
  ),
  passed = uniqueN(attractors$seed_id) == 500L &&
    uniqueN(attractors$O2_pct) == 201L &&
    nrow(attractors) == 100500L &&
    nrow(auc_table) == 3618L,
  source = attractor_path,
  stringsAsFactors = FALSE
)
validation <- rbindlist(validation_rows, fill = TRUE)
fwrite(validation, file.path(data_dir, "validation_summary.tsv"), sep = "\t")
if (any(!validation$passed)) {
  stop(
    "At least one validation check failed. See ",
    file.path(data_dir, "validation_summary.tsv")
  )
}

png_path <- file.path(figure_dir, "fixed_o2_parameter_auc_heatmap.png")
pdf_path <- file.path(figure_dir, "fixed_o2_parameter_auc_heatmap.pdf")
svg_path <- file.path(figure_dir, "fixed_o2_parameter_auc_heatmap.svg")

if (!analysis_only) {
plot_data <- copy(auc_table)
plot_data <- merge(
  plot_data,
  parameter_meta[, c("parameter", "parameter_group", "parameter_order")],
  by = c("parameter", "parameter_group"),
  all.x = TRUE,
  sort = FALSE
)
plot_data[, parameter := factor(parameter, levels = rev(parameter_meta$parameter))]
plot_data[, parameter_group := factor(parameter_group, levels = group_levels)]
parameter_annotation <- unique(plot_data[, .(parameter, parameter_group)])
group_palette <- c(
  "Oxygen response and supply-demand" = "#4E79A7",
  "Growth, stress-associated death, and clearance" = "#D89C2B",
  "Chromosome instability and WGD" = "#B07AA1",
  "Post-missegregation survival" = "#7A8B3A"
)

heatmap <- ggplot(plot_data, aes(x = O2_pct, y = parameter, fill = signed_auc_for_color)) +
  geom_tile(width = 0.025, height = 0.92) +
  geom_hline(
    yintercept = c(3.5, 7.5, 13.5),
    linewidth = 0.45,
    color = "#505050"
  ) +
  geom_segment(
    data = parameter_annotation,
    aes(
      x = -0.074,
      xend = -0.016,
      y = parameter,
      yend = parameter,
      color = parameter_group
    ),
    inherit.aes = FALSE,
    linewidth = 5.2,
    lineend = "butt"
  ) +
  scale_x_continuous(
    name = expression(paste("Fixed ", O[2], " concentration (%)")),
    limits = c(-0.082, 5.0125),
    breaks = seq(0, 5, by = 0.5),
    expand = c(0, 0)
  ) +
  scale_y_discrete(name = NULL, drop = TRUE) +
  scale_color_manual(
    name = "Parameter function",
    values = group_palette,
    breaks = names(group_palette),
    guide = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 4)
    )
  ) +
  scale_fill_gradientn(
    name = "Discriminative AUC",
    colors = c("#2166AC", "#F7F7F7", "#F7F7F7", "#F7F7F7", "#B2182B"),
    values = rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1),
    breaks = c(-1, -0.75, 0, 0.75, 1),
    labels = c("1.00", "0.75", "0.50", "0.75", "1.00"),
    na.value = "#D9D9D9",
    oob = squish,
    guide = guide_colorbar(
      order = 2,
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(12, "cm"),
      barheight = grid::unit(0.45, "cm"),
      ticks.colour = "#303030",
      frame.colour = "#303030"
    )
  ) +
  labs(
    title = expression(paste("Oxygen-dependent parameter separation of dominant-ploidy outcomes")),
    subtitle = paste0(
      "500 separate in vivo fitted parameter sets. Blue: higher values associated with ",
      "dominant_mean_ploidy < 2; red: higher values associated with ",
      "dominant_mean_ploidy \u2265 2."
    ),
    caption = paste0(
      "Color intensity is max(AUC, 1\u2212AUC); legend labels show non-negative AUC magnitude. ",
      "Near-white indicates AUC \u2248 0.5. Gray indicates an undefined AUC because only one ",
      "dominant-mean-ploidy class was present."
    )
  ) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#202020"),
    plot.subtitle = element_text(size = 10.5, color = "#303030", margin = margin(b = 8)),
    plot.caption = element_text(size = 8.5, color = "#4A4A4A", hjust = 0, margin = margin(t = 8)),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    axis.text.x = element_text(size = 9, color = "#202020"),
    axis.text.y = element_text(size = 9.5, color = "#202020"),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#707070", linewidth = 0.35),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y = grid::unit(2, "mm"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 16, 10, 10)
  )

figure_width_in <- 16
# At fixed title/legend sizes, 8.05 inches reduces the physical heatmap-row
# height to approximately two thirds of the previous 10.5-inch export.
figure_height_in <- 8.05
ggsave(
  png_path, heatmap,
  width = figure_width_in, height = figure_height_in,
  units = "in", dpi = 300, bg = "white"
)
ggsave(
  pdf_path, heatmap,
  width = figure_width_in, height = figure_height_in,
  units = "in", device = cairo_pdf, bg = "white"
)
ggsave(
  svg_path, heatmap,
  width = figure_width_in, height = figure_height_in,
  units = "in", device = svglite::svglite, bg = "white"
)
}

git_sha <- tryCatch(
  system2("git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)
run_arguments <- data.frame(
  argument = c(
    "fit_root", "fixed_o2_script", "output_dir", "n_seed", "o2_min", "o2_max",
    "o2_n", "o2_step", "n_parameter", "positive_class", "negative_class",
    "auc_definition", "signed_color_definition", "n_core", "overwrite",
    "analysis_only", "repository_head", "run_timestamp"
  ),
  value = c(
    fit_root, fixed_o2_script, output_dir, "500", "0", "5", "201", "0.025", "18",
    "dominant_mean_ploidy \u2265 2", "dominant_mean_ploidy < 2",
    "rank-based AUC with average ranks for ties",
    "sign(raw_auc-0.5) * max(raw_auc,1-raw_auc); exact raw_auc=0.5 is encoded as 0",
    as.character(n_core), as.character(overwrite), as.character(analysis_only),
    git_sha[[1L]],
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ),
  stringsAsFactors = FALSE
)
fwrite(run_arguments, file.path(data_dir, "run_arguments.tsv"), sep = "\t")

provenance_paths <- c(
  fit_root = fit_root,
  fixed_o2_engine = fixed_o2_script,
  legacy_dense_attractor_table = legacy_attractor_path,
  script = script_path
)
provenance <- data.frame(
  source = names(provenance_paths),
  path = unname(provenance_paths),
  exists = file.exists(provenance_paths) | dir.exists(provenance_paths),
  md5 = vapply(
    provenance_paths,
    function(path) {
      info <- file.info(path)
      if (file.exists(path) && !isTRUE(info$isdir[[1L]])) {
        unname(tools::md5sum(path)[[1L]])
      } else {
        NA_character_
      }
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)
fwrite(provenance, file.path(data_dir, "source_file_provenance.tsv"), sep = "\t")

message("Completed dense fixed-O2 AUC heatmap workflow.")
message("  Attractor rows: ", nrow(attractors))
message("  AUC rows: ", nrow(auc_table))
if (!analysis_only) message("  Figure: ", png_path)
