# Legacy fixed-O2 pipeline orchestration functions.
# Kept only to preserve the historical all-in-one CLI while domain work lives
# in separate modules.

fo2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}


fixo2_attractor_o2_grid <- function(args) {
  sort(unique(fixo2_as_num_vec(args$attractor_o2_grid, seq(0, 5, by = 0.05))))
}


fixo2_format_o2_list <- function(x, max_n = 18L) {
  x <- sort(unique(as.numeric(x)))
  labs <- format(x, scientific = FALSE, trim = TRUE)
  if (length(labs) > max_n) {
    labs <- c(labs[seq_len(max_n)], paste0("... (", length(x), " total)"))
  }
  paste(labs, collapse = ",")
}








fo2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs", "report")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}


fo2_metric_map <- function(path) {
  if (!file.exists(path)) return(list())
  tab <- tryCatch(fo2_read_tsv(path), error = function(e) data.frame())
  if (!all(c("metric", "value") %in% names(tab))) return(list())
  vals <- as.list(tab$value)
  names(vals) <- tab$metric
  vals
}


fo2_map_num <- function(map, keys) {
  for (key in keys) {
    val <- suppressWarnings(as.numeric(map[[key]]))
    if (length(val) && is.finite(val[[1]])) return(val[[1]])
  }
  NA_real_
}


fo2_map_chr <- function(map, keys) {
  for (key in keys) {
    val <- map[[key]]
    if (!is.null(val) && length(val) && !is.na(val[[1]]) && nzchar(as.character(val[[1]]))) return(as.character(val[[1]]))
  }
  NA_character_
}


fo2_seed_manifest_extras <- function(manifest) {
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    seed <- manifest$seed_id[[i]]
    seed_dir <- manifest$seed_dir[[i]]
    m <- if (!is.na(seed_dir) && nzchar(seed_dir)) fo2_metric_map(file.path(seed_dir, "fit_summary.tsv")) else list()
    data.frame(
      seed_id = seed,
      optimizer_iter_completed = fo2_map_num(m, c("optimizer_iter_completed", "deoptim_iter_completed", "itermax")),
      optimizer_stop_reason = fo2_map_chr(m, c("deoptim_stop_reason", "optimizer_stop_reason", "optimizer_local_convergence")),
      optimizer_local_objective = fo2_map_num(m, c("optimizer_local_objective")),
      optimizer_deoptim_objective = fo2_map_num(m, c("optimizer_deoptim_objective")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


fo2_write_spectral_gap_outputs <- function(attractors, out_dir, generate_figures = TRUE) {
  gap_by_seed <- fo2_spectral_gap_by_seed(attractors)
  gap_summary <- fo2_spectral_gap_summary(gap_by_seed)
  composite <- fo2_ploidy_gap_reliability_composite_table(gap_by_seed)
  violin <- fo2_ploidy_gap_reliability_violin_table(gap_by_seed)
  fo2_write_tsv(gap_by_seed, file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_by_seed.tsv"))
  fo2_write_tsv(gap_summary, file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_regime_summary.tsv"))
  fo2_write_tsv(composite, file.path(out_dir, "tables", "fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv"))
  fo2_write_tsv(violin, file.path(out_dir, "tables", "fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv"))
  if (isTRUE(generate_figures)) {
    dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
    fo2_plot_spectral_gap_outputs(gap_by_seed, gap_summary, file.path(out_dir, "figures"))
    fo2_plot_ploidy_gap_reliability_composite(composite, gap_summary, file.path(out_dir, "figures"))
    fo2_plot_ploidy_gap_reliability_violin(violin, file.path(out_dir, "figures"))
  }
  invisible(list(by_seed = gap_by_seed, summary = gap_summary))
}


cf2_write_tsv <- fo2_write_tsv
cf2_read_tsv <- fo2_read_tsv


cf2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs", "report")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}


cf2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}


cf2_default_time_grid <- vf2_default_time_grid


vf2_write_tsv <- fo2_write_tsv


vf2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}


vf2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}


vf2_as_seed_vec <- function(x, default) {
  vals <- o2ipa_split_csv(x, default)
  vals <- unique(o2ipa_norm_seed(vals))
  if (length(vals)) vals else default
}


fixo2_arg_nonempty <- function(x) {
  !is.null(x) && length(x) && !is.na(x[[1]]) && nzchar(trimws(as.character(x[[1]])))
}


vf2_seed_modes <- function(seed_ids, mode_arg = NULL) {
  modes <- o2ipa_split_csv(mode_arg, character())
  if (length(modes) == length(seed_ids)) {
    mode_map <- setNames(modes, seed_ids)
  } else {
    mode_map <- character()
  }
  out <- rep("", length(seed_ids))
  names(out) <- seed_ids
  has_mode <- seed_ids %in% names(mode_map) & nzchar(mode_map[seed_ids])
  out[has_mode] <- unname(mode_map[seed_ids[has_mode]])
  out
}


vf2_seed_labels <- function(seed_ids, mode_arg = NULL) {
  seed_modes <- vf2_seed_modes(seed_ids, mode_arg)
  labels <- seed_ids
  has_mode <- nzchar(seed_modes)
  labels[has_mode] <- paste0(seed_ids[has_mode], " (", seed_modes[has_mode], ")")
  labels
}


vf2_as_int_vec <- function(x, default) {
  vals <- suppressWarnings(as.integer(o2ipa_split_csv(x, as.character(default))))
  vals <- vals[!is.na(vals) & vals > 0L]
  if (length(vals)) unique(vals) else default
}


vf2_num_path_tag <- num_path_tag


fixo2_write_tsv <- fo2_write_tsv
fixo2_read_tsv <- fo2_read_tsv


fixo2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}


fixo2_default_run_dir <- function(repo_root) {
  default_hpc_run_dir <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed"
  local_run_dir <- file.path(repo_root, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
  if (dir.exists(default_hpc_run_dir)) default_hpc_run_dir else local_run_dir
}


fixo2_default_label_file <- function(repo_root) {
  ""
}


fixo2_resolve_repo_path <- function(path, repo_root, mustWork = FALSE) {
  if (is.null(path) || !length(path) || is.na(path[[1]]) || !nzchar(as.character(path[[1]]))) return(path)
  path <- as.character(path[[1]])
  if (identical(path, "~")) {
    path <- repo_root
  } else if (startsWith(path, "~/")) {
    path <- file.path(repo_root, substring(path, 3L))
  } else if (!grepl("^/", path)) {
    path <- file.path(repo_root, path)
  }
  normalizePath(path, mustWork = mustWork)
}


fixo2_default_simulation_dir <- function(repo_root) {
  file.path(repo_root, "oxygen", "results", "O2_fixed_simulation")
}


fixo2_default_html_report_script <- function(repo_root) {
  file.path(repo_root, "oxygen", "code", "O2_supply_demand_MAP", "report", "render_fixo2_invivo_report.R")
}


# Legacy cross-stage orchestration retained only for compatibility entrypoints.
fixo2_prepare_dirs <- function(out_dir) {
  dirs <- c(
    out_dir,
    file.path(out_dir, "tables"),
    file.path(out_dir, "logs"),
    file.path(out_dir, "report"),
    file.path(out_dir, "attractors"),
    file.path(out_dir, "counterfactual_trajectories"),
    file.path(out_dir, "simulation"),
    file.path(out_dir, "simulation", "analytical_simulation_agreement")
  )
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}


fixo2_write_output_manifest <- function(out_dir) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) {
    manifest <- data.frame()
  } else {
    rel <- substring(files, nchar(normalizePath(out_dir, mustWork = FALSE)) + 2L)
    manifest <- data.frame(
      section = ifelse(grepl("^attractors/", rel), "attractors",
        ifelse(grepl("^counterfactual_trajectories/", rel), "counterfactual_trajectories",
          ifelse(grepl("^simulation/analytical_simulation_agreement/", rel), "analytical_simulation_agreement",
            ifelse(grepl("^simulation/", rel), "simulation", "top_level")))),
      relative_path = rel,
      absolute_path = normalizePath(files, mustWork = FALSE),
      size_bytes = as.numeric(file.info(files)$size),
      stringsAsFactors = FALSE
    )
  }
  fixo2_write_tsv(manifest, file.path(out_dir, "tables", "FixO2_invivo_output_manifest.tsv"))
}


fixo2_normalize_parts <- function(x) {
  parts <- tolower(o2ipa_split_csv(x, "all"))
  if ("all" %in% parts) return(c("attractors", "counterfactual_trajectories", "simulation", "analytical_simulation_agreement"))
  parts <- sub("^counterfactual$", "counterfactual_trajectories", parts)
  parts <- sub("^trajectories$", "counterfactual_trajectories", parts)
  parts <- sub("^validation$", "simulation", parts)
  parts <- sub("^agreement$", "analytical_simulation_agreement", parts)
  parts <- sub("^analytical_agreement$", "analytical_simulation_agreement", parts)
  parts <- sub("^analytical_simulation$", "analytical_simulation_agreement", parts)
  parts <- sub("^scatter$", "analytical_simulation_agreement", parts)
  parts <- sub("^scatters$", "analytical_simulation_agreement", parts)
  valid <- c("attractors", "counterfactual_trajectories", "simulation", "analytical_simulation_agreement")
  bad <- setdiff(parts, valid)
  if (length(bad)) stop("Unknown run_parts value(s): ", paste(bad, collapse = ","))
  unique(parts)
}


fixo2_run_attractors <- function(args, repo_root, run_dir, label_file, out_dir) {
  if (!dir.exists(run_dir)) stop("run_dir does not exist: ", run_dir)
  fo2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "FixO2_invivo_attractors.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  objective_source <- o2ipa_as_chr(args$objective_source, "auto")
  o2_grid <- fixo2_attractor_o2_grid(args)
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, o2_grid)
  mode_o2_grid <- o2_grid
  max_seeds <- o2ipa_as_int(args$max_seeds, NA_integer_)
  generate_figures <- o2ipa_as_bool(args$generate_figures, TRUE)
  random_seed <- o2ipa_as_int(args$random_seed, 20260623L)
  set.seed(random_seed)

  run_args <- data.frame(
    argument = c("run_dir", "out_dir", "optional_label_file", "mode_source", "mode_rule", "mode_reference_o2", "objective_source", "attractor_o2_grid", "mode_o2_grid", "max_seeds", "generate_figures", "random_seed"),
    value = c(
      run_dir, out_dir, label_file,
      "fixed_o2_attractor_dominant_ploidy",
      "dominant_mean_ploidy >= 2 => mode1; dominant_mean_ploidy < 2 => mode2",
      mode_reference_o2, objective_source, paste(o2_grid, collapse = ","), paste(mode_o2_grid, collapse = ","),
      max_seeds, generate_figures, random_seed
    ),
    stringsAsFactors = FALSE
  )
  fo2_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  message("optional_label_file: ", if (nzchar(label_file)) label_file else "<none>")
  message("mode_reference_o2: ", mode_reference_o2)

  message("Collecting seed inputs.")
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = objective_source)
  if (!is.na(max_seeds) && max_seeds > 0L) {
    valid <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
    valid <- valid[order(o2ipa_seed_number(valid$seed_id)), , drop = FALSE]
    keep <- valid$seed_id[seq_len(min(max_seeds, nrow(valid)))]
    inputs$manifest <- inputs$manifest[inputs$manifest$seed_id %in% keep, , drop = FALSE]
    inputs$params_long <- inputs$params_long[inputs$params_long$seed_id %in% keep, , drop = FALSE]
    message("Using max_seeds subset: ", paste(keep, collapse = ", "))
  }
  extras <- fo2_seed_manifest_extras(inputs$manifest)
  manifest <- merge(inputs$manifest, extras, by = "seed_id", all.x = TRUE, sort = FALSE)
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  optional_metadata <- fo2_load_optional_seed_metadata(label_file, manifest)
  if (nrow(optional_metadata)) {
    manifest <- merge(manifest, optional_metadata, by = "seed_id", all.x = TRUE, sort = FALSE)
  }
  fo2_write_tsv(manifest, file.path(out_dir, "tables", "fixo2_seed_metadata.tsv"))
  fo2_write_tsv(manifest, file.path(out_dir, "tables", "seed_manifest_with_labels.tsv"))
  fo2_write_tsv(inputs$params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))

  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  cfg_meta <- o2pr_cfg_metadata(cfg)
  fo2_write_tsv(cfg_meta, file.path(out_dir, "tables", "fit_config_schema_summary.tsv"))
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  model_env <- o2ipa_source_model(SCRIPT_DIR)

  message("Computing fixed-O2 attractors.")
  rows <- list()
  valid_manifest <- manifest[manifest$fit_success %in% TRUE & manifest$seed_id %in% rownames(param_mat), , drop = FALSE]
  valid_manifest <- valid_manifest[order(o2ipa_seed_number(valid_manifest$seed_id)), , drop = FALSE]
  counter <- 0L
  for (i in seq_len(nrow(valid_manifest))) {
    seed <- valid_manifest$seed_id[[i]]
    if (i %% 25L == 0L || i == 1L) message("Processing ", seed, " (", i, "/", nrow(valid_manifest), ")")
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    for (O2 in mode_o2_grid) {
      counter <- counter + 1L
      rows[[counter]] <- fo2_dominant_attractor_one(seed, run_params, model_env, cfg, O2)
      rows[[counter]]$in_attractor_o2_grid <- any(abs(O2 - o2_grid) < 1e-9)
      rows[[counter]]$is_mode_reference_o2 <- abs(O2 - mode_reference_o2) < 1e-9
    }
  }
  all_mode_attractors <- do.call(rbind, rows)
  all_mode_attractors <- merge(all_mode_attractors, manifest[, intersect(c(
    "seed_id", "objective", "delta_objective", "mean_terminal_ploidy",
    "mean_late_drop_amplitude", "source_trajectory_regime",
    "source_mode_label", "source_mode_reason"
  ), names(manifest)), drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  all_mode_attractors <- fixo2_assign_attractor_modes(all_mode_attractors, "dominant_mean_ploidy")
  mode_by_seed_o2 <- fixo2_attractor_mode_table(all_mode_attractors)
  mode_reference_by_seed <- fixo2_reference_mode_table(mode_by_seed_o2, mode_reference_o2)
  fo2_write_tsv(mode_by_seed_o2, file.path(out_dir, "tables", "fixed_o2_attractor_mode_by_seed_o2.tsv"))
  fo2_write_tsv(mode_reference_by_seed, file.path(out_dir, "tables", "fixed_o2_attractor_mode_reference_by_seed.tsv"))
  fo2_write_tsv(mode_reference_by_seed, file.path(out_dir, "tables", "fixed_o2_attractor_mode_by_seed.tsv"))
  fo2_write_tsv(fixo2_attractor_mode_summary_by_seed(mode_by_seed_o2), file.path(out_dir, "tables", "fixed_o2_attractor_mode_summary_by_seed.tsv"))

  attractors <- all_mode_attractors[all_mode_attractors$in_attractor_o2_grid %in% TRUE, , drop = FALSE]
  attractors <- fixo2_assign_attractor_modes(attractors, "dominant_mean_ploidy")
  fo2_write_tsv(attractors, file.path(out_dir, "tables", "fixed_o2_attractors_by_seed.tsv"))
  fo2_write_tsv(fo2_mode_seed_stack_table(attractors), file.path(out_dir, "tables", "fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv"))
  fo2_write_mode_comparison_tables(attractors, out_dir)
  fo2_write_spectral_gap_outputs(attractors, out_dir, generate_figures = generate_figures)

  regime_summary <- do.call(rbind, lapply(split(attractors, list(attractors$O2_pct, attractors$trajectory_regime), drop = TRUE), function(d) {
    data.frame(
      O2_pct = d$O2_pct[[1]],
      trajectory_regime = d$trajectory_regime[[1]],
      n_seed = length(unique(d$seed_id)),
      dominant_mean_ploidy_median = stats::median(d$dominant_mean_ploidy, na.rm = TRUE),
      dominant_mean_ploidy_q25 = as.numeric(stats::quantile(d$dominant_mean_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
      dominant_mean_ploidy_q75 = as.numeric(stats::quantile(d$dominant_mean_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
      fraction_N_le_25_median = stats::median(d$dominant_fraction_N_le_25, na.rm = TRUE),
      fraction_N_below_44_median = stats::median(d$dominant_fraction_N_below_44, na.rm = TRUE),
      selection_22_vs_44_median = stats::median(d$selection_22_vs_44, na.rm = TRUE),
      selection_44_vs_88_median = stats::median(d$selection_44_vs_88, na.rm = TRUE),
      spectral_gap_median = stats::median(d$spectral_gap, na.rm = TRUE),
      eigenvector_nonnegative_fraction = mean(d$eigenvector_nonnegative %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  regime_summary <- regime_summary[order(regime_summary$O2_pct, regime_summary$trajectory_regime), , drop = FALSE]
  fo2_write_tsv(regime_summary, file.path(out_dir, "tables", "fixed_o2_attractor_regime_summary.tsv"))

  tests <- fo2_attractor_regime_tests(attractors)
  fo2_write_tsv(tests, file.path(out_dir, "tables", "fixed_o2_attractor_regime_tests.tsv"))

  corr <- fo2_parameter_correlations(attractors, inputs$params_long)
  fo2_write_tsv(corr[order(corr$O2_pct, corr$metric, -corr$abs_spearman_rho), , drop = FALSE], file.path(out_dir, "tables", "parameter_attractor_correlations.tsv"))

  if (generate_figures) {
    fo2_plot_outputs(fo2_mode_seed_stack_table(attractors), out_dir)
  }
  fo2_write_report(out_dir, run_dir, label_file, attractors, regime_summary, tests, corr)
  message("Completed fixed-O2 attractor analysis: ", out_dir)
  invisible(list(out_dir = out_dir, o2_grid = o2_grid, n_seed = nrow(valid_manifest)))
}


fixo2_run_counterfactual <- function(args, repo_root, run_dir, label_file, out_dir, o2_grid) {
  cf2_mkdirs(out_dir)
  time_grid <- sort(unique(cf2_as_num_vec(args$time_grid, cf2_default_time_grid())))
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  max_seeds <- o2ipa_as_int(args$max_seeds, 0L)
  generate_figures <- o2ipa_as_bool(args$generate_figures, TRUE)
  max_time <- max(time_grid, na.rm = TRUE)

  log_file <- file.path(out_dir, "logs", "FixO2_invivo_counterfactual_trajectories.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
  message("run_dir: ", run_dir)
  message("optional_label_file ignored for FixO2 mode assignment: ", if (nzchar(label_file)) label_file else "<none>")
  message("out_dir: ", out_dir)
  message("O2 grid: ", paste(o2_grid, collapse = ","))
  message("mode_reference_o2: ", mode_reference_o2)
  message("time grid length: ", length(time_grid), "; max_time=", max_time)

  seed_inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  manifest <- seed_inputs$manifest
  params_long <- seed_inputs$params_long
  if (max_seeds > 0L && nrow(manifest) > max_seeds) {
    manifest <- manifest[seq_len(max_seeds), , drop = FALSE]
    params_long <- params_long[params_long$seed_id %in% manifest$seed_id, , drop = FALSE]
  }
  cfg <- o2pr_first_seed_cfg(manifest)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  param_mat <- o2ipa_params_wide(params_long, "value")

  traj_rows <- list()
  summary_rows <- list()
  matrix_rows <- list()
  mode_reference_rows <- list()
  counter_traj <- 0L
  counter_summary <- 0L
  counter_matrix <- 0L
  counter_mode_reference <- 0L
  init_specs <- data.frame(
    initial_condition = c("init_2N", "init_4N"),
    requested_N = c(44, 88),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(manifest))) {
    seed <- manifest$seed_id[[i]]
    if (i == 1L || i %% 25L == 0L) message("Processing ", seed, " (", i, "/", nrow(manifest), ")")
    if (!seed %in% rownames(param_mat)) next
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    ref_fm <- tryCatch(cf2_fixed_matrix(model_env, cfg, run_params, mode_reference_o2), error = function(e) NULL)
    if (is.null(ref_fm)) {
      ref_dom <- data.frame(
        status = "matrix_error",
        dominant_growth_rate = NA_real_,
        second_growth_rate = NA_real_,
        spectral_gap = NA_real_,
        dominant_mean_N = NA_real_,
        dominant_mean_ploidy = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      ref_dom <- cf2_dominant_one(ref_fm$M, ref_fm$ngrid, n_unit)
    }
    mode_fields <- fixo2_mode_fields(ref_dom$dominant_mean_ploidy[[1]])
    mode_fields$mode_source <- "fixed_o2_attractor_dominant_ploidy_at_reference_o2"
    mode_fields$mode_rule <- paste0(
      "dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " >= 2 => mode1; dominant_mean_ploidy at fixed O2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " < 2 => mode2"
    )
    mode_fields$mode_reference_o2_pct <- mode_reference_o2
    mode_fields$mode_reference_o2_key <- fixo2_o2_key(mode_reference_o2)
    mode_fields$mode_reference_dominant_mean_ploidy <- ref_dom$dominant_mean_ploidy[[1]]
    mode_fields$mode_reference_status <- ref_dom$status[[1]]
    mode_fields$mode_reference_dominant_growth_rate <- ref_dom$dominant_growth_rate[[1]]
    mode_fields$mode_reference_spectral_gap <- ref_dom$spectral_gap[[1]]
    counter_mode_reference <- counter_mode_reference + 1L
    mode_reference_rows[[counter_mode_reference]] <- data.frame(seed_id = seed, mode_fields, stringsAsFactors = FALSE)
    for (O2 in o2_grid) {
      fm <- tryCatch(cf2_fixed_matrix(model_env, cfg, run_params, O2), error = function(e) NULL)
      if (is.null(fm)) next
      dom <- cf2_dominant_one(fm$M, fm$ngrid, n_unit)
      counter_matrix <- counter_matrix + 1L
      matrix_rows[[counter_matrix]] <- data.frame(seed_id = seed, O2_pct = O2, dom, mode_fields, stringsAsFactors = FALSE)
      for (j in seq_len(nrow(init_specs))) {
        init <- cf2_init_vector(fm$ngrid, init_specs$requested_N[[j]])
        sim <- cf2_eigen_trajectory(fm$M, fm$ngrid, init$vector, time_grid, n_unit)
        tr <- sim$trajectory
        if (nrow(tr)) {
          tr$seed_id <- seed
          tr$O2_pct <- O2
          tr$initial_condition <- init_specs$initial_condition[[j]]
          tr$requested_initial_N <- init_specs$requested_N[[j]]
          tr$used_initial_N <- init$used_N
          tr$status <- sim$status
          tr$trajectory_regime <- mode_fields$trajectory_regime[[1]]
          tr$mode_label <- mode_fields$mode_label[[1]]
          tr$mode_source <- mode_fields$mode_source[[1]]
          tr$mode_rule <- mode_fields$mode_rule[[1]]
          tr$mode_threshold_dominant_ploidy <- mode_fields$mode_threshold_dominant_ploidy[[1]]
          tr$mode_reference_o2_pct <- mode_fields$mode_reference_o2_pct[[1]]
          tr$mode_reference_o2_key <- mode_fields$mode_reference_o2_key[[1]]
          tr$mode_reference_dominant_mean_ploidy <- mode_fields$mode_reference_dominant_mean_ploidy[[1]]
          counter_traj <- counter_traj + 1L
          traj_rows[[counter_traj]] <- tr[, c(
            "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
            "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
            "mode_reference_dominant_mean_ploidy", "O2_pct", "initial_condition",
            "requested_initial_N", "used_initial_N", "status", "day",
            "mean_N", "mean_ploidy", "fraction_N_le_25", "fraction_N_below_44",
            "fraction_N_ge_44", "fraction_N_ge_66", "fraction_N_ge_88"
          )]
        }
        sm <- cf2_summarize_trajectory(tr, max_time)
        counter_summary <- counter_summary + 1L
        summary_rows[[counter_summary]] <- data.frame(
          seed_id = seed,
          mode_fields,
          O2_pct = O2,
          initial_condition = init_specs$initial_condition[[j]],
          requested_initial_N = init_specs$requested_N[[j]],
          used_initial_N = init$used_N,
          status = sim$status,
          sm,
          dom,
          terminal_minus_dominant_ploidy = sm$terminal_mean_ploidy - dom$dominant_mean_ploidy,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  trajectory <- if (length(traj_rows)) do.call(rbind, traj_rows) else data.frame()
  summary_by_seed <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
  dominant_by_seed <- if (length(matrix_rows)) do.call(rbind, matrix_rows) else data.frame()
  mode_reference_by_seed <- if (length(mode_reference_rows)) do.call(rbind, mode_reference_rows) else data.frame()
  counterfactual_mode_by_seed_o2 <- if (nrow(dominant_by_seed)) {
    dominant_by_seed$O2_key <- fixo2_o2_key(dominant_by_seed$O2_pct)
    dominant_by_seed[, intersect(c(
      "seed_id", "O2_pct", "O2_key", "dominant_mean_ploidy", "trajectory_regime",
      "mode_label", "mode_source", "mode_rule", "mode_threshold_dominant_ploidy",
      "mode_reference_o2_pct", "mode_reference_o2_key", "mode_reference_dominant_mean_ploidy",
      "dominant_growth_rate", "spectral_gap"
    ), names(dominant_by_seed)), drop = FALSE]
  } else {
    data.frame()
  }
  regime_summary <- cf2_regime_summary(summary_by_seed)
  tests <- cf2_regime_tests(summary_by_seed)
  correlations <- cf2_parameter_correlations(summary_by_seed, params_long)

  cf2_write_tsv(data.frame(
    argument = c("run_dir", "optional_label_file", "mode_source", "mode_rule", "mode_reference_o2", "out_dir", "o2_grid", "time_grid", "max_seeds", "generate_figures"),
    value = c(
      run_dir, label_file,
      "fixed_o2_attractor_dominant_ploidy_at_reference_o2",
      paste0(
        "dominant_mean_ploidy at fixed O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        " >= 2 => mode1; dominant_mean_ploidy at fixed O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        " < 2 => mode2"
      ),
      mode_reference_o2,
      out_dir, paste(o2_grid, collapse = ","), paste(time_grid, collapse = ","),
      as.character(max_seeds), as.character(generate_figures)
    ),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  cf2_write_tsv(trajectory, file.path(out_dir, "tables", "fixed_o2_counterfactual_trajectories.tsv"))
  cf2_write_tsv(summary_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_summary_by_seed.tsv"))
  cf2_write_tsv(dominant_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_dominant_by_seed.tsv"))
  cf2_write_tsv(counterfactual_mode_by_seed_o2, file.path(out_dir, "tables", "fixed_o2_counterfactual_mode_by_seed_o2.tsv"))
  cf2_write_tsv(mode_reference_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_mode_reference_by_seed.tsv"))
  cf2_write_tsv(regime_summary, file.path(out_dir, "tables", "fixed_o2_counterfactual_regime_summary.tsv"))
  cf2_write_tsv(tests, file.path(out_dir, "tables", "fixed_o2_counterfactual_regime_tests.tsv"))
  cf2_write_tsv(correlations, file.path(out_dir, "tables", "fixed_o2_counterfactual_parameter_correlations.tsv"))
  cf2_write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  if (isTRUE(generate_figures)) cf2_plot(trajectory, summary_by_seed, file.path(out_dir, "figures"))
  cf2_report(out_dir, args, regime_summary, tests, correlations)
  message("Completed fixed-O2 counterfactual trajectory analysis: ", out_dir)
  invisible(list(out_dir = out_dir, o2_grid = o2_grid, n_seed = nrow(manifest)))
}


fixo2_run_simulation_validation <- function(args, repo_root, run_dir, out_dir, simulation_dir, o2_grid) {
  simulation_mode <- o2ipa_as_chr(args$simulation_mode, "invivo")
  simulation_ids <- vf2_as_int_vec(args$simulation_ids, seq_len(10L))
  include_simulation <- o2ipa_as_bool(args$include_simulation, TRUE)
  generate_missing <- o2ipa_as_bool(args$generate_missing_simulation, TRUE)
  simulation_n_core <- fixo2_simulation_n_core(args)
  simulation_worker_threads <- fixo2_simulation_worker_threads(args)
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  seed_selection_n_workers <- o2ipa_as_int(args$simulation_seed_selection_n_workers %||% args$seed_selection_n_workers, simulation_n_core)
  if (!is.finite(seed_selection_n_workers) || seed_selection_n_workers < 1L) seed_selection_n_workers <- 1L
  seed_selection_n_workers <- as.integer(seed_selection_n_workers)
  manual_seed_selection <- fixo2_arg_nonempty(args$seed_ids) || fixo2_arg_nonempty(args$seed_id)
  simulation_analysis_dir <- normalizePath(file.path(out_dir, ".."), mustWork = FALSE)
  seed_mode_objective_metadata <- fixo2_simulation_mode_objective_metadata(
    analysis_dir = simulation_analysis_dir,
    run_dir = run_dir,
    o2_grid = o2_grid,
    mode_reference_o2 = mode_reference_o2,
    seed_ids = NULL,
    n_workers = seed_selection_n_workers
  )
  if (manual_seed_selection) {
    seed_ids <- vf2_as_seed_vec(args$seed_ids %||% args$seed_id, character())
    if (!length(seed_ids)) stop("Manual simulation seed selection was requested, but no valid seed_ids were supplied.")
    seed_mode_arg <- args$seed_modes %||% args$seed_labels
    seed_modes <- vf2_seed_modes(seed_ids, seed_mode_arg)
    seed_selection <- fixo2_manual_simulation_seed_table(
      seed_ids = seed_ids,
      seed_modes = seed_modes,
      metadata = seed_mode_objective_metadata,
      mode_reference_o2 = mode_reference_o2
    )
  } else {
    seed_selection <- fixo2_select_best_objective_seed_by_mode(
      metadata = seed_mode_objective_metadata,
      mode_reference_o2 = mode_reference_o2
    )
    seed_ids <- seed_selection$seed_id
  }
  seed_modes <- seed_selection$mode_label
  seed_modes[is.na(seed_modes)] <- ""
  seed_labels <- seed_ids
  has_mode <- nzchar(seed_modes)
  seed_labels[has_mode] <- paste0(seed_ids[has_mode], " (", seed_modes[has_mode], ")")
  seed_mode_map <- stats::setNames(seed_modes, seed_ids)
  seed_label_map <- stats::setNames(seed_labels, seed_ids)
  seed_selection_mode <- if (nrow(seed_selection)) seed_selection$selection_mode[[1]] else NA_character_
  default_time_grid <- if (include_simulation) {
    seq(0, o2ipa_as_num(args$time_days, 1000), by = o2ipa_as_num(args$time_step_days, 1))
  } else {
    vf2_default_time_grid()
  }
  time_grid <- sort(unique(vf2_as_num_vec(args$time_grid, default_time_grid)))
  dt_grid <- sort(unique(vf2_as_num_vec(args$dt_grid %||% args$dt, c(0.05))))
  plot_dt <- o2ipa_as_num(args$plot_dt, min(dt_grid))

  vf2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "FixO2_invivo_simulation_validation.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  message("simulation_dir: ", simulation_dir)
  message("include_simulation: ", include_simulation)
  message("generate_missing_simulation: ", generate_missing)
  message("simulation_n_core: ", simulation_n_core)
  message("simulation_worker_threads: ", simulation_worker_threads)
  message("mode_reference_o2: ", mode_reference_o2)
  message("seed_selection_mode: ", seed_selection_mode)
  message("seed_selection_n_workers: ", seed_selection_n_workers)
  if (include_simulation) message("simulation_ids: ", paste(simulation_ids, collapse = ","))
  message("seed_ids: ", paste(seed_ids, collapse = ","))
  message("seed_modes: ", paste(seed_modes, collapse = ","))
  message("seed_labels: ", paste(seed_labels, collapse = ","))
  message("O2 grid: ", paste(o2_grid, collapse = ","))
  message("dt grid: ", paste(dt_grid, collapse = ","))

  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  init_specs <- data.frame(
    initial_condition = c("init_2N", "init_4N"),
    requested_N = c(44, 88),
    initial_ploidy = c(2, 4),
    stringsAsFactors = FALSE
  )
  vf2_write_tsv(seed_selection, file.path(out_dir, "tables", "fixed_o2_simulation_representative_seeds.tsv"))

  if (include_simulation) {
    fixo2_ensure_simulation(
      args = args,
      run_dir = run_dir,
      simulation_dir = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      o2_grid = o2_grid,
      init_specs = init_specs,
      simulation_ids = simulation_ids,
      out_dir = out_dir
    )
  }

  trajectory_rows <- list()
  state_error_rows <- list()
  summary_rows <- list()
  tc <- 0L
  ec <- 0L
  sc <- 0L
  for (seed_id in seed_ids) {
    message("Processing seed=", seed_id, " label=", seed_label_map[[seed_id]])
    manifest <- inputs$manifest[inputs$manifest$seed_id == seed_id, , drop = FALSE]
    if (!nrow(manifest)) stop("seed_id not found in run_dir: ", seed_id)
    params_long <- inputs$params_long[inputs$params_long$seed_id == seed_id, , drop = FALSE]
    if (!nrow(params_long)) stop("No parameters found for ", seed_id)
    cfg <- o2pr_first_seed_cfg(manifest)
    n_unit <- as.numeric(cfg$N_UNIT %||% 22)
    if (!seed_id %in% rownames(param_mat)) stop("seed_id missing from parameter matrix: ", seed_id)
    pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    for (O2 in o2_grid) {
      message("  O2=", O2)
      fm <- vf2_fixed_matrix(model_env, cfg, run_params, O2)
      for (j in seq_len(nrow(init_specs))) {
        init <- vf2_init_vector(fm$ngrid, init_specs$requested_N[[j]])
        eig <- vf2_eigen_states(fm$M, init$vector, time_grid)
        expm <- vf2_expm_states(fm$M, init$vector, time_grid)
        for (dt in dt_grid) {
          message("    ", init_specs$initial_condition[[j]], ", dt=", dt)
          eu <- vf2_euler_states(fm$M, init$vector, time_grid, dt)
          comp <- vf2_compare_states(eig$states, expm, eu, fm$ngrid, n_unit, time_grid)
          comp$seed_id <- seed_id
          comp$seed_mode <- seed_mode_map[[seed_id]]
          comp$seed_label <- seed_label_map[[seed_id]]
          comp$O2_pct <- O2
          comp$initial_condition <- init_specs$initial_condition[[j]]
          comp$requested_initial_N <- init_specs$requested_N[[j]]
          comp$initial_ploidy <- init_specs$initial_ploidy[[j]]
          comp$used_initial_N <- init$used_N
          comp$dt <- dt
          comp$dominant_growth_rate <- eig$lambda_ref
          comp <- comp[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "requested_initial_N", "used_initial_N",
            "dt", "day", "dominant_growth_rate",
            "eigen_mean_N", "expm_mean_N", "euler_mean_N",
            "diff_euler_expm_mean_N", "abs_diff_euler_expm_mean_N",
            "diff_eigen_expm_mean_N", "abs_diff_eigen_expm_mean_N",
            "eigen_mean_ploidy", "expm_mean_ploidy", "euler_mean_ploidy",
            "diff_euler_expm_mean_ploidy", "abs_diff_euler_expm_mean_ploidy",
            "diff_eigen_expm_mean_ploidy", "abs_diff_eigen_expm_mean_ploidy",
            "eigen_sd_ploidy", "expm_sd_ploidy", "euler_sd_ploidy",
            "diff_euler_expm_sd_ploidy", "abs_diff_euler_expm_sd_ploidy",
            "diff_eigen_expm_sd_ploidy", "abs_diff_eigen_expm_sd_ploidy",
            "eigen_fraction_N_le_25", "expm_fraction_N_le_25", "euler_fraction_N_le_25",
            "eigen_fraction_N_below_44", "expm_fraction_N_below_44", "euler_fraction_N_below_44",
            "eigen_fraction_N_ge_44", "expm_fraction_N_ge_44", "euler_fraction_N_ge_44",
            "eigen_fraction_N_ge_66", "expm_fraction_N_ge_66", "euler_fraction_N_ge_66",
            "eigen_fraction_N_ge_88", "expm_fraction_N_ge_88", "euler_fraction_N_ge_88",
            "max_abs_state_diff_euler_expm", "l1_state_diff_euler_expm",
            "max_abs_state_diff_eigen_expm", "l1_state_diff_eigen_expm"
          )]
          tc <- tc + 1L
          trajectory_rows[[tc]] <- comp
          serr <- comp[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "dt", "day",
            "max_abs_state_diff_euler_expm", "l1_state_diff_euler_expm",
            "max_abs_state_diff_eigen_expm", "l1_state_diff_eigen_expm"
          ), drop = FALSE]
          ec <- ec + 1L
          state_error_rows[[ec]] <- serr
          sm <- vf2_error_summary_one(comp)
          sm$seed_id <- seed_id
          sm$seed_mode <- seed_mode_map[[seed_id]]
          sm$seed_label <- seed_label_map[[seed_id]]
          sm$O2_pct <- O2
          sm$initial_condition <- init_specs$initial_condition[[j]]
          sm$initial_ploidy <- init_specs$initial_ploidy[[j]]
          sm$requested_initial_N <- init_specs$requested_N[[j]]
          sm$used_initial_N <- init$used_N
          sm$dt <- dt
          sm$dominant_growth_rate <- eig$lambda_ref
          sm <- sm[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "requested_initial_N", "used_initial_N",
            "dt", "dominant_growth_rate", "n_timepoint",
            "max_abs_diff_euler_expm_mean_ploidy", "terminal_abs_diff_euler_expm_mean_ploidy",
            "max_abs_diff_eigen_expm_mean_ploidy", "terminal_abs_diff_eigen_expm_mean_ploidy",
            "max_abs_diff_euler_expm_sd_ploidy", "terminal_abs_diff_euler_expm_sd_ploidy",
            "max_abs_diff_eigen_expm_sd_ploidy", "terminal_abs_diff_eigen_expm_sd_ploidy",
            "max_abs_diff_euler_expm_mean_N", "terminal_abs_diff_euler_expm_mean_N",
            "max_abs_state_diff_euler_expm", "terminal_max_abs_state_diff_euler_expm",
            "max_l1_state_diff_euler_expm", "terminal_l1_state_diff_euler_expm"
          )]
          sc <- sc + 1L
          summary_rows[[sc]] <- sm
        }
      }
    }
  }

  trajectories <- do.call(rbind, trajectory_rows)
  state_errors <- do.call(rbind, state_error_rows)
  summaries <- do.call(rbind, summary_rows)
  summaries <- summaries[order(match(summaries$seed_id, seed_ids), summaries$O2_pct, summaries$initial_condition, summaries$dt), , drop = FALSE]
  seed_slug <- paste(seed_ids, collapse = "_")

  simulation_traj <- data.frame()
  simulation_summary <- data.frame()
  expm_simulation_comparison <- data.frame()
  expm_simulation_comparison_summary <- data.frame()
  eigen_simulation_comparison <- data.frame()
  eigen_simulation_comparison_summary <- data.frame()
  expm_phase_plane <- data.frame()
  eigen_phase_plane <- data.frame()
  simulation_phase_plane <- data.frame()
  if (include_simulation) {
    message("Reading fixed-O2 simulation trajectories")
    simulation_traj <- vf2_read_simulation_trajectories(
      sim_root = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      seed_label_map = seed_label_map,
      seed_mode_map = seed_mode_map,
      o2_grid = o2_grid,
      init_specs = init_specs,
      simulation_ids = simulation_ids
    )
    simulation_summary <- vf2_simulation_summary(simulation_traj)
    expm_simulation_comparison <- vf2_compare_solution_to_simulation(trajectories, simulation_summary, "expm")
    expm_simulation_comparison_summary <- vf2_simulation_comparison_summary(expm_simulation_comparison)
    eigen_simulation_comparison <- vf2_compare_solution_to_simulation(trajectories, simulation_summary, "eigen")
    eigen_simulation_comparison_summary <- vf2_simulation_comparison_summary(eigen_simulation_comparison)
  }
  expm_phase_plane <- vf2_solution_phase_plane(trajectories, plot_dt, "expm", seed_mode_map)
  eigen_phase_plane <- vf2_solution_phase_plane(trajectories, plot_dt, "eigen", seed_mode_map)
  if (include_simulation) {
    simulation_phase_plane <- vf2_simulation_phase_plane(simulation_traj, seed_mode_map)
  }

  run_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "simulation_dir", "simulation_mode", "include_simulation",
      "generate_missing_simulation", "simulation_n_core", "simulation_worker_threads",
      "mode_reference_o2", "seed_selection_mode", "seed_selection_n_workers",
      "simulation_ids", "seed_ids", "seed_modes", "seed_labels", "o2_grid", "time_grid", "dt_grid", "plot_dt"
    ),
    value = c(
      run_dir, out_dir, simulation_dir, simulation_mode, as.character(include_simulation),
      as.character(generate_missing), as.character(simulation_n_core), as.character(simulation_worker_threads),
      as.character(mode_reference_o2), seed_selection_mode, as.character(seed_selection_n_workers),
      paste(simulation_ids, collapse = ","), paste(seed_ids, collapse = ","), paste(seed_modes, collapse = ","), paste(seed_labels, collapse = ","),
      paste(o2_grid, collapse = ","), paste(time_grid, collapse = ","), paste(dt_grid, collapse = ","), as.character(plot_dt)
    ),
    stringsAsFactors = FALSE
  )
  vf2_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  vf2_write_tsv(trajectories, file.path(out_dir, "tables", "eigen_vs_euler_trajectories.tsv"))
  vf2_write_tsv(state_errors, file.path(out_dir, "tables", "eigen_vs_euler_state_errors.tsv"))
  vf2_write_tsv(summaries, file.path(out_dir, "tables", "eigen_vs_euler_error_summary.tsv"))
  vf2_write_tsv(expm_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_expm_analytical.tsv"))
  vf2_write_tsv(eigen_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_eigen_analytical.tsv"))
  if (include_simulation) {
    vf2_write_tsv(simulation_traj, file.path(out_dir, "tables", "simulation_replicate_mean_ploidy_trajectories.tsv"))
    vf2_write_tsv(simulation_summary, file.path(out_dir, "tables", "simulation_summary_mean_ploidy_trajectories.tsv"))
    vf2_write_tsv(expm_simulation_comparison, file.path(out_dir, "tables", "expm_vs_simulation_mean_ploidy_comparison.tsv"))
    vf2_write_tsv(expm_simulation_comparison_summary, file.path(out_dir, "tables", "expm_vs_simulation_error_summary.tsv"))
    vf2_write_tsv(eigen_simulation_comparison, file.path(out_dir, "tables", "eigen_vs_simulation_mean_ploidy_comparison.tsv"))
    vf2_write_tsv(eigen_simulation_comparison_summary, file.path(out_dir, "tables", "eigen_vs_simulation_error_summary.tsv"))
    vf2_write_tsv(simulation_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_simulation_replicates.tsv"))
    vf2_write_tsv(rbind(expm_phase_plane, simulation_phase_plane), file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_expm_vs_simulation.tsv"))
    vf2_write_tsv(rbind(eigen_phase_plane, simulation_phase_plane), file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_eigen_vs_simulation.tsv"))
  }
  vf2_plot(trajectories, file.path(out_dir, "figures", sprintf("expm_vs_euler_mean_ploidy_%s.pdf", seed_slug)), plot_dt = plot_dt)
  if (include_simulation) {
    vf2_plot_solution_vs_simulation(
      trajectories,
      simulation_traj,
      file.path(out_dir, "figures", sprintf("expm_vs_simulation_replicates_mean_ploidy_%s.pdf", seed_slug)),
      plot_dt = plot_dt,
      solution_name = "expm",
      solution_label = "Expm analytical",
      title = "Fixed-O2 expm analytical vs simulation replicates"
    )
    vf2_plot_solution_vs_simulation(
      trajectories,
      simulation_traj,
      file.path(out_dir, "figures", sprintf("eigen_vs_simulation_replicates_mean_ploidy_%s.pdf", seed_slug)),
      plot_dt = plot_dt,
      solution_name = "eigen",
      solution_label = "Eigen analytical",
      title = "Fixed-O2 eigen analytical vs simulation replicates"
    )
    vf2_plot_phase_plane_solution_vs_simulation(
      expm_phase_plane,
      simulation_phase_plane,
      file.path(out_dir, "figures", sprintf("expm_vs_simulation_phase_plane_mean_sd_ploidy_%s.pdf", seed_slug)),
      solution_label = "Expm analytical",
      title = "Fixed-O2 phase plane: expm analytical vs simulation replicates"
    )
    vf2_plot_phase_plane_solution_vs_simulation(
      eigen_phase_plane,
      simulation_phase_plane,
      file.path(out_dir, "figures", sprintf("eigen_vs_simulation_phase_plane_mean_sd_ploidy_%s.pdf", seed_slug)),
      solution_label = "Eigen analytical",
      title = "Fixed-O2 phase plane: eigen analytical vs simulation replicates"
    )
  }
  message("Completed simulation validation: ", out_dir)
  invisible(list(out_dir = out_dir, o2_grid = o2_grid, seed_ids = seed_ids, simulation_dir = simulation_dir))
}


path_results_suffix <- function(path) {
  path <- normalizePath(as.character(path), mustWork = FALSE)
  suffix <- sub("^.*(/oxygen/results/.*)$", "\\1", path)
  ifelse(identical(suffix, path), path, suffix)
}


split_csv <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) {
    return(default)
  }
  vals <- trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}


normalize_seed_ids <- function(x) {
  if (is.null(x) || !length(x)) return(character())
  o2ipa_norm_seed(x)
}


as_num_vec <- function(x, default) {
  vals <- suppressWarnings(as.numeric(split_csv(x, as.character(default))))
  vals <- vals[is.finite(vals)]
  if (length(vals)) unique(vals) else default
}


as_int_vec <- function(x, default) {
  vals <- suppressWarnings(as.integer(split_csv(x, as.character(default))))
  vals <- vals[!is.na(vals)]
  if (length(vals)) unique(vals) else default
}


as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  tolower(as.character(x[[1]])) %in% c("true", "1", "yes", "y", "on")
}


num_key <- function(x) {
  vapply(as.numeric(x), function(val) {
    if (!is.finite(val)) return(NA_character_)
    format(signif(val, 12), scientific = FALSE, trim = TRUE)
  }, character(1))
}




read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing file: ", path)
    return(data.frame())
  }
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (is_gz) {
      gzip <- Sys.which("gzip")
      if (!nzchar(gzip)) stop("Reading gzip-compressed tables requires gzip on PATH: ", path)
      return(as.data.frame(data.table::fread(cmd = paste(shQuote(gzip), "-cd", shQuote(path)), sep = "\t", data.table = FALSE, showProgress = FALSE)))
    }
    return(as.data.frame(data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE)))
  }
  if (is_gz) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(utils::read.table(
      con,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ))
  }
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}


write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    wrote <- FALSE
    if (requireNamespace("data.table", quietly = TRUE)) {
      wrote <- tryCatch({
        data.table::fwrite(x, file = path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
        TRUE
      }, error = function(e) FALSE)
    }
    if (!wrote) {
      con <- gzfile(path, open = "wt")
      on.exit(close(con), add = TRUE)
      utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    }
  } else {
    utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  invisible(path)
}


filter_by_values <- function(df, col, values) {
  if (!nrow(df) || !(col %in% names(df))) return(df)
  key <- num_key(df[[col]])
  keep <- key %in% num_key(values)
  df[keep, , drop = FALSE]
}


fixo2_arg_value <- function(args, primary, fallback = NULL, default = NULL) {
  for (key in c(primary, fallback)) {
    if (!is.null(key) && !is.null(args[[key]]) && length(args[[key]]) && !is.na(args[[key]][[1]])) {
      return(args[[key]])
    }
  }
  default
}


fixo2_agreement_cache_read_path <- function(path, legacy_path, recompute = FALSE) {
  if (isTRUE(recompute) || file.exists(path)) return(path)
  legacy_path <- legacy_path[nzchar(legacy_path)]
  hit <- legacy_path[file.exists(legacy_path)]
  if (length(hit)) hit[[1]] else path
}


fixo2_run_analytical_simulation_agreement <- function(args, repo_root, run_dir, analysis_dir, out_dir, simulation_dir, o2_grid) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package is required for analytical-simulation agreement plots.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(out_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  log_dir <- file.path(out_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(log_dir, "FixO2_invivo_analytical_simulation_agreement.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  fit_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_fit_dir", "fit_dir", run_dir), run_dir),
    repo_root,
    mustWork = FALSE
  )
  agreement_run_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_run_dir", "run_dir", fit_dir), fit_dir),
    repo_root,
    mustWork = FALSE
  )
  agreement_analysis_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_analysis_dir", "analysis_dir", analysis_dir), analysis_dir),
    repo_root,
    mustWork = FALSE
  )
  agreement_simulation_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "agreement_simulation_dir", "simulation_dir", simulation_dir), simulation_dir),
    repo_root,
    mustWork = FALSE
  )
  simulation_mode <- o2ipa_as_chr(args$simulation_mode, "invivo")
  time_points <- sort(as_num_vec(fixo2_arg_value(args, "agreement_time_points", "time_points", "25,50,100,200,300,500,700,1000"), c(25, 50, 100, 200, 300, 500, 700, 1000)))
  agreement_o2 <- sort(as_num_vec(fixo2_arg_value(args, "agreement_o2_values", "o2_values", paste(o2_grid, collapse = ",")), o2_grid))
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  initial_ploidy_values <- sort(as_num_vec(fixo2_arg_value(args, "agreement_initial_ploidy_values", "initial_ploidy_values", "2,4"), c(2, 4)))
  simulation_ids <- sort(as_int_vec(fixo2_arg_value(args, "agreement_simulation_ids", "simulation_ids", "1,2,3"), c(1L, 2L, 3L)))
  objective_transform <- o2ipa_as_chr(fixo2_arg_value(args, "agreement_objective_transform", "objective_transform", "identity"), "identity")
  if (!objective_transform %in% c("identity", "log10")) stop("--agreement_objective_transform must be identity or log10.")
  recompute <- o2ipa_as_bool(fixo2_arg_value(args, "agreement_recompute", "recompute", FALSE), FALSE)
  cache_all_times <- o2ipa_as_bool(fixo2_arg_value(args, "agreement_cache_all_times", "cache_all_times", TRUE), TRUE)
  recompute_analytical <- o2ipa_as_bool(fixo2_arg_value(args, "agreement_recompute_analytical", "recompute_analytical", recompute), recompute)
  progress_every <- o2ipa_as_int(fixo2_arg_value(args, "agreement_progress_every", "progress_every", 100L), 100L)
  n_workers <- o2ipa_as_int(
    fixo2_arg_value(args, "agreement_n_workers", "n_workers", Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
    1L
  )
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  analytical_methods <- split_csv(fixo2_arg_value(args, "agreement_analytical_methods", "analytical_methods", "eigen,expm"), default = c("eigen", "expm"))
  analytical_methods <- unique(vapply(analytical_methods, method_slug, character(1)))
  seed_ids <- normalize_seed_ids(split_csv(fixo2_arg_value(args, "agreement_seed_ids", "seed_ids", ""), default = character()))

  legacy_table_dir <- file.path(agreement_analysis_dir, "simulation", "scatters", "tables")
  analytical_cache_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_analytical_cache_table", "analytical_cache_table", file.path(table_dir, "agreement_generated_analytical_trajectories.tsv.gz")),
      file.path(table_dir, "agreement_generated_analytical_trajectories.tsv.gz")
    ),
    repo_root,
    mustWork = FALSE
  )
  all_time_sim_metric_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_all_time_simulation_metric_table", "all_time_simulation_metric_table", file.path(table_dir, "agreement_simulation_all_time_metrics_by_replicate.tsv.gz")),
      file.path(table_dir, "agreement_simulation_all_time_metrics_by_replicate.tsv.gz")
    ),
    repo_root,
    mustWork = FALSE
  )
  sim_metric_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_simulation_metric_table", "simulation_metric_table", file.path(table_dir, "agreement_simulation_selected_time_metrics_by_replicate.tsv")),
      file.path(table_dir, "agreement_simulation_selected_time_metrics_by_replicate.tsv")
    ),
    repo_root,
    mustWork = FALSE
  )
  sim_summary_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_simulation_summary_table", "simulation_summary_table", file.path(table_dir, "agreement_simulation_selected_time_metrics.tsv")),
      file.path(table_dir, "agreement_simulation_selected_time_metrics.tsv")
    ),
    repo_root,
    mustWork = FALSE
  )
  agreement_data_path <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "agreement_data_table", "scatter_data_table", file.path(table_dir, "agreement_analytical_vs_simulation_data.tsv")),
      file.path(table_dir, "agreement_analytical_vs_simulation_data.tsv")
    ),
    repo_root,
    mustWork = FALSE
  )

  analytical_cache_read_path <- fixo2_agreement_cache_read_path(
    analytical_cache_path,
    c(
      file.path(table_dir, "scatter_generated_analytical_trajectories.tsv.gz"),
      file.path(legacy_table_dir, "scatter_generated_analytical_trajectories.tsv.gz")
    ),
    recompute_analytical
  )
  all_time_sim_metric_read_path <- fixo2_agreement_cache_read_path(
    all_time_sim_metric_path,
    c(
      file.path(table_dir, "scatter_simulation_all_time_metrics_by_replicate.tsv.gz"),
      file.path(legacy_table_dir, "scatter_simulation_all_time_metrics_by_replicate.tsv.gz")
    ),
    recompute
  )
  sim_metric_read_path <- fixo2_agreement_cache_read_path(
    sim_metric_path,
    c(
      file.path(table_dir, "scatter_simulation_selected_time_metrics_by_replicate.tsv"),
      file.path(legacy_table_dir, "scatter_simulation_selected_time_metrics_by_replicate.tsv")
    ),
    recompute
  )
  sim_summary_read_path <- fixo2_agreement_cache_read_path(
    sim_summary_path,
    c(
      file.path(table_dir, "scatter_simulation_selected_time_metrics.tsv"),
      file.path(legacy_table_dir, "scatter_simulation_selected_time_metrics.tsv")
    ),
    recompute
  )

  message("run_dir: ", agreement_run_dir)
  message("analysis_dir: ", agreement_analysis_dir)
  message("out_dir: ", out_dir)
  message("simulation_dir: ", agreement_simulation_dir)
  message("simulation_mode: ", simulation_mode)
  message("time_points: ", paste(time_points, collapse = ","))
  message("O2 values: ", paste(agreement_o2, collapse = ","))
  message("mode_reference_o2: ", mode_reference_o2)
  message("initial_ploidy_values: ", paste(initial_ploidy_values, collapse = ","))
  message("simulation_ids: ", paste(simulation_ids, collapse = ","))
  message("analytical_methods: ", paste(analytical_methods, collapse = ","))
  message("n_workers: ", n_workers)

  expected_seed_ids <- expected_analytical_seed_ids(agreement_run_dir, seed_ids)
  message("Generating or reading analytical trajectories.")
  if (!recompute_analytical && file.exists(analytical_cache_read_path)) {
    message("Reading cached analytical trajectories: ", analytical_cache_read_path)
    analytical <- read_tsv(analytical_cache_read_path)
    missing_keys <- analytical_cache_missing_keys(
      analytical = analytical,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      expected_seed_ids = expected_seed_ids,
      expected_run_dir = agreement_run_dir
    )
    analytical <- filter_analytical_trajectories(
      analytical = analytical,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      seed_ids = seed_ids
    )
    missing_methods <- setdiff(analytical_methods, unique(analytical$analytical_method))
    if (length(missing_methods) || length(missing_keys)) analytical <- data.frame()
  } else {
    analytical <- data.frame()
  }
  if (!nrow(analytical)) {
    analytical <- generate_analytical_trajectories(
      run_dir = agreement_run_dir,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      n_workers = n_workers,
      seed_ids = seed_ids
    )
    write_tsv(analytical, analytical_cache_path)
    analytical <- filter_analytical_trajectories(
      analytical = analytical,
      time_points = time_points,
      o2_values = agreement_o2,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      seed_ids = seed_ids
    )
  } else if (!identical(normalizePath(analytical_cache_read_path, mustWork = FALSE), normalizePath(analytical_cache_path, mustWork = FALSE))) {
    write_tsv(analytical, analytical_cache_path)
  }
  if (!nrow(analytical)) stop("No analytical trajectory rows were found for the requested O2/time grid.")

  message("Reading seed objective values.")
  objectives <- read_seed_objectives(
    analysis_dir = agreement_analysis_dir,
    fit_dir = fit_dir,
    o2_values = agreement_o2,
    seed_ids = seed_ids,
    n_workers = n_workers,
    mode_reference_o2 = mode_reference_o2
  )

  if (!cache_all_times && !recompute && file.exists(sim_summary_read_path)) {
    message("Reading cached simulation summary: ", sim_summary_read_path)
    sim_summary <- read_tsv(sim_summary_read_path)
    if (!identical(normalizePath(sim_summary_read_path, mustWork = FALSE), normalizePath(sim_summary_path, mustWork = FALSE))) {
      write_tsv(sim_summary, sim_summary_path)
    }
  } else {
    if (cache_all_times) {
      if (!recompute && file.exists(all_time_sim_metric_read_path)) {
        message("Reading cached all-time simulation metrics: ", all_time_sim_metric_read_path)
        all_time_sim_metrics <- read_tsv(all_time_sim_metric_read_path)
      } else {
        if (!dir.exists(agreement_simulation_dir)) stop("Simulation directory is required when all-time cache is missing: ", agreement_simulation_dir)
        tasks <- read_simulation_tasks(
          simulation_dir = agreement_simulation_dir,
          simulation_mode = simulation_mode,
          analytical = analytical,
          o2_values = agreement_o2,
          initial_ploidy_values = initial_ploidy_values,
          simulation_ids = simulation_ids
        )
        if (!nrow(tasks)) stop("No simulation tasks matched the requested seed/O2/initial ploidy/simulation id filters.")
        message("Matched simulation tasks: ", nrow(tasks))
        all_time_sim_metrics <- read_simulation_metrics(
          tasks = tasks,
          time_points = NULL,
          progress_every = progress_every,
          n_workers = n_workers
        )
        if (!nrow(all_time_sim_metrics)) stop("No all-time simulation metrics were read from state trajectories.")
        write_tsv(all_time_sim_metrics, all_time_sim_metric_path)
      }
      sim_metrics <- filter_simulation_metrics(
        sim_metrics = all_time_sim_metrics,
        time_points = time_points,
        o2_values = agreement_o2,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids,
        seed_ids = unique(analytical$seed_id)
      )
    } else if (!recompute && file.exists(sim_metric_read_path)) {
      message("Reading cached selected-time simulation metrics: ", sim_metric_read_path)
      sim_metrics <- read_tsv(sim_metric_read_path)
      sim_metrics <- filter_simulation_metrics(
        sim_metrics = sim_metrics,
        time_points = time_points,
        o2_values = agreement_o2,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids,
        seed_ids = unique(analytical$seed_id)
      )
    } else {
      if (!dir.exists(agreement_simulation_dir)) stop("Simulation directory is required when selected-time cache is missing: ", agreement_simulation_dir)
      tasks <- read_simulation_tasks(
        simulation_dir = agreement_simulation_dir,
        simulation_mode = simulation_mode,
        analytical = analytical,
        o2_values = agreement_o2,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids
      )
      if (!nrow(tasks)) stop("No simulation tasks matched the requested seed/O2/initial ploidy/simulation id filters.")
      message("Matched simulation tasks: ", nrow(tasks))
      sim_metrics <- read_simulation_metrics(
        tasks = tasks,
        time_points = time_points,
        progress_every = progress_every,
        n_workers = n_workers
      )
    }
    if (!nrow(sim_metrics)) stop("No simulation metrics were read from state trajectories.")
    write_tsv(sim_metrics, sim_metric_path)
    sim_summary <- aggregate_replicates(sim_metrics)
    write_tsv(sim_summary, sim_summary_path)
  }

  message("Merging analytical and simulation summaries.")
  agreement_rows <- list()
  fig_dirs <- character()
  plot_dirs <- character()
  for (method in sort(unique(analytical$analytical_method))) {
    method_analytical <- analytical[analytical$analytical_method %in% method, , drop = FALSE]
    plot_dirs[[method]] <- make_analytical_solution_outputs(
      analytical = method_analytical,
      out_dir = out_dir,
      analytical_method = method,
      seed_metadata = objectives
    )
    agreement_data <- merge_scatter_data(method_analytical, sim_summary, objectives)
    if (!nrow(agreement_data)) {
      warning("No merged analytical-vs-simulation rows were produced for method: ", method)
      next
    }
    method_table_dir <- file.path(table_dir, method_slug(method))
    dir.create(method_table_dir, recursive = TRUE, showWarnings = FALSE)
    method_agreement_data_path <- file.path(method_table_dir, "agreement_analytical_vs_simulation_data.tsv")
    write_tsv(agreement_data, method_agreement_data_path)
    agreement_rows[[method]] <- agreement_data
    message("Drawing analytical-simulation agreement plots for method: ", method)
    agreement_plot_data <- fixo2_add_agreement_limits(agreement_data)
    fig_dirs[[method]] <- make_scatter_outputs(
      dat = agreement_plot_data,
      out_dir = out_dir,
      objective_transform = objective_transform,
      analytical_method = method
    )
  }
  if (!length(fig_dirs)) stop("No analytical-simulation agreement figures were produced for any analytical method.")
  combined_agreement_data <- do.call(rbind, agreement_rows)
  if (!is.null(combined_agreement_data) && nrow(combined_agreement_data)) write_tsv(combined_agreement_data, agreement_data_path)

  manifest <- data.frame(
    field = c(
      "simulation_dir", "analysis_dir", "fit_dir", "run_dir", "out_dir", "simulation_mode",
      "time_points", "o2_values", "mode_reference_o2", "initial_ploidy_values", "simulation_ids",
      "objective_transform", "cache_all_times", "n_workers", "seed_ids",
      "analytical_methods", "recompute", "recompute_analytical", "analytical_cache_table",
      "all_time_simulation_metric_table", "simulation_metric_table", "simulation_summary_table",
      "agreement_data_table", "figure_dirs", "analytical_solution_plot_dirs"
    ),
    value = c(
      agreement_simulation_dir, agreement_analysis_dir, fit_dir, agreement_run_dir, out_dir, simulation_mode,
      paste(time_points, collapse = ","), paste(agreement_o2, collapse = ","), mode_reference_o2,
      paste(initial_ploidy_values, collapse = ","), paste(simulation_ids, collapse = ","),
      objective_transform, as.character(cache_all_times), as.character(n_workers), paste(seed_ids, collapse = ","),
      paste(analytical_methods, collapse = ","), as.character(recompute), as.character(recompute_analytical), analytical_cache_path,
      all_time_sim_metric_path, sim_metric_path, sim_summary_path, agreement_data_path, paste(fig_dirs, collapse = ","), paste(plot_dirs, collapse = ",")
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(table_dir, "analytical_simulation_agreement_run_manifest.tsv"))
  message("Completed analytical-simulation agreement analysis: ", out_dir)
  invisible(list(out_dir = out_dir, figure_dirs = fig_dirs, analytical_methods = analytical_methods))
}


fixo2_main <- function(args = o2ipa_parse_args()) {
  repo_root <- o2ipa_repo_root(SCRIPT_DIR)
  run_dir <- fixo2_resolve_repo_path(o2ipa_as_chr(args$run_dir, fixo2_default_run_dir(repo_root)), repo_root, mustWork = FALSE)
  out_dir <- fixo2_resolve_repo_path(o2ipa_as_chr(args$out_dir, file.path(repo_root, "oxygen", "results", "analysis", "FixO2_invivo_500seed")), repo_root, mustWork = FALSE)
  label_file <- fixo2_resolve_repo_path(o2ipa_as_chr(args$label_file, fixo2_default_label_file(repo_root)), repo_root, mustWork = FALSE)
  simulation_dir <- fixo2_resolve_repo_path(o2ipa_as_chr(args$simulation_dir, fixo2_default_simulation_dir(repo_root)), repo_root, mustWork = FALSE)
  o2_grid <- sort(unique(fixo2_as_num_vec(args$o2_grid, c(0, 0.1, 0.5, 1, 2, 5))))
  mode_reference_o2 <- fixo2_mode_reference_o2(args)
  fixo2_validate_mode_reference_o2(mode_reference_o2, fixo2_attractor_o2_grid(args))
  parts <- fixo2_normalize_parts(args$run_parts)
  fixo2_prepare_dirs(out_dir)

  render_html_report <- o2ipa_as_bool(fixo2_arg_value(args, "render_html_report", "generate_html_report", TRUE), TRUE)
  html_report_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "html_report_dir", "report_out_dir", file.path(out_dir, "html_report")), file.path(out_dir, "html_report")),
    repo_root,
    mustWork = FALSE
  )
  html_report_basename <- o2ipa_as_chr(fixo2_arg_value(args, "html_report_basename", "report_basename", "index"), "index")

  top_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "optional_label_file", "mode_source", "mode_rule", "mode_reference_o2",
      "simulation_dir", "o2_grid", "run_parts", "simulation_n_core", "simulation_worker_threads",
      "render_html_report", "html_report_dir", "html_report_basename"
    ),
    value = c(
      run_dir, out_dir, label_file,
      "fixed_o2_attractor_dominant_ploidy_at_reference_o2",
      paste0(
        "downstream analyses use the seed-level mode assigned at fixed O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        "; the complete attractor output still reports O2-dependent modes for all attractor O2 values"
      ),
      mode_reference_o2,
      simulation_dir, paste(o2_grid, collapse = ","), paste(parts, collapse = ","),
      as.character(fixo2_simulation_n_core(args)), as.character(fixo2_simulation_worker_threads(args)),
      as.character(render_html_report), html_report_dir, html_report_basename
    ),
    stringsAsFactors = FALSE
  )
  fixo2_write_tsv(top_args, file.path(out_dir, "tables", "FixO2_invivo_run_arguments.tsv"))

  results <- list()
  if ("attractors" %in% parts) {
    results$attractors <- fixo2_run_attractors(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      label_file = label_file,
      out_dir = file.path(out_dir, "attractors")
    )
  }
  if ("counterfactual_trajectories" %in% parts) {
    results$counterfactual_trajectories <- fixo2_run_counterfactual(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      label_file = label_file,
      out_dir = file.path(out_dir, "counterfactual_trajectories"),
      o2_grid = o2_grid
    )
  }
  if ("simulation" %in% parts) {
    results$simulation <- fixo2_run_simulation_validation(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      out_dir = file.path(out_dir, "simulation"),
      simulation_dir = simulation_dir,
      o2_grid = o2_grid
    )
  }
  if ("analytical_simulation_agreement" %in% parts) {
    results$analytical_simulation_agreement <- fixo2_run_analytical_simulation_agreement(
      args = args,
      repo_root = repo_root,
      run_dir = run_dir,
      analysis_dir = out_dir,
      out_dir = file.path(out_dir, "simulation", "analytical_simulation_agreement"),
      simulation_dir = simulation_dir,
      o2_grid = o2_grid
    )
  }

  html_report <- fixo2_render_html_report(args, repo_root = repo_root, out_dir = out_dir)
  if (!is.null(html_report)) results$html_report <- html_report
  fixo2_write_output_manifest(out_dir)
  message("Completed FixO2 in vivo workflow: ", out_dir)
  invisible(results)
}
