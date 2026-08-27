#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
if (identical(Sys.getenv("FIGURE5_LOCAL_ONLY"), "1")) {
  WORKSPACE_ROOT <- normalizePath(
    Sys.getenv(
      "FIGURE_WORKSPACE_ROOT",
      unset = file.path(script_dir, "..", "..")
    ),
    mustWork = TRUE
  )
  DATA_ROOT <- file.path(WORKSPACE_ROOT, "data", "Figures")
  require_files <- function(paths, label = "required file") {
    missing <- paths[!file.exists(paths)]
    if (length(missing)) {
      stop("Missing ", label, "(s):\n", paste(missing, collapse = "\n"))
    }
    invisible(normalizePath(paths, mustWork = TRUE))
  }
  write_intermediate_tsv <- function(x, path) {
    data_root <- normalizePath(DATA_ROOT, mustWork = TRUE)
    normalized <- normalizePath(path, mustWork = FALSE)
    if (!startsWith(normalized, paste0(data_root, .Platform$file.sep))) {
      stop("Local Figure 5 output is outside iteration4 data: ", normalized)
    }
    dir.create(dirname(normalized), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      x,
      normalized,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
    normalizePath(normalized, mustWork = TRUE)
  }
} else {
  source(file.path(script_dir, "util", "runtime", "process_runner.R"))
}

figure5_select_winners <- function() {
  manifest <- read_joint_warmup_manifest()
  pair_dirs <- file.path(JOINT_RESULT_ROOT, manifest$joint_run_prefix)
  require_files(pair_dirs, "Figure 5 C01-C06 joint result directory")
  rows <- lapply(seq_along(pair_dirs), function(pair_index) {
    pair_dir <- pair_dirs[[pair_index]]
    seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    candidates <- lapply(seed_dirs, function(seed_dir) {
      summary_path <- file.path(seed_dir, "fit_summary.tsv")
      if (!file.exists(summary_path)) return(NULL)
      metrics <- read_metric_map(summary_path)
      objective <- suppressWarnings(as.numeric(metrics[["objective"]]))
      if (!is.finite(objective) ||
          !identical(as.character(metrics[["fit_status"]]), "ok")) {
        return(NULL)
      }
      data.frame(
        seed_dir = normalizePath(seed_dir, mustWork = TRUE),
        selected_seed = basename(seed_dir),
        objective = objective,
        objective_invivo =
          suppressWarnings(as.numeric(metrics[["objective_invivo"]])),
        objective_invitro =
          suppressWarnings(as.numeric(metrics[["objective_invitro"]])),
        objective_soft_coupling =
          suppressWarnings(as.numeric(metrics[["objective_soft_coupling"]])),
        stringsAsFactors = FALSE
      )
    })
    candidates <- Filter(Negate(is.null), candidates)
    candidates <- do.call(rbind, candidates)
    candidates <- candidates[
      order(
        candidates$objective,
        as.integer(sub("^seed", "", candidates$selected_seed))
      ),
      ,
      drop = FALSE
    ]
    winner <- candidates[1L, , drop = FALSE]
    pair_id <- basename(pair_dir)
    warmup_label <- sub("^fit_joint_", "", pair_id)
    invivo_seed <- manifest$invivo_seed[[pair_index]]
    invitro_seed <- manifest$invitro_seed[[pair_index]]
    data.frame(
      record_type = "joint_pair_best",
      warmup_label = warmup_label,
      invivo_seed = invivo_seed,
      invitro_seed = invitro_seed,
      selected_seed = winner$selected_seed,
      objective_metric = "objective",
      objective = winner$objective,
      objective_invivo = winner$objective_invivo,
      objective_invitro = winner$objective_invitro,
      objective_soft_coupling = winner$objective_soft_coupling,
      source_dir = winner$seed_dir,
      bundle_dir = file.path(
        "figure5_frozen_inputs", "winners", warmup_label
      ),
      stringsAsFactors = FALSE
    )
  })
  selected <- do.call(rbind, rows)
  selected$family <- joint_family_from_label(selected$warmup_label)
  selected <- selected[
    order(match(selected$family, JOINT_FAMILY_LEVELS)), , drop = FALSE
  ]
  if (!identical(selected$family, JOINT_FAMILY_LEVELS)) {
    stop("Figure 5 joint winners do not cover C01-C06 exactly once.")
  }
  rownames(selected) <- NULL

  add_separate <- function(type, seed, metric, root) {
    metrics <- read_metric_map(
      file.path(root, seed, "fit_summary.tsv")
    )
    value <- suppressWarnings(as.numeric(metrics[[metric]]))
    data.frame(
      record_type = type,
      warmup_label = "",
      invivo_seed = if (type == "separate_invivo_global_best") seed else "",
      invitro_seed = if (type == "separate_invitro_global_best") seed else "",
      selected_seed = seed,
      objective_metric = metric,
      objective = value,
      objective_invivo = NA_real_,
      objective_invitro = NA_real_,
      objective_soft_coupling = NA_real_,
      source_dir = normalizePath(file.path(root, seed), mustWork = TRUE),
      bundle_dir = file.path(
        "separate_fits",
        if (type == "separate_invivo_global_best") "invivo" else "invitro",
        "selected_seeds", seed
      ),
      family = "",
      stringsAsFactors = FALSE
    )
  }
  rbind(
    selected,
    add_separate(
      "separate_invivo_global_best", INVIVO_VISUALIZATION_SEED, "objective",
      INVIVO_RESULT_ROOT
    ),
    add_separate(
      "separate_invitro_global_best", INVITRO_VISUALIZATION_SEED, "objective_total",
      INVITRO_RESULT_ROOT
    )
  )
}

figure5_rebuild_lineage <- function(winner_dir) {
  suppressPackageStartupMessages(library(dplyr))
  source(file.path(
    MODEL_CODE_ROOT, "util", "o2_supply_demand_map_shared.R"
  ))
  source(file.path(
    MODEL_CODE_ROOT, "util", "o2_supply_demand_map_common_semantics.R"
  ))
  source(file.path(
    MODEL_CODE_ROOT, "model", "model_O2_supply_demand_MAP.R"
  ))
  helper_root <- file.path(CODE_ROOT, "util", "invitro")
  for (helper in c(
    "fit_io.R",
    "lineage_adapter.R",
    "simulation.R",
    "summaries.R",
    "objective.R"
  )) {
    source(file.path(helper_root, helper))
  }

  fit <- readRDS(file.path(
    INVITRO_RESULT_ROOT, INVITRO_VISUALIZATION_SEED, "fit_result.rds"
  ))
  fit_data <- list()
  for (run_name in c("run_2N", "run_4N")) {
    adapter <- fit$best_components[[run_name]]$adapter
    for (segment in adapter$segments) {
      for (observed in segment$observed) {
        flow_entry <- if (is.data.frame(observed$observed_flow)) {
          list(
            g0g1_ploidy_density = observed$observed_flow,
            sample_name_canonical = observed$observed_flow_sample_name
          )
        } else {
          NULL
        }
        fit_data[[observed$passage_id]] <- list(
          g = observed$observed_growth,
          passage_duration = observed$passage_duration,
          initial_cells = observed$initial_cells,
          final_cells = observed$final_cells,
          kary = observed$observed_kary,
          flow = flow_entry
        )
      }
    }
  }
  effective <- utils::read.delim(
    file.path(winner_dir, "invitro_effective_params.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  run_params <- fit$best_params
  for (i in seq_len(nrow(effective))) {
    run_params[[effective$parameter[[i]]]] <-
      as.numeric(effective$value[[i]])
  }
  model_core <- build_model_core(cfg = fit$cfg)
  run_2N <- ivt_run_lineage(
    fit$best_components$run_2N$adapter,
    cfg = fit$cfg, run_params = run_params, model_core = model_core
  )
  run_4N <- ivt_run_lineage(
    fit$best_components$run_4N$adapter,
    cfg = fit$cfg, run_params = run_params, model_core = model_core
  )
  summary_df <- dplyr::bind_rows(
    ivt_collect_lineage_summary(run_2N, fit_data),
    ivt_collect_lineage_summary(run_4N, fit_data)
  )
  build_path_map <- function(adapter) {
    segments <- adapter$segments
    ids <- vapply(segments, `[[`, character(1L), "segment_id")
    parents <- vapply(
      segments,
      function(x) {
        parent <- x$parent_segment_id
        if (is.null(parent) || !length(parent) || is.na(parent)) {
          ""
        } else {
          as.character(parent)
        }
      },
      character(1L)
    )
    terminal_keys <- setdiff(ids, parents[nzchar(parents)])
    dplyr::bind_rows(lapply(terminal_keys, function(terminal_key) {
      path <- character()
      current <- terminal_key
      repeat {
        path <- c(current, path)
        idx <- match(current, ids)
        parent <- parents[[idx]]
        if (!nzchar(parent)) break
        current <- parent
      }
      pieces <- strsplit(terminal_key, "_", fixed = TRUE)[[1L]]
      data.frame(
        cohort = segments[[1L]]$cohort,
        lineage_terminal_key = terminal_key,
        lineage_label = if (all(pieces == "20.5")) "control" else "deprived",
        segment_id = path,
        lineage_passage_index = seq_along(path),
        stringsAsFactors = FALSE
      )
    }))
  }
  path_map <- dplyr::bind_rows(
    build_path_map(fit$best_components$run_2N$adapter),
    build_path_map(fit$best_components$run_4N$adapter)
  )
  dplyr::left_join(
    summary_df, path_map, by = c("segment_id", "cohort"),
    relationship = "many-to-many"
  )
}

figure5_rebuild_functional_curves <- function(
    local_winner,
    invitro_fit_result_path
) {
  source(file.path(
    MODEL_CODE_ROOT, "util", "o2_supply_demand_map_common_semantics.R"
  ))
  Sys.setenv(
    MININGCLONEID_OXYGEN_CODE_DIR =
      file.path(MODEL_CODE_ROOT, "model")
  )
  source(file.path(
    MODEL_CODE_ROOT, "model", "model_O2_supply_demand_MAP.R"
  ))
  .first_non_null_local <- function(...) {
    values <- list(...)
    for (value in values) if (!is.null(value)) return(value)
    NULL
  }
  source(file.path(
    CODE_ROOT, "util", "oxygen", "adaptive_grids.R"
  ))
  source(file.path(
    CODE_ROOT, "util", "oxygen", "functional_response.R"
  ))

  best <- utils::read.delim(
    file.path(local_winner, "best_params.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  invivo_cfg <- normalize_sim_cfg_common(
    readRDS(file.path(local_winner, "fit_config.rds")),
    context = "viz"
  )
  invivo_params <- as.list(stats::setNames(
    as.numeric(best$value), best$parameter
  ))
  invivo_params <- normalize_run_params_common(
    invivo_params, cfg = invivo_cfg
  )

  invitro_fit <- readRDS(invitro_fit_result_path)
  invitro_cfg <- normalize_sim_cfg_common(
    invitro_fit$cfg, context = "viz"
  )
  invitro_params <- invitro_fit$best_params
  effective <- utils::read.delim(
    file.path(local_winner, "invitro_effective_params.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  for (i in seq_len(nrow(effective))) {
    invitro_params[[effective$parameter[[i]]]] <-
      as.numeric(effective$value[[i]])
  }
  invitro_params <- normalize_run_params_common(
    invitro_params, cfg = invitro_cfg
  )

  context_tables <- list(
    invivo = generate_invivo_functional_response_tables(
      invivo_params, invivo_cfg
    ),
    invitro = generate_invivo_functional_response_tables(
      invitro_params, invitro_cfg
    )
  )
  table_names <- c(
    "functional_curve_oxygen_multi_ploidy",
    "functional_curve_ploidy"
  )
  for (context in names(context_tables)) {
    context_dir <- file.path(local_winner, "viz", context)
    dir.create(context_dir, recursive = TRUE, showWarnings = FALSE)
    for (table_name in table_names) {
      write_intermediate_tsv(
        context_tables[[context]][[table_name]],
        file.path(context_dir, paste0(table_name, ".tsv"))
      )
    }
  }
  invisible(NULL)
}

# Archived reconstruction of a penalty-induced reference distribution. This
# is retained only for provenance and is not called by the active Figure 5
# pipeline; it is not the DE initial-population distribution used in panel F.
figure5_build_archived_induced_prior_solution_distributions <- function(
    workspace_root,
    master,
    compact,
    parameter_meta,
    n_prior_draws = 100000L,
    primary_seed = 20260813L,
    secondary_seed = 20260913L,
    density_points = 401L
) {
  figure5_dir <- file.path(
    workspace_root, "data", "Figures", "Figure5"
  )
  winner_root <- file.path(
    figure5_dir, "figure5_frozen_inputs", "winners"
  )
  config_paths <- sort(list.files(
    winner_root,
    pattern = "^fit_config\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  ))
  if (length(config_paths) != 6L) {
    stop(
      "Figure 5 prior reconstruction requires six local fit configs; found ",
      length(config_paths)
    )
  }
  configs <- lapply(config_paths, readRDS)
  prior_fields <- c(
    "use_soft_prior", "lambda_prior",
    "prior_center_buffer_smax", "prior_sd_buffer_smax",
    "prior_center_log10_buffer_beta", "prior_sd_log10_buffer_beta",
    "prior_center_log10_buffer_n_exp", "prior_sd_log10_buffer_n_exp",
    "prior_center_log10_mu_hp", "prior_sd_log10_mu_hp",
    "prior_center_gamma_mu", "prior_sd_gamma_mu"
  )
  for (field in prior_fields) {
    values <- vapply(
      configs,
      function(cfg) paste(cfg[[field]], collapse = ","),
      character(1)
    )
    if (length(unique(values)) != 1L) {
      stop("Figure 5 winner configs disagree for prior field: ", field)
    }
  }
  cfg <- configs[[1L]]
  if (!isTRUE(cfg$use_soft_prior) ||
      !is.finite(cfg$lambda_prior) || cfg$lambda_prior <= 0) {
    stop("Figure 5 requires the declared positive in-vivo soft-prior weight")
  }

  coupled_parameters <- c(
    "lam_max", "p_mis_base", "p_wgd",
    "p_misseg", "k_o_mis",
    "O2_crit", "n_O",
    "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  )
  required_master_fields <- c(
    "parameter", "joint_union_lower_transformed",
    "joint_union_upper_transformed", "joint_union_lower_bound",
    "joint_union_upper_bound", "regularization_sigma", "welsch_c"
  )
  if (!all(required_master_fields %in% names(master))) {
    stop(
      "Figure 5 prior reconstruction lacks fields: ",
      paste(setdiff(required_master_fields, names(master)), collapse = ", ")
    )
  }
  bounds <- unique(master[, required_master_fields, drop = FALSE])
  bounds <- bounds[bounds$parameter %in% coupled_parameters, , drop = FALSE]
  if (nrow(bounds) != 14L || anyDuplicated(bounds$parameter)) {
    stop("Figure 5 requires one joint-bound record per coupled parameter")
  }
  bounds <- merge(
    bounds,
    parameter_meta[, c(
      "parameter", "parameter_group", "parameter_order", "transformation"
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  bounds <- bounds[match(coupled_parameters, bounds$parameter), , drop = FALSE]
  if (anyNA(bounds$transformation) ||
      !all(bounds$transformation %in% c("log10", "identity"))) {
    stop("Figure 5 has an unsupported or missing optimizer transformation")
  }
  if (any(!is.finite(bounds$joint_union_lower_transformed)) ||
      any(!is.finite(bounds$joint_union_upper_transformed)) ||
      any(bounds$joint_union_upper_transformed <=
            bounds$joint_union_lower_transformed) ||
      any(bounds$regularization_sigma <= 0) ||
      any(bounds$welsch_c <= 0)) {
    stop("Figure 5 has invalid prior-support or Welsch configuration")
  }

  soft_prior <- data.frame(
    parameter = c(
      "mu_hp", "gamma_mu", "buffer_smax",
      "buffer_beta", "buffer_n_exp"
    ),
    soft_prior_mu = c(
      cfg$prior_center_log10_mu_hp,
      cfg$prior_center_gamma_mu,
      cfg$prior_center_buffer_smax,
      cfg$prior_center_log10_buffer_beta,
      cfg$prior_center_log10_buffer_n_exp
    ),
    soft_prior_sd = c(
      cfg$prior_sd_log10_mu_hp,
      cfg$prior_sd_gamma_mu,
      cfg$prior_sd_buffer_smax,
      cfg$prior_sd_log10_buffer_beta,
      cfg$prior_sd_log10_buffer_n_exp
    ),
    stringsAsFactors = FALSE
  )
  bounds <- merge(
    bounds, soft_prior, by = "parameter", all.x = TRUE, sort = FALSE
  )
  bounds <- bounds[match(coupled_parameters, bounds$parameter), , drop = FALSE]
  bounds$soft_prior_included <- !is.na(bounds$soft_prior_mu)
  bounds$soft_prior_mu[!bounds$soft_prior_included] <- NA_real_
  bounds$soft_prior_sd[!bounds$soft_prior_included] <- NA_real_
  if (any(bounds$soft_prior_included &
          (!is.finite(bounds$soft_prior_sd) | bounds$soft_prior_sd <= 0))) {
    stop("Figure 5 has an invalid configured soft-prior scale")
  }

  sample_one_prior <- function(row, n, seed) {
    set.seed(seed)
    accepted_values <- numeric()
    n_candidates <- 0L
    n_accepted_total <- 0L
    lower <- row$joint_union_lower_transformed
    upper <- row$joint_union_upper_transformed
    sigma <- row$regularization_sigma
    welsch_c <- row$welsch_c
    while (length(accepted_values) < n) {
      remaining <- n - length(accepted_values)
      batch_n <- max(5000L, as.integer(ceiling(remaining / 0.60)))
      vivo_z <- stats::runif(batch_n, lower, upper)
      vitro_z <- stats::runif(batch_n, lower, upper)
      delta <- vivo_z - vitro_z
      penalty <- (welsch_c^2 / 2) * (
        1 - exp(-((delta / (welsch_c * sigma))^2))
      )
      if (isTRUE(row$soft_prior_included)) {
        penalty <- penalty + cfg$lambda_prior * 0.5 * (
          (vivo_z - row$soft_prior_mu) / row$soft_prior_sd
        )^2
      }
      keep <- stats::runif(batch_n) <= exp(-penalty)
      n_candidates <- n_candidates + batch_n
      n_accepted_total <- n_accepted_total + sum(keep)
      if (any(keep)) {
        vivo_keep <- vivo_z[keep]
        vitro_keep <- vitro_z[keep]
        log2_ratio <- if (identical(row$transformation, "log10")) {
          (vivo_keep - vitro_keep) * log2(10)
        } else {
          log2(vivo_keep / vitro_keep)
        }
        log2_ratio <- log2_ratio[is.finite(log2_ratio)]
        accepted_values <- c(accepted_values, log2_ratio)
      }
    }
    list(
      values = accepted_values[seq_len(n)],
      n_candidates = n_candidates,
      n_accepted_total = n_accepted_total,
      acceptance_rate = n_accepted_total / n_candidates
    )
  }

  quantile_probs <- c(0.001, 0.005, 0.05, 0.25, 0.50, 0.75, 0.95, 0.995, 0.999)
  quantile_names <- c(
    "q001", "q005", "q05", "q25", "median",
    "q75", "q95", "q995", "q999"
  )
  prior_draw_rows <- vector("list", nrow(bounds))
  config_rows <- vector("list", nrow(bounds))
  secondary_quantiles <- vector("list", nrow(bounds))
  for (i in seq_len(nrow(bounds))) {
    row <- bounds[i, , drop = FALSE]
    primary <- sample_one_prior(
      row, n_prior_draws, primary_seed + row$parameter_order
    )
    secondary <- sample_one_prior(
      row, n_prior_draws, secondary_seed + row$parameter_order
    )
    primary_q <- stats::quantile(
      primary$values,
      probs = quantile_probs,
      names = FALSE,
      type = 8
    )
    secondary_q <- stats::quantile(
      secondary$values,
      probs = quantile_probs,
      names = FALSE,
      type = 8
    )
    q_scale <- max(1, primary_q[[7L]] - primary_q[[3L]])
    max_q_difference <- max(abs(primary_q - secondary_q))
    prior_draw_rows[[i]] <- data.frame(
      parameter = row$parameter,
      parameter_order = row$parameter_order,
      draw_id = seq_len(n_prior_draws),
      log2_ratio = primary$values,
      stringsAsFactors = FALSE
    )
    config_rows[[i]] <- data.frame(
      parameter = row$parameter,
      parameter_group = row$parameter_group,
      parameter_order = row$parameter_order,
      transformation = row$transformation,
      joint_lower_transformed = row$joint_union_lower_transformed,
      joint_upper_transformed = row$joint_union_upper_transformed,
      joint_lower_natural = row$joint_union_lower_bound,
      joint_upper_natural = row$joint_union_upper_bound,
      base_measure = paste0(
        "uniform context values on joint transformed bounds; ",
        "pushed forward to log2(in vivo/in vitro)"
      ),
      prior_definition = if (row$soft_prior_included) {
        "joint bounds + Welsch coupling + in-vivo Gaussian soft prior"
      } else {
        "joint bounds + Welsch coupling; no parameter-specific Gaussian prior"
      },
      soft_prior_included = row$soft_prior_included,
      soft_prior_mu = row$soft_prior_mu,
      soft_prior_sd = row$soft_prior_sd,
      lambda_prior = cfg$lambda_prior,
      regularization_sigma = row$regularization_sigma,
      welsch_c = row$welsch_c,
      primary_seed = primary_seed + row$parameter_order,
      secondary_seed = secondary_seed + row$parameter_order,
      retained_draws = n_prior_draws,
      primary_candidates = primary$n_candidates,
      primary_acceptance_rate = primary$acceptance_rate,
      secondary_candidates = secondary$n_candidates,
      secondary_acceptance_rate = secondary$acceptance_rate,
      max_replicate_quantile_difference = max_q_difference,
      scaled_replicate_quantile_difference = max_q_difference / q_scale,
      stringsAsFactors = FALSE
    )
    secondary_quantiles[[i]] <- data.frame(
      parameter = row$parameter,
      quantile = quantile_names,
      primary = primary_q,
      secondary = secondary_q,
      absolute_difference = abs(primary_q - secondary_q),
      stringsAsFactors = FALSE
    )
  }
  prior_draws <- do.call(rbind, prior_draw_rows)
  prior_config <- do.call(rbind, config_rows)
  prior_replicate_quantiles <- do.call(rbind, secondary_quantiles)
  rownames(prior_draws) <- NULL
  rownames(prior_config) <- NULL
  rownames(prior_replicate_quantiles) <- NULL

  summarize_values <- function(values) {
    q <- stats::quantile(
      values,
      probs = quantile_probs,
      names = FALSE,
      type = 8
    )
    out <- as.list(q)
    names(out) <- quantile_names
    out$width90 <- q[[7L]] - q[[3L]]
    out$n <- length(values)
    as.data.frame(out, stringsAsFactors = FALSE)
  }
  prior_summary_rows <- lapply(coupled_parameters, function(parameter) {
    cbind(
      data.frame(
        parameter = parameter,
        distribution = "Induced prior/reference",
        distribution_role = "prior",
        family = "",
        stringsAsFactors = FALSE
      ),
      summarize_values(
        prior_draws$log2_ratio[prior_draws$parameter == parameter]
      )
    )
  })
  solution_summary_rows <- lapply(coupled_parameters, function(parameter) {
    values <- compact[compact$parameter == parameter, , drop = FALSE]
    pooled <- cbind(
      data.frame(
        parameter = parameter,
        distribution = "All optimizer solutions",
        distribution_role = "solution_pooled",
        family = "All",
        stringsAsFactors = FALSE
      ),
      summarize_values(values$log2_ratio)
    )
    families <- do.call(rbind, lapply(
      c("C01", "C02", "C03"),
      function(family) {
        family_values <- values$log2_ratio[values$family == family]
        cbind(
          data.frame(
            parameter = parameter,
            distribution = family,
            distribution_role = "solution_family",
            family = family,
            stringsAsFactors = FALSE
          ),
          summarize_values(family_values)
        )
      }
    ))
    rbind(pooled, families)
  })
  prior_solution_summary <- rbind(
    do.call(rbind, prior_summary_rows),
    do.call(rbind, solution_summary_rows)
  )
  prior_solution_summary <- merge(
    prior_solution_summary,
    parameter_meta[, c(
      "parameter", "parameter_group", "parameter_order"
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  prior_width <- prior_solution_summary$width90[
    match(
      prior_solution_summary$parameter,
      prior_solution_summary$parameter[
        prior_solution_summary$distribution_role == "prior"
      ]
    )
  ]
  prior_width_lookup <- setNames(
    prior_solution_summary$width90[
      prior_solution_summary$distribution_role == "prior"
    ],
    prior_solution_summary$parameter[
      prior_solution_summary$distribution_role == "prior"
    ]
  )
  prior_solution_summary$width90_relative_to_prior <-
    prior_solution_summary$width90 /
    prior_width_lookup[prior_solution_summary$parameter]
  prior_solution_summary <- prior_solution_summary[order(
    prior_solution_summary$parameter_order,
    match(
      prior_solution_summary$distribution_role,
      c("prior", "solution_pooled", "solution_family")
    ),
    match(prior_solution_summary$family, c("", "All", "C01", "C02", "C03"))
  ), , drop = FALSE]
  rownames(prior_solution_summary) <- NULL

  seed_groups <- split(
    compact,
    interaction(compact$pair_id, compact$seed_number, drop = TRUE)
  )
  endpoint_signatures <- do.call(rbind, lapply(seed_groups, function(group) {
    group <- group[order(group$parameter_order), , drop = FALSE]
    data.frame(
      pair_id = group$pair_id[[1L]],
      family = group$family[[1L]],
      seed_number = group$seed_number[[1L]],
      objective_rank = group$objective_rank[[1L]],
      objective = group$objective[[1L]],
      signature = paste(
        c(
          sprintf("%.17g", group$vivo_natural),
          sprintf("%.17g", group$vitro_natural)
        ),
        collapse = "|"
      ),
      stringsAsFactors = FALSE
    )
  }))
  endpoint_signatures <- endpoint_signatures[order(
    endpoint_signatures$pair_id,
    endpoint_signatures$objective,
    endpoint_signatures$seed_number
  ), , drop = FALSE]
  endpoint_signatures$pair_signature <- paste(
    endpoint_signatures$pair_id,
    endpoint_signatures$signature,
    sep = "::"
  )
  unique_endpoint_keys <- endpoint_signatures[
    !duplicated(endpoint_signatures$pair_signature),
    c("pair_id", "seed_number"),
    drop = FALSE
  ]
  subset_keys <- list(
    all_endpoints = endpoint_signatures[, c("pair_id", "seed_number")],
    lowest_10pct_per_pair = endpoint_signatures[
      endpoint_signatures$objective_rank <= 50L,
      c("pair_id", "seed_number"),
      drop = FALSE
    ],
    unique_endpoints = unique_endpoint_keys
  )
  sensitivity_rows <- list()
  sensitivity_index <- 1L
  for (subset_name in names(subset_keys)) {
    keys <- subset_keys[[subset_name]]
    keys$key <- paste(keys$pair_id, keys$seed_number, sep = "::")
    compact$key <- paste(compact$pair_id, compact$seed_number, sep = "::")
    subset_data <- compact[compact$key %in% keys$key, , drop = FALSE]
    for (parameter in coupled_parameters) {
      parameter_data <- subset_data[
        subset_data$parameter == parameter,
        ,
        drop = FALSE
      ]
      pooled <- summarize_values(parameter_data$log2_ratio)
      family_summary <- do.call(rbind, lapply(
        c("C01", "C02", "C03"),
        function(family) {
          q <- stats::quantile(
            parameter_data$log2_ratio[parameter_data$family == family],
            probs = c(0.05, 0.95),
            names = FALSE,
            type = 8
          )
          data.frame(family = family, q05 = q[[1L]], q95 = q[[2L]])
        }
      ))
      direction <- if (all(family_summary$q05 > 0)) {
        "higher in vivo"
      } else if (all(family_summary$q95 < 0)) {
        "lower in vivo"
      } else {
        "family-dependent"
      }
      sensitivity_rows[[sensitivity_index]] <- data.frame(
        subset = subset_name,
        parameter = parameter,
        n_endpoints = length(unique(parameter_data$key)),
        q05 = pooled$q05,
        median = pooled$median,
        q95 = pooled$q95,
        width90 = pooled$width90,
        direction = direction,
        width90_relative_to_prior =
          pooled$width90 / prior_width_lookup[[parameter]],
        stringsAsFactors = FALSE
      )
      sensitivity_index <- sensitivity_index + 1L
    }
  }
  sensitivity <- do.call(rbind, sensitivity_rows)
  sensitivity <- merge(
    sensitivity,
    parameter_meta[, c(
      "parameter", "parameter_group", "parameter_order"
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  sensitivity <- sensitivity[order(
    sensitivity$parameter_order,
    match(
      sensitivity$subset,
      c("all_endpoints", "lowest_10pct_per_pair", "unique_endpoints")
    )
  ), , drop = FALSE]
  direction_sets <- split(sensitivity$direction, sensitivity$parameter)
  direction_stability <- vapply(
    direction_sets,
    function(x) length(unique(x)) == 1L,
    logical(1)
  )
  sensitivity$direction_stable_across_sensitivities <-
    direction_stability[sensitivity$parameter]
  rownames(sensitivity) <- NULL

  reflected_density <- function(values, from, to, support_lower, support_upper) {
    values <- values[is.finite(values)]
    if (length(values) < 2L || length(unique(values)) < 2L) {
      x <- seq(from, to, length.out = density_points)
      y <- rep(0, length(x))
      y[[which.min(abs(x - values[[1L]]))]] <- 1
      return(data.frame(x = x, density = y, density_scaled = y))
    }
    bw <- stats::bw.nrd0(values)
    if (!is.finite(bw) || bw <= 0) {
      bw <- max((to - from) / 100, .Machine$double.eps^0.25)
    }
    augmented <- values
    multiplier <- 1L
    if (is.finite(support_lower)) {
      augmented <- c(augmented, 2 * support_lower - values)
      multiplier <- multiplier + 1L
    }
    if (is.finite(support_upper)) {
      augmented <- c(augmented, 2 * support_upper - values)
      multiplier <- multiplier + 1L
    }
    d <- stats::density(
      augmented,
      bw = bw,
      from = from,
      to = to,
      n = density_points,
      cut = 0
    )
    y <- d$y * multiplier
    y_scaled <- if (max(y) > 0) y / max(y) else y
    data.frame(x = d$x, density = y, density_scaled = y_scaled)
  }

  density_rows <- list()
  density_index <- 1L
  for (i in seq_len(nrow(bounds))) {
    row <- bounds[i, , drop = FALSE]
    parameter <- row$parameter
    prior_values <- prior_draws$log2_ratio[
      prior_draws$parameter == parameter
    ]
    solution_data <- compact[compact$parameter == parameter, , drop = FALSE]
    support <- if (identical(row$transformation, "log10")) {
      width <- (
        row$joint_union_upper_transformed -
          row$joint_union_lower_transformed
      ) * log2(10)
      c(-width, width)
    } else if (row$joint_union_lower_transformed > 0) {
      width <- log2(
        row$joint_union_upper_transformed /
          row$joint_union_lower_transformed
      )
      c(-width, width)
    } else {
      c(-Inf, Inf)
    }
    prior_q <- stats::quantile(
      prior_values, c(0.001, 0.999), names = FALSE, type = 8
    )
    display_from <- min(prior_q[[1L]], min(solution_data$log2_ratio), 0)
    display_to <- max(prior_q[[2L]], max(solution_data$log2_ratio), 0)
    padding <- max((display_to - display_from) * 0.025, 0.02)
    display_from <- max(display_from - padding, support[[1L]])
    display_to <- min(display_to + padding, support[[2L]])
    distributions <- c(
      list(
        list(
          distribution = "Induced prior/reference",
          distribution_role = "prior",
          family = "",
          values = prior_values
        ),
        list(
          distribution = "All optimizer solutions",
          distribution_role = "solution_pooled",
          family = "All",
          values = solution_data$log2_ratio
        )
      ),
      lapply(c("C01", "C02", "C03"), function(family) {
        list(
          distribution = family,
          distribution_role = "solution_family",
          family = family,
          values = solution_data$log2_ratio[
            solution_data$family == family
          ]
        )
      })
    )
    for (distribution in distributions) {
      density <- reflected_density(
        distribution$values,
        display_from,
        display_to,
        support[[1L]],
        support[[2L]]
      )
      density_rows[[density_index]] <- cbind(
        data.frame(
          parameter = parameter,
          parameter_group = row$parameter_group,
          parameter_order = row$parameter_order,
          distribution = distribution$distribution,
          distribution_role = distribution$distribution_role,
          family = distribution$family,
          display_from = display_from,
          display_to = display_to,
          support_lower = support[[1L]],
          support_upper = support[[2L]],
          displayed_tail_fraction = mean(
            distribution$values < display_from |
              distribution$values > display_to
          ),
          stringsAsFactors = FALSE
        ),
        density
      )
      density_index <- density_index + 1L
    }
  }
  density_table <- do.call(rbind, density_rows)
  rownames(density_table) <- NULL

  if (nrow(prior_draws) != 14L * n_prior_draws ||
      nrow(prior_config) != 14L ||
      nrow(prior_solution_summary) != 70L ||
      nrow(sensitivity) != 42L ||
      nrow(density_table) != 14L * 5L * density_points ||
      any(prior_config$primary_acceptance_rate <= 0.50) ||
      any(prior_config$scaled_replicate_quantile_difference > 0.05) ||
      any(density_table$displayed_tail_fraction > 0.002) ||
      any(!is.finite(density_table$density_scaled))) {
    stop("Figure 5 prior/solution distribution validation failed")
  }

  write_intermediate_tsv(
    prior_draws,
    file.path(figure5_dir, "figure5f_induced_prior_draws.tsv")
  )
  write_intermediate_tsv(
    density_table,
    file.path(figure5_dir, "figure5f_prior_solution_density.tsv")
  )
  write_intermediate_tsv(
    prior_solution_summary,
    file.path(figure5_dir, "figure5f_prior_solution_summary.tsv")
  )
  write_intermediate_tsv(
    prior_config,
    file.path(figure5_dir, "figure5f_prior_sampling_config.tsv")
  )
  write_intermediate_tsv(
    prior_replicate_quantiles,
    file.path(figure5_dir, "figure5f_prior_replicate_quantiles.tsv")
  )
  write_intermediate_tsv(
    sensitivity,
    file.path(figure5_dir, "figure5f_sensitivity_validation.tsv")
  )
  invisible(list(
    draws = prior_draws,
    density = density_table,
    summary = prior_solution_summary,
    config = prior_config,
    sensitivity = sensitivity,
    config_paths = config_paths
  ))
}

# Figure 5F no longer displays an induced prior/reference. This helper builds
# only the optimizer-solution products used by the nested C01-C06 panels.
# Within a parameter, all six families use the same horizontal grid, KDE
# bandwidth, and density normalization constant so their locations, widths,
# and peak heights remain directly comparable. Different parameters retain
# parameter-specific horizontal ranges.
figure5_build_family_distribution_products <- function(
    workspace_root,
    compact,
    parameter_meta,
    density_points = 401L
) {
  figure5_dir <- file.path(
    workspace_root, "data", "Figures", "Figure5"
  )
  coupled_parameters <- c(
    "lam_max", "p_mis_base", "p_wgd",
    "p_misseg", "k_o_mis",
    "O2_crit", "n_O",
    "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  )
  family_levels <- JOINT_FAMILY_LEVELS

  summarize_values <- function(values) {
    q <- stats::quantile(
      values,
      probs = c(0.05, 0.50, 0.95),
      names = FALSE,
      type = 8
    )
    data.frame(
      q05 = q[[1L]],
      median = q[[2L]],
      q95 = q[[3L]],
      width90 = q[[3L]] - q[[1L]],
      stringsAsFactors = FALSE
    )
  }

  seed_groups <- split(
    compact,
    interaction(compact$pair_id, compact$seed_number, drop = TRUE)
  )
  endpoint_signatures <- do.call(rbind, lapply(seed_groups, function(group) {
    group <- group[order(group$parameter_order), , drop = FALSE]
    data.frame(
      pair_id = group$pair_id[[1L]],
      family = group$family[[1L]],
      seed_number = group$seed_number[[1L]],
      objective_rank = group$objective_rank[[1L]],
      objective = group$objective[[1L]],
      signature = paste(
        c(
          sprintf("%.17g", group$vivo_natural),
          sprintf("%.17g", group$vitro_natural)
        ),
        collapse = "|"
      ),
      stringsAsFactors = FALSE
    )
  }))
  endpoint_signatures <- endpoint_signatures[order(
    endpoint_signatures$pair_id,
    endpoint_signatures$objective,
    endpoint_signatures$seed_number
  ), , drop = FALSE]
  endpoint_signatures$pair_signature <- paste(
    endpoint_signatures$pair_id,
    endpoint_signatures$signature,
    sep = "::"
  )
  unique_endpoint_keys <- endpoint_signatures[
    !duplicated(endpoint_signatures$pair_signature),
    c("pair_id", "seed_number"),
    drop = FALSE
  ]
  subset_keys <- list(
    all_endpoints = endpoint_signatures[, c("pair_id", "seed_number")],
    lowest_10pct_per_pair = endpoint_signatures[
      endpoint_signatures$objective_rank <= 50L,
      c("pair_id", "seed_number"),
      drop = FALSE
    ],
    unique_endpoints = unique_endpoint_keys
  )
  compact$key <- paste(compact$pair_id, compact$seed_number, sep = "::")
  sensitivity_rows <- list()
  sensitivity_index <- 1L
  for (subset_name in names(subset_keys)) {
    keys <- subset_keys[[subset_name]]
    keys$key <- paste(keys$pair_id, keys$seed_number, sep = "::")
    subset_data <- compact[compact$key %in% keys$key, , drop = FALSE]
    for (parameter in coupled_parameters) {
      parameter_data <- subset_data[
        subset_data$parameter == parameter,
        ,
        drop = FALSE
      ]
      pooled <- summarize_values(parameter_data$log2_ratio)
      family_intervals <- do.call(rbind, lapply(
        family_levels,
        function(family) {
          family_values <- parameter_data$log2_ratio[
            parameter_data$family == family
          ]
          q <- stats::quantile(
            family_values,
            probs = c(0.05, 0.95),
            names = FALSE,
            type = 8
          )
          data.frame(
            family = family,
            q05 = q[[1L]],
            q95 = q[[2L]],
            width90 = q[[2L]] - q[[1L]],
            stringsAsFactors = FALSE
          )
        }
      ))
      direction <- if (all(family_intervals$q05 > 0)) {
        "higher in vivo"
      } else if (all(family_intervals$q95 < 0)) {
        "lower in vivo"
      } else {
        "family-dependent"
      }
      sensitivity_rows[[sensitivity_index]] <- data.frame(
        subset = subset_name,
        parameter = parameter,
        n_endpoints = length(unique(parameter_data$key)),
        q05 = pooled$q05,
        median = pooled$median,
        q95 = pooled$q95,
        width90 = pooled$width90,
        min_family_width90 = min(family_intervals$width90),
        max_family_width90 = max(family_intervals$width90),
        direction = direction,
        stringsAsFactors = FALSE
      )
      sensitivity_index <- sensitivity_index + 1L
    }
  }
  sensitivity <- do.call(rbind, sensitivity_rows)
  sensitivity <- merge(
    sensitivity,
    parameter_meta[, c(
      "parameter", "parameter_group", "parameter_order"
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  sensitivity <- sensitivity[order(
    sensitivity$parameter_order,
    match(
      sensitivity$subset,
      c("all_endpoints", "lowest_10pct_per_pair", "unique_endpoints")
    )
  ), , drop = FALSE]
  direction_sets <- split(sensitivity$direction, sensitivity$parameter)
  direction_stability <- vapply(
    direction_sets,
    function(x) length(unique(x)) == 1L,
    logical(1)
  )
  sensitivity$direction_stable_across_sensitivities <-
    direction_stability[sensitivity$parameter]
  rownames(sensitivity) <- NULL

  density_rows <- list()
  density_index <- 1L
  for (parameter in coupled_parameters) {
    parameter_data <- compact[
      compact$parameter == parameter,
      ,
      drop = FALSE
    ]
    values <- parameter_data$log2_ratio
    common_bandwidth <- stats::bw.nrd0(values)
    raw_span <- diff(range(values))
    if (!is.finite(common_bandwidth) || common_bandwidth <= 0) {
      common_bandwidth <- max(raw_span / 100, 0.02)
    }
    display_core <- range(c(values, 0))
    display_span <- diff(display_core)
    padding <- max(
      display_span * 0.04,
      3 * common_bandwidth,
      0.02
    )
    display_from <- display_core[[1L]] - padding
    display_to <- display_core[[2L]] + padding
    for (family in family_levels) {
      family_values <- parameter_data$log2_ratio[
        parameter_data$family == family
      ]
      if (length(family_values) != 500L) {
        stop(
          "Figure 5F requires 500 endpoints for ",
          parameter, " / ", family
        )
      }
      density <- stats::density(
        family_values,
        bw = common_bandwidth,
        from = display_from,
        to = display_to,
        n = density_points,
        cut = 0
      )
      density_rows[[density_index]] <- data.frame(
        parameter = parameter,
        parameter_group = parameter_data$parameter_group[[1L]],
        parameter_order = parameter_data$parameter_order[[1L]],
        family = family,
        n_solutions = length(family_values),
        common_bandwidth = common_bandwidth,
        display_from = display_from,
        display_to = display_to,
        grid_index = seq_along(density$x),
        x = density$x,
        density = density$y,
        stringsAsFactors = FALSE
      )
      density_index <- density_index + 1L
    }
  }
  density_table <- do.call(rbind, density_rows)
  parameter_density_max <- tapply(
    density_table$density,
    density_table$parameter,
    max
  )
  density_table$density_scaled_parameter <-
    density_table$density /
    parameter_density_max[density_table$parameter]
  density_table <- density_table[order(
    density_table$parameter_order,
    match(density_table$family, family_levels),
    density_table$grid_index
  ), , drop = FALSE]
  rownames(density_table) <- NULL

  density_counts <- table(
    density_table$parameter,
    density_table$family
  )
  common_scale_check <- vapply(
    split(density_table, density_table$parameter),
    function(group) {
      length(unique(group$display_from)) == 1L &&
        length(unique(group$display_to)) == 1L &&
        length(unique(group$common_bandwidth)) == 1L &&
        all(vapply(
          split(group$x, group$family),
          function(x) isTRUE(all.equal(x, group$x[group$family == "C01"])),
          logical(1)
        ))
    },
    logical(1)
  )
  if (nrow(density_table) != 14L * 6L * density_points ||
      any(density_counts != density_points) ||
      !all(common_scale_check) ||
      any(!is.finite(density_table$density_scaled_parameter)) ||
      any(density_table$density_scaled_parameter < 0) ||
      any(density_table$density_scaled_parameter > 1 + 1e-12) ||
      nrow(sensitivity) != 42L ||
      anyNA(sensitivity$direction_stable_across_sensitivities)) {
    stop("Figure 5F family-distribution validation failed")
  }

  write_intermediate_tsv(
    density_table,
    file.path(figure5_dir, "figure5f_family_density.tsv")
  )
  write_intermediate_tsv(
    sensitivity,
    file.path(figure5_dir, "figure5f_sensitivity_validation.tsv")
  )
  invisible(list(
    density = density_table,
    sensitivity = sensitivity
  ))
}

figure5_build_local_solution_ensemble <- function(
    workspace_root = WORKSPACE_ROOT
) {
  figure5_dir <- file.path(
    workspace_root, "data", "Figures", "Figure5"
  )
  supp5_1_master_path <- file.path(
    workspace_root, "data", "Figures", "Supp_Figure5_1",
    "soft_coupling_master_long.tsv"
  )
  supp5_2_joint_root <- file.path(
    workspace_root, "data", "Figures", "Supp_Figure5_2", "joint"
  )
  parameter_meta_path <- file.path(
    figure5_dir, "parameter_function_groups.tsv"
  )
  require_files(
    c(supp5_1_master_path, parameter_meta_path),
    "Figure 5 local optimizer-ensemble input"
  )
  if (!dir.exists(supp5_2_joint_root)) {
    stop("Missing Figure 5 local objective root: ", supp5_2_joint_root)
  }

  master <- utils::read.delim(
    supp5_1_master_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_master_fields <- c(
    "parameter", "vivo_natural", "vitro_natural",
    "ratio_vivo_to_vitro", "joint_union_lower_bound",
    "joint_union_upper_bound", "feasible_at_solution",
    "projection_applied", "seed_number", "pair_id"
  )
  if (!all(required_master_fields %in% names(master))) {
    stop(
      "Figure 5 optimizer ensemble lacks fields: ",
      paste(
        setdiff(required_master_fields, names(master)),
        collapse = ", "
      )
    )
  }

  objective_paths <- sort(list.files(
    supp5_2_joint_root,
    pattern = "^seed_objective_simple\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  ))
  if (length(objective_paths) != 6L) {
    stop(
      "Figure 5 requires six local joint objective tables; found ",
      length(objective_paths)
    )
  }
  objective_rows <- lapply(objective_paths, function(path) {
    objective <- utils::read.delim(
      path, check.names = FALSE, stringsAsFactors = FALSE
    )
    if (!all(c("seed", "objective_rank", "objective") %in% names(objective))) {
      stop("Malformed joint objective table: ", path)
    }
    objective$seed_number <- as.integer(sub("^seed", "", objective$seed))
    objective$warmup_label <- basename(dirname(path))
    objective$pair_id <- paste0("fit_joint_", objective$warmup_label)
    objective$pair_objective_min <- min(objective$objective, na.rm = TRUE)
    objective$delta_objective <-
      objective$objective - objective$pair_objective_min
    objective[, c(
      "pair_id", "seed_number", "objective_rank", "objective",
      "pair_objective_min", "delta_objective"
    )]
  })
  objectives <- do.call(rbind, objective_rows)
  if (nrow(objectives) != 3000L ||
      any(table(objectives$pair_id) != 500L)) {
    stop("Figure 5 local objective universe must contain 500 seeds per pair")
  }

  master$warmup_label <- sub("^fit_joint_", "", master$pair_id)
  master$family <- regmatches(
    master$pair_id, regexpr("C[0-9]{2}", master$pair_id)
  )
  master$subcluster <- ""
  master$log2_ratio <- log2(master$ratio_vivo_to_vitro)
  master <- merge(
    master,
    objectives,
    by = c("pair_id", "seed_number"),
    all.x = TRUE,
    sort = FALSE
  )
  if (any(!is.finite(master$objective)) ||
      any(!is.finite(master$log2_ratio))) {
    stop("Figure 5 optimizer ensemble has unmatched objectives or ratios")
  }

  parameter_meta <- utils::read.delim(
    parameter_meta_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  master <- merge(
    master,
    parameter_meta[, c(
      "parameter", "parameter_group", "parameter_order"
    )],
    by = "parameter",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(master$parameter_group)) {
    stop("Figure 5 optimizer ensemble has unmapped parameters")
  }

  bound_scale <- pmax(
    1,
    abs(master$joint_union_lower_bound),
    abs(master$joint_union_upper_bound),
    abs(
      master$joint_union_upper_bound -
        master$joint_union_lower_bound
    ),
    na.rm = TRUE
  )
  bound_tolerance <- 1e-8 * bound_scale
  master$vivo_at_active_bound <-
    abs(master$vivo_natural - master$joint_union_lower_bound) <=
      bound_tolerance |
    abs(master$vivo_natural - master$joint_union_upper_bound) <=
      bound_tolerance
  master$vitro_at_active_bound <-
    abs(master$vitro_natural - master$joint_union_lower_bound) <=
      bound_tolerance |
    abs(master$vitro_natural - master$joint_union_upper_bound) <=
      bound_tolerance
  master$either_context_at_active_bound <-
    master$vivo_at_active_bound | master$vitro_at_active_bound

  master <- master[order(
    master$parameter_order,
    master$family,
    master$subcluster,
    master$seed_number
  ), , drop = FALSE]
  compact <- master[, c(
    "family", "subcluster", "pair_id", "warmup_label",
    "seed_number", "objective_rank", "objective",
    "pair_objective_min", "delta_objective", "parameter",
    "parameter_group", "parameter_order", "vivo_natural",
    "vitro_natural", "ratio_vivo_to_vitro", "log2_ratio",
    "feasible_at_solution", "projection_applied",
    "vivo_at_active_bound", "vitro_at_active_bound",
    "either_context_at_active_bound"
  )]

  summarize_solution_group <- function(data) {
    values <- data$log2_ratio
    quantiles <- stats::quantile(
      values,
      probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
      names = FALSE,
      type = 8
    )
    data.frame(
      n_solutions = nrow(data),
      q05 = quantiles[[1L]],
      q25 = quantiles[[2L]],
      median = quantiles[[3L]],
      q75 = quantiles[[4L]],
      q95 = quantiles[[5L]],
      fraction_lower_in_vivo = mean(values < 0),
      fraction_equal = mean(abs(values) <= 1e-12),
      fraction_higher_in_vivo = mean(values > 0),
      active_bound_fraction = mean(
        data$either_context_at_active_bound
      ),
      objective_min = min(data$objective),
      objective_max = max(data$objective),
      objective_span = max(data$objective) - min(data$objective),
      stringsAsFactors = FALSE
    )
  }
  summarize_by <- function(data, keys) {
    groups <- split(
      data,
      interaction(data[, keys, drop = FALSE], drop = TRUE)
    )
    do.call(rbind, lapply(groups, function(group) {
      cbind(
        group[1L, keys, drop = FALSE],
        summarize_solution_group(group)
      )
    }))
  }

  pair_summary <- summarize_by(
    compact,
    c(
      "family", "subcluster", "pair_id", "parameter",
      "parameter_group", "parameter_order"
    )
  )
  pair_summary <- pair_summary[order(
    pair_summary$parameter_order,
    pair_summary$family,
    pair_summary$subcluster
  ), , drop = FALSE]
  rownames(pair_summary) <- NULL

  family_summary <- summarize_by(
    compact,
    c(
      "family", "parameter", "parameter_group",
      "parameter_order"
    )
  )
  family_pair_groups <- split(
    pair_summary,
    interaction(
      pair_summary$family,
      pair_summary$parameter,
      drop = TRUE
    )
  )
  family_pair_metrics <- do.call(rbind, lapply(
    family_pair_groups,
    function(group) {
      data.frame(
        family = group$family[[1L]],
        parameter = group$parameter[[1L]],
        n_pairs = nrow(group),
        pair_median_min = min(group$median),
        pair_median_max = max(group$median),
        pair_median_span = max(group$median) - min(group$median),
        all_pair_medians_lower = all(group$median < 0),
        all_pair_medians_higher = all(group$median > 0),
        stringsAsFactors = FALSE
      )
    }
  ))
  family_summary <- merge(
    family_summary,
    family_pair_metrics,
    by = c("family", "parameter"),
    all.x = TRUE,
    sort = FALSE
  )
  family_summary$family_direction <- ifelse(
    family_summary$q05 > 0,
    "higher in vivo",
    ifelse(
      family_summary$q95 < 0,
      "lower in vivo",
      "crosses 1"
    )
  )
  family_summary <- family_summary[order(
    family_summary$parameter_order,
    family_summary$family
  ), , drop = FALSE]
  rownames(family_summary) <- NULL

  cross_groups <- split(family_summary, family_summary$parameter)
  cross_family <- do.call(rbind, lapply(cross_groups, function(group) {
    parameter_pair <- pair_summary[
      pair_summary$parameter == group$parameter[[1L]],
      ,
      drop = FALSE
    ]
    stable_higher <- all(group$q05 > 0)
    stable_lower <- all(group$q95 < 0)
    data.frame(
      parameter = group$parameter[[1L]],
      parameter_group = group$parameter_group[[1L]],
      parameter_order = group$parameter_order[[1L]],
      n_families = nrow(group),
      n_pairs = nrow(parameter_pair),
      n_family_higher = sum(group$family_direction == "higher in vivo"),
      n_family_lower = sum(group$family_direction == "lower in vivo"),
      n_family_crosses_1 = sum(group$family_direction == "crosses 1"),
      n_pair_medians_higher = sum(parameter_pair$median > 0),
      n_pair_medians_lower = sum(parameter_pair$median < 0),
      directionally_stable = stable_higher || stable_lower,
      direction = if (stable_higher) {
        "higher in vivo"
      } else if (stable_lower) {
        "lower in vivo"
      } else {
        "family-dependent"
      },
      cross_family_median_min = min(group$median),
      cross_family_median_max = max(group$median),
      cross_family_median_span = max(group$median) - min(group$median),
      max_family_90_width = max(group$q95 - group$q05),
      max_active_bound_fraction = max(group$active_bound_fraction),
      stringsAsFactors = FALSE
    )
  }))
  cross_family <- cross_family[order(
    cross_family$parameter_order
  ), , drop = FALSE]
  rownames(cross_family) <- NULL

  if (nrow(compact) != 42000L ||
      length(unique(compact$pair_id)) != 6L ||
      length(unique(compact$parameter)) != 14L ||
      length(unique(compact$family)) != 6L ||
      any(table(compact$pair_id, compact$parameter) != 500L) ||
      any(!compact$feasible_at_solution) ||
      any(compact$projection_applied)) {
    stop("Figure 5 optimizer-ensemble validation failed")
  }

  write_intermediate_tsv(
    compact,
    file.path(figure5_dir, "figure5f_optimizer_solutions.tsv")
  )
  write_intermediate_tsv(
    pair_summary,
    file.path(figure5_dir, "figure5f_pair_summary.tsv")
  )
  write_intermediate_tsv(
    family_summary,
    file.path(figure5_dir, "figure5f_family_summary.tsv")
  )
  write_intermediate_tsv(
    cross_family,
    file.path(figure5_dir, "figure5f_cross_family_summary.tsv")
  )
  family_distribution <- figure5_build_family_distribution_products(
    workspace_root = workspace_root,
    compact = compact,
    parameter_meta = parameter_meta
  )
  invisible(list(
    solutions = compact,
    pair = pair_summary,
    family = family_summary,
    cross_family = cross_family,
    family_distribution = family_distribution,
    input_paths = c(supp5_1_master_path, objective_paths)
  ))
}

data_Figure5 <- function() {
  data_dir <- file.path(DATA_ROOT, "Figure5")
  frozen_root <- file.path(data_dir, "figure5_frozen_inputs")
  winners_root <- file.path(frozen_root, "winners")
  dir.create(winners_root, recursive = TRUE, showWarnings = FALSE)

  selection <- figure5_select_winners()
  selection_path <- file.path(frozen_root, "selected_results.tsv")
  write_intermediate_tsv(selection, selection_path)
  joint_selection <- selection[
    selection$record_type == "joint_pair_best", , drop = FALSE
  ]
  selected_pair_inputs <- joint_selection
  selected_pair_inputs$subcluster <- ""
  selected_pair_inputs$pair_id <- paste0(
    "fit_joint_", selected_pair_inputs$warmup_label
  )
  selected_pair_inputs$separate_invivo_objective <- vapply(
    selected_pair_inputs$invivo_seed,
    function(seed) {
      metrics <- read_metric_map(file.path(
        INVIVO_RESULT_ROOT, seed, "fit_summary.tsv"
      ))
      suppressWarnings(as.numeric(metrics[["objective"]]))
    },
    numeric(1L)
  )
  selected_pair_inputs$separate_invivo_objective_rank <- 1L
  selected_pair_inputs$family_selection_rank <- 1L
  selected_pair_inputs$selected_for_figure5f <- TRUE
  selected_pair_inputs$selection_rule <- paste0(
    "one objective-minimizing representative from each primary in-vivo ",
    "cluster C01-C06, paired to global in-vitro best ",
    INVITRO_VISUALIZATION_SEED
  )
  selected_pair_inputs <- selected_pair_inputs[
    order(match(selected_pair_inputs$family, JOINT_FAMILY_LEVELS)), ,
    drop = FALSE
  ]
  write_intermediate_tsv(
    selected_pair_inputs,
    file.path(data_dir, "figure5f_selected_pair_inputs.tsv")
  )

  copied_sources <- character()
  copied_destinations <- character()
  copied_roles <- character()
  table_names <- c(
    "invitro_growth_loglik.tsv",
    "invivo_burden_fit.tsv",
    "invivo_terminal_ploidy_fit.tsv",
    "joint_soft_coupling.tsv",
    "invitro_effective_params.tsv",
    "best_params.tsv",
    "fit_config.rds"
  )
  invitro_fit_result_path <- file.path(
    INVITRO_RESULT_ROOT, INVITRO_VISUALIZATION_SEED, "fit_result.rds"
  )
  for (i in seq_len(nrow(joint_selection))) {
    row <- joint_selection[i, , drop = FALSE]
    winner_dir <- row$source_dir[[1L]]
    local_winner <- file.path(winners_root, row$warmup_label[[1L]])
    dir.create(local_winner, recursive = TRUE, showWarnings = FALSE)
    sources <- file.path(winner_dir, table_names)
    destinations <- file.path(local_winner, table_names)
    Map(copy_input, sources, destinations)
    lineage <- figure5_rebuild_lineage(winner_dir)
    write_intermediate_tsv(
      lineage, file.path(local_winner, "invitro_lineage_summary.tsv")
    )
    figure5_rebuild_functional_curves(
      local_winner = local_winner,
      invitro_fit_result_path = invitro_fit_result_path
    )
    copied_sources <- c(copied_sources, sources)
    copied_destinations <- c(copied_destinations, destinations)
    copied_roles <- c(
      copied_roles,
      rep("selected joint-winner table/input", length(sources))
    )
  }

  tsne_dir <- file.path(
    JOINT_RESULT_ROOT,
    "joint_primary_clusters", "pooled_invivo_invitro",
    "full_data_in_vivo_clustring", "Tables"
  )
  tsne_names <- c(
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_full_coordinates.csv",
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
  )
  tsne_sources <- file.path(tsne_dir, tsne_names)
  tsne_destinations <- file.path(data_dir, tsne_names)
  Map(copy_input, tsne_sources, tsne_destinations)

  config_sources <- file.path(
    AUDIT_ROOT, "parameters",
    c("parameter_function_groups.tsv", "parameter_function_group_palette.tsv")
  )
  config_destinations <- file.path(data_dir, basename(config_sources))
  for (i in seq_along(config_sources)) {
    ok <- file.copy(
      config_sources[[i]], config_destinations[[i]], overwrite = TRUE
    )
    if (!ok) stop("Failed to stage Figure 5 config: ", config_sources[[i]])
  }

  all_joint_summary_files <- unlist(lapply(
    sort(list.dirs(JOINT_RESULT_ROOT, recursive = FALSE, full.names = TRUE)),
    function(pair_dir) {
      if (!grepl("^fit_joint_", basename(pair_dir))) return(character())
      seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
      seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
      file.path(seed_dirs, "fit_summary.tsv")
    }
  ), use.names = FALSE)
  require_files(
    all_joint_summary_files, "Figure 5 joint winner-selection input"
  )
  contract_sources <- c(
    copied_sources,
    tsne_sources,
    invitro_fit_result_path,
    all_joint_summary_files,
    file.path(
      MODEL_CODE_ROOT,
      c(
        "model/model_O2_supply_demand_MAP.R",
        "model/model_O2_supply_demand_MAP.cpp",
        "util/o2_supply_demand_map_shared.R",
        "util/o2_supply_demand_map_common_semantics.R"
      )
    )
  )
  contract <- data.frame(
    role = c(
      copied_roles,
      rep("pooled t-SNE coordinate input", length(tsne_sources)),
      paste(INVITRO_VISUALIZATION_SEED,
            "in-vitro lineage reconstruction state"),
      rep("joint winner-selection fit summary", length(all_joint_summary_files)),
      rep("external mechanistic model source", 4L)
    ),
    source = normalizePath(contract_sources, mustWork = TRUE),
    local_file = c(
      copied_destinations,
      normalizePath(tsne_destinations, mustWork = TRUE),
      NA_character_,
      rep(NA_character_, length(all_joint_summary_files)),
      rep(NA_character_, 4L)
    ),
    source_md5 = unname(tools::md5sum(contract_sources)),
    local_md5 = c(
      unname(tools::md5sum(copied_destinations)),
      unname(tools::md5sum(tsne_destinations)),
      NA_character_,
      rep(NA_character_, length(all_joint_summary_files)),
      rep(NA_character_, 4L)
    ),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure5", contract)
  invisible(selection)
}

if (sys.nframe() == 0L) data_Figure5()
