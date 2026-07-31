#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

figure5_select_winners <- function() {
  pair_dirs <- sort(list.dirs(
    JOINT_RESULT_ROOT, recursive = FALSE, full.names = TRUE
  ))
  pair_dirs <- pair_dirs[grepl("^fit_joint_", basename(pair_dirs))]
  if (length(pair_dirs) != 6L) {
    stop("Figure 5 requires six fit_joint_* result directories.")
  }
  rows <- lapply(pair_dirs, function(pair_dir) {
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
    invivo_seed <- sub("^.*_vi_(seed[0-9]+)_.*$", "\\1", pair_id)
    invitro_seed <- sub("^.*_vt_(seed[0-9]+)$", "\\1", pair_id)
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
  selected <- selected[order(
    match(
      sub("^.*_(C[0-9]+Sc[0-9]+)_.*$", "\\1", selected$warmup_label),
      c("C01Sc01", "C01Sc02", "C02Sc01", "C02Sc02", "C03Sc01", "C03Sc02")
    )
  ), , drop = FALSE]
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
      stringsAsFactors = FALSE
    )
  }
  rbind(
    selected,
    add_separate(
      "separate_invivo_global_best", "seed25", "objective",
      INVIVO_RESULT_ROOT
    ),
    add_separate(
      "separate_invitro_global_best", "seed10", "objective_total",
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

  fit <- readRDS(file.path(INVITRO_RESULT_ROOT, "seed10", "fit_result.rds"))
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
    seed10_fit_result_path
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

  seed10_fit <- readRDS(seed10_fit_result_path)
  invitro_cfg <- normalize_sim_cfg_common(
    seed10_fit$cfg, context = "viz"
  )
  invitro_params <- seed10_fit$best_params
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

  copied_sources <- character()
  copied_destinations <- character()
  table_names <- c(
    "invitro_growth_loglik.tsv",
    "invivo_burden_fit.tsv",
    "invivo_terminal_ploidy_fit.tsv",
    "joint_soft_coupling.tsv",
    "invitro_effective_params.tsv",
    "best_params.tsv",
    "fit_config.rds"
  )
  seed10_fit_result_path <- file.path(
    INVITRO_RESULT_ROOT, "seed10", "fit_result.rds"
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
      seed10_fit_result_path = seed10_fit_result_path
    )
    copied_sources <- c(copied_sources, sources)
    copied_destinations <- c(copied_destinations, destinations)
  }

  tsne_dir <- file.path(
    JOINT_RESULT_ROOT,
    "landscape_subcluster", "pooled_invivo_invitro",
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
    seed10_fit_result_path,
    all_joint_summary_files
  )
  contract <- data.frame(
    role = c(
      rep("selected joint-winner table/input", length(copied_sources)),
      rep("pooled t-SNE coordinate input", length(tsne_sources)),
      "seed10 in-vitro lineage reconstruction state",
      rep("joint winner-selection fit summary", length(all_joint_summary_files))
    ),
    source = normalizePath(contract_sources, mustWork = TRUE),
    local_file = c(
      copied_destinations,
      normalizePath(tsne_destinations, mustWork = TRUE),
      NA_character_,
      rep(NA_character_, length(all_joint_summary_files))
    ),
    source_md5 = unname(tools::md5sum(contract_sources)),
    local_md5 = c(
      unname(tools::md5sum(copied_destinations)),
      unname(tools::md5sum(tsne_destinations)),
      NA_character_,
      rep(NA_character_, length(all_joint_summary_files))
    ),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure5", contract)
  invisible(selection)
}

if (sys.nframe() == 0L) data_Figure5()
