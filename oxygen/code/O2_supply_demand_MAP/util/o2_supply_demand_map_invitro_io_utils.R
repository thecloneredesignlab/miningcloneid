# Canonical in-vitro input, configuration, and fitted-parameter helpers.

ivt_locate_repo_root <- function(start = getwd()) {
  candidates <- unique(normalizePath(
    c(start, file.path(start, ".."), file.path(start, "..", "..")),
    mustWork = FALSE
  ))
  hit <- candidates[
    file.exists(file.path(candidates, "code", "O2_supply_demand_MAP", "model", "model_O2_supply_demand_MAP.R")) |
      file.exists(file.path(candidates, "code", "O2_supply_demand_MAP", "model_O2_supply_demand_MAP.R"))
  ][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop("Could not locate repository root from: ", start)
  }
  hit
}

ivt_source_map_model <- function(repo_root) {
  code_dir <- file.path(repo_root, "code", "O2_supply_demand_MAP")
  util_dir <- file.path(code_dir, "util")
  model_dir <- file.path(code_dir, "model")
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = code_dir)
  shared_path <- file.path(util_dir, "o2_supply_demand_map_shared.R")
  common_path <- file.path(util_dir, "o2_supply_demand_map_common_semantics.R")
  model_path <- file.path(model_dir, "model_O2_supply_demand_MAP.R")
  legacy_shared_path <- file.path(code_dir, "o2_supply_demand_map_shared.R")
  legacy_common_path <- file.path(code_dir, "o2_supply_demand_map_common_semantics.R")
  legacy_model_path <- file.path(code_dir, "model_O2_supply_demand_MAP.R")

  if (!file.exists(shared_path)) shared_path <- legacy_shared_path
  if (!file.exists(common_path)) common_path <- legacy_common_path
  if (!file.exists(model_path)) model_path <- legacy_model_path

  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(shared_path))
  source(shared_path, local = FALSE)
  source(common_path, local = FALSE)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
  source(model_path, local = FALSE)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = code_dir)
  invisible(code_dir)
}

ivt_canonical_flow_sample_name <- function(x) {
  x <- as.character(x)
  x <- sub("^Sample_", "", x)
  x <- sub("^SUM159", "SUM-159_NLS", x)
  x
}

ivt_normalize_flow_slot <- function(flow_entry) {
  if (is.null(flow_entry)) return(NULL)

  if (is.list(flow_entry) && !is.data.frame(flow_entry)) {
    if (!is.null(flow_entry$g0g1_ploidy_density)) {
      g0g1 <- flow_entry$g0g1_ploidy_density
      g0g1$density <- as.numeric(g0g1$density)
      g0g1$log_density <- as.numeric(g0g1$log_density)
      g0g1$ploidy <- as.numeric(g0g1$ploidy)
      g0g1$grid_index <- as.integer(g0g1$grid_index)
      flow_entry$g0g1_ploidy_density <- g0g1
    }
    return(flow_entry)
  }

  if (is.data.frame(flow_entry)) {
    return(list(
      type = "legacy_histogram",
      raw_histogram = flow_entry
    ))
  }

  list(
    type = "unknown",
    raw_value = flow_entry
  )
}

ivt_attach_g0g1_flow_data <- function(fit_data,
                                      repo_root,
                                      csv_path = file.path(repo_root, "data", "g0g1_ploidy_density_grid.csv")) {
  if (!file.exists(csv_path)) {
    return(fit_data)
  }

  flow_csv <- read.csv(csv_path, stringsAsFactors = FALSE)
  required_cols <- c("sample_name", "grid_index", "ploidy", "log_density")
  missing_cols <- setdiff(required_cols, names(flow_csv))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in G0/G1 flow CSV: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  flow_csv$sample_name_raw <- as.character(flow_csv$sample_name)
  flow_csv$sample_name <- ivt_canonical_flow_sample_name(flow_csv$sample_name_raw)
  flow_csv$grid_index <- as.integer(flow_csv$grid_index)
  flow_csv$ploidy <- as.numeric(flow_csv$ploidy)
  flow_csv$log_density <- as.numeric(flow_csv$log_density)
  flow_csv$density <- exp(flow_csv$log_density)

  flow_by_sample <- split(flow_csv, flow_csv$sample_name)
  for (sample_name in names(flow_by_sample)) {
    flow_df <- flow_by_sample[[sample_name]]
    flow_df <- flow_df[order(flow_df$grid_index, flow_df$ploidy), , drop = FALSE]
    keep_cols <- c("sample_name_raw", "sample_name", "grid_index", "ploidy", "log_density", "density")
    flow_by_sample[[sample_name]] <- flow_df[, keep_cols, drop = FALSE]
  }

  fit_ids <- names(fit_data)
  sample_hits <- lapply(fit_ids, function(pid) {
    hits <- names(flow_by_sample)[vapply(names(flow_by_sample), function(sample_name) {
      grepl(sample_name, pid, fixed = TRUE)
    }, logical(1))]
    if (length(hits) > 1L) {
      hits <- hits[order(nchar(hits), decreasing = TRUE)]
      hits <- hits[[1]]
    }
    if (length(hits) == 0L) NA_character_ else hits[[1]]
  })
  names(sample_hits) <- fit_ids

  for (pid in fit_ids) {
    fit_data[[pid]]$flow <- ivt_normalize_flow_slot(fit_data[[pid]]$flow)
    sample_name <- sample_hits[[pid]]
    if (is.na(sample_name) || !nzchar(sample_name)) next
    fit_data[[pid]]$flow$type <- "g0g1_ploidy_density"
    fit_data[[pid]]$flow$sample_name_canonical <- sample_name
    fit_data[[pid]]$flow$g0g1_ploidy_density <- flow_by_sample[[sample_name]]
  }

  fit_data
}

ivt_load_death_data <- function(path) {
  if (is.null(path) || !length(path) || is.na(path[[1]]) || !nzchar(trimws(path[[1]]))) {
    return(data.frame())
  }
  path <- normalizePath(path[[1]], mustWork = FALSE)
  if (!file.exists(path)) {
    stop("In vitro death-likelihood data not found: ", path)
  }
  death_data <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c(
    "observation_id", "cohort", "lineage_id", "scenario_id",
    "model_passage_id", "model_segment_id", "lineage_passage_index",
    "likelihood_observation_day", "include_in_current_endpoint_likelihood",
    "dead_count", "eligible_denominator", "observed_dead_fraction"
  )
  missing <- setdiff(required, names(death_data))
  if (length(missing)) {
    stop(
      "In vitro death-likelihood data is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  include <- tolower(trimws(as.character(
    death_data$include_in_current_endpoint_likelihood
  ))) %in% c("true", "t", "1", "yes", "y")
  death_data <- death_data[include, , drop = FALSE]
  numeric_columns <- c(
    "lineage_passage_index", "likelihood_observation_day", "dead_count",
    "eligible_denominator", "observed_dead_fraction"
  )
  death_data[numeric_columns] <- lapply(
    death_data[numeric_columns],
    function(x) suppressWarnings(as.numeric(x))
  )
  if (nrow(death_data) == 0L) {
    stop("In vitro death-likelihood data has no enabled endpoint rows: ", path)
  }
  if (any(!is.finite(death_data$dead_count)) ||
      any(!is.finite(death_data$eligible_denominator)) ||
      any(death_data$dead_count < 0) ||
      any(death_data$eligible_denominator <= 0) ||
      any(death_data$dead_count > death_data$eligible_denominator)) {
    stop("In vitro death-likelihood counts must satisfy 0 <= dead_count <= eligible_denominator.")
  }
  recomputed_fraction <- death_data$dead_count / death_data$eligible_denominator
  if (any(!is.finite(death_data$observed_dead_fraction)) ||
      any(abs(death_data$observed_dead_fraction - recomputed_fraction) > 1e-7)) {
    stop("observed_dead_fraction does not match dead_count / eligible_denominator.")
  }
  if (any(!is.finite(death_data$likelihood_observation_day)) ||
      any(death_data$likelihood_observation_day < 0)) {
    stop("likelihood_observation_day must contain finite non-negative values.")
  }
  identity_columns <- c(
    "observation_id", "cohort", "lineage_id", "scenario_id",
    "model_passage_id", "model_segment_id"
  )
  if (any(vapply(death_data[identity_columns], function(x) {
    any(is.na(x) | !nzchar(trimws(as.character(x))))
  }, logical(1)))) {
    stop("In vitro death-likelihood identity columns cannot contain missing or empty values.")
  }
  if (anyDuplicated(as.character(death_data$observation_id))) {
    stop("In vitro death-likelihood observation_id values must be unique.")
  }
  if (anyDuplicated(as.character(death_data$model_passage_id))) {
    stop("In vitro death-likelihood model_passage_id values must be unique.")
  }
  death_data$death_data_path <- path
  death_data
}

.ivt_growth_timecourse_from_endpoint <- function(fit_entry, passage_id) {
  duration <- suppressWarnings(as.numeric(fit_entry$passage_duration))
  initial_cells <- suppressWarnings(as.numeric(fit_entry$initial_cells))
  final_cells <- suppressWarnings(as.numeric(fit_entry$final_cells))
  if (length(duration) != 1L || !is.finite(duration) || duration <= 0 ||
      length(initial_cells) != 1L || !is.finite(initial_cells) || initial_cells <= 0 ||
      length(final_cells) != 1L || !is.finite(final_cells) || final_cells <= 0) {
    return(NULL)
  }
  data.frame(
    growth_observation_id = c(
      paste0(passage_id, "__day0"),
      paste0(passage_id, "__endpoint")
    ),
    observation_day = c(0, duration),
    observed_live_cells = c(initial_cells, final_cells),
    growth_data_source = "fit_data_endpoint_fallback",
    stringsAsFactors = FALSE
  )
}

.ivt_attach_growth_timecourse_data <- function(fit_data, repo_root) {
  metadata_path <- file.path(repo_root, "data", "metadata.csv")
  metadata <- NULL
  if (file.exists(metadata_path)) {
    metadata <- utils::read.csv(
      metadata_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    required <- c("passage_id", "num_date", "correctedCount")
    missing <- setdiff(required, names(metadata))
    if (length(missing) > 0L) {
      stop(
        "In vitro growth metadata is missing columns: ",
        paste(missing, collapse = ", ")
      )
    }
    metadata$passage_id <- trimws(as.character(metadata$passage_id))
    metadata$observation_day <- suppressWarnings(as.numeric(metadata$num_date))
    metadata$observed_live_cells <- suppressWarnings(as.numeric(metadata$correctedCount))
    if (any(is.na(metadata$passage_id) | !nzchar(metadata$passage_id)) ||
        any(!is.finite(metadata$observation_day)) ||
        any(metadata$observation_day < 0) ||
        any(!is.finite(metadata$observed_live_cells)) ||
        any(metadata$observed_live_cells <= 0)) {
      stop(
        "In vitro growth metadata requires non-empty passage_id values, ",
        "finite non-negative num_date values, and positive correctedCount values."
      )
    }
    if ("id" %in% names(metadata)) {
      metadata$growth_observation_id <- trimws(as.character(metadata$id))
    } else {
      metadata$growth_observation_id <- paste0(
        metadata$passage_id,
        "__day",
        format(metadata$observation_day, trim = TRUE, scientific = FALSE)
      )
    }
    if (any(is.na(metadata$growth_observation_id) |
            !nzchar(metadata$growth_observation_id)) ||
        anyDuplicated(metadata$growth_observation_id)) {
      stop("In vitro growth metadata observation identifiers must be non-empty and unique.")
    }
    metadata$growth_data_source <- normalizePath(metadata_path, mustWork = TRUE)
    metadata <- metadata[
      order(metadata$passage_id, metadata$observation_day, metadata$growth_observation_id),
      c(
        "passage_id", "growth_observation_id", "observation_day",
        "observed_live_cells", "growth_data_source"
      ),
      drop = FALSE
    ]
  }

  for (passage_id in names(fit_data)) {
    entry <- fit_data[[passage_id]]
    timecourse <- NULL
    if (!is.null(metadata)) {
      timecourse <- metadata[
        metadata$passage_id == passage_id,
        setdiff(names(metadata), "passage_id"),
        drop = FALSE
      ]
      if (nrow(timecourse) == 0L) timecourse <- NULL
    }
    if (!is.null(timecourse)) {
      day_zero <- which(abs(timecourse$observation_day) <= 1e-10)
      if (length(day_zero) != 1L) {
        stop("Growth metadata must contain exactly one day-0 count for passage ", passage_id, ".")
      }
      duration <- suppressWarnings(as.numeric(entry$passage_duration))
      initial_cells <- suppressWarnings(as.numeric(entry$initial_cells))
      final_cells <- suppressWarnings(as.numeric(entry$final_cells))
      last_idx <- which.max(timecourse$observation_day)
      nearly_equal <- function(x, y) {
        is.finite(x) && is.finite(y) &&
          abs(x - y) <= 1e-7 * max(1, abs(x), abs(y))
      }
      if (length(duration) != 1L ||
          !nearly_equal(timecourse$observation_day[[last_idx]], duration) ||
          length(initial_cells) != 1L ||
          !nearly_equal(timecourse$observed_live_cells[[day_zero]], initial_cells) ||
          length(final_cells) != 1L ||
          !nearly_equal(timecourse$observed_live_cells[[last_idx]], final_cells)) {
        stop(
          "Growth metadata day-0/endpoint values do not match fit_data for passage ",
          passage_id, "."
        )
      }
    } else {
      timecourse <- .ivt_growth_timecourse_from_endpoint(entry, passage_id)
    }
    entry$growth_timecourse <- timecourse
    fit_data[[passage_id]] <- entry
  }
  fit_data
}

ivt_load_fit_objects <- function(repo_root,
                                 fit_objects_dir = file.path(repo_root, "ploidyOxygen", "data", "fit_objects"),
                                 flow_csv_path = file.path(repo_root, "data", "g0g1_ploidy_density_grid.csv"),
                                 death_data_path = NULL) {
  base_dir <- normalizePath(fit_objects_dir, mustWork = FALSE)
  required_files <- c("fit_data.Rds", "jobs_2N.Rds", "jobs_4N.Rds")
  missing_files <- required_files[!file.exists(file.path(base_dir, required_files))]
  if (length(missing_files) > 0L) {
    stop(
      "Missing in vitro fit-object files under ", base_dir, ": ",
      paste(missing_files, collapse = ", ")
    )
  }
  fit_data <- readRDS(file.path(base_dir, "fit_data.Rds"))
  fit_data <- ivt_attach_g0g1_flow_data(
    fit_data = fit_data,
    repo_root = repo_root,
    csv_path = flow_csv_path
  )
  fit_data <- .ivt_attach_growth_timecourse_data(
    fit_data = fit_data,
    repo_root = repo_root
  )
  list(
    fit_data = fit_data,
    jobs_2N = readRDS(file.path(base_dir, "jobs_2N.Rds")),
    jobs_4N = readRDS(file.path(base_dir, "jobs_4N.Rds")),
    death_data = ivt_load_death_data(death_data_path)
  )
}

ivt_build_default_cfg <- function(repo_root,
                                  dt = 0.05,
                                  init_total_size = 1e6,
                                  o2_upper_bound = 21,
                                  fixed_oxygen = TRUE) {
  cfg <- list(
    parameter_table = file.path(repo_root, "data", "O2_supply_demand", "parameter_table_invitro.csv"),
    N_UNIT = 22L,
    N_MIN = 22L,
    N_MAX = 154L,
    start_with = "chr_number",
    DT = dt,
    init_total_size = init_total_size,
    o2_Nref = init_total_size,
    o2_min = 0,
    o2_burden_feedback = !isTRUE(fixed_oxygen),
    O2_growth = TRUE,
    Crowding = TRUE,
    crowding = "logistic",
    K = 1.5e7,
    dose_ref = 30,
    tx_mult_min = 0.05,
    min_pop = 1e-12,
    fit_treatment = FALSE,
    o2_cache_profile = FALSE,
    burden_log_eps = 1e-12,
    ploidy_O2_death = "ploidy_related",
    o2_S0_upper_bound = o2_upper_bound
  )
  normalize_sim_cfg_common(cfg, context = "viz")
}

ivt_parameter_table_path <- function(repo_root) {
  fname <- "parameter_table_invitro_buffering.csv"
  file.path(repo_root, "data", "O2_supply_demand", fname)
}

ivt_load_default_run_params <- function(cfg) {
  parameter_table <- read.csv(cfg$parameter_table, stringsAsFactors = FALSE)
  run_params <- setNames(
    lapply(parameter_table$init_value, as.numeric),
    parameter_table$param_symbol
  )
  normalize_run_params_common(run_params, cfg = cfg)
}

ivt_default_best_param_path <- function(repo_root) {
  file.path(repo_root, "data", "O2_supply_demand", "invitro_best_parameters.tsv")
}

ivt_load_run_params_from_row <- function(best_row, cfg) {
  run_params <- ivt_load_default_run_params(cfg)
  nat_cols <- grep("^nat__", names(best_row), value = TRUE)
  for (col in nat_cols) {
    nm <- sub("^nat__", "", col)
    run_params[[nm]] <- as.numeric(best_row[[col]][[1]])
  }
  if ("p_mis_base" %in% names(best_row)) {
    run_params$p_mis_base <- as.numeric(best_row$p_mis_base[[1]])
  }

  normalize_run_params_common(run_params, cfg = cfg)
}

ivt_resolve_sampled_objective_tsv <- function(repo_root,
                                              tsv_path = file.path(repo_root, "workflow", "sampled_objective_draws", "invitro_objective_array.tsv")) {
  if (length(tsv_path) != 1L || is.na(tsv_path) || !nzchar(trimws(tsv_path))) {
    stop("Sampled objective TSV choice must be a non-empty string.")
  }

  tsv_path <- trimws(tsv_path)
  out_dir <- file.path(repo_root, "workflow", "sampled_objective_draws")
  if (!identical(tolower(tsv_path), "latest") && dirname(tsv_path) == ".") {
    tsv_path <- file.path(out_dir, tsv_path)
  }
  if (file.exists(tsv_path)) return(tsv_path)

  if (!identical(tolower(tsv_path), "latest")) {
    out_dir <- dirname(tsv_path)
  }
  if (!dir.exists(out_dir)) {
    stop("Sampled objective directory not found: ", out_dir)
  }

  candidates <- list.files(
    out_dir,
    pattern = "^invitro_objective_array_[0-9_]+\\.tsv$",
    full.names = TRUE
  )
  if (length(candidates) == 0L) {
    stop("Sampled objective TSV not found: ", tsv_path)
  }

  info <- file.info(candidates)
  candidates[[order(info$mtime, decreasing = TRUE)[[1]]]]
}

ivt_load_committed_best_run_params <- function(cfg,
                                               repo_root,
                                               best_param_path = ivt_default_best_param_path(repo_root)) {
  if (!file.exists(best_param_path)) {
    stop("Committed best-parameter TSV not found: ", best_param_path)
  }

  tab <- read.delim(best_param_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(tab) == 0L) stop("Committed best-parameter TSV is empty: ", best_param_path)

  best <- tab[1, , drop = FALSE]
  list(
    run_params = ivt_load_run_params_from_row(best_row = best, cfg = cfg),
    best_row = best,
    source_tsv = best_param_path,
    source_type = "committed_best_parameters"
  )
}

ivt_load_best_sampled_run_params <- function(cfg,
                                             repo_root,
                                             tsv_path = file.path(repo_root, "workflow", "sampled_objective_draws", "invitro_objective_array.tsv"),
                                             best_param_path = ivt_default_best_param_path(repo_root)) {
  sampled_try <- tryCatch({
    tsv_use <- ivt_resolve_sampled_objective_tsv(repo_root = repo_root, tsv_path = tsv_path)
    tab <- read.delim(tsv_use, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(tab) == 0L) stop("Sampled objective TSV is empty: ", tsv_use)
    if (!"objective" %in% names(tab)) stop("Sampled objective TSV missing 'objective' column.")
    ord <- order(as.numeric(tab$objective))
    best <- tab[ord[[1]], , drop = FALSE]

    list(
      run_params = ivt_load_run_params_from_row(best_row = best, cfg = cfg),
      best_row = best,
      source_tsv = tsv_use,
      source_type = "sampled_objective_tsv"
    )
  }, error = function(e) {
    NULL
  })

  if (!is.null(sampled_try)) {
    return(sampled_try)
  }

  committed_try <- tryCatch(
    ivt_load_committed_best_run_params(cfg = cfg, repo_root = repo_root, best_param_path = best_param_path),
    error = function(e) e
  )
  if (!inherits(committed_try, "error")) {
    return(committed_try)
  }

  stop(
    "Could not load sampled-objective results or committed best-parameter fallback. ",
    "Sampled path tried: ", tsv_path,
    ". Fallback path tried: ", best_param_path,
    ". Fallback error: ", conditionMessage(committed_try)
  )
}
