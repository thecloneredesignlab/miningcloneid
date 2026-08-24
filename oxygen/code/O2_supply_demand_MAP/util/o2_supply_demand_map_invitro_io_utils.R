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

ivt_load_fit_objects <- local({
  death_filename <- "sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"
  death_expected_md5 <- "f346d7803596e15968b90b8213f15d95"

  resolve_death_path <- function(repo_root) {
    root_use <- normalizePath(repo_root, mustWork = FALSE)
    checkout_root <- if (dir.exists(file.path(root_use, "code", "O2_supply_demand_MAP"))) {
      dirname(root_use)
    } else {
      root_use
    }
    file.path(checkout_root, "data", "InVitroData", death_filename)
  }

  load_death_data <- function(repo_root) {
    death_path <- resolve_death_path(repo_root)
    if (!file.exists(death_path)) {
      return(list(
        data = data.frame(),
        path = normalizePath(death_path, mustWork = FALSE),
        md5 = NA_character_,
        n_file_rows = 0L
      ))
    }
    death_path <- normalizePath(death_path, mustWork = TRUE)
    death_md5 <- unname(as.character(tools::md5sum(death_path)))
    if (!identical(tolower(death_md5), death_expected_md5)) {
      stop(
        "In-vitro Death likelihood input MD5 mismatch for ", death_path,
        ". Expected ", death_expected_md5, ", found ", death_md5, "."
      )
    }

    raw <- utils::read.delim(
      death_path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character",
      quote = "",
      comment.char = ""
    )
    required_cols <- c(
      "observation_id", "cohort", "lineage_id", "model_passage_id",
      "include_in_current_endpoint_likelihood", "dead_count",
      "eligible_denominator", "observed_dead_fraction"
    )
    missing_cols <- setdiff(required_cols, names(raw))
    if (length(missing_cols) > 0L) {
      stop(
        "Missing required columns in in-vitro Death likelihood input: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    include_text <- toupper(trimws(raw$include_in_current_endpoint_likelihood))
    valid_include <- include_text %in% c("TRUE", "FALSE", "T", "F", "1", "0")
    if (any(!valid_include)) {
      stop("Invalid include_in_current_endpoint_likelihood value(s) in Death likelihood input.")
    }
    include <- include_text %in% c("TRUE", "T", "1")
    if (sum(include) != 90L) {
      stop(
        "Death likelihood input must contain exactly 90 included observations; found ",
        sum(include), "."
      )
    }
    death <- raw[include, , drop = FALSE]
    death$include_in_current_endpoint_likelihood <- TRUE

    for (field in c("observation_id", "cohort", "lineage_id", "model_passage_id")) {
      death[[field]] <- trimws(as.character(death[[field]]))
      if (any(is.na(death[[field]]) | !nzchar(death[[field]]))) {
        stop("Death likelihood input contains missing or empty ", field, " value(s).")
      }
    }
    if (anyDuplicated(death$observation_id)) {
      stop("Death likelihood observation_id values must be unique.")
    }
    if (anyDuplicated(death$model_passage_id)) {
      stop("Death likelihood model_passage_id values must be unique.")
    }

    numeric_field <- function(field) {
      value <- suppressWarnings(as.numeric(death[[field]]))
      if (any(!is.finite(value))) {
        stop("Death likelihood input contains non-finite ", field, " value(s).")
      }
      value
    }
    death$dead_count <- numeric_field("dead_count")
    death$eligible_denominator <- numeric_field("eligible_denominator")
    death$observed_dead_fraction <- numeric_field("observed_dead_fraction")
    if (any(death$dead_count < 0 | death$dead_count > death$eligible_denominator)) {
      stop("Death likelihood input violates 0 <= dead_count <= eligible_denominator.")
    }
    death$observed_live_count <- death$eligible_denominator - death$dead_count
    if (any(!is.finite(death$observed_live_count) | death$observed_live_count <= 0)) {
      stop("Death likelihood observed live counts must be strictly positive.")
    }

    fraction_text <- trimws(raw$observed_dead_fraction[include])
    fraction_digits <- ifelse(
      grepl(".", fraction_text, fixed = TRUE),
      nchar(sub("^[^.]*[.]", "", fraction_text)),
      0L
    )
    fraction_tolerance <- 0.5 * 10^(-fraction_digits) + 4 * .Machine$double.eps
    fraction_from_counts <- death$dead_count / death$eligible_denominator
    if (any(abs(death$observed_dead_fraction - fraction_from_counts) > fraction_tolerance)) {
      stop("observed_dead_fraction is inconsistent with dead_count/eligible_denominator at file precision.")
    }

    expected_groups <- c("2N-O1" = 23L, "2N-O2" = 23L, "4N-O1" = 22L, "4N-O2" = 22L)
    observed_groups <- table(paste(death$cohort, death$lineage_id, sep = "-"))
    if (!setequal(names(observed_groups), names(expected_groups)) ||
        any(as.integer(observed_groups[names(expected_groups)]) != unname(expected_groups))) {
      stop(
        "Death likelihood cohort/lineage counts must be ",
        paste(names(expected_groups), expected_groups, sep = "=", collapse = ", "), "."
      )
    }

    list(
      data = death,
      path = death_path,
      md5 = death_md5,
      n_file_rows = nrow(raw)
    )
  }

  function(repo_root,
           fit_objects_dir = file.path(repo_root, "ploidyOxygen", "data", "fit_objects"),
           flow_csv_path = file.path(repo_root, "data", "g0g1_ploidy_density_grid.csv")) {
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
    list(
      fit_data = fit_data,
      jobs_2N = readRDS(file.path(base_dir, "jobs_2N.Rds")),
      jobs_4N = readRDS(file.path(base_dir, "jobs_4N.Rds")),
      death_enabled = FALSE,
      death_data = data.frame(),
      death_data_path = NA_character_,
      death_data_md5 = NA_character_,
      death_data_n_file_rows = 0L
    )
  }
})

INVITRO_PASSAGE_IMPLEMENTATION <- "v2"

ivt_reject_removed_passage_mode <- function(value, source = "input") {
  if (!is.null(value) && length(value) > 0L) {
    stop(
      "passage_mode has been removed from ", source,
      "; in vitro always uses the fixed v2 passage implementation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
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
    K = 1e12,
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
