ivt_locate_repo_root <- function(start = getwd()) {
  candidates <- unique(normalizePath(
    c(start, file.path(start, ".."), file.path(start, "..", "..")),
    mustWork = FALSE
  ))
  hit <- candidates[
    file.exists(file.path(candidates, "code", "O2_supply_demand_MAP", "model_O2_supply_demand_MAP.R"))
  ][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop("Could not locate repository root from: ", start)
  }
  hit
}

ivt_source_map_model <- function(repo_root) {
  code_dir <- file.path(repo_root, "code", "O2_supply_demand_MAP")
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = code_dir)
  old_wd <- getwd()
  setwd(code_dir)
  on.exit(setwd(old_wd), add = TRUE)
  source("o2_supply_demand_map_shared.R")
  source("o2_supply_demand_map_common_semantics.R")
  source("model_O2_supply_demand_MAP.R")
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

ivt_load_fit_objects <- function(repo_root) {
  base_dir <- file.path(repo_root, "ploidyOxygen", "data", "fit_objects")
  fit_data <- readRDS(file.path(base_dir, "fit_data.Rds"))
  fit_data <- ivt_attach_g0g1_flow_data(fit_data = fit_data, repo_root = repo_root)
  list(
    fit_data = fit_data,
    jobs_2N = readRDS(file.path(base_dir, "jobs_2N.Rds")),
    jobs_4N = readRDS(file.path(base_dir, "jobs_4N.Rds"))
  )
}

ivt_build_default_cfg <- function(repo_root,
                                  dt = 0.1,
                                  init_total_size = 1e6,
                                  o2_upper_bound = 21,
                                  fixed_oxygen = TRUE) {
  cfg <- list(
    parameter_table = file.path(repo_root, "data", "O2_supply_demand", "parameter_table_invitro.csv"),
    N_UNIT = 22L,
    N_MIN = 22L,
    N_MAX = 154L,
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
    ploidy_O2_death = "uniform",
    misseg_loss_survival = "nullisomy",
    o2_S0_upper_bound = o2_upper_bound
  )
  normalize_sim_cfg_common(cfg, context = "viz")
}

ivt_parameter_table_for_loss_mode <- function(repo_root, loss_mode = "nullisomy") {
  loss_mode <- canonical_misseg_loss_survival_mode(loss_mode, "nullisomy")
  fname <- if (identical(loss_mode, "buffering")) {
    "parameter_table_invitro_buffering.csv"
  } else {
    "parameter_table_invitro.csv"
  }
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
  if ("misseg_loss_survival" %in% names(best_row)) {
    run_params$misseg_loss_survival <- as.character(best_row$misseg_loss_survival[[1]])
  }
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
