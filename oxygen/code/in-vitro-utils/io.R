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

ivt_load_fit_objects <- function(repo_root) {
  base_dir <- file.path(repo_root, "ploidyOxygen", "data", "fit_objects")
  list(
    fit_data = readRDS(file.path(base_dir, "fit_data.Rds")),
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
    o2_S0_upper_bound = o2_upper_bound
  )
  normalize_sim_cfg_common(cfg, context = "viz")
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
  if (file.exists(tsv_path)) return(tsv_path)

  out_dir <- dirname(tsv_path)
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
