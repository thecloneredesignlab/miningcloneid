#!/usr/bin/env Rscript

# Exploratory standalone nullisomy-only fitter.
# This script reuses the active O2_supply_demand_MAP backend, but it is not yet
# integrated into the main optimizer/config path used by the primary invivo
# fitting workflow.

suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

repo_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
workflow_root <- file.path(repo_root, "oxygen", "code", "O2_supply_demand_MAP")
model_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
shared_path <- file.path(workflow_root, "util", "o2_supply_demand_map_shared.R")
if (!file.exists(model_path)) stop("Missing active model backend: ", model_path)
if (!file.exists(shared_path)) stop("Missing active shared utils: ", shared_path)
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
source(shared_path, local = TRUE)
source(model_path, local = TRUE)

parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0L) return(out)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
    a <- substring(a, 3L)
    if (grepl("=", a, fixed = TRUE)) {
      kv <- strsplit(a, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- paste(kv[-1], collapse = "=")
      i <- i + 1L
    } else {
      key <- a
      if (i == length(argv)) stop("Missing value for --", key)
      out[[key]] <- argv[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) default else v
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) default else v
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

default_data_path <- function() {
  file.path(
    repo_root,
    "data",
    "DrugResponseData_PloidyJumps",
    "DrugResponseData",
    "SNU668_doxorubicin_response_and_karyotype_ploidy.json"
  )
}

default_data_paths <- function() {
  c(default_data_path())
}

dataset_ids_from_cfg <- function(cfg) {
  ids <- cfg$dataset_ids
  if (is.null(ids) || length(ids) == 0L) "dataset1" else as.character(ids)
}

dataset_id_from_json_path <- function(json_path) {
  base <- basename(json_path)
  sub("\\.json$", "", base, ignore.case = TRUE)
}

timestamp_string <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

make_output_dir <- function(out_dir = NULL) {
  path <- if (!is.null(out_dir) && nzchar(out_dir)) {
    out_dir
  } else {
    file.path(getwd(), "doxo_outputs", paste0("doxo_nullisomy_", timestamp_string()))
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

load_doxo_data <- function(json_path) {
  if (!file.exists(json_path)) stop("JSON not found: ", json_path)
  dat <- jsonlite::fromJSON(json_path, simplifyDataFrame = TRUE)
  need <- c(
    "drug_response_summary_by_dose",
    "drug_response_replicate_level",
    "cell_ploidy",
    "sample_ploidy_summary"
  )
  miss <- setdiff(need, names(dat))
  if (length(miss) > 0L) stop("JSON missing fields: ", paste(miss, collapse = ", "))
  dat
}

load_doxo_datasets <- function(json_paths, dataset_ids = NULL) {
  json_paths <- as.character(json_paths)
  if (length(json_paths) == 0L) stop("No JSON paths provided.")
  if (is.null(dataset_ids)) {
    dataset_ids <- vapply(json_paths, dataset_id_from_json_path, character(1))
  }
  dataset_ids <- as.character(dataset_ids)
  if (length(dataset_ids) != length(json_paths)) stop("dataset_ids length must match json_paths length.")
  if (anyDuplicated(dataset_ids)) stop("Dataset ids must be unique.")
  out <- setNames(vector("list", length(json_paths)), dataset_ids)
  for (i in seq_along(json_paths)) {
    out[[dataset_ids[[i]]]] <- list(
      dataset_id = dataset_ids[[i]],
      json_path = json_paths[[i]],
      data_obj = load_doxo_data(json_paths[[i]])
    )
  }
  out
}

canonical_response_data_mode <- function(x) {
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("auto", "prefer_smoothed", "smoothed_if_available", "prefer_hill_smoothed")) return("prefer_hill_smoothed")
  if (s %in% c("summary", "summary_by_dose", "raw_summary")) return("summary_by_dose")
  if (s %in% c("smoothed", "hill_smoothed", "hill_smoothed_by_dose")) return("hill_smoothed_by_dose")
  stop("Unsupported response_data_mode: ", x)
}

normalize_response_dose_table <- function(tab, value_col, sd_col = NULL, source_name = "dose table") {
  dose_col <- if ("dose_uM" %in% names(tab)) {
    "dose_uM"
  } else if ("dose_nM" %in% names(tab)) {
    "dose_nM"
  } else {
    stop(source_name, " missing dose column: expected dose_uM or dose_nM.")
  }
  if (!("sample" %in% names(tab))) {
    stop(source_name, " missing required column: sample.")
  }
  need <- c("sample", dose_col, value_col)
  if (!is.null(sd_col)) need <- c(need, sd_col)
  miss <- setdiff(need, names(tab))
  if (length(miss) > 0L) {
    stop(source_name, " missing fields: ", paste(miss, collapse = ", "))
  }
  out <- tab[, need, drop = FALSE]
  names(out)[names(out) == dose_col] <- "dose_uM"
  names(out)[names(out) == value_col] <- "mean_normalized_to_0uM"
  if (!is.null(sd_col)) {
    names(out)[names(out) == sd_col] <- "sd_normalized_to_0uM"
  }
  if (identical(dose_col, "dose_nM")) {
    out$dose_uM <- as.numeric(out$dose_uM) / 1000
  } else {
    out$dose_uM <- as.numeric(out$dose_uM)
  }
  out
}

response_table_from_data <- function(data_obj, cfg) {
  mode <- canonical_response_data_mode(cfg$response_data_mode)
  dr_summary <- data_obj$drug_response_summary_by_dose
  has_smoothed <- "drug_response_hill_smoothed_by_dose" %in% names(data_obj)
  if (identical(mode, "summary_by_dose") || (identical(mode, "prefer_hill_smoothed") && !has_smoothed)) {
    return(normalize_response_dose_table(
      dr_summary,
      value_col = if ("mean_normalized_to_0uM" %in% names(dr_summary)) "mean_normalized_to_0uM" else "mean_normalized_to_0nM",
      sd_col = if ("sd_normalized_to_0uM" %in% names(dr_summary)) "sd_normalized_to_0uM" else "sd_normalized_to_0nM",
      source_name = "drug_response_summary_by_dose"
    ))
  }

  if (!has_smoothed) {
    stop("response_data_mode=hill_smoothed_by_dose requested but JSON lacks drug_response_hill_smoothed_by_dose.")
  }
  dr_smooth <- data_obj$drug_response_hill_smoothed_by_dose
  out <- normalize_response_dose_table(
    dr_smooth,
    value_col = if ("predicted_normalized_to_0uM" %in% names(dr_smooth)) "predicted_normalized_to_0uM" else "predicted_normalized_to_0nM",
    sd_col = NULL,
    source_name = "drug_response_hill_smoothed_by_dose"
  )
  sd_tab <- normalize_response_dose_table(
    dr_summary,
    value_col = if ("mean_normalized_to_0uM" %in% names(dr_summary)) "mean_normalized_to_0uM" else "mean_normalized_to_0nM",
    sd_col = if ("sd_normalized_to_0uM" %in% names(dr_summary)) "sd_normalized_to_0uM" else "sd_normalized_to_0nM",
    source_name = "drug_response_summary_by_dose"
  )[, c("sample", "dose_uM", "sd_normalized_to_0uM"), drop = FALSE]
  out <- merge(out, sd_tab, by = c("sample", "dose_uM"), all.x = TRUE, sort = FALSE)
  out <- out[, c("sample", "dose_uM", "mean_normalized_to_0uM", "sd_normalized_to_0uM"), drop = FALSE]
  out
}

extract_sample_chromosome_totals <- function(cell_ploidy, sample_name) {
  if (!("chromosome_total" %in% names(cell_ploidy))) {
    stop("cell_ploidy is missing required column 'chromosome_total'.")
  }
  if (!("sample" %in% names(cell_ploidy))) {
    stop("cell_ploidy is missing required column 'sample'.")
  }
  sub <- cell_ploidy[cell_ploidy$sample == sample_name, , drop = FALSE]
  if (nrow(sub) == 0L) stop("No cell_ploidy rows found for sample: ", sample_name)
  vals <- suppressWarnings(as.numeric(sub$chromosome_total))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    stop("No finite chromosome_total values found for sample: ", sample_name)
  }
  vals
}

validate_doxo_data <- function(data_obj, cfg = default_cfg()) {
  dr <- response_table_from_data(data_obj, cfg)
  cp <- data_obj$cell_ploidy

  if (!("sample" %in% names(cp))) stop("cell_ploidy missing field: sample")
  if (!("chromosome_total" %in% names(cp))) stop("cell_ploidy missing field: chromosome_total")

  if (any(!is.finite(suppressWarnings(as.numeric(dr$dose_uM))))) {
    stop("Response table contains non-finite dose_uM values.")
  }
  if (any(!is.finite(suppressWarnings(as.numeric(dr$mean_normalized_to_0uM))))) {
    stop("Response table contains non-finite mean_normalized_to_0uM values.")
  }
  sd_num <- suppressWarnings(as.numeric(dr$sd_normalized_to_0uM))
  if (any(!is.finite(sd_num))) {
    stop("Response table contains non-finite sd_normalized_to_0uM values.")
  }
  if (any(sd_num < 0)) {
    stop("Response table contains negative sd_normalized_to_0uM values.")
  }

  dup_key <- paste(dr$sample, dr$dose_uM, sep = "\r")
  if (anyDuplicated(dup_key)) {
    stop("Response table contains duplicated (sample, dose_uM) rows.")
  }

  samples <- unique(as.character(dr$sample))
  samples <- samples[!is.na(samples) & nzchar(samples)]
  if (length(samples) == 0L) stop("No samples found in response table.")

  samples_ploidy <- unique(as.character(cp$sample))
  samples_ploidy <- samples_ploidy[!is.na(samples_ploidy) & nzchar(samples_ploidy)]
  missing_ploidy <- setdiff(samples, samples_ploidy)
  if (length(missing_ploidy) > 0L) {
    stop("Missing cell_ploidy entries for samples: ", paste(missing_ploidy, collapse = ", "))
  }

  for (sample_name in samples) {
    extract_sample_chromosome_totals(cp, sample_name)
    sub <- dr[dr$sample == sample_name, , drop = FALSE]
    zero_idx <- which(as.numeric(sub$dose_uM) == 0)
    if (length(zero_idx) != 1L) {
      stop("Sample ", sample_name, " must have exactly one dose_uM == 0 row; found ", length(zero_idx), ".")
    }
  }

  invisible(TRUE)
}

default_cfg <- function() {
  list(
    mode = "doxorubicin_nullisomy",
    N_UNIT = 22L,
    N_MIN = 22L,
    N_MAX = 154L,
    DT = 0.25,
    start_with = "chr_number",
    chr_lengths_bp = default_chr_lengths_bp_1to22(),
    Crowding = FALSE,
    crowding = "logistic",
    K = 1e12,
    min_pop = 1e-12,
    fit_treatment = FALSE,
    dose_ref = 1.0,
    tx_mult_min = 1.0,
    o2_burden_feedback = FALSE,
    O2_growth = FALSE,
    o2_cache_bin_pct = 0.01,
    o2_cache_hysteresis_pct = 0.005,
    o2_cache_profile = FALSE,
    p_mis_base = 1e-5,
    gamma_loss = 1.0,
    p_wgd = 1e-8,
    rho_2N = 42332.020717,
    o2_fixed = 5.0,
    o2_S0 = 5.0,
    o2_min = 5.0,
    kappa_O = 1.0,
    tau_O2 = 1.0,
    o2_Nref = 1e6,
    assay_N0 = 3000,
    k_o = 1.0,
    alpha_o2 = 0.5,
    gamma_growth = 2.0,
    O2_crit = 1.0,
    n_O = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    k_clear = 1e-3,
    burden_log_eps = 1e-12,
    sigma_log_response_floor = 0.05,
    sigma_seed = 0.0,
    observable_mode = "cell_count",
    nullisomy_hidden_copy_mode = "dirichlet_multinomial",
    nullisomy_dirichlet_mc_samples = 10000L,
    nullisomy_dirichlet_seed = 12345L,
    fixed_params = list(),
    fixed_params_by_variant = list(),
    response_data_mode = "prefer_hill_smoothed",
    ec50_dataset_mode = "free",
    hafner_ec50_min_denom = 1e-3,
    hafner_ec50_penalty = 1e12,
    p_mis_cap = 0.95,
    p_mis_implausible_threshold = 0.3,
    N_2N_ref = 44,
    N_ref_2N_state = 44,
    N_ref_high_state = 95,
    high_ploidy_threshold = 66,
    assay_days_fixed = 4.0,
    optim_maxit = 12,
    fit_models_default = c("shared_per_copy", "continuous_ploidy_amplified", "categorical_high_ploidy"),
    seeds_default = c(1L, 2L, 3L, 4L, 5L),
    dataset_ids = "dataset1"
  )
}

canonical_observable_mode <- function(x) {
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("cell_count", "count", "cells")) return("cell_count")
  if (s %in% c("volume_weighted_burden", "volume", "burden")) return("volume_weighted_burden")
  stop("Unsupported observable_mode: ", x)
}

canonical_ec50_dataset_mode <- function(x) {
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("free", "dataset_specific", "dataset-specific")) return("free")
  if (s %in% c("shared", "global")) return("shared")
  if (s %in% c("hafner", "hafner_lambda", "lambda_hafner")) return("hafner_lambda")
  stop("Unsupported ec50_dataset_mode: ", x)
}

hill_pmis <- function(dose_uM, p_mis_base, p_mis_doxo_max, EC50, hill_n) {
  dose <- pmax(as.numeric(dose_uM), 0)
  p_base <- max(as.numeric(p_mis_base), 0)
  p_amp <- max(as.numeric(p_mis_doxo_max), 0)
  ec50 <- max(as.numeric(EC50), 1e-12)
  n <- max(as.numeric(hill_n), 1e-12)
  frac <- dose^n / (ec50^n + dose^n)
  clip01(p_base + p_amp * frac)
}

expected_mis_copies <- function(N, p) {
  as.numeric(N) * as.numeric(p)
}

sample_initial_state <- function(cell_ploidy, sample_name, cfg) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  vals <- extract_sample_chromosome_totals(cell_ploidy, sample_name)
  vals <- pmin(pmax(round(vals), cfg$N_MIN), cfg$N_MAX)
  tab <- table(factor(vals, levels = grid))
  prob <- as.numeric(tab)
  prob <- prob / sum(prob)
  names(prob) <- as.character(grid)
  prob
}

sample_initial_counts <- function(cell_ploidy, sample_name, cfg, assay_N0 = NULL) {
  if (is.null(assay_N0)) assay_N0 <- cfg$assay_N0
  sample_initial_state(cell_ploidy, sample_name, cfg) * as.numeric(assay_N0)
}

sample_names_from_data <- function(data_obj, cfg = default_cfg()) {
  dr <- response_table_from_data(data_obj, cfg)
  samples_resp <- unique(as.character(dr$sample))
  samples_resp <- samples_resp[!is.na(samples_resp) & nzchar(samples_resp)]
  if (length(samples_resp) == 0L) stop("No samples found in response table.")

  samples_ploidy <- unique(as.character(data_obj$cell_ploidy$sample))
  samples_ploidy <- samples_ploidy[!is.na(samples_ploidy) & nzchar(samples_ploidy)]
  missing_ploidy <- setdiff(samples_resp, samples_ploidy)
  if (length(missing_ploidy) > 0L) {
    stop("Missing cell_ploidy entries for samples: ", paste(missing_ploidy, collapse = ", "))
  }
  samples_resp
}

pick_reference_samples <- function(data_obj, cfg) {
  samples <- sample_names_from_data(data_obj, cfg)
  stats <- lapply(samples, function(sample_name) {
    init_prob <- sample_initial_state(data_obj$cell_ploidy, sample_name, cfg)
    grid <- as.numeric(names(init_prob))
    mean_chr <- sum(grid * init_prob)
    data.frame(sample = sample_name, mean_chromosome_total = mean_chr, stringsAsFactors = FALSE)
  })
  stats <- bind_rows(stats)
  ord <- order(stats$mean_chromosome_total, stats$sample)
  list(
    low = stats$sample[[ord[[1L]]]],
    high = stats$sample[[ord[[length(ord)]]]]
  )
}

make_run_params <- function(theta, cfg, dataset_id = NULL) {
  lambda_eff <- max(lambda_eff_for_dataset(theta, cfg, dataset_id), 1e-6)
  list(
    lam_min = lambda_eff,
    lam_max = lambda_eff,
    k_o = cfg$k_o,
    o2_min = cfg$o2_min,
    p_mis_base = cfg$p_mis_base,
    p_misseg = NULL,
    k_o_mis = 1.0,
    gamma_loss = cfg$gamma_loss,
    p_wgd = cfg$p_wgd,
    o2_S0 = cfg$o2_S0,
    kappa_O = cfg$kappa_O,
    eta_o2 = 1.0,
    rho_2N = cfg$rho_2N,
    alpha_o2 = cfg$alpha_o2,
    gamma_growth = cfg$gamma_growth,
    mu_hp = 0.0,
    gamma_mu = cfg$gamma_mu,
    O2_crit = cfg$O2_crit,
    n_O = cfg$n_O,
    k_clear = cfg$k_clear,
    tau_O2 = cfg$tau_O2,
    alpha = 0.0,
    gamma = 1.0,
    ploidy_O2_death = "uniform",
    start_with = "chr_number",
    boundary = "drop"
  )
}

lambda_eff_for_dataset <- function(theta, cfg, dataset_id = NULL) {
  dataset_ids <- dataset_ids_from_cfg(cfg)
  dataset_id <- if (is.null(dataset_id)) dataset_ids[[1L]] else as.character(dataset_id)
  nm <- paste0("lambda_eff__", dataset_id)
  if (!is.null(theta[[nm]])) return(as.numeric(theta[[nm]]))
  if (!is.null(theta$lambda_eff_by_dataset[[dataset_id]])) return(as.numeric(theta$lambda_eff_by_dataset[[dataset_id]]))
  if (!is.null(theta$lambda_eff)) return(as.numeric(theta$lambda_eff))
  stop("Missing lambda_eff for dataset_id: ", dataset_id)
}

lambda_eff_by_dataset_from_theta <- function(theta, cfg) {
  dataset_ids <- dataset_ids_from_cfg(cfg)
  setNames(
    as.list(vapply(dataset_ids, function(id) {
      nm <- paste0("lambda_eff__", id)
      if (!is.null(theta[[nm]])) return(as.numeric(theta[[nm]]))
      if (!is.null(theta$lambda_eff_by_dataset[[id]])) return(as.numeric(theta$lambda_eff_by_dataset[[id]]))
      if (!is.null(theta$lambda_eff) && length(dataset_ids) == 1L) return(as.numeric(theta$lambda_eff))
      stop("Missing lambda_eff for dataset_id: ", id)
    }, numeric(1))),
    dataset_ids
  )
}

hafner_ec50_terms <- function(theta, cfg) {
  dataset_ids <- dataset_ids_from_cfg(cfg)
  lambda_vals <- unlist(lambda_eff_by_dataset_from_theta(theta, cfg), use.names = FALSE)
  names(lambda_vals) <- dataset_ids
  Td <- log(2) / pmax(lambda_vals, 1e-12)
  denom <- cfg$assay_days_fixed * theta$S_M - Td
  ratio <- Td / denom
  EC50 <- rep(NA_real_, length(ratio))
  valid <- is.finite(ratio) & ratio > 0 & is.finite(theta$S50) & is.finite(theta$hill_n) & theta$S50 > 0 & theta$hill_n > 0
  EC50[valid] <- theta$S50 * ratio[valid]^(1 / theta$hill_n)
  list(dataset_ids = dataset_ids, lambda_eff = lambda_vals, doubling_time = Td, denom = denom, EC50 = EC50)
}

hafner_ec50_is_valid <- function(theta, cfg) {
  if (!identical(canonical_ec50_dataset_mode(cfg$ec50_dataset_mode), "hafner_lambda")) return(TRUE)
  terms <- hafner_ec50_terms(theta, cfg)
  all(is.finite(terms$EC50)) &&
    all(is.finite(terms$denom)) &&
    all(terms$denom > cfg$hafner_ec50_min_denom) &&
    is.finite(theta$S50) &&
    is.finite(theta$S_M) &&
    theta$S50 > 0 &&
    theta$S_M > 0
}

ec50_by_dataset_from_theta <- function(theta, cfg, fail_invalid = TRUE) {
  dataset_ids <- dataset_ids_from_cfg(cfg)
  mode <- canonical_ec50_dataset_mode(cfg$ec50_dataset_mode)
  if (identical(mode, "hafner_lambda")) {
    terms <- hafner_ec50_terms(theta, cfg)
    valid <- all(is.finite(terms$EC50)) &&
      all(is.finite(terms$denom)) &&
      all(terms$denom > cfg$hafner_ec50_min_denom)
    if (!valid && fail_invalid) {
      stop(
        "Invalid Hafner EC50 correction: require assay_days_fixed * S_M - doubling_time > ",
        cfg$hafner_ec50_min_denom,
        " for every dataset."
      )
    }
    out <- as.list(terms$EC50)
    names(out) <- terms$dataset_ids
    return(out)
  }
  if (identical(mode, "shared")) {
    if (is.null(theta$EC50)) stop("Missing shared EC50.")
    return(setNames(as.list(rep(as.numeric(theta$EC50), length(dataset_ids))), dataset_ids))
  }
  if (length(dataset_ids) == 1L) {
    if (!is.null(theta$EC50)) return(setNames(list(as.numeric(theta$EC50)), dataset_ids))
  }
  setNames(
    as.list(vapply(dataset_ids, function(id) {
      nm <- paste0("EC50__", id)
      if (!is.null(theta[[nm]])) return(as.numeric(theta[[nm]]))
      if (!is.null(theta$EC50_by_dataset[[id]])) return(as.numeric(theta$EC50_by_dataset[[id]]))
      if (!is.null(theta$EC50) && length(dataset_ids) == 1L) return(as.numeric(theta$EC50))
      stop("Missing EC50 for dataset_id: ", id)
    }, numeric(1))),
    dataset_ids
  )
}

refresh_dataset_parameter_maps <- function(theta, cfg, fail_invalid_ec50 = TRUE) {
  theta$lambda_eff_by_dataset <- lambda_eff_by_dataset_from_theta(theta, cfg)
  theta$EC50_by_dataset <- ec50_by_dataset_from_theta(theta, cfg, fail_invalid = fail_invalid_ec50)
  dataset_ids <- dataset_ids_from_cfg(cfg)
  if (length(dataset_ids) == 1L) {
    theta$lambda_eff <- theta$lambda_eff_by_dataset[[dataset_ids[[1L]]]]
    if (is.null(theta$EC50)) theta$EC50 <- theta$EC50_by_dataset[[dataset_ids[[1L]]]]
  }
  theta
}

ensure_valid_hafner_initial_theta <- function(theta, cfg) {
  if (!identical(canonical_ec50_dataset_mode(cfg$ec50_dataset_mode), "hafner_lambda")) return(theta)
  dataset_ids <- dataset_ids_from_cfg(cfg)
  lambda_vals <- unlist(lambda_eff_by_dataset_from_theta(theta, cfg), use.names = FALSE)
  Td <- log(2) / pmax(lambda_vals, 1e-12)
  min_SM <- (max(Td, na.rm = TRUE) + max(cfg$hafner_ec50_min_denom, 1e-3)) / cfg$assay_days_fixed
  if (is.finite(min_SM) && (!is.finite(theta$S_M) || theta$S_M <= min_SM)) {
    theta$S_M <- min(max(min_SM * 1.25, min_SM + 0.05), 5.0)
  }
  refresh_dataset_parameter_maps(theta, cfg, fail_invalid_ec50 = FALSE)
}

ec50_for_dataset <- function(theta, cfg, dataset_id = NULL) {
  dataset_ids <- dataset_ids_from_cfg(cfg)
  dataset_id <- if (is.null(dataset_id)) dataset_ids[[1L]] else as.character(dataset_id)
  vals <- ec50_by_dataset_from_theta(theta, cfg, fail_invalid = TRUE)
  if (!is.null(vals[[dataset_id]])) return(as.numeric(vals[[dataset_id]]))
  stop("Missing EC50 for dataset_id: ", dataset_id)
}

canonical_model_variant <- function(x) {
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("shared", "shared_per_copy", "null")) return("shared_per_copy")
  if (s %in% c("continuous", "continuous_ploidy_amplified", "alpha")) return("continuous_ploidy_amplified")
  if (s %in% c("categorical", "categorical_high_ploidy", "a4n")) return("categorical_high_ploidy")
  stop("Unsupported model variant: ", x)
}

model_variant_label <- function(x) {
  x <- canonical_model_variant(x)
  switch(
    x,
    shared_per_copy = "Shared per-copy",
    continuous_ploidy_amplified = "Continuous ploidy-amplified",
    categorical_high_ploidy = "Categorical high-ploidy amplified"
  )
}

parameter_spec <- function(variant, cfg) {
  variant <- canonical_model_variant(variant)
  dataset_ids <- dataset_ids_from_cfg(cfg)
  ec50_mode <- canonical_ec50_dataset_mode(cfg$ec50_dataset_mode)
  ec50_names <- switch(
    ec50_mode,
    free = if (length(dataset_ids) == 1L) "EC50" else paste0("EC50__", dataset_ids),
    shared = "EC50",
    hafner_lambda = c("S50", "S_M")
  )
  ec50_defaults <- if (identical(ec50_mode, "hafner_lambda")) c(0.2, 1.0) else rep(0.2, length(ec50_names))
  ec50_lowers <- if (identical(ec50_mode, "hafner_lambda")) c(1e-4, 0.05) else rep(1e-4, length(ec50_names))
  ec50_uppers <- if (identical(ec50_mode, "hafner_lambda")) c(100, 5) else rep(100, length(ec50_names))
  lambda_eff_lower <- if (identical(ec50_mode, "hafner_lambda")) 0.1 else 1e-3
  rows <- data.frame(
    name = c("p_mis_doxo_max", ec50_names, "hill_n", paste0("lambda_eff__", dataset_ids), "response_top", "nullisomy_dirichlet_alpha"),
    transform = c("logit_headroom", rep("log10", length(ec50_names)), "log10", rep("log10", length(dataset_ids)), "log10", "log10"),
    lower = c(NA_real_, ec50_lowers, 0.2, rep(lambda_eff_lower, length(dataset_ids)), 1e-3, 1e-4),
    upper = c(NA_real_, ec50_uppers, 8, rep(2, length(dataset_ids)), 10, 100),
    default = c(0.05, ec50_defaults, 1.5, rep(0.2, length(dataset_ids)), 1.0, 10.0),
    stringsAsFactors = FALSE
  )
  if (identical(variant, "continuous_ploidy_amplified")) {
    rows <- bind_rows(
      rows,
      data.frame(
        name = c("ec50_alpha", "hilln_alpha", "pmis_alpha"),
        transform = c("identity", "identity", "identity"),
        lower = c(-2.0, -2.0, -2.0),
        upper = c(2.0, 2.0, 2.0),
        default = c(0.0, 0.0, 0.0),
        stringsAsFactors = FALSE
      )
    )
  } else if (identical(variant, "categorical_high_ploidy")) {
    rows <- bind_rows(
      rows,
      data.frame(
        name = c("a_4N_ec50", "a_4N_hilln", "a_4N_pmis"),
        transform = c("log10", "log10", "log10"),
        lower = c(0.25, 0.25, 0.25),
        upper = c(8.0, 8.0, 8.0),
        default = c(1.0, 1.0, 1.0),
        stringsAsFactors = FALSE
      )
    )
  }
  rows
}

parameter_names_for_variant <- function(variant, cfg) {
  parameter_spec(variant, cfg)$name
}

parse_fixed_params <- function(x) {
  if (is.null(x) || !nzchar(trimws(as.character(x)))) return(list())
  parts <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) return(list())
  keys <- character(length(parts))
  vals <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    kv <- strsplit(parts[[i]], "=", fixed = TRUE)[[1]]
    if (length(kv) < 2L) stop("Invalid fixed_params entry: ", parts[[i]])
    key <- trimws(kv[[1]])
    val_txt <- trimws(paste(kv[-1], collapse = "="))
    if (!nzchar(key) || !nzchar(val_txt)) stop("Invalid fixed_params entry: ", parts[[i]])
    val_num <- suppressWarnings(as.numeric(val_txt))
    if (!is.finite(val_num)) stop("Fixed parameter must be numeric: ", key, "=", val_txt)
    keys[[i]] <- key
    vals[[i]] <- val_num
  }
  if (anyDuplicated(keys)) stop("Duplicate fixed parameter names in --fixed_params.")
  names(vals) <- keys
  vals
}

parse_variant_scoped_fixed_params <- function(args) {
  out <- list()
  key_map <- list(
    fixed_params_shared_per_copy = "shared_per_copy",
    fixed_params_shared = "shared_per_copy",
    fixed_params_continuous_ploidy_amplified = "continuous_ploidy_amplified",
    fixed_params_continuous = "continuous_ploidy_amplified",
    fixed_params_categorical_high_ploidy = "categorical_high_ploidy",
    fixed_params_categorical = "categorical_high_ploidy"
  )
  for (arg_name in names(key_map)) {
    if (!is.null(args[[arg_name]])) {
      out[[key_map[[arg_name]]]] <- parse_fixed_params(args[[arg_name]])
    }
  }
  out
}

natural_to_opt <- function(name, value, cfg, transform = NULL) {
  if (is.null(transform)) {
    spec <- parameter_spec("categorical_high_ploidy", cfg)
    transform <- spec$transform[match(name, spec$name)]
  }
  if (identical(transform, "log10")) return(log10(value))
  if (identical(transform, "identity")) return(as.numeric(value))
  if (identical(transform, "logit_headroom")) {
    p_headroom <- max(cfg$p_mis_cap - cfg$p_mis_base - 1e-8, 1e-8)
    frac <- clip01(as.numeric(value) / p_headroom)
    return(qlogis(pmin(pmax(frac, 1e-8), 1 - 1e-8)))
  }
  stop("Unsupported transform: ", transform)
}

opt_to_natural <- function(name, value, cfg, transform = NULL) {
  if (is.null(transform)) {
    spec <- parameter_spec("categorical_high_ploidy", cfg)
    transform <- spec$transform[match(name, spec$name)]
  }
  if (identical(transform, "log10")) return(10^value)
  if (identical(transform, "identity")) return(as.numeric(value))
  if (identical(transform, "logit_headroom")) {
    p_headroom <- max(cfg$p_mis_cap - cfg$p_mis_base - 1e-8, 1e-8)
    return(p_headroom * plogis(value))
  }
  stop("Unsupported transform: ", transform)
}

validate_fixed_params <- function(fixed_params, variant, cfg) {
  variant <- canonical_model_variant(variant)
  if (length(fixed_params) == 0L) return(invisible(TRUE))
  spec <- parameter_spec(variant, cfg)
  bad_names <- setdiff(names(fixed_params), spec$name)
  if (length(bad_names) > 0L) {
    stop("Fixed parameters not valid for ", variant, ": ", paste(bad_names, collapse = ", "))
  }
  for (nm in names(fixed_params)) {
    val <- as.numeric(fixed_params[[nm]])
    if (!is.finite(val)) stop("Fixed parameter must be finite numeric: ", nm)
    row <- spec[spec$name == nm, , drop = FALSE]
    lo <- row$lower[[1]]
    hi <- row$upper[[1]]
    if (is.finite(lo) && val < lo) stop("Fixed parameter ", nm, " below lower bound.")
    if (is.finite(hi) && val > hi) stop("Fixed parameter ", nm, " above upper bound.")
  }
  invisible(TRUE)
}

fixed_params_for_variant <- function(cfg, variant) {
  variant <- canonical_model_variant(variant)
  global_fixed <- cfg$fixed_params
  scoped_fixed <- cfg$fixed_params_by_variant[[variant]]
  if (is.null(global_fixed)) global_fixed <- list()
  if (is.null(scoped_fixed)) scoped_fixed <- list()
  merged <- global_fixed
  if (length(scoped_fixed) > 0L) {
    for (nm in names(scoped_fixed)) merged[[nm]] <- scoped_fixed[[nm]]
  }
  validate_fixed_params(merged, variant, cfg)
  merged
}

fixed_param_names <- function(cfg, variant) {
  names(fixed_params_for_variant(cfg, variant))
}

default_theta_from_spec <- function(variant, cfg) {
  spec <- parameter_spec(variant, cfg)
  out <- as.list(spec$default)
  names(out) <- spec$name
  refresh_dataset_parameter_maps(out, cfg, fail_invalid_ec50 = FALSE)
}

variant_param_count <- function(variant, cfg) {
  variant <- canonical_model_variant(variant)
  length(parameter_names_for_variant(variant, cfg)) - length(fixed_param_names(cfg, variant))
}

theta_from_par <- function(par, cfg, variant) {
  variant <- canonical_model_variant(variant)
  spec <- parameter_spec(variant, cfg)
  out <- default_theta_from_spec(variant, cfg)
  fixed <- fixed_params_for_variant(cfg, variant)
  if (length(fixed) > 0L) {
    for (nm in names(fixed)) out[[nm]] <- as.numeric(fixed[[nm]])
  }
  free_spec <- spec[!(spec$name %in% fixed_param_names(cfg, variant)), , drop = FALSE]
  stopifnot(length(par) == nrow(free_spec))
  for (i in seq_len(nrow(free_spec))) {
    out[[free_spec$name[[i]]]] <- opt_to_natural(free_spec$name[[i]], par[[i]], cfg, free_spec$transform[[i]])
  }
  refresh_dataset_parameter_maps(out, cfg, fail_invalid_ec50 = FALSE)
}

par_from_theta <- function(theta, cfg, variant) {
  variant <- canonical_model_variant(variant)
  spec <- parameter_spec(variant, cfg)
  free_spec <- spec[!(spec$name %in% fixed_param_names(cfg, variant)), , drop = FALSE]
  vapply(seq_len(nrow(free_spec)), function(i) {
    nm <- free_spec$name[[i]]
    natural_to_opt(nm, theta[[nm]], cfg, free_spec$transform[[i]])
  }, numeric(1))
}

default_theta <- function(cfg, variant) {
  variant <- canonical_model_variant(variant)
  out <- default_theta_from_spec(variant, cfg)
  fixed <- fixed_params_for_variant(cfg, variant)
  if (length(fixed) > 0L) {
    for (nm in names(fixed)) {
      out[[nm]] <- as.numeric(fixed[[nm]])
    }
  }
  refresh_dataset_parameter_maps(out, cfg, fail_invalid_ec50 = FALSE)
}

fit_bounds <- function(variant, cfg) {
  variant <- canonical_model_variant(variant)
  spec <- parameter_spec(variant, cfg)
  free_spec <- spec[!(spec$name %in% fixed_param_names(cfg, variant)), , drop = FALSE]
  lower <- vapply(seq_len(nrow(free_spec)), function(i) {
    natural_to_opt(free_spec$name[[i]], free_spec$lower[[i]], cfg, free_spec$transform[[i]])
  }, numeric(1))
  upper <- vapply(seq_len(nrow(free_spec)), function(i) {
    natural_to_opt(free_spec$name[[i]], free_spec$upper[[i]], cfg, free_spec$transform[[i]])
  }, numeric(1))
  list(lower = lower, upper = upper)
}

ec50_by_state <- function(N, theta, cfg, variant, dataset_id = NULL) {
  variant <- canonical_model_variant(variant)
  N_use <- pmax(as.numeric(N), 1)
  EC50 <- ec50_for_dataset(theta, cfg, dataset_id)
  if (identical(variant, "shared_per_copy")) return(rep(EC50, length(N_use)))
  if (identical(variant, "continuous_ploidy_amplified")) {
    return(EC50 * (N_use / cfg$N_2N_ref)^(-theta$ec50_alpha))
  }
  ifelse(N_use >= cfg$high_ploidy_threshold, EC50 / theta$a_4N_ec50, EC50)
}

hilln_by_state <- function(N, theta, cfg, variant) {
  variant <- canonical_model_variant(variant)
  N_use <- pmax(as.numeric(N), 1)
  if (identical(variant, "shared_per_copy")) return(rep(theta$hill_n, length(N_use)))
  if (identical(variant, "continuous_ploidy_amplified")) {
    return(theta$hill_n * (N_use / cfg$N_2N_ref)^(theta$hilln_alpha))
  }
  ifelse(N_use >= cfg$high_ploidy_threshold, theta$hill_n * theta$a_4N_hilln, theta$hill_n)
}

pmis_multiplier_by_state <- function(N, theta, cfg, variant) {
  variant <- canonical_model_variant(variant)
  N_use <- pmax(as.numeric(N), 1)
  if (identical(variant, "shared_per_copy")) return(rep(1.0, length(N_use)))
  if (identical(variant, "continuous_ploidy_amplified")) {
    return((N_use / cfg$N_2N_ref)^(theta$pmis_alpha))
  }
  ifelse(N_use >= cfg$high_ploidy_threshold, theta$a_4N_pmis, 1.0)
}

effective_pmis_by_state <- function(dose_uM, N, theta, cfg, variant, dataset_id = NULL) {
  variant <- canonical_model_variant(variant)
  N_use <- as.numeric(N)
  ec50_use <- ec50_by_state(N_use, theta, cfg, variant, dataset_id = dataset_id)
  hilln_use <- hilln_by_state(N_use, theta, cfg, variant)
  pmis_mult <- pmis_multiplier_by_state(N_use, theta, cfg, variant)
  out <- hill_pmis(
    dose_uM = dose_uM,
    p_mis_base = cfg$p_mis_base,
    p_mis_doxo_max = theta$p_mis_doxo_max,
    EC50 = ec50_use,
    hill_n = hilln_use
  ) * pmis_mult
  pmin(pmax(out, 0), cfg$p_mis_cap)
}

build_state_generator <- function(theta, dose_uM, cfg, variant, dataset_id = NULL) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  p_vec <- effective_pmis_by_state(dose_uM, grid, theta, cfg, variant, dataset_id = dataset_id)
  lambda_eff <- lambda_eff_for_dataset(theta, cfg, dataset_id)
  G <- .build_G_with_WGD(
    N0min = cfg$N_MIN,
    N0max = cfg$N_MAX,
    lambda0_vec = rep(lambda_eff, length(grid)),
    p0_vec = p_vec,
    wgd_prob_vec = rep(cfg$p_wgd, length(grid)),
    boundary = "drop",
    eps_tail = 1e-8,
    gamma_loss = cfg$gamma_loss,
    N_unit = cfg$N_UNIT,
    nullisomy_hidden_copy_mode = cfg$nullisomy_hidden_copy_mode,
    nullisomy_dirichlet_alpha = theta$nullisomy_dirichlet_alpha,
    nullisomy_dirichlet_mc_samples = cfg$nullisomy_dirichlet_mc_samples,
    nullisomy_dirichlet_seed = cfg$nullisomy_dirichlet_seed
  )
  list(G = G, p_vec = p_vec, grid = grid)
}

simulate_one_doxo <- function(init_state_prob, dose_uM, theta, cfg, variant, dataset_id = NULL) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  init_total <- cfg$assay_N0
  init_state <- as.numeric(init_state_prob) * init_total
  run_params <- make_run_params(theta, cfg, dataset_id)
  obs_steps <- as.integer(round(cfg$assay_days_fixed / cfg$DT))
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid, run_params = run_params, cfg = cfg))
  built <- build_state_generator(theta, dose_uM, cfg, variant, dataset_id = dataset_id)
  x_live <- as.numeric(init_state)
  for (step in seq_len(obs_steps)) {
    x_live <- step_dt(built$G, x_live, cfg$DT, steps = 1L, normalize = FALSE)
    x_live[!is.finite(x_live) | x_live < 0] <- 0
  }
  live_observable <- switch(
    canonical_observable_mode(cfg$observable_mode),
    cell_count = sum(x_live),
    volume_weighted_burden = sum(x_live * vol_by_N)
  )
  p_ref_2N <- effective_pmis_by_state(dose_uM, cfg$N_ref_2N_state, theta, cfg, variant, dataset_id = dataset_id)[[1]]
  p_ref_high <- effective_pmis_by_state(dose_uM, cfg$N_ref_high_state, theta, cfg, variant, dataset_id = dataset_id)[[1]]
  list(
    p_const = hill_pmis(dose_uM, cfg$p_mis_base, theta$p_mis_doxo_max, ec50_for_dataset(theta, cfg, dataset_id), theta$hill_n),
    p_ref_2N = p_ref_2N,
    p_ref_high = p_ref_high,
    expected_mis_2N = expected_mis_copies(cfg$N_ref_2N_state, p_ref_2N),
    expected_mis_high = expected_mis_copies(cfg$N_ref_high_state, p_ref_high),
    live_observable = live_observable,
    frac_N = x_live / max(sum(x_live), cfg$burden_log_eps)
  )
}

predict_dataset_single <- function(theta, data_obj, cfg, variant, dataset_id = NULL) {
  variant <- canonical_model_variant(variant)
  dataset_id <- if (is.null(dataset_id)) dataset_ids_from_cfg(cfg)[[1L]] else as.character(dataset_id)
  dr <- response_table_from_data(data_obj, cfg)
  samples <- sample_names_from_data(data_obj, cfg)
  init_states <- setNames(
    lapply(samples, function(sample_name) sample_initial_state(data_obj$cell_ploidy, sample_name, cfg)),
    samples
  )
  rows <- vector("list", length = 0L)
  for (sample_name in samples) {
    sub <- dr[dr$sample == sample_name, , drop = FALSE]
    sub <- sub[order(sub$dose_uM), , drop = FALSE]
    sims <- lapply(sub$dose_uM, function(dose) simulate_one_doxo(init_states[[sample_name]], dose, theta, cfg, variant, dataset_id = dataset_id))
    zero_idx <- which(as.numeric(sub$dose_uM) == 0)
    if (length(zero_idx) != 1L) {
      stop("Sample ", sample_name, " must have exactly one dose_uM == 0 row for normalization; found ", length(zero_idx), ".")
    }
    base_live <- sims[[zero_idx[[1L]]]]$live_observable
    base_live <- max(base_live, cfg$burden_log_eps)
    for (i in seq_len(nrow(sub))) {
        rows[[length(rows) + 1L]] <- data.frame(
          dataset_id = dataset_id,
          model = variant,
          model_label = model_variant_label(variant),
          sample = sample_name,
          dose_uM = sub$dose_uM[[i]],
          observed = sub$mean_normalized_to_0uM[[i]],
          observed_sd = sub$sd_normalized_to_0uM[[i]],
          predicted = theta$response_top * (sims[[i]]$live_observable / base_live),
          p_const = sims[[i]]$p_const,
          p_eff_2N = sims[[i]]$p_ref_2N,
          p_eff_high = sims[[i]]$p_ref_high,
          expected_mis_2N = sims[[i]]$expected_mis_2N,
        expected_mis_high = sims[[i]]$expected_mis_high,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

is_dataset_bundle <- function(x) {
  is.list(x) && length(x) > 0L && is.list(x[[1L]]) && all(c("dataset_id", "data_obj") %in% names(x[[1L]]))
}

predict_dataset <- function(theta, data_obj, cfg, variant) {
  if (is_dataset_bundle(data_obj)) {
    return(bind_rows(lapply(data_obj, function(ds) predict_dataset_single(theta, ds$data_obj, cfg, variant, dataset_id = ds$dataset_id))))
  }
  predict_dataset_single(theta, data_obj, cfg, variant, dataset_id = dataset_ids_from_cfg(cfg)[[1L]])
}

response_sigma <- function(obs, obs_sd, sigma_floor = 0.05) {
  obs_use <- pmax(as.numeric(obs), 1e-8)
  sd_use <- as.numeric(obs_sd)
  sigma <- ifelse(is.finite(sd_use) & sd_use > 0, sd_use / obs_use, sigma_floor)
  pmax(sigma, sigma_floor)
}

objective_from_table <- function(pred_tab, cfg, sigma_seed = cfg$sigma_seed) {
  obs <- pmax(as.numeric(pred_tab$observed), cfg$burden_log_eps)
  pred <- pmax(as.numeric(pred_tab$predicted), cfg$burden_log_eps)
  sigma_obs <- response_sigma(obs, pred_tab$observed_sd, cfg$sigma_log_response_floor)
  sigma <- sqrt(sigma_obs^2 + 2 * as.numeric(sigma_seed)^2)
  z <- (log(pred) - log(obs)) / sigma
  0.5 * sum(z^2 + log(2 * pi * sigma^2), na.rm = TRUE)
}

fit_one_seed <- function(seed, data_obj, cfg, variant) {
  variant <- canonical_model_variant(variant)
  base <- default_theta(cfg, variant)
  jittered <- jitter_theta(base, cfg, variant, seed)
  init_par <- par_from_theta(jittered, cfg, variant)
  bounds <- fit_bounds(variant, cfg)
  best_valid_par <- NULL
  best_valid_value <- Inf
  fn <- function(par) {
    theta <- theta_from_par(par, cfg, variant)
    if (!hafner_ec50_is_valid(theta, cfg)) return(cfg$hafner_ec50_penalty)
    pred <- tryCatch(
      predict_dataset(theta, data_obj, cfg, variant),
      error = function(e) NULL
    )
    if (is.null(pred)) return(cfg$hafner_ec50_penalty)
    value <- objective_from_table(pred, cfg, sigma_seed = 0.0)
    if (is.finite(value) && value < best_valid_value) {
      best_valid_value <<- value
      best_valid_par <<- par
    }
    value
  }
  opt <- optim(
    par = init_par,
    fn = fn,
    method = "L-BFGS-B",
    lower = bounds$lower,
    upper = bounds$upper,
    control = list(maxit = cfg$optim_maxit)
  )
  final_par <- opt$par
  final_value <- opt$value
  theta_hat <- theta_from_par(final_par, cfg, variant)
  if (!hafner_ec50_is_valid(theta_hat, cfg)) {
    if (is.null(best_valid_par)) {
      stop("Optimizer did not evaluate a valid Hafner EC50 parameter set; try different seeds or bounds.")
    }
    final_par <- best_valid_par
    final_value <- best_valid_value
    theta_hat <- theta_from_par(final_par, cfg, variant)
  }
  pred_hat <- predict_dataset(theta_hat, data_obj, cfg, variant)
  list(
    model = variant,
    seed = seed,
    objective = final_value,
    convergence = opt$convergence,
    message = opt$message,
    theta = theta_hat,
    predictions = pred_hat
  )
}

diagnostic_curves <- function(cfg, out_dir) {
  doses <- c(0, exp(seq(log(1e-4), log(10), length.out = 200L)))
  doses_plot <- pmax(doses, 1e-4)
  p_base <- cfg$p_mis_base
  pmaxs <- c(0.01, 0.05, 0.15)
  ec50s <- c(0.05, 0.2, 1.0)
  hill_ns <- c(1, 2)
  curve_rows <- vector("list", length = 0L)
  for (pmax in pmaxs) {
    for (ec50 in ec50s) {
      for (hn in hill_ns) {
        curve_rows[[length(curve_rows) + 1L]] <- data.frame(
          dose_uM = doses,
          dose_uM_plot = doses_plot,
          p_mis = hill_pmis(doses, p_base, pmax, ec50, hn),
          spec = sprintf("pmax=%.3f EC50=%.3f n=%.1f", pmax, ec50, hn),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  curve_tab <- bind_rows(curve_rows)
  g1 <- ggplot(curve_tab, aes(dose_uM_plot, p_mis, color = spec)) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    labs(title = "p_mis vs doxorubicin dose", x = "Dose (uM)", y = "Per-copy p_mis")
  ggsave(file.path(out_dir, "diag_pmis_vs_dose.png"), g1, width = 8, height = 5, dpi = 150)

  states <- c(cfg$N_ref_2N_state, 53, 90, cfg$N_ref_high_state)
  expected_tab <- bind_rows(lapply(states, function(N) {
    data.frame(
      dose_uM = doses,
      dose_uM_plot = doses_plot,
      expected_mis = expected_mis_copies(N, hill_pmis(doses, p_base, 0.05, 0.2, 2)),
      N = as.factor(N),
      stringsAsFactors = FALSE
    )
  }))
  g2 <- ggplot(expected_tab, aes(dose_uM_plot, expected_mis, color = N)) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    labs(title = "Expected missegregated copies per division", x = "Dose (uM)", y = "N * p_mis(C)")
  ggsave(file.path(out_dir, "diag_expected_mis_copies.png"), g2, width = 8, height = 5, dpi = 150)

  p_grid <- exp(seq(log(1e-5), log(0.5), length.out = 120L))
  nulli_tab <- bind_rows(lapply(states, function(N) {
    dead_prob <- vapply(p_grid, function(p) {
      attr(
        .pr_delta_vec(
          N, p,
          gamma_loss = cfg$gamma_loss,
          N_unit = cfg$N_UNIT,
          nullisomy_hidden_copy_mode = cfg$nullisomy_hidden_copy_mode,
          nullisomy_dirichlet_alpha = 10,
          nullisomy_dirichlet_mc_samples = cfg$nullisomy_dirichlet_mc_samples,
          nullisomy_dirichlet_seed = cfg$nullisomy_dirichlet_seed
        ),
        "mass_dropped"
      )
    }, numeric(1))
    data.frame(p_mis = p_grid, dead_prob = dead_prob, N = as.factor(N), stringsAsFactors = FALSE)
  }))
  g3 <- ggplot(nulli_tab, aes(p_mis, dead_prob, color = N)) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    labs(title = "Nullisomy/nonviable daughter mass vs p_mis", x = "Per-copy p_mis", y = "Dropped daughter mass")
  ggsave(file.path(out_dir, "diag_nullisomy_vs_pmis.png"), g3, width = 8, height = 5, dpi = 150)

  invisible(list(curve_tab = curve_tab, expected_tab = expected_tab, nulli_tab = nulli_tab))
}

initial_scale_diagnostics <- function(data_obj, cfg, out_dir) {
  if (is_dataset_bundle(data_obj)) {
    rows <- bind_rows(lapply(data_obj, function(ds) {
      tmp_dir <- tempfile("init_scale_")
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
      tab <- initial_scale_diagnostics(ds$data_obj, cfg, out_dir = tmp_dir)
      tab$dataset_id <- ds$dataset_id
      tab
    }))
    utils::write.table(rows, file.path(out_dir, "initial_scale_diagnostic.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    g <- ggplot(rows, aes(pmax(dose_uM, 1e-4), abs_diff, color = sample)) +
      geom_point(size = 1.8) +
      geom_line(linewidth = 0.9) +
      scale_x_log10() +
      facet_wrap(~ dataset_id) +
      labs(title = "Normalized response sensitivity to initialization scale", x = "Dose (uM)", y = "|response(3000) - response(1e6)|")
    ggsave(file.path(out_dir, "initial_scale_diagnostic.png"), g, width = 10, height = 6, dpi = 150)
    return(invisible(rows))
  }
  theta <- default_theta(cfg, "shared_per_copy")
  dr <- response_table_from_data(data_obj, cfg)
  samples <- unique(as.character(dr$sample))
  rows <- vector("list", length = 0L)
  scales <- c(cfg$assay_N0, cfg$o2_Nref)
  scale_labels <- c("assay_N0", "legacy_o2_Nref")
  for (sample_name in samples) {
    init_prob <- sample_initial_state(data_obj$cell_ploidy, sample_name, cfg)
    sub <- dr[dr$sample == sample_name, , drop = FALSE]
    sub <- sub[order(sub$dose_uM), , drop = FALSE]
    for (k in seq_along(scales)) {
      cfg_k <- cfg
      cfg_k$assay_N0 <- scales[[k]]
      sims <- lapply(sub$dose_uM, function(dose) simulate_one_doxo(init_prob, dose, theta, cfg_k, "shared_per_copy"))
      zero_idx <- which(as.numeric(sub$dose_uM) == 0)
      if (length(zero_idx) != 1L) {
        stop("Sample ", sample_name, " must have exactly one dose_uM == 0 row for initial-scale diagnostics; found ", length(zero_idx), ".")
      }
      base_live <- max(sims[[zero_idx[[1L]]]]$live_observable, cfg$burden_log_eps)
      for (i in seq_len(nrow(sub))) {
        rows[[length(rows) + 1L]] <- data.frame(
          sample = sample_name,
          dose_uM = sub$dose_uM[[i]],
          init_scale = scale_labels[[k]],
          init_cells = scales[[k]],
          predicted = sims[[i]]$live_observable / base_live,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  tab <- bind_rows(rows) %>%
    select(sample, dose_uM, init_scale, predicted) %>%
    tidyr::pivot_wider(names_from = init_scale, values_from = predicted) %>%
    mutate(abs_diff = abs(assay_N0 - legacy_o2_Nref))
  utils::write.table(tab, file.path(out_dir, "initial_scale_diagnostic.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  g <- ggplot(tab, aes(pmax(dose_uM, 1e-4), abs_diff, color = sample)) +
    geom_point(size = 1.8) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    labs(title = "Normalized response sensitivity to initialization scale", x = "Dose (uM)", y = "|response(3000) - response(1e6)|")
  ggsave(file.path(out_dir, "initial_scale_diagnostic.png"), g, width = 8, height = 5, dpi = 150)
  invisible(tab)
}

run_tests <- function(data_obj, cfg, out_dir) {
  dr_used <- response_table_from_data(data_obj, cfg)
  if ("drug_response_hill_smoothed_by_dose" %in% names(data_obj)) {
    sm <- normalize_response_dose_table(
      data_obj$drug_response_hill_smoothed_by_dose,
      value_col = if ("predicted_normalized_to_0uM" %in% names(data_obj$drug_response_hill_smoothed_by_dose)) "predicted_normalized_to_0uM" else "predicted_normalized_to_0nM",
      sd_col = NULL,
      source_name = "drug_response_hill_smoothed_by_dose"
    )
    key <- c("sample", "dose_uM")
    merged_sm <- merge(
      dr_used,
      sm[, c("sample", "dose_uM", "mean_normalized_to_0uM"), drop = FALSE],
      by = key,
      sort = FALSE
    )
    stopifnot(nrow(merged_sm) == nrow(dr_used))
    stopifnot(all(abs(merged_sm$mean_normalized_to_0uM.x - merged_sm$mean_normalized_to_0uM.y) < 1e-12))
  }

  data_nm <- data_obj
  data_nm$drug_response_summary_by_dose <- data_nm$drug_response_summary_by_dose %>%
    mutate(
      dose_nM = dose_uM * 1000,
      mean_normalized_to_0nM = mean_normalized_to_0uM,
      sd_normalized_to_0nM = sd_normalized_to_0uM
    ) %>%
    select(-dose_uM, -mean_normalized_to_0uM, -sd_normalized_to_0uM)
  dr_nm <- response_table_from_data(data_nm, cfg)
  stopifnot(nrow(dr_nm) == nrow(dr_used))
  stopifnot(max(abs(dr_nm$dose_uM - dr_used$dose_uM)) < 1e-12)
  stopifnot(max(abs(dr_nm$mean_normalized_to_0uM - dr_used$mean_normalized_to_0uM)) < 1e-12)
  stopifnot(max(abs(dr_nm$sd_normalized_to_0uM - dr_used$sd_normalized_to_0uM)) < 1e-12)

  doses <- c(0, 0.01, 0.1, 1.0, 10.0)
  p <- hill_pmis(doses, cfg$p_mis_base, 0.05, 0.2, 2.0)
  stopifnot(all(diff(p) >= -1e-12))
  stopifnot(abs(p[[1]] - cfg$p_mis_base) < 1e-12)
  stopifnot(all(p <= cfg$p_mis_base + 0.05 + 1e-12))

  p_a <- hill_pmis(c(0.01, 1.0), cfg$p_mis_base, 0.05, 0.2, 2.0)
  stopifnot(expected_mis_copies(cfg$N_ref_2N_state, p_a[[2]]) > expected_mis_copies(cfg$N_ref_2N_state, p_a[[1]]))
  stopifnot(expected_mis_copies(cfg$N_ref_high_state, 0.02) > expected_mis_copies(cfg$N_ref_2N_state, 0.02))

  surv_low <- .loss_survival_nullisomy(
    cfg$N_ref_2N_state, m_loss = 1, gamma_loss = cfg$gamma_loss, N_unit = cfg$N_UNIT,
    nullisomy_hidden_copy_mode = cfg$nullisomy_hidden_copy_mode,
    nullisomy_dirichlet_alpha = 10,
    nullisomy_dirichlet_mc_samples = cfg$nullisomy_dirichlet_mc_samples,
    nullisomy_dirichlet_seed = cfg$nullisomy_dirichlet_seed
  )
  surv_high <- .loss_survival_nullisomy(
    95, m_loss = 1, gamma_loss = cfg$gamma_loss, N_unit = cfg$N_UNIT,
    nullisomy_hidden_copy_mode = cfg$nullisomy_hidden_copy_mode,
    nullisomy_dirichlet_alpha = 10,
    nullisomy_dirichlet_mc_samples = cfg$nullisomy_dirichlet_mc_samples,
    nullisomy_dirichlet_seed = cfg$nullisomy_dirichlet_seed
  )
  stopifnot(surv_high >= surv_low)

  samples <- sample_names_from_data(data_obj, cfg)
  init_counts <- lapply(samples, function(s) sample_initial_counts(data_obj$cell_ploidy, s, cfg))
  stopifnot(all(vapply(init_counts, function(x) abs(sum(x) - cfg$assay_N0) < 1e-8, logical(1))))

  theta <- default_theta(cfg, "shared_per_copy")
  refs <- pick_reference_samples(data_obj, cfg)
  init_high <- sample_initial_state(data_obj$cell_ploidy, refs$high, cfg)
  sim_shared <- simulate_one_doxo(init_high, dose_uM = 1.0, theta = theta, cfg = cfg, variant = "shared_per_copy")
  theta_ec50 <- modifyList(default_theta(cfg, "continuous_ploidy_amplified"), list(ec50_alpha = 0.5, hilln_alpha = 0.0, pmis_alpha = 0.0))
  stopifnot(ec50_by_state(cfg$N_ref_high_state, theta_ec50, cfg, "continuous_ploidy_amplified") <
              ec50_by_state(cfg$N_ref_2N_state, theta_ec50, cfg, "continuous_ploidy_amplified"))
  stopifnot(effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta_ec50, cfg, "continuous_ploidy_amplified") >
              effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta, cfg, "shared_per_copy"))

  theta_hilln <- modifyList(default_theta(cfg, "continuous_ploidy_amplified"), list(ec50_alpha = 0.0, hilln_alpha = 0.5, pmis_alpha = 0.0))
  stopifnot(hilln_by_state(cfg$N_ref_high_state, theta_hilln, cfg, "continuous_ploidy_amplified") >
              hilln_by_state(cfg$N_ref_2N_state, theta_hilln, cfg, "continuous_ploidy_amplified"))

  theta_pmis <- modifyList(default_theta(cfg, "continuous_ploidy_amplified"), list(ec50_alpha = 0.0, hilln_alpha = 0.0, pmis_alpha = 0.5))
  stopifnot(pmis_multiplier_by_state(cfg$N_ref_high_state, theta_pmis, cfg, "continuous_ploidy_amplified") >
              pmis_multiplier_by_state(cfg$N_ref_2N_state, theta_pmis, cfg, "continuous_ploidy_amplified"))
  stopifnot(effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta_pmis, cfg, "continuous_ploidy_amplified") >
              effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta, cfg, "shared_per_copy"))

  theta_neutral_cont <- modifyList(default_theta(cfg, "continuous_ploidy_amplified"), list(ec50_alpha = 0.0, hilln_alpha = 0.0, pmis_alpha = 0.0))
  stopifnot(abs(effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta_neutral_cont, cfg, "continuous_ploidy_amplified") -
                  effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta, cfg, "shared_per_copy")) < 1e-12)

  theta_cat <- modifyList(default_theta(cfg, "categorical_high_ploidy"), list(a_4N_ec50 = 2.0, a_4N_hilln = 1.5, a_4N_pmis = 1.0))
  stopifnot(ec50_by_state(cfg$N_ref_high_state, theta_cat, cfg, "categorical_high_ploidy") < theta_cat$EC50)
  stopifnot(hilln_by_state(cfg$N_ref_high_state, theta_cat, cfg, "categorical_high_ploidy") > theta_cat$hill_n)

  theta_cat_pmis <- modifyList(default_theta(cfg, "categorical_high_ploidy"), list(a_4N_ec50 = 1.0, a_4N_hilln = 1.0, a_4N_pmis = 2.0))
  stopifnot(pmis_multiplier_by_state(cfg$N_ref_high_state, theta_cat_pmis, cfg, "categorical_high_ploidy") >
              pmis_multiplier_by_state(cfg$N_ref_2N_state, theta_cat_pmis, cfg, "categorical_high_ploidy"))
  stopifnot(effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta_cat_pmis, cfg, "categorical_high_ploidy") >
              effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta, cfg, "shared_per_copy"))

  theta_neutral_cat <- modifyList(default_theta(cfg, "categorical_high_ploidy"), list(a_4N_ec50 = 1.0, a_4N_hilln = 1.0, a_4N_pmis = 1.0))
  stopifnot(abs(effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta_neutral_cat, cfg, "categorical_high_ploidy") -
                  effective_pmis_by_state(0.1, cfg$N_ref_high_state, theta, cfg, "shared_per_copy")) < 1e-12)

  cfg_fix_shared <- cfg
  cfg_fix_shared$fixed_params <- list(EC50 = 0.03)
  theta_fix_shared <- default_theta(cfg_fix_shared, "shared_per_copy")
  stopifnot(abs(theta_fix_shared$EC50 - 0.03) < 1e-12)
  stopifnot(length(par_from_theta(theta_fix_shared, cfg_fix_shared, "shared_per_copy")) ==
              variant_param_count("shared_per_copy", cfg_fix_shared))

  cfg_fix_variant <- cfg
  cfg_fix_variant$fixed_params <- list(a_4N_pmis = 1.0, nullisomy_dirichlet_alpha = 0.25)
  theta_fix_variant <- default_theta(cfg_fix_variant, "categorical_high_ploidy")
  stopifnot(abs(theta_fix_variant$a_4N_pmis - 1.0) < 1e-12)
  stopifnot(abs(theta_fix_variant$nullisomy_dirichlet_alpha - 0.25) < 1e-12)
  par_fix_variant <- par_from_theta(theta_fix_variant, cfg_fix_variant, "categorical_high_ploidy")
  stopifnot(length(par_fix_variant) == variant_param_count("categorical_high_ploidy", cfg_fix_variant))
  stopifnot(length(fit_bounds("categorical_high_ploidy", cfg_fix_variant)$lower) == length(par_fix_variant))
  theta_reconstructed <- theta_from_par(par_fix_variant, cfg_fix_variant, "categorical_high_ploidy")
  stopifnot(abs(theta_reconstructed$a_4N_pmis - 1.0) < 1e-12)
  stopifnot(abs(theta_reconstructed$nullisomy_dirichlet_alpha - 0.25) < 1e-12)
  jittered_fix_variant <- jitter_theta(theta_fix_variant, cfg_fix_variant, "categorical_high_ploidy", seed = 123)
  stopifnot(abs(jittered_fix_variant$a_4N_pmis - 1.0) < 1e-12)
  stopifnot(abs(jittered_fix_variant$nullisomy_dirichlet_alpha - 0.25) < 1e-12)

  cfg_fix_multi <- cfg
  cfg_fix_multi$fixed_params <- parse_fixed_params("EC50=0.03,nullisomy_dirichlet_alpha=0.25")
  theta_fix_multi <- default_theta(cfg_fix_multi, "shared_per_copy")
  pred_fix_multi <- predict_dataset(theta_fix_multi, data_obj, cfg_fix_multi, "shared_per_copy")
  stopifnot(all(is.finite(pred_fix_multi$predicted)))

  cfg_fix_scoped <- cfg
  cfg_fix_scoped$fixed_params <- list(nullisomy_dirichlet_alpha = 0.25)
  cfg_fix_scoped$fixed_params_by_variant <- list(
    continuous_ploidy_amplified = list(hilln_alpha = 0),
    categorical_high_ploidy = list(a_4N_ec50 = 1, a_4N_hilln = 1)
  )
  theta_fix_scoped_shared <- default_theta(cfg_fix_scoped, "shared_per_copy")
  theta_fix_scoped_cont <- default_theta(cfg_fix_scoped, "continuous_ploidy_amplified")
  theta_fix_scoped_cat <- default_theta(cfg_fix_scoped, "categorical_high_ploidy")
  stopifnot(is.null(theta_fix_scoped_shared$hilln_alpha))
  stopifnot(abs(theta_fix_scoped_shared$nullisomy_dirichlet_alpha - 0.25) < 1e-12)
  stopifnot(abs(theta_fix_scoped_cont$hilln_alpha - 0) < 1e-12)
  stopifnot(abs(theta_fix_scoped_cat$a_4N_ec50 - 1) < 1e-12)
  stopifnot(abs(theta_fix_scoped_cat$a_4N_hilln - 1) < 1e-12)

  cfg_joint <- cfg
  cfg_joint$dataset_ids <- c("A", "B")
  theta_joint <- default_theta(cfg_joint, "shared_per_copy")
  stopifnot(all(c("EC50__A", "EC50__B") %in% names(theta_joint)))
  stopifnot(all(c("lambda_eff__A", "lambda_eff__B") %in% names(theta_joint)))
  stopifnot(length(par_from_theta(theta_joint, cfg_joint, "shared_per_copy")) ==
              variant_param_count("shared_per_copy", cfg_joint))
  datasets_joint <- list(
    A = list(dataset_id = "A", data_obj = data_obj),
    B = list(dataset_id = "B", data_obj = data_obj)
  )
  pred_A_ref <- predict_dataset_single(theta_joint, data_obj, cfg_joint, "shared_per_copy", dataset_id = "A")
  pred_joint_ref <- predict_dataset(theta_joint, datasets_joint, cfg_joint, "shared_per_copy")
  pred_joint_A <- pred_joint_ref[pred_joint_ref$dataset_id == "A", c("sample", "dose_uM", "predicted")]
  merged_joint <- merge(pred_A_ref[, c("sample", "dose_uM", "predicted")], pred_joint_A, by = c("sample", "dose_uM"), suffixes = c("_single", "_joint"))
  stopifnot(all(abs(merged_joint$predicted_single - merged_joint$predicted_joint) < 1e-12))
  theta_joint_shift <- theta_joint
  theta_joint_shift$lambda_eff__A <- theta_joint_shift$lambda_eff__A * 1.2
  theta_joint_shift$lambda_eff_by_dataset$A <- theta_joint_shift$lambda_eff__A
  pred_joint_shift <- predict_dataset(theta_joint_shift, datasets_joint, cfg_joint, "shared_per_copy")
  pred_B_ref <- pred_joint_ref[pred_joint_ref$dataset_id == "B", c("sample", "dose_uM", "predicted")]
  pred_B_shift <- pred_joint_shift[pred_joint_shift$dataset_id == "B", c("sample", "dose_uM", "predicted")]
  merged_B <- merge(pred_B_ref, pred_B_shift, by = c("sample", "dose_uM"), suffixes = c("_ref", "_shift"))
  stopifnot(all(abs(merged_B$predicted_ref - merged_B$predicted_shift) < 1e-12))

  theta_joint_ec50_shift <- theta_joint
  theta_joint_ec50_shift$EC50__A <- theta_joint_ec50_shift$EC50__A * 0.5
  theta_joint_ec50_shift$EC50_by_dataset$A <- theta_joint_ec50_shift$EC50__A
  pred_joint_ec50_shift <- predict_dataset(theta_joint_ec50_shift, datasets_joint, cfg_joint, "shared_per_copy")
  pred_ec50_B_ref <- pred_joint_ref[pred_joint_ref$dataset_id == "B", c("sample", "dose_uM", "predicted")]
  pred_ec50_B_shift <- pred_joint_ec50_shift[pred_joint_ec50_shift$dataset_id == "B", c("sample", "dose_uM", "predicted")]
  merged_ec50_B <- merge(pred_ec50_B_ref, pred_ec50_B_shift, by = c("sample", "dose_uM"), suffixes = c("_ref", "_shift"))
  stopifnot(all(abs(merged_ec50_B$predicted_ref - merged_ec50_B$predicted_shift) < 1e-12))
  pred_ec50_A_ref <- pred_joint_ref[pred_joint_ref$dataset_id == "A", c("sample", "dose_uM", "predicted")]
  pred_ec50_A_shift <- pred_joint_ec50_shift[pred_joint_ec50_shift$dataset_id == "A", c("sample", "dose_uM", "predicted")]
  merged_ec50_A <- merge(pred_ec50_A_ref, pred_ec50_A_shift, by = c("sample", "dose_uM"), suffixes = c("_ref", "_shift"))
  stopifnot(any(abs(merged_ec50_A$predicted_ref - merged_ec50_A$predicted_shift) > 1e-12))

  cfg_shared_ec50 <- cfg_joint
  cfg_shared_ec50$ec50_dataset_mode <- "shared"
  theta_shared_ec50 <- default_theta(cfg_shared_ec50, "shared_per_copy")
  stopifnot("EC50" %in% names(theta_shared_ec50))
  stopifnot(!any(c("EC50__A", "EC50__B") %in% names(theta_shared_ec50)))
  stopifnot(abs(theta_shared_ec50$EC50_by_dataset$A - theta_shared_ec50$EC50_by_dataset$B) < 1e-12)

  cfg_hafner <- cfg_joint
  cfg_hafner$ec50_dataset_mode <- "hafner_lambda"
  theta_hafner <- default_theta(cfg_hafner, "shared_per_copy")
  stopifnot(all(c("S50", "S_M") %in% names(theta_hafner)))
  stopifnot(!any(c("EC50__A", "EC50__B", "EC50") %in% names(theta_hafner)))
  free_bounds <- fit_bounds("shared_per_copy", cfg_joint)
  hafner_bounds <- fit_bounds("shared_per_copy", cfg_hafner)
  free_lambda_idx <- which(parameter_spec("shared_per_copy", cfg_joint)$name == "lambda_eff__A")
  hafner_lambda_idx <- which(parameter_spec("shared_per_copy", cfg_hafner)$name == "lambda_eff__A")
  stopifnot(abs(opt_to_natural("lambda_eff__A", free_bounds$lower[[free_lambda_idx]], cfg_joint, "log10") - 1e-3) < 1e-12)
  stopifnot(abs(opt_to_natural("lambda_eff__A", hafner_bounds$lower[[hafner_lambda_idx]], cfg_hafner, "log10") - 0.1) < 1e-12)
  theta_hafner$lambda_eff__A <- 0.25
  theta_hafner$lambda_eff__B <- 0.50
  theta_hafner <- refresh_dataset_parameter_maps(theta_hafner, cfg_hafner, fail_invalid_ec50 = TRUE)
  stopifnot(theta_hafner$EC50_by_dataset$A > theta_hafner$EC50_by_dataset$B)
  pred_hafner_ref <- predict_dataset(theta_hafner, datasets_joint, cfg_hafner, "shared_per_copy")
  theta_hafner_shift <- theta_hafner
  theta_hafner_shift$lambda_eff__A <- 0.30
  theta_hafner_shift <- refresh_dataset_parameter_maps(theta_hafner_shift, cfg_hafner, fail_invalid_ec50 = TRUE)
  pred_hafner_shift <- predict_dataset(theta_hafner_shift, datasets_joint, cfg_hafner, "shared_per_copy")
  pred_hafner_B_ref <- pred_hafner_ref[pred_hafner_ref$dataset_id == "B", c("sample", "dose_uM", "predicted")]
  pred_hafner_B_shift <- pred_hafner_shift[pred_hafner_shift$dataset_id == "B", c("sample", "dose_uM", "predicted")]
  merged_hafner_B <- merge(pred_hafner_B_ref, pred_hafner_B_shift, by = c("sample", "dose_uM"), suffixes = c("_ref", "_shift"))
  stopifnot(all(abs(merged_hafner_B$predicted_ref - merged_hafner_B$predicted_shift) < 1e-12))
  pred_hafner_A_ref <- pred_hafner_ref[pred_hafner_ref$dataset_id == "A", c("sample", "dose_uM", "predicted")]
  pred_hafner_A_shift <- pred_hafner_shift[pred_hafner_shift$dataset_id == "A", c("sample", "dose_uM", "predicted")]
  merged_hafner_A <- merge(pred_hafner_A_ref, pred_hafner_A_shift, by = c("sample", "dose_uM"), suffixes = c("_ref", "_shift"))
  stopifnot(any(abs(merged_hafner_A$predicted_ref - merged_hafner_A$predicted_shift) > 1e-12))
  theta_hafner_bad <- theta_hafner
  theta_hafner_bad$S_M <- 0.1
  stopifnot(!hafner_ec50_is_valid(theta_hafner_bad, cfg_hafner))
  cfg_hafner_fixed <- cfg_hafner
  cfg_hafner_fixed$fixed_params <- list(S50 = 0.2, S_M = 1.0)
  theta_hafner_fixed <- default_theta(cfg_hafner_fixed, "shared_per_copy")
  stopifnot(length(par_from_theta(theta_hafner_fixed, cfg_hafner_fixed, "shared_per_copy")) ==
              variant_param_count("shared_per_copy", cfg_hafner_fixed))
  cfg_hafner_bad_fix <- cfg_hafner
  cfg_hafner_bad_fix$fixed_params <- list(EC50__A = 0.2)
  stopifnot(inherits(try(default_theta(cfg_hafner_bad_fix, "shared_per_copy"), silent = TRUE), "try-error"))

  ploidy_tests <- data.frame(
    check = c(
      "response_table_uses_hill_smoothed_values_when_available",
      "response_table_accepts_nM_schema_and_converts_to_uM",
      "fixed_shared_parameter_removed_from_optimizer_vector",
      "fixed_variant_parameters_reconstruct_theta_correctly",
      "fit_bounds_match_free_parameter_vector_with_fixed_params",
      "fixed_parameters_remain_unchanged_after_jitter",
      "predict_dataset_runs_with_multiple_fixed_params",
      "variant_scoped_fixed_params_apply_only_to_their_variant",
      "joint_theta_includes_one_EC50_per_dataset",
      "joint_theta_includes_one_lambda_eff_per_dataset",
      "joint_prediction_matches_single_dataset_prediction_when_inputs_match",
      "changing_one_dataset_lambda_eff_does_not_move_the_other_dataset",
      "changing_one_dataset_EC50_does_not_move_the_other_dataset",
      "changing_one_dataset_EC50_does_move_that_dataset",
      "shared_ec50_mode_has_one_EC50_for_all_datasets",
      "hafner_ec50_mode_fits_S50_and_SM_not_EC50",
      "hafner_mode_uses_lambda_eff_lower_bound_0p1",
      "hafner_lower_lambda_eff_increases_derived_EC50",
      "changing_one_dataset_lambda_eff_in_hafner_mode_does_not_move_other_dataset",
      "changing_one_dataset_lambda_eff_in_hafner_mode_does_move_that_dataset",
      "hafner_invalid_denominator_is_detected",
      "hafner_fixed_S50_and_SM_removed_from_optimizer_vector",
      "hafner_mode_rejects_fixed_free_EC50",
      "response_top_controls_zero_dose_level",
      "zero_dose_is_not_forced_to_one_when_response_top_changes",
      "continuous_positive_ec50_alpha_lowers_high_state_EC50",
      "continuous_positive_hilln_alpha_raises_high_state_hill_n",
      "continuous_positive_pmis_alpha_raises_high_state_multiplier",
      "continuous_neutral_triplet_matches_shared_model",
      "categorical_a4n_ec50_gt1_lowers_high_state_EC50",
      "categorical_a4n_hilln_gt1_raises_high_state_hill_n",
      "categorical_a4n_pmis_gt1_raises_high_state_multiplier",
      "categorical_neutral_triplet_matches_shared_model",
      "continuous_positive_ec50_alpha_increases_high_state_pmis_at_mid_dose"
    ),
    passed = rep(TRUE, 34L),
    stringsAsFactors = FALSE
  )
  utils::write.table(ploidy_tests, file.path(out_dir, "ploidy_specific_hill_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  pred_shared <- predict_dataset(theta, data_obj, cfg, "shared_per_copy")
  stopifnot(all(abs(pred_shared$predicted[pred_shared$dose_uM == 0] - theta$response_top) < 1e-10))

  theta_top <- default_theta(cfg, "shared_per_copy")
  theta_top$response_top <- 0.8
  pred_top <- predict_dataset(theta_top, data_obj, cfg, "shared_per_copy")
  stopifnot(all(abs(pred_top$predicted[pred_top$dose_uM == 0] - 0.8) < 1e-10))
  stopifnot(any(abs(pred_top$predicted[pred_top$dose_uM == 0] - pred_shared$predicted[pred_shared$dose_uM == 0]) > 1e-6))

  theta_cli <- theta_from_cli_args(
    list(
      model = "categorical_high_ploidy",
      p_mis_doxo_max = "0.12",
      EC50_uM = "0.03",
      hill_n = "2.4",
      lambda_eff = "0.4",
      response_top = "0.85",
      a_4N_ec50 = "2.5",
      a_4N_hilln = "1.2",
      a_4N_pmis = "1.8"
    ),
    cfg,
    "shared_per_copy"
  )
  stopifnot(identical(canonical_model_variant("categorical_high_ploidy"), "categorical_high_ploidy"))
  stopifnot(abs(theta_cli$p_mis_doxo_max - 0.12) < 1e-12)
  stopifnot(abs(theta_cli$EC50 - 0.03) < 1e-12)
  stopifnot(abs(theta_cli$hill_n - 2.4) < 1e-12)
  stopifnot(abs(theta_cli$response_top - 0.85) < 1e-12)
  stopifnot(abs(theta_cli$a_4N_ec50 - 2.5) < 1e-12)
  stopifnot(abs(theta_cli$a_4N_hilln - 1.2) < 1e-12)
  stopifnot(abs(theta_cli$a_4N_pmis - 1.8) < 1e-12)

  predict_mode_dir <- file.path(out_dir, "predict_mode_smoke")
  dir.create(predict_mode_dir, recursive = TRUE, showWarnings = FALSE)
  predict_res <- run_predict_mode(
    data_obj,
    cfg,
    predict_mode_dir,
    list(
      model = "shared_per_copy",
      p_mis_doxo_max = "0.08",
      EC50 = "0.04",
      hill_n = "1.9",
      lambda_eff = "0.3",
      response_top = "0.77"
    )
  )
  zero_pred <- predict_res$fit$predictions$predicted[predict_res$fit$predictions$dose_uM == 0]
  stopifnot(length(zero_pred) >= 1L)
  stopifnot(all(abs(zero_pred - 0.77) < 1e-10))
  stopifnot(file.exists(file.path(predict_mode_dir, "fit_observed_vs_predicted_all_models.png")))

  cfg_large <- cfg
  cfg_large$assay_N0 <- cfg$o2_Nref
  pred_large <- predict_dataset(theta, data_obj, cfg_large, "shared_per_copy")
  merged_pred <- merge(
    pred_shared[, c("sample", "dose_uM", "predicted")],
    pred_large[, c("sample", "dose_uM", "predicted")],
    by = c("sample", "dose_uM"),
    suffixes = c("_3000", "_1e6"),
    sort = FALSE
  )
  stopifnot(all(abs(merged_pred$predicted_3000 - merged_pred$predicted_1e6) < 1e-8))

  sigma_zero <- objective_from_table(pred_shared, cfg, sigma_seed = 0.0)
  sigma_nonzero <- objective_from_table(pred_shared, cfg, sigma_seed = 0.2)
  stopifnot(is.finite(sigma_zero), is.finite(sigma_nonzero))
  stopifnot(abs(mean(pred_shared$predicted) - mean(pred_shared$predicted)) < 1e-12)
  stopifnot(sigma_nonzero < sigma_zero)

  cfg_volume <- cfg
  cfg_volume$observable_mode <- "volume_weighted_burden"
  pred_volume <- predict_dataset(theta, data_obj, cfg_volume, "shared_per_copy")
  stopifnot(nrow(pred_volume) == nrow(pred_shared))
  stopifnot(all(is.finite(pred_volume$predicted)))

  fake_fit <- list(model = "shared_per_copy", seed = 2L, objective = 5.0)
  fake_summary <- data.frame(
    model = c("shared_per_copy", "shared_per_copy"),
    seed = c(1L, 2L),
    objective = c(10.0, 5.0),
    convergence = c(0L, 1L),
    n_params = c(4L, 4L),
    aic = c(28.0, 18.0),
    stringsAsFactors = FALSE
  )
  fake_row <- fit_summary_row_for_fit(fake_summary, fake_fit)
  stopifnot(fake_row$seed[[1]] == 2L)
  stopifnot(abs(fake_row$aic[[1]] - 18.0) < 1e-12)
  stopifnot(fake_row$convergence[[1]] == 1L)

  writeLines("All doxorubicin-nullisomy checks passed.", con = file.path(out_dir, "tests_ok.txt"))
  invisible(TRUE)
}

theta_to_named_rows <- function(theta, variant, cfg) {
  variant <- canonical_model_variant(variant)
  dataset_ids <- names(theta$lambda_eff_by_dataset)
  ec50_rows <- data.frame(
    parameter = if (length(theta$EC50_by_dataset) == 1L) "EC50_uM" else paste0("EC50_uM_", names(theta$EC50_by_dataset)),
    value = unlist(theta$EC50_by_dataset, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  lambda_rows <- data.frame(
    parameter = if (length(theta$lambda_eff_by_dataset) == 1L) "lambda_eff" else paste0("lambda_eff_", dataset_ids),
    value = unlist(theta$lambda_eff_by_dataset, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  rows <- data.frame(
    parameter = c("p_mis_doxo_max", "hill_n", "response_top", "nullisomy_dirichlet_alpha"),
    value = c(theta$p_mis_doxo_max, theta$hill_n, theta$response_top, theta$nullisomy_dirichlet_alpha),
    stringsAsFactors = FALSE
  )
  rows <- bind_rows(rows[1, , drop = FALSE], ec50_rows, rows[2, , drop = FALSE], lambda_rows, rows[3:4, , drop = FALSE])
  if (identical(canonical_ec50_dataset_mode(cfg$ec50_dataset_mode), "hafner_lambda")) {
    rows <- bind_rows(
      rows,
      data.frame(parameter = "S50_uM", value = theta$S50, stringsAsFactors = FALSE),
      data.frame(parameter = "S_M", value = theta$S_M, stringsAsFactors = FALSE)
    )
  }
  if (identical(variant, "continuous_ploidy_amplified")) {
    rows <- bind_rows(
      rows,
      data.frame(parameter = "ec50_alpha", value = theta$ec50_alpha, stringsAsFactors = FALSE),
      data.frame(parameter = "hilln_alpha", value = theta$hilln_alpha, stringsAsFactors = FALSE),
      data.frame(parameter = "pmis_alpha", value = theta$pmis_alpha, stringsAsFactors = FALSE)
    )
  } else if (identical(variant, "categorical_high_ploidy")) {
    rows <- bind_rows(
      rows,
      data.frame(parameter = "a_4N_ec50", value = theta$a_4N_ec50, stringsAsFactors = FALSE),
      data.frame(parameter = "a_4N_hilln", value = theta$a_4N_hilln, stringsAsFactors = FALSE),
      data.frame(parameter = "a_4N_pmis", value = theta$a_4N_pmis, stringsAsFactors = FALSE)
    )
  }
  rows
}

param_row_num <- function(parameter, value) {
  data.frame(
    parameter = as.character(parameter),
    value_num = as.numeric(value),
    value_text = NA_character_,
    stringsAsFactors = FALSE
  )
}

param_row_text <- function(parameter, value) {
  data.frame(
    parameter = as.character(parameter),
    value_num = NA_real_,
    value_text = as.character(value),
    stringsAsFactors = FALSE
  )
}

theta_to_param_rows <- function(theta, variant, cfg) {
  num_rows <- theta_to_named_rows(theta, variant, cfg)
  data.frame(
    parameter = as.character(num_rows$parameter),
    value_num = as.numeric(num_rows$value),
    value_text = NA_character_,
    stringsAsFactors = FALSE
  )
}

theta_from_cli_args <- function(args, cfg, variant) {
  variant <- canonical_model_variant(first_non_null(args$model, args$variant, variant))
  theta <- default_theta(cfg, variant)
  dataset_ids <- dataset_ids_from_cfg(cfg)
  ec50_mode <- canonical_ec50_dataset_mode(cfg$ec50_dataset_mode)

  theta$p_mis_doxo_max <- as_num(args$p_mis_doxo_max, theta$p_mis_doxo_max)
  if (identical(ec50_mode, "hafner_lambda")) {
    theta$S50 <- as_num(first_non_null(args$S50, args$S50_uM), theta$S50)
    theta$S_M <- as_num(args$S_M, theta$S_M)
  } else if (identical(ec50_mode, "shared")) {
    theta$EC50 <- as_num(first_non_null(args$EC50, args$EC50_uM), theta$EC50)
  } else if (length(dataset_ids) == 1L) {
    theta$EC50 <- as_num(first_non_null(args$EC50, args$EC50_uM), theta$EC50)
  } else {
    for (dataset_id in dataset_ids) {
      nm <- paste0("EC50__", dataset_id)
      theta[[nm]] <- as_num(first_non_null(args[[nm]], args[[paste0("EC50_uM__", dataset_id)]]), theta[[nm]])
    }
  }
  theta$hill_n <- as_num(args$hill_n, theta$hill_n)
  theta$response_top <- as_num(args$response_top, theta$response_top)
  theta$nullisomy_dirichlet_alpha <- as_num(args$nullisomy_dirichlet_alpha, theta$nullisomy_dirichlet_alpha)
  for (dataset_id in dataset_ids) {
    nm <- paste0("lambda_eff__", dataset_id)
    theta[[nm]] <- as_num(first_non_null(args[[nm]], if (length(dataset_ids) == 1L) args$lambda_eff else NULL), theta[[nm]])
  }
  theta$lambda_eff_by_dataset <- setNames(
    as.list(vapply(dataset_ids, function(id) theta[[paste0("lambda_eff__", id)]], numeric(1))),
    dataset_ids
  )
  if (length(dataset_ids) == 1L) theta$lambda_eff <- theta[[paste0("lambda_eff__", dataset_ids[[1]])]]
  theta <- refresh_dataset_parameter_maps(theta, cfg, fail_invalid_ec50 = TRUE)

  if (identical(variant, "continuous_ploidy_amplified")) {
    theta$ec50_alpha <- as_num(args$ec50_alpha, theta$ec50_alpha)
    theta$hilln_alpha <- as_num(args$hilln_alpha, theta$hilln_alpha)
    theta$pmis_alpha <- as_num(args$pmis_alpha, theta$pmis_alpha)
  } else if (identical(variant, "categorical_high_ploidy")) {
    theta$a_4N_ec50 <- as_num(args$a_4N_ec50, theta$a_4N_ec50)
    theta$a_4N_hilln <- as_num(args$a_4N_hilln, theta$a_4N_hilln)
    theta$a_4N_pmis <- as_num(args$a_4N_pmis, theta$a_4N_pmis)
  }

  theta
}

jitter_theta <- function(theta, cfg, variant, seed) {
  variant <- canonical_model_variant(variant)
  set.seed(seed)
  fixed <- fixed_param_names(cfg, variant)
  maybe_set <- function(name, value) {
    if (!(name %in% fixed)) theta[[name]] <<- value
  }
  maybe_set("p_mis_doxo_max", clip01(theta$p_mis_doxo_max * 10^runif(1, -0.5, 0.5)))
  ec50_mode <- canonical_ec50_dataset_mode(cfg$ec50_dataset_mode)
  if (identical(ec50_mode, "hafner_lambda")) {
    maybe_set("S50", theta$S50 * 10^runif(1, -0.75, 0.75))
    maybe_set("S_M", theta$S_M * 10^runif(1, -0.3, 0.3))
  } else if (identical(ec50_mode, "shared")) {
    maybe_set("EC50", theta$EC50 * 10^runif(1, -0.75, 0.75))
  } else {
    for (dataset_id in dataset_ids_from_cfg(cfg)) {
      nm <- if (length(dataset_ids_from_cfg(cfg)) == 1L) "EC50" else paste0("EC50__", dataset_id)
      maybe_set(nm, theta[[nm]] * 10^runif(1, -0.75, 0.75))
    }
  }
  maybe_set("hill_n", theta$hill_n * 10^runif(1, -0.5, 0.5))
  for (dataset_id in dataset_ids_from_cfg(cfg)) {
    nm <- paste0("lambda_eff__", dataset_id)
    maybe_set(nm, theta[[nm]] * 10^runif(1, -0.4, 0.4))
  }
  maybe_set("response_top", theta$response_top * 10^runif(1, -0.2, 0.2))
  maybe_set("nullisomy_dirichlet_alpha", 10^runif(1, log10(0.5), log10(100)))

  if (identical(variant, "continuous_ploidy_amplified")) {
    strong_ec50_family <- as.logical(rbinom(1, 1, 0.5))
    if (strong_ec50_family) {
      maybe_set("ec50_alpha", runif(1, 0.4, 1.5))
      maybe_set("hilln_alpha", runif(1, -0.5, 0.2))
      maybe_set("pmis_alpha", runif(1, 0.5, 2.0))
    } else {
      maybe_set("ec50_alpha", runif(1, -1.0, 1.5))
      maybe_set("hilln_alpha", runif(1, -0.5, 0.5))
      maybe_set("pmis_alpha", runif(1, -0.5, 2.0))
    }
  } else if (identical(variant, "categorical_high_ploidy")) {
    strong_ec50_family <- as.logical(rbinom(1, 1, 0.5))
    if (strong_ec50_family) {
      maybe_set("a_4N_ec50", 10^runif(1, 0.4, 0.9))
      maybe_set("a_4N_hilln", 10^runif(1, -0.1, 0.15))
      maybe_set("a_4N_pmis", 10^runif(1, 0.4, 0.9))
    } else {
      maybe_set("a_4N_ec50", 10^runif(1, -0.2, 0.2))
      maybe_set("a_4N_hilln", 10^runif(1, -0.2, 0.2))
      maybe_set("a_4N_pmis", 10^runif(1, -0.2, 0.2))
    }
  }
  theta <- refresh_dataset_parameter_maps(theta, cfg, fail_invalid_ec50 = FALSE)
  ensure_valid_hafner_initial_theta(theta, cfg)
}

implausibility_note_for_predictions <- function(pred_tab, cfg) {
  max_p <- max(c(pred_tab$p_eff_2N, pred_tab$p_eff_high), na.rm = TRUE)
  if (!is.finite(max_p)) return("Could not evaluate p_mis plausibility.")
  if (max_p > cfg$p_mis_implausible_threshold) {
    paste0("Effective per-copy p_mis reaches ", signif(max_p, 4), ", which is likely biologically implausible.")
  } else {
    paste0("Effective per-copy p_mis stays below ", signif(cfg$p_mis_implausible_threshold, 4), ".")
  }
}

assess_identifiability <- function(fit_rows, variant) {
  variant <- canonical_model_variant(variant)
  if (identical(variant, "shared_per_copy")) return("Not applicable.")
  conv <- fit_rows[fit_rows$convergence == 0, , drop = FALSE]
  if (nrow(conv) < 2L) return("Insufficient converged seeds to assess identifiability.")
  target_cols <- if (identical(variant, "continuous_ploidy_amplified")) {
    c("ec50_alpha", "hilln_alpha", "pmis_alpha")
  } else {
    c("a_4N_ec50", "a_4N_hilln", "a_4N_pmis")
  }
  unstable <- vapply(target_cols, function(col) {
    target <- as.numeric(conv[[col]])
    any(!is.finite(target)) || sd(target) > max(0.15, 0.5 * abs(mean(target)))
  }, logical(1))
  if (any(unstable)) {
    return("Weakly identifiable across seeds.")
  }
  "Reasonably stable across seeds."
}

fit_summary_row_for_fit <- function(fit_summary, fit) {
  rows <- fit_summary[fit_summary$model == fit$model, , drop = FALSE]
  if (nrow(rows) == 0L) stop("No fit_summary row found for model: ", fit$model)
  seed_match <- if ("seed" %in% names(rows) && !is.null(fit$seed)) which(rows$seed == fit$seed) else integer(0)
  obj_match <- which(abs(as.numeric(rows$objective) - as.numeric(fit$objective)) < 1e-8)
  idx <- intersect(seed_match, obj_match)
  if (length(idx) == 0L) idx <- obj_match
  if (length(idx) == 0L) idx <- which.min(abs(as.numeric(rows$objective) - as.numeric(fit$objective)))
  rows[idx[[1L]], , drop = FALSE]
}

write_fit_outputs <- function(best_fits, fit_summary, data_obj, cfg, out_dir) {
  pred_all <- bind_rows(lapply(best_fits, `[[`, "predictions"))
  pred_all$dose_uM_plot <- pmax(pred_all$dose_uM, 1e-4)
  residual <- transform(
    pred_all,
    residual_log = log(pmax(predicted, cfg$burden_log_eps)) - log(pmax(observed, cfg$burden_log_eps))
  )

  params_all <- bind_rows(lapply(best_fits, function(fit) {
    fit_fixed <- fixed_params_for_variant(cfg, fit$model)
    bind_rows(
      mutate(param_row_num("p_mis_base", cfg$p_mis_base), is_fixed = FALSE),
      mutate(param_row_num("assay_N0_cells", cfg$assay_N0), is_fixed = FALSE),
      mutate(param_row_text("observable_mode", cfg$observable_mode), is_fixed = FALSE),
      mutate(param_row_text("nullisomy_hidden_copy_mode", cfg$nullisomy_hidden_copy_mode), is_fixed = FALSE),
      mutate(theta_to_param_rows(fit$theta, fit$model, cfg), is_fixed = as.character(theta_to_param_rows(fit$theta, fit$model, cfg)$parameter) %in% names(fit_fixed)),
      mutate(param_row_num("assay_days_fixed", cfg$assay_days_fixed), is_fixed = FALSE),
      mutate(param_row_num("objective", fit$objective), is_fixed = FALSE)
    ) %>%
      mutate(model = fit$model, model_label = model_variant_label(fit$model), .before = 1)
  }))
  utils::write.table(params_all, file.path(out_dir, "best_params_by_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  fixed_params_rows <- bind_rows(lapply(names(best_fits), function(variant) {
    fp <- fixed_params_for_variant(cfg, variant)
    if (length(fp) == 0L) return(data.frame(model = character(0), parameter = character(0), value = numeric(0), stringsAsFactors = FALSE))
    data.frame(model = variant, parameter = names(fp), value = vapply(fp, as.numeric, numeric(1)), stringsAsFactors = FALSE)
  }))
  fixed_params_tab <- if (nrow(fixed_params_rows) == 0L) data.frame(model = character(0), parameter = character(0), value = numeric(0), stringsAsFactors = FALSE) else fixed_params_rows
  utils::write.table(fixed_params_tab, file.path(out_dir, "fixed_params_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(residual, file.path(out_dir, "predictions_all_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  fit_summary$aic <- 2 * fit_summary$objective + 2 * fit_summary$n_params
  utils::write.table(fit_summary, file.path(out_dir, "model_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  g_fit <- ggplot(pred_all, aes(dose_uM_plot, observed, color = sample)) +
    geom_point(size = 1.8) +
    geom_line(aes(y = predicted), linewidth = 0.9) +
    geom_errorbar(aes(ymin = pmax(observed - observed_sd, 1e-6), ymax = observed + observed_sd), width = 0.02, alpha = 0.5) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Observed vs predicted normalized doxorubicin response", x = "Dose (uM)", y = "Normalized viable burden")
  ggsave(file.path(out_dir, "fit_observed_vs_predicted_all_models.png"), g_fit, width = 11, height = 7, dpi = 150)

  g_res <- ggplot(residual, aes(dose_uM_plot, residual_log, color = sample)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 1.8) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Residuals by model, sample, and dose", x = "Dose (uM)", y = "log(pred) - log(obs)")
  ggsave(file.path(out_dir, "fit_residuals_all_models.png"), g_res, width = 11, height = 7, dpi = 150)

  p_curve <- pred_all %>%
    distinct(model, model_label, dose_uM, p_eff_2N, p_eff_high) %>%
    mutate(dose_uM_plot = pmax(dose_uM, 1e-4)) %>%
    tidyr::pivot_longer(c(p_eff_2N, p_eff_high), names_to = "state", values_to = "p_eff")
  utils::write.table(p_curve, file.path(out_dir, "effective_pmis_by_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  g_p <- ggplot(p_curve, aes(dose_uM_plot, p_eff, color = state)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Effective per-copy p_mis by model", x = "Dose (uM)", y = "Per-copy p_mis")
  ggsave(file.path(out_dir, "effective_pmis_by_model.png"), g_p, width = 11, height = 7, dpi = 150)

  mis_curve <- pred_all %>%
    distinct(model, model_label, dose_uM, expected_mis_2N, expected_mis_high) %>%
    mutate(dose_uM_plot = pmax(dose_uM, 1e-4)) %>%
    tidyr::pivot_longer(c(expected_mis_2N, expected_mis_high), names_to = "state", values_to = "expected_mis")
  utils::write.table(mis_curve, file.path(out_dir, "expected_missegregated_chromosomes_by_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  g_m <- ggplot(mis_curve, aes(dose_uM_plot, expected_mis, color = state)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Expected missegregated chromosomes per division", x = "Dose (uM)", y = "Expected missegregated chromosomes")
  ggsave(file.path(out_dir, "expected_missegregated_chromosomes_by_model.png"), g_m, width = 11, height = 7, dpi = 150)

  report_lines <- c("Doxorubicin-nullisomy model comparison", "")
  for (fit in best_fits) {
    variant <- fit$model
    fit_tab <- fit_summary_row_for_fit(fit_summary, fit)
    pred_tab <- fit$predictions
    convergence_note <- if ("convergence" %in% names(fit_tab) && is.finite(as.numeric(fit_tab$convergence[[1]])) && as.numeric(fit_tab$convergence[[1]]) != 0) {
      paste("Convergence warning: selected best fit has convergence =", fit_tab$convergence[[1]], "(optim did not report successful convergence).")
    } else {
      paste("Convergence:", if ("convergence" %in% names(fit_tab)) fit_tab$convergence[[1]] else "not reported")
    }
    report_lines <- c(
      report_lines,
      paste0("[", model_variant_label(variant), "]"),
      paste("Objective:", signif(fit$objective, 6)),
      paste("AIC:", signif(fit_tab$aic[[1]], 6)),
      if ("seed" %in% names(fit_tab)) paste("Selected seed:", fit_tab$seed[[1]]) else NULL,
      convergence_note,
      paste("Parameters:", paste(paste(theta_to_named_rows(fit$theta, variant, cfg)$parameter, signif(theta_to_named_rows(fit$theta, variant, cfg)$value, 6), sep = "="), collapse = ", ")),
      paste("Identifiability:", assess_identifiability(fit_summary[fit_summary$model == variant, , drop = FALSE], variant)),
      implausibility_note_for_predictions(pred_tab, cfg),
      ""
    )
  }
  best_fit_rows <- bind_rows(lapply(best_fits, function(fit) fit_summary_row_for_fit(fit_summary, fit)))
  best_variant <- best_fit_rows$model[[which.min(best_fit_rows$aic)]]
  report_lines <- c(
    report_lines,
    "Assay assumptions:",
    paste("Cells are seeded at", signif(cfg$assay_N0, 6), "cells per well and treatment begins at drug addition (t = 0)."),
    paste("Treatment duration is fixed at", signif(cfg$assay_days_fixed, 6), "days (96 h in the assay protocol)."),
    if (nrow(fixed_params_tab) > 0L) {
      paste("Fixed parameters:", paste(paste(paste0(fixed_params_tab$model, "::", fixed_params_tab$parameter), signif(fixed_params_tab$value, 6), sep = "="), collapse = ", "))
    } else {
      "Fixed parameters: none."
    },
    if (identical(canonical_observable_mode(cfg$observable_mode), "cell_count")) {
      "CellTiter-Glo signal is approximated as proportional to viable cell number."
    } else {
      "CellTiter-Glo signal is approximated as proportional to viable volume-weighted burden."
    },
    paste("Nullisomy risk is computed using", cfg$nullisomy_hidden_copy_mode, "hidden-copy buffering rather than the balanced hidden-copy baseline."),
    "Seeding variation is currently fixed to zero in the observation model, so normalized-response uncertainty comes only from the assay-derived sigma_obs term.",
    "In this deterministic density-independent mode, absolute starting cell number cancels in normalized response except for numerical roundoff.",
    "",
    paste("Best AIC model:", model_variant_label(best_variant)),
    "Amplification was tested without adding any direct death term."
  )
  writeLines(report_lines, con = file.path(out_dir, "technical_report.txt"))
}

run_fit <- function(data_obj, cfg, out_dir, seeds, variants = NULL) {
  if (is.null(variants)) variants <- cfg$fit_models_default
  variants <- vapply(variants, canonical_model_variant, character(1))
  all_fits <- unlist(
    lapply(variants, function(variant) lapply(seeds, function(seed) fit_one_seed(seed, data_obj, cfg, variant))),
    recursive = FALSE
  )
  fit_summary <- bind_rows(lapply(all_fits, function(x) {
    row <- data.frame(
      model = x$model,
      model_label = model_variant_label(x$model),
      seed = x$seed,
      objective = x$objective,
      convergence = x$convergence,
      n_params = variant_param_count(x$model, cfg),
      p_mis_doxo_max = x$theta$p_mis_doxo_max,
      hill_n = x$theta$hill_n,
      response_top = x$theta$response_top,
      nullisomy_dirichlet_alpha = x$theta$nullisomy_dirichlet_alpha,
      assay_days_fixed = cfg$assay_days_fixed,
      stringsAsFactors = FALSE
    )
    if (length(x$theta$EC50_by_dataset) == 1L) {
      row$EC50_uM <- unlist(x$theta$EC50_by_dataset, use.names = FALSE)[[1]]
    } else {
      for (dataset_id in names(x$theta$EC50_by_dataset)) {
        row[[paste0("EC50_uM_", dataset_id)]] <- x$theta$EC50_by_dataset[[dataset_id]]
      }
    }
    if (!is.null(x$theta$S50)) row$S50_uM <- x$theta$S50
    if (!is.null(x$theta$S_M)) row$S_M <- x$theta$S_M
    if (length(x$theta$lambda_eff_by_dataset) == 1L) {
      row$lambda_eff <- unlist(x$theta$lambda_eff_by_dataset, use.names = FALSE)[[1]]
    } else {
      for (dataset_id in names(x$theta$lambda_eff_by_dataset)) {
        row[[paste0("lambda_eff_", dataset_id)]] <- x$theta$lambda_eff_by_dataset[[dataset_id]]
      }
    }
    if (!is.null(x$theta$ec50_alpha)) row$ec50_alpha <- x$theta$ec50_alpha
    if (!is.null(x$theta$hilln_alpha)) row$hilln_alpha <- x$theta$hilln_alpha
    if (!is.null(x$theta$pmis_alpha)) row$pmis_alpha <- x$theta$pmis_alpha
    if (!is.null(x$theta$a_4N_ec50)) row$a_4N_ec50 <- x$theta$a_4N_ec50
    if (!is.null(x$theta$a_4N_hilln)) row$a_4N_hilln <- x$theta$a_4N_hilln
    if (!is.null(x$theta$a_4N_pmis)) row$a_4N_pmis <- x$theta$a_4N_pmis
    row
  }))
  utils::write.table(fit_summary, file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  best_fits <- lapply(variants, function(variant) {
    idx <- which(fit_summary$model == variant)
    all_fits[[idx[which.min(fit_summary$objective[idx])]]]
  })
  names(best_fits) <- variants
  write_fit_outputs(best_fits, fit_summary, data_obj, cfg, out_dir)
  invisible(list(best_fits = best_fits, fit_summary = fit_summary))
}

run_predict_mode <- function(data_obj, cfg, out_dir, args) {
  variant <- canonical_model_variant(first_non_null(args$model, args$variant, "shared_per_copy"))
  theta <- theta_from_cli_args(args, cfg, variant)
  pred <- predict_dataset(theta, data_obj, cfg, variant)
  fit <- list(
    model = variant,
    theta = theta,
    predictions = pred,
    objective = 0
  )
  fit_summary <- data.frame(
    model = variant,
    model_label = model_variant_label(variant),
    seed = NA_integer_,
    objective = 0,
    convergence = 0,
    n_params = variant_param_count(variant, cfg),
    p_mis_doxo_max = theta$p_mis_doxo_max,
    hill_n = theta$hill_n,
    response_top = theta$response_top,
    nullisomy_dirichlet_alpha = theta$nullisomy_dirichlet_alpha,
    assay_days_fixed = cfg$assay_days_fixed,
    stringsAsFactors = FALSE
  )
  if (length(theta$EC50_by_dataset) == 1L) {
    fit_summary$EC50_uM <- unlist(theta$EC50_by_dataset, use.names = FALSE)[[1]]
  } else {
    for (dataset_id in names(theta$EC50_by_dataset)) {
      fit_summary[[paste0("EC50_uM_", dataset_id)]] <- theta$EC50_by_dataset[[dataset_id]]
    }
  }
  if (!is.null(theta$S50)) fit_summary$S50_uM <- theta$S50
  if (!is.null(theta$S_M)) fit_summary$S_M <- theta$S_M
  if (length(theta$lambda_eff_by_dataset) == 1L) {
    fit_summary$lambda_eff <- unlist(theta$lambda_eff_by_dataset, use.names = FALSE)[[1]]
  } else {
    for (dataset_id in names(theta$lambda_eff_by_dataset)) {
      fit_summary[[paste0("lambda_eff_", dataset_id)]] <- theta$lambda_eff_by_dataset[[dataset_id]]
    }
  }
  if (!is.null(theta$ec50_alpha)) fit_summary$ec50_alpha <- theta$ec50_alpha
  if (!is.null(theta$hilln_alpha)) fit_summary$hilln_alpha <- theta$hilln_alpha
  if (!is.null(theta$pmis_alpha)) fit_summary$pmis_alpha <- theta$pmis_alpha
  if (!is.null(theta$a_4N_ec50)) fit_summary$a_4N_ec50 <- theta$a_4N_ec50
  if (!is.null(theta$a_4N_hilln)) fit_summary$a_4N_hilln <- theta$a_4N_hilln
  if (!is.null(theta$a_4N_pmis)) fit_summary$a_4N_pmis <- theta$a_4N_pmis

  write_fit_outputs(setNames(list(fit), variant), fit_summary, data_obj, cfg, out_dir)
  invisible(list(fit = fit, fit_summary = fit_summary))
}

main <- function(argv) {
  args <- parse_args(argv)
  mode <- tolower(trimws(as.character(if (is.null(args$mode)) "fit" else args$mode)))
  cfg <- default_cfg()
  if (!is.null(args$p_mis_base)) cfg$p_mis_base <- as_num(args$p_mis_base, cfg$p_mis_base)
  if (!is.null(args$gamma_loss)) cfg$gamma_loss <- as_num(args$gamma_loss, cfg$gamma_loss)
  if (!is.null(args$dt)) cfg$DT <- as_num(args$dt, cfg$DT)
  if (!is.null(args$assay_n0)) cfg$assay_N0 <- as_num(args$assay_n0, cfg$assay_N0)
  if (!is.null(args$observable_mode)) cfg$observable_mode <- canonical_observable_mode(args$observable_mode)
  if (!is.null(args$response_data_mode)) cfg$response_data_mode <- canonical_response_data_mode(args$response_data_mode)
  if (!is.null(args$ec50_dataset_mode)) cfg$ec50_dataset_mode <- canonical_ec50_dataset_mode(args$ec50_dataset_mode)
  if (!is.null(args$hafner_ec50_min_denom)) cfg$hafner_ec50_min_denom <- as_num(args$hafner_ec50_min_denom, cfg$hafner_ec50_min_denom)
  if (!is.null(args$optim_maxit)) cfg$optim_maxit <- as_int(args$optim_maxit, cfg$optim_maxit)
  cfg$fixed_params <- parse_fixed_params(args$fixed_params)
  cfg$fixed_params_by_variant <- parse_variant_scoped_fixed_params(args)
  if (!is.null(args$fixed_nullisomy_dirichlet_alpha)) {
    cfg$fixed_params$nullisomy_dirichlet_alpha <- as_num(args$fixed_nullisomy_dirichlet_alpha, NA_real_)
  }
  json_paths <- if (!is.null(args$jsons)) {
    strsplit(args$jsons, ",", fixed = TRUE)[[1]]
  } else if (!is.null(args$json)) {
    c(args$json)
  } else {
    default_data_paths()
  }
  json_paths <- trimws(as.character(json_paths))
  json_paths <- json_paths[nzchar(json_paths)]
  dataset_ids <- if (!is.null(args$dataset_ids)) strsplit(args$dataset_ids, ",", fixed = TRUE)[[1]] else NULL
  out_dir <- make_output_dir(args$out_dir)
  datasets <- load_doxo_datasets(json_paths, dataset_ids = dataset_ids)
  cfg$dataset_ids <- names(datasets)
  for (ds in datasets) validate_doxo_data(ds$data_obj, cfg)
  diagnostic_curves(cfg, out_dir)
  initial_scale_diagnostics(datasets, cfg, out_dir)
  if (mode == "diagnostics") return(invisible(TRUE))
  if (mode == "test") return(invisible(run_tests(datasets[[1]]$data_obj, cfg, out_dir)))
  if (mode == "predict") return(invisible(run_predict_mode(datasets, cfg, out_dir, args)))

  seeds <- if (!is.null(args$seeds_csv)) {
    as.integer(strsplit(args$seeds_csv, ",", fixed = TRUE)[[1]])
  } else {
    cfg$seeds_default
  }
  best_fit <- run_fit(datasets, cfg, out_dir, seeds)
  message("Wrote doxorubicin-nullisomy outputs to: ", out_dir)
  invisible(best_fit)
}

if (sys.nframe() == 0L) {
  main(commandArgs(trailingOnly = TRUE))
}
