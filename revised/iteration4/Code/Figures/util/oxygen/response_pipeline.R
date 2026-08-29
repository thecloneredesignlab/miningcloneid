#!/usr/bin/env Rscript

# Standalone Oxygen response generator.
#
# The complete analytical model and curve-classification functions required by
# Oxygen response are defined below as one public workspace utility implementation.
# No external R or C++ source file is loaded at runtime.
# The complete dependency closure is preserved below in this one file.

suppressPackageStartupMessages({
  library(Matrix)
  library(Rcpp)
})

# ---- embedded generic and semantic helpers -------------------------------
# -----------------------------------------------------------------------------
# Shared helpers for O2_supply_demand_MAP scripts.
# This file consolidates generic utilities and shared runtime helpers used by
# fit/viz/report/model code. It intentionally excludes configuration semantics,
# which stay in o2_supply_demand_map_common_semantics.R.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Basic Utils
# -----------------------------------------------------------------------------
o2sd_null_coalesce <- function(a, b) {
  if (is.null(a)) b else a
}

o2sd_first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

o2sd_parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0L) return(out)
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    body <- sub("^--", "", arg)
    parts <- strsplit(body, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    val <- if (length(parts) > 1L) paste(parts[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

o2sd_as_num <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) return(as.numeric(default))
  y <- suppressWarnings(as.numeric(x))
  if (length(y) == 1L && !is.finite(y)) return(as.numeric(default))
  y
}

o2sd_as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || !length(x)) return(as.integer(default))
  y <- suppressWarnings(as.integer(x))
  if (length(y) == 1L && is.na(y)) return(as.integer(default))
  y
}

o2sd_as_num_vec <- function(x, default = numeric(0)) {
  if (is.null(x) || !length(x)) return(as.numeric(default))
  s <- trimws(as.character(x[[1]]))
  if (!nzchar(s)) return(as.numeric(default))
  parts <- trimws(unlist(strsplit(s, "[,;]", perl = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) return(as.numeric(default))
  vals <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(vals))) stop("Invalid numeric vector argument: ", s)
  vals
}

o2sd_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(isTRUE(default))
  if (is.logical(x)) return(isTRUE(x[[1]]))
  s <- tolower(trimws(as.character(x[[1]])))
  if (!nzchar(s)) return(isTRUE(default))
  if (s %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (s %in% c("false", "f", "0", "no", "n")) return(FALSE)
  isTRUE(default)
}

o2sd_as_bool_scalar <- o2sd_as_bool

o2sd_clip <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}

o2sd_get_script_dir <- function(default = getwd()) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }

  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }

  normalizePath(default, mustWork = FALSE)
}

# -----------------------------------------------------------------------------
# Parameter Description Helpers
# -----------------------------------------------------------------------------

o2sd_parameter_natural_name <- function(x) {
  out <- trimws(as.character(x))
  out <- sub("^ivt__", "", out)
  out <- sub("^log10_", "", out)
  out <- sub("^logit01_", "", out)
  out <- sub("^logit_", "", out)
  out
}

o2sd_read_parameter_description_rows <- function(parameter_tables = character()) {
  rows <- list()
  paths <- unique(as.character(if (is.null(parameter_tables)) character(0) else parameter_tables))
  paths <- paths[nzchar(paths) & file.exists(paths)]
  if (!length(paths)) {
    return(data.frame(parameter = character(0), parameter_description = character(0), stringsAsFactors = FALSE))
  }
  for (path in paths) {
    tab <- tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (!is.data.frame(tab) || !nrow(tab)) next
    param_col <- intersect(c("param_symbol", "param_prototype", "parameter", "transformed_parameter", "param_name"), names(tab))
    desc_col <- intersect(c("parameter_description", "description", "note"), names(tab))
    if (!length(param_col) || !length(desc_col)) next
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = o2sd_parameter_natural_name(tab[[param_col[[1L]]]]),
      parameter_description = trimws(as.character(tab[[desc_col[[1L]]]])),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame(parameter = character(0), parameter_description = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

o2sd_parameter_description_map <- function(parameter_tables = character()) {
  rows <- o2sd_read_parameter_description_rows(parameter_tables)
  rows$parameter <- trimws(as.character(rows$parameter))
  rows$parameter_description <- trimws(as.character(rows$parameter_description))
  keep <- nzchar(rows$parameter) & nzchar(rows$parameter_description)
  rows <- rows[keep, , drop = FALSE]
  rows <- rows[!duplicated(rows$parameter), , drop = FALSE]
  stats::setNames(rows$parameter_description, rows$parameter)
}

o2sd_parameter_column <- function(tab) {
  if (!is.data.frame(tab) || !ncol(tab)) return(NULL)
  hits <- intersect(c("parameter", "param_prototype", "param_symbol", "transformed_parameter", "param_name"), names(tab))
  if (!length(hits)) NULL else hits[[1L]]
}

o2sd_add_parameter_descriptions <- function(tab,
                                           parameter_tables = character(),
                                           description_col = "parameter_description") {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(tab)
  param_col <- o2sd_parameter_column(tab)
  if (is.null(param_col)) return(tab)

  desc_map <- o2sd_parameter_description_map(parameter_tables)
  param_names <- trimws(as.character(tab[[param_col]]))
  natural_names <- o2sd_parameter_natural_name(param_names)
  desc <- unname(desc_map[natural_names])
  missing_desc <- is.na(desc) | !nzchar(desc)
  if (any(missing_desc)) {
    desc[missing_desc] <- unname(desc_map[param_names[missing_desc]])
  }
  if (description_col %in% names(tab)) {
    existing <- trimws(as.character(tab[[description_col]]))
    keep_existing <- !is.na(existing) & nzchar(existing)
    desc[keep_existing] <- existing[keep_existing]
  }
  missing_desc <- is.na(desc) | !nzchar(desc)
  if ("description" %in% names(tab)) {
    fallback <- trimws(as.character(tab$description))
    desc[missing_desc & nzchar(fallback)] <- fallback[missing_desc & nzchar(fallback)]
  }
  missing_desc <- is.na(desc) | !nzchar(desc)
  if ("note" %in% names(tab)) {
    fallback <- trimws(as.character(tab$note))
    desc[missing_desc & nzchar(fallback)] <- fallback[missing_desc & nzchar(fallback)]
  }
  desc[is.na(desc)] <- ""

  out <- tab
  out[[description_col]] <- desc
  ordered_names <- names(out)
  ordered_names <- ordered_names[ordered_names != description_col]
  insert_after <- match(param_col, ordered_names)
  ordered_names <- append(ordered_names, description_col, after = insert_after)
  out[, ordered_names, drop = FALSE]
}

# -----------------------------------------------------------------------------
# Shared Runtime Aliases
# -----------------------------------------------------------------------------
o2sd_runtime_as_num <- o2sd_as_num
o2sd_runtime_first_non_null <- o2sd_first_non_null

# Keep compatibility aliases for existing fit/viz/model code paths.
.first_non_null_local <- o2sd_first_non_null

harvest_init_log_param_name <- function(harvest) {
  paste0("log_init_mult_", as.character(harvest))
}

harvest_init_natural_param_name <- function(harvest) {
  paste0("init_mult_", as.character(harvest))
}

# -----------------------------------------------------------------------------
# Shared Runtime Helpers
# -----------------------------------------------------------------------------
default_rho_2N_prior_bounds <- function(cfg = NULL) {
  lo <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$rho_2N_min else NULL, 3.2e4))
  hi <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$rho_2N_max else NULL, 5.6e4))
  if (!is.finite(lo) || lo <= 0) lo <- 3.2e4
  if (!is.finite(hi) || hi <= 0) hi <- 5.6e4
  if (lo > hi) {
    tmp <- lo
    lo <- hi
    hi <- tmp
  }
  c(rho_2N_min = lo, rho_2N_max = hi)
}

default_rho_2N_prior_center <- function(cfg = NULL) {
  b <- default_rho_2N_prior_bounds(cfg)
  sqrt(b[["rho_2N_min"]] * b[["rho_2N_max"]])
}

cell_volume_mm3_by_ploidy <- function(ploidy, run_params, cfg) {
  p <- as.numeric(ploidy)
  rho_2N <- suppressWarnings(as.numeric(run_params$rho_2N))
  rho_2N <- if (length(rho_2N) > 0) rho_2N[[1]] else NA_real_
  if (is.na(rho_2N) || !is.finite(rho_2N) || rho_2N <= 0) rho_2N <- default_rho_2N_prior_center(cfg)
  rep(1 / rho_2N, length(p))
}

cell_volume_mm3_by_chr_number <- function(N, run_params, cfg, N_dip = 44) {
  N_use <- as.numeric(N)
  rho_2N <- suppressWarnings(as.numeric(run_params$rho_2N))
  rho_2N <- if (length(rho_2N) > 0) rho_2N[[1]] else NA_real_
  if (is.na(rho_2N) || !is.finite(rho_2N) || rho_2N <= 0) rho_2N <- default_rho_2N_prior_center(cfg)
  rep(1 / rho_2N, length(N_use))
}

cell_volume_mm3_by_N <- function(N, run_params, cfg) {
  start_with_mode <- assert_canonical_start_with_mode(
    .first_non_null_local(if (!is.null(cfg)) cfg$start_with else NULL, "ploidy")
  )
  if (identical(start_with_mode, "chr_number")) {
    return(cell_volume_mm3_by_chr_number(N, run_params = run_params, cfg = cfg, N_dip = 44))
  }
  p_weighted <- weighted_ploidy_from_total_N(
    N_total = as.numeric(N),
    chr_lengths_bp = cfg$chr_lengths_bp
  )
  cell_volume_mm3_by_ploidy(p_weighted, run_params = run_params, cfg = cfg)
}

resolve_terminal_ploidy_path <- function(data_dir) {
  candidates <- c(
    file.path(data_dir, "all_ploidy.csv"),
    file.path(data_dir, "all_ploidy.tsv")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("No terminal single-cell table found. Tried: ", paste(candidates, collapse = ", "))
  }
  normalizePath(hit[[1]], mustWork = TRUE)
}

read_terminal_ploidy_table <- function(ploidy_path) {
  ext <- tolower(tools::file_ext(ploidy_path))
  if (identical(ext, "csv")) {
    tab <- utils::read.csv(ploidy_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(tab) == 1L && length(names(tab)) == 1L && grepl("\t", names(tab)[[1]], fixed = TRUE)) {
      return(utils::read.delim(ploidy_path, check.names = FALSE, stringsAsFactors = FALSE))
    }
    tab
  } else {
    utils::read.delim(ploidy_path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

collapse_layers_to_ploidy <- function(dists_df, N_UNIT = 22L) {
  if ("passage" %in% names(dists_df)) {
    time_col <- "passage"
  } else if ("day" %in% names(dists_df)) {
    time_col <- "day"
  } else {
    stop("dists_df must contain either a 'passage' column or a 'day' column.")
  }

  tmp <- dplyr::as_tibble(dists_df)
  tmp$time <- tmp[[time_col]]
  tmp %>%
    dplyr::group_by(time, N) %>%
    dplyr::summarise(fraction = sum(fraction), .groups = "drop") %>%
    dplyr::mutate(P = N / N_UNIT)
}

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  stopifnot(length(x) == length(w))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

summarize_ploidy_timecourse <- function(dP) {
  dplyr::as_tibble(dP) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      mean_ploidy = sum(P * fraction),
      q25_ploidy = weighted_quantile(P, fraction, probs = 0.25)[1],
      median_ploidy = weighted_quantile(P, fraction, probs = 0.50)[1],
      q75_ploidy = weighted_quantile(P, fraction, probs = 0.75)[1],
      frac_gt_2_5 = sum(fraction[P > 2.5]),
      frac_gt_3_0 = sum(fraction[P > 3.0]),
      .groups = "drop"
    )
}

read_necrosis_observations <- function(path, use_necrosis_loss = FALSE) {
  empty <- data.frame(
    harvest = character(0),
    obs_necrosis_fraction = numeric(0),
    n_necrosis_obs = integer(0),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(use_necrosis_loss)) return(empty)
  if (is.null(path) || !nzchar(trimws(as.character(path)))) return(empty)
  if (!file.exists(path)) {
    warning("Necrosis mapping CSV not found; necrosis loss will be skipped: ", path)
    return(empty)
  }

  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("dt_harvest", "percent_necrosis")
  missing_cols <- setdiff(required, names(tab))
  if (length(missing_cols) > 0L) {
    stop("Necrosis mapping CSV is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  harvest <- trimws(as.character(tab$dt_harvest))
  percent <- suppressWarnings(as.numeric(tab$percent_necrosis))
  keep <- nzchar(harvest) & is.finite(percent)
  if ("mapping_status" %in% names(tab)) {
    keep <- keep & tolower(trimws(as.character(tab$mapping_status))) == "mapped"
  }
  if (any(keep & (percent < 0 | percent > 100), na.rm = TRUE)) {
    bad <- unique(harvest[keep & (percent < 0 | percent > 100)])
    stop("percent_necrosis must be within [0, 100] for harvest(s): ", paste(bad, collapse = ", "))
  }
  if (!any(keep)) return(empty)

  split_frac <- split(percent[keep] / 100, harvest[keep])
  out <- data.frame(
    harvest = names(split_frac),
    obs_necrosis_fraction = vapply(split_frac, mean, numeric(1), na.rm = TRUE),
    n_necrosis_obs = vapply(split_frac, function(x) sum(is.finite(x)), integer(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out[order(out$harvest), , drop = FALSE]
}

necrosis_observation_for_harvest <- function(necrosis_tab, harvest) {
  empty <- list(obs_necrosis_fraction = NA_real_, n_necrosis_obs = 0L)
  if (!is.data.frame(necrosis_tab) || nrow(necrosis_tab) == 0L) return(empty)
  if (!all(c("harvest", "obs_necrosis_fraction", "n_necrosis_obs") %in% names(necrosis_tab))) {
    return(empty)
  }
  h <- trimws(as.character(harvest))
  if (!nzchar(h)) return(empty)
  idx <- match(h, as.character(necrosis_tab$harvest), nomatch = 0L)
  if (idx == 0L) return(empty)
  list(
    obs_necrosis_fraction = as.numeric(necrosis_tab$obs_necrosis_fraction[[idx]]),
    n_necrosis_obs = as.integer(necrosis_tab$n_necrosis_obs[[idx]])
  )
}

prepare_data <- function(dt_path, ploidy_path, cfg) {
  if (!file.exists(dt_path)) stop("Tumor-burden xlsx not found: ", dt_path)
  if (!file.exists(ploidy_path)) stop("Terminal single-cell table not found: ", ploidy_path)

  dt <- readxl::read_excel(dt_path)
  required <- c("harvest", "Dose", "Day of 1st treatment")
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in tumor-burden sheet: ", paste(missing_cols, collapse = ", "))
  }

  day_cols <- grep("^Day_", names(dt), value = TRUE)
  if (length(day_cols) == 0) stop("No Day_* columns found in tumor-burden sheet.")
  day_vals <- as.numeric(sub("^Day_", "", day_cols))

  day_num_df <- lapply(day_cols, function(col) suppressWarnings(as.numeric(dt[[col]])))
  names(day_num_df) <- day_cols
  day_num_df <- as.data.frame(day_num_df, stringsAsFactors = FALSE)

  pl <- read_terminal_ploidy_table(ploidy_path)
  required_pl_cols <- c("file", "ploidy")
  missing_pl_cols <- setdiff(required_pl_cols, names(pl))
  if (length(missing_pl_cols) > 0L) {
    stop("Terminal single-cell table must include columns: ", paste(required_pl_cols, collapse = ", "))
  }
  if (!("total_chromosomes" %in% names(pl))) {
    pl$total_chromosomes <- NA_real_
  }
  pl$harvest <- sub(".sps.cbs", "", pl$file, fixed = TRUE)
  pl_split <- split(pl[, c("ploidy", "total_chromosomes"), drop = FALSE], pl$harvest)
  necrosis_tab <- read_necrosis_observations(
    path = o2sd_runtime_first_non_null(cfg$necrosis_mapping_csv, NULL),
    use_necrosis_loss = isTRUE(o2sd_runtime_first_non_null(cfg$use_necrosis_loss, FALSE))
  )

  scenarios <- vector("list", nrow(dt))
  keep <- logical(nrow(dt))
  n_endpoint_obs_ploidy <- 0L
  n_endpoint_obs_chr_number <- 0L
  for (i in seq_len(nrow(dt))) {
    h <- as.character(dt$harvest[[i]])
    if (!nzchar(h)) next

    cohort <- if (grepl("2N", h, fixed = TRUE)) "2N" else if (grepl("4N", h, fixed = TRUE)) "4N" else NA_character_
    if (!is.finite(o2sd_runtime_as_num(dt$Dose[[i]], NA_real_))) next
    dose <- o2sd_runtime_as_num(dt$Dose[[i]], NA_real_)
    if (isTRUE(cfg$dose_zero_only) && dose != 0) next
    treat_day <- o2sd_runtime_as_num(dt[["Day of 1st treatment"]][[i]], Inf)
    if (!is.finite(treat_day)) treat_day <- Inf

    y <- as.numeric(day_num_df[i, ])
    idx <- which(is.finite(y))
    if (length(idx) < 2) next

    full_days <- day_vals[idx]
    full_burden <- y[idx]

    if (any(diff(idx) > 1)) next

    if (isTRUE(cfg$truncate_at_treatment)) {
      keep_pre <- full_days <= treat_day
      obs_days <- full_days[keep_pre]
      obs_burden <- full_burden[keep_pre]
    } else {
      obs_days <- full_days
      obs_burden <- full_burden
    }
    if (length(obs_days) < 2) next

    necrosis_obs <- necrosis_observation_for_harvest(necrosis_tab, h)
    obs_necrosis_fraction <- necrosis_obs$obs_necrosis_fraction
    n_necrosis_obs <- necrosis_obs$n_necrosis_obs

    obs_pl <- pl_split[[h]]
    if (is.null(obs_pl)) {
      obs_ploidy_z <- numeric(0)
      obs_chr_number <- numeric(0)
      endpoint_obs_z <- numeric(0)
    } else {
      obs_ploidy_raw <- suppressWarnings(as.numeric(obs_pl$ploidy))
      obs_chr_raw <- suppressWarnings(as.numeric(obs_pl$total_chromosomes))
      obs_ploidy_z <- obs_ploidy_raw[is.finite(obs_ploidy_raw)] * as.numeric(cfg$N_UNIT)
      obs_chr_number <- obs_chr_raw[is.finite(obs_chr_raw)]
      endpoint_obs_z <- if (identical(assert_canonical_start_with_mode(cfg$start_with), "chr_number")) {
        obs_chr_number
      } else {
        obs_ploidy_z
      }
      if (length(obs_ploidy_z) > 0L) n_endpoint_obs_ploidy <- n_endpoint_obs_ploidy + 1L
      if (length(obs_chr_number) > 0L) n_endpoint_obs_chr_number <- n_endpoint_obs_chr_number + 1L
    }

    harvest_day <- max(full_days)
    # Necrosis fractions are harvest endpoint observations.
    sim_end_day <- if (isTRUE(cfg$ploidy_at_harvest) ||
        isTRUE(o2sd_runtime_first_non_null(cfg$use_necrosis_loss, FALSE))) {
      harvest_day
    } else {
      max(obs_days)
    }

    scenarios[[i]] <- list(
      harvest = h,
      cohort = cohort,
      dose = dose,
      treat_day = treat_day,
      obs_days = obs_days,
      obs_burden = obs_burden,
      sim_end_day = sim_end_day,
      harvest_day = harvest_day,
      obs_necrosis_fraction = obs_necrosis_fraction,
      n_necrosis_obs = n_necrosis_obs,
      necrosis_day = harvest_day,
      ploidy_obs_z = obs_ploidy_z,
      chr_number_obs = obs_chr_number,
      endpoint_obs_z = endpoint_obs_z,
      log_init_mult_param = harvest_init_log_param_name(h),
      init_mult_param = harvest_init_natural_param_name(h)
    )
    keep[[i]] <- TRUE
  }

  scenarios <- scenarios[keep]
  if (length(scenarios) == 0) stop("No valid scenarios after preprocessing.")
  paired_only <- isTRUE(o2sd_runtime_first_non_null(cfg$paired_only, TRUE))
  n_before_pair_filter <- length(scenarios)
  n_ploidy_before_pair_filter <- sum(vapply(scenarios, function(s) length(s$ploidy_obs_z) > 0, logical(1)))

  if (paired_only) {
    scenarios <- scenarios[vapply(scenarios, function(s) length(s$endpoint_obs_z) > 0, logical(1))]
    if (length(scenarios) == 0) {
      stop("paired_only=TRUE but no scenarios have both burden and terminal endpoint data.")
    }
    if (length(scenarios) < n_before_pair_filter) {
      message(
        "paired_only=TRUE: retained ", length(scenarios), "/", n_before_pair_filter,
        " scenarios with both burden+endpoint data (dropped ", n_before_pair_filter - length(scenarios), ")."
      )
    }
  }

  if (is.finite(cfg$max_scenarios) && cfg$max_scenarios > 0) {
    scenarios <- scenarios[seq_len(min(length(scenarios), as.integer(cfg$max_scenarios)))]
  }

  matched_ploidy <- sum(vapply(scenarios, function(s) length(s$endpoint_obs_z) > 0, logical(1)))
  message(
    "Prepared scenarios: ", length(scenarios),
    " (with terminal endpoint data: ", matched_ploidy,
    "; paired_only=", paired_only,
    "; pre_pair_filter_ploidy=", n_ploidy_before_pair_filter, "/", n_before_pair_filter, ")"
  )
  message(
    "Endpoint observation mode: start_with=", cfg$start_with,
    " (ploidy_rows=", n_endpoint_obs_ploidy,
    ", chr_number_rows=", n_endpoint_obs_chr_number,
    ", N_UNIT=", cfg$N_UNIT, ")."
  )
  scenarios
}

prepare_cpp_scenarios <- function(scenarios, cfg) {
  n <- length(scenarios)
  cohort_code <- integer(n)
  dose_vec <- numeric(n)
  treat_day_vec <- numeric(n)
  obs_steps_list <- vector("list", n)
  sim_end_step_vec <- integer(n)
  obs_burden_list <- vector("list", n)
  keep_burden_list <- vector("list", n)
  ploidy_z_list <- vector("list", n)
  obs_necrosis_fraction <- rep(NA_real_, n)
  keep_necrosis <- rep(FALSE, n)

  for (i in seq_len(n)) {
    sc <- scenarios[[i]]
    cohort_code[[i]] <- if (identical(sc$cohort, "2N")) 0L else 1L
    dose_vec[[i]] <- as.numeric(sc$dose)
    treat_day_vec[[i]] <- as.numeric(sc$treat_day)
    obs_steps_list[[i]] <- as.integer(round(as.numeric(sc$obs_days) / cfg$DT))
    sim_end_day <- suppressWarnings(as.numeric(sc$sim_end_day))
    sim_end_day <- if (length(sim_end_day) > 0L) sim_end_day[[1]] else NA_real_
    endpoint_days <- suppressWarnings(as.numeric(c(sc$harvest_day, sc$necrosis_day, sc$sim_end_day)))
    endpoint_days <- endpoint_days[is.finite(endpoint_days)]
    endpoint_day <- if (length(endpoint_days) > 0L) endpoint_days[[1]] else NA_real_
    if (isTRUE(o2sd_runtime_first_non_null(cfg$use_necrosis_loss, FALSE)) && is.finite(endpoint_day)) {
      sim_end_day <- endpoint_day
    }
    sim_end_step_vec[[i]] <- as.integer(round(sim_end_day / cfg$DT))
    obs_burden <- as.numeric(sc$obs_burden)
    obs_burden_list[[i]] <- obs_burden

    day_obs <- as.numeric(sc$obs_days)
    keep_day <- rep(TRUE, length(obs_burden))
    if (isTRUE(o2sd_runtime_first_non_null(cfg$burden_exclude_day0, TRUE)) &&
        length(day_obs) == length(obs_burden)) {
      keep_day <- is.finite(day_obs) & (day_obs > 0)
    }
    keep_burden_list[[i]] <- as.logical(keep_day)

    z <- as.numeric(sc$endpoint_obs_z)
    ploidy_z_list[[i]] <- z[is.finite(z)]

    necrosis_fraction <- suppressWarnings(as.numeric(o2sd_runtime_first_non_null(sc$obs_necrosis_fraction, NA_real_)))
    necrosis_day <- suppressWarnings(as.numeric(sc$necrosis_day))
    necrosis_day <- if (length(necrosis_day) > 0L) necrosis_day[[1]] else NA_real_
    if (!is.finite(necrosis_day)) necrosis_day <- endpoint_day
    obs_necrosis_fraction[[i]] <- necrosis_fraction
    keep_necrosis[[i]] <- isTRUE(o2sd_runtime_first_non_null(cfg$use_necrosis_loss, FALSE)) &&
      is.finite(necrosis_fraction) && is.finite(necrosis_day)
  }

  list(
    cohort_code = cohort_code,
    dose = dose_vec,
    treat_day = treat_day_vec,
    obs_steps = obs_steps_list,
    sim_end_step = sim_end_step_vec,
    obs_burden = obs_burden_list,
    keep_burden = keep_burden_list,
    ploidy_z = ploidy_z_list,
    obs_necrosis_fraction = obs_necrosis_fraction,
    keep_necrosis = keep_necrosis
  )
}

build_model_core <- function(run_params = NULL, cfg) {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  R0 <- length(grid_pre)

  init_state_2N <- make_init_state(
    grid_pre = grid_pre,
    ploidy = 2,
    N_UNIT = cfg$N_UNIT,
    total_size = cfg$init_total_size,
    chr_lengths_bp = cfg$chr_lengths_bp
  )
  init_state_4N <- make_init_state(
    grid_pre = grid_pre,
    ploidy = 4,
    N_UNIT = cfg$N_UNIT,
    total_size = cfg$init_total_size,
    chr_lengths_bp = cfg$chr_lengths_bp
  )

  list(
    grid_pre = grid_pre,
    R0 = R0,
    init_state_2N = init_state_2N,
    init_state_4N = init_state_4N
  )
}

# ---- embedded configuration semantics -----------------------------------

# -----------------------------------------------------------------------------
# Function: canonical_ploidy_o2_death_mode
# Purpose: Canonicalize ploidy_O2_death from config/CLI aliases to one mode.
# -----------------------------------------------------------------------------
canonical_ploidy_o2_death_mode <- function(x, default = "diploid_NULL") {
  val <- o2sd_first_non_null(x, default)
  if (is.logical(val) && length(val) > 0L && !is.na(val[[1]])) {
    return(if (isTRUE(val[[1]])) "diploid_NULL" else "uniform")
  }
  s <- tolower(trimws(as.character(val[[1]])))
  if (!nzchar(s)) s <- tolower(trimws(as.character(default[[1]])))
  if (s %in% c("uniform", "false", "f", "0", "no", "n")) return("uniform")
  if (s %in% c("diploid_null", "diploid-null", "diploidnull", "true", "t", "1", "yes", "y")) {
    return("diploid_NULL")
  }
  if (s %in% c("ploidy_related", "ploidy-related", "ploidyrelated")) return("ploidy_related")
  stop(
    "Invalid ploidy_O2_death mode: '", as.character(val[[1]]),
    "'. Allowed values: uniform, diploid_NULL, ploidy_related."
  )
}

# -----------------------------------------------------------------------------
# Function: canonical_start_with_mode
# Purpose: Canonicalize endpoint observation mode for fit/viz dispatch.
# -----------------------------------------------------------------------------
canonical_start_with_mode <- function(x, default = "ploidy") {
  val <- o2sd_first_non_null(x, default)
  s <- tolower(trimws(as.character(val[[1]])))
  if (!nzchar(s)) s <- tolower(trimws(as.character(default[[1]])))
  if (s %in% c("ploidy", "p")) return("ploidy")
  if (s %in% c("chr_number", "chr-number", "chrnumber", "chromosome_number", "chromosome-number", "chromosomenumber", "chromosome", "chromosomes", "n")) {
    return("chr_number")
  }
  stop(
    "Invalid start_with mode: '", as.character(val[[1]]),
    "'. Allowed values: ploidy, chr_number."
  )
}

# -----------------------------------------------------------------------------
# Function: assert_canonical_start_with_mode
# Purpose: Enforce that runtime start_with is already canonical.
# -----------------------------------------------------------------------------
assert_canonical_start_with_mode <- function(x) {
  val <- o2sd_first_non_null(x, NA_character_)
  s <- trimws(as.character(val[[1]]))
  if (!nzchar(s)) {
    stop("start_with must be provided as one of: ploidy, chr_number.")
  }
  if (identical(s, "ploidy") || identical(s, "chr_number")) {
    return(s)
  }
  stop(
    "start_with must already be canonical before runtime dispatch. ",
    "Allowed values: ploidy, chr_number; got '", s, "'."
  )
}

# -----------------------------------------------------------------------------
# Function: assert_canonical_ploidy_o2_death_mode
# Purpose: Enforce that runtime ploidy_O2_death is already canonical.
# -----------------------------------------------------------------------------
assert_canonical_ploidy_o2_death_mode <- function(x) {
  val <- o2sd_first_non_null(x, NA_character_)
  s <- trimws(as.character(val[[1]]))
  if (!nzchar(s)) {
    stop("ploidy_O2_death must be provided as one of: uniform, diploid_NULL, ploidy_related.")
  }
  if (identical(s, "uniform") || identical(s, "diploid_NULL") || identical(s, "ploidy_related")) {
    return(s)
  }
  stop(
    "ploidy_O2_death must already be canonical before runtime dispatch. ",
    "Allowed values: uniform, diploid_NULL, ploidy_related; got '", s, "'."
  )
}

# -----------------------------------------------------------------------------
# Function: filter_family_specific_run_params_for_output_common
# Purpose: Remove inactive family-specific parameters from natural-scale output
#   tables so outputs reflect the active model family only.
# -----------------------------------------------------------------------------
filter_family_specific_run_params_for_output_common <- function(run_params) {
  rp <- as.list(run_params)
  rp
}

# -----------------------------------------------------------------------------
# Function: filter_fit_summary_metrics_for_output_common
# Purpose: Remove inactive family-specific parameter summary rows from
#   fit_summary.tsv outputs.
# -----------------------------------------------------------------------------
filter_fit_summary_metrics_for_output_common <- function(summary_df) {
  if (!is.data.frame(summary_df) || !"metric" %in% names(summary_df)) {
    return(summary_df)
  }

  drop_metrics <- character(0)
  summary_df[!(summary_df$metric %in% unique(drop_metrics)), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# Function: read_param_table_prototype_slot_common
# Purpose: Read one natural-scale prototype slot from parameter_table.csv.
# -----------------------------------------------------------------------------
read_param_table_prototype_slot_common <- function(path, param_prototype, slot = c("init", "lower", "upper")) {
  slot <- match.arg(slot)
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NA_real_)
  tab <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab)) return(NA_real_)

  if (all(c("param_symbol", "init_value", "lower_bound", "upper_bound") %in% names(tab))) {
    col_map <- c(
      init = "init_value",
      lower = "lower_bound",
      upper = "upper_bound"
    )
    proto <- trimws(as.character(tab$param_symbol))
  } else {
    col_map <- c(
      init = "prototype_init_value",
      lower = "prototype_lower_bound",
      upper = "prototype_upper_bound"
    )
    if (!all(c("param_prototype", unname(col_map)) %in% names(tab))) return(NA_real_)
    proto <- trimws(as.character(tab$param_prototype))
  }

  idx <- match(param_prototype, proto)
  if (is.na(idx)) return(NA_real_)
  suppressWarnings(as.numeric(tab[[col_map[[slot]]]][idx]))
}

# -----------------------------------------------------------------------------
# Function: read_o2_S0_natural_upper_bound_common
# Purpose: Read natural-scale o2_S0 upper bound from parameter_table.csv.
# -----------------------------------------------------------------------------
read_o2_S0_natural_upper_bound_common <- function(path, fallback = 5.0) {
  ub <- read_param_table_prototype_slot_common(path, "o2_S0", slot = "upper")
  if (!is.finite(ub) || ub <= 0) ub <- as.numeric(fallback)
  if (!is.finite(ub) || ub <= 0) ub <- 5.0
  ub
}

# -----------------------------------------------------------------------------
# Function: default_o2_parameter_table_path_common
# Purpose: Resolve the default natural-scale parameter table path.
# -----------------------------------------------------------------------------
default_o2_parameter_table_path_common <- function(script_dir, must_exist = FALSE) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  file_name <- "parameter_table_O2.csv"
  path <- normalizePath(
    file.path(workflow_root, "..", "..", "data", "O2_supply_demand", file_name),
    mustWork = FALSE
  )
  if (isTRUE(must_exist) && !file.exists(path)) {
    stop(
      "Default parameter table not found: ",
      path
    )
  }
  path
}

# -----------------------------------------------------------------------------
# Function: normalize_sim_cfg_common
# Purpose: Canonicalize shared simulation config semantics for fit/viz.
# -----------------------------------------------------------------------------
normalize_sim_cfg_common <- function(cfg, context = c("fit", "viz")) {
  context <- match.arg(context)
  if (is.null(cfg)) cfg <- list()

  cfg$N_UNIT <- as.integer(o2sd_first_non_null(cfg$N_UNIT, 22L))
  cfg$N_MIN <- as.integer(o2sd_first_non_null(cfg$N_MIN, 22L))
  cfg$N_MAX <- as.integer(o2sd_first_non_null(cfg$N_MAX, 154L))
  cfg$DT <- as.numeric(o2sd_first_non_null(cfg$DT, 0.5))

  cfg$o2_S0_upper_bound <- as.numeric(o2sd_first_non_null(
    cfg$o2_S0_upper_bound,
    read_o2_S0_natural_upper_bound_common(cfg$parameter_table, fallback = 5.0)
  ))
  if (!is.finite(cfg$o2_S0_upper_bound) || cfg$o2_S0_upper_bound <= 0) {
    stop(context, "_config o2_S0_upper_bound must be > 0.")
  }

  cfg$o2_min <- as.numeric(o2sd_first_non_null(cfg$o2_min, 0.0))
  if (!is.finite(cfg$o2_min) || cfg$o2_min < 0) cfg$o2_min <- 0.0
  cfg$o2_min <- min(max(cfg$o2_min, 0), cfg$o2_S0_upper_bound)

  cfg$init_total_size <- as.numeric(o2sd_first_non_null(cfg$init_total_size, 1e6))
  cfg$o2_Nref <- as.numeric(o2sd_first_non_null(cfg$o2_Nref, cfg$init_total_size, 1e6))
  if (!is.finite(cfg$o2_Nref) || cfg$o2_Nref <= 0) cfg$o2_Nref <- 1e6

  cfg$tau_O2_init <- as.numeric(o2sd_first_non_null(cfg$tau_O2_init, 2.0))
  cfg$tau_O2 <- as.numeric(o2sd_first_non_null(cfg$tau_O2, cfg$tau_O2_init))
  if (!is.finite(cfg$tau_O2) || cfg$tau_O2 <= 0) cfg$tau_O2 <- cfg$tau_O2_init

  cfg$K <- as.numeric(o2sd_first_non_null(cfg$K, 1e12))
  if (!is.finite(cfg$K) || cfg$K <= 0) cfg$K <- 1e12
  cfg$crowding <- as.character(o2sd_first_non_null(cfg$crowding, "logistic"))
  cfg$Crowding <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$Crowding, TRUE), TRUE)
  if (!cfg$crowding %in% c("logistic", "gompertz")) {
    stop(context, "_config crowding must be logistic or gompertz.")
  }

  cfg$dose_ref <- as.numeric(o2sd_first_non_null(cfg$dose_ref, 30))
  cfg$tx_mult_min <- as.numeric(o2sd_first_non_null(cfg$tx_mult_min, 0.05))
  cfg$min_pop <- as.numeric(o2sd_first_non_null(cfg$min_pop, 1e-12))
  cfg$alpha_o2_init <- as.numeric(o2sd_first_non_null(cfg$alpha_o2_init, 0.5))
  cfg$gamma_growth_init <- as.numeric(o2sd_first_non_null(cfg$gamma_growth_init, 2.0))
  cfg$mu_hp_init <- as.numeric(o2sd_first_non_null(cfg$mu_hp_init, 1e-3))
  cfg$gamma_mu_init <- as.numeric(o2sd_first_non_null(cfg$gamma_mu_init, 1.0))
  cfg$o2_crit_init <- as.numeric(o2sd_first_non_null(cfg$o2_crit_init, 1.0))
  cfg$n_O_init <- as.numeric(o2sd_first_non_null(cfg$n_O_init, 1.0))
  cfg$k_clear_init <- as.numeric(o2sd_first_non_null(cfg$k_clear_init, 1e-3))
  cfg$p_wgd_init <- as.numeric(o2sd_first_non_null(cfg$p_wgd_init, 1e-4))
  cfg$harvest_init_multiplier <- o2sd_as_bool_scalar(
    o2sd_first_non_null(cfg$harvest_init_multiplier, FALSE),
    FALSE
  )
  cfg$use_necrosis_loss <- o2sd_as_bool_scalar(
    o2sd_first_non_null(cfg$use_necrosis_loss, FALSE),
    FALSE
  )
  cfg$necrosis_mapping_csv <- o2sd_first_non_null(cfg$necrosis_mapping_csv, NULL)
  cfg$sigma_necrosis_logit <- as.numeric(o2sd_first_non_null(cfg$sigma_necrosis_logit, 0.75))
  cfg$lambda_necrosis <- as.numeric(o2sd_first_non_null(cfg$lambda_necrosis, 1.0))
  cfg$necrosis_fraction_eps <- as.numeric(o2sd_first_non_null(cfg$necrosis_fraction_eps, 1e-4))
  cfg$prior_center_log_init_mult <- as.numeric(o2sd_first_non_null(cfg$prior_center_log_init_mult, 0.0))
  cfg$prior_sd_log_init_mult <- as.numeric(o2sd_first_non_null(cfg$prior_sd_log_init_mult, 0.35))
  cfg$log_init_mult_lower <- as.numeric(o2sd_first_non_null(cfg$log_init_mult_lower, -1.0))
  cfg$log_init_mult_upper <- as.numeric(o2sd_first_non_null(cfg$log_init_mult_upper, 1.0))
  cfg$ploidy_O2_death <- canonical_ploidy_o2_death_mode(
    o2sd_first_non_null(cfg$ploidy_O2_death, "diploid_NULL"),
    default = "diploid_NULL"
  )
  cfg$buffer_smax_init <- as.numeric(o2sd_first_non_null(cfg$buffer_smax_init, 0.8))
  cfg$buffer_beta_init <- as.numeric(o2sd_first_non_null(cfg$buffer_beta_init, 1.0))
  cfg$buffer_n_exp_init <- as.numeric(o2sd_first_non_null(cfg$buffer_n_exp_init, 1.0))
  cfg$prior_center_buffer_smax <- as.numeric(o2sd_first_non_null(cfg$prior_center_buffer_smax, cfg$buffer_smax_init, 0.8))
  cfg$prior_sd_buffer_smax <- as.numeric(o2sd_first_non_null(cfg$prior_sd_buffer_smax, 0.25))
  cfg$prior_center_log10_buffer_beta <- as.numeric(o2sd_first_non_null(cfg$prior_center_log10_buffer_beta, log10(max(cfg$buffer_beta_init, 1e-8))))
  cfg$prior_sd_log10_buffer_beta <- as.numeric(o2sd_first_non_null(cfg$prior_sd_log10_buffer_beta, 0.75))
  cfg$prior_center_log10_buffer_n_exp <- as.numeric(o2sd_first_non_null(cfg$prior_center_log10_buffer_n_exp, log10(max(cfg$buffer_n_exp_init, 1e-8))))
  cfg$prior_sd_log10_buffer_n_exp <- as.numeric(o2sd_first_non_null(cfg$prior_sd_log10_buffer_n_exp, 0.75))
  cfg$start_with <- canonical_start_with_mode(
    o2sd_first_non_null(cfg$start_with, "ploidy"),
    default = "ploidy"
  )

  if (!is.finite(cfg$mu_hp_init) || cfg$mu_hp_init <= 0) cfg$mu_hp_init <- 1e-3
  if (!is.finite(cfg$gamma_mu_init) || cfg$gamma_mu_init <= 0) cfg$gamma_mu_init <- 1.0
  if (!is.finite(cfg$o2_crit_init) || cfg$o2_crit_init < 0) cfg$o2_crit_init <- 1.0
  if (!is.finite(cfg$n_O_init) || cfg$n_O_init < 0) cfg$n_O_init <- 1.0
  if (!is.finite(cfg$k_clear_init) || cfg$k_clear_init <= 0) cfg$k_clear_init <- 1e-3
  if (!is.finite(cfg$p_wgd_init) || cfg$p_wgd_init <= 0) cfg$p_wgd_init <- 1e-4
  if (!is.finite(cfg$buffer_smax_init)) cfg$buffer_smax_init <- 0.8
  cfg$buffer_smax_init <- min(max(cfg$buffer_smax_init, 0), 1)
  if (!is.finite(cfg$buffer_beta_init) || cfg$buffer_beta_init < 0) cfg$buffer_beta_init <- 1.0
  if (!is.finite(cfg$buffer_n_exp_init) || cfg$buffer_n_exp_init < 0) cfg$buffer_n_exp_init <- 1.0
  if (!is.finite(cfg$prior_center_buffer_smax)) cfg$prior_center_buffer_smax <- cfg$buffer_smax_init
  cfg$prior_center_buffer_smax <- min(max(cfg$prior_center_buffer_smax, 0), 1)
  if (!is.finite(cfg$sigma_necrosis_logit) || cfg$sigma_necrosis_logit <= 0) cfg$sigma_necrosis_logit <- 0.75
  if (!is.finite(cfg$lambda_necrosis) || cfg$lambda_necrosis < 0) cfg$lambda_necrosis <- 1.0
  if (!is.finite(cfg$necrosis_fraction_eps) || cfg$necrosis_fraction_eps <= 0 || cfg$necrosis_fraction_eps >= 0.5) {
    cfg$necrosis_fraction_eps <- 1e-4
  }
  if (!is.finite(cfg$prior_sd_buffer_smax) || cfg$prior_sd_buffer_smax <= 0) cfg$prior_sd_buffer_smax <- 0.25
  if (!is.finite(cfg$prior_center_log10_buffer_beta)) cfg$prior_center_log10_buffer_beta <- log10(max(cfg$buffer_beta_init, 1e-8))
  if (!is.finite(cfg$prior_sd_log10_buffer_beta) || cfg$prior_sd_log10_buffer_beta <= 0) cfg$prior_sd_log10_buffer_beta <- 0.75
  if (!is.finite(cfg$prior_center_log10_buffer_n_exp)) cfg$prior_center_log10_buffer_n_exp <- log10(max(cfg$buffer_n_exp_init, 1e-8))
  if (!is.finite(cfg$prior_sd_log10_buffer_n_exp) || cfg$prior_sd_log10_buffer_n_exp <= 0) cfg$prior_sd_log10_buffer_n_exp <- 0.75
  if (!is.finite(cfg$prior_center_log_init_mult)) cfg$prior_center_log_init_mult <- 0.0
  if (!is.finite(cfg$prior_sd_log_init_mult) || cfg$prior_sd_log_init_mult <= 0) cfg$prior_sd_log_init_mult <- 0.35
  if (!is.finite(cfg$log_init_mult_lower)) cfg$log_init_mult_lower <- -1.0
  if (!is.finite(cfg$log_init_mult_upper)) cfg$log_init_mult_upper <- 1.0
  if (cfg$log_init_mult_upper < cfg$log_init_mult_lower) {
    tmp <- cfg$log_init_mult_lower
    cfg$log_init_mult_lower <- cfg$log_init_mult_upper
    cfg$log_init_mult_upper <- tmp
  }

  cfg$dose_zero_only <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$dose_zero_only, TRUE), TRUE)
  cfg$fit_treatment <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$fit_treatment, FALSE), FALSE)
  cfg$max_scenarios <- as.numeric(o2sd_first_non_null(cfg$max_scenarios, Inf))
  cfg$o2_burden_feedback <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$o2_burden_feedback, TRUE), TRUE)
  cfg$O2_growth <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$O2_growth, TRUE), TRUE)
  cfg$o2_cache_profile <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$o2_cache_profile, FALSE), FALSE)

  if (is.null(cfg$truncate_at_treatment)) {
    cfg$truncate_at_treatment <- o2sd_as_bool_scalar(o2sd_first_non_null(cfg$pretreat_only, FALSE), FALSE)
  }
  cfg$truncate_at_treatment <- o2sd_as_bool_scalar(cfg$truncate_at_treatment, FALSE)

  if (is.null(cfg$ploidy_at_harvest)) cfg$ploidy_at_harvest <- TRUE
  cfg$ploidy_at_harvest <- o2sd_as_bool_scalar(cfg$ploidy_at_harvest, TRUE)

  cfg
}

# -----------------------------------------------------------------------------
# Function: normalize_run_params_common
# Purpose: Fill shared run_params defaults/canonical values after reconstruction.
# -----------------------------------------------------------------------------
normalize_run_params_common <- function(run_params, cfg = NULL) {
  if (is.null(run_params)) run_params <- list()
  if (is.null(cfg)) cfg <- list()

  run_params$p_mis_base <- as.numeric(o2sd_first_non_null(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5))
  run_params$harvest_init_multiplier <- o2sd_as_bool_scalar(
    o2sd_first_non_null(run_params$harvest_init_multiplier, cfg$harvest_init_multiplier, FALSE),
    FALSE
  )
  run_params$ploidy_O2_death <- canonical_ploidy_o2_death_mode(
    o2sd_first_non_null(run_params$ploidy_O2_death, cfg$ploidy_O2_death, "diploid_NULL"),
    default = "diploid_NULL"
  )
  run_params$start_with <- canonical_start_with_mode(
    o2sd_first_non_null(run_params$start_with, cfg$start_with, "ploidy"),
    default = "ploidy"
  )
  run_params$buffer_smax <- as.numeric(o2sd_first_non_null(run_params$buffer_smax, cfg$buffer_smax_init, 0.8))
  if (!is.finite(run_params$buffer_smax)) run_params$buffer_smax <- 0.8
  run_params$buffer_smax <- min(max(run_params$buffer_smax, 0), 1)
  run_params$buffer_beta <- as.numeric(o2sd_first_non_null(run_params$buffer_beta, cfg$buffer_beta_init, 1.0))
  if (!is.finite(run_params$buffer_beta) || run_params$buffer_beta < 0) run_params$buffer_beta <- 1.0
  run_params$buffer_n_exp <- as.numeric(o2sd_first_non_null(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1.0))
  if (!is.finite(run_params$buffer_n_exp) || run_params$buffer_n_exp < 0) run_params$buffer_n_exp <- 1.0

  run_params$o2_S0_upper_bound <- as.numeric(o2sd_first_non_null(
    run_params$o2_S0_upper_bound,
    cfg$o2_S0_upper_bound,
    read_o2_S0_natural_upper_bound_common(cfg$parameter_table, fallback = 5.0)
  ))
  if (!is.finite(run_params$o2_S0_upper_bound) || run_params$o2_S0_upper_bound <= 0) {
    run_params$o2_S0_upper_bound <- 5.0
  }

  run_params$o2_min <- as.numeric(o2sd_first_non_null(run_params$o2_min, cfg$o2_min, 0.0))
  if (!is.finite(run_params$o2_min) || run_params$o2_min < 0) run_params$o2_min <- 0.0
  run_params$o2_min <- min(max(run_params$o2_min, 0), run_params$o2_S0_upper_bound)

  run_params$o2_Nref <- as.numeric(o2sd_first_non_null(run_params$o2_Nref, cfg$o2_Nref, cfg$init_total_size, 1e6))
  if (!is.finite(run_params$o2_Nref) || run_params$o2_Nref <= 0) run_params$o2_Nref <- 1e6

  fit_treatment_use <- isTRUE(o2sd_first_non_null(cfg$fit_treatment, FALSE))
  if (fit_treatment_use) {
    if (is.null(run_params$alpha) || !is.finite(as.numeric(run_params$alpha))) {
      stop("run_params$alpha must be present when fit_treatment=TRUE.")
    }
    if (is.null(run_params$gamma) || !is.finite(as.numeric(run_params$gamma))) {
      stop("run_params$gamma must be present when fit_treatment=TRUE.")
    }
    run_params$alpha <- as.numeric(run_params$alpha)
    run_params$gamma <- as.numeric(run_params$gamma)
  } else {
    run_params$alpha <- as.numeric(o2sd_first_non_null(run_params$alpha, 0))
    run_params$gamma <- as.numeric(o2sd_first_non_null(run_params$gamma, 1))
  }

  tau_use <- as.numeric(run_params$tau_O2)
  if (!is.finite(tau_use) || tau_use <= 0) {
    tau_use <- as.numeric(o2sd_first_non_null(cfg$tau_O2, cfg$tau_O2_init, 2.0))
  }
  if (!is.finite(tau_use) || tau_use <= 0) tau_use <- 2.0
  run_params$tau_O2 <- tau_use

  run_params$p_wgd <- as.numeric(o2sd_first_non_null(run_params$p_wgd, cfg$p_wgd_init, 1e-4))
  if (!is.finite(run_params$p_wgd) || run_params$p_wgd < 0) run_params$p_wgd <- 0.0

  run_params
}

# ---- embedded fixed-O2 shared helpers -----------------------------------
#!/usr/bin/env Rscript

# Shared fixed-O2 fitting-parameter and model-loading helpers used by
# simulation/fix_o2_simulation.R and analyses that source it.
suppressPackageStartupMessages({
  if (requireNamespace("Matrix", quietly = TRUE)) library(Matrix)
})

o2ipa_null_coalesce <- function(x, y) {
    if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x)))
        y
    else x
}

o2ipa_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        arg <- args[[i]]
        if (!grepl("^--", arg)) {
            i <- i + 1L
            next
        }
        kv <- sub("^--", "", arg)
        eq <- regexpr("=", kv, fixed = TRUE)
        if (eq > 0L) {
            key <- substr(kv, 1L, eq - 1L)
            val <- substr(kv, eq + 1L, nchar(kv))
            out[[key]] <- val
            i <- i + 1L
        } else {
            key <- kv
            if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
                out[[key]] <- args[[i + 1L]]
                i <- i + 2L
            } else {
                out[[key]] <- TRUE
                i <- i + 1L
            }
        }
    }
    out
}

o2ipa_as_chr <- function(x, default = "") {
    val <- o2ipa_null_coalesce(x, default)
    val <- as.character(val[[1]])
    if (!nzchar(val))
        default
    else val
}

o2ipa_as_num <- function(x, default = NA_real_) {
    val <- suppressWarnings(as.numeric(o2ipa_null_coalesce(x, default)[[1]]))
    if (is.finite(val))
        val
    else default
}

o2ipa_as_int <- function(x, default = NA_integer_) {
    val <- suppressWarnings(as.integer(o2ipa_null_coalesce(x, default)[[1]]))
    if (!is.na(val))
        val
    else default
}

o2ipa_as_bool <- function(x, default = FALSE) {
    if (is.null(x) || !length(x) || is.na(x[[1]]))
        return(isTRUE(default))
    if (is.logical(x[[1]]))
        return(isTRUE(x[[1]]))
    tolower(trimws(as.character(x[[1]]))) %in% c("1", "true", "t", "yes", "y", "on")
}

o2ipa_split_csv <- function(x, default = character()) {
    txt <- trimws(o2ipa_as_chr(x, paste(default, collapse = ",")))
    if (!nzchar(txt))
        return(default)
    vals <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
    vals[nzchar(vals)]
}

o2ipa_script_dir <- function() {
    frame_files <- Filter(
        nzchar,
        vapply(sys.frames(), function(env) {
            ofile <- env$ofile
            if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
        }, character(1))
    )
    if (length(frame_files)) {
        return(dirname(frame_files[[length(frame_files)]]))
    }
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
    }
    normalizePath(getwd(), mustWork = FALSE)
}

o2ipa_find_workflow_root <- function(path = o2ipa_script_dir()) {
    explicit <- trimws(Sys.getenv("FIGURE_MODEL_CODE_ROOT", unset = ""))
    if (nzchar(explicit) &&
        file.exists(file.path(explicit, "util", "o2_supply_demand_map_shared.R")) &&
        file.exists(file.path(explicit, "model", "model_O2_supply_demand_MAP.R"))) {
        return(normalizePath(explicit, mustWork = TRUE))
    }
    cur <- normalizePath(path, mustWork = FALSE)
    if (file.exists(cur) && !dir.exists(cur))
        cur <- dirname(cur)
    for (i in seq_len(8L)) {
        if (file.exists(file.path(cur, "util", "o2_supply_demand_map_shared.R")) && file.exists(file.path(cur, "model", "model_O2_supply_demand_MAP.R"))) {
            return(normalizePath(cur, mustWork = FALSE))
        }
        repo_candidate <- file.path(cur, "oxygen", "code", "O2_supply_demand_MAP")
        if (file.exists(file.path(repo_candidate, "util", "o2_supply_demand_map_shared.R")) && file.exists(file.path(repo_candidate, "model", "model_O2_supply_demand_MAP.R"))) {
            return(normalizePath(repo_candidate, mustWork = FALSE))
        }
        parent <- dirname(cur)
        if (identical(parent, cur))
            break
        cur <- parent
    }
    normalizePath(file.path(path, "..", ".."), mustWork = FALSE)
}

o2ipa_workflow_root <- function(script_dir = o2ipa_script_dir()) {
    o2ipa_find_workflow_root(script_dir)
}

o2ipa_repo_root <- function(script_dir = o2ipa_script_dir()) {
    normalizePath(file.path(o2ipa_workflow_root(script_dir), "..", "..", ".."), mustWork = FALSE)
}

o2ipa_read_tsv <- function(path) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

o2ipa_seed_number <- function(seed_id) {
    suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

o2ipa_norm_seed <- function(x) {
    s <- as.character(x)
    s <- trimws(s)
    out <- ifelse(grepl("^seed[0-9]+$", s), s, ifelse(grepl("^[0-9]+$", s), paste0("seed", s), s))
    out
}

o2ipa_order_seeds <- function(seed_id) {
    n <- o2ipa_seed_number(seed_id)
    order(ifelse(is.na(n), Inf, n), seed_id)
}

o2ipa_target_params <- function() {
    c("O2_crit", "alpha_o2", "gamma_growth", "mu_hp", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta", "buffer_n_exp",
        "n_O", "lam_max", "p_mis_base", "p_wgd", "gamma_mu", "o2_S0", "kappa_O", "eta_o2", "rho_2N")
}

o2ipa_param_aliases <- function() {
    list(O2_crit = c("O2_crit", "o2_crit", "O2crit", "o2crit"), alpha_o2 = c("alpha_o2", "alpha_O2"), gamma_growth = c("gamma_growth"),
        mu_hp = c("mu_hp"), p_misseg = c("p_misseg", "p_mis", "p_missegregation"), k_o_mis = c("k_o_mis", "ko_mis", "k_O_mis"),
        buffer_smax = c("buffer_smax", "buffer_s_max"), buffer_beta = c("buffer_beta"), buffer_n_exp = c("buffer_n_exp",
            "buffer_n"), n_O = c("n_O", "n_o", "nO"), lam_max = c("lam_max", "lambda_max"), p_mis_base = c("p_mis_base",
            "p_misseg_base"), p_wgd = c("p_wgd", "wgd_prob"), gamma_mu = c("gamma_mu"), o2_S0 = c("o2_S0", "O2_S0", "S0_o2"),
        kappa_O = c("kappa_O", "kappa_o"), eta_o2 = c("eta_o2", "eta_O2"), rho_2N = c("rho_2N", "rho2N"))
}

o2ipa_optimizer_transform <- function(parameter) {
    log10_params <- c("O2_crit", "alpha_o2", "mu_hp", "p_misseg", "k_o_mis", "buffer_beta", "buffer_n_exp", "lam_max", "p_mis_base",
        "p_wgd", "o2_S0", "kappa_O", "eta_o2", "rho_2N")
    if (parameter %in% log10_params)
        "log10"
    else "identity"
}

o2ipa_transform_parameter_value <- function(parameter, value, epsilon = 1e-12) {
    v <- as.numeric(value)
    if (!is.finite(v))
        return(NA_real_)
    tr <- o2ipa_optimizer_transform(parameter)
    if (identical(tr, "log10")) {
        if (v <= 0)
            return(NA_real_)
        return(log10(v))
    }
    if (identical(tr, "logit")) {
        vv <- min(max(v, epsilon), 1 - epsilon)
        return(stats::qlogis(vv))
    }
    v
}

o2ipa_find_extra <- function(run_dir, file) {
    path <- file.path(run_dir, "extra_results", file)
    if (file.exists(path))
        path
    else NA_character_
}

o2ipa_read_extra_summary <- function(run_dir) {
    path <- o2ipa_find_extra(run_dir, "seed_summary.tsv")
    if (is.na(path))
        return(NULL)
    tab <- o2ipa_read_tsv(path)
    if (!"seed" %in% names(tab))
        return(NULL)
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    tab
}

o2ipa_read_param_matrix <- function(run_dir) {
    parent <- dirname(normalizePath(run_dir, mustWork = FALSE))
    base <- basename(normalizePath(run_dir, mustWork = FALSE))
    candidates <- c(file.path(run_dir, "param_matrix_with_seed.tsv"), file.path(run_dir, "parameter_matrix_with_seed.tsv"),
        file.path(parent, paste0(base, "_param_matrix_with_seed.tsv")), file.path(parent, paste0(base, "_param_matrix.tsv")))
    candidates <- candidates[file.exists(candidates)]
    if (!length(candidates))
        return(NULL)
    tab <- o2ipa_read_tsv(candidates[[1]])
    if (!"seed" %in% names(tab))
        tab$seed <- seq_len(nrow(tab))
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    attr(tab, "source_file") <- candidates[[1]]
    tab
}

o2ipa_discover_seeds <- function(run_dir) {
    dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[grepl("^seed[^/]*[0-9]+$", basename(dirs))]
    seed_ids <- o2ipa_norm_seed(basename(dirs))
    seed_dir_map <- setNames(normalizePath(dirs, mustWork = FALSE), seed_ids)
    summary_tab <- o2ipa_read_extra_summary(run_dir)
    matrix_tab <- o2ipa_read_param_matrix(run_dir)
    all_ids <- unique(c(seed_ids, if (!is.null(summary_tab)) summary_tab$seed_id, if (!is.null(matrix_tab)) matrix_tab$seed_id))
    all_ids <- all_ids[o2ipa_order_seeds(all_ids)]
    data.frame(seed_id = all_ids, seed_dir = unname(seed_dir_map[all_ids]), stringsAsFactors = FALSE)
}

o2ipa_metric_map <- function(path) {
    if (is.null(path) || is.na(path) || !file.exists(path))
        return(list())
    tab <- o2ipa_read_tsv(path)
    if (!all(c("metric", "value") %in% names(tab)))
        return(list())
    vals <- as.list(tab$value)
    names(vals) <- tab$metric
    vals
}

o2ipa_num_from_map <- function(map, key) {
    if (is.null(map[[key]]))
        return(NA_real_)
    suppressWarnings(as.numeric(map[[key]]))
}

o2ipa_chr_from_map <- function(map, key) {
    if (is.null(map[[key]]))
        return(NA_character_)
    as.character(map[[key]])
}

o2ipa_read_best_params <- function(path) {
    if (!file.exists(path))
        return(list(values = setNames(numeric(0), character(0)), aliases = data.frame()))
    tab <- o2ipa_read_tsv(path)
    if (!all(c("parameter", "value") %in% names(tab))) {
        stop("best_params.tsv must contain parameter and value columns: ", path)
    }
    vals <- suppressWarnings(as.numeric(tab$value))
    names(vals) <- trimws(as.character(tab$parameter))
    list(values = vals, aliases = data.frame())
}

o2ipa_extract_param <- function(seed_id, parameter, best_vals, summary_row, matrix_row) {
    aliases <- o2ipa_param_aliases()[[parameter]]
    sources <- list()
    if (length(best_vals)) {
        sources$best_params <- best_vals
    }
    if (!is.null(summary_row) && nrow(summary_row) == 1L) {
        vals <- numeric(0)
        for (a in aliases) {
            col <- paste0("value__", a)
            if (col %in% names(summary_row))
                vals[a] <- suppressWarnings(as.numeric(summary_row[[col]][[1]]))
        }
        sources$seed_summary_value_cols <- vals
    }
    if (!is.null(matrix_row) && nrow(matrix_row) == 1L) {
        vals <- numeric(0)
        for (a in aliases) {
            for (col in c(paste0("final_", a), a)) {
                if (col %in% names(matrix_row))
                  vals[a] <- suppressWarnings(as.numeric(matrix_row[[col]][[1]]))
            }
        }
        sources$param_matrix <- vals
    }
    for (src in names(sources)) {
        vals <- sources[[src]]
        for (a in aliases) {
            if (a %in% names(vals)) {
                v <- suppressWarnings(as.numeric(vals[[a]]))
                if (is.finite(v)) {
                  return(list(value = v, source = src, alias = a))
                }
            }
        }
    }
    list(value = NA_real_, source = NA_character_, alias = NA_character_)
}

o2ipa_extract_all_params <- function(seed_id, seed_dir, summary_tab, matrix_tab) {
    best_path <- if (!is.na(seed_dir) && nzchar(seed_dir))
        file.path(seed_dir, "best_params.tsv")
    else NA_character_
    best_vals <- if (!is.na(best_path) && file.exists(best_path))
        o2ipa_read_best_params(best_path)$values
    else setNames(numeric(0), character(0))
    summary_row <- if (!is.null(summary_tab))
        summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE]
    else NULL
    if (!is.null(summary_row) && nrow(summary_row) > 1L)
        summary_row <- summary_row[1, , drop = FALSE]
    matrix_row <- if (!is.null(matrix_tab))
        matrix_tab[matrix_tab$seed_id == seed_id, , drop = FALSE]
    else NULL
    if (!is.null(matrix_row) && nrow(matrix_row) > 1L)
        matrix_row <- matrix_row[1, , drop = FALSE]
    params <- o2ipa_target_params()
    rows <- lapply(params, function(p) {
        hit <- o2ipa_extract_param(seed_id, p, best_vals, summary_row, matrix_row)
        data.frame(seed_id = seed_id, parameter = p, value = hit$value, parameter_source = hit$source, matched_alias = hit$alias,
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

o2ipa_choose_objective <- function(summary_vals, objective_source = "auto") {
    raw <- NA_real_
    b <- o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw")
    p <- o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw")
    if (is.finite(b) && is.finite(p))
        raw <- b + p
    data <- o2ipa_num_from_map(summary_vals, "objective_data")
    total <- o2ipa_num_from_map(summary_vals, "objective_total")
    if (!is.finite(total))
        total <- o2ipa_num_from_map(summary_vals, "objective")
    source <- match.arg(objective_source, c("auto", "raw_likelihood", "data", "total"))
    if (identical(source, "raw_likelihood"))
        return(list(value = raw, source = "raw_likelihood"))
    if (identical(source, "data"))
        return(list(value = data, source = "objective_data"))
    if (identical(source, "total"))
        return(list(value = total, source = if (is.finite(o2ipa_num_from_map(summary_vals, "objective_total"))) "objective_total" else "objective"))
    if (is.finite(raw))
        return(list(value = raw, source = "raw_likelihood"))
    if (is.finite(data))
        return(list(value = data, source = "objective_data"))
    list(value = total, source = if (is.finite(total)) "objective_total_or_objective" else NA_character_)
}

o2ipa_collect_seed_inputs <- function(run_dir, objective_source = "auto") {
    manifest0 <- o2ipa_discover_seeds(run_dir)
    summary_tab <- o2ipa_read_extra_summary(run_dir)
    matrix_tab <- o2ipa_read_param_matrix(run_dir)
    boundary_path <- o2ipa_find_extra(run_dir, "parameter_boundary_long.tsv")
    boundary_long <- if (!is.na(boundary_path))
        o2ipa_read_tsv(boundary_path)
    else NULL
    if (!is.null(boundary_long) && "seed" %in% names(boundary_long)) {
        boundary_long$seed_id <- o2ipa_norm_seed(boundary_long$seed)
    }
    param_rows <- list()
    manifest_rows <- vector("list", nrow(manifest0))
    for (i in seq_len(nrow(manifest0))) {
        seed_id <- manifest0$seed_id[[i]]
        seed_dir <- manifest0$seed_dir[[i]]
        if (is.na(seed_dir))
            seed_dir <- ""
        fit_summary_path <- if (nzchar(seed_dir))
            file.path(seed_dir, "fit_summary.tsv")
        else NA_character_
        summary_vals <- if (!is.na(fit_summary_path) && file.exists(fit_summary_path)) {
            o2ipa_metric_map(fit_summary_path)
        }
        else {
            row <- if (!is.null(summary_tab))
                summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE]
            else NULL
            if (!is.null(row) && nrow(row))
                as.list(row[1, , drop = TRUE])
            else list()
        }
        obj <- o2ipa_choose_objective(summary_vals, objective_source = objective_source)
        params_long <- o2ipa_extract_all_params(seed_id, seed_dir, summary_tab, matrix_tab)
        param_rows[[i]] <- params_long
        missing_params <- params_long$parameter[!is.finite(params_long$value)]
        bseed <- if (!is.null(boundary_long))
            boundary_long[boundary_long$seed_id == seed_id, , drop = FALSE]
        else NULL
        target_boundary <- if (!is.null(bseed) && nrow(bseed)) {
            pname <- if ("param_prototype" %in% names(bseed))
                bseed$param_prototype
            else if ("parameter" %in% names(bseed))
                bseed$parameter
            else bseed$param_name
            bseed[pname %in% o2ipa_target_params(), , drop = FALSE]
        }
        else {
            data.frame()
        }
        boundary_risk <- FALSE
        n_near <- 0L
        if (nrow(target_boundary)) {
            status_col <- if ("bound_status" %in% names(target_boundary))
                "bound_status"
            else if ("status" %in% names(target_boundary))
                "status"
            else NA_character_
            if (!is.na(status_col)) {
                boundary_risk <- any(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior",
                  ""))
                n_near <- sum(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior", ""))
            }
        }
        fit_success <- length(missing_params) == 0L && is.finite(obj$value)
        convergence <- o2ipa_chr_from_map(summary_vals, "deoptim_stop_reason")
        if (is.na(convergence))
            convergence <- o2ipa_chr_from_map(summary_vals, "optimizer_local_convergence")
        failure_parts <- character(0)
        if (!is.finite(obj$value))
            failure_parts <- c(failure_parts, "missing_objective")
        if (length(missing_params))
            failure_parts <- c(failure_parts, paste0("missing_params:", paste(missing_params, collapse = ",")))
        manifest_rows[[i]] <- data.frame(seed_id = seed_id, seed_dir = if (nzchar(seed_dir))
            normalizePath(seed_dir, mustWork = FALSE)
        else NA_character_, fit_success = fit_success, convergence_status = convergence, objective = obj$value, objective_source = obj$source,
            objective_total = o2ipa_num_from_map(summary_vals, "objective_total"), objective_data = o2ipa_num_from_map(summary_vals,
                "objective_data"), objective_burden = o2ipa_num_from_map(summary_vals, "objective_burden"), objective_ploidy = o2ipa_num_from_map(summary_vals,
                "objective_ploidy"), burden_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw"),
            ploidy_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw"), runtime = o2ipa_num_from_map(summary_vals,
                "runtime"), parameter_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir,
                "best_params.tsv")))
                file.path(seed_dir, "best_params.tsv")
            else NA_character_, config_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir,
                "fit_config.rds")))
                file.path(seed_dir, "fit_config.rds")
            else NA_character_, visualization_available = !is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir,
                "viz_status.log")), failure_reason = if (fit_success)
                NA_character_
            else paste(failure_parts, collapse = ";"), boundary_risk = boundary_risk, number_of_parameters_near_boundary = n_near,
            stringsAsFactors = FALSE)
    }
    manifest <- do.call(rbind, manifest_rows)
    params_long <- do.call(rbind, param_rows)
    list(manifest = manifest, params_long = params_long, boundary_long = boundary_long, summary_tab = summary_tab, matrix_tab = matrix_tab)
}

o2ipa_params_wide <- function(params_long, value_col = "value") {
    seeds <- unique(params_long$seed_id)
    params <- o2ipa_target_params()
    mat <- matrix(NA_real_, nrow = length(seeds), ncol = length(params), dimnames = list(seeds, params))
    for (i in seq_len(nrow(params_long))) {
        mat[params_long$seed_id[[i]], params_long$parameter[[i]]] <- params_long[[value_col]][[i]]
    }
    as.data.frame(mat, check.names = FALSE)
}


o2ipa_call_model <- function(model_env, fn, ...) {
    if (!exists(fn, envir = model_env, inherits = TRUE))
        stop("Model helper missing: ", fn)
    get(fn, envir = model_env, inherits = TRUE)(...)
}

`%||%` <- function(x, y) {
    if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x)))
        y
    else x
}

o2pr_first_seed_cfg <- function(manifest) {
    for (p in manifest$config_file) {
        if (!is.na(p) && file.exists(p)) {
            cfg <- tryCatch(readRDS(p), error = function(e) NULL)
            if (!is.null(cfg))
                return(cfg)
        }
    }
    list(N_MIN = 22L, N_MAX = 154L, N_UNIT = 22L, DT = 0.05, start_with = "chr_number")
}

o2pr_cfg_metadata <- function(cfg) {
    data.frame(metric = c("Nmin", "Nmax", "N_UNIT", "DT", "start_with", "boundary_default", "o2_burden_feedback", "O2_growth",
        "ploidy_O2_death", "o2_S0_upper_bound", "o2_Nref", "trajectory_value_semantics"), value = c(cfg$N_MIN %||% NA, cfg$N_MAX %||%
        NA, cfg$N_UNIT %||% NA, cfg$DT %||% NA, cfg$start_with %||% NA, cfg$boundary %||% "drop", cfg$o2_burden_feedback %||%
        NA, cfg$O2_growth %||% NA, cfg$ploidy_O2_death %||% NA, cfg$o2_S0_upper_bound %||% NA, cfg$o2_Nref %||% NA, if (identical(as.character(cfg$start_with %||%
        "ploidy"), "chr_number")) {
        "weighted mean chromosome number N"
    } else {
        "weighted mean ploidy converted to N when needed"
    }), stringsAsFactors = FALSE)
}

o2pr_build_G <- function(model_env, cfg, run_params, O2) {
    fn <- get("cpp_o2simps_build_G_for_o2_triplet", envir = model_env, inherits = TRUE)
    tri <- fn(O2 = as.numeric(O2), O2_crit = as.numeric(run_params$O2_crit %||% cfg$o2_crit_init %||% 1), N0min = as.integer(cfg$N_MIN %||%
        22L), N0max = as.integer(cfg$N_MAX %||% 154L), N1min = as.integer(cfg$N_MIN %||% 22L), N1max = as.integer(cfg$N_MAX %||%
        154L), lam_max = as.numeric(run_params$lam_max), p_mis_base = as.numeric(run_params$p_mis_base %||% cfg$p_mis_base %||%
        1e-05), p_misseg = as.numeric(run_params$p_misseg %||% 0), k_o_mis = as.numeric(run_params$k_o_mis %||% 50), p_wgd = as.numeric(run_params$p_wgd %||%
        0), boundary = as.character(run_params$boundary %||% "drop"), eps_tail = 1e-08, buffer_smax = as.numeric(run_params$buffer_smax %||%
        1), buffer_beta = as.numeric(run_params$buffer_beta %||% 0), buffer_n_exp = as.numeric(run_params$buffer_n_exp %||%
        1), N_unit = as.integer(cfg$N_UNIT %||% 22L), beta_size = 0, O2_growth = isTRUE(run_params$O2_growth %||% cfg$O2_growth %||%
        TRUE), alpha_o2 = as.numeric(run_params$alpha_o2 %||% 0), gamma_growth = as.numeric(run_params$gamma_growth %||%
        1), mu_hp = as.numeric(run_params$mu_hp %||% 0), gamma_mu = as.numeric(run_params$gamma_mu %||% 1), n_O = as.numeric(run_params$n_O %||%
        1), ploidy_O2_death = as.character(run_params$ploidy_O2_death %||% cfg$ploidy_O2_death %||% "diploid_NULL"))
    G <- Matrix::sparseMatrix(i = as.integer(tri$i), j = as.integer(tri$j), x = as.numeric(tri$x), dims = c(as.integer(tri$nrow),
        as.integer(tri$ncol)), repr = "C")
    attr(G, "triplet") <- tri
    G
}

o2pr_run_params_from_vec <- function(vec, cfg) {
    rp <- as.list(vec)
    rp$o2_min <- if ("o2_min" %in% names(vec))
        vec[["o2_min"]]
    else (cfg$o2_min %||% 0)
    rp$o2_S0_upper_bound <- cfg$o2_S0_upper_bound %||% 5
    rp$o2_Nref <- cfg$o2_Nref %||% 1e+06
    rp$O2_growth <- cfg$O2_growth %||% TRUE
    rp$ploidy_O2_death <- cfg$ploidy_O2_death %||% "ploidy_related"
    rp$boundary <- cfg$boundary %||% "drop"
    rp
}

# ---- embedded C++ analytical backend -------------------------------------
.RESPONSE_CPP_CACHE <- file.path(tempdir(), "hypoxia_ltee_response_workflow_sourcecpp")
dir.create(.RESPONSE_CPP_CACHE, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  code = r"---(
#include <RcppEigen.h>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

namespace {

constexpr double kNDip = 44.0;
constexpr int kPloidyDeathUniform = 0;
constexpr int kPloidyDeathDiploidNull = 1;
constexpr int kPloidyDeathPloidyRelated = 2;
constexpr int kStartWithPloidy = 0;
constexpr int kStartWithChrNumber = 1;

// -----------------------------------------------------------------------------
// Function: trim_lower_ascii_cpp
// Purpose: Normalize mode strings for robust parsing.
// Parameters:
//   - x: Raw mode string.
// Returns:
//   std::string return value containing lowercase-trimmed ASCII text.
// -----------------------------------------------------------------------------
inline std::string trim_lower_ascii_cpp(const std::string& x) {
  size_t b = 0;
  while (b < x.size() && std::isspace(static_cast<unsigned char>(x[b]))) ++b;
  size_t e = x.size();
  while (e > b && std::isspace(static_cast<unsigned char>(x[e - 1]))) --e;
  std::string out = x.substr(b, e - b);
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return out;
}

// -----------------------------------------------------------------------------
// Function: canonical_ploidy_o2_death_mode_cpp
// Purpose: Parse canonical ploidy_O2_death mode.
// Parameters:
//   - mode_raw: Requested mode string.
// Returns:
//   int return value containing one of:
//     0=uniform, 1=diploid_NULL, 2=ploidy_related.
// -----------------------------------------------------------------------------
inline int canonical_ploidy_o2_death_mode_cpp(const std::string& mode_raw) {
  const std::string s = trim_lower_ascii_cpp(mode_raw);
  if (s.empty()) {
    stop(
      "ploidy_O2_death must be supplied as a canonical mode string. ",
      "Allowed values are: uniform, diploid_NULL, ploidy_related."
    );
  }
  if (s == "ploidy_related") {
    return kPloidyDeathPloidyRelated;
  }
  if (s == "uniform") {
    return kPloidyDeathUniform;
  }
  if (s == "diploid_null") {
    return kPloidyDeathDiploidNull;
  }
  stop(
    "Invalid ploidy_O2_death mode: '", mode_raw,
    "'. Allowed canonical values are: uniform, diploid_NULL, ploidy_related."
  );
}

// -----------------------------------------------------------------------------
// Function: ploidy_o2_death_mode_name_cpp
// Purpose: Return canonical mode name for logging/cache consistency.
// Parameters:
//   - mode_code: Integer mode code.
// Returns:
//   std::string return value containing canonical mode name.
// -----------------------------------------------------------------------------
inline std::string ploidy_o2_death_mode_name_cpp(int mode_code) {
  if (mode_code == kPloidyDeathUniform) return "uniform";
  if (mode_code == kPloidyDeathDiploidNull) return "diploid_NULL";
  return "ploidy_related";
}

// -----------------------------------------------------------------------------
// Function: canonical_start_with_mode_cpp
// Purpose: Parse canonical start_with mode.
// Parameters:
//   - mode_raw: Requested mode string.
// Returns:
//   int return value containing one of:
//     0=ploidy, 1=chr_number.
// -----------------------------------------------------------------------------
inline int canonical_start_with_mode_cpp(const std::string& mode_raw) {
  const std::string s = trim_lower_ascii_cpp(mode_raw);
  if (s.empty()) {
    stop(
      "start_with must be supplied as a canonical mode string. ",
      "Allowed values are: ploidy, chr_number."
    );
  }
  if (s == "ploidy") return kStartWithPloidy;
  if (s == "chr_number") return kStartWithChrNumber;
  stop(
    "Invalid start_with mode: '", mode_raw,
    "'. Allowed canonical values are: ploidy, chr_number."
  );
}

// -----------------------------------------------------------------------------
// Function: clamp01
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - x: Input value or vector to process.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

inline double clamp_prob_eps(double x, double eps) {
  const double eps_use = (std::isfinite(eps) && eps > 0.0 && eps < 0.5) ? eps : 1e-4;
  if (!std::isfinite(x)) return eps_use;
  if (x < eps_use) return eps_use;
  if (x > 1.0 - eps_use) return 1.0 - eps_use;
  return x;
}

inline double logit_prob_eps(double x, double eps) {
  const double p = clamp_prob_eps(x, eps);
  return std::log(p / (1.0 - p));
}

// -----------------------------------------------------------------------------
// Clamp a percent-O2 scalar to the physical [0, 100] interval. Use
// clamp_o2_pct_to_upper() for active-model O2_S0/o2_min/target semantics.
// -----------------------------------------------------------------------------
inline double clamp_o2_pct(double x) {
  if (x < 0.0) return 0.0;
  if (x > 100.0) return 100.0;
  return x;
}

inline double o2_upper_bound_use_cpp(double o2_upper) {
  if (!std::isfinite(o2_upper) || o2_upper <= 0.0) {
    stop("o2_S0_upper_bound must be finite and > 0.");
  }
  return std::min(o2_upper, 100.0);
}

// Clamp model oxygen percentages to the active model cap. R wrappers perform
// the same normalization early; C++ remains the final guard for direct calls.
inline double clamp_o2_pct_to_upper(double x, double o2_upper) {
  const double upper = o2_upper_bound_use_cpp(o2_upper);
  if (!std::isfinite(x)) return 0.0;
  if (x < 0.0) return 0.0;
  if (x > upper) return upper;
  return x;
}

// -----------------------------------------------------------------------------
// Function: hypoxia_weight_cpp
// Purpose: Compute Hill-type hypoxia weight used by growth/death modules.
// Parameters:
//   - O2_use: Oxygen level used by model rate functions.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double hypoxia_weight_cpp(double O2_use, double O2_crit_use, double n_O_use) {
  double o2_crit = (std::isfinite(O2_crit_use) && O2_crit_use >= 0.0) ? O2_crit_use : 1.0;
  o2_crit = std::max(o2_crit, 1e-12);
  const double o2 = clamp_o2_pct(O2_use);
  const double n_O = (std::isfinite(n_O_use) && n_O_use >= 0.0) ? n_O_use : 1.0;
  const double num = std::pow(o2_crit, n_O);
  const double den = num + std::pow(o2, n_O);
  if (!std::isfinite(den) || den <= 0.0) return 0.0;
  const double h = num / den;
  if (!std::isfinite(h)) return 0.0;
  return clamp01(h);
}

// -----------------------------------------------------------------------------
// Function: resource_stress_cpp
// Purpose: Compute the oxygen-only resource stress used by growth/death modules.
// Parameters:
//   - O2_use: Oxygen level used by model rate functions.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double resource_stress_cpp(
    double O2_use,
    double O2_crit_use,
    double n_O_use
) {
  return hypoxia_weight_cpp(O2_use, O2_crit_use, n_O_use);
}

// -----------------------------------------------------------------------------
// Function: constant_p_wgd_cpp
// Purpose: Compute the constant per-division WGD probability.
// Parameters:
//   - p_wgd: Per-division whole-genome doubling probability.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double constant_p_wgd_cpp(double p_wgd) {
  if (!std::isfinite(p_wgd) || p_wgd <= 0.0) return 0.0;
  return clamp01(p_wgd);
}

// -----------------------------------------------------------------------------
// Function: lambda_base_from_o2_cpp
// Purpose: Compute baseline proliferation as the maximal growth rate.
// Parameters:
//   - lam_max: Maximal proliferation rate.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_base_from_o2_cpp(double lam_max) {
  double lam_base = std::isfinite(lam_max) ? lam_max : 0.0;
  if (!std::isfinite(lam_base) || lam_base < 0.0) lam_base = 0.0;
  return lam_base;
}

// -----------------------------------------------------------------------------
// Function: lambda_base_from_resource_cpp
// Purpose: Compute baseline proliferation as the maximal growth rate.
// Parameters:
//   - lam_max: Maximal proliferation rate.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_base_from_resource_cpp(double lam_max) {
  double lam_base = std::isfinite(lam_max) ? lam_max : 0.0;
  if (!std::isfinite(lam_base) || lam_base < 0.0) lam_base = 0.0;
  return lam_base;
}

// -----------------------------------------------------------------------------
// Function: lambda_eff_soft_o2_only_cpp
// Purpose: Compute effective growth rate with oxygen-stress growth damping.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - lam_max: Maximal proliferation rate.
//   - alpha_o2: Oxygen-mediated growth-penalty strength.
//   - gamma_growth: Exponent for oxygen-mediated ploidy growth penalty.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_eff_soft_o2_only_cpp(
    int N_state,
    double O2_use,
    double lam_max,
    double alpha_o2,
    double gamma_growth,
    double O2_crit_use,
    double n_O_use
) {
  const double lam_base = lambda_base_from_o2_cpp(lam_max);
  if (lam_base <= 0.0) return 0.0;
  const double alpha_use = (std::isfinite(alpha_o2) && alpha_o2 > 0.0) ? alpha_o2 : 0.0;
  const double gamma_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double N_ratio = std::max(static_cast<double>(N_state) / kNDip, 0.0);
  const double h_o2 = hypoxia_weight_cpp(O2_use, O2_crit_use, n_O_use);
  const double denom = 1.0 + alpha_use * h_o2 * std::pow(N_ratio, gamma_use);
  if (!std::isfinite(denom) || denom <= 0.0) return 0.0;
  const double lam_eff = lam_base / denom;
  if (!std::isfinite(lam_eff) || lam_eff < 0.0) return 0.0;
  return lam_eff;
}

// -----------------------------------------------------------------------------
// Function: lambda_eff_soft_cpp
// Purpose: Compute effective growth rate with soft resource/ploidy penalty.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - lam_max: Maximal proliferation rate.
//   - alpha_o2: Resource-mediated growth-penalty strength.
//   - gamma_growth: Exponent for resource-mediated ploidy growth penalty.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_eff_soft_cpp(
    int N_state,
    double O2_use,
    double lam_max,
    double alpha_o2,
    double gamma_growth,
    double O2_crit_use,
    double n_O_use
) {
  const double lam_base = lambda_base_from_resource_cpp(lam_max);
  if (lam_base <= 0.0) return 0.0;
  const double alpha_use = (std::isfinite(alpha_o2) && alpha_o2 > 0.0) ? alpha_o2 : 0.0;
  const double gamma_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double N_ratio = std::max(static_cast<double>(N_state) / kNDip, 0.0);
  const double h_resource = resource_stress_cpp(
    O2_use,
    O2_crit_use,
    n_O_use
  );
  const double denom = 1.0 + alpha_use * h_resource * std::pow(N_ratio, gamma_use);
  if (!std::isfinite(denom) || denom <= 0.0) return 0.0;
  const double lam_eff = lam_base / denom;
  if (!std::isfinite(lam_eff) || lam_eff < 0.0) return 0.0;
  return lam_eff;
}

// -----------------------------------------------------------------------------
// Function: lambda_eff_runtime_cpp
// Purpose: Dispatch proliferation rate according to the O2_growth runtime switch.
// -----------------------------------------------------------------------------
inline double lambda_eff_runtime_cpp(
    int N_state,
    double O2_use,
    double lam_max,
    bool o2_growth,
    double alpha_o2,
    double gamma_growth,
    double O2_crit_use,
    double n_O_use
) {
  if (!o2_growth) {
    return lambda_base_from_o2_cpp(lam_max);
  }
  return lambda_eff_soft_o2_only_cpp(
    N_state,
    O2_use,
    lam_max,
    alpha_o2,
    gamma_growth,
    O2_crit_use,
    n_O_use
  );
}

// -----------------------------------------------------------------------------
// Function: mu_eff_soft_cpp
// Purpose: Compute effective hypoxia-linked death rate with optional ploidy modulation.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - mu_hp: Hypoxia-linked high-ploidy death strength.
//   - gamma_mu: Exponent for high-ploidy hypoxia death above diploid reference.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
//   - ploidy_O2_death_mode: Mode code parsed from ploidy_O2_death.
//     Allowed values:
//       uniform       -> mu_eff = mu_hp * h(O2)
//       diploid_NULL  -> mu_eff = mu_hp * h(O2) * (1 + max(N/N_dip - 1, 0)^gamma_mu)
//       ploidy_related-> mu_eff = mu_hp * h(O2) * (N/N_dip)^gamma_mu
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double mu_eff_soft_cpp(
    int N_state,
    double O2_use,
    double mu_hp,
    double gamma_mu,
    double O2_crit_use,
    double n_O_use,
    int ploidy_O2_death_mode
) {
  const double mu_hp_use = (std::isfinite(mu_hp) && mu_hp > 0.0) ? mu_hp : 0.0;
  if (mu_hp_use <= 0.0) return 0.0;
  const double h_o2 = hypoxia_weight_cpp(O2_use, O2_crit_use, n_O_use);
  if (h_o2 <= 0.0) return 0.0;
  if (ploidy_O2_death_mode == kPloidyDeathUniform) {
    const double mu_eff = mu_hp_use * h_o2;
    if (!std::isfinite(mu_eff) || mu_eff < 0.0) return 0.0;
    return mu_eff;
  }
  const double gamma_mu_use = (std::isfinite(gamma_mu) && gamma_mu > 0.0) ? gamma_mu : 1.0;
  const double n_ratio = std::max(static_cast<double>(N_state) / kNDip, 0.0);
  if (ploidy_O2_death_mode == kPloidyDeathDiploidNull) {
    const double above_dip = std::max(n_ratio - 1.0, 0.0);
    const double mu_eff = mu_hp_use * h_o2 * (1.0 + std::pow(above_dip, gamma_mu_use));
    if (!std::isfinite(mu_eff) || mu_eff < 0.0) return 0.0;
    return mu_eff;
  }
  const double mu_eff = mu_hp_use * h_o2 * std::pow(n_ratio, gamma_mu_use);
  if (!std::isfinite(mu_eff) || mu_eff < 0.0) return 0.0;
  return mu_eff;
}

// -----------------------------------------------------------------------------
// Function: sigmoid01
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - z: Function-specific input argument.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double sigmoid01(double z) {
  return 1.0 / (1.0 + std::exp(-z));
}

// -----------------------------------------------------------------------------
// Function: quantize_o2_key
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - o2_pct: Oxygen level in percent scale (0-100).
//   - bin_pct: Function-specific input argument.
// Returns:
//   int return value containing the computed result.
// -----------------------------------------------------------------------------
inline int quantize_o2_key(double o2_pct, double bin_pct) {
  const double o2_use = clamp_o2_pct(o2_pct);
  const double bin_use = (std::isfinite(bin_pct) && bin_pct > 0.0) ? bin_pct : 1e-3;
  const double raw = o2_use / bin_use;
  const double cap = static_cast<double>(std::numeric_limits<int>::max() / 4);
  const double clamped = std::min(std::max(raw, 0.0), cap);
  return static_cast<int>(std::llround(clamped));
}

// -----------------------------------------------------------------------------
// Function: boundary_mode
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
// Returns:
//   int return value containing the computed result.
// -----------------------------------------------------------------------------
inline int boundary_mode(const std::string& boundary) {
  if (boundary == "drop") return 0;
  if (boundary == "absorb_minmax") return 1;
  stop("boundary must be one of: drop, absorb_minmax");
}

// -----------------------------------------------------------------------------
// Function: append_with_boundary
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - Np: Function-specific input argument.
//   - row_min: Minimum valid row state index for sparse insertion.
//   - row_max: Maximum valid row state index for sparse insertion.
//   - col_1based: 1-based sparse column index for sparse insertion.
//   - value: Numeric transition value to append.
//   - bmode: Encoded boundary mode used internally by C++ helpers.
//   - ii: Sparse-triplet row index accumulator.
//   - jj: Sparse-triplet column index accumulator.
//   - xx: Sparse-triplet value accumulator.
//   - dropped_value: Optional accumulator for out-of-grid dropped transition mass.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
inline void append_with_boundary(
    int Np,
    int row_min,
    int row_max,
    int col_1based,
    double value,
    int bmode,
    std::vector<int>& ii,
    std::vector<int>& jj,
    std::vector<double>& xx,
    double* dropped_value = nullptr
) {
  if (Np < row_min || Np > row_max) {
    if (bmode == 1) {
      const int Np2 = std::max(std::min(Np, row_max), row_min);
      ii.push_back(Np2 - row_min + 1);
      jj.push_back(col_1based);
      xx.push_back(value);
    } else if (dropped_value != nullptr) {
      // Under boundary=drop, out-of-grid offspring are not written to live
      // states; caller is responsible for routing this mass into dead buffer.
      *dropped_value += value;
    }
    return;
  }
  ii.push_back(Np - row_min + 1);
  jj.push_back(col_1based);
  xx.push_back(value);
}

// -----------------------------------------------------------------------------
// Function: binom_prob_int
// Purpose: Numerically robust binomial PMF evaluator for integer n in [0, N].
// Parameters:
//   - n: Number of successes.
//   - N: Number of Bernoulli trials.
//   - p: Per-trial success probability.
// Returns:
//   double return value containing PMF value in [0,1].
// -----------------------------------------------------------------------------
inline double binom_prob_int(int n, int N, double p) {
  if (N < 0) return 0.0;
  if (n < 0 || n > N) return 0.0;
  const double p_use = clamp01(p);
  const double v = R::dbinom(n, N, p_use, false);
  if (!std::isfinite(v) || v < 0.0) return 0.0;
  return v;
}

// -----------------------------------------------------------------------------
// Function: buffering_survival_modifier
// Purpose: Compute buffering-model survival for any missegregated
//   daughter with m affected chromosome copies.
// Parameters:
//   - q: Mother chromosome count state.
//   - m_misseg: Number of missegregated chromosome copies.
//   - buffer_smax: Maximum per-copy survival factor, constrained to [0,1].
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - n_chr: Number of modeled chromosome classes.
// Returns:
//   double return value containing survival modifier in [0,1].
// -----------------------------------------------------------------------------
inline double buffering_survival_modifier(
    int q,
    int m_misseg,
    double buffer_smax,
    double buffer_beta,
    double buffer_n_exp,
    int n_chr
) {
  if (m_misseg <= 0) return 1.0;
  if (q <= 0) return 1.0;
  const double smax_use = clamp01(std::isfinite(buffer_smax) ? buffer_smax : 1.0);
  const double beta_use = (std::isfinite(buffer_beta) && buffer_beta >= 0.0) ? buffer_beta : 0.0;
  const double n_exp_use = (std::isfinite(buffer_n_exp) && buffer_n_exp >= 0.0) ? buffer_n_exp : 1.0;
  const double n_chr_use = (n_chr > 0) ? static_cast<double>(n_chr) : 22.0;
  const double ratio = (2.0 * n_chr_use) / static_cast<double>(q);
  double sN = smax_use * std::exp(-beta_use * std::pow(std::max(ratio, 0.0), n_exp_use));
  sN = clamp01(std::isfinite(sN) ? sN : 0.0);
  if (sN <= 0.0) return 0.0;
  const double log_survival = static_cast<double>(m_misseg) * std::log(sN);
  if (!std::isfinite(log_survival)) return 0.0;
  const double survival = std::exp(log_survival);
  if (!std::isfinite(survival)) return 0.0;
  return clamp01(survival);
}

// -----------------------------------------------------------------------------
// Function: finalize_pr_delta_internal
// Purpose: Finalize coarse ploidy transition weights and dropped daughter mass.
// -----------------------------------------------------------------------------
inline void finalize_pr_delta_internal(
    int shift_offset,
    const std::vector<double>& shift_mass,
    double survivors_total,
    std::vector<int>& ts_out,
    std::vector<double>& prob_out,
    double& mass_dropped
) {
  if (!std::isfinite(survivors_total) || survivors_total < 0.0) survivors_total = 0.0;

  const double dead_daughters = std::max(0.0, 2.0 - survivors_total);
  mass_dropped = std::max(0.0, std::min(1.0, dead_daughters / 2.0));

  ts_out.clear();
  prob_out.clear();
  ts_out.reserve(shift_mass.size());
  prob_out.reserve(shift_mass.size());
  for (int t = -shift_offset; t <= shift_offset; ++t) {
    const double w = shift_mass[static_cast<size_t>(t + shift_offset)];
    if (!std::isfinite(w) || w <= 0.0) continue;
    ts_out.push_back(t);
    prob_out.push_back(w);
  }

  if (ts_out.empty()) {
    ts_out.assign(1, 0);
    prob_out.assign(1, 0.0);
  }
}

// -----------------------------------------------------------------------------
// Function: o2simps_pr_delta_internal
// Purpose: Build symmetric buffering transition weights, matching the legacy O2
//   buffering model: gain and loss daughters share the same buffering survival.
// -----------------------------------------------------------------------------
void o2simps_pr_delta_internal(
    int N,
    double p,
    double eps_tail,
    double buffer_smax,
    double buffer_beta,
    double buffer_n_exp,
    int N_unit,
    std::vector<int>& ts_out,
    std::vector<double>& prob_out,
    double& mass_dropped
) {
  const int N_use = std::max(0, N);
  const double p_use = clamp01(p);
  const double eps_use = (std::isfinite(eps_tail) && eps_tail > 0.0) ? eps_tail : 0.0;
  const int shift_offset = N_use;
  std::vector<double> shift_mass(static_cast<size_t>(2 * shift_offset + 1), 0.0);

  double survivors_total = 0.0;
  for (int n = 0; n <= N_use; ++n) {
    const double pn = binom_prob_int(n, N_use, p_use);
    if (pn <= 0.0) continue;
    if (eps_use > 0.0 && pn < eps_use) continue;
    const int delta_gain = n;
    const int delta_loss = -n;
    const double s_buf = buffering_survival_modifier(
      N_use,
      n,
      buffer_smax,
      buffer_beta,
      buffer_n_exp,
      N_unit
    );
    const double w_gain = pn * s_buf;
    const double w_loss = pn * s_buf;
    if (w_gain > 0.0) shift_mass[static_cast<size_t>(shift_offset + delta_gain)] += w_gain;
    if (w_loss > 0.0) shift_mass[static_cast<size_t>(shift_offset + delta_loss)] += w_loss;
    survivors_total += (w_gain + w_loss);
  }

  finalize_pr_delta_internal(shift_offset, shift_mass, survivors_total, ts_out, prob_out, mass_dropped);
}

} // namespace

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_pr_delta_vec
// Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
// Parameters:
//   - N: Ploidy state value or chromosome-copy count.
//   - p: Missegregation probability parameter.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_pr_delta_vec(
    int N,
    double p,
    double eps_tail = 1e-8,
    double buffer_smax = 1.0,
    double buffer_beta = 0.0,
    double buffer_n_exp = 1.0,
    int N_unit = 22
) {
  std::vector<int> ts;
  std::vector<double> prob;
  double mass_dropped = 0.0;

  o2simps_pr_delta_internal(
    N,
    p,
    eps_tail,
    buffer_smax,
    buffer_beta,
    buffer_n_exp,
    N_unit,
    ts,
    prob,
    mass_dropped
  );

  return List::create(
    _["ts"] = IntegerVector(ts.begin(), ts.end()),
    _["prob"] = NumericVector(prob.begin(), prob.end()),
    _["mass_dropped"] = mass_dropped
  );
}

// -----------------------------------------------------------------------------
// Inputs:
// - Ntot: nonnegative live/effective oxygen-demand proxy.
// - o2_S0, o2_min: natural-scale percent O2 values.
// - kappa_O, o2_Nref: log-demand slope and demand normalization scale.
// - o2_S0_upper_bound: active model cap; C++ enforces this for direct callers.
//
// Returns:
// - target percent O2 values clamped to [0, min(o2_S0_upper_bound, 100)].
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cpp_o2simps_o2_window_supply(
    NumericVector Ntot,
    double o2_S0 = 0.5,
    double kappa_O = 1.0,
    double o2_Nref = 1e6,
    double o2_min = 0.0,
    double o2_S0_upper_bound = 5.0
) {
  const int n = Ntot.size();
  NumericVector out(n);
  const double o2_upper_use = o2_upper_bound_use_cpp(o2_S0_upper_bound);
  const double o2_S0_use = clamp_o2_pct_to_upper((std::isfinite(o2_S0) && o2_S0 >= 0.0) ? o2_S0 : 0.5, o2_upper_use);
  const double kappa_use = (std::isfinite(kappa_O) && kappa_O > 0.0) ? kappa_O : 1.0;
  const double Nref_use = (std::isfinite(o2_Nref) && o2_Nref > 0.0) ? o2_Nref : 1e6;
  const double O2_min_use = clamp_o2_pct_to_upper((std::isfinite(o2_min) && o2_min >= 0.0) ? o2_min : 0.0, o2_upper_use);

  for (int i = 0; i < n; ++i) {
    const double Nlive = (std::isfinite(Ntot[i]) && Ntot[i] > 0.0) ? Ntot[i] : 0.0;
    const double burden_ratio = Nlive / Nref_use;
    double o2_target = o2_S0_use - kappa_use * std::log1p(burden_ratio);
    if (!std::isfinite(o2_target)) o2_target = O2_min_use;
    o2_target = std::max(O2_min_use, o2_target);
    out[i] = clamp_o2_pct_to_upper(o2_target, o2_upper_use);
  }

  return out;
}

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_build_B_total_triplet
// Purpose: Build total missegregation transition operator on ploidy grid.
// Parameters:
//   - Nmin: Minimum ploidy state on source grid.
//   - Nmax: Maximum ploidy state on source grid.
//   - p_vec: State-specific missegregation probability vector.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_build_B_total_triplet(
    int Nmin,
    int Nmax,
    NumericVector p_vec,
    std::string boundary = "drop",
    double eps_tail = 1e-8,
    double buffer_smax = 1.0,
    double buffer_beta = 0.0,
    double buffer_n_exp = 1.0,
    int N_unit = 22
) {
  const int R = Nmax - Nmin + 1;
  if (R <= 0) stop("Nmax must be >= Nmin");

  const int p_len = p_vec.size();
  if (!(p_len == 1 || p_len == R)) stop("p_vec length must be 1 or R");

  const int bmode = boundary_mode(boundary);

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  ii.reserve(static_cast<size_t>(R) * 12);
  jj.reserve(static_cast<size_t>(R) * 12);
  xx.reserve(static_cast<size_t>(R) * 12);

  for (int col = 0; col < R; ++col) {
    const int N = Nmin + col;
    double pN = (p_len == 1) ? p_vec[0] : p_vec[col];
    pN = std::max(0.0, std::min(1.0, pN));

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    o2simps_pr_delta_internal(
      N,
      pN,
      eps_tail,
      buffer_smax,
      buffer_beta,
      buffer_n_exp,
      N_unit,
      ts,
      pr,
      mass_dropped
    );

    const int col_1based = col + 1;
    const int K = static_cast<int>(ts.size());
    for (int k = 0; k < K; ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;
      // Signed-shift contract: each (t, w) already encodes final daughter
      // displacement and mass, so write exactly once to N + t.
      append_with_boundary(
        N + t,
        Nmin,
        Nmax,
        col_1based,
        w,
        bmode,
        ii,
        jj,
        xx
      );
    }
  }

  return List::create(
    _["i"] = IntegerVector(ii.begin(), ii.end()),
    _["j"] = IntegerVector(jj.begin(), jj.end()),
    _["x"] = NumericVector(xx.begin(), xx.end()),
    _["nrow"] = R,
    _["ncol"] = R
  );
}

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_build_B_WGD_triplet
// Purpose: Build WGD transition operator between source and doubled-ploidy grids.
// Parameters:
//   - N0min: Minimum ploidy state on source grid.
//   - N0max: Maximum ploidy state on source grid.
//   - N1min: Minimum ploidy state on doubled-state target grid.
//   - N1max: Maximum ploidy state on doubled-state target grid.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - wgd_value: Function-specific input argument.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_build_B_WGD_triplet(
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    std::string boundary = "drop",
    double wgd_value = 1.0
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;
  if (R0 <= 0 || R1 <= 0) stop("Nmax must be >= Nmin for source and target grids");

  const int bmode = boundary_mode(boundary);
  const double val = wgd_value;

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  ii.reserve(static_cast<size_t>(R0));
  jj.reserve(static_cast<size_t>(R0));
  xx.reserve(static_cast<size_t>(R0));

  for (int col = 0; col < R0; ++col) {
    const int N0 = N0min + col;
    const int Np = 2 * N0;
    const int col_1based = col + 1;
    append_with_boundary(
      Np,
      N1min,
      N1max,
      col_1based,
      val,
      bmode,
      ii,
      jj,
      xx
    );
  }

  return List::create(
    _["i"] = IntegerVector(ii.begin(), ii.end()),
    _["j"] = IntegerVector(jj.begin(), jj.end()),
    _["x"] = NumericVector(xx.begin(), xx.end()),
    _["nrow"] = R1,
    _["ncol"] = R0
  );
}

namespace {

// -----------------------------------------------------------------------------
// Function: append_block_with_boundary
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - Np: Function-specific input argument.
//   - row_min: Minimum valid row state index for sparse insertion.
//   - row_max: Maximum valid row state index for sparse insertion.
//   - row_offset_1based: Function-specific input argument.
//   - col_1based: 1-based sparse column index for sparse insertion.
//   - value: Numeric transition value to append.
//   - bmode: Encoded boundary mode used internally by C++ helpers.
//   - ii: Sparse-triplet row index accumulator.
//   - jj: Sparse-triplet column index accumulator.
//   - xx: Sparse-triplet value accumulator.
//   - dropped_value: Optional accumulator for out-of-grid dropped transition mass.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
inline void append_block_with_boundary(
    int Np,
    int row_min,
    int row_max,
    int row_offset_1based,
    int col_1based,
    double value,
    int bmode,
    std::vector<int>& ii,
    std::vector<int>& jj,
    std::vector<double>& xx,
    double* dropped_value = nullptr
) {
  if (value == 0.0) return;

  int Np_use = Np;
  if (Np_use < row_min || Np_use > row_max) {
    if (bmode == 0) {
      if (dropped_value != nullptr) {
        // Under boundary=drop, out-of-grid offspring mass is accumulated here
        // so caller can add it to dead_buffer_rate.
        *dropped_value += value;
      }
      return;
    }
    Np_use = std::max(std::min(Np_use, row_max), row_min);
  }

  const int row_local_1based = Np_use - row_min + 1;
  ii.push_back(row_offset_1based + row_local_1based - 1);
  jj.push_back(col_1based);
  xx.push_back(value);
}

// -----------------------------------------------------------------------------
// Function: resolve_pmis_for_death
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - mu_eff: Effective death rate for the current ploidy state.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale for death-linked missegregation.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double resolve_pmis_for_death(
    double mu_eff,
    double p_mis_base,
    double p_misseg,
    double k_o_mis
) {
  const double p_base = clamp01(std::isfinite(p_mis_base) ? p_mis_base : 1e-5);
  const double p_amp = (std::isfinite(p_misseg) && p_misseg > 0.0) ? p_misseg : 0.0;
  const double k_use = (std::isfinite(k_o_mis) && k_o_mis > 0.0) ? k_o_mis : 1e-12;
  const double mu_use = std::max(mu_eff, 0.0);
  const double frac = mu_use / (mu_use + k_use);
  const double delta_p = p_amp * frac;
  return clamp01(p_base + delta_p);
}

} // namespace

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_build_G_for_o2_triplet
// Purpose: Build division-related live-state generator at the current oxygen/burden condition.
// Parameters:
//   - O2: Oxygen level used by model rate functions.
//   - N0min: Minimum ploidy state on the single chromosome-count grid.
//   - N0max: Maximum ploidy state on the single chromosome-count grid.
//   - N1min: Legacy argument kept for interface stability (unused).
//   - N1max: Legacy argument kept for interface stability (unused).
//   - lam_max: Maximal proliferation rate.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale for death-linked missegregation (mu_eff scale).
//   - p_wgd: Constant per-division WGD probability.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_build_G_for_o2_triplet(
    double O2,
    double O2_crit,
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double lam_max,
    double p_mis_base,
    double p_misseg,
    double k_o_mis,
    double p_wgd = 0.0,
    std::string boundary = "drop",
    double eps_tail = 1e-8,
    double buffer_smax = 1.0,
    double buffer_beta = 0.0,
    double buffer_n_exp = 1.0,
    int N_unit = 22,
    double beta_size = 0.0,
    bool O2_growth = true,
    double alpha_o2 = 0.0,
    double gamma_growth = 1.0,
    double mu_hp = 0.0,
    double gamma_mu = 1.0,
    double n_O = 1.0,
    std::string ploidy_O2_death = "diploid_NULL"
) {
  const int R = N0max - N0min + 1;
  if (R <= 0) stop("Nmax must be >= Nmin");
  (void)N1min;
  (void)N1max;

  const int bmode = boundary_mode(boundary);

  const double O2_use = clamp_o2_pct(O2);
  const double o2_crit_use = (std::isfinite(O2_crit) && O2_crit >= 0.0) ? O2_crit : 1.0;
  const bool o2_growth_use = O2_growth;
  const double alpha_o2_use = (std::isfinite(alpha_o2) && alpha_o2 > 0.0) ? alpha_o2 : 0.0;
  const double gamma_growth_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double mu_hp_use = (std::isfinite(mu_hp) && mu_hp > 0.0) ? mu_hp : 0.0;
  const double gamma_mu_use = (std::isfinite(gamma_mu) && gamma_mu > 0.0) ? gamma_mu : 1.0;
  const double n_O_use = (std::isfinite(n_O) && n_O >= 0.0) ? n_O : 1.0;
  const int ploidy_O2_death_mode_use = canonical_ploidy_o2_death_mode_cpp(ploidy_O2_death);
  (void)beta_size;
  auto lam_for_N = [&](int N_state) -> double {
    return lambda_eff_runtime_cpp(
      N_state,
      O2_use,
      lam_max,
      o2_growth_use,
      alpha_o2_use,
      gamma_growth_use,
      o2_crit_use,
      n_O_use
    );
  };
  auto mu_for_N = [&](int N_state) -> double {
    return mu_eff_soft_cpp(
      N_state,
      O2_use,
      mu_hp_use,
      gamma_mu_use,
      o2_crit_use,
      n_O_use,
      ploidy_O2_death_mode_use
    );
  };

  const double pw = constant_p_wgd_cpp(p_wgd);

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  std::vector<double> dead_buffer_rate(static_cast<size_t>(R), 0.0);
  std::vector<double> misseg_nonviable_rate(static_cast<size_t>(R), 0.0);
  std::vector<double> boundary_dropped_rate_vec(static_cast<size_t>(R), 0.0);
  std::vector<double> misseg_nonviable_division_prob(static_cast<size_t>(R), 0.0);
  std::vector<double> misseg_nonviable_daughters_per_division(static_cast<size_t>(R), 0.0);
  ii.reserve(static_cast<size_t>(R) * 20);
  jj.reserve(static_cast<size_t>(R) * 20);
  xx.reserve(static_cast<size_t>(R) * 20);

  for (int c = 0; c < R; ++c) {
    const int N = N0min + c;
    const double lam_N = lam_for_N(N);
    const double mu_N = mu_for_N(N);
    const double p_mis_N = resolve_pmis_for_death(
      mu_N,
      p_mis_base,
      p_misseg,
      k_o_mis
    );
    const int col_1based = c + 1;

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    o2simps_pr_delta_internal(
      N,
      p_mis_N,
      eps_tail,
      buffer_smax,
      buffer_beta,
      buffer_n_exp,
      N_unit,
      ts,
      pr,
      mass_dropped
    );

    const double scale_mis = lam_N * (1.0 - pw);
    const double scale_wgd = lam_N * pw;
    double boundary_dropped_rate = 0.0;
    const double dead_daughters_per_division = std::max(0.0, std::min(2.0, 2.0 * mass_dropped));
    {
      // Expected nonviable offspring inflow from missegregation-linked loss.
      double nonviable_daughters_per_division = (1.0 - pw) * dead_daughters_per_division;
      if (!std::isfinite(nonviable_daughters_per_division) || nonviable_daughters_per_division < 0.0) {
        nonviable_daughters_per_division = 0.0;
      }
      if (nonviable_daughters_per_division > 2.0) nonviable_daughters_per_division = 2.0;
      double nonviable_rate = lam_N * nonviable_daughters_per_division;
      if (!std::isfinite(nonviable_rate) || nonviable_rate < 0.0) nonviable_rate = 0.0;
      misseg_nonviable_division_prob[static_cast<size_t>(col_1based - 1)] = std::min(nonviable_daughters_per_division, 1.0);
      misseg_nonviable_daughters_per_division[static_cast<size_t>(col_1based - 1)] = nonviable_daughters_per_division;
      misseg_nonviable_rate[static_cast<size_t>(col_1based - 1)] = nonviable_rate;
      dead_buffer_rate[static_cast<size_t>(col_1based - 1)] = nonviable_rate;
    }
    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;
      // Signed-shift contract: each (t, w) already encodes final daughter
      // displacement and mass, so write exactly once to N + t.
      append_block_with_boundary(
        N + t,
        N0min,
        N0max,
        1,
        col_1based,
        scale_mis * w,
        bmode,
        ii,
        jj,
        xx,
        &boundary_dropped_rate
      );
    }

    ii.push_back(col_1based);
    jj.push_back(col_1based);
    xx.push_back(-lam_N);

    append_block_with_boundary(
      2 * N,
      N0min,
      N0max,
      1,
      col_1based,
      scale_wgd,
      bmode,
      ii,
      jj,
      xx,
      &boundary_dropped_rate
    );
    if (!std::isfinite(boundary_dropped_rate) || boundary_dropped_rate < 0.0) {
      boundary_dropped_rate = 0.0;
    }
    boundary_dropped_rate_vec[static_cast<size_t>(col_1based - 1)] = boundary_dropped_rate;
    // Boundary-drop losses are also routed to the dead buffer so mass does not
    // disappear from the represented system.
    dead_buffer_rate[static_cast<size_t>(col_1based - 1)] += boundary_dropped_rate;
  }

  return List::create(
    _["i"] = IntegerVector(ii.begin(), ii.end()),
    _["j"] = IntegerVector(jj.begin(), jj.end()),
    _["x"] = NumericVector(xx.begin(), xx.end()),
    _["nrow"] = R,
    _["ncol"] = R,
    _["dead_buffer_rate"] = NumericVector(dead_buffer_rate.begin(), dead_buffer_rate.end()),
    _["misseg_nonviable_rate"] = NumericVector(misseg_nonviable_rate.begin(), misseg_nonviable_rate.end()),
    _["boundary_dropped_rate"] = NumericVector(boundary_dropped_rate_vec.begin(), boundary_dropped_rate_vec.end()),
    _["misseg_nonviable_division_prob"] = NumericVector(
      misseg_nonviable_division_prob.begin(),
      misseg_nonviable_division_prob.end()
    ),
    _["misseg_nonviable_daughters_per_division"] = NumericVector(
      misseg_nonviable_daughters_per_division.begin(),
      misseg_nonviable_daughters_per_division.end()
    )
  );
}

namespace {

using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

struct SparseCacheEntry {
  SpMat mat;
  std::vector<double> dead_buffer_rate;
};

template <typename T>
// -----------------------------------------------------------------------------
// Function: hash_combine_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - seed: Random-seed value used for reproducible initialization.
//   - value: Numeric transition value to append.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
inline void hash_combine_cpp(std::size_t& seed, const T& value) {
  seed ^= std::hash<T>{}(value) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

// -----------------------------------------------------------------------------
// Function: bits_of_double_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - x: Input value or vector to process.
// Returns:
//   std::uint64_t return value containing the computed result.
// -----------------------------------------------------------------------------
inline std::uint64_t bits_of_double_cpp(double x) {
  std::uint64_t out = 0ULL;
  std::memcpy(&out, &x, sizeof(double));
  return out;
}

// -----------------------------------------------------------------------------
// Function: g_cache_signature_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - N0min: Minimum ploidy state on source grid.
//   - N0max: Maximum ploidy state on source grid.
//   - N1min: Legacy interface argument (unused in single-layer dynamics).
//   - N1max: Legacy interface argument (unused in single-layer dynamics).
//   - lam_max: Maximal proliferation rate.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale for death-linked missegregation (mu_eff scale).
//   - p_wgd: Constant per-division WGD probability.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
// Returns:
//   std::size_t return value containing the computed result.
// -----------------------------------------------------------------------------
inline std::size_t g_cache_signature_cpp(
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double lam_max,
    double p_mis_base,
    double p_misseg,
    double k_o_mis,
    double p_wgd,
    const std::string& boundary,
    double eps_tail,
    double buffer_smax,
    double buffer_beta,
    double buffer_n_exp,
    double beta_size,
    bool o2_growth,
    double alpha_o2,
    double O2_crit,
    double gamma_growth,
    double mu_hp,
    double gamma_mu,
    double n_O,
    int ploidy_O2_death_mode,
    int N_unit
) {
  std::size_t seed = 0ULL;
  hash_combine_cpp(seed, N0min);
  hash_combine_cpp(seed, N0max);
  hash_combine_cpp(seed, N1min);
  hash_combine_cpp(seed, N1max);
  hash_combine_cpp(seed, bits_of_double_cpp(lam_max));
  hash_combine_cpp(seed, bits_of_double_cpp(p_mis_base));
  hash_combine_cpp(seed, bits_of_double_cpp(p_misseg));
  hash_combine_cpp(seed, bits_of_double_cpp(k_o_mis));
  hash_combine_cpp(seed, bits_of_double_cpp(p_wgd));
  hash_combine_cpp(seed, boundary);
  hash_combine_cpp(seed, bits_of_double_cpp(eps_tail));
  hash_combine_cpp(seed, bits_of_double_cpp(buffer_smax));
  hash_combine_cpp(seed, bits_of_double_cpp(buffer_beta));
  hash_combine_cpp(seed, bits_of_double_cpp(buffer_n_exp));
  hash_combine_cpp(seed, bits_of_double_cpp(beta_size));
  hash_combine_cpp(seed, o2_growth ? 1 : 0);
  hash_combine_cpp(seed, bits_of_double_cpp(alpha_o2));
  hash_combine_cpp(seed, bits_of_double_cpp(O2_crit));
  hash_combine_cpp(seed, bits_of_double_cpp(gamma_growth));
  hash_combine_cpp(seed, bits_of_double_cpp(mu_hp));
  hash_combine_cpp(seed, bits_of_double_cpp(gamma_mu));
  hash_combine_cpp(seed, bits_of_double_cpp(n_O));
  hash_combine_cpp(seed, ploidy_O2_death_mode);
  hash_combine_cpp(seed, N_unit);
  return seed;
}

// -----------------------------------------------------------------------------
// Function: vector_sum_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - x: Input value or vector to process.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double vector_sum_cpp(const std::vector<double>& x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

// -----------------------------------------------------------------------------
// Function: sparse_mv_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - G: Generator or transition matrix used for state propagation.
//   - x: Input value or vector to process.
//   - y: Function-specific input argument.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
inline void sparse_mv_cpp(
    const SparseCacheEntry& G,
    const std::vector<double>& x,
    std::vector<double>& y
) {
  if (static_cast<int>(x.size()) != G.mat.cols() || static_cast<int>(y.size()) != G.mat.rows()) {
    stop("Sparse matvec dimension mismatch.");
  }
  Eigen::Map<const Eigen::VectorXd> xmap(x.data(), static_cast<int>(x.size()));
  Eigen::Map<Eigen::VectorXd> ymap(y.data(), static_cast<int>(y.size()));
  ymap.noalias() = G.mat * xmap;
}

// -----------------------------------------------------------------------------
// Function: o2_window_supply_scalar_cpp
// Purpose: Compute oxygen target from viable burden using a logarithmic
//   supply-demand form with lower oxygen floor.
// Parameters:
//   - Ntot: Total predicted cell count (or burden proxy) at current time.
//   - o2_S0: Baseline oxygen supply level at near-zero burden (%).
//   - kappa_O: Function-specific input argument.
//   - o2_Nref: Fixed viable-cell scaling constant for demand normalization.
//   - o2_min: Lower floor for oxygen target in logarithmic supply-demand model (%).
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double o2_window_supply_scalar_cpp(
    double Ntot,
    double o2_S0,
    double kappa_O,
    double o2_Nref,
    double o2_min,
    double o2_S0_upper_bound
) {
  const double o2_upper_use = o2_upper_bound_use_cpp(o2_S0_upper_bound);
  const double o2_S0_use = clamp_o2_pct_to_upper((std::isfinite(o2_S0) && o2_S0 >= 0.0) ? o2_S0 : 0.5, o2_upper_use);
  const double kappa_use = (std::isfinite(kappa_O) && kappa_O > 0.0) ? kappa_O : 1.0;
  const double Nref_use = (std::isfinite(o2_Nref) && o2_Nref > 0.0) ? o2_Nref : 1e6;
  const double O2_min_use = clamp_o2_pct_to_upper((std::isfinite(o2_min) && o2_min >= 0.0) ? o2_min : 0.0, o2_upper_use);
  const double Nlive = (std::isfinite(Ntot) && Ntot > 0.0) ? Ntot : 0.0;
  const double burden_ratio = Nlive / Nref_use;
  double o2_target = o2_S0_use - kappa_use * std::log1p(burden_ratio);
  if (!std::isfinite(o2_target)) o2_target = O2_min_use;
  o2_target = std::max(O2_min_use, o2_target);
  return clamp_o2_pct_to_upper(o2_target, o2_upper_use);
}

// -----------------------------------------------------------------------------
// Function: death_rate_for_N_cpp
// Purpose: Compute live->dead transfer death rate with optional ploidy modulation.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - O2_crit_use: Hill critical oxygen scale.
//   - mu_hp_use: Hypoxia-linked high-ploidy death strength.
//   - gamma_mu_use: Exponent for high-ploidy hypoxia death above diploid reference.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
//   - ploidy_O2_death_mode: Parsed mode code for hypoxia-death ploidy dependence.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double death_rate_for_N_cpp(
    int N_state,
    double O2_use,
    double O2_crit_use,
    double mu_hp_use,
    double gamma_mu_use,
    double n_O_use,
    int ploidy_O2_death_mode
) {
  return mu_eff_soft_cpp(
    N_state,
    O2_use,
    mu_hp_use,
    gamma_mu_use,
    O2_crit_use,
    n_O_use,
    ploidy_O2_death_mode
  );
}

// -----------------------------------------------------------------------------
// Function: build_sparse_cache_entry_from_triplet
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - tri: Function-specific input argument.
// Returns:
//   SparseCacheEntry return value containing the computed result.
// -----------------------------------------------------------------------------
inline SparseCacheEntry build_sparse_cache_entry_from_triplet(const List& tri) {
  IntegerVector ii = tri["i"];
  IntegerVector jj = tri["j"];
  NumericVector xx = tri["x"];
  const int nrow = as<int>(tri["nrow"]);
  const int ncol = as<int>(tri["ncol"]);
  const int n = xx.size();
  if (ii.size() != n || jj.size() != n) {
    stop("Triplet i/j/x length mismatch.");
  }
  SparseCacheEntry out;
  out.mat.resize(nrow, ncol);
  std::vector<Eigen::Triplet<double, int>> triplets;
  triplets.reserve(static_cast<size_t>(n));
  for (int e = 0; e < n; ++e) {
    const int r = ii[e] - 1;
    const int c = jj[e] - 1;
    if (r < 0 || r >= nrow || c < 0 || c >= ncol) {
      stop("Triplet index out of bounds.");
    }
    triplets.emplace_back(r, c, xx[e]);
  }
  out.mat.setFromTriplets(triplets.begin(), triplets.end());
  out.mat.makeCompressed();
  out.dead_buffer_rate.assign(static_cast<size_t>(ncol), 0.0);
  if (tri.containsElementNamed("dead_buffer_rate")) {
    NumericVector db = tri["dead_buffer_rate"];
    if (db.size() != ncol) {
      stop("dead_buffer_rate length mismatch.");
    }
    for (int i = 0; i < ncol; ++i) {
      double v = db[i];
      if (!std::isfinite(v) || v < 0.0) v = 0.0;
      out.dead_buffer_rate[static_cast<size_t>(i)] = v;
    }
  }
  return out;
}

} // namespace

// -----------------------------------------------------------------------------
// Inputs:
// - sim_args: natural-scale single-scenario simulation contract assembled by R.
//   init/live/dead states are on the N0min:N0max chromosome grid; observation
//   steps are integer DT steps.
// - oxygen fields are percent O2; o2_S0_upper_bound is enforced here as the
//   final guard for o2_S0, o2_min, O2 target, and O2 state.
// - event parameters encode the tested generator semantics: lambda is a parent
//   division rate, p_wgd diverts a division to one 2N transition, mu_hp feeds
//   only the hypoxia-death compartment, and nonviable daughters feed dead buffer.
//
// Returns:
// - observation/terminal live, dead-hypoxia, dead-buffer, total burden summaries;
//   optional per-observation state matrices and O2 trajectories.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_simulate_one(List sim_args) {
  NumericVector init_state = as<NumericVector>(sim_args["init_state"]);
  const bool has_init_dead_hypoxia_state = sim_args.containsElementNamed("init_dead_hypoxia_state");
  const bool has_init_dead_buffer_state = sim_args.containsElementNamed("init_dead_buffer_state");
  NumericVector init_dead_hypoxia_state;
  NumericVector init_dead_buffer_state;
  if (has_init_dead_hypoxia_state) {
    init_dead_hypoxia_state = as<NumericVector>(sim_args["init_dead_hypoxia_state"]);
  }
  if (has_init_dead_buffer_state) {
    init_dead_buffer_state = as<NumericVector>(sim_args["init_dead_buffer_state"]);
  }
  int N0min = as<int>(sim_args["N0min"]);
  int N0max = as<int>(sim_args["N0max"]);
  int N1min = as<int>(sim_args["N1min"]);
  int N1max = as<int>(sim_args["N1max"]);
  IntegerVector obs_steps = as<IntegerVector>(sim_args["obs_steps"]);
  int sim_end_step = as<int>(sim_args["sim_end_step"]);
  double DT = as<double>(sim_args["DT"]);
  double dose = as<double>(sim_args["dose"]);
  double dose_ref = as<double>(sim_args["dose_ref"]);
  double treat_day = as<double>(sim_args["treat_day"]);
  bool fit_treatment = as<bool>(sim_args["fit_treatment"]);
  double alpha = as<double>(sim_args["alpha"]);
  double gamma = as<double>(sim_args["gamma"]);
  double tx_mult_min = as<double>(sim_args["tx_mult_min"]);
  bool crowding_enabled = as<bool>(sim_args["crowding_enabled"]);
  std::string crowding = as<std::string>(sim_args["crowding"]);
  double K = as<double>(sim_args["K"]);
  double min_pop = as<double>(sim_args["min_pop"]);
  double O2_crit = as<double>(sim_args["O2_crit"]);
  bool o2_feedback = as<bool>(sim_args["o2_feedback"]);
  double o2_S0 = as<double>(sim_args["o2_S0"]);
  double kappa_O = as<double>(sim_args["kappa_O"]);
  double tau_O2 = as<double>(sim_args["tau_O2"]);
  double o2_Nref = as<double>(sim_args["o2_Nref"]);
  double o2_min = as<double>(sim_args["o2_min"]);
  double o2_S0_upper_bound = sim_args.containsElementNamed("o2_S0_upper_bound")
    ? as<double>(sim_args["o2_S0_upper_bound"])
    : 5.0;
  double eta_o2 = as<double>(sim_args["eta_o2"]);
  double o2_cache_bin_pct = as<double>(sim_args["o2_cache_bin_pct"]);
  double o2_cache_hysteresis_pct = as<double>(sim_args["o2_cache_hysteresis_pct"]);
  bool o2_cache_profile = as<bool>(sim_args["o2_cache_profile"]);
  double lam_max = as<double>(sim_args["lam_max"]);
  double p_mis_base = as<double>(sim_args["p_mis_base"]);
  double p_misseg = as<double>(sim_args["p_misseg"]);
  double k_o_mis = as<double>(sim_args["k_o_mis"]);
  double p_wgd = as<double>(sim_args["p_wgd"]);
  std::string boundary = as<std::string>(sim_args["boundary"]);
  double eps_tail = as<double>(sim_args["eps_tail"]);
  double buffer_smax = as<double>(sim_args["buffer_smax"]);
  double buffer_beta = as<double>(sim_args["buffer_beta"]);
  double buffer_n_exp = as<double>(sim_args["buffer_n_exp"]);
  int N_unit = as<int>(sim_args["N_unit"]);
  double beta_size = as<double>(sim_args["beta_size"]);
  bool O2_growth = as<bool>(sim_args["O2_growth"]);
  double alpha_o2 = as<double>(sim_args["alpha_o2"]);
  double gamma_growth = as<double>(sim_args["gamma_growth"]);
  double mu_hp = as<double>(sim_args["mu_hp"]);
  double gamma_mu = as<double>(sim_args["gamma_mu"]);
  double n_O = as<double>(sim_args["n_O"]);
  std::string ploidy_O2_death = as<std::string>(sim_args["ploidy_O2_death"]);
  std::string start_with = as<std::string>(sim_args["start_with"]);
  double k_clear = as<double>(sim_args["k_clear"]);
  NumericVector vol_by_N = as<NumericVector>(sim_args["vol_by_N"]);
  double burden_floor = as<double>(sim_args["burden_floor"]);
  bool return_full_trajectory = as<bool>(sim_args["return_full_trajectory"]);

  const int R = N0max - N0min + 1;
  if (R <= 0) stop("Nmax must be >= Nmin.");
  (void)N1min;
  (void)N1max;
  const int D = R;
  if (init_state.size() != D) {
    stop("init_state length mismatch: expected N0max - N0min + 1 live-state entries.");
  }
  if (has_init_dead_hypoxia_state && init_dead_hypoxia_state.size() != D) {
    stop("init_dead_hypoxia_state length mismatch: expected N0max - N0min + 1 entries.");
  }
  if (has_init_dead_buffer_state && init_dead_buffer_state.size() != D) {
    stop("init_dead_buffer_state length mismatch: expected N0max - N0min + 1 entries.");
  }
  if (vol_by_N.size() != R) stop("vol_by_N length mismatch.");

  if (!(crowding == "logistic" || crowding == "gompertz")) {
    stop("crowding must be logistic or gompertz.");
  }

  const double DT_use = (std::isfinite(DT) && DT > 0.0) ? DT : 0.5;
  const double K_use = (std::isfinite(K) && K > 0.0) ? K : 1e12;
  const double min_pop_use = (std::isfinite(min_pop) && min_pop > 0.0) ? min_pop : 1e-12;
  const double dose_ref_use = (std::isfinite(dose_ref) && dose_ref > 0.0) ? dose_ref : 30.0;
  const double dose_use = (std::isfinite(dose) && dose > 0.0) ? dose : 0.0;
  const double tx_min_use = clamp01(tx_mult_min);
  const double burden_floor_use = (std::isfinite(burden_floor) && burden_floor >= 0.0) ? burden_floor : 0.0;

  std::vector<int> obs_v(obs_steps.size());
  for (int i = 0; i < obs_steps.size(); ++i) obs_v[static_cast<size_t>(i)] = obs_steps[i];
  std::vector<int> step_unique = obs_v;
  std::sort(step_unique.begin(), step_unique.end());
  step_unique.erase(std::unique(step_unique.begin(), step_unique.end()), step_unique.end());
  const int max_obs_step = step_unique.empty() ? 0 : step_unique.back();
  const int final_step = std::max(sim_end_step, max_obs_step);

  std::unordered_map<int, int> step_to_idx;
  step_to_idx.reserve(step_unique.size() * 2 + 1);
  for (size_t i = 0; i < step_unique.size(); ++i) {
    step_to_idx[step_unique[i]] = static_cast<int>(i);
  }

  std::vector<double> Ntot_live_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_dead_hypoxia_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_dead_buffer_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_dead_total_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_total_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_live_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_dead_hypoxia_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_dead_buffer_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_dead_total_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_total_at_step(step_unique.size(), NA_REAL);
  std::vector<double> O2_target_at_step(step_unique.size(), NA_REAL);
  std::vector<double> O2_eff_at_step(step_unique.size(), NA_REAL);
  NumericMatrix live_state_at_step(
    return_full_trajectory ? static_cast<int>(step_unique.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_hypoxia_state_at_step(
    return_full_trajectory ? static_cast<int>(step_unique.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_buffer_state_at_step(
    return_full_trajectory ? static_cast<int>(step_unique.size()) : 0,
    return_full_trajectory ? R : 0
  );

  std::vector<double> v_live(static_cast<size_t>(D), 0.0);
  std::vector<double> v_dead_hypoxia(static_cast<size_t>(D), 0.0);
  std::vector<double> v_dead_buffer(static_cast<size_t>(D), 0.0);
  std::copy(init_state.begin(), init_state.end(), v_live.begin());
  if (has_init_dead_hypoxia_state) {
    std::copy(init_dead_hypoxia_state.begin(), init_dead_hypoxia_state.end(), v_dead_hypoxia.begin());
  }
  if (has_init_dead_buffer_state) {
    std::copy(init_dead_buffer_state.begin(), init_dead_buffer_state.end(), v_dead_buffer.begin());
  }
  std::vector<double> growth(static_cast<size_t>(D), 0.0);
  std::vector<double> death_flow_hypoxia(static_cast<size_t>(D), 0.0);
  std::vector<double> death_flow_buffer(static_cast<size_t>(D), 0.0);

  // Shared across scenario calls in the same worker process.
  // We keep one active parameter signature at a time so cache is reused
  // within one objective (same params), then reset when params change.
  static std::size_t active_sig = std::numeric_limits<std::size_t>::max();
  static std::unordered_map<std::size_t, SparseCacheEntry> shared_G_cache;

  const int ploidy_O2_death_mode_use = canonical_ploidy_o2_death_mode_cpp(ploidy_O2_death);
  const int start_with_mode_use = canonical_start_with_mode_cpp(start_with);
  const double n_O_use = (std::isfinite(n_O) && n_O >= 0.0) ? n_O : 1.0;
  const std::size_t cur_sig = g_cache_signature_cpp(
    N0min,
    N0max,
    N1min,
    N1max,
    lam_max,
    p_mis_base,
    p_misseg,
    k_o_mis,
    p_wgd,
    boundary,
    eps_tail,
    buffer_smax,
    buffer_beta,
    buffer_n_exp,
    beta_size,
    O2_growth,
    alpha_o2,
    O2_crit,
    gamma_growth,
    mu_hp,
    gamma_mu,
    n_O_use,
    ploidy_O2_death_mode_use,
    N_unit
  );
  if (cur_sig != active_sig) {
    shared_G_cache.clear();
    shared_G_cache.reserve(256);
    active_sig = cur_sig;
  }

  const double O2_crit_use = (std::isfinite(O2_crit) && O2_crit >= 0.0) ? O2_crit : 1.0;
  const double o2_upper_use = o2_upper_bound_use_cpp(o2_S0_upper_bound);
  const double o2_S0_use = clamp_o2_pct_to_upper((std::isfinite(o2_S0) ? o2_S0 : 0.5), o2_upper_use);
  const double kappa_O_use = (std::isfinite(kappa_O) ? kappa_O : 1.0);
  const double tau_use = (std::isfinite(tau_O2) && tau_O2 > 0.0) ? tau_O2 : 2.0;
  const double alpha_tau = 1.0 - std::exp(-DT_use / tau_use);
  const double o2_Nref_use = (std::isfinite(o2_Nref) && o2_Nref > 0.0) ? o2_Nref : 1e6;
  const double o2_min_use = clamp_o2_pct_to_upper((std::isfinite(o2_min) && o2_min >= 0.0) ? o2_min : 0.0, o2_upper_use);
  const double eta_o2_use = (std::isfinite(eta_o2) && eta_o2 >= 0.0) ? eta_o2 : 1.0;
  const double N_unit_use = (N_unit > 0) ? static_cast<double>(N_unit) : 22.0;
  const double o2_bin_use = (std::isfinite(o2_cache_bin_pct) && o2_cache_bin_pct > 0.0) ? o2_cache_bin_pct : 1e-3;
  const double o2_hyst_use = (std::isfinite(o2_cache_hysteresis_pct) && o2_cache_hysteresis_pct >= 0.0) ? o2_cache_hysteresis_pct : 0.0;
  const double k_clear_use = (std::isfinite(k_clear) && k_clear >= 0.0) ? k_clear : 0.0;
  const std::string ploidy_O2_death_mode_name = ploidy_o2_death_mode_name_cpp(ploidy_O2_death_mode_use);
  (void) o2_cache_profile;
  int cache_g_build = 0;
  int cache_g_hit = 0;
  int cache_g_hysteresis = 0;
  bool has_last_key = false;
  std::size_t last_key = 0ULL;
  double last_o2_eff = 0.0;
  std::vector<double> o2_demand_weight(static_cast<size_t>(D), 1.0);
  for (int i = 0; i < D; ++i) {
    const double N_state = static_cast<double>(N0min + i);
    const double ratio = (start_with_mode_use == kStartWithChrNumber)
      ? std::max(0.0, N_state / kNDip)
      : std::max(0.0, (N_state / N_unit_use) / 2.0);
    double w = std::pow(ratio, eta_o2_use);
    if (!std::isfinite(w) || w < 0.0) w = 1.0;
    o2_demand_weight[static_cast<size_t>(i)] = w;
  }
  auto compute_o2_demand_eff = [&](const std::vector<double>& live_state) -> double {
    double s = 0.0;
    for (int i = 0; i < D; ++i) {
      s += live_state[static_cast<size_t>(i)] * o2_demand_weight[static_cast<size_t>(i)];
    }
    if (!std::isfinite(s) || s < 0.0) s = 0.0;
    return s;
  };
  double O2_state = clamp_o2_pct_to_upper(o2_S0_use, o2_upper_use);
  if (o2_feedback) {
    O2_state = o2_window_supply_scalar_cpp(
      compute_o2_demand_eff(v_live),
      o2_S0_use,
      kappa_O_use,
      o2_Nref_use,
      o2_min_use,
      o2_upper_use
    );
    O2_state = clamp_o2_pct_to_upper(O2_state, o2_upper_use);
  }

  for (int step = 0; step <= final_step; ++step) {
    const double Ntot_live_now = vector_sum_cpp(v_live);
    const double Ntot_live_eff_for_o2_now = compute_o2_demand_eff(v_live);
    double O2_target_now = clamp_o2_pct_to_upper(o2_S0_use, o2_upper_use);
    if (o2_feedback) {
      O2_target_now = o2_window_supply_scalar_cpp(
        Ntot_live_eff_for_o2_now,
        o2_S0_use,
        kappa_O_use,
        o2_Nref_use,
        o2_min_use,
        o2_upper_use
      );
    }
    O2_target_now = clamp_o2_pct_to_upper(O2_target_now, o2_upper_use);
    const double O2_eff_now = clamp_o2_pct_to_upper(O2_state, o2_upper_use);

    auto it_obs = step_to_idx.find(step);
    if (it_obs != step_to_idx.end()) {
      const int idx = it_obs->second;
      const double Ntot_dead_h_now = vector_sum_cpp(v_dead_hypoxia);
      const double Ntot_dead_b_now = vector_sum_cpp(v_dead_buffer);
      const double Ntot_dead_now = Ntot_dead_h_now + Ntot_dead_b_now;
      const double Ntot_total_now = Ntot_live_now + Ntot_dead_now;
      Ntot_live_at_step[static_cast<size_t>(idx)] = Ntot_live_now;
      Ntot_dead_hypoxia_at_step[static_cast<size_t>(idx)] = Ntot_dead_h_now;
      Ntot_dead_buffer_at_step[static_cast<size_t>(idx)] = Ntot_dead_b_now;
      Ntot_dead_total_at_step[static_cast<size_t>(idx)] = Ntot_dead_now;
      Ntot_total_at_step[static_cast<size_t>(idx)] = Ntot_total_now;
      double burden_live_now = 0.0;
      double burden_dead_h_now = 0.0;
      double burden_dead_b_now = 0.0;
      double burden_dead_now = 0.0;
      double burden_total_now = 0.0;
      for (int i = 0; i < R; ++i) {
        const size_t idx_i = static_cast<size_t>(i);
        const double n_live_i = v_live[idx_i];
        const double n_dead_h_i = v_dead_hypoxia[idx_i];
        const double n_dead_b_i = v_dead_buffer[idx_i];
        const double n_dead_i = n_dead_h_i + n_dead_b_i;
        const double n_total_i = n_live_i + n_dead_i;
        burden_live_now += n_live_i * vol_by_N[i];
        burden_dead_h_now += n_dead_h_i * vol_by_N[i];
        burden_dead_b_now += n_dead_b_i * vol_by_N[i];
        burden_dead_now += n_dead_i * vol_by_N[i];
        burden_total_now += n_total_i * vol_by_N[i];
      }
      Vmm3_live_at_step[static_cast<size_t>(idx)] = burden_live_now;
      Vmm3_dead_hypoxia_at_step[static_cast<size_t>(idx)] = burden_dead_h_now;
      Vmm3_dead_buffer_at_step[static_cast<size_t>(idx)] = burden_dead_b_now;
      Vmm3_dead_total_at_step[static_cast<size_t>(idx)] = burden_dead_now;
      Vmm3_total_at_step[static_cast<size_t>(idx)] = burden_total_now;
      O2_target_at_step[static_cast<size_t>(idx)] = O2_target_now;
      O2_eff_at_step[static_cast<size_t>(idx)] = O2_eff_now;
      if (return_full_trajectory) {
        for (int i = 0; i < R; ++i) {
          live_state_at_step(idx, i) = v_live[static_cast<size_t>(i)];
          dead_hypoxia_state_at_step(idx, i) = v_dead_hypoxia[static_cast<size_t>(i)];
          dead_buffer_state_at_step(idx, i) = v_dead_buffer[static_cast<size_t>(i)];
        }
      }
    }
    if (step >= final_step) break;

    const double t = static_cast<double>(step) * DT_use;
    double tx_mult = 1.0;
    if (fit_treatment) {
      if (!(t < treat_day) && dose_use > 0.0) {
        double dose_scaled = dose_use / dose_ref_use;
        if (!std::isfinite(dose_scaled) || dose_scaled < 0.0) dose_scaled = 0.0;
        tx_mult = std::exp(-alpha * std::pow(dose_scaled, gamma));
      } else {
        tx_mult = 1.0;
      }
      if (!std::isfinite(tx_mult)) tx_mult = tx_min_use;
      if (tx_mult < tx_min_use) tx_mult = tx_min_use;
      if (tx_mult > 1.0) tx_mult = 1.0;
    }

    O2_state = O2_state + alpha_tau * (O2_target_now - O2_state);
    double O2_eff = clamp_o2_pct_to_upper(O2_state, o2_upper_use);

    const int o2_key = quantize_o2_key(O2_eff, o2_bin_use);
    std::size_t gkey = 0ULL;
    hash_combine_cpp(gkey, o2_key);
    if (o2_hyst_use > 0.0 && has_last_key &&
        std::abs(O2_eff - last_o2_eff) <= o2_hyst_use) {
      gkey = last_key;
      ++cache_g_hysteresis;
    }
    auto itG = shared_G_cache.find(gkey);
    if (itG == shared_G_cache.end()) {
      const List tri = cpp_o2simps_build_G_for_o2_triplet(
        O2_eff,
        O2_crit,
        N0min,
        N0max,
        N1min,
        N1max,
        lam_max,
        p_mis_base,
        p_misseg,
        k_o_mis,
        p_wgd,
        boundary,
        eps_tail,
        buffer_smax,
        buffer_beta,
        buffer_n_exp,
        N_unit,
        beta_size,
        O2_growth,
        alpha_o2,
        gamma_growth,
        mu_hp,
        gamma_mu,
        n_O_use,
        ploidy_O2_death_mode_name
      );
      SparseCacheEntry entry = build_sparse_cache_entry_from_triplet(tri);
      auto insert_res = shared_G_cache.emplace(gkey, std::move(entry));
      itG = insert_res.first;
      ++cache_g_build;
    } else {
      ++cache_g_hit;
    }
    has_last_key = true;
    last_key = gkey;
    last_o2_eff = O2_eff;

    sparse_mv_cpp(itG->second, v_live, growth);
    // Crowding scales division-linked activity when the config switch is enabled.
    double crowd_mult = 1.0;
    if (crowding_enabled) {
      if (crowding == "logistic") {
        crowd_mult = std::max(0.0, 1.0 - (Ntot_live_now / K_use));
      } else {
        crowd_mult = std::exp(-(Ntot_live_now / K_use));
      }
    }
    if (!std::isfinite(crowd_mult) || crowd_mult < 0.0) crowd_mult = 0.0;
    const double scalar = DT_use * crowd_mult * tx_mult;
    for (int i = 0; i < D; ++i) {
      const int N_state = N0min + i;
      const double mu_i = death_rate_for_N_cpp(
        N_state,
        O2_eff,
        O2_crit_use,
        mu_hp,
        gamma_mu,
        n_O_use,
        ploidy_O2_death_mode_use
      );
      const double src_live = v_live[static_cast<size_t>(i)];
      // Hypoxia death flow is independent of crowding/treatment scaling.
      double flow_h_i = DT_use * mu_i * src_live;
      if (!std::isfinite(flow_h_i) || flow_h_i < 0.0) flow_h_i = 0.0;
      death_flow_hypoxia[static_cast<size_t>(i)] = flow_h_i;
      double db_rate_i = 0.0;
      if (static_cast<size_t>(i) < itG->second.dead_buffer_rate.size()) {
        db_rate_i = itG->second.dead_buffer_rate[static_cast<size_t>(i)];
      }
      if (!std::isfinite(db_rate_i) || db_rate_i < 0.0) db_rate_i = 0.0;
      // Dead-buffer inflow from mitotic nonviability + boundary-dropped offspring
      // (not a continuous death hazard).
      double flow_b_i = scalar * db_rate_i * src_live;
      if (!std::isfinite(flow_b_i) || flow_b_i < 0.0) flow_b_i = 0.0;
      death_flow_buffer[static_cast<size_t>(i)] = flow_b_i;
    }
    for (size_t i = 0; i < v_live.size(); ++i) {
      const double next_v = v_live[i] + scalar * growth[i] - death_flow_hypoxia[i];
      if (!std::isfinite(next_v) || next_v < 0.0) {
        v_live[i] = 0.0;
      } else {
        v_live[i] = next_v;
      }
    }
    for (size_t i = 0; i < v_dead_hypoxia.size(); ++i) {
      const double dead_h_prev = v_dead_hypoxia[i];
      const double dead_h_next = dead_h_prev + death_flow_hypoxia[i] - DT_use * k_clear_use * dead_h_prev;
      if (!std::isfinite(dead_h_next) || dead_h_next < 0.0) {
        v_dead_hypoxia[i] = 0.0;
      } else {
        v_dead_hypoxia[i] = dead_h_next;
      }
      const double dead_b_prev = v_dead_buffer[i];
      const double dead_b_next = dead_b_prev + death_flow_buffer[i] - DT_use * k_clear_use * dead_b_prev;
      if (!std::isfinite(dead_b_next) || dead_b_next < 0.0) {
        v_dead_buffer[i] = 0.0;
      } else {
        v_dead_buffer[i] = dead_b_next;
      }
    }
    if (!return_full_trajectory &&
        vector_sum_cpp(v_live) <= min_pop_use &&
        (vector_sum_cpp(v_dead_hypoxia) + vector_sum_cpp(v_dead_buffer)) <= min_pop_use) break;
  }

  NumericVector Ntot_live_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_dead_hypoxia_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_dead_buffer_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_dead_total_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_total_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_live_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_dead_hypoxia_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_dead_buffer_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_dead_total_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_total_obs(obs_v.size(), NA_REAL);
  NumericVector O2_target_obs(obs_v.size(), NA_REAL);
  NumericVector O2_eff_obs(obs_v.size(), NA_REAL);
  NumericMatrix live_state_obs(
    return_full_trajectory ? static_cast<int>(obs_v.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_hypoxia_state_obs(
    return_full_trajectory ? static_cast<int>(obs_v.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_buffer_state_obs(
    return_full_trajectory ? static_cast<int>(obs_v.size()) : 0,
    return_full_trajectory ? R : 0
  );
  for (int i = 0; i < static_cast<int>(obs_v.size()); ++i) {
    auto it = step_to_idx.find(obs_v[static_cast<size_t>(i)]);
    if (it == step_to_idx.end()) {
      Ntot_live_obs[i] = min_pop_use;
      Ntot_dead_hypoxia_obs[i] = 0.0;
      Ntot_dead_buffer_obs[i] = 0.0;
      Ntot_dead_total_obs[i] = 0.0;
      Ntot_total_obs[i] = min_pop_use;
      Vmm3_live_obs[i] = burden_floor_use;
      Vmm3_dead_hypoxia_obs[i] = 0.0;
      Vmm3_dead_buffer_obs[i] = 0.0;
      Vmm3_dead_total_obs[i] = 0.0;
      Vmm3_total_obs[i] = burden_floor_use;
      O2_target_obs[i] = o2_window_supply_scalar_cpp(
        0.0,
        o2_S0_use,
        kappa_O_use,
        o2_Nref_use,
        o2_min_use,
        o2_upper_use
      );
      O2_eff_obs[i] = O2_target_obs[i];
      if (return_full_trajectory) {
        for (int j = 0; j < R; ++j) {
          live_state_obs(i, j) = 0.0;
          dead_hypoxia_state_obs(i, j) = 0.0;
          dead_buffer_state_obs(i, j) = 0.0;
        }
      }
      continue;
    }
    const int idx = it->second;
    double nv_live = Ntot_live_at_step[static_cast<size_t>(idx)];
    double nv_dead_h = Ntot_dead_hypoxia_at_step[static_cast<size_t>(idx)];
    double nv_dead_b = Ntot_dead_buffer_at_step[static_cast<size_t>(idx)];
    double nv_dead = Ntot_dead_total_at_step[static_cast<size_t>(idx)];
    double nv_total = Ntot_total_at_step[static_cast<size_t>(idx)];
    double bv_live = Vmm3_live_at_step[static_cast<size_t>(idx)];
    double bv_dead_h = Vmm3_dead_hypoxia_at_step[static_cast<size_t>(idx)];
    double bv_dead_b = Vmm3_dead_buffer_at_step[static_cast<size_t>(idx)];
    double bv_dead = Vmm3_dead_total_at_step[static_cast<size_t>(idx)];
    double bv_total = Vmm3_total_at_step[static_cast<size_t>(idx)];
    double o2_target_val = O2_target_at_step[static_cast<size_t>(idx)];
    double o2_eff_val = O2_eff_at_step[static_cast<size_t>(idx)];
    if (!std::isfinite(nv_live)) nv_live = min_pop_use;
    if (!std::isfinite(nv_dead_h) || nv_dead_h < 0.0) nv_dead_h = 0.0;
    if (!std::isfinite(nv_dead_b) || nv_dead_b < 0.0) nv_dead_b = 0.0;
    if (!std::isfinite(nv_dead) || nv_dead < 0.0) nv_dead = (nv_dead_h + nv_dead_b);
    if (!std::isfinite(nv_total)) nv_total = min_pop_use;
    if (!std::isfinite(bv_live)) bv_live = burden_floor_use;
    if (!std::isfinite(bv_dead_h) || bv_dead_h < 0.0) bv_dead_h = 0.0;
    if (!std::isfinite(bv_dead_b) || bv_dead_b < 0.0) bv_dead_b = 0.0;
    if (!std::isfinite(bv_dead) || bv_dead < 0.0) bv_dead = (bv_dead_h + bv_dead_b);
    if (!std::isfinite(bv_total)) bv_total = burden_floor_use;
    if (!std::isfinite(o2_target_val)) {
      o2_target_val = o2_window_supply_scalar_cpp(
        nv_live,
        o2_S0_use,
        kappa_O_use,
        o2_Nref_use,
        o2_min_use,
        o2_upper_use
      );
    }
    if (!std::isfinite(o2_eff_val)) o2_eff_val = o2_target_val;
    Ntot_live_obs[i] = nv_live;
    Ntot_dead_hypoxia_obs[i] = nv_dead_h;
    Ntot_dead_buffer_obs[i] = nv_dead_b;
    Ntot_dead_total_obs[i] = nv_dead;
    Ntot_total_obs[i] = nv_total;
    Vmm3_live_obs[i] = bv_live;
    Vmm3_dead_hypoxia_obs[i] = bv_dead_h;
    Vmm3_dead_buffer_obs[i] = bv_dead_b;
    Vmm3_dead_total_obs[i] = bv_dead;
    Vmm3_total_obs[i] = bv_total;
    O2_target_obs[i] = o2_target_val;
    O2_eff_obs[i] = o2_eff_val;
    if (return_full_trajectory) {
      for (int j = 0; j < R; ++j) {
        live_state_obs(i, j) = live_state_at_step(idx, j);
        dead_hypoxia_state_obs(i, j) = dead_hypoxia_state_at_step(idx, j);
        dead_buffer_state_obs(i, j) = dead_buffer_state_at_step(idx, j);
      }
    }
  }

  NumericVector frac_N_live(R, 0.0);
  double total_frac = 0.0;
  for (int i = 0; i < R; ++i) {
    const double val = v_live[static_cast<size_t>(i)];
    frac_N_live[i] = val;
    total_frac += val;
  }
  if (total_frac > 0.0 && std::isfinite(total_frac)) {
    for (int i = 0; i < R; ++i) frac_N_live[i] = frac_N_live[i] / total_frac;
  } else {
    const double u = 1.0 / static_cast<double>(R);
    for (int i = 0; i < R; ++i) frac_N_live[i] = u;
  }

  const double Ntot_live_terminal = vector_sum_cpp(v_live);
  const double Ntot_dead_hypoxia_terminal = vector_sum_cpp(v_dead_hypoxia);
  const double Ntot_dead_buffer_terminal = vector_sum_cpp(v_dead_buffer);
  const double Ntot_dead_total_terminal = Ntot_dead_hypoxia_terminal + Ntot_dead_buffer_terminal;
  const double Ntot_total_terminal = Ntot_live_terminal + Ntot_dead_total_terminal;
  double Vmm3_live_terminal = 0.0;
  double Vmm3_dead_hypoxia_terminal = 0.0;
  double Vmm3_dead_buffer_terminal = 0.0;
  double Vmm3_dead_total_terminal = 0.0;
  double Vmm3_total_terminal = 0.0;
  for (int i = 0; i < R; ++i) {
    const size_t idx_i = static_cast<size_t>(i);
    const double n_live_i = v_live[idx_i];
    const double n_dead_h_i = v_dead_hypoxia[idx_i];
    const double n_dead_b_i = v_dead_buffer[idx_i];
    const double n_dead_i = n_dead_h_i + n_dead_b_i;
    const double n_total_i = n_live_i + n_dead_i;
    Vmm3_live_terminal += n_live_i * vol_by_N[i];
    Vmm3_dead_hypoxia_terminal += n_dead_h_i * vol_by_N[i];
    Vmm3_dead_buffer_terminal += n_dead_b_i * vol_by_N[i];
    Vmm3_dead_total_terminal += n_dead_i * vol_by_N[i];
    Vmm3_total_terminal += n_total_i * vol_by_N[i];
  }

  return List::create(
    // Backward-compatible aliases:
    _["Ntot_obs"] = Ntot_total_obs,
    _["Vmm3_obs"] = Vmm3_total_obs,
    _["frac_N"] = frac_N_live,
    // Explicit live/dead-consistent outputs:
    _["Ntot_live_obs"] = Ntot_live_obs,
    _["Ntot_dead_hypoxia_obs"] = Ntot_dead_hypoxia_obs,
    _["Ntot_dead_buffer_obs"] = Ntot_dead_buffer_obs,
    _["Ntot_dead_total_obs"] = Ntot_dead_total_obs,
    _["Ntot_total_obs"] = Ntot_total_obs,
    _["Vmm3_live_obs"] = Vmm3_live_obs,
    _["Vmm3_dead_hypoxia_obs"] = Vmm3_dead_hypoxia_obs,
    _["Vmm3_dead_buffer_obs"] = Vmm3_dead_buffer_obs,
    _["Vmm3_dead_total_obs"] = Vmm3_dead_total_obs,
    _["Vmm3_total_obs"] = Vmm3_total_obs,
    _["Ntot_live_terminal"] = Ntot_live_terminal,
    _["Ntot_dead_hypoxia_terminal"] = Ntot_dead_hypoxia_terminal,
    _["Ntot_dead_buffer_terminal"] = Ntot_dead_buffer_terminal,
    _["Ntot_dead_total_terminal"] = Ntot_dead_total_terminal,
    _["Ntot_total_terminal"] = Ntot_total_terminal,
    _["Vmm3_live_terminal"] = Vmm3_live_terminal,
    _["Vmm3_dead_hypoxia_terminal"] = Vmm3_dead_hypoxia_terminal,
    _["Vmm3_dead_buffer_terminal"] = Vmm3_dead_buffer_terminal,
    _["Vmm3_dead_total_terminal"] = Vmm3_dead_total_terminal,
    _["Vmm3_total_terminal"] = Vmm3_total_terminal,
    _["O2_target_obs"] = O2_target_obs,
    _["O2_eff_obs"] = O2_eff_obs,
    _["frac_N_live"] = frac_N_live,
    _["live_state_obs"] = live_state_obs,
    _["dead_hypoxia_state_obs"] = dead_hypoxia_state_obs,
    _["dead_buffer_state_obs"] = dead_buffer_state_obs,
    _["cache_g_build"] = cache_g_build,
    _["cache_g_hit"] = cache_g_hit,
    _["cache_g_hysteresis"] = cache_g_hysteresis,
    _["cache_bin_pct"] = o2_bin_use,
    _["cache_hysteresis_pct"] = o2_hyst_use
  );
}

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_objective_components_map
// Purpose: Compute MAP objective components using log-normal burden likelihood
//   and continuous single-cell ploidy mixture likelihood with balanced
//   2N/4N tumor-group aggregation for ploidy loss.
// Parameters:
//   - ploidy_z_list: Per-tumor continuous single-cell ploidy observations.
//   - mu_by_N: Representative ploidy value for each modeled N state.
//   - sigma_burden: Log-normal observation SD for burden.
//   - sigma_ploidy: Gaussian observation SD for single-cell ploidy.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale on mu_eff for missegregation amplification.
//   - o2_min: Lower floor for oxygen target in logarithmic supply-demand model (%).
// Returns:
//   List return value containing per-modality mean NLL components.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_objective_components_map(
    List scenario_data,
    List objective_data,
    List state_data,
    List sim_args
) {
  IntegerVector cohort_code = as<IntegerVector>(scenario_data["cohort_code"]);
  NumericVector dose_vec = as<NumericVector>(scenario_data["dose_vec"]);
  NumericVector treat_day_vec = as<NumericVector>(scenario_data["treat_day_vec"]);
  List obs_steps_list = as<List>(scenario_data["obs_steps_list"]);
  IntegerVector sim_end_step_vec = as<IntegerVector>(scenario_data["sim_end_step_vec"]);
  List obs_burden_list = as<List>(scenario_data["obs_burden_list"]);
  List keep_burden_list = as<List>(scenario_data["keep_burden_list"]);
  List ploidy_z_list = as<List>(scenario_data["ploidy_z_list"]);
  NumericVector init_mult_vec = scenario_data.containsElementNamed("init_mult_vec")
    ? as<NumericVector>(scenario_data["init_mult_vec"])
    : NumericVector(cohort_code.size(), 1.0);
  NumericVector obs_necrosis_fraction_vec = scenario_data.containsElementNamed("obs_necrosis_fraction_vec")
    ? as<NumericVector>(scenario_data["obs_necrosis_fraction_vec"])
    : NumericVector(cohort_code.size(), NA_REAL);
  LogicalVector keep_necrosis_vec = scenario_data.containsElementNamed("keep_necrosis_vec")
    ? as<LogicalVector>(scenario_data["keep_necrosis_vec"])
    : LogicalVector(cohort_code.size(), false);

  NumericVector mu_by_N = as<NumericVector>(objective_data["mu_by_N"]);
  double sigma_burden = as<double>(objective_data["sigma_burden"]);
  double sigma_ploidy = as<double>(objective_data["sigma_ploidy"]);
  bool use_necrosis_loss = objective_data.containsElementNamed("use_necrosis_loss")
    ? as<bool>(objective_data["use_necrosis_loss"])
    : false;
  double sigma_necrosis_logit = objective_data.containsElementNamed("sigma_necrosis_logit")
    ? as<double>(objective_data["sigma_necrosis_logit"])
    : 0.75;
  double necrosis_fraction_eps = objective_data.containsElementNamed("necrosis_fraction_eps")
    ? as<double>(objective_data["necrosis_fraction_eps"])
    : 1e-4;

  NumericVector init_state_2N = as<NumericVector>(state_data["init_state_2N"]);
  NumericVector init_state_4N = as<NumericVector>(state_data["init_state_4N"]);
  int N0min = as<int>(state_data["N0min"]);
  int N0max = as<int>(state_data["N0max"]);
  int N1min = as<int>(state_data["N1min"]);
  int N1max = as<int>(state_data["N1max"]);
  int N_unit = as<int>(state_data["N_unit"]);
  NumericVector vol_by_N = as<NumericVector>(state_data["vol_by_N"]);

  double DT = as<double>(sim_args["DT"]);
  double dose_ref = as<double>(sim_args["dose_ref"]);
  bool fit_treatment = as<bool>(sim_args["fit_treatment"]);
  double alpha = as<double>(sim_args["alpha"]);
  double gamma = as<double>(sim_args["gamma"]);
  double tx_mult_min = as<double>(sim_args["tx_mult_min"]);
  bool crowding_enabled = as<bool>(sim_args["crowding_enabled"]);
  std::string crowding = as<std::string>(sim_args["crowding"]);
  double K = as<double>(sim_args["K"]);
  double min_pop = as<double>(sim_args["min_pop"]);
  double O2_crit = as<double>(sim_args["O2_crit"]);
  bool o2_feedback = as<bool>(sim_args["o2_feedback"]);
  double o2_S0 = as<double>(sim_args["o2_S0"]);
  double kappa_O = as<double>(sim_args["kappa_O"]);
  double tau_O2 = as<double>(sim_args["tau_O2"]);
  double o2_Nref = as<double>(sim_args["o2_Nref"]);
  double o2_min = as<double>(sim_args["o2_min"]);
  double o2_S0_upper_bound = sim_args.containsElementNamed("o2_S0_upper_bound")
    ? as<double>(sim_args["o2_S0_upper_bound"])
    : 5.0;
  double eta_o2 = as<double>(sim_args["eta_o2"]);
  double o2_cache_bin_pct = as<double>(sim_args["o2_cache_bin_pct"]);
  double o2_cache_hysteresis_pct = as<double>(sim_args["o2_cache_hysteresis_pct"]);
  bool o2_cache_profile = as<bool>(sim_args["o2_cache_profile"]);
  double lam_max = as<double>(sim_args["lam_max"]);
  double p_mis_base = as<double>(sim_args["p_mis_base"]);
  double p_misseg = as<double>(sim_args["p_misseg"]);
  double k_o_mis = as<double>(sim_args["k_o_mis"]);
  double p_wgd = as<double>(sim_args["p_wgd"]);
  std::string boundary = as<std::string>(sim_args["boundary"]);
  double eps_tail = as<double>(sim_args["eps_tail"]);
  double buffer_smax = as<double>(sim_args["buffer_smax"]);
  double buffer_beta = as<double>(sim_args["buffer_beta"]);
  double buffer_n_exp = as<double>(sim_args["buffer_n_exp"]);
  double beta_size = as<double>(sim_args["beta_size"]);
  double alpha_o2 = as<double>(sim_args["alpha_o2"]);
  double gamma_growth = as<double>(sim_args["gamma_growth"]);
  double mu_hp = as<double>(sim_args["mu_hp"]);
  double gamma_mu = as<double>(sim_args["gamma_mu"]);
  double n_O = as<double>(sim_args["n_O"]);
  std::string ploidy_O2_death = as<std::string>(sim_args["ploidy_O2_death"]);
  std::string start_with = as<std::string>(sim_args["start_with"]);
  double k_clear = as<double>(sim_args["k_clear"]);
  double burden_log_eps = as<double>(sim_args["burden_log_eps"]);

  const int n_sc = cohort_code.size();
  if (dose_vec.size() != n_sc || treat_day_vec.size() != n_sc ||
      obs_steps_list.size() != n_sc || sim_end_step_vec.size() != n_sc ||
      obs_burden_list.size() != n_sc || keep_burden_list.size() != n_sc ||
      ploidy_z_list.size() != n_sc || init_mult_vec.size() != n_sc ||
      obs_necrosis_fraction_vec.size() != n_sc || keep_necrosis_vec.size() != n_sc) {
    stop("Scenario containers must have consistent length.");
  }

  const double log_eps_use =
    (std::isfinite(burden_log_eps) && burden_log_eps > 0.0) ? burden_log_eps : 1e-15;
  const double sigma_b_use =
    (std::isfinite(sigma_burden) && sigma_burden > 0.0) ? sigma_burden : 0.35;
  const double sigma_p_use =
    (std::isfinite(sigma_ploidy) && sigma_ploidy > 0.0) ? sigma_ploidy : 0.08;
  const double sigma_n_use =
    (std::isfinite(sigma_necrosis_logit) && sigma_necrosis_logit > 0.0) ? sigma_necrosis_logit : 0.75;
  const double necrosis_eps_use =
    (std::isfinite(necrosis_fraction_eps) && necrosis_fraction_eps > 0.0 && necrosis_fraction_eps < 0.5)
      ? necrosis_fraction_eps
      : 1e-4;
  const double prob_eps = 1e-300;
  const bool o2_growth_use = !(std::isfinite(alpha_o2) && alpha_o2 < 0.0);
  const double alpha_o2_use = std::fabs(alpha_o2);
  std::vector<double> burden_losses;
  std::vector<double> ploidy_losses_2N;
  std::vector<double> ploidy_losses_4N;
  std::vector<double> necrosis_losses;
  burden_losses.reserve(static_cast<size_t>(n_sc));
  ploidy_losses_2N.reserve(static_cast<size_t>(n_sc));
  ploidy_losses_4N.reserve(static_cast<size_t>(n_sc));
  necrosis_losses.reserve(static_cast<size_t>(n_sc));
  double burden_nll_total = 0.0;
  double ploidy_nll_total = 0.0;
  double necrosis_nll_total = 0.0;
  int burden_obs_total = 0;
  int ploidy_obs_total = 0;
  int necrosis_obs_total = 0;

  int cache_g_build = 0;
  int cache_g_hit = 0;
  int cache_g_hysteresis = 0;

  for (int i = 0; i < n_sc; ++i) {
    const int cohort = cohort_code[i];
    NumericVector init_state = (cohort == 0) ? init_state_2N : init_state_4N;
    const double init_mult = (std::isfinite(init_mult_vec[i]) && init_mult_vec[i] > 0.0) ? init_mult_vec[i] : 1.0;
    if (init_mult != 1.0) {
      init_state = clone(init_state);
      for (int j = 0; j < init_state.size(); ++j) {
        init_state[j] = init_state[j] * init_mult;
      }
    }
    IntegerVector obs_steps = as<IntegerVector>(obs_steps_list[i]);
    NumericVector obs_burden = as<NumericVector>(obs_burden_list[i]);
    LogicalVector keep_day = as<LogicalVector>(keep_burden_list[i]);
    NumericVector ploidy_z = as<NumericVector>(ploidy_z_list[i]);

    List sim_one_args = clone(sim_args);
    sim_one_args["init_state"] = init_state;
    sim_one_args["N0min"] = N0min;
    sim_one_args["N0max"] = N0max;
    sim_one_args["N1min"] = N1min;
    sim_one_args["N1max"] = N1max;
    sim_one_args["obs_steps"] = obs_steps;
    sim_one_args["sim_end_step"] = sim_end_step_vec[i];
    sim_one_args["dose"] = dose_vec[i];
    sim_one_args["treat_day"] = treat_day_vec[i];
    sim_one_args["fit_treatment"] = fit_treatment;
    sim_one_args["alpha"] = alpha;
    sim_one_args["gamma"] = gamma;
    sim_one_args["tx_mult_min"] = tx_mult_min;
    sim_one_args["crowding_enabled"] = crowding_enabled;
    sim_one_args["crowding"] = crowding;
    sim_one_args["K"] = K;
    sim_one_args["min_pop"] = min_pop;
    sim_one_args["O2_crit"] = O2_crit;
    sim_one_args["o2_feedback"] = o2_feedback;
    sim_one_args["o2_S0"] = o2_S0;
    sim_one_args["kappa_O"] = kappa_O;
    sim_one_args["tau_O2"] = tau_O2;
    sim_one_args["o2_Nref"] = o2_Nref;
    sim_one_args["o2_min"] = o2_min;
    sim_one_args["o2_S0_upper_bound"] = o2_S0_upper_bound;
    sim_one_args["eta_o2"] = eta_o2;
    sim_one_args["o2_cache_bin_pct"] = o2_cache_bin_pct;
    sim_one_args["o2_cache_hysteresis_pct"] = o2_cache_hysteresis_pct;
    sim_one_args["o2_cache_profile"] = o2_cache_profile;
    sim_one_args["lam_max"] = lam_max;
    sim_one_args["p_mis_base"] = p_mis_base;
    sim_one_args["p_misseg"] = p_misseg;
    sim_one_args["k_o_mis"] = k_o_mis;
    sim_one_args["p_wgd"] = p_wgd;
    sim_one_args["boundary"] = boundary;
    sim_one_args["eps_tail"] = eps_tail;
    sim_one_args["buffer_smax"] = buffer_smax;
    sim_one_args["buffer_beta"] = buffer_beta;
    sim_one_args["buffer_n_exp"] = buffer_n_exp;
    sim_one_args["N_unit"] = N_unit;
    sim_one_args["beta_size"] = beta_size;
    sim_one_args["O2_growth"] = o2_growth_use;
    sim_one_args["alpha_o2"] = alpha_o2_use;
    sim_one_args["gamma_growth"] = gamma_growth;
    sim_one_args["mu_hp"] = mu_hp;
    sim_one_args["gamma_mu"] = gamma_mu;
    sim_one_args["n_O"] = n_O;
    sim_one_args["ploidy_O2_death"] = ploidy_O2_death;
    sim_one_args["start_with"] = start_with;
    sim_one_args["k_clear"] = k_clear;
    sim_one_args["vol_by_N"] = vol_by_N;
    sim_one_args["burden_floor"] = log_eps_use;
    sim_one_args["return_full_trajectory"] = false;

    List sim = cpp_o2simps_simulate_one(sim_one_args);

    NumericVector pred_burden = sim["Vmm3_total_obs"];
    NumericVector frac_N = sim["frac_N_live"];
    cache_g_build += as<int>(sim["cache_g_build"]);
    cache_g_hit += as<int>(sim["cache_g_hit"]);
    cache_g_hysteresis += as<int>(sim["cache_g_hysteresis"]);

    if (mu_by_N.size() != frac_N.size()) {
      stop("mu_by_N length must equal simulated terminal state vector length.");
    }

    // Burden log-normal NLL per tumor (mean across available time points).
    const int nb = std::min(obs_burden.size(), pred_burden.size());
    double burden_nll_sum = 0.0;
    int burden_n = 0;
    for (int j = 0; j < nb; ++j) {
      const bool keepj = (keep_day.size() == nb) ? static_cast<bool>(keep_day[j]) : true;
      if (!keepj) continue;
      const double obs = obs_burden[j];
      const double pred = pred_burden[j];
      if (!std::isfinite(obs) || !std::isfinite(pred) || obs <= 0.0 || pred <= 0.0) continue;
      const double resid = std::log(std::max(obs, log_eps_use)) - std::log(std::max(pred, log_eps_use));
      const double z = resid / sigma_b_use;
      burden_nll_sum += std::log(sigma_b_use) + 0.5 * std::log(2.0 * 3.14159265358979323846) + 0.5 * z * z;
      ++burden_n;
    }
    if (burden_n > 0) {
      burden_nll_total += burden_nll_sum;
      burden_obs_total += burden_n;
      burden_losses.push_back(burden_nll_sum / static_cast<double>(burden_n));
    }

    if (use_necrosis_loss && static_cast<bool>(keep_necrosis_vec[i])) {
      const double obs_necrosis = obs_necrosis_fraction_vec[i];
      const double pred_dead = as<double>(sim["Vmm3_dead_total_terminal"]);
      const double pred_total = as<double>(sim["Vmm3_total_terminal"]);
      if (std::isfinite(obs_necrosis) && std::isfinite(pred_dead) &&
          std::isfinite(pred_total) && pred_total > 0.0) {
        const double pred_necrosis = pred_dead / pred_total;
        const double resid =
          logit_prob_eps(pred_necrosis, necrosis_eps_use) -
          logit_prob_eps(obs_necrosis, necrosis_eps_use);
        const double z = resid / sigma_n_use;
        if (std::isfinite(z)) {
          const double necrosis_loss = z * z;
          necrosis_nll_total += necrosis_loss;
          ++necrosis_obs_total;
          necrosis_losses.push_back(necrosis_loss);
        }
      }
    }

    // Continuous single-cell ploidy NLL per tumor:
    // p(z) = sum_j pi_j * Normal(z; mu_j, sigma_ploidy^2), then average -log p(z).
    double p_sum = 0.0;
    for (int j = 0; j < frac_N.size(); ++j) {
      const double pv = frac_N[j];
      if (std::isfinite(pv) && pv > 0.0) p_sum += pv;
    }
    if (ploidy_z.size() > 0 && p_sum > 0.0) {
      double ploidy_nll_sum = 0.0;
      int ploidy_n = 0;
      for (int c = 0; c < ploidy_z.size(); ++c) {
        const double z_obs = ploidy_z[c];
        if (!std::isfinite(z_obs)) continue;
        double max_log = -std::numeric_limits<double>::infinity();
        std::vector<double> comp_log;
        comp_log.reserve(static_cast<size_t>(frac_N.size()));
        for (int j = 0; j < frac_N.size(); ++j) {
          const double pv = frac_N[j];
          if (!std::isfinite(pv) || pv <= 0.0) continue;
          const double pj = pv / p_sum;
          const double mu_j = mu_by_N[j];
          if (!std::isfinite(mu_j)) continue;
          const double log_comp =
            std::log(std::max(pj, prob_eps)) + R::dnorm4(z_obs, mu_j, sigma_p_use, /*give_log=*/1);
          comp_log.push_back(log_comp);
          if (log_comp > max_log) max_log = log_comp;
        }
        if (comp_log.empty() || !std::isfinite(max_log)) continue;
        double sum_exp = 0.0;
        for (size_t t = 0; t < comp_log.size(); ++t) {
          sum_exp += std::exp(comp_log[t] - max_log);
        }
        if (!std::isfinite(sum_exp) || sum_exp <= 0.0) continue;
        const double log_mix = max_log + std::log(sum_exp);
        ploidy_nll_sum += -log_mix;
        ++ploidy_n;
      }
      if (ploidy_n > 0) {
        ploidy_nll_total += ploidy_nll_sum;
        ploidy_obs_total += ploidy_n;
        const double tumor_ploidy_loss = ploidy_nll_sum / static_cast<double>(ploidy_n);
        if (cohort == 0) {
          ploidy_losses_2N.push_back(tumor_ploidy_loss);
        } else {
          ploidy_losses_4N.push_back(tumor_ploidy_loss);
        }
      }
    }
  }

  const double L_b = burden_losses.empty()
    ? 0.0
    : std::accumulate(burden_losses.begin(), burden_losses.end(), 0.0) /
        static_cast<double>(burden_losses.size());
  const double L_n = necrosis_losses.empty()
    ? 0.0
    : std::accumulate(necrosis_losses.begin(), necrosis_losses.end(), 0.0) /
        static_cast<double>(necrosis_losses.size());
  const bool has_2N = !ploidy_losses_2N.empty();
  const bool has_4N = !ploidy_losses_4N.empty();
  const double L_p_2N = has_2N
    ? std::accumulate(ploidy_losses_2N.begin(), ploidy_losses_2N.end(), 0.0) /
        static_cast<double>(ploidy_losses_2N.size())
    : 0.0;
  const double L_p_4N = has_4N
    ? std::accumulate(ploidy_losses_4N.begin(), ploidy_losses_4N.end(), 0.0) /
        static_cast<double>(ploidy_losses_4N.size())
    : 0.0;
  const double L_p = (has_2N && has_4N)
    ? (0.5 * L_p_2N + 0.5 * L_p_4N)
    : (has_2N ? L_p_2N : (has_4N ? L_p_4N : 0.0));
  const int n_ploidy_total = static_cast<int>(ploidy_losses_2N.size() + ploidy_losses_4N.size());
  const double objective_burden_neg2loglik_raw = 2.0 * burden_nll_total;
  const double objective_ploidy_neg2loglik_raw = 2.0 * ploidy_nll_total;
  const double objective_necrosis_neg2loglik_raw = 2.0 * necrosis_nll_total;

  return List::create(
    _["L_b"] = L_b,
    _["L_p"] = L_p,
    _["L_n"] = L_n,
    _["burden_nll_total"] = burden_nll_total,
    _["ploidy_nll_total"] = ploidy_nll_total,
    _["necrosis_nll_total"] = necrosis_nll_total,
    _["objective_burden_neg2loglik_raw"] = objective_burden_neg2loglik_raw,
    _["objective_ploidy_neg2loglik_raw"] = objective_ploidy_neg2loglik_raw,
    _["objective_necrosis_neg2loglik_raw"] = objective_necrosis_neg2loglik_raw,
    _["n_burden"] = static_cast<int>(burden_losses.size()),
    _["n_burden_obs_total"] = burden_obs_total,
    _["n_ploidy"] = n_ploidy_total,
    _["n_ploidy_obs_total"] = ploidy_obs_total,
    _["n_necrosis"] = static_cast<int>(necrosis_losses.size()),
    _["n_necrosis_obs_total"] = necrosis_obs_total,
    _["n_ploidy_2N"] = static_cast<int>(ploidy_losses_2N.size()),
    _["n_ploidy_4N"] = static_cast<int>(ploidy_losses_4N.size()),
    _["L_p_2N"] = L_p_2N,
    _["L_p_4N"] = L_p_4N,
    _["cache_g_build"] = cache_g_build,
    _["cache_g_hit"] = cache_g_hit,
    _["cache_g_hysteresis"] = cache_g_hysteresis
  );
}

)---",
  env = globalenv(), rebuild = FALSE, showOutput = FALSE, verbose = FALSE,
  cacheDir = .RESPONSE_CPP_CACHE
)
.USE_CPP_O2SIMPS_BACKEND <- TRUE

# ---- embedded analytical model functions --------------------------------
o2simps_cpp_dll_info <- function() {
  if (!isTRUE(.USE_CPP_O2SIMPS_BACKEND)) {
    stop("C++ backend is not initialized.")
  }
  if (!exists("cpp_o2simps_simulate_one", mode = "function", inherits = TRUE)) {
    stop("Required wrapper function missing: cpp_o2simps_simulate_one")
  }

  loaded <- getLoadedDLLs()
  if (length(loaded) == 0L) stop("No loaded DLLs found after sourceCpp initialization.")
  dll_names <- names(loaded)
  dll_paths <- vapply(loaded, function(x) as.character(x[["path"]]), FUN.VALUE = character(1), USE.NAMES = FALSE)
  valid <- nzchar(dll_paths) & file.exists(dll_paths)

  # Prefer DLLs from this model's dedicated sourceCpp cache.
  cache_pat <- ".rcpp_cache_o2_supply_demand_map"
  in_cache <- valid & grepl(cache_pat, dll_paths, fixed = TRUE)
  candidate <- if (any(in_cache)) in_cache else (valid & grepl("sourceCpp", dll_names, fixed = TRUE))
  if (!any(candidate)) {
    stop("Unable to resolve sourceCpp DLL path from loaded DLL list.")
  }

  cand_idx <- which(candidate)
  mt <- suppressWarnings(file.info(dll_paths[cand_idx])$mtime)
  best <- cand_idx[[1]]
  if (length(cand_idx) > 1L && any(is.finite(as.numeric(mt)))) {
    ord <- order(mt, decreasing = TRUE, na.last = TRUE)
    best <- cand_idx[[ord[[1]]]]
  }
  dll_path <- normalizePath(dll_paths[[best]], mustWork = TRUE)
  wrapper_candidates <- list.files(
    dirname(dll_path),
    pattern = "\\.cpp\\.R$",
    full.names = TRUE
  )
  wrapper_path <- ""
  if (length(wrapper_candidates) > 0L) {
    if (length(wrapper_candidates) == 1L) {
      wrapper_path <- wrapper_candidates[[1]]
    } else {
      hit <- vapply(
        wrapper_candidates,
        function(p) {
          txt <- tryCatch(readLines(p, warn = FALSE), error = function(e) character(0))
          any(grepl("cpp_o2simps_simulate_one", txt, fixed = TRUE))
        },
        logical(1)
      )
      if (any(hit)) wrapper_path <- wrapper_candidates[which(hit)[1]]
    }
  }
  if (!nzchar(wrapper_path) || !file.exists(wrapper_path)) {
    stop("Unable to resolve sourceCpp wrapper file (*.cpp.R) for O2_supply_demand_MAP backend.")
  }

  list(
    name = as.character(dll_names[[best]]),
    path = dll_path,
    wrapper_path = normalizePath(wrapper_path, mustWork = TRUE)
  )
}

# -----------------------------------------------------------------------------
# Function: .first_non_null
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - ...: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

# Euler stepper for generator matrix dynamics.
# -----------------------------------------------------------------------------
# Function: step_dt
# Purpose: Advance state vector by repeated generator-matrix time steps.
# Parameters:
#   - G: Generator or transition matrix used for state propagation.
#   - x: Input value or vector to process.
#   - dt: Time-step length in days.
#   - steps: Number of repeated matrix-step applications.
#   - normalize: Logical flag to normalize resulting state distribution.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
step_dt <- function(G, x, dt, steps = 1L, normalize = FALSE) {
  I <- Diagonal(n = nrow(G))
  A <- I + dt * G
  v <- as.numeric(x)
  for (k in seq_len(steps)) {
    v <- as.numeric(A %*% v)
    if (normalize) v <- v / sum(v)
  }
  v
}

# Chromosome-length-weighted ploidy mapping on fixed autosomes 1..22.
# -----------------------------------------------------------------------------
# Function: default_chr_lengths_bp_1to22
# Purpose: Return fixed chromosome lengths (bp) for autosomes 1..22.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Numeric vector of chromosome lengths in base pairs.
# -----------------------------------------------------------------------------
default_chr_lengths_bp_1to22 <- function() {
  c(
    `1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
    `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
    `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
    `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
    `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
    `21` = 46709983, `22` = 50818468
  )
}

# Precompute immutable defaults once; these are reused throughout the fit/eval path.
.chr_lengths_default_bp_1to22 <- default_chr_lengths_bp_1to22()
.chr_lengths_default_ord_desc <- order(.chr_lengths_default_bp_1to22, decreasing = TRUE)
.chr_lengths_default_denom <- sum(.chr_lengths_default_bp_1to22)

# -----------------------------------------------------------------------------
# Function: normalize_chr_lengths_bp_1to22
# Purpose: Validate and normalize chromosome-length vector for autosomes 1..22.
# Parameters:
#   - chr_lengths_bp: Optional chromosome-length vector.
# Returns:
#   Numeric vector of length 22 with positive finite values.
# -----------------------------------------------------------------------------
normalize_chr_lengths_bp_1to22 <- function(chr_lengths_bp = NULL) {
  w <- if (is.null(chr_lengths_bp)) default_chr_lengths_bp_1to22() else as.numeric(chr_lengths_bp)
  if (length(w) != 22L) stop("chr_lengths_bp must have length 22 (autosomes 1..22).")
  if (any(!is.finite(w) | w <= 0)) stop("chr_lengths_bp must be all positive finite values.")
  names(w) <- as.character(seq_len(22L))
  w
}

# -----------------------------------------------------------------------------
# Function: weighted_ploidy_from_total_N
# Purpose: Map total copy-number state N to chromosome-length-weighted ploidy.
# Parameters:
#   - N_total: Total copy number state(s) on the integer grid.
#   - chr_lengths_bp: Optional chromosome-length vector.
# Returns:
#   Numeric weighted ploidy values.
# -----------------------------------------------------------------------------
weighted_ploidy_from_total_N <- function(N_total, chr_lengths_bp = NULL) {
  if (is.null(chr_lengths_bp)) {
    w <- .chr_lengths_default_bp_1to22
    ord <- .chr_lengths_default_ord_desc
    n_chr <- 22L
    denom <- .chr_lengths_default_denom
  } else {
    w <- normalize_chr_lengths_bp_1to22(chr_lengths_bp)
    ord <- order(w, decreasing = TRUE)
    n_chr <- length(w)
    denom <- sum(w)
  }
  Nv <- as.numeric(N_total)
  vapply(Nv, function(nn) {
    if (!is.finite(nn)) return(NA_real_)
    n_int <- as.integer(round(nn))
    if (n_int < 0L) n_int <- 0L
    base <- n_int %/% n_chr
    rem <- n_int %% n_chr
    cn <- rep.int(base, n_chr)
    if (rem > 0L) cn[ord[seq_len(rem)]] <- cn[ord[seq_len(rem)]] + 1L
    sum(cn * w) / denom
  }, numeric(1))
}

# -----------------------------------------------------------------------------
# Function: map_ploidy_to_N_by_chrlen
# Purpose: Map ploidy value(s) to nearest integer N state under weighted mapping.
# Parameters:
#   - ploidy_values: Observed ploidy value(s).
#   - N_grid: Integer N grid.
#   - chr_lengths_bp: Optional chromosome-length vector.
# Returns:
#   Integer N states on N_grid.
# -----------------------------------------------------------------------------
map_ploidy_to_N_by_chrlen <- function(ploidy_values, N_grid, chr_lengths_bp = NULL) {
  grid <- as.integer(sort(unique(N_grid)))
  p_grid <- weighted_ploidy_from_total_N(grid, chr_lengths_bp = chr_lengths_bp)
  pv <- as.numeric(ploidy_values)
  vapply(pv, function(p) {
    if (!is.finite(p)) return(NA_integer_)
    d <- abs(p_grid - p)
    k <- which.min(d)
    as.integer(grid[[k]])
  }, integer(1))
}

# Build a normalized initial distribution on an integer N grid from ploidy values.
# -----------------------------------------------------------------------------
# Function: create_initial_dist
# Purpose: Construct initial ploidy-state distribution/state vector for simulation.
# Parameters:
#   - ploidy_values: Function-specific input argument.
#   - N_grid: Function-specific input argument.
#   - N_unit: Number of modeled chromosome classes.
#   - chr_lengths_bp: Optional chromosome-length vector for weighted ploidy mapping.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
create_initial_dist <- function(ploidy_values, N_grid, N_unit = 22L, chr_lengths_bp = NULL) {
  N_values <- map_ploidy_to_N_by_chrlen(
    ploidy_values = ploidy_values,
    N_grid = N_grid,
    chr_lengths_bp = chr_lengths_bp
  )
  N_counts <- table(N_values)
  N_fracs <- as.numeric(N_counts) / sum(N_counts)
  names(N_fracs) <- names(N_counts)
  x_vec <- rep(0, length(N_grid))
  names(x_vec) <- N_grid
  valid_names <- names(N_fracs)[names(N_fracs) %in% names(x_vec)]
  x_vec[valid_names] <- N_fracs[valid_names]
  x_vec
}

# Build explicit initial state vectors on a single chromosome-count grid.
# - Either pass ploidy or a custom initial fraction vector.
# - Returns absolute cell counts with total mass = total_size.
# -----------------------------------------------------------------------------
# Function: make_init_state
# Purpose: Construct initial ploidy-state distribution/state vector for simulation.
# Parameters:
#   - grid_pre: Chromosome-count grid.
#   - ploidy: Initial ploidy label used to choose nearest N state.
#   - init_dist: Optional named initial fractions on grid_pre.
#   - N_UNIT: Function-specific input argument.
#   - total_size: Function-specific input argument.
#   - chr_lengths_bp: Optional chromosome-length vector for weighted ploidy mapping.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
make_init_state <- function(grid_pre,
                            ploidy = c(2, 4),
                            init_dist = NULL,
                            N_UNIT = 22L, total_size = 1e6, chr_lengths_bp = NULL) {
  # Single-layer initialization: all cohorts are initialized on one N grid.
  ploidy <- match.arg(as.character(ploidy), choices = c("2", "4"))
  Pnum <- as.numeric(ploidy)
  grid_N <- as.integer(grid_pre)
  x <- rep(0, length(grid_N))
  names(x) <- grid_N

  if (!is.null(init_dist)) {
    nm <- intersect(names(init_dist), names(x))
    x[nm] <- x[nm] + as.numeric(init_dist[nm])
  } else {
    N_delta <- as.integer(map_ploidy_to_N_by_chrlen(
      ploidy_values = Pnum,
      N_grid = grid_N,
      chr_lengths_bp = chr_lengths_bp
    ))
    stopifnot(N_delta %in% grid_N)
    x[as.character(N_delta)] <- 1
  }

  s <- sum(x)
  if (s <= 0) stop("Init mass is zero.")
  x / s * total_size
}

# -----------------------------------------------------------------------------
# Function: .clip01
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.clip01 <- function(x) pmin(pmax(x, 0), 1)
# -----------------------------------------------------------------------------
# Function: .clip_o2pct
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.clip_o2pct <- function(x) pmin(pmax(x, 0), 100)
# -----------------------------------------------------------------------------
# Function: .assert_o2_pct
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
#   - label: Text label used for logging and progress messages.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.assert_o2_pct <- function(x, label = "O2") {
  x_num <- as.numeric(x)
  if (any(!is.finite(x_num))) stop(label, " must be finite.")
  if (any(x_num < 0 | x_num > 100)) stop(label, " must be in percent scale [0, 100].")
  x_num
}

# -----------------------------------------------------------------------------
# Function: .require_cpp_o2simps_fn
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - fn_name: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.require_cpp_o2simps_fn <- function(fn_name) {
  if (!(isTRUE(.USE_CPP_O2SIMPS_BACKEND) &&
        exists(fn_name, mode = "function", inherits = TRUE))) {
    stop("Required C++ backend function is unavailable: ", fn_name)
  }
}

# O2 target from a logarithmic supply-demand model on effective demand:
#   O2_target = max(o2_min, o2_S0 - kappa_O * log(1 + N_eff / o2_Nref)),
# where N_eff can encode ploidy-weighted demand.
# Inputs:
# - Ntot: nonnegative effective oxygen-demand proxy, often ploidy-weighted live burden
# - run_params: natural-scale oxygen parameters, including the active o2_S0_upper_bound
# - o2_Nref: demand normalization constant on the same scale as Ntot
#
# Returns:
# - oxygen target percentages after the same active-cap normalization enforced in C++
.o2_supply_demand_from_burden <- function(Ntot,
                                          run_params,
                                          o2_Nref = 1e6) {
  Ntot_use <- pmax(as.numeric(Ntot), 0)
  o2_s0_upper_use <- as.numeric(.first_non_null(run_params$o2_S0_upper_bound, run_params$o2_S0_max, 5.0))
  if (!is.finite(o2_s0_upper_use) || o2_s0_upper_use <= 0) o2_s0_upper_use <- 5.0
  o2_S0 <- as.numeric(.first_non_null(run_params$o2_S0, 0.5))
  kappa_O <- as.numeric(.first_non_null(run_params$kappa_O, 1.0))
  o2_min <- as.numeric(.first_non_null(run_params$o2_min, 0.0))
  if (!is.finite(kappa_O) || kappa_O <= 0) kappa_O <- 1.0
  if (!is.finite(o2_min) || o2_min < 0) o2_min <- 0.0
  if (!is.finite(o2_Nref) || o2_Nref <= 0) o2_Nref <- 1e6
  o2_S0 <- max(0, min(o2_s0_upper_use, o2_S0))
  o2_min <- max(0, min(o2_s0_upper_use, o2_min))

  .require_cpp_o2simps_fn("cpp_o2simps_o2_window_supply")
  return(as.numeric(cpp_o2simps_o2_window_supply(
    Ntot = as.numeric(Ntot_use),
    o2_S0 = as.numeric(o2_S0),
    kappa_O = as.numeric(kappa_O),
    o2_Nref = as.numeric(o2_Nref),
    o2_min = as.numeric(o2_min),
    o2_S0_upper_bound = as.numeric(o2_s0_upper_use)
  )))
}

# Richard buffering.R-style lambda: maximal baseline proliferation rate.
# -----------------------------------------------------------------------------
# Function: growth_lambda
# Purpose: Compute baseline proliferation rate for a given ploidy state.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - N: Ploidy state value or chromosome-copy count.
#   - lam_max: Maximal proliferation rate.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
growth_lambda <- function(O2, N, lam_max) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  invisible(O2_use)
  lam_max_use <- as.numeric(lam_max)
  if (!is.finite(lam_max_use)) stop("lam_max must be finite.")
  rep(pmax(lam_max_use, 0), length(N))
}

# Main-path oxygen-only death helper (aligned with C++).
# -----------------------------------------------------------------------------
# Function: .mu_eff_of_O2
# Purpose: Compute state-specific hypoxia death rate under the main model.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.resource_stress_of_O2 <- function(O2, run_params, O2_crit = NULL) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  O2_crit_use <- as.numeric(.first_non_null(O2_crit, run_params$O2_crit, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  n_O <- as.numeric(.first_non_null(run_params$n_O, 1.0))
  if (!is.finite(n_O) || n_O < 0) {
    stop("run_params$n_O must be finite and >= 0.")
  }
  o2_c <- pmax(O2_crit_use, 1e-12)
  h_o2 <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(as.numeric(O2_use), 0)^n_O))
  .clip01(h_o2)
}

# Main-path proliferation helper aligned with the current C++ runtime dispatch.
# -----------------------------------------------------------------------------
# Function: .lambda_eff_of_O2
# Purpose: Compute state-specific effective proliferation under the active model.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
#   - O2_growth: Optional runtime override for the growth-penalty switch.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.lambda_eff_of_O2 <- function(O2, run_params, N = 44, O2_crit = NULL, O2_growth = TRUE) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  N_use <- as.numeric(N)
  if (any(!is.finite(N_use))) stop("N must be finite.")
  n_out <- max(length(O2_use), length(N_use))
  if (!(length(O2_use) %in% c(1L, n_out) &&
        length(N_use) %in% c(1L, n_out))) {
    stop("O2 and N must have compatible lengths.")
  }
  O2_vec <- rep_len(as.numeric(O2_use), n_out)
  N_vec <- rep_len(N_use, n_out)
  O2_crit_use <- as.numeric(.first_non_null(O2_crit, run_params$O2_crit, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  n_O <- as.numeric(.first_non_null(run_params$n_O, 1.0))
  if (!is.finite(n_O) || n_O < 0) stop("run_params$n_O must be finite and >= 0.")
  lam_max_use <- as.numeric(.first_non_null(run_params$lam_max, 0.0))
  if (!is.finite(lam_max_use)) lam_max_use <- 0.0
  lam_base <- rep(pmax(lam_max_use, 0), n_out)
  alpha_o2_use <- pmax(as.numeric(.first_non_null(run_params$alpha_o2, 0.0)), 0)
  gamma_growth_use <- pmax(as.numeric(.first_non_null(run_params$gamma_growth, 1.0)), 1e-12)
  o2_c <- pmax(O2_crit_use, 1e-12)
  h_o2 <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(O2_vec, 0)^n_O))
  h_o2 <- .clip01(h_o2)

  if (!isTRUE(O2_growth)) return(pmax(lam_base, 0))
  denom <- 1 + alpha_o2_use * h_o2 * ((pmax(N_vec, 0) / 44)^gamma_growth_use)
  pmax(lam_base / pmax(denom, 1e-12), 0)
}

# Main-path baseline-plus-increment missegregation helper (aligned with C++).
# -----------------------------------------------------------------------------
# Function: .mu_eff_of_O2
# Purpose: Compute state-specific hypoxia death rate under the main model.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.mu_eff_of_O2 <- function(O2, run_params, N = 44, O2_crit = NULL) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  N_use <- as.numeric(N)
  if (any(!is.finite(N_use))) stop("N must be finite.")
  n_out <- max(length(O2_use), length(N_use))
  if (!(length(O2_use) %in% c(1L, n_out) &&
        length(N_use) %in% c(1L, n_out))) {
    stop("O2 and N must have compatible lengths.")
  }
  O2_vec <- rep_len(as.numeric(O2_use), n_out)
  N_vec <- rep_len(N_use, n_out)

  O2_crit_use <- as.numeric(.first_non_null(O2_crit, run_params$O2_crit, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  n_O <- as.numeric(.first_non_null(run_params$n_O, 1.0))
  if (!is.finite(n_O) || n_O < 0) stop("run_params$n_O must be finite and >= 0.")
  o2_c <- pmax(O2_crit_use, 1e-12)
  h_o2 <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(O2_vec, 0)^n_O))
  h_o2 <- .clip01(h_o2)

  mu_hp_use <- as.numeric(.first_non_null(run_params$mu_hp, 0.0))
  if (!is.finite(mu_hp_use) || mu_hp_use < 0) mu_hp_use <- 0.0
  gamma_mu_use <- as.numeric(.first_non_null(run_params$gamma_mu, 1.0))
  if (!is.finite(gamma_mu_use) || gamma_mu_use <= 0) gamma_mu_use <- 1.0
  ploidy_O2_death_mode <- assert_canonical_ploidy_o2_death_mode(
    .first_non_null(run_params$ploidy_O2_death, "diploid_NULL")
  )

  if (identical(ploidy_O2_death_mode, "uniform")) {
    mu_eff <- mu_hp_use * h_o2
  } else if (identical(ploidy_O2_death_mode, "diploid_NULL")) {
    above_dip <- pmax(N_vec / 44.0 - 1.0, 0.0)
    mu_eff <- mu_hp_use * h_o2 * (1.0 + (above_dip^gamma_mu_use))
  } else {
    ratio <- pmax(N_vec / 44.0, 0.0)
    mu_eff <- mu_hp_use * h_o2 * (ratio^gamma_mu_use)
  }
  pmax(mu_eff, 0.0)
}

# Main-path death-linked O2/ploidy missegregation helper (aligned with C++).
# -----------------------------------------------------------------------------
# Function: .pmisseg_of_O2
# Purpose: Compute state-specific missegregation probability (baseline + death-linked increment).
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.pmisseg_of_O2 <- function(O2, run_params, N = 44, O2_crit = NULL) {
  mu_eff <- .mu_eff_of_O2(
    O2 = O2,
    run_params = run_params,
    N = N,
    O2_crit = O2_crit
  )

  p_base <- as.numeric(.first_non_null(run_params$p_mis_base, 1e-5))
  if (!is.finite(p_base) || p_base < 0) p_base <- 1e-5
  p_amp <- as.numeric(.first_non_null(run_params$p_misseg, 0.0))
  if (!is.finite(p_amp) || p_amp < 0) p_amp <- 0.0
  k_o_mis <- as.numeric(.first_non_null(run_params$k_o_mis, 50.0))
  if (!is.finite(k_o_mis) || k_o_mis <= 0) k_o_mis <- 1e-12
  frac <- mu_eff / (mu_eff + k_o_mis)
  delta_p <- p_amp * frac
  .clip01(p_base + delta_p)
}

# -----------------------------------------------------------------------------
# Function: .p_wgd_of_O2
# Purpose: Return the constant per-division WGD probability under the main model.
# Parameters:
#   - O2: Oxygen level used only to match diagnostic vector lengths.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
# Returns:
#   Numeric vector used by downstream diagnostics and plotting helpers.
# -----------------------------------------------------------------------------
.p_wgd_of_O2 <- function(O2, run_params) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  p_wgd_use <- as.numeric(.first_non_null(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd_use) || p_wgd_use < 0) p_wgd_use <- 0.0
  rep(.clip01(p_wgd_use), length(O2_use))
}

# Intrinsic-buffer delta weight formula (aligned with asymmetric_intrinsic_buffer).
# -----------------------------------------------------------------------------
# Function: .pr_delta_vec
# Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
# Parameters:
#   - N: Ploidy state value or chromosome-copy count.
#   - p: Missegregation probability parameter.
#   - eps_tail: Small truncation threshold for tail probabilities.
#   - buffer_smax: Maximum per-copy survival factor.
#   - buffer_beta: Ploidy-buffering strength.
#   - buffer_n_exp: Ploidy-buffering exponent.
#   - N_unit: Number of modeled chromosome classes for buffering scale.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.pr_delta_vec <- function(N, p, eps_tail = 1e-8, mr_lethality = 0.9,
                          N_unit = 22L,
                          buffer_smax = 1.0, buffer_beta = 0.0,
                          buffer_n_exp = 1.0) {
  .require_cpp_o2simps_fn("cpp_o2simps_pr_delta_vec")
  res <- cpp_o2simps_pr_delta_vec(
    as.integer(N),
    as.numeric(p),
    eps_tail = as.numeric(eps_tail),
    buffer_smax = as.numeric(buffer_smax),
    buffer_beta = as.numeric(buffer_beta),
    buffer_n_exp = as.numeric(buffer_n_exp),
    N_unit = as.integer(N_unit)
  )
  out <- as.numeric(res$prob)
  names(out) <- as.character(res$ts)
  attr(out, "mass_dropped") <- as.numeric(res$mass_dropped)
  return(out)
}

# -----------------------------------------------------------------------------
# Function: .build_B_total
# Purpose: Build total missegregation transition operator on ploidy grid.
# Parameters:
#   - Nmin: Minimum ploidy state on source grid.
#   - Nmax: Maximum ploidy state on source grid.
#   - p_vec: State-specific missegregation probability vector.
#   - mr_lethality: Probability of lethal outcome after severe missegregation.
#   - boundary: Boundary handling mode when transitions leave the ploidy grid.
#   - eps_tail: Small truncation threshold for tail probabilities.
#   - return_sparse: Function-specific input argument.
#   - buffer_smax: Maximum per-copy survival factor.
#   - buffer_beta: Ploidy-buffering strength.
#   - buffer_n_exp: Ploidy-buffering exponent.
#   - N_unit: Number of modeled chromosome classes for buffering scale.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_B_total <- function(Nmin, Nmax, p_vec, mr_lethality = 0.9,
                           boundary = c("drop", "absorb_minmax"),
                           eps_tail = 1e-8, return_sparse = TRUE,
                           N_unit = 22L,
                           buffer_smax = 1.0, buffer_beta = 0.0,
                           buffer_n_exp = 1.0) {
  boundary <- match.arg(boundary)
  R <- Nmax - Nmin + 1L
  if (length(p_vec) == 1L) p_vec <- rep(p_vec, R)

  .require_cpp_o2simps_fn("cpp_o2simps_build_B_total_triplet")
  tri <- cpp_o2simps_build_B_total_triplet(
    as.integer(Nmin),
    as.integer(Nmax),
    as.numeric(p_vec),
    boundary = boundary,
    eps_tail = as.numeric(eps_tail),
    buffer_smax = as.numeric(buffer_smax),
    buffer_beta = as.numeric(buffer_beta),
    buffer_n_exp = as.numeric(buffer_n_exp),
    N_unit = as.integer(N_unit)
  )
  B <- sparseMatrix(
    i = as.integer(tri$i),
    j = as.integer(tri$j),
    x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
    repr = "C"
  )
  return(if (isTRUE(return_sparse)) B else as.matrix(B))
}

# -----------------------------------------------------------------------------
# Function: .build_B_WGD
# Purpose: Build WGD transition operator between source and doubled-ploidy grids.
# Parameters:
#   - N0min: Minimum ploidy state on source grid.
#   - N0max: Maximum ploidy state on source grid.
#   - N1min: Minimum ploidy state on target grid for doubled states.
#   - N1max: Maximum ploidy state on target grid for doubled states.
#   - boundary: Boundary handling mode when transitions leave the ploidy grid.
#   - return_sparse: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_B_WGD <- function(N0min, N0max, N1min, N1max,
                         boundary = c("drop", "absorb_minmax"),
                         return_sparse = TRUE) {
  boundary <- match.arg(boundary)
  R0 <- N0max - N0min + 1L
  R1 <- N1max - N1min + 1L

  .require_cpp_o2simps_fn("cpp_o2simps_build_B_WGD_triplet")
  tri <- cpp_o2simps_build_B_WGD_triplet(
    as.integer(N0min),
    as.integer(N0max),
    as.integer(N1min),
    as.integer(N1max),
    boundary = boundary,
    wgd_value = 1.0
  )
  B <- sparseMatrix(
    i = as.integer(tri$i),
    j = as.integer(tri$j),
    x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
    repr = "C"
  )
  return(if (isTRUE(return_sparse)) B else as.matrix(B))
}

# -----------------------------------------------------------------------------
# Function: .build_G_with_WGD
# Purpose: Build generator matrix at the current oxygen/burden condition.
# Parameters:
#   - N0min: Minimum ploidy state on source grid.
#   - N0max: Maximum ploidy state on source grid.
#   - lambda0_vec: Function-specific input argument.
#   - p0_vec: Function-specific input argument.
#   - wgd_prob_vec: Function-specific input argument.
#   - mr_lethality0: Function-specific input argument.
#   - mr_buffer_by_ploidy: Function-specific input argument.
#   - N_unit: Ploidy scaling unit used to map integer states to N values.
#   - P_low: Function-specific input argument.
#   - P_high: Function-specific input argument.
#   - boundary: Boundary handling mode when transitions leave the ploidy grid.
#   - eps_tail: Small truncation threshold for tail probabilities.
#   - buffer_smax: Maximum per-copy survival factor.
#   - buffer_beta: Ploidy-buffering strength.
#   - buffer_n_exp: Ploidy-buffering exponent.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_G_with_WGD <- function(
    N0min, N0max, lambda0_vec, p0_vec, wgd_prob_vec,
    mr_lethality0 = 0.9, mr_lethality1 = 0.9,
    mr_buffer_by_ploidy = TRUE, N_unit = 22L, P_low = 2.0, P_high = 4.0,
    boundary = "drop", eps_tail = 1e-8,
    buffer_smax = 1.0, buffer_beta = 0.0,
    buffer_n_exp = 1.0
) {
  R0 <- N0max - N0min + 1L
  if (length(lambda0_vec) == 1L) lambda0_vec <- rep(lambda0_vec, R0)
  if (length(p0_vec) == 1L) p0_vec <- rep(p0_vec, R0)
  if (length(wgd_prob_vec) == 1L) wgd_prob_vec <- rep(wgd_prob_vec, R0)
  wgd_prob_vec <- .clip01(wgd_prob_vec)

  B0 <- .build_B_total(
    N0min, N0max, p_vec = p0_vec,
    boundary = boundary, eps_tail = eps_tail,
    N_unit = N_unit,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp
  )
  BW <- .build_B_WGD(N0min, N0max, N0min, N0max, boundary = boundary)
  L0 <- Diagonal(x = lambda0_vec)
  S0 <- Diagonal(x = (1 - wgd_prob_vec))
  SW <- Diagonal(x = wgd_prob_vec)
  ((B0 %*% S0) %*% L0) + ((BW %*% SW) %*% L0) - L0
}

# -----------------------------------------------------------------------------
# Function: run_all_sims
# Purpose: Run forward simulations for all configured scenarios and collect outputs.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
run_all_sims <- function(run_params) {
  all_results_list <- list()
  passage_times <- list()

  init_P_2N <- x$P[x$passage == 0 & x$ploidy == "2N"]
  init_P_4N <- x$P[x$passage == 0 & x$ploidy == "4N"]

  buffer_smax <- as.numeric(.first_non_null(run_params$buffer_smax, 1.0))
  buffer_beta <- as.numeric(.first_non_null(run_params$buffer_beta, 0.0))
  buffer_n_exp <- as.numeric(.first_non_null(run_params$buffer_n_exp, 1.0))
  boundary_mode <- as.character(.first_non_null(run_params$boundary, "drop"))
  O2_crit_use <- as.numeric(.first_non_null(run_params$O2_crit, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0

  .require_cpp_o2simps_fn("cpp_o2simps_build_G_for_o2_triplet")
  lam_max_use <- as.numeric(.first_non_null(run_params$lam_max, 0.0))
  mu_hp_use <- as.numeric(.first_non_null(run_params$mu_hp, 0.0))
  gamma_mu_use <- as.numeric(.first_non_null(run_params$gamma_mu, 1.0))
  n_O_use <- as.numeric(.first_non_null(run_params$n_O, 1.0))
  o2_Nref_use <- as.numeric(.first_non_null(run_params$o2_Nref, if (exists("cfg", inherits = TRUE)) get("cfg", inherits = TRUE)$o2_Nref else NULL, 1e6))
  cfg_o2_growth <- if (exists("cfg", inherits = TRUE)) get("cfg", inherits = TRUE)$O2_growth else NULL
  o2_growth_use <- isTRUE(.first_non_null(run_params$O2_growth, cfg_o2_growth, TRUE))
  ploidy_O2_death_mode_use <- assert_canonical_ploidy_o2_death_mode(
    .first_non_null(run_params$ploidy_O2_death, "diploid_NULL")
  )
  if (!is.finite(mu_hp_use) || mu_hp_use < 0) mu_hp_use <- 0.0
  if (!is.finite(gamma_mu_use) || gamma_mu_use <= 0) gamma_mu_use <- 1.0
  if (!is.finite(n_O_use) || n_O_use < 0) stop("run_params$n_O must be finite and >= 0.")
  if (!is.finite(o2_Nref_use) || o2_Nref_use <= 0) o2_Nref_use <- 1e6

  G_cache <- new.env(parent = emptyenv())
  build_G_for_resources <- function(O2) {
    O2_use <- .assert_o2_pct(as.numeric(O2), label = "O2")
    key <- sprintf("%.3f", O2_use)
    if (!exists(key, envir = G_cache, inherits = FALSE)) {
      tri <- cpp_o2simps_build_G_for_o2_triplet(
        O2 = as.numeric(O2_use),
        O2_crit = as.numeric(O2_crit_use),
        N0min = as.integer(N_MIN),
        N0max = as.integer(N_MAX),
        N1min = as.integer(N_MIN),
        N1max = as.integer(N_MAX),
        lam_max = as.numeric(lam_max_use),
        p_mis_base = as.numeric(.first_non_null(run_params$p_mis_base, 1e-5)),
        p_misseg = as.numeric(.first_non_null(run_params$p_misseg, 0.0)),
        k_o_mis = as.numeric(.first_non_null(run_params$k_o_mis, 50.0)),
        p_wgd = as.numeric(.first_non_null(run_params$p_wgd, 0.0)),
        boundary = as.character(boundary_mode),
        eps_tail = as.numeric(1e-8),
        buffer_smax = as.numeric(buffer_smax),
        buffer_beta = as.numeric(buffer_beta),
        buffer_n_exp = as.numeric(buffer_n_exp),
        N_unit = as.integer(N_UNIT),
        beta_size = 0.0,
        O2_growth = isTRUE(o2_growth_use),
        alpha_o2 = as.numeric(.first_non_null(run_params$alpha_o2, 0.0)),
        gamma_growth = as.numeric(.first_non_null(run_params$gamma_growth, 1.0)),
        mu_hp = as.numeric(mu_hp_use),
        gamma_mu = as.numeric(gamma_mu_use),
        n_O = as.numeric(n_O_use),
        ploidy_O2_death = as.character(ploidy_O2_death_mode_use)
      )
      G <- sparseMatrix(
        i = as.integer(tri$i),
        j = as.integer(tri$j),
        x = as.numeric(tri$x),
        dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
        repr = "C"
      )
      assign(key, G, envir = G_cache)
    }
    get(key, envir = G_cache, inherits = FALSE)
  }

  for (sim in sim_configs) {
    O2_LEVEL <- .assert_o2_pct(as.numeric(sim$O2), label = "sim$O2")

    init_P_values <- if (sim$init_ploidy == "2N") init_P_2N else init_P_4N
    x_current <- create_initial_dist(
      init_P_values,
      grid_pre,
      N_UNIT,
      chr_lengths_bp = default_chr_lengths_bp_1to22()
    )
    x_current <- x_current / sum(x_current)

    sim_passage_times <- numeric(PASSAGES_TO_RUN)

    for (p in 1:PASSAGES_TO_RUN) {
      pop_start <- sum(x_current)
      pop_target <- pop_start * POP_GROWTH_FACTOR
      time_in_passage <- 0.0

      while (sum(x_current) < pop_target) {
        x_prev <- as.numeric(x_current)
        mu_vec_step <- as.numeric(.mu_eff_of_O2(
          O2 = O2_LEVEL,
          run_params = run_params,
          N = grid_pre,
          O2_crit = O2_crit_use
        ))
        G <- build_G_for_resources(O2_LEVEL)
        x_div <- step_dt(G, x_prev, DT, 1L)
        x_next <- x_div - DT * mu_vec_step * x_prev
        x_next[!is.finite(x_next) | x_next < 0] <- 0
        x_current <- x_next
        time_in_passage <- time_in_passage + DT
        if (sum(x_current) < pop_start * 1e-3 || time_in_passage > 1000) {
          break
        }
      }
      sim_passage_times[p] <- time_in_passage

      if (p %in% REPORT_PASSAGES) {
        pop_total <- sum(x_current)
        dist_df <- data.frame(
          sim_id = sim$id,
          passage = p,
          layer = "single",
          N = grid_pre,
          fraction = x_current / pop_total
        )
        all_results_list[[length(all_results_list) + 1]] <- dist_df
      }
      x_current <- x_current / sum(x_current) * pop_start
    }

    passage_times[[sim$id]] <- data.frame(
      sim_id = sim$id,
      passage = 1:PASSAGES_TO_RUN,
      duration = sim_passage_times
    )
  }

  all_dists <- do.call(rbind, all_results_list)
  all_passage_times <- do.call(rbind, passage_times)

  list(all_dists = all_dists, all_passage_times = all_passage_times)
}

# -----------------------------------------------------------------------------
# Function: plot_misseg_interp
# Purpose: Plot missegregation response curve under current parameterization.
# Parameters:
#   - par: Function-specific input argument.
#   - o2_ref: Reference oxygen level used for plotting interpolation curves.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
plot_misseg_interp <- function(par, o2_ref = 20.5) {
  O2 <- seq(0, 100, length.out = 401)
  p <- .pmisseg_of_O2(
    O2 = O2,
    run_params = par,
    N = 44,
    O2_crit = as.numeric(.first_non_null(par$O2_crit, 1.0))
  )
  df <- data.frame(O2_pct = O2, p = as.numeric(p))
  ggplot(df, aes(O2_pct, p)) +
    geom_line(linewidth = 1, color = "black") +
    geom_point(data = df[c(1, nrow(df)), ], size = 2, color = "red") +
    labs(x = "Oxygen (%)", y = "Missegregation rate") +
    theme_bw()
}

# ---- embedded fixed-O2 simulation/attractor functions -------------------
`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
.first_non_null_local <- o2sd_first_non_null

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript fix_o2_simulation.R --fit_dir=/path/to/seed --simulation=invivo --initial_ploidy=3 --time_days=30 --n_sim=1 --o2=2",
      "",
      "Required:",
      "  --fit_dir, --run_dir, or --best_params",
      "                                   Seed dir, parent run dir containing seedXX dirs, or best_params.tsv.",
      "  --simulation                    invivo, invitro, or joint.",
      "  --initial_ploidy or --initial_ploidy_values",
      "                                   Initial population ploidy value(s), comma-separated for batch.",
      "  --time_days                     Simulation horizon in days.",
      "  --n_sim                         Number of repeated deterministic trajectories.",
      "  --o2 or --o2_values             Fixed effective O2 percentage(s), comma-separated for batch.",
      "",
      "Optional:",
      "  --out_dir                       Single-run output dir, or batch output root; defaults to repo oxygen/results/O2_fixed_simulation.",
      "  --seeds                         Seed filter, e.g. 1:500, 25,322, or seed25,seed322.",
      "  --n_core                        Parallel worker count for batch tasks; defaults to 1.",
      "  --build_task_list_only          TRUE/FALSE; write batch task_list.tsv without running tasks.",
      "  --dt                            Model integration step in days.",
      "  --save_every_days               Output interval in days; defaults to 1.",
      "  --report_dt                     Backward-compatible alias for --save_every_days.",
      "  --initial_cells                 Initial live cell count; defaults to 1.",
      "  --seed_id                       Optional seed label when --best_params has no seed directory.",
      "  --joint_scope                   shared_invivo or invitro_effective; defaults to shared_invivo.",
      "  --Crowding                      TRUE/FALSE; defaults to fit config or FALSE.",
      "",
      sep = "\n"
    )
  )
}

resolve_path_value <- function(path_value, base_dir = getwd()) {
  if (is.null(path_value) || !length(path_value)) return(NULL)
  txt <- trimws(as.character(path_value[[1]]))
  if (!nzchar(txt)) return(NULL)
  if (startsWith(txt, "~")) return(normalizePath(path.expand(txt), mustWork = FALSE))
  if (grepl("^(/|[A-Za-z]:[/\\\\])", txt)) return(normalizePath(txt, mustWork = FALSE))
  normalizePath(file.path(base_dir, txt), mustWork = FALSE)
}

normalize_simulation_mode <- function(x) {
  s <- tolower(gsub("[[:space:]_-]+", "", trimws(as.character(x[[1]] %||% ""))))
  if (s %in% c("invivo", "vivo")) return("invivo")
  if (s %in% c("invitro", "vitro", "invivtro")) return("invitro")
  if (s %in% c("joint", "jointfit", "jointfitting")) return("joint")
  stop("Invalid --simulation. Expected one of: invivo, invitro, joint.")
}

safe_id <- function(...) {
  x <- paste(..., sep = "_")
  x <- gsub("\\+", "", x)
  x <- gsub("-", "m", x)
  x <- gsub("\\.", "p", x)
  gsub("[^A-Za-z0-9_]+", "", x)
}

num_path_tag <- function(x) {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) != 1L || !is.finite(val)) return("NA")
  txt <- format(val, scientific = FALSE, trim = TRUE)
  txt <- sub("^-", "m", txt)
  txt <- gsub("\\.", "p", txt)
  txt <- gsub("[^A-Za-z0-9]+", "", txt)
  if (!nzchar(txt)) "NA" else txt
}

resolve_seed_id <- function(seed_arg = NULL, fit_dir = NULL, best_params_path = NULL) {
  if (!is.null(seed_arg) && length(seed_arg) > 0L) {
    txt <- trimws(as.character(seed_arg[[1]]))
    if (nzchar(txt)) return(safe_id(txt))
  }
  candidates <- character(0)
  if (!is.null(fit_dir) && nzchar(as.character(fit_dir))) {
    candidates <- c(candidates, basename(normalizePath(fit_dir, mustWork = FALSE)))
  }
  if (!is.null(best_params_path) && nzchar(as.character(best_params_path))) {
    candidates <- c(candidates, basename(dirname(normalizePath(best_params_path, mustWork = FALSE))))
  }
  hit <- candidates[grepl("^seed[0-9A-Za-z_.-]+$", candidates)]
  if (length(hit) > 0L) return(safe_id(hit[[1]]))
  "seed_unknown"
}

default_output_dir <- function(simulation, fixed_o2, initial_ploidy, seed_id) {
  root <- normalizePath(
    file.path(REPO_ROOT, "oxygen", "results", "O2_fixed_simulation"),
    mustWork = FALSE
  )
  o2_tag <- num_path_tag(fixed_o2)
  ploidy_tag <- num_path_tag(initial_ploidy)
  file.path(
    root,
    safe_id(simulation),
    paste0("O2_", o2_tag),
    paste0("ploidy", ploidy_tag),
    seed_id
  )
}

resolve_save_every_days <- function(argv) {
  save_every <- as_num(argv$save_every_days, NA_real_)
  if (is.finite(save_every) && save_every > 0) return(save_every)
  legacy_report_dt <- as_num(argv$report_dt, NA_real_)
  if (is.finite(legacy_report_dt) && legacy_report_dt > 0) return(legacy_report_dt)
  1.0
}

resolve_simulation_ids <- function(argv, n_sim) {
  if (!is.null(argv$simulation_id) && length(argv$simulation_id) > 0L) {
    raw_ids <- split_cli_values(argv$simulation_id)
    ids <- suppressWarnings(as.integer(raw_ids))
    if (!length(ids) || any(is.na(ids)) || any(ids < 1L)) {
      stop("--simulation_id must contain positive integer value(s).")
    }
  } else {
    ids <- seq_len(n_sim)
  }
  if (length(ids) != n_sim) {
    stop("--simulation_id count must match --n_sim. For one worker task use --n_sim=1 --simulation_id=<id>.")
  }
  ids
}

split_cli_values <- function(x) {
  if (is.null(x) || !length(x)) return(character(0))
  txt <- trimws(as.character(x[[1]]))
  if (!nzchar(txt)) return(character(0))
  vals <- trimws(unlist(strsplit(txt, "[,;[:space:]]+", perl = TRUE)))
  vals[nzchar(vals)]
}

parse_o2_values <- function(argv) {
  raw <- split_cli_values(.first_non_null_local(argv$o2_values, argv$o2))
  vals <- suppressWarnings(as.numeric(raw))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) stop("Provide --o2 or --o2_values with at least one finite O2 value.")
  bad <- vals < 0 | vals > 100
  if (any(bad)) stop("O2 values must be between 0 and 100: ", paste(vals[bad], collapse = ","))
  vals
}

parse_initial_ploidy_values <- function(argv) {
  raw <- split_cli_values(.first_non_null_local(argv$initial_ploidy_values, argv$initial_ploidy))
  vals <- suppressWarnings(as.numeric(raw))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    stop("Provide --initial_ploidy or --initial_ploidy_values with at least one finite ploidy value.")
  }
  bad <- vals <= 0
  if (any(bad)) stop("Initial ploidy values must be > 0: ", paste(vals[bad], collapse = ","))
  vals
}

normalize_seed_label <- function(x) {
  s <- trimws(as.character(x))
  s <- sub("^seed", "", s, ignore.case = TRUE)
  if (!nzchar(s)) return(NA_character_)
  safe_id(paste0("seed", s))
}

parse_seed_values <- function(x) {
  raw <- split_cli_values(x)
  expanded <- unlist(lapply(raw, function(token) {
    token <- trimws(as.character(token))
    range_match <- regexec("^(?:seed)?([0-9]+)\\s*:\\s*(?:seed)?([0-9]+)$", token, ignore.case = TRUE)
    hit <- regmatches(token, range_match)[[1]]
    if (length(hit) == 3L) {
      start <- as.integer(hit[[2]])
      end <- as.integer(hit[[3]])
      return(as.character(seq.int(start, end)))
    }
    token
  }), use.names = FALSE)
  out <- vapply(expanded, normalize_seed_label, character(1))
  out <- out[!is.na(out) & nzchar(out)]
  unique(out)
}

seed_number_for_order <- function(seed_id) {
  n <- suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
  ifelse(is.na(n), .Machine$integer.max, n)
}

resolve_batch_run_dir <- function(argv, fit_dir = NULL) {
  run_dir <- resolve_path_value(argv$run_dir, getwd())
  if (!is.null(run_dir)) return(run_dir)
  if (!is.null(fit_dir) && dir.exists(fit_dir) && !file.exists(file.path(fit_dir, "best_params.tsv"))) {
    seed_dirs <- list.dirs(fit_dir, recursive = FALSE, full.names = TRUE)
    if (any(file.exists(file.path(seed_dirs, "best_params.tsv")))) {
      return(fit_dir)
    }
  }
  NULL
}

discover_seed_dirs <- function(run_dir, seed_ids = character(0)) {
  if (is.null(run_dir) || !dir.exists(run_dir)) {
    stop("--run_dir does not exist: ", run_dir)
  }
  if (length(seed_ids) > 0L) {
    dirs <- file.path(run_dir, seed_ids)
  } else {
    dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[grepl("^seed[0-9A-Za-z_.-]+$", basename(dirs))]
  }
  ids <- basename(dirs)
  keep <- dir.exists(dirs) & file.exists(file.path(dirs, "best_params.tsv"))
  missing <- ids[!keep]
  if (length(missing) > 0L) {
    warning("Skipping seed dirs without best_params.tsv: ", paste(missing, collapse = ","))
  }
  dirs <- normalizePath(dirs[keep], mustWork = TRUE)
  ids <- safe_id(ids[keep])
  if (!length(dirs)) stop("No eligible seed directories found in --run_dir: ", run_dir)
  ord <- order(seed_number_for_order(ids), ids)
  data.frame(seed_id = ids[ord], seed_dir = dirs[ord], stringsAsFactors = FALSE)
}

batch_output_root <- function(argv) {
  out_dir <- resolve_path_value(argv$out_dir, getwd())
  if (!is.null(out_dir)) return(out_dir)
  normalizePath(file.path(REPO_ROOT, "oxygen", "results", "O2_fixed_simulation"), mustWork = FALSE)
}

batch_leaf_output_dir <- function(root, simulation, fixed_o2, initial_ploidy, seed_id) {
  o2_tag <- num_path_tag(fixed_o2)
  ploidy_tag <- num_path_tag(initial_ploidy)
  file.path(
    root,
    safe_id(simulation),
    paste0("O2_", o2_tag),
    paste0("ploidy", ploidy_tag),
    seed_id
  )
}

build_task_list <- function(argv, simulation, run_dir, seed_tab, o2_values, initial_ploidy_values, n_sim) {
  root <- batch_output_root(argv)
  rows <- list()
  task_id <- 0L
  for (i in seq_len(nrow(seed_tab))) {
    for (o2 in o2_values) {
      for (initial_ploidy in initial_ploidy_values) {
        leaf_dir <- batch_leaf_output_dir(
          root = root,
          simulation = simulation,
          fixed_o2 = o2,
          initial_ploidy = initial_ploidy,
          seed_id = seed_tab$seed_id[[i]]
        )
        for (sim_id in seq_len(n_sim)) {
          task_id <- task_id + 1L
          rows[[task_id]] <- data.frame(
            task_id = task_id,
            seed_id = seed_tab$seed_id[[i]],
            seed_dir = seed_tab$seed_dir[[i]],
            fixed_o2_pct = as.numeric(o2),
            initial_ploidy = as.numeric(initial_ploidy),
            simulation_id = as.integer(sim_id),
            result_dir = leaf_dir,
            output_dir = file.path(leaf_dir, paste0("sim", sim_id)),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  tasks <- do.call(rbind, rows)
  attr(tasks, "output_root") <- root
  attr(tasks, "run_dir") <- run_dir
  tasks
}

optional_forward_args <- function(argv) {
  keys <- c(
    "dt", "save_every_days", "report_dt", "initial_cells",
    "joint_scope", "Crowding", "O2_growth", "start_with", "ploidy_O2_death"
  )
  out <- character(0)
  for (key in keys) {
    if (!is.null(argv[[key]]) && length(argv[[key]]) > 0L) {
      out <- c(out, paste0("--", key, "=", as.character(argv[[key]][[1]])))
    }
  }
  out
}

task_script_args <- function(task, argv, simulation, time_days) {
  c(
    normalizePath(file.path(SCRIPT_DIR, "fix_o2_simulation.R"), mustWork = TRUE),
    "--worker_task=TRUE",
    paste0("--fit_dir=", task$seed_dir),
    paste0("--simulation=", simulation),
    paste0("--initial_ploidy=", task$initial_ploidy),
    paste0("--time_days=", time_days),
    "--n_sim=1",
    paste0("--simulation_id=", task$simulation_id),
    paste0("--o2=", task$fixed_o2_pct),
    paste0("--seed_id=", task$seed_id),
    paste0("--out_dir=", task$output_dir),
    optional_forward_args(argv)
  )
}

run_task_process <- function(task, argv, simulation, time_days) {
  dir.create(task$output_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(task$output_dir, "task.log")
  args <- task_script_args(task, argv, simulation, time_days)
  status <- system2(file.path(R.home("bin"), "Rscript"), args = args, stdout = log_path, stderr = log_path)
  status <- if (is.null(status)) 0L else as.integer(status)
  data.frame(
    task_id = task$task_id,
    seed_id = task$seed_id,
    fixed_o2_pct = task$fixed_o2_pct,
    initial_ploidy = task$initial_ploidy,
    simulation_id = task$simulation_id,
    status = status,
    output_dir = task$output_dir,
    log_file = log_path,
    stringsAsFactors = FALSE
  )
}

fixo2_simulation_output_paths <- function(out_dir) {
  list(
    run_config = file.path(out_dir, "run_config.tsv"),
    parameters_used = file.path(out_dir, "parameters_used.tsv"),
    population = file.path(out_dir, "population_trajectory.tsv"),
    rate = file.path(out_dir, "rate_trajectory.tsv.gz"),
    state = file.path(out_dir, "state_trajectory.tsv.gz")
  )
}

fixo2_existing_file_ok <- function(path) {
  file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
}

fixo2_simulation_output_complete <- function(out_dir,
                                             require_rate = TRUE,
                                             require_state = TRUE) {
  paths <- fixo2_simulation_output_paths(out_dir)
  required <- c(paths$run_config, paths$parameters_used, paths$population)
  if (isTRUE(require_rate)) required <- c(required, paths$rate)
  if (isTRUE(require_state)) required <- c(required, paths$state)
  all(vapply(required, fixo2_existing_file_ok, logical(1)))
}

fixo2_annotate_task_completion <- function(tasks) {
  if (!nrow(tasks)) {
    tasks$complete <- logical(0)
    return(tasks)
  }
  tasks$complete <- vapply(tasks$output_dir, fixo2_simulation_output_complete, logical(1))
  tasks
}

fixo2_filter_missing_simulation_tasks <- function(tasks, force = FALSE) {
  tasks <- fixo2_annotate_task_completion(tasks)
  if (isTRUE(force)) {
    tasks$skip_reason <- ifelse(tasks$complete, "complete_but_force_requested", "")
    return(tasks)
  }
  tasks$skip_reason <- ifelse(tasks$complete, "complete", "")
  tasks[!tasks$complete, , drop = FALSE]
}

run_batch <- function(argv, simulation, fit_dir, best_params_path, o2_values, initial_ploidy_values, time_days, n_sim) {
  if (!is.null(best_params_path)) {
    stop("Batch mode requires --run_dir or parent --fit_dir, not --best_params.")
  }
  run_dir <- resolve_batch_run_dir(argv, fit_dir = fit_dir)
  if (is.null(run_dir)) {
    stop("Batch mode requires --run_dir or a parent --fit_dir containing seedXX subdirectories.")
  }
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  seed_ids <- parse_seed_values(argv$seeds)
  seed_tab <- discover_seed_dirs(run_dir, seed_ids = seed_ids)
  n_core <- as_int(argv$n_core, 1L)
  if (!is.finite(n_core) || is.na(n_core) || n_core < 1L) n_core <- 1L
  tasks <- build_task_list(
    argv = argv,
    simulation = simulation,
    run_dir = run_dir,
    seed_tab = seed_tab,
    o2_values = o2_values,
    initial_ploidy_values = initial_ploidy_values,
    n_sim = n_sim
  )
  n_core <- min(as.integer(n_core), nrow(tasks))
  task_root <- file.path(attr(tasks, "output_root"), safe_id(simulation))
  dir.create(task_root, recursive = TRUE, showWarnings = FALSE)
  force <- as_bool(argv$force, FALSE)
  tasks_annotated <- fixo2_annotate_task_completion(tasks)
  pending_tasks <- fixo2_filter_missing_simulation_tasks(tasks, force = force)
  task_list_path <- file.path(task_root, "task_list.tsv")
  pending_task_list_path <- file.path(task_root, "pending_task_list.tsv")
  write_tsv(tasks_annotated, task_list_path)
  write_tsv(pending_tasks, pending_task_list_path)

  message("Fixed-O2 batch simulation")
  message("  run_dir: ", run_dir)
  message("  seeds: ", paste(seed_tab$seed_id, collapse = ","))
  message("  o2_values: ", paste(o2_values, collapse = ","))
  message("  initial_ploidy_values: ", paste(initial_ploidy_values, collapse = ","))
  message("  n_sim: ", n_sim)
  message("  task_count: ", nrow(tasks), " total; ", nrow(pending_tasks), " pending")
  message("  n_core: ", n_core)
  message("  task_list: ", task_list_path)
  message("  pending_task_list: ", pending_task_list_path)
  if (as_bool(argv$build_task_list_only, FALSE)) {
    message("  build_task_list_only: TRUE")
    return(invisible(pending_tasks))
  }
  if (!nrow(pending_tasks)) {
    message("  all requested fixed-O2 simulation outputs are already complete; no tasks were run.")
    status <- tasks_annotated
    status$status <- 0L
    status$log_file <- NA_character_
    status_path <- file.path(task_root, "task_status.tsv")
    write_tsv(status, status_path)
    return(invisible(status))
  }

  n_core <- min(as.integer(n_core), nrow(pending_tasks))
  task_rows <- split(pending_tasks, seq_len(nrow(pending_tasks)))
  runner <- function(task_df) {
    run_task_process(task_df[1, , drop = FALSE], argv, simulation, time_days)
  }
  if (n_core > 1L) {
    if (.Platform$OS.type == "windows") {
      warning("Parallel batch execution with n_core > 1 is not supported on Windows in this script; falling back to sequential.")
      status_list <- lapply(task_rows, runner)
    } else {
      status_list <- parallel::mclapply(task_rows, runner, mc.cores = n_core, mc.preschedule = FALSE)
    }
  } else {
    status_list <- lapply(task_rows, runner)
  }
  status <- do.call(rbind, status_list)
  status_path <- file.path(task_root, "task_status.tsv")
  write_tsv(status, status_path)
  message("  task_status: ", status_path)
  if (any(status$status != 0L)) {
    failed <- status[status$status != 0L, , drop = FALSE]
    stop(
      "Batch failed for task(s): ",
      paste(failed$task_id, collapse = ","),
      ". Check task_status.tsv and task logs."
    )
  }
  message("Batch done.")
  invisible(status)
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required file was not found: ", path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

read_best_params_table <- function(path) {
  tab <- read_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain columns: parameter, value")
  }
  params <- trimws(as.character(tab$parameter))
  vals <- suppressWarnings(as.numeric(tab$value))
  bad <- !nzchar(params) | !is.finite(vals)
  if (any(bad)) {
    stop("Non-finite or unnamed parameter values in ", path, ": ", paste(params[bad], collapse = ", "))
  }
  if (any(duplicated(params))) {
    stop("Duplicate parameters in ", path, ": ", paste(unique(params[duplicated(params)]), collapse = ", "))
  }
  stats::setNames(vals, params)
}

read_joint_scope_params <- function(path, scope) {
  tab <- read_tsv(path)
  if (!all(c("scope", "parameter", "value") %in% names(tab))) {
    stop("joint_best_params_long.tsv must contain columns: scope, parameter, value")
  }
  scope <- trimws(as.character(scope[[1]]))
  sub <- tab[trimws(as.character(tab$scope)) == scope, , drop = FALSE]
  if (!nrow(sub)) {
    stop("No rows found for joint scope '", scope, "' in ", path)
  }
  params <- trimws(as.character(sub$parameter))
  vals <- suppressWarnings(as.numeric(sub$value))
  bad <- !nzchar(params) | !is.finite(vals)
  if (any(bad)) {
    stop("Non-finite or unnamed joint parameters in ", path, ": ", paste(params[bad], collapse = ", "))
  }
  if (any(duplicated(params))) {
    stop("Duplicate joint parameters in ", path, " for scope ", scope, ": ", paste(unique(params[duplicated(params)]), collapse = ", "))
  }
  stats::setNames(vals, params)
}

resolve_parameter_values <- function(fit_dir, best_params_path, simulation, joint_scope) {
  explicit_best <- !is.null(best_params_path)
  if (explicit_best) {
    path <- normalizePath(best_params_path, mustWork = TRUE)
    return(list(values = read_best_params_table(path), path = path, source = "best_params"))
  }

  if (is.null(fit_dir)) {
    stop("Provide --fit_dir or --best_params.")
  }
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  joint_long <- file.path(fit_dir, "joint_best_params_long.tsv")
  if (file.exists(joint_long) && simulation %in% c("joint", "invitro")) {
    scope <- if (identical(simulation, "invitro")) {
      as.character(.first_non_null_local(joint_scope, "invitro_effective"))
    } else {
      as.character(.first_non_null_local(joint_scope, "shared_invivo"))
    }
    path <- normalizePath(joint_long, mustWork = TRUE)
    return(list(values = read_joint_scope_params(path, scope), path = path, source = paste0("joint_best_params_long:", scope)))
  }

  path <- file.path(fit_dir, "best_params.tsv")
  if (!file.exists(path)) stop("Cannot find best_params.tsv in --fit_dir: ", fit_dir)
  path <- normalizePath(path, mustWork = TRUE)
  list(values = read_best_params_table(path), path = path, source = "best_params")
}

read_fit_config <- function(fit_dir) {
  if (is.null(fit_dir)) return(list())
  path <- file.path(fit_dir, "fit_config.rds")
  if (!file.exists(path)) return(list())
  cfg <- readRDS(path)
  if (!is.list(cfg)) stop("fit_config.rds did not contain a list: ", path)
  cfg
}

validate_range <- function(rp, name, lower = -Inf, upper = Inf, strict_lower = FALSE, strict_upper = FALSE) {
  val <- suppressWarnings(as.numeric(rp[[name]]))
  ok <- length(val) == 1L && is.finite(val)
  if (ok && strict_lower) ok <- val > lower else if (ok) ok <- val >= lower
  if (ok && strict_upper) ok <- val < upper else if (ok) ok <- val <= upper
  if (!ok) {
    lo <- if (is.finite(lower)) paste0(if (strict_lower) ">" else ">=", lower) else ""
    hi <- if (is.finite(upper)) paste0(if (strict_upper) "<" else "<=", upper) else ""
    stop("Invalid parameter ", name, ": expected ", paste(c(lo, hi)[nzchar(c(lo, hi))], collapse = " and "), ".")
  }
  invisible(TRUE)
}

prepare_run_params <- function(param_values, simulation, cfg, fixed_o2) {
  core_required <- c(
    "lam_max", "p_mis_base", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp",
    "p_wgd", "alpha_o2", "gamma_growth", "mu_hp",
    "gamma_mu", "O2_crit", "n_O"
  )
  required <- switch(
    simulation,
    invivo = c(core_required, "rho_2N"),
    invitro = core_required,
    joint = c(core_required, "rho_2N")
  )
  missing <- setdiff(required, names(param_values))
  if (length(missing) > 0L) {
    stop("Missing required ", simulation, " parameter(s): ", paste(missing, collapse = ", "))
  }

  rp <- as.list(param_values)
  rp$boundary <- as.character(.first_non_null_local(rp$boundary, cfg$boundary, "drop"))
  rp$o2_S0 <- as.numeric(fixed_o2)
  rp$o2_min <- as.numeric(.first_non_null_local(rp$o2_min, cfg$o2_min, 0.0))
  rp$o2_Nref <- as.numeric(.first_non_null_local(rp$o2_Nref, cfg$o2_Nref, cfg$init_total_size, 1e6))
  rp$o2_S0_upper_bound <- max(
    as.numeric(.first_non_null_local(rp$o2_S0_upper_bound, cfg$o2_S0_upper_bound, fixed_o2, 5.0)),
    as.numeric(fixed_o2)
  )
  rp$kappa_O <- as.numeric(.first_non_null_local(rp$kappa_O, cfg$kappa_O_init, 0.0))
  rp$eta_o2 <- as.numeric(.first_non_null_local(rp$eta_o2, cfg$eta_o2_init, 1.0))
  rp$tau_O2 <- as.numeric(.first_non_null_local(rp$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 0.1))
  rp$rho_2N <- as.numeric(.first_non_null_local(rp$rho_2N, default_rho_2N_prior_center(cfg)))
  rp$k_clear <- as.numeric(.first_non_null_local(rp$k_clear, cfg$k_clear_init, 1e-4))
  rp$alpha <- as.numeric(.first_non_null_local(rp$alpha, 0.0))
  rp$gamma <- as.numeric(.first_non_null_local(rp$gamma, 1.0))
  rp$O2_growth <- isTRUE(.first_non_null_local(rp$O2_growth, cfg$O2_growth, TRUE))
  rp$ploidy_O2_death <- assert_canonical_ploidy_o2_death_mode(
    .first_non_null_local(rp$ploidy_O2_death, cfg$ploidy_O2_death, "ploidy_related")
  )

  validate_range(rp, "lam_max", lower = 0)
  validate_range(rp, "p_mis_base", lower = 0, upper = 1)
  validate_range(rp, "p_misseg", lower = 0, upper = 1)
  validate_range(rp, "p_wgd", lower = 0, upper = 1)
  validate_range(rp, "k_o_mis", lower = 0, strict_lower = TRUE)
  validate_range(rp, "buffer_smax", lower = 0, upper = 1 + 1e-8)
  validate_range(rp, "buffer_beta", lower = 0)
  validate_range(rp, "buffer_n_exp", lower = 0, strict_lower = TRUE)
  validate_range(rp, "alpha_o2", lower = 0)
  validate_range(rp, "gamma_growth", lower = 0, strict_lower = TRUE)
  validate_range(rp, "mu_hp", lower = 0)
  validate_range(rp, "gamma_mu", lower = 0, strict_lower = TRUE)
  validate_range(rp, "O2_crit", lower = 0)
  validate_range(rp, "n_O", lower = 0)
  validate_range(rp, "rho_2N", lower = 0, strict_lower = TRUE)
  validate_range(rp, "o2_S0", lower = 0, upper = 100)
  validate_range(rp, "o2_S0_upper_bound", lower = 0, upper = 100, strict_lower = TRUE)
  validate_range(rp, "o2_Nref", lower = 0, strict_lower = TRUE)
  validate_range(rp, "k_clear", lower = 0)
  validate_range(rp, "tau_O2", lower = 0, strict_lower = TRUE)
  validate_range(rp, "eta_o2", lower = 0)

  rp
}

prepare_sim_cfg <- function(cfg, argv, fixed_o2, run_params) {
  cfg$N_UNIT <- as_int(cfg$N_UNIT, 22L)
  cfg$N_MIN <- as_int(cfg$N_MIN, 22L)
  cfg$N_MAX <- as_int(cfg$N_MAX, 154L)
  cfg$DT <- as_num(argv$dt, as.numeric(.first_non_null_local(cfg$DT, 0.05)))
  if (!is.finite(cfg$DT) || cfg$DT <= 0) stop("dt must be finite and > 0.")
  cfg$chr_lengths_bp <- .first_non_null_local(cfg$chr_lengths_bp, default_chr_lengths_bp_1to22())
  cfg$start_with <- assert_canonical_start_with_mode(
    canonical_start_with_mode(.first_non_null_local(argv$start_with, cfg$start_with, "chr_number"))
  )
  cfg$ploidy_O2_death <- assert_canonical_ploidy_o2_death_mode(
    canonical_ploidy_o2_death_mode(.first_non_null_local(argv$ploidy_O2_death, cfg$ploidy_O2_death, run_params$ploidy_O2_death, "ploidy_related"))
  )
  cfg$o2_burden_feedback <- FALSE
  cfg$O2_growth <- as_bool(argv$O2_growth, isTRUE(.first_non_null_local(cfg$O2_growth, TRUE)))
  cfg$o2_S0_upper_bound <- max(
    as.numeric(.first_non_null_local(run_params$o2_S0_upper_bound, cfg$o2_S0_upper_bound, fixed_o2, 5.0)),
    as.numeric(fixed_o2)
  )
  cfg$o2_Nref <- as.numeric(.first_non_null_local(run_params$o2_Nref, cfg$o2_Nref, cfg$init_total_size, 1e6))
  cfg$o2_min <- as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.0))
  cfg$Crowding <- as_bool(argv$Crowding, isTRUE(.first_non_null_local(cfg$Crowding, FALSE)))
  cfg$K <- as.numeric(.first_non_null_local(cfg$K, 1e12))
  cfg$crowding <- as.character(.first_non_null_local(cfg$crowding, "logistic"))
  cfg$min_pop <- as.numeric(.first_non_null_local(cfg$min_pop, 1e-12))
  cfg$dose_ref <- as.numeric(.first_non_null_local(cfg$dose_ref, 30))
  cfg$tx_mult_min <- as.numeric(.first_non_null_local(cfg$tx_mult_min, 0.05))
  cfg$o2_cache_bin_pct <- as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01))
  cfg$o2_cache_hysteresis_pct <- as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005))
  cfg$o2_cache_profile <- FALSE
  cfg$burden_log_eps <- as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12))
  if (!is.finite(cfg$o2_S0_upper_bound) || cfg$o2_S0_upper_bound <= 0) stop("o2_S0_upper_bound must be > 0.")
  cfg
}

make_initial_state <- function(cfg, initial_ploidy, initial_cells) {
  initial_ploidy <- as.numeric(initial_ploidy)
  initial_cells <- as.numeric(initial_cells)
  if (!is.finite(initial_ploidy) || initial_ploidy <= 0) stop("initial_ploidy must be finite and > 0.")
  if (!is.finite(initial_cells) || initial_cells <= 0) stop("initial_cells must be finite and > 0.")
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  ploidy_grid <- as.numeric(grid_pre) / as.numeric(cfg$N_UNIT)
  min_ploidy <- min(ploidy_grid)
  max_ploidy <- max(ploidy_grid)
  if (initial_ploidy < min_ploidy || initial_ploidy > max_ploidy) {
    stop(
      "initial_ploidy=", initial_ploidy, " is outside model grid [",
      signif(min_ploidy, 6), ", ", signif(max_ploidy, 6), "]."
    )
  }
  idx <- which.min(abs(ploidy_grid - initial_ploidy))
  state <- numeric(length(grid_pre))
  state[[idx]] <- initial_cells
  list(
    state = state,
    requested_ploidy = initial_ploidy,
    used_N = as.integer(grid_pre[[idx]]),
    used_ploidy = as.numeric(ploidy_grid[[idx]])
  )
}

fixo2_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22L), as.integer(cfg$N_MAX %||% 154L))
  G <- o2pr_build_G(model_env, cfg, run_params, O2)
  mu_all <- as.numeric(o2ipa_call_model(
    model_env,
    ".mu_eff_of_O2",
    O2 = rep(O2, length(ngrid)),
    run_params = run_params,
    N = ngrid
  ))
  M <- as.matrix(G - Matrix::Diagonal(x = mu_all))
  list(M = M, G = G, mu_all = mu_all, ngrid = ngrid)
}

fixo2_init_vector <- function(ngrid, init_N, n_unit = 22) {
  idx <- which.min(abs(ngrid - init_N))
  v <- numeric(length(ngrid))
  v[[idx]] <- 1
  list(vector = v, used_N = ngrid[[idx]], used_ploidy = ngrid[[idx]] / n_unit)
}

fixo2_normalize_state_matrix <- function(state_mat) {
  state_mat <- Re(state_mat)
  state_mat[!is.finite(state_mat)] <- NA_real_
  state_mat <- pmax(state_mat, 0)
  sums <- colSums(state_mat, na.rm = TRUE)
  valid <- is.finite(sums) & sums > 0
  if (any(valid)) {
    state_mat[, valid] <- sweep(state_mat[, valid, drop = FALSE], 2L, sums[valid], "/")
  }
  if (any(!valid)) state_mat[, !valid] <- NA_real_
  state_mat
}

fixo2_normalize_state <- function(x) {
  fixo2_normalize_state_matrix(matrix(x, ncol = 1L))[, 1L]
}

fixo2_trajectory_from_state_matrix <- function(state_mat, ngrid, time_grid, n_unit) {
  mean_N <- as.numeric(crossprod(ngrid, state_mat))
  data.frame(
    day = time_grid,
    mean_N = mean_N,
    mean_ploidy = mean_N / n_unit,
    fraction_N_le_25 = colSums(state_mat[ngrid <= 25, , drop = FALSE], na.rm = TRUE),
    fraction_N_below_44 = colSums(state_mat[ngrid < 44, , drop = FALSE], na.rm = TRUE),
    fraction_N_ge_44 = colSums(state_mat[ngrid >= 44, , drop = FALSE], na.rm = TRUE),
    fraction_N_ge_66 = colSums(state_mat[ngrid >= 66, , drop = FALSE], na.rm = TRUE),
    fraction_N_ge_88 = colSums(state_mat[ngrid >= 88, , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

fixo2_eigen_trajectory_cached <- function(eig, ngrid, init, time_grid, n_unit) {
  coef <- tryCatch(solve(eig$vectors, init), error = function(e) NULL)
  if (is.null(coef)) return(list(status = "eigen_solve_failed", trajectory = data.frame()))
  lambda_ref <- max(Re(eig$values), na.rm = TRUE)
  weight_mat <- exp(outer(eig$values - lambda_ref, time_grid, `*`)) *
    matrix(coef, nrow = length(coef), ncol = length(time_grid))
  state_mat <- fixo2_normalize_state_matrix(eig$vectors %*% weight_mat)
  out <- fixo2_trajectory_from_state_matrix(state_mat, ngrid, time_grid, n_unit)
  status <- if (any(!is.finite(out$mean_ploidy))) "nonfinite_state" else "ok"
  list(status = status, trajectory = out)
}

fixo2_eigen_states <- function(M, init, time_grid) {
  eig <- tryCatch(eigen(as.matrix(M), only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) stop("eigen decomposition failed")
  coef <- tryCatch(solve(eig$vectors, init), error = function(e) NULL)
  if (is.null(coef)) stop("eigen coefficient solve failed")
  lambda_ref <- max(Re(eig$values), na.rm = TRUE)
  lapply(time_grid, function(tt) {
    weights <- exp((eig$values - lambda_ref) * tt) * coef
    fixo2_normalize_state(eig$vectors %*% weights)
  })
}

fixo2_eigen_trajectory <- function(M, ngrid, init, time_grid, n_unit) {
  eig <- tryCatch(eigen(M, only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) return(list(status = "eigen_failed", trajectory = data.frame()))
  fixo2_eigen_trajectory_cached(eig, ngrid, init, time_grid, n_unit)
}

fixo2_expm_states <- function(M, init, time_grid) {
  time_grid <- sort(unique(as.numeric(time_grid)))
  x <- as.numeric(init)
  t_now <- 0
  states <- vector("list", length(time_grid))
  expm_cache <- new.env(parent = emptyenv())
  get_step_expm <- function(delta) {
    key <- format(signif(delta, 15), scientific = FALSE, trim = TRUE)
    if (!exists(key, envir = expm_cache, inherits = FALSE)) {
      assign(key, Matrix::expm(M * delta), envir = expm_cache)
    }
    get(key, envir = expm_cache, inherits = FALSE)
  }
  for (i in seq_along(time_grid)) {
    target <- time_grid[[i]]
    delta <- target - t_now
    if (delta > 1e-12) {
      x <- as.numeric(get_step_expm(delta) %*% x)
      scale <- suppressWarnings(max(abs(x), na.rm = TRUE))
      if (!is.finite(scale) || scale <= 0) {
        x[] <- NA_real_
      } else {
        x <- x / scale
      }
      t_now <- target
    }
    states[[i]] <- fixo2_normalize_state(x)
  }
  states
}

fixo2_expm_trajectory <- function(M, ngrid, init, time_grid, n_unit) {
  time_grid <- sort(unique(as.numeric(time_grid)))
  states <- fixo2_expm_states(M, init, time_grid)
  rows <- lapply(seq_along(time_grid), function(i) {
    state <- states[[i]]
    target <- time_grid[[i]]
    mean_N <- sum(ngrid * state, na.rm = TRUE)
    data.frame(
      day = target,
      mean_N = mean_N,
      mean_ploidy = mean_N / n_unit,
      fraction_N_le_25 = sum(state[ngrid <= 25], na.rm = TRUE),
      fraction_N_below_44 = sum(state[ngrid < 44], na.rm = TRUE),
      fraction_N_ge_44 = sum(state[ngrid >= 44], na.rm = TRUE),
      fraction_N_ge_66 = sum(state[ngrid >= 66], na.rm = TRUE),
      fraction_N_ge_88 = sum(state[ngrid >= 88], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  status <- if (any(!is.finite(out$mean_ploidy))) "nonfinite_state" else "ok"
  list(status = status, trajectory = out)
}

fixo2_trajectory_with_fallback <- function(M, eig, ngrid, init, time_grid, n_unit) {
  sim <- fixo2_eigen_trajectory_cached(eig, ngrid, init, time_grid, n_unit)
  method <- "eigen_cached"
  if (!identical(sim$status, "ok")) {
    fallback <- tryCatch(fixo2_expm_trajectory(M, ngrid, init, time_grid, n_unit), error = function(e) NULL)
    if (!is.null(fallback) && identical(fallback$status, "ok")) {
      sim <- fallback
      method <- "expm_fallback"
    } else if (!is.null(fallback)) {
      sim <- fallback
      method <- "expm_fallback_failed"
    }
  }
  sim$trajectory_method <- method
  sim
}

fixo2_dominant_from_eig <- function(eig, ngrid, n_unit) {
  idx <- which.max(Re(eig$values))
  v <- Re(eig$vectors[, idx])
  if (sum(v, na.rm = TRUE) < 0) v <- -v
  v <- fixo2_normalize_state(v)
  lambda1 <- Re(eig$values[[idx]])
  lambda2 <- sort(Re(eig$values), decreasing = TRUE)[min(2L, length(eig$values))]
  dominant_mean_N <- sum(ngrid * v, na.rm = TRUE)
  spectral_gap <- lambda1 - lambda2
  data.frame(
    dominant_mean_N = dominant_mean_N,
    dominant_mean_ploidy = dominant_mean_N / n_unit,
    dominant_fraction_N_le_25 = sum(v[ngrid <= 25], na.rm = TRUE),
    dominant_fraction_N_below_44 = sum(v[ngrid < 44], na.rm = TRUE),
    dominant_fraction_N_ge_44 = sum(v[ngrid >= 44], na.rm = TRUE),
    dominant_growth_rate = lambda1,
    second_growth_rate = lambda2,
    spectral_gap = spectral_gap,
    relative_spectral_gap = spectral_gap / pmax(abs(lambda1), .Machine$double.eps),
    relax_time_days = ifelse(spectral_gap > 0, 1 / spectral_gap, NA_real_),
    time_to_10x_days = ifelse(spectral_gap > 0, log(10) / spectral_gap, NA_real_),
    time_to_100x_days = ifelse(spectral_gap > 0, log(100) / spectral_gap, NA_real_),
    log10_advantage_1000d = spectral_gap * 1000 / log(10),
    dominance_class = ifelse(
      !is.finite(spectral_gap) | spectral_gap <= 0, "non_positive",
      ifelse(
        spectral_gap < 0.001, "very_weak",
        ifelse(spectral_gap < 0.005, "weak", ifelse(spectral_gap < 0.01, "moderate", "strong"))
      )
    ),
    stringsAsFactors = FALSE
  )
}

fixo2_dominant_one <- function(M, ngrid, n_unit) {
  eig <- tryCatch(eigen(M, only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) {
    return(data.frame(
      dominant_mean_N = NA_real_,
      dominant_mean_ploidy = NA_real_,
      dominant_fraction_N_le_25 = NA_real_,
      dominant_fraction_N_below_44 = NA_real_,
      dominant_fraction_N_ge_44 = NA_real_,
      dominant_growth_rate = NA_real_,
      second_growth_rate = NA_real_,
      spectral_gap = NA_real_,
      relative_spectral_gap = NA_real_,
      relax_time_days = NA_real_,
      time_to_10x_days = NA_real_,
      time_to_100x_days = NA_real_,
      log10_advantage_1000d = NA_real_,
      dominance_class = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  fixo2_dominant_from_eig(eig, ngrid, n_unit)
}

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

fixo2_mode_fields <- function(dominant_ploidy) {
  dominant_ploidy <- suppressWarnings(as.numeric(dominant_ploidy))
  threshold <- fixo2_mode_threshold()
  regime <- ifelse(
    !is.finite(dominant_ploidy),
    NA_character_,
    ifelse(dominant_ploidy >= threshold, fixo2_mode_regimes()[["mode1"]], fixo2_mode_regimes()[["mode2"]])
  )
  data.frame(
    trajectory_regime = regime,
    mode_label = names(fixo2_mode_regimes())[match(regime, unname(fixo2_mode_regimes()))],
    mode_source = "fixed_o2_attractor_dominant_ploidy",
    mode_rule = "dominant_mean_ploidy >= 2 => mode1; dominant_mean_ploidy < 2 => mode2",
    mode_threshold_dominant_ploidy = threshold,
    stringsAsFactors = FALSE
  )
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
    "seed_id", "O2_pct", "O2_key", "dominant_mean_ploidy", "population_average_p_misseg", "trajectory_regime",
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

fixo2_format_o2_list <- function(x, max_n = 18L) {
  x <- sort(unique(as.numeric(x)))
  labs <- format(x, scientific = FALSE, trim = TRUE)
  if (length(labs) > max_n) labs <- c(labs[seq_len(max_n)], paste0("... (", length(x), " total)"))
  paste(labs, collapse = ",")
}

fixo2_validate_mode_reference_o2 <- function(mode_reference_o2, attractor_o2_grid) {
  mode_reference_o2 <- suppressWarnings(as.numeric(mode_reference_o2))
  attractor_o2_grid <- sort(unique(suppressWarnings(as.numeric(attractor_o2_grid))))
  attractor_o2_grid <- attractor_o2_grid[is.finite(attractor_o2_grid)]
  if (!is.finite(mode_reference_o2)) {
    stop("--mode_reference_o2 must be a finite numeric O2 value.", call. = FALSE)
  }
  if (!length(attractor_o2_grid)) {
    stop("--attractor_o2_grid must contain at least one finite numeric O2 value.", call. = FALSE)
  }
  if (!any(abs(attractor_o2_grid - mode_reference_o2) < 1e-9)) {
    stop(
      "--mode_reference_o2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " is invalid. It must exactly match one value in --attractor_o2_grid. Available attractor O2 values: ",
      fixo2_format_o2_list(attractor_o2_grid),
      call. = FALSE
    )
  }
  invisible(mode_reference_o2)
}

fixo2_population_average_p_misseg <- function(state, ngrid, O2, run_params, model_env) {
  state <- as.numeric(state)
  ngrid <- as.numeric(ngrid)
  if (length(state) != length(ngrid)) {
    stop("state and ngrid must have the same length.", call. = FALSE)
  }
  if (!length(state) || any(!is.finite(ngrid)) || any(!is.finite(state))) {
    return(NA_real_)
  }
  state_sum <- sum(state)
  if (!is.finite(state_sum) || state_sum <= 0) return(NA_real_)
  state <- state / state_sum
  p_state <- as.numeric(o2ipa_call_model(
    model_env,
    ".pmisseg_of_O2",
    O2 = rep(O2, length(ngrid)),
    run_params = run_params,
    N = ngrid
  ))
  if (length(p_state) != length(state) || any(!is.finite(p_state))) {
    return(NA_real_)
  }
  weighted <- sum(state * p_state)
  if (!is.finite(weighted)) NA_real_ else min(max(weighted, 0), 1)
}

fixo2_dominant_attractor_one <- function(seed_id, run_params, model_env, cfg, O2) {
  fm <- tryCatch(fixo2_fixed_matrix(model_env, cfg, run_params, O2), error = function(e) e)
  if (inherits(fm, "error")) {
    return(data.frame(
      seed_id = seed_id,
      O2_pct = O2,
      status = paste0("matrix_failed:", conditionMessage(fm)),
      dominant_mean_N = NA_real_,
      dominant_mean_ploidy = NA_real_,
      population_average_p_misseg = NA_real_,
      dominant_fraction_N_le_25 = NA_real_,
      dominant_fraction_N_below_44 = NA_real_,
      dominant_fraction_N_ge_44 = NA_real_,
      dominant_growth_rate = NA_real_,
      spectral_gap = NA_real_,
      eigenvector_nonnegative = NA,
      selection_22_vs_44 = NA_real_,
      selection_44_vs_88 = NA_real_,
      eff_growth_N22 = NA_real_,
      eff_growth_N44 = NA_real_,
      eff_growth_N88 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  eig <- tryCatch(eigen(as.matrix(fm$M), only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) {
    return(data.frame(
      seed_id = seed_id,
      O2_pct = O2,
      status = "eigen_failed",
      dominant_mean_N = NA_real_,
      dominant_mean_ploidy = NA_real_,
      population_average_p_misseg = NA_real_,
      dominant_fraction_N_le_25 = NA_real_,
      dominant_fraction_N_below_44 = NA_real_,
      dominant_fraction_N_ge_44 = NA_real_,
      dominant_growth_rate = NA_real_,
      spectral_gap = NA_real_,
      eigenvector_nonnegative = NA,
      selection_22_vs_44 = NA_real_,
      selection_44_vs_88 = NA_real_,
      eff_growth_N22 = NA_real_,
      eff_growth_N44 = NA_real_,
      eff_growth_N88 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  idx <- which.max(Re(eig$values))
  v <- Re(eig$vectors[, idx])
  if (sum(v, na.rm = TRUE) < 0) v <- -v
  nonneg <- all(v >= -1e-8, na.rm = TRUE)
  v <- fixo2_normalize_state(v)
  lambda1 <- Re(eig$values[[idx]])
  lambda2 <- sort(Re(eig$values), decreasing = TRUE)[min(2L, length(eig$values))]
  eff <- vapply(c(22L, 44L, 88L), function(N) {
    col <- as.integer(N - (cfg$N_MIN %||% 22L) + 1L)
    if (col < 1L || col > ncol(fm$G)) return(NA_real_)
    sum(fm$G[, col]) - fm$mu_all[[col]]
  }, numeric(1))
  names(eff) <- c("22", "44", "88")
  data.frame(
    seed_id = seed_id,
    O2_pct = O2,
    status = if (all(is.na(v))) "empty_dominant_vector_after_truncation" else "ok",
    dominant_mean_N = sum(fm$ngrid * v, na.rm = TRUE),
    dominant_mean_ploidy = sum(fm$ngrid * v, na.rm = TRUE) / as.numeric(cfg$N_UNIT %||% 22),
    population_average_p_misseg = fixo2_population_average_p_misseg(
      state = v,
      ngrid = fm$ngrid,
      O2 = O2,
      run_params = run_params,
      model_env = model_env
    ),
    dominant_fraction_N_le_25 = sum(v[fm$ngrid <= 25], na.rm = TRUE),
    dominant_fraction_N_below_44 = sum(v[fm$ngrid < 44], na.rm = TRUE),
    dominant_fraction_N_ge_44 = sum(v[fm$ngrid >= 44], na.rm = TRUE),
    dominant_growth_rate = lambda1,
    spectral_gap = lambda1 - lambda2,
    eigenvector_nonnegative = nonneg,
    selection_22_vs_44 = eff[["22"]] - eff[["44"]],
    selection_44_vs_88 = eff[["44"]] - eff[["88"]],
    eff_growth_N22 = eff[["22"]],
    eff_growth_N44 = eff[["44"]],
    eff_growth_N88 = eff[["88"]],
    stringsAsFactors = FALSE
  )
}

generate_fixo2_attractor_mode_table <- function(run_dir, o2_values, seed_ids = NULL, n_workers = 1L) {
  if (is.null(run_dir) || !nzchar(run_dir) || !dir.exists(run_dir)) {
    stop("run_dir is required to generate FixO2 attractor mode table: ", run_dir)
  }
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  seeds <- if (is.null(seed_ids) || !length(seed_ids)) rownames(param_mat) else intersect(o2ipa_norm_seed(seed_ids), rownames(param_mat))
  if (!length(seeds)) stop("No seed parameters were found for FixO2 attractor mode generation.")
  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  o2_values <- sort(unique(as.numeric(o2_values)))
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, length(seeds))
  message("Generating FixO2 attractor mode table: ", length(seeds), " seeds, ", length(o2_values), " O2 values, workers=", n_workers)
  worker <- function(seed) {
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    rows <- lapply(o2_values, function(O2) fixo2_dominant_attractor_one(seed, run_params, model_env, cfg, O2))
    do.call(rbind, rows)
  }
  rows <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    parallel::mclapply(seeds, worker, mc.cores = n_workers)
  } else {
    lapply(seeds, worker)
  }
  attractors <- do.call(rbind, rows[vapply(rows, nrow, integer(1)) > 0L])
  if (is.null(attractors) || !nrow(attractors)) return(data.frame())
  manifest <- inputs$manifest
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  attractors <- merge(
    attractors,
    manifest[, intersect(c("seed_id", "objective", "delta_objective"), names(manifest)), drop = FALSE],
    by = "seed_id",
    all.x = TRUE,
    sort = FALSE
  )
  attractors <- fixo2_assign_attractor_modes(attractors, "dominant_mean_ploidy")
  fixo2_attractor_mode_table(attractors)
}

# Backward-compatible aliases used by fixed_o2 and dense-grid analyses.
cf2_fixed_matrix <- fixo2_fixed_matrix
cf2_init_vector <- fixo2_init_vector
cf2_normalize_state <- fixo2_normalize_state
cf2_eigen_trajectory <- fixo2_eigen_trajectory
cf2_dominant_one <- fixo2_dominant_one
fo2_dominant_attractor_one <- fixo2_dominant_attractor_one

make_time_grid <- function(time_days, dt, report_dt) {
  time_days <- as.numeric(time_days)
  dt <- as.numeric(dt)
  report_dt <- as.numeric(report_dt)
  if (!is.finite(time_days) || time_days < 0) stop("time_days must be finite and >= 0.")
  if (!is.finite(dt) || dt <= 0) stop("dt must be finite and > 0.")
  if (!is.finite(report_dt) || report_dt <= 0) stop("report_dt must be finite and > 0.")
  step_max <- as.integer(round(time_days / dt))
  keep_steps <- unique(as.integer(round(seq(0, time_days, by = report_dt) / dt)))
  keep_steps <- sort(unique(c(0L, keep_steps, step_max)))
  keep_steps <- keep_steps[keep_steps >= 0L & keep_steps <= step_max]
  list(
    step_max = step_max,
    keep_steps = keep_steps,
    days = as.numeric(keep_steps) * dt,
    requested_time_days = time_days,
    actual_end_day = as.numeric(step_max) * dt
  )
}

build_rate_table <- function(cfg, run_params, fixed_o2) {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  O2_vec <- rep(as.numeric(fixed_o2), length(grid_pre))
  lam <- as.numeric(.lambda_eff_of_O2(
    O2 = O2_vec,
    run_params = run_params,
    N = grid_pre,
    O2_growth = isTRUE(run_params$O2_growth)
  ))
  mu <- as.numeric(.mu_eff_of_O2(O2 = O2_vec, run_params = run_params, N = grid_pre))
  pmiss <- as.numeric(.pmisseg_of_O2(O2 = O2_vec, run_params = run_params, N = grid_pre))
  pwgd <- as.numeric(.p_wgd_of_O2(O2 = O2_vec, run_params = run_params))

  tri <- cpp_o2simps_build_G_for_o2_triplet(
    O2 = as.numeric(fixed_o2),
    O2_crit = as.numeric(.first_non_null_local(run_params$O2_crit, 1.0)),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    lam_max = as.numeric(run_params$lam_max),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = 1e-8,
    buffer_smax = as.numeric(.first_non_null_local(run_params$buffer_smax, 1.0)),
    buffer_beta = as.numeric(.first_non_null_local(run_params$buffer_beta, 0.0)),
    buffer_n_exp = as.numeric(.first_non_null_local(run_params$buffer_n_exp, 1.0)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(run_params$beta_size, 0.0)),
    O2_growth = isTRUE(run_params$O2_growth),
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, 0.0)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, 1.0)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, 0.0)),
    gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu, 1.0)),
    n_O = as.numeric(.first_non_null_local(run_params$n_O, 1.0)),
    ploidy_O2_death = as.character(.first_non_null_local(run_params$ploidy_O2_death, cfg$ploidy_O2_death))
  )
  G <- Matrix::sparseMatrix(
    i = as.integer(tri$i),
    j = as.integer(tri$j),
    x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
    repr = "C"
  )
  net_live_rate <- as.numeric(Matrix::colSums(G)) - mu

  data.frame(
    N = as.integer(grid_pre),
    ploidy = as.numeric(grid_pre) / as.numeric(cfg$N_UNIT),
    fixed_o2_pct = as.numeric(fixed_o2),
    lambda_growth_rate = lam,
    mu_hypoxia_death_rate = mu,
    p_miss = pmiss,
    p_wgd = pwgd,
    dead_buffer_rate = as.numeric(tri$dead_buffer_rate),
    misseg_nonviable_rate = as.numeric(tri$misseg_nonviable_rate),
    boundary_dropped_rate = as.numeric(tri$boundary_dropped_rate),
    net_live_rate_approx = net_live_rate,
    stringsAsFactors = FALSE
  )
}

make_sim_args <- function(init_state, cfg, run_params, time_grid, fixed_o2) {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg))
  list(
    init_state = as.numeric(init_state),
    init_dead_hypoxia_state = rep(0, length(grid_pre)),
    init_dead_buffer_state = rep(0, length(grid_pre)),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(time_grid$keep_steps),
    sim_end_step = as.integer(time_grid$step_max),
    DT = as.numeric(cfg$DT),
    dose = 0.0,
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = as.numeric(time_grid$step_max * cfg$DT + cfg$DT),
    fit_treatment = FALSE,
    alpha = 0.0,
    gamma = 1.0,
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding_enabled = isTRUE(cfg$Crowding),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(.first_non_null_local(run_params$O2_crit, 1.0)),
    o2_feedback = FALSE,
    o2_S0 = as.numeric(fixed_o2),
    kappa_O = as.numeric(.first_non_null_local(run_params$kappa_O, 0.0)),
    tau_O2 = as.numeric(.first_non_null_local(run_params$tau_O2, 0.1)),
    o2_Nref = as.numeric(.first_non_null_local(run_params$o2_Nref, cfg$o2_Nref, 1e6)),
    o2_min = as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.0)),
    o2_S0_upper_bound = as.numeric(cfg$o2_S0_upper_bound),
    eta_o2 = as.numeric(.first_non_null_local(run_params$eta_o2, 1.0)),
    o2_cache_bin_pct = as.numeric(cfg$o2_cache_bin_pct),
    o2_cache_hysteresis_pct = as.numeric(cfg$o2_cache_hysteresis_pct),
    o2_cache_profile = FALSE,
    O2_growth = isTRUE(run_params$O2_growth),
    lam_max = as.numeric(run_params$lam_max),
    p_mis_base = as.numeric(run_params$p_mis_base),
    p_misseg = as.numeric(run_params$p_misseg),
    k_o_mis = as.numeric(run_params$k_o_mis),
    p_wgd = as.numeric(run_params$p_wgd),
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = 1e-8,
    buffer_smax = as.numeric(run_params$buffer_smax),
    buffer_beta = as.numeric(run_params$buffer_beta),
    buffer_n_exp = as.numeric(run_params$buffer_n_exp),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(run_params$beta_size, 0.0)),
    alpha_o2 = as.numeric(run_params$alpha_o2),
    gamma_growth = as.numeric(run_params$gamma_growth),
    mu_hp = as.numeric(run_params$mu_hp),
    gamma_mu = as.numeric(run_params$gamma_mu),
    n_O = as.numeric(run_params$n_O),
    ploidy_O2_death = as.character(.first_non_null_local(run_params$ploidy_O2_death, cfg$ploidy_O2_death)),
    start_with = as.character(cfg$start_with),
    k_clear = as.numeric(run_params$k_clear),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(cfg$burden_log_eps),
    return_full_trajectory = TRUE
  )
}

matrix_long <- function(mat, sim_id, time_grid, grid_pre, cfg, status, rate_table, live_by_time) {
  n_time <- nrow(mat)
  n_state <- length(grid_pre)
  out <- data.frame(
    simulation_id = as.integer(sim_id),
    step = rep(as.integer(time_grid$keep_steps), each = n_state),
    day = rep(as.numeric(time_grid$days), each = n_state),
    N = rep(as.integer(grid_pre), times = n_time),
    ploidy = rep(as.numeric(grid_pre) / as.numeric(cfg$N_UNIT), times = n_time),
    status = status,
    cell_count = as.numeric(t(mat)),
    stringsAsFactors = FALSE
  )
  live_tot <- rep(as.numeric(live_by_time), each = n_state)
  if (identical(status, "live")) {
    out$fraction_of_live_cells <- out$cell_count / pmax(live_tot, 1e-300)
  } else {
    out$fraction_of_live_cells <- NA_real_
  }
  out$total_live_cells <- live_tot
  rate_rep <- rate_table[rep(seq_len(nrow(rate_table)), times = n_time), , drop = FALSE]
  cbind(out, rate_rep[, setdiff(names(rate_rep), c("N", "ploidy")), drop = FALSE])
}

simulate_one <- function(sim_id, init_state, cfg, run_params, time_grid, fixed_o2, rate_table) {
  sim_cpp <- cpp_o2simps_simulate_one(make_sim_args(
    init_state = init_state,
    cfg = cfg,
    run_params = run_params,
    time_grid = time_grid,
    fixed_o2 = fixed_o2
  ))
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  live_state <- as.matrix(sim_cpp$live_state_obs)
  dead_hypoxia_state <- as.matrix(sim_cpp$dead_hypoxia_state_obs)
  dead_buffer_state <- as.matrix(sim_cpp$dead_buffer_state_obs)
  expected_dim <- c(length(time_grid$keep_steps), length(grid_pre))
  if (!identical(dim(live_state), expected_dim)) stop("Unexpected live_state_obs shape.")
  if (!identical(dim(dead_hypoxia_state), expected_dim)) stop("Unexpected dead_hypoxia_state_obs shape.")
  if (!identical(dim(dead_buffer_state), expected_dim)) stop("Unexpected dead_buffer_state_obs shape.")

  population <- data.frame(
    simulation_id = as.integer(sim_id),
    step = as.integer(time_grid$keep_steps),
    day = as.numeric(time_grid$days),
    requested_time_days = as.numeric(time_grid$requested_time_days),
    actual_end_day = as.numeric(time_grid$actual_end_day),
    fixed_o2_pct = as.numeric(fixed_o2),
    o2_target_pct = as.numeric(sim_cpp$O2_target_obs),
    o2_eff_pct = as.numeric(sim_cpp$O2_eff_obs),
    live_cells = as.numeric(sim_cpp$Ntot_live_obs),
    dead_hypoxia_cells = as.numeric(sim_cpp$Ntot_dead_hypoxia_obs),
    dead_buffer_cells = as.numeric(sim_cpp$Ntot_dead_buffer_obs),
    dead_total_cells = as.numeric(sim_cpp$Ntot_dead_total_obs),
    total_cells = as.numeric(sim_cpp$Ntot_total_obs),
    live_volume_mm3 = as.numeric(sim_cpp$Vmm3_live_obs),
    dead_hypoxia_volume_mm3 = as.numeric(sim_cpp$Vmm3_dead_hypoxia_obs),
    dead_buffer_volume_mm3 = as.numeric(sim_cpp$Vmm3_dead_buffer_obs),
    dead_total_volume_mm3 = as.numeric(sim_cpp$Vmm3_dead_total_obs),
    total_volume_mm3 = as.numeric(sim_cpp$Vmm3_total_obs),
    stringsAsFactors = FALSE
  )

  list(
    population = population,
    live_state = live_state,
    dead_hypoxia_state = dead_hypoxia_state,
    dead_buffer_state = dead_buffer_state,
    live_by_time = rowSums(live_state, na.rm = TRUE)
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
  invisible(path)
}

append_tsv <- function(x, path, append = file.exists(path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = !isTRUE(append),
    append = isTRUE(append),
    na = "NA"
  )
  invisible(path)
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
  invisible(path)
}

append_tsv_gz <- function(x, path, append = file.exists(path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = if (isTRUE(append)) "at" else "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = !isTRUE(append),
    na = "NA"
  )
  invisible(path)
}

make_rate_trajectory <- function(rate_table, time_grid, sim_id) {
  cbind(
    data.frame(
      simulation_id = as.integer(sim_id),
      step = rep(as.integer(time_grid$keep_steps), each = nrow(rate_table)),
      day = rep(as.numeric(time_grid$days), each = nrow(rate_table)),
      stringsAsFactors = FALSE
    ),
    rate_table[rep(seq_len(nrow(rate_table)), times = length(time_grid$keep_steps)), , drop = FALSE]
  )
}

params_to_table <- function(run_params, input_names, fixed_o2, source_label) {
  nms <- sort(names(run_params))
  value <- vapply(run_params[nms], function(x) paste(as.character(x), collapse = ","), character(1))
  src <- ifelse(nms %in% input_names, source_label, "default_or_config")
  src[nms %in% c("o2_S0")] <- "fixed_o2_override"
  data.frame(parameter = nms, value = value, source = src, stringsAsFactors = FALSE)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  if (as_bool(argv$help, FALSE) || as_bool(argv$h, FALSE)) {
    usage()
    return(invisible(NULL))
  }
  simulation <- normalize_simulation_mode(argv$simulation %||% "")
  o2_values <- parse_o2_values(argv)
  initial_ploidy_values <- parse_initial_ploidy_values(argv)
  time_days <- as_num(argv$time_days, NA_real_)
  n_sim <- as_int(argv$n_sim, NA_integer_)
  if (!is.finite(n_sim) || is.na(n_sim) || n_sim < 1L) stop("--n_sim must be an integer >= 1.")
  if (!is.finite(time_days) || time_days < 0) stop("--time_days must be finite and >= 0.")

  fit_dir <- resolve_path_value(argv$fit_dir, getwd())
  best_params_path <- resolve_path_value(argv$best_params, getwd())
  run_dir <- resolve_batch_run_dir(argv, fit_dir = fit_dir)
  worker_task <- as_bool(argv$worker_task, FALSE)
  n_core_requested <- as_int(argv$n_core, 1L)
  n_core_batch_requested <- is.finite(n_core_requested) && !is.na(n_core_requested) && n_core_requested > 1L
  batch_requested <- !isTRUE(worker_task) && (
    !is.null(run_dir) ||
      length(o2_values) > 1L ||
      length(initial_ploidy_values) > 1L ||
      length(parse_seed_values(argv$seeds)) > 0L ||
      n_core_batch_requested
  )
  if (batch_requested) {
    return(run_batch(
      argv = argv,
      simulation = simulation,
      fit_dir = fit_dir,
      best_params_path = best_params_path,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      time_days = time_days,
      n_sim = n_sim
    ))
  }

  if (length(o2_values) != 1L) {
    stop("Single-run mode requires exactly one O2 value. Use --run_dir or parent --fit_dir for batch mode.")
  }
  if (length(initial_ploidy_values) != 1L) {
    stop("Single-run mode requires exactly one initial ploidy value. Use --run_dir or parent --fit_dir for batch mode.")
  }
  fixed_o2 <- o2_values[[1]]
  initial_ploidy <- initial_ploidy_values[[1]]
  simulation_ids <- resolve_simulation_ids(argv, n_sim)
  param_info <- resolve_parameter_values(
    fit_dir = fit_dir,
    best_params_path = best_params_path,
    simulation = simulation,
    joint_scope = argv$joint_scope
  )
  cfg_raw <- read_fit_config(fit_dir)
  run_params <- prepare_run_params(
    param_values = param_info$values,
    simulation = simulation,
    cfg = cfg_raw,
    fixed_o2 = fixed_o2
  )
  cfg <- prepare_sim_cfg(cfg_raw, argv, fixed_o2 = fixed_o2, run_params = run_params)
  run_params$O2_growth <- isTRUE(cfg$O2_growth)
  run_params$ploidy_O2_death <- cfg$ploidy_O2_death

  initial_cells <- as_num(argv$initial_cells, 1.0)
  init <- make_initial_state(cfg, initial_ploidy = initial_ploidy, initial_cells = initial_cells)
  save_every_days <- resolve_save_every_days(argv)
  time_grid <- make_time_grid(time_days = time_days, dt = cfg$DT, report_dt = save_every_days)
  rate_table_base <- build_rate_table(cfg = cfg, run_params = run_params, fixed_o2 = fixed_o2)

  seed_id <- resolve_seed_id(argv$seed_id, fit_dir = fit_dir, best_params_path = best_params_path)
  out_dir <- resolve_path_value(argv$out_dir, getwd())
  if (is.null(out_dir)) {
    out_dir <- default_output_dir(
      simulation = simulation,
      fixed_o2 = fixed_o2,
      initial_ploidy = initial_ploidy,
      seed_id = seed_id
    )
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("Fixed-O2 simulation")
  message("  simulation: ", simulation)
  message("  parameter_source: ", param_info$source, " (", param_info$path, ")")
  message("  fixed_o2_pct: ", fixed_o2)
  message("  initial_ploidy: ", init$requested_ploidy, " -> grid N=", init$used_N, " ploidy=", signif(init$used_ploidy, 8))
  message("  seed_id: ", seed_id)
  message("  time_days: ", time_days, " dt=", cfg$DT, " save_every_days=", save_every_days, " n_sim=", n_sim)
  message("  output: ", out_dir)

  output_paths <- fixo2_simulation_output_paths(out_dir)
  if (!as_bool(argv$force, FALSE) && fixo2_simulation_output_complete(out_dir)) {
    message("  existing fixed-O2 simulation output is complete; skipping. Use --force=TRUE to overwrite.")
    return(invisible(out_dir))
  }
  unlink(unlist(output_paths), force = TRUE)

  run_config <- data.frame(
    field = c(
      "simulation", "fit_dir", "parameter_file", "parameter_source",
      "fixed_o2_pct", "o2_feedback", "initial_ploidy_requested",
      "initial_N_used", "initial_ploidy_used", "initial_cells",
      "time_days_requested", "dt", "save_every_days", "report_dt", "n_sim",
      "simulation_ids", "actual_end_day",
      "seed_id", "output_dir",
      "N_MIN", "N_MAX", "N_UNIT", "start_with", "ploidy_O2_death",
      "O2_growth", "Crowding"
    ),
    value = as.character(c(
      simulation, fit_dir %||% "", param_info$path, param_info$source,
      fixed_o2, FALSE, init$requested_ploidy,
      init$used_N, init$used_ploidy, initial_cells,
      time_days, cfg$DT, save_every_days, save_every_days, n_sim,
      paste(simulation_ids, collapse = ","), time_grid$actual_end_day,
      seed_id, normalizePath(out_dir, mustWork = FALSE),
      cfg$N_MIN, cfg$N_MAX, cfg$N_UNIT, cfg$start_with, cfg$ploidy_O2_death,
      cfg$O2_growth, cfg$Crowding
    )),
    stringsAsFactors = FALSE
  )
  parameters_used <- params_to_table(
    run_params = run_params,
    input_names = names(param_info$values),
    fixed_o2 = fixed_o2,
    source_label = param_info$source
  )
  write_tsv(run_config, output_paths$run_config)
  write_tsv(parameters_used, output_paths$parameters_used)

  population_n <- 0L
  state_n <- 0
  rate_n <- 0
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  for (sim_index in seq_along(simulation_ids)) {
    sim_id <- simulation_ids[[sim_index]]
    message("  replicate ", sim_id, " (", sim_index, "/", length(simulation_ids), ")")
    res <- simulate_one(
      sim_id = sim_id,
      init_state = init$state,
      cfg = cfg,
      run_params = run_params,
      time_grid = time_grid,
      fixed_o2 = fixed_o2,
      rate_table = rate_table_base
    )

    if (any(abs(res$population$o2_eff_pct - fixed_o2) > 1e-10, na.rm = TRUE) ||
        any(abs(res$population$o2_target_pct - fixed_o2) > 1e-10, na.rm = TRUE)) {
      stop("Internal fixed-O2 check failed: O2 target/effective trajectory is not constant at requested O2.")
    }
    append_tsv(res$population, output_paths$population)
    population_n <- population_n + nrow(res$population)

    rate_chunk <- make_rate_trajectory(rate_table_base, time_grid, sim_id)
    append_tsv_gz(rate_chunk, output_paths$rate)
    rate_n <- rate_n + nrow(rate_chunk)
    rm(rate_chunk)
    gc(verbose = FALSE)

    state_chunk <- matrix_long(res$live_state, sim_id, time_grid, grid_pre, cfg, "live", rate_table_base, res$live_by_time)
    append_tsv_gz(state_chunk, output_paths$state)
    state_n <- state_n + nrow(state_chunk)
    rm(state_chunk)
    gc(verbose = FALSE)

    state_chunk <- matrix_long(res$dead_hypoxia_state, sim_id, time_grid, grid_pre, cfg, "dead_hypoxia", rate_table_base, res$live_by_time)
    append_tsv_gz(state_chunk, output_paths$state)
    state_n <- state_n + nrow(state_chunk)
    rm(state_chunk)
    gc(verbose = FALSE)

    state_chunk <- matrix_long(res$dead_buffer_state, sim_id, time_grid, grid_pre, cfg, "dead_buffer", rate_table_base, res$live_by_time)
    append_tsv_gz(state_chunk, output_paths$state)
    state_n <- state_n + nrow(state_chunk)
    rm(state_chunk, res)
    gc(verbose = FALSE)
  }

  message("Done.")
  message("  population_trajectory.tsv rows: ", population_n)
  message("  rate_trajectory.tsv.gz rows: ", rate_n)
  message("  state_trajectory.tsv.gz rows: ", state_n)
  invisible(out_dir)
}


# ---- embedded regression-smoothed curve classifier ----------------------
curve_classification_rule_version <- function() {
  "shape_first_v3"
}

curve_step_epsilon <- function(ploidy_range,
                               step_epsilon_abs = 1e-6,
                               step_epsilon_fraction = 1e-4) {
  scale <- if (is.finite(ploidy_range) && ploidy_range > 0) ploidy_range else 1
  max(step_epsilon_abs, step_epsilon_fraction * scale)
}

collapse_signs <- function(signs) {
  s <- signs[is.finite(signs) & !is.na(signs)]
  nz <- s[s != 0L]
  if (!length(nz)) return(integer(0))
  rle(nz)$values
}

finite_diff_curve <- function(curve,
                              step_epsilon = NULL,
                              value_col = "dominant_mean_ploidy",
                              x_col = "O2_pct",
                              id_col = "seed_id",
                              step_epsilon_abs = 1e-6,
                              step_epsilon_fraction = 1e-4) {
  if (!value_col %in% names(curve)) stop("Missing value column: ", value_col, call. = FALSE)
  if (!x_col %in% names(curve)) stop("Missing x column: ", x_col, call. = FALSE)
  curve <- curve[order(curve[[x_col]]), , drop = FALSE]
  y <- suppressWarnings(as.numeric(curve[[value_col]]))
  pr <- suppressWarnings(max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
  if (!is.finite(pr)) pr <- NA_real_
  if (is.null(step_epsilon) || !is.finite(step_epsilon)) {
    step_epsilon <- curve_step_epsilon(pr, step_epsilon_abs, step_epsilon_fraction)
  }
  dy <- c(diff(y), NA_real_)
  sign <- rep(NA_integer_, length(dy))
  sign[is.finite(dy) & dy > step_epsilon] <- 1L
  sign[is.finite(dy) & dy < -step_epsilon] <- -1L
  sign[is.finite(dy) & abs(dy) <= step_epsilon] <- 0L
  out <- data.frame(
    O2_pct = curve[[x_col]],
    finite_difference_next = dy,
    local_slope_sign = sign,
    step_epsilon = rep(step_epsilon, length(dy)),
    stringsAsFactors = FALSE
  )
  if (!is.null(id_col) && id_col %in% names(curve)) {
    out <- cbind(data.frame(seed_id = curve[[id_col]], stringsAsFactors = FALSE), out)
  }
  out
}

terminal_plateau_from_signs <- function(signs, plateau_min_points = 3L) {
  signs <- signs[!is.na(signs)]
  if (!length(signs)) return(FALSE)
  last_nz <- suppressWarnings(max(which(signs != 0L)))
  if (!is.finite(last_nz)) return(FALSE)
  (length(signs) - last_nz) >= plateau_min_points
}

classify_o2_ploidy_curve <- function(curve,
                                     value_col = "dominant_mean_ploidy",
                                     x_col = "O2_pct",
                                     id_col = "seed_id",
                                     flat_range_threshold = 0.05,
                                     step_epsilon_abs = 1e-6,
                                     step_epsilon_fraction = 1e-4,
                                     reverse_fraction_tolerance = 0.05,
                                     plateau_min_points = 3L) {
  if (!value_col %in% names(curve)) stop("Missing value column: ", value_col, call. = FALSE)
  if (!x_col %in% names(curve)) stop("Missing x column: ", x_col, call. = FALSE)
  curve <- curve[order(curve[[x_col]]), , drop = FALSE]
  y <- suppressWarnings(as.numeric(curve[[value_col]]))
  finite_y <- y[is.finite(y)]
  if (length(finite_y) < 2L) {
    step_epsilon <- curve_step_epsilon(NA_real_, step_epsilon_abs, step_epsilon_fraction)
    diffs <- finite_diff_curve(curve, step_epsilon, value_col, x_col, id_col, step_epsilon_abs, step_epsilon_fraction)
    summary <- data.frame(
      curve_class = "insufficient_data",
      sign_sequence = "",
      n_sign_changes = NA_integer_,
      step_epsilon = step_epsilon,
      slope_epsilon = step_epsilon,
      flat_range_threshold = flat_range_threshold,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      ploidy_range = NA_real_,
      net_ploidy_change = NA_real_,
      max_positive_step = NA_real_,
      max_negative_step = NA_real_,
      positive_step_total = NA_real_,
      negative_step_total = NA_real_,
      fraction_positive_steps = NA_real_,
      fraction_negative_steps = NA_real_,
      fraction_zero_steps = NA_real_,
      low_amplitude_curve = NA,
      terminal_plateau = NA,
      terminal_plateau_for_class = NA,
      classification_rule_version = curve_classification_rule_version(),
      stringsAsFactors = FALSE
    )
    return(list(summary = summary, differences = diffs, collapsed_signs = integer(0)))
  }

  ploidy_range <- max(finite_y, na.rm = TRUE) - min(finite_y, na.rm = TRUE)
  finite_idx <- which(is.finite(y))
  net_change <- y[finite_idx[[length(finite_idx)]]] - y[finite_idx[[1L]]]
  step_epsilon <- curve_step_epsilon(ploidy_range, step_epsilon_abs, step_epsilon_fraction)
  diffs <- finite_diff_curve(curve, step_epsilon, value_col, x_col, id_col, step_epsilon_abs, step_epsilon_fraction)
  dy <- diffs$finite_difference_next
  signs <- diffs$local_slope_sign
  step_ok <- is.finite(dy) & !is.na(signs)
  dy <- dy[step_ok]
  signs <- signs[step_ok]
  collapsed <- collapse_signs(signs)
  n_changes <- if (length(collapsed) <= 1L) 0L else length(collapsed) - 1L
  pos_total <- sum(dy[signs > 0], na.rm = TRUE)
  neg_total <- sum(abs(dy[signs < 0]), na.rm = TRUE)
  frac_pos <- if (length(signs)) mean(signs > 0) else NA_real_
  frac_neg <- if (length(signs)) mean(signs < 0) else NA_real_
  frac_zero <- if (length(signs)) mean(signs == 0) else NA_real_
  max_pos <- suppressWarnings(max(dy[signs > 0], na.rm = TRUE))
  max_neg <- suppressWarnings(max(abs(dy[signs < 0]), na.rm = TRUE))
  if (!is.finite(max_pos)) max_pos <- 0
  if (!is.finite(max_neg)) max_neg <- 0
  terminal_plateau <- terminal_plateau_from_signs(signs, plateau_min_points)
  low_amplitude_curve <- is.finite(ploidy_range) && ploidy_range <= flat_range_threshold
  terminal_plateau_for_class <- terminal_plateau && !low_amplitude_curve

  reverse_ok_for_increase <- neg_total <= reverse_fraction_tolerance * max(pos_total, ploidy_range, .Machine$double.eps)
  reverse_ok_for_decrease <- pos_total <= reverse_fraction_tolerance * max(neg_total, ploidy_range, .Machine$double.eps)

  curve_class <- "complex_nonmonotone"
  if (!length(collapsed)) {
    curve_class <- if (is.finite(net_change) && net_change > 0) {
      "monotone_increasing"
    } else if (is.finite(net_change) && net_change < 0) {
      "monotone_decreasing"
    } else {
      "approximately_flat"
    }
  } else if (length(collapsed) == 1L && collapsed[[1L]] == 1L) {
    curve_class <- if (terminal_plateau_for_class) "single_transition_increase_then_plateau" else "monotone_increasing"
  } else if (length(collapsed) == 1L && collapsed[[1L]] == -1L) {
    curve_class <- if (terminal_plateau_for_class) "single_transition_decrease_then_plateau" else "monotone_decreasing"
  } else if (length(collapsed) == 2L && identical(as.integer(collapsed), c(-1L, 1L))) {
    curve_class <- "u_shaped"
  } else if (length(collapsed) == 2L && identical(as.integer(collapsed), c(1L, -1L))) {
    curve_class <- "inverted_u_shaped"
  } else if (is.finite(net_change) && net_change > 0 && pos_total > 0 && reverse_ok_for_increase) {
    curve_class <- "monotone_increasing"
  } else if (is.finite(net_change) && net_change < 0 && neg_total > 0 && reverse_ok_for_decrease) {
    curve_class <- "monotone_decreasing"
  }

  summary <- data.frame(
    curve_class = curve_class,
    sign_sequence = if (length(collapsed)) paste(collapsed, collapse = ",") else "",
    n_sign_changes = n_changes,
    step_epsilon = step_epsilon,
    slope_epsilon = step_epsilon,
    flat_range_threshold = flat_range_threshold,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    ploidy_range = ploidy_range,
    net_ploidy_change = net_change,
    max_positive_step = max_pos,
    max_negative_step = max_neg,
    positive_step_total = pos_total,
    negative_step_total = neg_total,
    fraction_positive_steps = frac_pos,
    fraction_negative_steps = frac_neg,
    fraction_zero_steps = frac_zero,
    low_amplitude_curve = low_amplitude_curve,
    terminal_plateau = terminal_plateau,
    terminal_plateau_for_class = terminal_plateau_for_class,
    classification_rule_version = curve_classification_rule_version(),
    stringsAsFactors = FALSE
  )
  list(summary = summary, differences = diffs, collapsed_signs = collapsed)
}

run_curve_classification_validation <- function() {
  mk <- function(seed, y, gap = 0.02) {
    data.frame(
      seed_id = seed,
      O2_pct = seq_along(y),
      dominant_mean_ploidy = y,
      spectral_gap = rep(gap, length(y)),
      stringsAsFactors = FALSE
    )
  }
  cases <- list(
    flat = list(y = rep(2, 201), expected = "approximately_flat"),
    tiny_increasing = list(y = seq(2, 2.0002, length.out = 201), expected = "monotone_increasing"),
    tiny_u = list(y = c(seq(2.0004, 2, length.out = 101), seq(2.000004, 2.0002, length.out = 100)), expected = "u_shaped"),
    small_increasing = list(y = seq(2, 2.2, length.out = 201), expected = "monotone_increasing"),
    small_decreasing = list(y = seq(2.2, 2, length.out = 201), expected = "monotone_decreasing"),
    increasing = list(y = seq(1, 3, length.out = 201), expected = "monotone_increasing"),
    decreasing = list(y = seq(3, 1, length.out = 201), expected = "monotone_decreasing"),
    increase_plateau = list(y = c(seq(1, 2.2, length.out = 120), rep(2.2, 81)), expected = "single_transition_increase_then_plateau"),
    decrease_plateau = list(y = c(seq(3, 1.8, length.out = 120), rep(1.8, 81)), expected = "single_transition_decrease_then_plateau"),
    u = list(y = c(seq(3, 1.7, length.out = 101), seq(1.72, 3, length.out = 100)), expected = "u_shaped"),
    inverted_u = list(y = c(seq(1, 2.8, length.out = 101), seq(2.78, 1, length.out = 100)), expected = "inverted_u_shaped"),
    complex = list(y = c(seq(1, 2, length.out = 67), seq(1.98, 1.4, length.out = 67), seq(1.42, 2.4, length.out = 67)), expected = "complex_nonmonotone")
  )
  rows <- lapply(names(cases), function(nm) {
    z <- mk(paste0("seed_", nm), cases[[nm]]$y)
    obs <- classify_o2_ploidy_curve(z)$summary$curve_class[[1L]]
    data.frame(
      test_case = nm,
      expected_class = cases[[nm]]$expected,
      observed_class = obs,
      passed = identical(cases[[nm]]$expected, obs),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

smooth_curve_classification_rule_version <- function() {
  "loess_persistent_v1"
}

smooth_curve_values <- function(x, y,
                                span = 0.20,
                                degree = 2L,
                                family = "gaussian",
                                spline_spar = 0.65) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  out <- rep(NA_real_, length(y))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 3L) {
    out[ok] <- y[ok]
    return(out)
  }

  dat <- data.frame(x = x[ok], y = y[ok])
  span <- min(max(as.numeric(span), 0.05), 1)
  degree <- as.integer(degree)
  if (!degree %in% c(0L, 1L, 2L)) degree <- 2L
  family <- as.character(family)
  if (!family %in% c("gaussian", "symmetric")) family <- "gaussian"

  plausible_prediction <- function(pred) {
    if (length(pred) != sum(ok) || !all(is.finite(pred))) return(FALSE)
    raw_range <- max(dat$y, na.rm = TRUE) - min(dat$y, na.rm = TRUE)
    if (!is.finite(raw_range)) raw_range <- 0
    margin <- max(0.05, 0.10 * raw_range)
    pred_min <- min(pred, na.rm = TRUE)
    pred_max <- max(pred, na.rm = TRUE)
    pred_range <- pred_max - pred_min
    is.finite(pred_min) && is.finite(pred_max) &&
      pred_min >= min(dat$y, na.rm = TRUE) - margin &&
      pred_max <= max(dat$y, na.rm = TRUE) + margin &&
      pred_range <= raw_range + 2 * margin
  }

  loess_predict <- function(fam) {
    fit <- tryCatch(
      suppressWarnings(stats::loess(
        y ~ x,
        data = dat,
        span = span,
        degree = degree,
        family = fam,
        control = stats::loess.control(surface = "direct", trace.hat = "approximate")
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) return(rep(NA_real_, sum(ok)))
    tryCatch(
      suppressWarnings(as.numeric(stats::predict(fit, newdata = data.frame(x = x[ok])))),
      error = function(e) rep(NA_real_, sum(ok))
    )
  }

  families <- unique(c(family, "gaussian"))
  for (fam in families) {
    pred <- loess_predict(fam)
    if (plausible_prediction(pred)) {
      out[ok] <- pred
      break
    }
  }

  if (any(!is.finite(out[ok]))) {
    spl <- tryCatch(
      stats::smooth.spline(dat$x, dat$y, spar = spline_spar),
      error = function(e) NULL
    )
    if (!is.null(spl)) {
      pred <- tryCatch(
        as.numeric(stats::predict(spl, x[ok])$y),
        error = function(e) rep(NA_real_, sum(ok))
      )
      out[ok] <- pred
    }
  }

  bad <- which(ok & !is.finite(out))
  if (length(bad)) out[bad] <- y[bad]
  out
}

persistent_sign_segments <- function(x, y,
                                     signs,
                                     min_segment_points,
                                     min_segment_span,
                                     min_segment_amplitude) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  signs <- as.integer(signs)
  n <- length(x)
  if (n < 2L) {
    return(data.frame())
  }
  signs <- signs[seq_len(n - 1L)]
  signs[is.na(signs)] <- 0L
  runs <- rle(signs)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  rows <- list()
  for (i in seq_along(runs$values)) {
    s <- runs$values[[i]]
    if (!s %in% c(-1L, 1L)) next
    st <- starts[[i]]
    en <- ends[[i]]
    if (en + 1L > n) next
    x_span <- x[[en + 1L]] - x[[st]]
    signed_amplitude <- y[[en + 1L]] - y[[st]]
    amplitude <- abs(signed_amplitude)
    persistent <- is.finite(x_span) && is.finite(amplitude) &&
      (en - st + 1L) >= min_segment_points &&
      x_span >= min_segment_span &&
      amplitude >= min_segment_amplitude
    rows[[length(rows) + 1L]] <- data.frame(
      segment_index = length(rows) + 1L,
      sign = s,
      start_interval = st,
      end_interval = en,
      n_steps = en - st + 1L,
      x_start = x[[st]],
      x_end = x[[en + 1L]],
      x_span = x_span,
      signed_amplitude = signed_amplitude,
      amplitude = amplitude,
      persistent = persistent,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$kept_index <- NA_integer_
  keep <- which(out$persistent)
  if (length(keep)) out$kept_index[keep] <- seq_along(keep)
  out
}

classify_o2_ploidy_curve_smooth <- function(curve,
                                            value_col = "dominant_mean_ploidy",
                                            x_col = "O2_pct",
                                            id_col = "seed_id",
                                            flat_range_threshold = 0.05,
                                            step_epsilon_abs = 1e-6,
                                            step_epsilon_fraction = 1e-4,
                                            reverse_fraction_tolerance = 0.05,
                                            smooth_span = 0.20,
                                            smooth_degree = 2L,
                                            smooth_family = "gaussian",
                                            min_segment_span_fraction = 0.02,
                                            min_segment_amplitude_abs = 0.01,
                                            min_segment_amplitude_fraction = 0.03,
                                            min_segment_points = 3L,
                                            terminal_plateau_span_fraction = 0.10,
                                            terminal_plateau_amplitude_fraction = 0.03) {
  if (!value_col %in% names(curve)) stop("Missing value column: ", value_col, call. = FALSE)
  if (!x_col %in% names(curve)) stop("Missing x column: ", x_col, call. = FALSE)
  curve <- curve[order(curve[[x_col]]), , drop = FALSE]
  x <- suppressWarnings(as.numeric(curve[[x_col]]))
  y <- suppressWarnings(as.numeric(curve[[value_col]]))
  finite_y <- y[is.finite(y)]
  if (length(finite_y) < 2L) {
    step_epsilon <- curve_step_epsilon(NA_real_, step_epsilon_abs, step_epsilon_fraction)
    diffs <- data.frame(
      O2_pct = x,
      raw_value = y,
      fitted_value = y,
      finite_difference_next = c(diff(y), NA_real_),
      fitted_difference_next = c(diff(y), NA_real_),
      local_slope_sign = rep(NA_integer_, length(y)),
      fitted_local_slope_sign = rep(NA_integer_, length(y)),
      step_epsilon = rep(step_epsilon, length(y)),
      stringsAsFactors = FALSE
    )
    if (!is.null(id_col) && id_col %in% names(curve)) {
      diffs <- cbind(data.frame(seed_id = curve[[id_col]], stringsAsFactors = FALSE), diffs)
    }
    summary <- data.frame(
      curve_class = "insufficient_data",
      sign_sequence = "",
      n_sign_changes = NA_integer_,
      classification_rule_version = smooth_curve_classification_rule_version(),
      raw_ploidy_range = NA_real_,
      fitted_ploidy_range = NA_real_,
      net_ploidy_change = NA_real_,
      step_epsilon = step_epsilon,
      flat_range_threshold = flat_range_threshold,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      smooth_span = smooth_span,
      smooth_degree = smooth_degree,
      min_segment_points = min_segment_points,
      min_segment_span = NA_real_,
      min_segment_amplitude = NA_real_,
      terminal_plateau = NA,
      terminal_plateau_for_class = NA,
      positive_persistent_total = NA_real_,
      negative_persistent_total = NA_real_,
      n_persistent_segments = NA_integer_,
      stringsAsFactors = FALSE
    )
    return(list(summary = summary, differences = diffs, segments = data.frame()))
  }

  fitted <- smooth_curve_values(
    x = x,
    y = y,
    span = smooth_span,
    degree = smooth_degree,
    family = smooth_family
  )
  raw_range <- max(finite_y, na.rm = TRUE) - min(finite_y, na.rm = TRUE)
  finite_fitted <- fitted[is.finite(fitted)]
  fitted_range <- max(finite_fitted, na.rm = TRUE) - min(finite_fitted, na.rm = TRUE)
  if (!is.finite(fitted_range)) fitted_range <- NA_real_
  finite_idx <- which(is.finite(fitted))
  net_change <- fitted[finite_idx[[length(finite_idx)]]] - fitted[finite_idx[[1L]]]
  step_epsilon <- curve_step_epsilon(fitted_range, step_epsilon_abs, step_epsilon_fraction)
  raw_dy <- c(diff(y), NA_real_)
  fitted_dy <- c(diff(fitted), NA_real_)
  raw_sign <- rep(NA_integer_, length(raw_dy))
  fitted_sign <- rep(NA_integer_, length(fitted_dy))
  raw_sign[is.finite(raw_dy) & raw_dy > step_epsilon] <- 1L
  raw_sign[is.finite(raw_dy) & raw_dy < -step_epsilon] <- -1L
  raw_sign[is.finite(raw_dy) & abs(raw_dy) <= step_epsilon] <- 0L
  fitted_sign[is.finite(fitted_dy) & fitted_dy > step_epsilon] <- 1L
  fitted_sign[is.finite(fitted_dy) & fitted_dy < -step_epsilon] <- -1L
  fitted_sign[is.finite(fitted_dy) & abs(fitted_dy) <= step_epsilon] <- 0L

  x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (!is.finite(x_range) || x_range <= 0) x_range <- 1
  min_segment_points_eff <- max(as.integer(min_segment_points), ceiling(min_segment_span_fraction * max(length(x) - 1L, 1L)))
  min_segment_span <- min_segment_span_fraction * x_range
  min_segment_amplitude <- max(min_segment_amplitude_abs, min_segment_amplitude_fraction * fitted_range)
  if (!is.finite(min_segment_amplitude)) min_segment_amplitude <- min_segment_amplitude_abs

  seg <- persistent_sign_segments(
    x = x,
    y = fitted,
    signs = fitted_sign,
    min_segment_points = min_segment_points_eff,
    min_segment_span = min_segment_span,
    min_segment_amplitude = min_segment_amplitude
  )
  kept <- if (nrow(seg)) seg[seg$persistent, , drop = FALSE] else data.frame()
  kept_signs <- if (nrow(kept)) as.integer(kept$sign) else integer(0)
  collapsed <- if (length(kept_signs)) rle(kept_signs)$values else integer(0)
  n_changes <- if (length(collapsed) <= 1L) 0L else length(collapsed) - 1L
  pos_total <- if (nrow(kept)) sum(kept$amplitude[kept$sign > 0], na.rm = TRUE) else 0
  neg_total <- if (nrow(kept)) sum(kept$amplitude[kept$sign < 0], na.rm = TRUE) else 0

  terminal_plateau <- FALSE
  if (nrow(kept)) {
    last_end <- kept$end_interval[[nrow(kept)]]
    anchor <- min(last_end + 1L, length(x))
    trailing_span <- x[[length(x)]] - x[[anchor]]
    trailing_amplitude <- abs(fitted[[length(fitted)]] - fitted[[anchor]])
    terminal_plateau <- is.finite(trailing_span) && is.finite(trailing_amplitude) &&
      trailing_span >= terminal_plateau_span_fraction * x_range &&
      trailing_amplitude <= max(min_segment_amplitude, terminal_plateau_amplitude_fraction * fitted_range)
  }
  low_amplitude_curve <- is.finite(fitted_range) && fitted_range <= flat_range_threshold
  terminal_plateau_for_class <- terminal_plateau && !low_amplitude_curve
  reverse_ok_for_increase <- neg_total <= reverse_fraction_tolerance * max(pos_total, fitted_range, .Machine$double.eps)
  reverse_ok_for_decrease <- pos_total <= reverse_fraction_tolerance * max(neg_total, fitted_range, .Machine$double.eps)

  curve_class <- "complex_nonmonotone"
  if (low_amplitude_curve) {
    curve_class <- "approximately_flat"
  } else if (!length(collapsed)) {
    curve_class <- if (is.finite(net_change) && net_change > flat_range_threshold) {
      "monotone_increasing"
    } else if (is.finite(net_change) && net_change < -flat_range_threshold) {
      "monotone_decreasing"
    } else {
      "approximately_flat"
    }
  } else if (length(collapsed) == 1L && collapsed[[1L]] == 1L) {
    curve_class <- if (terminal_plateau_for_class) "single_transition_increase_then_plateau" else "monotone_increasing"
  } else if (length(collapsed) == 1L && collapsed[[1L]] == -1L) {
    curve_class <- if (terminal_plateau_for_class) "single_transition_decrease_then_plateau" else "monotone_decreasing"
  } else if (length(collapsed) == 2L && identical(as.integer(collapsed), c(-1L, 1L))) {
    curve_class <- "u_shaped"
  } else if (length(collapsed) == 2L && identical(as.integer(collapsed), c(1L, -1L))) {
    curve_class <- "inverted_u_shaped"
  } else if (is.finite(net_change) && net_change > 0 && pos_total > 0 && reverse_ok_for_increase) {
    curve_class <- "monotone_increasing"
  } else if (is.finite(net_change) && net_change < 0 && neg_total > 0 && reverse_ok_for_decrease) {
    curve_class <- "monotone_decreasing"
  }

  diffs <- data.frame(
    O2_pct = x,
    raw_value = y,
    fitted_value = fitted,
    finite_difference_next = raw_dy,
    fitted_difference_next = fitted_dy,
    local_slope_sign = raw_sign,
    fitted_local_slope_sign = fitted_sign,
    step_epsilon = rep(step_epsilon, length(y)),
    stringsAsFactors = FALSE
  )
  if (!is.null(id_col) && id_col %in% names(curve)) {
    diffs <- cbind(data.frame(seed_id = curve[[id_col]], stringsAsFactors = FALSE), diffs)
  }
  summary <- data.frame(
    curve_class = curve_class,
    sign_sequence = if (length(collapsed)) paste(collapsed, collapse = ",") else "",
    n_sign_changes = n_changes,
    classification_rule_version = smooth_curve_classification_rule_version(),
    raw_ploidy_range = raw_range,
    fitted_ploidy_range = fitted_range,
    net_ploidy_change = net_change,
    step_epsilon = step_epsilon,
    flat_range_threshold = flat_range_threshold,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    smooth_span = smooth_span,
    smooth_degree = smooth_degree,
    min_segment_points = min_segment_points_eff,
    min_segment_span = min_segment_span,
    min_segment_amplitude = min_segment_amplitude,
    terminal_plateau = terminal_plateau,
    terminal_plateau_for_class = terminal_plateau_for_class,
    positive_persistent_total = pos_total,
    negative_persistent_total = neg_total,
    n_persistent_segments = nrow(kept),
    stringsAsFactors = FALSE
  )
  list(summary = summary, differences = diffs, segments = seg)
}

run_smooth_curve_classification_validation <- function() {
  mk <- function(seed, y) {
    data.frame(
      seed_id = seed,
      O2_pct = seq(0, 5, length.out = length(y)),
      dominant_mean_ploidy = y,
      stringsAsFactors = FALSE
    )
  }
  y_inc_with_spike <- seq(1, 3, length.out = 201)
  y_inc_with_spike[100] <- y_inc_with_spike[100] - 0.35
  cases <- list(
    flat_noise = list(y = 2 + 0.004 * sin(seq(0, 8 * pi, length.out = 201)), expected = "approximately_flat"),
    increasing_spike = list(y = y_inc_with_spike, expected = "monotone_increasing"),
    increasing = list(y = seq(1, 3, length.out = 201), expected = "monotone_increasing"),
    decreasing = list(y = seq(3, 1, length.out = 201), expected = "monotone_decreasing"),
    u = list(y = c(seq(3, 1.6, length.out = 101), seq(1.62, 3, length.out = 100)), expected = "u_shaped"),
    inverted_u = list(y = c(seq(1, 2.8, length.out = 101), seq(2.78, 1, length.out = 100)), expected = "inverted_u_shaped"),
    increase_plateau = list(y = c(seq(1, 2.2, length.out = 115), rep(2.2, 86)), expected = "single_transition_increase_then_plateau")
  )
  rows <- lapply(names(cases), function(nm) {
    z <- mk(paste0("seed_", nm), cases[[nm]]$y)
    obs <- classify_o2_ploidy_curve_smooth(z)$summary$curve_class[[1L]]
    data.frame(
      test_case = nm,
      expected_class = cases[[nm]]$expected,
      observed_class = obs,
      passed = identical(cases[[nm]]$expected, obs),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# =============================================================================
# Oxygen response standalone workflow
# This section follows the embedded analytical model and curve
# classification functions. It contains every Oxygen response-specific calculation,
# validation, plotting, assembly, and provenance operation.
# =============================================================================

RESPONSE_SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

response_parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!grepl("^--[^=]+=", arg)) next
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    out[[gsub("-", "_", key, fixed = TRUE)]] <- sub("^[^=]+=", "", arg)
  }
  out
}

response_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  val <- tolower(trimws(as.character(x[[1L]])))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

response_as_int <- function(x, default = 1L) {
  out <- suppressWarnings(as.integer(if (is.null(x)) default else x[[1L]]))
  if (!is.finite(out) || out < 1L) default else out
}

response_read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  data.table::fread(path, data.table = FALSE, check.names = FALSE)
}

response_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  normalizePath(path, mustWork = TRUE)
}

response_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)[[1L]])
}

response_seed_number <- function(x) {
  suppressWarnings(as.integer(sub("^seed", "", basename(as.character(x)))))
}

response_read_metric_map <- function(path) {
  tab <- response_read_tsv(path)
  if (!all(c("metric", "value") %in% names(tab))) {
    stop("Invalid fit_summary.tsv: ", path, call. = FALSE)
  }
  tab <- tab[!duplicated(tab$metric), , drop = FALSE]
  stats::setNames(as.character(tab$value), as.character(tab$metric))
}

response_metric_number <- function(metric_map, metric, path) {
  value <- suppressWarnings(as.numeric(metric_map[[metric]]))
  if (length(value) != 1L || !is.finite(value)) {
    stop("Missing/non-finite metric '", metric, "' in ", path, call. = FALSE)
  }
  value
}

response_save_plot_pair <- function(plot, base_path, width, height) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  png_path <- paste0(base_path, ".png")
  pdf_path <- paste0(base_path, ".pdf")
  ggplot2::ggsave(
    png_path, plot, width = width, height = height, units = "in",
    dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    pdf_path, plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
  c(png = normalizePath(png_path, mustWork = TRUE),
    pdf = normalizePath(pdf_path, mustWork = TRUE))
}

response_curve_class_order <- c(
  "complex_nonmonotone",
  "inverted_u_shaped",
  "monotone_increasing",
  "u_shaped",
  "approximately_flat",
  "single_transition_increase_then_plateau",
  "single_transition_decrease_then_plateau",
  "monotone_decreasing"
)

response_curve_class_colors <- c(
  complex_nonmonotone = "#6A5ACD",
  inverted_u_shaped = "#D55E00",
  monotone_increasing = "#009E73",
  u_shaped = "#E69F00",
  approximately_flat = "#666666",
  single_transition_increase_then_plateau = "#CC79A7",
  single_transition_decrease_then_plateau = "#56B4E9",
  monotone_decreasing = "#0072B2"
)

response_region_colors <- c(
  C01 = "#C99700",
  C02 = "#6A3D9A",
  C03 = "#006D2C",
  C04 = "#0072B2",
  C05 = "#D55E00",
  C06 = "#009E73"
)
response_region_label_outline_color <- "#111111"
response_region_label_outline_radius <- 0.025

response_font_family <- "Helvetica"
response_font_sizes <- c(
  panel_title = 10.5,
  subtitle = 7.5,
  axis_title = 9,
  axis_text = 8,
  legend_title = 8,
  legend_text = 7.2,
  panel_tag = 14
)

response_curve_class_labels <- c(
  complex_nonmonotone = "Complex nonmonotone",
  inverted_u_shaped = "Inverted U",
  monotone_increasing = "Monotone increasing",
  u_shaped = "U-shaped",
  approximately_flat = "Approximately flat",
  single_transition_increase_then_plateau = "Increase then plateau",
  single_transition_decrease_then_plateau = "Decrease then plateau",
  monotone_decreasing = "Monotone decreasing"
)

response_expected_class_counts <- c(
  complex_nonmonotone = 316L,
  inverted_u_shaped = 50L,
  monotone_increasing = 50L,
  u_shaped = 47L,
  approximately_flat = 26L,
  single_transition_increase_then_plateau = 5L,
  single_transition_decrease_then_plateau = 4L,
  monotone_decreasing = 2L
)

# -----------------------------------------------------------------------------
# Oxygen responseA: regression-smoothed fixed-O2 response classes
# -----------------------------------------------------------------------------

response_classify_all_curves <- function(curves) {
  groups <- split(curves, curves$seed_id)
  summaries <- vector("list", length(groups))
  smooth_rows <- vector("list", length(groups))
  segment_rows <- vector("list", length(groups))
  group_names <- names(groups)

  for (i in seq_along(groups)) {
    curve <- groups[[i]]
    curve <- curve[order(curve$O2_pct), , drop = FALSE]
    result <- classify_o2_ploidy_curve_smooth(
      curve,
      value_col = "dominant_mean_ploidy",
      x_col = "O2_pct",
      id_col = "seed_id",
      flat_range_threshold = 0.05,
      step_epsilon_abs = 1e-6,
      step_epsilon_fraction = 1e-4,
      reverse_fraction_tolerance = 0.05,
      smooth_span = 0.20,
      smooth_degree = 2L,
      smooth_family = "gaussian",
      min_segment_span_fraction = 0.02,
      min_segment_amplitude_abs = 0.01,
      min_segment_amplitude_fraction = 0.03,
      min_segment_points = 3L,
      terminal_plateau_span_fraction = 0.10,
      terminal_plateau_amplitude_fraction = 0.03
    )

    summary_row <- result$summary
    summary_row$seed_id <- group_names[[i]]
    summary_row$seed_number <- response_seed_number(group_names[[i]])
    summaries[[i]] <- summary_row

    diff_row <- result$differences
    diff_row$seed_id <- group_names[[i]]
    diff_row$seed_number <- response_seed_number(group_names[[i]])
    smooth_rows[[i]] <- diff_row

    if (nrow(result$segments)) {
      seg <- result$segments
      seg$seed_id <- group_names[[i]]
      seg$seed_number <- response_seed_number(group_names[[i]])
      segment_rows[[i]] <- seg
    }
  }

  summary <- data.table::rbindlist(summaries, fill = TRUE)
  smooth <- data.table::rbindlist(smooth_rows, fill = TRUE)
  segments <- data.table::rbindlist(segment_rows, fill = TRUE)
  summary <- as.data.frame(summary)
  smooth <- as.data.frame(smooth)
  segments <- as.data.frame(segments)

  if (!"fitted_value" %in% names(smooth)) {
    stop("Curve classifier did not return fitted_value.", call. = FALSE)
  }
  names(smooth)[names(smooth) == "fitted_value"] <- "smoothed_dominant_mean_ploidy"
  names(summary)[names(summary) == "curve_class"] <- "smooth_curve_class"

  summary <- summary[order(summary$seed_number), , drop = FALSE]
  smooth <- smooth[order(smooth$seed_number, smooth$O2_pct), , drop = FALSE]
  list(summary = summary, smooth = smooth, segments = segments)
}

response_adaptive_lwd <- function(n) {
  max(0.18, min(1.05, 2.8 / sqrt(max(n, 1))))
}

response_panel_y_range <- function(y, min_span = 0.2, pad_frac = 0.08) {
  y <- y[is.finite(y)]
  if (!length(y)) return(c(0, 1))
  yr <- range(y)
  span <- diff(yr)
  if (!is.finite(span) || span < min_span) {
    center <- mean(yr)
    span <- min_span
    yr <- center + c(-0.5, 0.5) * span
  }
  pad <- max(span * pad_frac, .Machine$double.eps)
  yr + c(-pad, pad)
}

response_draw_panel_a <- function(smooth, by_seed, png_path, pdf_path) {
  merged <- merge(
    smooth,
    by_seed[, c("seed_id", "smooth_curve_class"), drop = FALSE],
    by = "seed_id", all.x = TRUE, sort = FALSE
  )
  class_counts <- table(factor(by_seed$smooth_curve_class, levels = response_curve_class_order))
  classes <- response_curve_class_order[class_counts > 0]

  draw <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(
      mfrow = c(4, 2),
      mar = c(0.65, 3.15, 0.42, 0.45),
      oma = c(2.35, 2.8, 1.30, 0.15),
      family = response_font_family
    )
    for (panel_index in seq_along(classes)) {
      class_name <- classes[[panel_index]]
      z <- merged[merged$smooth_curve_class == class_name, , drop = FALSE]
      seeds <- unique(z$seed_id)
      y_range <- response_panel_y_range(z$smoothed_dominant_mean_ploidy)
      graphics::plot(
        range(z$O2_pct, na.rm = TRUE), y_range, type = "n",
        xlab = "", ylab = "", main = "",
        xaxt = if (panel_index > 6L) "s" else "n",
        cex.axis = response_font_sizes[["axis_text"]] / 10,
        mgp = c(1.65, 0.45, 0),
        tcl = -0.22
      )
      graphics::grid(col = "#DDDDDD", lty = 3)
      lwd <- if (identical(class_name, "approximately_flat")) {
        max(0.85, response_adaptive_lwd(length(seeds)))
      } else {
        response_adaptive_lwd(length(seeds))
      }
      alpha <- if (identical(class_name, "approximately_flat")) {
        0.75
      } else {
        max(0.18, min(0.55, 10 / sqrt(max(length(seeds), 1))))
      }
      curve_color <- grDevices::adjustcolor(
        response_curve_class_colors[[class_name]], alpha.f = alpha
      )
      for (seed in seeds) {
        zz <- z[z$seed_id == seed, , drop = FALSE]
        zz <- zz[order(zz$O2_pct), , drop = FALSE]
        graphics::lines(
          zz$O2_pct, zz$smoothed_dominant_mean_ploidy,
          col = curve_color, lwd = lwd
        )
      }
      median_curve <- stats::aggregate(
        smoothed_dominant_mean_ploidy ~ O2_pct, z, stats::median, na.rm = TRUE
      )
      graphics::lines(
        median_curve$O2_pct, median_curve$smoothed_dominant_mean_ploidy,
        col = "#111111", lwd = 1.2
      )
      title_text <- paste0(
        response_curve_class_labels[[class_name]],
        " (n=", length(seeds), ")"
      )
      graphics::legend(
        "topright", legend = title_text,
        bty = "n", bg = grDevices::adjustcolor("white", alpha.f = 0.82),
        text.font = 2, text.col = "#222222",
        cex = response_font_sizes[["panel_title"]] / 10,
        inset = c(0.006, 0.012), x.intersp = 0.3, y.intersp = 0.8
      )
    }
    graphics::mtext(
      "Fixed oxygen (%)", side = 1, outer = TRUE,
      line = 1.15, cex = response_font_sizes[["axis_title"]] / 10
    )
    graphics::mtext(
      "Smoothed dominant mean ploidy", side = 2, outer = TRUE,
      line = 1.3, cex = response_font_sizes[["axis_title"]] / 10
    )
    graphics::mtext(
      "A. Oxygen-ploidy response classes",
      side = 3, outer = TRUE, at = 0.055, adj = 0,
      line = 0.08, font = 2,
      cex = response_font_sizes[["panel_title"]] / 10
    )
  }

  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(
    png_path, width = 8.6, height = 6.5,
    units = "in", res = 300, pointsize = 10, bg = "white"
  )
  draw()
  grDevices::dev.off()
  grDevices::cairo_pdf(
    pdf_path, width = 8.6, height = 6.5,
    pointsize = 10, family = response_font_family, bg = "white"
  )
  draw()
  grDevices::dev.off()
  c(png = normalizePath(png_path, mustWork = TRUE),
    pdf = normalizePath(pdf_path, mustWork = TRUE))
}

response_build_panel_a <- function(input_path, data_dir, figure_dir) {
  curves <- response_read_tsv(input_path)
  required <- c("seed_id", "O2_pct", "dominant_mean_ploidy")
  if (!all(required %in% names(curves))) {
    stop("Oxygen responseA input is missing: ",
         paste(setdiff(required, names(curves)), collapse = ", "), call. = FALSE)
  }
  curves$seed_id <- as.character(curves$seed_id)
  curves$O2_pct <- as.numeric(curves$O2_pct)
  curves$dominant_mean_ploidy <- as.numeric(curves$dominant_mean_ploidy)

  if (nrow(curves) != 100500L ||
      length(unique(curves$seed_id)) != 500L ||
      length(unique(curves$O2_pct)) != 201L ||
      min(curves$O2_pct) != 0 || max(curves$O2_pct) != 5 ||
      anyDuplicated(curves[, c("seed_id", "O2_pct")]) ||
      any(!is.finite(curves$dominant_mean_ploidy))) {
    stop("Oxygen responseA requires 500 seeds x 201 finite O2 points from 0 to 5.",
         call. = FALSE)
  }

  classified <- response_classify_all_curves(curves)
  class_table <- classified$summary
  class_counts <- as.data.frame(
    table(factor(class_table$smooth_curve_class, levels = response_curve_class_order)),
    stringsAsFactors = FALSE
  )
  names(class_counts) <- c("smooth_curve_class", "n_seed")
  class_counts <- class_counts[class_counts$n_seed > 0, , drop = FALSE]
  class_counts$fraction_seed <- class_counts$n_seed / sum(class_counts$n_seed)
  observed_counts <- stats::setNames(
    as.integer(class_counts$n_seed), as.character(class_counts$smooth_curve_class)
  )
  if (!identical(observed_counts[names(response_expected_class_counts)],
                 response_expected_class_counts)) {
    stop(
      "Oxygen responseA class counts differ from the validated target: ",
      paste(names(observed_counts), observed_counts, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }

  paths <- c(
    smoothed_curves = response_write_tsv(
      classified$smooth,
      file.path(data_dir, "response_class_smoothed_curves.tsv")
    ),
    class_by_seed = response_write_tsv(
      class_table,
      file.path(data_dir, "response_class_curve_class_by_seed.tsv")
    ),
    class_counts = response_write_tsv(
      class_counts,
      file.path(data_dir, "response_class_class_counts.tsv")
    ),
    persistent_segments = response_write_tsv(
      classified$segments,
      file.path(data_dir, "response_class_persistent_segments.tsv")
    )
  )
  figure_paths <- response_draw_panel_a(
    classified$smooth, class_table,
    file.path(figure_dir, "response_class_fixed_o2_response_classes.png"),
    file.path(figure_dir, "response_class_fixed_o2_response_classes.pdf")
  )
  list(
    curves = curves, class_by_seed = class_table, class_counts = class_counts,
    data_paths = paths, figure_paths = figure_paths
  )
}

# -----------------------------------------------------------------------------
# Oxygen responseB: pair-specific complete O2 x effective-p_misseg grids
# -----------------------------------------------------------------------------

response_collect_pair_valid_seeds <- function(pair_dir) {
  seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
  rows <- lapply(seed_dirs, function(seed_dir) {
    summary_path <- file.path(seed_dir, "fit_summary.tsv")
    if (!file.exists(summary_path)) return(NULL)
    metric_map <- response_read_metric_map(summary_path)
    status <- as.character(metric_map[["fit_status"]])
    objective <- suppressWarnings(as.numeric(metric_map[["objective"]]))
    parameter_path <- file.path(seed_dir, "joint_best_params_long.tsv")
    if (!identical(status, "ok") || !is.finite(objective) || !file.exists(parameter_path)) {
      return(NULL)
    }
    data.frame(
      pair_id = basename(pair_dir),
      seed = response_seed_number(seed_dir),
      seed_id = basename(seed_dir),
      objective = objective,
      fit_status = status,
      seed_dir = normalizePath(seed_dir, mustWork = TRUE),
      fit_summary_path = normalizePath(summary_path, mustWork = TRUE),
      parameter_path = normalizePath(parameter_path, mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- data.table::rbindlist(rows, fill = TRUE)
  if (!nrow(out)) stop("No valid joint seeds under ", pair_dir, call. = FALSE)
  out <- as.data.frame(out)
  out[order(out$objective, out$seed), , drop = FALSE]
}

response_select_joint_pair_best <- function(joint_root) {
  pair_dirs <- list.dirs(joint_root, recursive = FALSE, full.names = TRUE)
  pair_dirs <- sort(pair_dirs[grepl("^fit_joint_", basename(pair_dirs))])
  if (length(pair_dirs) != 6L) {
    stop("Expected six fit_joint_* pair directories; found ", length(pair_dirs),
         call. = FALSE)
  }
  selected <- do.call(rbind, lapply(pair_dirs, function(pair_dir) {
    response_collect_pair_valid_seeds(pair_dir)[1L, , drop = FALSE]
  }))
  selected <- selected[order(selected$pair_id), , drop = FALSE]
  rownames(selected) <- NULL

  observed <- stats::setNames(
    selected$seed,
    sub("^.*_(C[0-9]{2})_.*$", "\\1", selected$pair_id)
  )
  expected <- c(
    C01 = 111L, C02 = 468L, C03 = 194L,
    C04 = 372L, C05 = 252L, C06 = 270L
  )
  if (!identical(observed[names(expected)], expected)) {
    stop(
      "Best joint seeds differ from validated selection: ",
      paste(names(observed), observed, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  selected
}

response_pair_metadata <- function(pair_id, joint_seed) {
  invivo_seed <- suppressWarnings(as.integer(
    sub("^.*_vi_seed([0-9]+)_.*$", "\\1", pair_id)
  ))
  class_id <- sub("^.*_(C[0-9]{2})_.*$", "\\1", pair_id)
  invitro_seed <- suppressWarnings(as.integer(
    sub("^.*_vt_seed([0-9]+)$", "\\1", pair_id)
  ))
  if (!is.finite(invivo_seed) || identical(class_id, pair_id) ||
      !is.finite(invitro_seed) || !is.finite(joint_seed)) {
    stop("Cannot parse canonical pair metadata from ", pair_id, call. = FALSE)
  }
  data.frame(
    pair_id = pair_id,
    class_id = class_id,
    clone = class_id,
    subclone = "",
    invivo_seed = invivo_seed,
    invitro_seed = invitro_seed,
    joint_seed = as.integer(joint_seed),
    stringsAsFactors = FALSE
  )
}

response_force_effective_p_misseg <- function(run_params, effective_p_misseg) {
  value <- as.numeric(effective_p_misseg)
  if (length(value) != 1L || !is.finite(value) || value <= 0 || value >= 1) {
    stop("effective_p_misseg must be within (0,1).", call. = FALSE)
  }
  out <- run_params
  out$p_mis_base <- value
  out$p_misseg <- 0
  out
}

response_prepare_pair_model <- function(best_row, o2_values) {
  fit_dir <- best_row$seed_dir[[1L]]
  parameter_info <- resolve_parameter_values(
    fit_dir = fit_dir,
    best_params_path = NULL,
    simulation = "joint",
    joint_scope = "shared_invivo"
  )
  config_raw <- read_fit_config(fit_dir)
  run_params <- prepare_run_params(
    param_values = parameter_info$values,
    simulation = "joint",
    cfg = config_raw,
    fixed_o2 = max(o2_values)
  )
  config <- prepare_sim_cfg(
    config_raw, argv = list(), fixed_o2 = max(o2_values),
    run_params = run_params
  )
  run_params$O2_growth <- isTRUE(config$O2_growth)
  run_params$ploidy_O2_death <- config$ploidy_O2_death
  list(
    fit_dir = fit_dir, parameter_info = parameter_info,
    config = config, run_params = run_params
  )
}

response_compute_pair <- function(best_row, o2_values, cin_values) {
  prepared <- response_prepare_pair_model(best_row, o2_values)
  model_env <- globalenv()

  overlay_rows <- lapply(o2_values, function(o2) {
    fixo2_dominant_attractor_one(
      seed_id = best_row$seed_id[[1L]],
      run_params = prepared$run_params,
      model_env = model_env,
      cfg = prepared$config,
      O2 = o2
    )
  })
  overlay <- do.call(rbind, overlay_rows)
  overlay$pair_id <- best_row$pair_id[[1L]]
  overlay$joint_seed <- best_row$seed[[1L]]
  overlay$objective <- best_row$objective[[1L]]
  overlay$source_seed_dir <- prepared$fit_dir
  overlay$source_parameter_file <- prepared$parameter_info$path

  surface_rows <- lapply(cin_values, function(cin_value) {
    forced_params <- response_force_effective_p_misseg(
      prepared$run_params, cin_value
    )
    do.call(rbind, lapply(o2_values, function(o2) {
      row <- fixo2_dominant_attractor_one(
        seed_id = best_row$seed_id[[1L]],
        run_params = forced_params,
        model_env = model_env,
        cfg = prepared$config,
        O2 = o2
      )
      row$actual_effective_p_misseg <- row$population_average_p_misseg
      row$effective_p_misseg <- cin_value
      row$log10_effective_p_misseg <- log10(cin_value)
      row$population_average_p_misseg <- cin_value
      row$log10_population_average_p_misseg <- log10(cin_value)
      row
    }))
  })
  surface <- do.call(rbind, surface_rows)
  surface$pair_id <- best_row$pair_id[[1L]]
  surface$joint_seed <- best_row$seed[[1L]]
  surface$objective <- best_row$objective[[1L]]
  surface$source_seed_dir <- prepared$fit_dir
  surface$source_parameter_file <- prepared$parameter_info$path

  list(overlay = overlay, surface = surface)
}

response_validate_panel_b_data <- function(overlay, surface, selected) {
  if (nrow(overlay) != 1206L ||
      length(unique(overlay$pair_id)) != 6L ||
      anyDuplicated(overlay[, c("pair_id", "O2_pct")]) ||
      any(overlay$status != "ok") ||
      any(!is.finite(overlay$population_average_p_misseg)) ||
      any(!is.finite(overlay$dominant_mean_ploidy)) ||
      any(!is.finite(overlay$spectral_gap))) {
    stop("Oxygen responseB overlay validation failed.", call. = FALSE)
  }

  if (nrow(surface) != 72360L ||
      length(unique(surface$pair_id)) != 6L ||
      anyDuplicated(surface[, c("pair_id", "O2_pct", "effective_p_misseg")]) ||
      any(surface$status != "ok") ||
      any(!is.finite(surface$dominant_mean_ploidy)) ||
      any(!is.finite(surface$spectral_gap))) {
    stop("Oxygen responseB surface validation failed.", call. = FALSE)
  }
  per_pair <- table(surface$pair_id)
  if (any(per_pair != 12060L)) {
    stop("Each Oxygen responseB pair must contain 12,060 surface points.",
         call. = FALSE)
  }
  if (max(abs(surface$actual_effective_p_misseg -
              surface$effective_p_misseg), na.rm = TRUE) > 1e-8) {
    stop("Oxygen responseB forced effective p_misseg audit failed.",
         call. = FALSE)
  }
  if (min(surface$O2_pct) != 0 || max(surface$O2_pct) != 5 ||
      length(unique(surface$O2_pct)) != 201L ||
      length(unique(surface$effective_p_misseg)) != 60L ||
      abs(min(surface$effective_p_misseg) - 0.005) > 1e-12 ||
      abs(max(surface$effective_p_misseg) - 0.5) > 1e-12) {
    stop("Oxygen responseB O2/effective-p_misseg grid validation failed.",
         call. = FALSE)
  }

  pair_check <- merge(
    unique(overlay[, c("pair_id", "joint_seed", "objective")]),
    selected[, c("pair_id", "seed", "objective")],
    by = "pair_id", suffixes = c("_overlay", "_selected")
  )
  if (nrow(pair_check) != 6L ||
      any(pair_check$joint_seed != pair_check$seed) ||
      any(abs(pair_check$objective_overlay -
              pair_check$objective_selected) > 1e-10)) {
    stop("Oxygen responseB best-seed/objective provenance validation failed.",
         call. = FALSE)
  }
  invisible(TRUE)
}

response_scientific_labels <- function(log10_values) {
  formatC(10^log10_values, format = "e", digits = 0)
}

response_surface_plot <- function(surface, overlay, metadata, fill_limits) {
  surface <- surface[order(surface$O2_pct,
                           surface$log10_effective_p_misseg), , drop = FALSE]
  overlay <- overlay[order(overlay$O2_pct), , drop = FALSE]
  overlay$log10_population_average_p_misseg <-
    log10(overlay$population_average_p_misseg)

  tile_width <- min(diff(sort(unique(surface$O2_pct))))
  tile_height <- min(diff(sort(unique(
    surface$log10_effective_p_misseg
  ))))
  y_breaks <- pretty(range(surface$log10_effective_p_misseg), n = 5)

  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = surface,
      ggplot2::aes(
        x = O2_pct, y = log10_effective_p_misseg,
        fill = dominant_mean_ploidy
      ),
      width = tile_width, height = tile_height
    ) +
    ggplot2::geom_contour(
      data = surface,
      ggplot2::aes(
        x = O2_pct, y = log10_effective_p_misseg,
        z = as.numeric(spectral_gap < 0.005)
      ),
      breaks = 0.5, color = "white", linetype = "dotted",
      linewidth = 0.42, show.legend = FALSE
    ) +
    ggplot2::geom_path(
      data = overlay,
      ggplot2::aes(
        x = O2_pct, y = log10_population_average_p_misseg
      ),
      color = "#111111", linewidth = 0.75
    ) +
    ggplot2::geom_point(
      data = overlay[c(1L, nrow(overlay)), , drop = FALSE],
      ggplot2::aes(
        x = O2_pct, y = log10_population_average_p_misseg
      ),
      shape = 21, fill = "white", color = "#111111",
      stroke = 0.6, size = 1.8
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#2C7BB6", "#74ADD1", "#FEE08B", "#F46D43", "#A50026"),
      limits = fill_limits, oob = scales::squish,
      name = "Dominant\nmean ploidy"
    ) +
    ggplot2::scale_x_continuous(breaks = 0:5, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks, labels = response_scientific_labels,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = range(surface$log10_effective_p_misseg), expand = FALSE
    ) +
    ggplot2::labs(
      title = paste(metadata$clone, metadata$subclone),
      subtitle = paste0(
        "Complete O2-effective p_misseg grid\n",
        "vi seed ", metadata$invivo_seed,
        "; vt seed ", metadata$invitro_seed,
        "; best joint seed ", metadata$joint_seed
      ),
      x = "Fixed oxygen (%)",
      y = "Effective per-chromosome\nmissegregation probability p_misseg"
    ) +
    ggplot2::theme_classic(
      base_size = response_font_sizes[["axis_text"]],
      base_family = response_font_family
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = response_font_sizes[["panel_title"]],
        color = "#111111", hjust = 0,
        margin = ggplot2::margin(b = 2)
      ),
      plot.title.position = "panel",
      plot.subtitle = ggplot2::element_text(
        size = response_font_sizes[["subtitle"]],
        color = "#555555", lineheight = 0.9
      ),
      axis.title = ggplot2::element_text(
        size = response_font_sizes[["axis_title"]]
      ),
      axis.text = ggplot2::element_text(
        size = response_font_sizes[["axis_text"]]
      ),
      legend.position = "right",
      legend.title = ggplot2::element_text(
        size = response_font_sizes[["legend_title"]]
      ),
      legend.text = ggplot2::element_text(
        size = response_font_sizes[["legend_text"]]
      ),
      plot.margin = grid::unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
      aspect.ratio = 1
    )
}

response_draw_panel_b <- function(overlay, surface, selected, figure_dir) {
  metadata <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
    response_pair_metadata(selected$pair_id[[i]], selected$seed[[i]])
  }))
  metadata$row <- 1L
  metadata$col <- match(metadata$clone, sprintf("C%02d", 1:6))
  metadata <- metadata[order(metadata$row, metadata$col), , drop = FALSE]
  if (nrow(metadata) != 6L || anyNA(metadata$row) || anyNA(metadata$col) ||
      anyDuplicated(metadata[, c("row", "col")])) {
    stop("Oxygen responseB requires one complete C01-C06 primary-region row.",
         call. = FALSE)
  }

  fill_limits <- range(surface$dominant_mean_ploidy, c(1, 4), na.rm = TRUE)
  plots <- lapply(seq_len(nrow(metadata)), function(i) {
    meta <- metadata[i, , drop = FALSE]
    plot <- response_surface_plot(
      surface[surface$pair_id == meta$pair_id, , drop = FALSE],
      overlay[overlay$pair_id == meta$pair_id, , drop = FALSE],
      meta, fill_limits
    )
    if (meta$col != 1L) {
      plot <- plot + ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
    if (i == 1L) {
      plot <- plot +
        ggplot2::labs(tag = "B") +
        ggplot2::theme(
          plot.tag = ggplot2::element_text(
            face = "bold", size = response_font_sizes[["panel_tag"]]
          )
        )
    }
    plot
  })

  composite <- patchwork::wrap_plots(
    plots, nrow = 1, ncol = 6, byrow = TRUE, guides = "collect"
  ) +
    patchwork::plot_annotation(
      caption = paste0(
        "Pair-specific complete O2-effective p_misseg grids from the best ",
        "joint MAP parameter vectors.\n",
        "Black curves: unmodified fitted trajectory; white dotted contours: ",
        "spectral gap <0.005."
      )
    )

  outputs <- response_save_plot_pair(
    composite,
    file.path(figure_dir, "pair_surface_o2_effective_p_misseg_six_pair_grid"),
    width = 17.2, height = 3.0
  )
  list(plot = composite, figure_paths = outputs, metadata = metadata)
}

response_build_panel_b <- function(joint_root, data_dir, figure_dir,
                               workers = 1L, rebuild = FALSE) {
  selected_path <- file.path(data_dir, "pair_surface_best_joint_seed_summary.tsv")
  overlay_path <- file.path(data_dir, "pair_surface_best_trajectory.tsv")
  surface_path <- file.path(data_dir, "pair_surface_complete_surface.tsv")

  selected <- response_select_joint_pair_best(joint_root)
  use_cache <- !rebuild &&
    file.exists(selected_path) && file.exists(overlay_path) &&
    file.exists(surface_path)

  if (use_cache) {
    cached_selected <- response_read_tsv(selected_path)
    current_key <- selected[, c("pair_id", "seed", "objective"), drop = FALSE]
    cached_key <- cached_selected[, c("pair_id", "seed", "objective"), drop = FALSE]
    current_key <- current_key[order(current_key$pair_id), , drop = FALSE]
    cached_key <- cached_key[order(cached_key$pair_id), , drop = FALSE]
    use_cache <- identical(current_key$pair_id, cached_key$pair_id) &&
      identical(as.integer(current_key$seed), as.integer(cached_key$seed)) &&
      max(abs(current_key$objective - cached_key$objective)) < 1e-12
  }

  if (use_cache) {
    message("Oxygen responseB: reusing script-generated validated grid cache.")
    overlay <- response_read_tsv(overlay_path)
    surface <- response_read_tsv(surface_path)
  } else {
    message("Oxygen responseB: computing six complete 201 x 60 grids from raw joint fits.")
    o2_values <- seq(0, 5, length.out = 201L)
    cin_values <- 10^seq(log10(0.005), log10(0.5), length.out = 60L)
    compute_one <- function(i) {
      message(
        "  pair ", i, "/", nrow(selected), ": ",
        selected$pair_id[[i]], " seed", selected$seed[[i]]
      )
      response_compute_pair(selected[i, , drop = FALSE], o2_values, cin_values)
    }
    rows <- if (workers > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(
        seq_len(nrow(selected)), compute_one,
        mc.cores = min(workers, nrow(selected))
      )
    } else {
      lapply(seq_len(nrow(selected)), compute_one)
    }
    overlay <- do.call(rbind, lapply(rows, `[[`, "overlay"))
    surface <- do.call(rbind, lapply(rows, `[[`, "surface"))
    overlay <- overlay[order(overlay$pair_id, overlay$O2_pct), , drop = FALSE]
    surface <- surface[
      order(surface$pair_id, surface$O2_pct, surface$effective_p_misseg),
      , drop = FALSE
    ]
    response_validate_panel_b_data(overlay, surface, selected)
    response_write_tsv(selected, selected_path)
    response_write_tsv(overlay, overlay_path)
    response_write_tsv(surface, surface_path)
  }

  overlay_keep <- c(
    "pair_id", "joint_seed", "seed_id", "objective", "O2_pct", "status",
    "population_average_p_misseg", "dominant_mean_ploidy", "spectral_gap",
    "dominant_growth_rate"
  )
  surface_keep <- c(
    "pair_id", "joint_seed", "seed_id", "objective", "O2_pct", "status",
    "effective_p_misseg", "log10_effective_p_misseg",
    "population_average_p_misseg", "log10_population_average_p_misseg",
    "actual_effective_p_misseg", "dominant_mean_ploidy", "spectral_gap",
    "dominant_growth_rate"
  )
  if (!all(overlay_keep %in% names(overlay)) ||
      !all(surface_keep %in% names(surface))) {
    stop("Oxygen responseB compact output schema is incomplete.", call. = FALSE)
  }
  overlay_was_compact <- identical(names(overlay), overlay_keep)
  surface_was_compact <- identical(names(surface), surface_keep)
  overlay <- overlay[, overlay_keep, drop = FALSE]
  surface <- surface[, surface_keep, drop = FALSE]
  if (!overlay_was_compact) response_write_tsv(overlay, overlay_path)
  if (!surface_was_compact) response_write_tsv(surface, surface_path)

  response_validate_panel_b_data(overlay, surface, selected)
  drawn <- response_draw_panel_b(overlay, surface, selected, figure_dir)

  qc_rows <- lapply(split(surface, surface$pair_id), function(z) {
    data.frame(
      pair_id = z$pair_id[[1L]],
      joint_seed = unique(z$joint_seed)[[1L]],
      objective = unique(z$objective)[[1L]],
      n_o2 = length(unique(z$O2_pct)),
      n_effective_p_misseg = length(unique(z$effective_p_misseg)),
      n_grid_points = nrow(z),
      n_status_ok = sum(z$status == "ok"),
      max_abs_actual_minus_requested_p_misseg =
        max(abs(z$actual_effective_p_misseg - z$effective_p_misseg)),
      dominant_mean_ploidy_min = min(z$dominant_mean_ploidy),
      dominant_mean_ploidy_max = max(z$dominant_mean_ploidy),
      spectral_gap_min = min(z$spectral_gap),
      spectral_gap_median = stats::median(z$spectral_gap),
      fraction_spectral_gap_below_0p005 = mean(z$spectral_gap < 0.005),
      stringsAsFactors = FALSE
    )
  })
  qc <- do.call(rbind, qc_rows)
  qc_path <- response_write_tsv(
    qc, file.path(data_dir, "pair_surface_pair_qc_summary.tsv")
  )
  list(
    selected = selected, overlay = overlay, surface = surface,
    data_paths = c(
      selected = normalizePath(selected_path, mustWork = TRUE),
      overlay = normalizePath(overlay_path, mustWork = TRUE),
      surface = normalizePath(surface_path, mustWork = TRUE),
      qc = qc_path
    ),
    figure_paths = drawn$figure_paths,
    metadata = drawn$metadata
  )
}

# -----------------------------------------------------------------------------
# Oxygen responseC: response classes, parameter t-SNE, and full MAP fit quality
# -----------------------------------------------------------------------------

response_build_panel_c <- function(
    tsne_path, class_by_seed, data_dir, figure_dir, panel_tag = "C") {
  landscape <- utils::read.csv(
    tsne_path, check.names = FALSE, stringsAsFactors = FALSE
  )
  required <- c(
    "tSNE1", "tSNE2", "dataset", "point_type", "seed", "objective",
    "cluster_id"
  )
  if (!all(required %in% names(landscape))) {
    stop("Oxygen responseC t-SNE input is missing: ",
         paste(setdiff(required, names(landscape)), collapse = ", "),
         call. = FALSE)
  }
  landscape <- landscape[
    landscape$dataset == "invivo" & landscape$point_type == "best",
    required, drop = FALSE
  ]
  landscape$seed <- as.integer(landscape$seed)
  landscape$map_objective <- as.numeric(landscape$objective)
  landscape$objective <- NULL
  classes <- class_by_seed[, c(
    "seed_number", "smooth_curve_class"
  ), drop = FALSE]
  names(classes) <- c("seed", "curve_class")
  data <- merge(landscape, classes, by = "seed", all.x = TRUE, sort = FALSE)
  data <- data[order(data$seed), , drop = FALSE]

  if (nrow(data) != 500L || anyDuplicated(data$seed) ||
      any(!is.finite(data$tSNE1)) || any(!is.finite(data$tSNE2)) ||
      any(!is.finite(data$map_objective)) || anyNA(data$curve_class)) {
    stop("Oxygen responseC requires 500 unique, finite, classified in-vivo best points.",
         call. = FALSE)
  }
  unknown <- setdiff(unique(data$curve_class), response_curve_class_order)
  if (length(unknown)) {
    stop("Unknown Oxygen responseC response class: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  data$curve_class <- factor(
    data$curve_class, levels = response_curve_class_order
  )
  data$warm_start_region <- factor(
    sub("^vi_", "", as.character(data$cluster_id)),
    levels = names(response_region_colors)
  )
  expected_region_counts <- c(
    C01 = 84L, C02 = 18L, C03 = 165L,
    C04 = 46L, C05 = 160L, C06 = 27L
  )
  observed_region_counts <- table(data$warm_start_region)
  if (anyNA(data$warm_start_region) ||
      !identical(
        as.integer(observed_region_counts[names(expected_region_counts)]),
        as.integer(expected_region_counts)
      )) {
    stop(
      "Oxygen responseC warm-start-region validation failed: ",
      paste(
        names(observed_region_counts),
        as.integer(observed_region_counts),
        sep = "=", collapse = ", "
      ),
      call. = FALSE
    )
  }
  cutoff <- as.numeric(stats::quantile(
    data$map_objective, 0.10, names = FALSE, type = 7
  ))
  data$lowest_objective_decile <- data$map_objective <= cutoff
  if (sum(data$lowest_objective_decile) != 50L ||
      abs(cutoff - 14.3326611383349) > 1e-10) {
    stop(
      "Oxygen responseC full-MAP lowest-decile validation failed: cutoff=",
      signif(cutoff, 15), ", n=", sum(data$lowest_objective_decile),
      call. = FALSE
    )
  }
  class_best_index <- unlist(lapply(
    response_curve_class_order,
    function(class_name) {
      index <- which(as.character(data$curve_class) == class_name)
      index[order(data$map_objective[index], data$seed[index])][[1L]]
    }
  ), use.names = FALSE)
  data$class_best_objective <- FALSE
  data$class_best_objective[class_best_index] <- TRUE
  if (sum(data$class_best_objective) != length(response_curve_class_order) ||
      any(table(data$curve_class[data$class_best_objective]) != 1L)) {
    stop("Oxygen responseC requires one objective-minimum fit per response class.",
         call. = FALSE)
  }
  class_best_seed_by_class <- stats::setNames(
    data$seed[data$class_best_objective],
    as.character(data$curve_class[data$class_best_objective])
  )
  class_legend_labels <- stats::setNames(
    paste0(
      response_curve_class_labels[response_curve_class_order],
      " (best seed ",
      class_best_seed_by_class[response_curve_class_order],
      ")"
    ),
    response_curve_class_order
  )

  cluster_hulls <- do.call(rbind, lapply(
    split(data, data$warm_start_region),
    function(d) {
      xy <- unique(d[, c("tSNE1", "tSNE2"), drop = FALSE])
      index <- grDevices::chull(xy$tSNE1, xy$tSNE2)
      hull <- xy[c(index, index[[1L]]), , drop = FALSE]
      center <- colMeans(hull)
      hull$tSNE1 <- center[[1L]] + 1.035 * (hull$tSNE1 - center[[1L]])
      hull$tSNE2 <- center[[2L]] + 1.035 * (hull$tSNE2 - center[[2L]])
      hull$warm_start_region <- as.character(d$warm_start_region[[1L]])
      hull
    }
  ))
  cluster_labels <- do.call(rbind, lapply(
    split(data, data$warm_start_region),
    function(d) {
      x_range <- range(d$tSNE1, finite = TRUE)
      y_range <- range(d$tSNE2, finite = TRUE)
      data.frame(
        tSNE1 = x_range[[2L]] - 0.08 * diff(x_range),
        tSNE2 = y_range[[2L]] - 0.08 * diff(y_range),
        warm_start_region = as.character(d$warm_start_region[[1L]]),
        stringsAsFactors = FALSE
      )
    }
  ))

  class_n <- table(data$curve_class)
  class_axis_labels <- stats::setNames(
    c(
      "Complex\nnonmono.", "Inverted\nU", "Monotone\ninc.", "U-shaped",
      "Approx.\nflat", "Increase\nplateau", "Decrease\nplateau",
      "Monotone\ndec."
    ),
    names(class_n)
  )
  objective_span <- diff(range(data$map_objective, finite = TRUE))
  class_count_labels <- do.call(rbind, lapply(
    response_curve_class_order,
    function(class_name) {
      d <- data[as.character(data$curve_class) == class_name, , drop = FALSE]
      data.frame(
        curve_class = factor(class_name, levels = response_curve_class_order),
        class_x = match(class_name, response_curve_class_order),
        label_y = max(d$map_objective, na.rm = TRUE) + 0.035 * objective_span,
        label = paste0("n=", nrow(d)),
        stringsAsFactors = FALSE
      )
    }
  ))
  half_violin <- do.call(rbind, lapply(
    response_curve_class_order,
    function(class_name) {
      d <- data[as.character(data$curve_class) == class_name, , drop = FALSE]
      class_x <- match(class_name, response_curve_class_order)
      density_fit <- stats::density(
        d$map_objective, n = 256L, adjust = 0.85,
        from = min(d$map_objective), to = max(d$map_objective)
      )
      density_scale <- if (max(density_fit$y) > 0) {
        0.38 * density_fit$y / max(density_fit$y)
      } else {
        rep(0, length(density_fit$y))
      }
      data.frame(
        curve_class = factor(class_name, levels = response_curve_class_order),
        violin_x = c(class_x, class_x - density_scale, class_x),
        map_objective = c(
          density_fit$x[[1L]], density_fit$x, density_fit$x[[length(density_fit$x)]]
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  point_data <- data
  set.seed(20260729L)
  point_data$point_x <- (
    as.integer(point_data$curve_class) +
      stats::runif(nrow(point_data), min = 0.07, max = 0.36)
  )
  set.seed(20260729L)
  p_landscape <- ggplot2::ggplot(
    data, ggplot2::aes(x = tSNE1, y = tSNE2, color = curve_class)
  )
  for (region_name in names(response_region_colors)) {
    p_landscape <- p_landscape + ggplot2::geom_path(
      data = cluster_hulls[
        cluster_hulls$warm_start_region == region_name, , drop = FALSE
      ],
      ggplot2::aes(x = tSNE1, y = tSNE2),
      inherit.aes = FALSE,
      color = response_region_colors[[region_name]],
      linewidth = 0.58, linetype = "22", show.legend = FALSE
    )
  }
  p_landscape <- p_landscape +
    ggplot2::geom_point(size = 1.55, alpha = 0.82) +
    ggplot2::geom_point(
      data = data[data$class_best_objective, , drop = FALSE],
      shape = 21, fill = NA, color = "#111111",
      stroke = 0.85, size = 3.0, show.legend = FALSE
    )
  for (region_name in names(response_region_colors)) {
    region_label_data <- cluster_labels[
        cluster_labels$warm_start_region == region_name, , drop = FALSE
      ]
    p_landscape <- p_landscape +
      ggplot2::geom_label(
        data = region_label_data,
        ggplot2::aes(x = tSNE1, y = tSNE2, label = warm_start_region),
        inherit.aes = FALSE,
        color = response_region_colors[[region_name]], fill = "white",
        family = response_font_family, fontface = "bold",
        size = response_font_sizes[["axis_text"]] / ggplot2::.pt,
        linewidth = 0.32,
        label.padding = grid::unit(0.08, "lines"),
        label.r = grid::unit(0.05, "lines"),
        show.legend = FALSE
      ) +
      shadowtext::geom_shadowtext(
        data = region_label_data,
        ggplot2::aes(x = tSNE1, y = tSNE2, label = warm_start_region),
        inherit.aes = FALSE,
        color = response_region_colors[[region_name]],
        bg.colour = response_region_label_outline_color,
        bg.r = response_region_label_outline_radius,
        family = response_font_family, fontface = "bold",
        size = response_font_sizes[["axis_text"]] / ggplot2::.pt,
        show.legend = FALSE
      )
  }
  p_landscape <- p_landscape +
    ggplot2::scale_color_manual(
      values = response_curve_class_colors,
      labels = class_legend_labels, drop = FALSE
    ) +
    ggplot2::labs(
      x = "Pooled 14-parameter\nt-SNE coordinate 1",
      y = "Pooled 14-parameter\nt-SNE coordinate 2",
      title = paste0(panel_tag, ". Response classes in parameter space"),
      subtitle = paste0(
        "Black rings: class-best fits; dashed outlines: ",
        paste(sort(unique(data$warm_start_region)), collapse = ", ")
      ),
      color = "O2-ploidy response class"
    ) +
    ggplot2::theme_classic(
      base_size = response_font_sizes[["axis_text"]],
      base_family = response_font_family
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = response_font_sizes[["panel_title"]]
      ),
      plot.subtitle = ggplot2::element_text(
        size = response_font_sizes[["subtitle"]], color = "#555555"
      ),
      axis.title = ggplot2::element_text(
        size = response_font_sizes[["axis_title"]]
      ),
      axis.text = ggplot2::element_text(
        size = response_font_sizes[["axis_text"]]
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(
        size = response_font_sizes[["legend_title"]]
      ),
      legend.text = ggplot2::element_text(
        size = response_font_sizes[["legend_text"]]
      ),
      legend.key.height = grid::unit(0.3, "cm"),
      aspect.ratio = 1
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        ncol = 3, override.aes = list(size = 2.5, alpha = 1)
      )
    )

  set.seed(20260729L)
  p_objective <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = cutoff, color = "#333333",
      linetype = "dotted", linewidth = 0.45
    ) +
    ggplot2::geom_polygon(
      data = half_violin,
      ggplot2::aes(
        x = violin_x, y = map_objective,
        group = curve_class, fill = curve_class, color = curve_class
      ),
      alpha = 0.20, linewidth = 0.42,
      show.legend = FALSE
    ) +
    ggplot2::geom_boxplot(
      data = data,
      ggplot2::aes(
        x = as.integer(curve_class) - 0.02,
        y = map_objective,
        group = curve_class, color = curve_class
      ),
      inherit.aes = FALSE,
      width = 0.14, outlier.shape = NA,
      fill = "white", linewidth = 0.42,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(
        x = point_x, y = map_objective, color = curve_class
      ),
      inherit.aes = FALSE,
      size = 0.78, alpha = 0.30,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = class_count_labels,
      ggplot2::aes(x = class_x, y = label_y, label = label),
      inherit.aes = FALSE,
      size = response_font_sizes[["subtitle"]] / ggplot2::.pt,
      family = response_font_family, fontface = "bold",
      color = "#333333", vjust = -0.05
    ) +
    ggplot2::scale_color_manual(
      values = response_curve_class_colors, guide = "none", drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = response_curve_class_colors, guide = "none", drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(response_curve_class_order),
      labels = unname(class_axis_labels[response_curve_class_order]),
      limits = c(0.52, 8.55),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.03, 0.10))
    ) +
    ggplot2::labs(
      x = NULL, y = "Full in vivo MAP objective",
      title = "C. Full MAP fit quality across classes",
      subtitle = "Dotted line: lowest full-MAP decile cutoff; lower is better"
    ) +
    ggplot2::theme_classic(
      base_size = response_font_sizes[["axis_text"]],
      base_family = response_font_family
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = response_font_sizes[["panel_title"]]
      ),
      plot.subtitle = ggplot2::element_text(
        size = response_font_sizes[["subtitle"]], color = "#555555"
      ),
      axis.title = ggplot2::element_text(
        size = response_font_sizes[["axis_title"]]
      ),
      axis.text.y = ggplot2::element_text(
        size = response_font_sizes[["axis_text"]]
      ),
      axis.text.x = ggplot2::element_text(
        size = response_font_sizes[["axis_text"]],
        angle = 48, hjust = 1, vjust = 1
      ),
      aspect.ratio = 1
    )

  c_legend <- cowplot::get_legend(
    p_landscape + ggplot2::theme(legend.position = "bottom")
  )
  c_plot_row <- cowplot::plot_grid(
    p_landscape + ggplot2::theme(legend.position = "none"),
    p_objective,
    nrow = 1, align = "hv", axis = "tblr",
    rel_widths = c(1, 1)
  )
  composite <- cowplot::plot_grid(
    c_plot_row, c_legend,
    ncol = 1, rel_heights = c(1, 0.19)
  )
  figure_paths <- response_save_plot_pair(
    composite,
    file.path(figure_dir, "selection_diagnostic"),
    width = 8.6, height = 5.3
  )
  data_path <- response_write_tsv(
    data, file.path(data_dir, "selection_diagnostic_selection_data.tsv")
  )
  low_class <- as.data.frame(
    table(data$curve_class[data$lowest_objective_decile]),
    stringsAsFactors = FALSE
  )
  names(low_class) <- c("curve_class", "n_lowest_decile")
  class_best <- data[data$class_best_objective, , drop = FALSE]
  class_best <- class_best[
    order(match(as.character(class_best$curve_class), response_curve_class_order)),
    , drop = FALSE
  ]
  validation <- data.frame(
    metric = c(
      "n_seed", "n_class_best", "n_lowest_decile",
      "full_map_lowest_decile_cutoff",
      "map_objective_min", "map_objective_max",
      paste0("class_best_seed_", as.character(class_best$curve_class)),
      paste0(
        "class_best_map_objective_", as.character(class_best$curve_class)
      ),
      paste0("lowest_decile_", low_class$curve_class),
      paste0("warm_start_region_n_", names(expected_region_counts))
    ),
    value = c(
      nrow(data), sum(data$class_best_objective),
      sum(data$lowest_objective_decile), cutoff,
      min(data$map_objective), max(data$map_objective),
      class_best$seed, class_best$map_objective,
      low_class$n_lowest_decile, as.integer(observed_region_counts)
    ),
    stringsAsFactors = FALSE
  )
  validation_path <- response_write_tsv(
    validation, file.path(data_dir, "selection_diagnostic_objective_validation.tsv")
  )
  list(
    data = data, cutoff = cutoff,
    data_paths = c(selection = data_path, validation = validation_path),
    figure_paths = figure_paths
  )
}

# -----------------------------------------------------------------------------
# Assembly, provenance, validation, and orchestration
# -----------------------------------------------------------------------------

response_assemble_main <- function(
    panel_a_png, panel_b_png, panel_c_png, output_png, output_pdf) {
  image_a <- magick::image_read(panel_a_png)
  image_b <- magick::image_read(panel_b_png)
  image_c <- magick::image_read(panel_c_png)

  right_width <- max(
    magick::image_info(image_b)$width,
    magick::image_info(image_c)$width
  )
  image_b <- magick::image_resize(image_b, paste0(right_width, "x"))
  image_c <- magick::image_resize(image_c, paste0(right_width, "x"))
  right_stack <- magick::image_append(c(image_b, image_c), stack = TRUE)

  right_height <- magick::image_info(right_stack)$height
  image_a <- magick::image_resize(image_a, paste0("x", right_height))
  assembled <- magick::image_append(c(image_a, right_stack), stack = FALSE)
  magick::image_write(assembled, path = output_png, format = "png", density = 300)
  magick::image_write(assembled, path = output_pdf, format = "pdf", density = 300)
  c(
    png = normalizePath(output_png, mustWork = TRUE),
    pdf = normalizePath(output_pdf, mustWork = TRUE)
  )
}

response_copy_deliverable <- function(source, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE)
  if (!ok) stop("Failed to copy ", source, " -> ", destination, call. = FALSE)
  normalizePath(destination, mustWork = TRUE)
}

response_validation_rows <- function(
    panel_a, panel_b, panel_c, assembled_paths, script_path) {
  png_info <- function(path) {
    info <- magick::image_info(magick::image_read(path))
    c(width = as.integer(info$width[[1L]]), height = as.integer(info$height[[1L]]))
  }
  dim_a <- png_info(panel_a$figure_paths[["png"]])
  dim_b <- png_info(panel_b$figure_paths[["png"]])
  dim_c <- png_info(panel_c$figure_paths[["png"]])
  dim_main <- png_info(assembled_paths[["png"]])
  panel_a_single_column <- dim_a[["height"]] / dim_a[["width"]] >= 2.5
  abc_integrated <- (
    dim_b[["width"]] == dim_c[["width"]] &&
      dim_a[["height"]] == dim_b[["height"]] + dim_c[["height"]] &&
      dim_main[["height"]] == dim_a[["height"]] &&
      dim_main[["width"]] == dim_a[["width"]] + dim_b[["width"]]
  )

  data.frame(
    check = c(
      "standalone_no_source_calls",
      "response_class_single_column_layout",
      "response_class_n_seed",
      "response_class_n_o2",
      "response_class_n_curve_rows",
      "response_class_class_counts_match",
      "pair_surface_n_pair",
      "pair_surface_n_overlay_rows",
      "pair_surface_n_surface_rows",
      "pair_surface_all_status_ok",
      "pair_surface_max_abs_actual_minus_requested_p_misseg",
      "selection_diagnostic_n_seed",
      "selection_diagnostic_n_class_best",
      "selection_diagnostic_n_lowest_decile",
      "selection_diagnostic_cutoff",
      "selection_diagnostic_class_counts_match",
      "selection_diagnostic_warm_start_region_counts_match",
      "response_workflow_abc_integrated_layout"
    ),
    observed = c(
      !any(grepl(
        "(^|[^[:alnum:]_])source[[:space:]]*\\(",
        readLines(script_path, warn = FALSE)
      )),
      paste(dim_a, collapse = "x"),
      nrow(panel_a$class_by_seed),
      length(unique(panel_a$curves$O2_pct)),
      nrow(panel_a$curves),
      identical(
        stats::setNames(
          as.integer(panel_a$class_counts$n_seed),
          panel_a$class_counts$smooth_curve_class
        )[names(response_expected_class_counts)],
        response_expected_class_counts
      ),
      length(unique(panel_b$surface$pair_id)),
      nrow(panel_b$overlay),
      nrow(panel_b$surface),
      all(panel_b$surface$status == "ok") &&
        all(panel_b$overlay$status == "ok"),
      max(abs(
        panel_b$surface$actual_effective_p_misseg -
          panel_b$surface$effective_p_misseg
      )),
      nrow(panel_c$data),
      sum(panel_c$data$class_best_objective),
      sum(panel_c$data$lowest_objective_decile),
      panel_c$cutoff,
      identical(
        as.integer(table(panel_c$data$curve_class)),
        as.integer(response_expected_class_counts)
      ),
      identical(
        as.integer(table(panel_c$data$warm_start_region)),
        as.integer(c(
          C01 = 84L, C02 = 18L, C03 = 165L,
          C04 = 46L, C05 = 160L, C06 = 27L
        ))
      ),
      paste(
        paste0("A=", paste(dim_a, collapse = "x")),
        paste0("B=", paste(dim_b, collapse = "x")),
        paste0("C=", paste(dim_c, collapse = "x")),
        paste0("main=", paste(dim_main, collapse = "x")),
        sep = "; "
      )
    ),
    expected = c(
      "TRUE", "height/width>=2.5", "500", "201", "100500", "TRUE",
      "6", "1206", "72360", "TRUE", "<=1e-8",
      "500", "8", "50", "14.3592516529645", "TRUE", "TRUE",
      "A left spanning B+C; B upper-right; C lower-right"
    ),
    passed = c(
      !any(grepl(
        "(^|[^[:alnum:]_])source[[:space:]]*\\(",
        readLines(script_path, warn = FALSE)
      )),
      panel_a_single_column,
      nrow(panel_a$class_by_seed) == 500L,
      length(unique(panel_a$curves$O2_pct)) == 201L,
      nrow(panel_a$curves) == 100500L,
      identical(
        stats::setNames(
          as.integer(panel_a$class_counts$n_seed),
          panel_a$class_counts$smooth_curve_class
        )[names(response_expected_class_counts)],
        response_expected_class_counts
      ),
      length(unique(panel_b$surface$pair_id)) == 6L,
      nrow(panel_b$overlay) == 1206L,
      nrow(panel_b$surface) == 72360L,
      all(panel_b$surface$status == "ok") &&
        all(panel_b$overlay$status == "ok"),
      max(abs(
        panel_b$surface$actual_effective_p_misseg -
          panel_b$surface$effective_p_misseg
      )) <= 1e-8,
      nrow(panel_c$data) == 500L,
      sum(panel_c$data$class_best_objective) == 8L,
      sum(panel_c$data$lowest_objective_decile) == 50L,
      abs(panel_c$cutoff - 14.3592516529645) <= 1e-10,
      identical(
        as.integer(table(panel_c$data$curve_class)),
        as.integer(response_expected_class_counts)
      ),
      identical(
        as.integer(table(panel_c$data$warm_start_region)),
        as.integer(c(
          C01 = 84L, C02 = 18L, C03 = 165L,
          C04 = 46L, C05 = 160L, C06 = 27L
        ))
      ),
      abc_integrated
    ),
    stringsAsFactors = FALSE
  )
}

response_main <- function(argv = response_parse_args()) {
  suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(grid)
    library(shadowtext)
  })

  workspace_dir <- normalizePath(
    argv$workspace_dir %||% RESPONSE_SCRIPT_DIR,
    mustWork = FALSE
  )
  dir.create(workspace_dir, recursive = TRUE, showWarnings = FALSE)
  repo_root <- normalizePath(
    Sys.getenv(
      "HYPOXIA_REPO_ROOT",
      unset = file.path(RESPONSE_SCRIPT_DIR, "..", "..", "..", "..", "..", "..")
    ),
    mustWork = TRUE
  )
  if (is.null(argv$joint_root) ||
      !nzchar(trimws(as.character(argv$joint_root)))) {
    stop("--joint_root=PATH is required.")
  }
  joint_root <- normalizePath(argv$joint_root, mustWork = TRUE)
  data_dir <- normalizePath(
    argv$data_dir %||% file.path(workspace_dir, "data", "response_workflow"),
    mustWork = FALSE
  )
  figure_dir <- normalizePath(
    argv$figure_dir %||% file.path(workspace_dir, "figures", "panels"),
    mustWork = FALSE
  )
  output_dir <- normalizePath(
    argv$output_dir %||% file.path(workspace_dir, "figures"),
    mustWork = FALSE
  )
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  workers_default <- max(1L, min(6L, parallel::detectCores(logical = FALSE) - 1L))
  workers <- response_as_int(argv$workers, workers_default)
  rebuild <- response_as_bool(argv$rebuild, FALSE)

  panel_a_input <- file.path(
    joint_root, "cross_validation", "best_fit_parameter_feature",
    "03_dense-grid_monotonicity_classification",
    "monotonicity_classification",
    "dense-grid_monotonicity_classification", "tables",
    "fixed_o2_ploidy_monotonicity_curves.tsv"
  )
  panel_c_input <- file.path(
    joint_root, "landscape_subcluster", "pooled_invivo_invitro",
    "full_data_in_vivo_clustring", "Tables",
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
  )

  message("Building response-class analysis.")
  panel_a <- response_build_panel_a(panel_a_input, data_dir, figure_dir)
  message("Building pair-surface analysis.")
  panel_b <- response_build_panel_b(
    joint_root, data_dir, figure_dir, workers = workers, rebuild = rebuild
  )
  message("Building selection diagnostics.")
  panel_c <- response_build_panel_c(
    panel_c_input, panel_a$class_by_seed, data_dir, figure_dir
  )

  assembled_paths <- response_assemble_main(
    panel_a$figure_paths[["png"]],
    panel_b$figure_paths[["png"]],
    panel_c$figure_paths[["png"]],
    file.path(output_dir, "assembled_oxygen_response.png"),
    file.path(output_dir, "assembled_oxygen_response.pdf")
  )

  delivered <- c(
    revised_assembled_png = assembled_paths[["png"]],
    revised_assembled_pdf = assembled_paths[["pdf"]]
  )

  script_path <- normalizePath(
    file.path(RESPONSE_SCRIPT_DIR, "response_pipeline.R"), mustWork = TRUE
  )
  validation <- response_validation_rows(
    panel_a, panel_b, panel_c, assembled_paths, script_path
  )
  validation_path <- response_write_tsv(
    validation, file.path(data_dir, "response_workflow_validation.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Oxygen response validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", "),
      call. = FALSE
    )
  }

  source_files <- unique(c(
    script_path, panel_a_input, panel_c_input,
    panel_b$selected$fit_summary_path,
    panel_b$selected$parameter_path,
    file.path(panel_b$selected$seed_dir, "fit_config.rds")
  ))
  provenance <- data.frame(
    source_file = source_files,
    exists = file.exists(source_files),
    md5 = vapply(source_files, response_md5, character(1L)),
    role = c(
      "standalone_generator",
      "response_class_500seed_201o2_curves",
      "selection_diagnostic_pooled_full_tsne_and_map_objective",
      rep("pair_surface_selected_joint_fit_summary", nrow(panel_b$selected)),
      rep("pair_surface_selected_joint_parameter_file", nrow(panel_b$selected)),
      rep("pair_surface_selected_joint_fit_config", nrow(panel_b$selected))
    ),
    stringsAsFactors = FALSE
  )
  provenance_path <- response_write_tsv(
    provenance, file.path(data_dir, "response_workflow_source_file_provenance.tsv")
  )

  output_files <- unique(c(
    panel_a$data_paths, panel_a$figure_paths,
    panel_b$data_paths, panel_b$figure_paths,
    panel_c$data_paths, panel_c$figure_paths,
    assembled_paths, delivered, validation_path, provenance_path
  ))
  manifest <- data.frame(
    output_file = output_files,
    exists = file.exists(output_files),
    size_bytes = file.info(output_files)$size,
    md5 = vapply(output_files, response_md5, character(1L)),
    generating_script = script_path,
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(data_dir, "response_workflow_output_manifest.tsv")
  response_write_tsv(manifest, manifest_path)
  message("Oxygen response complete: ", assembled_paths[["png"]])
  invisible(list(
    panel_a = panel_a, panel_b = panel_b, panel_c = panel_c,
    assembled = assembled_paths, validation = validation,
    manifest = manifest_path
  ))
}

if (sys.nframe() == 0L) {
  response_main()
}
