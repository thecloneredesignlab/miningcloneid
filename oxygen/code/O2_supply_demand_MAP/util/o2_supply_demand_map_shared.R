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

    scenarios[[i]] <- list(
      harvest = h,
      cohort = cohort,
      dose = dose,
      treat_day = treat_day,
      obs_days = obs_days,
      obs_burden = obs_burden,
      sim_end_day = if (isTRUE(cfg$ploidy_at_harvest)) max(full_days) else max(obs_days),
      harvest_day = max(full_days),
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

  for (i in seq_len(n)) {
    sc <- scenarios[[i]]
    cohort_code[[i]] <- if (identical(sc$cohort, "2N")) 0L else 1L
    dose_vec[[i]] <- as.numeric(sc$dose)
    treat_day_vec[[i]] <- as.numeric(sc$treat_day)
    obs_steps_list[[i]] <- as.integer(round(as.numeric(sc$obs_days) / cfg$DT))
    sim_end_step_vec[[i]] <- as.integer(round(as.numeric(sc$sim_end_day) / cfg$DT))
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
  }

  list(
    cohort_code = cohort_code,
    dose = dose_vec,
    treat_day = treat_day_vec,
    obs_steps = obs_steps_list,
    sim_end_step = sim_end_step_vec,
    obs_burden = obs_burden_list,
    keep_burden = keep_burden_list,
    ploidy_z = ploidy_z_list
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
