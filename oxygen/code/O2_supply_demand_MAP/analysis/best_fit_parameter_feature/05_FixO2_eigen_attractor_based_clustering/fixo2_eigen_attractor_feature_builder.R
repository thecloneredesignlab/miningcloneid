#!/usr/bin/env Rscript

fixo2ea_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1L))
  )
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else normalizePath(getwd(), mustWork = FALSE)
}

FIXO2EA_CODE_DIR <- fixo2ea_script_dir()
FIXO2EA_ANALYSIS_DIR <- normalizePath(file.path(FIXO2EA_CODE_DIR, ".."), mustWork = FALSE)
source(file.path(FIXO2EA_ANALYSIS_DIR, "util", "parameter_landscape_shared_utils.R"))

fixo2ea_default_result_root <- function() {
  file.path(
    default_oxygen_dir(),
    "results", "analysis", "best_fit_parameter_feature",
    "05_FixO2_eigen_attractor_based_clustering"
  )
}

fixo2ea_default_parameter_source_root <- function() {
  file.path(
    default_oxygen_dir(),
    "results", "analysis", "best_fit_parameter_feature",
    "02_parameter_landscape_clustering"
  )
}

fixo2ea_default_best_csv <- function(dataset, source_root = fixo2ea_default_parameter_source_root()) {
  file.path(source_root, paste0(normalize_dataset(dataset), "_best_params_by_seed.csv"))
}

fixo2ea_default_initial_csv <- function(dataset, source_root = fixo2ea_default_parameter_source_root()) {
  file.path(source_root, paste0(normalize_dataset(dataset), "_deoptim_initial_population.csv"))
}

fixo2ea_load_fixo2_env <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    env <- new.env(parent = globalenv())
    sim_path <- file.path(default_oxygen_dir(), "code", "O2_supply_demand_MAP", "simulation", "fix_o2_simulation.R")
    if (!file.exists(sim_path)) stop("Missing FixO2 simulation script: ", sim_path, call. = FALSE)
    source(sim_path, local = env)
    cache <<- env
    env
  }
})

fixo2ea_o2_grid <- function(o2_values = NULL, o2_min = 0, o2_max = 5, o2_n = 201L) {
  vals <- as_num_vec(o2_values, numeric())
  if (length(vals)) return(sort(unique(vals)))
  o2_min <- as_num(o2_min, 0)
  o2_max <- as_num(o2_max, 5)
  o2_n <- as_int(o2_n, 201L)
  if (!is.finite(o2_min) || !is.finite(o2_max) || o2_max <= o2_min) {
    stop("Invalid O2 range.", call. = FALSE)
  }
  if (!is.finite(o2_n) || is.na(o2_n) || o2_n < 2L) stop("o2_n must be >= 2.", call. = FALSE)
  seq(o2_min, o2_max, length.out = o2_n)
}

fixo2ea_o2_key <- function(x) {
  txt <- formatC(as.numeric(x), format = "f", digits = 3)
  txt <- sub("^-", "m", txt)
  gsub("\\.", "p", txt)
}

fixo2ea_feature_cols <- function(o2_values) {
  paste0("fixo2_eigen_ploidy_o2_", fixo2ea_o2_key(o2_values))
}

fixo2ea_parse_seed_selector <- function(seeds = NULL) {
  if (is.null(seeds) || !length(seeds) || all(is.na(seeds))) return(integer())
  txt <- paste(as.character(seeds), collapse = ",")
  parts <- trimws(strsplit(txt, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  out <- integer()
  for (part in parts) {
    part <- sub("^seed", "", part)
    if (grepl("^[0-9]+:[0-9]+$", part)) {
      bounds <- as.integer(strsplit(part, ":", fixed = TRUE)[[1L]])
      out <- c(out, seq.int(bounds[[1L]], bounds[[2L]]))
    } else {
      val <- suppressWarnings(as.integer(part))
      if (is.finite(val) && !is.na(val)) out <- c(out, val)
    }
  }
  sort(unique(out))
}

fixo2ea_seed_filter <- function(seed_values, seeds = NULL, max_seeds = NA_integer_) {
  seed_values <- sort(unique(as.integer(seed_values)))
  selected <- fixo2ea_parse_seed_selector(seeds)
  if (length(selected)) seed_values <- intersect(seed_values, selected)
  max_seeds <- as_int(max_seeds, NA_integer_)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seed_values <- seed_values[seq_len(min(length(seed_values), max_seeds))]
  }
  seed_values
}

fixo2ea_seed_ids <- function(seed_numbers) {
  paste0("seed", as.integer(seed_numbers))
}

fixo2ea_run_dir_seed_numbers <- function(input_dir, seeds = NULL, max_seeds = NA_integer_) {
  seed_dirs <- list_seed_dirs(input_dir)
  if (!length(seed_dirs)) stop("No seed directories found under: ", input_dir, call. = FALSE)
  fixo2ea_seed_filter(vapply(seed_dirs, seed_from_dir, integer(1L)), seeds = seeds, max_seeds = max_seeds)
}

fixo2ea_read_best_metadata <- function(dataset, best_csv, input_dir, seed_numbers) {
  dataset <- normalize_dataset(dataset)
  if (file.exists(best_csv)) {
    best <- attach_objective(read_csv_plain(best_csv), objective_seed_dir = input_dir)
  } else {
    rows <- lapply(seed_numbers, function(seed) {
      seed_dir <- file.path(input_dir, paste0("seed", as.integer(seed)))
      best_names <- names(read_best_params(seed_dir))
      best_params_row(seed_dir, seed, best_names, dataset = dataset)
    })
    best <- rbind_fill_plain(rows)
  }
  if (!"seed" %in% names(best)) stop("Best parameter table is missing seed column: ", best_csv, call. = FALSE)
  best$seed <- as.integer(best$seed)
  best <- best[best$seed %in% seed_numbers, , drop = FALSE]
  best <- best[match(seed_numbers, best$seed), , drop = FALSE]
  if (any(is.na(best$seed))) stop("Best metadata missing selected seeds.", call. = FALSE)
  best$dataset <- dataset
  best$point_type <- "best"
  best$source_group <- paste0(dataset, "_best")
  best$point_id <- paste(dataset, "best", paste0("seed", best$seed), sep = "__")
  best$source_row_id <- seq_len(nrow(best))
  best$initial_index <- NA_integer_
  best
}

fixo2ea_read_initial_metadata <- function(dataset,
                                          initial_csv,
                                          input_dir,
                                          seed_numbers,
                                          best_names = NULL) {
  dataset <- normalize_dataset(dataset)
  if (file.exists(initial_csv)) {
    initial <- read_csv_plain(initial_csv)
  } else {
    if (is.null(best_names)) {
      first_seed_dir <- file.path(input_dir, paste0("seed", seed_numbers[[1L]]))
      best_names <- names(read_best_params(first_seed_dir))
    }
    rows <- lapply(seed_numbers, function(seed) {
      seed_dir <- file.path(input_dir, paste0("seed", as.integer(seed)))
      metrics <- metric_map(file.path(seed_dir, "fit_summary.tsv"))
      initial_population_natural(seed_dir, seed, metrics, best_names, dataset = dataset)
    })
    initial <- rbind_fill_plain(rows)
  }
  if (!"seed" %in% names(initial)) stop("Initial parameter table is missing seed column: ", initial_csv, call. = FALSE)
  initial$seed <- as.integer(initial$seed)
  initial$source_row_id <- seq_len(nrow(initial))
  initial <- initial[initial$seed %in% seed_numbers, , drop = FALSE]
  initial <- initial[order(match(initial$seed, seed_numbers), initial$source_row_id), , drop = FALSE]
  initial$initial_index <- ave(seq_len(nrow(initial)), initial$seed, FUN = seq_along)
  initial$dataset <- dataset
  initial$point_type <- "initial"
  initial$source_group <- paste0(dataset, "_initial")
  initial$objective <- NA_real_
  initial$point_id <- paste(dataset, "initial", paste0("seed", initial$seed), paste0("init", initial$initial_index), sep = "__")
  initial
}

fixo2ea_extract_cfg_from_fit_result <- function(fit_result) {
  if (!is.list(fit_result)) return(NULL)
  candidates <- list(
    fit_result$cfg,
    fit_result$ctx$invitro_cfg,
    fit_result$ctx$invitro$cfg,
    fit_result$ctx$cfg
  )
  for (cfg in candidates) {
    if (is.list(cfg) && length(cfg)) return(cfg)
  }
  NULL
}

fixo2ea_read_seed_cfg <- function(seed_dir) {
  if (is.na(seed_dir) || !nzchar(seed_dir)) return(NULL)
  fit_config_path <- file.path(seed_dir, "fit_config.rds")
  if (file.exists(fit_config_path)) {
    cfg <- tryCatch(readRDS(fit_config_path), error = function(e) NULL)
    if (is.list(cfg) && length(cfg)) {
      return(list(cfg = cfg, source = fit_config_path))
    }
  }

  fit_result_path <- file.path(seed_dir, "fit_result.rds")
  if (file.exists(fit_result_path)) {
    fit_result <- tryCatch(readRDS(fit_result_path), error = function(e) NULL)
    cfg <- fixo2ea_extract_cfg_from_fit_result(fit_result)
    if (is.list(cfg) && length(cfg)) {
      return(list(cfg = cfg, source = paste0(fit_result_path, "::cfg")))
    }
  }
  NULL
}

fixo2ea_first_available_cfg <- function(input_dir, manifest, first_cfg) {
  seed_dirs <- character()
  if (!is.null(manifest) && "seed_dir" %in% names(manifest)) {
    seed_dirs <- unique(as.character(manifest$seed_dir))
  }
  if (!length(seed_dirs)) {
    seed_dirs <- list_seed_dirs(input_dir)
  }
  seed_dirs <- seed_dirs[!is.na(seed_dirs) & nzchar(seed_dirs)]

  for (seed_dir in seed_dirs) {
    item <- fixo2ea_read_seed_cfg(seed_dir)
    if (!is.null(item)) {
      message("Using FixO2 cfg from ", item$source)
      return(item$cfg)
    }
  }

  warning("No fit_config.rds or cfg-containing fit_result.rds found under ", input_dir, "; falling back to FixO2 default cfg.", call. = FALSE)
  first_cfg(manifest)
}

fixo2ea_dataset_context <- function(input_dir, fixo2_env = fixo2ea_load_fixo2_env()) {
  collect_inputs <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)
  first_cfg <- get("o2pr_first_seed_cfg", envir = fixo2_env, inherits = TRUE)
  source_model <- get("o2ipa_source_model", envir = fixo2_env, inherits = TRUE)
  inputs <- collect_inputs(input_dir, objective_source = "auto")
  cfg <- fixo2ea_first_available_cfg(input_dir, inputs$manifest, first_cfg)
  model_env <- source_model(get("SCRIPT_DIR", envir = fixo2_env, inherits = TRUE))
  list(cfg = cfg, model_env = model_env)
}

fixo2ea_param_columns <- function(df) {
  meta <- c(
    "task_id", "dataset", "point_type", "source_group", "point_id", "seed", "seed_id",
    "source_row_id", "initial_index", "objective", "objective_ploidy",
    "objective_burden", "objective_growth", "objective_flow",
    "deoptim_iter_completed", "pred1000_both_gt44",
    "pred1000_ploidy_ratio_2N", "pred1000_ploidy_ratio_4N"
  )
  cols <- setdiff(names(df), meta)
  cols[vapply(df[cols], function(x) {
    vals <- suppressWarnings(as.numeric(x))
    any(is.finite(vals))
  }, logical(1L))]
}

fixo2ea_meta_columns <- function(df) {
  intersect(
    c(
      "task_id", "dataset", "point_type", "source_group", "point_id", "seed",
      "source_row_id", "initial_index", "objective", "objective_ploidy",
      "objective_burden", "objective_growth", "objective_flow",
      "deoptim_iter_completed"
    ),
    names(df)
  )
}

fixo2ea_compute_one_point <- function(point_row,
                                      param_cols,
                                      o2_values,
                                      context,
                                      fixo2_env,
                                      point_label) {
  pvec <- suppressWarnings(as.numeric(point_row[1L, param_cols, drop = TRUE]))
  names(pvec) <- param_cols
  finite <- is.finite(pvec)
  pvec <- pvec[finite]
  run_params <- get("o2pr_run_params_from_vec", envir = fixo2_env, inherits = TRUE)(pvec, context$cfg)
  dominant_one <- get("fixo2_dominant_attractor_one", envir = fixo2_env, inherits = TRUE)
  rows <- lapply(o2_values, function(o2) {
    dominant_one(point_label, run_params, context$model_env, context$cfg, o2)
  })
  out <- do.call(rbind, rows)
  out$O2_pct <- as.numeric(out$O2_pct)
  out
}

fixo2ea_point_summary <- function(long, o2_values) {
  long <- long[match(o2_values, long$O2_pct), , drop = FALSE]
  vals <- suppressWarnings(as.numeric(long$dominant_mean_ploidy))
  gaps <- suppressWarnings(as.numeric(long$spectral_gap))
  status <- as.character(long$status)
  data.frame(
    n_o2 = length(o2_values),
    n_o2_ok = sum(status == "ok", na.rm = TRUE),
    n_o2_failed = sum(status != "ok" | is.na(status), na.rm = TRUE),
    min_dominant_mean_ploidy = suppressWarnings(min(vals, na.rm = TRUE)),
    max_dominant_mean_ploidy = suppressWarnings(max(vals, na.rm = TRUE)),
    min_spectral_gap = suppressWarnings(min(gaps, na.rm = TRUE)),
    median_spectral_gap = suppressWarnings(stats::median(gaps, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

fixo2ea_fix_infinite_summary <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
  for (col in num_cols) {
    bad <- !is.finite(df[[col]])
    if (any(bad)) df[[col]][bad] <- NA_real_
  }
  df
}

fixo2ea_wide_from_long <- function(long,
                                   metadata,
                                   o2_values,
                                   feature_cols = fixo2ea_feature_cols(o2_values)) {
  rows <- vector("list", nrow(metadata))
  for (i in seq_len(nrow(metadata))) {
    point_id <- metadata$point_id[[i]]
    d <- long[long$point_id == point_id, , drop = FALSE]
    d <- d[match(o2_values, d$O2_pct), , drop = FALSE]
    vals <- suppressWarnings(as.numeric(d$dominant_mean_ploidy))
    names(vals) <- feature_cols
    summary <- fixo2ea_point_summary(d, o2_values)
    rows[[i]] <- cbind(metadata[i, fixo2ea_meta_columns(metadata), drop = FALSE], summary, as.data.frame(as.list(vals), check.names = FALSE))
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  fixo2ea_fix_infinite_summary(out)
}

fixo2ea_status_summary <- function(wide) {
  aggregate(
    cbind(n_points = wide$n_o2, n_o2_ok = wide$n_o2_ok, n_o2_failed = wide$n_o2_failed),
    by = list(dataset = wide$dataset, point_type = wide$point_type),
    FUN = sum,
    na.rm = TRUE
  )
}

fixo2ea_write_table <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, quote = FALSE, row.names = FALSE)
  message("Wrote ", path, " [", nrow(df), " x ", ncol(df), "]")
  invisible(path)
}

fixo2ea_append_csv <- function(df, path, append = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    df,
    file = path,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = !append,
    append = append,
    na = "NA"
  )
}

fixo2ea_build_best_feature_tables <- function(dataset,
                                              input_dir = default_dataset_input_dir(dataset),
                                              best_csv = fixo2ea_default_best_csv(dataset),
                                              result_root = fixo2ea_default_result_root(),
                                              o2_values = fixo2ea_o2_grid(),
                                              seeds = NULL,
                                              max_seeds = NA_integer_,
                                              n_workers = 1L,
                                              force = FALSE) {
  dataset <- normalize_dataset(dataset)
  tables_dir <- file.path(result_root, "FixO2EigenAttractors", "Tables")
  long_path <- file.path(tables_dir, paste0(dataset, "_best_fixo2_eigen_attractor_long.tsv"))
  wide_path <- file.path(tables_dir, paste0(dataset, "_best_fixo2_eigen_attractor_wide.csv"))
  status_path <- file.path(tables_dir, paste0(dataset, "_best_fixo2_eigen_attractor_status_summary.csv"))
  if (!isTRUE(force) && file.exists(wide_path)) {
    message("Reusing existing best FixO2 eigen wide table: ", wide_path)
    return(list(long = long_path, wide = wide_path, status = status_path))
  }

  input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
  seed_numbers <- fixo2ea_run_dir_seed_numbers(input_dir, seeds = seeds, max_seeds = max_seeds)
  seed_ids <- fixo2ea_seed_ids(seed_numbers)
  fixo2_env <- fixo2ea_load_fixo2_env()
  generator <- get("generate_fixo2_attractor_mode_table", envir = fixo2_env, inherits = TRUE)
  message("Building ", dataset, " best FixO2 eigen attractors via generate_fixo2_attractor_mode_table().")
  long <- generator(input_dir, o2_values = o2_values, seed_ids = seed_ids, n_workers = n_workers)
  if (!nrow(long)) stop("No best FixO2 eigen attractor rows were generated for ", dataset, ".", call. = FALSE)
  long$dataset <- dataset
  long$point_type <- "best"
  long$seed <- as.integer(sub("^seed", "", long$seed_id))
  long$source_group <- paste0(dataset, "_best")
  long$point_id <- paste(dataset, "best", long$seed_id, sep = "__")
  long <- long[order(match(long$seed, seed_numbers), long$O2_pct), , drop = FALSE]
  dir.create(dirname(long_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(long, long_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote ", long_path, " [", nrow(long), " x ", ncol(long), "]")

  metadata <- fixo2ea_read_best_metadata(dataset, best_csv, input_dir, seed_numbers)
  wide <- fixo2ea_wide_from_long(long, metadata, o2_values)
  fixo2ea_write_table(wide, wide_path)
  fixo2ea_write_table(fixo2ea_status_summary(wide), status_path)
  list(long = long_path, wide = wide_path, status = status_path)
}

fixo2ea_build_initial_feature_tables <- function(dataset,
                                                 input_dir = default_dataset_input_dir(dataset),
                                                 initial_csv = fixo2ea_default_initial_csv(dataset),
                                                 best_csv = fixo2ea_default_best_csv(dataset),
                                                 result_root = fixo2ea_default_result_root(),
                                                 o2_values = fixo2ea_o2_grid(),
                                                 seeds = NULL,
                                                 max_seeds = NA_integer_,
                                                 chunk_size = 250L,
                                                 n_workers = 1L,
                                                 write_long = FALSE,
                                                 force = FALSE) {
  dataset <- normalize_dataset(dataset)
  tables_dir <- file.path(result_root, "FixO2EigenAttractors", "Tables")
  wide_path <- file.path(tables_dir, paste0(dataset, "_initial_fixo2_eigen_attractor_wide.csv"))
  status_path <- file.path(tables_dir, paste0(dataset, "_initial_fixo2_eigen_attractor_status_summary.csv"))
  long_path <- file.path(tables_dir, paste0(dataset, "_initial_fixo2_eigen_attractor_long.tsv"))
  if (!isTRUE(force) && file.exists(wide_path)) {
    message("Reusing existing initial FixO2 eigen wide table: ", wide_path)
    return(list(long = if (file.exists(long_path)) long_path else NA_character_, wide = wide_path, status = status_path))
  }

  input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
  seed_numbers <- fixo2ea_run_dir_seed_numbers(input_dir, seeds = seeds, max_seeds = max_seeds)
  best_names <- names(read_best_params(file.path(input_dir, paste0("seed", seed_numbers[[1L]]))))
  initial <- fixo2ea_read_initial_metadata(
    dataset = dataset,
    initial_csv = initial_csv,
    input_dir = input_dir,
    seed_numbers = seed_numbers,
    best_names = best_names
  )
  if (!nrow(initial)) stop("No initial rows selected for ", dataset, ".", call. = FALSE)
  param_cols <- intersect(fixo2ea_param_columns(initial), best_names)
  if (!length(param_cols)) stop("No parameter columns found in initial table for ", dataset, ".", call. = FALSE)
  message("Building ", dataset, " full initial FixO2 eigen attractors: ", nrow(initial), " points x ", length(o2_values), " O2 values.")

  fixo2_env <- fixo2ea_load_fixo2_env()
  context <- fixo2ea_dataset_context(input_dir, fixo2_env = fixo2_env)
  feature_cols <- fixo2ea_feature_cols(o2_values)
  chunk_size <- as_int(chunk_size, 250L)
  if (!is.finite(chunk_size) || is.na(chunk_size) || chunk_size < 1L) chunk_size <- 250L
  n_workers <- as_int(n_workers, 1L)
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L

  if (file.exists(wide_path)) unlink(wide_path)
  if (file.exists(long_path)) unlink(long_path)
  chunk_starts <- seq.int(1L, nrow(initial), by = chunk_size)
  append_wide <- FALSE
  append_long <- FALSE
  all_status <- list()

  compute_row <- function(i) {
    row <- initial[i, , drop = FALSE]
    long <- fixo2ea_compute_one_point(
      point_row = row,
      param_cols = param_cols,
      o2_values = o2_values,
      context = context,
      fixo2_env = fixo2_env,
      point_label = row$point_id[[1L]]
    )
    long$dataset <- dataset
    long$point_type <- "initial"
    long$seed <- row$seed[[1L]]
    long$source_group <- row$source_group[[1L]]
    long$point_id <- row$point_id[[1L]]
    long$source_row_id <- row$source_row_id[[1L]]
    long$initial_index <- row$initial_index[[1L]]
    vals <- suppressWarnings(as.numeric(long$dominant_mean_ploidy[match(o2_values, long$O2_pct)]))
    names(vals) <- feature_cols
    summary <- fixo2ea_point_summary(long, o2_values)
    wide <- cbind(row[, fixo2ea_meta_columns(row), drop = FALSE], summary, as.data.frame(as.list(vals), check.names = FALSE))
    list(wide = fixo2ea_fix_infinite_summary(wide), long = long)
  }

  for (chunk_i in seq_along(chunk_starts)) {
    idx <- chunk_starts[[chunk_i]]:min(nrow(initial), chunk_starts[[chunk_i]] + chunk_size - 1L)
    rows <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
      parallel::mclapply(idx, compute_row, mc.cores = min(n_workers, length(idx)))
    } else {
      lapply(idx, compute_row)
    }
    wide_chunk <- do.call(rbind, lapply(rows, `[[`, "wide"))
    fixo2ea_append_csv(wide_chunk, wide_path, append = append_wide)
    append_wide <- TRUE
    all_status[[chunk_i]] <- fixo2ea_status_summary(wide_chunk)
    if (isTRUE(write_long)) {
      long_chunk <- do.call(rbind, lapply(rows, `[[`, "long"))
      dir.create(dirname(long_path), recursive = TRUE, showWarnings = FALSE)
      utils::write.table(
        long_chunk,
        file = long_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = !append_long,
        append = append_long,
        na = "NA"
      )
      append_long <- TRUE
    }
    message("Processed ", dataset, " initial FixO2 eigen chunk ", chunk_i, "/", length(chunk_starts), " [rows ", min(idx), "-", max(idx), "].")
  }
  status <- Reduce(function(a, b) {
    both <- rbind(a, b)
    aggregate(cbind(n_points, n_o2_ok, n_o2_failed) ~ dataset + point_type, both, sum, na.rm = TRUE)
  }, all_status)
  fixo2ea_write_table(status, status_path)
  list(long = if (isTRUE(write_long)) long_path else NA_character_, wide = wide_path, status = status_path)
}

fixo2ea_build_feature_tables <- function(result_root = fixo2ea_default_result_root(),
                                         source_root = fixo2ea_default_parameter_source_root(),
                                         invivo_input = default_dataset_input_dir("invivo"),
                                         invitro_input = default_dataset_input_dir("invitro"),
                                         invivo_best_csv = fixo2ea_default_best_csv("invivo", source_root),
                                         invitro_best_csv = fixo2ea_default_best_csv("invitro", source_root),
                                         invivo_initial_csv = fixo2ea_default_initial_csv("invivo", source_root),
                                         invitro_initial_csv = fixo2ea_default_initial_csv("invitro", source_root),
                                         o2_values = fixo2ea_o2_grid(),
                                         datasets = c("invivo", "invitro"),
                                         point_types = c("best", "initial"),
                                         seeds = NULL,
                                         max_seeds = NA_integer_,
                                         n_workers = 1L,
                                         chunk_size = 250L,
                                         write_initial_long = FALSE,
                                         force = FALSE) {
  result_root <- normalizePath(path.expand(result_root), mustWork = FALSE)
  dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
  datasets <- unique(vapply(as_char_vec(datasets, c("invivo", "invitro")), normalize_dataset, character(1L)))
  point_types <- unique(tolower(as_char_vec(point_types, c("best", "initial"))))
  invalid <- setdiff(point_types, c("best", "initial"))
  if (length(invalid)) stop("point_types must contain only best and/or initial: ", paste(invalid, collapse = ", "), call. = FALSE)

  paths <- list()
  for (dataset in datasets) {
    input_dir <- if (identical(dataset, "invivo")) invivo_input else invitro_input
    best_csv <- if (identical(dataset, "invivo")) invivo_best_csv else invitro_best_csv
    initial_csv <- if (identical(dataset, "invivo")) invivo_initial_csv else invitro_initial_csv
    if ("best" %in% point_types) {
      paths[[paste(dataset, "best", sep = "_")]] <- fixo2ea_build_best_feature_tables(
        dataset = dataset,
        input_dir = input_dir,
        best_csv = best_csv,
        result_root = result_root,
        o2_values = o2_values,
        seeds = seeds,
        max_seeds = max_seeds,
        n_workers = n_workers,
        force = force
      )
    }
    if ("initial" %in% point_types) {
      paths[[paste(dataset, "initial", sep = "_")]] <- fixo2ea_build_initial_feature_tables(
        dataset = dataset,
        input_dir = input_dir,
        initial_csv = initial_csv,
        best_csv = best_csv,
        result_root = result_root,
        o2_values = o2_values,
        seeds = seeds,
        max_seeds = max_seeds,
        n_workers = n_workers,
        chunk_size = chunk_size,
        write_long = write_initial_long,
        force = force
      )
    }
  }
  manifest <- do.call(rbind, lapply(names(paths), function(name) {
    item <- paths[[name]]
    data.frame(
      table = name,
      long_path = as.character(item$long %||% NA_character_),
      wide_path = as.character(item$wide),
      status_path = as.character(item$status),
      stringsAsFactors = FALSE
    )
  }))
  manifest_path <- file.path(result_root, "FixO2EigenAttractors", "Tables", "fixo2_eigen_feature_manifest.csv")
  fixo2ea_write_table(manifest, manifest_path)
  invisible(paths)
}

fixo2ea_task_dir <- function(result_root = fixo2ea_default_result_root()) {
  file.path(result_root, "FixO2EigenAttractors", "HPC")
}

fixo2ea_task_table_path <- function(result_root = fixo2ea_default_result_root()) {
  file.path(fixo2ea_task_dir(result_root), "fixo2_eigen_parameter_tasks.tsv")
}

fixo2ea_task_manifest_path <- function(result_root = fixo2ea_default_result_root()) {
  file.path(fixo2ea_task_dir(result_root), "fixo2_eigen_parameter_task_manifest.tsv")
}

fixo2ea_task_result_root <- function(result_root = fixo2ea_default_result_root()) {
  file.path(fixo2ea_task_dir(result_root), "task_rows")
}

fixo2ea_task_result_path <- function(result_root, task_id) {
  task_id <- as.integer(task_id)
  shard <- sprintf("shard_%04d", as.integer((task_id - 1L) %/% 1000L))
  file.path(fixo2ea_task_result_root(result_root), shard, sprintf("task_%06d.csv", task_id))
}

fixo2ea_write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote ", path, " [", nrow(df), " x ", ncol(df), "]")
  invisible(path)
}

fixo2ea_read_tsv_plain <- function(path) {
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

fixo2ea_dataset_task_rows <- function(dataset,
                                      point_type,
                                      input_dir,
                                      best_csv,
                                      initial_csv,
                                      seed_numbers) {
  dataset <- normalize_dataset(dataset)
  point_type <- match.arg(point_type, c("best", "initial"))
  best_names <- names(read_best_params(file.path(input_dir, paste0("seed", seed_numbers[[1L]]))))
  rows <- if (identical(point_type, "best")) {
    fixo2ea_read_best_metadata(dataset, best_csv, input_dir, seed_numbers)
  } else {
    fixo2ea_read_initial_metadata(dataset, initial_csv, input_dir, seed_numbers, best_names = best_names)
  }
  param_cols <- intersect(fixo2ea_param_columns(rows), best_names)
  if (!length(param_cols)) {
    stop("No parameter columns found for ", dataset, " ", point_type, " task rows.", call. = FALSE)
  }
  keep_cols <- unique(c(fixo2ea_meta_columns(rows), param_cols))
  rows <- rows[, keep_cols, drop = FALSE]
  rows$dataset <- dataset
  rows$point_type <- point_type
  rows
}

fixo2ea_build_task_table <- function(result_root = fixo2ea_default_result_root(),
                                     source_root = fixo2ea_default_parameter_source_root(),
                                     invivo_input = default_dataset_input_dir("invivo"),
                                     invitro_input = default_dataset_input_dir("invitro"),
                                     invivo_best_csv = fixo2ea_default_best_csv("invivo", source_root),
                                     invitro_best_csv = fixo2ea_default_best_csv("invitro", source_root),
                                     invivo_initial_csv = fixo2ea_default_initial_csv("invivo", source_root),
                                     invitro_initial_csv = fixo2ea_default_initial_csv("invitro", source_root),
                                     datasets = c("invivo", "invitro"),
                                     point_types = c("best", "initial"),
                                     seeds = NULL,
                                     max_seeds = NA_integer_,
                                     force = FALSE) {
  result_root <- normalizePath(path.expand(result_root), mustWork = FALSE)
  task_table <- fixo2ea_task_table_path(result_root)
  task_manifest <- fixo2ea_task_manifest_path(result_root)
  if (!isTRUE(force) && file.exists(task_table) && file.exists(task_manifest)) {
    message("Reusing existing FixO2 eigen task table: ", task_table)
    return(list(task_table = task_table, task_manifest = task_manifest))
  }
  datasets <- unique(vapply(as_char_vec(datasets, c("invivo", "invitro")), normalize_dataset, character(1L)))
  point_types <- unique(tolower(as_char_vec(point_types, c("best", "initial"))))
  invalid <- setdiff(point_types, c("best", "initial"))
  if (length(invalid)) stop("point_types must contain only best and/or initial: ", paste(invalid, collapse = ", "), call. = FALSE)

  pieces <- list()
  for (dataset in datasets) {
    input_dir <- normalizePath(path.expand(if (identical(dataset, "invivo")) invivo_input else invitro_input), mustWork = FALSE)
    best_csv <- if (identical(dataset, "invivo")) invivo_best_csv else invitro_best_csv
    initial_csv <- if (identical(dataset, "invivo")) invivo_initial_csv else invitro_initial_csv
    seed_numbers <- fixo2ea_run_dir_seed_numbers(input_dir, seeds = seeds, max_seeds = max_seeds)
    for (point_type in point_types) {
      message("Building task rows for ", dataset, " ", point_type, ".")
      rows <- fixo2ea_dataset_task_rows(
        dataset = dataset,
        point_type = point_type,
        input_dir = input_dir,
        best_csv = best_csv,
        initial_csv = initial_csv,
        seed_numbers = seed_numbers
      )
      pieces[[paste(dataset, point_type, sep = "_")]] <- rows
    }
  }
  tasks <- rbind_fill_plain(pieces)
  tasks$task_id <- seq_len(nrow(tasks))
  tasks <- tasks[, c("task_id", setdiff(names(tasks), "task_id")), drop = FALSE]
  fixo2ea_write_tsv(tasks, task_table)
  manifest <- aggregate(
    task_id ~ dataset + point_type,
    tasks,
    FUN = function(x) paste0(min(x), "-", max(x), " (n=", length(x), ")")
  )
  manifest$total_tasks <- nrow(tasks)
  fixo2ea_write_tsv(manifest, task_manifest)
  list(task_table = task_table, task_manifest = task_manifest)
}

fixo2ea_read_task_rows <- function(task_table, task_ids) {
  tasks <- fixo2ea_read_tsv_plain(task_table)
  task_ids <- as.integer(task_ids)
  rows <- tasks[tasks$task_id %in% task_ids, , drop = FALSE]
  rows <- rows[match(task_ids, rows$task_id), , drop = FALSE]
  if (!nrow(rows) || any(is.na(rows$task_id))) {
    stop("Task table is missing requested task id(s): ", paste(task_ids[is.na(rows$task_id)], collapse = ","), call. = FALSE)
  }
  rows
}

fixo2ea_array_task_ids <- function(array_task_id, points_per_task = 1L, max_task_id = NA_integer_) {
  array_task_id <- as_int(array_task_id, NA_integer_)
  points_per_task <- as_int(points_per_task, 1L)
  max_task_id <- as_int(max_task_id, NA_integer_)
  if (!is.finite(array_task_id) || is.na(array_task_id) || array_task_id < 1L) {
    stop("array_task_id must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(points_per_task) || is.na(points_per_task) || points_per_task < 1L) points_per_task <- 1L
  start <- (array_task_id - 1L) * points_per_task + 1L
  end <- start + points_per_task - 1L
  if (is.finite(max_task_id) && !is.na(max_task_id)) end <- min(end, max_task_id)
  if (end < start) return(integer())
  seq.int(start, end)
}

fixo2ea_compute_task_row <- function(task_row, result_root, o2_values, contexts, fixo2_env) {
  dataset <- normalize_dataset(task_row$dataset[[1L]])
  context <- contexts[[dataset]]
  if (is.null(context)) stop("Missing FixO2 context for dataset: ", dataset, call. = FALSE)
  param_cols <- fixo2ea_param_columns(task_row)
  long <- fixo2ea_compute_one_point(
    point_row = task_row,
    param_cols = param_cols,
    o2_values = o2_values,
    context = context,
    fixo2_env = fixo2_env,
    point_label = task_row$point_id[[1L]]
  )
  vals <- suppressWarnings(as.numeric(long$dominant_mean_ploidy[match(o2_values, long$O2_pct)]))
  names(vals) <- fixo2ea_feature_cols(o2_values)
  summary <- fixo2ea_point_summary(long, o2_values)
  wide <- cbind(
    task_row[, fixo2ea_meta_columns(task_row), drop = FALSE],
    summary,
    as.data.frame(as.list(vals), check.names = FALSE)
  )
  wide <- fixo2ea_fix_infinite_summary(wide)
  out_path <- fixo2ea_task_result_path(result_root, task_row$task_id[[1L]])
  fixo2ea_write_table(wide, out_path)
  out_path
}

fixo2ea_run_task_ids <- function(task_table,
                                 result_root = fixo2ea_default_result_root(),
                                 task_ids,
                                 invivo_input = default_dataset_input_dir("invivo"),
                                 invitro_input = default_dataset_input_dir("invitro"),
                                 o2_values = fixo2ea_o2_grid(),
                                 force = FALSE) {
  result_root <- normalizePath(path.expand(result_root), mustWork = FALSE)
  rows <- fixo2ea_read_task_rows(task_table, task_ids)
  if (!isTRUE(force)) {
    keep <- !file.exists(vapply(rows$task_id, function(id) fixo2ea_task_result_path(result_root, id), character(1L)))
    if (!any(keep)) {
      message("All requested task row outputs already exist: ", paste(task_ids, collapse = ","))
      return(invisible(character()))
    }
    rows <- rows[keep, , drop = FALSE]
  }
  fixo2_env <- fixo2ea_load_fixo2_env()
  datasets <- unique(as.character(rows$dataset))
  contexts <- list()
  if ("invivo" %in% datasets) contexts$invivo <- fixo2ea_dataset_context(invivo_input, fixo2_env = fixo2_env)
  if ("invitro" %in% datasets) contexts$invitro <- fixo2ea_dataset_context(invitro_input, fixo2_env = fixo2_env)
  out <- vapply(seq_len(nrow(rows)), function(i) {
    fixo2ea_compute_task_row(rows[i, , drop = FALSE], result_root, o2_values, contexts, fixo2_env)
  }, character(1L))
  invisible(out)
}

fixo2ea_merge_task_results <- function(result_root = fixo2ea_default_result_root(),
                                       task_table = fixo2ea_task_table_path(result_root),
                                       force = FALSE) {
  result_root <- normalizePath(path.expand(result_root), mustWork = FALSE)
  tasks <- fixo2ea_read_tsv_plain(task_table)
  result_paths <- vapply(tasks$task_id, function(id) fixo2ea_task_result_path(result_root, id), character(1L))
  missing <- result_paths[!file.exists(result_paths)]
  if (length(missing)) {
    stop("Cannot merge FixO2 eigen task rows; missing ", length(missing), " task output(s). First missing: ", missing[[1L]], call. = FALSE)
  }
  tables_dir <- file.path(result_root, "FixO2EigenAttractors", "Tables")
  grouped <- split(seq_along(result_paths), paste(tasks$dataset, tasks$point_type, sep = "_"))
  manifest_rows <- list()
  for (group in names(grouped)) {
    idx <- grouped[[group]]
    first <- TRUE
    out_path <- file.path(tables_dir, paste0(group, "_fixo2_eigen_attractor_wide.csv"))
    if (file.exists(out_path)) {
      if (isTRUE(force)) unlink(out_path) else stop("Output exists; use force=TRUE to overwrite: ", out_path, call. = FALSE)
    }
    status_rows <- list()
    for (path in result_paths[idx]) {
      row <- read_csv_plain(path)
      fixo2ea_append_csv(row, out_path, append = !first)
      first <- FALSE
      status_rows[[length(status_rows) + 1L]] <- row[, intersect(c("dataset", "point_type", "n_o2", "n_o2_ok", "n_o2_failed"), names(row)), drop = FALSE]
    }
    wide <- read_csv_plain(out_path)
    status_path <- file.path(tables_dir, paste0(group, "_fixo2_eigen_attractor_status_summary.csv"))
    fixo2ea_write_table(fixo2ea_status_summary(wide), status_path)
    manifest_rows[[group]] <- data.frame(
      table = group,
      long_path = NA_character_,
      wide_path = out_path,
      status_path = status_path,
      stringsAsFactors = FALSE
    )
    message("Merged ", length(idx), " task row(s) into ", out_path)
  }
  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(tables_dir, "fixo2_eigen_feature_manifest.csv")
  fixo2ea_write_table(manifest, manifest_path)
  invisible(manifest)
}
