#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(
      nzchar,
      vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      }, character(1))
    )
    if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
  }
})

source(file.path(SCRIPT_DIR, "process_fingerprint_utils.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "ploidy_regime_utils.R"), local = TRUE)

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required")
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    if (!grepl("=", arg, fixed = TRUE)) {
      out[[arg]] <- TRUE
    } else {
      key <- sub("=.*$", "", arg)
      val <- sub("^[^=]*=", "", arg)
      out[[key]] <- val
    }
  }
  out
}

repo_root <- function() {
  normalizePath(file.path(SCRIPT_DIR, "..", "..", "..", ".."), mustWork = FALSE)
}

resolve_repo_path <- function(path, root = repo_root(), mustWork = FALSE) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path)) return(path)
  path <- as.character(path[[1]])
  if (startsWith(path, "~/")) {
    path <- file.path(root, substring(path, 3L))
  } else if (identical(path, "~")) {
    path <- root
  } else if (!grepl("^/", path)) {
    path <- file.path(root, path)
  }
  normalizePath(path, mustWork = mustWork)
}

split_csv <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) {
    return(default)
  }
  vals <- trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

normalize_seed_ids <- function(x) {
  if (is.null(x) || !length(x)) return(character())
  o2ipa_norm_seed(x)
}

as_num_vec <- function(x, default) {
  vals <- suppressWarnings(as.numeric(split_csv(x, as.character(default))))
  vals <- vals[is.finite(vals)]
  if (length(vals)) unique(vals) else default
}

as_int_vec <- function(x, default) {
  vals <- suppressWarnings(as.integer(split_csv(x, as.character(default))))
  vals <- vals[!is.na(vals)]
  if (length(vals)) unique(vals) else default
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  tolower(as.character(x[[1]])) %in% c("true", "1", "yes", "y", "on")
}

num_key <- function(x) {
  vapply(as.numeric(x), function(val) {
    if (!is.finite(val)) return(NA_character_)
    format(signif(val, 12), scientific = FALSE, trim = TRUE)
  }, character(1))
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

read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing file: ", path)
    return(data.frame())
  }
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (is_gz) {
      gzip <- Sys.which("gzip")
      if (!nzchar(gzip)) stop("Reading gzip-compressed tables requires gzip on PATH: ", path)
      return(as.data.frame(data.table::fread(cmd = paste(shQuote(gzip), "-cd", shQuote(path)), sep = "\t", data.table = FALSE, showProgress = FALSE)))
    }
    return(as.data.frame(data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE)))
  }
  if (is_gz) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(utils::read.table(
      con,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ))
  }
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    wrote <- FALSE
    if (requireNamespace("data.table", quietly = TRUE)) {
      wrote <- tryCatch({
        data.table::fwrite(x, file = path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
        TRUE
      }, error = function(e) FALSE)
    }
    if (!wrote) {
      con <- gzfile(path, open = "wt")
      on.exit(close(con), add = TRUE)
      utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    }
  } else {
    utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  invisible(path)
}

format_o2_label <- function(x) {
  paste0("O2 = ", format(as.numeric(x), scientific = FALSE, trim = TRUE))
}

initial_condition_from_ploidy <- function(x) {
  paste0("init_", format(as.numeric(x), scientific = FALSE, trim = TRUE), "N")
}

filter_by_values <- function(df, col, values) {
  if (!nrow(df) || !(col %in% names(df))) return(df)
  key <- num_key(df[[col]])
  keep <- key %in% num_key(values)
  df[keep, , drop = FALSE]
}

method_slug <- function(x) {
  x <- tolower(trimws(as.character(x[[1]])))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "analytical" else x
}

method_label <- function(x) {
  x <- method_slug(x)
  if (identical(x, "eigen")) return("Eigen analytical")
  if (identical(x, "expm")) return("Expm analytical")
  paste0(x, " analytical")
}

analytical_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22L), as.integer(cfg$N_MAX %||% 154L))
  G <- o2pr_build_G(model_env, cfg, run_params, O2)
  mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
  M <- G - Matrix::Diagonal(x = mu_all)
  list(M = M, ngrid = ngrid)
}

analytical_init_vector <- function(ngrid, init_N) {
  idx <- which.min(abs(ngrid - init_N))
  v <- numeric(length(ngrid))
  v[[idx]] <- 1
  list(vector = v, used_N = ngrid[[idx]], used_ploidy = ngrid[[idx]] / 22)
}

analytical_normalize_state <- function(x) {
  x <- as.numeric(Re(x))
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  x <- pmax(x, 0)
  s <- sum(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  x / s
}

analytical_state_metrics <- function(state, ngrid, n_unit) {
  w <- analytical_normalize_state(state)
  if (all(is.na(w))) {
    return(data.frame(
      analytical_mean_N = NA_real_,
      analytical_mean_ploidy = NA_real_,
      analytical_sd_ploidy = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  ploidy_grid <- ngrid / n_unit
  mean_N <- sum(ngrid * w, na.rm = TRUE)
  mean_ploidy <- sum(ploidy_grid * w, na.rm = TRUE)
  second_ploidy <- sum(ploidy_grid^2 * w, na.rm = TRUE)
  data.frame(
    analytical_mean_N = mean_N,
    analytical_mean_ploidy = mean_ploidy,
    analytical_sd_ploidy = sqrt(max(0, second_ploidy - mean_ploidy^2)),
    stringsAsFactors = FALSE
  )
}

analytical_eigen_states <- function(M, init, time_grid) {
  Mdense <- as.matrix(M)
  eig <- tryCatch(eigen(Mdense, only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) stop("eigen decomposition failed")
  coef <- tryCatch(solve(eig$vectors, init), error = function(e) NULL)
  if (is.null(coef)) stop("eigen coefficient solve failed")
  lambda_ref <- max(Re(eig$values), na.rm = TRUE)
  lapply(time_grid, function(tt) {
    weights <- exp((eig$values - lambda_ref) * tt) * coef
    analytical_normalize_state(eig$vectors %*% weights)
  })
}

analytical_expm_states <- function(M, init, time_grid) {
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
      scale <- max(abs(x), na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) {
        x[] <- NA_real_
      } else {
        x <- x / scale
      }
      t_now <- target
    }
    states[[i]] <- analytical_normalize_state(x)
  }
  states
}

generate_seed_analytical_trajectories <- function(seed_id, inputs, param_mat, model_env,
                                                  time_points, o2_values, initial_ploidy_values,
                                                  analytical_methods) {
  manifest <- inputs$manifest[inputs$manifest$seed_id == seed_id, , drop = FALSE]
  if (!nrow(manifest) || !seed_id %in% rownames(param_mat)) return(data.frame())
  cfg <- o2pr_first_seed_cfg(manifest)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
  names(pvec) <- colnames(param_mat)
  run_params <- o2pr_run_params_from_vec(pvec, cfg)
  init_specs <- data.frame(
    initial_condition = paste0("init_", format(initial_ploidy_values, scientific = FALSE, trim = TRUE), "N"),
    initial_ploidy = initial_ploidy_values,
    requested_N = initial_ploidy_values * n_unit,
    stringsAsFactors = FALSE
  )

  rows <- list()
  k <- 0L
  for (O2 in o2_values) {
    fm <- analytical_fixed_matrix(model_env, cfg, run_params, O2)
    for (j in seq_len(nrow(init_specs))) {
      init <- analytical_init_vector(fm$ngrid, init_specs$requested_N[[j]])
      states_by_method <- list()
      if ("eigen" %in% analytical_methods) {
        states_by_method$eigen <- tryCatch(analytical_eigen_states(fm$M, init$vector, time_points), error = function(e) NULL)
      }
      if ("expm" %in% analytical_methods) {
        states_by_method$expm <- tryCatch(analytical_expm_states(fm$M, init$vector, time_points), error = function(e) NULL)
      }
      for (method in names(states_by_method)) {
        states <- states_by_method[[method]]
        if (is.null(states)) next
        for (i in seq_along(time_points)) {
          met <- analytical_state_metrics(states[[i]], fm$ngrid, n_unit)
          k <- k + 1L
          rows[[k]] <- data.frame(
            seed_id = seed_id,
            O2_pct = O2,
            O2_key = num_key(O2),
            initial_condition = init_specs$initial_condition[[j]],
            initial_ploidy = init_specs$initial_ploidy[[j]],
            requested_initial_N = init_specs$requested_N[[j]],
            used_initial_N = init$used_N,
            day = time_points[[i]],
            day_key = num_key(time_points[[i]]),
            analytical_method = method,
            analytical_method_label = method_label(method),
            met,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame()
  out
}

generate_analytical_trajectories <- function(run_dir, time_points, o2_values, initial_ploidy_values,
                                             analytical_methods, n_workers = 1L, seed_ids = NULL) {
  if (!dir.exists(run_dir)) stop("run_dir does not exist and is required to generate analytical trajectories: ", run_dir)
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  seeds <- if (is.null(seed_ids) || !length(seed_ids)) rownames(param_mat) else intersect(seed_ids, rownames(param_mat))
  if (!length(seeds)) stop("No seed parameters were found for analytical trajectory generation.")
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, length(seeds))
  message("Generating analytical trajectories from fitted parameters: ", length(seeds), " seeds, methods=", paste(analytical_methods, collapse = ","), ", workers=", n_workers)
  worker <- function(seed_id) {
    generate_seed_analytical_trajectories(
      seed_id = seed_id,
      inputs = inputs,
      param_mat = param_mat,
      model_env = model_env,
      time_points = time_points,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods
    )
  }
  rows <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    parallel::mclapply(seeds, worker, mc.cores = n_workers)
  } else {
    lapply(seeds, worker)
  }
  out <- do.call(rbind, rows[vapply(rows, nrow, integer(1)) > 0L])
  if (is.null(out)) out <- data.frame()
  out
}

filter_analytical_trajectories <- function(analytical, time_points, o2_values, initial_ploidy_values, analytical_methods, seed_ids = NULL) {
  if (!nrow(analytical)) return(analytical)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "analytical_mean_N", "analytical_mean_ploidy", "analytical_sd_ploidy"),
    names(analytical)
  )
  for (col in numeric_cols) analytical[[col]] <- suppressWarnings(as.numeric(analytical[[col]]))
  analytical$analytical_method <- vapply(analytical$analytical_method, method_slug, character(1))
  out <- filter_by_values(analytical, "day", time_points)
  out <- filter_by_values(out, "O2_pct", o2_values)
  out <- filter_by_values(out, "initial_ploidy", initial_ploidy_values)
  out <- out[out$analytical_method %in% analytical_methods, , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids)) out <- out[out$seed_id %in% seed_ids, , drop = FALSE]
  out$O2_key <- num_key(out$O2_pct)
  out$day_key <- num_key(out$day)
  out
}

expected_analytical_seed_ids <- function(run_dir, seed_ids = NULL) {
  seed_ids <- normalize_seed_ids(seed_ids)
  if (length(seed_ids)) return(seed_ids)
  if (is.null(run_dir) || !nzchar(run_dir) || !dir.exists(run_dir)) return(character())
  seeds <- tryCatch(o2ipa_discover_seeds(run_dir)$seed_id, error = function(e) character())
  normalize_seed_ids(seeds)
}

analytical_cache_missing_keys <- function(analytical, time_points, o2_values, initial_ploidy_values,
                                          analytical_methods, expected_seed_ids = character()) {
  if (!nrow(analytical)) return("all")
  required_cols <- c("seed_id", "O2_pct", "initial_ploidy", "day", "analytical_method", "analytical_mean_ploidy")
  missing_cols <- setdiff(required_cols, names(analytical))
  if (length(missing_cols)) return(paste0("missing column: ", paste(missing_cols, collapse = ",")))

  seeds <- normalize_seed_ids(expected_seed_ids)
  if (!length(seeds)) seeds <- normalize_seed_ids(unique(analytical$seed_id))
  if (!length(seeds)) return("no seeds")

  cache <- data.frame(
    seed_id = normalize_seed_ids(analytical$seed_id),
    O2_key = num_key(analytical$O2_pct),
    initial_ploidy_key = num_key(analytical$initial_ploidy),
    day_key = num_key(analytical$day),
    analytical_method = vapply(analytical$analytical_method, method_slug, character(1)),
    analytical_mean_ploidy = suppressWarnings(as.numeric(analytical$analytical_mean_ploidy)),
    stringsAsFactors = FALSE
  )
  cache <- cache[is.finite(cache$analytical_mean_ploidy), , drop = FALSE]
  available_keys <- unique(paste(cache$seed_id, cache$O2_key, cache$initial_ploidy_key, cache$day_key, cache$analytical_method, sep = "\r"))
  required <- expand.grid(
    seed_id = seeds,
    O2_key = num_key(o2_values),
    initial_ploidy_key = num_key(initial_ploidy_values),
    day_key = num_key(time_points),
    analytical_method = analytical_methods,
    stringsAsFactors = FALSE
  )
  required_keys <- paste(required$seed_id, required$O2_key, required$initial_ploidy_key, required$day_key, required$analytical_method, sep = "\r")
  missing <- setdiff(required_keys, available_keys)
  if (!length(missing)) return(character())
  head(missing, 5L)
}

read_seed_objectives <- function(analysis_dir, fit_dir = NULL) {
  attractor_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_spectral_gap_by_seed.tsv"
  )
  if (file.exists(attractor_path)) {
    tab <- read_tsv(attractor_path)
    cols <- intersect(c("seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective"), names(tab))
    if (all(c("seed_id", "objective") %in% cols)) {
      tab <- tab[, cols, drop = FALSE]
      tab <- tab[!duplicated(tab$seed_id), , drop = FALSE]
      tab$objective <- as.numeric(tab$objective)
      return(tab)
    }
  }

  if (is.null(fit_dir) || !nzchar(fit_dir) || !dir.exists(fit_dir)) {
    warning("Seed objective values were not found in the analysis directory, and fit_dir is unavailable.")
    return(data.frame(seed_id = character(), objective = numeric(), stringsAsFactors = FALSE))
  }

  seed_dirs <- list.dirs(fit_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("seed[0-9]+$", basename(seed_dirs))]
  rows <- lapply(seed_dirs, function(seed_dir) {
    candidates <- c(
      file.path(seed_dir, "best_params.tsv"),
      file.path(seed_dir, "best_parameters.tsv"),
      file.path(seed_dir, "parameter_table_input.csv")
    )
    hits <- candidates[file.exists(candidates)]
    path <- if (length(hits)) hits[[1]] else NA_character_
    if (is.na(path)) return(NULL)
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    tab <- tryCatch(
      utils::read.table(path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) data.frame()
    )
    objective_cols <- intersect(
      c("objective", "optimizer_local_objective", "optimizer_deoptim_objective", "loss", "best_objective"),
      names(tab)
    )
    if (!length(objective_cols) || !nrow(tab)) return(NULL)
    data.frame(
      seed_id = basename(seed_dir),
      objective = suppressWarnings(as.numeric(tab[[objective_cols[[1]]]][[1]])),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) out <- data.frame(seed_id = character(), objective = numeric(), stringsAsFactors = FALSE)
  out
}

task_table_path <- function(simulation_dir, simulation_mode) {
  file.path(simulation_dir, simulation_mode, "task_list.tsv")
}

build_task_table_from_paths <- function(simulation_dir, simulation_mode, seed_ids, o2_values, initial_ploidy_values, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    for (O2 in o2_values) {
      for (initial_ploidy in initial_ploidy_values) {
        for (simulation_id in simulation_ids) {
          k <- k + 1L
          output_dir <- file.path(
            simulation_dir,
            simulation_mode,
            paste0("O2_", num_path_tag(O2)),
            paste0("ploidy", num_path_tag(initial_ploidy)),
            seed_id,
            paste0("sim", simulation_id)
          )
          rows[[k]] <- data.frame(
            task_id = k,
            seed_id = seed_id,
            fixed_o2_pct = O2,
            initial_ploidy = initial_ploidy,
            initial_condition = initial_condition_from_ploidy(initial_ploidy),
            simulation_id = simulation_id,
            output_dir = output_dir,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

read_simulation_tasks <- function(simulation_dir, simulation_mode, analytical, o2_values, initial_ploidy_values, simulation_ids) {
  path <- task_table_path(simulation_dir, simulation_mode)
  if (file.exists(path)) {
    tasks <- read_tsv(path)
    if (!"fixed_o2_pct" %in% names(tasks) && "O2_pct" %in% names(tasks)) {
      names(tasks)[names(tasks) == "O2_pct"] <- "fixed_o2_pct"
    }
    required <- c("seed_id", "fixed_o2_pct", "initial_ploidy", "simulation_id", "output_dir")
    missing <- setdiff(required, names(tasks))
    if (length(missing)) stop("Simulation task table is missing column(s): ", paste(missing, collapse = ", "))
    if (!"initial_condition" %in% names(tasks)) {
      tasks$initial_condition <- initial_condition_from_ploidy(tasks$initial_ploidy)
    }
  } else {
    seed_ids <- sort(unique(analytical$seed_id))
    tasks <- build_task_table_from_paths(
      simulation_dir = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      simulation_ids = simulation_ids
    )
  }
  tasks$fixed_o2_pct <- as.numeric(tasks$fixed_o2_pct)
  tasks$initial_ploidy <- as.numeric(tasks$initial_ploidy)
  tasks$simulation_id <- as.integer(tasks$simulation_id)
  expected_output_dir <- file.path(
    simulation_dir,
    simulation_mode,
    paste0("O2_", vapply(tasks$fixed_o2_pct, num_path_tag, character(1))),
    paste0("ploidy", vapply(tasks$initial_ploidy, num_path_tag, character(1))),
    tasks$seed_id,
    paste0("sim", tasks$simulation_id)
  )
  use_expected <- !dir.exists(tasks$output_dir) & dir.exists(expected_output_dir)
  tasks$output_dir[use_expected] <- expected_output_dir[use_expected]
  tasks <- filter_by_values(tasks, "fixed_o2_pct", o2_values)
  tasks <- filter_by_values(tasks, "initial_ploidy", initial_ploidy_values)
  tasks <- tasks[tasks$simulation_id %in% simulation_ids, , drop = FALSE]
  tasks <- tasks[tasks$seed_id %in% unique(analytical$seed_id), , drop = FALSE]
  tasks$state_file <- file.path(tasks$output_dir, "state_trajectory.tsv.gz")
  tasks
}

read_state_metrics_awk <- function(path, time_points = NULL) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    return(data.frame())
  }
  gzip <- Sys.which("gzip")
  awk <- Sys.which("awk")
  if (!nzchar(gzip) || !nzchar(awk)) {
    stop("Reading state trajectories requires gzip and awk on PATH.")
  }
  all_days <- is.null(time_points) || !length(time_points)
  days_arg <- if (all_days) "" else paste(num_key(time_points), collapse = ",")
  awk_script <- paste(
    'BEGIN {',
    '  keep_all = all_days + 0;',
    '  if (!keep_all) {',
    '    split(days, d, ",");',
    '    for (i in d) keep[sprintf("%.10g", d[i] + 0)] = 1;',
    '  }',
    '  OFS = "\t";',
    '}',
    'NR == 1 {',
    '  for (i = 1; i <= NF; i++) idx[$i] = i;',
    '  next;',
    '}',
    '{',
    '  day = sprintf("%.10g", $(idx["day"]) + 0);',
    '  if ((keep_all || (day in keep)) && $(idx["status"]) == "live") {',
    '    c = $(idx["cell_count"]) + 0;',
    '    n = $(idx["N"]) + 0;',
    '    p = $(idx["ploidy"]) + 0;',
    '    sumc[day] += c;',
    '    sumn[day] += n * c;',
    '    sump[day] += p * c;',
    '    sump2[day] += p * p * c;',
    '  }',
    '}',
    'END {',
    '  for (day in sumc) {',
    '    if (sumc[day] > 0) {',
    '      meanp = sump[day] / sumc[day];',
    '      varp = sump2[day] / sumc[day] - meanp * meanp;',
    '      if (varp < 0 && varp > -1e-9) varp = 0;',
    '      print day, sumn[day] / sumc[day], meanp, sqrt(varp), sumc[day];',
    '    }',
    '  }',
    '}',
    sep = "\n"
  )
  cmd <- paste(
    shQuote(gzip), "-cd", shQuote(path), "|",
    shQuote(awk),
    "-v", shQuote(paste0("days=", days_arg)),
    "-v", shQuote(paste0("all_days=", if (all_days) "1" else "0")),
    shQuote(awk_script)
  )
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
  if (!length(out)) return(data.frame())
  con <- textConnection(out)
  on.exit(close(con), add = TRUE)
  tab <- utils::read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(tab) <- c("day", "simulation_mean_N", "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells")
  tab$day <- as.numeric(tab$day)
  tab <- tab[order(tab$day), , drop = FALSE]
  tab
}

read_simulation_metrics_chunk <- function(tasks, idx, time_points = NULL, progress_every = 100L, worker_label = NULL, total_tasks = nrow(tasks)) {
  rows <- vector("list", length(idx))
  missing <- 0L
  for (j in seq_along(idx)) {
    i <- idx[[j]]
    if (progress_every > 0L && (j == 1L || j %% progress_every == 0L || j == length(idx))) {
      label <- if (is.null(worker_label) || !nzchar(worker_label)) "" else paste0(worker_label, ": ")
      message(label, "Reading simulation state metrics: ", j, "/", length(idx), " (task ", i, "/", total_tasks, ")")
    }
    task <- tasks[i, , drop = FALSE]
    metric <- read_state_metrics_awk(task$state_file[[1]], time_points)
    if (!nrow(metric)) {
      missing <- missing + 1L
      next
    }
    metric$seed_id <- task$seed_id[[1]]
    metric$O2_pct <- task$fixed_o2_pct[[1]]
    metric$O2_key <- num_key(metric$O2_pct)
    metric$initial_ploidy <- task$initial_ploidy[[1]]
    metric$initial_condition <- task$initial_condition[[1]]
    metric$simulation_id <- task$simulation_id[[1]]
    metric$day_key <- num_key(metric$day)
    rows[[j]] <- metric
  }
  if (missing) {
    label <- if (is.null(worker_label) || !nzchar(worker_label)) "" else paste0(worker_label, ": ")
    warning(label, "Missing or unreadable simulation state files: ", missing)
  }
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) out <- data.frame()
  out
}

read_simulation_metrics <- function(tasks, time_points = NULL, progress_every = 100L, n_workers = 1L) {
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, nrow(tasks))
  if (n_workers <= 1L || !identical(.Platform$OS.type, "unix")) {
    return(read_simulation_metrics_chunk(
      tasks = tasks,
      idx = seq_len(nrow(tasks)),
      time_points = time_points,
      progress_every = progress_every,
      total_tasks = nrow(tasks)
    ))
  }

  message("Reading simulation state metrics with ", n_workers, " workers: ", nrow(tasks), " tasks")
  idx <- seq_len(nrow(tasks))
  chunks <- split(idx, ((idx - 1L) %% n_workers) + 1L)
  res <- parallel::mclapply(
    seq_along(chunks),
    function(worker_id) {
      read_simulation_metrics_chunk(
        tasks = tasks,
        idx = chunks[[worker_id]],
        time_points = time_points,
        progress_every = progress_every,
        worker_label = sprintf("worker %02d", worker_id),
        total_tasks = nrow(tasks)
      )
    },
    mc.cores = n_workers
  )
  out <- do.call(rbind, res[vapply(res, nrow, integer(1)) > 0L])
  if (is.null(out)) out <- data.frame()
  out
}

filter_simulation_metrics <- function(sim_metrics, time_points, o2_values, initial_ploidy_values, simulation_ids, seed_ids = NULL) {
  if (!nrow(sim_metrics)) return(sim_metrics)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "simulation_id", "simulation_mean_N",
      "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells"),
    names(sim_metrics)
  )
  for (col in numeric_cols) sim_metrics[[col]] <- suppressWarnings(as.numeric(sim_metrics[[col]]))
  out <- filter_by_values(sim_metrics, "day", time_points)
  out <- filter_by_values(out, "O2_pct", o2_values)
  out <- filter_by_values(out, "initial_ploidy", initial_ploidy_values)
  if ("simulation_id" %in% names(out)) out <- out[out$simulation_id %in% simulation_ids, , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids) && "seed_id" %in% names(out)) {
    out <- out[out$seed_id %in% seed_ids, , drop = FALSE]
  }
  out$O2_key <- num_key(out$O2_pct)
  out$day_key <- num_key(out$day)
  out
}

aggregate_replicates <- function(sim_metrics) {
  if (!nrow(sim_metrics)) return(sim_metrics)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "simulation_id", "simulation_mean_N",
      "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells"),
    names(sim_metrics)
  )
  for (col in numeric_cols) sim_metrics[[col]] <- suppressWarnings(as.numeric(sim_metrics[[col]]))
  keys <- c("seed_id", "O2_pct", "O2_key", "initial_condition", "initial_ploidy", "day", "day_key")
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(sim_metrics)
    out <- dt[, .(
      simulation_n = data.table::uniqueN(simulation_id),
      simulation_mean_N = mean(simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(simulation_live_cells, na.rm = TRUE)
    ), by = keys]
    return(as.data.frame(out))
  }
  split_key <- interaction(sim_metrics[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sim_metrics, split_key), function(x) {
    data.frame(
      x[1, keys, drop = FALSE],
      simulation_n = length(unique(x$simulation_id)),
      simulation_mean_N = mean(x$simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(x$simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(x$simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(x$simulation_live_cells, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

merge_scatter_data <- function(analytical, sim_summary, objectives) {
  by_cols <- c("seed_id", "O2_key", "initial_condition", "day_key")
  dat <- merge(
    analytical,
    sim_summary,
    by = by_cols,
    all = FALSE,
    suffixes = c("_analytical", "_simulation")
  )
  if ("O2_pct_analytical" %in% names(dat)) dat$O2_pct <- dat$O2_pct_analytical
  if ("day_analytical" %in% names(dat)) dat$day <- dat$day_analytical
  dat <- merge(dat, objectives, by = "seed_id", all.x = TRUE, suffixes = c("", "_objective"))
  if ("mode_label_objective" %in% names(dat) && "mode_label" %in% names(dat)) {
    fill <- !nzchar(as.character(dat$mode_label)) | is.na(dat$mode_label)
    dat$mode_label[fill] <- dat$mode_label_objective[fill]
  }
  dat$O2_factor <- factor(format_o2_label(dat$O2_pct), levels = format_o2_label(sort(unique(dat$O2_pct))))
  dat$day_factor <- factor(paste0("Day ", format(as.numeric(dat$day), scientific = FALSE, trim = TRUE)),
                           levels = paste0("Day ", format(sort(unique(as.numeric(dat$day))), scientific = FALSE, trim = TRUE)))
  dat$initial_condition <- factor(dat$initial_condition, levels = sort(unique(dat$initial_condition)))
  dat$objective <- as.numeric(dat$objective)
  dat[is.finite(dat$analytical_mean_ploidy) & is.finite(dat$simulation_mean_ploidy), , drop = FALSE]
}

plot_limits <- function(dat) {
  vals <- c(dat$analytical_mean_ploidy, dat$simulation_mean_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(0, 1))
  rng <- range(vals)
  pad <- diff(rng) * 0.04
  if (!is.finite(pad) || pad <= 0) pad <- 0.05
  rng + c(-pad, pad)
}

objective_aesthetic <- function(dat, transform = "identity") {
  objective <- as.numeric(dat$objective)
  label <- "Objective"
  if (identical(transform, "log10")) {
    objective <- ifelse(is.finite(objective) & objective > 0, log10(objective), NA_real_)
    label <- "log10(objective)"
  }
  list(value = objective, label = label)
}

base_scatter <- function(dat, limits, point_size = 0.9, alpha = 0.55) {
  ggplot2::ggplot(dat, ggplot2::aes(x = analytical_mean_ploidy, y = simulation_mean_ploidy)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey45", linetype = 2, linewidth = 0.35) +
    ggplot2::coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    ggplot2::labs(
      x = "Analytical solution mean ploidy",
      y = "Simulation-inferred mean ploidy"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
}

save_plot <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", device = "pdf")
  invisible(path)
}

plot_time_facets_color_o2 <- function(dat, path, limits) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = O2_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::labs(fill = "Fixed O2", shape = "Initial condition")
  save_plot(p, path, width = 13, height = 7)
}

plot_time_facets_color_objective <- function(dat, path, limits, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition")
  save_plot(p, path, width = 13, height = 7)
}

plot_o2_facets_color_objective <- function(dat, path, limits, title = NULL, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition", title = title)
  save_plot(p, path, width = 11, height = 7)
}

make_scatter_outputs <- function(dat, out_dir, objective_transform = "identity", analytical_method = "analytical") {
  fig_dir <- file.path(out_dir, "simulation", "scatters", method_slug(analytical_method))
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  limits <- plot_limits(dat)

  plot_time_facets_color_o2(
    dat,
    file.path(fig_dir, "scatter_analytical_vs_simulation_by_time_color_o2.pdf"),
    limits = limits
  )
  plot_time_facets_color_objective(
    dat,
    file.path(fig_dir, "scatter_analytical_vs_simulation_by_time_color_objective.pdf"),
    limits = limits,
    objective_transform = objective_transform
  )
  plot_o2_facets_color_objective(
    dat,
    file.path(fig_dir, "scatter_analytical_vs_simulation_by_o2_color_objective_all_times.pdf"),
    limits = limits,
    title = "All selected time points",
    objective_transform = objective_transform
  )

  for (day in sort(unique(as.numeric(dat$day)))) {
    day_dat <- dat[abs(as.numeric(dat$day) - day) < 1e-9, , drop = FALSE]
    day_label <- format(day, scientific = FALSE, trim = TRUE)
    plot_o2_facets_color_objective(
      day_dat,
      file.path(fig_dir, paste0("scatter_analytical_vs_simulation_by_o2_color_objective_day", day_label, ".pdf")),
      limits = limits,
      title = paste0("Day ", day_label),
      objective_transform = objective_transform
    )
  }
  invisible(fig_dir)
}

main <- function(argv = parse_args()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package is required for plotting.")

  root <- resolve_repo_path(argv$repo_root %||% "~", root = repo_root(), mustWork = TRUE)
  simulation_dir <- resolve_repo_path(argv$simulation_dir %||% "~/oxygen/results/O2_fixed_simulation", root = root, mustWork = FALSE)
  analysis_dir <- resolve_repo_path(argv$analysis_dir %||% "~/oxygen/results/analysis/FixO2_invivo_500seed", root = root, mustWork = TRUE)
  fit_dir <- resolve_repo_path(argv$fit_dir %||% "~/oxygen/results/fit_invitro_O2_buffering_500seed", root = root, mustWork = FALSE)
  run_dir <- resolve_repo_path(argv$run_dir %||% fit_dir, root = root, mustWork = FALSE)
  out_dir <- resolve_repo_path(argv$out_dir %||% "~/oxygen/results/analysis/FixO2_invivo_500seed", root = root, mustWork = FALSE)
  simulation_mode <- argv$simulation_mode %||% "invivo"
  time_points <- sort(as_num_vec(argv$time_points, c(25, 50, 100, 200, 300, 500, 700, 1000)))
  o2_values <- sort(as_num_vec(argv$o2_values, c(0, 0.1, 0.5, 1, 2, 5)))
  initial_ploidy_values <- sort(as_num_vec(argv$initial_ploidy_values, c(2, 4)))
  simulation_ids <- sort(as_int_vec(argv$simulation_ids, c(1L, 2L, 3L)))
  objective_transform <- argv$objective_transform %||% "identity"
  if (!objective_transform %in% c("identity", "log10")) stop("--objective_transform must be identity or log10.")
  recompute <- as_bool(argv$recompute, FALSE)
  cache_all_times <- as_bool(argv$cache_all_times, TRUE)
  n_workers <- suppressWarnings(as.integer(argv$n_workers %||% Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  analytical_methods <- split_csv(argv$analytical_methods %||% "eigen,expm", default = c("eigen", "expm"))
  analytical_methods <- unique(vapply(analytical_methods, method_slug, character(1)))
  seed_ids <- normalize_seed_ids(split_csv(argv$seed_ids %||% "", default = character()))
  recompute_analytical <- as_bool(argv$recompute_analytical, recompute)

  table_dir <- file.path(out_dir, "simulation", "scatters", "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  analytical_cache_path <- resolve_repo_path(argv$analytical_cache_table %||% file.path(table_dir, "scatter_generated_analytical_trajectories.tsv.gz"), root = root, mustWork = FALSE)
  all_time_sim_metric_path <- resolve_repo_path(argv$all_time_simulation_metric_table %||% file.path(table_dir, "scatter_simulation_all_time_metrics_by_replicate.tsv.gz"), root = root, mustWork = FALSE)
  sim_metric_path <- resolve_repo_path(argv$simulation_metric_table %||% file.path(table_dir, "scatter_simulation_selected_time_metrics_by_replicate.tsv"), root = root, mustWork = FALSE)
  sim_summary_path <- resolve_repo_path(argv$simulation_summary_table %||% file.path(table_dir, "scatter_simulation_selected_time_metrics.tsv"), root = root, mustWork = FALSE)
  scatter_data_path <- resolve_repo_path(argv$scatter_data_table %||% file.path(table_dir, "scatter_analytical_vs_simulation_data.tsv"), root = root, mustWork = FALSE)

  expected_seed_ids <- expected_analytical_seed_ids(run_dir, seed_ids)

  message("Generating analytical trajectories.")
  if (!recompute_analytical && file.exists(analytical_cache_path)) {
    message("Reading cached generated analytical trajectories: ", analytical_cache_path)
    analytical <- read_tsv(analytical_cache_path)
    missing_keys <- analytical_cache_missing_keys(
      analytical = analytical,
      time_points = time_points,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      expected_seed_ids = expected_seed_ids
    )
    analytical <- filter_analytical_trajectories(
      analytical = analytical,
      time_points = time_points,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      seed_ids = seed_ids
    )
    missing_methods <- setdiff(analytical_methods, unique(analytical$analytical_method))
    if (length(missing_methods)) {
      message("Generated analytical cache is missing method(s), rebuilding: ", paste(missing_methods, collapse = ","))
      analytical <- data.frame()
    }
    if (length(missing_keys)) {
      message("Generated analytical cache does not cover the requested grid, rebuilding. Example missing key(s): ", paste(missing_keys, collapse = " | "))
      analytical <- data.frame()
    }
  } else {
    analytical <- data.frame()
  }
  if (!nrow(analytical)) {
    analytical <- generate_analytical_trajectories(
      run_dir = run_dir,
      time_points = time_points,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      n_workers = n_workers,
      seed_ids = seed_ids
    )
    write_tsv(analytical, analytical_cache_path)
    analytical <- filter_analytical_trajectories(
      analytical = analytical,
      time_points = time_points,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods,
      seed_ids = seed_ids
    )
  }
  if (!nrow(analytical)) stop("No analytical trajectory rows were found for the requested O2/time grid.")
  message("Analytical methods available: ", paste(sort(unique(analytical$analytical_method)), collapse = ", "))

  message("Reading seed objective values.")
  objectives <- read_seed_objectives(analysis_dir = analysis_dir, fit_dir = fit_dir)

  if (!cache_all_times && !recompute && file.exists(sim_summary_path)) {
    message("Reading cached simulation summary: ", sim_summary_path)
    sim_summary <- read_tsv(sim_summary_path)
  } else {
    if (cache_all_times) {
      if (!recompute && file.exists(all_time_sim_metric_path)) {
        message("Reading cached all-time simulation metrics: ", all_time_sim_metric_path)
        all_time_sim_metrics <- read_tsv(all_time_sim_metric_path)
      } else {
        if (!dir.exists(simulation_dir)) stop("Simulation directory is required when all-time cache is missing: ", simulation_dir)
        message("Reading simulation task table.")
        tasks <- read_simulation_tasks(
          simulation_dir = simulation_dir,
          simulation_mode = simulation_mode,
          analytical = analytical,
          o2_values = o2_values,
          initial_ploidy_values = initial_ploidy_values,
          simulation_ids = simulation_ids
        )
        if (!nrow(tasks)) stop("No simulation tasks matched the requested seed/O2/initial ploidy/simulation id filters.")
        message("Matched simulation tasks: ", nrow(tasks))
        message("Reading all-time simulation state metrics.")
        all_time_sim_metrics <- read_simulation_metrics(
          tasks = tasks,
          time_points = NULL,
          progress_every = as.integer(argv$progress_every %||% 100L),
          n_workers = n_workers
        )
        if (!nrow(all_time_sim_metrics)) stop("No all-time simulation metrics were read from state trajectories.")
        write_tsv(all_time_sim_metrics, all_time_sim_metric_path)
      }
      sim_metrics <- filter_simulation_metrics(
        sim_metrics = all_time_sim_metrics,
        time_points = time_points,
        o2_values = o2_values,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids,
        seed_ids = unique(analytical$seed_id)
      )
    } else {
      if (!dir.exists(simulation_dir)) stop("Simulation directory is required when --cache_all_times=FALSE: ", simulation_dir)
      message("Reading simulation task table.")
      tasks <- read_simulation_tasks(
        simulation_dir = simulation_dir,
        simulation_mode = simulation_mode,
        analytical = analytical,
        o2_values = o2_values,
        initial_ploidy_values = initial_ploidy_values,
        simulation_ids = simulation_ids
      )
      if (!nrow(tasks)) stop("No simulation tasks matched the requested seed/O2/initial ploidy/simulation id filters.")
      message("Matched simulation tasks: ", nrow(tasks))
      sim_metrics <- read_simulation_metrics(
        tasks = tasks,
        time_points = time_points,
        progress_every = as.integer(argv$progress_every %||% 100L),
        n_workers = n_workers
      )
    }
    if (!nrow(sim_metrics)) stop("No simulation metrics were read from state trajectories.")
    write_tsv(sim_metrics, sim_metric_path)

    sim_summary <- aggregate_replicates(sim_metrics)
    write_tsv(sim_summary, sim_summary_path)
  }

  message("Merging analytical and simulation summaries.")
  scatter_rows <- list()
  fig_dirs <- character()
  for (method in sort(unique(analytical$analytical_method))) {
    method_analytical <- analytical[analytical$analytical_method %in% method, , drop = FALSE]
    scatter_data <- merge_scatter_data(method_analytical, sim_summary, objectives)
    if (!nrow(scatter_data)) {
      warning("No merged analytical-vs-simulation rows were produced for method: ", method)
      next
    }
    method_table_dir <- file.path(table_dir, method_slug(method))
    dir.create(method_table_dir, recursive = TRUE, showWarnings = FALSE)
    method_scatter_data_path <- file.path(method_table_dir, "scatter_analytical_vs_simulation_data.tsv")
    write_tsv(scatter_data, method_scatter_data_path)
    scatter_rows[[method]] <- scatter_data

    message("Drawing scatter plots for method: ", method)
    fig_dirs[[method]] <- make_scatter_outputs(
      dat = scatter_data,
      out_dir = out_dir,
      objective_transform = objective_transform,
      analytical_method = method
    )
  }
  if (!length(fig_dirs)) stop("No scatter figures were produced for any analytical method.")
  combined_scatter_data <- do.call(rbind, scatter_rows)
  if (!is.null(combined_scatter_data) && nrow(combined_scatter_data)) write_tsv(combined_scatter_data, scatter_data_path)

  manifest <- data.frame(
    field = c(
      "simulation_dir", "analysis_dir", "fit_dir", "out_dir", "simulation_mode",
      "time_points", "o2_values", "initial_ploidy_values", "simulation_ids",
      "objective_transform", "cache_all_times", "n_workers", "run_dir", "seed_ids",
      "analytical_methods", "recompute_analytical", "analytical_cache_table",
      "all_time_simulation_metric_table", "simulation_metric_table", "simulation_summary_table",
      "scatter_data_table", "figure_dirs"
    ),
    value = c(
      simulation_dir, analysis_dir, fit_dir, out_dir, simulation_mode,
      paste(time_points, collapse = ","), paste(o2_values, collapse = ","),
      paste(initial_ploidy_values, collapse = ","), paste(simulation_ids, collapse = ","),
      objective_transform, as.character(cache_all_times), as.character(n_workers), run_dir, paste(seed_ids, collapse = ","),
      paste(analytical_methods, collapse = ","), as.character(recompute_analytical), analytical_cache_path,
      all_time_sim_metric_path, sim_metric_path, sim_summary_path, scatter_data_path, paste(fig_dirs, collapse = ",")
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(table_dir, "scatter_run_manifest.tsv"))
  message("Done. Scatter figures written to: ", paste(fig_dirs, collapse = ", "))
  invisible(fig_dirs)
}

if (identical(environment(), globalenv())) {
  main()
}
