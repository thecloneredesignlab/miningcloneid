#!/usr/bin/env Rscript

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!grepl("^--[^=]+=", arg)) next
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    out[[key]] <- sub("^[^=]+=", "", arg)
  }
  out
}

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})
workflow_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
simulation_path <- file.path(workflow_root, "simulation", "fix_o2_simulation.R")
sim_env <- new.env(parent = globalenv())
source(simulation_path, local = sim_env, chdir = TRUE)

read_metric <- function(path, metric) {
  if (!file.exists(path)) return(NA_character_)
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (!all(c("metric", "value") %in% names(tab))) return(NA_character_)
  hit <- tab$value[as.character(tab$metric) == metric]
  if (length(hit)) as.character(hit[[1L]]) else NA_character_
}

seed_number <- function(path) suppressWarnings(as.integer(sub("^seed", "", basename(path))))

collect_pair_valid_seeds <- function(pair_dir) {
  seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
  rows <- lapply(seed_dirs, function(seed_dir) {
    summary_path <- file.path(seed_dir, "fit_summary.tsv")
    status <- read_metric(summary_path, "fit_status")
    objective <- suppressWarnings(as.numeric(read_metric(summary_path, "objective")))
    data.frame(
      pair_id = basename(pair_dir),
      seed = seed_number(seed_dir),
      seed_id = basename(seed_dir),
      objective = objective,
      fit_status = status,
      seed_dir = normalizePath(seed_dir, mustWork = FALSE),
      fit_summary_path = normalizePath(summary_path, mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  ok <- tab$fit_status == "ok" & is.finite(tab$objective) & file.exists(file.path(tab$seed_dir, "joint_best_params_long.tsv"))
  tab <- tab[ok, , drop = FALSE]
  if (!nrow(tab)) stop("No valid joint seeds found under: ", pair_dir, call. = FALSE)
  tab[order(tab$seed), , drop = FALSE]
}

select_pair_best <- function(pair_dir) {
  tab <- collect_pair_valid_seeds(pair_dir)
  tab[order(tab$objective, tab$seed), , drop = FALSE][1L, , drop = FALSE]
}

compute_pair_curve <- function(best_row, o2_values, model_env) {
  fit_dir <- best_row$seed_dir[[1L]]
  param_info <- sim_env$resolve_parameter_values(
    fit_dir = fit_dir,
    best_params_path = NULL,
    simulation = "joint",
    joint_scope = "shared_invivo"
  )
  cfg_raw <- sim_env$read_fit_config(fit_dir)
  run_params <- sim_env$prepare_run_params(
    param_values = param_info$values,
    simulation = "joint",
    cfg = cfg_raw,
    fixed_o2 = max(o2_values)
  )
  cfg <- sim_env$prepare_sim_cfg(cfg_raw, argv = list(), fixed_o2 = max(o2_values), run_params = run_params)
  run_params$O2_growth <- isTRUE(cfg$O2_growth)
  run_params$ploidy_O2_death <- cfg$ploidy_O2_death
  rows <- lapply(o2_values, function(o2) {
    sim_env$fixo2_dominant_attractor_one(
      seed_id = best_row$seed_id[[1L]],
      run_params = run_params,
      model_env = model_env,
      cfg = cfg,
      O2 = o2
    )
  })
  out <- do.call(rbind, rows)
  out$pair_id <- best_row$pair_id[[1L]]
  out$joint_seed <- best_row$seed[[1L]]
  out$objective <- best_row$objective[[1L]]
  out$source_seed_dir <- fit_dir
  out$source_parameter_file <- param_info$path
  out
}

force_effective_p_misseg <- function(run_params, effective_p_misseg) {
  effective_p_misseg <- as.numeric(effective_p_misseg)
  if (length(effective_p_misseg) != 1L || !is.finite(effective_p_misseg) ||
      effective_p_misseg <= 0 || effective_p_misseg >= 1) {
    stop("effective_p_misseg must be finite and within (0, 1).", call. = FALSE)
  }
  rp <- run_params
  rp$p_mis_base <- effective_p_misseg
  rp$p_misseg <- 0
  rp
}

compute_pair_cin_grid <- function(best_row, o2_values, effective_p_misseg_values, model_env) {
  fit_dir <- best_row$seed_dir[[1L]]
  param_info <- sim_env$resolve_parameter_values(
    fit_dir = fit_dir,
    best_params_path = NULL,
    simulation = "joint",
    joint_scope = "shared_invivo"
  )
  cfg_raw <- sim_env$read_fit_config(fit_dir)
  run_params <- sim_env$prepare_run_params(
    param_values = param_info$values,
    simulation = "joint",
    cfg = cfg_raw,
    fixed_o2 = max(o2_values)
  )
  cfg <- sim_env$prepare_sim_cfg(cfg_raw, argv = list(), fixed_o2 = max(o2_values), run_params = run_params)
  run_params$O2_growth <- isTRUE(cfg$O2_growth)
  run_params$ploidy_O2_death <- cfg$ploidy_O2_death
  rows <- lapply(effective_p_misseg_values, function(effective_p_misseg) {
    cin_run_params <- force_effective_p_misseg(run_params, effective_p_misseg)
    do.call(rbind, lapply(o2_values, function(o2) {
      row <- sim_env$fixo2_dominant_attractor_one(
        seed_id = best_row$seed_id[[1L]],
        run_params = cin_run_params,
        model_env = model_env,
        cfg = cfg,
        O2 = o2
      )
      row$actual_effective_p_misseg <- row$population_average_p_misseg
      row$effective_p_misseg <- effective_p_misseg
      row$log10_effective_p_misseg <- log10(effective_p_misseg)
      row$population_average_p_misseg <- effective_p_misseg
      row$log10_population_average_p_misseg <- log10(effective_p_misseg)
      row
    }))
  })
  out <- do.call(rbind, rows)
  out$pair_id <- best_row$pair_id[[1L]]
  out$joint_seed <- best_row$seed[[1L]]
  out$objective <- best_row$objective[[1L]]
  out$source_seed_dir <- fit_dir
  out$source_parameter_file <- param_info$path
  out
}

effective_p_misseg_grid_values <- function(p_min, p_max, n) {
  p_min <- as.numeric(p_min)
  p_max <- as.numeric(p_max)
  n <- as.integer(n)
  if (!is.finite(p_min) || !is.finite(p_max) || p_min <= 0 || p_max >= 1 || p_max <= p_min) {
    stop("--effective_p_misseg_max must be greater than --effective_p_misseg_min, with both values within (0, 1).", call. = FALSE)
  }
  if (!is.finite(n) || n < 2L) stop("--cin_n must be an integer >= 2.", call. = FALSE)
  10^seq(log10(p_min), log10(p_max), length.out = n)
}

main <- function(argv = parse_args()) {
  joint_root <- argv$joint_root
  out_path <- argv$out
  summary_out <- argv$summary_out
  if (is.null(joint_root) || is.null(out_path)) {
    stop(
      "Usage: generate_figure6_joint_overlay.R --joint_root=/abs/path --out=/abs/path.tsv ",
      "[--summary_out=/abs/path.tsv] [--selection_mode=best|all] [--surface_mode=curve|cin_grid] ",
      "[--pair_ids=id1,id2] [--workers=1] [--effective_p_misseg_min=0.005] [--effective_p_misseg_max=0.5] [--cin_n=60]",
      call. = FALSE
    )
  }
  joint_root <- normalizePath(path.expand(joint_root), mustWork = TRUE)
  out_path <- normalizePath(path.expand(out_path), mustWork = FALSE)
  if (is.null(summary_out) || !nzchar(summary_out)) {
    summary_out <- file.path(dirname(out_path), "joint_best_seed_summary.tsv")
  }
  summary_out <- normalizePath(path.expand(summary_out), mustWork = FALSE)
  o2_min <- suppressWarnings(as.numeric(argv$o2_min %||% 0))
  o2_max <- suppressWarnings(as.numeric(argv$o2_max %||% 5))
  o2_n <- suppressWarnings(as.integer(argv$o2_n %||% 201L))
  selection_mode <- match.arg(argv$selection_mode %||% "best", c("best", "all"))
  surface_mode <- match.arg(argv$surface_mode %||% "curve", c("curve", "cin_grid"))
  if (identical(surface_mode, "cin_grid") && !identical(selection_mode, "best")) {
    stop("--surface_mode=cin_grid currently supports --selection_mode=best only.", call. = FALSE)
  }
  workers <- suppressWarnings(as.integer(argv$workers %||% 1L))
  if (!is.finite(workers) || workers < 1L) stop("--workers must be a positive integer.", call. = FALSE)
  o2_values <- seq(o2_min, o2_max, length.out = o2_n)
  effective_p_misseg_min <- if (!is.null(argv$effective_p_misseg_min) && nzchar(argv$effective_p_misseg_min)) {
    argv$effective_p_misseg_min
  } else {
    10^as.numeric(argv$cin_log10_min %||% log10(0.005))
  }
  effective_p_misseg_max <- if (!is.null(argv$effective_p_misseg_max) && nzchar(argv$effective_p_misseg_max)) {
    argv$effective_p_misseg_max
  } else {
    10^as.numeric(argv$cin_log10_max %||% log10(0.5))
  }
  cin_values <- effective_p_misseg_grid_values(
    effective_p_misseg_min,
    effective_p_misseg_max,
    argv$cin_n %||% 60L
  )

  pair_dirs <- list.dirs(joint_root, recursive = FALSE, full.names = TRUE)
  pair_dirs <- sort(pair_dirs[grepl("^fit_joint_", basename(pair_dirs))])
  if (!is.null(argv$pair_ids) && nzchar(argv$pair_ids)) {
    requested_pair_ids <- strsplit(argv$pair_ids, ",", fixed = TRUE)[[1L]]
    missing_pair_ids <- setdiff(requested_pair_ids, basename(pair_dirs))
    if (length(missing_pair_ids)) {
      stop("Requested pair directory/directories not found: ", paste(missing_pair_ids, collapse = ", "), call. = FALSE)
    }
    pair_dirs <- pair_dirs[basename(pair_dirs) %in% requested_pair_ids]
  }
  if (!length(pair_dirs)) stop("No joint pair directories found under: ", joint_root, call. = FALSE)
  selected <- if (identical(selection_mode, "best")) {
    do.call(rbind, lapply(pair_dirs, select_pair_best))
  } else {
    do.call(rbind, lapply(pair_dirs, collect_pair_valid_seeds))
  }
  selected <- selected[order(selected$pair_id, selected$seed), , drop = FALSE]
  rownames(selected) <- NULL
  model_env <- sim_env$o2ipa_source_model(sim_env$SCRIPT_DIR)
  compute_one <- function(i) {
    if (identical(selection_mode, "best") || i == 1L || i %% 50L == 0L || i == nrow(selected)) {
      message(
        "Computing ", i, "/", nrow(selected), ": ",
        selected$pair_id[[i]], " ", selected$seed_id[[i]],
        " across ", length(o2_values), " O2 values",
        if (identical(surface_mode, "cin_grid")) paste0(" x ", length(cin_values), " effective p_misseg values.") else "."
      )
    }
    if (identical(surface_mode, "cin_grid")) {
      compute_pair_cin_grid(selected[i, , drop = FALSE], o2_values, cin_values, model_env)
    } else {
      compute_pair_curve(selected[i, , drop = FALSE], o2_values, model_env)
    }
  }
  curve_rows <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(selected)), compute_one, mc.cores = workers)
  } else {
    lapply(seq_len(nrow(selected)), compute_one)
  }
  curves <- do.call(rbind, curve_rows)
  curves <- if (identical(surface_mode, "cin_grid")) {
    curves[order(curves$pair_id, curves$joint_seed, curves$O2_pct, curves$effective_p_misseg), , drop = FALSE]
  } else {
    curves[order(curves$pair_id, curves$joint_seed, curves$O2_pct), , drop = FALSE]
  }
  if (identical(surface_mode, "cin_grid")) {
    compact_columns <- c(
      "pair_id", "joint_seed", "seed_id", "objective", "O2_pct", "status",
      "effective_p_misseg", "log10_effective_p_misseg",
      "population_average_p_misseg", "log10_population_average_p_misseg",
      "actual_effective_p_misseg", "dominant_mean_ploidy", "spectral_gap", "dominant_growth_rate"
    )
    if (!all(compact_columns %in% names(curves))) {
      stop("CIN-grid fixed-O2 output is missing: ", paste(setdiff(compact_columns, names(curves)), collapse = ", "), call. = FALSE)
    }
    curves <- curves[, compact_columns, drop = FALSE]
  } else if (identical(selection_mode, "all")) {
    compact_columns <- c(
      "pair_id", "joint_seed", "seed_id", "objective", "O2_pct", "status",
      "population_average_p_misseg", "dominant_mean_ploidy", "spectral_gap"
    )
    if (!all(compact_columns %in% names(curves))) {
      stop("All-seed fixed-O2 output is missing: ", paste(setdiff(compact_columns, names(curves)), collapse = ", "), call. = FALSE)
    }
    curves <- curves[, compact_columns, drop = FALSE]
  }
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(curves, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  selected$selection_mode <- selection_mode
  selected$surface_mode <- surface_mode
  if (identical(surface_mode, "cin_grid")) {
    selected$effective_p_misseg_min <- min(cin_values)
    selected$effective_p_misseg_max <- max(cin_values)
    selected$log10_effective_p_misseg_min <- min(log10(cin_values))
    selected$log10_effective_p_misseg_max <- max(log10(cin_values))
    selected$cin_n <- length(cin_values)
  }
  utils::write.table(selected, summary_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote ", out_path, " [", nrow(curves), " x ", ncol(curves), "]")
  message("Wrote ", summary_out, " [", nrow(selected), " x ", ncol(selected), "]")
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) y else x

if (sys.nframe() == 0L) main()
