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

select_pair_best <- function(pair_dir) {
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

main <- function(argv = parse_args()) {
  joint_root <- argv$joint_root
  out_path <- argv$out
  summary_out <- argv$summary_out
  if (is.null(joint_root) || is.null(out_path)) {
    stop("Usage: generate_figure6_joint_overlay.R --joint_root=/abs/path --out=/abs/path.tsv [--summary_out=/abs/path.tsv]", call. = FALSE)
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
  o2_values <- seq(o2_min, o2_max, length.out = o2_n)

  pair_dirs <- list.dirs(joint_root, recursive = FALSE, full.names = TRUE)
  pair_dirs <- sort(pair_dirs[grepl("^fit_joint_", basename(pair_dirs))])
  if (!length(pair_dirs)) stop("No joint pair directories found under: ", joint_root, call. = FALSE)
  best <- do.call(rbind, lapply(pair_dirs, select_pair_best))
  model_env <- sim_env$o2ipa_source_model(sim_env$SCRIPT_DIR)
  curves <- do.call(rbind, lapply(seq_len(nrow(best)), function(i) {
    message("Computing ", best$pair_id[[i]], " ", best$seed_id[[i]], " across ", length(o2_values), " O2 values.")
    compute_pair_curve(best[i, , drop = FALSE], o2_values, model_env)
  }))
  curves <- curves[order(curves$pair_id, curves$O2_pct), , drop = FALSE]
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(curves, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  utils::write.table(best, summary_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote ", out_path, " [", nrow(curves), " x ", ncol(curves), "]")
  message("Wrote ", summary_out, " [", nrow(best), " x ", ncol(best), "]")
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) y else x

if (sys.nframe() == 0L) main()
