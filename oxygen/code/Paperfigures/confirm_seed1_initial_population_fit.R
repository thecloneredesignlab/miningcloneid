#!/usr/bin/env Rscript

options(digits = 17)

DEFAULT_ORIGINAL_SEED_DIR <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed/seed1"

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (length(x) > 1L) return(x)
  if (is.na(x) || !nzchar(trimws(as.character(x)))) y else x
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1L) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

script_dir <- function() {
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
  if (length(frame_files) > 0L) return(dirname(frame_files[[length(frame_files)]]))
  normalizePath(getwd(), mustWork = FALSE)
}

read_tsv <- function(path) {
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

read_csv_plain <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

metric_map <- function(path) {
  tab <- read_tsv(path)
  if (!all(c("metric", "value") %in% names(tab))) {
    stop("Summary file is not a metric/value table: ", path)
  }
  stats::setNames(as.character(tab$value), as.character(tab$metric))
}

metric_value <- function(metrics, key, default = NA_character_) {
  if (!key %in% names(metrics)) return(default)
  val <- metrics[[key]]
  if (length(val) == 0L || is.na(val) || !nzchar(trimws(as.character(val)))) default else val
}

read_best_params <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  if (!file.exists(path)) stop("Missing best_params.tsv: ", path)
  tab <- read_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain parameter/value columns: ", path)
  }
  stats::setNames(suppressWarnings(as.numeric(tab$value)), as.character(tab$parameter))
}

read_fit_args_from_seed <- function(seed_dir) {
  path <- file.path(seed_dir, "run_effective_args.tsv")
  if (!file.exists(path)) stop("Missing run_effective_args.tsv: ", path)
  tab <- read_tsv(path)
  req <- c("source", "key", "value")
  if (!all(req %in% names(tab))) stop("Invalid run_effective_args.tsv: ", path)
  tab <- tab[tab$source == "fit_command", , drop = FALSE]
  vals <- as.list(as.character(tab$value))
  names(vals) <- as.character(tab$key)
  vals
}

read_parameter_table <- function(seed_dir) {
  candidates <- c(
    file.path(seed_dir, "parameter_table.csv"),
    file.path(dirname(seed_dir), "parameter_table.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) stop("Missing parameter_table.csv for ", seed_dir)
  tab <- read_csv_plain(existing[[1L]])
  req <- c("param_name", "estimate", "param_prototype")
  if (!all(req %in% names(tab))) {
    stop("parameter_table.csv missing columns: ", paste(setdiff(req, names(tab)), collapse = ", "))
  }
  tab$estimate <- tolower(trimws(as.character(tab$estimate))) %in% c("true", "t", "1", "yes", "y")
  tab
}

natural_to_transformed_initialpop <- function(initial_csv, seed, original_seed_dir) {
  if (!file.exists(initial_csv)) stop("Missing initial population CSV: ", initial_csv)
  pop <- read_csv_plain(initial_csv)
  if (!"seed" %in% names(pop)) stop("Initial population CSV has no seed column: ", initial_csv)
  pop <- pop[as.integer(pop$seed) == as.integer(seed), , drop = FALSE]
  if (!nrow(pop)) stop("No rows found for seed ", seed, " in ", initial_csv)

  param_tab <- read_parameter_table(original_seed_dir)
  opt_tab <- param_tab[param_tab$estimate, , drop = FALSE]
  transformed <- matrix(NA_real_, nrow = nrow(pop), ncol = nrow(opt_tab))
  colnames(transformed) <- as.character(opt_tab$param_name)
  for (j in seq_len(nrow(opt_tab))) {
    natural_name <- as.character(opt_tab$param_prototype[[j]])
    transformed_name <- as.character(opt_tab$param_name[[j]])
    if (!natural_name %in% names(pop)) {
      stop("Initial population CSV is missing natural parameter column: ", natural_name)
    }
    vals <- suppressWarnings(as.numeric(pop[[natural_name]]))
    if (any(!is.finite(vals))) {
      stop("Initial population has non-finite values for ", natural_name)
    }
    if (startsWith(transformed_name, "log10_")) {
      if (any(vals <= 0)) stop("Cannot log10-transform non-positive values for ", natural_name)
      vals <- log10(vals)
    }
    transformed[, j] <- vals
  }
  transformed
}

make_comparison_table <- function(original_seed_dir, confirm_seed_dir, output_dir) {
  original_params <- read_best_params(original_seed_dir)
  confirm_params <- read_best_params(confirm_seed_dir)
  param_names <- c(names(original_params), setdiff(names(confirm_params), names(original_params)))

  rows <- vector("list", length(param_names))
  for (i in seq_along(param_names)) {
    nm <- param_names[[i]]
    old <- unname(original_params[[nm]])
    new <- unname(confirm_params[[nm]])
    rows[[i]] <- data.frame(
      parameter_or_metric = nm,
      original = old,
      specified_initial_population = new,
      difference = new - old,
      stringsAsFactors = FALSE
    )
  }

  original_summary <- metric_map(file.path(original_seed_dir, "fit_summary.tsv"))
  confirm_summary <- metric_map(file.path(confirm_seed_dir, "fit_summary.tsv"))
  metric_names <- c(
    "objective",
    "optimizer_deoptim_objective",
    "optimizer_local_objective",
    "optimizer_iter_completed",
    "optimizer_iter_target",
    "optimizer_local_accepted",
    "optimizer_local_convergence"
  )
  metric_rows <- lapply(metric_names, function(nm) {
    old_raw <- metric_value(original_summary, nm)
    new_raw <- metric_value(confirm_summary, nm)
    old_num <- suppressWarnings(as.numeric(old_raw))
    new_num <- suppressWarnings(as.numeric(new_raw))
    numeric_pair <- is.finite(old_num) && is.finite(new_num)
    data.frame(
      parameter_or_metric = nm,
      original = if (numeric_pair) old_num else old_raw,
      specified_initial_population = if (numeric_pair) new_num else new_raw,
      difference = if (numeric_pair) new_num - old_num else if (identical(old_raw, new_raw)) "0" else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, c(rows, metric_rows))
  out_path <- file.path(output_dir, "seed1_initial_population_reproducibility_comparison.csv")
  utils::write.csv(out, out_path, quote = FALSE, row.names = FALSE)
  message("Wrote comparison table: ", out_path)
  invisible(out)
}

main <- function() {
  raw_command_args <- base::commandArgs
  wrapper_script_dir <- script_dir()
  assign(
    "commandArgs",
    function(trailingOnly = FALSE) {
      args <- raw_command_args(trailingOnly = trailingOnly)
      if (!isTRUE(trailingOnly)) {
        args <- args[!grepl("^--file=", args)]
      }
      args
    },
    envir = .GlobalEnv
  )
  on.exit({
    if (exists("commandArgs", envir = .GlobalEnv, inherits = FALSE)) {
      rm("commandArgs", envir = .GlobalEnv)
    }
  }, add = TRUE)

  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  seed <- as.integer(argv$seed %||% "1")
  default_paperfigures_dir <- file.path(dirname(dirname(wrapper_script_dir)), "results", "PaperFigures")
  default_initial_population_csv <- file.path(default_paperfigures_dir, "invivo_deoptim_initial_population.csv")
  default_output_dir <- file.path(default_paperfigures_dir, "ConfirmFits")
  original_seed_dir <- normalizePath(argv$original_seed_dir %||% DEFAULT_ORIGINAL_SEED_DIR, mustWork = FALSE)
  initial_csv <- normalizePath(path.expand(argv$initial_population_csv %||% default_initial_population_csv), mustWork = FALSE)
  output_dir <- normalizePath(path.expand(argv$output_dir %||% default_output_dir), mustWork = FALSE)
  confirm_seed_dir <- file.path(output_dir, "seed1_from_initial_population")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(confirm_seed_dir)) unlink(confirm_seed_dir, recursive = TRUE, force = TRUE)
  dir.create(confirm_seed_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(original_seed_dir)) stop("Original seed directory does not exist: ", original_seed_dir)

  this_script_dir <- wrapper_script_dir
  workflow_root <- normalizePath(file.path(dirname(this_script_dir), "O2_supply_demand_MAP"), mustWork = FALSE)
  optimizer_dir <- file.path(workflow_root, "optimizer")
  backend_path <- normalizePath(
    argv$backend %||% file.path(workflow_root, "util", "o2_supply_demand_map_fit_invivo_backend.R"),
    mustWork = FALSE
  )
  if (!file.exists(backend_path)) stop("Backend script not found: ", backend_path)

  initialpop_t <- natural_to_transformed_initialpop(
    initial_csv = initial_csv,
    seed = seed,
    original_seed_dir = original_seed_dir
  )

  fit_args <- read_fit_args_from_seed(original_seed_dir)
  fit_args$seed <- as.character(seed)
  fit_args$out_dir <- confirm_seed_dir
  fit_args$mode <- "fit_seed"
  fit_args$fit_invivo <- "TRUE"

  backend_env <- new.env(parent = globalenv())
  backend_env$commandArgs <- function(trailingOnly = FALSE) {
    args <- raw_command_args(trailingOnly = trailingOnly)
    if (!isTRUE(trailingOnly)) {
      args <- args[!grepl("^--file=", args)]
    }
    args
  }
  sys.source(backend_path, envir = backend_env, chdir = TRUE)
  backend_get_script_dir <- function(default = getwd()) {
    normalizePath(optimizer_dir, mustWork = FALSE)
  }
  assign("get_script_dir", backend_get_script_dir, envir = backend_env)
  assign("o2sd_get_script_dir", backend_get_script_dir, envir = backend_env)
  environment(backend_env$main_fit_single_seed)$get_script_dir <- backend_get_script_dir
  environment(backend_env$main)$get_script_dir <- backend_get_script_dir
  message("Backend path: ", backend_path)
  message("Backend optimizer dir: ", backend_get_script_dir())

  original_build_de_initialpop <- get(".build_de_initialpop", envir = backend_env, inherits = FALSE)
  assign(
    ".build_de_initialpop",
    function(np, lower, upper, init_use = NULL, mode = "hybrid", uniform_frac = 0.3, sigma_frac = 0.1) {
      regenerated_plan <- original_build_de_initialpop(
        np = np,
        lower = lower,
        upper = upper,
        init_use = init_use,
        mode = mode,
        uniform_frac = uniform_frac,
        sigma_frac = sigma_frac
      )
      regenerated <- if (is.list(regenerated_plan) && !is.null(regenerated_plan$pop)) {
        regenerated_plan$pop
      } else {
        regenerated_plan
      }
      requested <- initialpop_t[, names(lower), drop = FALSE]
      if (nrow(requested) != as.integer(np)) {
        stop("Initial population row count mismatch: CSV has ", nrow(requested), ", NP is ", np)
      }
      if (any(!is.finite(requested))) stop("Initial population contains non-finite transformed values")
      if (any(requested < matrix(lower, nrow = nrow(requested), ncol = length(lower), byrow = TRUE) - 1e-8) ||
          any(requested > matrix(upper, nrow = nrow(requested), ncol = length(upper), byrow = TRUE) + 1e-8)) {
        stop("Initial population contains values outside transformed optimizer bounds.")
      }
      if (is.matrix(regenerated) && all(names(lower) %in% colnames(regenerated))) {
        max_abs_delta <- max(abs(requested - regenerated[, names(lower), drop = FALSE]))
        message("Using supplied initial population; max abs delta vs regenerated population = ", signif(max_abs_delta, 17))
      } else {
        message("Using supplied initial population; regenerated population was not matrix-like for delta check.")
      }
      if (is.list(regenerated_plan)) {
        regenerated_plan$pop <- requested
        regenerated_plan$mode_effective <- paste0("specified_", regenerated_plan$mode_effective %||% "initial_population")
        return(regenerated_plan)
      }
      list(
        pop = requested,
        mode_effective = "specified_initial_population",
        warm_start_used = !is.null(init_use),
        n_local = NA_integer_,
        n_uniform = NA_integer_
      )
    },
    envir = backend_env
  )

  message("Running confirm fit into: ", confirm_seed_dir)
  backend_env$main_fit_single_seed(fit_args)

  comp <- make_comparison_table(
    original_seed_dir = original_seed_dir,
    confirm_seed_dir = confirm_seed_dir,
    output_dir = output_dir
  )
  diff_num <- suppressWarnings(as.numeric(comp$difference))
  finite_diff <- diff_num[is.finite(diff_num)]
  if (length(finite_diff)) {
    message("Max absolute numeric difference: ", signif(max(abs(finite_diff)), 17))
  }
}

main()
