#!/usr/bin/env Rscript

# Fixed-O2 wrapper for the O2_supply_demand_MAP karyotype model.

suppressPackageStartupMessages(library(Matrix))

.fix_o2_script_dir <- local({
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
  getwd()
})

SCRIPT_DIR <- normalizePath(.fix_o2_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, "../../../.."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())

MODEL_PATH <- file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R")
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(MODEL_PATH))
source(MODEL_PATH, local = environment())
rm(.fix_o2_script_dir)

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
  task_list_path <- file.path(task_root, "task_list.tsv")
  write_tsv(tasks, task_list_path)

  message("Fixed-O2 batch simulation")
  message("  run_dir: ", run_dir)
  message("  seeds: ", paste(seed_tab$seed_id, collapse = ","))
  message("  o2_values: ", paste(o2_values, collapse = ","))
  message("  initial_ploidy_values: ", paste(initial_ploidy_values, collapse = ","))
  message("  n_sim: ", n_sim)
  message("  task_count: ", nrow(tasks))
  message("  n_core: ", n_core)
  message("  task_list: ", task_list_path)
  if (as_bool(argv$build_task_list_only, FALSE)) {
    message("  build_task_list_only: TRUE")
    return(invisible(tasks))
  }

  task_rows <- split(tasks, seq_len(nrow(tasks)))
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

  output_paths <- list(
    run_config = file.path(out_dir, "run_config.tsv"),
    parameters_used = file.path(out_dir, "parameters_used.tsv"),
    population = file.path(out_dir, "population_trajectory.tsv"),
    rate = file.path(out_dir, "rate_trajectory.tsv.gz"),
    state = file.path(out_dir, "state_trajectory.tsv.gz")
  )
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

if (sys.nframe() == 0) {
  main()
}
