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
PROCESS_FINGERPRINT_DIR <- file.path(WORKFLOW_ROOT, "analysis", "process_fingerprints")
source(file.path(PROCESS_FINGERPRINT_DIR, "process_fingerprint_utils.R"), local = environment())
source(file.path(PROCESS_FINGERPRINT_DIR, "ploidy_regime_utils.R"), local = environment())

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
    "seed_id", "O2_pct", "O2_key", "dominant_mean_ploidy", "trajectory_regime",
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

fixo2_dominant_attractor_one <- function(seed_id, run_params, model_env, cfg, O2) {
  fm <- tryCatch(fixo2_fixed_matrix(model_env, cfg, run_params, O2), error = function(e) e)
  if (inherits(fm, "error")) {
    return(data.frame(
      seed_id = seed_id,
      O2_pct = O2,
      status = paste0("matrix_failed:", conditionMessage(fm)),
      dominant_mean_N = NA_real_,
      dominant_mean_ploidy = NA_real_,
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

if (sys.nframe() == 0) {
  main()
}
