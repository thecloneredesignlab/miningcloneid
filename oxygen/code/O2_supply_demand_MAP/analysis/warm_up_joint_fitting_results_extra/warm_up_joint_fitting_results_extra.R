#!/usr/bin/env Rscript

# Extra analyses for joint multi-warmup fitting results.
# This script regenerates the joint DEoptim initial population from the same
# joint backend code and random seeds used during fitting.

suppressPackageStartupMessages({
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
})

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  normalizePath(getwd(), mustWork = FALSE)
}

find_repo_root <- function(start = local_script_dir()) {
  cur <- normalizePath(start, mustWork = FALSE)
  if (file.exists(cur) && !dir.exists(cur)) cur <- dirname(cur)
  for (i in seq_len(10L)) {
    if (file.exists(file.path(cur, "oxygen", "code", "O2_supply_demand_MAP", "util", "o2_supply_demand_map_shared.R"))) {
      return(normalizePath(cur, mustWork = FALSE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not locate repository root from: ", start, call. = FALSE)
}

SCRIPT_DIR <- local_script_dir()
REPO_ROOT <- find_repo_root(SCRIPT_DIR)
WORKFLOW_DIR <- file.path(REPO_ROOT, "oxygen", "code", "O2_supply_demand_MAP")
PL_UTILS <- file.path(WORKFLOW_DIR, "analysis", "best_fit_parameter_feature", "util", "parameter_landscape_shared_utils.R")
JOINT_BACKEND <- file.path(WORKFLOW_DIR, "util", "o2_supply_demand_map_fit_joint_backend.R")
DENSE_DIR <- file.path(WORKFLOW_DIR, "analysis", "best_fit_parameter_feature", "03_dense-grid_monotonicity_classification")

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[[1]], fixed = TRUE)
    val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  if (is.logical(x[[1]])) return(isTRUE(x[[1]]))
  val <- tolower(trimws(as.character(x[[1]])))
  if (val %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  isTRUE(default)
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  val <- suppressWarnings(as.integer(x[[1]]))
  if (length(val) && !is.na(val)) val else default
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  val <- suppressWarnings(as.numeric(x[[1]]))
  if (length(val) && is.finite(val)) val else default
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) return(default)
  vals <- trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])
  vals <- suppressWarnings(as.numeric(vals))
  vals[is.finite(vals)]
}

is_finite_int <- function(x) {
  is.finite(x) && !is.na(x)
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

write_table <- function(x, path, sep = "\t") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = sep, quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

write_tsv <- function(x, path) write_table(x, path, sep = "\t")
write_csv <- function(x, path) write_table(x, path, sep = ",")

rbind_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(d) {
    missing <- setdiff(cols, names(d))
    for (col in missing) d[[col]] <- NA
    d[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

read_metric_map <- function(path) {
  if (!file.exists(path)) return(list())
  tab <- read_tsv(path)
  if (!all(c("metric", "value") %in% names(tab))) return(list())
  vals <- as.list(tab$value)
  names(vals) <- tab$metric
  vals
}

metric_value <- function(metrics, key, default = NA_character_) {
  val <- metrics[[key]]
  if (is.null(val) || !length(val)) default else as.character(val[[1]])
}

metric_num <- function(metrics, key, default = NA_real_) {
  suppressWarnings(as.numeric(metric_value(metrics, key, default)))
}

seed_number_from_dir <- function(path) {
  suppressWarnings(as.integer(sub("^seed", "", basename(path))))
}

order_seed_dirs <- function(seed_dirs) {
  seed_dirs[order(seed_number_from_dir(seed_dirs), basename(seed_dirs))]
}

safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  x
}

default_input_root <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_200seed_20260708_173447")
}

default_output_root <- function(input_root) {
  file.path(
    REPO_ROOT, "oxygen", "results", "analysis", "warm_up_joint_fitting_results_extra",
    basename(normalizePath(input_root, mustWork = FALSE))
  )
}

stage_values <- function(x) {
  vals <- trimws(strsplit(as.character(x %||% "all"), ",", fixed = TRUE)[[1]])
  vals <- vals[nzchar(vals)]
  vals <- gsub("-", "_", tolower(vals), fixed = TRUE)
  if ("all" %in% vals) return(c("prepare", "embedding", "curve", "summary"))
  match.arg(vals, c("prepare", "embedding", "curve", "summary"), several.ok = TRUE)
}

source_env <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  file_command_args <- function(trailingOnly = FALSE) {
    if (isTRUE(trailingOnly)) character(0) else paste0("--file=", path)
  }
  env$commandArgs <- file_command_args
  had_global_command_args <- exists("commandArgs", envir = globalenv(), inherits = FALSE)
  old_global_command_args <- if (had_global_command_args) get("commandArgs", envir = globalenv(), inherits = FALSE) else NULL
  assign("commandArgs", file_command_args, envir = globalenv())
  on.exit({
    if (had_global_command_args) {
      assign("commandArgs", old_global_command_args, envir = globalenv())
    } else {
      rm("commandArgs", envir = globalenv())
    }
  }, add = TRUE)
  sys.source(path, envir = env, chdir = TRUE)
  env
}

load_parameter_landscape_utils <- function() {
  if (!file.exists(PL_UTILS)) stop("Missing parameter landscape utility script: ", PL_UTILS, call. = FALSE)
  source(PL_UTILS, local = globalenv())
  invisible(TRUE)
}

load_joint_backend_env <- function() {
  if (!file.exists(JOINT_BACKEND)) stop("Missing joint backend: ", JOINT_BACKEND, call. = FALSE)
  source_env(JOINT_BACKEND)
}

remap_path_value <- function(x, from, to) {
  x <- as.character(x)
  if (!nzchar(from) || !nzchar(to)) return(x)
  ifelse(startsWith(x, from), paste0(to, substring(x, nchar(from) + 1L)), x)
}

remap_argv_paths <- function(argv, path_map_from, path_map_to) {
  if (!length(argv)) return(argv)
  for (key in names(argv)) {
    val <- argv[[key]]
    if (is.null(val) || !length(val)) next
    argv[[key]] <- remap_path_value(val[[1]], path_map_from, path_map_to)
  }
  argv
}

read_effective_argv <- function(seed_dir, path_map_from, path_map_to) {
  path <- file.path(seed_dir, "run_effective_args.tsv")
  if (!file.exists(path)) stop("Missing run_effective_args.tsv: ", path, call. = FALSE)
  tab <- read_tsv(path)
  if (!all(c("source", "key", "value") %in% names(tab))) {
    stop("run_effective_args.tsv must contain source/key/value: ", path, call. = FALSE)
  }
  tab <- tab[as.character(tab$source) == "fit_command", , drop = FALSE]
  argv <- as.list(tab$value)
  names(argv) <- tab$key
  remap_argv_paths(argv, path_map_from, path_map_to)
}

discover_pair_dirs <- function(input_root, max_pairs = NA_integer_) {
  if (!dir.exists(input_root)) stop("Input root does not exist: ", input_root, call. = FALSE)
  dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^fit_joint_", basename(dirs))]
  dirs <- dirs[vapply(dirs, function(d) {
    any(file.exists(file.path(list.dirs(d, recursive = FALSE, full.names = TRUE), "best_params.tsv")))
  }, logical(1))]
  dirs <- dirs[order(basename(dirs))]
  if (is_finite_int(max_pairs) && max_pairs > 0L) dirs <- dirs[seq_len(min(length(dirs), max_pairs))]
  if (!length(dirs)) stop("No fit_joint_* pair directories found under: ", input_root, call. = FALSE)
  dirs
}

parse_pair_label <- function(pair_dir) {
  run_prefix <- basename(pair_dir)
  warmup_label <- sub("^fit_joint_", "", run_prefix)
  out <- list(
    joint_run_prefix = run_prefix,
    warmup_label = warmup_label,
    method = NA_character_,
    invivo_warmup_seed = NA_integer_,
    invivo_cluster = NA_character_,
    invivo_subcluster = NA_character_,
    invitro_warmup_seed = NA_integer_
  )
  m <- regexec("^([A-Za-z0-9]+)_vi_seed([0-9]+)_([A-Za-z0-9]+)(Sc[0-9]+)_vt_seed([0-9]+)$", warmup_label)
  hit <- regmatches(warmup_label, m)[[1]]
  if (length(hit) == 6L) {
    out$method <- hit[[2]]
    out$invivo_warmup_seed <- as.integer(hit[[3]])
    out$invivo_cluster <- hit[[4]]
    out$invivo_subcluster <- hit[[5]]
    out$invitro_warmup_seed <- as.integer(hit[[6]])
  }
  out
}

read_best_param_vector <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  if (!file.exists(path)) stop("Missing best_params.tsv: ", path, call. = FALSE)
  tab <- read_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain parameter/value columns: ", path, call. = FALSE)
  }
  vals <- suppressWarnings(as.numeric(tab$value))
  names(vals) <- as.character(tab$parameter)
  vals[!duplicated(names(vals))]
}

vector_row <- function(vals, cols) {
  out <- rep(NA_real_, length(cols))
  names(out) <- cols
  hit <- intersect(cols, names(vals))
  out[hit] <- as.numeric(vals[hit])
  as.data.frame(as.list(out), check.names = FALSE)
}

decode_invivo_from_joint_vector <- function(par_t, ctx, joint_env, output_params) {
  context_vectors <- joint_env$joint_build_context_specific_transformed_vectors(par_t, ctx)
  params <- joint_env$INVIVO_ENV$decode_params(
    context_vectors$invivo_par,
    fit_treatment = ctx$invivo$cfg$fit_treatment,
    fit_tau_O2 = ctx$invivo$cfg$fit_tau_O2,
    cfg = ctx$invivo$cfg
  )
  params <- params[vapply(params, is.numeric, logical(1))]
  params <- joint_env$filter_family_specific_run_params_for_output_common(params)
  vals <- as.numeric(unlist(params, use.names = FALSE))
  names(vals) <- names(params)
  vals[output_params]
}

joint_initial_population_to_invivo <- function(pop, ctx, joint_env, output_params) {
  rows <- vector("list", nrow(pop))
  for (i in seq_len(nrow(pop))) {
    vals <- decode_invivo_from_joint_vector(pop[i, ], ctx, joint_env, output_params)
    rows[[i]] <- as.data.frame(as.list(vals), check.names = FALSE)
  }
  do.call(rbind, rows)
}

collect_output_param_names <- function(best_rows) {
  cols <- unique(unlist(lapply(best_rows, names), use.names = FALSE))
  drop <- c("point_type", "synthetic_seed_id", "synthetic_seed_number", "pair_id", "joint_seed")
  setdiff(cols, drop)
}

build_prepare_outputs <- function(args) {
  load_parameter_landscape_utils()
  joint_env <- load_joint_backend_env()

  input_root <- normalizePath(path.expand(args$input_root %||% default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% default_output_root(input_root)), mustWork = FALSE)
  path_map_from <- as.character(args$path_map_from %||% "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling")
  path_map_to <- normalizePath(path.expand(args$path_map_to %||% REPO_ROOT), mustWork = FALSE)
  max_pairs <- as_int(args$max_pairs, NA_integer_)
  max_joint_seeds <- as_int(args$max_joint_seeds, NA_integer_)
  max_initial_members <- as_int(args$max_initial_members, NA_integer_)
  prepare_workers <- as_int(args$prepare_workers %||% args$n_workers, 1L)
  prepare_workers <- max(1L, prepare_workers)
  overwrite <- as_bool(args$overwrite, TRUE)
  if (prepare_workers > 1L && !identical(.Platform$OS.type, "unix")) {
    stop("--prepare_workers > 1 requires a Unix-like platform with fork support.", call. = FALSE)
  }

  tables_dir <- file.path(out_root, "tables")
  synthetic_run_dir <- file.path(out_root, "curve_classification", "joint_invivo_seed_run")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(synthetic_run_dir, recursive = TRUE, showWarnings = FALSE)

  initial_path <- file.path(tables_dir, "joint_initial_population_invivo_params.tsv")
  best_path <- file.path(tables_dir, "joint_best_invivo_params.tsv")
  manifest_path <- file.path(tables_dir, "joint_seed_manifest.tsv")
  synthetic_manifest_path <- file.path(tables_dir, "joint_best_curve_synthetic_run_manifest.tsv")
  if (!overwrite && all(file.exists(c(initial_path, best_path, manifest_path, synthetic_manifest_path)))) {
    message("Reusing prepare outputs under: ", out_root)
    return(invisible(list(
      initial = initial_path,
      best = best_path,
      manifest = manifest_path,
      synthetic_manifest = synthetic_manifest_path,
      synthetic_run_dir = synthetic_run_dir
    )))
  }

  pair_dirs <- discover_pair_dirs(input_root, max_pairs = max_pairs)
  best_vectors <- list()
  seed_records <- list()
  initial_rows <- list()
  best_rows <- list()
  synthetic_rows <- list()
  synthetic_counter <- 0L

  message(
    "Preparing joint initial and best tables from ", length(pair_dirs),
    " pair directories with prepare_workers=", prepare_workers, "."
  )
  for (pair_index in seq_along(pair_dirs)) {
    pair_dir <- pair_dirs[[pair_index]]
    pair_info <- parse_pair_label(pair_dir)
    pair_id <- pair_info$warmup_label
    seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    seed_dirs <- seed_dirs[file.exists(file.path(seed_dirs, "best_params.tsv"))]
    seed_dirs <- order_seed_dirs(seed_dirs)
    if (is_finite_int(max_joint_seeds) && max_joint_seeds > 0L) {
      seed_dirs <- seed_dirs[seq_len(min(length(seed_dirs), max_joint_seeds))]
    }
    if (!length(seed_dirs)) next

    first_argv <- read_effective_argv(seed_dirs[[1L]], path_map_from, path_map_to)
    ctx <- joint_env$build_joint_context(first_argv)
    NP_use <- max(ctx$NP, ctx$joint_np_min_factor * length(ctx$init))
    member_count <- NP_use
    if (is_finite_int(max_initial_members) && max_initial_members > 0L) {
      member_count <- min(member_count, max_initial_members)
    }

    pair_best_vectors <- lapply(seed_dirs, read_best_param_vector)
    best_vectors <- c(best_vectors, pair_best_vectors)
    output_params <- collect_output_param_names(pair_best_vectors)
    output_params <- unique(c(umap_parameter_set("invivo"), output_params))

    pair_workers <- min(prepare_workers, length(seed_dirs))
    message(
      "Pair ", pair_index, "/", length(pair_dirs), ": ", pair_id,
      " (", length(seed_dirs), " seeds; NP_used=", NP_use,
      "; exported initial members=", member_count,
      "; workers=", pair_workers, ")"
    )

    pair_seed_offset <- synthetic_counter
    seed_worker <- function(seed_i) {
      seed_dir <- seed_dirs[[seed_i]]
      joint_seed <- seed_number_from_dir(seed_dir)
      if (!is_finite_int(joint_seed)) stop("Could not parse joint seed from: ", seed_dir, call. = FALSE)
      metrics <- read_metric_map(file.path(seed_dir, "fit_summary.tsv"))
      synthetic_seed_number <- pair_seed_offset + seed_i
      synthetic_seed_id <- paste0("seed", synthetic_seed_number)
      seed_ctx <- ctx
      seed_ctx$seed <- as.integer(joint_seed)
      seed_ctx$raw$seed <- as.character(joint_seed)
      set.seed(seed_ctx$seed)
      pop <- joint_env$joint_deoptim_initial_population(seed_ctx, NP_use)
      if (member_count < nrow(pop)) pop <- pop[seq_len(member_count), , drop = FALSE]

      invivo_init <- joint_initial_population_to_invivo(pop, seed_ctx, joint_env, output_params)
      init_meta <- data.frame(
        point_type = "initial",
        synthetic_seed_id = synthetic_seed_id,
        synthetic_seed_number = synthetic_seed_number,
        pair_id = pair_id,
        joint_run_prefix = pair_info$joint_run_prefix,
        method = pair_info$method,
        invivo_warmup_seed = pair_info$invivo_warmup_seed,
        invivo_cluster = pair_info$invivo_cluster,
        invivo_subcluster = pair_info$invivo_subcluster,
        invitro_warmup_seed = pair_info$invitro_warmup_seed,
        joint_seed = joint_seed,
        initial_member = seq_len(nrow(invivo_init)),
        NP_used = NP_use,
        source_seed_dir = normalizePath(seed_dir, mustWork = FALSE),
        stringsAsFactors = FALSE
      )
      initial_row <- data.frame(init_meta, invivo_init, check.names = FALSE)

      best_vec <- pair_best_vectors[[seed_i]]
      best_param_row <- vector_row(best_vec, output_params)
      best_meta <- data.frame(
        point_type = "best",
        synthetic_seed_id = synthetic_seed_id,
        synthetic_seed_number = synthetic_seed_number,
        pair_id = pair_id,
        joint_run_prefix = pair_info$joint_run_prefix,
        method = pair_info$method,
        invivo_warmup_seed = pair_info$invivo_warmup_seed,
        invivo_cluster = pair_info$invivo_cluster,
        invivo_subcluster = pair_info$invivo_subcluster,
        invitro_warmup_seed = pair_info$invitro_warmup_seed,
        joint_seed = joint_seed,
        objective = metric_num(metrics, "objective"),
        objective_invivo = metric_num(metrics, "objective_invivo"),
        objective_invitro = metric_num(metrics, "objective_invitro"),
        objective_soft_coupling = metric_num(metrics, "objective_soft_coupling"),
        objective_constraints = metric_num(metrics, "objective_constraints"),
        optimizer_start_objective = metric_num(metrics, "optimizer_start_objective"),
        optimizer_start_used = metric_value(metrics, "optimizer_start_used"),
        optimizer_deoptim_objective = metric_num(metrics, "optimizer_deoptim_objective"),
        optimizer_local_objective = metric_num(metrics, "optimizer_local_objective"),
        optimizer_local_accepted = metric_value(metrics, "optimizer_local_accepted"),
        optimizer_iter_completed = metric_num(metrics, "optimizer_iter_completed"),
        deoptim_stop_reason = metric_value(metrics, "deoptim_stop_reason"),
        NP_used = metric_num(metrics, "NP_used"),
        source_seed_dir = normalizePath(seed_dir, mustWork = FALSE),
        stringsAsFactors = FALSE
      )
      best_row <- data.frame(best_meta, best_param_row, check.names = FALSE)

      synthetic_seed_dir <- file.path(synthetic_run_dir, synthetic_seed_id)
      dir.create(synthetic_seed_dir, recursive = TRUE, showWarnings = FALSE)
      best_param_table <- data.frame(parameter = names(best_vec), value = as.numeric(best_vec), stringsAsFactors = FALSE)
      write_tsv(best_param_table, file.path(synthetic_seed_dir, "best_params.tsv"))
      summary_min <- data.frame(
        metric = c(
          "fit_status", "fit_mode", "objective", "objective_total", "objective_data",
          "objective_burden", "objective_ploidy", "deoptim_stop_reason"
        ),
        value = c(
          metric_value(metrics, "fit_status", "ok"),
          "fit_joint_invivo_projection",
          as.character(metric_num(metrics, "objective")),
          as.character(metric_num(metrics, "objective")),
          as.character(metric_num(metrics, "objective_invivo_data")),
          as.character(metric_num(metrics, "objective_invivo_burden")),
          as.character(metric_num(metrics, "objective_invivo_ploidy")),
          metric_value(metrics, "deoptim_stop_reason")
        ),
        stringsAsFactors = FALSE
      )
      write_tsv(summary_min, file.path(synthetic_seed_dir, "fit_summary.tsv"))
      cfg_path <- file.path(seed_dir, "fit_config.rds")
      if (file.exists(cfg_path)) {
        cfg <- readRDS(cfg_path)
        remap_cfg <- function(x) {
          if (is.character(x)) return(remap_path_value(x, path_map_from, path_map_to))
          if (is.list(x)) return(lapply(x, remap_cfg))
          x
        }
        saveRDS(remap_cfg(cfg), file.path(synthetic_seed_dir, "fit_config.rds"))
      }
      synthetic_row <- data.frame(
        synthetic_seed_id = synthetic_seed_id,
        synthetic_seed_number = synthetic_seed_number,
        synthetic_seed_dir = normalizePath(synthetic_seed_dir, mustWork = FALSE),
        pair_id = pair_id,
        joint_seed = joint_seed,
        source_seed_dir = normalizePath(seed_dir, mustWork = FALSE),
        stringsAsFactors = FALSE
      )

      list(
        initial = initial_row,
        best = best_row,
        synthetic = synthetic_row,
        seed_record = best_meta[, setdiff(names(best_meta), output_params), drop = FALSE]
      )
    }

    seed_results <- if (pair_workers > 1L) {
      parallel::mclapply(
        seq_along(seed_dirs),
        seed_worker,
        mc.cores = pair_workers,
        mc.set.seed = FALSE,
        mc.preschedule = FALSE
      )
    } else {
      lapply(seq_along(seed_dirs), seed_worker)
    }
    failed <- vapply(seed_results, inherits, logical(1), "try-error")
    if (any(failed)) {
      bad <- which(failed)[[1L]]
      stop(
        "Prepare worker failed for pair ", pair_id,
        ", seed_dir=", seed_dirs[[bad]], ": ", as.character(seed_results[[bad]]),
        call. = FALSE
      )
    }
    initial_rows <- c(initial_rows, lapply(seed_results, `[[`, "initial"))
    best_rows <- c(best_rows, lapply(seed_results, `[[`, "best"))
    synthetic_rows <- c(synthetic_rows, lapply(seed_results, `[[`, "synthetic"))
    seed_records <- c(seed_records, lapply(seed_results, `[[`, "seed_record"))
    synthetic_counter <- synthetic_counter + length(seed_dirs)
  }

  initial_df <- rbind_fill(initial_rows)
  best_df <- rbind_fill(best_rows)
  manifest_df <- rbind_fill(seed_records)
  synthetic_df <- rbind_fill(synthetic_rows)

  write_tsv(initial_df, initial_path)
  write_tsv(best_df, best_path)
  write_tsv(manifest_df, manifest_path)
  write_tsv(synthetic_df, synthetic_manifest_path)
  message("Wrote initial population table: ", initial_path, " (", nrow(initial_df), " rows)")
  message("Wrote joint best table: ", best_path, " (", nrow(best_df), " rows)")
  message("Wrote synthetic fixed-O2 run dir: ", synthetic_run_dir)

  invisible(list(
    initial = initial_path,
    best = best_path,
    manifest = manifest_path,
    synthetic_manifest = synthetic_manifest_path,
    synthetic_run_dir = synthetic_run_dir
  ))
}

parse_embedding_param_sets <- function(x) {
  vals <- trimws(strsplit(as.character(x %||% "invivo"), ",", fixed = TRUE)[[1]])
  vals <- vals[nzchar(vals)]
  vals <- gsub("-", "_", tolower(vals), fixed = TRUE)
  vals[vals %in% c("all", "both")] <- "invivo,shared"
  vals <- unlist(strsplit(vals, ",", fixed = TRUE), use.names = FALSE)
  vals <- vals[nzchar(vals)]
  vals <- unique(vals)
  match.arg(vals, c("invivo", "shared"), several.ok = TRUE)
}

embedding_param_config <- function(param_set) {
  param_set <- match.arg(param_set, c("invivo", "shared"))
  if (identical(param_set, "shared")) {
    return(list(
      token = "shared",
      label = "shared in vivo/in vitro parameters",
      params = pooled_umap_parameter_set(),
      log10_params = pooled_umap_log10_parameter_set()
    ))
  }
  list(
    token = "invivo_full",
    label = "full in vivo parameters",
    params = umap_parameter_set("invivo"),
    log10_params = umap_log10_parameter_set("invivo")
  )
}

embedding_plot <- function(coords, reduction, x_col, y_col, out_prefix, param_label = "parameters") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping plot: ", out_prefix, call. = FALSE)
    return(invisible(NULL))
  }
  d <- coords
  d$point_type <- factor(d$point_type, levels = c("initial", "best"))
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = d[d$point_type == "initial", , drop = FALSE],
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = pair_id),
      alpha = 0.22, size = 0.36, stroke = 0
    ) +
    ggplot2::geom_point(
      data = d[d$point_type == "best", , drop = FALSE],
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = pair_id),
      shape = 8, alpha = 0.95, size = 2.2, stroke = 0.65
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::labs(
      title = paste0("Joint in vivo ", reduction),
      subtitle = paste0("Regenerated initial population vs best; ", param_label),
      x = x_col, y = y_col, color = "Pair"
    )
  ggplot2::ggsave(paste0(out_prefix, ".pdf"), p, width = 7.2, height = 7.2, useDingbats = FALSE)
  ggplot2::ggsave(paste0(out_prefix, ".png"), p, width = 7.2, height = 7.2, dpi = 300)

  best <- d[d$point_type == "best", , drop = FALSE]
  p_best <- ggplot2::ggplot(best, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = pair_id)) +
    ggplot2::geom_point(shape = 8, alpha = 0.95, size = 2.2, stroke = 0.65) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::labs(
      title = paste0("Joint in vivo ", reduction, ": best only"),
      subtitle = param_label,
      x = x_col, y = y_col, color = "Pair"
    )
  ggplot2::ggsave(paste0(out_prefix, "_best_only.pdf"), p_best, width = 7.2, height = 7.2, useDingbats = FALSE)
  ggplot2::ggsave(paste0(out_prefix, "_best_only.png"), p_best, width = 7.2, height = 7.2, dpi = 300)
  invisible(NULL)
}

build_embedding_outputs <- function(args) {
  load_parameter_landscape_utils()
  input_root <- normalizePath(path.expand(args$input_root %||% default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% default_output_root(input_root)), mustWork = FALSE)
  tables_dir <- file.path(out_root, "tables")
  emb_root_dir <- file.path(out_root, "embeddings")
  fig_root_dir <- file.path(out_root, "figures", "embeddings")
  dir.create(emb_root_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_root_dir, recursive = TRUE, showWarnings = FALSE)

  initial_path <- file.path(tables_dir, "joint_initial_population_invivo_params.tsv")
  best_path <- file.path(tables_dir, "joint_best_invivo_params.tsv")
  if (!file.exists(initial_path) || !file.exists(best_path)) {
    stop("Missing prepare outputs. Run --stage=prepare first.", call. = FALSE)
  }
  initial_df <- read_tsv(initial_path)
  best_df <- read_tsv(best_path)
  meta_cols <- unique(c(
    "point_type", "synthetic_seed_id", "synthetic_seed_number", "pair_id", "joint_run_prefix",
    "method", "invivo_warmup_seed", "invivo_cluster", "invivo_subcluster",
    "invitro_warmup_seed", "joint_seed", "initial_member", "objective",
    "objective_invivo", "objective_invitro", "source_seed_dir"
  ))
  feature_meta <- rbind_fill(list(
    initial_df[, intersect(meta_cols, names(initial_df)), drop = FALSE],
    best_df[, intersect(meta_cols, names(best_df)), drop = FALSE]
  ))
  feature_meta$.row_id <- seq_len(nrow(feature_meta))
  param_sets <- parse_embedding_param_sets(args$embedding_param_set %||% args$embedding_param_sets)

  reductions <- trimws(strsplit(as.character(args$reductions %||% "pca,umap,tsne"), ",", fixed = TRUE)[[1]])
  reductions <- reductions[nzchar(reductions)]
  reductions <- gsub("-", "_", tolower(reductions), fixed = TRUE)
  reductions <- match.arg(reductions, c("pca", "umap", "tsne"), several.ok = TRUE)
  overwrite <- as_bool(args$overwrite, TRUE)
  umap_seed <- as_int(args$umap_seed, 123L)
  umap_neighbors <- as_int(args$umap_neighbors, 80L)
  umap_min_dist <- as_num(args$umap_min_dist, 0.1)
  umap_threads <- as_int(args$umap_threads, 1L)
  tsne_seed <- as_int(args$tsne_seed, 123L)
  tsne_perplexity <- as_num(args$tsne_perplexity, 30)
  tsne_theta <- as_num(args$tsne_theta, 0.5)
  tsne_max_iter <- as_int(args$tsne_max_iter, 1000L)
  run_full_tsne <- as_bool(args$run_full_tsne, FALSE)
  tsne_initial_sample <- as_int(args$tsne_initial_sample, 100000L)

  for (param_set in param_sets) {
    cfg <- embedding_param_config(param_set)
    emb_dir <- file.path(emb_root_dir, cfg$token)
    fig_dir <- file.path(fig_root_dir, cfg$token)
    dir.create(emb_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    initial_feature <- transform_umap_features(initial_df, cfg$params, cfg$log10_params)
    best_feature <- transform_umap_features(best_df, cfg$params, cfg$log10_params)
    all_feature <- rbind(initial_feature, best_feature)
    prep <- prepare_feature_matrix(all_feature, preprocess_mode = "zscore")
    write_tsv(prep$metadata, file.path(emb_dir, paste0("joint_", cfg$token, "_embedding_zscore_metadata.tsv")))

    for (reduction in reductions) {
      table_path <- file.path(emb_dir, paste0("joint_", cfg$token, "_", reduction, "_coordinates.tsv"))
      if (!overwrite && file.exists(table_path)) {
        message("Reusing embedding table: ", table_path)
        next
      }
      mat <- prep$mat
      meta <- feature_meta
      label_suffix <- ""
      if (identical(reduction, "tsne") && !run_full_tsne) {
        initial_idx <- which(meta$point_type == "initial")
        best_idx <- which(meta$point_type == "best")
        if (is_finite_int(tsne_initial_sample) && tsne_initial_sample > 0L && length(initial_idx) > tsne_initial_sample) {
          set.seed(tsne_seed)
          initial_idx <- sort(sample(initial_idx, tsne_initial_sample))
        }
        keep <- sort(c(initial_idx, best_idx))
        mat <- mat[keep, , drop = FALSE]
        meta <- meta[keep, , drop = FALSE]
        label_suffix <- paste0("_sampled_initial", length(initial_idx))
      }
      message("Running ", cfg$token, " ", reduction, " embedding with ", nrow(mat), " rows.")
      if (identical(reduction, "pca")) {
        emb <- run_pca_embedding(
          mat,
          label = paste("joint in vivo initial+best", cfg$label),
          variance_path = file.path(emb_dir, paste0("joint_", cfg$token, "_pca_variance.tsv")),
          variance_figure_prefix = file.path(fig_dir, paste0("joint_", cfg$token, "_pca_variance"))
        )
        coord_cols <- c("PCA1", "PCA2")
      } else if (identical(reduction, "umap")) {
        emb <- run_umap_embedding(
          mat,
          label = paste("joint in vivo initial+best", cfg$label),
          umap_seed = umap_seed,
          n_neighbors = umap_neighbors,
          min_dist = umap_min_dist,
          n_threads = umap_threads
        )
        coord_cols <- c("UMAP1", "UMAP2")
      } else {
        emb <- run_tsne_embedding(
          mat,
          label = paste("joint in vivo initial+best", cfg$label),
          tsne_seed = tsne_seed,
          perplexity = tsne_perplexity,
          theta = tsne_theta,
          max_iter = tsne_max_iter
        )
        coord_cols <- c("tSNE1", "tSNE2")
        table_path <- file.path(emb_dir, paste0("joint_", cfg$token, "_tsne", label_suffix, "_coordinates.tsv"))
      }
      coord_df <- cbind(meta, as.data.frame(emb, check.names = FALSE))
      coord_df$embedding_param_set <- cfg$token
      write_tsv(coord_df, table_path)
      plot_prefix <- file.path(fig_dir, sub("_coordinates.tsv$", "", basename(table_path)))
      embedding_plot(coord_df, toupper(reduction), coord_cols[[1]], coord_cols[[2]], plot_prefix, param_label = cfg$label)
      message("Wrote embedding: ", table_path)
    }
  }
  invisible(TRUE)
}

run_curve_classification <- function(args) {
  input_root <- normalizePath(path.expand(args$input_root %||% default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% default_output_root(input_root)), mustWork = FALSE)
  synthetic_run_dir <- file.path(out_root, "curve_classification", "joint_invivo_seed_run")
  if (!dir.exists(synthetic_run_dir)) stop("Missing synthetic run dir. Run --stage=prepare first: ", synthetic_run_dir, call. = FALSE)

  mono_script <- file.path(DENSE_DIR, "fixed_o2_ploidy_monotonicity.R")
  reg_script <- file.path(DENSE_DIR, "fixed_o2_ploidy_monotonicity_regression_classification.R")
  mono_env <- source_env(mono_script)
  reg_env <- source_env(reg_script)
  dense_out <- file.path(out_root, "curve_classification", "dense-grid_monotonicity_classification")
  reg_out <- file.path(out_root, "curve_classification", "dense-grid_monotonicity_regression_classification")
  o2_grid <- as_num_vec(args$o2_grid, seq(0, 5, by = 0.025))
  reporting_o2 <- as_num_vec(args$reporting_o2, c(0, 0.1, 0.5, 1, 2, 5))
  n_workers <- as_int(args$n_workers, 8L)
  max_seeds <- as_int(args$curve_max_seeds, NA_integer_)
  generate_figures <- as_bool(args$generate_figures, TRUE)
  overwrite <- as_bool(args$overwrite, TRUE)

  mono_args <- list(
    run_dir = synthetic_run_dir,
    out_dir = dense_out,
    o2_grid = paste(format(o2_grid, scientific = FALSE, trim = TRUE), collapse = ","),
    reporting_o2 = paste(format(reporting_o2, scientific = FALSE, trim = TRUE), collapse = ","),
    n_workers = as.character(n_workers),
    max_seeds = as.character(max_seeds),
    overwrite = as.character(overwrite),
    generate_figures = as.character(generate_figures),
    run_validation = as.character(as_bool(args$run_validation, TRUE))
  )
  mono_env$generate_outputs(mono_args)
  reg_args <- list(
    input_dir = dense_out,
    out_dir = reg_out,
    overwrite = as.character(overwrite),
    generate_figures = as.character(generate_figures)
  )
  reg_env$generate_outputs(reg_args)
  invisible(list(pointwise = dense_out, regression = reg_out))
}

plot_pair_curve_summary <- function(summary_df, out_root) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !nrow(summary_df)) return(invisible(NULL))
  fig_dir <- file.path(out_root, "figures", "curve_classification")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = pair_id, y = fraction, fill = curve_class)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = "Pair", y = "Fraction of joint seeds", fill = "Curve class", title = "Joint in vivo curve class composition by pair")
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_fraction_by_pair.pdf"), p, width = 8.4, height = 4.8, useDingbats = FALSE)
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_fraction_by_pair.png"), p, width = 8.4, height = 4.8, dpi = 300)

  p2 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = curve_class, y = pair_id, fill = enrichment)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b83232", midpoint = 1, na.value = "grey90") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = "Curve class", y = "Pair", fill = "Observed / expected", title = "Curve class enrichment by pair")
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_enrichment_heatmap.pdf"), p2, width = 8.4, height = 4.8, useDingbats = FALSE)
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_enrichment_heatmap.png"), p2, width = 8.4, height = 4.8, dpi = 300)
  invisible(NULL)
}

plot_embedding_by_curve_class <- function(master, out_root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  emb_dir <- file.path(out_root, "embeddings")
  fig_dir <- file.path(out_root, "figures", "embedding_curve_class")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  coord_paths <- list.files(emb_dir, pattern = "_coordinates\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (!length(coord_paths)) return(invisible(NULL))
  for (path in coord_paths) {
    coords <- read_tsv(path)
    best <- coords[coords$point_type == "best", , drop = FALSE]
    if (!nrow(best)) next
    best <- merge(
      best,
      master[, c("synthetic_seed_id", "curve_class", "smooth_curve_class", "pair_id"), drop = FALSE],
      by = c("synthetic_seed_id", "pair_id"),
      all.x = TRUE,
      sort = FALSE
    )
    coord_cols <- intersect(c("PCA1", "PCA2", "UMAP1", "UMAP2", "tSNE1", "tSNE2"), names(best))
    if (length(coord_cols) < 2L) next
    x_col <- coord_cols[[1L]]
    y_col <- coord_cols[[2L]]
    emb_base <- paste0(normalizePath(emb_dir, mustWork = FALSE), "/")
    path_norm <- normalizePath(path, mustWork = FALSE)
    rel_path <- if (startsWith(path_norm, emb_base)) substring(path_norm, nchar(emb_base) + 1L) else basename(path_norm)
    stub <- safe_id(sub("_coordinates\\.tsv$", "", gsub("/", "_", rel_path, fixed = TRUE)))
    p <- ggplot2::ggplot(best, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = curve_class)) +
      ggplot2::geom_point(shape = 8, alpha = 0.95, size = 2.2, stroke = 0.65) +
      ggplot2::coord_fixed() +
      ggplot2::facet_wrap(~ pair_id) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::labs(x = x_col, y = y_col, color = "Curve class", title = paste0(stub, ": best points by curve class"))
    ggplot2::ggsave(file.path(fig_dir, paste0(stub, "_best_by_curve_class.pdf")), p, width = 8.2, height = 8.2, useDingbats = FALSE)
    ggplot2::ggsave(file.path(fig_dir, paste0(stub, "_best_by_curve_class.png")), p, width = 8.2, height = 8.2, dpi = 300)
  }
  invisible(NULL)
}

build_summary_outputs <- function(args) {
  input_root <- normalizePath(path.expand(args$input_root %||% default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% default_output_root(input_root)), mustWork = FALSE)
  tables_dir <- file.path(out_root, "tables")
  synthetic_manifest_path <- file.path(tables_dir, "joint_best_curve_synthetic_run_manifest.tsv")
  best_path <- file.path(tables_dir, "joint_best_invivo_params.tsv")
  if (!file.exists(synthetic_manifest_path) || !file.exists(best_path)) {
    stop("Missing prepare outputs. Run --stage=prepare first.", call. = FALSE)
  }
  synthetic_manifest <- read_tsv(synthetic_manifest_path)
  best_df <- read_tsv(best_path)

  reg_by_seed_path <- file.path(
    out_root, "curve_classification", "dense-grid_monotonicity_regression_classification",
    "tables", "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"
  )
  point_by_seed_path <- file.path(
    out_root, "curve_classification", "dense-grid_monotonicity_classification",
    "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv"
  )
  if (file.exists(reg_by_seed_path)) {
    curve_seed <- read_tsv(reg_by_seed_path)
    curve_seed$curve_class <- curve_seed$smooth_curve_class
    curve_source <- "regression_smoothed"
  } else if (file.exists(point_by_seed_path)) {
    curve_seed <- read_tsv(point_by_seed_path)
    curve_source <- "pointwise"
  } else {
    stop("Missing curve classification table. Run --stage=curve first.", call. = FALSE)
  }
  master <- merge(
    best_df,
    curve_seed,
    by.x = "synthetic_seed_id",
    by.y = "seed_id",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("", "_curve")
  )
  master$curve_class_source <- curve_source
  write_tsv(master, file.path(tables_dir, "joint_best_master_table.tsv"))

  tab <- as.data.frame(table(master$pair_id, master$curve_class), stringsAsFactors = FALSE)
  names(tab) <- c("pair_id", "curve_class", "n_seed")
  tab <- tab[tab$n_seed > 0L, , drop = FALSE]
  pair_tot <- aggregate(n_seed ~ pair_id, tab, sum)
  names(pair_tot)[[2L]] <- "pair_total"
  class_tot <- aggregate(n_seed ~ curve_class, tab, sum)
  names(class_tot)[[2L]] <- "class_total"
  total <- sum(tab$n_seed)
  summary_df <- merge(tab, pair_tot, by = "pair_id", all.x = TRUE, sort = FALSE)
  summary_df <- merge(summary_df, class_tot, by = "curve_class", all.x = TRUE, sort = FALSE)
  summary_df$total_seed <- total
  summary_df$fraction <- summary_df$n_seed / summary_df$pair_total
  summary_df$overall_fraction <- summary_df$class_total / total
  summary_df$expected_n <- summary_df$pair_total * summary_df$class_total / total
  summary_df$enrichment <- summary_df$n_seed / summary_df$expected_n
  summary_df$standardized_residual <- (summary_df$n_seed - summary_df$expected_n) / sqrt(pmax(summary_df$expected_n, .Machine$double.eps))
  summary_df <- summary_df[order(summary_df$pair_id, summary_df$curve_class), , drop = FALSE]
  write_tsv(summary_df, file.path(tables_dir, "pair_curve_class_summary.tsv"))

  assoc <- data.frame(
    test = character(0),
    statistic = numeric(0),
    p_value = numeric(0),
    method = character(0),
    stringsAsFactors = FALSE
  )
  xt <- xtabs(n_seed ~ pair_id + curve_class, summary_df)
  if (nrow(xt) > 1L && ncol(xt) > 1L) {
    chi <- suppressWarnings(stats::chisq.test(xt))
    assoc <- data.frame(
      test = "pair_by_curve_class",
      statistic = unname(chi$statistic),
      p_value = chi$p.value,
      method = chi$method,
      stringsAsFactors = FALSE
    )
  }
  write_tsv(assoc, file.path(tables_dir, "pair_curve_class_association_test.tsv"))
  plot_pair_curve_summary(summary_df, out_root)
  plot_embedding_by_curve_class(master, out_root)

  curves_path <- file.path(
    out_root, "curve_classification", "dense-grid_monotonicity_regression_classification",
    "tables", "fixed_o2_ploidy_monotonicity_regression_curves.tsv"
  )
  if (file.exists(curves_path) && requireNamespace("ggplot2", quietly = TRUE)) {
    curves <- read_tsv(curves_path)
    curves <- merge(
      curves,
      master[, c("synthetic_seed_id", "pair_id", "joint_seed", "curve_class"), drop = FALSE],
      by.x = "seed_id",
      by.y = "synthetic_seed_id",
      all.x = TRUE,
      sort = FALSE
    )
    curve_class_col <- "curve_class"
    fig_dir <- file.path(out_root, "figures", "curve_classification")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    p <- ggplot2::ggplot(curves, ggplot2::aes(x = O2_pct, y = smoothed_dominant_mean_ploidy, group = seed_id, color = .data[[curve_class_col]])) +
      ggplot2::geom_line(alpha = 0.22, linewidth = 0.25) +
      ggplot2::facet_wrap(~ pair_id) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(x = "O2 (%)", y = "Smoothed dominant mean ploidy", color = "Curve class", title = "Joint in vivo fixed-O2 curves by pair")
    ggplot2::ggsave(file.path(fig_dir, "joint_invivo_smoothed_curves_by_pair_and_class.pdf"), p, width = 10, height = 6.5, useDingbats = FALSE)
    ggplot2::ggsave(file.path(fig_dir, "joint_invivo_smoothed_curves_by_pair_and_class.png"), p, width = 10, height = 6.5, dpi = 300)
  }

  message("Wrote master table and pair curve-class summaries under: ", tables_dir)
  invisible(list(master = file.path(tables_dir, "joint_best_master_table.tsv")))
}

main <- function(argv = parse_args()) {
  stages <- stage_values(argv$stage %||% "all")
  if ("prepare" %in% stages) build_prepare_outputs(argv)
  if ("embedding" %in% stages) build_embedding_outputs(argv)
  if ("curve" %in% stages) run_curve_classification(argv)
  if ("summary" %in% stages) build_summary_outputs(argv)
  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  main(parse_args())
}
