#!/usr/bin/env Rscript

# Table-only extra analyses for joint multi-warmup fitting results. This module
# materializes prepare, embedding, and summary tables; plotting and cross-stage
# curve orchestration live in vis/ and runner/ respectively.

local_script_dir <- function() {
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own_file <- frame_files[
    basename(frame_files) == "analyze_warm_up_joint_results.R"
  ]
  if (length(own_file)) return(dirname(own_file[[length(own_file)]]))
  workflow_frames <- frame_files[
    grepl("/O2_supply_demand_MAP/", frame_files, fixed = TRUE)
  ]
  if (length(workflow_frames)) {
    root <- dirname(workflow_frames[[length(workflow_frames)]])
    while (!identical(basename(root), "O2_supply_demand_MAP")) {
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (identical(basename(root), "O2_supply_demand_MAP")) {
      return(file.path(root, "analysis", "multi_warmup"))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[[1]]),
      mustWork = FALSE
    )))
  }
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
WORKFLOW_DIR <- normalizePath(
  Sys.getenv(
    "FIGURE_MODEL_CODE_ROOT",
    unset = file.path(SCRIPT_DIR, "..", "..")
  ),
  mustWork = TRUE
)
REPO_ROOT <- normalizePath(
  Sys.getenv("FIGURE_WORKSPACE_ROOT", unset = WORKFLOW_DIR),
  mustWork = TRUE
)
PL_UTILS <- file.path(WORKFLOW_DIR, "analysis", "parameter_landscape_clustering", "parameter_landscape_analysis_utils.R")
JOINT_BACKEND <- file.path(WORKFLOW_DIR, "util", "o2_supply_demand_map_fit_joint_backend.R")
DENSE_DIR <- file.path(WORKFLOW_DIR, "analysis", "best_fit_parameter_feature", "03_dense-grid_monotonicity_classification")

source(
  file.path(WORKFLOW_DIR, "util", "o2_supply_demand_map_multi_warmup_utils.R"),
  local = environment()
)

is_finite_int <- function(x) {
  is.finite(x) && !is.na(x)
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

seed_number_from_dir <- o2sd_seed_from_dir

order_seed_dirs <- function(seed_dirs) {
  seed_dirs[order(seed_number_from_dir(seed_dirs), basename(seed_dirs))]
}

safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  x
}

default_input_root <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540")
}

default_output_root <- function(input_root) {
  file.path(
    REPO_ROOT, "oxygen", "results", "analysis", "warm_up_joint_fitting_results_extra",
    basename(normalizePath(input_root, mustWork = FALSE))
  )
}

default_invivo_curve_class_table <- function() {
  file.path(
    REPO_ROOT, "oxygen", "results", "analysis", "best_fit_parameter_feature",
    "03_dense-grid_monotonicity_classification", "monotonicity_classification",
    "dense-grid_monotonicity_classification", "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv"
  )
}

legacy_invivo_curve_class_table <- function() {
  file.path(
    REPO_ROOT, "oxygen", "results", "analysis", "monotonicity_classification",
    "dense-grid_monotonicity_classification", "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv"
  )
}

stage_values <- function(x) {
  vals <- trimws(strsplit(as.character(x %||% "all"), ",", fixed = TRUE)[[1]])
  vals <- vals[nzchar(vals)]
  vals <- gsub("-", "_", tolower(vals), fixed = TRUE)
  if ("all" %in% vals) return(c("prepare", "embedding", "curve", "summary"))
  match.arg(vals, c("prepare", "embedding", "curve", "curve_regression", "summary"), several.ok = TRUE)
}

cli_args_from_list <- function(args, stage = NULL) {
  args <- args %||% list()
  if (!is.null(stage)) args$stage <- stage
  keys <- names(args)
  if (is.null(keys)) return(character())
  keys <- keys[nzchar(keys)]
  out <- character(0)
  for (key in keys) {
    val <- args[[key]]
    if (is.null(val) || !length(val) || (length(val) == 1L && is.na(val))) next
    out <- c(out, paste0("--", key, "=", as.character(val[[1]])))
  }
  out
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
  source(PL_UTILS, local = globalenv(), chdir = TRUE)
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

effective_args_file_for_seed <- function(seed_dir) {
  candidates <- c(
    file.path(seed_dir, "run_effective_args.tsv"),
    file.path(dirname(seed_dir), "run_effective_args.tsv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[[1L]], mustWork = FALSE) else NA_character_
}

read_effective_argv <- function(seed_dir, path_map_from, path_map_to) {
  path <- effective_args_file_for_seed(seed_dir)
  if (is.na(path)) stop("Missing run_effective_args.tsv in seed or pair directory: ", seed_dir, call. = FALSE)
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

parse_pair_id_label <- function(pair_id) {
  parse_pair_label(paste0("fit_joint_", as.character(pair_id)))
}

seed_id_from_value <- function(x) {
  x_chr <- trimws(as.character(x))
  out <- ifelse(grepl("^seed[0-9]+$", x_chr), x_chr, ifelse(grepl("^[0-9]+$", x_chr), paste0("seed", x_chr), NA_character_))
  out[!nzchar(x_chr) | is.na(x_chr)] <- NA_character_
  out
}

resolve_invivo_curve_class_table <- function(args, required = FALSE) {
  explicit <- args$invivo_curve_class_table %||%
    args$invivo_curve_class_path %||%
    args$source_invivo_curve_class_table %||%
    args$source_invivo_curve_table
  if (!is.null(explicit) && length(explicit) && !is.na(explicit[[1]]) && nzchar(as.character(explicit[[1]]))) {
    path <- normalizePath(path.expand(as.character(explicit[[1]])), mustWork = FALSE)
    if (!file.exists(path)) stop("Missing --invivo_curve_class_table: ", path, call. = FALSE)
    return(path)
  }
  candidates <- c(default_invivo_curve_class_table(), legacy_invivo_curve_class_table())
  hits <- candidates[file.exists(candidates)]
  if (length(hits)) return(normalizePath(hits[[1L]], mustWork = FALSE))
  if (required) {
    stop(
      "Could not find the source in vivo curve-class table. Pass --invivo_curve_class_table=PATH. Tried: ",
      paste(candidates, collapse = "; "),
      call. = FALSE
    )
  }
  NA_character_
}

fill_missing_pair_metadata <- function(df) {
  if (!"pair_id" %in% names(df)) stop("joint best table is missing pair_id.", call. = FALSE)
  defaults <- list(
    method = NA_character_,
    invivo_warmup_seed = NA_integer_,
    invivo_cluster = NA_character_,
    invivo_subcluster = NA_character_,
    invitro_warmup_seed = NA_integer_
  )
  for (nm in names(defaults)) {
    if (!nm %in% names(df)) df[[nm]] <- defaults[[nm]]
  }
  for (i in seq_len(nrow(df))) {
    parsed <- parse_pair_id_label(df$pair_id[[i]])
    for (nm in names(defaults)) {
      val <- df[[nm]][[i]]
      if ((is.na(val) || !nzchar(as.character(val))) && !is.null(parsed[[nm]]) && !is.na(parsed[[nm]])) {
        df[[nm]][[i]] <- parsed[[nm]]
      }
    }
  }
  df
}

export_pair_invivo_warmup_curve_class <- function(best_df, out_root, args, required = FALSE) {
  if (!nrow(best_df)) return(invisible(NULL))
  curve_path <- resolve_invivo_curve_class_table(args, required = required)
  if (is.na(curve_path) || !nzchar(curve_path)) {
    message("Skipping pair in vivo warm-up curve-class export; source curve-class table was not found.")
    return(invisible(NULL))
  }
  curve_df <- read_tsv(curve_path)
  if (!"seed_id" %in% names(curve_df)) stop("Source curve-class table is missing seed_id: ", curve_path, call. = FALSE)

  best_meta <- fill_missing_pair_metadata(best_df)
  best_meta$invivo_seed_id <- seed_id_from_value(best_meta$invivo_warmup_seed)
  group_cols <- c(
    "pair_id", "method", "invivo_warmup_seed", "invivo_seed_id",
    "invivo_cluster", "invivo_subcluster", "invitro_warmup_seed"
  )
  pair_df <- best_meta[, group_cols, drop = FALSE]
  pair_df$.row_count <- 1L
  pair_summary <- stats::aggregate(.row_count ~ ., pair_df, FUN = sum, na.action = stats::na.pass)
  names(pair_summary)[names(pair_summary) == ".row_count"] <- "n_joint_seeds"

  curve_cols <- intersect(
    c(
      "seed_id", "curve_class", "final_interpretation_class", "monotonicity_reliability",
      "classification_rule_version", "min_spectral_gap", "median_spectral_gap",
      "fraction_o2_gap_below_0p005", "fraction_o2_gap_below_0p01",
      "ploidy_range", "net_ploidy_change", "objective"
    ),
    names(curve_df)
  )
  curve_keep <- curve_df[, curve_cols, drop = FALSE]
  names(curve_keep)[names(curve_keep) == "seed_id"] <- "invivo_seed_id"
  names(curve_keep)[names(curve_keep) == "objective"] <- "objective_original_invivo_seed"

  pair_summary$.order <- seq_len(nrow(pair_summary))
  out <- merge(pair_summary, curve_keep, by = "invivo_seed_id", all.x = TRUE, sort = FALSE)
  out <- out[order(out$.order), , drop = FALSE]
  out$.order <- NULL
  out$source_curve_table <- normalizePath(curve_path, mustWork = FALSE)
  ordered_cols <- intersect(
    c(
      "pair_id", "method", "invivo_warmup_seed", "invivo_seed_id",
      "invivo_cluster", "invivo_subcluster", "invitro_warmup_seed", "n_joint_seeds",
      "curve_class", "final_interpretation_class", "monotonicity_reliability",
      "classification_rule_version", "min_spectral_gap", "median_spectral_gap",
      "fraction_o2_gap_below_0p005", "fraction_o2_gap_below_0p01",
      "ploidy_range", "net_ploidy_change", "objective_original_invivo_seed",
      "source_curve_table"
    ),
    names(out)
  )
  out <- out[, ordered_cols, drop = FALSE]

  tables_dir <- file.path(out_root, "tables")
  tsv_path <- file.path(tables_dir, "joint_pair_invivo_warmup_seed_curve_class.tsv")
  tables_csv_path <- file.path(tables_dir, "joint_pair_invivo_warmup_seed_curve_class.csv")
  root_csv_path <- file.path(out_root, "joint_pair_invivo_warmup_seed_curve_class.csv")
  write_tsv(out, tsv_path)
  write_csv(out, tables_csv_path)
  write_csv(out, root_csv_path)
  message("Wrote pair in vivo warm-up seed curve-class table: ", root_csv_path)
  invisible(list(tsv = tsv_path, tables_csv = tables_csv_path, csv = root_csv_path))
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
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  initial_path <- file.path(tables_dir, "joint_initial_population_invivo_params.tsv")
  best_path <- file.path(tables_dir, "joint_best_invivo_params.tsv")
  manifest_path <- file.path(tables_dir, "joint_seed_manifest.tsv")
  curve_manifest_path <- file.path(tables_dir, "joint_best_curve_seed_manifest.tsv")
  if (!overwrite && all(file.exists(c(initial_path, best_path, manifest_path, curve_manifest_path)))) {
    message("Reusing prepare outputs under: ", out_root)
    return(invisible(list(
      initial = initial_path,
      best = best_path,
      manifest = manifest_path,
      curve_seed_manifest = curve_manifest_path
    )))
  }

  pair_dirs <- discover_pair_dirs(input_root, max_pairs = max_pairs)
  best_vectors <- list()
  seed_records <- list()
  initial_rows <- list()
  best_rows <- list()
  curve_manifest_rows <- list()
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

      cfg_path <- file.path(seed_dir, "fit_config.rds")
      curve_manifest_row <- data.frame(
        seed_id = synthetic_seed_id,
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
        source_seed_dir = normalizePath(seed_dir, mustWork = FALSE),
        parameter_file = normalizePath(file.path(seed_dir, "best_params.tsv"), mustWork = FALSE),
        fit_summary_file = normalizePath(file.path(seed_dir, "fit_summary.tsv"), mustWork = FALSE),
        config_file = if (file.exists(cfg_path)) normalizePath(cfg_path, mustWork = FALSE) else NA_character_,
        run_effective_args_file = effective_args_file_for_seed(seed_dir),
        fit_status = metric_value(metrics, "fit_status", "ok"),
        fit_mode = metric_value(metrics, "fit_mode", "fit_joint"),
        objective = metric_num(metrics, "objective"),
        objective_total = metric_num(metrics, "objective"),
        objective_data = metric_num(metrics, "objective_invivo_data"),
        objective_burden = metric_num(metrics, "objective_invivo_burden"),
        objective_ploidy = metric_num(metrics, "objective_invivo_ploidy"),
        deoptim_stop_reason = metric_value(metrics, "deoptim_stop_reason"),
        stringsAsFactors = FALSE
      )

      list(
        initial = initial_row,
        best = best_row,
        curve_manifest = curve_manifest_row,
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
    curve_manifest_rows <- c(curve_manifest_rows, lapply(seed_results, `[[`, "curve_manifest"))
    seed_records <- c(seed_records, lapply(seed_results, `[[`, "seed_record"))
    synthetic_counter <- synthetic_counter + length(seed_dirs)
  }

  initial_df <- rbind_fill(initial_rows)
  best_df <- rbind_fill(best_rows)
  manifest_df <- rbind_fill(seed_records)
  curve_manifest_df <- rbind_fill(curve_manifest_rows)

  write_tsv(initial_df, initial_path)
  write_tsv(best_df, best_path)
  write_tsv(manifest_df, manifest_path)
  write_tsv(curve_manifest_df, curve_manifest_path)
  export_pair_invivo_warmup_curve_class(best_df, out_root, args, required = FALSE)
  message("Wrote initial population table: ", initial_path, " (", nrow(initial_df), " rows)")
  message("Wrote joint best table: ", best_path, " (", nrow(best_df), " rows)")
  message("Wrote curve seed manifest: ", curve_manifest_path, " (", nrow(curve_manifest_df), " rows)")

  invisible(list(
    initial = initial_path,
    best = best_path,
    manifest = manifest_path,
    curve_seed_manifest = curve_manifest_path
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




build_embedding_outputs <- function(args) {
  load_parameter_landscape_utils()
  input_root <- normalizePath(path.expand(args$input_root %||% default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% default_output_root(input_root)), mustWork = FALSE)
  tables_dir <- file.path(out_root, "tables")
  emb_root_dir <- file.path(out_root, "embeddings")
  dir.create(emb_root_dir, recursive = TRUE, showWarnings = FALSE)

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
  tsne_threads <- as_int(args$tsne_threads, 1L)
  run_full_tsne <- as_bool(args$run_full_tsne, FALSE)
  tsne_initial_sample <- as_int(args$tsne_initial_sample, 100000L)

  for (param_set in param_sets) {
    cfg <- embedding_param_config(param_set)
    emb_dir <- file.path(emb_root_dir, cfg$token)
    dir.create(emb_dir, recursive = TRUE, showWarnings = FALSE)

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
          variance_figure_prefix = NULL
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
          max_iter = tsne_max_iter,
          num_threads = tsne_threads
        )
        coord_cols <- c("tSNE1", "tSNE2")
        table_path <- file.path(emb_dir, paste0("joint_", cfg$token, "_tsne", label_suffix, "_coordinates.tsv"))
      }
      coord_df <- cbind(meta, as.data.frame(emb, check.names = FALSE))
      coord_df$embedding_param_set <- cfg$token
      write_tsv(coord_df, table_path)
      message("Wrote embedding: ", table_path)
    }
  }
  invisible(TRUE)
}





build_summary_outputs <- function(args) {
  input_root <- normalizePath(path.expand(args$input_root %||% default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% default_output_root(input_root)), mustWork = FALSE)
  tables_dir <- file.path(out_root, "tables")
  curve_manifest_path <- file.path(tables_dir, "joint_best_curve_seed_manifest.tsv")
  best_path <- file.path(tables_dir, "joint_best_invivo_params.tsv")
  if (!file.exists(curve_manifest_path) || !file.exists(best_path)) {
    stop("Missing prepare outputs. Run --stage=prepare first.", call. = FALSE)
  }
  best_df <- read_tsv(best_path)
  export_pair_invivo_warmup_curve_class(best_df, out_root, args, required = FALSE)

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

  message("Wrote master table and pair curve-class summaries under: ", tables_dir)
  invisible(list(master = file.path(tables_dir, "joint_best_master_table.tsv")))
}

main <- function(argv = parse_args()) {
  stages <- stage_values(argv$stage %||% "all")
  if ("prepare" %in% stages) build_prepare_outputs(argv)
  if ("embedding" %in% stages) build_embedding_outputs(argv)
  if ("summary" %in% stages) build_summary_outputs(argv)
  invisible(TRUE)
}
