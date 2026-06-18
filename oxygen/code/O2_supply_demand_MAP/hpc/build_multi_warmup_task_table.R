#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) {
      out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv))
    } else {
      out[[kv]] <- TRUE
    }
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x
}

as_chr <- function(x, default = "") as.character((x %||% default)[[1]])

as_int <- function(x, default) {
  out <- suppressWarnings(as.integer(as_chr(x, as.character(default))))
  if (!is.finite(out) || out <= 0L) as.integer(default) else out
}

as_bool <- function(x, default = FALSE) {
  val <- tolower(trimws(as_chr(x, if (default) "TRUE" else "FALSE")))
  val %in% c("true", "t", "1", "yes", "y", "on")
}

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript build_multi_warmup_task_table.R --multi_warmup_root=DIR [options]\n\n",
    "Options:\n",
    "  --multi_warmup_root=DIR        Root containing multi_warmup_manifest.tsv.\n",
    "  --manifest=FILE               Default: <root>/multi_warmup_manifest.tsv.\n",
    "  --jobs=FILE                   Default: <root>/multi_warmup_jobs.tsv when present.\n",
    "  --out=FILE                    Default: <root>/multi_warmup_tasks.tsv.\n",
    "  --project_root=DIR            Default: inferred from this script path or task settings.\n",
    "  --seeds_per_pair=N            Joint seeds to run for each warm-up pair.\n",
    "  --total_seeds=N               Legacy alias for --seeds_per_pair.\n",
    "  --array_tasks=N               Override pair-array task count.\n",
    "  --seeds_per_task=N            Override seeds per task.\n",
    "  --config_path=FILE            Joint config path written to every task row.\n",
    "  --runner_script=FILE          Joint runner path written to every task row.\n",
    "  --parameter_table=FILE        In vitro parameter table written to every task row.\n",
    "  --fit_objects_dir=DIR         Fit objects directory written to every task row.\n",
    "  --flow_density_path=FILE      Flow-density grid written to every task row.\n",
    "  --joint_n_cores=N             Joint runner core count written to every task row.\n",
    "  --itermax=N, --de_reltol=NUM, --de_steptol=N, --NP=N, --auto_viz=TRUE|FALSE\n",
    "  --joint_soft_coupling_sigma_default=NUM\n",
    "                                Sigma value written to every task row.\n",
    "  --joint_soft_coupling_huber_k=NUM\n",
    "                                Huber k value written to every task row.\n",
    "  --order=round_robin|pair_major Default: round_robin.\n",
    "  --refresh_status=TRUE|FALSE   Default: TRUE.\n",
    "  --log_root=DIR                Default: <project_root>/oxygen/results/log.\n",
    "  --help\n\n",
    "The output has one row per warm-up pair x joint seed and can be consumed by\n",
    "submit_multi_warmup_task_table.sh.\n",
    sep = ""
  )
}

read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Missing file: ", path, call. = FALSE)
    return(NULL)
  }
  utils::read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    na.strings = character()
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

abs_path <- function(path, mustWork = FALSE) {
  normalizePath(path, mustWork = mustWork)
}

script_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) return("")
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
  # hpc/build_multi_warmup_task_table.R -> O2_supply_demand_MAP -> code -> oxygen -> project
  normalizePath(file.path(dirname(script_path), "../../../.."), mustWork = FALSE)
}

kv_map <- function(path) {
  tab <- read_tsv(path, required = FALSE)
  if (is.null(tab) || !all(c("key", "value") %in% names(tab))) {
    return(setNames(character(), character()))
  }
  vals <- as.character(tab$value)
  vals[is.na(vals)] <- ""
  setNames(vals, as.character(tab$key))
}

first_existing_pair_manifest <- function(root, joint_run_prefix) {
  paths <- file.path(root, joint_run_prefix, "fit_array_manifest.tsv")
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1]] else ""
}

setting <- function(settings, name, default = "") {
  if (is.null(settings) || !length(settings) || !(name %in% names(settings))) return(default)
  val <- settings[[name]]
  if (is.null(val) || is.na(val) || !nzchar(val)) default else val
}

arg_value <- function(argv, names) {
  for (name in names) {
    if (!is.null(argv[[name]]) && length(argv[[name]]) && !is.na(argv[[name]][[1]]) && nzchar(as.character(argv[[name]][[1]]))) {
      return(as.character(argv[[name]][[1]]))
    }
  }
  NULL
}

apply_cli_settings <- function(settings, argv) {
  set_from_arg <- function(key, names) {
    val <- arg_value(argv, names)
    if (!is.null(val)) settings[[key]] <<- val
  }
  set_from_arg("config_path", "config_path")
  set_from_arg("runner_script", "runner_script")
  set_from_arg("invitro_parameter_table", c("parameter_table", "invitro_parameter_table", "parameter_table_invitro"))
  set_from_arg("fit_objects_dir", "fit_objects_dir")
  set_from_arg("flow_density_path", "flow_density_path")
  set_from_arg("flow_density_enabled", "flow_density_enabled")
  set_from_arg("n_cores", c("joint_n_cores", "n_cores"))
  set_from_arg("r_module", "r_module")
  set_from_arg("auto_viz", "auto_viz")
  set_from_arg("joint_fitting_mode", "joint_fitting_mode")
  set_from_arg("itermax", "itermax")
  set_from_arg("de_reltol", "de_reltol")
  set_from_arg("de_steptol", "de_steptol")
  set_from_arg("NP", c("NP", "np"))
  settings
}

spec_end <- function(spec) {
  spec <- as.character(spec)
  if (!length(spec) || is.na(spec) || !nzchar(spec)) return(NA_integer_)
  spec <- sub("%.*$", "", spec)
  pieces <- unlist(strsplit(spec, ",", fixed = TRUE), use.names = FALSE)
  ends <- suppressWarnings(as.integer(sub("^.*-", "", pieces)))
  ends <- ends[is.finite(ends)]
  if (!length(ends)) NA_integer_ else max(ends)
}

require_columns <- function(tab, cols, label) {
  missing <- setdiff(cols, names(tab))
  if (length(missing)) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

make_task_status <- function(seed_dir, fit_summary_path, best_params_path, refresh_status = TRUE) {
  if (!refresh_status) return("")
  if (file.exists(fit_summary_path) && file.exists(best_params_path)) return("complete")
  if (dir.exists(seed_dir)) return("seed_dir_present_incomplete_or_running")
  "not_started"
}

row_value <- function(row, col, default = "") {
  if (!col %in% names(row)) return(default)
  val <- row[[col]][[1]]
  if (is.na(val)) default else as.character(val)
}

build_tasks <- function(manifest,
                        jobs,
                        root,
                        project_root,
                        settings,
                        total_seeds,
                        array_tasks,
                        seeds_per_task,
                        order,
                        refresh_status,
                        log_root,
                        joint_soft_coupling_sigma_default,
                        joint_soft_coupling_huber_k,
                        joint_warmup_sigmaN) {
  total_pairs <- nrow(manifest)
  pair_indices <- seq_len(total_pairs)
  seed_indices <- seq_len(total_seeds)
  if (identical(order, "pair_major")) {
    grid <- expand.grid(pair_index = pair_indices, joint_seed = seed_indices)
    grid <- grid[order(grid$pair_index, grid$joint_seed), , drop = FALSE]
  } else {
    grid <- expand.grid(pair_index = pair_indices, joint_seed = seed_indices)
    grid <- grid[order(grid$joint_seed, grid$pair_index), , drop = FALSE]
  }
  row.names(grid) <- NULL

  array_script <- file.path(project_root, "oxygen/code/O2_supply_demand_MAP/hpc/submit_fit_seed_array_joint_buffering.sub")
  postprocess_script <- file.path(project_root, "oxygen/code/O2_supply_demand_MAP/hpc/postprocess_extra_results.sh")
  runner_script_default <- file.path(project_root, "oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh")
  config_default <- file.path(project_root, "oxygen/config/O2_supply_demand.yaml")
  parameter_default <- file.path(project_root, "oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv")
  fit_objects_default <- file.path(project_root, "oxygen/ploidyOxygen/data/fit_objects")
  flow_density_default <- file.path(project_root, "oxygen/data/g0g1_ploidy_density_grid.csv")

  job_match <- function(label) {
    if (is.null(jobs) || !nrow(jobs) || !"warmup_label" %in% names(jobs)) return(NULL)
    hit <- match(label, jobs$warmup_label)
    if (!is.finite(hit)) return(NULL)
    jobs[hit, , drop = FALSE]
  }

  rows <- vector("list", nrow(grid))
  for (idx in seq_len(nrow(grid))) {
    pair_i <- grid$pair_index[[idx]]
    seed <- grid$joint_seed[[idx]]
    m <- manifest[pair_i, , drop = FALSE]
    warmup_label <- row_value(m, "warmup_label")
    run_prefix <- row_value(m, "joint_run_prefix", paste0("fit_joint_", warmup_label))
    joint_run_dir <- file.path(root, run_prefix)
    output_seed_dir <- file.path(joint_run_dir, paste0("seed", seed))
    fit_summary_path <- file.path(output_seed_dir, "fit_summary.tsv")
    best_params_path <- file.path(output_seed_dir, "best_params.tsv")
    pair_array_task_id <- ceiling(seed / seeds_per_task)
    pair_seed_start <- (pair_array_task_id - 1L) * seeds_per_task + 1L
    pair_seed_end <- pair_array_task_id * seeds_per_task
    pair_job <- job_match(warmup_label)

    pair_job_id <- if (!is.null(pair_job) && "job_id" %in% names(pair_job)) as.character(pair_job$job_id[[1]]) else ""
    post_id <- if (!is.null(pair_job) && "postprocess_job_id" %in% names(pair_job)) as.character(pair_job$postprocess_job_id[[1]]) else ""
    qos <- if (!is.null(pair_job) && "qos" %in% names(pair_job)) as.character(pair_job$qos[[1]]) else "xxlarge"
    walltime <- if (!is.null(pair_job) && "walltime" %in% names(pair_job)) as.character(pair_job$walltime[[1]]) else "12:00:00"

    rows[[idx]] <- data.frame(
      global_task_id = idx,
      recommended_sbatch_array_index = idx,
      pair_index = pair_i,
      total_pairs = total_pairs,
      warmup_label = warmup_label,
      phase = row_value(m, "phase"),
      invivo_family = row_value(m, "invivo_family"),
      invivo_rank = row_value(m, "invivo_rank"),
      invivo_seed = row_value(m, "invivo_seed"),
      invivo_seed_dir = row_value(m, "invivo_seed_dir"),
      invitro_family = row_value(m, "invitro_family"),
      invitro_rank = row_value(m, "invitro_rank"),
      invitro_seed = row_value(m, "invitro_seed"),
      invitro_seed_dir = row_value(m, "invitro_seed_dir"),
      selection_reason = row_value(m, "selection_reason"),
      joint_seed = seed,
      joint_seed_label = paste0("seed", seed),
      pair_array_task_id = pair_array_task_id,
      pair_seed_start = pair_seed_start,
      pair_seed_end = pair_seed_end,
      seeds_csv = as.character(seed),
      joint_run_prefix = run_prefix,
      joint_run_dir = joint_run_dir,
      out_root = root,
      output_seed_dir = output_seed_dir,
      fit_summary_path = fit_summary_path,
      best_params_path = best_params_path,
      joint_soft_coupling_parameters_table = row_value(m, "joint_soft_coupling_parameters_table"),
      config_path = setting(settings, "config_path", config_default),
      runner_script = setting(settings, "runner_script", runner_script_default),
      array_script = array_script,
      postprocess_script = postprocess_script,
      parameter_table = setting(settings, "invitro_parameter_table", parameter_default),
      fit_objects_dir = setting(settings, "fit_objects_dir", fit_objects_default),
      flow_density_path = setting(settings, "flow_density_path", flow_density_default),
      flow_density_enabled = setting(settings, "flow_density_enabled", "TRUE"),
      total_seeds = total_seeds,
      array_tasks = array_tasks,
      seeds_per_task = seeds_per_task,
      n_cores = setting(settings, "n_cores", "22"),
      r_module = setting(settings, "r_module", "R/4.4"),
      auto_viz = setting(settings, "auto_viz", "TRUE"),
      joint_fitting_mode = setting(settings, "joint_fitting_mode", "DIRECT"),
      itermax = setting(settings, "itermax", "500"),
      de_reltol = setting(settings, "de_reltol", "1e-4"),
      de_steptol = setting(settings, "de_steptol", "25"),
      NP = setting(settings, "NP", "80"),
      joint_warmup_enable = "TRUE",
      joint_warmup_seed_label = warmup_label,
      joint_warmup_invivo_seed_dir = row_value(m, "invivo_seed_dir"),
      joint_warmup_invitro_seed_dir = row_value(m, "invitro_seed_dir"),
      joint_warmup_sigmaN = joint_warmup_sigmaN,
      joint_soft_coupling_sigma_default = joint_soft_coupling_sigma_default,
      joint_soft_coupling_huber_k = joint_soft_coupling_huber_k,
      qos = qos,
      walltime = walltime,
      previous_pair_job_id = pair_job_id,
      previous_postprocess_job_id = post_id,
      slurm_output_pattern = file.path(log_root, paste0("o2mw_global_%A_", idx, ".out")),
      slurm_error_pattern = file.path(log_root, paste0("o2mw_global_%A_", idx, ".err")),
      seed_dir_exists = if (refresh_status) dir.exists(output_seed_dir) else NA,
      fit_summary_exists = if (refresh_status) file.exists(fit_summary_path) else NA,
      best_params_exists = if (refresh_status) file.exists(best_params_path) else NA,
      task_status = make_task_status(output_seed_dir, fit_summary_path, best_params_path, refresh_status),
      task_key = paste(warmup_label, paste0("seed", seed), sep = "__"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  do.call(rbind, rows)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  if (isTRUE(argv$help) || isTRUE(argv$h)) {
    usage()
    return(invisible(NULL))
  }

  root <- as_chr(argv$multi_warmup_root, as_chr(argv$root))
  if (!nzchar(root)) {
    usage()
    stop("--multi_warmup_root is required.", call. = FALSE)
  }
  root <- abs_path(root, mustWork = TRUE)

  manifest_path <- abs_path(as_chr(argv$manifest, file.path(root, "multi_warmup_manifest.tsv")), mustWork = TRUE)
  jobs_path <- as_chr(argv$jobs, file.path(root, "multi_warmup_jobs.tsv"))
  jobs <- read_tsv(jobs_path, required = FALSE)
  manifest <- read_tsv(manifest_path, required = TRUE)
  require_columns(
    manifest,
    c(
      "warmup_label", "invivo_seed_dir", "invitro_seed_dir",
      "joint_run_prefix", "joint_soft_coupling_parameters_table"
    ),
    "multi_warmup_manifest.tsv"
  )
  if (!nrow(manifest)) stop("Manifest has no rows: ", manifest_path, call. = FALSE)

  first_manifest <- first_existing_pair_manifest(root, manifest$joint_run_prefix)
  settings <- kv_map(first_manifest)
  settings <- apply_cli_settings(settings, argv)

  inferred_project_root <- setting(settings, "project_root", script_project_root())
  project_root <- abs_path(as_chr(argv$project_root, inferred_project_root), mustWork = FALSE)

  spec_n <- NA_integer_
  if (!is.null(jobs) && nrow(jobs) && "array_spec" %in% names(jobs)) {
    spec_n <- spec_end(jobs$array_spec[[1]])
  }
  total_default <- if (is.finite(spec_n)) spec_n else as_int(setting(settings, "total_seeds", ""), 500L)
  seeds_per_pair <- as_int(argv$seeds_per_pair, as_int(argv$joint_seeds_per_pair, NA_integer_))
  total_seeds <- as_int(argv$total_seeds, as_int(setting(settings, "total_seeds", ""), total_default))
  if (is.finite(seeds_per_pair) && seeds_per_pair > 0L) {
    total_seeds <- seeds_per_pair
  }
  array_tasks <- as_int(argv$array_tasks, as_int(setting(settings, "array_tasks", ""), total_seeds))
  seeds_per_task <- as_int(argv$seeds_per_task, as_int(setting(settings, "seeds_per_task", ""), 1L))
  if (is.finite(seeds_per_pair) && seeds_per_pair > 0L && is.null(argv$array_tasks) && is.null(argv$seeds_per_task)) {
    array_tasks <- total_seeds
    seeds_per_task <- 1L
  }
  if (array_tasks * seeds_per_task != total_seeds) {
    stop(
      "array_tasks * seeds_per_task must equal total_seeds. Got array_tasks=",
      array_tasks, ", seeds_per_task=", seeds_per_task, ", total_seeds=", total_seeds,
      call. = FALSE
    )
  }

  order <- as_chr(argv$order, "round_robin")
  if (!order %in% c("round_robin", "pair_major")) {
    stop("--order must be round_robin or pair_major.", call. = FALSE)
  }
  refresh_status <- as_bool(argv$refresh_status, TRUE)
  log_root <- abs_path(as_chr(argv$log_root, file.path(project_root, "oxygen/results/log")), mustWork = FALSE)
  out_path <- abs_path(as_chr(argv$out, file.path(root, "multi_warmup_tasks.tsv")), mustWork = FALSE)
  sigma_default <- as_chr(
    argv$joint_soft_coupling_sigma_default,
    setting(settings, "joint_soft_coupling_sigma_default", "")
  )
  huber_k <- as_chr(
    argv$joint_soft_coupling_huber_k,
    setting(settings, "joint_soft_coupling_huber_k", "")
  )
  warmup_sigmaN <- as_chr(
    argv$joint_warmup_sigmaN,
    setting(settings, "joint_warmup_sigmaN", "")
  )

  tasks <- build_tasks(
    manifest = manifest,
    jobs = jobs,
    root = root,
    project_root = project_root,
    settings = settings,
    total_seeds = total_seeds,
    array_tasks = array_tasks,
    seeds_per_task = seeds_per_task,
    order = order,
    refresh_status = refresh_status,
    log_root = log_root,
    joint_soft_coupling_sigma_default = sigma_default,
    joint_soft_coupling_huber_k = huber_k,
    joint_warmup_sigmaN = warmup_sigmaN
  )

  write_tsv(tasks, out_path)
  message("Wrote ", nrow(tasks), " tasks to ", out_path)
  if (refresh_status && "task_status" %in% names(tasks)) {
    print(table(tasks$task_status, useNA = "ifany"))
  }
}

if (identical(environment(), globalenv())) {
  main()
}
