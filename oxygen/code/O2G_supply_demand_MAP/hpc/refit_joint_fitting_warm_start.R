#!/usr/bin/env Rscript
# Submit constrained O2G joint fitting using warm starts from completed
# separate in vivo and in vitro 500-seed runs.
#
# Run on the HPC login node:
#   Rscript oxygen/code/O2G_supply_demand_MAP/hpc/refit_joint_fitting_warm_start.R
#
# It can also be submitted as a lightweight manager job:
#   sbatch oxygen/code/O2G_supply_demand_MAP/hpc/refit_joint_fitting_warm_start.R

#SBATCH --job-name=o2g_refit_joint_ws
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --output=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/o2g_refit_joint_ws_manager_%j.out
#SBATCH --error=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/o2g_refit_joint_ws_manager_%j.err

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    body <- sub("^--", "", arg)
    pos <- regexpr("=", body, fixed = TRUE)
    if (pos < 0) {
      out[[body]] <- "TRUE"
    } else {
      key <- substr(body, 1L, pos - 1L)
      val <- substr(body, pos + 1L, nchar(body))
      out[[key]] <- val
    }
  }
  out
}

first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

as_chr <- function(x, default = "") {
  val <- first_non_null(x, default)
  if (is.null(val) || !length(val)) return(default)
  val <- as.character(val[[1L]])
  if (!nzchar(trimws(val))) default else val
}

as_bool <- function(x, default = FALSE) {
  val <- first_non_null(x, default)
  if (is.logical(val)) return(isTRUE(val[[1L]]))
  txt <- tolower(trimws(as.character(val[[1L]])))
  if (!nzchar(txt)) return(isTRUE(default))
  if (txt %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (txt %in% c("false", "f", "0", "no", "n")) return(FALSE)
  isTRUE(default)
}

as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(first_non_null(x, default)))
  if (!is.finite(val) || is.na(val)) default else val
}

as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(first_non_null(x, default)))
  if (!is.finite(val) || is.na(val)) default else val
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

normalize_path <- function(path, must_work = FALSE) {
  normalizePath(path.expand(path), mustWork = must_work)
}

resolve_relative_path <- function(path, base_dir) {
  txt <- as_chr(path, "")
  if (!nzchar(txt)) return(txt)
  if (grepl("^(/|~|[A-Za-z]:[/\\\\])", txt)) {
    return(normalize_path(txt, must_work = FALSE))
  }
  normalize_path(file.path(base_dir, txt), must_work = FALSE)
}

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript oxygen/code/O2G_supply_demand_MAP/hpc/refit_joint_fitting_warm_start.R [options]\n\n",
    "Main defaults:\n",
    "  --project_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid\n",
    "  --invivo_run_dir=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invivo_O2G_buffering_500seed\n",
    "  --invitro_run_dir=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invitro_O2G_buffering_500seed\n",
    "  --joint_run_prefix=refit_joint_O2G_buffering_warmstart_bioconstraints_500seed\n\n",
    "Useful options:\n",
    "  --dry_run=TRUE                         Build summaries/config/candidates, print sbatch command only.\n",
    "  --force_extra_results=TRUE             Regenerate seed_summary.tsv before selecting best seeds.\n",
    "  --enable_constraints=FALSE             Submit an unconstrained joint warm-start run.\n",
    "  --total_seeds=500 --array_tasks=500 --seeds_per_task=1\n",
    "  --joint_n_cores=22 --joint_qos=xxlarge --joint_time_limit=12:00:00\n",
    sep = ""
  )
}

stop_missing <- function(label, path, is_dir = FALSE) {
  ok <- if (is_dir) dir.exists(path) else file.exists(path)
  if (!ok) {
    stop("Missing ", label, ": ", path, call. = FALSE)
  }
  invisible(path)
}

run_command <- function(command, args, label) {
  message(label)
  message("  ", command, " ", paste(shQuote(args), collapse = " "))
  status <- system2(command, args = args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop(label, " failed with exit status ", status, call. = FALSE)
  }
  invisible(TRUE)
}

run_rscript <- function(script, args, label) {
  rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) {
    stop("Rscript not found in PATH.", call. = FALSE)
  }
  run_command(rscript, c(script, args), label)
}

read_best_dir <- function(path) {
  stop_missing("best-dir file", path)
  txt <- trimws(readLines(path, warn = FALSE))
  txt <- txt[nzchar(txt)]
  if (!length(txt)) {
    stop("Best-dir file is empty: ", path, call. = FALSE)
  }
  normalize_path(txt[[1L]], must_work = FALSE)
}

validate_best_dir <- function(label, seed_dir) {
  stop_missing(paste(label, "best seed directory"), seed_dir, is_dir = TRUE)
  candidates <- file.path(seed_dir, c("best_params_transformed.tsv", "fit_parameter_stages.tsv"))
  if (!any(file.exists(candidates))) {
    stop(
      label,
      " best seed directory lacks transformed parameters: ",
      seed_dir,
      "\nExpected best_params_transformed.tsv or fit_parameter_stages.tsv.",
      call. = FALSE
    )
  }
  invisible(seed_dir)
}

ensure_extra_results <- function(label,
                                 run_dir,
                                 extra_results_script,
                                 out_dir,
                                 summary_path,
                                 near_thresh,
                                 force_extra_results) {
  if (file.exists(summary_path) && !isTRUE(force_extra_results)) {
    message("Reusing existing ", label, " seed summary: ", summary_path)
    return(invisible(summary_path))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_rscript(
    extra_results_script,
    c(
      paste0("--run_dir=", run_dir),
      paste0("--out_dir=", out_dir),
      paste0("--near_thresh=", near_thresh)
    ),
    paste0("Running ", label, " extra_results")
  )
  stop_missing(paste(label, "seed_summary.tsv"), summary_path)
}

write_manifest <- function(path, values) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tab <- data.frame(
    key = names(values),
    value = vapply(values, as.character, character(1L)),
    stringsAsFactors = FALSE
  )
  utils::write.table(tab, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

write_joint_overlay_config <- function(config_path,
                                       out_path,
                                       config_overrides,
                                       original_config_dir) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to write the joint overlay config.", call. = FALSE)
  }
  cfg <- yaml::read_yaml(config_path)
  if (is.null(cfg)) cfg <- list()
  if (!is.list(cfg) || is.null(names(cfg))) {
    stop("Config must be a named YAML mapping: ", config_path, call. = FALSE)
  }

  path_keys <- c(
    "out_root", "data_dir", "seeds_file", "parameter_table", "parameters",
    "init_params_tsv", "joint_invivo_best_dir", "joint_invitro_best_dir",
    "joint_init_candidates_tsv", "fit_objects_dir", "flow_density_path",
    "invitro_parameter_table", "parameter_table_invitro"
  )
  for (key in intersect(path_keys, names(cfg))) {
    if (is.null(cfg[[key]]) || length(cfg[[key]]) != 1L || is.na(cfg[[key]])) next
    txt <- trimws(as.character(cfg[[key]]))
    if (!nzchar(txt) || tolower(txt) == "null") next
    cfg[[key]] <- resolve_relative_path(txt, original_config_dir)
  }

  for (key in names(config_overrides)) {
    cfg[[key]] <- config_overrides[[key]]
  }

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(cfg, out_path)
  normalize_path(out_path, must_work = FALSE)
}

submit_sbatch <- function(args, dry_run) {
  sbatch <- Sys.which("sbatch")
  if (!nzchar(sbatch) && !isTRUE(dry_run)) {
    stop("sbatch not found in PATH; cannot submit joint fitting.", call. = FALSE)
  }
  message("Submitting joint warm-start array")
  message("  sbatch ", paste(shQuote(args), collapse = " "))
  if (isTRUE(dry_run)) {
    message("DRY_RUN=TRUE; not submitting.")
    return(NA_character_)
  }
  out <- system2(sbatch, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    cat(paste(out, collapse = "\n"), "\n", sep = "")
    stop("sbatch failed with exit status ", status, call. = FALSE)
  }
  job_id <- trimws(out[[length(out)]])
  message("Submitted joint warm-start array job: ", job_id)
  job_id
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  if (as_bool(first_non_null(argv$help, argv$h), FALSE)) {
    usage()
    return(invisible(NULL))
  }

  default_project_root <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid"
  project_root <- normalize_path(as_chr(argv$project_root, default_project_root), must_work = FALSE)
  config_path <- normalize_path(
    as_chr(argv$config, file.path(project_root, "oxygen", "config", "O2G_supply_demand.yaml")),
    must_work = FALSE
  )
  out_root <- normalize_path(
    as_chr(argv$out_root, file.path(project_root, "oxygen", "results")),
    must_work = FALSE
  )
  invivo_run_dir <- normalize_path(
    as_chr(
      argv$invivo_run_dir,
      file.path(out_root, "fit_invivo_O2G_buffering_500seed")
    ),
    must_work = FALSE
  )
  invitro_run_dir <- normalize_path(
    as_chr(
      argv$invitro_run_dir,
      file.path(out_root, "fit_invitro_O2G_buffering_500seed")
    ),
    must_work = FALSE
  )
  joint_run_prefix <- as_chr(
    argv$joint_run_prefix,
    "refit_joint_O2G_buffering_warmstart_bioconstraints_500seed"
  )
  joint_run_dir <- file.path(out_root, joint_run_prefix)

  workflow_root <- normalize_path(file.path(project_root, "oxygen", "code", "O2G_supply_demand_MAP"), must_work = FALSE)
  extra_results_script <- normalize_path(
    as_chr(argv$extra_results_script, file.path(workflow_root, "analysis", "extra_results.R")),
    must_work = FALSE
  )
  invivo_selector_script <- normalize_path(
    as_chr(argv$invivo_selector_script, file.path(workflow_root, "analysis", "select_invivo_best_long_ploidy_seed.R")),
    must_work = FALSE
  )
  invitro_selector_script <- normalize_path(
    as_chr(argv$invitro_selector_script, file.path(workflow_root, "analysis", "select_best_seed_from_summary.R")),
    must_work = FALSE
  )
  build_script <- normalize_path(
    as_chr(argv$build_script, file.path(workflow_root, "analysis", "build_joint_init_candidates.R")),
    must_work = FALSE
  )
  joint_sub_script <- normalize_path(
    as_chr(argv$joint_sub_script, file.path(workflow_root, "hpc", "submit_fit_seed_array_joint_buffering_warmstart.sub")),
    must_work = FALSE
  )

  total_seeds <- as_int(argv$total_seeds, 500L)
  array_tasks <- as_int(argv$array_tasks, 500L)
  seeds_per_task <- as_int(argv$seeds_per_task, 1L)
  joint_n_cores <- as_int(first_non_null(argv$joint_n_cores, argv$n_cores), 22L)
  joint_qos <- as_chr(argv$joint_qos, "xxlarge")
  joint_time_limit <- as_chr(argv$joint_time_limit, "12:00:00")
  r_module <- as_chr(argv$r_module, "R/4.4")
  auto_viz <- as_bool(argv$auto_viz, TRUE)
  glucose <- as_bool(argv$glucose, TRUE)
  dry_run <- as_bool(argv$dry_run, FALSE)
  force_extra_results <- as_bool(argv$force_extra_results, FALSE)
  near_thresh <- as_num(argv$near_thresh, 0.05)

  selection_horizon <- as_num(argv$selection_horizon, 1000)
  selection_threshold_n <- as_num(argv$selection_threshold_N, 44)
  selection_cohort <- as_chr(argv$selection_cohort, "2N")
  selection_dose <- as_chr(argv$selection_dose, "ALL")
  invitro_objective_columns <- as_chr(argv$invitro_objective_columns, "objective,objective_total")

  enable_constraints <- as_bool(argv$enable_constraints, TRUE)
  require_invivo_pred1000 <- as_bool(argv$require_invivo_pred1000_ploidy_gt2, enable_constraints)
  require_invitro_growth <- as_bool(argv$require_invitro_growth_nonnegative, enable_constraints)
  require_invitro_ploidy <- as_bool(argv$require_invitro_ploidy_phenotype, enable_constraints)

  if (total_seeds <= 0L || array_tasks <= 0L || seeds_per_task <= 0L || joint_n_cores <= 0L) {
    stop("total_seeds, array_tasks, seeds_per_task, and joint_n_cores must all be positive.", call. = FALSE)
  }
  if (array_tasks * seeds_per_task != total_seeds) {
    stop(
      "array_tasks * seeds_per_task must equal total_seeds. Got ",
      array_tasks, " * ", seeds_per_task, " != ", total_seeds,
      call. = FALSE
    )
  }

  stop_missing("project root", project_root, is_dir = TRUE)
  stop_missing("config", config_path)
  stop_missing("in vivo run directory", invivo_run_dir, is_dir = TRUE)
  stop_missing("in vitro run directory", invitro_run_dir, is_dir = TRUE)
  for (path in c(extra_results_script, invivo_selector_script, invitro_selector_script, build_script, joint_sub_script)) {
    stop_missing("required helper script", path)
  }
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(joint_run_dir, recursive = TRUE, showWarnings = FALSE)

  invivo_extra_dir <- normalize_path(
    as_chr(argv$invivo_extra_results_dir, file.path(invivo_run_dir, "extra_results")),
    must_work = FALSE
  )
  invitro_extra_dir <- normalize_path(
    as_chr(argv$invitro_extra_results_dir, file.path(invitro_run_dir, "extra_results")),
    must_work = FALSE
  )
  invivo_summary <- normalize_path(
    as_chr(argv$invivo_seed_summary, file.path(invivo_extra_dir, "seed_summary.tsv")),
    must_work = FALSE
  )
  invitro_summary <- normalize_path(
    as_chr(argv$invitro_seed_summary, file.path(invitro_extra_dir, "seed_summary.tsv")),
    must_work = FALSE
  )
  invivo_selection_tsv <- normalize_path(
    as_chr(argv$invivo_selection_tsv, file.path(invivo_run_dir, "best_long_ploidy_gt2_seed.tsv")),
    must_work = FALSE
  )
  invivo_best_dir_file <- normalize_path(
    as_chr(argv$invivo_best_dir_file, file.path(invivo_run_dir, "best_long_ploidy_gt2_seed.dir")),
    must_work = FALSE
  )
  invitro_selection_tsv <- normalize_path(
    as_chr(argv$invitro_selection_tsv, file.path(invitro_run_dir, "best_seed_from_summary.tsv")),
    must_work = FALSE
  )
  invitro_best_dir_file <- normalize_path(
    as_chr(argv$invitro_best_dir_file, file.path(invitro_run_dir, "best_seed_from_summary.dir")),
    must_work = FALSE
  )
  joint_init_candidates_tsv <- normalize_path(
    as_chr(argv$joint_init_candidates_tsv, file.path(joint_run_dir, "joint_init_candidates.tsv")),
    must_work = FALSE
  )
  overlay_config_path <- normalize_path(
    as_chr(argv$joint_overlay_config, file.path(joint_run_dir, "config.refit_joint_warm_start.yaml")),
    must_work = FALSE
  )

  message("Starting O2G joint warm-start refit manager")
  message("  in vivo run dir: ", invivo_run_dir)
  message("  in vitro run dir: ", invitro_run_dir)
  message("  joint run dir: ", joint_run_dir)
  message("  constraints enabled: ", enable_constraints)

  ensure_extra_results(
    label = "in vivo",
    run_dir = invivo_run_dir,
    extra_results_script = extra_results_script,
    out_dir = invivo_extra_dir,
    summary_path = invivo_summary,
    near_thresh = near_thresh,
    force_extra_results = force_extra_results
  )

  run_rscript(
    invivo_selector_script,
    c(
      paste0("--invivo_dir=", invivo_run_dir),
      paste0("--out_tsv=", invivo_selection_tsv),
      paste0("--best_dir_file=", invivo_best_dir_file),
      paste0("--horizon=", selection_horizon),
      paste0("--threshold_N=", selection_threshold_n),
      paste0("--cohort=", selection_cohort),
      paste0("--dose=", selection_dose)
    ),
    "Selecting best in vivo seed"
  )
  invivo_best_dir <- read_best_dir(invivo_best_dir_file)
  validate_best_dir("in vivo", invivo_best_dir)

  ensure_extra_results(
    label = "in vitro",
    run_dir = invitro_run_dir,
    extra_results_script = extra_results_script,
    out_dir = invitro_extra_dir,
    summary_path = invitro_summary,
    near_thresh = near_thresh,
    force_extra_results = force_extra_results
  )

  run_rscript(
    invitro_selector_script,
    c(
      paste0("--run_dir=", invitro_run_dir),
      paste0("--summary_tsv=", invitro_summary),
      paste0("--out_tsv=", invitro_selection_tsv),
      paste0("--best_dir_file=", invitro_best_dir_file),
      paste0("--objective_columns=", invitro_objective_columns),
      "--required_files=best_params_transformed.tsv"
    ),
    "Selecting best in vitro seed"
  )
  invitro_best_dir <- read_best_dir(invitro_best_dir_file)
  validate_best_dir("in vitro", invitro_best_dir)

  run_rscript(
    build_script,
    c(
      paste0("--config=", config_path),
      paste0("--out_tsv=", joint_init_candidates_tsv),
      paste0("--glucose=", if (glucose) "TRUE" else "FALSE"),
      paste0("--invivo_best_dir=", invivo_best_dir),
      paste0("--invitro_best_dir=", invitro_best_dir)
    ),
    "Building joint warm-start candidates"
  )
  stop_missing("joint init candidates TSV", joint_init_candidates_tsv)

  config_overrides <- list(
    out_root = out_root,
    run_prefix = joint_run_prefix,
    append_run_prefix_timestamp = FALSE,
    auto_viz = auto_viz,
    glucose = glucose,
    joint_invivo_best_dir = invivo_best_dir,
    joint_invitro_best_dir = invitro_best_dir,
    joint_init_candidates_tsv = joint_init_candidates_tsv,
    joint_restriction = enable_constraints,
    joint_require_biological_constraints = enable_constraints,
    joint_constraint_penalty = as_num(argv$joint_constraint_penalty, 1e9),
    joint_require_invivo_pred1000_ploidy_gt2 = require_invivo_pred1000,
    joint_invivo_ploidy_horizon_day = as_num(argv$joint_invivo_ploidy_horizon_day, 1000),
    joint_invivo_min_ploidy_fold = as_num(argv$joint_invivo_min_ploidy_fold, 2),
    joint_require_invitro_growth_nonnegative = require_invitro_growth,
    joint_require_invitro_ploidy_phenotype = require_invitro_ploidy,
    joint_invitro_2N_wgd_min_N = as_num(argv$joint_invitro_2N_wgd_min_N, 80),
    joint_invitro_2N_wgd_min_fraction = as_num(argv$joint_invitro_2N_wgd_min_fraction, 0.01),
    joint_invitro_4N_min_chr_drop = as_num(argv$joint_invitro_4N_min_chr_drop, 2)
  )
  overlay_config_path <- write_joint_overlay_config(
    config_path = config_path,
    out_path = overlay_config_path,
    config_overrides = config_overrides,
    original_config_dir = dirname(config_path)
  )
  message("Wrote joint overlay config: ", overlay_config_path)

  manifest_path <- file.path(joint_run_dir, "refit_joint_warm_start_manifest.tsv")
  write_manifest(
    manifest_path,
    c(
      project_root = project_root,
      config_path = config_path,
      overlay_config_path = overlay_config_path,
      invivo_run_dir = invivo_run_dir,
      invitro_run_dir = invitro_run_dir,
      invivo_best_dir = invivo_best_dir,
      invitro_best_dir = invitro_best_dir,
      joint_run_dir = joint_run_dir,
      joint_init_candidates_tsv = joint_init_candidates_tsv,
      total_seeds = total_seeds,
      array_tasks = array_tasks,
      seeds_per_task = seeds_per_task,
      joint_n_cores = joint_n_cores,
      joint_restriction = enable_constraints,
      constraints_enabled = enable_constraints,
      require_invivo_pred1000_ploidy_gt2 = require_invivo_pred1000,
      require_invitro_growth_nonnegative = require_invitro_growth,
      require_invitro_ploidy_phenotype = require_invitro_ploidy
    )
  )
  message("Wrote manifest: ", manifest_path)

  joint_export <- paste(
    c(
      "ALL",
      paste0("PROJECT_ROOT=", project_root),
      paste0("CONFIG_PATH=", overlay_config_path),
      paste0("OUT_ROOT=", out_root),
      paste0("RUN_PREFIX=", joint_run_prefix),
      paste0("TOTAL_SEEDS=", total_seeds),
      paste0("ARRAY_TASKS=", array_tasks),
      paste0("SEEDS_PER_TASK=", seeds_per_task),
      paste0("N_CORES=", joint_n_cores),
      paste0("AUTO_VIZ=", if (auto_viz) "TRUE" else "FALSE"),
      paste0("GLUCOSE=", if (glucose) "TRUE" else "FALSE"),
      paste0("JOINT_INIT_CANDIDATES_TSV=", joint_init_candidates_tsv),
      paste0("R_MODULE=", r_module)
    ),
    collapse = ","
  )
  sbatch_args <- c(
    "--parsable",
    "--job-name=o2gj_refit_WS",
    paste0("--qos=", joint_qos),
    paste0("--time=", joint_time_limit),
    paste0("--cpus-per-task=", joint_n_cores),
    paste0("--array=1-", array_tasks),
    paste0("--output=", file.path(out_root, "o2g_refit_joint_ws_%A_%a.out")),
    paste0("--error=", file.path(out_root, "o2g_refit_joint_ws_%A_%a.err")),
    paste0("--export=", joint_export),
    joint_sub_script
  )
  joint_job_id <- submit_sbatch(sbatch_args, dry_run = dry_run)

  if (!is.na(joint_job_id)) {
    write_manifest(file.path(joint_run_dir, "refit_joint_warm_start_submission.tsv"), c(job_id = joint_job_id))
  }

  invisible(list(
    joint_job_id = joint_job_id,
    joint_run_dir = joint_run_dir,
    overlay_config_path = overlay_config_path,
    joint_init_candidates_tsv = joint_init_candidates_tsv,
    invivo_best_dir = invivo_best_dir,
    invitro_best_dir = invitro_best_dir
  ))
}

main()
