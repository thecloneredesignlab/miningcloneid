#!/usr/bin/env Rscript

.o2g_joint_init_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2g_joint_init_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = FALSE)

source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2g_joint_init_script_dir)

parse_args <- o2sd_parse_args
as_bool <- o2sd_as_bool
.first_non_null_local <- o2sd_first_non_null
clip <- o2sd_clip

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript build_joint_init_candidates.R \\\n",
    "    --config=/path/to/O2G_supply_demand.yaml \\\n",
    "    --out_tsv=/path/to/joint_init_candidates.tsv [joint options]\n\n",
    "Best seed directories are read from config keys joint_invivo_best_dir and\n",
    "joint_invitro_best_dir unless --invivo_best_dir/--invitro_best_dir are provided.\n\n",
    "Key joint options passed through to the current joint context:\n",
    "  --glucose=TRUE|FALSE\n",
    "  --data_dir=PATH\n",
    "  --fit_objects_dir=PATH\n",
    "  --flow_density_path=PATH\n",
    sep = ""
  )
}

trim_cli_scalar <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  txt <- trimws(as.character(x[[1]]))
  if (!nzchar(txt)) return(NULL)
  txt
}

require_arg <- function(argv, key) {
  val <- trim_cli_scalar(argv[[key]])
  if (is.null(val)) {
    stop("Missing required argument --", key, call. = FALSE)
  }
  val
}

resolve_best_dir <- function(path, config_path, from_cli = FALSE) {
  path <- trim_cli_scalar(path)
  if (is.null(path)) return(NULL)
  if (grepl("^(/|~)", path)) {
    return(normalizePath(path, mustWork = FALSE))
  }
  base_dir <- if (isTRUE(from_cli)) {
    getwd()
  } else {
    dirname(normalizePath(config_path, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, path), mustWork = FALSE)
}

sigmoid_safe <- function(log10_prob, label) {
  p <- 10^as.numeric(log10_prob)
  if (!is.finite(p)) {
    stop("Cannot convert non-finite legacy probability transform for ", label, call. = FALSE)
  }
  qlogis(clip(p, 1e-12, 1 - 1e-12))
}

augment_legacy_transforms <- function(vals) {
  if (!"logit_p_misseg" %in% names(vals) && "log10_p_misseg" %in% names(vals)) {
    vals[["logit_p_misseg"]] <- sigmoid_safe(vals[["log10_p_misseg"]], "p_misseg")
  }
  if (!"logit_p_wgd" %in% names(vals) && "log10_p_wgd" %in% names(vals)) {
    vals[["logit_p_wgd"]] <- sigmoid_safe(vals[["log10_p_wgd"]], "p_wgd")
  }
  if (!"logit_p_wgd_max" %in% names(vals) && "log10_p_wgd_max" %in% names(vals)) {
    vals[["logit_p_wgd_max"]] <- sigmoid_safe(vals[["log10_p_wgd_max"]], "p_wgd_max")
  }
  if (!"delta_lam" %in% names(vals) &&
      all(c("log10_lam_min", "log10_lam_max") %in% names(vals))) {
    lam_min <- 10^as.numeric(vals[["log10_lam_min"]])
    lam_max <- 10^as.numeric(vals[["log10_lam_max"]])
    lam_gap <- lam_max - lam_min
    if (is.finite(lam_gap) && lam_gap > 0) {
      vals[["delta_lam"]] <- log(lam_gap)
    }
  }
  if (!"log10_lam_max" %in% names(vals) &&
      all(c("log10_lam_min", "delta_lam") %in% names(vals))) {
    lam_min <- 10^as.numeric(vals[["log10_lam_min"]])
    lam_max <- lam_min + exp(as.numeric(vals[["delta_lam"]]))
    if (is.finite(lam_max) && lam_max > 0) {
      vals[["log10_lam_max"]] <- log10(lam_max)
    }
  }
  vals
}

read_best_transformed <- function(fit_dir, label) {
  candidates <- c(
    file.path(fit_dir, "best_params_transformed.tsv"),
    file.path(fit_dir, "fit_parameter_stages.tsv")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop(
      label,
      " transformed best-parameter table not found. Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  path <- hit[[1]]
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("transformed_parameter", "transformed_value") %in% names(tab))) {
    stop(
      label,
      " transformed best-parameter table must contain transformed_parameter and transformed_value columns: ",
      path,
      call. = FALSE
    )
  }
  vals <- setNames(
    suppressWarnings(as.numeric(tab$transformed_value)),
    as.character(tab$transformed_parameter)
  )
  vals <- vals[nzchar(names(vals))]
  vals <- vals[is.finite(vals)]
  augment_legacy_transforms(vals)
}

value_from_unprefixed <- function(vals, parameter) {
  raw <- sub("^ivt__", "", parameter)
  if (raw %in% names(vals)) vals[[raw]] else NA_real_
}

make_candidate_row <- function(candidate, vals, source, ctx, input_vals = vals) {
  full_names <- names(ctx$init)
  input_vals <- as.numeric(input_vals)
  names(input_vals) <- full_names
  vals <- as.numeric(vals)
  names(vals) <- full_names
  clipped <- clip(vals, ctx$lower, ctx$upper)
  data.frame(
    candidate = candidate,
    parameter = full_names,
    transformed_value = as.numeric(clipped),
    transformed_value_input = as.numeric(input_vals),
    lower = as.numeric(ctx$lower),
    upper = as.numeric(ctx$upper),
    clipped = as.logical(clipped != vals),
    source = source,
    stringsAsFactors = FALSE
  )
}

build_candidate_invivo_best <- function(ctx, invivo_vals, invitro_vals) {
  vals <- ctx$init
  source <- rep("joint_default", length(vals))
  names(source) <- names(vals)
  for (nm in names(vals)) {
    if (startsWith(nm, "ivt__")) {
      v <- value_from_unprefixed(invitro_vals, nm)
      if (is.finite(v)) {
        vals[[nm]] <- v
        source[[nm]] <- "invitro_best"
      }
    } else if (nm %in% names(invivo_vals) && is.finite(invivo_vals[[nm]])) {
      vals[[nm]] <- invivo_vals[[nm]]
      source[[nm]] <- "invivo_best"
    }
  }
  make_candidate_row("invivo_best_mapped", vals, source, ctx)
}

build_candidate_invitro_best <- function(ctx, invivo_vals, invitro_vals) {
  vals <- ctx$init
  source <- rep("joint_default", length(vals))
  names(source) <- names(vals)
  for (nm in names(vals)) {
    if (startsWith(nm, "ivt__")) {
      v <- value_from_unprefixed(invitro_vals, nm)
      if (is.finite(v)) {
        vals[[nm]] <- v
        source[[nm]] <- "invitro_best"
      }
    } else if (nm %in% names(invitro_vals) && is.finite(invitro_vals[[nm]])) {
      vals[[nm]] <- invitro_vals[[nm]]
      source[[nm]] <- "invitro_best"
    } else if (nm %in% names(invivo_vals) && is.finite(invivo_vals[[nm]])) {
      vals[[nm]] <- invivo_vals[[nm]]
      source[[nm]] <- "invivo_best"
    }
  }
  make_candidate_row("invitro_best_mapped", vals, source, ctx)
}

build_candidate_midpoint <- function(ctx, invivo_candidate, invitro_candidate) {
  a <- setNames(invivo_candidate$transformed_value, invivo_candidate$parameter)[names(ctx$init)]
  b <- setNames(invitro_candidate$transformed_value, invitro_candidate$parameter)[names(ctx$init)]
  vals <- (as.numeric(a) + as.numeric(b)) / 2
  names(vals) <- names(ctx$init)
  source <- rep("midpoint_transformed", length(vals))
  names(source) <- names(vals)
  make_candidate_row("midpoint_transformed", vals, source, ctx)
}

build_candidate_default <- function(ctx) {
  vals <- ctx$init
  source <- rep("joint_default", length(vals))
  names(source) <- names(vals)
  make_candidate_row("default_table_init", vals, source, ctx)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  if (isTRUE(as_bool(.first_non_null_local(argv$help, argv$h), FALSE))) {
    usage()
    return(invisible(NULL))
  }

  out_tsv <- normalizePath(require_arg(argv, "out_tsv"), mustWork = FALSE)
  config_path <- trim_cli_scalar(argv$config)
  if (is.null(config_path)) {
    config_path <- file.path(OXYGEN_ROOT, "config", "O2G_supply_demand.yaml")
  }

  joint_env <- new.env(parent = globalenv())
  joint_backend <- file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_fit_joint_backend.R")
  sys.source(joint_backend, envir = joint_env, chdir = TRUE)

  ctx_argv <- argv
  ctx_argv$config <- config_path
  ctx_argv$out_dir <- tempfile("joint_init_context_")
  ctx_argv$seed <- .first_non_null_local(ctx_argv$seed, 1L)
  ctx_argv$joint_init_candidates_tsv <- NULL
  ctx <- joint_env$build_joint_context(ctx_argv)

  invivo_best_from_cli <- !is.null(trim_cli_scalar(argv$invivo_best_dir))
  invitro_best_from_cli <- !is.null(trim_cli_scalar(argv$invitro_best_dir))
  invivo_best_dir <- trim_cli_scalar(.first_non_null_local(
    argv$invivo_best_dir,
    ctx$raw$joint_invivo_best_dir
  ))
  invitro_best_dir <- trim_cli_scalar(.first_non_null_local(
    argv$invitro_best_dir,
    ctx$raw$joint_invitro_best_dir
  ))
  if (is.null(invivo_best_dir)) {
    stop(
      "Missing in vivo best directory. Provide --invivo_best_dir or set ",
      "joint_invivo_best_dir in the config.",
      call. = FALSE
    )
  }
  if (is.null(invitro_best_dir)) {
    stop(
      "Missing in vitro best directory. Provide --invitro_best_dir or set ",
      "joint_invitro_best_dir in the config.",
      call. = FALSE
    )
  }
  invivo_best_dir <- resolve_best_dir(invivo_best_dir, config_path = config_path, from_cli = invivo_best_from_cli)
  invitro_best_dir <- resolve_best_dir(invitro_best_dir, config_path = config_path, from_cli = invitro_best_from_cli)

  invivo_vals <- read_best_transformed(invivo_best_dir, "in vivo")
  invitro_vals <- read_best_transformed(invitro_best_dir, "in vitro")

  invivo_candidate <- build_candidate_invivo_best(ctx, invivo_vals, invitro_vals)
  invitro_candidate <- build_candidate_invitro_best(ctx, invivo_vals, invitro_vals)
  midpoint_candidate <- build_candidate_midpoint(ctx, invivo_candidate, invitro_candidate)
  default_candidate <- build_candidate_default(ctx)

  out <- dplyr::bind_rows(
    invivo_candidate,
    invitro_candidate,
    midpoint_candidate,
    default_candidate
  )

  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  write.table(out, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  clipped <- out[out$clipped, , drop = FALSE]
  clipped_tsv <- sub("\\.tsv$", ".clipped.tsv", out_tsv)
  if (identical(clipped_tsv, out_tsv)) {
    clipped_tsv <- paste0(out_tsv, ".clipped.tsv")
  }
  write.table(clipped, file = clipped_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Wrote joint init candidates: ", normalizePath(out_tsv, mustWork = FALSE))
  message("Candidates: ", paste(unique(out$candidate), collapse = ", "))
  message("Parameters per candidate: ", length(ctx$init))
  message("Clipped rows: ", nrow(clipped), " (", normalizePath(clipped_tsv, mustWork = FALSE), ")")
  invisible(out_tsv)
}

if (sys.nframe() == 0) {
  main()
}
