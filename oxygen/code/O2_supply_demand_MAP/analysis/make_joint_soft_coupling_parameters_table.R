#!/usr/bin/env Rscript

# Generate the joint soft-coupling start table from separate seed fits.
#
# The script reads best_params.tsv from one in vivo seed result directory and
# one in vitro seed result directory, then writes a seed-labelled CSV such as:
#
#   oxygen/data/O2_supply_demand/joint_soft_coupling_parameters_table__invivo_seed50__invitro_seed350.csv
#
# Values are written on the optimizer scale expected by the joint backend:
#
#   center = (transformed_in_vivo + transformed_in_vitro) / 2
#   delta  = transformed_in_vivo - transformed_in_vitro
#
# Example:
#
#   Rscript oxygen/code/O2_supply_demand_MAP/analysis/make_joint_soft_coupling_parameters_table.R \
#     --invivo-seed-dir oxygen/results/fit_invivo_O2_buffering_500seed/seed50 \
#     --invitro-seed-dir oxygen/results/fit_invitro_O2_buffering_500seed/seed350 \
#     --seed-label invivo_seed50__invitro_seed350

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!grepl("^--", arg)) {
      i <- i + 1L
      next
    }
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) {
      key <- substr(kv, 1L, pos - 1L)
      val <- substr(kv, pos + 1L, nchar(kv))
      out[[key]] <- val
      i <- i + 1L
    } else if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
      out[[kv]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      out[[kv]] <- TRUE
      i <- i + 1L
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
  val <- as.character(first_non_null(x, default))
  if (!length(val) || !nzchar(val[[1]])) default else val[[1]]
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  }
  normalizePath("oxygen/code/O2_supply_demand_MAP/analysis/make_joint_soft_coupling_parameters_table.R", mustWork = FALSE)
}

project_root <- function() {
  script_dir <- dirname(script_path())
  dirname(dirname(dirname(dirname(script_dir))))
}

sanitize_label <- function(x, fallback = "seed") {
  x <- trimws(as.character(first_non_null(x, fallback)))
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) fallback else x
}

label_from_seed_dir <- function(seed_dir, prefix) {
  paste0(prefix, "_", sanitize_label(basename(normalizePath(seed_dir, mustWork = FALSE))))
}

infer_seed_label <- function(invivo_seed_dir, invitro_seed_dir) {
  paste(
    label_from_seed_dir(invivo_seed_dir, "invivo"),
    label_from_seed_dir(invitro_seed_dir, "invitro"),
    sep = "__"
  )
}

default_output_path <- function(seed_label) {
  filename <- paste0(
    "joint_soft_coupling_parameters_table__",
    sanitize_label(seed_label),
    ".csv"
  )
  file.path(project_root(), "oxygen/data/O2_supply_demand", filename)
}

param_specs <- data.frame(
  parameter = c(
    "O2_crit", "alpha_o2", "gamma_growth", "mu_hp",
    "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
    "buffer_n_exp", "n_O", "lam_max", "p_mis_base",
    "p_wgd", "gamma_mu"
  ),
  param_name = c(
    "log10_O2_crit", "log10_alpha_o2", "gamma_growth", "log10_mu_hp",
    "log10_p_misseg", "log10_k_o_mis", "buffer_smax", "log10_buffer_beta",
    "log10_buffer_n_exp", "n_O", "log10_lam_max", "log10_p_mis_base",
    "log10_p_wgd", "gamma_mu"
  ),
  scale = c(
    "log10", "log10", "identity", "log10",
    "log10", "log10", "identity", "log10",
    "log10", "identity", "log10", "log10",
    "log10", "identity"
  ),
  stringsAsFactors = FALSE
)

default_delta_params <- c(
  "O2_crit", "mu_hp", "p_misseg", "k_o_mis",
  "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O",
  "lam_max", "p_mis_base", "p_wgd", "gamma_mu"
)

usage <- function(default_out) {
  paste0(
    "Usage: Rscript oxygen/code/O2_supply_demand_MAP/analysis/make_joint_soft_coupling_parameters_table.R \\\n",
    "  --invivo-seed-dir <dir> --invitro-seed-dir <dir> [--seed-label <label>] [--out <csv>] [--digits 8]\n\n",
    "Default output: ", default_out, "\n",
    "--seed-label: default is invivo_<basename(invivo-seed-dir)>__invitro_<basename(invitro-seed-dir)>.\n",
    "--delta-params: default, all, none, or comma-separated natural parameter names.\n"
  )
}

read_best_params <- function(seed_dir, label) {
  path <- file.path(seed_dir, "best_params.tsv")
  if (!file.exists(path)) {
    stop("Missing best_params.tsv for ", label, ": ", path, call. = FALSE)
  }
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(c("parameter", "value"), names(tab))
  if (length(missing)) {
    stop(path, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  values <- suppressWarnings(as.numeric(tab$value))
  bad <- !is.finite(values)
  if (any(bad)) {
    stop(
      "Non-numeric value in ", path, " for: ",
      paste(tab$parameter[bad], collapse = ", "),
      call. = FALSE
    )
  }
  out <- values
  names(out) <- trimws(as.character(tab$parameter))
  out[nzchar(names(out))]
}

transform_value <- function(parameter, value, scale) {
  if (!is.finite(value)) {
    stop("Non-finite value for parameter ", parameter, ": ", value, call. = FALSE)
  }
  if (identical(scale, "identity")) return(value)
  if (identical(scale, "log10")) {
    if (value <= 0) {
      stop("Parameter ", parameter, " must be > 0 for log10 transform, got ", value, call. = FALSE)
    }
    return(log10(value))
  }
  stop("Unsupported scale for parameter ", parameter, ": ", scale, call. = FALSE)
}

require_param <- function(params, parameter, label) {
  if (!parameter %in% names(params)) {
    stop("Missing parameter ", parameter, " in ", label, " best_params.tsv", call. = FALSE)
  }
  unname(params[[parameter]])
}

parse_delta_params <- function(spec) {
  spec <- trimws(as_chr(spec, "default"))
  all_params <- param_specs$parameter
  if (identical(spec, "default")) return(default_delta_params)
  if (identical(spec, "all")) return(all_params)
  if (identical(spec, "none")) return(character())
  out <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  out <- out[nzchar(out)]
  unknown <- setdiff(out, all_params)
  if (length(unknown)) {
    stop("Unknown --delta-params entries: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  out
}

build_rows <- function(invivo_params,
                       invitro_params,
                       digits,
                       delta_params,
                       seed_label,
                       invivo_seed_label,
                       invitro_seed_label) {
  rows <- vector("list", 0L)
  row_df <- function(param_name, value, scale) {
    data.frame(
      param_name = param_name,
      value = sprintf(paste0("%.", digits, "f"), value),
      scale = scale,
      seed_label = seed_label,
      invivo_seed_label = invivo_seed_label,
      invitro_seed_label = invitro_seed_label,
      stringsAsFactors = FALSE
    )
  }
  for (i in seq_len(nrow(param_specs))) {
    parameter <- param_specs$parameter[[i]]
    param_name <- param_specs$param_name[[i]]
    scale <- param_specs$scale[[i]]

    vivo_t <- transform_value(parameter, require_param(invivo_params, parameter, "in vivo"), scale)
    vitro_t <- transform_value(parameter, require_param(invitro_params, parameter, "in vitro"), scale)
    center <- (vivo_t + vitro_t) / 2
    delta <- vivo_t - vitro_t

    rows[[length(rows) + 1L]] <- row_df(param_name, center, scale)
    if (parameter %in% delta_params) {
      rows[[length(rows) + 1L]] <- row_df(paste0("delta__", param_name), delta, scale)
    }
  }
  do.call(rbind, rows)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  default_out <- file.path(
    project_root(),
    "oxygen/data/O2_supply_demand/joint_soft_coupling_parameters_table__<seed_label>.csv"
  )

  if (isTRUE(argv$help) || isTRUE(argv$h)) {
    cat(usage(default_out))
    quit(status = 0L)
  }

  invivo_seed_dir <- as_chr(argv[["invivo-seed-dir"]])
  invitro_seed_dir <- as_chr(argv[["invitro-seed-dir"]])
  if (!nzchar(invivo_seed_dir) || !nzchar(invitro_seed_dir)) {
    stop(usage(default_out), call. = FALSE)
  }
  invivo_seed_dir <- normalizePath(invivo_seed_dir, mustWork = FALSE)
  invitro_seed_dir <- normalizePath(invitro_seed_dir, mustWork = FALSE)
  if (!dir.exists(invivo_seed_dir)) {
    stop("invivo seed dir does not exist: ", invivo_seed_dir, call. = FALSE)
  }
  if (!dir.exists(invitro_seed_dir)) {
    stop("invitro seed dir does not exist: ", invitro_seed_dir, call. = FALSE)
  }

  invivo_seed_label <- sanitize_label(as_chr(first_non_null(argv[["invivo-label"]], argv[["invivo_label"]]), basename(invivo_seed_dir)))
  invitro_seed_label <- sanitize_label(as_chr(first_non_null(argv[["invitro-label"]], argv[["invitro_label"]]), basename(invitro_seed_dir)))
  seed_label <- sanitize_label(as_chr(
    first_non_null(argv[["seed-label"]], argv[["seed_label"]], argv[["joint-warmup-seed-label"]], argv[["joint_warmup_seed_label"]]),
    infer_seed_label(invivo_seed_dir, invitro_seed_dir)
  ))

  out_arg <- as_chr(first_non_null(argv$out, argv[["out-file"]], argv[["out_file"]]), "")
  out_dir <- as_chr(first_non_null(argv[["out-dir"]], argv[["out_dir"]]), "")
  if (!nzchar(out_arg)) {
    out <- if (nzchar(out_dir)) {
      file.path(normalizePath(out_dir, mustWork = FALSE), basename(default_output_path(seed_label)))
    } else {
      default_output_path(seed_label)
    }
  } else if (dir.exists(out_arg)) {
    out <- file.path(normalizePath(out_arg, mustWork = FALSE), basename(default_output_path(seed_label)))
  } else {
    out <- out_arg
  }
  out <- normalizePath(out, mustWork = FALSE)
  digits <- suppressWarnings(as.integer(as_chr(argv$digits, "8")))
  if (!is.finite(digits) || digits < 0L) {
    stop("--digits must be a non-negative integer.", call. = FALSE)
  }
  delta_params <- parse_delta_params(argv[["delta-params"]])

  invivo_params <- read_best_params(invivo_seed_dir, "in vivo")
  invitro_params <- read_best_params(invitro_seed_dir, "in vitro")
  rows <- build_rows(
    invivo_params,
    invitro_params,
    digits,
    delta_params,
    seed_label,
    invivo_seed_label,
    invitro_seed_label
  )

  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, file = out, quote = FALSE, row.names = FALSE)

  cat("wrote=", out, "\n", sep = "")
  cat("seed_label=", seed_label, "\n", sep = "")
  cat("invivo_seed_label=", invivo_seed_label, "\n", sep = "")
  cat("invitro_seed_label=", invitro_seed_label, "\n", sep = "")
  cat("invivo_seed_dir=", invivo_seed_dir, "\n", sep = "")
  cat("invitro_seed_dir=", invitro_seed_dir, "\n", sep = "")
  cat("n_rows=", nrow(rows), "\n", sep = "")
}

if (identical(environment(), globalenv())) {
  main()
}
