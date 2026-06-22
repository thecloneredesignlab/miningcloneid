#!/usr/bin/env Rscript

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

default_paperfigures_dir <- function() {
  code_dir <- dirname(script_dir())
  oxygen_dir <- dirname(code_dir)
  file.path(oxygen_dir, "results", "PaperFigures")
}

DEFAULT_INPUT_DIR <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed"
DEFAULT_OUTPUT_DIR <- default_paperfigures_dir()

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

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  y <- tolower(trimws(as.character(x[[1]])))
  if (!nzchar(y)) return(default)
  if (y %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (y %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  y <- suppressWarnings(as.integer(x[[1]]))
  if (!is.finite(y) || is.na(y)) default else y
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
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

metric_map <- function(path) {
  tab <- read_tsv(path)
  if (!all(c("metric", "value") %in% names(tab))) {
    stop("Summary file is not a metric/value table: ", path)
  }
  stats::setNames(as.character(tab$value), as.character(tab$metric))
}

metric_value <- function(metrics, keys, default = NA_character_) {
  hit <- keys[keys %in% names(metrics)]
  if (!length(hit)) return(default)
  val <- metrics[[hit[[1]]]]
  if (length(val) == 0L || is.na(val) || !nzchar(trimws(as.character(val)))) default else val
}

seed_from_dir <- function(path) {
  as.integer(sub("^seed", "", basename(path)))
}

list_seed_dirs <- function(input_dir) {
  dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("^seed[0-9]+$", basename(dirs))]
  dirs[order(seed_from_dir(dirs))]
}

read_best_params <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  if (!file.exists(path)) stop("Missing best_params.tsv: ", path)
  tab <- read_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain parameter/value columns: ", path)
  }
  vals <- suppressWarnings(as.numeric(tab$value))
  stats::setNames(vals, as.character(tab$parameter))
}

read_param_table <- function(seed_dir) {
  candidates <- c(
    file.path(seed_dir, "parameter_table.csv"),
    file.path(dirname(seed_dir), "parameter_table.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) stop("Missing parameter_table.csv for ", seed_dir)
  path <- existing[[1L]]
  tab <- read_csv_plain(path)
  req <- c("param_name", "estimate", "init_value", "lower_bound", "upper_bound", "param_prototype")
  missing <- setdiff(req, names(tab))
  if (length(missing)) {
    stop("parameter_table.csv missing columns: ", paste(missing, collapse = ", "), " in ", path)
  }
  tab$estimate <- vapply(tab$estimate, as_bool, logical(1), default = FALSE)
  tab$init_value <- suppressWarnings(as.numeric(tab$init_value))
  tab$lower_bound <- suppressWarnings(as.numeric(tab$lower_bound))
  tab$upper_bound <- suppressWarnings(as.numeric(tab$upper_bound))
  tab
}

read_np <- function(metrics, param_table) {
  np <- as_int(metric_value(metrics, "NP"), NA_integer_)
  if (!is.finite(np) || is.na(np) || np < 1L) {
    np <- max(10L * sum(param_table$estimate), 1L)
    warning("Could not read NP from fit_summary.tsv; using ", np)
  }
  np
}

sample_uniform_box <- function(n, lower, upper) {
  n <- as.integer(n)
  d <- length(lower)
  if (is.na(n) || n <= 0L || d <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  u <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  span <- as.numeric(upper - lower)
  out <- sweep(u, 2L, span, `*`)
  out <- sweep(out, 2L, as.numeric(lower), `+`)
  colnames(out) <- names(lower)
  out
}

sample_truncnorm_box <- function(n, center, lower, upper, sigma_frac = 0.1) {
  n <- as.integer(n)
  d <- length(lower)
  if (is.na(n) || n <= 0L || d <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  sigma_frac <- as.numeric(sigma_frac)
  if (!is.finite(sigma_frac) || sigma_frac <= 0) sigma_frac <- 0.1
  sd_vec <- pmax((as.numeric(upper) - as.numeric(lower)) * sigma_frac, 1e-12)
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  out <- sweep(z, 2L, sd_vec, `*`)
  out <- sweep(out, 2L, as.numeric(center), `+`)
  out <- sweep(out, 2L, as.numeric(lower), pmax)
  out <- sweep(out, 2L, as.numeric(upper), pmin)
  colnames(out) <- names(lower)
  out
}

build_de_initialpop <- function(np,
                                lower,
                                upper,
                                init_use,
                                mode = "uniform",
                                uniform_frac = 0.3,
                                sigma_frac = 0.1) {
  np <- as.integer(np)
  if (is.na(np) || np < 1L) stop("NP must be >= 1")
  mode <- tolower(trimws(as.character(mode)))
  if (!mode %in% c("uniform", "hybrid")) stop("de_init_mode must be uniform or hybrid")

  pop <- sample_uniform_box(np, lower, upper)
  init_vec <- as.numeric(init_use[names(lower)])
  names(init_vec) <- names(lower)
  if (any(!is.finite(init_vec))) stop("Initial vector has missing/non-finite values")
  init_vec <- pmin(pmax(init_vec, lower), upper)
  pop[1L, ] <- init_vec

  rem <- as.integer(np - 1L)
  if (rem <= 0L) return(pop)

  if (identical(mode, "uniform")) {
    pop[2L:np, ] <- sample_uniform_box(rem, lower, upper)
    return(pop)
  }

  uniform_frac <- as.numeric(uniform_frac)
  if (!is.finite(uniform_frac) || uniform_frac < 0 || uniform_frac > 1) {
    stop("de_init_uniform_frac must be in [0,1]")
  }
  n_uniform <- as.integer(round(rem * uniform_frac))
  n_uniform <- max(0L, min(rem, n_uniform))
  n_local <- as.integer(rem - n_uniform)
  cursor <- 2L
  if (n_local > 0L) {
    idx <- cursor:(cursor + n_local - 1L)
    pop[idx, ] <- sample_truncnorm_box(n_local, init_vec, lower, upper, sigma_frac)
    cursor <- cursor + n_local
  }
  if (n_uniform > 0L) {
    idx <- cursor:(cursor + n_uniform - 1L)
    pop[idx, ] <- sample_uniform_box(n_uniform, lower, upper)
  }
  if (rem > 1L) {
    ord <- sample.int(rem, size = rem, replace = FALSE)
    pop[2L:np, ] <- pop[1L + ord, , drop = FALSE]
  }
  pop
}

inverse_transform_column <- function(param_name, x) {
  if (startsWith(param_name, "log10_")) {
    return(10^as.numeric(x))
  }
  as.numeric(x)
}

initial_population_natural <- function(seed_dir, seed, metrics, best_names) {
  param_table <- read_param_table(seed_dir)
  opt_tab <- param_table[param_table$estimate, , drop = FALSE]
  if (!nrow(opt_tab)) stop("No estimated parameters in parameter table for ", seed_dir)

  param_names <- as.character(opt_tab$param_name)
  init <- stats::setNames(as.numeric(opt_tab$init_value), param_names)
  lower <- stats::setNames(as.numeric(opt_tab$lower_bound), param_names)
  upper <- stats::setNames(as.numeric(opt_tab$upper_bound), param_names)
  np <- read_np(metrics, param_table)

  set.seed(as.integer(seed))
  pop_t <- build_de_initialpop(
    np = np,
    lower = lower,
    upper = upper,
    init_use = init,
    mode = "uniform"
  )

  out <- data.frame(matrix(nrow = nrow(pop_t), ncol = 0L))
  for (j in seq_along(param_names)) {
    natural_name <- as.character(opt_tab$param_prototype[[j]])
    out[[natural_name]] <- inverse_transform_column(param_names[[j]], pop_t[, j])
  }

  if ("c_vol_2N_eff_mm3" %in% best_names && "rho_2N" %in% names(out)) {
    out[["c_vol_2N_eff_mm3"]] <- 1 / out[["rho_2N"]]
  }

  best_vals <- read_best_params(seed_dir)
  for (nm in setdiff(best_names, names(out))) {
    if (nm %in% names(best_vals)) {
      out[[nm]] <- as.numeric(best_vals[[nm]])
    } else {
      out[[nm]] <- NA_real_
    }
  }
  out <- out[, best_names, drop = FALSE]
  out$seed <- as.integer(seed)
  out
}

best_params_row <- function(seed_dir, seed, best_names) {
  vals <- read_best_params(seed_dir)
  out <- as.data.frame(as.list(stats::setNames(rep(NA_real_, length(best_names)), best_names)))
  for (nm in intersect(best_names, names(vals))) {
    out[[nm]] <- as.numeric(vals[[nm]])
  }
  metrics <- metric_map(file.path(seed_dir, "fit_summary.tsv"))
  iter_completed <- suppressWarnings(as.integer(metric_value(
    metrics,
    c("optimizer_iter_completed", "deoptim_iter_completed", "iter_completed"),
    default = NA_character_
  )))
  if (!is.finite(iter_completed) || is.na(iter_completed)) {
    iter_completed <- suppressWarnings(as.integer(metric_value(metrics, "itermax", default = NA_character_)))
    warning("Using itermax fallback for DEoptim completed iteration in seed", seed)
  }
  out$deoptim_iter_completed <- iter_completed
  out$seed <- as.integer(seed)
  out
}

write_csv <- function(df, path) {
  utils::write.csv(df, file = path, quote = FALSE, row.names = FALSE)
  message("Wrote ", path, " [", nrow(df), " x ", ncol(df), "]")
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  input_dir <- normalizePath(argv$input_dir %||% DEFAULT_INPUT_DIR, mustWork = FALSE)
  output_dir <- normalizePath(path.expand(argv$output_dir %||% DEFAULT_OUTPUT_DIR), mustWork = FALSE)
  max_seeds <- as_int(argv$max_seeds, NA_integer_)

  if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  seed_dirs <- list_seed_dirs(input_dir)
  if (!length(seed_dirs)) stop("No seed directories found under: ", input_dir)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seed_dirs <- seed_dirs[seq_len(min(length(seed_dirs), max_seeds))]
  }

  best_vectors <- list()
  for (seed_dir in seed_dirs) {
    seed <- seed_from_dir(seed_dir)
    best_vectors[[as.character(seed)]] <- read_best_params(seed_dir)
  }
  best_names <- names(best_vectors[[1L]])
  for (vals in best_vectors[-1L]) {
    best_names <- c(best_names, setdiff(names(vals), best_names))
  }

  message("Found ", length(seed_dirs), " seed directories.")
  message("Using ", length(best_names), " natural parameter columns.")

  best_rows <- vector("list", length(seed_dirs))
  init_rows <- vector("list", length(seed_dirs))
  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    seed <- seed_from_dir(seed_dir)
    summary_path <- file.path(seed_dir, "fit_summary.tsv")
    if (!file.exists(summary_path)) stop("Missing fit_summary.tsv: ", summary_path)
    metrics <- metric_map(summary_path)

    best_rows[[i]] <- best_params_row(seed_dir, seed, best_names)
    init_rows[[i]] <- initial_population_natural(seed_dir, seed, metrics, best_names)

    if (i %% 25L == 0L || i == length(seed_dirs)) {
      message("Processed ", i, "/", length(seed_dirs), " seeds.")
    }
  }

  best_df <- do.call(rbind, best_rows)
  init_df <- do.call(rbind, init_rows)

  best_df <- best_df[, c(best_names, "deoptim_iter_completed", "seed"), drop = FALSE]
  init_df <- init_df[, c(best_names, "seed"), drop = FALSE]

  write_csv(
    init_df,
    file.path(output_dir, "invivo_deoptim_initial_population.csv")
  )
  write_csv(
    best_df,
    file.path(output_dir, "invivo_best_params_by_seed.csv")
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(trimws(as.character(x[[1]])))) y else x
}

main()
