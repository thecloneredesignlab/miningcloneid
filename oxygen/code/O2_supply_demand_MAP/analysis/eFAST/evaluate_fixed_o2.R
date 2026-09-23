#!/usr/bin/env Rscript

# Evaluate the canonical fixed-O2 operator for a SALib FAST sample, preserving row order.
args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) < 2L || !startsWith(arg, "--")) stop("Use --key=value arguments: ", arg)
    out[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
  }
  out
}
opt <- parse_args(args)
required <- c("samples", "metadata", "fit_root", "out")
missing <- required[!vapply(required, function(x) nzchar(opt[[x]] %||% ""), logical(1L))]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
fixed_script <- normalizePath(file.path(script_dir, "..", "..", "simulation", "o2", "fixed_o2",
                                        "run_fixed_o2_simulation.R"), mustWork = TRUE)
source(fixed_script, local = TRUE)
suppressPackageStartupMessages(library(jsonlite))

samples <- read.delim(gzfile(opt$samples), check.names = FALSE)
metadata <- jsonlite::fromJSON(opt$metadata)
oxygen <- as.numeric(metadata$oxygen_pct)
if (nrow(samples) != metadata$n_samples ||
    !identical(as.integer(samples$sample_id), seq_len(nrow(samples)))) {
  stop("FAST sample row order or count is invalid")
}
if (!identical(as.character(metadata$method), "eFAST")) stop("Not an eFAST design")
workers <- as.integer(opt$workers %||% "1")
if (!is.finite(workers) || workers < 1L) stop("workers must be positive")
workers <- min(workers, nrow(samples))
if (workers > 1L && .Platform$OS.type != "unix") stop("Parallel evaluation requires Unix")

seed_dir <- file.path(opt$fit_root, "seed25")
cfg_raw <- readRDS(file.path(seed_dir, "fit_config.rds"))
fit_params <- read.delim(file.path(seed_dir, "best_params.tsv"), check.names = FALSE)
base_params <- as.list(stats::setNames(as.numeric(fit_params$value), fit_params$parameter))
cfg <- prepare_sim_cfg(cfg_raw, list(), 0, prepare_run_params(base_params, "invivo", cfg_raw, 0))
if (cfg$N_MIN != 22L || cfg$N_MAX != 154L || cfg$N_UNIT != 22L ||
    !isTRUE(cfg$O2_growth) || !identical(cfg$ploidy_O2_death, "ploidy_related")) {
  stop("Unexpected in-vivo model configuration")
}

model_env <- globalenv()
evaluate_one <- function(i) {
  p <- base_params
  for (name in setdiff(names(samples), "sample_id")) p[[name]] <- as.numeric(samples[[name]][[i]])
  rp <- prepare_run_params(p, "invivo", cfg, 0)
  rows <- lapply(oxygen, function(o2) {
    result <- fixo2_dominant_attractor_one(as.character(i), rp, model_env, cfg, o2)
    data.frame(sample_id = i, O2_pct = o2,
               dominant_mean_ploidy = result$dominant_mean_ploidy,
               dominant_growth_rate = result$dominant_growth_rate,
               spectral_gap = result$spectral_gap,
               eigenvector_nonnegative = result$eigenvector_nonnegative,
               status = result$status, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

if (identical(opt$validate, "TRUE")) {
  reference <- read.delim(file.path(opt$figure4_dir,
                                    "fixed_o2_dominant_ploidy_201grid.tsv"))
  rp <- prepare_run_params(base_params, "invivo", cfg, 0)
  for (o2 in c(0, 2.5, 5)) {
    actual <- fixo2_dominant_attractor_one("seed25", rp, model_env, cfg, o2)
    expected <- reference[reference$seed_id == "seed25" & abs(reference$O2_pct - o2) < 1e-12, ]
    if (nrow(expected) != 1L || actual$status != "ok" ||
        abs(actual$dominant_mean_ploidy - expected$dominant_mean_ploidy) > 1e-8 ||
        abs(actual$dominant_growth_rate - expected$dominant_growth_rate) > 1e-8) {
      stop("Figure 4 reference mismatch at O2=", o2)
    }
    message("Reference matched at O2=", o2, ": ploidy=", signif(actual$dominant_mean_ploidy, 8),
            ", growth=", signif(actual$dominant_growth_rate, 8))
  }
}

start <- Sys.time()
indices <- split(seq_len(nrow(samples)), rep(seq_len(workers), length.out = nrow(samples)))
part_paths <- sprintf("%s.part%02d.gz", opt$out, seq_along(indices))
worker <- function(part) {
  results <- lapply(indices[[part]], evaluate_one)
  table <- do.call(rbind, results)
  con <- gzfile(part_paths[[part]], "wt")
  on.exit(close(con))
  utils::write.table(table, con, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = part == 1L, na = "NA")
  c(rows = nrow(table), failures = sum(table$status != "ok" |
                                     is.na(table$eigenvector_nonnegative) |
                                     !table$eigenvector_nonnegative))
}
counts <- if (workers == 1L) list(worker(1L)) else
  parallel::mclapply(seq_along(indices), worker, mc.cores = workers, mc.preschedule = FALSE)
if (any(vapply(counts, inherits, logical(1L), "try-error"))) stop("Parallel worker failed")
totals <- colSums(do.call(rbind, counts))
if (totals[["rows"]] != nrow(samples) * length(oxygen) || totals[["failures"]] != 0L) {
  stop("Incomplete/nonvalid model evaluation: rows=", totals[["rows"]],
       ", failures=", totals[["failures"]])
}
combined_path <- paste0(opt$out, ".incomplete")
if (!file.copy(part_paths[[1L]], combined_path, overwrite = TRUE)) stop("Cannot create output")
if (length(part_paths) > 1L && !all(file.append(combined_path, part_paths[-1L]))) {
  stop("Cannot concatenate worker outputs")
}
if (!file.rename(combined_path, opt$out)) stop("Cannot finalize output")
unlink(part_paths)
message("Evaluated ", totals[["rows"]], " fixed-O2 operators in ",
        round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1), " seconds")
