#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(argv) {
  out <- list()
  if (!length(argv)) return(out)
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    body <- sub("^--", "", arg)
    parts <- strsplit(body, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    val <- if (length(parts) > 1L) paste(parts[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  y <- suppressWarnings(as.numeric(x))
  if (length(y) != 1L || !is.finite(y)) default else y
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  y <- suppressWarnings(as.integer(x))
  if (length(y) != 1L || is.na(y)) default else y
}

resolve_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- args_full[grepl("^--file=", args_full)]
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
}

script_dir <- normalizePath(resolve_script_dir(), mustWork = FALSE)
workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
cpp_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.cpp")
cache_dir <- file.path(tempdir(), ".rcpp_cache_o2_supply_demand_map_diag")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(file = cpp_path, rebuild = FALSE, showOutput = FALSE, verbose = FALSE, cacheDir = cache_dir)

balanced_copy_vector <- function(N, n_chr = 22L) {
  base <- N %/% n_chr
  rem <- N %% n_chr
  out <- rep(base, n_chr)
  if (rem > 0L) out[seq_len(rem)] <- out[seq_len(rem)] + 1L
  out
}

sample_dirichlet_multinomial_with_floor <- function(N, alpha, n_chr = 22L) {
  if (N < n_chr) return(rep(1L, n_chr))
  excess <- N - n_chr
  if (excess <= 0L) return(rep(1L, n_chr))
  w <- rgamma(n_chr, shape = alpha, rate = 1)
  if (!all(is.finite(w)) || sum(w) <= 0) w <- rep(1, n_chr)
  probs <- w / sum(w)
  draws <- as.integer(rmultinom(1L, size = excess, prob = probs)[, 1])
  1L + draws
}

nullisomy_risk_curve <- function(N, mode, alpha = 100, mc_samples = 10000L, seed = 12345L, N_unit = 22L) {
  m_vals <- 0:N
  surv <- vapply(
    m_vals,
    function(m) cpp_o2simps_loss_survival_nullisomy(
      as.integer(N),
      as.integer(m),
      gamma_loss = 1.0,
      N_unit = as.integer(N_unit),
      nullisomy_hidden_copy_mode = mode,
      nullisomy_dirichlet_alpha = as.numeric(alpha),
      nullisomy_dirichlet_mc_samples = as.integer(mc_samples),
      nullisomy_dirichlet_seed = as.integer(seed)
    )[[1]],
    numeric(1)
  )
  data.frame(m_loss = m_vals, risk = pmin(pmax(1 - surv, 0), 1))
}

build_risk_table <- function(N_vals, specs, mc_samples, seed, N_unit = 22L) {
  rows <- vector("list", length = 0L)
  for (N in N_vals) {
    for (i in seq_len(nrow(specs))) {
      sp <- specs[i, , drop = FALSE]
      tab <- nullisomy_risk_curve(
        N = N,
        mode = sp$mode[[1]],
        alpha = sp$alpha[[1]],
        mc_samples = mc_samples,
        seed = seed,
        N_unit = N_unit
      )
      rows[[length(rows) + 1L]] <- transform(
        tab,
        N = N,
        curve_id = sp$curve_id[[1]],
        label = sp$label[[1]],
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

build_variance_table <- function(N_vals, alpha_vals, variance_samples, variance_seed, n_chr = 22L) {
  set.seed(variance_seed)
  rows <- vector("list", length = 0L)
  for (N in N_vals) {
    balanced_var <- stats::var(balanced_copy_vector(N, n_chr))
    rows[[length(rows) + 1L]] <- data.frame(
      N = N,
      curve_id = "balanced",
      label = "balanced / maximal buffering",
      mean_within_class_var = balanced_var,
      q10 = balanced_var,
      q90 = balanced_var,
      stringsAsFactors = FALSE
    )
    for (alpha in alpha_vals) {
      samples <- replicate(
        variance_samples,
        stats::var(sample_dirichlet_multinomial_with_floor(N, alpha = alpha, n_chr = n_chr)),
        simplify = TRUE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        N = N,
        curve_id = paste0("alpha_", alpha),
        label = paste0("Dirichlet-multinomial alpha=", alpha),
        mean_within_class_var = mean(samples),
        q10 = as.numeric(stats::quantile(samples, 0.10, names = FALSE)),
        q90 = as.numeric(stats::quantile(samples, 0.90, names = FALSE)),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

main <- function(argv) {
  args <- parse_args(argv)
  out_dir <- normalizePath(args$out_dir %||% file.path(getwd(), "nullisomy_hidden_copy_diagnostics"), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  N_vals <- c(44L, 53L, 88L, 90L, 95L)
  mc_samples <- as_int(args$mc_samples, 10000L)
  variance_samples <- as_int(args$variance_samples, 5000L)
  seed <- as_int(args$seed, 12345L)
  variance_seed <- as_int(args$variance_seed, 24680L)

  specs <- data.frame(
    curve_id = c("balanced", "alpha_1000", "alpha_100", "alpha_10", "alpha_1", "alpha_0.5"),
    mode = c("balanced", rep("dirichlet_multinomial", 5L)),
    alpha = c(NA_real_, 1000, 100, 10, 1, 0.5),
    label = c(
      "balanced / maximal buffering",
      "Dirichlet-multinomial alpha=1000",
      "Dirichlet-multinomial alpha=100",
      "Dirichlet-multinomial alpha=10",
      "Dirichlet-multinomial alpha=1",
      "Dirichlet-multinomial alpha=0.5"
    ),
    stringsAsFactors = FALSE
  )

  risk_tab <- build_risk_table(N_vals = N_vals, specs = specs, mc_samples = mc_samples, seed = seed)
  variance_tab <- build_variance_table(N_vals = N_vals, alpha_vals = c(1000, 100, 10, 1, 0.5), variance_samples = variance_samples, variance_seed = variance_seed)

  risk_tab$label <- factor(risk_tab$label, levels = specs$label)
  variance_tab$label <- factor(variance_tab$label, levels = c("balanced / maximal buffering", specs$label[-1]))

  palette_vals <- c(
    "balanced / maximal buffering" = "#17324c",
    "Dirichlet-multinomial alpha=1000" = "#4c6b8a",
    "Dirichlet-multinomial alpha=100" = "#6f8f3a",
    "Dirichlet-multinomial alpha=10" = "#cf8b2a",
    "Dirichlet-multinomial alpha=1" = "#c44e52",
    "Dirichlet-multinomial alpha=0.5" = "#7a3e9d"
  )

  utils::write.table(risk_tab, file.path(out_dir, "nullisomy_risk_curves.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(variance_tab, file.path(out_dir, "copy_variance_by_alpha.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  risk_caption <- paste(
    "Balanced is shown as a separate deterministic maximal-buffering baseline.",
    "Under the one-copy-floor Dirichlet-multinomial prior, alpha -> infinity yields an equal-probability multinomial allocation of excess copies, not the deterministic balanced vector."
  )
  g_risk <- ggplot(risk_tab, aes(m_loss, risk, color = label)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ N, scales = "free_x") +
    scale_color_manual(values = palette_vals) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Nullisomy Risk Curves by Hidden-Copy Approximation",
      subtitle = "Balanced is a maximal-buffering baseline; high-alpha Dirichlet approaches an equal-probability multinomial-like limit, not balanced.",
      x = "Chromosome copies lost (m)",
      y = "Nullisomy risk",
      color = NULL,
      caption = risk_caption
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  var_caption <- paste(
    "Within-sample copy-number variance across chromosome classes.",
    "Balanced has the minimum deterministic variance for a given total N; Dirichlet-multinomial samples remain heterogeneous even at large alpha."
  )
  g_var <- ggplot(variance_tab, aes(label, mean_within_class_var, color = label, group = label)) +
    geom_linerange(aes(ymin = q10, ymax = q90), linewidth = 0.7, alpha = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ N, scales = "free_y") +
    scale_color_manual(values = palette_vals) +
    labs(
      title = "Hidden Copy-Vector Variance by Alpha",
      subtitle = "High-alpha Dirichlet remains more variable than the balanced maximal-buffering baseline.",
      x = NULL,
      y = "Variance across chromosome classes",
      color = NULL,
      caption = var_caption
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(out_dir, "nullisomy_risk_curves_hidden_copy.pdf"), g_risk, width = 11.5, height = 8.5)
  ggsave(file.path(out_dir, "nullisomy_risk_curves_hidden_copy.png"), g_risk, width = 11.5, height = 8.5, dpi = 160)
  ggsave(file.path(out_dir, "copy_variance_hidden_copy.pdf"), g_var, width = 11.5, height = 8.5)
  ggsave(file.path(out_dir, "copy_variance_hidden_copy.png"), g_var, width = 11.5, height = 8.5, dpi = 160)

  note_lines <- c(
    "# Nullisomy Hidden-Copy Diagnostics",
    "",
    "Diagnostic framing:",
    "- `balanced / maximal buffering` is a separate deterministic baseline, not the `alpha -> infinity` limit of the current one-copy-floor Dirichlet-multinomial model.",
    "- Under the one-copy-floor prior, if `E = N - M` and `(e_1,...,e_M) ~ DirichletMultinomial(E; alpha,...,alpha)`, then as `alpha -> infinity`, `(e_1,...,e_M)` approaches `Multinomial(E; 1/M,...,1/M)`.",
    "- That limit remains stochastic and therefore does not collapse to the deterministic balanced copy vector.",
    "",
    "Outputs:",
    paste0("- Risk curves: `", file.path(out_dir, "nullisomy_risk_curves_hidden_copy.pdf"), "` and `.png`"),
    paste0("- Copy-variance diagnostic: `", file.path(out_dir, "copy_variance_hidden_copy.pdf"), "` and `.png`"),
    paste0("- Tables: `", file.path(out_dir, "nullisomy_risk_curves.tsv"), "` and `", file.path(out_dir, "copy_variance_by_alpha.tsv"), "`"),
    "",
    paste("Monte Carlo samples for risk curves:", mc_samples),
    paste("Monte Carlo samples for copy-variance diagnostic:", variance_samples)
  )
  writeLines(note_lines, con = file.path(out_dir, "diagnostic_note.md"))
  message("Wrote nullisomy hidden-copy diagnostics to: ", out_dir)
}

main(commandArgs(trailingOnly = TRUE))
