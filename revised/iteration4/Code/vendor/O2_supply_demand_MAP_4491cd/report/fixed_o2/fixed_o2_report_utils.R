# Fixed-O2 report assembly helpers.
# These functions summarize already materialized tables and invoke the standalone HTML assembler.

fo2_write_report <- function(out_dir, run_dir, label_file, attractors, regime_summary, tests, correlations) {
  counts <- if (nrow(attractors)) {
    rows <- lapply(split(attractors, attractors$trajectory_regime, drop = TRUE), function(d) {
      data.frame(
        trajectory_regime = d$trajectory_regime[[1]],
        n_seed_o2 = nrow(d),
        n_seed = length(unique(d$seed_id)),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  } else {
    data.frame()
  }
  low_o2 <- regime_summary[regime_summary$O2_pct %in% c(0, 0.01, 0.05, 0.1), , drop = FALSE]
  top_tests <- tests[is.finite(tests$BH_adjusted_p_value), , drop = FALSE]
  top_tests <- top_tests[order(top_tests$BH_adjusted_p_value, top_tests$wilcox_p_value), , drop = FALSE]
  top_corr <- correlations[is.finite(correlations$abs_spearman_rho), , drop = FALSE]
  top_corr <- top_corr[order(-top_corr$abs_spearman_rho), , drop = FALSE]
  gap_summary_path <- file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_regime_summary.tsv")
  gap_summary <- if (file.exists(gap_summary_path)) fo2_read_tsv(gap_summary_path) else data.frame()
  key_gap <- gap_summary[
    nrow(gap_summary) > 0 &
      gap_summary$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2") &
      gap_summary$O2_pct %in% c(0, 0.5, 1, 2, 5),
    ,
    drop = FALSE
  ]
  key_gap <- key_gap[, intersect(c(
    "O2_pct", "trajectory_regime", "n_seed", "spectral_gap_median",
    "time_to_10x_days_median", "time_to_100x_days_median",
    "log10_advantage_1000d_median", "fraction_gap_ge_0p005",
    "fraction_gap_ge_0p01"
  ), names(key_gap)), drop = FALSE]
  lines <- c(
    "# In vivo fixed-O2 ploidy attractor analysis",
    "",
    paste0("- run_dir: ", run_dir),
    paste0("- optional_label_file: ", if (nzchar(label_file)) label_file else "<none>"),
    "- FixO2 mode rule: dominant_mean_ploidy >= 2 is mode1; dominant_mean_ploidy < 2 is mode2.",
    paste0("- analyzed seeds: ", length(unique(attractors$seed_id))),
    "",
    "## Regime counts",
    "",
    paste(utils::capture.output(print(counts, row.names = FALSE)), collapse = "\n"),
    "",
    "## Low-O2 Attractor Summary",
    "",
    paste(utils::capture.output(print(low_o2, row.names = FALSE)), collapse = "\n"),
    "",
    "## Strongest Mode1-vs-Mode2 Attractor Differences",
    "",
    paste(utils::capture.output(print(utils::head(top_tests, 20), row.names = FALSE)), collapse = "\n"),
    "",
    "## Strongest Parameter-Attractor Correlations",
    "",
    paste(utils::capture.output(print(utils::head(top_corr, 30), row.names = FALSE)), collapse = "\n"),
    "",
    "## Seed-Level Mode1/Mode2 Stack Outputs",
    "",
    "- Table: `tables/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv`",
    "- Figure: `figures/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.pdf`",
    "",
    "## Spectral Gap Reliability Outputs",
    "",
    "- Seed table: `tables/fixed_o2_attractor_spectral_gap_by_seed.tsv`",
    "- Regime summary: `tables/fixed_o2_attractor_spectral_gap_regime_summary.tsv`",
    "- Figures: `figures/fixed_o2_spectral_gap_all_seed_stack_mode1_mode2.pdf`, `figures/fixed_o2_time_to_10x_all_seed_stack_mode1_mode2.pdf`, `figures/fixed_o2_gap_reliability_fraction_mode1_mode2.pdf`",
    "- Composite table: `tables/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv`",
    "- Composite figure: `figures/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.pdf`",
    "- Violin table: `tables/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv`",
    "- Violin figure: `figures/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.pdf`",
    "",
    paste(utils::capture.output(print(key_gap, row.names = FALSE)), collapse = "\n"),
    "",
    "## Interpretation Boundary",
    "",
    "These are fixed-O2 attractors from the fitted model parameters. They test model-internal mechanism consistency under standardized O2; they do not by themselves prove biological causality."
  )
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}


cf2_report <- function(out_dir, args, regime_summary, tests, correlations) {
  key_summary <- regime_summary[regime_summary$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2") &
                                  regime_summary$O2_pct %in% c(1, 2, 5), , drop = FALSE]
  key_tests <- tests[tests$metric %in% c("terminal_mean_ploidy", "time_crossing_ploidy_1p5_down_censored") &
                       tests$O2_pct %in% c(1, 2, 5), , drop = FALSE]
  key_tests <- key_tests[order(key_tests$BH_adjusted_p_value), , drop = FALSE]
  top_corr <- data.frame()
  if (is.data.frame(correlations) && nrow(correlations) &&
      all(c("correlation_scope", "O2_pct", "abs_spearman_rho", "BH_adjusted_p_value") %in% names(correlations))) {
    top_corr <- correlations[correlations$correlation_scope == "mode1_mode2_only" & correlations$O2_pct %in% c(1, 2, 5), , drop = FALSE]
    if (nrow(top_corr)) {
      top_corr <- head(top_corr[order(-top_corr$abs_spearman_rho, top_corr$BH_adjusted_p_value), , drop = FALSE], 30L)
    }
  }
  lines <- c(
    "# In vivo fixed-O2 trajectory counterfactual analysis",
    "",
    paste0("- run_dir: `", o2ipa_as_chr(args$run_dir, ""), "`"),
    paste0("- optional_label_file: `", o2ipa_as_chr(args$label_file, ""), "`"),
    "- FixO2 mode rule: dominant_mean_ploidy >= 2 is mode1; dominant_mean_ploidy < 2 is mode2.",
    paste0("- O2 grid: `", o2ipa_as_chr(args$o2_grid, "0,0.5,1,2,5"), "`"),
    paste0("- time grid: default dense early grid unless `--time_grid` was supplied"),
    "",
    "## Key high-O2 terminal summaries",
    "",
    paste(capture.output(print(key_summary, row.names = FALSE)), collapse = "\n"),
    "",
    "## Key mode1-vs-mode2 tests",
    "",
    paste(capture.output(print(head(key_tests, 30L), row.names = FALSE)), collapse = "\n"),
    "",
    "## Top mode1/mode2 parameter correlations at O2=1/2/5",
    "",
    paste(capture.output(print(top_corr, row.names = FALSE)), collapse = "\n"),
    "",
    "## Notes",
    "",
    "- This counterfactual uses the same fixed-O2 matrix as the attractor analysis.",
    "- Initial conditions are point masses at N=44 (`init_2N`) and N=88 (`init_4N`) unless clipped by the configured N grid.",
    "- State vectors are propagated through the exact eigen-decomposition of the fixed linear system and normalized to composition at each time point.",
    "- `time_crossing_ploidy_1p5_down_censored` is set to max_time + 1 when the trajectory never crosses 1.5."
  )
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}


fixo2_render_html_report <- function(args, repo_root, out_dir) {
  render_html <- o2ipa_as_bool(fixo2_arg_value(args, "render_html_report", "generate_html_report", TRUE), TRUE)
  if (!isTRUE(render_html)) {
    message("Skipping FixO2 HTML report rendering because --render_html_report=FALSE.")
    return(invisible(NULL))
  }

  render_script <- fixo2_resolve_repo_path(
    o2ipa_as_chr(
      fixo2_arg_value(args, "html_report_script", "report_script", fixo2_default_html_report_script(repo_root)),
      fixo2_default_html_report_script(repo_root)
    ),
    repo_root,
    mustWork = FALSE
  )
  if (!file.exists(render_script)) {
    stop("FixO2 HTML report script does not exist: ", render_script, call. = FALSE)
  }

  html_out_dir <- fixo2_resolve_repo_path(
    o2ipa_as_chr(fixo2_arg_value(args, "html_report_dir", "report_out_dir", file.path(out_dir, "html_report")), file.path(out_dir, "html_report")),
    repo_root,
    mustWork = FALSE
  )
  report_basename <- o2ipa_as_chr(fixo2_arg_value(args, "html_report_basename", "report_basename", "index"), "index")
  rscript <- file.path(R.home("bin"), "Rscript")
  if (!file.exists(rscript)) rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) stop("Rscript executable was not found for FixO2 HTML report rendering.", call. = FALSE)

  cmd_args <- c(
    render_script,
    paste0("--analysis_dir=", out_dir),
    paste0("--out_dir=", html_out_dir),
    paste0("--report_basename=", report_basename)
  )
  message("Rendering FixO2 HTML report: ", file.path(html_out_dir, paste0(report_basename, ".html")))
  status <- system2(rscript, cmd_args)
  if (!identical(as.integer(status), 0L)) {
    stop("FixO2 HTML report rendering failed with exit status ", status, call. = FALSE)
  }
  out_path <- file.path(html_out_dir, paste0(report_basename, ".html"))
  if (!file.exists(out_path) || is.na(file.info(out_path)$size) || file.info(out_path)$size <= 0) {
    stop("FixO2 HTML report renderer completed but did not create a non-empty report: ", out_path, call. = FALSE)
  }
  message("Completed FixO2 HTML report: ", normalizePath(out_path, mustWork = FALSE))
  invisible(out_path)
}
