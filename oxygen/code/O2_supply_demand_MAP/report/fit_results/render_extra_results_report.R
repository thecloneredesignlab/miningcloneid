#!/usr/bin/env Rscript

.o2sd_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2sd_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_html_utils.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_report_utils.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "report", "run_provenance_report.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

escape_html <- o2sd_html_escape_full

make_figure_spec <- function(extra_results_dir, filename, title, legend) {
  path <- file.path(extra_results_dir, filename)
  if (!file.exists(path)) {
    stop("Missing required figure for extra_results report: ", path)
  }
  list(
    filename = filename,
    path = normalizePath(path, mustWork = TRUE),
    title = title,
    legend = legend
  )
}

make_figure_spec_optional <- function(extra_results_dir, filename, title, legend) {
  path <- file.path(extra_results_dir, filename)
  if (!file.exists(path)) {
    return(NULL)
  }
  list(
    filename = filename,
    path = normalizePath(path, mustWork = TRUE),
    title = title,
    legend = legend
  )
}

with_figure_label <- function(fig, label) {
  if (is.null(fig)) return(NULL)
  fig$label <- as.character(label)
  fig
}

figure_display_label <- function(fig, fallback_index) {
  label <- fig$display_label %||% fig$label %||% fallback_index
  as.character(label)
}

build_prediction_figure_specs <- function(extra_results_dir) {
  figs <- list(
    make_figure_spec_optional(
      extra_results_dir,
      "predict1000_ploidy_mean_ci_2N_4N.pdf",
      "Cross-seed 1000-day Ploidy Prediction: 2N and 4N",
      "Mean trajectories and 95% confidence intervals across seed-level scenario means. Dashed lines show the cross-seed min/max envelope; the right y-axis reports ploidy."
    )
  )
  for (cohort in c("2N", "4N")) {
    figs <- c(figs, list(
      make_figure_spec_optional(
        extra_results_dir,
        paste0("predict1000_burden_total_log_seed_mean_", cohort, ".pdf"),
        paste0("Cross-seed 1000-day Total Burden Prediction: ", cohort),
        "All seed-level scenario-mean total-burden trajectories are drawn as thin grey lines on the log10 scale. The colored solid line is the cross-seed mean."
      )
    ))
  }
  for (cohort in c("2N", "4N")) {
    figs <- c(figs, list(
      make_figure_spec_optional(
        extra_results_dir,
        paste0("predict1000_ploidy_seed_curves_", cohort, ".pdf"),
        paste0("Cross-seed 1000-day Ploidy Seed Trajectories: ", cohort),
        "All seed-level scenario-mean trajectories are drawn as thin grey lines. The colored solid line is the cross-seed mean and the colored dashed line is the cross-seed median."
      )
    ))
  }
  Filter(Negate(is.null), figs)
}

read_report_run_mode <- function(extra_results_dir) {
  seed_summary_path <- file.path(extra_results_dir, "seed_summary.tsv")
  if (!file.exists(seed_summary_path)) {
    return("unknown")
  }
  seed_summary <- utils::read.delim(seed_summary_path, check.names = FALSE, stringsAsFactors = FALSE, nrows = 1000)
  fit_mode <- if ("fit_mode" %in% names(seed_summary)) unique(as.character(seed_summary$fit_mode)) else character(0)
  if (any(fit_mode == "fit_joint", na.rm = TRUE) ||
      ("objective_invivo" %in% names(seed_summary) && any(is.finite(suppressWarnings(as.numeric(seed_summary$objective_invivo)))))) {
    return("fit_joint")
  }
  if (any(fit_mode == "fit_invitro", na.rm = TRUE) ||
      ("objective_total" %in% names(seed_summary) && any(is.finite(suppressWarnings(as.numeric(seed_summary$objective_total))))) ||
      ("growth_loglik" %in% names(seed_summary) && any(is.finite(suppressWarnings(as.numeric(seed_summary$growth_loglik)))))) {
    return("fit_invitro")
  }
  "fit_invivo"
}

pred1000_top_seed_legend <- function(extra_results_dir, n = 3L) {
  seed_summary_path <- file.path(extra_results_dir, "seed_summary.tsv")
  base_legend <- "Boundary forest shows all seeds; the highlighted seeds are the top eligible seeds. Eligibility requires all final 2N and 4N 1000-day predictions to be above 44; eligible seeds are ordered by raw objective."
  if (!file.exists(seed_summary_path)) {
    return(base_legend)
  }
  seed_summary <- tryCatch(
    utils::read.delim(seed_summary_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(seed_summary) || !nrow(seed_summary) || !"seed" %in% names(seed_summary)) {
    return(base_legend)
  }
  if (!"pred1000_both_gt44" %in% names(seed_summary)) {
    return(base_legend)
  }
  eligible <- tolower(trimws(as.character(seed_summary$pred1000_both_gt44))) %in% c("true", "t", "1", "yes", "y", "on")
  seed_summary <- seed_summary[eligible, , drop = FALSE]
  if (!nrow(seed_summary)) {
    return(paste(base_legend, "No eligible seeds were found."))
  }
  rank_col <- if ("objective" %in% names(seed_summary)) {
    "objective"
  } else {
    NA_character_
  }
  if (!is.na(rank_col)) {
    seed_summary[[rank_col]] <- suppressWarnings(as.numeric(seed_summary[[rank_col]]))
    seed_summary <- seed_summary[order(seed_summary[[rank_col]], seed_summary$seed, na.last = TRUE), , drop = FALSE]
  }
  top_seeds <- utils::head(as.character(seed_summary$seed), as.integer(n))
  top_text <- paste(paste0("Top ", seq_along(top_seeds), ": ", top_seeds), collapse = "; ")
  paste(base_legend, top_text)
}

build_invivo_figure_specs <- function(extra_results_dir) {
  c(Filter(Negate(is.null), list(
    make_figure_spec_optional(
      extra_results_dir,
      "objective_vs_boundary_risk.pdf",
      "Objective vs Boundary Risk",
      "Objective score against minimum relative distance to the nearest fitted parameter bound."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "objective_components_violin.pdf",
      "Objective Components Violin",
      "Across-seed distributions of total objective, burden objective, and ploidy objective."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invivo_top20_seed_parameter_umap.pdf",
      "In Vivo Top 20 Seed Parameter UMAP",
      "UMAP embedding of the top 20 in vivo seeds by objective, using transformed values for the full parameter vector. Source coordinates are saved as invivo_top20_seed_parameter_umap_coords.tsv."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "joint_objective_components.pdf",
      "Joint Objective Components",
      "Across-seed violin and boxplot distributions of joint total objective, in vivo objective, and in vitro objective."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "joint_objective_tradeoff.pdf",
      "Joint In Vivo vs In Vitro Objective Tradeoff",
      "Seed-level in vivo objective against in vitro objective; color shows total joint objective."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "joint_soft_coupling_delta_magnitude.pdf",
      "Joint Soft-Coupling Delta Magnitudes",
      "Absolute context splits on the transformed optimizer scale for soft-coupled parameters."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "joint_soft_coupling_penalty_ranking.pdf",
      "Joint Soft-Coupling Penalty Ranking",
      "Soft-coupling penalty contributions by parameter and seed."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "joint_soft_coupling_vivo_vitro_pairs.pdf",
      "Joint Soft-Coupled In Vivo vs In Vitro Parameters",
      "Natural-scale paired parameter values derived from center plus or minus delta over two."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_objective_components.pdf",
      "In Vitro Objective Components",
      "Seed-level total in vitro objective and growth, ploidy, and flow negative log-likelihood contributions."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "parameter_boundary_forest.pdf",
      "Parameter Boundary Forest",
      "Relative fitted positions of active parameters within their transformed bounds across all seeds. Highlighted seeds are the top 3 by raw objective."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "parameter_boundary_forest_pred1000_gt44_top3.pdf",
      "Parameter Boundary Forest (Pred1000 > 44 Top 3)",
      pred1000_top_seed_legend(extra_results_dir)
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "parameter_boundary_forest_log_x.pdf",
      "Parameter Boundary Forest (Log X Scale)",
      "Original fitted parameter values across all seeds on a log10 x-axis. Highlighted seeds are the top 3 by raw objective. Horizontal lines span original lower-to-upper parameter bounds; non-positive raw values or bounds are floored and labeled 0."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "parameter_boundary_forest_pred1000_gt44_top3_log_x.pdf",
      "Parameter Boundary Forest (Pred1000 > 44 Top 3, Log X Scale)",
      paste(
        pred1000_top_seed_legend(extra_results_dir),
        "The x-axis shows original fitted parameter values on a log10 scale; non-positive raw values or bounds are floored and labeled 0."
      )
    )
  )), build_prediction_figure_specs(extra_results_dir))
}

build_invitro_figure_specs <- function(extra_results_dir) {
  Filter(Negate(is.null), list(
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_top20_seed_parameter_umap.pdf",
      "In Vitro Top 20 Seed Parameter UMAP",
      "UMAP embedding of the top 20 in vitro seeds by objective, using transformed values for the full parameter vector. Source coordinates are saved as invitro_top20_seed_parameter_umap_coords.tsv."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "optimization_diagnostics.pdf",
      "Optimization Diagnostics",
      "All-seed in vitro optimization diagnostics based on the current run outputs. Source data are saved in this extra_results directory as seed_summary.tsv and parameter_boundary_long.tsv."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_parameter_positions_log_x.pdf",
      "Fitted Parameter Positions (Log X Scale)",
      "All-seed in vitro fitted-parameter positions shown as original parameter values on a log10 x-axis. This is the parameter panel from Optimization Diagnostics only; non-positive raw values or bounds are floored and labeled 0."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "best_fit_ploidy_likelihood_comparison.pdf",
      "Best-Fit Ploidy Likelihood Comparison",
      "Passage-level ploidy likelihood comparison for the best available seed-level parameter sets. Source data are saved in best_fit_likelihood_comparison.tsv and best_fit_ploidy_likelihood_comparison_long.tsv."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "best_fit_flow_likelihood_comparison.pdf",
      "Best-Fit Flow-Density Likelihood Comparison",
      "Passage-level processed flow-density likelihood comparison for the best available seed-level parameter sets. Source data are saved in best_fit_likelihood_comparison.tsv and best_fit_flow_likelihood_comparison_long.tsv."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_objective_components.pdf",
      "In Vitro Objective Components",
      "Seed-level total in vitro objective and growth, ploidy, and flow negative log-likelihood contributions. Lower total objective is better."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "objective_components_violin.pdf",
      "In Vitro Objective Component Distributions",
      "Across-seed distributions of the total in vitro objective and the in vitro growth, ploidy, and flow negative log-likelihood components."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "objective_vs_boundary_risk.pdf",
      "In Vitro Objective vs Boundary Risk",
      "Total in vitro objective against minimum relative distance to the nearest fitted parameter bound."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_karyotype_quantiles_multiseed.pdf",
      "Multi-Seed Predicted Chromosome-Count Quantiles vs Observed Cells",
      "Across-seed version of the predicted chromosome-count quantile timecourse shown in simulation_basics.html. Panels are split by 2N/4N and control/deprived; green curves show predicted quantiles by passage and orange points show observed single cells."
    )
  ))
}

build_joint_summary_figure_specs <- function(extra_results_dir) {
  Filter(Negate(is.null), list(
    make_figure_spec_optional(
      extra_results_dir,
      "objective_vs_boundary_risk.pdf",
      "Objective vs Boundary Risk",
      "Objective score against minimum relative distance to the nearest fitted parameter bound."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "objective_components_violin.pdf",
      "Objective Components Violin",
      "Across-seed distributions of joint total objective, in vivo objective, and in vitro objective."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "joint_objective_components.pdf",
      "Joint Objective Components",
      "Across-seed violin and boxplot distributions of joint total objective, in vivo objective, and in vitro objective."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "joint_objective_tradeoff.pdf",
      "Joint In Vivo vs In Vitro Objective Tradeoff",
      "Seed-level in vivo objective against in vitro objective; color shows total joint objective."
    ),
      make_figure_spec_optional(
        extra_results_dir,
        "parameter_boundary_forest.pdf",
        "Parameter Boundary Forest",
        "Relative fitted positions of active parameters within their transformed bounds across all seeds. Highlighted seeds are the top 3 by raw objective."
      ),
      with_figure_label(
        make_figure_spec_optional(
          extra_results_dir,
          "joint_soft_coupling_context_boundary_forest.pdf",
          "Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter",
          "Soft-coupled parameter positions for in vivo and in vitro context-specific values, computed as center plus or minus delta over two within the joint union transformed bounds. Colors distinguish context; lines connect paired seed-parameter values; highlighted seeds match the top 3 by raw objective."
        ),
        "5b"
      ),
      make_figure_spec_optional(
        extra_results_dir,
        "parameter_boundary_forest_pred1000_gt44_top3.pdf",
        "Parameter Boundary Forest (Pred1000 > 44 Top 3)",
        pred1000_top_seed_legend(extra_results_dir)
      ),
      with_figure_label(
        make_figure_spec_optional(
          extra_results_dir,
          "joint_soft_coupling_context_boundary_forest_pred1000_gt44_top3.pdf",
          "Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Pred1000 > 44 Top 3)",
          paste(
            pred1000_top_seed_legend(extra_results_dir),
            "Soft-coupled in vivo and in vitro context-specific values are computed as center plus or minus delta over two and plotted within the joint union transformed bounds. Lines connect paired seed-parameter values."
          )
        ),
        "6b"
      ),
      make_figure_spec_optional(
        extra_results_dir,
        "parameter_boundary_forest_log_x.pdf",
        "Parameter Boundary Forest (Log X Scale)",
        "Original fitted parameter values across all seeds on a log10 x-axis. Highlighted seeds are the top 3 by raw objective. Horizontal lines span original lower-to-upper parameter bounds; non-positive raw values or bounds are floored and labeled 0."
      ),
      with_figure_label(
        make_figure_spec_optional(
          extra_results_dir,
          "joint_soft_coupling_context_boundary_forest_log_x.pdf",
          "Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Log X Scale)",
          "Natural-scale in vivo and in vitro context-specific soft-coupled parameter values across all seeds on a log10 x-axis. Values are computed as center plus or minus delta over two; paired seed-parameter values are connected; horizontal lines span natural joint union bounds."
        ),
        "7b"
      ),
      make_figure_spec_optional(
        extra_results_dir,
        "parameter_boundary_forest_pred1000_gt44_top3_log_x.pdf",
        "Parameter Boundary Forest (Pred1000 > 44 Top 3, Log X Scale)",
        paste(
          pred1000_top_seed_legend(extra_results_dir),
          "The x-axis shows original fitted parameter values on a log10 scale; non-positive raw values or bounds are floored and labeled 0."
        )
      ),
      with_figure_label(
        make_figure_spec_optional(
          extra_results_dir,
          "joint_soft_coupling_context_boundary_forest_pred1000_gt44_top3_log_x.pdf",
          "Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Pred1000 > 44 Top 3, Log X Scale)",
          paste(
            pred1000_top_seed_legend(extra_results_dir),
            "The x-axis shows natural-scale soft-coupled in vivo and in vitro context-specific values on a log10 scale; paired seed-parameter values are connected; non-positive raw values or bounds are floored and labeled 0."
          )
        ),
        "8b"
      )
    ))
  }

build_joint_invivo_figure_specs <- function(extra_results_dir) {
  build_prediction_figure_specs(extra_results_dir)
}

build_joint_invitro_figure_specs <- function(extra_results_dir) {
  Filter(Negate(is.null), list(
    make_figure_spec_optional(
      extra_results_dir,
      "optimization_diagnostics.pdf",
      "Optimization Diagnostics",
      "All-seed in vitro optimization diagnostics using the joint-fit in vitro objective for seed ordering."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_parameter_positions_log_x.pdf",
      "Fitted Parameter Positions (Log X Scale)",
      "All-seed joint-fit in vitro fitted-parameter positions shown as original parameter values on a log10 x-axis. This is the parameter panel from Optimization Diagnostics only; non-positive raw values or bounds are floored and labeled 0."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "best_fit_ploidy_likelihood_comparison.pdf",
      "Best-Fit Ploidy Likelihood Comparison",
      "Passage-level ploidy likelihood comparison for the best available seed-level parameter sets."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "best_fit_flow_likelihood_comparison.pdf",
      "Best-Fit Flow-Density Likelihood Comparison",
      "Passage-level processed flow-density likelihood comparison for the best available seed-level parameter sets."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_objective_components.pdf",
      "In Vitro Objective Components",
      "Seed-level total in vitro objective and growth, ploidy, and flow negative log-likelihood contributions. Lower total objective is better."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_objective_components_violin.pdf",
      "In Vitro Objective Component Distributions",
      "Across-seed distributions of the total in vitro objective and the in vitro growth, ploidy, and flow negative log-likelihood components."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_objective_vs_boundary_risk.pdf",
      "In Vitro Objective vs Boundary Risk",
      "Total in vitro objective against minimum relative distance to the nearest fitted parameter bound."
    ),
    make_figure_spec_optional(
      extra_results_dir,
      "invitro_karyotype_quantiles_multiseed.pdf",
      "Multi-Seed Predicted Chromosome-Count Quantiles vs Observed Cells",
      "Across-seed version of the predicted chromosome-count quantile timecourse shown in simulation_basics.html. Panels are split by 2N/4N and control/deprived; green curves show predicted quantiles by passage and orange points show observed single cells."
    )
  ))
}

build_joint_figure_chapters <- function(extra_results_dir) {
  chapters <- list(
    list(title = "Joint Summary", figures = build_joint_summary_figure_specs(extra_results_dir)),
    list(title = "In Vivo", figures = build_joint_invivo_figure_specs(extra_results_dir)),
    list(title = "In Vitro", figures = build_joint_invitro_figure_specs(extra_results_dir))
  )
  chapters <- Filter(function(chapter) length(chapter$figures) > 0L, chapters)
  if (!length(chapters)) {
    stop("No supported joint summary figures were found in extra_results directory: ", extra_results_dir)
  }
  figure_idx <- 0L
  figure_label_idx <- 0L
  for (chapter_i in seq_along(chapters)) {
    if (!length(chapters[[chapter_i]]$figures)) next
    for (figure_i in seq_along(chapters[[chapter_i]]$figures)) {
      figure_idx <- figure_idx + 1L
      chapters[[chapter_i]]$figures[[figure_i]]$index <- figure_idx
      if (is.null(chapters[[chapter_i]]$figures[[figure_i]]$label)) {
        figure_label_idx <- figure_label_idx + 1L
        chapters[[chapter_i]]$figures[[figure_i]]$display_label <- as.character(figure_label_idx)
      } else {
        chapters[[chapter_i]]$figures[[figure_i]]$display_label <- chapters[[chapter_i]]$figures[[figure_i]]$label
      }
    }
  }
  chapters
}

build_figure_specs <- function(extra_results_dir) {
  run_mode <- read_report_run_mode(extra_results_dir)
  figs <- if (identical(run_mode, "fit_invitro")) {
    build_invitro_figure_specs(extra_results_dir)
  } else {
    build_invivo_figure_specs(extra_results_dir)
  }
  if (!length(figs)) {
    stop("No supported figures were found in extra_results directory: ", extra_results_dir)
  }
  figs
}

infer_run_label <- function(extra_results_dir) {
  base <- basename(extra_results_dir)
  if (identical(base, "extra_results")) {
    return(basename(dirname(extra_results_dir)))
  }
  base
}

report_magick_available <- o2sd_report_magick_available
report_gs_available <- o2sd_report_gs_available
report_base64enc_available <- o2sd_report_base64enc_available
file_to_data_uri <- o2sd_report_file_to_data_uri
render_pdf_preview_png_gs <- o2sd_report_render_pdf_preview_png_gs

pdf_to_data_uri <- function(pdf_path, density = 180) {
  if (!report_magick_available() && !report_gs_available()) {
    return(file_to_data_uri(pdf_path, mime = "application/pdf"))
  }
  png_path <- tempfile("o2_extra_results_", fileext = ".png")
  on.exit(unlink(png_path, force = TRUE), add = TRUE)
  if (report_magick_available()) {
    img <- magick::image_read(pdf_path, density = density)
    if (length(img) > 1L) {
      img <- img[1]
    }
    magick::image_write(img, path = png_path, format = "png")
  } else {
    render_pdf_preview_png_gs(pdf_path, png_path, density = density)
  }
  file_to_data_uri(png_path, mime = "image/png")
}

figure_media_html <- function(fig) {
  data_uri <- pdf_to_data_uri(fig$path)
  if (report_magick_available() || report_gs_available()) {
    sprintf(
      '<div class="report-figure"><img src="%s" alt="%s" class="report-figure-image"/></div>',
      data_uri,
      escape_html(fig$title)
    )
  } else {
    sprintf(
      paste0(
        '<div class="report-figure">',
        '<object data="%s" type="application/pdf" class="report-figure-object">',
        '<div class="report-figure-fallback"><a href="%s">Open PDF figure</a></div>',
        '</object></div>'
      ),
      data_uri,
      escape_html(fig$filename)
    )
  }
}

figure_pair_match <- function(left_title, right_title, left_pattern, right_pattern) {
  grepl(left_pattern, left_title) && grepl(right_pattern, right_title)
}

figure_pair_match_any_order <- function(left_title, right_title, pattern_a, pattern_b) {
  figure_pair_match(left_title, right_title, pattern_a, pattern_b) ||
    figure_pair_match(left_title, right_title, pattern_b, pattern_a)
}

figure_half_width_placeholder <- function(fig) {
  title <- fig$title %||% ""
  grepl("Top 20 Seed Parameter UMAP", title, fixed = TRUE)
}

figure_layout_groups <- function(figure_specs) {
  n <- length(figure_specs)
  groups <- list()
  i <- 1L
  while (i <= n) {
    if (figure_half_width_placeholder(figure_specs[[i]])) {
      groups <- c(groups, list(c(i, NA_integer_)))
      i <- i + 1L
      next
    }
    if (i < n) {
      left_title <- figure_specs[[i]]$title %||% ""
      right_title <- figure_specs[[i + 1L]]$title %||% ""
      if (
        figure_pair_match_any_order(left_title, right_title, "Objective vs Boundary Risk", "Objective Components Violin") ||
          figure_pair_match_any_order(left_title, right_title, "Joint Objective Components", "Joint In Vivo vs In Vitro Objective Tradeoff") ||
          figure_pair_match_any_order(left_title, right_title, "Best-Fit Ploidy Likelihood Comparison", "Best-Fit Flow-Density Likelihood Comparison") ||
          figure_pair_match_any_order(left_title, right_title, "In Vitro Objective Component Distributions", "In Vitro Objective vs Boundary Risk")
      ) {
        groups <- c(groups, list(seq.int(i, i + 1L)))
        i <- i + 2L
        next
      }
    }
    groups <- c(groups, list(i))
    i <- i + 1L
  }
  groups
}

blank_figure_block <- function() {
  paste0(
    '<section class="report-section report-section--blank">',
    '<div class="report-figure report-figure-blank"></div>',
    '</section>'
  )
}

read_table_optional <- o2sd_report_read_table_optional
report_truthy <- o2sd_report_truthy

finite_or_na <- function(x, fn) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  fn(x)
}

first_nonempty <- function(x, default = "") {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) default else x[[1L]]
}

format_report_number <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (is.na(x)) return("")
  if (!is.finite(x)) return(as.character(x))
  if (abs(x - round(x)) <= 1e-12 * max(1, abs(x))) {
    return(format(round(x), scientific = FALSE, trim = TRUE))
  }
  rounded <- round(x, digits = 3)
  if (rounded == 0) {
    return(format(signif(x, 3), scientific = TRUE, trim = TRUE))
  }
  formatC(rounded, format = "f", digits = 3)
}

format_report_value <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  if (is.numeric(x) || is.integer(x)) {
    return(vapply(as.numeric(x), format_report_number, character(1), USE.NAMES = FALSE))
  }
  out <- as.character(x)
  x_trim <- trimws(out)
  numeric_like <- grepl("^[-+]?((\\d+\\.?\\d*)|(\\.\\d+))([eE][-+]?\\d+)?$", x_trim)
  numeric_like[is.na(numeric_like)] <- FALSE
  numeric_like <- numeric_like & nzchar(x_trim)
  if (any(numeric_like)) {
    num <- suppressWarnings(as.numeric(x_trim[numeric_like]))
    out[numeric_like] <- vapply(num, format_report_number, character(1), USE.NAMES = FALSE)
  }
  out[is.na(out)] <- ""
  out
}

table_to_html <- function(df, max_rows = 120) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return("<p class=\"report-empty\">No rows available.</p>")
  }
  if (nrow(df) > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
  }
  df[] <- lapply(df, format_report_value)
  header <- paste0("<th>", escape_html(names(df)), "</th>", collapse = "")
  body <- apply(
    df,
    1L,
    function(row) paste0("<tr>", paste0("<td>", escape_html(row), "</td>", collapse = ""), "</tr>")
  )
  paste0(
    '<table class="report-table"><thead><tr>',
    header,
    "</tr></thead><tbody>",
    paste0(body, collapse = "\n"),
    "</tbody></table>"
  )
}

build_parameter_summary_table <- function(extra_results_dir) {
  path <- file.path(extra_results_dir, "parameter_boundary_long.tsv")
  tab <- read_table_optional(path, sep = "\t")
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab) || !"param_prototype" %in% names(tab)) {
    return(NULL)
  }
  tab <- o2sd_add_parameter_descriptions(tab)
  active <- if ("active_in_fit" %in% names(tab)) report_truthy(tab$active_in_fit) else rep(FALSE, nrow(tab))
  estimate <- if ("estimate" %in% names(tab)) report_truthy(tab$estimate) else rep(FALSE, nrow(tab))
  current_param <- active | estimate
  tab <- tab[current_param, , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  params <- unique(trimws(as.character(tab$param_prototype)))
  params <- params[nzchar(params)]
  rows <- lapply(params, function(param) {
    part <- tab[trimws(as.character(tab$param_prototype)) == param, , drop = FALSE]
    active <- if ("active_in_fit" %in% names(part)) report_truthy(part$active_in_fit) else rep(FALSE, nrow(part))
    estimate <- if ("estimate" %in% names(part)) report_truthy(part$estimate) else rep(FALSE, nrow(part))
    value <- if ("prototype_value" %in% names(part)) part$prototype_value else NA_real_
    rel_dist <- if ("rel_dist_to_nearest" %in% names(part)) part$rel_dist_to_nearest else NA_real_
    bound_status <- if ("bound_status" %in% names(part)) {
      paste(unique(as.character(part$bound_status[active & !is.na(part$bound_status)])), collapse = ", ")
    } else {
      ""
    }
    if (!nzchar(bound_status)) bound_status <- "fixed_or_unavailable"
    data.frame(
      parameter = param,
      parameter_description = first_nonempty(part$parameter_description),
      estimated = any(estimate, na.rm = TRUE),
      active_seed_count = sum(active, na.rm = TRUE),
      seed_count = length(unique(as.character(part$seed))),
      min_value = finite_or_na(value, min),
      median_value = finite_or_na(value, stats::median),
      max_value = finite_or_na(value, max),
      min_rel_dist_to_bound = finite_or_na(rel_dist[active], min),
      boundary_status = bound_status,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(!out$estimated, out$parameter), , drop = FALSE]
  row.names(out) <- NULL
  out
}

report_deoptim_converged_flags <- function(seed_summary) {
  n <- nrow(seed_summary)
  if (!n) return(logical(0))
  de_obj <- if ("optimizer_deoptim_objective" %in% names(seed_summary)) {
    suppressWarnings(as.numeric(seed_summary$optimizer_deoptim_objective))
  } else {
    rep(NA_real_, n)
  }
  if ("objective" %in% names(seed_summary)) {
    objective <- suppressWarnings(as.numeric(seed_summary$objective))
    fallback_idx <- !is.finite(de_obj) & is.finite(objective)
    de_obj[fallback_idx] <- objective[fallback_idx]
  }
  has_objective <- is.finite(de_obj)
  converged <- rep(FALSE, n)
  has_convergence_evidence <- rep(FALSE, n)
  if ("optimizer_interrupted" %in% names(seed_summary)) {
    interrupted <- report_truthy(seed_summary$optimizer_interrupted)
  } else {
    interrupted <- rep(FALSE, n)
  }
  if (length(interrupted) != n) interrupted <- rep(FALSE, n)
  if ("deoptim_stop_reason" %in% names(seed_summary)) {
    reason <- trimws(as.character(seed_summary$deoptim_stop_reason))
    has_reason <- !is.na(reason) & nzchar(reason) & !tolower(reason) %in% c("na", "nan", "null")
    reason_lower <- tolower(reason)
    converged[has_reason] <- reason_lower[has_reason] %in% c("early_stop_reltol_or_steptol")
    has_convergence_evidence[has_reason] <- TRUE
  }

  iter_completed <- if ("optimizer_iter_completed" %in% names(seed_summary)) {
    suppressWarnings(as.numeric(seed_summary$optimizer_iter_completed))
  } else {
    rep(NA_real_, n)
  }
  iter_target <- if ("optimizer_iter_target" %in% names(seed_summary)) {
    suppressWarnings(as.numeric(seed_summary$optimizer_iter_target))
  } else {
    rep(NA_real_, n)
  }
  if ("itermax" %in% names(seed_summary)) {
    itermax <- suppressWarnings(as.numeric(seed_summary$itermax))
    fill_itermax <- !is.finite(iter_target) & is.finite(itermax)
    iter_target[fill_itermax] <- itermax[fill_itermax]
  }
  legacy_target_idx <- !is.finite(iter_target) & is.finite(iter_completed)
  iter_target[legacy_target_idx] <- 500

  has_iter_evidence <- is.finite(iter_completed) & is.finite(iter_target)
  fill_from_iter <- has_iter_evidence & !has_convergence_evidence
  converged[fill_from_iter] <- iter_completed[fill_from_iter] < iter_target[fill_from_iter]
  has_convergence_evidence[fill_from_iter] <- TRUE
  max_iter_reached <- has_iter_evidence & iter_completed >= iter_target
  converged[max_iter_reached] <- FALSE
  has_convergence_evidence[max_iter_reached] <- TRUE

  converged <- converged & has_convergence_evidence & has_objective & !interrupted
  converged[is.na(converged)] <- FALSE
  converged
}

report_fit_label_for_summary <- function(seed_summary) {
  fit_mode <- if ("fit_mode" %in% names(seed_summary)) unique(as.character(seed_summary$fit_mode)) else character(0)
  fit_mode <- fit_mode[!is.na(fit_mode) & nzchar(trimws(fit_mode))]
  if (any(fit_mode == "fit_joint")) return("joint fitting")
  if (any(fit_mode == "fit_invitro")) return("in vitro")
  "in vivo"
}

normalize_convergence_summary_columns <- function(tab) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(tab)
  fit_label <- if ("Fit" %in% names(tab)) trimws(as.character(tab$Fit[[1]])) else ""
  if (!identical(fit_label, "joint fitting")) {
    keep <- intersect(c("Fit", "Total seeds", "DEoptim converged"), names(tab))
    return(tab[, keep, drop = FALSE])
  }
  tab
}

convergence_summary_legend <- function(run_mode) {
  deoptim_text <- "DEoptim converged requires a finite objective plus explicit convergence evidence; reaching the iteration target is counted as not converged."
  if (identical(run_mode, "fit_joint")) {
    return(paste(
      deoptim_text,
      "L-BFGS-B accepted counts seed-level local refinements marked optimizer_local_accepted."
    ))
  }
  deoptim_text
}

report_missing_summary_value <- o2sd_missing_summary_value

supplement_report_optimizer_fields <- function(seed_summary, extra_results_dir) {
  if (is.null(seed_summary) || !is.data.frame(seed_summary) || !nrow(seed_summary)) return(seed_summary)
  if (!("seed" %in% names(seed_summary))) return(seed_summary)
  refinement_path <- file.path(dirname(extra_results_dir), "lbfgsb_refinement_accepted_seeds.csv")
  if (!file.exists(refinement_path)) return(seed_summary)
  refinement <- tryCatch(
    utils::read.csv(refinement_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(refinement) || !is.data.frame(refinement) || !nrow(refinement) || !("seed" %in% names(refinement))) {
    return(seed_summary)
  }
  optimizer_cols <- intersect(
    c(
      "optimizer_deoptim_objective",
      "optimizer_local_objective",
      "optimizer_local_attempted",
      "optimizer_local_accepted",
      "optimizer_local_convergence",
      "optimizer_local_maxit"
    ),
    names(refinement)
  )
  if (!length(optimizer_cols)) return(seed_summary)
  match_idx <- match(as.character(seed_summary$seed), as.character(refinement$seed))
  has_match <- !is.na(match_idx)
  if (!any(has_match)) return(seed_summary)
  for (col in optimizer_cols) {
    if (!(col %in% names(seed_summary))) seed_summary[[col]] <- NA
    replacement <- rep(NA, nrow(seed_summary))
    replacement[has_match] <- refinement[[col]][match_idx[has_match]]
    fill_idx <- has_match &
      report_missing_summary_value(seed_summary[[col]]) &
      !report_missing_summary_value(replacement)
    if (any(fill_idx)) seed_summary[[col]][fill_idx] <- replacement[fill_idx]
  }
  seed_summary
}

build_convergence_summary_table <- function(extra_results_dir) {
  summary_path <- file.path(extra_results_dir, "convergence_summary.tsv")
  summary_tab <- read_table_optional(summary_path, sep = "\t")
  if (is.data.frame(summary_tab) && nrow(summary_tab)) {
    return(normalize_convergence_summary_columns(summary_tab))
  }

  seed_summary <- read_table_optional(file.path(extra_results_dir, "seed_summary.tsv"), sep = "\t")
  if (is.null(seed_summary) || !is.data.frame(seed_summary) || !nrow(seed_summary)) {
    return(NULL)
  }
  seed_summary <- supplement_report_optimizer_fields(seed_summary, extra_results_dir)
  de_converged <- report_deoptim_converged_flags(seed_summary)
  fit_label <- report_fit_label_for_summary(seed_summary)
  if (!identical(fit_label, "joint fitting")) {
    return(data.frame(
      Fit = fit_label,
      `Total seeds` = nrow(seed_summary),
      `DEoptim converged` = sum(de_converged, na.rm = TRUE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  local_accepted <- if ("optimizer_local_accepted" %in% names(seed_summary)) {
    report_truthy(seed_summary$optimizer_local_accepted)
  } else {
    rep(FALSE, nrow(seed_summary))
  }
  both <- de_converged & local_accepted
  data.frame(
    Fit = fit_label,
    `Total seeds` = nrow(seed_summary),
    `DEoptim converged` = sum(de_converged, na.rm = TRUE),
    `L-BFGS-B accepted` = sum(local_accepted, na.rm = TRUE),
    `Converged and accepted` = sum(both, na.rm = TRUE),
    `Converged only` = sum(de_converged & !local_accepted, na.rm = TRUE),
    `Accepted only` = sum(!de_converged & local_accepted, na.rm = TRUE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

convergence_venn_html <- function(extra_results_dir, run_mode = "unknown") {
  if (!identical(run_mode, "fit_joint")) return("")
  venn_fig <- make_figure_spec_optional(
    extra_results_dir,
    "convergence_venn.pdf",
    "Convergence Venn Diagram",
    "Overlap between DEoptim converged seeds and L-BFGS-B accepted seeds."
  )
  if (is.null(venn_fig)) return("")
  iteration_fig <- make_figure_spec_optional(
    extra_results_dir,
    "convergence_iteration_histogram.pdf",
    "DEoptim Iteration Frequency",
    "Frequency distribution of seed-level DEoptim iterations completed before convergence or stopping."
  )
  iteration_html <- if (!is.null(iteration_fig)) {
    paste0(
      '<figure class="report-convergence-figure">',
      figure_media_html(iteration_fig),
      '<figcaption class="report-figure-legend">Frequency distribution of seed-level DEoptim iterations completed before convergence or stopping.</figcaption>',
      '</figure>'
    )
  } else {
    '<div class="report-convergence-blank"></div>'
  }
  paste0(
    '<div class="report-convergence-grid">',
    '<figure class="report-convergence-figure">',
    figure_media_html(venn_fig),
    '<figcaption class="report-figure-legend">Overlap between DEoptim converged seeds and L-BFGS-B accepted seeds.</figcaption>',
    '</figure>',
    iteration_html,
    '</div>'
  )
}

report_title_for_mode <- function(run_mode) {
  if (identical(run_mode, "fit_invitro")) return("In Vitro Extra Results Report")
  if (identical(run_mode, "fit_joint")) return("Joint Fit Extra Results Report")
  "Extra Results Report"
}

joint_chapter_figure_blocks <- function(figure_specs) {
  figure_groups <- figure_layout_groups(figure_specs)
  vapply(figure_groups, function(group) {
      blocks <- vapply(group, function(i) {
        if (is.na(i)) return(blank_figure_block())
        fig <- figure_specs[[i]]
        figure_index <- fig$index %||% i
        figure_label <- figure_display_label(fig, figure_index)
        sprintf(
          paste0(
            '<section class="report-section" id="figure-%d">',
            '%s',
            '<h2 class="report-figure-title">Figure %s. %s</h2>',
            '<p class="report-figure-legend">%s</p>',
            '<p class="report-figure-file"><code>%s</code></p>',
            '</section>'
        ),
          figure_index,
          figure_media_html(fig),
          escape_html(figure_label),
          escape_html(fig$title),
          escape_html(fig$legend),
          escape_html(fig$filename)
      )
    }, character(1))
    if (length(group) > 1L) {
      paste0('<div class="report-figure-grid report-figure-grid--', length(group), '">', paste(blocks, collapse = ""), '</div>')
    } else {
      blocks
    }
  }, character(1))
}

build_report_html_joint <- function(extra_results_dir, run_mode = "fit_joint") {
  run_label <- infer_run_label(extra_results_dir)
  report_title <- report_title_for_mode(run_mode)
  convergence_summary <- build_convergence_summary_table(extra_results_dir)
  has_convergence_summary <- is.data.frame(convergence_summary) && nrow(convergence_summary) > 0L
  parameter_summary <- build_parameter_summary_table(extra_results_dir)
  has_parameter_summary <- is.data.frame(parameter_summary) && nrow(parameter_summary) > 0L
  chapters <- build_joint_figure_chapters(extra_results_dir)

  nav_chapters <- vapply(seq_along(chapters), function(chapter_i) {
    chapter <- chapters[[chapter_i]]
    summary_nav <- character(0)
    if (chapter_i == 1L && isTRUE(has_convergence_summary)) {
      summary_nav <- c(
        summary_nav,
        '<li class="report-nav-item"><a class="report-nav-link report-nav-figure-link" href="#convergence-summary">Convergence Summary</a></li>'
      )
    }
    if (chapter_i == 1L && isTRUE(has_parameter_summary)) {
      summary_nav <- c(
        summary_nav,
        '<li class="report-nav-item"><a class="report-nav-link report-nav-figure-link" href="#parameter-summary">Parameter Summary</a></li>'
      )
    }
      fig_nav <- vapply(chapter$figures, function(fig) {
        figure_label <- figure_display_label(fig, fig$index)
        sprintf(
          '<li class="report-nav-item"><a class="report-nav-link report-nav-figure-link" href="#figure-%d">Figure %s %s</a></li>',
          fig$index,
          escape_html(figure_label),
          escape_html(fig$title)
        )
      }, character(1))
    sprintf(
      paste0(
        '<li class="report-nav-chapter">',
        '<a class="report-nav-link report-nav-chapter-link" href="#chapter-%d">%d. %s</a>',
        '<ul class="report-nav-sublist">%s</ul>',
        '</li>'
      ),
      chapter_i,
      chapter_i,
      escape_html(chapter$title),
      paste(c(summary_nav, fig_nav), collapse = "")
    )
  }, character(1))
  nav_chapters <- c(
    nav_chapters,
    '<li class="report-nav-chapter"><a class="report-nav-link report-nav-chapter-link" href="#run-provenance">4. Run Provenance</a></li>'
  )

  convergence_section <- if (isTRUE(has_convergence_summary)) {
    paste0(
      '<section class="report-section" id="convergence-summary">',
      '<h2 class="report-figure-title">Convergence Summary</h2>',
      '<p class="report-figure-legend">', escape_html(convergence_summary_legend(run_mode)), '</p>',
      table_to_html(convergence_summary, max_rows = 20),
      convergence_venn_html(extra_results_dir, run_mode = run_mode),
      '</section>'
    )
  } else {
    ""
  }

  parameter_section <- if (isTRUE(has_parameter_summary)) {
    paste0(
      '<section class="report-section" id="parameter-summary">',
      '<h2 class="report-figure-title">Parameter Summary</h2>',
      '<p class="report-figure-legend">One row per current active or estimated parameter across all seeds. The parameter_description column comes from the run parameter table snapshot.</p>',
      table_to_html(parameter_summary, max_rows = 200),
      '</section>'
    )
  } else {
    ""
  }
  summary_sections <- paste0(convergence_section, parameter_section)
  provenance_block <- o2sd_run_provenance_html(dirname(extra_results_dir), title = "4. Run Provenance", section_id = "run-provenance")

  chapter_blocks <- vapply(seq_along(chapters), function(chapter_i) {
    chapter <- chapters[[chapter_i]]
    figure_blocks <- paste(joint_chapter_figure_blocks(chapter$figures), collapse = "")
    chapter_intro <- switch(
      chapter$title,
      "Joint Summary" = "Joint objective, boundary risk, and parameter-boundary diagnostics.",
      "In Vivo" = "Cross-seed in vivo prediction figures generated from the joint-fit seed outputs.",
      "In Vitro" = "In vitro summary figures generated from the joint-fit in vitro components and seed outputs.",
      ""
    )
    paste0(
      '<section class="report-chapter" id="chapter-', chapter_i, '">',
      '<div class="report-chapter-heading">',
      '<h2>', chapter_i, '. ', escape_html(chapter$title), '</h2>',
      if (nzchar(chapter_intro)) paste0('<p>', escape_html(chapter_intro), '</p>') else "",
      '</div>',
      if (chapter_i == 1L) summary_sections else "",
      figure_blocks,
      '</section>'
    )
  }, character(1))

  paste0(
    '<!DOCTYPE html>',
    '<html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    '<title>Extra Results Report</title>',
    '<style>',
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    '.report-shell{display:flex;gap:28px;max-width:1600px;margin:0 auto;padding:24px;}',
    '.report-sidebar{position:sticky;top:24px;align-self:flex-start;width:300px;border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);overflow:hidden;}',
    '.report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}',
    '.report-kicker{font-size:11px;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;opacity:0.78;}',
    '.report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}',
    '.report-subtitle{margin-top:4px;font-size:12px;opacity:0.85;}',
    '.report-nav{padding:10px 8px 12px 8px;}',
    '.report-nav-list,.report-nav-sublist{margin:0;padding:0;list-style:none;}',
    '.report-nav-chapter{margin:5px 0 10px 0;}',
    '.report-nav-item{margin:3px 0;}',
    '.report-nav-link{display:block;border-radius:8px;text-decoration:none;color:#17324c;line-height:1.35;}',
    '.report-nav-chapter-link{padding:9px 12px;background:#eef3f8;font-size:14px;font-weight:700;}',
    '.report-nav-figure-link{padding:6px 10px 6px 16px;font-size:12px;font-weight:600;}',
    '.report-nav-link:hover{background:rgba(47,110,164,0.08);}',
    '.report-main{flex:1;min-width:0;max-width:none;}',
    '.report-hero{margin-bottom:20px;padding:18px 20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.report-hero h1{margin:0 0 6px 0;font-size:28px;line-height:1.15;}',
    '.report-meta{margin:0;color:#516274;font-size:14px;}',
    '.report-chapter{margin-bottom:30px;}',
    '.report-chapter-heading{margin:0 0 12px 0;padding:14px 16px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.report-chapter-heading h2{margin:0 0 4px 0;font-size:21px;line-height:1.2;}',
    '.report-chapter-heading p{margin:0;color:#516274;font-size:13px;}',
    '.report-section{margin-bottom:24px;padding:14px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.report-section--blank{background:transparent;border-color:transparent;box-shadow:none;}',
    '.report-figure-grid{display:grid;gap:14px;margin-bottom:24px;align-items:stretch;}',
    '.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}',
    '.report-convergence-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin-top:14px;align-items:start;}',
    '.report-convergence-figure{margin:0;min-width:0;}',
    '.report-convergence-blank{min-height:1px;}',
    '.report-figure-grid .report-section{margin-bottom:0;min-width:0;display:flex;flex-direction:column;}',
    '.report-figure-grid .report-figure{display:block;overflow:visible;}',
    '.report-figure-grid .report-figure-image{width:100%;height:auto;object-fit:contain;}',
    '.report-figure{margin:0 0 8px 0;}',
    '.report-figure-blank{min-height:1px;}',
    '.report-figure-image{display:block;width:100%;max-width:100%;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-figure-object{display:block;width:100%;min-height:680px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-figure-fallback{padding:18px;text-align:center;}',
    '.report-figure-fallback a{color:#2f6ea4;font-weight:600;text-decoration:none;}',
    '.report-figure-fallback a:hover{text-decoration:underline;}',
    '.report-figure-title{margin:0 0 5px 0;font-size:15px;line-height:1.18;}',
    '.report-figure-legend{margin:0 0 5px 0;font-size:11px;line-height:1.35;color:#425365;}',
    '.report-figure-file{margin:0;color:#5f7082;font-size:10px;}',
    '.report-table{width:100%;border-collapse:collapse;font-size:13px;background:#fff;margin-top:12px;}',
    '.report-table th,.report-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}',
    '.report-table th{background:#f7f9fb;font-weight:700;}',
    '.report-empty{color:#657789;font-style:italic;}',
    'code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}',
    '@media (max-width: 991px){.report-shell{display:block;padding:16px;}.report-sidebar{position:relative;top:auto;width:auto;margin-bottom:16px;}.report-main{max-width:none;}.report-figure-grid--2,.report-convergence-grid{grid-template-columns:1fr;}}',
    '</style></head><body>',
    '<div class="report-shell">',
    '<aside class="report-sidebar" aria-label="Figure navigation">',
    '<div class="report-sidebar-header">',
    '<div class="report-kicker">Navigation</div>',
    '<div class="report-title">', escape_html(report_title), '</div>',
    '<div class="report-subtitle">Figure guide for ', escape_html(run_label), '</div>',
    '</div>',
    '<nav class="report-nav"><ul class="report-nav-list">', paste(nav_chapters, collapse = ""), '</ul></nav>',
    '</aside>',
    '<main class="report-main">',
    '<section class="report-hero">',
    '<h1>', escape_html(report_title), '</h1>',
    '<p class="report-meta"><strong>Run:</strong> ', escape_html(run_label), '<br/>',
    '<strong>Detected mode:</strong> ', escape_html(run_mode), '<br/>',
    '<strong>Source directory:</strong> <code>', escape_html(extra_results_dir), '</code><br/>',
    '<strong>Generated at:</strong> ', escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), '</p>',
    '</section>',
    paste(chapter_blocks, collapse = ""),
    provenance_block,
    '</main></div></body></html>'
  )
}

build_report_html <- function(extra_results_dir, figure_specs, run_mode = "unknown") {
  if (identical(run_mode, "fit_joint")) {
    return(build_report_html_joint(extra_results_dir = extra_results_dir, run_mode = run_mode))
  }
  run_label <- infer_run_label(extra_results_dir)
  report_title <- report_title_for_mode(run_mode)
  convergence_summary <- build_convergence_summary_table(extra_results_dir)
  has_convergence_summary <- is.data.frame(convergence_summary) && nrow(convergence_summary) > 0L
  parameter_summary <- build_parameter_summary_table(extra_results_dir)
  has_parameter_summary <- is.data.frame(parameter_summary) && nrow(parameter_summary) > 0L
  summary_nav <- character(0)
  if (isTRUE(has_convergence_summary)) {
    summary_nav <- c(summary_nav, '<li class="report-nav-item"><a class="report-nav-link" href="#convergence-summary">Convergence Summary</a></li>')
  }
  if (isTRUE(has_parameter_summary)) {
    summary_nav <- c(summary_nav, '<li class="report-nav-item"><a class="report-nav-link" href="#parameter-summary">Parameter Summary</a></li>')
  }
  nav_items <- vapply(seq_along(figure_specs), function(i) {
    fig <- figure_specs[[i]]
    figure_label <- figure_display_label(fig, i)
    sprintf(
      '<li class="report-nav-item"><a class="report-nav-link" href="#figure-%d">Figure %s %s</a></li>',
      i,
      escape_html(figure_label),
      escape_html(fig$title)
    )
  }, character(1))
  nav_items <- c(nav_items, '<li class="report-nav-item"><a class="report-nav-link" href="#run-provenance">Run Provenance</a></li>')

  figure_groups <- figure_layout_groups(figure_specs)
  figure_blocks <- vapply(figure_groups, function(group) {
    blocks <- vapply(group, function(i) {
      if (is.na(i)) return(blank_figure_block())
      fig <- figure_specs[[i]]
      figure_label <- figure_display_label(fig, i)
      sprintf(
        paste0(
          '<section class="report-section" id="figure-%d">',
          '%s',
          '<h2 class="report-figure-title">Figure %s. %s</h2>',
          '<p class="report-figure-legend">%s</p>',
          '<p class="report-figure-file"><code>%s</code></p>',
          '</section>'
        ),
        i,
        figure_media_html(fig),
        escape_html(figure_label),
        escape_html(fig$title),
        escape_html(fig$legend),
        escape_html(fig$filename)
      )
    }, character(1))
    if (length(group) > 1L) {
      paste0('<div class="report-figure-grid report-figure-grid--', length(group), '">', paste(blocks, collapse = ""), '</div>')
    } else {
      blocks
    }
  }, character(1))

  convergence_section <- if (isTRUE(has_convergence_summary)) {
    paste0(
      '<section class="report-section" id="convergence-summary">',
      '<h2 class="report-figure-title">Convergence Summary</h2>',
      '<p class="report-figure-legend">', escape_html(convergence_summary_legend(run_mode)), '</p>',
      table_to_html(convergence_summary, max_rows = 20),
      convergence_venn_html(extra_results_dir, run_mode = run_mode),
      '</section>'
    )
  } else {
    ""
  }

  parameter_section <- if (isTRUE(has_parameter_summary)) {
    paste0(
      '<section class="report-section" id="parameter-summary">',
      '<h2 class="report-figure-title">Parameter Summary</h2>',
      '<p class="report-figure-legend">One row per current active or estimated parameter across all seeds. The parameter_description column comes from the run parameter table snapshot.</p>',
      table_to_html(parameter_summary, max_rows = 200),
      '</section>'
    )
  } else {
    ""
  }
  provenance_block <- o2sd_run_provenance_html(dirname(extra_results_dir), title = "Run Provenance", section_id = "run-provenance")

  paste0(
    '<!DOCTYPE html>',
    '<html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    '<title>Extra Results Report</title>',
    '<style>',
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    '.report-shell{display:flex;gap:28px;max-width:1600px;margin:0 auto;padding:24px;}',
    '.report-sidebar{position:sticky;top:24px;align-self:flex-start;width:280px;border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);overflow:hidden;}',
    '.report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}',
    '.report-kicker{font-size:11px;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;opacity:0.78;}',
    '.report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}',
    '.report-subtitle{margin-top:4px;font-size:12px;opacity:0.85;}',
    '.report-nav{padding:10px 8px 12px 8px;}',
    '.report-nav-list{margin:0;padding:0;list-style:none;}',
    '.report-nav-item{margin:4px 0;}',
    '.report-nav-link{display:block;padding:10px 12px;border-radius:8px;text-decoration:none;color:#17324c;font-size:14px;font-weight:600;line-height:1.35;}',
    '.report-nav-link:hover{background:rgba(47,110,164,0.08);}',
    '.report-main{flex:1;min-width:0;max-width:none;}',
    '.report-hero{margin-bottom:20px;padding:18px 20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.report-hero h1{margin:0 0 6px 0;font-size:28px;line-height:1.15;}',
    '.report-meta{margin:0;color:#516274;font-size:14px;}',
    '.report-section{margin-bottom:24px;padding:14px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.report-section--blank{background:transparent;border-color:transparent;box-shadow:none;}',
    '.report-figure-grid{display:grid;gap:14px;margin-bottom:24px;align-items:stretch;}',
    '.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}',
    '.report-convergence-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin-top:14px;align-items:start;}',
    '.report-convergence-figure{margin:0;min-width:0;}',
    '.report-convergence-blank{min-height:1px;}',
    '.report-figure-grid .report-section{margin-bottom:0;min-width:0;display:flex;flex-direction:column;}',
    '.report-figure-grid .report-figure{display:block;overflow:visible;}',
    '.report-figure-grid .report-figure-image{width:100%;height:auto;object-fit:contain;}',
    '.report-figure{margin:0 0 8px 0;}',
    '.report-figure-blank{min-height:1px;}',
    '.report-figure-image{display:block;width:100%;max-width:100%;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-figure-object{display:block;width:100%;min-height:680px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-figure-fallback{padding:18px;text-align:center;}',
    '.report-figure-fallback a{color:#2f6ea4;font-weight:600;text-decoration:none;}',
    '.report-figure-fallback a:hover{text-decoration:underline;}',
    '.report-figure-title{margin:0 0 5px 0;font-size:15px;line-height:1.18;}',
    '.report-figure-legend{margin:0 0 5px 0;font-size:11px;line-height:1.35;color:#425365;}',
    '.report-figure-file{margin:0;color:#5f7082;font-size:10px;}',
    '.report-table{width:100%;border-collapse:collapse;font-size:13px;background:#fff;margin-top:12px;}',
    '.report-table th,.report-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}',
    '.report-table th{background:#f7f9fb;font-weight:700;}',
    '.report-empty{color:#657789;font-style:italic;}',
    'code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}',
    '@media (max-width: 991px){.report-shell{display:block;padding:16px;}.report-sidebar{position:relative;top:auto;width:auto;margin-bottom:16px;}.report-main{max-width:none;}.report-figure-grid--2,.report-convergence-grid{grid-template-columns:1fr;}}',
    '</style></head><body>',
    '<div class="report-shell">',
    '<aside class="report-sidebar" aria-label="Figure navigation">',
    '<div class="report-sidebar-header">',
    '<div class="report-kicker">Navigation</div>',
    '<div class="report-title">', escape_html(report_title), '</div>',
    '<div class="report-subtitle">Figure guide for ', escape_html(run_label), '</div>',
    '</div>',
    '<nav class="report-nav"><ul class="report-nav-list">', paste(c(summary_nav, nav_items), collapse = ""), '</ul></nav>',
    '</aside>',
    '<main class="report-main">',
    '<section class="report-hero">',
    '<h1>', escape_html(report_title), '</h1>',
    '<p class="report-meta"><strong>Run:</strong> ', escape_html(run_label), '<br/>',
    '<strong>Detected mode:</strong> ', escape_html(run_mode), '<br/>',
    '<strong>Source directory:</strong> <code>', escape_html(extra_results_dir), '</code><br/>',
    '<strong>Generated at:</strong> ', escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), '</p>',
    '</section>',
    convergence_section,
    parameter_section,
    paste(figure_blocks, collapse = ""),
    provenance_block,
    '</main></div></body></html>'
  )
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  extra_results_dir <- argv$extra_results_dir
  staged_inputs <- FALSE
  if (is.null(extra_results_dir) || !nzchar(trimws(extra_results_dir))) {
    analysis_dir <- normalizePath(argv$analysis_dir, mustWork = TRUE)
    simulation_dir <- normalizePath(argv$simulation_dir, mustWork = TRUE)
    viz_dir <- normalizePath(argv$viz_dir, mustWork = TRUE)
    required <- c(
      file.path(analysis_dir, "analysis_manifest.tsv"),
      file.path(simulation_dir, "simulation_manifest.tsv"),
      file.path(viz_dir, "visualization_manifest.tsv")
    )
    if (any(!file.exists(required))) stop("Extra-results report requires completed simulation, analysis, and visualization manifests.", call. = FALSE)
    extra_results_dir <- tempfile("extra_results_report_inputs_")
    dir.create(extra_results_dir, recursive = TRUE, showWarnings = FALSE)
    staged_inputs <- TRUE
    on.exit(unlink(extra_results_dir, recursive = TRUE, force = TRUE), add = TRUE)
    for (source_dir in unique(c(simulation_dir, analysis_dir, viz_dir))) {
      files <- list.files(source_dir, full.names = TRUE, recursive = FALSE)
      files <- files[file.info(files)$isdir %in% FALSE]
      if (length(files)) file.copy(files, extra_results_dir, overwrite = TRUE)
    }
  }
  extra_results_dir <- normalizePath(extra_results_dir, mustWork = TRUE)
  report_dir <- normalizePath(argv$report_dir %||% dirname(argv$out_path %||% file.path(extra_results_dir, "extra_results_report.html")), mustWork = FALSE)
  out_path <- normalizePath(
    argv$out_path %||% file.path(report_dir, "extra_results_report.html"),
    mustWork = FALSE
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  run_mode <- read_report_run_mode(extra_results_dir)
  figure_specs <- build_figure_specs(extra_results_dir)
  html <- build_report_html(extra_results_dir = extra_results_dir, figure_specs = figure_specs, run_mode = run_mode)
  writeLines(html, con = out_path, useBytes = TRUE)
  o2sd_inject_report_image_lightbox(out_path)
  utils::write.table(
    data.frame(stage = "report", file = basename(out_path), staged_inputs = staged_inputs, stringsAsFactors = FALSE),
    file.path(dirname(out_path), "report_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )

  message("Wrote extra results report: ", out_path)
}

if (sys.nframe() == 0) {
  main()
}
