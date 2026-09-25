#!/usr/bin/env Rscript

# Figure-only fixed-O2 entrypoint. It reads materialized numerical and analysis
# tables and writes PDF figures; it never estimates trajectories or parameters.

.fixo2_vis_script_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      path <- env$ofile
      if (is.null(path)) "" else normalizePath(path, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[basename(frame_files) == "render_fixed_o2_figures.R"]
  if (length(own)) {
    dirname(own[[length(own)]])
  } else {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) {
      dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
    } else {
      normalizePath(getwd(), mustWork = FALSE)
    }
  }
})

.fixo2_vis_workflow_root <- normalizePath(
  file.path(.fixo2_vis_script_dir, "..", ".."),
  mustWork = TRUE
)
sys.source(
  file.path(
    .fixo2_vis_workflow_root,
    "util",
    "o2_supply_demand_map_fixed_o2_utils.R"
  ),
  envir = environment(),
  chdir = TRUE
)
sys.source(
  file.path(
    .fixo2_vis_workflow_root,
    "util",
    "o2_supply_demand_map_fixed_o2_format_utils.R"
  ),
  envir = environment(),
  chdir = TRUE
)
sys.source(
  file.path(
    .fixo2_vis_workflow_root,
    "util",
    "o2_supply_demand_map_fixed_o2_table_utils.R"
  ),
  envir = environment(),
  chdir = TRUE
)
sys.source(
  file.path(.fixo2_vis_script_dir, "fixed_o2_plot_utils.R"),
  envir = environment(),
  chdir = TRUE
)

fixo2_vis_parse_args <- o2ipa_parse_args

fixo2_vis_usage <- function() {
  paste(
    "Usage:",
    "  Rscript render_fixed_o2_figures.R --simulation_dir DIR --analysis_dir DIR [options]",
    "",
    "Required:",
    "  --simulation_dir DIR  Root containing materialized numerical tables.",
    "  --analysis_dir DIR    Root produced by run_fixed_o2_analysis.R.",
    "",
    "Options:",
    "  --vis_dir DIR         Figure output root (default: analysis_dir/figures).",
    "  --parts CSV           attractors,counterfactual,validation,agreement (default: all).",
    "  --plot_dt N           Validation dt to display (default: smallest available).",
    "  --objective_transform identity|log10 (default: identity).",
    "  --help                Print this message.",
    sep = "\n"
  )
}

fixo2_vis_read_tsv <- fo2_read_tsv
fixo2_vis_write_tsv <- fo2_write_tsv

fixo2_vis_find <- function(root, candidates, required = TRUE) {
  paths <- normalizePath(file.path(root, candidates), mustWork = FALSE)
  hit <- paths[file.exists(paths)]
  if (length(hit)) return(hit[[1L]])
  if (required) {
    stop(
      "Required table was not found. Checked: ",
      paste(paths, collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

fixo2_vis_parts <- function(x) {
  parts <- tolower(trimws(strsplit(as.character(x %||% "all"), ",", fixed = TRUE)[[1L]]))
  if ("all" %in% parts) {
    return(c("attractors", "counterfactual", "validation", "agreement"))
  }
  parts <- sub("^counterfactual_trajectories$", "counterfactual", parts)
  allowed <- c("attractors", "counterfactual", "validation", "agreement")
  unknown <- setdiff(parts, allowed)
  if (length(unknown)) {
    stop("Unknown --parts value(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  unique(parts)
}

fixo2_vis_attractors <- function(simulation_dir, analysis_dir, vis_dir) {
  plot_data_path <- fixo2_vis_find(
    analysis_dir,
    c(
      "attractors/tables/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv",
      "tables/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv"
    )
  )
  gap_path <- fixo2_vis_find(
    analysis_dir,
    c(
      "attractors/tables/fixed_o2_attractor_spectral_gap_by_seed.tsv",
      "tables/fixed_o2_attractor_spectral_gap_by_seed.tsv"
    )
  )
  gap_summary_path <- fixo2_vis_find(
    analysis_dir,
    c(
      "attractors/tables/fixed_o2_attractor_spectral_gap_regime_summary.tsv",
      "tables/fixed_o2_attractor_spectral_gap_regime_summary.tsv"
    )
  )
  composite_path <- fixo2_vis_find(
    analysis_dir,
    c(
      "attractors/tables/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv",
      "tables/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv"
    )
  )
  violin_path <- fixo2_vis_find(
    analysis_dir,
    c(
      "attractors/tables/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv",
      "tables/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv"
    )
  )
  attractor_plot_data <- fixo2_vis_read_tsv(plot_data_path)
  gap <- fixo2_vis_read_tsv(gap_path)
  gap_summary <- fixo2_vis_read_tsv(gap_summary_path)
  composite <- fixo2_vis_read_tsv(composite_path)
  violin <- fixo2_vis_read_tsv(violin_path)
  out_dir <- file.path(vis_dir, "attractors")
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  fo2_plot_outputs(attractor_plot_data, out_dir)
  fo2_plot_spectral_gap_outputs(gap, gap_summary, fig_dir)
  fo2_plot_ploidy_gap_reliability_composite(composite, gap_summary, fig_dir)
  fo2_plot_ploidy_gap_reliability_violin(violin, fig_dir)
  data.frame(
    part = "attractors",
    numerical_table = NA_character_,
    analysis_table = plot_data_path,
    figure_dir = normalizePath(fig_dir, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

fixo2_vis_counterfactual <- function(simulation_dir, vis_dir) {
  trajectory_path <- fixo2_vis_find(
    simulation_dir,
    c(
      "counterfactual_trajectories/tables/fixed_o2_counterfactual_trajectories.tsv",
      "counterfactual/tables/fixed_o2_counterfactual_trajectories.tsv",
      "tables/fixed_o2_counterfactual_trajectories.tsv"
    )
  )
  summary_path <- fixo2_vis_find(
    simulation_dir,
    c(
      "counterfactual_trajectories/tables/fixed_o2_counterfactual_summary_by_seed.tsv",
      "counterfactual/tables/fixed_o2_counterfactual_summary_by_seed.tsv",
      "tables/fixed_o2_counterfactual_summary_by_seed.tsv"
    )
  )
  trajectory <- fixo2_vis_read_tsv(trajectory_path)
  summary <- fixo2_vis_read_tsv(summary_path)
  fig_dir <- file.path(vis_dir, "counterfactual_trajectories", "figures")
  cf2_plot(trajectory, summary, fig_dir)
  data.frame(
    part = "counterfactual",
    numerical_table = trajectory_path,
    analysis_table = summary_path,
    figure_dir = normalizePath(fig_dir, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

fixo2_vis_validation <- function(simulation_dir, vis_dir, plot_dt = NA_real_) {
  table_root <- fixo2_vis_find(
    simulation_dir,
    c("simulation/tables", "validation/tables", "tables"),
    required = TRUE
  )
  trajectory_path <- file.path(table_root, "eigen_vs_euler_trajectories.tsv")
  if (!file.exists(trajectory_path)) {
    stop("Validation trajectory table was not found: ", trajectory_path, call. = FALSE)
  }
  trajectories <- fixo2_vis_read_tsv(trajectory_path)
  if (!is.finite(plot_dt)) {
    values <- suppressWarnings(as.numeric(trajectories$dt))
    plot_dt <- min(values[is.finite(values)], na.rm = TRUE)
  }
  fig_dir <- file.path(vis_dir, "validation", "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  vf2_plot(
    trajectories,
    file.path(fig_dir, "expm_vs_euler_mean_ploidy.pdf"),
    plot_dt = plot_dt
  )

  simulation_path <- file.path(table_root, "simulation_replicate_mean_ploidy_trajectories.tsv")
  if (file.exists(simulation_path)) {
    simulation <- fixo2_vis_read_tsv(simulation_path)
    for (method in c("expm", "eigen")) {
      vf2_plot_solution_vs_simulation(
        solution_traj = trajectories,
        simulation_traj = simulation,
        out_path = file.path(
          fig_dir,
          paste0(method, "_vs_simulation_replicates_mean_ploidy.pdf")
        ),
        plot_dt = plot_dt,
        solution_name = method,
        solution_label = paste(method, "analytical"),
        title = paste("Fixed-O2", method, "analytical vs simulation")
      )
    }
  }
  data.frame(
    part = "validation",
    numerical_table = trajectory_path,
    analysis_table = NA_character_,
    figure_dir = normalizePath(fig_dir, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

fixo2_vis_agreement <- function(analysis_dir, vis_dir, objective_transform) {
  agreement_path <- fixo2_vis_find(
    analysis_dir,
    c(
      "analytical_simulation_agreement/tables/agreement_augmented_data.tsv",
      "tables/agreement_augmented_data.tsv"
    )
  )
  agreement <- fixo2_vis_read_tsv(agreement_path)
  out_dir <- file.path(vis_dir, "analytical_simulation_agreement")
  methods <- if ("analytical_method" %in% names(agreement)) {
    sort(unique(as.character(agreement$analytical_method)))
  } else {
    "analytical"
  }
  for (method in methods) {
    method_data <- if ("analytical_method" %in% names(agreement)) {
      agreement[agreement$analytical_method == method, , drop = FALSE]
    } else {
      agreement
    }
    make_scatter_outputs(
      method_data,
      out_dir,
      objective_transform = objective_transform,
      analytical_method = method
    )
    make_analytical_solution_outputs(
      method_data,
      out_dir,
      analytical_method = method,
      seed_metadata = unique(method_data[, intersect(
        c("seed_id", "mode_label", "trajectory_regime"),
        names(method_data)
      ), drop = FALSE])
    )
  }
  data.frame(
    part = "agreement",
    numerical_table = NA_character_,
    analysis_table = agreement_path,
    figure_dir = normalizePath(out_dir, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

fixo2_vis_main <- function(args = fixo2_vis_parse_args()) {
  if (isTRUE(args$help)) {
    cat(fixo2_vis_usage(), "\n")
    return(invisible(NULL))
  }
  simulation_dir <- as.character(args$simulation_dir %||% "")
  analysis_dir <- as.character(args$analysis_dir %||% "")
  if (!nzchar(simulation_dir) || !nzchar(analysis_dir)) {
    stop(
      "--simulation_dir and --analysis_dir are required.\n",
      fixo2_vis_usage(),
      call. = FALSE
    )
  }
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  vis_dir <- normalizePath(
    as.character(args$vis_dir %||% file.path(analysis_dir, "figures")),
    mustWork = FALSE
  )
  dir.create(vis_dir, recursive = TRUE, showWarnings = FALSE)
  parts <- fixo2_vis_parts(args$parts)
  plot_dt <- o2ipa_as_num(args$plot_dt, NA_real_)
  objective_transform <- tolower(as.character(args$objective_transform %||% "identity"))
  if (!objective_transform %in% c("identity", "log10")) {
    stop("--objective_transform must be identity or log10", call. = FALSE)
  }
  rows <- list()
  if ("attractors" %in% parts) {
    rows[[length(rows) + 1L]] <- fixo2_vis_attractors(
      simulation_dir,
      analysis_dir,
      vis_dir
    )
  }
  if ("counterfactual" %in% parts) {
    rows[[length(rows) + 1L]] <- fixo2_vis_counterfactual(
      simulation_dir,
      vis_dir
    )
  }
  if ("validation" %in% parts) {
    rows[[length(rows) + 1L]] <- fixo2_vis_validation(
      simulation_dir,
      vis_dir,
      plot_dt
    )
  }
  if ("agreement" %in% parts) {
    rows[[length(rows) + 1L]] <- fixo2_vis_agreement(
      analysis_dir,
      vis_dir,
      objective_transform
    )
  }
  manifest <- if (length(rows)) do.call(rbind, rows) else data.frame()
  fixo2_vis_write_tsv(manifest, file.path(vis_dir, "fixed_o2_figure_manifest.tsv"))
  message("Completed fixed-O2 figure rendering: ", vis_dir)
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  fixo2_vis_main()
}
