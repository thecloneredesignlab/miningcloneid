#!/usr/bin/env Rscript

.joint_vis_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
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
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  getwd()
})

SCRIPT_DIR <- normalizePath(.joint_vis_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"),
  local = environment()
)
rm(.joint_vis_script_dir)

`%||%` <- o2sd_null_coalesce

read_comparison_tsv <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing simulation table: ",
      path,
      ". Run the in-vivo and in-vitro simulation producers first.",
      call. = FALSE
    )
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

comparison_finite_rows <- function(df, cols) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(data.frame())
  ok <- rep(TRUE, nrow(df))
  for (col in cols) {
    if (!col %in% names(df)) return(data.frame())
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    ok <- ok & is.finite(df[[col]])
  }
  df[ok, , drop = FALSE]
}

comparison_context_bind <- function(invivo_df, invitro_df, invivo_label, invitro_label) {
  if (is.null(invivo_df) || is.null(invitro_df) || !nrow(invivo_df) || !nrow(invitro_df)) {
    return(data.frame())
  }
  invivo_df$context <- invivo_label
  invitro_df$context <- invitro_label
  common_cols <- intersect(names(invivo_df), names(invitro_df))
  out <- rbind(
    invivo_df[, common_cols, drop = FALSE],
    invitro_df[, common_cols, drop = FALSE]
  )
  out$context <- factor(out$context, levels = c(invivo_label, invitro_label))
  out
}

comparison_palette <- function(levels) {
  levels <- as.character(levels)
  known <- c("2N" = "#1f77b4", "4N" = "#d62728")
  if (all(levels %in% names(known))) return(known[levels])
  stats::setNames(grDevices::hcl.colors(length(levels), palette = "Dark 3"), levels)
}

comparison_state_axis_label <- function(df) {
  endpoint <- suppressWarnings(as.numeric(df$endpoint_value))
  if (length(endpoint) && any(is.finite(endpoint) & endpoint > 10)) {
    "Chromosome Number (N)"
  } else {
    "Ploidy"
  }
}

save_comparison_plot_pair <- function(plot, out_dir, basename, width = 14, height = 7) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(basename, ".pdf"))
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  ggplot2::ggsave(
    pdf_path,
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 180,
    bg = "white"
  )
  c(pdf_path, png_path)
}

generate_invivo_invitro_comparison_figures <- function(invivo_dir,
                                                        invitro_dir,
                                                        out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for joint visualization.", call. = FALSE)
  }
  if (!dir.exists(invivo_dir)) {
    stop("Missing in-vivo simulation directory: ", invivo_dir, call. = FALSE)
  }
  if (!dir.exists(invitro_dir)) {
    stop("Missing in-vitro simulation directory: ", invitro_dir, call. = FALSE)
  }

  read_pair <- function(invivo_file, invitro_file, invivo_label, invitro_label) {
    comparison_context_bind(
      read_comparison_tsv(file.path(invivo_dir, invivo_file)),
      read_comparison_tsv(file.path(invitro_dir, invitro_file)),
      invivo_label,
      invitro_label
    )
  }

  context_basic <- c("In vivo", "In vitro")
  generated <- character(0)

  oxygen_curve <- read_pair(
    "functional_curve_oxygen.tsv",
    "functional_curve_oxygen.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  oxygen_curve <- comparison_finite_rows(oxygen_curve, c("death_rate", "ms_rate"))
  if (nrow(oxygen_curve) > 0L && "cohort" %in% names(oxygen_curve)) {
    oxygen_curve$cohort <- factor(
      oxygen_curve$cohort,
      levels = unique(oxygen_curve$cohort)
    )
    p <- ggplot2::ggplot(
      oxygen_curve,
      ggplot2::aes(x = death_rate, y = ms_rate, color = cohort, group = cohort)
    ) +
      ggplot2::geom_path(linewidth = 1) +
      ggplot2::facet_grid(. ~ context) +
      ggplot2::scale_color_manual(
        values = comparison_palette(levels(oxygen_curve$cohort)),
        drop = FALSE
      ) +
      ggplot2::labs(
        title = paste(
          "In Vivo vs In Vitro:",
          "Stress-Associated Death Rate vs Missegregation Rate"
        ),
        subtitle = "Left: in vivo; right: in vitro.",
        x = "Stress-associated death rate",
        y = "Missegregation rate",
        color = "Cohort"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80")
      )
    generated <- c(
      generated,
      save_comparison_plot_pair(
        p,
        out_dir,
        "invivo_vs_invitro_death_rate_vs_missegregation_rate"
      )
    )
  }

  multi_curve_basic <- read_pair(
    "functional_curve_oxygen_multi_ploidy.tsv",
    "functional_curve_oxygen_multi_ploidy.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  multi_curve_basic <- comparison_finite_rows(
    multi_curve_basic,
    c(
      "ms_rate",
      "misseg_nonviable_daughter_fraction",
      "misseg_nonviable_division_prob"
    )
  )
  if (nrow(multi_curve_basic) > 0L && "cohort" %in% names(multi_curve_basic)) {
    multi_curve_basic$cohort <- factor(
      multi_curve_basic$cohort,
      levels = unique(multi_curve_basic$cohort)
    )
    colors <- comparison_palette(levels(multi_curve_basic$cohort))
    make_ms_plot <- function(value_col, title, y_label, basename) {
      plot_df <- multi_curve_basic
      plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
      plot_df <- comparison_finite_rows(plot_df, c("ms_rate", "value"))
      if (!nrow(plot_df)) return(character(0))
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = ms_rate, y = value, color = cohort)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_grid(. ~ context) +
        ggplot2::scale_color_manual(values = colors, drop = FALSE) +
        ggplot2::labs(
          title = paste0("In Vivo vs In Vitro: ", title),
          subtitle = "Left: in vivo; right: in vitro.",
          x = "Missegregation rate",
          y = y_label,
          color = "Reference state"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          strip.background = ggplot2::element_rect(
            fill = "grey95",
            color = "grey80"
          )
        )
      save_comparison_plot_pair(p, out_dir, basename)
    }
    generated <- c(
      generated,
      make_ms_plot(
        "misseg_nonviable_daughter_fraction",
        "Nonviable Daughter Fraction vs Missegregation Rate",
        "Nonviable daughters / all daughters",
        "invivo_vs_invitro_ms_rate_vs_nonviable_daughter_fraction"
      ),
      make_ms_plot(
        "misseg_nonviable_division_prob",
        "Capped Nonviable Daughter Burden vs Missegregation Rate",
        "min(expected nonviable daughters / division, 1)",
        "invivo_vs_invitro_ms_rate_vs_nonviable_division_probability"
      )
    )
  }

  viability_curve <- read_pair(
    "functional_curve_ploidy.tsv",
    "functional_curve_ploidy.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  viability_curve <- comparison_finite_rows(
    viability_curve,
    c("endpoint_value", "viability_after_ms")
  )
  if (nrow(viability_curve) > 0L) {
    p <- ggplot2::ggplot(
      viability_curve,
      ggplot2::aes(x = endpoint_value, y = viability_after_ms)
    ) +
      ggplot2::geom_line(color = "#2ca02c", linewidth = 1) +
      ggplot2::facet_grid(. ~ context) +
      ggplot2::labs(
        title = paste(
          "In Vivo vs In Vitro:",
          "Ploidy-Dependent Post-Missegregation Survival"
        ),
        subtitle = "Left: in vivo; right: in vitro.",
        x = comparison_state_axis_label(viability_curve),
        y = "Post-missegregation survival"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80")
      )
    generated <- c(
      generated,
      save_comparison_plot_pair(
        p,
        out_dir,
        "invivo_vs_invitro_ploidy_vs_viability_after_ms"
      )
    )
  }

  ploidy_resource <- read_pair(
    "functional_curve_ploidy_by_o2.tsv",
    "functional_curve_ploidy_by_o2.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  ploidy_resource <- comparison_finite_rows(
    ploidy_resource,
    c("endpoint_value", "oxygen_pct", "death_rate", "proliferation_rate")
  )
  if (nrow(ploidy_resource) > 0L) {
    make_ploidy_o2_plot <- function(value_col, title, y_label, basename) {
      plot_df <- ploidy_resource
      plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
      plot_df <- comparison_finite_rows(
        plot_df,
        c("endpoint_value", "oxygen_pct", "value")
      )
      if (!nrow(plot_df)) return(character(0))
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = endpoint_value, y = value, color = oxygen_pct)
      ) +
        ggplot2::geom_point(shape = 15, size = 1.8, alpha = 0.95) +
        ggplot2::facet_grid(. ~ context) +
        ggplot2::scale_color_gradient(
          low = "#2C7BB6",
          high = "#F28E2B",
          name = "Oxygen level"
        ) +
        ggplot2::labs(
          title = paste0("In Vivo vs In Vitro: ", title),
          subtitle = "Left: in vivo; right: in vitro.",
          x = comparison_state_axis_label(plot_df),
          y = y_label
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          strip.background = ggplot2::element_rect(
            fill = "grey95",
            color = "grey80"
          )
        )
      save_comparison_plot_pair(p, out_dir, basename)
    }
    generated <- c(
      generated,
      make_ploidy_o2_plot(
        "death_rate",
        "Ploidy vs Stress-Associated Death Rate by Oxygen Level",
        "Stress-associated death rate",
        "invivo_vs_invitro_ploidy_vs_death_rate_by_o2"
      ),
      make_ploidy_o2_plot(
        "proliferation_rate",
        "Ploidy vs Proliferation Rate by Oxygen Level",
        "Proliferation rate",
        "invivo_vs_invitro_ploidy_vs_proliferation_rate_by_o2"
      )
    )
  }

  multi_curve_resource <- read_pair(
    "functional_curve_oxygen_multi_ploidy.tsv",
    "functional_curve_oxygen_multi_ploidy.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  multi_curve_resource <- comparison_finite_rows(
    multi_curve_resource,
    c("oxygen_pct", "ms_rate", "proliferation_rate", "death_rate")
  )
  if (nrow(multi_curve_resource) > 0L && "cohort" %in% names(multi_curve_resource)) {
    multi_curve_resource$cohort <- factor(
      multi_curve_resource$cohort,
      levels = unique(multi_curve_resource$cohort)
    )
    colors <- comparison_palette(levels(multi_curve_resource$cohort))
    make_o2_plot <- function(value_col, title, y_label, basename) {
      plot_df <- multi_curve_resource
      plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
      plot_df <- comparison_finite_rows(plot_df, c("oxygen_pct", "value"))
      if (!nrow(plot_df)) return(character(0))
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = oxygen_pct, y = value, color = cohort)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_grid(. ~ context) +
        ggplot2::scale_color_manual(values = colors, drop = FALSE) +
        ggplot2::labs(
          title = paste0("In Vivo vs In Vitro: ", title),
          subtitle = "Left: in vivo; right: in vitro.",
          x = "Oxygen level (%)",
          y = y_label,
          color = "Reference state"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          strip.background = ggplot2::element_rect(
            fill = "grey95",
            color = "grey80"
          )
        )
      save_comparison_plot_pair(p, out_dir, basename)
    }
    generated <- c(
      generated,
      make_o2_plot(
        "ms_rate",
        "Oxygen Level vs Missegregation Rate",
        "Missegregation rate",
        "invivo_vs_invitro_oxygen_vs_missegregation_rate_multi_ploidy"
      ),
      make_o2_plot(
        "proliferation_rate",
        "Oxygen Level vs Proliferation Rate",
        "Proliferation rate",
        "invivo_vs_invitro_oxygen_vs_proliferation_rate"
      ),
      make_o2_plot(
        "death_rate",
        "Oxygen Level vs Stress-Associated Death Rate",
        "Stress-associated death rate",
        "invivo_vs_invitro_oxygen_vs_death_rate"
      )
    )
  }

  generated <- unique(normalizePath(generated, mustWork = TRUE))
  if (!length(generated)) {
    stop("No joint comparison figure could be generated from the input tables.", call. = FALSE)
  }
  generated
}

write_joint_visualization_manifest <- function(paths,
                                               out_dir,
                                               invivo_dir,
                                               invitro_dir) {
  info <- file.info(paths)
  manifest <- data.frame(
    artifact = basename(paths),
    relative_path = basename(paths),
    format = tolower(tools::file_ext(paths)),
    bytes = as.numeric(info$size),
    invivo_simulation_dir = normalizePath(invivo_dir, mustWork = TRUE),
    invitro_simulation_dir = normalizePath(invitro_dir, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$artifact), , drop = FALSE]
  manifest_path <- file.path(out_dir, "visualization_manifest.tsv")
  utils::write.table(
    manifest,
    manifest_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  manifest_path
}

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir
  if (!is.null(fit_dir)) {
    fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  }

  invivo_dir <- argv$invivo_simulation_dir %||% (
    if (!is.null(fit_dir)) file.path(fit_dir, "simulation", "invivo") else NULL
  )
  invitro_dir <- argv$invitro_simulation_dir %||% (
    if (!is.null(fit_dir)) file.path(fit_dir, "simulation", "invitro") else NULL
  )
  out_dir <- argv$out_dir %||% argv$viz_dir %||% (
    if (!is.null(fit_dir)) file.path(fit_dir, "viz", "invivo_vs_invitro") else NULL
  )

  if (is.null(invivo_dir) || is.null(invitro_dir) || is.null(out_dir)) {
    stop(
      paste(
        "Usage: plot_invivo_invitro_comparison.R",
        "--fit_dir=/abs/seed",
        "[--invivo_simulation_dir=/abs/path]",
        "[--invitro_simulation_dir=/abs/path]",
        "[--out_dir=/abs/path]"
      ),
      call. = FALSE
    )
  }
  invivo_dir <- normalizePath(invivo_dir, mustWork = FALSE)
  invitro_dir <- normalizePath(invitro_dir, mustWork = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- generate_invivo_invitro_comparison_figures(
    invivo_dir,
    invitro_dir,
    out_dir
  )
  manifest <- write_joint_visualization_manifest(
    paths,
    out_dir,
    invivo_dir,
    invitro_dir
  )
  message(
    "[joint-vis] wrote ",
    length(paths),
    " figure files and manifest ",
    manifest
  )
  invisible(c(paths, manifest))
}

if (sys.nframe() == 0L) {
  main()
}
