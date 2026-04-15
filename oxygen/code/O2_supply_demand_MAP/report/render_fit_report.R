#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(magick))

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

get_script_dir <- function() {
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
}

parse_args <- function(args) {
  out <- list()
  if (length(args) == 0) return(out)
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

is_fit_dir <- function(path) {
  dir.exists(path) &&
    file.exists(file.path(path, "fit_summary.tsv")) &&
    file.exists(file.path(path, "best_params.tsv")) &&
    dir.exists(file.path(path, "viz"))
}

find_fit_dirs_under <- function(root_dir) {
  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  sub_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  fit_sub_dirs <- sub_dirs[vapply(sub_dirs, is_fit_dir, logical(1))]
  if (is_fit_dir(root_dir)) {
    fit_sub_dirs <- c(root_dir, fit_sub_dirs)
  }
  sort(unique(normalizePath(fit_sub_dirs, mustWork = TRUE)))
}

read_fit_summary_selected <- function(fit_dir) {
  path <- file.path(fit_dir, "fit_summary.tsv")
  tab <- read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c("metric" = "character", "value" = "character")
  )
  wanted <- c(
    "optimizer_method",
    "objective",
    "objective_data",
    "objective_prior_raw",
    "objective_prior",
    "objective_burden",
    "objective_ploidy",
    "objective_burden_neg2loglik_raw",
    "objective_ploidy_neg2loglik_raw"
  )
  idx <- match(wanted, tab$metric)
  data.frame(
    metric = wanted,
    value = ifelse(is.na(idx), NA_character_, tab$value[idx]),
    stringsAsFactors = FALSE
  )
}

read_best_params <- function(fit_dir) {
  path <- file.path(fit_dir, "best_params.tsv")
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c("parameter" = "character", "value" = "character")
  )
}

extract_horizon_day <- function(path) {
  base <- basename(path)
  as.integer(sub(".*_0_([0-9]+)day\\.pdf$", "\\1", base))
}

sort_paths_by_horizon <- function(paths) {
  if (length(paths) == 0) return(paths)
  paths[order(vapply(paths, extract_horizon_day, integer(1)))]
}

make_figure_spec <- function(path, title, legend) {
  list(
    src = normalizePath(path, mustWork = TRUE),
    title = title,
    legend = legend
  )
}

render_pdf_preview_png <- function(src_pdf, dest_png, density = 180) {
  img <- magick::image_read(src_pdf, density = density)
  if (length(img) > 1L) {
    img <- img[1]
  }
  magick::image_write(img, path = dest_png, format = "png")
  normalizePath(dest_png, mustWork = TRUE)
}

build_section_specs <- function(fit_dir) {
  viz_dir <- file.path(fit_dir, "viz")

  burden_predict <- sort_paths_by_horizon(list.files(
    viz_dir,
    pattern = "^predict_burden_live_dead_decomposition_0_[0-9]+day\\.pdf$",
    full.names = TRUE
  ))
  ploidy_predict <- sort_paths_by_horizon(list.files(
    viz_dir,
    pattern = "^predict_curves_0_[0-9]+day\\.pdf$",
    full.names = TRUE
  ))

  burden_figs <- Filter(Negate(is.null), c(
    if (file.exists(file.path(viz_dir, "burden_trend_absolute(real_scale).pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "burden_trend_absolute(real_scale).pdf"),
        "Burden Trend Absolute (Real Scale)",
        "Absolute tumor burden shown on the real scale across observation times."
      ))
    },
    if (file.exists(file.path(viz_dir, "burden_live_dead_decomposition.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "burden_live_dead_decomposition.pdf"),
        "Burden Live/Dead Decomposition",
        "Observed and fitted burden decomposed into live, hypoxia-dead, and buffer-dead components."
      ))
    },
    lapply(burden_predict, function(path) {
      hz <- extract_horizon_day(path)
      make_figure_spec(
        path,
        sprintf("Predicted Burden Live/Dead Decomposition (0-%s day)", hz),
        sprintf("Forward simulation from day 0 to %s showing the predicted live/dead burden decomposition.", hz)
      )
    })
  ))

  ploidy_figs <- Filter(Negate(is.null), c(
    if (file.exists(file.path(viz_dir, "ploidy_weighted_mean_over_time.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "ploidy_weighted_mean_over_time.pdf"),
        "Weighted Mean Ploidy Over Time",
        "Weighted mean chromosome-state trajectory over time under the fitted model."
      ))
    },
    if (file.exists(file.path(viz_dir, "terminal_ploidy_observed_vs_predicted_violin.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "terminal_ploidy_observed_vs_predicted_violin.pdf"),
        "Terminal Ploidy Observed vs Predicted",
        "Observed and predicted terminal ploidy/chromosome-number distributions at harvest."
      ))
    },
    lapply(ploidy_predict, function(path) {
      hz <- extract_horizon_day(path)
      make_figure_spec(
        path,
        sprintf("Predicted Curves (0-%s day)", hz),
        sprintf("Forward simulation from day 0 to %s summarizing predicted burden and ploidy trajectories.", hz)
      )
    })
  ))

  missegregation_figs <- Filter(Negate(is.null), c(
    if (file.exists(file.path(viz_dir, "ms_rate_vs_nonviable_division_probability.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "ms_rate_vs_nonviable_division_probability.pdf"),
        "Probability of >=1 Nonviable Daughter vs MS Rate",
        "Nullisomy-only per-division probability of producing at least one nonviable daughter, shown against missegregation rate."
      ))
    },
    if (file.exists(file.path(viz_dir, "ploidy_vs_viability_after_ms.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "ploidy_vs_viability_after_ms.pdf"),
        "Ploidy vs Viability After MS",
        "Viability modifier after a one-copy-loss missegregation event across the ploidy grid."
      ))
    }
  ))

  oxygen_figs <- Filter(Negate(is.null), c(
    if (file.exists(file.path(viz_dir, "ploidy_vs_death_rate_by_o2.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "ploidy_vs_death_rate_by_o2.pdf"),
        "Ploidy vs Death Rate by O2",
        "Death rate across the ploidy grid, colored by oxygen level."
      ))
    },
    if (file.exists(file.path(viz_dir, "ploidy_vs_proliferation_rate_by_o2.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "ploidy_vs_proliferation_rate_by_o2.pdf"),
        "Ploidy vs Proliferation Rate by O2",
        "Proliferation rate across the ploidy grid, colored by oxygen level."
      ))
    },
    if (file.exists(file.path(viz_dir, "oxygen_vs_proliferation_rate.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "oxygen_vs_proliferation_rate.pdf"),
        "Oxygen vs Proliferation Rate",
        "Oxygen-response curve for the fitted proliferation rate."
      ))
    },
    if (file.exists(file.path(viz_dir, "oxygen_vs_missegregation_rate_multi_ploidy.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "oxygen_vs_missegregation_rate_multi_ploidy.pdf"),
        "Oxygen vs Missegregation Rate Across Reference Ploidy States",
        "Oxygen-response curve for missegregation rate across multiple reference ploidy states."
      ))
    },
    if (file.exists(file.path(viz_dir, "oxygen_vs_death_rate.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "oxygen_vs_death_rate.pdf"),
        "Oxygen vs Death Rate",
        "Oxygen-response curve for the fitted death rate."
      ))
    }
  ))

  list(
    list(name = "Burden", figures = burden_figs),
    list(name = "Ploidy", figures = ploidy_figs),
    list(name = "Missegregation", figures = missegregation_figs),
    list(name = "Oxygen / O2", figures = oxygen_figs)
  )
}

stage_assets <- function(section_specs) {
  assets_dir <- file.path(tempdir(), paste0("o2sd_report_assets_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(section_specs)) {
    if (length(section_specs[[i]]$figures) == 0) next
    for (j in seq_along(section_specs[[i]]$figures)) {
      src <- section_specs[[i]]$figures[[j]]$src
      pdf_stage <- file.path(assets_dir, basename(src))
      if (!file.copy(src, pdf_stage, overwrite = TRUE)) {
        stop("Failed to stage PDF asset: ", src)
      }
      png_stage <- sub("\\.pdf$", ".png", pdf_stage, ignore.case = TRUE)
      png_stage <- render_pdf_preview_png(pdf_stage, png_stage)
      section_specs[[i]]$figures[[j]]$pdf_asset_abs <- normalizePath(pdf_stage, mustWork = TRUE)
      section_specs[[i]]$figures[[j]]$html_asset_abs <- normalizePath(png_stage, mustWork = TRUE)
    }
  }
  section_specs
}

render_one_fit <- function(fit_dir, template_path, out_subdir = "reprot", report_basename = "fit_report") {
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  out_dir <- file.path(fit_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  section_specs <- stage_assets(build_section_specs(fit_dir))
  params <- list(
    fit_label = basename(fit_dir),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    selected_summary = read_fit_summary_selected(fit_dir),
    best_params = read_best_params(fit_dir),
    sections = section_specs
  )

  html_out <- render(
    input = template_path,
    output_format = "html_document",
    output_file = paste0(report_basename, ".html"),
    output_dir = out_dir,
    params = params,
    envir = new.env(parent = baseenv()),
    quiet = TRUE
  )
  pdf_out <- render(
    input = template_path,
    output_format = "pdf_document",
    output_file = paste0(report_basename, ".pdf"),
    output_dir = out_dir,
    params = params,
    envir = new.env(parent = baseenv()),
    quiet = TRUE
  )

  html_files_dir <- file.path(out_dir, paste0(report_basename, "_files"))
  if (dir.exists(html_files_dir)) {
    unlink(html_files_dir, recursive = TRUE, force = TRUE)
  }
  reprot_assets_dir <- file.path(out_dir, "assets")
  if (dir.exists(reprot_assets_dir)) {
    unlink(reprot_assets_dir, recursive = TRUE, force = TRUE)
  }

  list(
    fit_dir = fit_dir,
    out_dir = normalizePath(out_dir, mustWork = TRUE),
    html = normalizePath(html_out, mustWork = TRUE),
    pdf = normalizePath(pdf_out, mustWork = TRUE)
  )
}

main <- function() {
  script_dir <- normalizePath(get_script_dir(), mustWork = TRUE)
  template_path <- file.path(script_dir, "fit_result_report.Rmd")
  if (!file.exists(template_path)) {
    stop("Missing report template: ", template_path)
  }

  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  fit_root <- argv$fit_dir %||% stop("Usage: render_fit_report.R --fit_dir=/abs/path/to/fit_or_root [--out_subdir=reprot] [--report_basename=fit_report]")
  fit_root <- normalizePath(fit_root, mustWork = TRUE)
  out_subdir <- argv$out_subdir %||% "reprot"
  report_basename <- argv$report_basename %||% "fit_report"

  fit_dirs <- find_fit_dirs_under(fit_root)
  if (length(fit_dirs) == 0) {
    stop("No fit directories found under: ", fit_root)
  }

  message("Found ", length(fit_dirs), " fit director", if (length(fit_dirs) == 1) "y" else "ies", " under: ", fit_root)
  for (i in seq_along(fit_dirs)) {
    fit_dir <- fit_dirs[[i]]
    message("[", i, "/", length(fit_dirs), "] Rendering report for: ", fit_dir)
    res <- render_one_fit(
      fit_dir = fit_dir,
      template_path = template_path,
      out_subdir = out_subdir,
      report_basename = report_basename
    )
    message("  HTML: ", res$html)
    message("  PDF : ", res$pdf)
  }
}

main()
