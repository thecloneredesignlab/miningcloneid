#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

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
  default_title <- tools::file_path_sans_ext(basename(path))
  title_use <- if (is.null(title) || !length(title) || is.na(title[[1]]) || !nzchar(trimws(title[[1]]))) {
    default_title
  } else {
    as.character(title[[1]])
  }
  legend_use <- if (is.null(legend) || !length(legend) || is.na(legend[[1]]) || !nzchar(trimws(legend[[1]]))) {
    paste0("Figure source: ", basename(path), ".")
  } else {
    as.character(legend[[1]])
  }
  list(
    src = normalizePath(path, mustWork = TRUE),
    title = title_use,
    legend = legend_use
  )
}

optional_figure <- function(viz_dir, filename, title, legend) {
  path <- file.path(viz_dir, filename)
  if (!file.exists(path)) return(list())
  list(make_figure_spec(path, title, legend))
}

optional_series_figures <- function(paths, title_tpl, legend_tpl) {
  if (length(paths) == 0L) return(list())
  lapply(paths, function(path) {
    hz <- extract_horizon_day(path)
    make_figure_spec(
      path,
      sprintf(title_tpl, hz),
      sprintf(legend_tpl, hz)
    )
  })
}

render_pdf_preview_png <- function(src_pdf, dest_png, density = 180) {
  img <- magick::image_read(src_pdf, density = density)
  if (length(img) > 1L) {
    img <- img[1]
  }
  magick::image_write(img, path = dest_png, format = "png")
  normalizePath(dest_png, mustWork = TRUE)
}

report_magick_available <- function() {
  if (identical(Sys.getenv("O2G_REPORT_FORCE_NO_MAGICK", unset = ""), "TRUE")) {
    return(FALSE)
  }
  requireNamespace("magick", quietly = TRUE)
}

report_gs_available <- function() {
  nzchar(Sys.which("gs"))
}

report_base64enc_available <- function() {
  requireNamespace("base64enc", quietly = TRUE)
}

render_pdf_preview_png_gs <- function(src_pdf, dest_png, density = 180) {
  gs_bin <- Sys.which("gs")
  if (!nzchar(gs_bin)) {
    stop("Ghostscript ('gs') was requested for PDF preview rendering but is not available in PATH.")
  }
  src_pdf_use <- normalizePath(src_pdf, mustWork = TRUE)
  dest_png_use <- normalizePath(dest_png, mustWork = FALSE)
  density_use <- suppressWarnings(as.integer(density))
  if (!is.finite(density_use) || density_use <= 0L) density_use <- 180L
  args <- c(
    "-dSAFER",
    "-dBATCH",
    "-dNOPAUSE",
    "-sDEVICE=pngalpha",
    sprintf("-r%d", density_use),
    sprintf("-sOutputFile=%s", shQuote(dest_png_use)),
    shQuote(src_pdf_use)
  )
  out <- suppressWarnings(system2(gs_bin, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "Ghostscript failed while rendering PDF preview for ", src_pdf, ": ",
      paste(out, collapse = "\n")
    )
  }
  if (!file.exists(dest_png_use)) {
    stop("Ghostscript did not create expected PNG preview: ", dest_png_use)
  }
  normalizePath(dest_png_use, mustWork = TRUE)
}

file_to_data_uri <- function(path, mime) {
  if (report_base64enc_available()) {
    return(base64enc::dataURI(file = path, mime = mime))
  }
  base64_bin <- Sys.which("base64")
  if (nzchar(base64_bin)) {
    enc <- tryCatch(
      suppressWarnings(system2(base64_bin, c("-w", "0", path), stdout = TRUE, stderr = TRUE)),
      error = function(e) character()
    )
    if (!length(enc)) {
      enc <- tryCatch(
        suppressWarnings(system2(base64_bin, path, stdout = TRUE, stderr = TRUE)),
        error = function(e) character()
      )
    }
    if (length(enc) > 0L) {
      return(sprintf("data:%s;base64,%s", mime, paste(enc, collapse = "")))
    }
  }
  stop(
    "HTML report fallback requires either the R package 'base64enc' or a system 'base64' command ",
    "when 'magick' is unavailable."
  )
}

pdf_to_data_uri <- function(pdf_path) {
  file_to_data_uri(pdf_path, mime = "application/pdf")
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
  oxygen_predict <- sort_paths_by_horizon(list.files(
    viz_dir,
    pattern = "^predict_o2_timecourse_0_[0-9]+day\\.pdf$",
    full.names = TRUE
  ))
  oxygen_predict <- Filter(function(path) !(extract_horizon_day(path) %in% c(100L, 300L)), oxygen_predict)
  glucose_predict <- sort_paths_by_horizon(list.files(
    viz_dir,
    pattern = "^predict_g_timecourse_0_[0-9]+day\\.pdf$",
    full.names = TRUE
  ))
  glucose_predict <- Filter(function(path) !(extract_horizon_day(path) %in% c(100L, 300L)), glucose_predict)

  burden_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "burden_trend_absolute(real_scale).pdf",
      "Burden Trend Absolute (Real Scale)",
      "Absolute tumor burden shown on the real scale across observation times."
    ),
    optional_figure(
      viz_dir,
      "burden_live_dead_decomposition.pdf",
      "Burden Live/Dead Decomposition",
      "Observed and fitted burden decomposed into live, hypoxia-dead, and buffer-dead components."
    ),
    optional_series_figures(
      burden_predict,
      "Predicted Burden Live/Dead Decomposition (0-%s day)",
      "Forward simulation from day 0 to %s showing the predicted live/dead burden decomposition."
    )
  ))

  ploidy_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "ploidy_weighted_mean_over_time.pdf",
      "Weighted Mean Ploidy Over Time",
      "Weighted mean chromosome-state trajectory over time under the fitted model."
    ),
    optional_figure(
      viz_dir,
      "terminal_ploidy_observed_vs_predicted_violin.pdf",
      "Terminal Ploidy Observed vs Predicted",
      "Observed and predicted terminal ploidy or chromosome-number distributions at harvest."
    ),
    optional_series_figures(
      ploidy_predict,
      "Predicted Curves (0-%s day)",
      "Forward simulation from day 0 to %s summarizing predicted burden and ploidy trajectories."
    )
  ))

  missegregation_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "ms_rate_vs_nonviable_daughter_fraction.pdf",
      "Nonviable Daughter Fraction vs MS Rate",
      "Fraction of all daughter cells that are nonviable because of nullisomy, shown against missegregation rate."
    ),
    optional_figure(
      viz_dir,
      "ms_rate_vs_nonviable_division_probability.pdf",
      "Probability of >=1 Nonviable Daughter vs MS Rate",
      "Per-division probability of producing at least one nonviable daughter because of nullisomy, shown against missegregation rate."
    ),
    optional_figure(
      viz_dir,
      "ploidy_vs_viability_after_ms.pdf",
      "Ploidy vs Viability After MS",
      "Viability modifier after a one-copy-loss missegregation event across the ploidy grid."
    )
  ))

  oxygen_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "o2_target_vs_eff_timecourse.pdf",
      "O2 Target vs Effective Timecourse",
      "Timecourse comparison between the oxygen target and the lagged effective oxygen state used by the model."
    ),
    optional_figure(
      viz_dir,
      "predict_burden_vs_o2.pdf",
      "Predicted Burden vs O2",
      "Forward-simulation burden trajectories plotted against the effective oxygen state."
    ),
    optional_figure(
      viz_dir,
      "ploidy_vs_death_rate_by_o2.pdf",
      "Ploidy vs Death Rate by O2",
      "Death rate across the ploidy grid, colored by oxygen level."
    ),
    optional_figure(
      viz_dir,
      "ploidy_vs_proliferation_rate_by_o2.pdf",
      "Ploidy vs Proliferation Rate by O2",
      "Proliferation rate across the ploidy grid, colored by oxygen level."
    ),
    optional_figure(
      viz_dir,
      "oxygen_vs_missegregation_rate.pdf",
      "Oxygen vs Missegregation Rate",
      "Oxygen-response curve for missegregation rate at the reference ploidy state."
    ),
    optional_figure(
      viz_dir,
      "oxygen_vs_missegregation_rate_multi_ploidy.pdf",
      "Oxygen vs Missegregation Rate Across Reference Ploidy States",
      "Oxygen-response curve for missegregation rate across multiple reference ploidy states."
    ),
    optional_figure(
      viz_dir,
      "oxygen_vs_proliferation_rate.pdf",
      "Oxygen vs Proliferation Rate Across Reference Ploidy States",
      "Oxygen-response curve for the fitted proliferation rate across multiple reference ploidy states."
    ),
    optional_figure(
      viz_dir,
      "oxygen_vs_death_rate.pdf",
      "Oxygen vs Death Rate Across Reference Ploidy States",
      "Oxygen-response curve for the fitted death rate across multiple reference ploidy states."
    ),
    optional_series_figures(
      oxygen_predict,
      "Predicted O2 Timecourse (0-%s day)",
      "Forward simulation from day 0 to %s showing the predicted oxygen target and effective state trajectories."
    )
  ))

  glucose_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "g_target_vs_eff_timecourse.pdf",
      "G Target vs Effective Timecourse",
      "Timecourse comparison between the glucose target and the lagged effective glucose state used by the model."
    ),
    optional_figure(
      viz_dir,
      "g_lag_gap_timecourse.pdf",
      "G Lag Gap Timecourse",
      "Difference between glucose target and effective glucose state over time."
    ),
    optional_figure(
      viz_dir,
      "predict_burden_vs_g.pdf",
      "Predicted Burden vs G",
      "Forward-simulation burden trajectories plotted against the effective glucose state."
    ),
    optional_figure(
      viz_dir,
      "glucose_vs_missegregation_rate.pdf",
      "Glucose vs Missegregation Rate",
      "Glucose-response curve for missegregation rate at the reference ploidy state."
    ),
    optional_figure(
      viz_dir,
      "glucose_vs_missegregation_rate_multi_ploidy.pdf",
      "Glucose vs Missegregation Rate Across Reference Ploidy States",
      "Glucose-response curve for missegregation rate across multiple reference ploidy states."
    ),
    optional_figure(
      viz_dir,
      "glucose_vs_proliferation_rate.pdf",
      "Glucose vs Proliferation Rate Across Reference Ploidy States",
      "Glucose-response curve for the fitted proliferation rate across multiple reference ploidy states."
    ),
    optional_figure(
      viz_dir,
      "glucose_vs_death_rate.pdf",
      "Glucose vs Death Rate Across Reference Ploidy States",
      "Glucose-response curve for the fitted death rate across multiple reference ploidy states."
    ),
    optional_series_figures(
      glucose_predict,
      "Predicted G Timecourse (0-%s day)",
      "Forward simulation from day 0 to %s showing the predicted glucose target and effective state trajectories."
    )
  ))

  sections <- list(
    list(name = "Burden", figures = burden_figs),
    list(name = "Ploidy", figures = ploidy_figs),
    list(name = "Missegregation", figures = missegregation_figs),
    list(name = "Oxygen / O2", figures = oxygen_figs),
    list(name = "Glucose / G", figures = glucose_figs)
  )
  Filter(function(section) length(section$figures) > 0L, sections)
}

stage_assets <- function(section_specs) {
  assets_dir <- file.path(
    tempdir(),
    paste0("o2g_report_assets_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  )
  dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
  use_magick <- report_magick_available()
  use_gs <- !use_magick && report_gs_available()
  for (i in seq_along(section_specs)) {
    if (length(section_specs[[i]]$figures) == 0) next
    for (j in seq_along(section_specs[[i]]$figures)) {
      src <- section_specs[[i]]$figures[[j]]$src
      pdf_stage <- file.path(assets_dir, basename(src))
      if (!file.copy(src, pdf_stage, overwrite = TRUE)) {
        stop("Failed to stage PDF asset: ", src)
      }
      section_specs[[i]]$figures[[j]]$pdf_asset_abs <- normalizePath(pdf_stage, mustWork = TRUE)
      if (use_magick || use_gs) {
        png_stage <- sub("\\.pdf$", ".png", pdf_stage, ignore.case = TRUE)
        png_stage <- if (use_magick) {
          render_pdf_preview_png(pdf_stage, png_stage)
        } else {
          render_pdf_preview_png_gs(pdf_stage, png_stage)
        }
        section_specs[[i]]$figures[[j]]$html_embed_kind <- "img"
        section_specs[[i]]$figures[[j]]$html_asset_uri <- normalizePath(png_stage, mustWork = TRUE)
      } else {
        section_specs[[i]]$figures[[j]]$html_embed_kind <- "pdf_object"
        section_specs[[i]]$figures[[j]]$html_asset_uri <- pdf_to_data_uri(pdf_stage)
      }
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

if (sys.nframe() == 0) {
  main()
}
