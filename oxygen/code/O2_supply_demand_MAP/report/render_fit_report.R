#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

safe_getwd <- function(fallback = tempdir()) {
  tryCatch(getwd(), error = function(e) fallback)
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
  safe_getwd()
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

report_magick_available <- function() {
  if (identical(Sys.getenv("O2SD_REPORT_FORCE_NO_MAGICK", unset = ""), "TRUE")) {
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

report_pandoc_available <- function() {
  if (identical(Sys.getenv("O2SD_REPORT_FORCE_NO_PANDOC", unset = ""), "TRUE")) {
    return(FALSE)
  }
  isTRUE(rmarkdown::pandoc_available("1.12.3"))
}

report_pdflatex_available <- function() {
  if (identical(Sys.getenv("O2SD_REPORT_FORCE_NO_PDFLATEX", unset = ""), "TRUE")) {
    return(FALSE)
  }
  nzchar(Sys.which("pdflatex"))
}

report_pdf_enabled <- function(argv = list()) {
  arg_value <- argv$render_pdf %||% Sys.getenv("O2SD_RENDER_PDF", unset = "")
  if (!nzchar(arg_value)) {
    return(FALSE)
  }
  tolower(arg_value) %in% c("true", "t", "1", "yes", "y")
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

escape_html <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

fit_report_table_html <- function(tab, class_name = "report-table") {
  if (is.null(tab) || !NROW(tab)) {
    return('<p class="report-empty">No data available.</p>')
  }
  headers <- paste(sprintf("<th>%s</th>", escape_html(names(tab))), collapse = "")
  rows <- vapply(seq_len(nrow(tab)), function(i) {
    vals <- vapply(tab[i, , drop = FALSE], function(x) escape_html(x[[1]]), character(1))
    paste0("<tr>", paste(sprintf("<td>%s</td>", vals), collapse = ""), "</tr>")
  }, character(1))
  paste0(
    '<table class="', class_name, '"><thead><tr>', headers,
    '</tr></thead><tbody>', paste(rows, collapse = ""), "</tbody></table>"
  )
}

fit_report_figure_data_uri <- function(fig) {
  embed_kind <- fig$html_embed_kind %||% "img"
  if (identical(embed_kind, "img")) {
    img_path <- fig$html_asset_abs %||% fig$html_asset_uri
    if (is.character(img_path) && length(img_path) && startsWith(img_path[[1]], "data:")) {
      return(img_path[[1]])
    }
    return(file_to_data_uri(img_path[[1]], mime = "image/png"))
  }
  pdf_uri <- fig$html_asset_uri %||% ""
  if (is.character(pdf_uri) && length(pdf_uri) && nzchar(pdf_uri[[1]]) && startsWith(pdf_uri[[1]], "data:")) {
    return(pdf_uri[[1]])
  }
  pdf_to_data_uri(fig$pdf_asset_abs[[1]])
}

build_fit_report_html <- function(params) {
  sections <- params$sections %||% list()
  nav_items <- c(
    '<li class="report-nav-item"><a class="report-nav-link" href="#report-metadata">Report Metadata</a></li>',
    '<li class="report-nav-item"><a class="report-nav-link" href="#fit-summary">1. Fit Summary</a></li>',
    '<li class="report-nav-item"><a class="report-nav-link" href="#best-parameters">2. Best Parameters</a></li>',
    '<li class="report-nav-item"><a class="report-nav-link" href="#figures">3. Figures</a></li>'
  )
  if (length(sections) > 0L) {
    section_nav <- vapply(seq_along(sections), function(i) {
      slug <- paste0("section-", i)
      sprintf(
        '<li class="report-nav-item"><a class="report-nav-link report-nav-sub" href="#%s">3.%d %s</a></li>',
        slug, i, escape_html(sections[[i]]$name)
      )
    }, character(1))
    nav_items <- c(nav_items, section_nav)
  }

  section_blocks <- if (length(sections) == 0L) {
    '<p class="report-empty">No figures found for this fit.</p>'
  } else {
    paste(vapply(seq_along(sections), function(i) {
      section <- sections[[i]]
      slug <- paste0("section-", i)
      fig_blocks <- paste(vapply(seq_along(section$figures), function(j) {
        fig <- section$figures[[j]]
        fig_title <- fig$title %||% tools::file_path_sans_ext(basename(fig$pdf_asset_abs[[1]]))
        fig_legend <- fig$legend %||% paste0("Figure source: ", basename(fig$pdf_asset_abs[[1]]), ".")
        fig_label <- sprintf("Figure %d.%d %s", i, j, fig_title)
        media <- if (identical(fig$html_embed_kind %||% "img", "pdf_object")) {
          sprintf(
            paste0(
              '<div class="report-figure">',
              '<object data="%s" type="application/pdf" class="report-figure-object">',
              '<div class="report-figure-fallback"><a href="%s">Open PDF figure</a></div>',
              '</object></div>'
            ),
            fit_report_figure_data_uri(fig),
            escape_html(basename(fig$pdf_asset_abs[[1]]))
          )
        } else {
          sprintf(
            '<div class="report-figure"><img src="%s" alt="%s" class="report-figure-image"/></div>',
            fit_report_figure_data_uri(fig),
            escape_html(fig_label)
          )
        }
        paste0(
          '<article class="report-figure-card">',
          media,
          '<p class="report-figure-title"><strong>', escape_html(fig_label), "</strong></p>",
          '<p class="report-figure-legend">', escape_html(fig_legend), "</p>",
          "</article>"
        )
      }, character(1)), collapse = "")
      paste0(
        '<section class="report-section" id="', slug, '">',
        '<h2>3.', i, " ", escape_html(section$name), "</h2>",
        fig_blocks,
        "</section>"
      )
    }, character(1)), collapse = "")
  }

  paste0(
    '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    "<title>O2 Supply-Demand MAP Fit Report</title>",
    "<style>",
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    ".report-shell{display:flex;gap:28px;max-width:1680px;margin:0 auto;padding:24px;}",
    ".report-sidebar{position:sticky;top:24px;align-self:flex-start;width:300px;border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);overflow:hidden;}",
    ".report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}",
    ".report-kicker{font-size:11px;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;opacity:0.78;}",
    ".report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}",
    ".report-subtitle{margin-top:4px;font-size:12px;opacity:0.85;}",
    ".report-nav{padding:10px 8px 12px 8px;}.report-nav-list{margin:0;padding:0;list-style:none;}.report-nav-item{margin:4px 0;}",
    ".report-nav-link{display:block;padding:10px 12px;border-radius:8px;text-decoration:none;color:#17324c;font-size:14px;font-weight:600;line-height:1.35;}",
    ".report-nav-link:hover{background:rgba(47,110,164,0.08);}.report-nav-sub{padding-left:22px;font-size:13px;font-weight:500;}",
    ".report-main{flex:1;min-width:0;max-width:1180px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}",
    ".report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-card h2,.report-section h2{margin-top:0;}",
    ".report-meta{margin:0;color:#516274;font-size:14px;}",
    ".report-table{width:100%;border-collapse:collapse;font-size:14px;border-top:2px solid #1b2a38;border-bottom:2px solid #1b2a38;}",
    ".report-table thead tr{border-bottom:1px solid #1b2a38;}.report-table th,.report-table td{padding:10px 12px;border:0;text-align:left;vertical-align:top;}",
    ".report-table th{background:transparent;font-weight:700;}.report-empty{color:#657789;font-style:italic;}",
    ".report-figure-card{margin-bottom:26px;}.report-figure{margin:0 0 12px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}",
    ".report-figure-image{width:100%;height:auto;display:block;border-radius:6px;}.report-figure-object{width:100%;min-height:680px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}",
    ".report-figure-fallback{padding:18px;text-align:center;}.report-figure-fallback a{color:#2f6ea4;font-weight:600;text-decoration:none;}.report-figure-title,.report-figure-legend{margin:8px 0 0 0;}",
    "@media (max-width: 1100px){.report-shell{display:block;}.report-sidebar{position:static;width:auto;margin-bottom:20px;}}",
    "</style></head><body>",
    '<div class="report-shell">',
    '<aside class="report-sidebar"><div class="report-sidebar-header">',
    '<div class="report-kicker">O2 Supply-Demand MAP</div>',
    '<div class="report-title">Fit Report</div>',
    '<div class="report-subtitle">', escape_html(params$fit_label), "</div></div>",
    '<nav class="report-nav"><ul class="report-nav-list">', paste(nav_items, collapse = ""), "</ul></nav></aside>",
    '<main class="report-main">',
    '<section class="report-card" id="report-metadata"><h1>O2 Supply-Demand MAP Fit Report</h1>',
    '<p class="report-meta"><strong>Fit Label:</strong> ', escape_html(params$fit_label), "<br/>",
    "<strong>Generated At:</strong> ", escape_html(params$generated_at), "</p></section>",
    '<section class="report-card" id="fit-summary"><h2>1. Fit Summary</h2>',
    fit_report_table_html(params$selected_summary),
    "</section>",
    '<section class="report-card" id="best-parameters"><h2>2. Best Parameters</h2>',
    fit_report_table_html(params$best_params),
    "</section>",
    '<section class="report-card" id="figures"><h2>3. Figures</h2>',
    section_blocks,
    "</section></main></div></body></html>"
  )
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
    if (file.exists(file.path(viz_dir, "death_rate_vs_ms_rate_by_ploidy_buffering.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "death_rate_vs_ms_rate_by_ploidy_buffering.pdf"),
        "MS Rate vs Death Rate Across Ploidy States",
        "Buffering-mode relationship between fitted death rate and missegregation rate across 1N-5N reference ploidy states in 0.5N increments."
      ))
    } else if (file.exists(file.path(viz_dir, "ms_rate_vs_death_rate.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "ms_rate_vs_death_rate.pdf"),
        "MS Rate vs Death Rate Across Ploidy States",
        "Relationship between fitted death rate and missegregation rate for 2N and 4N cohorts/reference ploidy states."
      ))
    },
    if (file.exists(file.path(viz_dir, "missegregation_survival_modifier_by_N.pdf"))) {
      list(make_figure_spec(
        file.path(viz_dir, "missegregation_survival_modifier_by_N.pdf"),
        "Missegregation Survival Modifier",
        "Buffering-mode survival modifier across mother ploidy state N for missegregation sizes 1-3."
      ))
    },
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
        section_specs[[i]]$figures[[j]]$html_asset_abs <- normalizePath(png_stage, mustWork = TRUE)
        section_specs[[i]]$figures[[j]]$html_asset_uri <- normalizePath(png_stage, mustWork = TRUE)
      } else {
        section_specs[[i]]$figures[[j]]$html_embed_kind <- "pdf_object"
        section_specs[[i]]$figures[[j]]$html_asset_uri <- pdf_to_data_uri(pdf_stage)
      }
    }
  }
  section_specs
}

render_one_fit <- function(fit_dir, template_path, out_subdir = "reprot", report_basename = "fit_report", render_pdf = FALSE) {
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

  html_out <- file.path(out_dir, paste0(report_basename, ".html"))
  if (report_pandoc_available()) {
    html_out <- render(
      input = template_path,
      output_format = "html_document",
      output_file = paste0(report_basename, ".html"),
      output_dir = out_dir,
      params = params,
      envir = new.env(parent = baseenv()),
      quiet = TRUE
    )
  } else {
    writeLines(build_fit_report_html(params), con = html_out, useBytes = TRUE)
    message("  pandoc unavailable: generated HTML-only fallback report.")
  }
  pdf_out <- NA_character_
  pdf_status <- "disabled"
  if (isTRUE(render_pdf) && report_pandoc_available() && report_pdflatex_available()) {
    pdf_status <- "rendered"
    pdf_out <- render(
      input = template_path,
      output_format = "pdf_document",
      output_file = paste0(report_basename, ".pdf"),
      output_dir = out_dir,
      params = params,
      envir = new.env(parent = baseenv()),
      quiet = TRUE
    )
  } else if (isTRUE(render_pdf) && !report_pandoc_available()) {
    pdf_status <- "skipped_pandoc_unavailable"
    message("  pandoc unavailable: generated HTML report and skipped PDF.")
  } else if (isTRUE(render_pdf)) {
    pdf_status <- "skipped_pdflatex_unavailable"
    message("  pdflatex unavailable: generated HTML report and skipped PDF.")
  } else {
    message("  PDF disabled: generated HTML report only.")
  }

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
    pdf = if (is.na(pdf_out)) NA_character_ else normalizePath(pdf_out, mustWork = TRUE),
    pdf_status = pdf_status
  )
}

main <- function() {
  script_dir <- normalizePath(get_script_dir(), mustWork = TRUE)
  template_path <- file.path(script_dir, "fit_result_report.Rmd")
  if (!file.exists(template_path)) {
    stop("Missing report template: ", template_path)
  }

  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  fit_root <- argv$fit_dir %||% stop("Usage: render_fit_report.R --fit_dir=/abs/path/to/fit_or_root [--out_subdir=reprot] [--report_basename=fit_report] [--render_pdf=TRUE]")
  fit_root <- normalizePath(fit_root, mustWork = TRUE)
  out_subdir <- argv$out_subdir %||% "reprot"
  report_basename <- argv$report_basename %||% "fit_report"
  render_pdf <- report_pdf_enabled(argv)

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
      report_basename = report_basename,
      render_pdf = render_pdf
    )
    message("  HTML: ", res$html)
    if (is.na(res$pdf)) {
      if (identical(res$pdf_status, "disabled")) {
        message("  PDF : skipped (disabled)")
      } else if (identical(res$pdf_status, "skipped_pdflatex_unavailable")) {
        message("  PDF : skipped (pdflatex unavailable)")
      } else {
        message("  PDF : skipped")
      }
    } else {
      message("  PDF : ", res$pdf)
    }
  }
}

main()
