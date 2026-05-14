#!/usr/bin/env Rscript

.o2g_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2g_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
rm(.o2g_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

read_table_optional <- function(path, sep = "\t") {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.table(
      path,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

format_numeric_like <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  out <- as.character(x)
  x_trim <- trimws(out)
  numeric_pattern <- "^[-+]?((\\d+\\.?\\d*)|(\\.\\d+))([eE][-+]?\\d+)?$"
  numeric_like <- !is.na(out) & nzchar(x_trim) & grepl(numeric_pattern, x_trim)
  if (!any(numeric_like)) return(out)

  num <- suppressWarnings(as.numeric(x_trim[numeric_like]))
  keep <- is.finite(num)
  if (!any(keep)) return(out)

  out_num <- as.character(num[keep])
  int_like <- abs(num[keep] - round(num[keep])) < 1e-9
  decimal_text <- formatC(num[keep], format = "f", digits = 3)
  sci_nonzero <- !int_like & num[keep] != 0 & grepl("\\.000$", decimal_text)
  decimal <- !int_like & !sci_nonzero

  out_num[int_like] <- format(round(num[keep][int_like]), scientific = FALSE, trim = TRUE, digits = 15)
  out_num[sci_nonzero] <- formatC(num[keep][sci_nonzero], format = "e", digits = 3)
  out_num[decimal] <- formatC(num[keep][decimal], format = "f", digits = 3)

  idx <- which(numeric_like)[keep]
  out[idx] <- out_num
  out
}

table_to_html <- function(df, max_rows = 80) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return("<p class=\"muted\">No rows available.</p>")
  }
  if (nrow(df) > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
  }
  df[] <- lapply(df, format_numeric_like)
  header <- paste0("<th>", html_escape(names(df)), "</th>", collapse = "")
  body <- apply(
    df,
    1L,
    function(row) paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  )
  paste0(
    "<table><thead><tr>", header, "</tr></thead><tbody>",
    paste0(body, collapse = "\n"),
    "</tbody></table>"
  )
}

parameter_description_table_paths <- function(fit_dir) {
  unique(c(
    file.path(fit_dir, "parameter_table_input.csv"),
    file.path(fit_dir, "parameter_table.csv")
  ))
}

annotate_parameter_table_for_report <- function(tab, fit_dir) {
  o2sd_add_parameter_descriptions(
    tab,
    parameter_tables = parameter_description_table_paths(fit_dir)
  )
}

report_truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y", "on")
}

transformed_param_to_natural <- function(x) {
  x <- as.character(x)
  x <- sub("^ivt__", "", x)
  x[x == "delta_lam"] <- "lam_max"
  x <- sub("^log10_", "", x)
  x <- sub("^logit_", "", x)
  x
}

read_current_invitro_parameter_table <- function(fit_dir) {
  candidates <- parameter_description_table_paths(fit_dir)
  for (path in candidates) {
    tab <- read_table_optional(path, sep = ",")
    if (is.data.frame(tab) && nrow(tab) > 0L && "param_symbol" %in% names(tab)) {
      return(tab)
    }
  }
  NULL
}

current_invitro_parameter_sets <- function(fit_dir) {
  tab <- read_current_invitro_parameter_table(fit_dir)
  if (!is.data.frame(tab) || nrow(tab) == 0L || !"param_symbol" %in% names(tab)) {
    return(list(current = character(0), fitted = character(0), fixed = character(0)))
  }
  current <- if ("use_invitro_fit" %in% names(tab)) report_truthy(tab$use_invitro_fit) else report_truthy(tab$estimate)
  fitted <- current & if ("estimate" %in% names(tab)) report_truthy(tab$estimate) else current
  symbols <- trimws(as.character(tab$param_symbol))
  list(
    current = unique(symbols[current & nzchar(symbols)]),
    fitted = unique(symbols[fitted & nzchar(symbols)]),
    fixed = unique(symbols[current & !fitted & nzchar(symbols)])
  )
}

split_parameter_table_for_current_mode <- function(tab,
                                                   fit_dir,
                                                   parameter_col,
                                                   transformed = FALSE) {
  empty <- list(fitted = NULL, fixed = NULL)
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0L || !(parameter_col %in% names(tab))) {
    return(empty)
  }
  sets <- current_invitro_parameter_sets(fit_dir)
  display_name <- if (isTRUE(transformed)) transformed_param_to_natural(tab[[parameter_col]]) else as.character(tab[[parameter_col]])
  if (length(sets$current) > 0L) {
    keep <- display_name %in% sets$current
    tab <- tab[keep, , drop = FALSE]
    display_name <- display_name[keep]
  }
  if (!nrow(tab)) return(empty)

  tab <- annotate_parameter_table_for_report(tab, fit_dir = fit_dir)
  list(
    fitted = tab[display_name %in% sets$fitted, , drop = FALSE],
    fixed = tab[display_name %in% sets$fixed, , drop = FALSE]
  )
}

parameter_section_html <- function(sections) {
  paste0(vapply(sections, function(section) {
    table <- section$table
    if (is.null(table) || !is.data.frame(table) || nrow(table) == 0L) {
      return("")
    }
    paste0("<h3>", html_escape(section$title), "</h3>", table_to_html(table, max_rows = 140))
  }, character(1)), collapse = "")
}

summary_subset <- function(summary_df) {
  if (is.null(summary_df) || !all(c("metric", "value") %in% names(summary_df))) {
    return(summary_df)
  }
  wanted <- c(
    "fit_mode",
    "objective_total",
    "total_loglik",
    "growth_loglik",
    "ploidy_loglik",
    "flow_loglik",
    "growth_loglik_sum",
    "ploidy_loglik_sum",
    "flow_loglik_sum",
    "n_growth",
    "n_growth_observed",
    "n_ploidy_passages",
    "n_kary_cells",
    "n_flow_passages",
    "n_flow_samples",
    "seed",
    "itermax",
    "itermax_requested",
    "itermax_max",
    "de_reltol",
    "de_steptol",
    "NP_requested",
    "NP_used",
    "n_cores_requested",
    "n_cores_used",
    "misseg_loss_survival",
    "parameter_table",
    "fit_objects_dir",
    "flow_density_path"
  )
  idx <- match(wanted, summary_df$metric)
  out <- data.frame(
    metric = wanted,
    value = ifelse(is.na(idx), NA_character_, as.character(summary_df$value[idx])),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$value), , drop = FALSE]
}

figure_title <- function(path) {
  title <- tools::file_path_sans_ext(basename(path))
  title <- gsub("^invitro_", "", title)
  title <- gsub("_", " ", title, fixed = TRUE)
  paste0(toupper(substr(title, 1, 1)), substr(title, 2, nchar(title)))
}

file_to_data_uri <- function(path, mime) {
  if (!file.exists(path)) return("")
  if (requireNamespace("base64enc", quietly = TRUE)) {
    return(base64enc::dataURI(file = path, mime = mime))
  }
  base64_bin <- Sys.which("base64")
  if (!nzchar(base64_bin)) {
    stop("Self-contained HTML report requires the R package 'base64enc' or a system 'base64' command.")
  }
  enc <- suppressWarnings(system2(base64_bin, path, stdout = TRUE, stderr = FALSE))
  paste0("data:", mime, ";base64,", paste(enc, collapse = ""))
}

make_figure_spec <- function(viz_dir, basename, title, legend) {
  pdf_path <- file.path(viz_dir, paste0(basename, ".pdf"))
  png_path <- file.path(viz_dir, paste0(basename, ".png"))
  if (!file.exists(pdf_path) && !file.exists(png_path)) return(NULL)
  list(
    basename = basename,
    title = title,
    legend = legend,
    pdf = if (file.exists(pdf_path)) normalizePath(pdf_path, mustWork = TRUE) else NULL,
    png = if (file.exists(png_path)) normalizePath(png_path, mustWork = TRUE) else NULL
  )
}

optional_figure <- function(viz_dir, basename, title, legend) {
  fig <- make_figure_spec(viz_dir, basename, title, legend)
  if (is.null(fig)) list() else list(fig)
}

build_invitro_section_specs <- function(viz_dir) {
  list(
    list(
      name = "Identifiability Diagnostics",
      figures = c(
        optional_figure(
          viz_dir,
          "invitro_identifiability_diagnostics",
          "Identifiability Diagnostics",
          "Local identifiability output for the current seed when available. If no Jacobian/Hessian was saved, the figure reports an optimizer-population proxy instead."
        )
      )
    ),
    list(
      name = "Rate-Function Diagnostics",
      figures = c(
        optional_figure(
          viz_dir,
          "invitro_constant_external_oxygen",
          "Constant External Oxygen",
          "Constant oxygen settings used by the in vitro runner for each passage condition."
        ),
        optional_figure(
          viz_dir,
          "invitro_rate_function_diagnostics",
          "Rate-Function Diagnostics",
          "Best-fit oxygen, ploidy, death, proliferation, and missegregation rate functions for the current seed."
        ),
        optional_figure(
          viz_dir,
          "invitro_daily_counts",
          "Daily Live-Cell Trajectories",
          "Predicted live-cell trajectories within each passage; selected propagation days are marked."
        ),
        optional_figure(
          viz_dir,
          "invitro_lineage_counts",
          "Selected-Day Live Cells",
          "Predicted live cells at the selected propagation day for each passage condition."
        ),
        optional_figure(
          viz_dir,
          "invitro_lineage_growth",
          "Growth Rate Fit",
          "Observed passage growth is overlaid with the fitted in vitro trajectory."
        ),
        optional_figure(
          viz_dir,
          "invitro_lineage_ploidy",
          "Mean Chromosome Count Fit",
          "Observed mean chromosome counts are compared with the fitted mean chromosome-count trajectory."
        ),
        optional_figure(
          viz_dir,
          "invitro_flow_density",
          "Flow-Density Fit",
          "Observed G0/G1 ploidy-density curves are overlaid with the fitted flow-density prediction."
        ),
        optional_figure(
          viz_dir,
          "invitro_distribution_heatmap",
          "Predicted Ploidy Distribution",
          "Full predicted chromosome-count distribution across in vitro passages."
        )
      )
    )
  )
}

figure_spec_html <- function(fig, section_index, figure_index) {
  title <- fig$title %||% figure_title(fig$basename)
  legend <- html_escape(fig$legend %||% paste0("Figure source: ", fig$basename, ".pdf."))
  label <- sprintf("Figure %d.%d %s", section_index, figure_index, title)
  img_html <- if (!is.null(fig$png) && file.exists(fig$png)) {
    paste0(
      "<img src=\"", file_to_data_uri(fig$png, "image/png"), "\" alt=\"",
      html_escape(label), "\">"
    )
  } else if (!is.null(fig$pdf) && file.exists(fig$pdf)) {
    paste0(
      "<object data=\"", file_to_data_uri(fig$pdf, "application/pdf"),
      "\" type=\"application/pdf\" width=\"100%\" height=\"620px\"></object>"
    )
  } else {
    "<p class=\"muted\">Figure asset missing.</p>"
  }
  pdf_link <- if (!is.null(fig$pdf) && file.exists(fig$pdf)) {
    paste0(
      "<p class=\"figlink\"><a download=\"", html_escape(basename(fig$pdf)),
      "\" href=\"", file_to_data_uri(fig$pdf, "application/pdf"),
      "\">Embedded PDF source</a></p>"
    )
  } else {
    ""
  }
  paste0(
    "<section class=\"figure\"><h3>", html_escape(label), "</h3>",
    img_html,
    "<p class=\"legend\">", legend, "</p>",
    pdf_link,
    "</section>"
  )
}

figure_html <- function(viz_dir, report_dir) {
  if (!dir.exists(viz_dir)) {
    return("<p class=\"muted\">No viz directory found.</p>")
  }
  sections <- build_invitro_section_specs(viz_dir)
  sections <- lapply(sections, function(section) {
    section$figures <- Filter(Negate(is.null), section$figures)
    section
  })
  sections <- Filter(function(section) length(section$figures) > 0L, sections)
  if (length(sections) == 0L) {
    return("<p class=\"muted\">No figures found.</p>")
  }
  paste0(vapply(seq_along(sections), function(i) {
    section <- sections[[i]]
    figures <- paste0(vapply(seq_along(section$figures), function(j) {
      figure_spec_html(section$figures[[j]], i, j)
    }, character(1)), collapse = "\n")
    paste0("<h2>3.", i, " ", html_escape(section$name), "</h2>", figures)
  }, character(1)), collapse = "\n")
}

write_html_report <- function(fit_dir, out_dir, report_basename = "fit_report") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  summary_df <- read_table_optional(file.path(fit_dir, "fit_summary.tsv"), sep = "\t")
  best_params <- split_parameter_table_for_current_mode(
    read_table_optional(file.path(fit_dir, "best_params.tsv"), sep = "\t"),
    fit_dir = fit_dir,
    parameter_col = "parameter",
    transformed = FALSE
  )
  best_params_t <- split_parameter_table_for_current_mode(
    read_table_optional(file.path(fit_dir, "best_params_transformed.tsv"), sep = "\t"),
    fit_dir = fit_dir,
    parameter_col = "transformed_parameter",
    transformed = TRUE
  )
  viz_dir <- file.path(fit_dir, "viz")

  html <- paste0(
    "<!doctype html>\n<html><head><meta charset=\"utf-8\">",
    "<title>In Vitro Fit Report</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:28px;color:#1f2933;background:#fbfbf8;}",
    "h1{font-size:28px;margin-bottom:4px;}h2{margin-top:30px;border-bottom:1px solid #d8d8d0;padding-bottom:6px;}h3{margin-bottom:10px;}",
    ".muted{color:#6b7280}.path{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:#475569;}",
    "table{border-collapse:collapse;width:100%;font-size:12px;background:white;margin:10px 0 22px 0;}th,td{border:1px solid #ddd;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#eef2f3;}",
    ".figure{background:white;border:1px solid #ddd;padding:14px;margin:18px 0;box-shadow:0 1px 2px rgba(0,0,0,0.04);}img{max-width:100%;height:auto;display:block;}.legend{font-size:13px;color:#4b5563;}.figlink{font-size:12px;margin:8px 0 0 0;}",
    "</style></head><body>",
    "<h1>In Vitro Fit Report</h1>",
    "<p class=\"path\">", html_escape(normalizePath(fit_dir, mustWork = FALSE)), "</p>",
    "<p class=\"muted\">Generated: ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p>",
    "<h2>Fit Summary</h2>",
    table_to_html(summary_subset(summary_df), max_rows = 120),
    "<h2>Best Parameters</h2>",
    parameter_section_html(list(
      list(title = "Fitted Parameters", table = best_params$fitted),
      list(title = "Fixed Parameters Used By Current Mode", table = best_params$fixed)
    )),
    "<h2>Transformed Parameters</h2>",
    parameter_section_html(list(
      list(title = "Fitted Parameters", table = best_params_t$fitted),
      list(title = "Fixed Parameters Used By Current Mode", table = best_params_t$fixed)
    )),
    "<h2>Figures</h2>",
    figure_html(viz_dir = viz_dir, report_dir = out_dir),
    "</body></html>"
  )
  out_path <- file.path(out_dir, paste0(report_basename, ".html"))
  writeLines(html, con = out_path, useBytes = TRUE)
  out_path
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% stop(
    "Usage: render_invitro_fit_report.R --fit_dir=/abs/path/to/seed_dir [--out_subdir=report] [--report_basename=fit_report]",
    call. = FALSE
  )
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  out_subdir <- argv$out_subdir %||% "report"
  report_basename <- argv$report_basename %||% "fit_report"
  out_dir <- file.path(fit_dir, out_subdir)
  out_path <- write_html_report(fit_dir, out_dir, report_basename)
  message("In vitro report written to: ", normalizePath(out_path, mustWork = FALSE))
  invisible(out_path)
}

if (sys.nframe() == 0) {
  main()
}
