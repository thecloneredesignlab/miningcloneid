#!/usr/bin/env Rscript

.o2_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(SCRIPT_DIR, "run_provenance_report.R"), local = environment())
rm(.o2_bootstrap_script_dir)

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

report_seed_token <- function(fit_dir) {
  base <- basename(normalizePath(fit_dir, mustWork = FALSE))
  hit <- regmatches(base, regexpr("seed[0-9]+", base, ignore.case = TRUE))
  if (length(hit) > 0L && nzchar(hit[[1]])) {
    return(tolower(hit[[1]]))
  }

  summary_df <- read_table_optional(file.path(fit_dir, "fit_summary.tsv"), sep = "\t")
  if (is.data.frame(summary_df) && all(c("metric", "value") %in% names(summary_df))) {
    seed_value <- summary_df$value[match("seed", summary_df$metric)]
    seed_value <- suppressWarnings(as.integer(seed_value[[1]] %||% NA_integer_))
    if (is.finite(seed_value)) return(paste0("seed", seed_value))
  }

  NA_character_
}

report_basename_with_seed <- function(report_basename, fit_dir) {
  seed_token <- report_seed_token(fit_dir)
  if (is.na(seed_token) || !nzchar(seed_token)) return(report_basename)
  if (grepl(seed_token, report_basename, fixed = TRUE)) return(report_basename)
  paste(report_basename, seed_token, sep = "_")
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

make_figure_spec <- function(viz_dir, basename, title, legend, display_index = NULL, layout_group = NULL) {
  pdf_path <- file.path(viz_dir, paste0(basename, ".pdf"))
  png_path <- file.path(viz_dir, paste0(basename, ".png"))
  if (!file.exists(pdf_path) && !file.exists(png_path)) return(NULL)
  list(
    basename = basename,
    title = title,
    legend = legend,
    display_index = display_index,
    layout_group = layout_group,
    pdf = if (file.exists(pdf_path)) normalizePath(pdf_path, mustWork = TRUE) else NULL,
    png = if (file.exists(png_path)) normalizePath(png_path, mustWork = TRUE) else NULL
  )
}

optional_figure <- function(viz_dir, basename, title, legend, display_index = NULL, layout_group = NULL) {
  fig <- make_figure_spec(viz_dir, basename, title, legend, display_index = display_index, layout_group = layout_group)
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
          "invitro_o2_selected_live_panels",
          "Assigned Fixed Oxygen and Selected-Day Viable Cells",
          "Branch-aware diagnostic for the in vitro runner. Each cohort is split into control/deprived lineage panels using the same branch-specific x-axis as the aligned growth/chromosome-number/burden composite; repeated lineage passages are not averaged across branches. The upper row shows assigned fixed oxygen and the lower row shows selected-day predicted viable cells.",
          display_index = "2.1/2.4"
        ),
        optional_figure(
          viz_dir,
          "invitro_rate_function_diagnostics",
          "Rate-Function Diagnostics",
          "Best-fit fixed-oxygen, chromosome-number, stress-associated death, proliferation, and missegregation rate functions for the current seed.",
          display_index = "2.2"
        ),
        optional_figure(
          viz_dir,
          "invitro_missegregation_probability_over_passage",
          "Mean Per-Chromosome Missegregation Probability Over Passage",
          "Viable-population-weighted mean per-chromosome missegregation probability across in vitro passage branches, computed from the fitted fixed-oxygen levels and selected-day chromosome-number distributions.",
          display_index = "2.2a"
        ),
        optional_figure(
          viz_dir,
          "invitro_daily_counts",
          "Daily Viable-Cell Trajectories",
          "Predicted viable-cell trajectories split into 2N/control, 2N/deprived, 4N/control, and 4N/deprived panels, with each lineage passage shown as an inset subplot; selected propagation days are marked.",
          display_index = "2.3"
        ),
        optional_figure(
          viz_dir,
          "invitro_growth_ploidy_burden_composite",
          "Aligned Growth, Chromosome-Number, and Burden Fit",
          "Composite in vitro fit view. The 2N and 4N cohort blocks are stacked vertically; each block contains growth rate, chromosome-number quantile, and burden-decomposition rows. Repeated lineage passages are split into branch-specific fixed-oxygen x-axis labels rather than averaged together; growth-rate and chromosome-number lines follow parent-child lineage links, so parallel p10 branches both connect to their shared p9 parent. Observed growth-rate points are drawn at their own branch positions. In vitro burden components are viable/dead fractions normalized by the displayed predicted cell components, so the burden row ranges from 0 to 1. Control and deprived panels use their own passage ranges, with rows aligned within each lineage panel.",
          display_index = "2.5"
        ),
        optional_figure(
          viz_dir,
          "invitro_flow_density",
          "Flow-Density Fit",
          "Observed G0/G1 chromosome-number density curves are overlaid with the fitted flow-density prediction.",
          display_index = "2.6",
          layout_group = "density-distribution"
        ),
        optional_figure(
          viz_dir,
          "invitro_distribution_heatmap",
          "Predicted Chromosome-Number Distribution",
          "Full predicted chromosome-number distribution across in vitro passages.",
          display_index = "2.7",
          layout_group = "density-distribution"
        )
      )
    ),
    list(
      name = "Missegregation-Linked Survival and Death Relationships",
      figures = c(
        optional_figure(
          viz_dir,
          "ms_rate_vs_nonviable_daughter_fraction",
          "Nonviable Daughter Fraction vs Missegregation Rate",
          "Per-division fraction of daughter cells that are nonviable because of missegregation-linked loss, shown against missegregation rate.",
          layout_group = "ms-primary"
        ),
        optional_figure(
          viz_dir,
          "death_rate_vs_missegregation_rate",
          "Stress-Associated Death Rate vs Missegregation Rate",
          "Missegregation-rate curve plotted against the fitted stress-associated death rate at the 2N and 4N reference ploidy states.",
          layout_group = "ms-primary"
        ),
        optional_figure(
          viz_dir,
          "ploidy_vs_viability_after_ms",
          "Ploidy-Dependent Post-Missegregation Survival",
          "Post-missegregation survival after a one-copy-loss event across the ploidy grid.",
          layout_group = "ms-primary"
        )
      )
    )
  )
}

figure_display_number <- function(fig, section_index, figure_index) {
  display_index <- fig$display_index %||% ""
  if (nzchar(display_index)) display_index else sprintf("%d.%d", section_index, figure_index)
}

figure_dom_id <- function(fig, section_index, figure_index) {
  display_number <- figure_display_number(fig, section_index, figure_index)
  paste0("figure-", gsub("[^A-Za-z0-9]+", "-", display_number), "-", gsub("[^A-Za-z0-9]+", "-", fig$basename))
}

figure_spec_html <- function(fig, section_index, figure_index) {
  title <- fig$title %||% figure_title(fig$basename)
  legend <- html_escape(fig$legend %||% paste0("Figure source: ", fig$basename, ".pdf."))
  label <- sprintf("Figure %s %s", figure_display_number(fig, section_index, figure_index), title)
  id <- figure_dom_id(fig, section_index, figure_index)
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
    "<section class=\"figure\" id=\"", html_escape(id), "\"><h3>", html_escape(label), "</h3>",
    img_html,
    "<p class=\"legend\">", legend, "</p>",
    pdf_link,
    "</section>"
  )
}

collect_invitro_sections <- function(viz_dir) {
  if (!dir.exists(viz_dir)) return(list())
  sections <- build_invitro_section_specs(viz_dir)
  sections <- lapply(sections, function(section) {
    section$figures <- Filter(Negate(is.null), section$figures)
    section
  })
  Filter(function(section) length(section$figures) > 0L, sections)
}

figure_layout_groups <- function(figures) {
  groups <- list()
  i <- 1L
  while (i <= length(figures)) {
    group <- figures[[i]]$layout_group %||% ""
    if (nzchar(group)) {
      j <- i
      while (j < length(figures) && identical(figures[[j + 1L]]$layout_group %||% "", group)) {
        j <- j + 1L
      }
      groups <- c(groups, list(seq.int(i, j)))
      i <- j + 1L
    } else {
      groups <- c(groups, list(i))
      i <- i + 1L
    }
  }
  groups
}

figure_html <- function(sections) {
  if (length(sections) == 0L) {
    return("<p class=\"muted\">No figures found.</p>")
  }
  paste0(vapply(seq_along(sections), function(i) {
    section <- sections[[i]]
    section_id <- paste0("figure-section-", i)
    figure_groups <- figure_layout_groups(section$figures)
    figures <- paste0(vapply(figure_groups, function(group) {
      blocks <- paste0(vapply(group, function(j) {
        figure_spec_html(section$figures[[j]], i, j)
      }, character(1)), collapse = "\n")
      if (length(group) > 1L) {
        paste0("<div class=\"figure-grid figure-grid--", length(group), "\">", blocks, "</div>")
      } else {
        blocks
      }
    }, character(1)), collapse = "\n")
    paste0("<h2 id=\"", html_escape(section_id), "\">3.", i, " ", html_escape(section$name), "</h2>", figures)
  }, character(1)), collapse = "\n")
}

report_nav_item <- function(id, label, class = "nav-h2") {
  sprintf(
    "<li><a class=\"%s\" href=\"#%s\">%s</a></li>",
    html_escape(class),
    html_escape(id),
    html_escape(label)
  )
}

build_sidebar_nav <- function(sections) {
  items <- c(
    report_nav_item("fit-summary", "Fit Summary", "nav-h2"),
    report_nav_item("best-parameters", "Best Parameters", "nav-h2"),
    report_nav_item("transformed-parameters", "Transformed Parameters", "nav-h2"),
    report_nav_item("figures", "Figures", "nav-h2"),
    report_nav_item("run-provenance", "Run Provenance", "nav-h2")
  )
  for (i in seq_along(sections)) {
    section <- sections[[i]]
    section_id <- paste0("figure-section-", i)
    items <- c(items, report_nav_item(section_id, paste0("3.", i, " ", section$name), "nav-h3"))
    for (j in seq_along(section$figures)) {
      fig <- section$figures[[j]]
      title <- fig$title %||% figure_title(fig$basename)
      items <- c(
        items,
        report_nav_item(
          figure_dom_id(fig, i, j),
          paste0("Figure ", figure_display_number(fig, i, j), " ", title),
          "nav-h4"
        )
      )
    }
  }
  paste0("<ul>", paste(items, collapse = ""), "</ul>")
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
  viz_dir <- file.path(fit_dir, "viz", "invitro")
  if (!dir.exists(viz_dir)) {
    viz_dir <- file.path(fit_dir, "viz")
  }
  figure_sections <- collect_invitro_sections(viz_dir)
  sidebar_nav <- build_sidebar_nav(figure_sections)
  provenance_block <- o2sd_run_provenance_html(fit_dir, title = "Run Provenance", section_id = "run-provenance")

  html <- paste0(
    "<!doctype html>\n<html><head><meta charset=\"utf-8\">",
    "<title>In Vitro Resource-Stress Fit Report</title>",
    "<style>",
    "html{scroll-behavior:smooth;}body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;color:#1f2933;background:#fbfbf8;}",
    ".report-shell{display:flex;gap:22px;align-items:flex-start;padding:24px;}",
    ".sidebar{position:sticky;top:24px;align-self:flex-start;width:292px;max-height:calc(100vh - 48px);overflow:auto;border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);}",
    ".sidebar-header{padding:14px;background:#24384d;color:#fff;}.sidebar-kicker{font-size:11px;text-transform:uppercase;letter-spacing:.06em;opacity:.8;}.sidebar-title{font-size:17px;font-weight:700;margin-top:3px;}",
    ".sidebar nav{padding:10px 8px 14px 8px;}.sidebar ul{margin:0;padding:0;list-style:none;}.sidebar li{margin:3px 0;}.sidebar a{display:block;border-radius:8px;padding:8px 10px;text-decoration:none;color:#17324c;font-size:13px;line-height:1.25;}.sidebar a:hover{background:rgba(47,110,164,0.09);}.sidebar .nav-h3{padding-left:20px;font-size:12px;color:#384c60;}.sidebar .nav-h4{padding-left:34px;font-size:11px;color:#5a6a78;}",
    ".content{flex:1;min-width:0;max-width:1320px;}",
    "h1{font-size:28px;margin:0 0 4px 0;}h2{margin-top:30px;border-bottom:1px solid #d8d8d0;padding-bottom:6px;scroll-margin-top:24px;}h3{margin:0 0 8px 0;font-size:15px;line-height:1.2;scroll-margin-top:24px;}",
    ".muted{color:#6b7280}.path{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:#475569;}",
    "table{border-collapse:collapse;width:100%;font-size:12px;background:white;margin:10px 0 22px 0;}th,td{border:1px solid #ddd;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#eef2f3;}",
    ".figure{background:white;border:1px solid #ddd;padding:12px;margin:16px 0;box-shadow:0 1px 2px rgba(0,0,0,0.04);scroll-margin-top:24px;}.figure-grid{display:grid;gap:14px;margin:16px 0;align-items:start;}.figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}.figure-grid .figure{margin:0;min-width:0;}img{width:100%;max-width:100%;height:auto;display:block;}.legend{font-size:11px;line-height:1.35;color:#4b5563;margin:8px 0 0 0;}.figlink{font-size:10px;margin:6px 0 0 0;}",
    "@media (max-width: 1100px){.report-shell{display:block;padding:16px;}.sidebar{position:relative;top:auto;width:auto;max-height:none;margin-bottom:18px;}.content{max-width:none;}.figure-grid--2,.figure-grid--3{grid-template-columns:1fr;}}",
    "</style></head><body><div class=\"report-shell\">",
    "<aside class=\"sidebar\"><div class=\"sidebar-header\"><div class=\"sidebar-kicker\">Navigation</div><div class=\"sidebar-title\">In Vitro Resource-Stress Report</div></div><nav>",
    sidebar_nav,
    "</nav></aside><main class=\"content\">",
    "<h1>In Vitro Resource-Stress Fit Report</h1>",
    "<p class=\"path\">", html_escape(normalizePath(fit_dir, mustWork = FALSE)), "</p>",
    "<p class=\"muted\">Generated: ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p>",
    "<h2 id=\"fit-summary\">Fit Summary</h2>",
    table_to_html(summary_subset(summary_df), max_rows = 120),
    "<h2 id=\"best-parameters\">Best Parameters</h2>",
    parameter_section_html(list(
      list(title = "Fitted Parameters", table = best_params$fitted),
      list(title = "Fixed Parameters Used By Current Mode", table = best_params$fixed)
    )),
    "<h2 id=\"transformed-parameters\">Transformed Parameters</h2>",
    parameter_section_html(list(
      list(title = "Fitted Parameters", table = best_params_t$fitted),
      list(title = "Fixed Parameters Used By Current Mode", table = best_params_t$fixed)
    )),
    "<h2 id=\"figures\">Figures</h2>",
    figure_html(figure_sections),
    provenance_block,
    "</main></div></body></html>"
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
  report_basename <- report_basename_with_seed(argv$report_basename %||% "fit_report", fit_dir)
  out_dir <- file.path(fit_dir, out_subdir)
  out_path <- write_html_report(fit_dir, out_dir, report_basename)
  message("In vitro report written to: ", normalizePath(out_path, mustWork = FALSE))
  invisible(out_path)
}

if (sys.nframe() == 0) {
  main()
}
