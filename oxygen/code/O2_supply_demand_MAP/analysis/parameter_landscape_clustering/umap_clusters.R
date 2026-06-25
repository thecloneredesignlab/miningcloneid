#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

source(file.path(local_script_dir(), "parameter_landscape_utils.R"))

html_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

html_id <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("(^-|-$)", "", x)
  if (!nzchar(x)) "section" else x
}

file_to_data_uri <- function(path, mime = "image/png") {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Required R package is not installed: base64enc")
  }
  if (!file.exists(path)) stop("Missing asset: ", path)
  paste0("data:", mime, ";base64,", base64enc::base64encode(path))
}

as_html_list <- function(items) {
  paste0("<ul>", paste0("<li>", html_escape(items), "</li>", collapse = ""), "</ul>")
}

table_to_html <- function(df, max_rows = 50L) {
  if (is.null(df) || !nrow(df)) return("<p class=\"report-empty\">No rows.</p>")
  df <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  vals <- as.data.frame(lapply(df, function(x) {
    if (is.numeric(x)) {
      ifelse(is.na(x), "", signif(x, 5))
    } else {
      ifelse(is.na(x), "", as.character(x))
    }
  }), stringsAsFactors = FALSE, check.names = FALSE)
  header <- paste0("<th>", html_escape(names(vals)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(vals)), function(i) {
    paste0("<tr>", paste0("<td>", html_escape(unlist(vals[i, , drop = TRUE], use.names = FALSE)), "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0("<table class=\"report-table\"><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

read_silhouette_selected <- function(root_dir, dataset, source_dir, prefix) {
  path <- file.path(root_dir, dataset, "TablesWclusters", source_dir, paste0(prefix, "_silhouette.csv"))
  if (!file.exists(path)) stop("Missing silhouette table: ", path)
  tab <- read_csv_plain(path)
  tab$selected <- coerce_logical_column(tab$selected, "selected")
  selected <- tab[tab$selected, , drop = FALSE]
  if (nrow(selected) != 1L) stop("Silhouette table must contain exactly one selected row: ", path)
  data.frame(
    `clustering basis` = if (identical(source_dir, "ByInputFeatures")) "standardized input variables" else "UMAP coordinates",
    `selected k` = selected$k,
    `mean silhouette width` = selected$average_silhouette,
    `silhouette sample size` = selected$sample_n,
    check.names = FALSE
  )
}

figure_paths <- function(root_dir, dataset, source_dir, prefix) {
  fig_dir <- file.path(root_dir, dataset, "FiguresWclusters", source_dir)
  list(
    png = file.path(fig_dir, paste0(prefix, ".png")),
    pdf = file.path(fig_dir, paste0(prefix, ".pdf"))
  )
}

figure_card <- function(fig, label, legend) {
  img <- if (file.exists(fig$png)) {
    paste0("<img class=\"report-figure-image\" src=\"", file_to_data_uri(fig$png), "\" alt=\"", html_escape(label), "\">")
  } else {
    paste0("<p class=\"report-empty\">Missing image: ", html_escape(fig$png), "</p>")
  }
  pdf_link <- if (file.exists(fig$pdf)) {
    paste0("<p class=\"report-figure-link\"><a download=\"", html_escape(basename(fig$pdf)), "\" href=\"", file_to_data_uri(fig$pdf, "application/pdf"), "\">Embedded PDF source</a></p>")
  } else {
    ""
  }
  paste0(
    "<article class=\"report-figure-card\">",
    "<div class=\"report-figure\">", img, "</div>",
    "<p class=\"report-figure-title\"><strong>", html_escape(label), "</strong></p>",
    "<p class=\"report-figure-legend\">", html_escape(legend), "</p>",
    pdf_link,
    "</article>"
  )
}

invivo_growth_turnover_o2_cluster_summary <- function(root_dir, prefix, source_dir) {
  coord_path <- file.path(root_dir, "invivo", "TablesWclusters", source_dir, paste0(prefix, "_coordinates.csv"))
  metric_path <- file.path(root_dir, "invivo", "Tables", "invivo_best_params_growth_turnover_100d.csv")
  if (!file.exists(coord_path)) stop("Missing clustered coordinate table: ", coord_path)
  if (!file.exists(metric_path)) stop("Missing growth/turnover support table: ", metric_path)

  coord <- read_csv_plain(coord_path)
  metrics <- read_csv_plain(metric_path)
  req_coord <- c("seed", "cluster_id")
  req_metrics <- c(
    "seed",
    "mean_net_growth_with_misseg_death_0_100d",
    "mean_turnover_rate_with_misseg_death_0_100d",
    "mean_O2_0_100d"
  )
  missing <- c(setdiff(req_coord, names(coord)), setdiff(req_metrics, names(metrics)))
  if (length(missing)) stop("Missing columns for cluster summary: ", paste(unique(missing), collapse = ", "))

  merged <- merge(coord[, req_coord, drop = FALSE], metrics[, req_metrics, drop = FALSE], by = "seed", all.x = TRUE)
  metric_cols <- setdiff(req_metrics, "seed")
  summary <- aggregate(
    merged[, metric_cols, drop = FALSE],
    list(cluster = merged$cluster_id),
    function(x) mean(as.numeric(x), na.rm = TRUE)
  )
  summary$n <- as.integer(table(merged$cluster_id)[summary$cluster])
  summary <- summary[order(summary$cluster), , drop = FALSE]
  data.frame(
    cluster = summary$cluster,
    n = summary$n,
    `mean net growth (0-100 d)` = summary$mean_net_growth_with_misseg_death_0_100d,
    `mean turnover rate (0-100 d)` = summary$mean_turnover_rate_with_misseg_death_0_100d,
    `mean O2 (0-100 d)` = summary$mean_O2_0_100d,
    check.names = FALSE
  )
}

extra_summary_block <- function(root_dir, dataset, spec) {
  summary_type <- spec$extra_summary %||% ""
  if (identical(summary_type, "invivo_growth_turnover_o2_by_input_cluster")) {
    input_tab <- invivo_growth_turnover_o2_cluster_summary(root_dir, spec$prefix, "ByInputFeatures")
    umap_tab <- invivo_growth_turnover_o2_cluster_summary(root_dir, spec$prefix, "ByUMAPCoordinates")
    return(paste0(
      "<h4>Cluster-level summaries</h4>",
      "<p>Values are cluster means computed from the same seed-level growth, turnover, and O2 summaries used in this UMAP.</p>",
      "<div class=\"report-table-grid report-table-grid--2\">",
      "<div><h5>Clusters estimated from the input variables</h5>",
      table_to_html(input_tab, max_rows = 10L),
      "</div><div><h5>Clusters estimated from the UMAP coordinates</h5>",
      table_to_html(umap_tab, max_rows = 10L),
      "</div></div>"
    ))
  }
  ""
}

pair_block <- function(root_dir, dataset, spec, figure_index) {
  source_dirs <- c("ByInputFeatures", "ByUMAPCoordinates")
  source_labels <- c("clusters estimated from the input variables", "clusters estimated from the UMAP coordinates")
  fig_labels <- paste0("Figure ", figure_index, LETTERS[1:2], ". ", spec$title, " - ", source_labels)
  figs <- lapply(source_dirs, function(src) figure_paths(root_dir, dataset, src, spec$prefix))
  silhouettes <- do.call(rbind, lapply(source_dirs, function(src) {
    read_silhouette_selected(root_dir, dataset, src, spec$prefix)
  }))
  legend_base <- paste0(
    spec$legend,
    " For best-fit seed results, color represents the objective value, with lower values shown in blue and higher values shown in yellow. The outline color identifies cluster membership."
  )
  paste0(
    "<section class=\"report-subsection\" id=\"", html_escape(spec$id), "\">",
    "<h3>", html_escape(paste0(figure_index, ". ", spec$title)), "</h3>",
    "<p>", html_escape(spec$description), "</p>",
    if (length(spec$bullets %||% character())) as_html_list(spec$bullets) else "",
    "<h4>Number of clusters selected by the silhouette criterion</h4>",
    table_to_html(silhouettes, max_rows = 10L),
    extra_summary_block(root_dir, dataset, spec),
    "<div class=\"report-figure-grid report-figure-grid--2\">",
    figure_card(figs[[1]], fig_labels[[1]], paste0(legend_base, " Clusters in this panel were estimated from the standardized variables used to construct the UMAP.")),
    figure_card(figs[[2]], fig_labels[[2]], paste0(legend_base, " Clusters in this panel were estimated from the two-dimensional UMAP coordinates.")),
    "</div></section>"
  )
}

method_block <- function(title, items, id) {
  paste0(
    "<section class=\"report-card\" id=\"", html_escape(id), "\">",
    "<h2>", html_escape(title), "</h2>",
    as_html_list(items),
    "</section>"
  )
}

parameter_table <- function(dataset) {
  params <- umap_parameter_set(dataset)
  log_params <- umap_log10_parameter_set(dataset)
  data.frame(
    parameter = params,
    log10_transformed = params %in% log_params,
    stringsAsFactors = FALSE
  )
}

build_nav <- function(sections) {
  paste0(
    "<ul class=\"report-nav-list\">",
    paste0(
      vapply(sections, function(x) {
        cls <- x$class %||% "report-nav-link"
        paste0("<li class=\"report-nav-item\"><a class=\"", cls, "\" href=\"#", html_escape(x$id), "\">", html_escape(x$title), "</a></li>")
      }, character(1)),
      collapse = ""
    ),
    "</ul>"
  )
}

build_report_html <- function(root_dir) {
  invivo_specs <- list(
    list(
      id = "invivo-18-params",
      prefix = "invivo_best_params_umap",
      title = "in vivo best-fit results: 18-parameter UMAP",
      description = "This analysis embeds one best-fit parameter vector from each seed using the 18 in vivo soft-coupling parameters included in the multi-paired analysis.",
      bullets = c(
        "Input rows are the 500 seed-level best parameter vectors.",
        "This report panel does not include the DEoptim initial population points.",
        "Positive-scale parameters are transformed as log10(x); the remaining parameters are retained on their original scale. Each variable is then standardized as z = (x - mean(x)) / sd(x) before UMAP."
      ),
      legend = "Each point represents one seed-level best-fit result. Circles indicate that at least one 1000-day ploidy prediction is at or below 44 chromosomes; diamonds indicate that both the 2N- and 4N-starting predictions exceed 44 chromosomes."
    ),
    list(
      id = "invivo-misseg-growth-turnover",
      prefix = "invivo_best_params_growth_turnover_with_misseg_death_100d_umap",
      title = "in vivo best-fit results: 18 parameters plus net growth and turnover",
      description = "This analysis extends the 18-parameter representation by adding the 0-100 day mean net growth rate and mean turnover rate calculated from each seed's best-fit parameter set.",
      bullets = c(
        "The 0-100 day net growth is the time-weighted mean of lambda_eff - mu_eff - d_misseg, where d_misseg is the chromosome missegregation-associated death rate.",
        "The 0-100 day turnover rate is the time-weighted mean of lambda_eff + mu_eff + d_misseg; therefore, it includes both oxygen-driven death and chromosome missegregation-associated death.",
        "The report intentionally excludes the turnover version that contains only O2-driven death."
      ),
      legend = "Each point represents one seed-level best-fit result. Circles indicate that at least one 1000-day ploidy prediction is at or below 44 chromosomes; diamonds indicate that both the 2N- and 4N-starting predictions exceed 44 chromosomes."
    ),
    list(
      id = "invivo-misseg-growth-turnover-o2",
      prefix = "invivo_best_params_growth_turnover_with_misseg_death_O2_100d_umap",
      title = "in vivo best-fit results: 18 parameters plus growth, turnover, and mean O2",
      description = "This analysis adds the 0-100 day mean O2 concentration to the parameter, net-growth, and turnover representation.",
      bullets = c(
        "The mean O2 concentration is averaged over the same 0-100 day simulation window as net growth and turnover.",
        "Turnover remains the version that includes missegregation-associated death.",
        "Mean O2 is included as an additional continuous variable and is standardized together with the other UMAP inputs."
      ),
      extra_summary = "invivo_growth_turnover_o2_by_input_cluster",
      legend = "Each point represents one seed-level best-fit result. Circles indicate that at least one 1000-day ploidy prediction is at or below 44 chromosomes; diamonds indicate that both the 2N- and 4N-starting predictions exceed 44 chromosomes."
    ),
    list(
      id = "invivo-misseg-growth-turnover-o2-ploidy-ratio",
      prefix = "invivo_best_params_growth_turnover_with_misseg_death_O2_ploidy_ratio_100d_umap",
      title = "in vivo best-fit results: 18 parameters plus growth, turnover, O2, and ploidy ratios",
      description = "This analysis adds two 1000-day ploidy-ratio variables to the parameter, net-growth, turnover, and mean-O2 representation.",
      bullets = c(
        "pred1000_ploidy_ratio_2N is the day-1000 predicted ploidy for the 2N-starting population divided by its initial ploidy.",
        "pred1000_ploidy_ratio_4N is the day-1000 predicted ploidy for the 4N-starting population divided by its initial ploidy.",
        "Both ratios are included as additional continuous variables and are standardized together with the other UMAP inputs."
      ),
      legend = "Each point represents one seed-level best-fit result. Circles indicate that at least one 1000-day ploidy prediction is at or below 44 chromosomes; diamonds indicate that both the 2N- and 4N-starting predictions exceed 44 chromosomes."
    )
  )

  invitro_specs <- list(
    list(
      id = "invitro-best-only",
      prefix = "invitro_best_params_umap",
      title = "in vitro best-fit results: 14-parameter UMAP",
      description = "This analysis embeds one in vitro best-fit parameter vector from each seed.",
      bullets = c(
        "Input rows are the 500 seed-level best parameter vectors.",
        "All points are seed-level best-fit results and are shown as circles.",
        "The 14 input parameters are transformed and standardized as described in the Methods."
      ),
      legend = "Each point represents one seed-level best-fit result. Point color encodes objective."
    )
  )

  nav <- build_nav(list(
    list(id = "report-metadata", title = "Overview"),
    list(id = "methods", title = "Methods"),
    list(id = "invivo", title = "in vivo"),
    list(id = "invivo-parameters", title = "in vivo parameters", class = "report-nav-link report-nav-h3"),
    list(id = "invivo-figures", title = "in vivo UMAP figures", class = "report-nav-link report-nav-h3"),
    list(id = "invitro", title = "in vitro"),
    list(id = "invitro-parameters", title = "in vitro parameters", class = "report-nav-link report-nav-h3"),
    list(id = "invitro-figures", title = "in vitro UMAP figures", class = "report-nav-link report-nav-h3")
  ))

  invivo_figures <- paste0(vapply(seq_along(invivo_specs), function(i) {
    pair_block(root_dir, "invivo", invivo_specs[[i]], i)
  }, character(1)), collapse = "")
  invitro_figures <- paste0(vapply(seq_along(invitro_specs), function(i) {
    pair_block(root_dir, "invitro", invitro_specs[[i]], length(invivo_specs) + i)
  }, character(1)), collapse = "")

  css <- paste0(
    "html{scroll-behavior:smooth;}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f4f7fa;color:#1b2a38;}",
    ".report-shell{display:flex;gap:28px;max-width:1680px;margin:0 auto;padding:24px;}.report-sidebar{position:sticky;top:24px;align-self:flex-start;width:310px;max-height:calc(100vh - 48px);border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);overflow:auto;scrollbar-gutter:stable;}",
    ".report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}.report-kicker{font-size:11px;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;opacity:0.78;}.report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}.report-subtitle{margin-top:4px;font-size:12px;opacity:0.85;}",
    ".report-nav{padding:10px 8px 12px 8px;}.report-nav-list{margin:0;padding:0;list-style:none;}.report-nav-item{margin:4px 0;}.report-nav-link{display:block;padding:10px 12px;border-radius:8px;text-decoration:none;color:#17324c;font-size:14px;font-weight:600;line-height:1.35;}.report-nav-link:hover{background:rgba(47,110,164,0.08);}.report-nav-h3{padding-left:26px;font-size:13px;font-weight:500;color:#536577;}",
    ".report-main{flex:1;min-width:0;max-width:1180px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}.report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-card h2,.report-section h2{margin-top:0;}.report-card h2,.report-section h2,.report-subsection h3{scroll-margin-top:24px;}",
    ".report-meta{margin:0;color:#516274;font-size:14px;}.report-empty{color:#657789;font-style:italic;}.report-table{width:100%;border-collapse:collapse;font-size:13px;margin:10px 0 18px 0;}.report-table th,.report-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}.report-table th{background:#f7f9fb;font-weight:700;}",
    ".report-subsection{margin:22px 0 34px 0;padding-top:2px;}.report-subsection h3{margin:0 0 8px 0;font-size:19px;color:#17324c;}.report-subsection h4{margin:16px 0 4px 0;font-size:14px;color:#284662;}ul{margin-top:8px;}.report-figure-grid{display:grid;gap:18px;margin:16px 0 24px 0;align-items:stretch;}.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}",
    ".report-table-grid{display:grid;gap:18px;margin:10px 0 18px 0;align-items:start;}.report-table-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-table-grid h5{margin:0 0 6px 0;font-size:13px;color:#284662;}",
    ".report-figure-card{min-width:0;}.report-figure{margin:0 0 10px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;aspect-ratio:1/1;display:flex;align-items:center;justify-content:center;overflow:hidden;}.report-figure-image{width:100%;height:100%;display:block;object-fit:contain;border-radius:6px;}.report-figure-title,.report-figure-legend{margin:8px 0 0 0;}.report-figure-title{font-size:13px;line-height:1.35;}.report-figure-legend{font-size:12px;line-height:1.4;color:#536577;}.report-figure-link{font-size:11px;margin:6px 0 0 0;}.report-figure-link a{color:#2f6ea4;text-decoration:none;font-weight:600;}",
    "@media (max-width: 1100px){.report-shell{display:block;}.report-sidebar{position:static;width:auto;margin-bottom:20px;}}@media (max-width: 900px){.report-figure-grid--2{grid-template-columns:1fr;}}"
  )

  methods <- c(
    "Parameters on positive multiplicative scales are transformed as log10(x) before analysis; all other variables are retained on their original scale.",
    "For each UMAP input column x_j, values are standardized as z_ij = (x_ij - mean(x_j)) / sd(x_j). This places all variables on a common scale before UMAP and before clustering in the input-variable space.",
    "The UMAPs shown in this report exclude all PCA-based UMAPs. For in vivo, O2-only turnover UMAPs are also excluded; only turnover that includes missegregation-associated death is reported.",
    "Each displayed UMAP is annotated with two clustering analyses. One estimates clusters from the standardized variables used to construct the UMAP, and the other estimates clusters from the final two-dimensional UMAP coordinates.",
    "For k = 2,...,8, k-means clustering is evaluated by the mean silhouette width. The selected k is the value with the largest mean silhouette width. For large data sets, this criterion is evaluated on a fixed random subset of at most 5000 observations, and the selected k is then used to assign clusters for all observations.",
    "Cluster color denotes cluster membership and does not modify the objective-value color scale."
  )

  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>UMAPs Report</title>",
    "<style>", css, "</style></head><body>",
    "<div class=\"report-shell\">",
    "<aside class=\"report-sidebar\"><div class=\"report-sidebar-header\"><div class=\"report-kicker\">parameter landscape clustering</div>",
    "<div class=\"report-title\">UMAPs Report</div><div class=\"report-subtitle\">in vivo and in vitro</div></div><nav class=\"report-nav\">",
    nav,
    "</nav></aside><main class=\"report-main\">",
    "<section class=\"report-card\" id=\"report-metadata\"><h1>UMAPs Report</h1>",
    "<p class=\"report-meta\"><strong>Results root:</strong> ", html_escape(normalizePath(root_dir, mustWork = FALSE)), "<br>",
    "<strong>Generated at:</strong> ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p></section>",
    method_block("Methods", methods, "methods"),
    "<section class=\"report-section\" id=\"invivo\"><h2>in vivo</h2>",
    "<section id=\"invivo-parameters\"><h3>in vivo UMAP input parameters</h3>",
    "<p>The base in vivo UMAP uses 18 multi-paired soft-coupling parameters. The table records whether each parameter is log10-transformed before standardization.</p>",
    table_to_html(parameter_table("invivo"), max_rows = 50L), "</section>",
    "<section id=\"invivo-figures\"><h3>in vivo clustered UMAPs</h3>", invivo_figures, "</section></section>",
    "<section class=\"report-section\" id=\"invitro\"><h2>in vitro</h2>",
    "<section id=\"invitro-parameters\"><h3>in vitro UMAP input parameters</h3>",
    "<p>The in vitro UMAPs use 14 soft-coupling parameters and do not include any ploidy trajectory-derived features.</p>",
    table_to_html(parameter_table("invitro"), max_rows = 50L), "</section>",
    "<section id=\"invitro-figures\"><h3>in vitro clustered UMAPs</h3>", invitro_figures, "</section></section>",
    "</main></div></body></html>"
  )
}

write_umap_cluster_report <- function(result_root = default_parameter_landscape_clustering_dir(),
                                      output_html = file.path(result_root, "parameter_landscape_clustering_umap_cluster_report.html")) {
  result_root <- normalizePath(path.expand(result_root), mustWork = FALSE)
  output_html <- normalizePath(path.expand(output_html), mustWork = FALSE)
  if (!dir.exists(result_root)) stop("Result root does not exist: ", result_root)
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  html <- build_report_html(result_root)
  writeLines(html, con = output_html, useBytes = TRUE)
  message("Wrote UMAP cluster report: ", output_html)
  invisible(output_html)
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
result_root <- argv$result_root %||% default_parameter_landscape_clustering_dir()
output_html <- argv$output_html %||% file.path(result_root, "parameter_landscape_clustering_umap_cluster_report.html")
write_umap_cluster_report(result_root = result_root, output_html = output_html)
