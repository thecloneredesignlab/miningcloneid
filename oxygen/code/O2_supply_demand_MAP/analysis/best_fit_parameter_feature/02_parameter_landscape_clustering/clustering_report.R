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

display_o2 <- function(x) {
  gsub("O2", "O\u2082", as.character(x), fixed = TRUE)
}

file_to_data_uri <- function(path, mime = "image/png") {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Required R package is not installed: base64enc")
  }
  if (!file.exists(path)) stop("Missing asset: ", path)
  paste0("data:", mime, ";base64,", base64enc::base64encode(path))
}

as_html_list <- function(items) {
  paste0("<ul>", paste0("<li>", html_escape(display_o2(items)), "</li>", collapse = ""), "</ul>")
}

report_tables_wclusters_dir <- function(root_dir, dataset, reduction = "umap") {
  reduction <- normalize_reduction(reduction)
  if (identical(dataset, "pooled_invivo_invitro")) {
    return(paper_pooled_reduction_tables_wclusters_dir(reduction = reduction, root_dir = root_dir))
  }
  paper_reduction_tables_wclusters_dir(dataset, reduction = reduction, root_dir = root_dir)
}

report_figures_wclusters_dir <- function(root_dir, dataset, reduction = "umap") {
  reduction <- normalize_reduction(reduction)
  if (identical(dataset, "pooled_invivo_invitro")) {
    return(paper_pooled_reduction_figures_wclusters_dir(reduction = reduction, root_dir = root_dir))
  }
  paper_reduction_figures_wclusters_dir(dataset, reduction = reduction, root_dir = root_dir)
}

report_figures_dir <- function(root_dir, dataset, reduction = "umap") {
  reduction <- normalize_reduction(reduction)
  if (identical(dataset, "pooled_invivo_invitro")) {
    return(paper_pooled_reduction_figures_dir(reduction = reduction, root_dir = root_dir))
  }
  paper_reduction_figures_dir(dataset, reduction = reduction, root_dir = root_dir)
}

report_output_subdirs <- function(dataset, prefix) {
  if (identical(dataset, "pooled_invivo_invitro")) {
    return(unique(c(pooled_prefix_output_subdir(prefix), "")))
  }
  subdir <- single_dataset_prefix_output_subdir(prefix)
  unique(c(if (nzchar(subdir)) subdir else character(), ""))
}

report_candidate_paths <- function(base_dir, dataset, prefix, ...) {
  components <- list(...)
  vapply(report_output_subdirs(dataset, prefix), function(subdir) {
    parts <- c(list(base_dir), if (nzchar(subdir)) list(subdir) else list(), components)
    do.call(file.path, parts)
  }, character(1), USE.NAMES = FALSE)
}

first_existing_path <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing)) existing[[1]] else paths[[1]]
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
  header <- paste0("<th>", html_escape(display_o2(names(vals))), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(vals)), function(i) {
    paste0("<tr>", paste0("<td>", html_escape(display_o2(unlist(vals[i, , drop = TRUE], use.names = FALSE))), "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0("<table class=\"report-table\"><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

section_title <- function(number, label) {
  number <- as.character(number %||% "")
  if (nzchar(number)) paste(number, label) else label
}

child_number <- function(parent_number, index) {
  parent_number <- as.character(parent_number %||% "")
  if (!nzchar(parent_number)) return(as.character(index))
  paste(parent_number, index, sep = ".")
}

spec_reduction <- function(spec) {
  normalize_reduction(spec$reduction %||% "umap")
}

spec_coordinate_source_dir <- function(spec) {
  spec$coordinate_source_dir %||% reduction_coordinate_cluster_dir(spec_reduction(spec))
}

spec_coordinate_source_label <- function(spec) {
  spec$coordinate_source_label %||% paste(reduction_axis_labels(spec_reduction(spec)), collapse = "/")
}

silhouette_path <- function(root_dir, dataset, source_dir, prefix, reduction = "umap") {
  first_existing_path(report_candidate_paths(
    report_tables_wclusters_dir(root_dir, dataset, reduction = reduction),
    dataset,
    prefix,
    source_dir,
    paste0(prefix, "_silhouette.csv")
  ))
}

read_silhouette_selected <- function(root_dir, dataset, source_dir, prefix, reduction = "umap") {
  path <- silhouette_path(root_dir, dataset, source_dir, prefix, reduction = reduction)
  if (!file.exists(path)) stop("Missing silhouette table: ", path)
  tab <- read_csv_plain(path)
  tab$selected <- coerce_logical_column(tab$selected, "selected")
  selected <- tab[tab$selected, , drop = FALSE]
  if (nrow(selected) != 1L) stop("Silhouette table must contain exactly one selected row: ", path)
  data.frame(
    `clustering basis` = if (identical(source_dir, "ByInputFeatures")) "input variables used for the embedding" else paste(reduction_axis_labels(reduction), collapse = "/"),
    `selected k` = selected$k,
    `mean silhouette width` = selected$average_silhouette,
    `silhouette sample size` = selected$sample_n,
    check.names = FALSE
  )
}

figure_paths <- function(root_dir, dataset, source_dir, prefix, reduction = "umap") {
  fig_dir <- report_figures_wclusters_dir(root_dir, dataset, reduction = reduction)
  list(
    png = first_existing_path(report_candidate_paths(fig_dir, dataset, prefix, source_dir, paste0(prefix, ".png")))
  )
}

base_figure_paths <- function(root_dir, dataset, prefix, reduction = "umap") {
  fig_dir <- report_figures_dir(root_dir, dataset, reduction = reduction)
  list(
    png = first_existing_path(report_candidate_paths(fig_dir, dataset, prefix, paste0(prefix, ".png")))
  )
}

cluster_output_pair_exists <- function(root_dir, dataset, spec) {
  reduction <- spec_reduction(spec)
  base_fig <- base_figure_paths(root_dir, dataset, spec$prefix, reduction = reduction)
  source_dirs <- c("ByInputFeatures", spec_coordinate_source_dir(spec))
  file.exists(base_fig$png) &&
    all(vapply(source_dirs, function(src) {
    fig <- figure_paths(root_dir, dataset, src, spec$prefix, reduction = reduction)
    file.exists(fig$png) &&
      file.exists(silhouette_path(root_dir, dataset, src, spec$prefix, reduction = reduction))
  }, logical(1)))
}

figure_card <- function(fig, label, legend) {
  img <- if (file.exists(fig$png)) {
    paste0("<img class=\"report-figure-image\" src=\"", file_to_data_uri(fig$png), "\" alt=\"", html_escape(label), "\">")
  } else {
    paste0("<p class=\"report-empty\">Missing image: ", html_escape(fig$png), "</p>")
  }
  paste0(
    "<article class=\"report-figure-card\">",
    "<div class=\"report-figure\">", img, "</div>",
    "<p class=\"report-figure-title\"><strong>", html_escape(label), "</strong></p>",
    "<p class=\"report-figure-legend\">", html_escape(legend), "</p>",
    "</article>"
  )
}

invivo_growth_turnover_o2_cluster_summary <- function(root_dir, prefix, source_dir) {
  coord_path <- first_existing_path(report_candidate_paths(
    paper_tables_wclusters_dir("invivo", root_dir = root_dir),
    "invivo",
    prefix,
    source_dir,
    paste0(prefix, "_coordinates.csv")
  ))
  metric_path <- file.path(paper_tables_dir("invivo", root_dir = root_dir), "invivo_best_params_growth_turnover_100d.csv")
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

extra_summary_block <- function(root_dir,
                                dataset,
                                spec,
                                section_number = NULL,
                                next_table_number = NULL,
                                add_table_entry = NULL) {
  summary_type <- spec$extra_summary %||% ""
  if (identical(summary_type, "invivo_growth_turnover_o2_by_input_cluster")) {
    input_tab <- invivo_growth_turnover_o2_cluster_summary(root_dir, spec$prefix, "ByInputFeatures")
    umap_tab <- invivo_growth_turnover_o2_cluster_summary(root_dir, spec$prefix, "ByUMAPCoordinates")
    next_num <- function() {
      if (is.function(next_table_number) && !is.null(section_number)) {
        next_table_number(section_number)
      } else {
        ""
      }
    }
    input_table_number <- next_num()
    umap_table_number <- next_num()
    input_table_id <- paste0(spec$id, "-input-cluster-summary")
    umap_table_id <- paste0(spec$id, "-coordinate-cluster-summary")
    if (is.function(add_table_entry)) {
      add_table_entry(input_table_id, paste0("Table ", input_table_number, ". ", spec$title, " input-cluster process summary"))
      add_table_entry(umap_table_id, paste0("Table ", umap_table_number, ". ", spec$title, " coordinate-cluster process summary"))
    }
    return(paste0(
      "<h4>Cluster-level summaries</h4>",
      "<p>Values are cluster means computed from the same seed-level growth, turnover, and O2 summaries used in this UMAP.</p>",
      "<div class=\"report-table-grid report-table-grid--2\">",
      "<div class=\"table-block\" id=\"", html_escape(input_table_id), "\" data-nav-target=\"", html_escape(input_table_id), "\">",
      "<h5>", html_escape(display_o2(paste0("Table ", input_table_number, ". Clusters estimated from the input variables."))), "</h5>",
      table_to_html(input_tab, max_rows = 10L),
      "</div><div class=\"table-block\" id=\"", html_escape(umap_table_id), "\" data-nav-target=\"", html_escape(umap_table_id), "\">",
      "<h5>", html_escape(display_o2(paste0("Table ", umap_table_number, ". Clusters estimated from the UMAP coordinates."))), "</h5>",
      table_to_html(umap_tab, max_rows = 10L),
      "</div></div>"
    ))
  }
  ""
}

pair_block <- function(root_dir,
                       dataset,
                       spec,
                       section_number,
                       figure_number,
                       heading_level = 5L,
                       silhouette_tables = NULL,
                       next_table_number = NULL,
                       add_table_entry = NULL) {
  reduction <- spec_reduction(spec)
  heading_level <- max(3L, min(6L, as.integer(heading_level)))
  heading_tag <- paste0("h", heading_level)
  source_dirs <- c("ByInputFeatures", spec_coordinate_source_dir(spec))
  coordinate_label <- spec_coordinate_source_label(spec)
  source_labels <- c(
    "no cluster overlay",
    "clusters estimated from the input variables",
    paste0("clusters estimated from the ", coordinate_label, " coordinates")
  )
  fig_labels <- paste0("Figure ", figure_number, LETTERS[1:3], ". ", spec$title, " - ", source_labels)
  base_fig <- base_figure_paths(root_dir, dataset, spec$prefix, reduction = reduction)
  figs <- c(
    list(base_fig),
    lapply(source_dirs, function(src) figure_paths(root_dir, dataset, src, spec$prefix, reduction = reduction))
  )
  input_silhouette <- read_silhouette_selected(root_dir, dataset, "ByInputFeatures", spec$prefix, reduction = reduction)
  coordinate_silhouette <- read_silhouette_selected(root_dir, dataset, spec_coordinate_source_dir(spec), spec$prefix, reduction = reduction)
  legend_base <- paste0(
    spec$legend,
    spec$legend_suffix %||% " For best-fit seed results, color represents the objective value, with lower values shown in blue and higher values shown in yellow. The outline color identifies cluster membership."
  )
  if (is.null(silhouette_tables)) {
    silhouette_tables <- list(
      input = list(
        id = paste0(spec$id, "-input-silhouette"),
        number = paste(section_number, "1", sep = "."),
        title = "Input-variable clustering silhouette selection"
      ),
      coordinate = list(
        id = paste0(spec$id, "-coordinate-silhouette"),
        number = paste(section_number, "2", sep = "."),
        title = paste0(reduction_label(reduction), "-coordinate clustering silhouette selection")
      )
    )
  }
  paste0(
    "<section class=\"report-subsection\" id=\"", html_escape(spec$id), "\" data-nav-target=\"", html_escape(spec$id), "\">",
    "<", heading_tag, ">", html_escape(display_o2(section_title(section_number, spec$title))), "</", heading_tag, ">",
    "<p>", html_escape(display_o2(spec$description)), "</p>",
    if (length(spec$bullets %||% character())) as_html_list(spec$bullets) else "",
    "<div class=\"report-table-grid report-table-grid--2\">",
    "<div class=\"table-block\" id=\"", html_escape(silhouette_tables$input$id), "\" data-nav-target=\"", html_escape(silhouette_tables$input$id), "\">",
    "<p class=\"table-caption\"><strong>",
    html_escape(display_o2(paste0("Table ", silhouette_tables$input$number, ". ", silhouette_tables$input$title, "."))),
    "</strong> Number of clusters selected by the silhouette criterion for clusters estimated from the input variables used to build this reduction.</p>",
    table_to_html(input_silhouette, max_rows = 10L),
    "</div>",
    "<div class=\"table-block\" id=\"", html_escape(silhouette_tables$coordinate$id), "\" data-nav-target=\"", html_escape(silhouette_tables$coordinate$id), "\">",
    "<p class=\"table-caption\"><strong>",
    html_escape(display_o2(paste0("Table ", silhouette_tables$coordinate$number, ". ", silhouette_tables$coordinate$title, "."))),
    "</strong> Number of clusters selected by the silhouette criterion for clusters estimated from the plotted ", html_escape(display_o2(coordinate_label)), " coordinate data.</p>",
    table_to_html(coordinate_silhouette, max_rows = 10L),
    "</div></div>",
    extra_summary_block(
      root_dir,
      dataset,
      spec,
      section_number = section_number,
      next_table_number = next_table_number,
      add_table_entry = add_table_entry
    ),
    "<div class=\"report-figure-grid report-figure-grid--3\">",
    figure_card(figs[[1]], fig_labels[[1]], paste0(legend_base, " This panel shows the original reduction without cluster labels or cluster contours.")),
    figure_card(figs[[2]], fig_labels[[2]], paste0(legend_base, " Clusters in this panel were estimated from the variables used to construct the embedding.")),
    figure_card(figs[[3]], fig_labels[[3]], paste0(legend_base, " Clusters in this panel were estimated from the two-dimensional ", coordinate_label, " coordinates.")),
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

pooled_parameter_table <- function() {
  params <- pooled_umap_parameter_set()
  log_params <- pooled_umap_log10_parameter_set()
  data.frame(
    parameter = params,
    log10_transformed = params %in% log_params,
    stringsAsFactors = FALSE
  )
}

preprocessing_guide_table <- function(dataset, reductions = c("pca", "umap", "tsne")) {
  reductions <- normalize_report_reductions(reductions)
  method_label <- reduction_slash_label(reductions)
  include_pca_umap <- "umap" %in% reductions
  if (identical(dataset, "pooled_invivo_invitro")) {
    rows <- data.frame(
      report_label = c(
        paste0("z-score shared-parameter ", method_label),
        paste0("context-prior-unit ", method_label),
        paste0("common-prior-unit ", method_label)
      ),
      data_rows = rep("Pooled in vivo and in vitro rows for the specific panel.", 3L),
      parameter_set = rep("The 14 parameters shared by in vivo and in vitro.", 3L),
      preprocessing = c(
        "Positive-scale parameters are log10-transformed where applicable, then each shared parameter is centered and scaled by its empirical mean and SD over the pooled rows in that panel.",
        "Each row is scaled within its own source context: in vivo rows use in vivo optimizer bounds; in vitro rows use in vitro optimizer bounds. Values are represented as prior-unit positions.",
        "A single common bound is built for each shared parameter using min(lower_in_vivo, lower_in_vitro) and max(upper_in_vivo, upper_in_vitro); both sources use that common prior-unit scale."
      ),
      interpretation = c(
        "Emphasizes empirical variation in the current pooled data; all parameters contribute with unit variance after z-scoring.",
        "Emphasizes where each source lies within its own optimizer search context; identical natural values can map differently if source-specific bounds differ.",
        "Emphasizes direct cross-source comparability on one shared optimizer scale; narrower source-specific ranges occupy only part of the common interval."
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (isTRUE(include_pca_umap)) {
      rows <- rbind(rows, data.frame(
        report_label = "PCA-UMAPs",
        data_rows = "Same rows and preprocessing as the matching UMAP variant.",
        parameter_set = "The same shared-parameter matrix after z-score, context-prior-unit, or common-prior-unit preprocessing.",
        preprocessing = "PCA is first run on the preprocessed matrix; retained PC scores are then embedded by UMAP. These are UMAP-coordinate plots, not true PCA-coordinate plots.",
        interpretation = "Tests whether the UMAP structure is stable after PCA compression of the same input representation.",
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
    }
    return(rbind(rows, data.frame(
      report_label = c("full data panels", "sampled panels"),
      data_rows = c(
        "All available DEoptim initial rows plus all seed-level best-fit rows.",
        "Down-sampled DEoptim initial rows plus all seed-level best-fit rows."
      ),
      parameter_set = rep("The same parameter set used by the corresponding reduction/preprocessing variant.", 2L),
      preprocessing = c(
        "No sampling is applied to initial-population rows.",
        "Sampling is applied only to initial-population rows; best-fit rows are retained."
      ),
      interpretation = c(
        "Shows the complete initial-vs-best landscape when full initial data were generated.",
        "Shows the same landscape with initial-row density reduced for visual balance and runtime."
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )))
  }

  label <- dataset_label(dataset)
  n_params <- if (identical(dataset, "invitro")) 14L else 18L
  data_desc <- paste0(label, " rows for the specific panel.")
  rows <- data.frame(
    report_label = c(
      paste0("z-score ", n_params, "-parameter ", method_label),
      paste0("prior-unit ", n_params, "-parameter ", method_label)
    ),
    data_rows = rep(data_desc, 2L),
    parameter_set = c(
      paste0("The base ", label, " parameter set used for this workflow."),
      paste0("The same ", label, " parameter set represented on optimizer prior-unit scales.")
    ),
    preprocessing = c(
      "Positive-scale parameters are log10-transformed where applicable, then each column is centered and scaled by its empirical mean and SD over the rows in that panel.",
      "Each parameter is mapped to its optimizer prior interval as (value - lower_bound) / (upper_bound - lower_bound)."
    ),
    interpretation = c(
      "Emphasizes empirical variation in the current data set; all parameters contribute with unit variance after z-scoring.",
      "Emphasizes where each seed lies within the optimizer search range for each parameter."
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (isTRUE(include_pca_umap)) {
    rows <- rbind(rows, data.frame(
      report_label = "PCA-UMAPs",
      data_rows = "Same rows and preprocessing as the matching UMAP variant.",
      parameter_set = "The same parameter matrix after z-score or prior-unit preprocessing.",
      preprocessing = "PCA is first run on the preprocessed matrix; retained PC scores are then embedded by UMAP. These are UMAP-coordinate plots, not true PCA-coordinate plots.",
      interpretation = "Tests whether UMAP structure is stable after PCA compression of the same input representation.",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  rbind(rows, data.frame(
    report_label = c("full data panels", "sampled panels", "best-fit only panels"),
    data_rows = c(
      "All available DEoptim initial rows plus all seed-level best-fit rows.",
      "Down-sampled DEoptim initial rows plus all seed-level best-fit rows.",
      "Seed-level best-fit rows only."
    ),
    parameter_set = c(
      "The same parameter set used by the corresponding reduction/preprocessing variant.",
      "The same parameter set used by the corresponding reduction/preprocessing variant.",
      paste0("The seed-level best-fit ", label, " parameter vectors and any explicitly added phenotype/process features for that panel.")
    ),
    preprocessing = c(
      "No sampling is applied to initial-population rows.",
      "Sampling is applied only to initial-population rows; best-fit rows are retained.",
      "Initial-population rows are not included."
    ),
    interpretation = c(
      "Shows the complete initial-vs-best landscape when full initial data were generated.",
      "Shows the same landscape with initial-row density reduced for visual balance and runtime.",
      "Focuses on optimized seed-level solutions without the surrounding DEoptim initial cloud."
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
}

best_parameter_reduction_specs <- function(dataset, legend) {
  dataset <- normalize_dataset(dataset)
  label <- dataset_label(dataset)
  n_params <- if (identical(dataset, "invitro")) 14L else 18L
  base <- paste0(dataset, "_best_params")
  list(
    list(
      id = paste0(dataset, "-best-prior-unit-umap"),
      prefix = paste0(base, "_prior_unit_umap"),
      title = paste0(label, " best-fit results: ", n_params, "-parameter prior-unit UMAP"),
      description = "This analysis repeats the best-fit parameter embedding after rescaling each transformed optimizer parameter to its prior interval [0, 1].",
      bullets = c(
        "Input rows are the seed-level best parameter vectors.",
        "Positive-scale parameters are transformed on the optimizer log10 scale before prior-unit rescaling.",
        "The original z-score UMAP is retained as the reference version."
      ),
      legend = legend,
      reduction = "umap"
    ),
    list(
      id = paste0(dataset, "-best-pca-umap"),
      prefix = paste0(base, "_pca_umap"),
      title = paste0(label, " best-fit results: PCA-to-UMAP sensitivity"),
      description = "This sensitivity analysis runs UMAP after retaining principal-component scores from the z-score input matrix.",
      bullets = c(
        "The output remains a UMAP coordinate plot.",
        "This panel is retained separately from the first-class PCA coordinate plots."
      ),
      legend = legend,
      reduction = "umap"
    ),
    list(
      id = paste0(dataset, "-best-prior-unit-pca-umap"),
      prefix = paste0(base, "_prior_unit_pca_umap"),
      title = paste0(label, " best-fit results: prior-unit PCA-to-UMAP sensitivity"),
      description = "This sensitivity analysis runs UMAP after PCA on prior-unit transformed parameters.",
      bullets = c(
        "The output remains a UMAP coordinate plot.",
        "PCA is centered for the prior-unit input before retained PCs are passed to UMAP."
      ),
      legend = legend,
      reduction = "umap"
    ),
    list(
      id = paste0(dataset, "-best-pca"),
      prefix = paste0(base, "_pca"),
      title = paste0(label, " best-fit results: ", n_params, "-parameter PCA"),
      description = "This analysis plots the first two principal-component coordinates for the best-fit parameter vectors.",
      bullets = c(
        "The x- and y-axes are PCA 1 and PCA 2.",
        "PCA variance tables are written next to the coordinate table."
      ),
      legend = legend,
      reduction = "pca"
    ),
    list(
      id = paste0(dataset, "-best-prior-unit-pca"),
      prefix = paste0(base, "_prior_unit_pca"),
      title = paste0(label, " best-fit results: prior-unit PCA"),
      description = "This analysis plots PCA coordinates after prior-unit transformed parameters are centered for PCA.",
      bullets = c(
        "The x- and y-axes are PCA 1 and PCA 2.",
        "Each parameter contributes on its optimizer prior scale rather than its empirical z-score scale."
      ),
      legend = legend,
      reduction = "pca"
    ),
    list(
      id = paste0(dataset, "-best-tsne"),
      prefix = paste0(base, "_tsne"),
      title = paste0(label, " best-fit results: ", n_params, "-parameter t-SNE"),
      description = "This analysis plots a two-dimensional t-SNE embedding of the z-score transformed best-fit parameter vectors.",
      bullets = c(
        "The x- and y-axes are t-SNE 1 and t-SNE 2.",
        "t-SNE is treated as a visualization companion to UMAP and PCA."
      ),
      legend = legend,
      reduction = "tsne"
    ),
    list(
      id = paste0(dataset, "-best-prior-unit-tsne"),
      prefix = paste0(base, "_prior_unit_tsne"),
      title = paste0(label, " best-fit results: prior-unit t-SNE"),
      description = "This analysis plots t-SNE after prior-unit transformed parameters.",
      bullets = c(
        "The x- and y-axes are t-SNE 1 and t-SNE 2.",
        "The embedding input is the transformed optimizer coordinate rescaled to [0, 1]."
      ),
      legend = legend,
      reduction = "tsne"
    )
  )
}

discover_pooled_sampled_prefixes <- function(root_dir) {
  tab_dir <- paper_pooled_tables_wclusters_dir(root_dir = root_dir)
  if (!dir.exists(tab_dir)) return(character(0))
  files <- list.files(
    tab_dir,
    pattern = "^pooled_invivo_invitro_initial_vs_best_sampled[0-9]+_umap_silhouette\\.csv$",
    recursive = TRUE
  )
  files <- files[basename(dirname(files)) == "ByInputFeatures"]
  sort(unique(sub("_silhouette\\.csv$", "", basename(files))))
}

sample_n_from_pooled_prefix <- function(prefix) {
  as.integer(sub("^.*_sampled([0-9]+)_umap$", "\\1", prefix))
}

discover_clustered_prefixes <- function(root_dir, dataset, reduction = "umap") {
  reduction <- normalize_reduction(reduction)
  tab_dir <- report_tables_wclusters_dir(root_dir, dataset, reduction = reduction)
  if (!dir.exists(tab_dir)) return(character(0))
  files <- list.files(tab_dir, pattern = "_silhouette\\.csv$", recursive = TRUE)
  files <- files[basename(dirname(files)) == "ByInputFeatures"]
  prefixes <- sort(unique(sub("_silhouette\\.csv$", "", basename(files))))
  prefixes[vapply(prefixes, function(prefix) {
    cluster_output_pair_exists(
      root_dir,
      dataset,
      list(prefix = prefix, reduction = reduction)
    )
  }, logical(1))]
}

scope_label <- function(dataset) {
  if (identical(dataset, "pooled_invivo_invitro")) return("pooled in vivo/in vitro")
  dataset_label(dataset)
}

reduction_label <- function(reduction) {
  switch(
    normalize_reduction(reduction),
    umap = "UMAP",
    pca = "PCA",
    tsne = "t-SNE"
  )
}

reduction_group_title <- function(reduction) {
  switch(
    normalize_reduction(reduction),
    umap = "UMAPs",
    pca = "PCAs",
    tsne = "t-SNEs"
  )
}

normalize_report_reductions <- function(reductions) {
  reductions <- as_char_vec(reductions, default = c("pca", "umap", "tsne"))
  reductions <- unique(vapply(reductions, normalize_reduction, character(1)))
  ordered <- c("pca", "umap", "tsne")
  reductions <- ordered[ordered %in% reductions]
  if (!length(reductions)) stop("No valid reductions were requested.")
  reductions
}

reduction_slash_label <- function(reductions) {
  paste(vapply(normalize_report_reductions(reductions), reduction_label, character(1)), collapse = "/")
}

reduction_sentence_label <- function(reductions) {
  labels <- vapply(normalize_report_reductions(reductions), reduction_label, character(1))
  if (length(labels) == 1L) return(labels[[1L]])
  if (length(labels) == 2L) return(paste(labels, collapse = " or "))
  paste0(paste(labels[-length(labels)], collapse = ", "), ", or ", labels[[length(labels)]])
}

cluster_report_filename <- function(reductions) {
  reductions <- normalize_report_reductions(reductions)
  if (identical(reductions, "pca")) return("parameter_landscape_clustering_pca_cluster_report.html")
  "parameter_landscape_clustering_umap_cluster_report.html"
}

preprocess_label_from_prefix <- function(prefix, pooled = FALSE) {
  if (grepl("common_prior_unit", prefix, fixed = TRUE)) return("common-prior-unit")
  if (grepl("context_prior_unit", prefix, fixed = TRUE)) return("context-prior-unit")
  if (grepl("prior_unit", prefix, fixed = TRUE)) return("prior-unit")
  if (isTRUE(pooled)) return("z-score shared-parameter")
  "z-score"
}

analysis_subject_from_prefix <- function(prefix, dataset) {
  if (identical(dataset, "pooled_invivo_invitro")) {
    if (grepl("sampled[0-9]+", prefix)) {
      sample_n <- suppressWarnings(as.integer(sub("^.*_sampled([0-9]+).*$", "\\1", prefix)))
      if (is.finite(sample_n)) return(paste0("sampled initial rows (sample_n = ", sample_n, ")"))
      return("sampled initial rows")
    }
    return("initial and best-fit results")
  }
  if (grepl("deoptim_initial_vs_best_sampled[0-9]+", prefix)) {
    sample_n <- suppressWarnings(as.integer(sub("^.*_sampled([0-9]+).*$", "\\1", prefix)))
    if (is.finite(sample_n)) return(paste0("DEoptim initial and best-fit results with sampled initial rows (sample_n = ", sample_n, ")"))
    return("DEoptim initial and best-fit results with sampled initial rows")
  }
  if (grepl("deoptim_initial_vs_best", prefix, fixed = TRUE)) return("DEoptim initial and best-fit results")
  if (grepl("best_params_growth_turnover_with_misseg_death_O2_attractor_dominant_ploidy", prefix, fixed = TRUE)) {
    return("best-fit results plus growth, turnover, O2, and fixed-O2 attractor ploidy")
  }
  if (grepl("best_params_growth_turnover_with_misseg_death_O2_ploidy_ratio", prefix, fixed = TRUE)) {
    return("best-fit results plus growth, turnover, O2, and ploidy ratios")
  }
  if (grepl("best_params_growth_turnover_with_misseg_death_O2", prefix, fixed = TRUE)) {
    return("best-fit results plus growth, turnover, and mean O2")
  }
  if (grepl("best_params_growth_turnover_with_misseg_death", prefix, fixed = TRUE)) {
    return("best-fit results plus net growth and turnover")
  }
  if (grepl("best_params_growth_turnover_O2_ploidy_ratio", prefix, fixed = TRUE)) {
    return("best-fit results plus growth, turnover, O2, and ploidy ratios")
  }
  if (grepl("best_params_growth_turnover_O2", prefix, fixed = TRUE)) {
    return("best-fit results plus growth, turnover, and mean O2")
  }
  if (grepl("best_params_growth_turnover", prefix, fixed = TRUE)) {
    return("best-fit results plus net growth and turnover")
  }
  "best-fit results"
}

generic_reduction_spec <- function(dataset, prefix, reduction, legend, legend_suffix = NULL) {
  reduction <- normalize_reduction(reduction)
  pooled <- identical(dataset, "pooled_invivo_invitro")
  method <- reduction_label(reduction)
  preprocess <- preprocess_label_from_prefix(prefix, pooled = pooled)
  subject <- analysis_subject_from_prefix(prefix, dataset)
  sensitivity <- identical(reduction, "umap") && grepl("_pca_umap$", prefix)
  title_method <- if (sensitivity) "PCA-to-UMAP sensitivity" else method
  axis_note <- switch(
    reduction,
    umap = if (sensitivity) {
      "The output remains a UMAP coordinate plot after PCA preprocessing."
    } else {
      "The x- and y-axes are UMAP 1 and UMAP 2."
    },
    pca = "The x- and y-axes are PCA 1 and PCA 2; PCA variance tables are written next to the coordinate table.",
    tsne = "The x- and y-axes are t-SNE 1 and t-SNE 2."
  )
  list(
    id = html_id(paste(dataset, prefix, reduction, sep = "-")),
    prefix = prefix,
    title = paste(scope_label(dataset), subject, "-", preprocess, title_method),
    description = paste0(
      "This panel is generated from the existing ", method,
      " output files for ", scope_label(dataset), " and is included so the report covers all current reduction methods."
    ),
    bullets = c(
      axis_note,
      "Cluster summaries compare k-means clusters estimated from the embedding-input variables against clusters estimated from the two-dimensional reduction coordinates.",
      "This entry was discovered from the latest result directory and included when both clustered figures and silhouette tables were present."
    ),
    legend = legend,
    legend_suffix = legend_suffix,
    reduction = reduction
  )
}

dedupe_specs <- function(specs) {
  if (!length(specs)) return(specs)
  keys <- vapply(specs, function(spec) paste(spec_reduction(spec), spec$prefix, sep = "::"), character(1))
  specs[!duplicated(keys)]
}

augment_specs_from_outputs <- function(root_dir, dataset, reduction, specs, legend, legend_suffix = NULL) {
  reduction <- normalize_reduction(reduction)
  specs <- specs[vapply(specs, function(spec) identical(spec_reduction(spec), reduction), logical(1))]
  existing_prefixes <- vapply(specs, function(spec) spec$prefix, character(1))
  discovered <- discover_clustered_prefixes(root_dir, dataset, reduction = reduction)
  missing <- setdiff(discovered, existing_prefixes)
  specs <- c(
    specs,
    lapply(missing, function(prefix) {
      generic_reduction_spec(
        dataset = dataset,
        prefix = prefix,
        reduction = reduction,
        legend = legend,
        legend_suffix = legend_suffix
      )
    })
  )
  dedupe_specs(available_specs(root_dir, dataset, specs))
}

specs_by_reduction <- function(root_dir, dataset, specs, legend, legend_suffix = NULL) {
  stats::setNames(
    lapply(c("umap", "pca", "tsne"), function(reduction) {
      augment_specs_from_outputs(
        root_dir = root_dir,
        dataset = dataset,
        reduction = reduction,
        specs = specs,
        legend = legend,
        legend_suffix = legend_suffix
      )
    }),
    c("umap", "pca", "tsne")
  )
}

available_specs <- function(root_dir, dataset, specs) {
  keep <- vapply(specs, function(spec) cluster_output_pair_exists(root_dir, dataset, spec), logical(1))
  specs[keep]
}

nav_link <- function(id, label, active = TRUE) {
  paste0(
    "<li><a href=\"#", html_escape(id), "\"",
    if (isTRUE(active)) {
      paste0(" data-target=\"", html_escape(id), "\"")
    } else {
      paste0(" data-quick-target=\"", html_escape(id), "\"")
    },
    ">",
    html_escape(display_o2(label)),
    "</a></li>"
  )
}

nav_entries_to_html <- function(entries, empty_text, active = TRUE) {
  if (!length(entries)) {
    return(paste0("<p class=\"nav-empty\">", html_escape(empty_text), "</p>"))
  }
  paste0(
    "<ul class=\"nav-list nav-sublist\">",
    paste0(vapply(entries, function(entry) nav_link(entry$id, entry$label, active = active), character(1)), collapse = ""),
    "</ul>"
  )
}

nav_group <- function(title, entries_html, count = NA_integer_, open = FALSE) {
  label <- if (is.finite(suppressWarnings(as.numeric(count)))) {
    paste0(title, " (", count, ")")
  } else {
    title
  }
  paste0(
    "<details class=\"nav-group\"", if (isTRUE(open)) " open" else "", ">",
    "<summary>", html_escape(display_o2(label)), "</summary>",
    entries_html,
    "</details>"
  )
}

filter_nav_entries <- function(entries, field, value) {
  entries[vapply(entries, function(entry) identical(entry[[field]], value), logical(1))]
}

toc_ref_list <- function(refs) {
  if (!length(refs)) return("<p class=\"nav-empty\">No table or figure references in this subsection.</p>")
  paste0(
    "<ul class=\"nav-list nav-tertiary\">",
    paste0(vapply(refs, function(ref) {
      nav_link(ref$id, paste0(if (identical(ref$type, "table")) "Table: " else "Figure: ", ref$label), active = TRUE)
    }, character(1)), collapse = ""),
    "</ul>"
  )
}

toc_heading_node <- function(heading, toc_refs) {
  refs <- filter_nav_entries(toc_refs, "heading_id", heading$id)
  paste0(
    "<details class=\"nav-branch nav-heading\">",
    "<summary><a href=\"#", html_escape(heading$id), "\" data-target=\"", html_escape(heading$id), "\">",
    html_escape(display_o2(heading$label)),
    "</a></summary>",
    toc_ref_list(refs),
    "</details>"
  )
}

toc_section_node <- function(section_id, section_label, toc_headings, toc_refs, child_html = "", open = FALSE, depth = 1L) {
  headings <- filter_nav_entries(toc_headings, "section_id", section_id)
  heading_html <- if (length(headings)) {
    paste0(vapply(headings, function(heading) toc_heading_node(heading, toc_refs), character(1)), collapse = "")
  } else {
    ""
  }
  body_html <- paste0(heading_html, child_html)
  if (!nzchar(body_html)) body_html <- "<p class=\"nav-empty\">No subsections available.</p>"
  depth <- max(1L, min(8L, as.integer(depth)))
  paste0(
    "<details class=\"nav-branch nav-section nav-depth-", depth, "\"", if (isTRUE(open)) " open" else "", ">",
    "<summary><a href=\"#", html_escape(section_id), "\" data-target=\"", html_escape(section_id), "\">",
    html_escape(display_o2(section_label)),
    "</a></summary>",
    body_html,
    "</details>"
  )
}

build_toc_html <- function(section_items, toc_headings, toc_refs) {
  if (!is.data.frame(section_items) || !nrow(section_items)) {
    return("<p class=\"nav-empty\">No sections available.</p>")
  }
  if (!"parent_id" %in% names(section_items)) section_items$parent_id <- ""
  section_items$parent_id[is.na(section_items$parent_id)] <- ""
  render_item <- function(i, open = FALSE, depth = 1L) {
    section_id <- section_items$id[[i]]
    child_idx <- which(section_items$parent_id == section_id)
    child_html <- if (length(child_idx)) {
      paste0(vapply(child_idx, function(j) render_item(j, open = FALSE, depth = depth + 1L), character(1)), collapse = "")
    } else {
      ""
    }
    toc_section_node(
      section_items$id[[i]],
      section_items$label[[i]],
      toc_headings,
      toc_refs,
      child_html = child_html,
      open = open,
      depth = depth
    )
  }
  top_idx <- which(!nzchar(section_items$parent_id))
  paste0(
    "<div class=\"toc-tree\">",
    paste0(vapply(seq_along(top_idx), function(k) render_item(top_idx[[k]], open = k == 1L, depth = 1L), character(1)), collapse = ""),
    "</div>"
  )
}

build_sidebar_nav <- function(section_items, toc_headings, toc_refs, table_entries, figure_entries) {
  toc_html <- build_toc_html(section_items, toc_headings, toc_refs)
  top_count <- if ("parent_id" %in% names(section_items)) {
    sum(!nzchar(ifelse(is.na(section_items$parent_id), "", section_items$parent_id)))
  } else {
    nrow(section_items)
  }
  paste0(
    nav_group("Contents", toc_html, count = top_count, open = TRUE),
    nav_group("Tables", nav_entries_to_html(table_entries, "No tables available.", active = FALSE), count = length(table_entries), open = FALSE),
    nav_group("Figures", nav_entries_to_html(figure_entries, "No figures available.", active = FALSE), count = length(figure_entries), open = FALSE)
  )
}

report_nav_script <- function() {
  paste0(
    "<script>",
    "(function(){",
    "var links=Array.prototype.slice.call(document.querySelectorAll('.sidebar a[data-target]'));",
    "var linkById={};links.forEach(function(a){linkById[a.getAttribute('data-target')]=a;});",
    "var targets=Array.prototype.slice.call(document.querySelectorAll('[data-nav-target]')).filter(function(el){return !!linkById[el.id];});",
    "function activeDetails(active){var keep=[];var group=active?active.closest('details'):null;while(group){keep.push(group);group=group.parentElement?group.parentElement.closest('details'):null;}return keep;}",
    "function collapseInactive(active){var keep=activeDetails(active);document.querySelectorAll('.toc-tree details.nav-branch').forEach(function(group){if(keep.indexOf(group)===-1){group.open=false;}});keep.forEach(function(group){group.open=true;});}",
    "function setActive(id){",
    "links.forEach(function(a){a.classList.toggle('active',a.getAttribute('data-target')===id);});",
    "var active=linkById[id];",
    "if(active){collapseInactive(active);var item=active.closest('li,details');if(item&&item.scrollIntoView){item.scrollIntoView({block:'nearest'});}}",
    "}",
    "links.forEach(function(a){a.addEventListener('click',function(){setActive(a.getAttribute('data-target'));});});",
    "if('IntersectionObserver' in window){",
    "var visible={};",
    "var observer=new IntersectionObserver(function(entries){",
    "entries.forEach(function(entry){if(entry.isIntersecting){visible[entry.target.id]=entry.boundingClientRect.top;}else{delete visible[entry.target.id];}});",
    "var best=null;Object.keys(visible).forEach(function(id){if(best===null||Math.abs(visible[id])<Math.abs(visible[best])){best=id;}});",
    "if(best!==null){setActive(best);}",
    "},{rootMargin:'-18% 0px -68% 0px',threshold:[0,0.05,0.25]});",
    "targets.forEach(function(el){observer.observe(el);});",
    "}else if(targets.length){setActive(targets[0].id);}",
    "})();",
    "</script>"
  )
}

build_report_html <- function(root_dir,
                              reductions = c("pca", "umap", "tsne"),
                              report_title = NULL) {
  reductions <- normalize_report_reductions(reductions)
  reduction_phrase <- reduction_sentence_label(reductions)
  report_title <- report_title %||% if (identical(reductions, "pca")) {
    "Parameter Landscape PCA Report"
  } else {
    "Parameter Landscape Reduction Report"
  }
  side_title <- if (identical(reductions, "pca")) "PCA Report" else "Reduction Report"
  invivo_mode_legend <- "Each point represents one seed-level best-fit result. Point shape indicates the fixed-O2 attractor mode at O2 = 2: circles denote Mode 1 and triangles denote Mode 2."
  invivo_specs_all <- list(
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
      legend = invivo_mode_legend
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
      legend = invivo_mode_legend
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
      legend = invivo_mode_legend
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
      legend = invivo_mode_legend
    ),
    list(
      id = "invivo-misseg-growth-turnover-o2-attractor-dominant-ploidy",
      prefix = "invivo_best_params_growth_turnover_with_misseg_death_O2_attractor_dominant_ploidy_100d_umap",
      title = "in vivo best-fit results: 18 parameters plus growth, turnover, O2, and fixed-O2 attractor ploidy",
      description = "This analysis replaces the two day-1000 ploidy-ratio variables with fixed-O2 attractor dominant-ploidy values evaluated over the selected O2 concentrations.",
      bullets = c(
        "The growth, turnover, and mean-O2 variables are the same variables used in the preceding O2-augmented representation.",
        "The additional attractor variables are the dominant mean ploidy values from the fixed-O2 analytical solution.",
        "These attractor variables are included as continuous variables and are standardized together with the other UMAP inputs."
      ),
      legend = invivo_mode_legend
    )
  )
  invivo_specs_all <- c(invivo_specs_all, best_parameter_reduction_specs("invivo", invivo_mode_legend))
  invivo_specs <- specs_by_reduction(
    root_dir = root_dir,
    dataset = "invivo",
    specs = invivo_specs_all,
    legend = invivo_mode_legend
  )

  invitro_legend <- "Each point represents one seed-level best-fit result. Point color encodes objective."
  invitro_specs_all <- c(list(
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
      legend = invitro_legend
    )
  ), best_parameter_reduction_specs("invitro", invitro_legend))
  invitro_specs <- specs_by_reduction(
    root_dir = root_dir,
    dataset = "invitro",
    specs = invitro_specs_all,
    legend = invitro_legend
  )

  pooled_legend_suffix <- paste0(
    " Initial points are gray. Best-fit in vivo points use a blue-to-yellow objective scale, ",
    "and best-fit in vitro points use a green-to-red objective scale; each color scale uses the original objective values for that data set. ",
    "Circles denote in vivo rows and triangles denote in vitro rows. Cluster membership is shown by the cluster contour and label annotations."
  )
  pooled_sampled_prefixes <- discover_pooled_sampled_prefixes(root_dir)
  pooled_specs_all <- c(
    list(
      list(
        id = "pooled-full",
        prefix = "pooled_invivo_invitro_initial_vs_best_umap",
        title = "pooled in vivo/in vitro initial and best-fit results",
        description = "This analysis embeds the pooled in vivo and in vitro initial-population rows together with the seed-level best-fit parameter vectors.",
        bullets = c(
          "Input rows combine in vivo and in vitro sources but keep source labels, point type, seed, and objective metadata.",
          "The parameter-table initial row is removed from each seed in both the in vivo and in vitro initial populations by default before pooling.",
          "The embedding uses the 14 parameters shared by the in vivo and in vitro UMAP representations."
        ),
        legend = "Each point is one pooled UMAP row.",
        legend_suffix = pooled_legend_suffix
      )
    ),
    lapply(pooled_sampled_prefixes, function(prefix) {
      sample_n <- sample_n_from_pooled_prefix(prefix)
      list(
        id = html_id(paste0("pooled-", prefix)),
        prefix = prefix,
        title = paste0("pooled in vivo/in vitro with sampled initial rows", if (!is.na(sample_n)) paste0(" (sample_n = ", sample_n, ")") else ""),
        description = "This analysis repeats the pooled UMAP after down-sampling initial-population rows separately within in vivo and in vitro, while retaining all seed-level best-fit parameter vectors.",
        bullets = c(
          "Down-sampling is applied to initial rows only and is stratified by source data set.",
          "Best-fit rows from both data sets are retained and colored by their original objective values.",
          "The same 14 shared parameters and transformation rules are used as in the full pooled UMAP."
        ),
        legend = "Each point is one pooled sampled UMAP row.",
        legend_suffix = pooled_legend_suffix
      )
    })
  )
  pooled_specs <- specs_by_reduction(
    root_dir = root_dir,
    dataset = "pooled_invivo_invitro",
    specs = pooled_specs_all,
    legend = "Each point is one pooled reduction row.",
    legend_suffix = pooled_legend_suffix
  )

  section_items <- list()
  table_entries <- list()
  figure_entries <- list()
  table_section_counts <- new.env(parent = emptyenv())
  figure_section_counts <- new.env(parent = emptyenv())

  section_title <- function(number, label) {
    number <- as.character(number %||% "")
    if (nzchar(number)) paste(number, label) else label
  }

  child_number <- function(parent_number, index) {
    parent_number <- as.character(parent_number %||% "")
    if (!nzchar(parent_number)) return(as.character(index))
    paste(parent_number, index, sep = ".")
  }

  next_section_item_number <- function(section_number, counter_env) {
    section_number <- as.character(section_number %||% "")
    if (!nzchar(section_number)) section_number <- "0"
    current <- if (exists(section_number, envir = counter_env, inherits = FALSE)) {
      get(section_number, envir = counter_env, inherits = FALSE)
    } else {
      0L
    }
    current <- current + 1L
    assign(section_number, current, envir = counter_env)
    paste(section_number, current, sep = ".")
  }

  next_table_number <- function(section_number) {
    next_section_item_number(section_number, table_section_counts)
  }

  next_figure_number <- function(section_number) {
    next_section_item_number(section_number, figure_section_counts)
  }

  add_section <- function(id, label, parent_id = "", number = "") {
    existing <- vapply(section_items, function(x) x$id, character(1))
    if (id %in% existing) return(invisible(FALSE))
    section_items[[length(section_items) + 1L]] <<- list(
      id = id,
      label = section_title(number, label),
      number = as.character(number %||% ""),
      parent_id = parent_id
    )
    invisible(TRUE)
  }

  add_table_entry <- function(id, label) {
    table_entries[[length(table_entries) + 1L]] <<- list(id = id, label = label)
    invisible(TRUE)
  }

  add_figure_entry <- function(id, label) {
    figure_entries[[length(figure_entries) + 1L]] <<- list(id = id, label = label)
    invisible(TRUE)
  }

  reduction_order <- reductions
  variant_order <- c("full", "sampled", "best")
  umap_subtype_order <- c("umap", "pca_umap")

  variant_label <- function(variant) {
    switch(
      variant,
      full = "full data",
      sampled = "sampled",
      best = "best-fit only",
      variant
    )
  }

  variant_note <- function(variant) {
    switch(
      variant,
      full = "These panels use the full initial-population rows together with seed-level best-fit rows.",
      sampled = "These panels use down-sampled initial-population rows together with all seed-level best-fit rows.",
      best = "These panels contain seed-level best-fit rows only; no DEoptim initial-population rows are included.",
      ""
    )
  }

  spec_variant <- function(spec, dataset) {
    prefix <- spec$prefix %||% ""
    if (grepl("sampled[0-9]+", prefix)) return("sampled")
    if (identical(dataset, "pooled_invivo_invitro")) return("full")
    if (grepl("deoptim_initial_vs_best", prefix, fixed = TRUE)) return("full")
    "best"
  }

  umap_subtype <- function(spec) {
    if (grepl("_pca_umap$", spec$prefix %||% "")) "pca_umap" else "umap"
  }

  umap_subtype_label <- function(subtype) {
    switch(
      subtype,
      umap = "UMAPs",
      pca_umap = "PCA-UMAPs",
      subtype
    )
  }

  split_specs_by_variant <- function(specs, dataset) {
    out <- stats::setNames(vector("list", length(variant_order)), variant_order)
    for (spec in specs) {
      key <- spec_variant(spec, dataset)
      out[[key]] <- c(out[[key]], list(spec))
    }
    out
  }

  split_umap_specs <- function(specs) {
    out <- stats::setNames(vector("list", length(umap_subtype_order)), umap_subtype_order)
    for (spec in specs) {
      key <- umap_subtype(spec)
      out[[key]] <- c(out[[key]], list(spec))
    }
    out
  }

  render_parameter_block <- function(scope_id, title, parameter_text, parameter_df, section_number) {
    param_id <- paste0(scope_id, "-parameters")
    table_number <- next_table_number(section_number)
    add_table_entry(param_id, paste0("Table ", table_number, ". ", title, " reduction input parameters"))
    paste0(
      "<div class=\"report-parameter-block\" id=\"", html_escape(param_id), "\" data-nav-target=\"", html_escape(param_id), "\">",
      "<p class=\"table-caption\"><strong>", html_escape(display_o2(paste0("Table ", table_number, ". Reduction input parameters."))), "</strong> ",
      html_escape(display_o2(parameter_text)), "</p>",
      table_to_html(parameter_df, max_rows = 50L),
      "</div>"
    )
  }

  render_preprocessing_guide_block <- function(scope_id, dataset, title, section_number) {
    guide_id <- paste0(scope_id, "-preprocessing-guide")
    table_number <- next_table_number(section_number)
    add_table_entry(guide_id, paste0("Table ", table_number, ". ", title, " preprocessing and data-use guide"))
    paste0(
      "<div class=\"report-parameter-block\" id=\"", html_escape(guide_id), "\" data-nav-target=\"", html_escape(guide_id), "\">",
      "<p class=\"table-caption\"><strong>", html_escape(display_o2(paste0("Table ", table_number, ". Preprocessing and data-use guide."))), "</strong> ",
      "This table defines how the report variants differ before ", html_escape(reduction_phrase), " is calculated.</p>",
      table_to_html(preprocessing_guide_table(dataset, reductions = reduction_order), max_rows = 20L),
      "</div>"
    )
  }

  render_spec_block <- function(dataset, spec, parent_id, parent_number, spec_index, heading_level) {
    section_number <- child_number(parent_number, spec_index)
    figure_number <- next_figure_number(section_number)
    input_table_number <- next_table_number(section_number)
    coordinate_table_number <- next_table_number(section_number)
    input_table_id <- paste0(spec$id, "-input-silhouette")
    coordinate_table_id <- paste0(spec$id, "-coordinate-silhouette")
    coordinate_title <- paste0(reduction_label(spec_reduction(spec)), "-coordinate clustering silhouette selection")
    add_section(spec$id, spec$title, parent_id, number = section_number)
    add_table_entry(input_table_id, paste0("Table ", input_table_number, ". ", spec$title, " input-variable clustering silhouette selection"))
    add_table_entry(coordinate_table_id, paste0("Table ", coordinate_table_number, ". ", spec$title, " ", coordinate_title))
    add_figure_entry(spec$id, paste0("Figure ", figure_number, ". ", spec$title))
    pair_block(
      root_dir,
      dataset,
      spec,
      section_number = section_number,
      figure_number = figure_number,
      heading_level = heading_level,
      silhouette_tables = list(
        input = list(
          id = input_table_id,
          number = input_table_number,
          title = "Input-variable clustering silhouette selection"
        ),
        coordinate = list(
          id = coordinate_table_id,
          number = coordinate_table_number,
          title = coordinate_title
        )
      ),
      next_table_number = next_table_number,
      add_table_entry = add_table_entry
    )
  }

  render_variant_group <- function(parent_id,
                                   parent_number,
                                   dataset,
                                   variant,
                                   variant_index,
                                   specs,
                                   heading_level = 4L,
                                   pair_heading_level = 5L,
                                   force = FALSE) {
    if (!length(specs) && !isTRUE(force)) return("")
    variant_id <- paste0(parent_id, "-", gsub("_", "-", variant, fixed = TRUE))
    variant_number <- child_number(parent_number, variant_index)
    heading_tag <- paste0("h", heading_level)
    add_section(variant_id, variant_label(variant), parent_id, number = variant_number)
    body <- if (length(specs)) {
      paste0(vapply(seq_along(specs), function(i) {
        render_spec_block(dataset, specs[[i]], variant_id, parent_number = variant_number, spec_index = i, heading_level = pair_heading_level)
      }, character(1)), collapse = "")
    } else {
      "<p class=\"report-empty\">No clustered outputs were found for this data subset.</p>"
    }
    paste0(
      "<section class=\"report-variant-section\" id=\"", html_escape(variant_id), "\" data-nav-target=\"", html_escape(variant_id), "\">",
      "<", heading_tag, ">", html_escape(display_o2(section_title(variant_number, variant_label(variant)))), "</", heading_tag, ">",
      "<p>", html_escape(display_o2(variant_note(variant))), "</p>",
      body,
      "</section>"
    )
  }

  render_standard_variants <- function(parent_id,
                                       parent_number,
                                       dataset,
                                       specs,
                                       heading_level = 4L,
                                       pair_heading_level = 5L) {
    grouped <- split_specs_by_variant(specs, dataset)
    paste0(
      render_variant_group(parent_id, parent_number, dataset, "full", 1L, grouped$full, heading_level, pair_heading_level, force = TRUE),
      render_variant_group(parent_id, parent_number, dataset, "sampled", 2L, grouped$sampled, heading_level, pair_heading_level, force = TRUE),
      render_variant_group(parent_id, parent_number, dataset, "best", 3L, grouped$best, heading_level, pair_heading_level, force = length(grouped$best) > 0L)
    )
  }

  render_umap_subtype_group <- function(scope_id, parent_number, dataset, subtype, subtype_index, specs) {
    subtype_id <- paste0(scope_id, "-umaps-", gsub("_", "-", subtype, fixed = TRUE))
    subtype_number <- child_number(parent_number, subtype_index)
    add_section(subtype_id, umap_subtype_label(subtype), paste0(scope_id, "-umaps"), number = subtype_number)
    paste0(
      "<section class=\"report-subtype-section\" id=\"", html_escape(subtype_id), "\" data-nav-target=\"", html_escape(subtype_id), "\">",
      "<h4>", html_escape(display_o2(section_title(subtype_number, umap_subtype_label(subtype)))), "</h4>",
      render_standard_variants(subtype_id, subtype_number, dataset, specs, heading_level = 5L, pair_heading_level = 6L),
      "</section>"
    )
  }

  render_reduction_group <- function(scope_id, scope_number, dataset, reduction, reduction_index, specs) {
    group_id <- paste0(scope_id, "-", reduction, "s")
    title <- reduction_group_title(reduction)
    group_number <- child_number(scope_number, reduction_index)
    add_section(group_id, title, scope_id, number = group_number)
    body <- if (identical(reduction, "umap")) {
      grouped <- split_umap_specs(specs)
      paste0(vapply(seq_along(umap_subtype_order), function(i) {
        subtype <- umap_subtype_order[[i]]
        render_umap_subtype_group(scope_id, group_number, dataset, subtype, i, grouped[[subtype]])
      }, character(1)), collapse = "")
    } else {
      render_standard_variants(group_id, group_number, dataset, specs, heading_level = 4L, pair_heading_level = 5L)
    }
    paste0(
      "<section class=\"report-reduction-section\" id=\"", html_escape(group_id), "\" data-nav-target=\"", html_escape(group_id), "\">",
      "<h3>", html_escape(display_o2(section_title(group_number, title))), "</h3>",
      body,
      "</section>"
    )
  }

  render_scope_section <- function(scope_id,
                                   dataset,
                                   title,
                                   parameter_heading,
                                   parameter_text,
                                   parameter_df,
                                   specs,
                                   scope_number) {
    add_section(scope_id, title, "", number = scope_number)
    paste0(
      "<section class=\"report-section\" id=\"", html_escape(scope_id), "\" data-nav-target=\"", html_escape(scope_id), "\"><h2>", html_escape(display_o2(section_title(scope_number, title))), "</h2>",
      render_parameter_block(scope_id, title, parameter_text, parameter_df, section_number = scope_number),
      render_preprocessing_guide_block(scope_id, dataset, title, section_number = scope_number),
      paste0(vapply(seq_along(reduction_order), function(i) {
        reduction <- reduction_order[[i]]
        render_reduction_group(scope_id, scope_number, dataset, reduction, i, specs[[reduction]])
      }, character(1)), collapse = ""),
      "</section>"
    )
	  }

  report_body <- paste0(
    render_scope_section(
      scope_id = "pooled",
      dataset = "pooled_invivo_invitro",
      title = "pooled in vivo/in vitro",
      parameter_heading = "Pooled reduction input parameters",
      parameter_text = "The pooled reductions use the intersection of the in vivo and in vitro parameter sets.",
      parameter_df = pooled_parameter_table(),
      specs = pooled_specs,
      scope_number = "1"
    ),
    render_scope_section(
      scope_id = "invivo",
      dataset = "invivo",
      title = "in vivo",
      parameter_heading = "in vivo reduction input parameters",
      parameter_text = "The base in vivo reductions use 18 multi-paired soft-coupling parameters. The table records whether each parameter is log10-transformed before standardization.",
      parameter_df = parameter_table("invivo"),
      specs = invivo_specs,
      scope_number = "2"
    ),
    render_scope_section(
      scope_id = "invitro",
      dataset = "invitro",
      title = "in vitro",
      parameter_heading = "in vitro reduction input parameters",
      parameter_text = "The in vitro reductions use 14 soft-coupling parameters and do not include any ploidy trajectory-derived features.",
      parameter_df = parameter_table("invitro"),
      specs = invitro_specs,
      scope_number = "3"
    )
  )

  section_df <- if (length(section_items)) {
    do.call(rbind, lapply(section_items, function(x) {
      data.frame(
        id = x$id,
        label = x$label,
        parent_id = x$parent_id %||% "",
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(id = character(), label = character(), parent_id = character(), stringsAsFactors = FALSE)
  }
  nav <- build_sidebar_nav(section_df, list(), list(), table_entries, figure_entries)

  css <- paste0(
    "html{scroll-behavior:smooth;}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f4f7fa;color:#1b2a38;}",
    ".shell{display:flex;gap:26px;max-width:1760px;margin:0 auto;padding:22px;}.sidebar{position:sticky;top:22px;align-self:flex-start;width:320px;max-height:calc(100vh - 44px);overflow:auto;border:1px solid #d8e0e8;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.07);scrollbar-gutter:stable;}",
    ".side-head{padding:16px;background:#26394d;color:#fff;}.kicker{font-size:11px;letter-spacing:.08em;text-transform:uppercase;opacity:.78;font-weight:700;}.side-title{font-size:18px;font-weight:700;line-height:1.2;margin-top:4px;}.side-subtitle{font-size:12px;line-height:1.35;margin-top:4px;opacity:.84;}",
    ".nav-group{border-bottom:1px solid #e4ebf2;}.nav-group summary{cursor:pointer;list-style:none;padding:10px 14px;font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#334e68;background:#f7fafc;}.nav-group summary::-webkit-details-marker,.nav-branch summary::-webkit-details-marker{display:none;}.nav-group summary:before,.nav-branch summary:before{content:'+';display:inline-block;width:16px;color:#6b7c8f;}.nav-group[open]>summary:before,.nav-branch[open]>summary:before{content:'-';}.toc-tree{padding:8px 8px 12px;}.nav-branch{margin:3px 0;border-radius:7px;}.nav-depth-1{margin-left:0;}.nav-depth-2{margin-left:8px;}.nav-depth-3{margin-left:16px;}.nav-depth-4{margin-left:24px;}.nav-depth-5{margin-left:32px;}.nav-depth-6,.nav-depth-7,.nav-depth-8{margin-left:40px;}.nav-branch summary{cursor:pointer;list-style:none;padding:7px 8px;font-size:12px;line-height:1.3;color:#21384d;background:#fff;border-radius:7px;}.nav-depth-1>summary{font-weight:700;}.nav-depth-2>summary{font-weight:650;}.nav-depth-3>summary{font-size:11.5px;}.nav-depth-4>summary,.nav-depth-5>summary,.nav-depth-6>summary{font-size:11px;color:#334e68;}.nav-branch summary:hover{background:#eef4fa;}.nav-branch summary a{text-decoration:none;color:inherit;}.nav-section{border:1px solid #edf2f7;background:#fbfdff;}.nav-heading{margin-left:10px;background:#fff;}.nav-heading summary{font-size:11px;color:#334e68;}.nav-list{list-style:none;margin:0;padding:8px 10px 10px;}.nav-sublist{max-height:52vh;overflow:auto;}.nav-tertiary{padding:4px 4px 8px 22px;}.nav-list li{margin:2px 0;}.nav-list a{display:block;padding:7px 10px;border-radius:7px;text-decoration:none;color:#21384d;font-size:12px;line-height:1.3;overflow:hidden;text-overflow:ellipsis;}.nav-tertiary a{font-size:11px;padding:5px 8px;color:#334e68;}.nav-list a:hover{background:#eef4fa;}.nav-list a.active,.nav-branch summary a.active{background:#dcecf8;color:#102a43;font-weight:700;box-shadow:inset 3px 0 0 #2f80c0;}.nav-branch summary a.active{display:inline-block;border-radius:6px;padding:3px 6px;}.nav-empty{font-size:12px;color:#6b7c8f;margin:8px 14px 12px;}",
    ".main{flex:1;min-width:0;max-width:1280px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}.report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-card h2,.report-section h2{margin-top:0;}.report-card,.report-section,.report-reduction-section,.report-subtype-section,.report-variant-section,.report-subsection,.report-parameter-block{scroll-margin-top:24px;}",
    ".report-meta{margin:0;color:#516274;font-size:14px;}.report-empty{color:#657789;font-style:italic;}.table-block{scroll-margin-top:24px;}.table-caption{font-size:12px;line-height:1.45;color:#334e68;margin:8px 0 6px;}.report-table{width:100%;border-collapse:collapse;font-size:13px;margin:10px 0 18px 0;}.report-table th,.report-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}.report-table th{background:#f7f9fb;font-weight:700;}",
    ".report-parameter-block{margin:8px 0 22px;padding:12px 14px;border:1px solid #e1e8f0;border-radius:8px;background:#fbfdff;}.report-reduction-section{margin:26px 0 30px;padding-top:2px;}.report-reduction-section h3{margin:0 0 12px;font-size:21px;color:#17324c;border-bottom:1px solid #e4ebf2;padding-bottom:6px;}.report-subtype-section{margin:18px 0 24px;}.report-subtype-section h4{margin:0 0 10px;font-size:17px;color:#243b52;}.report-variant-section{margin:18px 0 24px;padding:14px;border:1px solid #edf2f7;border-radius:8px;background:#fff;}.report-variant-section h4,.report-variant-section h5{margin:0 0 6px;font-size:15px;color:#334e68;}.report-subsection{margin:18px 0 30px;padding-top:2px;}.report-subsection h5{margin:0 0 8px;font-size:15px;color:#17324c;}.report-subsection h6{margin:0 0 8px;font-size:14px;color:#17324c;}ul{margin-top:8px;}.report-figure-grid{display:grid;gap:18px;margin:16px 0 24px 0;align-items:stretch;}.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}",
    ".report-table-grid{display:grid;gap:18px;margin:10px 0 18px 0;align-items:start;}.report-table-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-table-grid h5{margin:0 0 6px 0;font-size:13px;color:#284662;}",
    ".report-figure-card{min-width:0;}.report-figure{margin:0 0 10px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;aspect-ratio:1/1;display:flex;align-items:center;justify-content:center;overflow:hidden;}.report-figure-image{width:100%;height:100%;display:block;object-fit:contain;border-radius:6px;}.report-figure-title,.report-figure-legend{margin:8px 0 0 0;}.report-figure-title{font-size:13px;line-height:1.35;}.report-figure-legend{font-size:12px;line-height:1.4;color:#536577;}",
    "@media (max-width: 1100px){.shell{display:block;}.sidebar{position:static;width:auto;margin-bottom:20px;}}@media (max-width: 900px){.report-figure-grid--2,.report-figure-grid--3{grid-template-columns:1fr;}}"
  )

  methods <- c(
    "Parameters on positive multiplicative scales are transformed as log10(x) before analysis; all other variables are retained on their original scale.",
    "For each reduction input column x_j, values are standardized as z_ij = (x_ij - mean(x_j)) / sd(x_j), unless a prior-unit variant is requested. Prior-unit variants rescale transformed optimizer coordinates to their prior interval before the reduction.",
    if (length(reduction_order) == 1L && identical(reduction_order, "pca")) {
      "This PCA-only report includes true PCA coordinate plots. The x- and y-axes are PCA 1 and PCA 2."
    } else {
      "This report distinguishes direct UMAPs, true PCA coordinate plots, and true t-SNE embeddings. PCA-to-UMAP sensitivity outputs remain under the UMAP group because their final displayed coordinates are UMAP coordinates."
    },
    paste0(
      "Each displayed reduction is annotated with two clustering analyses. One estimates clusters from the variables used to construct the reduction, and the other estimates clusters from the final two-dimensional ",
      reduction_phrase,
      " coordinates."
    ),
    "For k = 2,...,8, k-means clustering is evaluated by the mean silhouette width. The selected k is the value with the largest mean silhouette width. For large data sets, this criterion is evaluated on a fixed random subset of at most 5000 observations, and the selected k is then used to assign clusters for all observations.",
    "Cluster color denotes cluster membership and does not modify the objective-value color scale.",
    "For in vivo reductions, point shape denotes the fixed-O2 attractor mode at O2 = 2. Mode 1 is shown with circles and Mode 2 is shown with triangles when that annotation is available.",
    "The pooled in vivo/in vitro reductions use the 14 parameters shared by the two data-set-specific representations and keep the in vivo and in vitro objective values on their original scales."
  )

  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
	    "<title>", html_escape(report_title), "</title>",
	    "<style>", css, "</style></head><body>",
	    "<div class=\"shell\">",
	    "<aside class=\"sidebar\"><div class=\"side-head\"><div class=\"kicker\">parameter landscape clustering</div>",
	    "<div class=\"side-title\">", html_escape(side_title), "</div><div class=\"side-subtitle\">pooled / in vivo / in vitro</div></div><nav>",
	    nav,
	    "</nav></aside><main class=\"main\">",
	    "<section class=\"report-card\" id=\"report-metadata\"><h1>", html_escape(report_title), "</h1>",
	    "<p class=\"report-meta\">Clustered ", html_escape(reduction_phrase), " summaries are organized as pooled, in vivo, and in vitro sections.</p>",
	    "<p class=\"report-meta\"><strong>Results root:</strong> ", html_escape(normalizePath(root_dir, mustWork = FALSE)), "<br>",
	    "<strong>Generated at:</strong> ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p></section>",
	    method_block("Methods", methods, "methods"),
	    report_body,
	    "</main></div>",
	    report_nav_script(),
	    "</body></html>"
	  )
	}

write_umap_cluster_report <- function(result_root = default_parameter_landscape_clustering_dir(),
                                      output_html = NULL,
                                      reductions = c("pca", "umap", "tsne"),
                                      report_title = NULL) {
  reductions <- normalize_report_reductions(reductions)
  result_root <- normalizePath(path.expand(result_root), mustWork = FALSE)
  output_html <- output_html %||% file.path(result_root, cluster_report_filename(reductions))
  output_html <- normalizePath(path.expand(output_html), mustWork = FALSE)
  if (!dir.exists(result_root)) stop("Result root does not exist: ", result_root)
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  html <- build_report_html(result_root, reductions = reductions, report_title = report_title)
  writeLines(html, con = output_html, useBytes = TRUE)
  message("Wrote cluster report: ", output_html)
  invisible(output_html)
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
result_root <- argv$result_root %||% default_parameter_landscape_clustering_dir()
reductions <- normalize_report_reductions(argv$reductions %||% argv$reduction %||% c("pca", "umap", "tsne"))
output_html <- argv$output_html %||% file.path(result_root, cluster_report_filename(reductions))
write_umap_cluster_report(
  result_root = result_root,
  output_html = output_html,
  reductions = reductions,
  report_title = argv$report_title
)
