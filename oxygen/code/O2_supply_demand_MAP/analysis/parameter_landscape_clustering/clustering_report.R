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
  if (!identical(dataset, "pooled_invivo_invitro")) return("")
  unique(c(pooled_prefix_output_subdir(prefix), ""))
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
  header <- paste0("<th>", html_escape(names(vals)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(vals)), function(i) {
    paste0("<tr>", paste0("<td>", html_escape(unlist(vals[i, , drop = TRUE], use.names = FALSE)), "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0("<table class=\"report-table\"><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
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
    png = first_existing_path(report_candidate_paths(fig_dir, dataset, prefix, source_dir, paste0(prefix, ".png"))),
    pdf = first_existing_path(report_candidate_paths(fig_dir, dataset, prefix, source_dir, paste0(prefix, ".pdf")))
  )
}

base_figure_paths <- function(root_dir, dataset, prefix, reduction = "umap") {
  fig_dir <- report_figures_dir(root_dir, dataset, reduction = reduction)
  list(
    png = first_existing_path(report_candidate_paths(fig_dir, dataset, prefix, paste0(prefix, ".png"))),
    pdf = first_existing_path(report_candidate_paths(fig_dir, dataset, prefix, paste0(prefix, ".pdf")))
  )
}

cluster_output_pair_exists <- function(root_dir, dataset, spec) {
  reduction <- spec_reduction(spec)
  base_fig <- base_figure_paths(root_dir, dataset, spec$prefix, reduction = reduction)
  source_dirs <- c("ByInputFeatures", spec_coordinate_source_dir(spec))
  file.exists(base_fig$png) &&
    file.exists(base_fig$pdf) &&
    all(vapply(source_dirs, function(src) {
    fig <- figure_paths(root_dir, dataset, src, spec$prefix, reduction = reduction)
    file.exists(fig$png) &&
      file.exists(fig$pdf) &&
      file.exists(silhouette_path(root_dir, dataset, src, spec$prefix, reduction = reduction))
  }, logical(1)))
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
  coord_path <- file.path(paper_tables_wclusters_dir("invivo", root_dir = root_dir), source_dir, paste0(prefix, "_coordinates.csv"))
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
  reduction <- spec_reduction(spec)
  source_dirs <- c("ByInputFeatures", spec_coordinate_source_dir(spec))
  coordinate_label <- spec_coordinate_source_label(spec)
  source_labels <- c(
    "no cluster overlay",
    "clusters estimated from the input variables",
    paste0("clusters estimated from the ", coordinate_label, " coordinates")
  )
  fig_labels <- paste0("Figure ", figure_index, LETTERS[1:3], ". ", spec$title, " - ", source_labels)
  base_fig <- base_figure_paths(root_dir, dataset, spec$prefix, reduction = reduction)
  figs <- c(
    list(base_fig),
    lapply(source_dirs, function(src) figure_paths(root_dir, dataset, src, spec$prefix, reduction = reduction))
  )
  silhouettes <- do.call(rbind, lapply(source_dirs, function(src) {
    read_silhouette_selected(root_dir, dataset, src, spec$prefix, reduction = reduction)
  }))
  legend_base <- paste0(
    spec$legend,
    spec$legend_suffix %||% " For best-fit seed results, color represents the objective value, with lower values shown in blue and higher values shown in yellow. The outline color identifies cluster membership."
  )
  paste0(
    "<section class=\"report-subsection\" id=\"", html_escape(spec$id), "\">",
    "<h3>", html_escape(paste0(figure_index, ". ", spec$title)), "</h3>",
    "<p>", html_escape(spec$description), "</p>",
    if (length(spec$bullets %||% character())) as_html_list(spec$bullets) else "",
    "<h4>Number of clusters selected by the silhouette criterion</h4>",
    table_to_html(silhouettes, max_rows = 10L),
    extra_summary_block(root_dir, dataset, spec),
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

  figure_counter <- 0L
  next_figure_index <- function() {
    figure_counter <<- figure_counter + 1L
    figure_counter
  }

  render_reduction_group <- function(scope_id, dataset, reduction, specs) {
    group_id <- paste0(scope_id, "-", reduction, "s")
    title <- reduction_group_title(reduction)
    if (!length(specs)) {
      return(paste0(
        "<section id=\"", html_escape(group_id), "\"><h3>", html_escape(title), "</h3>",
        "<p class=\"report-empty\">No clustered ", html_escape(title), " outputs were found for this section.</p></section>"
      ))
    }
    paste0(
      "<section id=\"", html_escape(group_id), "\"><h3>", html_escape(title), "</h3>",
      paste0(vapply(specs, function(spec) {
        pair_block(root_dir, dataset, spec, next_figure_index())
      }, character(1)), collapse = ""),
      "</section>"
    )
  }

  render_scope_section <- function(scope_id,
                                   dataset,
                                   title,
                                   parameter_heading,
                                   parameter_text,
                                   parameter_df,
                                   specs) {
    paste0(
      "<section class=\"report-section\" id=\"", html_escape(scope_id), "\"><h2>", html_escape(title), "</h2>",
      "<section id=\"", html_escape(paste0(scope_id, "-parameters")), "\"><h3>", html_escape(parameter_heading), "</h3>",
      "<p>", html_escape(parameter_text), "</p>",
      table_to_html(parameter_df, max_rows = 50L), "</section>",
      paste0(vapply(c("umap", "pca", "tsne"), function(reduction) {
        render_reduction_group(scope_id, dataset, reduction, specs[[reduction]])
      }, character(1)), collapse = ""),
      "</section>"
    )
  }

  scope_nav <- function(scope_id, title, specs) {
    c(
      list(list(id = scope_id, title = title)),
      list(list(id = paste0(scope_id, "-parameters"), title = paste0(title, " parameters"), class = "report-nav-link report-nav-h3")),
      lapply(c("umap", "pca", "tsne"), function(reduction) {
        count <- length(specs[[reduction]])
        list(
          id = paste0(scope_id, "-", reduction, "s"),
          title = paste0(title, " ", reduction_group_title(reduction), " (", count, ")"),
          class = "report-nav-link report-nav-h3"
        )
      })
    )
  }

  nav <- build_nav(c(
    list(
      list(id = "report-metadata", title = "Overview"),
      list(id = "methods", title = "Methods")
    ),
    scope_nav("pooled", "pooled", pooled_specs),
    scope_nav("invivo", "in vivo", invivo_specs),
    scope_nav("invitro", "in vitro", invitro_specs)
  ))

  css <- paste0(
    "html{scroll-behavior:smooth;}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f4f7fa;color:#1b2a38;}",
    ".report-shell{display:flex;gap:28px;max-width:1680px;margin:0 auto;padding:24px;}.report-sidebar{position:sticky;top:24px;align-self:flex-start;width:310px;max-height:calc(100vh - 48px);border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);overflow:auto;scrollbar-gutter:stable;}",
    ".report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}.report-kicker{font-size:11px;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;opacity:0.78;}.report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}.report-subtitle{margin-top:4px;font-size:12px;opacity:0.85;}",
    ".report-nav{padding:10px 8px 12px 8px;}.report-nav-list{margin:0;padding:0;list-style:none;}.report-nav-item{margin:4px 0;}.report-nav-link{display:block;padding:10px 12px;border-radius:8px;text-decoration:none;color:#17324c;font-size:14px;font-weight:600;line-height:1.35;}.report-nav-link:hover{background:rgba(47,110,164,0.08);}.report-nav-h3{padding-left:26px;font-size:13px;font-weight:500;color:#536577;}",
    ".report-main{flex:1;min-width:0;max-width:1180px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}.report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-card h2,.report-section h2{margin-top:0;}.report-card h2,.report-section h2,.report-subsection h3{scroll-margin-top:24px;}",
    ".report-meta{margin:0;color:#516274;font-size:14px;}.report-empty{color:#657789;font-style:italic;}.report-table{width:100%;border-collapse:collapse;font-size:13px;margin:10px 0 18px 0;}.report-table th,.report-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}.report-table th{background:#f7f9fb;font-weight:700;}",
    ".report-subsection{margin:22px 0 34px 0;padding-top:2px;}.report-subsection h3{margin:0 0 8px 0;font-size:19px;color:#17324c;}.report-subsection h4{margin:16px 0 4px 0;font-size:14px;color:#284662;}ul{margin-top:8px;}.report-figure-grid{display:grid;gap:18px;margin:16px 0 24px 0;align-items:stretch;}.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}",
    ".report-table-grid{display:grid;gap:18px;margin:10px 0 18px 0;align-items:start;}.report-table-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-table-grid h5{margin:0 0 6px 0;font-size:13px;color:#284662;}",
    ".report-figure-card{min-width:0;}.report-figure{margin:0 0 10px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;aspect-ratio:1/1;display:flex;align-items:center;justify-content:center;overflow:hidden;}.report-figure-image{width:100%;height:100%;display:block;object-fit:contain;border-radius:6px;}.report-figure-title,.report-figure-legend{margin:8px 0 0 0;}.report-figure-title{font-size:13px;line-height:1.35;}.report-figure-legend{font-size:12px;line-height:1.4;color:#536577;}.report-figure-link{font-size:11px;margin:6px 0 0 0;}.report-figure-link a{color:#2f6ea4;text-decoration:none;font-weight:600;}",
    "@media (max-width: 1100px){.report-shell{display:block;}.report-sidebar{position:static;width:auto;margin-bottom:20px;}}@media (max-width: 900px){.report-figure-grid--2,.report-figure-grid--3{grid-template-columns:1fr;}}"
  )

  methods <- c(
    "Parameters on positive multiplicative scales are transformed as log10(x) before analysis; all other variables are retained on their original scale.",
    "For each reduction input column x_j, values are standardized as z_ij = (x_ij - mean(x_j)) / sd(x_j), unless a prior-unit variant is requested. Prior-unit variants rescale transformed optimizer coordinates to their prior interval before the reduction.",
    "This report distinguishes direct UMAPs, true PCA coordinate plots, and true t-SNE embeddings. PCA-to-UMAP sensitivity outputs remain under the UMAP group because their final displayed coordinates are UMAP coordinates.",
    "Each displayed reduction is annotated with two clustering analyses. One estimates clusters from the variables used to construct the reduction, and the other estimates clusters from the final two-dimensional UMAP, PCA, or t-SNE coordinates.",
    "For k = 2,...,8, k-means clustering is evaluated by the mean silhouette width. The selected k is the value with the largest mean silhouette width. For large data sets, this criterion is evaluated on a fixed random subset of at most 5000 observations, and the selected k is then used to assign clusters for all observations.",
    "Cluster color denotes cluster membership and does not modify the objective-value color scale.",
    "For in vivo reductions, point shape denotes the fixed-O2 attractor mode at O2 = 2. Mode 1 is shown with circles and Mode 2 is shown with triangles when that annotation is available.",
    "The pooled in vivo/in vitro reductions use the 14 parameters shared by the two data-set-specific representations and keep the in vivo and in vitro objective values on their original scales."
  )

  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Parameter Landscape Reduction Report</title>",
    "<style>", css, "</style></head><body>",
    "<div class=\"report-shell\">",
    "<aside class=\"report-sidebar\"><div class=\"report-sidebar-header\"><div class=\"report-kicker\">parameter landscape clustering</div>",
    "<div class=\"report-title\">Reduction Report</div><div class=\"report-subtitle\">pooled, in vivo, and in vitro</div></div><nav class=\"report-nav\">",
    nav,
    "</nav></aside><main class=\"report-main\">",
    "<section class=\"report-card\" id=\"report-metadata\"><h1>Parameter Landscape Reduction Report</h1>",
    "<p class=\"report-meta\">Clustered UMAP, PCA, and t-SNE summaries are organized as pooled, in vivo, and in vitro sections.</p>",
    "<p class=\"report-meta\"><strong>Results root:</strong> ", html_escape(normalizePath(root_dir, mustWork = FALSE)), "<br>",
    "<strong>Generated at:</strong> ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p></section>",
    method_block("Methods", methods, "methods"),
    render_scope_section(
      scope_id = "pooled",
      dataset = "pooled_invivo_invitro",
      title = "pooled in vivo/in vitro",
      parameter_heading = "Pooled reduction input parameters",
      parameter_text = "The pooled reductions use the intersection of the in vivo and in vitro parameter sets.",
      parameter_df = pooled_parameter_table(),
      specs = pooled_specs
    ),
    render_scope_section(
      scope_id = "invivo",
      dataset = "invivo",
      title = "in vivo",
      parameter_heading = "in vivo reduction input parameters",
      parameter_text = "The base in vivo reductions use 18 multi-paired soft-coupling parameters. The table records whether each parameter is log10-transformed before standardization.",
      parameter_df = parameter_table("invivo"),
      specs = invivo_specs
    ),
    render_scope_section(
      scope_id = "invitro",
      dataset = "invitro",
      title = "in vitro",
      parameter_heading = "in vitro reduction input parameters",
      parameter_text = "The in vitro reductions use 14 soft-coupling parameters and do not include any ploidy trajectory-derived features.",
      parameter_df = parameter_table("invitro"),
      specs = invitro_specs
    ),
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
