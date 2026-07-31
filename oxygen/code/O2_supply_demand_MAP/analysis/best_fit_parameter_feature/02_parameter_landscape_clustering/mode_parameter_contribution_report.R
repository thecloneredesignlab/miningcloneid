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

format_report_number <- function(x, digits = 3L, integer_tol = 1e-9) {
  x <- suppressWarnings(as.numeric(x))
  vapply(x, function(v) {
    if (!is.finite(v)) return("")
    if (abs(v - round(v)) < integer_tol) return(format(round(v), scientific = FALSE, trim = TRUE))
    trimws(formatC(v, digits = digits, format = "fg"))
  }, character(1))
}

fmt_num <- function(x, digits = 3L) {
  format_report_number(x, digits = digits)
}

fmt_o2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- vapply(x, function(v) {
    if (!is.finite(v)) return("")
    txt <- sprintf("%.10f", v)
    txt <- sub("0+$", "", txt)
    txt <- sub("\\.$", "", txt)
    if (identical(txt, "-0")) txt <- "0"
    txt
  }, character(1))
  out
}

display_o2 <- function(x) {
  gsub("O2", "O\u2082", as.character(x), fixed = TRUE)
}

fmt_o2_label <- function(x) {
  paste0("O\u2082=", fmt_o2(x))
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read_csv_plain(path)
}

read_tsv_if_exists <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read_tsv(path)
}

table_to_html <- function(df, max_rows = 30L, class = "report-table") {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    return("<p class=\"report-empty\">No rows available.</p>")
  }
  df <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  vals <- as.data.frame(lapply(df, function(x) {
    if (is.numeric(x)) {
      format_report_number(x, digits = 4L)
    } else {
      display_o2(ifelse(is.na(x), "", as.character(x)))
    }
  }), stringsAsFactors = FALSE, check.names = FALSE)
  names(vals) <- display_o2(names(vals))
  header <- paste0("<th>", html_escape(names(vals)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(vals)), function(i) {
    cells <- paste0("<td>", html_escape(unlist(vals[i, , drop = TRUE], use.names = FALSE)), "</td>", collapse = "")
    paste0("<tr>", cells, "</tr>")
  }, character(1))
  paste0(
    "<table class=\"", class, "\"><thead><tr>", header, "</tr></thead><tbody>",
    paste(rows, collapse = ""), "</tbody></table>"
  )
}

new_report_context <- function() {
  ctx <- new.env(parent = emptyenv())
  ctx$table_counter <- 0L
  ctx$figure_counter <- 0L
  ctx$table_section_counts <- new.env(parent = emptyenv())
  ctx$figure_section_counts <- new.env(parent = emptyenv())
  ctx$table_nav <- list()
  ctx$figure_nav <- list()
  ctx$toc_headings <- list()
  ctx$toc_refs <- list()
  ctx$current_section_id <- ""
  ctx$current_section_number <- ""
  ctx$current_heading_id <- ""
  ctx
}

kv_data <- function(keys, values) {
  data.frame(item = keys, value = values, check.names = FALSE, stringsAsFactors = FALSE)
}

next_section_item_number <- function(ctx, counter_env) {
  section_number <- as.character(ctx$current_section_number %||% "")
  if (!nzchar(section_number)) section_number <- "0"
  current <- if (exists(section_number, envir = counter_env, inherits = FALSE)) {
    get(section_number, envir = counter_env, inherits = FALSE)
  } else {
    0L
  }
  current <- current + 1L
  assign(section_number, current, envir = counter_env)
  paste0(section_number, ".", current)
}

table_block <- function(ctx, df, title, caption, conclusion, max_rows = 30L, class = "report-table") {
  ctx$table_counter <- ctx$table_counter + 1L
  table_number <- next_section_item_number(ctx, ctx$table_section_counts)
  table_id <- paste0("table-", ctx$table_counter)
  ctx$table_nav[[length(ctx$table_nav) + 1L]] <- list(
    id = table_id,
    label = paste0("Table ", table_number, ". ", title)
  )
  ctx$toc_refs[[length(ctx$toc_refs) + 1L]] <- list(
    id = table_id,
    label = paste0("Table ", table_number, ". ", title),
    section_id = ctx$current_section_id,
    heading_id = ctx$current_heading_id,
    type = "table"
  )
  shown_note <- ""
  if (is.data.frame(df) && nrow(df) > max_rows) {
    shown_note <- sprintf(" The HTML table displays the first %d of %d rows.", max_rows, nrow(df))
  }
  table_title <- paste0("Table ", table_number, ". ", title, ".")
  paste0(
    "<div class=\"table-block\" id=\"", html_escape(table_id), "\" data-nav-target=\"", html_escape(table_id), "\">",
    "<p class=\"table-caption\"><strong>", html_escape(display_o2(table_title)), "</strong> ", html_escape(display_o2(paste0(caption, shown_note))), "</p>",
    table_to_html(df, max_rows = max_rows, class = class),
    if (nzchar(conclusion)) {
      paste0("<p class=\"table-conclusion\">", html_escape(display_o2(conclusion)), "</p>")
    } else {
      ""
    },
    "</div>"
  )
}

file_to_data_uri <- function(path, mime = "image/png") {
  if (!file.exists(path)) return(NA_character_)
  if (!requireNamespace("base64enc", quietly = TRUE)) return(NA_character_)
  paste0("data:", mime, ";base64,", base64enc::base64encode(path))
}

asset_src <- function(path, output_html, embed_assets = TRUE) {
  if (!file.exists(path)) return(NA_character_)
  if (isTRUE(embed_assets)) {
    uri <- file_to_data_uri(path)
    if (!is.na(uri)) return(uri)
  }
  rel <- tryCatch(
    utils::URLencode(normalizePath(path, mustWork = FALSE), reserved = TRUE),
    error = function(e) path
  )
  paste0("file://", rel)
}

image_card <- function(ctx, path, title, caption = "", conclusion = "", output_html, embed_assets = TRUE) {
  ctx$figure_counter <- ctx$figure_counter + 1L
  figure_number <- next_section_item_number(ctx, ctx$figure_section_counts)
  figure_id <- paste0("figure-", ctx$figure_counter)
  ctx$figure_nav[[length(ctx$figure_nav) + 1L]] <- list(
    id = figure_id,
    label = paste0("Figure ", figure_number, ". ", title)
  )
  ctx$toc_refs[[length(ctx$toc_refs) + 1L]] <- list(
    id = figure_id,
    label = paste0("Figure ", figure_number, ". ", title),
    section_id = ctx$current_section_id,
    heading_id = ctx$current_heading_id,
    type = "figure"
  )
  figure_title <- paste0("Figure ", figure_number, ". ", title, ".")
  caption_html <- paste0(
    "<p class=\"figure-caption\"><strong>", html_escape(display_o2(figure_title)), "</strong>",
    if (nzchar(caption)) paste0(" ", html_escape(display_o2(caption))) else "",
    "</p>"
  )
  interpretation_html <- if (nzchar(conclusion)) {
    paste0("<p class=\"figure-interpretation\">", html_escape(display_o2(conclusion)), "</p>")
  } else {
    ""
  }
  if (!file.exists(path)) {
    return(paste0(
      "<article class=\"figure-card missing\" id=\"", html_escape(figure_id), "\" data-nav-target=\"", html_escape(figure_id), "\">",
      caption_html,
      "<p>Missing image: ", html_escape(path), "</p>",
      interpretation_html,
      "</article>"
    ))
  }
  src <- asset_src(path, output_html = output_html, embed_assets = embed_assets)
  paste0(
    "<article class=\"figure-card\" id=\"", html_escape(figure_id), "\" data-nav-target=\"", html_escape(figure_id), "\">",
    "<div class=\"figure-frame\"><img src=\"", html_escape(src), "\" alt=\"", html_escape(display_o2(figure_title)), "\"></div>",
    caption_html,
    interpretation_html,
    "</article>"
  )
}

section <- function(id, number, title, body) {
  paste0(
    "<section class=\"report-section\" id=\"", html_escape(id), "\" data-nav-target=\"", html_escape(id), "\">",
    "<h2>", html_escape(number), ". ", html_escape(display_o2(title)), "</h2>",
    body,
    "</section>"
  )
}

subsection <- function(id, number, title, body) {
  paste0(
    "<section class=\"report-subsection\" id=\"", html_escape(id), "\" data-nav-target=\"", html_escape(id), "\">",
    "<h3>", html_escape(number), " ", html_escape(display_o2(title)), "</h3>",
    body,
    "</section>"
  )
}

numbered_heading <- function(ctx, level, number, title) {
  heading_id <- html_id(paste0("heading-", number, "-", title))
  ctx$toc_headings[[length(ctx$toc_headings) + 1L]] <- list(
    id = heading_id,
    label = paste0(number, " ", title),
    section_id = ctx$current_section_id
  )
  ctx$current_heading_id <- heading_id
  paste0(
    "<h", level, " id=\"", html_escape(heading_id), "\" data-nav-target=\"", html_escape(heading_id), "\">",
    html_escape(number), " ", html_escape(display_o2(title)),
    "</h", level, ">"
  )
}

o2_slug <- function(x) {
  fixed_o2_o2_slug(suppressWarnings(as.numeric(x)))
}

reference_dir <- function(mode_dir, o2) {
  file.path(mode_dir, paste0("reference_o2_", o2_slug(o2)))
}

reference_paths <- function(mode_dir, o2) {
  ref_dir <- reference_dir(mode_dir, o2)
  list(
    ref_dir = ref_dir,
    summary = file.path(ref_dir, "mode_parameter_contribution_summary.tsv"),
    status = file.path(ref_dir, "mode_parameter_elastic_net_status.tsv"),
    feature_importance = file.path(ref_dir, "mode_parameter_feature_importance.csv"),
    auc_by_feature = file.path(ref_dir, "auc_plots", "mode_parameter_auc_by_feature.csv"),
    seed_labels = file.path(ref_dir, "mode_parameter_seed_labels.csv"),
    decision_rules = file.path(ref_dir, "mode_parameter_decision_rules.csv"),
    top3_joint = file.path(ref_dir, "top3_joint", "top3_joint_auc.csv"),
    main_auc_png = file.path(ref_dir, "auc_plots", "main_feature_auc_bar.png"),
    top30_auc_png = file.path(ref_dir, "auc_plots", "combined_feature_auc_bar_top30.png"),
    top_feature_roc_png = file.path(ref_dir, "top3_feature_roc", paste0("top", 1:3, "_feature_roc.png")),
    top3_roc_png = file.path(ref_dir, "top3_joint", "top3_joint_roc.png"),
    cumulative_png = file.path(ref_dir, "cumulative_auc", "cumulative_feature_auc.png"),
    cumulative_csv = file.path(ref_dir, "cumulative_auc", "cumulative_feature_auc.csv")
  )
}

fixed_o2_distribution_paths <- function(mode_dir) {
  fig_dir <- file.path(dirname(normalizePath(mode_dir, mustWork = FALSE)), "FixO2Modes", "Figures")
  list(
    all_seeds_png = file.path(fig_dir, "fixed_o2_analytical_solution_all_seeds_violin_boxplot.png"),
    by_mode_png = file.path(fig_dir, "fixed_o2_analytical_solution_by_mode_violin_boxplot.png")
  )
}

safe_col <- function(df, col, default = NA) {
  if (!is.data.frame(df) || !col %in% names(df)) return(rep(default, max(1L, nrow(df))))
  df[[col]]
}

top_feature_table <- function(feature_importance, auc_by_feature, n = 10L) {
  if (!nrow(feature_importance)) return(data.frame())
  keep <- feature_importance[order(suppressWarnings(as.numeric(feature_importance$rank))), , drop = FALSE]
  keep <- keep[seq_len(min(nrow(keep), n)), , drop = FALSE]
  auc_cols <- c("feature_name", "auc_discriminative", "direction", "auc_mode1_higher", "auc_rank")
  auc_cols <- intersect(auc_cols, names(auc_by_feature))
  if (length(auc_cols) && nrow(auc_by_feature)) {
    keep <- merge(keep, auc_by_feature[, auc_cols, drop = FALSE], by = "feature_name", all.x = TRUE, sort = FALSE)
  }
  data.frame(
    rank = suppressWarnings(as.integer(keep$rank)),
    feature = keep$feature_name,
    type = keep$feature_type,
    parameter_a = keep$parameter_a,
    parameter_b = keep$parameter_b,
    elastic_net_coef = suppressWarnings(as.numeric(keep$elastic_net_coef)),
    abs_coef = suppressWarnings(as.numeric(keep$abs_elastic_net_coef)),
    selection_frequency = suppressWarnings(as.numeric(keep$selection_frequency)),
    auc = suppressWarnings(as.numeric(keep$auc_discriminative)),
    direction = keep$direction,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

feature_recurrence_table <- function(top_features, top_n = 10L) {
  if (!nrow(top_features)) return(data.frame())
  tf <- top_features[suppressWarnings(as.integer(top_features$rank)) <= top_n, , drop = FALSE]
  if (!nrow(tf)) return(data.frame())
  rows <- lapply(split(tf, tf$feature_name), function(df) {
    refs <- sort(unique(suppressWarnings(as.numeric(df$mode_reference_o2))))
    data.frame(
      feature = df$feature_name[[1L]],
      type = df$feature_type[[1L]],
      times_in_top_n = length(refs),
      reference_o2 = paste(fmt_o2(refs), collapse = ", "),
      median_rank = stats::median(suppressWarnings(as.numeric(df$rank)), na.rm = TRUE),
      median_selection_frequency = stats::median(suppressWarnings(as.numeric(df$selection_frequency)), na.rm = TRUE),
      median_abs_coef = stats::median(abs(suppressWarnings(as.numeric(df$elastic_net_coef))), na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$times_in_top_n, out$median_rank, -out$median_selection_frequency, out$feature), , drop = FALSE]
}

top3_summary_table <- function(index) {
  if (!nrow(index)) return(data.frame())
  data.frame(
    reference_O2 = fmt_o2(index$mode_reference_o2),
    n_seed = index$n_seed,
    mode1 = index$n_mode1,
    mode2 = index$n_mode2,
    top1 = index$top1_feature,
    top1_AUC = suppressWarnings(as.numeric(index$top1_auc)),
    top2 = index$top2_feature,
    top2_AUC = suppressWarnings(as.numeric(index$top2_auc)),
    top3 = index$top3_feature,
    top3_AUC = suppressWarnings(as.numeric(index$top3_auc)),
    top3_joint_CV_AUC = suppressWarnings(as.numeric(index$top3_joint_cv_auc)),
    top3_joint_apparent_AUC = suppressWarnings(as.numeric(index$top3_joint_apparent_auc)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mode_balance_table <- function(index) {
  if (!nrow(index)) return(data.frame())
  data.frame(
    reference_O2 = fmt_o2(index$mode_reference_o2),
    n_seed = index$n_seed,
    mode1 = index$n_mode1,
    mode2 = index$n_mode2,
    mode1_fraction = suppressWarnings(as.numeric(index$n_mode1)) / suppressWarnings(as.numeric(index$n_seed)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mode_definition_table <- function() {
  data.frame(
    label = c("Mode1", "Mode2"),
    encoded_response = c(1L, 0L),
    fixed_o2_attractor_rule = c(
      "dominant_mean_ploidy at the reference O2 >= 2",
      "dominant_mean_ploidy at the reference O2 < 2"
    ),
    interpretation = c(
      "The fixed-O2 attractor is dominated by ploidy at or above the diploid threshold.",
      "The fixed-O2 attractor is dominated by ploidy below the diploid threshold."
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

cumulative_highlight_table <- function(cumulative) {
  if (!nrow(cumulative)) return(data.frame())
  ks <- c(1L, 3L, 5L, 10L, 30L)
  out <- cumulative[suppressWarnings(as.integer(cumulative$k)) %in% ks, , drop = FALSE]
  if (!nrow(out)) return(data.frame())
  data.frame(
    k = out$k,
    added_feature = out$added_feature,
    CV_AUC = suppressWarnings(as.numeric(out$cv_auc)),
    CV_AUC_gain = suppressWarnings(as.numeric(out$cv_auc_gain)),
    apparent_AUC = suppressWarnings(as.numeric(out$apparent_auc)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

index_interpretation <- function(index) {
  if (!nrow(index)) return(character())
  best_joint <- index[which.max(suppressWarnings(as.numeric(index$top3_joint_cv_auc))), , drop = FALSE]
  worst_joint <- index[which.min(suppressWarnings(as.numeric(index$top3_joint_cv_auc))), , drop = FALSE]
  top1_features <- paste(unique(index$top1_feature), collapse = ", ")
  c(
    sprintf(
      "Across %d reference O2 values, the analysis used %d seeds per reference O2 with both mode labels present.",
      nrow(index),
      suppressWarnings(as.integer(stats::median(index$n_seed, na.rm = TRUE)))
    ),
    sprintf(
      "The top single retained features across reference O2 values were: %s.",
      top1_features
    ),
    sprintf(
      "The strongest top3 joint cross-validated AUC was %s at reference O2=%s; the weakest was %s at reference O2=%s.",
      fmt_num(best_joint$top3_joint_cv_auc[[1L]], digits = 3L),
      fmt_o2(best_joint$mode_reference_o2[[1L]]),
      fmt_num(worst_joint$top3_joint_cv_auc[[1L]], digits = 3L),
      fmt_o2(worst_joint$mode_reference_o2[[1L]])
    )
  )
}

mode_direction_text <- function(x) {
  x <- as.character(x %||% "")
  if (identical(x, "higher_in_mode1")) return("higher values are associated with Mode1")
  if (identical(x, "higher_in_mode2")) return("higher values are associated with Mode2")
  if (nzchar(x)) return(x)
  "direction was not available"
}

auc_strength_text <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return("unavailable")
  if (x >= 0.9) return("very strong")
  if (x >= 0.8) return("strong")
  if (x >= 0.7) return("moderate")
  "limited"
}

report_metadata_conclusion <- function(index, reference_o2) {
  sprintf(
    "The report covers %d reference-O2 analyses (%s), all generated from the same mode-contribution result directory.",
    length(reference_o2),
    paste(fmt_o2(reference_o2), collapse = ", ")
  )
}

run_args_conclusion <- function(run_args) {
  if (!nrow(run_args) || !"argument" %in% names(run_args) || !"value" %in% names(run_args)) {
    return("Run arguments were not available, so this report should be interpreted as a summary of the existing output files only.")
  }
  vals <- stats::setNames(as.character(run_args$value), as.character(run_args$argument))
  sprintf(
    "The model settings use elastic-net alpha=%s, %s stability resamples, sample fraction=%s, and random seed=%s, which determine sparsity, stability-selection variability, and reproducibility.",
    vals[["elastic_net_alpha"]] %||% "NA",
    vals[["n_bootstrap"]] %||% "NA",
    vals[["sample_fraction"]] %||% "NA",
    vals[["random_seed"]] %||% "NA"
  )
}

mode_balance_conclusion <- function(index) {
  if (!nrow(index)) return("Mode balance was not available.")
  frac <- suppressWarnings(as.numeric(index$n_mode1) / as.numeric(index$n_seed))
  min_i <- which.min(frac)
  max_i <- which.max(frac)
  sprintf(
    "Class balance is usable across reference O2 values: Mode1 fraction ranges from %s at O2=%s to %s at O2=%s, with %d seeds per reference analysis.",
    fmt_num(frac[[min_i]], digits = 3L),
    fmt_o2(index$mode_reference_o2[[min_i]]),
    fmt_num(frac[[max_i]], digits = 3L),
    fmt_o2(index$mode_reference_o2[[max_i]]),
    suppressWarnings(as.integer(stats::median(index$n_seed, na.rm = TRUE)))
  )
}

mode_definition_conclusion <- function() {
  "Mode labels are assigned separately for each reference O2 value using the fixed-O2 attractor dominant mean ploidy; the threshold value 2 is the ploidy threshold, not the O2 concentration."
}

top3_summary_conclusion <- function(index) {
  if (!nrow(index)) return("Top3 feature summary was not available.")
  best <- index[which.max(suppressWarnings(as.numeric(index$top3_joint_cv_auc))), , drop = FALSE]
  worst <- index[which.min(suppressWarnings(as.numeric(index$top3_joint_cv_auc))), , drop = FALSE]
  sprintf(
    "The top3 joint classifier performs best at O2=%s (CV AUC %s) and weakest at O2=%s (CV AUC %s), showing that mode separability depends on the chosen reference O2.",
    fmt_o2(best$mode_reference_o2[[1L]]),
    fmt_num(best$top3_joint_cv_auc[[1L]], digits = 3L),
    fmt_o2(worst$mode_reference_o2[[1L]]),
    fmt_num(worst$top3_joint_cv_auc[[1L]], digits = 3L)
  )
}

feature_recurrence_conclusion <- function(recurrence) {
  if (!nrow(recurrence)) return("Feature recurrence was not available.")
  top <- recurrence[order(-suppressWarnings(as.numeric(recurrence$times_in_top_n)), recurrence$median_rank), , drop = FALSE]
  sprintf(
    "%s is the most recurrent retained feature, appearing in the Top10 for %s reference-O2 analyses; recurrent features are the most stable candidates for describing common mode regimes.",
    top$feature[[1L]],
    fmt_num(top$times_in_top_n[[1L]], digits = 3L)
  )
}

auc_summary_conclusion <- function(auc_summary) {
  if (!nrow(auc_summary) || !"auc" %in% names(auc_summary)) return("AUC summary was not available.")
  best <- auc_summary[which.max(suppressWarnings(as.numeric(auc_summary$auc))), , drop = FALSE]
  metric <- if ("metric_label" %in% names(best)) best$metric_label[[1L]] else best$metric[[1L]]
  sprintf(
    "The highest AUC in this comparison is %s for %s at O2=%s, indicating the strongest mode separation among the summarized AUC metrics.",
    fmt_num(best$auc[[1L]], digits = 3L),
    metric,
    fmt_o2(best$mode_reference_o2[[1L]])
  )
}

joint_summary_conclusion <- function(joint_summary) {
  if (!nrow(joint_summary) || !"cv_auc" %in% names(joint_summary)) return("Top3 joint AUC summary was not available.")
  best <- joint_summary[which.max(suppressWarnings(as.numeric(joint_summary$cv_auc))), , drop = FALSE]
  sprintf(
    "The best joint Top3 model reaches CV AUC %s at O2=%s using %s, %s, and %s, which is %s discrimination of the fixed-O2 mode label.",
    fmt_num(best$cv_auc[[1L]], digits = 3L),
    fmt_o2(best$mode_reference_o2[[1L]]),
    best$top1_feature[[1L]],
    best$top2_feature[[1L]],
    best$top3_feature[[1L]],
    auc_strength_text(best$cv_auc[[1L]])
  )
}

top3_auc_conclusion <- function(top3_auc) {
  if (!nrow(top3_auc) || !"auc" %in% names(top3_auc)) return("Top3 single-feature AUC table was not available.")
  best <- top3_auc[which.max(suppressWarnings(as.numeric(top3_auc$auc))), , drop = FALSE]
  sprintf(
    "The strongest retained single-feature AUC is %s for %s at O2=%s; single-feature AUCs show which individual retained variables are most discriminative before joint modeling.",
    fmt_num(best$auc[[1L]], digits = 3L),
    best$feature_name[[1L]],
    fmt_o2(best$mode_reference_o2[[1L]])
  )
}

global_top3_figure_conclusion <- function(index) {
  if (!nrow(index)) return("The grouped bar figure could not be summarized because the index table was missing.")
  best <- index[which.max(suppressWarnings(as.numeric(index$top1_auc))), , drop = FALSE]
  sprintf(
    "Among the retained single features, the largest Top1 AUC is %s for %s at O2=%s.",
    fmt_num(best$top1_auc[[1L]], digits = 3L),
    best$top1_feature[[1L]],
    fmt_o2(best$mode_reference_o2[[1L]])
  )
}

global_cumulative_figure_conclusion <- function(cumulative_summary) {
  if (!nrow(cumulative_summary) || !"cv_auc" %in% names(cumulative_summary)) {
    return("The cumulative AUC figure could not be summarized because the cumulative table was missing.")
  }
  max_k <- max(suppressWarnings(as.integer(cumulative_summary$k)), na.rm = TRUE)
  final <- cumulative_summary[suppressWarnings(as.integer(cumulative_summary$k)) == max_k, , drop = FALSE]
  best <- final[which.max(suppressWarnings(as.numeric(final$cv_auc))), , drop = FALSE]
  sprintf(
    "At the largest cumulative feature count (Top%d), O2=%s has the highest CV AUC (%s), indicating the best cumulative separability among reference-O2 settings.",
    max_k,
    fmt_o2(best$mode_reference_o2[[1L]]),
    fmt_num(best$cv_auc[[1L]], digits = 3L)
  )
}

status_conclusion <- function(status, o2) {
  if (!nrow(status)) return(sprintf("No fit-status rows were available for O2=%s.", fmt_o2(o2)))
  status_val <- if ("fit_status" %in% names(status)) status$fit_status[[1L]] else "available"
  sprintf(
    "The elastic-net fit status for O2=%s is %s; feature rankings should be interpreted only for reference analyses with successful model fitting.",
    fmt_o2(o2),
    status_val
  )
}

top_feature_conclusion <- function(top_table, o2) {
  if (!nrow(top_table)) return(sprintf("No retained features were available for O2=%s.", fmt_o2(o2)))
  first <- top_table[1L, , drop = FALSE]
  sprintf(
    "At O2=%s, %s is the top retained feature (AUC %s, selection frequency %s); %s.",
    fmt_o2(o2),
    first$feature[[1L]],
    fmt_num(first$auc[[1L]], digits = 3L),
    fmt_num(first$selection_frequency[[1L]], digits = 3L),
    mode_direction_text(first$direction[[1L]])
  )
}

top_feature_roc_conclusion <- function(top_table, o2, rank) {
  if (!nrow(top_table) || rank > nrow(top_table)) {
    return(sprintf("No Top%d retained feature ROC was available for O2=%s.", rank, fmt_o2(o2)))
  }
  row <- top_table[rank, , drop = FALSE]
  sprintf(
    "The Top%d retained feature at O2=%s is %s, with single-feature AUC %s; %s.",
    rank,
    fmt_o2(o2),
    row$feature[[1L]],
    fmt_num(row$auc[[1L]], digits = 3L),
    mode_direction_text(row$direction[[1L]])
  )
}

joint_table_conclusion <- function(joint_table, o2) {
  if (!nrow(joint_table)) return(sprintf("No Top3 joint model row was available for O2=%s.", fmt_o2(o2)))
  sprintf(
    "The Top3 joint model at O2=%s reaches CV AUC %s using %s, %s, and %s, giving %s supervised discrimination.",
    fmt_o2(o2),
    fmt_num(joint_table$CV_AUC[[1L]], digits = 3L),
    joint_table$top1[[1L]],
    joint_table$top2[[1L]],
    joint_table$top3[[1L]],
    auc_strength_text(joint_table$CV_AUC[[1L]])
  )
}

cumulative_conclusion <- function(cumulative, o2) {
  if (!nrow(cumulative)) return(sprintf("No cumulative AUC rows were available for O2=%s.", fmt_o2(o2)))
  k <- suppressWarnings(as.integer(cumulative$k))
  auc <- suppressWarnings(as.numeric(cumulative$cv_auc))
  top1_i <- which(k == 1L)[1L]
  top3_i <- which(k == 3L)[1L]
  final_i <- which.max(k)
  if (is.na(top1_i) || is.na(top3_i) || is.na(final_i)) {
    return(sprintf("The cumulative table records how model AUC changes as ranked features are added for O2=%s.", fmt_o2(o2)))
  }
  sprintf(
    "For O2=%s, CV AUC increases from %s at Top1 to %s at Top3 and %s at Top%d, showing the incremental value of adding retained features.",
    fmt_o2(o2),
    fmt_num(auc[[top1_i]], digits = 3L),
    fmt_num(auc[[top3_i]], digits = 3L),
    fmt_num(auc[[final_i]], digits = 3L),
    k[[final_i]]
  )
}

rules_conclusion <- function(rules, o2) {
  if (!nrow(rules)) return(sprintf("No presentation rules were available for O2=%s.", fmt_o2(o2)))
  sprintf(
    "The rule table gives %d compact presentation rules for O2=%s; these rules summarize visible regimes but do not replace the elastic-net ranking.",
    nrow(rules),
    fmt_o2(o2)
  )
}

output_files_conclusion <- function() {
  "These files provide the audit trail for the report: global summaries, per-reference ranked features, joint AUCs, cumulative AUCs, and the two global summary figures embedded in the HTML."
}

figure_grid <- function(cards, columns = 2L) {
  paste0("<div class=\"figure-grid figure-grid--", columns, "\">", paste(cards, collapse = ""), "</div>")
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

toc_section_node <- function(section_id, section_label, toc_headings, toc_refs, child_html = "", open = FALSE) {
  headings <- filter_nav_entries(toc_headings, "section_id", section_id)
  heading_html <- if (length(headings)) {
    paste0(vapply(headings, function(heading) toc_heading_node(heading, toc_refs), character(1)), collapse = "")
  } else {
    ""
  }
  body_html <- paste0(heading_html, child_html)
  if (!nzchar(body_html)) body_html <- "<p class=\"nav-empty\">No subsections available.</p>"
  paste0(
    "<details class=\"nav-branch nav-section\"", if (isTRUE(open)) " open" else "", ">",
    "<summary><a href=\"#", html_escape(section_id), "\" data-target=\"", html_escape(section_id), "\">",
    html_escape(display_o2(section_label)),
    "</a></summary>",
    body_html,
    "</details>"
  )
}

build_toc_html <- function(section_items, toc_headings, toc_refs) {
  if (!"parent_id" %in% names(section_items)) section_items$parent_id <- ""
  section_items$parent_id[is.na(section_items$parent_id)] <- ""
  render_item <- function(i, open = FALSE) {
    section_id <- section_items$id[[i]]
    child_idx <- which(section_items$parent_id == section_id)
    child_html <- if (length(child_idx)) {
      paste0(vapply(child_idx, function(j) render_item(j, open = FALSE), character(1)), collapse = "")
    } else {
      ""
    }
    toc_section_node(
      section_items$id[[i]],
      section_items$label[[i]],
      toc_headings,
      toc_refs,
      child_html = child_html,
      open = open
    )
  }
  top_idx <- which(!nzchar(section_items$parent_id))
  paste0(
    "<div class=\"toc-tree\">",
    paste0(vapply(seq_along(top_idx), function(k) render_item(top_idx[[k]], open = k == 1L), character(1)), collapse = ""),
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

build_reference_section <- function(ctx, mode_dir, output_html, o2, section_number, embed_assets = TRUE, top_n = 10L) {
  section_id <- paste0("reference-o2-", o2_slug(o2))
  ctx$current_section_id <- section_id
  ctx$current_section_number <- section_number
  ctx$current_heading_id <- html_id(paste0("heading-", section_number, ".6-Reference-specific figures"))
  paths <- reference_paths(mode_dir, o2)
  status <- read_tsv_if_exists(paths$status)
  features <- read_csv_if_exists(paths$feature_importance)
  auc <- read_csv_if_exists(paths$auc_by_feature)
  joint <- read_csv_if_exists(paths$top3_joint)
  cumulative <- read_csv_if_exists(paths$cumulative_csv)
  rules <- read_csv_if_exists(paths$decision_rules)

  top_table <- top_feature_table(features, auc, n = top_n)
  joint_table <- if (nrow(joint)) {
    data.frame(
      top1 = joint$top1_feature,
      top2 = joint$top2_feature,
      top3 = joint$top3_feature,
      apparent_AUC = suppressWarnings(as.numeric(joint$apparent_auc)),
      CV_AUC = suppressWarnings(as.numeric(joint$cv_auc)),
      CV_success_folds = joint$cv_success_folds,
      fit_status = joint$fit_status,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    data.frame()
  }

  top30_card <- image_card(
    ctx,
    paths$top30_auc_png,
    paste0("Top30 retained-feature AUC at reference O2=", fmt_o2(o2)),
    "Top-ranked main and pairwise interaction features by elastic-net/stability ranking.",
    top_feature_conclusion(top_table, o2),
    output_html,
    embed_assets
  )
  top_feature_roc_cards <- vapply(
    seq_len(3L),
    function(rank) {
      image_card(
        ctx,
        paths$top_feature_roc_png[[rank]],
        paste0("Top", rank, " retained-feature ROC at reference O2=", fmt_o2(o2)),
        "Single-feature ROC for one of the top three elastic-net/stability-retained features; the plotted AUC is shown on the curve.",
        top_feature_roc_conclusion(top_table, o2, rank),
        output_html,
        embed_assets
      )
    },
    character(1)
  )
  joint_card <- image_card(
    ctx,
    paths$top3_roc_png,
    paste0("Top3 joint ROC at reference O2=", fmt_o2(o2)),
    "Apparent and cross-validated ROC curves for a logistic model using the top three retained features jointly.",
    joint_table_conclusion(joint_table, o2),
    output_html,
    embed_assets
  )
  cumulative_card <- image_card(
    ctx,
    paths$cumulative_png,
    paste0("Cumulative feature AUC at reference O2=", fmt_o2(o2)),
    "AUC after adding ranked features cumulatively from Top1 to Top30.",
    cumulative_conclusion(cumulative, o2),
    output_html,
    embed_assets
  )

  body <- paste0(
    "<p class=\"lead\">Reference O\u2082=", html_escape(fmt_o2(o2)), " uses the fixed-O\u2082 mode label at this O\u2082 concentration as the supervised response.</p>",
    numbered_heading(ctx, 4L, paste0(section_number, ".1"), "Elastic-net fit status"),
    table_block(
      ctx,
      status,
      paste0("Elastic-net fit status for reference O2=", fmt_o2(o2)),
      "Model diagnostics for the penalized logistic regression fit at this reference O2.",
      status_conclusion(status, o2),
      max_rows = max(1L, nrow(status))
    ),
    numbered_heading(ctx, 4L, paste0(section_number, ".2"), "Top retained features"),
    table_block(
      ctx,
      top_table,
      paste0("Top retained features for reference O2=", fmt_o2(o2)),
      "Elastic-net/stability-ranked features with coefficient direction, selection frequency, and single-feature AUC.",
      top_feature_conclusion(top_table, o2),
      max_rows = top_n
    ),
    numbered_heading(ctx, 4L, paste0(section_number, ".3"), "Top3 joint model"),
    table_block(
      ctx,
      joint_table,
      paste0("Top3 joint classifier for reference O2=", fmt_o2(o2)),
      "Logistic model performance when the top three retained features are used jointly.",
      joint_table_conclusion(joint_table, o2),
      max_rows = 1L
    ),
    numbered_heading(ctx, 4L, paste0(section_number, ".4"), "Cumulative AUC checkpoints"),
    table_block(
      ctx,
      cumulative_highlight_table(cumulative),
      paste0("Cumulative AUC checkpoints for reference O2=", fmt_o2(o2)),
      "Selected checkpoints showing how cross-validated AUC changes as more ranked features are added.",
      cumulative_conclusion(cumulative, o2),
      max_rows = 10L
    ),
    numbered_heading(ctx, 4L, paste0(section_number, ".5"), "Presentation rule tree"),
    "<p>The shallow rule tree is not used to rank features; it is a compact presentation layer fitted on main effects only.</p>",
    table_block(
      ctx,
      rules,
      paste0("Presentation rule tree for reference O2=", fmt_o2(o2)),
      "Shallow decision-tree rules summarizing mode-associated parameter regimes.",
      rules_conclusion(rules, o2),
      max_rows = 10L
    ),
    numbered_heading(ctx, 4L, paste0(section_number, ".6"), "Reference-specific figures"),
    figure_grid(top30_card, columns = 1L),
    figure_grid(top_feature_roc_cards, columns = 3L),
    figure_grid(c(joint_card, cumulative_card), columns = 2L)
  )
  subsection(section_id, section_number, paste0("Reference O2=", fmt_o2(o2)), body)
}

build_report_html <- function(mode_dir, output_html, embed_assets = TRUE, top_n = 10L) {
  ctx <- new_report_context()
  mode_dir <- normalizePath(path.expand(mode_dir), mustWork = FALSE)
  index <- read_csv_if_exists(file.path(mode_dir, "mode_parameter_contribution_index.csv"))
  if (!nrow(index)) stop("Missing or empty mode_parameter_contribution_index.csv under: ", mode_dir)
  index <- index[order(suppressWarnings(as.numeric(index$mode_reference_o2))), , drop = FALSE]
  reference_o2 <- suppressWarnings(as.numeric(index$mode_reference_o2))

  run_args <- read_tsv_if_exists(file.path(mode_dir, "mode_parameter_contribution_run_arguments.tsv"))
  top_features <- read_csv_if_exists(file.path(mode_dir, "mode_parameter_top_features_across_reference_o2.csv"))
  top3_auc <- read_csv_if_exists(file.path(mode_dir, "summary_tables", "mode_parameter_top3_auc_by_reference_o2.csv"))
  auc_summary <- read_csv_if_exists(file.path(mode_dir, "summary_tables", "o2_group_auc_summary.csv"))
  joint_summary <- read_csv_if_exists(file.path(mode_dir, "summary_tables", "top3_joint_auc_across_reference_o2.csv"))
  cumulative_summary <- read_csv_if_exists(file.path(mode_dir, "summary_tables", "cumulative_feature_auc_across_reference_o2.csv"))

  reference_numbers <- paste0("3.", seq_along(reference_o2) + 1L)
  main_section_nav_items <- c(
    "1. Methods" = "method",
    "2. Global Results" = "global-results",
    "3. Reference O2 Sections" = "reference-sections",
    "4. Output Files And Provenance" = "outputs"
  )
  reference_nav_items <- stats::setNames(
    paste0("reference-o2-", vapply(reference_o2, o2_slug, character(1))),
    paste0(reference_numbers, " O2=", fmt_o2(reference_o2))
  )
  section_items <- rbind(
    data.frame(
      label = names(main_section_nav_items),
      id = unname(main_section_nav_items),
      parent_id = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    data.frame(
      label = names(reference_nav_items),
      id = unname(reference_nav_items),
      parent_id = "reference-sections",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
  fixed_o2_figs <- fixed_o2_distribution_paths(mode_dir)

  ctx$current_section_id <- "method"
  ctx$current_section_number <- "1"
  ctx$current_heading_id <- ""
  method <- section(
    "method",
    "1",
    "Methods",
    paste0(
      numbered_heading(ctx, 3L, "1.1", "Analysis design"),
      "<p>The response variable was the fixed-O\u2082 attractor mode label for each seed at a chosen reference O\u2082 concentration. ",
      "Mode1 was encoded as 1 and Mode2 as 0. Mode1 means dominant mean ploidy at the reference O\u2082 is at least 2; Mode2 means it is below 2. ",
      "A separate supervised model was fit for each reference O\u2082 value, so feature importance is reference-O\u2082 specific.</p>",
      "<ol>",
      "<li><strong>Input rows.</strong> Each row is one best-fit in vivo seed with a valid Mode1/Mode2 label at the reference O\u2082.</li>",
      "<li><strong>Mode labels.</strong> Labels come from the fixed-O\u2082 attractor dominant mean ploidy rule: <code>dominant_mean_ploidy &gt;= 2</code> gives Mode1, and <code>dominant_mean_ploidy &lt; 2</code> gives Mode2.</li>",
      "<li><strong>Main effects.</strong> The 18 fitted parameters were transformed using the same log10 transform set used by the UMAP workflow where appropriate, then standardized to z-scores within the analysis matrix.</li>",
      "<li><strong>Pairwise interactions.</strong> All pairwise products of standardized main effects were added. With 18 main effects this gives 153 interactions and 171 candidate features total.</li>",
      "<li><strong>Elastic-net classifier.</strong> A binomial logistic regression model was fit with glmnet. The fitted probability is Pr(Mode1), and coefficients indicate whether larger feature values push the model toward Mode1 or Mode2 after accounting for other selected features.</li>",
      "<li><strong>Penalty.</strong> The elastic-net penalty combines L1 sparsity and L2 shrinkage: lambda * [alpha * |beta|_1 + (1-alpha) * |beta|_2^2 / 2]. The current run used alpha from the run arguments, and lambda.1se from cross-validation for a conservative sparse model.</li>",
      "<li><strong>Stability selection layer.</strong> The same lambda.1se was reused across stratified subsamples. Selection frequency is the fraction of successful bootstrap/subsample fits where a feature had nonzero coefficient.</li>",
      "<li><strong>Ranking.</strong> Features were ranked primarily by selection frequency, then by absolute elastic-net coefficient, then by fallback contribution scores from single-feature or interaction tests.</li>",
      "<li><strong>Interpretability outputs.</strong> The report also includes retained-feature AUC bars, Top1-Top3 single-feature ROC curves, top3 joint ROC/AUC, cumulative AUC curves, 2D phase plots for top interactions, and shallow decision rules as a presentation model.</li>",
      "</ol>",
      numbered_heading(ctx, 3L, "1.2", "Mode definition"),
      "<p>For each selected reference O\u2082 concentration, the workflow computes the fixed-O\u2082 attractor dominant mean ploidy for every seed and converts that value into a binary mode label.</p>",
      table_block(
        ctx,
        mode_definition_table(),
        "Mode1 and Mode2 definitions",
        "Binary mode labels used as the supervised response in the elastic-net logistic regression.",
        mode_definition_conclusion(),
        max_rows = 2L
      ),
      numbered_heading(ctx, 3L, "1.3", "Mode balance by reference O2"),
      table_block(
        ctx,
        mode_balance_table(index),
        "Mode balance by reference O2",
        "Seed counts and Mode1/Mode2 balance for each supervised reference-O2 label set.",
        mode_balance_conclusion(index),
        max_rows = nrow(index)
      ),
      numbered_heading(ctx, 3L, "1.4", "Interpretation notes"),
      "<p><strong>Important interpretation.</strong> A selected interaction such as <code>a:b</code> is the product of standardized transformed parameters. It does not mean a mechanistic multiplication was explicitly simulated; it means the mode boundary is better captured when the two parameters are considered jointly.</p>",
      figure_grid(c(
        image_card(
          ctx,
          fixed_o2_figs$all_seeds_png,
          "Fixed-O2 analytical solution across all seeds",
          "Violin and boxplot distributions of the fixed-O2 analytical solution, measured as dominant mean ploidy, across all 500 seeds at each fixed O2 concentration.",
          "This plot shows the overall attractor dominant-ploidy distribution as fixed O2 changes, before splitting seeds by Mode.",
          output_html,
          embed_assets
        ),
        image_card(
          ctx,
          fixed_o2_figs$by_mode_png,
          "Fixed-O2 analytical solution by Mode",
          "The same dominant-mean-ploidy distributions, grouped by the Mode label assigned at each fixed O2 concentration; Mode1 is blue and Mode2 is orange.",
          "This plot shows how the analytical solution separates the two Mode classes within each fixed-O2 condition.",
          output_html,
          embed_assets
        )
      ), columns = 2L)
    )
  )

  ctx$current_section_id <- "global-results"
  ctx$current_section_number <- "2"
  ctx$current_heading_id <- html_id("heading-2.7-Global figures")
  interp <- paste0("<ul>", paste0("<li>", html_escape(display_o2(index_interpretation(index))), "</li>", collapse = ""), "</ul>")
  global_figs <- figure_grid(c(
    image_card(
      ctx,
      file.path(mode_dir, "summary_plots", "o2_group_top3_feature_auc.png"),
      "Top1, Top2, and Top3 feature AUC",
      "Grouped bars compare the retained Top3 single-feature AUC values at each reference O2.",
      global_top3_figure_conclusion(index),
      output_html,
      embed_assets
    ),
    image_card(
      ctx,
      file.path(mode_dir, "summary_plots", "cumulative_feature_auc_summary.png"),
      "Cumulative cross-validated AUC",
      "All reference O2 curves are drawn in one square panel, with color indicating O2 concentration.",
      global_cumulative_figure_conclusion(cumulative_summary),
      output_html,
      embed_assets
    )
  ), columns = 2L)

  ctx$current_section_id <- "global-results"
  ctx$current_section_number <- "2"
  ctx$current_heading_id <- ""
  recurrence <- feature_recurrence_table(top_features, top_n = 10L)
  global_results <- section(
    "global-results",
    "2",
    "Global Results",
    paste0(
      numbered_heading(ctx, 3L, "2.1", "Executive interpretation"), interp,
      numbered_heading(ctx, 3L, "2.2", "Top3 feature summary by reference O2"),
      table_block(
        ctx,
        top3_summary_table(index),
        "Top3 feature summary by reference O2",
        "Per-reference summary of mode balance, retained Top1-Top3 features, single-feature AUCs, and joint Top3 AUCs.",
        top3_summary_conclusion(index),
        max_rows = nrow(index)
      ),
      numbered_heading(ctx, 3L, "2.3", "Most recurrent features among Top10 retained features"),
      table_block(
        ctx,
        recurrence,
        "Most recurrent Top10 retained features",
        "Features repeatedly retained among the Top10 across reference-O2 analyses.",
        feature_recurrence_conclusion(recurrence),
        max_rows = 30L
      ),
      numbered_heading(ctx, 3L, "2.4", "AUC summary table"),
      table_block(
        ctx,
        auc_summary,
        "AUC summary across reference O2",
        "Comparison of best-feature AUC, Top3 mean single-feature AUC, and Top3 joint cross-validated AUC.",
        auc_summary_conclusion(auc_summary),
        max_rows = 30L
      ),
      numbered_heading(ctx, 3L, "2.5", "Top3 joint AUC table"),
      table_block(
        ctx,
        joint_summary,
        "Top3 joint AUC across reference O2",
        "Cross-validated and apparent AUC values for the joint classifier using the retained Top3 features.",
        joint_summary_conclusion(joint_summary),
        max_rows = 30L
      ),
      numbered_heading(ctx, 3L, "2.6", "Top3 single-feature table"),
      table_block(
        ctx,
        top3_auc,
        "Top3 single-feature AUC across reference O2",
        "Long-form table of Top1, Top2, and Top3 retained single-feature AUCs.",
        top3_auc_conclusion(top3_auc),
        max_rows = 30L
      ),
      numbered_heading(ctx, 3L, "2.7", "Global figures"),
      global_figs
    )
  )

  ctx$current_section_id <- "reference-sections"
  ctx$current_section_number <- "3"
  ctx$current_heading_id <- ""
  reference_sections <- paste0(
    section(
      "reference-sections",
      "3",
      "Reference O2 Sections",
      paste0(
        numbered_heading(ctx, 3L, "3.1", "Per-reference analysis layout"),
        "<p>Each subsection below reports the model fit, top retained features, joint Top3 classifier, cumulative AUC checkpoints, rule-list summary, and core figures for one reference O\u2082 value.</p>"
      )
    ),
    paste0(vapply(seq_along(reference_o2), function(i) {
      build_reference_section(ctx, mode_dir, output_html, reference_o2[[i]], reference_numbers[[i]], embed_assets, top_n)
    }, character(1)), collapse = "")
  )

  ctx$current_section_id <- "outputs"
  ctx$current_section_number <- "4"
  ctx$current_heading_id <- ""
  outputs <- section(
    "outputs",
    "4",
    "Output Files And Provenance",
    paste0(
      numbered_heading(ctx, 3L, "4.1", "Report metadata"),
      table_block(
        ctx,
        kv_data(
          c("Mode contribution directory", "Generated at", "Reference O2 values", "Number of reference analyses"),
          c(mode_dir, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), paste(fmt_o2(reference_o2), collapse = ", "), nrow(index))
        ),
        "Report metadata",
        "Location, generation time, and reference-O2 scope for this HTML report.",
        report_metadata_conclusion(index, reference_o2),
        max_rows = 4L,
        class = "report-table report-table--kv"
      ),
      numbered_heading(ctx, 3L, "4.2", "Run arguments"),
      table_block(
        ctx,
        run_args,
        "Run arguments",
        "Command-line settings recorded by the mode-contribution workflow.",
        run_args_conclusion(run_args),
        max_rows = 50L
      ),
      numbered_heading(ctx, 3L, "4.3", "Output file list"),
      "<p>The HTML was generated from existing mode-contribution output files. It does not refit models and does not rerun UMAP.</p>",
      table_block(
        ctx,
        data.frame(
          file = c(
            "mode_parameter_contribution_index.csv",
            "mode_parameter_contribution_run_arguments.tsv",
            "mode_parameter_top_features_across_reference_o2.csv",
            "summary_tables/mode_parameter_top3_auc_by_reference_o2.csv",
            "summary_tables/o2_group_auc_summary.csv",
            "summary_tables/top3_joint_auc_across_reference_o2.csv",
            "summary_tables/cumulative_feature_auc_across_reference_o2.csv",
            "../FixO2Modes/Figures/fixed_o2_analytical_solution_all_seeds_violin_boxplot.png",
            "../FixO2Modes/Figures/fixed_o2_analytical_solution_by_mode_violin_boxplot.png",
            "reference_o2_*/top3_feature_roc/top[1-3]_feature_roc.png",
            "summary_plots/o2_group_top3_feature_auc.png",
            "summary_plots/cumulative_feature_auc_summary.png"
          ),
          purpose = c(
            "Global per-reference summary and Top3 retained features",
            "Command-line settings recorded by the mode-contribution workflow",
            "Top ranked features across all reference O2 values",
            "Long-form Top1/Top2/Top3 feature AUC table",
            "Best-feature, Top3 mean, and Top3 joint AUC comparison",
            "Joint model AUC from the top three retained features",
            "Cumulative AUC after adding Top1 through Top30 features",
            "Overall fixed-O2 analytical-solution distribution across all seeds",
            "Fixed-O2 analytical-solution distribution split by local Mode label",
            "Per-reference Top1, Top2, and Top3 retained-feature ROC curves",
            "Grouped bar plot of Top1, Top2, and Top3 single-feature AUC values",
            "Global cumulative CV AUC figure"
          ),
          stringsAsFactors = FALSE
        ),
        "Output files used by this report",
        "Primary CSV, TSV, and PNG artifacts consumed by the HTML generator.",
        output_files_conclusion(),
        max_rows = 20L
      )
    )
  )

  nav <- build_sidebar_nav(section_items, ctx$toc_headings, ctx$toc_refs, ctx$table_nav, ctx$figure_nav)

  css <- paste0(
    "html{scroll-behavior:smooth;}body{margin:0;background:#f4f7fa;color:#1d2a35;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;}",
    ".shell{display:flex;gap:26px;max-width:1760px;margin:0 auto;padding:22px;}.sidebar{position:sticky;top:22px;align-self:flex-start;width:300px;max-height:calc(100vh - 44px);overflow:auto;border:1px solid #d8e0e8;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.07);}",
    ".side-head{padding:16px;background:#26394d;color:#fff;}.kicker{font-size:11px;letter-spacing:.08em;text-transform:uppercase;opacity:.78;font-weight:700;}.side-title{font-size:18px;font-weight:700;line-height:1.2;margin-top:4px;}.nav-group{border-bottom:1px solid #e4ebf2;}.nav-group summary{cursor:pointer;list-style:none;padding:10px 14px;font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#334e68;background:#f7fafc;}.nav-group summary::-webkit-details-marker,.nav-branch summary::-webkit-details-marker{display:none;}.nav-group summary:before,.nav-branch summary:before{content:'+';display:inline-block;width:16px;color:#6b7c8f;}.nav-group[open]>summary:before,.nav-branch[open]>summary:before{content:'-';}.toc-tree{padding:8px 8px 12px;}.nav-branch{margin:3px 0;border-radius:7px;}.nav-branch summary{cursor:pointer;list-style:none;padding:7px 8px;font-size:12px;line-height:1.3;color:#21384d;background:#fff;border-radius:7px;}.nav-branch summary:hover{background:#eef4fa;}.nav-branch summary a{text-decoration:none;color:inherit;}.nav-section{border:1px solid #edf2f7;background:#fbfdff;}.nav-heading{margin-left:10px;background:#fff;}.nav-heading summary{font-size:11px;color:#334e68;}.nav-list{list-style:none;margin:0;padding:8px 10px 10px;}.nav-sublist{max-height:52vh;overflow:auto;}.nav-tertiary{padding:4px 4px 8px 22px;}.nav-list li{margin:2px 0;}.nav-list a{display:block;padding:7px 10px;border-radius:7px;text-decoration:none;color:#21384d;font-size:12px;line-height:1.3;overflow:hidden;text-overflow:ellipsis;}.nav-tertiary a{font-size:11px;padding:5px 8px;color:#334e68;}.nav-list a:hover{background:#eef4fa;}.nav-list a.active,.nav-branch summary a.active{background:#dcecf8;color:#102a43;font-weight:700;box-shadow:inset 3px 0 0 #2f80c0;}.nav-branch summary a.active{display:inline-block;border-radius:6px;padding:3px 6px;}.nav-empty{font-size:12px;color:#6b7c8f;margin:8px 14px 12px;}",
    ".main{flex:1;min-width:0;max-width:1280px;}.report-section,.report-subsection{margin-bottom:22px;padding:20px;border:1px solid #d8e0e8;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.045);scroll-margin-top:22px;}.report-subsection{padding-top:16px;}.table-block,.figure-card{scroll-margin-top:22px;}h1{font-size:30px;margin:0 0 8px;}h2{font-size:24px;margin:0 0 12px;color:#172a3b;}h3{font-size:19px;margin:18px 0 8px;color:#243b52;}h4{font-size:15px;margin:16px 0 6px;color:#334e68;}.lead{font-size:15px;line-height:1.55;color:#334e68;}p,li{font-size:14px;line-height:1.55;}code{background:#eef2f7;border-radius:4px;padding:1px 4px;}",
    ".table-block{margin:12px 0 22px;}.table-caption{font-size:12px;line-height:1.45;color:#334e68;margin:8px 0 6px;}.table-conclusion{font-size:13px;line-height:1.5;color:#1d2a35;margin:6px 0 0;}.report-table{width:100%;border-collapse:collapse;font-size:12px;margin:10px 0 8px;}.report-table th,.report-table td{padding:7px 8px;border-bottom:1px solid #e4ebf2;text-align:left;vertical-align:top;}.report-table th{background:#f7fafc;font-weight:700;color:#243b52;position:sticky;top:0;}.report-table--kv td:first-child{font-weight:700;color:#334e68;width:230px;}.report-empty{font-style:italic;color:#6b7c8f;}",
    ".figure-caption{font-size:12px;line-height:1.45;color:#334e68;margin:8px 0 4px;}.figure-interpretation{font-size:13px;line-height:1.5;color:#1d2a35;margin:6px 0 0;}",
    ".figure-grid{display:grid;gap:18px;margin:16px 0 22px;}.figure-grid--1{grid-template-columns:1fr;}.figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}.figure-card{min-width:0;}.figure-frame{border:1px solid #d8e0e8;border-radius:8px;background:#fff;padding:8px;display:flex;align-items:center;justify-content:center;overflow:hidden;}.figure-frame img{max-width:100%;height:auto;display:block;}.missing{color:#8a4b00;background:#fff8eb;}",
    "@media(max-width:1100px){.shell{display:block}.sidebar{position:static;width:auto;margin-bottom:18px}.figure-grid--2,.figure-grid--3{grid-template-columns:1fr;}}"
  )

  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Mode Parameter Contribution Report</title><style>", css, "</style></head><body>",
    "<div class=\"shell\"><aside class=\"sidebar\"><div class=\"side-head\"><div class=\"kicker\">parameter landscape</div><div class=\"side-title\">Mode Contribution Report</div></div><nav>", nav, "</nav></aside>",
    "<main class=\"main\"><section class=\"report-section\"><h1>Mode Parameter Contribution Report</h1>",
    "<p class=\"lead\">Elastic-net logistic regression with main effects and pairwise interactions, using fixed-O\u2082 attractor modes as supervised labels.</p></section>",
    method, global_results, reference_sections, outputs,
    "</main></div>",
    report_nav_script(),
    "</body></html>"
  )
}

write_mode_parameter_contribution_report <- function(mode_contribution_dir,
                                                     output_html = file.path(mode_contribution_dir, "mode_parameter_contribution_report.html"),
                                                     embed_assets = TRUE,
                                                     top_n = 10L) {
  mode_contribution_dir <- normalizePath(path.expand(mode_contribution_dir), mustWork = FALSE)
  output_html <- normalizePath(path.expand(output_html), mustWork = FALSE)
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  html <- build_report_html(
    mode_dir = mode_contribution_dir,
    output_html = output_html,
    embed_assets = embed_assets,
    top_n = top_n
  )
  writeLines(html, con = output_html, useBytes = TRUE)
  message("Wrote mode parameter contribution report: ", output_html)
  invisible(output_html)
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  result_root <- normalizePath(
    path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()),
    mustWork = FALSE
  )
  mode_contribution_dir <- normalizePath(
    path.expand(argv$mode_contribution_dir %||% file.path(result_root, "mode_parameter_contribution")),
    mustWork = FALSE
  )
  output_html <- argv$output_html %||% file.path(mode_contribution_dir, "mode_parameter_contribution_report.html")
  embed_assets <- as_bool(argv$embed_assets, TRUE)
  top_n <- as_int(argv$top_n, 10L)
  if (!is.finite(top_n) || is.na(top_n) || top_n < 1L) stop("top_n must be a positive integer.")
  write_mode_parameter_contribution_report(
    mode_contribution_dir = mode_contribution_dir,
    output_html = output_html,
    embed_assets = embed_assets,
    top_n = top_n
  )
}

main()
