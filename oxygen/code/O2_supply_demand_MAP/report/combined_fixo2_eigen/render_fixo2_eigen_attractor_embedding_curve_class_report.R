#!/usr/bin/env Rscript

# Canonical consume-only HTML renderer for the combined FixO2 eigen-attractor workflow.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

WORKFLOW_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
UTIL_DIR <- file.path(WORKFLOW_DIR, "util")
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_path_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_cli_utils.R"))

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

display_o2 <- function(x) {
  gsub("O2", "O<sub>2</sub>", html_escape(x), fixed = TRUE)
}

fmt_num <- function(x, digits = 4L) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return("")
  if (is.numeric(x)) {
    return(format(signif(x, digits), scientific = FALSE, trim = TRUE))
  }
  as.character(x)
}

fmt_pct <- function(x, digits = 1L) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return("")
  paste0(format(round(as.numeric(x) * 100, digits), nsmall = digits, trim = TRUE), "%")
}

read_tsv <- function(path, required = TRUE) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Required table does not exist: ", path, call. = FALSE)
    return(data.frame())
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

split_path <- function(path) {
  strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1L]]
}

relative_path <- function(path, from_dir) {
  path <- normalizePath(path, mustWork = FALSE)
  from_dir <- normalizePath(from_dir, mustWork = FALSE)
  p <- split_path(path)
  f <- split_path(from_dir)
  n <- min(length(p), length(f))
  i <- 1L
  while (i <= n && identical(p[[i]], f[[i]])) i <- i + 1L
  up <- if (i <= length(f)) rep("..", length(f) - i + 1L) else character()
  down <- if (i <= length(p)) p[i:length(p)] else character()
  out <- paste(c(up, down), collapse = "/")
  if (nzchar(out)) out else basename(path)
}

path_href <- function(path, output_dir) {
  html_escape(relative_path(path, output_dir))
}

file_link <- function(path, output_dir, label = NULL) {
  if (is.null(path) || !length(path) || is.na(path[[1L]]) || !nzchar(path[[1L]])) return("")
  label <- label %||% basename(path)
  if (!file.exists(path)) {
    return(paste0("<span class=\"missing\">Missing: ", html_escape(label), "</span>"))
  }
  paste0("<a href=\"", path_href(path, output_dir), "\">", html_escape(label), "</a>")
}

image_data_uri <- function(path) {
  if (is.null(path) || !length(path) || is.na(path[[1L]]) || !nzchar(path[[1L]]) || !file.exists(path)) {
    return("")
  }
  ext <- tolower(tools::file_ext(path))
  mime <- switch(
    ext,
    png = "image/png",
    jpg = "image/jpeg",
    jpeg = "image/jpeg",
    gif = "image/gif",
    svg = "image/svg+xml",
    "application/octet-stream"
  )
  paste0("data:", mime, ";base64,", base64enc::base64encode(path))
}

normalize_id <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  x <- gsub("^-+|-+$", "", x)
  if (nzchar(x)) x else "section"
}

nav_link <- function(id, label, track = TRUE) {
  paste0(
    "<a href=\"#", html_escape(id), "\"",
    if (isTRUE(track)) paste0(" data-target=\"", html_escape(id), "\"") else "",
    ">",
    display_o2(label),
    "</a>"
  )
}

nav_entries_to_html <- function(entries, empty_text, track = FALSE) {
  if (!length(entries)) return(paste0("<p class=\"nav-empty\">", html_escape(empty_text), "</p>"))
  paste0(
    "<ul class=\"nav-list nav-sublist\">",
    paste0(vapply(entries, function(entry) {
      paste0("<li>", nav_link(entry$id, entry$label, track = track), "</li>")
    }, character(1)), collapse = ""),
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
    "<summary>", display_o2(label), "</summary>",
    entries_html,
    "</details>"
  )
}

filter_nav_entries <- function(entries, field, value) {
  entries[vapply(entries, function(entry) identical(entry[[field]], value), logical(1))]
}

toc_section_node <- function(section_id,
                             section_label,
                             child_html = "",
                             open = FALSE,
                             depth = 1L) {
  depth <- max(1L, min(8L, as.integer(depth)))
  child_block <- if (nzchar(child_html)) {
    paste0("<div class=\"nav-children nav-children-depth-", depth, "\">", child_html, "</div>")
  } else {
    ""
  }
  paste0(
    "<details class=\"nav-branch nav-section nav-depth-", depth, "\"",
    if (isTRUE(open)) " open" else "", ">",
    "<summary>", nav_link(section_id, section_label), "</summary>",
    child_block,
    "</details>"
  )
}

build_toc_html <- function(section_items) {
  if (!is.data.frame(section_items) || !nrow(section_items)) {
    return("<p class=\"nav-empty\">No sections available.</p>")
  }
  section_items$parent_id[is.na(section_items$parent_id)] <- ""
  render_item <- function(i, open = FALSE, depth = 1L) {
    section_id <- section_items$id[[i]]
    child_idx <- which(section_items$parent_id == section_id)
    child_html <- if (length(child_idx)) {
      paste0(vapply(child_idx, function(j) {
        render_item(j, open = FALSE, depth = depth + 1L)
      }, character(1)), collapse = "")
    } else {
      ""
    }
    toc_section_node(section_items$id[[i]], section_items$label[[i]], child_html, open, depth)
  }
  top_idx <- which(!nzchar(section_items$parent_id))
  paste0(
    "<div class=\"toc-tree\">",
    paste0(vapply(seq_along(top_idx), function(k) {
      render_item(top_idx[[k]], open = k == 1L, depth = 1L)
    }, character(1)), collapse = ""),
    "</div>"
  )
}

build_sidebar_nav <- function(section_items, table_entries, figure_entries) {
  top_count <- if (nrow(section_items)) sum(!nzchar(section_items$parent_id)) else 0L
  nav_group("Contents", build_toc_html(section_items), count = top_count, open = TRUE)
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
    "function setActive(id){links.forEach(function(a){a.classList.toggle('active',a.getAttribute('data-target')===id);});var active=linkById[id];if(active){collapseInactive(active);var item=active.closest('li,details');if(item&&item.scrollIntoView){item.scrollIntoView({block:'nearest'});}}}",
    "links.forEach(function(a){a.addEventListener('click',function(){setActive(a.getAttribute('data-target'));});});",
    "if('IntersectionObserver' in window){var visible={};var observer=new IntersectionObserver(function(entries){entries.forEach(function(entry){if(entry.isIntersecting){visible[entry.target.id]=entry.boundingClientRect.top;}else{delete visible[entry.target.id];}});var best=null;Object.keys(visible).forEach(function(id){if(best===null||Math.abs(visible[id])<Math.abs(visible[best])){best=id;}});if(best!==null){setActive(best);}},{rootMargin:'-18% 0px -68% 0px',threshold:[0,0.05,0.25]});targets.forEach(function(el){observer.observe(el);});}else if(targets.length){setActive(targets[0].id);}",
    "})();",
    "</script>"
  )
}

report_css <- function() {
  paste0(
    "html{scroll-behavior:smooth;}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f4f7fa;color:#1b2a38;}",
    ".shell{display:flex;gap:26px;max-width:1760px;margin:0 auto;padding:22px;}.sidebar{position:sticky;top:22px;align-self:flex-start;width:320px;max-height:calc(100vh - 44px);overflow:auto;border:1px solid #d8e0e8;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.07);scrollbar-gutter:stable;}",
    ".side-head{padding:16px;background:#26394d;color:#fff;}.kicker{font-size:11px;letter-spacing:.08em;text-transform:uppercase;opacity:.78;font-weight:700;}.side-title{font-size:18px;font-weight:700;line-height:1.2;margin-top:4px;}.side-subtitle{font-size:12px;line-height:1.35;margin-top:4px;opacity:.84;}",
    ".nav-group{border-bottom:1px solid #e4ebf2;}.nav-group summary{cursor:pointer;list-style:none;padding:10px 14px;font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#334e68;background:#f7fafc;}.nav-group summary::-webkit-details-marker,.nav-branch summary::-webkit-details-marker{display:none;}.nav-group summary:before,.nav-branch summary:before{content:'+';display:inline-block;width:16px;color:#6b7c8f;}.nav-group[open]>summary:before,.nav-branch[open]>summary:before{content:'-';}.toc-tree{padding:8px 8px 12px;}.nav-branch{margin:4px 0;border-radius:7px;}.nav-depth-1{margin-left:0;}.nav-depth-2{margin-left:0;}.nav-depth-3{margin-left:0;}.nav-depth-4{margin-left:0;}.nav-depth-5{margin-left:0;}.nav-depth-6,.nav-depth-7,.nav-depth-8{margin-left:0;}.nav-children{margin:4px 0 8px 14px;padding-left:10px;border-left:1px solid #dbe5ee;}.nav-children-depth-2{margin-left:12px;}.nav-children-depth-3,.nav-children-depth-4,.nav-children-depth-5{margin-left:10px;}.nav-branch summary{cursor:pointer;list-style:none;padding:7px 8px;font-size:12px;line-height:1.3;color:#21384d;background:#fff;border-radius:7px;}.nav-depth-1>summary{font-weight:700;}.nav-depth-2>summary{font-weight:650;}.nav-depth-3>summary{font-size:11.5px;}.nav-depth-4>summary,.nav-depth-5>summary,.nav-depth-6>summary{font-size:11px;color:#334e68;}.nav-branch summary:hover{background:#eef4fa;}.nav-branch summary a{text-decoration:none;color:inherit;}.nav-section{border:1px solid #edf2f7;background:#fbfdff;}.nav-depth-2.nav-section,.nav-depth-3.nav-section,.nav-depth-4.nav-section{border-color:#eef3f7;background:#fff;}.nav-list{list-style:none;margin:0;padding:8px 10px 10px;}.nav-sublist{max-height:52vh;overflow:auto;}.nav-list li{margin:2px 0;}.nav-list a{display:block;padding:7px 10px;border-radius:7px;text-decoration:none;color:#21384d;font-size:12px;line-height:1.3;overflow:hidden;text-overflow:ellipsis;}.nav-list a:hover{background:#eef4fa;}.nav-list a.active,.nav-branch summary a.active{background:#dcecf8;color:#102a43;font-weight:700;box-shadow:inset 3px 0 0 #2f80c0;}.nav-branch summary a.active{display:inline-block;border-radius:6px;padding:3px 6px;}.nav-empty{font-size:12px;color:#6b7c8f;margin:8px 14px 12px;}",
    ".main{flex:1;min-width:0;max-width:1280px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}.report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-card h2,.report-section h2{margin-top:0;}.report-card,.report-section,.report-subsection,.table-block,.report-figure-card{scroll-margin-top:24px;}",
    ".report-meta{margin:0 0 6px;color:#516274;font-size:14px;line-height:1.45;}.method-list li{margin:8px 0;line-height:1.55;}.summary-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin:14px 0 18px;}.summary-card{padding:12px 14px;border:1px solid #e1e8f0;border-radius:8px;background:#fbfdff;}.summary-card .label{font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:#627386;font-weight:700;}.summary-card .value{margin-top:4px;font-size:22px;font-weight:750;color:#17324c;}.summary-card .sub{font-size:12px;color:#657789;margin-top:3px;}",
    ".table-block{margin:18px 0 24px;}.table-caption{font-size:12px;line-height:1.45;color:#334e68;margin:8px 0 6px;}.table-source,.figure-links{font-size:12px;color:#536577;margin:8px 0;}.report-table{width:100%;border-collapse:collapse;font-size:13px;margin:10px 0 18px 0;}.report-table th,.report-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}.report-table th{background:#f7f9fb;font-weight:700;}.report-table-wrapper{overflow:auto;border:1px solid #edf2f7;border-radius:8px;}.report-table-wrapper .report-table{margin:0;}.missing{color:#a6422b;font-weight:650;}",
    ".report-subsection{margin:24px 0 30px;padding-top:2px;}.report-subsection h3{margin:0 0 12px;font-size:21px;color:#17324c;border-bottom:1px solid #e4ebf2;padding-bottom:6px;}.report-subsection h4{margin:18px 0 10px;font-size:17px;color:#243b52;}.report-subsection h5{margin:16px 0 8px;font-size:15px;color:#334e68;}.report-figure-grid{display:grid;gap:18px;margin:16px 0 24px 0;align-items:start;}.report-figure-grid--1{grid-template-columns:1fr;}.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}",
    ".report-figure-card{min-width:0;}.report-figure{margin:0 0 10px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;display:flex;align-items:center;justify-content:center;overflow:hidden;}.report-figure--square{aspect-ratio:1/1;}.report-figure--wide{aspect-ratio:2/1;}.report-figure--natural{display:block;aspect-ratio:auto;}.report-figure-image{width:100%;height:100%;display:block;object-fit:contain;border-radius:6px;}.report-figure--natural .report-figure-image{height:auto;object-fit:contain;}.report-figure-title{margin:8px 0 0;font-size:13px;line-height:1.35;font-weight:700;color:#243b52;}.report-figure-legend{margin:6px 0 0;font-size:12px;line-height:1.45;color:#536577;}",
    "@media (max-width: 1100px){.shell{display:block;}.sidebar{position:static;width:auto;margin-bottom:20px;}.summary-grid{grid-template-columns:repeat(2,minmax(0,1fr));}}@media (max-width: 900px){.report-figure-grid--2,.report-figure-grid--3{grid-template-columns:1fr;}.summary-grid{grid-template-columns:1fr;}}"
  )
}

table_to_html <- function(df, id, label, caption, source_path, output_dir, max_rows = Inf) {
  table_df <- df
  note <- ""
  if (is.finite(max_rows) && nrow(df) > max_rows) {
    table_df <- df[seq_len(max_rows), , drop = FALSE]
    note <- paste0(" Showing first ", max_rows, " of ", nrow(df), " rows.")
  }
  if (!nrow(table_df)) {
    body <- "<p class=\"report-meta\">No rows available.</p>"
  } else {
    header <- paste0("<tr>", paste0("<th>", display_o2(names(table_df)), "</th>", collapse = ""), "</tr>")
    rows <- apply(table_df, 1L, function(row) {
      paste0("<tr>", paste0("<td>", display_o2(row), "</td>", collapse = ""), "</tr>")
    })
    body <- paste0(
      "<div class=\"report-table-wrapper\"><table class=\"report-table\">",
      "<thead>", header, "</thead><tbody>",
      paste0(rows, collapse = ""),
      "</tbody></table></div>"
    )
  }
  paste0(
    "<div class=\"table-block\" id=\"", html_escape(id), "\" data-nav-target>",
    "<h4>", display_o2(label), "</h4>",
    "<p class=\"table-caption\"><strong>Description:</strong> ", display_o2(caption), html_escape(note), "</p>",
    body,
    "</div>"
  )
}

figure_card <- function(id,
                        label,
                        png_path,
                        pdf_path = NULL,
                        output_dir,
                        caption = "",
                        wide = FALSE,
                        natural = FALSE) {
  fig_class <- if (isTRUE(natural)) {
    "report-figure report-figure--natural"
  } else if (isTRUE(wide)) {
    "report-figure report-figure--wide"
  } else {
    "report-figure report-figure--square"
  }
  data_uri <- image_data_uri(png_path)
  image_html <- if (nzchar(data_uri)) {
    paste0(
      "<img class=\"report-figure-image\" src=\"", data_uri,
      "\" alt=\"", html_escape(label), "\">"
    )
  } else {
    paste0("<p class=\"missing\">Missing preview: ", html_escape(basename(png_path %||% "")), "</p>")
  }
  paste0(
    "<div class=\"report-figure-card\" id=\"", html_escape(id), "\" data-nav-target>",
    "<figure class=\"", fig_class, "\">", image_html, "</figure>",
    "<p class=\"report-figure-title\">", display_o2(label), "</p>",
    if (nzchar(caption)) paste0("<p class=\"report-figure-legend\"><strong>Description:</strong> ", display_o2(caption), "</p>") else "",
    "</div>"
  )
}

section_html <- function(id, title, body, level = 2L) {
  paste0(
    "<section class=\"report-section\" id=\"", html_escape(id), "\" data-nav-target>",
    "<h", level, ">", display_o2(title), "</h", level, ">",
    body,
    "</section>"
  )
}

subsection_html <- function(id, title, body, level = 3L) {
  paste0(
    "<section class=\"report-subsection\" id=\"", html_escape(id), "\" data-nav-target>",
    "<h", level, ">", display_o2(title), "</h", level, ">",
    body,
    "</section>"
  )
}

arg_value <- function(run_args, key, default = "") {
  if (!nrow(run_args) || !"argument" %in% names(run_args) || !"value" %in% names(run_args)) return(default)
  hit <- which(run_args$argument == key)
  if (length(hit)) as.character(run_args$value[[hit[[1L]]]]) else default
}

classification_method_html <- function(run_args) {
  span <- arg_value(run_args, "smooth_span", "0.2")
  degree <- arg_value(run_args, "smooth_degree", "2")
  family <- arg_value(run_args, "smooth_family", "gaussian")
  rule <- arg_value(run_args, "classification_rule_version", "loess_persistent_v1")
  min_span <- arg_value(run_args, "min_segment_span_fraction", "0.02")
  min_amp_abs <- arg_value(run_args, "min_segment_amplitude_abs", "0.01")
  min_amp_frac <- arg_value(run_args, "min_segment_amplitude_fraction", "0.03")
  min_points <- arg_value(run_args, "min_segment_points", "3")
  flat <- arg_value(run_args, "flat_range_threshold", "0.05")
  terminal <- arg_value(run_args, "terminal_plateau_span_fraction", "0.1")
  paste0(
    "<ol class=\"method-list\">",
    "<li>The input curve for each seed is the fixed-O2 dense-grid dominant-eigenvector summary: dominant mean ploidy evaluated across the fixed O2 grid. This report does not rerun stochastic simulation or model fitting.</li>",
    "<li>For each seed, the raw ploidy-vs-O2 curve is smoothed by <code>stats::loess</code> with span <code>", html_escape(span), "</code>, degree <code>", html_escape(degree), "</code>, and family <code>", html_escape(family), "</code>. If loess cannot return a usable curve, the script falls back to <code>smooth.spline</code>.</li>",
    "<li>The smoothed curve is converted into adjacent fitted differences along the ordered O2 grid. A local step is treated as meaningful only when it exceeds <code>max(step_epsilon_abs, step_epsilon_fraction * fitted_ploidy_range)</code>; otherwise it is treated as flat.</li>",
    "<li>Consecutive positive, negative, or flat local signs are compressed into sign segments. Only persistent segments can define the final class: a segment must span at least <code>", html_escape(min_span), "</code> of the O2 range, contain at least <code>", html_escape(min_points), "</code> grid points, and exceed both the absolute amplitude threshold <code>", html_escape(min_amp_abs), "</code> and the fitted-range-scaled amplitude threshold <code>", html_escape(min_amp_frac), "</code>.</li>",
    "<li>If the fitted ploidy range is below <code>", html_escape(flat), "</code>, the curve is classified as <code>approximately_flat</code>. Otherwise, the persistent sign sequence determines whether the curve is monotone increasing, monotone decreasing, U-shaped, inverted-U-shaped, a single transition into a terminal plateau, or complex nonmonotone. Terminal plateaus are recognized only near the end of the O2 grid using a terminal span fraction of <code>", html_escape(terminal), "</code>.</li>",
    "<li>The final class used downstream is the regression-smoothed class from rule version <code>", html_escape(rule), "</code>. The pointwise class is retained only as an audit field to show how smoothing and persistence filters changed the earlier dense-grid sign classification.</li>",
    "</ol>"
  )
}

summary_cards_html <- function(class_counts, by_seed, global_test) {
  n_seed <- if (nrow(class_counts) && "n_seed" %in% names(class_counts)) {
    sum(class_counts$n_seed, na.rm = TRUE)
  } else if (nrow(by_seed)) {
    nrow(by_seed)
  } else {
    NA_integer_
  }
  top_class <- if (nrow(class_counts)) {
    paste0(class_counts$smooth_curve_class[[1L]], " (", fmt_pct(class_counts$fraction_seed[[1L]]), ")")
  } else {
    ""
  }
  changed <- if (nrow(by_seed) && "class_changed" %in% names(by_seed)) {
    sum(as.logical(by_seed$class_changed), na.rm = TRUE)
  } else {
    NA_integer_
  }
  global_p <- if (nrow(global_test) && "p_value" %in% names(global_test)) {
    fmt_num(global_test$p_value[[1L]], digits = 3L)
  } else {
    "NA"
  }
  paste0(
    "<div class=\"summary-grid\">",
    "<div class=\"summary-card\"><div class=\"label\">Seeds</div><div class=\"value\">", html_escape(fmt_num(n_seed)), "</div><div class=\"sub\">classified best-fit seeds</div></div>",
    "<div class=\"summary-card\"><div class=\"label\">Largest class</div><div class=\"value\">", display_o2(top_class), "</div><div class=\"sub\">by smoothed curve class</div></div>",
    "<div class=\"summary-card\"><div class=\"label\">Changed</div><div class=\"value\">", html_escape(fmt_num(changed)), "</div><div class=\"sub\">pointwise to smoothed class changes</div></div>",
    "<div class=\"summary-card\"><div class=\"label\">Global objective test</div><div class=\"value\">p=", html_escape(global_p), "</div><div class=\"sub\">Kruskal-Wallis across classes</div></div>",
    "</div>"
  )
}

pretty_stub_label <- function(stub) {
  x <- sub("^pooled_invivo_invitro_initial_vs_best_", "", stub)
  x <- sub("^sampled500_", "sampled 500 ", x)
  x <- gsub("common_prior_unit", "common prior-unit", x, fixed = TRUE)
  x <- gsub("context_prior_unit", "context prior-unit", x, fixed = TRUE)
  x <- gsub("pca_umap", "PCA-to-UMAP", x, fixed = TRUE)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("\\bpca\\b", "PCA", x, ignore.case = TRUE)
  x <- gsub("\\bumap\\b", "UMAP", x, ignore.case = TRUE)
  x <- gsub("\\btsne\\b", "t-SNE", x, ignore.case = TRUE)
  x
}

classification_figure_specs <- function(classification_dir) {
  fig_dir <- file.path(classification_dir, "figures")
  data.frame(
    file = c(
      "fixed_o2_regression_smoothed_all_seed_curves_by_class.png",
      "fixed_o2_regression_smoothed_median_iqr_by_class.png",
      "fixed_o2_regression_pointwise_vs_smooth_class_transition.png",
      "fixed_o2_regression_objective_by_curve_class_boxplot.png"
    ),
    label = c(
      "Regression-smoothed all-seed curves by class",
      "Regression-smoothed median and IQR by class",
      "Pointwise vs regression-smoothed class transition",
      "Objective by O2-Ploidy curve class"
    ),
    stringsAsFactors = FALSE
  ) |>
    within({
      png <- file.path(fig_dir, file)
      pdf <- sub("\\.png$", ".pdf", png)
    })
}

classification_figure_description <- function(label) {
  switch(
    label,
    "Regression-smoothed all-seed curves by class" =
      "For each regression-smoothed curve class, all seed-level smoothed dominant-mean-ploidy curves are shown across the fixed-O2 grid, with the class median overlaid.",
    "Regression-smoothed median and IQR by class" =
      "Class-level median smoothed dominant mean ploidy and interquartile range summarize the typical O2-ploidy response for each curve class.",
    "Pointwise vs regression-smoothed class transition" =
      "The heatmap audits how the earlier pointwise dense-grid class assignments map to the final regression-smoothed persistent-segment classes.",
    "Objective by O2-Ploidy curve class" =
      "Objective values for in vivo best-fit seeds are grouped by final regression-smoothed O2-Ploidy curve class, with significance annotations from the objective-by-class tests.",
    "Generated by the regression-smoothed fixed-O2 monotonicity classification workflow."
  )
}

default_classification_dir <- function(repo_root = bpf_repo_root(SCRIPT_DIR)) {
  file.path(bpf_dense_grid_result_root(repo_root), "dense-grid_monotonicity_classification")
}

default_embedding_dir <- function(repo_root = bpf_repo_root(SCRIPT_DIR)) {
  file.path(bpf_combine_fixo2_eigen_attractor_result_dir(repo_root), "pooled_embedding_curve_class")
}

default_output_html <- function(embedding_dir) {
  file.path(embedding_dir, "fixo2_eigen_attractor_embedding_curve_class_report.html")
}

manifest_variants_for_reduction <- function(manifest, reduction) {
  rows <- manifest[manifest$reduction == reduction, , drop = FALSE]
  vals <- unique(as.character(rows$variant))
  vals <- vals[nzchar(vals)]
  preferred <- c("Full", "BestOnly", "Sampled")
  c(intersect(preferred, vals), setdiff(vals, preferred))
}

build_report_html <- function(classification_dir, embedding_dir, output_html) {
  output_dir <- dirname(output_html)
  class_table_dir <- file.path(classification_dir, "tables")
  class_report_dir <- file.path(classification_dir, "report")
  manifest_path <- file.path(embedding_dir, "tables", "pooled_embedding_curve_class_manifest.tsv")

  manifest <- read_tsv(manifest_path)
  run_args_path <- file.path(class_table_dir, "fixed_o2_ploidy_monotonicity_regression_run_arguments.tsv")
  by_seed_path <- file.path(class_table_dir, "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv")
  class_counts_path <- file.path(class_table_dir, "fixed_o2_ploidy_monotonicity_regression_class_counts.tsv")
  change_audit_path <- file.path(class_table_dir, "fixed_o2_ploidy_monotonicity_regression_class_change_audit.tsv")
  objective_summary_path <- file.path(class_table_dir, "fixed_o2_regression_objective_by_curve_class_boxplot_summary.tsv")
  objective_global_path <- file.path(class_table_dir, "fixed_o2_regression_objective_by_curve_class_boxplot_global_test.tsv")
  objective_pairwise_path <- file.path(class_table_dir, "fixed_o2_regression_objective_by_curve_class_boxplot_pairwise_tests.tsv")
  summary_md_path <- file.path(class_report_dir, "fixed_o2_ploidy_monotonicity_regression_summary.md")

  run_args <- read_tsv(run_args_path, required = FALSE)
  by_seed <- read_tsv(by_seed_path)
  class_counts <- read_tsv(class_counts_path)
  change_audit <- read_tsv(change_audit_path, required = FALSE)
  objective_summary <- read_tsv(objective_summary_path, required = FALSE)
  objective_global <- read_tsv(objective_global_path, required = FALSE)
  objective_pairwise <- read_tsv(objective_pairwise_path, required = FALSE)

  section_items <- data.frame(id = character(), label = character(), parent_id = character(), stringsAsFactors = FALSE)
  table_entries <- list()
  figure_entries <- list()

  add_section <- function(id, label, parent_id = "") {
    section_items <<- rbind(
      section_items,
      data.frame(id = id, label = label, parent_id = parent_id, stringsAsFactors = FALSE)
    )
  }
  add_table <- function(id, label) {
    table_entries[[length(table_entries) + 1L]] <<- list(id = id, label = label)
  }
  add_figure <- function(id, label) {
    figure_entries[[length(figure_entries) + 1L]] <<- list(id = id, label = label)
  }
  table_counter <- 0L
  figure_counter <- 0L
  render_table <- function(df, id, label, caption, source_path, max_rows = Inf) {
    table_counter <<- table_counter + 1L
    numbered_label <- paste0("Table ", table_counter, ". ", label)
    add_table(id, numbered_label)
    table_to_html(
      df = df,
      id = id,
      label = numbered_label,
      caption = caption,
      source_path = source_path,
      output_dir = output_dir,
      max_rows = max_rows
    )
  }
  render_figure <- function(id, label, png_path, pdf_path = NULL, caption = "", wide = FALSE, natural = FALSE) {
    figure_counter <<- figure_counter + 1L
    numbered_label <- paste0("Figure ", figure_counter, ". ", label)
    add_figure(id, numbered_label)
    figure_card(
      id = id,
      label = numbered_label,
      png_path = png_path,
      pdf_path = pdf_path,
      output_dir = output_dir,
      caption = caption,
      wide = wide,
      natural = natural
    )
  }

  add_section("classification", "1. Regression classification")
  add_section("classification-method", "1.1 Method", "classification")
  add_section("classification-results", "1.2 Classification results", "classification")
  add_section("classification-objective", "1.3 Objective by class", "classification")
  add_section("classification-figures", "1.4 Classification figures", "classification")
  add_section("scatter-integration", "2. Curve class in FixO2 eigen-attractor scatter")
  reductions <- c("PCAs", "UMAPs", "TSNEs")
  reduction_numbers <- c(PCAs = "2.1", UMAPs = "2.2", TSNEs = "2.3")
  for (reduction in reductions) {
    rid <- paste0("scatter-", normalize_id(reduction))
    rlabel <- switch(reduction, PCAs = "PCA", UMAPs = "UMAP", TSNEs = "t-SNE")
    rnum <- reduction_numbers[[reduction]]
    add_section(rid, paste(rnum, rlabel), "scatter-integration")
    variants <- manifest_variants_for_reduction(manifest, reduction)
    for (j in seq_along(variants)) {
      variant <- variants[[j]]
      add_section(paste0(rid, "-", normalize_id(variant)), paste0(rnum, ".", j, " ", variant), rid)
    }
  }

  class_counts_disp <- class_counts
  if (nrow(class_counts_disp) && "fraction_seed" %in% names(class_counts_disp)) {
    class_counts_disp$fraction_seed <- vapply(class_counts_disp$fraction_seed, fmt_pct, character(1))
  }
  change_disp <- change_audit
  if (nrow(change_disp)) {
    frac_cols <- grep("^fraction_", names(change_disp), value = TRUE)
    for (col in frac_cols) change_disp[[col]] <- vapply(change_disp[[col]], fmt_pct, character(1))
  }
  significant_pairwise <- objective_pairwise
  if (nrow(significant_pairwise) && "significant" %in% names(significant_pairwise)) {
    significant_pairwise <- significant_pairwise[as.logical(significant_pairwise$significant), , drop = FALSE]
  }

  classification_results <- paste0(
    summary_cards_html(class_counts, by_seed, objective_global),
    render_table(class_counts_disp, "table-class-counts", "Regression-smoothed class counts",
                 "Seed counts and fractions for the final smoothed curve class used in the scatter mapping.",
                 class_counts_path),
    render_table(change_disp, "table-class-change-audit", "Pointwise vs smoothed class-change audit",
                 "Audit of how the persistent regression-smoothed rule changes the earlier pointwise class assignment.",
                 change_audit_path, max_rows = 20L),
    ""
  )

  objective_results <- paste0(
    render_table(objective_summary, "table-objective-summary", "Objective summary by curve class",
                 "Distribution summary for the objective values grouped by regression-smoothed curve class.",
                 objective_summary_path),
    render_table(objective_global, "table-objective-global", "Global objective test",
                 "Kruskal-Wallis test across all regression-smoothed curve classes.",
                 objective_global_path),
    render_table(significant_pairwise, "table-objective-pairwise", "Significant pairwise objective tests",
                 "Pairwise tests after multiple-testing adjustment; only significant rows are shown.",
                 objective_pairwise_path, max_rows = 50L),
    render_table(run_args, "table-run-arguments", "Classification run arguments",
                 "Arguments and derived thresholds recorded by the regression-smoothed classification script.",
                 run_args_path)
  )

  class_fig_specs <- classification_figure_specs(classification_dir)
  class_fig_specs <- class_fig_specs[file.exists(class_fig_specs$png), , drop = FALSE]
  classification_figures <- if (nrow(class_fig_specs)) {
    full_width_labels <- c(
      "Regression-smoothed all-seed curves by class",
      "Regression-smoothed median and IQR by class"
    )
    full_width_idx <- which(class_fig_specs$label %in% full_width_labels)
    full_width_html <- if (length(full_width_idx)) {
      paste0(vapply(full_width_idx, function(i) {
      id <- paste0("fig-classification-", normalize_id(class_fig_specs$label[[i]]))
        paste0(
        "<div class=\"report-figure-grid report-figure-grid--1\">",
        render_figure(
          id,
          class_fig_specs$label[[i]],
          class_fig_specs$png[[i]],
          class_fig_specs$pdf[[i]],
          caption = classification_figure_description(class_fig_specs$label[[i]]),
          wide = TRUE,
          natural = identical(class_fig_specs$label[[i]], "Regression-smoothed all-seed curves by class")
        ),
        "</div>"
      )
      }, character(1)), collapse = "")
    } else {
      ""
    }
    other_idx <- setdiff(seq_len(nrow(class_fig_specs)), full_width_idx)
    other_html <- if (length(other_idx)) {
      paste0(
        "<div class=\"report-figure-grid report-figure-grid--2\">",
        paste0(vapply(other_idx, function(i) {
          id <- paste0("fig-classification-", normalize_id(class_fig_specs$label[[i]]))
          render_figure(
            id,
            class_fig_specs$label[[i]],
            class_fig_specs$png[[i]],
            class_fig_specs$pdf[[i]],
            caption = classification_figure_description(class_fig_specs$label[[i]]),
            wide = TRUE
          )
        }, character(1)), collapse = ""),
        "</div>"
      )
    } else {
      ""
    }
    paste0(full_width_html, other_html)
  } else {
    "<p class=\"report-meta\">No classification figures were found.</p>"
  }

  scatter_intro <- paste0(
    "<p class=\"report-meta\">The pooled FixO2 eigen-attractor coordinate tables are reused directly. ",
    "Each embedding entry first shows a 2x2 summary figure with the original pooled in vivo/in vitro scatter, ",
    "the in vivo best-fit O2-Ploidy curve-class scatter, the per-class 3x3 highlight panel, and the average-slope 3x3 highlight panel. ",
    "A standalone monotone-increasing detail figure is shown immediately after each 2x2 figure.</p>"
  )

  scatter_sections <- c(scatter_intro)
  for (reduction in reductions) {
    rid <- paste0("scatter-", normalize_id(reduction))
    rlabel <- switch(reduction, PCAs = "PCA", UMAPs = "UMAP", TSNEs = "t-SNE")
    rnum <- reduction_numbers[[reduction]]
    variant_blocks <- character()
    variants <- manifest_variants_for_reduction(manifest, reduction)
    for (j in seq_along(variants)) {
      variant <- variants[[j]]
      vid <- paste0(rid, "-", normalize_id(variant))
      rows <- manifest[manifest$reduction == reduction & manifest$variant == variant, , drop = FALSE]
      if (nrow(rows)) {
        fig_cards <- paste0(vapply(seq_len(nrow(rows)), function(i) {
          row <- rows[i, , drop = FALSE]
          label <- pretty_stub_label(row$stub[[1L]])
          id <- paste0("fig-", normalize_id(reduction), "-", normalize_id(variant), "-", normalize_id(row$stub[[1L]]))
          main_card <- render_figure(
            id,
            paste(rlabel, variant, label),
            row$extended_combined_png[[1L]],
            row$extended_combined_pdf[[1L]],
            caption = paste0(
              "2x2 pooled FixO2 eigen-attractor summary for ", rlabel, " ", variant, " data using the ", label,
              " preprocessing. The top-left panel shows the original pooled in vivo/in vitro scatter; the top-right panel maps the regression-smoothed O2-Ploidy curve class onto in vivo best-fit points; ",
              "the bottom-left panel shows per-class curve-class highlights; and the bottom-right panel shows the same per-class layout colored by average slope. ",
              "Rows: ", row$n_rows[[1L]], "; initial points: ", row$n_initial[[1L]],
              "; in vivo best fits with class: ", row$n_invivo_best_with_class[[1L]],
              "; missing class: ", row$n_invivo_best_missing_class[[1L]], "."
            ),
            natural = TRUE
          )
          monotone_card <- ""
          if (
            all(c("monotone_increasing_combined_png", "monotone_increasing_combined_pdf") %in% names(row)) &&
              !is.na(row$monotone_increasing_combined_png[[1L]]) &&
              nzchar(row$monotone_increasing_combined_png[[1L]])
          ) {
            monotone_id <- paste0(id, "-monotone-increasing")
            monotone_card <- render_figure(
              monotone_id,
              paste(rlabel, variant, label, "monotone increasing detail"),
              row$monotone_increasing_combined_png[[1L]],
              row$monotone_increasing_combined_pdf[[1L]],
              caption = paste0(
                "Standalone monotone increasing detail plot for ", rlabel, " ", variant, " data using the ", label,
                " preprocessing. The left panel highlights only in vivo best-fit points classified as monotone increasing; ",
                "the right panel colors those same points by average slope using the monotone increasing class-specific scale. ",
                "Monotone increasing points: ", row$monotone_increasing_n[[1L]],
                "; class-specific average-slope range: ",
                fmt_num(row$monotone_increasing_average_slope_min[[1L]]), " to ",
                fmt_num(row$monotone_increasing_average_slope_max[[1L]]), "."
              ),
              wide = TRUE
            )
          }
          paste0(main_card, monotone_card)
        }, character(1)), collapse = "")
        body <- paste0("<div class=\"report-figure-grid report-figure-grid--1\">", fig_cards, "</div>")
      } else {
        body <- "<p class=\"report-meta\">No figures available for this group.</p>"
      }
      variant_blocks <- c(variant_blocks, subsection_html(vid, paste0(rnum, ".", j, " ", variant), body, level = 4L))
    }
    scatter_sections <- c(scatter_sections, subsection_html(rid, paste(rnum, rlabel), paste0(variant_blocks, collapse = ""), level = 3L))
  }

  classification_body <- paste0(
    subsection_html("classification-method", "1.1 Regression classification method", classification_method_html(run_args), level = 3L),
    subsection_html("classification-results", "1.2 Regression classification results", classification_results, level = 3L),
    subsection_html("classification-objective", "1.3 Objective values by O2-Ploidy curve class", objective_results, level = 3L),
    subsection_html("classification-figures", "1.4 Regression classification figures", classification_figures, level = 3L)
  )

  report_body <- paste0(
    section_html("classification", "1. Regression classification", classification_body, level = 2L),
    section_html("scatter-integration", "2. Curve class mapped onto pooled FixO2 eigen-attractor scatter plots", paste0(scatter_sections, collapse = ""), level = 2L)
  )

  nav <- build_sidebar_nav(section_items, table_entries, figure_entries)
  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>O2-Ploidy Curve Class in FixO2 Eigen-Attractor Scatter</title>",
    "<style>", report_css(), "</style></head><body>",
    "<div class=\"shell\">",
    "<aside class=\"sidebar\"><div class=\"side-head\"><div class=\"kicker\">best-fit parameter feature</div>",
    "<div class=\"side-title\">O<sub>2</sub>-Ploidy Curve Class</div>",
    "<div class=\"side-subtitle\">classification and pooled FixO2 eigen-attractor scatter</div></div><nav>",
    nav,
    "</nav></aside><main class=\"main\">",
    "<section class=\"report-card\" id=\"report-metadata\" data-nav-target>",
    "<h1>O<sub>2</sub>-Ploidy Curve Class in FixO2 Eigen-Attractor Scatter</h1>",
    "<p class=\"report-meta\">This report combines regression-smoothed fixed-O2 curve classification with pooled in vivo/in vitro PCA, UMAP, and t-SNE FixO2 eigen-attractor scatter plots.</p>",
    "<p class=\"report-meta\"><strong>Classification root:</strong> ", html_escape(normalizePath(classification_dir, mustWork = FALSE)), "<br>",
    "<strong>Embedding root:</strong> ", html_escape(normalizePath(embedding_dir, mustWork = FALSE)), "<br>",
    "<strong>Generated at:</strong> ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p>",
    "</section>",
    report_body,
    "</main></div>",
    report_nav_script(),
    "</body></html>"
  )
}

validate_report_inputs <- function(classification_dir, embedding_dir) {
  required <- c(
    file.path(classification_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"),
    file.path(classification_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_class_counts.tsv"),
    file.path(embedding_dir, "tables", "pooled_embedding_curve_class_manifest.tsv")
  )
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Missing required report inputs:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

write_pooled_embedding_curve_class_report <- function(classification_dir = default_classification_dir(),
                                                      embedding_dir = default_embedding_dir(),
                                                      output_html = default_output_html(embedding_dir),
                                                      dry_run = FALSE) {
  classification_dir <- normalizePath(path.expand(classification_dir), mustWork = FALSE)
  embedding_dir <- normalizePath(path.expand(embedding_dir), mustWork = FALSE)
  output_html <- normalizePath(path.expand(output_html), mustWork = FALSE)
  validate_report_inputs(classification_dir, embedding_dir)
  message("Classification root: ", classification_dir)
  message("Embedding root: ", embedding_dir)
  message("Output HTML: ", output_html)
  if (isTRUE(dry_run)) return(invisible(output_html))
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  html <- build_report_html(classification_dir, embedding_dir, output_html)
  writeLines(html, con = output_html, useBytes = TRUE)
  if (!file.exists(output_html) || file.info(output_html)$size <= 0) {
    stop("Report renderer did not create a non-empty HTML file: ", output_html, call. = FALSE)
  }
  message("Wrote report: ", output_html)
  invisible(output_html)
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  args <- bpf_parse_args(raw_args)
  repo_root <- bpf_repo_root(SCRIPT_DIR)
  classification_dir <- bpf_resolve_repo_path(args$classification_dir %||% default_classification_dir(repo_root), repo_root)
  embedding_dir <- bpf_resolve_repo_path(args$embedding_dir %||% default_embedding_dir(repo_root), repo_root)
  output_html <- bpf_resolve_repo_path(args$output_html %||% default_output_html(embedding_dir), repo_root)
  dry_run <- bpf_as_bool(args$dry_run, FALSE)
  write_pooled_embedding_curve_class_report(
    classification_dir = classification_dir,
    embedding_dir = embedding_dir,
    output_html = output_html,
    dry_run = dry_run
  )
}

if (identical(environment(), globalenv())) {
  main()
}
