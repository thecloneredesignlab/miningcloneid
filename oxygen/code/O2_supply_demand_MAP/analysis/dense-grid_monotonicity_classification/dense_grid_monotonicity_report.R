#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
WORKFLOW_DIR <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)
REPO_ROOT <- normalizePath(file.path(WORKFLOW_DIR, "..", "..", ".."), mustWork = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[[1]], fixed = TRUE)
    val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = TRUE) {
  if (is.null(x) || !length(x) || is.na(x)) return(default)
  val <- tolower(as.character(x[[1]]))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

as_int <- function(x, default = 30L) {
  val <- suppressWarnings(as.integer(x[[1]] %||% default))
  if (length(val) && is.finite(val)) val else default
}

default_result_dir <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "analysis", "monotonicity_classification", "dense-grid_monotonicity_classification")
}

read_tsv_if_exists <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

html_escape <- function(x) {
  x <- as.character(x %||% "")
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

html_id <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("(^-|-$)", "", x)
  if (!nzchar(x)) "section" else x
}

display_o2 <- function(x) {
  gsub("O2", paste0("O", intToUtf8(0x2082)), as.character(x), fixed = TRUE)
}

fmt_num <- function(x, digits = 4L) {
  x <- suppressWarnings(as.numeric(x))
  vapply(x, function(v) {
    if (!is.finite(v)) return("")
    if (abs(v - round(v)) < 1e-9) return(format(round(v), scientific = FALSE, trim = TRUE))
    trimws(formatC(v, digits = digits, format = "fg"))
  }, character(1))
}

table_to_html <- function(df, max_rows = 30L, class = "report-table") {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    return("<p class=\"report-empty\">No rows available.</p>")
  }
  shown_note <- if (nrow(df) > max_rows) sprintf(" Showing first %d of %d rows.", max_rows, nrow(df)) else ""
  df <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  vals <- as.data.frame(lapply(df, function(x) {
    if (is.numeric(x)) fmt_num(x) else display_o2(ifelse(is.na(x), "", as.character(x)))
  }), stringsAsFactors = FALSE, check.names = FALSE)
  names(vals) <- display_o2(names(vals))
  header <- paste0("<th>", html_escape(names(vals)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(vals)), function(i) {
    cells <- paste0("<td>", html_escape(unlist(vals[i, , drop = TRUE], use.names = FALSE)), "</td>", collapse = "")
    paste0("<tr>", cells, "</tr>")
  }, character(1))
  paste0(
    "<table class=\"", class, "\" data-shown-note=\"", html_escape(shown_note), "\"><thead><tr>",
    header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>"
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

kv_data <- function(keys, values) {
  data.frame(item = keys, value = values, check.names = FALSE, stringsAsFactors = FALSE)
}

table_block <- function(ctx, df, title, caption, conclusion = "", max_rows = 30L, class = "report-table") {
  ctx$table_counter <- ctx$table_counter + 1L
  table_number <- next_section_item_number(ctx, ctx$table_section_counts)
  table_id <- paste0("table-", ctx$table_counter)
  shown_note <- if (is.data.frame(df) && nrow(df) > max_rows) sprintf(" The HTML table displays the first %d of %d rows.", max_rows, nrow(df)) else ""
  ctx$table_nav[[length(ctx$table_nav) + 1L]] <- list(id = table_id, label = paste0("Table ", table_number, ". ", title))
  ctx$toc_refs[[length(ctx$toc_refs) + 1L]] <- list(
    id = table_id,
    label = paste0("Table ", table_number, ". ", title),
    section_id = ctx$current_section_id,
    heading_id = ctx$current_heading_id,
    type = "table"
  )
  paste0(
    "<div class=\"table-block\" id=\"", html_escape(table_id), "\" data-nav-target=\"", html_escape(table_id), "\">",
    "<p class=\"table-caption\"><strong>", html_escape(display_o2(paste0("Table ", table_number, ". ", title, "."))), "</strong> ",
    html_escape(display_o2(paste0(caption, shown_note))), "</p>",
    table_to_html(df, max_rows = max_rows, class = class),
    if (nzchar(conclusion)) paste0("<p class=\"table-conclusion\">", html_escape(display_o2(conclusion)), "</p>") else "",
    "</div>"
  )
}

file_to_data_uri <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (!requireNamespace("base64enc", quietly = TRUE)) return(NA_character_)
  ext <- tolower(tools::file_ext(path))
  mime <- if (identical(ext, "jpg") || identical(ext, "jpeg")) "image/jpeg" else "image/png"
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
  ctx$figure_nav[[length(ctx$figure_nav) + 1L]] <- list(id = figure_id, label = paste0("Figure ", figure_number, ". ", title))
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
  interpretation_html <- if (nzchar(conclusion)) paste0("<p class=\"figure-interpretation\">", html_escape(display_o2(conclusion)), "</p>") else ""
  if (!file.exists(path)) {
    return(paste0(
      "<article class=\"figure-card missing\" id=\"", html_escape(figure_id), "\" data-nav-target=\"", html_escape(figure_id), "\">",
      caption_html, "<p>Missing image: ", html_escape(path), "</p>", interpretation_html, "</article>"
    ))
  }
  src <- asset_src(path, output_html = output_html, embed_assets = embed_assets)
  paste0(
    "<article class=\"figure-card\" id=\"", html_escape(figure_id), "\" data-nav-target=\"", html_escape(figure_id), "\">",
    "<div class=\"figure-frame\"><img src=\"", html_escape(src), "\" alt=\"", html_escape(display_o2(figure_title)), "\"></div>",
    caption_html, interpretation_html, "</article>"
  )
}

figure_grid <- function(cards, columns = 2L) {
  paste0("<div class=\"figure-grid figure-grid--", columns, "\">", paste(cards, collapse = ""), "</div>")
}

section <- function(id, number, title, body) {
  paste0(
    "<section class=\"report-section\" id=\"", html_escape(id), "\" data-nav-target=\"", html_escape(id), "\">",
    "<h2>", html_escape(number), ". ", html_escape(display_o2(title)), "</h2>", body, "</section>"
  )
}

subsection <- function(id, number, title, body) {
  paste0(
    "<section class=\"report-subsection\" id=\"", html_escape(id), "\" data-nav-target=\"", html_escape(id), "\">",
    "<h3>", html_escape(number), " ", html_escape(display_o2(title)), "</h3>", body, "</section>"
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
    html_escape(number), " ", html_escape(display_o2(title)), "</h", level, ">"
  )
}

nav_link <- function(id, label, active = TRUE) {
  paste0(
    "<li><a href=\"#", html_escape(id), "\"",
    if (isTRUE(active)) paste0(" data-target=\"", html_escape(id), "\"") else paste0(" data-quick-target=\"", html_escape(id), "\""),
    ">", html_escape(display_o2(label)), "</a></li>"
  )
}

nav_entries_to_html <- function(entries, empty_text, active = TRUE) {
  if (!length(entries)) return(paste0("<p class=\"nav-empty\">", html_escape(empty_text), "</p>"))
  paste0("<ul class=\"nav-list nav-sublist\">", paste0(vapply(entries, function(entry) nav_link(entry$id, entry$label, active = active), character(1)), collapse = ""), "</ul>")
}

nav_group <- function(title, entries_html, count = NA_integer_, open = FALSE) {
  label <- if (is.finite(suppressWarnings(as.numeric(count)))) paste0(title, " (", count, ")") else title
  paste0(
    "<details class=\"nav-group\"", if (isTRUE(open)) " open" else "", ">",
    "<summary>", html_escape(display_o2(label)), "</summary>", entries_html, "</details>"
  )
}

filter_nav_entries <- function(entries, field, value) {
  entries[vapply(entries, function(entry) identical(entry[[field]], value), logical(1))]
}

toc_ref_list <- function(refs) {
  if (!length(refs)) return("<p class=\"nav-empty\">No table or figure references in this subsection.</p>")
  paste0("<ul class=\"nav-list nav-tertiary\">", paste0(vapply(refs, function(ref) {
    nav_link(ref$id, paste0(if (identical(ref$type, "table")) "Table: " else "Figure: ", ref$label), active = TRUE)
  }, character(1)), collapse = ""), "</ul>")
}

toc_heading_node <- function(heading, toc_refs) {
  refs <- filter_nav_entries(toc_refs, "heading_id", heading$id)
  paste0(
    "<details class=\"nav-branch nav-heading\"><summary><a href=\"#", html_escape(heading$id),
    "\" data-target=\"", html_escape(heading$id), "\">", html_escape(display_o2(heading$label)),
    "</a></summary>", toc_ref_list(refs), "</details>"
  )
}

toc_section_node <- function(section_id, section_label, toc_headings, toc_refs, child_html = "", open = FALSE) {
  headings <- filter_nav_entries(toc_headings, "section_id", section_id)
  heading_html <- if (length(headings)) paste0(vapply(headings, function(heading) toc_heading_node(heading, toc_refs), character(1)), collapse = "") else ""
  body_html <- paste0(heading_html, child_html)
  if (!nzchar(body_html)) body_html <- "<p class=\"nav-empty\">No subsections available.</p>"
  paste0(
    "<details class=\"nav-branch nav-section\"", if (isTRUE(open)) " open" else "", ">",
    "<summary><a href=\"#", html_escape(section_id), "\" data-target=\"", html_escape(section_id), "\">",
    html_escape(display_o2(section_label)), "</a></summary>", body_html, "</details>"
  )
}

build_toc_html <- function(section_items, toc_headings, toc_refs) {
  if (!"parent_id" %in% names(section_items)) section_items$parent_id <- ""
  section_items$parent_id[is.na(section_items$parent_id)] <- ""
  render_item <- function(i, open = FALSE) {
    section_id <- section_items$id[[i]]
    child_idx <- which(section_items$parent_id == section_id)
    child_html <- if (length(child_idx)) paste0(vapply(child_idx, function(j) render_item(j, open = FALSE), character(1)), collapse = "") else ""
    toc_section_node(section_items$id[[i]], section_items$label[[i]], toc_headings, toc_refs, child_html = child_html, open = open)
  }
  top_idx <- which(!nzchar(section_items$parent_id))
  paste0("<div class=\"toc-tree\">", paste0(vapply(seq_along(top_idx), function(k) render_item(top_idx[[k]], open = k == 1L), character(1)), collapse = ""), "</div>")
}

build_sidebar_nav <- function(section_items, toc_headings, toc_refs, table_entries, figure_entries) {
  top_count <- if ("parent_id" %in% names(section_items)) sum(!nzchar(ifelse(is.na(section_items$parent_id), "", section_items$parent_id))) else nrow(section_items)
  paste0(
    nav_group("Contents", build_toc_html(section_items, toc_headings, toc_refs), count = top_count, open = TRUE),
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
    "function setActive(id){links.forEach(function(a){a.classList.toggle('active',a.getAttribute('data-target')===id);});var active=linkById[id];if(active){collapseInactive(active);var item=active.closest('li,details');if(item&&item.scrollIntoView){item.scrollIntoView({block:'nearest'});}}}",
    "links.forEach(function(a){a.addEventListener('click',function(){setActive(a.getAttribute('data-target'));});});",
    "if('IntersectionObserver' in window){var visible={};var observer=new IntersectionObserver(function(entries){entries.forEach(function(entry){if(entry.isIntersecting){visible[entry.target.id]=entry.boundingClientRect.top;}else{delete visible[entry.target.id];}});var best=null;Object.keys(visible).forEach(function(id){if(best===null||Math.abs(visible[id])<Math.abs(visible[best])){best=id;}});if(best!==null){setActive(best);}},{rootMargin:'-18% 0px -68% 0px',threshold:[0,0.05,0.25]});targets.forEach(function(el){observer.observe(el);});}else if(targets.length){setActive(targets[0].id);}",
    "})();",
    "</script>"
  )
}

report_css <- function() {
  paste0(
    "html{scroll-behavior:smooth;}body{margin:0;background:#f4f7fa;color:#1d2a35;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;}",
    ".shell{display:flex;gap:26px;max-width:1760px;margin:0 auto;padding:22px;}.sidebar{position:sticky;top:22px;align-self:flex-start;width:300px;max-height:calc(100vh - 44px);overflow:auto;border:1px solid #d8e0e8;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.07);}",
    ".side-head{padding:16px;background:#26394d;color:#fff;}.kicker{font-size:11px;letter-spacing:.08em;text-transform:uppercase;opacity:.78;font-weight:700;}.side-title{font-size:18px;font-weight:700;line-height:1.2;margin-top:4px;}.nav-group{border-bottom:1px solid #e4ebf2;}.nav-group summary{cursor:pointer;list-style:none;padding:10px 14px;font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#334e68;background:#f7fafc;}.nav-group summary::-webkit-details-marker,.nav-branch summary::-webkit-details-marker{display:none;}.nav-group summary:before,.nav-branch summary:before{content:'+';display:inline-block;width:16px;color:#6b7c8f;}.nav-group[open]>summary:before,.nav-branch[open]>summary:before{content:'-';}.toc-tree{padding:8px 8px 12px;}.nav-branch{margin:3px 0;border-radius:7px;}.nav-branch summary{cursor:pointer;list-style:none;padding:7px 8px;font-size:12px;line-height:1.3;color:#21384d;background:#fff;border-radius:7px;}.nav-branch summary:hover{background:#eef4fa;}.nav-branch summary a{text-decoration:none;color:inherit;}.nav-section{border:1px solid #edf2f7;background:#fbfdff;}.nav-heading{margin-left:10px;background:#fff;}.nav-heading summary{font-size:11px;color:#334e68;}.nav-list{list-style:none;margin:0;padding:8px 10px 10px;}.nav-sublist{max-height:52vh;overflow:auto;}.nav-tertiary{padding:4px 4px 8px 22px;}.nav-list li{margin:2px 0;}.nav-list a{display:block;padding:7px 10px;border-radius:7px;text-decoration:none;color:#21384d;font-size:12px;line-height:1.3;overflow:hidden;text-overflow:ellipsis;}.nav-tertiary a{font-size:11px;padding:5px 8px;color:#334e68;}.nav-list a:hover{background:#eef4fa;}.nav-list a.active,.nav-branch summary a.active{background:#dcecf8;color:#102a43;font-weight:700;box-shadow:inset 3px 0 0 #2f80c0;}.nav-branch summary a.active{display:inline-block;border-radius:6px;padding:3px 6px;}.nav-empty{font-size:12px;color:#6b7c8f;margin:8px 14px 12px;}",
    ".main{flex:1;min-width:0;max-width:1280px;}.report-section,.report-subsection{margin-bottom:22px;padding:20px;border:1px solid #d8e0e8;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.045);scroll-margin-top:22px;}.report-subsection{padding-top:16px;}.table-block,.figure-card{scroll-margin-top:22px;}h1{font-size:30px;margin:0 0 8px;}h2{font-size:24px;margin:0 0 12px;color:#172a3b;}h3{font-size:19px;margin:18px 0 8px;color:#243b52;}h4{font-size:15px;margin:16px 0 6px;color:#334e68;}.lead{font-size:15px;line-height:1.55;color:#334e68;}p,li{font-size:14px;line-height:1.55;}code{background:#eef2f7;border-radius:4px;padding:1px 4px;}",
    ".metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px;margin:14px 0 22px;}.metric-card{border:1px solid #d8e0e8;border-radius:8px;background:#fbfdff;padding:12px;}.metric-value{font-size:24px;font-weight:700;color:#172a3b;}.metric-label{font-size:12px;color:#607080;line-height:1.35;}",
    ".table-block{margin:12px 0 22px;}.table-caption{font-size:12px;line-height:1.45;color:#334e68;margin:8px 0 6px;}.table-conclusion{font-size:13px;line-height:1.5;color:#1d2a35;margin:6px 0 0;}.report-table{width:100%;border-collapse:collapse;font-size:12px;margin:10px 0 8px;}.report-table th,.report-table td{padding:7px 8px;border-bottom:1px solid #e4ebf2;text-align:left;vertical-align:top;}.report-table th{background:#f7fafc;font-weight:700;color:#243b52;position:sticky;top:0;}.report-table--kv td:first-child{font-weight:700;color:#334e68;width:230px;}.report-empty{font-style:italic;color:#6b7c8f;}",
    ".figure-caption{font-size:12px;line-height:1.45;color:#334e68;margin:8px 0 4px;}.figure-interpretation{font-size:13px;line-height:1.5;color:#1d2a35;margin:6px 0 0;}.figure-grid{display:grid;gap:18px;margin:16px 0 22px;}.figure-grid--1{grid-template-columns:1fr;}.figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}.figure-card{min-width:0;}.figure-frame{border:1px solid #d8e0e8;border-radius:8px;background:#fff;padding:8px;display:flex;align-items:center;justify-content:center;overflow:hidden;}.figure-frame img{max-width:100%;height:auto;display:block;}.missing{color:#8a4b00;background:#fff8eb;}",
    "@media(max-width:1100px){.shell{display:block}.sidebar{position:static;width:auto;margin-bottom:18px}.figure-grid--2,.figure-grid--3{grid-template-columns:1fr;}}"
  )
}

metric_card <- function(value, label) {
  paste0("<div class=\"metric-card\"><div class=\"metric-value\">", html_escape(value), "</div><div class=\"metric-label\">", html_escape(display_o2(label)), "</div></div>")
}

pretty_file_title <- function(path) {
  x <- tools::file_path_sans_ext(basename(path))
  x <- gsub("^fixed_o2_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(x)
}

figure_caption <- function(path) {
  name <- basename(path)
  if (grepl("all_seed_curves_by_class_all", name)) return("All seeds are overlaid within each curve class; line visibility is adjusted for class density.")
  if (grepl("all_seed_curves_by_class_stable", name)) return("Only stable/reliable seeds are overlaid within each curve class.")
  if (grepl("representative_curves", name)) return("Representative seeds are selected from each curve class to summarize the class shape.")
  if (grepl("heatmap", name)) return("Seeds are ordered by monotonicity class to show class-level structure across the fixed-O2 grid.")
  if (grepl("median_iqr", name)) return("Class-level median and interquartile ranges summarize dominant ploidy across O2.")
  if (grepl("spectral_gap", name)) return("Spectral gap diagnostics flag regions where the dominant eigenmode is weakly separated.")
  if (grepl("all_seed_curves", name)) return("All seed trajectories are colored by the assigned curve class.")
  "Generated figure from the dense-grid fixed-O2 monotonicity workflow."
}

figure_interpretation <- function(path) {
  name <- basename(path)
  if (grepl("all_seed_curves_by_class", name)) return("This panel is intended for visual review of within-class heterogeneity, including whether a nominal class contains dispersed or borderline trajectories.")
  if (grepl("representative_curves", name)) return("Representative curves should be treated as examples, not as the full within-class distribution.")
  if (grepl("heatmap", name)) return("The heatmap gives a compact audit of how dominant ploidy or spectral gap changes across the 201-point analytical O2 grid.")
  if (grepl("min_spectral_gap", name)) return("Seeds with small spectral gaps should be interpreted cautiously because the dominant eigenmode is less well separated.")
  "Use this figure together with the summary tables to audit the classification outcome."
}

expected_figure_paths <- function(figures_dir) {
  list(
    all_seed_class_curves = file.path(figures_dir, "fixed_o2_dominant_ploidy_all_seed_curves_by_class.png"),
    class_median_iqr = file.path(figures_dir, "fixed_o2_dominant_ploidy_median_iqr_by_class.png"),
    dominant_ploidy_heatmap = file.path(figures_dir, "fixed_o2_dominant_ploidy_heatmap_ordered_by_class.png"),
    mode_label_heatmap = file.path(figures_dir, "fixed_o2_mode_label_heatmap_ordered_by_class.png"),
    spectral_gap_heatmap = file.path(figures_dir, "fixed_o2_spectral_gap_heatmap_ordered_by_class.png"),
    spectral_gap_scatter = file.path(figures_dir, "fixed_o2_min_spectral_gap_vs_ploidy_range.png"),
    representative_curves = file.path(figures_dir, "fixed_o2_representative_curves_by_class.png"),
    all_seed_overlay_by_class = file.path(figures_dir, "all_seed_curves_by_class", "fixed_o2_all_seed_curves_by_class_all.png"),
    reliable_seed_overlay_by_class = file.path(figures_dir, "all_seed_curves_by_class", "fixed_o2_all_seed_curves_by_class_stable_reliable.png")
  )
}

figure_definitions <- function(figures_dir) {
  paths <- expected_figure_paths(figures_dir)
  data.frame(
    key = names(paths),
    path = unname(unlist(paths, use.names = FALSE)),
    title = c(
      "All seed dominant-ploidy curves by class",
      "Dominant-ploidy median and IQR by class",
      "Dominant-ploidy heatmap ordered by class",
      "Mode-label heatmap ordered by class",
      "Spectral-gap heatmap ordered by class",
      "Minimum spectral gap versus ploidy range",
      "Representative curves by class",
      "All seed curves overlaid within each class",
      "Stable/reliable seed curves overlaid within each class"
    ),
    caption = c(
      "All seed curves are colored by monotonicity class across the fixed-O2 grid.",
      "Class-level median and interquartile range summarize dominant mean ploidy across O2.",
      "Dominant mean ploidy for every seed, ordered by curve class.",
      "Discrete mode labels at reporting O2 points, ordered by dense-grid curve class.",
      "Spectral gaps across the dense grid, ordered by curve class.",
      "Seed-level reliability diagnostic comparing minimum spectral gap with total ploidy range.",
      "Representative seed curves selected from each curve class.",
      "Every seed within each class is drawn in the same class panel.",
      "Only stable/reliable seeds are drawn within each class panel."
    ),
    conclusion = c(
      "This is the broadest visual audit of class assignment across all seeds.",
      "This summary separates the central class trend from within-class spread.",
      "This heatmap audits whether class labels reflect coherent dominant-ploidy patterns.",
      "This heatmap connects dense-grid curve classes to the fixed-O2 mode labels used elsewhere.",
      "Small-gap bands identify regions where monotonicity calls should be interpreted cautiously.",
      "Seeds with low spectral gaps and large ploidy ranges are the highest-priority reliability checks.",
      "Representative curves are examples only; they do not replace the full overlay figures.",
      "This figure is the main within-class heterogeneity audit for all seeds.",
      "This figure focuses the same audit on stable/reliable seed classifications."
    ),
    stringsAsFactors = FALSE
  )
}

figure_inventory_table <- function(fig_defs, result_dir) {
  rel <- sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", normalizePath(result_dir, mustWork = FALSE)), "/?"), "", normalizePath(fig_defs$path, mustWork = FALSE))
  data.frame(
    figure_key = fig_defs$key,
    file = rel,
    exists = file.exists(fig_defs$path),
    stringsAsFactors = FALSE
  )
}

figure_card_by_key <- function(ctx, fig_defs, key, output_html, embed_assets = TRUE) {
  row <- fig_defs[fig_defs$key == key, , drop = FALSE]
  if (!nrow(row)) stop("Unknown figure key: ", key)
  image_card(
    ctx,
    row$path[[1L]],
    row$title[[1L]],
    row$caption[[1L]],
    row$conclusion[[1L]],
    output_html,
    embed_assets
  )
}

curve_class_definition_table <- function() {
  data.frame(
    curve_class = c(
      "monotone_increasing",
      "single_transition_increase_then_plateau",
      "monotone_decreasing",
      "inverted_u_shaped",
      "u_shaped",
      "complex_nonmonotone"
    ),
    definition = c(
      "Dominant mean ploidy increases across the dense O2 grid within the configured reverse-step tolerance.",
      "Dominant mean ploidy increases and then reaches a terminal plateau under the plateau rule.",
      "Dominant mean ploidy decreases across the dense O2 grid within the configured reverse-step tolerance.",
      "Dominant mean ploidy rises and later falls, producing a single-peaked curve.",
      "Dominant mean ploidy falls and later rises, producing a single-trough curve.",
      "The curve contains shape changes that do not fit the monotone or single-transition classes."
    ),
    stringsAsFactors = FALSE
  )
}

class_summary_table <- function(by_seed, reps) {
  if (!nrow(by_seed)) return(data.frame())
  rows <- lapply(split(by_seed, by_seed$curve_class), function(z) {
    cls <- z$curve_class[[1L]]
    rep_row <- reps[reps$curve_class == cls, , drop = FALSE]
    data.frame(
      curve_class = cls,
      n_seed = length(unique(z$seed_id)),
      reliable = sum(z$monotonicity_reliability == "reliable", na.rm = TRUE),
      caution_small_gap = sum(z$monotonicity_reliability == "caution_small_gap", na.rm = TRUE),
      unreliable_small_gap = sum(z$monotonicity_reliability == "unreliable_small_gap", na.rm = TRUE),
      median_ploidy_range = stats::median(suppressWarnings(as.numeric(z$ploidy_range)), na.rm = TRUE),
      median_net_ploidy_change = stats::median(suppressWarnings(as.numeric(z$net_ploidy_change)), na.rm = TRUE),
      representative_seed = if (nrow(rep_row)) rep_row$representative_seed_id[[1L]] else "",
      representative_reliability = if (nrow(rep_row)) rep_row$representative_reliability[[1L]] else "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$n_seed, out$curve_class), , drop = FALSE]
}

top_seed_examples <- function(by_seed, cls, n = 10L) {
  keep <- by_seed[by_seed$curve_class == cls, , drop = FALSE]
  if (!nrow(keep)) return(data.frame())
  keep <- keep[order(keep$monotonicity_reliability != "reliable", -suppressWarnings(as.numeric(keep$ploidy_range)), keep$seed_number), , drop = FALSE]
  cols <- intersect(
    c("seed_id", "monotonicity_reliability", "ploidy_range", "net_ploidy_change", "min_spectral_gap", "fraction_o2_gap_below_0p005", "n_crossings_ploidy_2", "objective_total"),
    names(keep)
  )
  keep[seq_len(min(nrow(keep), n)), cols, drop = FALSE]
}

build_class_section <- function(ctx, by_seed, reps, output_html, result_dir, cls, number, embed_assets = TRUE) {
  section_id <- paste0("class-", html_id(cls))
  ctx$current_section_id <- section_id
  ctx$current_section_number <- number
  ctx$current_heading_id <- ""
  z <- by_seed[by_seed$curve_class == cls, , drop = FALSE]
  rep_row <- reps[reps$curve_class == cls, , drop = FALSE]
  body <- paste0(
    numbered_heading(ctx, 4L, paste0(number, ".1"), "Class diagnostics"),
    table_block(
      ctx,
      top_seed_examples(by_seed, cls, n = 10L),
      paste0("Example seeds for ", cls),
      "Seeds are ordered to prioritize reliable classifications and larger ploidy ranges.",
      sprintf("The %s class contains %d seeds, including %d reliable seeds.", cls, nrow(z), sum(z$monotonicity_reliability == "reliable", na.rm = TRUE)),
      max_rows = 10L
    ),
    numbered_heading(ctx, 4L, paste0(number, ".2"), "Representative seed"),
    table_block(
      ctx,
      rep_row,
      paste0("Representative seed for ", cls),
      "Representative seed selected from the class-level curve summary.",
      if (nrow(rep_row)) sprintf("Representative seed for %s is %s.", cls, rep_row$representative_seed_id[[1L]]) else sprintf("No representative seed was recorded for %s.", cls),
      max_rows = 5L
    )
  )
  subsection(section_id, number, cls, body)
}

build_report_html <- function(result_dir, output_html, embed_assets = TRUE, top_n = 30L) {
  result_dir <- normalizePath(path.expand(result_dir), mustWork = FALSE)
  output_html <- normalizePath(path.expand(output_html), mustWork = FALSE)
  tables_dir <- file.path(result_dir, "tables")
  figures_dir <- file.path(result_dir, "figures")
  ctx <- new_report_context()

  curves <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_curves.tsv"))
  by_seed <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_by_seed.tsv"))
  class_counts <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_class_counts.tsv"))
  class_by_mode <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_class_by_reporting_o2_mode.tsv"))
  class_curves <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_class_curve_summary.tsv"))
  reps <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_representative_seeds.tsv"))
  objective_rank <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_class_by_objective_rank.tsv"))
  validation <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_validation.tsv"))
  run_args <- read_tsv_if_exists(file.path(tables_dir, "fixed_o2_ploidy_monotonicity_run_arguments.tsv"))

  if (!nrow(by_seed) || !nrow(curves)) {
    stop("Missing dense-grid monotonicity tables under: ", tables_dir)
  }

  reliability <- as.data.frame(table(by_seed$monotonicity_reliability), stringsAsFactors = FALSE)
  names(reliability) <- c("monotonicity_reliability", "n_seed")
  reliability$fraction_seed <- reliability$n_seed / sum(reliability$n_seed)
  reliability <- reliability[order(reliability$monotonicity_reliability), , drop = FALSE]
  passed <- if ("passed" %in% names(validation)) tolower(as.character(validation$passed)) %in% c("true", "1", "yes") else logical()
  validation_summary <- if (nrow(validation)) paste0(sum(passed), "/", nrow(validation)) else "NA"
  classes <- class_counts$curve_class %||% sort(unique(by_seed$curve_class))

  main_section_nav_items <- c(
    "1. Methods" = "method",
    "2. Global Results" = "global-results",
    "3. Curve Class Sections" = "class-sections",
    "4. Output Files And Provenance" = "outputs"
  )
  class_numbers <- paste0("3.", seq_along(classes) + 1L)
  class_nav_items <- stats::setNames(
    paste0("class-", vapply(classes, html_id, character(1))),
    paste0(class_numbers, " ", classes)
  )
  section_items <- rbind(
    data.frame(label = names(main_section_nav_items), id = unname(main_section_nav_items), parent_id = "", stringsAsFactors = FALSE, check.names = FALSE),
    data.frame(label = names(class_nav_items), id = unname(class_nav_items), parent_id = "class-sections", stringsAsFactors = FALSE, check.names = FALSE)
  )

  fig_defs <- figure_definitions(figures_dir)
  fig_inventory <- figure_inventory_table(fig_defs, result_dir)

  ctx$current_section_id <- "method"
  ctx$current_section_number <- "1"
  method <- section(
    "method",
    "1",
    "Methods",
    paste0(
      numbered_heading(ctx, 3L, "1.1", "Analysis design"),
      "<p>Each seed was evaluated on the fixed-O2 analytical grid using the generator dominant eigenvector. The report classifies the dominant mean ploidy curve across O2; this is analytical evaluation, not stochastic simulation.</p>",
      "<ol>",
      "<li><strong>Input.</strong> The analysis uses the 500-seed in vivo best-parameter set and evaluates each seed over the dense fixed-O2 grid.</li>",
      "<li><strong>Grid.</strong> The current output contains ", html_escape(length(unique(curves$O2_pct))), " O2 values and ", html_escape(nrow(curves)), " seed-by-O2 curve rows.</li>",
      "<li><strong>Classification.</strong> Curve classes are assigned from the numerical shape of dominant mean ploidy versus fixed O2, using tolerance settings recorded in the run arguments.</li>",
      "<li><strong>Reliability.</strong> Spectral-gap summaries flag seeds where the dominant eigenmode is weakly separated; these classifications should be interpreted cautiously.</li>",
      "</ol>",
      numbered_heading(ctx, 3L, "1.2", "Curve class definition"),
      table_block(
        ctx,
        curve_class_definition_table(),
        "Curve class definitions",
        "Shape labels used by the dense-grid monotonicity classifier.",
        "These labels summarize the numerical dominant-ploidy curve shape across the fixed-O2 grid.",
        max_rows = 10L
      ),
      numbered_heading(ctx, 3L, "1.3", "Curve class balance by seed"),
      table_block(
        ctx,
        class_counts,
        "Curve class counts",
        "Number and fraction of seeds assigned to each monotonicity class.",
        "This is the primary class-balance table for the dense-grid analysis.",
        max_rows = max(1L, nrow(class_counts))
      ),
      numbered_heading(ctx, 3L, "1.4", "Interpretation notes"),
      "<p><strong>Important interpretation.</strong> The classes describe the dominant-eigenvector analytical solution across fixed O2. A non-monotone class is not automatically a stochastic trajectory result; spectral-gap reliability should be checked before making biological conclusions.</p>",
      figure_grid(c(
        figure_card_by_key(ctx, fig_defs, "all_seed_class_curves", output_html, embed_assets),
        figure_card_by_key(ctx, fig_defs, "mode_label_heatmap", output_html, embed_assets)
      ), columns = 2L)
    )
  )

  ctx$current_section_id <- "global-results"
  ctx$current_section_number <- "2"
  ctx$current_heading_id <- ""
  global_cards <- paste0(
    "<div class=\"metric-grid\">",
    metric_card(length(unique(by_seed$seed_id)), "seeds classified"),
    metric_card(nrow(curves), "dense seed x O2 rows"),
    metric_card(length(unique(curves$O2_pct)), "fixed-O2 grid points"),
    metric_card(validation_summary, "validation checks passed"),
    "</div>"
  )
  global_results <- section(
    "global-results",
    "2",
    "Global Results",
    paste0(
      numbered_heading(ctx, 3L, "2.1", "Executive summary"),
      global_cards,
      table_block(
        ctx,
        reliability,
        "Spectral-gap reliability",
        "Seed counts by spectral-gap reliability category.",
        "Small spectral gaps indicate weaker dominant-eigenmode separation.",
        max_rows = max(1L, nrow(reliability))
      ),
      numbered_heading(ctx, 3L, "2.2", "Curve class summary by class"),
      table_block(
        ctx,
        class_summary_table(by_seed, reps),
        "Class-level diagnostics",
        "Reliability counts, ploidy-range summaries, and representative seeds for each curve class.",
        "This table is useful for deciding which classes are visually stable and which are driven by small-gap or dispersed trajectories.",
        max_rows = 20L
      ),
      numbered_heading(ctx, 3L, "2.3", "Class curve summaries"),
      table_block(
        ctx,
        class_curves,
        "Class curve summary by O2",
        "Class-level median/IQR dominant ploidy and median spectral gap for each O2 grid point.",
        "The HTML shows the first rows; the complete table is available in the result tables directory.",
        max_rows = 30L
      ),
      numbered_heading(ctx, 3L, "2.4", "Reporting-O2 mode crosswalk"),
      table_block(
        ctx,
        class_by_mode,
        "Curve class by reporting-O2 mode",
        "Cross-tabulation of curve classes against reporting O2 mode labels.",
        "This crosswalk connects the dense-curve classes to the discrete mode labels used elsewhere in the parameter-landscape workflow.",
        max_rows = max(30L, nrow(class_by_mode))
      ),
      numbered_heading(ctx, 3L, "2.5", "Objective-rank summary"),
      table_block(
        ctx,
        objective_rank,
        "Curve class by objective rank",
        "Distribution of curve classes across objective-score ranks.",
        "This table checks whether curve classes concentrate among better or worse fitting seeds.",
        max_rows = max(30L, nrow(objective_rank))
      ),
      numbered_heading(ctx, 3L, "2.6", "Global figures"),
      figure_grid(c(
        figure_card_by_key(ctx, fig_defs, "class_median_iqr", output_html, embed_assets),
        figure_card_by_key(ctx, fig_defs, "representative_curves", output_html, embed_assets)
      ), columns = 2L),
      figure_grid(c(
        figure_card_by_key(ctx, fig_defs, "spectral_gap_scatter", output_html, embed_assets)
      ), columns = 1L),
      numbered_heading(ctx, 3L, "2.7", "Heatmap and overlay figures"),
      figure_grid(c(
        figure_card_by_key(ctx, fig_defs, "dominant_ploidy_heatmap", output_html, embed_assets),
        figure_card_by_key(ctx, fig_defs, "spectral_gap_heatmap", output_html, embed_assets)
      ), columns = 2L),
      figure_grid(c(
        figure_card_by_key(ctx, fig_defs, "all_seed_overlay_by_class", output_html, embed_assets),
        figure_card_by_key(ctx, fig_defs, "reliable_seed_overlay_by_class", output_html, embed_assets)
      ), columns = 1L)
    )
  )

  ctx$current_section_id <- "class-sections"
  ctx$current_section_number <- "3"
  ctx$current_heading_id <- ""
  class_sections <- paste0(
    section(
      "class-sections",
      "3",
      "Curve Class Sections",
      paste0(
        numbered_heading(ctx, 3L, "3.1", "Per-class review layout"),
        "<p>Each class subsection lists example seeds and the selected representative seed. The overlay figures in Section 2.4 should be used to inspect the full within-class heterogeneity.</p>"
      )
    ),
    paste0(vapply(seq_along(classes), function(i) {
      build_class_section(ctx, by_seed, reps, output_html, result_dir, classes[[i]], class_numbers[[i]], embed_assets)
    }, character(1)), collapse = "")
  )

  ctx$current_section_id <- "outputs"
  ctx$current_section_number <- "4"
  ctx$current_heading_id <- ""
  table_files <- list.files(tables_dir, pattern = "\\.(tsv|csv)$", full.names = FALSE)
  table_files <- table_files[order(table_files)]
  outputs <- section(
    "outputs",
    "4",
    "Output Files And Provenance",
    paste0(
      numbered_heading(ctx, 3L, "4.1", "Report metadata"),
      table_block(
        ctx,
        kv_data(
          c("Result directory", "Generated at", "Seed count", "O2 grid points", "Curve rows", "Expected figure files", "Available expected figure files"),
          c(result_dir, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), length(unique(by_seed$seed_id)), length(unique(curves$O2_pct)), nrow(curves), nrow(fig_inventory), sum(fig_inventory$exists))
        ),
        "Report metadata",
        "Location, generation time, and dense-grid scope for this HTML report.",
        "The HTML is generated from existing dense-grid monotonicity result files; it does not recompute the analytical grid.",
        max_rows = 10L,
        class = "report-table report-table--kv"
      ),
      numbered_heading(ctx, 3L, "4.2", "Output file list"),
      "<p>The HTML was generated from existing dense-grid monotonicity output files. It does not refit models, rerun UMAP, or rerun fixed-O2 analytical grid evaluation.</p>",
      table_block(
        ctx,
        data.frame(
          file = c(
            "fixed_o2_ploidy_monotonicity_curves.tsv",
            "fixed_o2_ploidy_monotonicity_by_seed.tsv",
            "fixed_o2_ploidy_monotonicity_class_counts.tsv",
            "fixed_o2_ploidy_monotonicity_class_curve_summary.tsv",
            "fixed_o2_ploidy_monotonicity_mode_crosswalk.tsv",
            "fixed_o2_ploidy_monotonicity_representative_seeds.tsv",
            "fixed_o2_ploidy_monotonicity_validation.tsv",
            "figures/fixed_o2_dominant_ploidy_all_seed_curves_by_class.png",
            "figures/fixed_o2_dominant_ploidy_median_iqr_by_class.png",
            "figures/fixed_o2_dominant_ploidy_heatmap_ordered_by_class.png",
            "figures/fixed_o2_mode_label_heatmap_ordered_by_class.png",
            "figures/fixed_o2_spectral_gap_heatmap_ordered_by_class.png",
            "figures/fixed_o2_min_spectral_gap_vs_ploidy_range.png",
            "figures/fixed_o2_representative_curves_by_class.png",
            "figures/all_seed_curves_by_class/fixed_o2_all_seed_curves_by_class_all.png",
            "figures/all_seed_curves_by_class/fixed_o2_all_seed_curves_by_class_stable_reliable.png",
            "report/fixed_o2_ploidy_monotonicity_summary.md"
          ),
          purpose = c(
            "Dense seed-by-O2 analytical curve rows",
            "Seed-level monotonicity class, curve metrics, spectral-gap reliability, and reporting-O2 summaries",
            "Primary class-count summary",
            "Class-level median/IQR curve summaries",
            "Crosswalk between dense-grid class and reporting-O2 mode labels",
            "Representative seed selection for each class",
            "Post-analysis validation checks",
            "All-seed dominant-ploidy curves colored by class",
            "Class median/IQR dominant-ploidy summary",
            "Dominant-ploidy heatmap ordered by class",
            "Reporting-O2 mode-label heatmap ordered by class",
            "Spectral-gap heatmap ordered by class",
            "Reliability diagnostic scatter plot",
            "Representative class curves",
            "All seeds overlaid inside each class panel",
            "Stable/reliable seeds overlaid inside each class panel",
            "Plain markdown summary generated by the monotonicity workflow"
          ),
          stringsAsFactors = FALSE
        ),
        "Output files used by this report",
        "Primary TSV and PNG artifacts consumed by the HTML generator.",
        "These files provide the audit trail for the report tables and embedded figures.",
        max_rows = 20L
      ),
      numbered_heading(ctx, 3L, "4.3", "Expected figure file audit"),
      table_block(
        ctx,
        fig_inventory,
        "Expected figure file audit",
        "Explicit figure-file mapping used by this report.",
        "Every figure embedded in the HTML is selected by this fixed mapping rather than by a directory-wide PNG scan.",
        max_rows = max(1L, nrow(fig_inventory))
      ),
      numbered_heading(ctx, 3L, "4.4", "Available tables"),
      table_block(
        ctx,
        data.frame(file = table_files, path = file.path("tables", table_files), stringsAsFactors = FALSE),
        "Available table files",
        "All TSV/CSV tables found in the result table directory.",
        "These tables are the complete tabular output set available to this report.",
        max_rows = max(1L, top_n)
      )
    )
  )

  nav <- build_sidebar_nav(section_items, ctx$toc_headings, ctx$toc_refs, ctx$table_nav, ctx$figure_nav)
  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Dense-Grid Fixed-O2 Monotonicity Classification Report</title><style>", report_css(), "</style></head><body>",
    "<div class=\"shell\"><aside class=\"sidebar\"><div class=\"side-head\"><div class=\"kicker\">monotonicity classification</div><div class=\"side-title\">Dense-Grid Report</div></div><nav>", nav, "</nav></aside>",
    "<main class=\"main\"><section class=\"report-section\"><h1>Dense-Grid Fixed-O2 Monotonicity Classification Report</h1>",
    "<p class=\"lead\">Analytical fixed-O2 dense-grid classification of dominant mean ploidy curves across seeds, with spectral-gap reliability diagnostics.</p></section>",
    method, global_results, class_sections, outputs,
    "</main></div>",
    report_nav_script(),
    "</body></html>"
  )
}

write_dense_grid_monotonicity_report <- function(result_dir,
                                                 output_html = file.path(result_dir, "report", "index.html"),
                                                 embed_assets = TRUE,
                                                 top_n = 30L) {
  result_dir <- normalizePath(path.expand(result_dir), mustWork = FALSE)
  output_html <- normalizePath(path.expand(output_html), mustWork = FALSE)
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  html <- build_report_html(
    result_dir = result_dir,
    output_html = output_html,
    embed_assets = embed_assets,
    top_n = top_n
  )
  writeLines(html, con = output_html, useBytes = TRUE)
  message("Wrote dense-grid monotonicity report: ", output_html)
  invisible(output_html)
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  result_dir <- normalizePath(path.expand(argv$result_dir %||% argv$out_dir %||% default_result_dir()), mustWork = FALSE)
  output_html <- argv$output_html %||% file.path(result_dir, "report", "index.html")
  embed_assets <- as_bool(argv$embed_assets, TRUE)
  top_n <- as_int(argv$top_n, 30L)
  write_dense_grid_monotonicity_report(
    result_dir = result_dir,
    output_html = output_html,
    embed_assets = embed_assets,
    top_n = top_n
  )
}

if (identical(environment(), globalenv())) {
  main()
}
