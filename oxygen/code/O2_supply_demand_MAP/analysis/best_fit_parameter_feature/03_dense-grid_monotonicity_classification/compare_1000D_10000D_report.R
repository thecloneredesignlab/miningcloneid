#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  normalizePath(sys.frames()[[1L]]$ofile %||% "", mustWork = FALSE)
}

script_dir <- function() {
  p <- script_path()
  if (nzchar(p)) dirname(p) else getwd()
}

repo_root <- function() {
  cur <- normalizePath(script_dir(), mustWork = FALSE)
  for (i in seq_len(12L)) {
    if (dir.exists(file.path(cur, "oxygen", "code", "O2_supply_demand_MAP"))) {
      return(normalizePath(cur, mustWork = FALSE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not locate repository root from: ", script_dir(), call. = FALSE)
}

default_compare_root <- function() {
  file.path(
    repo_root(), "oxygen", "results", "analysis", "best_fit_parameter_feature",
    "03_dense-grid_monotonicity_classification", "monotonicity_classification", "compare"
  )
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    if (!grepl("=", arg, fixed = TRUE)) {
      out[[arg]] <- TRUE
    } else {
      key <- sub("=.*$", "", arg)
      val <- sub("^[^=]*=", "", arg)
      out[[key]] <- val
    }
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x)) return(default)
  val <- tolower(trimws(as.character(x[[1L]])))
  if (val %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (val %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

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
  x[!nzchar(x)] <- "section"
  x
}

file_to_data_uri <- function(path, mime = "image/png") {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Required R package is not installed: base64enc")
  }
  if (!file.exists(path)) stop("Missing asset: ", path)
  paste0("data:", mime, ";base64,", base64enc::base64encode(path))
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing table: ", path)
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
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
    paste0(
      "<tr>",
      paste0("<td>", html_escape(unlist(vals[i, , drop = TRUE], use.names = FALSE)), "</td>", collapse = ""),
      "</tr>"
    )
  }, character(1))
  paste0("<table class=\"report-table\"><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
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
    html_escape(label),
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
    "<summary>", html_escape(label), "</summary>",
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
    html_escape(heading$label),
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
    html_escape(section_label),
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

report_css <- function() {
  paste0(
    "html{scroll-behavior:smooth;}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f4f7fa;color:#1b2a38;}",
    ".shell{display:flex;gap:26px;max-width:1760px;margin:0 auto;padding:22px;}.sidebar{position:sticky;top:22px;align-self:flex-start;width:320px;max-height:calc(100vh - 44px);overflow:auto;border:1px solid #d8e0e8;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.07);scrollbar-gutter:stable;}",
    ".side-head{padding:16px;background:#26394d;color:#fff;}.kicker{font-size:11px;letter-spacing:.08em;text-transform:uppercase;opacity:.78;font-weight:700;}.side-title{font-size:18px;font-weight:700;line-height:1.2;margin-top:4px;}.side-subtitle{font-size:12px;line-height:1.35;margin-top:4px;opacity:.84;}",
    ".nav-group{border-bottom:1px solid #e4ebf2;}.nav-group summary{cursor:pointer;list-style:none;padding:10px 14px;font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#334e68;background:#f7fafc;}.nav-group summary::-webkit-details-marker,.nav-branch summary::-webkit-details-marker{display:none;}.nav-group summary:before,.nav-branch summary:before{content:'+';display:inline-block;width:16px;color:#6b7c8f;}.nav-group[open]>summary:before,.nav-branch[open]>summary:before{content:'-';}.toc-tree{padding:8px 8px 12px;}.nav-branch{margin:3px 0;border-radius:7px;}.nav-depth-1{margin-left:0;}.nav-depth-2{margin-left:8px;}.nav-depth-3{margin-left:16px;}.nav-depth-4{margin-left:24px;}.nav-depth-5{margin-left:32px;}.nav-depth-6,.nav-depth-7,.nav-depth-8{margin-left:40px;}.nav-branch summary{cursor:pointer;list-style:none;padding:7px 8px;font-size:12px;line-height:1.3;color:#21384d;background:#fff;border-radius:7px;}.nav-depth-1>summary{font-weight:700;}.nav-depth-2>summary{font-weight:650;}.nav-depth-3>summary{font-size:11.5px;}.nav-depth-4>summary,.nav-depth-5>summary,.nav-depth-6>summary{font-size:11px;color:#334e68;}.nav-branch summary:hover{background:#eef4fa;}.nav-branch summary a{text-decoration:none;color:inherit;}.nav-section{border:1px solid #edf2f7;background:#fbfdff;}.nav-heading{margin-left:10px;background:#fff;}.nav-heading summary{font-size:11px;color:#334e68;}.nav-list{list-style:none;margin:0;padding:8px 10px 10px;}.nav-sublist{max-height:52vh;overflow:auto;}.nav-tertiary{padding:4px 4px 8px 22px;}.nav-list li{margin:2px 0;}.nav-list a{display:block;padding:7px 10px;border-radius:7px;text-decoration:none;color:#21384d;font-size:12px;line-height:1.3;overflow:hidden;text-overflow:ellipsis;}.nav-tertiary a{font-size:11px;padding:5px 8px;color:#334e68;}.nav-list a:hover{background:#eef4fa;}.nav-list a.active,.nav-branch summary a.active{background:#dcecf8;color:#102a43;font-weight:700;box-shadow:inset 3px 0 0 #2f80c0;}.nav-branch summary a.active{display:inline-block;border-radius:6px;padding:3px 6px;}.nav-empty{font-size:12px;color:#6b7c8f;margin:8px 14px 12px;}",
    ".main{flex:1;min-width:0;max-width:1280px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}.report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-card h2,.report-section h2{margin-top:0;}.report-card,.report-section,.report-reduction-section,.report-subtype-section,.report-variant-section,.report-subsection,.report-parameter-block{scroll-margin-top:24px;}",
    ".report-meta{margin:0;color:#516274;font-size:14px;}.report-empty{color:#657789;font-style:italic;}.table-block{scroll-margin-top:24px;}.table-caption{font-size:12px;line-height:1.45;color:#334e68;margin:8px 0 6px;}.report-table{width:100%;border-collapse:collapse;font-size:13px;margin:10px 0 18px 0;}.report-table th,.report-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}.report-table th{background:#f7f9fb;font-weight:700;}",
    ".report-parameter-block{margin:8px 0 22px;padding:12px 14px;border:1px solid #e1e8f0;border-radius:8px;background:#fbfdff;}.report-reduction-section{margin:26px 0 30px;padding-top:2px;}.report-reduction-section h3{margin:0 0 12px;font-size:21px;color:#17324c;border-bottom:1px solid #e4ebf2;padding-bottom:6px;}.report-subtype-section{margin:18px 0 24px;}.report-subtype-section h4{margin:0 0 10px;font-size:17px;color:#243b52;}.report-variant-section{margin:18px 0 24px;padding:14px;border:1px solid #edf2f7;border-radius:8px;background:#fff;}.report-variant-section h4,.report-variant-section h5{margin:0 0 6px;font-size:15px;color:#334e68;}.report-subsection{margin:18px 0 30px;padding-top:2px;}.report-subsection h5{margin:0 0 8px;font-size:15px;color:#17324c;}.report-subsection h6{margin:0 0 8px;font-size:14px;color:#17324c;}ul{margin-top:8px;}.report-figure-grid{display:grid;gap:18px;margin:16px 0 24px 0;align-items:stretch;}.report-figure-grid--1{grid-template-columns:1fr;}.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}",
    ".report-table-grid{display:grid;gap:18px;margin:10px 0 18px 0;align-items:start;}.report-table-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-table-grid h5{margin:0 0 6px 0;font-size:13px;color:#284662;}",
    ".report-figure-card{min-width:0;scroll-margin-top:24px;}.report-figure{margin:0 0 10px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;aspect-ratio:1/1;display:flex;align-items:center;justify-content:center;overflow:hidden;}.report-figure--wide{aspect-ratio:auto;min-height:0;}.report-figure-image{width:100%;height:100%;display:block;object-fit:contain;border-radius:6px;}.report-figure--wide .report-figure-image{height:auto;max-height:none;}.report-figure-title,.report-figure-legend{margin:8px 0 0 0;}.report-figure-title{font-size:13px;line-height:1.35;}.report-figure-legend{font-size:12px;line-height:1.4;color:#536577;}.report-links{font-size:12px;margin-top:6px;}.report-links a{color:#1b5f93;text-decoration:none;}.report-links a:hover{text-decoration:underline;}",
    "@media (max-width: 1100px){.shell{display:block;}.sidebar{position:static;width:auto;margin-bottom:20px;}}@media (max-width: 900px){.report-figure-grid--2,.report-figure-grid--3{grid-template-columns:1fr;}}"
  )
}

report_paths <- function(compare_root, label) {
  terminal_day <- as.integer(sub("D$", "", label))
  base <- file.path(compare_root, label)
  list(
    label = label,
    terminal_day = terminal_day,
    base = base,
    tables = file.path(base, "tables"),
    figures = file.path(base, "figures"),
    validation = file.path(base, "tables", "validation.tsv"),
    run_args = file.path(base, "tables", "analysis_run_arguments.tsv"),
    consistency_summary = file.path(base, "tables", paste0("day", terminal_day, "_curve_class_consistency_summary.tsv")),
    gap_png = file.path(base, "figures", paste0("day", terminal_day, "_abs_delta_vs_eigen_spectral_gap.png")),
    gap_pdf = file.path(base, "figures", paste0("day", terminal_day, "_abs_delta_vs_eigen_spectral_gap.pdf")),
    consistency_png = file.path(base, "figures", paste0("day", terminal_day, "_curve_class_consistency_summary.png")),
    consistency_pdf = file.path(base, "figures", paste0("day", terminal_day, "_curve_class_consistency_summary.pdf")),
    seed_curves_png = file.path(base, "figures", "convergence_to_eigen_by_time_seed_curves_panels.png"),
    seed_curves_pdf = file.path(base, "figures", "convergence_to_eigen_by_time_seed_curves_panels.pdf"),
    o2_curves_png = file.path(base, "figures", "convergence_to_eigen_by_time_o2_seed_curves_panels.png"),
    o2_curves_pdf = file.path(base, "figures", "convergence_to_eigen_by_time_o2_seed_curves_panels.pdf"),
    by_class_dir = file.path(base, "figures", "convergence_to_eigen_by_time_o2_seed_curves_by_class")
  )
}

rel_path <- function(path, base_dir) {
  path <- normalizePath(path, mustWork = FALSE)
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  prefix <- paste0(base_dir, .Platform$file.sep)
  if (startsWith(path, prefix)) sub(prefix, "", path, fixed = TRUE) else path
}

asset_src <- function(path, output_dir, embed_assets = TRUE) {
  if (isTRUE(embed_assets)) file_to_data_uri(path, mime = "image/png") else rel_path(path, output_dir)
}

figure_spec <- function(id, label, png, pdf = NULL, legend = "") {
  list(id = id, label = label, png = png, pdf = pdf, legend = legend)
}

figure_card <- function(spec, output_dir, embed_assets = TRUE, wide = TRUE) {
  img <- if (file.exists(spec$png)) {
    paste0(
      "<img class=\"report-figure-image\" src=\"",
      html_escape(asset_src(spec$png, output_dir = output_dir, embed_assets = embed_assets)),
      "\" alt=\"", html_escape(spec$label), "\">"
    )
  } else {
    paste0("<p class=\"report-empty\">Missing image: ", html_escape(spec$png), "</p>")
  }
  pdf_link <- if (!is.null(spec$pdf) && file.exists(spec$pdf)) {
    paste0("<p class=\"report-links\"><a href=\"", html_escape(rel_path(spec$pdf, output_dir)), "\">Open PDF</a></p>")
  } else {
    ""
  }
  paste0(
    "<article class=\"report-figure-card\" id=\"", html_escape(spec$id), "\" data-nav-target=\"", html_escape(spec$id), "\">",
    "<div class=\"report-figure", if (isTRUE(wide)) " report-figure--wide" else "", "\">", img, "</div>",
    "<p class=\"report-figure-title\"><strong>", html_escape(spec$label), "</strong></p>",
    if (nzchar(spec$legend %||% "")) paste0("<p class=\"report-figure-legend\">", html_escape(spec$legend), "</p>") else "",
    pdf_link,
    "</article>"
  )
}

figure_grid <- function(specs, output_dir, embed_assets = TRUE, columns = 1L, wide = TRUE) {
  paste0(
    "<div class=\"report-figure-grid report-figure-grid--", as.integer(columns), "\">",
    paste0(vapply(specs, function(spec) figure_card(spec, output_dir, embed_assets = embed_assets, wide = wide), character(1)), collapse = ""),
    "</div>"
  )
}

table_block <- function(id, label, df, max_rows = 50L) {
  paste0(
    "<div class=\"table-block\" id=\"", html_escape(id), "\">",
    "<p class=\"table-caption\"><strong>", html_escape(label), "</strong></p>",
    table_to_html(df, max_rows = max_rows),
    "</div>"
  )
}

section_block <- function(id, title, body) {
  paste0(
    "<section class=\"report-section\" id=\"", html_escape(id), "\" data-nav-target=\"", html_escape(id), "\">",
    "<h2>", html_escape(title), "</h2>",
    body,
    "</section>"
  )
}

class_key_from_file <- function(path) {
  key <- tools::file_path_sans_ext(basename(path))
  sub("^convergence_to_eigen_by_time_o2_seed_curves_", "", key)
}

class_label <- function(key) {
  tools::toTitleCase(gsub("_", " ", key, fixed = TRUE))
}

common_class_keys <- function(paths_1000, paths_10000) {
  keys_1000 <- vapply(paths_1000, class_key_from_file, character(1))
  keys_10000 <- vapply(paths_10000, class_key_from_file, character(1))
  common <- intersect(keys_1000, keys_10000)
  preferred <- c(
    "monotone_increasing",
    "monotone_decreasing",
    "single_transition_increase_then_plateau",
    "u_shaped",
    "inverted_u_shaped",
    "complex_nonmonotone"
  )
  c(intersect(preferred, common), sort(setdiff(common, preferred)))
}

validation_summary <- function(label, validation_path) {
  tab <- read_tsv(validation_path)
  passed <- if ("passed" %in% names(tab)) {
    tolower(as.character(tab$passed)) %in% c("true", "t", "1", "yes", "pass")
  } else {
    character()
  }
  data.frame(
    analysis = label,
    n_tests = nrow(tab),
    n_pass = sum(passed, na.rm = TRUE),
    all_passed = if (length(passed)) all(passed) else NA,
    stringsAsFactors = FALSE
  )
}

selected_run_args <- function(label, run_arg_path) {
  tab <- read_tsv(run_arg_path)
  keep <- c(
    "analysis_max_day_requested",
    "terminal_day",
    "analysis_label",
    "n_seed",
    "n_o2",
    "selected_days",
    "plot_days",
    "n_comparison_rows",
    "steady_state_interpretation"
  )
  tab <- tab[tab$argument %in% keep, , drop = FALSE]
  data.frame(analysis = label, tab, stringsAsFactors = FALSE, check.names = FALSE)
}

build_report_html <- function(compare_root,
                              output_html = NULL,
                              labels = c("1000D", "10000D"),
                              embed_assets = TRUE,
                              report_title = "Initial-Ploidy Trajectory vs Eigenvector Comparison Report") {
  compare_root <- normalizePath(path.expand(compare_root), mustWork = TRUE)
  output_html <- normalizePath(path.expand(output_html %||% file.path(compare_root, "compare_1000D_vs_10000D_report.html")), mustWork = FALSE)
  output_dir <- dirname(output_html)
  paths <- lapply(labels, function(label) report_paths(compare_root, label))
  names(paths) <- labels
  if (length(paths) != 2L) stop("This report expects exactly two analysis labels.")

  required <- unlist(lapply(paths, function(p) {
    c(
      p$validation,
      p$run_args,
      p$consistency_summary,
      p$gap_png,
      p$gap_pdf,
      p$consistency_png,
      p$consistency_pdf,
      p$seed_curves_png,
      p$seed_curves_pdf,
      p$o2_curves_png,
      p$o2_curves_pdf
    )
  }), use.names = FALSE)
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing required report inputs:\n", paste(missing, collapse = "\n"))

  class_files <- lapply(paths, function(p) {
    list.files(p$by_class_dir, pattern = "\\.png$", full.names = TRUE)
  })
  class_keys <- common_class_keys(class_files[[1L]], class_files[[2L]])
  if (!length(class_keys)) stop("No matching by-class convergence PNG files were found.")

  class_png <- function(p, key) file.path(p$by_class_dir, paste0("convergence_to_eigen_by_time_o2_seed_curves_", key, ".png"))
  class_pdf <- function(p, key) file.path(p$by_class_dir, paste0("convergence_to_eigen_by_time_o2_seed_curves_", key, ".pdf"))
  class_required <- unlist(lapply(paths, function(p) {
    as.vector(rbind(vapply(class_keys, function(key) class_png(p, key), character(1)),
                    vapply(class_keys, function(key) class_pdf(p, key), character(1))))
  }), use.names = FALSE)
  missing_class <- class_required[!file.exists(class_required)]
  if (length(missing_class)) stop("Missing by-class report inputs:\n", paste(missing_class, collapse = "\n"))

  overview_validation <- do.call(rbind, Map(function(label, p) validation_summary(label, p$validation), names(paths), paths))
  overview_args <- do.call(rbind, Map(function(label, p) selected_run_args(label, p$run_args), names(paths), paths))
  consistency_tables <- do.call(rbind, Map(function(label, p) {
    tab <- read_tsv(p$consistency_summary)
    data.frame(analysis = label, tab, stringsAsFactors = FALSE, check.names = FALSE)
  }, names(paths), paths))

  section_items <- data.frame(
    id = c(
      "overview",
      "spectral-gap",
      "curve-class-consistency",
      "seed-level-convergence",
      "o2-level-convergence",
      "by-class-convergence",
      paste0("by-class-", html_id(class_keys))
    ),
    label = c(
      "1 Overview",
      "2 Spectral-gap comparison",
      "3 Curve-class consistency",
      "4 Seed-level convergence",
      "5 O2-level convergence",
      "6 By-class O2 convergence",
      paste0("6.", seq_along(class_keys), " ", class_label(class_keys))
    ),
    parent_id = c(rep("", 6L), rep("by-class-convergence", length(class_keys))),
    stringsAsFactors = FALSE
  )

  table_entries <- list(
    list(id = "table-validation-summary", label = "Validation summary"),
    list(id = "table-run-arguments", label = "Run arguments"),
    list(id = "table-class-consistency", label = "Curve-class consistency summary")
  )

  fig_gap <- Map(function(label, p) {
    figure_spec(
      id = paste0("fig-", tolower(label), "-spectral-gap"),
      label = paste0(label, " day ", p$terminal_day, " absolute delta vs eigen spectral gap"),
      png = p$gap_png,
      pdf = p$gap_pdf,
      legend = "Terminal-day finite-time trajectory differences are compared against the eigenvector dominant-attractor ploidy and spectral gap."
    )
  }, names(paths), paths)

  fig_consistency <- Map(function(label, p) {
    figure_spec(
      id = paste0("fig-", tolower(label), "-class-consistency"),
      label = paste0(label, " day ", p$terminal_day, " curve-class consistency"),
      png = p$consistency_png,
      pdf = p$consistency_pdf,
      legend = "Seed-level curve-class agreement between the eigenvector method and finite-time initial-ploidy trajectories at the terminal day."
    )
  }, names(paths), paths)

  fig_seed <- Map(function(label, p) {
    figure_spec(
      id = paste0("fig-", tolower(label), "-seed-convergence"),
      label = paste0(label, " seed-level convergence to eigenvector result"),
      png = p$seed_curves_png,
      pdf = p$seed_curves_pdf,
      legend = "All seeds are summarized across time by 2N minus eigen, 4N minus eigen, and 4N minus 2N differences."
    )
  }, names(paths), paths)

  fig_o2 <- Map(function(label, p) {
    figure_spec(
      id = paste0("fig-", tolower(label), "-o2-convergence"),
      label = paste0(label, " O2-level convergence to eigenvector result"),
      png = p$o2_curves_png,
      pdf = p$o2_curves_pdf,
      legend = "The same convergence summaries are split by selected fixed-O2 concentrations."
    )
  }, names(paths), paths)

  class_figs <- lapply(class_keys, function(key) {
    Map(function(label, p) {
      figure_spec(
        id = paste0("fig-", tolower(label), "-class-", html_id(key)),
        label = paste0(label, " ", class_label(key)),
        png = class_png(p, key),
        pdf = class_pdf(p, key),
        legend = "By-class panels use the eigen/dominant-attractor curve classification label, then compare finite-time trajectories against the eigenvector result."
      )
    }, names(paths), paths)
  })
  names(class_figs) <- class_keys

  figure_entries <- c(
    lapply(c(fig_gap, fig_consistency, fig_seed, fig_o2), function(spec) list(id = spec$id, label = spec$label)),
    unlist(lapply(class_figs, function(specs) {
      lapply(specs, function(spec) list(id = spec$id, label = spec$label))
    }), recursive = FALSE)
  )

  overview_body <- paste0(
    "<p class=\"report-meta\">This report compares the 1000D and 10000D compare outputs under <strong>",
    html_escape(compare_root),
    "</strong>. The eigenvector method is treated as the fixed-O2 dominant-attractor steady-state ploidy, and finite-time 2N/4N initial-ploidy trajectories are compared against it.</p>",
    table_block("table-validation-summary", "Validation summary", overview_validation),
    table_block("table-run-arguments", "Selected run arguments", overview_args, max_rows = 40L),
    table_block("table-class-consistency", "Terminal-day curve-class consistency summary", consistency_tables, max_rows = 20L)
  )

  by_class_body <- paste0(vapply(class_keys, function(key) {
    section_block(
      paste0("by-class-", html_id(key)),
      class_label(key),
      figure_grid(class_figs[[key]], output_dir = output_dir, embed_assets = embed_assets, columns = 1L, wide = TRUE)
    )
  }, character(1)), collapse = "")

  report_body <- paste0(
    section_block("overview", "1 Overview", overview_body),
    section_block(
      "spectral-gap",
      "2 Spectral-gap comparison",
      figure_grid(fig_gap, output_dir = output_dir, embed_assets = embed_assets, columns = 1L, wide = TRUE)
    ),
    section_block(
      "curve-class-consistency",
      "3 Curve-class consistency",
      figure_grid(fig_consistency, output_dir = output_dir, embed_assets = embed_assets, columns = 2L, wide = TRUE)
    ),
    section_block(
      "seed-level-convergence",
      "4 Seed-level convergence",
      figure_grid(fig_seed, output_dir = output_dir, embed_assets = embed_assets, columns = 2L, wide = TRUE)
    ),
    section_block(
      "o2-level-convergence",
      "5 O2-level convergence",
      paste0(
        figure_grid(list(fig_o2[[1L]]), output_dir = output_dir, embed_assets = embed_assets, columns = 1L, wide = TRUE),
        figure_grid(list(fig_o2[[2L]]), output_dir = output_dir, embed_assets = embed_assets, columns = 1L, wide = TRUE)
      )
    ),
    section_block(
      "by-class-convergence",
      "6 By-class O2 convergence",
      paste0(
        "<p class=\"report-meta\">Each subsection stacks the matching 1000D and 10000D class-specific O2 convergence panels in two rows.</p>",
        by_class_body
      )
    )
  )

  nav <- build_sidebar_nav(section_items, list(), list(), table_entries, figure_entries)
  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>", html_escape(report_title), "</title>",
    "<style>", report_css(), "</style></head><body>",
    "<div class=\"shell\">",
    "<aside class=\"sidebar\"><div class=\"side-head\"><div class=\"kicker\">monotonicity comparison</div>",
    "<div class=\"side-title\">1000D vs 10000D</div><div class=\"side-subtitle\">eigenvector / initial ploidy</div></div><nav>",
    nav,
    "</nav></aside><main class=\"main\">",
    "<section class=\"report-card\" id=\"report-metadata\"><h1>", html_escape(report_title), "</h1>",
    "<p class=\"report-meta\">Side-by-side report for terminal-day compare outputs. PNG previews are embedded; PDF links point to the original figure files.</p>",
    "<p class=\"report-meta\"><strong>Compare root:</strong> ", html_escape(compare_root), "<br>",
    "<strong>Generated at:</strong> ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p></section>",
    report_body,
    "</main></div>",
    report_nav_script(),
    "</body></html>"
  )
}

write_compare_report <- function(compare_root = default_compare_root(),
                                 output_html = NULL,
                                 labels = c("1000D", "10000D"),
                                 embed_assets = TRUE,
                                 report_title = "Initial-Ploidy Trajectory vs Eigenvector Comparison Report") {
  compare_root <- normalizePath(path.expand(compare_root), mustWork = FALSE)
  if (!dir.exists(compare_root)) stop("Compare root does not exist: ", compare_root)
  output_html <- normalizePath(path.expand(output_html %||% file.path(compare_root, "compare_1000D_vs_10000D_report.html")), mustWork = FALSE)
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  html <- build_report_html(
    compare_root = compare_root,
    output_html = output_html,
    labels = labels,
    embed_assets = embed_assets,
    report_title = report_title
  )
  writeLines(html, con = output_html, useBytes = TRUE)
  message("Wrote compare report: ", output_html)
  invisible(output_html)
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
compare_root <- argv$compare_root %||% default_compare_root()
output_html <- argv$output_html %||% file.path(compare_root, "compare_1000D_vs_10000D_report.html")
labels <- if (!is.null(argv$labels)) unlist(strsplit(argv$labels, ",", fixed = TRUE), use.names = FALSE) else c("1000D", "10000D")
labels <- trimws(labels)
write_compare_report(
  compare_root = compare_root,
  output_html = output_html,
  labels = labels,
  embed_assets = as_bool(argv$embed_assets %||% TRUE, default = TRUE),
  report_title = argv$report_title %||% "Initial-Ploidy Trajectory vs Eigenvector Comparison Report"
)
