#!/usr/bin/env Rscript

# Assemble the materialized joint-coupling tables and visualization outputs into
# one portable technical report. This report layer does not recalculate
# statistics or draw figures.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)) else getwd()
UTIL_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "util"), mustWork = TRUE)
source(file.path(UTIL_DIR, "o2_supply_demand_map_html_utils.R"), local = TRUE)
source(file.path(UTIL_DIR, "o2_supply_demand_map_report_utils.R"), local = TRUE)

parse_args <- o2sd_report_parse_equals_args
`%||%` <- o2sd_report_null_coalesce

escape_html <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  o2sd_html_escape_minimal(x)
}

escape_attr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  o2sd_html_escape_full(x)
}

format_table_value <- function(x) {
  if (length(x) == 0L || is.na(x)) return("NA")
  o2sd_report_format_numeric_like(x)
}

table_html <- function(path, n = Inf) {
  if (!file.exists(path)) stop("Missing report table: ", path)
  x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.finite(n)) x <- utils::head(x, n)
  header <- paste0("<tr>", paste0("<th scope='col'>", escape_html(names(x)), "</th>", collapse = ""), "</tr>")
  if (!nrow(x)) return(paste0("<div class='table-wrap'><table>", header, "</table></div>"))
  rows <- vapply(seq_len(nrow(x)), function(i) {
    values <- vapply(x[i, , drop = FALSE], format_table_value, character(1L))
    paste0("<tr>", paste0("<td>", escape_html(values), "</td>", collapse = ""), "</tr>")
  }, character(1L))
  paste0("<div class='table-wrap'><table>", header, paste(rows, collapse = "\n"), "</table></div>")
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

read_figure_catalog <- function(path) {
  catalog <- read_tsv(path)
  required <- c(
    "figure_stem", "major_section_key", "major_section_title", "major_section_order",
    "result_order", "result_id", "subsection_title", "figure_title",
    "legend_explanation", "result_interpretation", "limitation_note"
  )
  missing_columns <- setdiff(required, names(catalog))
  if (length(missing_columns)) stop("Figure catalog is missing columns: ", paste(missing_columns, collapse = ", "))
  blank <- vapply(catalog[required], function(x) any(is.na(x) | !nzchar(trimws(as.character(x)))), logical(1L))
  if (any(blank)) stop("Figure catalog contains blank values in: ", paste(names(blank)[blank], collapse = ", "))
  catalog$major_section_order <- as.integer(catalog$major_section_order)
  catalog$result_order <- as.integer(catalog$result_order)
  if (anyDuplicated(catalog$figure_stem)) stop("Figure catalog contains duplicate figure_stem values")
  if (anyDuplicated(catalog$result_id)) stop("Figure catalog contains duplicate result_id values")
  section_map <- unique(catalog[c("major_section_key", "major_section_title", "major_section_order")])
  if (anyDuplicated(section_map$major_section_key) || anyDuplicated(section_map$major_section_order)) {
    stop("Each major section key and order must have exactly one title")
  }
  catalog <- catalog[order(catalog$major_section_order, catalog$result_order), , drop = FALSE]
  rownames(catalog) <- NULL
  catalog$figure_order <- seq_len(nrow(catalog))
  catalog
}

relative_path <- function(path, root) {
  path <- normalizePath(path, mustWork = TRUE)
  root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(path, prefix)) stop("Path is not within figure root: ", path)
  substring(path, nchar(prefix) + 1L)
}

validate_figure_catalog <- function(catalog, figure_root) {
  pngs <- list.files(figure_root, pattern = "[.]png$", recursive = TRUE, full.names = TRUE)
  if (!length(pngs)) stop("No PNG figures found under: ", figure_root)
  stems <- tools::file_path_sans_ext(basename(pngs))
  if (anyDuplicated(stems)) stop("Figure stems must be unique across figure_root")
  missing_figures <- setdiff(catalog$figure_stem, stems)
  undocumented_figures <- setdiff(stems, catalog$figure_stem)
  if (length(missing_figures) || length(undocumented_figures)) {
    stop(
      "Figure catalog mismatch. Missing: ", paste(missing_figures, collapse = ", "),
      "; undocumented: ", paste(undocumented_figures, collapse = ", ")
    )
  }
  catalog$source_png <- pngs[match(catalog$figure_stem, stems)]
  catalog$source_pdf <- sub("[.]png$", ".pdf", catalog$source_png)
  missing_pdf <- catalog$source_pdf[!file.exists(catalog$source_pdf)]
  if (length(missing_pdf)) stop("Missing PDF companions: ", paste(missing_pdf, collapse = ", "))
  catalog$source_group <- dirname(vapply(catalog$source_png, relative_path, character(1L), root = figure_root))
  catalog
}

copy_report_assets <- function(catalog, out_dir, figure_root) {
  asset_dir <- file.path(out_dir, "figures")
  dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)
  png_relative <- vapply(catalog$source_png, relative_path, character(1L), root = figure_root)
  pdf_relative <- vapply(catalog$source_pdf, relative_path, character(1L), root = figure_root)
  catalog$png_asset <- paste0("figures/", gsub("[/\\\\]", "__", png_relative))
  catalog$pdf_asset <- paste0("figures/", gsub("[/\\\\]", "__", pdf_relative))
  source_files <- c(catalog$source_png, catalog$source_pdf)
  target_files <- file.path(out_dir, c(catalog$png_asset, catalog$pdf_asset))
  copied <- file.copy(source_files, target_files, overwrite = TRUE, copy.mode = TRUE)
  if (!all(copied)) stop("Failed to copy report assets: ", paste(source_files[!copied], collapse = ", "))
  catalog
}

read_visualization_manifests <- function(figure_root) {
  paths <- list.files(figure_root, pattern = "^visualization_manifest[.]tsv$", recursive = TRUE, full.names = TRUE)
  if (!length(paths)) return(data.frame())
  manifests <- lapply(paths, function(path) {
    x <- read_tsv(path)
    x$manifest_source <- relative_path(path, figure_root)
    x
  })
  columns <- unique(unlist(lapply(manifests, names), use.names = FALSE))
  manifests <- lapply(manifests, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA_character_
    x[columns]
  })
  do.call(rbind, manifests)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

required_args <- c("ratio_analysis_dir", "ploidy_analysis_dir", "figure_root", "out_dir")
missing_args <- required_args[!vapply(required_args, function(x) !is.null(args[[x]]) && nzchar(args[[x]]), logical(1L))]
if (length(missing_args)) stop("Missing required arguments: ", paste(paste0("--", missing_args), collapse = ", "))

ratio_dir <- normalizePath(args$ratio_analysis_dir, mustWork = TRUE)
ploidy_dir <- normalizePath(args$ploidy_analysis_dir, mustWork = TRUE)
figure_root <- normalizePath(args$figure_root, mustWork = TRUE)
out_dir <- normalizePath(args$out_dir, mustWork = FALSE)
catalog_path <- normalizePath(args$figure_catalog %||% file.path(SCRIPT_DIR, "joint_coupling_figure_catalog.tsv"), mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ratio_summary_path <- file.path(ratio_dir, "soft_coupling_report_summary.tsv")
analysis_config_path <- file.path(ratio_dir, "analysis_config.tsv")
input_quality_path <- file.path(ratio_dir, "input_quality_summary.tsv")
ploidy_summary_path <- file.path(ploidy_dir, "ploidy_coupling_report_summary.tsv")
pair_assignment_path <- file.path(ploidy_dir, "ploidy_pair_category_assignment.tsv")
estimability_path <- file.path(ploidy_dir, "ploidy_category_estimability.tsv")
ploidy_definition_path <- file.path(ploidy_dir, "ploidy_category_definition.tsv")

ratio_summary <- read_tsv(ratio_summary_path)
analysis_config <- read_tsv(analysis_config_path)
input_quality <- read_tsv(input_quality_path)
ploidy_summary <- read_tsv(ploidy_summary_path)
pair_assignment <- read_tsv(pair_assignment_path)
estimability <- read_tsv(estimability_path)
ploidy_definition <- read_tsv(ploidy_definition_path)

required_definition_columns <- c("setting", "value")
if (!all(required_definition_columns %in% names(ploidy_definition))) {
  stop("Ploidy category definition table must contain setting and value columns")
}

definition_value <- function(setting, fallback = NA_character_) {
  hit <- ploidy_definition$value[ploidy_definition$setting == setting]
  if (length(hit)) as.character(hit[[1L]]) else fallback
}

definition_number <- function(setting, fallback = NA_real_) {
  value <- suppressWarnings(as.numeric(definition_value(setting, as.character(fallback))))
  if (!is.finite(value)) stop("Invalid numeric ploidy definition for ", setting)
  value
}

config_value <- function(key, fallback = NA_character_) {
  hit <- analysis_config$value[analysis_config$key == key]
  if (length(hit)) as.character(hit[[1L]]) else fallback
}

n_pairs <- as.integer(config_value("n_pairs", nrow(pair_assignment)))
n_parameters <- nrow(ratio_summary)
seed_values <- unique(input_quality$n_seed)
seed_text <- if (length(seed_values) == 1L) as.character(seed_values) else paste(sort(seed_values), collapse = "/")
class_threshold <- as.numeric(config_value("class_threshold", "1.1"))
class_lower <- 1 / class_threshold
analysis_window <- definition_value("analysis_window_days", "0-1000")
high_floor <- definition_number("high_floor", 44)
high_tolerance <- definition_number("high_tolerance", 0.5)
cat_a_floor <- high_floor - high_tolerance
low_endpoint <- definition_number("low_endpoint", 30)
terminal_window_days <- definition_number("terminal_window_days", 50)
rolling_slope_threshold <- definition_number("rolling_slope_threshold_chr_per_day", 0.025)
plateau_min_days <- definition_number("plateau_min_days", 60)
plateau_slope_limit <- definition_number("plateau_abs_slope_limit_chr_per_day", 0.02)
two_transition_bic_cutoff <- definition_number("two_transition_bic_delta_cutoff", -10)
anchor_seeds <- unique(ratio_summary$invitro_anchor_seed)
anchor_text <- paste(anchor_seeds, collapse = ", ")
shared_anchor <- all(as.logical(ratio_summary$shared_invitro_anchor))
cat_levels <- c("CatA", "CatB", "CatC", "CatU")
cat_counts <- table(factor(pair_assignment$pair_ploidy_category, levels = cat_levels))
cat_count_text <- paste0(cat_levels, "=", as.integer(cat_counts), collapse = ", ")
within_pair_estimable <- any(as.logical(estimability$within_pair_category_comparison_estimable))
result_root_text <- config_value("result_root", "joint_coupling_analysis")
fit_label <- basename(result_root_text)
generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

catalog <- read_figure_catalog(catalog_path)
catalog <- validate_figure_catalog(catalog, figure_root)
catalog <- copy_report_assets(catalog, out_dir, figure_root)
if (nrow(catalog) != 29L) stop("The report contract requires exactly 29 figures; found ", nrow(catalog))

catalog_output <- file.path(out_dir, "report_figure_catalog.tsv")
catalog_export <- catalog[c(
  "figure_order", "figure_stem", "major_section_key", "major_section_title",
  "major_section_order", "result_order", "result_id", "subsection_title",
  "figure_title", "legend_explanation", "result_interpretation", "limitation_note",
  "source_group", "png_asset", "pdf_asset"
)]
utils::write.table(catalog_export, catalog_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

section_intros <- c(
  overview = "Start with the pair-by-parameter landscape. The first two results separate class agreement from continuous effect magnitude and establish the reading frame for all later analyses.",
  within_pair = "These results use each pair as the unit of interpretation and the fitting seeds as repeated optimization outcomes. Together they quantify dominant class, threshold crossing, direction within ClassB, and continuous spread.",
  between_pair = "These pair-balanced results ask whether a parameter retains the same ClassA, ClassB, or ClassC behavior across the six pairs. Seed counts are not allowed to turn into artificial cross-pair replication.",
  process = "The audited parameter-to-process map links ratio-class behavior to model mechanisms. Process summaries are descriptive aggregations and should be read with their parameter-level components.",
  ploidy_categories = "These diagnostics define and audit CatA, CatB, and CatC from the in-vivo mean-ploidy trajectories, including both 2N and 4N cohorts and the one-stage versus two-stage transition decision.",
  category_association = "These results compare soft-coupling behavior across the ploidy categories. Because each category contains two pairs and no pair changes category across seeds, all associations are descriptive and pair-confounded.",
  robustness = "The final figure set tests dependence on class thresholds, objective quality, parameter bounds, stable-core membership, and ploidy-category cutoffs. These checks qualify the strength of the preceding patterns."
)

figure_block <- function(row) {
  section_number <- paste0(row$major_section_order, ".", row$result_order)
  n <- row$figure_order
  paste0(
    "<section class='result-subsection report-card report-figure-card' id='", escape_attr(row$result_id), "' aria-labelledby='", escape_attr(row$result_id), "-title' data-figure-number='", n, "'>\n",
    "<h3 id='", escape_attr(row$result_id), "-title'><span class='subsection-number'>", section_number, "</span> ", escape_html(row$subsection_title), "</h3>\n",
    "<figure class='report-figure' id='figure-", n, "' aria-labelledby='figure-caption-", n, "'>",
    "<a class='figure-image-link' href='", escape_attr(row$pdf_asset), "' title='Open vector PDF for Figure ", n, "'>",
    "<img class='report-figure-image' loading='lazy' decoding='async' src='", escape_attr(row$png_asset), "' alt='", escape_attr(row$figure_title), "'></a>",
    "<figcaption class='report-figure-title' id='figure-caption-", n, "'><span class='figure-number'>Figure ", n, ".</span> ", escape_html(row$figure_title),
    " <span class='asset-link'>[<a href='", escape_attr(row$pdf_asset), "'>vector PDF</a>]</span></figcaption></figure>\n",
    "<div class='legend-note'><strong>Legend.</strong> ", escape_html(row$legend_explanation), "</div>\n",
    "<p class='interpretation'><strong>Interpretation.</strong> ", escape_html(row$result_interpretation), "</p>\n",
    "<p class='limitation-note'><strong>Limitation.</strong> ", escape_html(row$limitation_note), "</p>\n",
    "</section>"
  )
}

section_keys <- unique(catalog$major_section_key)
result_sections <- vapply(section_keys, function(key) {
  rows <- catalog[catalog$major_section_key == key, , drop = FALSE]
  paste0(
    "<section class='result-group' id='section-", escape_attr(key), "' aria-labelledby='section-", escape_attr(key), "-title'>\n",
    "<div class='report-section-heading'><h2 id='section-", escape_attr(key), "-title'>", rows$major_section_order[[1L]], ". ", escape_html(rows$major_section_title[[1L]]), "</h2>\n",
    "<p class='section-intro'>", escape_html(section_intros[[key]]), "</p>\n",
    "</div>",
    paste(vapply(seq_len(nrow(rows)), function(i) figure_block(rows[i, , drop = FALSE]), character(1L)), collapse = "\n"),
    "\n</section>"
  )
}, character(1L))

nav_groups <- vapply(section_keys, function(key) {
  rows <- catalog[catalog$major_section_key == key, , drop = FALSE]
  child_list_id <- paste0("nav-children-", key)
  children <- paste0(
    "<li class='report-nav-item'><a class='report-nav-link report-nav-h3' href='#", escape_attr(rows$result_id), "' data-target-id='", escape_attr(rows$result_id), "'><span class='nav-figure'>Fig. ", rows$figure_order,
    "</span><span>", escape_html(rows$subsection_title), "</span></a></li>", collapse = "\n"
  )
  paste0(
    "<li class='report-nav-item nav-group' data-nav-group='", escape_attr(key), "'>",
    "<div class='report-nav-group-row'><a class='report-nav-link report-nav-h2' href='#section-", escape_attr(key), "' data-target-id='section-", escape_attr(key), "'>",
    rows$major_section_order[[1L]], ". ", escape_html(rows$major_section_title[[1L]]), "</a>",
    "<button class='report-nav-toggle' type='button' aria-expanded='false' aria-controls='", escape_attr(child_list_id), "' aria-label='Toggle ", escape_attr(rows$major_section_title[[1L]]), "'><span aria-hidden='true'>›</span></button></div>",
    "<ol class='report-nav-list report-nav-children' id='", escape_attr(child_list_id), "' hidden>", children, "</ol></li>"
  )
}, character(1L))

css <- paste0(
  "*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f4f7fa;color:#1b2a38;line-height:1.55}",
  ".report-shell{display:flex;gap:28px;align-items:flex-start;max-width:1680px;margin:0 auto;padding:24px}.report-sidebar{position:sticky;top:24px;align-self:flex-start;flex:0 0 300px;width:300px;max-height:calc(100vh - 48px);overflow:auto;scrollbar-gutter:stable;border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,.08)}",
  ".report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff}.report-kicker{font-size:11px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;opacity:.78}.report-sidebar-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15}.report-sidebar-subtitle{margin-top:4px;font-size:12px;line-height:1.3;opacity:.85;overflow-wrap:anywhere}",
  "#report-toc>summary{display:none;cursor:pointer;padding:11px 14px;border-bottom:1px solid #d6dde6;color:#17324c;font-size:14px;font-weight:700}.report-nav{padding:10px 8px 12px}.report-nav-list{margin:0;padding:0;list-style:none}.report-nav-item{margin:4px 0}.report-nav-link{display:flex;gap:7px;padding:10px 12px;border-radius:8px;text-decoration:none;color:#17324c;font-size:14px;font-weight:600;line-height:1.35;transition:background-color .15s ease,color .15s ease}.report-nav-link:hover,.report-nav-link:focus-visible{background:rgba(47,110,164,.08);color:#17324c;outline:none}.report-nav-link.is-active{background:#2f6ea4;color:#fff;box-shadow:0 6px 14px rgba(47,110,164,.22)}",
  ".report-nav-group-row{display:grid;grid-template-columns:minmax(0,1fr) 32px;align-items:center}.report-nav-group-row>.report-nav-link{min-width:0}.report-nav-toggle{display:flex;align-items:center;justify-content:center;width:28px;height:28px;padding:0;border:0;border-radius:7px;background:transparent;color:#536577;cursor:pointer}.report-nav-toggle:hover,.report-nav-toggle:focus-visible{background:rgba(47,110,164,.08);color:#17324c;outline:none}.report-nav-toggle span{display:block;font-size:21px;line-height:1;transition:transform .16s ease}.nav-group.is-expanded .report-nav-toggle span{transform:rotate(90deg)}.report-nav-children{margin-left:8px;border-left:1px solid #dfe6ee;padding-left:2px}.report-nav-children[hidden]{display:none}.report-nav-h3{padding:7px 10px 7px 18px;font-size:12px;font-weight:500;color:#536577}.nav-figure{flex:0 0 34px;color:#2f6ea4;font-variant-numeric:tabular-nums}.report-nav-link.is-active .nav-figure{color:#fff}.sidebar-footer{margin:10px 8px 0;padding-top:8px;border-top:1px solid #dfe6ee}.sidebar-footer a{font-size:12px;color:#536577}",
  ".report-main{flex:1;min-width:0;max-width:1180px}.report-card{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.05)}.report-header-card h1{margin:0 0 8px;font-size:28px;line-height:1.15}.report-meta{margin:0;color:#516274;font-size:13px;line-height:1.5}.report-meta strong{color:#1b2a38}.subtitle{margin:12px 0 0;color:#516274;font-size:14px;max-width:920px}",
  "h2,h3,.result-subsection{scroll-margin-top:24px}h2{margin:0 0 12px;color:#17324c;font-size:23px;line-height:1.25}h3{margin:0 0 12px;color:#284662;font-size:18px;line-height:1.3}.report-section-heading{margin:34px 0 16px;padding:0 2px 8px;border-bottom:1px solid #d8dfe7}.report-section-heading h2{margin-bottom:6px}.section-intro{max-width:900px;margin:0;color:#536577;font-size:14px}.subsection-number{color:#2f6ea4;font-variant-numeric:tabular-nums}",
  ".summary{margin:14px 0 0;padding:11px 13px;border-left:4px solid #2f6ea4;background:#f2f6fa}.warning{margin:14px 0 0;padding:11px 13px;border-left:4px solid #d56a27;background:#fff7f1}.summary p{margin:4px 0}.warning p{margin:4px 0}.result-subsection{padding:20px}.report-figure{margin:0 0 10px;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff}.figure-image-link{display:block}.figure-image-link:focus-visible{outline:3px solid #2f6ea4;outline-offset:3px}.report-figure-image{display:block;width:100%;height:auto;border-radius:6px;background:#fff}",
  ".report-figure-title{margin:10px 0 0;color:#1b2a38;font-size:14px;line-height:1.45}.figure-number{font-weight:700;color:#2f6ea4}.asset-link{font-size:12px;white-space:nowrap}.asset-link a{color:#2f6ea4;text-decoration:none}.asset-link a:hover{text-decoration:underline}.legend-note{margin:7px 2px 0;color:#536577;font-size:13px;line-height:1.5}.interpretation{margin:12px 0 0;padding:10px 12px;border-left:4px solid #2f6ea4;background:#f2f6fa;color:#374151;font-size:13px;line-height:1.5}.limitation-note{margin:9px 0 0;padding:9px 12px;border-left:4px solid #d56a27;background:#fff7f1;color:#4b5563;font-size:12px;line-height:1.5}",
  ".definition-section{margin:20px 0 24px;scroll-margin-top:24px}.definition-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}.definition-panel{padding:16px;border:1px solid #d6dde6;border-radius:10px;background:#f9fbfd}.definition-panel h4{margin:0 0 8px;color:#17324c;font-size:16px}.definition-lead{margin:0 0 12px;color:#536577;font-size:13px}.definition-list{display:grid;gap:8px;margin:0}.definition-row{display:grid;grid-template-columns:76px minmax(0,1fr);gap:10px;align-items:start;padding:9px 10px;border:1px solid #e0e6ed;border-radius:8px;background:#fff}.definition-row dt,.definition-row dd{margin:0}.definition-row dd{font-size:13px;color:#374151}.definition-chip{display:inline-block;padding:3px 7px;border-radius:999px;font-size:12px;font-weight:700;text-align:center}.chip-class-a{background:#dceaf6;color:#1f557c}.chip-class-b{background:#fff0bd;color:#765a00}.chip-class-c{background:#fde1d2;color:#8b3f17}.chip-cat-a{background:#d8f0ec;color:#17665d}.chip-cat-b{background:#fde1d8;color:#8a3f2e}.chip-cat-c{background:#eadff4;color:#60417d}.chip-cat-u{background:#e8edf2;color:#4e5c69}.definition-note{margin:12px 0 0;padding-top:10px;border-top:1px solid #dfe6ee;color:#536577;font-size:12px}.definition-crosswalk{margin-top:14px;padding:11px 13px;border-left:4px solid #2f6ea4;background:#f2f6fa;color:#374151;font-size:13px}",
  ".table-wrap{overflow:auto;max-width:100%;margin:12px 0 20px;border:1px solid #d6dde6;border-radius:8px;background:#fff}table{width:100%;min-width:720px;border-collapse:collapse;font-size:12px}th,td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top}th{position:sticky;top:0;z-index:1;background:#f7f9fb;color:#17324c;font-weight:700}tr:last-child td{border-bottom:0}code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;background:#f1f5f9;padding:1px 4px;border-radius:4px}",
  "@media(max-width:1100px){.report-shell{display:block;padding:16px}.report-sidebar{position:relative;top:auto;width:auto;max-height:none;margin-bottom:20px}.report-main{max-width:none}#report-toc>summary{display:block}.report-nav{max-height:65vh;overflow:auto}}",
  "@media(max-width:900px){.report-shell{padding:12px}.report-card{padding:16px;border-radius:10px}.report-header-card h1{font-size:24px}.report-section-heading{margin-top:28px}h2{font-size:21px}h3{font-size:16px}.report-nav-link{font-size:13px}.report-nav-h3{font-size:12px}.definition-grid{grid-template-columns:1fr}.table-wrap{margin-left:0;margin-right:0}}",
  "@media print{body{background:#fff}.report-shell{display:block;max-width:none;padding:0}.report-sidebar{display:none!important}.report-main{max-width:none}.report-card{border:0;box-shadow:none;padding:0;margin-bottom:20px}.result-subsection,.report-figure,.report-figure-title,.legend-note{break-inside:avoid}.report-section-heading{break-after:avoid}.asset-link{display:none}.summary,.warning,.interpretation,.limitation-note{border:1px solid #aaa;background:#fff}}
  "
)

pair_table <- table_html(pair_assignment_path)
estimability_table <- table_html(estimability_path)
ratio_table <- table_html(ratio_summary_path)
ploidy_table <- table_html(ploidy_summary_path)

html <- paste0(
  "<!doctype html>\n<html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>Joint Soft-Coupling and Ploidy-Category Analysis</title><style>", css, "</style></head><body id='report-top'>\n",
  "<div class='report-shell'>\n<aside class='report-sidebar' aria-label='Report navigation'>",
  "<div class='report-sidebar-header'><div class='report-kicker'>Navigation</div><div class='report-sidebar-title'>Joint Coupling Analysis</div><div class='report-sidebar-subtitle'>", escape_html(fit_label), "</div></div>",
  "<details id='report-toc' open><summary>Report contents</summary>",
  "<nav class='report-nav' id='report-nav' aria-label='Report sections'><ol class='report-nav-list'>",
  "<li class='report-nav-item'><a class='report-nav-link report-nav-h2 is-active' href='#report-metadata' data-target-id='report-metadata'>Report metadata</a></li>",
  "<li class='report-nav-item'><a class='report-nav-link report-nav-h2' href='#technical-summary' data-target-id='technical-summary'>1. Technical summary</a></li>",
  "<li class='report-nav-item'><a class='report-nav-link report-nav-h2' href='#scope-definitions' data-target-id='scope-definitions'>2. Scope, data, and definitions</a></li>",
  "<li class='report-nav-item'><a class='report-nav-link report-nav-h3' href='#class-cat-definitions' data-target-id='class-cat-definitions'>2.1 Class and Cat definitions</a></li>",
  paste(nav_groups, collapse = "\n"),
  "<li class='report-nav-item'><a class='report-nav-link report-nav-h2' href='#parameter-lookup' data-target-id='parameter-lookup'>10. Parameter and process lookup</a></li>",
  "<li class='report-nav-item'><a class='report-nav-link report-nav-h2' href='#limits-next-steps' data-target-id='limits-next-steps'>11. Methodological limits and next steps</a></li>",
  "</ol><div class='sidebar-footer'><a href='#report-top'>Back to top</a></div></nav></details></aside>\n",
  "<main class='report-main'>",
  "<header class='report-card report-header-card' id='report-metadata'><h1>Joint Soft-Coupling and Ploidy-Category Analysis</h1>",
  "<p class='report-meta'><strong>Result set:</strong> ", escape_html(fit_label), "<br><strong>Generated at:</strong> ", escape_html(generated_at), "<br><strong>Analysis scale:</strong> ", n_pairs, " pairs × ", escape_html(seed_text), " fitting seeds × ", n_parameters, " soft-coupling parameters</p>",
  "<p class='subtitle'>A figure-complete technical report of within-pair stability, between-pair consistency, biological-process mapping, ploidy trajectory categories, and robustness checks.</p></header>\n",
  "<section class='report-card' id='technical-summary' aria-labelledby='technical-summary-title'><h2 id='technical-summary-title'>1. Technical summary</h2>",
  "<div class='summary'><p>The primary analysis classifies ", n_parameters, " in-vivo / in-vitro parameter ratios for ", escape_html(seed_text), " fitting seeds per pair across ", n_pairs, " pairs. ClassB is the coupling band <code>", format(class_lower, digits = 6), " ≤ ratio ≤ ", format(class_threshold, digits = 6), "</code>; ClassA is below and ClassC is above it.</p>",
  "<p>Within-pair stability, cross-pair stable cores, threshold sensitivity, objective-stratum sensitivity, optimizer-boundary diagnostics, and ploidy-category robustness are reported in Figures 1–", nrow(catalog), ". Pair-balanced summaries use pair—not seed—as the cross-pair unit.</p></div>",
  "<div class='warning'><strong>Estimability constraint.</strong> Pair assignments are ", escape_html(cat_count_text), ". ",
  if (within_pair_estimable) "At least one pair contains within-pair category variation." else "Each pair has only one observed ploidy category across its available seeds. Therefore CatA/B/C comparisons are descriptive between-pair patterns; within-pair category effects and pair-stratified Cat × Class inference are not estimable.",
  "</div></section>\n",
  "<section class='report-card' id='scope-definitions' aria-labelledby='scope-definitions-title'><h2 id='scope-definitions-title'>2. Scope, data, and definitions</h2>",
  "<p>", if (shared_anchor) paste0("All ", n_pairs, " pairs share the same in-vitro seed", escape_html(anchor_text), " anchor") else paste0("The pairs use in-vitro anchors ", escape_html(anchor_text)), ", so cross-pair stability is conditional on that reference. Ratio intervals, quantiles, category calls, and robustness summaries are materialized by the analysis layer before plotting; this report only assembles those outputs.</p>",
  "<section class='definition-section' id='class-cat-definitions' aria-labelledby='class-cat-definitions-title'><h3 id='class-cat-definitions-title'>2.1 Class and Cat definitions</h3>",
  "<div class='definition-grid'>",
  "<article class='definition-panel'><h4>Soft-coupling ratio classes (Class)</h4>",
  "<p class='definition-lead'><code>ratio = fitted in-vivo parameter value / fitted in-vitro parameter value</code>. Classification is performed for every parameter × pair × fitting seed.</p>",
  "<dl class='definition-list'>",
  "<div class='definition-row' data-definition='ClassA'><dt><span class='definition-chip chip-class-a'>ClassA</span></dt><dd><code>0 &lt; ratio &lt; ", format(class_lower, digits = 6), "</code>. The fitted parameter is clearly lower in vivo than in vitro.</dd></div>",
  "<div class='definition-row' data-definition='ClassB'><dt><span class='definition-chip chip-class-b'>ClassB</span></dt><dd><code>", format(class_lower, digits = 6), " ≤ ratio ≤ ", format(class_threshold, digits = 6), "</code>. The values are approximately equal under the prespecified 1.1-fold band; this is the operational soft-coupling class.</dd></div>",
  "<div class='definition-row' data-definition='ClassC'><dt><span class='definition-chip chip-class-c'>ClassC</span></dt><dd><code>ratio &gt; ", format(class_threshold, digits = 6), "</code>. The fitted parameter is clearly higher in vivo than in vitro.</dd></div>",
  "</dl><p class='definition-note'>Both boundary values belong to ClassB. Non-finite or non-positive ratios are invalid rather than assigned to a Class. ClassB denotes numerical agreement within the chosen tolerance; it does not by itself establish causal biological coupling.</p></article>",
  "<article class='definition-panel'><h4>In-vivo ploidy trajectory categories (Cat)</h4>",
  "<p class='definition-lead'>Mean chromosome trajectories are evaluated over days ", escape_html(analysis_window), " for the 2N and 4N cohorts separately, then combined into one seed-level category.</p>",
  "<dl class='definition-list'>",
  "<div class='definition-row' data-definition='CatA'><dt><span class='definition-chip chip-cat-a'>CatA</span></dt><dd>Both cohorts remain in the high-ploidy regime: each trajectory stays at or above ", format(cat_a_floor, digits = 6), " chromosomes (", format(high_floor, digits = 6), " with a ", format(high_tolerance, digits = 6), "-chromosome tolerance).</dd></div>",
  "<div class='definition-row' data-definition='CatB'><dt><span class='definition-chip chip-cat-b'>CatB</span></dt><dd>Both cohorts enter the low-ploidy regime with terminal mean chromosome count ≤ ", format(low_endpoint, digits = 6), ", at least one effective rapid drop, and no qualifying middle plateau: a one-stage decline toward the approximately 22-chromosome state.</dd></div>",
  "<div class='definition-row' data-definition='CatC'><dt><span class='definition-chip chip-cat-c'>CatC</span></dt><dd>A two-stage decline toward the approximately 22-chromosome state: terminal mean ≤ ", format(low_endpoint, digits = 6), ", at least two effective drops, a middle plateau lasting ≥ ", format(plateau_min_days, digits = 6), " days with absolute slope ≤ ", format(plateau_slope_limit, digits = 6), " chromosome/day, and two-transition support <code>BIC_two − BIC_one ≤ ", format(two_transition_bic_cutoff, digits = 6), "</code>. Seed-level CatC also includes the biologically complementary 2N CatB + 4N CatC pattern because 2N begins on the 4N middle plateau.</dd></div>",
  "<div class='definition-row' data-definition='CatU'><dt><span class='definition-chip chip-cat-u'>CatU</span></dt><dd>Unclassified: a missing or edge trajectory, an unsupported shape, or an incompatible high-versus-low cohort combination. CatU is retained for quality control and sensitivity analysis but excluded from prespecified CatA/B/C association tests.</dd></div>",
  "</dl><p class='definition-note'>Terminal behavior is summarized over the final ", format(terminal_window_days, digits = 6), " days. Effective drops use a ", format(rolling_slope_threshold, digits = 6), " chromosome/day rolling-slope threshold. These are operational classification rules, not claims that chromosome count is exactly 44 or 22 at every time point.</p></article>",
  "</div><div class='definition-crosswalk'><strong>How to read Class and Cat together.</strong> Class describes the direction and magnitude of one fitted parameter's in-vivo versus in-vitro difference. Cat describes the shape of the simulated in-vivo ploidy trajectory. The Cat × Class results ask whether these two distinct axes are associated. In this result set, every pair has one Cat across its seeds, so that association is descriptive and confounded with pair identity.</div></section>",
  "<h3>2.2 Pair-level ploidy assignments</h3>", pair_table,
  "<h3>2.3 Within-pair estimability</h3>", estimability_table, "</section>\n",
  paste(result_sections, collapse = "\n"),
  "<section class='report-card' id='parameter-lookup' aria-labelledby='parameter-lookup-title'><h2 id='parameter-lookup-title'>10. Parameter and process lookup</h2>",
  "<p>The tables retain exact values for audit and downstream use; figures are intended for comparison and pattern recognition.</p>",
  "<h3>10.1 Soft-coupling parameter summary</h3>", ratio_table,
  "<h3>10.2 Ploidy category × ratio class summary</h3>", ploidy_table, "</section>\n",
  "<section class='report-card' id='limits-next-steps' aria-labelledby='limits-next-steps-title'><h2 id='limits-next-steps-title'>11. Methodological limits and next steps</h2>",
  "<p>The analysis is model-derived and does not establish a causal biological relationship. Available pair counts are ", escape_html(cat_count_text), ", with no within-pair category variation; pair identity is therefore confounded with ploidy category. The strongest follow-up is to add pairs or fitting seeds that yield multiple ploidy categories within the same pair, then repeat pair-stratified association tests.</p>",
  "<p>For parameter prioritization, require agreement across the primary ClassB threshold, nearby-threshold sensitivity, objective-quality subsets, and boundary/projection diagnostics before interpreting a process as strongly coupled. The per-figure limitation notes identify the specific inferential boundary for each displayed result.</p></section>\n",
  "</main></div>\n",
  "<script>(function(){const toc=document.getElementById('report-toc');const sidebar=document.querySelector('.report-sidebar');const narrow=()=>window.matchMedia('(max-width: 1100px)').matches;if(narrow()){toc.open=false;}const links=Array.from(document.querySelectorAll('#report-nav a[href^=\"#\"]'));const byId=new Map(links.map(a=>[a.getAttribute('href').slice(1),a]));const groups=Array.from(document.querySelectorAll('.nav-group[data-nav-group]'));const targets=Array.from(document.querySelectorAll('main [id]')).filter(e=>byId.has(e.id));let activeId='';let framePending=false;function setGroupExpanded(group,expanded){const list=group.querySelector('.report-nav-children');const toggle=group.querySelector('.report-nav-toggle');group.classList.toggle('is-expanded',expanded);list.hidden=!expanded;toggle.setAttribute('aria-expanded',String(expanded));}function revealActiveLink(link){if(narrow()||!sidebar)return;const linkRect=link.getBoundingClientRect();const sideRect=sidebar.getBoundingClientRect();if(linkRect.top<sideRect.top+12||linkRect.bottom>sideRect.bottom-12){sidebar.scrollTo({top:Math.max(0,sidebar.scrollTop+linkRect.top-sideRect.top-sidebar.clientHeight/2+linkRect.height/2),behavior:'smooth'});}}function activate(id){if(!id||id===activeId)return;activeId=id;links.forEach(a=>a.classList.remove('is-active'));const link=byId.get(id);if(!link)return;link.classList.add('is-active');const activeGroup=link.closest('.nav-group');groups.forEach(group=>setGroupExpanded(group,group===activeGroup));revealActiveLink(link);}function locate(){framePending=false;const marker=Math.max(110,window.innerHeight*.2);let current=targets[0];targets.forEach(target=>{if(target.getBoundingClientRect().top<=marker){current=target;}});if(current)activate(current.id);}function scheduleLocate(){if(framePending)return;framePending=true;window.requestAnimationFrame(locate);}groups.forEach(group=>{const toggle=group.querySelector('.report-nav-toggle');toggle.addEventListener('click',()=>{const expand=!group.classList.contains('is-expanded');groups.forEach(item=>setGroupExpanded(item,false));if(expand)setGroupExpanded(group,true);});});links.forEach(a=>a.addEventListener('click',()=>{const group=a.closest('.nav-group');if(group){groups.forEach(item=>setGroupExpanded(item,item===group));}if(narrow()){toc.open=false;}}));if('IntersectionObserver' in window){const observer=new IntersectionObserver(scheduleLocate,{rootMargin:'-20% 0px -70% 0px',threshold:0});targets.forEach(t=>observer.observe(t));}window.addEventListener('scroll',scheduleLocate,{passive:true});window.addEventListener('resize',scheduleLocate);locate();})();</script>\n",
  "</body></html>"
)

html_path <- file.path(out_dir, "joint_coupling_analysis_report.html")
writeLines(html, html_path, useBytes = TRUE)
o2sd_inject_report_image_lightbox(html_path)

chart_map <- read_visualization_manifests(figure_root)
if (nrow(chart_map)) {
  match_index <- match(chart_map$figure_id, catalog$figure_stem)
  if (anyNA(match_index)) stop("Visualization manifest contains undocumented figure IDs")
  chart_map$report_figure_order <- catalog$figure_order[match_index]
  chart_map$major_section_key <- catalog$major_section_key[match_index]
  chart_map$major_section_title <- catalog$major_section_title[match_index]
  chart_map$result_id <- catalog$result_id[match_index]
  chart_map$figure_title <- catalog$figure_title[match_index]
  chart_map$legend_explanation <- catalog$legend_explanation[match_index]
  chart_map$result_interpretation <- catalog$result_interpretation[match_index]
  chart_map$limitation_note <- catalog$limitation_note[match_index]
  format_order <- match(chart_map$format, c("png", "pdf"))
  chart_map <- chart_map[order(chart_map$report_figure_order, format_order), , drop = FALSE]
}
chart_map_path <- file.path(out_dir, "chart_map.tsv")
utils::write.table(chart_map, chart_map_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

asset_rows <- rbind(
  data.frame(
    artifact_type = "figure_png", figure_number = catalog$figure_order,
    result_id = catalog$result_id, figure_title = catalog$figure_title,
    path = catalog$png_asset, stringsAsFactors = FALSE
  ),
  data.frame(
    artifact_type = "figure_pdf", figure_number = catalog$figure_order,
    result_id = catalog$result_id, figure_title = catalog$figure_title,
    path = catalog$pdf_asset, stringsAsFactors = FALSE
  )
)
report_manifest <- rbind(
  data.frame(artifact_type = "html_report", figure_number = NA_integer_, result_id = NA_character_, figure_title = NA_character_, path = basename(html_path), stringsAsFactors = FALSE),
  data.frame(artifact_type = "figure_catalog", figure_number = NA_integer_, result_id = NA_character_, figure_title = NA_character_, path = basename(catalog_output), stringsAsFactors = FALSE),
  data.frame(artifact_type = "chart_map", figure_number = NA_integer_, result_id = NA_character_, figure_title = NA_character_, path = basename(chart_map_path), stringsAsFactors = FALSE),
  asset_rows
)
utils::write.table(report_manifest, file.path(out_dir, "report_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

message("Wrote figure-complete report: ", html_path)
message("Validated and embedded ", nrow(catalog), " numbered figures with PNG and PDF assets")
