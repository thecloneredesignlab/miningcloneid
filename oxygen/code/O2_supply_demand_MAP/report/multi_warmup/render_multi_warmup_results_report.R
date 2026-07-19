#!/usr/bin/env Rscript

# Report-only HTML assembly. All figures and analytical tables must already
# exist before this module runs.

.o2mw_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_multi_warmup_results_report.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  workflow_frames <- frames[
    grepl("/O2_supply_demand_MAP/", frames, fixed = TRUE)
  ]
  if (length(workflow_frames)) {
    root <- dirname(workflow_frames[[length(workflow_frames)]])
    while (!identical(basename(root), "O2_supply_demand_MAP")) {
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (identical(basename(root), "O2_supply_demand_MAP")) {
      return(file.path(root, "report", "multi_warmup"))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(.o2mw_report_dir, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_html_utils.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_report_utils.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "report", "run_provenance_report.R"), local = environment())

escape_html <- o2sd_html_escape_standard

format_unique_values <- function(x, digits = 6L, scientific = FALSE) {
  if (is.null(x) || !length(x)) return("NA")
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return("NA")
  num <- suppressWarnings(as.numeric(x))
  if (all(is.finite(num))) {
    vals <- unique(signif(num, digits))
    vals <- vals[order(vals)]
    return(paste(format(vals, scientific = scientific, trim = TRUE), collapse = ", "))
  }
  paste(unique(as.character(x)), collapse = ", ")
}

format_count_summary <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  ux <- sort(unique(x))
  if (length(ux) == 1L) return(as.character(ux[[1]]))
  paste0(min(ux), "-", max(ux), " (mixed)")
}

strip_project_path_prefix <- function(x) {
  prefix <- "/share/lab_crd/lab_crd/taoli/Project"
  if (!length(x)) return(x)
  gsub(prefix, "", as.character(x), fixed = TRUE)
}

strip_project_path_prefix_table <- function(tab) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(tab)
  out <- tab
  chr_cols <- vapply(out, is.character, logical(1))
  for (col in names(out)[chr_cols]) {
    out[[col]] <- strip_project_path_prefix(out[[col]])
  }
  out
}

drop_report_columns <- function(tab, cols) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(tab)
  tab[, setdiff(names(tab), cols), drop = FALSE]
}

table_html <- function(tab, max_rows = 50L) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return('<p class="empty">No data available.</p>')
  tab <- utils::head(tab, max_rows)
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", escape_html(names(tab))), collapse = ""), "</tr>")
  rows <- apply(tab, 1L, function(row) paste0("<tr>", paste(sprintf("<td>%s</td>", escape_html(row)), collapse = ""), "</tr>"))
  paste0('<table class="report-table">', header, paste(rows, collapse = ""), "</table>")
}

file_to_data_uri <- function(path, mime = "application/pdf") {
  if (requireNamespace("base64enc", quietly = TRUE)) {
    encoded <- base64enc::base64encode(path)
  } else if (requireNamespace("openssl", quietly = TRUE)) {
    bytes <- readBin(path, what = "raw", n = file.info(path)$size)
    encoded <- openssl::base64_encode(bytes)
  } else {
    return(NULL)
  }
  paste0("data:", mime, ";base64,", encoded)
}

magick_available <- o2sd_report_magick_available
ghostscript_bin <- o2sd_report_ghostscript_bin
ghostscript_available <- o2sd_report_gs_available
render_pdf_preview_png_gs <- o2sd_report_render_pdf_preview_png_gs

pdf_to_image_data_uri <- function(pdf_path, density = 180) {
  if (!magick_available() && !ghostscript_available()) return(NULL)
  png_path <- tempfile("multi_warmup_pdf_preview_", fileext = ".png")
  on.exit(unlink(png_path, force = TRUE), add = TRUE)
  tryCatch({
    if (magick_available()) {
      img <- magick::image_read(pdf_path, density = density)
      if (length(img) > 1L) img <- img[1]
      magick::image_write(img, path = png_path, format = "png")
    } else {
      render_pdf_preview_png_gs(pdf_path, png_path, density = density)
    }
    file_to_data_uri(png_path, mime = "image/png")
  }, error = function(e) NULL)
}

figure_html <- function(root, filename, title, caption, section_id = NULL) {
  path <- file.path(root, filename)
  if (!file.exists(path)) return("")
  section_attr <- if (!is.null(section_id) && nzchar(as.character(section_id))) {
    paste0(' id="', escape_html(section_id), '"')
  } else {
    ""
  }
  figure_data <- pdf_to_image_data_uri(path)
  figure_body <- if (!is.null(figure_data) && nzchar(figure_data)) {
    paste0(
      '<img class="figure-image" src="', escape_html(figure_data), '" alt="', escape_html(title), '"/>'
    )
  } else {
    paste0('<p class="empty">PDF preview rendering is unavailable. Open <code>', escape_html(filename), '</code>.</p>')
  }
  paste0(
    '<section class="card"', section_attr, '>',
    '<h2>', escape_html(title), '</h2>',
    '<p>', escape_html(caption), '</p>',
    figure_body,
    '</section>'
  )
}

insert_after_first_regex <- function(x, pattern, insert) {
  if (!nzchar(insert)) return(x)
  hit <- regexpr(pattern, x, perl = TRUE)
  if (hit[[1]] < 0L) return(x)
  end <- hit[[1]] + attr(hit, "match.length") - 1L
  paste0(substr(x, 1L, end), insert, substr(x, end + 1L, nchar(x)))
}


integrated_tradeoff_by_invivo_warmup_html <- function(root) {
  png_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.png")
  pdf_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.pdf")
  if (!file.exists(png_path) && !file.exists(pdf_path)) return("")
  image_data <- if (file.exists(png_path)) {
    file_to_data_uri(png_path, mime = "image/png")
  } else {
    pdf_to_image_data_uri(pdf_path)
  }
  figure_body <- if (!is.null(image_data) && nzchar(image_data)) {
    paste0(
      '<div class="report-figure"><img src="',
      escape_html(image_data),
      '" alt="Joint In Vivo vs In Vitro Objective Tradeoff by In Vivo Warm-Up Family" class="report-figure-image"/></div>'
    )
  } else {
    '<p class="report-empty">Plot preview rendering is unavailable.</p>'
  }
  paste0(
    '<section class="report-section" id="integrated-joint-figure-4b">',
    figure_body,
    '<h2 class="report-figure-title">Figure 2.4b. Joint In Vivo vs In Vitro Objective Tradeoff by In Vivo Warm-Up Family</h2>',
    '<p class="report-figure-legend">Seed-level in vivo objective against in vitro objective; color shows Vxx in vivo warm-up source. Points are one third of the Figure 2.4 point size and are not text-labeled.</p>',
    '<p class="report-figure-file"><code>multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.pdf</code></p>',
    '</section>'
  )
}
















report_figure_section_html <- function(root, section_id, title, legend, filename_stem, section_style = NULL) {
  png_path <- file.path(root, paste0(filename_stem, ".png"))
  pdf_path <- file.path(root, paste0(filename_stem, ".pdf"))
  if (!file.exists(png_path) && !file.exists(pdf_path)) return("")
  style_attr <- if (!is.null(section_style) && nzchar(as.character(section_style))) {
    paste0(' style="', escape_html(section_style), '"')
  } else {
    ""
  }
  image_data <- if (file.exists(png_path)) {
    file_to_data_uri(png_path, mime = "image/png")
  } else {
    pdf_to_image_data_uri(pdf_path)
  }
  figure_body <- if (!is.null(image_data) && nzchar(image_data)) {
    paste0(
      '<div class="report-figure"><img src="',
      escape_html(image_data),
      '" alt="', escape_html(title), '" class="report-figure-image"/></div>'
    )
  } else {
    '<p class="report-empty">Plot preview rendering is unavailable.</p>'
  }
  paste0(
    '<section class="report-section" id="', escape_html(section_id), '"', style_attr, '>',
    figure_body,
    '<h2 class="report-figure-title">', title, '</h2>',
    '<p class="report-figure-legend">', escape_html(legend), '</p>',
    '<p class="report-figure-file"><code>', escape_html(paste0(filename_stem, ".pdf")), '</code></p>',
    '</section>'
  )
}

format_distance_report_table <- function(tab) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(tab)
  out <- tab
  for (col in names(out)) {
    if (is.numeric(out[[col]])) {
      out[[col]] <- ifelse(
        is.finite(out[[col]]),
        format(signif(out[[col]], 5), trim = TRUE, scientific = FALSE),
        ""
      )
    } else if (is.logical(out[[col]])) {
      out[[col]] <- ifelse(is.na(out[[col]]), "", ifelse(out[[col]], "TRUE", "FALSE"))
    }
  }
  out
}

integrated_distance_diagnostics_html <- function(root) {
  fig_distance <- report_figure_section_html(
    root,
    section_id = "integrated-joint-figure-4c",
    title = "Figure 2.4c. Initial-to-Final 18D Profile Distance by Joint Seed",
    legend = "Each point is one joint seed. The y-axis is the scaled 18D distance from the final fitted profile to that seed's own warm-up profile; color shows the final basin assigned from all completed final seed profiles.",
    filename_stem = "multi_warmup_initial_to_final_distance"
  )
  fig_heatmap <- report_figure_section_html(
    root,
    section_id = "integrated-joint-figure-4d",
    title = "Figure 2.4d. Final Basin Representative vs Warm-Up 18D Profile Distance",
    legend = "Rows are the best-objective representative seed from each final basin; columns are the initial warm-up profiles. Larger distances indicate final profiles farther from the original warm-up regimes.",
    filename_stem = "multi_warmup_final_to_warmup_distance_heatmap"
  )

  basin <- read_table_optional(file.path(root, "multi_warmup_final_basin_distance_summary.tsv"))
  basin_cols <- c(
    "final_basin_id", "n_seeds", "best_warmup_label", "best_joint_seed", "best_objective",
    "distance_to_own_warmup", "nearest_warmup_label", "distance_to_nearest_warmup",
    "distance_to_second_nearest_warmup", "warmup_distance_reference_q25",
    "warmup_distance_reference_median", "is_new_space_candidate"
  )
  if (!is.null(basin) && nrow(basin)) {
    basin <- basin[, intersect(basin_cols, names(basin)), drop = FALSE]
    if ("best_objective" %in% names(basin)) {
      basin <- basin[order(num_col(basin$best_objective), basin$final_basin_id), , drop = FALSE]
    }
  }
  basin_html <- paste0(
    '<section class="report-section" id="integrated-joint-table-4f" style="grid-column:1 / -1;">',
    '<h2 class="report-figure-title">Table 2.4f. Per-Basin Representative Distance Summary</h2>',
    '<p class="report-figure-legend">A basin is flagged as a new-space candidate when its best representative is farther from every warm-up profile than the lower quartile of pairwise warm-up-profile distances. Full table: <code>multi_warmup_final_basin_distance_summary.tsv</code>.</p>',
    table_html(format_distance_report_table(basin), max_rows = 200),
    '</section>'
  )

  paste0(
    fig_distance,
    '<div class="distance-diagnostic-grid">',
    fig_heatmap,
    '<div class="distance-diagnostic-empty" aria-hidden="true"></div>',
    '</div>',
    basin_html
  )
}

replace_report_section_by_id <- function(fragment, section_id, replacement) {
  if (!nzchar(replacement)) return(fragment)
  pattern <- paste0('<section class="report-section" id="', section_id, '">.*?</section>')
  sub(pattern, replacement, fragment, perl = TRUE)
}

replace_integrated_colored_seed_sections <- function(fragment, root) {
  replacements <- list(
    "integrated-joint-figure-5" = list(
      title = "Figure 2.5. Parameter Boundary Forest",
      legend = "Fitted active-parameter positions within bounds. Seed-level points are colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_by_invivo_warmup"
    ),
    "integrated-joint-figure-6" = list(
      title = "Figure 2.5b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter",
      legend = "Soft-coupled paired in vivo and in vitro context-specific parameter positions. In vivo points are colored by Vxx warm-up family; in vitro points are blue; lines connect paired seed-parameter values.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_by_invivo_warmup"
    ),
    "integrated-joint-figure-7" = list(
      title = "Figure 2.6. Parameter Boundary Forest (Pred1000 &gt; 44 Top 3)",
      legend = "Top three integrated seeds passing the 2N/4N 1000-day prediction gate, with seed-level points colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_by_invivo_warmup"
    ),
    "integrated-joint-figure-8" = list(
      title = "Figure 2.6b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Pred1000 &gt; 44 Top 3)",
      legend = "Soft-coupled paired context-specific parameter positions for all integrated seeds, with eligible prediction-gated objective top 3 seeds highlighted. In vivo points are colored by Vxx warm-up family; in vitro points are blue.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_by_invivo_warmup"
    ),
    "integrated-joint-figure-9" = list(
      title = "Figure 2.7. Parameter Boundary Forest (Log X Scale)",
      legend = "Fitted active-parameter positions on the original parameter scale. Seed-level points are colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-10" = list(
      title = "Figure 2.7b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Log X Scale)",
      legend = "Natural-scale soft-coupled paired in vivo and in vitro context-specific parameter positions on a log10 x-axis. In vivo points are colored by Vxx warm-up family; in vitro points are blue.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-11" = list(
      title = "Figure 2.8. Parameter Boundary Forest (Pred1000 &gt; 44 Top 3, Log X Scale)",
      legend = "Top three integrated seeds passing the 2N/4N 1000-day prediction gate on the original parameter scale, with seed-level points colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-12" = list(
      title = "Figure 2.8b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Pred1000 &gt; 44 Top 3, Log X Scale)",
      legend = "Natural-scale soft-coupled paired context-specific parameter positions for all integrated seeds, with eligible prediction-gated objective top 3 seeds highlighted. In vivo points are colored by Vxx warm-up family; in vitro points are blue.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-14" = list(
      title = "Figure 2.10. Cross-seed 1000-day Total Burden Prediction: 2N",
      legend = "Seed-level 2N burden trajectories colored by Vxx in vivo warm-up family; the solid line shows the cross-seed mean.",
      stem = "multi_warmup_integrated_burden_total_log_seed_mean_2N_by_invivo_warmup"
    ),
    "integrated-joint-figure-15" = list(
      title = "Figure 2.11. Cross-seed 1000-day Total Burden Prediction: 4N",
      legend = "Seed-level 4N burden trajectories colored by Vxx in vivo warm-up family; the solid line shows the cross-seed mean.",
      stem = "multi_warmup_integrated_burden_total_log_seed_mean_4N_by_invivo_warmup"
    ),
    "integrated-joint-figure-16" = list(
      title = "Figure 2.12. Cross-seed 1000-day Ploidy Seed Trajectories: 2N",
      legend = "Seed-level 2N ploidy trajectories colored by Vxx in vivo warm-up family; solid and dashed lines show cross-seed mean and median.",
      stem = "multi_warmup_integrated_ploidy_seed_curves_2N_by_invivo_warmup"
    ),
    "integrated-joint-figure-17" = list(
      title = "Figure 2.13. Cross-seed 1000-day Ploidy Seed Trajectories: 4N",
      legend = "Seed-level 4N ploidy trajectories colored by Vxx in vivo warm-up family; solid and dashed lines show cross-seed mean and median.",
      stem = "multi_warmup_integrated_ploidy_seed_curves_4N_by_invivo_warmup"
    )
  )
  for (section_id in names(replacements)) {
    item <- replacements[[section_id]]
    fragment <- replace_report_section_by_id(
      fragment,
      section_id,
      report_figure_section_html(root, section_id, item$title, item$legend, item$stem)
    )
  }
  fragment
}

nav_link <- function(href, label, class = "mw-nav-link") {
  paste0('<li><a class="', class, '" href="#', escape_html(href), '">', escape_html(label), '</a></li>')
}

sidebar_html <- function() {
  part1_items <- c(
    nav_link("part-1", "Part 1. Warm-Up Sweep Summary", "mw-nav-link mw-nav-part"),
    nav_link("run-status", "1.1 Run Status"),
    nav_link("umap-profile", "1.2 18D Warm-Start Profile UMAP"),
    nav_link("umap-cluster", "1.3 18D UMAP by In Vivo Cluster"),
    nav_link("objective-summary", "1.4 Objective Summary"),
    nav_link("objective-tradeoff", "1.5 In Vivo vs In Vitro Objective"),
    nav_link("final-ratios", "1.6 Final Parameter Ratios"),
    nav_link("invivo-only-params", "1.7 In Vivo-Only Warm-Start Parameters"),
    nav_link("deoptim-iterations", "1.8 DEoptim Iteration Frequency"),
    nav_link("warmup-manifest", "1.9 Warm-Up Manifest"),
    nav_link("invivo-cluster-summary", "1.10 In Vivo Cluster Summary"),
    nav_link("best-seed-summary", "1.11 Best Seed Summary"),
    nav_link("final-basin-assignments", "1.12 Final Basin Assignments")
  )
  part2_items <- c(
    nav_link("part-2", "Part 2. Integrated Joint Fitting Results", "mw-nav-link mw-nav-part"),
    nav_link("integrated-joint-chapter-1", "2.1 Joint Summary"),
    nav_link("integrated-joint-convergence-summary", "2.1.1 Convergence Summary"),
    nav_link("integrated-joint-parameter-summary", "2.1.2 Parameter Summary"),
    nav_link("integrated-joint-chapter-2", "2.2 In Vivo"),
    nav_link("integrated-joint-chapter-3", "2.3 In Vitro")
  )
  part3_items <- c(
    nav_link("run-provenance", "Part 3. Run Provenance", "mw-nav-link mw-nav-part")
  )
  paste0(
    '<aside class="mw-sidebar" aria-label="Multi-warm-up report navigation">',
    '<div class="mw-sidebar-header">',
    '<div class="mw-kicker">Navigation</div>',
    '<div class="mw-sidebar-title">Multi-Warm-Up Results</div>',
    '</div>',
    '<nav class="mw-nav"><ul>',
    paste(c(part1_items, part2_items, part3_items), collapse = ""),
    '</ul></nav>',
    '</aside>'
  )
}

regex_extract <- function(x, pattern) {
  m <- regexec(pattern, x, perl = TRUE)
  hit <- regmatches(x, m)[[1]]
  if (length(hit) >= 2L) hit[[2]] else ""
}

rewrite_integrated_joint_fragment <- function(fragment, root) {
  if (!nzchar(fragment)) return(fragment)
  fragment <- sub(
    '(?s)<section class="report-hero">.*?</section>',
    paste0(
      '<section class="report-hero" id="part-2">',
      '<h1>Part 2. Integrated Joint Fitting Results</h1>',
      '<p class="report-meta">',
      'This section reuses the standard joint <code>extra_results.R</code> report logic after combining all warm-up pairs into one prefixed-seed run. ',
      'Seed IDs are written as <code>warmup_label__seed</code> so each row and figure can still be traced back to its warm-up pair.<br/>',
      '<strong>Integrated source:</strong> <code>',
      escape_html(file.path(root, "multi_warmup_integrated_joint_run", "extra_results")),
      '</code>',
      '</p></section>'
    ),
    fragment,
    perl = TRUE
  )
  replacements <- c(
    "1. Joint Summary" = "2.1 Joint Summary",
    "2. In Vivo" = "2.2 In Vivo",
    "3. In Vitro" = "2.3 In Vitro"
  )
  for (old in names(replacements)) {
    fragment <- gsub(old, replacements[[old]], fragment, fixed = TRUE)
  }
  fragment <- gsub('id="chapter-([0-9]+)"', 'id="integrated-joint-chapter-\\1"', fragment, perl = TRUE)
  fragment <- gsub('id="figure-([0-9]+)"', 'id="integrated-joint-figure-\\1"', fragment, perl = TRUE)
  fragment <- gsub('id="convergence-summary"', 'id="integrated-joint-convergence-summary"', fragment, fixed = TRUE)
  fragment <- gsub('id="parameter-summary"', 'id="integrated-joint-parameter-summary"', fragment, fixed = TRUE)
  fragment <- gsub(
    '<h2 class="report-figure-title">Convergence Summary</h2>',
    '<h2 class="report-figure-title">2.1.1 Convergence Summary</h2>',
    fragment,
    fixed = TRUE
  )
  fragment <- gsub(
    '<h2 class="report-figure-title">Parameter Summary</h2>',
    '<h2 class="report-figure-title">2.1.2 Parameter Summary</h2>',
    fragment,
    fixed = TRUE
  )
  fragment <- gsub("Figure ([0-9A-Za-z.]+)\\.", "Figure 2.\\1.", fragment, perl = TRUE)
  fragment <- insert_after_first_regex(
    fragment,
    '<section class="report-section" id="integrated-joint-figure-4">.*?</section>',
    paste0(
      integrated_tradeoff_by_invivo_warmup_html(root),
      integrated_distance_diagnostics_html(root)
    )
  )
  replace_integrated_colored_seed_sections(fragment, root)
}

integrated_joint_results_fragment <- function(root) {
  report_path <- file.path(root, "multi_warmup_integrated_joint_run", "extra_results", "extra_results_report.html")
  if (!file.exists(report_path)) {
    return(list(
      style = "",
      html = paste0(
        '<section class="card">',
        '<h2>Part 2. Integrated Joint Fitting Results</h2>',
        '<p class="empty">Integrated joint extra_results report was not found. Run ',
        '<code>collect_multi_warmup_results.R</code> first to build ',
        '<code>multi_warmup_integrated_joint_run/extra_results</code>.</p>',
        '</section>'
      )
    ))
  }
  txt <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  style <- regex_extract(txt, "(?s)<style>(.*?)</style>")
  main <- regex_extract(txt, '(?s)<main class="report-main">(.*)</main>\\s*</div>\\s*</body>\\s*</html>')
  if (!nzchar(main)) {
    return(list(
      style = style,
      html = paste0(
        '<section class="card">',
        '<h2>Part 2. Integrated Joint Fitting Results</h2>',
        '<p class="empty">Integrated joint report exists but its main content could not be extracted: <code>',
        escape_html(report_path),
        '</code>.</p></section>'
      )
    ))
  }
  list(
    style = style,
    html = paste0('<div class="integrated-joint-report report-main">', rewrite_integrated_joint_fragment(main, root), '</div>')
  )
}

read_pair_seed_summaries <- function(root, manifest) {
  if (is.null(manifest) || !is.data.frame(manifest) || !nrow(manifest) || !("joint_run_prefix" %in% names(manifest))) {
    return(data.frame())
  }
  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    warmup_label <- if ("warmup_label" %in% names(manifest)) as.character(manifest$warmup_label[[i]]) else as.character(i)
    run_dir <- file.path(root, as.character(manifest$joint_run_prefix[[i]]))
    seed_summary <- read_table_optional(file.path(run_dir, "extra_results", "seed_summary.tsv"))
    if (is.null(seed_summary) || !is.data.frame(seed_summary) || !nrow(seed_summary)) next
    seed_summary$warmup_label <- warmup_label
    seed_summary$joint_run_prefix <- as.character(manifest$joint_run_prefix[[i]])
    rows[[length(rows) + 1L]] <- seed_summary
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

run_status_table <- function(root, summary_df, manifest) {
  tasks <- read_table_optional(file.path(root, "multi_warmup_tasks.tsv"))
  seed_summary_all <- read_pair_seed_summaries(root, manifest)
  total <- if (is.data.frame(manifest)) nrow(manifest) else 0L
  completed <- if (is.data.frame(summary_df)) nrow(summary_df) else 0L
  planned_seed_tasks <- if (is.data.frame(tasks) && nrow(tasks)) nrow(tasks) else NA_integer_
  planned_per_pair <- if (is.data.frame(tasks) && nrow(tasks) && all(c("warmup_label", "joint_seed") %in% names(tasks))) {
    stats::aggregate(joint_seed ~ warmup_label, data = tasks, FUN = function(x) length(unique(x)))
  } else {
    data.frame()
  }
  completed_per_pair <- if (is.data.frame(seed_summary_all) && nrow(seed_summary_all) && all(c("warmup_label", "seed") %in% names(seed_summary_all))) {
    stats::aggregate(seed ~ warmup_label, data = seed_summary_all, FUN = function(x) length(unique(x)))
  } else {
    data.frame()
  }
  sigma_values <- c(
    if (is.data.frame(tasks) && "joint_soft_coupling_sigma_default" %in% names(tasks)) tasks$joint_soft_coupling_sigma_default else NULL,
    if (is.data.frame(seed_summary_all) && "joint_soft_coupling_sigma_default" %in% names(seed_summary_all)) seed_summary_all$joint_soft_coupling_sigma_default else NULL
  )
  welsch_c_values <- c(
    if (is.data.frame(tasks) && "joint_soft_coupling_welsch_c" %in% names(tasks)) tasks$joint_soft_coupling_welsch_c else NULL,
    if (is.data.frame(seed_summary_all) && "joint_soft_coupling_welsch_c" %in% names(seed_summary_all)) seed_summary_all$joint_soft_coupling_welsch_c else NULL
  )
  itermax_values <- c(
    if (is.data.frame(tasks) && "itermax" %in% names(tasks)) tasks$itermax else NULL,
    if (is.data.frame(seed_summary_all) && "itermax" %in% names(seed_summary_all)) seed_summary_all$itermax else NULL
  )
  n_soft_params <- if (is.data.frame(seed_summary_all) && "joint_soft_coupling_n_params" %in% names(seed_summary_all)) {
    seed_summary_all$joint_soft_coupling_n_params
  } else {
    NULL
  }
  status <- data.frame(
    metric = c(
      "planned_pairs",
      "completed_pairs",
      "missing_or_failed_pairs",
      "planned_seed_tasks",
      "planned_seeds_per_pair",
      "completed_seed_results",
      "completed_seeds_per_pair",
      "soft_coupling_sigma",
      "soft_coupling_welsch_c",
      "soft_coupling_parameter_count",
      "DEoptim_itermax",
      "qos",
      "walltime"
    ),
    value = c(
      total,
      completed,
      max(total - completed, 0L),
      if (is.na(planned_seed_tasks)) "NA" else planned_seed_tasks,
      format_count_summary(if (nrow(planned_per_pair)) planned_per_pair$joint_seed else integer(0)),
      if (is.data.frame(seed_summary_all) && nrow(seed_summary_all)) nrow(seed_summary_all) else 0L,
      format_count_summary(if (nrow(completed_per_pair)) completed_per_pair$seed else integer(0)),
      format_unique_values(sigma_values),
      format_unique_values(welsch_c_values),
      format_unique_values(n_soft_params),
      format_unique_values(itermax_values),
      if (is.data.frame(tasks) && "qos" %in% names(tasks)) format_unique_values(tasks$qos) else "NA",
      if (is.data.frame(tasks) && "walltime" %in% names(tasks)) format_unique_values(tasks$walltime) else "NA"
    ),
    source = c(
      "multi_warmup_manifest.tsv",
      "multi_warmup_best_seed_summary.tsv",
      "manifest vs completed summary",
      "multi_warmup_tasks.tsv",
      "multi_warmup_tasks.tsv",
      "per-pair extra_results/seed_summary.tsv",
      "per-pair extra_results/seed_summary.tsv",
      "multi_warmup_tasks.tsv / seed_summary.tsv",
      "multi_warmup_tasks.tsv / seed_summary.tsv",
      "per-pair extra_results/seed_summary.tsv",
      "multi_warmup_tasks.tsv / seed_summary.tsv",
      "multi_warmup_tasks.tsv",
      "multi_warmup_tasks.tsv"
    ),
    stringsAsFactors = FALSE
  )
  status
}

status_counts <- function(summary_df, manifest) {
  total <- if (is.data.frame(manifest)) nrow(manifest) else 0L
  completed <- if (is.data.frame(summary_df)) nrow(summary_df) else 0L
  data.frame(
    metric = c("planned_pairs", "completed_pairs", "missing_or_failed_pairs"),
    value = c(total, completed, max(total - completed, 0L)),
    stringsAsFactors = FALSE
  )
}

decision_text <- function(summary_df) {
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) {
    return("No completed warm-up pair summaries were available at report time.")
  }
  if ("final_basin_id" %in% names(summary_df)) {
    basins <- unique(na.omit(summary_df$final_basin_id))
    if (length(basins) <= 1L) {
      return("All completed warm-up pairs assigned to one final parameter basin; this supports warm-start stability for the completed runs.")
    }
    return(paste0("Completed warm-up pairs assigned to ", length(basins), " final parameter basins; compare objective and biological diagnostics before selecting a final result."))
  }
  "Final basin assignments were not available; use objective and diagnostic summaries for the completed runs."
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root <- normalizePath(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir)), mustWork = TRUE)
  out_path <- normalizePath(as_chr(argv$out, file.path(root, "multi-warm-up_results.html")), mustWork = FALSE)
  manifest <- read_table_optional(file.path(root, "multi_warmup_manifest.tsv"))
  invivo_clusters <- read_table_optional(file.path(root, "multi_warmup_invivo_cluster_summary.tsv"))
  summary_df <- read_table_optional(file.path(root, "multi_warmup_best_seed_summary.tsv"))
  basins <- read_table_optional(file.path(root, "multi_warmup_final_basin_assignments.tsv"))

  manifest_display <- strip_project_path_prefix_table(manifest)
  hidden_summary_cols <- c(
    "seed_dir",
    "cluster_feature_set",
    "invivo_seed_dir",
    "invitro_seed_dir",
    "joint_soft_coupling_parameters_table",
    "joint_run_dir"
  )
  invivo_clusters_display <- drop_report_columns(invivo_clusters, hidden_summary_cols)
  summary_display <- drop_report_columns(summary_df, hidden_summary_cols)
  integrated_joint <- integrated_joint_results_fragment(root)
  provenance_block <- o2sd_run_provenance_html(root, title = "Part 3. Run Provenance", section_id = "run-provenance")

  html <- paste0(
    '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    '<title>Multi-Warm-Up Results</title>',
    '<style>',
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    '.mw-layout{display:flex;gap:24px;max-width:1780px;margin:0 auto;padding:24px;}',
    '.shell{flex:1;min-width:0;max-width:none;padding:0;}',
    '.mw-sidebar{position:sticky;top:24px;align-self:flex-start;width:305px;max-height:calc(100vh - 48px);overflow:auto;border:1px solid #d6dde6;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.05);}',
    '.mw-sidebar-header{padding:14px 16px;background:#1f3348;color:#fff;}',
    '.mw-kicker{font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.08em;opacity:.8;}',
    '.mw-sidebar-title{margin-top:4px;font-size:17px;font-weight:700;line-height:1.2;}',
    '.mw-nav{padding:10px;}.mw-nav ul{margin:0;padding:0;list-style:none;}.mw-nav li{margin:3px 0;}',
    '.mw-nav-link{display:block;border-radius:7px;padding:6px 9px;color:#17324c;text-decoration:none;font-size:12px;font-weight:600;line-height:1.3;}',
    '.mw-nav-part{margin-top:6px;background:#eef3f8;font-size:13px;font-weight:800;}',
    '.mw-nav-link:hover{background:rgba(47,110,164,.08);}',
    '.hero,.card{background:#fff;border:1px solid #d6dde6;border-radius:10px;box-shadow:0 8px 22px rgba(0,0,0,.05);padding:16px;margin-bottom:18px;}',
    '.report-section{background:#fff;border:1px solid #d6dde6;border-radius:10px;box-shadow:0 8px 22px rgba(0,0,0,.05);padding:16px;margin-bottom:18px;}',
    'h1{margin:0 0 8px;font-size:28px;}h2{margin:0 0 8px;font-size:18px;}p{color:#425365;line-height:1.45;}',
    '.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px;}',
    '.distance-diagnostic-grid{display:grid;grid-template-columns:minmax(0,2fr) minmax(0,1fr);gap:18px;align-items:start;}',
    '.distance-diagnostic-grid .report-section{min-width:0;margin-bottom:18px;}',
    '.distance-diagnostic-empty{min-height:1px;}',
    '.part-card h2{font-size:20px;}.part-card p{margin-bottom:0;}',
    '.integrated-joint-report{margin-top:18px;}',
    '.integrated-joint-report .report-hero h1{font-size:25px;}',
    '.figure-image{display:block;width:100%;height:auto;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-table{width:100%;border-collapse:collapse;font-size:12px;}.report-table th,.report-table td{padding:7px 9px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}',
    '.report-table th{background:#f7f9fb;font-weight:700;}.empty{font-style:italic;color:#657789;}code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}',
    '@media(max-width:1100px){.mw-layout{display:block;padding:14px}.mw-sidebar{position:relative;top:auto;width:auto;max-height:none;margin-bottom:16px}.grid{grid-template-columns:1fr}}',
    integrated_joint$style,
    '</style></head><body><div class="mw-layout">',
    sidebar_html(),
    '<main class="shell">',
    '<section class="hero"><h1>Multi-Warm-Up Results</h1>',
    '<p><strong>Root:</strong> <code>', escape_html(root), '</code><br/>',
    '<strong>Generated:</strong> ', escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), '<br/>',
    '<strong>Decision summary:</strong> ', escape_html(decision_text(summary_df)), '</p></section>',
    '<section class="card part-card" id="part-1"><h2>Part 1. Warm-Up Sweep Summary</h2>',
    '<p>Cross-pair warm-start selection, run status, objective summaries, paired-parameter ratios, and final basin diagnostics.</p></section>',
    '<section class="card" id="run-status"><h2>1.1 Run Status</h2>', table_html(run_status_table(root, summary_df, manifest), max_rows = 30), '</section>',
    '<div class="grid">',
    figure_html(root, "joint_soft_coupling_18d_profile_umap_500seed.pdf", "1.2 18D Warm-Start Profile UMAP", "Standalone paired 18D warm-start profile UMAP. In vivo rank is shown by color and in vitro rank by point shape.", section_id = "umap-profile"),
    figure_html(root, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_500seed.pdf", "1.3 18D Warm-Start Profile UMAP by In Vivo Cluster", "Same 18D UMAP with in vivo cluster overlay. Red labels mark selected representative in vivo ranks.", section_id = "umap-cluster"),
    '</div>',
    '<div class="grid">',
    figure_html(root, "multi_warmup_objective_summary.pdf", "1.4 Objective Summary", "Best total objective for each completed warm-up pair.", section_id = "objective-summary"),
    figure_html(root, "multi_warmup_invivo_invitro_objective_scatter.pdf", "1.5 In Vivo vs In Vitro Objective", "Best-seed in vivo and in vitro objective components by warm-up pair.", section_id = "objective-tradeoff"),
    figure_html(root, "multi_warmup_final_parameter_ratio_heatmap.pdf", "1.6 Final Parameter Ratios", "Best-seed final in vivo/in vitro parameter ratios.", section_id = "final-ratios"),
    figure_html(root, "multi_warmup_invivo_only_parameter_heatmap.pdf", "1.7 In Vivo-Only Warm-Start Parameters", "Source in vivo seed values for the four parameters not represented in paired soft-coupling ratios.", section_id = "invivo-only-params"),
    '</div>',
    figure_html(root, "multi_warmup_deoptim_iteration_histograms.pdf", "1.8 DEoptim Iteration Frequency by Warm-Up Pair", "Facet histogram of seed-level DEoptim iterations completed for each warm-up pair.", section_id = "deoptim-iterations"),
    '<section class="card" id="warmup-manifest"><h2>1.9 Warm-Up Manifest</h2>', table_html(manifest_display, max_rows = 200), '</section>',
    '<section class="card" id="invivo-cluster-summary"><h2>1.10 In Vivo Cluster Summary</h2>', table_html(invivo_clusters_display, max_rows = 200), '</section>',
    '<section class="card" id="best-seed-summary"><h2>1.11 Best Seed Summary</h2>', table_html(summary_display, max_rows = 200), '</section>',
    '<section class="card" id="final-basin-assignments"><h2>1.12 Final Basin Assignments</h2>', table_html(basins, max_rows = 200), '</section>',
    integrated_joint$html,
    provenance_block,
    '</main></div></body></html>'
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_path, useBytes = TRUE)
  message("Wrote multi-warm-up report: ", out_path)
}

if (identical(environment(), globalenv())) {
  main()
}
