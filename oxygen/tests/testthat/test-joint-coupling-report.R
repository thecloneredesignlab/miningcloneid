find_joint_report_repo_root <- function() {
  current <- normalizePath(getwd(), mustWork = FALSE)
  for (i in seq_len(8L)) {
    target <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP", "README.md")
    if (file.exists(target)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Cannot locate repository root")
}

write_joint_report_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

count_joint_report_pattern <- function(text, pattern) {
  hit <- gregexpr(pattern, text, perl = TRUE)[[1L]]
  if (identical(hit[[1L]], -1L)) 0L else length(hit)
}

testthat::test_that("joint coupling figure catalog defines all report results", {
  root <- find_joint_report_repo_root()
  report_dir <- file.path(root, "oxygen", "code", "O2_supply_demand_MAP", "report", "joint_coupling")
  catalog_path <- file.path(report_dir, "joint_coupling_figure_catalog.tsv")
  generator_path <- file.path(report_dir, "render_joint_coupling_report.R")
  catalog <- utils::read.delim(catalog_path, check.names = FALSE, stringsAsFactors = FALSE, quote = "")

  required <- c(
    "figure_stem", "major_section_key", "major_section_title", "major_section_order",
    "result_order", "result_id", "subsection_title", "figure_title",
    "legend_explanation", "result_interpretation", "limitation_note"
  )
  testthat::expect_equal(nrow(catalog), 29L)
  testthat::expect_true(all(required %in% names(catalog)))
  testthat::expect_false(anyDuplicated(catalog$figure_stem) > 0L)
  testthat::expect_false(anyDuplicated(catalog$result_id) > 0L)
  testthat::expect_false(any(vapply(catalog[required], function(x) any(is.na(x) | !nzchar(trimws(x))), logical(1L))))
  expected_counts <- c(overview = 2L, within_pair = 5L, between_pair = 4L, process = 2L,
                       ploidy_categories = 5L, category_association = 5L, robustness = 6L)
  testthat::expect_identical(as.integer(table(factor(catalog$major_section_key, levels = names(expected_counts)))), unname(expected_counts))

  generator_text <- paste(readLines(generator_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(generator_text, "validate_figure_catalog")
  testthat::expect_match(generator_text, "IntersectionObserver")
  testthat::expect_match(generator_text, "o2_supply_demand_map_html_utils[.]R")
  testthat::expect_match(generator_text, "o2_supply_demand_map_report_utils[.]R")
  testthat::expect_false(grepl("fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed", generator_text, fixed = TRUE))
})

testthat::test_that("joint coupling report renders 29 navigable, explained figures", {
  testthat::skip_if_not(nzchar(Sys.which("Rscript")), "Rscript is required")
  root <- find_joint_report_repo_root()
  report_dir <- file.path(root, "oxygen", "code", "O2_supply_demand_MAP", "report", "joint_coupling")
  generator <- file.path(report_dir, "render_joint_coupling_report.R")
  catalog <- utils::read.delim(file.path(report_dir, "joint_coupling_figure_catalog.tsv"), check.names = FALSE, stringsAsFactors = FALSE, quote = "")

  fixture <- tempfile("joint-coupling-report-")
  ratio_dir <- file.path(fixture, "tables", "soft_coupling")
  ploidy_dir <- file.path(fixture, "tables", "ploidy_coupling")
  figure_root <- file.path(fixture, "figures")
  out_dir <- file.path(fixture, "report")
  dir.create(ratio_dir, recursive = TRUE)
  dir.create(ploidy_dir, recursive = TRUE)
  dir.create(figure_root, recursive = TRUE)

  write_joint_report_tsv(
    data.frame(parameter = "p", primary_process = "process", cross_pair_dominant_class = "ClassB",
               cross_pair_dominant_fraction = 1, all_pairs_same_dominant_class = TRUE,
               pair_median_direction_consistency = 1, intraclass_correlation = 0.9,
               shared_invitro_anchor = TRUE, invitro_anchor_seed = 10,
               interpretation_scope = "test fixture"),
    file.path(ratio_dir, "soft_coupling_report_summary.tsv")
  )
  write_joint_report_tsv(
    data.frame(key = c("class_threshold", "n_pairs"), value = c("1.1", "6")),
    file.path(ratio_dir, "analysis_config.tsv")
  )
  write_joint_report_tsv(
    data.frame(pair_id = paste0("pair", seq_len(6L)), n_seed = 500, n_parameter = 14,
               n_rows = 7000, expected_rows = 7000, row_count_pass = TRUE),
    file.path(ratio_dir, "input_quality_summary.tsv")
  )
  cats <- rep(c("CatA", "CatB", "CatC"), each = 2L)
  write_joint_report_tsv(
    data.frame(parameter = "p", n = 3000, n_pairs = 6, cramers_v = 0.5,
               association_scope = "descriptive test fixture", primary_process = "process"),
    file.path(ploidy_dir, "ploidy_coupling_report_summary.tsv")
  )
  write_joint_report_tsv(
    data.frame(pair_id = paste0("pair", seq_len(6L)), pair_label = paste("Pair", seq_len(6L)),
               pair_ploidy_category = cats, n_seed = 500, dominant_fraction = 1,
               n_observed_categories = 1, within_pair_category_comparison_estimable = FALSE),
    file.path(ploidy_dir, "ploidy_pair_category_assignment.tsv")
  )
  write_joint_report_tsv(
    data.frame(pair_id = paste0("pair", seq_len(6L)), n_seed = 500, n_observed_categories = 1,
               n_CatA = ifelse(cats == "CatA", 500, 0), n_CatB = ifelse(cats == "CatB", 500, 0),
               n_CatC = ifelse(cats == "CatC", 500, 0),
               within_pair_category_comparison_estimable = FALSE),
    file.path(ploidy_dir, "ploidy_category_estimability.tsv")
  )
  write_joint_report_tsv(
    data.frame(
      setting = c(
        "analysis_window_days", "high_floor", "high_tolerance", "low_endpoint",
        "terminal_window_days", "rolling_slope_threshold_chr_per_day",
        "plateau_min_days", "plateau_abs_slope_limit_chr_per_day",
        "two_transition_bic_delta_cutoff"
      ),
      value = c("0-1000", "44", "0.5", "30", "50", "0.025", "60", "0.02", "-10")
    ),
    file.path(ploidy_dir, "ploidy_category_definition.tsv")
  )

  manifest_rows <- list()
  for (i in seq_len(nrow(catalog))) {
    group_dir <- file.path(figure_root, catalog$major_section_key[[i]])
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    png_path <- file.path(group_dir, paste0(catalog$figure_stem[[i]], ".png"))
    pdf_path <- file.path(group_dir, paste0(catalog$figure_stem[[i]], ".pdf"))
    grDevices::png(png_path, width = 80, height = 60)
    graphics::par(mar = rep(0, 4))
    graphics::plot.new()
    grDevices::dev.off()
    grDevices::pdf(pdf_path, width = 1, height = 1)
    graphics::par(mar = rep(0, 4))
    graphics::plot.new()
    grDevices::dev.off()
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      figure_id = catalog$figure_stem[[i]], format = c("png", "pdf"),
      path = c(png_path, pdf_path), input_tables = "fixture.tsv", chart_family = "Fixture",
      analytical_question = "Does the report render?", width_inches = 1, height_inches = 1
    )
  }
  write_joint_report_tsv(do.call(rbind, manifest_rows), file.path(figure_root, "visualization_manifest.tsv"))

  output <- system2(
    Sys.which("Rscript"),
    c(
      shQuote(generator),
      paste0("--ratio_analysis_dir=", shQuote(ratio_dir)),
      paste0("--ploidy_analysis_dir=", shQuote(ploidy_dir)),
      paste0("--figure_root=", shQuote(figure_root)),
      paste0("--out_dir=", shQuote(out_dir))
    ),
    stdout = TRUE, stderr = TRUE
  )
  testthat::expect_true(is.null(attr(output, "status")), info = paste(output, collapse = "\n"))

  html_path <- file.path(out_dir, "joint_coupling_analysis_report.html")
  html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  testthat::expect_equal(count_joint_report_pattern(html, "<figure\\b"), 29L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='[^']*result-subsection[^']*'"), 29L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='legend-note'"), 29L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='interpretation'"), 29L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='limitation-note'"), 29L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='figure-number'>Figure [0-9]+[.]"), 29L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='report-nav-link report-nav-h3'"), 30L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='report-nav-toggle'"), 7L)
  testthat::expect_equal(count_joint_report_pattern(html, "data-nav-group='[^']+'"), 7L)
  testthat::expect_equal(count_joint_report_pattern(html, "class='report-nav-list report-nav-children'[^>]+hidden"), 7L)
  testthat::expect_match(html, "Figure 1[.]")
  testthat::expect_match(html, "Figure 29[.]")
  testthat::expect_match(html, "class='report-sidebar'")
  testthat::expect_match(html, "class='report-sidebar-header'")
  testthat::expect_match(html, "class='report-main'")
  testthat::expect_match(html, "class='report-card report-header-card'")
  testthat::expect_match(html, "id='report-nav'")
  testthat::expect_match(html, "href='#class-cat-definitions'", fixed = TRUE)
  testthat::expect_match(html, "ratio = fitted in-vivo parameter value / fitted in-vitro parameter value", fixed = TRUE)
  testthat::expect_match(html, "0 &lt; ratio &lt; 0.909091", fixed = TRUE)
  testthat::expect_match(html, "0.909091 ≤ ratio ≤ 1.1", fixed = TRUE)
  testthat::expect_match(html, "data-definition='ClassA'", fixed = TRUE)
  testthat::expect_match(html, "data-definition='ClassB'", fixed = TRUE)
  testthat::expect_match(html, "data-definition='ClassC'", fixed = TRUE)
  testthat::expect_match(html, "data-definition='CatA'", fixed = TRUE)
  testthat::expect_match(html, "data-definition='CatB'", fixed = TRUE)
  testthat::expect_match(html, "data-definition='CatC'", fixed = TRUE)
  testthat::expect_match(html, "data-definition='CatU'", fixed = TRUE)
  testthat::expect_match(html, "stays at or above 43.5 chromosomes", fixed = TRUE)
  testthat::expect_match(html, "middle plateau lasting ≥ 60 days", fixed = TRUE)
  testthat::expect_match(html, "BIC_two − BIC_one ≤ -10", fixed = TRUE)
  testthat::expect_match(html, "two distinct axes", fixed = TRUE)
  testthat::expect_match(html, "IntersectionObserver")
  testthat::expect_match(html, "function setGroupExpanded")
  testthat::expect_match(html, "function revealActiveLink")
  testthat::expect_match(html, "sidebar[.]scrollTo")
  testthat::expect_match(html, "aria-expanded='false'", fixed = TRUE)
  testthat::expect_match(html, "background:#f4f7fa", fixed = TRUE)
  testthat::expect_match(html, "linear-gradient(180deg,#1f3348 0%,#284662 100%)", fixed = TRUE)
  testthat::expect_match(html, "border-radius:12px", fixed = TRUE)
  testthat::expect_match(html, "@media(max-width:1100px)", fixed = TRUE)
  testthat::expect_equal(
    count_joint_report_pattern(html, "id=\"o2sd-report-image-lightbox-script\""),
    1L
  )
  testthat::expect_match(html, "var minimum=Math.max(fitScale*.35,.025);", fixed = TRUE)
  testthat::expect_match(html, "zoom(.8);", fixed = TRUE)
  testthat::expect_false(grepl("prefers-color-scheme:dark", html, fixed = TRUE))
  testthat::expect_false(grepl("[.]png</figcaption>", html))
  testthat::expect_true(all(vapply(catalog$result_id, function(id) {
    grepl(paste0("href='#", id, "'"), html, fixed = TRUE) &&
      grepl(paste0("id='", id, "'"), html, fixed = TRUE)
  }, logical(1L))))
  testthat::expect_true(all(file.exists(file.path(out_dir, c(
    "report_figure_catalog.tsv", "report_manifest.tsv", "chart_map.tsv"
  )))))
  testthat::expect_equal(length(list.files(file.path(out_dir, "figures"), pattern = "[.]png$")), 29L)
  testthat::expect_equal(length(list.files(file.path(out_dir, "figures"), pattern = "[.]pdf$")), 29L)
})
