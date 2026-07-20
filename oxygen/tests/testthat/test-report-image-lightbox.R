find_report_lightbox_repo_root <- function() {
  current <- normalizePath(getwd(), mustWork = FALSE)
  for (i in seq_len(8L)) {
    marker <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP", "README.md")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Cannot locate repository root")
}

count_report_lightbox_pattern <- function(text, pattern, fixed = FALSE) {
  hit <- gregexpr(pattern, text, fixed = fixed, perl = !fixed)[[1L]]
  if (identical(hit[[1L]], -1L)) 0L else length(hit)
}

testthat::test_that("shared lightbox uses the corrected shrink contract", {
  root <- find_report_lightbox_repo_root()
  util <- file.path(
    root, "oxygen", "code", "O2_supply_demand_MAP", "util",
    "o2_supply_demand_map_html_utils.R"
  )
  env <- new.env(parent = baseenv())
  sys.source(util, envir = env)
  assets <- env$o2sd_report_image_lightbox_assets()

  testthat::expect_match(
    assets,
    "fitScale=Math.min(paddedWidth/imageWidth(),paddedHeight/imageHeight(),1);",
    fixed = TRUE
  )
  testthat::expect_match(assets, "var minimum=Math.max(fitScale*.35,.025);", fixed = TRUE)
  testthat::expect_match(assets, "zoom(.8);", fixed = TRUE)
  testthat::expect_match(assets, "1/1.18", fixed = TRUE)
  testthat::expect_match(assets, "var ratio=next/previous;", fixed = TRUE)
  testthat::expect_match(assets, "panX+=dx*(1-ratio);panY+=dy*(1-ratio)", fixed = TRUE)
  testthat::expect_match(assets, "reset(\"fit\")", fixed = TRUE)
  testthat::expect_match(assets, "reset(\"actual\")", fixed = TRUE)
})

testthat::test_that("lightbox injection is content-preserving and idempotent", {
  root <- find_report_lightbox_repo_root()
  util <- file.path(
    root, "oxygen", "code", "O2_supply_demand_MAP", "util",
    "o2_supply_demand_map_html_utils.R"
  )
  env <- new.env(parent = baseenv())
  sys.source(util, envir = env)
  fixture <- tempfile(fileext = ".html")
  original <- paste0(
    "<!doctype html><html><head><title>Fixture</title></head><body>",
    "<nav id='keep-navigation'>Navigation</nav>",
    "<main><figure><img src='data:image/png;base64,AA==' alt='Plot'>",
    "<figcaption>Caption remains unchanged</figcaption></figure></main>",
    "</body></html>"
  )
  writeChar(original, fixture, eos = NULL, useBytes = TRUE)

  first <- env$o2sd_inject_report_image_lightbox(fixture)
  after_first <- paste(readLines(fixture, warn = FALSE), collapse = "\n")
  second <- env$o2sd_inject_report_image_lightbox(fixture)
  after_second <- paste(readLines(fixture, warn = FALSE), collapse = "\n")

  testthat::expect_true(first$changed)
  testthat::expect_identical(first$status, "injected")
  testthat::expect_false(second$changed)
  testthat::expect_identical(second$status, "already_current")
  testthat::expect_identical(after_second, after_first)
  testthat::expect_equal(
    count_report_lightbox_pattern(after_first, "id=\"o2sd-report-image-lightbox-script\"", fixed = TRUE),
    1L
  )
  testthat::expect_match(after_first, "id='keep-navigation'>Navigation", fixed = TRUE)
  testthat::expect_match(after_first, "Caption remains unchanged", fixed = TRUE)
  testthat::expect_match(after_first, "data:image/png;base64,AA==", fixed = TRUE)

  no_image <- tempfile(fileext = ".html")
  writeLines("<!doctype html><html><body><p>No image</p></body></html>", no_image)
  skipped <- env$o2sd_inject_report_image_lightbox(no_image)
  testthat::expect_identical(skipped$status, "skipped_no_images")
  testthat::expect_false(skipped$changed)
})

testthat::test_that("existing-report migration classifies and records every HTML file", {
  root <- find_report_lightbox_repo_root()
  migration_script <- file.path(
    root, "oxygen", "code", "O2_supply_demand_MAP", "report",
    "migrate_existing_html_report_lightboxes.R"
  )
  env <- new.env(parent = globalenv())
  source(migration_script, local = env, chdir = TRUE)
  fixture <- tempfile("lightbox-migration-")
  dir.create(fixture)
  writeLines("<html><body><img src='plot.png'></body></html>", file.path(fixture, "with-image.html"))
  writeLines("<html><body><p>No image</p></body></html>", file.path(fixture, "without-image.html"))

  dry_manifest <- file.path(fixture, "dry-run.tsv")
  dry <- env$o2sd_migrate_existing_html_report_lightboxes(fixture, dry_manifest, dry_run = TRUE)
  testthat::expect_setequal(dry$status, c("would_inject", "skipped_no_images"))
  testthat::expect_true(file.exists(dry_manifest))

  live_manifest <- file.path(fixture, "live.tsv")
  live <- env$o2sd_migrate_existing_html_report_lightboxes(fixture, live_manifest, dry_run = FALSE)
  image_row <- live[live$report_path == "with-image.html", , drop = FALSE]
  testthat::expect_identical(image_row$status, "injected")
  testthat::expect_true(image_row$changed)
  testthat::expect_gt(image_row$bytes_after, image_row$bytes_before)

  rerun <- env$o2sd_migrate_existing_html_report_lightboxes(fixture, file.path(fixture, "rerun.tsv"))
  image_rerun <- rerun[rerun$report_path == "with-image.html", , drop = FALSE]
  testthat::expect_identical(image_rerun$status, "already_current")
  testthat::expect_false(image_rerun$changed)
})

testthat::test_that("all canonical image-report writers install the shared lightbox", {
  root <- find_report_lightbox_repo_root()
  report_root <- file.path(root, "oxygen", "code", "O2_supply_demand_MAP", "report")
  direct_writers <- c(
    "combined_fixo2_eigen/render_fixo2_eigen_attractor_embedding_curve_class_report.R",
    "combined_parameter_landscape/render_combined_parameter_landscape_report.R",
    "fit_results/render_extra_results_report.R",
    "fit_results/render_joint_sigma_soft_coupled_paired_seeds_report.R",
    "fixed_o2_eigen/render_fixo2_eigen_attractor_report.R",
    "joint_coupling/render_joint_coupling_report.R",
    "multi_warmup/render_multi_warmup_results_report.R",
    "parameter_landscape/parameter_landscape_report_utils.R",
    "render_fit_report.R",
    "render_fixo2_invivo_report.R",
    "render_invitro_fit_report.R"
  )
  for (relative in direct_writers) {
    text <- paste(readLines(file.path(report_root, relative), warn = FALSE), collapse = "\n")
    testthat::expect_match(
      text,
      "o2sd_inject_report_image_lightbox(",
      fixed = TRUE,
      info = paste("missing lightbox injection in", relative)
    )
  }

  shared_parameter_writer <- paste(
    readLines(file.path(report_root, "parameter_landscape", "parameter_landscape_report_utils.R"), warn = FALSE),
    collapse = "\n"
  )
  for (relative in c(
    "dense_grid_monotonicity/render_dense_grid_report.R",
    "parameter_landscape/clustering_report.R",
    "parameter_landscape/dominant_ploidy_parameter_contribution_report.R",
    "parameter_landscape/mode_parameter_contribution_report.R"
  )) {
    text <- paste(readLines(file.path(report_root, relative), warn = FALSE), collapse = "\n")
    testthat::expect_match(text, "o2pl_write_report_html(", fixed = TRUE, info = relative)
  }
  testthat::expect_match(shared_parameter_writer, "o2sd_inject_report_image_lightbox(", fixed = TRUE)
})
