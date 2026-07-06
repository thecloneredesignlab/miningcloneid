#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(
      nzchar,
      vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      }, character(1))
    )
    if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
  }
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
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

repo_root <- function() {
  normalizePath(file.path(SCRIPT_DIR, "..", "..", "..", ".."), mustWork = FALSE)
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

slugify <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("(^-|-$)", "", x)
  ifelse(nzchar(x), x, "section")
}

read_table_optional <- function(path, sep = "\t") {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.table(
      path,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
}

format_numeric_like <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  out <- as.character(x)
  x_trim <- trimws(out)
  numeric_pattern <- "^[-+]?((\\d+\\.?\\d*)|(\\.\\d+))([eE][-+]?\\d+)?$"
  numeric_like <- !is.na(out) & nzchar(x_trim) & grepl(numeric_pattern, x_trim)
  if (!any(numeric_like)) return(out)

  num <- suppressWarnings(as.numeric(x_trim[numeric_like]))
  keep <- is.finite(num)
  if (!any(keep)) return(out)

  out_num <- as.character(num[keep])
  int_like <- abs(num[keep] - round(num[keep])) < 1e-9
  decimal_text <- formatC(num[keep], format = "f", digits = 3)
  sci_nonzero <- !int_like & num[keep] != 0 & grepl("\\.000$", decimal_text)
  decimal <- !int_like & !sci_nonzero

  out_num[int_like] <- format(round(num[keep][int_like]), scientific = FALSE, trim = TRUE, digits = 15)
  out_num[sci_nonzero] <- formatC(num[keep][sci_nonzero], format = "e", digits = 3)
  out_num[decimal] <- formatC(num[keep][decimal], format = "f", digits = 3)

  idx <- which(numeric_like)[keep]
  out[idx] <- out_num
  out
}

table_to_html <- function(df, max_rows = 40, columns = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return("<p class=\"muted\">No rows available.</p>")
  }
  if (!is.null(columns)) {
    keep <- intersect(columns, names(df))
    if (length(keep)) df <- df[, keep, drop = FALSE]
  }
  truncated <- FALSE
  if (nrow(df) > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
    truncated <- TRUE
  }
  df[] <- lapply(df, format_numeric_like)
  header <- paste0("<th>", html_escape(names(df)), "</th>", collapse = "")
  body <- apply(
    df,
    1L,
    function(row) paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  )
  note <- if (truncated) {
    paste0("<p class=\"table-note\">Showing first ", nrow(df), " rows.</p>")
  } else {
    ""
  }
  paste0(
    "<div class=\"table-wrap\"><table><thead><tr>", header, "</tr></thead><tbody>",
    paste0(body, collapse = "\n"),
    "</tbody></table></div>",
    note
  )
}

table_block <- function(table_number, title, df, max_rows = 40, columns = NULL) {
  paste0(
    "<div class=\"table-block\">",
    "<div class=\"table-caption\"><strong>Table ", as.integer(table_number), ". ", html_escape(title), ".</strong></div>",
    table_to_html(df, max_rows = max_rows, columns = columns),
    "</div>"
  )
}

filter_numeric_values <- function(df, col, values, tol = 1e-9) {
  if (is.null(df) || !is.data.frame(df) || !(col %in% names(df))) return(df)
  x <- suppressWarnings(as.numeric(df[[col]]))
  keep <- vapply(x, function(val) any(abs(val - values) <= tol), logical(1))
  df[keep %in% TRUE, , drop = FALSE]
}

file_to_data_uri <- function(path, mime = "image/png") {
  if (!file.exists(path)) stop("File not found for embedding: ", path)
  if (requireNamespace("base64enc", quietly = TRUE)) {
    return(base64enc::dataURI(file = path, mime = mime))
  }
  b64 <- Sys.which("base64")
  if (!nzchar(b64)) stop("The base64enc R package or a base64 executable is required.")
  txt <- paste(system2(b64, path, stdout = TRUE), collapse = "")
  paste0("data:", mime, ";base64,", txt)
}

run_status <- function(expr) {
  out <- tryCatch(expr, error = function(e) {
    structure(character(), status = 1L, error = conditionMessage(e))
  })
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = out)
}

render_pdf_first_page <- function(pdf_path, png_path, dpi = 180) {
  if (!file.exists(pdf_path)) stop("PDF figure not found: ", pdf_path)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(png_path)) {
    png_info <- file.info(png_path)
    pdf_info <- file.info(pdf_path)
    if (!is.na(png_info$size) && png_info$size > 0 &&
        !is.na(png_info$mtime) && !is.na(pdf_info$mtime) &&
        png_info$mtime >= pdf_info$mtime) {
      return(png_path)
    }
  }

  prefix <- sub("\\.png$", "", png_path)
  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    res <- run_status(system2(
      pdftoppm,
      c("-png", "-singlefile", "-f", "1", "-l", "1", "-r", as.character(dpi), pdf_path, prefix),
      stdout = TRUE,
      stderr = TRUE
    ))
    if (identical(res$status, 0L) && file.exists(png_path) && file.info(png_path)$size > 0) {
      return(png_path)
    }
  }

  magick <- Sys.which("magick")
  if (nzchar(magick)) {
    res <- run_status(system2(
      magick,
      c("-density", as.character(dpi), paste0(pdf_path, "[0]"), "-quality", "95", png_path),
      stdout = TRUE,
      stderr = TRUE
    ))
    if (identical(res$status, 0L) && file.exists(png_path) && file.info(png_path)$size > 0) {
      return(png_path)
    }
  }

  stop("Could not render PDF to PNG. Install poppler pdftoppm or ImageMagick magick.")
}

safe_read_args <- function(path) {
  tab <- read_table_optional(path)
  if (!is.data.frame(tab) || !all(c("argument", "value") %in% names(tab))) return(list())
  vals <- as.list(tab$value)
  names(vals) <- tab$argument
  vals
}

mode_counts_table <- function(mode_table) {
  if (is.null(mode_table) || !is.data.frame(mode_table) || !"trajectory_regime" %in% names(mode_table)) {
    return(NULL)
  }
  regs <- mode_table$trajectory_regime
  regs[is.na(regs) | !nzchar(regs)] <- "unlabeled"
  split_rows <- split(mode_table, regs, drop = TRUE)
  rows <- lapply(names(split_rows), function(reg_name) {
    d <- split_rows[[reg_name]]
    mode_vals <- if ("mode_label" %in% names(d)) unique(d$mode_label[!is.na(d$mode_label) & nzchar(d$mode_label)]) else character()
    data.frame(
      trajectory_regime = reg_name,
      mode_label = if (length(mode_vals)) mode_vals[[1]] else NA_character_,
      n_seed_o2 = nrow(d),
      n_seed = if ("seed_id" %in% names(d)) length(unique(d$seed_id)) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$trajectory_regime), , drop = FALSE]
}

mode_summary_sentence <- function(counts, count_col = "n_seed", unit_label = "seeds") {
  if (is.null(counts) || !nrow(counts)) return("")
  if (!count_col %in% names(counts)) count_col <- if ("n_seed_o2" %in% names(counts)) "n_seed_o2" else "n_seed"
  n1 <- counts[[count_col]][match("mode1_attractor_dominant_ploidy_ge_2", counts$trajectory_regime)]
  n2 <- counts[[count_col]][match("mode2_attractor_dominant_ploidy_lt_2", counts$trajectory_regime)]
  parts <- character()
  if (length(n1) && !is.na(n1)) parts <- c(parts, paste0(n1, " mode1 ", unit_label))
  if (length(n2) && !is.na(n2)) parts <- c(parts, paste0(n2, " mode2 ", unit_label))
  if (!length(parts)) return("")
  paste0("The current run contains ", paste(parts, collapse = ", "), ".")
}

arg_value <- function(arg_map, key, fallback = "not recorded") {
  val <- arg_map[[key]]
  if (is.null(val) || !length(val)) return(fallback)
  val <- as.character(val[[1]])
  if (is.na(val) || !nzchar(trimws(val))) fallback else val
}

split_csv_values <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) {
    return(character())
  }
  vals <- trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

compact_csv <- function(x, max_items = 12) {
  vals <- split_csv_values(x)
  if (!length(vals)) return("not recorded")
  if (length(vals) <= max_items) return(paste(vals, collapse = ", "))
  paste0(paste(vals[seq_len(max_items)], collapse = ", "), ", ... (", length(vals), " values)")
}

time_grid_summary <- function(x) {
  vals <- split_csv_values(x)
  if (!length(vals)) return("not recorded")
  num <- suppressWarnings(as.numeric(vals))
  if (length(num) >= 2L && all(is.finite(num))) {
    d <- diff(num)
    regular <- length(d) == 1L || all(abs(d - d[[1]]) < 1e-9)
    if (regular) {
      step <- format_numeric_like(d[[1]])[[1]]
      return(paste0(
        format_numeric_like(min(num))[[1]], " to ",
        format_numeric_like(max(num))[[1]], " days, saved every ",
        step, " day", ifelse(abs(d[[1]] - 1) < 1e-9, "", "s"),
        " (", length(num), " saved time points)"
      ))
    }
  }
  compact_csv(x)
}

unique_value_summary <- function(df, col, fallback = "not recorded") {
  if (is.null(df) || !is.data.frame(df) || !(col %in% names(df))) return(fallback)
  vals <- unique(as.character(df[[col]]))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) return(fallback)
  nums <- suppressWarnings(as.numeric(vals))
  if (all(is.finite(nums))) {
    vals <- format_numeric_like(sort(unique(nums)))
  } else {
    vals <- sort(vals)
  }
  paste(vals, collapse = ", ")
}

simulation_seed_selection_sentence <- function(seed_table, sim_args_map) {
  mode <- arg_value(sim_args_map, "seed_selection_mode", "not recorded")
  ref_o2 <- arg_value(sim_args_map, "mode_reference_o2", "not recorded")
  if (!is.data.frame(seed_table) || !nrow(seed_table)) {
    return(paste0(
      "Representative seed selection was recorded as ",
      html_escape(mode),
      ". The detailed representative-seed table was not available in this result directory."
    ))
  }
  seed_table$final_objective <- if ("final_objective" %in% names(seed_table)) suppressWarnings(as.numeric(seed_table$final_objective)) else NA_real_
  label_one <- function(d) {
    seed <- if ("seed_id" %in% names(d)) d$seed_id[[1]] else "unknown seed"
    mode_label <- if ("mode_label" %in% names(d) && nzchar(as.character(d$mode_label[[1]]))) d$mode_label[[1]] else "unknown mode"
    obj <- if (is.finite(d$final_objective[[1]])) format(d$final_objective[[1]], scientific = TRUE, digits = 4) else "NA"
    paste0(seed, " for ", mode_label, " (final objective ", obj, ")")
  }
  pieces <- vapply(split(seed_table, seq_len(nrow(seed_table))), label_one, character(1))
  if (identical(mode, "best_final_objective_by_reference_mode")) {
    paste0(
      "Representative seeds were selected automatically by first assigning each seed to mode1 or mode2 using the fixed-O2 attractor at O2 = ",
      html_escape(ref_o2),
      ", then choosing the seed with the smallest final objective within each mode. The selected seeds were ",
      html_escape(paste(pieces, collapse = "; ")),
      "."
    )
  } else {
    paste0(
      "Representative seeds were supplied manually for this run. The workflow attached the FixO2 reference mode at O2 = ",
      html_escape(ref_o2),
      " and the final objective when available. The analyzed seeds were ",
      html_escape(paste(pieces, collapse = "; ")),
      "."
    )
  }
}

code_span <- function(x) {
  paste0("<code>", html_escape(x), "</code>")
}

figure_png_name <- function(section, pdf_path, relative_path = NULL) {
  source_path <- relative_path %||% pdf_path
  source_path <- gsub("\\\\", "/", source_path)
  paste0(slugify(section), "__", slugify(tools::file_path_sans_ext(source_path)), ".png")
}

figure_id <- function(spec) {
  pdf_name <- basename(spec$relative_path %||% spec$title)
  slugify(paste("figure", spec$figure_number %||% "", tools::file_path_sans_ext(pdf_name), sep = "-"))
}

number_figures <- function(specs) {
  lapply(seq_along(specs), function(i) {
    specs[[i]]$figure_number <- i
    specs[[i]]
  })
}

figure_block <- function(spec, analysis_dir, asset_dir) {
  if (identical(spec$type %||% "single", "paired")) {
    return(paired_figure_block(spec, analysis_dir, asset_dir))
  }
  pdf_path <- file.path(analysis_dir, spec$relative_path)
  id <- figure_id(spec)
  figure_label <- paste0("Figure ", as.integer(spec$figure_number), ". ", spec$title, ".")
  if (!file.exists(pdf_path)) {
    return(paste0(
      "<article class=\"figure\" id=\"", html_escape(id), "\">",
      "<h3>", html_escape(figure_label), "</h3>",
      "<p class=\"missing\">Missing figure: ", html_escape(pdf_path), "</p>",
      "</article>"
    ))
  }
  png_path <- file.path(asset_dir, figure_png_name(spec$section, pdf_path, spec$relative_path))
  png_path <- render_pdf_first_page(pdf_path, png_path)
  uri <- file_to_data_uri(png_path, "image/png")
  paste0(
    "<article class=\"figure\" id=\"", html_escape(id), "\">",
    "<img src=\"", uri, "\" alt=\"", html_escape(spec$title), "\">",
    "<figcaption class=\"figure-caption\"><strong>", html_escape(figure_label), "</strong> ",
    html_escape(spec$what_changes), " ", html_escape(spec$legend), "</figcaption>",
    "<p class=\"figure-result\">", html_escape(spec$result), "</p>",
    "</article>"
  )
}

figure_panel <- function(panel, spec, analysis_dir, asset_dir) {
  pdf_path <- file.path(analysis_dir, panel$relative_path)
  panel_label <- paste0("Figure ", as.integer(spec$figure_number), panel$suffix, ". ", panel$title, ".")
  if (!file.exists(pdf_path)) {
    return(paste0(
      "<div class=\"figure-panel missing-panel\">",
      "<h4>", html_escape(panel_label), "</h4>",
      "<p class=\"missing\">Missing figure: ", html_escape(pdf_path), "</p>",
      "</div>"
    ))
  }
  png_path <- file.path(asset_dir, figure_png_name(spec$section, pdf_path, panel$relative_path))
  png_path <- render_pdf_first_page(pdf_path, png_path)
  uri <- file_to_data_uri(png_path, "image/png")
  paste0(
    "<div class=\"figure-panel\">",
    "<img src=\"", uri, "\" alt=\"", html_escape(panel$title), "\">",
    "<figcaption class=\"figure-caption\"><strong>", html_escape(panel_label), "</strong> ",
    html_escape(panel$what_changes), " ", html_escape(panel$legend), "</figcaption>",
    "<p class=\"figure-result\">", html_escape(panel$result), "</p>",
    "</div>"
  )
}

paired_figure_block <- function(spec, analysis_dir, asset_dir) {
  id <- figure_id(spec)
  suffixes <- vapply(spec$panels, function(panel) panel$suffix %||% "", character(1))
  suffixes <- suffixes[nzchar(suffixes)]
  suffix_range <- if (length(suffixes) > 1L) paste0(suffixes[[1]], "-", suffixes[[length(suffixes)]]) else suffixes[[1]]
  paste0(
    "<article class=\"figure paired-figure\" id=\"", html_escape(id), "\">",
    "<h3>Figure ", as.integer(spec$figure_number), html_escape(suffix_range), ". ", html_escape(spec$title), ".</h3>",
    paste0(vapply(spec$panels, figure_panel, character(1), spec = spec, analysis_dir = analysis_dir, asset_dir = asset_dir), collapse = ""),
    "</article>"
  )
}

build_nav <- function(items) {
  paste0(
    "<ul>",
    paste0(vapply(items, function(item) {
      cls <- paste0("nav-", item$level)
      paste0("<li><a class=\"", cls, "\" href=\"#", html_escape(item$id), "\">", html_escape(item$title), "</a></li>")
    }, character(1)), collapse = ""),
    "</ul>"
  )
}

chapter_header <- function(id, title, body) {
  paste0(
    "<section id=\"", html_escape(id), "\" class=\"chapter\">",
    "<h2>", html_escape(title), "</h2>",
    body,
    "</section>"
  )
}

simulation_figure_spec <- function(filename) {
  key <- tools::file_path_sans_ext(filename)
  if (grepl("^expm_vs_euler_mean_ploidy_", key)) {
    spec <- list(
      title = "Expm Analytical Solution vs Euler Numerical Integration",
      what_changes = "Rows split the 2N and 4N initial conditions, columns split the representative seeds, and colored curves show fixed O2 levels.",
      legend = "Solid colored curves are matrix-exponential analytical trajectories. Dashed curves are Euler numerical integration with the plotted dt.",
      result = "This comparison evaluates whether the finite-step numerical integration reproduces the matrix-exponential solution of the same fixed-O2 linear system before the analytical trajectories are compared with simulation output."
    )
  } else if (grepl("^eigen_vs_simulation_replicates_mean_ploidy_", key)) {
    spec <- list(
      title = "Eigen Analytical Solution vs Simulation Replicates",
      what_changes = "Each panel fixes seed and initial ploidy. O2 changes by color, and replicate simulations are overlaid against the eigen-based analytical trajectory.",
      legend = "Solid translucent curves are eigen analytical trajectories. Dashed curves are individual simulation replicates.",
      result = "This comparison tests whether the modal representation of the fixed-O2 generator captures the replicate-level simulation trajectories across the representative seeds."
    )
  } else if (grepl("^expm_vs_simulation_replicates_mean_ploidy_", key)) {
    spec <- list(
      title = "Expm Analytical Solution vs Simulation Replicates",
      what_changes = "Each panel fixes seed and initial ploidy. O2 changes by color, and replicate simulations are overlaid against the matrix-exponential trajectory.",
      legend = "Solid translucent curves are expm analytical trajectories. Dashed curves are individual simulation replicates.",
      result = "This comparison uses the direct finite-time solution as the deterministic reference for the fixed-O2 simulation replicates."
    )
  } else if (grepl("^expm_vs_simulation_phase_plane_mean_sd_ploidy_", key)) {
    spec <- list(
      title = "Phase Plane: Expm Analytical Solution vs Simulation",
      what_changes = "The x-axis is mean ploidy and the y-axis is ploidy standard deviation. Curves trace the trajectory through this phase plane under each fixed O2 level.",
      legend = "Solid translucent paths are expm analytical trajectories. Dashed paths are simulation replicates.",
      result = "This phase-plane view evaluates whether the simulation replicates preserve the joint evolution of mean ploidy and ploidy dispersion predicted by the matrix-exponential solution."
    )
  } else if (grepl("^eigen_vs_simulation_phase_plane_mean_sd_ploidy_", key)) {
    spec <- list(
      title = "Phase Plane: Eigen Analytical Solution vs Simulation",
      what_changes = "The x-axis is mean ploidy and the y-axis is ploidy standard deviation. Curves trace the trajectory through this phase plane under each fixed O2 level.",
      legend = "Solid translucent paths are eigen analytical trajectories. Dashed paths are simulation replicates.",
      result = "This phase-plane view evaluates whether the modal solution reproduces the same mean-SD trajectory geometry observed in the fixed-O2 simulation replicates."
    )
  } else {
    spec <- list(
      title = gsub("_", " ", tools::file_path_sans_ext(filename), fixed = TRUE),
      what_changes = "The figure is one of the generated simulation validation panels.",
      legend = "See the labels in the source figure.",
      result = "This figure is included because it was present in simulation/figures."
    )
  }
  spec$section <- "Simulation"
  spec$relative_path <- file.path("simulation", "figures", filename)
  spec
}

simulation_figure_order <- function(filename) {
  key <- tools::file_path_sans_ext(filename)
  if (grepl("^expm_vs_euler_mean_ploidy_", key)) return(1L)
  if (grepl("^eigen_vs_simulation_replicates_mean_ploidy_", key)) return(2L)
  if (grepl("^expm_vs_simulation_replicates_mean_ploidy_", key)) return(3L)
  if (grepl("^expm_vs_simulation_phase_plane_mean_sd_ploidy_", key)) return(4L)
  if (grepl("^eigen_vs_simulation_phase_plane_mean_sd_ploidy_", key)) return(5L)
  99L
}

agreement_method_label <- function(method) {
  method <- tolower(as.character(method[[1]]))
  if (identical(method, "eigen")) return("Eigen analytical")
  if (identical(method, "expm")) return("Expm analytical")
  paste0(method, " analytical")
}

agreement_method_summary <- function(analysis_dir, method) {
  path <- file.path(
    analysis_dir,
    "simulation",
    "analytical_simulation_agreement",
    method,
    "tables",
    "agreement_analytical_vs_simulation_data.tsv"
  )
  dat <- read_table_optional(path)
  if (!is.data.frame(dat) || !nrow(dat)) {
    path <- file.path(
      analysis_dir,
      "simulation",
      "analytical_simulation_agreement",
      "tables",
      method,
      "agreement_analytical_vs_simulation_data.tsv"
    )
    dat <- read_table_optional(path)
  }
  if (!is.data.frame(dat) || !nrow(dat)) {
    return(list(
      n = NA_integer_,
      n_seed = NA_integer_,
      n_o2 = NA_integer_,
      n_day = NA_integer_,
      bias = NA_real_,
      median_abs = NA_real_,
      p95_abs = NA_real_,
      cor = NA_real_
    ))
  }
  analytical <- suppressWarnings(as.numeric(dat$analytical_mean_ploidy))
  simulation <- suppressWarnings(as.numeric(dat$simulation_mean_ploidy))
  residual <- simulation - analytical
  ok <- is.finite(analytical) & is.finite(simulation) & is.finite(residual)
  list(
    n = sum(ok),
    n_seed = if ("seed_id" %in% names(dat)) length(unique(dat$seed_id[ok])) else NA_integer_,
    n_o2 = if ("O2_pct" %in% names(dat)) length(unique(dat$O2_pct[ok])) else NA_integer_,
    n_day = if ("day" %in% names(dat)) length(unique(dat$day[ok])) else NA_integer_,
    bias = if (any(ok)) mean(residual[ok], na.rm = TRUE) else NA_real_,
    median_abs = if (any(ok)) stats::median(abs(residual[ok]), na.rm = TRUE) else NA_real_,
    p95_abs = if (any(ok)) as.numeric(stats::quantile(abs(residual[ok]), 0.95, na.rm = TRUE, names = FALSE)) else NA_real_,
    cor = if (sum(ok) >= 3L) suppressWarnings(stats::cor(analytical[ok], simulation[ok], method = "pearson")) else NA_real_
  )
}

agreement_summary_sentence <- function(summary) {
  if (is.null(summary) || !is.finite(summary$n)) {
    return("The agreement table was unavailable, so the result is interpreted qualitatively from the plotted point cloud and residual distribution.")
  }
  paste0(
    "Across ", format(summary$n, big.mark = ","), " matched seed-O2-initial-condition-time rows from ",
    format(summary$n_seed, big.mark = ","), " seeds, the mean simulation-minus-analytical bias was ",
    format_numeric_like(summary$bias)[[1]], " ploidy units, the median absolute residual was ",
    format_numeric_like(summary$median_abs)[[1]], ", the 95th percentile absolute residual was ",
    format_numeric_like(summary$p95_abs)[[1]], ", and the analytical-simulation Pearson correlation was ",
    format_numeric_like(summary$cor)[[1]], "."
  )
}

agreement_pair_spec <- function(method, key, title_suffix, scatter_filename, bland_filename, summary, focus) {
  label <- agreement_method_label(method)
  summary_sentence <- agreement_summary_sentence(summary)
  list(
    type = "paired",
    method = method,
    section = "Analytical-Simulation Agreement",
    title = paste0(label, " Agreement: ", title_suffix),
    panels = list(
      list(
        suffix = "A",
        relative_path = file.path("simulation", "analytical_simulation_agreement", method, "scatter", scatter_filename),
        title = paste0(label, " scatter: ", title_suffix),
        what_changes = "The x-axis is analytical mean ploidy and the y-axis is simulation-inferred mean ploidy. The diagonal reference line marks exact agreement.",
        legend = focus$scatter_legend,
        result = paste0(focus$scatter_result, " ", summary_sentence)
      ),
      list(
        suffix = "B",
        relative_path = file.path("simulation", "analytical_simulation_agreement", method, "bland_altman", bland_filename),
        title = paste0(label, " Bland-Altman residuals: ", title_suffix),
        what_changes = "The x-axis is the mean of analytical and simulation mean ploidy, and the y-axis is simulation minus analytical mean ploidy.",
        legend = "The solid horizontal line is the mean residual bias, and dashed horizontal lines show the approximate 95% limits of agreement defined as bias +/- 1.96 standard deviations of the residuals.",
        result = paste0(focus$bland_result, " ", summary_sentence)
      )
    )
  )
}

agreement_plot_spec <- function(method, summary) {
  label <- agreement_method_label(method)
  summary_sentence <- agreement_summary_sentence(summary)
  list(
    type = "paired",
    method = method,
    section = "Analytical-Simulation Agreement",
    title = paste0(label, " Solution Distributions Across Fixed O2"),
    panels = list(
      list(
        suffix = "A",
        relative_path = file.path("simulation", "analytical_simulation_agreement", method, "plots", "analytical_solution_vs_fixed_o2_by_time_init_half_violin_median.pdf"),
        title = paste0(label, " solution distribution across fixed O2"),
        what_changes = "Fixed O2 is shown as an equally spaced categorical x-axis, with one panel for each selected time point.",
        legend = "Half violins show the seed-level distribution of analytical mean ploidy for init_2N in green and init_4N in purple; the white point and short line mark the median.",
        result = paste0(
          "This distributional view summarizes how the deterministic analytical solution varies with fixed O2 before simulation residuals are considered. ",
          summary_sentence
        )
      ),
      list(
        suffix = "B",
        relative_path = file.path("simulation", "analytical_simulation_agreement", method, "plots", "analytical_solution_vs_fixed_o2_by_time_mode_init_half_violin_median.pdf"),
        title = paste0(label, " solution distribution by mode and initial condition"),
        what_changes = "For each selected time point, fixed O2 is expanded into alternating mode1 and mode2 groups. Mode1 groups have a pale blue background and mode2 groups have a pale orange background.",
        legend = "Within each mode group, half violins show the seed-level distribution of analytical mean ploidy for init_2N in green and init_4N in purple; the white point and short line mark the median.",
        result = paste0(
          "This mode-resolved distributional view compares the spread of analytical solutions across fixed O2, fitted trajectory mode, and initial condition. ",
          summary_sentence
        )
      )
    )
  )
}

agreement_figure_specs <- function(analysis_dir) {
  focuses <- list(
    list(
      title_suffix = "All Selected Time Points Colored by Objective",
      scatter = "scatter_analytical_vs_simulation_by_o2_color_objective_all_times.pdf",
      bland = "bland_altman_simulation_vs_analytical_by_o2_color_objective_all_times.pdf",
      scatter_legend = "Panels split fixed O2, point color encodes the final fitting objective, and point shape marks the initial condition.",
      bland_legend = "",
      scatter_result = "This view asks whether agreement changes across O2 levels or is concentrated among seeds with poorer final objective values.",
      bland_result = "The residual view tests whether the residual magnitude or direction shifts with the average ploidy scale and whether objective-colored points form outlying residual clusters."
    ),
    list(
      title_suffix = "Selected Time Points Colored by Fixed O2",
      scatter = "scatter_analytical_vs_simulation_by_time_color_o2.pdf",
      bland = "bland_altman_simulation_vs_analytical_by_time_color_o2.pdf",
      scatter_legend = "Panels split the eight selected time points, point color marks fixed O2, and point shape marks the initial condition.",
      bland_legend = "",
      scatter_result = "This view asks whether analytical-simulation agreement changes over time and whether particular O2 levels separate from the main agreement diagonal.",
      bland_result = "The residual view shows whether bias broadens or shifts at specific time points while keeping fixed O2 visible by color."
    ),
    list(
      title_suffix = "Selected Time Points Colored by Objective",
      scatter = "scatter_analytical_vs_simulation_by_time_color_objective.pdf",
      bland = "bland_altman_simulation_vs_analytical_by_time_color_objective.pdf",
      scatter_legend = "Panels split the eight selected time points, point color encodes the final fitting objective, and point shape marks the initial condition.",
      bland_legend = "",
      scatter_result = "This view asks whether seeds with larger final objective values are responsible for deviations from the analytical-simulation agreement diagonal.",
      bland_result = "The residual view tests whether fitting objective explains the spread around zero residual across time."
    ),
    list(
      title_suffix = "Selected Time Points Colored by Mode",
      scatter = "scatter_analytical_vs_simulation_by_time_color_mode.pdf",
      bland = "bland_altman_simulation_vs_analytical_by_time_color_mode.pdf",
      scatter_legend = "Panels split the eight selected time points, point color marks the FixO2 attractor-derived mode or unknown mode status, and point shape marks the initial condition.",
      bland_legend = "",
      scatter_result = "This view asks whether mode1 and mode2 occupy different analytical-simulation agreement regions.",
      bland_result = "The residual view tests whether residual bias differs by the two FixO2 attractor-derived modes."
    )
  )
  methods <- c("eigen", "expm")
  specs <- list()
  k <- 0L
  for (method in methods) {
    summary <- agreement_method_summary(analysis_dir, method)
    for (focus in focuses) {
      k <- k + 1L
      specs[[k]] <- agreement_pair_spec(
        method = method,
        key = focus$title_suffix,
        title_suffix = focus$title_suffix,
        scatter_filename = focus$scatter,
        bland_filename = focus$bland,
        summary = summary,
        focus = focus
      )
    }
    k <- k + 1L
    specs[[k]] <- agreement_plot_spec(method, summary)
  }
  specs
}

load_report_data <- function(analysis_dir) {
  attractor_dir <- file.path(analysis_dir, "attractors")
  counter_dir <- file.path(analysis_dir, "counterfactual_trajectories")
  simulation_dir <- file.path(analysis_dir, "simulation")
  list(
    top_args = read_table_optional(file.path(analysis_dir, "tables", "FixO2_invivo_run_arguments.tsv")),
    top_args_map = safe_read_args(file.path(analysis_dir, "tables", "FixO2_invivo_run_arguments.tsv")),
    manifest = read_table_optional(file.path(attractor_dir, "tables", "fixo2_seed_metadata.tsv")),
    attractor_mode_reference = read_table_optional(file.path(attractor_dir, "tables", "fixed_o2_attractor_mode_reference_by_seed.tsv")),
    attractor_mode = read_table_optional(file.path(attractor_dir, "tables", "fixed_o2_attractor_mode_by_seed_o2.tsv")),
    attractor_display = read_table_optional(file.path(attractor_dir, "tables", "dominant_mean_ploidy_summary_mode1_mode2_display.tsv")),
    attractor_gap = read_table_optional(file.path(attractor_dir, "tables", "fixed_o2_attractor_spectral_gap_regime_summary.tsv")),
    attractor_tests = read_table_optional(file.path(attractor_dir, "tables", "fixed_o2_attractor_regime_tests.tsv")),
    counter_args = read_table_optional(file.path(counter_dir, "tables", "analysis_run_arguments.tsv")),
    counter_summary = read_table_optional(file.path(counter_dir, "tables", "fixed_o2_counterfactual_regime_summary.tsv")),
    counter_tests = read_table_optional(file.path(counter_dir, "tables", "fixed_o2_counterfactual_regime_tests.tsv")),
    simulation_args = read_table_optional(file.path(simulation_dir, "tables", "analysis_run_arguments.tsv")),
    simulation_args_map = safe_read_args(file.path(simulation_dir, "tables", "analysis_run_arguments.tsv")),
    simulation_seed_selection = read_table_optional(file.path(simulation_dir, "tables", "fixed_o2_simulation_representative_seeds.tsv")),
    euler_summary = read_table_optional(file.path(simulation_dir, "tables", "eigen_vs_euler_error_summary.tsv")),
    expm_sim_summary = read_table_optional(file.path(simulation_dir, "tables", "expm_vs_simulation_error_summary.tsv")),
    eigen_sim_summary = read_table_optional(file.path(simulation_dir, "tables", "eigen_vs_simulation_error_summary.tsv")),
    sim_file_status = read_table_optional(file.path(simulation_dir, "tables", "fixed_o2_simulation_file_status.tsv"))
  )
}

write_html_report <- function(analysis_dir, out_dir, report_basename = "index") {
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  asset_dir <- file.path(out_dir, "assets", "figures")
  dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- load_report_data(analysis_dir)
  shared_o2 <- c(0, 0.1, 0.5, 1, 2, 5)
  mode_reference_o2 <- arg_value(dat$top_args_map, "mode_reference_o2", "2")
  mode_reference_counts <- mode_counts_table(dat$attractor_mode_reference)
  mode_all_counts <- mode_counts_table(dat$attractor_mode)
  mode_sentence <- mode_summary_sentence(mode_reference_counts, count_col = "n_seed", unit_label = "seeds")

  attractor_display <- filter_numeric_values(dat$attractor_display, "O2", shared_o2)
  attractor_gap <- filter_numeric_values(dat$attractor_gap, "O2_pct", shared_o2)
  if (is.data.frame(attractor_gap) && "trajectory_regime" %in% names(attractor_gap)) {
    attractor_gap <- attractor_gap[attractor_gap$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  }
  counter_summary <- filter_numeric_values(dat$counter_summary, "O2_pct", shared_o2)
  if (is.data.frame(counter_summary) && "trajectory_regime" %in% names(counter_summary)) {
    counter_summary <- counter_summary[counter_summary$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  }
  counter_key <- counter_summary
  if (is.data.frame(counter_key) && "O2_pct" %in% names(counter_key)) {
    counter_key <- counter_key[counter_key$O2_pct %in% c(1, 2, 5), , drop = FALSE]
  }

  attractor_figures <- list(
    list(
      section = "Attractors",
      relative_path = file.path("attractors", "figures", "fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.pdf"),
      title = "Ploidy, Spectral Gap, and Reliability Composite",
      what_changes = "Fixed O2 is swept along the x-axis. The panels summarize dominant ploidy, spectral gap, relaxation or growth reliability metrics, and mode differences.",
      legend = "For the attractor analysis, mode1 and mode2 are assigned at each fixed O2 value from the corresponding fixed-O2 attractor: mode1 has dominant mean ploidy >= 2, and mode2 has dominant mean ploidy < 2. Higher spectral gap means that the largest eigenvalue is more clearly separated from the next eigenvalue. Grey shaded bands mark O2 settings with weak dominance, defined by median spectral gap < 0.005 or median time to 10-fold dominance > 500 days.",
      result = "The comparison shows how fixed-O2 attractor ploidy and eigenvalue reliability vary after grouping each seed-O2 condition by the attractor-derived dominant ploidy threshold."
    ),
    list(
      section = "Attractors",
      relative_path = file.path("attractors", "figures", "fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.pdf"),
      title = "Seed-Level Distributions of Ploidy and Reliability Metrics",
      what_changes = "Each violin panel shows the seed-level distribution at selected fixed O2 concentrations for mode1 and mode2.",
      legend = "Violin widths show the distribution across seed-O2 conditions. Blue corresponds to mode1, defined by dominant mean ploidy >= 2; orange corresponds to mode2, defined by dominant mean ploidy < 2.",
      result = "These distributions show the extent to which the fixed-O2 attractor and reliability metrics are consistent across fitted seeds rather than being driven only by group medians."
    )
  )
  counter_figures <- list(
    list(
      section = "Counterfactual Trajectories",
      relative_path = file.path("counterfactual_trajectories", "figures", "fixed_o2_counterfactual_median_trajectories.pdf"),
      title = "Fixed-O2 Counterfactual Ploidy Trajectories",
      what_changes = "Each panel fixes one O2 concentration. Within a panel, trajectories are split by initial condition and colored by mode.",
      legend = paste0("Thin semi-transparent curves are seed-level trajectories. Blue is mode1 and orange is mode2, assigned once per seed from the fixed-O2 attractor at O2 = ", html_escape(mode_reference_o2), "."),
      result = "The panel layout compares how fixed oxygen modifies the temporal evolution of mean ploidy from 2N and 4N initial conditions across the two FixO2 attractor-derived modes."
    )
  )
  sim_dir <- file.path(analysis_dir, "simulation", "figures")
  sim_pdf_files <- if (dir.exists(sim_dir)) basename(list.files(sim_dir, pattern = "\\.pdf$", full.names = TRUE)) else character()
  sim_pdf_files <- sim_pdf_files[order(vapply(sim_pdf_files, simulation_figure_order, integer(1)), sim_pdf_files)]
  simulation_figures <- lapply(sim_pdf_files, simulation_figure_spec)
  agreement_figures <- agreement_figure_specs(analysis_dir)
  all_figures <- number_figures(c(attractor_figures, counter_figures, simulation_figures, agreement_figures))
  n_attractor_figures <- length(attractor_figures)
  n_counter_figures <- length(counter_figures)
  n_simulation_figures <- length(simulation_figures)
  n_agreement_figures <- length(agreement_figures)
  pos <- 0L
  attractor_figures <- if (n_attractor_figures) all_figures[seq_len(n_attractor_figures)] else list()
  pos <- pos + n_attractor_figures
  counter_figures <- if (n_counter_figures) {
    all_figures[(pos + 1L):(pos + n_counter_figures)]
  } else {
    list()
  }
  pos <- pos + n_counter_figures
  simulation_figures <- if (n_simulation_figures) {
    all_figures[(pos + 1L):(pos + n_simulation_figures)]
  } else {
    list()
  }
  pos <- pos + n_simulation_figures
  agreement_figures <- if (n_agreement_figures) {
    all_figures[(pos + 1L):(pos + n_agreement_figures)]
  } else {
    list()
  }

  render_figures <- function(figs) {
    if (!length(figs)) return("")
    paste0(vapply(figs, figure_block, character(1), analysis_dir = analysis_dir, asset_dir = asset_dir), collapse = "")
  }

  agreement_methods <- c("eigen", "expm")
  agreement_method_numbers <- stats::setNames(paste0("4.", seq_along(agreement_methods)), agreement_methods)

  render_agreement_method <- function(method, figs) {
    if (!length(figs)) return("")
    summary <- agreement_method_summary(analysis_dir, method)
    label <- agreement_method_label(method)
    numbered_label <- paste(agreement_method_numbers[[method]], label)
    paste0(
      "<section class=\"agreement-method\" id=\"analytical-simulation-agreement-", html_escape(method), "\">",
      "<h3 class=\"method-title\">", html_escape(numbered_label), "</h3>",
      "<p>", html_escape(label), " is shown as a method-specific subsection. The first four figures compare analytical and simulation outputs with scatter and Bland-Altman views, and the final paired figure summarizes fixed-O2 analytical-solution distributions without seed-connecting lines. ",
      html_escape(agreement_summary_sentence(summary)), "</p>",
      render_figures(figs),
      "</section>"
    )
  }

  agreement_method_body <- paste0(vapply(agreement_methods, function(method) {
    figs <- Filter(function(spec) identical(spec$method %||% "", method), agreement_figures)
    render_agreement_method(method, figs)
  }, character(1)), collapse = "")

  sim_complete <- ""
  if (is.data.frame(dat$sim_file_status) && nrow(dat$sim_file_status) && "complete_after" %in% names(dat$sim_file_status)) {
    n_task <- nrow(dat$sim_file_status)
    n_complete <- sum(dat$sim_file_status$complete_after %in% TRUE | tolower(as.character(dat$sim_file_status$complete_after)) == "true", na.rm = TRUE)
    sim_complete <- paste0(n_complete, " of ", n_task, " requested fixed-O2 simulation trajectories were complete after the workflow check.")
  }
  sim_args_map <- dat$simulation_args_map
  sim_mode <- arg_value(sim_args_map, "simulation_mode")
  sim_seed_labels <- compact_csv(arg_value(sim_args_map, "seed_labels", arg_value(sim_args_map, "seed_ids")))
  sim_o2_grid <- compact_csv(arg_value(sim_args_map, "o2_grid"))
  sim_time_grid <- time_grid_summary(arg_value(sim_args_map, "time_grid"))
  sim_dt_grid <- arg_value(sim_args_map, "dt_grid")
  sim_plot_dt <- arg_value(sim_args_map, "plot_dt")
  sim_n_core <- arg_value(sim_args_map, "simulation_n_core")
  sim_worker_threads <- arg_value(sim_args_map, "simulation_worker_threads")
  sim_generate_missing <- arg_value(sim_args_map, "generate_missing_simulation")
  sim_initial_ploidy <- unique_value_summary(dat$sim_file_status, "initial_ploidy", "2, 4")
  sim_initial_condition <- unique_value_summary(dat$sim_file_status, "initial_condition", "init_2N, init_4N")
  sim_replicates <- unique_value_summary(dat$sim_file_status, "simulation_id", compact_csv(arg_value(sim_args_map, "simulation_ids")))
  sim_task_sentence <- ""
  if (is.data.frame(dat$sim_file_status) && nrow(dat$sim_file_status)) {
    sim_task_sentence <- paste0(
      "The expected fixed-O2 simulation set contains ",
      nrow(dat$sim_file_status),
      " seed-by-O2-by-initial-condition-by-replicate state trajectories."
    )
  }
  sim_mode_label <- if (tolower(sim_mode) == "invivo") "in vivo" else sim_mode
  sim_missing_sentence <- if (tolower(sim_generate_missing) %in% c("true", "1", "yes")) {
    paste0(
      "Completed trajectories were reused; absent trajectories were generated as needed and parallelized over ",
      sim_n_core, " cores with ", sim_worker_threads, " worker thread per process."
    )
  } else {
    "Completed trajectories were reused; absent trajectories were not generated automatically in this run."
  }
  sim_seed_selection_sentence <- simulation_seed_selection_sentence(dat$simulation_seed_selection, sim_args_map)

  figure_nav_item <- function(spec) {
    label <- if (identical(spec$type %||% "single", "paired")) {
      suffixes <- vapply(spec$panels, function(panel) panel$suffix %||% "", character(1))
      suffixes <- suffixes[nzchar(suffixes)]
      suffix_range <- if (length(suffixes) > 1L) paste0(suffixes[[1]], "-", suffixes[[length(suffixes)]]) else suffixes[[1]]
      paste0("Figure ", as.integer(spec$figure_number), suffix_range, ". ", spec$title)
    } else {
      paste0("Figure ", as.integer(spec$figure_number), ". ", spec$title)
    }
    list(
      level = "h3",
      id = figure_id(spec),
      title = label
    )
  }
  agreement_nav_items <- list()
  for (method in agreement_methods) {
    method_figs <- Filter(function(spec) identical(spec$method %||% "", method), agreement_figures)
    if (!length(method_figs)) next
    agreement_nav_items <- c(
      agreement_nav_items,
      list(list(level = "h3", id = paste0("analytical-simulation-agreement-", method), title = paste(agreement_method_numbers[[method]], agreement_method_label(method)))),
      lapply(method_figs, figure_nav_item)
    )
  }
  nav_items <- c(
    list(
      list(level = "h2", id = "overview", title = "Overview"),
      list(level = "h2", id = "attractors", title = "1. Attractors")
    ),
    lapply(attractor_figures, figure_nav_item),
    list(
      list(level = "h2", id = "counterfactual-trajectories", title = "2. Counterfactual Trajectories")
    ),
    lapply(counter_figures, figure_nav_item),
    list(
      list(level = "h2", id = "simulation", title = "3. Simulation")
    ),
    lapply(simulation_figures, figure_nav_item),
    list(
      list(level = "h2", id = "analytical-simulation-agreement", title = "4. Analytical-Simulation Agreement Across 500 Seeds")
    ),
    agreement_nav_items,
    list(
      list(level = "h2", id = "run-metadata", title = "Run Metadata")
    )
  )
  nav <- build_nav(nav_items)

  css <- paste0(
    "html{scroll-behavior:smooth;}body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;color:#1f2933;background:#fbfbf8;}",
    ".report-shell{display:flex;gap:22px;align-items:flex-start;padding:24px;}",
    ".sidebar{position:sticky;top:24px;align-self:flex-start;width:292px;max-height:calc(100vh - 48px);overflow:auto;border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);}",
    ".sidebar-header{padding:14px;background:#24384d;color:#fff;}.sidebar-kicker{font-size:11px;text-transform:uppercase;letter-spacing:.06em;opacity:.8;}.sidebar-title{font-size:17px;font-weight:700;margin-top:3px;}",
    ".sidebar nav{padding:10px 8px 14px 8px;}.sidebar ul{margin:0;padding:0;list-style:none;}.sidebar li{margin:3px 0;}.sidebar a{display:block;border-radius:8px;padding:8px 10px;text-decoration:none;color:#17324c;font-size:13px;line-height:1.25;}.sidebar a:hover{background:rgba(47,110,164,0.09);}.sidebar .nav-h3{padding-left:20px;font-size:12px;color:#384c60;}",
    ".content{flex:1;min-width:0;max-width:1320px;}h1{font-size:30px;margin:0 0 4px 0;}h2{margin-top:32px;border-bottom:1px solid #d8d8d0;padding-bottom:7px;scroll-margin-top:24px;}h3{margin:0 0 6px 0;font-size:16px;line-height:1.25;scroll-margin-top:24px;}h4{font-size:14px;margin:18px 0 6px 0;}",
    ".muted{color:#6b7280}.path{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:#475569;}.lead{font-size:14px;line-height:1.55;max-width:980px;}.chapter p,.chapter li{font-size:13px;line-height:1.55;}.chapter ul{padding-left:20px;}",
    ".callout{border-left:4px solid #2f6ea4;background:#f2f6fa;padding:10px 12px;margin:14px 0 18px 0;}.callout p{margin:3px 0;}",
    ".subchapter{margin-top:22px;scroll-margin-top:24px;}.agreement-method{margin-top:18px;}.figure{background:white;border:1px solid #ddd;padding:12px;margin:18px 0;box-shadow:0 1px 2px rgba(0,0,0,0.04);scroll-margin-top:24px;}.paired-figure h3{margin-bottom:10px;}.figure-panel{margin-top:14px;}.figure-panel:first-of-type{margin-top:4px;}.missing-panel{border:1px dashed #f2b8b5;background:#fff7f6;padding:10px;}",
    "img{width:100%;max-width:100%;height:auto;display:block;border:1px solid #eef2f3;}.figure-caption{font-size:12px;line-height:1.45;color:#374151;margin:10px 0 0 0;}.figure-result{font-size:12px;line-height:1.45;color:#374151;margin:7px 0 0 0;}.missing{color:#b91c1c;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}",
    ".table-block{margin:10px 0 20px 0;}.table-caption{font-size:12px;line-height:1.4;color:#374151;margin:0 0 6px 0;}.table-wrap{overflow-x:auto;margin:0 0 18px 0;}table{border-collapse:collapse;width:100%;font-size:12px;background:white;}th,td{border:1px solid #ddd;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#eef2f3;}.table-note{font-size:11px;color:#6b7280;margin-top:-12px;}",
    "code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;background:#f1f5f9;padding:1px 4px;border-radius:4px;}strong{font-weight:700;}",
    "@media (max-width: 1100px){.report-shell{display:block;padding:16px;}.sidebar{position:relative;top:auto;width:auto;max-height:none;margin-bottom:18px;}.content{max-width:none;}}"
  )

  overview <- paste0(
    "<section id=\"overview\" class=\"chapter\">",
    "<h2>Overview</h2>",
    "<p class=\"lead\">This report summarizes the unified fixed-O2 in vivo workflow across attractor analysis, counterfactual fixed-O2 trajectories, simulation validation, and analytical-simulation agreement across the full fitted seed ensemble.</p>",
    "<div class=\"callout\">",
    "<p><strong>Mode1</strong> is the FixO2 attractor-derived label <code>mode1_attractor_dominant_ploidy_ge_2</code>: at the selected reference fixed O2, the dominant attractor has mean ploidy >= 2.</p>",
    "<p><strong>Mode2</strong> is the FixO2 attractor-derived label <code>mode2_attractor_dominant_ploidy_lt_2</code>: at the selected reference fixed O2, the dominant attractor has mean ploidy < 2.</p>",
    "<p>The workflow also outputs the full O2-dependent mode table for every attractor O2 value. Downstream counterfactual and analytical-simulation agreement analyses use the seed-level mode assigned at fixed O2 = ", html_escape(mode_reference_o2), ". These labels are generated inside the FixO2 workflow and are not inherited from the in vivo trajectory classification. ", html_escape(mode_sentence), "</p>",
    "</div>",
    table_block(1, paste0("Reference attractor-derived mode counts by seed at fixed O2 = ", html_escape(mode_reference_o2)), mode_reference_counts, max_rows = 10),
    table_block(2, "Full O2-dependent attractor-derived mode counts by seed-O2 condition", mode_all_counts, max_rows = 10),
    "</section>"
  )

  attractor_body <- paste0(
    "<p>For each fitted seed and each fixed O2 value, the workflow builds a continuous-time linear system over the configured chromosome-count grid <code>N_MIN:N_MAX</code>. The model helper <code>o2pr_build_G()</code> creates the gain and transition operator <code>G(O2)</code>. Each column corresponds to a source ploidy state, and the column entries are the O2-dependent production or missegregation fluxes into destination ploidy states.</p>",
    "<p>The death or loss term is evaluated on the same grid with <code>.mu_eff_of_O2(O2, N)</code>, giving <code>mu_all</code>. The fixed-O2 generator is then <code>M(O2) = G(O2) - diag(mu_all)</code>. The unnormalized population mass follows <code>dx/dt = M x</code>.</p>",
    "<p>The dominant eigenvalue of <code>M</code> gives the asymptotic growth rate under that fixed O2. The corresponding dominant eigenvector is sign-corrected, truncated to nonnegative values, and normalized to sum to one; after normalization it is interpreted as the fixed-O2 attractor composition over ploidy states. Mean ploidy, low-ploidy fractions, selection contrasts between 22/44/88 chromosome states, and spectral gap are then summarized by seed and mode.</p>",
    table_block(3, "Dominant mean ploidy summary at the shared fixed O2 levels", attractor_display, max_rows = 20),
    table_block(
      4,
      "Spectral gap reliability summary for mode1 and mode2",
      attractor_gap,
      max_rows = 20,
      columns = c("O2_pct", "trajectory_regime", "n_seed", "dominant_mean_ploidy_median", "spectral_gap_median", "time_to_10x_days_median", "log10_advantage_1000d_median", "fraction_gap_ge_0p005", "fraction_gap_ge_0p01")
    ),
    render_figures(attractor_figures)
  )

  counter_body <- paste0(
    "<p>Counterfactual trajectory analysis uses the same fixed-O2 matrix <code>M</code> as the attractor analysis, but asks how an initial ploidy composition evolves through time if O2 is held constant. The workflow initializes a point mass at 2N (<code>init_2N</code>) or 4N (<code>init_4N</code>), performs an eigen-decomposition of <code>M</code>, propagates the state over the requested time grid, and normalizes the state to composition at each time point.</p>",
    "<p>For each seed, fixed O2, and initial condition, the script records mean ploidy, low-ploidy fractions, terminal ploidy, minimum ploidy, and crossing times such as the first downward crossing of mean ploidy 1.5. The figure uses the six shared fixed O2 levels: 0, 0.1, 0.5, 1, 2, and 5. Mode grouping in this downstream trajectory analysis is fixed at the seed level using the attractor-derived mode at O2 = ", html_escape(mode_reference_o2), ".</p>",
    table_block(
      5,
      "Summary of terminal counterfactual trajectories at high fixed-O2 levels",
      counter_key,
      max_rows = 24,
      columns = c("O2_pct", "initial_condition", "trajectory_regime", "n_seed", "terminal_mean_ploidy_median", "terminal_mean_ploidy_q25", "terminal_mean_ploidy_q75", "reaches_ploidy_1p5_fraction", "time_crossing_1p5_censored_median", "terminal_minus_dominant_median")
    ),
    render_figures(counter_figures)
  )

  euler_summary <- dat$euler_summary
  expm_sim_summary <- dat$expm_sim_summary
  eigen_sim_summary <- dat$eigen_sim_summary
  simulation_body <- paste0(
    "<p>Simulation validation used the fitted in vivo parameter sets without refitting. The validation cohort used the ", html_escape(sim_mode_label), " simulation mode, seeds ", html_escape(sim_seed_labels), ", fixed O2 values ", html_escape(sim_o2_grid), ", initial ploidies ", html_escape(sim_initial_ploidy), " represented by ", html_escape(sim_initial_condition), ", replicate IDs ", html_escape(sim_replicates), ", and saved observation times from ", html_escape(sim_time_grid), ". ", html_escape(sim_missing_sentence), " The Euler numerical comparison used a step size of ", html_escape(sim_dt_grid), ", and the diagnostic figures use the ", html_escape(sim_plot_dt), " dt setting for the deterministic overlays.</p>",
    "<p>", sim_seed_selection_sentence, "</p>",
    "<p>For each seed, fixed O2 value, initial ploidy, and replicate, the simulator initializes viable mass at the chromosome-count grid state nearest to the requested 2N or 4N ploidy. The in vivo parameterization includes the fitted growth, missegregation, WGD, stress-associated death, and post-missegregation survival parameters (<code>lam_max</code>, <code>p_mis_base</code>, <code>p_misseg</code>, <code>k_o_mis</code>, <code>p_wgd</code>, <code>alpha_o2</code>, <code>gamma_growth</code>, <code>mu_hp</code>, <code>gamma_mu</code>, <code>O2_crit</code>, <code>n_O</code>, <code>buffer_smax</code>, <code>buffer_beta</code>, <code>buffer_n_exp</code>, and the in vivo <code>rho_2N</code>). For the fixed-O2 intervention, the baseline effective oxygen supply parameter <code>o2_S0</code> is overwritten by the assigned fixed O2 value; this override is recorded in the parameter audit table as <code>fixed_o2_override</code>.</p>",
    "<p>In the resource-stress model, oxygen is an endogenous state variable driven by tumor burden, with target oxygen defined as <code>O2_target = max(o2_min, o2_S0 - kappa_O * log(1 + N_eff / o2_Nref))</code> and effective oxygen allowed to relax toward that target. The fixed-O2 simulation removes this feedback by disabling burden-dependent oxygen dynamics (<code>cfg$o2_burden_feedback = FALSE</code> and <code>o2_feedback = FALSE</code>) and setting both target and effective oxygen to the assigned fixed value at every recorded time point. Thus fixed O2 replaces the model-generated oxygen trajectory, while all O2-dependent biological rates are still evaluated at that constant oxygen level.</p>",
    "<p>At each saved time point, the simulator exports two complementary summaries. The population-level trajectory records fixed, target, and effective O2 together with viable, hypoxia-origin dead, CIN-associated dead, total cell counts, and the corresponding tumor volumes. The state-resolved trajectory records the viable, hypoxia-origin dead, and CIN-associated dead cell counts for every chromosome-count and ploidy state, along with the state-specific growth, stress-associated death, missegregation, WGD, CIN-associated non-viable daughter production, boundary loss, and net viable rates. The validation reads the viable state distribution to compute replicate-specific mean chromosome number, mean ploidy, ploidy standard deviation, and viable-cell totals; these replicate summaries are then averaged by seed, O2, initial condition, and day for comparison with the analytical solutions and for the mean-SD phase-plane plots. Completed state-resolved trajectories are reused, whereas missing trajectories are generated before the comparison step. ", html_escape(sim_task_sentence), "</p>",
    "<p><strong>Expm analytical solution:</strong> The matrix-exponential solution treats <code>M</code> as the full fixed-O2 generator and applies the finite-time propagator directly to the initial state. For each observation interval, the state is updated as <code>x(t_i) = exp(M * delta_t) x(t_{i-1})</code>, rescaled to avoid numerical overflow, and normalized to a composition before summary statistics are calculated. This approach preserves the combined effect of all transient modes without explicitly diagonalizing <code>M</code>, and it is used as the primary deterministic reference for finite-time trajectories.</p>",
    "<p><strong>Eigen analytical solution:</strong> The eigen solution first decomposes the generator into its modal representation, <code>M = V Lambda V^{-1}</code>. The initial state is projected onto the eigenvector basis by solving <code>c = V^{-1} x(0)</code>, each mode is propagated by <code>exp((lambda - lambda_ref) t)</code>, and the state is reconstructed from the weighted eigenvectors before taking the real, nonnegative, normalized composition. The subtraction of <code>lambda_ref</code>, the largest real eigenvalue, removes the common exponential growth scale so that the calculation focuses on composition dynamics. This representation is informative because it separates dominant attractor modes from transient relaxation modes, but it is more sensitive than expm to eigenvector conditioning, near-degenerate eigenvalues, complex components, and nonnegative truncation.</p>",
    "<p><strong>Euler numerical solution:</strong> Euler integration advances the same ODE with the finite-step update <code>x_{k+1} = x_k + dt M x_k</code>. It is not a separate biological model; it is a numerical check that the finite-step integration used in trajectory diagnostics follows the direct matrix-exponential solution.</p>",
    if (nzchar(sim_complete)) paste0("<div class=\"callout\"><p>", html_escape(sim_complete), "</p></div>") else "",
    table_block(
      6,
      "Representative seed selection for simulation validation",
      dat$simulation_seed_selection,
      max_rows = 8,
      columns = c("selection_rank", "selection_mode", "mode_label", "seed_id", "final_objective", "objective_source", "mode_reference_o2_pct", "mode_reference_dominant_mean_ploidy", "selection_rule")
    ),
    table_block(
      7,
      "Expm analytical solution versus Euler numerical integration error summary",
      euler_summary,
      max_rows = 24,
      columns = c("seed_id", "seed_mode", "O2_pct", "initial_condition", "dt", "max_abs_diff_euler_expm_mean_ploidy", "terminal_abs_diff_euler_expm_mean_ploidy", "max_abs_diff_eigen_expm_mean_ploidy", "terminal_abs_diff_eigen_expm_mean_ploidy")
    ),
    table_block(
      8,
      "Expm analytical solution versus simulation replicate error summary",
      expm_sim_summary,
      max_rows = 24,
      columns = c("seed_id", "seed_label", "O2_pct", "initial_condition", "simulation_n", "max_abs_diff_simulation_analytical_mean_ploidy", "terminal_abs_diff_simulation_analytical_mean_ploidy", "max_abs_diff_simulation_analytical_sd_ploidy")
    ),
    table_block(
      9,
      "Eigen analytical solution versus simulation replicate error summary",
      eigen_sim_summary,
      max_rows = 24,
      columns = c("seed_id", "seed_label", "O2_pct", "initial_condition", "simulation_n", "max_abs_diff_simulation_analytical_mean_ploidy", "terminal_abs_diff_simulation_analytical_mean_ploidy", "max_abs_diff_simulation_analytical_sd_ploidy")
    ),
    render_figures(simulation_figures)
  )

  agreement_body <- paste0(
    "<p>This analysis extends the representative seed simulation validation to the full 500-seed fitted ensemble. For each seed, fixed O2 value, initial condition, and selected time point, the deterministic analytical mean ploidy is paired with the mean ploidy inferred from fixed-O2 simulation replicates. The paired analysis is reported separately for the eigen and expm analytical solutions, with the eigen block shown first and the expm block below it.</p>",
    "<p>Mode-colored agreement panels use the same seed-level reference mode as the counterfactual trajectory analysis, assigned from the fixed-O2 attractor at O2 = ", html_escape(mode_reference_o2), ".</p>",
    "<p>The Bland-Altman panels are constructed from the same matched rows as the scatter panels. For each row, the x-axis is the average of the analytical and simulation mean ploidy values, <code>(analytical + simulation) / 2</code>, and the y-axis is the residual <code>simulation - analytical</code>. The horizontal bias line is the mean residual, and the dashed limits of agreement are <code>bias +/- 1.96 * SD(residual)</code>. These panels therefore emphasize systematic bias and residual spread even when the scatter points lie close to the equality diagonal.</p>",
    agreement_method_body
  )

  metadata_body <- paste0(
    table_block(10, "Unified workflow arguments", dat$top_args, max_rows = 30),
    table_block(11, "Counterfactual trajectory arguments", dat$counter_args, max_rows = 30),
    table_block(12, "Simulation validation arguments", dat$simulation_args, max_rows = 30)
  )

  html <- paste0(
    "<!doctype html>\n<html><head><meta charset=\"utf-8\">",
    "<title>FixO2 In Vivo Report</title>",
    "<style>", css, "</style>",
    "</head><body><div class=\"report-shell\">",
    "<aside class=\"sidebar\"><div class=\"sidebar-header\"><div class=\"sidebar-kicker\">Navigation</div><div class=\"sidebar-title\">FixO2 In Vivo</div></div><nav>",
    nav,
    "</nav></aside><main class=\"content\">",
    "<h1>FixO2 In Vivo Report</h1>",
    "<p class=\"path\">", html_escape(analysis_dir), "</p>",
    "<p class=\"muted\">Generated: ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p>",
    overview,
    chapter_header("attractors", "1. Attractors", attractor_body),
    chapter_header("counterfactual-trajectories", "2. Counterfactual Trajectories", counter_body),
    chapter_header("simulation", "3. Simulation", simulation_body),
    chapter_header("analytical-simulation-agreement", "4. Analytical-Simulation Agreement Across 500 Seeds", agreement_body),
    chapter_header("run-metadata", "Run Metadata", metadata_body),
    "</main></div></body></html>"
  )

  out_path <- file.path(out_dir, paste0(report_basename, ".html"))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, con = out_path, useBytes = TRUE)
  out_path
}

main <- function(argv = parse_args()) {
  root <- repo_root()
  analysis_dir <- argv$analysis_dir %||% argv$results_dir %||% file.path(root, "oxygen", "results", "analysis", "FixO2_invivo_500seed")
  analysis_dir <- normalizePath(path.expand(analysis_dir), mustWork = TRUE)
  out_dir <- argv$out_dir %||% file.path(analysis_dir, "html_report")
  report_basename <- argv$report_basename %||% "index"
  out_path <- write_html_report(analysis_dir = analysis_dir, out_dir = out_dir, report_basename = report_basename)
  message("FixO2 in vivo HTML report written to: ", normalizePath(out_path, mustWork = FALSE))
  invisible(out_path)
}

if (identical(environment(), globalenv())) {
  main()
}
