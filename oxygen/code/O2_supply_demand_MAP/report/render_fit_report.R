#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
}

parse_args <- function(args) {
  out <- list()
  if (length(args) == 0) return(out)
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

REPORT_SCRIPT_DIR <- normalizePath(get_script_dir(), mustWork = FALSE)
REPORT_WORKFLOW_ROOT <- normalizePath(file.path(REPORT_SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(REPORT_WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())

is_fit_dir <- function(path) {
  dir.exists(path) &&
    file.exists(file.path(path, "fit_summary.tsv")) &&
    file.exists(file.path(path, "best_params.tsv")) &&
    dir.exists(file.path(path, "viz"))
}

find_fit_dirs_under <- function(root_dir) {
  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  sub_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  fit_sub_dirs <- sub_dirs[vapply(sub_dirs, is_fit_dir, logical(1))]
  if (is_fit_dir(root_dir)) {
    fit_sub_dirs <- c(root_dir, fit_sub_dirs)
  }
  sort(unique(normalizePath(fit_sub_dirs, mustWork = TRUE)))
}

fit_report_seed_token <- function(fit_dir) {
  base <- basename(normalizePath(fit_dir, mustWork = FALSE))
  hit <- regmatches(base, regexpr("seed[0-9]+", base, ignore.case = TRUE))
  if (length(hit) > 0L && nzchar(hit[[1]])) {
    return(tolower(hit[[1]]))
  }

  summary_path <- file.path(fit_dir, "fit_summary.tsv")
  if (file.exists(summary_path)) {
    summary_tab <- tryCatch(
      utils::read.delim(
        summary_path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = c("metric" = "character", "value" = "character")
      ),
      error = function(e) NULL
    )
    if (is.data.frame(summary_tab) && all(c("metric", "value") %in% names(summary_tab))) {
      seed_value <- summary_tab$value[match("seed", summary_tab$metric)]
      seed_value <- suppressWarnings(as.integer(seed_value[[1]] %||% NA_integer_))
      if (is.finite(seed_value)) return(paste0("seed", seed_value))
    }
  }

  NA_character_
}

fit_report_basename_with_seed <- function(report_basename, fit_dir) {
  seed_token <- fit_report_seed_token(fit_dir)
  if (is.na(seed_token) || !nzchar(seed_token)) return(report_basename)
  if (grepl(seed_token, report_basename, fixed = TRUE)) return(report_basename)
  paste(report_basename, seed_token, sep = "_")
}

read_fit_summary_selected <- function(fit_dir) {
  path <- file.path(fit_dir, "fit_summary.tsv")
  tab <- read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c("metric" = "character", "value" = "character")
  )
  fit_mode <- tab$value[match("fit_mode", tab$metric)]
  wanted <- if (length(fit_mode) > 0L && identical(as.character(fit_mode[[1]]), "fit_joint")) {
    c(
      "fit_mode",
      "optimizer_method",
      "objective",
      "objective_invivo",
      "objective_invitro",
      "objective_soft_coupling",
      "objective_constraints",
      "joint_weight_invivo",
      "joint_weight_invitro",
      "joint_invitro_growth_weight",
      "joint_invitro_ploidy_weight",
      "joint_invitro_flow_weight",
      "joint_soft_coupling_enabled",
      "joint_soft_coupling_params",
      "n_cores_requested",
      "n_cores_used",
      "n_parameters"
    )
  } else {
    c(
      "optimizer_method",
      "objective",
      "objective_data",
      "objective_prior_raw",
      "objective_prior",
      "objective_burden",
      "objective_ploidy",
      "objective_burden_neg2loglik_raw",
      "objective_ploidy_neg2loglik_raw"
    )
  }
  idx <- match(wanted, tab$metric)
  data.frame(
    metric = wanted,
    value = ifelse(is.na(idx), NA_character_, tab$value[idx]),
    stringsAsFactors = FALSE
  )
}

read_fit_summary_map <- function(fit_dir) {
  path <- file.path(fit_dir, "fit_summary.tsv")
  tab <- read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c("metric" = "character", "value" = "character")
  )
  setNames(tab$value, tab$metric)
}

summary_flag_true <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  tolower(trimws(as.character(x[[1]]))) %in% c("true", "t", "1", "yes", "y", "on")
}

filter_best_params_for_report <- function(best_params, fit_summary_map) {
  drop_names <- character(0)
  best_params[!(best_params$parameter %in% unique(drop_names)), , drop = FALSE]
}

parameter_description_table_paths <- function(fit_dir) {
  unique(c(
    file.path(fit_dir, "parameter_table_input.csv"),
    file.path(fit_dir, "parameter_table.csv"),
    file.path(fit_dir, "parameter_table_input_invivo.csv"),
    file.path(fit_dir, "parameter_table_input_invitro.csv"),
    file.path(fit_dir, "parameter_table_invivo_transformed.csv"),
    file.path(fit_dir, "parameter_table_invitro_transformed.csv")
  ))
}

annotate_parameter_table_for_report <- function(tab, fit_dir) {
  o2sd_add_parameter_descriptions(
    tab,
    parameter_tables = parameter_description_table_paths(fit_dir)
  )
}

read_best_params <- function(fit_dir) {
  path <- file.path(fit_dir, "best_params.tsv")
  tab <- read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c("parameter" = "character", "value" = "character")
  )
  annotate_parameter_table_for_report(
    filter_best_params_for_report(tab, read_fit_summary_map(fit_dir)),
    fit_dir = fit_dir
  )
}

read_report_table_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (identical(tolower(tools::file_ext(path)), "csv")) {
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}

read_soft_coupling_table_for_report <- function(fit_dir) {
  tab <- read_report_table_optional(file.path(fit_dir, "joint_soft_coupling.tsv"))
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(NULL)
  keep <- c(
    "parameter",
    "center_natural",
    "vivo_natural",
    "vitro_natural",
    "delta_transformed",
    "delta_interpretable",
    "log10_ratio_vivo_to_vitro",
    "fold_change_vivo_to_vitro",
    "ratio_vivo_to_vitro",
    "logit_difference_vivo_to_vitro",
    "odds_ratio_vivo_to_vitro",
    "regularization_sigma",
    "penalty_paid",
    "joint_union_lower_transformed",
    "joint_union_upper_transformed",
    "feasible_at_solution"
  )
  keep <- intersect(keep, names(tab))
  tab[, keep, drop = FALSE]
}

report_truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y", "on")
}

transformed_param_to_natural <- function(x) {
  x <- as.character(x)
  x <- sub("^ivt__", "", x)
  x <- sub("^log10_", "", x)
  x <- sub("^logit_", "", x)
  x
}

read_invivo_fitted_natural_names <- function(fit_dir) {
  tab <- read_report_table_optional(file.path(fit_dir, "parameter_table_invivo_transformed.csv"))
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(character(0))
  if (all(c("estimate", "param_prototype") %in% names(tab))) {
    return(unique(as.character(tab$param_prototype[report_truthy(tab$estimate)])))
  }
  name_col <- intersect(c("param_name", "transformed_parameter"), names(tab))
  if (length(name_col) == 0L) return(character(0))
  unique(transformed_param_to_natural(tab[[name_col[[1]]]]))
}

read_invitro_fitted_natural_names <- function(fit_dir) {
  natural_tab <- read_report_table_optional(file.path(fit_dir, "parameter_table_input_invitro.csv"))
  if (is.data.frame(natural_tab) && nrow(natural_tab) > 0L &&
      all(c("param_symbol", "use_invitro_fit") %in% names(natural_tab))) {
    return(unique(as.character(natural_tab$param_symbol[report_truthy(natural_tab$use_invitro_fit)])))
  }

  transformed_tab <- read_report_table_optional(file.path(fit_dir, "parameter_table_invitro_transformed.csv"))
  if (!is.data.frame(transformed_tab) || nrow(transformed_tab) == 0L ||
      !"param_name" %in% names(transformed_tab)) {
    return(character(0))
  }
  unique(transformed_param_to_natural(transformed_tab$param_name))
}

split_joint_scope_parameters <- function(tab, fitted_names, shared_names) {
  if (!is.data.frame(tab) || nrow(tab) == 0L || !"parameter" %in% names(tab)) {
    return(list(fitted = NULL, fixed = NULL))
  }
  tab$parameter <- as.character(tab$parameter)
  tab <- tab[!(tab$parameter %in% shared_names), , drop = FALSE]
  fit_key <- sub("_(vivo|vitro)$", "", tab$parameter)
  fitted_names <- unique(as.character(fitted_names))
  list(
    fitted = tab[fit_key %in% fitted_names, , drop = FALSE],
    fixed = tab[!(fit_key %in% fitted_names), , drop = FALSE]
  )
}

joint_shared_value_for_scope <- function(shared_row, scope) {
  preferred <- if (identical(scope, "invivo")) {
    c("in_vivo_value", "shared_value", "value")
  } else {
    c("in_vitro_value", "shared_value", "value")
  }
  for (col in preferred) {
    if (col %in% names(shared_row)) return(shared_row[[col]][[1]])
  }
  NA
}

append_joint_only_parameter_row <- function(tab, parameter, value) {
  row <- data.frame(parameter = parameter, value = value, stringsAsFactors = FALSE)
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(row)
  if ("parameter" %in% names(tab) && parameter %in% as.character(tab$parameter)) return(tab)

  missing_in_row <- setdiff(names(tab), names(row))
  for (col in missing_in_row) row[[col]] <- NA
  missing_in_tab <- setdiff(names(row), names(tab))
  for (col in missing_in_tab) tab[[col]] <- NA
  row <- row[, names(tab), drop = FALSE]
  rbind(tab, row)
}

move_joint_o2_crit_from_shared <- function(shared_tab, invivo_tab, invitro_tab) {
  if (!is.data.frame(shared_tab) || !"parameter" %in% names(shared_tab)) {
    return(list(shared = shared_tab, invivo = invivo_tab, invitro = invitro_tab))
  }
  shared_tab$parameter <- as.character(shared_tab$parameter)
  o2_rows <- shared_tab[shared_tab$parameter == "O2_crit", , drop = FALSE]
  shared_tab <- shared_tab[shared_tab$parameter != "O2_crit", , drop = FALSE]
  if (nrow(o2_rows) > 0L) {
    o2_row <- o2_rows[1L, , drop = FALSE]
    invivo_tab <- append_joint_only_parameter_row(
      invivo_tab,
      "O2_crit_vivo",
      joint_shared_value_for_scope(o2_row, "invivo")
    )
    invitro_tab <- append_joint_only_parameter_row(
      invitro_tab,
      "O2_crit_vitro",
      joint_shared_value_for_scope(o2_row, "invitro")
    )
  }
  list(shared = shared_tab, invivo = invivo_tab, invitro = invitro_tab)
}

format_report_number <- function(x) {
  num <- as.numeric(x)
  out <- as.character(x)
  finite <- is.finite(num)
  if (!any(finite)) return(out)
  num_finite <- num[finite]
  int_like <- abs(num_finite - round(num_finite)) < 1e-9
  decimal_text <- formatC(num_finite, format = "f", digits = 3)
  sci_nonzero <- !int_like & num_finite != 0 & grepl("\\.000$", decimal_text)
  decimal <- !int_like & !sci_nonzero
  out_finite <- character(length(num_finite))
  out_finite[int_like] <- format(round(num_finite[int_like]), scientific = FALSE, trim = TRUE, digits = 15)
  out_finite[sci_nonzero] <- formatC(num_finite[sci_nonzero], format = "e", digits = 3)
  out_finite[decimal] <- formatC(num_finite[decimal], format = "f", digits = 3)
  out[finite] <- out_finite
  out
}

format_report_value <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  if (is.numeric(x) || is.integer(x)) {
    return(format_report_number(x))
  }
  out <- as.character(x)
  x_trim <- trimws(out)
  numeric_pattern <- "^[-+]?((\\d+\\.?\\d*)|(\\.\\d+))([eE][-+]?\\d+)?$"
  numeric_like <- !is.na(out) & nzchar(x_trim) & grepl(numeric_pattern, x_trim)
  if (any(numeric_like)) {
    num <- suppressWarnings(as.numeric(x_trim[numeric_like]))
    keep <- is.finite(num)
    if (any(keep)) {
      idx <- which(numeric_like)[keep]
      out[idx] <- format_report_number(num[keep])
    }
  }
  out
}

format_report_table <- function(tab) {
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0L) return(tab)
  out <- as.data.frame(lapply(tab, format_report_value), stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- names(tab)
  out
}

read_parameter_sections <- function(fit_dir) {
  fit_summary <- read_fit_summary_map(fit_dir)
  if (!identical(as.character(fit_summary[["fit_mode"]]), "fit_joint")) {
    return(list(list(name = "Best Parameters", table = read_best_params(fit_dir))))
  }
  shared_tab <- read_report_table_optional(file.path(fit_dir, "joint_params_shared.tsv"))
  soft_tab <- read_soft_coupling_table_for_report(fit_dir)
  invivo_tab <- read_report_table_optional(file.path(fit_dir, "joint_params_invivo_only.tsv"))
  invitro_tab <- read_report_table_optional(file.path(fit_dir, "joint_params_invitro_only.tsv"))
  moved_tabs <- move_joint_o2_crit_from_shared(shared_tab, invivo_tab, invitro_tab)
  shared_tab <- moved_tabs$shared
  invivo_tab <- moved_tabs$invivo
  invitro_tab <- moved_tabs$invitro
  shared_names <- if (is.data.frame(shared_tab) && "parameter" %in% names(shared_tab)) {
    unique(as.character(shared_tab$parameter))
  } else {
    character(0)
  }
  invivo_split <- split_joint_scope_parameters(
    invivo_tab,
    fitted_names = read_invivo_fitted_natural_names(fit_dir),
    shared_names = shared_names
  )
  invitro_split <- split_joint_scope_parameters(
    invitro_tab,
    fitted_names = read_invitro_fitted_natural_names(fit_dir),
    shared_names = shared_names
  )
  sections <- list(
    list(group = "Fitted Parameters", name = "Soft-Coupled Parameters", table = soft_tab),
    list(group = "Fitted Parameters", name = "Shared Parameters", table = annotate_parameter_table_for_report(shared_tab, fit_dir)),
    list(group = "Fitted Parameters", name = "In Vivo-Only Fitted Parameters", table = annotate_parameter_table_for_report(invivo_split$fitted, fit_dir)),
    list(group = "Fitted Parameters", name = "In Vitro-Only Fitted Parameters", table = annotate_parameter_table_for_report(invitro_split$fitted, fit_dir)),
    list(group = "Fixed Parameters", name = "In Vivo-Only Fixed Parameters", table = annotate_parameter_table_for_report(invivo_split$fixed, fit_dir)),
    list(group = "Fixed Parameters", name = "In Vitro-Only Fixed Parameters", table = annotate_parameter_table_for_report(invitro_split$fixed, fit_dir))
  )
  sections <- Filter(function(section) is.data.frame(section$table) && nrow(section$table) > 0L, sections)
  if (length(sections) == 0L) {
    sections <- list(list(name = "Best Parameters", table = read_best_params(fit_dir)))
  }
  sections
}

extract_horizon_day <- function(path) {
  base <- basename(path)
  as.integer(sub(".*_0_([0-9]+)day\\.pdf$", "\\1", base))
}

sort_paths_by_horizon <- function(paths) {
  if (length(paths) == 0) return(paths)
  paths[order(vapply(paths, extract_horizon_day, integer(1)))]
}

report_pandoc_available <- function() {
  isTRUE(rmarkdown::pandoc_available("1.12.3"))
}

report_pdflatex_available <- function() {
  if (identical(Sys.getenv("O2_REPORT_FORCE_NO_PDFLATEX", unset = ""), "TRUE")) {
    return(FALSE)
  }
  nzchar(Sys.which("pdflatex"))
}

report_pdf_enabled <- function(argv = list()) {
  arg_value <- argv$render_pdf %||% Sys.getenv("O2_RENDER_PDF", unset = "")
  if (!nzchar(arg_value)) {
    return(FALSE)
  }
  tolower(arg_value) %in% c("true", "t", "1", "yes", "y")
}

escape_html <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

html_id_slug <- function(x) {
  slug <- tolower(trimws(as.character(x %||% "section")))
  slug <- gsub("[^a-z0-9]+", "-", slug)
  slug <- gsub("(^-+|-+$)", "", slug)
  if (!nzchar(slug)) slug <- "section"
  slug
}

parameter_group_id <- function(group_index, group) {
  sprintf("params-group-%02d-%s", group_index, html_id_slug(group))
}

parameter_section_id <- function(index, name) {
  sprintf("params-table-%02d-%s", index, html_id_slug(name))
}

figure_part_id <- function(section_index, part) {
  part_index <- suppressWarnings(as.integer(part$part_index %||% NA_integer_))
  if (!is.finite(part_index) || part_index < 1L) part_index <- section_index
  part_title <- part$title %||% paste0("Part ", part_index)
  sprintf("figures-%02d-part-%02d-%s", section_index, part_index, html_id_slug(part_title))
}

figure_subpart_id <- function(section_index, part, subpart) {
  part_index <- suppressWarnings(as.integer(part$part_index %||% NA_integer_))
  if (!is.finite(part_index) || part_index < 1L) part_index <- section_index
  subpart_index <- suppressWarnings(as.integer(subpart$subpart_index %||% NA_integer_))
  if (!is.finite(subpart_index) || subpart_index < 1L) subpart_index <- 1L
  subpart_title <- subpart$title %||% paste0("Subpart ", subpart_index)
  sprintf(
    "figures-%02d-part-%02d-subpart-%02d-%s",
    section_index,
    part_index,
    subpart_index,
    html_id_slug(subpart_title)
  )
}

fresh_png_companion_for_pdf <- function(path) {
  png_companion <- sub("\\.pdf$", ".png", path, ignore.case = TRUE)
  if (!file.exists(png_companion)) return(NULL)
  src_info <- file.info(path)
  png_info <- file.info(png_companion)
  if (is.na(src_info$mtime[[1]]) || is.na(png_info$mtime[[1]])) return(NULL)
  if (png_info$mtime[[1]] < src_info$mtime[[1]]) return(NULL)
  normalizePath(png_companion, mustWork = TRUE)
}

make_figure_spec <- function(path, title, legend) {
  default_title <- tools::file_path_sans_ext(basename(path))
  title_use <- if (is.null(title) || !length(title) || is.na(title[[1]]) || !nzchar(trimws(title[[1]]))) {
    default_title
  } else {
    as.character(title[[1]])
  }
  legend_use <- if (is.null(legend) || !length(legend) || is.na(legend[[1]]) || !nzchar(trimws(legend[[1]]))) {
    paste0("Figure source: ", basename(path), ".")
  } else {
    as.character(legend[[1]])
  }
  list(
    src = normalizePath(path, mustWork = TRUE),
    html_src = fresh_png_companion_for_pdf(path),
    title = title_use,
    legend = legend_use
  )
}

optional_figure <- function(viz_dir, filename, title, legend) {
  path <- file.path(viz_dir, filename)
  if (!file.exists(path)) return(list())
  list(make_figure_spec(path, title, legend))
}

optional_series_figures <- function(paths, title_tpl, legend_tpl) {
  if (length(paths) == 0L) return(list())
  lapply(paths, function(path) {
    hz <- extract_horizon_day(path)
    make_figure_spec(
      path,
      sprintf(title_tpl, hz),
      sprintf(legend_tpl, hz)
    )
  })
}

optional_predicted_curve_series_figures <- function(curve_paths, death_language, annotated_path = NULL) {
  if (length(curve_paths) == 0L) return(list())
  figs <- list()
  for (curve_path in curve_paths) {
    hz <- extract_horizon_day(curve_path)
    figs <- c(figs, list(make_figure_spec(
      curve_path,
      sprintf("Predicted Curves (0-%s day)", hz),
      sprintf(
        paste0(
          "Forward simulation from day 0 to %s summarizing predicted burden on a log10 scale, ",
          "predicted live cells on a log10 scale, the %s death fraction among all deaths, ",
          "ploidy, chromosome-state density, and O2 trajectories."
        ),
        hz,
        death_language$figure_phrase
      )
    )))
    if (identical(hz, 1000L) && !is.null(annotated_path) && file.exists(annotated_path)) {
      figs <- c(figs, list(make_figure_spec(
        annotated_path,
        "Predicted (0-1000 day)",
        paste0(
          "Forward simulation from day 0 to 1000 using the chromosome-N probability heatmap as the main panel. ",
          "Column annotations show cohort-level burden, live cells, and O2% in 2N/4N rows; ",
          "the mean chromosome-number curve is aligned below the heatmap."
        )
      )))
    }
  }
  figs
}

render_pdf_preview_png <- function(src_pdf, dest_png, density = 180) {
  img <- magick::image_read(src_pdf, density = density)
  if (length(img) > 1L) {
    img <- img[1]
  }
  magick::image_write(img, path = dest_png, format = "png")
  normalizePath(dest_png, mustWork = TRUE)
}

report_magick_available <- function() {
  if (identical(Sys.getenv("O2_REPORT_FORCE_NO_MAGICK", unset = ""), "TRUE")) {
    return(FALSE)
  }
  requireNamespace("magick", quietly = TRUE)
}

report_gs_available <- function() {
  nzchar(Sys.which("gs"))
}

report_base64enc_available <- function() {
  requireNamespace("base64enc", quietly = TRUE)
}

render_pdf_preview_png_gs <- function(src_pdf, dest_png, density = 180) {
  gs_bin <- Sys.which("gs")
  if (!nzchar(gs_bin)) {
    stop("Ghostscript ('gs') was requested for PDF preview rendering but is not available in PATH.")
  }
  src_pdf_use <- normalizePath(src_pdf, mustWork = TRUE)
  dest_png_use <- normalizePath(dest_png, mustWork = FALSE)
  density_use <- suppressWarnings(as.integer(density))
  if (!is.finite(density_use) || density_use <= 0L) density_use <- 180L
  args <- c(
    "-dSAFER",
    "-dBATCH",
    "-dNOPAUSE",
    "-sDEVICE=pngalpha",
    sprintf("-r%d", density_use),
    sprintf("-sOutputFile=%s", shQuote(dest_png_use)),
    shQuote(src_pdf_use)
  )
  out <- suppressWarnings(system2(gs_bin, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "Ghostscript failed while rendering PDF preview for ", src_pdf, ": ",
      paste(out, collapse = "\n")
    )
  }
  if (!file.exists(dest_png_use)) {
    stop("Ghostscript did not create expected PNG preview: ", dest_png_use)
  }
  normalizePath(dest_png_use, mustWork = TRUE)
}

file_to_data_uri <- function(path, mime) {
  if (report_base64enc_available()) {
    return(base64enc::dataURI(file = path, mime = mime))
  }
  base64_bin <- Sys.which("base64")
  if (nzchar(base64_bin)) {
    enc <- tryCatch(
      suppressWarnings(system2(base64_bin, c("-w", "0", path), stdout = TRUE, stderr = TRUE)),
      error = function(e) character()
    )
    if (!length(enc)) {
      enc <- tryCatch(
        suppressWarnings(system2(base64_bin, path, stdout = TRUE, stderr = TRUE)),
        error = function(e) character()
      )
    }
    if (length(enc) > 0L) {
      return(sprintf("data:%s;base64,%s", mime, paste(enc, collapse = "")))
    }
  }
  stop(
    "HTML report fallback requires either the R package 'base64enc' or a system 'base64' command ",
    "when 'magick' is unavailable."
  )
}

pdf_to_data_uri <- function(pdf_path) {
  file_to_data_uri(pdf_path, mime = "application/pdf")
}

resource_death_language_report <- function() {
  list(
    report_phrase = "hypoxia-dead",
    figure_phrase = "hypoxia"
  )
}

is_joint_fit_dir_for_report <- function(fit_dir) {
  fit_summary <- read_fit_summary_map(fit_dir)
  identical(as.character(fit_summary[["fit_mode"]]), "fit_joint")
}

resolve_invivo_report_viz_dir <- function(fit_dir) {
  viz_root <- file.path(fit_dir, "viz")
  invivo_viz <- file.path(viz_root, "invivo")
  if (is_joint_fit_dir_for_report(fit_dir) && dir.exists(invivo_viz)) {
    return(invivo_viz)
  }
  viz_root
}

build_invivo_section_specs <- function(fit_dir) {
  viz_dir <- resolve_invivo_report_viz_dir(fit_dir)
  death_language <- resource_death_language_report()

  burden_predict <- sort_paths_by_horizon(list.files(
    viz_dir,
    pattern = "^predict_burden_live_dead_decomposition_0_[0-9]+day\\.pdf$",
    full.names = TRUE
  ))
  burden_predict_figs <- optional_figure(
    viz_dir,
    "predict_burden_live_dead_decomposition_combined.pdf",
    "Predicted Burden Live/Dead Decomposition (0-100/300/1000 day)",
    paste0(
      "Forward simulations from day 0 to 100, 300, and 1000 showing the predicted live/",
      death_language$figure_phrase,
      "/buffer burden decomposition with y-axis ticks shown as log10 burden values and a shared component legend."
    )
  )
  if (length(burden_predict_figs) == 0L && length(burden_predict) > 0L) {
    burden_predict_figs <- optional_series_figures(
      burden_predict,
      "Predicted Burden Live/Dead Decomposition (0-%s day)",
      paste0(
        "Forward simulation from day 0 to %s showing the predicted live/",
        death_language$figure_phrase,
        "/buffer burden decomposition with y-axis ticks shown as log10 burden values."
      )
    )
  }
  ploidy_predict <- sort_paths_by_horizon(list.files(
    viz_dir,
    pattern = "^predict_curves_0_[0-9]+day\\.pdf$",
    full.names = TRUE
  ))
  annotated_predict_1000 <- file.path(viz_dir, "predicted_0_1000day.pdf")

  burden_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "burden_trend_absolute(real_scale).pdf",
      "Burden Trend Absolute (Real Scale)",
      "Absolute tumor burden shown on the real scale across observation times."
    ),
    optional_figure(
      viz_dir,
      "burden_live_dead_decomposition.pdf",
      "Burden Live/Dead Decomposition",
      paste0(
        "Observed and fitted burden decomposed into live, ",
        death_language$report_phrase,
        ", and buffer-dead components."
      )
    ),
    burden_predict_figs
  ))

  ploidy_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "ploidy_weighted_mean_over_time.pdf",
      "Weighted Mean Ploidy Over Time",
      "Weighted mean chromosome-state trajectory over time under the fitted model."
    ),
    optional_figure(
      viz_dir,
      "terminal_ploidy_observed_vs_predicted_violin.pdf",
      "Terminal Ploidy Observed vs Predicted",
      "Observed and predicted terminal ploidy or chromosome-number distributions at harvest."
    ),
    optional_predicted_curve_series_figures(
      ploidy_predict,
      death_language,
      annotated_path = annotated_predict_1000
    )
  ))

  oxygen_dynamics_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "o2_target_vs_eff_timecourse.pdf",
      "O2 Target vs Effective Timecourse",
      "Timecourse comparison between the oxygen target and the lagged effective oxygen state used by the model."
    ),
    optional_figure(
      viz_dir,
      "predict_burden_vs_o2.pdf",
      "Predicted Burden vs O2",
      "Forward-simulation burden trajectories plotted against the effective oxygen state."
    )
  ))
  oxygen_ms_relationship_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "ms_rate_vs_nonviable_daughter_fraction.pdf",
      "Nonviable Daughter Fraction vs MS Rate",
      "Per-division fraction of daughter cells that are nonviable because of missegregation-linked loss, shown against missegregation rate."
    ),
    optional_figure(
      viz_dir,
      "ploidy_vs_viability_after_ms.pdf",
      "Ploidy vs Viability After MS",
      "Viability modifier after a one-copy-loss missegregation event across the ploidy grid."
    ),
    optional_figure(
      viz_dir,
      "death_rate_vs_missegregation_rate.pdf",
      "Death Rate vs Missegregation Rate",
      "Missegregation-rate curve plotted against the fitted death rate at the 2N and 4N reference ploidy states."
    )
  ))
  figure_index_range <- function(start, figs) {
    if (!length(figs)) return(integer(0))
    seq.int(start, length.out = length(figs))
  }
  oxygen_idx_start <- 1L
  oxygen_dynamics_idx <- figure_index_range(oxygen_idx_start, oxygen_dynamics_figs)
  oxygen_idx_start <- oxygen_idx_start + length(oxygen_dynamics_figs)
  oxygen_ms_relationship_idx <- figure_index_range(oxygen_idx_start, oxygen_ms_relationship_figs)
  oxygen_figs <- c(
    oxygen_dynamics_figs,
    oxygen_ms_relationship_figs
  )
  oxygen_figure_parts <- Filter(Negate(is.null), list(
    if (length(oxygen_ms_relationship_idx)) {
      list(
        part_index = 3L,
        title = "MS-Linked Viability and Death Relationships",
        description = "Missegregation-linked viability, nonviable daughter production, and death-rate relationships.",
        figure_indices = oxygen_ms_relationship_idx,
        cols = 3L
      )
    }
  ))

  sections <- list(
    list(
      name = "Burden",
      figures = burden_figs,
      layout_groups = list(list(indices = 1:2, cols = 2L))
    ),
    list(
      name = "Ploidy",
      figures = ploidy_figs,
      layout_groups = list(list(indices = 1:2, cols = 2L))
    ),
    list(
      name = "Oxygen / O2",
      figures = oxygen_figs,
      direct_figure_indices = oxygen_dynamics_idx,
      direct_layout_groups = list(list(indices = seq_along(oxygen_dynamics_idx), cols = min(2L, length(oxygen_dynamics_idx)))),
      figure_parts = oxygen_figure_parts
    )
  )
  Filter(function(section) length(section$figures) > 0L, sections)
}

flatten_section_figures <- function(sections) {
  unlist(lapply(sections, function(section) section$figures), recursive = FALSE)
}

optional_figure_with_layout <- function(viz_dir, filename, title, legend, layout_group = NULL, figure_index = NULL) {
  fig <- optional_figure(viz_dir, filename, title, legend)
  if (length(fig) > 0L) {
    fig[[1]]$layout_group <- layout_group
    if (!is.null(figure_index)) fig[[1]]$figure_index <- figure_index
  }
  fig
}

sequential_layout_groups_from_figures <- function(figures) {
  if (!length(figures)) return(list())
  groups <- list()
  i <- 1L
  while (i <= length(figures)) {
    group <- figures[[i]]$layout_group %||% ""
    if (nzchar(group)) {
      j <- i
      while (j < length(figures) && identical(figures[[j + 1L]]$layout_group %||% "", group)) {
        j <- j + 1L
      }
      run_indices <- seq.int(i, j)
      if (length(run_indices) > 1L) {
        if (length(run_indices) == 3L) {
          run_indices <- c(run_indices, NA_integer_)
        }
        groups <- c(groups, list(list(indices = run_indices, cols = 2L)))
        i <- j + 1L
      } else {
        groups <- c(groups, list(list(indices = i, cols = 1L)))
        i <- i + 1L
      }
    } else {
      groups <- c(groups, list(list(indices = i, cols = 1L)))
      i <- i + 1L
    }
  }
  groups
}

build_invitro_report_section_specs_for_joint <- function(fit_dir) {
  viz_dir <- file.path(fit_dir, "viz", "invitro")
  if (!dir.exists(viz_dir)) {
    return(list())
  }

  sections <- list(
    list(
      name = "Identifiability Diagnostics",
      figures = c(
        optional_figure(
          viz_dir,
          "invitro_identifiability_diagnostics.pdf",
          "Identifiability Diagnostics",
          "Local identifiability output for the current seed when available. If no Jacobian/Hessian was saved, the figure reports an optimizer-population proxy instead."
        )
      )
    ),
    list(
      name = "Rate-Function Diagnostics",
      figures = c(
        optional_figure_with_layout(
          viz_dir,
          "invitro_o2_selected_live_panels.pdf",
          "Constant External Oxygen and Selected-Day Live Cells",
          "Branch-aware diagnostic for the in vitro runner. Each cohort is split into control/deprived lineage panels using the same branch-specific x-axis as the aligned growth/chromosome/burden composite; repeated lineage passages are not averaged across branches. The upper row shows external oxygen and the lower row shows selected-day predicted live cells."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_rate_function_diagnostics.pdf",
          "Rate-Function Diagnostics",
          "Best-fit oxygen, ploidy, death, proliferation, and missegregation rate functions for the current seed."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_daily_counts.pdf",
          "Daily Live-Cell Trajectories",
          "Predicted live-cell trajectories split into 2N/control, 2N/deprived, 4N/control, and 4N/deprived panels, with each lineage passage shown as an inset subplot; selected propagation days are marked."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_growth_ploidy_burden_composite.pdf",
          "Aligned Growth, Chromosome Count, and Burden Fit",
          "Composite in vitro fit view. The 2N and 4N cohort blocks are stacked vertically; each block contains growth rate, chromosome-count quantile, and burden-decomposition rows. Repeated lineage passages are split into branch-specific O2 x-axis labels rather than averaged together; growth-rate and chromosome-count lines follow parent-child lineage links, so parallel p10 branches both connect to their shared p9 parent. Observed growth-rate points are drawn at their own branch positions. In vitro burden components are live/dead fractions normalized by the displayed predicted cell components, so the burden row ranges from 0 to 1. Control and deprived panels use their own passage ranges, with rows aligned within each lineage panel."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_flow_density.pdf",
          "Flow-Density Fit",
          "Observed G0/G1 ploidy-density curves are overlaid with the fitted flow-density prediction.",
          layout_group = "density-distribution"
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_distribution_heatmap.pdf",
          "Predicted Ploidy Distribution",
          "Full predicted chromosome-count distribution across in vitro passages.",
          layout_group = "density-distribution"
        )
      )
    ),
    list(
      name = "MS-Linked Viability and Death Relationships",
      figures = c(
        optional_figure_with_layout(
          viz_dir,
          "ms_rate_vs_nonviable_daughter_fraction.pdf",
          "Nonviable Daughter Fraction vs MS Rate",
          "Per-division fraction of daughter cells that are nonviable because of missegregation-linked loss, shown against missegregation rate.",
          figure_index = 3L
        ),
        optional_figure_with_layout(
          viz_dir,
          "death_rate_vs_missegregation_rate.pdf",
          "Death Rate vs Missegregation Rate",
          "Missegregation-rate curve plotted against the fitted death rate at the 2N and 4N reference ploidy states.",
          figure_index = 2L
        ),
        optional_figure_with_layout(
          viz_dir,
          "ploidy_vs_viability_after_ms.pdf",
          "Ploidy vs Viability After MS",
          "Viability modifier after a one-copy-loss missegregation event across the ploidy grid.",
          figure_index = 5L
        )
      ),
      layout_groups = list(
        list(indices = 1:3, cols = 3L)
      )
    )
  )
  sections <- lapply(sections, function(section) {
    section$figures <- Filter(Negate(is.null), section$figures)
    section$layout_groups <- section$layout_groups %||% sequential_layout_groups_from_figures(section$figures)
    section
  })
  Filter(function(section) length(section$figures) > 0L, sections)
}

joint_source_section_to_part <- function(source_section, part_index, global_start) {
  n_figs <- length(source_section$figures)
  if (!n_figs) return(NULL)
  part <- list(
    part_index = part_index,
    title = source_section$name %||% paste0("Part ", part_index),
    description = source_section$description %||% "",
    figure_indices = seq.int(global_start, length.out = n_figs)
  )

  source_parts <- source_section$figure_parts %||% list()
  direct_indices <- source_section$direct_figure_indices %||% integer(0)
  direct_indices <- direct_indices[direct_indices <= n_figs]
  if (length(source_parts) > 0L || length(direct_indices) > 0L) {
    subparts <- list()
    if (length(direct_indices) > 0L) {
      direct_layout_groups <- source_section$direct_layout_groups %||% list(
        list(indices = seq_along(direct_indices), cols = source_section$direct_cols %||% 1L)
      )
      subparts <- c(subparts, list(list(
        subpart_index = length(subparts) + 1L,
        title = if (identical(source_section$name %||% "", "Oxygen / O2")) "O2 Dynamics" else "Overview",
        description = source_section$direct_description %||% "",
        figure_indices = direct_indices,
        layout_groups = direct_layout_groups
      )))
    }
    for (source_part in source_parts) {
      source_part$subpart_index <- length(subparts) + 1L
      source_part$part_index <- NULL
      subparts <- c(subparts, list(source_part))
    }
    part$subparts <- subparts
    return(part)
  }

  layout_groups <- source_section$layout_groups %||% NULL
  if (length(layout_groups) > 0L) {
    part$layout_groups <- layout_groups
  } else {
    part$cols <- source_section$cols %||% 1L
  }
  part
}

build_joint_scope_section <- function(name, source_sections) {
  source_sections <- Filter(function(section) length(section$figures) > 0L, source_sections)
  if (!length(source_sections)) return(NULL)
  figures <- flatten_section_figures(source_sections)
  parts <- list()
  offset <- 0L
  for (i in seq_along(source_sections)) {
    part <- joint_source_section_to_part(source_sections[[i]], i, offset + 1L)
    if (!is.null(part)) {
      parts <- c(parts, list(part))
    }
    offset <- offset + length(source_sections[[i]]$figures)
  }
  list(name = name, figures = figures, figure_parts = parts)
}

read_comparison_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

comparison_finite_rows <- function(df, cols) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(data.frame())
  ok <- rep(TRUE, nrow(df))
  for (col in cols) {
    if (!col %in% names(df)) return(data.frame())
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    ok <- ok & is.finite(df[[col]])
  }
  df[ok, , drop = FALSE]
}

comparison_context_bind <- function(invivo_df, invitro_df, invivo_label, invitro_label) {
  if (is.null(invivo_df) || is.null(invitro_df) || !nrow(invivo_df) || !nrow(invitro_df)) {
    return(data.frame())
  }
  invivo_df$context <- invivo_label
  invitro_df$context <- invitro_label
  common_cols <- intersect(names(invivo_df), names(invitro_df))
  out <- rbind(invivo_df[, common_cols, drop = FALSE], invitro_df[, common_cols, drop = FALSE])
  out$context <- factor(out$context, levels = c(invivo_label, invitro_label))
  out
}

comparison_palette <- function(levels) {
  levels <- as.character(levels)
  known <- c("2N" = "#1f77b4", "4N" = "#d62728")
  if (all(levels %in% names(known))) return(known[levels])
  stats::setNames(grDevices::hcl.colors(length(levels), palette = "Dark 3"), levels)
}

comparison_state_axis_label <- function(df) {
  endpoint <- suppressWarnings(as.numeric(df$endpoint_value))
  if (length(endpoint) && any(is.finite(endpoint) & endpoint > 10)) {
    "Chromosome Number (N)"
  } else {
    "Ploidy"
  }
}

save_comparison_plot_pair <- function(plot, out_dir, basename, width = 14, height = 7) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(basename, ".pdf"))
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 180, bg = "white")
  TRUE
}

generate_invivo_invitro_comparison_figures <- function(fit_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
  invivo_dir <- file.path(fit_dir, "viz", "invivo")
  invitro_dir <- file.path(fit_dir, "viz", "invitro")
  out_dir <- file.path(fit_dir, "viz", "invivo_vs_invitro")
  if (!dir.exists(invivo_dir) || !dir.exists(invitro_dir)) return(FALSE)

  read_pair <- function(invivo_file, invitro_file, invivo_label, invitro_label) {
    comparison_context_bind(
      read_comparison_tsv(file.path(invivo_dir, invivo_file)),
      read_comparison_tsv(file.path(invitro_dir, invitro_file)),
      invivo_label,
      invitro_label
    )
  }

  context_basic <- c("In vivo", "In vitro")
  generated <- logical(0)

  oxygen_curve <- read_pair(
    "functional_curve_oxygen.tsv",
    "functional_curve_oxygen.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  oxygen_curve <- comparison_finite_rows(oxygen_curve, c("death_rate", "ms_rate"))
  if (nrow(oxygen_curve) > 0L && "cohort" %in% names(oxygen_curve)) {
    oxygen_curve$cohort <- factor(oxygen_curve$cohort, levels = unique(oxygen_curve$cohort))
    p <- ggplot2::ggplot(
      oxygen_curve,
      ggplot2::aes(x = death_rate, y = ms_rate, color = cohort, group = cohort)
    ) +
      ggplot2::geom_path(linewidth = 1) +
      ggplot2::facet_grid(. ~ context) +
      ggplot2::scale_color_manual(values = comparison_palette(levels(oxygen_curve$cohort)), drop = FALSE) +
      ggplot2::labs(
        title = "In Vivo vs In Vitro: Death Rate vs Missegregation Rate",
        subtitle = "Left: in vivo; right: in vitro.",
        x = "Death rate",
        y = "MS rate",
        color = "Cohort"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"))
    generated[["death_rate_vs_missegregation_rate"]] <- save_comparison_plot_pair(
      p,
      out_dir,
      "invivo_vs_invitro_death_rate_vs_missegregation_rate"
    )
  }

  multi_curve_basic <- read_pair(
    "functional_curve_oxygen_multi_ploidy.tsv",
    "functional_curve_oxygen_multi_ploidy.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  multi_curve_basic <- comparison_finite_rows(
    multi_curve_basic,
    c("ms_rate", "misseg_nonviable_daughter_fraction", "misseg_nonviable_division_prob")
  )
  if (nrow(multi_curve_basic) > 0L && "cohort" %in% names(multi_curve_basic)) {
    multi_curve_basic$cohort <- factor(multi_curve_basic$cohort, levels = unique(multi_curve_basic$cohort))
    colors <- comparison_palette(levels(multi_curve_basic$cohort))
    make_ms_plot <- function(value_col, title, y_label, basename) {
      plot_df <- multi_curve_basic
      plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
      plot_df <- comparison_finite_rows(plot_df, c("ms_rate", "value"))
      if (!nrow(plot_df)) return(FALSE)
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = ms_rate, y = value, color = cohort)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_grid(. ~ context) +
        ggplot2::scale_color_manual(values = colors, drop = FALSE) +
        ggplot2::labs(
          title = paste0("In Vivo vs In Vitro: ", title),
          subtitle = "Left: in vivo; right: in vitro.",
          x = "MS rate",
          y = y_label,
          color = "Reference state"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"))
      save_comparison_plot_pair(p, out_dir, basename)
    }
    generated[["ms_rate_vs_nonviable_daughter_fraction"]] <- make_ms_plot(
      "misseg_nonviable_daughter_fraction",
      "Nonviable Daughter Fraction vs MS Rate",
      "Nonviable daughters / all daughters",
      "invivo_vs_invitro_ms_rate_vs_nonviable_daughter_fraction"
    )
    generated[["ms_rate_vs_nonviable_division_probability"]] <- make_ms_plot(
      "misseg_nonviable_division_prob",
      "Capped Nonviable Daughter Burden vs MS Rate",
      "min(expected nonviable daughters / division, 1)",
      "invivo_vs_invitro_ms_rate_vs_nonviable_division_probability"
    )
  }

  viability_curve <- read_pair(
    "functional_curve_ploidy.tsv",
    "functional_curve_ploidy.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  viability_curve <- comparison_finite_rows(viability_curve, c("endpoint_value", "viability_after_ms"))
  if (nrow(viability_curve) > 0L) {
    p <- ggplot2::ggplot(
      viability_curve,
      ggplot2::aes(x = endpoint_value, y = viability_after_ms)
    ) +
      ggplot2::geom_line(color = "#2ca02c", linewidth = 1) +
      ggplot2::facet_grid(. ~ context) +
      ggplot2::labs(
        title = "In Vivo vs In Vitro: Ploidy vs Viability After MS",
        subtitle = "Left: in vivo; right: in vitro.",
        x = comparison_state_axis_label(viability_curve),
        y = "Viability after MS"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"))
    generated[["ploidy_vs_viability_after_ms"]] <- save_comparison_plot_pair(
      p,
      out_dir,
      "invivo_vs_invitro_ploidy_vs_viability_after_ms"
    )
  }

  ploidy_resource <- read_pair(
    "functional_curve_ploidy_by_o2.tsv",
    "functional_curve_ploidy_by_o2.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  ploidy_resource <- comparison_finite_rows(ploidy_resource, c("endpoint_value", "oxygen_pct", "death_rate", "proliferation_rate"))
  if (nrow(ploidy_resource) > 0L) {
    make_ploidy_o2_plot <- function(value_col, title, y_label, basename) {
      plot_df <- ploidy_resource
      plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
      plot_df <- comparison_finite_rows(plot_df, c("endpoint_value", "oxygen_pct", "value"))
      if (!nrow(plot_df)) return(FALSE)
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = endpoint_value, y = value, color = oxygen_pct)
      ) +
        ggplot2::geom_point(shape = 15, size = 1.8, alpha = 0.95) +
        ggplot2::facet_grid(. ~ context) +
        ggplot2::scale_color_gradient(low = "#2C7BB6", high = "#F28E2B", name = "O2 level") +
        ggplot2::labs(
          title = paste0("In Vivo vs In Vitro: ", title),
          subtitle = "Left: in vivo; right: in vitro.",
          x = comparison_state_axis_label(plot_df),
          y = y_label
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"))
      save_comparison_plot_pair(p, out_dir, basename)
    }
    generated[["ploidy_vs_death_rate_by_o2"]] <- make_ploidy_o2_plot(
      "death_rate",
      "Ploidy vs Death Rate by O2",
      "Death rate",
      "invivo_vs_invitro_ploidy_vs_death_rate_by_o2"
    )
    generated[["ploidy_vs_proliferation_rate_by_o2"]] <- make_ploidy_o2_plot(
      "proliferation_rate",
      "Ploidy vs Proliferation Rate by O2",
      "Proliferation rate",
      "invivo_vs_invitro_ploidy_vs_proliferation_rate_by_o2"
    )
  }

  multi_curve_resource <- read_pair(
    "functional_curve_oxygen_multi_ploidy.tsv",
    "functional_curve_oxygen_multi_ploidy.tsv",
    context_basic[[1]],
    context_basic[[2]]
  )
  multi_curve_resource <- comparison_finite_rows(
    multi_curve_resource,
    c("oxygen_pct", "ms_rate", "proliferation_rate", "death_rate")
  )
  if (nrow(multi_curve_resource) > 0L && "cohort" %in% names(multi_curve_resource)) {
    multi_curve_resource$cohort <- factor(multi_curve_resource$cohort, levels = unique(multi_curve_resource$cohort))
    colors <- comparison_palette(levels(multi_curve_resource$cohort))
    make_o2_plot <- function(value_col, title, y_label, basename) {
      plot_df <- multi_curve_resource
      plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
      plot_df <- comparison_finite_rows(plot_df, c("oxygen_pct", "value"))
      if (!nrow(plot_df)) return(FALSE)
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = oxygen_pct, y = value, color = cohort)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_grid(. ~ context) +
        ggplot2::scale_color_manual(values = colors, drop = FALSE) +
        ggplot2::labs(
          title = paste0("In Vivo vs In Vitro: ", title),
          subtitle = "Left: in vivo; right: in vitro.",
          x = "Oxygen (%)",
          y = y_label,
          color = "Reference state"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"))
      save_comparison_plot_pair(p, out_dir, basename)
    }
    generated[["oxygen_vs_missegregation_rate_multi_ploidy"]] <- make_o2_plot(
      "ms_rate",
      "Oxygen vs Missegregation Rate",
      "MS rate",
      "invivo_vs_invitro_oxygen_vs_missegregation_rate_multi_ploidy"
    )
    generated[["oxygen_vs_proliferation_rate"]] <- make_o2_plot(
      "proliferation_rate",
      "Oxygen vs Proliferation Rate",
      "Proliferation rate",
      "invivo_vs_invitro_oxygen_vs_proliferation_rate"
    )
    generated[["oxygen_vs_death_rate"]] <- make_o2_plot(
      "death_rate",
      "Oxygen vs Death Rate",
      "Death rate",
      "invivo_vs_invitro_oxygen_vs_death_rate"
    )
  }

  any(generated, na.rm = TRUE)
}

build_invivo_invitro_comparison_section <- function(fit_dir) {
  generate_invivo_invitro_comparison_figures(fit_dir)
  viz_dir <- file.path(fit_dir, "viz", "invivo_vs_invitro")
  if (!dir.exists(viz_dir)) return(NULL)
  figures <- c(
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ms_rate_vs_nonviable_daughter_fraction.pdf",
      "Nonviable Daughter Fraction vs MS Rate",
      "Left panel uses the in vivo functional-response output; right panel uses the in vitro functional-response output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ploidy_vs_viability_after_ms.pdf",
      "Ploidy vs Viability After MS",
      "Left panel uses the in vivo viability curve; right panel uses the in vitro viability curve."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_death_rate_vs_missegregation_rate.pdf",
      "Death Rate vs Missegregation Rate",
      "Left panel uses the in vivo functional-response output; right panel uses the in vitro functional-response output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ploidy_vs_death_rate_by_o2.pdf",
      "Ploidy vs Death Rate by O2",
      "Left panel uses the in vivo O2-only resource output; right panel uses the in vitro O2-only resource output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ploidy_vs_proliferation_rate_by_o2.pdf",
      "Ploidy vs Proliferation Rate by O2",
      "Left panel uses the in vivo O2-only resource output; right panel uses the in vitro O2-only resource output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_oxygen_vs_missegregation_rate_multi_ploidy.pdf",
      "Oxygen vs Missegregation Rate",
      "Left panel uses the in vivo O2-only resource output; right panel uses the in vitro O2-only resource output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_oxygen_vs_proliferation_rate.pdf",
      "Oxygen vs Proliferation Rate",
      "Left panel uses the in vivo O2-only resource output; right panel uses the in vitro O2-only resource output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_oxygen_vs_death_rate.pdf",
      "Oxygen vs Death Rate",
      "Left panel uses the in vivo O2-only resource output; right panel uses the in vitro O2-only resource output."
    )
  )
  figures <- Filter(Negate(is.null), figures)
  if (!length(figures)) return(NULL)
  layout_groups <- lapply(seq(1L, length(figures), by = 2L), function(i) {
    idx <- if (i + 1L <= length(figures)) c(i, i + 1L) else c(i, NA_integer_)
    list(indices = idx, cols = 2L)
  })
  list(
    name = "In Vivo vs In Vitro",
    figures = figures,
    layout_groups = layout_groups
  )
}

build_section_specs <- function(fit_dir) {
  fit_summary <- read_fit_summary_map(fit_dir)
  invivo_sections <- build_invivo_section_specs(fit_dir)
  if (!identical(as.character(fit_summary[["fit_mode"]]), "fit_joint")) {
    return(invivo_sections)
  }
  Filter(Negate(is.null), list(
    build_joint_scope_section("In Vivo", invivo_sections),
    build_joint_scope_section("In Vitro", build_invitro_report_section_specs_for_joint(fit_dir)),
    build_invivo_invitro_comparison_section(fit_dir)
  ))
}

stage_assets <- function(section_specs) {
  assets_dir <- file.path(
    tempdir(),
    paste0("o2_report_assets_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  )
  dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
  use_magick <- report_magick_available()
  use_gs <- !use_magick && report_gs_available()
  for (i in seq_along(section_specs)) {
    if (length(section_specs[[i]]$figures) == 0) next
    for (j in seq_along(section_specs[[i]]$figures)) {
      src <- section_specs[[i]]$figures[[j]]$src
      stage_prefix <- sprintf("section%02d_figure%03d_", i, j)
      pdf_stage <- file.path(assets_dir, paste0(stage_prefix, basename(src)))
      if (!file.copy(src, pdf_stage, overwrite = TRUE)) {
        stop("Failed to stage PDF asset: ", src)
      }
      section_specs[[i]]$figures[[j]]$pdf_asset_abs <- normalizePath(pdf_stage, mustWork = TRUE)
      html_src <- section_specs[[i]]$figures[[j]]$html_src %||% NULL
      if (!is.null(html_src) && file.exists(html_src)) {
        png_stage <- file.path(assets_dir, paste0(stage_prefix, basename(html_src)))
        if (!file.copy(html_src, png_stage, overwrite = TRUE)) {
          stop("Failed to stage PNG asset: ", html_src)
        }
        section_specs[[i]]$figures[[j]]$html_embed_kind <- "img"
        section_specs[[i]]$figures[[j]]$html_asset_abs <- normalizePath(png_stage, mustWork = TRUE)
        section_specs[[i]]$figures[[j]]$html_asset_uri <- normalizePath(png_stage, mustWork = TRUE)
      } else if (use_magick || use_gs) {
        png_stage <- sub("\\.pdf$", ".png", pdf_stage, ignore.case = TRUE)
        png_stage <- if (use_magick) {
          render_pdf_preview_png(pdf_stage, png_stage)
        } else {
          render_pdf_preview_png_gs(pdf_stage, png_stage)
        }
        section_specs[[i]]$figures[[j]]$html_embed_kind <- "img"
        section_specs[[i]]$figures[[j]]$html_asset_abs <- normalizePath(png_stage, mustWork = TRUE)
        section_specs[[i]]$figures[[j]]$html_asset_uri <- normalizePath(png_stage, mustWork = TRUE)
      } else {
        section_specs[[i]]$figures[[j]]$html_embed_kind <- "pdf_object"
        section_specs[[i]]$figures[[j]]$html_asset_uri <- pdf_to_data_uri(pdf_stage)
      }
    }
  }
  section_specs
}

fit_report_fig_title <- function(fig) {
  title <- fig$title %||% ""
  if (!nzchar(trimws(title)) && !is.null(fig$pdf_asset_abs) && length(fig$pdf_asset_abs)) {
    title <- tools::file_path_sans_ext(basename(fig$pdf_asset_abs[[1]]))
  }
  if (!nzchar(trimws(title))) title <- "Untitled Figure"
  title
}

fit_report_fig_legend <- function(fig) {
  legend <- fig$legend %||% ""
  if (!nzchar(trimws(legend)) && !is.null(fig$pdf_asset_abs) && length(fig$pdf_asset_abs)) {
    legend <- paste0("Figure source: ", basename(fig$pdf_asset_abs[[1]]), ".")
  }
  if (!nzchar(trimws(legend))) legend <- "Figure legend unavailable."
  legend
}

fit_report_figure_data_uri <- function(fig) {
  embed_kind <- fig$html_embed_kind %||% "img"
  if (identical(embed_kind, "img")) {
    img_path <- fig$html_asset_abs %||% fig$html_asset_uri
    if (is.character(img_path) && length(img_path) && startsWith(img_path[[1]], "data:")) {
      return(img_path[[1]])
    }
    return(file_to_data_uri(img_path[[1]], mime = "image/png"))
  }
  pdf_uri <- fig$html_asset_uri %||% ""
  if (nzchar(pdf_uri) && startsWith(pdf_uri, "data:")) return(pdf_uri)
  pdf_to_data_uri(fig$pdf_asset_abs[[1]])
}

fit_report_figure_card_id <- function(section_index, figure_index) {
  sprintf("figure-%02d-%03d", section_index, figure_index)
}

fit_report_figure_label <- function(
    fig,
    section_index,
    figure_index,
    part_index = NULL,
    part_figure_index = NULL,
    subpart_index = NULL,
    subpart_figure_index = NULL) {
  fig_title <- fit_report_fig_title(fig)
  figure_index_use <- suppressWarnings(as.integer(fig$figure_index %||% figure_index))
  if (!is.finite(figure_index_use) || figure_index_use < 1L) figure_index_use <- figure_index
  part_index_use <- suppressWarnings(as.integer(part_index %||% NA_integer_))
  part_figure_index_use <- suppressWarnings(as.integer(fig$figure_index %||% part_figure_index %||% NA_integer_))
  subpart_index_use <- suppressWarnings(as.integer(subpart_index %||% NA_integer_))
  subpart_figure_index_use <- suppressWarnings(as.integer(fig$figure_index %||% subpart_figure_index %||% NA_integer_))
  if (is.finite(part_index_use) && part_index_use > 0L &&
      is.finite(subpart_index_use) && subpart_index_use > 0L &&
      is.finite(subpart_figure_index_use) && subpart_figure_index_use > 0L) {
    return(sprintf(
      "Figure 3.%d.%d.%d.%d %s",
      section_index,
      part_index_use,
      subpart_index_use,
      subpart_figure_index_use,
      fig_title
    ))
  }
  if (is.finite(part_index_use) && part_index_use > 0L &&
      is.finite(part_figure_index_use) && part_figure_index_use > 0L) {
    return(sprintf("Figure 3.%d.%d.%d %s", section_index, part_index_use, part_figure_index_use, fig_title))
  }
  sprintf("Figure %d.%d %s", section_index, figure_index_use, fig_title)
}

render_report_figure_card <- function(
    fig,
    section_index,
    figure_index,
    part_index = NULL,
    part_figure_index = NULL,
    subpart_index = NULL,
    subpart_figure_index = NULL) {
  fig_legend <- fit_report_fig_legend(fig)
  fig_label <- fit_report_figure_label(
    fig,
    section_index,
    figure_index,
    part_index = part_index,
    part_figure_index = part_figure_index,
    subpart_index = subpart_index,
    subpart_figure_index = subpart_figure_index
  )
  fig_id <- fit_report_figure_card_id(section_index, figure_index)
  media <- if (identical(fig$html_embed_kind %||% "img", "pdf_object")) {
    sprintf(
      paste0(
        '<div class="report-figure">',
        '<object data="%s" type="application/pdf" class="report-figure-object">',
        '<div class="report-figure-fallback"><a href="%s">Open PDF figure</a></div>',
        '</object></div>'
      ),
      fit_report_figure_data_uri(fig),
      escape_html(basename(fig$pdf_asset_abs[[1]]))
    )
  } else {
    sprintf(
      '<div class="report-figure"><img src="%s" alt="%s" class="report-figure-image"/></div>',
      fit_report_figure_data_uri(fig),
      escape_html(fig_label)
    )
  }
  paste0(
    '<article class="report-figure-card" id="', fig_id, '">',
    media,
    '<p class="report-figure-title"><strong>', escape_html(fig_label), "</strong></p>",
    '<p class="report-figure-legend">', escape_html(fig_legend), "</p>",
    "</article>"
  )
}

render_report_blank_figure_card <- function() {
  paste0(
    '<article class="report-figure-card report-figure-card--blank">',
    '<div class="report-figure report-figure-blank"></div>',
    "</article>"
  )
}

figure_layout_groups_from_specs <- function(n_figs, layout_groups = NULL, default_cols = 1L) {
  if (length(layout_groups) > 0L) {
    groups <- lapply(layout_groups, function(group) {
      idx <- group$indices[is.na(group$indices) | group$indices <= n_figs]
      if (!length(idx)) return(NULL)
      if (!any(!is.na(idx))) return(NULL)
      list(indices = idx, cols = group$cols %||% default_cols)
    })
    groups <- Filter(Negate(is.null), groups)
    used <- unlist(lapply(groups, function(group) group$indices[!is.na(group$indices)]), use.names = FALSE)
    extra <- setdiff(seq_len(n_figs), used)
    if (length(extra) > 0L) {
      groups <- c(groups, lapply(extra, function(i) list(indices = i, cols = default_cols)))
    }
    return(groups)
  }
  list(list(indices = seq_len(n_figs), cols = default_cols))
}

in_vivo_figure_layout_groups <- function(n_figs) {
  final_row <- c(12L, 13L, if (n_figs >= 15L) seq.int(15L, n_figs) else integer(0))
  requested <- list(
    list(indices = 1:2, cols = 2L),
    list(indices = 3:5, cols = 3L),
    list(indices = 6:7, cols = 2L),
    list(indices = 8:10, cols = 3L),
    list(indices = c(11L, 14L), cols = 2L),
    list(indices = final_row, cols = 3L)
  )
  groups <- lapply(requested, function(group) {
    idx <- group$indices[is.na(group$indices) | group$indices <= n_figs]
    if (!length(idx)) return(NULL)
    list(indices = idx, cols = group$cols)
  })
  groups <- Filter(Negate(is.null), groups)
  used <- unlist(lapply(groups, function(group) group$indices[!is.na(group$indices)]), use.names = FALSE)
  extra <- setdiff(seq_len(n_figs), used)
  if (length(extra) > 0L) {
    groups <- c(groups, lapply(extra, function(i) list(indices = i, cols = 1L)))
  }
  groups
}

in_vitro_figure_layout_groups <- function(n_figs) {
  requested <- list(
    list(indices = 1:2, cols = 2L),
    list(indices = 3:4, cols = 2L),
    list(indices = c(5L, 6L, 7L, NA_integer_), cols = 2L)
  )
  groups <- lapply(requested, function(group) {
    available <- !is.na(group$indices) & group$indices <= n_figs
    if (!any(available)) return(NULL)
    idx <- group$indices[is.na(group$indices) | available]
    list(indices = idx, cols = group$cols)
  })
  groups <- Filter(Negate(is.null), groups)
  used <- unlist(lapply(groups, function(group) group$indices[!is.na(group$indices)]), use.names = FALSE)
  extra <- setdiff(seq_len(n_figs), used)
  if (length(extra) > 0L) {
    groups <- c(groups, lapply(extra, function(i) list(indices = i, cols = 1L)))
  }
  groups
}

section_figure_layout_groups <- function(section, n_figs) {
  layout_groups <- section$layout_groups %||% NULL
  if (length(layout_groups) > 0L) {
    return(figure_layout_groups_from_specs(n_figs, layout_groups = layout_groups, default_cols = 1L))
  }
  if (identical(section$name %||% "", "In Vivo")) {
    return(in_vivo_figure_layout_groups(n_figs))
  }
  if (identical(section$name %||% "", "In Vitro")) {
    return(in_vitro_figure_layout_groups(n_figs))
  }
  lapply(seq_len(n_figs), function(i) list(indices = i, cols = 1L))
}

part_figure_layout_groups <- function(part, n_figs) {
  layout_groups <- part$layout_groups %||% NULL
  if (length(layout_groups) > 0L) {
    return(figure_layout_groups_from_specs(n_figs, layout_groups = layout_groups, default_cols = 1L))
  }
  cols <- suppressWarnings(as.integer(part$cols %||% 1L))
  if (!is.finite(cols) || cols < 1L) cols <- 1L
  list(list(indices = seq_len(n_figs), cols = cols))
}

render_section_figure_parts <- function(section, section_index) {
  parts <- section$figure_parts %||% list()
  if (length(parts) == 0L) return(NULL)
  paste(vapply(parts, function(part) {
    fig_indices <- part$figure_indices
    fig_indices <- fig_indices[fig_indices <= length(section$figures)]
    if (!length(fig_indices)) return("")
    part_figs <- section$figures[fig_indices]
    part_index <- suppressWarnings(as.integer(part$part_index %||% NA_integer_))
    if (!is.finite(part_index) || part_index < 1L) part_index <- NA_integer_
    part_title <- part$title %||% paste0("Part ", part_index)
    part_description <- part$description %||% ""
    subparts <- part$subparts %||% list()
    if (length(subparts) > 0L) {
      body <- paste(vapply(subparts, function(subpart) {
        subpart_indices <- subpart$figure_indices
        subpart_indices <- subpart_indices[subpart_indices <= length(part_figs)]
        if (!length(subpart_indices)) return("")
        subpart_figs <- part_figs[subpart_indices]
        subpart_global_indices <- fig_indices[subpart_indices]
        subpart_index <- suppressWarnings(as.integer(subpart$subpart_index %||% NA_integer_))
        if (!is.finite(subpart_index) || subpart_index < 1L) subpart_index <- 1L
        subpart_title <- subpart$title %||% paste0("Subpart ", subpart_index)
        subpart_description <- subpart$description %||% ""
        groups <- part_figure_layout_groups(subpart, length(subpart_figs))
        subpart_body <- paste(vapply(groups, function(group) {
          cols <- suppressWarnings(as.integer(group$cols[[1]]))
          if (!is.finite(cols) || cols < 1L) cols <- 1L
          cards <- paste(vapply(group$indices, function(j) {
            if (is.na(j)) return(render_report_blank_figure_card())
            render_report_figure_card(
              subpart_figs[[j]],
              section_index,
              subpart_global_indices[[j]],
              part_index = part_index,
              subpart_index = subpart_index,
              subpart_figure_index = j
            )
          }, character(1)), collapse = "")
          paste0(
            '<div class="report-figure-grid report-figure-grid--', cols, '">',
            cards,
            "</div>"
          )
        }, character(1)), collapse = "")
        paste0(
          '<section class="report-figure-subpart" id="', figure_subpart_id(section_index, part, subpart), '">',
          '<h5>', sprintf("3.%d.%d.%d ", section_index, part_index, subpart_index), escape_html(subpart_title), "</h5>",
          if (nzchar(trimws(subpart_description))) {
            paste0('<p class="report-part-description">', escape_html(subpart_description), "</p>")
          } else {
            ""
          },
          subpart_body,
          "</section>"
        )
      }, character(1)), collapse = "")
    } else {
      groups <- part_figure_layout_groups(part, length(part_figs))
      body <- paste(vapply(groups, function(group) {
        cols <- suppressWarnings(as.integer(group$cols[[1]]))
        if (!is.finite(cols) || cols < 1L) cols <- 1L
        cards <- paste(vapply(group$indices, function(j) {
          if (is.na(j)) return(render_report_blank_figure_card())
          render_report_figure_card(
            part_figs[[j]],
            section_index,
            fig_indices[[j]],
            part_index = part_index,
            part_figure_index = j
          )
        }, character(1)), collapse = "")
        paste0(
          '<div class="report-figure-grid report-figure-grid--', cols, '">',
          cards,
          "</div>"
        )
      }, character(1)), collapse = "")
    }
    paste0(
      '<section class="report-figure-part" id="', figure_part_id(section_index, part), '">',
      '<h4>', sprintf("3.%d.%d ", section_index, part_index), escape_html(part_title), "</h4>",
      if (nzchar(trimws(part_description))) {
        paste0('<p class="report-part-description">', escape_html(part_description), "</p>")
      } else {
        ""
      },
      body,
      "</section>"
    )
  }, character(1)), collapse = "")
}

render_section_direct_figure_blocks <- function(section, section_index) {
  direct_indices <- section$direct_figure_indices %||% integer(0)
  direct_indices <- direct_indices[direct_indices <= length(section$figures)]
  if (!length(direct_indices)) return("")
  direct_layout_groups <- section$direct_layout_groups %||% list(
    list(
      indices = seq_along(direct_indices),
      cols = section$direct_cols %||% 1L
    )
  )
  groups <- figure_layout_groups_from_specs(
    length(direct_indices),
    layout_groups = direct_layout_groups,
    default_cols = 1L
  )
  paste(vapply(groups, function(group) {
    cols <- suppressWarnings(as.integer(group$cols[[1]]))
    if (!is.finite(cols) || cols < 1L) cols <- 1L
    cards <- paste(vapply(group$indices, function(j) {
      global_idx <- direct_indices[[j]]
      render_report_figure_card(section$figures[[global_idx]], section_index, global_idx)
    }, character(1)), collapse = "")
    paste0(
      '<div class="report-figure-grid report-figure-grid--', cols, '">',
      cards,
      "</div>"
    )
  }, character(1)), collapse = "")
}

render_section_figure_blocks <- function(section, section_index) {
  n_figs <- length(section$figures)
  if (n_figs == 0L) return("")
  part_blocks <- render_section_figure_parts(section, section_index)
  direct_blocks <- render_section_direct_figure_blocks(section, section_index)
  if (nzchar(direct_blocks) && !is.null(part_blocks) && nzchar(part_blocks)) {
    return(paste0(direct_blocks, part_blocks))
  }
  if (nzchar(direct_blocks)) return(direct_blocks)
  if (!is.null(part_blocks) && nzchar(part_blocks)) return(part_blocks)
  groups <- section_figure_layout_groups(section, n_figs)
  paste(vapply(groups, function(group) {
    cols <- suppressWarnings(as.integer(group$cols[[1]]))
    if (!is.finite(cols) || cols < 1L) cols <- 1L
    cards <- paste(vapply(group$indices, function(j) {
      if (is.na(j)) return(render_report_blank_figure_card())
      render_report_figure_card(section$figures[[j]], section_index, j)
    }, character(1)), collapse = "")
    paste0(
      '<div class="report-figure-grid report-figure-grid--', cols, '">',
      cards,
      "</div>"
    )
  }, character(1)), collapse = "")
}

fit_report_table_html <- function(tab, class_name = "report-table") {
  if (is.null(tab) || !NROW(tab)) {
    return('<p class="report-empty">No data available.</p>')
  }
  tab <- format_report_table(tab)
  headers <- paste(sprintf("<th>%s</th>", escape_html(names(tab))), collapse = "")
  rows <- vapply(seq_len(nrow(tab)), function(i) {
    vals <- vapply(tab[i, , drop = FALSE], function(x) escape_html(x[[1]]), character(1))
    paste0("<tr>", paste(sprintf("<td>%s</td>", vals), collapse = ""), "</tr>")
  }, character(1))
  paste0(
    '<table class="', class_name, '"><thead><tr>', headers, '</tr></thead><tbody>',
    paste(rows, collapse = ""),
    "</tbody></table>"
  )
}

build_parameter_blocks_html <- function(parameter_sections) {
  if (length(parameter_sections) == 0L) {
    return('<p class="report-empty">No parameter table available.</p>')
  }
  current_group <- NULL
  group_index <- 0L
  blocks <- character(0)
  for (i in seq_along(parameter_sections)) {
    section <- parameter_sections[[i]]
    group <- as.character(section$group %||% "")
    if (nzchar(group) && !identical(group, current_group)) {
      group_index <- group_index + 1L
      blocks <- c(
        blocks,
        paste0(
          '<h3 class="report-param-group" id="',
          parameter_group_id(group_index, group),
          '">',
          escape_html(group),
          "</h3>"
        )
      )
      current_group <- group
    }
    heading_tag <- if (nzchar(group)) "h4" else "h3"
    section_name <- section$name %||% paste0("Parameter Table ", i)
    blocks <- c(
      blocks,
      paste0(
        "<", heading_tag, ' class="report-param-section" id="',
        parameter_section_id(i, section_name),
        '">',
        escape_html(section_name),
        "</", heading_tag, ">",
        fit_report_table_html(section$table)
      )
    )
  }
  paste(blocks, collapse = "")
}

report_nav_item <- function(id, label, level_class) {
  sprintf(
    '<li class="report-nav-item"><a class="report-nav-link %s" href="#%s" data-target-id="%s">%s</a></li>',
    level_class,
    id,
    id,
    escape_html(label)
  )
}

append_figure_nav_items <- function(
    nav_items,
    section,
    section_index,
    figure_indices,
    level_class,
    part_index = NULL,
    subpart_index = NULL) {
  figure_indices <- figure_indices[figure_indices <= length(section$figures)]
  if (!length(figure_indices)) return(nav_items)
  for (local_idx in seq_along(figure_indices)) {
    global_idx <- figure_indices[[local_idx]]
    if (is.na(global_idx) || global_idx < 1L || global_idx > length(section$figures)) next
    fig_label <- fit_report_figure_label(
      section$figures[[global_idx]],
      section_index,
      global_idx,
      part_index = part_index,
      part_figure_index = if (is.null(subpart_index)) local_idx else NULL,
      subpart_index = subpart_index,
      subpart_figure_index = if (!is.null(subpart_index)) local_idx else NULL
    )
    nav_items <- c(
      nav_items,
      report_nav_item(
        fit_report_figure_card_id(section_index, global_idx),
        fig_label,
        level_class
      )
    )
  }
  nav_items
}

build_report_nav_items <- function(sections, parameter_sections) {
  nav_items <- c(
    report_nav_item("report-metadata", "Report Metadata", "report-nav-h2"),
    report_nav_item("fit-summary", "1. Fit Summary", "report-nav-h2"),
    report_nav_item("best-parameters", "2. Parameters", "report-nav-h2")
  )
  if (length(parameter_sections) > 0L) {
    current_group <- NULL
    group_index <- 0L
    for (i in seq_along(parameter_sections)) {
      section <- parameter_sections[[i]]
      group <- as.character(section$group %||% "")
      if (nzchar(group) && !identical(group, current_group)) {
        group_index <- group_index + 1L
        group_id <- parameter_group_id(group_index, group)
        nav_items <- c(nav_items, report_nav_item(group_id, group, "report-nav-h3"))
        current_group <- group
      }
      section_name <- section$name %||% paste0("Parameter Table ", i)
      nav_class <- if (nzchar(group)) "report-nav-h4" else "report-nav-h3"
      nav_items <- c(nav_items, report_nav_item(parameter_section_id(i, section_name), section_name, nav_class))
    }
  }
  nav_items <- c(nav_items, report_nav_item("figures", "3. Figures", "report-nav-h2"))
  if (length(sections) > 0L) {
    for (i in seq_along(sections)) {
      slug <- paste0("section-", i)
      nav_items <- c(nav_items, report_nav_item(slug, sprintf("3.%d %s", i, sections[[i]]$name), "report-nav-h3"))
      parts <- sections[[i]]$figure_parts %||% list()
      direct_indices <- sections[[i]]$direct_figure_indices %||% integer(0)
      direct_indices <- direct_indices[direct_indices <= length(sections[[i]]$figures)]
      if (length(direct_indices) > 0L) {
        nav_items <- append_figure_nav_items(
          nav_items,
          sections[[i]],
          i,
          direct_indices,
          "report-nav-h4"
        )
      }
      if (length(parts) > 0L) {
        for (part in parts) {
          part_index <- suppressWarnings(as.integer(part$part_index %||% NA_integer_))
          if (!is.finite(part_index) || part_index < 1L) part_index <- i
          part_title <- part$title %||% paste0("Part ", part_index)
          nav_items <- c(
            nav_items,
            report_nav_item(
              figure_part_id(i, part),
              sprintf("3.%d.%d %s", i, part_index, part_title),
              "report-nav-h4"
            )
          )
          subparts <- part$subparts %||% list()
          part_figure_indices <- part$figure_indices
          part_figure_indices <- part_figure_indices[part_figure_indices <= length(sections[[i]]$figures)]
          if (length(subparts) > 0L) {
            for (subpart in subparts) {
              subpart_index <- suppressWarnings(as.integer(subpart$subpart_index %||% NA_integer_))
              if (!is.finite(subpart_index) || subpart_index < 1L) subpart_index <- 1L
              subpart_title <- subpart$title %||% paste0("Subpart ", subpart_index)
              nav_items <- c(
                nav_items,
                report_nav_item(
                  figure_subpart_id(i, part, subpart),
                  sprintf("3.%d.%d.%d %s", i, part_index, subpart_index, subpart_title),
                  "report-nav-h5"
                )
              )
              subpart_indices <- subpart$figure_indices
              subpart_indices <- subpart_indices[subpart_indices <= length(part_figure_indices)]
              if (length(subpart_indices) > 0L) {
                nav_items <- append_figure_nav_items(
                  nav_items,
                  sections[[i]],
                  i,
                  part_figure_indices[subpart_indices],
                  "report-nav-h6",
                  part_index = part_index,
                  subpart_index = subpart_index
                )
              }
            }
          } else if (length(part_figure_indices) > 0L) {
            nav_items <- append_figure_nav_items(
              nav_items,
              sections[[i]],
              i,
              part_figure_indices,
              "report-nav-h5",
              part_index = part_index
            )
          }
        }
      } else if (!length(direct_indices) && length(sections[[i]]$figures) > 0L) {
        nav_items <- append_figure_nav_items(
          nav_items,
          sections[[i]],
          i,
          seq_along(sections[[i]]$figures),
          "report-nav-h4"
        )
      }
    }
  }
  nav_items
}

build_fit_report_html <- function(params) {
  sections <- params$sections %||% list()
  parameter_sections <- params$parameter_sections %||% list(list(name = "Best Parameters", table = params$best_params))
  nav_items <- build_report_nav_items(sections, parameter_sections)
  parameter_blocks <- build_parameter_blocks_html(parameter_sections)

  section_blocks <- if (length(sections) == 0L) {
    '<p class="report-empty">No figures found for this fit.</p>'
  } else {
    paste(vapply(seq_along(sections), function(i) {
      section <- sections[[i]]
      slug <- paste0("section-", i)
      fig_blocks <- render_section_figure_blocks(section, i)
      paste0(
        '<section class="report-section" id="', slug, '">',
        '<h2>3.', i, " ", escape_html(section$name), "</h2>",
        fig_blocks,
        "</section>"
      )
    }, character(1)), collapse = "")
  }

  paste0(
    "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"/>",
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    "<title>O2 Supply-Demand MAP Fit Report</title>",
    "<style>",
    "body{margin:0;font-family:-apple-system,BlinkMacSystemFont,\"Segoe UI\",sans-serif;background:#f4f7fa;color:#1b2a38;}",
    ".report-shell{display:flex;gap:28px;max-width:1680px;margin:0 auto;padding:24px;}",
    ".report-sidebar{position:sticky;top:24px;align-self:flex-start;width:300px;max-height:calc(100vh - 48px);border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);overflow:auto;scrollbar-gutter:stable;}",
    ".report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}",
    ".report-kicker{font-size:11px;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;opacity:0.78;}",
    ".report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}",
    ".report-subtitle{margin-top:4px;font-size:12px;opacity:0.85;}",
    ".report-nav{padding:10px 8px 12px 8px;}.report-nav-list{margin:0;padding:0;list-style:none;}.report-nav-item{margin:4px 0;}",
    ".report-nav-link{display:block;padding:10px 12px;border-radius:8px;text-decoration:none;color:#17324c;font-size:14px;font-weight:600;line-height:1.35;}",
    ".report-nav-link:hover{background:rgba(47,110,164,0.08);}.report-nav-h3{padding-left:22px;font-size:13px;font-weight:500;}.report-nav-h4{padding-left:36px;font-size:12px;font-weight:500;color:#536577;}.report-nav-h5{padding-left:50px;font-size:12px;font-weight:500;color:#6a7a89;}.report-nav-h6{padding-left:64px;font-size:11px;font-weight:500;color:#7a8794;}",
    ".report-main{flex:1;min-width:0;max-width:1180px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}",
    ".report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-card h2,.report-section h2{margin-top:0;}.report-card h2,.report-card h3,.report-card h4,.report-card h5,.report-section h2,.report-section h3,.report-section h4,.report-section h5,.report-figure-part,.report-figure-subpart{scroll-margin-top:24px;}",
    ".report-meta{margin:0;color:#516274;font-size:14px;}",
    ".report-param-group{margin:22px 0 8px 0;padding-top:12px;border-top:1px solid #dfe6ee;color:#17324c;font-size:19px;}.report-param-section{margin:14px 0 8px 0;color:#284662;font-size:16px;}",
    ".report-table{width:100%;border-collapse:collapse;font-size:14px;}.report-table th,.report-table td{padding:10px 12px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}",
    ".report-table th{background:#f7f9fb;font-weight:700;}.report-empty{color:#657789;font-style:italic;}",
    ".report-figure-part{margin:22px 0 32px 0;}.report-figure-part h4{margin:0 0 6px 0;font-size:18px;}.report-figure-subpart{margin:18px 0 26px 0;}.report-figure-subpart h5{margin:0 0 6px 0;font-size:16px;color:#284662;}.report-part-description{margin:0 0 14px 0;color:#536577;}",
    ".report-figure-grid{display:grid;gap:18px;margin-bottom:26px;align-items:stretch;}.report-figure-grid--1{grid-template-columns:1fr;}.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}",
    ".report-figure-card{margin-bottom:26px;min-width:0;}.report-figure-grid .report-figure-card{margin-bottom:0;}.report-figure{margin:0 0 12px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}",
    ".report-figure-grid--2 .report-figure,.report-figure-grid--3 .report-figure{aspect-ratio:4/3;display:flex;align-items:center;justify-content:center;overflow:hidden;}.report-figure-image{width:100%;height:auto;display:block;border-radius:6px;}.report-figure-grid--2 .report-figure-image,.report-figure-grid--3 .report-figure-image{height:100%;object-fit:contain;}.report-figure-object{width:100%;min-height:680px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}.report-figure-blank{border:1px dashed #d7dee7;background:#f7f9fb;}",
    ".report-figure-grid--2 .report-figure-object,.report-figure-grid--3 .report-figure-object{height:100%;min-height:0;}.report-figure-fallback{padding:18px;text-align:center;}.report-figure-fallback a{color:#2f6ea4;font-weight:600;text-decoration:none;}.report-figure-title,.report-figure-legend{margin:8px 0 0 0;}.report-figure-grid .report-figure-title{font-size:13px;line-height:1.35;}.report-figure-grid .report-figure-legend{font-size:12px;line-height:1.4;color:#536577;}",
    "@media (max-width: 1100px){.report-shell{display:block;}.report-sidebar{position:static;width:auto;margin-bottom:20px;}}@media (max-width: 900px){.report-figure-grid--2,.report-figure-grid--3{grid-template-columns:1fr;}}",
    "</style></head><body>",
    '<div class="report-shell">',
    '<aside class="report-sidebar"><div class="report-sidebar-header">',
    '<div class="report-kicker">O2 Supply-Demand MAP</div>',
    '<div class="report-title">Fit Report</div>',
    '<div class="report-subtitle">', escape_html(params$fit_label), "</div></div>",
    '<nav class="report-nav"><ul class="report-nav-list">', paste(nav_items, collapse = ""), "</ul></nav></aside>",
    '<main class="report-main">',
    '<section class="report-card" id="report-metadata"><h1>O2 Supply-Demand MAP Fit Report</h1>',
    '<p class="report-meta"><strong>Fit Label:</strong> ', escape_html(params$fit_label), "<br/>",
    "<strong>Generated At:</strong> ", escape_html(params$generated_at), "</p></section>",
    '<section class="report-card" id="fit-summary"><h2>1. Fit Summary</h2>',
    fit_report_table_html(params$selected_summary),
    "</section>",
    '<section class="report-card" id="best-parameters"><h2>2. Parameters</h2>',
    parameter_blocks,
    "</section>",
    '<section class="report-card" id="figures"><h2>3. Figures</h2>',
    section_blocks,
    "</section></main></div></body></html>"
  )
}

render_one_fit <- function(fit_dir, template_path, out_subdir = "report", report_basename = "fit_report", render_pdf = FALSE) {
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  report_basename <- fit_report_basename_with_seed(report_basename, fit_dir)
  out_dir <- file.path(fit_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  section_specs <- stage_assets(build_section_specs(fit_dir))
  params <- list(
    fit_label = basename(fit_dir),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    selected_summary = read_fit_summary_selected(fit_dir),
    best_params = read_best_params(fit_dir),
    parameter_sections = read_parameter_sections(fit_dir),
    sections = section_specs
  )

  html_out <- file.path(out_dir, paste0(report_basename, ".html"))
  pdf_out <- NA_character_
  pdf_status <- "disabled"
  if (report_pandoc_available()) {
    html_out <- render(
      input = template_path,
      output_format = "html_document",
      output_file = paste0(report_basename, ".html"),
      output_dir = out_dir,
      params = params,
      envir = new.env(parent = baseenv()),
      quiet = TRUE
    )
    if (isTRUE(render_pdf) && report_pdflatex_available()) {
      pdf_status <- "rendered"
      pdf_out <- render(
        input = template_path,
        output_format = "pdf_document",
        output_file = paste0(report_basename, ".pdf"),
        output_dir = out_dir,
        params = params,
        envir = new.env(parent = baseenv()),
        quiet = TRUE
      )
    } else if (isTRUE(render_pdf)) {
      pdf_status <- "skipped_pdflatex_unavailable"
      message("  pdflatex unavailable: generated HTML report and skipped PDF.")
    } else {
      message("  PDF disabled: generated HTML report only.")
    }
  } else {
    pdf_status <- "skipped_pandoc_unavailable"
    writeLines(build_fit_report_html(params), con = html_out, useBytes = TRUE)
    message("  pandoc unavailable: generated HTML-only fallback report and skipped PDF.")
  }

  html_files_dir <- file.path(out_dir, paste0(report_basename, "_files"))
  if (dir.exists(html_files_dir)) {
    unlink(html_files_dir, recursive = TRUE, force = TRUE)
  }
  report_assets_dir <- file.path(out_dir, "assets")
  if (dir.exists(report_assets_dir)) {
    unlink(report_assets_dir, recursive = TRUE, force = TRUE)
  }

  list(
    fit_dir = fit_dir,
    out_dir = normalizePath(out_dir, mustWork = TRUE),
    html = normalizePath(html_out, mustWork = TRUE),
    pdf = if (is.na(pdf_out)) NA_character_ else normalizePath(pdf_out, mustWork = TRUE),
    pdf_status = pdf_status
  )
}

main <- function() {
  script_dir <- normalizePath(get_script_dir(), mustWork = TRUE)
  template_path <- file.path(script_dir, "fit_result_report.Rmd")
  if (!file.exists(template_path)) {
    stop("Missing report template: ", template_path)
  }

  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  fit_root <- argv$fit_dir %||% stop("Usage: render_fit_report.R --fit_dir=/abs/path/to/fit_or_root [--out_subdir=report] [--report_basename=fit_report] [--render_pdf=TRUE]")
  fit_root <- normalizePath(fit_root, mustWork = TRUE)
  out_subdir <- argv$out_subdir %||% "report"
  report_basename <- argv$report_basename %||% "fit_report"
  render_pdf <- report_pdf_enabled(argv)

  fit_dirs <- find_fit_dirs_under(fit_root)
  if (length(fit_dirs) == 0) {
    stop("No fit directories found under: ", fit_root)
  }

  message("Found ", length(fit_dirs), " fit director", if (length(fit_dirs) == 1) "y" else "ies", " under: ", fit_root)
  for (i in seq_along(fit_dirs)) {
    fit_dir <- fit_dirs[[i]]
    message("[", i, "/", length(fit_dirs), "] Rendering report for: ", fit_dir)
    res <- render_one_fit(
      fit_dir = fit_dir,
      template_path = template_path,
      out_subdir = out_subdir,
      report_basename = report_basename,
      render_pdf = render_pdf
    )
    message("  HTML: ", res$html)
    if (is.na(res$pdf)) {
      if (identical(res$pdf_status, "disabled")) {
        message("  PDF : skipped (disabled)")
      } else
      if (identical(res$pdf_status, "skipped_pdflatex_unavailable")) {
        message("  PDF : skipped (pdflatex unavailable)")
      } else {
        message("  PDF : skipped (pandoc unavailable)")
      }
    } else {
      message("  PDF : ", res$pdf)
    }
  }
}

if (sys.nframe() == 0) {
  main()
}
