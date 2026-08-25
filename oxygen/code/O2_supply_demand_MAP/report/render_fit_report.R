#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

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
source(file.path(REPORT_WORKFLOW_ROOT, "util", "o2_supply_demand_map_html_utils.R"), local = environment())
source(file.path(REPORT_WORKFLOW_ROOT, "util", "o2_supply_demand_map_report_utils.R"), local = environment())
source(file.path(REPORT_SCRIPT_DIR, "run_provenance_report.R"), local = environment())

`%||%` <- o2sd_report_null_coalesce_simple

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
      "objective_invivo_necrosis",
      "objective_invitro",
      "objective_soft_coupling",
      "objective_constraints",
      "joint_weight_invivo",
      "joint_weight_invitro",
      "joint_invitro_growth_weight",
      "joint_invitro_ploidy_weight",
      "joint_invitro_flow_weight",
      "joint_soft_coupling_enabled",
      "joint_soft_coupling_sigma_default",
      "joint_soft_coupling_welsch_c",
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
      "objective_necrosis",
      "objective_burden_neg2loglik_raw",
      "objective_ploidy_neg2loglik_raw",
      "objective_necrosis_neg2loglik_raw"
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

summary_flag_true <- o2sd_flag_true

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

annotate_parameter_table_for_report <- o2sd_report_annotate_parameter_table_for_report

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
  if ("ratio_vivo_to_vitro" %in% names(tab)) {
    if (!"fold_change_vivo_to_vitro" %in% names(tab)) {
      tab$fold_change_vivo_to_vitro <- tab$ratio_vivo_to_vitro
    } else {
      missing_fold_change <- is.na(suppressWarnings(as.numeric(tab$fold_change_vivo_to_vitro)))
      tab$fold_change_vivo_to_vitro[missing_fold_change] <- tab$ratio_vivo_to_vitro[missing_fold_change]
    }
  }
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
    "regularization_sigma",
    "penalty_type",
    "welsch_c",
    "welsch_transition_delta",
    "abs_delta_over_sigma",
    "welsch_saturation_fraction",
    "penalty_region",
    "penalty_paid",
    "joint_union_lower_transformed",
    "joint_union_upper_transformed",
    "feasible_at_solution"
  )
  keep <- intersect(keep, names(tab))
  tab[, keep, drop = FALSE]
}

report_truthy <- o2sd_report_truthy
transformed_param_to_natural <- o2sd_report_transformed_param_to_natural

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

report_pandoc_available <- o2sd_report_pandoc_available
report_pdflatex_available <- o2sd_report_pdflatex_available

report_pdf_enabled <- function(argv = list()) {
  arg_value <- argv$render_pdf %||% Sys.getenv("O2_RENDER_PDF", unset = "")
  if (!nzchar(arg_value)) {
    return(FALSE)
  }
  tolower(arg_value) %in% c("true", "t", "1", "yes", "y")
}

escape_html <- o2sd_html_escape_full

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

parameter_figure_id <- function(index, title) {
  sprintf("params-figure-%02d-%s", index, html_id_slug(title))
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

stage_parameter_figures <- function(parameter_figures) {
  if (!length(parameter_figures)) return(list())
  staged <- stage_assets(list(list(name = "Parameter Figures", figures = parameter_figures)))
  staged[[1]]$figures %||% list()
}

build_parameter_figure_specs <- function(fit_dir) {
  if (!is_joint_fit_dir_for_report(fit_dir)) return(list())
  viz_dir <- file.path(fit_dir, "viz", "joint_parameters")
  if (!dir.exists(viz_dir)) return(list())
  candidates <- list.files(
    viz_dir,
    pattern = "_ratio_vivo_to_vitro_all_soft[.]pdf$",
    full.names = TRUE
  )
  if (!length(candidates)) return(list())
  info <- file.info(candidates)
  pdf_path <- candidates[[order(info$mtime, decreasing = TRUE)[[1L]]]]
  list(make_figure_spec(
    pdf_path,
    "In Vivo vs In Vitro Parameter Ratios",
    paste0(
      "The upstream joint-parameter analysis materializes the in vivo/in vitro ",
      "ratios consumed by this figure. Point and segment color show the direction ",
      "of difference; label color shows mechanism class."
    )
  ))
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
          "Forward simulation from day 0 to %s summarizing predicted total burden on a log10 scale, ",
          "predicted viable cells on a log10 scale, the %s fraction among all dead biomass, ",
          "mean ploidy, chromosome-number density, and effective-oxygen trajectories."
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
          "Forward simulation from day 0 to 1000 using the chromosome-number probability heatmap as the main panel. ",
          "Column annotations show cohort-level burden, viable cells, and effective O2 (%) in 2N/4N rows; ",
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

report_magick_available <- o2sd_report_magick_available
report_gs_available <- o2sd_report_gs_available
report_base64enc_available <- o2sd_report_base64enc_available
render_pdf_preview_png_gs <- o2sd_report_render_pdf_preview_png_gs
file_to_data_uri <- o2sd_report_file_to_data_uri

pdf_to_data_uri <- function(pdf_path) {
  file_to_data_uri(pdf_path, mime = "application/pdf")
}

resource_death_language_report <- function() {
  list(
    report_phrase = "hypoxia-origin dead",
    figure_phrase = "hypoxia-origin dead",
    cin_report_phrase = "CIN-associated dead"
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
    "Predicted Total Burden Viable/Dead Decomposition (0-100/300/1000 day)",
    paste0(
      "Forward simulations from day 0 to 100, 300, and 1000 showing the predicted viable/",
      death_language$figure_phrase,
      "/",
      death_language$cin_report_phrase,
      " burden decomposition with y-axis ticks shown as log10 burden values and a shared component legend."
    )
  )
  if (length(burden_predict_figs) == 0L && length(burden_predict) > 0L) {
    burden_predict_figs <- optional_series_figures(
      burden_predict,
      "Predicted Total Burden Viable/Dead Decomposition (0-%s day)",
      paste0(
        "Forward simulation from day 0 to %s showing the predicted viable/",
        death_language$figure_phrase,
        "/",
        death_language$cin_report_phrase,
        " burden decomposition with y-axis ticks shown as log10 burden values."
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
      "Total Tumor Burden Viable/Dead Decomposition",
      paste0(
        "Observed and fitted tumor burden decomposed into viable, ",
        death_language$report_phrase,
        ", and ",
        death_language$cin_report_phrase,
        " components."
      )
    ),
    burden_predict_figs
  ))

  ploidy_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "ploidy_weighted_mean_over_time.pdf",
      "Weighted Mean Ploidy Over Time",
      "Weighted mean chromosome-number trajectory over time under the fitted model."
    ),
    optional_figure(
      viz_dir,
      "terminal_ploidy_observed_vs_predicted_violin.pdf",
      "Terminal Chromosome-Number Observed vs Predicted",
      "Observed and predicted terminal chromosome-number distributions at harvest."
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
      "Oxygen Target vs Effective Oxygen Time Course",
      "Time-course comparison between the oxygen supply-demand target and the lagged effective oxygen state used by the model."
    ),
    optional_figure(
      viz_dir,
      "missegregation_probability_over_time.pdf",
      "Mean Per-Chromosome Missegregation Probability Over Time",
      "Viable-population-weighted mean per-chromosome missegregation probability over time, computed from the fitted effective-O2 trajectory and predicted chromosome-number distribution."
    ),
    optional_figure(
      viz_dir,
      "predict_burden_vs_o2.pdf",
      "Predicted Tumor Burden vs Effective Oxygen",
      "Forward-simulation burden trajectories plotted against the effective oxygen state."
    )
  ))
  oxygen_observation_missegregation_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "population_average_cin_by_initial_cohort_day100.pdf",
      "0-100 Day Population-average CIN rate over time",
      "Population-average CIN rate over days 0-100, shown for canonical 2N-derived and 4N-derived in vivo trajectories."
    ),
    optional_figure(
      viz_dir,
      "population_average_cin_by_initial_cohort_day300.pdf",
      "0-300 Day Population-average CIN rate over time",
      "Population-average CIN rate over days 0-300, shown for canonical 2N-derived and 4N-derived in vivo trajectories."
    ),
    optional_figure(
      viz_dir,
      "population_average_cin_by_initial_cohort_day1000.pdf",
      "0-1000 Day Population-average CIN rate over time",
      "Population-average CIN rate over days 0-1000, shown for canonical 2N-derived and 4N-derived in vivo trajectories."
    )
  ))
  oxygen_ms_relationship_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "ms_rate_vs_nonviable_daughter_fraction.pdf",
      "Nonviable Daughter Fraction vs Missegregation Rate",
      "Per-division fraction of daughter cells that are nonviable because of missegregation-linked loss, shown against missegregation rate."
    ),
    optional_figure(
      viz_dir,
      "ploidy_vs_viability_after_ms.pdf",
      "Ploidy-Dependent Post-Missegregation Survival",
      "Post-missegregation survival after a one-copy-loss event across the ploidy grid."
    ),
    optional_figure(
      viz_dir,
      "death_rate_vs_missegregation_rate.pdf",
      "Stress-Associated Death Rate vs Missegregation Rate",
      "Missegregation-rate curve plotted against the fitted stress-associated death rate at the 2N and 4N reference ploidy states."
    )
  ))
  oxygen_in_vivo_relationship_figs <- Filter(Negate(is.null), c(
    optional_figure(
      viz_dir,
      "oxygen_vs_missegregation_rate_multi_ploidy.pdf",
      "Effective Oxygen vs Missegregation Rate",
      "In vivo missegregation-rate curves across effective oxygen levels for multiple reference ploidy states."
    ),
    optional_figure(
      viz_dir,
      "oxygen_vs_proliferation_rate.pdf",
      "Effective Oxygen vs Proliferation Rate",
      "In vivo proliferation-rate curves across effective oxygen levels for fitted reference states."
    ),
    optional_figure(
      viz_dir,
      "oxygen_vs_death_rate.pdf",
      "Effective Oxygen vs Stress-Associated Death Rate",
      "In vivo stress-associated death-rate curves across effective oxygen levels for fitted reference states."
    )
  ))
  figure_index_range <- function(start, figs) {
    if (!length(figs)) return(integer(0))
    seq.int(start, length.out = length(figs))
  }
  oxygen_idx_start <- 1L
  oxygen_dynamics_idx <- figure_index_range(oxygen_idx_start, oxygen_dynamics_figs)
  oxygen_idx_start <- oxygen_idx_start + length(oxygen_dynamics_figs)
  oxygen_observation_missegregation_idx <- figure_index_range(oxygen_idx_start, oxygen_observation_missegregation_figs)
  oxygen_idx_start <- oxygen_idx_start + length(oxygen_observation_missegregation_figs)
  oxygen_ms_relationship_idx <- figure_index_range(oxygen_idx_start, oxygen_ms_relationship_figs)
  oxygen_idx_start <- oxygen_idx_start + length(oxygen_ms_relationship_figs)
  oxygen_in_vivo_relationship_idx <- figure_index_range(oxygen_idx_start, oxygen_in_vivo_relationship_figs)
  oxygen_figs <- c(
    oxygen_dynamics_figs,
    oxygen_observation_missegregation_figs,
    oxygen_ms_relationship_figs,
    oxygen_in_vivo_relationship_figs
  )
  oxygen_figure_parts <- Filter(Negate(is.null), list(
    if (length(oxygen_observation_missegregation_idx)) {
      list(
        part_index = 3L,
        title = "Population-average CIN rate over time",
        description = "In vivo population-average CIN rate over time for the 0-100, 0-300, and 0-1000 day prediction horizons, restricted to canonical 2N-derived and 4N-derived trajectories.",
        figure_indices = oxygen_observation_missegregation_idx,
        cols = 3L
      )
    },
    if (length(oxygen_ms_relationship_idx)) {
      list(
        part_index = 4L,
        title = "Missegregation-Linked Survival and Death Relationships",
        description = "Missegregation-linked post-missegregation survival, nonviable daughter production, and stress-associated death-rate relationships.",
        figure_indices = oxygen_ms_relationship_idx,
        cols = 3L
      )
    },
    if (length(oxygen_in_vivo_relationship_idx)) {
      list(
        part_index = 5L,
        title = "In Vivo Effective-Oxygen Rate Relationships",
        description = "In vivo-only effective-oxygen rate relationships corresponding to the In Vivo vs In Vitro comparison views.",
        figure_indices = oxygen_in_vivo_relationship_idx,
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
      name = "Effective Oxygen / Resource Stress",
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
          "Assigned Fixed Oxygen and Selected-Day Viable Cells",
          "Branch-aware diagnostic for the in vitro runner. Each cohort is split into control/deprived lineage panels using the same branch-specific x-axis as the aligned growth/chromosome-number/burden composite; repeated lineage passages are not averaged across branches. The upper row shows assigned fixed oxygen and the lower row shows selected-day predicted viable cells."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_rate_function_diagnostics.pdf",
          "Rate-Function Diagnostics",
          "Best-fit fixed-oxygen, chromosome-number, stress-associated death, proliferation, and missegregation rate functions for the current seed."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_missegregation_probability_over_passage.pdf",
          "Mean Per-Chromosome Missegregation Probability Over Passage",
          "Viable-population-weighted mean per-chromosome missegregation probability across in vitro passage branches, computed from the fitted fixed-oxygen levels and selected-day chromosome-number distributions."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_daily_counts.pdf",
          "Daily Viable-Cell Trajectories",
          "Predicted viable-cell trajectories split into 2N/control, 2N/deprived, 4N/control, and 4N/deprived panels, with each lineage passage shown as an inset subplot; selected propagation days are marked."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_growth_ploidy_burden_composite.pdf",
          "Aligned Growth, Chromosome-Number, and Burden Fit",
          "Composite in vitro fit view. The 2N and 4N cohort blocks are stacked vertically; each block contains growth rate, chromosome-number quantile, and burden-decomposition rows. Repeated lineage passages are split into branch-specific fixed-oxygen x-axis labels rather than averaged together; growth-rate and chromosome-number lines follow parent-child lineage links, so parallel p10 branches both connect to their shared p9 parent. Observed growth-rate points are drawn at their own branch positions. In vitro burden components are viable/dead fractions normalized by the displayed predicted cell components, so the burden row ranges from 0 to 1. Control and deprived panels use their own passage ranges, with rows aligned within each lineage panel."
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_flow_density.pdf",
          "Flow-Density Fit",
          "Observed G0/G1 chromosome-number density curves are overlaid with the fitted flow-density prediction.",
          layout_group = "density-distribution"
        ),
        optional_figure_with_layout(
          viz_dir,
          "invitro_distribution_heatmap.pdf",
          "Predicted Chromosome-Number Distribution",
          "Full predicted chromosome-number distribution across in vitro passages.",
          layout_group = "density-distribution"
        )
      )
    ),
    list(
      name = "Missegregation-Linked Survival and Death Relationships",
      figures = c(
        optional_figure_with_layout(
          viz_dir,
          "ms_rate_vs_nonviable_daughter_fraction.pdf",
          "Nonviable Daughter Fraction vs Missegregation Rate",
          "Per-division fraction of daughter cells that are nonviable because of missegregation-linked loss, shown against missegregation rate.",
          figure_index = 3L
        ),
        optional_figure_with_layout(
          viz_dir,
          "death_rate_vs_missegregation_rate.pdf",
          "Stress-Associated Death Rate vs Missegregation Rate",
          "Missegregation-rate curve plotted against the fitted stress-associated death rate at the 2N and 4N reference ploidy states.",
          figure_index = 2L
        ),
        optional_figure_with_layout(
          viz_dir,
          "ploidy_vs_viability_after_ms.pdf",
          "Ploidy-Dependent Post-Missegregation Survival",
          "Post-missegregation survival after a one-copy-loss event across the ploidy grid.",
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
        title = if (identical(source_section$name %||% "", "Effective Oxygen / Resource Stress")) "Effective Oxygen Dynamics" else "Overview",
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

build_invivo_invitro_comparison_section <- function(fit_dir) {
  viz_dir <- file.path(fit_dir, "viz", "invivo_vs_invitro")
  if (!dir.exists(viz_dir)) return(NULL)
  figures <- c(
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ms_rate_vs_nonviable_daughter_fraction.pdf",
      "Nonviable Daughter Fraction vs Missegregation Rate",
      "Left panel uses the in vivo functional-response output; right panel uses the in vitro functional-response output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ploidy_vs_viability_after_ms.pdf",
      "Ploidy-Dependent Post-Missegregation Survival",
      "Left panel uses the in vivo viability curve; right panel uses the in vitro viability curve."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_death_rate_vs_missegregation_rate.pdf",
      "Stress-Associated Death Rate vs Missegregation Rate",
      "Left panel uses the in vivo functional-response output; right panel uses the in vitro functional-response output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ploidy_vs_death_rate_by_o2.pdf",
      "Ploidy vs Stress-Associated Death Rate by Oxygen Level",
      "Left panel uses the in vivo effective-oxygen resource-stress output; right panel uses the in vitro fixed-oxygen resource-stress output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_ploidy_vs_proliferation_rate_by_o2.pdf",
      "Ploidy vs Proliferation Rate by Oxygen Level",
      "Left panel uses the in vivo effective-oxygen resource-stress output; right panel uses the in vitro fixed-oxygen resource-stress output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_oxygen_vs_missegregation_rate_multi_ploidy.pdf",
      "Oxygen Level vs Missegregation Rate",
      "Left panel uses the in vivo effective-oxygen resource-stress output; right panel uses the in vitro fixed-oxygen resource-stress output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_oxygen_vs_proliferation_rate.pdf",
      "Oxygen Level vs Proliferation Rate",
      "Left panel uses the in vivo effective-oxygen resource-stress output; right panel uses the in vitro fixed-oxygen resource-stress output."
    ),
    optional_figure(
      viz_dir,
      "invivo_vs_invitro_oxygen_vs_death_rate.pdf",
      "Oxygen Level vs Stress-Associated Death Rate",
      "Left panel uses the in vivo effective-oxygen resource-stress output; right panel uses the in vitro fixed-oxygen resource-stress output."
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

build_parameter_figure_blocks_html <- function(parameter_figures) {
  if (!length(parameter_figures)) return("")
  paste(vapply(seq_along(parameter_figures), function(i) {
    fig <- parameter_figures[[i]]
    title <- fit_report_fig_title(fig)
    legend <- fit_report_fig_legend(fig)
    media <- if (identical(fig$html_embed_kind %||% "img", "pdf_object")) {
      sprintf(
        paste0(
          '<div class="report-parameter-figure">',
          '<object data="%s" type="application/pdf" class="report-parameter-figure-object">',
          '<div class="report-figure-fallback"><a href="%s">Open PDF figure</a></div>',
          '</object></div>'
        ),
        fit_report_figure_data_uri(fig),
        escape_html(basename(fig$pdf_asset_abs[[1]]))
      )
    } else {
      sprintf(
        '<div class="report-parameter-figure"><img src="%s" alt="%s" class="report-parameter-figure-image"/></div>',
        fit_report_figure_data_uri(fig),
        escape_html(title)
      )
    }
    paste0(
      '<article class="report-parameter-figure-card" id="', parameter_figure_id(i, title), '">',
      '<h3 class="report-param-figure-title">', escape_html(title), "</h3>",
      media,
      '<p class="report-figure-legend">', escape_html(legend), "</p>",
      "</article>"
    )
  }, character(1)), collapse = "")
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

build_report_nav_items <- function(sections, parameter_sections, parameter_figures = list()) {
  nav_items <- c(
    report_nav_item("report-metadata", "Report Metadata", "report-nav-h2"),
    report_nav_item("fit-summary", "1. Fit Summary", "report-nav-h2"),
    report_nav_item("best-parameters", "2. Parameters", "report-nav-h2")
  )
  if (length(parameter_figures) > 0L) {
    for (i in seq_along(parameter_figures)) {
      title <- fit_report_fig_title(parameter_figures[[i]])
      nav_items <- c(nav_items, report_nav_item(parameter_figure_id(i, title), title, "report-nav-h3"))
    }
  }
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
  parameter_figures <- params$parameter_figures %||% list()
  nav_items <- build_report_nav_items(sections, parameter_sections, parameter_figures = parameter_figures)
  nav_items <- c(nav_items, '<li class="report-nav-item"><a class="report-nav-link report-nav-h3" href="#run-provenance">4. Run Provenance</a></li>')
  parameter_figure_blocks <- build_parameter_figure_blocks_html(parameter_figures)
  parameter_blocks <- build_parameter_blocks_html(parameter_sections)
  provenance_block <- o2sd_run_provenance_html(params$fit_dir, title = "4. Run Provenance", section_id = "run-provenance")

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
    "<title>Resource-Stress Chromosome-Number Evolution Fit Report</title>",
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
    ".report-param-group{margin:22px 0 8px 0;padding-top:12px;border-top:1px solid #dfe6ee;color:#17324c;font-size:19px;}.report-param-section{margin:14px 0 8px 0;color:#284662;font-size:16px;}.report-param-figure-title{margin:0 0 10px 0;color:#284662;font-size:16px;}",
    ".report-table{width:100%;border-collapse:collapse;font-size:14px;}.report-table th,.report-table td{padding:10px 12px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}",
    ".report-table th{background:#f7f9fb;font-weight:700;}.report-empty{color:#657789;font-style:italic;}",
    ".report-figure-part{margin:22px 0 32px 0;}.report-figure-part h4{margin:0 0 6px 0;font-size:18px;}.report-figure-subpart{margin:18px 0 26px 0;}.report-figure-subpart h5{margin:0 0 6px 0;font-size:16px;color:#284662;}.report-part-description{margin:0 0 14px 0;color:#536577;}",
    ".report-figure-grid{display:grid;gap:18px;margin-bottom:26px;align-items:stretch;}.report-figure-grid--1{grid-template-columns:1fr;}.report-figure-grid--2{grid-template-columns:repeat(2,minmax(0,1fr));}.report-figure-grid--3{grid-template-columns:repeat(3,minmax(0,1fr));}",
    ".report-figure-card{margin-bottom:26px;min-width:0;}.report-figure-grid .report-figure-card{margin-bottom:0;}.report-figure{margin:0 0 12px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}",
    ".report-figure-grid--2 .report-figure,.report-figure-grid--3 .report-figure{aspect-ratio:4/3;display:flex;align-items:center;justify-content:center;overflow:hidden;}.report-figure-image{width:100%;height:auto;display:block;border-radius:6px;}.report-figure-grid--2 .report-figure-image,.report-figure-grid--3 .report-figure-image{height:100%;object-fit:contain;}.report-figure-object{width:100%;min-height:680px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}.report-figure-blank{border:1px dashed #d7dee7;background:#f7f9fb;}.report-parameter-figure-card{margin:2px 0 24px 0;}.report-parameter-figure{aspect-ratio:18/6.5;display:flex;align-items:center;justify-content:center;overflow:hidden;margin:0 0 10px 0;padding:10px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}.report-parameter-figure-image{width:100%;height:100%;display:block;object-fit:contain;border-radius:6px;}.report-parameter-figure-object{width:100%;height:100%;min-height:420px;border:0;background:#fff;}",
    ".report-figure-grid--2 .report-figure-object,.report-figure-grid--3 .report-figure-object{height:100%;min-height:0;}.report-figure-fallback{padding:18px;text-align:center;}.report-figure-fallback a{color:#2f6ea4;font-weight:600;text-decoration:none;}.report-figure-title,.report-figure-legend{margin:8px 0 0 0;}.report-figure-grid .report-figure-title{font-size:13px;line-height:1.35;}.report-figure-grid .report-figure-legend{font-size:12px;line-height:1.4;color:#536577;}",
    "@media (max-width: 1100px){.report-shell{display:block;}.report-sidebar{position:static;width:auto;margin-bottom:20px;}}@media (max-width: 900px){.report-figure-grid--2,.report-figure-grid--3{grid-template-columns:1fr;}}",
    "</style></head><body>",
    '<div class="report-shell">',
    '<aside class="report-sidebar"><div class="report-sidebar-header">',
    '<div class="report-kicker">Resource-Stress Model</div>',
    '<div class="report-title">Fit Report</div>',
    '<div class="report-subtitle">', escape_html(params$fit_label), "</div></div>",
    '<nav class="report-nav"><ul class="report-nav-list">', paste(nav_items, collapse = ""), "</ul></nav></aside>",
    '<main class="report-main">',
    '<section class="report-card" id="report-metadata"><h1>Resource-Stress Chromosome-Number Evolution Fit Report</h1>',
    '<p class="report-meta"><strong>Fit Label:</strong> ', escape_html(params$fit_label), "<br/>",
    "<strong>Generated At:</strong> ", escape_html(params$generated_at), "</p></section>",
    '<section class="report-card" id="fit-summary"><h2>1. Fit Summary</h2>',
    fit_report_table_html(params$selected_summary),
    "</section>",
    '<section class="report-card" id="best-parameters"><h2>2. Parameters</h2>',
    parameter_figure_blocks,
    parameter_blocks,
    "</section>",
    '<section class="report-card" id="figures"><h2>3. Figures</h2>',
    section_blocks,
    "</section>",
    provenance_block,
    "</main></div></body></html>"
  )
}

render_one_fit <- function(fit_dir, template_path, out_subdir = "report", report_basename = "fit_report", render_pdf = FALSE) {
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  report_basename <- fit_report_basename_with_seed(report_basename, fit_dir)
  out_dir <- file.path(fit_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  section_specs <- stage_assets(build_section_specs(fit_dir))
  parameter_figures <- stage_parameter_figures(build_parameter_figure_specs(fit_dir))
  params <- list(
    fit_label = basename(fit_dir),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    fit_dir = fit_dir,
    selected_summary = read_fit_summary_selected(fit_dir),
    best_params = read_best_params(fit_dir),
    parameter_sections = read_parameter_sections(fit_dir),
    parameter_figures = parameter_figures,
    sections = section_specs,
    provenance_html = o2sd_run_provenance_html(fit_dir, title = "4. Run Provenance", section_id = "run-provenance")
  )

  html_out <- file.path(out_dir, paste0(report_basename, ".html"))
  pdf_out <- NA_character_
  pdf_status <- "disabled"
  if (report_pandoc_available()) {
    intermediates_dir <- file.path(
      out_dir,
      paste0(".", report_basename, "_intermediates")
    )
    dir.create(intermediates_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(
      unlink(intermediates_dir, recursive = TRUE, force = TRUE),
      add = TRUE
    )
    html_out <- render(
      input = template_path,
      output_format = "html_document",
      output_file = paste0(report_basename, ".html"),
      output_dir = out_dir,
      intermediates_dir = intermediates_dir,
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
        intermediates_dir = intermediates_dir,
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

  o2sd_inject_report_image_lightbox(html_out)

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
