.joint_parameter_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[
    basename(frame_files) == "o2_supply_demand_map_joint_parameter_utils.R"
  ]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(.joint_parameter_utils_dir, "o2_supply_demand_map_common_semantics.R"),
  local = TRUE,
  chdir = TRUE
)
rm(.joint_parameter_utils_dir)

joint_parameter_ratio_null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

`%||%` <- joint_parameter_ratio_null_coalesce

joint_parameter_ratio_param_specs <- function() {
  list(
    O2_crit = list(transform = "log10"),
    mu_hp = list(transform = "log10"),
    p_misseg = list(transform = "log10"),
    k_o_mis = list(transform = "log10"),
    buffer_smax = list(transform = "identity"),
    buffer_beta = list(transform = "log10"),
    buffer_n_exp = list(transform = "log10"),
    n_O = list(transform = "identity"),
    lam_max = list(transform = "log10"),
    p_mis_base = list(transform = "log10"),
    p_wgd = list(transform = "log10"),
    gamma_mu = list(transform = "identity"),
    alpha_o2 = list(transform = "log10"),
    gamma_growth = list(transform = "identity"),
    o2_S0 = list(transform = "log10"),
    kappa_O = list(transform = "log10"),
    eta_o2 = list(transform = "log10"),
    rho_2N = list(transform = "log10"),
    k_clear = list(transform = "log10"),
    sigma_burden = list(transform = "log10")
  )
}

joint_parameter_ratio_plot_defaults <- o2sd_joint_default_soft_coupling_params

joint_parameter_ratio_mechanism_map <- function() {
  data.frame(
    parameter = c(
      "p_mis_base", "p_wgd", "k_o_mis", "p_misseg",
      "gamma_mu", "mu_hp",
      "O2_crit", "n_O",
      "alpha_o2", "gamma_growth", "lam_max",
      "buffer_beta", "buffer_smax", "buffer_n_exp"
    ),
    mechanism_class = c(
      rep("Baseline growth", 4),
      rep("MS / WGD generation", 2),
      rep("Oxygen-death response", 2),
      rep("Oxygen-growth response", 3),
      rep("Post-MS survival", 3)
    ),
    stringsAsFactors = FALSE
  )
}

joint_parameter_ratio_mechanism_levels <- function() {
  c(
    "Baseline growth",
    "MS / WGD generation",
    "Oxygen-death response",
    "Oxygen-growth response",
    "Post-MS survival"
  )
}

joint_parameter_ratio_mechanism_colors <- function() {
  c(
    "Baseline growth" = "#E7298A",
    "MS / WGD generation" = "#7570B3",
    "Oxygen-death response" = "#D95F02",
    "Oxygen-growth response" = "#1B9E77",
    "Post-MS survival" = "#66A61E"
  )
}

joint_parameter_ratio_direction_colors <- function() {
  c(
    "Higher in vitro" = "#2F75B5",
    "Near 1x" = "#737D87",
    "Higher in vivo" = "#C23B50"
  )
}

joint_parameter_ratio_read_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  tryCatch(
    {
      if (identical(ext, "csv")) {
        utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
      } else {
        utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
      }
    },
    error = function(e) NULL
  )
}

joint_parameter_ratio_read_summary <- function(fit_dir) {
  tab <- joint_parameter_ratio_read_table(file.path(fit_dir, "fit_summary.tsv"))
  if (!is.data.frame(tab) || !all(c("metric", "value") %in% names(tab))) return(character())
  vals <- as.character(tab$value)
  names(vals) <- as.character(tab$metric)
  vals
}

joint_parameter_ratio_summary_value <- function(summary_map, key, default = NA_character_) {
  if (is.null(summary_map) || !length(summary_map) || !(key %in% names(summary_map))) return(default)
  if (is.na(summary_map[[key]])) return(default)
  value <- as.character(summary_map[[key]])
  if (!nzchar(trimws(value))) default else value
}

joint_parameter_ratio_split_csv <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(trimws(x[[1]]))) return(character())
  out <- trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

joint_parameter_ratio_active_params <- function(summary_map) {
  params <- joint_parameter_ratio_split_csv(
    joint_parameter_ratio_summary_value(summary_map, "joint_soft_coupling_params", "")
  )
  if (!length(params)) params <- joint_parameter_ratio_plot_defaults()
  unique(params)
}

joint_parameter_ratio_sanitize_stub <- function(x) {
  x <- if (is.null(x) || !length(x) || is.na(x[[1]])) "fit_report" else as.character(x[[1]])
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "fit_report" else x
}

joint_parameter_ratio_seed_token <- function(fit_dir, summary_map = NULL) {
  seed_value <- joint_parameter_ratio_summary_value(summary_map, "seed", "")
  if (nzchar(seed_value)) {
    seed_num <- suppressWarnings(as.integer(seed_value))
    if (is.finite(seed_num)) return(paste0("seed", seed_num))
  }

  base <- basename(normalizePath(fit_dir, mustWork = FALSE))
  hit <- regmatches(base, regexpr("seed[0-9]+", base, ignore.case = TRUE))
  if (length(hit) && nzchar(hit[[1]])) return(tolower(hit[[1]]))
  joint_parameter_ratio_sanitize_stub(base)
}

joint_parameter_ratio_fit_report_label <- function(fit_dir, summary_map = NULL) {
  paste0("fit_report_", joint_parameter_ratio_seed_token(fit_dir, summary_map), ".html")
}

joint_parameter_ratio_output_basename <- function(fit_dir) {
  summary_map <- joint_parameter_ratio_read_summary(fit_dir)
  paste0(
    "fit_report_",
    joint_parameter_ratio_seed_token(fit_dir, summary_map),
    "_ratio_vivo_to_vitro_all_soft"
  )
}

joint_parameter_ratio_apply_transform <- function(x, transform) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return(NA_real_)
  if (identical(transform, "log10")) return(log10(x))
  if (identical(transform, "logit")) return(stats::qlogis(x))
  x
}

joint_parameter_ratio_inverse_transform <- function(x, transform) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return(NA_real_)
  if (identical(transform, "log10")) return(10^x)
  if (identical(transform, "logit")) return(stats::plogis(x))
  x
}

joint_parameter_ratio_natural_name <- function(x) {
  x <- as.character(x)
  x <- sub("^delta__", "", x)
  x <- sub("^ivt__", "", x)
  x <- sub("^log10_", "", x)
  x <- sub("^logit_", "", x)
  x
}

joint_parameter_ratio_infer_transform <- function(param_name, scale = NA_character_) {
  scale <- if (is.null(scale) || !length(scale) || is.na(scale[[1]])) NA_character_ else as.character(scale[[1]])
  if (!is.na(scale) && nzchar(scale)) return(scale)
  param_name <- if (is.null(param_name) || !length(param_name) || is.na(param_name[[1]])) "" else as.character(param_name[[1]])
  if (grepl("^delta__log10_|^log10_", param_name)) return("log10")
  if (grepl("^delta__logit_|^logit_", param_name)) return("logit")
  "identity"
}

joint_parameter_ratio_find_best_params_path <- function(path) {
  if (is.null(path) || !length(path) || is.na(path[[1]]) || !nzchar(path[[1]])) return(character())
  path <- normalizePath(as.character(path[[1]]), mustWork = FALSE)
  if (dir.exists(path)) return(file.path(path, "best_params.tsv"))
  unique(c(
    if (identical(basename(path), "best_params.tsv")) path else character(),
    file.path(dirname(path), "best_params.tsv")
  ))
}

joint_parameter_ratio_seed_candidates <- function(summary_map, kind) {
  if (identical(kind, "invivo")) {
    keys <- c("joint_warmup_invivo_seed_dir", "joint_warmup_invivo_best_seed_dir", "joint_warmup_invivo_source_path")
  } else {
    keys <- c("joint_warmup_invitro_seed_dir", "joint_warmup_invitro_best_seed_dir", "joint_warmup_invitro_source_path")
  }
  paths <- unlist(lapply(keys, function(key) {
    joint_parameter_ratio_find_best_params_path(
      joint_parameter_ratio_summary_value(summary_map, key, "")
    )
  }), use.names = FALSE)
  unique(paths[nzchar(paths)])
}

joint_parameter_ratio_read_best_values <- function(paths) {
  paths <- unique(as.character(paths))
  paths <- paths[nzchar(paths)]
  for (path in paths) {
    if (!file.exists(path)) next
    tab <- joint_parameter_ratio_read_table(path)
    if (!is.data.frame(tab) || !all(c("parameter", "value") %in% names(tab))) next
    params <- as.character(tab$parameter)
    vals <- suppressWarnings(as.numeric(tab$value))
    keep <- nzchar(params) & is.finite(vals)
    if (!any(keep)) next
    out <- vals[keep]
    names(out) <- params[keep]
    attr(out, "source_path") <- normalizePath(path, mustWork = FALSE)
    return(out)
  }
  numeric()
}

joint_parameter_ratio_rows_from_maps <- function(invivo_vals, invitro_vals, active_params) {
  if (!length(invivo_vals) || !length(invitro_vals)) return(data.frame())
  specs <- joint_parameter_ratio_param_specs()
  rows <- list()
  for (parameter in intersect(active_params, names(specs))) {
    vivo_natural <- suppressWarnings(as.numeric(unname(invivo_vals[[parameter]])))
    vitro_natural <- suppressWarnings(as.numeric(unname(invitro_vals[[parameter]])))
    if (!is.finite(vivo_natural) || !is.finite(vitro_natural) || vitro_natural == 0) next
    transform <- specs[[parameter]]$transform
    vivo_t <- joint_parameter_ratio_apply_transform(vivo_natural, transform)
    vitro_t <- joint_parameter_ratio_apply_transform(vitro_natural, transform)
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = parameter,
      transform = transform,
      vivo_natural = vivo_natural,
      vitro_natural = vitro_natural,
      ratio_vivo_to_vitro = vivo_natural / vitro_natural,
      center_transformed = (vivo_t + vitro_t) / 2,
      delta_transformed = vivo_t - vitro_t,
      source = "separate best_params.tsv",
      invivo_source = attr(invivo_vals, "source_path") %||% NA_character_,
      invitro_source = attr(invitro_vals, "source_path") %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

joint_parameter_ratio_rows_from_warmup_cache <- function(fit_dir, active_params) {
  tab <- joint_parameter_ratio_read_table(file.path(fit_dir, "joint_soft_coupling_parameters_table_input.csv"))
  if (!is.data.frame(tab) || !nrow(tab) || !all(c("param_name", "value") %in% names(tab))) {
    return(data.frame())
  }

  tab$param_name <- as.character(tab$param_name)
  tab$value <- suppressWarnings(as.numeric(tab$value))
  tab$scale <- as.character(tab$scale %||% NA_character_)
  tab <- tab[is.finite(tab$value) & nzchar(tab$param_name), , drop = FALSE]
  if (!nrow(tab)) return(data.frame())

  tab$is_delta <- grepl("^delta__", tab$param_name)
  tab$parameter <- joint_parameter_ratio_natural_name(tab$param_name)
  rows <- list()
  for (parameter in intersect(active_params, unique(tab$parameter))) {
    param_rows <- tab[tab$parameter == parameter, , drop = FALSE]
    center_rows <- param_rows[!param_rows$is_delta, , drop = FALSE]
    if (!nrow(center_rows)) next
    delta_rows <- param_rows[param_rows$is_delta, , drop = FALSE]
    center_row <- center_rows[1L, , drop = FALSE]
    delta_value <- if (nrow(delta_rows)) delta_rows$value[[1]] else 0
    transform <- joint_parameter_ratio_infer_transform(center_row$param_name[[1]], center_row$scale[[1]])
    center_t <- center_row$value[[1]]
    vivo_t <- center_t + delta_value / 2
    vitro_t <- center_t - delta_value / 2
    vivo_natural <- joint_parameter_ratio_inverse_transform(vivo_t, transform)
    vitro_natural <- joint_parameter_ratio_inverse_transform(vitro_t, transform)
    ratio <- if (is.finite(vivo_natural) && is.finite(vitro_natural) && vitro_natural != 0) {
      vivo_natural / vitro_natural
    } else {
      NA_real_
    }
    if (!is.finite(ratio) || ratio <= 0) next
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = parameter,
      transform = transform,
      vivo_natural = vivo_natural,
      vitro_natural = vitro_natural,
      ratio_vivo_to_vitro = ratio,
      center_transformed = center_t,
      delta_transformed = delta_value,
      source = "joint_soft_coupling_parameters_table_input.csv",
      invivo_source = NA_character_,
      invitro_source = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

joint_parameter_ratio_rows_from_soft_coupling <- function(fit_dir, active_params) {
  tab <- joint_parameter_ratio_read_table(file.path(fit_dir, "joint_soft_coupling.tsv"))
  required <- c("parameter", "vivo_natural", "vitro_natural", "ratio_vivo_to_vitro")
  if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab))) {
    return(data.frame())
  }

  tab$parameter <- as.character(tab$parameter)
  tab$vivo_natural <- suppressWarnings(as.numeric(tab$vivo_natural))
  tab$vitro_natural <- suppressWarnings(as.numeric(tab$vitro_natural))
  tab$ratio_vivo_to_vitro <- suppressWarnings(as.numeric(tab$ratio_vivo_to_vitro))
  tab <- tab[
    tab$parameter %in% active_params &
      is.finite(tab$vivo_natural) &
      is.finite(tab$vitro_natural) &
      is.finite(tab$ratio_vivo_to_vitro) &
      tab$ratio_vivo_to_vitro > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(tab)) return(data.frame())

  specs <- joint_parameter_ratio_param_specs()
  rows <- list()
  for (i in seq_len(nrow(tab))) {
    parameter <- tab$parameter[[i]]
    transform <- if (parameter %in% names(specs)) specs[[parameter]]$transform else NA_character_
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = parameter,
      transform = transform,
      vivo_natural = tab$vivo_natural[[i]],
      vitro_natural = tab$vitro_natural[[i]],
      ratio_vivo_to_vitro = tab$ratio_vivo_to_vitro[[i]],
      center_transformed = NA_real_,
      delta_transformed = NA_real_,
      source = "joint_soft_coupling.tsv",
      invivo_source = NA_character_,
      invitro_source = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

joint_parameter_ratio_build_data <- function(fit_dir) {
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  summary_map <- joint_parameter_ratio_read_summary(fit_dir)
  active_params <- joint_parameter_ratio_active_params(summary_map)

  plot_df <- joint_parameter_ratio_rows_from_soft_coupling(fit_dir, active_params)
  if (!nrow(plot_df)) {
    invivo_vals <- joint_parameter_ratio_read_best_values(
      joint_parameter_ratio_seed_candidates(summary_map, "invivo")
    )
    invitro_vals <- joint_parameter_ratio_read_best_values(
      joint_parameter_ratio_seed_candidates(summary_map, "invitro")
    )
    plot_df <- joint_parameter_ratio_rows_from_maps(invivo_vals, invitro_vals, active_params)
  }
  if (!nrow(plot_df)) {
    plot_df <- joint_parameter_ratio_rows_from_warmup_cache(fit_dir, active_params)
  }
  if (!nrow(plot_df)) return(data.frame())

  plot_df <- plot_df[plot_df$parameter %in% active_params, , drop = FALSE]
  mech <- joint_parameter_ratio_mechanism_map()
  plot_df <- merge(plot_df, mech, by = "parameter", all.x = TRUE, sort = FALSE)
  plot_df$mechanism_class[is.na(plot_df$mechanism_class)] <- "Baseline growth"
  plot_df$mechanism_class <- factor(
    plot_df$mechanism_class,
    levels = joint_parameter_ratio_mechanism_levels()
  )
  plot_df <- plot_df[is.finite(plot_df$ratio_vivo_to_vitro) & plot_df$ratio_vivo_to_vitro > 0, , drop = FALSE]
  plot_df[order(plot_df$ratio_vivo_to_vitro), , drop = FALSE]
}
