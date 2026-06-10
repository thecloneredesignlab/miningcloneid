joint_parameter_ratio_null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

`%||%` <- joint_parameter_ratio_null_coalesce

joint_parameter_ratio_soft_specs <- function() {
  list(
    O2_crit = list(transform = "log10", center_name = "log10_O2_crit", delta_name = "delta__log10_O2_crit"),
    alpha_o2 = list(transform = "log10", center_name = "log10_alpha_o2", delta_name = "delta__log10_alpha_o2"),
    mu_hp = list(transform = "log10", center_name = "log10_mu_hp", delta_name = "delta__log10_mu_hp"),
    p_misseg = list(transform = "log10", center_name = "log10_p_misseg", delta_name = "delta__log10_p_misseg"),
    k_o_mis = list(transform = "log10", center_name = "log10_k_o_mis", delta_name = "delta__log10_k_o_mis"),
    buffer_smax = list(transform = "identity", center_name = "buffer_smax", delta_name = "delta__buffer_smax"),
    buffer_beta = list(transform = "log10", center_name = "log10_buffer_beta", delta_name = "delta__log10_buffer_beta"),
    buffer_n_exp = list(transform = "log10", center_name = "log10_buffer_n_exp", delta_name = "delta__log10_buffer_n_exp"),
    n_O = list(transform = "identity", center_name = "n_O", delta_name = "delta__n_O"),
    gamma_growth = list(transform = "identity", center_name = "gamma_growth", delta_name = "delta__gamma_growth")
  )
}

joint_parameter_ratio_nonsoft_specs <- function() {
  list(
    lam_max = list(transform = "log10", joint_init_name = "log10_lam_max"),
    p_mis_base = list(transform = "log10", joint_init_name = "log10_p_mis_base"),
    p_wgd = list(transform = "log10", joint_init_name = "log10_p_wgd"),
    o2_S0 = list(transform = "log10", joint_init_name = "log10_o2_S0"),
    kappa_O = list(transform = "log10", joint_init_name = "log10_kappa_O"),
    eta_o2 = list(transform = "log10", joint_init_name = "log10_eta_o2"),
    rho_2N = list(transform = "log10", joint_init_name = "log10_rho_2N"),
    gamma_mu = list(transform = "identity", joint_init_name = "gamma_mu"),
    qc = list(transform = "identity", joint_init_name = "qc"),
    k_clear = list(transform = "log10", joint_init_name = "log10_k_clear"),
    sigma_burden = list(transform = "log10", joint_init_name = "log10_sigma_burden")
  )
}

joint_parameter_ratio_plot_defaults <- function() {
  c(
    "O2_crit", "lam_max", "buffer_beta", "p_mis_base", "n_O", "p_wgd",
    "buffer_smax", "gamma_mu", "buffer_n_exp", "mu_hp", "k_o_mis", "p_misseg"
  )
}

joint_parameter_ratio_mechanism_map <- function() {
  data.frame(
    parameter = c(
      "p_mis_base", "p_wgd", "k_o_mis", "p_misseg",
      "gamma_mu", "mu_hp",
      "O2_crit", "n_O",
      "lam_max",
      "buffer_beta", "buffer_smax", "buffer_n_exp"
    ),
    mechanism_class = c(
      rep("Baseline growth", 4),
      rep("MS / WGD generation", 2),
      rep("Oxygen-death response", 2),
      "Oxygen-growth response",
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
  specs <- c(joint_parameter_ratio_soft_specs(), joint_parameter_ratio_nonsoft_specs())
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

joint_parameter_ratio_build_data <- function(fit_dir) {
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  summary_map <- joint_parameter_ratio_read_summary(fit_dir)
  active_params <- joint_parameter_ratio_active_params(summary_map)

  invivo_vals <- joint_parameter_ratio_read_best_values(
    joint_parameter_ratio_seed_candidates(summary_map, "invivo")
  )
  invitro_vals <- joint_parameter_ratio_read_best_values(
    joint_parameter_ratio_seed_candidates(summary_map, "invitro")
  )
  plot_df <- joint_parameter_ratio_rows_from_maps(invivo_vals, invitro_vals, active_params)
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

joint_parameter_ratio_axis_label <- function(power) {
  if (power == 0) return("1x")
  if (power < 0) return(paste0("1/", 2^abs(power), "x"))
  paste0(2^power, "x")
}

joint_parameter_ratio_value_label <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return("")
  if (abs(x) >= 0.01 && abs(x) < 1000) return(sprintf("%.2fx", x))
  paste0(formatC(x, format = "e", digits = 2), "x")
}

joint_parameter_ratio_make <- function(plot_df, fit_label = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  if (!is.data.frame(plot_df) || !nrow(plot_df)) return(NULL)

  plot_df <- plot_df[order(plot_df$ratio_vivo_to_vitro), , drop = FALSE]
  plot_df$log2_ratio <- log2(plot_df$ratio_vivo_to_vitro)
  plot_df <- plot_df[is.finite(plot_df$log2_ratio), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)

  n <- nrow(plot_df)
  plot_df$y <- rev(seq_len(n))
  plot_df$segment_x <- pmin(0, plot_df$log2_ratio)
  plot_df$segment_xend <- pmax(0, plot_df$log2_ratio)
  plot_df$ratio_label <- vapply(plot_df$ratio_vivo_to_vitro, joint_parameter_ratio_value_label, character(1))
  plot_df$direction <- ifelse(
    plot_df$ratio_vivo_to_vitro < 1 / 1.10,
    "Higher in vitro",
    ifelse(plot_df$ratio_vivo_to_vitro > 1.10, "Higher in vivo", "Near 1x")
  )
  plot_df$direction <- factor(
    plot_df$direction,
    levels = c("Higher in vitro", "Near 1x", "Higher in vivo")
  )
  plot_df$ratio_label_x <- ifelse(plot_df$log2_ratio < 0, plot_df$log2_ratio - 0.08, plot_df$log2_ratio + 0.08)
  plot_df$ratio_hjust <- ifelse(plot_df$log2_ratio < 0, 1, 0)

  limit_power <- max(6, ceiling(max(abs(plot_df$log2_ratio), na.rm = TRUE)))
  limit_power <- min(max(limit_power, 1), 20)
  x_breaks <- seq(-limit_power, limit_power, by = 1)
  x_limits <- c(-limit_power - 0.75, limit_power + 0.75)
  parameter_label_x <- -limit_power + 0.10

  subtitle <- "Point/segment color = direction of difference | Label color = mechanism class"
  if (!is.null(fit_label) && length(fit_label) && nzchar(as.character(fit_label[[1]]))) {
    subtitle <- paste0(as.character(fit_label[[1]]), " | ", subtitle)
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(y = y)) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linetype = "22", linewidth = 0.7) +
    ggplot2::geom_segment(
      ggplot2::aes(x = segment_x, xend = segment_xend, yend = y, color = direction),
      linewidth = 1.35,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = log2_ratio, color = direction),
      size = 4.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = ratio_label_x,
        label = ratio_label,
        color = direction,
        hjust = ratio_hjust
      ),
      size = 4.6,
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      ggplot2::aes(
        x = parameter_label_x,
        label = parameter,
        fill = mechanism_class
      ),
      color = "black",
      fontface = "bold",
      size = 4.9,
      linewidth = 0.25,
      label.padding = grid::unit(0.16, "lines"),
      label.r = grid::unit(0.11, "lines"),
      show.legend = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = vapply(x_breaks, joint_parameter_ratio_axis_label, character(1)),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.55, n + 0.75),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_color_manual(
      name = "Direction",
      values = joint_parameter_ratio_direction_colors(),
      breaks = c("Higher in vitro", "Higher in vivo"),
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = joint_parameter_ratio_mechanism_colors(),
      breaks = joint_parameter_ratio_mechanism_levels(),
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(label = "a", color = "black", size = 4.4)
      ),
      color = ggplot2::guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(size = 4.5, linewidth = 1.5)
      )
    ) +
    ggplot2::labs(
      title = "In Vivo vs In Vitro Parameter Ratios",
      subtitle = subtitle,
      x = "ratio_vivo_to_vitro on a symmetric log2 scale",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.7),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "#555555", size = 12),
      axis.title.x = ggplot2::element_text(size = 15),
      plot.title = ggplot2::element_text(size = 23, face = "plain", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0, color = "black"),
      legend.position = "top",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.spacing.x = grid::unit(0.30, "cm"),
      legend.margin = ggplot2::margin(0, 0, 6, 0),
      plot.margin = ggplot2::margin(10, 22, 10, 10)
    )
}

plot_joint_parameter_ratio_figure <- function(
  fit_dir,
  out_dir = file.path(fit_dir, "viz", "joint_parameters"),
  basename = joint_parameter_ratio_output_basename(fit_dir),
  width = 18,
  height = 6.5,
  dpi = 180
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  plot_df <- joint_parameter_ratio_build_data(fit_dir)
  if (!nrow(plot_df)) return(NULL)
  summary_map <- joint_parameter_ratio_read_summary(fit_dir)
  p <- joint_parameter_ratio_make(
    plot_df,
    fit_label = joint_parameter_ratio_fit_report_label(fit_dir, summary_map)
  )
  if (is.null(p)) return(NULL)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(basename, ".pdf"))
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  ggplot2::ggsave(pdf_path, p, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_path, p, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  list(
    pdf = normalizePath(pdf_path, mustWork = TRUE),
    png = normalizePath(png_path, mustWork = TRUE),
    data = plot_df
  )
}
