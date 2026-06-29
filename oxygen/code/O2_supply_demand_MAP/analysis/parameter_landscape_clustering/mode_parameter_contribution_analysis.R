#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

source(file.path(local_script_dir(), "parameter_landscape_utils.R"))

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_+.-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- sub("^_", "", sub("_$", "", x))
  ifelse(nzchar(x), x, "value")
}

parameter_axis_label <- function(parameter, log10_params) {
  if (parameter %in% log10_params) paste0("log10(", parameter, ")") else parameter
}

format_o2_value <- function(o2) {
  vals <- suppressWarnings(as.numeric(o2))
  vapply(vals, function(x) {
    if (!is.finite(x)) return(NA_character_)
    txt <- sprintf("%.10f", x)
    txt <- sub("0+$", "", txt)
    txt <- sub("\\.$", "", txt)
    if (identical(txt, "-0")) txt <- "0"
    txt
  }, character(1))
}

reference_o2_plot_title <- function(prefix, o2) {
  bquote(.(prefix) * O[2] == .(format_o2_value(o2)[[1L]]))
}

normalize_mode_contribution_target <- function(target) {
  target <- tolower(trimws(as.character(target %||% "mode")))
  if (target %in% c("mode", "discrete", "classification", "logistic")) return("mode")
  if (target %in% c(
    "dominant_ploidy", "dominant_mean_ploidy", "mean_ploidy",
    "ploidy", "continuous", "regression"
  )) {
    return("dominant_ploidy")
  }
  stop("mode_contribution_target must be mode or dominant_ploidy.", call. = FALSE)
}

mode_contribution_default_output_dir <- function(root_dir, target) {
  target <- normalize_mode_contribution_target(target)
  if (identical(target, "dominant_ploidy")) {
    return(file.path(root_dir, "dominant_ploidy_parameter_contribution"))
  }
  file.path(root_dir, "mode_parameter_contribution")
}

mode_contribution_target_label <- function(target) {
  target <- normalize_mode_contribution_target(target)
  if (identical(target, "dominant_ploidy")) "dominant mean ploidy" else "Mode"
}

o2_discrete_label_expressions <- function(x) {
  vals <- sub("^O2=", "", as.character(x))
  vals[!nzchar(vals) | is.na(vals)] <- "NA"
  parse(text = paste0("O[2] == ", vals))
}

GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin",
  ggplot2::GeomViolin,
  draw_group = function(self, data, panel_params, coord, ..., draw_quantiles = NULL) {
    data <- data[order(data$y), , drop = FALSE]
    data$xminv <- data$x - data$violinwidth * (data$x - data$xmin)
    data$xmaxv <- data$x + data$violinwidth * (data$xmax - data$x)
    if (data$group[[1L]] %% 2L == 1L) {
      newdata <- transform(data, x = xminv)
      newdata <- newdata[order(newdata$y), , drop = FALSE]
    } else {
      newdata <- transform(data, x = xmaxv)
      newdata <- newdata[order(-newdata$y), , drop = FALSE]
    }
    newdata <- rbind(newdata[1L, , drop = FALSE], newdata, newdata[nrow(newdata), , drop = FALSE], newdata[1L, , drop = FALSE])
    newdata$x[c(1L, nrow(newdata) - 1L, nrow(newdata))] <- data$x[[1L]]
    ggplot2::GeomPolygon$draw_panel(newdata, panel_params, coord, ...)
  }
)

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              ...,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, na.rm = na.rm, ...)
  )
}

numeric_column_or_na <- function(df, column) {
  if (!column %in% names(df)) return(rep(NA_real_, nrow(df)))
  suppressWarnings(as.numeric(df[[column]]))
}

signed_auc <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- as.integer(y[keep])
  score <- as.numeric(score[keep])
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

safe_wilcox_p <- function(y, x) {
  out <- tryCatch(
    stats::wilcox.test(x[y == 1L], x[y == 0L], exact = FALSE)$p.value,
    error = function(e) NA_real_,
    warning = function(w) suppressWarnings(stats::wilcox.test(x[y == 1L], x[y == 0L], exact = FALSE)$p.value)
  )
  if (is.finite(out)) out else NA_real_
}

safe_glm_term <- function(df, formula, term, predicted_auc = FALSE) {
  fit <- tryCatch(
    suppressWarnings(stats::glm(formula, data = df, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(coef = NA_real_, z = NA_real_, p = NA_real_, auc = NA_real_))
  }
  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  coef <- z <- p <- NA_real_
  if (!is.null(co) && term %in% rownames(co)) {
    coef <- suppressWarnings(as.numeric(co[term, "Estimate"]))
    z <- suppressWarnings(as.numeric(co[term, "z value"]))
    p <- suppressWarnings(as.numeric(co[term, "Pr(>|z|)"]))
  }
  auc <- NA_real_
  if (isTRUE(predicted_auc)) {
    pred <- tryCatch(suppressWarnings(stats::predict(fit, type = "response")), error = function(e) NULL)
    if (!is.null(pred)) auc <- signed_auc(df$y, pred)
  }
  list(coef = coef, z = z, p = p, auc = auc)
}

standardize_matrix <- function(mat) {
  z <- scale(as.matrix(mat))
  scales <- attr(z, "scaled:scale")
  zero <- !is.finite(scales) | scales == 0
  if (any(zero)) {
    stop("Zero/non-finite feature SD for: ", paste(colnames(mat)[zero], collapse = ", "))
  }
  z
}

build_interaction_design <- function(z_main) {
  params <- colnames(z_main)
  pairs <- utils::combn(params, 2L, simplify = FALSE)
  inter <- vapply(
    pairs,
    function(pair) z_main[, pair[[1L]]] * z_main[, pair[[2L]]],
    numeric(nrow(z_main))
  )
  if (is.null(dim(inter))) inter <- matrix(inter, ncol = 1L)
  colnames(inter) <- vapply(pairs, function(pair) paste(pair, collapse = ":"), character(1))
  list(
    interaction_matrix = inter,
    interaction_pairs = data.frame(
      feature_name = colnames(inter),
      parameter_a = vapply(pairs, `[[`, character(1), 1L),
      parameter_b = vapply(pairs, `[[`, character(1), 2L),
      stringsAsFactors = FALSE
    )
  )
}

main_effect_table <- function(best_df, transformed_df, z_main, y, params, log10_params) {
  rows <- lapply(params, function(param) {
    x <- z_main[, param]
    xt <- transformed_df[[param]]
    xo <- suppressWarnings(as.numeric(best_df[[param]]))
    mode1 <- y == 1L
    mode2 <- y == 0L
    pooled_sd <- stats::sd(x, na.rm = TRUE)
    delta_z <- mean(x[mode1], na.rm = TRUE) - mean(x[mode2], na.rm = TRUE)
    auc <- signed_auc(y, x)
    glm <- safe_glm_term(data.frame(y = y, x = x), y ~ x, "x", predicted_auc = TRUE)
    data.frame(
      feature_name = param,
      feature_type = "main",
      parameter = param,
      parameter_a = param,
      parameter_b = NA_character_,
      transformed_label = parameter_axis_label(param, log10_params),
      n_seed = length(y),
      n_mode1 = sum(mode1),
      n_mode2 = sum(mode2),
      mean_mode1_original = mean(xo[mode1], na.rm = TRUE),
      mean_mode2_original = mean(xo[mode2], na.rm = TRUE),
      delta_mode1_minus_mode2_original = mean(xo[mode1], na.rm = TRUE) - mean(xo[mode2], na.rm = TRUE),
      mean_mode1_transformed = mean(xt[mode1], na.rm = TRUE),
      mean_mode2_transformed = mean(xt[mode2], na.rm = TRUE),
      delta_mode1_minus_mode2_transformed = mean(xt[mode1], na.rm = TRUE) - mean(xt[mode2], na.rm = TRUE),
      delta_mode1_minus_mode2_z = delta_z,
      cohen_d_z = if (is.finite(pooled_sd) && pooled_sd > 0) delta_z / pooled_sd else NA_real_,
      auc_mode1_higher = auc,
      auc_contribution = if (is.finite(auc)) abs(auc - 0.5) * 2 else NA_real_,
      wilcox_p = safe_wilcox_p(y, x),
      glm_coef_mode1 = glm$coef,
      glm_z = glm$z,
      glm_p = glm$p,
      glm_predicted_auc = glm$auc,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$wilcox_fdr <- stats::p.adjust(out$wilcox_p, method = "BH")
  out$glm_fdr <- stats::p.adjust(out$glm_p, method = "BH")
  out[order(-out$auc_contribution, out$wilcox_p, out$parameter), , drop = FALSE]
}

interaction_table <- function(z_main, y) {
  design <- build_interaction_design(z_main)
  pairs <- design$interaction_pairs
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pa <- pairs$parameter_a[[i]]
    pb <- pairs$parameter_b[[i]]
    xi <- design$interaction_matrix[, i]
    df <- data.frame(
      y = y,
      xa = z_main[, pa],
      xb = z_main[, pb],
      xi = xi
    )
    glm <- safe_glm_term(df, y ~ xa + xb + xi, "xi", predicted_auc = TRUE)
    auc <- signed_auc(y, xi)
    data.frame(
      feature_name = pairs$feature_name[[i]],
      feature_type = "interaction",
      parameter = NA_character_,
      parameter_a = pa,
      parameter_b = pb,
      n_seed = length(y),
      n_mode1 = sum(y == 1L),
      n_mode2 = sum(y == 0L),
      interaction_mean_mode1 = mean(xi[y == 1L], na.rm = TRUE),
      interaction_mean_mode2 = mean(xi[y == 0L], na.rm = TRUE),
      interaction_delta_mode1_minus_mode2 = mean(xi[y == 1L], na.rm = TRUE) - mean(xi[y == 0L], na.rm = TRUE),
      interaction_auc_mode1_higher = auc,
      interaction_auc_contribution = if (is.finite(auc)) abs(auc - 0.5) * 2 else NA_real_,
      interaction_glm_coef_mode1 = glm$coef,
      interaction_glm_z = glm$z,
      interaction_glm_p = glm$p,
      interaction_glm_fitted_auc = glm$auc,
      interaction_model_auc_contribution = if (is.finite(glm$auc)) abs(glm$auc - 0.5) * 2 else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$interaction_glm_fdr <- stats::p.adjust(out$interaction_glm_p, method = "BH")
  out[order(out$interaction_glm_p, -out$interaction_model_auc_contribution, out$feature_name), , drop = FALSE]
}

discriminative_auc <- function(auc) {
  auc <- suppressWarnings(as.numeric(auc))
  ifelse(is.finite(auc), pmax(auc, 1 - auc), NA_real_)
}

auc_direction <- function(auc) {
  auc <- suppressWarnings(as.numeric(auc))
  ifelse(
    is.finite(auc),
    ifelse(auc >= 0.5, "higher_in_mode1", "higher_in_mode2"),
    NA_character_
  )
}

feature_vector_from_row <- function(feature_row, z_main) {
  feature_type <- as.character(feature_row$feature_type[[1L]])
  feature_name <- as.character(feature_row$feature_name[[1L]])
  if (identical(feature_type, "main")) {
    param <- as.character(feature_row$parameter_a[[1L]] %||% feature_row$parameter[[1L]] %||% feature_name)
    if (!param %in% colnames(z_main)) stop("Unknown main feature parameter: ", param)
    return(as.numeric(z_main[, param]))
  }
  pa <- as.character(feature_row$parameter_a[[1L]])
  pb <- as.character(feature_row$parameter_b[[1L]])
  if (!pa %in% colnames(z_main) || !pb %in% colnames(z_main)) {
    stop("Unknown interaction feature parameters: ", feature_name)
  }
  as.numeric(z_main[, pa]) * as.numeric(z_main[, pb])
}

feature_auc_table <- function(combined, z_main, y) {
  if (!nrow(combined)) return(data.frame())
  rows <- lapply(seq_len(nrow(combined)), function(i) {
    row <- combined[i, , drop = FALSE]
    x <- feature_vector_from_row(row, z_main)
    auc <- signed_auc(y, x)
    data.frame(
      feature_name = row$feature_name,
      feature_type = row$feature_type,
      parameter = row$parameter %||% NA_character_,
      parameter_a = row$parameter_a %||% NA_character_,
      parameter_b = row$parameter_b %||% NA_character_,
      rank = suppressWarnings(as.integer(row$rank)),
      auc_mode1_higher = auc,
      auc_discriminative = discriminative_auc(auc),
      direction = auc_direction(auc),
      elastic_net_coef = suppressWarnings(as.numeric(row$elastic_net_coef %||% NA_real_)),
      selection_frequency = suppressWarnings(as.numeric(row$selection_frequency %||% NA_real_)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$auc_discriminative, out$rank, out$feature_name), , drop = FALSE]
  out$auc_rank <- seq_len(nrow(out))
  out
}

plot_feature_auc_bar <- function(auc_df,
                                 output_prefix,
                                 title,
                                 top_n = NULL,
                                 width = 7.8,
                                 height = NULL,
                                 flip = TRUE) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(character()))
  df <- auc_df[is.finite(auc_df$auc_discriminative), , drop = FALSE]
  df <- df[order(-df$auc_discriminative, df$feature_name), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n) && top_n > 0L) {
    df <- head(df, as.integer(top_n))
  }
  if (!nrow(df)) return(invisible(character()))
  df$feature_index <- seq_len(nrow(df))
  if (is.null(height)) {
    height <- max(4.2, min(30, 2.2 + 0.16 * nrow(df)))
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(fill = direction)) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = feature_index - 0.36,
        xmax = feature_index + 0.36,
        ymin = 0.5,
        ymax = auc_discriminative
      )
    )
  if (isTRUE(flip)) {
    p <- p + ggplot2::coord_flip()
  }
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = df$feature_index,
      labels = df$feature_name,
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
    ggplot2::scale_fill_manual(
      values = c(higher_in_mode1 = "#1f78b4", higher_in_mode2 = "#e31a1c"),
      na.value = "grey70",
      labels = c(higher_in_mode1 = "higher in Mode1", higher_in_mode2 = "higher in Mode2")
    ) +
    ggplot2::labs(
      title = title,
      x = if (isTRUE(flip)) NULL else "Feature",
      y = "AUC",
      fill = "Direction"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = if (isTRUE(flip)) {
        ggplot2::element_text()
      } else {
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
      },
      legend.position = if (isTRUE(flip)) "right" else "bottom"
    )
  save_plot_pair(p, output_prefix, width = width, height = height)
  invisible(paste0(output_prefix, c(".pdf", ".png")))
}

roc_curve_df <- function(y, score, curve = "score") {
  keep <- is.finite(y) & is.finite(score)
  y <- as.integer(y[keep])
  score <- as.numeric(score[keep])
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(data.frame())
  ord <- order(score, decreasing = TRUE)
  yy <- y[ord]
  tp <- cumsum(yy == 1L)
  fp <- cumsum(yy == 0L)
  data.frame(
    fpr = c(0, fp / n0),
    tpr = c(0, tp / n1),
    curve = curve,
    stringsAsFactors = FALSE
  )
}

stratified_fold_ids <- function(y, k = 10L, seed = 123L) {
  y <- as.integer(y)
  k <- max(2L, min(as.integer(k), min(table(y))))
  folds <- integer(length(y))
  set.seed(seed)
  for (cls in sort(unique(y))) {
    idx <- sample(which(y == cls))
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

fit_top3_joint_auc <- function(combined,
                               z_main,
                               y,
                               mode_reference_o2,
                               output_dir,
                               cv_k = 10L,
                               seed = 123L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  top3 <- head(combined[order(combined$rank), , drop = FALSE], 3L)
  if (nrow(top3) < 1L) return(NULL)
  x <- do.call(cbind, lapply(seq_len(nrow(top3)), function(i) feature_vector_from_row(top3[i, , drop = FALSE], z_main)))
  colnames(x) <- paste0("top", seq_len(ncol(x)))
  df <- data.frame(y = as.integer(y), x, check.names = FALSE)
  fit <- tryCatch(
    suppressWarnings(stats::glm(stats::reformulate(colnames(x), response = "y"), data = df, family = stats::binomial())),
    error = function(e) e
  )
  apparent_pred <- rep(NA_real_, length(y))
  apparent_auc <- NA_real_
  fit_status <- "ok"
  fit_reason <- ""
  if (inherits(fit, "error")) {
    fit_status <- "failed"
    fit_reason <- conditionMessage(fit)
  } else {
    apparent_pred <- suppressWarnings(stats::predict(fit, type = "response"))
    apparent_auc <- signed_auc(y, apparent_pred)
  }

  folds <- stratified_fold_ids(y, k = cv_k, seed = seed)
  cv_pred <- rep(NA_real_, length(y))
  cv_success <- 0L
  for (fold in sort(unique(folds))) {
    train <- folds != fold
    test <- folds == fold
    fold_fit <- tryCatch(
      suppressWarnings(stats::glm(stats::reformulate(colnames(x), response = "y"), data = df[train, , drop = FALSE], family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(fold_fit)) next
    cv_pred[test] <- suppressWarnings(stats::predict(fold_fit, newdata = df[test, , drop = FALSE], type = "response"))
    cv_success <- cv_success + 1L
  }
  cv_auc <- signed_auc(y, cv_pred)

  pred_tab <- data.frame(
    row_id = seq_along(y),
    mode_reference_o2 = mode_reference_o2,
    mode_label = ifelse(y == 1L, "mode1", "mode2"),
    apparent_probability_mode1 = apparent_pred,
    cv_probability_mode1 = cv_pred,
    fold_id = folds,
    stringsAsFactors = FALSE
  )
  write_csv(pred_tab, file.path(output_dir, "top3_joint_predictions.csv"))

  roc <- rbind_fill_plain(list(
    roc_curve_df(y, apparent_pred, "apparent"),
    roc_curve_df(y, cv_pred, "cross_validated")
  ))
  if (nrow(roc)) {
    write_csv(roc, file.path(output_dir, "top3_joint_roc_points.csv"))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      auc_label <- function(x) {
        if (!is.finite(x)) return("NA")
        if (abs(x - round(x)) < 1e-9) return(format(round(x), scientific = FALSE, trim = TRUE))
        formatC(x, digits = 3L, format = "fg")
      }
      curve_labels <- c(
        apparent = paste0("apparent AUC=", auc_label(apparent_auc)),
        cross_validated = paste0("cross-validated AUC=", auc_label(cv_auc))
      )
      p <- ggplot2::ggplot(roc, ggplot2::aes(x = fpr, y = tpr, color = curve)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
        ggplot2::geom_path(linewidth = 0.9) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        ggplot2::scale_color_manual(
          values = c(apparent = "#1f78b4", cross_validated = "#e31a1c"),
          breaks = names(curve_labels),
          labels = unname(curve_labels)
        ) +
        ggplot2::labs(
          title = reference_o2_plot_title("Top3 joint ROC, reference ", mode_reference_o2),
          x = "False positive rate",
          y = "True positive rate",
          color = "Evaluation"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
      save_plot_pair(p, file.path(output_dir, "top3_joint_roc"), width = 5.2, height = 4.6)
    }
  }

  top_names <- top3$feature_name
  top_names <- c(top_names, rep(NA_character_, 3L - length(top_names)))
  summary <- data.frame(
    mode_reference_o2 = mode_reference_o2,
    n_seed = length(y),
    n_mode1 = sum(y == 1L),
    n_mode2 = sum(y == 0L),
    top1_feature = top_names[[1L]],
    top2_feature = top_names[[2L]],
    top3_feature = top_names[[3L]],
    fit_status = fit_status,
    fit_reason = fit_reason,
    apparent_auc_mode1_higher = apparent_auc,
    apparent_auc = discriminative_auc(apparent_auc),
    cv_auc_mode1_higher = cv_auc,
    cv_auc = discriminative_auc(cv_auc),
    cv_k = length(unique(folds)),
    cv_success_folds = cv_success,
    stringsAsFactors = FALSE
  )
  write_csv(summary, file.path(output_dir, "top3_joint_auc.csv"))
  summary
}

plot_top_feature_rocs <- function(combined,
                                  z_main,
                                  y,
                                  mode_reference_o2,
                                  output_dir,
                                  n = 3L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(character()))
  if (!nrow(combined)) return(invisible(character()))
  top <- head(combined[order(combined$rank), , drop = FALSE], as.integer(n))
  if (!nrow(top)) return(invisible(character()))

  auc_label <- function(x) {
    if (!is.finite(x)) return("NA")
    if (abs(x - round(x)) < 1e-9) return(format(round(x), scientific = FALSE, trim = TRUE))
    formatC(x, digits = 3L, format = "fg")
  }
  direction_label <- function(direction) {
    if (identical(direction, "higher_in_mode1")) return("Higher values in Mode1")
    if (identical(direction, "higher_in_mode2")) return("Higher values in Mode2")
    "Direction unavailable"
  }

  roc_rows <- list()
  written <- character()
  for (i in seq_len(nrow(top))) {
    row <- top[i, , drop = FALSE]
    score <- feature_vector_from_row(row, z_main)
    signed <- signed_auc(y, score)
    discriminative <- discriminative_auc(signed)
    direction <- auc_direction(signed)
    roc_score <- if (is.finite(signed) && signed < 0.5) -score else score
    roc <- roc_curve_df(y, roc_score, paste0("top", i))
    if (!nrow(roc)) next
    roc$mode_reference_o2 <- mode_reference_o2
    roc$feature_rank <- i
    roc$feature_name <- row$feature_name[[1L]]
    roc$feature_type <- row$feature_type[[1L]]
    roc$auc_mode1_higher <- signed
    roc$auc_discriminative <- discriminative
    roc$direction <- direction
    roc_rows[[length(roc_rows) + 1L]] <- roc

    p <- ggplot2::ggplot(roc, ggplot2::aes(x = fpr, y = tpr)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
      ggplot2::geom_path(color = "#1f78b4", linewidth = 0.95) +
      ggplot2::annotate(
        "label",
        x = 0.54,
        y = 0.16,
        hjust = 0,
        label = paste0("AUC = ", auc_label(discriminative), "\n", direction_label(direction)),
        size = 3.1,
        label.size = 0.2
      ) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::labs(
        title = reference_o2_plot_title(paste0("Top", i, " feature ROC, reference "), mode_reference_o2),
        subtitle = row$feature_name[[1L]],
        x = "False positive rate",
        y = "True positive rate"
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_text(size = 8.5)
      )
    prefix <- file.path(output_dir, paste0("top", i, "_feature_roc"))
    save_plot_pair(p, prefix, width = 4.2, height = 4.2)
    written <- c(written, paste0(prefix, c(".pdf", ".png")))
  }
  if (length(roc_rows)) {
    write_csv(rbind_fill_plain(roc_rows), file.path(output_dir, "top3_single_feature_roc_points.csv"))
  }
  invisible(written)
}

cumulative_feature_matrix <- function(features, z_main) {
  x <- do.call(cbind, lapply(seq_len(nrow(features)), function(i) feature_vector_from_row(features[i, , drop = FALSE], z_main)))
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  colnames(x) <- paste0("feature_", seq_len(ncol(x)))
  x
}

fit_logistic_auc_for_matrix <- function(x, y, folds) {
  x <- as.matrix(x)
  colnames(x) <- paste0("feature_", seq_len(ncol(x)))
  df <- data.frame(y = as.integer(y), x, check.names = FALSE)
  formula <- stats::reformulate(colnames(x), response = "y")
  fit <- tryCatch(
    suppressWarnings(stats::glm(formula, data = df, family = stats::binomial())),
    error = function(e) e
  )
  apparent_pred <- rep(NA_real_, length(y))
  fit_status <- "ok"
  fit_reason <- ""
  if (inherits(fit, "error")) {
    fit_status <- "failed"
    fit_reason <- conditionMessage(fit)
  } else {
    apparent_pred <- suppressWarnings(stats::predict(fit, type = "response"))
  }

  cv_pred <- rep(NA_real_, length(y))
  cv_success <- 0L
  for (fold in sort(unique(folds))) {
    train <- folds != fold
    test <- folds == fold
    fold_fit <- tryCatch(
      suppressWarnings(stats::glm(formula, data = df[train, , drop = FALSE], family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(fold_fit)) next
    cv_pred[test] <- suppressWarnings(stats::predict(fold_fit, newdata = df[test, , drop = FALSE], type = "response"))
    cv_success <- cv_success + 1L
  }
  list(
    apparent_auc_mode1_higher = signed_auc(y, apparent_pred),
    apparent_auc = discriminative_auc(signed_auc(y, apparent_pred)),
    cv_auc_mode1_higher = signed_auc(y, cv_pred),
    cv_auc = discriminative_auc(signed_auc(y, cv_pred)),
    cv_success_folds = cv_success,
    fit_status = fit_status,
    fit_reason = fit_reason
  )
}

plot_cumulative_feature_auc <- function(cumulative_tab, output_prefix, title) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE) || !nrow(cumulative_tab)) return(invisible(character()))
  plot_df <- rbind_fill_plain(list(
    data.frame(
      k = cumulative_tab$k,
      auc = cumulative_tab$cv_auc,
      evaluation = "cross_validated",
      evaluation_label = "Cross-validated",
      stringsAsFactors = FALSE
    ),
    data.frame(
      k = cumulative_tab$k,
      auc = cumulative_tab$apparent_auc,
      evaluation = "apparent",
      evaluation_label = "Apparent",
      stringsAsFactors = FALSE
    )
  ))
  plot_df <- plot_df[is.finite(plot_df$k) & is.finite(plot_df$auc), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(character()))
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = k, y = auc, color = evaluation_label, linetype = evaluation_label)
  ) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_vline(xintercept = c(3, 5, 10), linetype = "dotted", color = "grey60") +
    ggplot2::scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
    ggplot2::scale_color_manual(values = c("Cross-validated" = "#1f78b4", "Apparent" = "grey45")) +
    ggplot2::scale_linetype_manual(values = c("Cross-validated" = "solid", "Apparent" = "dashed")) +
    ggplot2::labs(
      title = title,
      x = "Cumulative feature count",
      y = "AUC",
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  save_plot_pair(p, output_prefix, width = 6.4, height = 4.6)
  invisible(paste0(output_prefix, c(".pdf", ".png")))
}

fit_cumulative_feature_auc <- function(combined,
                                       z_main,
                                       y,
                                       mode_reference_o2,
                                       output_dir,
                                       top_n = 30L,
                                       cv_k = 10L,
                                       seed = 123L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  top_n <- max(1L, min(as.integer(top_n), nrow(combined)))
  ranked <- combined[order(combined$rank), , drop = FALSE]
  ranked <- ranked[seq_len(top_n), , drop = FALSE]
  folds <- stratified_fold_ids(y, k = cv_k, seed = seed)
  rows <- lapply(seq_len(nrow(ranked)), function(k) {
    selected <- ranked[seq_len(k), , drop = FALSE]
    x <- cumulative_feature_matrix(selected, z_main)
    auc <- fit_logistic_auc_for_matrix(x, y, folds)
    data.frame(
      mode_reference_o2 = mode_reference_o2,
      k = k,
      added_feature = selected$feature_name[[k]],
      added_feature_type = selected$feature_type[[k]],
      feature_set = paste(selected$feature_name, collapse = ";"),
      n_seed = length(y),
      n_mode1 = sum(y == 1L),
      n_mode2 = sum(y == 0L),
      cv_k = length(unique(folds)),
      cv_success_folds = auc$cv_success_folds,
      apparent_auc_mode1_higher = auc$apparent_auc_mode1_higher,
      apparent_auc = auc$apparent_auc,
      cv_auc_mode1_higher = auc$cv_auc_mode1_higher,
      cv_auc = auc$cv_auc,
      fit_status = auc$fit_status,
      fit_reason = auc$fit_reason,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$apparent_auc_gain <- c(NA_real_, diff(out$apparent_auc))
  out$cv_auc_gain <- c(NA_real_, diff(out$cv_auc))
  write_csv(out, file.path(output_dir, "cumulative_feature_auc.csv"))
  plot_cumulative_feature_auc(
    out,
    file.path(output_dir, "cumulative_feature_auc"),
    reference_o2_plot_title("Cumulative feature AUC, reference ", mode_reference_o2)
  )
  out
}

write_auc_outputs <- function(combined,
                              z_main,
                              y,
                              mode_reference_o2,
                              ref_dir) {
  auc_dir <- file.path(ref_dir, "auc_plots")
  dir.create(auc_dir, recursive = TRUE, showWarnings = FALSE)
  auc_tab <- feature_auc_table(combined, z_main, y)
  write_csv(auc_tab, file.path(auc_dir, "mode_parameter_auc_by_feature.csv"))
  main_auc <- auc_tab[auc_tab$feature_type == "main", , drop = FALSE]
  if (nrow(main_auc)) {
    plot_feature_auc_bar(
      main_auc,
      file.path(auc_dir, "main_feature_auc_bar"),
      reference_o2_plot_title("Main feature AUC, reference ", mode_reference_o2)
    )
  }
  if (nrow(auc_tab)) {
    plot_feature_auc_bar(
      auc_tab,
      file.path(auc_dir, "all_feature_auc_bar"),
      reference_o2_plot_title("All feature AUC, reference ", mode_reference_o2)
    )
    plot_feature_auc_bar(
      auc_tab,
      file.path(auc_dir, "combined_feature_auc_bar_top30"),
      reference_o2_plot_title("Top30 retained-feature AUC, reference ", mode_reference_o2),
      top_n = 30L,
      width = 10.5,
      height = 5.25,
      flip = FALSE
    )
    plot_top_feature_rocs(
      combined,
      z_main,
      y,
      mode_reference_o2,
      file.path(ref_dir, "top3_feature_roc"),
      n = 3L
    )
  }
  auc_tab
}

top_feature_summary_fields <- function(combined, auc_tab, n = 3L) {
  top <- head(combined[order(combined$rank), , drop = FALSE], n)
  out <- data.frame(.row = 1L)
  for (i in seq_len(n)) {
    prefix <- paste0("top", i)
    if (i <= nrow(top)) {
      feature <- top$feature_name[[i]]
      auc_row <- auc_tab[auc_tab$feature_name == feature, , drop = FALSE]
      out[[paste0(prefix, "_feature")]] <- feature
      out[[paste0(prefix, "_feature_type")]] <- top$feature_type[[i]]
      out[[paste0(prefix, "_auc")]] <- if (nrow(auc_row)) auc_row$auc_discriminative[[1L]] else NA_real_
      out[[paste0(prefix, "_auc_mode1_higher")]] <- if (nrow(auc_row)) auc_row$auc_mode1_higher[[1L]] else NA_real_
      out[[paste0(prefix, "_direction")]] <- if (nrow(auc_row)) auc_row$direction[[1L]] else NA_character_
      out[[paste0(prefix, "_elastic_net_coef")]] <- suppressWarnings(as.numeric(top$elastic_net_coef[[i]]))
      out[[paste0(prefix, "_selection_frequency")]] <- suppressWarnings(as.numeric(top$selection_frequency[[i]]))
    } else {
      out[[paste0(prefix, "_feature")]] <- NA_character_
      out[[paste0(prefix, "_feature_type")]] <- NA_character_
      out[[paste0(prefix, "_auc")]] <- NA_real_
      out[[paste0(prefix, "_auc_mode1_higher")]] <- NA_real_
      out[[paste0(prefix, "_direction")]] <- NA_character_
      out[[paste0(prefix, "_elastic_net_coef")]] <- NA_real_
      out[[paste0(prefix, "_selection_frequency")]] <- NA_real_
    }
  }
  out$.row <- NULL
  out
}

stratified_subsample <- function(y, fraction) {
  idx1 <- which(y == 1L)
  idx0 <- which(y == 0L)
  n1 <- max(1L, floor(length(idx1) * fraction))
  n0 <- max(1L, floor(length(idx0) * fraction))
  c(sample(idx1, n1), sample(idx0, n0))
}

fit_glmnet_importance <- function(z_main,
                                  y,
                                  n_bootstrap = 100L,
                                  sample_fraction = 0.75,
                                  alpha = 0.5,
                                  seed = 123L) {
  inter <- build_interaction_design(z_main)
  x <- cbind(z_main, inter$interaction_matrix)
  feature_names <- colnames(x)
  info <- data.frame(
    feature_name = feature_names,
    feature_type = ifelse(feature_names %in% colnames(z_main), "main", "interaction"),
    parameter = ifelse(feature_names %in% colnames(z_main), feature_names, NA_character_),
    parameter_a = NA_character_,
    parameter_b = NA_character_,
    stringsAsFactors = FALSE
  )
  info$parameter_a[info$feature_type == "main"] <- info$feature_name[info$feature_type == "main"]
  info$parameter_b[info$feature_type == "main"] <- NA_character_
  idx <- match(info$feature_name, inter$interaction_pairs$feature_name)
  hit <- !is.na(idx)
  info$parameter_a[hit] <- inter$interaction_pairs$parameter_a[idx[hit]]
  info$parameter_b[hit] <- inter$interaction_pairs$parameter_b[idx[hit]]

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    info$elastic_net_coef <- NA_real_
    info$selection_frequency <- NA_real_
    attr(info, "status") <- data.frame(
      model = "elastic_net_logistic",
      status = "unavailable",
      reason = "R package glmnet is not installed",
      stringsAsFactors = FALSE
    )
    return(info)
  }

  min_class <- min(table(y))
  nfolds <- max(3L, min(10L, as.integer(min_class)))
  if (nfolds < 3L) {
    info$elastic_net_coef <- NA_real_
    info$selection_frequency <- NA_real_
    attr(info, "status") <- data.frame(
      model = "elastic_net_logistic",
      status = "unavailable",
      reason = "fewer than three seeds in one mode",
      stringsAsFactors = FALSE
    )
    return(info)
  }

  set.seed(seed)
  cv <- tryCatch(
    glmnet::cv.glmnet(x, y, family = "binomial", alpha = alpha, standardize = FALSE, nfolds = nfolds),
    error = function(e) e
  )
  if (inherits(cv, "error")) {
    info$elastic_net_coef <- NA_real_
    info$selection_frequency <- NA_real_
    attr(info, "status") <- data.frame(
      model = "elastic_net_logistic",
      status = "failed",
      reason = conditionMessage(cv),
      stringsAsFactors = FALSE
    )
    return(info)
  }

  coef_mat <- as.matrix(stats::coef(cv, s = "lambda.1se"))
  coef_vec <- coef_mat[setdiff(rownames(coef_mat), "(Intercept)"), 1L]
  coef_vec <- coef_vec[feature_names]
  selected_counts <- stats::setNames(rep(0L, length(feature_names)), feature_names)
  n_success <- 0L
  n_bootstrap <- max(0L, as.integer(n_bootstrap))
  if (n_bootstrap > 0L) {
    set.seed(seed + 1L)
    for (b in seq_len(n_bootstrap)) {
      idx <- stratified_subsample(y, sample_fraction)
      fit <- tryCatch(
        glmnet::glmnet(
          x[idx, , drop = FALSE],
          y[idx],
          family = "binomial",
          alpha = alpha,
          lambda = cv$lambda.1se,
          standardize = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      co <- as.matrix(stats::coef(fit))
      co <- co[setdiff(rownames(co), "(Intercept)"), 1L]
      co <- co[feature_names]
      selected_counts[abs(co) > 0] <- selected_counts[abs(co) > 0] + 1L
      n_success <- n_success + 1L
    }
  }

  info$elastic_net_coef <- as.numeric(coef_vec)
  info$selection_frequency <- if (n_success > 0L) as.numeric(selected_counts) / n_success else NA_real_
  attr(info, "status") <- data.frame(
    model = "elastic_net_logistic",
    status = "ok",
    reason = "",
    alpha = alpha,
    lambda_min = cv$lambda.min,
    lambda_1se = cv$lambda.1se,
    nfolds = nfolds,
    n_bootstrap_requested = n_bootstrap,
    n_bootstrap_success = n_success,
    sample_fraction = sample_fraction,
    stringsAsFactors = FALSE
  )
  info
}

safe_cor_result <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(list(estimate = NA_real_, p = NA_real_))
  }
  out <- tryCatch(
    suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out)) return(list(estimate = NA_real_, p = NA_real_))
  list(
    estimate = suppressWarnings(as.numeric(out$estimate[[1L]])),
    p = suppressWarnings(as.numeric(out$p.value))
  )
}

regression_r2 <- function(y, pred) {
  keep <- is.finite(y) & is.finite(pred)
  y <- as.numeric(y[keep])
  pred <- as.numeric(pred[keep])
  if (length(y) < 2L) return(NA_real_)
  sst <- sum((y - mean(y))^2)
  if (!is.finite(sst) || sst <= 0) return(NA_real_)
  1 - sum((y - pred)^2) / sst
}

regression_rmse <- function(y, pred) {
  keep <- is.finite(y) & is.finite(pred)
  if (!any(keep)) return(NA_real_)
  sqrt(mean((as.numeric(y[keep]) - as.numeric(pred[keep]))^2))
}

safe_lm_term <- function(df, formula, term) {
  fit <- tryCatch(
    suppressWarnings(stats::lm(formula, data = df)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(coef = NA_real_, t = NA_real_, p = NA_real_, r2 = NA_real_, pred = rep(NA_real_, nrow(df))))
  }
  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  coef <- tval <- p <- NA_real_
  if (!is.null(co) && term %in% rownames(co)) {
    coef <- suppressWarnings(as.numeric(co[term, "Estimate"]))
    tval <- suppressWarnings(as.numeric(co[term, "t value"]))
    p <- suppressWarnings(as.numeric(co[term, "Pr(>|t|)"]))
  }
  pred <- tryCatch(suppressWarnings(stats::predict(fit)), error = function(e) rep(NA_real_, nrow(df)))
  list(coef = coef, t = tval, p = p, r2 = regression_r2(df$y, pred), pred = pred)
}

continuous_direction <- function(coef) {
  coef <- suppressWarnings(as.numeric(coef))
  ifelse(
    is.finite(coef),
    ifelse(coef >= 0, "higher_ploidy_with_higher_feature", "lower_ploidy_with_higher_feature"),
    NA_character_
  )
}

continuous_main_effect_table <- function(best_df, transformed_df, z_main, y, params, log10_params) {
  rows <- lapply(params, function(param) {
    x <- as.numeric(z_main[, param])
    xt <- transformed_df[[param]]
    xo <- suppressWarnings(as.numeric(best_df[[param]]))
    lm <- safe_lm_term(data.frame(y = y, x = x), y ~ x, "x")
    sp <- safe_cor_result(x, y, method = "spearman")
    pe <- safe_cor_result(x, y, method = "pearson")
    data.frame(
      feature_name = param,
      feature_type = "main",
      parameter = param,
      parameter_a = param,
      parameter_b = NA_character_,
      transformed_label = parameter_axis_label(param, log10_params),
      n_seed = length(y),
      response_mean = mean(y, na.rm = TRUE),
      response_sd = stats::sd(y, na.rm = TRUE),
      feature_mean_original = mean(xo, na.rm = TRUE),
      feature_sd_original = stats::sd(xo, na.rm = TRUE),
      feature_mean_transformed = mean(xt, na.rm = TRUE),
      feature_sd_transformed = stats::sd(xt, na.rm = TRUE),
      lm_coef_per_z = lm$coef,
      lm_t = lm$t,
      lm_p = lm$p,
      lm_r2 = lm$r2,
      spearman_rho = sp$estimate,
      spearman_p = sp$p,
      pearson_r = pe$estimate,
      pearson_p = pe$p,
      direction = continuous_direction(lm$coef),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$lm_fdr <- stats::p.adjust(out$lm_p, method = "BH")
  out$spearman_fdr <- stats::p.adjust(out$spearman_p, method = "BH")
  out$pearson_fdr <- stats::p.adjust(out$pearson_p, method = "BH")
  out[order(-out$lm_r2, out$lm_p, out$parameter), , drop = FALSE]
}

continuous_interaction_table <- function(z_main, y) {
  design <- build_interaction_design(z_main)
  pairs <- design$interaction_pairs
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pa <- pairs$parameter_a[[i]]
    pb <- pairs$parameter_b[[i]]
    xi <- as.numeric(design$interaction_matrix[, i])
    df <- data.frame(
      y = y,
      xa = as.numeric(z_main[, pa]),
      xb = as.numeric(z_main[, pb]),
      xi = xi
    )
    full <- safe_lm_term(df, y ~ xa + xb + xi, "xi")
    reduced_fit <- tryCatch(
      suppressWarnings(stats::lm(y ~ xa + xb, data = df)),
      error = function(e) NULL
    )
    reduced_pred <- if (is.null(reduced_fit)) rep(NA_real_, nrow(df)) else suppressWarnings(stats::predict(reduced_fit))
    reduced_r2 <- regression_r2(y, reduced_pred)
    single <- safe_lm_term(data.frame(y = y, xi = xi), y ~ xi, "xi")
    sp <- safe_cor_result(xi, y, method = "spearman")
    data.frame(
      feature_name = pairs$feature_name[[i]],
      feature_type = "interaction",
      parameter = NA_character_,
      parameter_a = pa,
      parameter_b = pb,
      n_seed = length(y),
      interaction_mean = mean(xi, na.rm = TRUE),
      interaction_sd = stats::sd(xi, na.rm = TRUE),
      interaction_lm_coef = full$coef,
      interaction_lm_t = full$t,
      interaction_lm_p = full$p,
      interaction_full_r2 = full$r2,
      main_effects_r2 = reduced_r2,
      interaction_delta_r2 = if (is.finite(full$r2) && is.finite(reduced_r2)) full$r2 - reduced_r2 else NA_real_,
      interaction_single_feature_r2 = single$r2,
      interaction_spearman_rho = sp$estimate,
      interaction_spearman_p = sp$p,
      direction = continuous_direction(full$coef),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$interaction_lm_fdr <- stats::p.adjust(out$interaction_lm_p, method = "BH")
  out$interaction_spearman_fdr <- stats::p.adjust(out$interaction_spearman_p, method = "BH")
  out[order(-out$interaction_delta_r2, out$interaction_lm_p, out$feature_name), , drop = FALSE]
}

random_subsample <- function(n, fraction) {
  n_sub <- max(3L, floor(n * fraction))
  sample(seq_len(n), min(n, n_sub))
}

fit_glmnet_regression_importance <- function(z_main,
                                             y,
                                             n_bootstrap = 100L,
                                             sample_fraction = 0.75,
                                             alpha = 0.5,
                                             seed = 123L) {
  inter <- build_interaction_design(z_main)
  x <- cbind(z_main, inter$interaction_matrix)
  feature_names <- colnames(x)
  info <- data.frame(
    feature_name = feature_names,
    feature_type = ifelse(feature_names %in% colnames(z_main), "main", "interaction"),
    parameter = ifelse(feature_names %in% colnames(z_main), feature_names, NA_character_),
    parameter_a = NA_character_,
    parameter_b = NA_character_,
    stringsAsFactors = FALSE
  )
  info$parameter_a[info$feature_type == "main"] <- info$feature_name[info$feature_type == "main"]
  info$parameter_b[info$feature_type == "main"] <- NA_character_
  idx <- match(info$feature_name, inter$interaction_pairs$feature_name)
  hit <- !is.na(idx)
  info$parameter_a[hit] <- inter$interaction_pairs$parameter_a[idx[hit]]
  info$parameter_b[hit] <- inter$interaction_pairs$parameter_b[idx[hit]]

  y <- as.numeric(y)
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    info$elastic_net_coef <- NA_real_
    info$selection_frequency <- NA_real_
    attr(info, "status") <- data.frame(
      model = "elastic_net_gaussian",
      status = "unavailable",
      reason = "R package glmnet is not installed",
      stringsAsFactors = FALSE
    )
    return(info)
  }

  nfolds <- max(3L, min(10L, length(y)))
  if (length(y) < 3L || stats::sd(y) == 0) {
    info$elastic_net_coef <- NA_real_
    info$selection_frequency <- NA_real_
    attr(info, "status") <- data.frame(
      model = "elastic_net_gaussian",
      status = "unavailable",
      reason = "fewer than three finite seeds or zero response variance",
      stringsAsFactors = FALSE
    )
    return(info)
  }

  set.seed(seed)
  cv <- tryCatch(
    glmnet::cv.glmnet(x, y, family = "gaussian", alpha = alpha, standardize = FALSE, nfolds = nfolds),
    error = function(e) e
  )
  if (inherits(cv, "error")) {
    info$elastic_net_coef <- NA_real_
    info$selection_frequency <- NA_real_
    attr(info, "status") <- data.frame(
      model = "elastic_net_gaussian",
      status = "failed",
      reason = conditionMessage(cv),
      stringsAsFactors = FALSE
    )
    return(info)
  }

  coef_mat <- as.matrix(stats::coef(cv, s = "lambda.1se"))
  coef_vec <- coef_mat[setdiff(rownames(coef_mat), "(Intercept)"), 1L]
  coef_vec <- coef_vec[feature_names]
  selected_counts <- stats::setNames(rep(0L, length(feature_names)), feature_names)
  n_success <- 0L
  n_bootstrap <- max(0L, as.integer(n_bootstrap))
  if (n_bootstrap > 0L) {
    set.seed(seed + 1L)
    for (b in seq_len(n_bootstrap)) {
      idx <- random_subsample(length(y), sample_fraction)
      fit <- tryCatch(
        glmnet::glmnet(
          x[idx, , drop = FALSE],
          y[idx],
          family = "gaussian",
          alpha = alpha,
          lambda = cv$lambda.1se,
          standardize = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      co <- as.matrix(stats::coef(fit))
      co <- co[setdiff(rownames(co), "(Intercept)"), 1L]
      co <- co[feature_names]
      selected_counts[abs(co) > 0] <- selected_counts[abs(co) > 0] + 1L
      n_success <- n_success + 1L
    }
  }

  info$elastic_net_coef <- as.numeric(coef_vec)
  info$selection_frequency <- if (n_success > 0L) as.numeric(selected_counts) / n_success else NA_real_
  attr(info, "status") <- data.frame(
    model = "elastic_net_gaussian",
    status = "ok",
    reason = "",
    alpha = alpha,
    lambda_min = cv$lambda.min,
    lambda_1se = cv$lambda.1se,
    nfolds = nfolds,
    n_bootstrap_requested = n_bootstrap,
    n_bootstrap_success = n_success,
    sample_fraction = sample_fraction,
    cvm_min = min(cv$cvm, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  info
}

continuous_fold_ids <- function(y, k = 10L, seed = 123L) {
  y <- as.numeric(y)
  n <- length(y)
  k <- max(2L, min(as.integer(k), n))
  set.seed(seed)
  ord <- order(y, stats::runif(n), na.last = TRUE)
  folds <- integer(n)
  folds[ord] <- rep(seq_len(k), length.out = n)
  folds
}

fit_regression_for_matrix <- function(x, y, folds) {
  x <- as.matrix(x)
  colnames(x) <- paste0("feature_", seq_len(ncol(x)))
  df <- data.frame(y = as.numeric(y), x, check.names = FALSE)
  formula <- stats::reformulate(colnames(x), response = "y")
  fit <- tryCatch(
    suppressWarnings(stats::lm(formula, data = df)),
    error = function(e) e
  )
  apparent_pred <- rep(NA_real_, length(y))
  fit_status <- "ok"
  fit_reason <- ""
  if (inherits(fit, "error")) {
    fit_status <- "failed"
    fit_reason <- conditionMessage(fit)
  } else {
    apparent_pred <- suppressWarnings(stats::predict(fit))
  }

  cv_pred <- rep(NA_real_, length(y))
  cv_success <- 0L
  for (fold in sort(unique(folds))) {
    train <- folds != fold
    test <- folds == fold
    fold_fit <- tryCatch(
      suppressWarnings(stats::lm(formula, data = df[train, , drop = FALSE])),
      error = function(e) NULL
    )
    if (is.null(fold_fit)) next
    cv_pred[test] <- suppressWarnings(stats::predict(fold_fit, newdata = df[test, , drop = FALSE]))
    cv_success <- cv_success + 1L
  }
  list(
    apparent_r2 = regression_r2(y, apparent_pred),
    apparent_rmse = regression_rmse(y, apparent_pred),
    cv_r2 = regression_r2(y, cv_pred),
    cv_rmse = regression_rmse(y, cv_pred),
    cv_success_folds = cv_success,
    apparent_pred = apparent_pred,
    cv_pred = cv_pred,
    fit_status = fit_status,
    fit_reason = fit_reason
  )
}

continuous_feature_r2_table <- function(combined, z_main, y) {
  if (!nrow(combined)) return(data.frame())
  rows <- lapply(seq_len(nrow(combined)), function(i) {
    row <- combined[i, , drop = FALSE]
    x <- feature_vector_from_row(row, z_main)
    lm <- safe_lm_term(data.frame(y = y, x = x), y ~ x, "x")
    data.frame(
      feature_name = row$feature_name,
      feature_type = row$feature_type,
      parameter = row$parameter %||% NA_character_,
      parameter_a = row$parameter_a %||% NA_character_,
      parameter_b = row$parameter_b %||% NA_character_,
      rank = suppressWarnings(as.integer(row$rank)),
      single_feature_coef = lm$coef,
      single_feature_r2 = lm$r2,
      direction = continuous_direction(lm$coef),
      elastic_net_coef = suppressWarnings(as.numeric(row$elastic_net_coef %||% NA_real_)),
      selection_frequency = suppressWarnings(as.numeric(row$selection_frequency %||% NA_real_)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$single_feature_r2, out$rank, out$feature_name), , drop = FALSE]
  out$r2_rank <- seq_len(nrow(out))
  out
}

plot_feature_r2_bar <- function(r2_df,
                                output_prefix,
                                title,
                                top_n = NULL,
                                width = 7.8,
                                height = NULL,
                                flip = TRUE) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(character()))
  df <- r2_df[is.finite(r2_df$single_feature_r2), , drop = FALSE]
  df <- df[order(-df$single_feature_r2, df$feature_name), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n) && top_n > 0L) {
    df <- head(df, as.integer(top_n))
  }
  if (!nrow(df)) return(invisible(character()))
  df$feature_index <- seq_len(nrow(df))
  if (is.null(height)) height <- max(4.2, min(30, 2.2 + 0.16 * nrow(df)))
  p <- ggplot2::ggplot(df, ggplot2::aes(fill = direction)) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = feature_index - 0.36,
        xmax = feature_index + 0.36,
        ymin = 0,
        ymax = single_feature_r2
      )
    )
  if (isTRUE(flip)) p <- p + ggplot2::coord_flip()
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = df$feature_index,
      labels = df$feature_name,
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, max(0.01, max(df$single_feature_r2, na.rm = TRUE) * 1.05))) +
    ggplot2::scale_fill_manual(
      values = c(
        higher_ploidy_with_higher_feature = "#1f78b4",
        lower_ploidy_with_higher_feature = "#e31a1c"
      ),
      na.value = "grey70",
      labels = c(
        higher_ploidy_with_higher_feature = "positive",
        lower_ploidy_with_higher_feature = "negative"
      )
    ) +
    ggplot2::labs(
      title = title,
      x = if (isTRUE(flip)) NULL else "Feature",
      y = expression(R^2),
      fill = "Slope"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = if (isTRUE(flip)) {
        ggplot2::element_text()
      } else {
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
      },
      legend.position = if (isTRUE(flip)) "right" else "bottom"
    )
  save_plot_pair(p, output_prefix, width = width, height = height)
  invisible(paste0(output_prefix, c(".pdf", ".png")))
}

plot_top_feature_continuous <- function(combined,
                                        z_main,
                                        y,
                                        mode_reference_o2,
                                        output_dir,
                                        n = 3L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE) || !nrow(combined)) return(invisible(character()))
  top <- head(combined[order(combined$rank), , drop = FALSE], as.integer(n))
  if (!nrow(top)) return(invisible(character()))
  written <- character()
  rows <- list()
  for (i in seq_len(nrow(top))) {
    row <- top[i, , drop = FALSE]
    score <- feature_vector_from_row(row, z_main)
    lm <- safe_lm_term(data.frame(y = y, x = score), y ~ x, "x")
    plot_df <- data.frame(
      feature_value = score,
      dominant_mean_ploidy = y,
      feature_rank = i,
      feature_name = row$feature_name[[1L]],
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- plot_df
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = feature_value, y = dominant_mean_ploidy)) +
      ggplot2::geom_point(size = 1.75, alpha = 0.78, color = "#1f78b4") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#222222", linewidth = 0.75) +
      ggplot2::geom_hline(yintercept = 2, linetype = "dashed", color = "#b2182b") +
      ggplot2::labs(
        title = reference_o2_plot_title(paste0("Top", i, " continuous feature, reference "), mode_reference_o2),
        subtitle = paste0(row$feature_name[[1L]], "; R2=", formatC(lm$r2, digits = 3L, format = "fg")),
        x = row$feature_name[[1L]],
        y = "Dominant mean ploidy"
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_text(size = 8.5)
      )
    prefix <- file.path(output_dir, paste0("top", i, "_feature_dominant_ploidy"))
    save_plot_pair(p, prefix, width = 4.8, height = 4.2)
    written <- c(written, paste0(prefix, c(".pdf", ".png")))
  }
  if (length(rows)) write_csv(rbind_fill_plain(rows), file.path(output_dir, "top3_single_feature_points.csv"))
  invisible(written)
}

write_continuous_r2_outputs <- function(combined,
                                        z_main,
                                        y,
                                        mode_reference_o2,
                                        ref_dir) {
  r2_dir <- file.path(ref_dir, "r2_plots")
  dir.create(r2_dir, recursive = TRUE, showWarnings = FALSE)
  r2_tab <- continuous_feature_r2_table(combined, z_main, y)
  write_csv(r2_tab, file.path(r2_dir, "dominant_ploidy_parameter_r2_by_feature.csv"))
  main_r2 <- r2_tab[r2_tab$feature_type == "main", , drop = FALSE]
  if (nrow(main_r2)) {
    plot_feature_r2_bar(
      main_r2,
      file.path(r2_dir, "main_feature_r2_bar"),
      reference_o2_plot_title("Main feature R2, reference ", mode_reference_o2)
    )
  }
  if (nrow(r2_tab)) {
    plot_feature_r2_bar(
      r2_tab,
      file.path(r2_dir, "all_feature_r2_bar"),
      reference_o2_plot_title("All feature R2, reference ", mode_reference_o2)
    )
    plot_feature_r2_bar(
      r2_tab,
      file.path(r2_dir, "combined_feature_r2_bar_top30"),
      reference_o2_plot_title("Top30 retained-feature R2, reference ", mode_reference_o2),
      top_n = 30L,
      width = 10.5,
      height = 5.25,
      flip = FALSE
    )
    plot_top_feature_continuous(
      combined,
      z_main,
      y,
      mode_reference_o2,
      file.path(ref_dir, "top3_feature_continuous"),
      n = 3L
    )
  }
  r2_tab
}

continuous_top_feature_summary_fields <- function(combined, r2_tab, n = 3L) {
  top <- head(combined[order(combined$rank), , drop = FALSE], n)
  out <- data.frame(.row = 1L)
  for (i in seq_len(n)) {
    prefix <- paste0("top", i)
    if (i <= nrow(top)) {
      feature <- top$feature_name[[i]]
      r2_row <- r2_tab[r2_tab$feature_name == feature, , drop = FALSE]
      out[[paste0(prefix, "_feature")]] <- feature
      out[[paste0(prefix, "_feature_type")]] <- top$feature_type[[i]]
      out[[paste0(prefix, "_r2")]] <- if (nrow(r2_row)) r2_row$single_feature_r2[[1L]] else NA_real_
      out[[paste0(prefix, "_direction")]] <- if (nrow(r2_row)) r2_row$direction[[1L]] else NA_character_
      out[[paste0(prefix, "_elastic_net_coef")]] <- suppressWarnings(as.numeric(top$elastic_net_coef[[i]]))
      out[[paste0(prefix, "_selection_frequency")]] <- suppressWarnings(as.numeric(top$selection_frequency[[i]]))
    } else {
      out[[paste0(prefix, "_feature")]] <- NA_character_
      out[[paste0(prefix, "_feature_type")]] <- NA_character_
      out[[paste0(prefix, "_r2")]] <- NA_real_
      out[[paste0(prefix, "_direction")]] <- NA_character_
      out[[paste0(prefix, "_elastic_net_coef")]] <- NA_real_
      out[[paste0(prefix, "_selection_frequency")]] <- NA_real_
    }
  }
  out$.row <- NULL
  out
}

fit_top3_joint_regression <- function(combined,
                                      z_main,
                                      y,
                                      mode_reference_o2,
                                      output_dir,
                                      cv_k = 10L,
                                      seed = 123L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  top3 <- head(combined[order(combined$rank), , drop = FALSE], 3L)
  if (nrow(top3) < 1L) return(NULL)
  x <- do.call(cbind, lapply(seq_len(nrow(top3)), function(i) feature_vector_from_row(top3[i, , drop = FALSE], z_main)))
  colnames(x) <- paste0("top", seq_len(ncol(x)))
  folds <- continuous_fold_ids(y, k = cv_k, seed = seed)
  metrics <- fit_regression_for_matrix(x, y, folds)
  pred_tab <- data.frame(
    row_id = seq_along(y),
    mode_reference_o2 = mode_reference_o2,
    dominant_mean_ploidy = y,
    apparent_predicted_dominant_mean_ploidy = metrics$apparent_pred,
    cv_predicted_dominant_mean_ploidy = metrics$cv_pred,
    fold_id = folds,
    stringsAsFactors = FALSE
  )
  write_csv(pred_tab, file.path(output_dir, "top3_joint_predictions.csv"))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_df <- rbind_fill_plain(list(
      data.frame(
        observed = y,
        predicted = metrics$apparent_pred,
        evaluation = "Apparent",
        stringsAsFactors = FALSE
      ),
      data.frame(
        observed = y,
        predicted = metrics$cv_pred,
        evaluation = "Cross-validated",
        stringsAsFactors = FALSE
      )
    ))
    plot_df <- plot_df[is.finite(plot_df$observed) & is.finite(plot_df$predicted), , drop = FALSE]
    if (nrow(plot_df)) {
      lim <- range(c(plot_df$observed, plot_df$predicted), na.rm = TRUE)
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = observed, y = predicted, color = evaluation)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
        ggplot2::geom_point(alpha = 0.72, size = 1.6) +
        ggplot2::coord_equal(xlim = lim, ylim = lim, expand = TRUE) +
        ggplot2::scale_color_manual(values = c("Apparent" = "grey45", "Cross-validated" = "#1f78b4")) +
        ggplot2::labs(
          title = reference_o2_plot_title("Top3 joint regression, reference ", mode_reference_o2),
          x = "Observed dominant mean ploidy",
          y = "Predicted dominant mean ploidy",
          color = NULL
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
      save_plot_pair(p, file.path(output_dir, "top3_joint_observed_vs_predicted"), width = 5.2, height = 4.6)
    }
  }

  top_names <- c(top3$feature_name, rep(NA_character_, 3L - nrow(top3)))
  summary <- data.frame(
    mode_reference_o2 = mode_reference_o2,
    n_seed = length(y),
    top1_feature = top_names[[1L]],
    top2_feature = top_names[[2L]],
    top3_feature = top_names[[3L]],
    fit_status = metrics$fit_status,
    fit_reason = metrics$fit_reason,
    apparent_r2 = metrics$apparent_r2,
    apparent_rmse = metrics$apparent_rmse,
    cv_r2 = metrics$cv_r2,
    cv_rmse = metrics$cv_rmse,
    cv_k = length(unique(folds)),
    cv_success_folds = metrics$cv_success_folds,
    stringsAsFactors = FALSE
  )
  write_csv(summary, file.path(output_dir, "top3_joint_regression_metrics.csv"))
  summary
}

fit_cumulative_feature_r2 <- function(combined,
                                      z_main,
                                      y,
                                      mode_reference_o2,
                                      output_dir,
                                      top_n = 30L,
                                      cv_k = 10L,
                                      seed = 123L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  top_n <- max(1L, min(as.integer(top_n), nrow(combined)))
  ranked <- combined[order(combined$rank), , drop = FALSE]
  ranked <- ranked[seq_len(top_n), , drop = FALSE]
  folds <- continuous_fold_ids(y, k = cv_k, seed = seed)
  rows <- lapply(seq_len(nrow(ranked)), function(k) {
    selected <- ranked[seq_len(k), , drop = FALSE]
    x <- cumulative_feature_matrix(selected, z_main)
    metrics <- fit_regression_for_matrix(x, y, folds)
    data.frame(
      mode_reference_o2 = mode_reference_o2,
      k = k,
      added_feature = selected$feature_name[[k]],
      added_feature_type = selected$feature_type[[k]],
      feature_set = paste(selected$feature_name, collapse = ";"),
      n_seed = length(y),
      cv_k = length(unique(folds)),
      cv_success_folds = metrics$cv_success_folds,
      apparent_r2 = metrics$apparent_r2,
      apparent_rmse = metrics$apparent_rmse,
      cv_r2 = metrics$cv_r2,
      cv_rmse = metrics$cv_rmse,
      fit_status = metrics$fit_status,
      fit_reason = metrics$fit_reason,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$apparent_r2_gain <- c(NA_real_, diff(out$apparent_r2))
  out$cv_r2_gain <- c(NA_real_, diff(out$cv_r2))
  write_csv(out, file.path(output_dir, "cumulative_feature_r2.csv"))
  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(out)) {
    plot_df <- rbind_fill_plain(list(
      data.frame(k = out$k, r2 = out$cv_r2, evaluation = "Cross-validated", stringsAsFactors = FALSE),
      data.frame(k = out$k, r2 = out$apparent_r2, evaluation = "Apparent", stringsAsFactors = FALSE)
    ))
    plot_df <- plot_df[is.finite(plot_df$k) & is.finite(plot_df$r2), , drop = FALSE]
    if (nrow(plot_df)) {
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = k, y = r2, color = evaluation, linetype = evaluation)) +
        ggplot2::geom_line(linewidth = 0.85) +
        ggplot2::geom_point(size = 1.6) +
        ggplot2::geom_vline(xintercept = c(3, 5, 10), linetype = "dotted", color = "grey60") +
        ggplot2::scale_color_manual(values = c("Cross-validated" = "#1f78b4", "Apparent" = "grey45")) +
        ggplot2::scale_linetype_manual(values = c("Cross-validated" = "solid", "Apparent" = "dashed")) +
        ggplot2::labs(
          title = reference_o2_plot_title("Cumulative feature R2, reference ", mode_reference_o2),
          x = "Cumulative feature count",
          y = expression(R^2),
          color = NULL,
          linetype = NULL
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
      save_plot_pair(p, file.path(output_dir, "cumulative_feature_r2"), width = 6.4, height = 4.6)
    }
  }
  out
}

plot_top_interactions_continuous <- function(top_interactions,
                                             transformed_df,
                                             y,
                                             output_dir,
                                             log10_params,
                                             top_n = 6L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    write_tsv_plain(
      data.frame(status = "unavailable", reason = "R package ggplot2 is not installed", stringsAsFactors = FALSE),
      file.path(output_dir, "phase_plot_status.tsv")
    )
    return(character())
  }
  top_interactions <- top_interactions[seq_len(min(nrow(top_interactions), top_n)), , drop = FALSE]
  written <- character()
  for (i in seq_len(nrow(top_interactions))) {
    pa <- top_interactions$parameter_a[[i]]
    pb <- top_interactions$parameter_b[[i]]
    plot_df <- data.frame(
      x = transformed_df[[pa]],
      y = transformed_df[[pb]],
      dominant_mean_ploidy = y,
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = dominant_mean_ploidy)) +
      ggplot2::geom_point(size = 1.8, alpha = 0.82) +
      ggplot2::scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 2, na.value = "grey70") +
      ggplot2::labs(
        title = paste0("Dominant-ploidy phase plot: ", pa, " x ", pb),
        subtitle = bquote("Best in vivo parameters colored by fixed-" * O[2] * " dominant mean ploidy"),
        x = parameter_axis_label(pa, log10_params),
        y = parameter_axis_label(pb, log10_params),
        color = "Dominant\nploidy"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    prefix <- file.path(output_dir, sprintf("%02d_%s__%s", i, clean_filename(pa), clean_filename(pb)))
    save_plot_pair(p, prefix, width = 5.6, height = 4.8)
    written <- c(written, paste0(prefix, c(".pdf", ".png")))
  }
  written
}

gini_impurity <- function(y) {
  if (!length(y)) return(0)
  p <- mean(y == 1L)
  2 * p * (1 - p)
}

best_tree_split <- function(x_df, y, rows, min_node) {
  parent <- gini_impurity(y[rows])
  best <- NULL
  for (param in names(x_df)) {
    x <- x_df[[param]][rows]
    finite <- is.finite(x)
    if (sum(finite) < 2L * min_node) next
    ux <- sort(unique(x[finite]))
    if (length(ux) < 2L) next
    thresholds <- (utils::head(ux, -1L) + utils::tail(ux, -1L)) / 2
    if (length(thresholds) > 80L) {
      thresholds <- unique(as.numeric(stats::quantile(thresholds, probs = seq(0.05, 0.95, length.out = 80L), na.rm = TRUE)))
    }
    for (thr in thresholds) {
      left <- rows[is.finite(x_df[[param]][rows]) & x_df[[param]][rows] <= thr]
      right <- rows[is.finite(x_df[[param]][rows]) & x_df[[param]][rows] > thr]
      if (length(left) < min_node || length(right) < min_node) next
      weighted <- (length(left) * gini_impurity(y[left]) + length(right) * gini_impurity(y[right])) / length(rows)
      gain <- parent - weighted
      if (is.null(best) || gain > best$gain) {
        best <- list(parameter = param, threshold = thr, gain = gain, left = left, right = right)
      }
    }
  }
  best
}

build_rule_tree <- function(x_df,
                            y,
                            log10_params,
                            max_depth = 3L,
                            min_node = 20L) {
  leaves <- list()
  recurse <- function(rows, depth, path) {
    if (depth >= max_depth || length(unique(y[rows])) < 2L || length(rows) < 2L * min_node) {
      leaves[[length(leaves) + 1L]] <<- list(rows = rows, path = path)
      return(invisible(NULL))
    }
    split <- best_tree_split(x_df, y, rows, min_node)
    if (is.null(split) || !is.finite(split$gain) || split$gain <= 0) {
      leaves[[length(leaves) + 1L]] <<- list(rows = rows, path = path)
      return(invisible(NULL))
    }
    base <- split$parameter
    label <- parameter_axis_label(base, log10_params)
    threshold_original <- if (base %in% log10_params) 10^split$threshold else split$threshold
    left_cond <- data.frame(
      parameter = base,
      transformed_label = label,
      operator = "<=",
      threshold_transformed = split$threshold,
      threshold_original = threshold_original,
      gini_gain = split$gain,
      stringsAsFactors = FALSE
    )
    right_cond <- left_cond
    right_cond$operator <- ">"
    recurse(split$left, depth + 1L, rbind(path, left_cond))
    recurse(split$right, depth + 1L, rbind(path, right_cond))
  }
  recurse(seq_along(y), 0L, data.frame())
  rows <- lapply(seq_along(leaves), function(i) {
    leaf <- leaves[[i]]
    yy <- y[leaf$rows]
    rule <- if (!nrow(leaf$path)) {
      "TRUE"
    } else {
      paste(
        sprintf(
          "%s %s %.6g",
          leaf$path$transformed_label,
          leaf$path$operator,
          leaf$path$threshold_transformed
        ),
        collapse = " & "
      )
    }
    data.frame(
      rule_id = i,
      rule = rule,
      n_seed = length(yy),
      n_mode1 = sum(yy == 1L),
      n_mode2 = sum(yy == 0L),
      fraction_mode1 = mean(yy == 1L),
      predicted_mode = if (mean(yy == 1L) >= 0.5) "mode1" else "mode2",
      depth = if (!nrow(leaf$path)) 0L else nrow(leaf$path),
      split_parameters = if (!nrow(leaf$path)) "" else paste(unique(leaf$path$parameter), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$n_seed, -abs(out$fraction_mode1 - 0.5), out$rule_id), , drop = FALSE]
}

plot_top_interactions <- function(top_interactions,
                                  transformed_df,
                                  mode_df,
                                  output_dir,
                                  log10_params,
                                  top_n = 6L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    write_tsv_plain(
      data.frame(status = "unavailable", reason = "R package ggplot2 is not installed", stringsAsFactors = FALSE),
      file.path(output_dir, "phase_plot_status.tsv")
    )
    return(character())
  }
  top_interactions <- top_interactions[seq_len(min(nrow(top_interactions), top_n)), , drop = FALSE]
  written <- character()
  for (i in seq_len(nrow(top_interactions))) {
    pa <- top_interactions$parameter_a[[i]]
    pb <- top_interactions$parameter_b[[i]]
    plot_df <- data.frame(
      x = transformed_df[[pa]],
      y = transformed_df[[pb]],
      mode_label = mode_df$mode_label,
      seed = mode_df$seed,
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = mode_label)) +
      ggplot2::geom_point(size = 1.8, alpha = 0.82) +
      ggplot2::scale_color_manual(values = c(mode1 = "#1f78b4", mode2 = "#e31a1c"), na.value = "grey60") +
      ggplot2::labs(
        title = paste0("Mode phase plot: ", pa, " x ", pb),
        subtitle = bquote("Best in vivo parameters colored by fixed-" * O[2] * " mode"),
        x = parameter_axis_label(pa, log10_params),
        y = parameter_axis_label(pb, log10_params),
        color = "Mode"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    prefix <- file.path(output_dir, sprintf("%02d_%s__%s", i, clean_filename(pa), clean_filename(pb)))
    save_plot_pair(p, prefix, width = 5.6, height = 4.8)
    written <- c(written, paste0(prefix, c(".pdf", ".png")))
  }
  written
}

reference_contribution_dir <- function(output_dir, mode_reference_o2) {
  file.path(output_dir, paste0("reference_o2_", fixed_o2_o2_slug(mode_reference_o2)))
}

reference_contribution_paths <- function(output_dir, mode_reference_o2) {
  ref_dir <- reference_contribution_dir(output_dir, mode_reference_o2)
  list(
    reference_dir = ref_dir,
    seed_labels = file.path(ref_dir, "mode_parameter_seed_labels.csv"),
    main_effects = file.path(ref_dir, "mode_parameter_main_effects.csv"),
    pairwise_interactions = file.path(ref_dir, "mode_parameter_pairwise_interactions.csv"),
    stability_selection = file.path(ref_dir, "mode_parameter_stability_selection.csv"),
    elastic_net_status = file.path(ref_dir, "mode_parameter_elastic_net_status.tsv"),
    feature_importance = file.path(ref_dir, "mode_parameter_feature_importance.csv"),
    decision_rules = file.path(ref_dir, "mode_parameter_decision_rules.csv"),
    auc_by_feature = file.path(ref_dir, "auc_plots", "mode_parameter_auc_by_feature.csv"),
    top3_joint_auc = file.path(ref_dir, "top3_joint", "top3_joint_auc.csv"),
    top_feature_rocs = file.path(ref_dir, "top3_feature_roc", paste0("top", 1:3, "_feature_roc.png")),
    cumulative_auc = file.path(ref_dir, "cumulative_auc", "cumulative_feature_auc.csv"),
    summary = file.path(ref_dir, "mode_parameter_contribution_summary.tsv")
  )
}

reference_contribution_complete <- function(output_dir, mode_reference_o2) {
  paths <- reference_contribution_paths(output_dir, mode_reference_o2)
  all(file.exists(c(
    paths$feature_importance,
    paths$summary,
    paths$auc_by_feature,
    paths$top3_joint_auc,
    paths$top_feature_rocs,
    paths$cumulative_auc
  )))
}

read_reference_contribution <- function(output_dir, mode_reference_o2) {
  paths <- reference_contribution_paths(output_dir, mode_reference_o2)
  if (!file.exists(paths$summary)) stop("Missing reference summary: ", paths$summary)
  if (!file.exists(paths$feature_importance)) stop("Missing reference feature importance: ", paths$feature_importance)
  summary <- read_tsv(paths$summary)
  top <- head(read_csv_plain(paths$feature_importance), 30L)
  list(summary = summary, top_features = top, output_dir = paths$reference_dir)
}

dominant_ploidy_contribution_paths <- function(output_dir, mode_reference_o2) {
  ref_dir <- reference_contribution_dir(output_dir, mode_reference_o2)
  list(
    reference_dir = ref_dir,
    seed_values = file.path(ref_dir, "dominant_ploidy_parameter_seed_values.csv"),
    main_effects = file.path(ref_dir, "dominant_ploidy_parameter_main_effects.csv"),
    pairwise_interactions = file.path(ref_dir, "dominant_ploidy_parameter_pairwise_interactions.csv"),
    stability_selection = file.path(ref_dir, "dominant_ploidy_parameter_stability_selection.csv"),
    elastic_net_status = file.path(ref_dir, "dominant_ploidy_parameter_elastic_net_status.tsv"),
    feature_importance = file.path(ref_dir, "dominant_ploidy_parameter_feature_importance.csv"),
    r2_by_feature = file.path(ref_dir, "r2_plots", "dominant_ploidy_parameter_r2_by_feature.csv"),
    top3_joint_metrics = file.path(ref_dir, "top3_joint", "top3_joint_regression_metrics.csv"),
    top_feature_plots = file.path(ref_dir, "top3_feature_continuous", paste0("top", 1:3, "_feature_dominant_ploidy.png")),
    cumulative_r2 = file.path(ref_dir, "cumulative_r2", "cumulative_feature_r2.csv"),
    summary = file.path(ref_dir, "dominant_ploidy_parameter_contribution_summary.tsv")
  )
}

dominant_ploidy_contribution_complete <- function(output_dir, mode_reference_o2) {
  paths <- dominant_ploidy_contribution_paths(output_dir, mode_reference_o2)
  all(file.exists(c(
    paths$feature_importance,
    paths$summary,
    paths$r2_by_feature,
    paths$top3_joint_metrics,
    paths$top_feature_plots,
    paths$cumulative_r2
  )))
}

read_dominant_ploidy_contribution <- function(output_dir, mode_reference_o2) {
  paths <- dominant_ploidy_contribution_paths(output_dir, mode_reference_o2)
  if (!file.exists(paths$summary)) stop("Missing reference summary: ", paths$summary)
  if (!file.exists(paths$feature_importance)) stop("Missing reference feature importance: ", paths$feature_importance)
  summary <- read_tsv(paths$summary)
  top <- head(read_csv_plain(paths$feature_importance), 30L)
  list(summary = summary, top_features = top, output_dir = paths$reference_dir)
}

mode_reference_o2_group <- function(o2) {
  o2 <- suppressWarnings(as.numeric(o2))
  group <- ifelse(o2 <= 0.1, "Low O2", ifelse(o2 <= 1, "Mid O2", "High O2"))
  factor(group, levels = c("Low O2", "Mid O2", "High O2"))
}

mode_reference_o2_label <- function(o2) {
  paste0("O2=", format_o2_value(o2))
}

write_o2_group_auc_summary <- function(output_dir, reference_o2, index) {
  tables_dir <- file.path(output_dir, "summary_tables")
  plots_dir <- file.path(output_dir, "summary_plots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  top3_rows <- lapply(seq_len(nrow(index)), function(i) {
    row <- index[i, , drop = FALSE]
    o2 <- suppressWarnings(as.numeric(row$mode_reference_o2[[1L]]))
    rows <- lapply(seq_len(3L), function(rank) {
      feature_col <- paste0("top", rank, "_feature")
      type_col <- paste0("top", rank, "_feature_type")
      auc_col <- paste0("top", rank, "_auc")
      direction_col <- paste0("top", rank, "_direction")
      coef_col <- paste0("top", rank, "_elastic_net_coef")
      stability_col <- paste0("top", rank, "_selection_frequency")
      data.frame(
        mode_reference_o2 = o2,
        o2_group = as.character(mode_reference_o2_group(o2)),
        feature_rank = rank,
        feature_name = if (feature_col %in% names(row)) as.character(row[[feature_col]][[1L]]) else NA_character_,
        feature_type = if (type_col %in% names(row)) as.character(row[[type_col]][[1L]]) else NA_character_,
        auc = if (auc_col %in% names(row)) suppressWarnings(as.numeric(row[[auc_col]][[1L]])) else NA_real_,
        direction = if (direction_col %in% names(row)) as.character(row[[direction_col]][[1L]]) else NA_character_,
        elastic_net_coef = if (coef_col %in% names(row)) suppressWarnings(as.numeric(row[[coef_col]][[1L]])) else NA_real_,
        selection_frequency = if (stability_col %in% names(row)) suppressWarnings(as.numeric(row[[stability_col]][[1L]])) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  top3_long <- rbind_fill_plain(top3_rows)
  top3_long <- top3_long[is.finite(top3_long$mode_reference_o2), , drop = FALSE]
  write_csv(top3_long, file.path(tables_dir, "mode_parameter_top3_auc_by_reference_o2.csv"))

  joint_rows <- lapply(reference_o2, function(o2) {
    paths <- reference_contribution_paths(output_dir, o2)
    if (!file.exists(paths$top3_joint_auc)) return(NULL)
    out <- read_csv_plain(paths$top3_joint_auc)
    out$mode_reference_o2 <- suppressWarnings(as.numeric(out$mode_reference_o2))
    out$o2_group <- as.character(mode_reference_o2_group(out$mode_reference_o2))
    out
  })
  joint <- rbind_fill_plain(joint_rows)
  write_csv(joint, file.path(tables_dir, "top3_joint_auc_across_reference_o2.csv"))

  plot_rows <- list()
  if (nrow(top3_long)) {
    top1 <- top3_long[top3_long$feature_rank == 1L, , drop = FALSE]
    if (nrow(top1)) {
      plot_rows[[length(plot_rows) + 1L]] <- data.frame(
        mode_reference_o2 = top1$mode_reference_o2,
        o2_group = top1$o2_group,
        metric = "best_feature_auc",
        metric_label = "Best feature AUC",
        auc = top1$auc,
        evaluation = "single_feature",
        stringsAsFactors = FALSE
      )
    }
    mean_rows <- stats::aggregate(
      auc ~ mode_reference_o2 + o2_group,
      data = top3_long[is.finite(top3_long$auc), , drop = FALSE],
      FUN = mean
    )
    if (nrow(mean_rows)) {
      mean_rows$metric <- "top3_mean_single_feature_auc"
      mean_rows$metric_label <- "Top3 mean feature AUC"
      mean_rows$evaluation <- "single_feature_mean"
      plot_rows[[length(plot_rows) + 1L]] <- mean_rows
    }
  }
  if (nrow(joint)) {
    joint_auc <- if ("cv_auc" %in% names(joint)) suppressWarnings(as.numeric(joint$cv_auc)) else rep(NA_real_, nrow(joint))
    apparent_auc <- if ("apparent_auc" %in% names(joint)) suppressWarnings(as.numeric(joint$apparent_auc)) else rep(NA_real_, nrow(joint))
    joint_plot <- data.frame(
      mode_reference_o2 = suppressWarnings(as.numeric(joint$mode_reference_o2)),
      o2_group = as.character(mode_reference_o2_group(joint$mode_reference_o2)),
      metric = "top3_joint_auc",
      metric_label = "Top3 joint AUC",
      auc = ifelse(is.finite(joint_auc), joint_auc, apparent_auc),
      evaluation = ifelse(is.finite(joint_auc), "cross_validated", "apparent"),
      stringsAsFactors = FALSE
    )
    plot_rows[[length(plot_rows) + 1L]] <- joint_plot
  }

  summary_long <- rbind_fill_plain(plot_rows)
  if (nrow(summary_long)) {
    summary_long <- summary_long[is.finite(summary_long$mode_reference_o2) & is.finite(summary_long$auc), , drop = FALSE]
  }
  write_csv(summary_long, file.path(tables_dir, "o2_group_auc_summary.csv"))
  if (nrow(summary_long) && requireNamespace("ggplot2", quietly = TRUE)) {
    ref_levels <- mode_reference_o2_label(reference_o2)
    summary_long$mode_reference_o2_label <- factor(mode_reference_o2_label(summary_long$mode_reference_o2), levels = ref_levels)
    summary_long$o2_group <- factor(summary_long$o2_group, levels = c("Low O2", "Mid O2", "High O2"))
    summary_long$metric_label <- factor(
      summary_long$metric_label,
      levels = c("Best feature AUC", "Top3 mean feature AUC", "Top3 joint AUC")
    )
    p <- ggplot2::ggplot(
      summary_long,
      ggplot2::aes(
        x = mode_reference_o2_label,
        y = auc,
        fill = metric_label
      )
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.76), width = 0.68) +
      ggplot2::coord_cartesian(ylim = c(0.5, 1)) +
      ggplot2::scale_y_continuous(breaks = seq(0.5, 1, by = 0.1)) +
      ggplot2::scale_x_discrete(labels = o2_discrete_label_expressions) +
      ggplot2::scale_fill_manual(values = c(
        "Best feature AUC" = "#1f78b4",
        "Top3 mean feature AUC" = "#33a02c",
        "Top3 joint AUC" = "#e31a1c"
      )) +
      ggplot2::labs(
        title = bquote("Mode contribution AUC across reference " * O[2]),
        x = expression("Reference " * O[2] * " (%)"),
        y = "AUC",
        fill = NULL
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
        legend.position = "bottom"
      )
    save_plot_pair(p, file.path(plots_dir, "o2_group_auc_summary"), width = 6.2, height = 6.2)

    top3_plot <- top3_long[is.finite(top3_long$auc), , drop = FALSE]
    if (nrow(top3_plot)) {
      top3_plot$mode_reference_o2_label <- factor(mode_reference_o2_label(top3_plot$mode_reference_o2), levels = ref_levels)
      top3_plot$o2_group <- factor(top3_plot$o2_group, levels = c("Low O2", "Mid O2", "High O2"))
      top3_plot$feature_rank_label <- factor(paste0("Top", top3_plot$feature_rank), levels = paste0("Top", 1:3))
      p_top3 <- ggplot2::ggplot(
        top3_plot,
        ggplot2::aes(
          x = mode_reference_o2_label,
          y = auc,
          fill = feature_rank_label
        )
      ) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.76), width = 0.68) +
        ggplot2::coord_cartesian(ylim = c(0.5, 1)) +
        ggplot2::scale_y_continuous(breaks = seq(0.5, 1, by = 0.1)) +
        ggplot2::scale_x_discrete(labels = o2_discrete_label_expressions) +
        ggplot2::scale_fill_manual(values = c(Top1 = "#1f78b4", Top2 = "#33a02c", Top3 = "#ff7f00")) +
        ggplot2::labs(
          title = bquote("Top3 retained feature AUC across reference " * O[2]),
          x = expression("Reference " * O[2] * " (%)"),
          y = "AUC",
          fill = NULL
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
          legend.position = "bottom"
        )
      save_plot_pair(p_top3, file.path(plots_dir, "o2_group_top3_feature_auc"), width = 6.2, height = 6.2)
    }
  }
  invisible(file.path(tables_dir, "o2_group_auc_summary.csv"))
}

write_cumulative_auc_global_summary <- function(output_dir, reference_o2) {
  tables_dir <- file.path(output_dir, "summary_tables")
  plots_dir <- file.path(output_dir, "summary_plots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- lapply(reference_o2, function(o2) {
    paths <- reference_contribution_paths(output_dir, o2)
    if (!file.exists(paths$cumulative_auc)) return(NULL)
    tab <- read_csv_plain(paths$cumulative_auc)
    tab$mode_reference_o2 <- suppressWarnings(as.numeric(tab$mode_reference_o2))
    tab$o2_group <- as.character(mode_reference_o2_group(tab$mode_reference_o2))
    tab$mode_reference_o2_label <- mode_reference_o2_label(tab$mode_reference_o2)
    tab
  })
  cumulative <- rbind_fill_plain(rows)
  write_csv(cumulative, file.path(tables_dir, "cumulative_feature_auc_across_reference_o2.csv"))
  if (!nrow(cumulative) || !requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(file.path(tables_dir, "cumulative_feature_auc_across_reference_o2.csv")))
  }

  plot_df <- cumulative[is.finite(cumulative$k) & is.finite(cumulative$cv_auc), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(file.path(tables_dir, "cumulative_feature_auc_across_reference_o2.csv")))
  ref_levels <- mode_reference_o2_label(reference_o2)
  plot_df$mode_reference_o2_label <- factor(plot_df$mode_reference_o2_label, levels = ref_levels)
  plot_df$o2_group <- factor(plot_df$o2_group, levels = c("Low O2", "Mid O2", "High O2"))
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = k, y = cv_auc, color = mode_reference_o2_label, group = mode_reference_o2_label)
  ) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::geom_point(size = 1.35, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = c(3, 5, 10), linetype = "dotted", color = "grey60") +
    ggplot2::scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
    ggplot2::scale_color_discrete(labels = o2_discrete_label_expressions) +
    ggplot2::labs(
      title = bquote("Cumulative cross-validated AUC across reference " * O[2]),
      x = "Cumulative feature count",
      y = "Cross-validated AUC",
      color = expression("Reference " * O[2])
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  save_plot_pair(p, file.path(plots_dir, "cumulative_feature_auc_summary"), width = 6.2, height = 6.2)
  invisible(file.path(tables_dir, "cumulative_feature_auc_across_reference_o2.csv"))
}

fixed_o2_attractor_distribution_paths <- function(mode_tables_dir) {
  fig_dir <- file.path(mode_tables_dir, "Figures")
  reference_grouped_dir <- file.path(fig_dir, "ReferenceO2ModeGroupedAnalyticalSolution")
  list(
    figure_dir = fig_dir,
    reference_grouped_dir = reference_grouped_dir,
    all_seeds_prefix = file.path(fig_dir, "fixed_o2_analytical_solution_all_seeds_violin_boxplot"),
    by_mode_prefix = file.path(fig_dir, "fixed_o2_analytical_solution_by_mode_violin_boxplot")
  )
}

plot_pair_files <- function(prefix) {
  paste0(prefix, c(".pdf", ".png"))
}

plot_pair_exists <- function(prefix) {
  all(file.exists(plot_pair_files(prefix)))
}

fixed_o2_reference_grouped_prefix <- function(paths, reference_o2) {
  file.path(
    paths$reference_grouped_dir,
    paste0("fixed_o2_analytical_solution_grouped_by_reference_o2_", fixed_o2_o2_slug(reference_o2))
  )
}

fixed_o2_transition_table_paths <- function(mode_tables_dir, reference_o2 = 5) {
  tables_dir <- file.path(mode_tables_dir, "Tables")
  prefix <- file.path(
    tables_dir,
    paste0("reference_o2_", fixed_o2_o2_slug(reference_o2), "_mode_transition")
  )
  list(
    tables_dir = tables_dir,
    counts_by_o2 = paste0(prefix, "_counts_by_o2.csv"),
    paired_seed_delta = paste0(prefix, "_paired_seed_delta.csv"),
    mean_ploidy_and_delta_summary = paste0(prefix, "_mean_ploidy_and_delta_summary.csv"),
    mean_ploidy_and_delta_comparison_wide = paste0(prefix, "_mean_ploidy_and_delta_comparison_wide.csv"),
    delta_direction_comparison = paste0(prefix, "_delta_direction_comparison.csv"),
    sankey_seed_paths = paste0(prefix, "_sankey_seed_paths.csv"),
    sankey_path_counts = paste0(prefix, "_sankey_path_counts.csv"),
    sankey_node_summary = paste0(prefix, "_sankey_node_summary.csv")
  )
}

fixed_o2_transition_figure_paths <- function(mode_tables_dir, reference_o2 = 5) {
  figures_dir <- file.path(
    mode_tables_dir,
    "Figures",
    paste0("ModeTransitionReferenceO2_", fixed_o2_o2_slug(reference_o2))
  )
  prefix <- file.path(
    figures_dir,
    paste0("reference_o2_", fixed_o2_o2_slug(reference_o2), "_mode_transition")
  )
  list(
    figures_dir = figures_dir,
    counts_by_o2 = paste0(prefix, "_counts_by_o2"),
    fraction_by_o2 = paste0(prefix, "_fraction_by_o2"),
    mean_ploidy_by_o2 = paste0(prefix, "_mean_ploidy_by_o2"),
    paired_delta_by_o2 = paste0(prefix, "_paired_delta_by_o2"),
    sankey_by_o2 = paste0(prefix, "_sankey_by_o2"),
    sankey_pairwise_combined = paste0(prefix, "_sankey_pairwise_from_reference_combined"),
    index = file.path(figures_dir, "reference_o2_mode_transition_figure_index.csv")
  )
}

fixed_o2_attractor_reference_o2_values <- function(mode_tables_dir) {
  path <- file.path(mode_tables_dir, "fixed_o2_attractor_mode_by_seed_o2.tsv")
  if (!file.exists(path)) return(numeric())
  tab <- read_tsv(path)
  if (!"O2_pct" %in% names(tab)) return(numeric())
  vals <- suppressWarnings(as.numeric(tab$O2_pct))
  sort(unique(vals[is.finite(vals)]))
}

fixed_o2_transition_tables_complete <- function(mode_tables_dir, reference_o2 = 5) {
  reference_values <- fixed_o2_attractor_reference_o2_values(mode_tables_dir)
  if (!length(reference_values) || !any(abs(reference_values - reference_o2) < 1e-9)) return(TRUE)
  paths <- fixed_o2_transition_table_paths(mode_tables_dir, reference_o2 = reference_o2)
  all(file.exists(unlist(paths[setdiff(names(paths), "tables_dir")], use.names = FALSE)))
}

fixed_o2_transition_figures_complete <- function(mode_tables_dir, reference_o2 = 5) {
  reference_values <- fixed_o2_attractor_reference_o2_values(mode_tables_dir)
  if (!length(reference_values) || !any(abs(reference_values - reference_o2) < 1e-9)) return(TRUE)
  paths <- fixed_o2_transition_figure_paths(mode_tables_dir, reference_o2 = reference_o2)
  all(file.exists(c(
    plot_pair_files(paths$counts_by_o2),
    plot_pair_files(paths$fraction_by_o2),
    plot_pair_files(paths$mean_ploidy_by_o2),
    plot_pair_files(paths$paired_delta_by_o2),
    plot_pair_files(paths$sankey_by_o2),
    plot_pair_files(paths$sankey_pairwise_combined),
    paths$index
  )))
}

fixed_o2_attractor_distribution_complete <- function(mode_tables_dir) {
  paths <- fixed_o2_attractor_distribution_paths(mode_tables_dir)
  reference_o2 <- fixed_o2_attractor_reference_o2_values(mode_tables_dir)
  reference_files <- unlist(lapply(reference_o2, function(o2) {
    plot_pair_files(fixed_o2_reference_grouped_prefix(paths, o2))
  }), use.names = FALSE)
  all(file.exists(c(
    plot_pair_files(paths$all_seeds_prefix),
    plot_pair_files(paths$by_mode_prefix),
    reference_files
  ))) &&
    fixed_o2_transition_tables_complete(mode_tables_dir, reference_o2 = 5) &&
    fixed_o2_transition_figures_complete(mode_tables_dir, reference_o2 = 5)
}

prepare_fixed_o2_attractor_distribution_data <- function(mode_tables_dir) {
  path <- file.path(mode_tables_dir, "fixed_o2_attractor_mode_by_seed_o2.tsv")
  if (!file.exists(path)) stop("Missing fixed-O2 attractor table: ", path)
  tab <- read_tsv(path)
  required <- c("seed_id", "O2_pct", "dominant_mean_ploidy", "mode_label")
  missing <- setdiff(required, names(tab))
  if (length(missing)) stop("Missing columns in fixed-O2 attractor table: ", paste(missing, collapse = ", "))
  tab$fixed_o2 <- suppressWarnings(as.numeric(tab$O2_pct))
  tab$analytical_solution <- suppressWarnings(as.numeric(tab$dominant_mean_ploidy))
  tab <- tab[is.finite(tab$fixed_o2) & is.finite(tab$analytical_solution), , drop = FALSE]
  if (!nrow(tab)) stop("No finite fixed-O2 attractor analytical solution rows were available.")
  o2_levels <- sort(unique(tab$fixed_o2))
  tab$fixed_o2_label <- factor(mode_reference_o2_label(tab$fixed_o2), levels = mode_reference_o2_label(o2_levels))
  tab$fixed_o2_index <- match(tab$fixed_o2_label, levels(tab$fixed_o2_label))
  tab$mode_label_display <- ifelse(tab$mode_label == "mode1", "Mode1", ifelse(tab$mode_label == "mode2", "Mode2", NA_character_))
  tab$mode_label_display <- factor(tab$mode_label_display, levels = c("Mode1", "Mode2"))
  tab
}

summarise_numeric_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(data.frame(
      mean = NA_real_,
      median = NA_real_,
      sd = NA_real_,
      min = NA_real_,
      max = NA_real_
    ))
  }
  data.frame(
    mean = mean(x),
    median = stats::median(x),
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
    min = min(x),
    max = max(x)
  )
}

fixed_o2_sankey_mode_col <- function(o2) {
  paste0("mode_at_o2_", fixed_o2_o2_slug(o2))
}

fixed_o2_sankey_ploidy_col <- function(o2) {
  paste0("mean_ploidy_at_o2_", fixed_o2_o2_slug(o2))
}

fixed_o2_sankey_state_col <- function(o2) {
  paste0("reference_comparison_at_o2_", fixed_o2_o2_slug(o2))
}

fixed_o2_mode_display <- function(mode_label) {
  out <- as.character(mode_label)
  out[out == "mode1"] <- "Mode1"
  out[out == "mode2"] <- "Mode2"
  out[!(out %in% c("Mode1", "Mode2"))] <- NA_character_
  out
}

fixed_o2_reference_comparison_state <- function(reference_mode, current_mode, is_reference_axis = FALSE) {
  ref <- as.character(reference_mode)
  cur <- as.character(current_mode)
  out <- if (isTRUE(is_reference_axis)) {
    paste0(ref, "_reference")
  } else {
    paste0(ref, "_to_", cur)
  }
  valid <- (ref %in% c("mode1", "mode2")) & (cur %in% c("mode1", "mode2"))
  out[!valid] <- NA_character_
  out
}

fixed_o2_reference_comparison_display <- function(state) {
  out <- as.character(state)
  out[out == "mode1_reference"] <- "Mode1"
  out[out == "mode2_reference"] <- "Mode2"
  out[out == "mode1_to_mode1"] <- "Mode1->Mode1"
  out[out == "mode1_to_mode2"] <- "Mode1->Mode2"
  out[out == "mode2_to_mode1"] <- "Mode2->Mode1"
  out[out == "mode2_to_mode2"] <- "Mode2->Mode2"
  out[!(out %in% c("Mode1", "Mode2", "Mode1->Mode1", "Mode1->Mode2", "Mode2->Mode1", "Mode2->Mode2"))] <- NA_character_
  out
}

fixed_o2_reference_comparison_current_mode <- function(state) {
  state <- as.character(state)
  out <- rep(NA_character_, length(state))
  out[state %in% c("mode1_reference", "mode1_to_mode1", "mode2_to_mode1")] <- "Mode1"
  out[state %in% c("mode2_reference", "mode1_to_mode2", "mode2_to_mode2")] <- "Mode2"
  out
}

fixed_o2_reference_comparison_fill_key <- function(state) {
  out <- as.character(state)
  out[out == "mode1_reference"] <- "Mode1"
  out[out == "mode2_reference"] <- "Mode2"
  out[out == "mode1_to_mode1"] <- "M1->M1"
  out[out == "mode2_to_mode2"] <- "M2->M2"
  out[out == "mode1_to_mode2"] <- "M1->M2"
  out[out == "mode2_to_mode1"] <- "M2->M1"
  out[!(out %in% c("Mode1", "Mode2", "M1->M1", "M2->M2", "M1->M2", "M2->M1"))] <- NA_character_
  out
}

fixed_o2_reference_comparison_fill_colors <- function() {
  c(
    Mode1 = "#1f78b4",
    Mode2 = "#ff7f00",
    "M1->M1" = "#1f78b4",
    "M2->M2" = "#ff7f00",
    "M1->M2" = "#2ca25f",
    "M2->M1" = "#7b3294"
  )
}

fixed_o2_reference_comparison_state_order <- function(is_reference_axis = FALSE) {
  if (isTRUE(is_reference_axis)) {
    return(c("mode2_reference", "mode1_reference"))
  }
  c("mode2_to_mode2", "mode1_to_mode2", "mode2_to_mode1", "mode1_to_mode1")
}

fixed_o2_reference_sankey_o2_order <- function(o2_values, reference_o2 = 5) {
  numeric_o2 <- suppressWarnings(as.numeric(o2_values))
  vals <- sort(unique(numeric_o2[is.finite(numeric_o2)]), decreasing = TRUE)
  if (!length(vals) || !any(abs(vals - reference_o2) < 1e-9)) return(numeric())
  reference_value <- vals[which.min(abs(vals - reference_o2))]
  lower_values <- vals[vals < reference_o2 - 1e-9]
  c(reference_value, lower_values)
}

build_fixed_o2_reference_sankey_tables <- function(plot_df, reference_o2 = 5) {
  o2_order <- fixed_o2_reference_sankey_o2_order(plot_df$fixed_o2, reference_o2 = reference_o2)
  if (!length(o2_order)) return(NULL)
  base_cols <- c("seed_id", "mode_label", "analytical_solution")
  seed_paths <- NULL
  mode_cols <- character()
  ploidy_cols <- character()
  state_cols <- character()
  for (current_o2 in o2_order) {
    cur <- plot_df[
      is.finite(plot_df$fixed_o2) &
        abs(plot_df$fixed_o2 - current_o2) < 1e-9 &
        plot_df$mode_label %in% c("mode1", "mode2") &
        is.finite(plot_df$analytical_solution),
      base_cols,
      drop = FALSE
    ]
    cur <- cur[!duplicated(cur$seed_id), , drop = FALSE]
    mode_col <- fixed_o2_sankey_mode_col(current_o2)
    ploidy_col <- fixed_o2_sankey_ploidy_col(current_o2)
    names(cur) <- c("seed_id", mode_col, ploidy_col)
    seed_paths <- if (is.null(seed_paths)) {
      cur
    } else {
      merge(seed_paths, cur, by = "seed_id", all = FALSE, sort = FALSE)
    }
    mode_cols <- c(mode_cols, mode_col)
    ploidy_cols <- c(ploidy_cols, ploidy_col)
    state_cols <- c(state_cols, fixed_o2_sankey_state_col(current_o2))
  }
  if (is.null(seed_paths) || !nrow(seed_paths)) return(NULL)
  complete_modes <- Reduce(`&`, lapply(mode_cols, function(col) seed_paths[[col]] %in% c("mode1", "mode2")))
  seed_paths <- seed_paths[complete_modes, , drop = FALSE]
  if (!nrow(seed_paths)) return(NULL)
  seed_paths$reference_o2 <- reference_o2
  seed_paths$reference_mode_label <- seed_paths[[mode_cols[[1L]]]]
  seed_paths$reference_mode_display <- fixed_o2_mode_display(seed_paths$reference_mode_label)
  for (i in seq_along(mode_cols)) {
    seed_paths[[state_cols[[i]]]] <- fixed_o2_reference_comparison_state(
      seed_paths$reference_mode_label,
      seed_paths[[mode_cols[[i]]]],
      is_reference_axis = i == 1L
    )
  }
  seed_paths$mode_path_key <- apply(seed_paths[, mode_cols, drop = FALSE], 1L, paste, collapse = "|")
  seed_paths$mode_path_display <- apply(
    seed_paths[, mode_cols, drop = FALSE],
    1L,
    function(x) paste(fixed_o2_mode_display(x), collapse = " -> ")
  )
  seed_paths$reference_comparison_path_key <- apply(seed_paths[, state_cols, drop = FALSE], 1L, paste, collapse = "|")
  seed_paths$reference_comparison_path_display <- apply(
    seed_paths[, state_cols, drop = FALSE],
    1L,
    function(x) paste(fixed_o2_reference_comparison_display(x), collapse = " -> ")
  )
  seed_paths <- seed_paths[, c(
    "reference_o2",
    "seed_id",
    "reference_mode_label",
    "reference_mode_display",
    "mode_path_key",
    "mode_path_display",
    "reference_comparison_path_key",
    "reference_comparison_path_display",
    mode_cols,
    state_cols,
    ploidy_cols
  ), drop = FALSE]

  split_idx <- split(seq_len(nrow(seed_paths)), seed_paths$reference_comparison_path_key)
  path_rows <- lapply(names(split_idx), function(reference_comparison_path_key) {
    idx <- split_idx[[reference_comparison_path_key]]
    first_row <- seed_paths[idx[[1L]], , drop = FALSE]
    out <- data.frame(
      reference_o2 = reference_o2,
      n_seeds = length(idx),
      fraction_seeds = length(idx) / nrow(seed_paths),
      reference_mode_label = first_row$reference_mode_label[[1L]],
      reference_mode_display = first_row$reference_mode_display[[1L]],
      reference_comparison_path_key = reference_comparison_path_key,
      reference_comparison_path_display = first_row$reference_comparison_path_display[[1L]],
      mode_path_key = first_row$mode_path_key[[1L]],
      mode_path_display = first_row$mode_path_display[[1L]],
      stringsAsFactors = FALSE
    )
    for (col in mode_cols) out[[col]] <- first_row[[col]][[1L]]
    for (col in state_cols) out[[col]] <- first_row[[col]][[1L]]
    out
  })
  path_counts <- rbind_fill_plain(path_rows)
  if (nrow(path_counts)) {
    path_counts <- path_counts[order(
      match(path_counts$reference_mode_label, c("mode1", "mode2")),
      -path_counts$n_seeds,
      path_counts$reference_comparison_path_key
    ), , drop = FALSE]
    path_counts$path_rank <- seq_len(nrow(path_counts))
    path_counts <- path_counts[, c(
      "reference_o2",
      "path_rank",
      "n_seeds",
      "fraction_seeds",
      "reference_mode_label",
      "reference_mode_display",
      "reference_comparison_path_key",
      "reference_comparison_path_display",
      "mode_path_key",
      "mode_path_display",
      mode_cols,
      state_cols
    ), drop = FALSE]
  }

  reference_ploidy_col <- ploidy_cols[[1L]]
  node_summary_rows <- list()
  for (axis_idx in seq_along(state_cols)) {
    state_levels <- fixed_o2_reference_comparison_state_order(is_reference_axis = axis_idx == 1L)
    for (state in state_levels) {
      seed_idx <- seed_paths[[state_cols[[axis_idx]]]] == state
      n_state <- sum(seed_idx, na.rm = TRUE)
      current_ploidy <- suppressWarnings(as.numeric(seed_paths[[ploidy_cols[[axis_idx]]]][seed_idx]))
      reference_ploidy <- suppressWarnings(as.numeric(seed_paths[[reference_ploidy_col]][seed_idx]))
      current_mean <- mean(current_ploidy, na.rm = TRUE)
      reference_mean <- mean(reference_ploidy, na.rm = TRUE)
      paired_delta <- current_ploidy - reference_ploidy
      delta_mean <- mean(paired_delta, na.rm = TRUE)
      if (!is.finite(current_mean)) current_mean <- NA_real_
      if (!is.finite(reference_mean)) reference_mean <- NA_real_
      if (!is.finite(delta_mean)) delta_mean <- NA_real_
      node_summary_rows[[length(node_summary_rows) + 1L]] <- data.frame(
        reference_o2 = reference_o2,
        axis_idx = axis_idx,
        current_o2 = o2_order[[axis_idx]],
        comparison_state = state,
        comparison_display = fixed_o2_reference_comparison_display(state),
        mode_display = fixed_o2_reference_comparison_current_mode(state),
        is_mode_transition = axis_idx > 1L && state %in% c("mode1_to_mode2", "mode2_to_mode1"),
        n_seeds = n_state,
        mean_ploidy_current_o2 = current_mean,
        mean_ploidy_at_reference_o2_same_seeds = reference_mean,
        mean_paired_delta_current_minus_reference_o2 = delta_mean,
        stringsAsFactors = FALSE
      )
    }
  }
  node_summary <- rbind_fill_plain(node_summary_rows)

  list(
    seed_paths = seed_paths,
    path_counts = path_counts,
    node_summary = node_summary,
    o2_order = o2_order,
    mode_cols = mode_cols,
    state_cols = state_cols,
    ploidy_cols = ploidy_cols
  )
}

write_fixed_o2_reference_transition_tables <- function(plot_df,
                                                       mode_tables_dir,
                                                       reference_o2 = 5,
                                                       overwrite = FALSE) {
  o2_values <- sort(unique(plot_df$fixed_o2[is.finite(plot_df$fixed_o2)]))
  if (!length(o2_values) || !any(abs(o2_values - reference_o2) < 1e-9)) {
    message("Skipping fixed-O2 transition tables; reference O2=", format_o2_value(reference_o2), " is not available.")
    return(invisible(NULL))
  }
  paths <- fixed_o2_transition_table_paths(mode_tables_dir, reference_o2 = reference_o2)
  if (!overwrite && fixed_o2_transition_tables_complete(mode_tables_dir, reference_o2 = reference_o2)) {
    message("Skipping fixed-O2 transition tables; existing outputs found under: ", paths$tables_dir)
    return(invisible(paths))
  }
  dir.create(paths$tables_dir, recursive = TRUE, showWarnings = FALSE)

  base_cols <- c("seed_id", "fixed_o2", "mode_label", "analytical_solution")
  ref <- plot_df[
    is.finite(plot_df$fixed_o2) &
      abs(plot_df$fixed_o2 - reference_o2) < 1e-9 &
      plot_df$mode_label %in% c("mode1", "mode2") &
      is.finite(plot_df$analytical_solution),
    base_cols,
    drop = FALSE
  ]
  ref <- ref[!duplicated(ref$seed_id), , drop = FALSE]
  names(ref) <- c("seed_id", "reference_o2", "mode_label_reference_o2", "mean_ploidy_reference_o2")
  if (!nrow(ref)) {
    message("Skipping fixed-O2 transition tables; no mode1/mode2 reference seeds at O2=", format_o2_value(reference_o2), ".")
    return(invisible(NULL))
  }

  transition_levels <- c("mode2_at_reference_o2_to_mode1", "mode1_at_reference_o2_to_mode2")
  count_rows <- list()
  seed_rows <- list()
  for (current_o2 in o2_values) {
    cur <- plot_df[
      is.finite(plot_df$fixed_o2) &
        abs(plot_df$fixed_o2 - current_o2) < 1e-9 &
        plot_df$mode_label %in% c("mode1", "mode2") &
        is.finite(plot_df$analytical_solution),
      base_cols,
      drop = FALSE
    ]
    cur <- cur[!duplicated(cur$seed_id), , drop = FALSE]
    names(cur) <- c("seed_id", "current_o2", "mode_label_current_o2", "mean_ploidy_current_o2")
    joined <- merge(ref, cur, by = "seed_id", all.x = TRUE)
    ref_mode2 <- joined$mode_label_reference_o2 == "mode2"
    ref_mode1 <- joined$mode_label_reference_o2 == "mode1"
    mode2_to_mode1 <- ref_mode2 & joined$mode_label_current_o2 == "mode1"
    mode1_to_mode2 <- ref_mode1 & joined$mode_label_current_o2 == "mode2"
    count_rows[[length(count_rows) + 1L]] <- data.frame(
      reference_o2 = reference_o2,
      current_o2 = current_o2,
      n_reference_mode2 = sum(ref_mode2, na.rm = TRUE),
      mode2_at_reference_o2_to_mode1 = sum(mode2_to_mode1, na.rm = TRUE),
      mode2_to_mode1_fraction = sum(mode2_to_mode1, na.rm = TRUE) / sum(ref_mode2, na.rm = TRUE),
      mode2_at_reference_o2_still_mode2 = sum(ref_mode2 & joined$mode_label_current_o2 == "mode2", na.rm = TRUE),
      n_reference_mode1 = sum(ref_mode1, na.rm = TRUE),
      mode1_at_reference_o2_to_mode2 = sum(mode1_to_mode2, na.rm = TRUE),
      mode1_to_mode2_fraction = sum(mode1_to_mode2, na.rm = TRUE) / sum(ref_mode1, na.rm = TRUE),
      mode1_at_reference_o2_still_mode1 = sum(ref_mode1 & joined$mode_label_current_o2 == "mode1", na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    changed <- joined[mode2_to_mode1 | mode1_to_mode2, , drop = FALSE]
    if (nrow(changed)) {
      changed$reference_o2 <- reference_o2
      changed$current_o2 <- current_o2
      changed$transition <- ifelse(
        changed$mode_label_reference_o2 == "mode2" & changed$mode_label_current_o2 == "mode1",
        transition_levels[[1L]],
        transition_levels[[2L]]
      )
      changed$paired_delta_current_minus_reference_o2 <- changed$mean_ploidy_current_o2 - changed$mean_ploidy_reference_o2
      changed$paired_abs_delta <- abs(changed$paired_delta_current_minus_reference_o2)
      seed_rows[[length(seed_rows) + 1L]] <- changed[, c(
        "reference_o2",
        "current_o2",
        "seed_id",
        "transition",
        "mode_label_reference_o2",
        "mode_label_current_o2",
        "mean_ploidy_reference_o2",
        "mean_ploidy_current_o2",
        "paired_delta_current_minus_reference_o2",
        "paired_abs_delta"
      ), drop = FALSE]
    }
  }
  counts <- rbind_fill_plain(count_rows)
  seed_delta <- rbind_fill_plain(seed_rows)
  if (!nrow(seed_delta)) {
    seed_delta <- data.frame(
      reference_o2 = numeric(),
      current_o2 = numeric(),
      seed_id = character(),
      transition = character(),
      mode_label_reference_o2 = character(),
      mode_label_current_o2 = character(),
      mean_ploidy_reference_o2 = numeric(),
      mean_ploidy_current_o2 = numeric(),
      paired_delta_current_minus_reference_o2 = numeric(),
      paired_abs_delta = numeric(),
      stringsAsFactors = FALSE
    )
  }

  summary_rows <- list()
  for (current_o2 in o2_values) {
    for (transition in transition_levels) {
      df <- seed_delta[
        is.finite(seed_delta$current_o2) &
          abs(seed_delta$current_o2 - current_o2) < 1e-9 &
          seed_delta$transition == transition,
        ,
        drop = FALSE
      ]
      delta_stats <- summarise_numeric_vector(df$paired_delta_current_minus_reference_o2)
      abs_delta_stats <- summarise_numeric_vector(df$paired_abs_delta)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        reference_o2 = reference_o2,
        current_o2 = current_o2,
        transition = transition,
        n_converted = nrow(df),
        mean_ploidy_at_reference_o2_same_seeds = if (nrow(df)) mean(df$mean_ploidy_reference_o2, na.rm = TRUE) else NA_real_,
        mean_ploidy_at_current_o2_same_seeds = if (nrow(df)) mean(df$mean_ploidy_current_o2, na.rm = TRUE) else NA_real_,
        mean_paired_delta_current_minus_reference_o2 = delta_stats$mean,
        median_paired_delta_current_minus_reference_o2 = delta_stats$median,
        sd_paired_delta_current_minus_reference_o2 = delta_stats$sd,
        min_paired_delta_current_minus_reference_o2 = delta_stats$min,
        max_paired_delta_current_minus_reference_o2 = delta_stats$max,
        mean_abs_paired_delta = abs_delta_stats$mean,
        median_abs_paired_delta = abs_delta_stats$median,
        stringsAsFactors = FALSE
      )
    }
  }
  summary_tab <- rbind_fill_plain(summary_rows)

  direction_rows <- lapply(o2_values, function(current_o2) {
    a <- seed_delta[
      is.finite(seed_delta$current_o2) &
        abs(seed_delta$current_o2 - current_o2) < 1e-9 &
        seed_delta$transition == transition_levels[[1L]],
      ,
      drop = FALSE
    ]
    b <- seed_delta[
      is.finite(seed_delta$current_o2) &
        abs(seed_delta$current_o2 - current_o2) < 1e-9 &
        seed_delta$transition == transition_levels[[2L]],
      ,
      drop = FALSE
    ]
    mean_a <- if (nrow(a)) mean(a$paired_delta_current_minus_reference_o2, na.rm = TRUE) else NA_real_
    mean_b <- if (nrow(b)) mean(b$paired_delta_current_minus_reference_o2, na.rm = TRUE) else NA_real_
    abs_a <- if (nrow(a)) mean(a$paired_abs_delta, na.rm = TRUE) else NA_real_
    abs_b <- if (nrow(b)) mean(b$paired_abs_delta, na.rm = TRUE) else NA_real_
    p_delta <- if (nrow(a) > 0L && nrow(b) > 0L) {
      suppressWarnings(stats::wilcox.test(
        a$paired_delta_current_minus_reference_o2,
        b$paired_delta_current_minus_reference_o2,
        exact = FALSE
      )$p.value)
    } else {
      NA_real_
    }
    p_abs <- if (nrow(a) > 0L && nrow(b) > 0L) {
      suppressWarnings(stats::wilcox.test(a$paired_abs_delta, b$paired_abs_delta, exact = FALSE)$p.value)
    } else {
      NA_real_
    }
    data.frame(
      reference_o2 = reference_o2,
      current_o2 = current_o2,
      n_mode2_to_mode1 = nrow(a),
      n_mode1_to_mode2 = nrow(b),
      mean_paired_delta_mode2_to_mode1 = mean_a,
      mean_paired_delta_mode1_to_mode2 = mean_b,
      difference_of_mean_paired_delta_m21_minus_m12 = mean_a - mean_b,
      mean_abs_paired_delta_mode2_to_mode1 = abs_a,
      mean_abs_paired_delta_mode1_to_mode2 = abs_b,
      difference_of_mean_abs_delta_m21_minus_m12 = abs_a - abs_b,
      wilcox_p_delta_between_transition_directions = p_delta,
      wilcox_p_abs_delta_between_transition_directions = p_abs,
      stringsAsFactors = FALSE
    )
  })
  direction_comparison <- rbind_fill_plain(direction_rows)

  wide_rows <- lapply(o2_values, function(current_o2) {
    a <- summary_tab[
      is.finite(summary_tab$current_o2) &
        abs(summary_tab$current_o2 - current_o2) < 1e-9 &
        summary_tab$transition == transition_levels[[1L]],
      ,
      drop = FALSE
    ]
    b <- summary_tab[
      is.finite(summary_tab$current_o2) &
        abs(summary_tab$current_o2 - current_o2) < 1e-9 &
        summary_tab$transition == transition_levels[[2L]],
      ,
      drop = FALSE
    ]
    data.frame(
      reference_o2 = reference_o2,
      current_o2 = current_o2,
      mode2_to_mode1_n = a$n_converted,
      mode2_to_mode1_mean_ploidy_at_reference_o2 = a$mean_ploidy_at_reference_o2_same_seeds,
      mode2_to_mode1_mean_ploidy_at_current_o2 = a$mean_ploidy_at_current_o2_same_seeds,
      mode2_to_mode1_mean_paired_delta = a$mean_paired_delta_current_minus_reference_o2,
      mode1_to_mode2_n = b$n_converted,
      mode1_to_mode2_mean_ploidy_at_reference_o2 = b$mean_ploidy_at_reference_o2_same_seeds,
      mode1_to_mode2_mean_ploidy_at_current_o2 = b$mean_ploidy_at_current_o2_same_seeds,
      mode1_to_mode2_mean_paired_delta = b$mean_paired_delta_current_minus_reference_o2,
      difference_between_directions_mean_ploidy_at_reference_o2 =
        a$mean_ploidy_at_reference_o2_same_seeds - b$mean_ploidy_at_reference_o2_same_seeds,
      difference_between_directions_mean_ploidy_at_current_o2 =
        a$mean_ploidy_at_current_o2_same_seeds - b$mean_ploidy_at_current_o2_same_seeds,
      difference_between_directions_mean_paired_delta =
        a$mean_paired_delta_current_minus_reference_o2 - b$mean_paired_delta_current_minus_reference_o2,
      stringsAsFactors = FALSE
    )
  })
  wide <- rbind_fill_plain(wide_rows)

  write_csv(counts, paths$counts_by_o2)
  write_csv(seed_delta, paths$paired_seed_delta)
  write_csv(summary_tab, paths$mean_ploidy_and_delta_summary)
  write_csv(wide, paths$mean_ploidy_and_delta_comparison_wide)
  write_csv(direction_comparison, paths$delta_direction_comparison)
  sankey_tables <- build_fixed_o2_reference_sankey_tables(plot_df, reference_o2 = reference_o2)
  if (!is.null(sankey_tables)) {
    write_csv(sankey_tables$seed_paths, paths$sankey_seed_paths)
    write_csv(sankey_tables$path_counts, paths$sankey_path_counts)
    write_csv(sankey_tables$node_summary, paths$sankey_node_summary)
  }
  invisible(paths)
}

transition_display_label <- function(transition) {
  out <- as.character(transition)
  out[out == "mode2_at_reference_o2_to_mode1"] <- "Mode2 at reference O2 -> Mode1"
  out[out == "mode1_at_reference_o2_to_mode2"] <- "Mode1 at reference O2 -> Mode2"
  out
}

transition_display_label_expr <- function(transition) {
  labels <- transition_display_label(transition)
  labels <- gsub("reference O2", "reference O[2]", labels, fixed = TRUE)
  parse(text = paste0("\"", labels, "\""))
}

write_fixed_o2_reference_transition_plots <- function(mode_tables_dir, reference_o2 = 5, overwrite = FALSE) {
  table_paths <- fixed_o2_transition_table_paths(mode_tables_dir, reference_o2 = reference_o2)
  figure_paths <- fixed_o2_transition_figure_paths(mode_tables_dir, reference_o2 = reference_o2)
  if (!file.exists(table_paths$counts_by_o2) || !file.exists(table_paths$mean_ploidy_and_delta_summary)) {
    message("Skipping fixed-O2 transition plots; transition tables are not available under: ", table_paths$tables_dir)
    return(invisible(NULL))
  }
  if (!overwrite && fixed_o2_transition_figures_complete(mode_tables_dir, reference_o2 = reference_o2)) {
    message("Skipping fixed-O2 transition plots; existing outputs found under: ", figure_paths$figures_dir)
    return(invisible(figure_paths))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping fixed-O2 transition plots; ggplot2 is not available.")
    return(invisible(NULL))
  }
  dir.create(figure_paths$figures_dir, recursive = TRUE, showWarnings = FALSE)

  counts <- read_csv_plain(table_paths$counts_by_o2)
  summary_tab <- read_csv_plain(table_paths$mean_ploidy_and_delta_summary)
  counts$current_o2 <- suppressWarnings(as.numeric(counts$current_o2))
  summary_tab$current_o2 <- suppressWarnings(as.numeric(summary_tab$current_o2))
  o2_values <- sort(unique(counts$current_o2[is.finite(counts$current_o2)]))
  o2_levels <- mode_reference_o2_label(o2_values)
  counts$current_o2_label <- factor(mode_reference_o2_label(counts$current_o2), levels = o2_levels)
  summary_tab$current_o2_label <- factor(mode_reference_o2_label(summary_tab$current_o2), levels = o2_levels)
  transition_levels <- c("mode2_at_reference_o2_to_mode1", "mode1_at_reference_o2_to_mode2")
  transition_labels <- c(
    "Mode2 at reference O2 -> Mode1",
    "Mode1 at reference O2 -> Mode2"
  )
  transition_colors <- c(
    "Mode2 at reference O2 -> Mode1" = "#1f78b4",
    "Mode1 at reference O2 -> Mode2" = "#ff7f00"
  )

  counts_long <- rbind(
    data.frame(
      current_o2 = counts$current_o2,
      current_o2_label = counts$current_o2_label,
      transition = transition_labels[[1L]],
      n_converted = suppressWarnings(as.numeric(counts$mode2_at_reference_o2_to_mode1)),
      fraction_converted = suppressWarnings(as.numeric(counts$mode2_to_mode1_fraction)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      current_o2 = counts$current_o2,
      current_o2_label = counts$current_o2_label,
      transition = transition_labels[[2L]],
      n_converted = suppressWarnings(as.numeric(counts$mode1_at_reference_o2_to_mode2)),
      fraction_converted = suppressWarnings(as.numeric(counts$mode1_to_mode2_fraction)),
      stringsAsFactors = FALSE
    )
  )
  counts_long <- counts_long[is.finite(counts_long$current_o2) & counts_long$current_o2 != reference_o2, , drop = FALSE]
  counts_long$transition <- factor(counts_long$transition, levels = transition_labels)
  counts_long$count_label <- ifelse(
    is.finite(counts_long$fraction_converted),
    paste0(counts_long$n_converted, "\n", sprintf("%.1f%%", 100 * counts_long$fraction_converted)),
    as.character(counts_long$n_converted)
  )

  base_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )

  p_counts <- ggplot2::ggplot(
    counts_long,
    ggplot2::aes(x = current_o2_label, y = n_converted, fill = transition)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.76), width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(label = count_label),
      position = ggplot2::position_dodge(width = 0.76),
      vjust = -0.18,
      size = 2.8,
      lineheight = 0.92
    ) +
    ggplot2::scale_x_discrete(labels = o2_discrete_label_expressions) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
    ggplot2::scale_fill_manual(values = transition_colors, drop = FALSE) +
    ggplot2::labs(
      title = bquote("Mode transition counts from reference " * O[2] == .(format_o2_value(reference_o2))),
      x = expression("Current fixed " * O[2] * " (%)"),
      y = "Converted seeds",
      fill = "Transition"
    ) +
    base_theme
  save_plot_pair(p_counts, figure_paths$counts_by_o2, width = 7.2, height = 5.0)

  p_fraction <- ggplot2::ggplot(
    counts_long,
    ggplot2::aes(x = current_o2_label, y = fraction_converted, fill = transition)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.76), width = 0.68) +
    ggplot2::scale_x_discrete(labels = o2_discrete_label_expressions) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      limits = c(0, max(c(0.6, counts_long$fraction_converted), na.rm = TRUE)),
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::scale_fill_manual(values = transition_colors, drop = FALSE) +
    ggplot2::labs(
      title = bquote("Mode transition fractions from reference " * O[2] == .(format_o2_value(reference_o2))),
      x = expression("Current fixed " * O[2] * " (%)"),
      y = "Converted fraction",
      fill = "Transition"
    ) +
    base_theme
  save_plot_pair(p_fraction, figure_paths$fraction_by_o2, width = 7.2, height = 5.0)

  plot_summary <- summary_tab[
    summary_tab$transition %in% transition_levels &
      is.finite(summary_tab$current_o2) &
      summary_tab$current_o2 != reference_o2 &
      is.finite(summary_tab$n_converted) &
      summary_tab$n_converted > 0,
    ,
    drop = FALSE
  ]
  plot_summary$transition_label <- factor(transition_display_label(plot_summary$transition), levels = transition_labels)
  ploidy_long <- rbind(
    data.frame(
      current_o2 = plot_summary$current_o2,
      current_o2_label = plot_summary$current_o2_label,
      transition_label = plot_summary$transition_label,
      ploidy_context = "Reference O2=5",
      mean_ploidy = suppressWarnings(as.numeric(plot_summary$mean_ploidy_at_reference_o2_same_seeds)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      current_o2 = plot_summary$current_o2,
      current_o2_label = plot_summary$current_o2_label,
      transition_label = plot_summary$transition_label,
      ploidy_context = "Current O2",
      mean_ploidy = suppressWarnings(as.numeric(plot_summary$mean_ploidy_at_current_o2_same_seeds)),
      stringsAsFactors = FALSE
    )
  )
  ploidy_long$ploidy_context <- factor(ploidy_long$ploidy_context, levels = c("Reference O2=5", "Current O2"))
  facet_labels <- c(
    "Mode2 at reference O2 -> Mode1" = "Mode2~at~reference~O[2]~'->'~Mode1",
    "Mode1 at reference O2 -> Mode2" = "Mode1~at~reference~O[2]~'->'~Mode2"
  )
  p_ploidy <- ggplot2::ggplot(
    ploidy_long,
    ggplot2::aes(x = current_o2_label, y = mean_ploidy, color = ploidy_context, group = ploidy_context)
  ) +
    ggplot2::geom_hline(yintercept = 2, linetype = "dashed", color = "#4b5563", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::geom_point(size = 2.1) +
    ggplot2::facet_wrap(
      ~ transition_label,
      ncol = 1,
      labeller = ggplot2::labeller(transition_label = ggplot2::as_labeller(facet_labels, default = ggplot2::label_parsed))
    ) +
    ggplot2::scale_x_discrete(labels = o2_discrete_label_expressions) +
    ggplot2::scale_color_manual(
      values = c("Reference O2=5" = "#5b6778", "Current O2" = "#1f78b4"),
      labels = c(expression("Reference " * O[2] * "=5"), expression("Current " * O[2])),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = bquote("Mean ploidy of switching seeds from reference " * O[2] == .(format_o2_value(reference_o2))),
      x = expression("Current fixed " * O[2] * " (%)"),
      y = "Mean dominant ploidy",
      color = NULL
    ) +
    base_theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#eef3f8", color = "#d8e0e8"),
      strip.text = ggplot2::element_text(face = "bold")
    )
  save_plot_pair(p_ploidy, figure_paths$mean_ploidy_by_o2, width = 7.2, height = 6.2)

  delta_df <- plot_summary
  delta_df$transition_label <- factor(transition_display_label(delta_df$transition), levels = transition_labels)
  delta_df$se <- with(delta_df, sd_paired_delta_current_minus_reference_o2 / sqrt(pmax(n_converted, 1)))
  p_delta <- ggplot2::ggplot(
    delta_df,
    ggplot2::aes(
      x = current_o2_label,
      y = mean_paired_delta_current_minus_reference_o2,
      fill = transition_label
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#3b4652", linewidth = 0.35) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.76), width = 0.68) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_paired_delta_current_minus_reference_o2 - se,
        ymax = mean_paired_delta_current_minus_reference_o2 + se
      ),
      position = ggplot2::position_dodge(width = 0.76),
      width = 0.18,
      linewidth = 0.35
    ) +
    ggplot2::scale_x_discrete(labels = o2_discrete_label_expressions) +
    ggplot2::scale_fill_manual(values = transition_colors, drop = FALSE) +
    ggplot2::labs(
      title = bquote("Paired ploidy change in switching seeds from reference " * O[2] == .(format_o2_value(reference_o2))),
      x = expression("Current fixed " * O[2] * " (%)"),
      y = expression("Mean paired " * Delta * " ploidy (current " * O[2] * " - reference " * O[2] * ")"),
      fill = "Transition"
    ) +
    base_theme
  save_plot_pair(p_delta, figure_paths$paired_delta_by_o2, width = 7.2, height = 5.0)

  index <- data.frame(
    figure = c(
      "Mode transition counts by O2",
      "Mode transition fractions by O2",
      "Mean ploidy of switching seeds by O2",
      "Paired ploidy delta of switching seeds by O2"
    ),
    pdf = paste0(c(
      figure_paths$counts_by_o2,
      figure_paths$fraction_by_o2,
      figure_paths$mean_ploidy_by_o2,
      figure_paths$paired_delta_by_o2
    ), ".pdf"),
    png = paste0(c(
      figure_paths$counts_by_o2,
      figure_paths$fraction_by_o2,
      figure_paths$mean_ploidy_by_o2,
      figure_paths$paired_delta_by_o2
    ), ".png"),
    stringsAsFactors = FALSE
  )
  write_csv(index, figure_paths$index)
  invisible(figure_paths)
}

append_fixed_o2_transition_figure_index <- function(figure_paths, row) {
  existing <- if (file.exists(figure_paths$index)) {
    read_csv_plain(figure_paths$index)
  } else {
    data.frame(figure = character(), pdf = character(), png = character(), stringsAsFactors = FALSE)
  }
  if ("figure" %in% names(existing)) {
    existing <- existing[existing$figure != row$figure[[1L]], , drop = FALSE]
  }
  write_csv(rbind_fill_plain(list(existing, row)), figure_paths$index)
  invisible(figure_paths$index)
}

write_fixed_o2_reference_mode_sankey <- function(plot_df, mode_tables_dir, reference_o2 = 5, overwrite = FALSE) {
  table_paths <- fixed_o2_transition_table_paths(mode_tables_dir, reference_o2 = reference_o2)
  figure_paths <- fixed_o2_transition_figure_paths(mode_tables_dir, reference_o2 = reference_o2)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping fixed-O2 mode Sankey plot; ggplot2 is not available.")
    return(invisible(NULL))
  }
  sankey_tables <- build_fixed_o2_reference_sankey_tables(plot_df, reference_o2 = reference_o2)
  if (is.null(sankey_tables) || !nrow(sankey_tables$seed_paths) || !nrow(sankey_tables$path_counts)) {
    message("Skipping fixed-O2 mode Sankey plot; no complete seed-matched mode paths were available.")
    return(invisible(NULL))
  }
  dir.create(table_paths$tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
  if (overwrite || !file.exists(table_paths$sankey_seed_paths)) {
    write_csv(sankey_tables$seed_paths, table_paths$sankey_seed_paths)
  }
  if (overwrite || !file.exists(table_paths$sankey_path_counts)) {
    write_csv(sankey_tables$path_counts, table_paths$sankey_path_counts)
  }
  if (overwrite || !file.exists(table_paths$sankey_node_summary)) {
    write_csv(sankey_tables$node_summary, table_paths$sankey_node_summary)
  }

  index_row <- data.frame(
    figure = "Mode label Sankey by O2",
    pdf = paste0(figure_paths$sankey_by_o2, ".pdf"),
    png = paste0(figure_paths$sankey_by_o2, ".png"),
    stringsAsFactors = FALSE
  )
  if (!overwrite && plot_pair_exists(figure_paths$sankey_by_o2)) {
    append_fixed_o2_transition_figure_index(figure_paths, index_row)
    message("Skipping fixed-O2 mode Sankey plot; existing outputs found: ", figure_paths$sankey_by_o2)
    return(invisible(figure_paths))
  }

  groups <- sankey_tables$path_counts
  groups$n_seeds <- suppressWarnings(as.numeric(groups$n_seeds))
  groups$group_id <- seq_len(nrow(groups))
  o2_order <- sankey_tables$o2_order
  state_cols <- sankey_tables$state_cols
  ploidy_cols <- sankey_tables$ploidy_cols
  seed_paths <- sankey_tables$seed_paths
  reference_ploidy_col <- ploidy_cols[[1L]]
  mode_colors <- fixed_o2_reference_comparison_fill_colors()
  node_gap <- max(10, 0.028 * sum(groups$n_seeds, na.rm = TRUE))
  node_width <- 0.22

  node_rows <- list()
  for (axis_idx in seq_along(state_cols)) {
    y0 <- 0
    state_levels <- fixed_o2_reference_comparison_state_order(is_reference_axis = axis_idx == 1L)
    for (state in state_levels) {
      n_state <- sum(groups$n_seeds[groups[[state_cols[[axis_idx]]]] == state], na.rm = TRUE)
      seed_idx <- seed_paths[[state_cols[[axis_idx]]]] == state
      node_mean_ploidy <- mean(suppressWarnings(as.numeric(seed_paths[[ploidy_cols[[axis_idx]]]][seed_idx])), na.rm = TRUE)
      node_reference_mean_ploidy <- mean(suppressWarnings(as.numeric(seed_paths[[reference_ploidy_col]][seed_idx])), na.rm = TRUE)
      if (!is.finite(node_mean_ploidy)) node_mean_ploidy <- NA_real_
      if (!is.finite(node_reference_mean_ploidy)) node_reference_mean_ploidy <- NA_real_
      node_rows[[length(node_rows) + 1L]] <- data.frame(
        axis_idx = axis_idx,
        o2_value = o2_order[[axis_idx]],
        comparison_state = state,
        comparison_display = fixed_o2_reference_comparison_display(state),
        mode_display = fixed_o2_reference_comparison_current_mode(state),
        fill_key = fixed_o2_reference_comparison_fill_key(state),
        is_mode_transition = axis_idx > 1L && state %in% c("mode1_to_mode2", "mode2_to_mode1"),
        n_seeds = n_state,
        mean_ploidy = node_mean_ploidy,
        mean_ploidy_at_reference_o2_same_seeds = node_reference_mean_ploidy,
        ymin = y0,
        ymax = y0 + n_state,
        stringsAsFactors = FALSE
      )
      y0 <- y0 + n_state + node_gap
    }
  }
  nodes <- rbind_fill_plain(node_rows)
  nodes <- nodes[nodes$n_seeds > 0, , drop = FALSE]
  nodes$xmin <- nodes$axis_idx - node_width
  nodes$xmax <- nodes$axis_idx + node_width
  nodes$comparison_display_short <- gsub("Mode", "M", nodes$comparison_display, fixed = TRUE)
  nodes$mean_ploidy_label <- ifelse(is.finite(nodes$mean_ploidy), sprintf("%.2f", nodes$mean_ploidy), "NA")
  nodes$reference_mean_ploidy_label <- ifelse(
    is.finite(nodes$mean_ploidy_at_reference_o2_same_seeds),
    sprintf("%.2f", nodes$mean_ploidy_at_reference_o2_same_seeds),
    "NA"
  )
  base_node_label <- paste0("n=", nodes$n_seeds, "\nMean=", nodes$mean_ploidy_label)
  nodes$label <- ifelse(
    nodes$is_mode_transition,
    paste0(base_node_label, "\nRef=", nodes$reference_mean_ploidy_label),
    base_node_label
  )
  nodes$text_size <- ifelse(
    nodes$is_mode_transition & nodes$n_seeds < 35,
    1.25,
    ifelse(nodes$n_seeds < 35, 1.55, ifelse(nodes$n_seeds < 80, 1.95, 2.35))
  )
  nodes$label_inside <- ifelse(nodes$n_seeds < 35, "", nodes$label)
  small_nodes <- nodes[nodes$n_seeds < 35, , drop = FALSE]
  if (nrow(small_nodes)) {
    small_nodes$label_x <- ifelse(small_nodes$axis_idx == length(state_cols), small_nodes$xmin - 0.06, small_nodes$xmax + 0.06)
    small_nodes$segment_x <- ifelse(small_nodes$axis_idx == length(state_cols), small_nodes$xmin, small_nodes$xmax)
    small_nodes$hjust <- ifelse(small_nodes$axis_idx == length(state_cols), 1, 0)
  }

  rank_groups_for_axis <- function(df) {
    ranks <- lapply(seq_along(state_cols), function(i) {
      match(df[[state_cols[[i]]]], fixed_o2_reference_comparison_state_order(is_reference_axis = i == 1L))
    })
    do.call(order, c(ranks, list(-suppressWarnings(as.numeric(df$n_seeds)), df$reference_comparison_path_key)))
  }

  allocation_rows <- list()
  for (axis_idx in seq_along(state_cols)) {
    col <- state_cols[[axis_idx]]
    state_levels <- fixed_o2_reference_comparison_state_order(is_reference_axis = axis_idx == 1L)
    for (state in state_levels) {
      idx <- which(groups[[col]] == state)
      if (!length(idx)) next
      node <- nodes[nodes$axis_idx == axis_idx & nodes$comparison_state == state, , drop = FALSE]
      if (!nrow(node)) next
      node_groups <- groups[idx, , drop = FALSE]
      node_groups <- node_groups[rank_groups_for_axis(node_groups), , drop = FALSE]
      ymin <- node$ymin[[1L]]
      for (i in seq_len(nrow(node_groups))) {
        n <- node_groups$n_seeds[[i]]
        allocation_rows[[length(allocation_rows) + 1L]] <- data.frame(
          group_id = node_groups$group_id[[i]],
          axis_idx = axis_idx,
          comparison_state = state,
          ymin = ymin,
          ymax = ymin + n,
          stringsAsFactors = FALSE
        )
        ymin <- ymin + n
      }
    }
  }
  allocation <- rbind_fill_plain(allocation_rows)

  ribbon_rows <- list()
  curve_t <- seq(0, 1, length.out = 28)
  smooth_t <- curve_t * curve_t * (3 - 2 * curve_t)
  if (length(state_cols) > 1L) {
    for (group_id in groups$group_id) {
      group_row <- groups[groups$group_id == group_id, , drop = FALSE]
      for (axis_idx in seq_len(length(state_cols) - 1L)) {
        left <- allocation[allocation$group_id == group_id & allocation$axis_idx == axis_idx, , drop = FALSE]
        right <- allocation[allocation$group_id == group_id & allocation$axis_idx == axis_idx + 1L, , drop = FALSE]
        if (!nrow(left) || !nrow(right)) next
        right_state <- group_row[[state_cols[[axis_idx + 1L]]]][[1L]]
        x_left <- axis_idx + node_width
        x_right <- axis_idx + 1L - node_width
        x <- x_left + (x_right - x_left) * curve_t
        top <- left$ymax[[1L]] + (right$ymax[[1L]] - left$ymax[[1L]]) * smooth_t
        bottom <- left$ymin[[1L]] + (right$ymin[[1L]] - left$ymin[[1L]]) * smooth_t
        ribbon_rows[[length(ribbon_rows) + 1L]] <- data.frame(
          segment_id = paste(group_id, axis_idx, sep = "_"),
          x = c(x, rev(x)),
          y = c(top, rev(bottom)),
          fill_key = fixed_o2_reference_comparison_fill_key(right_state),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  ribbons <- rbind_fill_plain(ribbon_rows)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = ribbons,
      ggplot2::aes(x = x, y = y, group = segment_id, fill = fill_key),
      color = NA,
      alpha = 0.42
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_key),
      color = "white",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = axis_idx, y = (ymin + ymax) / 2, label = label_inside, size = text_size),
      color = "white",
      lineheight = 0.78,
      fontface = "bold"
    ) +
    ggplot2::geom_segment(
      data = small_nodes,
      ggplot2::aes(
        x = segment_x,
        xend = label_x,
        y = (ymin + ymax) / 2,
        yend = (ymin + ymax) / 2
      ),
      inherit.aes = FALSE,
      color = "#111827",
      linewidth = 0.25
    ) +
    ggplot2::geom_label(
      data = small_nodes,
      ggplot2::aes(x = label_x, y = (ymin + ymax) / 2, label = label, hjust = hjust),
      inherit.aes = FALSE,
      color = "#111827",
      fill = "white",
      linewidth = 0.15,
      label.padding = grid::unit(0.08, "lines"),
      size = 2.0,
      lineheight = 0.78,
      fontface = "bold"
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(
      breaks = seq_along(o2_order),
      labels = o2_discrete_label_expressions(mode_reference_o2_label(o2_order)),
      expand = ggplot2::expansion(add = 0.35)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.06))) +
    ggplot2::scale_fill_manual(values = mode_colors, drop = FALSE) +
    ggplot2::labs(
      title = bquote("Seed-matched Mode-label changes relative to " * O[2] == .(format_o2_value(reference_o2))),
      subtitle = bquote("Each downstream node is " * O[2] == .(format_o2_value(reference_o2)) * " Mode -> current " * O[2] * " Mode; green marks Mode1->Mode2 and purple marks Mode2->Mode1."),
      x = expression("Fixed " * O[2] * " (%)"),
      y = NULL,
      fill = "Mode change"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(10, 18, 10, 10)
    )
  save_plot_pair(p, figure_paths$sankey_by_o2, width = 9.4, height = 6.0)
  append_fixed_o2_transition_figure_index(figure_paths, index_row)
  invisible(figure_paths)
}

write_fixed_o2_reference_pairwise_sankey_combined <- function(plot_df,
                                                              mode_tables_dir,
                                                              reference_o2 = 5,
                                                              overwrite = FALSE) {
  table_paths <- fixed_o2_transition_table_paths(mode_tables_dir, reference_o2 = reference_o2)
  figure_paths <- fixed_o2_transition_figure_paths(mode_tables_dir, reference_o2 = reference_o2)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping pairwise fixed-O2 mode Sankey plot; ggplot2 is not available.")
    return(invisible(NULL))
  }
  sankey_tables <- build_fixed_o2_reference_sankey_tables(plot_df, reference_o2 = reference_o2)
  if (is.null(sankey_tables) || !nrow(sankey_tables$seed_paths) || length(sankey_tables$state_cols) < 2L) {
    message("Skipping pairwise fixed-O2 mode Sankey plot; no complete seed-matched mode paths were available.")
    return(invisible(NULL))
  }
  dir.create(table_paths$tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
  if (overwrite || !file.exists(table_paths$sankey_seed_paths)) {
    write_csv(sankey_tables$seed_paths, table_paths$sankey_seed_paths)
  }
  if (overwrite || !file.exists(table_paths$sankey_path_counts)) {
    write_csv(sankey_tables$path_counts, table_paths$sankey_path_counts)
  }
  if (overwrite || !file.exists(table_paths$sankey_node_summary)) {
    write_csv(sankey_tables$node_summary, table_paths$sankey_node_summary)
  }

  index_row <- data.frame(
    figure = "Pairwise mode label Sankey from reference O2",
    pdf = paste0(figure_paths$sankey_pairwise_combined, ".pdf"),
    png = paste0(figure_paths$sankey_pairwise_combined, ".png"),
    stringsAsFactors = FALSE
  )
  if (!overwrite && plot_pair_exists(figure_paths$sankey_pairwise_combined)) {
    append_fixed_o2_transition_figure_index(figure_paths, index_row)
    message("Skipping pairwise fixed-O2 mode Sankey plot; existing outputs found: ", figure_paths$sankey_pairwise_combined)
    return(invisible(figure_paths))
  }

  seed_paths <- sankey_tables$seed_paths
  state_cols <- sankey_tables$state_cols
  ploidy_cols <- sankey_tables$ploidy_cols
  o2_order <- sankey_tables$o2_order
  reference_state_col <- state_cols[[1L]]
  reference_ploidy_col <- ploidy_cols[[1L]]
  target_indices <- seq.int(2L, length(state_cols))
  mode_colors <- fixed_o2_reference_comparison_fill_colors()
  node_gap <- max(10, 0.028 * nrow(seed_paths))
  node_width <- 0.20
  curve_t <- seq(0, 1, length.out = 28)
  smooth_t <- curve_t * curve_t * (3 - 2 * curve_t)

  pair_label_for_o2 <- function(target_o2) {
    paste0(
      "O[2]==",
      format_o2_value(reference_o2),
      "~'vs'~O[2]==",
      format_o2_value(target_o2)
    )
  }

  node_rows <- list()
  flow_rows <- list()
  allocation_rows <- list()
  ribbon_rows <- list()
  axis_label_rows <- list()

  for (target_idx in target_indices) {
    target_o2 <- o2_order[[target_idx]]
    pair_label <- pair_label_for_o2(target_o2)
    target_state_col <- state_cols[[target_idx]]
    target_ploidy_col <- ploidy_cols[[target_idx]]

    panel_nodes <- list()
    for (axis_idx in c(1L, 2L)) {
      y0 <- 0
      state_levels <- fixed_o2_reference_comparison_state_order(is_reference_axis = axis_idx == 1L)
      state_col <- if (axis_idx == 1L) reference_state_col else target_state_col
      ploidy_col <- if (axis_idx == 1L) reference_ploidy_col else target_ploidy_col
      for (state in state_levels) {
        seed_idx <- seed_paths[[state_col]] == state
        n_state <- sum(seed_idx, na.rm = TRUE)
        mean_ploidy <- mean(suppressWarnings(as.numeric(seed_paths[[ploidy_col]][seed_idx])), na.rm = TRUE)
        reference_mean_ploidy <- mean(suppressWarnings(as.numeric(seed_paths[[reference_ploidy_col]][seed_idx])), na.rm = TRUE)
        if (!is.finite(mean_ploidy)) mean_ploidy <- NA_real_
        if (!is.finite(reference_mean_ploidy)) reference_mean_ploidy <- NA_real_
        panel_nodes[[length(panel_nodes) + 1L]] <- data.frame(
          pair_label = pair_label,
          target_o2 = target_o2,
          axis_idx = axis_idx,
          comparison_state = state,
          comparison_display = fixed_o2_reference_comparison_display(state),
          mode_display = fixed_o2_reference_comparison_current_mode(state),
          fill_key = fixed_o2_reference_comparison_fill_key(state),
          is_mode_transition = axis_idx > 1L && state %in% c("mode1_to_mode2", "mode2_to_mode1"),
          n_seeds = n_state,
          mean_ploidy = mean_ploidy,
          mean_ploidy_at_reference_o2_same_seeds = reference_mean_ploidy,
          ymin = y0,
          ymax = y0 + n_state,
          stringsAsFactors = FALSE
        )
        y0 <- y0 + n_state + node_gap
      }
    }
    panel_nodes <- rbind_fill_plain(panel_nodes)
    panel_nodes <- panel_nodes[panel_nodes$n_seeds > 0, , drop = FALSE]
    node_rows[[length(node_rows) + 1L]] <- panel_nodes

    target_state_levels <- fixed_o2_reference_comparison_state_order(is_reference_axis = FALSE)
    panel_flows <- lapply(target_state_levels, function(target_state) {
      seed_idx <- seed_paths[[target_state_col]] == target_state
      n_flow <- sum(seed_idx, na.rm = TRUE)
      if (!n_flow) return(NULL)
      reference_state <- unique(seed_paths[[reference_state_col]][seed_idx])
      reference_state <- reference_state[reference_state %in% fixed_o2_reference_comparison_state_order(is_reference_axis = TRUE)]
      if (!length(reference_state)) return(NULL)
      data.frame(
        pair_label = pair_label,
        target_o2 = target_o2,
        flow_id = paste(pair_label, target_state, sep = "_"),
        reference_state = reference_state[[1L]],
        target_state = target_state,
        reference_mode_display = fixed_o2_reference_comparison_current_mode(reference_state[[1L]]),
        fill_key = fixed_o2_reference_comparison_fill_key(target_state),
        n_seeds = n_flow,
        stringsAsFactors = FALSE
      )
    })
    panel_flows <- rbind_fill_plain(panel_flows)
    flow_rows[[length(flow_rows) + 1L]] <- panel_flows

    for (reference_state in fixed_o2_reference_comparison_state_order(is_reference_axis = TRUE)) {
      node <- panel_nodes[panel_nodes$axis_idx == 1L & panel_nodes$comparison_state == reference_state, , drop = FALSE]
      if (!nrow(node)) next
      ref_flows <- panel_flows[panel_flows$reference_state == reference_state, , drop = FALSE]
      if (!nrow(ref_flows)) next
      ref_flows <- ref_flows[order(match(ref_flows$target_state, target_state_levels)), , drop = FALSE]
      ymin <- node$ymin[[1L]]
      for (i in seq_len(nrow(ref_flows))) {
        allocation_rows[[length(allocation_rows) + 1L]] <- data.frame(
          pair_label = pair_label,
          flow_id = ref_flows$flow_id[[i]],
          axis_idx = 1L,
          ymin = ymin,
          ymax = ymin + ref_flows$n_seeds[[i]],
          stringsAsFactors = FALSE
        )
        ymin <- ymin + ref_flows$n_seeds[[i]]
      }
    }
    for (target_state in target_state_levels) {
      node <- panel_nodes[panel_nodes$axis_idx == 2L & panel_nodes$comparison_state == target_state, , drop = FALSE]
      flow <- panel_flows[panel_flows$target_state == target_state, , drop = FALSE]
      if (!nrow(node) || !nrow(flow)) next
      allocation_rows[[length(allocation_rows) + 1L]] <- data.frame(
        pair_label = pair_label,
        flow_id = flow$flow_id[[1L]],
        axis_idx = 2L,
        ymin = node$ymin[[1L]],
        ymax = node$ymax[[1L]],
        stringsAsFactors = FALSE
      )
    }

    axis_label_rows[[length(axis_label_rows) + 1L]] <- data.frame(
      pair_label = pair_label,
      x = c(1, 2),
      y = -1.55 * node_gap,
      label = c(
        paste0("O[2]==", format_o2_value(reference_o2)),
        paste0("O[2]==", format_o2_value(target_o2))
      ),
      stringsAsFactors = FALSE
    )
  }

  nodes <- rbind_fill_plain(node_rows)
  flows <- rbind_fill_plain(flow_rows)
  allocation <- rbind_fill_plain(allocation_rows)
  axis_labels <- rbind_fill_plain(axis_label_rows)
  if (!nrow(nodes) || !nrow(flows) || !nrow(allocation)) {
    message("Skipping pairwise fixed-O2 mode Sankey plot; no drawable pairwise flows were available.")
    return(invisible(NULL))
  }

  nodes$xmin <- nodes$axis_idx - node_width
  nodes$xmax <- nodes$axis_idx + node_width
  nodes$comparison_display_short <- gsub("Mode", "M", nodes$comparison_display, fixed = TRUE)
  nodes$mean_ploidy_label <- ifelse(is.finite(nodes$mean_ploidy), sprintf("%.2f", nodes$mean_ploidy), "NA")
  nodes$reference_mean_ploidy_label <- ifelse(
    is.finite(nodes$mean_ploidy_at_reference_o2_same_seeds),
    sprintf("%.2f", nodes$mean_ploidy_at_reference_o2_same_seeds),
    "NA"
  )
  base_node_label <- paste0("n=", nodes$n_seeds, "\nMean=", nodes$mean_ploidy_label)
  nodes$label <- ifelse(
    nodes$is_mode_transition,
    paste0(base_node_label, "\nRef=", nodes$reference_mean_ploidy_label),
    base_node_label
  )
  nodes$text_size <- ifelse(nodes$n_seeds < 35, 1.25, ifelse(nodes$n_seeds < 80, 1.45, 1.65))
  nodes$label_inside <- ifelse(nodes$n_seeds < 35, "", nodes$label)
  small_nodes <- nodes[nodes$n_seeds < 35, , drop = FALSE]
  if (nrow(small_nodes)) {
    small_nodes$label_x <- ifelse(small_nodes$axis_idx == 2L, small_nodes$xmax + 0.06, small_nodes$xmin - 0.06)
    small_nodes$segment_x <- ifelse(small_nodes$axis_idx == 2L, small_nodes$xmax, small_nodes$xmin)
    small_nodes$hjust <- ifelse(small_nodes$axis_idx == 2L, 0, 1)
  }

  for (i in seq_len(nrow(flows))) {
    flow <- flows[i, , drop = FALSE]
    left <- allocation[allocation$flow_id == flow$flow_id[[1L]] & allocation$axis_idx == 1L, , drop = FALSE]
    right <- allocation[allocation$flow_id == flow$flow_id[[1L]] & allocation$axis_idx == 2L, , drop = FALSE]
    if (!nrow(left) || !nrow(right)) next
    x_left <- 1 + node_width
    x_right <- 2 - node_width
    x <- x_left + (x_right - x_left) * curve_t
    top <- left$ymax[[1L]] + (right$ymax[[1L]] - left$ymax[[1L]]) * smooth_t
    bottom <- left$ymin[[1L]] + (right$ymin[[1L]] - left$ymin[[1L]]) * smooth_t
    ribbon_rows[[length(ribbon_rows) + 1L]] <- data.frame(
      pair_label = flow$pair_label[[1L]],
      segment_id = flow$flow_id[[1L]],
      x = c(x, rev(x)),
      y = c(top, rev(bottom)),
      fill_key = flow$fill_key[[1L]],
      stringsAsFactors = FALSE
    )
  }
  ribbons <- rbind_fill_plain(ribbon_rows)
  if (!nrow(ribbons)) {
    message("Skipping pairwise fixed-O2 mode Sankey plot; no drawable ribbons were available.")
    return(invisible(NULL))
  }

  pair_levels <- unique(vapply(o2_order[target_indices], pair_label_for_o2, character(1)))
  nodes$pair_label <- factor(nodes$pair_label, levels = pair_levels)
  ribbons$pair_label <- factor(ribbons$pair_label, levels = pair_levels)
  axis_labels$pair_label <- factor(axis_labels$pair_label, levels = pair_levels)
  if (nrow(small_nodes)) {
    small_nodes$pair_label <- factor(small_nodes$pair_label, levels = pair_levels)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = ribbons,
      ggplot2::aes(x = x, y = y, group = segment_id, fill = fill_key),
      color = NA,
      alpha = 0.42
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_key),
      color = "white",
      linewidth = 0.28
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = axis_idx, y = (ymin + ymax) / 2, label = label_inside, size = text_size),
      color = "white",
      lineheight = 0.78,
      fontface = "bold"
    ) +
    ggplot2::geom_segment(
      data = small_nodes,
      ggplot2::aes(
        x = segment_x,
        xend = label_x,
        y = (ymin + ymax) / 2,
        yend = (ymin + ymax) / 2
      ),
      inherit.aes = FALSE,
      color = "#111827",
      linewidth = 0.25
    ) +
    ggplot2::geom_label(
      data = small_nodes,
      ggplot2::aes(x = label_x, y = (ymin + ymax) / 2, label = label, hjust = hjust),
      inherit.aes = FALSE,
      color = "#111827",
      fill = "white",
      linewidth = 0.15,
      label.padding = grid::unit(0.08, "lines"),
      size = 1.7,
      lineheight = 0.78,
      fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = axis_labels,
      ggplot2::aes(x = x, y = y, label = label),
      parse = TRUE,
      inherit.aes = FALSE,
      size = 2.4,
      color = "#3b4652"
    ) +
    ggplot2::facet_wrap(
      ~ pair_label,
      nrow = 1,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(breaks = NULL, limits = c(0.68, 2.58), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.04))) +
    ggplot2::scale_fill_manual(values = mode_colors, drop = FALSE) +
    ggplot2::labs(
      title = bquote("Pairwise Mode-label changes from reference " * O[2] == .(format_o2_value(reference_o2))),
      subtitle = bquote("Each panel compares " * O[2] == .(format_o2_value(reference_o2)) * " directly with one target " * O[2] * "; green marks Mode1->Mode2 and purple marks Mode2->Mode1."),
      x = NULL,
      y = NULL,
      fill = "Mode change"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#eef3f8", color = "#d8e0e8"),
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )
  save_plot_pair(p, figure_paths$sankey_pairwise_combined, width = 15.5, height = 4.9)
  append_fixed_o2_transition_figure_index(figure_paths, index_row)
  invisible(figure_paths)
}

write_fixed_o2_reference_grouped_distribution_plots <- function(plot_df, paths, base_theme, overwrite = FALSE) {
  reference_o2 <- sort(unique(plot_df$fixed_o2[is.finite(plot_df$fixed_o2)]))
  if (!length(reference_o2)) return(invisible(NULL))
  dir.create(paths$reference_grouped_dir, recursive = TRUE, showWarnings = FALSE)
  index_rows <- lapply(reference_o2, function(ref_o2) {
    prefix <- fixed_o2_reference_grouped_prefix(paths, ref_o2)
    ref_rows <- plot_df[is.finite(plot_df$fixed_o2) & abs(plot_df$fixed_o2 - ref_o2) < 1e-9, , drop = FALSE]
    ref_rows <- ref_rows[ref_rows$mode_label %in% c("mode1", "mode2"), c("seed_id", "mode_label"), drop = FALSE]
    ref_rows <- ref_rows[!duplicated(ref_rows$seed_id), , drop = FALSE]
    if (!nrow(ref_rows)) {
      message("Skipping reference-grouped fixed-O2 plot for reference O2=", format_o2_value(ref_o2), "; no mode1/mode2 reference labels.")
      return(NULL)
    }
    idx <- match(plot_df$seed_id, ref_rows$seed_id)
    ref_mode <- ref_rows$mode_label[idx]
    ref_df <- plot_df[!is.na(idx) & ref_mode %in% c("mode1", "mode2"), , drop = FALSE]
    ref_df$reference_mode_label <- ref_mode[!is.na(idx) & ref_mode %in% c("mode1", "mode2")]
    ref_df$reference_mode_label_display <- ifelse(
      ref_df$reference_mode_label == "mode1",
      "Mode1",
      ifelse(ref_df$reference_mode_label == "mode2", "Mode2", NA_character_)
    )
    ref_df$reference_mode_label_display <- factor(ref_df$reference_mode_label_display, levels = c("Mode1", "Mode2"))
    ref_df <- ref_df[!is.na(ref_df$reference_mode_label_display), , drop = FALSE]
    if (!nrow(ref_df)) {
      message("Skipping reference-grouped fixed-O2 plot for reference O2=", format_o2_value(ref_o2), "; no joined analytical-solution rows.")
      return(NULL)
    }

    ref_df$split_group <- (ref_df$fixed_o2_index - 1L) * 2L + ifelse(ref_df$reference_mode_label_display == "Mode1", 1L, 2L)
    median_df <- stats::aggregate(
      analytical_solution ~ fixed_o2_index + fixed_o2_label + reference_mode_label_display,
      data = ref_df,
      FUN = stats::median
    )
    names(median_df)[names(median_df) == "analytical_solution"] <- "median_analytical_solution"
    median_df$x_start <- median_df$fixed_o2_index + ifelse(median_df$reference_mode_label_display == "Mode1", -0.32, 0.06)
    median_df$x_end <- median_df$fixed_o2_index + ifelse(median_df$reference_mode_label_display == "Mode1", -0.06, 0.32)
    o2_levels <- levels(ref_df$fixed_o2_label)
    p_ref <- ggplot2::ggplot(
      ref_df,
      ggplot2::aes(
        x = fixed_o2_index,
        y = analytical_solution,
        fill = reference_mode_label_display,
        group = split_group
      )
    ) +
      geom_split_violin(color = "#38485a", linewidth = 0.25, trim = FALSE, scale = "width", adjust = 0.75, width = 0.86, alpha = 0.95) +
      ggplot2::geom_hline(yintercept = 2, linetype = "dashed", color = "#4b5563", linewidth = 0.45) +
      ggplot2::geom_segment(
        data = median_df,
        ggplot2::aes(
          x = x_start,
          xend = x_end,
          y = median_analytical_solution,
          yend = median_analytical_solution,
          color = reference_mode_label_display
        ),
        inherit.aes = FALSE,
        linewidth = 1.15,
        lineend = "round"
      ) +
      ggplot2::geom_point(
        data = median_df,
        ggplot2::aes(x = fixed_o2_index, y = median_analytical_solution, color = reference_mode_label_display),
        inherit.aes = FALSE,
        size = 1.7,
        stroke = 0.2
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(o2_levels),
        labels = o2_discrete_label_expressions(o2_levels)
      ) +
      ggplot2::scale_fill_manual(values = c(Mode1 = "#1f78b4", Mode2 = "#ff7f00"), drop = FALSE) +
      ggplot2::scale_color_manual(values = c(Mode1 = "#0d4f8b", Mode2 = "#b85b00"), guide = "none", drop = FALSE) +
      ggplot2::labs(
        title = bquote("Fixed-" * O[2] * " analytical solution grouped by reference " * O[2] == .(format_o2_value(ref_o2))),
        x = expression("Fixed " * O[2] * " (%)"),
        y = "Analytical solution (dominant mean ploidy)",
        fill = "Reference Mode"
      ) +
      base_theme +
      ggplot2::theme(
        legend.position = c(0.985, 0.985),
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = ggplot2::element_rect(fill = "white", color = "#d8e0e8", linewidth = 0.25),
        legend.margin = ggplot2::margin(3, 5, 3, 5),
        legend.key.size = grid::unit(0.36, "cm")
      )
    if (overwrite || !plot_pair_exists(prefix)) {
      save_plot_pair(p_ref, prefix, width = 6.2, height = 4.8)
    } else {
      message("Skipping reference-grouped fixed-O2 plot; existing outputs found: ", prefix)
    }

    mode_counts <- table(factor(ref_rows$mode_label, levels = c("mode1", "mode2")))
    data.frame(
      reference_o2 = ref_o2,
      reference_o2_label = mode_reference_o2_label(ref_o2),
      n_reference_mode1 = unname(mode_counts[["mode1"]]),
      n_reference_mode2 = unname(mode_counts[["mode2"]]),
      pdf = paste0(prefix, ".pdf"),
      png = paste0(prefix, ".png"),
      stringsAsFactors = FALSE
    )
  })
  index <- rbind_fill_plain(index_rows)
  if (nrow(index)) {
    write_csv(index, file.path(paths$reference_grouped_dir, "reference_o2_grouped_plot_index.csv"))
  }
  invisible(index)
}

write_fixed_o2_attractor_distribution_plots <- function(mode_tables_dir, overwrite = FALSE) {
  if (!dir.exists(mode_tables_dir)) {
    message("Skipping fixed-O2 attractor distribution plots; mode table directory does not exist: ", mode_tables_dir)
    return(invisible(NULL))
  }
  paths <- fixed_o2_attractor_distribution_paths(mode_tables_dir)
  if (!overwrite && fixed_o2_attractor_distribution_complete(mode_tables_dir)) {
    message("Skipping fixed-O2 attractor distribution plots; existing outputs found under: ", paths$figure_dir)
    return(invisible(paths))
  }
  plot_df <- prepare_fixed_o2_attractor_distribution_data(mode_tables_dir)
  write_fixed_o2_reference_transition_tables(plot_df, mode_tables_dir, reference_o2 = 5, overwrite = overwrite)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping fixed-O2 attractor distribution plots; ggplot2 is not available.")
    return(invisible(NULL))
  }
  dir.create(paths$figure_dir, recursive = TRUE, showWarnings = FALSE)
  write_fixed_o2_reference_transition_plots(mode_tables_dir, reference_o2 = 5, overwrite = overwrite)
  write_fixed_o2_reference_mode_sankey(plot_df, mode_tables_dir, reference_o2 = 5, overwrite = overwrite)
  write_fixed_o2_reference_pairwise_sankey_combined(plot_df, mode_tables_dir, reference_o2 = 5, overwrite = overwrite)

  base_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )

  p_all <- ggplot2::ggplot(plot_df, ggplot2::aes(x = fixed_o2_label, y = analytical_solution)) +
    ggplot2::geom_violin(fill = "#d9dee7", color = "#657789", linewidth = 0.3, trim = TRUE, scale = "width", adjust = 0.55) +
    ggplot2::geom_boxplot(width = 0.16, fill = "white", color = "#1d2a35", linewidth = 0.35, outlier.size = 0.45, outlier.alpha = 0.35) +
    ggplot2::geom_hline(yintercept = 2, linetype = "dashed", color = "#4b5563", linewidth = 0.45) +
    ggplot2::scale_x_discrete(labels = o2_discrete_label_expressions) +
    ggplot2::labs(
      title = bquote("Fixed-" * O[2] * " analytical solution across all seeds"),
      x = expression("Fixed " * O[2] * " (%)"),
      y = "Analytical solution (dominant mean ploidy)"
    ) +
    base_theme
  if (overwrite || !plot_pair_exists(paths$all_seeds_prefix)) {
    save_plot_pair(p_all, paths$all_seeds_prefix, width = 6.2, height = 4.8)
  } else {
    message("Skipping fixed-O2 all-seed plot; existing outputs found: ", paths$all_seeds_prefix)
  }

  mode_df <- plot_df[!is.na(plot_df$mode_label_display), , drop = FALSE]
  if (nrow(mode_df)) {
    mode_df$split_group <- (mode_df$fixed_o2_index - 1L) * 2L + ifelse(mode_df$mode_label_display == "Mode1", 1L, 2L)
    median_df <- stats::aggregate(
      analytical_solution ~ fixed_o2_index + fixed_o2_label + mode_label_display,
      data = mode_df,
      FUN = stats::median
    )
    names(median_df)[names(median_df) == "analytical_solution"] <- "median_analytical_solution"
    median_df$x_start <- median_df$fixed_o2_index + ifelse(median_df$mode_label_display == "Mode1", -0.32, 0.06)
    median_df$x_end <- median_df$fixed_o2_index + ifelse(median_df$mode_label_display == "Mode1", -0.06, 0.32)
    o2_levels <- levels(mode_df$fixed_o2_label)
    p_mode <- ggplot2::ggplot(
      mode_df,
      ggplot2::aes(x = fixed_o2_index, y = analytical_solution, fill = mode_label_display, group = split_group)
    ) +
      geom_split_violin(color = "#38485a", linewidth = 0.25, trim = TRUE, scale = "width", adjust = 0.75, width = 0.86, alpha = 0.95) +
      ggplot2::geom_hline(yintercept = 2, linetype = "dashed", color = "#4b5563", linewidth = 0.45) +
      ggplot2::geom_segment(
        data = median_df,
        ggplot2::aes(
          x = x_start,
          xend = x_end,
          y = median_analytical_solution,
          yend = median_analytical_solution,
          color = mode_label_display
        ),
        inherit.aes = FALSE,
        linewidth = 1.15,
        lineend = "round"
      ) +
      ggplot2::geom_point(
        data = median_df,
        ggplot2::aes(x = fixed_o2_index, y = median_analytical_solution, color = mode_label_display),
        inherit.aes = FALSE,
        size = 1.7,
        stroke = 0.2
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(o2_levels),
        labels = o2_discrete_label_expressions(o2_levels)
      ) +
      ggplot2::scale_fill_manual(values = c(Mode1 = "#1f78b4", Mode2 = "#ff7f00"), drop = FALSE) +
      ggplot2::scale_color_manual(values = c(Mode1 = "#0d4f8b", Mode2 = "#b85b00"), guide = "none", drop = FALSE) +
      ggplot2::labs(
        title = bquote("Fixed-" * O[2] * " analytical solution by Mode"),
        x = expression("Fixed " * O[2] * " (%)"),
        y = "Analytical solution (dominant mean ploidy)",
        fill = "Mode"
      ) +
      base_theme +
      ggplot2::theme(
        legend.position = c(0.985, 0.985),
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = ggplot2::element_rect(fill = "white", color = "#d8e0e8", linewidth = 0.25),
        legend.margin = ggplot2::margin(3, 5, 3, 5),
        legend.key.size = grid::unit(0.36, "cm")
      )
    if (overwrite || !plot_pair_exists(paths$by_mode_prefix)) {
      save_plot_pair(p_mode, paths$by_mode_prefix, width = 6.2, height = 4.8)
    } else {
      message("Skipping fixed-O2 row-wise mode plot; existing outputs found: ", paths$by_mode_prefix)
    }
  }
  write_fixed_o2_reference_grouped_distribution_plots(plot_df, paths, base_theme, overwrite = overwrite)
  invisible(paths)
}

write_mode_contribution_global_summary <- function(output_dir, reference_o2) {
  runs <- lapply(reference_o2, function(o2) read_reference_contribution(output_dir, o2))
  index <- do.call(rbind, lapply(runs, `[[`, "summary"))
  top <- do.call(rbind, lapply(seq_along(runs), function(i) {
    df <- runs[[i]]$top_features
    if (!nrow(df)) return(NULL)
    df$mode_reference_o2 <- reference_o2[[i]]
    df
  }))
  index_path <- file.path(output_dir, "mode_parameter_contribution_index.csv")
  write_csv(index, index_path)
  write_csv(top, file.path(output_dir, "mode_parameter_top_features_across_reference_o2.csv"))
  write_o2_group_auc_summary(output_dir, reference_o2, index)
  write_cumulative_auc_global_summary(output_dir, reference_o2)
  invisible(index_path)
}

write_dominant_ploidy_contribution_global_summary <- function(output_dir, reference_o2) {
  runs <- lapply(reference_o2, function(o2) read_dominant_ploidy_contribution(output_dir, o2))
  index <- do.call(rbind, lapply(runs, `[[`, "summary"))
  top <- do.call(rbind, lapply(seq_along(runs), function(i) {
    df <- runs[[i]]$top_features
    if (!nrow(df)) return(NULL)
    df$mode_reference_o2 <- reference_o2[[i]]
    df
  }))
  index_path <- file.path(output_dir, "dominant_ploidy_parameter_contribution_index.csv")
  write_csv(index, index_path)
  write_csv(top, file.path(output_dir, "dominant_ploidy_parameter_top_features_across_reference_o2.csv"))

  plots_dir <- file.path(output_dir, "summary_plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(index)) {
    index$mode_reference_o2 <- suppressWarnings(as.numeric(index$mode_reference_o2))
    index$mode_reference_o2_label <- factor(
      mode_reference_o2_label(index$mode_reference_o2),
      levels = mode_reference_o2_label(reference_o2)
    )
    plot_df <- index[is.finite(index$mode_reference_o2), , drop = FALSE]
    if (nrow(plot_df) && all(c("top3_joint_cv_r2", "top3_joint_apparent_r2") %in% names(plot_df))) {
      r2_long <- rbind_fill_plain(list(
        data.frame(
          mode_reference_o2_label = plot_df$mode_reference_o2_label,
          r2 = suppressWarnings(as.numeric(plot_df$top3_joint_cv_r2)),
          evaluation = "Cross-validated",
          stringsAsFactors = FALSE
        ),
        data.frame(
          mode_reference_o2_label = plot_df$mode_reference_o2_label,
          r2 = suppressWarnings(as.numeric(plot_df$top3_joint_apparent_r2)),
          evaluation = "Apparent",
          stringsAsFactors = FALSE
        )
      ))
      r2_long <- r2_long[is.finite(r2_long$r2), , drop = FALSE]
      if (nrow(r2_long)) {
        p <- ggplot2::ggplot(r2_long, ggplot2::aes(x = mode_reference_o2_label, y = r2, color = evaluation, group = evaluation))
        if (length(unique(r2_long$mode_reference_o2_label)) > 1L) {
          p <- p + ggplot2::geom_line(linewidth = 0.85)
        }
        p <- p +
          ggplot2::geom_point(size = 1.8) +
          ggplot2::scale_color_manual(values = c("Cross-validated" = "#1f78b4", "Apparent" = "grey45")) +
          ggplot2::labs(
            title = expression("Top3 joint regression " * R^2 * " across reference " * O[2]),
            x = expression("Reference fixed " * O[2] * " (%)"),
            y = expression(R^2),
            color = NULL
          ) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
        save_plot_pair(p, file.path(plots_dir, "top3_joint_r2_across_reference_o2"), width = 6.4, height = 4.2)
      }
    }
  }
  invisible(index_path)
}

run_reference_contribution <- function(best_df,
                                       best_csv,
                                       mode_tables_dir,
                                       mode_reference_o2,
                                       output_dir,
                                       n_bootstrap,
                                       sample_fraction,
                                       alpha,
                                       random_seed,
                                       top_n_interactions,
                                       tree_depth,
                                       tree_min_node,
                                       cumulative_auc_top_n) {
  params <- umap_parameter_set("invivo")
  log10_params <- umap_log10_parameter_set("invivo")
  mode_tab <- read_fixed_o2_reference_mode_table(
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2
  )
  joined <- add_seed_id_if_needed(best_df)
  idx <- match(joined$seed_id, mode_tab$seed_id)
  joined <- joined[!is.na(idx), , drop = FALSE]
  mode_tab <- mode_tab[idx[!is.na(idx)], , drop = FALSE]
  keep_mode <- mode_tab$mode_label %in% c("mode1", "mode2")
  joined <- joined[keep_mode, , drop = FALSE]
  mode_tab <- mode_tab[keep_mode, , drop = FALSE]
  if (!nrow(joined)) stop("No seeds with mode1/mode2 labels for reference O2=", mode_reference_o2)
  y <- ifelse(mode_tab$mode_label == "mode1", 1L, 0L)
  if (length(unique(y)) < 2L) {
    stop("Only one mode is present for reference O2=", mode_reference_o2, "; contribution model requires both modes.")
  }

  transformed <- transform_umap_features(joined, params, log10_params)
  complete <- stats::complete.cases(transformed) & !is.na(y)
  if (!all(complete)) {
    joined <- joined[complete, , drop = FALSE]
    mode_tab <- mode_tab[complete, , drop = FALSE]
    transformed <- transformed[complete, , drop = FALSE]
    y <- y[complete]
  }
  if (length(unique(y)) < 2L) {
    stop("Only one mode remains after complete-case filtering for reference O2=", mode_reference_o2, ".")
  }
  z_main <- standardize_matrix(transformed)
  mode_df <- data.frame(
    seed = joined$seed,
    seed_id = joined$seed_id,
    mode_label = mode_tab$mode_label,
    trajectory_regime = mode_tab$trajectory_regime,
    mode_reference_o2 = mode_reference_o2,
    mode_reference_dominant_mean_ploidy = numeric_column_or_na(mode_tab, "mode_reference_dominant_mean_ploidy"),
    stringsAsFactors = FALSE
  )

  ref_dir <- reference_contribution_dir(output_dir, mode_reference_o2)
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

  main_tab <- main_effect_table(joined, transformed, z_main, y, params, log10_params)
  pair_tab <- interaction_table(z_main, y)
  glmnet_tab <- fit_glmnet_importance(
    z_main = z_main,
    y = y,
    n_bootstrap = n_bootstrap,
    sample_fraction = sample_fraction,
    alpha = alpha,
    seed = random_seed
  )
  glmnet_status <- attr(glmnet_tab, "status")

  combined <- rbind_fill_plain(list(
    main_tab[, intersect(names(main_tab), c(
      "feature_name", "feature_type", "parameter", "parameter_a", "parameter_b",
      "auc_contribution", "wilcox_p", "wilcox_fdr", "glm_coef_mode1", "glm_p", "glm_fdr"
    )), drop = FALSE],
    pair_tab[, intersect(names(pair_tab), c(
      "feature_name", "feature_type", "parameter", "parameter_a", "parameter_b",
      "interaction_auc_contribution", "interaction_glm_coef_mode1", "interaction_glm_p",
      "interaction_glm_fdr", "interaction_model_auc_contribution"
    )), drop = FALSE]
  ))
  combined <- merge(
    combined,
    glmnet_tab[, c("feature_name", "elastic_net_coef", "selection_frequency"), drop = FALSE],
    by = "feature_name",
    all.x = TRUE,
    sort = FALSE
  )
  combined$contribution_score <- ifelse(
    is.finite(combined$selection_frequency),
    combined$selection_frequency,
    pmax(
      suppressWarnings(as.numeric(combined$auc_contribution)),
      suppressWarnings(as.numeric(combined$interaction_model_auc_contribution)),
      suppressWarnings(as.numeric(combined$interaction_auc_contribution)),
      na.rm = TRUE
    )
  )
  combined$contribution_score[!is.finite(combined$contribution_score)] <- NA_real_
  combined$abs_elastic_net_coef <- abs(suppressWarnings(as.numeric(combined$elastic_net_coef)))
  combined <- combined[order(
    -ifelse(is.na(combined$selection_frequency), -Inf, combined$selection_frequency),
    -ifelse(is.na(combined$abs_elastic_net_coef), -Inf, combined$abs_elastic_net_coef),
    -ifelse(is.na(combined$contribution_score), -Inf, combined$contribution_score),
    combined$feature_name
  ), , drop = FALSE]
  combined$rank <- seq_len(nrow(combined))
  auc_tab <- write_auc_outputs(
    combined = combined,
    z_main = z_main,
    y = y,
    mode_reference_o2 = mode_reference_o2,
    ref_dir = ref_dir
  )
  top3_joint <- fit_top3_joint_auc(
    combined = combined,
    z_main = z_main,
    y = y,
    mode_reference_o2 = mode_reference_o2,
    output_dir = file.path(ref_dir, "top3_joint"),
    cv_k = 10L,
    seed = random_seed
  )
  cumulative_auc <- fit_cumulative_feature_auc(
    combined = combined,
    z_main = z_main,
    y = y,
    mode_reference_o2 = mode_reference_o2,
    output_dir = file.path(ref_dir, "cumulative_auc"),
    top_n = cumulative_auc_top_n,
    cv_k = 10L,
    seed = random_seed
  )
  top3_fields <- top_feature_summary_fields(combined, auc_tab, n = 3L)

  rule_tab <- build_rule_tree(
    x_df = transformed,
    y = y,
    log10_params = log10_params,
    max_depth = tree_depth,
    min_node = tree_min_node
  )

  top_interactions <- combined[combined$feature_type == "interaction", , drop = FALSE]
  top_interactions <- top_interactions[!is.na(top_interactions$parameter_a) & !is.na(top_interactions$parameter_b), , drop = FALSE]
  phase_paths <- character()
  if (nrow(top_interactions)) {
    phase_paths <- plot_top_interactions(
      top_interactions = top_interactions,
      transformed_df = transformed,
      mode_df = mode_df,
      output_dir = file.path(ref_dir, "phase_plots"),
      log10_params = log10_params,
      top_n = top_n_interactions
    )
  }

  summary_base <- data.frame(
    mode_reference_o2 = mode_reference_o2,
    n_seed = length(y),
    n_mode1 = sum(y == 1L),
    n_mode2 = sum(y == 0L),
    best_csv = normalizePath(best_csv, mustWork = FALSE),
    mode_tables_dir = normalizePath(mode_tables_dir, mustWork = FALSE),
    output_dir = normalizePath(ref_dir, mustWork = FALSE),
    top_main_feature = if (nrow(main_tab)) main_tab$feature_name[[1L]] else NA_character_,
    top_interaction_feature = if (nrow(pair_tab)) pair_tab$feature_name[[1L]] else NA_character_,
    top_combined_feature = if (nrow(combined)) combined$feature_name[[1L]] else NA_character_,
    cumulative_auc_top_n = if (!is.null(cumulative_auc)) nrow(cumulative_auc) else NA_integer_,
    n_phase_plot_files = length(phase_paths),
    stringsAsFactors = FALSE
  )
  top3_joint_fields <- data.frame(
    top3_joint_fit_status = if (!is.null(top3_joint) && "fit_status" %in% names(top3_joint)) top3_joint$fit_status[[1L]] else NA_character_,
    top3_joint_apparent_auc_mode1_higher = if (!is.null(top3_joint) && "apparent_auc_mode1_higher" %in% names(top3_joint)) top3_joint$apparent_auc_mode1_higher[[1L]] else NA_real_,
    top3_joint_apparent_auc = if (!is.null(top3_joint) && "apparent_auc" %in% names(top3_joint)) top3_joint$apparent_auc[[1L]] else NA_real_,
    top3_joint_cv_auc_mode1_higher = if (!is.null(top3_joint) && "cv_auc_mode1_higher" %in% names(top3_joint)) top3_joint$cv_auc_mode1_higher[[1L]] else NA_real_,
    top3_joint_cv_auc = if (!is.null(top3_joint) && "cv_auc" %in% names(top3_joint)) top3_joint$cv_auc[[1L]] else NA_real_,
    top3_joint_cv_k = if (!is.null(top3_joint) && "cv_k" %in% names(top3_joint)) top3_joint$cv_k[[1L]] else NA_integer_,
    top3_joint_cv_success_folds = if (!is.null(top3_joint) && "cv_success_folds" %in% names(top3_joint)) top3_joint$cv_success_folds[[1L]] else NA_integer_,
    stringsAsFactors = FALSE
  )
  summary <- cbind(summary_base, top3_fields, top3_joint_fields)

  write_csv(mode_df, file.path(ref_dir, "mode_parameter_seed_labels.csv"))
  write_csv(main_tab, file.path(ref_dir, "mode_parameter_main_effects.csv"))
  write_csv(pair_tab, file.path(ref_dir, "mode_parameter_pairwise_interactions.csv"))
  write_csv(glmnet_tab, file.path(ref_dir, "mode_parameter_stability_selection.csv"))
  write_tsv_plain(glmnet_status, file.path(ref_dir, "mode_parameter_elastic_net_status.tsv"))
  write_csv(combined, file.path(ref_dir, "mode_parameter_feature_importance.csv"))
  write_csv(rule_tab, file.path(ref_dir, "mode_parameter_decision_rules.csv"))
  write_tsv_plain(summary, file.path(ref_dir, "mode_parameter_contribution_summary.tsv"))

  list(
    summary = summary,
    top_features = head(combined, 30L),
    output_dir = ref_dir
  )
}

run_reference_dominant_ploidy_contribution <- function(best_df,
                                                       best_csv,
                                                       mode_tables_dir,
                                                       mode_reference_o2,
                                                       output_dir,
                                                       n_bootstrap,
                                                       sample_fraction,
                                                       alpha,
                                                       random_seed,
                                                       top_n_interactions,
                                                       cumulative_auc_top_n) {
  params <- umap_parameter_set("invivo")
  log10_params <- umap_log10_parameter_set("invivo")
  mode_tab <- read_fixed_o2_reference_mode_table(
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2
  )
  joined <- add_seed_id_if_needed(best_df)
  idx <- match(joined$seed_id, mode_tab$seed_id)
  joined <- joined[!is.na(idx), , drop = FALSE]
  mode_tab <- mode_tab[idx[!is.na(idx)], , drop = FALSE]
  y_raw <- numeric_column_or_na(mode_tab, "mode_reference_dominant_mean_ploidy")
  keep_response <- is.finite(y_raw)
  joined <- joined[keep_response, , drop = FALSE]
  mode_tab <- mode_tab[keep_response, , drop = FALSE]
  y_raw <- y_raw[keep_response]
  if (length(y_raw) < 3L) stop("Fewer than three finite dominant-ploidy values for reference O2=", mode_reference_o2)
  if (!is.finite(stats::sd(y_raw)) || stats::sd(y_raw) <= 0) {
    stop("Dominant-ploidy response has zero/non-finite variance for reference O2=", mode_reference_o2)
  }

  transformed <- transform_umap_features(joined, params, log10_params)
  complete <- stats::complete.cases(transformed) & is.finite(y_raw)
  if (!all(complete)) {
    joined <- joined[complete, , drop = FALSE]
    mode_tab <- mode_tab[complete, , drop = FALSE]
    transformed <- transformed[complete, , drop = FALSE]
    y_raw <- y_raw[complete]
  }
  if (length(y_raw) < 3L || !is.finite(stats::sd(y_raw)) || stats::sd(y_raw) <= 0) {
    stop("Dominant-ploidy response is not estimable after complete-case filtering for reference O2=", mode_reference_o2)
  }

  z_main <- standardize_matrix(transformed)
  y <- y_raw - 2
  seed_df <- data.frame(
    seed = joined$seed,
    seed_id = joined$seed_id,
    mode_label = mode_tab$mode_label,
    trajectory_regime = mode_tab$trajectory_regime,
    mode_reference_o2 = mode_reference_o2,
    mode_reference_dominant_mean_ploidy = y_raw,
    dominant_mean_ploidy_minus_2 = y,
    mode_reference_status = if ("mode_reference_status" %in% names(mode_tab)) mode_tab$mode_reference_status else NA_character_,
    mode_reference_dominant_growth_rate = numeric_column_or_na(mode_tab, "mode_reference_dominant_growth_rate"),
    mode_reference_spectral_gap = numeric_column_or_na(mode_tab, "mode_reference_spectral_gap"),
    stringsAsFactors = FALSE
  )

  ref_dir <- reference_contribution_dir(output_dir, mode_reference_o2)
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

  main_tab <- continuous_main_effect_table(joined, transformed, z_main, y, params, log10_params)
  pair_tab <- continuous_interaction_table(z_main, y)
  glmnet_tab <- fit_glmnet_regression_importance(
    z_main = z_main,
    y = y,
    n_bootstrap = n_bootstrap,
    sample_fraction = sample_fraction,
    alpha = alpha,
    seed = random_seed
  )
  glmnet_status <- attr(glmnet_tab, "status")

  combined <- rbind_fill_plain(list(
    main_tab[, intersect(names(main_tab), c(
      "feature_name", "feature_type", "parameter", "parameter_a", "parameter_b",
      "lm_coef_per_z", "lm_t", "lm_p", "lm_fdr", "lm_r2", "spearman_rho",
      "spearman_p", "spearman_fdr", "pearson_r", "pearson_p", "pearson_fdr",
      "direction"
    )), drop = FALSE],
    pair_tab[, intersect(names(pair_tab), c(
      "feature_name", "feature_type", "parameter", "parameter_a", "parameter_b",
      "interaction_lm_coef", "interaction_lm_t", "interaction_lm_p", "interaction_lm_fdr",
      "interaction_full_r2", "main_effects_r2", "interaction_delta_r2",
      "interaction_single_feature_r2", "interaction_spearman_rho",
      "interaction_spearman_p", "interaction_spearman_fdr", "direction"
    )), drop = FALSE]
  ))
  combined <- merge(
    combined,
    glmnet_tab[, c("feature_name", "elastic_net_coef", "selection_frequency"), drop = FALSE],
    by = "feature_name",
    all.x = TRUE,
    sort = FALSE
  )
  combined$contribution_score <- ifelse(
    is.finite(combined$selection_frequency),
    combined$selection_frequency,
    pmax(
      suppressWarnings(as.numeric(combined$lm_r2)),
      suppressWarnings(as.numeric(combined$interaction_delta_r2)),
      suppressWarnings(as.numeric(combined$interaction_single_feature_r2)),
      na.rm = TRUE
    )
  )
  combined$contribution_score[!is.finite(combined$contribution_score)] <- NA_real_
  combined$abs_elastic_net_coef <- abs(suppressWarnings(as.numeric(combined$elastic_net_coef)))
  combined <- combined[order(
    -ifelse(is.na(combined$selection_frequency), -Inf, combined$selection_frequency),
    -ifelse(is.na(combined$abs_elastic_net_coef), -Inf, combined$abs_elastic_net_coef),
    -ifelse(is.na(combined$contribution_score), -Inf, combined$contribution_score),
    combined$feature_name
  ), , drop = FALSE]
  combined$rank <- seq_len(nrow(combined))

  r2_tab <- write_continuous_r2_outputs(
    combined = combined,
    z_main = z_main,
    y = y_raw,
    mode_reference_o2 = mode_reference_o2,
    ref_dir = ref_dir
  )
  top3_joint <- fit_top3_joint_regression(
    combined = combined,
    z_main = z_main,
    y = y_raw,
    mode_reference_o2 = mode_reference_o2,
    output_dir = file.path(ref_dir, "top3_joint"),
    cv_k = 10L,
    seed = random_seed
  )
  cumulative_r2 <- fit_cumulative_feature_r2(
    combined = combined,
    z_main = z_main,
    y = y_raw,
    mode_reference_o2 = mode_reference_o2,
    output_dir = file.path(ref_dir, "cumulative_r2"),
    top_n = cumulative_auc_top_n,
    cv_k = 10L,
    seed = random_seed
  )
  top3_fields <- continuous_top_feature_summary_fields(combined, r2_tab, n = 3L)

  top_interactions <- combined[combined$feature_type == "interaction", , drop = FALSE]
  top_interactions <- top_interactions[!is.na(top_interactions$parameter_a) & !is.na(top_interactions$parameter_b), , drop = FALSE]
  phase_paths <- character()
  if (nrow(top_interactions)) {
    phase_paths <- plot_top_interactions_continuous(
      top_interactions = top_interactions,
      transformed_df = transformed,
      y = y_raw,
      output_dir = file.path(ref_dir, "phase_plots"),
      log10_params = log10_params,
      top_n = top_n_interactions
    )
  }

  summary_base <- data.frame(
    mode_reference_o2 = mode_reference_o2,
    n_seed = length(y_raw),
    response_target = "mode_reference_dominant_mean_ploidy",
    response_model_target = "dominant_mean_ploidy_minus_2",
    response_mean = mean(y_raw, na.rm = TRUE),
    response_sd = stats::sd(y_raw, na.rm = TRUE),
    response_min = min(y_raw, na.rm = TRUE),
    response_q25 = as.numeric(stats::quantile(y_raw, 0.25, na.rm = TRUE)),
    response_median = stats::median(y_raw, na.rm = TRUE),
    response_q75 = as.numeric(stats::quantile(y_raw, 0.75, na.rm = TRUE)),
    response_max = max(y_raw, na.rm = TRUE),
    n_above_or_equal_2 = sum(y_raw >= 2, na.rm = TRUE),
    n_below_2 = sum(y_raw < 2, na.rm = TRUE),
    best_csv = normalizePath(best_csv, mustWork = FALSE),
    mode_tables_dir = normalizePath(mode_tables_dir, mustWork = FALSE),
    output_dir = normalizePath(ref_dir, mustWork = FALSE),
    top_main_feature = if (nrow(main_tab)) main_tab$feature_name[[1L]] else NA_character_,
    top_interaction_feature = if (nrow(pair_tab)) pair_tab$feature_name[[1L]] else NA_character_,
    top_combined_feature = if (nrow(combined)) combined$feature_name[[1L]] else NA_character_,
    cumulative_r2_top_n = if (!is.null(cumulative_r2)) nrow(cumulative_r2) else NA_integer_,
    n_phase_plot_files = length(phase_paths),
    stringsAsFactors = FALSE
  )
  top3_joint_fields <- data.frame(
    top3_joint_fit_status = if (!is.null(top3_joint) && "fit_status" %in% names(top3_joint)) top3_joint$fit_status[[1L]] else NA_character_,
    top3_joint_apparent_r2 = if (!is.null(top3_joint) && "apparent_r2" %in% names(top3_joint)) top3_joint$apparent_r2[[1L]] else NA_real_,
    top3_joint_apparent_rmse = if (!is.null(top3_joint) && "apparent_rmse" %in% names(top3_joint)) top3_joint$apparent_rmse[[1L]] else NA_real_,
    top3_joint_cv_r2 = if (!is.null(top3_joint) && "cv_r2" %in% names(top3_joint)) top3_joint$cv_r2[[1L]] else NA_real_,
    top3_joint_cv_rmse = if (!is.null(top3_joint) && "cv_rmse" %in% names(top3_joint)) top3_joint$cv_rmse[[1L]] else NA_real_,
    top3_joint_cv_k = if (!is.null(top3_joint) && "cv_k" %in% names(top3_joint)) top3_joint$cv_k[[1L]] else NA_integer_,
    top3_joint_cv_success_folds = if (!is.null(top3_joint) && "cv_success_folds" %in% names(top3_joint)) top3_joint$cv_success_folds[[1L]] else NA_integer_,
    stringsAsFactors = FALSE
  )
  summary <- cbind(summary_base, top3_fields, top3_joint_fields)

  write_csv(seed_df, file.path(ref_dir, "dominant_ploidy_parameter_seed_values.csv"))
  write_csv(main_tab, file.path(ref_dir, "dominant_ploidy_parameter_main_effects.csv"))
  write_csv(pair_tab, file.path(ref_dir, "dominant_ploidy_parameter_pairwise_interactions.csv"))
  write_csv(glmnet_tab, file.path(ref_dir, "dominant_ploidy_parameter_stability_selection.csv"))
  write_tsv_plain(glmnet_status, file.path(ref_dir, "dominant_ploidy_parameter_elastic_net_status.tsv"))
  write_csv(combined, file.path(ref_dir, "dominant_ploidy_parameter_feature_importance.csv"))
  write_tsv_plain(summary, file.path(ref_dir, "dominant_ploidy_parameter_contribution_summary.tsv"))

  list(
    summary = summary,
    top_features = head(combined, 30L),
    output_dir = ref_dir
  )
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE)), forced_target = NULL) {
  root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
  mode_contribution_target <- normalize_mode_contribution_target(
    forced_target %||% argv$mode_contribution_target %||% argv$contribution_target %||% argv$response_target %||% "mode"
  )
  best_csv <- normalizePath(path.expand(argv$best_csv %||% paper_best_params_csv("invivo", root_dir = root_dir)), mustWork = FALSE)
  mode_tables_dir <- normalizePath(path.expand(argv$mode_tables_dir %||% paper_fixo2_mode_tables_dir(root_dir = root_dir)), mustWork = FALSE)
  output_dir <- normalizePath(path.expand(argv$output_dir %||% mode_contribution_default_output_dir(root_dir, mode_contribution_target)), mustWork = FALSE)
  mode_reference_o2 <- as_num(argv$mode_reference_o2, 2)
  reference_o2 <- parse_mode_reference_o2_values(mode_reference_o2, argv$mode_reference_o2_values)
  overwrite <- as_bool(argv$overwrite, FALSE)
  write_global_summary <- as_bool(argv$write_global_summary, TRUE)
  merge_only <- as_bool(argv$merge_only, FALSE)
  n_bootstrap <- as_int(argv$n_bootstrap, 100L)
  sample_fraction <- as_num(argv$sample_fraction, 0.75)
  alpha <- as_num(argv$elastic_net_alpha, 0.5)
  random_seed <- as_int(argv$random_seed, 123L)
  top_n_interactions <- as_int(argv$top_n_interactions, 6L)
  tree_depth <- as_int(argv$tree_depth, 3L)
  tree_min_node <- as_int(argv$tree_min_node, 20L)
  cumulative_auc_top_n <- as_int(argv$cumulative_auc_top_n, 30L)

  index_path <- if (identical(mode_contribution_target, "dominant_ploidy")) {
    file.path(output_dir, "dominant_ploidy_parameter_contribution_index.csv")
  } else {
    file.path(output_dir, "mode_parameter_contribution_index.csv")
  }
  reference_outputs_complete <- all(vapply(reference_o2, function(o2) {
    if (identical(mode_contribution_target, "dominant_ploidy")) {
      dominant_ploidy_contribution_complete(output_dir, o2)
    } else {
      reference_contribution_complete(output_dir, o2)
    }
  }, logical(1)))
  if (!overwrite && !merge_only && isTRUE(write_global_summary) && file.exists(index_path) && reference_outputs_complete) {
    write_fixed_o2_attractor_distribution_plots(mode_tables_dir, overwrite = FALSE)
    message("Skipping ", mode_contribution_target_label(mode_contribution_target), " parameter contribution; existing index and per-reference outputs found: ", index_path)
    return(invisible(index_path))
  }
  if (isTRUE(merge_only)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    out <- if (identical(mode_contribution_target, "dominant_ploidy")) {
      write_dominant_ploidy_contribution_global_summary(output_dir, reference_o2)
    } else {
      write_mode_contribution_global_summary(output_dir, reference_o2)
    }
    write_fixed_o2_attractor_distribution_plots(mode_tables_dir, overwrite = overwrite)
    message("Merged ", mode_contribution_target_label(mode_contribution_target), " parameter contribution outputs: ", out)
    return(invisible(out))
  }
  if (!file.exists(best_csv)) stop("Missing best-parameter CSV: ", best_csv)
  if (!dir.exists(mode_tables_dir)) stop("Missing fixed-O2 mode table directory: ", mode_tables_dir)
  if (!is.finite(sample_fraction) || sample_fraction <= 0 || sample_fraction > 1) {
    stop("sample_fraction must be in (0, 1].")
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) stop("elastic_net_alpha must be in [0, 1].")
  if (!is.finite(cumulative_auc_top_n) || is.na(cumulative_auc_top_n) || cumulative_auc_top_n < 1L) {
    stop("cumulative_auc_top_n must be a positive integer.")
  }

  best_df <- read_csv_plain(best_csv)
  max_seeds <- as_int(argv$max_seeds, NA_integer_)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L && "seed" %in% names(best_df)) {
    best_df <- best_df[order(suppressWarnings(as.integer(best_df$seed))), , drop = FALSE]
    best_df <- best_df[seq_len(min(nrow(best_df), max_seeds)), , drop = FALSE]
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  args_df <- data.frame(
    argument = c(
      "result_root", "best_csv", "mode_tables_dir", "output_dir", "mode_contribution_target", "mode_reference_o2_values",
      "n_bootstrap", "sample_fraction", "elastic_net_alpha", "random_seed",
      "top_n_interactions", "tree_depth", "tree_min_node", "cumulative_auc_top_n", "max_seeds"
    ),
    value = c(
      root_dir, best_csv, mode_tables_dir, output_dir, mode_contribution_target,
      paste(format(reference_o2, scientific = FALSE, trim = TRUE), collapse = ","),
      as.character(n_bootstrap), as.character(sample_fraction), as.character(alpha), as.character(random_seed),
      as.character(top_n_interactions), as.character(tree_depth), as.character(tree_min_node),
      as.character(cumulative_auc_top_n), as.character(max_seeds)
    ),
    stringsAsFactors = FALSE
  )
  run_args_file <- if (identical(mode_contribution_target, "dominant_ploidy")) {
    "dominant_ploidy_parameter_contribution_run_arguments.tsv"
  } else {
    "mode_parameter_contribution_run_arguments.tsv"
  }
  write_tsv_plain(args_df, file.path(output_dir, run_args_file))

  runs <- lapply(reference_o2, function(o2) {
    if (identical(mode_contribution_target, "dominant_ploidy")) {
      if (!overwrite && dominant_ploidy_contribution_complete(output_dir, o2)) {
        message("Skipping reference O2=", format(o2, scientific = FALSE, trim = TRUE), "; existing dominant-ploidy contribution outputs are complete.")
        return(read_dominant_ploidy_contribution(output_dir, o2))
      }
      return(run_reference_dominant_ploidy_contribution(
        best_df = best_df,
        best_csv = best_csv,
        mode_tables_dir = mode_tables_dir,
        mode_reference_o2 = o2,
        output_dir = output_dir,
        n_bootstrap = n_bootstrap,
        sample_fraction = sample_fraction,
        alpha = alpha,
        random_seed = random_seed,
        top_n_interactions = top_n_interactions,
        cumulative_auc_top_n = cumulative_auc_top_n
      ))
    }
    if (!overwrite && reference_contribution_complete(output_dir, o2)) {
      message("Skipping reference O2=", format(o2, scientific = FALSE, trim = TRUE), "; existing mode contribution outputs are complete.")
      return(read_reference_contribution(output_dir, o2))
    }
    run_reference_contribution(
      best_df = best_df,
      best_csv = best_csv,
      mode_tables_dir = mode_tables_dir,
      mode_reference_o2 = o2,
      output_dir = output_dir,
      n_bootstrap = n_bootstrap,
      sample_fraction = sample_fraction,
      alpha = alpha,
      random_seed = random_seed,
      top_n_interactions = top_n_interactions,
      tree_depth = tree_depth,
      tree_min_node = tree_min_node,
      cumulative_auc_top_n = cumulative_auc_top_n
    )
  })

  if (isTRUE(write_global_summary)) {
    if (identical(mode_contribution_target, "dominant_ploidy")) {
      write_dominant_ploidy_contribution_global_summary(output_dir, reference_o2)
    } else {
      write_mode_contribution_global_summary(output_dir, reference_o2)
    }
  }
  write_fixed_o2_attractor_distribution_plots(mode_tables_dir, overwrite = overwrite)
  message(mode_contribution_target_label(mode_contribution_target), " parameter contribution analysis complete: ", output_dir)
}

if (sys.nframe() == 0) {
  main()
}
