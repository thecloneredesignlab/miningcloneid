#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
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
                                 width = 7.8) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(character()))
  df <- auc_df[is.finite(auc_df$auc_discriminative), , drop = FALSE]
  df <- df[order(-df$auc_discriminative, df$feature_name), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n) && top_n > 0L) {
    df <- head(df, as.integer(top_n))
  }
  if (!nrow(df)) return(invisible(character()))
  df$feature_index <- seq_len(nrow(df))
  height <- max(4.2, min(30, 2.2 + 0.16 * nrow(df)))
  p <- ggplot2::ggplot(df, ggplot2::aes(fill = direction)) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = feature_index - 0.36,
        xmax = feature_index + 0.36,
        ymin = 0.5,
        ymax = auc_discriminative
      )
    ) +
    ggplot2::coord_flip() +
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
      x = NULL,
      y = "AUC",
      fill = "Direction"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
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
      p <- ggplot2::ggplot(roc, ggplot2::aes(x = fpr, y = tpr, color = curve)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
        ggplot2::geom_path(linewidth = 0.9) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        ggplot2::scale_color_manual(values = c(apparent = "#1f78b4", cross_validated = "#e31a1c")) +
        ggplot2::labs(
          title = paste0("Top3 joint ROC, reference O2=", format(mode_reference_o2, scientific = FALSE, trim = TRUE)),
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
      paste0("Main feature AUC, reference O2=", format(mode_reference_o2, scientific = FALSE, trim = TRUE))
    )
  }
  if (nrow(auc_tab)) {
    plot_feature_auc_bar(
      auc_tab,
      file.path(auc_dir, "all_feature_auc_bar"),
      paste0("All feature AUC, reference O2=", format(mode_reference_o2, scientific = FALSE, trim = TRUE))
    )
    plot_feature_auc_bar(
      auc_tab,
      file.path(auc_dir, "combined_feature_auc_bar_top30"),
      paste0("Top30 feature AUC, reference O2=", format(mode_reference_o2, scientific = FALSE, trim = TRUE)),
      top_n = 30L
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
        subtitle = "Best in vivo parameters colored by fixed-O2 mode",
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
    summary = file.path(ref_dir, "mode_parameter_contribution_summary.tsv")
  )
}

reference_contribution_complete <- function(output_dir, mode_reference_o2) {
  paths <- reference_contribution_paths(output_dir, mode_reference_o2)
  all(file.exists(c(paths$feature_importance, paths$summary, paths$auc_by_feature, paths$top3_joint_auc)))
}

read_reference_contribution <- function(output_dir, mode_reference_o2) {
  paths <- reference_contribution_paths(output_dir, mode_reference_o2)
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
  paste0("O2=", format(suppressWarnings(as.numeric(o2)), scientific = FALSE, trim = TRUE))
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
        color = metric_label,
        shape = metric_label,
        group = metric_label
      )
    ) +
      ggplot2::geom_line(linewidth = 0.75, alpha = 0.9) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_grid(. ~ o2_group, scales = "free_x", space = "free_x") +
      ggplot2::scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
      ggplot2::scale_color_manual(values = c(
        "Best feature AUC" = "#1f78b4",
        "Top3 mean feature AUC" = "#33a02c",
        "Top3 joint AUC" = "#e31a1c"
      )) +
      ggplot2::labs(
        title = "Mode contribution AUC across reference O2 groups",
        x = "Reference O2 (%)",
        y = "AUC",
        color = NULL,
        shape = NULL
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
        legend.position = "bottom"
      )
    save_plot_pair(p, file.path(plots_dir, "o2_group_auc_summary"), width = 8.4, height = 4.8)

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
          color = feature_rank_label,
          shape = feature_rank_label,
          group = feature_rank_label
        )
      ) +
        ggplot2::geom_line(linewidth = 0.75, alpha = 0.9) +
        ggplot2::geom_point(size = 2.4) +
        ggplot2::facet_grid(. ~ o2_group, scales = "free_x", space = "free_x") +
        ggplot2::scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
        ggplot2::scale_color_manual(values = c(Top1 = "#1f78b4", Top2 = "#33a02c", Top3 = "#ff7f00")) +
        ggplot2::labs(
          title = "Top3 retained feature AUC across reference O2 groups",
          x = "Reference O2 (%)",
          y = "AUC",
          color = NULL,
          shape = NULL
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
          legend.position = "bottom"
        )
      save_plot_pair(p_top3, file.path(plots_dir, "o2_group_top3_feature_auc"), width = 8.4, height = 4.8)
    }
  }
  invisible(file.path(tables_dir, "o2_group_auc_summary.csv"))
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
                                       tree_min_node) {
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

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
  tables_dir <- argv$tables_dir %||% paper_tables_dir("invivo", root_dir = root_dir)
  best_csv <- normalizePath(path.expand(argv$best_csv %||% file.path(tables_dir, "invivo_best_params_by_seed.csv")), mustWork = FALSE)
  mode_tables_dir <- normalizePath(path.expand(argv$mode_tables_dir %||% paper_fixo2_mode_tables_dir(root_dir = root_dir)), mustWork = FALSE)
  output_dir <- normalizePath(path.expand(argv$output_dir %||% file.path(root_dir, "mode_parameter_contribution")), mustWork = FALSE)
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

  index_path <- file.path(output_dir, "mode_parameter_contribution_index.csv")
  if (!overwrite && !merge_only && isTRUE(write_global_summary) && file.exists(index_path)) {
    message("Skipping mode parameter contribution; existing index found: ", index_path)
    return(invisible(index_path))
  }
  if (isTRUE(merge_only)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    out <- write_mode_contribution_global_summary(output_dir, reference_o2)
    message("Merged mode parameter contribution outputs: ", out)
    return(invisible(out))
  }
  if (!file.exists(best_csv)) stop("Missing best-parameter CSV: ", best_csv)
  if (!dir.exists(mode_tables_dir)) stop("Missing fixed-O2 mode table directory: ", mode_tables_dir)
  if (!is.finite(sample_fraction) || sample_fraction <= 0 || sample_fraction > 1) {
    stop("sample_fraction must be in (0, 1].")
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) stop("elastic_net_alpha must be in [0, 1].")

  best_df <- read_csv_plain(best_csv)
  max_seeds <- as_int(argv$max_seeds, NA_integer_)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L && "seed" %in% names(best_df)) {
    best_df <- best_df[order(suppressWarnings(as.integer(best_df$seed))), , drop = FALSE]
    best_df <- best_df[seq_len(min(nrow(best_df), max_seeds)), , drop = FALSE]
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  args_df <- data.frame(
    argument = c(
      "result_root", "best_csv", "mode_tables_dir", "output_dir", "mode_reference_o2_values",
      "n_bootstrap", "sample_fraction", "elastic_net_alpha", "random_seed",
      "top_n_interactions", "tree_depth", "tree_min_node", "max_seeds"
    ),
    value = c(
      root_dir, best_csv, mode_tables_dir, output_dir,
      paste(format(reference_o2, scientific = FALSE, trim = TRUE), collapse = ","),
      as.character(n_bootstrap), as.character(sample_fraction), as.character(alpha), as.character(random_seed),
      as.character(top_n_interactions), as.character(tree_depth), as.character(tree_min_node), as.character(max_seeds)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_plain(args_df, file.path(output_dir, "mode_parameter_contribution_run_arguments.tsv"))

  runs <- lapply(reference_o2, function(o2) {
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
      tree_min_node = tree_min_node
    )
  })

  if (isTRUE(write_global_summary)) {
    write_mode_contribution_global_summary(output_dir, reference_o2)
  }
  message("Mode parameter contribution analysis complete: ", output_dir)
}

main()
