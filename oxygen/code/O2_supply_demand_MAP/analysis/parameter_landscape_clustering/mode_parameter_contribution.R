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
    summary = file.path(ref_dir, "mode_parameter_contribution_summary.tsv")
  )
}

reference_contribution_complete <- function(output_dir, mode_reference_o2) {
  paths <- reference_contribution_paths(output_dir, mode_reference_o2)
  all(file.exists(c(paths$feature_importance, paths$summary)))
}

read_reference_contribution <- function(output_dir, mode_reference_o2) {
  paths <- reference_contribution_paths(output_dir, mode_reference_o2)
  if (!file.exists(paths$summary)) stop("Missing reference summary: ", paths$summary)
  if (!file.exists(paths$feature_importance)) stop("Missing reference feature importance: ", paths$feature_importance)
  summary <- read_tsv(paths$summary)
  top <- head(read_csv_plain(paths$feature_importance), 30L)
  list(summary = summary, top_features = top, output_dir = paths$reference_dir)
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

  summary <- data.frame(
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
