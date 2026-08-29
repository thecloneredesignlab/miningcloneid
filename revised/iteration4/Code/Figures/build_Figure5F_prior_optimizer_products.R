#!/usr/bin/env Rscript

# Build the active Figure 5F products from the reconstructed differential-
# evolution (DE) initial populations and the optimizer endpoints of one
# selected joint pair per declared primary family. Neither distribution is a
# Bayesian posterior, a confidence distribution, or a biological replicate
# distribution.

options(stringsAsFactors = FALSE, warn = 1)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE))
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
})
source(file.path(script_dir, "util", "runtime", "workspace_paths.R"))
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")

paths <- list(
  selection = file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv"),
  initial_population = file.path(
    figure5_dir, "figure5f_de_initial_population_log2_ratios.rds"
  ),
  initial_context = file.path(
    figure5_dir, "figure5f_de_initial_population_context_values.rds"
  ),
  initial_config = file.path(
    figure5_dir, "figure5f_de_initial_population_config.tsv"
  ),
  initial_readiness = file.path(
    figure5_dir, "figure5f_de_initial_population_readiness.tsv"
  ),
  optimizer = file.path(figure5_dir, "figure5f_optimizer_solutions.tsv"),
  parameter_meta = file.path(figure5_dir, "parameter_function_groups.tsv")
)
missing <- unlist(paths, use.names = FALSE)
missing <- missing[!file.exists(missing)]
if (length(missing)) {
  stop(
    "Missing Figure 5F DE-initial/optimizer input(s):\n",
    paste(missing, collapse = "\n")
  )
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
write_tsv <- function(x, path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  allowed <- normalizePath(figure5_dir, mustWork = TRUE)
  if (!startsWith(normalized, paste0(allowed, .Platform$file.sep))) {
    stop("Refusing to write outside iteration4 Figure5 data: ", normalized)
  }
  utils::write.table(
    x, normalized, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
  invisible(normalized)
}

parameters <- c(
  "lam_max", "p_mis_base", "p_wgd",
  "p_misseg", "k_o_mis",
  "O2_crit", "n_O",
  "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
families <- JOINT_FAMILY_LEVELS

selection <- read_tsv(paths$selection)
initial_population <- readRDS(paths$initial_population)
initial_context <- readRDS(paths$initial_context)
initial_config <- read_tsv(paths$initial_config)
initial_readiness <- read_tsv(paths$initial_readiness)
optimizer_all <- read_tsv(paths$optimizer)
parameter_meta <- read_tsv(paths$parameter_meta)

required_selection <- c(
  "warmup_label", "invivo_seed", "invitro_seed", "selected_seed",
  "objective", "objective_invivo", "objective_invitro", "family",
  "subcluster", "separate_invivo_objective", "selected_for_figure5f",
  "selection_rule"
)
required_initial <- c(
  "family", "warmup_label", "joint_seed", "initial_member",
  "exact_warm_start", parameters
)
required_initial_context <- c(
  "family", "warmup_label", "joint_seed", "initial_member",
  "exact_warm_start",
  paste0("vivo__", parameters), paste0("vitro__", parameters)
)
required_initial_config <- c(
  "family", "warmup_label", "parameter", "transformation",
  "center_name", "delta_name", "center_init", "delta_init",
  "center_sampling_sd", "delta_sampling_sd", "joint_warmup_sigmaN",
  "n_joint_seeds", "NP_used", "n_initial_population_values",
  "sampling_definition", "displayed_quantity", "ratio_formula"
)
required_optimizer <- c(
  "family", "subcluster", "pair_id", "warmup_label", "seed_number",
  "objective", "parameter", "parameter_group", "parameter_order",
  "ratio_vivo_to_vitro", "log2_ratio", "feasible_at_solution",
  "projection_applied", "either_context_at_active_bound"
)
if (!all(required_selection %in% names(selection)) ||
    !all(required_initial %in% names(initial_population)) ||
    !all(required_initial_context %in% names(initial_context)) ||
    !all(required_initial_config %in% names(initial_config)) ||
    !all(required_optimizer %in% names(optimizer_all)) ||
    !all(c("check", "passed") %in% names(initial_readiness))) {
  stop("Figure 5F DE-initial/optimizer input schema is incomplete.")
}
if (!all(as.logical(initial_readiness$passed))) {
  stop(
    "Figure 5F DE initial-population release checks failed upstream: ",
    paste(initial_readiness$check[!as.logical(initial_readiness$passed)],
          collapse = ", ")
  )
}

selected <- selection[as.logical(selection$selected_for_figure5f), , drop = FALSE]
selected <- selected[order(match(selected$family, families)), , drop = FALSE]
expected_selection <- stats::setNames(
  selection$warmup_label[match(families, selection$family)], families
)
observed_selection <- setNames(selected$warmup_label, selected$family)
if (nrow(selected) != length(families) ||
    !identical(selected$family, families) ||
    !identical(observed_selection[families], expected_selection) ||
    any(selected$invitro_seed != INVITRO_VISUALIZATION_SEED)) {
  stop("Figure 5F selected-pair contract does not match the approved inputs.")
}

optimizer <- optimizer_all[
  optimizer_all$warmup_label %in% selected$warmup_label,
  , drop = FALSE
]
initial_counts <- table(initial_population$family)
initial_seed_counts <- table(
  initial_population$family, initial_population$joint_seed
)
optimizer_counts <- table(optimizer$family, optimizer$parameter)
if (!identical(sort(unique(initial_population$family)), sort(families)) ||
    any(initial_counts[families] != 200000L) ||
    any(initial_seed_counts[families, , drop = FALSE] != 400L) ||
    any(!vapply(
      initial_population[, parameters, drop = FALSE],
      function(x) all(is.finite(x)), logical(1)
    ))) {
  stop("Figure 5F reconstructed DE initial-population coverage is invalid.")
}
if (nrow(initial_context) != length(families) * 500L * 400L ||
    !identical(initial_context[, c(
      "family", "warmup_label", "joint_seed", "initial_member",
      "exact_warm_start"
    )], initial_population[, c(
      "family", "warmup_label", "joint_seed", "initial_member",
      "exact_warm_start"
    )]) ||
    any(!vapply(
      initial_context[, c(
        paste0("vivo__", parameters), paste0("vitro__", parameters)
      ), drop = FALSE],
      function(x) all(is.finite(x) & x > 0), logical(1)
    ))) {
  stop("Figure 5F reconstructed DE context-value coverage is invalid.")
}
if (!setequal(unique(optimizer$family), families) ||
    !setequal(unique(optimizer$parameter), parameters) ||
    length(optimizer_counts) != 14L * length(families) ||
    any(optimizer_counts != 500L) ||
    any(!is.finite(optimizer$log2_ratio)) ||
    any(!is.finite(optimizer$ratio_vivo_to_vitro)) ||
    any(optimizer$ratio_vivo_to_vitro <= 0) ||
    any(!optimizer$feasible_at_solution) ||
    any(optimizer$projection_applied)) {
  stop("Figure 5F selected optimizer-endpoint coverage is invalid.")
}

meta <- parameter_meta[
  match(parameters, parameter_meta$parameter),
  c("parameter", "parameter_group", "parameter_order", "transformation"),
  drop = FALSE
]
initial_config <- initial_config[order(
  match(initial_config$parameter, parameters),
  match(initial_config$family, families)
), , drop = FALSE]
if (anyNA(meta$parameter) || !identical(meta$parameter, parameters) ||
    nrow(initial_config) != 14L * length(families) ||
    any(table(initial_config$parameter, initial_config$family) != 1L) ||
    any(initial_config$NP_used != 400L) ||
    any(initial_config$n_joint_seeds != 500L) ||
    any(initial_config$n_initial_population_values != 200000L) ||
    any(abs(initial_config$joint_warmup_sigmaN - 0.1216) > 1e-12)) {
  stop("Figure 5F parameter/DE-initialization metadata is not aligned.")
}

quantile_row <- function(values, prefix) {
  q <- stats::quantile(
    values, c(0.001, 0.05, 0.25, 0.50, 0.75, 0.95, 0.999),
    names = FALSE, type = 8
  )
  out <- data.frame(
    q001 = q[[1L]], q05 = q[[2L]], q25 = q[[3L]], median = q[[4L]],
    q75 = q[[5L]], q95 = q[[6L]], q999 = q[[7L]],
    width90 = q[[6L]] - q[[2L]], stringsAsFactors = FALSE
  )
  names(out) <- paste0(prefix, names(out))
  out
}

initial_summary <- do.call(rbind, lapply(parameters, function(parameter) {
  do.call(rbind, lapply(families, function(family) {
    values <- initial_population[[parameter]][
      initial_population$family == family
    ]
    config_row <- initial_config[
      initial_config$parameter == parameter &
        initial_config$family == family,
      , drop = FALSE
    ]
    cbind(
      data.frame(
        parameter = parameter,
        family = family,
        n_de_initial_values = length(values),
        n_de_seeds = config_row$n_joint_seeds,
        de_population_per_seed = config_row$NP_used,
        joint_warmup_sigmaN = config_row$joint_warmup_sigmaN,
        de_center_init = config_row$center_init,
        de_delta_init = config_row$delta_init,
        de_center_sampling_sd = config_row$center_sampling_sd,
        de_delta_sampling_sd = config_row$delta_sampling_sd,
        de_sampling_definition = config_row$sampling_definition,
        de_ratio_formula = config_row$ratio_formula,
        stringsAsFactors = FALSE
      ),
      meta[
        meta$parameter == parameter,
        c("parameter_group", "parameter_order", "transformation"),
        drop = FALSE
      ],
      quantile_row(values, "de_initial_")
    )
  }))
}))
rownames(initial_summary) <- NULL

optimizer_groups <- split(
  optimizer,
  interaction(optimizer$parameter, optimizer$family, drop = TRUE)
)
optimizer_summary <- do.call(rbind, lapply(optimizer_groups, function(group) {
  parameter <- group$parameter[[1L]]
  family <- group$family[[1L]]
  selected_row <- selected[selected$family == family, , drop = FALSE]
  cbind(
    data.frame(
      parameter = parameter,
      family = family,
      warmup_label = selected_row$warmup_label,
      invivo_seed = selected_row$invivo_seed,
      invitro_seed = selected_row$invitro_seed,
      selected_joint_seed = selected_row$selected_seed,
      joint_objective = selected_row$objective,
      separate_invivo_objective = selected_row$separate_invivo_objective,
      n_optimizer_endpoints = nrow(group),
      optimizer_active_bound_fraction =
        mean(group$either_context_at_active_bound),
      optimizer_fraction_higher_invivo = mean(group$log2_ratio > 0),
      stringsAsFactors = FALSE
    ),
    quantile_row(group$log2_ratio, "optimizer_")
  )
}))
rownames(optimizer_summary) <- NULL
optimizer_summary <- merge(
  optimizer_summary,
  initial_summary,
  by = c("parameter", "family"), all.x = TRUE, sort = FALSE
)
optimizer_summary <- optimizer_summary[order(
  match(optimizer_summary$parameter, parameters),
  match(optimizer_summary$family, families)
), , drop = FALSE]
optimizer_summary$endpoint_to_initial_width90_ratio <-
  optimizer_summary$optimizer_width90 /
  optimizer_summary$de_initial_width90
optimizer_summary$optimizer_direction <- ifelse(
  optimizer_summary$optimizer_q05 > 0,
  "higher in vivo",
  ifelse(
    optimizer_summary$optimizer_q95 < 0,
    "lower in vivo",
    "overlaps equality"
  )
)
optimizer_summary$interval_semantics <- paste0(
  "descriptive 5th-95th percentile spans across DE initial-population ",
  "members and optimizer endpoints"
)

cross_family <- do.call(rbind, lapply(parameters, function(parameter) {
  rows <- optimizer_summary[
    optimizer_summary$parameter == parameter, , drop = FALSE
  ]
  directions <- rows$optimizer_direction
  all_same_direction <- length(unique(directions)) == 1L &&
    directions[[1L]] != "overlaps equality"
  width_ratios <- rows$endpoint_to_initial_width90_ratio
  data.frame(
    parameter = parameter,
    parameter_group = rows$parameter_group[[1L]],
    parameter_order = rows$parameter_order[[1L]],
    optimizer_direction_C01 = directions[rows$family == "C01"],
    optimizer_direction_C02 = directions[rows$family == "C02"],
    optimizer_direction_C03 = directions[rows$family == "C03"],
    optimizer_direction_C04 = directions[rows$family == "C04"],
    optimizer_direction_C05 = directions[rows$family == "C05"],
    optimizer_direction_C06 = directions[rows$family == "C06"],
    cross_family_direction_agreement = all_same_direction,
    cross_family_direction = if (all_same_direction) {
      directions[[1L]]
    } else {
      "family-dependent or overlaps equality"
    },
    optimizer_median_min = min(rows$optimizer_median),
    optimizer_median_max = max(rows$optimizer_median),
    optimizer_median_span = diff(range(rows$optimizer_median)),
    endpoint_to_initial_width90_ratio_min = min(width_ratios),
    endpoint_to_initial_width90_ratio_max = max(width_ratios),
    optimizer_active_bound_fraction_max =
      max(rows$optimizer_active_bound_fraction),
    all_families_narrower_than_initial = all(width_ratios < 1),
    all_families_narrower_than_half_initial = all(width_ratios <= 0.5),
    fitted_solution_concentration = if (all(width_ratios <= 0.5)) {
      "narrower than half the DE-initial width in all families"
    } else if (all(width_ratios < 1)) {
      "narrower than the DE-initial width in all families"
    } else {
      "not consistently narrower than the DE-initial width"
    },
    inference_limit = paste0(
      "descriptive optimizer-solution concentration relative to the DE ",
      "initial population; not posterior uncertainty or proof of ",
      "structural identifiability"
    ),
    stringsAsFactors = FALSE
  )
}))
rownames(cross_family) <- NULL

safe_density <- function(values, grid_from, grid_to) {
  span <- grid_to - grid_from
  raw_bw <- suppressWarnings(stats::bw.nrd0(values))
  bw <- max(raw_bw, span / 800, 1e-6)
  estimate <- stats::density(
    values, bw = bw, from = grid_from, to = grid_to,
    n = 401, cut = 0
  )
  data.frame(
    x = estimate$x,
    density = estimate$y,
    density_scaled_distribution = estimate$y / max(estimate$y),
    bandwidth = bw,
    stringsAsFactors = FALSE
  )
}

density_rows <- list()
counter <- 0L
for (parameter in parameters) {
  initial_values_all <- initial_population[[parameter]]
  endpoint_values_all <- optimizer$log2_ratio[
    optimizer$parameter == parameter
  ]
  display_quantiles <- stats::quantile(
    c(initial_values_all, endpoint_values_all), c(0.0005, 0.9995),
    names = FALSE, type = 8
  )
  display_span <- diff(display_quantiles)
  if (!is.finite(display_span) || display_span <= 0) display_span <- 1
  display_from <- min(display_quantiles[[1L]], 0) - 0.04 * display_span
  display_to <- max(display_quantiles[[2L]], 0) + 0.04 * display_span

  for (family in families) {
    initial_values <- initial_population[[parameter]][
      initial_population$family == family
    ]
    endpoint_values <- optimizer$log2_ratio[
      optimizer$parameter == parameter & optimizer$family == family
    ]
    initial_density <- safe_density(
      initial_values, display_from, display_to
    )
    endpoint_density <- safe_density(
      endpoint_values, display_from, display_to
    )
    for (role in c("de_initial_population", "optimizer_endpoints")) {
      values <- if (role == "de_initial_population") {
        initial_values
      } else {
        endpoint_values
      }
      curve <- if (role == "de_initial_population") {
        initial_density
      } else {
        endpoint_density
      }
      counter <- counter + 1L
      density_rows[[counter]] <- cbind(
        data.frame(
          parameter = parameter,
          family = family,
          distribution_role = role,
          display_from = display_from,
          display_to = display_to,
          displayed_tail_fraction = mean(
            values < display_from | values > display_to
          ),
          n_values = length(values),
          stringsAsFactors = FALSE
        ),
        curve
      )
    }
  }
}
density_table <- do.call(rbind, density_rows)
rownames(density_table) <- NULL
density_table <- merge(
  density_table, meta,
  by = "parameter", all.x = TRUE, sort = FALSE
)
density_table <- density_table[order(
  match(density_table$parameter, parameters),
  match(density_table$family, families),
  match(
    density_table$distribution_role,
    c("de_initial_population", "optimizer_endpoints")
  ),
  density_table$x
), , drop = FALSE]

# Direct context-specific distributions for the active Figure 5F display.
# Density estimation is performed on the parameter's configured optimizer
# scale (log10 or identity), while summaries remain on the natural scale.
contexts <- c("in vivo", "in vitro")
context_initial_column <- c(
  "in vivo" = "vivo__",
  "in vitro" = "vitro__"
)
context_endpoint_column <- c(
  "in vivo" = "vivo_natural",
  "in vitro" = "vitro_natural"
)
context_bound_column <- c(
  "in vivo" = "vivo_at_active_bound",
  "in vitro" = "vitro_at_active_bound"
)

to_display_scale <- function(values, transformation) {
  if (identical(transformation, "log10")) log10(values) else values
}

context_summary_rows <- list()
context_density_rows <- list()
context_summary_index <- 0L
context_density_index <- 0L
for (parameter in parameters) {
  meta_row <- meta[meta$parameter == parameter, , drop = FALSE]
  transformation <- meta_row$transformation[[1L]]
  all_natural_values <- c(
    initial_context[[paste0("vivo__", parameter)]],
    initial_context[[paste0("vitro__", parameter)]],
    optimizer$vivo_natural[optimizer$parameter == parameter],
    optimizer$vitro_natural[optimizer$parameter == parameter]
  )
  all_display_values <- to_display_scale(all_natural_values, transformation)
  display_quantiles <- stats::quantile(
    all_display_values, c(0.0005, 0.9995), names = FALSE, type = 8
  )
  display_span <- diff(display_quantiles)
  if (!is.finite(display_span) || display_span <= 0) display_span <- 1
  display_from <- display_quantiles[[1L]] - 0.04 * display_span
  display_to <- display_quantiles[[2L]] + 0.04 * display_span

  for (family in families) {
    selected_row <- selected[selected$family == family, , drop = FALSE]
    initial_family <- initial_context$family == family
    endpoint_family_parameter <-
      optimizer$family == family & optimizer$parameter == parameter
    endpoint_group <- optimizer[endpoint_family_parameter, , drop = FALSE]

    for (context_name in contexts) {
      initial_column <- paste0(
        context_initial_column[[context_name]], parameter
      )
      endpoint_column <- context_endpoint_column[[context_name]]
      bound_column <- context_bound_column[[context_name]]
      initial_values <- initial_context[[initial_column]][initial_family]
      endpoint_values <- endpoint_group[[endpoint_column]]
      if (length(initial_values) != 200000L ||
          length(endpoint_values) != 500L ||
          any(!is.finite(initial_values) | initial_values <= 0) ||
          any(!is.finite(endpoint_values) | endpoint_values <= 0)) {
        stop(
          "Invalid Figure 5F context distribution for ", parameter,
          ", ", family, ", ", context_name
        )
      }

      context_summary_index <- context_summary_index + 1L
      context_summary_rows[[context_summary_index]] <- cbind(
        data.frame(
          parameter = parameter,
          family = family,
          context = context_name,
          warmup_label = selected_row$warmup_label,
          invivo_seed = selected_row$invivo_seed,
          invitro_seed = selected_row$invitro_seed,
          selected_joint_seed = selected_row$selected_seed,
          joint_objective = selected_row$objective,
          n_de_initial_values = length(initial_values),
          n_optimizer_endpoints = length(endpoint_values),
          optimizer_active_bound_fraction = mean(endpoint_group[[bound_column]]),
          parameter_group = meta_row$parameter_group,
          parameter_order = meta_row$parameter_order,
          transformation = transformation,
          summary_scale = "natural parameter scale",
          stringsAsFactors = FALSE
        ),
        quantile_row(initial_values, "de_initial_"),
        quantile_row(endpoint_values, "optimizer_")
      )

      role_values <- list(
        de_initial_population = initial_values,
        optimizer_endpoints = endpoint_values
      )
      for (role in names(role_values)) {
        natural_values <- role_values[[role]]
        display_values <- to_display_scale(
          natural_values, transformation
        )
        curve <- safe_density(
          display_values, display_from, display_to
        )
        context_density_index <- context_density_index + 1L
        context_density_rows[[context_density_index]] <- cbind(
          data.frame(
            parameter = parameter,
            family = family,
            context = context_name,
            distribution_role = role,
            transformation = transformation,
            display_from = display_from,
            display_to = display_to,
            displayed_tail_fraction = mean(
              display_values < display_from | display_values > display_to
            ),
            n_values = length(display_values),
            stringsAsFactors = FALSE
          ),
          curve
        )
      }
    }
  }
}

context_summary <- do.call(rbind, context_summary_rows)
rownames(context_summary) <- NULL
context_summary <- context_summary[order(
  match(context_summary$parameter, parameters),
  match(context_summary$family, families),
  match(context_summary$context, contexts)
), , drop = FALSE]

context_density <- do.call(rbind, context_density_rows)
rownames(context_density) <- NULL
context_density$x_natural <- ifelse(
  context_density$transformation == "log10",
  10^context_density$x,
  context_density$x
)
context_density <- merge(
  context_density,
  meta[, c("parameter", "parameter_group", "parameter_order")],
  by = "parameter", all.x = TRUE, sort = FALSE
)
context_density <- context_density[order(
  match(context_density$parameter, parameters),
  match(context_density$family, families),
  match(context_density$context, contexts),
  match(
    context_density$distribution_role,
    c("de_initial_population", "optimizer_endpoints")
  ),
  context_density$x
), , drop = FALSE]

density_cell_counts <- table(
  density_table$parameter,
  density_table$family,
  density_table$distribution_role
)
density_peaks <- tapply(
  density_table$density_scaled_distribution,
  interaction(
    density_table$parameter,
    density_table$family,
    density_table$distribution_role,
    drop = TRUE
  ),
  max
)
context_density_cell_counts <- table(
  context_density$parameter,
  context_density$family,
  context_density$context,
  context_density$distribution_role
)
context_density_peaks <- tapply(
  context_density$density_scaled_distribution,
  interaction(
    context_density$parameter,
    context_density$family,
    context_density$context,
    context_density$distribution_role,
    drop = TRUE
  ),
  max
)

checks <- data.frame(
  check = c(
    "selected_one_pair_per_family",
    "selected_pairs_match_approved_contract",
    "de_initial_values_per_family_parameter",
    "de_initial_seeds_per_family",
    "de_initial_population_members_per_seed",
    "optimizer_endpoints_per_family_parameter",
    "optimizer_endpoints_feasible",
    "optimizer_endpoints_unprojected",
    "optimizer_ratios_positive_and_finite",
    "density_rows",
    "density_grid_points_per_cell",
    "density_peak_normalized_separately",
    "displayed_tail_fraction_at_most_0p002",
    "summary_rows",
    "cross_family_rows",
    "de_initial_context_values_per_family_parameter_context",
    "optimizer_endpoint_context_values_per_family_parameter_context",
    "context_density_rows",
    "context_density_grid_points_per_cell",
    "context_density_peak_normalized_separately",
    "context_displayed_tail_fraction_at_most_0p002",
    "context_summary_rows",
    "active_panel_uses_direct_context_values",
    "posterior_claim_absent"
  ),
  observed = c(
    nrow(selected),
    paste(observed_selection[families], collapse = ";"),
    paste(sort(unique(initial_summary$n_de_initial_values)), collapse = ","),
    paste(sort(unique(initial_summary$n_de_seeds)), collapse = ","),
    paste(sort(unique(initial_summary$de_population_per_seed)), collapse = ","),
    paste(sort(unique(as.vector(optimizer_counts))), collapse = ","),
    all(optimizer$feasible_at_solution),
    !any(optimizer$projection_applied),
    all(is.finite(optimizer$ratio_vivo_to_vitro) &
          optimizer$ratio_vivo_to_vitro > 0),
    nrow(density_table),
    paste(sort(unique(as.vector(density_cell_counts))), collapse = ","),
    all(abs(density_peaks - 1) <= 1e-12),
    max(density_table$displayed_tail_fraction),
    nrow(optimizer_summary),
    nrow(cross_family),
    paste(sort(unique(context_summary$n_de_initial_values)), collapse = ","),
    paste(sort(unique(context_summary$n_optimizer_endpoints)), collapse = ","),
    nrow(context_density),
    paste(sort(unique(as.vector(context_density_cell_counts))), collapse = ","),
    all(abs(context_density_peaks - 1) <= 1e-12),
    max(context_density$displayed_tail_fraction),
    nrow(context_summary),
    TRUE,
    TRUE
  ),
  expected = c(
    "6",
    paste(expected_selection[families], collapse = ";"),
    "200000",
    "500",
    "400",
    "500",
    "TRUE",
    "TRUE",
    "TRUE",
    as.character(14L * length(families) * 2L * 401L),
    "401",
    "TRUE",
    "0.002",
    "84",
    "14",
    "200000",
    "500",
    as.character(14L * length(families) * 2L * 2L * 401L),
    "401",
    "TRUE",
    "0.002",
    "168",
    "TRUE",
    "TRUE"
  ),
  stringsAsFactors = FALSE
)
checks$passed <- mapply(function(observed, expected, check) {
  if (grepl("displayed_tail_fraction_at_most_0p002$", check)) {
    return(as.numeric(observed) <= as.numeric(expected))
  }
  identical(as.character(observed), as.character(expected))
}, checks$observed, checks$expected, checks$check)

provenance <- data.frame(
  role = c(
    "selected pair audit",
    "reconstructed DE initial-population ratios",
    "reconstructed DE initial-population natural context values",
    "DE initial-population configuration",
    "DE initial-population upstream readiness",
    "optimizer endpoints for all declared pairs",
    "parameter taxonomy"
  ),
  path = normalizePath(unlist(paths, use.names = FALSE), mustWork = TRUE),
  md5 = unname(tools::md5sum(unlist(paths, use.names = FALSE))),
  stringsAsFactors = FALSE
)

outputs <- list(
  density = file.path(figure5_dir, "figure5f_prior_optimizer_density.tsv"),
  summary = file.path(figure5_dir, "figure5f_prior_optimizer_summary.tsv"),
  cross = file.path(figure5_dir, "figure5f_prior_optimizer_cross_family.tsv"),
  readiness = file.path(figure5_dir, "figure5f_prior_optimizer_readiness.tsv"),
  provenance = file.path(figure5_dir, "figure5f_prior_optimizer_provenance.tsv")
)
outputs$context_density <- file.path(
  figure5_dir, "figure5f_context_initial_optimizer_density.tsv"
)
outputs$context_summary <- file.path(
  figure5_dir, "figure5f_context_initial_optimizer_summary.tsv"
)
write_tsv(density_table, outputs$density)
write_tsv(optimizer_summary, outputs$summary)
write_tsv(cross_family, outputs$cross)
write_tsv(checks, outputs$readiness)
write_tsv(provenance, outputs$provenance)
write_tsv(context_density, outputs$context_density)
write_tsv(context_summary, outputs$context_summary)

chart_contract <- c(
  "# Figure 5D chart contract",
  "",
  paste0(
    "The active panel compares the exact differential-evolution initial populations with optimizer endpoints for the 14 soft-coupled parameters in ",
    paste(families, collapse = ", "), "."
  ),
  "",
  "- One row is shown per parameter and one same-size distribution column per primary joint-fit family.",
  "- The upper half is in vitro and the lower half is in vivo.",
  "- Gray dashed densities are DE initial populations; solid family-colored outlines are optimizer endpoints.",
  "- Endpoint fill and median markers encode context; red edge bands mark the outer 5% of the complete joint-union optimizer interval.",
  "- KDE is evaluated on the configured optimizer scale and displayed with natural-scale tick labels.",
  "- Each family contains 500 joint seeds and each reconstructed DE population contains 400 members per seed.",
  sprintf("- Direct-context density table: %s rows.", nrow(context_density)),
  sprintf("- Cross-family endpoint direction agrees for %s of %s parameters.", sum(cross_family$cross_family_direction_agreement), nrow(cross_family)),
  "",
  "These distributions describe numerical initialization and optimization endpoints; they are not posterior distributions, confidence distributions, or biological replicates."
)
writeLines(
  chart_contract,
  file.path(figure5_dir, "figure5f_chart_contract.md"),
  useBytes = TRUE
)

if (any(!checks$passed)) {
  stop(
    "Figure 5F DE-initial/optimizer release checks failed: ",
    paste(checks$check[!checks$passed], collapse = ", ")
  )
}

cat("Built Figure 5F DE-initial/optimizer products.\n")
cat("Selected families:", paste(families, collapse = ", "), "\n")
cat("Density rows:", nrow(density_table), "\n")
cat("Direct-context density rows:", nrow(context_density), "\n")
cat(
  "Cross-family endpoint directional agreement:",
  sum(cross_family$cross_family_direction_agreement), "of", nrow(cross_family),
  "parameters\n"
)
