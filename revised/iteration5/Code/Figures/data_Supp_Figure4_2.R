#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))
source(resolve_utility_file("analysis/optimizer_diagnostics.R"))

data_Supp_Figure4_2 <- function() {
  output_root <- file.path(DATA_ROOT, "Supp_Figure4_2")
  metrics <- optimizer_collect_seed_metrics(INVIVO_RESULT_ROOT, "objective")
  selected <- optimizer_write_bundle(
    metrics = metrics,
    output_dir = output_root,
    fit_label = "in vivo",
    selected_seed = "seed25",
    run_dir = INVIVO_RESULT_ROOT,
    include_all_boundaries = TRUE
  )

  best_objective <- min(metrics$objective, na.rm = TRUE)
  relative_excess <- metrics$objective / best_objective - 1
  near_optimal <- data.frame(
    threshold_type = c("rank", "relative", "relative", "relative"),
    threshold = c(1, 0.01, 0.05, 0.10),
    n_seeds = c(
      1L,
      sum(relative_excess <= 0.01),
      sum(relative_excess <= 0.05),
      sum(relative_excess <= 0.10)
    ),
    total_seeds = nrow(metrics),
    selected_seed = "seed25",
    selected_objective = best_objective,
    stringsAsFactors = FALSE
  )
  write_intermediate_tsv(
    near_optimal,
    file.path(output_root, "objective_near_optimal_summary.tsv")
  )

  stop_reason <- ifelse(
    is.na(metrics$deoptim_stop_reason) |
      !nzchar(metrics$deoptim_stop_reason) |
      metrics$deoptim_stop_reason == "NA",
    "not recorded",
    metrics$deoptim_stop_reason
  )
  termination <- rbind(
    data.frame(
      diagnostic = "DEoptim stop reason",
      category = names(table(stop_reason)),
      n_seeds = as.integer(table(stop_reason)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      diagnostic = "Local refinement accepted",
      category = ifelse(c(FALSE, TRUE), "accepted", "not accepted"),
      n_seeds = c(
        sum(!metrics$optimizer_local_accepted, na.rm = TRUE),
        sum(metrics$optimizer_local_accepted, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      diagnostic = "Local convergence code",
      category = as.character(sort(unique(
        metrics$optimizer_local_convergence[
          is.finite(metrics$optimizer_local_convergence)
        ]
      ))),
      n_seeds = as.integer(table(factor(
        metrics$optimizer_local_convergence[
          is.finite(metrics$optimizer_local_convergence)
        ],
        levels = sort(unique(metrics$optimizer_local_convergence[
          is.finite(metrics$optimizer_local_convergence)
        ]))
      ))),
      stringsAsFactors = FALSE
    )
  )
  write_intermediate_tsv(
    termination,
    file.path(output_root, "termination_summary.tsv")
  )

  sources <- attr(selected, "source_files")
  contract <- data.frame(
    role = "in-vivo optimizer diagnostic input",
    source = sources,
    local_file = NA_character_,
    source_md5 = unname(tools::md5sum(sources)),
    local_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  write_data_contract("Supp_Figure4_2", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) data_Supp_Figure4_2()
