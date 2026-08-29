#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

run_all_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(
    n_core = 8L,
    recompute_fixed_o2 = FALSE,
    recompute_invivo_tsne = FALSE,
    rebuild_figure6_grid = FALSE,
    model_dependent_only = FALSE,
    first_main_figure = 1L,
    resume_after_figure5f_de = FALSE
  )
  for (arg in args) {
    if (grepl("^--n-core=", arg)) {
      out$n_core <- as.integer(sub("^--n-core=", "", arg))
    } else if (grepl("^--recompute-fixed-o2=", arg)) {
      out$recompute_fixed_o2 <- as_boolean(
        sub("^--recompute-fixed-o2=", "", arg)
      )
    } else if (grepl("^--recompute-invivo-tsne=", arg)) {
      out$recompute_invivo_tsne <- as_boolean(
        sub("^--recompute-invivo-tsne=", "", arg)
      )
    } else if (grepl("^--rebuild-figure6-grid=", arg)) {
      out$rebuild_figure6_grid <- as_boolean(
        sub("^--rebuild-figure6-grid=", "", arg)
      )
    } else if (grepl("^--model-dependent-only=", arg)) {
      out$model_dependent_only <- as_boolean(
        sub("^--model-dependent-only=", "", arg)
      )
    } else if (grepl("^--first-main-figure=", arg)) {
      out$first_main_figure <- as.integer(
        sub("^--first-main-figure=", "", arg)
      )
    } else if (grepl("^--resume-after-figure5f-de=", arg)) {
      out$resume_after_figure5f_de <- as_boolean(
        sub("^--resume-after-figure5f-de=", "", arg)
      )
    } else if (is_runtime_path_argument(arg)) {
      next
    } else {
      stop("Unknown run_all_figures.R option: ", arg)
    }
  }
  if (!is.finite(out$n_core) || out$n_core < 1L) {
    stop("--n-core must be a positive integer.")
  }
  if (!out$first_main_figure %in% c(1L, 3L, 4L)) {
    stop("--first-main-figure must be one of 1, 3, or 4.")
  }
  out
}

run_figure_entry <- function(script, args = character()) {
  path <- file.path(CODE_ROOT, script)
  require_files(path, "figure entry")
  message("\n==> ", script)
  run_process(
    "Rscript",
    args = c(path, runtime_path_arguments(), args)
  )
}

run_all_figures <- function(
    n_core = 8L,
    recompute_fixed_o2 = FALSE,
    recompute_invivo_tsne = FALSE,
    rebuild_figure6_grid = FALSE,
    model_dependent_only = FALSE,
    first_main_figure = 1L,
    resume_after_figure5f_de = FALSE
) {
  ensure_workspace_directories()

  if (!isTRUE(resume_after_figure5f_de)) {
    if (!isTRUE(model_dependent_only) && first_main_figure == 1L) {
      run_figure_entry("data_Figure1.R")
      run_figure_entry("draw_Figure1.R")
      run_figure_entry("data_Figure2.R")
      run_figure_entry("draw_Figure2.R")
    }
    if ((!isTRUE(model_dependent_only) && first_main_figure <= 3L) ||
        first_main_figure == 3L) {
      run_figure_entry("data_Figure3.R")
      run_figure_entry("draw_Figure3.R")
    } else {
      message("\nPreserving existing earlier figures; rebuilding Figure 4-6.")
    }
    run_figure_entry(
      "data_Figure4.R",
      c(
        paste0("--n-core=", n_core),
        paste0(
          "--recompute-fixed-o2=",
          if (isTRUE(recompute_fixed_o2)) "TRUE" else "FALSE"
        ),
        paste0(
          "--recompute-tsne=",
          if (isTRUE(recompute_invivo_tsne)) "TRUE" else "FALSE"
        )
      )
    )
    run_figure_entry("draw_Figure4.R")
    run_figure_entry("data_Supp_Figure4_1.R")
    run_figure_entry("draw_Supp_Figure4_1.R")
    run_figure_entry("data_Supp_Figure4_2.R")
    run_figure_entry("draw_Supp_Figure4_2.R")
    run_figure_entry("data_Figure5.R")
    run_figure_entry("data_Supp_Figure5_1.R")
    run_figure_entry("data_Supp_Figure5_2.R")
    run_figure_entry("finalize_Figure5_optimizer_ensemble.R")
    run_figure_entry("prepare_Figure5F_de_initial_population.R")
  } else {
    message(
      "\nResuming after validated Figure 5F DE initial-population products."
    )
  }
  run_figure_entry("audit_Figure5F_prior_optimizer_inputs.R")
  run_figure_entry("build_Figure5F_prior_optimizer_products.R")
  run_figure_entry("build_Figure5F_supplementary_table.R")
  run_figure_entry("draw_Figure5.R")
  run_figure_entry("draw_Supp_Figure5_1.R")
  run_figure_entry("draw_Supp_Figure5_2.R")
  run_figure_entry(
    "data_Figure6.R",
    c(
      paste0("--n-core=", n_core),
      paste0(
        "--rebuild=",
        if (isTRUE(rebuild_figure6_grid)) "TRUE" else "FALSE"
      )
    )
  )
  run_figure_entry("build_model_dependent_manuscript_tables.R")
  run_figure_entry("draw_Figure6.R")
  run_figure_entry(
    "data_Supp_Figure6_1.R",
    c(
      paste0("--n-core=", n_core),
      paste0(
        "--rebuild=",
        if (isTRUE(rebuild_figure6_grid)) "TRUE" else "FALSE"
      )
    )
  )
  run_figure_entry("draw_Supp_Figure6_1.R")
  run_figure_entry("data_Supp_Figure6_2.R")
  run_figure_entry("draw_Supp_Figure6_2.R")
  run_figure_entry(
    "data_Supp_Figure6_3.R",
    c(
      paste0("--n-core=", n_core),
      paste0(
        "--rebuild=",
        if (isTRUE(rebuild_figure6_grid)) "TRUE" else "FALSE"
      )
    )
  )
  run_figure_entry("draw_Supp_Figure6_3.R")
  run_figure_entry(
    "data_Supp_Figure6_4.R",
    c(
      paste0("--n-core=", n_core),
      paste0(
        "--rebuild=",
        if (isTRUE(rebuild_figure6_grid)) "TRUE" else "FALSE"
      )
    )
  )
  run_figure_entry("draw_Supp_Figure6_4.R")
  run_figure_entry("build_figure6_manuscript_tables.R")

  expected <- file.path(OUTPUT_ROOT, c(
    paste0("assembled_fig", 1:6, ".png"),
    "supp_fig4-1_all18_cluster_prior_violins.png",
    "supp_fig4-2_invivo_optimizer_diagnostics.png",
    "supp_fig5-1_joint_parameter_stability.png",
    "supp_fig5-2_joint_fit_optimizer_diagnostics.png",
    "supp_fig6-1_response_class_diagnostics.png",
    "supp_fig6-2_joint_ensemble_robustness.png",
    "supp_fig6-3_weak_gap_regime_robustness.png",
    "supp_fig6-4_extended_invitro_o2_response.png"
  ))
  require_files(expected, "complete figure set")
  message("\nAll Figure 1-6 and parent-indexed supplementary outputs are complete.")
  invisible(expected)
}

if (sys.nframe() == 0L) {
  options <- run_all_parse_args()
  run_all_figures(
    n_core = options$n_core,
    recompute_fixed_o2 = options$recompute_fixed_o2,
    recompute_invivo_tsne = options$recompute_invivo_tsne,
    rebuild_figure6_grid = options$rebuild_figure6_grid,
    model_dependent_only = options$model_dependent_only,
    first_main_figure = options$first_main_figure,
    resume_after_figure5f_de = options$resume_after_figure5f_de
  )
}
