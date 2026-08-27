#!/usr/bin/env Rscript

.postfit_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
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
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  getwd()
})

SCRIPT_DIR <- normalizePath(.postfit_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"),
  local = environment()
)
rm(.postfit_script_dir)

parse_args <- o2sd_parse_args
as_bool <- o2sd_as_bool
`%||%` <- o2sd_null_coalesce

postfit_command_text <- function(script, args) {
  paste(
    c(
      "Rscript",
      shQuote(normalizePath(script, mustWork = FALSE), type = "sh"),
      vapply(args, shQuote, character(1), type = "sh")
    ),
    collapse = " "
  )
}

run_postfit_step <- function(label,
                             script,
                             args,
                             log_path,
                             dry_run = FALSE) {
  if (!file.exists(script)) {
    stop("Missing ", label, " script: ", script, call. = FALSE)
  }
  command_text <- postfit_command_text(script, args)
  message("[postfit] ", label, ": ", command_text)
  if (isTRUE(dry_run)) return(invisible(command_text))

  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript",
    args = shQuote(c(normalizePath(script, mustWork = TRUE), args)),
    stdout = log_path,
    stderr = log_path
  )
  if (!identical(status, 0L)) {
    tail_text <- if (file.exists(log_path)) {
      paste(utils::tail(readLines(log_path, warn = FALSE), 40L), collapse = "\n")
    } else {
      ""
    }
    stop(
      label,
      " failed with exit status ",
      status,
      ". See ",
      log_path,
      if (nzchar(tail_text)) paste0("\nLast log lines:\n", tail_text) else "",
      call. = FALSE
    )
  }
  invisible(status)
}

append_optional_arg <- function(args, key, value) {
  if (is.null(value) || !length(value)) return(args)
  value <- trimws(as.character(value[[1]]))
  if (!nzchar(value)) return(args)
  c(args, paste0("--", key, "=", value))
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% stop(
    "Usage: run_postfit_pipeline.R --fit_dir=/abs/seed --scope=invivo|invitro|joint",
    call. = FALSE
  )
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  scope <- tolower(trimws(as.character(argv$scope %||% argv$mode %||% "")))
  if (!scope %in% c("invivo", "invitro", "joint")) {
    stop("--scope must be one of: invivo, invitro, joint.", call. = FALSE)
  }

  dry_run <- as_bool(argv$dry_run %||% FALSE, FALSE)
  run_simulation <- as_bool(argv$run_simulation %||% TRUE, TRUE)
  run_analysis <- as_bool(argv$run_analysis %||% TRUE, TRUE)
  run_vis <- as_bool(argv$run_vis %||% TRUE, TRUE)
  run_report <- as_bool(argv$run_report %||% TRUE, TRUE)
  data_dir_arg <- argv$data_dir
  if (!is.null(data_dir_arg) && length(data_dir_arg)) {
    data_dir_arg <- normalizePath(
      as.character(data_dir_arg[[1]]),
      mustWork = FALSE
    )
  }

  invivo_sim_dir <- normalizePath(
    argv$invivo_simulation_dir %||% file.path(fit_dir, "simulation", "invivo"),
    mustWork = FALSE
  )
  invitro_sim_dir <- normalizePath(
    argv$invitro_simulation_dir %||% file.path(fit_dir, "simulation", "invitro"),
    mustWork = FALSE
  )
  invitro_analysis_dir <- normalizePath(
    argv$invitro_analysis_dir %||% file.path(fit_dir, "analysis", "invitro"),
    mustWork = FALSE
  )
  joint_analysis_dir <- normalizePath(
    argv$joint_analysis_dir %||%
      file.path(fit_dir, "analysis", "joint_parameters"),
    mustWork = FALSE
  )
  invivo_viz_dir <- normalizePath(
    argv$invivo_viz_dir %||% file.path(fit_dir, "viz", "invivo"),
    mustWork = FALSE
  )
  invitro_viz_dir <- normalizePath(
    argv$invitro_viz_dir %||% file.path(fit_dir, "viz", "invitro"),
    mustWork = FALSE
  )

  scripts <- list(
    invivo_sim = file.path(
      WORKFLOW_ROOT,
      "simulation",
      "invivo",
      "generate_invivo_simulation_outputs.R"
    ),
    invitro_sim = file.path(
      WORKFLOW_ROOT,
      "simulation",
      "invitro",
      "run_invitro_simulation_outputs.R"
    ),
    invitro_analysis = file.path(
      WORKFLOW_ROOT,
      "analysis",
      "fit_diagnostics",
      "run_invitro_fit_diagnostics.R"
    ),
    joint_analysis = file.path(
      WORKFLOW_ROOT,
      "analysis",
      "fit_diagnostics",
      "run_joint_parameter_diagnostics.R"
    ),
    invivo_vis = file.path(
      WORKFLOW_ROOT,
      "vis",
      "viz_invivo_model_O2_supply_demand_MAP_results.R"
    ),
    invitro_vis = file.path(
      WORKFLOW_ROOT,
      "vis",
      "viz_invitro_model_O2_supply_demand_MAP_results.R"
    ),
    joint_vis = file.path(
      WORKFLOW_ROOT,
      "vis",
      "joint",
      "plot_invivo_invitro_comparison.R"
    ),
    joint_parameter_vis = file.path(
      WORKFLOW_ROOT,
      "vis",
      "joint",
      "plot_joint_parameter_ratios.R"
    ),
    fit_report = file.path(WORKFLOW_ROOT, "report", "render_fit_report.R"),
    invitro_report = file.path(
      WORKFLOW_ROOT,
      "report",
      "render_invitro_fit_report.R"
    )
  )

  wants_invivo <- scope %in% c("invivo", "joint")
  wants_invitro <- scope %in% c("invitro", "joint")

  if (run_simulation && wants_invivo) {
    args <- c(
      paste0("--fit_dir=", fit_dir),
      paste0("--simulation_dir=", invivo_sim_dir)
    )
    for (key in c(
      "data_dir",
      "report_dt",
      "predict_horizons",
      "predict_report_dt",
      "dose_zero_only",
      "truncate_at_treatment",
      "ploidy_at_harvest",
      "max_scenarios"
    )) {
      value <- if (identical(key, "data_dir")) data_dir_arg else argv[[key]]
      args <- append_optional_arg(args, key, value)
    }
    run_postfit_step(
      "in-vivo simulation materialization",
      scripts$invivo_sim,
      args,
      file.path(fit_dir, "invivo_simulation_status.log"),
      dry_run
    )
  }

  if (run_simulation && wants_invitro) {
    run_postfit_step(
      "in-vitro simulation materialization",
      scripts$invitro_sim,
      c(
        paste0("--fit_dir=", fit_dir),
        paste0("--simulation_dir=", invitro_sim_dir)
      ),
      file.path(fit_dir, "invitro_simulation_status.log"),
      dry_run
    )
  }

  if (run_analysis && wants_invitro) {
    run_postfit_step(
      "in-vitro fit diagnostics",
      scripts$invitro_analysis,
      c(
        paste0("--fit_dir=", fit_dir),
        paste0("--simulation_dir=", invitro_sim_dir),
        paste0("--analysis_dir=", invitro_analysis_dir)
      ),
      file.path(fit_dir, "invitro_analysis_status.log"),
      dry_run
    )
  }

  if (run_analysis && identical(scope, "joint")) {
    run_postfit_step(
      "joint parameter diagnostics",
      scripts$joint_analysis,
      c(
        paste0("--fit_dir=", fit_dir),
        paste0("--analysis_dir=", joint_analysis_dir)
      ),
      file.path(fit_dir, "joint_parameter_analysis_status.log"),
      dry_run
    )
  }

  if (run_vis && wants_invivo) {
    run_postfit_step(
      "in-vivo visualization",
      scripts$invivo_vis,
      c(
        paste0("--fit_dir=", fit_dir),
        paste0("--simulation_dir=", invivo_sim_dir),
        paste0("--out_dir=", invivo_viz_dir)
      ),
      file.path(fit_dir, "viz_status.log"),
      dry_run
    )
  }

  if (run_vis && wants_invitro) {
    run_postfit_step(
      "in-vitro visualization",
      scripts$invitro_vis,
      c(
        paste0("--fit_dir=", fit_dir),
        paste0("--simulation_dir=", invitro_sim_dir),
        paste0("--analysis_dir=", invitro_analysis_dir),
        paste0("--out_dir=", invitro_viz_dir)
      ),
      file.path(fit_dir, "invitro_viz_status.log"),
      dry_run
    )
  }

  if (run_vis && identical(scope, "joint")) {
    run_postfit_step(
      "joint in-vivo/in-vitro visualization",
      scripts$joint_vis,
      c(
        paste0("--fit_dir=", fit_dir),
        paste0("--invivo_simulation_dir=", invivo_sim_dir),
        paste0("--invitro_simulation_dir=", invitro_sim_dir),
        paste0(
          "--out_dir=",
          normalizePath(
            argv$joint_viz_dir %||% file.path(fit_dir, "viz", "invivo_vs_invitro"),
            mustWork = FALSE
          )
        )
      ),
      file.path(fit_dir, "joint_viz_status.log"),
      dry_run
    )
    run_postfit_step(
      "joint parameter visualization",
      scripts$joint_parameter_vis,
      c(
        paste0("--fit_dir=", fit_dir),
        paste0("--analysis_dir=", joint_analysis_dir),
        paste0(
          "--out_dir=",
          normalizePath(
            argv$joint_parameter_viz_dir %||%
              file.path(fit_dir, "viz", "joint_parameters"),
            mustWork = FALSE
          )
        )
      ),
      file.path(fit_dir, "joint_parameter_viz_status.log"),
      dry_run
    )
  }

  if (run_report) {
    report_script <- if (identical(scope, "invitro")) {
      scripts$invitro_report
    } else {
      scripts$fit_report
    }
    report_args <- paste0("--fit_dir=", fit_dir)
    for (key in c("render_pdf", "out_subdir", "report_basename")) {
      report_args <- append_optional_arg(report_args, key, argv[[key]])
    }
    run_postfit_step(
      "fit report",
      report_script,
      report_args,
      file.path(fit_dir, "report_status.log"),
      dry_run
    )
  }

  message("[postfit] completed scope=", scope, " fit_dir=", fit_dir)
  invisible(fit_dir)
}

if (sys.nframe() == 0L) {
  main()
}
