#!/usr/bin/env Rscript

# Generate the per-source-file README registry and its machine-readable TSV.
#
# Usage:
#   Rscript runner/documentation/generate_code_file_registry.R \
#     --workflow_root=/path/to/O2_supply_demand_MAP

parse_args <- function(argv) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[bits[[1L]]]] <- if (length(bits) > 1L) {
      paste(bits[-1L], collapse = "=")
    } else {
      "TRUE"
    }
  }
  out
}

script_dir <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_arg)) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[[1L]]),
      mustWork = FALSE
    )))
  }
  getwd()
}

read_lines <- function(path) {
  tryCatch(readLines(path, warn = FALSE), error = function(e) character())
}

first_descriptive_comment <- function(lines) {
  lines <- head(lines, 100L)
  text <- trimws(lines)
  is_comment <- grepl("^#", text) | grepl("^//", text)
  text <- text[is_comment]
  text <- trimws(sub("^(#+|//+)\\s*", "", text))
  text <- text[
    nzchar(text) &
      !grepl("^!|^SBATCH|^[-=]{3,}$", text) &
      !grepl(
        "^(usage|example|arguments?|inputs?|outputs?|optional|called by|test):?\\s*$",
        text,
        ignore.case = TRUE
      )
  ]
  if (!length(text)) return("")
  out <- text[[1L]]
  if (nchar(out) > 180L) out <- paste0(substr(out, 1L, 177L), "...")
  out
}

humanize <- function(path) {
  stem <- sub("\\.[^.]+$", "", basename(path))
  stem <- gsub("[-_]+", " ", stem)
  stem <- gsub("\\bviz\\b", "visualization", stem, ignore.case = TRUE)
  stem <- gsub("\\bo2\\b", "O2", stem, ignore.case = TRUE)
  paste0(toupper(substr(stem, 1L, 1L)), substring(stem, 2L))
}

infer_kind <- function(path, text) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("sh", "sub")) {
    if (grepl("sbatch|slurm", text, ignore.case = TRUE)) {
      return("HPC shell entrypoint")
    }
    return("shell entrypoint")
  }
  if (ext == "py") return("Python entrypoint")
  if (ext == "cpp") return("C++ model source")
  if (ext == "rmd") return("report template")
  if (grepl(
    "if\\s*\\([^\\n]*(sys\\.nframe|identical\\(environment\\(\\),\\s*globalenv)",
    text,
    perl = TRUE
  )) {
    return("R entrypoint")
  }
  "R library/module"
}

purpose_override <- function(relative_path) {
  overrides <- c(
    "model/model_O2_supply_demand_MAP.cpp" =
      "Implements the compiled O2 supply-demand state transitions and numerical kernels.",
    "model/model_O2_supply_demand_MAP.R" =
      "Loads and validates the protected R/Rcpp O2 supply-demand model interface.",
    "optimizer/auto_calibrate_boundary_params.R" =
      "Calibrates fitting parameter boundaries from completed fit results.",
    "optimizer/fit_model_O2_supply_demand_MAP.R" =
      "Protected fitting dispatcher for in-vivo, in-vitro, and joint modes.",
    "optimizer/profile_likelihood_O2_supply_demand_MAP.R" =
      "Runs the protected profile-likelihood fitting workflow."
  )
  unname(overrides[relative_path])
}

infer_purpose <- function(path, relative_path, comment) {
  override <- purpose_override(relative_path)
  if (length(override) && !is.na(override) && nzchar(override)) return(override)
  if (nzchar(comment)) return(comment)
  b <- basename(path)
  action <- if (grepl("^(generate|simulate|materialize|run_fixed)", b)) {
    "Materializes fitted-parameter simulation products for"
  } else if (grepl("^(plot|viz|visualize)", b)) {
    "Draws figures from materialized tables for"
  } else if (grepl("^(render|.*_report)", b)) {
    "Assembles a report from existing tables and figures for"
  } else if (grepl(
    "^(analyze|analyse|compare|collect|select|calculate|build)",
    b
  )) {
    "Computes analysis tables for"
  } else if (grepl("^(submit|run_.*array)", b)) {
    "Launches or coordinates HPC work for"
  } else if (grepl("(util|helper|shared)", b)) {
    "Provides reusable functions for"
  } else {
    "Implements"
  }
  paste(action, humanize(path))
}

infer_inputs <- function(text, layer) {
  vals <- character()
  if (grepl("best_params", text, fixed = TRUE)) vals <- c(vals, "best parameters")
  if (grepl("fit_result", text, fixed = TRUE)) vals <- c(vals, "fit result")
  if (grepl("fit_summary", text, fixed = TRUE)) vals <- c(vals, "fit summary")
  if (grepl("parameter_table", text, fixed = TRUE)) vals <- c(vals, "parameter table")
  if (grepl("simulation_manifest|simulation_dir", text, perl = TRUE)) {
    vals <- c(vals, "simulation tables/manifest")
  }
  if (grepl("analysis_manifest|analysis_dir", text, perl = TRUE)) {
    vals <- c(vals, "analysis tables/manifest")
  }
  if (grepl("\\.ya?ml|config", text, ignore.case = TRUE)) {
    vals <- c(vals, "configuration")
  }
  if (!length(vals)) {
    vals <- switch(
      layer,
      model = "model state/configuration",
      optimizer = "fit configuration and observations",
      hpc = "CLI/environment settings",
      runner = "CLI paths and workflow settings",
      util = "caller-supplied R objects",
      "materialized files or caller-supplied R objects"
    )
  }
  paste(unique(vals), collapse = "; ")
}

infer_outputs <- function(text, layer) {
  vals <- character()
  if (grepl(
    "write\\.table|write\\.csv|write_tsv|saveRDS",
    text,
    perl = TRUE
  )) {
    vals <- c(vals, "tables/data files")
  }
  if (grepl(
    "ggsave|grDevices::pdf|\\bpdf\\s*\\(|\\bpng\\s*\\(",
    text,
    perl = TRUE
  )) {
    vals <- c(vals, "figures")
  }
  if (grepl("rmarkdown::render|\\.html", text, perl = TRUE)) {
    vals <- c(vals, "report/HTML")
  }
  if (grepl("\\bsbatch\\b|#SBATCH", text, perl = TRUE)) {
    vals <- c(vals, "Slurm jobs/logs")
  }
  if (!length(vals)) {
    vals <- switch(
      layer,
      model = "in-memory model functions/state transitions",
      optimizer = "fit artifacts",
      util = "R values/functions; no workflow output when sourced",
      "in-memory values or layer-specific artifacts"
    )
  }
  paste(unique(vals), collapse = "; ")
}

escape_md <- function(x) {
  x <- gsub("\\|", "\\\\|", x)
  gsub("\n", " ", x, fixed = TRUE)
}

display_functions <- function(x, n, max_n = 8L) {
  if (!n) return("none")
  vals <- strsplit(x, ", ", fixed = TRUE)[[1L]]
  shown <- head(vals, max_n)
  suffix <- if (length(vals) > max_n) {
    paste0(" (+", length(vals) - max_n, " more)")
  } else {
    ""
  }
  paste0(paste0("`", shown, "`", collapse = ", "), suffix)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  default_root <- normalizePath(
    file.path(script_dir(), "..", ".."),
    mustWork = FALSE
  )
  workflow_root <- normalizePath(
    argv$workflow_root %||% default_root,
    mustWork = TRUE
  )
  out_md <- normalizePath(
    argv$out_md %||% file.path(workflow_root, "docs", "CODE_FILE_REGISTRY.md"),
    mustWork = FALSE
  )
  out_tsv <- normalizePath(
    argv$out_tsv %||% file.path(workflow_root, "docs", "code_file_registry.tsv"),
    mustWork = FALSE
  )

  all_files <- list.files(workflow_root, recursive = TRUE, full.names = TRUE)
  all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
  all_files <- all_files[
    grepl("(\\.R|\\.Rmd|\\.py|\\.sh|\\.sub|\\.cpp)$", all_files) &
      !grepl("/\\.rcpp_cache", all_files)
  ]
  all_files <- sort(normalizePath(all_files, mustWork = TRUE))
  relative_paths <- substring(all_files, nchar(workflow_root) + 2L)

  test_root <- normalizePath(
    file.path(workflow_root, "..", "..", "tests"),
    mustWork = FALSE
  )
  test_files <- if (dir.exists(test_root)) {
    list.files(
      test_root,
      recursive = TRUE,
      full.names = TRUE,
      pattern = "\\.[Rr]$"
    )
  } else {
    character()
  }
  test_text <- lapply(test_files, read_lines)

  rows <- lapply(seq_along(all_files), function(i) {
    path <- all_files[[i]]
    relative_path <- relative_paths[[i]]
    lines <- read_lines(path)
    text <- paste(lines, collapse = "\n")
    layer <- strsplit(relative_path, "/", fixed = TRUE)[[1L]][[1L]]
    function_lines <- lines[grepl(
      "^\\s*[.A-Za-z][A-Za-z0-9._]*\\s*<-\\s*function",
      lines
    )]
    functions <- sub(
      "^\\s*([.A-Za-z][A-Za-z0-9._]*)\\s*<-\\s*function.*$",
      "\\1",
      function_lines
    )
    refs <- if (length(test_files)) {
      hit <- vapply(
        test_text,
        function(x) any(grepl(basename(path), x, fixed = TRUE)),
        logical(1)
      )
      basename(test_files[hit])
    } else {
      character()
    }
    data.frame(
      file = relative_path,
      layer = layer,
      kind = infer_kind(path, text),
      purpose = infer_purpose(
        path,
        relative_path,
        first_descriptive_comment(lines)
      ),
      inputs = infer_inputs(text, layer),
      outputs = infer_outputs(text, layer),
      function_count = length(functions),
      functions = paste(functions, collapse = ", "),
      tests = paste(refs, collapse = ", "),
      lines = length(lines),
      stringsAsFactors = FALSE
    )
  })
  registry <- do.call(rbind, rows)

  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    registry,
    out_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  layer_order <- c(
    "model",
    "optimizer",
    "util",
    "simulation",
    "analysis",
    "vis",
    "report",
    "runner",
    "hpc",
    "docs"
  )
  layers <- c(
    intersect(layer_order, unique(registry$layer)),
    setdiff(unique(registry$layer), layer_order)
  )
  markdown <- c(
    "# O2 supply-demand MAP per-file registry",
    "",
    paste0(
      "This generated registry documents every R, R Markdown, Python, shell, ",
      "Slurm, and C++ source file under `",
      basename(workflow_root),
      "`. Generated caches are excluded. Each row records its owner, ",
      "responsibility, input/output contract, top-level function surface, and ",
      "direct test references."
    ),
    "",
    paste0("Total documented source files: **", nrow(registry), "**."),
    ""
  )
  for (layer in layers) {
    tab <- registry[registry$layer == layer, , drop = FALSE]
    markdown <- c(
      markdown,
      paste0("## ", layer),
      "",
      "| File | Type and responsibility | Inputs -> outputs | Functions | Direct tests |",
      "|---|---|---|---|---|"
    )
    for (i in seq_len(nrow(tab))) {
      tests <- if (nzchar(tab$tests[[i]])) {
        paste0("`", gsub(", ", "`, `", tab$tests[[i]], fixed = TRUE), "`")
      } else {
        "covered by parse/static/full regression"
      }
      markdown <- c(
        markdown,
        paste0(
          "| `",
          escape_md(tab$file[[i]]),
          "` | ",
          escape_md(tab$kind[[i]]),
          ". ",
          escape_md(tab$purpose[[i]]),
          " | ",
          escape_md(tab$inputs[[i]]),
          " -> ",
          escape_md(tab$outputs[[i]]),
          " | ",
          display_functions(tab$functions[[i]], tab$function_count[[i]]),
          " | ",
          tests,
          " |"
        )
      )
    }
    markdown <- c(markdown, "")
  }
  writeLines(markdown, out_md, useBytes = TRUE)

  cat("FILES=", nrow(registry), "\n", sep = "")
  cat("MARKDOWN=", normalizePath(out_md, mustWork = TRUE), "\n", sep = "")
  cat("TSV=", normalizePath(out_tsv, mustWork = TRUE), "\n", sep = "")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

if (sys.nframe() == 0L) {
  main()
}
