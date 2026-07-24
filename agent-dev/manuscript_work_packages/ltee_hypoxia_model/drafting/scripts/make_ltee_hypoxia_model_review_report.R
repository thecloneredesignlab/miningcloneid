#!/usr/bin/env Rscript

# Copy this file into a drafting package's scripts/ directory and adapt the
# CONFIG block. It writes a self-contained review_report.html with embedded PNGs.

# CONFIG -----------------------------------------------------------------------

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "scripts/make_ltee_hypoxia_model_review_report.R"
default_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
root <- Sys.getenv("REPORT_ROOT", unset = default_root)
out_path <- file.path(root, "review_report.html")
repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = "")

report_title <- "LTee Hypoxia Model — Gate B Manuscript-Figure Review"
reviewer_action <- paste(
  "Review the five recommended figures in order, focusing on claim support,",
  "panel density, and legibility at the displayed manuscript dimensions.",
  "The supplementary diagnostic and all intermediate drafts follow."
)

directive_table <- data.frame(
  directive_id = sprintf("FD%02d", 1:18),
  source = rep("feedback_manager_context.md; feedback_intake.md", 18),
  directive = c(
    "Add chromosome trajectories to the matched-design figure.",
    "Retain one integrated model overview.",
    "Limit Figure 3A to growth rate over lineage passage.",
    "Retain the historical control-versus-deprived comparison.",
    "Do not draft Figure 3E or run a restricted-model comparison.",
    "Lead Figures 3–5 with direct fit-quality evidence.",
    "Treat historical Figure 4A–E as in-vivo analyses.",
    "Keep the fixed-O2 triptych and solution landscape in the main figure.",
    "Use exact low/high-O2 parameter definitions and noncausal AUC language.",
    "Use the six July joint-pair winners throughout Figure 5.",
    "Order Figure 5 as fit quality, ratios, proliferation, missegregation, survival.",
    "Retain a full all-six parameter-ratio view.",
    "Leave Figure 6 out.",
    "Defer journal-specific constraints until polishing.",
    "Retain vi_C01, vi_C02, and vi_C03.",
    "Preserve the pooled in-vivo/in-vitro embedding and saved geometry.",
    "Keep main fit blocks compact; place optimizer diagnostics supplementary.",
    "Draft only within the approved five-figure/exclusion boundary."
  ),
  status = c(rep("addressed", 13), "deferred", rep("addressed", 4)),
  output = c(
    "final_figures/recommended/figure1.png",
    "final_figures/recommended/figure2.png",
    "final_figures/recommended/figure3.png",
    "final_figures/recommended/figure3.png",
    "not_drafted.md; no Figure 3E artifact",
    "final_figures/recommended/figure3.png; final_figures/recommended/figure4.png; final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions.png",
    "final_figures/recommended/figure4.png; subagent_notes/figure4.md",
    "final_figures/recommended/figure4.png",
    "final_figures/recommended/figure4.png; subagent_notes/figure4_legend.md",
    "final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions.png",
    "final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions.png",
    "final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions.png",
    "not_drafted.md; no Figure 6 artifact",
    "not_drafted.md; repository style contract applied",
    "final_figures/recommended/figure4.png",
    "final_figures/recommended/figure4.png; initial_subpanels/F4D/figure4D_pooled_embedding_preserved_source.png",
    "final_figures/recommended/figure3.png; final_figures/recommended/figure4.png; final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions.png; supplementary/fit_quality_optimizer_diagnostics.png",
    "final_figures/recommended/; not_drafted.md"
  ),
  report_section = c(
    "Figure 1", "Figure 2", "Figure 3", "Figure 3", "Scope Boundaries",
    "Figures 3–5", "Figure 4", "Figure 4", "Figure 4",
    "Figure 5", "Figure 5", "Figure 5", "Scope Boundaries",
    "Scope Boundaries", "Figure 4", "Figure 4",
    "Fit Diagnostics", "Scope Boundaries"
  ),
  stringsAsFactors = FALSE
)

cards <- list(
  list(
    section = "Recommended Final Drafts",
    title = "Figure 1 — Matched design and chromosome trajectories",
    badge = "recommended",
    image = file.path("final_figures", "recommended", "figure1.png"),
    legend = file.path("subagent_notes", "figure1_legend.md"),
    legend_heading = NULL,
    summary = "The study design is paired directly with the chromosome observations that motivate the modelling claims.",
    feedback = c(
      "Adds the requested in-vitro chromosome trajectories.",
      "Shows the in-vivo starting reference and terminal tumor distributions."
    ),
    caveats = c(
      "This is descriptive experimental context; no model fitting is shown."
    ),
    files = c(
      "scripts/make_figure1.R",
      "subagent_notes/figure1.md",
      "source_tables/figure1_displayed_summaries.tsv"
    ),
    directive_ids = c("FD01"),
    gallery = c(
      "refined_subpanels/figure1_recommended.png",
      "initial_subpanels/F1A/matched_design.png",
      "initial_subpanels/F1B/invitro_chromosome_trajectories.png",
      "initial_subpanels/F1C/invivo_start_terminal_distributions.png"
    ),
    gallery_open = FALSE,
    represented_images = character(),
    image_class = "landscape"
  ),
  list(
    section = "Recommended Final Drafts",
    title = "Figure 2 — Integrated model overview",
    badge = "recommended",
    image = "final_figures/recommended/figure2.png",
    legend = "subagent_notes/figure2_legend.md",
    legend_heading = NULL,
    summary = "One schematic connects resource state, fitted growth and death, chromosome-state transitions, survival filtering, and context-specific observables.",
    feedback = c("Retains one integrated model overview rather than splitting the mechanism across figures."),
    caveats = c("Arrows encode model structure and fitted links, not direct causal identification from the experiments."),
    files = c("scripts/make_figure2.R", "subagent_notes/figure2.md"),
    directive_ids = c("FD02"),
    gallery = c(
      "refined_subpanels/figure2_recommended.png",
      "initial_subpanels/F2/model_overview.png"
    ),
    gallery_open = FALSE,
    represented_images = character(),
    image_class = "landscape"
  ),
  list(
    section = "Recommended Final Drafts",
    title = "Figure 3 — In-vitro fit adequacy and fitted functions",
    badge = "recommended",
    image = "final_figures/recommended/figure3.png",
    legend = "subagent_notes/figure3_legend.md",
    legend_heading = NULL,
    summary = "Direct growth-fit evidence leads the figure; controls remain visible before conditional fitted survival and nonviable-fraction summaries.",
    feedback = c(
      "Panel A is restricted to growth over passage.",
      "The left-hand control conditions are retained.",
      "The rejected restricted-model panel is absent."
    ),
    caveats = c("Fit comparisons are in-sample summaries of saved outputs; they are not held-out validation."),
    files = c("scripts/make_figure3.R", "subagent_notes/figure3.md"),
    directive_ids = c("FD03", "FD04", "FD05", "FD06", "FD17", "FD18"),
    gallery = c(
      "refined_subpanels/figure3_refined.png",
      "initial_subpanels/F3A/figure3A_initial.png",
      "initial_subpanels/F3B/figure3B_initial.png",
      "initial_subpanels/F3C/figure3C_initial.png",
      "initial_subpanels/F3D/figure3D_initial.png"
    ),
    gallery_open = FALSE,
    represented_images = character(),
    image_class = "portrait"
  ),
  list(
    section = "Recommended Final Drafts",
    title = "Figure 4 — In-vivo fit adequacy and solution landscape",
    badge = "recommended",
    image = "final_figures/recommended/figure4.png",
    legend = "subagent_notes/figure4_legend.md",
    legend_heading = NULL,
    summary = "Fit adequacy precedes latent oxygen, fixed-O2 separators, the unchanged pooled embedding, and parameter distributions for all three in-vivo regions.",
    feedback = c(
      "The pooled in-vivo/in-vitro embedding and saved geometry are preserved.",
      "All vi_C01, vi_C02, and vi_C03 regions remain visible.",
      "AUCs are labelled as univariate separators, not causal effects."
    ),
    caveats = c(
      "Effective O2 is a fitted latent state, not a tumor oxygen measurement.",
      "Necrosis is excluded because the prediction export remains unavailable."
    ),
    files = c("scripts/make_figure4.R", "subagent_notes/figure4.md"),
    directive_ids = c("FD06", "FD07", "FD08", "FD09", "FD15", "FD16", "FD17"),
    gallery = c(
      "refined_subpanels/figure4_refined.png",
      "initial_subpanels/F4A/figure4A_fit_adequacy.png",
      "initial_subpanels/F4B/figure4B_latent_effective_o2.png",
      "initial_subpanels/F4C/figure4C_univariate_separators.png",
      "initial_subpanels/F4D/figure4D_pooled_embedding.png",
      "initial_subpanels/F4D/figure4D_pooled_embedding_layout_source.png",
      "initial_subpanels/F4D/figure4D_pooled_embedding_preserved_source.png",
      "initial_subpanels/F4E/figure4E_region_parameters.png"
    ),
    gallery_open = FALSE,
    represented_images = character(),
    image_class = "portrait"
  ),
  list(
    section = "Recommended Final Drafts",
    title = "Figure 5 — Joint-fit adequacy and context functions",
    badge = "recommended",
    image = "final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions.png",
    legend = "subagent_notes/figure5_legend.md",
    legend_heading = NULL,
    summary = "The common six-winner joint-fit universe is used consistently from direct fit checks through ratios and context-specific functions.",
    feedback = c(
      "Fit quality is shown before biological interpretation.",
      "All parameter ratios and all six joint winners are retained.",
      "The approved proliferation, missegregation, and survival order is used."
    ),
    caveats = c("Displayed curves summarize saved fitted parameterizations; no optimizer refit was run."),
    files = c("scripts/make_figure5.R", "subagent_notes/figure5.md"),
    directive_ids = c("FD06", "FD10", "FD11", "FD12", "FD17", "FD18"),
    gallery = c(
      "refined_subpanels/figure5_joint_fit_adequacy_and_context_functions.png",
      "initial_subpanels/F5A/figure5A_fit_adequacy.png",
      "initial_subpanels/F5B/figure5B_parameter_ratios.png",
      "initial_subpanels/F5C/figure5C_proliferation.png",
      "initial_subpanels/F5D/figure5D_missegregation.png",
      "initial_subpanels/F5E/figure5E_survival.png"
    ),
    gallery_open = FALSE,
    represented_images = character(),
    image_class = "portrait"
  ),
  list(
    section = "Other Review Outputs",
    title = "Supplementary fit and optimizer diagnostics",
    badge = "supplementary",
    image = "supplementary/fit_quality_optimizer_diagnostics.png",
    legend = "subagent_notes/fit_diagnostics_legend.md",
    legend_heading = NULL,
    summary = "Saved objective distributions and optimizer-status diagnostics support the compact direct-fit blocks without crowding the main figures.",
    feedback = c("Implements the approved separation of compact main fit blocks and supplementary optimizer diagnostics."),
    caveats = c("Diagnostics report saved optimization records as-is and do not constitute new optimization."),
    files = c(
      "scripts/make_fit_diagnostics.R",
      "supplementary/fit_quality_optimizer_diagnostics.tsv",
      "subagent_notes/fit_diagnostics.md"
    ),
    directive_ids = c("FD17"),
    gallery = character(),
    gallery_open = FALSE,
    represented_images = character(),
    image_class = "portrait"
  ),
  list(
    section = "Other Review Outputs",
    title = "Approved scope boundaries",
    badge = "scope",
    image = character(),
    legend = "not_drafted.md",
    legend_heading = NULL,
    summary = "The package stops at the approved Gate B boundary.",
    feedback = c(
      "No Figure 3E or restricted-model comparison was drafted.",
      "Figure 6 remains deferred.",
      "No optimizer refits or production manuscript edits were performed."
    ),
    caveats = c("Journal-specific polishing remains a later gate."),
    files = c("not_drafted.md", "feedback_intake.md"),
    directive_ids = c("FD05", "FD13", "FD14", "FD18"),
    gallery = character(),
    gallery_open = FALSE,
    represented_images = character(),
    image_class = "compact"
  )
)

section_order <- c(
  "Recommended Final Drafts",
  "Blockers And Missing Required Outputs",
  "Conservative Alternatives",
  "Exploratory Alternatives",
  "Refined And Initial Subpanels",
  "Other Review Outputs"
)
directive_status_levels <- c("addressed", "partially_addressed", "blocked", "dropped", "deferred")
prior_code_fidelity_path <- file.path(root, "prior_code_fidelity.csv")
prior_code_required_columns <- c(
  "panel_id",
  "inheritance_class",
  "prior_png_path",
  "prior_code_path",
  "copied_baseline_code_path",
  "active_local_script_path",
  "diff_path",
  "allowed_change_directive_ids",
  "fidelity_status",
  "blocker"
)
prior_code_inherited_classes <- c("inherited_preserve", "inherited_targeted_fix", "inherited_move")
prior_code_allowed_classes <- c(
  prior_code_inherited_classes,
  "inherited_replace",
  "novel_no_prior",
  "blocked_prior_missing"
)
png_scan_dirs <- c("initial_subpanels", "refined_subpanels", "final_figures", "supplementary")
exclude_from_report <- character()
report_manifest_path <- file.path(root, "report_manifest.csv")

appendix_files <- c(
  "feedback_manager_context.md",
  "feedback_intake.md",
  "report_manifest.csv",
  "source_manifest.csv",
  "source_tables/frozen_input_manifest.csv",
  "source_context/README.md",
  "final_figures/recommended/legend.md",
  "validation_report.md",
  "prior_code_fidelity.csv",
  "drafting_panels.md",
  "prior_panel_disposition.csv",
  "not_drafted.md"
)

# HELPERS ----------------------------------------------------------------------

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

read_text <- function(path) {
  if (!file.exists(path)) return("")
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

png_data_uri <- function(path) {
  if (!file.exists(path)) stop("Missing PNG for report: ", path)
  encoded <- system2("base64", c("-i", path), stdout = TRUE)
  paste0("data:image/png;base64,", paste(encoded, collapse = ""))
}

extract_md_section <- function(path, heading_fragment = NULL) {
  text <- read_text(path)
  if (!nzchar(text) || is.null(heading_fragment) || !nzchar(heading_fragment)) {
    return(text)
  }
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  start <- grep(paste0("^## .*", heading_fragment), lines)
  if (!length(start)) return(text)
  start <- start[[1]]
  next_heading <- grep("^## ", lines)
  next_heading <- next_heading[next_heading > start]
  end <- if (length(next_heading)) next_heading[[1]] - 1L else length(lines)
  paste(lines[start:end], collapse = "\n")
}

paragraphs_html <- function(text) {
  if (!nzchar(text)) return("<p>No legend text found.</p>")
  text <- gsub("^## +", "", text)
  blocks <- unlist(strsplit(text, "\n[ \t]*\n"))
  blocks <- blocks[nzchar(trimws(blocks))]
  inline_html <- function(x) {
    x <- html_escape(x)
    x <- gsub("`([^`]+)`", "<code>\\1</code>", x)
    gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", x)
  }
  paste0(
    vapply(blocks, function(block) {
      block <- trimws(block)
      if (grepl("^# +", block)) {
        return(paste0("<h4>", inline_html(sub("^# +", "", block)), "</h4>"))
      }
      block <- gsub("\n", " ", block)
      paste0("<p>", inline_html(block), "</p>")
    }, character(1)),
    collapse = "\n"
  )
}

list_html <- function(items) {
  items <- items[nzchar(items)]
  if (!length(items)) return("<p>None recorded.</p>")
  paste0("<ul>", paste0("<li>", html_escape(items), "</li>", collapse = ""), "</ul>")
}

or_default <- function(x, default) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(x[[1]])) default else x
}

is_true <- function(x) {
  !is.null(x) && length(x) && isTRUE(x[[1]])
}

collapse_nonempty <- function(x, sep = "; ") {
  x <- unique(x[nzchar(x)])
  if (!length(x)) "" else paste(x, collapse = sep)
}

norm_rel <- function(path) {
  path <- gsub("\\\\", "/", path)
  path <- sub("^\\./", "", path)
  path
}

status_class <- function(status) {
  paste0("status-", gsub("[^a-z0-9]+", "-", tolower(status)))
}

validate_directive_table <- function(directives) {
  required <- c("directive_id", "source", "directive", "status", "output", "report_section")
  missing <- setdiff(required, names(directives))
  if (length(missing)) {
    stop("directive_table is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  directives[] <- lapply(directives, as.character)
  for (field in required) {
    bad <- !nzchar(trimws(directives[[field]]))
    if (any(bad)) {
      stop(
        "directive_table has blank values in column ", field, " for row(s): ",
        paste(which(bad), collapse = ", "),
        call. = FALSE
      )
    }
  }
  duplicate_ids <- directives$directive_id[duplicated(directives$directive_id)]
  if (length(duplicate_ids)) {
    stop("directive_table has duplicate directive_id values: ", paste(unique(duplicate_ids), collapse = ", "), call. = FALSE)
  }
  invalid_status <- setdiff(unique(directives$status), directive_status_levels)
  if (length(invalid_status)) {
    stop(
      "directive_table has invalid status values: ", paste(invalid_status, collapse = ", "),
      ". Use one of: ", paste(directive_status_levels, collapse = ", "),
      call. = FALSE
    )
  }
  directives
}

split_recorded_paths <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(character())
  paths <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  trimws(paths[nzchar(trimws(paths))])
}

path_exists_from_root <- function(path) {
  if (!nzchar(path)) return(FALSE)
  if (grepl("^/", path)) {
    return(file.exists(path))
  }
  if (file.exists(file.path(root, path))) return(TRUE)
  nzchar(repo_root) && file.exists(file.path(repo_root, path))
}

all_recorded_paths_exist <- function(value) {
  paths <- split_recorded_paths(value)
  length(paths) > 0 && all(vapply(paths, path_exists_from_root, logical(1)))
}

validate_prior_code_fidelity <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing required prior_code_fidelity.csv. Create the prior-code fidelity manifest ",
      "before generating review_report.html.",
      call. = FALSE
    )
  }
  fidelity <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(prior_code_required_columns, names(fidelity))
  if (length(missing)) {
    stop(
      "prior_code_fidelity.csv is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  fidelity[prior_code_required_columns] <- lapply(fidelity[prior_code_required_columns], as.character)
  required_nonblank <- c("panel_id", "inheritance_class", "fidelity_status")
  for (field in required_nonblank) {
    bad <- !nzchar(trimws(fidelity[[field]]))
    if (any(bad)) {
      stop(
        "prior_code_fidelity.csv has blank values in column ", field, " for row(s): ",
        paste(which(bad), collapse = ", "),
        call. = FALSE
      )
    }
  }
  duplicate_ids <- fidelity$panel_id[duplicated(fidelity$panel_id)]
  if (length(duplicate_ids)) {
    stop(
      "prior_code_fidelity.csv has duplicate panel_id values: ",
      paste(unique(duplicate_ids), collapse = ", "),
      call. = FALSE
    )
  }
  invalid_class <- setdiff(unique(fidelity$inheritance_class), prior_code_allowed_classes)
  if (length(invalid_class)) {
    stop(
      "prior_code_fidelity.csv has invalid inheritance_class values: ",
      paste(invalid_class, collapse = ", "),
      ". Use one of: ", paste(prior_code_allowed_classes, collapse = ", "),
      call. = FALSE
    )
  }

  inherited <- fidelity$inheritance_class %in% prior_code_inherited_classes
  inherited_fields <- c(
    "prior_png_path",
    "prior_code_path",
    "copied_baseline_code_path",
    "active_local_script_path",
    "diff_path"
  )
  for (field in inherited_fields) {
    bad <- inherited & !nzchar(trimws(fidelity[[field]]))
    if (any(bad)) {
      stop(
        "Inherited rows in prior_code_fidelity.csv need non-empty ", field,
        " values. Bad row(s): ", paste(which(bad), collapse = ", "),
        call. = FALSE
      )
    }
  }
  for (field in c("copied_baseline_code_path", "active_local_script_path", "diff_path")) {
    bad <- inherited & !vapply(fidelity[[field]], all_recorded_paths_exist, logical(1))
    if (any(bad)) {
      stop(
        "Inherited rows in prior_code_fidelity.csv reference missing local ", field,
        " path(s). Bad row(s): ", paste(which(bad), collapse = ", "),
        call. = FALSE
      )
    }
  }

  blocked_missing <- fidelity$inheritance_class == "blocked_prior_missing"
  bad_blocker <- blocked_missing & !nzchar(trimws(fidelity$blocker))
  if (any(bad_blocker)) {
    stop(
      "blocked_prior_missing rows in prior_code_fidelity.csv need a blocker. Bad row(s): ",
      paste(which(bad_blocker), collapse = ", "),
      call. = FALSE
    )
  }

  replaced <- fidelity$inheritance_class == "inherited_replace"
  bad_replace <- replaced &
    !nzchar(trimws(fidelity$allowed_change_directive_ids)) &
    !nzchar(trimws(fidelity$blocker))
  if (any(bad_replace)) {
    stop(
      "inherited_replace rows in prior_code_fidelity.csv need directive ids or a blocker. Bad row(s): ",
      paste(which(bad_replace), collapse = ", "),
      call. = FALSE
    )
  }

  fidelity
}

status_counts_html <- function(directives) {
  counts <- table(factor(directives$status, levels = directive_status_levels))
  items <- vapply(names(counts), function(status) {
    paste0(
      "<span class=\"status-pill ", status_class(status), "\">",
      html_escape(status), ": ", counts[[status]], "</span>"
    )
  }, character(1))
  paste0("<div class=\"status-row\">", paste(items, collapse = ""), "</div>")
}

directive_table_html <- function(directives) {
  rows <- vapply(seq_len(nrow(directives)), function(i) {
    row <- directives[i, , drop = FALSE]
    paste0(
      "<tr>",
      "<td><code>", html_escape(row$directive_id), "</code></td>",
      "<td>", html_escape(row$source), "</td>",
      "<td>", html_escape(row$directive), "</td>",
      "<td><span class=\"status-pill ", status_class(row$status), "\">", html_escape(row$status), "</span></td>",
      "<td>", html_escape(row$output), "</td>",
      "<td>", html_escape(row$report_section), "</td>",
      "</tr>"
    )
  }, character(1))
  paste0(
    "<section class=\"directive-block\"><h2>Feedback Coverage</h2>",
    status_counts_html(directives),
    "<div class=\"table-wrap\"><table><thead><tr>",
    "<th>ID</th><th>Source</th><th>Directive</th><th>Status</th><th>Output Or Disposition</th><th>First Report Section</th>",
    "</tr></thead><tbody>",
    paste(rows, collapse = "\n"),
    "</tbody></table></div></section>"
  )
}

prior_code_fidelity_html <- function(fidelity) {
  cols <- c(
    "panel_id",
    "inheritance_class",
    "prior_png_path",
    "prior_code_path",
    "copied_baseline_code_path",
    "active_local_script_path",
    "diff_path",
    "allowed_change_directive_ids",
    "fidelity_status",
    "blocker"
  )
  rows <- vapply(seq_len(nrow(fidelity)), function(i) {
    row <- fidelity[i, cols, drop = FALSE]
    paste0(
      "<tr>",
      "<td><code>", html_escape(row$panel_id), "</code></td>",
      "<td>", html_escape(row$inheritance_class), "</td>",
      "<td>", html_escape(row$prior_png_path), "</td>",
      "<td>", html_escape(row$prior_code_path), "</td>",
      "<td>", html_escape(row$copied_baseline_code_path), "</td>",
      "<td>", html_escape(row$active_local_script_path), "</td>",
      "<td>", html_escape(row$diff_path), "</td>",
      "<td>", html_escape(row$allowed_change_directive_ids), "</td>",
      "<td><span class=\"status-pill ", status_class(row$fidelity_status), "\">",
      html_escape(row$fidelity_status), "</span></td>",
      "<td>", html_escape(row$blocker), "</td>",
      "</tr>"
    )
  }, character(1))
  paste0(
    "<section class=\"directive-block\"><h2>Prior-Code Fidelity</h2>",
    "<p>Inherited panels must start from copied prior reviewed code, with local active scripts and diffs visible here.</p>",
    "<div class=\"table-wrap\"><table><thead><tr>",
    "<th>Panel</th><th>Class</th><th>Prior PNG</th><th>Prior Code</th>",
    "<th>Copied Baseline</th><th>Active Local Script</th><th>Diff</th>",
    "<th>Allowed Directives</th><th>Status</th><th>Blocker</th>",
    "</tr></thead><tbody>",
    paste(rows, collapse = "\n"),
    "</tbody></table></div></section>"
  )
}

directive_subset_html <- function(ids, directives) {
  ids <- unique(ids[nzchar(ids)])
  if (!length(ids)) return("")
  matched <- directives[directives$directive_id %in% ids, , drop = FALSE]
  if (!nrow(matched)) {
    return(paste0("<p>Unknown directive id(s): <code>", html_escape(paste(ids, collapse = ", ")), "</code></p>"))
  }
  rows <- vapply(seq_len(nrow(matched)), function(i) {
    row <- matched[i, , drop = FALSE]
    paste0(
      "<li><code>", html_escape(row$directive_id), "</code> ",
      "<span class=\"status-pill ", status_class(row$status), "\">", html_escape(row$status), "</span> ",
      html_escape(row$directive), "</li>"
    )
  }, character(1))
  paste0("<ul class=\"directive-mini\">", paste(rows, collapse = ""), "</ul>")
}

card_section <- function(card) {
  if (!is.null(card$section) && length(card$section) && nzchar(card$section[[1]])) {
    return(card$section[[1]])
  }
  badge <- if (!is.null(card$badge) && length(card$badge)) tolower(card$badge[[1]]) else ""
  if (identical(badge, "recommended")) return("Recommended Final Drafts")
  if (identical(badge, "conservative")) return("Conservative Alternatives")
  if (identical(badge, "exploratory")) return("Exploratory Alternatives")
  "Other Review Outputs"
}

card_directive_ids <- function(card) {
  if (is.null(card$directive_ids)) character() else as.character(card$directive_ids)
}

all_draft_pngs <- function(dirs) {
  paths <- unlist(lapply(dirs, function(dir) {
    abs_dir <- file.path(root, dir)
    if (!dir.exists(abs_dir)) return(character())
    file.path(dir, list.files(abs_dir, pattern = "\\.png$", recursive = TRUE, full.names = FALSE))
  }), use.names = FALSE)
  sort(unique(norm_rel(paths)))
}

card_display_rows <- function(card) {
  primary <- if (!is.null(card$image)) norm_rel(card$image) else character()
  gallery <- if (!is.null(card$gallery)) norm_rel(card$gallery) else character()
  represented <- if (!is.null(card$represented_images)) norm_rel(card$represented_images) else character()
  paths <- c(primary, gallery, represented)
  paths <- paths[nzchar(paths)]
  directive_ids <- collapse_nonempty(card_directive_ids(card))
  modes <- c(
    rep("primary", length(primary)),
    rep("gallery", length(gallery)),
    rep("represented_by_contact_sheet", length(represented))
  )
  modes <- modes[nzchar(c(primary, gallery, represented))]
  if (!length(paths)) {
    return(data.frame(relative_path = character(), report_section = character(), display_mode = character(), directive_id = character()))
  }
  data.frame(
    relative_path = paths,
    report_section = card$title,
    display_mode = modes,
    directive_id = rep(directive_ids, length(paths)),
    stringsAsFactors = FALSE
  )
}

display_manifest <- function(cards) {
  rows <- do.call(rbind, lapply(cards, card_display_rows))
  if (is.null(rows) || !nrow(rows)) {
    return(data.frame(relative_path = character(), report_section = character(), display_mode = character(), directive_id = character()))
  }
  rows[!duplicated(rows$relative_path), , drop = FALSE]
}

infer_kind <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (!length(parts)) return("")
  parts[[1]]
}

infer_variant <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (length(parts) >= 2 && identical(parts[[1]], "final_figures")) return(parts[[2]])
  ""
}

infer_figure_or_subpanel <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (length(parts) >= 2 && parts[[1]] %in% c("initial_subpanels", "refined_subpanels")) return(parts[[2]])
  tools::file_path_sans_ext(basename(path))
}

write_coverage_manifest <- function(draft_pngs, display_rows, exclude_paths) {
  report_section <- rep("", length(draft_pngs))
  display_mode <- rep("", length(draft_pngs))
  directive_id <- rep("", length(draft_pngs))
  matched <- match(draft_pngs, display_rows$relative_path)
  has_match <- !is.na(matched)
  report_section[has_match] <- display_rows$report_section[matched[has_match]]
  display_mode[has_match] <- display_rows$display_mode[matched[has_match]]
  directive_id[has_match] <- display_rows$directive_id[matched[has_match]]
  include_status <- ifelse(draft_pngs %in% exclude_paths, "not_reviewable_internal", ifelse(has_match, "included", "missing"))
  rationale <- ifelse(include_status == "not_reviewable_internal", "Excluded by template CONFIG.", "")
  manifest <- data.frame(
    relative_path = draft_pngs,
    kind = vapply(draft_pngs, infer_kind, character(1)),
    figure_or_subpanel = vapply(draft_pngs, infer_figure_or_subpanel, character(1)),
    variant = vapply(draft_pngs, infer_variant, character(1)),
    directive_id = directive_id,
    feedback_basis = "",
    report_section = report_section,
    display_mode = display_mode,
    include_status = include_status,
    rationale = rationale,
    stringsAsFactors = FALSE
  )
  write.csv(manifest, report_manifest_path, row.names = FALSE, quote = TRUE)
  manifest
}

gallery_html <- function(paths, open = FALSE, label = "Additional Draft Outputs") {
  paths <- norm_rel(paths)
  paths <- paths[nzchar(paths)]
  if (!length(paths)) return("")
  items <- vapply(paths, function(path) {
    image_uri <- png_data_uri(file.path(root, path))
    paste0(
      "<figure class=\"gallery-item\"><img src=\"", image_uri, "\" alt=\"", html_escape(path), "\">",
      "<figcaption>", html_escape(path), "</figcaption></figure>"
    )
  }, character(1))
  paste0(
    "<details class=\"gallery-block\"", if (isTRUE(open)) " open" else "", "><summary><strong>",
    html_escape(label), "</strong></summary>",
    "<div class=\"gallery-grid\">", paste(items, collapse = "\n"), "</div></details>"
  )
}

card_html <- function(card) {
  image <- if (!is.null(card$image)) norm_rel(card$image) else character()
  image <- image[nzchar(image)]
  legend <- if (!is.null(card$legend)) card$legend else character()
  legend_path <- if (length(legend) && nzchar(legend[[1]])) file.path(root, legend[[1]]) else ""
  legend_heading <- if (!is.null(card$legend_heading)) card$legend_heading else NULL
  legend_text <- extract_md_section(legend_path, legend_heading)
  image_class <- if (!is.null(card$image_class)) card$image_class else "landscape"
  gallery_label <- if (!is.null(card$gallery_label)) card$gallery_label else "Additional Draft Outputs"
  directive_ids <- card_directive_ids(card)
  summary <- or_default(card$summary, "")
  primary_image_html <- if (length(image)) {
    image_uri <- png_data_uri(file.path(root, image[[1]]))
    paste0(
      "<img class=\"review-image\" src=\"", image_uri, "\" alt=\"", html_escape(card$title), "\">\n",
      "<details class=\"full-res\"><summary>Show full-resolution image</summary>",
      "<img src=\"", image_uri, "\" alt=\"", html_escape(paste(card$title, "full resolution")), "\">",
      "</details>\n"
    )
  } else {
    "<div class=\"no-image\">No primary image for this report card.</div>\n"
  }

  paste0(
    "<section class=\"review-card ", html_escape(image_class), "\">\n",
    "<div class=\"card-heading\"><h2>", html_escape(card$title), "</h2>",
    "<span>", html_escape(card$badge), "</span></div>\n",
    "<div class=\"card-grid\">\n",
    "<div class=\"image-pane\">\n",
    primary_image_html,
    gallery_html(card$gallery, open = is_true(card$gallery_open), label = gallery_label),
    "</div>\n",
    "<aside class=\"text-pane\">\n",
    "<div class=\"summary-note\">", html_escape(summary[[1]]), "</div>\n",
    "<h3>Feedback Status</h3>\n", directive_subset_html(directive_ids, directive_table), "\n",
    "<h3>Legend And Panel Guide</h3>\n", paragraphs_html(legend_text), "\n",
    "<h3>Feedback Or Prompt Items Addressed</h3>\n", list_html(card$feedback), "\n",
    "<h3>Caveats</h3>\n", list_html(card$caveats), "\n",
    "<h3>Review Files</h3>\n", list_html(card$files), "\n",
    "</aside>\n",
    "</div>\n",
    "</section>\n"
  )
}

coverage_html <- function(manifest) {
  included <- sum(manifest$include_status == "included")
  internal <- sum(manifest$include_status == "not_reviewable_internal")
  total <- nrow(manifest)
  paste0(
    "<div class=\"coverage\"><strong>Report coverage:</strong> ",
    included, " of ", total, " drafted PNGs included",
    if (internal > 0) paste0("; ", internal, " internal PNGs excluded") else "",
    ". See the embedded appendix for <code>report_manifest.csv</code>.</div>"
  )
}

appendix_html <- function(paths) {
  blocks <- vapply(paths, function(rel_path) {
    abs_path <- file.path(root, rel_path)
    paste0(
      "<section><h3>", html_escape(rel_path), "</h3>",
      "<pre>", html_escape(read_text(abs_path)), "</pre></section>"
    )
  }, character(1))
  paste0(
    "<details class=\"appendix-block\"><summary><strong>Appendix: full notes and feedback handoff context</strong></summary>",
    paste(blocks, collapse = "\n"),
    "</details>"
  )
}

validate_cards <- function(cards, draft_pngs, directives) {
  if (!length(cards)) stop("cards must contain at least one report card.", call. = FALSE)
  sections <- vapply(cards, card_section, character(1))
  has_recommended_outputs <- any(grepl("^final_figures/recommended/", draft_pngs))
  if (has_recommended_outputs && !any(sections == "Recommended Final Drafts")) {
    stop(
      "Recommended final PNGs exist, but no card is in section 'Recommended Final Drafts'. ",
      "Show recommended drafts first.",
      call. = FALSE
    )
  }
  titles <- vapply(cards, function(card) or_default(card$title, "")[[1]], character(1))
  if (any(!nzchar(titles))) {
    stop("Every report card needs a non-empty title.", call. = FALSE)
  }
  unknown_ids <- setdiff(unique(unlist(lapply(cards, card_directive_ids), use.names = FALSE)), directives$directive_id)
  unknown_ids <- unknown_ids[nzchar(unknown_ids)]
  if (length(unknown_ids)) {
    stop("cards reference unknown directive_ids: ", paste(unknown_ids, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

cards_for_section <- function(cards, section_name) {
  cards[vapply(cards, card_section, character(1)) == section_name]
}

sectioned_cards_html <- function(cards) {
  sections <- vapply(cards, card_section, character(1))
  ordered_sections <- c(section_order, setdiff(unique(sections), section_order))
  ordered_sections <- ordered_sections[ordered_sections %in% sections]
  paste(vapply(ordered_sections, function(section_name) {
    section_cards <- cards_for_section(cards, section_name)
    paste0(
      "<section class=\"card-section\"><h2>", html_escape(section_name), "</h2>",
      paste(vapply(section_cards, card_html, character(1)), collapse = "\n"),
      "</section>"
    )
  }, character(1)), collapse = "\n")
}

# RENDER -----------------------------------------------------------------------

draft_pngs <- all_draft_pngs(png_scan_dirs)
directive_table <- validate_directive_table(directive_table)
prior_code_fidelity <- validate_prior_code_fidelity(prior_code_fidelity_path)
validate_cards(cards, draft_pngs, directive_table)
display_rows <- display_manifest(cards)
exclude_from_report <- norm_rel(exclude_from_report)
coverage_manifest <- write_coverage_manifest(draft_pngs, display_rows, exclude_from_report)
missing_pngs <- coverage_manifest$relative_path[coverage_manifest$include_status == "missing"]
if (length(missing_pngs)) {
  stop(
    "Draft PNGs are missing from review_report.html coverage. Add them to cards, gallery, ",
    "or represented_images, or mark true internals in exclude_from_report: ",
    paste(missing_pngs, collapse = ", ")
  )
}

html <- paste0(
  "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
  "<title>", html_escape(report_title), "</title>\n",
  "<style>\n",
  ":root{--bg:#f4f6f8;--panel:#fff;--ink:#1f2933;--muted:#57606a;--line:#d0d7de;--soft:#f6f8fa;--accent:#264f78;--ok:#137333;--warn:#9a6700;--bad:#b42318;--drop:#6e7781;}\n",
  "*{box-sizing:border-box;} body{font-family:Arial,Helvetica,sans-serif;margin:0;color:var(--ink);line-height:1.42;background:var(--bg);} .shell{max-width:1680px;margin:0 auto;padding:24px;}\n",
  "h1{font-size:25px;margin:0 0 8px;} h2{font-size:19px;margin:0 0 12px;} h3{font-size:13px;text-transform:uppercase;letter-spacing:.04em;color:var(--muted);margin:18px 0 7px;} p{margin:0 0 10px;} code{background:#eef2f6;padding:1px 4px;border-radius:3px;}\n",
  ".summary{background:var(--panel);border:1px solid var(--line);border-left:5px solid var(--accent);padding:14px 16px;margin-bottom:18px;max-width:1120px;}\n",
  ".directive-block{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:18px 0;} .table-wrap{overflow-x:auto;} table{border-collapse:collapse;width:100%;font-size:13px;} th,td{border:1px solid var(--line);padding:8px 9px;text-align:left;vertical-align:top;} th{background:#f6f8fa;color:#44515f;} .status-row{display:flex;flex-wrap:wrap;gap:8px;margin:4px 0 12px;}\n",
  ".status-pill{display:inline-block;border-radius:999px;padding:2px 8px;font-size:12px;font-weight:bold;white-space:nowrap;background:#eef2f6;color:#44515f;} .status-addressed{background:#dafbe1;color:var(--ok);} .status-partially-addressed{background:#fff8c5;color:var(--warn);} .status-blocked{background:#ffebe9;color:var(--bad);} .status-dropped,.status-deferred{background:#eaeef2;color:var(--drop);} .directive-mini{padding-left:18px;margin:6px 0 10px;} .directive-mini li{margin-bottom:6px;}\n",
  ".card-section{margin:24px 0 32px;} .card-section>h2{border-bottom:1px solid var(--line);padding-bottom:8px;}\n",
  ".review-card{background:var(--panel);border:1px solid var(--line);border-radius:8px;margin:18px 0 26px;box-shadow:0 1px 2px rgba(15,23,42,.06);overflow:hidden;}\n",
  ".card-heading{display:flex;align-items:center;justify-content:space-between;gap:16px;border-bottom:1px solid var(--line);padding:13px 16px;background:#fbfcfd;} .card-heading span{font-size:12px;color:#fff;background:var(--accent);border-radius:999px;padding:4px 9px;white-space:nowrap;}\n",
  ".card-grid{display:grid;grid-template-columns:minmax(620px,1fr) minmax(360px,440px);gap:0;} .image-pane{padding:14px;background:#fff;border-right:1px solid var(--line);} .text-pane{padding:14px 16px;background:#fcfcfd;}\n",
  ".review-image{display:block;width:100%;height:auto;object-fit:contain;border:1px solid var(--line);background:white;margin:0 auto;} .landscape .review-image{max-height:76vh;} .portrait .review-image{max-height:84vh;width:auto;max-width:100%;} .compact .review-image{max-height:74vh;} .subpanel .review-image{max-height:66vh;} .no-image{border:1px dashed var(--line);background:#f6f8fa;padding:24px;text-align:center;color:var(--muted);}\n",
  ".summary-note{font-weight:bold;background:#eef4fa;border-left:4px solid var(--accent);padding:9px 10px;margin-bottom:12px;} ul{padding-left:18px;margin:6px 0 10px;} li{margin-bottom:5px;} .full-res{margin-top:10px;} .full-res summary{cursor:pointer;color:var(--accent);font-weight:bold;} .full-res img{width:100%;height:auto;border:1px solid var(--line);margin-top:8px;background:white;}\n",
  ".gallery-block{margin-top:14px;border-top:1px solid var(--line);padding-top:10px;} .gallery-block summary{cursor:pointer;color:var(--accent);} .gallery-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px;margin-top:10px;} .gallery-item{margin:0;border:1px solid var(--line);background:white;padding:8px;} .gallery-item img{display:block;width:100%;height:auto;border:1px solid #eef2f6;} figcaption{font-size:11px;color:var(--muted);margin-top:6px;overflow-wrap:anywhere;} .coverage{background:#eef4fa;border:1px solid #c7d8ea;padding:10px 12px;margin:0 0 18px;max-width:1120px;}\n",
  ".appendix-block{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:18px 0;} pre{white-space:pre-wrap;background:var(--soft);border:1px solid var(--line);padding:12px;overflow-x:auto;font-size:12px;}\n",
  "@media (max-width:1100px){.shell{padding:14px;}.card-grid{grid-template-columns:1fr;}.image-pane{border-right:0;border-bottom:1px solid var(--line);}.text-pane{max-height:none;}.portrait .review-image{width:100%;max-height:none;}}\n",
  "</style>\n</head>\n<body><div class=\"shell\">\n",
  "<h1>", html_escape(report_title), "</h1>\n",
  "<div class=\"summary\"><p><strong>Reviewer action:</strong> ", html_escape(reviewer_action), "</p>",
  "<p>Each figure is paired with its legend, feedback or prompt responses, caveats, and local source files. Recommended final drafts appear before alternatives; raw refined and initial galleries should be collapsed unless they are the primary review item.</p></div>\n",
  directive_table_html(directive_table),
  prior_code_fidelity_html(prior_code_fidelity),
  coverage_html(coverage_manifest),
  sectioned_cards_html(cards),
  appendix_html(appendix_files),
  "</div></body>\n</html>\n"
)

writeLines(html, out_path)
cat("Wrote ", out_path, "\n", sep = "")
