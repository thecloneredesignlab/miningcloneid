#!/usr/bin/env Rscript

# Reclassify an existing validated joint soft-coupling master table under a
# different ClassA/B/C interval. Raw ratios and fitting diagnostics are kept
# unchanged; only classification fields and downstream summaries are rebuilt.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

reclassify_joint_soft_coupling_master_table <- function(
    source_analysis_dir,
    out_dir,
    class_threshold = 1.1,
    class_lower_bound = NA_real_,
    class_upper_bound = NA_real_,
    class_boundary_rule = "classb_inclusive",
    analysis_label = basename(dirname(dirname(out_dir)))) {
  source_analysis_dir <- normalizePath(source_analysis_dir, mustWork = TRUE)
  if (o2jca_path_is_within(out_dir, source_analysis_dir)) {
    stop("Reclassification output must be outside source_analysis_dir", call. = FALSE)
  }
  spec <- o2jca_classification_spec(
    class_threshold, class_lower_bound, class_upper_bound, class_boundary_rule
  )
  source_master_path <- file.path(source_analysis_dir, "soft_coupling_master_long.tsv")
  master <- o2jca_read_tsv(source_master_path)
  o2jca_assert_columns(
    master,
    c("pair_id", "seed", "parameter", "ratio_vivo_to_vitro", "log2_ratio_vivo_to_vitro"),
    "source joint soft-coupling master table"
  )
  o2jca_assert_unique_key(master, c("pair_id", "seed", "parameter"), "source joint soft-coupling master table")

  master$ratio_class <- as.character(o2jca_ratio_class(
    master$ratio_vivo_to_vitro,
    threshold = spec$threshold,
    lower_bound = spec$lower_bound,
    upper_bound = spec$upper_bound,
    boundary_rule = spec$boundary_rule
  ))
  master$class_threshold <- spec$threshold
  master$class_lower_bound <- spec$lower_bound
  master$class_upper_bound <- spec$upper_bound
  master$class_boundary_rule <- spec$boundary_rule
  master$class_scheme <- spec$scheme
  master$class_rule <- spec$class_rule

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  carry_names <- c("input_quality_summary.tsv", "pair_metadata.tsv", "parameter_process_map.tsv")
  carry <- lapply(carry_names, function(name) {
    source_path <- file.path(source_analysis_dir, name)
    if (!file.exists(source_path)) stop("Missing source analysis table: ", source_path, call. = FALSE)
    o2jca_write_tsv(o2jca_read_tsv(source_path), file.path(out_dir, name))
  })
  source_config <- o2jca_read_tsv(file.path(source_analysis_dir, "analysis_config.tsv"))
  config_value <- function(key, fallback = NA_character_) {
    hit <- source_config$value[source_config$key == key]
    if (length(hit)) as.character(hit[[1L]]) else fallback
  }
  config <- data.frame(
    key = c(
      "result_root", "analysis_label", "source_analysis_dir", "source_master_md5",
      "reclassification_only", "class_threshold", "class_lower_bound",
      "class_upper_bound", "class_boundary_rule", "class_scheme", "class_rule",
      "max_pairs", "max_seeds", "n_pairs", "n_rows"
    ),
    value = c(
      config_value("result_root"), analysis_label, source_analysis_dir,
      unname(tools::md5sum(source_master_path)), "TRUE", spec$threshold,
      spec$lower_bound, spec$upper_bound, spec$boundary_rule, spec$scheme,
      spec$class_rule, config_value("max_pairs"), config_value("max_seeds"),
      length(unique(master$pair_id)), nrow(master)
    ),
    stringsAsFactors = FALSE
  )
  files <- c(
    o2jca_write_tsv(master, file.path(out_dir, "soft_coupling_master_long.tsv")),
    o2jca_write_tsv(master[master$ratio_class == "Invalid", , drop = FALSE], file.path(out_dir, "excluded_ratio_rows.tsv")),
    unlist(carry, use.names = FALSE),
    o2jca_write_tsv(config, file.path(out_dir, "analysis_config.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "master_table_manifest.tsv"))
  invisible(master)
}

main <- function() {
  args <- o2jca_parse_args()
  source_analysis_dir <- args$source_analysis_dir %||% stop("--source_analysis_dir is required", call. = FALSE)
  out_dir <- args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  reclassify_joint_soft_coupling_master_table(
    source_analysis_dir = source_analysis_dir,
    out_dir = out_dir,
    class_threshold = o2jca_as_num(args$class_threshold, 1.1),
    class_lower_bound = o2jca_as_num(args$class_lower_bound, NA_real_),
    class_upper_bound = o2jca_as_num(args$class_upper_bound, NA_real_),
    class_boundary_rule = args$class_boundary_rule %||% "classb_inclusive",
    analysis_label = args$analysis_label %||% basename(dirname(dirname(out_dir)))
  )
}

if (identical(environment(), globalenv())) main()
