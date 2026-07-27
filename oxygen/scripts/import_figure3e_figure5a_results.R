#!/usr/bin/env Rscript

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
soft_root <- "/Users/4482173/Documents/GitHub/soft_coupling"
hpc_root <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"
result_dir <- file.path(repo_root, "oxygen", "results", "figure3e_figure5a")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

invitro_root <- file.path(soft_root, "oxygen/results/fit_invitro_O2_buffering_500seed")
invitro_daily_path <- file.path(invitro_root, "seed10/invitro_daily_counts.tsv")
invitro_distribution_path <- file.path(invitro_root, "seed10/invitro_distribution_summary.tsv")

joint_run_name <- "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540"
joint_root <- file.path(soft_root, "oxygen/results", joint_run_name)
joint_analysis_root <- file.path(soft_root, "oxygen/results/analysis/joint_coupling", joint_run_name)
joint_task_runs_path <- file.path(joint_root, ".multi_warmup_task_table_runs.tsv")
joint_pair_metadata_path <- file.path(joint_analysis_root, "tables/soft_coupling/pair_metadata.tsv")
joint_args_path <- file.path(joint_root, "fit_joint_tsne_vi_seed138_C03Sc01_vt_seed10/run_effective_args.tsv")
joint_process_map_path <- file.path(joint_analysis_root, "tables/soft_coupling/parameter_process_map.tsv")

required_sources <- c(
  invitro_daily_path, invitro_distribution_path, joint_task_runs_path,
  joint_pair_metadata_path, joint_args_path, joint_process_map_path
)
missing_sources <- required_sources[!file.exists(required_sources)]
if (length(missing_sources)) {
  stop("Missing required source files: ", paste(missing_sources, collapse = ", "), call. = FALSE)
}

path_tokens <- function(x) strsplit(as.character(x), "_", fixed = TRUE)
path_length <- function(x) lengths(path_tokens(x))
path_has_low_o2 <- function(x) {
  vapply(path_tokens(x), function(z) any(as.numeric(z) < 20.5), logical(1L))
}

daily <- utils::read.delim(invitro_daily_path, check.names = FALSE, stringsAsFactors = FALSE)
distribution <- utils::read.delim(invitro_distribution_path, check.names = FALSE, stringsAsFactors = FALSE)

daily_endpoint <- daily[
  daily$cohort == "4N" & daily$day == daily$selected_day,
  , drop = FALSE
]
candidate_paths <- unique(daily_endpoint$segment_id[path_has_low_o2(daily_endpoint$segment_id)])
if (!length(candidate_paths)) stop("No oxygen-deprived 4N lineage was found.", call. = FALSE)
terminal_path <- candidate_paths[[which.max(path_length(candidate_paths))]]
terminal_tokens <- path_tokens(terminal_path)[[1L]]
terminal_prefixes <- vapply(
  seq_along(terminal_tokens),
  function(i) paste(terminal_tokens[seq_len(i)], collapse = "_"),
  character(1L)
)

daily_branch <- daily_endpoint[daily_endpoint$segment_id %in% terminal_prefixes, , drop = FALSE]
distribution_branch <- distribution[
  distribution$cohort == "4N" & distribution$segment_id %in% terminal_prefixes,
  , drop = FALSE
]
weighted <- stats::aggregate(
  cbind(weighted_chromosome = distribution_branch$N * distribution_branch$fraction,
        distribution_mass = distribution_branch$fraction),
  by = list(segment_id = distribution_branch$segment_id),
  FUN = sum
)
weighted$predicted_mean_chromosome <- weighted$weighted_chromosome / weighted$distribution_mass

figure3e <- merge(
  daily_branch,
  weighted[, c("segment_id", "predicted_mean_chromosome")],
  by = "segment_id",
  all = FALSE,
  sort = FALSE
)
figure3e$lineage_passage <- path_length(figure3e$segment_id)
figure3e$live_burden_fraction <- figure3e$burden_live / figure3e$burden_total
figure3e$direct_hypoxia_death_burden_fraction <- figure3e$burden_dead_hypoxia / figure3e$burden_total
figure3e$nonviable_daughter_burden_fraction <- figure3e$burden_dead_buffer / figure3e$burden_total
figure3e <- figure3e[order(figure3e$lineage_passage), c(
  "segment_id", "cohort", "lineage_passage", "passage_index", "oxygen_pct", "selected_day",
  "predicted_mean_chromosome", "live_burden_fraction",
  "direct_hypoxia_death_burden_fraction", "nonviable_daughter_burden_fraction"
)]
if (nrow(figure3e) != length(terminal_prefixes) || anyDuplicated(figure3e$lineage_passage) ||
    any(!is.finite(as.matrix(figure3e[, 7:10]))) ||
    tail(figure3e$predicted_mean_chromosome, 1L) >= figure3e$predicted_mean_chromosome[[1L]]) {
  stop("Figure 3E branch summary failed lineage, finiteness, or direction checks.", call. = FALSE)
}
figure3e_path <- file.path(result_dir, "figure3e_4n_deprived_passage_summary.csv")
utils::write.csv(figure3e, figure3e_path, row.names = FALSE, quote = TRUE)

pair_metadata <- utils::read.delim(joint_pair_metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
task_runs <- utils::read.delim(
  joint_task_runs_path, header = FALSE, col.names = c("warmup_label", "hpc_run_directory"),
  check.names = FALSE, stringsAsFactors = FALSE
)
args <- utils::read.delim(joint_args_path, check.names = FALSE, stringsAsFactors = FALSE)
arg_value <- function(key) {
  value <- args$value[args$key == key]
  if (length(value) != 1L) stop("Expected one effective argument for ", key, call. = FALSE)
  value[[1L]]
}

pair_metadata$warmup_label <- sub("^fit_joint_", "", pair_metadata$pair_id)
joint_pairs <- merge(pair_metadata, task_runs, by = "warmup_label", all = FALSE, sort = FALSE)
joint_pairs$optimizer_seed_count <- vapply(joint_pairs$pair_id, function(pair_id) {
  pair_dir <- file.path(joint_root, pair_id)
  seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = FALSE)
  sum(grepl("^seed[0-9]+$", seed_dirs))
}, integer(1L))
joint_pairs$joint_run_name <- joint_run_name
joint_pairs$soft_coupling_sigma_default <- as.numeric(arg_value("joint_soft_coupling_sigma_default"))
joint_pairs$soft_coupling_welsch_c <- as.numeric(arg_value("joint_soft_coupling_welsch_c"))
joint_pairs$warmup_sigmaN <- as.numeric(arg_value("joint_warmup_sigmaN"))
joint_pairs <- joint_pairs[order(joint_pairs$invivo_cluster, joint_pairs$invivo_subcluster), c(
  "joint_run_name", "pair_id", "warmup_label", "warmup_method", "invivo_seed",
  "invivo_cluster", "invivo_subcluster", "invitro_seed", "optimizer_seed_count",
  "warmup_sigmaN", "soft_coupling_sigma_default", "soft_coupling_welsch_c", "hpc_run_directory"
)]
if (nrow(joint_pairs) != 6L || length(unique(joint_pairs$invivo_cluster)) != 3L ||
    !all(table(joint_pairs$invivo_cluster) == 2L) || !identical(unique(joint_pairs$invitro_seed), 10L) ||
    !all(joint_pairs$optimizer_seed_count == 500L)) {
  stop("Figure 5A warm-up pairing does not match the expected 3 x 2, seed10-anchored, 500-seed design.", call. = FALSE)
}
joint_pairs_path <- file.path(result_dir, "figure5a_joint_warmup_pairs.csv")
utils::write.csv(joint_pairs, joint_pairs_path, row.names = FALSE, quote = TRUE)

soft_parameters <- strsplit(arg_value("joint_soft_coupling_params"), ",", fixed = TRUE)[[1L]]
process_map <- utils::read.delim(joint_process_map_path, check.names = FALSE, stringsAsFactors = FALSE)
soft_setup <- process_map[match(soft_parameters, process_map$parameter), c(
  "parameter", "primary_process", "downstream_processes", "effect_direction"
)]
soft_setup$soft_coupling_sigma_default <- as.numeric(arg_value("joint_soft_coupling_sigma_default"))
soft_setup$soft_coupling_welsch_c <- as.numeric(arg_value("joint_soft_coupling_welsch_c"))
if (nrow(soft_setup) != 14L || anyNA(soft_setup$primary_process) || anyDuplicated(soft_setup$parameter)) {
  stop("Figure 5A requires 14 uniquely mapped soft-coupled parameters.", call. = FALSE)
}
soft_setup_path <- file.path(result_dir, "figure5a_soft_coupled_parameters.csv")
utils::write.csv(soft_setup, soft_setup_path, row.names = FALSE, quote = TRUE)

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  if (!length(output)) return(NA_character_)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
to_hpc_path <- function(path) sub(paste0("^", soft_root), hpc_root, path)
generating_script <- normalizePath(file.path(script_dir, "import_figure3e_figure5a_results.R"), mustWork = TRUE)
provenance_row <- function(destination, source, purpose, source_git_commit = NA_character_) {
  data.frame(
    destination_file = basename(destination),
    destination_sha256 = sha256(destination),
    destination_bytes = as.numeric(file.info(destination)$size),
    source_host = "red.moffitt.org (local synchronized copy)",
    origin_directory = dirname(to_hpc_path(source)),
    origin_file = basename(source),
    origin_path = to_hpc_path(source),
    local_source_path = source,
    source_sha256 = sha256(source),
    transfer_type = "derived_filtered_local",
    purpose = purpose,
    source_git_commit = source_git_commit,
    generating_script = generating_script,
    recorded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

provenance <- do.call(rbind, list(
  provenance_row(figure3e_path, invitro_daily_path, "Figure 3E passage-end live/dead burden components"),
  provenance_row(figure3e_path, invitro_distribution_path, "Figure 3E passage-end mean chromosome number"),
  provenance_row(joint_pairs_path, joint_pair_metadata_path, "Figure 5A in vivo cluster/subcluster warm-start identities", "f8abc65b46d3bd2dfaf95ae284123f07313bc326"),
  provenance_row(joint_pairs_path, joint_task_runs_path, "Figure 5A exact HPC joint run directories", "3d4985997c54f34c3b1a30cc72462a282cff4817"),
  provenance_row(joint_pairs_path, joint_args_path, "Figure 5A optimizer-seed and soft-coupling configuration", "3d4985997c54f34c3b1a30cc72462a282cff4817"),
  provenance_row(soft_setup_path, joint_args_path, "Figure 5A active soft-coupled parameter list and penalty settings", "3d4985997c54f34c3b1a30cc72462a282cff4817"),
  provenance_row(soft_setup_path, joint_process_map_path, "Figure 5A audited biological-process grouping", "f8abc65b46d3bd2dfaf95ae284123f07313bc326")
))
provenance <- provenance[order(provenance$destination_file, provenance$origin_path), , drop = FALSE]
provenance_path <- file.path(result_dir, "source_file_provenance.csv")
utils::write.csv(provenance, provenance_path, row.names = FALSE, quote = TRUE, na = "")

message("Wrote ", figure3e_path, " [", nrow(figure3e), " lineage passages]")
message("Wrote ", joint_pairs_path, " [", nrow(joint_pairs), " joint warm-up pairs]")
message("Wrote ", soft_setup_path, " [", nrow(soft_setup), " soft-coupled parameters]")
message("Wrote ", provenance_path, " [", nrow(provenance), " source links]")
