# Figure7 population-mean reporting policy; no simulation or source-label update.
f7pv_flags <- function(z,relative_min=1e-8) {
  source <- z$group%in%c("direct_first_generation","historical_lineage")
  # Initially empty labels, pure roots and a complete H-to-L block determine
  # exact empty compartments independently of floating subtraction residuals.
  structural <- (source&(z$q==1|z$day==0)) |
    (z$initial_ploidy==4 & (z$group=="never_HL_lineage"|z$q==1|z$day==0))
  pos <- is.finite(z$mass_fraction)&z$mass_fraction>0
  log_cells <- rep(NA_real_,nrow(z))
  count_known <- is.finite(z$log_total_cells)
  counted <- pos&count_known
  log_cells[counted] <- log(z$mass_fraction[counted])+z$log_total_cells[counted]
  log_cells[is.finite(z$mass_fraction)&z$mass_fraction==0&count_known] <- -Inf
  below_one <- ifelse(is.na(log_cells),NA,log_cells<0)
  count_status <- ifelse(is.na(log_cells),"unavailable",
    ifelse(log_cells== -Inf,"zero",ifelse(log_cells<0,"below_one","at_least_one")))
  identified <- !is.na(z$ancestry_identifiable)&z$ancestry_identifiable
  finite <- is.finite(z$mean_low_birth_rate)&is.finite(z$mean_division_rate)&is.finite(z$mean_ploidy)
  rel <- pos & z$mass_fraction>=relative_min
  # Neither the absolute count nor its availability participates in this mask.
  valid <- !structural & identified & finite & rel
  reason <- ifelse(structural,"structurally_empty",
    ifelse(!identified,"unidentified_ancestry",ifelse(!pos,"zero_or_invalid_mass",
    ifelse(!finite,"nonfinite_mean",ifelse(!rel,"relative_below_threshold","valid")))))
  data.table::data.table(valid=valid,reason=reason,log_expected_cells=log_cells,
    expected_cells_lt1=below_one,expected_cells_status=count_status,
    relative_pass=rel,finite_means=finite,structurally_empty=structural)
}

# Authenticate a separate derived layer. Original simulation tables and receipts
# retain their original hashes; this receipt never re-signs model calculations.
f7pv_render_preflight <- function(ctx) {
  root <- file.path(ctx$data_root,"Figure7")
  folder <- file.path(root,"postprocess/population_validity")
  receipt <- file.path(folder,"completion.json")
  if(!file.exists(receipt))stop("supp_fig7_8 requires cached-telemetry postprocessing: run postprocess_Figure7_population_validity.R; do not simulate")
  done <- jsonlite::fromJSON(receipt)
  if(!identical(done$status,"postprocessed_validated")||done$schema!=2L||done$relative_mass_min!=1e-8||
    !isFALSE(done$absolute_count_filter)||!isTRUE(done$expected_cells_lt1_flag_only)||
    !isFALSE(done$count_availability_filters_means)||
    !isTRUE(done$all_protocol_routes_required)||done$simulation_calls!=0||done$job_submissions!=0)
    stop("Invalid Figure7 conditional-mean reporting receipt")
  record <- ctx$current[["Figure7"]]
  if(is.null(record)||done$fingerprint!=record$receipt$fingerprint)stop("Conditional-mean layer belongs to a different numerical dataset")
  manifest <- file.path(folder,"manifest.tsv")
  if(f78ro_hash(manifest)!=done$manifest_sha256)stop("Conditional-mean layer manifest changed")
  m <- data.table::fread(manifest)
  if(anyDuplicated(paste(m$scope,m$relative_path))||any(grepl("(^|/)\\.\\.(/|$)",m$relative_path)))stop("Invalid postprocess manifest paths")
  paths <- vapply(seq_len(nrow(m)),function(i)file.path(switch(m$scope[i],
    workspace=ctx$workspace,figure7=root,overlay=folder,stop("Unknown postprocess scope")),m$relative_path[i]),"")
  f78ro_no_archive(paths)
  if(!all(file.exists(paths))||!identical(f78ro_hash(paths),m$sha256))stop("Stale or modified conditional-mean inputs/code/results; repeat cached postprocessing only")
  check_path <- file.path(folder,"checks.tsv")
  checks <- data.table::fread(check_path)
  if(!nrow(checks)||!all(checks$passed))stop("Conditional-mean postprocessing checks are not complete")
  verified <- data.frame(path=c(receipt,manifest,paths),sha256=f78ro_hash(c(receipt,manifest,paths)))
  ctx$authenticated <- unique(rbind(ctx$authenticated,verified))
  rel <- ifelse(startsWith(verified$path,paste0(ctx$data_root,"/")),
    substring(verified$path,nchar(ctx$data_root)+2L),
    paste0("workspace:",substring(verified$path,nchar(ctx$workspace)+2L)))
  extra <- data.frame(path=verified$path,relative_path=rel,sha256=verified$sha256)
  ctx$checked <- unique(rbind(ctx$checked,extra))
  ctx$population_validity_table <- file.path(folder,"low_population_time_summary.tsv")
  message("PASS: supp_fig7_8 relative mass >=1e-8, identifiable finite means, every route; expected count <1 is flag only")
  invisible(ctx$population_validity_table)
}
