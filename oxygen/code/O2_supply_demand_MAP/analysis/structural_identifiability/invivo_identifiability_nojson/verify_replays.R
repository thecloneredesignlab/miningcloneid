# Inspect saved evidence only; no model loading or evaluation.
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==1L)
study <- normalizePath(args[1]); pkg <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(pkg,"adapter.R"))
numeric_same <- function(a,b) {
  stopifnot(identical(dim(a),dim(b)),identical(is.na(a),is.na(b)))
  ok <- !is.na(a); stopifnot(all(abs(a[ok]-b[ok])<=1e-8+1e-8*abs(a[ok])))
}
z <- lapply(1:3,function(i) readRDS(file.path(study,"tests",paste0("repair1_rep",i),"replay.rds")))
telemetry <- c("cache_g_build","cache_g_hit","cache_g_hysteresis")
for(k in 2:3) {
 stopifnot(identical(component_values(z[[1]]$components),component_values(z[[k]]$components)),identical(z[[1]]$predictions$tables,z[[k]]$predictions$tables))
 for(i in seq_along(z[[1]]$predictions$native)) {
  a <- z[[1]]$predictions$native[[i]];b <- z[[k]]$predictions$native[[i]]
  stopifnot(identical(a$arguments,b$arguments))
  for(n in setdiff(names(a$result),telemetry)) numeric_same(a$result[[n]],b$result[[n]])
 }
}
write_kv(list(pass=TRUE,fresh_processes=3L,objective=z[[1]]$components$L,difference=z[[1]]$components$L-readRDS(file.path(study,"manifests/contract.rds"))$expected_score,components_identical=TRUE,predictions_identical=TRUE,native_states_and_oxygen_identical=TRUE,cache_counters_excluded=TRUE),file.path(study,"tables/replay_gate.tsv"))
ctx <- readRDS(file.path(study,"manifests/contract.rds"))
pl <- read.delim(file.path(study,"source/inputs/all_ploidy.csv"),stringsAsFactors=FALSE)
pl$raw_data_row <- seq_len(nrow(pl));pl$harvest <- sub(".sps.cbs","",pl$file,fixed=TRUE)
rows <- lapply(ctx$scenarios,function(s) {
 p <- pl[pl$harvest==s$harvest & is.finite(pl$total_chromosomes),]
 stopifnot(nrow(p)==length(s$endpoint_obs_z),isTRUE(all.equal(as.numeric(p$total_chromosomes),s$endpoint_obs_z,tolerance=1e-12)))
 data.frame(observation_id=paste(s$harvest,"chromosome",seq_len(nrow(p)),sep=":"),harvest=s$harvest,source_file=p$file,cell_id=p$cell_id,raw_data_row=p$raw_data_row,total_chromosomes=p$total_chromosomes)
})
rows <- do.call(rbind,rows);stopifnot(!anyDuplicated(rows$observation_id),!anyDuplicated(paste(rows$source_file,rows$cell_id)))
write_tsv(rows,file.path(study,"tables/chromosome_source_rows.tsv"))
cat("Saved-state/oxygen replay and physical cell identity checks passed.\n")
