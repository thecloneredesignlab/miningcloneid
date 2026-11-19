source(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])),"common.R"))
preflight <- length(args)>2L && args[3]=="local_preflight"
out <- file.path(study,if(preflight)"preflight/local_metadata" else "results/metadata")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!preflight) {
  for(n in c("Rcpp","digest","Matrix","dplyr","DEoptim")) library(n,character.only=TRUE)
  record_environment(out)
  writeLines(system2(file.path(R.home("bin"),"R"),c("CMD","config","CXX17"),stdout=TRUE),file.path(out,"compiler.txt"))
}
load_defs()
saved <- readRDS(file.path(wp0,"source/seed144/fit_result.rds"))
pt <- file.path(wp0,"source/seed144/parameter_table_input.csv")
map <- o2prof_map(pt,saved$best_params,"zero_extended")
role <- tab(file.path(contract,"parameter_roles.tsv"))
stopifnot(identical(map$table$parameter,role$parameter),all(map$table$lower==role$lower),all(map$table$upper==role$upper),identical(map$table$transform,role$transform),all(map$table$scale==role$scale))
checks <- list()
initial <- c("init_mean_2N","init_sd_2N","init_mean_4N","init_sd_4N")
errors <- c("sigma_growth","sigma_kary")
fixed <- list(list(),saved$best_params[c(initial,errors)],saved$best_params[errors],saved$best_params[initial],list(p_misseg=0),list(buffer_beta=0),list(p_misseg=0,buffer_beta=0))
old <- tab(file.path(contract,"parameter_map_checks.tsv"))
for(i in seq_along(fixed)) {
  m <- o2prof_map(pt,saved$best_params,"zero_extended",fixed[[i]])
  p <- o2prof_project_start(m,saved$best_params)$params
  a <- o2prof_to_free(m,p);b <- o2prof_from_free(m,a);aa <- o2prof_to_free(m,b)
  gap <- max(abs(a-aa));stopifnot(length(a)==old$expected[i],gap<=1e-12)
  for(n in names(fixed[[i]]))stopifnot(b[[n]]==fixed[[i]][[n]])
  stopifnot(identical(paste(m$table$parameter[m$table$status=="inactive"],collapse=";"),if(is.na(old$inactive[i]))"" else old$inactive[i]))
  checks[[length(checks)+1L]] <- data.frame(case=old$case[i],parameter="all",free_coordinates=length(a),maximum_error=gap,pass=TRUE)
}
base <- o2prof_to_free(map,saved$best_params);spec <- o2prof_spec(map,saved$best_params)
for(i in seq_along(base))for(kind in c("single_coordinate","lower","upper")) {
  a <- base;unit <- (base[i]-spec$lower[i])/(spec$upper[i]-spec$lower[i])
  a[i] <- if(kind=="lower")spec$lower[i] else if(kind=="upper")spec$upper[i] else spec$lower[i]+(if(unit>.99)unit-.01 else unit+.01)*(spec$upper[i]-spec$lower[i])
  p <- o2prof_from_free(map,a);b <- o2prof_to_free(map,p)
  gap <- max(abs((a-b)/(spec$upper-spec$lower)))
  stopifnot(gap<=64*.Machine$double.eps)
  if(kind=="single_coordinate")stopifnot(all(b[-i]==base[-i]))
  checks[[length(checks)+1L]] <- data.frame(case=kind,parameter=map$table$parameter[i],free_coordinates=length(a),maximum_error=gap,pass=TRUE)
}
for(n in c("sigma_growth","sigma_kary")) {
  caught <- tryCatch({o2prof_map(pt,saved$best_params,"zero_extended",setNames(list(0),n));FALSE},error=function(e)TRUE)
  stopifnot(caught)
  checks[[length(checks)+1L]] <- data.frame(case="invalid_zero_error_rejected",parameter=n,free_coordinates=NA_integer_,maximum_error=NA_real_,pass=TRUE)
}
emit(do.call(rbind,checks),file.path(out,"mapping_checks.tsv"))
emit(role,file.path(out,"parameter_role_comparison.tsv"))
# Reconstruct observed summaries from original staged counts; no model or optimizer.
f <- readRDS(file.path(wp0,"source/inputs/oxygen/ploidyOxygen/data/fit_objects/fit_data.Rds"))
m <- read.csv(file.path(study,"inputs/metadata.csv"))
db <- tab(file.path(study,"inputs/cloneid_passaging_sum159_snapshot_20260731.tsv"))
g <- list()
for(id in names(f))if(is.finite(f[[id]]$g)) {
  z <- m[m$passage_id==id,];origin <- "metadata"
  if(nrow(z)<2L) {
    d <- db[db$id==id|db$passaged_from_id1==id,];d <- d[!is.na(d$id)&is.finite(as.numeric(d$correctedCount)),]
    t <- as.numeric(as.Date(substr(d$date,1,10)));t <- t-min(t)
    z <- data.frame(num_date=t,correctedCount=as.numeric(d$correctedCount));origin <- "archived_database_snapshot"
  }
  t <- z$num_date;y <- log(pmax(1,z$correctedCount))
  slope <- sum((t-mean(t))*(y-mean(y)))/sum((t-mean(t))^2)
  reference <- tab(file.path(contract,"growth_definition_audit.tsv"));j <- match(id,reference$observation_id)
  delta <- slope-reference$reconstructed_ols_g[j]
  stopifnot(is.finite(delta),abs(delta)<=1e-12,abs(slope-f[[id]]$g)<=1e-12)
  g[[length(g)+1L]] <- data.frame(observation_id=id,source=origin,raw_records=nrow(z),reconstructed_g=slope,local_difference=delta,saved_difference=slope-f[[id]]$g,pass=TRUE)
}
stopifnot(length(g)==114L)
emit(do.call(rbind,g),file.path(out,"growth_definition_comparison.tsv"))
# Real observation units and names; a density grid is never an event count.
hist <- list();kary <- list()
for(id in names(f)) {
  x <- f[[id]]
  if(length(x$kary)) kary[[length(kary)+1L]] <- data.frame(id=id,n=sum(is.finite(x$kary)))
  if(is.data.frame(x$flow)) {v<-x$flow;stopifnot(nrow(v)==500L,all(v$count>=0&v$count==floor(v$count)));hist[[length(hist)+1L]]<-data.frame(observation_id=id,bins=nrow(v),count_sum=sum(v$count),pass=TRUE)}
}
stopifnot(length(kary)==12L,sum(vapply(kary,function(x)x$n,numeric(1)))==220,length(hist)==20L)
emit(do.call(rbind,hist),file.path(out,"flow_histogram_checks.tsv"))
flow <- read.csv(file.path(wp0,"source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv"));dens <- list()
for(n in unique(flow$sample_name)) {
  v <- flow[flow$sample_name==n,];v <- v[order(v$ploidy),];value <- sum(exp(v$log_density))*median(diff(v$ploidy));stopifnot(nrow(v)==200L,abs(value-1)<=1e-10)
  dens[[length(dens)+1L]]<-data.frame(sample_name=n,grid_points=nrow(v),integral=value,pass=TRUE)
}
stopifnot(length(dens)==20L);emit(do.call(rbind,dens),file.path(out,"flow_density_checks.tsv"))
d <- tab(file.path(study,"inputs/sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"))
stopifnot(nrow(d)==90L,all(d$model_passage_id%in%names(f)),!anyDuplicated(d$model_passage_id))
stopifnot(unname(tools::md5sum(file.path(study,"inputs/sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv")))=="f346d7803596e15968b90b8213f15d95",
 all(d$dead_count>=0&d$dead_count<=d$eligible_denominator&d$eligible_denominator>0),
 all(d$dead_count==floor(d$dead_count)&d$eligible_denominator==floor(d$eligible_denominator)),
 all(abs(d$observed_dead_fraction-d$dead_count/d$eligible_denominator)<=1e-8))
original <- tab(file.path(contract,"dead_cell_sample_audit.tsv"));timing <- list()
for(i in seq_len(nrow(d))) {
  r <- original[match(d$model_passage_id[i],original$model_passage_id),]
  delta <- d$elapsed_days_since_seed[i]-r$archived_seed144_selected_day
  stopifnot(abs(delta-r$time_offset_from_archived_selected_day)<=1e-12)
  timing[[i]]<-data.frame(observation_id=d$observation_id[i],denominator=d$eligible_denominator[i],actual_day=d$elapsed_days_since_seed[i],archived_selected_day=r$archived_seed144_selected_day,time_difference=delta,partial=d$classification_partial[i],death_score_enabled=FALSE,assay_mapping_verified=FALSE)
}
emit(do.call(rbind,timing),file.path(out,"live_dead_timing_comparison.tsv"))
guard(out,preflight)
kv(list(pass=TRUE,container_runtime_verified=!preflight,map_cases=length(checks),growth_rows=length(g),growth_max_error=max(abs(do.call(rbind,g)$saved_difference)),karyotype_samples=length(kary),karyotype_cells=220L,flow_samples=length(hist),dead_rows=nrow(d),dead_max_time_difference=max(abs(do.call(rbind,timing)$time_difference)),model_propagations=0L,optimizer_calls=0L,synthetic_datasets=0L,random_draws=0L),file.path(out,"validation.tsv"))
writeLines("Metadata validation complete; no model propagation or optimization",file.path(out,"COMPLETE"))
