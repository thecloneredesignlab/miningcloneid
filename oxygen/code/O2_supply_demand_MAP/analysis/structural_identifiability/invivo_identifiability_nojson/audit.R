# Source/configuration/data inspection only: no simulator or optimizer loading.
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==1L)
study <- normalizePath(args[1],mustWork=TRUE)
source(file.path(dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))),"adapter.R"))
e <- load_backend(study,FALSE)
cfg <- readRDS(file.path(study,"source/seed25/fit_config.rds"))
paper <- readRDS(file.path(study,"source/paper/source_config.rds"))
runtime <- readRDS(file.path(study,"source/paper/runtime_config.rds"))
stopifnot(identical(cfg,paper),identical(setdiff(names(cfg),"necrosis_mapping_csv"),setdiff(names(runtime),"necrosis_mapping_csv")),
 identical(cfg[setdiff(names(cfg),"necrosis_mapping_csv")],runtime[setdiff(names(runtime),"necrosis_mapping_csv")]))
stopifnot(unname(tools::md5sum(file.path(study,"source/numerical/model/model_O2_supply_demand_MAP.R")))=="231419cf66c4b5a2e058d08b5a24bb21",
 unname(tools::md5sum(file.path(study,"source/numerical/model/model_O2_supply_demand_MAP.cpp")))=="6e037440a184e5bc19cd6df6d3296bf7")
cfg$necrosis_mapping_csv <- file.path(study,"source/inputs/histology_to_dt_Gem_VT_20260209_v5_mapping.csv")
sc <- e$prepare_data(file.path(study,"source/inputs/dt_Gem_VT_20260209_v5.xlsx"),file.path(study,"source/inputs/all_ploidy.csv"),cfg)
rebuilt <- e$prepare_cpp_scenarios(sc,cfg)
comparison <- lapply(names(cfg$scenario_cpp),function(k) data.frame(field=k,pass=isTRUE(all.equal(cfg$scenario_cpp[[k]],rebuilt[[k]],tolerance=1e-12,check.attributes=TRUE))))
write_tsv(do.call(rbind,comparison),file.path(study,"tables/scenario_array_audit.tsv"))
if(!all(vapply(comparison,function(z) z$pass,logical(1)))) stop("Rebuilt observations differ from archived scoring arrays")
stopifnot(identical(vapply(sc,`[[`,character(1),"harvest"),cfg$harvest_param_ids))
fit <- readRDS(file.path(study,"source/seed25/deoptim_result.rds"))
p <- fit$optim$bestmem
stopifnot(length(p)==20L,identical(names(p),e$get_param_names(FALSE,FALSE,FALSE)),is.finite(fit$optim$bestval))
tab <- read.csv(file.path(study,"source/seed25/parameter_table.csv"),stringsAsFactors=FALSE)
ind <- match(names(p),tab$param_name); stopifnot(!anyNA(ind),all(tab$estimate[ind]))
lower <- setNames(tab$lower_bound[ind],names(p)); upper <- setNames(tab$upper_bound[ind],names(p))
stopifnot(all(p>=lower-1e-12),all(p<=upper+1e-12),all(lower<upper))
rp <- e$decode_params(p,FALSE,FALSE,cfg)
rt <- e$encode_params(rp,FALSE,FALSE,cfg)
stopifnot(max(abs(rt[names(p)]-p))<1e-12)
rows <- list(); j <- 0L
for(i in seq_along(sc)) {
 s <- sc[[i]]
 for(k in seq_along(s$obs_days)) { j<-j+1L; rows[[j]]<-data.frame(observation_id=paste(s$harvest,"burden",k,sep=":"),harvest=s$harvest,cohort=s$cohort,modality="burden",within_sample_index=k,day=s$obs_days[k],value=s$obs_burden[k],scored=rebuilt$keep_burden[[i]][k] && is.finite(s$obs_burden[k]) && s$obs_burden[k]>0,time_basis="recorded Day column; native score excludes nonpositive observations") }
 for(k in seq_along(s$endpoint_obs_z)) { j<-j+1L; rows[[j]]<-data.frame(observation_id=paste(s$harvest,"chromosome",k,sep=":"),harvest=s$harvest,cohort=s$cohort,modality="chromosome",within_sample_index=k,day=s$harvest_day,value=s$endpoint_obs_z[k],scored=TRUE,time_basis="production max finite burden day; assay date audit pending") }
 j<-j+1L; rows[[j]]<-data.frame(observation_id=paste(s$harvest,"necrosis",sep=":"),harvest=s$harvest,cohort=s$cohort,modality="necrosis",within_sample_index=1L,day=s$necrosis_day,value=s$obs_necrosis_fraction,scored=rebuilt$keep_necrosis[i],time_basis="production max finite burden day; assay date audit pending")
}
obs <- do.call(rbind,rows); stopifnot(!anyDuplicated(obs$observation_id))
write_tsv(obs,file.path(study,"manifests/observations.tsv"))
counts <- list(burden=sum(obs$modality=="burden" & obs$scored),ploidy=sum(lengths(rebuilt$ploidy_z)),necrosis=sum(rebuilt$keep_necrosis))
nec <- e$read_necrosis_observations(cfg$necrosis_mapping_csv,TRUE)
nec$included_in_scored_scenarios <- nec$harvest %in% vapply(sc,`[[`,character(1),"harvest")
nec$reason <- ifelse(nec$included_in_scored_scenarios,"included","not retained by production scenario selection")
write_tsv(nec,file.path(study,"tables/necrosis_mapping_audit.tsv"))
pm <- data.frame(parameter=names(p),role="active",lower=lower,upper=upper,endpoint=p,transform=ifelse(startsWith(names(p),"log10_"),"log10","identity"))
write_tsv(pm,file.path(study,"manifests/parameters.tsv"))
conditioning <- data.frame(item=c("initial_abundance","oxygen_demand_reference","harvest_initial_multiplier","crowding","beta_size","harvest_time"),status=c("fixed model conditioning","fixed reference convention","disabled","disabled","inactive","production last finite burden day"),value=c(cfg$init_total_size,cfg$o2_Nref,cfg$harvest_init_multiplier,cfg$Crowding,0,"not an independently verified assay date"))
write_tsv(conditioning,file.path(study,"manifests/conditioning.tsv"))
atomic_rds(list(cfg=cfg,scenarios=sc,endpoint=p,lower=lower,upper=upper,expected_score=fit$optim$bestval,expected_counts=counts),file.path(study,"manifests/contract.rds"))
write_kv(c(list(pass=TRUE,model_evaluations=0L,scenarios=length(sc),necrosis_raw_usable=nrow(nec),necrosis_scored=sum(rebuilt$keep_necrosis)),counts),file.path(study,"tables/preflight.tsv"))
write_tsv(data.frame(column=names(readxl::read_excel(file.path(study,"source/inputs/dt_Gem_VT_20260209_v5.xlsx")))),file.path(study,"tables/burden_source_columns.tsv"))
assert_namespaces(); print(nec); print(counts)
