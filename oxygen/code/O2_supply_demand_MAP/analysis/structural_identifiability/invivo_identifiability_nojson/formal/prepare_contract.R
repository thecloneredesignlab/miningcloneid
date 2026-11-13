# Formal contract preparation: metadata and native parameter transformations only.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"));e <- load_backend(study,FALSE)
contract <- readRDS(file.path(study,"reference/pilot_contract.rds"))
direct <- e$build_transformed_parameter_table(file.path(study,"source/seed25/parameter_table_input.csv"),FALSE,FALSE,TRUE,FALSE)
stopifnot(identical(names(direct$optimizer$lower),names(contract$endpoint)))
old <- contract
contract$lower <- direct$optimizer$lower;contract$upper <- direct$optimizer$upper
stopifnot(max(abs(contract$lower-old$lower),abs(contract$upper-old$upper))<1e-12,
 all(contract$lower<contract$upper))
contract$endpoint <- pmin(contract$upper,pmax(contract$lower,contract$endpoint))
contract$rho_grid <- seq(contract$lower[['log10_rho_2N']],contract$upper[['log10_rho_2N']],length.out=7L)
atomic_rds(contract,file.path(study,"manifests/contract.rds"))
write_tsv(data.frame(parameter=names(contract$endpoint),lower=contract$lower,upper=contract$upper,
 pilot_endpoint=old$endpoint,formal_reference=contract$endpoint,reference_projection=contract$endpoint-old$endpoint,
 lower_difference=contract$lower-old$lower,upper_difference=contract$upper-old$upper),file.path(study,"tables/formal_bound_audit.tsv"))
write_tsv(direct$transformed_output,file.path(study,"manifests/transformed_parameters.tsv"))
write_tsv(data.frame(index=1:7,log10_density=contract$rho_grid,density=10^contract$rho_grid),file.path(study,"manifests/density_grid.tsv"))
specs <- list()
add <- function(task,role,mode,policy="production",family="reference",init_mode="hybrid",seed=1L,start_seed=25L,
 source_task="",required_check="",selector="seed",density_index=0L,fixture=FALSE) {
 stopifnot(is.null(specs[[task]]))
 specs[[task]] <<- list(task=task,role=role,mode=mode,policy=policy,family=family,init_mode=init_mode,
 seed=as.integer(seed),start_seed=as.integer(start_seed),source_task=source_task,required_check=required_check,
 selector=selector,density_index=as.integer(density_index),fixed=if(density_index>0L) setNames(contract$rho_grid[density_index],"log10_rho_2N") else NULL,
 fixture=fixture,NP=if(fixture) 40L else 256L,itermax=if(fixture) 3L else 1000L,reltol=1e-6,steptol=100L,
 local_maxit=if(fixture) 2L else 15000L)
}
add("engineering_parallel","engineering_parallel","engineering_parallel",seed=28100L,fixture=TRUE)
for(policy in c("production","data_only")) {
 add(paste0("engineering_serial_",policy),"engineering_serial","engineering_serial",policy,seed=28101L,fixture=TRUE)
 add(paste0("reference_",policy),"source_reference_checks","reference",policy)
 for(k in 1:4) {
  chain <- c("h25","h392","u1","u2")[k];family <- if(k<=2) "hybrid" else "uniform"
  prefix <- paste0("b_",policy,"_",chain)
  add(paste0(prefix,"_global"),"baseline_global","global",policy,family,family,seed=29000L+100L*match(policy,c("production","data_only"))+k,
   start_seed=if(k==1) 25L else if(k==2) 392L else 0L)
  add(paste0(prefix,"_local"),"baseline_local","local",policy,family,source_task=paste0(prefix,"_global"),selector="task")
  add(paste0(prefix,"_fresh"),"baseline_fresh","fresh",policy,family,source_task=paste0(prefix,"_local"),selector="task")
  add(paste0(prefix,"_directions"),"baseline_directional","directions",policy,family,
   source_task=paste0(prefix,"_local"),required_check=paste0(prefix,"_fresh"),selector="task")
 }
 add(paste0("b_",policy,"_reserve_global"),"reserve_baseline_global","global",policy,"uniform","uniform",seed=31000L+match(policy,c("production","data_only")),start_seed=0L)
 for(kind in c("reserve_g","reserve_probe")) {
  prefix <- paste0("b_",policy,"_",kind)
  add(paste0(prefix,"_local"),"reserve_baseline_local","local",policy,"restart",
   source_task=if(kind=="reserve_g") paste0("b_",policy,"_reserve_global") else "",
   selector=if(kind=="reserve_g") "task" else "baseline_probe")
  add(paste0(prefix,"_fresh"),"reserve_baseline_fresh","fresh",policy,"restart",source_task=paste0(prefix,"_local"),selector="task")
  add(paste0(prefix,"_directions"),"reserve_baseline_directional","directions",policy,"restart",source_task=paste0(prefix,"_local"),required_check=paste0(prefix,"_fresh"),selector="task")
 }
 for(i in 1:7) {
  prefix <- paste0("d_",policy,"_",i)
  add(paste0(prefix,"_global"),"density_global","global",policy,"global","hybrid",seed=32000L+100L*match(policy,c("production","data_only"))+i,
   selector="baseline_best",density_index=i)
  add(paste0(prefix,"_local"),"density_global_local","local",policy,"global",source_task=paste0(prefix,"_global"),selector="task",density_index=i)
  for(direction in c("asc","desc")) {
   first <- if(direction=="asc") i==1L else i==7L;previous <- if(direction=="asc") i-1L else i+1L
   add(paste0(prefix,"_",direction),"density_continuation_local","local",policy,direction,
    source_task=if(first) "" else paste0("d_",policy,"_",previous,"_",direction),
    selector=if(first) "baseline_best" else "task",density_index=i)
  }
  add(paste0(prefix,"_fresh"),"density_fresh","profile_fresh",policy,"selected",selector="density_selected",density_index=i)
  add(paste0(prefix,"_reserve_global"),"reserve_density_global","global",policy,"repair","hybrid",seed=34000L+100L*match(policy,c("production","data_only"))+i,
   selector="density_selected",density_index=i)
  add(paste0(prefix,"_reserve_local"),"reserve_density_local","local",policy,"repair",selector="density_selected",density_index=i)
  add(paste0(prefix,"_reserve_fresh"),"reserve_density_fresh","profile_fresh",policy,"repair",selector="density_selected",density_index=i)
  add(paste0(prefix,"_reserve_directions"),"reserve_density_directional","directions",policy,"repair",selector="density_selected",density_index=i)
 }
}
atomic_rds(specs,file.path(study,"manifests/task_specs.rds"))
fields <- c("task","role","mode","policy","family","init_mode","seed","start_seed","source_task","required_check","selector","density_index","fixture")
write_tsv(do.call(rbind,lapply(specs,function(s) as.data.frame(s[fields],stringsAsFactors=FALSE))),file.path(study,"manifests/task_plan.tsv"))
atomic_rds(list(prior_core_seconds=30491L,formal_cap_core_seconds=648000L,cumulative_cap_core_seconds=691200L,
 score_tolerance=1e-8,prediction_absolute_tolerance=1e-8,prediction_relative_tolerance=1e-8,
 baseline_recovery=0.002,direction_improvement=0.002,profile_recovery=0.01,steps=c(1e-3,1e-4)),file.path(study,"manifests/protocol.rds"))
write_kv(list(pass=TRUE,model_evaluations=0L,active_parameters=20L,density_anchors=7L,
 formal_tasks_in_registry=length(specs),maximum_coordinate_projection=max(abs(contract$endpoint-old$endpoint)),
 maximum_bound_difference=max(abs(contract$lower-old$lower),abs(contract$upper-old$upper))),file.path(study,"tables/metadata_preflight.tsv"))
assert_namespaces();cat("Formal metadata prepared;",length(specs),"predeclared task identities, bounded by role attempt caps. Zero model calls.\n")
