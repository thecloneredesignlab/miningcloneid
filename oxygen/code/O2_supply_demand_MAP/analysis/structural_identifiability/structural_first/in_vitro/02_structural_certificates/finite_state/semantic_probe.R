# Production-SIF metadata and rate-function dispatch only. No trajectories.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
project <- normalizePath(args[1]); study <- normalizePath(args[2])
out <- file.path(study,'verification/B0');dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(out,'semantic_gate.tsv')),getRversion()=='4.4.2')
writeLines(as.character(Sys.getpid()),file.path(out,'worker_pid.txt'))
emit <- function(x,name)write.table(x,file.path(out,name),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
env <- new.env(parent=baseenv())
extract <- function(file,names) {
  code <- parse(file); found <- character()
  for(e in code) {
    if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.name(e[[2]])&&as.character(e[[2]])%in%names) {
      stopifnot(is.call(e[[3]]),identical(e[[3]][[1]],as.name('function')))
      eval(e,envir=env);found<-c(found,as.character(e[[2]]))
    }
  }
  stopifnot(setequal(found,names),!anyDuplicated(found))
  emit(data.frame(source=basename(file),function_name=found),'last_extracted_functions.tsv')
}
src<-file.path(study,'inputs/production')
extract(file.path(src,'o2_supply_demand_map_shared.R'),'o2sd_first_non_null')
extract(file.path(src,'o2_supply_demand_map_common_semantics.R'),c('canonical_ploidy_o2_death_mode','assert_canonical_ploidy_o2_death_mode'))
extract(file.path(src,'model_O2_supply_demand_MAP.R'),c('.first_non_null','.clip01','.assert_o2_pct','.mu_eff_of_O2','.lambda_eff_of_O2','.pmisseg_of_O2'))
fitpath<-file.path(project,'oxygen/results/structural_identifiability/20261109_seed144_wp0/source/seed144/fit_result.rds')
fit<-readRDS(fitpath);cfg<-fit$cfg;p<-fit$best_params
stopifnot(identical(env$canonical_ploidy_o2_death_mode(cfg$ploidy_O2_death),'ploidy_related'),
 identical(env$assert_canonical_ploidy_o2_death_mode(p$ploidy_O2_death),'ploidy_related'),
 isTRUE(cfg$O2_growth),isTRUE(cfg$Crowding),cfg$crowding=='logistic',cfg$K==1e12,
 !isTRUE(cfg$o2_burden_feedback),cfg$N_MIN==22,cfg$N_MAX==154,cfg$N_UNIT==22,cfg$DT==0.05,p$k_clear==1e-4)
checks<-list();j<-0L
for(O in c(0.2,1,20.5))for(N in c(22,44,88,154)) {
 h<-p$O2_crit^p$n_O/(p$O2_crit^p$n_O+O^p$n_O)
 actual<-env$.mu_eff_of_O2(O,p,N=N)
 expected<-p$mu_hp*h*(N/44)^p$gamma_mu
 other<-p$mu_hp*h*(1+pmax(N/44-1,0)^p$gamma_mu)
 stopifnot(abs(actual-expected)<=1e-13*max(1,abs(expected)))
 if(N!=44)stopifnot(abs(expected-other)>1e-12)
 j<-j+1L;checks[[j]]<-data.frame(oxygen=O,N=N,actual=actual,ploidy_related=expected,diploid_NULL=other,residual=actual-expected)
}
emit(do.call(rbind,checks),'death_mode_checks.tsv')
emit(data.frame(key=c('cfg_mode','run_mode','R_version','grid','K','k_clear','probability_epsilon','continuous_state','native_generator_calls','forward_calls'),
 value=c(cfg$ploidy_O2_death,p$ploidy_O2_death,as.character(getRversion()),'22:154','1e12','1e-4','1e-8','133 live plus two passive retained totals','0','0')),'semantic_settings.tsv')
loaded<-sort(loadedNamespaces())
stopifnot(!length(intersect(loaded,c('rJava','V8','jsonlite','rjson','RJSONIO'))))
emit(data.frame(namespace=loaded),'loaded_namespaces.tsv')
emit(data.frame(status='PASS',death_mode='ploidy_related',source_loading='allowlisted unchanged function expressions only',native_generator_calls=0,forward_calls=0),'semantic_gate.tsv')
cat('PRODUCTION_SIF_SEMANTIC_GATE_PASSED\n')
