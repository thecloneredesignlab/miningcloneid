# Up to nine instantaneous operator evaluations; no trajectory/optimizer API.
args<-commandArgs(TRUE);stopifnot(length(args)==3L)
project<-normalizePath(args[1]);study<-normalizePath(args[2]);phase<-args[3]
out<-file.path(study,'verification',phase);dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(getRversion()=='4.4.2',!file.exists(file.path(out,'native_gate.tsv')))
emit<-function(x,name){p<-file.path(out,name);tmp<-tempfile('.pending_',out);old<-options(digits=17);on.exit(options(old));write.table(x,tmp,sep='\t',quote=TRUE,row.names=FALSE);stopifnot(file.rename(tmp,p))}
env<-new.env(parent=baseenv())
extract<-function(file,names){
  found<-character()
  for(e in parse(file))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.name(e[[2]])&&as.character(e[[2]])%in%names){
    stopifnot(is.call(e[[3]]),identical(e[[3]][[1]],as.name('function')));eval(e,envir=env);found<-c(found,as.character(e[[2]]))}
  stopifnot(setequal(found,names),!anyDuplicated(found))
}
src<-file.path(study,'inputs/production')
extract(file.path(src,'o2_supply_demand_map_shared.R'),'o2sd_first_non_null')
extract(file.path(src,'o2_supply_demand_map_common_semantics.R'),c('canonical_ploidy_o2_death_mode','assert_canonical_ploidy_o2_death_mode'))
extract(file.path(src,'model_O2_supply_demand_MAP.R'),c('.first_non_null','.clip01','.assert_o2_pct','.mu_eff_of_O2','.lambda_eff_of_O2','.pmisseg_of_O2'))
fit<-readRDS(file.path(project,'oxygen/results/structural_identifiability/20261109_seed144_wp0/source/seed144/fit_result.rds'))
stopifnot(fit$cfg$ploidy_O2_death=='ploidy_related',fit$best_params$ploidy_O2_death=='ploidy_related')
cache<-file.path(out,'tmp/cpp_cache');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(src,'model_O2_supply_demand_MAP.cpp'),env=env,cacheDir=cache,rebuild=FALSE,showOutput=FALSE,verbose=FALSE)
env$cpp_o2simps_simulate_one<-function(...)stop('Production trajectory calls are not released')
env$cpp_o2simps_objective_components_map<-function(...)stop('Objective calls are not released')
params<-read.delim(file.path(study,'verification/B2/native_witness_values.tsv'),stringsAsFactors=FALSE)
calls<-data.frame(ordinal=integer(),witness=character(),oxygen=double(),label=character(),status=character())
emit(calls,'kernel_calls.tsv')
for(wid in c('W1','W2','W3'))for(e in seq_len(3L)){
  O<-c(0.2,1,20.5)[e];label<-paste0(wid,'_O',e);p<-fit$best_params
  values<-params[params$witness==wid,,drop=FALSE];stopifnot(nrow(values)==20L,!anyDuplicated(values$parameter))
  for(i in seq_len(nrow(values)))p[[values$parameter[i]]]<-as.numeric(values$value[i])
  stopifnot(nrow(calls)<9L)
  calls<-rbind(calls,data.frame(ordinal=nrow(calls)+1L,witness=wid,oxygen=O,label=label,status='started'))
  emit(calls,'kernel_calls.tsv')
  tri<-env$cpp_o2simps_build_G_for_o2_triplet(O2=O,O2_crit=p$O2_crit,N0min=22L,N0max=154L,N1min=22L,N1max=154L,
    lam_max=p$lam_max,p_mis_base=p$p_mis_base,p_misseg=p$p_misseg,k_o_mis=p$k_o_mis,p_wgd=p$p_wgd,boundary='drop',eps_tail=1e-8,
    buffer_smax=p$buffer_smax,buffer_beta=p$buffer_beta,buffer_n_exp=p$buffer_n_exp,N_unit=22L,beta_size=p$beta_size,
    O2_growth=TRUE,alpha_o2=p$alpha_o2,gamma_growth=p$gamma_growth,mu_hp=p$mu_hp,gamma_mu=p$gamma_mu,n_O=p$n_O,ploidy_O2_death=p$ploidy_O2_death)
  stopifnot(tri$nrow==133L,tri$ncol==133L)
  G<-matrix(0,133,133)
  for(k in seq_along(tri$x))G[tri$i[k],tri$j[k]]<-G[tri$i[k],tri$j[k]]+tri$x[k]
  index<-expand.grid(i=seq_len(133),j=seq_len(133))
  emit(data.frame(i=index$i,j=index$j,value=as.numeric(G)),paste0(label,'_G.tsv'))
  N<-22:154
  emit(data.frame(N=N,lambda=env$.lambda_eff_of_O2(O,p,N=N),mu=env$.mu_eff_of_O2(O,p,N=N),
      probability=env$.pmisseg_of_O2(O,p,N=N),buffer=as.numeric(tri$dead_buffer_rate),
      nonviable=as.numeric(tri$misseg_nonviable_rate),boundary=as.numeric(tri$boundary_dropped_rate)),paste0(label,'_rates.tsv'))
  calls$status[nrow(calls)]<-'returned';emit(calls,'kernel_calls.tsv')
  cat('NATIVE_INSTANTANEOUS_KERNEL_RETURNED ',label,'\n',sep='');flush.console()
}
loaded<-sort(loadedNamespaces());stopifnot(!length(intersect(loaded,c('rJava','V8','jsonlite','rjson','RJSONIO'))))
emit(data.frame(namespace=loaded),'loaded_namespaces.tsv')
emit(data.frame(status='PASS',native_calls=nrow(calls),forward_calls=0,death_mode='ploidy_related'),'native_gate.tsv')
cat('NINE_PRODUCTION_KERNELS_PASSED\n')
