# New forward-only adapter. No optimizer, synthetic observation or legacy counter.
sf_tab <- function(p) read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)
sf_emit <- function(x,p) {
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  old<-options(digits=17);on.exit(options(old));q<-tempfile('.pending_',dirname(p))
  write.table(x,q,sep='\t',quote=TRUE,row.names=FALSE,na='NA');stopifnot(file.rename(q,p))
}
sf_rds <- function(x,p) {q<-tempfile('.pending_',dirname(p));saveRDS(x,q);stopifnot(file.rename(q,p))}
sf_hash <- function(x) digest::digest(x,algo='sha256')
sf_inventory <- function(out) {
  ns<-sort(loadedNamespaces());stopifnot(!length(intersect(ns,c('rJava','V8','jsonlite','rjson','RJSONIO','DEoptim','nloptr'))))
  sf_emit(data.frame(namespace=ns,version=vapply(ns,function(n)as.character(getNamespaceVersion(n)),character(1))),file.path(out,'loaded_namespaces.tsv'))
}
sf_backend <- function(study,out,native=FALSE) {
  stopifnot(getRversion()=='4.4.2',startsWith(normalizePath(tempdir()),paste0(normalizePath(file.path(out,'tmp')),'/')))
  suppressPackageStartupMessages(library(digest));suppressPackageStartupMessages(library(dplyr))
  env<-new.env(parent=.GlobalEnv);records<-list();src<-file.path(study,'inputs/production')
  files<-c('o2_supply_demand_map_shared.R','o2_supply_demand_map_common_semantics.R','model_O2_supply_demand_MAP.R',
    'o2_supply_demand_map_invitro_io_utils.R','o2_supply_demand_map_invitro_lineage_utils.R','o2_supply_demand_map_invitro_summary_utils.R',
    'o2_supply_demand_map_invitro_lineage_simulation_utils.R','o2_supply_demand_map_invitro_objective_utils.R','parameter_helpers.R')
  constants<-c('.chr_lengths_default_bp_1to22','.chr_lengths_default_ord_desc','.chr_lengths_default_denom','INVITRO_PASSAGE_IMPLEMENTATION','o2prof_symbols')
  for(file in files)for(e in parse(file.path(src,file)))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])) {
    name<-as.character(e[[2]]);rhs<-e[[3]]
    func<-is.call(rhs)&&identical(rhs[[1]],as.name('function'))
    closure<-name%in%c('ivt_load_fit_objects','ivt_run_segment_fixed_o2')
    if(func||closure||name%in%constants) {
      eval(e,envir=env);records[[length(records)+1L]]<-data.frame(file=file,name=name,kind=if(func)'unchanged_function' else if(closure)'unchanged_pure_closure_definition' else 'unchanged_constant')
    }
  }
  env$.first_non_null_local<-env$o2sd_first_non_null
  stopifnot(identical(env$INVITRO_PASSAGE_IMPLEMENTATION,'v2'),!exists('.init_cpp_o2simps_backend',env,inherits=FALSE))
  sf_emit(do.call(rbind,records),file.path(out,'loaded_source_definitions.tsv'))
  if(native) {
    cache<-file.path(out,'tmp/native_cache');dir.create(cache,recursive=TRUE)
    Rcpp::sourceCpp(file.path(src,'model_O2_supply_demand_MAP.cpp'),env=env,cacheDir=cache,rebuild=FALSE,showOutput=FALSE,verbose=FALSE)
    stopifnot(exists('cpp_o2simps_simulate_one',env,inherits=FALSE));env$.USE_CPP_O2SIMPS_BACKEND<-TRUE
    sf_emit(data.frame(export=c('cpp_o2simps_simulate_one','cpp_o2simps_build_G_for_o2_triplet'),loading='direct sourceCpp into private cache; original timed bootstrap not evaluated'),file.path(out,'native_interface.tsv'))
  }
  env
}
sf_context <- function(study,env) {
  fit<-readRDS(file.path(study,'inputs/data/fit_result.rds'));cfg<-fit$cfg
  stopifnot(cfg$ploidy_O2_death=='ploidy_related',fit$best_params$ploidy_O2_death=='ploidy_related')
  map<-env$o2prof_map(file.path(study,'inputs/data/parameter_table_input.csv'),fit$best_params,'zero_extended')
  projected<-env$o2prof_project_start(map,fit$best_params);p<-projected$params
  oldp<-readRDS(file.path(study,'inputs/evidence/d0_truth_parameters.rds'))
  stopifnot(isTRUE(all.equal(p,oldp,tolerance=0,check.attributes=TRUE)))
  objects<-env$ivt_load_fit_objects(repo_root=file.path(study,'inputs/data'),fit_objects_dir=file.path(study,'inputs/data'),flow_csv_path=file.path(study,'inputs/data/g0g1_ploidy_density_grid.csv'))
  stopifnot(!isTRUE(objects$death_enabled))
  adapters<-lapply(c('2N','4N'),function(c)env$ivt_build_job_table_adapter(objects[[paste0('jobs_',c)]],objects$fit_data,c));names(adapters)<-c('2N','4N')
  list(p=p,cfg=cfg,map=map,projection=projected$projection,core=env$build_model_core(cfg=cfg),objects=objects,adapters=adapters,flow_sd=fit$best_components$sigma_flow_ploidy)
}
sf_unit <- function(env,b,p) {
  spec<-env$o2prof_spec(b$map,b$p)
  setNames((unname(env$o2prof_to_free(b$map,p))-spec$lower)/(spec$upper-spec$lower),b$map$table$parameter)
}
sf_params <- function(env,b,u) {
  spec<-env$o2prof_spec(b$map,b$p)
  env$o2prof_from_free(b$map,setNames(spec$lower+unname(u)*(spec$upper-spec$lower),spec$param_name))
}
sf_predictions <- function(env,b,runs,p) {
  g<-numeric();k<-list();f<-list();selection<-list();states<-list()
  for(cohort in names(runs))for(seg in runs[[cohort]]$segment_results) {
    frac<-as.numeric(seg$selection$selected_frac)
    stopifnot(all(is.finite(frac)),min(frac)>=0,abs(sum(frac)-1)<1e-12)
    rate<-env$ivt_log_growth_rate(seg$sim$Ntot_live_obs[1],seg$selection$selected_live_cells,seg$selection$selected_day)
    kp<-env$ivt_smooth_kary_distribution(runs[[cohort]]$grid_pre,frac,p$sigma_kary)
    for(id in seg$segment$data_ids) {
      o<-env$ivt_observed_passage_summary(b$objects$fit_data[[id]])
      if(is.finite(o$observed_growth))g[id]<-rate
      if(length(o$observed_kary))k[[id]]<-kp
      if(!is.null(o$observed_flow)&&nrow(o$observed_flow))f[[id]]<-env$ivt_predicted_flow_density(o$observed_flow$ploidy,runs[[cohort]]$grid_pre,frac,b$flow_sd,b$cfg$N_UNIT)
    }
    parent<-seg$parent_boundary_selection
    selection[[length(selection)+1L]]<-data.frame(cohort=cohort,segment_id=seg$segment$segment_id,index=seg$selection$selected_index,day=seg$selection$selected_day,
      parent_index=if(is.null(parent))NA_integer_ else parent$selected_index,parent_day=if(is.null(parent))NA_real_ else parent$selected_day)
    states[[paste(cohort,seg$segment$segment_id)]]<-seg$sim$live_state_obs
  }
  stopifnot(length(g)==114,length(k)==12,length(f)==20)
  vector<-g/b$p$sigma_growth
  names(vector)<-paste0('growth:',names(g))
  for(id in names(k)) {v<-k[[id]];stopifnot(abs(sum(v)-1)<1e-12);vector<-c(vector,setNames(v[-length(v)],paste0('kary:',id,':',seq_len(length(v)-1))))}
  for(id in names(f)) {o<-env$ivt_observed_passage_summary(b$objects$fit_data[[id]])$observed_flow;step<-median(diff(sort(unique(o$ploidy))));v<-f[[id]]*step;stopifnot(abs(sum(v)-1)<1e-12);vector<-c(vector,setNames(v[-length(v)],paste0('flow:',id,':',seq_len(length(v)-1))))}
  list(growth=g,kary=k,flow=f,selection=do.call(rbind,selection),states=states,vector=vector)
}
sf_score <- function(env,b,runs,p,pred) {
  summary<-dplyr::bind_rows(lapply(runs,function(r)env$ivt_collect_lineage_summary(r,b$objects$fit_data)))
  a<-env$ivt_growth_loglik_df(summary,p$sigma_growth)$loglik
  kd<-dplyr::bind_rows(lapply(runs,function(r)env$ivt_ploidy_loglik_df(r,b$objects$fit_data,p$sigma_kary)))
  fd<-dplyr::bind_rows(lapply(runs,function(r)env$ivt_flow_loglik_df(r,b$objects$fit_data,b$cfg$N_UNIT,b$flow_sd)))
  q<--sum(c(mean(a),mean(kd$mean_loglik[is.finite(kd$mean_loglik)]),mean(fd$mean_loglik[is.finite(fd$mean_loglik)])))
  ag<-vapply(names(pred$growth),function(id)dnorm(env$ivt_observed_passage_summary(b$objects$fit_data[[id]])$observed_growth,pred$growth[id],p$sigma_growth,log=TRUE),numeric(1))
  ak<-vapply(names(pred$kary),function(id){obs<-env$ivt_observed_passage_summary(b$objects$fit_data[[id]])$observed_kary;mean(log(pmax(pred$kary[[id]][match(obs,b$core$grid_pre)],1e-12)))},numeric(1))
  af<-vapply(names(pred$flow),function(id){o<-env$ivt_observed_passage_summary(b$objects$fit_data[[id]])$observed_flow;w<-o$density*median(diff(sort(unique(o$ploidy))));w<-w/sum(w);sum(w*log(pmax(pred$flow[[id]],1e-12)))},numeric(1))
  independent<--sum(c(mean(ag),mean(ak),mean(af)))
  stopifnot(is.finite(q),abs(q-independent)<=1e-8)
  list(objective=q,independent=independent,difference=q-independent,growth=mean(a),kary=mean(ak),flow=mean(af))
}
sf_branches <- function(env,b,runs,p,pred) {
  oxy<-sort(unique(unlist(lapply(b$adapters,function(a)vapply(a$segments,function(s)s$oxygen_pct,numeric(1))))))
  masks<-vapply(oxy,function(o){probs<-env$.pmisseg_of_O2(o,p,N=b$core$grid_pre);paste(vapply(seq_along(probs),function(j){n<-b$core$grid_pre[j];paste(as.integer(dbinom(0:n,n,probs[j])>=1e-8),collapse='')},character(1)),collapse=':')},character(1))
  negative_margin<-1-b$cfg$DT*(p$lam_max+p$mu_hp*(b$cfg$N_MAX/44)^p$gamma_mu)
  stock_margin<-1-b$cfg$DT*p$k_clear
  crowd_bounds<-unlist(lapply(runs,function(run)vapply(run$segment_results,function(s){gap<-max(diff(s$segment$obs_days_local));max(s$sim$Ntot_live_obs)*exp(p$lam_max*gap)/b$cfg$K},numeric(1))))
  list(selection_sha=sf_hash(pred$selection),pruning_sha=sf_hash(masks),initial_zero_sha=sf_hash(lapply(runs,function(run)lapply(run$segment_results,function(s)s$sim$live_state_obs[1,]==0))),
    score_floor_sha=sf_hash(list(lapply(pred$kary,function(v)v<1e-12),lapply(pred$flow,function(v)v<1e-12))),
    negative_live_margin=negative_margin,retained_margin=stock_margin,crowding_upper_ratio=max(crowd_bounds),
    clamp_branch_certified=negative_margin>0 && stock_margin>0 && max(crowd_bounds)<1)
}
sf_run <- function(env,b,p,out,id,ledger) {
  stopifnot(!id%in%ledger$eval_id,nrow(ledger)<74)
  all_ledgers<-list.files(file.path(dirname(dirname(out)),'verification'),pattern='^forward_ledger[.]tsv$',recursive=TRUE,full.names=TRUE)
  used<-sum(vapply(all_ledgers,function(f)nrow(sf_tab(f)),integer(1)));stopifnot(used<74L)
  ledger<-rbind(ledger,data.frame(eval_id=id,status='started',parameter_sha=sf_hash(p),seconds=NA_real_));sf_emit(ledger,file.path(out,'forward_ledger.tsv'))
  start<-proc.time()[3]
  error<-NULL
  runs<-tryCatch({z<-lapply(c('2N','4N'),function(c)env$ivt_run_lineage(b$adapters[[c]],cfg=b$cfg,run_params=p,model_core=b$core));names(z)<-c('2N','4N');z},invitro_protocol_infeasible=function(e){error<<-conditionMessage(e);NULL})
  elapsed<-proc.time()[3]-start
  if(is.null(runs)) {
    ledger$status[nrow(ledger)]<-'protocol_infeasible';ledger$seconds[nrow(ledger)]<-elapsed;sf_emit(ledger,file.path(out,'forward_ledger.tsv'))
    return(list(ledger=ledger,result=list(status='protocol_infeasible',error=error)))
  }
  pred<-sf_predictions(env,b,runs,p);score<-sf_score(env,b,runs,p,pred);branch<-sf_branches(env,b,runs,p,pred)
  ledger$status[nrow(ledger)]<-'returned';ledger$seconds[nrow(ledger)]<-elapsed;sf_emit(ledger,file.path(out,'forward_ledger.tsv'))
  list(ledger=ledger,result=list(status='returned',params=p,pred=pred,score=score,branch=branch,runs=runs,seconds=elapsed))
}
