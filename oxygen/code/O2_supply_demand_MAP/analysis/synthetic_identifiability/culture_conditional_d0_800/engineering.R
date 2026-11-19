args<-commandArgs(TRUE);stopifnot(length(args)==3L)
s<-normalizePath(args[1]);task<-args[2];mode<-args[3];release<-Sys.getenv('D0_RELEASE')
source(file.path(release,'common.R'))
out<-file.path(s,'engineering',task);stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
started<-proc.time()[3]
tryCatch({
if(mode=='prepare') {
 source(file.path(release,'unit_loss.R'));d0_unit_loss(out)
 c<-d0_backend(s,task);fit_path<-file.path(c$paths$wp0,'source/seed144/fit_result.rds');fit<-readRDS(fit_path)
 cfg<-fit$cfg;source_table<-file.path(c$paths$wp0,'source/seed144/parameter_table_input.csv')
 allmap<-o2prof_map(source_table,fit$best_params,'zero_extended')
 projected<-o2prof_project_start(allmap,fit$best_params);p<-projected$params
 nuisance<-c('sigma_growth','sigma_kary','init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N')
 map<-o2prof_map(source_table,p,'zero_extended',p[nuisance])
 map$table$reason[map$table$status=='active']<-'Free mechanism in known-nuisance conditional recovery'
 map$table$reason[map$table$status=='fixed']<-'Known generating nuisance; not estimated in D0'
 # Remove every active generating value from the bundle used by blind fits.
 blind_map<-map;active<-map$table[map$table$status=='active',]
 mid<-o2prof_transform((o2prof_transform(active$lower,active)+o2prof_transform(active$upper,active))/2,active,inverse=TRUE)
 for(i in seq_len(nrow(active)))blind_map$base_params[[active$parameter[i]]]<-mid[i]
 tab<-read.csv(source_table,stringsAsFactors=FALSE)
 for(i in seq_len(nrow(active)))tab$init_value[tab$param_symbol==active$parameter[i]]<-mid[i]
 write.csv(tab,file.path(s,'design/neutral_parameter_table.csv'),row.names=FALSE)
 cfg$parameter_table<-file.path(s,'design/neutral_parameter_table.csv')
 # This registration precedes any biological simulation or expected outcome.
 d0_rds(fit$best_params,file.path(s,'truth/raw_generating_parameters.rds'))
 d0_rds(p,file.path(s,'truth/truth_parameters.rds'));d0_rds(cfg,file.path(s,'truth/truth_configuration.rds'))
 d0_emit(projected$projection,file.path(s,'truth/requested_actual_truth.tsv'))
 d0_emit(map$table,file.path(s,'truth/parameter_roles.tsv'))
 d0_emit(data.frame(parameter=names(p),raw=vapply(fit$best_params,as.character,character(1)),effective=vapply(p,as.character,character(1)),role=ifelse(names(p)%in%nuisance,'known_nuisance',ifelse(names(p)%in%active$parameter,'free_mechanism','fixed_runtime'))),file.path(s,'truth/truth_parameters.tsv'))
 d0_write(list(truth_id='T01_seed144',selection='Prespecified validated production seed144, not selected by recovery',source_sha256=d0_hash(fit_path),effective_truth_sha256=d0_hash(file.path(s,'truth/truth_parameters.rds')),global_optimum_claim=FALSE,free_mechanisms=14L,known_nuisance=6L,registered_before_simulation=TRUE),file.path(s,'truth/truth_registry.tsv'))
 d0_emit(data.frame(key=names(cfg),value=vapply(cfg,as.character,character(1))),file.path(s,'truth/truth_configuration.tsv'))
 inputs<-file.path(c$paths$wp0,'source/inputs/oxygen')
 objects<-c$backend$ivt_load_fit_objects(repo_root=inputs,fit_objects_dir=file.path(inputs,'ploidyOxygen/data/fit_objects'),flow_csv_path=file.path(inputs,'data/g0g1_ploidy_density_grid.csv'))
 stopifnot(!isTRUE(objects$death_enabled));observations<-list();design<-list()
 for(id in names(objects$fit_data)) {
  z<-c$backend$ivt_observed_passage_summary(objects$fit_data[[id]]);flow<-!is.null(z$observed_flow)&&nrow(z$observed_flow)>0L
  observations[[id]]<-list(growth=is.finite(z$observed_growth),kary=length(z$observed_kary)>0L,n_kary=length(z$observed_kary),flow=flow,flow_grid=if(flow)as.numeric(z$observed_flow$ploidy)else numeric(),step=if(flow)median(diff(sort(unique(z$observed_flow$ploidy))))else NA_real_,sample_name=z$observed_flow_sample_name)
  design[[id]]<-list(g=NA_real_,passage_duration=z$passage_duration,initial_cells=z$initial_cells,final_cells=z$final_cells,kary=numeric(),flow=NULL)
 }
 stopifnot(length(observations)==131L,sum(vapply(observations,function(z)z$growth,logical(1)))==114L,sum(vapply(observations,function(z)z$kary,logical(1)))==12L,sum(vapply(observations,function(z)z$flow,logical(1)))==20L,sum(vapply(observations,function(z)z$n_kary,integer(1)))==220L)
 adapters<-lapply(c('2N','4N'),function(cohort)c$backend$ivt_build_job_table_adapter(objects[[paste0('jobs_',cohort)]],design,cohort));names(adapters)<-c('2N','4N')
 bundle<-list(cfg=cfg,map=blind_map,adapters=adapters,observations=observations,flow_sd=fit$best_components$sigma_flow_ploidy,scope='Production-conditioned expected observations; physical finite sampling remains held',synthetic_dataset_id=basename(s),fixed_nuisance=nuisance,source_fit_sha256=d0_hash(fit_path))
 d0_rds(bundle,file.path(s,'design/blind_bundle.rds'))
 d0_emit(do.call(rbind,lapply(names(observations),function(id){o<-observations[[id]];data.frame(passage_id=id,growth=o$growth,kary=o$kary,n_kary_design_only=o$n_kary,flow=o$flow,flow_grid_points=length(o$flow_grid),step=o$step,live_dead_scored=FALSE)})),file.path(s,'design/observation_manifest.tsv'))
 c$bundle<-bundle;c$core<-c$backend$build_model_core(cfg=cfg);c$engineering<-TRUE;c$study<-s
 t<-proc.time()[3];runs<-d0_forward(c,p,task);first<-proc.time()[3]-t;identity<-list()
 for(cohort in c('2N','4N')) {
  actual<-runs[[cohort]]$segment_results;original<-fit$best_components[[paste0('run_',cohort)]]$segment_results;stopifnot(length(actual)==length(original))
  for(i in seq_along(actual)) {
   stopifnot(isTRUE(all.equal(actual[[i]]$selection,original[[i]]$selection,tolerance=1e-12,check.attributes=TRUE)))
   for(n in intersect(names(actual[[i]]$sim),names(original[[i]]$sim))) {
    if(n%in%c('cache_g_build','cache_g_hit','cache_g_hysteresis'))next
    a<-actual[[i]]$sim[[n]];b<-original[[i]]$sim[[n]];if(!is.numeric(a)||!is.numeric(b))next
    stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=TRUE)))
    finite<-is.finite(a)&is.finite(b);identity[[length(identity)+1L]]<-data.frame(cohort=cohort,segment=i,quantity=n,max_scaled_gap=if(any(finite))max(abs(a[finite]-b[finite])/pmax(1,abs(b[finite])))else 0)
   }
  }
 }
 d0_emit(do.call(rbind,identity),file.path(out,'truth_native_identity.tsv'))
 pred<-d0_predict(c,runs,p);data<-pred[c('growth','kary','flow')];data$dataset_id<-basename(s);data$kind<-'deterministic_expected_not_finite_sample'
 d0_rds(data,file.path(s,'data/expected_observations.rds'));c$dataset<-data;c$data_sha<-d0_hash(file.path(s,'data/expected_observations.rds'));base<-d0_loss(c,pred)
 d0_rds(list(params=p,runs=runs,predictions=pred,loss=base),file.path(s,'truth/truth_native.rds'));d0_emit(pred$selection,file.path(s,'truth/truth_selected_days.tsv'))
 # The weighted empirical-cell bridge must reproduce the frozen cell mean.
 bridge<-lapply(names(data$kary),function(id){obs<-c$backend$ivt_observed_passage_summary(objects$fit_data[[id]])$observed_kary;grid<-cfg$N_MIN:cfg$N_MAX;freq<-tabulate(match(obs,grid),nbins=length(grid))/length(obs);pr<-data$kary[[id]];native<-mean(log(pmax(pr[match(obs,grid)],1e-12)));weighted<-sum(freq*log(pmax(pr,1e-12)));stopifnot(abs(native-weighted)<=1e-12);data.frame(passage_id=id,n=length(obs),native=native,weighted=weighted,difference=native-weighted)})
 d0_emit(do.call(rbind,bridge),file.path(out,'empirical_kary_bridge.tsv'))
 changed<-data;changed$growth[1]<-changed$growth[1]+map$fixed$sigma_growth
 sensitivity<-d0_loss(c,pred,changed)$objective-base$objective;stopifnot(abs(sensitivity-1/(2*114))<=1e-8)
 d0_write(list(pass=TRUE,synthetic_response_sensitivity=sensitivity,expected_sensitivity=1/(2*114),real_responses_in_blind_bundle=FALSE),file.path(out,'dataset_consumption.tsv'))
 source(file.path(c$paths$archive,'implementation/diagnostics.R'),local=.GlobalEnv)
 oxygen<-sort(unique(unlist(lapply(adapters,function(a)vapply(a$segments,function(z)z$oxygen_pct,numeric(1))))))
 ledger<-lapply(oxygen,function(o)pilotdiag_rates(c,p,cfg,o)$rates);d0_emit(do.call(rbind,ledger),file.path(s,'truth/latent_event_audit.tsv'))
 generated<-unlist(lapply(c('design','data','truth'),function(name)list.files(file.path(s,name),recursive=TRUE,full.names=TRUE)));generated<-generated[!dir.exists(generated)&basename(generated)!='manifest.tsv']
 d0_emit(data.frame(path=substring(generated,nchar(s)+2L),sha256=vapply(generated,d0_hash,character(1))),file.path(s,'truth/manifest.tsv'))
 d0_write(list(pass=TRUE,data_sha256=c$data_sha,first_forward_seconds=first,model_forward_calls=1L,truth_expected_objective=base$objective,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
} else if(mode%in%c('replay','replay_final')) {
 c<-d0_context(s,task,TRUE);truth<-readRDS(file.path(s,'truth/truth_native.rds'));p<-truth$params
 runs<-d0_forward(c,p,task);pred<-d0_predict(c,runs,p);loss<-d0_loss(c,pred)
 for(cohort in c('2N','4N'))for(i in seq_along(runs[[cohort]]$segment_results)) {
  a<-runs[[cohort]]$segment_results[[i]];b<-truth$runs[[cohort]]$segment_results[[i]]
  stopifnot(isTRUE(all.equal(a$selection,b$selection,tolerance=1e-12,check.attributes=TRUE)))
  for(n in intersect(names(a$sim),names(b$sim)))if(!n%in%c('cache_g_build','cache_g_hit','cache_g_hysteresis')&&is.numeric(a$sim[[n]]))stopifnot(isTRUE(all.equal(a$sim[[n]],b$sim[[n]],tolerance=1e-12,check.attributes=TRUE)))
 }
 stopifnot(abs(loss$objective-truth$loss$objective)<=1e-8,isTRUE(all.equal(pred,truth$predictions,tolerance=1e-12,check.attributes=TRUE)))
 d0_write(c(loss,list(pass=TRUE,data_sha256=c$data_sha,model_forward_calls=1L,wall_seconds=proc.time()[3]-started)),file.path(out,'result.tsv'))
} else if(mode=='parallel') {
 c<-d0_context(s,task,TRUE);cl<-d0_cluster(s,task,TRUE,as.numeric(Sys.time())+1200,2L);waves<-list()
 for(i in 1:2) {
  t<-proc.time()[3];v<-parallel::clusterCall(cl,function(s){p<-readRDS(file.path(s,'truth/truth_parameters.rds'));q<-d0_ev$fn(d0_unit(d0_ev$context,p));list(pid=Sys.getpid(),q=q,calls=d0_ev$state$n)},s);elapsed<-proc.time()[3]-t
  stopifnot(length(v)==7L,max(vapply(v,function(z)z$q,numeric(1)))-min(vapply(v,function(z)z$q,numeric(1)))<=1e-8)
  waves[[i]]<-data.frame(wave=i,calls=7L,wall_seconds=elapsed,allocated_core_seconds=8*elapsed,throughput=7/elapsed)
 }
 parallel::clusterCall(cl,function(){.GlobalEnv$d0_smoke<-function(x)sum((x-.5)^2);TRUE})
 set.seed(2026111300L);population<-matrix(runif(256*14),nrow=256)
 smoke_fn<-function(x).GlobalEnv$d0_smoke(x)
 smoke<-DEoptim::DEoptim(smoke_fn,rep(0,14),rep(1,14),control=list(NP=256L,itermax=1L,trace=FALSE,cluster=cl,initialpop=population));stopifnot(is.finite(smoke$optim$bestval),smoke$optim$nfeval>=256L)
 d0_rds(smoke,file.path(out,'analytic_de_fixture.rds'));parallel::stopCluster(cl)
 d0_emit(do.call(rbind,waves),file.path(out,'parallel_timing.tsv'))
 d0_write(list(pass=TRUE,workers=7L,allocated_cpus=8L,model_forward_calls=14L,data_sha256=c$data_sha,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
} else if(mode=='engineering_check') {
 c<-d0_context(s,task,TRUE);p<-readRDS(file.path(s,'truth/truth_parameters.rds'));u<-d0_unit(c,p)
 ev<-d0_evaluator(c,out,task,as.numeric(Sys.time())+6500,29L,save_native=FALSE);q<-ev$fn(u);rows<-list()
 for(i in seq_along(u))for(h in c(-.02,.02)){v<-u;v[i]<-min(1,max(0,v[i]+h));value<-ev$fn(v);rows[[length(rows)+1L]]<-data.frame(parameter=c$spec$param_name[i],requested_step=h,actual_step=v[i]-u[i],objective=value)}
 d0_emit(do.call(rbind,rows),file.path(out,'coordinate_perturbations.tsv'))
 stopifnot(ev$state$n==29L)
 logs<-d0_tab(file.path(out,paste0('evaluations_',Sys.getpid(),'.tsv')));good<-logs$seconds[logs$feasible]
 measured<-max(as.numeric(quantile(good,.95)),1e-3)
 d0_write(list(pass=TRUE,cost_fits=is.finite(measured)&&measured*256*301/7<41400,measured_forward_p95_seconds=measured,minimum300_forecast_core_hours=8*measured*256*301/7/3600,truth_expected_objective=q,wall_seconds=proc.time()[3]-started),file.path(s,'tables/engineering_gate.tsv'))
 d0_write(list(pass=TRUE,model_forward_calls=29L,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
} else if(mode=='local_api_check') {
 source(file.path(release,'pattern.R'));d0_emit(d0_pattern_fixtures(),file.path(out,'analytic_local_stopping_fixtures.tsv'))
 d0_write(list(pass=TRUE,model_forward_calls=0L,registered_directional_tolerance=1e-7,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
} else stop('Unsupported engineering mode')
d0_guard(s,out);d0_seal(out)
},error=function(e){d0_write(list(failed=TRUE,message=conditionMessage(e),wall_seconds=proc.time()[3]-started),file.path(out,'FAILED.tsv'));stop(e)})
