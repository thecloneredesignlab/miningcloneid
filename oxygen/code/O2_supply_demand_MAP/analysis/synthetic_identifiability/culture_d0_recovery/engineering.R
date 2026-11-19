args<-commandArgs(TRUE);stopifnot(length(args)==3L)
s<-normalizePath(args[1]);task<-args[2];mode<-args[3];release<-Sys.getenv('D0_RELEASE')
source(file.path(release,'common.R'))
out<-file.path(s,'engineering',task);stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
started<-proc.time()[3]
tryCatch({
d0_guard(s,out)
if(mode=='prepare') {
 source(file.path(release,'unit_loss.R'));d0_unit_loss(out)
 c<-d0_backend(s,task);fit<-readRDS(file.path(c$paths$wp0,'source/seed144/fit_result.rds'))
 cfg<-fit$cfg;cfg$parameter_table<-file.path(c$paths$wp0,'source/seed144/parameter_table_input.csv')
 allmap<-o2prof_map(cfg$parameter_table,fit$best_params,'zero_extended')
 projected<-o2prof_project_start(allmap,fit$best_params);p<-projected$params
 nuisance<-c('sigma_growth','sigma_kary','init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N')
 map<-o2prof_map(cfg$parameter_table,p,'zero_extended',p[nuisance])
 inputs<-file.path(c$paths$wp0,'source/inputs/oxygen')
 objects<-c$backend$ivt_load_fit_objects(repo_root=inputs,fit_objects_dir=file.path(inputs,'ploidyOxygen/data/fit_objects'),flow_csv_path=file.path(inputs,'data/g0g1_ploidy_density_grid.csv'))
 stopifnot(!isTRUE(objects$death_enabled));observations<-list();design<-list()
 for(id in names(objects$fit_data)) {
  z<-c$backend$ivt_observed_passage_summary(objects$fit_data[[id]])
  flow<-!is.null(z$observed_flow)&&nrow(z$observed_flow)>0L
  observations[[id]]<-list(growth=is.finite(z$observed_growth),kary=length(z$observed_kary)>0L,n_kary=length(z$observed_kary),flow=flow,
     flow_grid=if(flow)as.numeric(z$observed_flow$ploidy)else numeric(),
     step=if(flow)median(diff(sort(unique(z$observed_flow$ploidy))))else NA_real_,sample_name=z$observed_flow_sample_name)
  design[[id]]<-list(g=NA_real_,passage_duration=z$passage_duration,initial_cells=z$initial_cells,final_cells=z$final_cells,kary=numeric(),flow=NULL)
 }
 stopifnot(sum(vapply(observations,function(z)z$growth,logical(1)))==114L,
   sum(vapply(observations,function(z)z$kary,logical(1)))==12L,
   sum(vapply(observations,function(z)z$flow,logical(1)))==20L,
   sum(vapply(observations,function(z)z$n_kary,integer(1)))==220L)
 adapters<-lapply(c('2N','4N'),function(cohort)c$backend$ivt_build_job_table_adapter(objects[[paste0('jobs_',cohort)]],design,cohort));names(adapters)<-c('2N','4N')
 blind_map<-map;active<-blind_map$table[blind_map$table$status=='active',]
 mid<-o2prof_transform((o2prof_transform(active$lower,active)+o2prof_transform(active$upper,active))/2,active,inverse=TRUE)
 for(i in seq_len(nrow(active)))blind_map$base_params[[active$parameter[i]]]<-mid[i]
 bundle<-list(cfg=cfg,map=blind_map,adapters=adapters,observations=observations,
   flow_sd=fit$best_components$sigma_flow_ploidy,scope='Conditional expected engineering; no physical sampling claim',
   synthetic_dataset_id=basename(s),fixed_nuisance=nuisance,source_fit_sha256=d0_hash(file.path(c$paths$wp0,'source/seed144/fit_result.rds')))
 d0_rds(bundle,file.path(s,'design/blind_bundle.rds'));d0_rds(p,file.path(s,'truth/truth_parameters.rds'))
 d0_emit(projected$projection,file.path(out,'requested_actual_truth.tsv'));d0_emit(map$table,file.path(out,'parameter_roles.tsv'))
 d0_emit(data.frame(parameter=o2prof_symbols,value=vapply(p[o2prof_symbols],as.numeric,numeric(1)),role=ifelse(o2prof_symbols%in%nuisance,'known_nuisance','mechanism')),file.path(s,'truth/truth_parameters.tsv'))
 d0_write(list(truth_id='T01_seed144',selection='Prespecified validated all-active production reference; not synthetic outcome',source_sha256=bundle$source_fit_sha256,global_optimum_claim=FALSE,free_mechanisms=14L,known_nuisance=6L),file.path(s,'truth/truth_registry.tsv'))
 c$bundle<-bundle;c$core<-c$backend$build_model_core(cfg=cfg);c$engineering<-TRUE;c$study<-s
 t<-proc.time()[3];runs<-d0_forward(c,p,task);first<-proc.time()[3]-t
 identity<-list()
 for(cohort in c('2N','4N')) {
  actual<-runs[[cohort]]$segment_results;original<-fit$best_components[[paste0('run_',cohort)]]$segment_results
  stopifnot(length(actual)==length(original))
  for(i in seq_along(actual)) {
   stopifnot(isTRUE(all.equal(actual[[i]]$selection,original[[i]]$selection,tolerance=1e-12,check.attributes=TRUE)))
   for(n in intersect(names(actual[[i]]$sim),names(original[[i]]$sim))) {
    if(n%in%c('cache_g_build','cache_g_hit','cache_g_hysteresis'))next
    a<-actual[[i]]$sim[[n]];b<-original[[i]]$sim[[n]]
    if(!is.numeric(a)||!is.numeric(b))next
    stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=TRUE)))
    finite<-is.finite(a)&is.finite(b);identity[[length(identity)+1L]]<-data.frame(cohort=cohort,segment=i,quantity=n,max_scaled_gap=if(any(finite))max(abs(a[finite]-b[finite])/pmax(1,abs(b[finite])))else 0)
   }
  }
 }
 d0_emit(do.call(rbind,identity),file.path(out,'truth_native_identity.tsv'))
 pred<-d0_predict(c,runs,p);data<-pred[c('growth','kary','flow')];data$dataset_id<-basename(s);data$kind<-'deterministic_expected_not_finite_sample'
 d0_rds(data,file.path(s,'data/expected_observations.rds'));c$dataset<-data;c$data_sha<-d0_hash(file.path(s,'data/expected_observations.rds'))
 base<-d0_loss(c,pred);d0_rds(list(params=p,runs=runs,predictions=pred,loss=base),file.path(s,'truth/truth_native.rds'))
 repeat_rows<-list()
 for(i in 1:2){t<-proc.time()[3];z<-d0_predict(c,d0_forward(c,p,task),p);v<-d0_loss(c,z);stopifnot(abs(v$objective-base$objective)<=1e-8,v$max_growth_error<=1e-10,v$max_kary_tv<=1e-10,v$max_flow_tv<=1e-10);repeat_rows[[i]]<-data.frame(repeat_id=i,seconds=proc.time()[3]-t,objective=v$objective,difference=v$objective-base$objective)}
 d0_emit(do.call(rbind,repeat_rows),file.path(out,'serial_replays.tsv'))
 changed<-data;changed$growth[1]<-changed$growth[1]+map$fixed$sigma_growth
 sensitivity<-d0_loss(c,pred,changed)$objective-base$objective
 stopifnot(abs(sensitivity-1/(2*114))<=1e-8)
 d0_write(list(pass=TRUE,synthetic_response_sensitivity=sensitivity,expected_sensitivity=1/(2*114),old_real_responses_in_blind_bundle=FALSE,model_calls_for_sensitivity=0L),file.path(out,'dataset_consumption.tsv'))
 generated<-unlist(lapply(c('design','data','truth'),function(name)list.files(file.path(s,name),recursive=TRUE,full.names=TRUE)))
 generated<-generated[!dir.exists(generated)&basename(generated)!='manifest.tsv']
 d0_emit(data.frame(path=substring(generated,nchar(s)+2L),sha256=vapply(generated,d0_hash,character(1))),file.path(s,'truth/manifest.tsv'))
 d0_write(list(generator_manifest_sha256=d0_hash(file.path(s,'truth/manifest.tsv')),data_sha256=c$data_sha),file.path(out,'generator_binding.tsv'))
 d0_write(list(pass=TRUE,data_sha256=c$data_sha,first_forward_seconds=first,model_forward_calls=3L,free_mechanisms=14L,known_nuisance=6L,truth_expected_objective=base$objective,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
} else if(mode=='parallel') {
 c<-d0_context(s,task,TRUE);deadline<-as.numeric(Sys.time())+1200
 cl<-d0_cluster(s,task,TRUE,deadline,2L)
 waves<-list()
 for(i in 1:2){t<-proc.time()[3];v<-parallel::clusterCall(cl,function(s){p<-readRDS(file.path(s,'truth/truth_parameters.rds'));q<-d0_ev$fn(d0_unit(d0_ev$context,p));list(pid=Sys.getpid(),q=q,calls=d0_ev$state$n)},s);elapsed<-proc.time()[3]-t
  stopifnot(length(v)==8L,max(vapply(v,function(z)z$q,numeric(1)))-min(vapply(v,function(z)z$q,numeric(1)))<=1e-8)
  waves[[i]]<-data.frame(wave=i,calls=8L,wall_seconds=elapsed,allocated_core_seconds=8*elapsed,throughput=8/elapsed)
 }
 # Exercise the actual DEoptim/cluster API on an analytic function only.
 parallel::clusterCall(cl,function(){.GlobalEnv$d0_smoke<-function(x)sum((x-.5)^2);TRUE})
 set.seed(2026111300L);smoke_population<-matrix(runif(256*14),nrow=256)
 smoke_fn<-function(x).GlobalEnv$d0_smoke(x)
 smoke<-DEoptim::DEoptim(smoke_fn,rep(0,14),rep(1,14),control=list(NP=256L,itermax=1L,trace=FALSE,cluster=cl,initialpop=smoke_population))
 stopifnot(is.finite(smoke$optim$bestval),smoke$optim$nfeval>=256L)
 d0_write(list(pass=TRUE,model_calls=0L,analytic_callbacks=smoke$optim$nfeval,native_optimizer='DEoptim',parallel_workers=8L),file.path(out,'native_optimizer_api_fixture.tsv'))
 parallel::stopCluster(cl);d0_emit(do.call(rbind,waves),file.path(out,'parallel_timing.tsv'))
 d0_write(list(pass=TRUE,workers=8L,model_forward_calls=16L,data_sha256=c$data_sha,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
} else if(mode=='engineering_check') {
 c<-d0_context(s,task,TRUE);p<-readRDS(file.path(s,'truth/truth_parameters.rds'));u<-d0_unit(c,p)
 ev<-d0_evaluator(c,out,task,as.numeric(Sys.time())+9000,45L)
 q<-ev$fn(u);truth<-readRDS(file.path(s,'truth/truth_native.rds'));stopifnot(abs(q-truth$loss$objective)<=1e-8)
 proposals<-list(u)
 for(i in seq_along(u))for(h in c(-.02,.02)){v<-u;v[i]<-min(1,max(0,v[i]+h));ev$fn(v);proposals[[length(proposals)+1L]]<-v}
 set.seed(2026111301L);random_bank<-matrix(runif(16L*14L),nrow=16L,ncol=14L);d0_rds(random_bank,file.path(out,'timing_uniform_bank.rds'))
 for(i in 1:16)ev$fn(random_bank[i,])
 stopifnot(ev$state$n==45L,nrow(d0_tab(file.path(s,'accounting/engineering_forward_slots.tsv')))==64L)
 logs<-d0_tab(file.path(out,paste0('evaluations_',Sys.getpid(),'.tsv')));good<-logs$seconds[logs$feasible]
 sources<-d0_tab(file.path(s,'tables/engineering_sources.tsv'))
 prepare_task<-sources$task[match('prepare',sources$mode)];parallel_task<-sources$task[match('parallel',sources$mode)]
 timing<-d0_tab(file.path(s,'engineering',parallel_task,'parallel_timing.tsv'));warm<-timing$throughput[timing$wave==2]
 conservative_seconds<-max(as.numeric(quantile(good,.95)),8/warm)
 forecast_global_hours<-256*1001*conservative_seconds*1.2/3600
 cost_fits<-is.finite(forecast_global_hours)&&forecast_global_hours<=8*(7-300/3600)
 # Register native max generations/callback limits before any blind start.
 global_calls<-256L*1001L;local_calls<-min(20000L,as.integer(floor((7200-300)/conservative_seconds)))
 controls<-list(NP=256L,global_max_generations=1000L,global_max_callbacks=global_calls,global_worker_max_callbacks=ceiling(global_calls/8)+512L,
    global_soft_seconds=24900L,local_max_callbacks=local_calls,local_soft_seconds=6900L,
    generation_seed_A=2026111311L,generation_seed_B=2026111322L,
    max_growth_standardized_error=.01,max_kary_tv=.01,max_flow_tv=.01,truth_selection_by_distance=FALSE)
 d0_write(controls,file.path(s,'tables/search_controls.tsv'))
 d0_write(list(pass=TRUE,cost_fits=cost_fits,data_sha256=c$data_sha,engineering_full_forward_slots=64L,
     measured_feasible_p95_seconds=as.numeric(quantile(good,.95)),parallel_warm_calls_per_second=warm,
     forecast_global_core_hours=forecast_global_hours,forecast_complete_chain_core_hours=forecast_global_hours+2,
     chain_core_hour_cap=58L,wall_seconds=proc.time()[3]-started),file.path(s,'tables/engineering_gate.tsv'))
 evidence<-c('tables/search_controls.tsv','truth/manifest.tsv','data/expected_observations.rds','design/blind_bundle.rds',paste0('engineering/',prepare_task,'/payload_manifest.tsv'),paste0('engineering/',parallel_task,'/payload_manifest.tsv'))
 d0_emit(data.frame(path=evidence,sha256=vapply(file.path(s,evidence),d0_hash,character(1))),file.path(s,'tables/engineering_evidence.tsv'))
 d0_write(list(pass=TRUE,cost_fits=cost_fits,model_forward_calls=45L,forecast_global_core_hours=forecast_global_hours,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
} else stop('Unsupported engineering mode')
d0_guard(s,out);d0_seal(out)
},error=function(e){d0_write(list(failed=TRUE,message=conditionMessage(e),wall_seconds=proc.time()[3]-started),file.path(out,'FAILED.tsv'));stop(e)})
