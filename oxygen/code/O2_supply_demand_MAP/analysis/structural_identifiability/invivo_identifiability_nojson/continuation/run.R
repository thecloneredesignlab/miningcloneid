args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3L)
study<-normalizePath(args[1]);mode<-args[2];task<-args[3]
entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
for(file in c('adapter.R','helpers.R','machine.R','optimize.R','observation.R','reference_checks.R','control.R','global.R','test_machine.R'))source(file.path(entry,file))
dispatch<-file.path(study,'dispatch',task);m<-read_tsv(file.path(dispatch,'manifest.tsv'))
stopifnot(all(vapply(file.path(dispatch,m$path),hash_file,character(1))==m$sha256))
spec<-read_kv(file.path(dispatch,'spec.tsv'));spec$task<-task;stopifnot(spec$mode==mode)
contract<-readRDS(file.path(study,'manifests/contract.rds'));spec$fixed<-if(spec$density_index>0L)setNames(contract$rho_grid[spec$density_index],'log10_rho_2N') else NULL
spec$family<-spec$root_family
stopifnot(spec$input_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')),
 spec$release_manifest_sha256==hash_file(file.path(entry,'manifest.tsv')),
 as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))==spec$cores)
deadline<-as.numeric(Sys.getenv('IVI_ALLOCATION_START'))+spec$wall_seconds-spec$margin_seconds
options(ivi.job.deadline=deadline,ivi.root.task=task,ivi.entry=entry,ivi.provenance=spec)
out<-file.path(study,if(mode%in%c('local','global'))'fits' else 'tests',task)
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
tryCatch({
 if(mode=='engineering_parallel')Sys.setenv(MININGCLONEID_RCPP_REBUILD='TRUE')
 ctx<-formal_context(study,spec);Sys.setenv(MININGCLONEID_RCPP_REBUILD='FALSE');reference<-validate_point(ctx,ctx$endpoint)
 old<-readRDS(file.path(study,'reference/production_replay.rds'));expected<-if(spec$policy=='production')old$components$L else old$components$L_data
 stopifnot(abs(reference$components$L-expected)<1e-8)
 write_kv(list(pass=TRUE,difference=reference$components$L-expected),file.path(out,'reference_recheck.tsv'))
 atomic_rds(spec,file.path(out,'dispatch_spec.rds'))
 if(mode=='local')run_local_window(study,task,spec)
 else if(mode=='global')run_global_window(study,entry,task,spec)
 else if(mode%in%c('fresh','directions','profile_fresh'))run_acceptance(study,task,spec,mode!='fresh')
 else if(mode=='replay') {
  roots<-read_tsv(file.path(study,'manifests/roots.tsv'));roots<-roots[roots$policy==spec$policy,,drop=FALSE];gaps<-numeric()
  for(i in seq_len(nrow(roots))) {
   s<-spec;s$root_id<-roots$root_id[i];s$root_family<-roots$root_family[i];s$import_id<-roots$import_id[i];s$selector<-'import';s$task<-paste0(task,'/points/',s$import_id)
   run_import_replay(study,s$task,s);z<-read_kv(file.path(study,'tests',s$task,'result.tsv'));gaps<-c(gaps,z$difference,z$history_repeat_difference)
  }
  stopifnot(max(abs(gaps))<1e-11);write_kv(list(pass=TRUE,points=nrow(roots),maximum_repeat_gap=max(abs(gaps))),file.path(out,'result.tsv'));finalize_payload(out)
 } else if(mode=='engineering_serial') {
  s<-spec;s$max_calls<-160L;s$task<-paste0(task,'/split_first');run_local_window(study,s$task,s,TRUE)
  s$selector<-'checkpoint';s$source_task<-s$task;s$task<-paste0(task,'/split_second');s$max_calls<-200L;run_local_window(study,s$task,s,TRUE)
  second<-readRDS(file.path(study,'tests',s$task,'checkpoint.rds'))$state
  s<-spec;s$task<-paste0(task,'/uninterrupted');s$max_calls<-359L;run_local_window(study,s$task,s,TRUE)
  full<-readRDS(file.path(study,'tests',s$task,'checkpoint.rds'))$state
  stopifnot(identical(full,second));write_kv(list(pass=TRUE,real_checkpoint_sequence_identity=TRUE,method_calls=full$total_calls,free_coordinates=20L-length(spec$fixed)),file.path(out,'result.tsv'));finalize_payload(out)
 } else if(mode=='engineering_parallel') {
  unit<-test_machine();write_kv(unit,file.path(out,'analytic_fixture.tsv'))
  # Real aggregate token exhaustion across eight concurrent workers.
  token<-file.path(out,'token_fixture.rds');atomic_rds(list(issued=0L,limit=37L),token)
  cl<-parallel::makePSOCKcluster(8L);parallel::clusterCall(cl,function(entry){source(file.path(entry,'adapter.R'));source(file.path(entry,'control.R'));TRUE},entry)
  counts<-parallel::clusterCall(cl,function(path){n<-0L;repeat{v<-try(take_token(path),silent=TRUE);if(inherits(v,'try-error'))break;n<-n+1L};n},token)
  parallel::stopCluster(cl);stopifnot(sum(unlist(counts))==37L,readRDS(token)$issued==37L)
  k<-0L
  for(policy in c('production','data_only'))for(init in c('hybrid','uniform'))for(fixed in c(FALSE,TRUE)) {
   k<-k+1L;s<-spec;s$policy<-policy;s$init_mode<-init;s$root_id<-paste0('engineering_',k);s$root_family<-'reference';s$selector<-'reference';s$seed<-52000L+k;s$max_calls<-160L
   s$fixed<-if(fixed)setNames(contract$rho_grid[4],'log10_rho_2N') else NULL;s$task<-paste0(task,'/cases/',policy,'_',init,'_',if(fixed)'fixed' else 'free')
   run_global_window(study,entry,s$task,s,TRUE)
  }
  # One complete and one interrupted native block, followed by a new window.
  s<-spec;s$policy<-'production';s$fixed<-NULL;s$root_id<-'engineering_partial';s$root_family<-'reference';s$selector<-'reference';s$init_mode<-'hybrid';s$max_calls<-131L;s$task<-paste0(task,'/partial_first')
  run_global_window(study,entry,s$task,s,TRUE)
  previous<-s$task;s$selector<-'global_population';s$source_task<-previous;s$max_calls<-80L;s$task<-paste0(task,'/partial_resumed')
  run_global_window(study,entry,s$task,s,TRUE)
  q<-read_kv(file.path(study,'tests',s$task,'result.tsv'));stopifnot(q$root_callbacks==211L)
  requested<-readRDS(file.path(study,'tests',s$task,'blocks/block_000001/actual_initial_population.rds'))$free
  center<-readRDS(file.path(study,'tests',previous,'accepted_candidate.rds'))$record$parameters
  stopifnot(identical(as.numeric(requested[1,]),as.numeric(center)))
  # Exhaust a quota with fewer than NP calls left for the final native block.
  s<-spec;s$policy<-'production';s$fixed<-NULL;s$root_id<-'engineering_short_partial';s$root_family<-'reference';s$selector<-'reference';s$init_mode<-'hybrid';s$max_calls<-97L;s$task<-paste0(task,'/short_partial_first')
  run_global_window(study,entry,s$task,s,TRUE)
  short<-read_kv(file.path(study,'tests',s$task,'blocks/block_000002/initial_population_validation.tsv'))
  stopifnot(!short$initial_population_complete,!short$exact_population_multiset,short$partial_population_subset,short$native_initial_evaluations==17L)
  previous<-s$task;s$selector<-'global_population';s$source_task<-previous;s$max_calls<-80L;s$task<-paste0(task,'/short_partial_resumed')
  run_global_window(study,entry,s$task,s,TRUE)
  q<-read_kv(file.path(study,'tests',s$task,'result.tsv'));stopifnot(q$root_callbacks==177L)
  previous<-s$task;s$selector<-'task';s$source_task<-previous;s$fixed<-setNames(contract$rho_grid[4],'log10_rho_2N');s$max_calls<-80L;s$task<-paste0(task,'/projected_global')
  run_global_window(study,entry,s$task,s,TRUE)
  projected<-readRDS(file.path(study,'tests',s$task,'actual_start.rds'));point<-readRDS(file.path(study,'tests',s$task,'accepted_candidate.rds'))
  stopifnot(is.null(projected$expected),!identical(projected$record$parameters,projected$parameters),point$record$parameters['log10_rho_2N']==s$fixed)
  previous<-s$task;s$source_task<-previous;s$max_calls<-1L;s$task<-paste0(task,'/tiny_population')
  run_global_window(study,entry,s$task,s,TRUE)
  tiny<-read_kv(file.path(study,'tests',s$task,'blocks/block_000001/initial_population_validation.tsv'));tiny_point<-readRDS(file.path(study,'tests',s$task,'accepted_candidate.rds'))
  stopifnot(tiny$native_initial_evaluations==1L,tiny$partial_population_subset,!tiny$initial_population_complete,tiny_point$checked$L<=point$checked$L)
  write_kv(list(pass=TRUE,actual_DE_cases=k,aggregate_token_calls=37L,analytic_beyond_25=TRUE,partial_window_calls=131L,resumed_root_calls=211L,partial_incumbent_retained=TRUE,short_partial_calls=97L,short_population_evaluations=17L,short_resumed_root_calls=177L,projected_global_fixed_residual=0,tiny_partial_callbacks=1L,tiny_incumbent_retained=TRUE),file.path(out,'result.tsv'));finalize_payload(out)
 } else stop('Unknown mode')
},error=function(e){write_kv(list(failed=TRUE,message=conditionMessage(e)),file.path(out,'FAILED.tsv'));stop(e)})
