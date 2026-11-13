# Population warm-restart blocks are serialization intervals, not stopping caps.
run_global_window<-function(study,entry,task,spec,fixture=FALSE) {
 out<-file.path(study,if(fixture)'tests' else 'fits',task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 ctx<-formal_context(study,spec);start<-resolve_start(ctx,spec);token_artifact<-file.path(out,'callback_budget.rds')
 token_directory<-if(dir.exists('/dev/shm')&&file.access('/dev/shm',2)==0)'/dev/shm' else tempdir()
 token<-tempfile(pattern=paste0('ivc_tokens_',Sys.getpid(),'_'),tmpdir=token_directory);sealed<-FALSE
 atomic_rds(list(issued=0L,limit=spec$max_calls),token);Sys.setenv(CONT_TOKEN_PATH=token)
 on.exit({Sys.unsetenv('CONT_TOKEN_PATH');if(!sealed&&file.exists(token))atomic_rds(readRDS(token),token_artifact);unlink(c(token,paste0(token,'.lock')),recursive=TRUE)},add=TRUE)
 population<-NULL;rng<-NULL;old_blocks<-0L;old_calls<-0L
 if(!is.null(start$checkpoint)) {
  # Global checkpoints have a different schema from exact local state.
  stop('Use global_population selector for native population warm restarts')
 }
 if(spec$selector=='global_population') {
  parent<-task_root(study,spec$source_task);cp<-readRDS(file.path(parent,'global_checkpoint.rds'))
  stopifnot(cp$root_id==spec$root_id,cp$policy==spec$policy);population<-cp$population;rng<-cp$rng
  previous<-read_kv(file.path(parent,'result.tsv'));old_calls<-previous$root_callbacks;old_blocks<-previous$root_blocks
 }
 atomic_rds(start,file.path(out,'actual_start.rds'));atomic_rds(spec,file.path(out,'actual_spec.rds'))
 best_start<-start;block<-0L;total_calls<-0L;generation_seconds<-NA_real_;reason<-'time_limit';t<-proc.time()[[3]];summaries<-list()
 repeat {
  if(as.numeric(Sys.time())>=getOption('ivi.job.deadline')){reason<-'time_limit';break}
  budget<-readRDS(token);if(budget$issued>=budget$limit){reason<-'call_limit';break}
  block<-block+1L
  generations<-if(fixture)1L else if(is.na(generation_seconds))4L else max(1L,as.integer(min(600,max(30,getOption('ivi.job.deadline')-as.numeric(Sys.time())))/generation_seconds))
  child<-paste0(task,'/blocks/block_',sprintf('%06d',block));s<-spec;s$task<-child;s$block_generations<-generations
  if(!is.null(population)&&!is.null(best_start$record)) {
   free<-setdiff(names(ctx$endpoint),names(spec$fixed));population[1L,]<-best_start$parameters[free];colnames(population)<-free
  }
  s$initial_population<-population;s$rng_state<-rng;s$start_override<-best_start
  s$reltol<-0;s$steptol<-generations+1L
  before<-readRDS(token)$issued;bt<-proc.time()[[3]]
  run_global(study,entry,child,s,fixture)
  root<-file.path(study,if(fixture)'tests' else 'fits',child);result<-readRDS(file.path(root,'native_result.rds'));a<-readRDS(file.path(root,'accepted_candidate.rds'))
  counts<-readRDS(token)$issued-before;total_calls<-total_calls+counts
  if(is.null(best_start$expected)||a$checked$L<best_start$expected)best_start<-list(parameters=a$record$parameters,record=a$record,expected=a$checked$L,source=child,source_sha256=hash_file(file.path(root,'payload_manifest.tsv')),root_id=spec$root_id)
  complete<-!is.null(result)
  if(complete) {
   population<-result$member$pop;rng<-readRDS(file.path(root,'final_rng_and_counters.rds'))$master
   generation_seconds<-(proc.time()[[3]]-bt)/max(1L,result$optim$iter)
   atomic_rds(list(population=population,rng=rng,root_id=spec$root_id,policy=spec$policy,
    total_blocks=old_blocks+block,total_callbacks=old_calls+readRDS(token)$issued,semantics='population_warm_restart'),file.path(out,'global_checkpoint.rds'))
  }
  summaries[[block]]<-data.frame(block=block,generations_requested=generations,complete=complete,callbacks=counts,objective=a$checked$L,wall_seconds=proc.time()[[3]]-bt)
  if(!complete){reason<-if(readRDS(token)$issued>=spec$max_calls)'call_limit' else 'time_limit';break}
 }
 paths<-list.files(out,pattern='^improvement_.*[.]rds$',recursive=TRUE,full.names=TRUE)
 own_record<-if(!is.null(best_start$record)&&identical(best_start$record$parameters,best_start$parameters))best_start$record else NULL
 best<-if(length(paths))incumbents(out,ctx) else if(!is.null(own_record))own_record else stop('No complete feasible incumbent')
 if(!is.null(own_record)&&own_record$components$L<best$components$L)best<-own_record
 point<-validate_point(ctx,best$parameters,record=best)
 stopifnot(readRDS(token)$issued==total_calls,total_calls<=spec$max_calls)
 atomic_rds(list(record=best,checked=point$components,policy=spec$policy,fixed=spec$fixed,root_id=spec$root_id,root_family=spec$root_family),file.path(out,'accepted_candidate.rds'))
 if(length(summaries))write_tsv(do.call(rbind,summaries),file.path(out,'blocks.tsv'))
 write_kv(list(pass=TRUE,policy=spec$policy,method='DEoptim_population_warm_restart',termination=reason,evaluations=total_calls,
  root_callbacks=old_calls+total_calls,objective=best$components$L,start_objective=if(is.null(start$expected))NA_real_ else start$expected,
  window_gain=if(is.null(start$expected))NA_real_ else start$expected-best$components$L,
  blocks=block,root_blocks=old_blocks+block,wall_seconds=proc.time()[[3]]-t,budget_limited=reason=='time_limit',call_limited=reason=='call_limit',
  free_coordinates=length(best$free_names),fixed_residual=if(length(spec$fixed))max(abs(best$parameters[names(spec$fixed)]-spec$fixed)) else 0),file.path(out,'result.tsv'))
 atomic_rds(readRDS(token),token_artifact)
 write_kv(list(token_filesystem=token_directory,aggregate_callbacks=readRDS(token)$issued),file.path(out,'token_accounting.tsv'))
 finalize_payload(out);sealed<-TRUE
}
