# Owned continuation control; the numerical model/observation functions are unchanged.
take_token<-function(path) {
 lock<-paste0(path,'.lock')
 while(!dir.create(lock,showWarnings=FALSE))Sys.sleep(.001)
 on.exit(unlink(lock,recursive=TRUE))
 z<-readRDS(path);if(z$issued>=z$limit)stop('IVI_CALL_LIMIT')
 z$issued<-z$issued+1L;atomic_rds(z,path);invisible(z$issued)
}
formal_context<-function(study,spec) {
 ctx<-load_contract(study,TRUE,spec$policy);ctx$provenance<-getOption('ivi.provenance',list())
 ctx$provenance$task<-spec$task;ctx$provenance$role<-spec$role;ctx$provenance$family<-spec$root_family
 ctx$provenance$root_id<-spec$root_id;ctx$provenance$root_family<-spec$root_family;ctx$provenance$seed<-spec$seed
 ctx
}
task_root<-function(study,task) {
 a<-file.path(study,'fits',task);b<-file.path(study,'tests',task)
 if(dir.exists(a))a else b
}
resolve_start<-function(ctx,spec) {
 if(!is.null(spec$start_override))return(spec$start_override)
 if(spec$selector=='reference')return(list(parameters=ctx$endpoint,record=NULL,expected=NULL,source='canonical_reference',source_sha256='',root_id=spec$root_id))
 checkpoint<-NULL
 if(spec$selector=='import') {
  path<-file.path(ctx$study,'reference/imports',paste0(spec$import_id,'.rds'));im<-readRDS(path)
  stopifnot(im$target_policy==spec$policy,im$root_id==spec$root_id);r<-im$original_record;expected<-im$expected_objective;origin<-paste0('import:',spec$import_id)
 } else if(spec$selector=='selection') {
  path<-file.path(ctx$study,'dispatch',getOption('ivi.root.task'),'selection.rds');a<-readRDS(path)
  stopifnot(a$policy==spec$policy);r<-a$record;expected<-r$components$L;origin<-a$source_task
  if(spec$root_family!='profile')stopifnot(a$root_id==spec$root_id)
 } else {
  parent<-task_root(ctx$study,spec$source_task);verify_payload_R(parent)
  a<-readRDS(file.path(parent,'accepted_candidate.rds'));r<-a$record;expected<-r$components$L;origin<-spec$source_task;path<-file.path(parent,'payload_manifest.tsv')
  stopifnot(a$policy==spec$policy,r$source_manifest_sha256==hash_file(file.path(ctx$study,'manifests/input_manifest.tsv')))
  if(spec$root_family!='profile')stopifnot(r$provenance$root_id==spec$root_id)
  if(spec$selector=='checkpoint') {
   checkpoint<-readRDS(file.path(parent,'checkpoint.rds'))
   stopifnot(checkpoint$root_id==spec$root_id,checkpoint$policy==spec$policy,
    checkpoint$machine_sha256==hash_file(file.path(getOption('ivi.entry'),'machine.R')),
    identical(checkpoint$state$best_p,r$parameters))
  }
 }
 p<-r$parameters;stopifnot(identical(names(p),names(ctx$endpoint)),all(p>=ctx$lower),all(p<=ctx$upper))
 original<-p;if(length(spec$fixed))p[names(spec$fixed)]<-spec$fixed
 if(!identical(p,original)){expected<-NULL;checkpoint<-NULL}
 list(parameters=p,original_parameters=original,record=r,expected=expected,source=origin,source_sha256=hash_file(path),root_id=spec$root_id,checkpoint=checkpoint)
}
run_local_window<-function(study,task,spec,fixture=FALSE) {
 out<-file.path(study,if(fixture)'tests' else 'fits',task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 ctx<-formal_context(study,spec);start<-resolve_start(ctx,spec);p<-start$parameters;fixed<-spec$fixed;free<-setdiff(names(p),names(fixed))
 atomic_rds(spec,file.path(out,'actual_spec.rds'));atomic_rds(list(parameters=p,resolution=start),file.path(out,'actual_local_start.rds'))
 ctx$provenance$start_sha256<-hash_file(file.path(out,'actual_local_start.rds'))
 fn<-new_evaluator(ctx,out,fixed,'cyclic_block_pattern',getOption('ivi.job.deadline'))
 t<-proc.time()[[3]];first<-fn(p[free]);if(!is.null(start$expected))stopifnot(abs(first-start$expected)<1e-8)
 write_kv(list(pass=TRUE,source_score=if(is.null(start$expected))NA_real_ else start$expected,first_score=first,
  difference=if(is.null(start$expected))NA_real_ else first-start$expected,projection=any(p!=start$original_parameters)),file.path(out,'start_identity.tsv'))
 st<-if(is.null(start$checkpoint))machine_new(p,first,ctx$lower,ctx$upper,fixed) else start$checkpoint$state
 stopifnot(abs(st$best_value-first)<1e-8,identical(st$best_p,p))
 old_calls<-if(is.null(start$checkpoint))0L else start$checkpoint$root_callbacks
 checkpoint_index<-0L;proposals<-list()
 save<-function(state) {
  checkpoint_index<<-checkpoint_index+1L
  z<-list(state=state,record=.ivi_state$best_record,root_id=spec$root_id,policy=spec$policy,
   root_callbacks=old_calls+.ivi_state$n,machine_sha256=hash_file(file.path(getOption('ivi.entry'),'machine.R')),
   parent_task=spec$source_task,task=task,input_manifest_sha256=hash_file(file.path(study,'manifests/input_manifest.tsv')))
  path<-file.path(out,'checkpoints',sprintf('checkpoint_%08d.rds',checkpoint_index));atomic_rds(z,path)
  write_kv(list(path=basename(path),sha256=hash_file(path)),file.path(out,'checkpoints/CURRENT.tsv'))
 }
 evaluate<-function(q,kind,step,direction,sequence) {
  value<-fn(q[free]);if(kind=='observation')proposals[[length(proposals)+1L]]<<-list(parameters=q,actual=value,sequence=sequence)
  value
 }
 propose<-function(p,v){stopifnot(max(abs(p-.ivi_state$best_record$parameters))<1e-12);observation_proposal(ctx,.ivi_state$best_record,fixed)}
 done<-function(poll)atomic_rds(poll,file.path(out,'polls',sprintf('poll_%08d.rds',poll$poll)))
 result<-machine_run(st,evaluate,propose,max(0L,spec$max_calls-1L),getOption('ivi.job.deadline'),save,done)
 best<-.ivi_state$best_record;point<-validate_point(ctx,best$parameters,record=best)
 stopifnot(abs(result$state$best_value-best$components$L)<1e-8,all(best$parameters[names(fixed)]==fixed))
 latest<-read_kv(file.path(out,'checkpoints/CURRENT.tsv'));cp<-readRDS(file.path(out,'checkpoints',latest$path));atomic_rds(cp,file.path(out,'checkpoint.rds'))
 atomic_rds(list(record=best,checked=point$components,policy=spec$policy,fixed=fixed,root_id=spec$root_id,root_family=spec$root_family),file.path(out,'accepted_candidate.rds'))
 atomic_rds(proposals,file.path(out,'observation_proposals.rds'));atomic_rds(result$state,file.path(out,'native_result.rds'))
 write_kv(list(pass=TRUE,policy=spec$policy,method='cyclic_block_pattern',termination=result$reason,objective=best$components$L,
  start_objective=first,window_gain=first-best$components$L,wall_seconds=proc.time()[[3]]-t,evaluations=.ivi_state$n,
  root_callbacks=old_calls+.ivi_state$n,method_calls=result$state$total_calls,visits_max=max(result$state$visits),
  completed_polls=result$state$polls,cycle=result$state$cycle,improving=result$improved,
  fixed_residual=if(length(fixed))max(abs(best$parameters[names(fixed)]-fixed)) else 0,
  free_coordinates=length(free),budget_limited=result$reason=='time_limit',call_limited=result$reason=='call_limit',native_convergence=NA_integer_),file.path(out,'result.tsv'))
 finalize_payload(out)
}
