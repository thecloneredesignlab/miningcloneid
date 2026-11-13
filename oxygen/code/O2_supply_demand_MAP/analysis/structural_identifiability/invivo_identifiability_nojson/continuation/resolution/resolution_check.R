# Explicit fixed-center full multiscale/coupled/observation checks, not an optimizer.
resolution_probes<-function(p,lower,upper,fixed){
 steps<-c(.02,.005,.001,.0001)
 probes<-make_direction_probes(p,lower,upper,fixed,steps)
 for(i in seq_along(probes))probes[[i]]$kind<-'axis'
 active<-setdiff(names(p),c(names(fixed),'log10_rho_2N','log10_sigma_burden'))
 dirs<-source_directions(active);dirs<-dirs[vapply(dirs,function(v)sum(v!=0)>1L,logical(1))]
 for(step in steps)for(name in names(dirs)){
  q<-p;q[active]<-p[active]+step*(upper-lower)[active]*dirs[[name]]
  if(any(q<lower|q>upper))next
  probes[[length(probes)+1L]]<-list(parameter=name,step=step,direction=0L,parameters=q,kind='coupled')
 }
 probes
}
run_resolution_check<-function(study,task,spec){
 out<-file.path(study,'tests',task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 ctx<-formal_context(study,spec);start<-resolve_start(ctx,spec);p<-start$parameters
 atomic_rds(start,file.path(out,'original_source.rds'))
 ctx$provenance$start_sha256<-hash_file(file.path(out,'original_source.rds'))
 fn<-new_evaluator(ctx,out,spec$fixed,'resolution_fixed_center_check',getOption('ivi.job.deadline'))
 free<-setdiff(names(p),names(spec$fixed));center_value<-fn(p[free]);center_record<-.ivi_state$best_record
 expected<-start$expected
 point<-validate_point(ctx,p,expected=expected,record=if(identical(p,start$original_parameters))start$record else center_record)
 stopifnot(abs(point$components$L-center_value)<1e-8,all(p[names(spec$fixed)]==spec$fixed))
 atomic_rds(list(parameters=p,record=center_record,expected=center_value,source=start$source,
  source_sha256=start$source_sha256,original_parameters=start$original_parameters),file.path(out,'actual_start.rds'))
 atomic_rds(list(parameters=p,components=point$components,predictions=point$predictions,
  source_task=start$source,policy=spec$policy),file.path(out,'fresh_acceptance.rds'))
 center<-list(record=center_record,checked=point$components,policy=spec$policy,
  fixed=spec$fixed,root_id=spec$root_id,root_family=spec$root_family,source_task=task,source_scope='tests')
 atomic_rds(center,file.path(out,'checked_center.rds'))
 probes<-resolution_probes(p,ctx$lower,ctx$upper,spec$fixed)
 stopifnot(spec$max_calls>=2L*length(probes)+2L)
 scores<-list();updates<-list();started<-proc.time()[[3]]
 observe<-function(record,label){
  proposal<-observation_proposal(ctx,record,spec$fixed)
  if(max(abs(proposal$parameters-record$parameters))<1e-14)return(invisible(NULL))
  value<-fn(proposal$parameters[free]);stopifnot(abs(value-proposal$expected)<1e-8)
  updates[[length(updates)+1L]]<<-list(parent=label,parameters=proposal$parameters,
    expected=proposal$expected,objective=value)
 }
 observe(center_record,'center')
 for(i in seq_along(probes)){
  q<-probes[[i]]$parameters;value<-fn(q[free]);components<-.ivi_state$last_components
  scores[[i]]<-c(probes[[i]],list(objective=value))
  if(value<center_value&&probes[[i]]$parameter%in%setdiff(names(p),c('log10_rho_2N','log10_sigma_burden'))||
     value<center_value&&probes[[i]]$kind=='coupled'){
   record<-list(parameters=q,components=components,predictions=predictions(ctx,q,TRUE))
   observe(record,paste0('probe_',i))
  }
 }
 best<-.ivi_state$best_record;checked<-validate_point(ctx,best$parameters,record=best)
 atomic_rds(list(record=best,checked=checked$components,policy=spec$policy,fixed=spec$fixed,
  root_id=spec$root_id,root_family=spec$root_family,source_task=task,source_scope='tests'),file.path(out,'accepted_candidate.rds'))
 atomic_rds(scores,file.path(out,'resolution_probes.rds'));atomic_rds(updates,file.path(out,'observation_updates.rds'))
 write_tsv(do.call(rbind,lapply(scores,function(z)data.frame(parameter=z$parameter,kind=z$kind,
  step=z$step,direction=z$direction,objective=z$objective,improvement=center_value-z$objective))),file.path(out,'resolution_scores.tsv'))
 write_kv(list(pass=TRUE,mode=spec$mode,policy_version='resolution01',policy=spec$policy,
  scope='fixed_center_all_four_scales_axis_coupled_observation',complete=TRUE,objective=center_value,
  best_objective=best$components$L,difference=point$components$L-center_value,
  source_projection=any(p!=start$original_parameters),probes=length(probes),
  conditional_observation_updates=length(updates),evaluations=.ivi_state$n,
  directional_improvement=center_value-best$components$L,
  representative_gap=center_value-best$components$L,free_coordinates=length(free),
  fixed_residual=if(length(spec$fixed))max(abs(p[names(spec$fixed)]-spec$fixed)) else 0,
  wall_seconds=proc.time()[[3]]-started),file.path(out,'result.tsv'))
 finalize_payload(out)
}
