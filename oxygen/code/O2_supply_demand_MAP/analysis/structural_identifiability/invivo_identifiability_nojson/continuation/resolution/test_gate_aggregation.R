# Synthetic metadata aggregation fixture. No simulator/objective is loaded or called.
args<-commandArgs(trailingOnly=TRUE);entry<-normalizePath(if(length(args)>1L)args[2] else args[1])
adapter<-if(file.exists(file.path(entry,'adapter.R')))file.path(entry,'adapter.R') else file.path(dirname(entry),'adapter.R')
source(adapter);source(file.path(entry,'resolution_rules.R'))
# Repair work after a completed batch cannot retroactively earn batch two.
frozen<-resolution_completed_batch(14,13.9981,13,13)
stopifnot(!frozen$eligible,
 !resolution_completed_batch(14,13.99,13,12.99,frozen)$eligible,
 resolution_completed_batch(14,13.9979,13,13)$eligible,
 !resolution_completed_batch(14,13.998,13,13)$eligible)
cat('Completed batch outcome is frozen across later repair gains: PASS\n')
for(case in c('unresolved','recovered','candidate_representative_edge')){
 recovered<-case!='unresolved';edge<-case=='candidate_representative_edge'
 study<-tempfile('resolution_gate_fixture_');dir.create(study);dir.create(file.path(study,'tables'));dir.create(file.path(study,'manifests'))
 writeLines('synthetic input identity only',file.path(study,'manifests/input_manifest.tsv'))
 jobs<-list();files<-character();file_values<-list();evidence<-character()
 root_for<-function(task)file.path(study,'payloads',task)
 closed<-function(task)task%in%jobs$task
 candidate<-function(task)readRDS(file.path(root_for(task),'accepted_candidate.rds'))
 payload<-function(task)root_for(task)
 prediction_proof<-function(task)invisible(TRUE)
 last_emitted<-NULL;emit<-function(name,pass,extra=list())last_emitted<<-pass
 for(policy in c('production','data_only'))for(i in seq_len(if(policy=='production')3L else 4L)){
  family<-if(i==1L)'archive' else 'uniform';score<-if(policy=='production')14 else 13
  score<-score+if(i==1L)0 else if(edge){if(i==2L).00203 else i/10} else if(recovered).001 else i/10
  task<-paste0('fit_',policy,'_',i);check<-paste0('check_',task);rid<-paste0(policy,'_',family,'_',i)
  record<-list(policy=policy,fixed=NULL,parameters=setNames(rep(.5,20),paste0('p',1:20)),
    components=list(L=score),provenance=list(task=task,root_id=rid,root_family=family))
  point<-list(record=record,checked=list(L=score),policy=policy,root_id=rid,root_family=family,source_task=task,source_scope='fits')
  for(t in c(task,check)){
   dir.create(root_for(t),recursive=TRUE);atomic_rds(point,file.path(root_for(t),'accepted_candidate.rds'));writeLines('fixture',file.path(root_for(t),'payload_manifest.tsv'))
  }
  if(edge&&i==2L){
   improved<-point;improved$record$components$L<-improved$checked$L<-score-.00005
   atomic_rds(improved,file.path(root_for(check),'accepted_candidate.rds'))
  }
  write_kv(list(termination='call_limit',window_gain=.0025,method_calls=20000L),file.path(root_for(task),'result.tsv'))
  write_kv(list(pass=TRUE,complete=TRUE,scope='fixed_center_all_four_scales_axis_coupled_observation',
    difference=0,fixed_residual=0,free_coordinates=20L,directional_improvement=if(edge&&i==2L).00005 else 0),file.path(root_for(check),'result.tsv'))
  atomic_rds(point,file.path(root_for(check),'checked_center.rds'))
  dir.create(file.path(study,'dispatch',check),recursive=TRUE);write_kv(list(checked_fit=task),file.path(study,'dispatch',check,'spec.tsv'))
  for(t in c(task,check))jobs[[t]]<-data.frame(task=t,mode=if(t==task)'local' else 'resolution_check',
    state='COMPLETED',exit_code='0:0',kind='baseline',policy=policy,root_id=rid,root_family=family)
  f<-file.path(root_for(task),'improvement_00000001.rds');atomic_rds(record,f);files<-c(files,f)
  file_values[[f]]<-if(policy=='production')c(production=score,data_only=score+.5) else c(production=score+2,data_only=score)
 }
 jobs<-do.call(rbind,jobs);floors<-c(production=14,data_only=13)
 fragment<-readLines(file.path(entry,'gate_tail.R'))
 fragment[startsWith(fragment,"} else stop('Unsupported metadata gate')")]<-'}'
 eval(parse(text=c('if(TRUE){',fragment)),envir=environment())
 p<-read_tsv(file.path(study,'tables/baseline_policy_gates.tsv'))
 stopifnot(nrow(p)==2L,all(p$pass==recovered),identical(last_emitted,recovered),
   nrow(read_tsv(file.path(study,'tables/resolution_local_decisions.tsv')))==7L,
   nrow(read_tsv(file.path(study,'tables/resolution_anchors.tsv')))==2L)
 if(edge)stopifnot(all(abs(p$uniform_gap-.00198)<1e-12),all(p$uniform_representative_gap>.002))
 unlink(study,recursive=TRUE)
 cat('Synthetic unequal3+4rows, actualminimum/representative/anchor case=',case,': PASS\n')
}
cat('Zero model calls.\n')
