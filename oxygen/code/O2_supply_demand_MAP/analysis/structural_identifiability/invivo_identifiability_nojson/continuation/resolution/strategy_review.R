# Finite policy decisions over validated records; no model calls.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'))
source(file.path(entry,'resolution_rules.R'))
jobs<-read_tsv(file.path(study,'tables/jobs.tsv'));slots<-read_tsv(file.path(study,'tables/resolution_assignments.tsv'))
policies<-read_tsv(file.path(study,'tables/baseline_policy_gates.tsv'))
decisions_path<-file.path(study,'tables/resolution_local_decisions.tsv');d<-if(file.exists(decisions_path))read_tsv(decisions_path) else NULL
if(!is.null(d))d<-d[!duplicated(d$task,fromLast=TRUE),,drop=FALSE]
bpath<-file.path(study,'tables/resolution_batches.tsv');batches<-if(file.exists(bpath))read_tsv(bpath) else NULL
outcome_path<-file.path(study,'tables/resolution_batch_outcomes.tsv')
outcomes<-if(file.exists(outcome_path))read_tsv(outcome_path) else NULL
closed<-function(j)sub('[+ ].*$','',j$state)%in%c('COMPLETED','FAILED','TIMEOUT','CANCELLED','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')
freeze_completed_outcome<-function(policy,b,u,p){
 batch<-as.integer(b$batch)
 saved<-if(is.null(outcomes))NULL else outcomes[outcomes$policy==policy&outcomes$batch==batch,,drop=FALSE]
 if(!is.null(saved)&&!nrow(saved))saved<-NULL
 outcome<-resolution_completed_batch(b$before_uniform,min(u$verified_best),b$before_feasible,p$known_feasible_bound,saved)
 if(is.null(saved)){
  stopifnot(length(args)>=2L,grepl('^c_[0-9a-f]+$',args[2]))
  collection<-args[2];collection_manifest<-paste0('collections/',collection,'/manifest.tsv')
  directory<-file.path(study,'tables/resolution_completed_batches',paste0(policy,'_c',sprintf('%04d',batch)))
  stopifnot(!dir.exists(directory));dir.create(directory,recursive=TRUE)
  own<-u[which.min(u$verified_best),,drop=FALSE]
  for(pair in list(c(paste0('baseline_',policy,'.rds'),'feasible_candidate.rds'),c(paste0('root_',own$root_id,'_selected.rds'),'uniform_candidate.rds')))
   stopifnot(file.copy(file.path(study,'tables',pair[1]),file.path(directory,pair[2])))
  evidence<-c(substring(file.path(directory,c('feasible_candidate.rds','uniform_candidate.rds')),nchar(study)+2L),collection_manifest)
  evidence_path<-paste0(substring(directory,nchar(study)+2L),'/evidence.tsv')
  write_tsv(data.frame(path=evidence,sha256=vapply(file.path(study,evidence),hash_file,character(1))),file.path(study,evidence_path))
  outcome<-cbind(data.frame(policy=policy,batch=batch),outcome,
   data.frame(uniform_root=own$root_id,collection=collection,collection_manifest_sha256=hash_file(file.path(study,collection_manifest)),
    evidence_manifest=evidence_path,recorded_utc=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC')))
  outcomes<<-if(is.null(outcomes))outcome else rbind(outcomes,outcome)
  write_tsv(outcomes,outcome_path)
 }
 proof<-read_tsv(file.path(study,outcome$evidence_manifest))
 stopifnot(all(vapply(file.path(study,proof$path),hash_file,character(1))==proof$sha256))
 outcome
}
result<-list()
for(policy in c('production','data_only')){
 p<-policies[policies$policy==policy,,drop=FALSE];choice<-'await_active_parents';reason<-'Active fitting/checking or missing complete endpoint evidence';batch<-0L
 ix<-match(slots$task,jobs$task);a<-slots[!is.na(ix)&jobs$policy[ix]==policy&slots$phase=='baseline',,drop=FALSE]
 fit<-a[a$slot_class%in%c('existing_local','baseline_global','baseline_polish','baseline_repair'),,drop=FALSE]
 active<-any(!closed(jobs[match(fit$task,jobs$task),,drop=FALSE]))
 local<-a[a$slot_class=='existing_local',,drop=FALSE];latest<-local[!duplicated(local$slot_key,fromLast=TRUE),,drop=FALSE]
 latest_decision<-if(is.null(d))NULL else d[match(latest$task,d$task),,drop=FALSE]
 latest_jobs<-jobs[match(latest$task,jobs$task),,drop=FALSE]
 latest_failed<-closed(latest_jobs)&!(latest_jobs$state=='COMPLETED'&latest_jobs$exit_code=='0:0')
 incomplete<-is.null(d)||any(is.na(latest_decision$task)&!latest_failed)
 eligible<-!incomplete&&any(vapply(seq_len(nrow(latest)),function(i){
  isTRUE(latest_decision$decision[i]=='continue')&&sum(local$slot_key==latest$slot_key[i])<3L},logical(1)))
 pb<-if(is.null(batches))NULL else batches[batches$policy==policy,,drop=FALSE]
 # Acceptance cannot skip recording a completed batch: a later external
 # counterexample may invalidate the baseline, but cannot rewrite this outcome.
 if(p$pass&&!is.null(pb)&&nrow(pb)&&p$known_feasible_pass){
  done_b<-pb[which.max(pb$batch),,drop=FALSE];done_prefix<-paste0(policy,':',done_b$batch,':')
  done_g<-a[a$slot_class=='baseline_global'&startsWith(a$slot_key,done_prefix),,drop=FALSE]
  done_polish<-a[a$slot_class=='baseline_polish'&startsWith(a$slot_key,done_prefix),,drop=FALSE]
  if(nrow(done_g)==2L&&nrow(done_polish)==2L&&!is.null(d)&&all(c(done_g$task,done_polish$task)%in%d$task)){
   done_chains<-read_tsv(file.path(study,'tables/baseline_chain_summary.tsv'))
   freeze_completed_outcome(policy,done_b,done_chains[done_chains$policy==policy&done_chains$family=='uniform',,drop=FALSE],p)
  }
 }
 if(p$pass){choice<-'accepted';reason<-'Actual minimum, checked representative and genuine independent witnesses pass'}
 else if(active){choice<-'await_active_parents'}
 else if(incomplete){choice<-'complete_independent_checks'}
 else if(eligible){choice<-'continue_local';reason<-'Whole-window/full-check criterion and remaining existing-root slots justify continuation'}
 else if(is.null(pb)||!nrow(pb)){choice<-'initial_global_batch';batch<-1L;reason<-'Finite available current-root work stopped or exhausted; independent recovery remains unresolved'}
 else{
  b<-pb[which.max(pb$batch),,drop=FALSE];batch<-as.integer(b$batch)
  prefix<-paste0(policy,':',batch,':');g<-a[a$slot_class=='baseline_global'&startsWith(a$slot_key,prefix),,drop=FALSE]
  polish<-a[a$slot_class=='baseline_polish'&startsWith(a$slot_key,prefix),,drop=FALSE]
  checked<-if(is.null(d))character() else d$task
  if(nrow(g)!=2L){choice<-'reconcile_global_batch';reason<-'Both branch receipts must be reconciled without duplication'}
  else if(nrow(polish)!=sum(jobs$state[match(g$task,jobs$task)]=='COMPLETED'&jobs$exit_code[match(g$task,jobs$task)]=='0:0')){choice<-'polish_global_endpoints';reason<-'Each successful original branch has one declared polish allowance'}
  else if(any(!polish$task%in%checked&jobs$state[match(polish$task,jobs$task)]=='COMPLETED'&jobs$exit_code[match(polish$task,jobs$task)]=='0:0')){choice<-'complete_independent_checks'}
  else{
   chains<-read_tsv(file.path(study,'tables/baseline_chain_summary.tsv'));u<-chains[chains$policy==policy&chains$family=='uniform',,drop=FALSE]
   fully_complete<-nrow(polish)==2L&&all(c(g$task,polish$task)%in%checked)
   # Freeze the outcome at the first complete, target-map-verified batch
   # review. Later repair gains must never be credited to that earlier batch.
   eligible_batch<-FALSE
   if(fully_complete&&p$known_feasible_pass){
    outcome<-freeze_completed_outcome(policy,b,u,p)
    eligible_batch<-isTRUE(outcome$eligible)
   }
   if(batch==1L&&fully_complete&&eligible_batch){choice<-'next_global_batch';batch<-2L;reason<-'Frozen complete first-batch outcome advanced a protected minimum by more than0.002'}
   else{
    repair<-a[a$slot_class=='baseline_repair',,drop=FALSE]
    if(nrow(repair)>=2L){choice<-'resolution_unresolved';reason<-'Finite existing-root, global and repair allowances ended without independent recovery; no extra windows authorized'}
    else if(!nrow(repair)){choice<-'finite_baseline_repair';reason<-'A separately unresolved recovery/check condition remains; first finite repair slot available'}
    else{
     prior<-jobs[jobs$task==tail(repair$task,1L),,drop=FALSE]
     proof<-if(is.null(d))NULL else d[d$task==prior$task,,drop=FALSE]
     if(prior$state!='COMPLETED'||prior$exit_code!='0:0'){choice<-'resolution_unresolved';reason<-'The allocated repair failed; no verified improvement supports another local repair'}
     else if(is.null(proof)||!nrow(proof)){choice<-'complete_independent_checks';reason<-'Await the completed repair endpoint full-check evidence'}
     else if(isTRUE(proof$complete)&&proof$decision=='continue'){choice<-'finite_baseline_repair';reason<-'The completed repair meets the whole-window/full-check improvement rule and one repair slot remains'}
     else{choice<-'resolution_unresolved';reason<-'The completed repair is below the improvement resolution; unused allowance is not another fitting entitlement'}
    }
   }
  }
 }
 result[[policy]]<-data.frame(policy=policy,decision=choice,batch=batch,reason=reason,model_evaluations=0L)
}
write_tsv(do.call(rbind,result),file.path(study,'tables/strategy_review.tsv'));print(do.call(rbind,result),row.names=FALSE)
assert_namespaces()
