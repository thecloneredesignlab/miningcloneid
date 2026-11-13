# Evidence-driven baseline strategy decision; no model calls, no resource ceiling.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character');tables<-file.path(study,'tables')
terminal<-c('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')
closed<-function(j)sub('[+ ].*$','',j$state)%in%terminal
successful<-function(j)all(j$state=='COMPLETED'&j$exit_code=='0:0')
optional<-function(name)if(file.exists(file.path(tables,name)))read_tsv(file.path(tables,name)) else NULL
policies<-optional('baseline_policy_gates.tsv');chains<-optional('baseline_chain_summary.tsv');batches<-optional('search_batches.tsv');result<-list()
reviewed<-optional('reviewed_failures.tsv');reviewed_tasks<-if(is.null(reviewed))character() else reviewed$task
for(policy in c('production','data_only')){
 x<-jobs[jobs$kind=='baseline'&jobs$policy==policy,,drop=FALSE];p<-if(is.null(policies))NULL else policies[policies$policy==policy,,drop=FALSE]
 decision<-'await_checked_local';reason<-'No complete independently checked local cohort';batch<-0L;uniform_gain<-feasible_gain<-NA_real_
 checked_policy<-if(is.null(chains))NULL else chains[chains$policy==policy,,drop=FALSE]
 previous_policy<-if(is.null(batches))NULL else batches[batches$policy==policy,,drop=FALSE]
 initial_ready<-!is.null(p)&&nrow(p)&&!p$pass&&!is.null(checked_policy)&&any(checked_policy$termination=='poll_stationary')&&any(checked_policy$family=='uniform')&&(is.null(previous_policy)||!nrow(previous_policy))&&!any(closed(x)&(x$state!='COMPLETED'|x$exit_code!='0:0')&!x$task%in%reviewed_tasks)
 if(!is.null(p)&&nrow(p)&&p$pass){decision<-'accepted';reason<-'Policy scientific recovery gates passed; wait for all baseline parents before profile release'}
 else if(initial_ready){decision<-'initial_global_batch';reason<-'A fully checked local root is stationary and a verified independent uniform chain still fails recovery; other ready local windows may continue concurrently';batch<-1L}
 else if(nrow(x)&&any(!closed(x))){decision<-'await_active_parents';reason<-'Current fitting/check windows still active; no stagnation inference'}
 else if(nrow(x)&&any((x$state!='COMPLETED'|x$exit_code!='0:0')&!x$task%in%reviewed_tasks)){decision<-'review_execution_failure';reason<-'A baseline allocation failed; repair/review before numerical strategy decisions'}
 else if(!is.null(chains)){
  locals<-x[x$mode=='local',,drop=FALSE];component<-sub('_w[0-9]+$','',locals$task);latest<-locals[!duplicated(component,fromLast=TRUE),,drop=FALSE]
  checked<-chains[match(latest$task,chains$task),,drop=FALSE]
  if(any(is.na(checked$task))){decision<-'complete_independent_checks';reason<-'Every latest local endpoint needs fresh and original-coordinate checks'}
  else if(any(checked$termination%in%c('time_limit','call_limit'))){decision<-'continue_local';reason<-'Incomplete resource-window polls must resume their exact state or a material own-root probe'}
  else if(nrow(checked)&&all(checked$termination=='poll_stationary')){
   previous<-if(is.null(batches))NULL else batches[batches$policy==policy,,drop=FALSE]
   if(is.null(previous)||!nrow(previous)){decision<-'initial_global_batch';reason<-'All current local components are checked/poll-stationary and independent recovery still fails';batch<-1L}
   else {
    last<-previous[which.max(previous$batch),,drop=FALSE];batch<-as.integer(last$batch);prefix<-paste0('g_',policy,'_c',sprintf('%04d',batch),'_')
    globals<-x[x$mode=='global'&startsWith(x$task,prefix),,drop=FALSE];components<-sub('_w[0-9]+$','',globals$task);latest_global<-globals[!duplicated(components,fromLast=TRUE),,drop=FALSE]
    expected<-paste0(prefix,c('hybrid','uniform'))
    if(nrow(latest_global)!=2L||!setequal(sub('_w[0-9]+$','',latest_global$task),expected)){decision<-'reconcile_global_batch';reason<-'Both original exploration branches must be reconciled, not duplicated'}
    else if(any(!paste0('b_',latest_global$task,'_polish_w0001')%in%locals$task)){decision<-'polish_global_endpoints';reason<-'A closed latest population endpoint still needs its own local branch and checks'}
    else {
     gains<-vapply(seq_len(nrow(latest_global)),function(i){task<-latest_global$task[i];z<-read_kv(file.path(study,'fits',task,'result.tsv'));if(is.na(z$window_gain)){blocks<-read_tsv(file.path(study,'fits',task,'blocks.tsv'));blocks$objective[1]-z$objective}else z$window_gain},numeric(1))
     uniform_now<-min(chains$verified_best[chains$policy==policy&chains$family=='uniform']);feasible_now<-min(p$best_objective,p$known_feasible_bound)
     uniform_gain<-last$before_uniform-uniform_now;feasible_gain<-last$before_feasible-feasible_now
     if(any(gains>1e-10)){decision<-'continue_global';reason<-'At least one latest native population window still gives effective improvement; retain its own population and all partial calls'}
     else if(max(uniform_gain,feasible_gain)>1e-10){decision<-'next_global_batch';reason<-'The completed batch advanced a protected minimum; a further deficit-driven batch remains justified';batch<-batch+1L}
     else {decision<-'strategy_stagnation';reason<-'Complete two-branch batch plus all local polishes gave no effective protected-minimum gain; all relevant local polls are stationary and recovery still fails'}
    }
   }
  }
 }
 result[[policy]]<-data.frame(policy=policy,decision=decision,batch=batch,uniform_gain=uniform_gain,feasible_gain=feasible_gain,reason=reason,model_evaluations=0L)
}
write_tsv(do.call(rbind,result),file.path(tables,'strategy_review.tsv'));print(do.call(rbind,result),row.names=FALSE);assert_namespaces()
