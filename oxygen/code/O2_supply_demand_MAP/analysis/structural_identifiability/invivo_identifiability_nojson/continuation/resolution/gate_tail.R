 # The actual minimum, stopping representative and repair anchor are distinct.
 tables<-file.path(study,'tables')
 minimums<-list();replays<-list();ledger<-list()
 for(policy in c('production','data_only')){
  order_files<-files[order(vapply(file_values,function(v)unname(v[policy]),numeric(1)))]
  chosen<-NULL;lowest_source<-NULL
  for(path in order_files){
   r<-readRDS(path)
   if(is.null(lowest_source))lowest_source<-list(record=r,path=path,objective=unname(file_values[[path]][policy]))
   if(r$policy==policy&&!length(r$fixed)){
    owner<-sub('/.*$','',r$provenance$task)
    chosen<-list(record=r,checked=r$components,policy=policy,root_id=r$provenance$root_id,
      root_family=r$provenance$root_family,source_task=owner,
      source_scope=if(jobs$mode[match(owner,jobs$task)]%in%c('local','global'))'fits' else 'tests',
      source_record=substring(path,nchar(study)+2L));break
   }
  }
  stopifnot(!is.null(chosen));minimums[[policy]]<-chosen
  atomic_rds(chosen,file.path(tables,paste0('baseline_',policy,'.rds')))
  evidence<-c(evidence,paste0('tables/baseline_',policy,'.rds'))
  ledger[[policy]]<-data.frame(policy=policy,candidate=chosen$checked$L,feasible_floor=unname(floors[policy]),
    candidate_root=chosen$root_id,candidate_source=chosen$source_record,policy_version='resolution01')
  if(chosen$checked$L>floors[policy]+1e-8){
   r<-lowest_source$record;rid<-paste0('informed_',policy,'_',substr(hash_file(lowest_source$path),1,12))
   a<-list(record=r,policy=policy,root_id=rid,root_family='informed',source_task=r$provenance$task,
       expected_target_objective=lowest_source$objective,original_policy=r$policy)
   name<-paste0('resolution_target_replay_',policy);atomic_rds(a,file.path(tables,paste0(name,'.rds')))
   replays[[policy]]<-data.frame(policy=policy,selection=name,root_id=rid,expected=lowest_source$objective,
        record_sha256=hash_file(lowest_source$path))
  }
 }
 if(length(replays))write_tsv(do.call(rbind,replays),file.path(tables,'resolution_target_replays.tsv'))
 else write_tsv(data.frame(policy=character(),selection=character(),root_id=character(),expected=numeric(),record_sha256=character()),file.path(tables,'resolution_target_replays.tsv'))
 write_tsv(do.call(rbind,ledger),file.path(tables,'resolution_candidates.tsv'))
 evidence<-c(evidence,'tables/resolution_target_replays.tsv',if(length(replays))paste0('tables/',vapply(replays,function(r)r$selection,character(1)),'.rds'))
 checks<-jobs[jobs$mode=='resolution_check'&jobs$state=='COMPLETED'&jobs$exit_code=='0:0',,drop=FALSE]
 rows<-list();representatives<-list()
 for(i in seq_len(nrow(checks))){
  j<-checks[i,];spec<-read_kv(file.path(study,'dispatch',j$task,'spec.tsv'))
  task<-spec$checked_fit;if(is.null(task)||!closed(task))next
  fit<-jobs[jobs$task==task,,drop=FALSE];if(fit$kind!='baseline')next
  d<-payload(j$task);a<-candidate(task);b<-candidate(j$task);z<-read_kv(file.path(root_for(task),'result.tsv'))
  dr<-read_kv(file.path(d,'result.tsv'));center<-readRDS(file.path(d,'checked_center.rds'));prediction_proof(j$task)
  stopifnot(dr$pass,dr$complete,dr$scope=='fixed_center_all_four_scales_axis_coupled_observation',
     abs(dr$difference)<1e-8,dr$fixed_residual==0,center$root_id==fit$root_id,dr$free_coordinates==20L)
  stopifnot(a$record$provenance$root_id==fit$root_id,b$record$provenance$root_id==fit$root_id,
    center$record$provenance$root_id==fit$root_id,a$record$provenance$root_family==fit$root_family,
    b$record$provenance$root_family==fit$root_family,center$record$provenance$root_family==fit$root_family)
  best<-min(a$checked$L,b$checked$L);gap<-center$checked$L-best
  valid<-gap<=1e-4&&dr$directional_improvement<=.002
  rows[[j$task]]<-data.frame(task=task,check_task=j$task,policy=j$policy,root_id=j$root_id,
    family=fit$root_family,objective=a$checked$L,verified_best=best,representative_objective=center$checked$L,
    representative_gap=gap,representative_pass=valid,probe_gain=dr$directional_improvement,
    termination=z$termination,window_gain=z$window_gain,method_calls=if(is.null(z$method_calls))NA else z$method_calls)
  representatives[[j$task]]<-center
 }
 tab<-if(length(rows))do.call(rbind,rows) else data.frame()
 policies<-list();decisions<-list();anchors_path<-file.path(tables,'resolution_anchors.tsv')
 anchors<-if(file.exists(anchors_path))read_tsv(anchors_path) else data.frame(policy=character(),objective=numeric(),source=character(),status=character())
 for(policy in c('production','data_only')){
  x<-if(nrow(tab))tab[tab$policy==policy,,drop=FALSE] else tab;floor<-unname(floors[policy])
  minrec<-minimums[[policy]];known<-minrec$checked$L<=floor+1e-8
  rep<-if(nrow(x))x[x$representative_pass&x$representative_objective<=floor+1e-4,,drop=FALSE] else x
  good<-if(nrow(x))x[x$representative_pass&x$verified_best<=floor+.002,,drop=FALSE] else x
  recovery<-nrow(good)>0L&&any(vapply(unique(good$root_id[good$family=='archive']),function(a)any(good$family=='uniform'&good$root_id!=a),logical(1)))
  uniform_gap<-if(nrow(x)&&any(x$family=='uniform'))min(x$verified_best[x$family=='uniform'])-floor else Inf
  uniform_rep_gap<-if(nrow(x)&&any(x$family=='uniform'))min(x$representative_objective[x$family=='uniform'])-floor else Inf
  rep_pass<-nrow(rep)>0L;gain<-if(rep_pass)min(rep$probe_gain) else Inf
  if(rep_pass){
   bestrep<-rep[which.min(rep$representative_objective),,drop=FALSE];rep_record<-representatives[[bestrep$check_task]]
   atomic_rds(rep_record,file.path(tables,paste0('baseline_stopping_rep_',policy,'.rds')))
   evidence<-c(evidence,paste0('tables/baseline_stopping_rep_',policy,'.rds'))
   if(!policy%in%anchors$policy)anchors<-rbind(anchors,data.frame(policy=policy,objective=bestrep$representative_objective,source=bestrep$check_task,status='initial_completed_reference'))
  }
  cumulative<-if(policy%in%anchors$policy)anchors$objective[match(policy,anchors$policy)]-floor else NA_real_
  policies[[policy]]<-data.frame(policy=policy,pass=known&&recovery&&rep_pass,
    family_recovery=recovery,directional_pass=rep_pass,known_feasible_pass=known,
    best_objective=minrec$checked$L,known_feasible_bound=floor,uniform_gap=uniform_gap,uniform_representative_gap=uniform_rep_gap,
    directional_gain=gain,representative_available=rep_pass,cumulative_gain=cumulative,
    material_repair_trigger=!is.na(cumulative)&&isTRUE(resolution_compare(cumulative,.002)>0L))
 }
 # Preserve the original completed anchor vector even after later small gains.
 if(!'selection'%in%names(anchors))anchors$selection<-rep('',nrow(anchors))
 if(!'selection_sha256'%in%names(anchors))anchors$selection_sha256<-rep('',nrow(anchors))
 for(i in seq_len(nrow(anchors))){
  name<-paste0('resolution_anchor_',anchors$policy[i],'.rds');target<-file.path(tables,name)
  source<-file.path(payload(anchors$source[i]),'checked_center.rds');point<-readRDS(source)
  stopifnot(abs(point$checked$L-anchors$objective[i])<1e-8,point$policy==anchors$policy[i])
  if(!file.exists(target))atomic_rds(point,target)
  stopifnot(abs(readRDS(target)$checked$L-anchors$objective[i])<1e-8)
  anchors$selection[i]<-name;anchors$selection_sha256[i]<-hash_file(target)
  evidence<-c(evidence,paste0('tables/',name))
 }
 write_tsv(anchors,anchors_path);write_tsv(do.call(rbind,policies),file.path(tables,'baseline_policy_gates.tsv'))
 if(nrow(tab)){
  tab$gap<-tab$verified_best-vapply(tab$policy,function(p)unname(floors[p]),numeric(1))
  write_tsv(tab,file.path(tables,'baseline_chain_summary.tsv'))
  for(i in seq_len(nrow(tab))){
   x<-tab[i,];pass<-policies[[x$policy]]$pass
   decision<-resolution_local_rule(x$window_gain,x$probe_gain,TRUE,pass,x$representative_pass)
   decisions[[x$check_task]]<-data.frame(task=x$task,check_task=x$check_task,policy=x$policy,
      root_id=x$root_id,window_gain=x$window_gain,full_check_gain=x$probe_gain,
      representative_gap=x$representative_gap,complete=TRUE,decision=decision,
      selector=if(isTRUE(resolution_compare(x$probe_gain,.002)>0L))'task' else 'checkpoint',
      continuation_source=if(isTRUE(resolution_compare(x$probe_gain,.002)>0L))x$check_task else x$task,policy_version='resolution01')
  }
  write_tsv(do.call(rbind,decisions),file.path(tables,'resolution_local_decisions.tsv'))
 }
 # Root selection retains genuine ancestry and never injects another root.
 for(root in unique(jobs$root_id[jobs$root_family%in%c('archive','uniform','informed')])){
  candidates<-jobs$task[jobs$root_id==root&jobs$state=='COMPLETED'&jobs$exit_code=='0:0'&jobs$kind!='engineering'&jobs$mode%in%c('local','global','resolution_check','directions','reference_reconstruct')]
  if(!length(candidates))next
  choices<-lapply(candidates,candidate);best<-choices[[which.min(vapply(choices,function(a)a$checked$L,numeric(1)))]]
  name<-paste0('root_',root,'_selected');atomic_rds(best,file.path(tables,paste0(name,'.rds')))
  write_kv(list(root_id=root,policy=best$policy,root_family=best$root_family,selection_sha256=hash_file(file.path(tables,paste0(name,'.rds')))),file.path(tables,paste0(name,'.tsv')))
  evidence<-c(evidence,paste0('tables/',name,'.rds'),paste0('tables/',name,'.tsv'))
 }
 evidence<-unique(c(evidence,'tables/resolution_candidates.tsv','tables/baseline_policy_gates.tsv',
   'tables/resolution_anchors.tsv',if(nrow(tab))c('tables/resolution_local_decisions.tsv','tables/baseline_chain_summary.tsv')))
 review_evidence<-'tables/resolution_review_evidence.tsv'
 write_tsv(data.frame(path=evidence,sha256=vapply(file.path(study,evidence),hash_file,character(1))),file.path(study,review_evidence))
 write_kv(list(pass=TRUE,policy_version='resolution01',input_manifest_sha256=hash_file(file.path(study,'manifests/input_manifest.tsv')),
   evidence_manifest=review_evidence,model_evaluations=0L),file.path(tables,'resolution_review_gate.tsv'))
 pending<-any(jobs$kind=='baseline'&!sub('[+ ].*$','',jobs$state)%in%c('COMPLETED','FAILED','TIMEOUT','CANCELLED','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'))
 emit('baseline_gate.tsv',all(vapply(policies,function(p)p$pass,logical(1)))&&!pending,list(policy_version='resolution01'))
 print(do.call(rbind,policies),row.names=FALSE)
} else stop('Unsupported metadata gate')
assert_namespaces()
