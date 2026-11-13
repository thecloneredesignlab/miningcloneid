# Closed-record arithmetic and lineage gates only; never calls the tumor model.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);action<-args[2]
entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'))
source(file.path(entry,'audit_math.R'));fast_audit<-install_validated_audit_math(study,entry)
source(file.path(entry,'metadata_stage.R'))
source(file.path(entry,'resolution_rules.R'))
verify_inputs(study)
jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character');contract<-readRDS(file.path(study,'manifests/contract.rds'));backend<-load_backend(study,FALSE)
release<-if(action=='engineering'&&length(args)>2L)args[3L] else read.delim(file.path(study,'tables/execution_release.tsv'),stringsAsFactors=FALSE)$value[1];evidence<-character()
if(fast_audit)evidence<-c(evidence,'tables/audit_math_validation.tsv','tables/audit_math_validation_cases.tsv')
root_for<-function(task){j<-jobs[jobs$task==task,,drop=FALSE];file.path(study,if(j$mode%in%c('local','global'))'fits' else 'tests',task)}
closed<-function(task){j<-jobs[jobs$task==task,,drop=FALSE];nrow(j)==1L&&j$state=='COMPLETED'&&j$exit_code=='0:0'}
payload_verified<-new.env(hash=TRUE,parent=emptyenv())
payload<-function(task){
 stopifnot(closed(task));p<-root_for(task);fingerprint<-hash_file(file.path(p,'payload_manifest.tsv'))
 if(!exists(task,envir=payload_verified,inherits=FALSE)){verify_payload_R(p);assign(task,fingerprint,envir=payload_verified)}
 else stopifnot(identical(get(task,envir=payload_verified),fingerprint))
 evidence<<-unique(c(evidence,paste0(substring(p,nchar(study)+2L),'/payload_manifest.tsv')));p
}
prediction_proof<-function(task){
 path<-paste0('tables/prediction_audits/',task,'_summary.tsv');z<-read_kv(file.path(study,path))
 stopifnot(z$pass,z$payload_sha256==hash_file(file.path(root_for(task),'payload_manifest.tsv')));evidence<<-unique(c(evidence,path))
}
audit<-function(r){
 stopifnot(r$source_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')),r$contract_sha256==hash_file(file.path(study,'manifests/contract.rds')),all(r$parameters>=contract$lower),all(r$parameters<=contract$upper),r$components$n_burden_obs_total==92L,r$components$n_ploidy_obs_total==6502L,r$components$n_necrosis_obs_total==6L)
 ctx<-contract;ctx$backend<-backend;ctx$cfg$parameter_table<-file.path(study,'source/seed25/parameter_table_input.csv');if(r$policy=='data_only')ctx$cfg$use_soft_prior<-FALSE
 stopifnot(max(abs(independent_components(ctx,r$parameters,r$predictions)-unlist(r$components[c('L_b','L_p','L_n')])))<1e-8)
 invisible(TRUE)
}
candidate<-function(task){p<-payload(task);a<-readRDS(file.path(p,'accepted_candidate.rds'));audit(a$record);stopifnot(abs(a$record$components$L-a$checked$L)<1e-8);a$source_task<-task;a$source_scope<-if(jobs$mode[match(task,jobs$task)]%in%c('local','global'))'fits' else 'tests';a$root_id<-jobs$root_id[match(task,jobs$task)];a$root_family<-jobs$root_family[match(task,jobs$task)];a}
emit<-function(name,pass,extra=list()){
 if(name=='baseline_gate.tsv'&&isTRUE(pass))evidence<-freeze_baseline_evidence(study,evidence)
 ep<-paste0('tables/',sub('[.]tsv$','',name),'_evidence.tsv');write_tsv(data.frame(path=unique(evidence),sha256=vapply(file.path(study,unique(evidence)),hash_file,character(1))),file.path(study,ep))
 write_kv(c(list(pass=pass,release=release,input_manifest_sha256=hash_file(file.path(study,'manifests/input_manifest.tsv')),release_manifest_sha256=hash_file(file.path(study,'releases',release,'manifest.tsv')),evidence_manifest=ep,model_evaluations=0L),extra),file.path(study,'tables',name))
}
if(action=='engineering'){
 parallel_tasks<-jobs$task[jobs$kind=='engineering'&jobs$mode=='engineering_parallel'&jobs$release==release]
 serial_tasks<-jobs$task[jobs$kind=='engineering'&jobs$release=='cont03'&jobs$mode%in%c('engineering_serial','replay')]
 tasks<-c(parallel_tasks,serial_tasks);stopifnot(length(parallel_tasks)==1L,length(serial_tasks)==6L,all(vapply(tasks,closed,logical(1))))
 # Later parallel-only repairs retain the real cont03 local/replay proofs by
 # explicit byte identity of every serial runtime method and branch.
 if(release!='cont03') {
  for(file in c('adapter.R','helpers.R','machine.R','observation.R','reference_checks.R','control.R')) {
   oldpath<-file.path(study,'releases/cont03',file);newpath<-file.path(study,'releases',release,file)
   stopifnot(hash_file(oldpath)==hash_file(newpath));evidence<-c(evidence,paste0('releases/cont03/',file),paste0('releases/',release,'/',file))
  }
  evaluator_section<-function(path)strsplit(paste(readLines(path),collapse='\n'),'run_global <- function',fixed=TRUE)[[1]][1L]
  stopifnot(identical(evaluator_section(file.path(study,'releases/cont03/optimize.R')),evaluator_section(file.path(study,'releases',release,'optimize.R'))))
  evidence<-c(evidence,'releases/cont03/optimize.R',paste0('releases/',release,'/optimize.R'),paste0('releases/',release,'/global.R'))
  section<-function(path){z<-paste(readLines(path),collapse='\n');a<-regexpr(" if(mode=='local')",z,fixed=TRUE)[1];b<-regexpr(" } else if(mode=='engineering_parallel')",z,fixed=TRUE)[1];substr(z,a,b-1L)}
  stopifnot(identical(section(file.path(study,'releases/cont03/run.R')),section(file.path(study,'releases',release,'run.R'))))
  evidence<-c(evidence,'releases/cont03/run.R',paste0('releases/',release,'/run.R'))
 }
 for(task in tasks){p<-payload(task);z<-read_kv(file.path(p,'result.tsv'));stopifnot(z$pass)}
 parallel<-tasks[jobs$mode[match(tasks,jobs$task)]=='engineering_parallel'];p<-root_for(parallel)
 unit<-read_kv(file.path(p,'analytic_fixture.tsv'));stopifnot(unit$pass,unit$maximum_visits>25L,unit$partial_poll_sequence_identity,unit$observation_boundary_identity)
 partial<-read_kv(file.path(p,'result.tsv'));stopifnot(partial$partial_incumbent_retained,partial$resumed_root_calls==211L,partial$partial_window_calls==131L)
 if(!is.null(partial$short_partial_calls)){
  short<-read_kv(file.path(p,'short_partial_first/blocks/block_000002/initial_population_validation.tsv'))
  stopifnot(partial$short_partial_calls==97L,partial$short_resumed_root_calls==177L,short$partial_population_subset,!short$initial_population_complete,!short$exact_population_multiset,short$native_initial_evaluations==17L)
  stopifnot(partial$projected_global_fixed_residual==0,partial$tiny_partial_callbacks==1L,partial$tiny_incumbent_retained)
 }
 windows<-list.dirs(file.path(p,'cases'),recursive=FALSE,full.names=TRUE);stopifnot(length(windows)==8L)
 for(window in windows){verify_payload_R(window);a<-readRDS(file.path(window,'accepted_candidate.rds'));audit(a$record);z<-read_kv(file.path(window,'result.tsv'));stopifnot(z$evaluations==160L,z$blocks>=2L,file.exists(file.path(window,'global_checkpoint.rds')))
  for(block in list.dirs(file.path(window,'blocks'),recursive=FALSE,full.names=TRUE)){stopifnot(read_kv(file.path(block,'initial_population_validation.tsv'))$exact_population_multiset,read_kv(file.path(block,'worker_shutdown.tsv'))$workers_stopped)}
 }
 for(task in tasks[jobs$mode[match(tasks,jobs$task)]=='engineering_serial']){p<-root_for(task);z<-read_kv(file.path(p,'result.tsv'));stopifnot(z$real_checkpoint_sequence_identity)
  for(name in c('split_first','split_second','uninterrupted')){a<-readRDS(file.path(p,name,'accepted_candidate.rds'));audit(a$record);q<-read_kv(file.path(p,name,'result.tsv'));stopifnot(q$fixed_residual==0)}
 }
 for(task in tasks[jobs$mode[match(tasks,jobs$task)]=='replay']){stopifnot(read_kv(file.path(root_for(task),'result.tsv'))$maximum_repeat_gap<1e-11);prediction_proof(task)}
 geometry<-paste0('tables/engineering_population_geometry_',release,'.tsv')
 if(file.exists(file.path(study,geometry))){
  gp<-root_for(parallel);g<-read_tsv(file.path(study,geometry));stopifnot(all(g$pass),all(g$model_evaluations==0L),nrow(g)==length(list.files(gp,pattern='^actual_initial_population[.]rds$',recursive=TRUE)))
  stopifnot(all(vapply(file.path(gp,g$block,'actual_initial_population.rds'),hash_file,character(1))==g$population_sha256));evidence<-c(evidence,geometry)
 }
 if(file.exists(file.path(study,'tables/cache_source_manifest.tsv')))evidence<-c(evidence,'tables/cache_source_manifest.tsv')
 emit('engineering_gate.tsv',TRUE,list(owned_parents=7L,DE_cases=8L,real_checkpoint_cases=4L,source_replays=6L,beyond_25_visits=TRUE));cat('Engineering gate passed.\n')
} else if(action=='baseline'){
 roots<-read_tsv(file.path(study,'manifests/roots.tsv'));rows<-list();points<-list();floors<-setNames(vapply(c('production','data_only'),function(policy)min(roots$expected_objective[roots$policy==policy]),numeric(1)),c('production','data_only'))
 # Every canonically exported improvement is a feasible floor, also under the
 # other policy's unchanged prior. It does not create another recovery lineage.
 files<-character()
 for(task in jobs$task[jobs$state=='COMPLETED'&jobs$exit_code=='0:0']){
  p<-payload(task);files<-c(files,list.files(p,pattern='^improvement_[0-9]+[.]rds$',recursive=TRUE,full.names=TRUE))
 }
 # An obsolete/failed engineering parent can contain fully closed native
 # children. Infrastructure failure does not erase their verified feasible
 # points. Unfinished children remain archived and require a fresh replay.
 reviewed<-read_tsv(file.path(study,'tables/reviewed_failures.tsv'))$task
 for(task in jobs$task[jobs$task%in%reviewed&jobs$mode=='engineering_parallel'&jobs$state!='COMPLETED']){
  p<-root_for(task);cases<-if(dir.exists(file.path(p,'cases')))list.dirs(file.path(p,'cases'),recursive=FALSE,full.names=TRUE) else character()
  candidates<-c(cases,file.path(p,c('partial_first','partial_resumed','short_partial_first','short_partial_resumed','projected_global','tiny_population')))
  for(case in cases)if(dir.exists(file.path(case,'blocks')))candidates<-c(candidates,list.dirs(file.path(case,'blocks'),recursive=FALSE,full.names=TRUE))
  for(child in candidates[file.exists(file.path(candidates,'COMPLETE'))]){
   verify_payload_R(child);evidence<-c(evidence,paste0(substring(child,nchar(study)+2L),'/payload_manifest.tsv'))
   files<-c(files,list.files(child,pattern='^improvement_[0-9]+[.]rds$',recursive=TRUE,full.names=TRUE))
  }
 }
 files<-unique(files)
 cfg<-contract$cfg;cfg$parameter_table<-file.path(study,'source/seed25/parameter_table_input.csv');cfg$use_soft_prior<-TRUE
 # Cache only results of the same independent arithmetic audit over exactly
 # hashed immutable records and frozen inputs. Every miss still receives the
 # complete original audit; a cache hit never changes the feasible-point pool.
 cache_path<-file.path(study,'tables/feasible_floor_cache.rds')
 cache_key<-list(input=hash_file(file.path(study,'manifests/input_manifest.tsv')),contract=hash_file(file.path(study,'manifests/contract.rds')),adapter=hash_file(file.path(entry,'adapter.R')),
  audit=digest::digest(paste(deparse(body(audit)),collapse='\n'),serialize=FALSE,algo='sha256'),component_method=hash_file(file.path(entry,'audit_math.R')),fast_audit=fast_audit)
 cache<-new.env(hash=TRUE,parent=emptyenv());hits<-misses<-0L
 if(file.exists(cache_path)){saved<-readRDS(cache_path);if(identical(saved$key,cache_key))list2env(saved$records,envir=cache)}
 file_values<-list()
 for(path in files){
  key<-hash_file(path)
  if(exists(key,envir=cache,inherits=FALSE)){value<-get(key,envir=cache);hits<-hits+1L}
  else{
   r<-readRDS(path);audit(r)
   value<-c(data_only=r$components$L_data,production=r$components$L_data+cfg$lambda_prior*backend$compute_soft_prior_penalty(r$parameters,cfg)$L_prior_raw)
   assign(key,value,envir=cache);misses<-misses+1L
  }
  stopifnot(identical(sort(names(value)),c('data_only','production')),all(is.finite(value)))
  file_values[[path]]<-value
  floors[names(value)]<-pmin(floors[names(value)],value)
 }
 atomic_rds(list(key=cache_key,records=as.list(cache)),cache_path);evidence<-c(evidence,'tables/feasible_floor_cache.rds')
 write_kv(list(cache_hits=hits,independent_audit_misses=misses,verified_records=length(files),input_manifest_sha256=cache_key$input,audit_function_sha256=cache_key$audit,cache_sha256=hash_file(cache_path),model_evaluations=0L),file.path(study,'tables/feasible_cache_usage.tsv'))
 evidence<-c(evidence,'tables/feasible_cache_usage.tsv')
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
