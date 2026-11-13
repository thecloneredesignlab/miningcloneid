# Closed-payload acceptance only; no model is loaded or evaluated.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L)
study<-normalizePath(args[1]);action<-args[2]
entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'))
specs<-readRDS(file.path(study,'manifests/task_specs.rds'));contract<-readRDS(file.path(study,'manifests/contract.rds'))
protocol<-readRDS(file.path(study,'manifests/protocol.rds'));jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character')
backend<-load_backend(study,FALSE);evidence<-character()
release<-if(file.exists(file.path(study,'tables/execution_release.tsv')))read_kv(file.path(study,'tables/execution_release.tsv'))$release else 'repair01'
if(action=='engineering')release<-'repair01'
if(release=='repair02') {
 bridge<-read_tsv(file.path(study,'tables/release_bridge_repair02.tsv'))
 stopifnot(all(vapply(file.path(study,bridge$path),hash_file,character(1))==bridge$sha256),
  isTRUE(read_kv(file.path(study,'tables/engineering_gate_repair02.tsv'))$pass))
}
root_for<-function(task)file.path(study,if(specs[[task]]$mode%in%c('global','local'))'fits' else 'tests',task)
closed<-function(task) {
 j<-jobs[jobs$task==task,,drop=FALSE]
 nrow(j)==1L && j$state=='COMPLETED' && j$exit_code=='0:0'
}
payload<-function(task) {
 stopifnot(closed(task));root<-root_for(task);verify_payload_R(root)
 evidence<<-unique(c(evidence,paste0(substring(root,nchar(study)+2L),'/payload_manifest.tsv')));root
}
audit<-function(record) {
 stopifnot(record$source_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')),
  record$contract_sha256==hash_file(file.path(study,'manifests/contract.rds')),
  all(record$parameters>=contract$lower),all(record$parameters<=contract$upper))
 ctx<-contract;ctx$backend<-backend;ctx$cfg$parameter_table<-file.path(study,'source/seed25/parameter_table_input.csv')
 if(record$policy=='data_only')ctx$cfg$use_soft_prior<-FALSE
 ind<-independent_components(ctx,record$parameters,record$predictions)
 stopifnot(max(abs(ind-unlist(record$components[names(ind)])))<1e-8)
 prior<-backend$compute_soft_prior_penalty(record$parameters,ctx$cfg)
 stopifnot(abs(prior$L_prior_raw-record$components$L_prior_raw)<1e-8)
 invisible(TRUE)
}
candidate<-function(task) {
 root<-payload(task);a<-readRDS(file.path(root,'accepted_candidate.rds'));audit(a$record)
 stopifnot(abs(a$checked$L-a$record$components$L)<1e-8)
 actual<-readRDS(file.path(root,'dispatch_spec.rds'))
 if(actual$root_family!='profile')stopifnot(identical(a$record$provenance$root_id,actual$root_id))
 if(actual$mode=='local' && actual$root_family!='profile') {
  start<-readRDS(file.path(root,'actual_local_start.rds'))$resolution
  if(actual$selector=='import') {
   im<-readRDS(file.path(study,'reference/imports',paste0(actual$import_id,'.rds')))
   stopifnot(im$root_id==actual$root_id,start$source_sha256==hash_file(file.path(study,'reference/imports',paste0(actual$import_id,'.rds'))))
  } else {
   stopifnot(!is.null(start$record),identical(start$record$provenance$root_id,actual$root_id))
   if(actual$selector=='task')payload(actual$source_task)
  }
 }
 list(record=a$record,policy=a$policy,root_id=actual$root_id,root_family=actual$root_family,
  source_task=task,source_scope=if(specs[[task]]$mode%in%c('global','local'))'fits' else 'tests',payload_sha256=hash_file(file.path(root,'payload_manifest.tsv')))
}
emit_gate<-function(name,pass,extra=list()) {
 evidence_path<-paste0('tables/',sub('[.]tsv$','',name),'_evidence.tsv')
 write_tsv(data.frame(path=evidence,sha256=vapply(file.path(study,evidence),hash_file,character(1))),file.path(study,evidence_path))
 write_kv(c(list(pass=pass,release=release,input_manifest_sha256=hash_file(file.path(study,'manifests/input_manifest.tsv')),
  release_manifest_sha256=hash_file(file.path(study,'releases',release,'manifest.tsv')),evidence_manifest=evidence_path,model_evaluations=0L),extra),file.path(study,'tables',name))
}
if(action=='engineering') {
 parents<-names(specs)[vapply(specs,function(s)s$stage=='engineering',logical(1))]
 stopifnot(length(parents)==9L,all(vapply(parents,closed,logical(1))))
 for(task in parents){root<-payload(task);stopifnot(isTRUE(read_kv(file.path(root,'result.tsv'))$pass))}
 unit<-file.path(study,'tests/engineering_parallel_unit');verify_payload_R(unit)
 evidence<-c(evidence,'tests/engineering_parallel_unit/payload_manifest.tsv')
 extra<-file.path(study,'tests/arithmetic_extra')
 if(dir.exists(extra)) {verify_payload_R(extra);stopifnot(isTRUE(read_kv(file.path(extra,'result.tsv'))$pass));evidence<-c(evidence,'tests/arithmetic_extra/payload_manifest.tsv')}
 for(policy in c('production','data_only'))for(init in c('hybrid','uniform'))for(fixed in c(FALSE,TRUE)) {
  name<-paste0('engineering_parallel_',policy,'_',init,'_',if(fixed)'fixed' else 'free');root<-file.path(study,'tests',name);verify_payload_R(root)
  evidence<-c(evidence,paste0('tests/',name,'/payload_manifest.tsv'))
  a<-readRDS(file.path(root,'accepted_candidate.rds'));audit(a$record)
  stopifnot(!a$budget_limited,isTRUE(read_kv(file.path(root,'initial_population_validation.tsv'))$exact_population_multiset),
   isTRUE(read_kv(file.path(root,'worker_shutdown.tsv'))$workers_stopped),length(a$record$free_names)==if(fixed)19L else 20L)
 }
 for(policy in c('production','data_only'))for(fixed in c(FALSE,TRUE)) {
  task<-paste0('engineering_',policy,'_',if(fixed)'fixed' else 'free');a<-candidate(task)
  result<-read_kv(file.path(root_for(task),'result.tsv'))
  stopifnot(!result$budget_limited,!result$call_limited,result$free_coordinates==if(fixed)19L else 20L,result$dynamic_coordinates==18L)
  pp<-readRDS(file.path(root_for(task),'observation_proposals.rds'))
  stopifnot(all(vapply(pp,function(z)abs(z$difference)<1e-8,logical(1))))
  if(!fixed)stopifnot(length(pp)>0L)
  if(fixed && !length(pp))stopifnot(abs(a$record$parameters[['log10_sigma_burden']]-contract$upper[['log10_sigma_burden']])<1e-12)
 }
 emit_gate('engineering_gate.tsv',TRUE,list(parent_allocations=9L,actual_de_cases=8L,pattern_cases=4L,import_cases=4L))
 cat('Repair engineering accepted: nine allocations, eight DE cases, four pattern cases and four provenance/observation replays.\n')
} else if(action=='baseline') {
 # Every engineering point is a real feasible point, including fixed-density
 # fixtures and canonical checks of observation proposals. Keep their floor.
 known<-list()
 for(task in names(specs)[vapply(specs,function(s)s$stage=='engineering' && s$mode%in%c('reference','engineering_serial'),logical(1))]) {
  root<-payload(task);a<-readRDS(file.path(root,'accepted_candidate.rds'));audit(a$record)
  known[[length(known)+1L]]<-data.frame(policy=specs[[task]]$policy,source=task,objective=a$checked$L)
  obs<-file.path(root,'observation_validation.rds')
  if(file.exists(obs))known[[length(known)+1L]]<-data.frame(policy=specs[[task]]$policy,source=paste0(task,':observation'),objective=readRDS(obs)$checked$components$L)
 }
 for(policy in c('production','data_only'))for(init in c('hybrid','uniform'))for(fixed in c(FALSE,TRUE)) {
  task<-paste0('engineering_parallel_',policy,'_',init,'_',if(fixed)'fixed' else 'free');root<-file.path(study,'tests',task)
  verify_payload_R(root);evidence<-c(evidence,paste0('tests/',task,'/payload_manifest.tsv'))
  a<-readRDS(file.path(root,'accepted_candidate.rds'));known[[length(known)+1L]]<-data.frame(policy=policy,source=task,objective=a$checked$L)
 }
 known<-do.call(rbind,known);write_tsv(known,file.path(study,'tables/engineering_feasible_points.tsv'))
 # Cross-policy feasible values are arithmetic comparisons of saved real
 # trajectories, never independent fits or transplanted root endpoints.
 cross_floor<-setNames(rep(Inf,2),c('production','data_only'));cross_rows<-list()
 audit_files<-list.files(file.path(study,'tables/record_audits'),pattern='_records[.]tsv$',full.names=TRUE)
 if(length(audit_files)) {
  ix<-do.call(rbind,lapply(audit_files,read_tsv));cfg<-contract$cfg;cfg$use_soft_prior<-TRUE
  cfg$parameter_table<-file.path(study,'source/seed25/parameter_table_input.csv')
  for(task in unique(ix$task)) {
   root<-payload(task);summary<-file.path(study,'tables/record_audits',paste0(task,'_summary.tsv'))
   proof<-read_kv(summary);stopifnot(isTRUE(proof$pass),proof$payload_sha256==hash_file(file.path(root,'payload_manifest.tsv')))
   evidence<-c(evidence,paste0('tables/record_audits/',task,'_summary.tsv'),paste0('tables/record_audits/',task,'_records.tsv'))
  }
  data_idx<-which.min(ix$data);cross_floor['data_only']<-ix$data[data_idx]
  for(k in seq_len(nrow(ix))) {
   row<-ix[k,];sp<-specs[[row$task]]
   if(sp$policy=='production')value<-row$objective
   else {
    if(row$data>=cross_floor['production'])next
    record<-readRDS(file.path(root_for(row$task),'evaluators',row$evaluator,row$record))
    value<-record$components$L_data+cfg$lambda_prior*backend$compute_soft_prior_penalty(record$parameters,cfg)$L_prior_raw
   }
   if(value<cross_floor['production']) {
    cross_floor['production']<-value
    cross_rows$production<-data.frame(target_policy='production',original_policy=sp$policy,source_task=row$task,evaluator=row$evaluator,record=row$record,objective=value)
   }
  }
  r<-ix[data_idx,];cross_rows$data_only<-data.frame(target_policy='data_only',original_policy=specs[[r$task]]$policy,source_task=r$task,evaluator=r$evaluator,record=r$record,objective=cross_floor['data_only'])
  write_tsv(do.call(rbind,cross_rows),file.path(study,'tables/cross_policy_feasible_floors.tsv'))
 }
 counters<-list()
 if(file.exists(file.path(study,'tables/profile_counterexamples.rds'))) {
  counters<-readRDS(file.path(study,'tables/profile_counterexamples.rds'))
  for(point in counters){payload(point$source_task);audit(point$record)}
 }
 summaries<-list();policy_rows<-list();allpass<-TRUE
 for(policy in c('production','data_only')) {
  endpoints<-list();probes<-list();rows<-list()
  tasks<-names(specs)[vapply(specs,function(s)s$mode=='local' && s$policy==policy && s$role%in%c('incumbent_local','uniform_local','reserve_local'),logical(1))]
  for(task in tasks) {
   fresh<-sub('_local$','_fresh',task);dir<-sub('_local$','_directions',task)
   if(!all(vapply(c(task,fresh,dir),closed,logical(1))))next
   a<-candidate(task);fr<-payload(fresh);pr<-candidate(dir);d<-read_kv(file.path(root_for(dir),'result.tsv'))
   f<-read_kv(file.path(fr,'result.tsv'));stopifnot(isTRUE(f$pass),abs(f$difference)<1e-8,isTRUE(d$pass))
   endpoints[[task]]<-a;probes[[task]]<-pr
   result<-read_kv(file.path(root_for(task),'result.tsv'))
   rows[[task]]<-data.frame(policy=policy,task=task,root_id=a$root_id,family=a$root_family,
    objective=a$record$components$L,probe_best=pr$record$components$L,method=result$method,
    native_code=result$native_convergence,budget_limited=result$budget_limited,independent=specs[[task]]$role!='reserve_local')
  }
  cheap<-paste0('r_',policy,'_',c('anchor','uniform_old'),'_local')
  if(!all(cheap%in%names(rows))){cat('WAIT_BASELINE',policy,'cheap endpoints not complete\n');next}
  tab<-do.call(rbind,rows);tab$verified_best<-pmin(tab$objective,tab$probe_best)
  idx<-which.min(tab$verified_best);best<-tab$verified_best[idx];tab$gap<-tab$verified_best-best
  recovered<-tab$gap<=protocol$baseline_recovery
  families<-unique(tab$family[recovered]);root_pass<-all(c('archive','uniform')%in%families) && length(unique(tab$root_id[recovered]))>=2L
  probe_best<-min(tab$probe_best);local_best<-min(tab$objective)
  direction_gain<-local_best-probe_best;direction_pass<-direction_gain<=protocol$direction_improvement
  import_table<-read_tsv(file.path(study,'reference/imports.tsv'));floor<-min(import_table$expected_objective[import_table$target_policy==policy],known$objective[known$policy==policy],cross_floor[policy])
  if(!is.null(counters[[policy]]))floor<-min(floor,counters[[policy]]$record$components$L)
  floor_pass<-best<=floor+protocol$score_tolerance
  pass<-root_pass && direction_pass && floor_pass;allpass<-allpass && pass
  selected<-if(tab$probe_best[idx]<=tab$objective[idx])probes[[tab$task[idx]]] else endpoints[[tab$task[idx]]]
  selected$diagnostic_center_task<-tab$task[which.min(tab$objective)]
  selected$diagnostic_center_objective<-local_best
  selected$directional_improvement_from_optimized_endpoint<-direction_gain
  selected$selection_kind<-'lowest_canonical_verified_chain_point'
  atomic_rds(selected,file.path(study,'tables',paste0('baseline_',policy,'.rds')))
  evidence<-c(evidence,paste0('tables/baseline_',policy,'.rds'))
  atomic_rds(probes[[tab$task[which.min(tab$probe_best)]]],file.path(study,'tables',paste0('baseline_',policy,'_probe_best.rds')))
  for(root in unique(tab$root_id)) {
   take<-which(tab$root_id==root);k<-take[which.min(tab$probe_best[take])];point<-probes[[tab$task[k]]]
   name<-paste0('root_',root,'_selected');atomic_rds(point,file.path(study,'tables',paste0(name,'.rds')))
   write_kv(list(policy=policy,root_id=root,root_family=point$root_family,selection_sha256=hash_file(file.path(study,'tables',paste0(name,'.rds')))),file.path(study,'tables',paste0(name,'.tsv')))
  }
  summaries[[policy]]<-tab;policy_rows[[policy]]<-data.frame(policy=policy,pass=pass,family_recovery=root_pass,directional_pass=direction_pass,
   known_feasible_pass=floor_pass,known_feasible_bound=floor,known_feasible_gap=best-floor,best_objective=best,best_probe=probe_best,
   uniform_gap=min(tab$gap[tab$family=='uniform']),best_optimized_endpoint=local_best,directional_improvement=direction_gain,best_task=selected$source_task)
  emit_gate(paste0('extended_',policy,'_gate.tsv'),!pass,list(reason=if(pass)'Recovery already accepted; additional globals disabled' else 'Cheap roots closed; recovery deficit remains'))
 }
 if(length(policy_rows)) {
  write_tsv(do.call(rbind,summaries),file.path(study,'tables/baseline_chain_summary.tsv'));write_tsv(do.call(rbind,policy_rows),file.path(study,'tables/baseline_policy_gates.tsv'))
  evidence<-c(evidence,'reference/imports.tsv')
  budget<-read_tsv(file.path(study,'manifests/budget.tsv'));repair_roles<-budget$role[budget$phase=='repair']
  pending_repair<-any(jobs$role%in%repair_roles & !jobs$state%in%c('COMPLETED','FAILED','TIMEOUT','CANCELLED','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'))
  emit_gate('baseline_gate.tsv',length(policy_rows)==2L && allpass && !pending_repair,list(production_pass=if(!is.null(policy_rows$production))policy_rows$production$pass else FALSE,
   data_only_pass=if(!is.null(policy_rows$data_only))policy_rows$data_only$pass else FALSE))
  print(do.call(rbind,policy_rows),row.names=FALSE)
 }
} else if(action%in%c('profile_select','profile_review')) {
 branch_for<-function(task,seen=character()) {
  stopifnot(!task%in%seen);sp<-specs[[task]]
  if(sp$role=='density_global')return('global')
  if(sp$role=='density_global_local')return('global')
  if(sp$role=='density_continuation_local')return(if(grepl('_asc$',task))'asc' else 'desc')
  root<-root_for(task)
  if(sp$mode=='local')source<-readRDS(file.path(root,'actual_local_start.rds'))$resolution$source
  else source<-readRDS(file.path(root,'actual_start.rds'))$source
  stopifnot(source%in%names(specs));branch_for(source,c(seen,task))
 }
 reviews<-list();counterexamples<-list();complete<-TRUE;overall<-TRUE
 for(policy in c('production','data_only'))for(i in 1:7) {
  prefix<-paste0('d_',policy,'_',i);primary<-paste0(prefix,c('_local','_asc','_desc'))
  if(!all(vapply(primary,closed,logical(1)))){complete<-FALSE;next}
  choices<-primary
  for(task in paste0(prefix,c('_reserve_global','_reserve_local','_fresh','_reserve_fresh','_reserve_directions')))
   if(closed(task))choices<-c(choices,task)
  points<-lapply(choices,candidate)
  for(a in points)stopifnot(a$record$parameters[['log10_rho_2N']]==contract$rho_grid[i])
  vals<-vapply(points,function(a)a$record$components$L,numeric(1));idx<-which.min(vals);selected<-points[[idx]]
  selected$branch<-branch_for(choices[idx])
  atomic_rds(selected,file.path(study,'tables',paste0('density_',policy,'_',i,'_selected.rds')))
  if(action=='profile_select')next
  checked_tasks<-paste0(prefix,c('_fresh','_reserve_fresh'));checked_tasks<-checked_tasks[vapply(checked_tasks,closed,logical(1))]
  if(!length(checked_tasks)){complete<-FALSE;next}
  checks<-lapply(checked_tasks,candidate);ck<-which.min(vapply(checks,function(a)a$record$components$L,numeric(1)))
  best<-checks[[ck]];value<-best$record$components$L;result<-read_kv(file.path(payload(checked_tasks[ck]),'result.tsv'))
  scores<-setNames(rep(Inf,3),c('global','asc','desc'))
  for(k in seq_along(points)) {branch<-branch_for(choices[k]);scores[branch]<-min(scores[branch],points[[k]]$record$components$L)}
  selected_verified<-value<=selected$record$components$L+1e-8
  baseline<-readRDS(file.path(study,'tables',paste0('baseline_',policy,'.rds')));delta<-value-baseline$record$components$L
  recovered<-sum(scores-value<=protocol$profile_recovery)
  direction_pass<-result$directional_improvement<=protocol$direction_improvement
  consistent<-delta>=-protocol$baseline_recovery
  pass<-selected_verified && recovered>=2L && direction_pass && consistent
  overall<-overall && pass
  if(!consistent) {
   old<-counterexamples[[policy]]
   if(is.null(old)||value<old$record$components$L)counterexamples[[policy]]<-best
  }
  reviews[[prefix]]<-data.frame(policy=policy,density_index=i,log10_density=contract$rho_grid[i],density=10^contract$rho_grid[i],
   objective=value,data=best$record$components$L_data,prior=best$record$components$L_prior,
   baseline_objective=baseline$record$components$L,raw_delta=delta,recovered_branches=recovered,
   global_branch=unname(scores['global']),ascending_branch=unname(scores['asc']),descending_branch=unname(scores['desc']),
   directional_improvement=result$directional_improvement,selected_verified=selected_verified,baseline_consistent=consistent,pass=pass,fresh_task=checked_tasks[ck])
  atomic_rds(best,file.path(study,'tables',paste0('profile_',policy,'_',i,'_verified.rds')))
 }
 if(action=='profile_select')cat('Selected closed three-branch profile candidates.\n')
 else {
  if(length(reviews))write_tsv(do.call(rbind,reviews),file.path(study,'tables/profile_point_review.tsv'))
  if(length(counterexamples)) {
   oldfile<-file.path(study,'tables/profile_counterexamples.rds')
   old<-if(file.exists(oldfile))readRDS(oldfile) else list()
   for(policy in names(counterexamples))if(is.null(old[[policy]])||counterexamples[[policy]]$record$components$L<old[[policy]]$record$components$L)old[[policy]]<-counterexamples[[policy]]
   atomic_rds(old,oldfile)
   history<-file.path(study,'tables/gate_history');dir.create(history,showWarnings=FALSE)
   gate<-file.path(study,'tables/baseline_gate.tsv');save<-file.path(history,paste0('baseline_before_profile_',hash_file(gate),'.tsv'))
   if(!file.exists(save))stopifnot(file.copy(gate,save))
   evidence<-c(evidence,'tables/profile_counterexamples.rds')
   emit_gate('baseline_gate.tsv',FALSE,list(reason='A constrained canonical candidate beats the free baseline beyond0.002'))
  }
  emit_gate('profile_gate.tsv',complete && length(reviews)==14L && overall,list(complete=complete,reviewed_anchors=length(reviews),baseline_consistent=!length(counterexamples)))
  if(length(reviews))print(do.call(rbind,reviews)[,c('policy','density_index','raw_delta','recovered_branches','directional_improvement','pass')],row.names=FALSE)
 }
} else stop('Unsupported metadata gate action')
assert_namespaces()
