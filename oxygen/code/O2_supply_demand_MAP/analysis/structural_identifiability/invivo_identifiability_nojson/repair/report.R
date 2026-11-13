# Reproducibility, root recovery, timing and integer allocation report; no model calls.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study<-normalizePath(args[1]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'))
source(file.path(entry,'pattern.R'))
jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character');budget<-read_tsv(file.path(study,'manifests/budget.tsv'))
specs<-readRDS(file.path(study,'manifests/task_specs.rds'));contract<-readRDS(file.path(study,'manifests/contract.rds'))
closed<-jobs$state%in%c('COMPLETED','FAILED','TIMEOUT','CANCELLED','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')
charge<-as.integer(jobs$reserved_core_seconds);have<-closed & nzchar(jobs$actual_core_seconds);charge[have]<-as.integer(jobs$actual_core_seconds[have])
actual<-sum(as.integer(jobs$actual_core_seconds[have]));stopifnot(192497+sum(charge)<=648000,222988+sum(charge)<=691200)
roles<-do.call(rbind,lapply(seq_len(nrow(budget)),function(i) {
 b<-budget[i,];take<-jobs$role==b$role
 data.frame(role=b$role,phase=b$phase,attempts=sum(take),attempt_cap=b$attempts,closed=sum(take & closed),
  actual_core_seconds=sum(as.integer(jobs$actual_core_seconds[take & have])),used_plus_live_core_seconds=sum(charge[take]),cap_core_seconds=b$cap_core_seconds)
}))
stopifnot(all(roles$attempts<=roles$attempt_cap),all(roles$used_plus_live_core_seconds<=roles$cap_core_seconds),
 sum(roles$used_plus_live_core_seconds[roles$phase=='repair'])<=129600,sum(roles$used_plus_live_core_seconds[roles$phase=='profile'])<=316800)
write_tsv(roles,file.path(study,'tables/role_accounting.tsv'))
fit<-params<-acceptance<-meters<-schedules<-list()
for(task in jobs$task[closed & jobs$state=='COMPLETED']) {
 sp<-specs[[task]];scope<-if(sp$mode%in%c('global','local'))'fits' else 'tests';root<-file.path(study,scope,task)
 if(!file.exists(file.path(root,'COMPLETE')))next
 verify_payload_R(root);j<-jobs[jobs$task==task,,drop=FALSE];result<-read_kv(file.path(root,'result.tsv'))
 if(sp$mode%in%c('global','local')) {
  a<-readRDS(file.path(root,'accepted_candidate.rds'));z<-a$checked;actual_spec<-readRDS(file.path(root,'actual_spec.rds'))
  logged_generation<-NA_integer_
  if(sp$mode=='global') {
   log<-file.path(study,'logs',paste0(task,'_',j$job_id,'.log'))
   if(file.exists(log)) {
    iterations<-grep('^Iteration: [0-9]+ ',readLines(log,warn=FALSE),value=TRUE)
    if(length(iterations))logged_generation<-as.integer(sub('^Iteration: ([0-9]+).*','\\1',tail(iterations,1)))
   }
  }
  fit[[task]]<-data.frame(task=task,policy=sp$policy,role=sp$role,phase=budget$phase[match(sp$role,budget$role)],
   root_id=actual_spec$root_id,root_family=actual_spec$root_family,release=j$release,job_id=j$job_id,
   parent_seconds=as.integer(j$actual_core_seconds)/as.integer(j$cores),allocated_core_seconds=as.integer(j$actual_core_seconds),stage_seconds=result$wall_seconds,
   objective=z$L,data=z$L_data,prior=z$L_prior,burden=z$L_b,ploidy=z$L_p,necrosis=z$L_n,
   calls=result$evaluations,budget_limited=result$budget_limited,call_limited=if(is.null(result$call_limited))FALSE else result$call_limited,
   native_generations=if(sp$mode=='global')result$iterations else NA_integer_,last_logged_generation=logged_generation,method=if(sp$mode=='global')'DEoptim' else result$method,
   termination=if(sp$mode=='global')if(result$budget_limited)'soft_time_limit' else 'native_stop' else result$termination,
   fixed_residual=result$fixed_residual,density_index=sp$density_index)
  p<-a$record$parameters;natural<-p;use<-startsWith(names(p),'log10_');natural[use]<-10^p[use];position<-(p-contract$lower)/(contract$upper-contract$lower)
  params[[task]]<-data.frame(task=task,policy=sp$policy,parameter=names(p),transformed=unname(p),natural=unname(natural),
   normalized_position=unname(position),near_lower_bound=unname(position)<=1e-6,near_upper_bound=unname(position)>=1-1e-6,
   bound_flag_normalized_tolerance=1e-6,fixed=names(p)%in%names(sp$fixed))
  if(sp$mode=='local') {
   native<-readRDS(file.path(root,'native_result.rds'));events<-readRDS(file.path(root,'polish_evaluations.rds'))
   dynamic<-setdiff(a$record$free_names,c('log10_rho_2N','log10_sigma_burden'))
   directions<-pattern_directions(dynamic);span<-contract$upper-contract$lower
   groups<-list();last_center<-NULL;last_step<-NA_real_;best_so_far<-result$start_objective
   for(event in events) {
    if(event$kind=='poll') {
     center<-event$parameters[dynamic]-event$step*span[dynamic]*directions[[event$direction]]
     new_group<-is.null(last_center)||event$step!=last_step||max(abs(center-last_center))>1e-12
     if(new_group) {
      groups[[length(groups)+1L]]<-list(step=event$step,before=best_so_far,best=event$objective)
      last_center<-center;last_step<-event$step
     } else groups[[length(groups)]]$best<-min(groups[[length(groups)]]$best,event$objective)
    }
    best_so_far<-min(best_so_far,event$objective)
   }
   rows<-lapply(unique(vapply(groups,`[[`,numeric(1),'step')),function(step) {
    take<-groups[vapply(groups,function(z)z$step==step,logical(1))];last<-tail(take,1)[[1]]
    improving<-last$best<last$before-1e-10
    data.frame(task=task,step=step,sweeps=length(take),sweep_cap=native$max_sweeps_per_scale,
     last_recorded_sweep_improved=improving,scale_completed=step%in%native$completed_scales,
     improving_at_sweep_cap=length(take)==native$max_sweeps_per_scale && improving && step%in%native$completed_scales)
   })
   if(length(rows))schedules[[task]]<-do.call(rbind,rows)
  }
 }
 if(sp$mode%in%c('fresh','directions','profile_fresh'))acceptance[[task]]<-data.frame(task=task,policy=sp$policy,mode=sp$mode,
  role=sp$role,job_id=j$job_id,transport_pass=result$pass,objective=result$objective,score_difference=result$difference,
  fixed_residual=result$fixed_residual,probes=result$probes,directional_improvement=result$directional_improvement,allocated_core_seconds=as.integer(j$actual_core_seconds))
 meter<-file.path(root,'runtime_meter_master.tsv')
 if(file.exists(meter)) {
  all<-read_kv(meter);children<-if(sp$mode=='global')list(root) else if(sp$mode=='engineering_parallel')as.list(list.dirs(file.path(study,'tests'),recursive=FALSE,full.names=TRUE)[startsWith(basename(list.dirs(file.path(study,'tests'),recursive=FALSE,full.names=TRUE)),paste0(task,'_'))]) else list()
  for(child in children) {
   p<-file.path(child,'worker_runtime_meters.rds')
   if(file.exists(p))for(m in readRDS(p))for(key in names(all))all[[key]]<-all[[key]]+m[[key]]
  }
  meters[[task]]<-cbind(data.frame(task=task,mode=sp$mode),as.data.frame(all))
 }
}
if(length(fit)){write_tsv(do.call(rbind,fit),file.path(study,'tables/fit_summary.tsv'));write_tsv(do.call(rbind,params),file.path(study,'tables/fit_parameters.tsv'))}
if(length(acceptance))write_tsv(do.call(rbind,acceptance),file.path(study,'tables/acceptance_summary.tsv'))
if(length(meters))write_tsv(do.call(rbind,meters),file.path(study,'tables/runtime_model_counts.tsv'))
if(length(schedules))write_tsv(do.call(rbind,schedules),file.path(study,'tables/local_poll_schedule.tsv'))
gate<-function(name,missing) {
 p<-file.path(study,'tables',name);if(!file.exists(p))return(missing)
 if(isTRUE(read_kv(p)$pass))'PASS' else 'NOT ACCEPTED'
}
lines<-c('# In vivo baseline repair and conditional profile report','',paste0('Updated ',format(Sys.time(),tz='UTC',usetz=TRUE),'.'),'',
 paste0(nrow(jobs),' owned allocations; ',sum(closed),' closed; ',sum(!closed),' active/uncertain.'),'',
 paste0('New closed actual ',actual,' core-seconds (',sprintf('%.8f',actual/3600),' core-hours); actual plus live ',sum(charge),'. Prior lifetime debit 222988 and formal debit 192497 remain charged. Lifetime actual plus live ',222988+sum(charge),' (',sprintf('%.8f',(222988+sum(charge))/3600),' hours); formal actual plus live ',192497+sum(charge),' (',sprintf('%.8f',(192497+sum(charge))/3600),' hours). Repair/profile phase limits 129600 / 316800; formal/lifetime 648000 / 691200.'),'',
 paste0('Engineering: ',gate('engineering_gate.tsv','PENDING'),'. Both unrestricted baselines: ',gate('baseline_gate.tsv','PENDING'),'. Density profile: ',gate('profile_gate.tsv','NOT STARTED'),'.'),'',
 'All numerical work uses the exact SHA256-pinned R4.4.2 SIF and unchanged five-file model, raw observations, weights, priors, initial stocks and demand reference. Clean-parent fresh children score each point; exported predictions are separately prepared canonically. Culture, shared source/runner/image/libraries and manuscript assets are protected. No commit/push.','',
 'repair02 is a pure lineage-metadata correction of repair01: a dynamic reserve parent is validated using its actual dispatch root. All other ten execution/method files are byte-identical. Its proof combines the original SIF engineering gate and a labeled pure metadata fixture; no additional numerical engineering allocation is claimed. Existing repair01 jobs are unchanged.','',
 'Independent recovery requires archive and uniform roots within 0.002 and feasible original-coordinate probes improving no more than 0.002. Uniform roots cannot receive the shared archive winner. A root extension is not a new independent start. The protected known-feasible pool includes all verified engineering points, observation proposals and any later constrained counterexamples. Exact score/prediction identity checks are distinct from search recovery.','',
 'Primary polishing uses normalized bounded coordinate/coupled-block polls and canonical-verified observation proposals. All 20 coordinates remain free in baselines; density is held exactly only for 19-coordinate profile branches. Call/time truncation and completed finite polls do not establish global optimality. Profile deltas are retained without clipping, and no confidence interval or likelihood-ratio threshold is calibrated.','',
 'Allocation is charged once from Slurm parent AllocCPUS times ElapsedRaw, including failures and launch/compilation. Point counters distinguish score children and prediction children; completed trajectory calls include the canonical score preparation used before prediction export. Stage counters alone do not represent all charged work. Private caches and native outputs remain archived.')
if(file.exists(file.path(study,'tables/baseline_policy_gates.tsv'))) {
 x<-read_tsv(file.path(study,'tables/baseline_policy_gates.tsv'))
 lines<-c(lines,'','## Baseline recovery','','The selected incumbent is the lowest canonically verified point, including a better coordinate probe. The independent direction-gain test retains the improvement measured from the optimized endpoint; selection does not reset that gain to zero. Diagnostic centers, all evaluated vectors and selected-point checks are retained.','',
  '| Policy | Best verified Q | Known feasible bound | Uniform gap | Family recovery | Direction check | Accepted |','|---|---:|---:|---:|---|---|---|')
 for(i in seq_len(nrow(x)))lines<-c(lines,paste0('| ',paste(c(x$policy[i],sprintf('%.10f',x$best_objective[i]),sprintf('%.10f',x$known_feasible_bound[i]),sprintf('%.10f',x$uniform_gap[i]),x$family_recovery[i],x$directional_pass[i],x$pass[i]),collapse=' | '),' |'))
 imports<-read_tsv(file.path(study,'reference/imports.tsv'))
 for(i in seq_len(nrow(x))) {
  original<-min(imports$expected_objective[imports$target_policy==x$policy[i]])
  lines<-c(lines,'',paste0(x$policy[i],': the freshly reproduced imported feasible anchor was ',sprintf('%.15f',original),
   '; the selected verified point is ',sprintf('%.15f',x$best_objective[i]),', an improvement of ',sprintf('%.15f',original-x$best_objective[i]),
   '. Thus the higher-than-known-point defect is repaired, while independent root recovery remains a separate requirement.'))
 }
}
if(length(acceptance)) {
 x<-do.call(rbind,acceptance);fr<-x[x$mode=='fresh',]
 lines<-c(lines,'','## Separate acceptance checks','',paste0(nrow(fr),' fresh checks; maximum absolute fresh gap ',if(nrow(fr))format(max(abs(fr$score_difference)),digits=17) else 'not available','; ',sum(x$mode!='fresh'),' direction/anchor checks. Exact fixed-target residuals and counts are in tables/acceptance_summary.tsv.'))
}
if(length(fit)) {
 x<-do.call(rbind,fit)
 lines<-c(lines,'','## Optimization timing','',
  'Parent wall time includes launch and SIF/cache preparation. Stage time and callback calls are method diagnostics; the charged resource measure is the parent allocation. DE time/call truncation remains visible. A completed bounded pattern schedule is a finite search result, not a convergence certificate.','',
  '| Task | Parent minutes | Stage minutes | Callback calls | Q | Termination |',
  '|---|---:|---:|---:|---:|---|')
 for(i in seq_len(nrow(x)))lines<-c(lines,paste0('| ',x$task[i],' | ',sprintf('%.3f',x$parent_seconds[i]/60),
  ' | ',sprintf('%.3f',x$stage_seconds[i]/60),' | ',x$calls[i],' | ',sprintf('%.10f',x$objective[i]),' | ',x$termination[i],' |'))
 local<-x$task[x$method=='bounded_block_pattern'];gaps<-numeric();count<-0L
 for(task in local) {
  proposals<-readRDS(file.path(study,'fits',task,'observation_proposals.rds'))
  count<-count+length(proposals)
  gaps<-c(gaps,vapply(proposals,function(z)abs(z$difference),numeric(1)))
 }
 lines<-c(lines,'',paste0(count,' observation-block proposals in completed optimization locals; maximum absolute full canonical score discrepancy ',
  if(length(gaps))format(max(gaps),digits=17) else 'not applicable','. The unchanged canonical scorer evaluates every proposed full vector before acceptance.'))
 if(length(schedules)) {
  schedule<-do.call(rbind,schedules);limited<-unique(schedule$task[schedule$improving_at_sweep_cap])
  lines<-c(lines,'',paste0(length(limited),' completed local schedules had at least one scale still improving on its 25th and final allowed sweep. These are finite sweep-cap stops, even when termination is recorded as complete. Exact scale/sweep counts are reconstructed from the archived poll centers in tables/local_poll_schedule.tsv; they do not add model calls.'))
 }
}
lines<-c(lines,'','## Inspectable outputs','',
 'The 36-core-hour repair allowance is a reservation ceiling with separate frozen attempt caps. Jobs release unused wall time on their scheduled finite poll or soft-time stops. Unused hours are not transferred to additional attempts. All submitted repair parents are closed before any profile release; the original 88-core-hour profile allowance remains held when either baseline fails.',
 '',
 '[Recovery figure](../figures/baseline_recovery.png); [fit timing and scores](../tables/fit_summary.tsv); [local poll caps](../tables/local_poll_schedule.tsv); [model-call meters](../tables/runtime_model_counts.tsv); [fresh and direction checks](../tables/acceptance_summary.tsv); [selected score components](../tables/selected_candidate_components.tsv); [selected parameters](../tables/selected_candidate_parameters.tsv).',
 '', 'Final integrity and accounting: [closure audit](../tables/closure_audit.tsv), [delivery manifest](../tables/final_delivery_manifest.tsv), [remote transport verification](../tables/delivery_transport_verification.tsv). These closure artifacts are available only after all submitted jobs close.',
 '','## Role accounting','','| Role | Attempts / cap | Actual hours | Used plus live hours | Cap hours |','|---|---:|---:|---:|---:|')
for(i in seq_len(nrow(roles)))lines<-c(lines,paste0('| ',roles$role[i],' | ',roles$attempts[i],' / ',roles$attempt_cap[i],' | ',sprintf('%.8f',roles$actual_core_seconds[i]/3600),' | ',sprintf('%.8f',roles$used_plus_live_core_seconds[i]/3600),' | ',sprintf('%.8f',roles$cap_core_seconds[i]/3600),' |'))
audits<-list.files(file.path(study,'tables/record_audits'),pattern='_summary[.]tsv$',full.names=TRUE)
if(length(audits)) {
 a<-lapply(audits,read_kv);lines<-c(lines,'','## Archived improvement audit','',paste0(sum(vapply(a,`[[`,integer(1),'records')),' full-precision improvements audited across ',length(a),' payloads; maximum independent component gap ',format(max(vapply(a,`[[`,numeric(1),'independent_component_max_gap')),digits=17),'. No model calls used in this metadata audit.'))
}
writeLines(lines,file.path(study,'reports/repair_execution_report.md'));assert_namespaces();cat('Repair metadata report updated; zero model calls.\n')
