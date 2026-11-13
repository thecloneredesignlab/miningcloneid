# Closed payload / allocation / checkpoint diagnostics; no model evaluation.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'))
jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character');terminal<-c('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')
closed<-sub('[+ ].*$','',jobs$state)%in%terminal;have<-closed&nzchar(jobs$actual_core_seconds)
actual<-sum(as.numeric(jobs$actual_core_seconds[have]));used<-sum(ifelse(have,as.numeric(jobs$actual_core_seconds),as.numeric(jobs$reserved_core_seconds)),na.rm=TRUE)
fit<-checks<-meters<-list()
for(i in which(jobs$state=='COMPLETED'&jobs$exit_code=='0:0')){
 j<-jobs[i,];root<-file.path(study,if(j$mode%in%c('local','global'))'fits' else 'tests',j$task);verify_payload_R(root)
 z<-read_kv(file.path(root,'result.tsv'))
 if(j$mode%in%c('local','global')){
  a<-readRDS(file.path(root,'accepted_candidate.rds'))
  fit[[j$task]]<-data.frame(task=j$task,policy=j$policy,root_id=j$root_id,family=j$root_family,kind=j$kind,mode=j$mode,release=j$release,job_id=j$job_id,
   parent_seconds=as.numeric(j$actual_core_seconds)/as.integer(j$cores),allocated_core_seconds=as.numeric(j$actual_core_seconds),
   objective=a$checked$L,data=a$checked$L_data,prior=a$checked$L_prior,burden=a$checked$L_b,ploidy=a$checked$L_p,necrosis=a$checked$L_n,
   calls=z$evaluations,root_callbacks=z$root_callbacks,window_gain=z$window_gain,termination=z$termination,
   visits_max=if(is.null(z$visits_max))NA_integer_ else z$visits_max,completed_polls=if(is.null(z$completed_polls))NA_integer_ else z$completed_polls,
   stage_seconds=z$wall_seconds,source_task=j$source_task,selector=j$selector,density_index=as.integer(j$density_index))
 }
 if(j$mode%in%c('fresh','directions','profile_fresh'))checks[[j$task]]<-data.frame(task=j$task,policy=j$policy,root_id=j$root_id,mode=j$mode,pass=z$pass,objective=z$objective,score_difference=z$difference,directional_gain=z$directional_improvement,probes=z$probes,fixed_residual=z$fixed_residual)
 meter<-read_kv(file.path(root,'runtime_meter_master.tsv'))
 # Parent counters already include all child/master work; add worker meters once.
 for(path in list.files(root,pattern='^worker_runtime_meters[.]rds$',recursive=TRUE,full.names=TRUE))for(m in readRDS(path))for(name in names(meter))meter[[name]]<-meter[[name]]+m[[name]]
 meters[[j$task]]<-cbind(data.frame(task=j$task,mode=j$mode),as.data.frame(meter))
}
if(length(fit)){
 x<-do.call(rbind,fit);write_tsv(x,file.path(study,'tables/fit_summary.tsv'))
 root_totals<-do.call(rbind,lapply(unique(x$root_id),function(root){z<-x[x$root_id==root,,drop=FALSE];data.frame(root_id=root,policy=z$policy[1],family=z$family[1],closed_windows=nrow(z),optimizer_callbacks=sum(z$calls),allocated_core_seconds=sum(z$allocated_core_seconds),best_closed_objective=min(z$objective))}))
 write_tsv(root_totals,file.path(study,'tables/root_accounting.tsv'))
}
if(length(checks))write_tsv(do.call(rbind,checks),file.path(study,'tables/acceptance_summary.tsv'))
if(length(meters))write_tsv(do.call(rbind,meters),file.path(study,'tables/runtime_model_counts.tsv'))
roles<-do.call(rbind,lapply(unique(jobs$role),function(role){x<-jobs[jobs$role==role,,drop=FALSE];data.frame(role=role,attempts=nrow(x),actual_core_seconds=sum(as.numeric(x$actual_core_seconds[nzchar(x$actual_core_seconds)])),live_core_seconds=sum(as.numeric(x$reserved_core_seconds[!sub('[+ ].*$','',x$state)%in%terminal])))}))
write_tsv(roles,file.path(study,'tables/role_accounting.tsv'))
flag<-function(name){p<-file.path(study,'tables',name);if(!file.exists(p))'PENDING' else if(read_kv(p)$pass)'PASS' else 'NOT ACCEPTED'}
origin_path<-file.path(study,'tables/analysis_snapshot_origin.tsv');origin<-if(file.exists(origin_path))read.delim(origin_path,stringsAsFactors=FALSE) else NULL
snapshot_note<-if(is.null(origin))'The tables below describe the collected Slurm snapshot.' else paste0('Scientific/timing tables use frozen Slurm generation ',origin$value[origin$key=='state_generation'],' captured ',origin$value[origin$key=='captured_utc'],'. The enclosing validation allocation is billed after it closes; the current owned ledger/STATUS records later accounting.')
runtime_release<-read.delim(file.path(study,'tables/execution_release.tsv'),stringsAsFactors=FALSE)$value[1]
lines<-c('# Improvement-driven in vivo continuation report','',paste('Updated',format(Sys.time(),tz='UTC',usetz=TRUE)),'',snapshot_note,
 '',paste0(nrow(jobs),' owned allocations; ',sum(closed),' closed; ',sum(!closed),' active/uncertain.'),
 '',paste0('New closed actual ',actual,' core-seconds (',sprintf('%.8f',actual/3600),' core-hours). Actual plus live ',used,'. Historical lifetime301846/formal271355remain charged. Lifetime actual ',301846+actual,' (',sprintf('%.8f',(301846+actual)/3600),' hours); historical formal plus new actual ',271355+actual,' (',sprintf('%.8f',(271355+actual)/3600),' hours).'),
 '', 'No aggregate core-hour or role/attempt ceiling applies. Per-window call/time limits are interruptions, not convergence. Further windows depend on observed improvement, unresolved checks and scientific gates. Estimates22/74/171.5hours are planning values only. All failed/cancelled/partial parents remain charged.',
 '', 'Per-window root_callbacks fields describe the selected continuous optimizer ancestry; an explicit endpoint restart/fork creates a new component. The authoritative non-resetting root totals sum distinct owned window callbacks in tables/root_accounting.tsv, without double-counting inherited path counters. Canonical replay/export score children and completed native trajectories are separately metered.',
 '',paste0('Engineering ',flag('engineering_gate.tsv'),'; both unrestricted baselines ',flag('baseline_gate.tsv'),'; profiles ',flag('profile_gate.tsv'),'.'),
 '', 'All actual model evaluations use the exact SHA256-pinned R4.4.2 SIF, unchanged numerical5/raw3inputs, original masks/weights/priors/bounds, fixed N0=Nref=1e6 and clean-parent fresh-child isolation. Private runtime/cache sources are owned by this study. Culture, shared model/runner/image/libraries/cache sources and closed prior studies remain unchanged.',
 '', 'The uncapped multiscale state machine cycles across four scales without a fixed sweep/cycle/window count. Analytic43-visit fixtures and real split/uninterrupted tests verify continuation. A poll_stationary label denotes finite checks with no effective improvement, not global optimality. Partial native DE blocks retain all best points/calls; population continuation is a warm restart, not bitwise uninterrupted DE.',
 '',paste0('Current science execution release: ',runtime_release,'. Six cont03serial/replay proofs carry by explicit unchanged local-method/branch identity. cont06parallel validation covers eight policy/initialization/map cases, aggregate37tokens,131+80=211 and97+80=177 partial-population continuations, projected-density incumbent exclusion and a1-callback incumbent-retention case. Earlier cont01dispatch-float failures, obsolete cont02NFS-token validation cancellation, cont03worker cache-load failure and superseded cont05validation remain in the reviewed allocation ledger. No numerical model/domain/acceptance rule was changed.'),
 '', 'Full native/checkpoint/record payloads are authoritative on HPC. The local lean mirror retains source/control manifests, reports, typed summaries, per-parent remote payload fingerprints and verified metadata packets. Local storage exhaustion was repaired by evicting verified current-study download replicas. No canonical HPC result or unrelated study was removed.',
 '', 'Archived-output validation runs in allocated single-core jobs in the same pinned SIF, with score/prediction calls explicitly forbidden. These jobs add no optimizer/model calls and their parent allocation time remains charged in the same ledger. Actual fitting, fresh model replays and coordinate checks remain separately allocated.',
 '', 'Independent baseline recovery requires archive and genuinely uniform roots within0.002, canonical identity1e-8, known-feasible retention and independent original-coordinate direction gain<=0.002. Root extensions count once. No study/archive winner is injected into uniform histories. Density profiles remain conditional on BOTH baselines, keep19nuisancesfree and requiretwo-of-threebranchagreement0.01. No calibrated likelihood-ratio interval is claimed.')
if(file.exists(file.path(study,'tables/baseline_policy_gates.tsv'))){
 x<-read_tsv(file.path(study,'tables/baseline_policy_gates.tsv'));lines<-c(lines,'','## Independent recovery','','| Policy | Best verified Q | Known floor | Uniform gap | Accepted |','|---|---:|---:|---:|---|')
 for(i in seq_len(nrow(x)))lines<-c(lines,paste0('| ',x$policy[i],' | ',sprintf('%.12f',x$best_objective[i]),' | ',sprintf('%.12f',x$known_feasible_bound[i]),' | ',sprintf('%.12f',x$uniform_gap[i]),' | ',x$pass[i],' |'))
}
if(length(fit)){
 x<-do.call(rbind,fit);lines<-c(lines,'','## Window timing and improvement','','| Task | Parent minutes | Calls | Window gain | Q | Termination | Maximum scale visits |','|---|---:|---:|---:|---:|---|---:|')
 for(i in seq_len(nrow(x)))lines<-c(lines,paste0('| ',x$task[i],' | ',sprintf('%.3f',x$parent_seconds[i]/60),' | ',x$calls[i],' | ',sprintf('%.12g',x$window_gain[i]),' | ',sprintf('%.12f',x$objective[i]),' | ',x$termination[i],' | ',x$visits_max[i],' |'))
}
lines<-c(lines,'','## Inspectable evidence','', '[Fit windows](../tables/fit_summary.tsv), [fresh/direction checks](../tables/acceptance_summary.tsv), [model counters](../tables/runtime_model_counts.tsv), [roles/accounting](../tables/role_accounting.tsv), [owned job ledger](../tables/jobs.tsv). Native predictions, complete improvements, checkpoint generations and polls are retained in closed fit/test payloads. Job parent allocation is the accounting authority; child elapsed timers are diagnostic and are not added again as charged core-hours.')
writeLines(lines,file.path(study,'reports/continuation_report.md'));assert_namespaces();cat('Closed metadata report updated; zero model calls.\n')
