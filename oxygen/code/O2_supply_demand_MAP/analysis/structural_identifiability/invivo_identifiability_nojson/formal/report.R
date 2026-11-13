# Generate source-backed interim/final tables and static figures without model calls.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"));source(file.path(entry,"helpers.R"))
jobs <- read.delim(file.path(study,"tables/jobs.tsv"),colClasses="character")
raw <- read.delim(file.path(study,"tables/raw_accounting.txt"),sep="|",header=FALSE,quote="",fill=TRUE,colClasses="character")
contract <- readRDS(file.path(study,"manifests/contract.rds"));specs <- readRDS(file.path(study,"manifests/task_specs.rds"))
closed <- jobs$state%in%c("COMPLETED","FAILED","TIMEOUT","CANCELLED","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE","REVOKED")
charge <- as.integer(jobs$reserved_core_seconds)
charge[closed & nzchar(jobs$actual_core_seconds)] <- as.integer(jobs$actual_core_seconds[closed & nzchar(jobs$actual_core_seconds)])
actual <- sum(as.integer(jobs$actual_core_seconds[closed & nzchar(jobs$actual_core_seconds)]))
stopifnot(sum(charge)<=648000,30491+sum(charge)<=691200)
budget <- read_tsv(file.path(study,"manifests/budget.tsv"))
stopifnot(all(jobs$role%in%budget$role))
role <- do.call(rbind,lapply(budget$role,function(r) {
 take <- jobs$role==r
 cap <- budget[budget$role==r,,drop=FALSE]
 data.frame(role=r,attempts=sum(take),closed=sum(take & closed),
  attempt_cap=cap$attempts,cap_core_seconds=cap$cap_core_seconds,
  actual_core_seconds=sum(as.integer(jobs$actual_core_seconds[take & closed & nzchar(jobs$actual_core_seconds)])),
  used_plus_committed_core_seconds=sum(charge[take]))
}))
stopifnot(sum(role$actual_core_seconds)==actual,sum(role$used_plus_committed_core_seconds)==sum(charge),
 all(role$attempts<=role$attempt_cap),all(role$used_plus_committed_core_seconds<=role$cap_core_seconds))
write_tsv(role,file.path(study,"tables/role_accounting.tsv"))
fit_rows <- list();parameter_rows <- list()
for(task in jobs$task[closed & jobs$state=="COMPLETED"]) {
 if(!specs[[task]]$mode%in%c("global","local")) next
 root <- file.path(study,"fits",task)
 if(!file.exists(file.path(root,"COMPLETE"))) next
 verify_payload_R(root)
 result <- read_kv(file.path(root,"result.tsv"));accepted <- readRDS(file.path(root,"accepted_candidate.rds"))
 j <- jobs[jobs$task==task,,drop=FALSE];parent <- raw[raw$V1==j$job_id,,drop=FALSE];batch <- raw[raw$V1==paste0(j$job_id,".batch"),,drop=FALSE]
 stopifnot(nrow(parent)==1L)
 c <- accepted$checked;spec <- readRDS(file.path(root,"actual_spec.rds"))
 fit_rows[[task]] <- data.frame(task=task,policy=spec$policy,role=spec$role,family=spec$family,
  job_id=j$job_id,cores=as.integer(j$cores),parent_seconds=as.integer(parent$V5),
  allocated_core_seconds=as.integer(j$cores)*as.integer(parent$V5),stage_seconds=result$wall_seconds,
  objective=c$L,data=c$L_data,prior=c$L_prior,burden=c$L_b,ploidy=c$L_p,necrosis=c$L_n,
  calls=result$evaluations,budget_limited=result$budget_limited,
  native_generations=if(spec$mode=="global") result$iterations else NA_integer_,
  native_local_code=if(spec$mode=="local") result$native_convergence else NA_integer_,
  density_index=spec$density_index,fixed_residual=result$fixed_residual,
  slurm_batch_maxrss_reported=if(nrow(batch)==1L) batch$V6 else "")
 p <- accepted$record$parameters;natural <- p;use <- startsWith(names(p),"log10_");natural[use] <- 10^p[use]
 position <- unname((p-contract$lower)/(contract$upper-contract$lower))
 parameter_rows[[task]] <- data.frame(task=task,policy=spec$policy,parameter=names(p),transformed=unname(p),natural=unname(natural),
  transformed_lower=unname(contract$lower),transformed_upper=unname(contract$upper),
  normalized_position=position,near_lower_bound=position<=1e-6,near_upper_bound=position>=1-1e-6,
  bound_flag_normalized_tolerance=1e-6,fixed=names(p)%in%names(spec$fixed))
}
if(length(fit_rows)) {
 write_tsv(do.call(rbind,fit_rows),file.path(study,"tables/fit_summary.tsv"))
 write_tsv(do.call(rbind,parameter_rows),file.path(study,"tables/fit_parameters.tsv"))
}
acceptance <- list()
for(task in jobs$task[closed & jobs$state=="COMPLETED" & jobs$mode%in%c("fresh","directions","profile_fresh")]) {
 root <- file.path(study,"tests",task)
 if(!file.exists(file.path(root,"COMPLETE"))) next
 verify_payload_R(root);z <- read_kv(file.path(root,"result.tsv"));j <- jobs[jobs$task==task,,drop=FALSE]
 acceptance[[task]] <- data.frame(task=task,role=j$role,mode=j$mode,policy=z$policy,job_id=j$job_id,
  transport_pass=z$pass,objective=z$objective,score_difference=z$difference,fixed_residual=z$fixed_residual,
  probes=z$probes,directional_improvement=z$directional_improvement,
  allocated_core_seconds=as.integer(j$actual_core_seconds),source_task=z$source_task)
}
if(length(acceptance)) write_tsv(do.call(rbind,acceptance),file.path(study,"tables/acceptance_summary.tsv"))
baseline_pass <- FALSE;profile_pass <- FALSE
baseline_label <- "PENDING";profile_label <- "NOT STARTED"
if(file.exists(file.path(study,"tables/baseline_gate.tsv"))) {
 baseline_pass <- isTRUE(read_kv(file.path(study,"tables/baseline_gate.tsv"))$pass)
 baseline_label <- if(baseline_pass) "PASS" else "NOT ACCEPTED"
}
if(file.exists(file.path(study,"tables/profile_gate.tsv"))) {
 g <- read_kv(file.path(study,"tables/profile_gate.tsv"));profile_pass <- isTRUE(g$pass)
 profile_label <- if(profile_pass) "PASS" else if(isTRUE(g$complete)) "NOT ACCEPTED" else "INCOMPLETE"
}
lines <- c("# Formal tumor identifiability execution report","",
 paste0("Updated ",format(Sys.time(),tz="UTC",usetz=TRUE),"."),"",
 paste0("Owned allocations: ",nrow(jobs),"; closed: ",sum(closed),"; active/uncertain: ",sum(!closed),"."),
 paste0("Formal actual allocation: ",actual," integer core-seconds (",sprintf("%.8f",actual/3600)," core-hours). Formal actual plus live reservations: ",sum(charge)," core-seconds. Prior pilot debit: 30491. Lifetime actual plus live reservations: ",30491+sum(charge)," core-seconds (",sprintf("%.8f",(30491+sum(charge))/3600)," core-hours). Caps remain 648000 new / 691200 lifetime core-seconds."),"",
 paste0("Both-policy unrestricted baseline gate: ",baseline_label,". Complete density-profile gate: ",profile_label,"."),"",
 "The exact five-file numerical source, observation masks, initial stocks, demand reference and pinned R4.4.2 SIF remain unchanged. Formal bounds were transformed directly from the original natural parameter input and frozen in RDS. Scores use clean-parent child evaluations and the declared within-objective scenario order. Production and data-only policies retain separate optimization and baseline references.","",
 "Every reported search endpoint is a best observed feasible archived incumbent. Native local codes, DE generation limits, budget-limited flags and independent score acceptance are separate fields. Source, repeatability, fixed-target and complete-prediction mismatches stop dependent work. A finite recovery/poll gate is not a proof of global optimality.","",
 "Engineering: eight actual DE cases, four serial cases, ten reference points, four real direction fixtures and exact typed/analytic checks passed in the pinned SIF. Five engineering parents used 2222 core-seconds. Control-plane mocks verified minute-resolution time requests, dual budget caps, duplicate-callback idempotence, stale payload rejection and reservation retention after a lost receipt. Mock tests launched no real Slurm work.","",
 "The baseline design has four independent starts per policy (seed25/392 hybrid neighborhoods and two full-box uniform populations without an endpoint insertion). Full twenty-coordinate search uses NP=256, strategy=2, maximum 1000 generations, reltol=1e-6 and patience=100, followed by separate native L-BFGS-B polishing and independent fresh/direction checks. Cross-family recovery must be within 0.002; feasible coordinate probes at normalized spans 1e-3 and 1e-4 must not improve the best by more than 0.002. Exact score identity remains 1e-8. All formally replayed reference points and the compatible corrected pilot constrained points remain known feasible objective bounds. Agreement at a higher score cannot override a better known feasible point.","",
 "Profiles may start only after both matching unrestricted baselines pass. The seven positive-density anchors cover 50000-500000 cells/mm^3; all nineteen active nuisance coordinates are refitted with a global branch and ascending/descending local traversals. Each anchor requires two branches within 0.01 of the best verified score and no nuisance-poll improvement above 0.002. A constrained score beating its baseline by more than 0.002 invalidates the baseline gate and blocks further profile submission. Raw deltas are retained, including negative values.","",
 "Density is conditional on fixed initial live abundance 1e6, initial chromosome diracs 44/88 and fixed demand reference 1e6. Effective founding abundance is not established by nominal injected cells. Profile objectives are weighted/composite diagnostics; no automatic likelihood-ratio cutoff, p value or confidence interval is produced. No initial-size extension or joint culture/tumor fit is included. Parameter bound flags mean within 1e-6 of a transformed-box span; these are search-domain flags, not measured constraints.","",
 "Allocation is charged once from Slurm parents, including failed/partial work; raw batch/extern accounting is retained without double charging. Stage timers include in-stage acceptance/bookkeeping, while parent elapsed includes launch/compilation. Slurm batch MaxRSS is exported verbatim, without claiming aggregate peak worker memory. Culture source/jobs/budget, shared model/libraries/image/runner and manuscript assets were not changed. No commit or push.")
if(file.exists(file.path(study,"tables/baseline_policy_gates.tsv"))) {
 x <- read_tsv(file.path(study,"tables/baseline_policy_gates.tsv"))
 lines <- c(lines,"","## Unrestricted recovery","","| Policy | Best objective | Best probe | Known feasible bound | Cross-family recovery | Direction gate | Accepted |","|---|---:|---:|---:|---|---|---|")
 for(i in seq_len(nrow(x))) lines <- c(lines,paste0("| ",paste(c(x$policy[i],sprintf("%.10f",x$best_objective[i]),sprintf("%.10f",x$best_probe[i]),sprintf("%.10f",x$known_feasible_bound[i]),x$family_recovery[i],x$directional_pass[i],x$pass[i]),collapse=" | ")," |"))
}
if(length(acceptance)) {
 checks <- do.call(rbind,acceptance);fresh <- checks[checks$mode=="fresh",,drop=FALSE]
 lines <- c(lines,"","## Independent acceptance evidence","",
  paste0(nrow(fresh)," separate fresh-check allocations are closed; their maximum absolute score gap is ",
   if(nrow(fresh)) format(max(abs(fresh$score_difference)),digits=17) else "not available",". ",
   sum(checks$mode%in%c("directions","profile_fresh"))," separate directional/anchor checks are closed. ",
   "Their transport, score differences, exact fixed-target residuals, probe counts and allocations are retained in tables/acceptance_summary.tsv. A successful transport check does not replace the search-recovery gate."))
}
density_jobs <- sum(startsWith(jobs$role,"density_") | startsWith(jobs$role,"reserve_density_"))
if(all(closed) && !baseline_pass && density_jobs==0L) lines <- c(lines,"",
 "All owned allocations are closed. The unrestricted baseline quality gate remains unresolved, so the approved Stage A stopping rule leaves Stage B unlaunched. The selected reserve attempts remain charged. Unused profile capacity is not reassigned to more baseline attempts.")
lines <- c(lines,"","## Allocation by approved role","",
 "| Role | Attempts / cap | Actual core-hours | Actual plus live core-hours | Role cap core-hours |",
 "|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(role))) lines <- c(lines,paste0("| ",role$role[i]," | ",role$attempts[i]," / ",role$attempt_cap[i]," | ",
 sprintf("%.8f",role$actual_core_seconds[i]/3600)," | ",sprintf("%.8f",role$used_plus_committed_core_seconds[i]/3600)," | ",
 sprintf("%.8f",role$cap_core_seconds[i]/3600)," |"))
if(length(fit_rows)) {
 fits <- do.call(rbind,fit_rows)
 lines <- c(lines,"","## Numerical stage timing","",
  "| Stage | Closed fits | Stage seconds, minimum | Stage seconds, maximum | Objective calls | Budget-limited fits |",
  "|---|---:|---:|---:|---:|---:|")
 for(mode in c("global","local")) {
  keep <- vapply(fits$task,function(task) specs[[task]]$mode==mode,logical(1));z <- fits[keep,,drop=FALSE]
  if(nrow(z)) lines <- c(lines,paste0("| ",mode," | ",nrow(z)," | ",sprintf("%.3f",min(z$stage_seconds))," | ",
   sprintf("%.3f",max(z$stage_seconds))," | ",sum(z$calls)," | ",sum(z$budget_limited)," |"))
 }
 lines <- c(lines,"","Per-fit native codes, exact full-precision score components, parent elapsed, stage timers, calls and reported batch MaxRSS are in tables/fit_summary.tsv. Per-role allocation includes launch and all acceptance costs; the stage table alone is not the full budget.")
}
audits <- list.files(file.path(study,"tables/record_audits"),pattern="_summary[.]tsv$",full.names=TRUE)
if(length(audits)) {
 a <- lapply(audits,read_kv)
 lines <- c(lines,"","## Archived improvement audit","",
  paste0(sum(vapply(a,`[[`,integer(1),"records"))," improvement records across ",length(a)," fit/direction payloads were audited from archived arrays; maximum independent component gap ",
   format(max(vapply(a,`[[`,numeric(1),"independent_component_max_gap")),digits=17),". No model evaluations were used for this audit."))
}
writeLines(lines,file.path(study,"reports/formal_execution_report.md"))
assert_namespaces();cat("Metadata report updated; zero model evaluations.\n")
