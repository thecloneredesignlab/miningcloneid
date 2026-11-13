# Summarize collected immutable outputs and closed parent accounting only.
# No additional simulation or objective evaluation is performed.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);source(file.path(study,"releases/fresh_pipeline01/adapter.R"))
jobs <- read.delim(file.path(study,"tables/jobs.tsv"),sep="\t",quote="\"",colClasses="character",check.names=FALSE)
jobs$cores <- as.integer(jobs$cores)
ledger_charge <- vapply(jobs$actual_ch,parse_exported_double,numeric(1))
stopifnot(all(jobs$state%in%c("COMPLETED","FAILED")),all(is.finite(ledger_charge)))
raw <- read.delim(file.path(study,"tables/raw_accounting.txt"),header=FALSE,sep="|",quote="",fill=TRUE,colClasses="character")
parents <- raw[match(jobs$job_id,raw$V1),,drop=FALSE]
stopifnot(identical(parents$V1,jobs$job_id),all(as.integer(parents$V4)==jobs$cores))
jobs$core_seconds <- jobs$cores*as.integer(parents$V5)
jobs$actual_ch <- jobs$core_seconds/3600
stopifnot(all(abs(jobs$actual_ch-ledger_charge)<1e-12))
total_charge <- sum(jobs$core_seconds)/3600
tasks <- c("fresh_global_production","fresh_local_production","fresh_global_data_only","fresh_local_data_only")
rows <- parameters <- list()
contract <- readRDS(file.path(study,"manifests/contract.rds"))
for(i in seq_along(tasks)) {
 task <- tasks[i];root <- file.path(study,"fits",task)
 result <- read_kv(file.path(root,"result.tsv"));accepted <- readRDS(file.path(root,"accepted_candidate.rds"))
 check <- readRDS(file.path(study,"tests",paste0("check_",task),"fresh_acceptance.rds"))
 j <- jobs[jobs$task==task,,drop=FALSE]
 stopifnot(nrow(j)==1L,j$state=="COMPLETED",j$exit_code=="0:0",isTRUE(result$pass),
   abs(check$components$L-accepted$record$components$L)<1e-8)
 parent <- raw[raw$V1==as.character(j$job_id),,drop=FALSE]
 batch <- raw[raw$V1==paste0(j$job_id,".batch"),,drop=FALSE]
 stopifnot(nrow(parent)==1L,nrow(batch)==1L)
 c <- check$components
 local <- grepl("_local_",task,fixed=TRUE)
 native <- readRDS(file.path(root,"native_result.rds"))
 rows[[i]] <- data.frame(task=task,policy=accepted$policy,stage=if(local) "local" else "global",
   job_id=j$job_id,allocated_cores=j$cores,parent_elapsed_seconds=as.integer(parent$V5),
   allocated_core_hours=j$actual_ch,stage_seconds=result$wall_seconds,
   objective_calls=result$evaluations,aggregate_score_seconds=result$objective_seconds,
   native_generations=if(local) NA_integer_ else result$iterations,
   native_local_code=if(local) result$native_convergence else NA_integer_,
   native_message=if(local && !is.null(native$message)) native$message else "",
   budget_limited=result$budget_limited,objective=c$L,data=c$L_data,prior=c$L_prior,
   burden=c$L_b,ploidy=c$L_p,necrosis=c$L_n,
   burden_observations=c$n_burden_obs_total,chromosome_values=c$n_ploidy_obs_total,necrosis_harvests=c$n_necrosis_obs_total,
   fixed_log10_density=result$fixed_log10_density,
   fresh_score_gap=check$components$L-accepted$record$components$L,
   slurm_batch_maxrss_reported=batch$V6,
   improvement_records=length(list.files(file.path(root,"evaluators"),pattern="^improvement_.*[.]rds$",recursive=TRUE)))
 p <- check$parameters;span <- contract$upper-contract$lower
 natural <- p;islog <- startsWith(names(p),"log10_");natural[islog] <- 10^p[islog]
 parameters[[i]] <- data.frame(task=task,policy=accepted$policy,parameter=names(p),
   transformed=unname(p),natural=unname(natural),lower=unname(contract$lower),upper=unname(contract$upper),
   relative_position=unname((p-contract$lower)/span),
   near_lower=unname(abs(p-contract$lower)<=1e-6*span),near_upper=unname(abs(p-contract$upper)<=1e-6*span),
   fixed=names(p)%in%names(accepted$fixed))
}
summary <- do.call(rbind,rows)
write_tsv(summary,file.path(study,"tables/corrected_timing_summary.tsv"))
write_tsv(do.call(rbind,parameters),file.path(study,"tables/corrected_candidate_parameters.tsv"))
role <- do.call(rbind,lapply(unique(jobs$role),function(r) {
 x <- jobs[jobs$role==r,,drop=FALSE]
 data.frame(role=r,attempts=nrow(x),failed=sum(x$state!="COMPLETED"),allocated_core_seconds=sum(x$core_seconds),actual_core_hours=sum(x$core_seconds)/3600)
}))
write_tsv(role,file.path(study,"tables/lifetime_role_accounting.tsv"))
budget <- read_tsv(file.path(study,"manifests/budget.tsv"))
stopifnot(total_charge<=192)
for(r in role$role) {
 b <- budget[budget$role==r,,drop=FALSE];x <- role[role$role==r,,drop=FALSE]
 stopifnot(x$attempts<=b$attempts,x$actual_core_hours<=b$cap_core_hours)
}
f <- function(x,digits=9) format(x,digits=digits,scientific=FALSE,trim=TRUE)
lines <- c("# Corrected constrained timing results", "",
 "These are warm/reference timing chains with density fixed at log10(rho)=5.198970004336019, approximately 158113.883 cells/mm^3. They are not unrestricted baselines, independent-start recovery or density profile results.", "",
 "| Policy | Stage | Objective | Calls | Stage seconds | Parent seconds | Allocated core-hours | Native result | Fresh gap |",
 "|---|---|---:|---:|---:|---:|---:|---|---:|")
for(i in seq_len(nrow(summary))) {
 x <- summary[i,]
 native_label <- if(x$stage=="global") paste(x$native_generations,"generations") else paste("code",x$native_local_code)
 lines <- c(lines,paste0("| ",paste(c(x$policy,x$stage,f(x$objective),x$objective_calls,f(x$stage_seconds),x$parent_elapsed_seconds,f(x$allocated_core_hours),native_label,f(x$fresh_score_gap)),collapse=" | ")," |"))
}
lines <- c(lines,"",paste0("Total lifetime actual allocated usage, including every failed replay, failed parallel fixture, legacy timing run, repair and independent acceptance: **",f(total_charge,17)," core-hours** (",sum(jobs$core_seconds)," integer core-seconds). All ",nrow(jobs)," owned allocations are closed. The separately approved ceiling remains 192; remaining unspent capacity is ",f(192-total_charge,17)," core-hours. Role and attempt caps were checked against the amended budget."),"",
 "Stage wall time includes optimization and in-stage acceptance/bookkeeping. Parent elapsed time includes launch and compilation. Aggregate score seconds sum evaluator timers and exclude improvement-prediction exports. Core-hours come only from closed Slurm parent AllocCPUS x ElapsedRaw / 3600; batch/extern steps are retained as evidence and are not charged again. Slurm batch MaxRSS is exported verbatim in the TSV; it is not substituted for an aggregate-memory claim.","",
 "Every new score evaluation uses the validated clean-parent child adapter with the pinned SIF and unchanged numerical source. NP=256, a maximum of 1000 DE generations, relative tolerance 1e-4 and patience 25 were retained for this pilot. Local runs retain native L-BFGS-B and its maximum 15000 iterations. Reported candidate scores refer to the best archived feasible incumbent, with native routine termination reported separately. Budget-limit flags, native termination, fresh acceptance and scientific search adequacy are distinct. A fresh-valid candidate with a nonzero native local code does not establish local stability or global optimality.","",
 "Both global candidates were independently accepted before separate local allocations. Each final local candidate was independently accepted in another allocation, with the same fixed density, complete prediction tables and independent burden/ploidy/necrosis arithmetic. Full-precision improvement vectors, component scores, native live/dead state arrays and oxygen outputs remain in immutable RDS records. The companion archive audit checks every saved improvement without new model calls.","",
 "See source_and_observation_audit.md, engineering_validation_amendment02.md and legacy_timing_findings.md for source identity, observation exclusions, initial-scale conditioning, history-cache diagnosis and preserved failed attempts. No formal baseline, bulk density profile, initial-size extension or culture work was launched by this pilot.")
writeLines(lines,file.path(study,"reports/corrected_timing_results.md"))
print(summary[,c("task","objective","objective_calls","stage_seconds","allocated_core_hours","native_local_code","fresh_score_gap")],row.names=FALSE)
cat("Lifetime allocated core-hours:",format(total_charge,digits=17),"\n")
